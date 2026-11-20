"""Chapman-Enskog binary diffusion and a closed diffusion cell (page A4.6).

Chapman and Cowling (1970): the first and second approximations to the binary
diffusion coefficient, eqs. (14.2, 1), (14.2, 4) and (14.21, 1-6), evaluated
from tabulated collision integrals; the 12,6 collision integrals recomputed by
quadrature of the classical deflection angle (no table used); and a closed
one-dimensional diffusion cell (Loschmidt tube) with a composition-dependent
diffusivity, solved with pymrm by backward Euler and Newton.

Units are c.g.s. throughout, as in the book: cm, g, s, erg.
"""
import math

import numpy as np
import scipy.sparse as sps
from scipy.constants import N_A, atm
from scipy.constants import k as KB_SI
from scipy.interpolate import CubicSpline
from scipy.sparse.linalg import splu
from pymrm import (construct_boundary_value_matrices, construct_coefficient_matrix,
                   construct_div, construct_grad, newton)
from pymrm.checks import find_roots

KB = KB_SI * 1e7                 # erg/K
T_STP = 273.15                   # K, section 14.32: "reduced to 0 degC"
P_ATM = atm * 10.0               # dyn/cm^2
N_STP = P_ATM / (KB * T_STP)     # molecules/cm^3 at 1 atm and 0 degC

# zero flux through a closed end, outward normal n: dx1/dn = 0 -> 1*dx1/dn + 0*x1 = 0
CLOSED_END = {"a": 1.0, "b": 0.0, "d": 0.0}


def sigma12_from_d(d12, m1, m2, T, n):
    """Invert eq. (14.2, 1): the rigid-sphere sigma_12 a measured D12 implies.

    Parameters
    ----------
    d12 : array_like
        Diffusion coefficients, cm^2/s.
    m1, m2 : array_like
        Molecular masses, g.
    T : float
        Temperature, K.
    n : float
        Number density, 1/cm^3.

    Returns
    -------
    numpy.ndarray
        sigma_12 in units of 1e-8 cm.
    """
    return np.sqrt(3.0 / (8.0 * n * np.asarray(d12, float))
                   * np.sqrt(KB * T * (m1 + m2) / (2.0 * np.pi * m1 * m2))) * 1e8


def sigma_from_viscosity(mu, M, T=T_STP, corr=1.016):
    """Eq. (12.1, 6): the rigid-sphere diameter from a viscosity.

    ``mu = corr * (5/16) sqrt(k m T / pi) / sigma^2``, with the 1.016 of the
    exact rigid-sphere value (section 12.1).

    Parameters
    ----------
    mu : array_like
        Viscosity, poise.
    M : array_like
        Molar mass, g/mol.
    T : float
        Temperature, K.
    corr : float
        Correction factor of the first approximation.

    Returns
    -------
    numpy.ndarray
        sigma in units of 1e-8 cm.
    """
    m = np.asarray(M, float) / N_A
    return np.sqrt(corr * (5.0 / 16.0) * np.sqrt(KB * m * T / np.pi) / np.asarray(mu, float)) * 1e8


def delta_generic(x1, M1, M2, P1, P2, A, B, C):
    """Eqs. (14.21, 2-6): Delta with ``[D12]_2 = [D12]_1 / (1 - Delta)``.

    Parameters
    ----------
    x1 : array_like
        Mole fraction of gas 1.
    M1, M2 : float
        Mass fractions ``m_i / (m_1 + m_2)``.
    P1, P2 : float
        ``M_i E / [mu_i]_1``.
    A, B, C : float
        Ratios of collision integrals, eq. (9.8, 7).

    Returns
    -------
    numpy.ndarray or float
        Delta at ``x1``.
    """
    num, den = delta_polynomials(M1, M2, P1, P2, A, B, C)
    return num(x1) / den(x1)


def delta_polynomials(M1, M2, P1, P2, A, B, C):
    """Numerator and denominator of Delta as polynomials in ``x1`` (``x2 = 1 - x1``).

    Returns
    -------
    tuple of numpy.polynomial.Polynomial
        ``(numerator, denominator)``.
    """
    P = np.polynomial.Polynomial
    x1, x2 = P([0.0, 1.0]), P([1.0, -1.0])
    P12 = 3 * (M1 - M2) ** 2 + 4 * M1 * M2 * A
    Q1 = P1 * (6 * M2 ** 2 + 5 * M1 ** 2 - 4 * M1 ** 2 * B + 8 * M1 * M2 * A)
    Q2 = P2 * (6 * M1 ** 2 + 5 * M2 ** 2 - 4 * M2 ** 2 * B + 8 * M1 * M2 * A)
    Q12 = 3 * (M1 - M2) ** 2 * (5 - 4 * B) + 4 * M1 * M2 * A * (11 - 4 * B) + 2 * P1 * P2
    num = 5 * C ** 2 * (M1 ** 2 * P1 * x1 ** 2 + M2 ** 2 * P2 * x2 ** 2 + P12 * x1 * x2)
    den = Q1 * x1 ** 2 + Q2 * x2 ** 2 + Q12 * x1 * x2
    return num, den


class CollisionTable:
    """Table 6 of the book as cubic splines in ``ln(kT/eps_12)``.

    Parameters
    ----------
    t_star, w11, a, b, c : array_like
        The columns kT/eps_12, W^(1)_12(1), A, B and C.
    """

    def __init__(self, t_star, w11, a, b, c):
        self.lt = np.log(np.asarray(t_star, float))
        self.spl_w = CubicSpline(self.lt, np.log(np.asarray(w11, float)))  # log-log
        self.spl_a = CubicSpline(self.lt, np.asarray(a, float))
        self.spl_b = CubicSpline(self.lt, np.asarray(b, float))
        self.spl_c = CubicSpline(self.lt, np.asarray(c, float))

    def w11(self, t_star):
        """W^(1)_12(1) at ``t_star``."""
        return np.exp(self.spl_w(np.log(t_star)))

    def abc(self, t_star):
        """A, B, C at ``t_star``."""
        lt = np.log(t_star)
        return float(self.spl_a(lt)), float(self.spl_b(lt)), float(self.spl_c(lt))


class ChapmanEnskog:
    """First and second approximations to D12 for named gases.

    Parameters
    ----------
    molar_mass, mu : dict
        Molar mass (g/mol) and viscosity at S.T.P. (poise) per formula.
    eps, sig : dict
        12,6 force constants per formula: eps/k (K) and sigma (1e-8 cm).
    table : CollisionTable
        The tabulated collision integrals.
    """

    def __init__(self, molar_mass, mu, eps, sig, table):
        self.M, self.mu, self.eps, self.sig, self.table = molar_mass, mu, eps, sig, table

    def pair_constants(self, g1, g2, eps=None, sig=None, eps_rule="geometric",
                       sig_rule="arithmetic", sig_scale=1.0):
        """Masses (g), eps_12/k (K) and sigma_12 (cm) of a pair by the combination rules.

        ``eps_rule`` and ``sig_rule`` are ``"geometric"`` or ``"arithmetic"``;
        the book's are a geometric eps_12 (p. 262) and an arithmetic sigma_12
        (eq. 14.2, 1). ``sig_scale`` multiplies sigma_12.
        """
        eps = self.eps if eps is None else eps
        sig = self.sig if sig is None else sig
        rule = {"geometric": lambda a, b: math.sqrt(a * b),
                "arithmetic": lambda a, b: 0.5 * (a + b)}
        return (self.M[g1] / N_A, self.M[g2] / N_A, rule[eps_rule](eps[g1], eps[g2]),
                sig_scale * rule[sig_rule](sig[g1], sig[g2]) * 1e-8)

    def t_star(self, g1, g2, T=T_STP, **kw):
        """kT/eps_12 of the pair."""
        return T / self.pair_constants(g1, g2, **kw)[2]

    def d12_first(self, g1, g2, T=T_STP, n=N_STP, W=None, **kw):
        """Eq. (14.2, 4); ``W = 1`` with Table 11 diameters is eq. (14.2, 1).

        Parameters
        ----------
        g1, g2 : str
            Formulas of the two gases.
        T : float
            Temperature, K.
        n : float
            Number density, 1/cm^3.
        W : float, optional
            W^(1)_12(1) to use instead of the table.
        **kw
            Passed to :meth:`pair_constants`.

        Returns
        -------
        float
            [D12]_1 in cm^2/s.
        """
        m1, m2, e12, s12 = self.pair_constants(g1, g2, **kw)
        w = self.table.w11(T / e12) if W is None else W
        return float(3.0 / (8.0 * n * s12 ** 2 * w)
                     * np.sqrt(KB * T * (m1 + m2) / (2.0 * np.pi * m1 * m2)))

    def delta_args(self, g1, g2, T=T_STP, A=None, B=None, C=None, swap_masses=False):
        """The arguments ``(M1, M2, P1, P2, A, B, C)`` of :func:`delta_generic` for a pair.

        ``E = kT / (8 M1 M2 Omega^(1)_12(1))`` (eq. 9.8, 8) with Omega
        eliminated through eq. (14.2, 4). A, B, C default to Table 6.
        """
        m1, m2, e12, s12 = self.pair_constants(g1, g2)
        if swap_masses:
            m1, m2 = m2, m1
        ts = T / e12
        w = float(self.table.w11(ts))
        a, b, c = self.table.abc(ts)
        m0 = m1 + m2
        E = (m0 / (4.0 * s12 ** 2 * w)) * np.sqrt(KB * T * m0 / (2.0 * np.pi * m1 * m2))
        M1, M2 = m1 / m0, m2 / m0
        return (M1, M2, M1 * E / self.mu[g1], M2 * E / self.mu[g2],
                a if A is None else A, b if B is None else B, c if C is None else C)

    def delta(self, g1, g2, x1, **kw):
        """Delta of eqs. (14.21, 2-6) at mole fraction ``x1`` of gas 1."""
        return delta_generic(np.asarray(x1, float), *self.delta_args(g1, g2, **kw))

    def d12_second(self, g1, g2, x1, **kw):
        """[D12]_2 = [D12]_1 / (1 - Delta), eq. (14.21, 1), in cm^2/s."""
        return self.d12_first(g1, g2) / (1.0 - self.delta(g1, g2, x1, **kw))


# ---- 12,6 collision integrals by quadrature --------------------------------

def turning_point(e_star, b):
    """Outermost root r0 of ``1 - b^2/r^2 - 4 (r^-12 - r^-6)/E* = 0`` (reduced units).

    Multiplied by ``r^12`` the condition is the sextic
    ``s^6 - b^2 s^5 + (4/E*) s^3 - 4/E* = 0`` in ``s = r^2``. Its roots are the
    eigenvalues of the companion matrix; the largest real one is the turning
    point (beyond it the polynomial is positive), polished by Newton steps.

    Parameters
    ----------
    e_star, b : numpy.ndarray
        Reduced energies and impact parameters, broadcast together.

    Returns
    -------
    r0 : numpy.ndarray
        Turning points.
    n_pos : numpy.ndarray
        Number of positive real roots of the sextic (1 or 3, by Descartes).
    """
    e_star, b = np.broadcast_arrays(np.asarray(e_star, float), np.asarray(b, float))
    shape = e_star.shape
    e, bb = e_star.ravel(), b.ravel()
    q = 4.0 / e
    coef = np.zeros((e.size, 6))          # monic: s^6 + c5 s^5 + ... + c0
    coef[:, 0] = -q                       # s^0
    coef[:, 3] = q                        # s^3
    coef[:, 5] = -bb ** 2                 # s^5
    comp = np.zeros((e.size, 6, 6))
    comp[:, 1:, :-1] = np.eye(5)
    comp[:, :, -1] = -coef
    lam = np.linalg.eigvals(comp)
    real = np.abs(lam.imag) <= 1e-7 * np.maximum(np.abs(lam), 1.0)
    pos = real & (lam.real > 0.0)
    s = np.where(pos, lam.real, -np.inf).max(axis=1)
    for _ in range(3):                    # Newton polish on the sextic
        p = s ** 6 - bb ** 2 * s ** 5 + q * s ** 3 - q
        dp = 6 * s ** 5 - 5 * bb ** 2 * s ** 4 + 3 * q * s ** 2
        s = s - p / dp
    return np.sqrt(s).reshape(shape), pos.sum(axis=1).reshape(shape)


def deflection_angle(e_star, b, n_g=96):
    """Classical deflection angle chi for the reduced 12,6 potential.

    ``chi = pi - 2 b int_{r0}^inf dr / (r^2 sqrt(f))``; the substitution
    ``u = r0/r = 1 - s^2`` removes the inverse-square-root singularity at r0,
    and Gauss-Legendre with ``n_g`` points integrates in s.

    Parameters
    ----------
    e_star, b : numpy.ndarray
        Reduced energies and impact parameters (same shape, 1-D).
    n_g : int
        Gauss points of the radial integral.

    Returns
    -------
    numpy.ndarray
        chi, same shape as ``b``.
    """
    r0, _ = turning_point(e_star, b)
    xs, ws = np.polynomial.legendre.leggauss(n_g)
    s, w = 0.5 * (xs + 1.0), 0.5 * ws
    rr = r0[:, None] / (1.0 - s * s)[None, :]
    fv = 1.0 - (b[:, None] / rr) ** 2 - 4.0 * (rr ** -12 - rr ** -6) / e_star[:, None]
    fv = np.maximum(fv, 1e-300)
    return np.pi - 2.0 * b / r0 * np.sum(w * 2.0 * s / np.sqrt(fv), axis=1)


def omega_11_22(t_star, n_E=72, x_max=40.0, n_b=120, n_panel=6, b_max=6.0, n_g=96):
    """Reduced collision integrals Omega^(1,1)* and Omega^(2,2)* (1 for rigid spheres).

    The transport cross-sections are integrals over the impact parameter on
    ``n_panel`` panels of ``n_b`` Gauss points up to ``b_max`` (panels resolve
    the orbiting singularity at low energy); the energy average uses
    Gauss-Legendre with ``n_E`` points on ``[0, x_max]`` of ``x = E/kT``.
    No table is consulted.

    Parameters
    ----------
    t_star : float
        kT/eps_12.

    Returns
    -------
    tuple of float
        ``(Omega11*, Omega22*)``.
    """
    xe, we = np.polynomial.legendre.leggauss(n_E)
    x = 0.5 * x_max * (xe + 1.0)
    wx = 0.5 * x_max * we
    edges = np.linspace(0.0, b_max, n_panel + 1)
    xb, wb = np.polynomial.legendre.leggauss(n_b)
    bs = np.concatenate([0.5 * (hi - lo) * (xb + 1.0) + lo for lo, hi in zip(edges[:-1], edges[1:])])
    wbs = np.concatenate([0.5 * (hi - lo) * wb for lo, hi in zip(edges[:-1], edges[1:])])
    E = np.repeat(x * t_star, bs.size)
    B = np.tile(bs, x.size)
    c = np.cos(deflection_angle(E, B, n_g=n_g)).reshape(x.size, bs.size)
    q1 = 2.0 * np.sum(wbs * (1.0 - c) * bs, axis=1)              # denominator 1 for l = 1
    q2 = 2.0 * np.sum(wbs * (1.0 - c ** 2) * bs, axis=1) / (2.0 / 3.0)
    e = np.exp(-x)
    return (float(np.sum(wx * e * x ** 2 * q1) / 2.0),
            float(np.sum(wx * e * x ** 3 * q2) / 6.0))


# ---- the closed diffusion cell ----------------------------------------------

def analytic_split(D, L, t, n_terms=2000):
    """Top-minus-bottom mean composition of a closed tube with constant D.

    The cosine series of the initial step at ``L/2``:
    ``sum over odd k of 8/(pi k)^2 exp(-D (k pi/L)^2 t)``, ``n_terms`` terms.
    """
    k = np.arange(1, 2 * n_terms, 2)
    return float(np.sum(8.0 / (np.pi * k) ** 2 * np.exp(-D * (k * np.pi / L) ** 2 * t)))


def split_series_tail(D, L, t, n_terms=2000):
    """Bound on the series tail after ``n_terms`` terms, relative to the sum.

    The terms fall with a ratio that itself falls, so the tail is at most
    ``t_{N+1} / (1 - t_{N+2}/t_{N+1})``.
    """
    k1, k2 = 2 * n_terms + 1, 2 * n_terms + 3
    t1 = 8.0 / (np.pi * k1) ** 2 * np.exp(-D * (k1 * np.pi / L) ** 2 * t)
    t2 = 8.0 / (np.pi * k2) ** 2 * np.exp(-D * (k2 * np.pi / L) ** 2 * t)
    tail = t1 / (1.0 - t2 / t1) if t1 > 0.0 else 0.0
    return tail / analytic_split(D, L, t, n_terms)


def invert_split(split, L, t, lo, hi, n_terms=2000):
    """The constant D whose closed form gives ``split``: one root on ``[lo, hi]``.

    Raises unless :func:`pymrm.checks.find_roots` returns exactly one root.
    """
    r = find_roots(lambda D: analytic_split(D, L, t, n_terms) - split, lo, hi, xtol=1e-15)
    if len(r["roots"]) != 1:
        raise RuntimeError(f"apparent D: {len(r['roots'])} roots on [{lo:.3g}, {hi:.3g}]")
    return r["roots"][0]


class DiffusionCell:
    """A closed tube, gas 1 above ``z = L/2`` and gas 2 below at ``t = 0``.

    ``dx1/dt = d/dz (D(x1) dx1/dz)`` with zero flux at both ends, on ``n_z``
    equal cells, by the theta method (``theta = 1``: backward Euler). The face
    diffusivity is ``D`` at the face composition: boundary faces from
    ``construct_boundary_value_matrices``, interior faces from the mean of the
    two neighbours (``face_rule="mean"``) or the left cell (``"left"``).

    Parameters
    ----------
    L : float
        Tube length, cm.
    n_z : int
        Number of cells.
    d : float or callable
        A constant diffusivity (the problem is then linear and solved by one
        factorisation), or ``d(x1)`` in cm^2/s (solved by Newton per step).
    nu : int
        Geometry index of ``construct_div`` (0: slab).
    face_rule : str
        ``"mean"`` or ``"left"``.
    theta : float
        Implicitness of the time step.
    step_offset : float
        Initial step at ``L/2 + step_offset * dz``.
    tol : float
        Newton step tolerance (unknowns are mole fractions of order one).
    """

    def __init__(self, L, n_z, d, nu=0, face_rule="mean", theta=1.0, step_offset=0.0,
                 tol=1e-12):
        self.L, self.n_z, self.d, self.theta, self.tol = L, n_z, d, theta, tol
        self.shape = (n_z, 1)
        self.x_f = np.linspace(0.0, L, n_z + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        bc = (CLOSED_END, CLOSED_END)
        # constant operators, assembled once
        self.grad, gb = construct_grad(self.shape, self.x_f, self.x_c, bc)
        self.grad_bc = gb.toarray().reshape((-1, 1))
        self.div = construct_div(self.shape, self.x_f, nu=nu)
        bl, _ = construct_boundary_value_matrices(self.shape, self.x_f, self.x_c, bc[0], bound_id=0)
        br, _ = construct_boundary_value_matrices(self.shape, self.x_f, self.x_c, bc[1], bound_id=1)
        if face_rule == "mean":
            w = 0.5 * (sps.eye(n_z - 1, n_z) + sps.eye(n_z - 1, n_z, k=1))
        elif face_rule == "left":
            w = sps.eye(n_z - 1, n_z)
        else:
            raise ValueError(f"face_rule {face_rule!r}")
        self.face = sps.vstack([bl, w, br], format="csc")   # zero-flux faces: no bc vector
        self.eye = sps.eye(n_z, format="csc")
        self.x0 = np.where(self.x_c < 0.5 * L + step_offset * (L / n_z), 1.0, 0.0).reshape(self.shape)

    def _flux(self, x):
        """Face fluxes -D(x_f) dx1/dz and the face values and gradients."""
        xf = self.face @ x
        g = self.grad @ x + self.grad_bc
        return -self.d(xf) * g, xf, g

    def _residual(self, x_old, dt):
        div_old = self.div @ self._flux(x_old)[0] if self.theta < 1.0 else 0.0

        def fun(y):
            y = np.reshape(y, (-1, 1))
            flux, xf, g = self._flux(y)
            r = (y - x_old) / dt + self.theta * (self.div @ flux) + (1.0 - self.theta) * div_old
            h = 1e-6                                   # D varies by a few % over x1 in [0, 1]
            dd = (self.d(xf + h) - self.d(xf - h)) / (2.0 * h)
            dflux = -(construct_coefficient_matrix(self.d(xf)) @ self.grad
                      + construct_coefficient_matrix(dd * g) @ self.face)
            return r, self.eye / dt + self.theta * (self.div @ dflux)
        return fun

    def solve(self, t_end, n_t):
        """Mole-fraction profile at ``t_end`` after ``n_t`` equal steps; shape ``(n_z, 1)``.

        Raises unless every Newton solve reports success.
        """
        dt = t_end / n_t
        x = self.x0.copy()
        if not callable(self.d):                      # constant D: linear, one factorisation
            lap = self.div @ (-self.d * self.grad)
            lu = splu(self.eye / dt + self.theta * lap)
            expl = self.eye / dt - (1.0 - self.theta) * lap
            for _ in range(n_t):
                x = lu.solve(expl @ x).reshape(self.shape)
            return x
        for _ in range(n_t):
            res = newton(self._residual(x, dt), x, tol=self.tol, stagnation=True, maxfev=60)
            if not res.success:
                raise RuntimeError(f"cell step failed: {res.message}")
            x = np.reshape(res.x, self.shape)
        return x

    def split(self, x):
        """Mean composition of the top half minus that of the bottom half."""
        lo = self.x_c < 0.5 * self.L
        x = np.ravel(x)
        return float(x[lo].mean() - x[~lo].mean())
