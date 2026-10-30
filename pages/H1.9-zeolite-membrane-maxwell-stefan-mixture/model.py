"""Maxwell-Stefan mixture permeation across a zeolite membrane (Krishna & Baur 2003).

Model code of gallery page H1.9. The equation numbers are those of the review,
Krishna, R. & Baur, R., Separation and Purification Technology 33(3) 213-254
(2003).

Contents
--------
DSL
    Dual-site Langmuir isotherm (eq. 32) with its surface potential (eq. 64).
IAST
    Ideal adsorbed solution theory, eqs. (62)-(65), forward (pressures to
    loadings) and inverse (loadings to pressures), and the matrix of
    thermodynamic factors (eq. 38) by implicit differentiation.
B_matrix
    The friction matrix of eq. (44) with the exchange coefficients of eq. (43).
gamma_numerical
    The thermodynamic-factor matrix by central differences of the IAST inverse,
    the independent route for ``IAST.inverse_and_gamma``.
friction_solve
    The friction system of eq. (42) solved matrix-free by GMRES, the independent
    route for ``B_matrix``.
Membrane
    Eqs. (57)-(59) with the flux relations of eq. (45) on a pymrm slab grid:
    steady state by ``pymrm.newton`` (with ``pymrm.continuation`` in the upstream
    pressure as the fallback route), and the transient by backward Euler.
selectivity, locate_peak
    The permeation selectivity of eq. (60), and the maximum of a sampled curve.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, gmres, spsolve

from pymrm import (NewtonBreakdown, NumJac, construct_coefficient_matrix, construct_div,
                   construct_grad, continuation, newton)

__all__ = ["DSL", "IAST", "B_matrix", "gamma_numerical", "friction_solve", "Membrane",
           "selectivity", "locate_peak"]


class DSL:
    """Dual-site Langmuir isotherm, eq. (32), and its surface potential, eq. (64).

    ``Theta(P) = Th_sat_A b_A P / (1 + b_A P) + Th_sat_B b_B P / (1 + b_B P)``;
    ``Th_sat_B = 0`` gives the single-site Langmuir.

    Parameters
    ----------
    b_A, b_B : float
        Site affinities (1/Pa).
    Th_sat_A, Th_sat_B : float
        Site saturation capacities (molecules per unit cell).

    Attributes
    ----------
    Th_sat : float
        Total saturation capacity.
    H : float
        Henry coefficient ``dTheta/dP`` at ``P -> 0`` (molecules per unit cell per Pa).
    """

    def __init__(self, b_A, Th_sat_A, b_B=0.0, Th_sat_B=0.0):
        self.bA, self.tA = float(b_A), float(Th_sat_A)
        self.bB, self.tB = float(b_B), float(Th_sat_B)
        self.Th_sat = self.tA + self.tB
        self.H = self.tA * self.bA + self.tB * self.bB
        self._lp = np.linspace(-45.0, 45.0, 9001)      # seed table for the inverse
        self._tab = self.psi(10.0 ** self._lp)

    def theta(self, P):
        """Loading at pressure ``P`` (Pa)."""
        r = self.tA * self.bA * P / (1 + self.bA * P)
        return r + (self.tB * self.bB * P / (1 + self.bB * P) if self.tB > 0 else 0.0)

    def dtheta(self, P):
        """``dTheta/dP``."""
        r = self.tA * self.bA / (1 + self.bA * P) ** 2
        return r + (self.tB * self.bB / (1 + self.bB * P) ** 2 if self.tB > 0 else 0.0)

    def psi(self, P):
        """Reduced surface potential ``int_0^P Theta(P') dP'/P'`` (eq. 64 over rho k_B T)."""
        r = self.tA * np.log1p(self.bA * P)
        return r + (self.tB * np.log1p(self.bB * P) if self.tB > 0 else 0.0)

    def P_of_psi(self, ps):
        """Pressure at which the pure component has surface potential ``ps``.

        ``psi`` is increasing in ``P`` with ``dpsi/dlnP = Theta``: a table seed in
        ``log10 P`` and six Newton steps in ``ln P``.
        """
        ps = np.asarray(ps, float)
        P = 10.0 ** np.interp(ps, self._tab, self._lp)
        for _ in range(6):
            st = np.clip(-(self.psi(P) - ps) / np.maximum(self.theta(P), 1e-300), -50.0, 50.0)
            P = np.clip(P * np.exp(st), 1e-300, 1e300)
        return P


def _solve_psi(g, dg_dlnpsi, m, incr, seed=None, rtol=1e-14, maxit=50):
    """Monotone scalar root of ``g(psi)`` at each of ``m`` points.

    The bracket is widened until it brackets. Without a seed: 30 log-bisections
    on ``[1e-14, 1e2]``, then five Newton steps in ``ln psi`` kept inside the
    bracket. With a seed (the solution of a previous call at nearby states): a
    bracket of six decades around it and safeguarded Newton in ``ln psi`` from
    the seed, the bracket shrinking with every sign of ``g`` seen and a
    bisection taken whenever Newton would leave it, until the step is below
    ``rtol``. A seeded solve that has not met ``rtol`` after ``maxit`` steps
    raises ``FloatingPointError``.
    """
    lo = np.full(m, 1e-14) if seed is None else np.maximum(seed * 1e-3, 1e-300)
    hi = np.full(m, 1e2) if seed is None else seed * 1e3
    for _ in range(40):
        bad = (g(lo) > 0) if incr else (g(lo) < 0)
        if not np.any(bad):
            break
        lo = np.where(bad, lo / 100, lo)
    for _ in range(40):
        bad = (g(hi) < 0) if incr else (g(hi) > 0)
        if not np.any(bad):
            break
        hi = np.where(bad, hi * 100, hi)
    if seed is None:
        for _ in range(30):
            mid = np.sqrt(lo * hi)
            gm = g(mid)
            up = (gm < 0) if incr else (gm > 0)
            lo = np.where(up, mid, lo)
            hi = np.where(up, hi, mid)
        psi = np.sqrt(lo * hi)
        for _ in range(5):
            psi = np.clip(psi * np.exp(np.clip(-g(psi) / dg_dlnpsi(psi), -2.0, 2.0)), lo, hi)
        return psi
    psi = np.clip(seed, lo, hi)
    for _ in range(maxit):
        gv = g(psi)
        up = (gv < 0) if incr else (gv > 0)          # the root lies above psi
        lo = np.where(up, psi, lo)
        hi = np.where(up, hi, psi)
        step = -gv / dg_dlnpsi(psi)
        new = psi * np.exp(np.clip(step, -2.0, 2.0))
        out = ~((new > lo) & (new < hi))
        new = np.where(out & (gv != 0), np.sqrt(lo * hi), np.where(gv == 0, psi, new))
        done = np.max(np.abs(np.log(new / psi))) < rtol
        psi = new
        if done:
            return psi
    raise FloatingPointError("IAST surface potential not converged from its seed")


class IAST:
    """Ideal adsorbed solution theory, eqs. (62)-(65), in both directions.

    Parameters
    ----------
    isos : list of DSL
        Pure-component isotherms, one per species (anything with ``theta``,
        ``dtheta``, ``psi`` and ``P_of_psi``).
    """

    def __init__(self, isos):
        self.isos = isos
        self.n = len(isos)
        self.psi_last = None

    def _P0(self, psi):
        return np.stack([self.isos[i].P_of_psi(psi) for i in range(self.n)], -1)

    def _T0(self, P0):
        return np.stack([self.isos[i].theta(P0[..., i]) for i in range(self.n)], -1)

    def forward(self, p, seed=None):
        """Loadings from partial pressures; eq. (62) closes as ``sum p_i/P_i^0 = 1``.

        Parameters
        ----------
        p : numpy.ndarray
            Partial pressures (Pa), shape ``(m, n)`` or ``(n,)``.

        Returns
        -------
        numpy.ndarray
            Loadings (molecules per unit cell), shape ``(m, n)``.
        """
        p = np.atleast_2d(np.asarray(p, float))

        def g(psi):
            return np.sum(p / self._P0(psi), -1) - 1.0        # decreasing in psi

        def dg(psi):
            P0 = self._P0(psi)
            return -psi * np.sum(p / (P0 * self._T0(P0)), -1)

        psi = _solve_psi(g, dg, p.shape[0], incr=False, seed=seed)
        self.psi_last = psi
        P0 = self._P0(psi)
        x = p / P0
        return x * (1.0 / np.sum(x / self._T0(P0), -1))[..., None]     # eq. (65)

    def _psi_of(self, Th, seed=None):
        """Surface potential from loadings. The adsorbed mole fractions and the
        total loading are known from ``Th``, so eq. (65) is a scalar monotone
        equation in the potential."""
        Th = np.atleast_2d(np.asarray(Th, float))
        Tt = np.sum(Th, -1)
        x = Th / Tt[..., None]

        def g(psi):
            return 1.0 / np.sum(x / self._T0(self._P0(psi)), -1) - Tt    # increasing

        def dg(psi):
            P0 = self._P0(psi)
            T0 = self._T0(P0)
            dT0 = np.stack([self.isos[i].dtheta(P0[..., i]) * P0[..., i]
                            for i in range(self.n)], -1) / T0
            S = np.sum(x / T0, -1)
            return psi * np.sum(x * dT0 / T0 ** 2, -1) / S ** 2

        psi = _solve_psi(g, dg, Th.shape[0], incr=True, seed=seed)
        self.psi_last = psi
        return psi, x

    def inverse(self, Th, seed=None):
        """Partial pressures (Pa) in equilibrium with loadings ``Th``, shape ``(m, n)``."""
        psi, x = self._psi_of(Th, seed)
        return x * self._P0(psi)

    def inverse_and_gamma(self, Th, Th_sat, seed=None, drop_mix_term=False):
        """Fugacities and the thermodynamic-factor matrix of eq. (38) from one solve.

        ``ln f_i = ln x_i + ln P_i^0(psi)``, so
        ``dln f_i/dTheta_j = delta_ij/Theta_i - 1/Theta_mix + (1/Theta_i^0) dpsi/dTheta_j``,
        with ``dpsi/dTheta_j`` from differentiating eq. (65).

        Parameters
        ----------
        Th : numpy.ndarray
            Loadings, shape ``(m, n)``.
        Th_sat : numpy.ndarray
            Saturation capacities, shape ``(n,)``.
        seed : numpy.ndarray, optional
            Surface potentials of a previous call, to narrow the bracket.
        drop_mix_term : bool, optional
            A defect for the page's break table: drop ``-1/Theta_mix``.

        Returns
        -------
        f : numpy.ndarray
            Fugacities (Pa), ``(m, n)``.
        Gamma : numpy.ndarray
            ``Gamma_ij = (theta_i/theta_j) dln f_i/dln theta_j``, ``(m, n, n)``.
        """
        Th = np.atleast_2d(np.asarray(Th, float))
        n = self.n
        psi, x = self._psi_of(Th, seed)
        P0 = self._P0(psi)
        T0 = self._T0(P0)
        dT0 = np.stack([self.isos[i].dtheta(P0[..., i]) for i in range(n)], -1) * P0 / T0
        Tt = np.sum(Th, -1)
        S = np.sum(x / T0, -1)
        dF_dpsi = np.sum(x * dT0 / T0 ** 2, -1) / S ** 2
        inv = 1.0 / T0
        dF_dT = -(inv / Tt[:, None] - np.sum(x * inv, -1)[:, None] / Tt[:, None]) \
            / S[:, None] ** 2 - 1.0
        dpsi_dT = -dF_dT / dF_dpsi[:, None]
        mix = 0.0 if drop_mix_term else 1.0 / Tt[:, None, None]
        dlnf = (np.eye(n)[None] / Th[:, :, None] - mix
                + dpsi_dT[:, None, :] / T0[:, :, None])
        return x * P0, (Th / Th_sat)[:, :, None] * dlnf * Th_sat[None, None, :]


def B_matrix(th, D_i, mode="vignes", swap_exponents=False):
    """The friction matrix of eq. (44), exchange coefficients from eq. (43).

    ``B_ii = 1/D_i + sum_(j != i) theta_j/D_ij``, ``B_ij = -theta_i/D_ij``,
    ``D_ij = D_i^(theta_i/(theta_i+theta_j)) D_j^(theta_j/(theta_i+theta_j))``.

    Parameters
    ----------
    th : numpy.ndarray
        Fractional occupancies, ``(m, n)``.
    D_i : numpy.ndarray
        Maxwell-Stefan diffusivities, ``(m, n)``.
    mode : {"vignes", "inf"}
        ``"inf"`` takes ``D_ij -> infinity``, eq. (47) (Habgood's model).
    swap_exponents : bool, optional
        A defect for the page's break table: the two exponents of eq. (43)
        exchanged.

    Returns
    -------
    numpy.ndarray
        ``(m, n, n)``.
    """
    m, n = th.shape
    B = np.zeros((m, n, n))
    for i in range(n):
        B[:, i, i] = 1.0 / D_i[:, i]
        if mode == "inf":
            continue
        for j in range(n):
            if j == i:
                continue
            xi = (th[:, j] if swap_exponents else th[:, i]) / np.maximum(th[:, i] + th[:, j],
                                                                         1e-300)
            Dij = D_i[:, i] ** xi * D_i[:, j] ** (1 - xi)
            B[:, i, i] += th[:, j] / Dij
            B[:, i, j] = -th[:, i] / Dij
    return B


def gamma_numerical(iast, Th, Th_sat, h=1e-5):
    """Thermodynamic-factor matrix by central differences of the IAST inverse in ``ln Theta``.

    Parameters
    ----------
    iast : IAST
    Th : numpy.ndarray
        Loadings, ``(m, n)``.
    Th_sat : numpy.ndarray
        Saturation capacities, ``(n,)``.
    h : float, optional
        Relative perturbation of each loading.

    Returns
    -------
    numpy.ndarray
        ``(m, n, n)``.
    """
    Th = np.atleast_2d(np.asarray(Th, float))
    m, n = Th.shape
    th = Th / Th_sat
    G = np.zeros((m, n, n))
    for j in range(n):
        Tp = Th.copy()
        Tp[:, j] *= (1 + h)
        Tm = Th.copy()
        Tm[:, j] *= (1 - h)
        d = (np.log(iast.inverse(Tp)) - np.log(iast.inverse(Tm))) \
            / (np.log1p(h) - np.log1p(-h))
        for i in range(n):
            G[:, i, j] = th[:, i] / th[:, j] * d[:, i]
    return G


def friction_solve(th, D_i, d, mode="vignes"):
    """Solve the friction system of eq. (42) for the scaled fluxes without forming ``[B]``.

    The matrix-vector product is the right-hand side of eq. (42) as printed,
    so this route does not share the ``B_ii``/``B_ij`` formulas of eq. (44).

    Parameters
    ----------
    th, D_i, d : numpy.ndarray
        Occupancies, diffusivities and driving forces, each ``(n,)``.
    mode : {"vignes", "inf"}

    Returns
    -------
    numpy.ndarray
        ``u_i = N_i/(rho Theta_i,sat)``, ``(n,)``.
    """
    n = len(th)

    def Dij(i, j):
        xi = th[i] / max(th[i] + th[j], 1e-300)
        return D_i[i] ** xi * D_i[j] ** (1 - xi)

    def mv(u):
        r = np.empty(n)
        for i in range(n):
            acc = u[i] / D_i[i]
            if mode != "inf":
                for j in range(n):
                    if j != i:
                        acc += (th[j] * u[i] - th[i] * u[j]) / Dij(i, j)
            r[i] = acc
        return r

    u, info = gmres(LinearOperator((n, n), matvec=mv), d, rtol=1e-14, atol=0.0)
    if info != 0:
        raise RuntimeError(f"GMRES did not converge (info = {info})")
    return u


class Membrane:
    """Multicomponent Maxwell-Stefan permeation across a zeolite slab, eqs. (57)-(59).

    ``dTheta_i/dt = -(1/rho) dN_i/dz`` with ``(N) = -rho [Theta_sat][B]^-1[Gamma] (dtheta/dz)``
    (eq. 45), and loadings fixed at both faces by equilibrium with the upstream
    and downstream partial pressures. State ``Theta``, shape ``(n_z, n_c)``:
    spatial axis first, species last. The face matrices ``[B]^-1[Gamma]`` depend
    on the loadings of the two neighbouring cells, so the steady Jacobian is the
    full-residual ``NumJac`` with a tridiagonal stencil in space and full
    coupling of the species axis. Fluxes are ``N/rho`` in molecules per unit cell
    times ``m/s`` on a grid of thickness ``delta``.

    Parameters
    ----------
    isos : list of DSL
        Pure-component isotherms.
    Th_sat : array_like
        Saturation capacities (molecules per unit cell).
    D0 : array_like
        Zero-loading Maxwell-Stefan diffusivities (m2/s, or any common unit).
    n_z : int, optional
        Number of cells across the membrane.
    Dij : {"vignes", "inf"}, optional
        Exchange coefficients by eq. (43), or infinite (eq. 47).
    loading_dep : callable, optional
        ``D_i(theta)/D_i(0)``; the default is eq. (17), constant.
    iso_model : {"iast", "langmuir"}, optional
        Mixture isotherm: IAST (eqs. 62-65) or the multicomponent Langmuir (eq. 39).
    delta : float, optional
        Membrane thickness (m); 1 for the thickness-free selectivities.
    pd_frac : float, optional
        Downstream partial pressures as this fraction of the upstream ones
        (the review's "vanishing values").
    nu, gamma_mode, swap_exponents
        Defects for the page's break table, not modelling options: ``nu=1``
        treats the slab as a cylinder, ``gamma_mode`` ``"identity"`` or
        ``"diagonal"`` replaces ``[Gamma]``, ``swap_exponents`` exchanges the
        exponents of eq. (43).
    """

    def __init__(self, isos, Th_sat, D0, n_z=30, Dij="vignes", loading_dep=None,
                 iso_model="iast", delta=1.0, pd_frac=1e-6, nu=0, gamma_mode="full",
                 swap_exponents=False):
        self.iast = IAST(isos)
        self.isos = isos
        self.n = len(isos)
        self.Th_sat = np.asarray(Th_sat, float)
        self.D0 = np.asarray(D0, float)
        self.n_z, self.Dij, self.iso_model, self.delta = n_z, Dij, iso_model, delta
        self.pd_frac, self.gamma_mode, self.swap_exponents = pd_frac, gamma_mode, swap_exponents
        self.loading_dep = loading_dep or (lambda th: np.ones_like(th))
        self.z_f = np.linspace(0.0, delta, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.shape = (n_z, self.n)                        # (space, species)
        self._seed = None
        self.nu = nu
        self._build_operators(nu)

    def _build_operators(self, nu):
        n = self.n
        # Outward normal, a dTheta/dn + b Theta = d; both faces Dirichlet (eqs. 57, 58):
        #   z = 0:     Theta_i = Theta_i,0  (equilibrium with the upstream gas)
        #   z = delta: Theta_i = Theta_i,d  (equilibrium with the swept downstream gas)
        # The face loadings change with the feed, so d = 1 multiplies an external
        # vector (shapes_d) and the operator is assembled once.
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
        self.grad, self.grad_bc_in, self.grad_bc_out = construct_grad(
            self.shape, self.z_f, self.z_c, bc, axis=0, shapes_d=((1, n), (1, n)))
        self.div = construct_div(self.shape, self.z_f, nu=nu, axis=0)   # nu = 0: slab
        # [B]^-1[Gamma] couples the species at each face: an (n, n) block per
        # face, placed by construct_coefficient_matrix with rows (face, i) and
        # columns (face, j).
        nf = self.n_z + 1
        self._face_shapes = ((nf, n, 1), (nf, 1, n))
        # the residual reads neighbouring cells through the face matrices
        self.numjac = NumJac(self.shape, axes_diagonals=[0])

    # -- isotherm: IAST (eqs. 62-65) or multicomponent Langmuir (eq. 39) --------
    def load(self, p):
        """Equilibrium loadings for partial pressures ``p``, ``(m, n)``."""
        p = np.atleast_2d(np.asarray(p, float))
        if self.iso_model == "iast":
            return self.iast.forward(p)
        den = 1 + np.sum([self.isos[i].bA * p[:, i] for i in range(self.n)], 0)
        return np.stack([self.Th_sat[i] * self.isos[i].bA * p[:, i] / den
                         for i in range(self.n)], -1)

    def gamma(self, Th):
        """Thermodynamic-factor matrices at loadings ``Th``, ``(m, n, n)``."""
        Th = np.atleast_2d(Th)
        if self.gamma_mode == "identity":
            return np.repeat(np.eye(self.n)[None], Th.shape[0], 0)
        if self.iso_model == "iast":
            _, G = self.iast.inverse_and_gamma(Th, self.Th_sat, seed=self._seed)
            self._seed = self.iast.psi_last
        else:
            th = Th / self.Th_sat                                       # eq. (39)
            G = np.eye(self.n)[None] + th[:, :, None] / (1 - th.sum(-1))[:, None, None]
        if self.gamma_mode == "diagonal":
            G = G * np.eye(self.n)[None]
        return G

    def boundary_loadings(self, p0, pd=None):
        """Face loadings upstream and downstream, each ``(n,)``."""
        p0 = np.asarray(p0, float)
        Th0 = self.load(p0[None, :])[0]
        pdv = np.asarray(pd, float) if pd is not None else self.pd_frac * p0
        Thd = self.load(pdv[None, :])[0]
        return Th0, Thd

    def _gbc(self, Th0, Thd):
        return (self.grad_bc_in @ Th0 + self.grad_bc_out @ Thd).ravel()

    def _M(self, Th, Th0, Thd):
        """Face-wise ``[Theta_sat][B]^-1[Gamma][Theta_sat]^-1`` as a block-diagonal matrix."""
        Thf = np.empty((self.n_z + 1, self.n))
        Thf[1:-1] = 0.5 * (Th[:-1] + Th[1:])
        Thf[0], Thf[-1] = Th0, Thd
        th = Thf / self.Th_sat
        G = self.gamma(Thf)
        B = B_matrix(th, self.D0 * self.loading_dep(th), self.Dij, self.swap_exponents)
        M = np.linalg.solve(B, G) * self.Th_sat[None, :, None] / self.Th_sat[None, None, :]
        return construct_coefficient_matrix(M, shape=self._face_shapes, format="csr")

    def flux(self, Th, Th0, Thd):
        """Face fluxes ``N/rho``, ``(n_z + 1, n)``."""
        gbc = self._gbc(Th0, Thd)
        return -(self._M(Th, Th0, Thd) @ (self.grad @ Th.reshape(-1) + gbc)).reshape(-1, self.n)

    def residual_log(self, L, Th0, Thd):
        """Steady residual ``-div N / Theta_sat`` in ``L = ln Theta``, and its Jacobian.

        Working in the logarithm keeps the loadings positive without a bound and
        makes the ``NumJac`` perturbation relative (the loadings span decades
        across the membrane). IAST has no solution at or above saturation, so a
        state with ``sum theta >= 1`` has a non-finite residual: ``newton``'s
        damping then shortens the step, and a ``continuation`` prediction there
        fails and halves the parameter step. (A raised exception would do for the
        damping, which evaluates trial points guarded, but not for a prediction,
        which ``newton`` evaluates as its initial guess.)
        """
        gbc = self._gbc(Th0, Thd)

        def res(Lf):
            Th = np.exp(Lf)
            if np.max(np.sum(Th / self.Th_sat, -1)) >= 1 - 1e-9:
                return np.full(self.shape, np.nan)       # outside the IAST domain
            r = self.div @ (self._M(Th, Th0, Thd) @ (self.grad @ Th.reshape(-1) + gbc))
            return r.reshape(self.shape) / self.Th_sat

        g, jac = self.numjac(res, L)
        return g.reshape(-1), jac

    def guess(self, Th0, Thd):
        """Initial profile, linear in the loadings between the two faces (exact in
        the Henry limit, where the flux law is Fick's law with a constant
        diffusivity)."""
        w = (self.z_c / self.delta)[:, None]
        return (1 - w) * Th0 + w * Thd

    def steady(self, p0, pd=None, Th_start=None, tol=1e-9, maxfev=30, accept_tol=1e-7,
               n_cont=4, force_continuation=False, damping="monotone"):
        """Steady state, with ``div N = 0`` checked on the result.

        Newton (``pymrm.newton``, monotone damping) in ``ln Theta`` from
        ``Th_start`` or ``guess``.
        If that does not converge (``success`` false, or ``NewtonBreakdown``), or
        ``force_continuation`` is set, the solve is
        repeated by ``pymrm.continuation`` in the upstream (and downstream)
        pressure, over ``n_cont`` factors from 0.01 to 1. A result is returned
        only when Newton reports success and the face fluxes are
        ``z``-independent to ``accept_tol``; otherwise ``RuntimeError``.

        Returns
        -------
        Th : numpy.ndarray
            Loadings, ``(n_z, n)``.
        J : numpy.ndarray
            Face fluxes ``N/rho``, ``(n_z + 1, n)``.
        """
        p0 = np.asarray(p0, float)
        # No bounds: the admissible set, sum theta < 1, is not a box. The
        # residual refuses states outside it and the damping shortens the step.
        nk = dict(tol=tol, maxfev=maxfev, damping=damping)
        Th0, Thd = self.boundary_loadings(p0, pd)
        self.continued = False
        ok = False
        if not force_continuation:
            L0 = np.log(Th_start if Th_start is not None else self.guess(Th0, Thd))
            try:
                res = newton(self.residual_log, L0, args=(Th0, Thd), **nk)
                ok = bool(res.success)
                self.result = res
            except NewtonBreakdown:          # statuses 4 to 6: no Newton step formed
                ok = False
        if not ok:
            self.continued = True
            facs = np.logspace(-2, 0, n_cont)

            def bnd(f):
                return self.boundary_loadings(p0 * f, None if pd is None else
                                              np.asarray(pd, float) * f)

            def fun(L, f):
                return self.residual_log(L, *bnd(f))

            L0 = np.log(self.guess(*bnd(facs[0])))
            cont = continuation(fun, L0, facs, newton_kwargs=nk)
            if not cont.success:
                raise RuntimeError(f"steady solve failed: {cont.message}")
            res = cont.results[-1]
            self.result = res
        Th = np.exp(res.x)
        J = self.flux(Th, Th0, Thd)
        self.newton_its = int(res.nit)
        # div N = 0 on the slab: the face flux is the same at every face. (With
        # the nu != 0 defect of the break table only newton's success is checked.)
        self.flux_spread = float(np.max(np.abs(J / J.mean(0) - 1.0)))
        if self.nu == 0 and self.flux_spread > accept_tol:
            raise RuntimeError(f"div N = 0 not met: flux spread {self.flux_spread:.1e}")
        return Th, J

    def transient(self, p0, t_out, pd=None, dt_frac=0.09, dt0=1e-6, tol=1e-11, maxit=40):
        """Permeation from a membrane equilibrated with the downstream gas.

        Backward Euler with steps ``dt_frac * t`` (at least ``dt0``), cut at each
        output time. Within a step the face matrices are frozen and the linear
        system re-solved until the loadings change by less than ``tol``
        (relative to the largest loading); a step whose inner iteration has not
        converged after ``maxit`` passes is cut by four, never accepted.
        Iterates are held above ``1e-10`` and below saturation; the counts of
        accepted steps whose final iterate needed either guard are kept in
        ``guard_accepted``.

        Returns
        -------
        list of tuple
            ``(t, Theta, N/rho at the faces)`` at each output time.
        """
        Th0, Thd = self.boundary_loadings(p0, pd)
        gbc = self._gbc(Th0, Thd)
        eye = sp.eye_array(self.n_z * self.n, format="csr")
        Th = np.tile(Thd, (self.n_z, 1))
        t = 0.0
        rec = []
        self.nsteps = self.maxinner = self.cutbacks = self.guard_accepted = 0
        for target in np.atleast_1d(np.asarray(t_out, float)):
            while t < target - 1e-15:
                step = min(max(dt0, dt_frac * max(t, dt0)), target - t)
                Th_old = Th.copy()
                while True:
                    Th_it = Th_old.copy()
                    ok = False
                    for k in range(1, maxit + 1):
                        Mf = self._M(Th_it, Th0, Thd)
                        A = (eye / step - self.div @ Mf @ self.grad).tocsc()
                        b = Th_old.reshape(-1) / step + (self.div @ Mf) @ gbc
                        Tn = spsolve(A, b).reshape(self.n_z, self.n)
                        guard = bool(np.any(Tn < 1e-10))
                        Tn = np.maximum(Tn, 1e-10)
                        fs = np.sum(Tn / self.Th_sat, -1)
                        bd = fs > 1 - 1e-9
                        if np.any(bd):
                            guard = True
                            Tn[bd] *= (1 - 1e-9) / fs[bd, None]
                        d = np.max(np.abs(Tn - Th_it)) / max(1.0, np.max(Th_it))
                        Th_it = Tn
                        if d < tol:
                            ok = True
                            break
                    if ok:
                        break
                    step *= 0.25
                    self.cutbacks += 1
                    if step < 1e-14 * max(1.0, t):
                        raise RuntimeError("inner iteration failed")
                Th = Th_it
                self.guard_accepted += int(guard)
                self.maxinner = max(self.maxinner, k)
                self.nsteps += 1
                t += step
            rec.append((t, Th.copy(), self.flux(Th, Th0, Thd)))
        return rec


def selectivity(J, p0):
    """Permeation selectivity of eq. (60), ``(N_2/N_1)/(p_2,0/p_1,0)``, from face fluxes."""
    Jm = J.mean(0)
    return float((Jm[1] / Jm[0]) / (p0[1] / p0[0]))


def locate_peak(t, y):
    """Maximum of a sampled curve: a parabola in ``ln t`` through the largest sample
    and its two neighbours, vertex kept inside that bracket.

    Returns
    -------
    t_peak, y_peak : float
    i : int
        Index of the largest sample.
    """
    i = int(np.argmax(y))
    if i == 0 or i == len(y) - 1:
        raise ValueError("the largest sample is at an end; the maximum is not bracketed")
    lt = np.log(t[i - 1:i + 2])
    c = np.polyfit(lt, y[i - 1:i + 2], 2)
    tp = float(np.clip(np.exp(-c[1] / (2 * c[0])), t[i - 1], t[i + 1]))
    return tp, float(np.polyval(c, np.log(tp))), i
