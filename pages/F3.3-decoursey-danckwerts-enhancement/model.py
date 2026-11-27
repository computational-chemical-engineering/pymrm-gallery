"""DeCoursey's enhancement factor and the exact surface-renewal model (page F3.3).

Gas A dissolves into an element of liquid that holds a non-volatile reactant B,
``A + zB -> products``, first order in each, for the age ``t`` of the element.
In ``theta = k2 B0 t``, ``zeta = x sqrt(k2 B0/D_A)``, ``a = A/A_i``,
``b = B/B0``, ``q = B0/(z A_i)`` and ``r = D_B/D_A``::

    a_theta = a_zetazeta - a b,      b_theta = r b_zetazeta - a b / q,
    a(0) = 1, b_zeta(0) = 0,  a(inf) = 0, b(inf) = 1,  a = 0, b = 1 at theta = 0.

The flux is ``F(theta) = -a_zeta(0, theta)``. Averaged over the Danckwerts
age distribution and over the Higbie exposure time ``theta_e = 4M/pi``::

    E_D(M) = M^(-1/2) int_0^inf F e^(-theta/M) dtheta,
    E_H(M) = sqrt(pi)/(2 sqrt(theta_e)) int_0^theta_e F dtheta.

The march works in ``eta = zeta/(2 sqrt(theta))`` and ``sigma = ln(theta)`` on a
fixed domain ``[0, eta_max]``::

    a_sigma = a_etaeta/4 + (eta/2) a_eta - theta a b,
    b_sigma = (r/4) b_etaeta + (eta/2) b_eta - (theta/q) a b,

so one march gives ``F theta = -a_eta(0) sqrt(theta)/2`` at every age, and
hence both averages at every ``M``.

Contents:

- the closed forms of DeCoursey (1974): :func:`e21` (Eq. 21), :func:`e20_residual`,
  :func:`e22_vkh` (Eq. 22), :func:`e23_hikita_asai` (Eq. 23, misprint M1 as an
  argument), the pseudo-first-order averages, and the instantaneous asymptote
  :func:`q_from_ea`;
- :class:`RenewalMarch`, the exact penetration problem by finite volumes
  (pymrm operators, ``pymrm.newton`` with ``NumJac`` per BDF2 step in sigma);
- :class:`History`, the flux and interface history of one march, with the
  Danckwerts and Higbie averages and their derivatives in ``M``;
- :func:`unique_root` and :func:`extremum_log`, the root and extremum helpers.
"""
import math

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import erf, erfc, erfcx, erfinv
from scipy.sparse import identity
from pymrm import (NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, construct_grad, newton)
from pymrm.checks import find_roots


# ------------------------------------------------------------- root helpers
def unique_root(f, lo, hi, n_scan=16, log=False, xtol=1e-13):
    """The single root of ``f`` on ``[lo, hi]``.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        The range, scanned in full by ``pymrm.checks.find_roots``.
    n_scan : int
        Scan points.
    log : bool
        Scan on a logarithmic grid.
    xtol : float
        Absolute tolerance of the refined root.

    Returns
    -------
    float

    Raises
    ------
    RuntimeError
        Unless the scan finds exactly one root and evaluates ``f`` everywhere.
    """
    found = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if found["failed"] or len(found["roots"]) != 1:
        raise RuntimeError(f"expected one root on [{lo:g}, {hi:g}]: {found['message']}")
    return float(found["roots"][0])


def extremum_log(f, df, lo, hi, n_scan=41, sign=1.0):
    """The largest ``sign * f`` on ``[lo, hi]``, refined from a log scan.

    A log-spaced scan of the whole window brackets the largest ``sign * f``;
    the extremum is the single root of the analytic derivative ``df`` between
    the scan neighbours. A scan maximum on a window end is returned as that
    end, flagged not interior.

    Parameters
    ----------
    f, df : callable
        The function of ``M`` and its derivative ``d f / d ln M``.
    lo, hi : float
        The window.
    n_scan : int
        Scan points.
    sign : float
        +1 for the maximum, -1 for the minimum.

    Returns
    -------
    tuple
        ``(M, f(M), interior)``.
    """
    xs = np.geomspace(lo, hi, n_scan)
    fs = np.array([sign * f(x) for x in xs])
    i = int(np.argmax(fs))
    if i in (0, n_scan - 1):
        return float(xs[i]), float(f(xs[i])), False
    lx = unique_root(lambda u: df(math.exp(u)), math.log(xs[i - 1]), math.log(xs[i + 1]), n_scan=8)
    return math.exp(lx), float(f(math.exp(lx))), True


# ------------------------------------------------------ DeCoursey's closed forms
def e21(M, Ea, q_power=2):
    """DeCoursey's Eq. (21), the positive root of Eq. (20).

    ``E = -M/(2(Ea-1)) + sqrt(M^2/(4(Ea-1)^2) + Ea M/(Ea-1) + 1)``.

    Parameters
    ----------
    M : float or ndarray
        ``D_A k2 B0 / k_L*^2`` (``k2 B0/s`` for Danckwerts).
    Ea : float
        Instantaneous asymptote of the enhancement factor.
    q_power : int
        Power of ``Ea - 1`` under ``M^2/4`` (2 as printed; 1 is a mistranscription).

    Returns
    -------
    float or ndarray
    """
    M = np.asarray(M, dtype=float)
    return -M / (2 * (Ea - 1)) + np.sqrt(M**2 / (4 * (Ea - 1) ** q_power) + Ea * M / (Ea - 1) + 1)


def de21_dlnM(M, Ea):
    """``d E21 / d ln M``, differentiated in closed form."""
    s = np.sqrt(M**2 / (4 * (Ea - 1) ** 2) + Ea * M / (Ea - 1) + 1)
    return M * (-1 / (2 * (Ea - 1)) + (M / (2 * (Ea - 1) ** 2) + Ea / (Ea - 1)) / (2 * s))


def e20_residual(E, M, Ea):
    """Residual of Eq. (20): ``E - sqrt(1 + M (Ea - E)/(Ea - 1))``."""
    return E - math.sqrt(1 + M * (Ea - E) / (Ea - 1))


def e22_vkh(M, Ea):
    """Van Krevelen-Hoftijzer, Eq. (22): ``E = sqrt(M f)/tanh sqrt(M f)``, ``f = (Ea-E)/(Ea-1)``."""
    def g(E):
        s = math.sqrt(max(M * (Ea - E) / (Ea - 1), 0.0))
        return E - (s / math.tanh(s) if s > 1e-12 else 1.0)
    return unique_root(g, 1.0, Ea)


def e23_hikita_asai(M, Ea, fn="erf"):
    """Hikita-Asai, Eq. (23), solved for ``E`` on ``[1, Ea]`` (or wider for ``erfc``).

    ``E = (sqrt(Mf) + (pi/8)/sqrt(Mf)) fn(sqrt(4 Mf/pi)) + exp(-4 Mf/pi)/2``.

    Parameters
    ----------
    M, Ea : float
    fn : {"erf", "erfc"}
        ``"erf"`` as used (misprint M1 corrected); ``"erfc"`` as printed.

    Returns
    -------
    float
    """
    phi = erf if fn == "erf" else erfc

    def g(E):
        m = M * (Ea - E) / (Ea - 1)
        if m <= 1e-300:
            return E - 0.5 - (0.0 if fn == "erf" else np.inf)
        return E - ((math.sqrt(m) + math.pi / 8 / math.sqrt(m)) * phi(math.sqrt(4 * m / math.pi))
                    + 0.5 * math.exp(-4 * m / math.pi))
    return unique_root(g, 1.0, Ea * (1 - 1e-15), n_scan=32)


def e_danckwerts_pfo(M):
    """Danckwerts pseudo-first-order enhancement, ``sqrt(1 + M)``."""
    return np.sqrt(1 + np.asarray(M, dtype=float))


def e_higbie_pfo(M):
    """Higbie pseudo-first-order enhancement (Eq. 23 with f = 1 and erf)."""
    M = np.asarray(M, dtype=float)
    return (np.sqrt(M) + np.pi / (8 * np.sqrt(M))) * erf(np.sqrt(4 * M / np.pi)) + 0.5 * np.exp(-4 * M / np.pi)


def q_from_ea(Ea, r, rule="penetration"):
    """``q = B0/(z A_i)`` that gives the instantaneous asymptote ``Ea``.

    The penetration relation (exact for Higbie and Danckwerts): with
    ``erf(lam) = 1/Ea``, ``q = e^(-lam^2) erfc(lam/sqrt r) / (sqrt(r) e^(-lam^2/r) erf(lam))``.
    At ``r = 1`` it is ``Ea - 1``, Eq. (14).

    Parameters
    ----------
    Ea, r : float
    rule : {"penetration", "eq14"}
        ``"eq14"`` returns ``Ea - 1`` at every ``r`` (misprint M2 as a defect).

    Returns
    -------
    float
    """
    if rule == "eq14":
        return Ea - 1.0
    lam = float(erfinv(1.0 / Ea))
    sr = math.sqrt(r)
    return math.exp(-lam**2) * erfcx(lam / sr) / (sr * erf(lam))   # erfc(x) = erfcx(x) e^(-x^2)


def ea_from_q(q, r):
    """The inverse of :func:`q_from_ea`: ``Ea = 1/erf(lam)`` with ``lam`` the root of the relation."""
    sr = math.sqrt(r)
    lam = unique_root(lambda L: math.exp(-L**2) * erfcx(L / sr) - sr * q * erf(L),
                      1e-12, 10.0, n_scan=64)
    return 1.0 / erf(lam)


# ------------------------------------------------------- the exact surface renewal
def sigma_steps(theta0, theta_end, dsigma, early_theta=1e-3, early_factor=4.0):
    """The step times in ``sigma = ln(theta)``.

    Steps of ``early_factor * dsigma`` up to ``early_theta`` (the reaction-free
    phase, steady in similarity variables) and ``dsigma`` after it; each part
    is fitted to a whole number of equal steps, so halving ``dsigma`` halves
    every step.

    Returns
    -------
    ndarray
        ``sigma_0 .. sigma_end``.
    """
    s0, s1, s2 = math.log(theta0), math.log(max(early_theta, theta0)), math.log(theta_end)
    n_e = max(int(math.ceil((s1 - s0) / (early_factor * dsigma) - 1e-9)), 0)
    n_l = int(math.ceil((s2 - s1) / dsigma - 1e-9))
    early = np.linspace(s0, s1, n_e + 1) if n_e else np.array([s0])
    return np.concatenate([early, np.linspace(s1, s2, n_l + 1)[1:]])


class RenewalMarch:
    """The exact penetration problem with a second-order reaction, marched in age.

    Finite volumes on a graded grid in ``eta`` (fields ``a``, ``b`` on the last
    axis), central fluxes for the ``(eta/2) d/d eta`` term in conservative form
    (flux ``-(eta/2) u`` plus the source ``-u/2``), and BDF2 steps in
    ``sigma = ln(theta)``, each solved with ``pymrm.newton``.

    Parameters
    ----------
    q : float
        ``B0/(z A_i)``; ``numpy.inf`` freezes ``b = 1`` (pseudo-first order).
    r : float
        ``D_B/D_A``.
    n : int
        Cells in ``eta``.
    dsigma : float
        Step in ``sigma`` after the reaction-free phase.
    theta0, theta_end : float
        First and last age of the march.
    eta_max : float, optional
        Domain end; default ``6 max(1, sqrt r)``.
    b_wall : {"neumann", "dirichlet"}
        ``b_eta = 0`` at the interface (B non-volatile), or ``b = 1`` there (a defect).
    stretch : float
        Faces at ``eta_max (e^(stretch xi) - 1)/(e^stretch - 1)``, ``xi`` uniform:
        cells cluster at the interface, where the reaction layer and the
        reaction plane sit (0 gives a uniform grid).
    conv : {"central", "upwind"}
        Face value of the ``eta/2`` flux: the mean of the two cells, or the
        upwind cell (first order, a defect).
    scheme : {"bdf2", "euler"}
        Step in ``sigma``: variable-step BDF2, or backward Euler (a defect).
    """

    def __init__(self, q, r, n=800, dsigma=0.02, theta0=1e-7, theta_end=5e6, eta_max=None,
                 b_wall="neumann", stretch=4.0, conv="central", scheme="bdf2"):
        self.q, self.r, self.n = float(q), float(r), int(n)
        self.theta0, self.theta_end, self.dsigma = float(theta0), float(theta_end), float(dsigma)
        self.eta_max = float(eta_max if eta_max is not None else 6.0 * max(1.0, math.sqrt(r)))
        self.conv, self.scheme = conv, scheme
        self.shape = (self.n, 2)                       # cells x fields (a, b)
        xi = np.linspace(0.0, 1.0, self.n + 1)          # the map is fixed, so refinement is nested
        self.x_f = (self.eta_max * np.expm1(stretch * xi) / np.expm1(stretch) if stretch > 0
                    else self.eta_max * xi)
        self.x_c = 0.5 * (self.x_f[1:] + self.x_f[:-1])
        # eta = 0, outward normal -eta:  a = 1 (Dirichlet);  -b_eta = 0 (B non-volatile)
        left = ({"a": [0.0, 1.0], "b": [1.0, 0.0], "d": [1.0, 0.0]} if b_wall == "neumann"
                else {"a": [0.0, 0.0], "b": [1.0, 1.0], "d": [1.0, 1.0]})   # defect: b = 1
        # eta = eta_max, outward normal +eta:  a = 0, b = 1 (undisturbed liquid)
        right = {"a": [0.0, 0.0], "b": [1.0, 1.0], "d": [0.0, 1.0]}
        self.bc = (left, right)
        self._build_operators()

    def _build_operators(self):
        """Assemble the constant transport operator once: diffusion plus the eta/2 term."""
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc)
        diff = construct_coefficient_matrix(np.array([[0.25, 0.25 * self.r]]), shape=(self.n + 1, 2))
        w = (0.5 * self.x_f).reshape((-1, 1))          # eta/2 at the faces
        up, up_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c, self.bc, v=w)
        dn, dn_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c, self.bc, v=-w)
        # flux J = -(1/4, r/4) u_eta - (eta/2) u_face, u_face the central (mean) face value
        if self.conv == "central":
            conv, conv_bc = 0.5 * (dn - up), 0.5 * (dn_bc - up_bc)
        else:                                          # the flow -eta/2 points to the interface
            conv, conv_bc = dn, dn_bc
        flux = -diff @ grad + conv
        flux_bc = -diff @ grad_bc + conv_bc
        div = construct_div(self.shape, self.x_f, nu=0)   # nu=0: planar, the penetration layer
        eye = identity(2 * self.n, format="csc")
        self.lin = -div @ flux - 0.5 * eye             # u_sigma = -dJ/deta - u/2 - reaction
        self.g = -(div @ flux_bc).toarray().ravel()
        self.eye = eye
        self.numjac = NumJac(self.shape)

    def _source(self, theta):
        k_b = 0.0 if not np.isfinite(self.q) else theta / self.q
        def src(u):
            ab = u[..., 0] * u[..., 1]
            return np.stack([-theta * ab, -k_b * ab], axis=-1)
        return src

    def _read(self, u, theta):
        """``F theta`` and ``b(0)`` from the interface face under the same bc."""
        val, grad, _, _ = compute_boundary_values(u, self.x_f, self.x_c, self.bc)
        val, grad = np.ravel(val), np.ravel(grad)
        return -grad[0] * math.sqrt(theta) / 2.0, val[1], -grad[0] / 2.0

    def run(self, tol=1e-10):
        """March from ``theta0`` to ``theta_end``.

        Returns
        -------
        History

        Raises
        ------
        RuntimeError
            If a Newton solve does not report success.
        """
        from scipy.sparse.linalg import spsolve
        sig = sigma_steps(self.theta0, self.theta_end, self.dsigma)
        u = spsolve(self.lin, -self.g).reshape(self.shape)   # reaction-free state at theta0
        ft, b0, c0 = self._read(u, math.exp(sig[0]))
        F, B, u_prev, h_prev = [ft], [b0], None, None
        for k in range(1, sig.size):
            h, theta = sig[k] - sig[k - 1], math.exp(sig[k])
            if u_prev is None or self.scheme == "euler":
                c = (1.0, 1.0, 0.0)                    # backward Euler, first step
            else:
                om = h / h_prev                         # variable-step BDF2
                c = ((1 + 2 * om) / (1 + om), 1 + om, om**2 / (1 + om))
            hist = (c[1] * u - (c[2] * u_prev if c[2] else 0.0)).ravel()
            src = self._source(theta)

            def fun(x, hist=hist, c0_=c[0], h=h, src=src):
                gs, js = self.numjac(src, x.reshape(self.shape))
                res = (c0_ * x.ravel() - hist) / h - (self.lin @ x.ravel() + self.g + gs.ravel())
                return res, (c0_ / h) * self.eye - self.lin - js

            guess = u if u_prev is None else (u + (u - u_prev) * (h / h_prev))
            sol = newton(fun, np.clip(guess, 0.0, 1.0), tol=tol, maxfev=50, stagnation=True)
            if not sol.success:
                raise RuntimeError(f"newton failed at theta = {theta:.3g}: {sol.message}")
            u_prev, u, h_prev = u, sol.x.reshape(self.shape), h
            ft, b0, _ = self._read(u, theta)
            F.append(ft)
            B.append(b0)
        return History(sig, np.array(F), np.array(B), c0, self.theta0, u, self)


class History:
    """The flux and interface history of one march, and its age averages.

    Parameters
    ----------
    sigma, ftheta, b0 : ndarray
        Step times ``ln(theta)``, ``F theta`` and ``b(0)`` at each step.
    c_start : float
        ``F sqrt(theta)`` of the reaction-free start, so ``int_0^theta0 F = 2 c sqrt(theta0)``.
    theta0 : float
    u_end : ndarray
        The last profile.
    march : RenewalMarch
    n_gauss : int
        Gauss-Legendre points per march interval for the Danckwerts integrals.
    """

    def __init__(self, sigma, ftheta, b0, c_start, theta0, u_end, march, n_gauss=4):
        self.sigma, self.ftheta, self.b0 = sigma, ftheta, b0
        self.c, self.theta0, self.u_end, self.march = c_start, theta0, u_end, march
        self.theta_end = math.exp(sigma[-1])
        self.S = CubicSpline(sigma, ftheta)
        self.Sb = CubicSpline(sigma, b0)
        self.S_int = self.S.antiderivative()
        self.set_quadrature(n_gauss)

    def set_quadrature(self, n_gauss, sub=1):
        """Gauss-Legendre nodes: ``n_gauss`` per subinterval, ``sub`` subintervals per step."""
        self.n_gauss = n_gauss
        x, w = np.polynomial.legendre.leggauss(n_gauss)
        edges = np.interp(np.arange(0, (self.sigma.size - 1) * sub + 1) / sub,
                          np.arange(self.sigma.size), self.sigma)
        mid, half = 0.5 * (edges[1:] + edges[:-1]), 0.5 * np.diff(edges)
        self.sq = (mid[:, None] + half[:, None] * x[None, :]).ravel()
        self.wq = (half[:, None] * w[None, :]).ravel()
        self.Fq, self.Bq, self.tq = self.S(self.sq), self.Sb(self.sq), np.exp(self.sq)

    def _check_M(self, M, factor=30.0):
        if not (M * factor <= self.theta_end):
            raise ValueError(f"M = {M:g}: the march ends at theta = {self.theta_end:g}, "
                             f"below {factor:g} M")

    def danckwerts(self, M, weight="renewal"):
        """``E_D(M)``; ``weight="unit"`` uses ``e^(-theta)`` (the 1/M dropped, a defect)."""
        self._check_M(M)
        k = 1.0 / M if weight == "renewal" else 1.0
        start = self.c * math.sqrt(math.pi / k) * erf(math.sqrt(self.theta0 * k))
        return (start + float(np.sum(self.wq * self.Fq * np.exp(-self.tq * k)))) / math.sqrt(M)

    def danckwerts_dlnM(self, M):
        """``d E_D / d ln M`` by differentiating under the integral."""
        self._check_M(M)
        e = np.exp(-self.tq / M)
        J = self.c * math.sqrt(math.pi * M) * erf(math.sqrt(self.theta0 / M)) + np.sum(self.wq * self.Fq * e)
        dJ = (self.c * math.sqrt(math.pi) / (2 * math.sqrt(M)) * erf(math.sqrt(self.theta0 / M))
              - self.c * math.sqrt(self.theta0) * math.exp(-self.theta0 / M) / M
              + np.sum(self.wq * self.Fq * e * self.tq) / M**2)
        return M * (-J / (2 * M**1.5) + dJ / math.sqrt(M))

    def mean_b0(self, M):
        """Danckwerts mean of the interfacial ``b``: ``(1/M) int b(0) e^(-theta/M) dtheta``."""
        self._check_M(M)
        start = 1.0 - math.exp(-self.theta0 / M)
        return start + float(np.sum(self.wq * self.Bq * self.tq * np.exp(-self.tq / M))) / M

    def _flux_integral(self, theta_e):
        if not (self.theta0 <= theta_e <= self.theta_end):
            raise ValueError(f"theta_e = {theta_e:g} outside the march")
        return 2 * self.c * math.sqrt(self.theta0) + float(self.S_int(math.log(theta_e)) - self.S_int(self.sigma[0]))

    def higbie(self, M, exposure="4M/pi"):
        """``E_H(M)`` with ``theta_e = 4M/pi`` (Brian et al.'s ``M``); ``"M"`` is a defect."""
        te = 4 * M / math.pi if exposure == "4M/pi" else M
        return math.sqrt(math.pi) / (2 * math.sqrt(te)) * self._flux_integral(te)

    def higbie_dlnM(self, M, exposure="4M/pi"):
        """``d E_H / d ln M``: the flux at ``theta_e`` minus the dilution of the mean."""
        te = 4 * M / math.pi if exposure == "4M/pi" else M
        I, F = self._flux_integral(te), float(self.S(math.log(te))) / te
        return te * math.sqrt(math.pi) / 2 * (F / math.sqrt(te) - I / (2 * te**1.5))
