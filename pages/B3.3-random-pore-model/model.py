"""The random pore model of Bhatia and Perlmutter (1980, Part I) and its rivals.

Model code of gallery page B3.3 (kinetic control). Dimensionless groups of the
paper's NOTATION list: ``psi = 4 pi L_o (1 - eps_o) / S_o**2``,
``sigma = R_o S_o / (1 - eps_o)``, ``tau = k_s C**n S_o t / (1 - eps_o)``.
Lengths in cm.

Contents
--------
psi_of, X_of_tau, SS0_of_X, dSS0_dX, rate37, XM_eq33, X_grain
    The random pore model: eqs. (32), (31), its derivative, (37), (33) and the
    grain model (39)/(40).
pet_S, pet_eps, u_of_eps, G_of_eps, cubic7_roots, X_petersen, X_pet_via_eps,
SS0_petersen, SS0_pet_via_eq6, petersen_audit, curve_D_length
    Petersen (1957) as the paper reprints it, eqs. (3)-(7).
Fig7
    Curves A (random pore) and B (Petersen) of Figure 7 and the comparison of
    the two: crossings, extrema and parting conversion, all root-found.
tau_at_X_block, eps_o_opt_block, march37
    The pymrm parts: blocks of independent scalar inversions with
    ``NumJac((K, 1))`` and ``newton``, and the rate law (37) marched in time.
J_gam, I_closed, I_energy, I_nrms, I_gl, dI_dm, m_star, m_star_brent, psi_m_star_zero
    The paper's squared-deviation integral I between eq. (37) and the grain
    model, in closed form and by Gauss-Legendre quadrature, and its minimiser.
S_of_eps_o, stat35
    Eq. (27) at fixed final porosity and the stationarity condition (35).
overlap_gap, overlap_gap_slope
    Avrami's overlap law against Petersen's, as functions of u.
yagi_kunii_eq6, eta_sphere
    Eq. (6) of Yagi and Kunii (page B3.1) and the sphere effectiveness factor.
unique_root
    ``pymrm.checks.find_roots`` that raises unless exactly one root is found.

Several functions take an argument whose default is the printed model and
whose alternative is a defect the page's break table injects.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar
from scipy.special import gamma as gamma_fn
from scipy.special import expit, gammaincc

from pymrm import NumJac, newton
from pymrm.checks import find_roots


# ------------------------------------------------------------------ helpers
def unique_root(f, lo, hi, n_scan=64, log=False, xtol=1e-15):
    """The single root of ``f`` on ``[lo, hi]``, found by a scan and refined.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        Search interval.
    n_scan : int
        Scan points of ``pymrm.checks.find_roots``.
    log : bool
        Logarithmic scan.
    xtol : float
        Absolute tolerance of the refined root.

    Returns
    -------
    float
        The root.

    Raises
    ------
    RuntimeError
        When the scan finds no root or more than one.
    """
    res = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(res["roots"]) != 1:
        raise RuntimeError(f"expected one root on [{lo:g}, {hi:g}], found {res['roots']}: "
                           f"{res['message']}")
    return float(res["roots"][0])


# --------------------------------------------------- the random pore model
def psi_of(So, Lo, eps, coef=4 * np.pi, keep_porosity=True):
    """Pore structure parameter of the NOTATION list (p. 386).

    Parameters
    ----------
    So, Lo, eps : float
        Initial surface (cm2/cm3), pore length (cm/cm3) and porosity.
    coef : float
        ``4 pi`` as printed; another value is a defect.
    keep_porosity : bool
        ``False`` drops the factor ``(1 - eps_o)`` (a defect).

    Returns
    -------
    float
        ``psi = 4 pi L_o (1 - eps_o) / S_o**2``.
    """
    return coef * Lo * ((1 - eps) if keep_porosity else 1.0) / So**2


def X_of_tau(tau, psi, sigma=np.inf, quarter=4.0):
    """Conversion, eq. (32).

    Parameters
    ----------
    tau : array_like
        Dimensionless time.
    psi : float or array_like
        Pore structure parameter.
    sigma : float
        Particle size parameter; ``inf`` removes the receding external surface.
    quarter : float
        The 4 of ``psi tau / 4``; another value is a defect.

    Returns
    -------
    numpy.ndarray
        ``X = 1 - (1 - tau/sigma)**3 exp[-tau (1 + psi tau / 4)]``.
    """
    tau = np.asarray(tau, float)
    f = 1.0 if not np.isfinite(sigma) else (1 - tau / sigma) ** 3
    return 1 - f * np.exp(-tau * (1 + psi * tau / quarter))


def SS0_of_X(X, psi):
    """Reaction surface over its initial value, eq. (31) at ``sigma -> inf``."""
    w = 1 - np.asarray(X, float)
    return w * np.sqrt(1 - psi * np.log(w))


def dSS0_dX(X, psi):
    """Derivative of eq. (31) with respect to ``X``, written out by hand.

    With ``w = 1 - X`` and ``r = 1 - psi ln w``:
    ``dS/dX = -sqrt(r) + psi / (2 sqrt(r))``.
    """
    r = 1 - psi * np.log(1 - np.asarray(X, float))
    return -np.sqrt(r) + psi / (2 * np.sqrt(r))


def rate37(X, psi):
    """Rate law eq. (37), ``dX/dtau = (1 - X) sqrt(1 - psi ln(1 - X))``."""
    w = np.maximum(1 - np.asarray(X, float), 1e-300)
    return w * np.sqrt(1 - psi * np.log(w))


def XM_eq33(psi, denom=2.0):
    """Conversion at the surface maximum, eq. (33).

    ``denom = 2`` is the printed ``(2 - psi)/(2 psi)``; ``denom = 1`` is the
    misgrouped ``(2 - psi)/psi`` (a defect).
    """
    return 1 - np.exp((2 - psi) / (denom * psi))


def X_grain(tau, m):
    """Grain model, eq. (39) integrated; ``m = 2/3`` is eq. (40), ``m = 1`` its branch."""
    tau = np.asarray(tau, float)
    if m == 1:
        return 1 - np.exp(-tau)
    return 1 - (1 - (1 - m) * tau) ** (1 / (1 - m))


# --------------------------------------------- Petersen (1957), as reprinted
def pet_S(r, L):
    """Petersen's surface, eq. (3): ``2 pi r L (1 - r sqrt(pi L / 3))``."""
    return 2 * np.pi * r * L * (1 - r * np.sqrt(np.pi * L / 3))


def pet_eps(r, L, c=2 / 3):
    """Petersen's porosity, eq. (4): ``pi r**2 L (1 - c r sqrt(pi L / 3))``.

    ``c = 2/3`` as printed; ``c = 1/3`` is the misread eq. (4) (a defect).
    """
    return np.pi * r**2 * L * (1 - c * r * np.sqrt(np.pi * L / 3))


def u_of_eps(eps, c=2 / 3):
    """Eq. (4) in ``u = r sqrt(pi L / 3)``: ``eps = 3 u**2 - 3 c u**3``, for ``u`` in (0, 2/(3c)]."""
    return unique_root(lambda u: 3 * u**2 - 3 * c * u**3 - eps, 0.0, 2 / (3 * c), n_scan=16)


def G_of_eps(eps):
    """The root of eq. (7) that eqs. (3), (4) and (6) select: ``G = 3/(2 u_o)``."""
    return 1.5 / u_of_eps(eps)


def cubic7_roots(eps, coef=4 / 27):
    """The three roots of eq. (7), ``(4/27) eps G**3 - G + 1 = 0``, sorted.

    ``coef = 4/27`` as printed; ``2/27`` is a misread (a defect).
    """
    r = np.roots([coef * eps, 0.0, -1.0, 1.0])
    return np.sort(r[np.isreal(r)].real)


def X_petersen(s, eps, G=None):
    """Petersen conversion, eq. (6), with ``s = k_s C**n t / r_po``.

    ``G`` defaults to the selected root of eq. (7); another root is a defect.
    """
    G = G_of_eps(eps) if G is None else G
    s = np.asarray(s, float)
    return eps / (1 - eps) * ((1 + s) ** 2 * (G - 1 - s) / (G - 1) - 1)


def dX_petersen_ds(s, eps, G=None):
    """Derivative of eq. (6) with respect to ``s``."""
    G = G_of_eps(eps) if G is None else G
    s = np.asarray(s, float)
    return eps / (1 - eps) * (2 * (1 + s) * (G - 1 - s) - (1 + s) ** 2) / (G - 1)


def X_pet_via_eps(s, eps):
    """Petersen conversion by the porosity route: ``u = u_o (1 + s)``, eq. (4) in ``u``."""
    u = u_of_eps(eps) * (1 + np.asarray(s, float))
    return (3 * u**2 - 2 * u**3 - eps) / (1 - eps)


def SS0_petersen(X, eps):
    """Petersen ``S/S_o`` along the conversion history, from eq. (3): ``S`` is proportional to ``u (1 - u)``."""
    u0 = u_of_eps(eps)
    uu = np.array([u_of_eps(eps + x * (1 - eps)) for x in np.atleast_1d(X)])
    return (uu * (1 - uu) / (u0 * (1 - u0))).reshape(np.shape(X))


def SS0_pet_via_eq6(X, eps, G=None):
    """Petersen ``S/S_o`` by differentiating eq. (6) in ``s`` (normalised at ``s = 0``)."""
    G = G_of_eps(eps) if G is None else G
    s = u_of_eps(eps + X * (1 - eps)) / u_of_eps(eps) - 1
    return (2 * (1 + s) * (G - 1 - s) - (1 + s) ** 2) / (2 * G - 3)


def petersen_audit(So, Lo, eps, root="physical", c=2 / 3):
    """Check a printed triple ``(S_o, L_o, eps_o)`` against eqs. (3) and (4), both ways.

    Parameters
    ----------
    So, Lo, eps : float
        The printed triple.
    root : {"physical", "upper"}
        Which root of eq. (3) in ``r_p``: below the surface maximum (physical)
        or above it (a defect).
    c : float
        Coefficient of eq. (4), see :func:`pet_eps`.

    Returns
    -------
    dict
        ``Smax`` (largest surface at length ``Lo``), ``feasible``, ``rp_from_S``,
        ``eps_pred`` and ``eps_rel`` (eq. 3 then eq. 4), ``S_pred`` and ``S_rel``
        (eq. 4 then eq. 3), ``psi``.
    """
    b = np.sqrt(np.pi * Lo / 3)
    rmax = 1 / (2 * b)
    out = dict(So=So, Lo=Lo, eps=eps, Smax=pet_S(rmax, Lo), psi=psi_of(So, Lo, eps))
    out["feasible"] = So <= out["Smax"]
    if out["feasible"]:
        lo, hi = (1e-12 * rmax, rmax) if root == "physical" else (rmax, 1 / b)
        rp = unique_root(lambda z: pet_S(z, Lo) / So - 1, lo, hi, n_scan=32, xtol=1e-20)
        out.update(rp_from_S=rp, eps_pred=pet_eps(rp, Lo, c))
        out["eps_rel"] = abs(out["eps_pred"] - eps) / eps
    else:
        out.update(rp_from_S=np.nan, eps_pred=np.nan, eps_rel=np.nan)
    r2 = u_of_eps(eps, c) / b
    out.update(rp_from_eps=r2, S_pred=pet_S(r2, Lo))
    out["S_rel"] = abs(out["S_pred"] - So) / So
    return out


def curve_D_length(So, eps):
    """Pore length that makes eqs. (3) and (4) hold together at ``(S_o, eps_o)``.

    Returns
    -------
    tuple of float
        ``(L, r_p)`` in cm/cm3 and cm.
    """
    u0 = u_of_eps(eps)
    L = So**2 / (4 * u0**2 * (1 - u0) ** 2 * 3 * np.pi)
    return L, u0 / np.sqrt(np.pi * L / 3)


# ------------------------------------------------------------- Figure 7
class Fig7:
    """Curves A (random pore, eq. 32) and B (Petersen, eq. 6) on one time axis.

    The time axis is ``k_s C**n t`` (cm). Curve A uses ``psi`` from
    ``(S_o, L_o, eps_o)``; curve B is the Petersen structure fixed by one
    pairing of the over-determined triple.

    Parameters
    ----------
    So, Lo, eps : float
        Curves A and B's printed parameters.
    pairing : {"eps_L", "eps_S", "S_L"}
        Which two of the three fix curve B: ``"eps_L"`` is the paper's (p. 384).
    G : float, optional
        Root of eq. (7) for curve B; default the selected one (another is a defect).
    quarter : float
        Passed to :func:`X_of_tau` for curve A (4 as printed).
    """

    def __init__(self, So, Lo, eps, pairing="eps_L", G=None, quarter=4.0):
        self.So, self.Lo, self.eps_A, self.quarter = So, Lo, eps, quarter
        self.psi = psi_of(So, Lo, eps)
        if pairing == "eps_L":
            self.eps_B, self.L_B = eps, Lo
            self.rpo = u_of_eps(eps) / np.sqrt(np.pi * Lo / 3)
        elif pairing == "eps_S":
            self.eps_B = eps
            self.L_B, self.rpo = curve_D_length(So, eps)
        elif pairing == "S_L":
            aud = petersen_audit(So, Lo, eps)
            self.eps_B, self.L_B, self.rpo = aud["eps_pred"], Lo, aud["rp_from_S"]
        else:
            raise ValueError(pairing)
        self.G = G_of_eps(self.eps_B) if G is None else G
        self.kt_end = self.rpo * (2 * self.G / 3 - 1)   # complete gasification of curve B

    def XA(self, kt):
        """Curve A, eq. (32) at ``sigma -> inf``."""
        return X_of_tau(np.asarray(kt) * self.So / (1 - self.eps_A), self.psi, quarter=self.quarter)

    def XB(self, kt):
        """Curve B, eq. (6)."""
        return X_petersen(np.asarray(kt) / self.rpo, self.eps_B, self.G)

    def dAB(self, kt):
        """``X_A - X_B``."""
        return self.XA(kt) - self.XB(kt)

    def ddAB(self, kt):
        """Derivative of ``X_A - X_B`` in ``k_s C**n t``, analytic."""
        tau = np.asarray(kt) * self.So / (1 - self.eps_A)
        dA = (1 - self.XA(kt)) * (1 + 2 * self.psi * tau / self.quarter) * self.So / (1 - self.eps_A)
        dB = dX_petersen_ds(np.asarray(kt) / self.rpo, self.eps_B, self.G) / self.rpo
        return dA - dB

    def kt_at_XA(self, X):
        """Time at which curve A reaches ``X`` (closed-form inverse of eq. 32)."""
        tau = (-1 + np.sqrt(1 - self.psi / self.quarter * 4 * np.log(1 - X))) * self.quarter / (2 * self.psi)
        return tau * (1 - self.eps_A) / self.So

    def compare(self, X_upto=0.75, gap=0.02, n_scan=400):
        """Crossings, extrema up to ``X_A = X_upto`` and the first ``|X_A - X_B| = gap``.

        Every location is a root: crossings of ``X_A - X_B``, extrema as roots
        of its analytic derivative, the parting conversion as the first root of
        ``|X_A - X_B| - gap``, each on a logarithmic scan of ``n_scan`` points
        (``pymrm.checks.find_roots``).

        Returns
        -------
        dict
            ``crossings`` and ``extrema`` (lists of ``k_s C**n t``), ``kt_upto``,
            ``max_upto`` (largest ``|X_A - X_B|`` for ``X_A <= X_upto``) and
            ``X_part`` (``X_A`` where the curves first part by ``gap``).
        """
        hi = self.kt_at_XA(0.999999)
        if self.kt_end > 0:          # curve B ends at complete gasification
            hi = min(hi, self.kt_end)
        lo = 1e-12 * hi
        kt_up = self.kt_at_XA(X_upto)
        kw = dict(n_scan=n_scan, log=True, xtol=1e-24)
        cross = find_roots(self.dAB, lo, hi, **kw)["roots"]
        ext = find_roots(self.ddAB, lo, hi, **kw)["roots"]
        cand = [kt_up] + [k for k in ext if k <= kt_up]
        part = find_roots(lambda z: abs(self.dAB(z)) - gap, lo, hi, **kw)["roots"]
        return dict(crossings=list(cross), extrema=list(ext), kt_upto=kt_up,
                    max_upto=float(max(abs(self.dAB(k)) for k in cand)),
                    X_part=float(self.XA(part[0])) if part else np.nan)


# ------------------------------------------------------ pymrm: blocks, march
def tau_at_X_block(Xt, psis, sigmas, tol=1e-13, residual="log"):
    """Invert eq. (32) for ``tau`` over a block of ``(psi, sigma)`` at once.

    Each entry is an independent scalar equation, so the state has shape
    ``(K, 1)`` and ``NumJac((K, 1))`` gives a diagonal Jacobian. ``tau`` is
    carried as ``sigma / (1 + exp(-w))`` so the iterate cannot leave
    ``(0, sigma)``; Newton starts at ``w = 0``, ``tau = sigma / 2``.

    Parameters
    ----------
    Xt : float or array_like
        Target conversion.
    psis, sigmas : array_like
        Pore structure and particle size parameters (broadcast together).
    tol : float
        ``newton`` step tolerance on ``w``.
    residual : {"log", "direct"}
        ``"log"`` solves ``3 ln(1 - tau/sigma) - tau (1 + psi tau/4) - ln(1 - X) = 0``;
        ``"direct"`` solves ``X(tau) - X = 0``, whose exponential underflows
        at the start for large ``psi`` and ``sigma``.

    Returns
    -------
    numpy.ndarray
        ``tau`` for each entry.

    Raises
    ------
    RuntimeError
        When ``newton`` does not converge.
    """
    ps = np.atleast_1d(np.asarray(psis, float)).ravel()
    sg = np.broadcast_to(np.asarray(sigmas, float), ps.shape).ravel()
    xt = np.broadcast_to(np.asarray(Xt, float), ps.shape).ravel()
    numjac = NumJac((ps.size, 1))   # K independent scalars: (K, 1), never (K,)

    def res(w):
        g = expit(w[:, 0])   # tau / sigma in (0, 1)
        t = sg * g
        if residual == "log":
            return (3 * np.log1p(-g) - t * (1 + ps * t / 4) - np.log1p(-xt))[:, None]
        return (1 - (1 - g) ** 3 * np.exp(-t * (1 + ps * t / 4)) - xt)[:, None]

    sol = newton(lambda w: numjac(res, w), np.zeros((ps.size, 1)), tol=tol, maxfev=300)
    if not sol.success:
        raise RuntimeError(f"tau block: newton did not converge ({sol.message})")
    return sg * expit(sol.x[:, 0])


def eps_o_opt_block(eps_final, Ghat, tol=1e-14):
    """Solve eq. (35) for the initial porosity that maximises the final surface.

    ``eps_o = eps sigmoid(w)`` keeps ``0 < eps_o < eps``, the physical branch.

    Parameters
    ----------
    eps_final : float or array_like
        Final porosity.
    Ghat : float or array_like
        The group ``S_o**2 / (2 pi L_o)``.
    tol : float
        ``newton`` step tolerance on ``w``.

    Returns
    -------
    numpy.ndarray
        The optimal ``eps_o`` per entry.

    Raises
    ------
    RuntimeError
        When ``newton`` does not converge.
    """
    ef = np.atleast_1d(np.asarray(eps_final, float)).ravel()
    gh = np.broadcast_to(np.asarray(Ghat, float), ef.shape).ravel()
    numjac = NumJac((ef.size, 1))

    def res(w):
        e0 = ef * expit(w[:, 0])
        return stat35(e0, ef, gh)[:, None]

    sol = newton(lambda w: numjac(res, w), -2.0 * np.ones((ef.size, 1)), tol=tol, maxfev=300)
    if not sol.success:
        raise RuntimeError(f"eq. (35) block: newton did not converge ({sol.message})")
    return ef * expit(sol.x[:, 0])


def march37(psis, tau_end, nsteps, theta=0.5):
    """March the rate law eq. (37) from ``X = 0`` with the theta method.

    One ``newton`` solve per step for all ``psi`` at once, state ``(K, 1)``
    and ``NumJac((K, 1))``. ``theta = 0.5`` is Crank-Nicolson, ``1`` implicit
    Euler.

    Parameters
    ----------
    psis : array_like
        Pore structure parameters.
    tau_end : float
        Final time.
    nsteps : int
        Number of equal steps.
    theta : float
        Implicitness.

    Returns
    -------
    numpy.ndarray
        ``X(tau_end)`` per ``psi``.

    Raises
    ------
    RuntimeError
        When a step's ``newton`` solve does not converge.
    """
    ps = np.asarray(psis, float).ravel()
    dt = tau_end / nsteps
    X = np.zeros((ps.size, 1))
    numjac = NumJac((ps.size, 1))
    for _ in range(nsteps):
        X_old = X.copy()
        r_old = rate37(X_old[:, 0], ps)

        def res(Xn):
            return Xn - X_old - dt * (theta * rate37(Xn[:, 0], ps) + (1 - theta) * r_old)[:, None]

        sol = newton(lambda z: numjac(res, z), X.copy(), tol=1e-14, maxfev=100)
        if not sol.success:
            raise RuntimeError(f"march37: newton did not converge ({sol.message})")
        X = sol.x.reshape(X.shape)
    return X[:, 0]


# ------------------------------------------- the squared-deviation integral I
@lru_cache(maxsize=None)
def _gl(n):
    return leggauss(n)


def J_gam(b, psi):
    """``J(b) = int_0^1 w**b sqrt(1 - psi ln w) dw`` in closed form (incomplete gamma)."""
    if psi == 0:
        return 1.0 / (b + 1)
    c = (b + 1) / psi
    return np.exp(c) / psi * (psi / (b + 1)) ** 1.5 * gammaincc(1.5, c) * gamma_fn(1.5)


def I_energy(psi, p=1):
    """``int_0^1 [w**p sqrt(1 - psi ln w)]**2 dw``, the scale of I."""
    return 1 / (2 * p + 1) + psi / (2 * p + 1) ** 2


def I_closed(psi, m, p=1):
    """The paper's I between ``w**p sqrt(1 - psi ln w)`` and ``w**m``, closed form.

    ``p = 1`` is the printed model term ``(1 - X) sqrt(...)``; ``p = 0`` drops
    the factor ``(1 - X)`` (a defect).
    """
    return I_energy(psi, p) - 2 * J_gam(m + p, psi) + 1 / (2 * m + 1)


def I_nrms(psi, m, p=1):
    """Normalised rms misfit ``sqrt(I / I_energy)``."""
    return np.sqrt(max(I_closed(psi, m, p), 0.0) / I_energy(psi, p))


def I_gl(psi, m, n=400, p=1, v_end=30.0):
    """I by Gauss-Legendre quadrature in ``v = -ln w`` on ``[0, v_end]``."""
    x, wt = _gl(n)
    v = 0.5 * v_end * (x + 1)
    W = 0.5 * v_end * wt
    return float(np.sum(W * (np.exp(-p * v) * np.sqrt(1 + psi * v) - np.exp(-m * v)) ** 2 * np.exp(-v)))


def dI_dm(psi, m, n=600, p=1, X_max=None):
    """``dI/dm`` by Gauss-Legendre quadrature.

    Full range (``X_max=None``): ``-2 dJ/db(m + p) - 2/(2m + 1)**2`` with
    ``dJ/db`` by ``n`` nodes on ``v`` in ``[0, 80]``. ``X_max < 1`` truncates
    the integral to ``X <= X_max`` (a defect) and integrates every term.
    """
    x, wt = _gl(n)
    if X_max is None:
        v, W = 40 * (x + 1), 40 * wt
        dJ = float(np.sum(W * (-v) * np.exp(-(m + p + 1) * v) * np.sqrt(1 + psi * v)))
        return -2 * dJ - 2 / (2 * m + 1) ** 2
    v_end = -np.log(1 - X_max)
    v, W = 0.5 * v_end * (x + 1), 0.5 * v_end * wt
    e = np.exp(-v)
    return float(np.sum(W * 2 * (e**m - e**p * np.sqrt(1 + psi * v)) * e**m * (-v) * e))


@lru_cache(maxsize=None)
def m_star(psi, n=600, p=1, X_max=None):
    """Best-fit grain shape factor: the single root of ``dI/dm`` on ``(-0.45, 3)``."""
    return unique_root(lambda z: dI_dm(psi, z, n, p, X_max), -0.45, 3.0, n_scan=48)


def m_star_brent(psi):
    """Best-fit ``m`` by bounded minimisation of the closed-form I (a second route)."""
    lo = max(m_star(psi) - 0.25, -0.4)
    return minimize_scalar(lambda z: I_closed(psi, z), bounds=(lo, lo + 0.5),
                           method="bounded", options=dict(xatol=1e-13)).x


def psi_m_star_zero(X_max=None):
    """The ``psi`` at which the best-fit ``m`` is zero: the root of ``dI/dm(psi, 0)``."""
    return unique_root(lambda pv: dI_dm(pv, 0.0, X_max=X_max), 1.0, 100.0, n_scan=48, log=True,
                       xtol=1e-13)


# ----------------------------------------------------- optimal porosity (35)
def S_of_eps_o(e0, ef, Ghat):
    """Eq. (27) with ``V = eps`` and ``V_o = eps_o``: ``S/S_o`` at final porosity ``ef``."""
    return (1 - ef) / (1 - e0) * np.sqrt(1 - (2 * (1 - e0) / Ghat) * np.log((1 - ef) / (1 - e0)))


def stat35(e0, ef, Ghat):
    """Residual of eq. (35): ``Ghat/(1 - eps_o) - 1 - ln[(1 - eps)/(1 - eps_o)]``."""
    return Ghat / (1 - e0) - 1 - np.log((1 - ef) / (1 - e0))


# ----------------------------------------------------------- overlap laws
def overlap_gap(u, c=2.0, avrami=1.0):
    """Avrami's ``1 - exp(-avrami V_E)`` minus Petersen's ``3u**2 - c u**3``, with ``V_E = 3u**2``.

    ``c = 2`` is eq. (4) as printed (``c = 1`` the misread); ``avrami = 1`` is eq. (20).
    """
    u = np.asarray(u, float)
    return (1 - np.exp(-avrami * 3 * u**2)) - (3 * u**2 - c * u**3)


def overlap_gap_slope(u, c=2.0):
    """``d(overlap_gap)/du`` divided by ``6u``: ``exp(-3u**2) - 1 + (c/2) u``."""
    u = np.asarray(u, float)
    return np.exp(-3 * u**2) - 1 + 0.5 * c * u


# ----------------------------------------------------------- other pages
def yagi_kunii_eq6(x, omega, gamma, x_power=1.0, gamma_coef=1.0):
    """Yagi and Kunii eq. (6) as transcribed on page B3.1: ``theta/theta_B`` against ``x = r/R``.

    ``x_power`` and ``gamma_coef`` change the first numerator term
    ``gamma x`` into ``gamma_coef gamma x**x_power`` (defects).
    """
    x = np.asarray(x, float)
    num = gamma_coef * gamma * x**x_power + 3.0 * omega * gamma * x**2 + omega * (1.0 - 2.0 * gamma) * x**3
    return 1.0 - num / (omega + omega * gamma + gamma)


def eta_sphere(phi):
    """Effectiveness factor of a first-order sphere, ``3 (phi coth phi - 1)/phi**2``."""
    return 3 * (phi / np.tanh(phi) - 1) / phi**2
