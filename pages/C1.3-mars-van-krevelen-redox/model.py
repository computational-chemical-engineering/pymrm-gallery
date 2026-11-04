"""Mars and van Krevelen (1954): the SO2 chain of the redox rate law.

Model code of gallery page C1.3. Every printed constant the functions use sits
in one dictionary, ``PAR``, which every function takes, so a defect in one
constant is one changed entry. The tables themselves are passed in as data
frames; nothing here loads data.

Contents
--------
PAR
    The printed constants: eq. (8), the oxygen content of air, the rate form
    of eq. (15), the converter targets and the temperature window.
k_eq, alpha_eq
    Eq. (8) and the equilibrium conversion it implies for a feed, one root of
    the equilibrium condition on the physical range.
f14
    The right-hand side of eq. (14), as printed or as the exact integral of
    eq. (13a), with or without the retardation term.
fit_neumann, fit_kuster
    The scale (and SO2 fraction, or alpha_eq) that reproduce the authors' k'
    columns, scanned over the whole physical window and refined.
conversion, fit_block
    Conversion against ``k/V`` for the three equal-parameter rate laws, and
    the one rate constant per temperature block from a root of the analytic
    derivative of the sum of squares.
oxygen_from_so2, fit_nu
    Burner-gas O2 from SO2 for a stoichiometry ``nu`` (linear in ``nu``), and
    the closed-form least-squares ``nu``.
OptimumCurve
    The optimum temperature curve ``T_opt(alpha)`` (a root of the analytic
    temperature derivative of the log rate) and the rate along it.
split_quad
    The catalyst split by adaptive quadrature of ``d(alpha)/r``.
OptimumProfileConverter
    The same split from a pymrm plug-flow solve: upwind convection at unit
    velocity plus the rate as a source, ``newton``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.special import beta as beta_fn, betainc, betaincinv, lambertw

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)
from pymrm.checks import find_roots, residual_check

__all__ = ["R_GAS", "CAL", "PAR", "k_eq", "alpha_eq", "f14", "kprime_pred", "fit_neumann",
           "fit_kuster", "conversion", "fit_block", "oxygen_from_so2", "fit_nu",
           "OptimumCurve", "split_quad", "OptimumProfileConverter", "one_root"]

R_GAS = 8.314     # J/(mol K)
CAL = 4184.0      # J per kcal

PAR = {
    "slope": 9880.0,       # eq. (8): 10log K = slope / T - intercept, K in 1/atm
    "intercept": 9.34,
    "x_o2_air": 0.209,     # O2 mole fraction of air (Neumann's gas, Table 4's burner gas)
    "ap_power": 0.5,       # rate (1 - a') / a'^p: p = 1/2 eq. (15), p = 1 eq. (13a)
    "E_fac": 1.0,          # multiplies the activation energy of Table 3's two-point line
    "alpha_mid": 0.81,     # converter 1 outlet, journal page 55
    "alpha_end": 0.97,     # converter 2 outlet
    "T_lo": 450.0,         # K, the window searched for the optimum temperature
    "T_hi": 1400.0,
    "T_cap": None,         # K, an optional ceiling on the optimum curve
    "boreskov_power": 0.8,  # Boreskov-Sokolova (p_SO2 / p_SO3)^0.8, journal page 51
}


def one_root(f, lo, hi, what, **kw):
    """The single root of ``f`` on ``[lo, hi]`` (``pymrm.checks.find_roots``).

    Raises ``RuntimeError`` unless the scan finds exactly one sign change.
    """
    res = find_roots(f, lo, hi, **kw)
    if len(res["roots"]) != 1:
        raise RuntimeError(f"{what}: {res['message']} on [{lo:g}, {hi:g}]")
    return res["roots"][0]


def k_eq(T, par=PAR):
    """Eq. (8): ``K = p_SO3^2 / (p_SO2^2 p_O2)`` in 1/atm at ``T`` (K)."""
    return 10.0 ** (par["slope"] / T - par["intercept"])


def alpha_eq(T, s, w, par=PAR, P=1.0):
    """Equilibrium SO2 conversion of a feed with mole fractions ``s`` (SO2), ``w`` (O2).

    The condition ``(x/(1-x))^2 (1 - s x/2) / ((w - s x/2) P) = K`` rises
    monotonically from 0 to infinity on the physical range ``0 < x <
    min(1, 2 w / s)``; its one root is found there.
    """
    ln_k = math.log(10.0) * (par["slope"] / T - par["intercept"])
    hi = min(1.0, 2.0 * w / s)

    def g(x):
        return (2.0 * math.log(x / (1.0 - x)) + math.log(1.0 - s * x / 2.0)
                - math.log((w - s * x / 2.0) * P) - ln_k)

    return one_root(g, 1e-12 * hi, hi * (1.0 - 1e-12), f"alpha_eq at {T:.2f} K", n_scan=24,
                    xtol=1e-15)


def f14(alpha, aeq, power=1, retard=True):
    """Right-hand side of eq. (14), ``aeq^power {ln 1/(1 - a') - a'}`` with ``a' = alpha/aeq``.

    ``power=1`` as printed, ``power=2`` the exact integral of eq. (13a);
    ``retard=False`` drops the ``-a'`` term (first order to equilibrium).
    Infinite at or beyond equilibrium.
    """
    ap = np.asarray(alpha, dtype=float) / aeq
    with np.errstate(divide="ignore", invalid="ignore"):
        base = -np.log1p(-ap) - (ap if retard else 0.0)
    return np.where(ap < 1.0, aeq ** power * base, np.inf)


def kprime_pred(rows, C, s, par=PAR, power=1, retard=True, aeq_one=False):
    """Recomputed k' of Table 2 rows: ``C V f14(alpha, alpha_eq)``.

    Neumann's gas is taken as SO2 in air, SO2 fraction ``s``, so the O2
    fraction is ``x_o2_air (1 - s)``. ``aeq_one`` sets alpha_eq to 1.
    """
    w = par["x_o2_air"] * (1.0 - s)
    temps = rows.temperature_C.to_numpy()
    aeq = {T: 1.0 if aeq_one else alpha_eq(T + 273.15, s, w, par) for T in np.unique(temps)}
    a = rows.conversion_pct.to_numpy() / 100.0
    f = np.array([f14(ai, aeq[T], power, retard) for ai, T in zip(a, temps)])
    return C * rows.gas_velocity.to_numpy() * f


def _log_fit(obs, f):
    """Closed-form scale of ``C f`` to ``obs`` in log space, and the loss."""
    if not np.all(np.isfinite(f)) or np.any(f <= 0):
        return np.inf, np.nan
    ln_c = np.mean(np.log(obs) - np.log(f))
    return float(np.sum((ln_c + np.log(f) - np.log(obs)) ** 2)), float(np.exp(ln_c))


def _scan_refine(loss, grid):
    """Minimum of ``loss`` scanned on ``grid``, refined between the neighbours of the best point."""
    vals = np.array([loss(x) for x in grid])
    i = int(np.argmin(vals))
    inner = np.isfinite(vals[1:-1]) & (vals[1:-1] < vals[:-2]) & (vals[1:-1] < vals[2:])
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    x = float(minimize_scalar(loss, bounds=(lo, hi), method="bounded",
                              options={"xatol": 1e-12}).x)
    return x, {"grid": np.asarray(grid), "loss": vals, "i_best": i,
               "n_local_minima": int(inner.sum()), "at_edge": i in (0, len(grid) - 1)}


def fit_neumann(rows, par=PAR, power=1, retard=True, s_fixed=None, aeq_one=False, n_scan=120):
    """Fit the scale ``C`` and the SO2 fraction ``s`` to the authors' k' (log least squares).

    ``C`` is the closed-form minimiser at each ``s``. ``s`` is scanned over its
    whole physical window, from trace SO2 to the stoichiometric mixture
    ``s = 2 x_air / (1 + 2 x_air)`` (where ``2 w = s``), and refined around the
    best scan point; ``s_fixed`` or ``aeq_one`` leave only ``C``.

    Returns
    -------
    dict
        ``s``, ``C`` and ``scan`` (the scan record, None when ``s`` is not fitted).
    """
    obs = rows.k_prime.to_numpy()

    def loss(s):
        return _log_fit(obs, kprime_pred(rows, 1.0, s, par, power, retard, aeq_one))

    if aeq_one or s_fixed is not None:
        s = 0.0 if aeq_one else float(s_fixed)
        return {"s": s, "C": loss(s)[1], "scan": None}
    x = par["x_o2_air"]
    s_max = 2.0 * x / (1.0 + 2.0 * x)
    s, scan = _scan_refine(lambda s: loss(s)[0], np.linspace(1e-4, s_max * (1 - 1e-6), n_scan))
    return {"s": s, "C": loss(s)[1], "scan": scan}


def fit_kuster(rows, power=1, retard=True, aeq_one=False, n_scan=400):
    """Fit alpha_eq and the scale to Kuster's k' (Table 1), ``k' = C f14(alpha, aeq) / t``.

    alpha_eq is scanned over its physical window, from the largest measured
    conversion to 1, and refined. Returns ``aeq``, ``C``, ``scan``, ``dev``
    (relative deviation per row).
    """
    obs = rows.k_prime.to_numpy()
    a = rows.conversion_pct.to_numpy() / 100.0
    t = rows.time_s.to_numpy()

    def pred(aeq):
        return f14(a, aeq, power, retard) / t

    if aeq_one:
        aeq, scan = 1.0, None
    else:
        lo = a.max()
        aeq, scan = _scan_refine(lambda q: _log_fit(obs, pred(q))[0],
                                 np.linspace(lo + 1e-9, 1.0, n_scan))
    C = _log_fit(obs, pred(aeq))[1]
    return {"aeq": aeq, "C": C, "scan": scan, "dev": C * pred(aeq) / obs - 1.0}


def conversion(y, aeq, law, par=PAR):
    """Conversion and its derivative in ``y = k/V`` for one of three rate laws.

    ``law``: ``"mvk"`` eq. (14) as printed, ``aeq {ln 1/(1-a') - a'} = y``
    (closed form through the Lambert W function); ``"first_order"``,
    ``aeq ln 1/(1-a') = y``; ``"boreskov"``, ``d alpha/dt = k ((1-alpha)/alpha)^0.8``
    (an incomplete beta function), capped at alpha_eq since the law has no
    equilibrium term.

    Returns
    -------
    tuple of ndarray
        ``alpha`` and ``d alpha / d y``.
    """
    y = np.asarray(y, dtype=float)
    if law == "mvk":
        # 1 - a' = q with q - ln q = 1 + y/aeq, so -q = W0(-exp(-1 - y/aeq))
        q = -np.real(lambertw(-np.exp(-1.0 - y / aeq), 0))
        return aeq * (1.0 - q), q / (1.0 - q)
    if law == "first_order":
        q = np.exp(-y / aeq)
        return aeq * (1.0 - q), q
    if law == "boreskov":
        p = par["boreskov_power"]
        a, b = 1.0 + p, 1.0 - p
        total = beta_fn(a, b)
        y_cap = total * betainc(a, b, aeq)
        al = np.where(y < y_cap, betaincinv(a, b, np.minimum(y, y_cap) / total), aeq)
        with np.errstate(divide="ignore"):
            d = np.where(y < y_cap, ((1.0 - al) / al) ** p, 0.0)
        return al, d
    raise ValueError(f"unknown rate law {law!r}")


def fit_block(rows, aeq, law, par=PAR, k_max=None):
    """One rate constant ``k`` for a temperature block, least squares in conversion points.

    ``k`` is the root of the analytic derivative of the sum of squares, on a
    logarithmic scan from a hundredth of the smallest to thirty times the
    largest velocity (for Boreskov's law, up to where every row is capped), or
    up to ``k_max``.

    Returns
    -------
    dict
        ``k``, ``rms`` (conversion points), ``at_edge`` (True when the sum of
        squares still falls at the end of the window) and ``n_stationary`` (sign
        changes of the derivative; the lowest of them is taken).
    """
    V = rows.gas_velocity.to_numpy()
    cv = rows.conversion_pct.to_numpy()

    def sse(k):
        return float(np.sum((100.0 * conversion(k / V, aeq, law, par)[0] - cv) ** 2))

    def dsse(k):
        al, d = conversion(k / V, aeq, law, par)
        return float(np.sum(2.0 * (100.0 * al - cv) * 100.0 * d / V))

    lo, hi = 0.01 * V.min(), 30.0 * V.max()
    if law == "boreskov":
        p = par["boreskov_power"]
        hi = min(hi, V.max() * beta_fn(1 + p, 1 - p) * betainc(1 + p, 1 - p, aeq) * (1 - 1e-9))
    if k_max is not None:
        hi = min(hi, k_max)
    if dsse(hi) < 0.0:
        return {"k": hi, "rms": math.sqrt(sse(hi) / len(V)), "at_edge": True, "n_stationary": 0}
    roots = find_roots(dsse, lo, hi, log=True)["roots"]
    if not roots:
        raise RuntimeError(f"fit_block {law}: no stationary point on [{lo:g}, {hi:g}]")
    # Boreskov's cap makes the sum of squares kinked: keep the lowest stationary point
    k = min(roots, key=sse)
    return {"k": k, "rms": math.sqrt(sse(k) / len(V)), "at_edge": False,
            "n_stationary": len(roots)}


def oxygen_from_so2(pct_so2, nu, par=PAR):
    """%O2 of burner gas from air, ``nu`` mol O2 consumed per mol SO2 formed.

    Per 100 mol air burnt, SO2 formed ``n`` and O2 left ``100 x_air - nu n`` in
    ``100 + (1 - nu) n`` mol; eliminating ``n`` with %SO2 gives
    ``%O2 = 100 x_air - %SO2 (x_air + (1 - x_air) nu)``, linear in ``nu``.
    """
    x = par["x_o2_air"]
    return 100.0 * x - np.asarray(pct_so2) * (x + (1.0 - x) * nu)


def fit_nu(pct_o2, pct_so2, par=PAR):
    """Least-squares ``nu`` of :func:`oxygen_from_so2` (closed form, the model is linear in ``nu``)."""
    x = par["x_o2_air"]
    p, o = np.asarray(pct_so2, float), np.asarray(pct_o2, float)
    return float(np.sum(p * (100.0 * x - x * p - o)) / ((1.0 - x) * np.sum(p * p)))


class OptimumCurve:
    """The optimum temperature curve and the rate along it, for one feed.

    The rate is eq. (15)'s diffusion-limited form ``A(T) (1 - a') / a'^p`` with
    ``a' = alpha / alpha_eq(T)``; ``A(T)`` is the Arrhenius line through
    Table 3's two rate groups. At each conversion the optimum temperature is
    the root of ``d ln r / dT = E/(R T^2) + (1/(1 - a') + p/a') a' alpha_eq'/alpha_eq``
    on the window from ``T_lo`` to where ``alpha_eq(T) = alpha`` (or ``T_hi``,
    or ``T_cap``). ``alpha_eq(T)`` is a cubic spline through ``n_T`` exact
    values, and the rate on the curve is ``g(alpha) / alpha^p`` with ``g`` a
    cubic spline through ``n_alpha`` conversions; both routes to the split use
    these splines.

    Parameters
    ----------
    group : tuple
        ``((T1, A1), (T2, A2))``: Table 3's temperatures (K) and rate groups (atm/hr).
    s, w : float
        Feed mole fractions of SO2 and O2.
    par : dict, optional
        The printed constants (``PAR``).
    n_T, n_alpha : int, optional
        Spline nodes in temperature and in conversion.
    """

    def __init__(self, group, s, w, par=PAR, n_T=1200, n_alpha=350):
        (T1, A1), (T2, A2) = group
        self.par, self.s, self.w = par, s, w
        self.e_over_r = math.log(A2 / A1) / (1.0 / T1 - 1.0 / T2) * par["E_fac"]
        self.T1, self.A1 = T1, A1
        self.p = par["ap_power"]
        T_grid = np.linspace(par["T_lo"], par["T_hi"], n_T)
        self.aeq_spline = CubicSpline(T_grid, [alpha_eq(T, s, w, par) for T in T_grid])
        self.daeq_spline = self.aeq_spline.derivative()
        self.a_nodes = np.linspace(1e-4, 0.9745, n_alpha)
        self.T_nodes = np.array([self.T_opt(a) for a in self.a_nodes])
        self.n_at_ceiling = int(np.sum(self.T_nodes >= self._ceiling() - 1e-9))
        # the rate grows as alpha^-p towards the inlet: spline the smooth r alpha^p
        self.g_spline = CubicSpline(self.a_nodes, [self.rate(a, T) * a ** self.p for a, T in
                                                   zip(self.a_nodes, self.T_nodes)])

    def A(self, T):
        """Table 3's rate group on the two-point Arrhenius line (atm/hr)."""
        return self.A1 * math.exp(-self.e_over_r * (1.0 / T - 1.0 / self.T1))

    def rate(self, alpha, T):
        """``A(T) (1 - a') / a'^p``; zero at or beyond equilibrium."""
        ap = alpha / float(self.aeq_spline(T))
        return 0.0 if ap >= 1.0 else self.A(T) * (1.0 - ap) / ap ** self.p

    def dlnr_dT(self, alpha, T):
        """Analytic temperature derivative of ``ln rate`` at fixed conversion."""
        aeq = float(self.aeq_spline(T))
        ap = alpha / aeq
        return (self.e_over_r / T ** 2
                + (1.0 / (1.0 - ap) + self.p / ap) * ap * float(self.daeq_spline(T)) / aeq)

    def _ceiling(self):
        cap = self.par["T_cap"]
        return self.par["T_hi"] if cap is None else min(self.par["T_hi"], cap)

    def T_eq(self, alpha):
        """Temperature at which ``alpha`` is the equilibrium conversion, or None above ``T_hi``."""
        if float(self.aeq_spline(self.par["T_hi"])) > alpha:
            return None
        return one_root(lambda T: float(self.aeq_spline(T)) - alpha, self.par["T_lo"],
                        self.par["T_hi"], f"T_eq at alpha {alpha:.4f}")

    def T_opt(self, alpha):
        """Optimum temperature (K) at conversion ``alpha``; a window end when the rate peaks there."""
        T_eq = self.T_eq(alpha)
        hi = self._ceiling() if T_eq is None else min(T_eq - 1e-6, self._ceiling())
        lo = self.par["T_lo"]
        if self.dlnr_dT(alpha, hi) > 0.0:
            return hi
        if self.dlnr_dT(alpha, lo) < 0.0:
            return lo
        return one_root(lambda T: self.dlnr_dT(alpha, T), lo, hi, f"T_opt at alpha {alpha:.4f}")

    def R(self, alpha):
        """Rate on the optimum curve (atm/hr), with ``alpha`` clipped to the node range."""
        a = np.clip(alpha, self.a_nodes[0], self.a_nodes[-1])
        return self.g_spline(a) / a ** self.p


def split_quad(curve, par=PAR, epsrel=1e-10, scale=None):
    """Catalyst split along the optimum curve by adaptive quadrature of ``d(alpha)/r``.

    The integral is taken in ``u = sqrt(alpha)``, where the integrand is smooth.

    ``scale(alpha)`` multiplies the rate (a uniform factor must leave the split
    unchanged). Returns the converter-1 share (%) and the total contact time.
    """
    sc = (lambda a: 1.0) if scale is None else scale

    def inv(u):
        # alpha = u^2 removes the square-root behaviour of 1/r at the inlet
        a = u * u
        return 2.0 * u / (sc(a) * float(curve.R(a)))

    kw = dict(limit=400, epsabs=0.0, epsrel=epsrel)
    u_mid, u_end = math.sqrt(par["alpha_mid"]), math.sqrt(par["alpha_end"])
    w1 = quad(inv, 0.0, u_mid, **kw)[0]
    w2 = quad(inv, u_mid, u_end, **kw)[0]
    return 100.0 * w1 / (w1 + w2), w1 + w2


class OptimumProfileConverter:
    """Steady plug flow ``d(alpha)/d(tau) = R(alpha)`` along the optimum temperature curve.

    ``tau`` is contact time (proportional to catalyst amount); the velocity is
    1 by construction, so upwind convection of ``alpha`` at ``v = 1`` plus a
    pointwise source is the whole model. State layout ``(n_tau, 1)``.

    Parameters
    ----------
    R : callable
        Rate on the optimum curve, ``R(alpha)`` elementwise.
    tau_end : float
        Bed length in contact time.
    n_tau : int
        Number of cells.
    """

    def __init__(self, R, tau_end, n_tau):
        self.R = R
        self.shape = (int(n_tau), 1)
        self.tau_f = np.linspace(0.0, tau_end, self.shape[0] + 1)
        self.tau_c = 0.5 * (self.tau_f[:-1] + self.tau_f[1:])
        # tau = 0, outward normal -tau: the feed is unconverted, alpha = 0 -> a = 0, b = 1, d = 0
        # tau = tau_end, outward normal +tau: a pure convective outlet carries the
        #     last cell's value out; no condition is imposed there
        self.bc = ({"a": 0.0, "b": 1.0, "d": 0.0}, {"outflow": True})
        self.alpha = np.zeros(self.shape)
        self.result = None
        self._build_operators()

    def _build_operators(self):
        # nu = 0: contact time is a Cartesian coordinate
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.tau_f, self.tau_c,
                                                      self.bc, v=1.0, axis=0)
        div_mat = construct_div(self.shape, self.tau_f, nu=0, axis=0)
        self.jac_const = (div_mat @ conv_mat).tocsc()
        self.g_const = (div_mat @ conv_bc).toarray().reshape(-1, 1)
        # the source at a cell depends on that cell's conversion only
        self.numjac = NumJac(self.shape)

    def reaction(self, u):
        """Pointwise source ``R(alpha)``, shape ``(n_tau, 1)``."""
        return self.R(u)

    def residual(self, u):
        """Steady residual ``d(v alpha)/d tau - R(alpha)`` and its Jacobian."""
        g_rxn, jac_rxn = self.numjac(self.reaction, u.reshape(self.shape))
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self, maxfev=100, strict=True):
        """Newton from a linear profile; raises ``RuntimeError`` unless it converges."""
        u0 = np.linspace(0.05, 0.97, self.shape[0]).reshape(self.shape)
        result = newton(self.residual, u0, maxfev=maxfev)
        if strict and not result.success:
            raise RuntimeError(f"OptimumProfileConverter.solve: {result.message}")
        self.alpha, self.result = result.x.reshape(self.shape), result
        return self

    def backward_error(self):
        """Componentwise backward error of the current state (``pymrm.checks``)."""
        return residual_check(self.residual, self.alpha.reshape(-1))["backward_error"]

    def outlet(self):
        """Conversion on the outlet face, read with ``compute_boundary_values``."""
        _, _, a_out, _ = compute_boundary_values(self.alpha, self.tau_f, self.tau_c, self.bc,
                                                 axis=0)
        return float(np.asarray(a_out).ravel()[0])

    def profile(self, tau, read="face"):
        """Conversion at ``tau``, linear between face values or between cell centres.

        First-order upwind at ``v = 1`` is backward Euler in ``tau``: each cell
        value approximates the profile at the cell's downstream face, which
        ``"face"`` uses (with the inlet value at ``tau = 0``); ``"centre"``
        places it at the cell centre; ``"step"`` takes the value of the cell
        that contains ``tau``, with no interpolation.
        """
        a = self.alpha[:, 0]
        if read == "face":
            return np.interp(tau, self.tau_f, np.concatenate([[0.0], a]))
        if read == "centre":
            return np.interp(tau, self.tau_c, a)
        if read == "step":
            i = np.searchsorted(self.tau_f, tau, side="right") - 1
            return a[np.clip(i, 0, a.size - 1)]
        raise ValueError(f"unknown read {read!r}")

    def crossing(self, target, read="face"):
        """Contact time at which the profile reaches ``target``, the one root on the bed."""
        lo = 0.0 if read == "face" else self.tau_c[0]
        hi = self.tau_f[-1] if read == "face" else self.tau_c[-1]
        if self.profile(hi, read) < target:
            raise RuntimeError(f"crossing: the bed ends at alpha {self.profile(hi, read):.5f}, "
                               f"below {target}; lengthen it")
        return one_root(lambda t: self.profile(t, read) - target, lo, hi,
                        f"crossing alpha = {target}", n_scan=200)
