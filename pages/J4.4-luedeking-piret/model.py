"""Luedeking and Piret (1959): product formation tied to growth and to biomass (page J4.4).

Luedeking and Piret, J. Biochem. Microbiol. Technol. Eng. 1, 393 (1959). The
relation is eq. (2), ``dP/dt = alpha dN/dt + beta N``, or eq. (3) after division
by ``N``, ``(1/N) dP/dt = alpha k + beta`` with ``k = (1/N) dN/dt``.

- the three least-squares routes to ``(alpha, beta)`` from the printed tables
  and the two one-term nulls the paper rejects (:func:`refit_all`,
  :func:`null_split`, :func:`ols`), and the transcription test of Table II's
  ``k`` column (:func:`digit_detection`);
- the closed form of eq. (2) with ``k(t)`` prescribed from the printed ``k``
  column (:class:`ClosedForm`);
- a batch fermenter marched by backward Euler with ``newton`` on
  ``NumJac((1, 2))`` (:func:`batch_march`) and a steady plug-flow fermenter
  built from ``construct_convflux_upwind`` and ``construct_div``
  (:class:`PlugFlowFermenter`);
- crossings found by ``pymrm.checks.find_roots`` on the whole physical range,
  raising unless there is exactly one (:func:`unique_root`).

Every defect the page's break table injects enters through an argument of these
functions. No data are loaded here; the tables are passed in as DataFrames.
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator, make_interp_spline
from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)
from pymrm.checks import find_roots

K_LOG = 0.45          # 1/h, the gap in the printed k column (0.40 to 0.46): not a fit
QUAD_RTOL = 1e-13     # relative tolerance of the closed-form integral of N


# ------------------------------------------------------------------ helpers
def ols(X, y):
    """Ordinary least squares coefficients of ``y ~ X``.

    Parameters
    ----------
    X : ndarray, shape (m, p)
        Regressors, one column per coefficient.
    y : ndarray, shape (m,)
        Observations.

    Returns
    -------
    ndarray, shape (p,)
        The least-squares coefficients.
    """
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return c


def rms(r):
    """Root mean square of a residual vector.

    Parameters
    ----------
    r : array_like
        Residuals.

    Returns
    -------
    float
        ``sqrt(mean(r**2))``.
    """
    return float(np.sqrt((np.asarray(r, dtype=float) ** 2).mean()))


def interpolant(t, v, kind="pchip"):
    """Interpolant of a printed column, with an ``antiderivative`` method.

    Parameters
    ----------
    t, v : array_like
        Printed abscissae (increasing) and values.
    kind : {"pchip", "linear"}
        Monotone cubic (the page's choice) or piecewise linear (a break row).

    Returns
    -------
    PchipInterpolator or BSpline
        Callable on ``[t[0], t[-1]]``.
    """
    t, v = np.asarray(t, dtype=float), np.asarray(v, dtype=float)
    if kind == "pchip":
        return PchipInterpolator(t, v)
    if kind == "linear":
        return make_interp_spline(t, v, k=1)
    raise ValueError(f"unknown interpolant kind {kind!r}")


def unique_root(f, lo, hi, what, n_scan=256):
    """The single root of ``f`` on ``[lo, hi]``.

    Parameters
    ----------
    f : callable
        Scalar function of one variable.
    lo, hi : float
        The physical range scanned in full.
    what : str
        Name used in the error message.
    n_scan : int
        Number of scan intervals handed to ``find_roots``.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        Unless ``find_roots`` reports exactly one root, so a second crossing or
        none is never replaced silently by "the first root".
    """
    found = find_roots(lambda x: float(f(x)), lo, hi, n_scan=n_scan, xtol=1e-13)
    if len(found["roots"]) != 1 or found["failed"]:
        raise ValueError(f"{what}: expected one root on [{lo:g}, {hi:g}], "
                         f"got {found['message']}")
    return float(found["roots"][0])


# -------------------------------------------------------------- closed form
class ClosedForm:
    """Eq. (2) integrated with ``k(t)`` prescribed from the printed ``k`` column.

    ``N(t) = N0 exp(int k dt)`` and ``P(t) = P0 + alpha (N - N0) + beta int N dt``,
    with the integral of ``N`` by adaptive quadrature broken at the printed
    times (relative tolerance :data:`QUAD_RTOL`).

    Parameters
    ----------
    t_k, k : array_like
        Printed times (h) and specific growth rates (1/h).
    n0 : float
        ``N`` at the first printed ``k`` time, U.O.D./ml.
    kind : {"pchip", "linear"}
        Interpolant of the ``k`` column.
    """

    def __init__(self, t_k, k, n0, kind="pchip"):
        self.t_k = np.asarray(t_k, dtype=float)
        self.kt = interpolant(self.t_k, k, kind)
        self._ka = self.kt.antiderivative()
        self.t_lo, self.t_hi, self.n0 = float(self.t_k[0]), float(self.t_k[-1]), float(n0)

    def k(self, t):
        """Prescribed specific growth rate at ``t``, clipped to the printed range."""
        return self.kt(np.clip(t, self.t_lo, self.t_hi))

    def N(self, t):
        """Bacterial density from the integrated printed ``k`` column."""
        return self.n0 * np.exp(self._ka(t) - self._ka(self.t_lo))

    def int_N(self, t):
        """``int_{t_lo}^{t} N dt``, with its quadrature error estimate.

        Returns
        -------
        tuple of float
            The integral and ``quad``'s absolute error estimate.
        """
        pts = self.t_k[(self.t_k > self.t_lo) & (self.t_k < t)]
        val, err = quad(self.N, self.t_lo, t, points=pts, epsabs=0.0,
                        epsrel=QUAD_RTOL, limit=400)
        return float(val), float(err)

    def P(self, t, a, b, p0=0.0):
        """Acid concentration from eq. (2) integrated, mg/ml."""
        return float(p0 + a * (self.N(t) - self.n0) + b * self.int_N(t)[0])


# ------------------------------------------------------------ pymrm models
def batch_march(nt, a, b, cf, tol=1e-13):
    """Backward Euler on ``(N, P)`` with ``newton`` and ``NumJac((1, 2))``.

    ``dN/dt = k N`` and ``dP/dt = alpha k N + beta N``, with ``k`` taken at the
    end of each step, ``t_lo + (i + 1) dt``. The step is linear in ``(N, P)``,
    so newton's first iterate already solves it and the next confirms it.
    ``NumJac((1, 2))`` has the field index as its last axis, so its default
    stencil couples the two fields of the one cell in full.

    Parameters
    ----------
    nt : int
        Number of steps over ``[cf.t_lo, cf.t_hi]``.
    a, b : float
        alpha (mg/U.O.D.) and beta (mg/U.O.D./h).
    cf : ClosedForm
        Supplies ``k(t)`` and the start state ``(N0, 0)``.
    tol : float
        newton step tolerance (with ``stagnation=True``).

    Returns
    -------
    ndarray, shape (2,)
        ``(N, P)`` at ``cf.t_hi``.

    Raises
    ------
    RuntimeError
        When a newton solve does not report success.
    """
    dt = (cf.t_hi - cf.t_lo) / nt
    jac = NumJac((1, 2))
    y = np.array([[cf.n0, 0.0]])
    for i in range(nt):
        kc = float(cf.k(cf.t_lo + (i + 1) * dt))
        y_old = y.copy()

        def src(c, kc=kc):
            N = c[..., 0:1]
            return np.concatenate([kc * N, a * kc * N + b * N], axis=-1)

        def resid(c, y_old=y_old):
            s, js = jac(src, c.reshape((1, 2)))
            r = (c.reshape((1, 2)) - y_old) / dt - s
            return r.reshape((-1, 1)), np.eye(2) / dt - js.toarray()

        sol = newton(resid, y.reshape((1, 2)), tol=tol, stagnation=True, maxfev=50)
        if not sol.success:
            raise RuntimeError(f"batch step {i} of {nt}: newton failed ({sol.message})")
        y = sol.x.reshape((1, 2))
    return y[0]


class PlugFlowFermenter:
    """Steady plug flow carrying ``(N, P)`` with the Luedeking-Piret source.

    ``v dc/dz = s(c)`` on ``0 < z < v tau``, where residence time is
    ``t = t_lo + z / v``, so ``k(z)`` is a constant vector assembled once. The
    operators are ``construct_convflux_upwind`` and ``construct_div`` (``nu=0``,
    Cartesian), assembled once in ``__init__``; the source is pointwise in the
    two fields, ``NumJac((n, 2))``.

    Parameters
    ----------
    ncell : int
        Number of cells.
    a, b : float
        alpha and beta.
    cf : ClosedForm
        Supplies ``k(t)`` and the feed ``(N0, 0)``.
    v : float
        Velocity; the result depends on ``z / v`` only.
    outlet : {"outflow", "zero_gradient"}
        Outlet condition: ``{"outflow": True}`` (the face carries the last cell's
        value, the first-order convection-only outlet) or ``{"a": 1, "b": 0,
        "d": 0}`` (a break row).
    tol : float
        newton step tolerance (with ``stagnation=True``).
    """

    def __init__(self, ncell, a, b, cf, v=1.0, outlet="outflow", tol=1e-12):
        self.shape, self.a, self.b, self.tol = (ncell, 2), a, b, tol
        self.z_f = np.linspace(0.0, v * (cf.t_hi - cf.t_lo), ncell + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        d0 = np.array([[cf.n0, 0.0]])
        # outward normal: a dc/dn + b c = d
        #   z = 0: c = (N0, 0), a Dirichlet feed
        #   z = L: pure outflow, the exit face carries the last cell's value
        out = {"outflow": True} if outlet == "outflow" else \
            {"a": 1.0, "b": 0.0, "d": 0.0}           # dc/dz = 0 at z = L
        if outlet not in ("outflow", "zero_gradient"):
            raise ValueError(f"unknown outlet {outlet!r}")
        self.bc = ({"a": 0.0, "b": 1.0, "d": d0}, out)
        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                  self.bc, v=v, axis=0)
        div = construct_div(self.shape, self.z_f, nu=0, axis=0)  # nu=0: Cartesian
        self.A, self.rhs = div @ conv, div @ conv_bc
        self.numjac = NumJac(self.shape)          # last axis = the 2 fields
        self.k_z = cf.k(cf.t_lo + self.z_c / v)[:, None]
        self.c0 = np.tile(d0, (ncell, 1))

    def source(self, c):
        """``(k N, alpha k N + beta N)`` in every cell."""
        N = c[..., 0:1]
        return np.concatenate([self.k_z * N, self.a * self.k_z * N + self.b * N], axis=-1)

    def solve(self):
        """Solve the steady state; raises unless newton reports success."""
        def resid(c):
            s, js = self.numjac(self.source, c)
            return (self.rhs + self.A @ c.reshape((-1, 1)) - s.reshape((-1, 1)),
                    self.A - js)
        r = newton(resid, self.c0, tol=self.tol, stagnation=True, maxfev=50)
        if not r.success:
            raise RuntimeError(f"plug-flow solve, {self.shape[0]} cells: newton failed "
                               f"({r.message})")
        self.c = r.x.reshape(self.shape)
        return self

    def outlet(self):
        """Outlet ``(N, P)`` from ``compute_boundary_values`` under the model's bc."""
        return np.asarray(compute_boundary_values(
            self.c, self.z_f, self.z_c, self.bc, axis=0)[2]).reshape(-1)[:2]

    def outlet_last_cell(self):
        """``(N, P)`` in the last cell, for the comparison with :meth:`outlet`."""
        return self.c[self.shape[0] - 1].copy()


# --------------------------------------------------------------- the fits
def digit_detection(rows, tol):
    """Single-digit substitutions in the ``(N, dN/dt, k)`` triple the k identity catches.

    Parameters
    ----------
    rows : DataFrame
        Table II rows with ``N_UOD_per_ml``, ``dNdt_UOD_per_ml_h`` and ``k_per_h``.
    tol : float
        A substitution is caught when ``|dN/dt / N / k - 1|`` exceeds ``tol``.

    Returns
    -------
    tuple of int
        ``(caught, total)``.
    """
    cols = ["N_UOD_per_ml", "dNdt_UOD_per_ml_h", "k_per_h"]
    total = caught = 0
    for _, r in rows.iterrows():
        for c in cols:
            s = f"{r[c]:g}"
            for j, ch in enumerate(s):
                if not ch.isdigit():
                    continue
                for d in "0123456789":
                    if d == ch:
                        continue
                    v = float(s[:j] + d + s[j + 1:])
                    if v <= 0:
                        continue
                    trip = {k: float(r[k]) for k in cols}
                    trip[c] = v
                    total += 1
                    if abs(trip["dNdt_UOD_per_ml_h"] / trip["N_UOD_per_ml"]
                           / trip["k_per_h"] - 1.0) > tol:
                        caught += 1
    return caught, total


def null_split(df, k_min=None, k_max=None):
    """The two-term eq. (2) fit and its two one-term nulls on a window of ``k``.

    Parameters
    ----------
    df : DataFrame
        Table II rows with ``N``, ``dN/dt``, ``dP/dt`` and ``k``.
    k_min, k_max : float, optional
        Keep rows with ``k >= k_min`` and ``k < k_max``.

    Returns
    -------
    tuple
        ``(rows, two-term rms, growth-only rms ratio, density-only rms ratio,
        alpha, beta)``.
    """
    f = df
    if k_min is not None:
        f = f[f.k_per_h >= k_min]
    if k_max is not None:
        f = f[f.k_per_h < k_max]
    x = np.column_stack([f.dNdt_UOD_per_ml_h.to_numpy(), f.N_UOD_per_ml.to_numpy()])
    y = f.dPdt_mg_per_ml_h.to_numpy()
    c = ols(x, y)
    r2 = rms(y - x @ c)
    a1 = float(ols(x[:, [0]], y)[0])
    b1 = float(ols(x[:, [1]], y)[0])
    return (len(f), r2, rms(y - a1 * x[:, 0]) / r2, rms(y - b1 * x[:, 1]) / r2,
            float(c[0]), float(c[1]))


def refit_all(t1, t2, t3, kinterp="pchip", quad_rule="trapz", use_quotients=False,
              drop_first_P_row=False, k_min=None, k_log=K_LOG, ph=6.0, ph_text=None):
    """Every table-derived quantity of the page, from (possibly corrupted) tables.

    Parameters
    ----------
    t1, t2, t3 : DataFrame
        Tables I, II and III as transcribed.
    kinterp : {"pchip", "linear"}
        Interpolant of the printed ``N`` and ``k`` columns.
    quad_rule : {"trapz", "rect"}
        Quadrature of the integrated acid split over the printed rows.
    use_quotients : bool
        Fit eq. (3) on ``dP/dt / N`` instead of the printed ``(1/N) dP/dt``.
    drop_first_P_row : bool
        Start the integral form at the second printed ``P``.
    k_min : float, optional
        Fit eqs. (2) and (3) on rows with ``k >= k_min`` only.
    k_log : float
        Split between the logarithmic plateau and the rest.
    ph : float
        The Table III row compared with the refit.
    ph_text : sequence of float, optional
        The pH levels as printed a second time in the running text; when given,
        ``ph_column_matches_text`` says whether Table III's column equals them.

    Returns
    -------
    dict
        Name to value, the names of the page's claims.
    """
    out = {}
    f3 = t2.dropna(subset=["k_per_h", "invN_dPdt_mg_per_UOD_h", "N_UOD_per_ml",
                           "dPdt_mg_per_ml_h"])
    if k_min is not None:
        f3 = f3[f3.k_per_h >= k_min]
    k = f3.k_per_h.to_numpy()
    y = (f3.dPdt_mg_per_ml_h / f3.N_UOD_per_ml).to_numpy() if use_quotients \
        else f3.invN_dPdt_mg_per_UOD_h.to_numpy()
    X = np.column_stack([k, np.ones_like(k)])
    c = ols(X, y)
    r = y - X @ c
    cov = (r ** 2).sum() / (len(k) - 2) * np.linalg.inv(X.T @ X)
    out["n_rows_eq3"] = len(k)
    out["alpha_refit_eq3"], out["beta_refit_eq3"] = float(c[0]), float(c[1])
    out["eq3_refit_rms"] = rms(r)
    out["eq3_R2"] = float(1 - (r ** 2).sum() / ((y - y.mean()) ** 2).sum())
    out["alpha_se_eq3"] = float(np.sqrt(cov[0, 0]))
    out["beta_se_eq3"] = float(np.sqrt(cov[1, 1]))
    out["alpha_beta_corr_eq3"] = float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]))
    r3 = t3.loc[t3.pH == ph].iloc[0]
    ap, bp = float(r3.alpha_mg_lactic_per_UOD), float(r3.beta_mg_lactic_per_UOD_h)
    out["alpha_printed"], out["beta_printed"] = ap, bp
    out["alpha_refit_eq3_vs_printed"] = float(c[0] / ap - 1)
    out["beta_refit_eq3_vs_printed"] = float(c[1] / bp - 1)
    out["printed_alpha_in_se"] = float(abs(c[0] - ap) / np.sqrt(cov[0, 0]))
    out["printed_beta_in_se"] = float(abs(c[1] - bp) / np.sqrt(cov[1, 1]))
    out["beta_given_printed_alpha"] = float(np.mean(y - ap * k))
    out["alpha_given_printed_beta"] = float(np.sum(k * (y - bp)) / np.sum(k * k))
    out["eq3_printed_pair_rms"] = rms(y - ap * k - bp)
    out["eq3_printed_pair_rms_penalty"] = float(rms(y - ap * k - bp) / rms(r) - 1)
    out["k_star_refit"] = float(c[1] / c[0])

    f2 = t2.dropna(subset=["N_UOD_per_ml", "dNdt_UOD_per_ml_h", "dPdt_mg_per_ml_h"])
    if k_min is not None:
        f2 = f2[f2.k_per_h >= k_min]
    x2 = np.column_stack([f2.dNdt_UOD_per_ml_h.to_numpy(), f2.N_UOD_per_ml.to_numpy()])
    dp = f2.dPdt_mg_per_ml_h.to_numpy()
    c2 = ols(x2, dp)
    rms2 = rms(dp - x2 @ c2)
    out["alpha_refit_eq2"], out["beta_refit_eq2"] = float(c2[0]), float(c2[1])
    out["eq2_refit_rms"] = rms2
    a1 = float(ols(x2[:, [0]], dp)[0])
    b1 = float(ols(x2[:, [1]], dp)[0])
    out["eq2_null_growth_only_rms"] = rms(dp - a1 * x2[:, 0])
    out["eq2_null_density_only_rms"] = rms(dp - b1 * x2[:, 1])
    out["eq2_null_growth_only_ratio"] = out["eq2_null_growth_only_rms"] / rms2
    out["eq2_null_density_only_ratio"] = out["eq2_null_density_only_rms"] / rms2

    # both halves of the paper's qualifier, on the unfiltered table so that the
    # log-phase-only row cannot empty one of the two windows
    fsp = t2.dropna(subset=["N_UOD_per_ml", "dNdt_UOD_per_ml_h", "dPdt_mg_per_ml_h",
                            "k_per_h"])
    n_lg, r_lg, ga_lg, gb_lg, a_lg, b_lg = null_split(fsp, k_min=k_log)
    n_ot, r_ot, ga_ot, gb_ot, _, _ = null_split(fsp, k_max=k_log)
    out["n_rows_logphase"], out["n_rows_outside_log"] = n_lg, n_ot
    out["eq2_logphase_refit_rms"] = r_lg
    out["eq2_logphase_null_growth_only_ratio"] = ga_lg
    out["eq2_logphase_null_density_only_ratio"] = gb_lg
    out["alpha_logphase_eq2"], out["beta_logphase_eq2"] = a_lg, b_lg
    out["eq2_outside_log_refit_rms"] = r_ot
    out["eq2_outside_log_null_growth_only_ratio"] = ga_ot
    out["eq2_outside_log_null_density_only_ratio"] = gb_ot

    # integral form: P itself, no differentiated column
    kk = t2.dropna(subset=["k_per_h"])
    t_lo, t_hi = float(kk.time_h.min()), float(kk.time_h.max())
    nn = t2.dropna(subset=["N_UOD_per_ml"])
    nt = interpolant(nn.time_h, nn.N_UOD_per_ml, kinterp)
    nti = nt.antiderivative()
    fp = t2.dropna(subset=["P_mg_per_ml"])
    fp = fp[fp.time_h <= t_hi]
    if drop_first_P_row:
        fp = fp.iloc[1:]
    tp, po = fp.time_h.to_numpy(), fp.P_mg_per_ml.to_numpy()
    dn, inn = nt(tp) - nt(tp[0]), nti(tp) - nti(tp[0])
    xc = np.column_stack([dn, inn])
    cc = ols(xc, po - po[0])
    rc = rms(po - po[0] - xc @ cc)
    out["alpha_refit_integral"], out["beta_refit_integral"] = float(cc[0]), float(cc[1])
    out["integral_refit_rms"] = rc
    ai = float(ols(xc[:, [0]], po - po[0])[0])
    bi = float(ols(xc[:, [1]], po - po[0])[0])
    out["integral_null_growth_only_ratio"] = rms(po - po[0] - ai * xc[:, 0]) / rc
    out["integral_null_density_only_ratio"] = rms(po - po[0] - bi * xc[:, 1]) / rc
    out["alpha_route_spread"] = float(max(c[0], c2[0], cc[0]) / min(c[0], c2[0], cc[0]) - 1)
    out["beta_route_spread"] = float(max(c[1], c2[1], cc[1]) / min(c[1], c2[1], cc[1]) - 1)
    pred = po[0] + ap * dn + bp * inn
    m = tp > tp[0]
    out["P_curve_rms_rel"] = rms(pred[m] / po[m] - 1)
    out["P_curve_max_rel"] = float(np.abs(pred[m] / po[m] - 1).max())
    out["P_final_predicted"] = float(pred[-1])
    out["P_final_vs_printed"] = float(pred[-1] / po[-1] - 1)

    # the crossover k = beta/alpha, on two column paths that share no cell
    kt = interpolant(kk.time_h, kk.k_per_h, kinterp)
    out["k_star_printed"] = float(bp / ap)
    out["crossover_time_h"] = unique_root(lambda t: kt(t) - bp / ap, t_lo, t_hi,
                                          "k(t) = beta/alpha")
    dd = t2.dropna(subset=["dNdt_UOD_per_ml_h", "N_UOD_per_ml"])
    dnt = interpolant(dd.time_h, dd.dNdt_UOD_per_ml_h, kinterp)
    nnt = interpolant(dd.time_h, dd.N_UOD_per_ml, kinterp)
    out["crossover_time_dNdt_route_h"] = unique_root(
        lambda t: ap * dnt(t) - bp * nnt(t), t_lo, t_hi, "alpha dN/dt = beta N")
    out["crossover_two_routes_rel"] = float(abs(out["crossover_time_h"]
                                                / out["crossover_time_dNdt_route_h"] - 1))
    kt_lin = interpolant(kk.time_h, kk.k_per_h, "linear")
    t_lin = unique_root(lambda t: kt_lin(t) - bp / ap, t_lo, t_hi, "linear k(t) = beta/alpha")
    out["crossover_interp_band_rel"] = float(abs(t_lin / out["crossover_time_h"] - 1))
    out["crossover_window_fraction"] = float((out["crossover_time_h"] - t_lo) / (t_hi - t_lo))

    # the acid split over the printed window
    tf2 = f2.time_h.to_numpy()
    if quad_rule == "trapz":
        ga, ng = ap * np.trapezoid(x2[:, 0], tf2), bp * np.trapezoid(x2[:, 1], tf2)
        tot = np.trapezoid(dp, tf2)
    elif quad_rule == "rect":
        w = np.diff(tf2)
        ga, ng = ap * np.sum(x2[:-1, 0] * w), bp * np.sum(x2[:-1, 1] * w)
        tot = np.sum(dp[:-1] * w)
    else:
        raise ValueError(f"unknown quadrature {quad_rule!r}")
    out["growth_assoc_share"] = float(ga / (ga + ng))
    out["integrated_closure_rel"] = float((ga + ng) / tot - 1)

    # Table II's own identities, and Table I against Table II's interpolated N
    kx = t2.dropna(subset=["k_per_h", "N_UOD_per_ml", "dNdt_UOD_per_ml_h"])
    out["k_identity_max_rel"] = float(np.abs(
        kx.dNdt_UOD_per_ml_h / kx.N_UOD_per_ml / kx.k_per_h - 1).max())
    qx = t2.dropna(subset=["invN_dPdt_mg_per_UOD_h", "N_UOD_per_ml", "dPdt_mg_per_ml_h"])
    out["invN_identity_max_rel"] = float(np.abs(
        qx.dPdt_mg_per_ml_h / qx.N_UOD_per_ml / qx.invN_dPdt_mg_per_UOD_h - 1).max())
    n1 = PchipInterpolator(t1.time_h.to_numpy(), t1.N_UOD_per_ml.to_numpy())
    ov = nn[(nn.time_h >= t1.time_h.min()) & (nn.time_h <= t1.time_h.max())]
    rr = np.abs(n1(ov.time_h.to_numpy()) / ov.N_UOD_per_ml.to_numpy() - 1)
    out["table1_vs_table2_N_max_rel"] = float(rr.max())
    out["table1_vs_table2_N_rms_rel"] = rms(rr)

    # Table III trends, in order of falling pH
    s3 = t3.sort_values("pH", ascending=False)
    al, be = s3.alpha_mg_lactic_per_UOD.to_numpy(), s3.beta_mg_lactic_per_UOD_h.to_numpy()
    ks = be / al
    out["beta_range_ph6_over_ph45"] = float(be[0] / be[-1])
    out["alpha_range_ph45_over_ph6"] = float(al[-1] / al[0])
    out["kstar_range_ph6_over_ph45"] = float(ks[0] / ks[-1])
    out["kstar_ph45_per_h"] = float(ks[-1])
    out["kstar_ph56_per_h"] = float(ks[1])
    out["beta_monotone_margin"] = float(np.min(-np.diff(be)))
    out["alpha_nondecreasing"] = bool(np.all(np.diff(al) >= 0))
    out["kstar_decreasing"] = bool(np.all(np.diff(ks) < 0))
    out["beta_sum_all_ph"] = float(be.sum())
    if ph_text is not None:
        out["ph_column_matches_text"] = bool(
            tuple(float(x) for x in s3.pH.to_numpy()) == tuple(float(x) for x in ph_text))
    return out
