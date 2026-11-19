"""Andrews (1968): Haldane substrate inhibition in batch, chemostat and plug flow (page J4.2).

Andrews, Biotechnol. Bioeng. 10, 707 (1968). Every function takes the
parameter dictionary ``P`` with the keys

- ``muh`` (muhat, 1/h), ``Ks`` and ``Ki`` (g/l): the inhibition function,
  eq. (1) or its rewritten form eq. (6);
- ``Y`` (yield, dimensionless);
- ``theta`` (h), ``S0`` and ``X0`` (g/l): the chemostat, eqs. (7)-(8);
- ``Si`` and ``Xi`` (g/l): the batch start, eqs. (4)-(5).

The module holds

- the growth laws and their derivative (:func:`mu`, :func:`mu_eq1`,
  :func:`mu_monod`, :func:`dmu_dS`) and the exact half-rate roots
  (:func:`half_rate_roots`);
- the batch time in closed form (:func:`t_batch`), the same by LSODA
  (:func:`batch_ivp`), the lag split of the batch time (:func:`lag_split`,
  :func:`lag_share`) and the acclimated-seed procedure (:func:`acclimate`);
- the chemostat steady states, their Jacobian, the outcome classifier and the
  saddle's stable manifold (:func:`steady_states`, :func:`jacobian`,
  :func:`fate`, :func:`separatrix_X_at`), and the fold of a seeded feed
  (:func:`fold_eq`, :func:`fold_point`, :func:`count_seeded`,
  :func:`seeded_steady_states`);
- two pymrm models: a backward-Euler chemostat (:class:`ChemostatMarch`) and
  a steady plug-flow fermenter (:class:`PlugFlowFermenter`), with the LSODA
  references they are refined against (:func:`chemostat_ivp`,
  :func:`batch_state`).

Every root is found by :func:`one_root`, a scan of a stated range with
``pymrm.checks.find_roots`` that raises unless exactly one root is found.
Every defect the page's break table injects enters through ``P`` or an
argument of these functions.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import identity as speye

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)
from pymrm.checks import find_roots

NEWTON_TOL = 1e-11    # newton step tolerance; unknowns are order-one concentrations


# --------------------------------------------------------------- roots
def one_root(f, lo, hi, what, n_scan=16, log=True, xtol=1e-16):
    """The single root of ``f`` on ``[lo, hi]``, or an error.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        The physically possible range, scanned in full.
    what : str
        Name of the root, for the error message.
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
        Unless the scan finds exactly one sign change and ``f`` could be
        evaluated at every scan point.
    """
    r = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(r["roots"]) != 1 or r["failed"]:
        raise RuntimeError(f"{what}: {r['message']} on [{lo:.6g}, {hi:.6g}]")
    return r["roots"][0]


# --------------------------------------------------------------- growth laws
def mu(S, P):
    """Andrews eq. (6), ``muhat S/(S^2/K_i + S + K_s)``.

    Parameters
    ----------
    S : float or numpy.ndarray
        Substrate concentration, g/l.
    P : dict
        Parameters (module docstring).

    Returns
    -------
    numpy.ndarray
        Specific growth rate, 1/h. Algebraically eq. (1), and finite at
        ``S = 0``.
    """
    S = np.asarray(S, float)
    return P["muh"] * S / (S * S / P["Ki"] + S + P["Ks"])


def mu_eq1(S, P):
    """Andrews eq. (1) as printed, ``muhat/(1 + K_s/S + S/K_i)``.

    Parameters
    ----------
    S : float or numpy.ndarray
        Substrate concentration, g/l.
    P : dict
        Parameters.

    Returns
    -------
    numpy.ndarray
        Specific growth rate, 1/h; ``0/0``-shaped at ``S = 0``.
    """
    S = np.asarray(S, float)
    with np.errstate(divide="ignore"):
        return P["muh"] / (1.0 + P["Ks"] / S + S / P["Ki"])


def mu_monod(S, P):
    """The ``K_i -> infinity`` limit, Monod at the same ``muhat`` and ``K_s``.

    Parameters
    ----------
    S : float or numpy.ndarray
        Substrate concentration, g/l.
    P : dict
        Parameters; ``Ki`` is not used.

    Returns
    -------
    numpy.ndarray
        Specific growth rate, 1/h.
    """
    S = np.asarray(S, float)
    return P["muh"] * S / (P["Ks"] + S)


def dmu_dS(S, P):
    """Analytic derivative of eq. (6), ``muhat (K_s - S^2/K_i)/d^2``.

    Parameters
    ----------
    S : float
        Substrate concentration, g/l.
    P : dict
        Parameters.

    Returns
    -------
    float
        ``d mu/d S``, l/(g h). Its single positive root is ``S_m``.
    """
    d = S * S / P["Ki"] + S + P["Ks"]
    return P["muh"] * (P["Ks"] - S * S / P["Ki"]) / (d * d)


def half_rate_roots(Ki, Ks):
    """Exact roots of eq. (1) equal to ``muhat/2``.

    Parameters
    ----------
    Ki, Ks : float
        Inhibition and saturation constants, g/l.

    Returns
    -------
    tuple of float
        ``K_i (1 -+ sqrt(1 - 4 K_s/K_i))/2``; they exist only for
        ``K_i >= 4 K_s``.
    """
    r = np.sqrt(1 - 4 * Ks / Ki)
    return Ki * (1 - r) / 2, Ki * (1 + r) / 2


# --------------------------------------------------------------- batch
def t_batch(Xi, Xf, Si, P, drop_constant_term=False):
    """Closed-form time for batch biomass to grow from ``Xi`` to ``Xf``.

    Constant yield makes eqs. (4)-(5) one autonomous ODE, and partial
    fractions integrate it:
    ``t = [(S_f - S_i)/K_i + (K_s/a) ln(S_i/S_f) + B ln(X_f/X_i)]/muhat``
    with ``a = S_i + X_i/Y``, ``S_f = a - X_f/Y`` and
    ``B = 1 + a/K_i + K_s/a``.

    Parameters
    ----------
    Xi, Xf, Si : float
        Initial and final biomass, initial substrate, g/l.
    P : dict
        Parameters (``muh``, ``Ks``, ``Ki``, ``Y``).
    drop_constant_term : bool
        Leave out the ``(S_f - S_i)/K_i`` term (a break-table defect).

    Returns
    -------
    float
        Time, h.
    """
    a = Si + Xi / P["Y"]
    Sf = a - Xf / P["Y"]
    B = 1.0 + a / P["Ki"] + P["Ks"] / a
    const = 0.0 if drop_constant_term else (Sf - Si) / P["Ki"]
    return (const + (P["Ks"] / a) * np.log(Si / Sf) + B * np.log(Xf / Xi)) / P["muh"]


def _batch_rhs(P):
    def rhs(t, y):
        S, X = y
        s = max(S, 0.0)
        m = P["muh"] * s / (s * s / P["Ki"] + s + P["Ks"])
        return [-m * X / P["Y"], m * X]
    return rhs


def batch_ivp(Xi, Xf, Si, P, rtol=1e-12):
    """Eqs. (4)-(5) by LSODA until ``X`` reaches ``Xf`` (a terminal event).

    Parameters
    ----------
    Xi, Xf, Si : float
        Initial and final biomass, initial substrate, g/l.
    P : dict
        Parameters.
    rtol : float
        LSODA relative tolerance; ``atol = 1e-4 rtol``.

    Returns
    -------
    float
        Time of the event, h.
    """
    hit = lambda t, y: y[1] - Xf            # noqa: E731
    hit.terminal, hit.direction = True, 1
    s = solve_ivp(_batch_rhs(P), [0, 1e5], [Si, Xi], events=hit, rtol=rtol,
                  atol=rtol * 1e-4, method="LSODA")
    if s.status < 0 or not len(s.t_events[0]):
        raise RuntimeError(f"batch never reached X = {Xf}: {s.message}")
    return float(s.t_events[0][0])


def lag_split(P, Ki=None, Xf=5.0, tangent="peak"):
    """Exact lag, exponential and tail split of the batch time, closed form.

    The lag is the intercept of the tangent to ``ln X`` at the steepest point
    of the run, ``S = S_m`` (eq. 3), or at the start of the run when
    ``S_m >= S_i``.

    Parameters
    ----------
    P : dict
        Parameters (``Si``, ``Xi`` set the run).
    Ki : float, optional
        Inhibition constant, default ``P["Ki"]``.
    Xf : float
        Target biomass, g/l.
    tangent : {"peak", "start"}
        ``"start"`` takes the tangent at ``t = 0`` (a break-table defect).

    Returns
    -------
    dict
        ``branch``, ``mu_max``, ``S_m``, ``X_peak``, ``t_peak``, ``lag``,
        ``expo``, ``tail`` and ``t_f``.
    """
    Q = dict(P, Ki=P["Ki"] if Ki is None else Ki)
    Si, Xi = Q["Si"], Q["Xi"]
    a = Si + Xi / Q["Y"]
    Sm, mum = np.sqrt(Q["Ks"] * Q["Ki"]), Q["muh"] / (1 + 2 * np.sqrt(Q["Ks"] / Q["Ki"]))
    branch = "interior peak"
    if Sm >= Si or tangent == "start":
        branch = "peak beyond S_i" if Sm >= Si else "tangent at t = 0"
        Sm, mum = Si, float(mu(Si, Q))
    Xs = Q["Y"] * (a - Sm)
    ts = t_batch(Xi, Xs, Si, Q)
    lam = ts - np.log(Xs / Xi) / mum
    tf = t_batch(Xi, Xf, Si, Q)
    expo = np.log(Xf / Xi) / mum
    return dict(Ki=Q["Ki"], branch=branch, mu_max=mum, S_m=Sm, X_peak=Xs, t_peak=ts,
                lag=lam, expo=expo, tail=tf - lam - expo, t_f=tf)


def lag_share(P, Ki=None, Xf=5.0, base_Ki=1e12, tangent="peak"):
    """Shares of the extra batch time that are lag, exponential and tail.

    Parameters
    ----------
    P : dict
        Parameters.
    Ki : float, optional
        Inhibition constant of the run, default ``P["Ki"]``.
    Xf : float
        Target biomass, g/l.
    base_Ki : float
        Inhibition constant of the no-inhibition control.
    tangent : {"peak", "start"}
        Passed to :func:`lag_split`.

    Returns
    -------
    dict
        ``dtf`` (extra time, h) and the shares ``lag``, ``expo``, ``tail``.
    """
    r = lag_split(P, Ki, Xf=Xf, tangent=tangent)
    b = lag_split(P, base_Ki, Xf=Xf, tangent=tangent)
    d = r["t_f"] - b["t_f"]
    return dict(dtf=d, lag=(r["lag"] - b["lag"]) / d, expo=(r["expo"] - b["expo"]) / d,
                tail=(r["tail"] - b["tail"]) / d)


def acclimate(mode, P, target=5.0, S_add=2.0, S_stop=0.02):
    """Andrews' acclimated seed, book p. 715: organisms kept, substrate topped up.

    Parameters
    ----------
    mode : {"increment", "reset", "ladder"}
        The three literal readings of the sentence: add ``S_add`` to what is
        left (the reading Fig. 7's annotation prints), refill to ``S_add``, or
        step the initial concentration up ``S_add, 2 S_add, ...``.
    P : dict
        Parameters (``Xi`` is the first inoculum).
    target : float
        Biomass to reach, g/l.
    S_add, S_stop : float
        Substrate added per cycle and left at its end, g/l.

    Returns
    -------
    tuple
        Total time (h), number of cycles, list of cycle times (h).
    """
    X, t, S_start, times = P["Xi"], 0.0, S_add, []
    for k in range(60):
        a = S_start + X / P["Y"]
        X_end = P["Y"] * (a - S_stop)
        if X_end >= target:
            dt = t_batch(X, target, S_start, P)
            times.append(dt)
            return t + dt, k + 1, times
        dt = t_batch(X, X_end, S_start, P)
        t += dt
        times.append(dt)
        X = X_end
        S_start = {"increment": S_stop + S_add, "reset": S_add,
                   "ladder": S_start + S_add}[mode]
    raise RuntimeError("acclimation never reached the target")


# --------------------------------------------------------------- chemostat
def steady_states(P, S0=None):
    """Nontrivial steady states of eqs. (7)-(8) with a sterile feed.

    The roots of ``mu(S) = 1/theta``: one on each limb of ``mu``, each found
    by :func:`one_root` on the whole limb (``mu`` is monotone there).

    Parameters
    ----------
    P : dict
        Parameters.
    S0 : float, optional
        Feed substrate, default ``P["S0"]``.

    Returns
    -------
    list of float
        ``[]`` past the fold, ``[S_low]`` when the upper root lies above the
        feed, else ``[S_low, S_high]``.
    """
    S0 = P["S0"] if S0 is None else S0
    Sm = np.sqrt(P["Ks"] * P["Ki"])
    f = lambda s: float(mu(s, P)) - 1 / P["theta"]      # noqa: E731
    if f(Sm) <= 0:
        return []
    lo = one_root(f, 1e-14, Sm, "lower steady state", xtol=1e-17)
    hi_cap = max(S0, Sm * (1 + 1e-12))
    if f(hi_cap) > 0:
        return [lo]
    hi = one_root(f, Sm, hi_cap, "upper steady state", xtol=1e-17)
    return [lo, hi]


def jacobian(S, X, P):
    """Jacobian of eqs. (7)-(8) in ``(X, S)``, with the analytic ``dmu/dS``.

    Parameters
    ----------
    S, X : float
        State, g/l.
    P : dict
        Parameters.

    Returns
    -------
    numpy.ndarray
        2 x 2 Jacobian, 1/h.
    """
    dmu = dmu_dS(S, P)
    m = float(mu(S, P))
    th, Y = P["theta"], P["Y"]
    return np.array([[-1 / th + m, dmu * X],
                     [-m / Y, -1 / th - dmu * X / Y]])


def fate(X0_, S0_ini, S0_feed, P, tend=800.0, rtol=1e-12, event=True):
    """Outcome of a chemostat run: +1 it reaches the operating state, -1 it washes out.

    Parameters
    ----------
    X0_, S0_ini : float
        Initial biomass and substrate, g/l.
    S0_feed : float or callable
        Feed substrate, g/l, or a function of time (h).
    P : dict
        Parameters (``theta``, ``Y``, ``muh``, ``Ks``, ``Ki``).
    tend : float
        Integration limit, h.
    rtol : float
        LSODA tolerance; ``atol = 1e-3 rtol``.
    event : bool
        Terminate on reaching 1.001 times the low steady substrate or on
        biomass falling through 1e-9 g/l. ``False`` classifies on the state
        at ``tend`` instead (a break-table defect).

    Returns
    -------
    int
        +1 or -1.
    """
    ss = steady_states(P, S0_feed if not callable(S0_feed) else S0_feed(tend))
    if not ss:
        return -1                     # past the fold: washout is the only state
    s_lo = ss[0]
    th, Y = P["theta"], P["Y"]

    def rhs(t, y):
        X, S = y
        S0t = S0_feed(t) if callable(S0_feed) else S0_feed
        m = float(mu(max(S, 0.0), P))
        return [-X / th + m * X, (S0t - S) / th - m * X / Y]

    if not event:
        s = solve_ivp(rhs, [0, tend], [X0_, S0_ini], rtol=rtol, atol=rtol * 1e-3,
                      method="LSODA")
        if s.status < 0:
            raise RuntimeError(f"chemostat integration failed: {s.message}")
        return 1 if s.y[1, -1] < s_lo * 2 else -1
    won = lambda t, y: y[1] - s_lo * 1.001         # noqa: E731
    won.terminal, won.direction = True, -1
    lost = lambda t, y: y[0] - 1e-9              # noqa: E731
    lost.terminal, lost.direction = True, -1
    s = solve_ivp(rhs, [0, tend], [X0_, S0_ini], events=[won, lost],
                  rtol=rtol, atol=rtol * 1e-3, method="LSODA")
    if s.status < 0:
        raise RuntimeError(f"chemostat integration failed: {s.message}")
    if len(s.t_events[0]):
        return 1
    if len(s.t_events[1]):
        return -1
    return 1 if s.y[1, -1] < s_lo * 2 else -1


def separatrix_X_at(S_query, S0_feed, P, eps=1e-8):
    """The saddle's stable manifold, integrated backward in time to ``S = S_query``.

    Shares only the right-hand side with :func:`fate`: no bisection, no
    outcome classifier, no event on ``X``.

    Parameters
    ----------
    S_query : float
        Substrate at which the manifold is read, g/l.
    S0_feed : float
        Feed substrate, g/l.
    P : dict
        Parameters.
    eps : float
        Offset from the saddle along the stable eigenvector.

    Returns
    -------
    float
        Biomass on the manifold at ``S_query``, g/l.
    """
    Sm = np.sqrt(P["Ks"] * P["Ki"])
    th, Y = P["theta"], P["Y"]
    s_hi = one_root(lambda s: float(mu(s, P)) - 1 / th, Sm, max(S0_feed, 1e4),
                    "saddle", xtol=1e-17)
    x_hi = Y * (S0_feed - s_hi)
    dmu = dmu_dS(s_hi, P)
    J = np.array([[0.0, dmu * x_hi],
                  [-float(mu(s_hi, P)) / Y, -1 / th - dmu * x_hi / Y]])
    w, v = np.linalg.eig(J)
    k_stable = 0 if w[0].real < w[1].real else 1  # the saddle's negative eigenvalue
    e = v[:, k_stable].real
    e = e / np.linalg.norm(e)
    if e[1] > 0:
        e = -e                                     # head down in S

    def rhs_back(t, y):                            # time reversed
        X, S = y
        m = float(mu(max(S, 0.0), P))
        return [X / th - m * X, -((S0_feed - S) / th - m * X / Y)]

    hit = lambda t, y: y[1] - S_query              # noqa: E731
    hit.terminal, hit.direction = True, -1
    s = solve_ivp(rhs_back, [0, 600], np.array([x_hi, s_hi]) + eps * e, events=[hit],
                  rtol=1e-13, atol=1e-18, method="LSODA")
    if s.status < 0 or not len(s.t_events[0]):
        raise RuntimeError("the stable manifold did not reach the query substrate")
    return float(s.y_events[0][0][0])


def fold_eq(S, P):
    """The fold of a seeded feed as one equation, ``mu (theta mu - 1) = mu' (S_0 - S)``.

    ``G(S) = (S_0 - S)(1 - theta mu) - theta mu X_0/Y`` has a double root
    there (``G = G' = 0``); eliminating ``X_0`` leaves this equation.

    Parameters
    ----------
    S : float
        Substrate, g/l.
    P : dict
        Parameters.

    Returns
    -------
    float
        Residual, 1/h^2.
    """
    m = float(mu(S, P))
    return m * (P["theta"] * m - 1) - dmu_dS(S, P) * (P["S0"] - S)


def fold_point(P):
    """Substrate at the fold and the feed biomass that puts the fold there.

    Parameters
    ----------
    P : dict
        Parameters.

    Returns
    -------
    tuple of float
        ``(S_fold, X0_fold)``, g/l: the root of :func:`fold_eq` on the upper
        limb ``[S_m, S_0]`` and ``Y (S_0 - S)(1 - theta mu)/(theta mu)``.
    """
    Sm = np.sqrt(P["Ks"] * P["Ki"])
    if Sm >= P["S0"]:
        raise RuntimeError("no upper limb below the feed: the fold does not exist")
    S = one_root(lambda s: fold_eq(s, P), Sm, P["S0"], "fold", log=False)
    m = float(mu(S, P))
    return S, P["Y"] * (P["S0"] - S) * (1 - P["theta"] * m) / (P["theta"] * m)


def _G(S, X0f, P):
    m = mu(S, P)
    return (P["S0"] - S) * (1 - P["theta"] * m) - P["theta"] * m * X0f / P["Y"]


def count_seeded(X0f, P, ngrid):
    """Number of steady states with feed biomass ``X0f``, counted on a scan.

    Counts the sign changes of the pole-free ``G(S)`` on ``ngrid`` geometric
    points of ``(0, S_0)``. Grid-limited by construction: two roots closer than
    one scan cell are not seen.

    Parameters
    ----------
    X0f : float
        Feed biomass, g/l.
    P : dict
        Parameters.
    ngrid : int
        Scan points.

    Returns
    -------
    int
        Number of sign changes.
    """
    grid = np.geomspace(1e-10, P["S0"] * (1 - 1e-12), ngrid)
    v = np.asarray(_G(grid, X0f, P), float)
    return int(np.count_nonzero(v[:-1] * v[1:] < 0))


def seeded_steady_states(X0f, P, n_scan=6001):
    """Steady substrate concentrations with feed biomass ``X0f``, root-found.

    Parameters
    ----------
    X0f : float
        Feed biomass, g/l (positive).
    P : dict
        Parameters.
    n_scan : int
        Scan points of ``find_roots`` on a geometric grid of ``(0, S_0)``.

    Returns
    -------
    list of float
        Every root of ``G``; their number is the result.
    """
    r = find_roots(lambda s: float(_G(s, X0f, P)), 1e-10, P["S0"] * (1 - 1e-12),
                   n_scan=n_scan, log=True, xtol=1e-16)
    if r["failed"]:
        raise RuntimeError(r["message"])
    return r["roots"]


def chemostat_ivp(c0, T, P, rtol=1e-12):
    """Eqs. (7)-(8) by LSODA over ``[0, T]``, the marcher's reference.

    Parameters
    ----------
    c0 : sequence of float
        Initial ``(S, X)``, g/l.
    T : float
        End time, h.
    P : dict
        Parameters.
    rtol : float
        LSODA tolerance; ``atol = 1e-3 rtol``.

    Returns
    -------
    numpy.ndarray
        ``(S, X)`` at ``T``.
    """
    th, Y = P["theta"], P["Y"]

    def rhs(t, y):
        S, X = y
        m = float(mu(max(S, 0.0), P))
        return [(P["S0"] - S) / th - m * X / Y, (P["X0"] - X) / th + m * X]
    s = solve_ivp(rhs, [0, T], c0, rtol=rtol, atol=rtol * 1e-3, method="LSODA")
    if s.status != 0:
        raise RuntimeError(f"chemostat reference failed: {s.message}")
    return s.y[:, -1]


def batch_state(S0i, X0i, tau, P, rtol=1e-13):
    """Batch state after ``tau`` hours by LSODA, the plug-flow outlet's reference.

    Parameters
    ----------
    S0i, X0i : float
        Initial substrate and biomass, g/l.
    tau : float
        Time, h.
    P : dict
        Parameters.
    rtol : float
        LSODA relative tolerance; ``atol = 1e-16``.

    Returns
    -------
    numpy.ndarray
        ``(S, X)`` at ``tau``.
    """
    s = solve_ivp(_batch_rhs(P), [0, tau], [S0i, X0i], rtol=rtol, atol=1e-16,
                  method="LSODA")
    if s.status != 0:
        raise RuntimeError(f"batch reference failed: {s.message}")
    return s.y[:, -1]


# --------------------------------------------------------------- pymrm models
def _growth_source(c, P):
    """Growth terms of ``(S, X)`` on a ``(..., 2)`` field: ``(-mu X/Y, mu X)``."""
    S = np.clip(c[..., 0:1], 0.0, None)
    m = P["muh"] * S / (S * S / P["Ki"] + S + P["Ks"]) * c[..., 1:2]
    return np.concatenate([-m / P["Y"], m], axis=-1)


class ChemostatMarch:
    """Backward Euler on Andrews eqs. (7)-(8) with pymrm's ``newton`` + ``NumJac``.

    The state has shape ``(1, 2)``: a completely mixed reactor is one cell and
    two fields ``(S, X)``, and the growth term is pointwise in the fields,
    which is ``NumJac``'s default stencil. The dilution terms are linear, so
    their operator (an identity) is assembled once in the constructor.

    Parameters
    ----------
    P : dict
        Parameters (``theta``, ``S0``, ``X0`` and the kinetics).
    nstep : int
        Number of time steps.
    tol : float
        ``newton`` step tolerance.
    dilution : bool
        Keep the ``-c/theta`` outflow term; ``False`` drops it (a break-table
        defect).
    theta_t : float
        Weight of the new time level: 1 is backward Euler, 1/2 the trapezoidal
        rule (a break-table defect that changes the order).
    """

    def __init__(self, P, nstep=600, tol=NEWTON_TOL, dilution=True, theta_t=1.0):
        self.P, self.nstep, self.tol, self.dilution = P, nstep, tol, dilution
        self.theta_t = theta_t
        self.shape = (1, 2)                       # (cells, fields) = (S, X)
        self.numjac = NumJac(self.shape)          # last axis coupled in full
        self.eye = speye(2, format="csc")         # constant, assembled once
        self.inflow = np.array([[P["S0"], P["X0"]]]) / P["theta"]

    def source(self, c):
        """Growth terms; the linear dilution is in the constant operator."""
        return _growth_source(c, self.P)

    def step(self, c_old, h):
        """One backward-Euler step of length ``h``.

        Parameters
        ----------
        c_old : numpy.ndarray
            State at the start of the step, shape ``(1, 2)``.
        h : float
            Step, h.

        Returns
        -------
        numpy.ndarray
            State at the end of the step.
        """
        out_rate = 1.0 / self.P["theta"] if self.dilution else 0.0
        w = self.theta_t
        rate_old = (self.inflow.reshape(-1, 1) - out_rate * c_old.reshape(-1, 1)
                    + self.source(c_old).reshape(-1, 1))

        def res(c):
            s, js = self.numjac(self.source, c)
            rate = self.inflow.reshape(-1, 1) - out_rate * c.reshape(-1, 1) + s.reshape(-1, 1)
            r = (c.reshape(-1, 1) - c_old.reshape(-1, 1)) / h - w * rate - (1 - w) * rate_old
            return r, self.eye / h + w * (out_rate * self.eye - js)
        out = newton(res, c_old, maxfev=100, tol=self.tol, stagnation=True)
        if not out.success:
            raise RuntimeError(f"chemostat step did not converge: {out.message}")
        return out.x.reshape(self.shape)

    def march(self, c0, T):
        """March from ``c0`` over ``[0, T]``.

        Parameters
        ----------
        c0 : sequence of float
            Initial ``(S, X)``, g/l.
        T : float
            End time, h.

        Returns
        -------
        numpy.ndarray
            State at ``T``, shape ``(1, 2)``.
        """
        c = np.asarray(c0, float).reshape(self.shape)
        h = T / self.nstep
        for _ in range(self.nstep):
            c = self.step(c, h)
        return c


class PlugFlowFermenter:
    """Steady one-dimensional plug flow with Andrews' kinetics, fields ``(S, X)``.

    ``v dc/dz = r(c)`` on ``0 < z < v tau``: first-order upwind convection
    (``construct_convflux_upwind``), ``construct_div`` with ``nu = 0``
    (Cartesian), the growth term through ``NumJac((n, 2))`` and ``newton``.
    The constant operators are assembled once in the constructor.

    Parameters
    ----------
    Sf, Xf : float
        Feed substrate and biomass, g/l.
    tau : float
        Residence time, h.
    P : dict
        Kinetic parameters.
    ncell : int
        Number of cells.
    v : float
        Velocity (only ``v tau`` matters).
    nu : int
        Geometry of ``construct_div``; 1 (cylindrical) is a break-table defect.
    inlet_bc : {"dirichlet", "neumann"}
        Inlet condition; ``"neumann"`` (zero gradient) is a break-table defect.
    init : {"feed", "batch"}
        Initial guess: the feed in every cell, or the batch trajectory sampled
        at ``z/v`` (a per-``tau`` guess with no continuation).
    tol : float
        ``newton`` step tolerance.
    """

    def __init__(self, Sf, Xf, tau, P, ncell=400, v=1.0, nu=0, inlet_bc="dirichlet",
                 init="feed", tol=NEWTON_TOL):
        self.P, self.shape, self.v, self.tau, self.tol = P, (ncell, 2), v, tau, tol
        self.z_f = np.linspace(0.0, v * tau, ncell + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        d0 = np.array([[Sf, Xf]])
        # Outward normal, a dc/dn + b c = d:
        #   z = 0 inlet:  c = (S_f, X_f), Dirichlet {a: 0, b: 1, d: feed}
        #   z = L outlet: convection only, the face carries the last cell value,
        #                 {"outflow": True}
        left = ({"a": 0.0, "b": 1.0, "d": d0} if inlet_bc == "dirichlet"
                else {"a": 1.0, "b": 0.0, "d": 0.0})     # dc/dz = 0 at the inlet (defect)
        self.bc = (left, {"outflow": True})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                  self.bc, v=v, axis=0)
        self.div = construct_div(self.shape, self.z_f, nu=nu, axis=0)  # nu=0 Cartesian
        self.A, self.b = self.div @ conv, self.div @ conv_bc
        self.numjac = NumJac(self.shape)          # pointwise in the 2 fields
        if init == "feed":
            self.c0 = np.tile(d0, (ncell, 1))
        else:
            sol = solve_ivp(_batch_rhs(P), [0, tau], [Sf, Xf], t_eval=self.z_c / v,
                            rtol=1e-10, atol=1e-14, method="LSODA")
            if sol.status != 0:
                raise RuntimeError(f"batch initial guess failed: {sol.message}")
            self.c0 = np.ascontiguousarray(sol.y.T)

    def source(self, c):
        """Growth terms of ``(S, X)``."""
        return _growth_source(c, self.P)

    def solve(self):
        """Steady solve; raises unless ``newton`` reports success.

        Returns
        -------
        PlugFlowFermenter
            ``self``, with the solution in ``self.c``, shape ``(ncell, 2)``.
        """
        def resid(c):
            s, js = self.numjac(self.source, c)
            return (self.b + self.A @ c.reshape((-1, 1)) - s.reshape((-1, 1)),
                    self.A - js)
        # bounds: S and X are concentrations and stay non-negative; S -> 0 at long tau
        r = newton(resid, self.c0, maxfev=300, tol=self.tol, stagnation=True, bounds=0.0)
        if not r.success:
            raise RuntimeError(f"plug-flow solve did not converge: {r.message}")
        self.c = r.x.reshape(self.shape)
        return self

    def outlet(self):
        """Outlet ``(S, X)`` from ``compute_boundary_values`` under the same bc.

        Returns
        -------
        numpy.ndarray
            ``(S, X)`` on the outlet face, g/l.
        """
        return np.asarray(compute_boundary_values(self.c, self.z_f, self.z_c,
                                                  self.bc, axis=0)[2]).reshape(-1)[:2]
