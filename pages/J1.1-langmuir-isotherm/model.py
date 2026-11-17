"""Langmuir (1918): the isotherm forms, the fits and the rate equation (page J1.1).

Langmuir, J. Am. Chem. Soc. 40, 1361 (1918). Everything here is a closed form
in one variable or a march of one ordinary differential equation:

- the isotherm forms the paper fits, eq. (31) (:func:`q_eq31`), eq. (33)
  (:func:`q_eq33`), eq. (19) (:func:`q_eq19`), and Freundlich's power law
  (:func:`q_freundlich`);
- Case I and Case VI in BET variables, ``x = sigma_2 mu``, ``c = sigma_1/sigma_2``
  (:func:`theta_caseI`, :func:`theta_caseVI`, :func:`departure`), and the
  occupancy ladder of eqs. (22)/(23) summed term by term (:func:`caseVI_series`);
- least squares with the linear amplitude profiled out: the stationarity condition
  root-found by ``newton`` and by Brent (:func:`profile_fit`, :func:`fit_langmuir`,
  :func:`fit_freundlich`), and a derivative-free scan plus Brent refinement
  (:func:`profile_rms`);
- eq. (14) marched in time (:func:`march_eq14`), the two branches of eq. (15)
  and the half-coverage time read off a marched trajectory
  (:func:`t_half_of_march`);
- the scatter-to-separation ratio of one table with its recovered ``p0`` scaled
  (:func:`sos_at_p0_factor`).
"""
import numpy as np
from scipy.optimize import minimize_scalar
from pymrm import NumJac, newton
from pymrm.checks import find_roots


# ---------------------------------------------------------- the isotherm forms
def q_eq31(p, a, b):
    """Eq. (31), the Langmuir isotherm in the paper's variables.

    Parameters
    ----------
    p : float or ndarray
        Pressure in bars (dyn/cm2).
    a, b : float
        The fitted constants of the table header; ``b`` in cu. mm.

    Returns
    -------
    float or ndarray
        Adsorbed amount ``a b p / (1 + a p)`` in cu. mm.
    """
    return a * b * p / (1.0 + a * p)


def q_eq33(p, b1, a2, b2):
    """Eq. (33): a site type saturated below the lowest pressure, plus eq. (31).

    Parameters
    ----------
    p : float or ndarray
        Pressure in bars.
    b1 : float
        Amount on the saturated sites, cu. mm.
    a2, b2 : float
        Constants of the second, eq. (31) term.

    Returns
    -------
    float or ndarray
        ``b1 + a2 b2 p / (1 + a2 p)``.
    """
    return b1 + a2 * b2 * p / (1.0 + a2 * p)


def q_eq19(p, pairs):
    """Eq. (19), Case II: a sum of Langmuir terms.

    Parameters
    ----------
    p : float or ndarray
        Pressure in bars.
    pairs : iterable of (float, float)
        ``(a, b)`` of each term.

    Returns
    -------
    float or ndarray
        The summed adsorbed amount.
    """
    return sum(a * b * p / (1.0 + a * p) for a, b in pairs)


def q_freundlich(p, amp, n):
    """Freundlich's ``q = A p^n``, the rival Langmuir prints in Table VII.

    Parameters
    ----------
    p : float or ndarray
        Pressure in bars.
    amp, n : float
        Amplitude ``A`` and exponent ``n``.

    Returns
    -------
    float or ndarray
        ``A p^n``.
    """
    return amp * p ** n


def lang_basis(p, a):
    """Shape of eq. (31) with the amplitude ``a b`` taken out: ``p / (1 + a p)``."""
    return p / (1.0 + a * p)


def freu_basis(p, n):
    """Shape of Freundlich's law with the amplitude taken out: ``p^n``."""
    return p ** n


def rms(r):
    """Root mean square of a residual vector."""
    return float(np.sqrt(np.mean(np.asarray(r, float) ** 2)))


# ------------------------------------------------- Case I and Case VI (BET)
def theta_caseI(x, c):
    """Case I in BET variables, ``x = sigma_2 mu`` and ``c = sigma_1/sigma_2``.

    Parameters
    ----------
    x, c : float or ndarray
        Relative pressure and the ratio of lives.

    Returns
    -------
    float or ndarray
        ``c x / (1 + c x)``.
    """
    return c * x / (1.0 + c * x)


def theta_caseVI(x, c):
    """Case VI, eqs. (29) and (30) with ``sigma_3 = sigma_4 = ... = sigma_2``.

    Written from the coefficients of eq. (30), not from the factored BET form,
    so that the identity with BET is a statement about the paper's own algebra.

    Parameters
    ----------
    x, c : float or ndarray
        Relative pressure and the ratio of lives.

    Returns
    -------
    float or ndarray
        Amount adsorbed in monolayers.
    """
    a30 = c - 2.0                      # (sigma_1 - 2 sigma_2)/sigma_2
    b30 = 4.0 - 3.0 - c                # sigma_2(4 sigma_2 - 3 sigma_2 - sigma_1)/sigma_2^2
    return c / (1.0 / x + a30 + b30 * x)


def departure(x, c):
    """Relative departure of Case VI from Case I, ``theta_VI/theta_I - 1``.

    Parameters
    ----------
    x, c : float or ndarray
        Relative pressure and the ratio of lives.

    Returns
    -------
    float or ndarray
        ``(1 + c x) / ((1 - x)(1 + (c - 1) x)) - 1``, strictly increasing in ``x``.
    """
    return (1.0 + c * x) / ((1.0 - x) * (1.0 + (c - 1.0) * x)) - 1.0


def caseVI_series(x, c, n_layers, bare=True):
    """Eqs. (22)/(23) summed term by term with ``sigma_3.. = sigma_2``.

    No closed form is used: the occupancies are ``theta_0 = theta_1/(c x)``
    and ``theta_k = theta_1 x^(k-1)``, summed to ``n_layers``.

    Parameters
    ----------
    x, c : float
        Relative pressure and the ratio of lives.
    n_layers : int
        Number of layers kept.
    bare : bool
        Count the bare fraction ``theta_0`` in the inventory of eq. (22);
        ``False`` is the defect of the break row that drops it.

    Returns
    -------
    float
        Amount adsorbed in monolayers.
    """
    th = [1.0 / (c * x)] + [np.prod([1.0] + [x] * (k - 1)) for k in range(1, n_layers + 1)]
    inv = sum(th) if bare else sum(th[1:])
    tot = sum(k * th[k] for k in range(1, n_layers + 1))
    return tot / inv


# ---------------------------------------------------------------- the fits
def profile_fit(basis, p, q, t0, lo, hi, h_rel=1e-6, newton_tol=1e-12,
                backward_tol=1e-8, xtol=1e-13, log=None, label=""):
    """Least squares of ``amp * basis(p, t)`` on ``q``, the amplitude profiled out.

    The amplitude enters linearly and is eliminated exactly, so one parameter
    ``t`` is left. Its stationarity condition ``dS/dt = 0`` (central difference of
    step ``h_rel * max(|t|, 1)``) is root-found twice: by ``newton`` on a
    ``NumJac((1, 1))`` Jacobian and by Brent on ``[lo, hi]``, after a root scan of the
    bracket that raises unless it holds exactly one root. Both get the same
    ``dS``, so their gap is a root-finder cross-check, not an independent fit.
    ``newton`` runs with ``stagnation=True``: a central-difference ``dS`` has a
    noise floor, and on two tables the step settles there without meeting
    ``newton_tol``.

    Parameters
    ----------
    basis : callable
        ``basis(p, t)``, the shape with the amplitude taken out.
    p, q : ndarray
        Pressures and adsorbed amounts.
    t0 : float
        Newton's starting value.
    lo, hi : float
        Brent's bracket.
    h_rel : float
        Relative central-difference step of ``dS``.
    newton_tol, backward_tol : float
        ``newton`` step and backward-error tolerances.
    xtol : float
        Brent's absolute tolerance.
    log : list, optional
        ``(label, |t_newton / t_brent - 1|)`` is appended to it.
    label : str
        Name of the fit in ``log``.

    Returns
    -------
    t, amp, rms : float
        Brent's shape parameter, the amplitude there and the RMS residual.

    Raises
    ------
    RuntimeError
        When ``newton`` does not report success, or ``dS`` has not exactly one root on
        ``[lo, hi]``.
    """
    def amp(t):
        b = basis(p, t)
        return float(np.dot(q, b) / np.dot(b, b))

    def S(t):
        return float(np.sum((amp(t) * basis(p, t) - q) ** 2))

    def dS(t):
        h = h_rel * max(abs(t), 1.0)
        return (S(t + h) - S(t - h)) / (2.0 * h)

    jac = NumJac((1, 1))
    res = newton(lambda z: jac(lambda zz: np.array([[dS(float(zz[0, 0]))]]), z),
                 np.array([[float(t0)]]), tol=newton_tol, backward_tol=backward_tol,
                 stagnation=True, maxfev=120)
    if not res.success:
        raise RuntimeError(f"newton failed on the fit {label!r}: {res.message}")
    t_newton = float(np.asarray(res.x, float).reshape(1)[0])
    t_brent = unique_root(dS, lo, hi, log=True, xtol=xtol)
    if log is not None:
        log.append((label, abs(t_newton / t_brent - 1.0)))
    return t_brent, amp(t_brent), rms(amp(t_brent) * basis(p, t_brent) - q)


def fit_langmuir(p, q, label="", **kw):
    """Best ``(a, b)`` of eq. (31) on ``q`` itself, not on the p/q linearisation.

    Parameters
    ----------
    p, q : ndarray
        Pressures and adsorbed amounts.
    label : str
        Name of the fit for the two-route log.
    **kw
        Passed to :func:`profile_fit`.

    Returns
    -------
    a, b, rms : float
    """
    slope, icept = np.polyfit(p, p / q, 1)
    a, ab, r = profile_fit(lang_basis, p, q, max(1e-4, slope / icept), 1e-6, 50.0,
                           label="langmuir " + label, **kw)
    return a, ab / a, r


def fit_freundlich(p, q, label="", **kw):
    """Best ``(A, n)`` of ``q = A p^n``, so the rival is beaten at its best.

    Parameters
    ----------
    p, q : ndarray
        Pressures and adsorbed amounts.
    label : str
        Name of the fit for the two-route log.
    **kw
        Passed to :func:`profile_fit`.

    Returns
    -------
    amp, n, rms : float
    """
    n0 = np.polyfit(np.log(p), np.log(q), 1)[0]
    n, amp, r = profile_fit(freu_basis, p, q, n0, 0.02, 2.0,
                            label="freundlich " + label, **kw)
    return amp, n, r


def profile_rms(basis, p, q, tlo, thi, n_grid=4000, xtol=1e-13, refine=True):
    """Best RMS of ``amp * basis(p, t)``, located without any derivative.

    The shape parameter is scanned on a log grid of ``n_grid`` points over
    ``[tlo, thi]`` and the grid minimum refined by Brent's bracketed
    minimisation in ``ln t``. Shares no code with :func:`profile_fit`.

    Parameters
    ----------
    basis : callable
        ``basis(p, t)``.
    p, q : ndarray
        Pressures and adsorbed amounts.
    tlo, thi : float
        The scanned window.
    n_grid : int
        Scan points.
    xtol : float
        Brent's tolerance in ``ln t``.
    refine : bool
        ``False`` returns the grid minimum, the defect of a break row.

    Returns
    -------
    rms, t, amp, rms_grid : float
        The refined optimum, its parameter and amplitude, and the grid minimum.

    Raises
    ------
    RuntimeError
        When the grid minimum sits on an end of the window.
    """
    ts = np.logspace(np.log10(tlo), np.log10(thi), n_grid)

    def obj(lt):
        b = basis(p, float(np.exp(lt)))
        return rms(float(np.dot(q, b) / np.dot(b, b)) * b - q)

    vals = np.array([obj(np.log(t)) for t in ts])
    i = int(np.argmin(vals))
    if not 0 < i < n_grid - 1:
        raise RuntimeError(f"optimum at or beyond the end of [{tlo:g}, {thi:g}]: widen it")
    if refine:
        r = minimize_scalar(obj, bracket=(np.log(ts[i - 1]), np.log(ts[i]), np.log(ts[i + 1])),
                            method="brent", options={"xtol": xtol})
        t_best, fun = float(np.exp(r.x)), float(r.fun)
    else:
        t_best, fun = float(ts[i]), float(vals[i])
    b = basis(p, t_best)
    return fun, t_best, float(np.dot(q, b) / np.dot(b, b)), float(vals.min())


# ------------------------------------------------------- the rate equation
def march_eq14(theta1, k_rate, t_end, n_steps, scheme="implicit", newton_tol=1e-14):
    """March eq. (14), ``dtheta'/dt = k (theta1 - theta')``, from ``theta' = 0``.

    ``k = (N/N0)(nu1 + alpha mu)`` and ``theta1 = alpha mu/(nu1 + alpha mu)``.
    The implicit Euler step is solved by ``newton`` on a ``NumJac((1, 1))``
    Jacobian, so the march shares no algebra with any closed form.

    Parameters
    ----------
    theta1 : float
        Equilibrium coverage.
    k_rate : float
        Rate constant of the approach to equilibrium.
    t_end : float
        End time.
    n_steps : int
        Number of equal steps.
    scheme : {"implicit", "midpoint"}
        Implicit Euler (first order), or the explicit midpoint rule (second
        order), the defect of a break row.
    newton_tol : float
        ``newton`` step tolerance of each implicit step.

    Returns
    -------
    t, theta : ndarray
        The ``n_steps + 1`` step times and coverages.

    Raises
    ------
    RuntimeError
        When a ``newton`` step does not report success.
    """
    dt = t_end / n_steps
    th = np.zeros((1, 1))
    jac = NumJac((1, 1))
    out = [0.0]
    for i in range(n_steps):
        if scheme == "midpoint":
            half = th + 0.5 * dt * k_rate * (theta1 - th)
            th = th + dt * k_rate * (theta1 - half)
        else:
            prev = th.copy()

            def res(z, prev=prev):
                return (z - prev) / dt - k_rate * (theta1 - z)

            sol = newton(lambda z: jac(res, z), th, tol=newton_tol, maxfev=60)
            if not sol.success:
                raise RuntimeError(f"implicit step {i + 1} of {n_steps} failed: {sol.message}")
            th = np.asarray(sol.x, float).reshape(1, 1)
        out.append(float(th[0, 0]))
    return np.linspace(0.0, t_end, n_steps + 1), np.array(out)


def t_eq15_printed(thp, theta1, k_rate):
    """Eq. (15) as printed: ``t = ln(theta' / (theta1 - theta')) / k``."""
    return np.log(thp / (theta1 - thp)) / k_rate


def t_eq15_corrected(thp, theta1, k_rate):
    """Eq. (15) with ``theta1`` in the numerator: ``t = ln(theta1 / (theta1 - theta')) / k``."""
    return np.log(theta1 / (theta1 - thp)) / k_rate


def unique_root(f, lo, hi, log=False, n_scan=64, xtol=1e-16):
    """The one root of ``f`` on ``[lo, hi]``, by ``pymrm.checks.find_roots``.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        The physical range scanned.
    log : bool
        Scan on a log grid.
    n_scan : int
        Scan points.
    xtol : float
        Brent tolerance of the refined root.

    Returns
    -------
    float

    Raises
    ------
    RuntimeError
        Unless exactly one root is found.
    """
    found = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(found["roots"]) != 1:
        raise RuntimeError(f"{len(found['roots'])} roots on [{lo:g}, {hi:g}]: "
                           f"{found['message']}")
    return float(found["roots"][0])


def t_half_of_march(t, th, theta1):
    """Time at which a marched trajectory reaches ``theta1 / 2``.

    The root of the piecewise-linear interpolant of the marched samples,
    located over the whole march by :func:`unique_root`; nothing analytic enters.

    Parameters
    ----------
    t, th : ndarray
        Step times and coverages of the march.
    theta1 : float
        Equilibrium coverage.

    Returns
    -------
    float
    """
    return unique_root(lambda tt: float(np.interp(tt, t, th)) - 0.5 * theta1,
                       float(t[0]), float(t[-1]), n_scan=4 * len(t))


# --------------------------------------------- the p0 a sceptic would attack
def sos_at_p0_factor(f, p, q, a, b, p0_ref):
    """Scatter over separation for one table when its recovered ``p0`` is scaled by ``f``.

    ``x = p / (f p0)`` and ``c = a f p0``, both as eq. (39) requires.

    Parameters
    ----------
    f : float
        Factor on ``p0``.
    p, q : ndarray
        Pressures and observed amounts of the table.
    a, b : float
        The table's printed eq. (31) constants.
    p0_ref : float
        The recovered ``p0`` in bars.

    Returns
    -------
    float
        RMS of ``q - q_eq31`` over the largest ``|departure * q_eq31|``.
    """
    qf = q_eq31(p, a, b)
    sep = np.abs(departure(p / (p0_ref * f), a * p0_ref * f) * qf).max()
    return rms(q - qf) / sep
