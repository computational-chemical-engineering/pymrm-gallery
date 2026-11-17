"""Brunauer, Emmett & Teller (1938): the BET forms, the two ladders and the rival (page J1.3).

Brunauer, Emmett and Teller, J. Am. Chem. Soc. 60, 309 (1938). Everything here
is a closed form in one variable or a finite algebraic system in a layer index:

- the BET isotherm, eq. (26) = (A) (:func:`v_over_vm_A`), and its n-layer form,
  eq. (B) (:func:`v_eqB`);
- the inflection point of eq. (26), root-found on the numerator of its second
  derivative (:func:`d2_numerator`, :func:`x_inflection`,
  :func:`inflection_ratio`), and the c that maximises it, root-found on a
  central-difference derivative with the step as a parameter
  (:func:`d_ratio_dc`, :func:`c_star_root`);
- the layer equilibria (10) to (12) solved for the occupancies with ``newton``
  and eq. (15) summed term by term (:func:`ladder_v_over_vm`), and the
  continuous n that reproduces a target value of eq. (B) (:func:`n_for_target`);
- the rival's polarization recursion (1a) solved as a tridiagonal system
  (:func:`dipole_ladder`) against its closed forms (1b) and (1c)
  (:func:`C_eq1b`, :func:`C_eq1c`);
- the arithmetic of the tables: cross-sections from v_m and a surface
  (:func:`sigma_A2`), the rounding envelope of their ratio
  (:func:`rounding_envelope`), spreads and bands (:func:`spread`, :func:`band`,
  :func:`to10`, :func:`to5`) and the single-point ratio (:func:`shortcut_ratio`);
- a root finder that raises unless it finds exactly one root
  (:func:`unique_root`).
"""
import numpy as np
from pymrm import NumJac, newton
from pymrm.checks import find_roots

X_EDGE = 1e-13   # the inflection scan stays this far inside (0, 1), where eq. (26) is finite


# ------------------------------------------------------------- root finding
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


# ----------------------------------------------------------------- BET forms
def v_over_vm_A(x, c):
    """Eq. (26) = eq. (A) = eq. (B) at n -> infinity.

    Parameters
    ----------
    x : float or ndarray
        Relative pressure p/p0, in (0, 1).
    c : float
        The BET constant.

    Returns
    -------
    float or ndarray
        ``v/v_m = c x / ((1 - x)(1 - x + c x))``.
    """
    x = np.asarray(x, float)
    return c * x / ((1.0 - x) * (1.0 - x + c * x))


def v_eqB(x, c, n, vm=1.0):
    """Eq. (B): at most ``n`` layers (``n`` may be non-integer, as in Table IV).

    Parameters
    ----------
    x : float or ndarray
        Relative pressure.
    c : float
        The BET constant.
    n : float
        Largest number of layers.
    vm : float
        Monolayer volume; 1 gives ``v/v_m``.

    Returns
    -------
    float or ndarray
        The adsorbed volume of eq. (B).
    """
    x = np.asarray(x, float)
    return (vm * c * x / (1.0 - x) * (1.0 - (n + 1) * x**n + n * x**(n + 1))
            / (1.0 + (c - 1.0) * x - c * x**(n + 1)))


def d2_numerator(x, c):
    """Sign-carrying numerator of d2(v/v_m)/dx2 for eq. (26).

    Parameters
    ----------
    x : float
        Relative pressure.
    c : float
        The BET constant.

    Returns
    -------
    float
        ``-2 c ((c-1)^2 x^3 + 3 (c-1) x - c + 2)``.
    """
    return -2.0 * c * ((c - 1.0)**2 * x**3 + 3.0 * (c - 1.0) * x - c + 2.0)


def x_inflection(c):
    """The inflection point of eq. (26) in (0, 1), root-found.

    The cubic of :func:`d2_numerator` is increasing in x for c > 1 and equals
    ``2 - c`` at x = 0, so eq. (26) has exactly one inflection for c > 2 and
    none below (it is concave everywhere, with no knee to read a point B off).

    Parameters
    ----------
    c : float
        The BET constant.

    Returns
    -------
    float
        The inflection abscissa, or nan for c <= 2.
    """
    if c <= 2.0:
        return float("nan")
    return unique_root(lambda t: d2_numerator(t, c), X_EDGE, 1.0 - X_EDGE)


def inflection_ratio(c):
    """``v/v_m`` of eq. (26) at its inflection point (nan for c <= 2).

    Parameters
    ----------
    c : float
        The BET constant.

    Returns
    -------
    float
    """
    xi = x_inflection(c)
    return float("nan") if not np.isfinite(xi) else float(v_over_vm_A(xi, c))


def d_ratio_dc(cv, h, central=True):
    """Finite-difference derivative of :func:`inflection_ratio` with respect to c.

    Parameters
    ----------
    cv : float
        The BET constant.
    h : float
        Relative step: the derivative is taken from ``c (1 +- h)``.
    central : bool
        Central difference (error of order h^2); False gives the forward
        difference (order h).

    Returns
    -------
    float
    """
    if central:
        return (inflection_ratio(cv * (1 + h)) - inflection_ratio(cv * (1 - h))) / (2 * h * cv)
    return (inflection_ratio(cv * (1 + h)) - inflection_ratio(cv)) / (h * cv)


def c_star_root(h, central=True, lo=20.0, hi=200.0, xtol=1e-12):
    """The c that maximises the inflection value: the one root of :func:`d_ratio_dc`.

    Parameters
    ----------
    h : float
        Relative finite-difference step of the derivative.
    central : bool
        Central (True) or forward difference.
    lo, hi : float
        The range of c scanned (log grid).
    xtol : float
        Brent tolerance on c.

    Returns
    -------
    float
    """
    return unique_root(lambda cv: d_ratio_dc(cv, h, central), lo, hi, log=True, xtol=xtol)


# ---------------------------------------------- the layer ladder, solved
def ladder_v_over_vm(x, c, n, tol=1e-12):
    """Solve eqs. (10) to (12) for the occupancies and form eq. (15)'s sums term by term.

    Shares no algebra with :func:`v_eqB`: no geometric series is used. The
    residual is ``r_1 = s_1 - c x s_0`` and ``r_i = s_i - x s_(i-1)`` with
    ``s_0 = 1``, on a field of shape ``(n, 1)``; the residual reads the
    neighbouring layer, so the stencil is ``axes_diagonals=[0]``.

    Parameters
    ----------
    x : float
        Relative pressure.
    c : float
        The BET constant.
    n : int
        Number of layers.
    tol : float
        ``newton`` step tolerance on the occupancies scaled by ``c x``.

    Returns
    -------
    float
        ``v/v_m``.
    OptimizeResult
        The ``newton`` result (``nit`` is reported).

    Raises
    ------
    RuntimeError
        When ``newton`` does not report success.
    """
    jac = NumJac((n, 1), axes_diagonals=[0]) if n > 1 else NumJac((n, 1))

    def residual(s):
        r = np.empty_like(s)
        r[0, 0] = s[0, 0] - c * x * 1.0          # s_1 = c x s_0, s_0 = 1  (eqs. 10, 18, 22)
        if n > 1:
            r[1:, 0] = s[1:, 0] - x * s[:-1, 0]  # s_i = x s_(i-1)         (eqs. 11, 12, 19)
        return r

    sol = newton(lambda s: jac(residual, s), np.full((n, 1), c * x), tol=tol,
                 scale=c * x, stagnation=True, maxfev=50)
    if not sol.success:
        raise RuntimeError(f"ladder newton failed at x={x}, c={c}, n={n}: {sol.message}")
    s = sol.x.ravel()
    i = np.arange(1, n + 1)
    return float((i * s).sum() / (1.0 + s.sum())), sol


def n_for_target(xx, cc, target, n0, tol=1e-13):
    """The continuous n at which eq. (B) (with v_m = 1) equals ``target``, by ``newton``.

    Parameters
    ----------
    xx : float
        Relative pressure.
    cc : float
        The BET constant.
    target : float
        Required ``v/v_m``.
    n0 : float
        Starting n.
    tol : float
        ``newton`` step tolerance on n.

    Returns
    -------
    float

    Raises
    ------
    RuntimeError
        When ``newton`` does not report success.
    """
    jac = NumJac((1, 1))
    res = lambda nn: np.atleast_2d(float(v_eqB(xx, cc, float(nn.ravel()[0]), 1.0)) - target)
    sol = newton(lambda nn: jac(res, nn), np.array([[float(n0)]]), tol=tol,
                 stagnation=True, maxfev=80)
    if not sol.success:
        raise RuntimeError(f"n_for_target newton failed at x={xx}, c={cc}: {sol.message}")
    return float(sol.x.ravel()[0])


# ------------------------------------------ the polarization ladder, solved
def dipole_ladder(k, N):
    """Eq. (1a), ``mu_i = k (mu_(i-1) + mu_(i+1))`` for i = 1..N, with mu_0 = 1, mu_(N+1) = 0.

    Parameters
    ----------
    k : float
        The polarization coupling.
    N : int
        Number of layers.

    Returns
    -------
    ndarray
        ``mu_1 .. mu_N`` relative to ``mu_0``.
    """
    A = np.zeros((N, N))
    b = np.zeros(N)
    for i in range(N):
        A[i, i] = 1.0
        if i > 0:
            A[i, i - 1] = -k
        else:
            b[i] += k
        if i < N - 1:
            A[i, i + 1] = -k
    return np.linalg.solve(A, b)


def C_eq1b(k):
    """Eq. (1b), ``C = (1 - sqrt(1 - 4 k^2)) / (2 k)``: real only for k <= 1/2 (nan above).

    Parameters
    ----------
    k : float

    Returns
    -------
    float
    """
    disc = 1.0 - 4.0 * k * k
    return (1.0 - np.sqrt(disc)) / (2.0 * k) if disc >= 0.0 else np.nan


def C_eq1c(k):
    """Eq. (1c), de Boer and Zwicker's approximation ``C = k / (1 - k^2)``, the one Bradley used.

    Parameters
    ----------
    k : float

    Returns
    -------
    float
    """
    return k / (1.0 - k * k)


# ------------------------------------------------- arithmetic of the tables
def sigma_A2(vm_ccg, S_m2g, v_stp, n_av):
    """Cross-section per molecule implied by a v_m and a specific surface.

    Parameters
    ----------
    vm_ccg : float or array_like
        Monolayer volume, cm3 STP per g.
    S_m2g : float or array_like
        Specific surface, m2/g.
    v_stp : float
        Molar volume at STP, cm3/mol.
    n_av : float
        Avogadro's number, 1/mol.

    Returns
    -------
    float or ndarray
        Cross-section in Angstrom^2.
    """
    return S_m2g * 1e4 / (np.asarray(vm_ccg, float) / v_stp * n_av) * 1e16


def rounding_envelope(vm_a, S_a, vm_b, S_b, du_S=0.5, du_vm=0.05):
    """Largest relative difference, in %, that the printed rounding allows between two sigma.

    ``v_stp`` and ``n_av`` cancel in the ratio, so the envelope follows from the
    printed quantisation alone: surfaces to ``du_S`` and v_m to ``du_vm``.

    Parameters
    ----------
    vm_a, S_a, vm_b, S_b : float or array_like
        The two (v_m, surface) pairs.
    du_S, du_vm : float
        Half a unit in the last printed place of the surface and of v_m.

    Returns
    -------
    float or ndarray
    """
    return (du_S / S_a + du_vm / vm_a + du_S / S_b + du_vm / vm_b) * 100


def spread(v):
    """Mean and largest deviation from it in %.

    Parameters
    ----------
    v : array_like

    Returns
    -------
    (float, float)
    """
    v = np.asarray(v, float)
    m = float(v.mean())
    return m, float(np.abs(v / m - 1).max() * 100)


def band(v):
    """Midrange and half-range of ``v``: the paper's "X +- Y".

    Parameters
    ----------
    v : array_like

    Returns
    -------
    (float, float)
    """
    v = np.asarray(v, float)
    return (v.min() + v.max()) / 2.0, (v.max() - v.min()) / 2.0


def to10(v):
    """``v`` rounded to the nearest 10 (half to even)."""
    return float(np.round(float(v) / 10.0) * 10.0)


def to5(v):
    """``v`` rounded to the nearest 5 (half to even)."""
    return float(np.round(float(v) / 5.0) * 5.0)


def shortcut_ratio(cv, xv):
    """``v (1 - x) / v_m`` on eq. (26): the single-point rule's estimate of v_m over v_m.

    Parameters
    ----------
    cv : float or ndarray
        The BET constant.
    xv : float
        Relative pressure of the single point.

    Returns
    -------
    float or ndarray
        ``c x / (1 + (c - 1) x)``, below 1 for every c and x in (0, 1).
    """
    return cv * xv / (1.0 + (cv - 1.0) * xv)
