"""The Weisz-Prater criterion: the pellet of Weisz and Hicks (1962) and its references.

Model code of gallery page B1.4. Two routes to the same steady pellet,

``Pellet``
    pymrm finite volume on ``u = r/R``, the Thiele modulus ``phi_0`` given,
    returning the profile, the effectiveness factor and the observable
    ``Phi = phi_0^2 eta``;
``shoot`` and ``Branch``
    an initial-value integration outward from a prescribed centre
    concentration, which reaches every state of a solution branch (upper,
    middle and lower) in one deterministic pass;

plus the closed-form isothermal effectiveness factor and the apparent
activation energy of Weisz and Prater's printed relation.

Contents
--------
Power, InhibitedLH, WeiszHicks
    Rate laws ``R(y)`` normalised to ``R(1) = 1``, with their exact derivative.
Pellet
    ``nabla^2 y = phi_0^2 R(y)`` by pymrm, exact Jacobian.
shoot, Branch, fold_depth, turning_points, max_log_slope, fold_onset
    The shooting reference and what is read off a traced branch.
phi0_eq10, Phi_c_eq10a
    The thresholds of Weisz and Hicks' eqs. (10) and (10a).
threshold_exact, eta_exact, dln_eta_dln_phi, q_ratio
    Closed forms of the isothermal first-order pellet.

Several functions take an argument whose default is the correct model and whose
alternative is a defect the page's break table injects.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import i0e, i1e

from pymrm import (NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_div, construct_grad, newton)
from pymrm.checks import find_roots

__all__ = ["NU", "Power", "InhibitedLH", "WeiszHicks", "SolveFailed", "Pellet", "shoot",
           "Branch", "fold_depth", "turning_points", "max_log_slope", "fold_onset",
           "phi0_eq10", "Phi_c_eq10a", "threshold_exact", "eta_exact", "dln_eta_dln_phi",
           "q_ratio"]

NU = {"slab": 0, "cylinder": 1, "sphere": 2}


# ------------------------------------------------------------------ rate laws
class Power:
    """Power law ``R = y^n``.

    Parameters
    ----------
    n : float
        Reaction order, ``n >= 0``.
    """

    def __init__(self, n):
        self.n = float(n)
        self.label = {0.0: "zero order", 0.5: "half order", 1.0: "first order",
                      2.0: "second order"}.get(self.n, f"order {n:g}")

    def r(self, y):
        """Normalised rate."""
        y = np.asarray(y, float)
        return np.ones_like(y) if self.n == 0.0 else y ** self.n

    def dr(self, y):
        """Derivative ``dR/dy``."""
        y = np.asarray(y, float)
        return np.zeros_like(y) if self.n == 0.0 else self.n * y ** (self.n - 1.0)


class InhibitedLH:
    """Self-inhibited Langmuir-Hinshelwood rate ``R = (1 + K)^2 y/(1 + K y)^2``.

    The rate falls with concentration above ``y = 1/K``.

    Parameters
    ----------
    K : float, optional
        Adsorption constant times surface concentration.
    """

    def __init__(self, K=5.0):
        self.K = float(K)
        self.label = f"inhibited LH (K={K:g})"

    def r(self, y):
        """Normalised rate."""
        y = np.asarray(y, float)
        return (1.0 + self.K) ** 2 * y / (1.0 + self.K * y) ** 2

    def dr(self, y):
        """Derivative ``dR/dy``."""
        y = np.asarray(y, float)
        return (1.0 + self.K) ** 2 * (1.0 - self.K * y) / (1.0 + self.K * y) ** 3


class WeiszHicks:
    """First order with the Prater-reduced Arrhenius factor, Weisz and Hicks eq. (9).

    ``R = scale * y * exp(sign * gamma beta (1 - y)/(1 + beta (1 - y)))``.

    Parameters
    ----------
    beta : float
        Prater number, positive for an exothermic reaction.
    gamma : float
        Arrhenius number ``Q/(R_g T_0)``.
    sign : float, optional
        Sign of the exponent; ``-1`` is a defect for the break table.
    scale : float, optional
        Factor on the rate; anything but 1 is a defect for the break table.
    denominator : bool, optional
        ``False`` drops the ``1 + beta (1 - y)`` of the exponent, the
        exponential approximation (a defect for the break table).
    """

    def __init__(self, beta, gamma, sign=1.0, scale=1.0, denominator=True):
        self.beta, self.gamma = float(beta), float(gamma)
        self.sign, self.scale = float(sign), float(scale)
        self.denominator = denominator
        self.label = f"beta = {beta:g}, gamma = {gamma:g}"

    def _expo(self, y):
        w = 1.0 + self.beta * (1.0 - y) if self.denominator else np.ones_like(y)
        return self.sign * self.gamma * self.beta * (1.0 - y) / w, w

    def r(self, y):
        """Normalised rate."""
        y = np.asarray(y, float)
        e, _ = self._expo(y)
        return self.scale * y * np.exp(e)

    def dr(self, y):
        """Derivative ``dR/dy``."""
        y = np.asarray(y, float)
        e, w = self._expo(y)
        de = -self.sign * self.gamma * self.beta / w ** 2
        return self.scale * np.exp(e) * (1.0 + y * de)


# ------------------------------------------------------------------ pymrm route
class SolveFailed(RuntimeError):
    """Raised when ``newton`` reports no success."""


class Pellet:
    """``nabla^2 y = phi_0^2 R(y)`` on ``0 < u < 1`` by pymrm finite volume.

    Symmetry at the centre, ``y = 1`` at the surface. Layout ``(n_u, 1)``:
    the spatial axis first, the single field last.

    Parameters
    ----------
    geom : {"slab", "cylinder", "sphere"}
        Geometry; sets ``nu``.
    n_u : int
        Number of cells (uniform).
    nu : int, optional
        Geometry index passed to ``construct_div`` instead of ``NU[geom]``
        (a wrong value is a defect for the break table).
    """

    def __init__(self, geom="sphere", n_u=400, nu=None):
        self.geom = geom
        self.nu = NU[geom] if nu is None else nu
        self.n_u, self.shape = n_u, (n_u, 1)
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])
        # outward normal: at u = 0 it points to -u, at u = 1 to +u
        #   centre  u = 0: symmetry, dy/dn = 0       -> a = 1, b = 0, d = 0
        #   surface u = 1: y = 1 (external value)    -> a = 0, b = 1, d = 1
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
        grad, grad_bc = construct_grad(self.shape, self.u_f, self.u_c, self.bc)
        div = construct_div(self.shape, self.u_f, nu=self.nu)   # 0 slab, 1 cylinder, 2 sphere
        self.lap, self.lap_bc = div @ grad, div @ grad_bc       # assembled once
        self.numjac = NumJac(self.shape)                       # only for jacobian="numjac"
        vf = self.u_f ** (self.nu + 1) / (self.nu + 1.0)
        self.dv, self.v_tot = np.diff(vf), vf[-1]

    def residual(self, y, phi, rate, jacobian="exact"):
        """Residual ``lap y - phi^2 R(y)`` and its Jacobian.

        Parameters
        ----------
        y : numpy.ndarray
            State, any shape with ``n_u`` entries.
        phi : float
            Thiele modulus.
        rate : Power, InhibitedLH or WeiszHicks
            Rate law.
        jacobian : {"exact", "numjac"}, optional
            ``"exact"`` puts ``phi^2 dR/dy`` on the diagonal with
            ``construct_coefficient_matrix``; ``"numjac"`` differences the
            source with ``NumJac``.

        Returns
        -------
        g : numpy.ndarray
            Residual, shape ``(n_u, 1)``.
        jac : scipy.sparse matrix
            Jacobian.
        """
        y = np.asarray(y, float).reshape(self.shape)
        p2 = phi * phi
        if jacobian == "exact":
            g_src = p2 * rate.r(y)
            jac_src = construct_coefficient_matrix(p2 * rate.dr(y))
        else:
            g_src, jac_src = self.numjac(lambda v: p2 * rate.r(v), y)
        g = self.lap @ y.reshape((-1, 1)) + self.lap_bc - np.asarray(g_src).reshape((-1, 1))
        return g, self.lap - jac_src

    def solve(self, phi, rate, y_init=None, maxfev=60, tol=1e-12, check=True,
              jacobian="exact", lin_solver_kwargs=None):
        """Steady profile at Thiele modulus ``phi``.

        ``newton`` with ``bounds=(0, 1)``: the maximum principle keeps ``y`` in
        ``[0, 1]`` for every rate law here, and above ``y = 1 + 1/beta`` the
        Arrhenius factor of eq. (9) is undefined. ``stagnation=True`` stops on a
        step that no longer decreases.

        Parameters
        ----------
        phi : float
            Thiele modulus.
        rate : Power, InhibitedLH or WeiszHicks
            Rate law.
        y_init : numpy.ndarray, optional
            Initial guess; default ``y = 1``.
        maxfev : int, optional
            Iteration limit.
        tol : float, optional
            Step tolerance of ``newton``.
        check : bool, optional
            Raise :class:`SolveFailed` when ``newton`` reports no success.
        jacobian : {"exact", "numjac"}, optional
            See :meth:`residual`.
        lin_solver_kwargs : dict, optional
            Passed to ``newton`` (``splu`` options such as ``permc_spec``).

        Returns
        -------
        y : numpy.ndarray
            Cell values, shape ``(n_u,)``.
        sol : OptimizeResult
            The ``newton`` result.
        """
        y0 = np.ones(self.shape) if y_init is None else np.asarray(y_init, float)
        sol = newton(lambda y: self.residual(y, phi, rate, jacobian), y0.reshape((-1, 1)),
                     tol=tol, maxfev=maxfev, bounds=(0.0, 1.0), stagnation=True,
                     lin_solver_kwargs=lin_solver_kwargs)
        if check and not sol.success:
            raise SolveFailed(f"phi_0 = {phi:g}, {rate.label}: {sol.message}")
        return np.asarray(sol.x).ravel(), sol

    def Phi(self, y, two_point=False):
        """Observable modulus ``(nu + 1) dy/du`` at the surface; needs no ``phi_0``.

        The surface gradient is the one the discrete operator transports, from
        ``compute_boundary_values`` (along ``+u``, the outward normal there).

        Parameters
        ----------
        y : numpy.ndarray
            Cell values.
        two_point : bool, optional
            ``True`` takes the gradient from the last cell centre to the
            surface, a first-order read (a defect for the break table).
        """
        y = np.asarray(y, float).reshape(self.shape)
        if two_point:
            g_s = (1.0 - y[-1, 0]) / (self.u_f[-1] - self.u_c[-1])
        else:
            _, g_s = compute_boundary_values(y, self.u_f, self.u_c, self.bc[1], bound_id=1)
        return (self.nu + 1.0) * float(np.ravel(g_s)[0])

    def eta(self, y, phi, two_point=False):
        """Effectiveness factor ``Phi/phi_0^2``."""
        return self.Phi(y, two_point) / phi ** 2

    def threshold(self, rate, target=1.0, lo=0.05, hi=3.0, n_scan=16):
        """Every Thiele modulus at which ``Phi = target``, from cold-started solves.

        A :func:`pymrm.checks.find_roots` scan in ``ln phi_0``; a solve that
        fails is reported by the scan, never taken as a sign.

        Returns
        -------
        dict
            The ``find_roots`` result, roots converted to ``phi_0``.
        """
        def f(lphi):
            return self.Phi(self.solve(float(np.exp(lphi)), rate)[0]) - target

        res = find_roots(f, np.log(lo), np.log(hi), n_scan=n_scan, xtol=1e-13)
        res["roots"] = [float(np.exp(r)) for r in res["roots"]]
        return res

    def Phi_from_volume(self, y, phi, rate):
        """``phi_0^2`` times the volume average of the rate.

        Summing the discrete residual over the cells, weighted by their
        volumes, telescopes the divergence to the surface flux, so
        ``Phi - Phi_from_volume`` is the volume-weighted residual: the same sum,
        not a second route.
        """
        return phi * phi * float(np.sum(rate.r(np.asarray(y).ravel()) * self.dv) / self.v_tot)


# ------------------------------------------------------------- shooting route
def shoot(nu, rate, t, rtol=1e-10, dense=False):
    """One branch state from its centre concentration ``y_c = 1 - exp(-t)``.

    Integrates ``y'' + (nu/s) y' = R(y)`` in ``s = phi_0 u`` outward from the
    centre until ``y = 1``; that radius is ``phi_0``. No grid, no operator, no
    Newton iteration: it shares only the rate law with :class:`Pellet`.

    Parameters
    ----------
    nu : int
        Geometry index.
    rate : Power, InhibitedLH or WeiszHicks
        Rate law.
    t : float
        Branch parameter, ``t > 0``.
    rtol : float, optional
        Relative tolerance of the DOP853 integration.
    dense : bool, optional
        Also return the ``solve_ivp`` result with dense output.

    Returns
    -------
    tuple or None
        ``(phi_0, eta, Phi)`` (plus the integration result when ``dense``), or
        ``None`` when ``y`` never reaches 1.
    """
    y_c = -np.expm1(-t)
    r0 = float(rate.r(y_c))
    if not (0.0 < y_c < 1.0) or r0 <= 0.0:
        return None
    # start on the series y = y_c + R s^2/(2 (nu + 1)), offset 1e-10 of the
    # distance from y_c to the nearer of 0 and 1
    s0 = min(1e-3, np.sqrt(2 * (nu + 1) * 1e-10 * min(y_c, np.exp(-t)) / r0))

    def rhs(s, v):
        return [v[1], float(rate.r(min(max(v[0], 0.0), 1.0))) - (nu / s) * v[1]]

    def hit(s, v):
        return v[0] - 1.0

    hit.terminal, hit.direction = True, 1
    v0 = [y_c + r0 * s0 ** 2 / (2 * (nu + 1)), r0 * s0 / (nu + 1)]
    sol = solve_ivp(rhs, (s0, 1e5), v0, events=hit, rtol=rtol, atol=1e-16,
                    method="DOP853", dense_output=dense)
    if not sol.t_events[0].size:
        return None
    s = float(sol.t_events[0][0])
    yp = float(sol.y_events[0][0][1])
    out = (s, (nu + 1.0) * yp / s, (nu + 1.0) * s * yp)
    return out + (sol,) if dense else out


class Branch:
    """One solution branch, traced once on a fixed ladder of centre concentrations.

    The branch parameter is ``tau = ln t`` with ``y_c = 1 - exp(-t)``, on
    ``n`` equally spaced values; the arrays run from the mild end
    (``y_c -> 1``) to the severe end. Every state is cached, so a root scan with
    :func:`pymrm.checks.find_roots` on the same ladder costs only its
    refinements.

    Parameters
    ----------
    nu : int
        Geometry index.
    rate : Power, InhibitedLH or WeiszHicks
        Rate law.
    n : int, optional
        Number of traced states.
    tmin, tmax : float, optional
        Range of ``t``.
    rtol : float, optional
        Tolerance of every integration.
    """

    def __init__(self, nu, rate, n=110, tmin=1e-9, tmax=24.0, rtol=1e-9):
        self.nu, self.rate, self.rtol, self.n = nu, rate, rtol, n
        self.lo, self.hi = float(np.log(tmin)), float(np.log(tmax))
        self._cache = {}
        taus = np.linspace(self.lo, self.hi, n)[::-1]
        rows = [(tau, self.state(tau)) for tau in taus]
        rows = [(tau, o) for tau, o in rows if o is not None and np.all(np.isfinite(o))]
        self.tau = np.array([r[0] for r in rows])
        self.phi = np.array([r[1][0] for r in rows])
        self.eta = np.array([r[1][1] for r in rows])
        self.Phi = np.array([r[1][2] for r in rows])

    def state(self, tau):
        """``(phi_0, eta, Phi)`` at ``tau = ln t``, or ``None``; cached."""
        key = float(tau)
        if key not in self._cache:
            self._cache[key] = shoot(self.nu, self.rate, float(np.exp(key)), self.rtol)
        return self._cache[key]

    def all_where(self, f):
        """Every state along the branch where ``f(phi_0, eta, Phi)`` changes sign.

        A :func:`pymrm.checks.find_roots` scan on the traced ladder, each root
        refined by Brent's method in ``tau`` to ``1e-12``. All roots are
        returned, mild end first.

        Returns
        -------
        list of tuple
            ``(tau, phi_0, eta, Phi)`` per root.
        """
        def g(tau):
            o = self.state(tau)
            return np.nan if o is None else f(*o)

        res = find_roots(g, self.lo, self.hi, n_scan=self.n, xtol=1e-12)
        return [(tau,) + self.state(tau) for tau in res["roots"][::-1]]

    def where(self, f):
        """The first such state (mild end first) as ``(phi_0, eta, Phi)``, or ``None``."""
        r = self.all_where(f)
        return r[0][1:] if r else None

    def profile(self, tau, u):
        """``(phi_0, y(u))`` of the state at ``tau``, ``u`` scaled on the radius."""
        o = shoot(self.nu, self.rate, float(np.exp(tau)), self.rtol, dense=True)
        phi, sol = o[0], o[3]
        return phi, np.clip(sol.sol(np.clip(u * phi, sol.t[0], phi))[0], 0.0, 1.0)


def fold_depth(x):
    """Largest relative backward excursion of a sequence traced in branch order.

    Zero means monotone: the quantity determines the state uniquely.
    """
    run = np.maximum.accumulate(x)
    return float(((run - x) / run).max())


def turning_points(br, which="Phi", refine=True):
    """Every local extremum of ``phi_0`` or ``Phi`` along a traced branch.

    Each extremum bracketed by the traced ladder is refined by a bounded Brent
    search in ``tau`` on the shooting solution.

    Parameters
    ----------
    br : Branch
        Traced branch.
    which : {"Phi", "phi"}, optional
        Quantity whose extrema are located.
    refine : bool, optional
        ``False`` returns the traced ladder point (a defect for the break table).

    Returns
    -------
    list of tuple
        ``("max" or "min", tau, phi_0, eta, Phi)`` in branch order.
    """
    x = br.Phi if which == "Phi" else br.phi
    k = 2 if which == "Phi" else 0
    out = []
    for i in range(1, len(x) - 1):
        if (x[i] - x[i - 1]) * (x[i + 1] - x[i]) >= 0:
            continue
        kind = "max" if x[i] > x[i - 1] else "min"
        if not refine:
            out.append((kind, float(br.tau[i]), br.phi[i], br.eta[i], br.Phi[i]))
            continue
        lo, hi = sorted((float(br.tau[i - 1]), float(br.tau[i + 1])))
        sgn = -1.0 if kind == "max" else 1.0
        r = minimize_scalar(lambda tt: sgn * shoot(br.nu, br.rate, np.exp(tt), br.rtol)[k],
                            bounds=(lo, hi), method="bounded",
                            options={"xatol": (hi - lo) * 1e-8})
        out.append((kind, float(r.x)) + shoot(br.nu, br.rate, np.exp(r.x), br.rtol))
    return out


def max_log_slope(br, which="Phi", h=1e-3):
    """Largest ``d ln x/d tau`` along a branch, ``x`` = ``Phi`` or ``phi_0``.

    Along an unfolded branch ``x`` falls as ``tau`` grows (towards the mild
    end), so the slope is negative everywhere; a fold makes it positive
    somewhere. The largest traced difference quotient is refined by a bounded
    Brent search on central differences of the shooting solution (half-step
    ``h`` in ``tau``).
    """
    x = np.log(br.Phi if which == "Phi" else br.phi)
    k = 2 if which == "Phi" else 0
    tau = br.tau
    s = (x[:-1] - x[1:]) / (tau[:-1] - tau[1:])
    i = int(np.argmax(s))
    j0, j1 = max(i - 1, 0), min(i + 2, len(tau) - 1)
    lo, hi = float(min(tau[j0], tau[j1])), float(max(tau[j0], tau[j1]))

    def slope(tt):
        a = shoot(br.nu, br.rate, np.exp(tt + h), br.rtol)
        b = shoot(br.nu, br.rate, np.exp(tt - h), br.rtol)
        return (np.log(a[k]) - np.log(b[k])) / (2.0 * h)

    r = minimize_scalar(lambda tt: -slope(tt), bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-6})
    return float(-r.fun)


def fold_onset(nu, rate_of, which, lo, hi, n=60, rtol=1e-11, xtol=1e-6):
    """Parameter value at which a fold in ``Phi`` or ``phi_0`` first appears.

    The root of :func:`max_log_slope` in the parameter ``p`` of
    ``rate_of(p)``: the cusp where the branch stops being monotone. No
    fold-depth threshold is involved.

    Parameters
    ----------
    nu : int
        Geometry index.
    rate_of : callable
        ``rate_of(p)`` returns the rate law at parameter ``p``.
    which : {"Phi", "phi"}
        Quantity whose fold is located.
    lo, hi : float
        Bracket: no fold at ``lo``, a fold at ``hi``.
    n : int, optional
        States traced per branch before refinement.
    rtol : float, optional
        Integration tolerance.
    xtol : float, optional
        Absolute tolerance on ``p``.
    """
    def f(p):
        return max_log_slope(Branch(nu, rate_of(p), n=n, rtol=rtol), which)

    return float(brentq(f, lo, hi, xtol=xtol))


def phi0_eq10(beta, gamma, denominator=True):
    """Largest safe Thiele modulus of eq. (10), ``exp(-gamma beta/(2 (1 + beta)))``.

    ``denominator=False`` drops the ``1 + beta`` (a defect for the break table).
    """
    return float(np.exp(-0.5 * gamma * beta / ((1.0 + beta) if denominator else 1.0)))


def Phi_c_eq10a(beta, gamma, denominator=True):
    """Threshold on the observable of eq. (10a), ``exp(-gamma beta/(1 + beta))``.

    ``denominator=False`` drops the ``1 + beta`` (a defect for the break table).
    """
    return float(np.exp(-gamma * beta / ((1.0 + beta) if denominator else 1.0)))


def threshold_exact(nu, target=1.0, length="R"):
    """Isothermal first-order state at which the observable equals ``target``.

    Parameters
    ----------
    nu : int
        Geometry index.
    target : float, optional
        Threshold on the observable.
    length : {"R", "V/S", "V/S unsquared"}, optional
        Length scale of the observable. On ``L = V/S = R/(nu + 1)`` the
        observable is ``Phi/(nu + 1)^2``; ``"V/S unsquared"`` divides by
        ``nu + 1`` only (a defect for the break table).

    Returns
    -------
    phi, eta : float
    """
    scale = {"R": 1.0, "V/S": (nu + 1.0) ** 2, "V/S unsquared": nu + 1.0}[length]
    res = find_roots(lambda p: p * p * eta_exact(nu, p) - target * scale, 1e-3, 1e3,
                     n_scan=32, log=True, xtol=1e-14)
    phi = res["roots"][0]
    return phi, eta_exact(nu, phi)


# ----------------------------------------------------------------- closed forms
def eta_exact(nu, phi):
    """Isothermal first-order effectiveness factor, closed form.

    Slab ``tanh(phi)/phi``, cylinder ``2 I1(phi)/(phi I0(phi))`` (scaled Bessel
    functions, no overflow), sphere ``3 (phi coth phi - 1)/phi^2``. For
    ``phi < 0.1`` the sphere uses its Taylor series to ``phi^10``, where the
    closed form loses digits to cancellation.
    """
    phi = float(phi)
    if nu == 0:
        return np.tanh(phi) / phi
    if nu == 1:
        return 2.0 * i1e(phi) / (phi * i0e(phi))
    if phi < 0.1:
        x = phi * phi
        return 1.0 - x / 15 + 2 * x ** 2 / 315 - x ** 3 / 1575 + 2 * x ** 4 / 31185 \
            - 1382 * x ** 5 / 212837625
    return (3.0 / phi) * (1.0 / np.tanh(phi) - 1.0 / phi)


def dln_eta_dln_phi(nu, phi, h=1e-4, stable=True):
    """``d ln eta / d ln phi`` of :func:`eta_exact` by a central difference in ``ln phi``.

    Parameters
    ----------
    nu : int
        Geometry index.
    phi : float
        Thiele modulus.
    h : float, optional
        Half-step in ``ln phi``.
    stable : bool, optional
        ``False`` evaluates the sphere's closed form at every ``phi``, without
        the series (a defect for the break table).
    """
    def e(p):
        if nu == 2 and not stable:
            return (3.0 / p) * (1.0 / np.tanh(p) - 1.0 / p)
        return eta_exact(nu, p)

    return (np.log(e(phi * np.exp(h))) - np.log(e(phi * np.exp(-h)))) / (2.0 * h)


def q_ratio(nu, phi, factor=0.5, h=1e-4, stable=True):
    """Apparent over true activation energy, ``Q'/Q = 1 + factor dln eta/dln phi_0``.

    Weisz and Prater's relation as printed by Weisz and Hicks (journal page 268)
    has ``factor = 1/2``.
    """
    return 1.0 + factor * dln_eta_dln_phi(nu, phi, h, stable)
