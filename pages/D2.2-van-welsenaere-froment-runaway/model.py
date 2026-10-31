"""Runaway criteria of Van Welsenaere and Froment (1970) and the 1-D tube they apply to.

Model code of gallery page D2.2. The closed-form criteria and the phase-plane
quadrature need only the groups ``A``, ``B``, ``C`` and the Arrhenius constants
``a``, ``b``; the reactor is a pymrm finite-volume model of the same equations.

Contents
--------
VWFCase
    The printed constants of the paper's base case and the groups ``A``, ``B``,
    ``C`` of its Eq. 4.
kexp
    The Arrhenius factor ``exp(-a/T + b)``.
maxima_curve, ps_curve
    The two loci of the p-T plane (Eqs. 7 and 18).
t_max_peak
    ``T_M``, the maximum of the maxima curve (Eq. 8), criterion 1.
t_tangency
    ``(T_i)_t`` from the tangency condition (Eq. 20), criterion 2.
inlet_limits
    Lower and upper inlet pressures (Eqs. 27 and 28) and their mean.
back_integrate
    The true critical inlet pressure: the trajectory through a critical point
    followed back to the inlet in the p-T plane.
critical_points
    Both criteria's critical points and their back-integrated inlet pressures.
RunawayTube
    The 1-D pseudo-homogeneous tube: first-order upwind convection in the
    contact time, pointwise source, ``pymrm.newton``.
critical_p0
    The inlet pressure whose reactor hot spot reaches a target temperature,
    by a root-finder on independent reactor solves.
critical_p0_bisect
    The same, by bisection.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from pymrm import NumJac, construct_convflux_upwind, construct_div, continuation, newton

__all__ = ["VWFCase", "kexp", "maxima_curve", "ps_curve", "t_max_peak", "t_tangency",
           "inlet_limits", "back_integrate", "critical_points", "RunawayTube",
           "critical_p0", "critical_p0_bisect"]


@dataclass(frozen=True)
class VWFCase:
    """Printed constants of the base case (Section 1), units kmol, m, hr, kcal, K, atm.

    Parameters
    ----------
    M : float
        Mean molecular weight (kg/kmol).
    P : float
        Total pressure (atm).
    rho_b, rho_g : float
        Catalyst bulk density and gas density (kg/m3).
    minus_dH : float
        Heat of reaction, positive for exothermic (kcal/kmol).
    c_p : float
        Heat capacity of the gas, printed VOLUMETRIC (kcal/(m3 K)).
    u : float
        Superficial velocity (m/hr).
    U : float
        Overall wall heat-transfer coefficient (kcal/(m2 hr K)).
    R : float
        Tube radius (m).
    p_B0 : float
        Partial pressure of the second reactant, in excess (atm).
    a, b : float
        Arrhenius constants of ``exp(-a/T + b)``: ``a = E/R_gas`` (K), ``b`` the
        logarithm of the frequency factor.
    """

    M: float
    P: float
    rho_b: float
    rho_g: float
    minus_dH: float
    c_p: float
    u: float
    U: float
    R: float
    p_B0: float
    a: float
    b: float

    @classmethod
    def from_table(cls, table, **changes):
        """Build the case from ``{symbol: value}`` as in the page's parameter CSV."""
        kw = {k: float(v) for k, v in dict(table).items()}
        kw.update(changes)
        return cls(**kw)

    @property
    def A(self):
        """Consumption group ``M P rho_b p_B0 / rho_g`` (1/hr per unit rate factor)."""
        return self.M * self.P * self.rho_b / self.rho_g * self.p_B0

    @property
    def B(self):
        """Heating group ``(-dH) rho_b p_B0 / c_p`` (K/(hr atm) per unit rate factor)."""
        return self.minus_dH * self.rho_b / self.c_p * self.p_B0

    @property
    def C(self):
        """Cooling group ``2 U / (c_p R)`` (1/hr). ``c_p`` is volumetric: no ``rho_g``."""
        return 2.0 * self.U / (self.c_p * self.R)

    @property
    def contact_time(self):
        """Contact time of a 1 m tube, ``1/u`` (hr)."""
        return 1.0 / self.u


def kexp(case, T):
    """Arrhenius factor ``exp(-a/T + b)`` for ``T > 0`` (K)."""
    return np.exp(-case.a / T + case.b)


def maxima_curve(case, T, Tw):
    """Eq. 7: the partial pressure at which ``dT/dz = 0`` at temperature ``T``."""
    return (T - Tw) / (case.B / case.C * kexp(case, T))


def ps_curve(case, T):
    """Eq. 18: the locus of the inflexion points before the maximum."""
    return T ** 2 / (case.a * case.B / case.C * kexp(case, T))


def t_max_peak(case, Tw):
    """Eq. 8: ``T_M``, the temperature at the maximum of the maxima curve."""
    return 0.5 * (case.a - np.sqrt(case.a * (case.a - 4.0 * Tw)))


def t_tangency(case, Tw, n_scan=4000):
    """``(T_i)_t`` from Eq. 20, solved directly for ``t = a / (T_i)_t``.

    The paper reaches it through its Figs. 6 and 7. Here the function is
    scanned for sign changes between ``0.8 a / T_M`` and ``t_w = a / T_w`` and
    the last bracket is refined with ``brentq``.

    Returns
    -------
    float
        ``(T_i)_t`` (K), or ``nan`` when Eq. 20 has no root in the scan.
    """
    K = case.A / case.C * np.exp(case.b - 20.0)
    tw = case.a / Tw

    def f(t):
        return (t - 2.0) / (t * np.exp(-t + 20.0)) * (1.0 - t * (1.0 - t / tw)) - K

    ts = np.linspace(case.a / t_max_peak(case, Tw) * 0.80, tw - 1e-9, n_scan)
    sign = np.nonzero(np.diff(np.sign(f(ts))) != 0)[0]
    if len(sign) == 0:
        return np.nan
    i = sign[-1]
    return case.a / brentq(f, ts[i], ts[i + 1])


def inlet_limits(case, T_cr, Tw, T0=None):
    """Eqs. 27 and 28 and their mean, the inlet pressures that bracket a criterion.

    Eq. 27 is the general lower limit; at the first criterion's ``T_M`` it
    coincides with Eq. 29, because Eq. 8 makes ``T_M^2 = a (T_M - T_w)``.

    Returns
    -------
    tuple of float
        ``(lower, upper, mean)`` (atm).
    """
    T0 = Tw if T0 is None else T0
    X = case.A / case.C * kexp(case, T_cr)
    lower = T_cr ** 2 / (case.a * case.B / case.C * kexp(case, T_cr)) \
        + case.A / case.B * (T_cr - T0)
    upper = case.A / case.B * (T_cr - Tw) * (1.0 / np.sqrt(X) + 1.0) ** 2
    return lower, upper, 0.5 * (lower + upper)


def back_integrate(case, T_cr, p_cr, Tw, p_top=0.3):
    """Inlet pressure of the trajectory through ``(T_cr, p_cr)``: the paper's exact value.

    Dividing the energy balance by the mass balance removes ``z`` and leaves
    ``dT/dp`` (their Eq. 5). The trajectory is followed from the critical point
    upstream, where ``p`` rises, until it reaches the inlet temperature
    ``T = T_w``. A phase-plane quadrature (adaptive Runge-Kutta, ``rtol``
    1e-11): no reactor grid, so it is independent of :class:`RunawayTube`.

    Returns
    -------
    float
        Inlet partial pressure (atm), or ``nan`` when ``T_w`` is not reached
        below ``p_top``.
    """
    A, B, C = case.A, case.B, case.C

    def dT_dp(p, T):
        return [-B / A + (C / A) * (T[0] - Tw) / (p * kexp(case, T[0]))]

    def reached_inlet(p, T):
        return T[0] - Tw
    reached_inlet.terminal, reached_inlet.direction = True, -1
    sol = solve_ivp(dT_dp, (p_cr, p_top), [T_cr], events=reached_inlet,
                    rtol=1e-11, atol=1e-12, max_step=p_cr / 500)
    return float(sol.t_events[0][0]) if len(sol.t_events[0]) else np.nan


def critical_points(case, Tw):
    """Critical points of both criteria and their back-integrated inlet pressures.

    Returns
    -------
    dict
        ``T1, p1, p0_1`` (criterion 1: the maximum of the maxima curve) and
        ``T2, p2, p0_2`` (criterion 2: tangency with the ``p_s`` curve).
    """
    TM = t_max_peak(case, Tw)
    p1 = maxima_curve(case, TM, Tw)
    Tt = t_tangency(case, Tw)
    p2 = ps_curve(case, Tt)
    return dict(T1=TM, p1=p1, p0_1=back_integrate(case, TM, p1, Tw),
                T2=Tt, p2=p2, p0_2=back_integrate(case, Tt, p2, Tw))


class RunawayTube:
    """1-D pseudo-homogeneous tube with a constant wall temperature.

    ``dp/dz = -A p k(T)``, ``dT/dz = B p k(T) - C (T - T_w)`` with
    ``k(T) = exp(-a/T + b)`` and ``z`` the contact time ``z'/u`` (hr), the
    paper's own coordinate. Written as a steady convection problem with unit
    velocity: first-order upwind fluxes, ``construct_div`` on a uniform grid,
    the source by ``NumJac`` per cell, solved by ``pymrm.newton``.

    State ``(n_z, 2)``: ``[..., 0] = p`` (atm), ``[..., 1] = T`` (K).

    Parameters
    ----------
    case : VWFCase
        Constants and groups.
    p0 : float
        Inlet partial pressure (atm).
    Tw : float
        Wall temperature (K).
    T0 : float, optional
        Inlet temperature (K); ``Tw`` by default, as in the paper.
    length : float, optional
        Tube length (m).
    n_z : int, optional
        Number of cells.

    Attributes
    ----------
    jac_const : scipy.sparse.csc_array
        Constant convection operator ``div @ conv``.
    g_inlet : scipy.sparse.csc_array
        Inlet contribution per unit inlet value, shape ``(2 n_z, 2)``:
        ``g_const = g_inlet @ (p_in, T_in)``.
    """

    #: ``newton`` bounds: ``T > 0``, where the Arrhenius factor is defined.
    BOUNDS = (np.array([-np.inf, 0.0]), None)

    def __init__(self, case, p0, Tw, T0=None, length=1.0, n_z=1500):
        self.case, self.Tw = case, float(Tw)
        self.shape = (int(n_z), 2)
        self._build_grid(length * case.contact_time)
        self._build_operators()
        self.set_inlet(p0, self.Tw if T0 is None else T0)
        self.u = np.tile(self.u_in, (self.shape[0], 1))

    def _build_grid(self, z_end):
        self.z_f = np.linspace(0.0, z_end, self.shape[0] + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])

    def _build_operators(self):
        # z = 0, outward normal -z: p = p_in, T = T_in       -> a = 0, b = 1, d = 1 x (p_in, T_in)
        #        (shapes_d: d is the coefficient on the inlet vector passed in set_inlet)
        # z = L, outward normal +z: dp/dz = dT/dz = 0 (outflow) -> a = 1, b = 0, d = 0
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        conv_mat, conv_bc_in, _ = construct_convflux_upwind(
            self.shape, self.z_f, self.z_c, bc, v=1.0, axis=0, shapes_d=((1, 2), None))
        # v = 1: z is the contact time, so d(v x)/dz = source. nu = 0: Cartesian.
        div_mat = construct_div(self.shape, self.z_f, nu=0, axis=0)
        self.jac_const = (div_mat @ conv_mat).tocsc()
        self.g_inlet = (div_mat @ conv_bc_in).tocsc()
        self.numjac = NumJac(self.shape)          # pointwise source: last axis in full

    def set_inlet(self, p0, T0=None):
        """Set the inlet values; the operators are not rebuilt."""
        T0 = self.u_in[1] if T0 is None else T0
        self.p0 = float(p0)
        self.u_in = np.array([self.p0, float(T0)])
        self.g_const = (self.g_inlet @ self.u_in).reshape((-1, 1))

    def reaction(self, u):
        """Source ``(-A r, B r - C (T - T_w))`` with ``r = p k(T)``."""
        p, T = u[..., 0], u[..., 1]
        rate = p * kexp(self.case, T)
        return np.stack([-self.case.A * rate,
                         self.case.B * rate - self.case.C * (T - self.Tw)], axis=-1)

    def residual(self, u):
        """Steady residual and Jacobian: convection minus source."""
        g_rxn, jac_rxn = self.numjac(self.reaction, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self, u0=None, tol=1e-10, maxfev=50):
        """Newton from ``u0`` (default: the inlet state everywhere).

        The step test is on the unknowns scaled by the inlet values (``scale``),
        so ``p`` and ``T`` are judged relative to their own size. ``T`` is
        bounded below by 0, where the Arrhenius factor is undefined.

        Raises
        ------
        RuntimeError
            When Newton does not converge (``NewtonBreakdown`` propagates).
        """
        u0 = np.tile(self.u_in, (self.shape[0], 1)) if u0 is None else u0
        result = newton(self.residual, np.asarray(u0, dtype=float).reshape(self.shape),
                        tol=tol, maxfev=maxfev, scale=np.abs(self.u_in), bounds=self.BOUNDS)
        if not result.success:
            raise RuntimeError(f"RunawayTube.solve at p0 = {self.p0:g}: {result.message}")
        self.u, self.result = result.x.reshape(self.shape), result
        return self

    def trace_inlet(self, p0s, max_change=None, tol=1e-10, maxfev=50):
        """Follow the solution over the inlet pressures ``p0s`` with ``pymrm.continuation``.

        Each point starts from a secant prediction of the two before it. The
        branch-jump check of ``continuation`` is off by default. Past the
        critical inlet pressure the hot spot moves upstream, and a secant
        prediction can be hundreds of kelvin off in the cells it passes, so the
        check rejects corrections that are the solution. Pass ``max_change``
        to switch it on.

        Uniqueness is a property of the differential equations (an
        initial-value problem in ``z``), not of the discrete system: on a
        coarse grid some backward-Euler cells near the inlet have further,
        ignited roots in the runaway region. Newton from the neighbouring
        solution follows the low root, the one that converges to the
        differential solution under refinement; check the traced states
        against independent solves (as the page does) before relying on them.

        On return the tube holds the last traced state and its inlet value.

        Returns
        -------
        scipy.optimize.OptimizeResult
            The result of ``pymrm.continuation``; ``xs`` holds the states.
        """
        p0s = np.asarray(p0s, dtype=float)

        def residual(u, p0):
            self.set_inlet(p0)
            return self.residual(u)

        self.set_inlet(p0s[0])
        start = self.solve(tol=tol, maxfev=maxfev).u
        scale = np.array([np.max(np.abs(p0s)), abs(self.u_in[1])])
        res = continuation(residual, start, p0s, max_change=max_change,
                           newton_kwargs=dict(tol=tol, maxfev=maxfev, scale=scale,
                                              bounds=self.BOUNDS))
        self.set_inlet(res.last_param if res.last_param is not None else p0s[0])
        self.u = np.asarray(res.x, dtype=float).reshape(self.shape)
        return res

    @property
    def p(self):
        """Partial pressure in the cells (atm)."""
        return self.u[:, 0]

    @property
    def T(self):
        """Temperature in the cells (K)."""
        return self.u[:, 1]

    def hot_spot(self):
        """The largest cell temperature (K)."""
        return float(self.T.max())


def _hot_spot(case, p0, Tw, **kw):
    return RunawayTube(case, p0, Tw, **kw).solve().hot_spot()


def critical_p0(case, Tw, lo, hi, target=None, xtol=1e-12, grow=1.01, max_grow=30, **kw):
    """Inlet pressure whose reactor hot spot equals ``target`` (default ``T_M``).

    The hot spot rises monotonically with the inlet pressure, so the critical
    value is the root of ``hot_spot(p0) - target``, found by ``brentq``; every
    evaluation is an independent solve from the inlet state. ``hi`` is an
    estimate (the paper's Eq. 30 mean is not a bound): while the hot spot at
    ``hi`` is still below the target, ``hi`` is raised by the factor ``grow``,
    at most ``max_grow`` times. Keep ``hi`` close: far above the critical value
    Newton from the inlet state does not converge (``RuntimeError``).

    Parameters
    ----------
    case : VWFCase
    Tw : float
        Wall (and inlet) temperature (K).
    lo, hi : float
        Lower end, below the root, and a first upper end (atm), for instance
        the paper's Eq. 27 and Eq. 30 values.
    target : float, optional
        Hot-spot temperature (K); ``T_M`` of Eq. 8 by default.
    xtol : float, optional
        Absolute tolerance on the inlet pressure (atm).
    grow, max_grow : float, int, optional
        Factor and number of steps by which ``hi`` may be raised.
    **kw
        Passed to :class:`RunawayTube` (``n_z``, ``length``).

    Raises
    ------
    ValueError
        When ``lo`` is not below the root, or ``hi`` is not above it after
        ``max_grow`` raises.
    """
    target = t_max_peak(case, Tw) if target is None else target

    def gap(p0):
        return _hot_spot(case, p0, Tw, **kw) - target

    if gap(lo) >= 0.0:
        raise ValueError(f"critical_p0: the hot spot at lo = {lo:g} atm already reaches "
                         f"{target:g} K")
    for _ in range(max_grow + 1):
        if gap(hi) > 0.0:
            break
        lo, hi = hi, hi * grow
    else:
        raise ValueError(f"critical_p0: no upper end found below {hi:g} atm "
                         f"({max_grow} raises by {grow:g})")
    return brentq(gap, lo, hi, xtol=xtol, rtol=4 * np.finfo(float).eps)


def critical_p0_bisect(case, Tw, lo, hi, target=None, tol=1e-7, **kw):
    """:func:`critical_p0` by bisection, to a bracket narrower than ``tol``.

    Returns
    -------
    float
        The midpoint of the final bracket (atm).
    """
    target = t_max_peak(case, Tw) if target is None else target
    while hi - lo >= tol:
        mid = 0.5 * (lo + hi)
        if _hot_spot(case, mid, Tw, **kw) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
