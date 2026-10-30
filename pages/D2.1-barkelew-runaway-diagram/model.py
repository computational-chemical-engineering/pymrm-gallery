"""Barkelew's runaway diagram, as Van Welsenaere and Froment (1970) reprint it.

Model code of gallery page D2.1. The reduced (dimensionless) tube has two
parameters, the cooling group ``N`` and the adiabatic group ``S``; its hot spot
is computed two independent ways, by a pymrm finite-volume solve and by a
phase-plane quadrature that transcribes the reduced model separately. The
dimensional tube of their Eqs. 3-4 is here too, for the comparison with their
Figs. 3 and 10.

Contents
--------
Case
    The printed constants of the base case and the groups ``A``, ``B``, ``C``.
kexp, t_max_peak
    The Arrhenius factor and ``T_M`` of their Eq. 8.
n_group, s_group, p0_of_s
    Barkelew's groups and the inverse map from ``S`` to an inlet pressure.
tau_max_phase
    Hot spot of the reduced tube by quadrature in the phase plane.
BarkelewTube, tau_max_pymrm
    The reduced tube as a pymrm model, and its hot spot.
u_at_tau, tangency_at_h, tangency_richardson, tangency
    Barkelew's critical point: the tangency of an ``S``-curve with the envelope.
s_crit_tangency, s_crit_level
    The critical ``S`` for a given ``N``, by a direct root in ``ln S``.
p0_crit_vwf
    Van Welsenaere and Froment's criterion 1 with the true rate law, by
    back-integration in the phase plane.
Tube
    The dimensional tube as a pymrm model, true or exponential rate law.
profile_ivp, blowup_z, linearisations
    The dimensional tube marched as an initial-value problem (in ``z``, or in
    the phase plane), and the rate-law linearisations the printed text permits.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from pymrm import NumJac, construct_convflux_upwind, construct_div, newton

__all__ = ["Case", "kexp", "t_max_peak", "n_group", "s_group", "p0_of_s",
           "tau_max_phase", "BarkelewTube", "tau_max_pymrm", "u_at_tau",
           "tangency_at_h", "tangency_richardson", "tangency", "s_crit_tangency",
           "s_crit_level", "p0_crit_vwf", "Tube",
           "profile_ivp", "blowup_z", "linearisations", "TAU_CAP", "H_TAN"]

#: Cap on the dimensionless hot spot in the quadrature (deep runaway only).
TAU_CAP = 20.0
#: Coarse member of the Richardson pair ``(h, h/2)`` of the tangency.
H_TAN = 0.01


@dataclass(frozen=True)
class Case:
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
        Arrhenius constants of ``exp(-a/T + b)``.
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
        """Consumption group ``M P rho_b p_B0 / rho_g``."""
        return self.M * self.P * self.rho_b / self.rho_g * self.p_B0

    @property
    def B(self):
        """Heating group ``(-dH) rho_b p_B0 / c_p``."""
        return self.minus_dH * self.rho_b / self.c_p * self.p_B0

    @property
    def C(self):
        """Cooling group ``2 U / (c_p R)``. ``c_p`` is volumetric: no ``rho_g``."""
        return 2.0 * self.U / (self.c_p * self.R)


def kexp(case, T):
    """The true Arrhenius factor ``exp(-a/T + b)`` of their Eqs. 3-4."""
    return np.exp(case.b - case.a / np.asarray(T, float))


def t_max_peak(case, Tw):
    """Their Eq. 8: ``T_M``, the critical hot spot of criterion 1 (K)."""
    return 0.5 * (case.a - np.sqrt(case.a * (case.a - 4.0 * Tw)))


def n_group(case, Tw):
    """``N = C / (A k_w)``: cooling capacity over reaction rate. Independent of ``p0``."""
    return case.C / (case.A * kexp(case, Tw))


def s_group(case, p0, Tw):
    """``S = a (dT)_ad / T_w^2``, with ``(dT)_ad = (B/A) p0`` from their Eq. 32."""
    return case.a * case.B * p0 / (case.A * Tw ** 2)


def p0_of_s(case, S, Tw):
    """Inverse of :func:`s_group`: the inlet pressure (atm) with group ``S``."""
    return S * case.A * Tw ** 2 / (case.a * case.B)


def tau_max_phase(N, S, cap=TAU_CAP, rtol=1e-9, expo=np.exp, cool=1.0, order=1.0):
    """Hot spot of the reduced tube from the phase-plane form.

    ``dtau/dy = -S + N tau / (y^n e^tau)``, integrated from ``y = 1`` (the
    inlet) towards ``y = 0`` by LSODA; the hot spot is where ``dtau/dy = 0``.
    The independent variable is the reactant fraction, so this route forms no
    grid at all, and its right-hand side is written here separately from
    :meth:`BarkelewTube.reaction`: the two share no code.

    Parameters
    ----------
    N, S : float
        Cooling and adiabatic groups.
    cap : float, optional
        Stop at ``tau = cap`` and return ``cap`` (deep in the runaway region).
    rtol : float, optional
        Relative tolerance of the quadrature.
    expo, cool, order : optional
        Rate law ``y^order expo(tau)`` and the sign of the cooling term; the
        model is ``np.exp``, ``+1`` and ``1``. Other values are defects for
        the page's break table.

    Returns
    -------
    float
        ``tau_m``.

    Raises
    ------
    RuntimeError
        When the integrator stops before an event or the end of the domain.
    """
    def rhs(y, t):
        return [-S + cool * N * t[0] / (y ** order * expo(t[0]))]

    def flat(y, t):                      # dtau/dy = 0: the hot spot
        return S * y ** order * expo(t[0]) - cool * N * t[0]

    def hot(y, t):                       # deep-runaway cap
        return t[0] - cap

    flat.direction, flat.terminal = -1.0, True
    hot.direction, hot.terminal = 1.0, True
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        sol = solve_ivp(rhs, [1.0, 1e-12], [0.0], rtol=rtol, atol=1e-12,
                        events=[flat, hot], method="LSODA")
    if sol.status < 0:
        raise RuntimeError(f"tau_max_phase at N = {N:g}, S = {S:g}: the quadrature "
                           f"stopped: {sol.message}")
    if sol.t_events[0].size:
        return float(sol.y_events[0][0][0])
    if sol.t_events[1].size:
        return cap
    return float(sol.y[0].max())


class BarkelewTube:
    """The reduced tube on a ``zeta`` grid, state ``(n_z, 2) = [y, tau]``.

    ``dy/dzeta = -y e^tau``, ``dtau/dzeta = S y e^tau - N tau``, ``y(0) = 1``,
    ``tau(0) = 0``, with ``zeta = A k_w z`` the scaled contact time. Written as a
    steady convection problem with unit velocity: first-order upwind fluxes,
    ``construct_div``, the source by ``NumJac`` per cell, ``pymrm.newton``.

    Parameters
    ----------
    N, S : float
        Cooling and adiabatic groups.
    zeta_end : float
        Length of the domain in ``zeta``.
    n_z : int
        Number of cells.
    nu : int, optional
        Geometry of ``construct_div``; ``0`` (Cartesian) for a contact time.
        Other values are a defect for the page's break table.
    cool, s_fac : float, optional
        Sign of the cooling term and a factor on ``S`` in :meth:`reaction`
        only; ``1`` in the model, other values are defects for the break table.
    """

    def __init__(self, N, S, zeta_end, n_z, nu=0, cool=1.0, s_fac=1.0):
        self.N, self.S = float(N), float(S)
        self.cool, self.s_fac = float(cool), float(s_fac)
        self.nu = nu
        self.shape = (int(n_z), 2)
        self._build_grid(zeta_end)
        self._build_operators()
        self.u = np.tile(np.array([1.0, 0.0]), (self.shape[0], 1))

    def _build_grid(self, zeta_end):
        self.z_f = np.linspace(0.0, zeta_end, self.shape[0] + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])

    def _build_operators(self):
        # zeta = 0, outward normal -zeta: y = 1, tau = 0         -> a = 0, b = 1, d = (1, 0)
        # zeta = end, outward normal +zeta: dy/dzeta = dtau/dzeta = 0 (outflow)
        #                                                          -> a = 1, b = 0, d = 0
        bc = ({"a": 0.0, "b": 1.0, "d": np.array([1.0, 0.0])},
              {"a": 1.0, "b": 0.0, "d": 0.0})
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                      bc, v=1.0, axis=0)
        # v = 1: zeta is a scaled contact time, so d(v x)/dzeta = source.
        # nu = 0: Cartesian (a contact time has no curvature).
        div_mat = construct_div(self.shape, self.z_f, nu=self.nu, axis=0)
        self.jac_const = (div_mat @ conv_mat).tocsc()
        self.g_const = div_mat @ conv_bc
        self.numjac = NumJac(self.shape)          # pointwise source: last axis in full

    def reaction(self, u):
        """Source ``(-r, S r - N tau)`` with ``r = y e^tau``."""
        y, t = u[..., 0], u[..., 1]
        with np.errstate(over="ignore"):          # an overflowing trial step is
            r = y * np.exp(t)                     # shortened by newton's damping
        return np.stack([-r, self.s_fac * self.S * r - self.cool * self.N * t], axis=-1)

    def residual(self, u):
        """Steady residual and Jacobian: convection minus source."""
        g_rxn, jac_rxn = self.numjac(self.reaction, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self, tol=1e-10, maxfev=100):
        """Newton from the inlet state everywhere; ``y`` and ``tau`` are of order one.

        Raises
        ------
        RuntimeError
            When Newton stops without converging (``NewtonBreakdown``, a
            singular Jacobian or a non-finite step, propagates as raised).
        """
        result = newton(self.residual, self.u, tol=tol, maxfev=maxfev)
        if not result.success:
            raise RuntimeError(f"BarkelewTube.solve at N = {self.N:g}, S = {self.S:g}: "
                               f"{result.message}")
        self.u, self.result = result.x.reshape(self.shape), result
        return self

    def hot_spot(self):
        """The largest cell value of ``tau``."""
        return float(self.u[:, 1].max())


def tau_max_pymrm(N, S, n_z=1600, nu=0, cool=1.0, s_fac=1.0, return_tube=False):
    """Hot spot of the reduced tube by pymrm, in two passes.

    A 200-cell solve locates the peak (the domain doubles while the peak sits
    in the last tenth of it); the ``n_z``-cell solve then covers three times the
    peak position. Both are independent solves from the inlet state; nothing is
    warm-started.

    Returns
    -------
    float or (float, BarkelewTube)
        ``tau_m``, and the fine tube when ``return_tube`` is true.

    Raises
    ------
    RuntimeError, pymrm.NewtonBreakdown
        When a solve does not converge.
    """
    kw = dict(nu=nu, cool=cool, s_fac=s_fac)
    zeta_end = 1.0
    for _ in range(8):
        coarse = BarkelewTube(N, S, zeta_end, 200, **kw).solve()
        i = int(np.argmax(coarse.u[:, 1]))
        if i < 180:
            break
        zeta_end *= 2.0
    fine = BarkelewTube(N, S, max(3.0 * coarse.z_c[i], 1e-3), n_z, **kw).solve()
    return (fine.hot_spot(), fine) if return_tube else fine.hot_spot()


def u_at_tau(S, target, tm=tau_max_phase, lo=1e-4, hi=40.0):
    """``N/S`` at which the hot spot equals ``target``; ``tau_m`` falls with ``N/S``."""
    return brentq(lambda u: tm(u * S, S) - target, lo, hi, xtol=1e-10)


def tangency_at_h(S, tm=tau_max_phase, h=0.04, bracket=None):
    """Barkelew's critical point at one step ``h`` in ``ln S``, unextrapolated.

    The envelope of a one-parameter family is where ``d/dS`` vanishes at fixed
    abscissa, so the tangency on curve ``S`` is the ``N/S`` at which the curves
    for ``S e^h`` and ``S e^-h`` cross: a centred difference in ``ln S`` with an
    ``O(h^2)`` truncation error. The root is bracketed by the hot spots
    ``tau_m = 1.7`` and ``0.7``, and the sign change is tested, not assumed.

    ``bracket`` replaces that default bracket on ``N/S`` (a narrower one is
    cheaper; the sign change is tested all the same).

    Returns
    -------
    (float, float) or (None, None)
        ``N/S`` and ``tau_m`` at the tangency; ``None`` when the bracket holds
        no sign change (no envelope).
    """
    Sp, Sm = S * np.exp(h), S * np.exp(-h)

    def f(u):
        return tm(u * Sp, Sp) / Sp - tm(u * Sm, Sm) / Sm

    lo, hi = (u_at_tau(S, 1.7, tm), u_at_tau(S, 0.7, tm)) if bracket is None else bracket
    if f(lo) * f(hi) >= 0:
        return None, None
    ur = brentq(f, lo, hi, xtol=1e-8)
    return ur, tm(ur * S, S)


def tangency_richardson(S, tm=tau_max_phase, h=H_TAN):
    """Richardson extrapolation of :func:`tangency_at_h` from ``(h, h/2)`` to ``h -> 0``."""
    u1, t1 = tangency_at_h(S, tm, h)
    if u1 is None:
        return None, None
    u2, t2 = tangency_at_h(S, tm, 0.5 * h)
    if u2 is None:
        return None, None
    return (4.0 * u2 - u1) / 3.0, (4.0 * t2 - t1) / 3.0


@lru_cache(maxsize=None)
def _tangency_cached(S, h):
    return tangency_richardson(S, tau_max_phase, h)


def tangency(S, tm=None, h=H_TAN):
    """The extrapolated tangency ``(N/S, tau_m)`` on the curve ``S``.

    With ``tm=None`` (the model's own quadrature) the value is a pure function
    of ``(S, h)`` and is memoised; another ``tm`` (a defect) is not cached.
    """
    if tm is None:
        return _tangency_cached(float(S), float(h))
    return tangency_richardson(S, tm, h)


def _root_in_log_s(f, S_lo, S_hi, xtol=1e-10, grow=2.0, max_grow=12):
    """Root of ``f(ln S)`` between ``S_lo`` and ``S_hi``, widened until it brackets."""
    a, b = np.log(S_lo), np.log(S_hi)
    fa, fb = f(a), f(b)
    for _ in range(max_grow):
        if fa * fb < 0:
            break
        a, b = a - np.log(grow), b + np.log(grow)
        fa, fb = f(a), f(b)
    else:
        raise ValueError(f"no sign change between S = {np.exp(a):g} and {np.exp(b):g}")
    return float(np.exp(brentq(f, a, b, xtol=xtol)))


def s_crit_tangency(N, S_lo, S_hi, tm=tau_max_phase, h=H_TAN, u_guess=None, xtol=1e-9):
    """Barkelew's critical ``S`` for a given ``N``: a direct root, no interpolation.

    Solves ``S u_t(S) = N`` in ``ln S``, with ``u_t(S)`` the extrapolated
    tangency ``N/S`` of :func:`tangency_richardson` (``N u_t`` rises with ``S``).
    ``S_lo, S_hi`` bracket the root (widened if needed), for instance two
    neighbouring curves of a table. ``u_guess(S)``, an estimate of ``u_t(S)``,
    only narrows the inner bracket of each tangency search; the bracket is
    tested and the full one used when it holds no sign change.
    """
    def u_t(S):
        def tan_h(hh):
            if u_guess is not None:
                g = u_guess(S)
                u, _ = tangency_at_h(S, tm, hh, bracket=(g * np.exp(-0.03), g * np.exp(0.03)))
                if u is not None:
                    return u
            return tangency_at_h(S, tm, hh)[0]
        return (4.0 * tan_h(0.5 * h) - tan_h(h)) / 3.0

    return _root_in_log_s(lambda lnS: np.exp(lnS) * u_t(np.exp(lnS)) - N, S_lo, S_hi,
                          xtol=xtol)


def s_crit_level(N, S_lo, S_hi, target=1.0, tm=tau_max_phase):
    """``S`` at which the hot spot equals ``target`` for a given ``N`` (criterion 1
    at ``target = 1``): a direct root in ``ln S``; ``tau_m`` rises with ``S``."""
    return _root_in_log_s(lambda lnS: tm(N, np.exp(lnS)) - target, S_lo, S_hi)


def p0_crit_vwf(case, Tw, Bs=1.0, Cs=1.0):
    """Criterion 1 with the true rate law: the critical inlet pressure (atm).

    Back-integrates their Eq. 5, ``dT/dp``, from the top of the maxima curve
    ``(T_M, p_M)`` until ``T`` falls to ``T_w``, the inlet. ``Bs`` and ``Cs``
    scale ``B`` and ``C`` (defects for the page's blindness table; 1 in the
    model).
    """
    A, B, C = case.A, case.B * Bs, case.C * Cs
    TM = t_max_peak(case, Tw)
    pM = (TM - Tw) / ((B / C) * kexp(case, TM))          # their Eq. 7 at T = T_M

    def dTdp(p, T):                                      # their Eq. 5
        return [-B / A + (C / A) * (T[0] - Tw) / (p * kexp(case, T[0]))]

    def hit_Tw(p, T):
        return T[0] - Tw
    hit_Tw.terminal, hit_Tw.direction = True, -1.0
    sol = solve_ivp(dTdp, [pM, pM + 1.0], [TM], rtol=1e-10, atol=1e-12,
                    events=hit_Tw, method="LSODA")
    return float(sol.t_events[0][0])


class Tube:
    """The dimensional tube of their Eqs. 3-4 on the contact time ``z = z'/u`` (hr).

    ``dp/dz = -A p k(T)``, ``dT/dz = B p k(T) - C (T - T_w)``, inlet
    ``p = p0``, ``T = T_w``. State ``(n_z, 2) = [p, T]``. ``simplified=True``
    replaces ``k`` by Barkelew's ``k_w exp(a (T - T_w) / T_w^2)``.

    Parameters
    ----------
    case : Case
    p0 : float
        Inlet partial pressure (atm).
    Tw : float
        Wall and inlet temperature (K).
    length : float, optional
        Tube length (m).
    n_z : int, optional
        Number of cells.
    simplified : bool, optional
        Barkelew's exponential rate law instead of the Arrhenius factor.
    """

    #: ``newton`` bounds: ``T > 0``, where the Arrhenius factor is defined.
    BOUNDS = (np.array([-np.inf, 0.0]), None)

    def __init__(self, case, p0, Tw, length=1.0, n_z=2000, simplified=False):
        self.case, self.p0, self.Tw, self.simplified = case, float(p0), float(Tw), simplified
        self.kw = float(kexp(case, Tw))
        self.shape = (int(n_z), 2)
        self.z_f = np.linspace(0.0, length / case.u, self.shape[0] + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self._build_operators()
        self.u = np.tile(np.array([self.p0, self.Tw]), (self.shape[0], 1))

    def _build_operators(self):
        # z = 0, outward normal -z: p = p0, T = T_w             -> a = 0, b = 1, d = (p0, T_w)
        # z = L, outward normal +z: dp/dz = dT/dz = 0 (outflow) -> a = 1, b = 0, d = 0
        bc = ({"a": 0.0, "b": 1.0, "d": np.array([self.p0, self.Tw])},
              {"a": 1.0, "b": 0.0, "d": 0.0})
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                      bc, v=1.0, axis=0)
        # v = 1 on the contact time; nu = 0: Cartesian.
        div_mat = construct_div(self.shape, self.z_f, nu=0, axis=0)
        self.jac_const = (div_mat @ conv_mat).tocsc()
        self.g_const = div_mat @ conv_bc
        self.numjac = NumJac(self.shape)

    def k(self, T):
        """Rate factor: Arrhenius, or Barkelew's exponential about ``T_w``."""
        with np.errstate(over="ignore"):
            if self.simplified:
                return self.kw * np.exp(self.case.a * (T - self.Tw) / self.Tw ** 2)
            return kexp(self.case, T)

    def reaction(self, u):
        """Source ``(-A r, B r - C (T - T_w))`` with ``r = p k(T)``."""
        p, T = u[..., 0], u[..., 1]
        r = p * self.k(T)
        return np.stack([-self.case.A * r, self.case.B * r - self.case.C * (T - self.Tw)],
                        axis=-1)

    def residual(self, u):
        """Steady residual and Jacobian: convection minus source."""
        g_rxn, jac_rxn = self.numjac(self.reaction, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self, tol=1e-10, maxfev=50):
        """Newton from the inlet state, step test on ``u / (p0, T_w)``, ``T > 0``.

        Raises
        ------
        RuntimeError
            When Newton does not converge (``NewtonBreakdown`` propagates).
        """
        result = newton(self.residual, self.u, tol=tol, maxfev=maxfev,
                        scale=np.array([self.p0, self.Tw]), bounds=self.BOUNDS)
        if not result.success:
            raise RuntimeError(f"Tube.solve at p0 = {self.p0:g}, Tw = {self.Tw:g}: "
                               f"{result.message}")
        self.u, self.result = result.x.reshape(self.shape), result
        return self

    def hot_spot(self):
        """The largest cell temperature (K)."""
        return float(self.u[:, 1].max())


def _march(case, p0, Tw, k_of_T, length, T_stop, method, dense, rtol=1e-10):
    def rhs(z, w):
        r = w[0] * k_of_T(w[1])
        return [-case.A * r, case.B * r - case.C * (w[1] - Tw)]

    def hot(z, w):
        return w[1] - T_stop
    hot.terminal, hot.direction = True, 1.0
    return solve_ivp(rhs, [0.0, length / case.u], [p0, Tw], rtol=rtol, atol=1e-12,
                     events=hot, method=method, dense_output=dense)


def profile_ivp(case, p0, Tw, simplified, length=1.0, T_stop=800.0):
    """The dimensional tube marched as an initial-value problem (LSODA).

    Stops where ``T`` reaches ``T_stop`` (the top of the paper's axis).

    Returns
    -------
    (numpy.ndarray, numpy.ndarray, float or None)
        Positions ``z'`` (m), temperatures (K), and the ``z'`` at which
        ``T_stop`` is reached (``None`` if it is not).
    """
    kw = float(kexp(case, Tw))
    if simplified:
        def k_of_T(T):
            return kw * np.exp(case.a * (T - Tw) / Tw ** 2)
    else:
        def k_of_T(T):
            return np.exp(case.b - case.a / T)
    sol = _march(case, p0, Tw, k_of_T, length, T_stop, "LSODA", True)
    if sol.status < 0:
        raise RuntimeError(f"profile_ivp at T_w = {Tw:g}: {sol.message}")
    zz = np.linspace(0.0, sol.t[-1], 800)
    z_stop = float(sol.t_events[0][0] * case.u) if sol.t_events[0].size else None
    return zz * case.u, sol.sol(zz)[1], z_stop


def blowup_z(case, Tw, k_wall, slope, p0, length=1.0, T_stop=800.0, rtol=1e-10,
             route="phase"):
    """``z'`` at which ``T`` reaches ``T_stop`` for ``k = k_wall exp(slope (T - T_w))``.

    ``route="phase"`` (default) integrates in the phase plane with the partial
    pressure as the independent variable, ``dT/dp = -B/A + C (T - T_w) / (A p k)``
    and ``dz/dp = -1/(A p k)``, from ``p0`` downwards (Radau). Where ``k`` grows
    without bound (a very steep exponent ignites at the inlet) both right-hand
    sides stay finite, while the march in ``z`` needs a step below the spacing
    of floating-point numbers and stops. ``route="z"`` marches in ``z`` (Radau)
    instead; it is the second route of the page's checks.

    Parameters
    ----------
    rtol : float, optional
        Relative tolerance of the integrator.
    route : {"phase", "z"}, optional
        Independent variable of the integration.

    Returns
    -------
    (float or None, float)
        The ``z'`` (m), ``None`` without runaway inside ``length``, and the
        highest ``T`` (K) reached.

    Raises
    ------
    RuntimeError
        When the integrator stops before an event or the end of the domain.
    """
    def k_of_T(T):
        with np.errstate(over="ignore"):
            return k_wall * np.exp(slope * (T - Tw))

    if route == "z":
        sol = _march(case, p0, Tw, k_of_T, length, T_stop, "Radau", False, rtol)
        if sol.status < 0:
            raise RuntimeError(f"blowup_z (z) at T_w = {Tw:g}: {sol.message}")
        z = float(sol.t_events[0][0] * case.u) if sol.t_events[0].size else None
        return z, float(sol.y[1].max())

    A, B, C, z_end = case.A, case.B, case.C, length / case.u

    def rhs(p, w):
        k = k_of_T(w[0])
        return [-B / A + C * (w[0] - Tw) / (A * p * k), -1.0 / (A * p * k)]

    def hot(p, w):
        return w[0] - T_stop

    def tube_end(p, w):
        return w[1] - z_end
    hot.terminal, hot.direction = True, 1.0
    tube_end.terminal, tube_end.direction = True, 1.0
    sol = solve_ivp(rhs, [p0, p0 * 1e-12], [Tw, 0.0], method="Radau", rtol=rtol,
                    atol=[1e-9, 1e-16], events=[hot, tube_end])
    if sol.status != 1:
        raise RuntimeError(f"blowup_z (phase) at T_w = {Tw:g}: {sol.message}")
    z = float(sol.y_events[0][0][1] * case.u) if sol.t_events[0].size else None
    return z, float(sol.y[0].max())


def linearisations(case, Tw):
    """The linearisations of ``ln k`` that every printed sentence permits.

    Each is ``ln k(T) ~ ln k(T_r) + (T - T_r) a / T_r^2`` written as
    ``k_w,eff e^tau`` with ``tau = a (T - T_w) / D``: all vanish at the wall and
    are linear in ``T - T_w``; only ``D`` differs. The last one drops the square
    and is a defect, not a permitted reading.

    Returns
    -------
    dict
        ``name -> (k_wall, slope)``.
    """
    TM = float(t_max_peak(case, Tw))
    kw = float(kexp(case, Tw))
    return {
        "tangent at T_w  (this page)": (kw, case.a / Tw ** 2),
        "tangent at T_M": (float(kexp(case, TM)) * np.exp(-case.a * (TM - Tw) / TM ** 2),
                           case.a / TM ** 2),
        "chord T_w to T_M": (kw, case.a / (Tw * TM)),
        "square dropped, a(T-T_w)/T_w": (kw, case.a / Tw),
    }
