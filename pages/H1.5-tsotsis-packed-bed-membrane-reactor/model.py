"""Tsotsis et al. (1992) packed-bed catalytic membrane reactor, isothermal version.

Model code of gallery page H1.5. Two cocurrent plug-flow channels on the
dimensionless length ``zeta`` in ``[0, 1]``: the tube (feed side, catalyst bed)
and the shell (permeate side, argon sweep), four molar flows each
(C2H6, C2H4, H2, Ar), scaled by the total feed flow. Every species permeates by
Knudsen diffusion through the permselective layer (the paper's eq. 1 without
reaction in the membrane), ethane dehydrogenation runs in the tube (eq. 3) and,
in a named variant, also in the shell (eq. 5). The paper prints no kinetics and
no permeances: the permeation number ``Pi`` and the Damkohler number ``Da`` are
the two free scales.

Contents
--------
Thermo, ln_K, delta_H
    Standard-state equilibrium constants from tabulated ideal-gas heat
    capacities and formation properties (closed-form integrals).
xeq, k_of_xeq, smr_equilibrium
    Equilibrium conversion of a feed (E11), its inverse, and the steam
    reforming equilibrium of Fig. 4 (E13).
Case
    Feed, sweep, pressures, equilibrium constant, permeance ratios, shell
    reaction switch.
source
    Pointwise right-hand side of E8 for the eight flows (shared by both routes).
conversion
    The conversion definitions D1, D2, D4 (E14).
MembraneReactor
    The two channels as one pymrm finite-volume system solved by ``newton``.
solve_ivp_reactor, drain_function_ivp
    The same model as an initial-value problem (Radau), an independent route.
LocalEquilibrium
    The Da -> infinity limit (E16) in the autonomous variable ``s = Pi zeta``.
packed_bed_x, exchange_tube_ethane, relaxation_factor, knudsen_d
    Closed forms and quadratures of E15, E17 and E10.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from pymrm import (BlockLayout, NumJac, compute_boundary_values,
                   construct_convflux_upwind, construct_div, newton)

__all__ = ["R", "PSI", "ATM", "BAR", "T0", "SPECIES", "NU", "molar_masses", "q_knudsen", "Thermo",
           "ln_K", "delta_H", "xeq", "k_of_xeq", "smr_equilibrium", "Case", "source", "conversion",
           "MembraneReactor", "solve_ivp_reactor", "drain_function_ivp", "LocalEquilibrium",
           "packed_bed_x", "exchange_tube_ethane", "relaxation_factor", "knudsen_d", "replace"]

R = 8.314462618          # J/(mol K)
PSI = 6894.757293168     # Pa per psi
ATM = 101325.0           # Pa
BAR = 1.0e5              # Pa, standard state of K_p
T0 = 298.15              # K
SPECIES = ("C2H6", "C2H4", "H2", "Ar")
NU = np.array([-1.0, 1.0, 1.0, 0.0])
DEHYD = {"C2H6": -1.0, "C2H4": 1.0, "H2": 1.0}
SMR = {"CH4": -1.0, "H2O": -1.0, "CO": 1.0, "H2": 3.0}
WGS = {"CO": -1.0, "H2O": -1.0, "CO2": 1.0, "H2": 1.0}


def molar_masses(M_C=12.011, M_H=1.008, M_Ar=39.948):
    """Molar masses (g/mol) of C2H6, C2H4, H2, Ar from atomic weights."""
    return np.array([2 * M_C + 6 * M_H, 2 * M_C + 4 * M_H, 2 * M_H, M_Ar])


def q_knudsen(M=None):
    """Knudsen permeance ratios ``q_j = (M_E/M_j)^(1/2)``, ethane = 1 (A2)."""
    M = molar_masses() if M is None else M
    return np.sqrt(M[0] / M)


@dataclass(frozen=True)
class Thermo:
    """Tabulated ideal-gas data: ``{species: (A, B, C, D, dHf298, dGf298)}``.

    ``cp=False`` drops the heat-capacity terms (``dH`` constant at its 298.15 K
    value): a defect for the break table, not a model.
    """

    table: dict
    cp: bool = True

    def coeffs(self, nu):
        """Reaction sums ``(dA, dB, dC, dD, dH0, dG0)`` for stoichiometry ``nu``."""
        s = np.zeros(6)
        for sp, n in nu.items():
            s += n * np.asarray(self.table[sp], dtype=float)
        if not self.cp:
            s[:4] = 0.0
        return s


def _idcph(T, dA, dB, dC, dD):
    """``int_T0^T dCp/R dT`` in closed form (K)."""
    return (dA * (T - T0) + dB / 2 * (T ** 2 - T0 ** 2) + dC / 3 * (T ** 3 - T0 ** 3)
            + dD * (1.0 / T0 - 1.0 / T))


def _idcps(T, dA, dB, dC, dD):
    """``int_T0^T dCp/(R T) dT`` in closed form (dimensionless)."""
    return (dA * np.log(T / T0) + dB * (T - T0) + dC / 2 * (T ** 2 - T0 ** 2)
            + dD / 2 * (1.0 / T0 ** 2 - 1.0 / T ** 2))


def delta_H(T, thermo, nu=DEHYD):
    """Standard reaction enthalpy at ``T`` (J/mol), E12."""
    dA, dB, dC, dD, dH0, _ = thermo.coeffs(nu)
    return dH0 + R * _idcph(T, dA, dB, dC, dD)


def ln_K(T, thermo, nu=DEHYD):
    """``ln K`` (standard state 1 bar) at ``T``, E12: ``dG = dH0 - T (dH0 - dG0)/T0 + R IDCPH - R T IDCPS``."""
    dA, dB, dC, dD, dH0, dG0 = thermo.coeffs(nu)
    dG = dH0 - T / T0 * (dH0 - dG0) + R * _idcph(T, dA, dB, dC, dD) - R * T * _idcps(T, dA, dB, dC, dD)
    return -dG / (R * T)


def xeq(K, feed, P, Pr=0.0):
    """Equilibrium conversion of the ethane in ``feed`` (E11), at total pressure ``P`` (Pa).

    ``Pr`` mol argon per mol feed is added (the sweep-diluted reference). With the
    extent ``xi = y_E0 x``: ``(1 + K') xi^2 + (p (y_Y + y_H) + K'(1 + Pr - y_E)) xi
    + p y_Y y_H - K y_E (1 + Pr) = 0`` with ``p = P/1 bar``; one positive root,
    taken in the cancellation-free form.
    """
    yE, yY, yH, _ = feed
    p = P / BAR
    a = p + K
    b = p * (yY + yH) + K * (1.0 + Pr - yE)
    c = p * yY * yH - K * yE * (1.0 + Pr)
    disc = np.sqrt(b * b - 4 * a * c)
    xi = (-b + disc) / (2 * a) if b < 0 else -2 * c / (b + disc)
    return float(xi / yE)


def k_of_xeq(x, feed, P, Pr=0.0):
    """Inverse of ``xeq``: the ``K`` (bar) whose equilibrium conversion is ``x``."""
    yE, yY, yH, _ = feed
    xi = yE * x
    return float((yY + xi) * (yH + xi) * P / BAR / ((yE - xi) * (1.0 + Pr + xi)))


def smr_equilibrium(T, P, thermo, feed=(1.0, 4.0, 0.20)):
    """Steam reforming equilibrium (E13) per mol CH4 fed: ``(X_CH4, X_CO2)``.

    Feed ``CH4 : H2O : H2``; extents ``a`` (CH4 + H2O = CO + 3 H2) and ``b``
    (CO + H2O = CO2 + H2). For a given ``a`` the shift equilibrium is a quadratic
    in ``b`` with one root in ``[0, a]``; the reforming residual then decreases
    strictly in ``a`` on ``(0, 1)`` and ``brentq`` finds the one root.
    """
    c0, w0, h0 = feed
    K1 = np.exp(ln_K(T, thermo, SMR))
    K2 = np.exp(ln_K(T, thermo, WGS))
    p = P / BAR

    def b_of(a):
        # K2 (a - b)(w0 - a - b) = b (h0 + 3a + b)
        qa = K2 - 1.0
        qb = -K2 * (a + w0 - a) - (h0 + 3 * a)
        qc = K2 * a * (w0 - a)
        if abs(qa) < 1e-14:
            return -qc / qb
        disc = np.sqrt(qb * qb - 4 * qa * qc)
        r1, r2 = (-qb - disc) / (2 * qa), (-qb + disc) / (2 * qa)
        return r1 if 0.0 <= r1 <= a else r2

    def f(a):
        b = b_of(a)
        n = np.array([c0 - a, w0 - a - b, a - b, b, h0 + 3 * a + b])
        y = n / n.sum()
        return np.log(K1) - (np.log(y[2]) + 3 * np.log(y[4]) - np.log(y[0]) - np.log(y[1]) + 2 * np.log(p))

    a = brentq(f, 1e-12, c0 - 1e-12, xtol=1e-15, rtol=1e-14)
    return a / c0, b_of(a) / c0


@dataclass(frozen=True)
class Case:
    """Operating state of the isothermal model (E8).

    Parameters
    ----------
    feed : array
        Tube feed mole fractions ``(y_E, y_Y, y_H, y_Ar)``.
    Pr : float
        Sweep ratio, mol argon per mol feed (A3).
    PF, PP : float
        Tube and shell pressures (Pa, absolute).
    K : float
        Equilibrium constant (bar).
    q : array
        Permeance ratios to ethane (A2: Knudsen).
    shell_rx : bool
        Reaction also in the shell at the tube's Da (variant, A6).
    """

    feed: np.ndarray
    Pr: float
    PF: float
    PP: float
    K: float
    q: np.ndarray = field(default_factory=q_knudsen)
    shell_rx: bool = False

    @property
    def beta(self):
        """``P^P / P^F``."""
        return self.PP / self.PF

    @property
    def Kred(self):
        """``K_p / P^F`` (dimensionless)."""
        return self.K * BAR / self.PF

    def inlet(self):
        """Inlet flows ``(tube 4, shell 4)`` scaled by the total feed."""
        return np.concatenate([np.asarray(self.feed, float), [0.0, 0.0, 0.0, self.Pr]])


def source(n, Pi, Da, case):
    """Right-hand side of E8, ``dn/dzeta``, for flows ``n[..., 8]`` (tube 0:4, shell 4:8).

    Permeation ``J_j = Pi q_j (y_j^F - beta y_j^P)`` leaves the tube and enters
    the shell; the tube rate is ``Da (y_E - y_Y y_H / K')`` with ``K' = K_p/P^F``;
    the shell rate (variant) is ``Da (beta y_E - beta^2 y_Y y_H / K')``. Both
    rates conserve C and H atoms exactly, so the flows close the atom balances.
    """
    nF, nP = n[..., :4], n[..., 4:]
    yF = nF / nF.sum(axis=-1, keepdims=True)
    yP = nP / nP.sum(axis=-1, keepdims=True)
    J = Pi * case.q * (yF - case.beta * yP)
    r = Da * (yF[..., 0] - yF[..., 1] * yF[..., 2] / case.Kred)
    dF = -J + NU * r[..., None]
    dP = J.copy()
    if case.shell_rx:
        b = case.beta
        rs = Da * (b * yP[..., 0] - b * b * yP[..., 1] * yP[..., 2] / case.Kred)
        dP = dP + NU * rs[..., None]
    return np.concatenate([dF, dP], axis=-1)


def conversion(out, case, kind="D1"):
    """Conversion of the outlet flows ``out[8]`` (E14).

    ``D1`` ethylene leaving both sides over ethane fed; ``D2`` the tube outlet
    composition ``n_Y/(n_E + n_Y)``; ``D4`` ``1 - n_E^F/n_E0``.
    """
    yE0 = case.feed[0]
    if kind == "D1":
        return (out[1] + out[5]) / yE0
    if kind == "D2":
        return out[1] / (out[0] + out[1])
    if kind == "D4":
        return 1.0 - out[0] / yE0
    raise ValueError(f"conversion kind {kind!r}")


class MembraneReactor:
    """The two cocurrent channels as one steady pymrm finite-volume system.

    State ``(n_z, 8)``: ``pymrm.BlockLayout`` with domains ``"tube"`` and
    ``"shell"``, ``(n_z, 4)`` each, concatenated along the field axis. Each
    channel carries its molar flows with a constant unit-velocity first-order
    upwind operator (the balances are in flow form, so the state-dependent gas
    velocity never enters the operator); the source couples the channels in
    every cell and is linearised by ``NumJac`` (pointwise: the last axis in
    full). The operators are assembled once.

    Parameters
    ----------
    case : Case
    Pi, Da : float
    n_z : int
    """

    def __init__(self, case, Pi, Da, n_z=200):
        self.case, self.Pi, self.Da = case, float(Pi), float(Da)
        self.z_f = np.linspace(0.0, 1.0, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.layout = BlockLayout({"tube": (n_z, 4), "shell": (n_z, 4)}, axis=-1)
        self.shape = self.layout.shape
        self._build_operators()

    def _build_operators(self):
        u0 = self.case.inlet()
        # Each channel flows +zeta. Outward normal, a dn/dn_out + b n = d:
        #   zeta = 0 (inlet):  n = n_in (Dirichlet)  -> a = 0, b = 1, d = n_in
        #   zeta = 1 (outlet): nothing imposed       -> {"outflow": True}
        self.bc = {"tube": ({"a": 0.0, "b": 1.0, "d": u0[:4]}, {"outflow": True}),
                   "shell": ({"a": 0.0, "b": 1.0, "d": u0[4:]}, {"outflow": True})}
        place = self.layout.place
        jac = None
        g = np.zeros(self.shape)
        for name in ("tube", "shell"):
            shp = self.layout.shapes[name]
            conv, conv_bc = construct_convflux_upwind(shp, self.z_f, self.z_c, self.bc[name], v=1.0, axis=0)
            div = construct_div(shp, self.z_f, nu=0, axis=0)     # nu = 0: zeta is Cartesian
            blk = place(div @ conv, name, name)
            jac = blk if jac is None else jac + blk
            self.layout.put(g, name, (div @ conv_bc).toarray().reshape(shp))
        self.jac_const = jac
        self.g_const = g.reshape((-1, 1))
        self.numjac = NumJac(self.shape)          # pointwise source: last axis in full

    def residual(self, n):
        """Steady residual ``dn/dzeta - S(n)`` and its Jacobian."""
        g_s, jac_s = self.numjac(lambda v: source(v, self.Pi, self.Da, self.case), n)
        g = self.g_const + self.jac_const @ n.reshape((-1, 1)) - g_s.reshape((-1, 1))
        return g, self.jac_const - jac_s

    def initial_guess(self):
        """The inlet state in every cell."""
        return np.tile(self.case.inlet(), (self.shape[0], 1))

    def solve(self, n0=None, tol=1e-12, maxfev=80, lin_solver_kwargs=None):
        """Solve with ``pymrm.newton``; raises ``RuntimeError`` unless ``success``.

        ``tol`` is on flows of order one, above their round-off floor; the lower
        bound 0 keeps every flow, and so every mole fraction, physical.
        """
        n0 = self.initial_guess() if n0 is None else np.asarray(n0, float).reshape(self.shape)
        kw = {} if lin_solver_kwargs is None else dict(lin_solver_kwargs=lin_solver_kwargs)
        with np.errstate(invalid="ignore", divide="ignore"):
            res = newton(self.residual, n0, tol=tol, maxfev=maxfev, bounds=(0.0, None),
                         stagnation=True, **kw)
        if not res.success:
            raise RuntimeError(f"MembraneReactor Pi={self.Pi:g} Da={self.Da:g} n_z={self.shape[0]}: "
                               f"{res.message}")
        self.n, self.result = res.x.reshape(self.shape), res
        return self

    def outlet(self, read="face"):
        """Outlet flows ``[8]`` read on the outlet faces with ``compute_boundary_values``.

        ``read="cell"`` takes the last cells instead (pitfall P8; with the
        outflow outlet the two coincide).
        """
        if read == "cell":
            return self.n[-1].copy()
        out = []
        for name in ("tube", "shell"):
            v = self.layout.view(self.n, name)
            _, _, v1, _ = compute_boundary_values(v, self.z_f, self.z_c, self.bc[name])
            out.append(np.ravel(v1))
        return np.concatenate(out)

    def conversion(self, kind="D1", read="face"):
        """Conversion of the outlet flows (``conversion``)."""
        return conversion(self.outlet(read), self.case, kind)

    def atom_balances(self):
        """Relative closure of C, H and Ar between inlet and outlet faces."""
        u0, u1 = self.case.inlet(), self.outlet()
        C = np.array([2, 2, 0, 0, 2, 2, 0, 0.0])
        H = np.array([6, 4, 2, 0, 6, 4, 2, 0.0])
        A = np.array([0, 0, 0, 1, 0, 0, 0, 1.0])
        return np.array([abs(w @ u1 - w @ u0) / (w @ u0) for w in (C, H, A)])


def solve_ivp_reactor(case, Pi, Da, rtol=1e-11, events=None, dense=False):
    """E8 as an initial-value problem in ``zeta`` (Radau): the route with no pymrm.

    Returns the ``solve_ivp`` result; raises ``RuntimeError`` on failure.
    """
    s = solve_ivp(lambda z, u: source(u, Pi, Da, case), (0.0, 1.0), case.inlet(), method="Radau",
                  rtol=rtol, atol=1e-3 * rtol, events=events, dense_output=dense)
    if s.status not in (0, 1):
        raise RuntimeError(f"solve_ivp_reactor Pi={Pi:g} Da={Da:g}: {s.message}")
    return s


def drain_function_ivp(case, Pi, Da, eps, rtol=1e-11):
    """Continuous drainage function of ``Pi`` (E18): tube outflow minus ``eps``.

    When the tube flow reaches ``eps`` at ``zeta_e < 1`` the integration stops
    there and the function returns ``-(1 - zeta_e)``; both branches vanish at
    ``zeta_e = 1``, so the function is continuous through the threshold.
    """
    def ev(z, u):
        return u[:4].sum() - eps
    ev.terminal, ev.direction = True, -1
    s = solve_ivp_reactor(case, Pi, Da, rtol=rtol, events=ev)
    if s.status == 1:
        return -(1.0 - s.t_events[0][0])
    return s.y[:4, -1].sum() - eps


def _eq_state(a, b, c, Kp):
    """Equilibrium flows ``(E, Y, H, Ar)`` from invariants ``a = n_E + n_Y``,
    ``b = n_H - n_Y``, ``c = n_Ar`` at ``Kp = K_p/P``.

    ``xi (b + xi) = Kp (a - xi)(S + xi)``, ``S = a + b + c``:
    ``(1 + Kp) xi^2 + (b + Kp (b + c)) xi - Kp a S = 0``, positive root.
    """
    S = a + b + c
    B = b + Kp * (b + c)
    C = Kp * a * S
    xi = 2 * C / (B + np.sqrt(B * B + 4 * (1 + Kp) * C)) if B > 0 else (-B + np.sqrt(B * B + 4 * (1 + Kp) * C)) / (2 * (1 + Kp))
    return np.array([a - xi, xi, b + xi, c])


class LocalEquilibrium:
    """The Da -> infinity limit (E16): the tube (and, with ``shell_rx``, the shell)
    at chemical equilibrium in every section.

    The model is autonomous in ``s = Pi zeta``: one trajectory gives the outlet
    for every ``Pi`` (``outlet(Pi)``), and drainage is the ``s`` at which the tube
    flow falls to ``eps``. State: tube invariants ``(a, b, c)``; shell flows (4)
    or shell invariants (3).

    Parameters
    ----------
    case : Case
    s_max : float
        End of the trajectory.
    eps : float or None
        Tube flow at which the trajectory stops (drainage threshold).
    rtol : float
    """

    def __init__(self, case, s_max=10.0, eps=None, rtol=1e-12):
        self.case = case
        c = case
        u0 = c.inlet()
        tube0 = [u0[0] + u0[1], u0[2] - u0[1], u0[3]]
        shell0 = [u0[4] + u0[5], u0[6] - u0[5], u0[7]] if c.shell_rx else list(u0[4:])
        self.KpF, self.KpP = c.K * BAR / c.PF, c.K * BAR / c.PP
        ev = None
        if eps is not None:
            def ev(s, v):
                return self.flows(v)[:4].sum() - eps
            ev.terminal, ev.direction = True, -1
        self.sol = solve_ivp(self.rhs, (0.0, s_max), np.array(tube0 + shell0), method="LSODA",
                             rtol=rtol, atol=1e-3 * rtol, events=ev, dense_output=True)
        if self.sol.status not in (0, 1):
            raise RuntimeError(f"LocalEquilibrium: {self.sol.message}")
        self.s_end = float(self.sol.t[-1])
        self.s_drain = float(self.sol.t_events[0][0]) if (eps is not None and self.sol.status == 1) else None

    def flows(self, v):
        """The eight flows from the state."""
        nF = _eq_state(v[0], v[1], v[2], self.KpF)
        nP = _eq_state(v[3], v[4], v[5], self.KpP) if self.case.shell_rx else v[3:7]
        return np.concatenate([nF, nP])

    def rhs(self, s, v):
        n = self.flows(v)
        c = self.case
        J = c.q * (n[:4] / n[:4].sum() - c.beta * n[4:] / n[4:].sum())
        dt = [-(J[0] + J[1]), -(J[2] - J[1]), -J[3]]
        dsh = [J[0] + J[1], J[2] - J[1], J[3]] if c.shell_rx else list(J)
        return np.array(dt + dsh)

    def outlet(self, Pi):
        """Outlet flows ``[8]`` at permeation number ``Pi`` (``Pi <= s_end``)."""
        if Pi > self.s_end * (1 + 1e-12):
            raise ValueError(f"Pi = {Pi:g} beyond the trajectory end {self.s_end:g}")
        return self.flows(self.sol.sol(min(Pi, self.s_end)))

    def conversion(self, Pi, kind="D1"):
        return conversion(self.outlet(Pi), self.case, kind)


def packed_bed_x(case, Da):
    """Pi = 0 (E17 i): the conversion at ``zeta = 1`` from ``zeta(xi) = int dxi/(Da g)``.

    ``g(xi)`` is the tube rate over ``Da`` along the extent; the integral grows
    without bound as ``xi`` approaches equilibrium, so ``zeta(xi) = 1`` has one
    root below it.
    """
    yE, yY, yH, _ = case.feed
    Kr = case.Kred

    def g(xi):
        return ((yE - xi) - (yY + xi) * (yH + xi) / Kr / (1 + xi)) / (1 + xi)
    xi_eq = xeq(case.K, case.feed, case.PF) * yE

    def zeta(xi):
        return quad(lambda t: 1.0 / (Da * g(t)), 0.0, xi, epsabs=0, epsrel=1e-13, limit=200)[0]
    k = 1
    while zeta(xi_eq * (1 - 10.0 ** -k)) <= 1.0:
        k += 1
        if k > 12:
            raise ValueError("the packed bed reaches equilibrium before zeta = 1")
    return brentq(lambda xi: zeta(xi) - 1.0, 0.0, xi_eq * (1 - 10.0 ** -k), xtol=1e-16, rtol=1e-15) / yE


def exchange_tube_ethane(yE, Pi, Pr):
    """E17 ii (q = 1, Da = 0, beta = 1): tube ethane at ``zeta = 1``."""
    return (yE + Pr * yE * np.exp(-Pi * (1 + 1 / Pr))) / (1 + Pr)


def relaxation_factor(case):
    """E15: ``-dg/dxi`` at the tube equilibrium of the feed (no co-fed ethylene)."""
    yE, _, yH, _ = case.feed
    xi = xeq(case.K, case.feed, case.PF) * yE
    d1 = -(1 + yE) / (1 + xi) ** 2
    d2 = -((yH + 2 * xi) * (1 + xi) - 2 * xi * (yH + xi)) / (case.Kred * (1 + xi) ** 3)
    return -(d1 + d2)


def knudsen_d(d_pore, T, M_kg):
    """Knudsen diffusivity ``(d_p/3) (8 R T/(pi M))^(1/2)`` (m2/s), E10."""
    return d_pore / 3.0 * np.sqrt(8 * R * T / (np.pi * M_kg))
