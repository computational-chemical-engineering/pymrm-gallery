"""Itoh (1987) palladium membrane reactor for cyclohexane dehydrogenation.

Model code of gallery page H1.4. Two plug-flow channels on dimensionless
length ``L`` in ``[0, 1]``: the reaction (tube) side carries cyclohexane ``u_C``
and hydrogen ``u_H``, the separation (shell) side carries the permeated
hydrogen ``v_H``. Benzene and both argon flows are algebraic. The paper's
Eqs. 1-8 as transcribed on the page; the equilibrium constant is the page's
reconstruction from the printed equilibrium conversion.

Contents
--------
ItohCase
    Printed constants at 473 K (Eqs. 1 and 7), the feed dilution and the
    reconstructed equilibrium constant.
kp_of_xeq, unconverted_at_equilibrium, xeq_of_kp
    Eq. 7's equilibrium condition for a membrane-free feed, and its inverse.
source
    Pointwise sources of ``(u_C, u_H, v_H)`` in mol/s: Eq. 7 and Eq. 1.
ceiling
    The kinetics-free co-current bound: infinite permeance, instantaneous
    reaction.
MembraneReactor
    The two channels as one pymrm finite-volume system, co-current or
    counter-current, solved by ``pymrm.newton``.
conversion_ivp
    The co-current model as an initial-value problem in ``L`` (Radau), an
    independent route.
root_log
    A root in a logarithmic variable whose bracket is widened until it
    brackets.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from pymrm import (BlockLayout, NumJac, compute_boundary_values,
                   construct_convflux_upwind, construct_div, newton)

__all__ = ["ItohCase", "kp_of_xeq", "unconverted_at_equilibrium", "xeq_of_kp", "source",
           "ceiling", "MembraneReactor", "conversion_ivp", "root_log"]


@dataclass(frozen=True)
class ItohCase:
    """Constants of Itoh (1987) at the one temperature of the model.

    Parameters
    ----------
    T : float
        Temperature (K).
    P0 : float
        Total pressure on both sides, ``P_Tr = P_Ts`` (Pa).
    y_C0 : float
        Feed mole fraction of cyclohexane in argon.
    alpha_H : float
        Permeation constant of Eq. 1 (mol/s), as printed.
    k_rate : float
        Rate constant of Eq. 7 at ``T`` (mol m-3 Pa-1 s-1).
    K_B : float
        Benzene adsorption constant of Eq. 7 at ``T`` (1/Pa).
    V_r : float
        Gross volume of the reaction section (m3).
    K_p : float
        Equilibrium constant at ``T`` (Pa3). The printed expression is not
        usable (page, misprint 2); pass the reconstruction from ``x_eq``.
    nu_H2 : float, optional
        Hydrogen formed per cyclohexane converted, 3 for the reaction.
    sqrt_eps : float, optional
        0 (default): Eq. 1 as printed. A positive value replaces ``sqrt(y)`` by
        ``y / sqrt(y + sqrt_eps)``, smooth and linear at ``y = 0``: a numerical
        check that the infinite slope of the square root at zero is not what
        stops a solve, not a model.
    """

    T: float
    P0: float
    y_C0: float
    alpha_H: float
    k_rate: float
    K_B: float
    V_r: float
    K_p: float
    nu_H2: float = 3.0
    sqrt_eps: float = 0.0

    @classmethod
    def printed(cls, T, P0, y_C0, K_p, alpha_H=4.47e-5, r_i=8.5e-3, l_0=0.14, **changes):
        """The case with Eq. 7's printed Arrhenius forms evaluated at ``T``.

        ``k = 0.221 exp(-4270/T)``, ``K_B = 2.03e-10 exp(6270/T)`` and
        ``V_r = pi r_i^2 l_0`` (the gross volume of the reaction section).
        """
        kw = dict(T=float(T), P0=float(P0), y_C0=float(y_C0), alpha_H=float(alpha_H),
                  k_rate=0.221 * np.exp(-4270.0 / T), K_B=2.03e-10 * np.exp(6270.0 / T),
                  V_r=np.pi * r_i ** 2 * l_0, K_p=float(K_p))
        kw.update(changes)
        return cls(**kw)

    @property
    def n_Ar(self):
        """Argon fed per mole of cyclohexane, ``(1 - y_C0)/y_C0``."""
        return (1.0 - self.y_C0) / self.y_C0


def kp_of_xeq(x, y_C0, P0):
    """Eq. 7's equilibrium condition ``K_p p_C / p_H^3 = p_B`` solved for ``K_p``.

    For one mole of cyclohexane with ``n_Ar = (1 - y_C0)/y_C0`` mol argon at
    conversion ``x`` and total pressure ``P0``:
    ``K_p = 27 x^4 P0^3 / ((1 - x)(1 + n_Ar + 3x)^3)``.
    """
    n_Ar = (1.0 - y_C0) / y_C0
    return 27.0 * x ** 4 * P0 ** 3 / ((1.0 - x) * (1.0 + n_Ar + 3.0 * x) ** 3)


def _log_kp_of_unconverted(log_y, y_C0, P0):
    """``ln kp_of_xeq`` written in the logarithm of the unconverted fraction ``y = 1 - x``."""
    n_Ar = (1.0 - y_C0) / y_C0
    y = np.exp(log_y)
    return (np.log(27.0) + 4.0 * np.log1p(-y) + 3.0 * np.log(P0) - log_y
            - 3.0 * np.log(4.0 + n_Ar - 3.0 * y))


def unconverted_at_equilibrium(K_p, y_C0, P0):
    """Membrane-free equilibrium ``1 - x_eq`` for a given ``K_p``.

    ``kp_of_xeq`` increases strictly on ``0 < x < 1``: its logarithmic
    derivative ``4/x + 1/(1 - x) - 9/(1 + n_Ar + 3x)`` is at least
    ``9 - 9/(1 + n_Ar) > 0``. The root is therefore unique, and it is found in
    ``ln y`` with ``y = 1 - x`` so that a conversion within 1e-16 of 1 is still
    resolved. The bracket ``y`` in ``[1e-300, 1 - 1e-15]`` covers every ``K_p``
    from about ``1e-50`` to ``1e300`` Pa3; outside it ``brentq`` raises.
    """
    def f(log_y):
        return _log_kp_of_unconverted(log_y, y_C0, P0) - np.log(K_p)
    return float(np.exp(brentq(f, np.log(1e-300), np.log1p(-1e-15), xtol=1e-14)))


def xeq_of_kp(K_p, y_C0, P0):
    """Membrane-free equilibrium conversion for ``K_p``; see ``unconverted_at_equilibrium``."""
    return 1.0 - unconverted_at_equilibrium(K_p, y_C0, P0)


def source(u, case, uC0, vA0):
    """Pointwise sources of ``(u_C, u_H, v_H)`` in mol/s per unit ``L``.

    Eq. 7 in its ``p_H^3``-multiplied form,
    ``r_C = -k (K_p p_C - p_B p_H^3) / (p_H^3 + K_B K_p p_C)``, which is smooth
    at the ``p_H = 0`` inlet (zero order far from equilibrium); the 1e-30 in
    the denominator only keeps ``0/0`` out of a state with neither cyclohexane
    nor hydrogen; it stays positive only for ``u_C >= 0``. Eq. 1 with
    ``P_Tr = P_Ts = P0``: square roots of the hydrogen mole fractions, defined
    for ``u_H >= 0`` and ``v_H >= 0``. The solver bounds all three below by 0.

    Parameters
    ----------
    u : numpy.ndarray
        State ``[..., (u_C, u_H, v_H)]`` in mol/s.
    case : ItohCase
    uC0, vA0 : float
        Cyclohexane feed and argon purge (mol/s).
    """
    uC, uH, vH = u[..., 0], u[..., 1], u[..., 2]
    uA = case.n_Ar * uC0                    # argon, constant (Eq. 6)
    uB = uC0 - uC                           # benzene closure (Eq. 4)
    su = uC + uB + uH + uA                  # total reaction-side flow
    sv = vA0 + vH                           # total separation-side flow
    P0 = case.P0
    pC, pB, pH = uC / su * P0, uB / su * P0, uH / su * P0
    rC = -case.k_rate * (case.K_p * pC - pB * pH ** 3) / (pH ** 3 + case.K_B * case.K_p * pC
                                                        + 1e-30)
    if case.sqrt_eps > 0.0:
        def root(y):
            return y / np.sqrt(y + case.sqrt_eps)
    else:
        root = np.sqrt
    perm = case.alpha_H * (root(uH / su) - root(vH / sv))
    return np.stack([rC * case.V_r, -case.nu_H2 * rC * case.V_r - perm, perm], axis=-1)


def ceiling(case, uC0, vA0, K_p=None):
    """Largest co-current conversion: infinite permeance, instantaneous reaction.

    At the outlet the hydrogen mole fraction ``y`` is then equal on both
    sides and the reaction is at equilibrium. With the hydrogen balance
    ``u_H + v_H = 3 x u_C0`` the equilibrium ``K_p p_C = p_B p_H^3`` reduces to
    ``K_p (1 - x) = x (y P0)^3``, ``y/(1 - y) = 3 x u_C0 / (u_C0 (1 + n_Ar) + v_A0)``.
    The left side falls and the right side rises with ``x``, so the root on
    ``(0, 1)`` is unique; ``brentq`` raises if the ends do not bracket it.
    """
    K_p = case.K_p if K_p is None else K_p
    F0 = uC0 * (1.0 + case.n_Ar)            # non-hydrogen reaction-side flow

    def f(x):
        q = 3.0 * x * uC0 / (F0 + vA0)      # y/(1 - y) for the common H2 fraction
        y = q / (1.0 + q)
        return K_p * (1.0 - x) - x * (y * case.P0) ** 3
    return brentq(f, 1e-12, 1.0 - 1e-12, xtol=1e-15)


def root_log(f, lo, hi, grow=10.0, max_grow=20, xtol=1e-12):
    """Root of ``f`` in ``ln p`` for a positive parameter ``p``, increasing ``f``.

    ``lo`` and ``hi`` are first guesses: while ``f(lo) > 0`` ``lo`` is divided
    by ``grow``, while ``f(hi) < 0`` ``hi`` is multiplied by it, at most
    ``max_grow`` times each. Raises ``ValueError`` when no bracket is found.
    """
    for _ in range(max_grow):
        if f(lo) <= 0.0:
            break
        lo /= grow
    else:
        raise ValueError(f"root_log: f > 0 down to {lo:g}")
    for _ in range(max_grow):
        if f(hi) >= 0.0:
            break
        hi *= grow
    else:
        raise ValueError(f"root_log: f < 0 up to {hi:g}")
    return float(np.exp(brentq(lambda s: f(np.exp(s)), np.log(lo), np.log(hi), xtol=xtol)))


class MembraneReactor:
    """Itoh's two channels as one steady pymrm finite-volume system.

    Unknowns are the flows scaled by the feed, ``w = (u_C, u_H, v_H)/u_C0``,
    order one. The state is ``(n_z, 3)``, the reaction-side fields
    ``(u_C, u_H)`` and the separation-side ``v_H`` concatenated along the
    field axis with ``pymrm.BlockLayout`` (domains ``"tube"`` ``(n_z, 2)`` and
    ``"shell"`` ``(n_z, 1)``). Each channel has its own first-order upwind
    convection operator, unit velocity along ``L`` (the balances are in
    flow-rate form) and its own boundary conditions; the blocks are placed
    once. The source couples the channels in every cell and is linearised by
    ``NumJac`` on the full state (pointwise: the last axis in full).

    Parameters
    ----------
    case : ItohCase
    uC0, vA0 : float
        Cyclohexane feed and argon purge (mol/s).
    n_z : int, optional
        Number of cells.
    flow : {"co", "counter"}, optional
        Purge entering at ``L = 0`` (the paper's arrangement) or at ``L = 1``.
    outlet : {"outflow", "zero_gradient"}, optional
        Boundary condition at each channel's outlet. ``"outflow"`` is the
        library's pure-outflow marker: the face carries the last cell's value,
        the scheme is exact backward Euler in ``L`` and the balances close.
        ``"zero_gradient"`` (``dw/dn = 0``) is kept only as a comparison: the
        upwind operator then extrapolates the outlet face from the last two
        cells.

    Attributes
    ----------
    layout : pymrm.BlockLayout
        Where the tube and shell fields sit in the state.
    jac_const : scipy.sparse.csc_array
        The constant convection operator of both channels.
    g_const : numpy.ndarray
        The inlet contribution, shape ``(n_z * 3, 1)``.
    """

    #: lower bound 0 on every flow: Eq. 1 takes square roots of u_H and v_H, and
    #: Eq. 7's denominator p_H^3 + K_B K_p p_C can reach zero for u_C < 0
    BOUNDS = (0.0, None)

    def __init__(self, case, uC0, vA0, n_z=400, flow="co", outlet="outflow"):
        if flow not in ("co", "counter") or outlet not in ("outflow", "zero_gradient"):
            raise ValueError(f"flow {flow!r} or outlet {outlet!r} not recognised")
        self.case, self.uC0, self.vA0 = case, float(uC0), float(vA0)
        self.flow, self.outlet = flow, outlet
        self.z_f = np.linspace(0.0, 1.0, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.layout = BlockLayout({"tube": (n_z, 2), "shell": (n_z, 1)}, axis=-1)
        self.shape = self.layout.shape
        self._build_operators()

    def _channel(self, shape, bc, v):
        """Convection operator of one channel: ``div @ conv`` and its inlet vector."""
        conv, conv_bc = construct_convflux_upwind(shape, self.z_f, self.z_c, bc, v=v, axis=0)
        div = construct_div(shape, self.z_f, nu=0, axis=0)   # nu = 0: L is Cartesian
        return div @ conv, (div @ conv_bc).toarray()

    def _build_operators(self):
        exit_bc = ({"outflow": True} if self.outlet == "outflow"
                   else {"a": 1.0, "b": 0.0, "d": 0.0})
        # tube, flowing +L. Outward normal, a dw/dn + b w = d:
        #   L = 0 (inlet):  u_C = u_C0, u_H = 0 (Eq. 8)    -> a = 0, b = 1, d = (1, 0)
        #   L = 1 (outlet): nothing imposed, pure outflow -> {"outflow": True}
        #                   (or dw/dL = 0, a = 1, b = 0, as the comparison)
        self.bc_tube = ({"a": 0.0, "b": 1.0, "d": np.array([1.0, 0.0])}, exit_bc)
        if self.flow == "co":
            # shell, flowing +L: L = 0: v_H = 0 (Eq. 8)       -> a = 0, b = 1, d = 0
            #                    L = 1: pure outflow
            self.bc_shell, v_shell = ({"a": 0.0, "b": 1.0, "d": 0.0}, exit_bc), 1.0
        else:
            # shell, flowing -L: L = 0: pure outflow (outward normal -L)
            #                    L = 1: v_H = 0, clean purge -> a = 0, b = 1, d = 0
            self.bc_shell, v_shell = (exit_bc, {"a": 0.0, "b": 1.0, "d": 0.0}), -1.0
        jac_t, g_t = self._channel(self.layout.shapes["tube"], self.bc_tube, 1.0)
        jac_s, g_s = self._channel(self.layout.shapes["shell"], self.bc_shell, v_shell)
        place = self.layout.place
        self.jac_const = (place(jac_t, "tube", "tube") + place(jac_s, "shell", "shell")).tocsc()
        g = np.zeros(self.shape)
        self.layout.put(g, "tube", g_t)
        self.layout.put(g, "shell", g_s)
        self.g_const = g.reshape((-1, 1))
        self.numjac = NumJac(self.shape)        # pointwise source: last axis in full

    def scaled_source(self, w):
        """``source`` for the scaled state ``w``, divided by ``u_C0``."""
        return source(w * self.uC0, self.case, self.uC0, self.vA0) / self.uC0

    def residual(self, w):
        """Steady residual ``d(w)/dL - S(w)`` and its Jacobian."""
        g_s, jac_s = self.numjac(self.scaled_source, w)
        g = self.g_const + self.jac_const @ w.reshape((-1, 1)) - g_s.reshape((-1, 1))
        return g, self.jac_const - jac_s

    def initial_guess(self):
        """The inlet state in every cell: cyclohexane fed, no hydrogen anywhere."""
        w0 = np.zeros(self.shape)
        w0[:, 0] = 1.0
        return w0

    def solve(self, w0=None, tol=1e-10, maxfev=100):
        """Solve with ``pymrm.newton`` from ``w0`` (default ``initial_guess``).

        ``damping="monotone"`` carries the iteration over the inlet cells, where
        the rate switches from zero order to equilibrium-pinned over a narrow
        range of ``p_H``; the lower bound 0 keeps ``u_H`` and ``v_H`` inside
        the domain of Eq. 1's square roots and ``u_C`` where Eq. 7's
        denominator is positive. ``tol`` is on the scaled flows,
        above their round-off floor.

        Raises
        ------
        RuntimeError
            When Newton does not converge (``NewtonBreakdown`` propagates).
        """
        w0 = self.initial_guess() if w0 is None else np.asarray(w0, dtype=float)
        # A damping trial may evaluate the residual beyond the bounds, where the
        # square roots are NaN; newton rejects a non-finite trial and shortens the
        # step, so the warning carries no information.
        with np.errstate(invalid="ignore"):
            result = newton(self.residual, w0.reshape(self.shape), tol=tol, maxfev=maxfev,
                            damping="monotone", bounds=self.BOUNDS)
        if not result.success:
            raise RuntimeError(f"MembraneReactor ({self.flow}) at vA0 = {self.vA0:g}, "
                               f"n_z = {self.shape[0]}: {result.message}")
        self.w, self.result = result.x.reshape(self.shape), result
        return self

    def boundary_values(self):
        """Face values of ``w`` at ``L = 0`` and ``L = 1``, each ``(u_C, u_H, v_H)``."""
        tube = self.layout.view(self.w, "tube")
        shell = self.layout.view(self.w, "shell")
        t0, _, t1, _ = compute_boundary_values(tube, self.z_f, self.z_c, self.bc_tube)
        s0, _, s1, _ = compute_boundary_values(shell, self.z_f, self.z_c, self.bc_shell)
        return (np.concatenate([np.ravel(t0), np.ravel(s0)]),
                np.concatenate([np.ravel(t1), np.ravel(s1)]))

    def outlets(self, read="face"):
        """Scaled outlet flows ``(u_C, u_H, v_H)``, each at its own channel's exit.

        ``read="cell"`` takes the last cell before each exit instead of the
        face (pitfall P8); with the ``"outflow"`` outlet the two are the same.
        """
        if read == "face":
            at0, at1 = self.boundary_values()
        elif read == "cell":
            at0, at1 = self.w[0], self.w[-1]
        else:
            raise ValueError(f"read {read!r} not recognised")
        shell = at1[2] if self.flow == "co" else at0[2]
        return np.array([at1[0], at1[1], shell])

    def conversion(self, read="face"):
        """Cyclohexane conversion ``1 - u_C(L = 1)/u_C0``."""
        return 1.0 - self.outlets(read)[0]

    def h_balance_error(self, read="face"):
        """Relative hydrogen closure, the corrected Eq. 5: ``u_H + v_H = 3 (u_C0 - u_C)``.

        Outlet values of each channel at its own exit. The source conserves
        hydrogen exactly (``nu_H2 = 3``), so with face values this is a discrete
        identity in either flow arrangement.
        """
        uC, uH, vH = self.outlets(read)
        return (uH + vH - 3.0 * (1.0 - uC)) / (3.0 * (1.0 - uC))

    def face_profiles(self):
        """Every field at the faces ``z_f``, in mol/s scaled by ``u_C0``.

        First-order upwind transports the upstream cell's value through each
        face, so an interior face carries the cell before it in the flow
        direction; the boundary faces come from ``compute_boundary_values``.
        """
        at0, at1 = self.boundary_values()
        n = self.shape[0]
        out = np.empty((n + 1, 3))
        out[0], out[-1] = at0, at1
        out[1:-1, :2] = self.w[:-1, :2]                      # tube flows +L
        out[1:-1, 2] = self.w[:-1, 2] if self.flow == "co" else self.w[1:, 2]
        return out


def conversion_ivp(case, uC0, vA0, rtol=1e-10, dense=False):
    """Co-current model as an initial-value problem in ``L`` (Radau): independent route.

    The paper integrated Eqs. 2-8 by Runge-Kutta-Gill; a stiff integrator plays
    that role here, with no pymrm.

    Raises
    ------
    RuntimeError
        When the integration does not reach ``L = 1``.
    """
    def rhs(L, y):
        return source(np.asarray(y), case, uC0, vA0)
    # atol far below rtol times the flows: with a fixed atol of 1e-14 mol/s, 1e2 to
    # 1e4 times rtol * u_C0, the result moved by 3e-7 between rtol 1e-10 and 1e-12
    # Radau's implicit stages may try a negative flow, where Eq. 1's square roots are
    # NaN; the integrator rejects such a stage and shortens the step.
    with np.errstate(invalid="ignore"):
        sol = solve_ivp(rhs, [0.0, 1.0], [uC0, 0.0, 0.0], method="Radau", rtol=rtol,
                        atol=1e-3 * rtol * uC0, dense_output=dense)
    if sol.status != 0:
        raise RuntimeError(f"conversion_ivp at vA0 = {vA0:g}: {sol.message}")
    X = 1.0 - sol.y[0, -1] / uC0
    return (X, sol) if dense else X
