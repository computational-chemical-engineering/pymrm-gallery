"""Voltz et al. (1973) CO and propylene oxidation kinetics and their integral reactor.

Units follow the paper (its eqs 2 to 6 and Table I): CO and O2 in mol %, the
hydrocarbon and NO in ppm, rates per unit catalyst volume in those units per
second, temperatures in degrees Fahrenheit with ``T + 460`` the absolute
(Rankine) temperature of the Arrhenius terms. Every published constant sits in
one dictionary ``par`` (keys ``k0_r1 ... k0_a4`` and ``E_r1 ... E_a4``, the
latter as E/R_g in degR), which every function takes.

No plotting and no data loading here; the page loads Table I and passes it in.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)
from pymrm.checks import find_roots

RANKINE_OFFSET = 460.0      # the paper's eqs 5 and 6: T_s + 460
O2_PER_CO = 0.5             # CO + 1/2 O2 -> CO2
O2_PER_C3H6 = 4.5           # C3H6 + 9/2 O2 -> 3 CO2 + 3 H2O


def rate_constants(par, T_F):
    """Arrhenius values of the six Table I parameters at ``T_F`` (degF).

    Returns
    -------
    dict
        ``kr1, kr2`` (1/(s mol% O2)) and ``ka1`` (1/mol% CO), ``ka2`` (1/ppm),
        ``ka3`` (1/(mol% ppm)^2), ``ka4`` (1/ppm^0.7).
    """
    T_R = T_F + RANKINE_OFFSET
    out = {}
    for name in ("r1", "r2", "a1", "a2", "a3", "a4"):
        out["k" + name] = par["k0_" + name] * np.exp(-par["E_" + name] / T_R)
    return out


def resistance(k, co, hc, no):
    """The inhibition term R(theta) of eq 4 (``hc`` in the rate's ppm basis)."""
    return ((1.0 + k["ka1"] * co + k["ka2"] * hc) ** 2
            * (1.0 + k["ka3"] * (co * hc) ** 2)
            * (1.0 + k["ka4"] * no ** 0.7))


def rates(k, co, hc, o2, no):
    """Disappearance rates of eqs 2 and 3: ``(r_CO in mol%/s, r_HC in ppm/s)``."""
    R = resistance(k, co, hc, no)
    return k["kr1"] * co * o2 / R, k["kr2"] * hc * o2 / R


def co_rate_maximum(par, T_F, hc, no=0.0, power=2, co_range=(1e-4, 50.0)):
    """CO concentration (mol %) at which r_CO of eq 2 is largest, at fixed HC, NO and O2.

    The maximum is the root of the analytic derivative of ln r_CO,
    ``1/CO - power*ka1/(1 + ka1 CO + ka2 HC) - 2 ka3 CO HC^2/(1 + ka3 CO^2 HC^2)``,
    found with :func:`pymrm.checks.find_roots` over ``co_range``. ``power`` is the
    exponent of the adsorption group (2 in the paper). Returns ``nan`` unless
    exactly one root lies in the range. O2 and the NO term are factors that do
    not move the maximum.
    """
    k = rate_constants(par, T_F)

    def dlnr(co):
        return (1.0 / co - power * k["ka1"] / (1.0 + k["ka1"] * co + k["ka2"] * hc)
                - 2.0 * k["ka3"] * co * hc ** 2 / (1.0 + k["ka3"] * (co * hc) ** 2))

    roots = find_roots(dlnr, *co_range, n_scan=200, log=True)["roots"]
    return roots[0] if len(roots) == 1 else np.nan


def co_rate_maximum_closed_form(par, T_F, hc):
    """The same maximum with the ka3 term dropped: ``CO* = (1 + ka2 HC) / ka1``."""
    k = rate_constants(par, T_F)
    return (1.0 + k["ka2"] * hc) / k["ka1"]


class IntegralReactor:
    """Isothermal integral (plug-flow) reactor of eq 1 on the fractional bed length.

    The unknowns are CO, HC and O2 scaled by their inlet values, ``u`` of shape
    ``(n, 3)`` on ``x`` in [0, 1]. Steady convection at unit velocity with the
    rates of eqs 2 and 3 multiplied by the space time ``tau = V/F`` (s):
    ``du/dx = -tau r(u) / Y_in``.

    Parameters
    ----------
    par : dict
        Table I.
    T_F : float
        Catalyst temperature (degF).
    inlet : sequence of 4 floats
        ``(CO mol %, HC ppm as propylene, O2 mol %, NO ppm)``, the legend of
        Figures 19 and 20.
    hc_basis : float
        Divisor taking propylene ppm to the ppm basis of the rate equations:
        2 for ppm as C6 (the paper's "(C3H6) = ppm of propylene (as C6 FID)"),
        1 for ppm of propylene.
    o2_consumed : bool
        Deplete O2 by the stoichiometry of both oxidations.
    n : int
        Number of cells.
    """

    def __init__(self, par, T_F, inlet, hc_basis=2.0, o2_consumed=True, n=400):
        co, hc_ppm, o2, no = (float(v) for v in inlet)
        self.k = rate_constants(par, T_F)
        self.no = no
        self.hc_basis = float(hc_basis)
        self.o2_factor = 1.0 if o2_consumed else 0.0
        self.y_in = np.array([co, hc_ppm / self.hc_basis, o2])
        self.shape = (int(n), 3)
        self.x_f = np.linspace(0.0, 1.0, self.shape[0] + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self._build_operators()

    def _build_operators(self):
        # x = 0, outward normal -x: u = 1 for every field (inlet)  -> a = 0, b = 1, d = 1
        # x = 1, outward normal +x: convective outflow, the face carries the cell value
        self.bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"outflow": True})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                  self.bc, v=1.0, axis=0)
        div = construct_div(self.shape, self.x_f, nu=0, axis=0)   # nu=0: Cartesian bed
        self.jac_const = div @ conv
        self.g_const = (div @ conv_bc).toarray().reshape(-1, 1)
        self.numjac = NumJac(self.shape)      # pointwise source: the 3 fields coupled

    def source(self, u):
        """``-dY/dtau / Y_in`` per cell for scaled unknowns ``u`` of shape ``(..., 3)``."""
        y = u * self.y_in
        r_co, r_hc = rates(self.k, y[..., 0], y[..., 1], y[..., 2], self.no)
        r_o2 = self.o2_factor * (O2_PER_CO * r_co
                                 + O2_PER_C3H6 * self.hc_basis * r_hc * 1e-4)
        return np.stack([r_co, r_hc, r_o2], axis=-1) / self.y_in

    def residual(self, u, tau):
        g_src, jac_src = self.numjac(lambda v: -tau * self.source(v), u)
        g = self.g_const + self.jac_const @ u.reshape(-1, 1) - g_src.reshape(-1, 1)
        return g, self.jac_const - jac_src

    def solve(self, tau, u0=None):
        """Solve at space time ``tau`` (s); returns outlet conversions (%) of CO and HC."""
        if u0 is None:
            u0 = np.ones(self.shape)
        # A trial step far above the inlet value can overflow R(theta); raising
        # makes newton treat it as a recoverable step and shorten it.
        with np.errstate(over="raise"):
            result = newton(lambda u: self.residual(u, tau), u0, tol=1e-12,
                            bounds=(0.0, None), stagnation=True, maxfev=60)
        if not result.success:
            raise RuntimeError(f"IntegralReactor.solve: {result.message}")
        self.u = result.x.reshape(self.shape)
        self.nit = result.nit
        _, _, u_out, _ = compute_boundary_values(self.u, self.x_f, self.x_c, self.bc)
        u_out = np.asarray(u_out).reshape(-1)
        return 100.0 * (1.0 - u_out[0]), 100.0 * (1.0 - u_out[1])


def march(par, T_F, inlet, tau, hc_basis=2.0, o2_consumed=True, rtol=1e-11):
    """Independent route: integrate eq 1 in tau with ``solve_ivp`` (LSODA).

    Returns outlet conversions (%) of CO and HC at ``tau`` (s).
    """
    rx = IntegralReactor(par, T_F, inlet, hc_basis=hc_basis,
                         o2_consumed=o2_consumed, n=1)

    def rhs(t, u):
        return -rx.source(u)

    sol = solve_ivp(rhs, (0.0, tau), np.ones(3), method="LSODA", rtol=rtol,
                    atol=1e-14)
    if sol.status != 0:
        raise RuntimeError(f"march: {sol.message}")
    u = sol.y[:, -1]
    return 100.0 * (1.0 - u[0]), 100.0 * (1.0 - u[1])


def tau_at_conversion(conv_of_tau, level, lo=1e-5, hi=10.0):
    """Space time (s) at which ``conv_of_tau(tau)`` equals ``level`` (%).

    Root of ``conv - level`` by :func:`pymrm.checks.find_roots` on a log scan of
    ``[lo, hi]``; ``nan`` unless exactly one root.
    """
    roots = find_roots(lambda t: conv_of_tau(t) - level, lo, hi, n_scan=16,
                       log=True, xtol=1e-14)["roots"]
    return roots[0] if len(roots) == 1 else np.nan
