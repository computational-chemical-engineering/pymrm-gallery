"""Model of page E1.1: a compressible gas column through a bubbling fluidised bed.

Toomey and Johnstone (1952) note that the gas must expand as the pressure falls
through the bed and never compute the effect. ``CompressibleBed`` solves for the
pressure profile and the expanded bed height with their own equation (1)
expansion law and their own measured bed weight, by a pymrm finite-volume Newton
solve, and gives two validation routes: the closed form (the equation
separates) and an adaptive Runge-Kutta integration with a terminal event.
"""

from __future__ import annotations

import numpy as np
from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, generate_grid, newton)
from pymrm.checks import find_roots
from scipy.integrate import quad, solve_ivp

PATM = 101325.0                     # Pa
M_AIR, R_GAS = 0.0289647, 8.314462618  # kg/mol, J/(mol K)


def mu_air(T):
    """Viscosity of air by Sutherland's law.

    Parameters
    ----------
    T : float
        Temperature, K.

    Returns
    -------
    float
        Dynamic viscosity, Pa s.
    """
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def rho_air(T, P=PATM):
    """Density of air as an ideal gas.

    Parameters
    ----------
    T : float
        Temperature, K.
    P : float, optional
        Pressure, Pa.

    Returns
    -------
    float
        Density, kg/m3.
    """
    return P * M_AIR / (R_GAS * T)


def unique_root(f, lo, hi, n_scan=8, xtol=1e-13):
    """The one root of ``f`` on ``[lo, hi]``, by a ``find_roots`` scan.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        Physically possible range.
    n_scan : int, optional
        Number of scan points.
    xtol : float, optional
        Absolute tolerance of the refined root.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        Unless the scan finds exactly one sign change and no failed point.
    """
    r = find_roots(f, lo, hi, n_scan=n_scan, xtol=xtol)
    if len(r["roots"]) != 1 or r["failed"]:
        raise ValueError(f"expected one root on [{lo:.6g}, {hi:.6g}]: {r['message']}")
    return float(r["roots"][0])


class CompressibleBed:
    """Isothermal ideal gas rising through a bubbling bed of fixed solids inventory.

    The depth ``s`` is measured downward from the bed surface, on the
    dimensionless coordinate ``zeta = s / L_f``. The gas mass flux is fixed, so
    ``u(P) = u_top P_top / P``, and the pressure gradient is the bed weight
    reduced by the bubble fraction,

        dP/ds = (1 - delta(P)) w_mf,   delta(P) = max(0, a (u(P) - u_e)),

    with ``P(0) = P_top``. The bed height ``L_f`` closes the solids inventory,
    ``int_0^{L_f} (1 - delta) ds = L_mf``.

    Parameters
    ----------
    w_mf : float
        Bed weight per unit volume at minimum fluidisation, dP_mf/L_mf, Pa/m.
    a : float
        Equation (1) coefficient, 0.065/D_p^0.5 in SI units, s/m.
    u_e : float
        Dense-phase superficial velocity, m/s.
    u_top : float
        Superficial gas velocity at the bed surface, m/s.
    L_mf : float
        Bed height at minimum fluidisation, m.
    P_top : float, optional
        Pressure at the bed surface, Pa.
    n : int, optional
        Number of finite-volume cells.
    v : float, optional
        Upwind direction of the pressure transport on ``zeta``; ``+1`` carries
        the surface value downward (a break row reverses it).
    P_inlet : float, optional
        Pressure imposed at ``zeta = 0``; defaults to ``P_top`` (a break row
        imposes another value).
    outlet_bc : dict, optional
        Boundary condition at the distributor, ``{"outflow": True}`` by default
        (a break row passes the zero-gradient dictionary).
    """

    def __init__(self, w_mf, a, u_e, u_top, L_mf, P_top=PATM, n=400, v=1.0, P_inlet=None,
                 outlet_bc=None):
        self.w, self.a, self.u_e, self.u_top = w_mf, a, u_e, u_top
        self.L_mf, self.P_top, self.n, self.v = L_mf, P_top, n, v
        self.P_inlet = P_top if P_inlet is None else P_inlet
        self.outlet_bc = {"outflow": True} if outlet_bc is None else outlet_bc
        self._build_grid()
        self._build_operators()

    def _build_grid(self):
        self.x_f, self.x_c = generate_grid(self.n, [0.0, 1.0], generate_x_c=True)

    def _build_operators(self):
        # Unknown p = P / P_top, of order one. Bc dictionaries use the OUTWARD
        # normal. Left (zeta = 0, bed surface): p = P_inlet / P_top (a=0, b=1).
        # Right (zeta = 1, distributor): pure outflow, the face carries the
        # adjacent cell value; with v > 0 nothing enters there.
        self.bc = ({"a": 0, "b": 1, "d": self.P_inlet / self.P_top}, self.outlet_bc)
        conv, conv_bc = construct_convflux_upwind((self.n, 1), self.x_f, self.x_c,
                                                  self.bc, v=self.v)
        div = construct_div((self.n, 1), self.x_f, nu=0)   # nu=0: plane column
        self.jac_const, self.g_const = div @ conv, div @ conv_bc
        self.numjac = NumJac((self.n, 1))                  # one field: (n, 1), never (n,)

    def delta(self, P):
        """Bubble fraction from equation (1) at pressure ``P`` (Pa)."""
        return np.maximum(0.0, self.a * (self.u_top * self.P_top / P - self.u_e))

    def source(self, p, Lf):
        """Weight term of d(v p)/dzeta = L_f w (1 - delta) / P_top."""
        return self.w * Lf / self.P_top * (1.0 - self.delta(p * self.P_top))

    def residual(self, p, Lf):
        """Residual and Jacobian of the steady pressure equation."""
        p = p.reshape(-1, 1)
        g_src, jac_src = self.numjac(lambda q: self.source(q, Lf), p)
        g = self.g_const + self.jac_const @ p - g_src.reshape(-1, 1)
        return g, self.jac_const - jac_src

    def solve(self, Lf):
        """Pressure profile for a given bed height.

        Parameters
        ----------
        Lf : float
            Expanded bed height, m.

        Returns
        -------
        dict
            ``P`` (cell values, Pa), ``P_bottom`` (distributor face, Pa),
            ``delta`` (cell values), ``nit`` and ``inventory_residual``
            (``L_f (1 - mean delta) - L_mf``, m).

        Raises
        ------
        RuntimeError
            When Newton does not converge.
        """
        p0 = np.full((self.n, 1), 1.0 + 0.5 * self.w * Lf / self.P_top)
        sol = newton(lambda p: self.residual(p, Lf), p0, tol=1e-12, maxfev=60)
        if not sol.success:
            raise RuntimeError(f"newton failed at L_f = {Lf:.6g} m: {sol.message}")
        p = sol.x.reshape(-1, 1)
        p_top, _, p_bot, _ = compute_boundary_values(p, self.x_f, self.x_c, self.bc)
        P = p.ravel() * self.P_top
        d = self.delta(P)
        return dict(P=P, P_surface=float(np.ravel(p_top)[0]) * self.P_top,
                    P_bottom=float(np.ravel(p_bot)[0]) * self.P_top, delta=d,
                    nit=int(sol.nit), inventory_residual=Lf * (1.0 - float(d.mean())) - self.L_mf)

    def excess_spread(self, sol, u_mf, at="faces"):
        """Spread of the excess velocity ``u - u_mf`` across the bed.

        ``(u(surface) - u(distributor)) / depth average of (u - u_mf)``.

        Parameters
        ----------
        sol : dict
            Result of :meth:`solve`.
        u_mf : float
            Minimum-fluidisation velocity, m/s.
        at : {"faces", "cells"}, optional
            ``"faces"``: the end values at the surface and distributor faces
            (``compute_boundary_values``) and the average by the trapezoidal rule
            over the face values the upwind scheme transports (face ``i + 1/2``
            carries cell ``i``). ``"cells"``: the first and last cell values and
            the cell mean, which read the surface half a cell to one cell short.

        Returns
        -------
        float
            The spread.
        """
        u = lambda P: self.u_top * self.P_top / np.asarray(P)
        if at == "cells":
            uc = u(sol["P"])
            return float((uc[0] - uc[-1]) / (uc.mean() - u_mf))
        uf = u(np.concatenate([[sol["P_surface"]], sol["P"][:-1], [sol["P_bottom"]]]))
        return float((uf[0] - uf[-1]) / (np.mean(0.5 * (uf[1:] + uf[:-1])) - u_mf))

    def Lf(self, xtol=1e-13):
        """Expanded bed height (m), the one root of the solids-inventory residual.

        The physical range is ``L_mf <= L_f <= L_mf / (1 - delta_max)``, with
        ``delta_max`` the bubble fraction at the surface (zeta = 0) pressure (delta falls as
        the pressure rises); the scan runs slightly beyond both ends. ``xtol``
        is the absolute tolerance of the root, m.
        """
        d_max = float(self.delta(self.P_inlet))
        if d_max >= 1.0:
            raise ValueError("equation (1) gives a bed of pure bubbles at the inlet")
        lo, hi = self.L_mf * (1 - 1e-3), self.L_mf / (1.0 - d_max) * (1 + 1e-3)
        return unique_root(lambda L: self.solve(L)["inventory_residual"], lo, hi, xtol=xtol)

    # -- validation route 1: the closed form ---------------------------------
    # Where the bed bubbles, dP/ds = w (b - c/P) with b = 1 + a u_e and
    # c = a u_top P_top, which integrates to [P + (c/b) ln(bP - c)]/b = w s + const.
    def _F(self, P):
        b = 1.0 + self.a * self.u_e
        c = self.a * self.u_top * self.P_top
        return (P + (c / b) * np.log(b * P - c)) / b

    def P_switch(self):
        """Pressure at which ``u(P) = u_e``: above it the bed does not bubble."""
        return self.u_top * self.P_top / self.u_e

    def s_analytic(self, P):
        """Depth (m) at which the pressure reaches ``P``, closed form."""
        Ps = self.P_switch()
        if self.P_top >= Ps:                          # never bubbles: linear throughout
            return (P - self.P_top) / self.w
        if P <= Ps:                                   # bubbling branch
            return (self._F(P) - self._F(self.P_top)) / self.w
        return (self._F(Ps) - self._F(self.P_top)) / self.w + (P - Ps) / self.w

    def P_bottom_exact(self):
        """Distributor pressure: integrating dP/ds over the bed gives
        ``P(L_f) - P_top = w int (1 - delta) ds = w L_mf`` exactly."""
        return self.P_top + self.w * self.L_mf

    def Lf_analytic(self):
        """Expanded bed height from the closed form."""
        return self.s_analytic(self.P_bottom_exact())

    def mean_excess_velocity(self, u_mf):
        """Depth average of ``u - u_mf`` over the bed, by quadrature in ``P``.

        ``int (u - u_mf) ds = int (u(P) - u_mf) / (w (1 - delta(P))) dP``.
        """
        Pb, Ps = self.P_bottom_exact(), self.P_switch()
        pts = [Ps] if self.P_top < Ps < Pb else None
        f = lambda P: (self.u_top * self.P_top / P - u_mf) / (self.w * (1.0 - self.delta(P)))
        val, _ = quad(f, self.P_top, Pb, points=pts, epsabs=0.0, epsrel=1e-13, limit=200)
        return val / self.Lf_analytic()

    # -- validation route 2: adaptive Runge-Kutta ----------------------------
    def Lf_ivp(self, rtol=1e-12):
        """Expanded bed height: integrate dP/ds from the surface until P reaches
        ``P_bottom_exact`` (a terminal event), with relative tolerance ``rtol``."""
        Pb = self.P_bottom_exact()
        ev = lambda s, P: P[0] - Pb
        ev.terminal, ev.direction = True, 1
        r = solve_ivp(lambda s, P: self.w * (1.0 - self.delta(np.asarray(P))),
                      [0.0, 60 * self.L_mf], [self.P_top],
                      rtol=rtol, atol=1e-8, events=ev)
        return float(r.t_events[0][0])
