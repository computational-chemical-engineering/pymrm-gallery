"""Larkins, White and Jeffrey (1961): two-phase concurrent downflow in a packed bed.

Model code of gallery page G1.1. Two pieces:

- the correlation chain of the paper: Eqs. (17)-(19) (Ergun regrouped, the
  Reynolds number and the particle diameter), the ideal-gas density, and the
  two-phase correlations Eqs. (22) and (23);
- the column: ``dP/dz = rho_m(P)/144 - delta_lg(P)`` with ``z`` measured
  downward, solved three ways that share no assembly: pymrm finite volumes
  (``ColumnFV``), ``z(P)`` by adaptive quadrature inverted by ``brentq``
  (``quad_outlet``), and the paper's own sectional design procedure
  (``recipe_outlet``). A hand-written backward-Euler march
  (``backward_euler_march``) checks what the pymrm outflow boundary does.

Units are the paper's: pressure in lb-force/sq.in., length in ft, density in
lb-mass/cu.ft., mass rates in lb-mass/(sq.ft. min), viscosity in centipoise.

Contents
--------
Larkins
    The printed constants and the reading of each correlation, with Eqs. (17)
    to (19), (22), (23) and the gas density as methods.
recover_mu_air
    The air viscosity behind Table 4's Reynolds column, inverted from Eq. (18).
Column
    One operating point: the single-phase losses and the right-hand side
    ``dP/dz`` of the column.
ColumnFV
    The column on a finite-volume grid (``construct_convflux_upwind``,
    ``construct_div``, ``NumJac``, ``newton``, ``compute_boundary_values``).
quad_outlet, choking_length
    The column by quadrature of ``dz/dP``.
recipe_outlet, recipe_error
    The paper's design steps 1-5 (book p. 239) on equal sections.
backward_euler_march
    A scalar backward-Euler march whose outlet value is its last cell.
threshold_length
    The bed length at which the one-section recipe reaches a given error.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)

__all__ = ["RecipeNotConverged", "P_ATM", "R_OFF", "CP_LB_FT_MIN", "SEC_PER_MIN2", "Larkins",
           "recover_mu_air", "Column", "ColumnFV", "quad_outlet", "choking_length",
           "recipe_outlet", "recipe_error", "backward_euler_march", "threshold_length"]

# Physical constants the paper does not print.
P_ATM = 14.696                      # lb/sq.in., standard atmosphere (Table 3 is gauge)
R_OFF = 459.67                      # degF -> degR
CP_LB_FT_MIN = 6.71968975e-4 * 60.0  # 1 centipoise in lb-mass/(ft min)
SEC_PER_MIN2 = 3600.0               # g_c is printed per s^2; the notation's mu is per min


class RecipeNotConverged(RuntimeError):
    """The section-average fixed point of ``recipe_outlet`` did not converge."""


@dataclass(frozen=True)
class Larkins:
    """The printed constants of the 3/8-in. Raschig-ring bed and each correlation.

    Parameters
    ----------
    alpha, beta : float
        Eq. (17) constants (266 and 2.33 for the Raschig rings, book p. 236).
    eps : float
        Bed porosity.
    s_surf : float
        Specific surface of the packing, sq.ft./cu.ft.
    gc : float
        The printed conversion factor, 32.17.
    rho_w : float
        Liquid density, lb/cu.ft.
    rho_air_ref, t_air_ref : float
        Table 2's air density (lb/cu.ft.) and the temperature it is quoted at (degF).
    mu_air_table2 : float
        Table 2's air viscosity at 80 degF, cP.
    eq22_a, eq22_b : float
        Eq. (22): ``log10(parameter) = a / ((log10 chi)^2 + b)``.
    eq23_c : tuple of float
        Eq. (23): ``log10 R_l = c0 + c1 log10 chi + c2 (log10 chi)^2``.
    dp_override : float, optional
        Particle diameter in ft instead of Eq. (19).
    p_atm, cp : float
        Standard atmosphere and the centipoise conversion (lb/(ft min)).
    mu_air_line : tuple of float, optional
        ``(slope, intercept)`` of the recovered air viscosity in cP against degF.
    mu_air_mode : {"recovered", "table2"}
        Air viscosity from the recovered line or Table 2's single value.
    eq22_ln : bool
        Read Eq. (22) with natural logarithms (a misreading, for the break table).
    chi_inverted : bool
        Use ``sqrt(delta_g/delta_l)`` for chi (a misreading).
    gravity : bool
        Keep the static-head term of Eq. (5) (downflow); False is horizontal flow.
    two_phase : bool
        Use Eqs. (22) and (23); False gives ``delta_lg = delta_l + delta_g``.
    """

    alpha: float
    beta: float
    eps: float
    s_surf: float
    gc: float
    rho_w: float
    rho_air_ref: float
    t_air_ref: float
    mu_air_table2: float
    eq22_a: float
    eq22_b: float
    eq23_c: tuple
    dp_override: float | None = None
    p_atm: float = P_ATM
    cp: float = CP_LB_FT_MIN
    mu_air_line: tuple | None = None
    mu_air_mode: str = "recovered"
    eq22_ln: bool = False
    chi_inverted: bool = False
    gravity: bool = True
    two_phase: bool = True

    @classmethod
    def from_constants(cls, kv, **changes):
        """Build from the page's constants table (``quantity -> value``)."""
        m = cls(alpha=kv["raschig_alpha"], beta=kv["raschig_beta"],
                eps=kv["raschig_bed_porosity"], s_surf=kv["raschig_specific_surface"],
                gc=kv["g_c"], rho_w=kv["water_density"],
                rho_air_ref=kv["air_density_1atm"],
                t_air_ref=kv["air_density_reference_T"],
                mu_air_table2=kv["air_viscosity"],
                eq22_a=kv["eq22_numerator"], eq22_b=kv["eq22_offset"],
                eq23_c=(kv["eq23_c0"], kv["eq23_c1"], kv["eq23_c2"]))
        return replace(m, **changes)

    @property
    def dp(self):
        """Eq. (19): ``D_p = 6 (1 - eps) / S``, unless overridden (ft)."""
        return (self.dp_override if self.dp_override is not None
                else 6.0 * (1.0 - self.eps) / self.s_surf)

    def reynolds(self, G, mu_cP):
        """Eq. (18), ``D_p G / (mu (1 - eps))``; equal to ``6 G / (S mu)``."""
        return self.dp * G / (mu_cP * self.cp * (1.0 - self.eps))

    def delta(self, Re, mu_cP, rho):
        """Eq. (17) solved for the friction loss, lb-force/sq.in. per ft."""
        mu = mu_cP * self.cp
        return (Re * (self.alpha + self.beta * Re) * mu**2
                / (self.gc * SEC_PER_MIN2 * rho * self.dp**3)
                * ((1.0 - self.eps) / self.eps)**3 / 144.0)

    def delta_primitive(self, G, mu_cP, rho, cp=None):
        """The same loss written out in ft-lb-second units, with no Reynolds group.

        ``cp`` overrides the centipoise conversion of this route only (a defect
        for the break table: the two routes then no longer share their inputs).
        """
        mu_s = mu_cP * ((self.cp if cp is None else cp) / 60.0)   # lb/(ft s)
        u = (G / 60.0) / rho                                     # ft/s, superficial
        e, dp = self.eps, self.dp
        lam = self.alpha * mu_s * u * (1 - e)**2 / (e**3 * dp**2)
        tur = self.beta * rho * u**2 * (1 - e) / (e**3 * dp)
        return (lam + tur) / self.gc / 144.0

    def rho_air(self, P_psia, TF):
        """Ideal gas from Table 2's single entry (0.0808 lb/cu.ft. at 1 atm, 32 degF)."""
        return (self.rho_air_ref * (P_psia / self.p_atm)
                * ((self.t_air_ref + R_OFF) / (TF + R_OFF)))

    def mu_air(self, TF):
        """Air viscosity in cP: the recovered line, or Table 2's single value."""
        TF = np.asarray(TF, float)
        if self.mu_air_mode == "table2":
            return np.full_like(TF, self.mu_air_table2)
        if self.mu_air_line is None:
            raise ValueError("no recovered air-viscosity line: set mu_air_line")
        s, i = self.mu_air_line
        return s * TF + i

    def f22(self, chi):
        """Eq. (22): the two-phase parameter ``delta_lg/(delta_l + delta_g)``.

        Returns 1 in both single-phase limits (chi -> 0 and chi -> infinity).
        """
        chi = np.asarray(chi, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.log(chi) if self.eq22_ln else np.log10(chi)
            v = self.eq22_a / (t**2 + self.eq22_b)
            out = np.exp(v) if self.eq22_ln else 10.0**v
        return np.where(np.isfinite(t), out, 1.0)

    def f23(self, chi):
        """Eq. (23): the liquid saturation, declared for 0.05 < chi < 30."""
        c0, c1, c2 = self.eq23_c
        t = np.log10(np.asarray(chi, float))
        return 10.0**(c0 + c1 * t + c2 * t * t)


def recover_mu_air(m, G, TF, Re_air):
    """Invert Eq. (18) on Table 4's air Reynolds column, row by row.

    Parameters
    ----------
    m : Larkins
        Constants (``D_p``, ``eps`` and ``cp`` enter).
    G, TF, Re_air : numpy.ndarray
        Air mass rate, temperature (degF) and printed air Reynolds number of
        the rows with air.

    Returns
    -------
    slope, intercept : float
        Straight line through the recovered viscosities, cP against degF.
    resid_pct : float
        Largest relative residual of the seven values against the line, %.
    mu : numpy.ndarray
        The recovered viscosities, cP.
    """
    mu = m.dp * G / (Re_air * m.cp * (1.0 - m.eps))
    slope, icept = np.polyfit(TF, mu, 1)
    resid = float(np.max(np.abs(np.polyval((slope, icept), TF) / mu - 1.0)) * 100.0)
    return float(slope), float(icept), resid, mu


class Column:
    """Two-phase downflow at one operating point, Larkins Eqs. (5), (17)-(23).

    Parameters
    ----------
    m : Larkins
        Constants and correlation readings.
    Gl, Gg : float
        Liquid and air mass rates, lb/(sq.ft. min).
    TF : float
        Temperature, degF.
    mu_l : float
        Liquid viscosity, cP.
    """

    def __init__(self, m, Gl, Gg, TF, mu_l):
        self.m, self.Gl, self.Gg, self.TF, self.mu_l = m, Gl, Gg, TF, mu_l
        self.Rel = m.reynolds(Gl, mu_l) if Gl > 0 else 0.0
        self.mua = float(m.mu_air(TF))
        self.Reg = m.reynolds(Gg, self.mua) if Gg > 0 else 0.0
        self.dl = m.delta(self.Rel, mu_l, m.rho_w) if Gl > 0 else 0.0

    def dg(self, P):
        """Eq. (17) for the gas alone; equal to ``K/P`` through the ideal gas."""
        return (self.m.delta(self.Reg, self.mua, self.m.rho_air(P, self.TF))
                if self.Gg > 0 else 0.0)

    def state(self, P):
        """``(delta_lg, rho_m, chi, R_l)`` at pressure ``P`` (psia)."""
        m, dg = self.m, self.dg(P)
        if not m.two_phase:
            rm = m.rho_w if self.Gl > 0 else m.rho_air(P, self.TF)
            return self.dl + dg, rm, np.nan, np.nan
        if dg <= 0:
            return self.dl, m.rho_w, np.inf, 1.0
        chi = np.sqrt(dg / self.dl) if m.chi_inverted else np.sqrt(self.dl / dg)
        dlg = float(m.f22(chi)) * (self.dl + dg)
        # a saturation is a fraction: the cap acts only for misread constants
        # (Eq. (23) never exceeds 0.73 as printed)
        Rl = float(np.clip(m.f23(chi), 0.0, 1.0))
        rm = m.rho_w * Rl + m.rho_air(P, self.TF) * (1 - Rl)
        return dlg, rm, float(chi), Rl

    def dPdz(self, P):
        """Eq. (5) with Eq. (20) as a differential equation, ``z`` measured downward."""
        dlg, rm, _, _ = self.state(P)
        return (rm / 144.0 if self.m.gravity else 0.0) - dlg


def quad_outlet(col, p_in, length):
    """Outlet pressure from ``z(P) = int dP / (dP/dz)``, inverted by ``brentq``.

    ``z(P)`` is monotone wherever ``dP/dz < 0``; the sign is checked at the
    inlet. A bed at or beyond the choking length has no root and raises.
    """
    if not col.dPdz(p_in) < 0:
        raise ValueError("pressure does not fall at the inlet: z(P) is not monotone")
    f = lambda P: 1.0 / col.dPdz(P)
    g = lambda Po: quad(f, p_in, Po, limit=400)[0] - length
    return float(brentq(g, p_in * 1e-4, p_in - 1e-12, xtol=1e-12, rtol=8.9e-16))


def choking_length(col, p_in):
    """Bed length at which the correlation drives the pressure to zero, ft."""
    return float(quad(lambda P: 1.0 / col.dPdz(P), p_in, p_in * 1e-8, limit=400)[0])


def recipe_outlet(col, p_in, length, n_sections=1, rule="arithmetic", itmax=300,
                  tol=1e-14):
    """The paper's design steps 1-5 (book p. 239) on ``n_sections`` equal sections.

    Step 1: assume the section's average pressure. Steps 2-4: delta_l, delta_g,
    chi, the parameter, R_l and rho_m at that pressure. Step 5: the pressure
    drop, and confirm the assumption (a fixed-point iteration to ``tol``).

    Parameters
    ----------
    rule : {"arithmetic", "geometric", "inlet"}
        The average of the section's terminal pressures; "inlet" never forms one.

    Raises
    ------
    RecipeNotConverged
        When the fixed-point iteration of a section does not converge.
    """
    P, dz, grav = p_in, length / n_sections, col.m.gravity
    for _ in range(n_sections):
        Pb = P
        for _ in range(itmax):
            dlg, rm, _, _ = col.state(Pb)
            Pn = P + ((rm / 144.0 if grav else 0.0) - dlg) * dz
            if rule == "arithmetic":
                nb = 0.5 * (P + Pn)
            elif rule == "geometric":
                nb = float(np.sqrt(max(P * Pn, 1e-12)))
            else:
                nb = P
            done = abs(nb - Pb) < tol * max(1.0, abs(Pb))
            Pb = nb
            if done:
                break
        else:
            raise RecipeNotConverged(f"section average did not converge in {itmax} iterations")
        dlg, rm, _, _ = col.state(Pb)
        P = P + ((rm / 144.0 if grav else 0.0) - dlg) * dz
    return float(P)


def recipe_error(col, p_in, length, n_sections=1, rule="arithmetic", tol=1e-14):
    """Relative error of the recipe in the pressure drop, against the quadrature."""
    Pe = quad_outlet(col, p_in, length)
    Pr = recipe_outlet(col, p_in, length, n_sections, rule, tol=tol)
    return abs((Pr - Pe) / (p_in - Pe))


def backward_euler_march(col, p_in, length, n, itmax=200, tol=1e-14):
    """``P_i = P_(i-1) + h dP/dz(P_i)`` cell by cell, each step by scalar Newton.

    Shares no code with ``ColumnFV``: an independent statement of the scheme
    the pymrm outflow boundary is meant to give.

    Raises
    ------
    RuntimeError
        When a step does not converge.
    """
    h, P = length / n, p_in
    for _ in range(n):
        x = P
        for _ in range(itmax):
            f = x - P - h * col.dPdz(x)
            xd = x * (1.0 + 1e-7)
            df = ((xd - P - h * col.dPdz(xd)) - f) / (xd - x)
            xn = x - f / df
            done = abs(xn - x) <= tol * max(1.0, abs(x))
            x = xn
            if done:
                break
        else:
            raise RuntimeError("backward-Euler step did not converge")
        P = x
    return float(P)


class ColumnFV:
    """The column ``dP/dz = rho_m/144 - delta_lg`` on a finite-volume grid.

    First-order upwind convection at ``v = 1`` (``construct_convflux_upwind``)
    and ``construct_div`` with ``nu = 0``; the pointwise source through
    ``NumJac((n, 1))``, so its Jacobian is diagonal. Operators are assembled
    once, in ``_build_operators``.

    Parameters
    ----------
    col : Column
        The operating point.
    p_in : float
        Inlet pressure, psia.
    length : float
        Bed length, ft.
    n : int
        Number of cells.
    outlet : {"outflow", "zero_gradient"}
        The outlet boundary. ``"outflow"`` (pymrm ``{"outflow": True}``) puts no
        condition there, which is right for a first-order problem: the face
        value is the last cell. ``"zero_gradient"`` imposes ``dP/dz = 0`` at
        ``z = L``, which the problem does not have; kept for comparison.
    """

    def __init__(self, col, p_in, length, n, outlet="outflow"):
        self.col, self.p_in, self.length, self.n, self.outlet = col, p_in, length, n, outlet
        self._build_grid()
        self._build_operators()

    def _build_grid(self):
        self.x_f = np.linspace(0.0, self.length, self.n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.shape = (self.n, 1)
        # bc on the OUTWARD normal, a dP/dn + b P = d
        #   inlet  z = 0: P = p_in                -> a = 0, b = 1, d = p_in
        #   outlet z = L: no condition (outflow)  -> {"outflow": True}
        #         or, for comparison, dP/dz = 0   -> a = 1, b = 0, d = 0
        right = ({"outflow": True} if self.outlet == "outflow"
                 else {"a": 1.0, "b": 0.0, "d": 0.0})
        self.bc = ({"a": 0.0, "b": 1.0, "d": self.p_in}, right)

    def _build_operators(self):
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                  self.bc, v=1.0)
        div = construct_div(self.shape, self.x_f, nu=0)     # nu = 0: Cartesian axial
        self.conv = conv
        self.jac_const = (div @ conv).tocsr()
        self.g_const = np.asarray((div @ conv_bc).todense()).reshape(-1, 1)
        self.numjac = NumJac(self.shape)                    # (n, 1): diagonal source
        self._src = np.vectorize(self.col.dPdz)

    def residual(self, P):
        """``d(v P)/dz - dP/dz(P)`` per cell, and its Jacobian."""
        g_src, jac_src = self.numjac(lambda p: -self._src(p), P.reshape(self.shape))
        g = self.jac_const @ P.reshape(-1, 1) + self.g_const + g_src.reshape(-1, 1)
        return g, self.jac_const + jac_src

    def solve(self, tol=1e-10):
        """Cell pressures, shape ``(n,)``; raises unless ``newton`` converged.

        ``tol`` is an absolute step on pressures of order 100 psia, above the
        round-off floor of the discrete solution; ``bounds`` keeps the
        iterates where the gas density and chi are defined (``P > 0``).
        """
        res = newton(self.residual, np.full(self.shape, self.p_in), tol=tol,
                     maxfev=50, bounds=(0.0, None))
        if not res.success:
            raise RuntimeError(f"newton did not converge: {res.message}")
        self.result = res
        return np.asarray(res.x).reshape(-1)

    def outlet_value(self, P, read="face"):
        """Outlet pressure: the face value the operator transports, or the last cell.

        ``read="cell"`` is pitfall P8 unless the outlet is ``"outflow"``.
        """
        if read == "cell":
            return float(P[-1])
        v = compute_boundary_values(np.reshape(P, self.shape), self.x_f, self.x_c,
                                    self.bc[1], axis=0, bound_id=1)[0]
        return float(np.asarray(v).reshape(-1)[0])

    def balance_gap(self, P, read="face"):
        """Integral balance ``P_out - P_in - sum h dP/dz(P_i)``, over ``P_in - P_out``.

        Zero to round-off for the transported face value whatever the source
        is: it checks the discretisation's bookkeeping, not the physics.
        """
        p_out = self.outlet_value(P, read)
        h = np.diff(self.x_f)
        integral = float(np.sum(h * self._src(P)))
        return (p_out - self.p_in - integral) / (self.p_in - p_out)


def threshold_length(col, p_in, target=0.01, rule="arithmetic", lo=0.05, n_check=12,
                     xtol=1e-10):
    """Bed length at which the one-section recipe's error reaches ``target``.

    The bracket starts at ``lo``, halved until the error is below target (the
    recipe is second order, so the error vanishes as the bed shortens), and its
    upper end is doubled from there until the error exceeds target; near the
    choking length the step halves the remaining distance instead, so the
    bracket never reaches it, and where the one-section fixed point does not
    converge (it can diverge close to choking) the step is halved back toward
    the last good point. Every evaluation is a converged recipe and a converged
    quadrature. The error is then checked to rise monotonically on ``n_check``
    points of the bracket (a check on those points, not a proof that the root
    is unique).

    Returns
    -------
    length : float
        The root, ft.
    bracket : tuple of float
        The bracket it was found on.

    Raises
    ------
    ValueError
        When the error does not reach ``target`` before the bed chokes, or is
        not monotone on the bracket.
    """
    Lch = choking_length(col, p_in)
    err = lambda L: recipe_error(col, p_in, L, rule=rule)
    while err(lo) >= target:
        lo *= 0.5
        if lo < 1e-6:
            raise ValueError("the recipe error does not vanish for short beds")
    good = lo
    while True:
        # double, or halve the remaining distance to choking once doubling would pass it
        hi = 2.0 * good if 2.0 * good < Lch else good + 0.5 * (Lch - good)
        if hi > (1.0 - 1e-6) * Lch:
            raise ValueError(f"the recipe error stays below {target} up to choking")
        for _ in range(40):
            try:
                e_hi = err(hi)
                break
            except RecipeNotConverged:
                # near choking the one-section fixed point can diverge: step back
                hi = 0.5 * (good + hi)
        else:
            raise RuntimeError("no converged recipe beyond the last bracket point")
        if e_hi >= target:
            break
        good = hi
    e = np.array([err(L) for L in np.geomspace(lo, hi, n_check)])
    if not np.all(np.diff(e) > 0):
        raise ValueError("the recipe error is not monotone on the bracket")
    return float(brentq(lambda L: err(L) - target, lo, hi, xtol=xtol, rtol=8.9e-16)), (lo, hi)
