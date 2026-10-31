"""Oh and Cavendish (1982): transient monolith channel, gas and wall in one Newton system.

The model of their eqs 1 to 18 for one channel of a monolith catalytic converter:
a quasi-steady gas phase (temperature and five mole fractions, first order in x)
convected past a wall that stores heat, conducts it axially and carries the
reaction at its surface (temperature and five surface mole fractions). The twelve
fields share one axial grid, so the state is ``(n_x, 12)`` with the field axis
last; the conduction, convection and divergence operators are constant and
assembled once, and ``NumJac`` on the last axis supplies the pointwise Jacobian of
reaction, film transfer and accumulation. Each backward-Euler time step is one
``newton`` solve.

No plotting and no data loading: the parameter table is read by the notebook
and handed to :meth:`OhCavendish.from_table` and :meth:`Feed.from_table` as a
lookup function.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse import diags_array

from pymrm import (NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, construct_grad,
                   generate_grid, newton)

R_GAS = 8.314          # J/(mol K)
M_GAS = 28.97          # g/mol, exhaust taken as air-like (not in the paper)
SPECIES = ("CO", "C3H6", "CH4", "H2", "O2")

# Heats of combustion, J/mol, from standard formation enthalpies with water as
# vapour. Not printed in the paper; standard thermochemistry.
DH = np.array([282.98e3, 1926.42e3, 802.29e3, 241.826e3])
# water formed per mole of CO, C3H6, CH4, H2 burnt, and the standard formation
# enthalpies of water vapour and liquid water, kJ/mol
N_H2O = np.array([0, 3, 2, 1])
HF_H2O = {"gas": -241.826, "liquid": -285.83}


def heats_of_combustion(water="gas"):
    """Heats of combustion ``(-dH)_i`` of CO, C3H6, CH4 and H2, J/mol.

    Parameters
    ----------
    water : {"gas", "liquid"}, optional
        State of the product water. The exhaust is far above its dew point,
        so ``"gas"`` (:data:`DH`) is the model's value; ``"liquid"`` adds the
        heat of condensation of the water formed.

    Returns
    -------
    numpy.ndarray
        Four heats of combustion, positive for an exothermic reaction.
    """
    return DH + 1e3 * N_H2O * (HF_H2O["gas"] - HF_H2O[water])

# Critical constants for the Slattery-Bird formula (Bird et al. 1960, p. 505).
TC = {"CO": 133.0, "C3H6": 365.0, "CH4": 190.7, "H2": 33.3, "O2": 154.4, "N2": 126.2}
PC = {"CO": 34.5, "C3H6": 45.6, "CH4": 45.8, "H2": 12.8, "O2": 49.7, "N2": 33.5}
MW = {"CO": 28.01, "C3H6": 42.08, "CH4": 16.04, "H2": 2.016, "O2": 32.00, "N2": 28.013}
SB_A, SB_B = 2.745e-4, 1.823          # non-polar pair, Bird et al. eq. 16.3-1

# field index in the last axis of the state
I_TS, I_TG = 0, 1
I_YG = slice(2, 7)
I_YS = slice(7, 12)
N_FIELD = 12


class StepFailed(RuntimeError):
    """A time step whose Newton solve did not converge; ``result`` holds it."""

    def __init__(self, message, result=None):
        super().__init__(message)
        self.result = result


def sb_coefficient(sp, b_sp="N2"):
    """Temperature-independent factor of the Slattery-Bird diffusivity.

    Parameters
    ----------
    sp, b_sp : str
        The diffusing species and the bath gas.

    Returns
    -------
    float
        ``D / T**SB_B`` in cm2/(s K**SB_B) at 1 atm.
    """
    tc = TC[sp] * TC[b_sp]
    return (SB_A * tc ** (-SB_B / 2) * (PC[sp] * PC[b_sp]) ** (1 / 3)
            * tc ** (5 / 12) * np.sqrt(1 / MW[sp] + 1 / MW[b_sp]))


SB_COEF = np.array([sb_coefficient(s) for s in SPECIES])


def d_slattery_bird(sp, T, p_atm=1.0):
    """Binary diffusivity of ``sp`` in N2, cm2/s (Bird et al. 1960, eq. 16.3-1).

    Parameters
    ----------
    sp : str
        Species name, one of :data:`SPECIES`.
    T : float or numpy.ndarray
        Temperature, K.
    p_atm : float, optional
        Pressure, atm.

    Returns
    -------
    float or numpy.ndarray
        Diffusivity, cm2/s.
    """
    return sb_coefficient(sp) * T ** SB_B / p_atm


@dataclass(frozen=True)
class OhCavendish:
    """Table I and the rate constants of eqs 5 and 6 (both cases).

    Attributes
    ----------
    p_tot : float
        Total pressure, Pa.
    w_g, area, length : float
        Exhaust mass flow (g/s), frontal area (cm2), monolith length (cm).
    r_h : float
        Hydraulic radius as the paper defines it, 2 * area / wetted perimeter
        (half the usual hydraulic diameter), cm.
    eps, rho_s, a_pt : float
        Void fraction, solid density (g/cm3), platinum area per reactor
        volume (cm2 Pt/cm3).
    lam_s, cp_g : float
        Solid conductivity J/(cm s K) and gas heat capacity J/(g K).
    cps : tuple of float
        ``C_ps = c0 + c1 T + c2 / T**2``, J/(g K).
    kin : tuple
        ``(A, E/R)`` of k1, k2, k3 (eq 5).
    ads : tuple
        ``(A, dH/R)`` of K1 to K4, ``K = A exp(-(dH/R)/T)`` (eq 6).
    lamg : tuple
        ``(A, n)`` of ``lambda_g = A T**n``, J/(cm s K).
    """

    p_tot: float
    w_g: float
    area: float
    length: float
    r_h: float
    eps: float
    rho_s: float
    a_pt: float
    lam_s: float
    cp_g: float
    cps: tuple
    kin: tuple
    ads: tuple
    lamg: tuple

    @classmethod
    def from_table(cls, value):
        """Build from a lookup ``value(symbol) -> float`` over the Table I symbols."""
        return cls(
            p_tot=value("P_tot") * 1e3, w_g=value("W_g"), area=value("A"),
            length=value("L"), r_h=value("R_h"), eps=value("eps"),
            rho_s=value("rho_s"), a_pt=value("a"), lam_s=value("lambda_s"),
            cp_g=value("C_pg"),
            cps=(value("C_ps_c0"), value("C_ps_c1"), value("C_ps_c2")),
            kin=tuple((value(f"k{i}_A"), value(f"k{i}_E_over_R")) for i in (1, 2, 3)),
            ads=tuple((value(f"K{i}_A"), value(f"K{i}_dH_over_R")) for i in (1, 2, 3, 4)),
            lamg=(value("lambda_g_A"), value("lambda_g_n")))

    def lam_gas(self, T):
        """Gas thermal conductivity, approximated by N2 (paper, p. 31), J/(cm s K)."""
        return self.lamg[0] * T ** self.lamg[1]

    def cp_solid(self, T):
        """Table I: solid heat capacity ``C_ps(T)``, J/(g K)."""
        return self.cps[0] + self.cps[1] * T + self.cps[2] / T ** 2

    @property
    def s_geo(self):
        """Geometric surface per reactor volume ``S = 2 eps / R_h``, cm2/cm3."""
        return 2.0 * self.eps / self.r_h

    @property
    def g_mass(self):
        """Mass flux ``W_g / A`` through the frontal area, g/(cm2 s): constant."""
        return self.w_g / self.area


@dataclass(frozen=True)
class Feed:
    """One case of Table I: feed after the step and the initial wall temperature.

    Attributes
    ----------
    y_in : tuple of float
        Inlet mole fractions of CO, C3H6, CH4, H2, O2.
    y_no : float
        Inlet NO mole fraction (enters only the inhibition term of eq 6).
    T_g_in, T_s0 : float
        Inlet gas temperature after the step and the uniform initial wall
        temperature, K.
    """

    y_in: tuple
    y_no: float
    T_g_in: float
    T_s0: float

    @classmethod
    def from_table(cls, value, case):
        """Build from a lookup ``value(symbol, case) -> float``."""
        return cls(y_in=tuple(value(f"y_{s}_in", case) for s in SPECIES),
                   y_no=value("y_NO_in", case), T_g_in=value("T_g_in", case),
                   T_s0=value("T_s0", case))


@dataclass
class March:
    """What :meth:`Monolith.march` returns.

    Attributes
    ----------
    t : numpy.ndarray
        Time levels, s, starting at 0.
    records : list
        ``record(model)`` at every level.
    nit : numpy.ndarray
        Newton iterations of every step.
    t_event : list
        Per event, the first time at which it reaches 0 from below, root-found
        on the length of the step that crosses it; None if it never does.
    """

    t: np.ndarray
    records: list
    nit: np.ndarray
    t_event: list


class Monolith:
    """One monolith channel, eqs 9 to 18, twelve fields monolithic implicit.

    Parameters
    ----------
    params : OhCavendish
        Table I and the kinetics.
    feed : Feed
        The case: inlet state and initial wall temperature.
    n_x : int, optional
        Number of axial cells.
    nu_nusselt, sh : float, optional
        Limiting Nusselt and Sherwood numbers (square duct, H1 boundary
        condition, by default; ``sh`` defaults to ``nu_nusselt``).
    enthalpy_form : bool, optional
        True: eq 12 as printed, ``d(C_ps T_s)/dt``. False: the alternative
        reading ``C_ps dT_s/dt``.
    transfer_length : float, optional
        Length in ``h = Nu lambda_g / transfer_length`` and
        ``k_m = Sh D / transfer_length`` (eqs 17 and 18); ``2 R_h`` by default.
    tol, maxfev : float, int, optional
        ``newton`` step tolerance on the scaled unknowns, and iteration limit.
    damping : str, optional
        ``newton`` damping.
    outlet : {"outflow", "zero_gradient"}, optional
        Outlet of the convected gas fields. ``"outflow"`` (default): the face
        carries the last cell's value, exact upwinding for an equation first
        order in x. ``"zero_gradient"``: ``dc/dx = 0`` imposed on the face, a
        condition the gas equations do not have (kept to show what it costs).
    scaled : bool, optional
        Apply the constant two-sided diagonal scaling (default) or solve the
        system in physical units.
    dh : numpy.ndarray, optional
        Heats of combustion of CO, C3H6, CH4, H2, J/mol; :data:`DH` by default.
    """

    def __init__(self, params, feed, n_x=80, nu_nusselt=3.608, sh=None,
                 enthalpy_form=True, transfer_length=None, tol=1e-9, maxfev=60,
                 damping="monotone", outlet="outflow", scaled=True, dh=None):
        self.p, self.feed = params, feed
        self.n_x, self.n_c = n_x, N_FIELD
        self.nu_nusselt = nu_nusselt
        self.sh = nu_nusselt if sh is None else sh
        self.enthalpy_form = enthalpy_form
        self.transfer_length = 2.0 * params.r_h if transfer_length is None else transfer_length
        self.tol, self.maxfev, self.damping = tol, maxfev, damping
        self.outlet_kind, self.scaled = outlet, scaled
        self.dh = DH if dh is None else np.asarray(dh, float)
        self.length = params.length
        self.y_in = np.asarray(feed.y_in, float)
        self.s_geo = params.s_geo
        self.g_mass = params.g_mass                 # g/(cm2 s), constant
        self.g_mol = self.g_mass / M_GAS            # mol/(cm2 s), constant
        self.dt = 1.0
        self._build_grid()
        self._init_state()
        self._build_operators()

    # ---------------------------------------------------------------- grid
    def _build_grid(self):
        self.x_f, self.x_c = generate_grid(self.n_x, [0.0, self.length], generate_x_c=True)
        self.dx = self.length / self.n_x

    def _init_state(self):
        u = np.zeros((self.n_x, self.n_c))
        u[:, I_TS] = self.feed.T_s0
        u[:, I_TG] = self.feed.T_g_in
        u[:, I_YG] = self.y_in
        u[:, I_YS] = self.y_in
        self.u = u
        self.t = 0.0

    # ----------------------------------------------------------- operators
    def _build_operators(self):
        shape = (self.n_x, self.n_c)
        shape_f = (self.n_x + 1, self.n_c)

        # eq 15, adiabatic solid at both ends, OUTWARD normal:
        #   x = 0:  -dTs/dx = 0  -> a = 1, b = 0, d = 0
        #   x = L:  +dTs/dx = 0  -> a = 1, b = 0, d = 0
        bc_diff = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        self.bc_diff = bc_diff
        grad_mat, grad_bc = construct_grad(shape, self.x_f, self.x_c, bc_diff, axis=0)
        lam = np.zeros(self.n_c)
        lam[I_TS] = self.p.lam_s * (1.0 - self.p.eps)   # conduction on T_s only
        lam_mat = construct_coefficient_matrix(lam, shape=shape_f, axis=-1)
        div_mat = construct_div(shape, self.x_f, nu=0, axis=0)   # nu = 0: Cartesian

        # constant carriers: mass flux for T_g, molar flux for the mole fractions;
        # zero on the wall fields, so the convection operator leaves them alone
        v = np.zeros(self.n_c)
        v[I_TG] = self.g_mass * self.p.cp_g
        v[I_YG] = self.g_mol
        d_in = np.zeros(self.n_c)
        d_in[I_TG] = self.feed.T_g_in                  # eq 14: T_g(0) = T_g,in
        d_in[I_YG] = self.y_in                         # eq 13: c_g,i(0) = c_g,i,in
        # inlet x = 0: Dirichlet b*c = d on every field (only the convected ones
        # use it); outlet x = L: the gas equations are first order in x and need
        # no condition there, so the face carries the last cell's value
        if self.outlet_kind == "outflow":
            bc_out = {"outflow": True}
        elif self.outlet_kind == "zero_gradient":     # x = L: dc/dx = 0
            bc_out = {"a": np.ones(self.n_c), "b": np.zeros(self.n_c), "d": np.zeros(self.n_c)}
        else:
            raise ValueError(f"outlet must be 'outflow' or 'zero_gradient', not {self.outlet_kind!r}")
        bc_conv = ({"a": np.zeros(self.n_c), "b": np.ones(self.n_c), "d": d_in}, bc_out)
        conv_mat, conv_bc = construct_convflux_upwind(shape, self.x_f, self.x_c,
                                                      bc_conv, v=v, axis=0)
        self.bc_conv = bc_conv

        self.jac_const = (-div_mat @ (lam_mat @ grad_mat) + div_mat @ conv_mat).tocsc()
        self.g_const = np.asarray(
            (-div_mat @ (lam_mat @ grad_bc) + div_mat @ conv_bc).todense()).reshape(-1, 1)
        self.numjac = NumJac(shape)        # last axis in full: every term is pointwise
        self._build_scaling(shape)

    def _build_scaling(self, shape):
        """Constant two-sided diagonal scaling of the Newton system.

        Columns: temperatures by 1000 K, mole fractions by their inlet value.
        Rows: each equation by its film-transfer coefficient at 600 K times the
        column scale, so every scaled row and unknown is of order one.
        """
        t_ref = 1000.0
        y_ref = np.maximum(self.y_in, 1e-5)
        self.col_scale = np.concatenate(([t_ref, t_ref], y_ref, y_ref))
        h_ref = self.nu_nusselt * self.p.lam_gas(600.0) / self.transfer_length * self.s_geo
        km_ref = (self.sh * SB_COEF * 600.0 ** SB_B / self.transfer_length * self.s_geo
                  * self.p.p_tot / (R_GAS * 600.0) / 1e6)
        self.row_scale = np.concatenate(([1 / (h_ref * t_ref)] * 2,
                                         1 / (km_ref * y_ref), 1 / (km_ref * y_ref)))
        if not self.scaled:
            self.col_scale = np.ones(self.n_c)
            self.row_scale = np.ones(self.n_c)
            t_ref, y_ref = 1.0, np.ones(5)
        # bounds where the residual stops being meaningful: T > 0 for exp(-E/T)
        # and T**1.823, and 0 <= y <= y_in, the range a consumed-only species
        # can take (below 0 the inhibition denominator G of eq 6 can vanish)
        self._lb = np.concatenate(([250 / t_ref] * 2, np.zeros(10)))
        self._ub = np.concatenate(([3000 / t_ref] * 2, self.y_in / y_ref, self.y_in / y_ref))
        cs = np.broadcast_to(self.col_scale, shape).reshape(-1)
        rs = np.broadcast_to(self.row_scale, shape).reshape(-1)
        self._Dc, self._Dr = diags_array(cs, format="csc"), diags_array(rs, format="csc")
        self._rs = rs.reshape(-1, 1)

    # ------------------------------------------------------------- physics
    def rates(self, y_s, T_s):
        """Specific rates of eqs 1 to 8, mol/(cm2 Pt s); order CO C3H6 CH4 H2 O2.

        Parameters
        ----------
        y_s : numpy.ndarray
            Surface mole fractions, last axis the five species.
        T_s : numpy.ndarray
            Wall temperature, K, the shape of ``y_s`` without its last axis.

        Returns
        -------
        numpy.ndarray
            Rates, the shape of ``y_s``.
        """
        k1, k2, k3 = (a * np.exp(-e / T_s) for a, e in self.p.kin)
        K1, K2, K3, K4 = (a * np.exp(-e / T_s) for a, e in self.p.ads)
        co, hc, ch4, h2, o2 = (y_s[..., i] for i in range(5))
        G = (T_s * (1 + K1 * co + K2 * hc) ** 2
             * (1 + K3 * co ** 2 * hc ** 2) * (1 + K4 * self.feed.y_no ** 0.7))
        r_co, r_hc = k1 * co * o2 / G, k2 * hc * o2 / G
        r_ch4, r_h2 = k3 * ch4 * o2 / G, k1 * h2 * o2 / G   # eq 7: H2 uses k1
        r_o2 = 0.5 * r_co + 4.5 * r_hc + 2.0 * r_ch4 + 0.5 * r_h2   # eq 8
        return np.stack([r_co, r_hc, r_ch4, r_h2, r_o2], axis=-1)

    def transport(self, Tg):
        """Film coefficients and molar density at gas temperature ``Tg``.

        Returns
        -------
        h : numpy.ndarray
            Heat-transfer coefficient, eq 17, J/(cm2 s K).
        km : numpy.ndarray
            Mass-transfer coefficients, eq 18, cm/s, last axis the species.
        c_tot : numpy.ndarray
            Molar density ``P / (R T)``, mol/cm3.
        """
        h = self.nu_nusselt * self.p.lam_gas(Tg) / self.transfer_length
        d_i = SB_COEF * (Tg ** SB_B)[..., None]
        km = self.sh * d_i / self.transfer_length
        c_tot = self.p.p_tot / (R_GAS * Tg) / 1e6
        return h, km, c_tot

    def enthalpy(self, T):
        """Solid enthalpy per gram, J/g, as it sits inside the time derivative.

        Eq 12 as printed puts ``C_ps`` inside ``d/dt``: ``C_ps(T) T``. The
        alternative reading, ``C_ps dT/dt``, integrates ``C_ps``.
        """
        c0, c1, c2 = self.p.cps
        if self.enthalpy_form:
            return self.p.cp_solid(T) * T
        return c0 * T + 0.5 * c1 * T ** 2 - c2 / T

    def source(self, u):
        """Every pointwise term of the residual: accumulation, film transfer, reaction."""
        Ts, Tg = u[..., I_TS], u[..., I_TG]
        y_g, y_s = u[..., I_YG], u[..., I_YS]
        h, km, c_tot = self.transport(Tg)
        rr = self.rates(y_s, Ts)
        exch = km * self.s_geo * c_tot[..., None] * (y_g - y_s)    # mol/(cm3 s)
        q = self.p.a_pt * (self.dh * rr[..., :4]).sum(axis=-1)           # J/(cm3 s)
        acc = (1 - self.p.eps) * self.p.rho_s * self.enthalpy(Ts)
        f = np.zeros_like(u)
        f[..., I_TS] = -acc / self.dt + h * self.s_geo * (Tg - Ts) + q   # eq 12
        f[..., I_TG] = h * self.s_geo * (Ts - Tg)                        # eq 10
        f[..., I_YG] = -exch                                             # eq 9
        f[..., I_YS] = -self.p.a_pt * rr + exch                          # eq 11
        return f

    def outlet(self, u=None):
        """Gas state on the outlet face, as the convection operator transports it."""
        u = self.u if u is None else u
        vals, _ = compute_boundary_values(u, self.x_f, self.x_c, bc=self.bc_conv[1],
                                          axis=0, bound_id=1)
        return np.asarray(vals).ravel()

    def wall_ends(self, u=None):
        """Wall temperature on the inlet and outlet faces, K, from the wall's own bc."""
        u = self.u if u is None else u
        t0, _, t1, _ = compute_boundary_values(u[:, I_TS:I_TS + 1], self.x_f, self.x_c,
                                               self.bc_diff, axis=0)
        return float(np.ravel(t0)[0]), float(np.ravel(t1)[0])

    def solid_enthalpy(self, u=None):
        """Wall enthalpy per frontal area, J/cm2."""
        u = self.u if u is None else u
        return float(((1 - self.p.eps) * self.p.rho_s
                      * self.enthalpy(u[:, I_TS])).sum() * self.dx)

    # -------------------------------------------------------------- solver
    def residual(self, w, g_old):
        """Scaled residual and Jacobian of one backward-Euler step.

        Parameters
        ----------
        w : numpy.ndarray
            Scaled state ``u / col_scale``, shape ``(n_x, 12)``.
        g_old : numpy.ndarray
            Old-time accumulation term, flat column.
        """
        u = w.reshape(self.u.shape) * self.col_scale
        f, jac_f = self.numjac(self.source, u)
        g = self.jac_const @ u.reshape(-1, 1) + self.g_const - f.reshape(-1, 1) - g_old
        return self._rs * g, (self._Dr @ (self.jac_const - jac_f) @ self._Dc).tocsc()

    def step_from(self, u_old, dt):
        """One backward-Euler step of length ``dt`` from ``u_old``, model state untouched.

        Returns
        -------
        u_new : numpy.ndarray
        result : scipy.optimize.OptimizeResult
            The ``newton`` result.

        Raises
        ------
        StepFailed
            When ``newton`` returns without converging.
        pymrm.NewtonBreakdown
            When no Newton step can be formed.
        """
        self.dt = dt
        g_old = np.zeros(u_old.shape)
        g_old[:, I_TS] = (1 - self.p.eps) * self.p.rho_s * self.enthalpy(u_old[:, I_TS]) / dt
        g_old = g_old.reshape(-1, 1)
        r = newton(lambda w: self.residual(w, g_old), u_old / self.col_scale,
                   tol=self.tol, maxfev=self.maxfev, damping=self.damping,
                   bounds=(self._lb, self._ub))
        if not r.success:
            raise StepFailed(f"time step dt = {dt:g} s from t = {self.t:g} s did not "
                             f"converge: status {r.status}, {r.message}", r)
        return r.x.reshape(u_old.shape) * self.col_scale, r

    def step(self, dt):
        """Advance the model state by one backward-Euler step; returns the result."""
        self.u, r = self.step_from(self.u, dt)
        self.t += dt
        return r

    def march(self, t_end, dt, record=None, events=()):
        """March to ``t_end`` with steps ``dt``, every step's convergence checked.

        Parameters
        ----------
        t_end, dt : float
            End time and step, s.
        record : callable, optional
            ``record(model)`` at every time level, including t = 0.
        events : sequence of callable, optional
            ``event(u) -> float`` each. The first crossing of each from below 0
            to 0 or above is located by ``brentq`` on the length of the step
            that crosses it, every trial a converged backward-Euler step from
            the level before.

        Returns
        -------
        March
            ``t_event`` is a list with one entry per event (None if it never
            crosses).
        """
        ts, out, nit = [self.t], [], []
        if record is not None:
            out.append(record(self))
        t_event = [None] * len(events)
        g_prev = [ev(self.u) for ev in events]
        for _ in range(int(round(t_end / dt))):
            u_prev, t_prev = self.u, self.t
            nit.append(self.step(dt).nit)
            ts.append(self.t)
            if record is not None:
                out.append(record(self))
            for i, ev in enumerate(events):
                if t_event[i] is not None:
                    continue
                g_now = ev(self.u)
                if g_prev[i] < 0.0 <= g_now:
                    t_event[i] = t_prev + self.crossing(u_prev, dt, ev, g_prev[i])
                    self.dt = dt
                g_prev[i] = g_now
        return March(np.array(ts), out, np.array(nit), t_event)

    def crossing(self, u_prev, dt, event, g_prev, xtol=1e-10):
        """Step length in ``[0, dt]`` at which ``event`` reaches 0, by ``brentq``.

        At length 0 the state is ``u_prev`` itself (``event`` there is
        ``g_prev < 0``); at ``dt`` it is the level where the crossing was seen.
        """
        def g(tau):
            return g_prev if tau == 0.0 else event(self.step_from(u_prev, tau)[0])
        return brentq(g, 0.0, dt, xtol=xtol)

    def peak(self, u_prev, dt, xatol=1e-7):
        """Refine a sampled maximum of the wall temperature in time.

        Maximises ``max_x T_s`` over the length ``tau`` in ``[0, 2 dt]`` of a
        backward-Euler step from ``u_prev``, the level before the sampled
        maximum (``tau = dt`` reproduces the sampled level).

        Returns
        -------
        tau : float
            Step length at the maximum, s.
        T_max : float
            Peak wall temperature, K.
        j : int
            Cell of the peak.
        """
        def neg(tau):
            return -self.step_from(u_prev, tau)[0][:, I_TS].max()
        r = minimize_scalar(neg, bounds=(0.0, 2.0 * dt), method="bounded",
                            options={"xatol": xatol})
        if not r.success:
            raise RuntimeError(f"peak refinement did not converge: {r.message}")
        if not (1e-3 * dt < r.x < 1.999 * dt):
            raise RuntimeError(f"peak refinement ended on its bracket, tau = {r.x:g}")
        u = self.step_from(u_prev, r.x)[0]
        return float(r.x), float(u[:, I_TS].max()), int(u[:, I_TS].argmax())
