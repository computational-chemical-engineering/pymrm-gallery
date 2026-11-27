"""Vanden Bussche and Froment (1996) methanol and RWGS kinetics in a plug-flow bed.

Model code of gallery page C2.5. The kinetics know nothing about the reactor
and the reactor knows nothing about the chemistry beyond a pointwise source,
so either piece is reusable on its own.

Contents
--------
R_GAS, SPECIES, NU, ATOMS, MOLAR_MASS, EQUIL
    Gas constant, species order, stoichiometry of the two reactions, element
    content, molar masses and the printed equilibrium correlations (p. 5).
table2, table1_to_table2
    Table 2 as a parameter dictionary; the Table 1 to Table 2 conversion.
Kinetics
    Eq. 3: methanol from CO2 and the reverse water gas shift, with the
    equilibrium constants K1* and K3*.
Thermo
    Heats of reaction and heat capacities (NIST WebBook Shomate and TRC
    tables), or the constant mass heat capacity of the test model.
PlugFlowBed
    Adiabatic or isothermal, isobaric plug flow in the reduced catalyst mass
    ``zeta = W/W_bed``: first-order upwind convection and ``pymrm.newton``.
reference_solution
    The same ODEs by ``scipy.integrate.solve_ivp`` (Radau), with dense output.
equilibrium, adiabatic_equilibrium
    Two-reaction equilibrium at fixed T and p, and at the adiabatic invariant
    of the constant-heat-capacity test model.
half_units
    Distance of a value from a printed entry in half-units of its last digit.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)
from pymrm.solve import NewtonBreakdown
from pymrm.checks import find_roots, residual_check

__all__ = ["R_GAS", "SPECIES", "NU", "ATOMS", "MOLAR_MASS", "EQUIL", "GROUP_KEYS", "table2",
           "table1_to_table2", "Kinetics", "Thermo", "PlugFlowBed", "reference_solution",
           "equilibrium", "adiabatic_equilibrium", "half_units"]

R_GAS = 8.314462618  # J/(mol K), CODATA 2018 (Table 1 to Table 2 needs 8.31413 to 8.31463)
SPECIES = ("CO", "CO2", "H2", "MeOH", "H2O", "inert")
# Columns: methanol synthesis CO2 + 3 H2 -> CH3OH + H2O; RWGS CO2 + H2 -> CO + H2O
NU = np.array([[0.0, 1.0],     # CO
               [-1.0, -1.0],   # CO2
               [-3.0, -1.0],   # H2
               [1.0, 0.0],     # CH3OH
               [1.0, 1.0],     # H2O
               [0.0, 0.0]])    # inert
ATOMS = np.array([[1, 1, 0, 1, 0, 0],    # C
                  [0, 0, 2, 4, 2, 0],    # H
                  [1, 2, 0, 1, 1, 0]],   # O
                 dtype=float)
MOLAR_MASS = {"CO": 0.0280101, "CO2": 0.0440095, "H2": 0.00201588, "MeOH": 0.03204186,
              "H2O": 0.01801528, "Ar": 0.039948, "N2": 0.0280134}   # kg/mol
# log10 K1* = a/T + b (bar^-2), log10 (1/K3*) = a/T + b (p. 5, from Graaf et al. 1986)
EQUIL = {"K1_slope": 3066.0, "K1_icpt": -10.592, "K3inv_slope": -2073.0, "K3inv_icpt": 2.029}
# Table rows (group column) to the keys used here
GROUP_KEYS = {"sqrt(K_H2)": "sqrtKH2", "K_H2O": "KH2O", "K_H2O/(K8 K9 K_H2)": "Kw",
              "k'5a K'2 K3 K4 K_H2": "k5", "k'1": "k1"}


def table2(df):
    """Table 2 as ``{"A_<key>": A, "B_<key>": B}`` (B in J/mol; the printed dash is 0)."""
    par = {}
    for r in df.itertuples():
        key = GROUP_KEYS[r.group]
        v = r.value_printed
        par[f"{r.param}_{key}"] = 0.0 if (not isinstance(v, str) or not v.strip()) else float(v)
    return par


def table1_to_table2(a_star, b, r_gas=R_GAS, t_av=501.57):
    """Frequency factor ``A = A* exp(-B/(R T_av))`` of the reparametrised Table 1 form."""
    return a_star * np.exp(-b / (r_gas * t_av))


def half_units(printed, value):
    """``|value - printed|`` in half-units of the last printed digit of ``printed`` (a string)."""
    s = printed.strip().lower().replace(",", "")
    mant, _, exp = s.partition("e")
    decimals = len(mant.split(".")[1]) if "." in mant else 0
    unit = 10.0 ** (-decimals + (int(exp) if exp else 0))
    return abs(value - float(s)) / (0.5 * unit)


class Kinetics:
    """Eq. 3 (p. 6): rates of methanol synthesis from CO2 and of the RWGS.

    Parameters
    ----------
    par : dict
        ``A_<key>`` and ``B_<key>`` for the keys of ``GROUP_KEYS``; ``kappa = A exp(B/(R T))``.
    r_gas : float, optional
        Gas constant (J/(mol K)) in the exponentials.
    factor : float, optional
        Multiplies both rates (the page's rate factor F).
    equil : dict, optional
        Coefficients of the two equilibrium correlations.
    drop : iterable of str, optional
        Denominator groups set to zero (``"KH2O"``, ``"Kw"``, ``"sqrtKH2"``).
    p_unit : float, optional
        Bar per pressure unit fed to the kinetic groups (1: bar, as printed).
    """

    # Table 2, the gas constant, the rate factor, the equilibrium correlations, removed groups
    # and the pressure unit are separate readings of the paper.
    # codehealth: allow params -- each argument is one reading of the paper
    def __init__(self, par, r_gas=R_GAS, factor=1.0, equil=EQUIL, drop=(), p_unit=1.0):
        self.par, self.r_gas, self.factor = dict(par), float(r_gas), float(factor)
        self.equil, self.drop, self.p_unit = dict(equil), frozenset(drop), float(p_unit)

    def constants(self, T):
        """The five groups of Table 2 and K1* (bar^-2), K3* (-) at temperature ``T`` (K)."""
        T = np.asarray(T, dtype=float)
        c = {k: self.par["A_" + k] * np.exp(self.par["B_" + k] / (self.r_gas * T))
             for k in GROUP_KEYS.values()}
        for k in self.drop:
            c[k] = 0.0 * c[k]
        e = self.equil
        c["K1"] = 10.0 ** (e["K1_slope"] / T + e["K1_icpt"])
        c["K3"] = 10.0 ** (-(e["K3inv_slope"] / T + e["K3inv_icpt"]))
        return c

    def denominator_terms(self, p, T):
        """The four terms of ``1/beta`` (eq. 3), shape ``(..., 4)``; ``p`` in bar, last axis species."""
        c = self.constants(T)
        q = p / self.p_unit
        return np.stack([np.ones_like(q[..., 2]), c["Kw"] * q[..., 4] / q[..., 2],
                         c["sqrtKH2"] * np.sqrt(q[..., 2]), c["KH2O"] * q[..., 4]], axis=-1)

    def approach(self, p, T):
        """Reaction quotient over equilibrium constant of the two reactions, shape ``(..., 2)``.

        Methanol: ``p_W p_M / (K1* p_H2^3 p_CO2)``; RWGS: ``K3* p_W p_CO / (p_CO2 p_H2)``.
        """
        c = self.constants(T)
        pco, pco2, ph2, pm, pw = (p[..., i] for i in range(5))
        return np.stack([pw * pm / (c["K1"] * ph2 ** 3 * pco2),
                         c["K3"] * pw * pco / (pco2 * ph2)], axis=-1)

    def rates(self, p, T):
        """``(r_MeOH, r_RWGS)`` in mol/(kg_cat s), shape ``(..., 2)``; ``p`` in bar.

        The approach terms are multiplied out, so both rates are smooth and exactly zero
        at ``p_CO2 = p_W = 0``.
        """
        c = self.constants(T)
        pco, pco2, ph2, pm, pw = (p[..., i] for i in range(5))
        q = p / self.p_unit
        den = 1.0 + c["Kw"] * q[..., 4] / q[..., 2] + c["sqrtKH2"] * np.sqrt(q[..., 2]) \
            + c["KH2O"] * q[..., 4]
        # k5 q_CO2 q_H2 (1 - approach) / den^3 with q = p / p_unit, the approach in bar
        r1 = c["k5"] * (pco2 * ph2 - pw * pm / (c["K1"] * ph2 ** 2)) / self.p_unit ** 2 / den ** 3
        r2 = c["k1"] * (pco2 - c["K3"] * pw * pco / ph2) / self.p_unit / den
        return self.factor * np.stack([r1, r2], axis=-1)


class Thermo:
    """Heats of reaction and the heat-capacity flow of the gas.

    Parameters
    ----------
    shomate : pandas.DataFrame
        NIST Shomate rows (species, T_min_K, T_max_K, A..H), read with ``dtype=str``.
    meoh_cp : pandas.DataFrame
        Methanol ideal-gas heat capacity table (T_K, Cp_J_per_molK).
    inert : {"Ar", "N2"}
        The unnamed inert of Table 3.
    heats : {"H298", "HvH"} or tuple
        ``"H298"``: standard heats at 298.15 K from the formation enthalpies;
        ``"HvH"``: the van 't Hoff heats of the p. 5 correlations; a tuple gives
        ``(dH_MeOH, dH_RWGS)`` in J/mol directly.
    cp_scale : float, optional
        Multiplies every molar heat capacity.
    cp_mass : float, optional
        Constant mass heat capacity (J/(kg K)) of the test model; then the
        heat-capacity flow is ``m c_p``, independent of T and composition.
    dhf_meoh : float, optional
        Formation enthalpy of gaseous methanol, J/mol.
    """

    # codehealth: allow params -- two data tables and the unprinted thermodynamic readings
    def __init__(self, shomate, meoh_cp, inert="Ar", heats="H298", cp_scale=1.0, cp_mass=None,
                 dhf_meoh=-201.9e3, r_gas=R_GAS):
        self.inert, self.cp_scale, self.cp_mass = inert, float(cp_scale), cp_mass
        self.rows = {}
        for r in shomate.itertuples():
            coef = np.array([float(getattr(r, k)) for k in "ABCDEFGH"])
            self.rows.setdefault(r.species, []).append((float(r.T_min_K), coef))
        self.meoh = CubicSpline(meoh_cp.T_K.astype(float), meoh_cp.Cp_J_per_molK.astype(float))
        hf = {s: 1e3 * self.rows[s][0][1][7] for s in ("CO", "CO2", "H2", "H2O")}
        hf["MeOH"] = float(dhf_meoh)
        if isinstance(heats, str) and heats == "H298":
            self.dH = np.array([hf["MeOH"] + hf["H2O"] - hf["CO2"] - 3 * hf["H2"],
                                hf["CO"] + hf["H2O"] - hf["CO2"] - hf["H2"]])
        elif isinstance(heats, str) and heats == "HvH":
            ln10 = np.log(10.0)
            self.dH = np.array([-r_gas * ln10 * EQUIL["K1_slope"],
                                r_gas * ln10 * (-EQUIL["K3inv_slope"])])
        else:
            self.dH = np.asarray(heats, dtype=float)
        names = ("CO", "CO2", "H2", "MeOH", "H2O", inert)
        self.m_feed_of = lambda y: float(np.dot(y, [MOLAR_MASS[s] for s in names]))

    def _shomate(self, species, T):
        rows = self.rows[species]
        t = T / 1000.0
        cp = np.zeros_like(t)
        for i, (tmin, c) in enumerate(rows):
            use = T >= tmin if i else np.ones_like(T, dtype=bool)
            cp = np.where(use, c[0] + c[1] * t + c[2] * t ** 2 + c[3] * t ** 3 + c[4] / t ** 2, cp)
        return cp

    def cp(self, T):
        """Molar heat capacities (J/(mol K)) in species order, shape ``T.shape + (6,)``."""
        T = np.asarray(T, dtype=float)
        cps = [self._shomate("CO", T), self._shomate("CO2", T), self._shomate("H2", T),
               self.meoh(T), self._shomate("H2O", T), self._shomate(self.inert, T)]
        return self.cp_scale * np.stack(cps, axis=-1)

    def capacity(self, u, T, m_feed):
        """Heat-capacity flow per mol/s of feed, J/(K mol_feed): ``sum u_i Cp_i(T)`` or ``M c_p``."""
        if self.cp_mass is not None:
            return m_feed * self.cp_mass * np.ones_like(np.asarray(T, dtype=float))
        return np.sum(u * self.cp(T), axis=-1)


class PlugFlowBed:
    """Adiabatic or isothermal, isobaric plug flow in ``zeta = W/W_bed``.

    ``d(v u)/dzeta = s(u)`` with ``v = 1``, ``u = (F_i/F0 for the six species, T/T0)``:
    ``du_i/dzeta = (W/F0) sum_j nu_ij r_j`` and, adiabatic,
    ``d(T/T0)/dzeta = (W/F0) sum_j (-dH_j) r_j / (C T0)`` with ``C`` the heat-capacity
    flow per mol/s of feed. State layout ``(n, 7)``: the coordinate first, fields last.

    Parameters
    ----------
    kin : Kinetics
    y0 : array_like
        Feed mole fractions in species order.
    T0, p : float
        Inlet temperature (K) and pressure (bar).
    w_f0 : float
        ``W/F0`` of the bed (kg s/mol).
    thermo : Thermo or None
        None: isothermal at ``T0``.
    n : int
        Number of cells.
    zeta_max : float
        Length of the computed bed in units of ``W`` (a longer bed for the
        rate-factor scaling, or a long bed for the equilibrium plateau).
    nu : array_like, optional
        Stoichiometric matrix (species x reactions).
    """

    # codehealth: allow params -- the bed: kinetics, feed, T0, p, W/F0, energy model, grid, stoichiometry
    def __init__(self, kin, y0, T0, p, w_f0, thermo=None, n=1000, zeta_max=1.0, nu=NU):
        self.kin, self.thermo, self.nu = kin, thermo, np.asarray(nu, dtype=float)
        self.y0, self.T0, self.p, self.w_f0 = np.asarray(y0, float), float(T0), float(p), float(w_f0)
        self.m_feed = thermo.m_feed_of(self.y0) if thermo is not None else np.nan
        self.shape = (int(n), 7)
        self.z_f = np.linspace(0.0, zeta_max, self.shape[0] + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.u_in = np.append(self.y0, 1.0)
        # zeta = 0, outward normal -zeta: the feed enters, u = u_in   -> a = 0, b = 1, d = u_in
        # zeta = zeta_max, outward normal +zeta: pure convective outlet, nothing imposed
        self.bc = ({"a": 0.0, "b": 1.0, "d": self.u_in.reshape(1, -1)}, {"outflow": True})
        self.u, self.result = np.tile(self.u_in, (self.shape[0], 1)), None
        self._build_operators()

    def _build_operators(self):
        # The convected quantity is u with unit velocity in zeta.
        # nu = 0: the catalyst-mass coordinate is Cartesian.
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c, self.bc,
                                                      v=1.0, axis=0)
        div_mat = construct_div(self.shape, self.z_f, nu=0, axis=0)
        self.jac_const = div_mat @ conv_mat
        self.g_const = (div_mat @ conv_bc).toarray().reshape(-1, 1)
        # The source at a cell depends on that cell's seven fields only: the default
        # stencil (last axis coupled in full) is exact, a 7 x 7 block per cell.
        self.numjac = NumJac(self.shape)

    def source(self, u):
        """Pointwise ``s(u)``, shape ``u.shape``.

        Newton's trial states can leave the physical domain (a flow near zero); the rates
        there overflow, and newton rejects the non-finite residual and shortens the step.
        Those floating-point warnings are silenced here, so they never reach the output.
        """
        with np.errstate(all="ignore"):
            return self._source(u)

    def _source(self, u):
        y = u[..., :6] / np.sum(u[..., :6], axis=-1, keepdims=True)
        T = self.T0 * u[..., 6]
        r = self.kin.rates(self.p * y, T)
        s_species = self.w_f0 * (r @ self.nu.T)
        if self.thermo is None:
            s_T = np.zeros_like(T)
        else:
            q = r @ (-self.thermo.dH)
            s_T = self.w_f0 * q / (self.thermo.capacity(u[..., :6], T, self.m_feed) * self.T0)
        return np.concatenate([s_species, s_T[..., None]], axis=-1)

    def residual(self, u):
        """Steady residual ``div(v u) - s(u)`` and its Jacobian."""
        g_s, jac_s = self.numjac(self.source, u.reshape(self.shape))
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_s.reshape((-1, 1))
        return g, self.jac_const - jac_s

    # codehealth: allow params -- Newton controls passed through to pymrm.newton
    def solve(self, guess=None, maxfev=60, tol=1e-10, damping="monotone", strict=True):
        """Newton from the feed state (or ``guess``); raises ``RuntimeError`` unless it converges.

        ``strict=False`` keeps whatever state one Newton run from the guess stopped at.
        Molar flows are bounded below by 0 and T/T0 by 0.5. When Newton from the flat
        feed state fails (a front far from it), the bed is first solved with W/F0
        divided by 2^k, k = 6 down to 1, each solution the guess of the next
        (continuation in the bed length, on the same grid).
        """
        lower = np.array([0.0] * 6 + [0.5])

        def attempt(x0, bed):
            try:
                res = newton(bed.residual, x0, maxfev=maxfev, tol=tol, stagnation=True,
                             bounds=(lower, None), damping=damping)
            except NewtonBreakdown:
                return None
            return res if res.success else None

        x0 = self.u if guess is None else np.asarray(guess, float).reshape(self.shape)
        if not strict:
            res = newton(self.residual, x0, maxfev=maxfev, tol=tol, stagnation=True,
                         bounds=(lower, None), damping=damping)
            self.u, self.result = res.x.reshape(self.shape), res
            return self
        result = attempt(x0, self)
        if result is None and guess is None:
            x = self.u.copy()
            for k in range(6, -1, -1):
                bed = self if k == 0 else PlugFlowBed(
                    self.kin, self.y0, self.T0, self.p, self.w_f0 / 2 ** k, self.thermo,
                    n=self.shape[0], zeta_max=self.z_f[-1], nu=self.nu)
                result = attempt(x, bed)
                if result is None:
                    break
                x = result.x.reshape(self.shape)
        if result is None:
            raise RuntimeError(f"PlugFlowBed.solve on {self.shape[0]} cells did not converge")
        self.u, self.result = result.x.reshape(self.shape), result
        return self

    def backward_error(self):
        """Componentwise backward error of the current state (``pymrm.checks``)."""
        return residual_check(self.residual, self.u.reshape(-1))["backward_error"]

    def outlet(self):
        """``u`` on the outlet face, read with ``compute_boundary_values``."""
        _, _, u_out, _ = compute_boundary_values(self.u, self.z_f, self.z_c, self.bc, axis=0)
        return np.asarray(u_out, dtype=float).ravel()

    def at(self, zeta):
        """``u`` at ``zeta``, linear between face values.

        First-order upwind is implicit Euler in ``zeta``: cell ``i`` holds the value at
        its downstream face, so the nodes are the faces, with the inlet value at the
        first face and the outlet face read by ``compute_boundary_values``. Read this
        way the error is the scheme's alone, first order without cancellation.
        """
        zeta = np.asarray(zeta, dtype=float)
        if np.any(zeta < 0.0) or np.any(zeta > self.z_f[-1]):
            raise ValueError("PlugFlowBed.at: zeta outside the computed bed")
        nodes = self.z_f
        vals = np.vstack([self.u_in, self.u[:-1], self.outlet()])
        return np.stack([np.interp(zeta, nodes, vals[:, j]) for j in range(7)], axis=-1)

    def mole_fractions(self, u):
        return u[..., :6] / np.sum(u[..., :6], axis=-1, keepdims=True)

    def rates_at(self, u):
        """``(r_MeOH, r_RWGS)`` at states ``u``."""
        return self.kin.rates(self.p * self.mole_fractions(u), self.T0 * u[..., 6])


def reference_solution(bed, zeta_end=None, rtol=1e-12, atol=1e-14):
    """``bed``'s ODEs by ``solve_ivp`` (Radau, dense output); raises unless it succeeds.

    Shares with the pymrm route only ``bed.source`` (rates, stoichiometry, energy balance).
    """
    zeta_end = bed.z_f[-1] if zeta_end is None else zeta_end

    def rhs(_, u):
        return bed.source(u.reshape(1, -1)).ravel()

    sol = solve_ivp(rhs, (0.0, zeta_end), bed.u_in, method="Radau", rtol=rtol, atol=atol,
                    dense_output=True)
    if sol.status != 0:
        raise RuntimeError(f"reference_solution: {sol.message}")
    return sol


# codehealth: allow params -- find_roots arguments and a label for the error message
def _one_root(f, lo, hi, what, n_scan=16, xtol=1e-14):
    res = find_roots(f, lo, hi, n_scan=n_scan, xtol=xtol)
    if len(res["roots"]) != 1:
        raise ValueError(f"{what}: {len(res['roots'])} roots on [{lo}, {hi}] ({res['message']})")
    return res["roots"][0]


# codehealth: allow params -- kinetics (for K1*, K3*), feed, T, p and the RWGS switch
def equilibrium(kin, y0, T, p, shift=True):
    """Two-reaction equilibrium at fixed ``T`` (K) and ``p`` (bar).

    Unknowns: extents ``a`` (methanol) and ``b`` (RWGS) per mol of feed. For each
    ``a`` the RWGS condition fixes ``b`` (monotone in ``b``); the methanol
    condition then fixes ``a``. Both by ``find_roots`` with exactly one root.
    ``shift=False`` freezes the RWGS (``b = 0``).

    Returns
    -------
    dict
        ``a``, ``b``, mole fractions ``y`` and moles ``n`` per mol of feed.
    """
    y0 = np.asarray(y0, float)
    c = kin.constants(T)
    ln_k1, ln_krwgs = np.log(c["K1"]), -np.log(c["K3"])
    eps = 1e-13

    def moles(a, b):
        return y0 + NU @ np.array([a, b])

    def b_of(a):
        if not shift:
            return 0.0
        lo = -min(y0[0], y0[4] + a) * (1 - eps) + 0.0
        hi = min(y0[1] - a, (y0[2] - 3 * a)) * (1 - eps)

        def g(b):
            n = moles(a, b)
            return np.log(n[0] * n[4]) - np.log(n[1] * n[2]) - ln_krwgs
        return _one_root(g, lo, hi, "RWGS extent", n_scan=8)

    def h(a):
        b = b_of(a)
        n = moles(a, b)
        y = n / n.sum()
        return np.log(y[3] * y[4]) - np.log(y[2] ** 3 * y[1] * p ** 2) - ln_k1

    a_hi = (min(y0[1] + y0[0], y0[2] / 2.0) if shift else min(y0[1], y0[2] / 3.0)) * (1 - 1e-9)
    # scan points outside the domain give nan (log of a negative amount), never a sign
    with np.errstate(all="ignore"):
        a = _one_root(h, 1e-12, a_hi, "methanol extent", n_scan=24)
        b = b_of(a)
    n = moles(a, b)
    return {"a": a, "b": b, "n": n, "y": n / n.sum()}


# codehealth: allow params -- kinetics, energy model, feed, inlet state and the T bracket
def adiabatic_equilibrium(kin, thermo, y0, T0, p, T_lo=None, T_hi=None):
    """Equilibrium of the constant-heat-capacity test model: ``T`` on the energy invariant.

    ``T - T0 = sum_j (-dH_j) xi_j / (M_feed c_p)`` with the extents ``xi`` of
    ``equilibrium(kin, y0, T, p)``; only valid with ``thermo.cp_mass`` set.
    """
    if thermo.cp_mass is None:
        raise ValueError("adiabatic_equilibrium needs the constant mass heat capacity")
    m_feed = thermo.m_feed_of(np.asarray(y0, float))

    def f(T):
        e = equilibrium(kin, y0, T, p)
        return T - T0 - (-thermo.dH) @ np.array([e["a"], e["b"]]) / (m_feed * thermo.cp_mass)

    T = _one_root(f, T0 + 1e-6 if T_lo is None else T_lo, T0 + 200.0 if T_hi is None else T_hi,
                  "adiabatic equilibrium temperature", n_scan=32)
    e = equilibrium(kin, y0, T, p)
    e["T"] = T
    return e
