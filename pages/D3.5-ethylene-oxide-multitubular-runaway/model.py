"""Borman and Westerterp (1995) model III kinetics in a wall-cooled plug-flow tube.

The rate law of model III (Table 5, caption form, constants of Table 6; the same reading as
gallery page C2.11), the heat release of epoxidation and combustion, the inlet of ethene in dry
air, and a steady 1-D pseudo-homogeneous tube solved as one sparse Newton system with pymrm:
first-order upwind convection with unit velocity in the space time ``w = W/m_dot`` (kg_cat s per
kg_gas), the pointwise source by ``NumJac``. The converged Jacobian gives the sensitivity of the
whole profile to the cooling ``h``, the coolant temperature ``T_c`` or the inlet ethene fraction
``y_E`` by one sparse solve, which is what Morbidelli and Varma's generalized criterion needs.

Units: amounts ``n_j`` in mol per kg of gas, ``P`` in bar, ``T`` in K, rates in mol/(kg_cat s),
heat release ``q`` in W/kg_cat, ``h`` in W/(kg_cat K), ``c_p`` in J/(kg_gas K). Field order of a
state ``u[..., k]``: E, O, C, EO, W, N (ethene, oxygen, carbon dioxide, ethene oxide, water,
nitrogen), then T. No plotting, no data loading.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.sparse.linalg import splu
from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind, construct_div,
                   continuation, newton)
from pymrm.checks import residual_check

SPECIES = ("E", "O", "C", "EO", "W", "N")
R_GAS = 8.314462618              # J/(mol K)
T_REF = 476.0                    # K, temperature of the printed adsorption constants (p. 54)
M_ETHENE = 28.0532e-3            # kg/mol
M_AIR = 28.9647e-3               # kg/mol, dry air
Y_O2_AIR = 0.20946               # O2 mole fraction of dry air


def table6_model3(table6):
    """The model III constants of Table 6 as floats.

    Parameters
    ----------
    table6 : pandas.DataFrame
        Columns ``model``, ``reaction``, ``parameter``, ``value`` (C2.11's CSV).

    Returns
    -------
    dict
        ``{"r1": p, "r2": p}``; each ``p`` has ``k_inf``, ``T_act``, ``K`` (``{j: K_j at
        476 K}``) and ``pair`` (``{j: (K_inf, T_ads)}`` where Table 6 prints a temperature
        dependence).
    """
    out = {}
    for rx in ("r1", "r2"):
        rows = table6[(table6.model == "III") & (table6.reaction == rx)]
        v = dict(zip(rows.parameter, rows.value.astype(float)))
        out[rx] = {"k_inf": v["k_inf"], "T_act": v["T_act"],
                   "K": {j: v[f"K_{j}"] for j in SPECIES[:5] if f"K_{j}" in v},
                   "pair": {j: (v[f"K_{j}_inf"], v[f"T_ads_{j}"]) for j in SPECIES[:5]
                            if f"K_{j}_inf" in v}}
    return out


def adsorption(p, T, inhibition=True):
    """``{j: K_j(T)}`` in 1/bar: ``K_inf exp(T_ads/T)`` where a pair is printed, the 476 K
    value otherwise (p. 54). ``inhibition=False`` sets the product constants to zero."""
    K = dict(p["K"])
    for j, (k_inf, t_ads) in p["pair"].items():
        K[j] = k_inf * np.exp(t_ads / T)
    if not inhibition:
        K.update(C=0.0, EO=0.0, W=0.0)
    return K


def rate(p, T, P, inhibition=True):
    """Model III, ``k_r K_E K_O^(1/2) P_E P_O^(1/2) / D^2`` with ``k_r = k_inf exp(-T_act/T)``
    and ``D = 1 + K_E P_E + (K_O P_O)^(1/2) + K_C P_C + K_EO P_EO + K_W P_W``."""
    K = adsorption(p, T, inhibition)
    D = (1.0 + K["E"] * P[..., 0] + np.sqrt(K["O"] * P[..., 1]) + K["C"] * P[..., 2]
         + K["EO"] * P[..., 3] + K["W"] * P[..., 4])
    return (p["k_inf"] * np.exp(-p["T_act"] / T) * K["E"] * np.sqrt(K["O"])
            * P[..., 0] * np.sqrt(P[..., 1]) / D ** 2)


def apparent_activation(p, T, P):
    """``-R dln r/d(1/T)`` in J/mol at fixed composition:
    ``R [T_act - T_ads,E + 2 sum_j T_ads,j K_j P_j / D]`` over the printed pairs."""
    K = adsorption(p, T)
    D = (1.0 + K["E"] * P[..., 0] + np.sqrt(K["O"] * P[..., 1]) + K["C"] * P[..., 2]
         + K["EO"] * P[..., 3] + K["W"] * P[..., 4])
    s = sum(t_ads * K[j] * P[..., SPECIES.index(j)] for j, (_, t_ads) in p["pair"].items())
    return R_GAS * (p["T_act"] - p["pair"]["E"][1] + 2.0 * s / D)


def inlet(y_E):
    """Inlet amounts ``n_j0`` (mol per kg of gas) for ethene fraction ``y_E`` in dry air."""
    M = y_E * M_ETHENE + (1.0 - y_E) * M_AIR
    return np.array([y_E, (1.0 - y_E) * Y_O2_AIR, 0.0, 0.0, 0.0,
                     (1.0 - y_E) * (1.0 - Y_O2_AIR)]) / M


def dinlet_dy(y_E):
    """``d n_j0 / d y_E`` (analytic)."""
    M = y_E * M_ETHENE + (1.0 - y_E) * M_AIR
    v = M * inlet(y_E)
    dv = np.array([1.0, -Y_O2_AIR, 0.0, 0.0, 0.0, -(1.0 - Y_O2_AIR)])
    return (dv * M - v * (M_ETHENE - M_AIR)) / M ** 2


class Kinetics:
    """Rates, species sources and heat release of the two parallel reactions.

    ``C2H4 + 1/2 O2 -> C2H4O`` (-105 kJ/mol) and ``C2H4 + 3 O2 -> 2 CO2 + 2 H2O``
    (-1323 kJ/mol ethene, p. 49). ``r2`` is the CO2 production rate (eq. 2), so the
    combustion burns ``burnt_per_co2 * r2`` ethene with ``burnt_per_co2 = 1/2``.

    Parameters
    ----------
    par : dict
        From :func:`table6_model3`.
    F : float
        Factor on the r2 rate constant (1 = Table 6 as printed).
    a : float
        Catalyst activity on both rates.
    P_total : float
        Total pressure (bar).
    dH : tuple
        ``(-dH1, -dH2)`` in J/mol ethene.
    burnt_per_co2, o2_per_burnt : float
        Ethene burnt per CO2 formed (1/2) and O2 per ethene burnt (3): break-row options.
    inhibition : bool
        Product inhibition on (Table 6) or off (a break-row option).
    factor_on : str
        ``"r2"`` (as C2.11) or ``"r1"`` (a break-row option).
    """

    # codehealth: allow params -- rate law, factor, activity and the stoichiometric break readings
    def __init__(self, par, F=1.0, a=1.0, P_total=5.0, dH=(105e3, 1323e3), burnt_per_co2=0.5,
                 o2_per_burnt=3.0, inhibition=True, factor_on="r2"):
        self.par, self.F, self.a, self.P_total, self.dH = par, F, a, P_total, dH
        self.nu, self.o2, self.inhibition, self.factor_on = (burnt_per_co2, o2_per_burnt,
                                                             inhibition, factor_on)

    def pressures(self, n):
        """Partial pressures ``[E, O, C, EO, W]`` (bar) of amounts ``n[..., 6]``."""
        return self.P_total * n[..., :5] / np.sum(n, axis=-1, keepdims=True)

    def rates(self, T, n):
        """``(r1, r2)`` in mol/(kg_cat s), activity and factor applied."""
        P = self.pressures(n)
        r1 = self.a * rate(self.par["r1"], T, P, self.inhibition)
        r2 = self.a * rate(self.par["r2"], T, P, self.inhibition)
        if self.factor_on == "r2":
            return r1, self.F * r2
        return self.F * r1, r2

    def sources(self, T, n):
        """``(dn/dw [..., 6], q)``: species sources and heat release (W/kg_cat)."""
        r1, r2 = self.rates(T, n)
        burnt = self.nu * r2
        dn = np.stack([-(r1 + burnt), -(0.5 * r1 + self.o2 * burnt), r2, r1, r2,
                       np.zeros_like(r1)], axis=-1)
        return dn, self.dH[0] * r1 + self.dH[1] * burnt

    def heat(self, T, n):
        """Heat release ``q`` (W/kg_cat)."""
        return self.sources(T, n)[1]

    def dheat_dT(self, T, n):
        """``dq/dT`` at fixed composition, from the apparent activation energies (E10)."""
        P = self.pressures(n)
        r1, r2 = self.rates(T, n)
        e1 = apparent_activation(self.par["r1"], T, P)
        e2 = apparent_activation(self.par["r2"], T, P)
        return (self.dH[0] * r1 * e1 + self.dH[1] * self.nu * r2 * e2) / (R_GAS * T ** 2)

    def selectivity(self, T, n):
        """Differential selectivity ``r1/(r1 + r2/2)``."""
        r1, r2 = self.rates(T, n)
        return r1 / (r1 + self.nu * r2)


class CooledTube:
    """Steady wall-cooled plug-flow tube, pseudo-homogeneous, isobaric (E3).

    ``dn/dw = a S(n, T)``, ``c_p dT/dw = q(n, T) - h (T - T_c)`` with ``T(0) = T_c``, written
    as a steady convection problem with unit velocity in ``w``: first-order upwind fluxes on
    a uniform grid, ``construct_div``, the source by ``NumJac`` per cell, ``pymrm.newton``.
    First-order upwind with ``v > 0`` is implicit Euler in ``w``: the value of cell ``i`` is
    the value at its downstream face, so the profile on the faces is ``[inlet, cells]``.

    Parameters
    ----------
    kin : Kinetics
        Rate law.
    y_E, T_c, h, c_p : float
        Inlet ethene fraction, coolant (and inlet) temperature, cooling coefficient
        ``4 U/(d_t rho_b)`` and gas heat capacity.
    w_end : float
        Tube length in space time (kg_cat s/kg_gas).
    n : int
        Number of cells.
    mode : str
        ``"full"``; ``"frozen"`` (rates held at the inlet composition, no species change:
        the Semenov limit); ``"isothermal"`` (no energy balance, ``T = T_c``).
    """

    #: ``newton`` lower bounds: ethene and oxygen amounts (square roots and the rate) and T.
    LOWER = np.array([0.0, 0.0, -np.inf, -np.inf, -np.inf, -np.inf, 0.0])

    # codehealth: allow params -- one operating point of the tube and its grid
    def __init__(self, kin, y_E, T_c, h, c_p, w_end, n, mode="full", eps_jac=1e-8):
        self.kin, self.h, self.c_p, self.mode = kin, float(h), float(c_p), mode
        self.eps_jac = eps_jac
        self.shape = (int(n), 7)
        self.w_f = np.linspace(0.0, float(w_end), self.shape[0] + 1)
        self.w_c = 0.5 * (self.w_f[:-1] + self.w_f[1:])
        self._build_operators()
        self.set_inlet(y_E, T_c)
        self.u = np.tile(self.u_in, (self.shape[0], 1))

    def _build_operators(self):
        # w = 0, outward normal -w: u = u_in            -> a = 0, b = 1, d = 1 x u_in (shapes_d)
        # w = w_end: first-order convective outlet       -> {"outflow": True}
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"outflow": True})
        conv, conv_in, _ = construct_convflux_upwind(self.shape, self.w_f, self.w_c, bc, v=1.0,
                                                     axis=0, shapes_d=((1, 7), None))
        # v = 1: w is the space time, so d(v u)/dw = source. nu = 0: Cartesian axis.
        div = construct_div(self.shape, self.w_f, nu=0, axis=0)
        self.jac_const = div @ conv
        self.g_inlet = div @ conv_in
        self.numjac = NumJac(self.shape, eps_jac=self.eps_jac)   # pointwise source: last axis in full

    def set_inlet(self, y_E=None, T_c=None):
        """Set the inlet (and coolant) state; the operators are not rebuilt."""
        self.y_E = self.y_E if y_E is None else float(y_E)
        self.T_c = self.T_c if T_c is None else float(T_c)
        self.n_in = inlet(self.y_E)
        self.u_in = np.append(self.n_in, self.T_c)
        self.g_const = (self.g_inlet @ self.u_in).reshape((-1, 1))
        self.scale = np.append(np.full(6, self.n_in[0]), self.T_c)

    def source(self, u):
        """Pointwise source ``(dn/dw, dT/dw)`` of a state ``u[..., 7]``."""
        n, T = u[..., :6], u[..., 6]
        if self.mode == "frozen":
            q = self.kin.heat(T, np.broadcast_to(self.n_in, n.shape))
            dn = np.zeros_like(n)
        else:
            dn, q = self.kin.sources(T, n)
        dT = np.zeros_like(T) if self.mode == "isothermal" else (
            q - self.h * (T - self.T_c)) / self.c_p
        return np.concatenate([dn, dT[..., None]], axis=-1)

    def residual(self, u):
        """Steady residual (convection minus source) and its Jacobian."""
        u = u.reshape(self.shape)
        g_src, jac_src = self.numjac(self.source, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_src.reshape((-1, 1))
        return g, self.jac_const - jac_src

    def solve(self, u0=None, tol=1e-11, maxfev=60):
        """Newton from ``u0`` (default: the current state); raises unless it converges.

        The step test is per field on ``scale`` (``n_E0`` for the amounts, ``T_c`` for T,
        pitfall P7); the bounds keep the ethene and oxygen amounts and T non-negative.
        """
        u0 = self.u if u0 is None else u0
        res = newton(self.residual, np.asarray(u0, dtype=float).reshape(self.shape), tol=tol,
                     maxfev=maxfev, scale=self.scale, stagnation=True,
                     bounds=(self.LOWER, None))
        if not res.success:
            raise RuntimeError(f"CooledTube.solve (y_E {self.y_E:g}, T_c {self.T_c:g}, "
                               f"h {self.h:g}): {res.message}")
        self.u, self.result = res.x.reshape(self.shape), res
        return self

    def check(self, tol=1e-9):
        """``pymrm.checks.residual_check`` on the converged state."""
        return residual_check(self.residual, self.u, tol=tol)

    def trace(self, name, values, max_change=0.2, tol=1e-11):
        """Follow the solution over values of ``name`` (``"h"``, ``"y_E"`` or ``"T_c"``)
        with ``pymrm.continuation``, from the current state; returns the result and leaves
        the tube at the last value reached. Newton follows the branch it starts on: check
        traced states against fresh solves before relying on them."""
        values = np.asarray(values, dtype=float)

        def residual(u, p):
            self._set(name, p)
            return self.residual(u)

        self._set(name, values[0])
        start = self.solve(tol=tol).u
        with warnings.catch_warnings():
            # a secant prediction past the bounds is projected onto them by newton, which
            # warns; the corrector then converges from the projected guess or is rejected
            warnings.filterwarnings("ignore", message="newton: .* projected onto them")
            res = continuation(residual, start, values, max_change=max_change,
                               newton_kwargs=dict(tol=tol, maxfev=60, scale=self.scale,
                                                  stagnation=True,
                                                  bounds=(self.LOWER, None)))
        if not res.success:
            raise RuntimeError(f"CooledTube.trace in {name}: {res.message}")
        self._set(name, res.last_param if res.last_param is not None else values[0])
        self.u = np.asarray(res.x, dtype=float).reshape(self.shape)
        return res

    def move(self, name, value, steps=4):
        """Continue the current state to ``name = value`` in ``steps`` equal steps."""
        start = {"h": self.h, "y_E": self.y_E, "T_c": self.T_c}[name]
        if start == value:
            return self.solve()
        self.trace(name, np.linspace(start, value, steps + 1))
        return self

    def _set(self, name, value):
        if name == "h":
            self.h = float(value)
        elif name == "y_E":
            self.set_inlet(y_E=value)
        elif name == "T_c":
            self.set_inlet(T_c=value)
        else:
            raise ValueError(f"unknown parameter {name!r}")

    # ---- profiles on the faces ------------------------------------------------------------

    @property
    def faces(self):
        """The state on the faces ``w_f``: the inlet, then the cells (downstream faces)."""
        return np.vstack([self.u_in, self.u])

    def outlet(self):
        """Outlet face state by ``compute_boundary_values`` under the outlet bc."""
        bc = ({"a": 0.0, "b": 1.0, "d": self.u_in[None, :]}, {"outflow": True})
        return compute_boundary_values(self.u, self.w_f, self.w_c, bc)[2].ravel()

    def invariants(self):
        """C, H, O invariants (E9) on the faces, and ``c_p T - dH1 n_EO - dH2 nu n_C``."""
        f = self.faces
        C = 2 * f[:, 0] + f[:, 2] + 2 * f[:, 3]
        H = 4 * f[:, 0] + 4 * f[:, 3] + 2 * f[:, 4]
        O = 2 * f[:, 1] + 2 * f[:, 2] + f[:, 3] + f[:, 4]
        E = (self.c_p * f[:, 6] - self.kin.dH[0] * f[:, 3]
             - self.kin.dH[1] * self.kin.nu * f[:, 2])
        return {"C": C, "H": H, "O": O, "energy": E}

    def hot_spot(self, dT=None):
        """Hot spot ``(w*, T*, dT*/dphi)`` from the parabola through the three faces around
        the largest face temperature; ``dT`` is a face sensitivity ``dT_f/dphi`` (optional).

        Raises when the largest face temperature is at either end of the tube."""
        T = self.faces[:, 6]
        k = int(np.argmax(T))
        if k == 0 or k == T.size - 1:
            raise RuntimeError(f"hot spot at a tube end (face {k} of {T.size})")
        a, b, c = T[k - 1], T[k], T[k + 1]
        den, d = a - 2.0 * b + c, a - c
        x = 0.5 * d / den
        dw = self.w_f[1] - self.w_f[0]
        Ts = b - 0.25 * d * x
        if dT is None:
            return self.w_f[k] + x * dw, Ts, None
        wa = -d / (4.0 * den) + d * d / (8.0 * den * den)
        wc = d / (4.0 * den) + d * d / (8.0 * den * den)
        wb = 1.0 - d * d / (4.0 * den * den)
        return self.w_f[k] + x * dw, Ts, wa * dT[k - 1] + wb * dT[k] + wc * dT[k + 1]

    def sensitivity(self, name):
        """Face sensitivity ``dT_f/dphi`` for ``phi`` = ``"h"``, ``"T_c"`` or ``"y_E"``:
        ``J du/dphi = -dg/dphi`` on the Jacobian at the converged state (E13)."""
        g, jac = self.residual(self.u)
        n = self.shape[0]
        dg = np.zeros(self.shape)
        e = np.zeros(7)
        if name == "h":
            dg[:, 6] = (self.u[:, 6] - self.T_c) / self.c_p
        elif name == "T_c":
            dg[:, 6] = -self.h / self.c_p if self.mode != "isothermal" else 0.0
            e[6] = 1.0
        elif name == "y_E":
            e[:6] = dinlet_dy(self.y_E)
        else:
            raise ValueError(f"unknown parameter {name!r}")
        rhs = -(dg.reshape((-1, 1)) + self.g_inlet @ e.reshape((-1, 1)))
        du = splu(jac).solve(rhs).reshape((n, 7))
        return np.concatenate([[e[6]], du[:, 6]])

    def mv(self, name, value):
        """Morbidelli-Varma normalised sensitivity ``S* = (phi/T*) dT*/dphi`` (eq. 34b, T*
        absolute) and the rise-normalised variant ``d ln(T* - T_c)/d ln phi``."""
        _, Ts, dTs = self.hot_spot(self.sensitivity(name))
        drise = dTs - (1.0 if name == "T_c" else 0.0)
        return value / Ts * dTs, value / (Ts - self.T_c) * drise

    def at_conversion(self, X=0.5):
        """``(w, state)`` where the ethene conversion crosses ``X``: the root of the quadratic
        through the three faces around the crossing; raises unless exactly one face interval
        holds a crossing."""
        from pymrm.checks import find_roots
        f = self.faces
        conv = 1.0 - f[:, 0] / self.n_in[0] - X
        k = np.nonzero(np.sign(conv[:-1]) != np.sign(conv[1:]))[0]
        if k.size != 1:
            raise RuntimeError(f"conversion {X} crossed in {k.size} face intervals")
        k = int(k[0])
        kk = min(max(k, 1), f.shape[0] - 2)            # the three faces around interval k
        w3 = self.w_f[kk - 1:kk + 2]

        def interp(w, col):
            L = [np.prod([(w - w3[j]) / (w3[i] - w3[j]) for j in range(3) if j != i])
                 for i in range(3)]
            return sum(L[i] * col[kk - 1 + i] for i in range(3))

        lo, hi = self.w_f[k], self.w_f[k + 1]
        r = find_roots(lambda w: interp(w, conv), lo, hi, n_scan=8, xtol=1e-13 * hi)
        if len(r["roots"]) != 1:
            raise RuntimeError(f"conversion root: {r['message']}")
        w = r["roots"][0]
        return w, np.array([interp(w, f[:, j]) for j in range(7)])


# ---- Morbidelli-Varma critical points ----------------------------------------------------------

# codehealth: allow params -- the march, its direction and its step control
def scan_sensitivity(tube, name, x_start, x_stop, dT=8.0, max_dlog=0.05, past=3,
                     min_peak=1.0, which=0):
    """March ``name`` from ``x_start`` toward ``x_stop`` by continuation and record the
    normalised sensitivity ``S*`` at every point.

    The step in ``ln x`` is ``min(max_dlog, dT/(T* |S*|))``, so the hot spot moves by about
    ``dT`` per step and a sharp extremum of ``|S*|`` is sampled by many points. The march
    stops at ``x_stop``, or ``past`` points beyond the largest ``|S*|`` once ``|S*|`` has
    fallen below half of it (and the largest exceeds ``min_peak``), or when continuation
    cannot follow the branch (a fold of the discrete system on the runaway side).
    ``which`` selects the criterion: 0 for eq. 34b (absolute ``T*``), 1 for the
    rise-normalised variant.

    Returns
    -------
    dict
        ``x``, ``S``, ``S_rise``, ``T``, ``w`` (arrays), ``states`` (list), ``stopped``
        (the reason).
    """
    sgn = 1.0 if x_stop > x_start else -1.0
    tube.move(name, x_start, steps=4)
    out = {k: [] for k in ("x", "S", "S_rise", "T", "w", "states")}
    x, stopped = float(x_start), "end of range"
    while True:
        w, Ts, _ = tube.hot_spot()
        s, sr = tube.mv(name, x)
        for k, v in zip(out, (x, s, sr, Ts, w, tube.u.copy())):
            out[k].append(v)
        a = np.abs(np.array(out["S_rise"] if which else out["S"]))
        j = int(np.argmax(a))
        if a[j] > min_peak and len(a) - 1 - j >= past and a[-1] < 0.5 * a[j]:
            stopped = "past the extremum"
            break
        if sgn * (x - x_stop) >= 0:
            break
        step = min(max_dlog, dT / (Ts * abs(s) + 1e-300))
        x_new = x * np.exp(sgn * step)
        x_new = min(x_new, x_stop) if sgn > 0 else max(x_new, x_stop)
        try:
            tube.move(name, x_new, steps=1)
        except RuntimeError as err:
            stopped = f"branch lost: {err}"
            break
        x = x_new
    res = {k: (np.array(v) if k != "states" else v) for k, v in out.items()}
    res["stopped"] = stopped
    return res


def _set_param(tube, name, x):
    tube._set(name, x)


def refine_extremum(tube, name, scan, xatol=3e-7, which=0):
    """Bounded maximisation of ``|S*|`` between the scan neighbours of its largest value.

    ``xatol`` is relative to the critical value. Near a flat extremum the function changes
    by about ``(xatol/width)^2``; keeping that far above the round-off noise of ``S*`` keeps
    the maximiser's steps, and so the hot spot at the critical point, the same on every
    BLAS kernel (a much smaller ``xatol`` lets it wander in that noise).

    Each evaluation starts from the scan state on the start side of ``x`` (the side the
    march came from) and continues to ``x`` in two steps. Raises when the largest ``|S*|``
    of the scan is at either end of the scan.

    Returns
    -------
    dict
        ``x`` (critical value), ``S``, ``S_rise``, ``T`` (hot spot), ``w``.
    """
    from scipy.optimize import minimize_scalar
    xs, a = scan["x"], np.abs(scan["S_rise"] if which else scan["S"])
    j = int(np.argmax(a))
    if j == 0 or j == xs.size - 1:
        raise RuntimeError(f"largest |S*| at the end of the scan ({name} = {xs[j]:g}; "
                           f"{scan['stopped']})")
    sgn = np.sign(xs[-1] - xs[0])

    def solve_at(x):
        i = int(np.nonzero(sgn * (xs - x) <= 0)[0][-1])     # last scan point before x
        _set_param(tube, name, xs[i])
        tube.u = scan["states"][i].copy()
        tube.move(name, x, steps=2)
        return tube

    def f(x):
        return -abs(solve_at(x).mv(name, x)[which])

    lo, hi = sorted((xs[j - 1], xs[j + 1]))
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded",
                        options=dict(xatol=xatol * abs(xs[j])))
    tube = solve_at(r.x)
    w, Ts, _ = tube.hot_spot()
    s, sr = tube.mv(name, r.x)
    if not lo < r.x < hi:
        raise RuntimeError("the refined extremum sits on its bracket")
    return {"x": float(r.x), "S": float(s), "S_rise": float(sr), "T": float(Ts), "w": float(w)}


def refined_grid(tube, factor=2):
    """A copy of ``tube`` on ``factor`` times as many cells, its state repeated per cell."""
    new = CooledTube(tube.kin, tube.y_E, tube.T_c, tube.h, tube.c_p, tube.w_f[-1],
                     factor * tube.shape[0], tube.mode, tube.eps_jac)
    new.u = np.repeat(tube.u, factor, axis=0)
    return new


# codehealth: allow params -- the operating point, the parameter scanned and the grid ladder
def critical_point(kin, y_E, T_c, h, c_p, name, x_start, x_stop, levels=(500, 1000),
                   w_end=600.0, restart=0.02, which=0, **scan_kw):
    """Morbidelli-Varma critical value of ``name`` on each grid of ``levels``.

    The coarsest grid scans from ``x_start`` (the stable side) toward ``x_stop`` and refines
    the largest ``|S*|``; each finer grid (twice the cells, state carried over) scans again
    from ``exp(restart)`` before the previous critical value on the stable side. Raises when
    the hot spot at the critical point lies in the last fifth of the tube.

    Returns
    -------
    list of dict
        One :func:`refine_extremum` result per level, with ``n`` and the coarse ``scan``.
    """
    base = {"h": h, "y_E": y_E, "T_c": T_c}
    base[name] = x_start
    tube = CooledTube(kin, base["y_E"], base["T_c"], base["h"], c_p, w_end, levels[0])
    tube.trace("h", np.geomspace(4.0 * base["h"] + 1.0, base["h"], 5))   # from strong cooling
    sgn = 1.0 if x_stop > x_start else -1.0
    out, scan = [], scan_sensitivity(tube, name, x_start, x_stop, which=which, **scan_kw)
    for k, n in enumerate(levels):
        if k:
            x0 = out[-1]["x"] * np.exp(-sgn * restart)
            before = np.nonzero(sgn * (scan["x"] - x0) <= 0)[0]
            i = int(before[-1]) if before.size else 0
            x0 = x0 if before.size else float(scan["x"][0])
            _set_param(tube, name, scan["x"][i])
            tube.u = scan["states"][i].copy()
            tube = refined_grid(tube, n // tube.shape[0])
            scan = scan_sensitivity(tube, name, x0, x_stop, which=which, **scan_kw)
        r = refine_extremum(tube, name, scan, which=which)
        if r["w"] > 0.8 * w_end:
            raise RuntimeError(f"hot spot at w = {r['w']:.4g} of {w_end:g}: lengthen the tube")
        r["n"] = n
        if k == 0:
            r["scan"] = scan
        out.append(r)
    return out


def richardson(values, order=1.0, ratio=2.0):
    """First Richardson step on the two finest values of a ladder."""
    a, b = float(values[-2]), float(values[-1])
    return b + (b - a) / (ratio ** order - 1.0)
