"""Xu and Froment (1989) Part II: the slab pellet and the industrial reformer tube.

The intrinsic rates are Part I's (page C2.1, passed in as its ``XuFroment``
object); this module adds the effective diffusivities (Fuller binaries, Wilke
mixture rule, Knudsen in each pore class), the slab pellet of eqs. 13 solved with
pymrm operators, and the plug-flow tube of eqs. 10 marched in z with a nested
Newton solve per implicit-Euler step. No data loading and no plotting here: the
notebook passes every table in.

Units: SI, except partial pressures in bar and rates in kmol/(kg_cat h) as
Part I prints them. Species order CH4, H2O, H2, CO, CO2, N2.
"""
import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy import sparse
import scipy.sparse.linalg  # noqa: F401  (spsolve in thiele_slab_eta)
from scipy.optimize import brentq
from pymrm import NumJac, compute_boundary_values, construct_div, construct_grad, newton

R = 8.314            # J/(mol K), the value C2.1 uses (R_KJ = 8.314e-3)
BAR = 1.0e5          # Pa per bar
KMOLH = 1000.0 / 3600.0   # kmol/h to mol/s
P_ATM = 1.01325      # bar per atm (Fuller's correlation takes atm)
SPECIES = ("CH4", "H2O", "H2", "CO", "CO2", "N2")
# change of moles per mol CH4 fed per unit (x, y): x total CH4 conversion, y conversion to CO2
ST = np.array([[-1.0, 0.0], [-1.0, -1.0], [3.0, 1.0], [1.0, -1.0], [0.0, 1.0], [0.0, 0.0]])
ATOMS = np.array([[1, 0, 0, 1, 1, 0], [4, 2, 2, 0, 0, 0], [0, 1, 0, 1, 2, 0]], dtype=float)


def _dense(a):
    return np.asarray(a.toarray() if sparse.issparse(a) else a).ravel()


@dataclass
class Case:
    """Inputs of the reformer: feed, geometry, catalyst, kinetics and property tables.

    Parameters
    ----------
    kin : object
        C2.1's ``XuFroment`` (rates through ``constants`` and ``equilibrium``).
    shomate : ndarray
        ``(6, 6)`` Shomate rows A..F (kJ/mol with formation) in species order.
    mw, vol : ndarray
        Molar masses (g/mol) and Fuller diffusion volumes of the six species.
    pore_d, pore_dv : sequence
        Pore class diameters (Angstrom) and volumes (cm3/g), Table 1.
    n0 : ndarray
        Moles per mol of equivalent CH4 fed (Table 2 ratios).
    """

    kin: object
    shomate: np.ndarray
    mw: np.ndarray
    vol: np.ndarray
    pore_d: tuple
    pore_dv: tuple
    n0: np.ndarray = field(default_factory=lambda: np.array([1.0, 3.358, 0.122, 0.0, 0.056, 0.164]))
    F0: float = 5.168 * KMOLH          # mol/s equivalent CH4 per tube
    T0: float = 793.15
    d_ti: float = 0.1016
    d_te: float = 0.1322
    L_heat: float = 11.12
    L_tot: float = 12.0
    rho_s: float = 2355.2
    tau: float = 3.54
    l_act: float = 0.002
    ring_factor: float = 1.071
    eps_B: float = 0.50
    eps_s: float = None                # None: V_g rho_s (rho_s read as pellet density)
    knudsen: bool = True
    act: float = 1.0
    K_H2O_T: float = None              # break-row knob: adsorption constant of H2O at a fixed T

    def __post_init__(self):
        if self.eps_s is None:
            self.eps_s = sum(self.pore_dv) * 1e-3 * self.rho_s

    def with_(self, **kw):
        """A copy with some fields replaced."""
        return replace(self, **kw)

    @property
    def omega(self):
        return math.pi * self.d_ti ** 2 / 4.0

    @property
    def rho_B(self):
        return self.rho_s * (1.0 - self.eps_B)

    # ---- composition (E1)
    def moles(self, x, y):
        return self.n0 + ST @ np.array([x, y])

    def state(self, x, y, pt):
        """Partial pressures (bar, five reacting species) and the six mole fractions."""
        n = self.moles(x, y)
        yf = n / n.sum()
        return pt * yf[:5], yf

    # ---- kinetics (E2, E3), partial-pressure form of C2.1's rates
    def rates(self, p, T, act=None):
        """``(r1, r2, r3)`` in kmol/(kg h) at partial pressures ``p`` (bar)."""
        pCH4, pH2O, pH2, pCO, pCO2 = p
        k = self.kin.constants(T)
        if self.K_H2O_T is not None:
            k["K_H2O"] = self.kin.constants(self.K_H2O_T)["K_H2O"]
        a = self.act if act is None else act
        K1, K2, K3 = self.kin.equilibrium(T)
        den = (1 + k["K_CO"] * pCO + k["K_H2"] * pH2 + k["K_CH4"] * pCH4
               + k["K_H2O"] * pH2O / pH2) ** 2
        r1 = a * k["k1"] / pH2 ** 2.5 * (pCH4 * pH2O - pH2 ** 3 * pCO / K1) / den
        r2 = a * k["k2"] / pH2 * (pCO * pH2O - pH2 * pCO2 / K2) / den
        r3 = a * k["k3"] / pH2 ** 3.5 * (pCH4 * pH2O ** 2 - pH2 ** 4 * pCO2 / K3) / den
        return r1, r2, r3

    def approach(self, p, T):
        """``V_i = Q_i / K_i`` of reactions I, II, III (eq. 14)."""
        pCH4, pH2O, pH2, pCO, pCO2 = p
        K1, K2, K3 = self.kin.equilibrium(T)
        return (pH2 ** 3 * pCO / (pCH4 * pH2O * K1), pH2 * pCO2 / (pCO * pH2O * K2),
                pH2 ** 4 * pCO2 / (pCH4 * pH2O ** 2 * K3))

    # ---- thermo (E12)
    def h_i(self, T):
        """Molar enthalpies (J/mol, formation included) of the six species."""
        t = T / 1000.0
        A, B, C, D, E, F = self.shomate.T
        return 1e3 * (A * t + B * t ** 2 / 2 + C * t ** 3 / 3 + D * t ** 4 / 4 - E / t + F)

    def cp_i(self, T):
        t = T / 1000.0
        A, B, C, D, E, _ = self.shomate.T
        return A + B * t + C * t ** 2 + D * t ** 3 + E / t ** 2

    def dH(self, T):
        """Reaction enthalpies I and II (J/mol) from the same h_i."""
        h = self.h_i(T)
        return h @ ST[:, 0], h @ ST[:, 1]

    def h_flow(self, x, y, T):
        """Enthalpy flow (W) of the gas at conversions x, y and temperature T."""
        return self.F0 * self.moles(x, y) @ self.h_i(T)

    # ---- effective diffusivities (E5-E7)
    def d_eff(self, yf, T, pt, pore_d=None):
        """Effective diffusivities (m2/s) of CH4, H2O, H2, CO, CO2 at bulk mole fractions."""
        pore_d = self.pore_d if pore_d is None else pore_d
        mw, vol = self.mw, self.vol
        Dij = 1e-7 * T ** 1.75 * np.sqrt(1 / mw[:, None] + 1 / mw[None, :]) / (
            pt / P_ATM * (vol[:, None] ** (1 / 3) + vol[None, :] ** (1 / 3)) ** 2)
        S = np.asarray(self.pore_dv) / sum(self.pore_dv)
        De = np.empty(5)
        for i in range(5):
            Dm = (1 - yf[i]) / sum(yf[j] / Dij[i, j] for j in range(6) if j != i)
            if self.knudsen:
                DK = (np.asarray(pore_d) * 1e-10 / 3.0
                      * math.sqrt(8 * R * T / (math.pi * mw[i] * 1e-3)))
                De[i] = S @ (1.0 / (1.0 / Dm + 1.0 / DK))
            else:
                De[i] = Dm
        return self.eps_s / self.tau * De

    # ---- equilibrium tube (E14)
    def equilibrium_state(self, H, pt):
        """Bulk equilibrium of I and II at enthalpy flow ``H`` (W) and ``pt`` (bar): (x, y, T)."""
        n0 = self.n0

        def at_T(T):
            K1, K2, _ = self.kin.equilibrium(T)

            def y_of_x(x):
                f = lambda y: (n0[4] + y) * (n0[2] + 3 * x + y) - K2 * (x - y) * (n0[1] - x - y)
                return brentq(f, -n0[4] + 1e-14, x - 1e-14, xtol=1e-15, rtol=1e-15)

            def g(x):
                p = self.state(x, y_of_x(x), pt)[0]
                return math.log(p[3] * p[2] ** 3 / (p[0] * p[1])) - math.log(K1)
            x = brentq(g, 1e-9, 1 - 1e-9, xtol=1e-15, rtol=1e-15)
            return x, y_of_x(x)
        T = brentq(lambda T: self.h_flow(*at_T(T), T) - H, 700.0, 1300.0, xtol=1e-12, rtol=1e-15)
        return (*at_T(T), T)


class Pellet:
    """Eqs. 13 in one active layer: p_CH4 and p_CO2 by finite volumes, the rest by eqs. 6e-g.

    Parameters
    ----------
    case : Case
    n : int
        Cells across the layer, graded to the surface by ``depth = l expm1(beta s) / expm1(beta)``
        with ``s`` uniform, so that doubling ``n`` nests the grids (Richardson).
    geom : {"slab", "outer", "inner"}
        Planar layer (the paper), or the outer or inner lateral layer of the ring as a
        cylindrical shell.
    cut : float, optional
        Depth (m) beyond which the rate is set to zero (the paper's innermost spline).
    """

    # codehealth: allow params -- one pellet: its grid, geometry, cut and ring dimensions
    def __init__(self, case, n, geom="slab", cut=None, beta=6.0, d_pe=0.0173, d_pi=0.0084):
        self.case, self.n, self.geom, self.cut = case, n, geom, cut
        L = case.l_act
        if cut is not None:              # a face exactly at the cut, at s = 1/2 on every level
            beta = 2.0 * math.log(L / cut - 1.0)
        s = np.linspace(0.0, 1.0, n + 1)
        depth = L * np.expm1(beta * s) / np.expm1(beta)
        if geom == "slab":
            x_f, nu, left = depth, 0, True
        elif geom == "outer":
            x_f, nu, left = d_pe / 2 - depth[::-1], 1, False
        else:
            x_f, nu, left = d_pi / 2 + depth, 1, True
        self.x_f, self.x_c, self.left = x_f, 0.5 * (x_f[1:] + x_f[:-1]), left
        dep_c = 0.5 * (depth[1:] + depth[:-1])
        if not left:
            dep_c = dep_c[::-1]
        self.depth_c = dep_c
        self.vol_c = np.diff(x_f) if nu == 0 else np.diff(x_f ** 2) / 2
        self.active = np.ones(n) if cut is None else (dep_c < cut).astype(float)
        dirichlet = {"a": 0.0, "b": 1.0, "d": 1.0}       # p = p_bulk at the surface
        noflux = {"a": 1.0, "b": 0.0, "d": 0.0}          # dp/dn = 0 at the inner face
        self.bc_unit = (dirichlet, noflux) if left else (noflux, dirichlet)
        shape = (n, 2)                                    # cells first, fields (CH4, CO2) last
        G, Gb = construct_grad(shape, x_f, self.x_c, self.bc_unit, axis=0)
        Dv = construct_div(shape, x_f, nu=nu, axis=0)    # nu = 0 planar, 1 cylindrical shell
        self.lap = (Dv @ G).tocsc()
        self.lap_b = _dense(Dv @ Gb).reshape(n, 2)
        self.numjac = NumJac(shape)

    def dependent(self, pc, pv, pb, De):
        """Eqs. 6e-6g: the five partial pressures from p_CH4 (pc) and p_CO2 (pv)."""
        dC, dV = pc - pb[0], pv - pb[4]
        pH2O = pb[1] - (De[4] / De[1]) * dV + (De[0] / De[1]) * dC
        pH2 = pb[2] + (De[4] / De[2]) * dV - 3 * (De[0] / De[2]) * dC
        pCO = pb[3] - (De[4] / De[3]) * dV - (De[0] / De[3]) * dC
        return pc, pH2O, pH2, pCO, pv

    # codehealth: allow params -- the bulk state the pellet sits in
    def cell_rates(self, U, pb, T, De, act=None):
        return self.case.rates(self.dependent(U[:, 0], U[:, 1], pb, De), T, act)

    # codehealth: allow params -- the bulk state, the diffusivities and the warm start
    def solve(self, pb, T, De, act=None, guess=None, tol=1e-12):
        """Converged cell values ``(n, 2)`` of p_CH4, p_CO2 (bar)."""
        case, n = self.case, self.n
        c = R * T * case.rho_s * KMOLH / BAR             # bar/m2 per kmol/(kg h), times 1/D
        dscale = np.array([De[0], De[4]]) / c
        Ad = sparse.diags(np.tile(dscale, n)) @ self.lap
        bvec = (self.lap_b * np.array([pb[0], pb[4]]) * dscale).ravel()
        a = self.active[:, None]

        def src(U):
            r1, r2, r3 = self.cell_rates(U, pb, T, De, act)
            return a * np.stack([-(r1 + r3), r2 + r3], axis=1)

        def fun(u):
            U = u.reshape(n, 2)
            s, J = self.numjac(src, U)
            return Ad @ u + bvec + s.ravel(), Ad + J

        U0 = np.tile([pb[0], pb[4]], (n, 1)) if guess is None else guess
        res = newton(fun, U0.ravel(), tol=tol, maxfev=60)
        if not res.success:
            raise RuntimeError("pellet Newton solve did not converge")
        return np.asarray(res.x).reshape(n, 2)

    # codehealth: allow params -- the bulk state and the continuation length
    def solve_from_scratch(self, pb, T, De, act=None, steps=25):
        """Activity continuation from 1e-6 to the target, for a first solve."""
        a_t = self.case.act if act is None else act
        U = None
        for a in np.geomspace(1e-6, a_t, steps):
            U = self.solve(pb, T, De, act=a, guess=U)
        return U

    # codehealth: allow params -- the bulk state the pellet sits in
    def observed(self, U, pb, T, De, act=None):
        """Volume-averaged r1, r2, r3 over the whole layer (kmol/(kg h), layer basis)."""
        r = np.array(self.cell_rates(U, pb, T, De, act))
        return (r * self.active) @ self.vol_c / self.vol_c.sum()

    def surface_flux_rate(self, U, pb, T, De):
        """r_CH4,obs from the surface gradient of p_CH4 (slab only), eq. 9's flux identity."""
        bc = ({"a": 0.0, "b": 1.0, "d": pb[0]}, {"a": 1.0, "b": 0.0, "d": 0.0})
        _, g_left, _, _ = compute_boundary_values(U[:, 0], self.x_f, self.x_c, bc)
        c = R * T * self.case.rho_s * KMOLH / BAR
        return float(-De[0] * np.ravel(g_left)[0] / (c * self.case.l_act))

    def prolong(self, U, other):
        """Interpolate a solution of ``other`` (a coarser pellet) onto this grid."""
        return np.stack([np.interp(self.x_c, other.x_c, U[:, k]) for k in range(2)], axis=1)


# codehealth: allow params -- the bulk state, the pellet grid and the Richardson ladder
def pellet_rates(case, pb, T, De, n=80, levels=2, geom="slab", cut=None, guesses=None,
                 act=None, return_all=False):
    """Richardson-extrapolated observed rates (layer basis) from n, 2n (, 4n) cells.

    Returns ``(r_ext, values)`` with ``values`` per level when ``return_all``; the
    extrapolation assumes the scheme's second order, checked by ``observed_orders``.
    """
    vals, sols, pels = [], [], []
    for k in range(levels):
        pel = Pellet(case, n * 2 ** k, geom=geom, cut=cut)
        g = None if guesses is None else guesses[k]
        if g is None and k > 0:
            g = pel.prolong(sols[-1], pels[-1])
        U = pel.solve_from_scratch(pb, T, De, act) if g is None else pel.solve(pb, T, De, act, g)
        sols.append(U)
        pels.append(pel)
        vals.append(pel.observed(U, pb, T, De, act))
    v = np.array(vals)
    ext = v[-1] + (v[-1] - v[-2]) / 3.0
    if return_all:
        return ext, v, sols, pels
    return ext


class Tube:
    """Plug-flow tube, eqs. 10a-c, marched in z by implicit Euler with a nested pellet solve.

    Parameters
    ----------
    case : Case
    pt : callable
        Total pressure (bar) against z (m).
    q : callable, optional
        Imposed inner-wall heat flux (W/m2) against z (zero beyond the heated length).
    Tw, U : callable, float, optional
        Imposed external wall temperature (K) against z and the overall coefficient
        (W/(m2 K), inner area): q = U (T_w - T) on the heated length.
    n_r : int
        Pellet cells on the coarse pellet level (two levels, Richardson).
    """

    # codehealth: allow params -- the two heat forms, the pressure profile and the pellet grid
    def __init__(self, case, pt, q=None, Tw=None, U=None, n_r=80, pellet_levels=2, cut=None,
                 intrinsic=False):
        self.case, self.pt, self.q_fun, self.Tw, self.U = case, pt, q, Tw, U
        self.intrinsic = intrinsic                    # True: eta = 1 (bulk rates, no pellet)
        self.n_r, self.levels = n_r, pellet_levels
        self.pels = [Pellet(case, n_r * 2 ** k, cut=cut) for k in range(pellet_levels)]
        self.guess = None

    def q(self, z, T):
        if z > self.case.L_heat + 1e-12:
            return 0.0
        if self.q_fun is not None:
            return self.q_fun(z)
        return self.U * (self.Tw(z) - T)

    def robs(self, x, y, T, z):
        """Observed rates on the ring basis (kmol/(kg h)) at the bulk state."""
        case = self.case
        pt = self.pt(z)
        pb, yf = case.state(x, y, pt)
        if self.intrinsic:
            return np.array(case.rates(pb, T))
        De = case.d_eff(yf, T, pt)
        if self.guess is None:
            g0 = self.pels[0].solve_from_scratch(pb, T, De)
            self.guess = [g0] + [self.pels[k].prolong(g0, self.pels[0])
                                 for k in range(1, self.levels)]
        vals = []
        for k, pel in enumerate(self.pels):
            self.guess[k] = pel.solve(pb, T, De, guess=self.guess[k])
            vals.append(pel.observed(self.guess[k], pb, T, De))
        v = np.array(vals)
        ext = v[-1] + (v[-1] - v[-2]) / 3.0 if self.levels > 1 else v[-1]
        return ext / case.ring_factor

    def march(self, z, tol=1e-11):
        """States (x, y, T) at the nodes ``z`` (z[0] = 0, the inlet) and the heat sum (W).

        Returns
        -------
        Y : ndarray ``(len(z), 3)``; Qsum : ndarray, the heat added up to each node by the
        same right-point rule the energy residual uses (for the invariant).
        """
        case = self.case
        k_ = case.omega * case.rho_B * KMOLH / case.F0
        Y = np.empty((len(z), 3))
        Y[0] = (0.0, 0.0, case.T0)
        Qs = np.zeros(len(z))
        Hs = 1e5                                      # W, scale of the energy residual
        perim = math.pi * case.d_ti
        J = None
        for i in range(1, len(z)):
            dz, zi = z[i] - z[i - 1], z[i]
            Hk = case.h_flow(*Y[i - 1])
            xk, yk = Y[i - 1, 0], Y[i - 1, 1]

            def res(u):
                x, y, T = u
                r = self.robs(x, y, T, zi)
                return np.array([x - xk - dz * k_ * (r[0] + r[2]),
                                 y - yk - dz * k_ * (r[1] + r[2]),
                                 (case.h_flow(x, y, T) - Hk - dz * perim * self.q(zi, T)) / Hs])

            def fun(u):
                nonlocal J
                f0 = res(u)
                if J is None:
                    J = np.empty((3, 3))
                    h = np.array([1e-7, 1e-7, 1e-4])
                    for j in range(3):
                        e = np.zeros(3)
                        e[j] = h[j]
                        J[:, j] = (res(u + e) - f0) / h[j]
                return f0, J
            u0 = Y[i - 1].copy()
            try:
                out = newton(fun, u0, tol=tol, maxfev=40)
                ok = bool(getattr(out, "success", True))
            except Exception:
                ok = False
            if not ok:                                # refresh the chord Jacobian once
                J = None
                out = newton(fun, u0, tol=tol, maxfev=60)
                if not out.success and not self.intrinsic:     # eta = 1 is the stiff break-row mode, kept as is
                    raise RuntimeError(f"tube Newton solve did not converge at z = {zi:.4g} m")
            Y[i] = np.asarray(out.x)
            Qs[i] = Qs[i - 1] + dz * perim * self.q(zi, Y[i, 2])
            if i % 10 == 0:
                J = None                              # chord Jacobian refreshed every 10 steps
        return Y, Qs


# codehealth: allow params -- picklable description of one tube for a process pool
def run_tube(case, heat, pt_z, pt_v, n_z, n_r=80, cut=None, intrinsic=False):
    """One tube march; picklable arguments only (for a process pool).

    Parameters
    ----------
    heat : tuple
        ``("q", coef)``: imposed inner-wall flux ``polyval(coef, z)`` (W/m2) on the heated
        length; ``("Tw", z_nodes, Tw_nodes, U)``: imposed external wall temperature
        (linear interpolation) with overall coefficient U (W/(m2 K)).
    pt_z, pt_v : array_like
        Total pressure (bar) at nodes, interpolated linearly.

    Returns
    -------
    z, Y, Qs : ndarray
        Nodes, states (x, y, T) and the heat added up to each node (W).
    """
    def pt(z):
        return float(np.interp(z, pt_z, pt_v))
    if heat[0] == "q":
        coef = np.asarray(heat[1], dtype=float)
        tb = Tube(case, pt, q=lambda z: float(np.polyval(coef, z)), n_r=n_r, cut=cut,
                  intrinsic=intrinsic)
    else:
        _, tz, tv, U = heat
        tb = Tube(case, pt, Tw=lambda z: float(np.interp(z, tz, tv)), U=U, n_r=n_r, cut=cut,
                  intrinsic=intrinsic)
    z = z_grid(n_z)
    Y, Qs = tb.march(z)
    return z, Y, Qs


def z_grid(n, nodes=(0.25, 4.0, 11.12, 12.0), weights=(0.25, 0.3, 0.3, 0.15), betas=(4.0, 0.0, 0.0, 4.0)):
    """Graded axial nodes with nodes at the given z for every ``n`` (a multiple of 20).

    Each segment maps a uniform ``s`` by ``expm1(beta s) / expm1(beta)`` (uniform for
    beta 0), graded to the inlet layer and to the end of heating; doubling ``n`` nests
    the grids, so the first-order march error extrapolates (Richardson).
    """
    parts, lo = [np.zeros(1)], 0.0
    for hi, w, b in zip(nodes, weights, betas):
        m = int(round(n * w))
        s = np.linspace(0.0, 1.0, m + 1)[1:]
        seg = np.expm1(b * s) / np.expm1(b) if b > 0 else s
        parts.append(lo + (hi - lo) * seg)
        lo = hi
    return np.concatenate(parts)


def thiele_slab_eta(phi, n, beta=6.0):
    """Isothermal first-order slab on the graded grid: eta from the volume average and from
    the surface flux (for the V2 check against tanh(phi)/phi)."""
    s = np.linspace(0.0, 1.0, n + 1)
    x_f = np.expm1(beta * s) / np.expm1(beta)
    x_c = 0.5 * (x_f[1:] + x_f[:-1])
    bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
    G, Gb = construct_grad((n, 1), x_f, x_c, bc, axis=0)
    Dv = construct_div((n, 1), x_f, nu=0, axis=0)    # planar
    A = (Dv @ G - phi ** 2 * sparse.eye(n)).tocsc()
    c = sparse.linalg.spsolve(A, -_dense(Dv @ Gb))
    eta_vol = c @ np.diff(x_f)
    _, g_left, _, _ = compute_boundary_values(c.reshape(n, 1), x_f, x_c, bc)
    eta_flux = -float(np.ravel(g_left)[0]) / phi ** 2
    return eta_vol, eta_flux
