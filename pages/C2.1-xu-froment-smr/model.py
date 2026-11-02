"""Xu and Froment (1989) steam-reforming kinetics and an integral plug-flow bed.

Model code of gallery page C2.1. The kinetics know nothing about the reactor
and the reactor knows nothing about the chemistry beyond a pointwise source of
its two conversions, so either piece is reusable on its own.

Contents
--------
STOICH, ATOMS
    Moles of CH4, H2O, H2, CO, CO2 per mol CH4 fed per unit of the two
    conversions, and the C, H, O content of each species.
XuFroment
    Table 6 preexponentials, Table 5 energies, the equilibrium correlation and
    the three rate equations (the paper's Eqs. 3).
partial_pressures
    Partial pressures from the two conversions for a CH4/H2O/H2 feed.
SteamReformer
    Isothermal, isobaric integral plug flow in the space time ``tau =
    W/F0_CH4``: first-order upwind convection and ``pymrm.newton``.
reference_profile, reference_solution
    The same two ODEs integrated by ``scipy.integrate.solve_ivp`` (Radau), the
    independent route for the discretisation error; ``reference_solution``
    keeps the dense output, for root-finding on the profile.
half_ulp, round_trip
    The Table 5 to Table 6 round trip, each row scored against its own printed
    precision.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)
from pymrm.checks import residual_check

__all__ = ["R_KJ", "SPECIES", "STOICH", "ATOMS", "XuFroment", "partial_pressures",
           "SteamReformer", "reference_profile", "reference_solution", "half_ulp", "round_trip"]

R_KJ = 8.314e-3  # gas constant, kJ/(mol K)
SPECIES = ("CH4", "H2O", "H2", "CO", "CO2")

# Moles per mol CH4 fed: n = n_feed + STOICH @ (x_CH4, x_CO2). x_CH4 is all
# methane converted, x_CO2 the part of it that ended up as CO2.
STOICH = np.array([[-1.0, 0.0],    # CH4  1 - x
                   [-1.0, -1.0],   # H2O  sr - x - y
                   [3.0, 1.0],     # H2   hr + 3x + y
                   [1.0, -1.0],    # CO   x - y
                   [0.0, 1.0]])    # CO2  y
ATOMS = np.array([[1, 0, 0, 1, 1],    # C
                  [4, 2, 2, 0, 0],    # H
                  [0, 1, 0, 1, 2]],   # O
                 dtype=float)


def partial_pressures(x_ch4, x_co2, p_t, sr, hr, stoich=STOICH):
    """Partial pressures (bar) of CH4, H2O, H2, CO, CO2 from the two conversions.

    Parameters
    ----------
    x_ch4, x_co2 : array_like
        Total methane conversion and conversion into CO2.
    p_t : float
        Total pressure (bar).
    sr, hr : float
        Feed ratios H2O/CH4 and H2/CH4.
    stoich : ndarray, optional
        ``(5, 2)`` change of each species per unit of the two conversions.

    Returns
    -------
    tuple of ndarray
        ``(p_CH4, p_H2O, p_H2, p_CO, p_CO2)``.
    """
    x_ch4, x_co2 = np.asarray(x_ch4, dtype=float), np.asarray(x_co2, dtype=float)
    feed = np.array([1.0, sr, hr, 0.0, 0.0])
    n = [feed[i] + stoich[i, 0] * x_ch4 + stoich[i, 1] * x_co2 for i in range(5)]
    n_t = sum(n)
    return tuple(p_t * n_i / n_t for n_i in n)


@dataclass(frozen=True)
class XuFroment:
    """Xu and Froment (1989) intrinsic kinetics, rates in kmol/(kg_cat h).

    Parameters
    ----------
    A : dict
        Preexponential factor per symbol (``k1, k2, k3, K_CO, K_H2, K_CH4,
        K_H2O``), Table 6.
    E : dict
        Activation energy or adsorption enthalpy per symbol (kJ/mol), Table 5.
    eq : tuple of float, optional
        ``(a1, b1, a2, b2)`` of the equilibrium correlation ``K1 = exp(-a1/T +
        b1)`` (bar^2) and ``K2 = exp(a2/T + b2)``, with ``K3 = K1 K2``. Not in
        the paper; the correlation usually quoted with these kinetics.
    eq_factor : float, optional
        Multiplies each of K1, K2, K3 (a sensitivity knob; 1 as published).
    activity : float, optional
        Multiplies k1, k2, k3 (1 for the steam-reforming reference level of
        Tables 5 and 6).
    """

    A: dict
    E: dict
    eq: tuple = (26830.0, 30.114, 4400.0, -4.036)
    eq_factor: float = 1.0
    activity: float = 1.0
    stoich: np.ndarray = field(default_factory=STOICH.copy, repr=False)

    @classmethod
    def from_table(cls, par, **changes):
        """Build from the parameter table (columns ``symbol, A, E_or_dH_kJ_per_mol``)."""
        A = {s: float(a) for s, a in zip(par["symbol"], par["A"])}
        E = {s: float(e) for s, e in zip(par["symbol"], par["E_or_dH_kJ_per_mol"])}
        return cls(A=A, E=E, **changes)

    def edited(self, symbol, A=None, E=None):
        """A copy with one parameter's preexponential or energy replaced."""
        new_A, new_E = dict(self.A), dict(self.E)
        if A is not None:
            new_A[symbol] = A
        if E is not None:
            new_E[symbol] = E
        return replace(self, A=new_A, E=new_E)

    def constants(self, T):
        """Rate and adsorption constants at ``T`` (K): ``A exp(-E / (R T))``."""
        kc = {s: self.A[s] * np.exp(-self.E[s] / (R_KJ * T)) for s in self.A}
        for s in ("k1", "k2", "k3"):
            kc[s] *= self.activity
        return kc

    def equilibrium(self, T):
        """``(K1, K2, K3)`` at ``T`` (K); K1 and K3 in bar^2, K2 dimensionless."""
        a1, b1, a2, b2 = self.eq
        K1, K2 = np.exp(-a1 / T + b1), np.exp(a2 / T + b2)
        f = self.eq_factor
        return f * K1, f * K2, f * K1 * K2

    def vant_hoff_enthalpies(self):
        """Reaction enthalpies (kJ/mol) implied by the correlation's slopes, I, II, III."""
        a1, _, a2, _ = self.eq
        return {"I": a1 * R_KJ, "II": -a2 * R_KJ, "III": (a1 - a2) * R_KJ}

    def rates(self, x_ch4, x_co2, T, p_t, sr, hr, kc=None):
        """``(r1, r2, r3)`` of the paper's Eqs. 3 at the given conversions."""
        kc = self.constants(T) if kc is None else kc
        K1, K2, K3 = self.equilibrium(T)
        p_ch4, p_h2o, p_h2, p_co, p_co2 = partial_pressures(x_ch4, x_co2, p_t, sr, hr,
                                                            self.stoich)
        den = (1.0 + kc["K_CO"] * p_co + kc["K_H2"] * p_h2
               + kc["K_CH4"] * p_ch4 + kc["K_H2O"] * p_h2o / p_h2) ** 2
        r1 = kc["k1"] / p_h2 ** 2.5 * (p_ch4 * p_h2o - p_h2 ** 3 * p_co / K1) / den
        r2 = kc["k2"] / p_h2 * (p_co * p_h2o - p_h2 * p_co2 / K2) / den
        r3 = kc["k3"] / p_h2 ** 3.5 * (p_ch4 * p_h2o ** 2 - p_h2 ** 4 * p_co2 / K3) / den
        return r1, r2, r3

    def beta(self, x_ch4, x_co2, T, p_t, sr, hr):
        """Approach to equilibrium ``Q_i / K_i`` of reactions I, II, III (1 = equilibrated)."""
        K1, K2, K3 = self.equilibrium(T)
        p_ch4, p_h2o, p_h2, p_co, p_co2 = partial_pressures(x_ch4, x_co2, p_t, sr, hr,
                                                            self.stoich)
        return (p_h2 ** 3 * p_co / (p_ch4 * p_h2o * K1),
                p_h2 * p_co2 / (p_co * p_h2o * K2),
                p_h2 ** 4 * p_co2 / (p_ch4 * p_h2o ** 2 * K3))


class SteamReformer:
    """Isothermal, isobaric integral plug flow in the space time ``tau = W/F0_CH4``.

    ``d(v u)/dtau = (r_CH4, r_CO2)`` with ``v = 1`` and ``u = (x_CH4, x_CO2)``,
    the paper's Eq. 1. State layout ``(n_tau, 2)``: space time first, the two
    conversions last.

    Parameters
    ----------
    kin : XuFroment
        The kinetics.
    T : float
        Temperature (K).
    p_t, sr, hr : float, optional
        Total pressure (bar) and feed ratios H2O/CH4 and H2/CH4.
    tau_max : float, optional
        Bed length in space time (g_cat h / mol_CH4).
    n_tau : int, optional
        Number of cells.
    """

    def __init__(self, kin, T, p_t=10.0, sr=3.0, hr=1.25, tau_max=0.45, n_tau=1000):
        self.kin, self.T, self.p_t, self.sr, self.hr = kin, float(T), p_t, sr, hr
        self.shape = (int(n_tau), 2)
        self.tau_f = np.linspace(0.0, tau_max, self.shape[0] + 1)
        self.tau_c = 0.5 * (self.tau_f[:-1] + self.tau_f[1:])
        # tau = 0, outward normal -tau: no methane converted yet, x = 0  -> a = 0, b = 1, d = 0
        # tau = tau_max, outward normal +tau: what leaves is what the last
        #     cell holds; no condition is imposed on a pure convective outlet
        self.bc = ({"a": 0.0, "b": 1.0, "d": 0.0}, {"outflow": True})
        # g/mol and kg/kmol are the same scale, so tau in g h/mol takes the
        # rates in kmol/(kg h) without a unit factor
        self.kc = kin.constants(self.T)
        self.u = np.zeros(self.shape)
        self.result = None
        self._build_operators()

    def _build_operators(self):
        # The flow coordinate is tau and the velocity is 1 by construction.
        # nu = 0: the space time is a Cartesian coordinate.
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.tau_f, self.tau_c,
                                                      self.bc, v=1.0, axis=0)
        div_mat = construct_div(self.shape, self.tau_f, nu=0, axis=0)
        self.jac_const = (div_mat @ conv_mat).tocsc()
        self.g_const = (div_mat @ conv_bc).toarray().reshape(-1, 1)
        # The source at a cell depends on that cell's two conversions only, so
        # the default stencil (last axis, the two fields, coupled in full) is exact
        # in structure: a 2 x 2 block per cell.
        self.numjac = NumJac(self.shape)

    def reaction(self, u):
        """Pointwise source ``(r_CH4, r_CO2) = (r1 + r3, r2 + r3)``, shape ``(n_tau, 2)``."""
        r1, r2, r3 = self.kin.rates(u[..., 0], u[..., 1], self.T, self.p_t, self.sr,
                                    self.hr, self.kc)
        return np.stack([r1 + r3, r2 + r3], axis=-1)

    def residual(self, u):
        """Steady residual ``div(v u) - r(u)`` and its Jacobian."""
        g_rxn, jac_rxn = self.numjac(self.reaction, u.reshape(self.shape))
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self, maxfev=60, strict=True):
        """Newton from ``x = 0``; raises ``RuntimeError`` unless it converges.

        ``strict=False`` keeps whatever state the iteration stopped at.
        """
        result = newton(self.residual, np.zeros(self.shape), maxfev=maxfev)
        if strict and not result.success:
            raise RuntimeError(f"SteamReformer.solve at {self.T} K: {result.message}")
        self.u, self.result = result.x.reshape(self.shape), result
        return self

    def backward_error(self):
        """Componentwise backward error of the current state (``pymrm.checks``)."""
        return residual_check(self.residual, self.u.reshape(-1))["backward_error"]

    def at(self, tau, field=0, read="centre"):
        """Conversion ``field`` (0: x_CH4, 1: x_CO2) at space time ``tau``.

        ``read="centre"`` interpolates linearly between cell centres;
        ``"face"`` between the upwind face values (the inlet value and each
        cell's value at its downstream face); ``"nearest"`` takes the nearest
        cell value. Past the last centre the outlet face value is used.
        """
        tau = np.asarray(tau, dtype=float)
        if np.any(tau < self.tau_c[0]) or np.any(tau > self.tau_f[-1]):
            raise ValueError("SteamReformer.at: tau outside the computed range")
        col = self.u[:, field]
        if read == "centre":
            return np.interp(tau, np.append(self.tau_c, self.tau_f[-1]),
                             np.append(col, self.outlet()[field]))
        if read == "face":
            return np.interp(tau, self.tau_f, np.concatenate([[0.0], col]))
        if read == "nearest":
            return col[np.abs(self.tau_c[None, :] - np.atleast_1d(tau)[:, None]).argmin(axis=1)]
        raise ValueError(f"unknown read {read!r}")

    def outlet(self):
        """``(x_CH4, x_CO2)`` on the outlet face, read with ``compute_boundary_values``."""
        _, _, u_out, _ = compute_boundary_values(self.u, self.tau_f, self.tau_c, self.bc,
                                                 axis=0)
        return np.asarray(u_out, dtype=float).ravel()

    def beta(self, tau):
        """``Q_i / K_i`` of reactions I, II, III at space time ``tau``."""
        return self.kin.beta(self.at(tau, 0), self.at(tau, 1), self.T, self.p_t,
                             self.sr, self.hr)


def reference_solution(kin, T, tau_max, p_t=10.0, sr=3.0, hr=1.25, rtol=1e-12, atol=1e-14):
    """The plug-flow ODEs integrated by ``solve_ivp`` (Radau) with dense output.

    Returns
    -------
    callable
        ``f(tau)`` giving ``(2, len(tau))``: x_CH4 and x_CO2 at the space times
        ``tau`` in ``[0, tau_max]``.
    """
    kc = kin.constants(T)

    def rhs(_, u):
        r1, r2, r3 = kin.rates(u[0], u[1], T, p_t, sr, hr, kc)
        return [r1 + r3, r2 + r3]

    sol = solve_ivp(rhs, (0.0, float(tau_max)), [0.0, 0.0], method="Radau",
                    rtol=rtol, atol=atol, dense_output=True)
    if not sol.success:
        raise RuntimeError(f"reference_solution at {T} K: {sol.message}")
    return sol.sol


def reference_profile(kin, T, taus, p_t=10.0, sr=3.0, hr=1.25, rtol=1e-12, atol=1e-14):
    """The same plug-flow ODEs integrated by ``solve_ivp`` (Radau).

    Returns
    -------
    ndarray
        ``(len(taus), 2)``: x_CH4 and x_CO2 at each space time in ``taus``.
    """
    taus = np.asarray(taus, dtype=float)
    return reference_solution(kin, T, taus.max(), p_t, sr, hr, rtol, atol)(taus).T


def half_ulp(s):
    """Half a unit in the last printed place of a number written as a string."""
    m, _, e = str(s).strip().lower().partition("e")
    dec = len(m.split(".")[1]) if "." in m else 0
    return 0.5 * 10.0 ** (-dec) * 10.0 ** (int(e) if e else 0)


def round_trip(rows, T_ref=None):
    """Table 5 to Table 6 round trip, each row against its own printed precision.

    ``A = value(T_ref) exp(E / (R T_ref))``. Half a unit in the last printed
    place of ``A``, of ``value(T_ref)`` and of ``E`` propagates to the bound
    ``u(A)/A + u(value)/value + u(E)/(R T_ref)``.

    Parameters
    ----------
    rows : iterable of dict
        Each with the printed strings ``A``, ``E``, ``val`` and the number
        ``T_ref`` (K), plus ``symbol``.
    T_ref : float, optional
        Force one reference temperature on every row.

    Returns
    -------
    list of dict
        Per row: ``symbol, T_ref, A, A_computed, rel_pct, bound_pct, ratio``.
    """
    out = []
    for r in rows:
        A, E, val = float(r["A"]), float(r["E"]), float(r["val"])
        Tr = float(r["T_ref"] if T_ref is None else T_ref)
        A_c = val * np.exp(E / (R_KJ * Tr))
        rel = abs(A_c - A) / A * 100
        bound = (half_ulp(r["A"]) / abs(A) + half_ulp(r["val"]) / abs(val)
                 + half_ulp(r["E"]) / (R_KJ * Tr)) * 100
        out.append(dict(symbol=r["symbol"], T_ref=Tr, A=A, A_computed=A_c, rel_pct=rel,
                        bound_pct=bound, ratio=rel / bound))
    return out
