"""Borman and Westerterp (1995): ethene epoxidation and combustion over silver, model III.

Rate law (Table 5, model III), selectivity, analytic derivatives, the gas composition of
the gradientless (Berty) reactor (eqs. 1, 2; Table 8; A-4, A-5), eq. 8 and a forward
solve of the Berty reactor with ``pymrm.newton``. Partial pressures in bar, rates in
mol/(kg s), temperature in K, catalyst mass in kg. No plotting, no data loading.

Species order of a pressure array ``P[..., j]``: ethene E, oxygen O, carbon dioxide C,
ethene oxide EO, water W (``SPECIES``).
"""
from __future__ import annotations

import numpy as np
from pymrm import NumJac, newton
from pymrm.checks import residual_check

SPECIES = ("E", "O", "C", "EO", "W")
R_GAS = 8.314462618              # J/(mol K), CODATA 2018 to ten significant digits
T_REF = 476.0                    # K, temperature of the printed K_j (p. 54)
PRODUCTS = ("C", "EO", "W")
M_ETHENE = 28.0532               # g/mol
M_AIR = 28.9647                  # g/mol, dry air
Y_O2_AIR = 0.20946               # O2 mole fraction of dry air


def table6_params(table6, model="III"):
    """The constants of one model of Table 6 as floats.

    Parameters
    ----------
    table6 : pandas.DataFrame
        Columns ``model``, ``reaction``, ``parameter``, ``value`` (strings, as printed).
    model : str
        ``"I"``, ``"II"`` or ``"III"``.

    Returns
    -------
    dict
        ``{"r1": p, "r2": p}``; each ``p`` has ``k_inf``, ``T_act``, ``k_lump``,
        ``T_lump``, ``K`` (``{j: K_j at 476 K}``) and ``pair`` (``{j: (K_inf, T_ads)}``
        for the constants Table 6 gives a temperature dependence), plus ``n``, ``f2``,
        ``f3`` where printed.
    """
    return {rx: _column(table6[(table6.model == model) & (table6.reaction == rx)])
            for rx in ("r1", "r2")}


def _column(rows):
    """One column (reaction) of Table 6 as floats."""
    v = dict(zip(rows.parameter, rows.value.astype(float)))
    p = {k: v[k] for k in ("k_inf", "T_act", "k_lump", "T_lump")}
    p["K"] = {j: v[f"K_{j}"] for j in SPECIES if f"K_{j}" in v}
    p["pair"] = {j: (v[f"K_{j}_inf"], v[f"T_ads_{j}"]) for j in SPECIES if f"K_{j}_inf" in v}
    p.update({e: v[e] for e in ("n", "f2", "f3") if e in v})
    return p


def adsorption(p, T, frozen=False):
    """``{j: K_j(T)}`` in 1/bar (caption of Table 5).

    A constant with a printed ``(K_inf, T_ads)`` pair is ``K_inf exp(T_ads/T)``; the
    others keep their printed 476 K value at every temperature (p. 54). ``frozen``
    evaluates every pair at 476 K instead (a break-row option).
    """
    K = dict(p["K"])
    for j, (k_inf, t_ads) in p["pair"].items():
        K[j] = k_inf * np.exp(t_ads / (T_REF if frozen else T))
    return K


def denominator(p, T, P, frozen=False):
    """``D = 1 + K_E P_E + (K_O P_O)^(1/2) + K_C P_C + K_W P_W + K_EO P_EO`` (model III)."""
    K = adsorption(p, T, frozen)
    P = np.asarray(P, dtype=float)
    return (1.0 + K["E"] * P[..., 0] + np.sqrt(K["O"] * P[..., 1]) + K["C"] * P[..., 2]
            + K["EO"] * P[..., 3] + K["W"] * P[..., 4])


# codehealth: allow params -- one rate law; form and frozen are the two readings the page tests
def rate(p, T, P, form="caption", frozen=False):
    """Model III rate ``r = k_r K_E K_O^(1/2) P_E P_O^(1/2) / D^2`` in mol/(kg s).

    ``form="caption"`` builds the numerator from ``k_inf``, ``T_act`` and ``K_E(T)``
    (Table 5 caption); ``form="lump"`` from the printed ``(k_r)_lump`` and
    ``(T_r)_lump``. The denominator is the same in both. ``frozen`` freezes the
    temperature dependence of the adsorption constants at 476 K in ``D`` only.
    """
    P = np.asarray(P, dtype=float)
    if form == "caption":
        K = adsorption(p, T)
        num = p["k_inf"] * np.exp(-p["T_act"] / T) * K["E"] * np.sqrt(K["O"])
    elif form == "lump":
        num = p["k_lump"] * np.exp(-p["T_lump"] / T)
    else:
        raise ValueError(f"unknown form {form!r}")
    return num * P[..., 0] * np.sqrt(P[..., 1]) / denominator(p, T, P, frozen) ** 2


# codehealth: allow params -- the rate law plus the r2 factor F the page is about
def rates(par, T, P, F=1.0, form="caption", frozen=False):
    """``(r1, r2)``: epoxidation (EO production) and combustion (CO2 production),
    the latter multiplied by the factor ``F`` on its rate constant."""
    return (rate(par["r1"], T, P, form, frozen), F * rate(par["r2"], T, P, form, frozen))


def selectivity(r1, r2, half=True):
    """Selectivity to EO, ``r1/(r1 + r2/2)`` (A-3 with eqs. 1, 2: two CO2 per ethene).

    ``half=False`` drops the factor 1/2 (a break-row option)."""
    return r1 / (r1 + (0.5 if half else 1.0) * r2)


def factor_for_selectivity(r1, r2, sigma, half=True):
    """The factor ``F`` on ``r2`` that gives selectivity ``sigma`` (closed form)."""
    return r1 * (1.0 - sigma) / (sigma * r2 * (0.5 if half else 1.0))


def apparent_activation(p, T, P, frozen=False):
    """Apparent activation energy ``-R dln r/d(1/T)`` in J/mol (caption form).

    ``R [T_act - T_ads,E + 2 sum_j T_ads,j K_j P_j / D]`` over the constants with a
    printed pair (the oxygen term would carry a factor 1/2; model III prints none).
    """
    P = np.asarray(P, dtype=float)
    K = adsorption(p, T, frozen)
    D = denominator(p, T, P, frozen)
    t_num = p["T_act"] - (p["pair"]["E"][1] if "E" in p["pair"] else 0.0)
    s = 0.0
    if not frozen:
        for j, (_, t_ads) in p["pair"].items():
            i = SPECIES.index(j)
            term = np.sqrt(K[j] * P[..., i]) / 2.0 if j == "O" else K[j] * P[..., i]
            s = s + t_ads * term / D
    return R_GAS * (t_num + 2.0 * s)


def orders(p, T, P):
    """Reaction orders ``(n_O, n_E)``: ``1/2 - (K_O P_O)^(1/2)/D`` and ``1 - 2 K_E P_E/D``."""
    P = np.asarray(P, dtype=float)
    K = adsorption(p, T)
    D = denominator(p, T, P)
    return 0.5 - np.sqrt(K["O"] * P[..., 1]) / D, 1.0 - 2.0 * K["E"] * P[..., 0] / D


def dln_ratio_dP(par, T, P, j):
    """``d ln(r1/r2)/dP_j`` (1/bar) for a product ``j``: ``-2 K_j^1/D_1 + 2 K_j^2/D_2``.

    Its sign is the sign of ``dsigma/dP_j`` whatever the factor on ``r2``."""
    out = 0.0
    for rx, sign in (("r1", -1.0), ("r2", 1.0)):
        out = out + sign * 2.0 * adsorption(par[rx], T)[j] / denominator(par[rx], T, P)
    return out


def dsigma_dT(par, T, P, F=1.0):
    """``dsigma/dT = -sigma (1 - sigma) (E_app,2 - E_app,1) / (R T^2)`` in 1/K."""
    s = selectivity(*rates(par, T, P, F))
    de = apparent_activation(par["r2"], T, P) - apparent_activation(par["r1"], T, P)
    return -s * (1.0 - s) * de / (R_GAS * T ** 2)


# codehealth: allow params -- the three printed constants of eq. 8, so a misreading is one argument
def r3_eq8(T, P_EO, k=1.7e5, T3=10500.0, n=0.4):
    """Eq. 8, ethene oxide oxidation measured without ethene, mol/(kg s)."""
    return k * np.exp(-T3 / T) * np.asarray(P_EO, dtype=float) ** n


# ---- gas composition of the gradientless reactor (Table 8, A-4, A-5) ----------------

def air_feed(x_E, y_O2=Y_O2_AIR):
    """Feed mole fractions for ethene in dry, CO2-free air; returns ``(x, M)``."""
    x = {"E": x_E, "O": (1.0 - x_E) * y_O2, "C": 0.0, "EO": 0.0, "W": 0.0}
    x["N"] = 1.0 - sum(x.values())
    return x, x_E * M_ETHENE + (1.0 - x_E) * M_AIR


def outlet_fractions(xI, zE, zC, o2_per_zC=2.5):
    """Table 8: outlet mole fractions after conversions ``zE`` (ethene) and ``zC``
    (ethene to CO2), divided by the total ``S = 1 - x_E^I (zE - zC)/2``.

    ``o2_per_zC`` is the O2 consumed per unit of ``zC`` beyond the 1/2 per unit of
    ``zE`` (3 - 1/2 = 5/2 by the stoichiometry of eq. E1; other values are break-row
    options)."""
    a = xI["E"]
    S = 1.0 - 0.5 * a * (zE - zC)
    x = {"E": a * (1.0 - zE),
         "O": xI["O"] - a * (0.5 * zE + o2_per_zC * zC),
         "C": xI["C"] + 2.0 * a * zC,
         "EO": xI["EO"] + a * (zE - zC),
         "W": xI["W"] + 2.0 * a * zC,
         "N": xI["N"]}
    return {k: v / S for k, v in x.items()}


def conversions_from_fractions(xI, x, contraction=True):
    """A-4 and A-5: ``(zE, zC)`` from the feed and outlet fractions of ethene and CO2.

    ``contraction=False`` ignores the change in total moles (``S = 1``), a break-row
    option."""
    a, c = xI["E"], xI["C"]
    if not contraction:
        return 1.0 - x["E"] / a, (x["C"] - c) / (2.0 * a)
    zE = (a - x["E"] - 0.25 * (x["C"] * a - x["E"] * c)) / (a * (1.0 - 0.25 * x["C"]
                                                                - 0.5 * x["E"]))
    zC = (x["C"] * (1.0 - 0.5 * a * zE) - c) / (2.0 * a * (1.0 - 0.25 * x["C"]))
    return zE, zC


def element_balances(xI, x):
    """Relative C, H and O balance gaps, outlet against feed, per mole of feed."""
    S = xI["N"] / x["N"]                  # moles out per mole of feed, from the inert
    atoms = {"C": {"E": 2, "C": 1, "EO": 2}, "H": {"E": 4, "EO": 4, "W": 2},
             "O": {"O": 2, "C": 2, "EO": 1, "W": 1}}
    gaps = {}
    for el, w in atoms.items():
        fin = sum(n * xI[k] for k, n in w.items())
        fout = sum(n * x[k] for k, n in w.items()) * S
        gaps[el] = fout / fin - 1.0
    return gaps


def pressures(x, P_total):
    """Partial pressures ``[E, O, C, EO, W]`` in bar."""
    return np.stack([np.asarray(x[j], dtype=float) * P_total for j in SPECIES], axis=-1)


# codehealth: allow params -- eqs. 1 and 2: two rates, catalyst mass, molar feed, feed fractions
def invert_berty(r1, r2, W, phi, xI):
    """Eqs. 1 and 2 solved for the conversions: ``zE x_E^I phi = (r1 + r2/2) W``,
    ``zC x_E^I phi = (r2/2) W`` (``W`` in kg, ``phi`` total molar feed in mol/s)."""
    n_E = xI["E"] * phi
    return (r1 + 0.5 * r2) * W / n_E, 0.5 * r2 * W / n_E


class BertyReactor:
    """The gradientless (Berty) reactor solved forward for its conversions (E12).

    Unknowns ``zE, zC`` scaled by ``scale`` (the expected conversions) so both are of
    order one (pitfall P7); the rates are model III at the outlet composition, with
    activity ``a1`` on ``r1``, ``a2`` on ``r2`` and the factor ``F`` on ``r2``.
    """

    # the operating condition of one run (T, P, feed, flow, mass), the two activities,
    # the r2 factor and the scale of the unknowns
    # codehealth: allow params -- one Berty run: operating condition, activities, factor
    def __init__(self, par, T, P_total, xI, phi, W, F=1.0, a1=1.0, a2=1.0, scale=(1e-2, 1e-2)):
        self.par, self.T, self.P_total, self.xI = par, T, P_total, xI
        self.phi, self.W, self.F, self.a1, self.a2 = phi, W, F, a1, a2
        self.scale = np.asarray(scale, dtype=float).reshape(1, 2)
        self.numjac = NumJac((1, 2))          # the two unknowns coupled in full

    def outlet_rates(self, zE, zC):
        """``(r1, r2)`` at the outlet composition of conversions ``zE, zC``."""
        P = pressures(outlet_fractions(self.xI, zE, zC), self.P_total)
        r1, r2 = rates(self.par, self.T, P, self.F)
        return self.a1 * r1, self.a2 * r2

    def _g(self, u):
        z = u * self.scale
        r1, r2 = self.outlet_rates(z[..., 0], z[..., 1])
        n_E = self.xI["E"] * self.phi
        g = np.stack([z[..., 0] * n_E - self.W * (r1 + 0.5 * r2),
                      z[..., 1] * n_E - self.W * 0.5 * r2], axis=-1)
        return g / (n_E * self.scale)

    def residual(self, u):
        g, jac = self.numjac(self._g, u.reshape(1, 2))
        return g.reshape(-1, 1), jac

    def solve(self, tol=1e-12, guess=(0.0, 0.0), backward_tol=1e-8):
        """Return ``(zE, zC, result, check)``; raises unless newton reports success.

        The default guess is the feed (no conversion); ``backward_tol`` goes to
        ``newton`` (``None`` stops on the step alone, a break-row option)."""
        u0 = np.asarray(guess, dtype=float).reshape(1, 2)
        result = newton(self.residual, u0, tol=tol, stagnation=True, maxfev=50,
                        backward_tol=backward_tol)
        if not result.success:
            raise RuntimeError(f"BertyReactor.solve: {result.message}")
        check = residual_check(self.residual, result.x, tol=1e-10)
        z = result.x.reshape(1, 2) * self.scale
        return float(z[0, 0]), float(z[0, 1]), result, check
