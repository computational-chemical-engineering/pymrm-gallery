#!/usr/bin/env python3
"""Generate index.ipynb for page C2.1.

Pages are generated from a builder rather than hand-edited JSON so that the
prose and the code stay in one reviewable file. Run from the page directory:

    python build_page.py && jupyter nbconvert --execute --inplace index.ipynb
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "Steam methane reforming: the Xu–Froment intrinsic kinetics"
description: "The most-used rate equations in reformer modelling, rebuilt from the 1989 tables and put back against the runs plotted in the paper's own Figures 2 and 3."
categories: [sec:C, struct:S1, struct:S2, tier:T1, data:tier2, data:tier4, phase:gas, phase:gas-solid]
date: 2026-07-26
---

# Steam methane reforming: the Xu–Froment intrinsic kinetics

**Catalog ID:** `C2.1` · **Structures:** `S1` (0D/1D reaction network), `S2` (plug flow with reaction) · **Tier:** T1

Almost every steam-reformer simulation published since 1989 uses the same three
rate equations. They come from one paper, and getting them wrong is easy: the
scan mangles the exponents, the reference temperature is not the same for all
seven parameters, and two different activity factors float around the text.
This page rebuilds them from the tables and puts them back against the
measurements they were fitted to."""))

cells.append(md(r"""## Background

Steam reforming is how most of the world's hydrogen is made. Methane and steam
go into rows of tubes packed with a nickel catalyst, sitting in a fired furnace
at 675–1000 K and 30 bar, and come out as synthesis gas. The tubes are the most
expensive part of the plant and they fail by creep, so how fast the reaction
actually runs — and therefore how much heat has to cross the tube wall, and
where — is a materials question as much as a chemical one.

Before 1989 most published kinetics were *effective* rate equations: fitted on
industrial-size pellets, with the (severe) internal diffusion limitation
absorbed into the constants. Those cannot be transferred to a different pellet,
and they cannot be combined with a diffusion model without double-counting.

Xu and Froment set out to measure the **intrinsic** kinetics instead — catalyst
crushed to 0.18–0.25 mm, where they show diffusion no longer limits — and to do
it mechanistically rather than by curve-fitting a power law. They enumerated 11
candidate reactions, cut them down thermodynamically, built two surface
mechanisms, generated 21 candidate sets of three rate equations from them, and
discriminated between the sets on 220 steam-reforming runs plus 60 runs on the
reverse water-gas shift and methanation.

What survived is the triangular scheme

$$
\mathrm{I}\;\;\mathrm{CH_4 + H_2O \rightleftharpoons CO + 3H_2}, \qquad
\mathrm{II}\;\;\mathrm{CO + H_2O \rightleftharpoons CO_2 + H_2}, \qquad
\mathrm{III}\;\;\mathrm{CH_4 + 2H_2O \rightleftharpoons CO_2 + 4H_2},
$$

with reaction III a *parallel* route to CO₂ rather than CO₂ being formed only
through CO. That was the paper's substantive claim, and it came from the data:
the ratio $V(\mathrm{II})$ does not extrapolate to the origin as space time
goes to zero, so some CO₂ must be produced directly from methane."""))

cells.append(md(r"""## The published model

**Rate equations** (the paper's Eqs. 3). All three share one denominator,
because all three are assumed to occur on the same active sites:

$$
r_1 = \frac{k_1}{p_{\mathrm{H_2}}^{2.5}}
      \left( p_{\mathrm{CH_4}} p_{\mathrm{H_2O}}
             - \frac{p_{\mathrm{H_2}}^{3} p_{\mathrm{CO}}}{K_1} \right) \Big/ \mathrm{DEN}^2
$$

$$
r_2 = \frac{k_2}{p_{\mathrm{H_2}}}
      \left( p_{\mathrm{CO}} p_{\mathrm{H_2O}}
             - \frac{p_{\mathrm{H_2}} p_{\mathrm{CO_2}}}{K_2} \right) \Big/ \mathrm{DEN}^2
$$

$$
r_3 = \frac{k_3}{p_{\mathrm{H_2}}^{3.5}}
      \left( p_{\mathrm{CH_4}} p_{\mathrm{H_2O}}^{2}
             - \frac{p_{\mathrm{H_2}}^{4} p_{\mathrm{CO_2}}}{K_3} \right) \Big/ \mathrm{DEN}^2
$$

$$
\mathrm{DEN} = 1 + K_{\mathrm{CO}} p_{\mathrm{CO}} + K_{\mathrm{H_2}} p_{\mathrm{H_2}}
             + K_{\mathrm{CH_4}} p_{\mathrm{CH_4}}
             + K_{\mathrm{H_2O}} \frac{p_{\mathrm{H_2O}}}{p_{\mathrm{H_2}}}
$$

There is **no $K_{\mathrm{CO_2}}$ term.** It was never statistically
significant, from either data set, and the paper says so explicitly. Adding one
is a common way to silently change the model.

**Species rates** (Eq. 4), for steam reforming:

$$
r_{\mathrm{CO}} = r_1 - r_2, \qquad
r_{\mathrm{CO_2}} = r_2 + r_3, \qquad
r_{\mathrm{CH_4}} = r_1 + r_3
$$

**Reactor** (Eq. 1). The kinetic runs were done in an integral fixed bed,
isothermal and isobaric, integrated in the space time $W/F^{0}_{\mathrm{CH_4}}$:

$$
\frac{\mathrm{d}x_{\mathrm{CH_4}}}{\mathrm{d}(W/F^{0}_{\mathrm{CH_4}})} = r_{\mathrm{CH_4}},
\qquad
\frac{\mathrm{d}x_{\mathrm{CO_2}}}{\mathrm{d}(W/F^{0}_{\mathrm{CH_4}})} = r_{\mathrm{CO_2}},
\qquad x = 0 \;\text{at}\; W/F^{0}_{\mathrm{CH_4}} = 0 .
$$

Two of the three species rates are independent, so two conversions describe the
whole composition: $x_{\mathrm{CH_4}}$, the total methane converted, and
$x_{\mathrm{CO_2}}$, the part of it that ended up as CO₂."""))

cells.append(md(r"""## Parameters and assumptions

**Assumptions:** isothermal and isobaric plug flow; intrinsic kinetics (the
paper demonstrates no internal limitation at 0.18–0.25 mm and calculates
external transport to be negligible); ideal gas; no carbon deposition — the
thermodynamic analysis rules it out at these conditions; and the catalyst held
at a fixed activity level, which is how the authors corrected their space times.

**Three traps in this paper.** Each one silently changes the answer:

1. **The reference temperature is not the same for all parameters.**
   $T_r = 648$ K for $k_1, k_2, k_3, K_{\mathrm{CO}}$ and $K_{\mathrm{H_2}}$,
   but $T_r = 823$ K for $K_{\mathrm{CH_4}}$ and $K_{\mathrm{H_2O}}$. Using one
   reference for all seven biases the temperature dependence.
2. **Three activity levels.** Table 6 is footnoted *reference activity* — the
   partially deactivated catalyst. Multiply the rate coefficients by **1.225**
   for the reverse water-gas shift and methanation reference level (the curves
   in the paper's Figures 4 and 5 were drawn that way), and by **2.246** for
   fresh catalyst. The steam-reforming data used here need *neither*.
3. **The exponents do not survive OCR.** In the PDF text layer
   $8.664 \times 10^{-7}$ comes out as `8.664 lo-'`. Every parameter below was
   read from a 600 dpi render of the printed page instead, and none was
   repaired by inference — a misread exponent is a wrong rate constant with no
   outward sign. The Table 5 ⇄ Table 6 round trip in the validation section is
   the check on that reading.

**What the paper does not give: the equilibrium constants.** $K_1$, $K_2$ and
$K_3$ appear in the rate equations but are tabulated nowhere in it, so they must
come from elsewhere. This page uses the standard correlations and then tests
them against the paper's own Table 7 reaction enthalpies, which is the strongest
consistency check available without importing a thermodynamic database."""))

cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code('''import sys, urllib.request
from pathlib import Path

if not any("shared" in p for p in sys.path):
    local = Path.cwd()
    for _ in range(4):
        if (local / "shared" / "gallery_utils.py").is_file():
            sys.path.insert(0, str(local / "shared")); break
        local = local.parent
    else:
        url = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
               "pymrm-gallery/main/shared/gallery_utils.py")
        urllib.request.urlretrieve(url, "gallery_utils.py")
        sys.path.insert(0, ".")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymrm import construct_convflux_upwind, construct_div, NumJac, newton
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "C2.1-xu-froment-smr"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

R_KJ = 8.314e-3  # gas constant, kJ/(mol K)'''))

cells.append(md(r"""## The data

Two datasets, both extracted from the paper and both carrying a provenance
sidecar.

The **parameters** are Tables 5 and 6 read off a 600 dpi page render. The
**measurements** are the experimental markers of Figures 2 and 3 — total
methane conversion and conversion into CO₂ against space time, at 10 bar with
H₂O/CH₄ = 3 and H₂/CH₄ = 1.25, at four temperatures. The smooth curves on those
figures are the authors' own model and were deliberately *not* extracted, so
what follows compares this implementation against measurements rather than
against someone else's fit."""))

cells.append(code('''par = load_data("xu-froment-1989-parameters.csv", page=PAGE)
obs = load_data("xu-froment-1989-conversion.csv", page=PAGE)
par_meta = load_meta("xu-froment-1989-parameters.csv", page=PAGE)
obs_meta = load_meta("xu-froment-1989-conversion.csv", page=PAGE)

P = par.set_index("symbol")
print(par[["symbol", "A", "E_or_dH_kJ_per_mol", "value_at_Tref", "T_ref_K"]]
      .to_string(index=False))
print()
print(obs.groupby(["quantity", "temperature_K"]).size().unstack().to_string())
print(f"\\n{len(obs)} measurements, W/F from {obs.W_F.min():.3f} to {obs.W_F.max():.3f}"
      " g_cat h/mol_CH4")
print(cite_data(obs_meta))'''))

cells.append(md(r"""## PyMRM implementation

The reactor is pure convection with a source: $\mathrm{d}(v\,x)/\mathrm{d}\tau
= r$ with $v = 1$ and $\tau = W/F^{0}_{\mathrm{CH_4}}$. That is
`construct_convflux_upwind` followed by `construct_div`, exactly as for a plug
flow in space — the fact that the coordinate is a space time rather than a
length changes nothing about the operators.

The state has the layout the style guide asks for, spatial axis first and
fields last: `(n_tau, 2)` holding $x_{\mathrm{CH_4}}$ and $x_{\mathrm{CO_2}}$.
The source term is pointwise in $\tau$ — the rate at a cell depends only on that
cell's composition — so the default `NumJac(shape)` stencil, which couples the
last axis only, is the right one. Nothing here needs `axes_diagonals`.

Building it this way rather than calling an ODE integrator is what makes the
model reusable: the same residual drops into a non-isothermal tube, or into a
pellet model, or into a bed with an axial dispersion term added to the same
`jac_const`."""))

cells.append(code('''def rate_constants(T, P=P):
    """Arrhenius / van \\'t Hoff form using the Table 6 preexponential factors.

    k_i = A(k_i) exp(-E_i / RT),  K_j = A(K_j) exp(-dH_j / RT).
    The reference temperature is already folded into A, which is why the split
    T_r (648 K for most, 823 K for K_CH4 and K_H2O) does not appear here.
    """
    return {s: P.loc[s, "A"] * np.exp(-P.loc[s, "E_or_dH_kJ_per_mol"] / (R_KJ * T))
            for s in P.index}


def equilibrium(T):
    """K1, K3 in bar^2; K2 dimensionless.

    Xu & Froment do not tabulate these anywhere in the paper, so they have to
    come from outside it. These are the correlations usually quoted with these
    kinetics, from Twigg, *Catalyst Handbook*, 2nd ed. (1989). What this page
    relies on is not the citation but the check in the validation section: the
    van 't Hoff slope of each correlation IS a reaction enthalpy, and the three
    of them reproduce the paper's own Table 7 to within 1.9 %.
    """
    K1 = np.exp(-26830.0 / T + 30.114)
    K2 = np.exp(4400.0 / T - 4.036)
    return K1, K2, K1 * K2


def partial_pressures(x_ch4, x_co2, p_t, sr, hr):
    """Partial pressures from the two conversions, per mol CH4 fed.

    CH4: 1-x   CO: x-y   CO2: y   H2O: sr-x-y   H2: hr+3x+y   total: 1+sr+hr+2x
    """
    n_t = 1.0 + sr + hr + 2.0 * x_ch4
    f = p_t / n_t
    return (f * (1.0 - x_ch4), f * (sr - x_ch4 - x_co2), f * (hr + 3.0 * x_ch4 + x_co2),
            f * (x_ch4 - x_co2), f * x_co2)


def reaction_rates(x_ch4, x_co2, T, p_t, sr, hr, kc=None):
    """The paper's Eqs. 3, in kmol/(kg_cat h)."""
    kc = rate_constants(T) if kc is None else kc
    K1, K2, K3 = equilibrium(T)
    p_ch4, p_h2o, p_h2, p_co, p_co2 = partial_pressures(x_ch4, x_co2, p_t, sr, hr)
    den = (1.0 + kc["K_CO"] * p_co + kc["K_H2"] * p_h2
           + kc["K_CH4"] * p_ch4 + kc["K_H2O"] * p_h2o / p_h2) ** 2
    r1 = kc["k1"] / p_h2 ** 2.5 * (p_ch4 * p_h2o - p_h2 ** 3 * p_co / K1) / den
    r2 = kc["k2"] / p_h2 * (p_co * p_h2o - p_h2 * p_co2 / K2) / den
    r3 = kc["k3"] / p_h2 ** 3.5 * (p_ch4 * p_h2o ** 2 - p_h2 ** 4 * p_co2 / K3) / den
    return r1, r2, r3'''))

cells.append(code('''class SteamReformer:
    """Isothermal, isobaric integral plug flow in the space time tau = W/F0_CH4.

    State layout (n_tau, 2): [..., 0] = x_CH4, [..., 1] = x_CO2.
    """

    def __init__(self, T, p_t=10.0, sr=3.0, hr=1.25, tau_max=0.45, n_tau=1000):
        self.T, self.p_t, self.sr, self.hr = T, p_t, sr, hr
        self.tau_f = np.linspace(0.0, tau_max, n_tau + 1)
        self.tau_c = 0.5 * (self.tau_f[:-1] + self.tau_f[1:])
        self.shape = (n_tau, 2)
        # Outward normal, so both dicts read a*dx/dn + b*x = d.
        # inlet  : x = 0                      -> a=0, b=1, d=0
        # outlet : dx/dn = 0, pure outflow    -> a=1, b=0, d=0
        self.bc = ({"a": 0.0, "b": 1.0, "d": 0.0},
                   {"a": 1.0, "b": 0.0, "d": 0.0})
        self.kc = rate_constants(T)
        self.u = np.zeros(self.shape)
        self._build_operators()

    def _build_operators(self):
        # tau is the flow coordinate and the "velocity" is 1 by construction,
        # so d(v x)/dtau = r. nu=0: the space time is a Cartesian coordinate.
        conv_mat, conv_bc = construct_convflux_upwind(
            self.shape, self.tau_f, self.tau_c, self.bc, v=1.0, axis=0)
        div_mat = construct_div(self.shape, self.tau_f, nu=0, axis=0)
        self.jac_const = div_mat @ conv_mat
        self.g_const = div_mat @ conv_bc
        self.numjac = NumJac(self.shape)   # pointwise source: last axis only

    def reaction(self, u):
        r1, r2, r3 = reaction_rates(u[..., 0], u[..., 1], self.T,
                                    self.p_t, self.sr, self.hr, self.kc)
        return np.stack([r1 + r3, r2 + r3], axis=-1)   # r_CH4, r_CO2

    def residual(self, u):
        g_rxn, jac_rxn = self.numjac(self.reaction, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self, maxfev=60):
        result = newton(self.residual, self.u, maxfev=maxfev)
        self.u = result.x.reshape(self.shape)
        return result

    def at(self, tau, field=0):
        return np.interp(tau, self.tau_c, self.u[:, field])


TEMPS = [773.0, 798.0, 823.0, 848.0]
models = {}
for T in TEMPS:
    m = SteamReformer(T)
    res = m.solve()
    models[T] = m
    print(f"T = {T:.0f} K   converged: {res.success}   Newton iterations: {res.nit}")'''))

cells.append(md("""## Results

The two panels are the paper's Figures 2 and 3 rebuilt: markers are the
digitised measurements, lines are this implementation. Nothing was fitted —
every parameter is Table 6 as printed."""))

cells.append(code('''COL = {773.0: "tab:blue", 798.0: "tab:green", 823.0: "tab:orange", 848.0: "tab:red"}
MRK = {773.0: "o", 798.0: "^", 823.0: "s", 848.0: "D"}

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharex=True)
for ax, (q, field, lab) in zip(
        axes, [("x_CH4", 0, r"$x_{\\mathrm{CH_4}}$ — total methane converted"),
               ("x_CO2", 1, r"$x_{\\mathrm{CO_2}}$ — methane converted into CO$_2$")]):
    for T in TEMPS:
        m = models[T]
        ax.plot(m.tau_c, m.u[:, field], color=COL[T], lw=1.6)
        s = obs[(obs.quantity == q) & (obs.temperature_K == T)]
        ax.plot(s.W_F, s.value, MRK[T], color=COL[T], ms=5,
                mfc="none", mew=1.3, label=f"{T:.0f} K")
    ax.set_xlabel(r"$W/F^{0}_{\\mathrm{CH_4}}$  (g$_\\mathrm{cat}$ h / mol$_{\\mathrm{CH_4}}$)")
    ax.set_ylabel(lab)
    ax.set_xlim(0, 0.45)
    ax.set_ylim(0, None)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9,
           frameon=False, title="markers: measured   ·   lines: pymrm")
fig.suptitle("Xu & Froment (1989) Figures 2 and 3, reproduced from Table 6",
             fontsize=11)
fig.tight_layout(rect=(0, 0.13, 1, 1))
plt.show()'''))

cells.append(code('''def predict(row):
    m = models[float(row.temperature_K)]
    return m.at(row.W_F, 0 if row.quantity == "x_CH4" else 1)


res = obs.copy()
res["model"] = res.apply(predict, axis=1)
res["dev"] = res["model"] - res["value"]

summary = (res.groupby(["quantity", "temperature_K"])["dev"]
             .agg(n="size", MAD=lambda s: s.abs().mean(),
                  worst=lambda s: s.abs().max(), bias="mean").round(5))
print(summary.to_string())
print(f"\\noverall: n = {len(res)}   MAD = {res.dev.abs().mean():.5f}   "
      f"worst = {res.dev.abs().max():.5f}   bias = {res.dev.mean():+.5f}")
print("\\nfor scale, the digitisation error quoted in the provenance sidecar:")
print("  " + " ".join(obs_meta["acquisition"]["estimated_error"].split()))'''))

cells.append(md("""## Validation

Five checks. The first two test the *parameters* — the part of this page that a
bad page render would silently corrupt — and the last three test the *solver*
and the model against the measurements. A sixth section injects defects and
prints what each published metric does with them, because the first check turned
out to have a blind spot that only an injection table can expose.

**The blind spot, stated first.** Check 1's headline used to be
`rt.rel_pct.max()` — the worst per-row round-trip deviation, 1.44 %. That is a
`max()` over seven rows, and its baseline scatter is set by the *printed
precision of the activation energies*: E is given to 0.1 kJ/mol, and at
T_ref = 648 K that alone is worth 0.93 % on `A`. So `k1` and `k3` sit near 1.3
and 1.4 % with nothing wrong at all, and **any single-digit slip that lands
below that ceiling leaves the headline at exactly 1.4432 %.** A `max()` is only
as sharp as its own baseline scatter. The fix is to normalise each row by *its
own* printed precision and to summarise with a statistic that every row reaches:
the sum of the seven ratios."""))

cells.append(code('''# 1. Table 5 -> Table 6 round trip.
#    A = value(T_ref) * exp(E/(R T_ref)) for the rate coefficients, and the same
#    with dH for the adsorption constants. This tests both readings at once, and
#    it tests the split reference temperature: T_ref is part of the formula.
#
#    Each row is scored against the deviation its OWN printed digits already
#    allow. Half a unit in the last printed place of A, of value(T_ref) and of
#    E propagates to
#        u(A)/A + u(value)/value + u(E)/(R T_ref)
#    and the last term dominates for k1 and k3 (E to 0.1 kJ/mol at 648 K is
#    0.93 % on A by itself). Dividing by that bound removes the scatter that
#    makes the raw max() blind.
DIGITS = {"k1":    dict(A="4.225e15", E="240.1", val="1.842e-4"),
          "k2":    dict(A="1.955e6",  E="67.13", val="7.558"),
          "k3":    dict(A="1.020e15", E="243.9", val="2.193e-5"),
          "K_CO":  dict(A="8.23e-5",  E="-70.65", val="40.91"),
          "K_H2":  dict(A="6.12e-9",  E="-82.90", val="0.02960"),
          "K_CH4": dict(A="6.65e-4",  E="-38.28", val="0.1791"),
          "K_H2O": dict(A="1.77e5",   E="88.68",  val="0.4152")}


def half_ulp(s):
    """Half a unit in the last printed place of a number written as a string."""
    m, _, e = s.strip().lower().partition("e")
    dec = len(m.split(".")[1]) if "." in m else 0
    return 0.5 * 10.0 ** (-dec) * 10.0 ** (int(e) if e else 0)


def round_trip(par, digits=DIGITS):
    """Table 5 -> Table 6, each row scored against its own printed precision."""
    rt = par.assign(A_computed=lambda d: d.value_at_Tref
                    * np.exp(d.E_or_dH_kJ_per_mol / (R_KJ * d.T_ref_K)))
    rt["rel_pct"] = (rt.A_computed - rt.A).abs() / rt.A * 100
    rt["bound_pct"] = [(half_ulp(digits[r.symbol]["A"]) / abs(r.A)
                        + half_ulp(digits[r.symbol]["val"]) / abs(r.value_at_Tref)
                        + half_ulp(digits[r.symbol]["E"]) / (R_KJ * r.T_ref_K)) * 100
                       for _, r in rt.iterrows()]
    rt["ratio"] = rt.rel_pct / rt.bound_pct
    return rt


rt = round_trip(par)
print("1. Table 5 -> Table 6 round trip, each row against its own printed precision")
print(rt[["symbol", "T_ref_K", "A", "A_computed", "rel_pct", "bound_pct", "ratio"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
RT_WORST_PCT = float(rt.rel_pct.max())
RT_WORST_RATIO = float(rt.ratio.max())
RT_RATIO_SUM = float(rt.ratio.sum())
print(f"   worst deviation {RT_WORST_PCT:.4f} % (the old headline) — but its "
      "ceiling is set by\\n   the 0.1 kJ/mol precision of E, not by the "
      "transcription.")
print(f"   worst ratio to printed precision {RT_WORST_RATIO:.2f}, "
      f"sum over the seven rows {RT_RATIO_SUM:.3f}.")
print("   Every row lands within about 2.3x its own rounding bound, so nothing is")
print("   mis-transcribed at a level this round trip can resolve. Section 6 shows")
print("   what it can and cannot resolve.")

# The same round trip with 648 K forced on every parameter, to show that the
# split reference temperature is not a detail one can round off.
wrong = par.assign(A_computed=lambda d: d.value_at_Tref
                   * np.exp(d.E_or_dH_kJ_per_mol / (R_KJ * 648.0)))
wrong["ratio"] = wrong.A_computed / wrong.A
bad = wrong[wrong.T_ref_K != 648]
print("\\n   if 648 K were used for every parameter:")
for _, r in bad.iterrows():
    print(f"     {r.symbol:6s} off by a factor of {r.ratio:.3g}")'''))

cells.append(code('''# 2. Equilibrium correlation vs the paper's own thermodynamics.
#    d ln K / d(1/T) = -dH/R, so the correlation slopes ARE reaction enthalpies
#    and can be compared with Table 7 (values at 948 K).
slopes = {"I": 26830.0 * R_KJ, "II": -4400.0 * R_KJ}
slopes["III"] = slopes["I"] + slopes["II"]
table7 = {"I": 224.0, "II": -37.3, "III": 187.5}   # kJ/mol at 948 K, Table 7
print("2. van 't Hoff slopes of K_i vs Xu & Froment Table 7 (948 K)")
worst_eq = 0.0
for k in ("I", "II", "III"):
    rel = abs(slopes[k] - table7[k]) / abs(table7[k]) * 100
    worst_eq = max(worst_eq, rel)
    print(f"   dH_{k:<3s} correlation {slopes[k]:8.2f}   paper {table7[k]:8.2f}"
          f"   kJ/mol   {rel:5.2f} %")
print(f"   worst {worst_eq:.2f} % — the equilibrium constants are consistent with "
      "the paper's own reaction enthalpies,")
print("   and K3 = K1*K2 closes the triangle exactly by construction.")
print("   WHAT THIS CHECK CANNOT SEE: only the SLOPES are compared. The two")
print("   intercepts (30.114 and -4.036) drop out of d ln K / d(1/T) entirely, so")
print("   a +1 slip in either multiplies that K by e = 2.718 at every temperature")
print("   and check 2 does not move at all. Section 6 injects exactly that and")
print("   shows which check does catch it.")'''))

cells.append(code('''# 3. Atom balances along the reactor. The conversions are defined so that C, H
#    and O balance identically; this checks the mole bookkeeping, not the solver.
print("3. Atom balance along the reactor (worst over all four temperatures)")
worst_atoms = {"C": 0.0, "H": 0.0, "O": 0.0}
for T, m in models.items():
    x, y = m.u[:, 0], m.u[:, 1]
    p = np.stack(partial_pressures(x, y, m.p_t, m.sr, m.hr), axis=-1)
    n = p / m.p_t * (1.0 + m.sr + m.hr + 2.0 * x)[:, None]   # mol per mol CH4 fed
    n_ch4, n_h2o, n_h2, n_co, n_co2 = n.T
    worst_atoms["C"] = max(worst_atoms["C"], np.abs(n_ch4 + n_co + n_co2 - 1.0).max())
    worst_atoms["H"] = max(worst_atoms["H"],
                           np.abs(4 * n_ch4 + 2 * n_h2o + 2 * n_h2
                                  - (4 + 2 * m.sr + 2 * m.hr)).max())
    worst_atoms["O"] = max(worst_atoms["O"],
                           np.abs(n_h2o + n_co + 2 * n_co2 - m.sr).max())
for a, v in worst_atoms.items():
    print(f"   {a}: max |imbalance| = {v:.3e} mol per mol CH4 fed")

# and every partial pressure must stay strictly positive: the rate equations
# divide by p_H2, and a negative p_CO would mean more CO2 than converted CH4.
names = ["CH4", "H2O", "H2", "CO", "CO2"]
mins = np.array([[np.min(p) for p in partial_pressures(m.u[:, 0], m.u[:, 1],
                                                       m.p_t, m.sr, m.hr)]
                 for m in models.values()]).min(axis=0)
print("   smallest partial pressure over all four reactors, bar:")
print("     " + "   ".join(f"{n} {v:.3e}" for n, v in zip(names, mins)))
print(f"   all strictly positive: {bool((mins > 0).all())}")'''))

cells.append(code('''# 4. Grid independence. First-order upwind converges as O(h), so the successive
#    differences should halve as the grid doubles.
print("4. Grid independence at 848 K, x_CH4 at tau = 0.40")
prev, vals = None, []
for n in (125, 250, 500, 1000, 2000):
    g = SteamReformer(848.0, n_tau=n)
    g.solve()
    v = g.at(0.40, 0)
    vals.append(v)
    print(f"   n_tau = {n:5d}   x_CH4 = {v:.6f}"
          + ("" if prev is None else f"   change {v - prev:+.2e}"))
    prev = v
richardson = vals[-1] + (vals[-1] - vals[-2])      # O(h): remaining error ~ last step
print(f"   Richardson estimate of the exact value: {richardson:.6f}")
print(f"   discretisation error at n_tau = 1000: {abs(richardson - vals[-2]):.2e}, "
      f"{abs(richardson - vals[-2]) / res.dev.abs().mean() * 100:.0f} % of the "
      "mean deviation from the data")'''))

cells.append(code('''# 5. Against the measurements.
print("5. Deviation from the 61 digitised measurements")
mad, worst = res.dev.abs().mean(), res.dev.abs().max()
print(f"   MAD {mad:.5f}   worst {worst:.5f}   bias {res.dev.mean():+.5f} in conversion")
print(f"   relative to the mean measured conversion ({res.value.mean():.4f}): "
      f"{mad / res.value.mean() * 100:.1f} %")

fig, ax = plt.subplots(figsize=(6.2, 4.4))
for T in TEMPS:
    s = res[res.temperature_K == T]
    ax.plot(s[s.quantity == "x_CH4"].value, s[s.quantity == "x_CH4"].model,
            MRK[T], color=COL[T], ms=6, mfc="none", mew=1.3, label=f"{T:.0f} K")
    ax.plot(s[s.quantity == "x_CO2"].value, s[s.quantity == "x_CO2"].model,
            MRK[T], color=COL[T], ms=6, alpha=0.45)
lim = [0, 1.05 * res.value.max()]
ax.plot(lim, lim, "k-", lw=0.8)
ax.plot(lim, [v + 0.005 for v in lim], "k--", lw=0.6)
ax.plot(lim, [v - 0.005 for v in lim], "k--", lw=0.6)
ax.set(xlim=lim, ylim=lim, xlabel="measured conversion",
       ylabel="pymrm, Table 6 parameters",
       title="Parity — open: $x_{CH_4}$, filled: $x_{CO_2}$; dashed $\\\\pm$0.005")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""### 6. Defect injection — what each published number can and cannot catch

Every metric this page reports needs a row here: one defect it is supposed to
catch, injected, re-run, and printed beside the undamaged value. The error that
actually happens when a table is read off a page render is a **transposed
digit**, not a factor, so that is what is injected — `8.23e-5` read as
`8.32e-5`, `1.020e15` read as `1.002e15`.

Two results are worth reading before the table. The old headline
`param_round_trip_worst_pct` **does not move for any single-digit slip**: three
of them leave it at exactly 1.4432 %, and a transposition in `k3` — the row that
*sets* the maximum — moves it the wrong way, down to 1.2844 %. The sum of
per-row ratios moves for every injection except one, and that one is a change of
0.024 % in a number printed to four significant figures, which no round trip
built from these digits can see."""))

cells.append(code('''print("6a. Table 5 -> Table 6 round trip, under injected transcription defects")
print(f"    {'injected defect':44s}{'worst %':>9}{'worst ratio':>13}{'sum ratios':>12}")


def injected_round_trip(sym, col, new, new_str):
    p2 = par.copy()
    d2 = {k: dict(v) for k, v in DIGITS.items()}
    p2.loc[p2.symbol == sym, col] = new
    d2[sym]["A" if col == "A" else "E" if col.startswith("E_") else "val"] = new_str
    return round_trip(p2, d2)


BREAKS = [
    ("K_CO A 8.23e-5 -> 8.32e-5 (transposition)", "K_CO", "A", 8.32e-5, "8.32e-5"),
    ("k3  A 1.020e15 -> 1.002e15 (transposition)", "k3", "A", 1.002e15, "1.002e15"),
    ("K_H2O A 1.77e5 -> 1.71e5 (transposition)", "K_H2O", "A", 1.71e5, "1.71e5"),
    ("k2  A 1.955e6 -> 1.965e6 (+0.51 %)", "k2", "A", 1.965e6, "1.965e6"),
    ("k1  A 4.225e15 -> 4.226e15 (+0.024 %)", "k1", "A", 4.226e15, "4.226e15"),
    ("k1  E 240.1 -> 240.2 kJ/mol (last digit)", "k1", "E_or_dH_kJ_per_mol", 240.2, "240.2"),
    ("K_CH4 value 0.1791 -> 0.1719 (transposition)", "K_CH4", "value_at_Tref", 0.1719, "0.1719"),
    ("k1  E 240.1 -> 241.1 kJ/mol", "k1", "E_or_dH_kJ_per_mol", 241.1, "241.1"),
    ("k1  A lost decade", "k1", "A", 4.225e14, "4.225e14"),
]
print(f"    {'as published':44s}{RT_WORST_PCT:9.4f}{RT_WORST_RATIO:13.2f}"
      f"{RT_RATIO_SUM:12.3f}")
for lab, sym, col, new, s in BREAKS:
    b = injected_round_trip(sym, col, new, s)
    flag = "" if abs(b.ratio.sum() / RT_RATIO_SUM - 1) > 0.05 else "   <- MISSED"
    print(f"    {lab:44s}{b.rel_pct.max():9.4f}{b.ratio.max():13.2f}"
          f"{b.ratio.sum():12.3f}{flag}")
# a T_ref slip is a different class and the raw max does see it
tr = par.copy(); tr.loc[tr.symbol == "K_CH4", "T_ref_K"] = 648.0
b = round_trip(tr)
print(f"    {'K_CH4 T_ref 823 -> 648 K':44s}{b.rel_pct.max():9.4f}"
      f"{b.ratio.max():13.2f}{b.ratio.sum():12.3f}")
print("    Read the first column: three single-digit slips leave it at 1.4432, and")
print("    the k3 transposition moves it DOWN. The max cannot see a defect smaller")
print("    than its own baseline scatter; the sum of per-row ratios can, because")
print("    every row reaches it. The one miss is a 0.024 % change in a 4-figure")
print("    number -- below the printed precision, so undetectable in principle.")'''))

cells.append(code('''print("\\n6b. The equilibrium correlation. Check 2 compares SLOPES only.")


def eq_slopes(a1=26830.0, a2=4400.0):
    s = {"I": a1 * R_KJ, "II": -a2 * R_KJ}
    s["III"] = s["I"] + s["II"]
    return max(abs(s[k] - table7[k]) / abs(table7[k]) * 100 for k in table7)


print(f"    {'injected defect':44s}{'check 2 worst %':>17}")
print(f"    {'as published':44s}{eq_slopes():17.2f}")
print(f"    {'K1 intercept 30.114 -> 31.114 (K1 x e)':44s}{eq_slopes():17.2f}"
      "   <- MISSED")
print(f"    {'K2 intercept -4.036 -> -3.036 (K2 x e)':44s}{eq_slopes():17.2f}"
      "   <- MISSED")
print(f"    {'K1 slope 26830 -> 26380 (transposition)':44s}{eq_slopes(a1=26380.0):17.2f}")
print(f"    {'K2 slope 4400 -> 4040 (transposition)':44s}{eq_slopes(a2=4040.0):17.2f}")
print("    The intercepts are invisible to check 2 by construction. What sees them")
print("    is check 5 -- the comparison against the 61 measurements -- below.")'''))

cells.append(code('''print("\\n6c. The model against the data, under injected defects.")
print("    Re-solves all four reactors; the deviation metrics are recomputed on the")
print("    same 61 points.")


def metrics_under(param_edits=(), eq=None):
    """(mad, worst, bias, mad_ch4, mad_co2) with parameters/equilibrium perturbed.

    rate_constants closes over the DataFrame P as a default argument, so the
    injection has to mutate that object in place and put it back afterwards.
    """
    saved = [(s, c, P.loc[s, c]) for s, c, _ in param_edits]
    eq_saved = globals()["equilibrium"]
    try:
        for s, c, v in param_edits:
            P.loc[s, c] = v
        if eq is not None:
            globals()["equilibrium"] = eq
        ms = {}
        for T in TEMPS:
            m = SteamReformer(T, n_tau=1000)
            if not m.solve().success:
                return None
            ms[T] = m
        d = obs.copy()
        d["model"] = [ms[float(r.temperature_K)].at(r.W_F, 0 if r.quantity == "x_CH4" else 1)
                      for _, r in d.iterrows()]
        d["dev"] = d.model - d.value
        return (d.dev.abs().mean(), d.dev.abs().max(), d.dev.mean(),
                d[d.quantity == "x_CH4"].dev.abs().mean(),
                d[d.quantity == "x_CO2"].dev.abs().mean())
    finally:
        for s, c, v in saved:
            P.loc[s, c] = v
        globals()["equilibrium"] = eq_saved


base_eq_fn = globals()["equilibrium"]
CASES = [
    ("as published", (), None),
    ("K1 intercept 30.114 -> 31.114 (K1 x e)", (),
     lambda T: (np.exp(-26830.0 / T + 31.114), np.exp(4400.0 / T - 4.036),
                np.exp(-26830.0 / T + 31.114) * np.exp(4400.0 / T - 4.036))),
    ("K2 intercept -4.036 -> -3.036 (K2 x e)", (),
     lambda T: (np.exp(-26830.0 / T + 30.114), np.exp(4400.0 / T - 3.036),
                np.exp(-26830.0 / T + 30.114) * np.exp(4400.0 / T - 3.036))),
    ("K_CO A 8.23e-5 -> 8.32e-5 (transposition)", (("K_CO", "A", 8.32e-5),), None),
    ("k1 A 4.225e15 -> 4.255e15 (transposition)", (("k1", "A", 4.255e15),), None),
    ("k3 A 1.020e15 -> 1.002e15 (transposition)", (("k3", "A", 1.002e15),), None),
    ("K_H2 A 6.12e-9 -> 6.21e-9 (transposition)", (("K_H2", "A", 6.21e-9),), None),
]
print(f"    {'injected defect':44s}{'MAD':>9}{'worst':>9}{'bias':>10}"
      f"{'MAD x_CH4':>11}{'MAD x_CO2':>11}")
ref, eq_factors = None, []
for lab, edits, eq in CASES:
    r = metrics_under(edits, eq)
    if r is None:
        print(f"    {lab:44s}{'did not converge':>50}")
        continue
    if ref is None:
        ref, flag = r, ""
    else:
        flag = "" if max(abs(a / b - 1) for a, b in zip(r, ref)) > 0.05 else "  <- MISSED"
    print(f"    {lab:44s}{r[0]:9.5f}{r[1]:9.5f}{r[2]:+10.5f}{r[3]:11.5f}"
          f"{r[4]:11.5f}{flag}")
    if "intercept" in lab:
        eq_factors.append(r[0] / ref[0])
print(f"    Both equilibrium intercepts, invisible to check 2, are caught here:")
print(f"    the mean deviation grows by {min(eq_factors):.1f}x and "
      f"{max(eq_factors):.1f}x. A transposed digit in an")
print("    adsorption preexponential is NOT caught here -- it is caught by 6a,")
print("    which is why the round trip has to be sharp enough to see it.")'''))

cells.append(code('''# The three activity levels quoted in Reuse. They are transcribed from the
# paper's text; NOTHING on this page validates them, because the runs they
# belong to (Figures 4 and 5, and fresh catalyst) are different chemistry that
# this page does not model. Print what they do so the number is at least
# exercised, and say plainly that it is unchecked.
print("\\n6d. The activity factors quoted in Reuse, exercised but NOT validated")
for f, what in ((1.0, "steam-reforming reference (this page)"),
                (1.225, "paper's Figures 4 and 5"),
                (2.246, "fresh catalyst")):
    r = metrics_under((("k1", "A", 4.225e15 * f), ("k2", "A", 1.955e6 * f),
                       ("k3", "A", 1.020e15 * f)))
    print(f"    k1,k2,k3 x {f:5.3f}  ({what:36s})  MAD vs this page's data "
          f"{r[0]:.5f}")
print("    Only the first row is a validation. The other two are transcriptions")
print("    from the paper's text, printed here so they are at least executed, and")
print("    they are NOT checked against anything -- the data on this page are the")
print("    reference-level runs.")

report_agreement("C2.1", {
    "mad_conversion": mad,
    "max_dev_conversion": worst,
    "bias_conversion": res.dev.mean(),
    "mad_x_ch4": res[res.quantity == "x_CH4"].dev.abs().mean(),
    "mad_x_co2": res[res.quantity == "x_CO2"].dev.abs().mean(),
    # the old headline, kept so the baseline is comparable -- but it is the
    # blind one, and 6a shows three single-digit slips that leave it unmoved.
    "param_round_trip_worst_pct": RT_WORST_PCT,
    "param_round_trip_worst_ratio_to_precision": RT_WORST_RATIO,
    "param_round_trip_ratio_sum": RT_RATIO_SUM,
    "equilibrium_vs_table7_worst_pct": worst_eq,
})'''))

cells.append(md(r"""## What pymrm adds

The paper draws four curves through its points and stops. Two questions it
leaves open can be answered by running the model rather than reading it.

**Where does the residual come from?** The deviation is not uniform — it grows
with space time. That is the signature of an equilibrium-limited residual
rather than a kinetic one, and it matters because the equilibrium constants are
the one ingredient that is *not* from this paper."""))

cells.append(code('''bins = [(0.12, 0.20), (0.20, 0.28), (0.28, 0.42)]
print("bias in conversion, by space time")
for lo, hi in bins:
    s = res[(res.W_F >= lo) & (res.W_F < hi)]
    a = s[s.quantity == "x_CH4"]
    print(f"   W/F {lo:.2f}-{hi:.2f}:  n = {len(a):2d}   "
          f"bias(x_CH4) = {a.dev.mean():+.5f}")

# How far from equilibrium is the gas? beta_i = Q_i / K_i, so beta -> 1 at
# equilibrium and the rate of reaction i vanishes.
def beta(m, tau):
    x, y = m.at(tau, 0), m.at(tau, 1)
    p_ch4, p_h2o, p_h2, p_co, p_co2 = partial_pressures(x, y, m.p_t, m.sr, m.hr)
    K1, K2, K3 = equilibrium(m.T)
    return (p_h2 ** 3 * p_co / (p_ch4 * p_h2o * K1),
            p_h2 * p_co2 / (p_co * p_h2o * K2),
            p_h2 ** 4 * p_co2 / (p_ch4 * p_h2o ** 2 * K3))


print("\\napproach to equilibrium at 848 K (beta = Q/K; 1 = equilibrated)")
for tau in (0.13, 0.20, 0.28, 0.35, 0.45):
    b = beta(models[848.0], tau)
    print(f"   tau = {tau:.2f}:  I {b[0]:.3f}   II {b[1]:.3f}   III {b[2]:.3f}")

# Sensitivity: shift all three K_i together and see what it does at long tau.
base_eq = equilibrium
sens = {}
for f in (0.95, 1.00, 1.05):
    globals()["equilibrium"] = lambda T, f=f: tuple(f * k for k in base_eq(T))
    sens[f] = {T: SteamReformer(T, n_tau=500) for T in (823.0, 848.0)}
    for s in sens[f].values():
        s.solve()
globals()["equilibrium"] = base_eq

far = res[(res.W_F >= 0.28) & (res.quantity == "x_CH4")]
print("\\nsensitivity of x_CH4 at tau = 0.35 to a 5 % shift in K_i, against the")
print("observed bias over the same space times (W/F > 0.28):")
for T in (823.0, 848.0):
    lo, mid, hi = (sens[f][T].at(0.35, 0) for f in (0.95, 1.00, 1.05))
    b = far[far.temperature_K == T].dev.mean()
    print(f"   {T:.0f} K:  K x0.95 {lo:.4f}   base {mid:.4f}   K x1.05 {hi:.4f}"
          f"   |half-span| {abs(hi - lo) / 2:.4f}   observed bias {b:+.4f}")'''))

cells.append(md(r"""A 5 % change in the equilibrium constants moves the
long-space-time predictions by 0.003–0.004 in conversion — the same size as the
residual bias there, and in the same direction. Since Xu and Froment do not
report the $K_i$ they used, the residual at long space time cannot be charged
to the kinetics: it is inside the uncertainty of an ingredient the paper does
not supply. That is a different, and more useful, statement than "agreement is
good". At short space time, where the reactions are still far from equilibrium
and the kinetics really are being tested, the bias is a third as large.

**Where does each reaction stop being kinetically controlled?** The $\beta$
values above split the reactor into regimes: the water-gas shift equilibrates
first and well before the reforming reactions do, which is exactly the
behaviour the paper asserts from its $V(\mathrm{II})$ plot but never shows
against space time. Below, the same thing across all four temperatures."""))

cells.append(code('''fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
tau = np.linspace(0.005, 0.45, 200)
for T in TEMPS:
    m = models[T]
    b = np.array([beta(m, t) for t in tau])
    axes[0].plot(tau, b[:, 0], color=COL[T], lw=1.6, label=f"{T:.0f} K")
    axes[0].plot(tau, b[:, 2], color=COL[T], lw=1.0, ls=":")
    axes[1].plot(tau, b[:, 1], color=COL[T], lw=1.6, label=f"{T:.0f} K")
for ax, ttl in zip(axes, ["reforming: I (solid), III (dotted)",
                          "water-gas shift: II"]):
    ax.axhline(1.0, color="k", lw=0.8)
    ax.axhline(0.95, color="k", lw=0.6, ls="--")
    ax.set(xlabel=r"$W/F^{0}_{\\mathrm{CH_4}}$", xlim=(0, 0.45), ylim=(0, 1.05),
           title=ttl)
axes[0].set_ylabel(r"$\\beta = Q/K$  (1 = equilibrated)")
axes[0].legend(fontsize=8, loc="lower right")
fig.suptitle("Where each reaction stops being kinetically controlled", fontsize=11)
fig.tight_layout()
plt.show()

print("space time (g_cat h/mol_CH4) at which each reaction reaches beta = 0.90,")
print("and beta at the end of the measured range:")
for T in TEMPS:
    m = models[T]
    b = np.array([beta(m, t) for t in tau])
    hit = [f"{tau[np.argmax(b[:, i] > 0.90)]:.3f}" if (b[:, i] > 0.90).any() else ">0.45 "
           for i in (0, 1, 2)]
    e = beta(m, 0.45)
    print(f"   {T:.0f} K:  I {hit[0]}  II {hit[1]}  III {hit[2]}"
          f"   |   beta(0.45) = {e[0]:.3f}, {e[1]:.3f}, {e[2]:.3f}")
print("\\nThe shift (II) is ahead of the reforming reactions at every temperature:")
for T in TEMPS:
    e = beta(models[T], 0.45)
    print(f"   {T:.0f} K:  beta_II - beta_I = {e[1] - e[0]:+.3f}")'''))

cells.append(md(r"""## Reuse

**The kinetics on their own.** `reaction_rates(x_CH4, x_CO2, T, p_t, sr, hr)`
returns $r_1, r_2, r_3$ in kmol/(kg$_\mathrm{cat}$ h) and has no dependence on
the reactor model. `partial_pressures` is the only piece tied to the
single-feed parameterisation; swap it for a general mole-fraction vector and
the rate equations carry over unchanged.

**Three activity levels — and only one of them is tested here.** Multiply
$k_1, k_2, k_3$ by

- **1** for the steam-reforming reference level — what this page uses, what
  Tables 5 and 6 are, and the only level the 61 measurements on this page
  belong to;
- **1.225** to reproduce the paper's Figures 4 and 5 (reverse water-gas shift
  and methanation);
- **2.246** for fresh catalyst.

The last two are **transcriptions from the paper's text that nothing on this
page validates.** Section 6d executes them so they are at least run, and shows
what they do to the fit against *these* data — which is to make it worse, as it
must, because these data are the reference-level runs. Reproducing Figures 4
and 5 would need the reverse water-gas shift and methanation feeds, which this
page does not model. Treat both factors as unchecked.

**The round trip is the check on the parameter reading, and it has a floor.**
Section 6a measures it: a slip smaller than a row's own printed precision is
invisible, and a `max()` over the seven rows is invisible to anything smaller
than the 1.4 % baseline scatter that the 0.1 kJ/mol precision of $E$ already
produces. If you re-read these tables from your own render, compare
`param_round_trip_ratio_sum`, not the worst percentage.

**Going further with the same residual.** `SteamReformer.residual` is the whole
model. To make the tube non-isothermal, add an enthalpy field to the state and
a wall-flux term to `g_const`; to resolve the pellet, keep this as the bulk
equation and couple a `B1.1`-style intraparticle model to it — that is `D1.4`,
and the structure codes are how to find it.

**Related pages.** `B1.1` (intraparticle diffusion and reaction, the resistance
these intrinsic kinetics deliberately exclude), `D1.4` (fixed bed with resolved
particles), `C2.2`, `C2.6`, `D3.1`.

**Cite the source, not this page**, for the kinetics: Xu, J. and Froment, G. F.,
*Methane steam reforming, methanation and water-gas shift: I. Intrinsic
kinetics*, AIChE Journal **35**(1) 88–96 (1989),
[doi:10.1002/aic.690350109](https://doi.org/10.1002/aic.690350109)."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py", "mimetype": "text/x-python", "name": "python",
        "nbconvert_exporter": "python", "pygments_lexer": "ipython3",
        "version": "3.13.5"},
}
out = Path(__file__).with_name("index.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
