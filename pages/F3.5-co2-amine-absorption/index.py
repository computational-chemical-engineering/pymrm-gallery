# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ---
# title: "CO2 absorption into amine-promoted carbonate solutions"
# description: "Seven species, two kinetic reactions and two acid-base equilibria decide how fast CO2 enters hot potash with an amine promoter. This page rebuilds the film model that settled the shuttle-versus-catalysis debate."
# categories: [sec:F, struct:S3, struct:S10, tier:T1, data:tier6, phase:gas-liquid]
# date: 2026-08-02
# ---
#
# # CO2 absorption into amine-promoted carbonate solutions
#
# **Catalog ID:** `F3.5` · **Structures:** `S3` (1D steady BVP) + `S10`
# (instantaneous-equilibrium constraints) · **Tier:** T1
#
# Add a few mol% of an amine to a hot potassium carbonate solution and CO2
# absorbs several times faster. For twenty years the field argued about *why*:
# a "shuttle" carrying CO2 into the bulk as carbamate, or homogeneous catalysis
# of the hydration reaction, each argued with an approximate model that
# captured one mechanism and excluded the other. Bosch, Versteeg and van Swaaij
# ended the argument by refusing to approximate: solve the diffusion-reaction
# film with **every** reaction in it, reversible, simultaneously, and ask the
# solution what it is doing.
#
# That model is a seven-species boundary-value problem with two kinetic
# reactions and two acid-base equilibria, the shape of problem pymrm's
# operators are built for. This page reproduces their Table 2 (twelve
# enhancement factors, twelve fluxes, eight promotion factors), their desorption
# promotion factor and their concentration profiles, and in doing so finds an
# inconsistency between the printed equilibrium constants and the ones their
# code must have used. One rescaled constant, fitted to the paper's own numbers,
# accounts for it; the page states exactly what that fit consumes, where the
# alternatives fail, and where the paper's own record argues *against* the
# reconstruction.

# %% [markdown]
# ## Background
#
# Hot potash (aqueous K2CO3) was the workhorse solvent for bulk CO2 removal:
# ammonia plants, hydrogen plants, natural gas sweetening. Absorption is
# enhanced by the hydration reaction CO2 + OH⁻ → HCO3⁻, but not by much, so
# small amounts of amine are added as **rate promoters**. By 1984 the
# experiments of Savage, Sartori and Astarita (Exxon) had shown that
# diethanolamine (DEA) raises the absorption flux about four-fold and a
# sterically hindered diamine ("HDA", structure never disclosed) about six-fold,
# and two mechanistic pictures competed:
#
# - **Shuttle** (Danckwerts & McNeil 1967; Astarita et al. 1981): the amine
#   reacts with CO2 near the interface, the carbamate diffuses into the bulk,
#   reverts, and the amine returns for more, a ferry for CO2.
# - **Homogeneous catalysis** (Savage et al. 1984): the amine merely catalyses
#   CO2 hydration, spent and regenerated at the same place.
#
# Each camp had an analytical model built on its own limiting assumptions.
# Bosch, Versteeg and van Swaaij applied the rigorous numerical method of
# Versteeg, Kuipers, van Beckum and van Swaaij (1989) to the full reversible
# network and found that **neither picture survives**: the computed profiles
# show shuttle transport *and* an elevated free-CO2 level that accelerates
# hydration, mechanisms that only coexist in a model that keeps all reactions.
# Their simulations bracket the measured promotion factors using independently
# determined parameters only, and explain why tertiary amines do not promote
# (no carbamate, so neither mechanism is available to them).
#
# **The reference.** The gallery catalogue cited this case as "Versteeg & van
# Swaaij (1988)", which is the same group's *kinetics* paper, the source of the
# k2 correlation this paper uses as its eq. (11). The model the case wants,
# absorption with parallel reversible reactions in CO2-amine systems, is
# carried by the paper reproduced here: Bosch, Versteeg & van Swaaij,
# *Gas-liquid mass transfer with parallel reversible reactions, II. Absorption
# of CO2 into amine-promoted carbonate solutions*, Chem. Eng. Sci. 44(11)
# 2735-2743 (1989), doi:10.1016/0009-2509(89)85216-9.

# %% [markdown]
# ## The published model
#
# **Reaction network** (their reactions 1 to 5), for a secondary amine R2NH:
#
# $$\mathrm{CO_2 + OH^- \rightleftharpoons HCO_3^-} \tag{1}$$
# $$\mathrm{CO_2 + 2\,R_2NH \rightleftharpoons R_2NCOO^- + R_2NH_2^+} \tag{2}$$
# $$\mathrm{CO_3^{2-} + H_2O \rightleftharpoons HCO_3^- + OH^-} \tag{3}$$
# $$\mathrm{R_2NH_2^+ + H_2O \rightleftharpoons R_2NH + H_3O^+} \tag{4}$$
# $$\mathrm{2\,H_2O \rightleftharpoons H_3O^+ + OH^-} \tag{5}$$
#
# Reactions (1) and (2) are kinetically controlled; the proton transfers (3)
# to (5) are instantaneous with respect to mass transfer and stay at
# equilibrium. Reaction (2) proceeds through the zwitterion mechanism with
# zwitterion formation rate-determining at low loading (their reading of Tseng
# et al. 1988), so its forward rate is $k_2[\mathrm{CO_2}][\mathrm{Am}]$; the
# reverse rate is fixed by the overall equilibrium constant, in the paper's
# notation
#
# $$K_{eq} = \frac{[\mathrm{AmCOO^-}][\mathrm{AmH^+}]}{[\mathrm{CO_2}][\mathrm{Am}]^2}
#        = \frac{K_c\,K_{c1}\,K_w}{K_p},$$
#
# with $K_c$ the carbamate stability constant, $K_{c1}$ the CO2-hydrolysis
# equilibrium, $K_w$ water's ion product and $K_p$ the amine protonation
# constant, all defined in their Notation and valued in their Table 1.
#
# **Film model.** A stagnant film of thickness $\delta = D_{CO_2}/k_L$
# separates the interface from a well-mixed bulk at chemical equilibrium:
#
# $$D_j \frac{\mathrm{d}^2 c_j}{\mathrm{d}x^2} + \sum_r \nu_{jr}\,R_r = 0,
# \qquad 0 \le x \le \delta,$$
#
# for $j \in$ {CO2, HCO3⁻, CO3²⁻, OH⁻, Am, AmH⁺, AmCOO⁻}. At the interface only
# CO2 crosses (ions and amine are non-volatile, zero flux); with the paper's
# "no gas phase resistance", the interfacial CO2 is $c_i = m\,p_{CO_2}/(RT)$.
# At $x = \delta$ every species takes its bulk equilibrium value. Reported
# quantities:
#
# $$J = -D_{CO_2}\left.\frac{\mathrm{d}c_{CO_2}}{\mathrm{d}x}\right|_{x=0},
# \qquad
# E = \frac{J}{k_L\,(c_i - c_b)},
# \qquad
# F = \frac{J_\text{promoted}}{J_\text{unpromoted}}.$$
#
# The paper argues (its section 5.3) that the approximate promotion-factor
# expressions, the shuttle ceiling (their eq. 14) and the
# instantaneous-equilibrium result of Bosch et al. 1989b (their eq. 15),
# $F = \sqrt{1 + [\mathrm{Am}]^2 K_{eq}/[\mathrm{AmH^+}]}$, both fail against
# the full solution; eq. (15) "gives much too high values". That claim is
# checked below.
#
# **What is simplified here, and why it is safe.** H3O⁺ is eliminated
# analytically: at the bulk pH of these solutions it is about ten orders below
# the major ions, so reactions (4)+(5) are combined into AmH⁺ + OH⁻ ⇌ Am + H2O
# with constant $K_w/K_p$, and (3) is used directly with constant
# $K_3 = K_w/K_{c2}$. The instantaneous equilibria are imposed as fast
# reversible kinetics with a rate constant $k_{inst}$, a model parameter whose
# value is the last rung of `model.K_LADDER`; the effect of a further decade is
# measured under Validation (`kinst_delta`).

# %% [markdown]
# ## Parameters and assumptions
#
# All parameters are Table 1 of the paper, read from a 600 dpi page render
# because the scan's OCR text layer mangles exponents. The transcription is
# cross-checked in the data sidecar: $K_w$, $K_{c2}$ and $k_{OH}$ agree with
# independent literature values at the paper's temperature, and $k_2$(HDA) is
# consistent with the paper's own correlation (11).
#
# **Assumptions carried over from the paper:** film model at the paper's
# temperature; ideal solution (their statement: "non-idealities in the liquid
# phase equilibria were not taken into account"); a single diffusion
# coefficient for all ions, which makes the migration term vanish identically,
# so Fick suffices (checked below); loading $\alpha$ = mol CO2 absorbed per mol
# K2CO3, so that $[\mathrm{HCO_3^-}] = 2\alpha c_{carb}$ and
# $[\mathrm{CO_3^{2-}}] = (1-\alpha)c_{carb}$ in the unpromoted bulk.
#
# **The paper is internally inconsistent about $\alpha$, and this page follows
# Table 2.** The convention just stated is confirmed by Figure 4, whose
# unpromoted bulk values agree with the speciation at the loading printed in
# its inset. Figures 2 and 3 do not: their unpromoted bulk is the composition
# at half the loading their insets print, and only that reading is
# electroneutral against $[\mathrm{K^+}] = 2c_{carb}$. Table 2 follows Figure 4:
# re-reading Table 2's loadings on the Figure 2/3 convention makes the
# consistency test of the next section collapse, and the scalar it then
# demands (printed below) is several times larger and far more scattered. The
# most economical explanation is a misprint in the Figure 2/3 insets; the
# figure readings behind this paragraph are set out, with their review status,
# further down the page.
#
# **Assumptions added here (both measured):** the gas-side Robin condition is
# imposed as a Dirichlet value, because the gas-side fraction of the total
# resistance is $m\,E\,k_L/k_g$ (from $1/K_g = 1/k_g + 1/(m E k_L)$), computed
# under Validation at the largest $E$ of the reconstructed run
# (`gas_side_fraction`); and H3O⁺ is eliminated as described above. The paper
# itself flags its softest inputs: the concentrations $c_{Am}$ and $c_{carb}$
# "had to be guessed", and everything about HDA's equilibrium behaviour is
# uncertain because its structural formula was never published.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.27.8"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "F3.5-co2-amine-absorption"  # equals the page directory

import os, sys, urllib.request
from pathlib import Path
_shared = next((p / "shared" for p in (Path.cwd(), *Path.cwd().parents)
                if (p / "shared" / "gallery_bootstrap.py").is_file()), None)
if _shared is None:
    _shared = (Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
               / "pymrm-gallery" / GALLERY_REF)
    _shared.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/computational-chemical-engineering/"
        f"pymrm-gallery/{GALLERY_REF}/shared/gallery_bootstrap.py",
        _shared / "gallery_bootstrap.py")
sys.path.insert(0, str(_shared))
sys.modules.pop("gallery_bootstrap", None)
from gallery_bootstrap import setup
setup(PAGE, PYMRM, GALLERY_REF);  # the semicolon stops the cell echoing the path

# %%
# Imports. Packages beyond pymrm and the Colab defaults (numpy, scipy, pandas,
# matplotlib, pyyaml) are installed here with subprocess, not %pip.
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymrm.checks import find_roots, observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims
import model as M          # the page's model.py, put on sys.path by setup()
from model import (SP, Film, Penetration, bulk, c_interface, run_case, params,
                   closures, ceiling, pfo_enhancement, vkh, second_order_dimensionless)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# fixed series colours, assigned once and used everywhere on this page
COL = {"u": "tab:gray", "DEA": "tab:blue", "HDA": "tab:orange"}

# %%
# ---- Table 1: every constant in the model, as printed ---------------------
t1 = load_data("bosch1989-table1-parameters.csv", page=PAGE)
print(cite_data(load_meta("bosch1989-table1-parameters.csv", page=PAGE)))
P = params(t1)      # one dictionary; every study takes it, so a defect is one entry
DELTA = P["DCO2"] / P["KL"]                    # film thickness delta = D_CO2/kL, m
K3 = P["KW"] / P["KC2"]                        # [HCO3][OH]/[CO3], reaction (3)
KEQ = {a: P["KCC"][a] * P["KC1"] * P["KW"] / P["KP"][a] for a in ("DEA", "HDA")}


def cgas(p_kpa, P=P):
    """CO2 concentration in the gas, p/(RT), mol/m3."""
    return p_kpa * 1e3 / (P["R_GAS"] * P["T"])


print(t1.to_string(index=False))
print(f"\nderived: delta = {DELTA*1e6:.2f} um   K3 = {K3:.3f} mol/m3   "
      f"Keq(DEA) = {KEQ['DEA']:.3e}, Keq(HDA) = {KEQ['HDA']:.3e} m3/mol")

# %% [markdown]
# ## The data
#
# **This page is provenance tier 6: reproduction, not validation.** The paper
# tabulates no measurements. Its Table 2 contains the *authors' own numerical
# results*, twelve (flux, enhancement) pairs plus eight promotion factors over
# four (loading, pressure) conditions, and its Table 3 compares the resulting
# promotion-factor ranges with the experimental values of Savage et al. (1984),
# whose glyphs print as **"±4"** (DEA) and **"±6"** (HDA) in this scan and are
# almost certainly "≈4" and "≈6"; they are used only in that weak sense
# throughout. The text adds two statements that are **not** in Table 2, the
# desorption promotion factor at the conditions of its Figure 4 and that the
# amine lowers the desorption driving force "by about one third", and those two
# are the only out-of-sample tests available. They are run, both ways, in the
# desorption section, where they turn out to argue *against* the
# reconstruction. Reproducing Table 2 tests this implementation against their
# numerics; the comparison with experiment is inherited from the paper and is
# only as sharp as "about 4" and "about 6".
#
# Two identities make the transcription self-checking, and the second does real
# work later:
#
# 1. each printed F must equal the ratio of printed fluxes, $F = J_p/J_u$;
# 2. each printed (J, E) pair implies the bulk dissolved CO2 the authors used,
#    $c_b = c_i - J/(E\,k_L)$, a model-independent statement, since it is just
#    the definition of E.

# %%
# ---- Table 2: the paper's numerical results + transcription checks --------
t2 = load_data("bosch1989-table2-results.csv", page=PAGE)
print(cite_data(load_meta("bosch1989-table2-results.csv", page=PAGE)))
print(t2.to_string(index=False), "\n")
CONDITIONS = [(float(a), float(p)) for a, p in zip(t2.alpha, t2.p_co2_kpa)]


def implied_cb(t2=t2, P=P):
    """The bulk CO2 each printed (J, E) pair implies, c_b = c_i - J/(E kL), mol/m3."""
    ci = c_interface(P, t2.p_co2_kpa.to_numpy())
    return {k: ci - t2[f"J_{k}"].to_numpy() / (t2[f"E_{k}"].to_numpy() * P["KL"])
            for k in ("u", "DEA", "HDA")}


def study_transcription(t2=t2):
    """Check 1: printed F against the ratio of printed fluxes."""
    f_dev = np.concatenate([t2.F_DEA / (t2.J_DEA / t2.J_u) - 1,
                            t2.F_HDA / (t2.J_HDA / t2.J_u) - 1])
    return {"table2_F_vs_Jratio_maxdev": float(np.max(np.abs(f_dev)))}


TR = study_transcription()
print(f"printed F vs printed J-ratios: max |dev| = "
      f"{TR['table2_F_vs_Jratio_maxdev']:.2%} (print rounding only)")
print(f"this ties {4 * 5} of the {4 * 8} printed numbers together; the "
      f"{4 * 3} it leaves untested\nare the twelve enhancement factors, which the "
      "reconstruction below leans on.")

# check 2: the implied bulk CO2
CB_IMPL = implied_cb()
# the two alpha = 0.4 rows overdetermine c_b without the definition of E,
# only that both rows share one bulk state:
r1, r2 = t2.iloc[1], t2.iloc[2]           # the two rows at alpha = 0.4
_k = (r2.J_u / r1.J_u) * (r1.E_u / r2.E_u)
_ci = c_interface(P, np.array([r1.p_co2_kpa, r2.p_co2_kpa]))
CB_TWO_ROW = float((_ci[1] - _k * _ci[0]) / (1 - _k))
print("\nimplied unpromoted bulk CO2 (mol/m3): "
      + ", ".join(f"alpha={a}: {c:.3f}" for (a, _), c in zip(CONDITIONS, CB_IMPL["u"])))
print(f"two-row consistency at alpha=0.4 (no definition of E used): c_b = {CB_TWO_ROW:.3f}")

# %% [markdown]
# ## PyMRM implementation
#
# Two solvers, sharing the constants above, both in `model.py`.
#
# **Bulk speciation** (`bulk`): the well-mixed liquid at loading α satisfies
# the four equilibria of reactions (1) to (5) plus three balances: total amine,
# total carbon $c_{carb}(1+\alpha)$, and electroneutrality with
# $[K^+] = 2c_{carb}$. Solved in log-concentrations with `scipy.optimize.root`.
#
# **Film** (`Film`): the seven-species BVP on ξ = x/δ with a geometrically
# stretched grid clustered at the interface (the reaction layer is about δ/Ha,
# a small fraction of the film for HDA). Operators are assembled **once**:
# `construct_grad` with the boundary conditions built in, `construct_div` with
# `nu=0` (slab). Both boundary conditions are written on the **outward
# normal** (`a·∂c/∂n + b·c = d`), each with its physical equation in a comment
# of `Film.__init__`:
#
# - interface, ξ=0 (n = −x): CO2 Dirichlet `{a:0, b:1, d: m c_g}` (the gas-side
#   resistance is negligible, measured under Validation); all other species
#   zero-flux `{a:1, b:0, d:0}`, non-volatile;
# - bulk, ξ=1 (n = +x): Dirichlet `{a:0, b:1, d: c_bulk,j}` for every species.
#
# The pointwise reaction source couples only the species axis, so
# `NumJac((n_x, n_c))`, spatial axis first, with its default last-axis stencil
# is the right Jacobian structure. The steady state is found with
# `pymrm.newton(..., bounds=...)` (positivity), with `success` checked. The
# instantaneous equilibria make the system stiff: from the flat bulk start
# Newton does not reach the solution at the model's $k_{inst}$ directly, so the
# solve is continued along a fixed ladder of $k_{inst}$ (`model.K_LADDER`),
# each rung starting from the converged profile of the one below. The flux is
# read from the face-flux operator on the interface face, boundary term
# included, so it is the flux the solve transports.

# %%
# ---- bulk speciation ------------------------------------------------------
b06, b02 = bulk(P, 0.6), bulk(P, 0.2)
print(f"unpromoted bulk at alpha=0.6: HCO3 = {b06['HCO3']:.0f}, CO3 = {b06['CO3']:.0f}"
      f" mol/m3 (compared with Figure 4 in the figure section)")
print(f"equilibrium residual {b06['res']:.1e}; bulk CO2 = {b06['CO2']:.3f} mol/m3")
print(f"unpromoted bulk at alpha=0.2: HCO3 = {b02['HCO3']:.0f}, CO3 = {b02['CO3']:.0f}"
      " mol/m3, the composition Figures 2 and 3\n  plot, although their insets "
      "print alpha = 0.4.")
print("film species order:", SP, "; k_inst ladder:", M.K_LADDER)

# %% [markdown]
# ## Results
#
# ### First pass: every constant exactly as printed
#
# Deviation convention here and everywhere below: **(model − paper)/paper**.

# %%
E_COLS = ["E_u", "E_DEA", "E_HDA"]


def table2_comparison(P=P, s=1.0, t2=t2, amines=(None, "DEA", "HDA")):
    """Film solves at the four Table 2 conditions; values and deviations."""
    rows, films = [], {}
    for alpha, p in CONDITIONS:
        fs = {am: run_case(P, alpha, p, am, s) for am in amines}
        films[(alpha, p)] = fs
        row = {"alpha": alpha, "p_kpa": p}
        for am, f in fs.items():
            row[f"E_{am or 'u'}"], row[f"J_{am or 'u'}"] = f.enhancement(), f.J_co2()
        for am in amines:
            if am:
                row[f"F_{am}"] = fs[am].J_co2() / fs[None].J_co2()
        rows.append(row)
    df = pd.DataFrame(rows)
    dev = df[["alpha", "p_kpa"]].copy()
    for col in df.columns[2:]:
        dev[col] = df[col].to_numpy() / t2[col].to_numpy() - 1.0
    return df, dev, films


def study_printed(P=P):
    """Table 2 with every constant as printed."""
    df, dev, films = table2_comparison(P, 1.0)
    return {"E_maxdev_printed": float(np.abs(dev[E_COLS].to_numpy()).max()),
            "_df": df, "_dev": dev}


PR = study_printed()
print("model (constants as printed):")
print(PR["_df"].round(3).to_string(index=False))
print("\n(model - paper)/paper:")
print(PR["_dev"].round(3).to_string(index=False))
_e = PR["_dev"][E_COLS].to_numpy()
_j = PR["_dev"][["J_u", "J_DEA", "J_HDA"]].to_numpy()
print(f"\nE deviations, printed constants: {int((_e > 0).sum())} of {_e.size} high, "
      f"by {_e.min():+.1%} to {_e.max():+.1%}\n(mean {_e.mean():+.1%}); "
      f"the fluxes are {_j.min():+.1%} to {_j.max():+.1%}.")

# %% [markdown]
# ### The offset is systematic, and the paper's own numbers locate it
#
# Every enhancement factor comes out high, both amines and none, by the margin
# the cell above prints. A single systematic cause is more likely than twelve
# coincidences, and the paper provides the instrument to find it: the implied
# bulk CO2 from identity (2) above.
#
# With the constants as printed, the speciation gives a bulk dissolved CO2
# well **below** what Table 2 itself implies, consistently at every loading,
# and definition-free at α = 0.4, where the two printed rows force the value
# printed under Table 2 against a much lower printed-constants value (the cell
# below prints both). The paper's own Figure 4 corroborates the higher value
# independently (its dashed, unpromoted CO2 profile, read in the figure section).
#
# A single equilibrium constant can carry the whole discrepancy. The bulk CO2
# at fixed loading scales as $c_b \propto 1/(K_{c1} K_3) = K_{c2}/(K_{c1}
# K_w)$; if the offending factor sits in $K_3 = K_w/K_{c2}$, the *same*
# rescaling also lowers the bulk OH⁻ (through $[\mathrm{OH^-}] =
# K_3[\mathrm{CO_3^{2-}}]/[\mathrm{HCO_3^-}]$) and thereby cuts the Hatta
# number of the hydration reaction by the square root of the scalar, the size
# of the observed enhancement offset. Placing the factor in $K_{c1}$ instead
# would fix $c_b$ but leave OH⁻, and the offset, untouched (that alternative is
# run below and misses $E_u$, `Eu_max_branch_Kc1`).
#
# **The location is a three-way choice, not a two-way one.** Because
# $K_3 = K_w/K_{c2}$, dividing $K_w$ by the scalar is indistinguishable from
# multiplying $K_{c2}$ by it for *every unpromoted* quantity. The two differ
# only through the amine equilibria, since $K_w$ also enters $K_B = K_w/K_p$
# and $K_{eq} = K_c K_{c1} K_w/K_p$. **The code implements the $K_{c2}$
# branch** ($K_B$ and $K_{eq}$ untouched); the $K_w$ branch is run below for
# comparison, and the data barely separate them. So: the *magnitude* of the
# correction is well determined; the *$K_3$-versus-$K_{c1}$* choice is decided
# by the sign and size of the E offset; the *$K_{c2}$-versus-$K_w$* split
# inside $K_3$ is essentially undetermined by these data.
#
# **Why $K_3$ differs from its printed parts is not determined by anything on
# this page.** An ionic-strength effect has the right sign (the ionic strength
# here is several molar, and the paper cites Savage et al. (1980) for exactly
# such a correction to $k_{OH}$), but the magnitude is not derived: the scalar
# is a small shift in $pK_{c2}$ for that ionic strength (the fit cell prints
# it), and the paper itself says non-idealities "were not taken into account".
# A competing clerical explanation fits at least as well: on the $K_w$ branch
# the required $pK_w$ (also printed) corresponds to a reference temperature
# some kelvin off on the standard $pK_w(T)$ curve. Both are candidate
# explanations; neither is established here.
#
# ### What the scalar is fitted to, stated exactly
#
# The scalar is fitted to the four unpromoted **implied** $c_b$ values: it is
# the root of the mean log ratio of implied to modelled $c_b$, found by
# `pymrm.checks.find_roots`. Those values are computed as
# $c_b = c_i - J/(E\,k_L)$, so they are built from the printed $J$ *and* the
# printed $E$: the fit consumes a particular combination of the two, and it
# would be false to say the enhancement factors play no part. What is true is
# narrower and is what the next cells show:
#
# - The identity does not admit one scalar exactly. Reconciling it condition by
#   condition demands four different scalars, printed below with their spread
#   (`k3_scale_spread_per_condition`, n = 4 conditions). One scalar is a
#   compromise, and the residual of that compromise has to land somewhere.
# - It lands almost entirely in the **flux**. $E$ carries the fitted $c_b$ in
#   its denominator, so an error in $c_b$ partially cancels between numerator
#   and denominator. The cell below prints the two local exponents,
#   $E \sim s^{a}$ and $J \sim s^{b}$, measured by re-running the model on
#   either side of the fit: $E$ is the less sensitive. The four $E_u$
#   deviations are therefore nearly one constant offset repeated four times,
#   not four independent agreements. **The reproduction metric is the flux,
#   and it is reported as the headline below.**
# - Refitting the same scalar to $E_u$ **alone** (the root of the mean $E_u$
#   deviation) gives a different number, also printed below. The two fitting
#   routes agreeing to about a percent is the real content of the sub-percent
#   $E_u$ figure, not sub-percent accuracy of the film model.
#
# Its *location* ($K_3$ rather than $K_{c1}$) is a discrete choice among the
# three named above, and making it does consult the sign and size of the E
# offset once, the same move as `F2.3`'s discrete choice between printed
# alternatives. What is then predictive, with nothing further adjustable: the
# variation of $E_u$ and $J_u$ across loading and pressure, all eight promoted
# enhancement and promotion factors, the desorption case, and the bulk
# compositions the three figures show.

# %%
# ---- fit the single reconstruction scalar s (K3 -> K3/s, i.e. K_c2 -> s K_c2)
S_RANGE = (1.0, 3.0)        # the scalar raises c_b, so s > 1; 3 is far beyond any fit


def unique_root(f, lo, hi, what, n_scan=9):
    """The single root of f on [lo, hi] by find_roots; raises otherwise."""
    r = find_roots(f, lo, hi, n_scan=n_scan)
    if len(r["roots"]) != 1:
        raise RuntimeError(f"{what}: {len(r['roots'])} roots on [{lo}, {hi}] "
                           f"({r['message']})")
    return float(r["roots"][0])


def k3_residual(s, P, implied):
    """Mean log ratio of implied to modelled unpromoted bulk CO2 at scalar s."""
    cb = np.array([bulk(P, a, None, k3_scale=s)["CO2"] for a, _ in CONDITIONS])
    return float(np.mean(np.log(implied / cb)))


def unpromoted_devs(P, s, t2=t2):
    """Mean deviation of E_u and of J_u from Table 2 at scalar s."""
    films = [run_case(P, a, p, None, s) for a, p in CONDITIONS]
    e = np.mean([f.enhancement() / v - 1 for f, v in zip(films, t2.E_u)])
    j = np.mean([f.J_co2() / v - 1 for f, v in zip(films, t2.J_u)])
    return float(e), float(j)


def study_fit(P=P, t2=t2):
    """The scalar from the c_b identity, its per-condition spread, and the E_u refit."""
    implied = implied_cb(t2, P)["u"]
    s = unique_root(lambda x: k3_residual(x, P, implied), *S_RANGE, "k3 scale")
    per = implied / np.array([bulk(P, a)["CO2"] for a, _ in CONDITIONS])
    s_eu = unique_root(lambda x: unpromoted_devs(P, x, t2)[0], *S_RANGE,
                       "E_u refit", n_scan=7)
    return {"k3_scale_fitted": s,
            "k3_scale_spread_per_condition": float(per.max() / per.min() - 1),
            "k3_scale_refitted_to_Eu": s_eu, "_per": per}


FIT = study_fit()
S_FIT = FIT["k3_scale_fitted"]
_res = k3_residual(S_FIT, P, CB_IMPL["u"])
print(f"fitted scale s = {S_FIT:.4f} (root of the mean log residual, which is "
      f"{_res:.1e} there)\n  K3_eff = {K3/S_FIT:.3f} mol/m3, equivalently "
      f"K_c2_eff = {P['KC2']*S_FIT:.2e} mol/m3 at fixed K_w;"
      f"\n  a shift of {np.log10(S_FIT):.3f} in pK_c2, or pK_w "
      f"{-np.log10(P['KW']):.2f} -> {-np.log10(P['KW']/S_FIT):.2f} on the K_w branch")
_per = FIT["_per"]
print("\nper-condition scalars the c_b identity demands:", np.round(_per, 4))
print(f"  span {_per.min():.3f}-{_per.max():.3f}, i.e. "
      f"{FIT['k3_scale_spread_per_condition']:.1%} spread over n = {_per.size} "
      "conditions. One scalar cannot satisfy\n  them exactly; the next cells show "
      "where the residual of the compromise lands.")
_cb1 = bulk(P, CONDITIONS[1][0])["CO2"]
print(f"\nat alpha = 0.4 the two printed rows force c_b = {CB_TWO_ROW:.2f} mol/m3; "
      f"the printed constants give {_cb1:.2f}")

sa, sb = S_FIT * 0.96, S_FIT * 1.03            # either side of the fit
(fa, ja), (fb, jb) = unpromoted_devs(P, sa), unpromoted_devs(P, sb)
_ls = np.log(sb / sa)
_pe, _pj = np.log((1 + fb) / (1 + fa)) / _ls, np.log((1 + jb) / (1 + ja)) / _ls
print(f"\nsensitivity of the two reproduction targets to s:")
print(f"  s = {sa:.3f}: mean E_u dev {fa:+.2%}, mean J_u dev {ja:+.2%}")
print(f"  s = {sb:.3f}: mean E_u dev {fb:+.2%}, mean J_u dev {jb:+.2%}")
print(f"  -> local exponents  E ~ s^{_pe:+.2f},  J ~ s^{_pj:+.2f}: E is "
      f"{abs(_pj/_pe):.1f}x less sensitive to the scalar than the flux is.")
S_FROM_EU = FIT["k3_scale_refitted_to_Eu"]
print(f"scalar refitted to E_u ALONE (root of the mean E_u deviation): "
      f"s = {S_FROM_EU:.4f} ({S_FROM_EU/S_FIT - 1:+.1%} from the c_b-fitted value)")

# corroboration on quantities NOT used in the fit: the promoted implied c_b
cb_check = pd.DataFrame([
    {"alpha": a, "p_kpa": p,
     "cb_DEA_recon": bulk(P, a, "DEA", S_FIT)["CO2"], "cb_DEA_implied": iD,
     "cb_HDA_recon": bulk(P, a, "HDA", S_FIT)["CO2"], "cb_HDA_implied": iH}
    for (a, p), iD, iH in zip(CONDITIONS, CB_IMPL["DEA"], CB_IMPL["HDA"])])
print("\npromoted bulk CO2, reconstructed vs implied (not fitted):")
print(cb_check.round(4).to_string(index=False))
_ci30, _cbD = c_interface(P, CONDITIONS[0][1]), float(cb_check.cb_DEA_implied.iloc[0])
print(f"\nThe implied promoted values inherit the 3-digit rounding of E and J; for "
      f"DEA at alpha={CONDITIONS[0][0]}\nc_i - J/(E kL) is {_cbD:.3f} out of "
      f"{_ci30:.3f}, so a 0.5% rounding in J or E moves it by "
      f"{0.005*(_ci30 - _cbD)/_cbD:.0%}.")

# %%
# ---- second pass: the reconstructed run -----------------------------------
def study_recon(P=P, s=None, t2=t2):
    """Table 2 with the fitted scalar on K_c2."""
    s = S_FIT if s is None else s
    df, dev, films = table2_comparison(P, s, t2=t2)
    a = lambda cols: np.abs(dev[cols].to_numpy())
    f_max = float(a(["F_DEA", "F_HDA"]).max())
    return {"Ju_maxdev_reconstructed": float(a(["J_u"]).max()),
            "Eu_maxdev_reconstructed": float(a(["E_u"]).max()),
            "E_promoted_maxdev_reconstructed": float(a(["E_DEA", "E_HDA"]).max()),
            "F_maxdev_reconstructed": f_max, "Fmax_branch_Kc2": f_max,
            "_df": df, "_dev": dev, "_films": films}


RC = study_recon()
df_rec, dev_rec = RC["_df"], RC["_dev"]
print("model (reconstructed K_c2; one scalar fitted to the paper's own c_b identity):")
print(df_rec.round(3).to_string(index=False))
print("\n(model - paper)/paper:")
print(dev_rec.round(3).to_string(index=False))
print(f"\nHEADLINE: unpromoted FLUX J_u max |dev| = {RC['Ju_maxdev_reconstructed']:.1%}")
print("  per condition, J_u: " + ", ".join(f"{d:+.1%}" for d in dev_rec["J_u"]))
print("  per condition, E_u: " + ", ".join(f"{d:+.1%}" for d in dev_rec["E_u"])
      + f"  (max {RC['Eu_maxdev_reconstructed']:.1%})")
print("  The E_u figures are one offset repeated, because the fitted c_b sits in")
print("  E's denominator and cancels most of the residual.")
print(f"\npromoted E: max |dev| = {RC['E_promoted_maxdev_reconstructed']:.1%}   "
      f"F: max |dev| = {RC['F_maxdev_reconstructed']:.1%}")
print(f"\nthe paper's own Table 2 gives F(DEA) down to {t2.F_DEA.min():.2f} and F(HDA) "
      f"to {t2.F_HDA.min():.2f}; its Table 3\nquotes the numerical ranges and the "
      "experimental '+-4' and '+-6' (read as ~4 and ~6).")
print(f"Here: F(DEA) {df_rec.F_DEA.min():.2f}-{df_rec.F_DEA.max():.2f}, "
      f"F(HDA) {df_rec.F_HDA.min():.2f}-{df_rec.F_HDA.max():.2f}: both amines "
      "promote, HDA the stronger, as the paper found.")

# %%
# ---- parity plot: all twelve enhancement factors --------------------------
fig, ax = plt.subplots(figsize=(5.6, 5.2))
paper_E = t2[E_COLS].to_numpy()
lims = np.array([0.5 * paper_E.min(), 1.6 * paper_E.max()])
ax.fill_between(lims, lims * 0.9, lims * 1.1, color="0.88", label="±10% band")
ax.plot(lims, lims, "-", color="0.4", lw=1)
for j, (key, label) in enumerate([("u", "unpromoted"), ("DEA", "DEA"), ("HDA", "HDA")]):
    ax.plot(paper_E[:, j], PR["_df"][E_COLS[j]], "o", mfc="none", color=COL[key],
            ms=7, label=f"{label}, printed constants")
    ax.plot(paper_E[:, j], df_rec[E_COLS[j]], "o", color=COL[key], ms=7,
            label=f"{label}, reconstructed")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("enhancement factor E, paper Table 2")
ax.set_ylabel("enhancement factor E, this page")
ax.legend(fontsize=8, loc="upper left")
ax.set_title("Twelve enhancement factors, before and after the\n"
             "one-scalar $K_{c2}$ reconstruction")
plt.tight_layout(); plt.show()

# %%
# ---- where does the scalar sit? three candidates, one data set ------------
# c_b at fixed loading scales as 1/(K_c1 K_3) with K_3 = K_w/K_c2, so raising
# c_b by s can be done by dividing K_c1 by s, by multiplying K_c2 by s (the
# reconstruction), or by dividing K_w by s. The last also divides
# K_B = K_w/K_p and K_eq = K_c K_c1 K_w/K_p, which only the promoted cases see.
def summarise(dev):
    a = lambda cols: np.abs(dev[cols].to_numpy())
    out = dict(Ju_max=a(["J_u"]).max(), Eu_max=a(["E_u"]).max())
    if "E_DEA" in dev:
        e_p, f_p = a(["E_DEA", "E_HDA"]), a(["F_DEA", "F_HDA"])
        out.update(Ep_max=e_p.max(), Ep_rms=float(np.sqrt((e_p**2).mean())),
                   F_max=f_p.max(), F_rms=float(np.sqrt((f_p**2).mean())))
    return out


def study_branch_kw(P=P, s=None):
    """The same scalar on K_w: K_3, K_B and K_eq all move."""
    s = S_FIT if s is None else s
    _, dev, _ = table2_comparison(dict(P, KW=P["KW"] / s), 1.0)
    return {"Fmax_branch_Kw": float(summarise(dev)["F_max"]), "_dev": dev}


def study_branch_kc1(P=P, s=None, amines=(None, "DEA", "HDA")):
    """The same scalar on K_c1: the bulk CO2 moves, the bulk OH- does not."""
    s = S_FIT if s is None else s
    _, dev, _ = table2_comparison(dict(P, KC1=P["KC1"] / s), 1.0, amines=amines)
    return {"Eu_max_branch_Kc1": float(summarise(dev)["Eu_max"]), "_dev": dev}


BKW, BKC1 = study_branch_kw(), study_branch_kc1()
branch = {"K_c2 (as coded)": summarise(dev_rec), "K_w": summarise(BKW["_dev"]),
          "K_c1": summarise(BKC1["_dev"])}
print("deviation (model - paper)/paper, same scalar in three places:\n")
print(f"{'branch':18s} {'J_u max':>8s} {'E_u max':>8s} {'E_p max':>8s} "
      f"{'E_p rms':>8s} {'F max':>7s} {'F rms':>7s}")
for k, v in branch.items():
    print(f"{k:18s} {v['Ju_max']:8.1%} {v['Eu_max']:8.1%} {v['Ep_max']:8.1%} "
          f"{v['Ep_rms']:8.1%} {v['F_max']:7.1%} {v['F_rms']:7.1%}")
print("\nK_c1 fixes the bulk CO2 and nothing else: it leaves the bulk OH- at its")
print("printed value, so the hydration Hatta number, and E_u, stay high. It is")
print("rejected. K_c2 and K_w are the same change on every unpromoted quantity and")
print("nearly so on the promoted enhancement factors; only the promotion factors")
print("prefer the K_c2 branch, at the level of the residual scatter. The page uses")
print("K_c2 and does not claim the split is resolved.")
# scaling c_carb instead is algebraically the same change as K_c1: OH- =
# K_3 [CO3]/[HCO3] is invariant under a uniform scaling of the carbonate
# concentrations, so it fails identically and is not run separately.

# %% [markdown]
# ### An independent route to the same scalar, using no enhancement factor
#
# The fit above consumes the printed $(J, E)$ pairs. The paper's three
# concentration-profile figures offer a route that consumes neither: their
# *unpromoted* (dashed) curves are a direct picture of the bulk speciation, so
# reading the bulk dissolved CO2 off a figure and dividing by what the printed
# constants give returns the scalar with no Table 2 in the calculation at all.
#
# **Review status: a reconstruction cross-check, pending maintainer visual
# review.** These readings are not a dataset, no CSV is shipped, and nothing
# on this page depends on them. The extraction script and its numbered
# overlays are in `queue_cases/F3.5/review/`; the overlays are drawn on the
# copyrighted page images and are git-ignored. The reading identifies two
# curves that are easy to confuse: Figure 2's *promoted* CO3²⁻ is the solid
# rule that leaves the plot box through the top at the right-hand edge, and
# Figure 3's two carbonate curves converge to within one line width at the
# frame.
#
# **Two falsifiable tests.** (1) The dashed curves must be electroneutral
# against $[\mathrm{K^+}] = 2c_{carb}$; swapping the HCO3⁻ and CO3²⁻
# assignments breaks that by tens of percent. (2) The solid curves must satisfy
# $[\mathrm{K^+}] + [\mathrm{AmH^+}] = [\mathrm{HCO_3^-}] +
# 2[\mathrm{CO_3^{2-}}] + [\mathrm{AmCOO^-}]$, and the amine must total
# $c_{Am}$. Test (2) cannot be written down without the promoted CO3²⁻ curve.
# Both tests are printed below with their worst deviations
# (`fig_charge_dev_unpromoted_max`, `fig_charge_dev_promoted_max`).
#
# Test (1) is also what shows Figures 2 and 3 to be at half the loading their
# insets print, the internal inconsistency flagged under Parameters; the
# promoted carbonates do not enter it.

# %%
# ---- figure readings (see queue_cases/F3.5/review/extract_figs.py) --------
# 600 dpi renders; axes calibrated on the printed ticks; curve centres tracked
# through the dashes and extrapolated to x/delta = 1. Hand reads made by this
# page, pending maintainer visual review: a cross-check, not data.
FIG_READ = {   # figure -> bulk (x/delta = 1) readings, mol/m3;
               # _u = dashed (unpromoted), _p = solid (promoted)
    2: dict(alpha=0.2, inset=0.4, CO2_u=0.1172, HCO3_u=805.5, CO3_u=1596.4,
            HCO3_p=693.3, CO3_p=1634.7, AmCOO_p=77.2, AmH_p=47.2, Am_p=484.4,
            Am_p_interface=369.0),
    3: dict(alpha=0.2, inset=0.4, CO2_u=0.1160, HCO3_u=803.8, CO3_u=1598.2,
            HCO3_p=579.8, CO3_p=1591.9, AmCOO_p=228.3, AmH_p=0.0, Am_p=375.9),
    4: dict(alpha=0.6, inset=0.6, CO2_u=2.1810, HCO3_u=2390.9, CO3_u=802.3,
            HCO3_p=2086.5, CO3_p=966.3, AmCOO_p=133.7, AmH_p=145.0, Am_p=314.6),
}


def study_figures(P=P):
    """Charge balances of the readings, and the scalar from the dashed bulk CO2."""
    k_plus = 2 * P["C_CARB"]
    chg_u = max(abs((f["HCO3_u"] + 2 * f["CO3_u"]) / k_plus - 1) for f in FIG_READ.values())
    chg_p = max(abs((f["HCO3_p"] + 2 * f["CO3_p"] + f["AmCOO_p"]) / (k_plus + f["AmH_p"]) - 1)
                for f in FIG_READ.values())
    s_fig = np.array([f["CO2_u"] / bulk(P, f["alpha"])["CO2"] for f in FIG_READ.values()])
    return {"s_from_figures_mean": float(s_fig.mean()),
            "fig_charge_dev_unpromoted_max": float(chg_u),
            "fig_charge_dev_promoted_max": float(chg_p), "_s_fig": s_fig}


FG = study_figures()
C_CARB, C_AM = P["C_CARB"], P["C_AM"]
print("(a) the DASHED curves against the potassium: [HCO3-] + 2[CO3 2-] = [K+]")
print("figure   HCO3-_u   CO3(2-)_u   charge sum   alpha(HCO3)  alpha(CO3)  inset")
for k, f in FIG_READ.items():
    print(f"  {k}    {f['HCO3_u']:8.1f}  {f['CO3_u']:9.1f}   "
          f"{f['HCO3_u'] + 2*f['CO3_u']:9.1f}   "
          f"{f['HCO3_u']/(2*C_CARB):10.2f}  {1 - f['CO3_u']/C_CARB:10.2f}  "
          f"{f['inset']:6.1f}")
print(f"\n[K+] = 2 c_carb = {2*C_CARB:.0f} mol/m3: worst deviation "
      f"{FG['fig_charge_dev_unpromoted_max']:.2%}, so the axis assignment is right;")
print("Figures 2 and 3 plot the alpha = 0.2 composition under a label that says 0.4.")
_swap = [(f["CO3_u"] + 2*f["HCO3_u"])/(2*C_CARB) - 1 for f in FIG_READ.values()]
print("with HCO3- and CO3(2-) swapped, the same test gives "
      + ", ".join(f"{d:+.0%}" for d in _swap) + ", so the assignment is not a free choice.")

print("\n(b) the SOLID curves: [K+] + [AmH+] = [HCO3-] + 2[CO3 2-] + [AmCOO-]")
print("    (Figure 3 draws DEAH+ on the x axis, so [AmH+] = 0 there)")
print("figure   cations    anions     dev   CO3(2-)_p read   demanded by charge   amine sum")
for k, f in FIG_READ.items():
    _cat = 2 * C_CARB + f["AmH_p"]
    _ani = f["HCO3_p"] + 2 * f["CO3_p"] + f["AmCOO_p"]
    print(f"  {k}    {_cat:8.1f}  {_ani:8.1f}  {_ani/_cat-1:+6.2%}   {f['CO3_p']:12.1f}   "
          f"{0.5*(_cat - f['HCO3_p'] - f['AmCOO_p']):16.1f}   "
          f"{f['Am_p'] + f['AmH_p'] + f['AmCOO_p']:9.1f}")
print(f"\nworst deviation {FG['fig_charge_dev_promoted_max']:.2%}; the amine sum is "
      f"to be compared with c_Am = {C_AM:.0f} mol/m3.\n")

# the price of reading Table 2 on the Figure 2/3 convention instead:
# a row labelled alpha would then carry the composition bulk() gives at alpha/2
S_ALT = CB_IMPL["u"] / np.array([bulk(P, a / 2)["CO2"] for a, _ in CONDITIONS])
print(f"if Table 2's loadings are read on the Figure 2/3 convention, the c_b identity "
      f"demands scalars {np.round(S_ALT, 1)},\nspan {S_ALT.max()/S_ALT.min()-1:.0%}, "
      f"against {FIT['k3_scale_spread_per_condition']:.0%} on the Table 2/Figure 4 "
      "convention. That is why this page follows Table 2 and Figure 4.\n")

print("unpromoted bulk CO2, deviation convention (model - figure)/figure:")
print("  fig  alpha    figure   printed    dev   reconstructed    dev   -> s")
for (k, f), s_k in zip(FIG_READ.items(), FG["_s_fig"]):
    cb_p, cb_r = bulk(P, f["alpha"])["CO2"], bulk(P, f["alpha"], None, S_FIT)["CO2"]
    print(f"   {k}   {f['alpha']:.1f}  {f['CO2_u']:8.4f}  {cb_p:8.4f} "
          f"{cb_p/f['CO2_u']-1:+6.1%}  {cb_r:13.4f} {cb_r/f['CO2_u']-1:+6.1%}  {s_k:5.3f}")
print(f"\nfigure-only scalars span {FG['_s_fig'].min():.3f}-{FG['_s_fig'].max():.3f} "
      f"(mean {FG['s_from_figures_mean']:.3f}); the Table 2 fit gave {S_FIT:.3f}.")
print(f"Figures 2 and 3 plot the SAME state, so their "
      f"{abs(FG['_s_fig'][0]/FG['_s_fig'][1]-1):.1%} difference measures the reading "
      "error.\nNo enhancement factor enters this route.")

# %%
# ---- concentration profiles at the conditions Figure 2 actually plots -----
_fs = RC["_films"][(0.2, 30.0)]      # alpha = 0.2 (what Fig. 2 shows), its 30 kPa
fu, fh = _fs[None], _fs["HDA"]
SPECIES_COL = {"CO2": "tab:blue", "HCO3": "tab:orange", "CO3": "tab:green",
               "Am": "tab:purple", "AmH": "tab:brown", "AmCOO": "tab:red"}
LABELS = {"HCO3": "HCO$_3^-$", "CO3": "CO$_3^{2-}$", "Am": "HDA",
          "AmH": "HDAH$^+$", "AmCOO": "HDACOO$^-$"}
fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
ax = axes[0]
ax.plot(fh.x_c, fh.c[:, 0], color=SPECIES_COL["CO2"], label="CO2, HDA-promoted")
ax.plot(fu.x_c, fu.c[:, 0], "--", color=SPECIES_COL["CO2"], label="CO2, unpromoted")
ax.set_ylabel("CO2 (mol m$^{-3}$)")
ax.legend(fontsize=8)
ax.set_title(r"Absorption at $\alpha$ = 0.2, $p_{CO_2}$ = 30 kPa"
             "\n(the composition the paper's Fig. 2 plots)")
ax = axes[1]
LABEL_POS = {"HCO3": (0.55, 6), "CO3": (0.55, 6), "Am": (0.55, 6),
             "AmH": (0.78, 6), "AmCOO": (0.25, -14)}
for sp in ["HCO3", "CO3", "Am", "AmH", "AmCOO"]:
    prof = fh.c[:, fh.idx[sp]]
    ax.plot(fh.x_c, prof, color=SPECIES_COL[sp])
    x_lab, dy = LABEL_POS[sp]
    ax.annotate(LABELS[sp], (x_lab, prof[np.searchsorted(fh.x_c, x_lab)]),
                textcoords="offset points", xytext=(0, dy), fontsize=8,
                color=SPECIES_COL[sp])
for sp in ["HCO3", "CO3"]:
    ax.plot(fu.x_c, fu.c[:, fu.idx[sp]], "--", color=SPECIES_COL[sp], lw=1)
ax.set_xlabel(r"$x/\delta$"); ax.set_ylabel("ions, amine (mol m$^{-3}$)")
ax.set_xlim(0, 1)
plt.tight_layout(); plt.show()

# everything quoted below is computed from the plotted solution
GRAD_RATIO = fh.J_co2() / fu.J_co2()
AM_DEPLETION = 1 - fh.c[0, fh.idx["Am"]] / fh.cb_dict["Am"]
print("What the profiles show: the HDACOO- concentration falls from interface to")
print("bulk, so carbamate carries CO2 into the bulk (the shuttle mechanism),")
print(f"while free HDA is depleted at the interface by {AM_DEPLETION:.0%}. The promoted")
print(f"CO2 flux at the interface is {GRAD_RATIO:.1f}x the unpromoted one, and HCO3-")
print("production in the film is enhanced (solid vs dashed): the hydration route")
print("is accelerated too. Both mechanisms operate at once, the paper's central")
print("point against the either-or models.")
print(f"\nAgainst Figure 2 itself: free HDA runs {fh.c[0, fh.idx['Am']]:.0f} at the "
      f"interface to {fh.cb_dict['Am']:.0f} mol/m3 in the bulk;\nthe figure reads "
      f"{FIG_READ[2]['Am_p_interface']:.0f} -> {FIG_READ[2]['Am_p']:.0f} (same review "
      "caveat as above).")

# %% [markdown]
# ### The one out-of-sample test in the paper, and it prefers the printed constants
#
# Everything compared so far is in Table 2, which the reconstruction was fitted
# against. The paper states exactly two numbers that are **not** in Table 2,
# both about the desorption case of its Figure 4: the promotion factor, and that
# adding the amine "lowers the driving force for desorption by about one
# third". A reconstruction fitted to Table 2 cannot have absorbed either.
#
# **Run both ways, the printed constants win on both counts.** The numbers are
# in the cell below. This is evidence *against* the reconstruction, and it is
# the only out-of-sample evidence the paper offers. The page does not explain
# it away.
#
# What can be said on the other side, and no more: the reconstruction is
# supported at the *same* loading by Figure 4's own bulk CO2 (on the very state
# this desorption case starts from; the cell below prints by how much each set
# of constants misses it), and by the twelve Table 2 numbers, which the printed
# constants miss systematically. So the paper's own record is internally
# inconsistent here: its Figure 4 argues one way and its section 4.2 text the
# other. A reader weighing this page should weigh that, rather than take either
# half alone.

# %%
# ---- desorption: the promotion factor stated in section 4.2 (Fig. 4) ------
F_DESORB_PAPER = 3.81          # section 4.2 text, not in Table 2 (this page's transcription)
DESORB_ALPHA, DESORB_P = 0.6, 0.030      # Fig. 4: alpha = 0.6, p = 30 Pa (kPa here)


def study_desorption(P=P, s=None):
    """F and the driving-force cut, with the printed constants and with the scalar."""
    s = S_FIT if s is None else s
    out = {}
    ci = c_interface(P, DESORB_P)
    for suffix, sc in (("_printed_constants", 1.0), ("", s)):
        fu_d = run_case(P, DESORB_ALPHA, DESORB_P, None, sc)
        fh_d = run_case(P, DESORB_ALPHA, DESORB_P, "HDA", sc)
        F = fh_d.J_co2() / fu_d.J_co2()
        out[f"desorption_F_dev{suffix}"] = float(F / F_DESORB_PAPER - 1)
        out[f"desorption_drive_cut{suffix}"] = float(
            1 - (fh_d.cb_dict["CO2"] - ci) / (fu_d.cb_dict["CO2"] - ci))
        out[f"_F{suffix}"] = F
    return out


DS = study_desorption()
for tag, suf in (("printed constants", "_printed_constants"), ("reconstructed", "")):
    print(f"{tag:18s}: F = {DS['_F' + suf]:5.2f} (paper {F_DESORB_PAPER}, dev "
          f"{DS['desorption_F_dev' + suf]:+5.1%});   driving force cut by "
          f"{DS['desorption_drive_cut' + suf]:4.1%} (paper: 'about one third')")
print("\nBOTH out-of-sample numbers come out better with the constants exactly as")
print("printed. This is reported as what it is: evidence against the reconstruction.")
_cb_fig4 = FIG_READ[4]["CO2_u"]
print(f"Its Figure 4, at this same alpha, shows a bulk CO2 of {_cb_fig4:.2f} mol/m3, "
      f"which the printed\nconstants miss by {bulk(P, 0.6)['CO2']/_cb_fig4 - 1:+.0%} "
      f"and the reconstruction by {bulk(P, 0.6, None, S_FIT)['CO2']/_cb_fig4 - 1:+.0%}.")

# the approximate eq. (15) the paper rejects, F = sqrt(1 + Am^2 Keq / AmH):
print("\ntheir eq. (15) vs the full model (absorption cases, reconstructed run):")
eq15_ratio = []
for (alpha, p) in CONDITIONS:
    for am in ("DEA", "HDA"):
        cbp = bulk(P, alpha, am, S_FIT)
        f15 = np.sqrt(1 + cbp["Am"]**2 * KEQ[am] / max(cbp["AmH"], 1e-12))
        fnum = df_rec.loc[(df_rec.alpha == alpha) & (df_rec.p_kpa == p), f"F_{am}"].iloc[0]
        eq15_ratio.append(f15 / fnum)
        print(f"  alpha={alpha} p={p:>5.0f} {am}: eq.(15) F = {f15:7.1f}   "
              f"full model F = {fnum:5.2f}   ratio {f15/fnum:5.1f}x")
print(f"eq. (15) overshoots by {min(eq15_ratio):.1f}x to {max(eq15_ratio):.0f}x, as the "
      "paper states ('much too\nhigh values'): the promoter is not at equilibrium near "
      "the interface.")

# %% [markdown]
# ## Validation
#
# The reproduction above is model against model. These checks are independent
# of the paper. They are **not equally strong**, and the list says which is
# which: some test the discretisation, others are consistency checks that are
# true by construction and can only catch a coding slip. Agreement at round-off
# is reported as a bound with its basis, never as evidence of accuracy.
#
# **The checks that test the numerics:**
#
# 1. **Pseudo-first-order limit.** With OH⁻ buffered and reaction (1)
#    irreversible, the film has the closed form
#    $E = \mathrm{Ha}(\cosh\mathrm{Ha} - c_b/c_i)/(\sinh\mathrm{Ha}\,(1 -
#    c_b/c_i))$, the `F3.1` result extended to a nonzero bulk concentration.
#    Swept over Hatta numbers from a reaction layer the size of the film down
#    to a small fraction of it. This is the check that can detect an
#    under-resolved boundary layer, and grid doubling (item 7) backs it.
# 2. **Van Krevelen-Hoftijzer.** An external published approximation for the
#    second-order irreversible film, matched at the Hatta number of the DEA
#    case. Nothing in this implementation was tuned to it.
#
# **Consistency checks, real but weaker than they look:**
#
# 3. **Physical absorption limit.** With all reactions off, E must equal 1.
#    Exact for *any* conservative scheme on a linear profile, so this tests unit
#    and dimensional bookkeeping, not accuracy.
# 4. **Dimensionless cross-assembly.** The amine-only irreversible
#    configuration of this dimensional seven-species code against an
#    independently written dimensionless two-species solve
#    (`second_order_dimensionless`). The two share the grid,
#    `construct_grad`/`construct_div`/`NumJac`/`newton` and the
#    discretisation; what is written twice is the non-dimensionalisation, the δ
#    and D scaling, the ν = 2 stoichiometry and $E_\infty$. So the agreement
#    tests the *dimensional bookkeeping* and cannot detect a discretisation
#    error common to both.
# 5. **Reversible to irreversible collapse.** Scaling $K_{eq}$ up drives the
#    reversible amine solve onto its irreversible limit. Both use the same
#    `Film` class, so this is self-consistency: it checks that the reverse term
#    vanishes as $K_{eq}\to\infty$.
# 6. **Conservation and electroneutrality.** Carbon-flux uniformity and the
#    amine-flux closure test the **stoichiometry matrix**, and nothing else.
#    Every reaction conserves carbon, so $\nabla\!\cdot\!\sum_j \gamma_j J_j = 0$
#    identically for any set of $D_j$: the break table halves
#    $D(\mathrm{HCO_3^-})$ and both closures stay inside their bounds. What
#    they catch is a wrong entry in the stoichiometry matrix (two break rows).
#    **The charge flux** vanishes for the same kind of reason: every reaction
#    conserves charge and only neutral CO2 crosses the interface, so at steady
#    state the charge flux is zero on every face for any set of $D_j$; it too
#    tests the stoichiometry. **Electroneutrality** of the concentrations is
#    what needs the equal diffusivities: with every charged species carrying
#    the same $D_{ion}$ and electroneutral Dirichlet data at the bulk edge, the
#    charge density satisfies a linear problem with zero source and zero
#    boundary data. The halved $D(\mathrm{HCO_3^-})$ row breaks it and leaves
#    the three flux checks inside their bounds.
# 7. **Physical bound and discretisation.** Every E below the
#    instantaneous-reaction ceiling from the total transportable base capacity
#    (a count with its margin); grid doubling and a further decade of $k_{inst}$
#    must leave E unchanged well below the reproduction tolerance.
# 8. **Transport-model sensitivity.** The same network solved as a Higbie
#    *penetration* model (a transient problem with the same operators plus an
#    accumulation term, `Penetration`) shows how much of the comparison depends
#    on the film idealisation. Its quadrature bias is calibrated on the
#    physical-absorption case and divided out; the size of that bias is
#    printed. The time integration is backward Euler on a geometric time grid
#    with `newton` per step, first order, so the number of time levels `n_t` is
#    refined at the condition the maximum is attained, the observed orders are
#    printed, and the result is Richardson-extrapolated from the finest triple.

# %%
# (3) physical absorption, (1) pseudo-first-order, (2)+(4)+(5) second order
P_SWEEP = 30.0                               # kPa, the alpha = 0.2 condition
HA_SWEEP = np.geomspace(0.1, 300, 13)        # reaction layer from the film to 1/300 of it


def study_physical(P=P, length_factor=1.0):
    """All reactions off: E - 1 on a film of length length_factor * delta."""
    f = Film(P, bulk(P, 0.2), cgas(P_SWEEP, P), None, enable_r1=False, enable_r34=False,
             length=length_factor * P["DCO2"] / P["KL"])
    f.solve(k_ladder=(0.0,))
    return {"physical_E1": abs(f.enhancement() - 1.0)}


def study_pfo(P=P, n_x=400):
    """Pseudo-first-order film against its closed form over HA_SWEEP."""
    cb02, ci = bulk(P, 0.2), c_interface(P, P_SWEEP)
    err = []
    for ha in HA_SWEEP:
        k1 = (ha * P["KL"]) ** 2 / P["DCO2"]          # first-order constant, 1/s
        fp = Film(P, cb02, cgas(P_SWEEP, P), None, enable_r34=False, n_x=n_x, pfo_k1=k1)
        fp.solve(k_ladder=(1.0,), tol=1e-8, bounds=None)   # linear: no floor
        err.append(abs(fp.enhancement() / pfo_enhancement(ha, cb02["CO2"] / ci) - 1))
    return {"pfo_max_rel_err": float(np.max(err)), "_err": err}


def study_second_order(P=P, nu_am=-2.0, keq_top=1e6):
    """Irreversible amine-only film: against the dimensionless assembly, VKH, and
    the reversible film as K_eq is scaled up to keq_top."""
    amb, ci = P["C_AM"], c_interface(P, P_SWEEP)
    cb_syn = {sp: 1e-8 for sp in SP}
    cb_syn["Am"] = amb
    kw = dict(enable_r1=False, enable_r34=False, nu_set={("Am", M.R2): nu_am})
    f_dim = Film(P, cb_syn, cgas(P_SWEEP, P), "DEA", irreversible2=True, **kw)
    f_dim.solve(k_ladder=(1.0,))
    e_dim = f_dim.enhancement()
    ha = np.sqrt(P["K2"]["DEA"] * amb * P["DCO2"]) / P["KL"]
    e_inf = 1.0 + P["DAM"] * amb / (2.0 * P["DCO2"] * ci)       # two amine per CO2
    e_dl, e_vkh = second_order_dimensionless(ha, e_inf - 1.0), vkh(ha, e_inf)
    collapse = []
    for ks in np.geomspace(1.0, keq_top, 4):
        f_rev = Film(P, cb_syn, cgas(P_SWEEP, P), "DEA", keq_scale=ks, **kw)
        f_rev.solve(k_ladder=(1.0,))
        collapse.append(f_rev.enhancement())
    return {"second_order_cross": abs(e_dim / e_dl - 1), "vkh_dev": abs(e_dim / e_vkh - 1),
            "collapse_residual": abs(collapse[-1] / e_dim - 1),
            "_e": (ha, e_inf, e_dim, e_dl, e_vkh), "_collapse": collapse}


PHY, PFO, SO = study_physical(), study_pfo(), study_second_order()
print(f"(3) physical absorption: |E - 1| at most 1e-9 (here {PHY['physical_E1']:.1e}) "
      "(exact for any conservative scheme: bookkeeping check)")
print(f"(1) pseudo-first-order vs closed form: max rel err {PFO['pfo_max_rel_err']:.2e}"
      f"\n    over Ha in [{HA_SWEEP[0]:g}, {HA_SWEEP[-1]:g}]; this one tests the "
      "boundary-layer resolution.")
ha, e_inf, e_dim, e_dl, e_vkh = SO["_e"]
print(f"(2) Ha = {ha:.1f}, E_inf = {e_inf:.0f}: this code {e_dim:.3f} vs the published"
      f"\n    Van Krevelen-Hoftijzer approximation {e_vkh:.3f}, dev {SO['vkh_dev']:.2e}")
print(f"(4) same case against the dimensionless assembly: {e_dl:.3f}, rel dev "
      f"at most 1e-6 (here {SO['second_order_cross']:.1e})\n    (shared grid and operators: tests the "
      "dimensional bookkeeping only)")
print("(5) E vs Keq scale " + ", ".join(f"{e:.2f}" for e in SO["_collapse"])
      + f" -> irreversible {e_dim:.2f};\n    residual {SO['collapse_residual']:.1e} "
      "(same Film class both sides)")

# %%
# (6) conservation + electroneutrality, (7) grid and k_inst, on HDA at alpha = 0.2
def study_film_checks(P=P, s=None, n_x=400, k_top=M.K_LADDER[-1], nu_set=None):
    """Closures on the stiffest case, and its change on grid doubling and on one
    more decade of k_inst."""
    s = S_FIT if s is None else s
    ladder = tuple(k for k in M.K_LADDER if k <= k_top)
    kw = dict(nu_set=nu_set)
    base = run_case(P, 0.2, P_SWEEP, "HDA", s, n_x=n_x, solve_kw=dict(k_ladder=ladder), **kw)
    fine = run_case(P, 0.2, P_SWEEP, "HDA", s, n_x=2 * n_x, solve_kw=dict(k_ladder=ladder), **kw)
    kup = run_case(P, 0.2, P_SWEEP, "HDA", s, n_x=n_x,
                   solve_kw=dict(k_ladder=ladder + (10 * ladder[-1],)), **kw)
    cl = closures(base)
    e0 = base.enhancement()
    return {"carbon_flux_closure": cl["carbon"], "amine_flux_closure": cl["amine"],
            "charge_flux": cl["charge"], "electroneutrality": cl["electroneutrality"],
            "grid_delta": abs(fine.enhancement() / e0 - 1),
            "kinst_delta": abs(kup.enhancement() / e0 - 1), "_E": e0}


def study_gas(P=P, films=None):
    """Gas-side share of the total resistance at the largest E of the run."""
    films = RC["_films"] if films is None else films
    e_max = max(f.enhancement() for fs in films.values() for f in fs.values())
    return {"gas_side_fraction": float(P["M_DIST"] * e_max * P["KL"] / P["KG"]),
            "_E_max": e_max}


def study_ceiling(films=None, with_carbonate=True):
    """How many of the twelve cases lie below their instantaneous ceiling."""
    films = RC["_films"] if films is None else films
    r = np.array([f.enhancement() / ceiling(f, with_carbonate)
                  for fs in films.values() for f in fs.values()])
    return {"cases_below_ceiling": int(np.sum(r < 1.0)), "_ratio": r}


FC, GS, CE = study_film_checks(), study_gas(), study_ceiling()
print(f"(6) carbon-flux uniformity (here {FC['carbon_flux_closure']:.1e}), amine flux/J "
      f"(here {FC['amine_flux_closure']:.1e}): tests of the stoichiometry matrix")
print(f"    charge flux/J (here {FC['charge_flux']:.1e}), stoichiometry again; "
      f"electroneutrality (here {FC['electroneutrality']:.1e}), the equal-ion-diffusivity premise")
print(f"(7) {CE['cases_below_ceiling']} of {CE['_ratio'].size} cases below their "
      f"instantaneous ceiling (largest E/ceiling {CE['_ratio'].max():.3f})")
print(f"    E change on grid doubling 400 -> 800: {FC['grid_delta']:.1e}; on k_inst "
      f"{M.K_LADDER[-1]:.0e} -> {10*M.K_LADDER[-1]:.0e}: {FC['kinst_delta']:.1e}")
print(f"    gas-side share of the total resistance at the largest E of the run "
      f"(E = {GS['_E_max']:.0f}):\n    m E k_L / k_g = {GS['gas_side_fraction']:.1e}, "
      "so the Dirichlet interface stands.")

# %%
# (8) transport-model sensitivity: the network as a Higbie penetration model
N_T = 640                         # time levels of the table below
PEN_LADDER = (160, 320, 640, 1280)
PEN_TOL = 1e-9                    # Newton step per time level, of the largest concentration
print(f"contact time theta = 4 D/(pi kL^2) = {M.contact_time(P)*1e3:.0f} ms")


def pen_bias(P, n_t, **kw):
    """Physical-absorption mean enhancement of the penetration solve (1 when exact)."""
    p0 = Penetration(P, bulk(P, 0.2), cgas(P_SWEEP, P), None, enable_r1=False,
                     enable_r34=False, theta_factor=kw.get("theta_factor", 1.0))
    p0.solve_transient(n_t=n_t, tol_rel=kw.get("tol_rel", PEN_TOL),
                       uniform=kw.get("uniform", False))
    return p0.enhancement_avg()


def pen_E(P, alpha, p, amine, s, n_t, bias, **kw):
    """Penetration enhancement with the quadrature bias divided out."""
    pen = Penetration(P, bulk(P, alpha, amine, s), cgas(p, P), amine, k3_scale=s,
                      theta_factor=kw.get("theta_factor", 1.0))
    pen.solve_transient(n_t=n_t, tol_rel=kw.get("tol_rel", PEN_TOL),
                        uniform=kw.get("uniform", False))
    return pen.enhancement_avg() / bias


def study_pen(P=P, s=None, n_t=N_T, films=None):
    """Penetration against film at every case; the largest relative difference."""
    s = S_FIT if s is None else s
    films = RC["_films"] if films is None else films
    b = pen_bias(P, n_t)
    rows = []
    for (alpha, p), fs in films.items():
        row = {"alpha": alpha, "p_kpa": p}
        for am, f in fs.items():
            row[f"E_{am or 'u'}_pen"] = pen_E(P, alpha, p, am, s, n_t, b)
            row[f"E_{am or 'u'}_film"] = f.enhancement()
        rows.append(row)
    df = pd.DataFrame(rows)
    tags = ("u", "DEA", "HDA")
    ratio = np.array([[r[f"E_{t}_pen"] / r[f"E_{t}_film"] for t in tags]
                      for _, r in df.iterrows()])
    # the max over the twelve printed conditions: a discrete set, so an argmax
    i_row, i_col = np.unravel_index(np.abs(ratio - 1).argmax(), ratio.shape)
    am = (None, "DEA", "HDA")[i_col]
    worst = (df.alpha.iloc[i_row], df.p_kpa.iloc[i_row], am,
             films[(df.alpha.iloc[i_row], df.p_kpa.iloc[i_row])][am].enhancement())
    return {"pen_vs_film_max": float(np.max(np.abs(ratio - 1))), "_df": df,
            "_worst": worst, "_bias": b}


PEN = study_pen()
WORST = PEN["_worst"]
print(f"physical-absorption quadrature bias {PEN['_bias'] - 1:+.2%} at n_t = {N_T} "
      "(divided out below)")
print(PEN["_df"].round(1).to_string(index=False))
print(f"penetration vs film: E differs by at most {PEN['pen_vs_film_max']:.1%}, at "
      f"alpha = {WORST[0]}, p = {WORST[1]:.0f} kPa, {WORST[2] or 'unpromoted'}")

# %%
# (8b) time refinement of the penetration result at the condition of its maximum
def study_pen_ladder(P=P, s=None, worst=None, ns=PEN_LADDER, known=None, **kw):
    """Observed orders of the film-penetration gap in n_t, and its Richardson limit
    from the finest triple, at one condition."""
    s = S_FIT if s is None else s
    alpha, p, am, e_film = WORST if worst is None else worst
    known = known or {}

    def gap(n_t):
        if n_t in known:
            return known[n_t]
        return abs(pen_E(P, alpha, p, am, s, n_t, pen_bias(P, n_t, **kw), **kw) / e_film - 1)

    st = observed_orders(gap, ns)
    vals = np.asarray(st["values"], float)
    orders = np.asarray(st["orders"], float)
    return {"pen_vs_film_max_richardson": float(st["extrapolated"]),
            "pen_richardson_error_estimate": float(st["error_estimate"]),
            "pen_time_order": float(orders[np.argmax(np.abs(orders - 1.0))]),
            "pen_ladder_over_newton_tol": float(abs(vals[-1] - vals[-2])
                                                / kw.get("tol_rel", PEN_TOL)),
            "_st": st}


LAD = study_pen_ladder(known={N_T: PEN["pen_vs_film_max"]})
_st = LAD["_st"]
print(f"(8b) n_t ladder at alpha = {WORST[0]}, p = {WORST[1]:.0f} kPa, "
      f"{WORST[2] or 'unpromoted'}:")
for n_t, v in zip(PEN_LADDER, _st["values"]):
    print(f"     n_t = {n_t:5d}   E_pen/E_film - 1 = {v:.6f}")
print(f"     observed orders {', '.join(f'{o:.3f}' for o in _st['orders'])} "
      "(backward Euler: 1)")
print(f"     Richardson limit from the finest triple {LAD['pen_vs_film_max_richardson']:.4f}"
      f", error of the finest level {LAD['pen_richardson_error_estimate']:.1e};")
print(f"     the n_t = {N_T} table value is "
      f"{abs(PEN['pen_vs_film_max']/LAD['pen_vs_film_max_richardson'] - 1):.1%} below it.")
print(f"     finest successive difference over the Newton tolerance "
      f"{LAD['pen_ladder_over_newton_tol']:.1e}: the orders measure the time step.")
print("     The penetration model gives the LARGER E, so the reproduction does not")
print("     hinge on the film idealisation, and the film variant, which the paper's")
print("     section 4.2 names, is the one that matches Table 2.")

# %% [markdown]
# ### Claims
#
# Every number the page states, with how it is checked. Reproduction results
# (deviations from Table 2, the fitted scalar and its refits, the desorption
# numbers, the penetration gap and its extrapolation) are values held to 5 % of
# the committed baseline. Transcription and figure-reading checks, the
# closed-form and cross-assembly gaps, closures, grid and $k_{inst}$ effects,
# the gas-side share, the time order and its round-off companion are bounds
# with their basis. The ceiling statement is a count with its margin.

# %%
BASE = {study_transcription: TR, study_fit: FIT, study_printed: PR, study_recon: RC,
        study_branch_kw: BKW, study_branch_kc1: BKC1, study_figures: FG,
        study_desorption: DS, study_physical: PHY, study_pfo: PFO,
        study_second_order: SO, study_film_checks: FC, study_gas: GS,
        study_ceiling: CE, study_pen: PEN, study_pen_ladder: LAD}
ALL = {k: v for s in BASE.values() for k, v in s.items() if not k.startswith("_")}
BOUNDS = {
    "table2_F_vs_Jratio_maxdev": dict(max=5e-3, basis="chosen",
        note="printed F against the ratio of printed fluxes, each printed to three "
             "significant digits: rounding alone reaches a few 1e-3"),
    "fig_charge_dev_unpromoted_max": dict(max=5e-3, basis="chosen",
        note="charge balance of the dashed bulk readings against [K+] = 2 c_carb; the "
             "reading resolution is a few mol/m3 on thousands"),
    "fig_charge_dev_promoted_max": dict(max=5e-3, basis="chosen",
        note="charge balance of the solid bulk readings, the same reading resolution"),
    "pfo_max_rel_err": dict(max=1e-4, basis="chosen",
        note="graded n_x = 400 film against the closed form over Ha 0.1 to 300; set by "
             "the second-order discretisation of the thinnest reaction layer"),
    "vkh_dev": dict(max=1e-3, basis="chosen",
        note="the film against the Van Krevelen-Hoftijzer approximation at the DEA "
             "Hatta number; the approximation error there is of order 1e-4"),
    "grid_delta": dict(max=1e-4, basis="chosen",
        note="E of HDA at alpha = 0.2 on 800 against 400 cells, held two decades "
             "below the reproduction deviations it supports"),
    "kinst_delta": dict(max=1e-4, basis="chosen",
        note="E of HDA at alpha = 0.2 with one more decade of k_inst, held two "
             "decades below the reproduction deviations"),
    "pen_time_order": dict(min=0.9, max=1.1, basis="theory",
        note="backward Euler; the order farthest from 1 over the two triples of the "
             "n_t ladder 160 to 1280"),
    "pen_ladder_over_newton_tol": dict(min=1e3, basis="theory",
        note="finest successive difference of the ladder over the Newton step "
             "tolerance per level: the orders measure the time step, not the solver"),
    "physical_E1": dict(max=1e-9, basis="identity",
        note="all reactions off: E = 1 for any conservative scheme on a linear profile"),
    "second_order_cross": dict(max=1e-6, basis="tolerance",
        note="two assemblies of the same problem on the same grid; their gap is the "
             "Newton tolerances, 1e-8 mol/m3 and 1e-11"),
    "collapse_residual": dict(max=1e-4, basis="theory",
        note="the reverse term of reaction (2) at K_eq scaled by 1e6 is of order "
             "1e-5 of the forward rate"),
    "carbon_flux_closure": dict(max=1e-8, basis="tolerance",
        note="every reaction conserves carbon, so the carbon flux is uniform up to the "
             "Newton residual"),
    "amine_flux_closure": dict(max=1e-8, basis="tolerance",
        note="every reaction conserves the amine count and the amine is non-volatile"),
    "charge_flux": dict(max=1e-8, basis="tolerance",
        note="charge-conserving reactions and only neutral CO2 through the interface: "
             "zero on every face up to the Newton residual, for any diffusivities"),
    "electroneutrality": dict(max=1e-9, basis="tolerance",
        note="equal ion diffusivities and electroneutral Dirichlet data: zero up to "
             "the Newton residual"),
    "gas_side_fraction": dict(max=1e-4, basis="chosen",
        note="m E k_L/k_g at the largest E of the reconstructed run: the share of the "
             "resistance the Dirichlet interface drops"),
}
for k, v in ALL.items():
    if k in BOUNDS:
        claims.bound(k, v, **BOUNDS[k])
    elif k == "cases_below_ceiling":
        claims.count(k, CE["_ratio"], 1.0, min_margin=1e-2, above=False,
                     note="E over its instantaneous ceiling, the twelve reconstructed cases")
    else:
        claims.value(k, v)
print(f"{len(claims.metrics)} claims registered")

# %% [markdown]
# ### The break table
#
# Every claim above needs a row that fails its own check. Each row injects one
# defect through a call argument or a physical parameter (a constant of Table 1,
# a printed value of Table 2, the reconstruction scalar, a diffusivity, an
# entry of the stoichiometry matrix, the grid, the $k_{inst}$ ladder, the film
# thickness, the time resolution or contact time of the penetration model) and
# reruns every study it names. A study that cannot be completed under a defect
# records `nan` for its claims, which fails them.
#
# Two rows are worth reading first. **"D(HCO3-) halved"** moves
# electroneutrality but leaves the carbon, amine and charge flux checks inside
# their bounds: those three test the stoichiometry, which the two stoichiometry
# rows show. **"penetration at n_t = 160 instead of 640"** moves the penetration gap by more than its
# tolerance: one time grid is not a converged result.

# %%
P_KW = dict(P, KW=0.9 * P["KW"])
P_KC1 = dict(P, KC1=1.1 * P["KC1"])
P_CC = dict(P, C_CARB=1.05 * P["C_CARB"])
P_D = dict(P, D=dict(P["D"], HCO3=0.5 * P["D"]["HCO3"]))
P_KG = dict(P, KG=1.0)
P_CAM = dict(P, C_AM=1.2 * P["C_AM"])
T2_JU = t2.copy()
T2_JU.loc[T2_JU.index[0], "J_u"] *= 1.1
U_ONLY = dict(amines=(None,))

ROWS = [
    ("K_w transcribed 10 % low (Table 1)",
     [(study_printed, dict(P=P_KW)), (study_desorption, dict(P=P_KW)),
      (study_fit, dict(P=P_KW)), (study_figures, dict(P=P_KW))]),
    ("K_c1 transcribed 10 % high (Table 1)",
     [(study_fit, dict(P=P_KC1)), (study_figures, dict(P=P_KC1)),
      (study_branch_kc1, dict(P=P_KC1, **U_ONLY))]),
    ("the reconstruction scalar left at s = 1",
     [(study_recon, dict(s=1.0)), (study_branch_kw, dict(s=1.0)),
      (study_branch_kc1, dict(s=1.0, **U_ONLY)), (study_desorption, dict(s=1.0))]),
    ("the first printed J_u read 10 % high (Table 2)",
     [(study_transcription, dict(t2=T2_JU)), (study_fit, dict(t2=T2_JU)),
      (study_recon, dict(t2=T2_JU))]),
    ("c_carb read 5 % high (Table 1)", [(study_figures, dict(P=P_CC))]),
    ("D(HCO3-) halved", [(study_film_checks, dict(P=P_D))]),
    ("reaction (3) does not make HCO3-",
     [(study_film_checks, dict(nu_set={("HCO3", M.R3): 0.0}))]),
    ("reaction (2) does not make AmCOO-",
     [(study_film_checks, dict(nu_set={("AmCOO", M.R2): 0.0}))]),
    ("the film on 50 cells instead of 400",
     [(study_pfo, dict(n_x=50)), (study_film_checks, dict(n_x=50))]),
    ("one amine per CO2 in reaction (2) instead of two",
     [(study_second_order, dict(nu_am=-1.0))]),
    ("the K_eq scan stopped at 1e2 instead of 1e6",
     [(study_second_order, dict(keq_top=1e2))]),
    ("the film 10 % thicker than D_CO2/k_L", [(study_physical, dict(length_factor=1.1))]),
    ("the k_inst ladder stopped at 1e4", [(study_film_checks, dict(k_top=1e4))]),
    ("k_g read as 1 m/s (Table 1)", [(study_gas, dict(P=P_KG))]),
    ("the instantaneous ceiling without the carbonate's capacity",
     [(study_ceiling, dict(with_carbonate=False))]),
    ("penetration at n_t = 160 instead of 640", [(study_pen, dict(n_t=160))]),
    ("Newton per time level stopped at 1e-3 of the largest concentration",
     [(study_pen_ladder, dict(tol_rel=1e-3))]),
    ("contact time D/k_L^2, the 4/pi of Higbie dropped",
     [(study_pen_ladder, dict(theta_factor=np.pi / 4))]),
    ("uniform time steps instead of the geometric time grid",
     [(study_pen_ladder, dict(uniform=True))]),
    ("c_Am read 20 % high (Table 1)", [(study_desorption, dict(P=P_CAM))]),
]
claims.define_rows([r[0] for r in ROWS])
BREAK = []
for row, studies in ROWS:
    for fn, kw in studies:
        keys = [k for k in BASE[fn] if not k.startswith("_")]
        try:
            out = fn(**kw)
        except Exception as exc:                 # noqa: BLE001 - a failed study is a result
            print(f"  {row}: {fn.__name__} raised {type(exc).__name__}:"
                  f" {str(exc)[:70]}; its claims record nan")
            out = {k: np.nan for k in keys}
        for k in keys:
            BREAK.append((row, k, out[k], claims.breaks(row, k, out[k])))

bt = pd.DataFrame(BREAK, columns=["row", "claim", "broken", "caught"])
caught = bt[bt.caught].groupby("row", sort=False).claim.apply(list)
print(f"{len(ROWS)} rows defined, {bt.row.nunique()} in the table; claims each row fails:")
for row, names in caught.items():
    print(f"  {row}:")
    print("      " + "\n      ".join(textwrap.wrap(", ".join(names), 96)))
idle = [r for r in bt.row.unique() if r not in caught.index]
print(f"\nrows that fail no claim: {idle if idle else 'none'}")
_kept = bt[(bt.row == "D(HCO3-) halved") & bt.claim.isin(
    ["carbon_flux_closure", "amine_flux_closure", "charge_flux"])]
print("D(HCO3-) halved, the three flux checks: "
      + ", ".join(f"{r.claim} (here {r.broken:.1e})" for r in _kept.itertuples())
      + "; all inside their bounds" * bool((~_kept.caught).all()))
missing = [n for n in claims.metrics if not claims.metrics[n].get("structural")
           and not bt[(bt.claim == n) & bt.caught].shape[0]]
print(f"claims no row fails: {missing if missing else 'none'}")

# %% [markdown]
# ## What pymrm adds
#
# This page reproduces a 1989 result; the model itself is unchanged, and no
# claim of improvement is made. What the reimplementation adds:
#
# - **A forensic reading of the paper's parameter set.** The printed
#   equilibrium constants are inconsistent with the paper's own Table 2 and its
#   three concentration-profile figures, and the structure of the model
#   localises the inconsistency to one constant, $K_3 = K_w/K_{c2}$. With a
#   single scalar on it, fitted to the paper's bulk-CO2 identity, the
#   unpromoted fluxes and the full Table 2 land within the deviations printed
#   above. It is a **reconstruction, not a closure**: the $K_{c2}$-versus-$K_w$
#   split is undetermined by the data, *why* the constant differs is not
#   established, and the paper's two out-of-sample numbers (the desorption
#   promotion factor and its driving-force statement) come out better with the
#   constants exactly as printed, reported above rather than set aside.
# - **The machinery, reusable.** In 1989 this calculation was a dedicated
#   research code. Here it is `model.py` on pymrm operators: `construct_grad`
#   with per-species boundary conditions on a graded grid, `construct_div`,
#   `NumJac` on the pointwise reaction network and `newton` with bounds,
#   continued in $k_{inst}$. The same `Film` class runs absorption and
#   desorption, any amine, any loading, reversible or not, and reduces to the
#   classical limits (`F3.1`'s Hatta and Van Krevelen-Hoftijzer results) on
#   demand; the flux is read from the face-flux operator.
# - **Model-form sensitivity the paper could not cheaply report.** The same
#   network solved as a Higbie penetration model, the same operators with an
#   accumulation term and `newton` per backward-Euler step, gives enhancement
#   factors larger than the film model's, by the extrapolated gap printed in
#   (8b), with `pymrm.checks.observed_orders` supplying the time orders and the
#   Richardson limit.
#
# What this page does **not** do: validate against experiment. The only
# measured numbers in the paper are "about 4" and "about 6", and both the
# paper's numerics and this reproduction bracket them.

# %% [markdown]
# ## Reuse
#
# - **Different amine:** add a row per constant to the Table 1 CSV ($K_p$,
#   $K_c$, $k_2$); nothing else changes. Tertiary amines: set $K_c = 0$ and
#   reaction (2) off; the model then predicts no promotion, the paper's
#   explanation of the Laddha-Danckwerts observation.
# - **Different solvent chemistry:** the `Film` class is generic in its
#   stoichiometry matrix and rate function; swap the network, keep the
#   operators, boundary conditions and the $k_{inst}$ continuation for any
#   instantaneous acid-base equilibria.
# - **The transient variant:** `Penetration.solve_transient` is backward Euler
#   on a geometric time grid, first order in `n_t`; refine `n_t` at the
#   condition the result is reported at, as (8b) does.
# - **Related pages:** `F3.1` (Hatta regimes, the enhancement-factor limits this
#   page reduces to), `F3.2` (the second-order film solved over the whole
#   chart), `A4.9` (multicomponent diffusion), `J4.8` (stoichiometry
#   conservation checks as validation).
#
# **Numerical points worth carrying away:** assemble the reaction source with
# the sign matched to `div(-grad)` (a wrong sign makes the coupled system
# anti-diffusive and Newton diverges); impose a numerically negligible
# gas-side resistance as a Dirichlet condition rather than a stiff Robin row
# (the ratio $k_g/(m\,k_L)$ is large and poisons the Jacobian conditioning);
# and continue the instantaneous-equilibrium rate constant up from a modest
# value rather than starting stiff.

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
# The README's literals, recorded so that a cut in prose precision shows in CI's
# loosening check (agreement.json "quoted").
SAYS = {
    'E_maxdev_printed': ['21.4%'],
    'k3_scale_fitted': ['1.4585'],
    'Ju_maxdev_reconstructed': ['3.8%'],
    'Eu_maxdev_reconstructed': ['0.6%'],
    'E_promoted_maxdev_reconstructed': ['7.6%'],
    'F_maxdev_reconstructed': ['5.6%'],
    'k3_scale_spread_per_condition': ['5.0%'],
    'k3_scale_refitted_to_Eu': ['1.4735'],
    'Eu_max_branch_Kc1': ['19.8%'],
    'Fmax_branch_Kw': ['8.6%'],
    'Fmax_branch_Kc2': ['5.6%'],
    's_from_figures_mean': ['1.445'],
    'fig_charge_dev_unpromoted_max': ['at most 5e-3'],
    'fig_charge_dev_promoted_max': ['at most 5e-3'],
    'desorption_F_dev_printed_constants': ['-2.3%'],
    'desorption_F_dev': ['-12.0%'],
    'desorption_drive_cut_printed_constants': ['34.7%'],
    'desorption_drive_cut': ['39.1%'],
    'pfo_max_rel_err': ['at most 1e-4'],
    'grid_delta': ['at most 1e-4'],
    'kinst_delta': ['at most 1e-4'],
    'vkh_dev': ['at most 1e-3'],
    'pen_vs_film_max': ['11.3%'],
    'pen_time_order': ['between 0.9 and 1.1'],
    'pen_ladder_over_newton_tol': ['at least 1e3'],
    'pen_vs_film_max_richardson': ['11.57%'],
    'pen_richardson_error_estimate': ['1.3e-3'],
    'physical_E1': ['at most 1e-9'],
    'second_order_cross': ['at most 1e-6'],
    'collapse_residual': ['at most 1e-4'],
    'carbon_flux_closure': ['at most 1e-8'],
    'amine_flux_closure': ['at most 1e-8'],
    'charge_flux': ['at most 1e-8'],
    'electroneutrality': ['at most 1e-9'],
    'cases_below_ceiling': ['12'],
    'gas_side_fraction': ['at most 1e-4'],
    'table2_F_vs_Jratio_maxdev': ['at most 5e-3'],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)
claims.report();   # the semicolon keeps the notebook from echoing the payload
