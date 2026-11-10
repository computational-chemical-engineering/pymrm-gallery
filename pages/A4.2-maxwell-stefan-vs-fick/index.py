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
# title: "Maxwell–Stefan vs Fick for multicomponent mixtures"
# description: "One physics, two bookkeepings: the Fick matrix [D] and the Maxwell–Stefan diffusivities describe the same mixture, and the page shows where the scalar effective diffusivity everyone actually uses falls apart."
# categories: [sec:A, struct:S9, tier:T0, data:tier4, phase:gas]
# date: 2026-07-31
# ---
#
# # Maxwell–Stefan vs Fick for multicomponent mixtures
#
# **Catalog ID:** `A4.2` · **Structure:** `S9` (implicit multicomponent flux) ·
# **Tier:** T0
#
# For a binary mixture there is one diffusion coefficient and everyone agrees
# what it means. For three or more species there are three descriptions in
# circulation, and they are not equally good:
#
# 1. the **Maxwell–Stefan** formulation: one drag coefficient
#    $\mathrm{D}_{ij}$ per pair of species;
# 2. the **generalized Fick matrix** $[D]$: $(n{-}1)^2$ coefficients with
#    composition-dependent off-diagonal terms;
# 3. the **scalar effective diffusivity** $D_{i,\mathrm{eff}}$: one number per
#    species, which is what most reactor models quietly assume.
#
# Krishna and Wesselingh's 1997 review is the standard argument that (1) and (2)
# are *exactly the same physics* in different bookkeeping, while (3) is a
# different, weaker claim that fails qualitatively as soon as the species
# gradients couple. This page reproduces the review's two worked matrix examples
# element by element, demonstrates the (1) ≡ (2) equivalence to machine precision
# in a pymrm transport solve, and then shows against the Duncan–Toor measurements
# exactly how (3) breaks: osmotic diffusion, reverse (uphill) diffusion and the
# diffusion barrier.

# %% [markdown]
# ## Background
#
# The review opens with four experimental demonstrations that a single
# composition-gradient-proportional flux law cannot describe multicomponent
# mixtures: the Duncan–Toor two-bulb experiment, the Vinograd–McBain mixed-ion
# cell, negative rejection in PEG/dextran ultrafiltration, and selectivity
# reversal in a zeolite membrane. This page works the first one, because it is
# an ideal gas at constant temperature and pressure; nothing exotic is present,
# and Fick's scalar law still fails.
#
# In that experiment (reproduced in full, against the measurements, on the
# related page [`A4.9`](../A4.9-duncan-toor/index.ipynb)), nitrogen sits between
# hydrogen and carbon dioxide, which exchange through a capillary. Reading the
# nitrogen trajectory in time, Toor's (1957) three named anomalies appear in
# sequence, and each corresponds to a statement about the *effective* Fick
# diffusivity $D_2 = -J_2/(c_t \nabla x_2)$ of nitrogen:
#
# | phenomenon | observation | required behaviour of $D_2$ |
# |---|---|---|
# | **osmotic diffusion** | flux without a gradient, $\nabla x_2 = 0,\ J_2 \neq 0$ | $D_2 \to \infty$ |
# | **reverse (uphill) diffusion** | flux against the gradient | $D_2 < 0$ |
# | **diffusion barrier** | no flux at the largest gradient, $J_2 = 0,\ \nabla x_2 \neq 0$ | $D_2 = 0$ |
#
# A parameter that must pass through $\pm\infty$, negative values and zero
# during one closed-cell experiment is not a transport property. The review's
# resolution is that the *pair* Maxwell–Stefan diffusivities stay positive,
# composition-independent and predictable, and all of the strangeness lives in
# the transformation to Fick form. This page makes that statement quantitative.
#
# **Scope.** The source is a 51-page review of the entire Maxwell–Stefan
# research programme. This page reproduces its bulk-fluid core: the
# $[B]$, $[\Gamma]$, $[D]$ formalism and its two worked examples, and the
# breakdown of the effective-diffusivity simplification. Later sections of the
# review (dusty gas, surface diffusion, electrolytes, membranes) belong to other
# catalog entries (`A4.3`, `H1.9`, `J5.x`).

# %% [markdown]
# ## The published model
#
# **Maxwell–Stefan** (review eqs. 16–17). At constant $T$ and $p$, the driving
# force on species $i$ is balanced by friction against every other species,
# with one inverse drag coefficient $\mathrm{D}_{ij}$ per pair:
#
# $$
# -\frac{x_i}{RT}\nabla_{T,p}\,\mu_i
# \;=\; \sum_{\substack{j=1\\ j\neq i}}^{n} \frac{x_j J_i - x_i J_j}{c_t\,\mathrm{D}_{ij}},
# \qquad i = 1,\dots,n
# $$
#
# (the friction term is frame-invariant, so molar diffusion fluxes $J_i$ can be
# used in place of the $N_i$ of the review's eq. 17). For an ideal gas the left
# side is just $-\nabla x_i$ (eq. 36).
#
# **Generalized Fick** (eqs. 24–28). Writing the driving force with the
# $(n{-}1)\times(n{-}1)$ thermodynamic-factor matrix
# $\Gamma_{ij} = \delta_{ij} + x_i\,\partial\ln\gamma_i/\partial x_j$ and
# eliminating species $n$ puts the same equations in matrix form,
#
# $$
# (J) = -c_t\,[B]^{-1}[\Gamma]\,(\nabla x) \;=\; -c_t\,[D]\,(\nabla x),
# \qquad [D] \equiv [B]^{-1}[\Gamma],
# $$
#
# with (eq. 26)
#
# $$
# B_{ii} = \frac{x_i}{\mathrm{D}_{in}} + \sum_{k\neq i}\frac{x_k}{\mathrm{D}_{ik}},
# \qquad
# B_{ij} = -x_i\!\left(\frac{1}{\mathrm{D}_{ij}} - \frac{1}{\mathrm{D}_{in}}\right)
# \quad (j \neq i).
# $$
#
# Note what "eliminating species $n$" costs: the $n{-}1$ retained species get
# explicit flux expressions and species $n$ is left to absorb the closure
# $\sum_i J_i = 0$. Which species is eliminated is pure bookkeeping: $[B]$,
# $[\Gamma]$ and $[D]$ all change with the choice, the fluxes do not. That
# invariance is tested on this page, and it is where the scalar simplification
# below parts company with Maxwell–Stefan.
#
# **This is an identity, not an approximation**: the two formulations are the
# same physics as long as $[D]$ keeps its full composition dependence. The
# review's point is about *parametrisation*: the $\mathrm{D}_{ij}$ are
# near-constant and predictable (for ideal gases, from kinetic theory; the
# second law bounds them, $\mathrm{D}_{ij} \geq 0$, eq. 23), while the four
# $D_{ij}$ elements swing with composition and $D_{12} \neq D_{21}$. In the
# binary limit the transformation collapses to $D = \mathrm{D}\,\Gamma$
# (eq. 13).
#
# **The simplification that is *not* an identity** (eqs. 29–32). Dropping the
# off-diagonal coupling gives one scalar per species,
# $J_i = -c_t D_{i,\mathrm{eff}} \nabla x_i$. The review lists when that is
# legitimate: all $\mathrm{D}_{ij}$ equal (eq. 29), or everything dilute in a
# solvent (eq. 30), or species $i$ diffusing through *stagnant* others: the
# condition behind Wilke's formula (eq. 32),
#
# $$
# D_{i,\mathrm{eff}} = \frac{1 - x_i}{\sum_{k \neq i} x_k / \mathrm{D}_{ik}} .
# $$
#
# The two-bulb experiment satisfies none of these (the pair diffusivities span a
# factor printed below the parameters, and every species counter-diffuses), which is what this page
# exploits.
#
# The scalar law is a statement about *individual* species, so it too has to be
# imposed on $n{-}1$ of them with the last one taking whatever the closure
# leaves, and unlike the matrix form, **the answer now depends on that
# choice**. The page measures the dependence.
#
# | Symbol | Code | Meaning |
# |---|---|---|
# | $\mathrm{D}_{ij}$ | `d_ms` | Maxwell–Stefan pair diffusivity |
# | $[B]$ | `build_b(...)` | inverted-drag matrix, eq. 26 |
# | $[\Gamma]$ | `gamma_p` | thermodynamic-factor matrix, eq. 24 |
# | $[D]$ | `fick_matrix(...)` | generalized Fick matrix, eq. 27 |
# | eqs. 16–17 direct | `friction_solve(...)` | the $n$ friction balances solved as they stand, no $[B]$ |
# | $D_{i,\mathrm{eff}}$ | `wilke_d(...)` | Wilke effective diffusivity, eq. 32 |
# | $J_i$ | `flux_*(...)` | molar diffusion flux |
# | $x_i$ | `x[:, i]` | mole fraction, capillary cell-centred |

# %% [markdown]
# ## Parameters and assumptions
#
# All numbers below were read off 600 dpi renders of the printed pages on
# 31 July 2026; the PDF is an Acrobat Capture OCR of a scan and its text layer
# cannot be trusted with signs or exponents.
#
# **Worked example 1** (p. 872, "Ideal ternary gas mixtures revisited"):
# H₂ (1)/N₂ (2)/CO₂ (3), $\mathrm{D}_{12} = 8.33$, $\mathrm{D}_{13} = 6.8$,
# $\mathrm{D}_{23} = 1.68 \times 10^{-5}$ m² s⁻¹ ("from the kinetic gas
# theory"; numerically identical to the values Duncan & Toor *state* for their
# own run at 35.2 °C and 1 atm). **Neither paper says these were measured**, and
# this page does not: Duncan & Toor (p. 40) call them "the best experimental
# values available (16) corrected to the thermostat temperature by the method
# suggested by Hirschfelder, Curtis, and Bird", i.e. literature values with a
# temperature correction applied, while Krishna & Wesselingh call the same three
# numbers kinetic-theory values. The `A4.9` page, which owns the dataset, makes
# the same distinction, and the provenance sidecar loaded below carries it.
# $[B]$ and $[D]$ are printed at the equilibrium composition stored in the CSV,
# together with a flux estimate for nitrogen.
#
# **Worked example 2** (p. 871, "Non-ideal ternary mixtures"): acetone
# (1)/benzene (2)/CCl₄ (3) at 25 °C and the composition stored in the CSV, with printed
# $[B]$, $[\Gamma]$ and $[D]$. One printed element is wrong, the sign of
# $B_{12}$, and the page proves it below from the paper's own numbers rather
# than repairing it silently.
#
# **Two-bulb demonstration**: geometry and conditions of the Duncan & Toor
# cell, taken from the published `A4.9` page where each value was verified
# against the original 1962 paper: bulbs 77.99 and 78.63 cm³, capillary 8.59 cm,
# 35.2 °C, 1 atm, effective $L/A = 258.1$ cm⁻¹. Initial compositions are the
# $t=0$ points of the digitised dataset. Note a bookkeeping trap: the review's
# eq. (3) numbers **swap the bulb labels** relative to Duncan & Toor's Fig. 2:
# in the review bulb 1 holds N₂ + CO₂, in the original (and on `A4.9`) bulb 1
# is the H₂ + N₂ charge. This page follows the original's labelling.
#
# **Assumptions** (all the review's): ideal gas, so $[\Gamma] = [I]$ for the
# transport solve; isothermal, isobaric; closed cell, so $\sum_i N_i = 0$
# (equimolar exchange) and $J_i = N_i$; capillary quasi-steady (its volume is
# small against the bulbs; the ratio is printed below); bulbs well mixed; no viscous flow, no thermal diffusion.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "main"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A4.2-maxwell-stefan-vs-fick"  # equals the page directory

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
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from pymrm.checks import find_roots, observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims

PAGE_A49 = "A4.9-duncan-toor"          # cross-page dataset source
claims = Claims(PAGE.split("-")[0])
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

# %%
# --- Conditions and cell geometry --------------------------------------------
# Geometry/conditions of the Duncan & Toor cell, from the published A4.9 page
# (verified there against the 1962 paper).
T      = 308.35        # K       35.2 C water bath
p      = 101.325e3     # Pa      1 atm
R_GAS  = 8.314462      # J/(mol K)
c_t    = p / (R_GAS * T)          # mol/m3

length = 85.9e-3       # m       capillary length, 8.59 cm
# BULB LABELLING. Everything with a bulb NUMBER on this page follows Duncan &
# Toor (bulb 1 = the H2 + N2 charge), which is also the labelling of the A4.9
# dataset plotted below. The review's eq. (3) numbers the bulbs the other way
# round; the one cell that uses them names its variables by CHARGE
# (x_h2n2 / x_n2co2) instead of by number, and no volume enters that cell.
vol_1  = 77.99e-6      # m3      bulb 1 = H2 + N2 charge  (Duncan & Toor)
vol_2  = 78.63e-6      # m3      bulb 2 = N2 + CO2 charge (Duncan & Toor)
LA_EFF = 258.1e2       # 1/m     effective L/A measured from seven binary runs
area   = length / LA_EFF          # m2, effective

# --- Species: 0 = H2, 1 = N2, 2 = CO2 (the review's 1, 2, 3) ----------------
SPECIES = ["H$_2$", "N$_2$", "CO$_2$"]
n_c, n_i = 3, 2                   # species, independent mole fractions

# Maxwell-Stefan pair diffusivities as printed in the review (p. 872), m2/s.
# Numerically identical to the 0.833 / 0.680 / 0.168 cm2/s Duncan & Toor STATE
# (p. 40) - which they describe as the best experimental values available in the
# literature, corrected to their thermostat temperature by the Hirschfelder-
# Curtis-Bird method, not as anything they measured. The review calls the same
# three numbers kinetic-theory values. Neither provenance is "measured here".
d_ms = np.zeros((3, 3))
d_ms[0, 1] = d_ms[1, 0] = 8.33e-5      # H2 - N2
d_ms[0, 2] = d_ms[2, 0] = 6.8e-5       # H2 - CO2
d_ms[1, 2] = d_ms[2, 1] = 1.68e-5      # N2 - CO2

t_end = 16 * 3600.0    # s, the span of Duncan & Toor's Fig. 2
n_z   = 40             # capillary cells

print(f"c_t = {c_t:.2f} mol/m3;  pair diffusivities span a factor "
      f"{d_ms[0,1]/d_ms[1,2]:.1f}; capillary volume / bulb volume = "
      f"{area * length / min(vol_1, vol_2):.2%}")

# %% [markdown]
# ## The data
#
# Two datasets, of different kinds, and the page keeps them apart:
#
# 1. **The review's own printed numbers** (tier 6, the authors' computations,
#    not measurements): the matrices $[B]$, $[\Gamma]$, $[D]$ of the two worked
#    examples plus the eq. (3) compositions, transcribed from 600 dpi renders
#    into `krishna-wesselingh-1997-worked-examples.csv`. One value is stored *as
#    printed* although it is provably a typo (the sign of $B_{12}$ in the
#    non-ideal example, settled below).
# 2. **The Duncan–Toor measurements** (tier 4, experimental): the digitised,
#    maintainer-reviewed bulb-composition dataset published with page `A4.9`,
#    loaded cross-page. Its provenance, extraction method and error estimate
#    live in that page's sidecar; nothing was re-digitised for this page. This
#    is the dataset every experimental number on this page is measured against,
#    so it is the one that sets the catalogue's `data.tier` for the page: 4,
#    with `own_tier: 6` recorded beside it (`docs/data-strategy.md` §3a).
#
#    What `A4.9` establishes about these rows, and whether it affects this page:
#
#    - bulb 1 is the H₂ + N₂ charge (the review numbers the bulbs the other way
#      round, see the note in *Parameters*): used throughout;
#    - the digitisation error stated in the sidecar (printed below) is the floor
#      under any agreement with these points: the Maxwell–Stefan agreement below
#      is reported as resolution-limited for that reason;
#    - the closure of the digitised compositions and the bulb-identity note in
#      the sidecar: the $t=0$ compositions are renormalised to sum to one before
#      they start the runs, and the $t=0$ points stay in the comparison;
#    - `A4.9` reads point O and the diffusion barrier off its step values and
#      refines its grid and time step: this page root-finds both features in
#      time and repeats both refinements (Validation), so its times are not
#      `A4.9`'s read-offs.

# %%
wx = load_data("krishna-wesselingh-1997-worked-examples.csv", page=PAGE)
wx_meta = load_meta("krishna-wesselingh-1997-worked-examples.csv", page=PAGE)
V = {(r.example, r.symbol): r.value for r in wx.itertuples()}   # value lookup

exp = load_data("duncan-toor-1962-run1.csv", page=PAGE_A49)
exp_meta = load_meta("duncan-toor-1962-run1.csv", page=PAGE_A49)

print("page dataset :", cite_data(wx_meta))
print("cross-page   :", cite_data(exp_meta))
print(f"\n{len(wx)} transcribed values; {len(exp)} experimental points "
      f"(error {exp_meta['acquisition']['estimated_error'].strip().splitlines()[0]}...)")

# %% [markdown]
# ## PyMRM implementation
#
# The transport model is the two-bulb quasi-steady capillary of the `A4.9` page
# (same grid, same operators, same Newton march) with the **flux closure made
# swappable**. Four closures share one skeleton:
#
# - **`flux_ms`**: Maxwell–Stefan: build $[B]$ at every capillary face from the
#   local composition and solve the $2\times2$ system
#   $(J) = -c_t [B]^{-1}(\nabla x)$ (structure `S9`: the flux is implicit).
# - **`flux_fick`**: generalized Fick: the *explicit closed-form* elements of
#   $[D](x) = [B]^{-1}$ (ideal gas, $[\Gamma]=[I]$), then
#   $(J) = -c_t[D](\nabla x)$.
# - **`flux_ms_friction`**: the same physics from the other end: the $n$
#   friction balances of eqs. (16)–(17) assembled **as they stand**, closed by
#   the cell constraint $\sum_k J_k = 0$, and solved for all $n$ fluxes at once.
#   No species is eliminated, no $[B]$ is ever formed, no code is shared with
#   `build_b`.
# - **`flux_wilke`**: scalar effective diffusivity, eq. (32), evaluated at the
#   local composition. Diagonal by construction: species $i$ responds to
#   $\nabla x_i$ only.
#
# **On what each comparison can and cannot show.** `flux_ms` and `flux_fick`
# call the *same* `build_b`; one applies `np.linalg.solve` to the $2\times2$
# matrix and the other its adjugate over its determinant. They must agree, and
# a mistake inside `build_b` (a wrong `d_ms` index, a dropped minus on an
# off-diagonal) would propagate identically into both. That comparison
# therefore confirms the closed-form inverse and nothing more, and the page
# labels it as the algebraic identity it is. **`flux_ms_friction` is the check
# with power**: it shares no line of code with `build_b`, so it does test
# eq. (26), the elimination of species $n$, and the transformation to
# $[D] = [B]^{-1}$.
#
# pymrm specifics: operators are assembled **once**, in `model.build_operators`.
# The bulb compositions enter through the `shapes_d` boundary matrices of
# `construct_grad`, and the end-face compositions the flux closures need come
# from `construct_boundary_value_matrices` with the same Dirichlet conditions, so
# no matrix changes during the march. `construct_div` uses `nu=0` (Cartesian: a straight uniform
# capillary). The `NumJac` stencil is `axes_diagonals=[0]` (each cell's residual
# depends on its neighbours through the face fluxes) plus `axes_blocks=[-1]`
# (all species couple pointwise through $[B]$); there is no second spatial axis,
# so a `[-2,-1]` block is not needed.

# %%
import model
from model import (full_x, build_b, fick_matrix, friction_solve, friction_residual,
                   fick_matrix_from_friction, wilke_d, grad_of, face_of, flux_ms,
                   flux_fick, flux_ms_friction, flux_wilke, solve_capillary, simulate)

# model.py keeps the page's parameters and operators as module state. The
# operators are assembled once in model.build_operators, whose docstring writes
# out both boundary conditions: x = x_bulb1 at z = 0 and x = x_bulb2 at z = L,
# a = 0, b = 1, d = 1 times the bulb composition passed through shapes_d.
model.configure(c_t=c_t, d_ms=d_ms, n_c=n_c, n_i=n_i, area=area, vol_1=vol_1,
                vol_2=vol_2, t_end=t_end, **model.build_operators(n_z, length))
print(f"{model.n_z} cells, field shape {model.shape}, operators assembled once")

# %% [markdown]
# ## Results
#
# ### The worked matrix examples, element by element
#
# **Example 1: ideal gas.** The review evaluates eq. (26) and $[D]=[B]^{-1}$
# for H₂/N₂/CO₂ at the equilibrium composition and prints both matrices. The
# same `build_b`/`fick_matrix` functions that drive the transport solve below
# must reproduce them, so this check tests the *implementation*, not only the
# transcription. The printed matrices are rounded to two or three significant
# figures and the review computed $[D]$ from its unrounded $[B]$, so the
# attainable agreement is set by that rounding; the most heavily rounded
# element, $B_{12}$, shows the largest deviation in the table below.

# %%
x_eq = np.array([V["ideal", "x1"], V["ideal", "x2"]])
b_print = np.array([[V["ideal", "B11"], V["ideal", "B12"]],
                    [V["ideal", "B21"], V["ideal", "B22"]]])
d_print = np.array([[V["ideal", "DF11"], V["ideal", "DF12"]],
                    [V["ideal", "DF21"], V["ideal", "DF22"]]])


def ideal_example(dms=d_ms):
    """Example 1 recomputed by eq. (26) and [D] = [B]^-1 from the pair diffusivities."""
    b = build_b(x_eq[None, :], dms)[0]
    d = fick_matrix(x_eq[None, :], dms)[0]
    return {"ex1_b_maxdev": float((np.abs(b - b_print) / np.abs(b_print)).max()),
            "ex1_fick_matrix_maxdev": float((np.abs(d - d_print) / np.abs(d_print)).max()),
            "ex1_d21_over_d22": float(d[1, 0] / d[1, 1]),
            "_b": b, "_d": d}


EX1 = ideal_example()
b_calc, d_calc = EX1["_b"], EX1["_d"]
print(f"equilibrium composition (CSV): x1 = {x_eq[0]}, x2 = {x_eq[1]}")
print("[B] (1e5 s/m2)      printed      recomputed    rel.dev")
print("[D] (1e-5 m2/s)")
for nm, Mp, Mc, s in [("B", b_print, b_calc, 1e5), ("D", d_print, d_calc, 1e-5)]:
    for i in range(2):
        for j in range(2):
            print(f"  {nm}{i+1}{j+1}   {Mp[i,j]/s:12.4f} {Mc[i,j]/s:12.4f}"
                  f"   {abs(Mc[i,j] - Mp[i,j]) / abs(Mp[i,j]):8.2%}")
print(f"\nmax relative deviation of the recomputed [B]: {EX1['ex1_b_maxdev']:.2%}, "
      f"of [D]: {EX1['ex1_fick_matrix_maxdev']:.2%}")
print(f"note the sign structure the review argues from: D21 = "
      f"{d_calc[1,0]*1e5:+.2f}e-5 is LARGER in magnitude than D22 = "
      f"{d_calc[1,1]*1e5:+.2f}e-5 (D21/D22 = {EX1['ex1_d21_over_d22']:+.2f}): nitrogen "
      f"answers hydrogen's gradient more strongly than its own.")

# %% [markdown]
# The review then estimates the initial nitrogen flux in its highlighted run
# from these matrices: $J_2 = -(c_t/\ell)\,(D_{21}\Delta x_1 + D_{22}\Delta x_2)$,
# with the bulb-minus-equilibrium differences stored in the CSV and printed
# below. Nitrogen's own difference $\Delta x_2$ is zero, so a scalar Fick law
# of *any* diffusivity gives $J_2 = 0$; the printed coefficients give a
# nitrogen flux of opposite sign to the hydrogen flux and of the size printed
# below: osmotic diffusion, from the cross-coefficient alone.

# %%
dx = np.array([V["ideal", "dx1"], V["ideal", "dx2"]])   # bulb - equilibrium
j_scaled_print = d_print @ dx        # J_i * (l/c_t), review's estimate, m2/s
j_scaled_calc = d_calc @ dx
ratio_print = float(j_scaled_print[1] / j_scaled_print[0])
print(f"bulb - equilibrium (CSV): dx1 = {dx[0]}, dx2 = {dx[1]}")
print(f"review's J2*(l/c_t) from the printed matrix : {abs(j_scaled_print[1])*1e5:.3f}e-5 m2/s")
print(f"same from the recomputed matrix             : {abs(j_scaled_calc[1])*1e5:.3f}e-5")
print(f"flux ratio J2/J1 at zero N2 driving force   : {ratio_print:+.3f}"
      f"  (recomputed {j_scaled_calc[1]/j_scaled_calc[0]:+.3f})")

# %% [markdown]
# **Example 2: non-ideal liquid, and a typo settled from the paper's own
# numbers.** For acetone/benzene/CCl₄ the review prints all three matrices of
# $[D] = [B]^{-1}[\Gamma]$, which makes the triple self-checking. The printed
# $B_{12}$ is negative, and it cannot be, on four independent grounds; the
# first never uses $B_{12}$ at all. Every number is printed by the cell below.
#
# 1. **Back-solve, $B_{12}$ never touched.** The *other three* printed elements
#    of $[B]$ determine all three diffusivities uniquely. Writing
#    $a = 1/\mathrm{D}_{12}$, $b = 1/\mathrm{D}_{13}$, $c = 1/\mathrm{D}_{23}$,
#    eq. (26) gives $B_{21} = -x_2(a-c)$, $B_{22} = x_1 a + (x_2+x_3)c$ and
#    $B_{11} = x_2 a + (x_1+x_3) b$, so
#    $$a = B_{22} - \frac{x_2+x_3}{x_2}B_{21},\qquad
#      c = a + \frac{B_{21}}{x_2},\qquad
#      b = \frac{B_{11} - x_2 a}{x_1 + x_3}.$$
#    The printed $B_{11}, B_{21}, B_{22}$ return the authors' *unrounded*
#    diffusivities, close to the printed two-figure ones. Feeding those back
#    into eq. (26) predicts a $B_{12}$ with the printed **magnitude** and the
#    opposite **sign**. Nothing in this route depends on any assumption about
#    $B_{12}$, or on the closure test of point 4.
# 2. eq. (26) with the printed (rounded) diffusivities gives
#    $B_{12} = -x_1(1/\mathrm{D}_{12} - 1/\mathrm{D}_{13})$, necessarily
#    positive because $\mathrm{D}_{12} > \mathrm{D}_{13}$.
# 3. **The review's other worked example disagrees with this one.** The ideal
#    gas example on p. 872 has the same ordering
#    $\mathrm{D}_{12} > \mathrm{D}_{13}$ and prints a positive $B_{12}$. Under
#    any single sign convention for eq. (26) the two printed examples are
#    mutually inconsistent, so one of them is a typo; and the ideal one is the
#    one this page reproduces element by element (above).
# 4. $[B]^{-1}[\Gamma]$ reproduces the printed $[D]$ with the $+$ sign, and
#    misses it with the printed $-$ sign. The rival hypothesis, that $B_{21}$
#    carries the typo instead, misses by more still.
#
# Nothing here involves anything outside the paper. The CSV stores the value as
# printed; the transformation below uses the corrected sign.

# %%
# Non-ideal example: acetone (1) / benzene (2) / CCl4 (3) at 25 C.
d_ms_liq = np.zeros((3, 3))
d_ms_liq[0, 1] = d_ms_liq[1, 0] = V["nonideal", "D12"]
d_ms_liq[0, 2] = d_ms_liq[2, 0] = V["nonideal", "D13"]
d_ms_liq[1, 2] = d_ms_liq[2, 1] = V["nonideal", "D23"]
x_liq = np.array([V["nonideal", "x1"], V["nonideal", "x2"], V["nonideal", "x3"]])

# eq. (26) recomputation from the printed diffusivities (independent of the printed [B]):
b_liq = np.zeros((2, 2))
for i in range(2):
    b_liq[i, i] = x_liq[i] / d_ms_liq[i, 2] + sum(
        x_liq[k] / d_ms_liq[i, k] for k in range(3) if k != i)
    for j in range(2):
        if j != i:
            b_liq[i, j] = -x_liq[i] * (1 / d_ms_liq[i, j] - 1 / d_ms_liq[i, 2])

b_liq_print = np.array([[V["nonideal", "B11"], V["nonideal", "B12"]],
                        [V["nonideal", "B21"], V["nonideal", "B22"]]])
gamma_p = np.array([[V["nonideal", "G11"], V["nonideal", "G12"]],
                    [V["nonideal", "G21"], V["nonideal", "G22"]]])
d_liq_print = np.array([[V["nonideal", "DF11"], V["nonideal", "DF12"]],
                        [V["nonideal", "DF21"], V["nonideal", "DF22"]]])


def gamma_dev(B):
    """Largest relative deviation of [B]^-1 [Gamma] from the printed [D]."""
    return float((np.abs(np.linalg.solve(B, gamma_p) - d_liq_print) / np.abs(d_liq_print)).max())


def nonideal_example(b_printed=b_liq_print):
    """Example 2 from a printed [B]: the back-solve from B11, B21, B22 (B12 unused),
    and [B]^-1 [Gamma] against the printed [D] under each sign hypothesis."""
    a_inv = b_printed[1, 1] - (x_liq[1] + x_liq[2]) * b_printed[1, 0] / x_liq[1]
    c_inv = a_inv + b_printed[1, 0] / x_liq[1]
    b_inv = (b_printed[0, 0] - x_liq[1] * a_inv) / (x_liq[0] + x_liq[2])
    b12_pred = -x_liq[0] * (a_inv - b_inv)
    fix = b_printed.copy(); fix[0, 1] = abs(fix[0, 1])            # this page's reading
    rival = b_printed.copy(); rival[1, 0] = -rival[1, 0]          # 'B21 is the typo'
    both = fix.copy(); both[1, 0] = -both[1, 0]                   # both negated
    return {"d12_backsolved": float(1 / a_inv), "d13_backsolved": float(1 / b_inv),
            "d23_backsolved": float(1 / c_inv), "b12_backsolved": float(b12_pred),
            "b12_backsolve_magnitude_dev": float(abs(abs(b12_pred) - abs(b_printed[0, 1]))
                                                 / abs(b_printed[0, 1])),
            "ex2_identity_maxdev": gamma_dev(fix),
            "ex2_as_printed_maxdev": gamma_dev(b_printed),
            "ex2_rival_b21_maxdev": gamma_dev(rival),
            "_both": gamma_dev(both), "_fix": fix, "_rival": rival}


EX2 = nonideal_example()
print("route 1 - back-solve the diffusivities from B11, B21, B22 ONLY (B12 unused):")
print(f"   D12, D13, D23 back-solved : {EX2['d12_backsolved']*1e9:.3f}, "
      f"{EX2['d13_backsolved']*1e9:.3f}, {EX2['d23_backsolved']*1e9:.3f} e-9 m2/s   "
      f"(printed {d_ms_liq[0,1]*1e9:.1f}, {d_ms_liq[0,2]*1e9:.1f}, {d_ms_liq[1,2]*1e9:.1f})")
print(f"   eq. (26) then predicts B12 = {EX2['b12_backsolved']/1e9:+.4f}e9 vs printed "
      f"{b_liq_print[0,1]/1e9:+.3f}e9 -> magnitude to "
      f"{EX2['b12_backsolve_magnitude_dev']:.1%}, sign opposite")
print(f"route 2 - eq. (26) with the printed diffusivities: B12 = "
      f"{b_liq[0,1]/1e9:+.3f}e9   (positive because D12 > D13)")
print(f"route 3 - the ideal example (p. 872), same D12 > D13 ordering, prints "
      f"B12 = {V['ideal','B12']/1e5:+.3f}e5 - positive\n")
print("route 4 - which hypothesis reproduces the printed [D]?")
for label, dev in [(f"as printed, B12 = {b_liq_print[0,1]/1e9:+.3f}e9", EX2["ex2_as_printed_maxdev"]),
                   (f"B12 = {EX2['_fix'][0,1]/1e9:+.3f}e9 (this page)", EX2["ex2_identity_maxdev"]),
                   ("[B] recomputed from eq. (26)", gamma_dev(b_liq)),
                   (f"rival: B21 = {EX2['_rival'][1,0]/1e9:+.3f}e9 instead", EX2["ex2_rival_b21_maxdev"]),
                   ("both off-diagonals negated", EX2["_both"])]:
    print(f"   {label:<32}: [B]^-1[Gamma] vs printed [D], max dev {dev:8.2%}")
print(f"\nOnly the corrected sign closes the printed triple, to "
      f"{EX2['ex2_identity_maxdev']:.2%} per element, i.e. [D] = [B]^-1 [Gamma] holds")
print("between the paper's own printed matrices. The residual is print rounding:")
print("the matrices carry two significant figures, and the back-solved "
      "diffusivities above show the authors used unrounded ones.")

# %% [markdown]
# ### The two-bulb experiment under four flux closures
#
# Now the same pymrm skeleton is run four times, Maxwell–Stefan via $[B]$,
# generalized Fick with the explicit composition-dependent $[D](x)$,
# Maxwell–Stefan from the raw friction balances, and the Wilke scalar, from the
# measured $t=0$ compositions. Lines are model; markers are the digitised
# measurements from `A4.9`. Nothing is fitted anywhere on this page.
#
# **The time march.** The capillary is quasi-steady, so time enters only
# through the two bulb balances, $V_1 c_t\,\mathrm{d}x_1/\mathrm{d}t = -A J(0)$
# and $V_2 c_t\,\mathrm{d}x_2/\mathrm{d}t = +A J(L)$. `simulate` steps them by
# backward Euler, with the end-of-step bulb compositions and the capillary
# profile coupled by three fixed-point iterations per step, and records the
# quasi-steady capillary flux at every step time, $t = 0$ included. Its order
# in time is measured in the Validation section.
#
# **Reading the model at the measurement times.** The model is known at the
# step times only; it is read at a measurement time by linear interpolation
# between the two neighbouring steps. The step is printed below; the two
# features read off in time (point O and the barrier) are roots of a cubic
# spline through the step values, found with `pymrm.checks.find_roots`. Backward
# Euler is first order, so each feature is read at 400 and 800 steps and
# extrapolated to zero step as $2r_{800} - r_{400}$; the Validation section
# measures the order of these readings.

# %%
# Initial compositions: the t=0 points of the digitised dataset (bulb 1 = H2+N2).
t0 = exp[exp.time_h == 0]
x1_0 = np.array([float(t0[(t0.bulb == 1) & (t0.species == s)].x.iloc[0])
                 for s in ["H2", "N2", "CO2"]])
x2_0 = np.array([float(t0[(t0.bulb == 2) & (t0.species == s)].x.iloc[0])
                 for s in ["H2", "N2", "CO2"]])
x1_0, x2_0 = x1_0 / x1_0.sum(), x2_0 / x2_0.sum()

idx = {"H2": 0, "N2": 1, "CO2": 2}
PERMS = {"CO2": (0, 1, 2), "H2": (1, 2, 0), "N2": (0, 2, 1)}   # last = eliminated
exp_n2 = exp[exp.species == "N2"].pivot(index="time_h", columns="bulb",
                                        values="x").dropna()
dx_exp = (exp_n2[1] - exp_n2[2]).values
T_END_H = t_end / 3600.0


def dev_vs_exp(h1, h2, t):
    """|model - measured| at every digitised point (model linearly interpolated
    between its step times), mole-fraction units."""
    d = [np.interp(r.time_h, t, (h1 if r.bulb == 1 else h2)[:, idx[r.species]]) - r.x
         for r in exp.itertuples()]
    return np.abs(np.asarray(d))


def n2_signal(h1, h2, t):
    """N2 bulb-difference trajectory, and |model - measured| at the paired points."""
    dxn = h1[:, 1] - h2[:, 1]
    return dxn, np.abs(np.interp(exp_n2.index, t, dxn) - dx_exp)


def run_ordering(perm, flux, dms=d_ms, reorder=True):
    """The whole transient with species perm[-1] eliminated. reorder=False keeps
    the pair diffusivities in the H2/N2/CO2 order (a bookkeeping defect)."""
    perm = np.asarray(perm)
    dms_p = dms[np.ix_(perm, perm)] if reorder else dms
    t, h1, h2, hJ = simulate(lambda x, a, b: flux(x, a, b, dms_p), x1_0[perm], x2_0[perm])
    p1, p2 = np.empty_like(h1), np.empty_like(h2)
    p1[:, perm], p2[:, perm] = h1, h2                # back to H2 / N2 / CO2 order
    return t, p1, p2, hJ                             # hJ stays in the permuted order


def one_root(f):
    """The single root of f on the run's time span; nan unless there is exactly one
    (a defect that removes the feature leaves no root, and the study goes on)."""
    r = find_roots(f, 0.0, T_END_H, n_scan=400)
    if len(r["roots"]) != 1 or r["failed"]:
        return float("nan")
    return float(r["roots"][0])


def bulb_features(t, h1, h2, hJ, dms):
    """The features of one march read in time: point O (zero N2 difference) and the
    barrier (zero N2 flux), roots of cubic splines through the step values, the N2
    difference at the barrier, and hydrogen's effective-diffusivity excursion outside
    its two pair values (flux and gradient at the same step time)."""
    dx = h1[:, 1] - h2[:, 1]
    s_dx, s_j2 = CubicSpline(t, dx), CubicSpline(t, hJ[:, 1])
    t_o, t_b = one_root(s_dx), one_root(s_j2)
    d1 = -hJ[:, 0] / (c_t * (h2[:, 0] - h1[:, 0]) / length)
    lo, hi = dms[0, 2], dms[0, 1]
    exc = max(max(lo - d1.min(), 0.0) / lo, max(d1.max() - hi, 0.0) / hi)
    return np.array([t_o, t_b, float(s_dx(t_b)) if np.isfinite(t_b) else np.nan, exc])


FEATURES = ("n2_point_o_time_h", "n2_barrier_time_h", "n2_barrier_dx", "h2_eff_excursion")


def paired(a, b):
    """Mean of b - a over paired points, its standard error (points taken as
    independent) and their ratio."""
    d = b - a
    se = d.std(ddof=1) / np.sqrt(len(d))
    if se <= 1e-12:                    # differences at round-off: no difference to rank
        return float(d.mean()), float(se), 0.0
    return float(d.mean()), float(se), float(d.mean() / se)


def closure_study(dms=d_ms, dms_check=None, nu=0, wilke_perm=(0, 1, 2)):
    """The two-bulb transient under the four closures and its comparison with the data.

    dms_check: pair diffusivities given to the explicit-[D] and friction routes
    only (a defect in those routes); nu: geometry of construct_div; wilke_perm:
    species order of the Wilke run, the last one eliminated."""
    dc = dms if dms_check is None else dms_check
    with model.operators(n_z, length, nu=nu):
        r = {"Maxwell-Stefan": simulate(lambda x, a, b: flux_ms(x, a, b, dms), x1_0, x2_0),
             "generalized Fick [D](x)": simulate(lambda x, a, b: flux_fick(x, a, b, dc),
                                                 x1_0, x2_0),
             "MS, friction form": simulate(lambda x, a, b: flux_ms_friction(x, a, b, dc),
                                           x1_0, x2_0),
             "Wilke effective D": run_ordering(wilke_perm, flux_wilke, dms)}

    def traj_diff(a, b):
        return float(max(np.abs(r[a][1] - r[b][1]).max(), np.abs(r[a][2] - r[b][2]).max()))

    t, h1, h2, hJ = r["Maxwell-Stefan"]
    _, w1, w2, _ = r["Wilke effective D"]
    out = {"ms_fick_algebraic_identity": traj_diff("Maxwell-Stefan", "generalized Fick [D](x)"),
           "ms_independent_equivalence": traj_diff("Maxwell-Stefan", "MS, friction form")}
    tot0, tot1 = vol_1 * x1_0 + vol_2 * x2_0, vol_1 * h1[-1] + vol_2 * h2[-1]
    out["mole_balance_error"] = float(np.abs(tot1 - tot0).max() / tot0.sum())
    dev_ms, dev_wk = dev_vs_exp(h1, h2, t), dev_vs_exp(w1, w2, t)
    dx_ms, n2_ms = n2_signal(h1, h2, t)
    dx_wk, n2_wk = n2_signal(w1, w2, t)
    # The six t = 0 points are the initial condition: their difference is exactly zero,
    # so the paired statistics run over the points after t = 0.
    after = (exp.time_h > 0).to_numpy()
    after_n2 = exp_n2.index.to_numpy() > 0
    out.update({"n_exp_points": len(dev_ms), "n_n2_points": len(n2_ms),
                "n_paired_points": int(after.sum()), "n_paired_n2_points": int(after_n2.sum()),
                "mean_abs_dev_ms": float(dev_ms.mean()), "mean_abs_dev_wilke": float(dev_wk.mean()),
                "n2_diff_err_ms": float(n2_ms.mean()), "n2_diff_err_wilke": float(n2_wk.mean())})
    for key, (a, b) in (("mean_abs_dev", (dev_ms[after], dev_wk[after])),
                        ("n2_diff_err", (n2_ms[after_n2], n2_wk[after_n2]))):
        m, se, z = paired(a, b)
        out.update({f"{key}_wilke_minus_ms": m, f"{key}_wilke_minus_ms_se": se,
                    f"{key}_wilke_minus_ms_z": z})
    out["wilke_n2_dx_above_zero"] = int((dx_wk > 0).sum())
    out["_wilke_n2_margin"] = float(np.abs(dx_wk).min())
    # Effective diffusivities backed out of the Maxwell-Stefan run (end-to-end gradient).
    d1_eff = -hJ[:, 0] / (c_t * (h2[:, 0] - h1[:, 0]) / length)
    d2_eff = -hJ[:, 1] / (c_t * (h2[:, 1] - h1[:, 1]) / length)
    # Point O, the barrier and hydrogen's excursion are step-dependent in a backward-Euler
    # march: read at 400 and 800 steps and extrapolated to zero step (first order).
    with model.operators(n_z, length, nu=nu):
        fine = simulate(lambda x, a, b: flux_ms(x, a, b, dms), x1_0, x2_0, n_steps=2 * (len(t) - 1))
    f_main, f_fine = bulb_features(t, h1, h2, hJ, dms), bulb_features(*fine, dms)
    out.update(dict(zip(FEATURES, 2 * f_fine - f_main)))
    out.update({"_runs": r, "_dx_ms": dx_ms, "_dx_wk": dx_wk, "_dev_ms": dev_ms,
                "_dev_wk": dev_wk, "_d1_eff": d1_eff, "_d2_eff": d2_eff, "_t": t})
    return out


def rnd(v, bound):
    """A round-off-level result printed as its bound (the digits vary between BLAS kernels)."""
    return f"below {bound:.0e} (round-off)" if v <= bound else f"{v:.2e} (ABOVE its bound {bound:.0e})"


CL = closure_study()
runs, t_h = CL["_runs"], CL["_t"]
print(f"time step {t_h[1] - t_h[0]:.3f} h, {len(t_h) - 1} steps; every capillary solve "
      "converged (solve_capillary raises otherwise)")
print(f"\nMS via [B] vs generalized Fick [D](x) = [B]^-1 : {rnd(CL['ms_fick_algebraic_identity'], 1e-10)}")
print("   ALGEBRAIC IDENTITY. Both call the same build_b; this confirms the")
print("   closed-form 2x2 inverse, and cannot detect an error inside [B].")
print(f"MS via [B] vs the raw friction system         : {rnd(CL['ms_independent_equivalence'], 1e-10)}")
print("   INDEPENDENT IMPLEMENTATION. No shared code, no [B], no species")
print("   eliminated: this is what tests eq. (26) and eq. (27).")

# %%
COL = {"H2": "tab:blue", "N2": "tab:red", "CO2": "tab:green"}
KEY = ["H2", "N2", "CO2"]
_, x1_ms, x2_ms, J_ms = runs["Maxwell-Stefan"]
_, x1_wk, x2_wk, _ = runs["Wilke effective D"]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
for ax, (ms, wk), b in zip(axes, [(x1_ms, x1_wk), (x2_ms, x2_wk)], (1, 2)):
    for i, sp in enumerate(KEY):
        ax.plot(t_h, ms[:, i], color=COL[sp], lw=2, label=SPECIES[i])
        ax.plot(t_h, wk[:, i], "--", color=COL[sp], lw=1.4, alpha=0.8,
                label="_nolegend_")
        sub = exp[(exp.bulb == b) & (exp.species == sp)]
        ax.plot(sub.time_h, sub.x, "o", color=COL[sp], ms=7, mfc="white",
                mew=1.8, zorder=5, label="_nolegend_")
    ax.set_xlabel("time, $t$ [h]")
    ax.set_title(f"Bulb {b}" + (": charged H$_2$ + N$_2$" if b == 1
                                else ": charged N$_2$ + CO$_2$"))
axes[0].set_ylabel("mole fraction, $x$ [-]")
axes[0].legend(frameon=False, loc="center right")
fig.suptitle("Maxwell–Stefan (solid) vs Wilke effective diffusivity (dashed) "
             "vs experiment (markers)")
fig.tight_layout()
plt.show()

# %%
# The nitrogen signal: its difference between the bulbs.
dx_ms, dx_wk = CL["_dx_ms"], CL["_dx_wk"]
t_zero, t_peak, dx_peak = CL["n2_point_o_time_h"], CL["n2_barrier_time_h"], CL["n2_barrier_dx"]

fig, ax = plt.subplots(figsize=(8.5, 4.8))
ax.plot(t_h, dx_ms, color="tab:red", lw=2.2, label="Maxwell–Stefan = Fick $[D](x)$")
ax.plot(t_h, dx_wk, "--", color="tab:gray", lw=2,
        label="Wilke effective $D$ (scalar law on H$_2$, N$_2$)")
ax.plot(exp_n2.index, dx_exp, "o", color="tab:red", ms=8, mfc="white",
        mew=1.8, zorder=5, label="experiment (A4.9 dataset)")
ax.axhline(0.0, color="k", lw=0.8)
ax.plot(t_zero, 0.0, "s", color="k", ms=7, zorder=6)
ax.plot(t_peak, dx_peak, "D", color="k", ms=7, zorder=6)
ax.annotate(f"point O: osmotic diffusion\nzero gradient, flux $\\neq$ 0\n"
            f"$t$ = {t_zero:.2f} h", xy=(t_zero, 0.0), xytext=(2.2, -0.055),
            arrowprops=dict(arrowstyle="->", lw=1.1), fontsize=9)
ax.annotate(f"diffusion barrier\nlargest gradient, flux = 0\n$t$ = {t_peak:.2f} h",
            xy=(t_peak, dx_peak), xytext=(9.2, 0.035),
            arrowprops=dict(arrowstyle="->", lw=1.1), fontsize=9)
ax.set_xlabel("time, $t$ [h]")
ax.set_ylabel(r"$x_{\mathrm{N_2},1} - x_{\mathrm{N_2},2}$ [-]")
ax.set_title("Nitrogen crosses zero; a scalar law imposed on nitrogen cannot")
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
plt.show()

print(f"N2 difference: starts {dx_ms[0]:+.4f}, crosses zero at t = {t_zero:.3f} h "
      f"(the one root on [0, {T_END_H:.0f}] h); the N2 flux vanishes at t = {t_peak:.3f} h, "
      f"where the difference is {dx_peak:+.4f} (roots read at {len(t_h) - 1} and "
      f"{2 * (len(t_h) - 1)} steps, extrapolated to zero step; the curve drawn is the "
      f"{len(t_h) - 1}-step march, whose own largest value is {dx_ms.max():+.4f})")
print(f"Wilke closure with the scalar law imposed on H2 and N2 (CO2 taking the "
      f"closure):\n   from {dx_wk[0]:+.4f} to {dx_wk[-1]:+.4f}, "
      f"{CL['wilke_n2_dx_above_zero']} of {len(dx_wk)} step values above zero "
      f"(closest approach to zero {CL['_wilke_n2_margin']:.4f}): no crossing, no barrier.")
print("Structural reason: a positive scalar diffusivity moves species i strictly "
      "down\nits OWN gradient, so the bulb difference of any species CARRYING the "
      "scalar law can\nonly decay towards zero. Which species carry it is a "
      "bookkeeping choice: next cell.")

# %% [markdown]
# ### Which species carries the scalar law, and why it matters
#
# The dashed curve above is *the* Wilke curve only for one bookkeeping choice.
# `n_i = 2` imposes eq. (32) on H₂ and N₂, and leaves CO₂ to absorb
# $\sum_i J_i = 0$. Nothing physical picks that ordering, and the Reuse section
# below invites you to change it, so the honest thing is to measure what it
# costs.
#
# Below, the whole calculation is repeated with each species in turn taking the
# dependent role. Same grid, same time steps, same data, same $t=0$
# compositions; only the labelling changes. Maxwell–Stefan should be exactly
# invariant (the elimination is undone by the closure), and the scalar law
# should not be.

# %%
def ordering_study(dms=d_ms, reorder=True):
    """Every closure run with each species in turn eliminated."""
    res = {}
    for elim, perm in PERMS.items():
        for tag, fl in (("MS", flux_ms), ("Wilke", flux_wilke)):
            t_p, p1, p2, _ = run_ordering(perm, fl, dms, reorder)
            dxn, n2e = n2_signal(p1, p2, t_p)
            res[(elim, tag)] = (float(dev_vs_exp(p1, p2, t_p).mean()), float(n2e.mean()), dxn, t_p)
    ms_means = [res[(e, "MS")][0] for e in PERMS]
    ms_n2 = [res[(e, "MS")][1] for e in PERMS]
    return {"ms_elimination_spread": float(max(max(ms_means) - min(ms_means),
                                                max(ms_n2) - min(ms_n2))),
            "mean_abs_dev_wilke_elim_h2": res[("H2", "Wilke")][0],
            "mean_abs_dev_wilke_elim_n2": res[("N2", "Wilke")][0],
            "n2_diff_err_wilke_elim_h2": res[("H2", "Wilke")][1],
            "n2_diff_err_wilke_elim_n2": res[("N2", "Wilke")][1],
            "_res": res}


OR = ordering_study()
perm_res = OR["_res"]


def _fmt(key):
    mean, n2e, dxn = perm_res[key][0] * 100, perm_res[key][1] * 100, perm_res[key][2]
    cross = "yes" if dxn.min() < 0 < dxn.max() else "no "
    return f"{mean:5.2f}  {n2e:6.2f}   {cross} {dxn.max():+.4f}"


print("deviations vs the A4.9 measurements, mole % (model - measured)\n")
print("eliminated |        Maxwell-Stefan          |         Wilke scalar")
print(" species   | mean   N2 diff  cross  max     | mean   N2 diff  cross  max")
print("-" * 74)
for elim in PERMS:
    star = " *" if elim == "CO2" else "  "
    print(f"  {elim:<4}{star}   | {_fmt((elim, 'MS'))}   | {_fmt((elim, 'Wilke'))}")
print("-" * 74)
print("* the ordering used everywhere else on this page; 'max' is the largest step value")

wk_means = [perm_res[(e, "Wilke")][0] * 100 for e in PERMS]
wk_n2 = [perm_res[(e, "Wilke")][1] * 100 for e in PERMS]
print(f"\nMaxwell-Stefan: the three orderings agree to {rnd(OR['ms_elimination_spread'], 1e-10)} in "
      "mole fraction: the\n   elimination is exact bookkeeping, as eq. (27) "
      "requires. Invariance is the claim,\n   and it is checked, not assumed.")
print(f"Wilke: mean deviation {min(wk_means):.2f} - {max(wk_means):.2f} mole %, "
      f"N2 signal {min(wk_n2):.2f} - {max(wk_n2):.2f} mole %, depending on nothing")
print("   but which species you chose to eliminate. And when N2 is the eliminated "
      "species\n   its flux is minus the sum of two scalar-law fluxes, so it CAN "
      "cross zero: the\n   'barrier' then costs nothing structurally, it is an "
      "artefact of the bookkeeping.")

# %%
fig, ax = plt.subplots(figsize=(8.5, 4.6))
STY = {"CO2": ("tab:gray", "--"), "H2": ("tab:orange", "-."), "N2": ("tab:purple", ":")}
for elim in PERMS:
    t_p = perm_res[(elim, "MS")][3]
    ax.plot(t_p, perm_res[(elim, "MS")][2], color="tab:red", lw=2.4, alpha=0.55,
            label="Maxwell–Stefan (all three orderings)" if elim == "CO2" else None)
    col, ls = STY[elim]
    ax.plot(t_p, perm_res[(elim, "Wilke")][2], ls, color=col, lw=1.8,
            label=f"Wilke, {elim} eliminated")
ax.plot(exp_n2.index, dx_exp, "o", color="k", ms=7, mfc="white", mew=1.6,
        zorder=5, label="experiment (A4.9 dataset)")
ax.axhline(0.0, color="k", lw=0.8)
ax.set_xlabel("time, $t$ [h]")
ax.set_ylabel(r"$x_{\mathrm{N_2},1} - x_{\mathrm{N_2},2}$ [-]")
ax.set_title("The scalar closure depends on which species is eliminated; "
             "Maxwell–Stefan does not")
ax.set_ylim(-0.09, 0.30)                       # headroom for the legend
ax.legend(frameon=False, fontsize=8.5, loc="upper center", ncol=2)
fig.tight_layout()
plt.show()

# %% [markdown]
# ### What the 'effective diffusivity of nitrogen' would have to do
#
# Back out $D_{2,\mathrm{eff}} = -J_2/(c_t\nabla x_2)$ from the Maxwell–Stefan
# solution (the quasi-steady capillary flux and the end-to-end gradient, both at
# the same step time). This is the quantity a Fickian fit to the experiment
# would report, and it is the review's Fig. 4 argument made quantitative: it
# diverges at the osmotic point, is negative through the reverse-diffusion
# window, and crosses zero at the barrier. Hydrogen's effective diffusivity,
# for contrast, stays close to its two pair values throughout (the excursion is
# printed below), which is why binary intuition survives for hydrogen and dies
# for nitrogen.

# %%
d1_eff, d2_eff = CL["_d1_eff"], CL["_d2_eff"]

fig, ax = plt.subplots(figsize=(8.5, 4.8))
ax.plot(t_h, d2_eff * 1e5, color="tab:red", lw=2, label=r"N$_2$: $-J_2/(c_t\nabla x_2)$")
ax.plot(t_h, d1_eff * 1e5, color="tab:blue", lw=2, label=r"H$_2$: $-J_1/(c_t\nabla x_1)$")
for v, lab in [(d_ms[0, 1] * 1e5, r"$\mathrm{D}_{12}$"),
               (d_ms[1, 2] * 1e5, r"$\mathrm{D}_{23}$")]:
    ax.axhline(v, color="k", lw=0.7, ls=":")
    ax.text(15.4, v * 1.15, lab, fontsize=9)
ax.axhline(0, color="k", lw=0.8)
ax.set_yscale("symlog", linthresh=1.0)
ax.set_ylim(-300, 300)
ax.axvline(t_zero, color="k", lw=0.7, ls="--")
ax.text(t_zero + 0.15, -120, "osmotic point:\n$D_{2}\\to\\pm\\infty$", fontsize=9)
ax.text(3.0, -25, "reverse diffusion:\n$D_{2}<0$", fontsize=9)
ax.annotate("barrier: $D_{2}=0$", xy=(t_peak, 0), xytext=(t_peak + 1.3, 30),
            arrowprops=dict(arrowstyle="->", lw=1.0), fontsize=9)
ax.set_xlabel("time, $t$ [h]")
ax.set_ylabel(r"effective Fick $D_i$ [$10^{-5}$ m$^2$s$^{-1}$], symlog")
ax.set_title("A 'diffusion coefficient' that visits $-\\infty$, 0 and $+\\infty$ "
             "in one experiment")
ax.legend(frameon=False, loc="upper right")
fig.tight_layout()
plt.show()

print(f"D1_eff (H2) range: {d1_eff.min()*1e5:.2f}e-5 ... {d1_eff.max()*1e5:.2f}e-5 m2/s, "
      f"positive and finite throughout\n   (this {len(t_h) - 1}-step march); extrapolated to zero step it strays "
      f"{CL['h2_eff_excursion']:.1%}\n   outside its own pair values ({d_ms[0,2]*1e5:.2f} to {d_ms[0,1]*1e5:.2f}e-5), so "
      "binary intuition survives for hydrogen")

# %% [markdown]
# ### The review's frozen-matrix flux estimate, tested
#
# The review's initial-flux estimate froze $[D]$ at the equilibrium composition.
# The full solve resolves $[D](x)$ along the capillary, where the composition
# swings across the whole ternary range. Comparing the *flux ratio*
# $J_2/J_1$ at $t=0$ (a convention-free number: the $\ell$ and $\Delta x$
# bookkeeping cancels) for the review's own eq. (3) run measures what the
# frozen-matrix shortcut costs. The same steady capillary carries the
# flux-uniformity check, the Wilke closure test and the grid ladder of the
# Validation section, so `steady_study` computes them together.

# %%
# The review's highlighted run, its eq. (3). CAREFUL: the review's "bulb 1" is
# the N2 + CO2 charge, i.e. Duncan & Toor's (and this page's) bulb 2. The two
# ends are therefore named by CHARGE here, not by number. Only the steady
# capillary enters, so no bulb volume is used; and the reported ratio J2/J1 is
# invariant under swapping the two ends anyway.
x_n2co2 = np.array([V["dt_run", "x1_bulb1"], V["dt_run", "x2_bulb1"],
                    V["dt_run", "x3_bulb1"]])      # review's bulb 1
x_h2n2 = np.array([V["dt_run", "x1_bulb2"], V["dt_run", "x2_bulb2"],
                   V["dt_run", "x3_bulb2"]])       # review's bulb 2
NZ = (10, 20, 40, 80, 160)                         # nested grids of the grid ladder


def capillary(xa, xb, flux, tol=1e-12):
    """Steady capillary between compositions xa (z = 0) and xb (z = L) on the
    current operators: cell profile and face fluxes."""
    x0 = np.linspace(0, 1, model.n_z)[:, None] * (xb[:n_i] - xa[:n_i]) + xa[:n_i]
    xs = solve_capillary(flux, x0, xa, xb, tol)
    return xs, flux(xs, xa, xb)


def steady_study(dms=d_ms, nu=0, face_rule="mean", tol=1e-12, dms_residual=None):
    """The eq. (3) run's steady capillary: flux ratio and the frozen-matrix error,
    flux uniformity (also with D(H2-N2) a decade low), the Wilke closure defect,
    the friction balance the solve never used, and the grid ladder."""
    fl = lambda x, a, b: flux_ms(x, a, b, dms)                       # noqa: E731
    d_low = dms.copy(); d_low[0, 1] = d_low[1, 0] = dms[0, 1] / 10.0
    with model.operators(n_z, length, nu=nu, face_rule=face_rule):
        xs, j = capillary(x_n2co2, x_h2n2, fl, tol)
        _, j_low = capillary(x_n2co2, x_h2n2, lambda x, a, b: flux_ms(x, a, b, d_low), tol)
        g0, xf0 = grad_of(xs, x_n2co2, x_h2n2), face_of(xs, x_n2co2, x_h2n2)
    jm, jl = j.mean(axis=0), j_low.mean(axis=0)
    out = {"j2_over_j1_full": float(jm[1] / jm[0]),
           "flux_nonuniformity": float(np.abs(j - jm).max() / np.abs(j).max()),
           "structural_blind_lost_decade": float(np.abs(j_low - jl).max() / np.abs(j_low).max()),
           "lost_decade_j2_over_j1": float(jl[1] / jl[0])}
    out["frozen_matrix_ratio_error"] = (abs(out["j2_over_j1_full"] - ratio_print)
                                        / abs(out["j2_over_j1_full"]))
    # The scalar law on ALL THREE species, on this profile (grad x3 = -(grad x1 + grad x2)):
    # an isobaric closed cell needs sum_i J_i = 0.
    grad_all = np.column_stack([g0, -(g0[:, 0] + g0[:, 1])])
    jw = -c_t * wilke_d(xf0, dms, which=range(n_c)) * grad_all
    out["wilke_closure_defect"] = float(np.abs(jw.sum(axis=1)).max() / np.abs(jw).max())
    # The friction balance the independent route never used must vanish.
    j_full = friction_solve(full_x(xf0), g0, dms)
    out["unused_friction_balance"] = float(friction_residual(
        full_x(xf0), g0, j_full, dms if dms_residual is None else dms_residual))
    # Grid ladder: J1 and J2/J1 on nested grids, successive differences.
    q = {}
    for n in NZ:
        with model.operators(n, length, nu=nu, face_rule=face_rule):
            _, j_n = capillary(x_n2co2, x_h2n2, fl, tol)
        q[n] = np.array([j_n.mean(axis=0)[0], j_n.mean(axis=0)[1] / j_n.mean(axis=0)[0]])
    st = observed_orders(lambda n: q[n], NZ)
    rel = np.abs(q[NZ[-1]] - q[NZ[-2]]) / np.abs(q[NZ[-1]])
    worst = max(st["orders"], key=lambda o: abs(o - 2.0))        # farthest from theory
    out.update({"grid_order": float(worst),
                "grid_ladder_over_roundoff": float(rel.min() / tol),
                "grid_error_ratio_nz40": float(abs(q[40][1] - st["extrapolated"][1])
                                               / abs(st["extrapolated"][1])),
                "_orders": st["orders"], "_j": j, "_j_low": j_low})
    return out


ST = steady_study()
print(f"J2/J1 at t=0, review's frozen-[D] estimate : {ratio_print:+.3f}")
print(f"J2/J1 at t=0, full composition-dependent   : {ST['j2_over_j1_full']:+.3f}")
print(f"frozen-matrix shortcut error in the ratio  : {ST['frozen_matrix_ratio_error']:.1%}")
print("\nWith the N2 driving force essentially zero, both give a nitrogen flux "
      "of opposite sign to\nthe hydrogen flux: osmotic diffusion. "
      f"The {ST['frozen_matrix_ratio_error']:.1%} shift is the composition dependence "
      "of [D]\nalong the capillary, which the review's estimate deliberately "
      "ignores and pymrm resolves.")

# %% [markdown]
# ## Validation
#
# Eight checks: two against the paper's printed numbers, five internal
# (3, 3b, 4, 5, 6), and one against the (cross-page) measurements. The
# deviation convention on this page is **model − measured**, reported in mole %
# (absolute mole-fraction difference × 100) for compositions and relative for
# matrix elements.
#
# Check 3 is the one that carries weight on the implementation: it compares
# $[B]$-based results against the raw friction system over the whole ternary
# composition triangle and over the whole transient, with no code in common.
# Check 3b is labelled for what it is: an algebraic identity between two ways
# of inverting the same $2\times2$ matrix, which would still agree if that
# matrix were wrong.
#
# **Checks 4 and 5's first number are one identity, and they test the solve,
# not the physics.** `solve_capillary` drives $\nabla\!\cdot\!J$ to its Newton
# tolerance with `construct_div(nu=0)`, so on a Cartesian grid the
# flux-uniformity number is the converged solve's residual restated; and with
# the flux uniform, $J_0 = J_L$, so the closed-cell mole balance closes
# **iff** the uniformity check passes. Neither can detect a wrong diffusivity, a
# wrong $(L/A)$ or a wrong bulb volume: the same uniformity computed with
# $\mathrm{D}_{\mathrm{H_2N_2}}$ a decade low stays at the round-off level
# while $J_2/J_1$ moves, both printed below. A wrong geometry (`nu = 1`) is
# what they catch; the defect table at the end shows it. Checks 3 and 6 carry
# the argument.
#
# **Numerical errors, split.** The grid order is measured on the steady
# capillary alone, where no time step enters; the time order on bulb-1 N₂ at
# three step times common to every level of a step ladder on the production
# grid, so the shared spatial error cancels in the successive differences. Each
# ladder carries a round-off companion: its finest successive difference over
# the Newton tolerance (accumulated over the steps for the time ladder).

# %%
def probe_study(dms_b=d_ms):
    """[D] by the [B] route (pair diffusivities dms_b) against the friction route
    (the configured ones) at 200 fixed random ternary compositions."""
    rng = np.random.default_rng(20260731)
    x_probe = rng.dirichlet(np.ones(n_c), size=200)[:, :n_i]
    d_b, d_f = fick_matrix(x_probe, dms_b), fick_matrix_from_friction(x_probe)
    return {"fick_matrix_independent_maxdev": float(np.max(np.abs(d_b - d_f) / np.abs(d_f)))}


NT = (50, 100, 200, 400, 800)                    # time steps over the 16 h run


def time_study(theta=1.0, tol=1e-12):
    """Bulb march on a step ladder: observed orders of bulb-1 N2 at t = 3.2, 8 and
    16 h (step times at every level), the N2 shift on halving the step, and the
    round-off companion."""
    h, feat = {}, {}
    for n in NT:
        t_n, h1_n, h2_n, hJ_n = simulate(flux_ms, x1_0, x2_0, n_steps=n, theta=theta, tol=tol)
        h[n] = h1_n[:, 1]
        feat[n] = bulb_features(t_n, h1_n, h2_n, hJ_n, d_ms)
    q = {n: h[n][[n // 5, n // 2, n]] for n in NT}
    st = observed_orders(lambda n: q[n], NT)
    # The features read in time on the last three levels: their order, and the
    # round-off companion (the finest difference over the steps times the tolerance).
    fo = observed_orders(lambda n: feat[n], NT[-3:])
    f_orders = np.asarray(fo["orders"]).ravel()
    f_diff = np.abs(feat[NT[-1]] - feat[NT[-2]])
    return {"feature_time_order": float(f_orders[np.argmax(np.abs(f_orders - 1.0))]),
            "feature_ladder_over_roundoff": float(f_diff.min() / (NT[-1] * tol)),
            "_feature_orders": f_orders,
            "time_order": float(max(st["orders"], key=lambda o: abs(o - 1.0))),
            "timestep_sensitivity": float(np.abs(h[800][::2] - h[400]).max()),
            "time_ladder_over_roundoff": float(np.abs(q[NT[-1]] - q[NT[-2]]).min()
                                               / (NT[-1] * tol)),
            "_orders": st["orders"]}


PR, TS = probe_study(), time_study()

print(f"1. ideal-gas [D] vs printed, max element dev   : {EX1['ex1_fick_matrix_maxdev']:.2%}")
print(f"2. non-ideal identity [B]^-1[Gamma] = [D]      : {EX2['ex2_identity_maxdev']:.2%} max")
print("3. independent friction route (no [B] anywhere):")
print(f"     [D] over 200 random ternary compositions   : "
      f"{rnd(PR['fick_matrix_independent_maxdev'], 1e-9)} rel")
print(f"     whole two-bulb transient                   : "
      f"{rnd(CL['ms_independent_equivalence'], 1e-10)} mole fraction")
print(f"     the friction balance NOT used in the solve : "
      f"{rnd(ST['unused_friction_balance'], 1e-12)} (must vanish)")
print(f"3b. algebraic identity, [B]-solve vs adjugate/det: {rnd(CL['ms_fick_algebraic_identity'], 1e-10)}")
print(f"4. species mole-balance error (closed cell)    : {rnd(CL['mole_balance_error'], 1e-10)}")
print(f"5. flux non-uniformity along capillary         : {rnd(ST['flux_nonuniformity'], 1e-10)}")
print(f"   same, D(H2-N2) a decade low                 : {rnd(ST['structural_blind_lost_decade'], 1e-10)}"
      f"   while J2/J1 goes {ST['j2_over_j1_full']:+.4f} -> {ST['lost_decade_j2_over_j1']:+.4f}")
print(f"   grid ladder n_z = {NZ}: observed orders "
      f"{', '.join(f'{p:.3f}' for p in ST['_orders'])}; J2/J1 at n_z = {n_z} is "
      f"{ST['grid_error_ratio_nz40']:.1e} (relative) from the Richardson limit; "
      f"finest difference / Newton tol = {ST['grid_ladder_over_roundoff']:.1e}")
print(f"   time ladder n_steps = {NT}: observed orders "
      f"{', '.join(f'{p:.3f}' for p in TS['_orders'])} (backward Euler: 1); "
      f"finest difference / (n_steps tol) = {TS['time_ladder_over_roundoff']:.1e}")
print(f"   point O, barrier time, barrier difference and hydrogen excursion read on the last three "
      f"levels {NT[-3:]}: observed orders {', '.join(f'{p:.2f}' for p in TS['_feature_orders'])} "
      f"(backward Euler: 1); finest difference / (n_steps tol) = "
      f"{TS['feature_ladder_over_roundoff']:.1e}")
print(f"   max N2 shift on halving dt at {len(t_h) - 1} steps     : {TS['timestep_sensitivity']:.2e}")
print(f"6. elimination invariance of Maxwell-Stefan    : {rnd(OR['ms_elimination_spread'], 1e-10)} "
      f"mole fraction over the three orderings\n"
      f"   same for the Wilke scalar closure           : "
      f"{(max(wk_means) - min(wk_means)) / 100:.1e}: the scalar answer is a bookkeeping choice")
print(f"7. vs experiment, mean |model - measured| over {CL['n_exp_points']} points: "
      f"Maxwell-Stefan {CL['mean_abs_dev_ms']*100:.2f} mole %, Wilke (CO2 eliminated) "
      f"{CL['mean_abs_dev_wilke']*100:.2f} mole %")
print(f"   paired statistics over the {CL['n_paired_points']} points after t = 0 (the t = 0 points "
      "are the initial condition, paired difference exactly 0)")
print(f"   paired difference Wilke - MS: {CL['mean_abs_dev_wilke_minus_ms']*100:.2f} "
      f"+/- {CL['mean_abs_dev_wilke_minus_ms_se']*100:.2f} mole % (standard error), "
      f"{CL['mean_abs_dev_wilke_minus_ms_z']:.1f} standard errors from zero")
print(f"   N2 bulb-difference error over {CL['n_n2_points']} paired times "
      f"({CL['n_paired_n2_points']} after t = 0 in the paired difference): "
      f"Maxwell-Stefan {CL['n2_diff_err_ms']*100:.2f} mole %, Wilke "
      f"{CL['n2_diff_err_wilke']*100:.2f} mole %; paired difference "
      f"{CL['n2_diff_err_wilke_minus_ms']*100:.2f} +/- {CL['n2_diff_err_wilke_minus_ms_se']*100:.2f}"
      f" mole %, {CL['n2_diff_err_wilke_minus_ms_z']:.1f} standard errors")
print("   The standard errors treat the points as independent; points along one")
print("   trajectory are correlated, so they probably understate the uncertainty.")
print("   This is the SAME computation A4.9 published (same solver, grid, step count and")
print("   dataset), so the Maxwell-Stefan mean is necessarily reproduced, not")
print("   independently corroborated. For scale, the digitisation error the A4.9")
print(f"   sidecar states: {exp_meta['acquisition']['estimated_error'].strip().splitlines()[0]}")
print("   So the Maxwell-Stefan agreement is resolution-limited and should not be read")
print("   as tighter than that.")
print(f"\n   Wilke closure defect: |sum_i J_i| reaches {ST['wilke_closure_defect']:.1%} of the "
      "largest flux: the scalar\n   law cannot even satisfy the equimolar "
      "constraint, which is why one species has to\n   be demoted to a dependent "
      "one in the first place.")

# %% [markdown]
# ### Claims
#
# Every number the page states, with how it is checked. Arithmetic on the
# printed matrices and every solve are held to `FIXED`: the solves are
# converged to a Newton step of 1e-12 in mole fraction, far inside it.
# Residuals, two-route gaps and orders are bounds with their basis.

# %%
FIXED = 1e-6
for k in ("ex1_fick_matrix_maxdev", "ex1_b_maxdev", "ex1_d21_over_d22"):
    claims.value(k, EX1[k], rtol=FIXED)
for k in ("ex2_identity_maxdev", "ex2_as_printed_maxdev", "ex2_rival_b21_maxdev",
          "b12_backsolve_magnitude_dev", "b12_backsolved", "d12_backsolved",
          "d13_backsolved", "d23_backsolved"):
    claims.value(k, EX2[k], rtol=FIXED)
claims.bound("fick_matrix_independent_maxdev", PR["fick_matrix_independent_maxdev"],
             max=1e-9, basis="identity",
             note="two linear-algebra routes to one matrix; round-off of order cond([B]) eps")
claims.bound("ms_independent_equivalence", CL["ms_independent_equivalence"], max=1e-10,
             basis="tolerance", note="two routes, each solve converged to a 1e-12 Newton step")
claims.bound("ms_fick_algebraic_identity", CL["ms_fick_algebraic_identity"], max=1e-10,
             basis="identity", note="[B] solve against its adjugate over determinant")
claims.bound("unused_friction_balance", ST["unused_friction_balance"], max=1e-12,
             basis="identity", note="sum_i grad x_i = 0 makes the n balances dependent")
claims.bound("mole_balance_error", CL["mole_balance_error"], max=1e-10, basis="tolerance",
             note="J(0) = J(L) to the capillary solve's residual closes the bulb balances")
for k in ("flux_nonuniformity", "structural_blind_lost_decade"):
    claims.bound(k, ST[k], max=1e-10, basis="tolerance",
                 note="div J = 0 makes J uniform; the spread is the Newton solve's residual")
for k in ("j2_over_j1_full", "lost_decade_j2_over_j1", "frozen_matrix_ratio_error",
          "wilke_closure_defect", "grid_error_ratio_nz40"):
    claims.value(k, ST[k], rtol=FIXED)
claims.bound("grid_order", ST["grid_order"], min=1.95, max=2.1, basis="theory",
             note="order farthest from 2 on the ladder; central face gradients and mean "
                  "face composition: second order")
claims.bound("grid_ladder_over_roundoff", ST["grid_ladder_over_roundoff"], min=1e3,
             basis="theory", note="finest relative successive difference over the Newton tolerance")
claims.bound("time_order", TS["time_order"], min=0.85, max=1.1, basis="theory",
             note="order farthest from 1 on the ladder; backward Euler: first order, "
                  "approached from below as the coarsest triple is pre-asymptotic")
claims.bound("time_ladder_over_roundoff", TS["time_ladder_over_roundoff"], min=1e3,
             basis="theory", note="finest successive difference over n_steps times the Newton tolerance")
claims.bound("feature_time_order", TS["feature_time_order"], min=0.85, max=1.15, basis="theory",
             note="order farthest from 1 over point O, barrier time, barrier difference and "
                  "hydrogen excursion read in time; backward Euler: first order")
claims.bound("feature_ladder_over_roundoff", TS["feature_ladder_over_roundoff"], min=1e3,
             basis="theory", note="finest successive difference of the four features over n_steps "
                                  "times the Newton tolerance")
claims.value("timestep_sensitivity", TS["timestep_sensitivity"], rtol=FIXED)
claims.bound("ms_elimination_spread", OR["ms_elimination_spread"], max=1e-10, basis="tolerance",
             note="three orderings of one physics, each solve converged to a 1e-12 Newton step")
for k in ("mean_abs_dev_wilke_elim_h2", "mean_abs_dev_wilke_elim_n2",
          "n2_diff_err_wilke_elim_h2", "n2_diff_err_wilke_elim_n2"):
    claims.value(k, OR[k], rtol=FIXED)
for k in ("n_exp_points", "n_n2_points", "n_paired_points", "n_paired_n2_points"):
    claims.exact(k, CL[k], structural=True, note="rows of the A4.9 CSV the comparison uses")
for k in ("mean_abs_dev_ms", "mean_abs_dev_wilke", "n2_diff_err_ms", "n2_diff_err_wilke",
          "mean_abs_dev_wilke_minus_ms", "mean_abs_dev_wilke_minus_ms_se",
          "n2_diff_err_wilke_minus_ms", "n2_diff_err_wilke_minus_ms_se",
          "n2_point_o_time_h", "n2_barrier_time_h", "n2_barrier_dx", "h2_eff_excursion"):
    claims.value(k, CL[k], rtol=FIXED)
for k in ("mean_abs_dev_wilke_minus_ms_z", "n2_diff_err_wilke_minus_ms_z"):
    claims.bound(k, CL[k], min=3.0, basis="chosen",
                 note="verdict: Wilke worse than Maxwell-Stefan by at least three standard errors")
claims.count("wilke_n2_dx_above_zero", CL["_dx_wk"], 0.0, above=True, min_margin=1e-3,
             note="Wilke N2 difference (CO2 eliminated) never reaches zero")
print(f"{len(claims.metrics)} claims registered")

# %% [markdown]
# ### The defect-injection table
#
# Every claim above needs a row that fails its own check. Each row injects one
# defect through a call argument (a pair diffusivity, a printed matrix element,
# the geometry of `construct_div`, the face-composition rule, the Newton
# tolerance, the time scheme, the species ordering) and reruns the whole study
# the claim comes from, so every claim the defect reaches is tested. A study
# that cannot be completed under a defect records `nan` for its claims, which
# fails them. The rows with all pair diffusivities equal make the Wilke scalar
# law exact by construction (the review's eq. 29).

# %%
BREAK, DEFINED = [], []


def brk(row, values):
    """Record one row: every claim in ``values`` against its own check."""
    DEFINED.append(row)
    for name, v in values.items():
        if name in claims.metrics:
            BREAK.append((row, name, claims.breaks(row, name, v)))


def run_row(row, fn, keys, **kw):
    """Rerun one study under a defect; nan for all its claims if it cannot finish."""
    try:
        out = fn(**kw)
    except Exception as exc:                     # noqa: BLE001 - a failed study is a result
        print(f"  {row}: study raised {type(exc).__name__}: {exc}; its claims record nan")
        out = {k: np.nan for k in keys}
    brk(row, {k: v for k, v in out.items() if not k.startswith("_")})


def keys(d):
    return [k for k in d if not k.startswith("_")]


def scaled(dms, i, j, f):
    """Pair diffusivities with D_ij (and D_ji) multiplied by f."""
    out = dms.copy(); out[i, j] = out[j, i] = dms[i, j] * f
    return out


d_swap = d_ms.copy()
d_swap[0, 2] = d_swap[2, 0] = d_ms[1, 2]; d_swap[1, 2] = d_swap[2, 1] = d_ms[0, 2]
d_equal = np.where(d_ms > 0, d_ms[0, 1], 0.0)
b_neg21 = b_liq_print.copy(); b_neg21[1, 0] = -b_neg21[1, 0]
WHOLE = {"steady capillary": (steady_study, keys(ST)),
         "two-bulb transient": (closure_study, keys(CL))}

run_row("ideal example with D(H2-CO2) and D(N2-CO2) swapped", ideal_example, keys(EX1),
        dms=d_swap)
run_row(f"non-ideal example with B21 read as {b_neg21[1, 0]:+.4g}", nonideal_example,
        keys(EX2), b_printed=b_neg21)
run_row("[B] route with D(H2-N2) x 1.1, friction route unchanged", probe_study, keys(PR),
        dms_b=scaled(d_ms, 0, 1, 1.1))
for row, kw in (("D(H2-N2) a decade low", dict(dms=scaled(d_ms, 0, 1, 0.1))),
                ("all pair diffusivities equal to D(H2-N2) (Wilke exact)", dict(dms=d_equal)),
                ("construct_div with nu = 1 (cylinder) instead of 0", dict(nu=1))):
    for what, (fn, ks) in WHOLE.items():
        run_row(f"{row} ({what})", fn, ks, **kw)
run_row("interior face composition from the left cell (first order)", steady_study,
        keys(ST), face_rule="left")
run_row("capillary Newton tolerance 1e-12 -> 1e-5 (grid ladder)", steady_study, keys(ST),
        tol=1e-5)
run_row("unused friction balance evaluated with D(N2-CO2) x 1.1", steady_study, keys(ST),
        dms_residual=scaled(d_ms, 1, 2, 1.1))
run_row("explicit [D] and friction routes with D(H2-N2) x 1.1, [B] route unchanged",
        closure_study, keys(CL), dms_check=scaled(d_ms, 0, 1, 1.1))
run_row("Wilke run with N2 eliminated (scalar law on H2 and CO2)", closure_study, keys(CL),
        wilke_perm=PERMS["N2"])
run_row("orderings run without reordering the pair diffusivities", ordering_study,
        keys(OR), reorder=False)
run_row("bulb march by the trapezoidal rule (theta = 0.5)", time_study, keys(TS), theta=0.5)
run_row("capillary Newton tolerance 1e-12 -> 1e-8 (time ladder)", time_study, keys(TS),
        tol=1e-8)

bt = pd.DataFrame(BREAK, columns=["row", "claim", "caught"])
caught = bt[bt.caught].groupby("row", sort=False).claim.apply(list)
assert len(DEFINED) == len(set(DEFINED)) == bt.row.nunique(), "every defined row is printed once"
print(f"{len(DEFINED)} rows defined, {bt.row.nunique()} in the table; claims each row fails:")
for row, names in caught.items():
    print(f"  {row}:")
    print("      " + "\n      ".join(textwrap.wrap(", ".join(names), 96)))
idle = [r for r in bt.row.unique() if r not in caught.index]
print(f"\nrows that fail no claim: {idle if idle else 'none'}")
missing = [n for n in claims.metrics
           if not claims.metrics[n].get("structural") and not bt[(bt.claim == n) & bt.caught].shape[0]]
print(f"claims no row fails: {missing if missing else 'none'}")

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises if any check fails.
claims.report();   # the semicolon keeps the notebook from echoing the payload

# %% [markdown]
# ## What pymrm adds
#
# **The transport solver is the one of the `A4.9` page**: same grid, operators
# and Newton march, so the Maxwell–Stefan deviation from the measurements is the
# number `A4.9` publishes, by construction. pymrm contributes the capillary
# solve (`construct_grad` with `shapes_d` for the moving bulb compositions,
# `construct_boundary_value_matrices` for the end-face compositions,
# `construct_div`, `NumJac`, `newton`) and the refinement and root tools of
# `pymrm.checks`. What this page adds beyond the solver sits at the formulation
# level; every number is printed above:
#
# - **The equivalence checked by a route that could have failed.** Eq. (27) says
#   Maxwell–Stefan and the generalized Fick matrix are the same physics. Coding
#   $[B]$ twice (a linear solve and its closed-form inverse) cannot show that,
#   because a mistake in $[B]$ survives both. So the page also solves the raw
#   friction balances with the $\sum_k J_k = 0$ bootstrap, eliminating nothing
#   and never forming $[B]$: it reproduces $[D](x)$ over the composition
#   triangle and the entire two-bulb transient to round-off. The often-repeated
#   "Maxwell–Stefan vs Fick" contest is really Maxwell–Stefan vs the *scalar*
#   shortcut.
# - **Why the scalar shortcut fails, not just by how much.** That it fails on
#   this experiment is `A4.9`'s published result; this page reaches it by a
#   different route, a full local-composition Wilke solve rather than a
#   frozen-$[D]$ analytic relaxation, and ranks the two closures by the paired
#   difference with its standard error. The new content is the structure of the
#   failure: (i) a positive scalar diffusivity moves a species strictly down its
#   own gradient, so the bulb difference of any species *carrying* the scalar
#   law decays monotonically and no barrier is possible; (ii) applied to all
#   three species at once the scalar fluxes violate the closed-cell constraint
#   $\sum_i J_i = 0$ by a large fraction of the largest flux, which is why one
#   species always has to be demoted; and (iii) the scalar answer *depends on
#   that demotion*, while Maxwell–Stefan is invariant to round-off. A model
#   whose prediction moves under a relabelling is not a transport model.
# - **The effective-diffusivity trace.** The review argues in a sketch (its
#   Fig. 4) that $D_{2}$ must visit $\pm\infty$, negative values and zero;
#   backing $-J_2/(c_t\nabla x_2)$ out of the pymrm solution turns that sketch
#   into a computed curve, with point O and the barrier root-found in time.
# - **A typo in the source settled from its own numbers, non-circularly.** The
#   printed $B_{12}$ of the non-ideal example is contradicted by a back-solve
#   that never uses it: the other three printed elements of $[B]$ fix the
#   authors' unrounded diffusivities, and eq. (26) then returns the printed
#   magnitude with the opposite sign. The review's *other* worked example prints
#   a positive $B_{12}$ under the same diffusivity ordering, the corrected sign
#   closes the printed triple, and the rival "$B_{21}$ is the typo" hypothesis
#   fails.
# - **The frozen-matrix estimate priced.** The review's constant-$[D]$
#   initial-flux estimate is compared with the composition-dependent solve in
#   the flux ratio; pymrm resolves the composition dependence along the
#   capillary that the estimate discards.
#
# What the page does **not** do: no new experimental data (the measurements are
# `A4.9`'s digitised dataset, used cross-page), no non-ideal *transport* solve
# (the $[\Gamma]$ example is verified as printed matrix algebra only), and
# nothing beyond the review's bulk-fluid sections: the dusty-gas and surface
# chapters belong to `A4.3` and `H1.9`.

# %% [markdown]
# ## Reuse
#
# **To adapt this page to your system**, change in order:
#
# | Change | Where |
# |---|---|
# | Different species or more of them | `SPECIES`, `d_ms`, `n_c`; `build_b`/`fick_matrix` generalize as written except the 2×2 closed-form inverse: use `np.linalg.inv` for $n>3$. `friction_solve` generalizes unchanged |
# | A different species order | free for Maxwell–Stefan, *not* for the scalar closure: `run_ordering` and the permutation table above show the cost. If you use `flux_wilke`, say which species carries the scalar law |
# | Non-ideal mixture | multiply the solved $[B]^{-1}$ by $[\Gamma](x)$ in `flux_ms` (the `nonideal` example's matrices show the pattern) |
# | Different geometry | `length`, `vol_1`, `vol_2`, `LA_EFF`, and `nu` in `model.build_operators` |
# | Porous medium | add Knudsen terms to $[B]$: the dusty gas model, page `A4.3` |
# | A quick screening decision | evaluate `fick_matrix` at your feed composition: if the off-diagonal elements rival the diagonal (here $D_{21}/D_{22}$ is printed in the first example), a scalar effective diffusivity will misbehave |
#
# The three `flux_*` closures against one skeleton are the reusable pattern:
# any `struct:S9` page can swap constitutive models this way and regression-test
# their agreement.
#
# **Related pages:** `A4.9` Duncan–Toor ternary diffusion (the experiment, and
# the source of the dataset used here) · `A4.3` dusty gas model · `A4.8` Stefan
# tube · `B1.11` multicomponent pellet · `H1.9` zeolite membrane permeation.
#
# ## References
#
# Krishna, R. and Wesselingh, J. A. (1997). *The Maxwell–Stefan approach to
# mass transfer.* Chemical Engineering Science **52**(6), 861–911.
# [doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7)
#
# Duncan, J. B. and Toor, H. L. (1962). *An experimental study of three
# component gas diffusion.* AIChE Journal **8**(1), 38–41.
# [doi:10.1002/aic.690080112](https://doi.org/10.1002/aic.690080112)
#
# Toor, H. L. (1957). *Diffusion in three-component gas mixtures.* AIChE
# Journal **3**(2), 198–207.
#
# Wilke, C. R. (1950). *Diffusional properties of multicomponent gases.*
# Chemical Engineering Progress **46**, 95–104.
#
# Taylor, R. and Krishna, R. (1993). *Multicomponent Mass Transfer.* Wiley.
