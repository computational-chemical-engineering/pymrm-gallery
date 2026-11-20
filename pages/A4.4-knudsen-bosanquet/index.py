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
# title: "Knudsen diffusion and the Bosanquet relation"
# description: "Where the Bosanquet interpolation 1/D = 1/D_K + 1/D_AB is right and where it is not: exact under Knudsen control, worst under bulk control, with the error set by the flux ratio."
# categories: [sec:A, struct:S3, tier:T0, data:tier6, phase:gas]
# date: 2026-08-01
# ---
#
# # Knudsen diffusion and the Bosanquet relation
#
# **Catalog ID:** `A4.4` · **Structure:** `S3` (scalar diffusion and reaction) ·
# **Tier:** T0
#
# Inside a catalyst pellet or a membrane, a molecule collides with other
# molecules *and* with the pore wall. Which collision dominates is set by the
# Knudsen number, mean free path over pore diameter, and the two limits have
# different physics: bulk diffusion couples the species to each other, Knudsen
# diffusion couples each species only to the wall.
#
# The universal engineering closure is the **Bosanquet relation**, adding the
# two resistances in series,
#
# $$
# \frac{1}{D_{\mathrm{eff}}} = \frac{1}{D_{K}} + \frac{1}{D_{AB}} .
# $$
#
# It is an interpolation, not an identity, and this page is about the
# difference. Krishna and Wesselingh's 1997 review derives it as a special case
# of the dusty gas model, names it, states the three conditions it needs, and
# says that it is "very restricted in its applicability", without saying how
# wrong it gets, or where.
#
# This page answers that. It reads the dusty-gas relations off the review's
# printed pages, implements them in pymrm, and measures the Bosanquet error as
# a function of Knudsen number and flux ratio. The answer is not the obvious
# one: **the error is worst in the bulk-diffusion limit and vanishes in the
# Knudsen limit**. The interpolation is exact precisely where the Knudsen term
# dominates, and wrong where it does not matter that it is there.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.20.3"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A4.4-knudsen-bosanquet"  # equals the page directory

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
# Imports. sympy is a Colab default.
import textwrap
import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from pymrm.checks import observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import (RGAS, Cell, d_groups, knudsen_D, pellet, isobaric_slab,
                   isobaric_flux_closed_form, eta_classical)

PAGE_A42 = "A4.2-maxwell-stefan-vs-fick"      # cross-page: the printed pair diffusivities
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# **The mechanisms.** The review distinguishes three transport mechanisms inside
# a porous particle (its Fig. 40): bulk or free-molecular diffusion, significant
# at large pore sizes and high pressure, "where molecule–molecule collisions
# dominate over molecule–wall collisions"; Knudsen diffusion, which "becomes
# predominant when the mean-free path of the molecular species is much larger
# than the pore diameter and hence molecule–wall collisions become important";
# and surface diffusion of adsorbed species, dominant in micropores. The review
# is explicit that bulk and Knudsen diffusion "occur together and it is prudent
# to take both mechanisms into account rather than assume that one or other
# mechanism is 'controlling'".
#
# **Why an interpolation is needed at all.** In the Knudsen limit each species
# moves independently: it exchanges momentum only with the wall, so its flux
# depends on its own partial-pressure gradient and on nothing else. In the bulk
# limit species exchange momentum with each other, so one species' flux depends
# on every other species' flux. Those are structurally different laws, and the
# transition between them is not a small correction: it is a change in what the
# flux is a function of. Bosanquet's relation asserts that the two resistances
# simply add. Whether that is right turns out to depend entirely on **what the
# other species is doing**.
#
# **Scope, and what belongs elsewhere.** This page is the *scalar* closure:
# `S3`, one effective diffusivity per species, the thing that goes into a Thiele
# modulus. The general $n$-component matrix form of the dusty gas model (the
# review's eqs. (86)–(87), its ternary He/Ne/Ar comparison against the Remick &
# Geankoplis capillary measurements in Fig. 44, and the viscous-flow term) is
# page [`A4.3`](../A4.3-dusty-gas-model/index.ipynb). Everything below is
# **binary**, which is the setting in which the Bosanquet relation is stated and
# used. Bulk multicomponent diffusion with no wall is page
# [`A4.2`](../A4.2-maxwell-stefan-vs-fick/index.ipynb), built from the same
# review; nothing on that page involves a pore.

# %% [markdown]
# ## The published model
#
# **Source, and which document each equation came from.** The Knudsen and
# Bosanquet results are pre-DOI: Knudsen's kinetic-theory treatment dates from
# 1909 and Bosanquet's interpolation from a 1944 British wartime report
# (BR-507), neither of which is on disk or reachable. **Neither was consulted.**
# Every equation below was read from **Krishna & Wesselingh (1997)**, which
# prints all of them with attribution and *names* the Bosanquet formula.
# Equation numbers throughout are that review's. Each was read off a 600 dpi
# render of the printed page: the PDF is an Acrobat 3.0 Capture OCR of a scan
# whose text layer is unreliable (it renders eq. (86)'s prefactor as $c_t/RT$
# where the page prints $1/RT$, and eq. (109)'s stoichiometric subscript $\nu_1$
# as $\nu_{1i}$).
#
# **The dusty gas model** (eq. 82). The pore wall is treated as a giant
# motionless pseudo-species, "dust", and the Maxwell–Stefan equations are
# written for the $(n{+}1)$-component mixture. Eliminating the dust gives, for
# $i = 1,\dots,n$,
#
# $$
# -\frac{1}{RT}\nabla p_i \;=\; \sum_{j=1}^{n}\frac{x_j N_i - x_i N_j}{\mathrm{D}^e_{ij}}
# \;+\; \frac{N_i}{\mathrm{D}^e_{iM}} .
# $$
#
# The first term is the Maxwell–Stefan friction against the other species; the
# second is friction against the wall. $\mathrm{D}^e_{ij}$ are the effective
# binary pair diffusivities, $\mathrm{D}^e_{iM}$ the effective Knudsen
# diffusivities. In matrix form (eqs. 86–87), $(N) = -[B^e]^{-1}\nabla(p)/RT$ with
# $B^e_{ii} = 1/\mathrm{D}^e_{iM} + \sum_{k\ne i} x_k/\mathrm{D}^e_{ik}$ and
# $B^e_{ij} = -x_i/\mathrm{D}^e_{ij}$.
#
# **The porous matrix** (eqs. 83–84): $\mathrm{D}^e_{ij} = (\varepsilon/\tau)\,
# \mathrm{D}_{ij}$, with $\tau = 1$ for a cylindrical pore.
#
# **The Knudsen diffusivity** (eq. 85), attributed to Jackson (1977) and Mason &
# Malinauskas (1983):
#
# $$
# \mathrm{D}^e_{iM} = \frac{\varepsilon}{\tau}\,\frac{d_0}{3}\sqrt{\frac{8RT}{\pi M_i}} ,
# $$
#
# with $d_0$ the pore **diameter** and, in the review's own words, "the
# square-root term represents the velocity of motion", the mean molecular
# speed $\bar v_i$. So $\mathrm{D}^e_{iM} = (\varepsilon/\tau)\,d_0\bar v_i/3$:
# proportional to pore size, independent of pressure, and scaling as
# $\sqrt{T}$ and $1/\sqrt{M_i}$.
#
# **What that implies** (stated on p. 887): the $\mathrm{D}^e_{ij}$ are
# inversely proportional to pressure and independent of pore size, so in the
# bulk-controlled regime the fluxes are independent of both; the
# $\mathrm{D}^e_{iM}$ are independent of pressure and proportional to pore size,
# so in the Knudsen-controlled regime the fluxes are directly proportional to
# both. Two testable claims.
#
# **Graham's law** (eq. 107). Imposing $\nabla p = 0$ on a gaseous mixture
# constrains the fluxes:
#
# $$
# \sum_{i=1}^{n} N_i \sqrt{M_i} = 0 \qquad (\nabla p = 0;\ \text{gaseous mixtures}).
# $$
#
# **The effective diffusivity** (eqs. 108–109). Defining $D_i$ by
# $N_i = -D_i(\,c_i\nabla_{T,p}\mu_i/RT + c_i \bar V_i \nabla p / RT + \dots)$,
# which for an ideal gas is $N_i = -D_i \nabla c_i$, and eliminating
# $\nabla p$ with the review's eq. (105) gives, for species 1,
#
# $$
# \frac{1}{D_1} = \frac{1}{\mathrm{D}^e_{1M}}
#  + \sum_{j=2}^{n}\frac{x_j}{\mathrm{D}^e_{1j}}\left(1 - \frac{x_1\nu_j}{x_j\nu_1}\right)
#  \;-\; \underbrace{\alpha'_1 x_1 \frac{\mathrm{D}_{\rm visc}}{\mathrm{D}_{1M}}
#    \frac{\sum_i (\nu_i/\nu_1)/\mathrm{D}^e_{iM}}
#         {1 + \mathrm{D}_{\rm visc}\sum_i \alpha'_i x_i/\mathrm{D}^e_{iM}}}_{\text{viscous flow}} .
# $$
#
# The $\nu_i$ are the stoichiometric coefficients: for a heterogeneous reaction
# the **flux ratios are fixed by the stoichiometry**, $N_j/N_1 = \nu_j/\nu_1$.
# The review's comment on this expression is the heart of the case: "The
# effective diffusivity thus defined will be a strong function of the
# composition and also the flux ratios of all the species participating in say a
# chemical reaction within the pellet." This page works at
# $\mathrm{D}_{\rm visc} = 0$ throughout (no viscous flow), for which the
# review's own remark is that "neglect of the viscous flow contribution … is not
# very serious; this result is typical (Haynes, 1978)".
#
# **The Bosanquet formula** (eq. 110). "For the special case of a (i) binary
# mixture, (ii) with no net change in the number of moles, and (iii) satisfying
# eq. (107), eq. (109) simplifies to"
#
# $$
# \frac{1}{D_1} = \frac{1}{\mathrm{D}^e_{1M}} + \frac{1}{\mathrm{D}^e_{12}},
# $$
#
# "a relation usually referred to as the Bosanquet formula. As noted above this
# formula is very restricted in its applicability."
#
# The review does not merely assert the interpolation: it prints the exact
# expression the interpolation approximates *and* the conditions under which the
# two coincide. Everything this page measures is the gap between eq. (109) and
# eq. (110).
#
# | Symbol | Code | Meaning |
# |---|---|---|
# | $\mathrm{D}^e_{iM}$ | `knudsen_D`, `dam`/`dbm` | effective Knudsen diffusivity, eq. 85 |
# | $\mathrm{D}^e_{12}$ | `dab` | effective binary pair diffusivity, eq. 83 |
# | $[B^e]$ | `Cell.b_matrix` | dusty-gas resistance matrix, eqs. 82/87 |
# | $D_1$ | `Cell.d_scalar(...)` | scalar effective diffusivity, eq. 109 |
# | $D_{\rm Bos}$ | `Cell.d_scalar(..., bosanquet=True)` | eq. 110 |
# | $\rho = -N_2/N_1$ | `rho` | the flux ratio eq. 109 depends on |
# | $\mathrm{Kn}$ | `kn` | $\mathrm{D}^e_{12}/\mathrm{D}^e_{1M}$, see below |
# | $\phi$ | `phi` | Thiele modulus of the pellet |

# %% [markdown]
# ## Parameters and assumptions
#
# **The Knudsen number used here.** The review never defines one. This page uses
#
# $$
# \mathrm{Kn} \;\equiv\; \frac{\mathrm{D}^e_{12}}{\mathrm{D}^e_{1M}},
# $$
#
# which is the ratio of the two resistances in the Bosanquet sum and is
# therefore the natural abscissa for everything below. It is *proportional* to
# the physical mean-free-path-over-diameter, and under the elementary
# kinetic-theory estimate $\mathrm{D}_{12} = \bar v\lambda/3$ it equals it,
# because eq. (85) reads $\mathrm{D}^e_{1M} = \bar v d_0/3$.
# $\mathrm{Kn} \to 0$ is bulk control, $\mathrm{Kn}\to\infty$ Knudsen control,
# and the two resistances are equal at $\mathrm{Kn} = 1$.
#
# **Dimensionless working.** Most results are reported in the groups
# $(\mathrm{Kn}, \rho, \phi)$ with concentrations scaled by the surface value
# and diffusivities by $D_{\rm Bos}$ (`d_groups`: $\mathrm{D}^e_{AM} = 1 +
# 1/\mathrm{Kn}$, $\mathrm{D}^e_{AB} = 1 + \mathrm{Kn}$). The review prints no
# pore sizes, porosities or rate constants for a reacting pellet, so inventing
# them would add nothing traceable. Where a dimensional number is wanted (the
# transition curve and the p. 887 scaling check), the page uses the review's
# **own printed** binary pair diffusivities from its worked example on p. 872 for
# the H₂/N₂/CO₂ pairs, which the review introduces as values that "can be
# estimated from the kinetic gas theory" for the Duncan & Toor (1962) system,
# loaded cross-page from the published `A4.2` dataset. The review prints **no
# temperature and no pressure** beside them; the 35.2 °C and 1 atm at which they
# are used here are Duncan & Toor's own experimental conditions, from page
# [`A4.9`](../A4.9-duncan-toor/index.ipynb)'s dataset sidecar, and they enter only
# through a mean molecular speed and a total concentration. Results item 5b
# measures what they are worth. Also used: $\varepsilon/\tau = 1$, which is
# eq. (84)'s printed $\tau = 1$ (a straight cylindrical pore) together with an
# **assumed** $\varepsilon = 1$, and standard molar masses.
#
# **Assumptions**, all of them the review's: isothermal ideal gas; no viscous
# flow ($\mathrm{D}_{\rm visc} = 0$, so the third term of eq. 109 is dropped);
# no surface diffusion; a single pore size; and, in the reacting pellet,
# first-order kinetics in A and a steady state, so the flux ratio is fixed by
# the stoichiometry.
#
# **The two settings.** The Bosanquet relation is used in two places, and its
# conditions fail differently in each:
#
# | setting | flux ratio $\rho = -N_2/N_1$ | which condition fails |
# |---|---|---|
# | isobaric binary counter-diffusion cell | $\sqrt{M_1/M_2}$ (eq. 107, Graham) | (ii), unless $M_1 = M_2$ |
# | pellet with $A \to \nu_B B$ | $\nu_B$ (stoichiometry) | (ii) **and** (iii), unless $\nu_B = 1$ |
#
# Conditions (ii) and (iii) are simultaneously satisfiable **only when the two
# species have the same molar mass**: equimolar counter-diffusion of unequal
# masses violates Graham's law and therefore generates a pressure gradient.
# That is the sense in which the formula is "very restricted".

# %%
# Numerical parameters (chosen; each is named beside the result it sets).
N_PELLET = 600                      # cells of the pellet sweep
TOL = 1e-12                         # Newton step tolerance relative to the surface concentration
KN_HI = 1e6                         # deep Knudsen control, for Jackson's limit
KN_TAIL = [1e5, 1e6, 1e7, 1e8]      # Kn ladder for the finite-Kn residue
NU_LIST = [1.0, 2.0, 3.0, 4.0]      # stoichiometries for Jackson's limit
PHI = 30.0                          # Thiele modulus of the sweep: the Bosanquet error depends on it
PHI_LIST = [10.0, 30.0, 100.0]      # Thiele moduli at which the bulk-limit error is quoted
KN_SWEEP = np.logspace(-3, 6, 19)   # Knudsen numbers of the pellet sweep
NS_HEAD = (300, 600, 1200, 2400, 4800)            # ladder for the quoted effectiveness errors
NS_ORDER = (50, 100, 200, 400, 800, 1600, 3200, 6400)   # ladder for the grid order
XA1, XA2 = 0.9, 0.1                 # mole fraction of A at the two faces of the slab
KN_GRID = np.array([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3])   # Kn of the isobaric slab
N_SLAB = 200                        # cells of the dimensionless slab
DELTA = 1.0e-3                      # m, thickness of the dimensional slab (sets the flux scale)
N_DIM = 120                         # cells of the dimensional slab
H_LOG = 1e-3                        # step in ln p and ln d0 for the local slopes

# %% [markdown]
# ## The data
#
# There is no dataset here in the experimental sense. **This page is tier 6.**
# Two collections of printed numbers are used, and neither is a measurement:
#
# 1. **`krishna-wesselingh-1997-porous-media.csv`**: the constants and stated
#    results of the review's dusty-gas sections (pp. 886–893), read off 600 dpi
#    renders: the $d_0/3$ and $M^{-1/2}$ of eq. (85), $\tau = 1$ of eq. (84),
#    the $\sqrt{M_i}$ of Graham's law, and the review's quotation of Jackson's
#    pellet-centre pressure with its "40 % increase" for $\nu_B = 2$.
# 2. **The review's own binary pair diffusivities** (p. 872), loaded cross-page
#    from the published `A4.2` dataset (rows `ideal`, `D12`, `D13`, `D23`), which
#    the review gives as estimates "from the kinetic gas theory" for the Duncan &
#    Toor (1962) system and prints with no temperature and no pressure. Nothing
#    was re-transcribed. What `A4.2` establishes about these rows: they are
#    kinetic-theory estimates, not measurements, and its recomputation of the
#    review's ideal-gas worked example from them flags none of the three. Its
#    findings concern the non-ideal example (`B12`'s printed sign), whose rows
#    this page does not use. This page makes no experimental claim on them.
#
# The review *does* contain a comparison with measurement, its Fig. 44: the
# dusty gas model against Remick & Geankoplis's (1974) He/Ne/Ar capillary data.
# That comparison is ternary, needs figure digitisation, and belongs to `A4.3`.
# The conditions printed inside that figure are recorded in the sidecar; the
# measured points are not, and this page makes no experimental claim.

# %%
pm = load_data("krishna-wesselingh-1997-porous-media.csv", page=PAGE)
pm_meta = load_meta("krishna-wesselingh-1997-porous-media.csv", page=PAGE)
P = {(r.group, r.symbol): r.value for r in pm.itertuples()}

wx = load_data("krishna-wesselingh-1997-worked-examples.csv", page=PAGE_A42)
wx_meta = load_meta("krishna-wesselingh-1997-worked-examples.csv", page=PAGE_A42)
W = {(r.example, r.symbol): r.value for r in wx.itertuples()}

print("page dataset :", cite_data(pm_meta))
print("cross-page   :", cite_data(wx_meta))

# The review prints the three pair diffusivities with no temperature and no
# pressure; 35.2 C and 1 atm are Duncan & Toor's conditions (page A4.9's
# sidecar). They set a mean molecular speed and a total concentration only.
T_REV = 308.35                       # K = 35.2 C  (Duncan & Toor's bath, via A4.9)
P_ATM = 101325.0                     # Pa          (Duncan & Toor's cell, via A4.9)
D_AB = {"H2-N2": W[("ideal", "D12")], "H2-CO2": W[("ideal", "D13")],
        "N2-CO2": W[("ideal", "D23")]}                       # m2/s, K&W p. 872
MOLAR = {"H2": 2.01588e-3, "N2": 28.0134e-3, "CO2": 44.0095e-3}   # kg/mol, IUPAC
TAU_CYL = P[("eq84", "tau_cylindrical")]     # eq. (84): tau = 1, cylindrical pore
EPS_POROSITY = 1.0                           # assumed, not printed: one straight open pore
EPS_TAU = EPS_POROSITY / TAU_CYL             # the eps/tau factor of eqs. (83), (85)
KN_DEN = P[("eq85", "knudsen_denominator")]                  # eq. (85): d_0/3
KN_EXP = P[("eq85", "knudsen_mass_exponent")]                # eq. (85): M^-1/2
GR_EXP = P[("eq107", "graham_mass_exponent")]                # eq. (107): sqrt(M)
NU_B = P[("jackson", "nu_b_example")]                        # nu_B = 2
DP_PRINTED = P[("jackson", "centre_pressure_increase")]      # "a 40% increase"
EQ85 = dict(eps_tau=EPS_TAU, den=KN_DEN, exp=KN_EXP)
MR = MOLAR["N2"] / MOLAR["H2"]                               # M_B/M_A for A = H2, B = N2
RHO_GRAHAM = (MOLAR["H2"] / MOLAR["N2"]) ** GR_EXP           # eq. (107) for the binary
for k, v in D_AB.items():
    print(f"  D_{k:8s} = {v:.3e} m2/s   (K&W p. 872)")
print(f"  eps/tau = {EPS_TAU:.1f}  (tau = {TAU_CYL:.0f} printed in eq. 84; "
      f"eps = {EPS_POROSITY:.0f} assumed)")

# %% [markdown]
# ## PyMRM implementation
#
# ### The Knudsen diffusivity, eq. (85)
#
# One line in `model.knudsen_D`, and the only place a molar mass enters. Its
# three scalings are checked against the exponents printed in eq. (85), and the
# pore size at which each pair's two resistances are equal ($\mathrm{Kn} = 1$)
# is printed: since eq. (85) gives $D_K = \bar v d_0/3$, it is also the
# elementary kinetic-theory mean free path $3\mathrm{D}_{AB}/\bar v$.

# %%
def eq85_study(mass_exp=KN_EXP, den=KN_DEN):
    """Ratios of D^e_iM for 4 T, 4 M and 4 d0 against sqrt(T), M^-1/2 and d0."""
    kw = dict(eps_tau=EPS_TAU, den=den, exp=mass_exp)
    base = knudsen_D(1e-7, T_REV, MOLAR["H2"], **kw)
    s = (knudsen_D(1e-7, 4 * T_REV, MOLAR["H2"], **kw) / base,
         knudsen_D(1e-7, T_REV, 4 * MOLAR["H2"], **kw) / base,
         knudsen_D(4e-7, T_REV, MOLAR["H2"], **kw) / base)
    return {"eq85_scaling_dev": max(abs(s[0] - 2.0), abs(s[1] - 0.5), abs(s[2] - 4.0)),
            "_s": s}


E85 = eq85_study()
D0_CROSS = {pair: KN_DEN * dab / np.sqrt(8 * RGAS * T_REV / (np.pi * MOLAR[pair.split("-")[0]]))
            for pair, dab in D_AB.items()}
for pair in D_AB:
    print(f"  Kn = 1 for {pair.split('-')[0]:3s} in {pair:8s}: d_0 = {D0_CROSS[pair]*1e9:5.1f} nm")
print("eq. (85) ratios for 4T, 4M, 4d0: " + ", ".join(f"{v:.12f}" for v in E85["_s"])
      + f"  (2, 0.5, 4); worst deviation {E85['eq85_scaling_dev']:.1e}")

d0s = np.logspace(-9, -4, 400)
fig, ax = plt.subplots(figsize=(7.0, 4.4))
for pair, dab in D_AB.items():
    dk = knudsen_D(d0s, T_REV, MOLAR[pair.split("-")[0]], **EQ85)
    ax.loglog(d0s * 1e9, 1.0 / (1.0 / dk + 1.0 / dab) / dab, lw=2,
              label=f"{pair.split('-')[0]} in {pair}")
    ax.plot([D0_CROSS[pair] * 1e9], [0.5], "o", ms=6, color="k", zorder=5)
ax.loglog(d0s * 1e9, np.minimum(knudsen_D(d0s, T_REV, MOLAR["H2"], **EQ85) / D_AB["H2-N2"], 1.0),
          "k--", lw=1, label="asymptotes (H$_2$ in H$_2$-N$_2$)")
ax.set_xlabel("pore diameter $d_0$ / nm")
ax.set_ylabel(r"$D_{\rm Bos}\,/\,\mathrm{D}^e_{12}$")
ax.set_title("Bosanquet interpolation, eq. (110), for K&W's p. 872 pair diffusivities\n"
             r"at 1 atm, 35.2 $^\circ$C (Duncan & Toor's conditions);  $\bullet$ = $\mathrm{Kn}=1$")
ax.legend(fontsize=8, loc="lower right"); ax.set_ylim(1e-3, 2)
plt.tight_layout(); plt.show()

# %% [markdown]
# ### Eq. (109) checked against eq. (82), symbolically
#
# Eq. (109) is the equation everything on this page rests on, it was read once
# off a scan, and its stoichiometric subscripts are what the OCR mangles. So it
# is derived independently from eq. (82) before being used.
#
# For a binary mixture with flux ratio $N_2/N_1 = -\rho$, eq. (82) for species 1
# reads $-\nabla c_1 = (x_2 N_1 - x_1 N_2)/\mathrm{D}^e_{12} + N_1/\mathrm{D}^e_{1M}$
# (ideal gas, so $\nabla p_1 = RT\nabla c_1$), and solving it for
# $D_1 = -N_1/\nabla c_1$ must reproduce eq. (109) with $\mathrm{D}_{\rm visc} = 0$.
# Setting $\rho = 1$ must then give eq. (110). Four residuals are formed: eq. (82)
# against eq. (109), eq. (109) at $\rho = 1$ against eq. (110), the closed form of
# the error derived below against the ratio of the two, and the root of the
# viscous factor against Graham's ratio. `nu_swap` reads eq. (109)'s ratio upside
# down ($x_1/\rho$ for $\rho x_1$), the misreading the scan invites, for the break
# table.

# %%
def symbolic_study(nu_swap=False):
    """Residuals of eq. (109) against eq. (82), eq. (110), the error form and Graham."""
    NA, xa, xb, dab, dam, dbm, rho, kn, ma, mb = sp.symbols(
        "N_A x_A x_B D_AB D_AM D_BM rho Kn M_A M_B", positive=True)
    row = (xb * NA - xa * (-rho * NA)) / dab + NA / dam      # eq. (82), species 1, N_2 = -rho N_1
    from82 = sp.simplify(sp.expand(row / NA))                # = 1/D_1
    eq109 = 1 / dam + (xb / dab) * (1 + (xa / (rho * xb) if nu_swap else rho * xa / xb))
    eq110 = sp.simplify(eq109.subs(rho, 1).subs(xb, 1 - xa))
    ratio = sp.simplify((eq109 / (1 / dam + 1 / dab)).subs(xb, 1 - xa)
                        .subs({dam: 1 + 1 / kn, dab: 1 + kn}))
    visc = sp.simplify(1 / dam - rho / dbm)                  # viscous factor of eq. (109)
    rho0 = sp.solve(sp.Eq(visc, 0), rho)[0]
    res = [sp.simplify(from82 - eq109),
           sp.simplify(eq110 - (1 / dam + 1 / dab)),
           sp.simplify(ratio - (1 + (rho - 1) * xa / (1 + kn))),
           sp.simplify(rho0.subs(dbm, dam * (mb / ma) ** sp.nsimplify(KN_EXP))
                       - (ma / mb) ** sp.nsimplify(GR_EXP))]
    return {"symbolic_nonzero_residuals": int(sum(r != 0 for r in res)), "_res": res,
            "_ratio": ratio, "_visc": visc, "_rho0": rho0, "_visc1": sp.simplify(visc.subs(rho, 1)),
            "_eq110": eq110, "_from82": from82}


SY = symbolic_study()
print("1/D_1 from eq. (82)       :", SY["_from82"])
print("residual against eq. (109):", SY["_res"][0])
print("eq. (109) at rho = 1      :", SY["_eq110"], "  residual against eq. (110):", SY["_res"][1])
print("D_Bos/D_1 (scaled)        :", SY["_ratio"], "  residual against 1 + (rho-1) x_A/(1+Kn):",
      SY["_res"][2])
print("viscous factor            :", SY["_visc"], "  vanishes at rho =", SY["_rho0"])
print("   with eq. (85), minus Graham's sqrt(M_A/M_B):", SY["_res"][3])
print("   its value at rho = 1   :", SY["_visc1"], "  (zero only if D_AM = D_BM, i.e. M_A = M_B)")
print(f"nonzero residuals: {SY['symbolic_nonzero_residuals']} of {len(SY['_res'])}")

# %% [markdown]
# **What that settles, and what it does not.** Eq. (109) as transcribed *is*
# eq. (82) rearranged, and it collapses to eq. (110) when $\rho = 1$, which is the
# review's condition (ii). It also gives a closed form for the whole question this
# page asks:
#
# $$
# \boxed{\;\frac{D_{\rm Bos}}{D_1} \;=\; 1 + \frac{(\rho - 1)\,x_1}{1 + \mathrm{Kn}}\;}
# $$
#
# The Bosanquet relation is exact when $\rho = 1$, for any Knudsen number; its
# error grows linearly in the flux-ratio departure $\rho - 1$ and in the mole
# fraction of the species itself; and it decays as $1/(1+\mathrm{Kn})$, so it is
# **largest in the bulk-diffusion limit and vanishes under Knudsen control**. That
# is a statement about the physics, not about the algebra: under Knudsen control
# species do not exchange momentum with each other, so what the other species is
# doing cannot matter, and any formula that ignores the flux ratio is right.
#
# The last symbolic lines identify the review's condition (iii). The viscous
# term of eq. (109) carries the factor
# $\sum_i(\nu_i/\nu_1)/\mathrm{D}^e_{iM} = 1/\mathrm{D}^e_{1M} -
# \rho/\mathrm{D}^e_{2M}$, and the $\rho$ that kills it is
# $\rho = \mathrm{D}^e_{2M}/\mathrm{D}^e_{1M} = \sqrt{M_1/M_2}$ by eq. (85), which
# is the flux ratio Graham's law imposes. **Condition (iii) is "the viscous term
# vanishes".** At $\rho = 1$ that factor is $1/\mathrm{D}^e_{1M} -
# 1/\mathrm{D}^e_{2M}$, so conditions (ii) and (iii) are two separate conditions
# and they close simultaneously only when $M_1 = M_2$.
#
# **This is algebra, and it is labelled as algebra.** It cannot detect an error
# in the pymrm assembly, in the boundary conditions, or in the discretisation.
# Everything below is there to test those.
#
# ### The solver
#
# One skeleton, three flux closures (`model.Cell`):
#
# - **`dgm`**: the dusty gas model as eq. (82) stands. The $2\times2$ matrix
#   $[B^e]$ of eq. (87) is built at every face from the local composition and
#   $(N) = -[B^e]^{-1}(\nabla c)$ is solved, which is eq. (86) for an ideal gas.
#   **Nothing about the flux ratio is assumed**: it comes out of the solve.
# - **`eq109`**: the scalar effective diffusivity of eq. (109), with the flux
#   ratio $\rho$ prescribed. Diagonal: each species responds to its own
#   gradient only.
# - **`bosanquet`**: eq. (110), the same with $\rho$ forced to 1, so the
#   diffusivities are constants.
#
# The three routes **differ only in the constitutive closure**: they share the
# face compositions, both operators, the residual assembly and the same `NumJac`,
# and only the body of `flux` changes. Their comparison therefore tests the
# presumed flux ratio; it cannot see a discretisation error, because a wrong grid
# moves both routes together.
#
# pymrm choices:
#
# - **Layout `(n, 2)`**: spatial axis first, species last, so `NumJac(shape,
#   axes_diagonals=[0], axes_blocks=[-1])` reads neighbouring cells through the
#   face fluxes and couples the two species pointwise through $[B^e]$.
# - **Boundary conditions on the outward normal.** At a pellet centre the outward
#   normal points inward and symmetry is `{a: 1, b: 0, d: 0}` ($dc/dn = 0$); a
#   fixed external composition is `{a: 0, b: 1, d: c_bnd}` ($c = c_{\rm bnd}$).
# - **Face compositions** from `construct_boundary_value_matrices` on both
#   boundary faces and the arithmetic mean of the neighbours inside; the same
#   matrices give the composition on the centre face, where the pellet-centre
#   pressure is read.
# - **`nu` in `construct_div` is geometry**: `0` slab, `1` cylinder, `2` sphere.
# - **Constant operators assembled once** in `Cell.__init__`; only the
#   composition-dependent flux is recomputed.
# - **`newton`** with bounds $c \ge 0$, a step tolerance of `TOL` relative to the
#   surface concentration and `stagnation=True`; `Cell.solve` raises unless it
#   reports success. The reaction is linear in $c_A$ and the problem is
#   nonlinear only through the composition dependence of $[B^e]$ and eq. (109).

# %%
p_demo, c_demo = pellet(1.0, NU_B, PHI, n=N_PELLET)
print(f"pellet, A -> {NU_B:.0f}B, Kn = 1, phi = {PHI:.0f}, {N_PELLET} cells: "
      f"p_0/p = {p_demo.centre_total(c_demo):.6f}, eta = {p_demo.eta(c_demo):.6f}")

# %% [markdown]
# ## Results
#
# ### 1. The isobaric counter-diffusion cell, where condition (ii) fails
#
# The textbook setting: a porous slab with fixed compositions on either side and
# **equal total pressure** at both ends. Everyone applies Bosanquet here. But
# $\nabla p = 0$ triggers eq. (107), so the flux ratio is not 1: it is
# $\sqrt{M_1/M_2}$, and Bosanquet's condition (ii) is violated by anything
# except two gases of the same molar mass.
#
# The `dgm` solver knows nothing about eq. (107). It is given equal total
# concentration at both ends and left to find the fluxes, so Graham's law is a
# *prediction* here. Because the flux ratio it predicts is constant, eq. (109)
# for species A can be integrated across the slab in closed form
# (`isobaric_flux_closed_form`: $N_A L = c_t\int_{x_2}^{x_1} D_A(x)\,dx$ with
# $1/D_A = 1/\mathrm{D}^e_{AM} + (1 + (\rho-1)x)/\mathrm{D}^e_{AB}$), and the
# Bosanquet flux is the same integral at $\rho = 1$. The quoted flux errors come
# from that closed form; the finite-volume `dgm` solve is checked against it.

# %%
def at_most(x, bound):
    """A round-off-level value printed as its bound (its digits vary with the BLAS kernel)."""
    return f"at most {bound:.0e}" if x <= bound else f"{x:.1e} (above {bound:.0e})"


def isobaric_study(mass_exp=KN_EXP, b_offdiag=1.0, n=N_SLAB):
    """H2 (A) against N2 (B), equal total concentration at both ends, over KN_GRID."""
    rho = MR ** mass_exp           # the flux ratio eq. (85) and grad p = 0 imply
    rows, gres, rdev, fv = {}, 0.0, 0.0, 0.0
    for kn in KN_GRID:
        s, c = isobaric_slab(kn, MR, xa=(XA1, XA2), n=n, mass_exp=mass_exp, b_offdiag=b_offdiag)
        N = s.flux(c)
        gl = N[:, 0] * MOLAR["H2"] ** GR_EXP + N[:, 1] * MOLAR["N2"] ** GR_EXP   # eq. (107)
        gres = max(gres, float(np.max(np.abs(gl)) / np.max(np.abs(N))))
        rdev = max(rdev, abs(float(-N[n // 2, 1] / N[n // 2, 0]) / RHO_GRAHAM - 1.0))
        dam, dab = d_groups(kn)
        n_dgm = isobaric_flux_closed_form(dam, dab, rho, XA1, XA2)
        n_bos = isobaric_flux_closed_form(dam, dab, 1.0, XA1, XA2)
        fv = max(fv, abs(float(N[n // 2, 0]) / n_dgm - 1.0))
        rows[kn] = (float(N[n // 2, 0]), n_dgm, n_bos, n_bos / n_dgm - 1.0)
    return {"graham_law_residual": gres, "graham_flux_ratio_dev": rdev,
            "isobaric_fv_vs_closed_form": fv,
            "bosanquet_flux_error_isobaric_bulk": rows[KN_GRID[0]][3],
            "bosanquet_flux_error_isobaric_knudsen": rows[KN_GRID[-1]][3],
            "isobaric_local_coefficient_error_bulk": (rho - 1.0) * XA1 / (1.0 + KN_GRID[0]),
            "_rows": rows}


ISO = isobaric_study()
print(f"Graham's law, eq. (107): N_B/N_A should be -sqrt(M_A/M_B) = -{RHO_GRAHAM:.6f} "
      f"for H2 (A) in N2 (B)")
print(f"  worst |sum N_i sqrt(M_i)| / max|N| over the solves : "
      f"{at_most(ISO['graham_law_residual'], 1e-10)}")
print(f"  worst deviation of the solved ratio from Graham    : "
      f"{at_most(ISO['graham_flux_ratio_dev'], 1e-10)}")
print(f"  finite-volume N_A ({N_SLAB} cells) against the closed form: "
      f"{ISO['isobaric_fv_vs_closed_form']:.1e} at worst\n")
print(f"{'Kn':>8} {'N_A (FV dgm)':>13} {'N_A (closed)':>13} {'N_A (Bos)':>11} {'Bosanquet error':>16}")
for kn, r in ISO["_rows"].items():
    print(f"{kn:8.3g} {r[0]:13.7f} {r[1]:13.7f} {r[2]:11.7f} {r[3]:15.2%}")

# %%
display(Markdown(rf"""
For hydrogen counter-diffusing against nitrogen at equal total pressure, the
Bosanquet formula underestimates the hydrogen flux by
**{-ISO['bosanquet_flux_error_isobaric_bulk']:.1%}** at $\mathrm{{Kn}} = {KN_GRID[0]:g}$ and by
**{-ISO['bosanquet_flux_error_isobaric_knudsen']:.2%}** at $\mathrm{{Kn}} = {KN_GRID[-1]:g}$. The
whole error is the flux ratio: Graham's law brings nitrogen back at only
{RHO_GRAHAM:.3f} times the hydrogen flux, so a net molar flow towards the
nitrogen side carries the hydrogen along, and eq. (110), which assumes the two
fluxes cancel, cannot see it.

The closed form predicts a *local* coefficient error of
$(\rho-1)x_A/(1+\mathrm{{Kn}})$, which at $x_A = {XA1}$ and
$\mathrm{{Kn}} = {KN_GRID[0]:g}$ is {ISO['isobaric_local_coefficient_error_bulk']:+.1%}. The flux
error is {ISO['bosanquet_flux_error_isobaric_bulk']:+.1%}, and the two differ because the
flux integrates $D_A(x_A)$ across the slab, where $x_A$ falls from {XA1} to {XA2}.
"""))

# %% [markdown]
# ### 2. The reacting pellet, where conditions (ii) and (iii) both fail
#
# A spherical pellet, pure A outside, first-order irreversible $A \to \nu_B B$
# on the pore walls. Stoichiometry fixes the flux ratio at $\rho = \nu_B$, so
# condition (ii) fails for any $\nu_B \neq 1$, and mass conservation fixes
# $M_B = M_A/\nu_B$: no free molar mass, and hence no way to satisfy (iii)
# either. A pressure gradient develops, and the review quotes Jackson's (1977)
# closed form for it: with complete conversion of A, the pressure at the pellet
# centre is $p_0 = \sqrt{\nu_B}\,p$.
#
# **That is the check with power on this page.** It is a printed result, it
# appears nowhere in the code, and it is not an interpolation: the solver is
# given eq. (82), eq. (85), the stoichiometry and the boundary conditions, and
# $\sqrt{\nu_B}$ has to come out. Jackson's result is the
# $\mathrm{Kn}\to\infty$ limit, so at a finite Knudsen number a departure from
# it remains; its scaling with Kn is measured on a fixed grid.

# %%
def jackson_study(mass_exp=KN_EXP, src_one=False, phi=PHI):
    """Centre pressure in the Knudsen limit for NU_LIST, and the finite-Kn residue."""
    kw = dict(n=N_PELLET, mass_exp=mass_exp, tol=TOL)
    rows = {}
    for nb in NU_LIST:
        p, c = pellet(KN_HI, nb, phi, nu_src=1.0 if src_one else None, **kw)
        p0 = p.centre_total(c)
        rows[nb] = (p0, abs(p0 / np.sqrt(nb) - 1.0), float(p.faces(c)[0][0, 0]))
    tail = []
    for k in KN_TAIL:
        p, c = pellet(k, NU_B, phi, nu_src=1.0 if src_one else None, **kw)
        dev = abs(p.centre_total(c) / np.sqrt(NU_B) - 1.0)
        tail.append((k, dev, dev * k))
    prods = [t[2] for t in tail]
    return {"jackson_centre_pressure_dev": max(r[1] for r in rows.values()),
            "jackson_nu2_solved_increase": rows[NU_B][0] - 1.0,
            "jackson_dev_x_kn_spread_1e5_to_1e8": max(prods) / min(prods) - 1.0,
            "_rows": rows, "_tail": tail}


JK = jackson_study()
print("Jackson's pellet-centre pressure, quoted by K&W p. 892, at Kn = "
      f"{KN_HI:.0e}, phi = {PHI:.0f}, {N_PELLET} cells")
print(f"{'nu_B':>6} {'p_0/p solved':>14} {'sqrt(nu_B)':>12} {'rel dev':>10} {'c_A(0)':>10}")
for nb, v in JK["_rows"].items():
    dev = f"{v[1]:10.2e}" if v[1] > 1e-13 else f"{'< 1e-13':>10}"   # nu_B = 1: round-off
    print(f"{nb:6.1f} {v[0]:14.7f} {np.sqrt(nb):12.7f} {dev} {v[2]:10.2e}")
print(f"\nthe departure against Kn (nu_B = {NU_B:.0f}, same grid):")
print(f"{'Kn':>10} {'rel dev':>12} {'dev x Kn':>12}")
for k, dev, prod in JK["_tail"]:
    print(f"{k:10.0e} {dev:12.3e} {prod:12.5f}")
print(f"dev x Kn constant to {JK['jackson_dev_x_kn_spread_1e5_to_1e8']:.2%} over the four values of Kn")

# %%
display(Markdown(rf"""
The review's sentence is *"Thus, for $\nu_B$ = 2, we have a
{DP_PRINTED:.0%} increase in pressure as we proceed towards the centre of the
pellet"*. The pymrm dusty-gas solve gives
**{JK['jackson_nu2_solved_increase']:.2%}**, and the review's own closed form
$p_0 = \sqrt{{\nu_B}}\,p$ gives $\sqrt{{2}} - 1$ = {np.sqrt(NU_B) - 1:.2%}. The
printed {DP_PRINTED:.0%} is the authors' rounding of their own expression; the
solve reproduces the expression for $\nu_B$ = 1, 2, 3 and 4 to
{JK['jackson_centre_pressure_dev']:.1e} relative at worst.

**What that last figure is.** Jackson's $p_0 = \sqrt{{\nu_B}}\,p$ is the
$\mathrm{{Kn}}\to\infty$ limit, and the solve runs at
$\mathrm{{Kn}} = 10^{{{np.log10(KN_HI):.0f}}}$, so a bulk-friction remnant of eq. (82)
is still present. The table above measures it: the departure times
$\mathrm{{Kn}}$ is constant to {JK['jackson_dev_x_kn_spread_1e5_to_1e8']:.2%} across
four decades of $\mathrm{{Kn}}$ on a fixed grid, so it falls as $1/\mathrm{{Kn}}$, the
finite-Kn term of the model. The centre pressure is read on the $r = 0$ face,
from the same boundary-value matrices the solver uses.
"""))

# %% [markdown]
# **Is that agreement worth anything?** Validation section V3 measures it
# directly. Two properties are worth noting first, and both are checked below:
# $p_0/\sqrt{\nu_B}\,p$ must be independent of pellet geometry, of the Thiele
# modulus and of the grid (Jackson's result contains none of them), and it is
# *not* a discriminator between the dusty gas model and Bosanquet, because in
# the Knudsen limit the two coincide.
#
# What it does test is narrower than it looks: the $M^{-1/2}$ of eq. (85) and
# the stoichiometric closure. It does **not** test the signs of the friction
# terms of eq. (82): V3 flips the off-diagonals of $[B^e]$ and the centre
# pressure moves by about as much as the finite-Kn residue this check already
# carries. There is no ideal-gas step to test either, because the solver works in
# concentrations from end to end and never forms a partial pressure. The checks
# that do catch a defect in the friction terms are V1 (the two closures against
# each other) and Graham's law; the break table runs both.
#
# ### 3. The Bosanquet error across the whole transition
#
# Sweep the Knudsen number with both closures and read off what the
# approximation costs in the two observables an engineer uses: the internal
# pressure and the effectiveness factor. The sweep runs on `N_PELLET` cells; the
# numbers quoted after it come from a refinement ladder at each quoted condition,
# extrapolated to zero cell size.

# %%
def sweep_study(nu_b=NU_B, nu_src=None):
    """Both closures over KN_SWEEP at PHI on N_PELLET cells."""
    rows = {}
    for kn in KN_SWEEP:
        r = {}
        for model in ("dgm", "bosanquet"):
            p, c = pellet(kn, nu_b, PHI, n=N_PELLET, model=model, nu_src=nu_src, tol=TOL)
            r[model] = (p.centre_total(c), p.eta(c))
        # constant D: D_B c_B + nu_B D_A c_A is harmonic with zero flux at the
        # centre, so the Bosanquet model's c_t(0) is nu_B D_A/D_B once c_A(0) -> 0
        r["closed"] = nu_b * (1.0 / (1.0 / p.dam + 1.0 / p.dab)) / (1.0 / (1.0 / p.dbm + 1.0 / p.dab))
        r["err"] = r["bosanquet"][1] / r["dgm"][1] - 1.0
        rows[kn] = r
    err = np.array([rows[k]["err"] for k in KN_SWEEP])
    drops = -np.diff(err)
    return {"bosanquet_p0_vs_closed_form": max(abs(r["bosanquet"][0] / r["closed"] - 1.0)
                                               for r in rows.values()),
            "bosanquet_eta_error_decreasing_in_kn": int(np.sum(drops > 0)),
            "_drops": drops, "_rows": rows}


SW = sweep_study()
kn = KN_SWEEP
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
axes[0].semilogx(kn, [SW["_rows"][k]["dgm"][0] for k in kn], "o-", lw=2,
                 label="dusty gas model, eqs. (82)/(87)")
axes[0].semilogx(kn, [SW["_rows"][k]["bosanquet"][0] for k in kn], "s--", lw=2,
                 label="Bosanquet, eq. (110)")
axes[0].axhline(np.sqrt(NU_B), color="k", ls=":", label=r"Jackson $\sqrt{\nu_B}$ (K&W p. 892)")
axes[0].axhline(1.0, color="0.6", lw=1)
axes[0].set_xlabel(r"$\mathrm{Kn} = \mathrm{D}^e_{AB}/\mathrm{D}^e_{AM}$")
axes[0].set_ylabel(r"$p_0/p$  (pellet centre)")
axes[0].set_title(f"internal pressure, $A\\to{NU_B:.0f}B$, $\\phi={PHI:.0f}$, {N_PELLET} cells")
axes[0].legend(fontsize=8)
axes[1].semilogx(kn, [SW["_rows"][k]["err"] * 100 for k in kn], "o-", lw=2, color="C3")
axes[1].axhline(0.0, color="0.6", lw=1)
axes[1].set_xlabel(r"$\mathrm{Kn} = \mathrm{D}^e_{AB}/\mathrm{D}^e_{AM}$")
axes[1].set_ylabel(r"$\eta_{\rm Bos}/\eta_{\rm DGM} - 1$  /  %")
axes[1].set_title("error in the effectiveness factor")
plt.tight_layout(); plt.show()

print(f"{'Kn':>9} {'p0/p DGM':>10} {'p0/p Bos':>10} {'eta DGM':>9} {'eta Bos':>9} {'eta error':>10}")
for k in KN_SWEEP[::3]:
    r = SW["_rows"][k]
    print(f"{k:9.3g} {r['dgm'][0]:10.5f} {r['bosanquet'][0]:10.5f} "
          f"{r['dgm'][1]:9.5f} {r['bosanquet'][1]:9.5f} {r['err']:9.2%}")
print(f"\nthe error falls between consecutive Kn at {SW['bosanquet_eta_error_decreasing_in_kn']} "
      f"of {len(KN_SWEEP) - 1} steps (smallest fall {SW['_drops'].min():.1e})")
print(f"Bosanquet centre pressure against its closed form nu_B D_A/D_B: "
      f"{SW['bosanquet_p0_vs_closed_form']:.1e} at worst")

# %%
HEAD_CASES = {"bosanquet_eta_error_bulk": (1e-3, 30.0),
              "bosanquet_eta_error_bulk_phi10": (1e-3, 10.0),
              "bosanquet_eta_error_bulk_phi100": (1e-3, 100.0),
              "bosanquet_eta_error_kn1": (1.0, 30.0),
              "bosanquet_eta_error_knudsen": (KN_HI, 30.0)}


def headline_study(nu_b=NU_B, mass_exp=KN_EXP, face_rule="mean", ns=NS_HEAD):
    """Effectiveness-factor error and centre pressures on a refinement ladder per case."""
    out, lad, rel = {}, {}, []
    for name, (kn_, phi_) in HEAD_CASES.items():
        def q(n, kn_=kn_, phi_=phi_):
            kw = dict(n=n, mass_exp=mass_exp, face_rule=face_rule, tol=TOL)
            a, ca = pellet(kn_, nu_b, phi_, **kw)
            b, cb = pellet(kn_, nu_b, phi_, model="bosanquet", **kw)
            return np.array([b.eta(cb) / a.eta(ca) - 1.0, a.centre_total(ca), b.centre_total(cb)])
        st = observed_orders(q, ns)
        vals, ext = np.array(st["values"]), np.array(st["extrapolated"])
        out[name] = float(ext[0])
        rel.append(abs(vals[-1][0] - ext[0]) / abs(ext[0]))
        d = np.abs(np.diff(vals[:, 0]))
        lad[name] = (vals, ext, np.log(d[:-1] / d[1:]) / np.log(ns[1] / ns[0]))
    bulk = lad["bosanquet_eta_error_bulk"]
    out["bosanquet_eta_error_bulk_grid_error"] = abs(float(bulk[0][list(ns).index(N_PELLET), 0]
                                                          - bulk[1][0])) if N_PELLET in ns else np.nan
    out["dgm_centre_pressure_bulk"] = float(bulk[1][1])
    out["bosanquet_centre_pressure_bulk"] = float(bulk[1][2])
    out["eta_extrapolation_worst"] = max(rel)
    phis = [out["bosanquet_eta_error_bulk_phi10"], out["bosanquet_eta_error_bulk"],
            out["bosanquet_eta_error_bulk_phi100"]]
    out["bosanquet_eta_error_grows_with_phi"] = int(np.sum(np.diff(phis) > 0))
    out["_rises"], out["_lad"] = np.diff(phis), lad
    return out


HD = headline_study()
print(f"{'case':<34} {'Kn':>7} {'phi':>5} {'600 cells':>11} {'extrapolated':>13}  orders on {NS_HEAD}")
for name, (kn_, phi_) in HEAD_CASES.items():
    vals, ext, ords = HD["_lad"][name]
    print(f"{name:<34} {kn_:7.0e} {phi_:5.0f} {vals[1][0]:11.4e} {ext[0]:13.6e}  "
          + " ".join(f"{o:.2f}" for o in ords))
print(f"\nfinest level against the extrapolated value, worst over the cases: "
      f"{HD['eta_extrapolation_worst']:.1e} (relative)")
print(f"the sweep's {N_PELLET} cells miss the bulk-limit value by "
      f"{HD['bosanquet_eta_error_bulk_grid_error']:.2e}")
print(f"centre pressure at Kn = 1e-3, extrapolated: DGM {HD['dgm_centre_pressure_bulk']:.6f}, "
      f"Bosanquet {HD['bosanquet_centre_pressure_bulk']:.6f}")

# %%
display(Markdown(rf"""
**Read the left panel from the right.** Deep in the Knudsen regime the two
closures sit on top of each other, both on Jackson's $\sqrt{{\nu_B}}$ line.
Coming left, the dusty gas model relaxes towards $p_0/p = 1$ (with bulk
diffusion available, the extra moles of B leave without building a pressure)
while the Bosanquet model runs the other way, to
{HD['bosanquet_centre_pressure_bulk']:.3f}, because constant diffusivities force the
total concentration at the centre to $\nu_B D_A/D_B$ regardless of how the
species couple. At $\mathrm{{Kn}} = 10^{{-3}}$ the dusty gas model puts the internal
pressure rise at {HD['dgm_centre_pressure_bulk'] - 1:.2%} and Bosanquet at
{HD['bosanquet_centre_pressure_bulk'] - 1:.0%}.

**The effectiveness factor** is the number that reaches a reactor model. For
this pellet ($A\to{NU_B:.0f}B$, $\phi = {PHI:.0f}$, extrapolated to zero cell size)
Bosanquet overestimates it by **{HD['bosanquet_eta_error_bulk']:.1%}** at
$\mathrm{{Kn}} = 10^{{-3}}$, {HD['bosanquet_eta_error_kn1']:.1%} at $\mathrm{{Kn}} = 1$, and
{HD['bosanquet_eta_error_knudsen']:.1e} (relative) at $\mathrm{{Kn}} = 10^{{6}}$.

**Those percentages are conditions-specific, and the condition is the Thiele
modulus.** At $\mathrm{{Kn}} = 10^{{-3}}$ the same error is
{HD['bosanquet_eta_error_bulk_phi10']:.1%} at $\phi = 10$ and
{HD['bosanquet_eta_error_bulk_phi100']:.1%} at $\phi = 100$: it grows with $\phi$
across the three values, because a thinner reacting shell samples the
composition where $x_A$ is closest to 1 and the closed form
$(\rho-1)x_A/(1+\mathrm{{Kn}})$ is largest. What is not conditions-specific is the
**shape**: on the sweep at $\phi = {PHI:.0f}$ the error falls at every step in
$\mathrm{{Kn}}$ and dies under Knudsen control.

So the interpolation whose purpose is to add the Knudsen resistance is accurate
when the Knudsen resistance dominates and inaccurate when it does not, and the
practical reading is the reverse of the usual instinct: *if your pores are small
enough for Knudsen diffusion to matter, Bosanquet is fine; if they are not, do
not reach for it just because the two terms are still comparable.*
"""))

# %% [markdown]
# ### 4. The closed form, against the solves
#
# The boxed expression $D_{\rm Bos}/D_1 = 1 + (\rho-1)x_A/(1+\mathrm{Kn})$ was
# derived above with sympy. Here it is evaluated against the diffusivities the
# solver used, at every face of every converged pellet of the sweep. This is an
# **algebraic identity** (both sides are built from `d_scalar`), so it confirms
# the derivation and the implementation of eq. (109), and nothing about the
# physics. It is shown because it is the map.

# %%
def closed_map_study(nu_swap=False):
    """Closed-form error map against the coefficients of the converged pellets."""
    dev = 0.0
    for kn_ in KN_SWEEP:
        p, c = pellet(kn_, NU_B, PHI, n=N_PELLET, nu_swap=nu_swap, tol=TOL)
        cf, _ = p.faces(c)
        xf = cf / cf.sum(axis=1, keepdims=True)
        got = p.d_scalar(xf, True)[:, 0] / p.d_scalar(xf, False)[:, 0]
        pred = 1.0 + (p.rho - 1.0) * xf[:, 0] / (1.0 + kn_)
        dev = max(dev, float(np.max(np.abs(got - pred))))
    return {"closed_form_error_map_dev": dev}


MAP = closed_map_study()
fig, ax = plt.subplots(figsize=(7.0, 4.2))
knc = np.logspace(-3, 3, 300)
for rho, lab in [(RHO_GRAHAM, rf"$\rho={RHO_GRAHAM:.3f}$  (H$_2$/N$_2$ isobaric, Graham)"),
                 (1.0, r"$\rho=1$  (Bosanquet exact)"),
                 (2.0, r"$\rho=2$  ($A\to2B$)"),
                 (4.0, r"$\rho=4$  ($A\to4B$)")]:
    ax.semilogx(knc, (rho - 1.0) / (1.0 + knc) * 100, lw=2, label=lab)
ax.axhline(0, color="0.6", lw=1); ax.axvline(1.0, color="0.6", ls=":", lw=1)
ax.set_xlabel(r"$\mathrm{Kn} = \mathrm{D}^e_{AB}/\mathrm{D}^e_{AM}$")
ax.set_ylabel(r"$D_{\rm Bos}/D_1 - 1$  at $x_A = 1$  /  %")
ax.set_title("How wrong the Bosanquet relation is, and where\n"
             r"(set by the flux ratio, decaying as $1/(1+\mathrm{Kn})$)")
ax.legend(fontsize=8); plt.tight_layout(); plt.show()
print(f"closed form against the coefficients the solver used, all faces, all Kn: "
      f"{MAP['closed_form_error_map_dev']:.1e}  (algebraic identity)")

# %% [markdown]
# ### 5. The review's two scaling claims
#
# p. 887 states them plainly: the pair diffusivities are inversely proportional
# to pressure and independent of pore size, so under bulk control the fluxes are
# independent of both; the Knudsen diffusivities are independent of pressure and
# proportional to pore size, so under Knudsen control the fluxes are
# proportional to both. What is measured is the **local** slope
# $\partial\ln N/\partial\ln p$ and $\partial\ln N/\partial\ln d_0$ of the
# dimensional isobaric H₂/N₂ slab (`N_DIM` cells, the review's $\mathrm{D}_{12}$
# at 1 atm scaled as $1/p$), by central differences in $\ln p$ and $\ln d_0$, at
# two named points: Knudsen control at $d_0 = 20$ nm, $p = 1$ kPa, and bulk
# control at $d_0 = 20$ µm, $p = 1$ MPa. Each departure from the stated exponent
# is a finite-Kn effect, so it is measured again a decade deeper into its limit
# (100 Pa and 10 MPa) and the rate $d\ln(\text{departure})/d\ln\mathrm{Kn}$
# reported: $-1$ under Knudsen control, $+1$ under bulk control, if the
# departures are the leading finite-Kn terms. The same slopes from
# `isobaric_flux_closed_form` check the solver.

# %%
SCAL_POINTS = {"knudsen": (20e-9, 1e3, 1e2), "bulk": (20e-6, 1e6, 1e7)}   # d0 (m), p, p deeper (Pa)
SCAL_TARGET = {"knudsen": 1.0, "bulk": 0.0}                                # p. 887 exponents


def scaling_study(T=None, d_scale=1.0, pair_const=False, dk_const=False, b_offdiag=1.0):
    """Local slopes of the slab flux in ln p and ln d0 at the two named points."""
    T = T_REV if T is None else T

    def groups(p, d0, d0_named):
        dab = d_scale * D_AB["H2-N2"] * (1.0 if pair_const else P_ATM / p)   # D_ij ~ 1/p
        dam = knudsen_D(d0_named if dk_const else d0, T, MOLAR["H2"], **EQ85)   # eq. (85)
        return dam, dab

    def ln_flux(p, d0, d0_named, closed):
        dam, dab = groups(p, d0, d0_named)
        ct = p / (RGAS * T)
        if closed:
            return np.log(isobaric_flux_closed_form(dam, dab, RHO_GRAHAM, XA1, XA2, ct=ct, L=DELTA))
        s, c = isobaric_slab(None, MR, xa=(XA1, XA2), n=N_DIM, dam=dam, dab=dab, ct=ct,
                             L=DELTA, b_offdiag=b_offdiag)
        return np.log(s.flux(c)[N_DIM // 2, 0])

    def slopes(d0, p, closed=False):
        e = np.exp(H_LOG)
        return np.array([(ln_flux(p * e, d0, d0, closed) - ln_flux(p / e, d0, d0, closed)) / (2 * H_LOG),
                         (ln_flux(p, d0 * e, d0, closed) - ln_flux(p, d0 / e, d0, closed)) / (2 * H_LOG)])

    out, gap, tab = {}, 0.0, {}
    for regime, (d0, p, p_deep) in SCAL_POINTS.items():
        s, s_deep = slopes(d0, p), slopes(d0, p_deep)
        dev = float(np.max(np.abs(s - SCAL_TARGET[regime])))
        dev_deep = float(np.max(np.abs(s_deep - SCAL_TARGET[regime])))
        kn_, kn_deep = (lambda g: g[1] / g[0])(groups(p, d0, d0)), (lambda g: g[1] / g[0])(groups(p_deep, d0, d0))
        out[f"scaling_slope_dev_{regime}"] = dev
        with np.errstate(divide="ignore", invalid="ignore"):   # nan when a defect freezes Kn
            out[f"scaling_rate_{regime}"] = float(np.log(dev_deep / dev) / np.log(kn_deep / kn_))
        gap = max(gap, float(np.max(np.abs(s - slopes(d0, p, closed=True)))))
        tab[regime] = (kn_, s, dev, kn_deep, dev_deep)
    out["scaling_slope_vs_closed_form"] = gap
    out["_tab"] = tab
    return out


SC = scaling_study()
print(f"{'regime':>8} {'Kn':>10} {'dlnN/dlnp':>11} {'dlnN/dlnd0':>11} {'departure':>10} "
      f"{'Kn deeper':>10} {'departure':>10}")
for regime, (kn_, s, dev, kn_deep, dev_deep) in SC["_tab"].items():
    print(f"{regime:>8} {kn_:10.3g} {s[0]:11.6f} {s[1]:11.6f} {dev:10.2e} {kn_deep:10.3g} {dev_deep:10.2e}")
print(f"\nrate d ln(departure)/d ln Kn: Knudsen {SC['scaling_rate_knudsen']:.3f} (expected -1), "
      f"bulk {SC['scaling_rate_bulk']:.3f} (expected +1)")
print(f"finite-volume slopes against the closed-form slopes: {SC['scaling_slope_vs_closed_form']:.1e}")

# %% [markdown]
# ### 5b. What the borrowed $T$, $p$ and diffusivities are worth
#
# The 35.2 °C, 1 atm and the three pair diffusivities are all inherited: the
# review prints no temperature and no pressure, and the diffusivities come
# cross-page from `A4.2`. Everything else on the page is computed in the
# dimensionless groups $(\mathrm{Kn}, \rho, \phi)$, which contain none of the three.
# The scaling check is the exception: its two points are named by a pore size and
# a pressure, and the Knudsen number there is set by exactly these numbers. So it
# is re-run at $T = 400$ K and with the borrowed $\mathrm{D}_{ij}$ scaled by 1.3.

# %%
T_BREAK, D_BREAK = 400.0, 1.3


def perturbed_study(pair_const=False):
    """The two slope departures at T = 400 K and with the pair diffusivities x 1.3."""
    a = scaling_study(T=T_BREAK, pair_const=pair_const)
    b = scaling_study(d_scale=D_BREAK, pair_const=pair_const)
    return {"scaling_slope_dev_knudsen_T400": a["scaling_slope_dev_knudsen"],
            "scaling_slope_dev_bulk_T400": a["scaling_slope_dev_bulk"],
            "scaling_slope_dev_knudsen_D1p3": b["scaling_slope_dev_knudsen"],
            "scaling_slope_dev_bulk_D1p3": b["scaling_slope_dev_bulk"]}


PT = perturbed_study()
print(f"{'substitution':>32}{'Knudsen departure':>19}{'bulk departure':>16}")
for lab, k, b in [(f"none (T = {T_REV:.2f} K, D as loaded)", SC["scaling_slope_dev_knudsen"],
                   SC["scaling_slope_dev_bulk"]),
                  (f"T -> {T_BREAK:.0f} K", PT["scaling_slope_dev_knudsen_T400"],
                   PT["scaling_slope_dev_bulk_T400"]),
                  (f"borrowed D x {D_BREAK}", PT["scaling_slope_dev_knudsen_D1p3"],
                   PT["scaling_slope_dev_bulk_D1p3"])]:
    print(f"{lab:>32}{k:>19.2e}{b:>16.2e}")

# %% [markdown]
# ## Validation
#
# Ranked, and **not of equal weight**. There is no measurement anywhere on this
# page: the review prints one dusty-gas/experiment comparison (Fig. 44) and it is
# ternary and belongs to `A4.3`. **This page is tier 6**: printed constants,
# printed closed forms, and internal identities.
#
# 1. **A stated result, reproduced.** Jackson's pellet-centre pressure
#    $p_0 = \sqrt{\nu_B}\,p$, quoted by the review on p. 892 with its "40 %
#    increase" for $\nu_B = 2$ (Results 2). The solver encodes it nowhere. V3
#    and V4 measure what that agreement is and is not sensitive to.
# 2. **A second stated result, reproduced.** Graham's law, eq. (107), emerges
#    from the isobaric slab, which is given equal total concentration at both
#    ends and nothing else (Results 1). It is an identity of eqs. (82)+(85) once
#    $\nabla p = 0$ is imposed, so it is not independent evidence for the model:
#    it checks the $M^{-1/2}$ exponent of eq. (85) and the friction terms of
#    eq. (82), and the break table shows that a defect in either moves it.
# 3. **The review's two scaling claims**, as local slopes with their rate of
#    approach (Results 5).
# 4. **The transcription of eq. (109)**, re-derived symbolically from eq. (82)
#    and collapsing to eq. (110) at $\rho = 1$. Algebra, labelled as algebra.
# 5. **Two closures that differ only in the constitutive relation** (V1).
# 6. **Grid convergence** (V2) and the quoted numbers extrapolated on their own
#    ladders (Results 3).
# 7. **The geometry index**, against the classical effectiveness factor of a
#    slab, a cylinder and a sphere (V5).
# 8. **Closed forms the solves must meet**: the isobaric flux (Results 1), the
#    Bosanquet model's centre pressure $\nu_B D_A/D_B$ (Results 3), and the
#    closed-form error map against the coefficients the solver used (Results 4,
#    algebraic).
#
# ### V1. The two closures against each other
#
# Where the flux ratio the scalar closure presumes is the one the dusty gas
# model produces, the two must agree. In the pellet the discrete divergence
# telescopes, so $N_B = -\nu_B N_A$ holds on every face of the `dgm` solve and
# eq. (109) with $\rho = \nu_B$ gives the same face fluxes: the agreement is an
# identity of the scheme, at the solve's tolerance. Where the presumed ratio is
# wrong (eq. 110's $\rho = 1$ against $\nu_B = 2$), they must not agree. Both
# directions are checked, because only the second one has any power; the break
# table runs the second at $\nu_B = 1$, where eq. (110) is exact.

# %%
def closure_study(b_offdiag=1.0, nu_swap=False, nu_wrong=NU_B):
    """dgm against eq. (109) with the stoichiometric ratio, and against eq. (110)."""
    agree = 0.0
    for kn_ in [1e-2, 1e0, 1e2]:
        for nb in [1.0, 2.0, 4.0]:
            _, ca = pellet(kn_, nb, 20.0, n=N_PELLET, b_offdiag=b_offdiag, tol=TOL)
            _, cb = pellet(kn_, nb, 20.0, n=N_PELLET, model="eq109", nu_swap=nu_swap, tol=TOL)
            agree = max(agree, float(np.max(np.abs(ca - cb)) / np.max(np.abs(ca))))
    _, cw = pellet(1e-2, nu_wrong, 20.0, n=N_PELLET, model="bosanquet", tol=TOL)
    _, ct = pellet(1e-2, nu_wrong, 20.0, n=N_PELLET, tol=TOL)
    return {"dgm_vs_eq109_correct_rho": agree,
            "dgm_vs_eq110_wrong_rho": float(np.max(np.abs(cw - ct)) / np.max(np.abs(ct)))}


V1 = closure_study()
print(f"dgm against eq. (109) with the stoichiometric flux ratio (9 pellets): "
      f"{at_most(V1['dgm_vs_eq109_correct_rho'], 1e-12)}")
print(f"dgm against eq. (110), flux ratio forced to 1, Kn = 0.01, nu_B = {NU_B:.0f}: "
      f"{V1['dgm_vs_eq110_wrong_rho']:.2e}")

# %% [markdown]
# ### V2. Grid convergence, with the round-off companion
#
# The order of the scheme on the effectiveness factor at $\mathrm{Kn} = 1$,
# $\phi = 10$, $A \to 2B$, from successive differences on a ladder of eight grids.
# The orders fall towards 2 along the ladder (the higher-order terms of the thin
# reacting shell fade as the cells shrink), so the claim is the order of the
# finest triple, with the whole ladder printed. The finest difference is compared
# with the Newton tolerance, so that the order is not read off round-off. The
# numbers quoted in Results 3 are extrapolated on their own ladders at their own
# conditions.

# %%
def grid_study(ns=NS_ORDER, face_rule="mean", tol=TOL):
    """Observed orders of eta on a refinement ladder, and the finest difference over tol."""
    def eta_n(n):
        p, c = pellet(1.0, NU_B, 10.0, n=n, face_rule=face_rule, tol=tol)
        return p.eta(c)
    st = observed_orders(eta_n, ns)
    v = np.asarray(st["values"])
    return {"grid_convergence_order": float(st["orders"][-1]),
            "grid_ladder_over_roundoff": float(abs(v[-1] - v[-2]) / (tol * abs(v[-1]))),
            "_st": st}


V2 = grid_study()
st = V2["_st"]
print(f"{'n':>6} {'eta':>14} {'|diff|':>12} {'order':>7}")
for i, n in enumerate(NS_ORDER):
    d = f"{abs(st['values'][i] - st['values'][i - 1]):12.3e}" if i else " " * 12
    o = f"{st['orders'][i - 2]:7.3f}" if i >= 2 else ""
    print(f"{n:6d} {st['values'][i]:14.10f} {d} {o}")
print(f"order of the finest triple: {V2['grid_convergence_order']:.3f}; finest difference over "
      f"the Newton tolerance: {V2['grid_ladder_over_roundoff']:.1e}")

# %% [markdown]
# ### V3. What the $\sqrt{\nu_B}$ agreement is sensitive to
#
# An agreement number that is algebraically guaranteed and presented as evidence
# is the most common defect of a page like this. So the headline check is broken
# on purpose, and the result decides what the page may claim.
#
# **What it catches.** A flipped mass exponent in eq. (85), a dropped mass
# dependence, a wrong stoichiometric coefficient in the source term, and an
# incomplete conversion ($\phi = 3$): each is a row of the break table below, and
# each moves the Jackson claims outside their tolerance.
#
# **What it cannot see**, measured here in the Knudsen limit where the check is
# made: a sign flip of the off-diagonals of $[B^e]$ (the $-x_i/\mathrm{D}^e_{ij}$
# of eq. 87, where species–species friction enters), a slab in place of the
# sphere, and a five-cell grid. Each moves the centre pressure by the shift
# printed below. The blindness is a property of the Knudsen limit, not of the
# defects: the break table re-measures the same three at $\mathrm{Kn} = 1$ and at
# $\phi = 3$, where they show. So the $\sqrt{\nu_B}$ agreement says **nothing**
# about the friction terms, the geometry or the grid; V1, Graham's law, V2 and V5
# do.

# %%
def blindness_study(kn=KN_HI, phi=PHI):
    """Shift of the centre pressure under three defects the Knudsen limit hides."""
    kw = dict(tol=TOL)
    p, c = pellet(kn, NU_B, phi, n=N_PELLET, **kw)
    ref = p.centre_total(c)
    out = {}
    for name, extra in [("jackson_shift_b_offdiag_flipped", dict(b_offdiag=-1.0, n=N_PELLET)),
                        ("jackson_shift_slab", dict(geom=0, n=N_PELLET)),
                        ("jackson_shift_five_cells", dict(n=5))]:
        q, d = pellet(kn, NU_B, phi, **extra, **kw)
        out[name] = abs(q.centre_total(d) / ref - 1.0)
    return out


V3 = blindness_study()
print(f"centre pressure at Kn = {KN_HI:.0e}, phi = {PHI:.0f}; relative shift under each defect:")
for k, v in V3.items():
    print(f"  {k:<34} {v:.1e}")
print(f"the finite-Kn departure from sqrt(nu_B) itself: {JK['_rows'][NU_B][1]:.1e}")

# %% [markdown]
# ### V4. Invariance of the headline result
#
# Jackson's result contains no geometry, no Thiele modulus and no grid, so the
# centre pressure in the Knudsen limit must be the same for every combination of
# three geometries, three Thiele moduli (all with complete conversion) and two
# grids.

# %%
def invariance_study(phis=(20.0, 40.0, 80.0)):
    """Spread of p_0/p over geometry x phi x grid at Kn = KN_HI."""
    vals = [pellet(KN_HI, NU_B, ph, n=n, geom=g, tol=TOL) for g in (0, 1, 2)
            for ph in phis for n in (100, 400)]
    v = np.array([p.centre_total(c) for p, c in vals])
    return {"jackson_invariance_spread": float(v.max() - v.min()), "_v": v}


V4 = invariance_study()
print(f"p_0/p over {len(V4['_v'])} combinations: min {V4['_v'].min():.8f}, max {V4['_v'].max():.8f}, "
      f"spread {V4['jackson_invariance_spread']:.1e}")

# %% [markdown]
# ### V5. The geometry index, against a closed form
#
# V3 shows that the $\sqrt{\nu_B}$ check is blind to a wrong `nu` in
# `construct_div`, and V4 sweeps geometry as an invariance. So the geometry is
# checked directly. Set $\nu_B = 1$: then $\rho = 1$, eq. (109) collapses to
# $1/D_A = 1/\mathrm{D}^e_{AM} + 1/\mathrm{D}^e_{AB}$, a **constant**, and the
# dusty-gas route gives the same thing, because the friction term
# $(x_A + x_B)/\mathrm{D}^e_{AB}$ loses its composition dependence. The pellet is
# then the textbook first-order problem, whose effectiveness factor is known in
# closed form and is *different for each geometry*: $\tanh\phi/\phi$ (slab),
# $2I_1(\phi)/\phi I_0(\phi)$ (cylinder), $(3/\phi)(\coth\phi - 1/\phi)$ (sphere).
# None of the three appears in the solver. The break table compares each solve
# with the next geometry's closed form, so that the agreement is not mistaken for
# an insensitive one.

# %%
GEOM_NAME = {0: "slab", 1: "cylinder", 2: "sphere"}


def geometry_study(ref_shift=0):
    """eta at nu_B = 1 against the classical closed form of geometry (geom + ref_shift) % 3."""
    rows, dev = [], 0.0
    for geom in (0, 1, 2):
        for phi in (1.0, 5.0):
            p, c = pellet(1.0, 1.0, phi, n=N_PELLET, geom=geom, tol=TOL)
            ref = eta_classical((geom + ref_shift) % 3, phi)
            rows.append((geom, phi, p.eta(c), ref))
            dev = max(dev, abs(p.eta(c) / ref - 1.0))
    return {"geometry_eta_vs_closed_form_dev": dev, "_rows": rows}


V5 = geometry_study()
print(f"{'geometry':>10} {'nu':>3} {'phi':>5} {'eta solved':>12} {'closed form':>12} {'rel dev':>9}")
for geom, phi, e, ref in V5["_rows"]:
    print(f"{GEOM_NAME[geom]:>10} {geom:3d} {phi:5.1f} {e:12.8f} {ref:12.8f} {abs(e / ref - 1):9.1e}")
print(f"worst relative deviation over the three geometries ({N_PELLET} cells): "
      f"{V5['geometry_eta_vs_closed_form_dev']:.1e}")

# %% [markdown]
# ### Claims
#
# Every number the page states, with how it is checked. Values from solves are
# held to `FIXED`: the solves are converged to a Newton step of 1e-12 of the
# surface concentration, far inside it; the departures that are small differences
# of converged quantities (Jackson's finite-Kn departure, the Knudsen-limit
# effectiveness error, the slope departures) to 1e-4. Residuals, identities, two-route gaps,
# orders, rates and the blindness shifts are bounds with their basis.

# %%
FIXED = 1e-6
STUDIES = {"eq85": E85, "symbolic": SY, "isobaric": ISO, "jackson": JK, "sweep": SW,
           "headline": HD, "map": MAP, "scaling": SC, "perturbed": PT, "closure": V1,
           "grid": V2, "blindness": V3, "invariance": V4, "geometry": V5}
ALL = {k: v for s in STUDIES.values() for k, v in s.items() if not k.startswith("_")}
BOUNDS = {
    "eq85_scaling_dev": dict(max=1e-14, basis="identity",
        note="ratios of eq. (85) at 4T, 4M and 4d0 against 2, 1/2 and 4"),
    "graham_law_residual": dict(max=1e-10, basis="tolerance",
        note="identity of eqs. (82)+(85) at grad p = 0; Newton tolerance"),
    "graham_flux_ratio_dev": dict(max=1e-10, basis="tolerance",
        note="solved N_B/N_A against -sqrt(M_A/M_B); Newton tolerance"),
    "isobaric_fv_vs_closed_form": dict(max=1e-5, basis="theory",
        note="second-order discretisation error of the 200-cell slab against eq. (109) integrated"),
    "jackson_dev_x_kn_spread_1e5_to_1e8": dict(max=1e-2, basis="theory",
        note="the departure from sqrt(nu_B) falls as 1/Kn: dev x Kn constant from Kn = 1e5 to 1e8"),
    "jackson_shift_b_offdiag_flipped": dict(max=1e-5, basis="chosen",
        note="of the order of the finite-Kn departure (4e-7); the defects the check catches move "
             "it by more than 1e-2"),
    "jackson_shift_slab": dict(max=1e-9, basis="chosen",
        note="geometry leaves the Knudsen-limit centre pressure unchanged with complete conversion"),
    "jackson_shift_five_cells": dict(max=1e-6, basis="chosen",
        note="a five-cell grid leaves the Knudsen-limit centre pressure within the finite-Kn departure"),
    "jackson_invariance_spread": dict(max=1e-6, basis="chosen",
        note="Jackson's result contains no geometry, phi or grid; 1e-6 is the finite-Kn departure's order"),
    "bosanquet_p0_vs_closed_form": dict(max=1e-10, basis="tolerance",
        note="nu_B D_A/D_B is exact once c_A(0) -> 0 (complete conversion at phi = 30); Newton tolerance"),
    "eta_extrapolation_worst": dict(max=1e-4, basis="theory",
        note="finest level (4800 cells) against the extrapolated value; second order"),
    "closed_form_error_map_dev": dict(max=1e-14, basis="identity",
        note="both sides built from d_scalar"),
    "scaling_rate_knudsen": dict(min=-1.02, max=-0.98, basis="theory",
        note="departure of the local slope from 1 falls as 1/Kn under Knudsen control"),
    "scaling_rate_bulk": dict(min=0.98, max=1.02, basis="theory",
        note="departure of the local slope from 0 falls as Kn under bulk control"),
    "scaling_slope_vs_closed_form": dict(max=1e-6, basis="theory",
        note="finite-volume slopes (120 cells) against eq. (109) integrated; discretisation"),
    "dgm_vs_eq109_correct_rho": dict(max=1e-12, basis="identity",
        note="N_B = -nu_B N_A on every face of the scheme; Newton tolerance"),
    "dgm_vs_eq110_wrong_rho": dict(min=0.1, basis="chosen",
        note="a wrong flux ratio is visible: 0.1 lies orders of magnitude above the identity gap"),
    "grid_convergence_order": dict(min=1.95, max=2.1, basis="theory",
        note="order of the finest triple (1600 to 6400 cells); orders fall towards 2 along the ladder"),
    "grid_ladder_over_roundoff": dict(min=1e3, basis="theory",
        note="finest successive difference over the Newton tolerance"),
    "geometry_eta_vs_closed_form_dev": dict(max=1e-5, basis="theory",
        note="second-order discretisation error of 600 cells against the classical closed forms"),
}
# small differences of converged quantities, where the Newton tolerance (1e-12 of
# the concentration) reaches about 1e-6 of the value: held to 1e-4
SLOW = {k: 1e-4 for k in ("jackson_centre_pressure_dev", "bosanquet_eta_error_knudsen",
                          "scaling_slope_dev_knudsen", "scaling_slope_dev_bulk",
                          "scaling_slope_dev_knudsen_T400", "scaling_slope_dev_bulk_T400",
                          "scaling_slope_dev_knudsen_D1p3", "scaling_slope_dev_bulk_D1p3")}
for k, v in ALL.items():
    if k in BOUNDS:
        claims.bound(k, v, **BOUNDS[k])
    elif k == "symbolic_nonzero_residuals":
        claims.exact(k, v, note="residuals of eq. (109) against eqs. (82), (110), the error form "
                                "and Graham's ratio")
    elif k == "bosanquet_eta_error_decreasing_in_kn":
        claims.count(k, SW["_drops"], 0.0, min_margin=1e-9, above=True,
                     note=f"steps between consecutive Kn of the sweep, of {len(KN_SWEEP) - 1}, "
                          "where the effectiveness error falls")
    elif k == "bosanquet_eta_error_grows_with_phi":
        claims.count(k, HD["_rises"], 0.0, min_margin=1e-4, above=True,
                     note="steps phi 10 -> 30 -> 100, of 2, where the bulk-limit error rises")
    else:
        claims.value(k, v, rtol=SLOW.get(k, FIXED))
print(f"{len(claims.metrics)} claims registered")

# %% [markdown]
# ### The defect-injection table
#
# Every claim above needs a row that fails its own check. Each row injects one
# defect through a call argument (eq. 85's mass exponent, eq. 87's off-diagonals,
# eq. 109's stoichiometric ratio, the source coefficient, the face-composition
# rule, the Newton tolerance, the pressure dependence of the pair diffusivities,
# the pore-size dependence of the Knudsen diffusivity, the closed form compared
# with) or changes a physical condition the claim depends on ($\phi$, $\nu_B$,
# $\mathrm{Kn}$), and reruns the whole study the claim comes from. A study that
# cannot be completed under a defect records `nan` for its claims.

# %%
ROWS = [
    ("eq. (85) mass exponent +1/2: scalings", eq85_study, dict(mass_exp=-KN_EXP)),
    ("eq. (109) read with x_A/rho for rho x_A: symbolic", symbolic_study, dict(nu_swap=True)),
    ("eq. (85) mass exponent +1/2: isobaric slab", isobaric_study, dict(mass_exp=-KN_EXP)),
    ("eq. (87) off-diagonals x 0.9: isobaric slab", isobaric_study, dict(b_offdiag=0.9)),
    ("eq. (85) mass exponent +1/2: Jackson limit", jackson_study, dict(mass_exp=-KN_EXP)),
    ("eq. (85) mass dependence dropped: Jackson limit", jackson_study, dict(mass_exp=0.0)),
    ("source term with 1 B per A, fluxes with nu_B: Jackson limit", jackson_study,
     dict(src_one=True)),
    ("conversion incomplete (phi = 3): Jackson limit", jackson_study, dict(phi=3.0)),
    ("the three hidden defects at Kn = 1", blindness_study, dict(kn=1.0)),
    ("the three hidden defects at phi = 3", blindness_study, dict(phi=3.0)),
    ("invariance over phi = 3, 40, 80", invariance_study, dict(phis=(3.0, 40.0, 80.0))),
    ("nu_B = 1/2: sweep", sweep_study, dict(nu_b=0.5)),
    ("source term with 1 B per A: sweep", sweep_study, dict(nu_src=1.0)),
    ("eq. (85) mass exponent +1/2: quoted errors", headline_study, dict(mass_exp=-KN_EXP)),
    ("nu_B = 1/2: quoted errors", headline_study, dict(nu_b=0.5)),
    ("interior face composition from the left cell: quoted errors", headline_study,
     dict(face_rule="left")),
    ("quoted errors on ladders of 20 to 320 cells", headline_study,
     dict(ns=(20, 40, 80, 160, 320))),
    ("eq. (109) read with x_A/rho: closed-form map", closed_map_study, dict(nu_swap=True)),
    ("pair diffusivities held at their 1 atm value: scaling", scaling_study,
     dict(pair_const=True)),
    ("Knudsen diffusivity held at the named pore size: scaling", scaling_study,
     dict(dk_const=True)),
    ("eq. (87) off-diagonals x 0.9 in the finite-volume slab: scaling", scaling_study,
     dict(b_offdiag=0.9)),
    ("pair diffusivities held at their 1 atm value: T and D substitutions", perturbed_study,
     dict(pair_const=True)),
    ("eq. (87) off-diagonals x 0.9: closures", closure_study, dict(b_offdiag=0.9)),
    ("eq. (109) read with x_A/rho: closures", closure_study, dict(nu_swap=True)),
    ("wrong-ratio test at nu_B = 1, where eq. (110) is exact: closures", closure_study,
     dict(nu_wrong=1.0)),
    ("interior face composition from the left cell: grid order", grid_study,
     dict(face_rule="left")),
    ("Newton tolerance 1e-12 -> 1e-6: grid ladder", grid_study, dict(tol=1e-6)),
    ("closed form of the next geometry: geometry", geometry_study, dict(ref_shift=1)),
]
BASE = {eq85_study: E85, symbolic_study: SY, isobaric_study: ISO, jackson_study: JK,
        blindness_study: V3, invariance_study: V4, sweep_study: SW, headline_study: HD,
        closed_map_study: MAP, scaling_study: SC, perturbed_study: PT, closure_study: V1,
        grid_study: V2, geometry_study: V5}       # each study's baseline result
claims.define_rows([r[0] for r in ROWS])
BREAK = []
for row, fn, kw in ROWS:
    keys = [k for k in BASE[fn] if not k.startswith("_")]
    try:
        out = fn(**kw)
    except Exception as exc:                     # noqa: BLE001 - a failed study is a result
        print(f"  {row}: study raised {type(exc).__name__}: {str(exc)[:80]}; its claims record nan")
        out = {k: np.nan for k in keys}
    for k in keys:
        BREAK.append((row, k, claims.breaks(row, k, out[k])))

bt = pd.DataFrame(BREAK, columns=["row", "claim", "caught"])
caught = bt[bt.caught].groupby("row", sort=False).claim.apply(list)
print(f"{len(ROWS)} rows defined, {bt.row.nunique()} in the table; claims each row fails:")
for row, names in caught.items():
    print(f"  {row}:")
    print("      " + "\n      ".join(textwrap.wrap(", ".join(names), 96)))
idle = [r for r in bt.row.unique() if r not in caught.index]
print(f"\nrows that fail no claim: {idle if idle else 'none'}")
missing = [n for n in claims.metrics if not claims.metrics[n].get("structural")
           and not bt[(bt.claim == n) & bt.caught].shape[0]]
print(f"claims no row fails: {missing if missing else 'none'}")

# %% [markdown]
# ## What pymrm adds
#
# pymrm supplies the transport skeleton of the binary solver: `construct_grad`
# for the face gradients, `construct_div` with `nu` for slab, cylinder or sphere,
# `construct_boundary_value_matrices` for the compositions on the boundary faces
# (the pellet-centre pressure is read there), `NumJac` with a cell-and-species
# stencil and `newton` with bounds; and `pymrm.checks.observed_orders` for the
# refinement ladders. The three flux closures, the closed forms and the studies
# are the page's own.
#
# What the page adds to the review, which states the Bosanquet formula, calls it
# "very restricted in its applicability" and prints no estimate of the error:
#
# - **A closed form for the error.** $D_{\rm Bos}/D_1 = 1 + (\rho-1)x_1/(1+\mathrm{Kn})$,
#   derived from the review's own eqs. (109) and (110): the error is set by the
#   flux ratio and decays as $1/(1+\mathrm{Kn})$.
# - **The consequence, measured.** Bosanquet is exact under Knudsen control and
#   worst under bulk control, on the effectiveness factor of a reacting pellet and
#   on the flux of an isobaric counter-diffusion cell.
# - **The two failure modes separated.** Condition (iii) is "the viscous term of
#   eq. (109) vanishes"; an isobaric counter-diffusion cell satisfies it and breaks
#   (ii) through Graham's law, while a reacting pellet with a net mole change
#   breaks both. The two close simultaneously only for equal molar masses.
# - **The internal pressure as an observable.** Under bulk control the dusty gas
#   model and the Bosanquet closure disagree about the pellet-centre pressure
#   (Results 3); under Knudsen control they agree. Neither number is in the review.
# - **Jackson's $\sqrt{\nu_B}$**, quoted in the review as an aside about mechanical
#   strength, used as the page's one non-algebraic test, with its sensitivity and
#   its blind spots measured.
#
# What the page does **not** do: no experimental comparison (the review's Fig. 44
# is ternary and belongs to `A4.3`); no viscous flow, no surface diffusion, no
# pore-size distribution; no multicomponent dusty gas model; and no dimensional
# catalyst, because the review prints none.

# %% [markdown]
# ## Reuse
#
# **To adapt this page to your system**, change in order:
#
# | Change | Where |
# |---|---|
# | Your own pore size, temperature and gas | `knudsen_D(d0, T, M)`; then `Kn = D_AB/D_K` decides whether any of this matters |
# | A different stoichiometry | `nu_b` in `pellet(...)`; mass conservation fixes $M_B = M_A/\nu_B$, so `dbm` follows |
# | A different flux ratio (not from a reaction) | `rho=` on `Cell`; the isobaric cell shows the Graham case |
# | A porosity and tortuosity | `EPS_POROSITY` and `TAU_CYL`, whose ratio is `EPS_TAU`; eq. (83) scales $\mathrm{D}^e_{ij}$ and eq. (85) $\mathrm{D}^e_{iM}$ by the same factor, so `Kn` is unchanged and every dimensionless result carries over |
# | More than two species | eq. (87) for general $n$ with `np.linalg.solve`: page `A4.3` (`DustyGas`) |
# | Viscous flow | the third term of eq. (109), which needs $B_0$, $\eta$ and the $\alpha'_i$: page `A4.3` |
# | A different geometry | `geom=` (0 slab, 1 cylinder, 2 sphere), which is `nu` in `construct_div` |
#
# **The reusable pattern** is several flux closures against one skeleton: a
# `flux(c)` method, `div @ flux`, and `NumJac` over the whole divergence. Any `S3`
# or `S9` page can swap constitutive models this way and test their agreement, and
# the wrong-flux-ratio comparison of V1, with its break row at $\nu_B = 1$, is the
# template for showing that such an agreement is not vacuous.
#
# **A screening rule that falls out of the page.** Before reaching for
# $1/D = 1/D_K + 1/D_{AB}$, compute $\mathrm{Kn} = \mathrm{D}_{AB}/D_K$ and the
# flux ratio $\rho$ your problem imposes. The relative error in the effective
# diffusivity is $(\rho-1)x/(1+\mathrm{Kn})$. If $\rho = 1$ (equimolar
# counter-diffusion of equal masses) it is exact. If $\mathrm{Kn} \gg |\rho-1|$
# it is negligible. Otherwise use eq. (109), which costs one extra term.
#
# **Related pages:** `A4.2` Maxwell–Stefan vs Fick (the same review, bulk
# mixtures, no pore) · `A4.3` dusty gas model, the $n$-component matrix form ·
# `A4.9` Duncan–Toor ternary diffusion · `B1.1` Thiele modulus and the
# effectiveness factor this page perturbs · `B1.4` Weisz–Prater · `H1.7`
# solution–diffusion membranes.
#
# ## References
#
# Krishna, R. and Wesselingh, J. A. (1997). *The Maxwell–Stefan approach to mass
# transfer.* Chemical Engineering Science **52**(6), 861–911.
# [doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7).
# The document read: every equation on this page is its eq. (82), (83), (84),
# (85), (105), (107), (109) or (110), read off 600 dpi renders.
#
# Bosanquet, C. H. (1944). *British TA Report BR-507.* The origin of eq. (110).
# **Not consulted**: a wartime report, pre-DOI, not on disk and not reachable.
# The relation is used here as the review prints and names it.
#
# Knudsen, M. (1909). *Die Gesetze der Molekularströmung und der inneren
# Reibungsströmung der Gase durch Röhren.* Annalen der Physik **28**, 75–130. The
# origin of free-molecule flow through a tube. **Not consulted**; eq. (85) is
# taken from the review, which attributes it to Jackson (1977) and Mason &
# Malinauskas (1983).
#
# Jackson, R. (1977). *Transport in Porous Catalysts.* Elsevier. The source of the
# pellet-centre pressure result and of the dusty-gas derivation. **Not
# consulted**; both reach this page through the review.
#
# Mason, E. A. and Malinauskas, A. P. (1983). *Gas Transport in Porous Media: The
# Dusty-Gas Model.* Elsevier. **Not consulted**; page `A4.3`.
#
# Remick, R. R. and Geankoplis, C. J. (1974). *Ternary diffusion of gases in
# capillaries in the transition region between Knudsen and molecular diffusion.*
# Chemical Engineering Science **29**, 1447–1455. The measurements in the review's
# Fig. 44. **Not consulted, and not used here.**

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises if any check fails.
claims.report();   # the semicolon keeps the notebook from echoing the payload
