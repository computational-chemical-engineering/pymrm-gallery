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
# title: "The dusty gas model"
# description: "The n-component form of pore diffusion: why one effective diffusivity per species cannot be enough, what the viscous term is worth, and how the composition window for uphill diffusion narrows with pore size."
# categories: [sec:A, struct:S9, tier:T1, data:tier6, phase:gas]
# date: 2026-08-01
# ---
#
# # The dusty gas model
#
# **Catalog ID:** `A4.3` · **Structure:** `S9` (coupled multicomponent transport) ·
# **Tier:** T1
#
# Inside a pore a molecule collides with other molecules and with the wall. The
# **dusty gas model** handles both at once by treating the wall as a giant
# motionless pseudo-species, "dust", and writing the Maxwell–Stefan equations
# for the $(n{+}1)$-component mixture. What comes out is a *matrix*: the flux of
# each species depends on every other species' flux, and on the pressure
# gradient.
#
# Page [`A4.4`](../A4.4-knudsen-bosanquet/index.ipynb) takes the **scalar** limit of
# that, one effective diffusivity per species (the Bosanquet interpolation), and
# measures its error. This page is the matrix the scalar is an approximation *of*,
# plus the term `A4.4` sets to zero: **viscous flow**.
#
# Three things fall out that the scalar form cannot express at all:
#
# - an **inert species with zero flux and a concentration that varies by tens of
#   percent** across a pellet; any scalar closure with a finite positive
#   diffusivity forces that variation to be exactly zero;
# - **uphill diffusion**, which no pore size rules out: what shrinks as the pore
#   shrinks is the *composition window* in which it happens, and this page
#   measures how;
# - a **reaction-generated internal pressure**, and how much of it viscous flow
#   bleeds away: under a percent at every pore size.
#
# Krishna and Wesselingh's 1997 review prints all of it and states the last two
# qualitatively. Everything below is that review's equations, implemented in
# pymrm and made numerical.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.27.7"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A4.3-dusty-gas-model"  # equals the page directory

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
from scipy.optimize import minimize_scalar
from IPython.display import Markdown, display
from pymrm.checks import find_roots, observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import (RGAS, SYMMETRY, DustyGas, dirichlet, knudsen_D, pair_matrix,
                   capillary as _capillary, pellet as _pellet, capillary_bvp as _capillary_bvp,
                   ms_bvp as _ms_bvp)

PAGE_A42 = "A4.2-maxwell-stefan-vs-fick"      # cross-page: the printed pair diffusivities
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# **Why a matrix and not a diffusivity.** In the bulk limit a species exchanges
# momentum only with the *other species*, so the Maxwell–Stefan friction terms
# couple all the fluxes and the resulting resistance matrix is **singular**:
# there is no absolute reference for velocity, only relative motion. In the
# Knudsen limit each species exchanges momentum only with the *wall*, so the
# fluxes decouple completely and the matrix is diagonal. The dusty gas model
# holds both terms simultaneously, and it is the wall term that makes the matrix
# invertible at all. That is why the model can be written explicitly as
# $(N) = -[B^e]^{-1}(\nabla p)$, and why the bulk limit has to be approached
# rather than taken.
#
# **And that is a theorem, not a numerical accident.** Delete the wall term from
# eq. (87) and sum a *column* $j$ of what is left:
#
# $$
# \sum_{i} B_{ij}\Big|_{\text{no wall}}
# = \underbrace{\sum_{k\neq j}\frac{x_k}{\mathrm{D}_{jk}}}_{B_{jj}}
# \;-\; \underbrace{\sum_{i\neq j}\frac{x_i}{\mathrm{D}_{ij}}}_{\text{off-diagonals}}
# \;=\; 0 ,
# $$
#
# because $\mathrm{D}_{ij} = \mathrm{D}_{ji}$. So $(1,\dots,1)$ is an **exact left
# null vector** of the wall-free matrix, for *any* composition and *any*
# $n \ge 2$, and the matrix is singular identically. Restore the wall term and
# every column sum becomes $1/\mathrm{D}^e_{jM} > 0$; since the off-diagonals are
# all of one sign, the matrix is then **strictly column-diagonally dominant**,
# hence nonsingular by Levy–Desplanques. The margin by which it is invertible *is*
# the Knudsen resistance. Both halves are checked against the assembled
# `b_matrix` below, and deleting the wall term in the solver leaves no solution.
#
# The physical reading of the same identity: summing eq. (82) over $i$ gives
# $\nabla p/RT = \sum_i N_i/\mathrm{D}^e_{iM}$, which without the wall term forces
# $\nabla p = 0$ (Graham's law), so the wall-free system is consistent only on
# that subspace and has nothing to say off it.
#
# **What "dust" means.** The wall is species $n{+}1$: infinitely heavy, uniformly
# distributed, motionless. Friction against it is then just another
# Maxwell–Stefan term, and eliminating it from the $(n{+}1)$-component equations
# leaves $n$ equations in the $n$ real fluxes. The review attributes the
# construction to Mason & Malinauskas (1983) and Jackson (1977).
#
# **Three mechanisms, and which two are here.** The review's Fig. 40 separates
# bulk diffusion, Knudsen diffusion and *surface* diffusion of adsorbed species.
# Only the first two are in the dusty gas model and only those two are on this
# page; the review is explicit that bulk and Knudsen "occur together and it is
# prudent to take both mechanisms into account rather than assume that one or
# other mechanism is 'controlling'".
#
# **Scope, and what belongs elsewhere.**
#
# | | page | what it owns |
# |---|---|---|
# | bulk multicomponent, no wall | [`A4.2`](../A4.2-maxwell-stefan-vs-fick/index.ipynb) | $[D] = [B]^{-1}[\Gamma]$, Maxwell–Stefan against generalized Fick |
# | the *scalar* pore closure | [`A4.4`](../A4.4-knudsen-bosanquet/index.ipynb) | eq. (109), the Bosanquet relation eq. (110), and its error |
# | **the $n$-component matrix + viscous flow** | **this page** | eqs. (82), (86), (87), (98) to (101), (103), (105) to (107) |
# | uphill diffusion in free space | [`A4.9`](../A4.9-duncan-toor/index.ipynb) | the Duncan & Toor two-bulb measurements |
#
# `A4.4` uses the binary $2\times2$ $[B^e]$ as the yardstick for its scalar
# closure and leaves the general $n$-component form, the viscous term and the
# review's Fig. 44 to this page. This page does not re-derive `A4.4`'s scalar
# result; it checks against it (`V5`).

# %% [markdown]
# ## The published model
#
# **Which document each equation comes from.** The dusty gas model originates with
# **Mason & Malinauskas (1983)** and **Jackson (1977)**; the Knudsen and Bosanquet
# results with **Knudsen (1909)** and a 1944 British wartime report. Every
# equation below is taken from **Krishna & Wesselingh (1997)**, which prints all
# of them with attribution and is the reference this page cites; equation numbers
# are that review's. They are read from renders of the printed page, not from the
# PDF text layer, which is an OCR of a scan that renders eq. (86)'s prefactor as
# $c_t/RT$ where the page prints $1/RT$ and turns eq. (109)'s stoichiometric
# subscript $\nu_1$ into $\nu_{1i}$.
#
# **The dusty gas equations** (eq. 82), for $i = 1,\dots,n$:
#
# $$
# -\frac{1}{RT}\nabla p_i \;=\; \sum_{j=1}^{n}\frac{x_j N_i - x_i N_j}{\mathrm{D}^e_{ij}}
# \;+\; \frac{N_i}{\mathrm{D}^e_{iM}} .
# $$
#
# **The porous matrix** (eqs. 83 and 84): $\mathrm{D}^e_{ij} = (\varepsilon/\tau)\,
# \mathrm{D}_{ij}$, with $\tau = 1$ for a cylindrical pore.
#
# **The Knudsen diffusivity** (eq. 85), attributed to Jackson (1977) and Mason &
# Malinauskas (1983), with $d_0$ the pore **diameter** and the square root "the
# velocity of motion":
#
# $$
# \mathrm{D}^e_{iM} = \frac{\varepsilon}{\tau}\,\frac{d_0}{3}\sqrt{\frac{8RT}{\pi M_i}} .
# $$
#
# **The matrix form** (eqs. 86 and 87), which is what this page is about:
#
# $$
# (N) = -\frac{1}{RT}\,[B^e]^{-1}(\nabla p), \qquad
# B^e_{ii} = \frac{1}{\mathrm{D}^e_{iM}} + \sum_{k\neq i}\frac{x_k}{\mathrm{D}^e_{ik}},
# \qquad
# B^e_{ij(i\neq j)} = -\frac{x_i}{\mathrm{D}^e_{ij}} .
# $$
#
# **Viscous flow** (eqs. 90, 91, 95, 98 to 101). Under a pressure gradient the
# fluid also flows bodily, $v = -(B_0/\eta)\nabla p$, with $B_0 = d_0^2/32$ for a
# cylindrical pore (eq. 91, Poiseuille). Allowing that flow to be *selective*
# through factors $\alpha_i$ (eq. 95) and folding it into the species velocities
# gives the review's final working form,
#
# $$
# -\frac{c_i}{RT}\nabla_{T,p}\mu_i - \frac{c_i}{RT}\bar V_i\nabla p
# - \alpha'_i c_i \frac{B_0}{\eta\,\mathrm{D}^e_{iM}}\nabla p
# - c_i z_i \frac{F}{RT}\nabla\Phi
# = \sum_{j=1}^{n}\frac{x_j N_i - x_i N_j}{\mathrm{D}^e_{ij}} + \frac{N_i}{\mathrm{D}^e_{iM}} ,
# $$
#
# which is eq. (101). For an ideal gas with no electric field the first two terms
# collapse to $\nabla p_i/RT$ exactly, and with the "viscous diffusivity"
# $\mathrm{D}_{\rm visc} \equiv c_t B_0 RT/\eta$ of eq. (106) the third becomes
# $\alpha'_i x_i (\mathrm{D}_{\rm visc}/\mathrm{D}^e_{iM})\nabla p / RT$. That is
# the only form used below, and eq. (82) is its $\mathrm{D}_{\rm visc} = 0$ case.
#
# **Three consequences the review prints, and this page uses as targets.**
#
# *Eq. (105)*, summing eq. (101) over all species:
#
# $$
# \nabla p \;=\; -\,\frac{\sum_i N_i/\mathrm{D}^e_{iM}}
# {\dfrac{1}{RT}\left(1 + \mathrm{D}_{\rm visc}\sum_i \alpha'_i x_i/\mathrm{D}^e_{iM}\right)} .
# $$
#
# *Eq. (107), Graham's law*, the same sum with $\nabla p = 0$, using eq. (85)'s
# $\mathrm{D}^e_{iM}\propto M_i^{-1/2}$:
#
# $$
# \sum_{i=1}^{n} N_i\sqrt{M_i} = 0 \qquad (\nabla p = 0;\ \text{gaseous mixtures}).
# $$
#
# *Eq. (103)*, a single component, for which $\alpha'_1 = 1$:
#
# $$
# N_1 = -\left(\frac{\mathrm{D}^e_{1M}}{RT} + c_1\frac{B_0}{\eta}\right)\nabla p ,
# $$
#
# Knudsen slip plus Poiseuille flow. Eqs. (105) and (103) are *derivable* from
# eq. (101) and are used below as checks on the transcription rather than as
# inputs.
#
# **What the review states in words and this page turns into numbers.**
#
# - p. 887: "the fluxes $N_i$ are independent of system pressure and pore size"
#   under bulk control, and "directly proportional to the system pressure and pore
#   size" under Knudsen control. Four exponents, 0 and 1.
# - p. 892: "Finite pressure gradients can be generated inside a porous catalyst
#   when there is a net change in the number of moles", and, from Jackson (1977),
#   for $A \to \nu_B B$ with complete conversion "the pressure at the centre of the
#   catalyst pellet … is $p_0 = \sqrt{\nu_B}\,p$ … for $\nu_B = 2$, we have a 40%
#   increase in pressure".
# - p. 893: "we also note that neglect of the viscous flow contribution … is not
#   very serious; this result is typical (Haynes, 1978)". **How not-very-serious
#   is not stated.**
# - p. 893: "Kaza and Jackson (1980) demonstrate the possibility of uphill
#   diffusion within a catalyst particle, a phenomena impossible to explain with
#   say the Bosanquet formula (110)". **No condition on it is given.**
#
# | Symbol | Code | Meaning |
# |---|---|---|
# | $[B^e]$ | `DustyGas.b_matrix` | eq. (87) |
# | $\mathrm{D}^e_{iM}$ | `knudsen_D`, `diM` | eq. (85) |
# | $\mathrm{D}^e_{ij}$ | `dij_ct`/$c_t$ | eq. (83), stored as the pressure-independent product $\mathrm{D}_{ij}c_t$ |
# | $\mathrm{D}_{\rm visc}$ | `d_visc` | eq. (106) |
# | $\alpha'_i$ | `alpha` | eq. (100), viscous selectivity; 1 throughout |
# | $B_0$ | `B0` | eq. (91), $d_0^2/32$ |

# %% [markdown]
# ## Parameters and assumptions
#
# **Where the numbers come from, and where they do not.** The review prints no
# binary pair diffusivities for any porous medium and no porosity for any of its
# worked configurations. So:
#
# - the **pair diffusivities** are the review's own printed values for the
#   H₂/N₂/CO₂ system on its p. 872, loaded cross-page from the
#   [`A4.2`](../A4.2-maxwell-stefan-vs-fick/index.ipynb) dataset, as `A4.4` does.
#   The review introduces them as values that "can be estimated from the kinetic
#   gas theory" and prints **no temperature and no pressure** beside them; the
#   35.2 °C and 1 atm used here are Duncan & Toor's own conditions, from
#   [`A4.9`](../A4.9-duncan-toor/index.ipynb)'s sidecar, and are labelled as such
#   wherever they appear;
# - $\varepsilon/\tau = 1$ is eq. (84)'s printed $\tau = 1$ for a cylindrical pore
#   together with an **assumed** $\varepsilon = 1$, one straight open pore. It
#   multiplies every diffusivity equally, so no dimensionless result on this page
#   depends on it;
# - the **pellet results are dimensionless**, in the groups
#   $(\mathrm{Kn}, \nu_B, \phi, \mathrm{D}_{\rm visc}/\mathrm{D}^e_{AM})$,
#   following `B3.1`: the review prints no pellet radius and no rate constant.
#
# **The Knudsen number** is `A4.4`'s, so that the two pages' abscissae are the
# same quantity: $\mathrm{Kn} \equiv \mathrm{D}^e_{ij}/\mathrm{D}^e_{iM}$; large is
# Knudsen control, small is bulk control.
#
# **Assumptions**, all the review's: isothermal ideal gas; no surface diffusion; a
# single pore size; no electric field; $\alpha'_i = 1$ (non-selective viscous
# flow, the review's own default outside ultrafiltration); and, where there is a
# reaction, first order in $A$ and at steady state.
#
# **What this page is not.** It is **data tier 6**: nothing on it is compared with
# a measurement. The review's one experimental comparison for the dusty gas model,
# its Fig. 44 against Remick & Geankoplis's (1974) He/Ne/Ar capillary data, is not
# used here.

# %% [markdown]
# ## The data
#
# **There is no dataset here in the experimental sense.** Two collections of
# printed numbers are used, and neither is a measurement.
#
# 1. **`krishna-wesselingh-1997-dgm-constants.csv`**: the structural constants of
#    the dusty gas model (the $d_0/3$ and $M^{-1/2}$ of eq. 85, the $\tau = 1$ of
#    eq. 84, the off-diagonal sign of eq. 87, the 32 of eq. 91, the $\sqrt{M_i}$
#    of eq. 107, the $\alpha'_1 = 1$ of eq. 103) together with the **four
#    exponents** the review states in words on p. 887. Those four are stored as
#    numbers so a computed log–log slope is compared with something that came off
#    the page rather than with something typed into this notebook. The molar
#    masses are IUPAC values, **not** printed in the review, and flagged as such
#    in the file.
# 2. **The review's own binary pair diffusivities** (p. 872), the `ideal` rows
#    `D12`, `D13`, `D23` of `A4.2`'s `krishna-wesselingh-1997-worked-examples.csv`.
#    What `A4.2` establishes about these rows: they are printed to two or three
#    significant figures, and recomputing the review's ideal-gas $[D]$ matrix from
#    them reproduces the printed elements to about one percent, the rounding of the
#    print. That rounding affects every dimensional number on this page in the
#    same way as a slightly different gas would, and none of the page's checks,
#    which compare routes handed the same three values. `A4.2`'s sign finding
#    concerns the non-ideal (liquid) example, whose rows are not used here.

# %%
K = load_data("krishna-wesselingh-1997-dgm-constants.csv", page=PAGE)
K_meta = load_meta("krishna-wesselingh-1997-dgm-constants.csv", page=PAGE)
P = {(r.group, r.symbol): r.value for r in K.itertuples()}

wx = load_data("krishna-wesselingh-1997-worked-examples.csv", page=PAGE_A42)
wx_meta = load_meta("krishna-wesselingh-1997-worked-examples.csv", page=PAGE_A42)
W = {(r.example, r.symbol): r.value for r in wx.itertuples()}

print("page dataset :", cite_data(K_meta))
print("cross-page   :", cite_data(wx_meta))

# --- structural constants ----------------------------------------------------
EPS_TAU   = P[("eq83", "eps_over_tau")] / P[("eq84", "tau_cylindrical")]  # eps assumed 1
KN_DEN    = P[("eq85", "knudsen_denominator")]        # eq. (85): d_0 / 3
KN_EXP    = P[("eq85", "knudsen_mass_exponent")]      # eq. (85): M^(-1/2)
B_OFF     = P[("eq87", "b_offdiag_sign")]             # eq. (87): B_ij = -x_i/D_ij
POIS_DEN  = P[("eq91", "poiseuille_denominator")]     # eq. (91): B_0 = d_0^2/32
ALPHA_DEF = P[("eq101", "viscous_selectivity_default")]
GR_EXP    = P[("eq107", "graham_mass_exponent")]      # eq. (107): sqrt(M)
EXP_KN_P  = P[("p887", "knudsen_flux_pressure_exponent")]    # stated 1
EXP_KN_D  = P[("p887", "knudsen_flux_poresize_exponent")]    # stated 1
EXP_BK_P  = P[("p887", "bulk_flux_pressure_exponent")]       # stated 0
EXP_BK_D  = P[("p887", "bulk_flux_poresize_exponent")]       # stated 0

# --- the review's own system, for the dimensional illustrations --------------
# The three pair diffusivities are printed on p. 872 with no temperature and no
# pressure; 35.2 C and 1 atm are Duncan & Toor's conditions (via A4.9).
T_REV, P_ATM = 308.35, 101325.0               # K, Pa  (via A4.9, not the review)
CT0 = P_ATM / (RGAS * T_REV)                  # mol/m3
SPEC = ["H2", "N2", "CO2"]
MOLAR = np.array([P[("molar", f"M_{s}")] for s in SPEC])          # kg/mol, IUPAC
D_PAIR = {(0, 1): W[("ideal", "D12")],        # H2-N2
          (0, 2): W[("ideal", "D13")],        # H2-CO2
          (1, 2): W[("ideal", "D23")]}        # N2-CO2

print(f"\n{len(K)} printed constants; {len(D_PAIR)} pair diffusivities reused from A4.2")
for (i, j), v in D_PAIR.items():
    print(f"  D_{SPEC[i]}-{SPEC[j]:4s} = {v:.3e} m2/s   (K&W p. 872)")
print(f"  eps/tau = {EPS_TAU:.1f}   (tau = 1 printed in eq. 84; eps = 1 assumed)")
print(f"  stated exponents, p. 887: Knudsen ({EXP_KN_P:.0f}, {EXP_KN_D:.0f}), "
      f"bulk ({EXP_BK_P:.0f}, {EXP_BK_D:.0f})")

# %% [markdown]
# ## PyMRM implementation
#
# ### Eq. (85), and the one place a molar mass enters
#
# `knudsen_D` in `model.py` is eq. (85) with $d_0$ the pore diameter;
# `pair_matrix` stores the pressure-independent products $\mathrm{D}_{ij}c_t$, so
# one object serves every pressure.

# %%
for s, M in zip(SPEC, MOLAR):
    print(f"  {s:4s} v_bar = {np.sqrt(8*RGAS*T_REV/(np.pi*M)):7.1f} m/s   "
          f"D^e_iM(d0 = 100 nm) = {knudsen_D(100e-9, T_REV, M, EPS_TAU, KN_DEN, KN_EXP):.3e} m2/s")

# %% [markdown]
# ### The solver
#
# `DustyGas` (in `model.py`) is one class covering every configuration on this
# page: a capillary or a pellet, any number of species, with or without the
# viscous term, with or without reaction.
#
# - **Layout `(n_cell, n_species)`**: spatial axis first, species last.
# - **Constant operators assembled once**: `construct_grad` for the face
#   gradients, `construct_div` for the balance (`nu` is the geometry: `0` slab,
#   `1` cylinder, `2` sphere), and one sparse face-composition matrix, the
#   boundary rows from `construct_boundary_value_matrices` stacked on the
#   interior face means. Only $[B^e]$, which depends on the local composition, is
#   rebuilt each iteration.
# - **`NumJac(shape, axes_diagonals=[0], axes_blocks=[-1])`**: the residual in a
#   cell reads its neighbours through the face fluxes (`axes_diagonals=[0]`) and
#   the species couple pointwise through $[B^e]$ (`axes_blocks=[-1]`).
# - **Boundary conditions on the outward normal**, $a\,\partial c/\partial n + b\,c = d$.
#   At a pellet centre (outward normal $-r$) symmetry $\partial c/\partial r = 0$
#   is `{a:1, b:0, d:0}`; an imposed composition $c = c_{\rm bnd}$ is
#   `{a:0, b:1, d:c_bnd}`.
# - **`newton`** with a step tolerance of $10^{-12}$ relative to the largest
#   concentration (round-off is near $10^{-15}$), `stagnation=True`, and
#   concentrations bounded below by zero; a solve that does not converge raises.
#   The reaction is first order in $A$, smooth through zero.
#
# Arguments exist so defects can be injected later: `b_offdiag` flips the sign of
# eq. (87)'s off-diagonals, `knudsen=False` deletes the wall term, `face_rule`
# takes the interior face composition from the left cell, and `scalar=True`
# replaces $[B^e]^{-1}$ by the inverse of its diagonal, a scalar closure used
# only for contrast.
#
# #### The wall term, and why $[B^e]$ is invertible
#
# The Background's claim, on the assembled matrix rather than on paper. For random
# compositions and $n = 2\dots6$ random *symmetric* pair matrices, the wall-free
# $[B^e]$ must have column sums zero (so $(1,\dots,1)$ is a left null vector and
# the matrix is singular), and the full $[B^e]$ must have column sums
# $1/\mathrm{D}^e_{jM}$, which is also its margin of strict column diagonal
# dominance, and hence of invertibility.

# %%
def bmatrix_study(b_offdiag=B_OFF, wall_in_free=False, wall_in_full=True, seed=0):
    """Column sums of the assembled [B^e] with and without the wall term."""
    rng = np.random.default_rng(seed)
    null, colsum, dom, rows = 0.0, 0.0, np.inf, []
    for n in [2, 3, 4, 5, 6]:
        dij = rng.random((n, n)) * 1e-4 + 1e-6
        dij = 0.5 * (dij + dij.T); np.fill_diagonal(dij, 0.0)      # symmetric, eq. (83)
        diM = rng.random(n) * 1e-6 + 1e-8
        cf = rng.random((40, n))                                   # 40 random compositions
        cf = CT0 * cf / cf.sum(axis=1, keepdims=True)
        bc = (dirichlet(np.full(n, CT0 / n)), dirichlet(np.full(n, CT0 / n)))
        Bf = DustyGas(2, 1.0, dij * CT0, diM, T_REV, bc=bc, b_offdiag=b_offdiag,
                      knudsen=wall_in_free).b_matrix(cf)
        Bw = DustyGas(2, 1.0, dij * CT0, diM, T_REV, bc=bc, b_offdiag=b_offdiag,
                      knudsen=wall_in_full).b_matrix(cf)
        # axis 0 of [B^e] is the face, axis 1 the row i, axis 2 the column j
        nv = float(np.max(np.abs(Bf.sum(axis=1)) / np.abs(Bf).sum(axis=1)))
        cs = float(np.max(np.abs(Bw.sum(axis=1) * diM[None, :] - 1.0)))
        diag = np.diagonal(Bw, axis1=1, axis2=2)
        # strict column diagonal dominance: B_jj - sum_{i != j} |B_ij|, relative to B_jj
        dm = float(np.min((2 * diag - np.abs(Bw).sum(axis=1)) / diag))
        null, colsum, dom = max(null, nv), max(colsum, cs), min(dom, dm)
        rows.append((n, nv, cs, dm))
    return {"bmatrix_wall_free_null_vector": null, "bmatrix_colsum_knudsen_dev": colsum,
            "bmatrix_column_dominance_margin": dom, "_rows": rows}


BM = bmatrix_study()
print(f"{'n':>3} {'|1^T B| / |B| (no wall)':>24} {'col sum vs 1/D_jM':>19} {'dominance / B_jj':>18}")
for n, nv, cs, dm in BM["_rows"]:
    print(f"{n:3d} {nv:24.2e} {cs:19.2e} {dm:18.3e}")

# %% [markdown]
# ### The two configurations
#
# A **capillary**, a slab with an imposed composition at each end and no reaction,
# and a **pellet**, a sphere with symmetry at the centre and an imposed
# composition at the surface, with $A \to \nu_B B$ first order in $A$
# (`capillary` and `pellet` in `model.py`; the wrappers below hand them this page's
# constants, and every keyword overrides one).
#
# The capillary is the review's Fig. 44 geometry with the H₂/N₂/CO₂ system
# substituted, because the review prints pair diffusivities for that system and
# not for the He/Ne/Ar one. The length is Duncan & Toor's capillary, 8.6 cm, via
# `A4.9`; it sets the flux scale and cancels out of every dimensionless result.
# Pellet centre and surface values are read on the boundary faces.

# %%
L_CAP = 0.086                       # m, Duncan & Toor's capillary, via A4.9
REF = dict(molar=MOLAR, dij_ct=pair_matrix(D_PAIR, CT0), T=T_REV, p=P_ATM, L=L_CAP,
           eps_tau=EPS_TAU, kn_den=KN_DEN, mass_exp=KN_EXP)   # shared by both capillary routes
XB_CAP = np.array([0.5, 0.5, 0.0])  # far end: H2/N2, no CO2
X_UP = np.array([0.0, 0.49, 0.51])  # near end, inside the uphill window


def capillary(d0, xa, xb, **kw):
    """Finite-volume capillary with this page's constants; returns (model, c)."""
    return _capillary(d0, xa, xb, **{**REF, "b_offdiag": B_OFF, "alpha": None, **kw})


def capillary_bvp(d0, xa, xb, **kw):
    """Collocation route (no pymrm operator) with this page's constants."""
    return _capillary_bvp(d0, xa, xb, **{**REF, **kw})


def pellet(nu_b, kn, phi, **kw):
    """Dimensionless pellet with this page's eq. (85) and (87) constants."""
    return _pellet(nu_b, kn, phi, **{"mass_exp": KN_EXP, "b_offdiag": B_OFF, **kw})


def mid_flux(d0, xa=X_UP, xb=XB_CAP, n_cell=100, **kw):
    """Capillary fluxes on the middle face (uniform at steady state)."""
    m, c = capillary(d0, xa, xb, n_cell=n_cell, **kw)
    return m.flux(c)[n_cell // 2]


def centre(m, c):
    """Composition on the r = 0 face of a pellet."""
    return m.faces(c)[0][0]


m0, c0 = capillary(1e-7, X_UP, XB_CAP)
print(f"capillary at d0 = 100 nm: N = {m0.flux(c0)[m0.n_cell // 2]} mol/(m2 s)")
mp, cp = pellet(2.0, 1e4, 60.0)
print(f"pellet nu_B = 2, Kn = 1e4, phi = 60: p_0/p = {centre(mp, cp).sum():.6f}")

# %% [markdown]
# ### The transcription, checked before it is used
#
# Eq. (101) is read once off a scan, and everything on this page rests on it. So
# it is checked against two *other* printed equations that must follow from it
# (the `B3.1` pattern: a transcription you can only read once is checked against
# something you can derive).
#
# Summing eq. (101) over all species makes the Maxwell–Stefan double sum vanish
# identically, because the term $(i,j)$ is the negative of the term $(j,i)$. What
# is left must be eq. (105); setting $\nabla p = 0$ in it must give eq. (107); and
# setting $n = 1$ must give eq. (103). None of these three is an input to the
# solver. `asym` makes $\mathrm{D}_{ji} \ne \mathrm{D}_{ij}$, for a defect row.

# %%
NS = 3


def transcription_study(asym=1):
    """Symbolic residuals of eqs. (105), (107), (103) derived from eq. (101)."""
    _N = sp.symbols("N1:4")                      # fluxes, sign matters, no assumptions
    _x = sp.symbols("x1:4", positive=True)
    _gpi = sp.symbols("gp1:4")                   # the individual grad p_i
    _gp, _RT = sp.symbols("gradp RT")
    _dM = sp.symbols("DM1:4", positive=True)
    _dv, _al = sp.symbols("D_visc alpha", positive=True)
    _dij = sp.Matrix(NS, NS, lambda i, j: 0 if i == j else
                     sp.Symbol(f"D{min(i, j) + 1}{max(i, j) + 1}", positive=True)
                     * (asym if i > j else 1))
    # eq. (101) for an ideal gas with no electric field:  lhs_i = rhs_i
    lhs = lambda i: -_gpi[i] / _RT - _al * _x[i] * (_dv / _dM[i]) * _gp / _RT
    rhs = lambda i: (sum((_x[j] * _N[i] - _x[i] * _N[j]) / _dij[i, j]
                         for j in range(NS) if j != i) + _N[i] / _dM[i])
    # eq. (105): sum over i; the friction double sum cancels and sum_i grad p_i = grad p
    total = sp.expand(sum(rhs(i) for i in range(NS)) - sum(lhs(i) for i in range(NS)))
    total = sp.simplify(total.subs(_gpi[2], _gp - _gpi[0] - _gpi[1]))
    friction_left = sp.simplify(total - sp.simplify(
        sum(_N[i] / _dM[i] for i in range(NS)) + _gp / _RT
        + _al * _dv * _gp / _RT * sum(_x[i] / _dM[i] for i in range(NS))))
    eq105_derived = sp.solve(sp.Eq(total, 0), _gp)[0]
    eq105_printed = -sum(_N[i] / _dM[i] for i in range(NS)) / (
        (1 / _RT) * (1 + _dv * sum(_al * _x[i] / _dM[i] for i in range(NS))))
    eq105 = sp.simplify(eq105_derived - eq105_printed)
    # eq. (107): Graham's law, from eq. (105) with grad p = 0, using eq. (85)
    _M = sp.symbols("M1:4", positive=True)
    _Kk = sp.Symbol("Kk", positive=True)         # (eps/tau)(d_0/3) sqrt(8RT/pi)
    knud = {_dM[i]: _Kk * _M[i] ** sp.Rational(-1, 2) for i in range(NS)}
    graham = sp.expand(sum(_N[i] / _dM[i] for i in range(NS)).subs(knud) * _Kk)
    eq107 = sp.simplify(graham - sum(_N[i] * sp.sqrt(_M[i]) for i in range(NS)))
    # eq. (103): one component, for which the review states alpha'_1 = 1
    _c1, _B0, _eta, _DM1 = sp.symbols("c1 B0 eta DM1", positive=True)
    one = sp.Eq(-_gp / _RT - ((_c1 * _B0 * _RT / _eta) / _DM1) * _gp / _RT, _N[0] / _DM1)
    eq103 = sp.simplify(sp.solve(one, _N[0])[0] + (_DM1 / _RT + _c1 * _B0 / _eta) * _gp)
    res = {"friction terms left by the sum": friction_left, "eq. (105)": eq105,
           "eq. (107)": eq107, "eq. (103)": eq103}
    return {"symbolic_nonzero_residuals": int(sum(r != 0 for r in res.values())),
            "_res": res, "_gp": sp.simplify(eq105_derived)}


TR = transcription_study()
for k, v in TR["_res"].items():
    print(f"  {k:32s} residual: {v}")
print("\nwhat summing eq. (101) gives for grad p:")
sp.pprint(TR["_gp"])

# %% [markdown]
# **What that settles, and what it does not.** Eqs. (105), (107) and (103) all
# drop out of eq. (101) with zero residual, so the three transcriptions are
# mutually consistent and the one equation the solver implements is the one on the
# page. The antisymmetry that kills the friction double sum is what makes Graham's
# law exact: it is a property of eq. (82)/(101), not an extra assumption, and it
# holds at *every* Knudsen number, not only in the free-molecule limit.
#
# **This is algebra.** It cannot detect an error in the pymrm assembly, in the
# boundary conditions or in the discretisation. The Validation section tests those.

# %% [markdown]
# ## Results
#
# ### 1. An inert species with zero flux and a concentration that is not flat
#
# Put an inert third species into a pellet running $A \to 2B$. At steady state it
# is not consumed and the centre is a symmetry point, so **its flux is zero
# everywhere**: that is conservation, true in any model, and not a prediction.
#
# The question is what its *concentration* does. A scalar effective-diffusivity
# closure says $N_C = -D_C\nabla c_C$ with a **finite positive** $D_C$, so
# $N_C = 0$ forces $\nabla c_C = 0$: a flat profile, necessarily. The review's
# own scalar closure eq. (109) does not satisfy that premise for an inert: with
# $\nu_C = 0$ its $x_C\nu_j/(x_j\nu_C)$ term diverges and $D_C \to 0$, so
# eq. (109) predicts *nothing* here rather than a flat profile. The statement
# below is about closures with a finite positive diffusivity, which is what "an
# effective diffusivity per species" means in practice; the diagonal of $[B^e]$
# is one of them. The dusty gas model for the inert reads
#
# $$
# -\nabla c_C = -x_C\sum_j \frac{N_j}{\mathrm{D}^e_{Cj}} \quad\text{when } N_C = 0,
# $$
#
# which is nonzero as soon as the reaction changes the number of moles. This is
# the **diffusion barrier** of `A4.9`'s triad, inside a pore. The variation is
# read between the centre and surface faces.

# %%
INERT_KN = {"kn0p1": 0.1, "kn1": 1.0, "kn10": 10.0}


def inert_study(contrast_scalar=True, **kw):
    """Inert C in an A -> 2B pellet: centre-to-surface variation of c_C."""
    out, var_s, flux, sols = {}, [], [], {}
    for tag, kn in INERT_KN.items():
        args = (2.0, kn, 20.0)
        opt = dict(x_as=0.5, n_extra=1, x_extra=0.5)
        m, c = pellet(*args, **opt, **kw)
        # the same pellet with [B^e]^-1 replaced by its diagonal: a scalar closure
        # with a finite positive D_C, run rather than asserted
        ms_, cs_ = pellet(*args, **opt, **{**kw, "scalar": contrast_scalar})
        cf, cfs = m.faces(c)[0], ms_.faces(cs_)[0]
        N = m.flux(c)
        out[f"inert_concentration_variation_{tag}"] = abs(cf[-1, 2] / cf[0, 2] - 1.0)
        var_s.append(abs(cfs[-1, 2] / cfs[0, 2] - 1.0))
        flux.append(np.abs(N[:, 2]).max() / np.abs(N[:, 0]).max())
        sols[kn] = (m, c, cf)
    out["inert_scalar_closure_variation"] = max(var_s)
    out["inert_flux_ratio"] = max(flux)
    out["_sols"] = sols
    return out


IN = inert_study()
print(f"{'Kn':>6} {'p_0/p':>9} {'c_C centre':>11} {'c_C surface':>12} {'variation':>10}")
for tag, kn in INERT_KN.items():
    m, c, cf = IN["_sols"][kn]
    print(f"{kn:6.1f} {cf[0].sum() / cf[-1].sum():9.5f} {cf[0, 2]:11.5f} {cf[-1, 2]:12.5f} "
          f"{IN[f'inert_concentration_variation_{tag}']:9.1%}")
print(f"scalar closure, variation below 1e-12: {IN['inert_scalar_closure_variation'] < 1e-12}; "
      f"|N_C|/|N_A| below 1e-10: {IN['inert_flux_ratio'] < 1e-10}")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
for kn in INERT_KN.values():
    m, c, _ = IN["_sols"][kn]
    axes[0].plot(m.x_c, c[:, 2], lw=2, label=rf"$\mathrm{{Kn}} = {kn:g}$")
    axes[1].plot(m.x_c, c.sum(axis=1), lw=2, label=rf"$\mathrm{{Kn}} = {kn:g}$")
axes[0].axhline(0.5, color="0.6", lw=1, ls="--")
axes[0].set_xlabel("r / R"); axes[0].set_ylabel(r"$c_C$ (inert) / $c_{t,\rm surface}$")
axes[0].set_title("The inert carries no flux, and a scalar closure with\n"
                  "$D_C > 0$ therefore makes these curves horizontal")
axes[1].set_xlabel("r / R"); axes[1].set_ylabel(r"$c_t$ / $c_{t,\rm surface}$")
axes[1].set_title(r"Total pressure, $A \rightarrow 2B$, $\phi = 20$")
for a in axes:
    a.legend(fontsize=8)
plt.tight_layout(); plt.show()

# %%
display(Markdown(rf"""
At $\mathrm{{Kn}} = 1$ the inert species' concentration varies by
**{IN['inert_concentration_variation_kn1']:.1%}** between the centre and the surface,
while its flux is below 1e-10 of the reactant's (here {IN['inert_flux_ratio']:.0e}),
zero to the solve's tolerance; that zero is conservation, which any model gives. Under bulk
control the effect is bigger ({IN['inert_concentration_variation_kn0p1']:.1%} at
$\mathrm{{Kn}} = 0.1$) and under Knudsen control it nearly disappears
({IN['inert_concentration_variation_kn10']:.1%} at $\mathrm{{Kn}} = 10$): the inert is
dragged by momentum exchange with the reacting species, and under Knudsen control
there is none to exchange.

**A scalar effective-diffusivity closure with a finite positive $D_C$ predicts no
variation here, at every Knudsen number**, and that is run, not asserted:
replacing $[B^e]^{{-1}}$ by its diagonal and re-solving the same three pellets gives
a variation below 1e-12 (here {IN['inert_scalar_closure_variation']:.0e}). Not a
small error: a structurally different answer.
"""))

# %% [markdown]
# ### 2. Uphill diffusion, and what pore size controls
#
# The review says Kaza and Jackson (1980) showed uphill diffusion is possible
# inside a catalyst particle, "a phenomena impossible to explain with say the
# Bosanquet formula (110)". It gives no condition for it. What follows is not
# Kaza and Jackson's example, but it makes the statement quantitative.
#
# Uphill diffusion is a *bulk* phenomenon: it needs one species' flux to be driven
# by another's, which is exactly the off-diagonal of $[B^e]$. Knudsen friction is
# diagonal. So the coupling weakens as the pore shrinks, and one might expect a
# pore size below which the effect dies.
#
# **It does not work like that.** Uphill diffusion needs two things at once: a
# driving species whose flux drags the second one, and a second species close
# enough to *balanced* that the drag beats its own gradient. The first weakens
# with pore size; the second is a condition on composition. What pore size
# controls is therefore not *whether* uphill diffusion happens but *how finely
# tuned the composition has to be*. Fixing the composition and solving for a
# critical pore diameter answers a question whose answer is whatever composition
# was chosen, so the composition is an axis below and not a constant.
#
# Every edge and threshold is a root found by `find_roots` on the whole physical
# range (near-end N₂ mole fraction 0 to 1, pore diameter 0.1 nm to 0.1 m), which
# raises unless there is exactly one.

# %%
D0_BULK = 1e-3                                  # 1 mm: bulk-controlled


def n2_mid(d0, x_n2, n_cell=100, **kw):
    """Fluxes and face gradients on the middle face, near end (0, x_N2, 1 - x_N2)."""
    m, c = capillary(d0, [0.0, x_n2, 1.0 - x_n2], XB_CAP, n_cell=n_cell, **kw)
    _, g = m.faces(c)
    k = n_cell // 2
    return m.flux(c)[k], g[k]


def one_root(f, lo, hi, what, log=False, n_scan=33):
    """The single root of f on [lo, hi]; raises for none, several or failed points."""
    r = find_roots(f, lo, hi, n_scan=n_scan, log=log)
    if len(r["roots"]) != 1 or r["failed"]:
        raise RuntimeError(f"{what}: {r['message']}")
    return r["roots"][0]


def window_1mm_study(**kw):
    """Edges of the uphill window at 1 mm: N_N2 = 0 (lower) and grad c_N2 = 0 (upper)."""
    lo = one_root(lambda a: n2_mid(D0_BULK, a, **kw)[0][1], 0.0, 1.0, "N_N2 = 0 at 1 mm")
    hi = one_root(lambda a: n2_mid(D0_BULK, a, **kw)[1][1], 0.0, 1.0, "grad c_N2 = 0 at 1 mm")
    return {"uphill_window_lower_1mm": lo, "uphill_window_upper_1mm": hi}


UW = window_1mm_study()
print(f"pore diameter {D0_BULK*1e3:g} mm (bulk-controlled), far end x_N2 = 0.5")
print(f"{'x_N2 at near end':>17} {'N_N2':>12} {'grad c_N2':>12}   direction")
for a in [0.515, 0.505, 0.500, 0.495, 0.490, 0.485, 0.480]:
    N, g = n2_mid(D0_BULK, a)
    print(f"{a:17.3f} {N[1]:12.4e} {g[1]:12.4e}   {'uphill' if N[1] * g[1] > 0 else 'down-gradient'}")
print(f"\nuphill window at 1 mm: x_N2(near end) in [{UW['uphill_window_lower_1mm']:.5f}, "
      f"{UW['uphill_window_upper_1mm']:.5f}]")

# %%
X_SCAN = {"0p483": 0.4830, "0p485": 0.4850, "0p49": 0.4900, "0p495": 0.4950,
          "0p498": 0.4980, "0p4995": 0.4995, "0p4999": 0.4999}


def kn_of(d0, i=1, kn_den=KN_DEN, mass_exp=KN_EXP):
    """Kn = D^e_ij / D^e_iM with D(H2-N2) as D_ij, the abscissa of A4.4; i = 1 is N2."""
    return EPS_TAU * D_PAIR[(0, 1)] / knudsen_D(d0, T_REV, MOLAR[i], EPS_TAU, kn_den, mass_exp)


def d0_star(x_n2, n_cell=100, **kw):
    """Pore diameter at which N_N2 changes sign, at a given near-end composition."""
    return 10 ** one_root(lambda l: n2_mid(10 ** l, x_n2, n_cell=n_cell, **kw)[0][1],
                          -10.0, -1.0, f"N_N2 = 0 at x_N2 = {x_n2}")


def threshold_study(**kw):
    """The 'critical pore size' as a function of the near-end composition."""
    d = {tag: d0_star(x, **kw) for tag, x in X_SCAN.items()}
    out = {f"uphill_threshold_nm_at_x_{t}": d[t] * 1e9 for t in ("0p483", "0p49", "0p4999")}
    out["uphill_threshold_range_over_window"] = max(d.values()) / min(d.values())
    out["kn_star_at_x_0p49"] = float(kn_of(d["0p49"], kn_den=kw.get("kn_den", KN_DEN)))
    out["_d"] = d
    return out


EQ_THR = 1e-6     # N2 flux over H2 flux counted as nonzero; Newton noise is far below


def equal_end_study(**kw):
    """Equal N2 mole fraction at both ends: does N2 still carry a flux?"""
    rows = []
    for d0 in np.geomspace(1e-9, 1e-3, 25):
        N, g = n2_mid(d0, 0.5, **kw)
        rows.append((d0, abs(N[1] / N[0])))
    ratios = np.array([r[1] for r in rows])
    return {"uphill_at_equal_end_compositions": int(np.sum(ratios > EQ_THR)),
            "_rows": rows, "_ratios": ratios}


TH = threshold_study()
print("near-end x_N2 (far end 0.5), threshold pore diameter, and Kn there:")
print(f"{'x_N2':>9} {'d_0*':>16} {'Kn*':>12}")
for tag, x in X_SCAN.items():
    print(f"{x:9.4f} {TH['_d'][tag]*1e9:13.4f} nm {kn_of(TH['_d'][tag]):12.4g}")
print(f"across x_N2 = {min(X_SCAN.values())} to {max(X_SCAN.values())}, all inside the "
      f"uphill window, the threshold moves by a factor {TH['uphill_threshold_range_over_window']:.0f}")

EQ = equal_end_study()
print("\nnear end at x_N2 = 0.5 exactly (the same as the far end):")
print(f"{'d_0':>12} {'|N_N2/N_H2|':>13}")
for d0, r in EQ["_rows"][::4]:
    print(f"{d0*1e9:11.1f}n {r:13.2e}")

# %%
# What pore size does control: the width of the composition window. Its lower
# edge is N_N2 = 0; the upper edge (grad c_N2 = 0) sits within 3e-4 of the far-end
# composition 0.5, so the window is reported by its half-width below 0.5.
WIN_D0 = 10.0 ** np.arange(-9.0, -2.9, 0.5)      # 1 nm to 1 mm, half decades
SMALL = WIN_D0 <= 1e-8


def window_study(**kw):
    """Lower window edge at each pore diameter, and the half-width 0.5 - x*."""
    lo = np.array([one_root(lambda a: n2_mid(d, a, **kw)[0][1], 0.0, 1.0,
                            f"N_N2 = 0 at d0 = {d:.1e}") for d in WIN_D0])
    w = 0.5 - lo
    return {"uphill_window_halfwidth_1nm": w[0],
            "uphill_window_halfwidth_1um": w[int(np.flatnonzero(WIN_D0 == 1e-6)[0])],
            "uphill_window_halfwidth_1mm": w[-1],
            "uphill_window_width_slope_in_d0":
                float(np.polyfit(np.log(WIN_D0[SMALL]), np.log(w[SMALL]), 1)[0]),
            "_lo": lo, "_w": w}


def grid_spread_study(ns=(50, 100, 200, 400), **kw):
    """Largest relative change of the threshold over grids of 50 to 400 cells."""
    spread = []
    for x in (0.4900, 0.4980):
        v = [d0_star(x, n_cell=nc, **kw) for nc in ns]
        spread.append(max(v) / min(v) - 1.0)
    return {"uphill_threshold_grid_spread": max(spread)}


WN = window_study()
GS = grid_spread_study()
print(f"{'d_0':>12} {'Kn':>10} {'lower edge x_N2':>17} {'window half-width':>19}")
for d, lo, w in zip(WIN_D0, WN["_lo"], WN["_w"]):
    print(f"{d*1e9:11.2f}n {kn_of(d):10.4g} {lo:17.9f} {w:19.4e}")
print(f"\nlog-log slope of the half-width against d_0 up to 10 nm: "
      f"{WN['uphill_window_width_slope_in_d0']:.3f}")
print(f"threshold over grids of 50 to 400 cells moves by at most {GS['uphill_threshold_grid_spread']:.1e}")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
dd = np.geomspace(1e-9, 1e-3, 22)
for a, col in zip([0.4850, 0.4900, 0.4950, 0.4999, 0.5000], ["C0", "C1", "C2", "C4", "C3"]):
    NN = np.array([n2_mid(d, a)[0][1] for d in dd])
    axes[0].plot(dd * 1e9, NN, lw=2, color=col, label=rf"$x_{{\rm N_2}} = {a:g}$")
axes[0].axhline(0.0, color="0.5", lw=1)
axes[0].set_xscale("log"); axes[0].set_yscale("symlog", linthresh=1e-9)
axes[0].set_xlabel(r"pore diameter $d_0$ / nm")
axes[0].set_ylabel(r"$N_{\rm N_2}$ / mol m$^{-2}$ s$^{-1}$")
axes[0].set_title("Where N$_2$ turns round depends on the composition\n"
                  r"(positive = uphill, since $\nabla c_{\rm N_2} > 0$)")
axes[0].legend(fontsize=8)
axes[1].loglog(WIN_D0 * 1e9, WN["_w"], "o-", lw=2, color="C0", label="computed")
axes[1].loglog(WIN_D0 * 1e9, WN["_w"][0] * (WIN_D0 / WIN_D0[0]), ":", color="0.5",
               label=r"slope 1 ($\propto d_0$)")
axes[1].set_xlabel(r"pore diameter $d_0$ / nm")
axes[1].set_ylabel(r"$0.5 - x_{\rm N_2}^*$, half-width of the uphill window")
axes[1].set_title("The window narrows with the pore, but never closes")
axes[1].legend(fontsize=8)
plt.tight_layout(); plt.show()

# %%
display(Markdown(rf"""
With the far end held at (H₂, N₂, CO₂) = (0.5, 0.5, 0) and the near end at
(0, 0.49, 0.51), nitrogen's concentration **rises** towards the far end and its flux
points that way too: it diffuses up its own gradient, dragged by the much larger
hydrogen flux going the other way. At a 1 mm pore the window in which it does so
is narrow, $x_{{\rm N_2}}$ between {UW['uphill_window_lower_1mm']:.4f} and
{UW['uphill_window_upper_1mm']:.4f} at the near end, which is why the phenomenon is
easy to miss.

**There is no single critical pore size.** The pore diameter at which
$N_{{\rm N_2}}$ changes sign is {TH['uphill_threshold_nm_at_x_0p49']:.1f} nm
($\mathrm{{Kn}}^* = {TH['kn_star_at_x_0p49']:.2f}$) at $x_{{\rm N_2}} = 0.49$, but
{TH['uphill_threshold_nm_at_x_0p483']:.0f} nm at 0.483 and
{TH['uphill_threshold_nm_at_x_0p4999']:.3f} nm at 0.4999, a factor
**{TH['uphill_threshold_range_over_window']:.0f}** across compositions that all sit
inside the window above. With the near end at $x_{{\rm N_2}} = 0.5$ exactly, N₂
still carries a flux at each of {EQ['uphill_at_equal_end_compositions']} pore
diameters from 1 nm to 1 mm (a quarter decade apart), at least
{EQ['_ratios'].min():.1e} of hydrogen's, while its mole fraction is the same at both
ends: its gradient averages to zero along the capillary, so over part of the
capillary N₂ runs up its own gradient at every one of those pore sizes. A threshold
quoted alone is a property of the composition it was
computed at. Over grids of 50 to 400 cells the threshold moves by at most
{GS['uphill_threshold_grid_spread']:.0e}: that measures the solver, not the
threshold, and it is the composition, not the discretisation, that decides what a
threshold is worth.

**What pore size does control** is the width of the composition window. Its lower
edge $x^*_{{\rm N_2}}(d_0)$ climbs towards the far-end composition as the pore
shrinks, and the half-width $0.5 - x^*$ falls from
{WN['uphill_window_halfwidth_1mm']:.2e} at 1 mm to
{WN['uphill_window_halfwidth_1um']:.2e} at 1 µm and
{WN['uphill_window_halfwidth_1nm']:.2e} at 1 nm, a log–log slope of
**{WN['uphill_window_width_slope_in_d0']:.2f}** up to 10 nm: nearly proportional to
$d_0$, as expected when the bulk coupling that drives the effect is set against a
Knudsen resistance that is not. It does not reach zero at any pore size computed.

**The screening rule.** Uphill diffusion is not a macropore phenomenon that stops
below some diameter. Its composition window shrinks in proportion to the pore.
Rare in a nanoporous catalyst, yes, but the dusty gas model does not forbid it,
and no scalar closure shows it at any pore size.
"""))

# %% [markdown]
# ### 3. Internal pressure, and what Jackson's $\sqrt{\nu_B}$ needs
#
# The review quotes Jackson (1977): for $A \to \nu_B B$ with complete conversion,
# the pellet-centre pressure is $p_0 = \sqrt{\nu_B}\,p$, "for $\nu_B = 2$ … a 40%
# increase". Page `A4.4` reproduces that number for a binary pellet. Running the
# $n$-component model over the conditions shows it is the corner of a wider result.
#
# In the Knudsen limit each species obeys $N_i = -\mathrm{D}^e_{iM}\nabla c_i$
# with $N_B = -\nu_B N_A$, so
# $\nabla c_B = -\nu_B(\mathrm{D}^e_{AM}/\mathrm{D}^e_{BM})\nabla c_A
# = -\nu_B\sqrt{M_B/M_A}\,\nabla c_A$ by eq. (85), and integrating from the surface
# to a fully converted centre gives
#
# $$
# \frac{p_0}{p} = 1 + \left(\nu_B\sqrt{M_B/M_A} - 1\right) x_{A,s} .
# $$
#
# That reduces to $\sqrt{\nu_B}$ only when **all three** of the following hold:
# $x_{A,s} = 1$ (pure $A$ at the surface), $M_A = \nu_B M_B$ (mass conserved by
# the stoichiometry, which is automatic for a real reaction) and Knudsen control.
# The third is the one that is easy to forget. The departure at a finite Knudsen
# number is extrapolated to zero cell size (200, 400 and 800 cells), so that the
# discretisation is removed before it is compared across $\mathrm{Kn}$ and $\phi$.

# %%
PHI_HI, KN_HI = 60.0, 1e4
KN_SWEEP = [1e6, 1e5, 1e4, 1e2, 1e1, 1.0, 1e-1, 1e-2]


def p0_of(*args, **kw):
    """Pellet-centre total pressure over the surface pressure."""
    return float(centre(*pellet(*args, **kw)).sum())


def jackson_study(phis=(PHI_HI, 200.0), **kw):
    """Centre pressure against sqrt(nu_B) and against the Knudsen-limit closed form."""
    jack = [(nu, p0_of(nu, KN_HI, PHI_HI, **kw)) for nu in (1.0, 2.0, 3.0, 4.0)]
    surf = [(x, p0_of(2.0, KN_HI, PHI_HI, **{**kw, "x_as": x}), 1 + (np.sqrt(2.0) - 1) * x)
            for x in (1.0, 0.75, 0.5, 0.25)]
    mass = [(mr, p0_of(2.0, KN_HI, PHI_HI, **{**kw, "mass_ratio": mr}), 2.0 * np.sqrt(mr))
            for mr in (0.5, 1.0, 2.0)]
    knsw = [(kn, p0_of(2.0, kn, PHI_HI, **kw)) for kn in KN_SWEEP]

    def dev_ext(kn, phi):    # departure from sqrt(2), extrapolated to zero cell size
        return observed_orders(lambda n: p0_of(2.0, kn, phi, **{**kw, "n_cell": n})
                               / np.sqrt(2.0) - 1.0, (200, 400, 800))["extrapolated"]

    d4, d6, d4phi = dev_ext(1e4, phis[0]), dev_ext(1e6, phis[0]), dev_ext(1e4, phis[1])
    return {"jackson_sqrt_nu_worst": max(abs(p / np.sqrt(nu) - 1) for nu, p in jack),
            "knudsen_closed_form_worst": max(abs(p / ref - 1) for _, p, ref in surf + mass),
            "jackson_p0_nu2": jack[1][1], "jackson_p0_half_inert": surf[2][1],
            "jackson_p0_equal_mass": mass[1][1], "jackson_p0_bulk_kn0p01": knsw[-1][1],
            "jackson_dev_times_kn": d4 * 1e4,
            "jackson_kn_rate": float(np.log(abs(d6 / d4)) / np.log(1e6 / 1e4)),
            "jackson_phi_effect": abs(d4phi / d4 - 1.0),
            "_jack": jack, "_surf": surf, "_mass": mass, "_knsw": knsw,
            "_d": (d4, d6, d4phi)}


JK = jackson_study()
print("Knudsen control (Kn = 1e4, phi = 60, 200 cells), pure A at the surface, mass conserved:")
for nu, p in JK["_jack"]:
    print(f"  nu_B = {nu:.0f}   p_0/p = {p:.6f}   sqrt(nu_B) = {np.sqrt(nu):.6f}   dev {p/np.sqrt(nu)-1:+.2e}")
print("surface not pure A  ->  1 + (sqrt(nu_B) - 1) x_A,s:")
for x, p, ref in JK["_surf"]:
    print(f"  x_A,s = {x:.2f}   p_0/p = {p:.6f}   closed form {ref:.6f}   dev {p/ref-1:+.2e}")
print("mass ratio free  ->  1 + (nu_B sqrt(M_B/M_A) - 1):")
for mr, p, ref in JK["_mass"]:
    print(f"  M_B/M_A = {mr:.1f}   p_0/p = {p:.6f}   closed form {ref:.6f}   dev {p/ref-1:+.2e}")
print("bulk friction spoils it:")
for kn, p in JK["_knsw"]:
    print(f"  Kn = {kn:8.4g}   p_0/p = {p:.6f}   dev from sqrt(2) {p/np.sqrt(2)-1:+.3e}")
d4, d6, d4phi = JK["_d"]
print(f"extrapolated departure from sqrt(2): {d4:+.4e} at Kn = 1e4, {d6:+.4e} at Kn = 1e6, "
      f"{d4phi:+.4e} at Kn = 1e4 and phi = 200")

# %%
display(Markdown(rf"""
Under Knudsen control, with pure $A$ at the surface and mass-conserving
stoichiometry, the solver returns $p_0/p = {JK['jackson_p0_nu2']:.5f}$ against Jackson's
$\sqrt{{2}}$; the worst departure over $\nu_B = 1, 2, 3, 4$ is
{JK['jackson_sqrt_nu_worst']:.1e}. Extrapolated to zero cell size that departure is
the finite Knudsen number of the run: $\mathrm{{Kn}}$ times the departure is
{JK['jackson_dev_times_kn']:.4f}, so the centre pressure falls short of
$\sqrt{{2}}$ by $1/(2\,\mathrm{{Kn}})$, and the departure falls as
$\mathrm{{Kn}}^{{{JK['jackson_kn_rate']:.3f}}}$ between $10^4$ and $10^6$. Pushing $\phi$
from 60 to 200 changes the extrapolated departure by
{JK['jackson_phi_effect']:.0e} of itself: conversion is complete at $\phi = 60$.

Relax any one of the three conditions and the number moves a long way. With half
the surface inert the rise is {JK['jackson_p0_half_inert']:.4f}. With $M_B = M_A$
instead of $M_A/2$ (the same $\nu_B = 2$, but mass not conserved) it is
{JK['jackson_p0_equal_mass']:.4f}. And under bulk control the pressure build-up
nearly vanishes: at $\mathrm{{Kn}} = 0.01$ the centre is only
{JK['jackson_p0_bulk_kn0p01']:.4f} times the surface pressure. The closed form above
reproduces the surface and mass-ratio cases to {JK['knudsen_closed_form_worst']:.1e}
at $\mathrm{{Kn}} = 10^4$.

**So $\sqrt{{\nu_B}}$ is a free-molecule result**, and quoting it for a
bulk-controlled pellet overstates the internal pressure by a factor approaching
$\sqrt{{\nu_B}}$ itself.
"""))

# %% [markdown]
# ### 4. What the viscous term is worth
#
# The review, on p. 893: "neglect of the viscous flow contribution … is not very
# serious; this result is typical (Haynes, 1978)". It does not say how serious.
#
# The viscous term enters eq. (101) as
# $\alpha'_i x_i (\mathrm{D}_{\rm visc}/\mathrm{D}^e_{iM})\nabla p/RT$, so the
# group that decides everything is $\mathrm{D}_{\rm visc}/\mathrm{D}^e_{iM}$.
# Substituting eq. (106) for $\mathrm{D}_{\rm visc}$, eq. (91) for $B_0$ and
# eq. (85) for $\mathrm{D}^e_{iM}$ gives a closed form in printed quantities only:
#
# $$
# \frac{\mathrm{D}_{\rm visc}}{\mathrm{D}^e_{iM}}
# = \frac{c_t B_0 RT/\eta}{(\varepsilon/\tau)(d_0/3)\bar v_i}
# = \frac{3\,p\,d_0}{32\,\eta\,\bar v_i}\cdot\frac{\tau}{\varepsilon},
# \qquad \bar v_i = \sqrt{\frac{8RT}{\pi M_i}} .
# $$
#
# Linear in pore size and in pressure. That is the screening rule, and it makes
# the review's remark checkable.
#
# **But it cannot be swept on its own.** The pellet's other dimensionless group,
# $\mathrm{Kn} = \mathrm{D}^e_{ij}/\mathrm{D}^e_{iM}$, is *also* a function of
# $d_0$, inversely, because eq. (85) puts $d_0$ in $\mathrm{D}^e_{iM}$ and
# eq. (91) puts $d_0^2$ in $B_0$. Holding Kn fixed while dialling
# $\mathrm{D}_{\rm visc}/\mathrm{D}^e_{iM}$ therefore describes no pore that
# exists: the two groups are one variable. The sweep below moves $d_0$ and lets
# both follow, and the cost of not doing that is measured at the end.

# %%
ETA_GAS = 2.0e-5                # Pa s, a representative gas viscosity, not from the
                                # review; it enters only this group


def d_visc_over_diM(d0, pois_den=POIS_DEN, kn_den=KN_DEN, p=P_ATM, T=T_REV, eta=ETA_GAS):
    """eq. (106) over eq. (85), with eq. (91) for B_0; one value per species."""
    ct = p / (RGAS * T)
    B0 = d0 ** 2 / pois_den
    return (ct * B0 * RGAS * T / eta) / knudsen_D(d0, T, MOLAR, EPS_TAU, kn_den, KN_EXP)


def d_visc_closed(d0, p=P_ATM, T=T_REV, eta=ETA_GAS):
    """The same group written out: 3 p d_0 / (32 eta v_bar_i) * tau/eps."""
    vbar = np.sqrt(8 * RGAS * T / (np.pi * MOLAR))
    return (KN_DEN * p * d0) / (POIS_DEN * eta * vbar) / EPS_TAU


def visc_group_study(pois_den=POIS_DEN, kn_den=KN_DEN):
    """The assembled screening group against its closed form, and Kn at two pores."""
    scan = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
    dev = max(float(np.max(np.abs(d_visc_over_diM(d, pois_den, kn_den) / d_visc_closed(d) - 1)))
              for d in scan)
    return {"visc_group_closed_form_dev": dev,
            "visc_group_1nm": float(d_visc_over_diM(1e-9, pois_den, kn_den)[0]),
            "visc_group_1um": float(d_visc_over_diM(1e-6, pois_den, kn_den)[0]),
            "kn_1nm": float(kn_of(1e-9, 0, kn_den)), "kn_1um": float(kn_of(1e-6, 0, kn_den))}


VG = visc_group_study()
print(f"{'d_0':>10}  " + "  ".join(f"D_visc/D^e_iM ({s})" for s in SPEC))
for d in [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]:
    print(f"{d*1e9:9.1f}nm  " + "  ".join(f"{x:18.4e}" for x in d_visc_over_diM(d)))
print(f"\nclosed form 3 p d_0/(32 eta v_bar) reproduces the assembled group to "
      f"{VG['visc_group_closed_form_dev']:.1e}")

# %%
# One variable, d_0: Kn and D_visc/D^e_AM both follow it (A = H2, B its product).
VISC_D0 = 10.0 ** np.arange(-9.0, -3.9, 0.25)     # 1 nm to 100 um, quarter decades
LADDER = (200, 400, 800)


def relief(d0, n_cell=200, pois_den=POIS_DEN, kn_den=KN_DEN, frozen_kn=None, **kw):
    """Relative change of p_0/p when the viscous term is kept, A -> 2B, phi = 60."""
    kn = float(kn_of(d0, 0, kn_den)) if frozen_kn is None else frozen_kn
    dv = float(d_visc_over_diM(d0, pois_den, kn_den)[0])
    a = p0_of(2.0, kn, PHI_HI, n_cell=n_cell, **kw)
    b = p0_of(2.0, kn, PHI_HI, n_cell=n_cell, dvisc=dv, **kw)
    return b / a - 1.0, kn, dv, a, b


def relief_ext(d0, **kw):
    """relief extrapolated to zero cell size from 200, 400 and 800 cells."""
    return observed_orders(lambda n: relief(d0, n_cell=n, **kw)[0], LADDER)["extrapolated"]


def viscous_study(**kw):
    """Viscous relief at named pores, at its maximum over pore size, and at frozen Kn."""
    scan = [relief(d, **kw) for d in VISC_D0]
    k = int(np.argmax([abs(s[0]) for s in scan]))       # brackets the maximum
    if k in (0, len(VISC_D0) - 1):
        raise RuntimeError("largest relief at the end of the scan: the maximum is at or beyond it")
    opt = minimize_scalar(lambda l: -abs(relief_ext(10 ** l, **kw)),
                          bounds=(np.log10(VISC_D0[k - 1]), np.log10(VISC_D0[k + 1])),
                          method="bounded", options={"xatol": 1e-4})
    d_max = 10 ** opt.x
    one_um = relief_ext(1e-6, **kw)
    frozen = relief_ext(1e-6, frozen_kn=KN_HI, **kw)
    kn_den = kw.get("kn_den", KN_DEN)
    return {"viscous_relief_1nm_pore": abs(relief_ext(1e-9, **kw)),
            "viscous_relief_1um_pore": abs(one_um),
            "viscous_relief_worst_over_pore_sizes": abs(relief_ext(d_max, **kw)),
            "viscous_relief_worst_at_kn": float(kn_of(d_max, 0, kn_den)),
            "viscous_relief_worst_d0_nm": d_max * 1e9,
            "viscous_relief_fixed_kn_artefact": abs(frozen),
            "viscous_fixed_kn_overstatement": frozen / one_um,
            "_scan": scan}


VS = viscous_study()
print("A -> 2B, pure A at the surface, Kn and the viscous group both set by d_0 (200 cells):")
print(f"{'d_0':>11} {'Kn':>10} {'D_visc/D^e_AM':>14} {'p_0/p no visc':>14} "
      f"{'p_0/p visc':>12} {'change':>10}")
for d0, (ch, kn, dv, a, b) in zip(VISC_D0, VS["_scan"]):
    print(f"{d0*1e9:10.2f}n {kn:10.4g} {dv:14.4g} {a:14.6f} {b:12.6f} {ch:+10.4%}")
print(f"\nextrapolated to zero cell size: {VS['viscous_relief_1nm_pore']:.3%} at 1 nm, "
      f"{VS['viscous_relief_1um_pore']:.3%} at 1 um; largest "
      f"{VS['viscous_relief_worst_over_pore_sizes']:.3%} at d_0 = "
      f"{VS['viscous_relief_worst_d0_nm']:.0f} nm (Kn = {VS['viscous_relief_worst_at_kn']:.3f})")
print(f"Kn frozen at {KN_HI:g} with the 1 um viscous group: "
      f"{VS['viscous_relief_fixed_kn_artefact']:.2%}")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
d0n = VISC_D0 * 1e9
axes[0].semilogx(d0n, [s[3] for s in VS["_scan"]], "o-", lw=2, color="C0", ms=4,
                 label="viscous term off")
axes[0].semilogx(d0n, [s[4] for s in VS["_scan"]], "s-", lw=2, color="C3", ms=4,
                 label="viscous term on")
axes[0].axhline(np.sqrt(2.0), color="0.5", ls=":", lw=1)
axes[0].annotate(r"$\sqrt{\nu_B}$", (d0n[0], np.sqrt(2.0)), fontsize=9, va="bottom", color="0.4")
axes[0].set_xlabel(r"pore diameter $d_0$ / nm"); axes[0].set_ylabel(r"$p_0/p$")
axes[0].set_title("There is barely any pressure rise left to bleed off\n"
                  "once the pore is big enough for viscous flow to matter")
axes[0].legend(fontsize=8)
axes[1].semilogx(d0n, [100 * abs(s[0]) for s in VS["_scan"]], "o-", lw=2, color="C0", ms=4)
axes[1].axvline(VS["viscous_relief_worst_d0_nm"], color="C3", ls="--", lw=1.5)
axes[1].set_xlabel(r"pore diameter $d_0$ / nm")
axes[1].set_ylabel(r"change in $p_0/p$ when viscous flow is kept / %")
axes[1].set_title("Viscous flow bleeds off under a percent, at every pore size\n"
                  r"($A\rightarrow 2B$, H$_2$/$B$ at 1 atm, $\phi = 60$)")
plt.tight_layout(); plt.show()

# %%
display(Markdown(rf"""
For hydrogen at 1 atm the screening group is
$\mathrm{{D}}_{{\rm visc}}/\mathrm{{D}}^e_{{AM}} = {VG['visc_group_1nm']:.3g}$ in a 1 nm
pore and ${VG['visc_group_1um']:.3g}$ in a 1 µm one, exactly the pore-size ratio,
because the group is linear in $d_0$. The Knudsen number runs the other way over the
same range, {VG['kn_1nm']:.1f} down to {VG['kn_1um']:.3f}, and that is what decides
the answer.

Sweeping $d_0$ and letting both follow, viscous flow changes the pellet-centre
pressure by **{VS['viscous_relief_1nm_pore']:.4%}** at 1 nm and
**{VS['viscous_relief_1um_pore']:.2%}** at 1 µm (both extrapolated to zero cell
size). Over the whole range 1 nm to 100 µm the change is largest at
{VS['viscous_relief_worst_d0_nm']:.0f} nm, where $\mathrm{{Kn}} =
{VS['viscous_relief_worst_at_kn']:.2f}$, and is **{VS['viscous_relief_worst_over_pore_sizes']:.2%}**
there.

**So "not very serious" holds at every pore size, not just for catalysts**, the
opposite of what a fixed-Kn sweep suggests. The reason is visible in the left
panel: the reaction-generated pressure rise that viscous flow could bleed off is
itself a Knudsen effect. In a small pore there is a large rise but almost no
viscous conductance; in a large pore there is plenty of viscous conductance but
almost nothing to relieve. The two factors multiply to a bounded, non-monotonic
curve.

Freezing $\mathrm{{Kn}}$ at $10^4$, a free-molecule pellet, while dialling the
viscous group to its 1 µm value gives {VS['viscous_relief_fixed_kn_artefact']:.1%},
a factor {VS['viscous_fixed_kn_overstatement']:.1f} more than a 1 µm pore shows: a
percentage of a pressure rise no 1 µm pore has. What this does *not* license is
dropping the term where the pressure gradient is imposed externally rather than
generated by reaction: eq. (103) is that case, and there viscous flow carries most
of the flux at large $d_0$.
"""))

# %% [markdown]
# ## Validation
#
# Ranked by how much each one can catch, strongest first. Every number below is a
# claim in the table at the end, with a **defect row**: a defect is injected that
# the check ought to see, and the claim has to fail.
#
# 1. **`V1`, an independent collocation solution of the same equations**, and the
#    grid convergence of the finite-volume route towards it. The only check here
#    that shares no *transport* code with the solver (it shares the two helpers that
#    evaluate eqs. 83 and 85 into numbers).
# 2. **`V2`, eq. (103), a printed closed form the code does not contain.** The one
#    check with power over the viscous term.
# 3. **`V3`, Graham's law, predicted rather than imposed.** An algebraic identity of
#    eq. (82); labelled as one.
# 4. **`V4`, the two asymptotes**, free-molecule and Maxwell–Stefan, with their rates.
# 5. **`V5`, the collapse to `A4.4`'s scalar result.**
# 6. **`V6`, the four exponents the review states on p. 887**, as log–log slopes.
#
# ### V1: an independent collocation solution
#
# The same physics is solved a second time by a different numerical method: eq. (82)
# is written out as an ODE system in the three concentrations with the three
# *constant* fluxes as unknown parameters, and handed to `scipy.integrate.solve_bvp`
# (`capillary_bvp` in `model.py`, which raises unless it converges).
#
# What they do **not** share: any line of the transport code. The collocation route
# never calls `b_matrix`, `faces` or `flux`, never touches `construct_grad`,
# `construct_div`, `NumJac` or `newton`, has no finite-volume grid, and solves for a
# different set of unknowns, the three fluxes, against the concentration field.
# Uniform pressure and Graham's law are imposed in neither.
#
# What they *do* share: `knudsen_D` and the pair products, which turn eqs. (83) and
# (85) into numbers, and the molar masses. A mistake inside those moves both routes
# together. Eq. (85)'s mass exponent is pinned separately by Graham's law (`V3`);
# its $d_0/3$ prefactor and eq. (83)'s $\varepsilon/\tau$ are pinned by nothing on
# this page, which is measured below ("What nothing here pins").
#
# The finite-volume flux on the middle face is refined from 25 to 800 cells; the
# order comes from successive differences (`observed_orders`), the extrapolated flux
# is compared with the collocation route, and the finest difference is compared with
# the Newton tolerance (the round-off companion).

# %%
V1_D0 = [1e-3, 1e-5, 1e-7, 1e-8, 1e-9]
NS_GRID = (25, 50, 100, 200, 400, 800)


def collocation_study(tol=1e-12, **kw):
    """Finite-volume route against collocation, and its grid ladder at 100 nm."""
    devs = [float(np.max(np.abs(mid_flux(d0, tol=tol, **kw) / capillary_bvp(d0, X_UP, XB_CAP) - 1)))
            for d0 in V1_D0]
    st = observed_orders(lambda n: mid_flux(1e-7, n_cell=n, tol=tol, **kw), NS_GRID)
    nb = capillary_bvp(1e-7, X_UP, XB_CAP)
    orders = np.asarray(st["orders"], float)
    vals = np.asarray(st["values"], float)
    finest = float(np.min(np.abs(vals[-1] - vals[-2]) / np.abs(vals[-1])))
    return {"collocation_worst_dev": max(devs),
            "collocation_extrapolated_dev": float(np.max(np.abs(st["extrapolated"] / nb - 1))),
            "grid_order": float(orders[np.argmax(np.abs(orders - 2.0))]),
            "grid_ladder_over_roundoff": finest / tol,
            "grid_error_n100": float(np.max(np.abs(vals[2] / nb - 1))),
            "_devs": devs, "_orders": orders, "_vals": vals}


V1 = collocation_study()
print(f"{'d_0':>12} {'finite volume (100 cells) vs collocation':>42}")
for d0, dev in zip(V1_D0, V1["_devs"]):
    print(f"{d0*1e9:11.1f}nm {dev:42.1e}")
print(f"\ngrid ladder at 100 nm: cells {NS_GRID}")
print(f"  orders from successive differences: {np.round(V1['_orders'], 4)}")
print(f"  extrapolated flux vs collocation within 1e-9: {V1['collocation_extrapolated_dev'] < 1e-9}; "
      f"100-cell error {V1['grid_error_n100']:.1e}; finest difference over the Newton "
      f"tolerance {V1['grid_ladder_over_roundoff']:.1e}")

# %%
display(Markdown(rf"""
The two routes agree to {V1['collocation_worst_dev']:.1e} at worst over five pore
sizes spanning six decades on the 100-cell grid, and the finite-volume flux converges
at order {V1['grid_order']:.3f} (the order farthest from 2 on the ladder). Extrapolated
to zero cell size it meets the collocation route to within 1e-9 (here
{V1['collocation_extrapolated_dev']:.0e}), the collocation tolerance: the remaining
deviation is the finite-volume discretisation error, not a disagreement about the
physics. The defect rows show what moves it: eq. (87)'s off-diagonals, the mass exponent
of eq. (85) in one route, the scalar closure, the face-composition rule; deleting
the Knudsen term leaves $[B^e]$ singular and the solve without an answer.
"""))

# %% [markdown]
# ### V2: eq. (103), which the code does not contain
#
# For a single component the review prints a closed form,
# $N_1 = -(\mathrm{D}^e_{1M}/RT + c_1 B_0/\eta)\nabla p$. The solver implements
# eq. (101), knows nothing about eq. (103), and has no special case for $n = 1$.
# Integrating the closed form across a slab with $c$ imposed at both ends gives
#
# $$
# N_1 L = \mathrm{D}^e_{1M}(c_A - c_B) + \frac{RTB_0}{\eta}\,\frac{c_A^2 - c_B^2}{2},
# $$
#
# with no free constant. With the arithmetic face mean, the discrete flux
# $-(\mathrm{D}^e_{1M} + (RTB_0/\eta)\,\bar c)\,\Delta c/h$ telescopes to the same
# expression, so the comparison is exact up to the solve's tolerance on any grid.

# %%
D0_103, L_103 = 1e-7, 1e-3
B0_103 = D0_103 ** 2 / POIS_DEN                       # eq. (91)
CA_103, CB_103 = CT0, 0.9 * CT0


def eq103_study(B0=B0_103, alpha=ALPHA_DEF, n_cell=200, kn_den_fv=KN_DEN,
                kn_den_ref=KN_DEN, **kw):
    """Single-component solver flux against eq. (103)'s closed form."""
    dm_fv = knudsen_D(D0_103, T_REV, MOLAR[:1], EPS_TAU, kn_den_fv, KN_EXP)
    dm_ref = knudsen_D(D0_103, T_REV, MOLAR[:1], EPS_TAU, kn_den_ref, KN_EXP)
    m = DustyGas(n_cell, L_103, np.zeros((1, 1)), dm_fv, T_REV, eps_tau=EPS_TAU,
                 bc=(dirichlet([CA_103]), dirichlet([CB_103])),
                 B0=B0, eta=ETA_GAS, alpha=[alpha], **kw)
    c = m.solve(np.linspace(CA_103, CB_103, n_cell).reshape(n_cell, 1))
    exact = (dm_ref[0] * (CA_103 - CB_103)
             + (RGAS * T_REV * B0_103 / ETA_GAS) * (CA_103 ** 2 - CB_103 ** 2) / 2) / L_103
    return {"eq103_viscous_dev": abs(m.flux(c)[n_cell // 2, 0] / exact - 1.0)}


V2 = eq103_study()
print(f"eq. (103) recovered by the general solver to {V2['eq103_viscous_dev']:.1e}")

# %% [markdown]
# ### V3: Graham's law, predicted and not imposed
#
# Nothing in `DustyGas` mentions eq. (107) or $\sqrt{M_i}$. The capillary is given
# equal total concentration at its two ends and left to find the fluxes; whether
# $\sum_i N_i\sqrt{M_i}$ vanishes, and whether the pressure stays uniform, are
# outputs. **This is an algebraic identity of eq. (82)**, derived symbolically
# above. What it tests is the implementation: eq. (85)'s mass exponent and
# eq. (87)'s antisymmetry.

# %%
def graham_study(**kw):
    """Graham residual and total-pressure spread over three pore sizes."""
    g, ps = [], []
    for d0 in (1e-3, 1e-6, 1e-8):
        m, c = capillary(d0, X_UP, XB_CAP, **kw)
        w = m.flux(c)[m.n_cell // 2] * MOLAR ** GR_EXP
        ct = m.faces(c)[0].sum(axis=1)
        g.append(float(abs(w.sum()) / np.abs(w).sum()))
        ps.append(float(np.ptp(ct) / ct.mean()))
    return {"graham_residual_worst": max(g), "pressure_uniformity": max(ps), "_g": g, "_p": ps}


V3 = graham_study()
print(f"at d_0 = 1 mm, 1 um, 10 nm: |sum N sqrt(M)| (relative) within 1e-9: "
      f"{V3['graham_residual_worst'] < 1e-9}; pressure spread within 1e-10: "
      f"{V3['pressure_uniformity'] < 1e-10}")

# %% [markdown]
# ### What nothing here pins: eq. (85)'s $d_0/3$
#
# Every check so far compares two routes, and a defect in one of them shows. Is
# there a defect that *every* route on this page would let through? Reading $d_0$
# as a radius rather than a diameter, $d_0/3 \to d_0/2$, injected into
# `knudsen_D` itself moves every route at once:
#
# - **`V1`** compares pymrm against collocation, and both call `knudsen_D`;
# - **`V2`**'s closed form $\mathrm{D}^e_{1M}\Delta c + (RTB_0/\eta)\Delta(c^2)/2$
#   calls `knudsen_D` too, so a wrong prefactor moves both sides together;
# - **`V3`** is $\sum_i N_i\sqrt{M_i} = 0$, out of which the prefactor **cancels
#   identically**: Graham's law pins the mass exponent and nothing else;
# - **`V4`**, **`V5`** and **`V6`** share eq. (85) with their references; `V6`
#   tests linearity in $d_0$, which $d_0/2$ preserves.

# %%
def prefactor_study(den_fv=2.0, den_ref=2.0, mass_exp=KN_EXP):
    """Eq. (85)'s prefactor broken (d0/2) in the finite-volume and reference routes."""
    base = float(np.max(np.abs(mid_flux(1e-7) / capillary_bvp(1e-7, X_UP, XB_CAP) - 1)))
    v1 = float(np.max(np.abs(mid_flux(1e-7, kn_den=den_fv, mass_exp=mass_exp)
                             / capillary_bvp(1e-7, X_UP, XB_CAP, kn_den=den_ref) - 1)))
    m, c = capillary(1e-6, X_UP, XB_CAP, kn_den=den_fv, mass_exp=mass_exp)
    w = m.flux(c)[m.n_cell // 2] * MOLAR ** GR_EXP
    return {"break_eq85_prefactor_v1": v1, "prefactor_v1_ratio": v1 / base,
            "break_eq85_prefactor_v2":
                eq103_study(kn_den_fv=den_fv, kn_den_ref=den_ref)["eq103_viscous_dev"],
            "break_eq85_prefactor_v3": float(abs(w.sum()) / np.abs(w).sum())}


PF = prefactor_study()
print(f"d_0/3 -> d_0/2 in every route: V1 {PF['break_eq85_prefactor_v1']:.2e} "
      f"(x{PF['prefactor_v1_ratio']:.1f} the printed-prefactor value), "
      f"V2 within 1e-10: {PF['break_eq85_prefactor_v2'] < 1e-10}, "
      f"V3 within 1e-9: {PF['break_eq85_prefactor_v3'] < 1e-9}")

# %%
display(Markdown(rf"""
Measured: $d_0/3 \to d_0/2$ in every route moves `V1` by a factor
{PF['prefactor_v1_ratio']:.1f}, which is the change of the 100-cell discretisation
error with the Knudsen number, not a disagreement, and leaves `V2` below 1e-10
(here {PF['break_eq85_prefactor_v2']:.0e}) and `V3` below 1e-9 (here
{PF['break_eq85_prefactor_v3']:.0e}), the solve's tolerance. **So eq. (85)'s $d_0/3$ is pinned by nothing on this page**,
and neither is eq. (83)'s $\varepsilon/\tau$. Had the transcription been wrong
there, every dimensional number on this page (the Knudsen numbers, the uphill
window, the viscous screening group) would be off by that factor while the
Validation section still passed. The defence is the transcription itself, with the
review's own sentence naming $d_0$ the pore *diameter*: a blind spot, stated as
one.
"""))

# %% [markdown]
# ### V4: the two asymptotes
#
# **Free-molecule.** As $d_0 \to 0$ the model must collapse to
# $N_i = -\mathrm{D}^e_{iM}\nabla c_i$, one species at a time, with a departure
# that falls in proportion to $d_0$. The reference is a one-line expression; it
# shares `knudsen_D` with the solver, so it tests the matrix assembly and the
# discretisation, not eq. (85).
#
# **Maxwell–Stefan.** As $d_0 \to \infty$ the wall term disappears and the model
# must collapse to Maxwell–Stefan with a departure that falls as $1/d_0$; but
# $[B^e]$ becomes **singular** in that limit, so it can only be approached. The
# reference is an independent collocation solve of the Maxwell–Stefan equations in
# mole fractions with Graham's law as the closure (`ms_bvp` in `model.py`). Near the
# limit the flux is sensitive to the total-pressure gradient, so the Newton
# tolerance on the concentrations enters the flux increasingly as d0 grows; its effect is measured by repeating the 1 mm point at a tenfold
# tighter tolerance. The departure is quoted at 1 mm, where that effect is small.
#
# Both use the finite-volume flux extrapolated to zero cell size.

# %%
FM_D0 = (1e-9, 1e-10, 1e-11)
MS_D0 = (1e-4, 1e-3, 1e-2)


def v4_study(tol=1e-12, ms_d0=MS_D0, **kw):
    """Departures from the free-molecule and Maxwell-Stefan limits, and their rates."""
    def ext(d0, ns, t):
        return observed_orders(lambda n: mid_flux(d0, n_cell=n, tol=t, **kw), ns)["extrapolated"]

    fm = []
    for d0 in FM_D0:
        n_knud = (knudsen_D(d0, T_REV, MOLAR, EPS_TAU, KN_DEN, KN_EXP)
                  * (CT0 * (X_UP - XB_CAP)) / L_CAP)
        fm.append(float(np.max(np.abs(ext(d0, (50, 100, 200), tol) / n_knud - 1))))
    n_ms = _ms_bvp(X_UP, XB_CAP, molar=MOLAR, dij_ct=REF["dij_ct"], T=T_REV, p=P_ATM,
                   L=L_CAP, eps_tau=EPS_TAU, gr_exp=GR_EXP)
    ms = [float(np.max(np.abs(ext(d0, (100, 200, 400), tol) / n_ms - 1))) for d0 in ms_d0]
    ms_tight = float(np.max(np.abs(ext(ms_d0[1], (100, 200, 400), tol / 10) / n_ms - 1)))
    return {"free_molecule_limit_dev": fm[-1],
            "free_molecule_rate": float(np.log(fm[-2] / fm[-1]) / np.log(FM_D0[-2] / FM_D0[-1])),
            "maxwell_stefan_limit_dev": ms_tight,
            "maxwell_stefan_rate": float(np.log(ms[-1] / ms[-2]) / np.log(ms_d0[-1] / ms_d0[-2])),
            "maxwell_stefan_tol_effect": abs(ms[1] / ms_tight - 1.0),
            "_fm": fm, "_ms": ms}


V4 = v4_study()
for d0, dev in zip(FM_D0, V4["_fm"]):
    print(f"  free molecule  d_0 = {d0*1e9:7.3f} nm   departure {dev:.4e}")
for d0, dev in zip(MS_D0, V4["_ms"]):
    print(f"  Maxwell-Stefan d_0 = {d0*1e3:7.1f} mm   departure {dev:.1e}")
print(f"  rates: free molecule {V4['free_molecule_rate']:.2f} in d_0, Maxwell-Stefan "
      f"{V4['maxwell_stefan_rate']:.2f}; at 1 mm a tenfold tighter Newton tolerance moves "
      f"the departure by less than 2 %: {V4['maxwell_stefan_tol_effect'] < 0.02}")

# %% [markdown]
# ### V5: the collapse to `A4.4`'s scalar result
#
# Page [`A4.4`](../A4.4-knudsen-bosanquet/index.ipynb) derives, for a binary with
# flux ratio $\rho = -N_2/N_1$,
#
# $$
# \frac{1}{D_1} = \frac{1}{\mathrm{D}^e_{1M}} + \frac{x_2 + \rho x_1}{\mathrm{D}^e_{12}}
# \qquad\text{(eq. 109 at } \mathrm{D}_{\rm visc} = 0)
# $$
#
# and the Bosanquet formula eq. (110) as its $\rho = 1$ case. A pellet running
# $A\to\nu_B B$ has $\rho = \nu_B$ fixed by stoichiometry, so the *local* effective
# diffusivity $D_A = -N_A/\nabla c_A$ that this page's matrix route produces must
# equal that expression, face by face. $\mathrm{D}^e_{ij}$ is inversely
# proportional to pressure and the pressure is *not* uniform in a reacting pellet,
# so the comparison uses the **local** $\mathrm{D}^e_{ij}$.
#
# **Shared:** the diffusivity values. **Not shared:** the assembly; one route
# inverts a $2\times 2$ matrix, the other evaluates a scalar formula. It tests
# eq. (87) and the inversion, and cannot see an error in eq. (85).

# %%
def eq109_study(**kw):
    """Matrix-route effective diffusivity of A against eq. (109), face by face."""
    worst, rows = 0.0, []
    for nu_b, kn in [(2.0, 1.0), (2.0, 10.0), (3.0, 0.5), (1.0, 2.0)]:
        m, c = pellet(nu_b, kn, 8.0, x_as=0.6, n_cell=200, **kw)
        cf, g = m.faces(c)
        N = m.flux(c)
        ok = np.abs(g[:, 0]) > 1e-6 * np.abs(g[:, 0]).max()     # skip the centre face
        d_matrix = -N[ok, 0] / g[ok, 0]
        # x_k / D^e_ik = c_k / (D_ik c_t): local pressure through c
        d_eq109 = 1.0 / (1.0 / m.diM[0] + (cf[ok, 1] + nu_b * cf[ok, 0]) / kn)
        d_eq110 = 1.0 / (1.0 / m.diM[0] + cf[ok].sum(axis=1) / kn)
        dev = float(np.max(np.abs(d_matrix / d_eq109 - 1)))
        rows.append((nu_b, kn, dev, float(np.max(np.abs(d_eq110 / d_eq109 - 1)))))
        worst = max(worst, dev)
    return {"eq109_collapse_worst": worst, "_rows": rows}


V5 = eq109_study()
print(f"matrix route against eq. (109), worst face, within 1e-10: {V5['eq109_collapse_worst'] < 1e-10}")
print(f"{'nu_B':>5} {'Kn':>6} {'eq.(110) vs eq.(109)':>22}")
for nu_b, kn, dev, bos in V5["_rows"]:
    print(f"{nu_b:5.1f} {kn:6.2f} {bos:22.3e}")

# %% [markdown]
# ### V6: the four exponents stated on page 887
#
# The review says the fluxes are "independent of system pressure and pore size"
# under bulk control and "directly proportional to" both under Knudsen control.
# Four exponents, read off the page into the dataset as 0, 0, 1, 1, and compared
# here with log–log slopes of the $n$-component solver over a decade in pressure
# and in pore size. `A4.4` checks the binary scalar version of the same claim; this
# is the matrix version.

# %%
def slopes(d0, p_range, d0_range, pair_const=False, **kw):
    """Log-log slopes of the three fluxes against pressure and against pore size."""
    ps, ds = np.geomspace(*p_range, 6), np.geomspace(*d0_range, 6)
    Np = []
    for pp in ps:
        extra = {"dij_ct": pair_matrix(D_PAIR, pp / (RGAS * T_REV))} if pair_const else {}
        Np.append(mid_flux(d0, p=pp, **extra, **kw))
    Nd = np.array([mid_flux(dd, **kw) for dd in ds])
    Np = np.array(Np)
    sp_ = [np.polyfit(np.log(ps), np.log(np.abs(Np[:, i])), 1)[0] for i in range(3)]
    sd_ = [np.polyfit(np.log(ds), np.log(np.abs(Nd[:, i])), 1)[0] for i in range(3)]
    return np.array(sp_), np.array(sd_)


def exponent_study(**kw):
    """Worst departure of the computed slopes from the four printed exponents."""
    sp_kn, sd_kn = slopes(5e-10, (1e3, 1e4), (3e-10, 1e-9), **kw)
    sp_bk, sd_bk = slopes(1e-3, (1e3, 1e4), (3e-4, 3e-3), **kw)
    dev = max(float(np.max(np.abs(sp_kn - EXP_KN_P))), float(np.max(np.abs(sd_kn - EXP_KN_D))),
              float(np.max(np.abs(sp_bk - EXP_BK_P))), float(np.max(np.abs(sd_bk - EXP_BK_D))))
    return {"p887_exponent_worst_dev": dev, "_s": (sp_kn, sd_kn, sp_bk, sd_bk)}


V6 = exponent_study()
sp_kn, sd_kn, sp_bk, sd_bk = V6["_s"]
print(f"{'regime':<26} {'d ln N / d ln p':>34} {'d ln N / d ln d_0':>34}")
print(f"{'Knudsen (d_0 = 0.5 nm)':<26} {str(np.round(sp_kn, 4)):>34} {str(np.round(sd_kn, 4)):>34}")
print(f"{'bulk (d_0 = 1 mm)':<26} {str(np.round(sp_bk, 4)):>34} {str(np.round(sd_bk, 4)):>34}")
print(f"\nstated on p. 887: Knudsen {EXP_KN_P:.0f} and {EXP_KN_D:.0f}; bulk {EXP_BK_P:.0f} and "
      f"{EXP_BK_D:.0f}; worst departure {V6['p887_exponent_worst_dev']:.4f}")

# %% [markdown]
# ### What each check can and cannot see
#
# | check | catches | blind to |
# |---|---|---|
# | **V1** collocation and grid ladder | the matrix assembly, the sign of eq. (87), the Knudsen term, the discretisation, the closure | a mistake inside `knudsen_D` or the pair products, which both routes use, and the diffusivity *values* |
# | **V2** eq. (103) | the viscous term's prefactor, eq. (91)'s 32, $\alpha'$, the face-composition rule | anything multicomponent, and eq. (85)'s $d_0/3$ (its reference calls `knudsen_D`) |
# | **V3** Graham | eq. (85)'s mass exponent, eq. (87)'s antisymmetry | the pair diffusivity values, any error that keeps antisymmetry, eq. (85)'s $d_0/3$ (it cancels); an algebraic identity |
# | **V4** asymptotes | that both limits are reached, at the right rate | eq. (85) itself, shared with the reference |
# | **V5** eq. (109) | eq. (87) and the $2\times2$ inversion | eq. (85), which both routes use |
# | **V6** p. 887 exponents | the pressure dependence of $\mathrm{D}^e_{ij}$ and the pore-size dependence of $\mathrm{D}^e_{iM}$ | anything that scales both regimes together, including a wrong constant in front of $d_0$ |
#
# Not pinned by anything on the page: eq. (85)'s $d_0/3$ (measured above), eq. (83)'s
# $\varepsilon/\tau$ (it multiplies every diffusivity equally), the pair diffusivity
# values (every route is handed the same three numbers), and any measurement (data
# tier 6).
#
# ### Claims
#
# Every number the page states, with how it is checked. Values from solves and root
# finds are held to `FIXED`: the solves are converged to a Newton step of 1e-12 of
# the largest concentration and the roots to 1e-12, far inside it; the location of
# the viscous maximum is held to the bounded maximiser's tolerance, and the two
# values taken near the bulk limit (1 mm), where the Newton tolerance enters the
# fluxes increasingly as d0 grows, and the 100-cell discretisation gaps derived from them,
# to 1e-3 (the Maxwell-Stefan departure there to a band); the rest to a hundred times their spread over three BLAS kernels. Residuals,
# identities, two-route gaps, tolerance effects, orders and rates are bounds with
# their basis.

# %%
FIXED = 1e-6
STUDIES = {"bmatrix": BM, "transcription": TR, "inert": IN, "window_1mm": UW,
           "threshold": TH, "equal_end": EQ, "window": WN, "grid_spread": GS,
           "jackson": JK, "visc_group": VG, "viscous": VS, "collocation": V1,
           "eq103": V2, "graham": V3, "prefactor": PF, "v4": V4, "eq109": V5,
           "exponents": V6}
ALL = {k: v for s in STUDIES.values() for k, v in s.items() if not k.startswith("_")}
BOUNDS = {
    "bmatrix_wall_free_null_vector": dict(max=1e-14, basis="identity",
        note="column sums of the wall-free [B^e] vanish by D_ij = D_ji; round-off"),
    "bmatrix_colsum_knudsen_dev": dict(max=1e-14, basis="identity",
        note="column sums of [B^e] equal 1/D^e_jM; round-off"),
    "bmatrix_column_dominance_margin": dict(min=1e-3, basis="chosen",
        note="column sum over B_jj is positive by eq. (87); 1e-3 lies far below every draw"),
    "inert_scalar_closure_variation": dict(max=1e-12, basis="identity",
        note="zero flux and a finite positive D_C force a flat profile"),
    "inert_flux_ratio": dict(max=1e-10, basis="tolerance",
        note="the inert's flux vanishes by conservation; the gap is the Newton tolerance"),
    "jackson_kn_rate": dict(min=-1.01, max=-0.99, basis="theory",
        note="departure from sqrt(nu_B) falls as 1/Kn between Kn = 1e4 and 1e6"),
    "jackson_phi_effect": dict(max=1e-3, basis="chosen",
        note="conversion complete at phi = 60: phi = 200 gives the same departure to "
             "the accuracy of the 800-cell extrapolation (about 2e-4)"),
    "visc_group_closed_form_dev": dict(max=1e-14, basis="identity",
        note="eqs. (85), (91), (106) assembled against 3 p d_0/(32 eta v_bar)"),
    "collocation_extrapolated_dev": dict(max=1e-9, basis="tolerance",
        note="extrapolated finite-volume flux against solve_bvp at tol 1e-10"),
    "grid_order": dict(min=1.95, max=2.05, basis="theory",
        note="order farthest from 2 on the 25-800 cell ladder; central faces: second order"),
    "grid_ladder_over_roundoff": dict(min=1e3, basis="theory",
        note="finest relative successive difference over the Newton tolerance"),
    "eq103_viscous_dev": dict(max=1e-10, basis="identity",
        note="the arithmetic face mean makes the discrete single-component flux exact"),
    "graham_residual_worst": dict(max=1e-9, basis="tolerance",
        note="identity of eq. (82); the Newton tolerance enters increasingly as d0 grows (1 mm measured)"),
    "pressure_uniformity": dict(max=1e-10, basis="tolerance",
        note="equal end pressures and Graham's law keep p uniform; Newton tolerance"),
    "break_eq85_prefactor_v2": dict(max=1e-10, basis="identity",
        note="d_0/2 in both sides of eq. (103) leaves the identity intact"),
    "break_eq85_prefactor_v3": dict(max=1e-9, basis="tolerance",
        note="the prefactor cancels from Graham's sum; Newton tolerance"),
    "free_molecule_rate": dict(min=0.99, max=1.01, basis="theory",
        note="departure from N_i = -D^e_iM grad c_i is first order in d_0 (0.01 to 0.1 nm)"),
    "maxwell_stefan_rate": dict(min=-1.01, max=-0.99, basis="theory",
        note="departure from Maxwell-Stefan falls as 1/d_0 (1 mm to 10 mm)"),
    "collocation_worst_dev": dict(max=1e-5, basis="chosen",
        note="two-route gap on 100 cells, the second-order discretisation error; it moves "
             "by 1e-3 of itself over three BLAS kernels through the Newton tolerance at 1 mm"),
    "grid_error_n100": dict(max=1e-6, basis="theory",
        note="second-order error of the 100-cell grid at 100 nm against the collocation route"),
    "prefactor_v1_ratio": dict(max=20.0, basis="chosen",
        note="d_0/2 in both routes moves the gap by the discretisation error's change with "
             "Kn; a defect V1 sees moves it by orders of magnitude"),
    "maxwell_stefan_limit_dev": dict(min=5.0e-4, max=5.8e-4, basis="chosen",
        note="departure at 1 mm, a physical value near 5.4e-4 whose spread over three BLAS "
             "kernels (2.2e-3 of itself) is too large for a value claim: held to a band"),
    "maxwell_stefan_tol_effect": dict(max=0.02, basis="tolerance",
        note="tenfold tighter Newton tolerance at 1 mm, relative to the departure"),
    "eq109_collapse_worst": dict(max=1e-10, basis="identity",
        note="2x2 inversion against eq. (109), face by face"),
}
SLOW = {"viscous_relief_worst_d0_nm": 1e-3, "viscous_relief_worst_at_kn": 1e-3,
        # near the bulk limit the Newton tolerance enters the flux increasingly as d0 grows
        "break_eq85_prefactor_v1": 1e-3,
        # measured spread over three BLAS kernels
        "uphill_threshold_grid_spread": 1e-4, "p887_exponent_worst_dev": 1e-5}
for k, v in ALL.items():
    if k in BOUNDS:
        claims.bound(k, v, **BOUNDS[k])
    elif k == "symbolic_nonzero_residuals":
        claims.exact(k, v, note="residuals of eqs. (105), (107), (103) and the friction sum")
    elif k == "uphill_at_equal_end_compositions":
        claims.count(k, EQ["_ratios"], EQ_THR, min_margin=10.0, above=True,
                     note="pore diameters, of 25, where |N_N2/N_H2| exceeds 1e-6 at equal ends")
    else:
        claims.value(k, v, rtol=SLOW.get(k, FIXED))
print(f"{len(claims.metrics)} claims registered")

# %% [markdown]
# ### The defect-injection table
#
# Every claim above needs a row that fails its own check. Each row injects one
# defect through a call argument (eq. 87's off-diagonals, the wall term, eq. 85's
# mass exponent or prefactor, a pair diffusivity, eq. 91's 32, $\alpha'$, the
# face-composition rule, the Newton tolerance, the scalar closure, the
# stoichiometric mass ratio, the grids or the conversion of a study) and reruns the
# whole study the claim comes from, so every claim the defect reaches is tested. A
# study that cannot be completed under a defect records `nan` for its claims, which
# fails them: deleting the wall term makes $[B^e]$ singular (shown by the null
# vector above), and the solve has no answer.

# %%
def scaled_pair(f, pair=(0, 1)):
    """Pair products with D(H2-N2) multiplied by f."""
    d = dict(D_PAIR); d[pair] = D_PAIR[pair] * f
    return pair_matrix(d, CT0)


SIGN, NOWALL = dict(b_offdiag=0.9 * B_OFF), dict(knudsen=False)
MASS, LEFT = dict(mass_exp=-KN_EXP), dict(face_rule="left")
ROWS = [
    ("eq. (87) off-diagonals x 0.9: [B^e] column sums", bmatrix_study, SIGN),
    ("eq. (87) off-diagonals x 0.9: collocation and grid", collocation_study, SIGN),
    ("eq. (87) off-diagonals x 0.9: Graham's law", graham_study, SIGN),
    ("eq. (87) off-diagonals x 0.9: eq. (109) collapse", eq109_study, SIGN),
    ("eq. (87) off-diagonals x 0.9: asymptotes", v4_study, SIGN),
    ("eq. (87) off-diagonals x 0.9: uphill window at 1 mm", window_1mm_study, SIGN),
    ("eq. (87) off-diagonals x 0.9: pellet pressures", jackson_study, SIGN),
    ("wall term deleted from the full [B^e]", bmatrix_study, dict(wall_in_full=False)),
    ("wall term kept in the wall-free [B^e]", bmatrix_study, dict(wall_in_free=True)),
    ("wall term deleted: collocation and grid", collocation_study, NOWALL),
    ("eq. (85) mass exponent +1/2, finite-volume route: collocation", collocation_study, MASS),
    ("eq. (85) mass exponent +1/2: Graham's law", graham_study, MASS),
    ("eq. (85) mass exponent +1/2, finite-volume route: asymptotes", v4_study, MASS),
    ("eq. (85) mass exponent +1/2: pellet pressures", jackson_study, MASS),
    ("eq. (85) mass exponent +1/2 in every route: prefactor study", prefactor_study, MASS),
    ("eq. (85) d0/2 in the finite-volume route only: prefactor study", prefactor_study,
     dict(den_ref=KN_DEN)),
    ("eq. (85) d0/2: uphill window at 1 mm", window_1mm_study, dict(kn_den=2.0)),
    ("eq. (85) d0/2: uphill window over pore sizes", window_study, dict(kn_den=2.0)),
    ("eq. (85) d0/2: Kn and the viscous group", visc_group_study, dict(kn_den=2.0)),
    ("eq. (85) d0/2: viscous sweep", viscous_study, dict(kn_den=2.0)),
    ("D(H2-N2) x 1.1: thresholds", threshold_study, dict(dij_ct=scaled_pair(1.1))),
    ("eq. (85) mass dependence dropped: uphill window over pore sizes", window_study,
     dict(mass_exp=0.0)),
    ("scalar closure (diagonal of [B^e]): inert pellets", inert_study, dict(scalar=True)),
    ("scalar closure: collocation and grid", collocation_study, dict(scalar=True)),
    ("scalar closure: equal end compositions", equal_end_study, dict(scalar=True)),
    ("inert contrast run with the full matrix", inert_study, dict(contrast_scalar=False)),
    ("inert consumed at 1 % of the rate of A", inert_study, dict(sink_extra=0.01)),
    ("D_ji = 2 D_ij in the symbolic sum", transcription_study, dict(asym=2)),
    ("interior face composition from the left cell: collocation and grid",
     collocation_study, LEFT),
    ("threshold grid spread over 10 to 80 cells instead of 50 to 400",
     grid_spread_study, dict(ns=(10, 20, 40, 80))),
    ("interior face composition from the left cell: eq. (103)", eq103_study, LEFT),
    ("Newton tolerance 1e-12 -> 1e-6: collocation and grid", collocation_study, dict(tol=1e-6)),
    ("Maxwell-Stefan limit at 30 mm with Newton tolerance 1e-11: asymptotes", v4_study,
     dict(tol=1e-11, ms_d0=(1e-3, 3e-2, 1e-1))),
    ("viscous term deleted (B0 = 0): eq. (103)", eq103_study, dict(B0=0.0)),
    ("eq. (91) 32 -> 8: eq. (103)", eq103_study, dict(B0=D0_103 ** 2 / 8)),
    ("alpha' = 0 instead of 1: eq. (103)", eq103_study, dict(alpha=0.0)),
    ("eq. (91) 32 -> 8 in the assembled group", visc_group_study, dict(pois_den=8.0)),
    ("eq. (91) 32 -> 8: viscous sweep", viscous_study, dict(pois_den=8.0)),
    ("M_B/M_A = 1 (mass not conserved): pellet pressures", jackson_study,
     dict(mass_ratio=1.0)),
    ("conversion incomplete (phi = 6 and 20): pellet pressures", jackson_study,
     dict(phis=(6.0, 20.0))),
    ("pair diffusivities held at their 1 atm value: exponents", exponent_study,
     dict(pair_const=True)),
]
BASE = {bmatrix_study: BM, transcription_study: TR, inert_study: IN,
        window_1mm_study: UW, threshold_study: TH, equal_end_study: EQ,
        window_study: WN, grid_spread_study: GS, jackson_study: JK,
        visc_group_study: VG, viscous_study: VS, collocation_study: V1,
        eq103_study: V2, graham_study: V3, prefactor_study: PF, v4_study: V4,
        eq109_study: V5, exponent_study: V6}       # each study's baseline result
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

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises if any check fails.
claims.report();   # the semicolon keeps the notebook from echoing the payload

# %% [markdown]
# ## What pymrm adds
#
# The review states the dusty gas model completely and then stops: it prints
# eqs. (82) to (87) and (98) to (110), shows one figure of computed fluxes, and makes
# four claims in words that it does not quantify. pymrm supplies the transport
# skeleton of the $n$-component solver (`construct_grad`, `construct_div`,
# `construct_boundary_value_matrices` for the face compositions, `NumJac` with a
# cell-and-species stencil, `newton`) and `pymrm.checks` for the root scans and the
# refinement ladders. What the page adds on top:
#
# - **An $n$-component dusty gas solver** covering a capillary or a pellet, any
#   number of species, with or without reaction and viscous flow, checked against
#   an independent collocation solution of the same equations that shares no
#   transport code, and against the review's own eqs. (103), (105), (107) and
#   `A4.4`'s eq. (109).
# - **"Not very serious" (Haynes 1978), quantified.** The screening group
#   $\mathrm{D}_{\rm visc}/\mathrm{D}^e_{iM} = 3pd_0/(32\eta\bar v_i)$ cannot be
#   swept alone, because eq. (85) makes the Knudsen number a function of the same
#   $d_0$. Tying both to $d_0$, viscous flow changes the pellet-centre pressure by
#   under a percent at every pore size from 1 nm to 100 µm; freezing Kn instead
#   overstates the 1 µm figure about elevenfold.
# - **Kaza & Jackson's uphill diffusion, given the right axis**: the threshold
#   pore diameter depends on composition, and what pore size controls is the width
#   of the composition window, nearly proportional to $d_0$.
# - **Jackson's $p_0 = \sqrt{\nu_B}\,p$, given its preconditions**, as the corner
#   case of a closed form the solver reproduces, with the departure at finite Kn
#   measured as $1/(2\,\mathrm{Kn})$.
# - **A demonstration that a scalar diffusivity is not merely inaccurate**: the inert
#   species in a reacting pellet.
#
# What this page does **not** do: it compares nothing with a measurement (data tier
# 6). It does not treat surface diffusion, pore-size distributions, viscous
# selectivity ($\alpha'_i \neq 1$), thermal diffusion or non-ideal mixtures, though
# the review discusses all of them. It does not re-derive `A4.4`'s Bosanquet error,
# only checks against it.

# %% [markdown]
# ## Reuse
#
# **To adapt this page to your system**, change in order:
#
# | Change | Where |
# |---|---|
# | the species and their molar masses | `SPEC`, `MOLAR` |
# | the binary pair diffusivities | `D_PAIR`, then `pair_matrix` |
# | pore size, porosity, tortuosity | the `d0` argument and `EPS_TAU` |
# | slab / cylinder / sphere | `geom` in `DustyGas` (`construct_div`'s `nu`) |
# | boundary conditions | `SYMMETRY`, `dirichlet`, both on the outward normal |
# | the reaction | `sink` and `kvol` in `DustyGas` |
# | viscous flow on or off | `B0` (eq. 91 gives $d_0^2/32$), `eta`, `alpha` |
#
# **The reusable pattern.** `DustyGas.flux` is a *constitutive closure that solves
# a small dense linear system at every face* (`np.linalg.solve` on a stacked
# `(n_faces, n, n)` array) while the transport skeleton stays pure pymrm. Any `S9`
# model whose flux law is implicit in the fluxes has this shape: Maxwell–Stefan
# (`A4.2`), Nernst–Planck with electroneutrality, or a Maxwell–Stefan surface
# diffusion model. Storing $\mathrm{D}_{ij}c_t$ rather than $\mathrm{D}_{ij}$ lets
# one object serve every pressure. Near the bulk limit the fluxes respond to the
# total-pressure gradient, so the Newton tolerance
# on the concentrations enters the fluxes increasingly as d0 grows: measure its
# effect there.
#
# **A screening rule that falls out, and a trap inside it.** Before reaching for
# the full matrix, form two ratios: $\mathrm{Kn} = \mathrm{D}^e_{ij}/\mathrm{D}^e_{iM}$
# and $3pd_0/(32\eta\bar v_i)$. If the first is small the species are strongly
# coupled and a scalar diffusivity will mislead; the second says how much viscous
# conductance the pore has. **They are not independent**: the first is inverse-linear
# in $d_0$ and the second linear, so a sweep that holds one fixed while varying the
# other describes no pore that exists. If two dimensionless groups share a physical
# variable, sweep the variable.
#
# **Related pages:**
# [`A4.4`](../A4.4-knudsen-bosanquet/index.ipynb) the scalar closure this is the
# matrix form of · [`A4.2`](../A4.2-maxwell-stefan-vs-fick/index.ipynb) the same
# matrix without a wall · [`A4.9`](../A4.9-duncan-toor/index.ipynb) uphill
# diffusion measured in free space ·
# [`B1.1`](../B1.1-thiele-weisz-hicks/index.ipynb) the pellet this replaces the
# effective diffusivity of.
#
# ## References
#
# - **Krishna, R. & Wesselingh, J. A. (1997)**, "The Maxwell–Stefan approach to
#   mass transfer", *Chemical Engineering Science* **52**(6) 861–911,
#   [doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7).
#   Every equation, every quoted sentence and every constant on this page comes
#   from it.
# - Mason, E. A. & Malinauskas, A. P. (1983), *Gas Transport in Porous Media: the
#   Dusty Gas Model*; Jackson, R. (1977), *Transport in Porous Catalysts*; Mason,
#   E. A. & del Castillo, L. F. (1985): the origins of the dusty gas model and of
#   eqs. (85), (99) to (101), cited by the review with attribution and cited here
#   through it.
# - Kaza, K. R. & Jackson, R. (1980), on uphill diffusion in a catalyst particle;
#   Haynes (1978), on the size of the viscous contribution; Remick, R. R. &
#   Geankoplis, C. J. (1974), the source of the review's Fig. 44 measurements.
#   Cited through the review's sentences about them; no result here is attributed
#   to any of them.
#
# *Cite the source, not this page.*
