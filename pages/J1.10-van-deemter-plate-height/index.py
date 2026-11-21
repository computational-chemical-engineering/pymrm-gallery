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
# title: "The van Deemter equation: the three terms, the constants they were read off, and what the Gaussian costs"
# description: "H = A + B/u + Cu is fitted to twelve points on one page of a 1956 paper. Rebuilding the chain that produced it recovers the film thickness, the particle diffusivity, and one printed equation that cannot be right."
# categories: [sec:J, struct:S4, tier:T0, data:tier6, phase:gas-liquid, phase:liquid-solid]
# date: 2026-08-09
# ---
#
# # The van Deemter equation: the three terms, the constants they were read off, and what the Gaussian costs
#
# **Catalog ID:** `J1.10` · **Structures:** `S4` (1-D transient PDE) · **Tier:** T0
#
# Every chromatography textbook prints
#
# $$H = A + \frac{B}{u} + C u,$$
#
# attributes it to van Deemter, Zuiderweg and Klinkenberg (1956), and moves on.
# The 1956 paper does not start there. It starts from a two-phase material
# balance, reduces it to a Gaussian, and only then, on page 286, from twelve
# measured points, writes down three numbers. This page rebuilds that chain,
# end to end, and reports what each step actually supports.

# %% [markdown]
# ## Background
#
# A band injected into a chromatographic column comes out wider than it went in.
# The plate theory of Martin and Synge describes *how much* wider with one
# empirical number per column, the height equivalent to a theoretical plate, and
# says nothing about where that number comes from. The rate theory of Lapidus and
# Amundson describes the mechanisms, longitudinal diffusion, finite mass
# transfer, but its solution is an integral over a Bessel function, which in
# 1956 was not something you evaluated for a design.
#
# The 1956 paper is the join. Its abstract lists the three causes it treats:
#
# > axial molecular diffusion, axial eddy diffusion, finiteness of transfer
# > coefficient
#
# and its contribution is that under conditions it states explicitly, the rate
# theory collapses to a Gaussian whose width is *additive* in those three
# mechanisms. That additivity is the van Deemter equation. It is why the
# familiar $A + B/u + Cu$ has one term per mechanism and no cross terms.
#
# Three things about the paper are easy to miss and matter here.
#
# 1. **Two experiments, not one.** Sections apply the theory to *ion exclusion*
#    (Simpson and Wheaton's data, restated in Table 2, reduced to plate heights in
#    Table 3) and to *gas–liquid partition chromatography* (the authors' own
#    Table 4). Only the second produces the famous three-term formula; the first
#    produces the particle diffusivity and the packing factor $\lambda$ that the
#    second then reuses.
# 2. **The constants were read off graphs.** Eq. (54) is a description of
#    Fig. 11, and the paper says so, "the values of the constants in (54) are
#    highly uncertain, owing to the relatively small number of experimental data
#    that were available".
# 3. **The paper criticises its own fit.** The eddy term "should obviously have
#    the same value for both components. The formulae (54) do not agree with
#    theory in this respect", 0.3 mm for *n*-butane against 0.8 mm for
#    iso-butane.

# %% [markdown]
# ## The published model
#
# **The balance (eqs. 29 and 30, p. 278),** per unit cross-sectional area, with
# $F_I$ and $F_{II}$ the fractional volumes of the mobile and immobile phases,
# $K = (c_I/c_{II})_{\text{equilibrium}}$ the distribution factor, $\alpha$ the
# mass transfer coefficient per unit volume of packing, and $u$ the **interstitial**
# velocity (the paper's own List of Symbols, p. 288):
#
# $$F_I\frac{\partial c_I}{\partial t}
#   = F_I D\frac{\partial^2 c_I}{\partial z^2}
#   - F_I u\frac{\partial c_I}{\partial z}
#   + \alpha\,(K c_{II} - c_I),
# \qquad
# F_{II}\frac{\partial c_{II}}{\partial t} = \alpha\,(c_I - K c_{II}).$$
#
# **The Gaussian (eqs. 33 and 34, p. 279).** For $z$ much larger than both
# $2D/u$ and $F_I u/\alpha$, the Lapidus–Amundson solution reduces to
#
# $$\frac{c_I}{c_0} = \frac{\beta t_0}{\sqrt{2\pi(\sigma_1^2+\sigma_2^2)}}
#    \exp\!\left[-\frac{(z/u-\beta t)^2}{2(\sigma_1^2+\sigma_2^2)}\right],
# \qquad
# \frac{1}{\beta} = 1 + \frac{F_{II}}{F_I K},\quad
# \sigma_1^2 = \frac{2Dz}{u^3},\quad
# \sigma_2^2 = 2\beta^2\frac{F_{II}^2 z}{\alpha F_I K^2 u}.$$
#
# **The plate height (eq. 38, p. 279).** Comparing that Gaussian with the plate
# theory's gives
#
# $$H = 2\frac{D}{u} + 2\,\frac{u F_I/\alpha}{\left(1 + K\dfrac{F_I}{F_{II}}\right)^2}.$$
#
# **The three terms (eqs. 47, 48 and 53, pp. 282–283).** Writing the longitudinal
# diffusivity as $D = \gamma D_I + \lambda u d_p$ and the transfer resistance as
# diffusion into a liquid film of thickness $d_f$ (eqs. 50–52):
#
# $$H = \underbrace{2\gamma\frac{D_I}{u}}_{B/u}
#     + \underbrace{2\lambda d_p}_{A}
#     + \underbrace{\frac{8}{\pi^2}\,
#       \frac{K F_I d_f^2}{\left(1+K\frac{F_I}{F_{II}}\right)^2 D_{II} F_{II}}\;u}_{Cu}.$$
#
# That is the van Deemter equation. Each of the three constants is a *product of
# printed quantities*, which is what makes the rest of this page possible: the
# paper inverts the relation to obtain $\gamma D_I$, $\lambda$ and $d_f$, and each
# inversion can be redone.
#
# **The particle form (eq. 49, p. 282)**, used for the ion-exclusion column, where
# the resistance is inside a solid particle rather than a film:
#
# $$\frac{H}{d_p} = 2\gamma\frac{D_I}{u d_p} + 2\lambda + C\,\frac{u d_p}{D_I},
# \qquad
# C = \frac{\dfrac{F_I^{2}}{75(1-F_I)^{2}}
#         + \dfrac{F_I K D_I}{2\pi^{2}(1-F_I)D_{II}}}
#         {\left(1+K\dfrac{F_I}{F_{II}}\right)^{2}}.$$

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.21"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "J1.10-van-deemter-plate-height"  # equals the page directory

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
# Imports. The page's model.py is put on sys.path by setup().
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.special import erfinv
import pymrm
from pymrm import vanleer
from pymrm.checks import find_roots, observed_orders
from gallery_utils import Claims, load_data, load_meta, cite_data
from model import (Column, cumulants, plate_height, kappas_per_z, one_root,
                   width_red, cmax, n_from_height, n_from_width, vd, refit, fit_eq49,
                   dii_from_slope)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])
SQ2 = np.sqrt(2.0)
EPS = np.finfo(float).eps


# Deviations are reported ONE way everywhere on this page:
#     dev = (computed - printed) / printed
def dev(computed, printed):
    return np.asarray(computed, float) / np.asarray(printed, float) - 1.0


# %% [markdown]
# ## Parameters and assumptions
#
# **What is printed, and what is not.** The gas–liquid columns are described
# completely enough to evaluate eq. (53): both are 1.80 m long and 6 mm in
# internal diameter, the Celite carrier averages 30 µm and the Sterchamol 240 µm,
# $K F_I/F_{II}$ is given for all three column/solute combinations, and $D_{II}$
# in the oil film is estimated from its viscosity as $0.3\times10^{-10}$ m² s⁻¹.
# What is **not** printed for those columns is $F_I$ itself. That matters, and the
# page returns to it.
#
# The ion-exclusion column is described the other way round: $F_I = F_{II} = 0.4$
# and $K = 1.3$ are printed (p. 282), but its **length is not**, anywhere. The
# length this page uses is recovered from Table 3's own columns - the value is
# computed and printed in *Table 3 checked against itself*, below - and is
# labelled an inference every time it appears.
#
# **Assumptions carried from the paper, unchanged.** Linear isotherm; constant
# $F_I$, $F_{II}$ and $u$; a single solute at low concentration; eddy diffusivity
# $E = \lambda u d_p$ (eq. 47); gas-phase transfer resistance neglected in the
# gas–liquid case, which the paper checks afterwards with its eq. (56).
#
# **One convention this page does not adopt silently.** The List of Symbols
# defines $u$ as the interstitial velocity and eqs. (29)/(30) use it that way, but
# Tables 2, 3 and 4 and Figs. 11–14 are all headed *superficial* velocity. Fig. 9
# converts (its abscissa is $u F_I d_p$); the gas–liquid section does not appear
# to. The page carries both readings explicitly wherever the difference changes a
# number.

# %% [markdown]
# ## The data
#
# Five files, all transcription, no digitisation. The paper has **four** printed
# tables, not two, and the ones that carry the plate height are Tables 3 and 4.
#
# | file | what it is | who measured it |
# |---|---|---|
# | `…table1-feed-volume` | Table 1, p. 277: elution-curve width and height against feed volume, 24 rows | the authors |
# | `…table2-simpson-wheaton` | Table 2, p. 280: ion-exclusion elution curves, 13 rows | **Simpson & Wheaton**, restated |
# | `…table3-hetp` | Table 3, p. 281: plate heights derived from Table 2, 10 rows (one printed twice) | derived by the authors |
# | `…table4-glc-hetp` | Table 4, p. 286: plate height against gas velocity, two columns, 24 rows | the authors |
# | `…printed-scalars` | 56 scalars from pp. 277–287, tagged `input` or `result` | (none) |
#
# **How they were read.** Every numeric cell was read on a 2× digit-scale crop of
# a 300 ppi render, the scan's native resolution, checked with `pdfimages -list`.
# The text layer was used for nothing numeric: it misreads the journal *volume* in
# the document's first line, printing `Vol 6` where the page prints `Vol. 5`. The
# 1956 Pergamon face sets the decimal point as a raised mid-dot which the scan
# drops to a space, so a cell reading `0 163 1100 0 27 0 30` is
# `0` `163` `1100` `0.27` `0.30`.
#
# **Five readings were settled by arithmetic constraints the tables impose on
# themselves, not by pixel shape**, and the sidecars name the constraint for each:
#
# * Table 2's velocity column is a **doubling series**, 0.046, 0.091, 0.182,
#   0.365, 0.73, each twice the last to within the printed rounding. That fixes
#   `0 365` as 0.365 and not 0.865.
# * The 23 ml base width on that same row is the only reading consistent with the
#   125–156 plates Table 3 prints for it.
# * Table 3's `11 -18` could be 11–18 or 1.1–1.8. $H$ divided by $H/d_p$ is the
#   particle diameter, and it must agree with the mesh: 0.56/11 = 0.051 mm is a
#   200–400 mesh particle, 0.56/1.1 = 0.51 mm is a 30–40 mesh one. The row is
#   200–400 mesh.
# * The same test fixes `13 -14` and `16 -20` as integer pairs.
# * Table 3's `0 56- 0 93` keeps its decimal points where $nH$, the column
#   length, which must be the same in every row, puts them.
#
# **Figs. 11 and 13 were consulted for three readings, spanning five cells of
# Table 4**: the Celite `3 4`/`3 5` pair, the Sterchamol `3 2`/`3 3` pair, and
# `31` against `81` (the 1956 face sets `3` with a nearly flat top). No point was
# taken off either figure; what was used is that neither carries a point between
# $u = 5$ and 15 mm s⁻¹, and that Fig. 13's rightmost point sits at $u \approx 31$
# on an axis stopping just past it. **All three are also forced from inside
# Table 4 alone**, the stronger statement, and the one checked in code below. The
# check is run on the cells *not* in dispute, because testing whether the
# transcribed column ascends would assume the reading it is offered as evidence
# for: the undisputed cells of each block ascend, and of the four candidate
# readings of the disputed pair only the printed one can be inserted at its two
# positions without breaking that order. Position does **not** decide the third
# cell, 81 is larger than everything above it either way, and the page does not
# claim it does; eq. (55), which the authors fitted to these very rows, overshoots
# the printed plate height by more than a factor two at $u = 81$ while landing
# within a few per cent of it at $u = 31$.
#
# Two of the 56 scalars, $n = 240$ and $n = 277$, are numbers **printed as curve
# labels** on Fig. 8 (p. 280). They are read as printed text, not as point
# coordinates, and are used in exactly one place: to identify which of Table 3's
# two printed endpoints came from which of the paper's two routes.
#
# **Simpson & Wheaton's paper is not on disk and was not consulted.** Table 2 is
# van Deemter et al.'s restatement of it. Nothing on this page claims anything
# about the original.

# %%
T1 = load_data("vandeemter-1956-table1-feed-volume.csv", page=PAGE)
T2 = load_data("vandeemter-1956-table2-simpson-wheaton.csv", page=PAGE)
T3 = load_data("vandeemter-1956-table3-hetp.csv", page=PAGE)
T4 = load_data("vandeemter-1956-table4-glc-hetp.csv", page=PAGE)
SC = load_data("vandeemter-1956-printed-scalars.csv", page=PAGE)
P = dict(zip(SC.quantity, SC.value))                 # every printed scalar, by key

print(cite_data(load_meta("vandeemter-1956-table3-hetp.csv", page=PAGE)))
print(f"\nTable 1 {len(T1):3d} rows   Table 2 {len(T2):3d} rows   "
      f"Table 3 {len(T3):3d} rows   Table 4 {len(T4):3d} rows")
print(f"printed scalars: {len(SC)} ({(SC.kind=='input').sum()} inputs, "
      f"{(SC.kind=='result').sum()} stated results), journal pp. "
      f"{int(SC.page.min())}-{int(SC.page.max())}")

# The Celite sieve analysis printed inside Table 4's carrier column.
sieve = {">40 um": 24.0, "40-20 um": 52.0, "20-10 um": 18.5, "10-6 um": 4.5, "<6 um": 1.0}
print(f"\nCelite sieve analysis as printed sums to {sum(sieve.values()):.1f} %")

# The two constraints INTERNAL to Table 4 that settle its three ambiguous
# readings without Figs. 11/13 (see the text above).
#
# The monotonicity test must be run on the cells NOT in dispute, otherwise it is
# circular: asking whether the transcribed column ascends, when the transcription
# already contains 3.4/3.5, can only ever answer True. So: establish the order on
# the undisputed cells, then insert each candidate reading at its printed
# position and see which candidates survive.
DISPUTED = {"celite": (3.4, 3.5), "sterchamol": (3.2, 3.3)}
ALTS = {"celite":     [(8.4, 8.5), (3.4, 8.5), (8.4, 3.5)],   # the leading 3 read
        "sterchamol": [(8.2, 8.3), (3.2, 8.3), (8.2, 3.3)]}   # as an 8, in either cell
ascends = lambda v: bool(np.all(np.diff(np.asarray(v, float)) >= 0))

print("\nTable 4's velocity column, tested on the cells that are NOT in dispute:")
BLOCKS = {}
for c, g in T4.groupby("solid_carrier"):
    v = g.superficial_velocity_mm_s.to_numpy(float)
    pos = [int(np.flatnonzero(np.isclose(v, x))[0]) for x in DISPUTED[c]]
    keep = np.delete(v, pos)
    BLOCKS[c] = (v, pos, keep)
    print(f"   {c:<11s}{len(keep):2d} undisputed cells "
          f"[{', '.join(f'{x:g}' for x in keep)}]  ascending: {ascends(keep)}")

POS_SETTLES = True
for c, (v, pos, keep) in BLOCKS.items():
    print(f"   {c}: inserting each candidate reading at its two printed positions")
    for cand in [DISPUTED[c]] + ALTS[c]:
        w = v.copy(); w[pos[0]], w[pos[1]] = cand
        ok = ascends(w)
        POS_SETTLES &= (ok == (cand == DISPUTED[c]))
        print(f"      {cand[0]:>4g} / {cand[1]:<4g}  whole column ascending: {ok}"
              f"{'   <- the reading used' if cand == DISPUTED[c] else ''}")
print(f"   -> position alone settles both pairs, and only the printed reading "
      f"survives: {POS_SETTLES}")

# The third ambiguous cell is NOT settled by position, and the page does not
# pretend otherwise: 81 is larger than every cell above it either way. Same
# undisputed column as above, with only its last cell varied.
v_31 = BLOCKS["sterchamol"][2]
v_81 = v_31.copy(); v_81[-1] = 81.0
print(f"\n   The last Sterchamol cell: reading it 31 leaves the undisputed column "
      f"ascending ({ascends(v_31)}),")
print(f"      and so does reading it 81 ({ascends(v_81)}) - position does NOT decide "
      f"this one.")
h55 = lambda u: P["eq55_B_sterchamol"] / u + P["eq55_C_sterchamol"] * u
h31 = float(T4[(T4.solid_carrier == "sterchamol")
               & (T4.superficial_velocity_mm_s == 31)].hetp_n_butane_mm.iloc[0])
print(f"   eq (55), fitted by the authors to these rows, gives {h55(31):.2f} mm at "
      f"u = 31 against the printed {h31} mm ({dev(h55(31), h31):+.0%}),")
print(f"      and {h55(81):.2f} mm at u = 81 ({dev(h55(81), h31):+.0%}) - which excludes "
      f"the 81 reading from inside the table.")

# %% [markdown]
# ### Table 3 checked against itself, before anything is fitted to it
#
# Table 3 prints $n$, $H$ and $H/d_p$ for the same experiment. Two of those three
# are enough to determine the third, so the table constrains itself twice over:
#
# * $nH$ is the column length, and must be the **same in every row**, even though
#   no column length is printed anywhere in the paper;
# * $H \div (H/d_p)$ is the particle diameter, and must be the **same for a given
#   mesh**, whatever the flow rate.
#
# Neither identity enters anything below. They are here because they are
# the only check available on a transcription of a table whose every cell is a
# range.

# %%
# Both printed ranges in a row run in OPPOSITE senses: the larger plate number
# pairs with the smaller plate height. That is what makes the product sharp.
Lhi = T3.n_plates_lo * T3.hetp_hi_mm
Llo = T3.n_plates_hi * T3.hetp_lo_mm
L_all = np.concatenate([Lhi.to_numpy(), Llo.to_numpy()])
L_SW = float(np.mean(L_all))

print("Recovered column length  l = n H   (mm)")
print(f"{'mesh':>8}{'u':>7}{'n_lo*H_hi':>11}{'n_hi*H_lo':>11}")
for (_, r), a, b in zip(T3.iterrows(), Lhi, Llo):
    print(f"{r.particle_mesh:>8}{r.superficial_velocity_mm_s:>7.3f}{a:>11.1f}{b:>11.1f}")
print(f"\n  mean {L_SW:.1f} mm, spread {L_all.min():.0f}-{L_all.max():.0f} mm "
      f"({(L_all.max()/L_all.min()-1)*100:.1f} % wide, against 2-sig-fig rounding "
      f"on H of up to {(0.05/T3.hetp_lo_mm).max()*100:.1f} %)")

dp_row = 0.5*(T3.hetp_lo_mm + T3.hetp_hi_mm) / (0.5*(T3.h_over_dp_lo + T3.h_over_dp_hi))
print("\nRecovered particle diameter  d_p = H / (H/d_p)   (mm)")
for mesh, g in dp_row.groupby(T3.particle_mesh):
    lo, hi = [float(x) for x in mesh.split("-")]
    print(f"  {mesh:>8} mesh : {'  '.join(f'{v:.4f}' for v in g)}"
          f"   (mean {g.mean():.4f}, {len(g)} row(s))")
uniq = ~T3.duplicated(subset=["particle_mesh", "superficial_velocity_mm_s"])
DP_50_100 = float(dp_row[(T3.particle_mesh == "50-100") & uniq].mean())
print(f"\n  50-100 mesh, the only series with a flow sweep: d_p = {DP_50_100:.4f} mm "
      f"(five distinct flow rates, row spread "
      f"{dp_row[(T3.particle_mesh=='50-100') & uniq].min():.4f}-"
      f"{dp_row[(T3.particle_mesh=='50-100') & uniq].max():.4f})")


# %% [markdown]
# ## PyMRM implementation
#
# Two things are built here. The first is a **symbolic** reconstruction of the
# paper's chain of equations, which needs no solver and is what most of this page
# rests on. The second is the **transient two-phase column** of eqs. (29)/(30),
# solved with pymrm, which is the only thing on the page that can test eq. (38)
# without using eq. (38).
#
# ### The column
#
# `Column` (in `model.py`) assembles the operators once, in `__init__`, and never
# inside a time step. Only $c_I$ carries a spatial operator; $c_{II}$ exchanges with it
# pointwise, so the two fields are assembled monolithically with `bmat` rather
# than through a shared axis, the immobile phase has no flux and no gradient, and
# giving it one would be a modelling error dressed as convenience.
#
# The inlet is a Danckwerts condition, written on the **outward** normal
# ($n = -z$ at the left face, so $\partial c/\partial n = -\partial c/\partial z$):
#
# $$u c_I - D\frac{\partial c_I}{\partial z}\bigg|_{0} = u\,c_{\text{in}}(t)
# \;\Longrightarrow\;
# D\frac{\partial c_I}{\partial n} + u\,c_I = u\,c_{\text{in}}(t)
# \;\Longrightarrow\; a = D,\; b = u,\; d = u\,c_{\text{in}}.$$
#
# Only $d$ varies in time, and the boundary vectors are linear in $d$, so they are
# built once for $c_{\text{in}} = 1$ and scaled: the matrices stay constant, which
# is the whole point of the convention.
#
# Convection is upwind plus a van Leer correction. The correction is the only
# nonlinear part of a step, and `Column.run` treats it as deferred: computed from
# the state at the start of the step (with that step's inlet value) and moved to
# the right-hand side, so each step of the Crank-Nicolson march is one solve with a
# matrix factorised once. `Column.run_newton` solves the same step with the limited
# flux implicit, through `ConvFluxTVD` and its exact Jacobian; the validation
# section measures what the deferral costs.

# %%
# Column, cumulants and plate_height live in model.py. The two boundary conditions
# as pymrm reads them, on a ten-cell copy of the column:
_demo = Column(u=0.23, D=0.17, alpha=1.0, K=1.3, F_I=0.4, F_II=0.4, L=700.0, n_z=10)
print(pymrm.describe_bc(_demo.bc, _demo.x_f))


# %% [markdown]
# ### The paper's chain, verified symbolically
#
# Four steps. For each of them the paper names the route, which equations are
# introduced into which, but never prints the manipulation, and none of the four
# is checked anywhere in the source. All four are exact.

# %%
s_, D_, u_, z_, a_, FI_, FII_, K_ = sp.symbols(
    "s D u z alpha F_I F_II K", positive=True)
dp_, DI_, DII_, g_, l_, df_ = sp.symbols(
    "d_p D_I D_II gamma lambda d_f", positive=True)

# --- what the paper prints -------------------------------------------------
beta_34 = 1 / (1 + FII_ / (FI_ * K_))                                     # (34)
sig1sq = 2 * D_ * z_ / u_ ** 3                                            # (34)
sig2sq = 2 * beta_34 ** 2 * FII_ ** 2 * z_ / (a_ * FI_ * K_ ** 2 * u_)    # (34)
H_38 = 2 * D_ / u_ + 2 * (u_ * FI_ / a_) / (1 + K_ * FI_ / FII_) ** 2     # (38)
C_49 = ((FI_ ** 2 / (75 * (1 - FI_) ** 2)
         + FI_ * K_ * DI_ / (2 * sp.pi ** 2 * (1 - FI_) * DII_))
        / (1 + K_ * FI_ / FII_) ** 2)                                     # (49)
H_49 = 2 * g_ * DI_ / (u_ * dp_) + 2 * l_ + C_49 * u_ * dp_ / DI_         # (49)
H_52 = (2 * D_ / u_ + sp.Rational(8, 1) / sp.pi ** 2 * K_ * FI_ * u_ * df_ ** 2
        / ((1 + K_ * FI_ / FII_) ** 2 * DII_ * FII_))                     # (52)

# --- 1. (38) is what (34) implies, through H = z sigma_t^2 / t_R^2 ---------
t_R = z_ / (u_ * beta_34)
H_from_34 = sp.simplify(z_ * ((sig1sq + sig2sq) / beta_34 ** 2) / t_R ** 2)
r1 = sp.simplify(H_from_34 - H_38)

# --- 2. (49) and its C are (38) with (45), (46) and (48) substituted -------
a_p = 6 * (1 - FI_) / dp_                                                 # (46)
inv_alpha = (sp.Rational(6, 25) * FI_ / (DI_ * a_p ** 2)
             + sp.Rational(3, 2) / sp.pi ** 2 * K_ * dp_ / (DII_ * a_p))  # (45)
D_48 = g_ * DI_ + l_ * u_ * dp_                                           # (47),(48)
r2 = sp.simplify(sp.simplify(H_38.subs({D_: D_48, a_: 1 / inv_alpha}) / dp_) - H_49)

# --- 3. (52) is (38) with the film resistance (41),(50),(51) --------------
alpha_II = sp.Rational(1, 4) * sp.pi ** 2 * DII_ * FII_ / df_ ** 2        # (51)
r3 = sp.simplify(H_38.subs(a_, alpha_II / K_) - H_52)   # 1/alpha = K/alpha_II

# --- 4. a route the paper never takes: Laplace transform of (29)/(30) ------
# c_I(z,s) ~ exp(-lambda(s) z) with  D lam^2 + u lam - Phi(s) = 0, so the
# cumulants of the residence time are exactly linear in z.  H = kappa2/kappa1^2.
Phi = s_ * (1 + (a_ / FI_) / (s_ + a_ * K_ / FII_))
lam_s = (-u_ + sp.sqrt(u_ ** 2 + 4 * D_ * Phi)) / (2 * D_)
ser = sp.expand(sp.series(lam_s, s_, 0, 5).removeO())
kap = [sp.simplify(ser.coeff(s_, m) * sp.factorial(m) * (-1) ** (m + 1))
       for m in range(1, 5)]                       # kappa_1..kappa_4 per unit length
r4 = sp.simplify(kap[1] / kap[0] ** 2 - H_38)

CHAIN = [("(34) -> (38)   H = z sigma_t^2 / t_R^2", r1),
         ("(38)+(45,46,48) -> (49) and its C", r2),
         ("(38)+(41,50,51) -> (52)", r3),
         ("Laplace cumulants of (29),(30) -> (38)", r4)]
for name, r in CHAIN:
    print(f"  {name:<42s} residual = {r}")
CHAIN_OK = int(sum(r == 0 for _, r in CHAIN))
print(f"\n  {CHAIN_OK} of {len(CHAIN)} exact")
print(f"\n  kappa_1 / z = {sp.simplify(kap[0])}")
print(f"  kappa_3 / z = {sp.factor(kap[2])}")
print("  kappa_4 / z comes out of the same expansion (taken to s^4); it is used in")
print("  *What pymrm adds*, so the excess kurtosis there rests on no PDE cumulant.")
kap_f = sp.lambdify((D_, u_, a_, FI_, FII_, K_), kap, "numpy")

# %% [markdown]
# ## Results
#
# ### 1. Eq. (54) against Table 4, a goodness of fit, and a null baseline beside it
#
# **The coefficients of eqs. (54) and (55) were fitted to these very rows.**
# Nothing below is a test of the van Deemter equation; it is a measure of how well
# three numbers describe twelve points, and it is reported that way. The null
# baseline is the best constant $H$, the model you get by ignoring velocity
# entirely, which is the thing the paper's Fig. 11 exists to refute.
#
# The minimum of each curve is **root-found** on $\mathrm{d}H/\mathrm{d}u = 0$ and
# then checked against the closed form $u_{\text{opt}} = \sqrt{B/C}$,
# $H_{\min} = A + 2\sqrt{BC}$. That agreement is an identity, not evidence, and is
# labelled as one; its only job is to guarantee no reported optimum is a sampled
# maximum of a sweep.

# %%
FITS = {
    "Celite, n-butane":    (P["eq54_B_n_butane"], P["eq54_A_n_butane"], P["eq54_C_n_butane"]),
    "Celite, iso-butane":  (P["eq54_B_iso_butane"], P["eq54_A_iso_butane"], P["eq54_C_iso_butane"]),
    "Sterchamol, n-butane": (P["eq55_B_sterchamol"], 0.0, P["eq55_C_sterchamol"]),
}
SERIES = {
    "Celite, n-butane": (T4[T4.solid_carrier == "celite"].superficial_velocity_mm_s.to_numpy(),
                         T4[T4.solid_carrier == "celite"].hetp_n_butane_mm.to_numpy()),
    "Celite, iso-butane": (T4[T4.solid_carrier == "celite"].superficial_velocity_mm_s.to_numpy(),
                           T4[T4.solid_carrier == "celite"].hetp_iso_butane_mm.to_numpy()),
    "Sterchamol, n-butane": (T4[T4.solid_carrier == "sterchamol"].superficial_velocity_mm_s.to_numpy(),
                             T4[T4.solid_carrier == "sterchamol"].hetp_n_butane_mm.to_numpy()),
}
fitstat = {}
print(f"{'series':<22}{'N':>3}{'rmse':>8}{'bias':>8}{'max|r|':>8}{'null rmse':>11}"
      f"{'R^2':>8}{'u_opt':>8}{'H_min':>8}")
for name, (B, A, C) in FITS.items():
    u, h = SERIES[name]
    r = vd(u, B, A, C) - h
    rmse = float(np.sqrt((r ** 2).mean()))
    null = float(np.sqrt(((h - h.mean()) ** 2).mean()))
    uo = float(np.sqrt(B / C))                  # dH/du = -B/u^2 + C = 0
    u_root = one_root(lambda v: -B / v ** 2 + C, 1e-3, 1e3, f"dH/du = 0, {name}",
                      log=True, xtol=1e-15)     # the same root, found on the physical range
    Hm = vd(uo, B, A, C)
    fitstat[name] = dict(rmse=rmse, bias=float(r.mean()), maxr=float(np.abs(r).max()),
                         null=null, R2=1 - (r ** 2).sum() / ((h - h.mean()) ** 2).sum(),
                         u_opt=uo, H_min=float(Hm), N=len(u), u_root=u_root,
                         H_closed=float(A + 2 * np.sqrt(B * C)))
    f = fitstat[name]
    print(f"{name:<22}{f['N']:>3}{rmse:>8.3f}{f['bias']:>+8.3f}{f['maxr']:>8.3f}"
          f"{null:>11.3f}{f['R2']:>8.3f}{uo:>8.3f}{Hm:>8.3f}")

u17 = P["sterchamol_hetp_at_17"]
h17 = vd(17.0, *FITS["Sterchamol, n-butane"])
rows17 = T4[(T4.solid_carrier == "sterchamol") & (T4.superficial_velocity_mm_s == 17)]
print(f"\n   The p. 286 footnote says the Sterchamol H.E.T.P. at ~17 mm/s is 'about "
      f"{u17} mm'.")
print(f"   Eq (55) at u = 17 gives {h17:.3f} mm; Table 4's four rows at that velocity "
      f"average {rows17.hetp_n_butane_mm.mean():.3f} mm.")
print(f"   The footnote quotes the measurement, not the formula "
      f"({dev(h17, u17):+.0%} apart).")

ROOT_ID = max(abs(f["u_root"] / f["u_opt"] - 1) for f in fitstat.values())
print(f"\n  root-found vs closed-form optimum: max relative difference below 1e-12 "
      f"(here {ROOT_ID:.1e}; the same three constants both ways)")
print(f"  every fit beats the best constant H by a factor "
      f"{min(f['null']/f['rmse'] for f in fitstat.values()):.1f} to "
      f"{max(f['null']/f['rmse'] for f in fitstat.values()):.1f} in rmse")


# %% [markdown]
# #### The one term the paper says it could not determine
#
# Eq. (55) has no constant: "the absence of a constant in the right hand member is
# not a definite proof that eddy diffusion does not in principle occur … the
# constant could not be determined accurately. It should not differ much from
# zero, however."
#
# That is a checkable statement, and the Sterchamol rows are enough to check it.
# Refitting all three constants freely gives a value and a standard error for $A$.

# %%
# refit (least squares with standard errors) is in model.py.
us_, hs_ = SERIES["Sterchamol, n-butane"]
b3, se3, rms3 = refit(us_, hs_)
b2, se2, rms2 = refit(us_, hs_, ("B", "C"))
print("Sterchamol, refitted freely (this page's fit, not the paper's):")
for k in ("B", "A", "C"):
    print(f"   {k} = {b3[k]:+9.4f} +/- {se3[k]:.4f}   (t = {b3[k]/se3[k]:+.2f})")
print(f"   rmse {rms3:.4f} mm, against {rms2:.4f} mm with A forced to zero, "
      f"and {fitstat['Sterchamol, n-butane']['rmse']:.4f} mm for eq (55) as printed")
A_STER_T = float(b3["A"] / se3["A"])
print(f"\n  A = {b3['A']:+.3f} +/- {se3['A']:.3f} mm: consistent with zero "
      f"(|t| = {abs(A_STER_T):.2f}), so the paper's statement stands. It is also "
      f"consistent with the {P['eq54_A_n_butane']:.1f} mm of the Celite column, "
      f"which is the comparison the paper actually wants to make and cannot.")
print(f"  Eq (55) as printed sits {fitstat['Sterchamol, n-butane']['bias']:+.3f} mm "
      f"above the rows on average - a systematic offset a free constant absorbs.")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3))
ufine = np.linspace(1.0, 33.0, 400)
cols = {"Celite, n-butane": "tab:blue", "Celite, iso-butane": "tab:red",
        "Sterchamol, n-butane": "tab:green"}
mk = {"Celite, n-butane": "x", "Celite, iso-butane": "o", "Sterchamol, n-butane": "s"}
for name, (B, A, C) in FITS.items():
    u, h = SERIES[name]
    axes[0].plot(u, h, mk[name], color=cols[name], ms=6, mfc="none", label=f"{name} (Table 4)")
    axes[0].plot(ufine, vd(ufine, B, A, C), color=cols[name], lw=1.6)
    f = fitstat[name]
    axes[0].plot(f["u_opt"], f["H_min"], "*", color=cols[name], ms=13)
axes[0].set(xlabel="superficial gas velocity $u$  (mm s$^{-1}$)",
            ylabel="H.E.T.P.  $H$  (mm)", ylim=(0, 7),
            title="Table 4 and eqs. (54)/(55): a FIT, not a test")
axes[0].legend(fontsize=8)

B, A, C = FITS["Celite, n-butane"]
axes[1].plot(ufine, np.full_like(ufine, A), "--", color="0.5", label=f"$A = {A}$ mm (eddy)")
axes[1].plot(ufine, B / ufine, ":", color="tab:purple", label=f"$B/u$, $B = {B:.0f}$ mm$^2$s$^{{-1}}$")
axes[1].plot(ufine, C * ufine, "-.", color="tab:orange", label=f"$Cu$, $C = {C}$ s")
axes[1].plot(ufine, vd(ufine, B, A, C), "k", lw=2, label="sum")
axes[1].plot(fitstat["Celite, n-butane"]["u_opt"], fitstat["Celite, n-butane"]["H_min"],
             "k*", ms=13)
axes[1].set(xlabel="superficial gas velocity $u$  (mm s$^{-1}$)", ylabel="contribution (mm)",
            ylim=(0, 7), title="Celite / $n$-butane: the three terms separately")
axes[1].legend(fontsize=8)
plt.tight_layout(); plt.show()

# %% [markdown]
# ### 2. The film thickness, three times
#
# The paper inverts the $C$ term of eq. (53) to get the effective liquid film
# thickness, and reports it twice, "It is found from both curves that
# $d_f \approx 10\ \mu$" for the Celite column, and $d_f = 9\ \mu$ for Sterchamol.
# Since $KF_I/F_{II}$ is printed for all three cases and $D_{II}$ for the oil is
# printed once, each inversion can be redone from printed numbers only:
#
# $$d_f^2 = \frac{\pi^2}{8}\,\frac{C\,(1+x)^2\,D_{II}}{x},
# \qquad x \equiv K\frac{F_I}{F_{II}}.$$
#
# The three cases share $D_{II}$ but have three different $x$ and two different
# $C$, so they are three genuinely different arithmetic routes to the same
# physical quantity.

# %%
DII_FILM = P["D_II_liquid_film"] * 1e6            # m2/s -> mm2/s
CASES = [("Celite, n-butane", P["eq54_C_n_butane"], P["KFI_over_FII_n_butane"], P["d_f_celite"]),
         ("Celite, iso-butane", P["eq54_C_iso_butane"], P["KFI_over_FII_iso_butane"], P["d_f_celite"]),
         ("Sterchamol, n-butane", P["eq55_C_sterchamol"], P["KFI_over_FII_sterchamol"], P["d_f_sterchamol"])]

df_um = {}
print(f"{'case':<22}{'C (s)':>8}{'K F_I/F_II':>12}{'d_f (um)':>10}{'printed':>9}{'dev':>9}")
for name, C, x, printed in CASES:
    df2 = np.pi ** 2 / 8 * C * (1 + x) ** 2 * DII_FILM / x
    d = np.sqrt(df2) * 1e3                       # mm -> um
    df_um[name] = float(d)
    print(f"{name:<22}{C:>8.2f}{x:>12.3f}{d:>10.3f}{printed:>9.0f}{dev(d, printed):>+9.1%}")

DF_MAX_DEV = max(abs(dev(df_um[n], p)) for n, _, _, p in CASES)
print(f"\n  worst deviation from the printed value: {DF_MAX_DEV:.1%}")
print(f"  the two Celite determinations differ from each other by "
      f"{abs(dev(df_um['Celite, iso-butane'], df_um['Celite, n-butane'])):.1%} - "
      f"which is what 'found from both curves that d_f approx 10 mu' means, "
      f"quantified.")

# %% [markdown]
# ### 3. Eq. (56): where the resistance sits, and the paper's own check on it
#
# Having assumed the gas-phase resistance away, the paper checks it on p. 287 with
#
# $$\frac{1}{100}\left(\frac{d_p}{d_f}\right)^2\frac{D_{II}}{D_I},$$
#
# reporting the coefficient as $D_{II}/(100 D_I) = 3.3\times10^{-8}$ and the two
# ratios as $3\times10^{-7}$ and $2.5\times10^{-5}$. All three come out of printed
# numbers with no freedom at all, so they are the cleanest transcription check on
# the page, and the Sterchamol one is sharper when this page's own $d_f$ is used
# than when the paper's rounded 9 µm is.

# %%
coef = P["D_II_liquid_film"] / (100 * P["D_I_gilliland"])
print(f"D_II / (100 D_I) = {coef:.4e}   printed {P['D_II_over_100_D_I']:.1e}"
      f"   dev {dev(coef, P['D_II_over_100_D_I']):+.2%}")

R56 = {}
for label, dp_um, df_used, printed in (
        ("Celite,     d_f = 10 um printed", P["celite_dp"], P["d_f_celite"], P["ratio56_celite"]),
        ("Celite,     d_f from eq (53)", P["celite_dp"], df_um["Celite, n-butane"], P["ratio56_celite"]),
        ("Sterchamol, d_f = 9 um printed", P["sterchamol_dp"], P["d_f_sterchamol"], P["ratio56_sterchamol"]),
        ("Sterchamol, d_f from eq (53)", P["sterchamol_dp"], df_um["Sterchamol, n-butane"], P["ratio56_sterchamol"])):
    r = coef * (dp_um / df_used) ** 2
    R56[label] = float(r)
    print(f"  {label:<34} ratio = {r:.3e}   printed {printed:.1e}   dev {dev(r, printed):+7.1%}")

d_eq = np.sqrt(coef)
print(f"\n  film thickness at which the two resistances are EQUAL: "
      f"d_f/d_p = {d_eq:.3e}")
print("  What that tests is an ORDER-OF-MAGNITUDE sentence quoted to one significant")
print("  figure: 'the two resistances become of equal order of magnitude in the case")
print("  of a film thickness smaller than about 0 0001 of the particle diameter'")
print(f"  (p. 287), i.e. {P['df_over_dp_equal_resistance']:.0e}. Eq (56) set to unity gives {d_eq:.2e}, a factor")
print(f"  {d_eq/P['df_over_dp_equal_resistance']:.2f} larger - the SAME order of magnitude, which is all the sentence")
print("  claims. Recorded as a rounding of a one-significant-figure statement, not as")
print(f"  an error: the conclusion it carries (both columns far from that regime, at")
print(f"  {P['ratio56_celite']:.0e} and {P['ratio56_sterchamol']:.1e}) holds on either number.")
EQ56_DEV = float(abs(dev(R56["Celite,     d_f from eq (53)"], P["ratio56_celite"])))

# %% [markdown]
# ### 4. The slope ratio, which is the paper's own consistency test
#
# The paper offers exactly one test of eq. (53) against eq. (54) that does not go
# through a fitted constant. The $C$ term depends on the component only through
# $x = KF_I/F_{II}$, as $x/(1+x)^2$, so the ratio of the two straight-part slopes
# of Fig. 11 is fixed by the two printed values of $x$ alone, no $d_f$, no
# $D_{II}$, no velocity convention. The paper says this "must be about 1.5" and
# that eq. (54) gives 1.6.

# %%
shape = lambda x: x / (1 + x) ** 2
ratio_theory = shape(P["KFI_over_FII_iso_butane"]) / shape(P["KFI_over_FII_n_butane"])
ratio_fit = P["eq54_C_iso_butane"] / P["eq54_C_n_butane"]
print(f"from eq (53) and the two printed x:   {ratio_theory:.4f}   "
      f"(paper: 'about {P['slope_ratio_predicted']}')  dev {dev(ratio_theory, P['slope_ratio_predicted']):+.1%}")
print(f"from the two slopes of eq (54):       {ratio_fit:.4f}   "
      f"(paper: {P['slope_ratio_from_eq54']})           dev {dev(ratio_fit, P['slope_ratio_from_eq54']):+.1%}")
print(f"\n  theory vs fit: {dev(ratio_fit, ratio_theory):+.1%}")
print("  Note what this test is immune to: x is a ratio of volumes times K, so")
print("  the superficial/interstitial question below cannot touch it, and neither")
print("  can d_f or D_II. It is the strongest thing eq (53) says about Table 4.")
SLOPE_RATIO_DEV = float(abs(dev(ratio_theory, P["slope_ratio_predicted"])))

# %% [markdown]
# ### 5. The velocity convention, and what it does to $\gamma$
#
# The paper reads $2\gamma D_I$ off the $B$ coefficients of eq. (54) and gets
# $8\times10^{-6}$ m² s⁻¹, then checks it: "This value is of the right order of
# magnitude as $\gamma$ lies between 0 5 and 1, and $D_I = 9\times10^{-6}$ m²/sec
# as estimated from the well-known formula of GILLILAND".
#
# Two facts about that step, both arithmetic:
#
# * $8$ mm² s⁻¹ **is** $8\times10^{-6}$ m² s⁻¹. The quoted $2\gamma D_I$ is the
#   iso-butane $B$ coefficient carried over unchanged, i.e. with $u$ taken as the
#   velocity Table 4 tabulates, the **superficial** one, and no factor $F_I$.
# * The List of Symbols defines $u$ as **interstitial**. Under that reading,
#   $B_{\text{fitted}} = 2\gamma D_I F_I$, so $2\gamma D_I = B/F_I$ and
#   $\gamma = 0.5/F_I$ for $n$-butane.
#
# $F_I$ is not printed for these columns, so the ambiguity cannot be closed. What
# *can* be done is to state the band of $F_I$ each reading requires, using the
# paper's own bound $0.5 \le \gamma \le 1$.

# %%
DI = P["D_I_gilliland"]
print(f"{'reading':<44}{'2 gamma D_I':>14}{'gamma':>9}")
for lab, B in (("superficial (F_I = 1), n-butane B = 9", P["eq54_B_n_butane"]),
               ("superficial (F_I = 1), iso-butane B = 8", P["eq54_B_iso_butane"])):
    twogD = B * 1e-6
    print(f"  {lab:<42}{twogD:>14.2e}{twogD/(2*DI):>9.3f}")
GAMMA_SUP = P["eq54_B_n_butane"] * 1e-6 / (2 * DI)
GAMMA_ISO = P["eq54_B_iso_butane"] * 1e-6 / (2 * DI)
print(f"\n  the paper's quoted 2 gamma D_I = {P['two_gamma_D_I']:.0e} m2/s is the "
      f"iso-butane coefficient read this way, to {dev(P['eq54_B_iso_butane']*1e-6, P['two_gamma_D_I']):+.0%}.")
print(f"  It implies gamma = {GAMMA_ISO:.3f}, which is BELOW the paper's own "
      f"stated range 0.5-1.0 by {(P['gamma_lower']/GAMMA_ISO - 1)*100:.0f} %.")
print(f"  The n-butane coefficient gives gamma = {GAMMA_SUP:.3f}, exactly the lower bound.")

FI_lo = GAMMA_SUP / P["gamma_upper"]
FI_hi = GAMMA_SUP / P["gamma_lower"]
print(f"\n  Under the interstitial reading, gamma = {GAMMA_SUP:.1f}/F_I, so the paper's")
print(f"  own bound 0.5 <= gamma <= 1 requires  {FI_lo:.2f} <= F_I <= {FI_hi:.2f}.")
print("  Fig. 10 puts BOTH the intergranular space and the gas-filled pores in F_I,")
print("  so that band is physically reachable for these columns and the reading")
print("  cannot be excluded. What can be said is which way each constant moves, and")
print("  the two move in OPPOSITE directions:")
print(f"    gamma  = B/(2 D_I F_I), so the interstitial reading makes it LARGER by")
print(f"             1/F_I - up to {1/FI_lo:.1f}x.")
print(f"    d_f^2 ~ C, and the fitted C multiplies u_sup = F_I u_int, so eq (53)'s C")
print(f"             is F_I times the fitted one: d_f^2 is SMALLER by F_I, down by up")
print(f"             to {1/FI_lo:.1f}x, i.e. {1/np.sqrt(FI_lo):.2f}x in d_f itself.")
print("    x = K F_I/F_II is printed as one number, so it carries no ambiguity.")
print("  The paper applies the conversion in Fig. 9 (abscissa u F_I d_p) but not in")
print("  eqs (53)-(55).")
print("\n  Unaffected by any of this: 2 lambda d_p, because the eddy term does not")
print(f"  contain u at all. lambda = {P['two_lambda_dp']/(2*P['celite_dp']*1e-3):.2f} "
      f"from the printed 0.5 mm and 30 um, against the printed 'lambda approx 8'.")
LAMBDA_CELITE = float(P["two_lambda_dp"] / (2 * P["celite_dp"] * 1e-3))

# %% [markdown]
# ### 6. Table 3 refitted: the particle diffusivity and the packing factor
#
# This is the paper's *other* extraction, and it is the one that can be redone
# without any inferred input. From p. 282:
#
# > The particle diffusion coefficient $D_{II}$ may be calculated from the slope of
# > the straight line in Fig. 9 With $F_I = F_{II} = 0.4$ and $K = 1.3$ it is found
# > that $D_{II} = 1.3\times10^{-10}$ m²/sec., a value which is about 10 times
# > smaller than the diffusion coefficient in the free liquid.
#
# and
#
# > The values found for $\lambda$ are 8 for 200-400 mesh, 3 for 50-100 mesh and
# > practically zero for 20-40 mesh.
#
# Fig. 9 plots $H/d_p$ against $u F_I d_p$, that is, against **superficial**
# velocity times particle diameter, since $u$ there is interstitial. Table 3 holds
# the same numbers, so the regression can be done on the table and the figure need
# never be touched. Five rows: the 50-100 mesh series at five flow rates.
#
# Two things make this a real extraction rather than an algebraic rearrangement.
# The intercept gives $\lambda$ and the slope gives $D_{II}$, and they are
# independent; and, because the first term of $C$ is small, the slope pins
# $D_{II}$ **almost without needing $D_I$**, which the paper never states and never
# prints for this system.

# %%
FLOW = T3[T3.block == "flow"].sort_values("superficial_velocity_mm_s")
Hmid = (0.5 * (FLOW.hetp_lo_mm + FLOW.hetp_hi_mm)).to_numpy()
Ymid = (0.5 * (FLOW.h_over_dp_lo + FLOW.h_over_dp_hi)).to_numpy()
u_sup = FLOW.superficial_velocity_mm_s.to_numpy()
X = u_sup * DP_50_100                                   # = u_interstitial F_I d_p

FI = P["F_I_ion_exclusion"]; FII = P["F_II_ion_exclusion"]; KK = P["K_ion_exclusion"]
xK = KK * FI / FII

(two_lam, S), rms2t = fit_eq49(X, Ymid)
b3t, rms3t = fit_eq49(X, Ymid, three_term=True)
print("H/d_p = 2 lambda + S (u_sup d_p)      [molecular term dropped, as the paper does]")
print(f"   2 lambda = {two_lam:.4f}  ->  lambda = {two_lam/2:.4f}   "
      f"(printed {P['lambda_50_100']:.0f}, dev {dev(two_lam/2, P['lambda_50_100']):+.2%})")
print(f"   slope    S = {S:.3f} (mm^2/s)^-1        rmse {rms2t:.4f}")
print(f"\n   with the 1/X term restored: B' = {b3t[0]:+.5f}, 2 lambda = {b3t[1]:.3f}, "
      f"S = {b3t[2]:.3f}, rmse {rms3t:.4f}")
print(f"   B' is NEGATIVE and {abs(b3t[0])/ (two_lam*X.mean()) * 100:.2g} % of the "
      f"intercept term at the mean X - i.e. the five rows carry no molecular-diffusion")
print("   signal at all, which is what the paper asserts on p. 282 without measuring.")

a_c = FI ** 2 / (75 * (1 - FI) ** 2)
b_c = FI * KK / (2 * np.pi ** 2 * (1 - FI))
def DII_from_slope(S, D_I=None, visc=75.0):
    """Eq (49)'s C inverted for D_II (model.dii_from_slope) on the ion-exclusion constants."""
    return dii_from_slope(S, FI, FII, KK, D_I=D_I, visc=visc)

D_I_MM = 1.3e-3            # = 10 x D_II, the paper's own "about 10 times smaller"
print(f"\n   a_c = {a_c:.5f} (viscous term),  b_c = {b_c:.5f} (particle term)")
print(f"   {'D_I assumed (mm2/s)':<26}{'D_II (m2/s)':>14}{'dev from printed':>18}")
DII_ROWS = {}
for lab, D_I in (("1.3e-3 (= 10 x D_II)", D_I_MM), ("1.0e-3", 1.0e-3),
                 ("2.0e-3", 2.0e-3), ("infinite (term dropped)", None)):
    d2 = DII_from_slope(S, D_I) * 1e-6                  # mm2/s -> m2/s
    DII_ROWS[lab] = float(d2)
    print(f"   {lab:<26}{d2:>14.4e}{dev(d2, P['D_II_ion_exclusion']):>+18.2%}")
DII_BEST = DII_ROWS["1.3e-3 (= 10 x D_II)"]
DII_DEV = float(abs(dev(DII_BEST, P["D_II_ion_exclusion"])))
DII_SPREAD = float(max(DII_ROWS.values()) / min(DII_ROWS.values()) - 1)
print(f"\n   D_II = {DII_BEST:.3e} m2/s against the printed "
      f"{P['D_II_ion_exclusion']:.1e} m2/s: {DII_DEV:.1%}.")
print(f"   Across a 2x change in the assumed D_I - and even with the term removed")
print(f"   entirely - D_II moves by only {DII_SPREAD:.1%}. The extraction is")
print("   essentially free of the one input the paper does not print.")
LAMBDA_DEV = float(abs(dev(two_lam / 2, P["lambda_50_100"])))

# %%
fig, ax = plt.subplots(figsize=(6.4, 4.4))
xx = np.linspace(0, X.max() * 1.12, 200)
ax.plot(X, Ymid, "o", color="tab:blue", ms=7, label="Table 3, 50-100 mesh (midpoints)")
for x0, lo, hi in zip(X, FLOW.h_over_dp_lo, FLOW.h_over_dp_hi):
    ax.plot([x0, x0], [lo, hi], color="tab:blue", lw=1.4, alpha=0.7)
ax.plot(xx, two_lam + S * xx, "k", lw=1.8,
        label=f"fit: $2\\lambda = {two_lam:.2f}$, slope ${S:.0f}$")
ax.axhline(2 * P["lambda_50_100"], ls="--", color="tab:red", lw=1.3,
           label=f"$2\\lambda$ from the printed $\\lambda = {P['lambda_50_100']:.0f}$")
ax.set(xlabel="$u_{sup}\\,d_p = u F_I d_p$   (mm$^2$ s$^{-1}$)", ylabel="$H/d_p$",
       xlim=(0, xx.max()), ylim=(0, 36),
       title="Eq. (49) on Table 3: the vertical bars are the two printed determinations")
ax.legend(fontsize=8.5, loc="upper left"); plt.tight_layout(); plt.show()

# %% [markdown]
# ### 7. Table 1's `calc` column, and where $n = 650$ comes from
#
# Table 1's last column is the authors' own calculation, footnoted "*calculated
# for $n = 650$". The chain that produces it is short and entirely printed:
# $\Delta S_0 = 4v\sqrt{n}$ (eq. 24), $S_{rf} = nv + 0.31\Delta S_0$ (eq. 26), and
# $(c_{I,n}/c_0)_{\max} = \mathrm{erf}\!\left(a/2\sqrt2\right)$ with
# $a = A/(v\sqrt n)$ (eq. 20). Eq. (20) itself follows from eq. (19) in one line:
# the maximum of the band sits at the centre of the injected slug, so the integral
# runs symmetrically over $\pm a/2$.
#
# The whole `calc` column depends on **one** number, $v\sqrt n = \Delta S_0/4$,
# not on $n$ and $v$ separately. So the eleven printed values pin $\Delta S_0$, and
# that is a sharper constraint than the paper's own "$\Delta S_0 = \sim 165$ ml".

# %%
# eq (20) derived from eq (19): the maximum is at the slug centre.
xs_ = sp.symbols("x", real=True); a_s = sp.symbols("a", positive=True)
lhs = sp.integrate(sp.exp(-xs_ ** 2 / 2) / sp.sqrt(2 * sp.pi), (xs_, -a_s / 2, a_s / 2))
print(f"   eq (19) at the band maximum  ->  {sp.simplify(lhs)}   =  eq (20)\n")

calc = T1.dropna(subset=["cmax_calc"]).groupby("feed_volume_ml").cmax_calc.first()
A_vals, printed = calc.index.to_numpy(float), calc.to_numpy(float)

n_from_24 = (4 * P["nv_glc"] / P["delta_S0_glc"]) ** 2
print(f"   n from eq (24) with the printed nv and dS0: (4 nv / dS0)^2 = "
      f"(4 x {P['nv_glc']:.0f} / {P['delta_S0_glc']:.0f})^2 = {n_from_24:.1f}"
      f"   (paper prints {P['n_glc']:.0f}, dev {dev(n_from_24, P['n_glc']):+.1%})")
S_rf_mean = float(T1.final_retention_volume_ml.mean())
nv_from_26 = S_rf_mean - 0.31 * P["delta_S0_glc"]
print(f"   nv from eq (26) with the mean printed S_rf ({S_rf_mean:.1f} ml): "
      f"{nv_from_26:.1f} ml  (paper: ~{P['nv_glc']:.0f})")

def n_match(w):
    return int(sum(abs(round(cmax(A / w), 2) - c) < 1e-9 for A, c in zip(A_vals, printed)))

def calc_window(cells):
    """The interval of v sqrt(n) on which every printed calc value is reproduced.

    round(cmax(A/w), 2) = c holds for cmax(A/w) in [c - 0.005, c + 0.005], i.e. for
    w between A / a(c + 0.005) and A / a(c - 0.005), a(y) = 2 sqrt(2) erfinv(y)
    (eq 20 inverted). The window is the intersection over the rows, in closed form.
    """
    a_of = lambda y: 2 * SQ2 * erfinv(y)
    lo = max(A / a_of(c + 0.005) for A, c in zip(A_vals, cells))
    hi = min(A / a_of(c - 0.005) for A, c in zip(A_vals, cells))
    return float(lo), float(hi)

W_LO, W_HI = calc_window(printed)
if not W_LO < W_HI:
    raise RuntimeError(f"no v sqrt(n) reproduces every printed calc value: [{W_LO}, {W_HI}]")
_inner = (W_LO * (1 + 1e-9), W_HI * (1 - 1e-9))
_outer = (W_LO * (1 - 1e-6), W_HI * (1 + 1e-6))
if min(n_match(w) for w in _inner) != len(A_vals) or max(n_match(w) for w in _outer) == len(A_vals):
    raise RuntimeError("the closed-form window edges do not bound the reproduced rows")
print(f"\n   {'v sqrt(n)':>10}{'= dS0/4 from':>26}{'rows of 11 reproduced':>24}")
for w, lab in ((P["delta_S0_glc"] / 4, "the printed dS0 = 165"),
               (P["nv_glc"] / P["n_glc"] * np.sqrt(P["n_glc"]), "the printed nv and n"),
               (0.5 * (W_LO + W_HI), "this page's window, centre")):
    print(f"   {w:>10.3f}{lab:>26}{n_match(w):>18} / {len(A_vals)}")
print(f"\n   All {len(A_vals)} printed calc values are reproduced for "
      f"v sqrt(n) in [{W_LO:.3f}, {W_HI:.3f}],")
print(f"   i.e. dS0 in [{4*W_LO:.1f}, {4*W_HI:.1f}] ml - {dev(4*W_LO, P['delta_S0_glc']):+.1%} to "
      f"{dev(4*W_HI, P['delta_S0_glc']):+.1%} on the paper's own '~165 ml'.")
print("   Nothing here is a defect: the paper writes both 165 and 1055 with a tilde.")
print("   What the calc column pins is dS0 alone, and it pins it to 0.9 %.")
CALC_N_OK = n_match(0.5 * (W_LO + W_HI))
CALC_FRAC = CALC_N_OK / len(A_vals)

print(f"\n   {'A (ml)':>8}{'printed':>9}{'this page':>11}{'a':>8}")
w0 = 0.5 * (W_LO + W_HI)
for A, c in zip(A_vals, printed):
    print(f"   {A:>8.0f}{c:>9.2f}{cmax(A/w0):>11.4f}{A/w0:>8.3f}")

# %% [markdown]
# ### 8. Table 2 → Table 3, both routes, and eq. (21) as printed cannot be right
#
# The paper says the plate number is obtainable "in two independent ways" (p. 280)
# - from the **height** of the elution curve through eq. (20) and Fig. 3, and from
# its **width** through eqs. (21)/(22) and Fig. 4. Table 3's printed ranges are
# those two determinations. With Table 2 and $nv = 70$ ml, both can be redone.
#
# The height route needs only eq. (20). The width route needs eq. (21), and
# **eq. (21) as printed does not reproduce the paper's own eqs. (24) and (25).**
# It reads
#
# $$\Delta s = a + 2\delta + \sqrt{2\pi}\;\frac{a+\delta}{\delta}\;e^{\frac12\delta^2}
#   \left(\mathrm{erf}\frac{a+\delta}{\sqrt2} - \mathrm{erf}\frac{\delta}{\sqrt2}\right)$$
#
# with $\delta$ from eq. (22). Pin what is *not* free first:
#
# * **Eq. (22) is not free.** It is exactly the condition that the second
#   derivative of eq. (19) vanish, the inflection point, which is what $\Delta S$
#   is defined from (p. 275, Fig. 2). Verified symbolically below.
# * **The $a\to0$ limit is not free.** Eq. (24) states $\Delta S_0 = 4v\sqrt n$,
#   i.e. $\Delta s \to 4$, and the paper says why: "For the Gaussian curve this
#   value of $\Delta S$ is just equal to 4 times the standard deviation".
# * **The $a\to\infty$ limit is not free.** Eq. (25) states
#   $\Delta S = A + v\sqrt{2\pi n}$ for $a>3$, i.e. $\Delta s \to a + \sqrt{2\pi}$.
#
# The only thing left free is the denominator of the third term, and both limits
# force it to be $a$, not $\delta$. Deriving $\Delta s$ from eq. (19) directly
# gives $a$; and $(a+\delta)/a$ is what makes the $a\to0$ limit finite at all,
# because the bracket vanishes linearly in $a$.

# %%
# --- eq (22) IS the inflection condition of eq (19) -------------------------
xv, av, dv = sp.symbols("x a delta", positive=True)
xr = sp.symbols("xr", real=True)
phi = sp.exp(-xr ** 2 / 2) / sp.sqrt(2 * sp.pi)
f19 = sp.integrate(phi, (xr, sp.Symbol("y", real=True) - av, sp.Symbol("y", real=True)))
y = sp.Symbol("y", real=True)
f19 = (sp.erf(y / sp.sqrt(2)) - sp.erf((y - av) / sp.sqrt(2))) / 2
infl = sp.simplify(sp.diff(f19, y, 2).subs(y, av + dv))
eq22 = dv * sp.exp(-dv ** 2 / 2) - (av + dv) * sp.exp(-(av + dv) ** 2 / 2)
ratio22 = sp.simplify(sp.simplify(infl) / eq22)
print(f"   d2(eq 19)/dy2 at y = a + delta, divided by eq (22): {sp.simplify(ratio22)}")
print("   -> eq (22) is exactly the inflection condition; it is not free.\n")

# --- and the corrected denominator is what eq (19) gives, not a guess ---------
# Delta s = 2 (y_R - a/2) + 2 f(y_R)/|f'(y_R)| with y_R = a + delta the right
# inflection point.  Write E = exp((a+delta)^2/2) and e = exp(delta^2/2); then
# f' = (1/E - 1/e)/sqrt(2 pi) and eq (22) reads delta E = (a + delta) e.
E_, e_, Q_ = sp.symbols("E e Q", positive=True)     # Q = the erf difference
tail_raw = 2 * (Q_ / 2) / ((1 / e_ - 1 / E_) / sp.sqrt(2 * sp.pi))
tail = sp.simplify(tail_raw.subs(E_, (av + dv) * e_ / dv))          # apply eq (22)
tail_with_a = sp.sqrt(2 * sp.pi) * (av + dv) / av * e_ * Q_
tail_with_delta = sp.sqrt(2 * sp.pi) * (av + dv) / dv * e_ * Q_
r21 = sp.simplify(tail - tail_with_a)
r21_printed = sp.simplify(tail - tail_with_delta)
print(f"   2 f / |f'| at the inflection, after eq (22):  {sp.simplify(tail)}")
print(f"     minus the third term with 'a'      underneath: {r21}")
print(f"     minus the third term with 'delta'  underneath: {sp.factor(r21_printed)}")
print("   -> the denominator required by eq (19) is a. Not a guess, and not")
print("      a choice made to fit the limits: the limits merely confirm it.")
print("")


# delta_of (eq 22, the unique root in (0, 1)) and width_red (eq 21) are in model.py.
print(f"   {'a':>10}{'as printed':>16}{'with a':>10}{'eq (24)/(25) requires':>24}")
for a, req in ((1e-8, 4.0), (0.01, 4.0), (1.0, None), (3.0, 3 + np.sqrt(2*np.pi)),
               (5.0, 5 + np.sqrt(2*np.pi)), (10.0, 10 + np.sqrt(2*np.pi))):
    pr, co = width_red(a, "delta"), width_red(a, "a")
    r = "-" if req is None else f"{req:.4f}"
    print(f"   {a:>10.2g}{pr:>16.4g}{co:>10.4f}{r:>24}")
EQ21_A0_PRINTED = float(width_red(1e-8, "delta"))
EQ21_A0_CORRECTED = float(width_red(1e-8, "a"))
EQ21_LARGE_A_DEV = float(abs(dev(width_red(5.0, "a"), 5 + np.sqrt(2 * np.pi))))
print(f"\n   As printed:  a -> 0 gives {EQ21_A0_PRINTED:.4f}, not the 4 of eq (24);")
print(f"                a -> large DIVERGES instead of approaching a + sqrt(2 pi).")
print(f"   With a:      both limits are met, the second to {EQ21_LARGE_A_DEV:.1e} at a = 5.")
print("   REPORTED, NOT REPAIRED. The corrected form is used for the width route")
print("   below and every number that depends on it is labelled; nothing else on")
print("   this page uses eq (21) at all.")
EQ21_DERIVED_OK = int(r21 == 0)

# %%
NV = P["nv_ion_exclusion"]

# n_from_height (eq 20) and n_from_width (eq 21, one root on [A/25, 1e4] ml) are in
# model.py; both take nv.
single = T2[T2.feed_volume_ml == 5]
print("Table 3's plate numbers, recomputed from Table 2 (A = 5 ml rows only)")
print(f"{'mesh':>8}{'u':>7}{'from height':>13}{'from width':>12} |{'printed lo':>11}"
      f"{'printed hi':>11}{'dev lo':>9}{'dev hi':>9}")


def table3_from_table2(nv):
    """Both routes for every A = 5 ml row of Table 2, against Table 3's printed pair.

    Returns (mesh, u, n_height, n_width, dev_lo, dev_hi, printed_lo, printed_hi)
    per row, the printed pair being matched by size (it is printed ascending).
    """
    out = []
    for _, r in single.iterrows():
        nh = n_from_height(r.feed_volume_ml, r.max_height_frac_c0, nv)
        nw = n_from_width(r.feed_volume_ml, r.width_at_base_ml, nv)
        m = T3[(T3.particle_mesh == r.particle_mesh)
               & (np.isclose(T3.superficial_velocity_mm_s, r.superficial_velocity_mm_s))].iloc[0]
        lo, hi = min(nh, nw), max(nh, nw)
        out.append((r.particle_mesh, r.superficial_velocity_mm_s, nh, nw,
                    dev(lo, m.n_plates_lo), dev(hi, m.n_plates_hi), m.n_plates_lo, m.n_plates_hi))
    return out


def route_medians(recs):
    """Median |dev| by route: the height route takes the endpoint its size gives it."""
    dh, dw = [], []
    for mesh, uu, nh, nw, dlo, dhi, *_ in recs:
        a, b = (dlo, dhi) if nh <= nw else (dhi, dlo)
        dh.append(abs(a)); dw.append(abs(b))
    return float(np.median(dh)), float(np.median(dw))


recs = table3_from_table2(NV)
for mesh, uu, nh, nw, dlo, dhi, plo, phi in recs:
    print(f"{mesh:>8}{uu:>7.3f}{nh:>13.0f}{nw:>12.0f} |"
          f"{plo:>11.0f}{phi:>11.0f}{dlo:>+9.1%}{dhi:>+9.1%}")

devs = np.abs([d for *_, dlo, dhi, _, _ in recs for d in (dlo, dhi)])
N_MED = float(np.median(devs)); N_MAX = float(devs.max())
print(f"\n   median |dev| {N_MED:.1%}, worst {N_MAX:.1%} over "
      f"{len(devs)} endpoints of {len(recs)} rows")

# Split by ROUTE, which is not the same as splitting by endpoint. The printed
# pair is in ascending order and the paper does not say which determination is
# which; Fig. 8's own printed labels settle it for one condition, and the same
# ascending convention assigns the rest. In most rows the height route gives the
# smaller n - but not in all of them, and where it inverts, so does the pairing.
fig8 = T3[(T3.particle_mesh == "50-100") & (T3.block == "particle")].iloc[0]
print(f"\n   Fig. 8's printed curve labels are n = {P['n_fig8_height']:.0f} (height) and "
      f"{P['n_fig8_width']:.0f} (width), which are exactly")
print(f"   this table's {fig8.n_plates_lo:.0f}-{fig8.n_plates_hi:.0f} for 50-100 mesh at "
      f"{fig8.superficial_velocity_mm_s} mm/s: low end = height, high end = width there.")
inverted = [f"{mesh} at {uu} mm/s" for mesh, uu, nh, nw, *_ in recs if nh > nw]
N_MED_H, N_MED_W = route_medians(recs)
print(f"   BY ROUTE, the height route is the accurate one: median {N_MED_H:.1%} against "
      f"{N_MED_W:.1%} for the width route.")
print(f"   The two routes invert in {len(inverted)} of {len(recs)} rows ({', '.join(inverted)}),")
print("   where the height determination is therefore the printed HIGH endpoint;")
print("   sorting by size instead of by route would report "
      f"{np.median([abs(r[4]) for r in recs]):.1%} and {np.median([abs(r[5]) for r in recs]):.1%}.")
print("   That asymmetry is expected: Table 2's widths are printed as integers, so")
print(f"   a +-0.5 ml read on the 200-400 mesh row's 10 ml is +-5 %, and n goes as")
print("   1/w^2. The height column carries two decimals.")
print("\n   The 50-100 / 0.091 row is the worst of the nine, and it is the one")
print("   condition carrying a five-point feed-volume sweep: the paper fitted n to")
print("   the WHOLE sweep in Fig. 8, not to the A = 5 ml pair the table above uses.")
print("   Over that sweep the height route alone gives:")
for _, r in T2[(T2.particle_mesh == "50-100")
               & (np.isclose(T2.superficial_velocity_mm_s, 0.091))].iterrows():
    nh = n_from_height(r.feed_volume_ml, r.max_height_frac_c0, NV)
    print(f"      A = {r.feed_volume_ml:>4.0f} ml, height {r.max_height_frac_c0:.2f}"
          f"  ->  n = {nh:>7.0f}" + ("   (saturated, no information)" if np.isnan(nh) else ""))
sw_row = T3[(T3.particle_mesh == "50-100") & (T3.block == "flow")
            & np.isclose(T3.superficial_velocity_mm_s, 0.091)].iloc[0]
print(f"   which brackets the printed {sw_row.n_plates_lo:.0f}-{sw_row.n_plates_hi:.0f}.")

# %% [markdown]
# ## Validation
#
# Eq. (38) is the hinge of the whole paper: everything after it, eq. (49),
# eq. (53), $\lambda$, $D_{II}$, $d_f$, is eq. (38) with substitutions. So the
# validation that matters is one that produces $H$ **without using eq. (38)**.
#
# The route taken here is to integrate eqs. (29)/(30) in time with pymrm, record
# $c_I(t)$ at two interior stations, and form
#
# $$H = (z_2-z_1)\,\frac{\sigma^2_{t}(z_2)-\sigma^2_{t}(z_1)}
#                       {\left[\bar t(z_2)-\bar t(z_1)\right]^2}.$$
#
# Taking **differences between two stations** is deliberate. For a linear system
# the cumulants of the response are exactly additive along the column, so every
# entrance effect, every exit effect and the entire shape of the injected pulse are
# constants that cancel in the difference. The number that comes out is the
# asymptotic plate height and nothing else, which is what eq. (38) claims to be.
#
# ### Conditions
#
# The ion-exclusion column, because it is the one whose $F_I$, $F_{II}$ and $K$ are
# printed. Everything else is printed too, except two quantities the paper itself
# declares immaterial; their contribution is measured rather than assumed.

# %%
DP = DP_50_100                                  # recovered from Table 3, an INFERENCE
U_SUP = 0.091                                   # printed, Tables 2 and 3
U_INT = U_SUP / FI                              # eqs (29)/(30) use interstitial u
DII_M = P["D_II_ion_exclusion"] * 1e6           # m2/s -> mm2/s, printed
D_I_ASSUMED = 10 * DII_M                        # "about 10 times smaller than in the free liquid"
GAMMA_ASSUMED = 0.5 * (P["gamma_lower"] + P["gamma_upper"])   # midpoint of the printed range
LAM = P["lambda_50_100"]                        # printed

D_LONG = GAMMA_ASSUMED * D_I_ASSUMED + LAM * U_INT * DP       # eqs (47),(48)
a_p = 6 * (1 - FI) / DP                                        # eq (46)
INV_ALPHA = (6 / 25) * FI / (D_I_ASSUMED * a_p ** 2) + (3 / (2 * np.pi ** 2)) * KK * DP / (DII_M * a_p)
ALPHA = 1 / INV_ALPHA                                          # eq (45)
BETA = 1 / (1 + FII / (FI * KK))                               # eq (34)

H38 = 2 * D_LONG / U_INT + 2 * (U_INT * FI * INV_ALPHA) / (1 + xK) ** 2
MT_SHARE = 2 * (U_INT * FI * INV_ALPHA) / (1 + xK) ** 2 / H38   # the transfer term's share of H
mol_share = 2 * GAMMA_ASSUMED * D_I_ASSUMED / U_INT / H38
visc_share = (6 / 25) * FI / (D_I_ASSUMED * a_p ** 2) / INV_ALPHA
VISC_SHARE_H = visc_share * MT_SHARE                            # ... and OF H, which is smaller

print(f"   d_p        = {DP:.4f} mm          (RECOVERED from Table 3, an inference)")
print(f"   u          = {U_INT:.4f} mm/s        (= {U_SUP} superficial / F_I)")
print(f"   D          = {D_LONG:.5f} mm2/s      (eq 48; molecular part is "
      f"{mol_share:.2%} of H)")
print(f"   1/alpha    = {INV_ALPHA:.3f} s           (eq 45; the D_I term is "
      f"{visc_share:.2%} of IT, and {VISC_SHARE_H:.3%} of H)")
print(f"   beta       = {BETA:.5f}            band velocity {BETA*U_INT:.5f} mm/s")
print(f"\n   H from eq (38) = {H38:.5f} mm")
ROW_091 = T3[(T3.particle_mesh == "50-100") & (T3.block == "flow")
             & np.isclose(T3.superficial_velocity_mm_s, U_SUP)].iloc[0]
print(f"   Table 3 prints {ROW_091.hetp_lo_mm}-{ROW_091.hetp_hi_mm} mm for this row "
      f"(50-100 mesh, u = {U_SUP} mm/s), i.e. n = {ROW_091.n_plates_lo:.0f}-{ROW_091.n_plates_hi:.0f}.")
print("   That agreement is NOT evidence: lambda and D_II were extracted FROM those")
print("   very rows in section 6. It is a closure check on the extraction only.")
print(f"\n   H splits {1-MT_SHARE:.2%} longitudinal / {MT_SHARE:.2%} mass transfer, so the")
print(f"   two assumed inputs carry {mol_share:.2%} of H (gamma D_I, directly) and")
print(f"   {VISC_SHARE_H:.3%} of H (the D_I term of 1/alpha, which is {visc_share:.2%} of 1/alpha")
print("   but only its transfer share of H). 'Longitudinal molecular diffusion was of")
print("   no importance in SIMPSON and WHEATON's experiments' (p. 282) is confirmed,")
print("   quantitatively, and by a wider margin than the 1/alpha share alone suggests.")

# %%
L_TOT, Z1, Z2, T_END, T_SIG = 700.0, 200.0, 600.0, 7000.0, 12.0

def solve(n_z=700, dt=2.0, limiter=vanleer, t_end=T_END, nu=0,
          alpha=ALPHA, K=KK, D=D_LONG, bc_left=None, sweeps=1, L=L_TOT,
          stations=(Z1, Z2), theta=0.5, newton=False, inlet_data=True):
    """One run of eqs (29)/(30): records at the stations and their cumulants.

    The step is linear except for the van Leer limiter. The default marches it
    as a deferred correction, one linear solve per step (Column.run); newton=True
    solves the same theta scheme with the exact Jacobian of ConvFluxTVD
    (Column.run_newton), which is what the deferred route is checked against.
    """
    col = Column(U_INT, D, alpha, K, FI, FII, L, n_z, limiter=limiter,
                 nu=nu, bc_left=bc_left)
    if newton:
        ts, rec, zs = col.run_newton(t_end, dt, T_SIG, list(stations), theta=theta)
    else:
        ts, rec, zs = col.run(t_end, dt, T_SIG, list(stations), sweeps=sweeps, theta=theta,
                            inlet_data=inlet_data)
    return ts, rec, zs, [cumulants(ts, rec[:, j]) for j in range(len(stations))]

ts_r, rec_r, zs_r, Kc_r = solve(n_z=1400, dt=1.0)
H_PYMRM = float(plate_height(Kc_r, zs_r))
H_REL = float(dev(H_PYMRM, H38))
print(f"   stations at z = {zs_r[0]:.4f} and {zs_r[1]:.4f} mm (cell centres)")
print(f"   H from the pymrm solution   = {H_PYMRM:.6f} mm")
print(f"   H from eq (38)              = {H38:.6f} mm")
print(f"   relative difference         = {H_REL:+.3e}")
print("\n   Nothing in this route uses eq (38), the Gaussian (33), or the plate")
print("   theory. It uses the balance (29)/(30) and two definitions.")

# %% [markdown]
# ### Refinement, on both axes, and a coefficient that is known in advance
#
# Bare first-order upwind adds a numerical dispersion of exactly $u h/2$, and
# eq. (38)'s first term is $2D/u$, so switching the limiter off must shift $H$ by
# **exactly $h$**, not "by order $h$". That turns the grid study into a
# prediction with a known answer, which is a stronger check on the assembly than
# an observed order alone. Orders come from successive differences of four
# nested grids (`observed_orders`), so no converged reference is needed.
#
# The time step is refined separately. The two-station $H$ is nearly
# *insensitive* to it: the two stations are sampled identically, so the leading
# time error cancels in the difference. The step is linear apart from the van
# Leer limiter, and the march treats the limiter as a deferred correction: one
# linear solve per step, the correction taken from the state at the start of the
# step. That lag makes the scheme formally first order in time. The same theta
# scheme with the limited flux solved exactly at every step (`ConvFluxTVD` with
# its Jacobian, `newton`) is the check on what the lag costs.
#
# The time order is measured on the linear scheme, the limiter switched off,
# where Crank-Nicolson's second order should show cleanly on $c_I$ at $z_2$ at two
# fixed times on the flanks of the band, which carry the time error pointwise.
# With the limiter on, the same ladder is printed: first order at t = 4800 s, irregular
# at 5200 s, and the page claims no order for it.

# %%
KAPZ = {}                    # (n_z, dt) -> (kappa_3/z, kappa_4/z), kept for later
NZS = (350, 700, 1400, 2800)

print("Grid, van Leer deferred correction (dt = 0.5 s throughout)")
GRID = {}
for n_z in NZS:
    _, _, zs, Kc = solve(n_z=n_z, dt=0.5)
    KAPZ[(n_z, 0.5)] = kappas_per_z(Kc, zs)
    GRID[n_z] = float(plate_height(Kc, zs))
    print(f"   n_z {n_z:>5d}  h = {L_TOT / n_z:.4f} mm   H = {GRID[n_z]:.7f}   "
          f"|dev| from eq (38) {abs(dev(GRID[n_z], H38)):.3e}")
OBS_GRID = observed_orders(lambda n: GRID[n], NZS)
ORD_GRID = float(max(OBS_GRID["orders"], key=lambda o: abs(o - 3.0)))
GRID_ROUNDOFF = float(min(abs(GRID[a] - GRID[b]) for a, b in zip(NZS[:-1], NZS[1:]))
                      / (int(round(T_END / 0.5)) * EPS * H38))
print(f"   observed orders (successive differences): "
      + ", ".join(f"{o:.3f}" for o in OBS_GRID["orders"]))
print(f"   smallest successive difference over the round-off floor (n_steps eps H): "
      f"above 1e3 (here {GRID_ROUNDOFF:.1e})")

print("\nGrid, BARE UPWIND - predicted H_num = H + h exactly")
up_rows = []
for n_z in NZS[:3]:
    _, _, zs, Kc = solve(n_z=n_z, dt=0.5, limiter=None)
    H = float(plate_height(Kc, zs)); h = L_TOT / n_z
    up_rows.append((h, (H - H38) / h, H))
    print(f"   n_z {n_z:>5d}  h = {h:.4f} mm   H = {H:.6f}   (H - H_38)/h = {(H-H38)/h:.5f}")
ORD_UPWIND = float(observed_orders(lambda n: dict(zip(NZS[:3], [r[2] for r in up_rows]))[n],
                                   NZS[:3])["orders"][0])
UPWIND_COEF = float(np.mean([c for _, c, _ in up_rows]))
print(f"   observed order {ORD_UPWIND:.2f} (first order, against {ORD_GRID:.2f} with the limiter)")
UPWIND_COEF_DEV = float(abs(UPWIND_COEF - 1.0))
print(f"   mean coefficient {UPWIND_COEF:.6f} against the predicted 1: "
      f"|error| below 1e-4 (here {UPWIND_COEF_DEV:.1e})")

# Time: the ladder is given as step counts, so observed_orders sees a ratio of 2.
# The observable is c_I at z2 at two fixed times on the falling and rising flanks of
# the band (whole steps on every ladder), which carries the time error pointwise.
DTS = (8.0, 4.0, 2.0, 1.0, 0.5)
STEPS = tuple(int(round(T_END / dt)) for dt in DTS)
T_OBS = (4800.0, 5200.0)


def dt_ladder(limiter=vanleer, theta=0.5, n_z=700, inlet_data=True):
    """H, mu1(z2) and c_I(z2, T_OBS) for every step of DTS, keyed by the step count."""
    out = {}
    for dt, m in zip(DTS, STEPS):
        ts, rec, zs, Kc = solve(n_z=n_z, dt=dt, limiter=limiter, theta=theta,
                               inlet_data=inlet_data)
        if limiter is vanleer and theta == 0.5 and inlet_data:
            KAPZ[(n_z, dt)] = kappas_per_z(Kc, zs)
        cz = np.array([rec[int(round(t / dt)), 1] for t in T_OBS])
        out[m] = (float(plate_height(Kc, zs)), float(Kc[1][1]), cz)
    return out


def dt_study(ladder):
    """Orders of c_I(z2, T_OBS), worst departure from 2, round-off margin, H swing."""
    orders = np.array([observed_orders(lambda m, j=j: ladder[m][2][j], STEPS)["orders"]
                       for j in range(len(T_OBS))]).T        # (triples, times)
    H_ = np.array([ladder[m][0] for m in STEPS])
    fin = float(np.min(np.abs(ladder[STEPS[-1]][2] - ladder[STEPS[-2]][2])
                       / (STEPS[-1] * EPS * np.abs(ladder[STEPS[-1]][2]))))
    return orders, float(np.max(np.abs(orders - 2.0))), fin, float(H_.max() / H_.min() - 1)


print(f"\nTime step (n_z = 700), dt = {', '.join(f'{d:g}' for d in DTS)} s")
LAD_VL = dt_ladder()
LAD_UP = dt_ladder(limiter=None)
O_VL, _, _, DT_SWING_H = dt_study(LAD_VL)
O_UP, DT_ORDER_DEV, DT_ROUNDOFF, DT_SWING_H_UP = dt_study(LAD_UP)
print(f"   {'dt (s)':>7}{'H, van Leer':>14}{'mu1(z2), van Leer':>20}"
      + "".join(f"{f'c_I(z2, {t:.0f} s)':>20}" for t in T_OBS) + "   (no limiter)")
for dt, m in zip(DTS, STEPS):
    print(f"   {dt:>7g}{LAD_VL[m][0]:>14.8f}{LAD_VL[m][1]:>20.6f}"
          + "".join(f"{x:>20.10f}" for x in LAD_UP[m][2]))
fmt_o = lambda o: "; ".join(", ".join(f"{x:.3f}" for x in col) for col in np.asarray(o).T)
def mu1_sign_changes(ladder):
    """Sign changes in the successive differences of mu1(z2) along the dt ladder."""
    d = np.sign(np.diff([ladder[m][1] for m in STEPS]))
    return int(np.sum(d[1:] != d[:-1]))


MU1_SIGN_CHANGES = mu1_sign_changes(LAD_VL)
print(f"\n   orders, no limiter (linear Crank-Nicolson), at t = {T_OBS[0]:.0f} s; {T_OBS[1]:.0f} s: "
      + fmt_o(O_UP))
print(f"   worst departure from 2: below 0.1 (here {DT_ORDER_DEV:.3f}); smallest finest difference "
      f"over the round-off floor (n_steps eps c): above 1e3 (here {DT_ROUNDOFF:.1e})")
print(f"   orders with the van Leer correction: " + fmt_o(O_VL)
      + f"  (mu1(z2) differences change sign {MU1_SIGN_CHANGES} time(s) on this ladder)"
      + "  (first order at t = 4800 s, irregular at 5200 s; none is claimed)")
print(f"   over the 16x change in dt the two-station H moves below 1e-6 relative "
      f"(here {DT_SWING_H:.1e}): it is grid-limited, not dt-limited.")

# The order is a claim about the SCHEME, so break the scheme: theta = 1 is
# backward Euler, first order in time, and nothing else changes.
LAD_BE = dt_ladder(limiter=None, theta=1.0)
O_BE, DT_ORDER_DEV_BE, _, _ = dt_study(LAD_BE)
LAD_BE_VL = dt_ladder(theta=1.0)
DT_SWING_H_BE = dt_study(LAD_BE_VL)[3]
print("\nSame ladders with BACKWARD EULER (theta = 1) in place of Crank-Nicolson:")
print(f"   orders, no limiter: " + fmt_o(O_BE) + "  (first order)")
print(f"   H with van Leer swings {DT_SWING_H_BE:.1e} over the same ladder, against "
      f"{DT_SWING_H:.1e} with Crank-Nicolson.")

# What the lag costs: the same theta scheme with the limited flux solved exactly.
N_NEWTON, DT_NEWTON = 350, 2.0
_, _, zs_n, Kc_n = solve(n_z=N_NEWTON, dt=DT_NEWTON, newton=True)
_, _, zs_d, Kc_d = solve(n_z=N_NEWTON, dt=DT_NEWTON)
H_NEWTON, H_DEFER = float(plate_height(Kc_n, zs_n)), float(plate_height(Kc_d, zs_d))
DEFER_VS_NEWTON = float(abs(dev(H_DEFER, H_NEWTON)))
MU1_DEFER_VS_NEWTON = float(abs(dev(Kc_d[1][1], Kc_n[1][1])))
MU1Z1_DEFER_VS_NEWTON = float(abs(dev(Kc_d[0][1], Kc_n[0][1])))
print(f"\nDeferred correction against ConvFluxTVD + newton (n_z = {N_NEWTON}, dt = {DT_NEWTON:g} s):")
print(f"   H {H_DEFER:.9f} against {H_NEWTON:.9f}: below 1e-5 relative (here {DEFER_VS_NEWTON:.1e});")
print(f"   the same grid is {abs(dev(H_NEWTON, H38)):.1e} from eq (38), so the lag is not what limits H")
print(f"   mu1(z1) {Kc_d[0][1]:.6f} against {Kc_n[0][1]:.6f} s: {MU1Z1_DEFER_VS_NEWTON:.1e} relative")
print(f"   mu1(z2) {Kc_d[1][1]:.6f} against {Kc_n[1][1]:.6f} s: {MU1_DEFER_VS_NEWTON:.1e} relative")
print(f"   The lag moves the first moment at both stations; the two-station H is closer"
      f" than either ({DEFER_VS_NEWTON:.1e} against {MU1Z1_DEFER_VS_NEWTON:.1e} and"
      f" {MU1_DEFER_VS_NEWTON:.1e}).")

# %% [markdown]
# ### Break table
#
# Every number this page reports needs something that moves it. The last row is
# a deliberate *non*-mover, included because the insensitivity is the design claim
# being made: the inlet condition moves the two-station $H$ by less than $10^{-5}$
# relative. That claim is a bound, not a value, because it is a difference of two
# nearly equal plate heights; what breaks it is reading $H$ from ONE station, where
# the same change of inlet condition is no longer cancelled.

# %%
brk = []
def add(label, value, note=""):
    """Every row is judged against the BASELINE run, which is what 'does it move'
    means; the deviation from eq (38) is shown too, because that is the metric."""
    base = brk[0][1] if brk else value
    brk.append((label, value, abs(dev(value, base)), abs(dev(value, H38)), note))

_, _, zsb, Kcb = solve(n_z=700, dt=2.0)
add("baseline: H from pymrm (n_z=700, dt=2, van Leer)", float(plate_height(Kcb, zsb)))
BRK1_KW = {"bare upwind (limiter=None)": dict(limiter=None),
           "coarse grid, n_z = 175": dict(n_z=175),
           "F_II halved (a different phase ratio)": dict(),
           "alpha halved": dict(alpha=ALPHA / 2),
           "K perturbed by +1 %": dict(K=KK * 1.01),
           "D perturbed by +1 %": dict(D=D_LONG * 1.01),
           "time window truncated to 60 %": dict(t_end=0.6 * T_END)}
for lab, kw in BRK1_KW.items():
    if lab.startswith("F_II halved"):
        col = Column(U_INT, D_LONG, ALPHA, KK, FI, FII / 2, L_TOT, 700)
        ts_, rec_, zs_ = col.run(T_END, 2.0, T_SIG, [Z1, Z2])
        Kc_ = [cumulants(ts_, rec_[:, j]) for j in range(2)]
    else:
        kw2 = dict(n_z=700, dt=2.0); kw2.update(kw)
        _, _, zs_, Kc_ = solve(**kw2)
    add(lab, float(plate_height(Kc_, zs_)))

# nu is load-bearing and is stated in a comment in Column.__init__; on a slab
# grid that starts at z = 0 the cylindrical divergence is singular there, so the
# wrong nu does not perturb the answer - it destroys the solve outright.
with np.errstate(invalid="ignore", divide="ignore"):
    _, _, zsN, KcN = solve(n_z=700, dt=2.0, nu=1)
    H_nu = float(plate_height(KcN, zsN))
add("construct_div with nu=1 (cylindrical)", H_nu, "singular at z=0: NaN, not a perturbation")

# a single station instead of two: the entrance effect is no longer cancelled
_, rec1, zs1, Kc1 = solve(n_z=700, dt=2.0)
mu1_in, mu2_in = 5 * T_SIG, T_SIG ** 2                  # exact for the Gaussian slug
H_single = float(zs1[1] * (Kc1[1][2] - mu2_in) / (Kc1[1][1] - mu1_in) ** 2)
add("ONE station (z2 = 600 mm) instead of two", H_single)
_, _, zsS, KcS = solve(n_z=700, dt=2.0, stations=(20.0, Z2))
H_single_short = float(zsS[0] * (KcS[0][2] - mu2_in) / (KcS[0][1] - mu1_in) ** 2)
add("ONE station at z = 20 mm (a short column)", H_single_short)

# a deliberate non-mover: the inlet condition, which the difference must cancel
_, _, zsD, KcD = solve(n_z=700, dt=2.0,
                       bc_left={"a": 0.0, "b": 1.0, "d": 1.0})   # Dirichlet inlet
H_DIR = float(plate_height(KcD, zsD))
add("Dirichlet inlet instead of Danckwerts", H_DIR, "DESIGNED non-mover")
BC_MOVE = float(brk[-1][2])                          # against the BASELINE run
# what breaks that insensitivity: the same inlet change read from ONE station
# (z2), where nothing cancels it. Same two runs, no new solve.
H_single_D = float(zsD[1] * (KcD[1][2] - mu2_in) / (KcD[1][1] - mu1_in) ** 2)
BC_MOVE_SINGLE = float(abs(dev(H_single_D, H_single)))

print(f"{'break':<44}{'H (mm)':>10}{'move vs base':>14}{'|dev| eq (38)':>15}  note")
for lab, v, dbase, d38, note in brk:
    print(f"{lab:<44}{v:>10.5f}{dbase:>14.2e}{d38:>15.2e}  {note}")
moves = [d for lab, v, d, d38, note in brk[1:-1] if np.isfinite(d)]
BREAK_MIN_MOVE = float(min(moves))
print(f"\n   smallest move among the {len(moves)} genuine breaks: {BREAK_MIN_MOVE:.1e}, "
      f"i.e. {BREAK_MIN_MOVE/abs(H_REL):.0f}x the baseline's own error against eq (38)")
print(f"   nu = 1 gives {H_nu} - the wrong geometry does not bias this model, it")
print("   removes it, which is a weaker kind of alarm and is labelled as one.")
print(f"\n   Changing the INLET CONDITION moves H by less than 1e-5 relative to the")
print(f"   baseline (here {BC_MOVE:.1e}) - the two-station difference cancels it, as")
print("   designed and claimed. Read from ONE station (z2) the same change moves H by")
print(f"   {BC_MOVE_SINGLE:.1e}: that is the row that breaks the claim.")

# %% [markdown]
# ## What pymrm adds
#
# **Not the plate height.** Eq. (38) is exact for this model, the Laplace route
# in the implementation section proves it, and the PDE reproduces it to
# $4\times10^{-5}$. Nothing here improves on it.
#
# What the paper leaves open is the *other* half of eqs. (33)/(34). The Gaussian
# is stated to hold "for locations $z$ much larger than both $2D/u$ and
# $F_I u/\alpha$", and the appendix restates the same two conditions as
# $2D/(uz) \ll 1$ and $\beta^2 F_{II}^2 \, 2u/(K^2 F_I \alpha z) \ll 1$. How much
# smaller than 1 is never said, and it cannot be read off Table 3 or Table 4:
# those record widths, not shapes.
#
# It can be computed. The **cumulants** of the response are exactly additive along
# the column, so $\kappa_2$, $\kappa_3$ and $\kappa_4$ *per unit length* fix the
# skewness and excess kurtosis of a column of any length in closed form:
#
# $$\gamma_1(z) = \frac{\kappa_3/z}{(\kappa_2/z)^{3/2}}\,\frac{1}{\sqrt z},
# \qquad
# \gamma_2(z) = \frac{\kappa_4/z}{(\kappa_2/z)^{2}}\,\frac{1}{z}.$$
#
# The Gaussian has $\gamma_1 = \gamma_2 = 0$, so these *are* the error the two
# conditions are supposed to control.
#
# **Where the four cumulants come from matters, so it is stated plainly.**
# $\kappa_1$ to $\kappa_4$ all come from the Laplace expansion of eqs. (29)/(30)
# in the implementation section, they are exact for this model, and every number
# quoted below is one of them. What the pymrm run adds is the *check*: the same
# two-station reading that produced $H$ also produces $\kappa_3$ and $\kappa_4$
# numerically, and the agreement is reported below with the caveat it needs. The
# PDE's $\kappa_3$ and $\kappa_4$ are used for **nothing but that check**, the
# comparison table and the two agreement metrics it feeds, and no shape number
# quoted anywhere on the page. Two printed numbers further down do use PDE
# moments, and are flagged where they appear: the peak-normalised
# $\max|\mathrm{curve}-\mathrm{Gaussian}|$ under the figure below overlays each
# computed curve with a Gaussian built from *that curve's own* mean and variance.
# Those two are a configuration-dependent scale, not agreement metrics, and
# nothing rests on them.

# %%
kap_num = kap_f(D_LONG, U_INT, ALPHA, FI, FII, KK)      # analytic, per unit length
dz_r = zs_r[1] - zs_r[0]
k2_pdes = (Kc_r[1][2] - Kc_r[0][2]) / dz_r
k3_pdes = (Kc_r[1][3] - Kc_r[0][3]) / dz_r
k4_pdes = (Kc_r[1][4] - Kc_r[0][4]) / dz_r

print(f"{'per unit length':<18}{'analytic (Laplace)':>22}{'pymrm PDE':>16}{'dev':>10}")
for nm, an, nu_ in (("kappa_1  (s/mm)", kap_num[0], (Kc_r[1][1]-Kc_r[0][1])/dz_r),
                    ("kappa_2 (s2/mm)", kap_num[1], k2_pdes),
                    ("kappa_3 (s3/mm)", kap_num[2], k3_pdes),
                    ("kappa_4 (s4/mm)", kap_num[3], k4_pdes)):
    print(f"{nm:<18}{float(an):>22.6g}{nu_:>16.6g}{dev(nu_, float(an)):>+10.2%}")
K3_DEV = float(abs(dev(k3_pdes, float(kap_num[2]))))
K4_DEV = float(abs(dev(k4_pdes, float(kap_num[3]))))
k2a, k3a, k4a = float(kap_num[1]), float(kap_num[2]), float(kap_num[3])

print(f"\n   Read those two deviations carefully: {K3_DEV:.1e} and {K4_DEV:.1e} are NOT")
print("   converged discretisation errors. The grid error on the higher cumulants is")
print("   POSITIVE and the time error NEGATIVE, so at fixed grid they cancel - and the")
print("   cancellation passes through zero inside the dt range used on this page:")
_, _, zsK, KcK = solve(n_z=1400, dt=2.0)
KAPZ[(1400, 2.0)] = kappas_per_z(KcK, zsK)
KAPZ[(1400, 1.0)] = (k3_pdes, k4_pdes)
print(f"   {'n_z':>6}{'dt (s)':>9}{'kappa_3 dev':>14}{'kappa_4 dev':>14}")
for dt in (0.5, 1.0, 2.0):
    k3_, k4_ = KAPZ[(1400, dt)]
    print(f"   {1400:>6}{dt:>9.1f}{dev(k3_, k3a):>+14.2e}{dev(k4_, k4a):>+14.2e}")
print(f"   So {K3_DEV:.1e} is what the PDE gives AT THESE SETTINGS, and is an upper bound")
print("   on nothing: it must not be read as 'the third cumulant agrees to 0.2 %' in")
print("   general. The quantity that does converge cleanly here is H itself")
print(f"   ({abs(H_REL):.1e}, grid order {ORD_GRID:.2f}, moving only {DT_SWING_H:.1e} over 16x in dt),")
print("   because the two-station difference cancels the time error and the higher")
print("   cumulants have no such protection. Everything quoted below is analytic.")

# Both quoted shape numbers use the ANALYTIC cumulants; the PDE's are the check.
skew = lambda z: k3a / (k2a ** 1.5 * np.sqrt(z))
exk = lambda z: k4a / (k2a ** 2 * z)
G1 = lambda z: 2 * D_LONG / (U_INT * z)                  # the paper's first condition
G2 = lambda z: FI * U_INT * INV_ALPHA / z                # transfer-unit height / z

print(f"\n{'z (mm)':>9}{'2D/(uz)':>11}{'F_I u/(alpha z)':>18}{'skewness':>11}{'ex. kurtosis':>14}")
for z in (5.0, 20.0, L_SW, 1800.0):
    print(f"{z:>9.0f}{G1(z):>11.4f}{G2(z):>18.4f}{skew(z):>11.4f}{exk(z):>14.4f}")
SKEW_SW = float(skew(L_SW))
SKEW_GLC = float(skew(1800.0))
print(f"\n   At the ion-exclusion column recovered from Table 3 ({L_SW:.0f} mm) the")
print(f"   paper's two groups are {G1(L_SW):.4f} and {G2(L_SW):.4f} - comfortably '<< 1' -")
print(f"   and the elution curve is still skewed by {SKEW_SW:.3f}.")

# skew(z) = (k3a / k2a^1.5) z^(-1/2) falls monotonically, so each level is crossed
# once, at a length in closed form; G1 = 2D/(uz) likewise.
z_at_skew = lambda level, k2=k2a, k3=k3a: (k3 / (k2 ** 1.5 * level)) ** 2
Z_SKEW = float(z_at_skew(0.05))
Z_G1 = 2 * D_LONG / (U_INT * 1e-2)                       # where 2D/(uz) = 1e-2
Z_SKEW_ROOT = one_root(lambda z: skew(z) - 0.05, 1.0, 1e9, "skewness = 0.05", log=True,
                       xtol=1e-9)
print(f"\n   In closed form: skewness falls to 0.05 only at z = {Z_SKEW:.0f} mm")
print(f"   = {Z_SKEW/L_SW:.1f} times the column actually used, where the paper's groups")
print(f"   are {G1(Z_SKEW):.2e} and {G2(Z_SKEW):.2e}.")
print("\n   The reason is in the exponents. Both of the paper's groups fall as 1/z;")
print("   the skewness falls only as 1/sqrt(z) and the excess kurtosis as 1/z, so")
print("   'the groups are small' and 'the curve is Gaussian' are not the same")
print("   statement and differ by a square root. Satisfying 2D/(uz) < 1e-2 leaves")
print(f"   a skewness of {skew(Z_G1):.3f}; reaching a")
print(f"   skewness of 0.05 needs it below {G1(Z_SKEW):.0e}.")
if abs(Z_SKEW_ROOT / Z_SKEW - 1) > 1e-9:
    raise RuntimeError(f"root-found length {Z_SKEW_ROOT} against the closed form {Z_SKEW}")
print(f"   (find_roots on 1 mm to 1 km finds the one crossing, at the closed-form length "
      f"to below 1e-9: here {abs(Z_SKEW_ROOT / Z_SKEW - 1):.0e})")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.2))
zz = np.logspace(0.3, 4.2, 300)
axes[0].loglog(zz, np.abs(skew(zz)), "k", lw=2, label=r"skewness $\gamma_1 \propto z^{-1/2}$")
axes[0].loglog(zz, np.abs(exk(zz)), color="tab:purple", lw=1.7,
               label=r"excess kurtosis $\gamma_2 \propto z^{-1}$")
axes[0].loglog(zz, G1(zz), "--", color="tab:blue", lw=1.4, label=r"$2D/(uz)$")
axes[0].loglog(zz, G2(zz), ":", color="tab:red", lw=1.6, label=r"$F_I u/(\alpha z)$")
axes[0].axvline(L_SW, color="0.4", lw=1.0)
axes[0].annotate(f"ion-exclusion column\n(recovered, {L_SW:.0f} mm)", (L_SW, 2e-3),
                 fontsize=7.5, ha="right", color="0.3")
axes[0].axhline(0.05, color="tab:green", lw=1.0, ls="-.")
axes[0].set(xlabel="column length $z$ (mm)", ylabel="magnitude",
            ylim=(1e-5, 20), title="The paper's two conditions fall faster than the error they bound")
axes[0].legend(fontsize=8, loc="lower left")

for j, z in enumerate((Z1, Z2)):
    y = rec_r[:, j] / np.trapezoid(rec_r[:, j], ts_r)
    axes[1].plot(ts_r, y / y.max(), lw=1.7, label=f"pymrm, $z = {zs_r[j]:.0f}$ mm")
    g = np.exp(-0.5 * ((ts_r - Kc_r[j][1]) / np.sqrt(Kc_r[j][2])) ** 2)
    axes[1].plot(ts_r, g, "--", lw=1.2, color="0.4",
                 label="Gaussian, eq. (33)" if j == 0 else None)
axes[1].set(xlabel="time (s)", ylabel="$c_I$, normalised to its own peak",
            title="Elution curves and the Gaussian they are claimed to be")
axes[1].legend(fontsize=8)
plt.tight_layout(); plt.show()

gauss_err = []
for j in range(2):
    y = rec_r[:, j]; y = y / y.max()
    g = np.exp(-0.5 * ((ts_r - Kc_r[j][1]) / np.sqrt(Kc_r[j][2])) ** 2)
    gauss_err.append(float(np.abs(y - g).max()))
GAUSS_ERR_Z2 = gauss_err[1]
print(f"   max |curve - Gaussian| (peak-normalised): {gauss_err[0]:.4f} at z = {zs_r[0]:.0f} mm,"
      f" {gauss_err[1]:.4f} at z = {zs_r[1]:.0f} mm")
print(f"   Each Gaussian uses that curve's OWN pymrm mean and variance - these two")
print("   numbers are the one place on the page where a PDE moment reaches print.")
print(f"   Both curves carry the Danckwerts inlet and the {T_SIG:.0f} s Gaussian slug, so")
print("   they are configuration-dependent in a way the cumulants are not; they are")
print("   shown for scale, and the skewness above is the transferable one.")

# %% [markdown]
# ### Break table, part 2: the algebra and the extractions
#
# The table above moves the PDE. This one moves everything else. Each row names
# the metric it guards, and the coverage printed with the claims is generated from
# the **measured** movers rather than asserted.

# %%
print("Breaks on the algebra and on the extractions - no PDE involved.")
print("Each row is a deliberate error; the last column is the size of the move.\n")
abrk = []

# --- the symbolic chain: 75 -> 80 in eq (49)'s viscous term ------------------
C_bad = ((FI_ ** 2 / (80 * (1 - FI_) ** 2)
          + FI_ * K_ * DI_ / (2 * sp.pi ** 2 * (1 - FI_) * DII_)) / (1 + K_ * FI_ / FII_) ** 2)
H49_bad = 2 * g_ * DI_ / (u_ * dp_) + 2 * l_ + C_bad * u_ * dp_ / DI_
r2_bad = sp.simplify(sp.simplify(H_38.subs({D_: D_48, a_: 1 / inv_alpha}) / dp_) - H49_bad)
CHAIN_BREAK_OK = int(r2_bad != 0)
abrk.append(("sym_chain_residual", "eq (49): 75 -> 80 in the viscous term",
             "0" if r2_bad == 0 else "non-zero", "0 -> non-zero"))

# --- D_II and lambda: the two inputs the regression actually depends on -------
dp_mesh = np.sqrt(0.297 * 0.149)                # US 50 and 100 mesh openings, mm
(twolam_mesh, S_mesh), _ = fit_eq49(u_sup * dp_mesh, Ymid)
DII_MESH = DII_from_slope(S_mesh, D_I_MM) * 1e-6
DII_MESH_DEV = float(abs(dev(DII_MESH, DII_BEST)))
abrk.append(("DII_from_table3_rel_dev", f"d_p = mesh geo-mean {dp_mesh:.4f} mm, not Table 3's",
             f"{DII_MESH:.3e} m2/s", f"{DII_MESH_DEV:.1%}"))
(twolam_drop, S_drop), _ = fit_eq49(X[:-1], Ymid[:-1])
DII_DROP = DII_from_slope(S_drop, D_I_MM) * 1e-6
DII_DROP_DEV = float(abs(dev(DII_DROP, DII_BEST)))
abrk.append(("DII_from_table3_rel_dev", f"drop the fastest row (u = {u_sup[-1]:.3f} mm/s)",
             f"{DII_DROP:.3e} m2/s", f"{DII_DROP_DEV:.1%}"))
def _spread(visc):
    """The D_I sensitivity of section 6 with the 75 of eq (49) replaced by visc."""
    vals = [DII_from_slope(S, d, visc) for d in (1.0e-3, 2.0e-3)]
    vals.append(DII_from_slope(S, None))
    return max(vals) / min(vals) - 1
abrk.append(("DII_D_I_sensitivity", "75 -> 7.5 in eq (49) (viscous term x10)",
             f"spread {_spread(7.5):.1%}", f"{abs(dev(_spread(7.5), DII_SPREAD)):.0%}"))
abrk.append(("lambda_from_table3_rel_dev",
             "restore the 1/X molecular term the paper drops",
             f"lambda = {b3t[1]/2:.4f}", f"{abs(dev(b3t[1]/2, two_lam/2)):.0%}"))
(twolam_d, _), _ = fit_eq49(X[1:], Ymid[1:])
abrk.append(("lambda_from_table3_rel_dev", "drop the slowest row",
             f"lambda = {twolam_d/2:.4f}", f"{abs(dev(twolam_d/2, two_lam/2)):.0%}"))

# --- d_f, eq (56), the slope ratio -------------------------------------------
df_bad = np.sqrt(np.pi ** 2 / 8 * P["eq54_C_n_butane"]
                 * (1 + P["KFI_over_FII_iso_butane"]) ** 2 * DII_FILM
                 / P["KFI_over_FII_iso_butane"]) * 1e3
abrk.append(("df_worst_rel_dev", "iso-butane x used on the n-butane C",
             f"{df_bad:.3f} um", f"{abs(dev(df_bad, df_um['Celite, n-butane'])):.1%}"))
r56_bad = coef * (P["sterchamol_dp"] / df_um["Celite, n-butane"]) ** 2
abrk.append(("eq56_celite_rel_dev", "Sterchamol d_p on the Celite d_f",
             f"{r56_bad:.3e}", f"{abs(dev(r56_bad, P['ratio56_celite'])):.0f}x"))
ratio_swapped = shape(P["KFI_over_FII_n_butane"]) / shape(P["KFI_over_FII_iso_butane"])
abrk.append(("slope_ratio_rel_dev", "the two printed x exchanged",
             f"{ratio_swapped:.4f}", f"{abs(dev(ratio_swapped, ratio_theory)):.0%}"))

# --- Table 1 and Table 3 ------------------------------------------------------
abrk.append(("table1_calc_frac_reproduced", "v sqrt(n) = 41.25 (dS0 = 165 exactly)",
             f"{n_match(P['delta_S0_glc']/4)} of {len(A_vals)}",
             f"{(CALC_N_OK-n_match(P['delta_S0_glc']/4))/len(A_vals):.0%} of rows lost"))
printed_bad = printed.copy(); printed_bad[3] = 0.88          # A = 122 ml: 0.86 -> 0.88
W_LO_BAD, W_HI_BAD = calc_window(printed_bad)
abrk.append(("table1_dS0_window_lo_ml", "one printed calc cell altered (A = 122: 0.86 -> 0.88)",
             f"window [{W_LO_BAD:.3f}, {W_HI_BAD:.3f}]",
             "window empties" if W_LO_BAD >= W_HI_BAD else "no move"))
L_bad = np.concatenate([(T3.n_plates_lo * T3.hetp_hi_mm).to_numpy(),
                        (T3.n_plates_hi * T3.hetp_lo_mm).to_numpy()])
L_bad[0] = T3.n_plates_lo.iloc[0] * 9.3          # 0.93 mis-read as 9.3, a plausible slip
abrk.append(("table3_recovered_length_mm", "Table 3 cell 0.93 mis-read as 9.3 mm",
             f"mean {L_bad.mean():.0f} mm ({L_bad.min():.0f}-{L_bad.max():.0f})",
             f"{abs(dev(L_bad.mean(), L_SW)):.0%}, and the band stops being a band"))
n_solvable = 0
for _, r in single.iterrows():                 # rows with a root of the printed eq (21)
    f_ = lambda w: width_red(r.feed_volume_ml / w, "delta") * w - r.width_at_base_ml
    n_solvable += int(len(find_roots(f_, r.feed_volume_ml / 25.0, 1e4, n_scan=400,
                                     log=True)["roots"]) > 0)
EQ21_N_SOLVABLE = int(n_solvable)
abrk.append(("table3_n_median_rel_dev", "eq (21) EXACTLY AS PRINTED for the width route",
             f"{EQ21_N_SOLVABLE} of {len(single)} rows solvable",
             "the width route ceases to exist"))

EQ21_LARGE_A_PRINTED = float(width_red(5.0, "delta"))
abrk.append(("eq21_corrected_large_a_dev", "eq (21) as printed at a = 5, vs eq (25)",
             f"{EQ21_LARGE_A_PRINTED:.4g} vs {5+np.sqrt(2*np.pi):.4f}",
             f"{abs(dev(EQ21_LARGE_A_PRINTED, 5+np.sqrt(2*np.pi))):.3g}x"))
abrk.append(("eq21_printed_a0_limit", "the corrected denominator, a -> 0",
             f"{EQ21_A0_CORRECTED:.4f} vs {EQ21_A0_PRINTED:.4f}",
             f"{abs(dev(EQ21_A0_CORRECTED, EQ21_A0_PRINTED)):.0%}"))

# --- the fit residual, and the Sterchamol constant ---------------------------
u_nb, h_nb = SERIES["Celite, n-butane"]
u_ib, h_ib = SERIES["Celite, iso-butane"]
rmse_bad = float(np.sqrt(((vd(u_ib, P["eq54_B_iso_butane"], P["eq54_A_iso_butane"],
                              1.5 * P["eq54_C_iso_butane"]) - h_ib) ** 2).mean()))
abrk.append(("eq54_fit_rmse_mm", "C for iso-butane (the max-rmse series) x1.5",
             f"rmse {rmse_bad:.4f} mm",
             f"{abs(dev(rmse_bad, fitstat['Celite, iso-butane']['rmse'])):.0%}"))
u_keep = u_nb > 2.5
null_bad = float(np.sqrt(((h_nb[u_keep] - h_nb[u_keep].mean()) ** 2).mean()))
abrk.append(("eq54_null_rmse_mm", "the two slowest Celite rows dropped",
             f"null rmse {null_bad:.4f} mm",
             f"{abs(dev(null_bad, fitstat['Celite, n-butane']['null'])):.0%}"))
b3_off, se3_off, _ = refit(us_, hs_ + 0.3)
abrk.append(("sterchamol_free_A_t_stat", "+0.3 mm added to every Sterchamol H",
             f"|t| = {abs(b3_off['A']/se3_off['A']):.2f}",
             f"from {abs(A_STER_T):.2f}"))

# --- the cumulant results -----------------------------------------------------
kap_half = kap_f(D_LONG, U_INT, ALPHA / 2, FI, FII, KK)
skew_half = float(kap_half[2]) / (float(kap_half[1]) ** 1.5 * np.sqrt(L_SW))
abrk.append(("skewness_at_recovered_column", "alpha halved (analytic kappa_3)",
             f"skew {skew_half:.4f}", f"{abs(dev(skew_half, SKEW_SW)):.0%}"))
Z_half = float(z_at_skew(0.05, float(kap_half[1]), float(kap_half[2])))
abrk.append(("z_for_skewness_0p05_mm", "same, on the root-found length",
             f"z = {Z_half:.0f} mm", f"{abs(dev(Z_half, Z_SKEW)):.0%}"))
_, _, zsT, KcT = solve(n_z=700, dt=2.0, t_end=0.6 * T_END)
k3_trunc, k4_trunc = kappas_per_z(KcT, zsT)
abrk.append(("pymrm_kappa3_vs_laplace_rel", "time window truncated to 60 %",
             f"kappa_3/z = {k3_trunc:.4g}", f"{abs(dev(k3_trunc, k3a)):.0%}"))
abrk.append(("pymrm_kappa4_vs_laplace_rel", "time window truncated to 60 %",
             f"kappa_4/z = {k4_trunc:.4g}", f"{abs(dev(k4_trunc, k4a)):.0%}"))

print(f"{'metric':<32}{'break':<50}{'value':>30}   move")
for m, lab, val, mv in abrk:
    print(f"{m:<32}{lab:<50}{val:>30}   {mv}")
print(f"\n   {len(abrk)} rows, covering "
      f"{len(set(m for m, *_ in abrk))} of the metrics reported below.")

# %% [markdown]
# ### Claims
#
# Every number above is registered once, with how it is checked
# (`docs/claims-design.md`): a **value** against the committed baseline within
# 5 %, a **bound** against its stated limits on every run and every machine, an
# **exact** result by equality. Residuals, identities, tolerance gaps and the
# fitted convergence orders are bounds, whatever their size.
#
# Some bounds need a reason. `pymrm_bc_insensitivity`, `pymrm_dt_swing_H` and
# `pymrm_deferred_vs_newton_H_rel` are differences of nearly equal plate heights,
# claimed only as below a limit. `pymrm_grid_order` is an observed order from
# successive differences, claimed as a range, and `pymrm_dt_order_upwind_dev` is
# the departure of the linear Crank-Nicolson march from second order; each order
# comes with a round-off margin.
#
# Two values need a different warning. `pymrm_kappa3_vs_laplace_rel` and
# `pymrm_kappa4_vs_laplace_rel` are **cancellation-limited**, not converged: the
# grid and time-step errors on the higher cumulants have opposite signs and the
# reported values sit near where they cancel; the table under *What pymrm adds*
# shows how they move with dt. They are honest at the stated settings and are regression
# guards like any other, but neither is evidence that the PDE reproduces a
# cumulant to that precision. The converging companion is `pymrm_H_vs_eq38_rel`,
# which is grid-limited and barely moves with $\Delta t$.
#
# Every claim has a break row that fails its own check: break table 1, break
# table 2, the bare-upwind and backward-Euler runs of the refinement study, and,
# for the identities, a deliberately wrong algebra (75 to 80 in eq. (49)) and
# a deliberately inverted closed form $u_{\text{opt}} = \sqrt{C/B}$. The rows are
# declared before they run, and the table below them is printed with every row.

# %%
# --- the second, independent computation of the headline ---------------------
claims.value("pymrm_H_vs_eq38_rel", abs(H_REL),
             note="discretisation error at n_z = 1400, dt = 1 s; grid-limited")
claims.bound("pymrm_grid_order", ORD_GRID, min=2.5, max=3.5, basis="chosen",
             note="observed order of the two-station H from successive differences on"
                  " n_z = 350, 700, 1400, 2800 (dt = 0.5 s), the one farther from 3; van"
                  " Leer is formally second order and H converges faster, so only the range"
                  " is claimed")
claims.bound("pymrm_grid_roundoff_margin", GRID_ROUNDOFF, min=1e3, basis="theory",
             note="smallest successive difference of the grid ladder over n_steps eps H")
claims.value("pymrm_upwind_dispersion_coeff", UPWIND_COEF,
             note="mean (H - H_38)/h with bare upwind; the prediction is exactly 1")
claims.bound("pymrm_upwind_dispersion_coeff_err", UPWIND_COEF_DEV, max=1e-4, basis="chosen",
             note="|mean (H - H_38)/h - 1| with bare upwind over n_z = 350, 700, 1400 (dt = 0.5 s):"
                  " the leading-order coefficient is exactly 1, so this is the O(h) remainder plus"
                  " the time error; the bound admits that remainder and fails any scheme whose"
                  " numerical dispersion is not u h / 2")
claims.bound("pymrm_dt_order_upwind_dev", DT_ORDER_DEV, max=0.1, basis="theory",
             note="worst |order - 2| of c_I(z2) at t = 4800 and 5200 s over dt = 8 to 0.5 s,"
                  " linear Crank-Nicolson (limiter off), successive differences")
claims.bound("pymrm_dt_roundoff_margin", DT_ROUNDOFF, min=1e3, basis="theory",
             note="smallest finest successive difference of c_I(z2) over n_steps eps c_I")
claims.bound("pymrm_dt_swing_H", DT_SWING_H, max=1e-6, basis="chosen",
             note="relative spread of the two-station H, van Leer, over dt = 8 to 0.5 s;"
                  " the bound keeps it 40x below the baseline's own error against eq (38)")
claims.bound("pymrm_deferred_vs_newton_H_rel", DEFER_VS_NEWTON, max=1e-5, basis="chosen",
             note="two-station H of the one-sweep deferred correction against the same theta"
                  " scheme solved by newton with the ConvFluxTVD Jacobian, n_z = 350, dt = 2 s;"
                  " the bound keeps it below the reported run's own error against eq (38)")
claims.exact("pymrm_dt_mu1_sign_changes", MU1_SIGN_CHANGES,
             note="sign changes in the successive differences of mu1(z2) on the van Leer dt ladder"
                  " (n_z = 700): the limited march has no regular time order there")
claims.bound("pymrm_mu1_deferred_vs_newton_rel", MU1_DEFER_VS_NEWTON, max=1e-4, basis="chosen",
             note="mu1(z2) of the one-sweep deferred correction against the same theta scheme"
                  " solved by newton, n_z = 350, dt = 2 s; the lag moves the absolute moment")
claims.bound("pymrm_mu1z1_deferred_vs_newton_rel", MU1Z1_DEFER_VS_NEWTON, max=1e-4, basis="chosen",
             note="mu1(z1) of the one-sweep deferred correction against the same theta scheme"
                  " solved by newton, n_z = 350, dt = 2 s; printed beside the mu1(z2) and H gaps")
claims.value("pymrm_kappa3_vs_laplace_rel", K3_DEV,
             note="cancellation-limited at n_z = 1400, dt = 1 s; a regression guard")
claims.value("pymrm_kappa4_vs_laplace_rel", K4_DEV,
             note="cancellation-limited at n_z = 1400, dt = 1 s; a regression guard")
claims.bound("pymrm_bc_insensitivity", BC_MOVE, max=1e-5, basis="chosen",
             note="relative move of the two-station H when the Danckwerts inlet is"
                  " replaced by a Dirichlet one; a difference of nearly equal plate"
                  " heights. The bound keeps it 4x below the baseline's own error against eq (38)")
# --- the paper's stated results, recomputed ----------------------------------
claims.value("DII_from_table3_m2s", DII_BEST, unit="m2/s")
claims.value("DII_from_table3_rel_dev", DII_DEV)
claims.value("DII_D_I_sensitivity", DII_SPREAD)
claims.value("lambda_from_table3", two_lam / 2)
claims.value("lambda_from_table3_rel_dev", LAMBDA_DEV)
claims.value("df_celite_n_butane_um", df_um["Celite, n-butane"], unit="um")
claims.value("df_celite_iso_butane_um", df_um["Celite, iso-butane"], unit="um")
claims.value("df_sterchamol_um", df_um["Sterchamol, n-butane"], unit="um")
claims.value("df_worst_rel_dev", DF_MAX_DEV)
claims.value("eq56_coefficient", coef, note="D_II / (100 D_I) from the printed constants")
claims.value("eq56_celite_rel_dev", EQ56_DEV)
claims.value("slope_ratio_theory", ratio_theory)
claims.value("slope_ratio_rel_dev", SLOPE_RATIO_DEV)
claims.value("gamma_iso_butane_superficial", GAMMA_ISO)
claims.value("table3_n_median_rel_dev", N_MED)
claims.value("table3_n_max_rel_dev", N_MAX)
claims.value("table3_n_median_height_route", N_MED_H)
claims.value("table3_n_median_width_route", N_MED_W)
claims.exact("table1_calc_frac_reproduced", CALC_FRAC,
             note="fraction of Table 1's eleven printed calc values reproduced")
claims.value("table1_dS0_window_lo_ml", 4 * W_LO, unit="ml")
claims.value("table1_dS0_window_hi_ml", 4 * W_HI, unit="ml")
claims.value("table3_recovered_length_mm", L_SW, unit="mm")
# --- fit, not test: labelled as such everywhere it appears -------------------
claims.value("eq54_fit_rmse_mm", max(f["rmse"] for f in fitstat.values()), unit="mm")
claims.value("eq54_null_rmse_mm", min(f["null"] for f in fitstat.values()), unit="mm")
claims.value("sterchamol_free_A_t_stat", abs(A_STER_T))
# --- the printed defect ------------------------------------------------------
claims.value("eq21_printed_a0_limit", EQ21_A0_PRINTED)
claims.value("eq21_corrected_large_a_dev", EQ21_LARGE_A_DEV,
             note="the corrected eq (21) at a = 5 against the a + sqrt(2 pi) of eq (25);"
                  " an asymptotic error, not round-off")
claims.exact("eq21_printed_solvable_rows", EQ21_N_SOLVABLE,
             note="Table 2 rows (A = 5 ml) for which the width route with eq (21) as printed"
                  " has a root on [A/25, 1e4] ml (find_roots, 400-point log scan)")
# --- what pymrm adds ---------------------------------------------------------
claims.value("group1_at_recovered_column", G1(L_SW), note="2D/(uz) at the recovered length")
claims.value("group2_at_recovered_column", G2(L_SW), note="F_I u/(alpha z) at the recovered length")
claims.value("skewness_at_recovered_column", SKEW_SW)
claims.value("z_for_skewness_0p05_mm", Z_SKEW, unit="mm")
claims.value("z_for_skewness_over_column", Z_SKEW / L_SW)
# --- identities --------------------------------------------------------------
claims.bound("sym_chain_residual", 0.0 if r2 == 0 else 1.0, max=0.0, basis="identity",
             note="sympy simplifies (38) with (45), (46), (48) substituted, minus (49), to 0;"
                  " recorded as 1 if it does not")
claims.exact("sym_chain_exact_count", CHAIN_OK,
             note="steps of the symbolic chain whose residual sympy simplifies to 0, of 4")
claims.bound("rootfind_vs_closed_form_optimum", ROOT_ID, max=1e-12, basis="identity",
             note="max relative gap between u_opt root-found on dH/du = -B/u^2 + C = 0"
                  " (find_roots, 1e-3 to 1e3 mm/s, xtol 1e-15) and the closed form sqrt(B/C):"
                  " an identity that tests the root finder on that expression only")

# --- break rows: the value each claim takes with the defect injected ----------
# Break table 1's defects are solved again at the claim's own setting
# (n_z = 1400, dt = 1 s), each row its own solve with the defect injected.
def h_defect(lab, n_z=1400, dt=1.0):
    """H with break table 1's defect `lab` at (n_z, dt): two-station, or one-station
    for the two single-station rows."""
    if lab.startswith("F_II halved"):
        col = Column(U_INT, D_LONG, ALPHA, KK, FI, FII / 2, L_TOT, n_z)
        ts_, rec_, zs_ = col.run(T_END, dt, T_SIG, [Z1, Z2])
        return float(plate_height([cumulants(ts_, rec_[:, j]) for j in range(2)], zs_))
    if lab.startswith("ONE station (z2"):
        _, _, z_, K_ = solve(n_z=n_z, dt=dt)
        return float(z_[1] * (K_[1][2] - mu2_in) / (K_[1][1] - mu1_in) ** 2)
    if lab.startswith("ONE station at z = 20"):
        _, _, z_, K_ = solve(n_z=n_z, dt=dt, stations=(20.0, Z2))
        return float(z_[0] * (K_[0][2] - mu2_in) / (K_[0][1] - mu1_in) ** 2)
    kw = dict(n_z=n_z, dt=dt)
    kw.update({**BRK1_KW, "construct_div with nu=1 (cylindrical)": dict(nu=1)}[lab])
    with np.errstate(invalid="ignore", divide="ignore"):
        _, _, z_, K_ = solve(**kw)
    return float(plate_height(K_, z_))


def df_of(C, x, dii=DII_FILM):
    """Film thickness (um) from eq (53)'s C, x = K F_I/F_II and the film D_II."""
    return float(np.sqrt(np.pi ** 2 / 8 * C * (1 + x) ** 2 * dii / x) * 1e3)


_STALE = {}
def _stale_cumulants():
    """Cumulants of the deferred route with stale inlet data, n_z = 350, dt = 2 s (one solve)."""
    if not _STALE:
        _STALE["K"] = solve(n_z=N_NEWTON, dt=DT_NEWTON, inlet_data=False)[3]
    return _STALE["K"]


def _defer_upwind_vs_newton():
    _, _, z_, K_ = solve(n_z=N_NEWTON, dt=DT_NEWTON, limiter=None)
    return float(abs(dev(plate_height(K_, z_), H_NEWTON)))


def _solvable(denominator):
    n = 0
    for _, r in single.iterrows():
        f_ = lambda w: width_red(r.feed_volume_ml / w, denominator) * w - r.width_at_base_ml
        n += int(len(find_roots(f_, r.feed_volume_ml / 25.0, 1e4, n_scan=400,
                                log=True)["roots"]) > 0)
    return n


def _chain_count_bad():
    return int(sum(r == 0 for r in (r1, r2_bad, r3, r4)))


def _dt_repeated():
    with np.errstate(divide="ignore", invalid="ignore"):
        return dt_study({**LAD_UP, STEPS[-1]: LAD_UP[STEPS[-2]]})[2]


recs_nv = table3_from_table2(1.1 * NV)
devs_nv = np.abs([d for *_, dlo, dhi, _, _ in recs_nv for d in (dlo, dhi)])
mh_nv, mw_nv = route_medians(recs_nv)
D_LAM8 = GAMMA_ASSUMED * D_I_ASSUMED + 8.0 * U_INT * DP
Y_LO = FLOW.h_over_dp_lo.to_numpy(float)
(two_lam_lo, S_lo), _ = fit_eq49(X, Y_LO)
BREAKS = [(f"break table 1: {lab}", "pymrm_H_vs_eq38_rel", lambda lab=lab: abs(dev(h_defect(lab), H38)))
          for lab, *_ in brk[1:-1]]
BREAKS += [
    ("ONE station: Dirichlet inlet against Danckwerts", "pymrm_bc_insensitivity", lambda: BC_MOVE_SINGLE),
    ("bare upwind in place of van Leer", "pymrm_grid_order", lambda: ORD_UPWIND),
    ("grid ladder with its finest level repeated", "pymrm_grid_roundoff_margin",
     lambda: float(min(abs(GRID[min(a, 1400)] - GRID[min(b, 1400)]) for a, b in zip(NZS[:-1], NZS[1:]))
                   / (int(round(T_END / 0.5)) * EPS * H38))),
    ("van Leer runs read as if bare upwind", "pymrm_upwind_dispersion_coeff",
     lambda: float(np.mean([(GRID[n] - H38) / (L_TOT / n) for n in NZS[:3]]))),
    ("van Leer runs read as if bare upwind", "pymrm_upwind_dispersion_coeff_err",
     lambda: abs(float(np.mean([(GRID[n] - H38) / (L_TOT / n) for n in NZS[:3]])) - 1.0)),
    ("backward Euler (theta = 1) in place of Crank-Nicolson", "pymrm_dt_order_upwind_dev",
     lambda: DT_ORDER_DEV_BE),
    ("backward Euler (theta = 1) in place of Crank-Nicolson", "pymrm_dt_swing_H",
     lambda: DT_SWING_H_BE),
    ("stale inlet data in the deferred correction", "pymrm_dt_mu1_sign_changes",
     lambda: mu1_sign_changes(dt_ladder(inlet_data=False))),
    ("stale inlet data in the deferred correction", "pymrm_mu1_deferred_vs_newton_rel",
     lambda: float(abs(dev(_stale_cumulants()[1][1], Kc_n[1][1])))),
    ("stale inlet data in the deferred correction", "pymrm_mu1z1_deferred_vs_newton_rel",
     lambda: float(abs(dev(_stale_cumulants()[0][1], Kc_n[0][1])))),
    ("dt ladder with its finest level repeated", "pymrm_dt_roundoff_margin", _dt_repeated),
    ("deferred route without the limiter, against newton", "pymrm_deferred_vs_newton_H_rel",
     _defer_upwind_vs_newton),
    ("time window truncated to 60 %", "pymrm_kappa3_vs_laplace_rel", lambda: abs(dev(k3_trunc, k3a))),
    ("time window truncated to 60 %", "pymrm_kappa4_vs_laplace_rel", lambda: abs(dev(k4_trunc, k4a))),
    ("d_p = mesh geometric mean", "DII_from_table3_rel_dev",
     lambda: abs(dev(DII_MESH, P["D_II_ion_exclusion"]))),
    ("drop the fastest row", "DII_from_table3_rel_dev", lambda: abs(dev(DII_DROP, P["D_II_ion_exclusion"]))),
    ("75 -> 7.5 in eq (49)", "DII_D_I_sensitivity", lambda: _spread(7.5)),
    ("restore the 1/X molecular term", "lambda_from_table3_rel_dev",
     lambda: abs(dev(b3t[1] / 2, P["lambda_50_100"]))),
    ("drop the slowest row", "lambda_from_table3_rel_dev",
     lambda: abs(dev(twolam_d / 2, P["lambda_50_100"]))),
    ("H/d_p at the low end of each printed range", "lambda_from_table3", lambda: two_lam_lo / 2),
    ("H/d_p at the low end of each printed range", "DII_from_table3_m2s",
     lambda: DII_from_slope(S_lo, D_I_MM) * 1e-6),
    ("iso-butane x used on the n-butane C", "df_worst_rel_dev",
     lambda: max(abs(dev(dict(df_um, **{"Celite, n-butane": df_bad})[n], p)) for n, _, _, p in CASES)),
    ("iso-butane x used on the n-butane C", "df_celite_n_butane_um", lambda: df_bad),
    ("film D_II doubled", "df_celite_iso_butane_um",
     lambda: df_of(P["eq54_C_iso_butane"], P["KFI_over_FII_iso_butane"], 2 * DII_FILM)),
    ("film D_II doubled", "df_sterchamol_um",
     lambda: df_of(P["eq55_C_sterchamol"], P["KFI_over_FII_sterchamol"], 2 * DII_FILM)),
    ("D_I doubled", "eq56_coefficient",
     lambda: P["D_II_liquid_film"] / (100 * 2 * P["D_I_gilliland"])),
    ("Sterchamol d_p on the Celite d_f", "eq56_celite_rel_dev", lambda: abs(dev(r56_bad, P["ratio56_celite"]))),
    ("the two printed x exchanged", "slope_ratio_rel_dev",
     lambda: abs(dev(ratio_swapped, P["slope_ratio_predicted"]))),
    ("the two printed x exchanged", "slope_ratio_theory", lambda: ratio_swapped),
    ("iso-butane B read as interstitial at F_I = 0.5", "gamma_iso_butane_superficial",
     lambda: P["eq54_B_iso_butane"] * 1e-6 / (2 * DI * 0.5)),
    ("eq (21) exactly as printed", "table3_n_median_rel_dev",
     lambda: N_MED if EQ21_N_SOLVABLE == len(single) else np.nan),
    ("eq (21) exactly as printed", "table3_n_max_rel_dev",
     lambda: N_MAX if EQ21_N_SOLVABLE == len(single) else np.nan),
    ("nv = 77 ml in place of 70 (10 % high)", "table3_n_median_height_route", lambda: mh_nv),
    ("nv = 77 ml in place of 70 (10 % high)", "table3_n_median_width_route", lambda: mw_nv),
    ("nv = 77 ml in place of 70 (10 % high)", "table3_n_median_rel_dev", lambda: float(np.median(devs_nv))),
    ("the corrected denominator in the width-route scan", "eq21_printed_solvable_rows",
     lambda: _solvable("a")),
    ("v sqrt(n) = 41.25 (dS0 = 165 exactly)", "table1_calc_frac_reproduced",
     lambda: n_match(P["delta_S0_glc"] / 4) / len(A_vals)),
    ("one printed calc cell altered", "table1_dS0_window_lo_ml",
     lambda: 4 * W_LO_BAD if W_LO_BAD < W_HI_BAD else np.nan),
    ("one printed calc cell altered", "table1_dS0_window_hi_ml",
     lambda: 4 * W_HI_BAD if W_LO_BAD < W_HI_BAD else np.nan),
    ("Table 3 cell 0.93 mis-read as 9.3", "table3_recovered_length_mm", lambda: float(L_bad.mean())),
    ("C for iso-butane x1.5", "eq54_fit_rmse_mm",
     lambda: max([rmse_bad] + [f["rmse"] for n, f in fitstat.items() if n != "Celite, iso-butane"])),
    ("the two slowest Celite rows dropped", "eq54_null_rmse_mm",
     lambda: min([null_bad] + [f["null"] for n, f in fitstat.items() if n != "Celite, n-butane"])),
    ("+0.3 mm on every Sterchamol H", "sterchamol_free_A_t_stat", lambda: abs(b3_off["A"] / se3_off["A"])),
    ("the corrected denominator, a -> 0", "eq21_printed_a0_limit", lambda: EQ21_A0_CORRECTED),
    ("eq (21) as printed at a = 5", "eq21_corrected_large_a_dev",
     lambda: abs(dev(EQ21_LARGE_A_PRINTED, 5 + np.sqrt(2 * np.pi)))),
    ("lambda = 8 (the 200-400 mesh value) on this column", "group1_at_recovered_column",
     lambda: 2 * D_LAM8 / (U_INT * L_SW)),
    ("alpha halved", "group2_at_recovered_column", lambda: FI * U_INT * 2 * INV_ALPHA / L_SW),
    ("alpha halved", "skewness_at_recovered_column", lambda: skew_half),
    ("alpha halved", "z_for_skewness_0p05_mm", lambda: Z_half),
    ("alpha halved", "z_for_skewness_over_column", lambda: Z_half / L_SW),
    ("eq (49): 75 -> 80", "sym_chain_residual", lambda: 0.0 if r2_bad == 0 else 1.0),
    ("eq (49): 75 -> 80", "sym_chain_exact_count", _chain_count_bad),
    ("closed form inverted, u_opt = sqrt(C/B)", "rootfind_vs_closed_form_optimum",
     lambda: max(abs(f["u_root"] / np.sqrt(FITS[n][2] / FITS[n][0]) - 1) for n, f in fitstat.items())),
]
claims.define_rows(list(dict.fromkeys(r for r, _, _ in BREAKS)))
BT = []
for row, name, fn in BREAKS:
    try:
        out = float(fn())
    except Exception as exc:                 # noqa: BLE001 - a failed study is a result
        print(f"  {row}: {name} raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
        out = np.nan
    BT.append((row, name, claims.breaks(row, name, out)))
    base_ = claims.metrics[name]["value"]
    print(f"  {row[:58]:<58} {name:<34} {base_:>11.4g} -> {out:>11.4g}  "
          f"{'caught' if BT[-1][2] else 'NOT caught'}")
_missed = [f"{r} -> {m}" for r, m, c in BT if not c]
_caught = {m for _, m, c in BT if c}
print(f"{len(BT)} break rows over {len(dict.fromkeys(r for r, _, _ in BT))} defects; "
      f"{len(_caught)} of {len(claims.metrics)} claims failed by at least one row")
print("rows that fail no claim: " + ("; ".join(_missed) or "none"))

# --- prose: every number the Markdown and the README quote for a claim ---------
SAYS = {
    "pymrm_bc_insensitivity": ["below 1e-5"],
    "pymrm_grid_order": ["between 2.5 and 3.5", "at least 2.5", "at most 3.5"],
    "pymrm_H_vs_eq38_rel": ["4e-5", "4.1e-5"],
    "skewness_at_recovered_column": ["0.14"],
    "pymrm_upwind_dispersion_coeff": ["0.99999"],
    "pymrm_dt_swing_H": ["below 1e-6"],
    "pymrm_dt_mu1_sign_changes": ["1"],
    "pymrm_dt_order_upwind_dev": ["below 0.1"],
    "pymrm_deferred_vs_newton_H_rel": ["below 1e-5"],
    "lambda_from_table3": ["2.99"],
    "lambda_from_table3_rel_dev": ["0.24 %"],
    "DII_from_table3_m2s": ["1.34e-10"],
    "DII_from_table3_rel_dev": ["3.7 %"],
    "DII_D_I_sensitivity": ["1.8 %"],
    "df_celite_n_butane_um": ["10.1", "10.06"],
    "df_celite_iso_butane_um": ["10.6", "10.57"],
    "df_sterchamol_um": ["8.7", "8.66"],
    "eq56_coefficient": ["3.333e-8"],
    "table1_dS0_window_lo_ml": ["166.0"],
    "table1_dS0_window_hi_ml": ["167.4"],
    "table3_n_median_rel_dev": ["2.1 %"],
    "table3_n_median_height_route": ["1.6 %"],
    "table3_n_median_width_route": ["2.2 %"],
    "eq21_printed_a0_limit": ["2"],
    "eq21_printed_solvable_rows": ["2"],
    "gamma_iso_butane_superficial": ["0.444"],
    "group1_at_recovered_column": ["2e-3", "2.4e-3"],
    "group2_at_recovered_column": ["3.4e-3"],
    "z_for_skewness_over_column": ["8.2"],
    "pymrm_kappa3_vs_laplace_rel": ["0.21 %"],
    "pymrm_kappa4_vs_laplace_rel": ["1.5 %"],
    "slope_ratio_theory": ["1.448"],
    "table3_recovered_length_mm": ["558"],
    "eq21_corrected_large_a_dev": ["1.1e-6"],
    "eq54_fit_rmse_mm": ["0.28"],
    "eq54_null_rmse_mm": ["0.92"],
}
for k, lits in SAYS.items():
    for lit in lits:
        claims.says(lit, k)

_kinds = {kd: [k for k, c in claims.metrics.items() if c["kind"] == kd]
          for kd in ("value", "bound", "exact")}
print(f"\n{len(claims.metrics)} claims: {len(_kinds['value'])} values, "
      f"{len(_kinds['bound'])} bounds, {len(_kinds['exact'])} exact. The bounds:")
for k in _kinds["bound"]:
    c = claims.metrics[k]
    lo, hi = c.get("min"), c.get("max")
    rng = f"[{lo:g}, {hi:g}]" if lo is not None and hi is not None else (
        f">= {lo:g}" if lo is not None else f"<= {hi:g}")
    print(f"  {k:<34} {rng:<12} [{c['basis']}]  here {c['value']:.3g}")
print(f"{sum(len(v) for v in SAYS.values())} prose literals checked with claims.says")

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload

# %% [markdown]
# ## Reuse
#
# **What this page establishes.**
#
# * The chain (34) → (38) → (49) → (52)/(53) is exact, all four steps, and an
#   entirely separate route, the Laplace-transform cumulants of eqs. (29)/(30),
#   reproduces eq. (38) symbolically as well. The 1956 algebra is sound.
# * Solving eqs. (29)/(30) with pymrm and taking moments reproduces eq. (38) to
#   $4\times10^{-5}$, with the bare-upwind numerical dispersion coming out at
#   0.99999 of the predicted $h$.
# * Table 3 refitted through eq. (49) returns $\lambda = 2.99$ against the printed
#   3, and $D_{II} = 1.34\times10^{-10}$ m² s⁻¹ against the printed
#   $1.3\times10^{-10}$, and the extraction turns out to be almost independent of
#   $D_I$, which the paper never prints for that system.
# * The film thickness comes back as 10.1, 10.6 and 8.7 µm against the printed
#   $\approx$ 10, $\approx$ 10 and 9 µm.
# * Eq. (21) **as printed cannot be right**: with the printed
#   $(a+\delta)/\delta$ it gives 2 in the limit its own eq. (24) fixes at 4, and it
#   diverges where its own eq. (25) fixes it at $a + \sqrt{2\pi}$. With $(a+\delta)/a$
#   both limits are met. Reported, not repaired; the corrected form is used and
#   labelled.
# * The Gaussian's shape error falls as $z^{-1/2}$ while both of the paper's own
#   sufficient conditions fall as $z^{-1}$, so satisfying them "$\ll 1$" is weaker
#   than it looks: at the column that produced Table 3, the groups are $2\times10^{-3}$
#   and the skewness is still 0.14 (both printed above).
#
# **What it cannot conclude.**
#
# * **Nothing about whether eq. (53) is right.** Its three constants were fitted to
#   the twelve points they are compared with. Every residual in section 1 is a
#   goodness of fit and is labelled as one; the closest thing to a test in the
#   paper is the slope ratio of section 4, which is one number.
# * **Nothing about $\gamma$ or $d_f$ as physical quantities.** They are the
#   paper's own inversions, redone, and the velocity-convention ambiguity of
#   section 5 cannot be closed from what is printed: under the interstitial
#   reading $\gamma$ is *larger* by $1/F_I$ and $d_f^2$ *smaller* by $F_I$, both by
#   up to the factor printed in section 5. $\lambda$ is the exception: the eddy
#   term carries no $u$, so $2\lambda d_p$ is the one constant here the ambiguity
#   cannot touch.
# * **Nothing about Simpson and Wheaton's measurements.** Their paper was not
#   consulted; Table 2 is van Deemter et al.'s restatement.
# * **Nothing about the ion-exclusion column's geometry.** Its length is printed
#   nowhere; the value used here is recovered from Table 3 (printed above) and is
#   an inference, not a measurement.
# * **Nothing about eq. (53) in a real gas chromatograph.** These columns run at
#   0.046 to 31 mm s⁻¹ with plate heights of 0.56 to 16 mm. Modern capillary columns
#   are three orders of magnitude away in both.
#
# **How to reuse the solver.** `Column` is eqs. (29)/(30) with nothing case-specific
# in it: give it $u$, $D$, $\alpha$, $K$, $F_I$, $F_{II}$ and a length. Any linear
# two-phase column with a first-order exchange is the same object, an adsorption
# breakthrough ([`J1.5`](../J1.5-ldf-breakthrough/)), an ion-exchange column, a
# simulated-moving-bed section. The two-station cumulant reading is the part worth
# copying: it gives the asymptotic plate height without ever needing to know what
# the inlet did, and it costs one extra recording station.
#
# **The single-station form is in the break table** rather than in the results.
# On this column it is only slightly out (break table 1 prints by how much),
# because a Danckwerts inlet happens to put very little of an entrance effect into
# the first two moments here; that is a property of *these* conditions and of
# *that* inlet condition, not a general one, and the two-station form costs one
# extra recording array.
#
# **The time march.** The step is linear except for the van Leer limiter, so
# `Column.run` takes the limiter as a deferred correction and factorises one
# matrix for the whole run. `Column.run_newton` solves the same theta scheme with
# the limited flux and its exact Jacobian (`ConvFluxTVD`, `newton` with monotone
# damping) at every step; on this page the two give the same two-station $H$
# (the bound is in the refinement study) at a fraction of the cost. For a quantity
# that is not differenced between stations, or a larger Courant number, measure
# the gap again before keeping the deferred route.
#
# **Related pages.** [`J1.5`](../J1.5-ldf-breakthrough/) (the linear driving force
# that replaces $\alpha$), `J1.6`, `J1.8` (the general rate model),
# [`A2.3`](../A2.3-taylor-aris-dispersion/) and
# [`A2.6`](../A2.6-gunn-dispersion-correlations/) (where $D$ comes from),
# [`A2.4`](../A2.4-tanks-in-series/) (plates and dispersion as the same statement).
#
# **Cite the source, not this page:** van Deemter, J. J., Zuiderweg, F. J. &
# Klinkenberg, A., *Longitudinal diffusion and resistance to mass transfer as
# causes of nonideality in chromatography*, Chem. Eng. Sci. **5**(6) 271–289
# (1956), [doi:10.1016/0009-2509(56)80003-1](https://doi.org/10.1016/0009-2509(56)80003-1).
