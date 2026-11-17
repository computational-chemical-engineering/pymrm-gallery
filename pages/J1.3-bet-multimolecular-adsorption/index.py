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
# title: "BET 1938: the theory that yields the surface, against the theory that has to be given it"
# description: "Brunauer, Emmett and Teller open by refuting a rival, the polarization theory of de Boer, Zwicker and Bradley, and only then derive their own isotherm. This page rebuilds both. The refutation is reproduced from its own constants: the inter-layer dipole ratio comes out -0.0102 where Bradley's fit needs 0.989, a factor 96.8, and Bradley's k lies outside the domain in which the polarization recursion has any decaying solution at all, confirmed by solving that recursion instead of quoting its closed form. The BET side is then separated into what was fitted and what was not: the seven-gas surface agreement is a real out-of-sample test that collapses a 46.9 % spread in v_m to 10.6 %, while the celebrated agreement between v_m and point B is shown to be near-guaranteed by the equation's own shape, the BET inflection point cannot exceed 2/sqrt(3) v_m for any c, a ceiling here derived in closed form and root-found independently. One printed statement does not survive: '840 +- 70 cal. for nitrogen on all twelve adsorbents' is the eleven-row number."
# categories: [sec:J, struct:S1, struct:S10, tier:T0, data:tier6, phase:gas-solid]
# date: 2026-08-07
# ---
#
# # BET 1938: the theory that yields the surface, against the theory that has to be given it
#
# **Catalog ID:** `J1.3` · **Structures:** `S1` (pointwise algebra / equilibrium
# ladder), `S10` (constrained algebraic system) · **Tier:** T0
#
# Almost everybody meets the BET equation as a line to draw through five points:
#
# $$\frac{p}{v\,(p_0-p)} \;=\; \frac{1}{v_\mathrm{m}c} \;+\; \frac{c-1}{v_\mathrm{m}c}\,\frac{p}{p_0}.
# \tag{A}$$
#
# Two parameters, fitted over $0.05 < p/p_0 < 0.35$, and the fit is good. That
# tells you nothing. A two-parameter fit over a fifth of an isotherm's range is
# supposed to be good, and the theory Brunauer, Emmett and Teller were arguing
# against could fit the same isotherms too, de Boer and Zwicker's polarization
# theory had been doing exactly that for nine years, and the paper says so on its
# first page.
#
# **So this page is organised the way the paper is, not the way the textbooks
# are.** Section I of the 1938 paper is a *refutation of a competitor*, carried out
# with numbers, before a single line of the BET derivation appears. That
# refutation is the page's spine, because it is the only part of the argument
# where two theories are made to disagree about something measurable.
#
# Three questions, in order.
#
# 1. **Can the rival theory be dismissed on its own constants?** The paper says
#    yes and prints the whole chain. Reproduced here: $\alpha/r^3 = 0.029$ for
#    argon from the refractive index and the solid density, a lattice constant
#    $d = -0.35$, and hence an inter-layer dipole ratio $C = -0.0102$, against
#    the $0.989$ that Bradley's own fitted $k = 0.615$ implies. That is a factor
#    **96.8** in $C$ and **9361** in the binding energy it transmits. And
#    Bradley's $k$ is not merely large: it is **outside the range where the
#    polarization recursion has a decaying solution at all**, which this page
#    shows by *solving the recursion* rather than by quoting its closed form.
# 2. **Which parts of the BET case are fits, and which are tests?** Kept strictly
#    apart. The isotherm agreement is a fit and is labelled one everywhere. The
#    thing that is *not* a fit is that seven different gases, each with its own
#    fitted $v_\mathrm{m}$, must give the *same* surface area once converted
#    through cross-sections computed from bulk densities, a conversion with no
#    adsorption in it. That collapses a **46.9 %** spread in $v_\mathrm{m}$ to
#    **10.6 %** in area on silica gel, and **18.6 %** to **7.9 %** on charcoal.
#    The polarization theory cannot even be asked this question: it has no
#    $v_\mathrm{m}$, and de Boer and Zwicker could not evaluate $K_1$ *because*
#    the surface was unknown.
# 3. **Is the famous $v_\mathrm{m} \approx$ point B agreement evidence?** Much
#    less than it looks. This page derives, in closed form and confirms by
#    root-finding, that the BET isotherm's inflection point can never sit above
#    $2/\sqrt{3} = 1.1547\,v_\mathrm{m}$ **for any $c$ whatever**, attained at
#    $c = 27+15\sqrt3$. An eye-read landmark at or below that inflection therefore
#    cannot exceed $v_\mathrm{m}$ by more than 15.5 %, **a one-sided bound**: it
#    can still sit arbitrarily far below, and on Table III's butane isotherm point
#    B is only 0.4828 of $v_\mathrm{m}$. The authors report the two "seldom differing
#    by as much as 10 %" over the twelve isotherms of Table I; the shape of their
#    own equation had already promised the upper half of that.
#
# Along the way one printed statement is shown to be wrong, in the F2.3 order,
# pin what is not free first. The paper's own $\pm$ notation is decoded exactly
# from Table II (four gases, four hits, centre to the nearest 10 cal, half-width
# to the nearest 5), and applied to Table I it shows that "$E_1-E_\mathrm{L}$ is
# uniformly $840 \pm 70$ cal. […] **for nitrogen on all twelve adsorbents**" is
# the number **eleven** of those twelve give. The twelve give $824.5 \pm 86.5$,
# which rounds to $820 \pm 85$ (or $820 \pm 90$ on the coarser half-width rule),
# and two rows fall outside the band as printed.

# %% [markdown]
# ## Background
#
# ### Two explanations of one curve
#
# By the mid-1930s it was settled that low-temperature van der Waals isotherms on
# most adsorbents are S-shaped: concave to the pressure axis at low $p$, convex at
# high $p$. It was not settled why. The paper's own opening sentence lists the two
# camps, capillary condensation, and the build-up of *multimolecular* layers,
# and then names the only quantitative account of the second that existed:
#
# > The isotherm equation which they, and later Bradley, derived on the basis of
# > this polarization theory is practically the only quantitative expression that
# > has been so far proposed to account for multimolecular adsorption.
#
# **The polarization theory.** De Boer and Zwicker supposed the adsorbent
# polarizes the first adsorbed layer, that layer polarizes the second, and so on.
# With $\mu_i$ the dipole moment in layer $i$:
#
# $$\mu_i = c_1 C^i, \tag{1}$$
#
# $$\phi_i = c_2 C^{2i}, \tag{2}$$
#
# binding energy going as the square of the moment. A Boltzmann equilibrium for
# the top layer, $p_n = c_3 e^{-\phi_n/RT}$ (3), then gives
#
# $$\ln\frac{p_n}{c_3} = -\frac{c_2}{RT}\,C^{2n}, \tag{4}$$
#
# or, equivalently,
#
# $$\ln\frac{p_n}{K_3 p_0} = K_2 K_1^{\,n}. \tag{4a}$$
#
# Substituting $n = v/v_\mathrm{m}$ makes (4a) an isotherm. In its logarithmic
# form with $K_3 = 1$ it is the familiar empirical statement that **a plot of
# $\log\log(p/p_0)$ against the amount adsorbed is straight**, which, the paper
# notes, "has often been found to be linear over a considerable range of $p/p_0$
# values".
#
# Two things make this a genuine rival rather than a straw man. It fits: "DeBoer
# in several papers showed that various experimental adsorption isotherms could be
# fitted by equation (4a)." And it had served to extract physics: Bradley
# applied it to argon on copper and aluminium sulfate, got $k \approx 0.6$, and
# concluded that very strong polarization explained adsorbed films of *more than
# 30 layers* at half the saturation pressure.
#
# **Why that is testable.** $C$ is not a free parameter. It is the ratio of
# induced dipoles in successive layers, and induced dipoles are computable from a
# polarizability and a lattice. Section I of the paper computes it. That is the
# whole refutation, and it is the part of the argument this page can check hardest,
# because every input is printed.
#
# **The BET answer.** Section II then generalises Langmuir's kinetic derivation
# from one layer to many, assuming only that the second and all higher layers
# behave like the *liquid*: same heat ($E_2 = E_3 = \dots = E_\mathrm{L}$), same
# evaporation–condensation ratio ($b_i/a_i = g$). Two constants survive,
# $v_\mathrm{m}$ and $c$, and the second turns out to carry the first-layer heat.
#
# ### What is on this page and what is not
#
# `J1.5` owns adsorption *breakthrough*; this page owns the *equilibrium relation*
# and, above all, the discrimination between it and the theory it displaced.
# **No curve is digitised anywhere and no figure is reproduced.** The paper's
# isotherms exist in it only as drawn curves, and everything used here is either a
# printed table, a printed sentence, or a block of constants typeset inside a
# figure frame. Where a figure's constants are used, the page says so.

# %% [markdown]
# ## The published model
#
# Transcribed from renders of the original at the scan's native 300 ppi, with
# every constant cropped and re-read at digit scale. Equation numbers are the
# paper's own; the four lettered equations (A), (B), (C), (E) are lettered in the
# paper too.
#
# ### Section I, the polarization theory, and its own arithmetic
#
# The polarization recursion, from footnote 3 on journal page 309:
#
# $$\mu_i = k\,(\mu_{i-1} + \mu_{i+1}), \tag{1a}$$
#
# which (1) solves if
#
# $$C = \frac{1-\sqrt{1-4k^2}}{2k}, \tag{1b}$$
#
# and de Boer and Zwicker's approximation
#
# $$C = \frac{k}{1-k^2}. \tag{1c}$$
#
# **Equation (1b) is real only for $k \le 1/2$**, and the footnote says so
# explicitly, immediately after recording that Bradley derived $k = 0.6075$ and
# $k = 0.615$ using (1c).
#
# The magnitude of $C$, journal page 310. Polarizability from the refractive
# index, nearest-neighbour distance from the solid molar volume in a
# face-centred-cubic lattice:
#
# $$2\pi\alpha = (n-1)\,v_\mathrm{g}, \tag{5}$$
#
# $$r^3 = \sqrt{2}\,v_\mathrm{s}, \tag{6}$$
#
# $$\frac{\alpha}{r^3} = \frac{n-1}{2^{3/2}\pi}\,\frac{v_\mathrm{g}}{v_\mathrm{s}}. \tag{7}$$
#
# For argon at 0 °C and 760 mm, $n-1 = 278\times10^{-6}$, and with the density of
# solid argon at 40 K,
#
# $$\alpha/r^3 = 0.029, \tag{8}$$
#
# $$\frac{\mu_{i+1}}{\mu_i} = C = d\,(\alpha/r^3), \tag{9}$$
#
# with $d = -0.35$ for close-packed spheres with the dipoles normal to the
# surface. Footnote 10 gives the lattice sums behind $d$ and the reduced recursion
#
# $$\mu_i = \frac{\alpha}{r^3}\Bigl\{-11.1\,\mu_i - 0.466\,(\mu_{i-1}+\mu_{i+1})\Bigr\}
# \;\;\longrightarrow\;\;
# \mu_i = -0.35\,\frac{\alpha}{r^3}\,(\mu_{i-1}+\mu_{i+1}),$$
#
# so that $k = 0.35\,\alpha/r^3$ and, since $k \ll 1$, $k \approx C$. The result on
# journal page 311: $C = -0.01$, and $K_1 = C^2 \approx 1\times10^{-4}$.
#
# ### Section II, the BET derivation
#
# $s_i$ is the surface area covered by exactly $i$ layers. Equilibrium of the bare
# surface with the first layer, and then of every layer with the next:
#
# $$a_1 p\,s_0 = b_1 s_1 e^{-E_1/RT}, \tag{10}$$
#
# $$a_2 p\,s_1 = b_2 s_2 e^{-E_2/RT}, \tag{11}$$
#
# $$a_i p\,s_{i-1} = b_i s_i e^{-E_i/RT}. \tag{12}$$
#
# Total area and total volume:
#
# $$A = \sum_{i=0}^{\infty} s_i, \tag{13}$$
#
# $$v = v_0 \sum_{i=0}^{\infty} i\,s_i, \tag{14}$$
#
# $$\frac{v}{v_\mathrm{m}} = \frac{\sum_{i=0}^{\infty} i s_i}{\sum_{i=0}^{\infty} s_i}. \tag{15}$$
#
# The two simplifying assumptions, the whole physical content of the model:
#
# $$E_2 = E_3 = \dots = E_i = E_\mathrm{L}, \tag{16}$$
#
# $$\frac{b_2}{a_2} = \frac{b_3}{a_3} = \dots = \frac{b_i}{a_i} = g. \tag{17}$$
#
# They give $s_1 = y s_0$ with $y = (a_1/b_1)p\,e^{E_1/RT}$ (18), $s_2 = x s_1$
# with $x = (p/g)e^{E_\mathrm{L}/RT}$ (19), and hence $s_i = x^{i-1}s_1 = c x^i s_0$
# (20)–(21) with
#
# $$c = \frac{y}{x} = \frac{a_1 g}{b_1}\,e^{(E_1-E_\mathrm{L})/RT}. \tag{22}$$
#
# Summing the two geometric series (23)–(25),
#
# $$\frac{v}{v_\mathrm{m}} = \frac{c\,x}{(1-x)(1-x+c\,x)}, \tag{26}$$
#
# $$x = p/p_0, \tag{27}$$
#
# $$v = \frac{v_\mathrm{m}\,c\,p}{(p_0-p)\bigl[1+(c-1)(p/p_0)\bigr]}, \tag{28}$$
#
# which is eq. (A) rearranged. Its low-pressure limit is Langmuir's (29). If no
# more than $n$ layers can build up:
#
# $$v = \frac{v_\mathrm{m}\,c\,x}{1-x}\;\frac{1-(n+1)x^n+n\,x^{n+1}}{1+(c-1)x-c\,x^{n+1}}, \tag{B}$$
#
# with the two limits the paper names: $n=1$ gives Langmuir, $n=\infty$ gives (A).
# At $n=1$ eq. (B) is also written as
#
# $$\frac{p}{v} = \frac{p_0}{c\,v_\mathrm{m}} + \frac{p}{v_\mathrm{m}}, \tag{E}$$
#
# the form used for charcoal, the one adsorbent in the paper that gave no S-shaped
# isotherm. Footnote 18 adds two further generalisations, (C) with a distinct
# first-layer packing and second-layer heat, and (D) with a distribution of
# capillary sizes; **neither is exercised numerically in the paper and neither is
# used here**, which the page states rather than leaving them to look supported.
#
# Finally, footnote 16 supplies the reading of $c$ used throughout the tables:
#
# $$c = \frac{a_1 b_2}{b_1 a_2}\,e^{(E_1-E_\mathrm{L})/RT}
# \;\approx\; e^{(E_1-E_\mathrm{L})/RT},
# \qquad E_1-E_\mathrm{L} = 2.303\,RT\log c,$$
#
# on the argument that $a_1b_2/b_1a_2$ "will not differ much from unity".
#
# **And the direction of that arrow is the whole story about $c$ on this page.**
# Journal page 313 says which quantity is derived from which: *"From them the
# values of $v_\mathrm{m}$ and $c$ can be evaluated as explained in the previous
# section. From $c$ one can obtain an approximate value for $E_1-E_\mathrm{L}$.$^{16}$"*
# Every $E_1-E_\mathrm{L}$ printed in this paper is therefore $2.303\,RT\log c$
# **with the prefactor already set to 1 by the authors**, so the prefactor is not
# a measurable quantity anywhere in the paper, and inverting footnote 16 returns
# the authors' own fitted $c$ rather than an estimate of it. Section 8 shows what
# *can* be measured in its place, and it is a genuine inconsistency.

# %% [markdown]
# ## Parameters and assumptions
#
# **Everything is pointwise algebra.** There is no grid, no time step and no
# transport anywhere on this page, so items 7 and 8 of the page contract (axis
# refinement, boundary reads) have nothing to act on. The two systems that *are*
# solved, the polarization dipole ladder and the BET layer-occupancy ladder, are
# finite algebraic systems in a *discrete layer index*, not discretisations of a
# continuum, so "refining" them means adding layers, and both are swept over layer
# count in the Validation section for exactly that reason.
#
# **Units, as printed.** Volumes adsorbed are cm³ at STP, per gram of adsorbent
# unless a figure says otherwise; energies are cal/mol; surfaces m²/g except
# Palmer and Clark's, which the paper gives in cm²/g. $R = 1.9872$ cal mol⁻¹ K⁻¹
# and $V_\mathrm{STP} = 22414$ cm³/mol are ours, not the paper's, it prints
# neither.
#
# **Temperatures.** Table I is at 90.1 K, stated. Elsewhere the paper prints
# Celsius, and −183 °C and −195.8 °C are converted with 273.16 K, giving 90.16 and
# 77.36 K; nothing on the page is sensitive at that level and the one place it
# could matter (the reconstruction of $c$) already carries a much larger admitted
# uncertainty.
#
# **The assumptions that carry the physics**, all the paper's:
#
# - eq. (16), every layer above the first has the heat of liquefaction;
# - eq. (17), and the same evaporation–condensation ratio;
# - $a_1$, $b_1$ and $E_1$ independent of how much is already in the first layer,
#   which the authors themselves say fails below $p/p_0 \approx 0.05$, "equation
#   (A) breaks down for the most active points on the surface";
# - for eq. (B), one single $n$ for the whole surface, relaxed only in eq. (D),
#   which is never evaluated.
#
# **Two reconstructions, both labelled as such wherever used.**
#
# 1. **$c$ from $E_1-E_\mathrm{L}$.** No table in the paper prints $c$. Where this
#    page needs $c$ it inverts footnote 16 with the prefactor set to 1, which is
#    the *same* operation the authors performed in the other direction, since
#    journal page 313 says the $E$ column was obtained *from* $c$. The inversion
#    therefore recovers the authors' own fitted $c$, and the only thing lost is
#    that the $E$ column is printed as integers: $\pm0.5$ cal at 90.1 K is
#    $\pm0.28$ % in $c$. What the inversion does **not** do is make $c$ a measured
#    quantity, it is the authors' fit, and section 8 shows the paper's own two
#    routes to one particular $c$ disagreeing by a factor 2.24. Every conclusion
#    drawn from a reconstructed $c$ is checked for sensitivity across *that*
#    factor, and section 8 reports which of them survives it and which does not.
# 2. **Molecular cross-sections.** Tables III and V print $v_\mathrm{m}$ and the
#    surfaces derived from it, but the cross-sections themselves are in Emmett and
#    Brunauer (1937), which is **not on disk and was not consulted**. They are
#    recovered here by inverting the two columns, and the recovery is then checked
#    across two independent tables.
#
# **What is deliberately not attempted.** Fig. 6 (ethyl chloride on charcoal at
# three temperatures) is a fully specified eq. (E) prediction, but reproducing it
# needs the vapour pressure of ethyl chloride at −15.3, 0 and 20 °C, which the
# paper does not print; it is scoped out rather than closed with an outside
# correlation. The Darco statement ($n = 2.2$ from $x = 0.02$ to $0.93$, no point
# off by 4 %) refers to an isotherm that appears nowhere in the paper, in a table
# or otherwise, and is likewise scoped out.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.17.4"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "J1.3-bet-multimolecular-adsorption"  # equals the page directory

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
# Imports. sympy is a Colab default; nothing else beyond pymrm is needed.
import itertools
from fractions import Fraction

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sp
from IPython.display import display

from gallery_utils import Claims, cite_data, load_data, load_meta
from pymrm.checks import observed_orders
from model import (C_eq1b, C_eq1c, X_EDGE, band, c_star_root, d2_numerator, dipole_ladder,
                   inflection_ratio, ladder_v_over_vm, n_for_target, rounding_envelope,
                   shortcut_ratio, sigma_A2, spread, to5, to10, unique_root, v_eqB,
                   v_over_vm_A, x_inflection)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
pd.set_option("display.width", 130)
claims = Claims(PAGE.split("-")[0])

R_CAL = 1.98720        # cal/(mol K), ours; the paper prints no value
V_STP = 22414.0        # cm3/mol at STP, ours
N_AV = 6.02214076e23   # 1/mol, ours
T_TABLE1 = 90.1        # K, printed in the Table I heading
KELVIN = 273.16        # used only to convert the printed Celsius temperatures

# %% [markdown]
# ## The data
#
# Six files, all transcribed for this page from renders of the original at its
# **native 300 ppi** (`pdfimages -list` reports every page as CCITT-G4 bilevel at
# 300 ppi, so a larger render only interpolates), each numeric cropped and re-read
# at digit scale. The PDF's text layer was not used for any digit. **No other
# page's dataset is loaded**, so none of the cross-page reconciliation obligations
# apply here.
#
# | file | what it is |
# |---|---|
# | `bet-1938-table1-nitrogen-90K.csv` | Table I, $v_\mathrm{m}$, point B and $E_1-E_\mathrm{L}$ for nitrogen at 90.1 K on twelve adsorbents |
# | `bet-1938-table2-e1el-gases.csv` | Table II, $E_1-E_\mathrm{L}$ for four gases on two adsorbents |
# | `bet-1938-table3-silica-gel.csv` | Table III, seven gases on one silica gel, with both surface columns |
# | `bet-1938-table4-so2-silica-gel.csv` | Table IV, SO₂ on silica gel at six temperatures (McGavack & Patrick's isotherms) |
# | `bet-1938-table5-charcoal.csv` | Table V, eight isotherms on charcoal |
# | `bet-1938-printed-claims.csv` | every scalar in prose, in a footnote or in a figure's constant block that this page checks |
#
# **Tier 6, and it matters here more than usual.** Not one number in any of these
# tables is a measurement. Every entry is a constant the authors *derived* from an
# isotherm, $v_\mathrm{m}$ and $c$ from their own eq. (A) plot, point B read by
# eye, $E_1-E_\mathrm{L}$ from $c$ through footnote 16, the surface columns from
# $v_\mathrm{m}$ and an external cross-section. Reproducing one of these numbers
# is reproduction, never validation, and the page says so at each such comparison.
#
# ### Fit or test: the split this page keeps
#
# | quantity | fitted to what | so agreement with it is |
# |---|---|---|
# | $v_\mathrm{m}$, $c$ | the isotherm itself, over $0.05<p/p_0<0.35$ | **a fit** |
# | the eq. (A) straight line inside that window | itself | **a fit**, and the page never reports it |
# | $n$ in eq. (B) | the isotherm above the window, by trial and error | **a fit** |
# | $E_1-E_\mathrm{L}$ | nothing further, a transform of $c$ | **the same fit, re-expressed** |
# | agreement of surface area **between gases** | nothing, cross-sections come from bulk densities | **a test** |
# | point B against $v_\mathrm{m}$ | nothing, point B is read off the raw isotherm | **a test, but a weak one** (§6) |
# | $v_\mathrm{m}(77.3\,\mathrm{K})$ predicted from the 90.1 K isotherm | nothing | **a test** |
# | BET area against Palmer & Clark's HF-dissolution area | nothing, a different measurement entirely | **a test** |
#
# The single most important line in that table is the third-from-last. It is the
# reason this page exists in the form it does: the polarization theory can fit an
# isotherm, and cannot be asked any of the questions in the "test" rows, because
# it contains no $v_\mathrm{m}$.
#
# **One residual risk in that row, stated rather than buried.** That the
# cross-sections come from bulk densities rests on footnote 17 and on Emmett and
# Brunauer (1937), **which is not on disk and was not consulted**. If those
# cross-sections had in fact been adjusted to make the gases agree, the test would
# be circular. Three things bound that risk without the 1937 paper. The numbers
# recovered below are *not* round or tuned-looking (section 5 prints them, on
# solid and on liquid packing). One single set reproduces two *different*
# adsorbents' tables, which §5 verifies by recovering it from each independently.
# And, the only one of the three that is a **number rather than an argument**,
# §5 builds the worst case explicitly: a σ set tuned to flatten the silica-gel
# areas exactly, applied to charcoal, and the resulting spread compared with what
# the printed set achieves on the same rows. The bound it gives is a factor of
# about two, which is real but modest, and it is quoted with the basis it holds on;
# on the butane-inclusive basis the test has no power at all, which §5 also
# prints.

# %%
t1 = load_data("bet-1938-table1-nitrogen-90K.csv", page=PAGE)
t2 = load_data("bet-1938-table2-e1el-gases.csv", page=PAGE)
t3 = load_data("bet-1938-table3-silica-gel.csv", page=PAGE)
t4 = load_data("bet-1938-table4-so2-silica-gel.csv", page=PAGE)
t5 = load_data("bet-1938-table5-charcoal.csv", page=PAGE)
printed = load_data("bet-1938-printed-claims.csv", page=PAGE)
P = dict(zip(printed.key, printed.printed.astype(float)))

print(cite_data(load_meta("bet-1938-table1-nitrogen-90K.csv", page=PAGE)))
print(f"\nTable I ({len(t1)} rows) | Table II ({len(t2)}) | Table III ({len(t3)}) | "
      f"Table IV ({len(t4)}) | Table V ({len(t5)}) | printed claims ({len(printed)})")
display(t1)

# %% [markdown]
# ### Five transcription checks the tables give for free
#
# Table II's silica-gel row repeats four values that also appear in Table III, and
# its catalyst-954 nitrogen entry also appears in Table I. Those five coincidences
# are not physics; they are a check that the transcription of three separate
# tables is mutually consistent, and they are the cheapest thing on this page.

# %%
xchk = []
for g, T in [("N2", -183), ("A", -183), ("CO2", -78), ("C4H10", 0)]:
    a = t2.loc[(t2.substance == "Silica gel") & (t2.gas == g), "E1_minus_EL_cal_per_mol"].item()
    b = t3.loc[(t3.gas == g) & (t3.temperature_C == T), "E1_minus_EL_cal_per_mol"].item()
    xchk.append(("Table II silica gel", f"Table III {g} at {T} C", a, b, a - b))
a = t2.loc[(t2.substance.str.contains("954")) & (t2.gas == "N2"), "E1_minus_EL_cal_per_mol"].item()
b = t1.loc[t1.substance == "Fe-Al2O3 catalyst 954", "E1_minus_EL_cal_per_mol"].item()
xchk.append(("Table II catalyst 954", "Table I catalyst 954", a, b, a - b))
xchk = pd.DataFrame(xchk, columns=["source A", "source B", "A", "B", "A-B"])
display(xchk)
TRANSCRIPTION_MAX_DEV = float(xchk["A-B"].abs().max())
print(f"worst disagreement across the five repeated cells: {TRANSCRIPTION_MAX_DEV:.0f} cal/mol")


# %% [markdown]
# ## PyMRM implementation
#
# Two algebraic ladders, and a set of root-finds. Nothing here is a PDE, and the
# page does not pretend otherwise, see *What pymrm adds*.
#
# **1. The BET layer-occupancy ladder** (eqs. 10–12, 18–21). The closed forms (26)
# and (B) come from summing two geometric series. To get a route to
# $v/v_\mathrm{m}$ that shares *no algebra at all* with them, the equilibria are
# solved as a system for the occupancies $\{s_i\}$ and the two sums in eq. (15)
# are then formed **term by term**:
#
# $$r_1 = s_1 - c\,x\,s_0, \qquad r_i = s_i - x\,s_{i-1}\;\;(i = 2\dots n), \qquad s_0 \equiv 1,$$
#
# solved with `newton` and a `NumJac((n, 1), axes_diagonals=[0])` Jacobian. The
# shape is `(n, 1)` and not `(n,)`: with a bare 1-D shape `NumJac` treats the layer
# index as a field index and builds a dense $n\times n$ Jacobian. `axes_diagonals=[0]`
# is correct here and is one of the rare cases where it is, the residual itself
# reads the neighbouring layer, which is exactly the condition the house rule
# names. The step tolerance is `1e-12` on the occupancies scaled by $c\,x$, with
# `stagnation=True`, and the call raises unless `newton` reports success; the
# continuous-$n$ solve of section 4 is checked the same way.
#
# **2. The polarization dipole ladder** (eq. 1a). Same idea applied to the
# *rival* theory: instead of quoting its closed form (1b), the recursion
# $\mu_i = k(\mu_{i-1}+\mu_{i+1})$ is assembled as a tridiagonal system with
# $\mu_0 = 1$ and $\mu_{N+1} = 0$ and solved directly. The decay ratio is then
# *measured* off the solution. For $k \le 1/2$ this must reproduce (1b); for
# $k > 1/2$ it must fail to settle on any ratio at all, and the page shows it
# doing so.
#
# **3. Root-finds, never sampled sweeps.** The inflection point of eq. (26), the
# $c$ at which the inflection value is largest, the $c$ at which the paper's
# single-point rule reaches exactly 5 % error, and the continuous $n$ that
# reproduces a printed percentage, every one is a root-find: `pymrm.checks.find_roots`
# on the physical range, raising unless it finds exactly one root. Two of them also
# have closed forms, derived symbolically here, and the closed form and the
# root-find are reported against each other. The maximum over $c$ is the root of a
# central-difference derivative; its step is a parameter, and section 6 measures
# what it costs.

# %%
# The forms and both ladders live in model.py (imported above); one call of each here.
_v, _sol = ladder_v_over_vm(0.3, 100.0, 6)
print(f"layer ladder, n = 6, x = 0.3, c = 100: v/v_m = {_v:.12f} after {_sol.nit} newton "
      f"step(s); eq. (B) gives {float(v_eqB(0.3, 100.0, 6)):.12f}")
print(f"dipole ladder, k = 0.01, 50 layers: mu_1/mu_0 = {dipole_ladder(0.01, 50)[0]:.12f}; "
      f"eq. (1b) gives {C_eq1b(0.01):.12f}")

# %% [markdown]
# ## Results
#
# ### 1. The rival theory, refuted on its own constants
#
# The paper's chain, reproduced input by input. Nothing here is fitted to
# anything: $\alpha/r^3$ comes from a refractive index and a solid density, $d$
# comes from three lattice sums, and $C$ follows.
#
# Two by-products worth printing. Eq. (7) is not stated to follow from (5) and (6);
# it is checked symbolically that it does. And the paper never prints the density
# of solid argon it used, only the reference, so inverting eq. (7) against the
# printed $\alpha/r^3 = 0.029$ says what density that was, which is a cheap test of
# whether the printed 0.029 is even self-consistent.

# %%
# --- eq. (7) follows from (5) and (6): symbolic, not asserted
n_ref, v_g, v_s = sp.symbols("n_ref v_g v_s", positive=True)
alpha_sym = (n_ref - 1) * v_g / (2 * sp.pi)               # eq. (5)
r3_sym    = sp.sqrt(2) * v_s                              # eq. (6)
EQ7_RESID = sp.simplify(alpha_sym / r3_sym
                        - (n_ref - 1) / (2**sp.Rational(3, 2) * sp.pi) * v_g / v_s)
print(f"eq. (7) from eqs. (5) and (6):  residual = {EQ7_RESID}  ->  identity: {EQ7_RESID == 0}")

# --- the printed chain
A_R3 = P["polar_alpha_over_r3"]
LATTICE_NET = P["polar_lattice_adj_quad"] - P["polar_lattice_adj_cube"]     # 8.357 - 8.823
D_CONST = LATTICE_NET / (1.0 + P["polar_lattice_self"] * A_R3)             # footnote 10
C_POL = D_CONST * A_R3                                                     # eq. (9)
K_POL = abs(C_POL)                                                         # k = 0.35 alpha/r^3
K1_POL = C_POL**2                                                          # eq. (4a)

chain = pd.DataFrame([
    ("8.357 - 8.823 -> the coefficient of (mu_i-1 + mu_i+1)", LATTICE_NET, -P["polar_net_adj"]),
    ("d = -0.466 / (1 + 11.1 alpha/r^3)",                     D_CONST,     P["polar_d"]),
    ("C = d (alpha/r^3)",                                     C_POL,       P["polar_C"]),
    ("K1 = C^2",                                              K1_POL,      P["polar_K1"]),
], columns=["step", "recomputed", "printed"])
display(chain)

# --- what density of solid argon the printed 0.029 implies
VG_OVER_VS = A_R3 * 2**1.5 * np.pi / P["polar_n_minus_1"]
V_SOLID = V_STP / VG_OVER_VS
RHO_SOLID = 39.948 / V_SOLID
R_NN = (np.sqrt(2.0) * V_SOLID / N_AV)**(1.0 / 3.0) * 1e8
print(f"\ninverting eq. (7) against the printed alpha/r^3 = {A_R3}:")
print(f"  v_g/v_s = {VG_OVER_VS:.4f}  ->  molar volume of solid argon {V_SOLID:.3f} cm3/mol")
print(f"  -> density {RHO_SOLID:.4f} g/cm3, nearest-neighbour distance {R_NN:.4f} Angstrom")
print("  (both are recovered, not printed; the paper cites Int. Crit. Tables I p. 103,")
print("   which is NOT on disk and was not consulted. They are stated here as the")
print("   inputs the printed 0.029 implies, and they are physically ordinary ones.)")

# %% [markdown]
# **A printed grouping that does not parse, proved from the paper's own
# arithmetic.** Footnote 10 sets the adjacent-layer lattice sum as
# `(8.357 - 8.823/r^3)`. As printed that is dimensionally inconsistent, a pure
# number minus a number over $r^3$, and it does not produce the coefficient used
# two lines later. The intended grouping is $(8.357-8.823)/r^3$, and the cell above
# shows $8.357-8.823$ is exactly $-0.466$, which is exactly the coefficient of the
# next display. **Reported, not repaired**, and the CSV stores the two lattice sums
# and the 0.466 as separate keys so the misprint can be proved rather than
# asserted.
#
# ### The comparison the whole refutation turns on
#
# Bradley obtained $k$ by fitting, and used eq. (1c) to turn it into $C$. Compare
# that $C$ with the one just computed.

# %%
k_brad = P["bradley_k_Al2SO4"]
C_BRAD = C_eq1c(k_brad)
C_FACTOR = C_BRAD / abs(C_POL)
E_FACTOR = C_FACTOR**2                       # binding energy goes as C^2, eq. (2)
K_LIMIT_RATIO = k_brad / P["bradley_k_upper_limit"]
K_LO, K_HI = 1e-12, 0.49999   # k on which eq. (1c) is scanned: up to the eq. (1b) limit 1/2
K_ADMISSIBLE = unique_root(lambda kk: C_eq1c(kk) - abs(C_POL), K_LO, K_HI)

cmp = pd.DataFrame([
    ("computed from the polarizability (eqs. 5-9)", K_POL,  abs(C_POL), K_POL**2),
    ("Bradley's fit, CuSO4 (eq. 1c)", P["bradley_k_CuSO4"], C_eq1c(P["bradley_k_CuSO4"]),
     C_eq1c(P["bradley_k_CuSO4"])**2),
    ("Bradley's fit, Al2(SO4)3 (eq. 1c)", k_brad, C_BRAD, C_BRAD**2),
], columns=["source of k", "k", "|C|", "C^2 = K1  (energy ratio per layer)"])
display(cmp)

print(f"C from Bradley's k = {k_brad} is {C_BRAD:.6f}; the paper prints {P['bradley_C_implied']:g}")
print(f"factor in the dipole ratio C  : {C_FACTOR:.4f}")
print(f"factor in the energy ratio C^2: {E_FACTOR:.1f}")
print(f"\nk that eq. (1c) needs to give the COMPUTED |C|: {K_ADMISSIBLE:.8f}")
print(f"  Bradley's k is {k_brad / K_ADMISSIBLE:.2f} times larger.")
print(f"\nand eq. (1b) is real only for k <= {P['bradley_k_upper_limit']:g}:")
print(f"  Bradley's two k exceed that limit by "
      f"{(P['bradley_k_CuSO4']/P['bradley_k_upper_limit']-1)*100:.1f} % and "
      f"{(K_LIMIT_RATIO-1)*100:.1f} %; 1 - 4k^2 = "
      f"{1-4*P['bradley_k_CuSO4']**2:+.4f} and {1-4*k_brad**2:+.4f}")

# %% [markdown]
# ### 2. The same conclusion without eq. (1b): solve the recursion
#
# Everything above rests on closed forms, (1b) and (1c), that the paper quotes.
# A closed form quoted from a rival's paper is exactly the kind of thing that
# should not be load-bearing, so the recursion (1a) itself is assembled and solved
# here, and the decay ratio *measured* off the solution.
#
# For $k \le 1/2$ this must return eq. (1b). For $k > 1/2$ the recursion has no
# decaying solution, so there is nothing for it to return, and what it does
# instead is the interesting part: the answer depends on where the ladder is
# truncated, and the moments change sign along it. **That is the strongest
# statement the page makes against Bradley's $k$**, and it needs neither (1b) nor
# (1c).

# %%
N_LADDER = 400
mu_phys = dipole_ladder(K_POL, N_LADDER)
LADDER_C = float(mu_phys[0])
LADDER_VS_1B = abs(LADDER_C / C_eq1b(K_POL) - 1.0)
print(f"physical k = {K_POL:.8f}")
print(f"  measured mu_1/mu_0     = {LADDER_C:.12f}")
print(f"  eq. (1b)               = {C_eq1b(K_POL):.12f}   relative deviation {LADDER_VS_1B:.3e}")
print(f"  geometric? mu_2/mu_1   = {mu_phys[1]/mu_phys[0]:.12f}")
print(f"  and k ~ C to within k^2: |C/k - 1| = {abs(C_eq1b(K_POL)/K_POL - 1):.4e}, "
      f"k^2 = {K_POL**2:.4e}")

N_SWEEP = (100, 200, 300, 400, 500, 600)   # ladder lengths of the truncation sweep


def trunc_spread(k, ns=N_SWEEP):
    """Relative spread of the measured mu_1/mu_0 over the ladder lengths ns."""
    m1 = np.array([dipole_ladder(k, N)[0] for N in ns])
    return float(np.ptp(m1) / abs(np.mean(m1)))


trunc = pd.DataFrame([(N, dipole_ladder(K_POL, N)[0], dipole_ladder(k_brad, N)[0]) for N in N_SWEEP],
                     columns=["layers N", "mu_1/mu_0 at the physical k", "mu_1/mu_0 at Bradley's k"])
display(trunc)
LADDER_SPREAD_PHYS = trunc_spread(K_POL)
LADDER_SPREAD_BRAD = trunc_spread(k_brad)
print(f"relative spread over N:  physical k {LADDER_SPREAD_PHYS:.3e}   Bradley's k {LADDER_SPREAD_BRAD:.4f}")
sgn = np.sign(dipole_ladder(k_brad, N_LADDER)[:8]).astype(int).tolist()
print(f"sign of mu_1..mu_8 at Bradley's k: {sgn}  -- an oscillation, not a decaying film")

# %%
fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
i = np.arange(1, 21)
ax[0].semilogy(i, np.abs(dipole_ladder(K_POL, N_LADDER)[:20]), "o-", ms=4,
               label=f"computed k = {K_POL:.5f}")
ax[0].semilogy(i, np.abs(C_eq1b(K_POL)**i), "k--", lw=1, label="eq. (1b): $C^i$")
ax[0].set(xlabel="layer $i$", ylabel=r"$|\mu_i / \mu_0|$",
          title="Polarization ladder at the computed $k$")
ax[0].legend(fontsize=8)
ax[1].plot(i, dipole_ladder(k_brad, N_LADDER)[:20], "s-", ms=4, color="C3",
           label=f"Bradley's k = {k_brad}")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set(xlabel="layer $i$", ylabel=r"$\mu_i / \mu_0$",
          title="Polarization ladder at Bradley's $k$ (no decaying solution)")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()
print("Left: at the k the polarizability gives, the ladder decays geometrically and the")
print("      measured ratio is eq. (1b) to round-off. Right: at Bradley's k it oscillates,")
print("      changes sign, and depends on where the ladder is cut off.")

# %% [markdown]
# **What this does and does not establish.** It establishes that *the
# polarization theory as Brunauer, Emmett and Teller restate it* cannot carry the
# binding energy Bradley's fit assigns it, using only inputs the paper prints. It
# does **not** establish anything about de Boer, Zwicker or Bradley's own papers:
# none of the three is on disk and none was consulted. The rival is being judged in
# the form its opponent gives it, and that limitation is real. What makes the
# judgement worth something anyway is that the inputs, a refractive index, a solid
# density, three lattice sums, are checkable physics rather than assertions, and
# two of them are reproduced above from independent directions.
#
# ### 3. The BET algebra: five identities, checked symbolically
#
# Every step from the layer equilibria to the lettered equations. These are
# identities, they are exactly zero, and the page treats them as transcription
# checks and nothing more, a wrong exponent anywhere would show up here, and
# nothing else on this page would notice.

# %%
x, c, n, vm, p, p0 = sp.symbols("x c n v_m p p_0", positive=True)
eq26 = c * x / ((1 - x) * (1 - x + c * x))
eqB  = c * x / (1 - x) * (1 - (n + 1) * x**n + n * x**(n + 1)) / (1 + (c - 1) * x - c * x**(n + 1))
eq28 = vm * c * p / ((p0 - p) * (1 + (c - 1) * (p / p0)))

ident = []
ident.append(("(28) rearranges to (A)",
              sp.simplify(p / (eq28 * (p0 - p)) - (1 / (vm * c) + (c - 1) / (vm * c) * p / p0))))
ident.append(("(B) at n -> infinity is (26)",
              sp.simplify(eqB.subs({x**n: 0, x**(n + 1): 0}) - eq26)))
ident.append(("(B) at n = 1 is the Langmuir form",
              sp.simplify(eqB.subs(n, 1) - c * x / (1 + c * x))))
v_n1 = (vm * eqB.subs(n, 1)).subs(x, p / p0)
ident.append(("(B) at n = 1 rearranges to (E)",
              sp.simplify(p / v_n1 - (p0 / (c * vm) + p / vm))))
ident.append(("(7) follows from (5) and (6)", EQ7_RESID))
ident = pd.DataFrame(ident, columns=["identity", "residual"])
display(ident)
SYMBOLIC_MAX_RESID = float(max(abs(float(r)) for r in ident.residual))
print(f"largest symbolic residual: {SYMBOLIC_MAX_RESID:g}  (exactly zero -- these are identities)")

# %% [markdown]
# ### 4. Eq. (B) on every branch, and the two percentages the paper prints
#
# The identities above are worthless as a check on the *numbers*, so eq. (B) is
# now evaluated the hard way: the layer equilibria (10)–(12) are solved for the
# occupancies and eq. (15)'s two sums formed term by term, with no geometric
# series anywhere. Exercised on **every branch the paper uses**, $n = 1$
# (charcoal), $n = 2$ (Darco, quoted), $n = 5,6,7$ (Fig. 3), $n = 3.5$ is
# non-integer and belongs to Table IV, and $n \to \infty$ (eq. A), and at four
# relative pressures, two inside the fitting window and two outside it.
#
# Then the printed claim itself. Journal page 316: with $n = 6$ against the $n$
# that actually fits, "the theoretical curve will be in error to the extent of
# $+5\%$ at $x = 0.58$, and $-7\%$ at $x = 0.72$". Fig. 3 prints $v_\mathrm{m}$ and
# $c$ for that very isotherm, so both percentages are computable.

# %%
N_BRANCH = (1, 2, 5, 6, 7, 20, 60)     # every branch the paper uses, and two long ladders
X_BRANCH = (0.05, 0.30, 0.58, 0.72)    # two inside the fitting window, two outside


def ladder_worst(closed=v_eqB, c_ladder=None, n_shift=0, ns=N_BRANCH, xs=X_BRANCH):
    """Worst relative gap between the ladder solve and a closed form over (n, x).

    closed: the closed form compared with; c_ladder: the c the ladder is run at
    (default Fig. 3's); n_shift: the closed form evaluated at n + n_shift."""
    cl = P["fig3_c"] if c_ladder is None else c_ladder
    return max(abs(ladder_v_over_vm(xx, cl, nn)[0] / float(closed(xx, P["fig3_c"], nn + n_shift)) - 1)
               for nn in ns for xx in xs)


rows = []
for nn in N_BRANCH:
    for xx in X_BRANCH:
        lad, sol = ladder_v_over_vm(xx, P["fig3_c"], nn)      # raises unless newton succeeds
        cf = float(v_eqB(xx, P["fig3_c"], nn))
        rows.append((nn, xx, lad, cf, abs(lad / cf - 1), sol.nit))
branch = pd.DataFrame(rows, columns=["n", "x", "ladder solve", "eq. (B) closed form",
                                     "rel. dev.", "Newton its"])
display(branch.round({"ladder solve": 9, "eq. (B) closed form": 9}))
LADDER_VS_EQB = float(branch["rel. dev."].max())
print(f"worst relative deviation over {len(branch)} (n, x) pairs: {LADDER_VS_EQB:.3e}, "
      f"every newton solve reporting success in {branch['Newton its'].min()} to "
      f"{branch['Newton its'].max()} steps")

# eq. (A) as the large-n limit -- the paper's own claim, quantified
print("\nthe paper: \"when x has a small value, and n is as large as 4 or 5, equation (A)")
print("becomes a very good approximation to (B)\".  Measured, with Fig. 3's c:")
for xx in (0.05, 0.35, 0.58, 0.72):
    for nn in (4, 5, 6):
        rel = float(v_eqB(xx, P["fig3_c"], nn) / v_over_vm_A(xx, P["fig3_c"]) - 1) * 100
        print(f"   x = {xx:<5} n = {nn}:  eq.(B)/eq.(A) - 1 = {rel:+9.4f} %")

# %%
VM3, C3 = P["fig3_vm"], P["fig3_c"]
ERR_058 = float(v_eqB(P["eqB_x_n5"], C3, 6, VM3) / v_eqB(P["eqB_x_n5"], C3, 5, VM3) - 1) * 100
ERR_072 = float(v_eqB(P["eqB_x_n7"], C3, 6, VM3) / v_eqB(P["eqB_x_n7"], C3, 7, VM3) - 1) * 100
pr = pd.DataFrame([
    (P["eqB_x_n5"], "n = 6 against n = 5", v_eqB(P["eqB_x_n5"], C3, 5, VM3),
     v_eqB(P["eqB_x_n5"], C3, 6, VM3), ERR_058, P["eqB_err_at_058"]),
    (P["eqB_x_n7"], "n = 6 against n = 7", v_eqB(P["eqB_x_n7"], C3, 7, VM3),
     v_eqB(P["eqB_x_n7"], C3, 6, VM3), ERR_072, P["eqB_err_at_072"]),
], columns=["x", "comparison", "reference v (cc)", "n = 6 v (cc)", "recomputed %", "printed %"])
display(pr.round(4))

# the continuous n each printed percentage implies: a newton root-find (n_for_target
# raises unless newton reports success)
def n_from_pct(xx, ref_n, pct, cc=None):
    """The n whose eq. (B) value is pct % off the ref_n branch at xx (c of Fig. 3)."""
    cc = C3 if cc is None else cc
    return n_for_target(xx, cc, (1 + pct / 100) * float(v_eqB(xx, cc, ref_n, 1.0)), 6.0)


N_FROM_058 = n_from_pct(P["eqB_x_n5"], 5, P["eqB_err_at_058"])
N_FROM_072 = n_from_pct(P["eqB_x_n7"], 7, P["eqB_err_at_072"])
print(f"n giving exactly +5 % over the n = 5 curve at x = 0.58 : {N_FROM_058:.5f}")
print(f"n giving exactly -7 % under the n = 7 curve at x = 0.72: {N_FROM_072:.5f}")
print(f"the two printed percentages are mutually consistent with a single n in "
      f"[{min(N_FROM_058, N_FROM_072):.3f}, {max(N_FROM_058, N_FROM_072):.3f}], i.e. the paper's")
print(f"'average value of n = 6' is right to {max(abs(N_FROM_058/6-1), abs(N_FROM_072/6-1))*100:.2f} %.")
print("Non-integer n is the paper's own device, not ours: Table IV's footer reads n = 3.5.")

# %%
xs = np.linspace(0.02, 0.80, 400)
fig, ax = plt.subplots(figsize=(6.6, 4.2))
for nn, ls in [(1, ":"), (5, "-"), (6, "-"), (7, "-")]:
    ax.plot(xs, v_eqB(xs, C3, nn, VM3), ls, lw=1.5, label=f"eq. (B), n = {nn}")
ax.plot(xs, v_over_vm_A(xs, C3) * VM3, "k--", lw=1.2, label=r"eq. (A), $n\to\infty$")
ax.axhline(VM3, color="0.5", lw=0.8)
ax.text(0.03, VM3 * 1.03, "$v_m$", fontsize=8, color="0.4")
for xx in (P["eqB_x_n5"], P["eqB_x_n7"]):
    ax.axvline(xx, color="0.75", lw=0.8, ls="-.")
ax.set(xlabel="$p/p_0$", ylabel="v adsorbed, cm$^3$ STP", ylim=(0, 520),
       title=f"Eq. (B) on its branches, at Fig. 3's own constants\n$v_m$ = {VM3} cc, c = {C3}")
ax.legend(fontsize=8, loc="upper left"); fig.tight_layout(); plt.show()
print("Constants from the typeset block inside the Fig. 3 frame. NO CURVE IS DIGITISED:")
print("this is eq. (B) evaluated, not the paper's drawing traced.")


# %% [markdown]
# ### 5. The evidence that is not a fit
#
# $v_\mathrm{m}$ is fitted, once per isotherm. What is *not* fitted is that seven
# different gases must give the same surface. Turning $v_\mathrm{m}$ into an area
# needs a cross-sectional area per molecule, and those come from the density of the
# solidified and of the liquefied gas, bulk properties, with no adsorption in
# them. Nothing was adjusted to make the seven agree.
#
# **First, recover the cross-sections and check them across two tables.** Tables
# III and V were computed with the same set, so inverting each independently must
# return the same numbers. This also fixes something the page would otherwise have
# to take on faith: Emmett and Brunauer (1937), where the cross-sections are
# tabulated, is not on disk and was not consulted.

# %%
def sigma_cross(t3_, t5_, v_stp=V_STP, du_S=0.5):
    """Cross-sections recovered from Tables III and V independently, compared pair by pair.

    v_stp enters the recovery; du_S is half a unit in the last printed place of the
    surfaces (integers: 0.5). The rounding envelope of a ratio of two sigma follows
    from the printed quantisation alone, since v_stp and N_A cancel in it."""
    sg = lambda t: t.assign(sig_solid=sigma_A2(t.v_m_cc_per_g, t.surface_solid_m2_per_g, v_stp, N_AV),
                            sig_liquid=sigma_A2(t.v_m_cc_per_g, t.surface_liquid_m2_per_g, v_stp, N_AV))
    m = sg(t3_).merge(sg(t5_), on=["gas", "temperature_C"], suffixes=("_silica", "_charcoal"))
    for pk in ("solid", "liquid"):
        m[f"dev_{pk}_%"] = (m[f"sig_{pk}_silica"] / m[f"sig_{pk}_charcoal"] - 1).abs() * 100
        m[f"env_{pk}_%"] = rounding_envelope(m.v_m_cc_per_g_silica, m[f"surface_{pk}_m2_per_g_silica"],
                                             m.v_m_cc_per_g_charcoal, m[f"surface_{pk}_m2_per_g_charcoal"],
                                             du_S=du_S)
        m[f"over_{pk}"] = m[f"dev_{pk}_%"] / m[f"env_{pk}_%"]
    return m


def sigma_stats(m):
    """Worst disagreement %, worst ratio to the envelope, pairs inside their envelope."""
    over = np.r_[m.over_solid.to_numpy(), m.over_liquid.to_numpy()]
    return {"sigma_cross_table_max_pct": float(max(m["dev_solid_%"].max(), m["dev_liquid_%"].max())),
            "sigma_cross_table_max_over_envelope": float(over.max()),
            "sigma_pairs_inside_envelope": int((over <= 1).sum())}


mg = sigma_cross(t3, t5)
display(mg[["gas", "temperature_C", "sig_solid_silica", "sig_solid_charcoal", "dev_solid_%",
            "sig_liquid_silica", "sig_liquid_charcoal", "dev_liquid_%"]].round(4))
SIG = sigma_stats(mg)
SIGMA_CROSS_MAX = SIG["sigma_cross_table_max_pct"]
SIGMA_OVER_ENVELOPE = SIG["sigma_cross_table_max_over_envelope"]
n_over = 2 * len(mg) - SIG["sigma_pairs_inside_envelope"]
SIG_WORST_ROW = mg.loc[mg["over_liquid"].idxmax()]
display(mg[["gas", "temperature_C", "dev_solid_%", "env_solid_%", "over_solid",
            "dev_liquid_%", "env_liquid_%", "over_liquid"]].round(4))
print(f"{len(mg)} isotherms are shared between Tables III and V, giving {2*len(mg)} independent")
print(f"pairs of recovered cross-sections. Worst disagreement: {SIGMA_CROSS_MAX:.4f} %.")
print("\nIS THAT THE ROUNDING? Mostly, but NOT on the row that sets the maximum -- and the")
print("page does not attribute it to rounding, because the arithmetic says it cannot be.")
print(f"  {2*len(mg)-n_over} of the {2*len(mg)} pairs sit inside their own rounding envelope.")
print(f"  {n_over} does not: {SIG_WORST_ROW.gas} at {SIG_WORST_ROW.temperature_C:g} C, "
      f"liquid packing --")
print(f"    {SIG_WORST_ROW['dev_liquid_%']:.4f} % observed against a "
      f"{SIG_WORST_ROW['env_liquid_%']:.4f} % envelope,")
print(f"    a factor {SIG_WORST_ROW['over_liquid']:.4f}. All four inputs "
      f"({SIG_WORST_ROW.v_m_cc_per_g_silica:g}, "
      f"{SIG_WORST_ROW.surface_liquid_m2_per_g_silica:.0f}, "
      f"{SIG_WORST_ROW.v_m_cc_per_g_charcoal:g}, "
      f"{SIG_WORST_ROW.surface_liquid_m2_per_g_charcoal:.0f})")
print("    were re-read at digit scale; this is in the paper, not in the CSV. CAUSE UNKNOWN,")
print("    reported and not repaired. Counted with the printed defects of section 9.")
print("\nWhat the check still establishes: the two tables were computed with ONE set of")
print("cross-sections, and this page has that set without the 1937 paper.")


# %% [markdown]
# **Now the test itself, with its null baseline beside it.** If the
# cross-sections were doing no work, if $v_\mathrm{m}$ were an arbitrary number
# per gas, the areas would scatter as widely as the $v_\mathrm{m}$ do. So the
# $v_\mathrm{m}$ spread *is* the null baseline, and it is printed next to the area
# spread every time.

# %%
k5 = (t5.gas != "C4H10")                     # the paper's own exclusion, journal page 318
sets = [
    ("silica gel (Table III), solid packing",  t3.surface_solid_m2_per_g,  P["t3_solid_mean"],  P["t3_solid_maxdev"]),
    ("silica gel (Table III), liquid packing", t3.surface_liquid_m2_per_g, P["t3_liquid_mean"], P["t3_liquid_maxdev"]),
    ("charcoal (Table V, 7 rows), solid",      t5.surface_solid_m2_per_g[k5],  P["t5_solid_mean"],  P["t5_solid_maxdev"]),
    ("charcoal (Table V, 7 rows), liquid",     t5.surface_liquid_m2_per_g[k5], P["t5_liquid_mean"], P["t5_liquid_maxdev"]),
]
rows = []
for lab, S, pm, pd_ in sets:
    m, d = spread(S)
    rows.append((lab, m, pm, d, pd_))
surf = pd.DataFrame(rows, columns=["set", "mean recomputed", "mean printed",
                                   "max deviation %", "printed %"])
display(surf.round(4))

T3_SOLID_MEAN, T3_SOLID_DEV = spread(t3.surface_solid_m2_per_g)
T3_LIQ_MEAN,   T3_LIQ_DEV   = spread(t3.surface_liquid_m2_per_g)
T5_SOLID_MEAN, T5_SOLID_DEV = spread(t5.surface_solid_m2_per_g[k5])
T5_LIQ_MEAN,   T5_LIQ_DEV   = spread(t5.surface_liquid_m2_per_g[k5])
_, T3_VM_DEV = spread(t3.v_m_cc_per_g)
_, T5_VM_DEV = spread(t5.v_m_cc_per_g[k5])
DISC_SILICA = T3_VM_DEV / T3_LIQ_DEV
DISC_CHARCOAL = T5_VM_DEV / T5_LIQ_DEV
print(f"\nNULL BASELINE, silica gel: the seven fitted v_m span {T3_VM_DEV:.4f} % about their mean.")
print(f"  Converted through cross-sections computed from bulk densities they span "
      f"{T3_LIQ_DEV:.4f} %.")
print(f"  Collapse factor {DISC_SILICA:.4f}.")
print(f"NULL BASELINE, charcoal:   {T5_VM_DEV:.4f} % -> {T5_LIQ_DEV:.4f} %, factor {DISC_CHARCOAL:.4f}.")
print("\nThis is the page's central discriminating number and it is NOT a goodness of fit:")
print("each v_m was fitted to its own isotherm and to nothing else, and the conversion")
print("contains no adsorption. The polarization theory has no v_m and cannot be asked.")

# --- HOW MUCH POWER DOES THE CHARCOAL CHECK ACTUALLY HAVE? Bound the circularity
#     risk with a number instead of a sentence: build the WORST-CASE tuned set --
#     sigma proportional to 1/v_m on silica gel, which makes the silica areas exactly
#     equal by construction -- and see what it does on charcoal. Six (gas, temperature)
#     pairs are shared and non-butane, so the comparison is like for like.
shared = t3.merge(t5, on=["gas", "temperature_C"], suffixes=("_s", "_c"))
sh6 = shared[shared.gas != "C4H10"]
vm_s = sh6.v_m_cc_per_g_s.to_numpy(float); vm_c = sh6.v_m_cc_per_g_c.to_numpy(float)


def tuned_spread(sh, power=1.0):
    """Charcoal area spread with sigma tuned as 1/v_m^power on silica gel (power 1
    flattens the silica areas exactly; the scale of sigma is irrelevant to a spread)."""
    return spread(sh.v_m_cc_per_g_c.to_numpy(float) / sh.v_m_cc_per_g_s.to_numpy(float)**power)[1]


SIGMA_TUNED_SILICA = spread(vm_s / vm_s)[1]
SIGMA_TUNED_CHARCOAL = tuned_spread(sh6)
SIGMA_PRINTED_CHARCOAL_6 = spread(sh6.surface_liquid_m2_per_g_c)[1]
SIGMA_TUNED_POWER = SIGMA_TUNED_CHARCOAL / SIGMA_PRINTED_CHARCOAL_6
print(f"\nCIRCULARITY BOUND, on the {len(sh6)} non-butane isotherms Tables III and V share.")
print(f"  A sigma set tuned to flatten silica gel exactly ({SIGMA_TUNED_SILICA:.2e} % spread)")
print(f"  leaves the charcoal areas spread by {SIGMA_TUNED_CHARCOAL:.4f} %, against "
      f"{SIGMA_PRINTED_CHARCOAL_6:.4f} % for the")
print(f"  printed set -- a factor {SIGMA_TUNED_POWER:.4f}. So the charcoal agreement does carry")
print("  independent information, but the margin is a factor of about two, not an order.")
print("  AND THE TEST HAS NO POWER AT ALL IF BUTANE IS KEPT IN: on all seven shared pairs")
sh7 = shared
SIGMA_TUNED_CHARCOAL_7 = tuned_spread(sh7)
SIGMA_PRINTED_CHARCOAL_7 = spread(sh7.surface_liquid_m2_per_g_c)[1]
print(f"  the tuned set gives {SIGMA_TUNED_CHARCOAL_7:.4f} % "
      f"and the PRINTED set {SIGMA_PRINTED_CHARCOAL_7:.4f} % -- indistinguishable,")
print("  because butane dominates both. The bound below is quoted on the six-row basis only.")
print("  This bounds, and does not remove, the risk that Emmett & Brunauer (1937)'s")
print("  cross-sections were themselves adjusted; that paper is not on disk.")

# %%
fig, ax = plt.subplots(figsize=(7.0, 3.6))
labels = ["silica gel\n(7 gases)", "charcoal\n(7 isotherms)"]
xpos = np.arange(2); w = 0.35
ax.bar(xpos - w/2, [T3_VM_DEV, T5_VM_DEV], w, label="spread in fitted $v_m$ (null baseline)", color="C3")
ax.bar(xpos + w/2, [T3_LIQ_DEV, T5_LIQ_DEV], w, label="spread in derived surface area", color="C0")
for i, (a, b) in enumerate([(T3_VM_DEV, T3_LIQ_DEV), (T5_VM_DEV, T5_LIQ_DEV)]):
    ax.text(i - w/2, a + 1, f"{a:.1f} %", ha="center", fontsize=9)
    ax.text(i + w/2, b + 1, f"{b:.1f} %", ha="center", fontsize=9)
ax.set_xticks(xpos); ax.set_xticklabels(labels)
ax.set(ylabel="maximum deviation from the mean, %", ylim=(0, 56),
       title="One surface from many gases: the test that is not a fit")
ax.legend(fontsize=8); fig.tight_layout(); plt.show()

# %% [markdown]
# **Two further things in the "test" column, both cheap and both real.**
#
# *A surface measured without adsorption at all.* Palmer and Clark measured the
# specific surface of a vitreous silica by the rate at which it dissolves in
# hydrofluoric acid. Brunauer, Emmett and Teller ran an eq. (A) plot on the same
# sample's acetone isotherm and compared. Two cross-sections were available for
# acetone and they bracket the HF result.
#
# *A $v_\mathrm{m}$ predicted at a temperature it was not fitted at.* Fig. 4's
# −195.8 °C curves are **calculated** from the −183 °C isotherms with
# $v_\mathrm{m} \propto d_\mathrm{L}^{2/3}$; Fig. 3 fits the −195.8 °C nitrogen
# isotherm independently. The two $v_\mathrm{m}$ can be compared, and the paper
# never does it.

# %%
PC_LIQUID = P["pc_area_from_liquid"] / P["pc_area_HF"]
PC_ADAM = P["pc_area_from_adam"] / P["pc_area_HF"]
pc_recalc_liquid = (P["pc_vm_acetone"] * 1e-6 / P["pc_sample_mass"] * N_AV
                    * P["pc_sigma_acetone_liquid"] * 1e-16)
pc_recalc_adam = (P["pc_vm_acetone"] * 1e-6 / P["pc_sample_mass"] * N_AV
                  * P["pc_sigma_acetone_adam"] * 1e-16)
pcp = pd.DataFrame([
    ("from the liquid-acetone density, 26.9 A^2", pc_recalc_liquid, P["pc_area_from_liquid"],
     PC_LIQUID, f"about {P['pc_larger']:.0f} % larger"),
    ("from Adam's close-packed films, 20.5 A^2", pc_recalc_adam, P["pc_area_from_adam"],
     PC_ADAM, f"about {P['pc_smaller']:.0f} % smaller"),
], columns=["cross-section used", "area recomputed cm2/g", "area printed cm2/g",
            "BET / HF-dissolution", "the paper says"])
display(pcp.round(4))
print(f"The two printed areas are recovered from v_m = {P['pc_vm_acetone']:g} micromol on "
      f"{P['pc_sample_mass']:g} g to "
      f"{(pc_recalc_liquid/P['pc_area_from_liquid']-1)*100:+.3f} % and "
      f"{(pc_recalc_adam/P['pc_area_from_adam']-1)*100:+.3f} %, using a modern Avogadro number;")
print("the residual is of the size the 1938 value of N_A would explain, and the two RATIOS")
print("against Palmer and Clark are exact to the printed words either way.")

VM_TEMP_DEV = (P["fig4_N2_1958_vm"] / P["fig3_vm"] - 1) * 100
DL_IMPLIED = (P["fig4_N2_1958_vm"] / P["fig4_N2_183_vm"])**1.5
print(f"\nv_m for N2 on catalyst 954 at 77.3 K:")
print(f"  calculated from the 90.1 K isotherm (Fig. 4) : {P['fig4_N2_1958_vm']:g} cc")
print(f"  fitted to the 77.3 K isotherm    (Fig. 3)    : {P['fig3_vm']:g} cc")
print(f"  the calculated value is {VM_TEMP_DEV:+.4f} % of the fitted one.")
print(f"  (the d_L^(2/3) rule with these two numbers implies a liquid-density ratio of "
      f"{DL_IMPLIED:.5f}\n   between 90.16 and 77.36 K -- recovered, not printed, and not checked")
print("   against any density table, none being on disk.)")

# %% [markdown]
# ### 6. The agreement that was mostly guaranteed: point B against $v_\mathrm{m}$
#
# The paper leans on this one hard. Point B is where the isotherm's approximately
# linear middle begins, read by eye off the raw curve, nothing about it is fitted
# , and it is claimed to agree with $v_\mathrm{m}$, "the two seldom differing by as
# much as 10 %".
#
# The claim is true. It is also worth much less than it looks, and the BET equation
# itself says why.
#
# **A scale-free landmark.** "Where the linear portion begins" is not scale-free,
# move the axes and it moves, but the isotherm's *inflection point*, where
# curvature changes sign, is. Any eye-read "beginning of the linear portion" lies
# at or below it. So compute where the inflection is, and how much is adsorbed
# there.
#
# **The result is a ceiling.** Maximising the inflection value over $c$ gives a
# number that does not depend on the adsorbent, the gas, the temperature or the
# fit:
#
# $$\max_{c>1}\;\frac{v(x_\mathrm{infl})}{v_\mathrm{m}} \;=\; \frac{2}{\sqrt3} \;=\; 1.1547\ldots
# \qquad\text{at}\quad c = 27+15\sqrt3, \quad x = 3\sqrt3-5.$$
#
# Derived below symbolically, and confirmed by a root-find that shares none of that
# algebra.
#
# **The guarantee is one-sided, and saying otherwise would be refuted by this
# paper's own Table III.** $2/\sqrt3$ bounds how far the inflection can sit
# *above* $v_\mathrm{m}$; it says nothing about how far *below* $v_\mathrm{m}$ an
# eye-read landmark can fall. So an eye-read knee cannot be **high** by more than
# 15.5 %, it can be arbitrarily low, and on Table III's butane isotherm point B
# is only 0.4828 of $v_\mathrm{m}$, printed two cells below. Butane is dropped from
# this page's eighteen-pair average because *the authors flag it empirically*, not
# because the ceiling excludes it. What the ceiling deflates is therefore the
# upper half of the agreement, which is the half the authors are claiming, since
# they report the two "seldom differing by as much as 10 %" and most of that was
# promised by the shape of their own equation before any data arrived.

# %%
# --- closed form: maximise v/v_m over c ON the inflection locus (Lagrange)
f_sym = c * x / ((1 - x) * (1 - x + c * x))
g_sym = sp.expand(sp.numer(sp.together(sp.diff(f_sym, x, 2))) / (-2 * c))     # inflection locus
lagr = sp.numer(sp.together(sp.simplify(sp.diff(f_sym, x) * sp.diff(g_sym, c)
                                        - sp.diff(f_sym, c) * sp.diff(g_sym, x))))
sol = [s for s in sp.solve([lagr, g_sym], [x, c], dict=True)]
sol = [s for s in sol if sp.im(sp.N(s[x])) == 0 and sp.N(s[x]) > 0][0]
X_STAR_SYM = sp.radsimp(sp.simplify(sol[x]))
C_STAR_SYM = sp.nsimplify(sp.simplify(sol[c]), [sp.sqrt(3)])
CEIL_SYM = sp.radsimp(sp.simplify(f_sym.subs(sol)))
print(f"symbolic:  x* = {X_STAR_SYM} = {float(sp.N(X_STAR_SYM)):.12f}")
print(f"           c* = {C_STAR_SYM} = {float(sp.N(C_STAR_SYM)):.10f}")
print(f"      ceiling = {CEIL_SYM} = {float(sp.N(CEIL_SYM, 20)):.15f}")

# --- SECOND, INDEPENDENT ROUTE: root-find the inflection, then root-find its maximum in c
# on a central-difference derivative with relative step H_FD (model.c_star_root)
H_FD = 1e-5                          # the page's step: truncation and round-off balance here
H_LADDER = (1e-2, 1e-3, 1e-4)        # steps of the truncation study, ratio 10
C_STAR_EXACT = float(sp.N(C_STAR_SYM, 20))
C_STAR_NUM = c_star_root(H_FD)
CEIL_NUM = inflection_ratio(C_STAR_NUM)
X_STAR_NUM = x_inflection(C_STAR_NUM)
CEIL_EXACT = 2.0 / np.sqrt(3.0)
CEIL_TWO_ROUTE = abs(CEIL_NUM / CEIL_EXACT - 1.0)
C_STAR_TWO_ROUTE = abs(C_STAR_NUM / C_STAR_EXACT - 1.0)
print(f"\nroot-found: x* = {X_STAR_NUM:.12f}   c* = {C_STAR_NUM:.10f}   ceiling = {CEIL_NUM:.15f}")
print(f"two routes on the ceiling agree below 1e-12 relative (here {CEIL_TWO_ROUTE:.1e}; the routes share no "
      f"algebra: one is a Lagrange")
print(" condition solved symbolically, the other a root-find on a differentiated root-find).")
print(f"c* to the closed form at the step h = {H_FD:g}: {C_STAR_TWO_ROUTE:.1e} relative.")


def c_star_study(h_ladder=H_LADDER, central=True):
    """Observed order of the c* error in the derivative step, from successive differences."""
    return observed_orders(lambda nh: c_star_root(1.0 / nh, central), [round(1.0 / h) for h in h_ladder])


FD = c_star_study()
C_STAR_FD_ORDER = float(FD["orders"][-1])
C_STAR_LADDER_OVER_FLOOR = abs(FD["values"][-2] - FD["values"][-1]) / abs(C_STAR_NUM - C_STAR_EXACT)
print("\nWHAT THE STEP COSTS. c* against the closed form on a ladder of steps:")
for h, cv in zip(H_LADDER, FD["values"]):
    print(f"   h = {h:g}: c* = {cv:.10f}, relative gap {abs(cv / C_STAR_EXACT - 1):.3e}, "
          f"gap / h^2 = {abs(cv / C_STAR_EXACT - 1) / h**2:.4f}")
print(f"observed order from successive differences {C_STAR_FD_ORDER:.4f} (central difference, theory 2);")
print(f"the finest successive difference is {C_STAR_LADDER_OVER_FLOOR:.0f} times the gap left at "
      f"h = {H_FD:g}, so the ladder is clear of that floor.")
print(f"At h = {H_FD:g} the h^2 term predicts {abs(FD['values'][-1] / C_STAR_EXACT - 1) / H_LADDER[-1]**2 * H_FD**2:.1e};"
      " the rest of the gap is round-off in the difference quotient, which grows as h shrinks.")
print("The ceiling is a quadratic maximum, so it carries the SQUARE of that error: its two")
print("routes agree to round-off whatever step the ladder uses.")

# %%
# c reconstructed from E1 - EL, footnote 16, prefactor set to 1 -- A RECONSTRUCTION
c_t1 = np.exp(t1.E1_minus_EL_cal_per_mol.to_numpy(float) / (R_CAL * T_TABLE1))
c_t3 = np.exp(t3.E1_minus_EL_cal_per_mol.to_numpy(float)
              / (R_CAL * (t3.temperature_C.to_numpy(float) + KELVIN)))

obs = np.r_[(t1.point_B_cc_per_g / t1.v_m_cc_per_g).to_numpy(float),
            (t3.point_B_cc_per_g / t3.v_m_cc_per_g).to_numpy(float)]
pred = np.array([inflection_ratio(cv) for cv in np.r_[c_t1, c_t3]])
lab = np.r_[t1.substance.to_numpy(),
            (t3.gas + " at " + t3.temperature_C.astype(str) + " C, silica gel").to_numpy()]
cvals = np.r_[c_t1, c_t3]
keep = ~np.char.startswith(lab.astype(str), "C4H10")     # the row the authors flag themselves

pb = pd.DataFrame({"isotherm": lab, "c (reconstructed)": cvals,
                   "point B / v_m observed": obs,
                   "inflection / v_m predicted": pred,
                   "predicted / observed": pred / obs})
display(pb.round(4))



def pointB_stats(t1_=t1, t3_=t3, c_scale=1.0, with_butane=False):
    """Point B against v_m and against the inflection, from the two tables as given.

    c_scale multiplies every reconstructed c; with_butane puts back the row the
    authors exclude."""
    ob = np.r_[(t1_.point_B_cc_per_g / t1_.v_m_cc_per_g).to_numpy(float),
               (t3_.point_B_cc_per_g / t3_.v_m_cc_per_g).to_numpy(float)]
    pr = np.array([inflection_ratio(cv) for cv in cvals * c_scale])
    kp = np.ones_like(keep) if with_butane else keep
    dev1 = np.abs(t1_.v_m_cc_per_g - t1_.point_B_cc_per_g) / t1_.point_B_cc_per_g * 100
    # the count in exact decimal arithmetic on the printed cells: one row sits at exactly 10 %
    tol = Fraction(str(P["t1_vm_vs_B_tolerance"])) / 100
    dev_exact = [abs(Fraction(str(v)) - Fraction(str(b))) / Fraction(str(b))
                 for v, b in zip(t1_.v_m_cc_per_g, t1_.point_B_cc_per_g)]
    n_over_tol = sum(d > tol for d in dev_exact)
    at_tol = [s_ for s_, d in zip(t1_.substance, dev_exact) if d == tol]
    ratio = float((pr[kp] / ob[kp]).mean())
    return {"pointB_over_vm_mean": float(ob[kp].mean()), "pointB_over_vm_sd": float(ob[kp].std(ddof=1)),
            "inflection_over_pointB_mean": ratio, "inflection_above_pointB_pct": (ratio - 1) * 100,
            "pointB_below_inflection_pct": (1 - 1 / ratio) * 100,
            "t1_worst_vm_vs_B_pct": float(dev1.max()), "t1_rows_vm_vs_B_over_tolerance": n_over_tol,
            "t1_rows_at_tolerance": at_tol, "t1_worst_row": str(t1_.substance.iloc[int(np.argmax(dev1.to_numpy()))]),
            "t3_butane_vm_over_B": float((t3_.v_m_cc_per_g / t3_.point_B_cc_per_g)[t3_.gas == "C4H10"].item()),
            "lowest_pointB_over_vm": float(ob.min())}


PBS = pointB_stats()
PB_MEAN, PB_SD = PBS["pointB_over_vm_mean"], PBS["pointB_over_vm_sd"]
PRED_MEAN = float(pred[keep].mean())
PRED_OVER_OBS = PBS["inflection_over_pointB_mean"]
T1_WORST_VM_B = PBS["t1_worst_vm_vs_B_pct"]
T1_ROWS_OVER_10 = PBS["t1_rows_vm_vs_B_over_tolerance"]
T3_BUTANE_RATIO = PBS["t3_butane_vm_over_B"]
print(f"18 isotherms (the butane row excluded, as the authors exclude it):")
print(f"  observed point B / v_m : mean {PB_MEAN:.6f}, sd {PB_SD:.6f}")
print(f"  predicted inflection   : mean {PRED_MEAN:.6f}, range over all 19 rows "
      f"{pred.min():.6f} to {pred.max():.6f}")
print(f"  every predicted value is under the ceiling {CEIL_EXACT:.6f}: "
      f"{bool((pred < CEIL_EXACT).all())}")
INFL_ABOVE_B = PBS["inflection_above_pointB_pct"]        # mean(infl/B) - 1
PB_BELOW_INFL = PBS["pointB_below_inflection_pct"]       # the same statement the other way up
PB_BELOW_INFL_ELEM = float((1 - obs[keep] / pred[keep]).mean()) * 100
print(f"  DIRECTION, stated once and used everywhere: mean(inflection / point B) = "
      f"{PRED_OVER_OBS:.6f},")
print(f"  i.e. the inflection sits {INFL_ABOVE_B:.3f} % ABOVE point B, equivalently point B sits")
print(f"  {PB_BELOW_INFL:.3f} % BELOW the inflection ({PB_BELOW_INFL_ELEM:.3f} % if the ratio is")
print("  averaged elementwise as 1 - B/inflection). The page quotes the two-thirds-place forms.")
print(f"\nTHE CEILING IS ONE-SIDED. Lowest point B / v_m over all 19 rows: "
      f"{obs.min():.4f} ({lab[int(obs.argmin())]}),")
print(f"  i.e. {(1-obs.min())*100:.2f} % BELOW v_m -- so 'inside 15.5 % of v_m' is false as a "
      f"two-sided claim.")
print("  Nothing bounds a landmark from below; the ceiling bounds it only from above.")
print(f"\nTable I, the paper's own claim ('seldom differing by as much as "
      f"{P['t1_vm_vs_B_tolerance']:.0f} %'):")
_tol = P["t1_vm_vs_B_tolerance"]
print(f"  worst |v_m - B|/B = {T1_WORST_VM_B:.4f} % ({PBS['t1_worst_row']}, where both entries")
print(f"  carry only two significant figures); rows above {_tol:g} % in exact arithmetic on the printed")
print(f"  cells: {T1_ROWS_OVER_10} of {len(t1)}; at exactly {_tol:g} %: {', '.join(PBS['t1_rows_at_tolerance'])}")
print("  (the paper's 'as much as' includes equality). The claim holds.")
print(f"Table III butane, which the authors flag: v_m / B = {T3_BUTANE_RATIO:.4f} "
      f"against their 'twice as large'.")

# %%
C_SSHAPE_MIN = 2.0          # the cubic's value at x = 0 is (2 - c): no inflection below c = 2
print(f"eq. (26) has an inflection -- is S-shaped at all -- only for c > {C_SSHAPE_MIN:g}.")
print(f"  check: inflection_ratio(1.9) = {inflection_ratio(1.9)}, "
      f"inflection_ratio(2.1) = {inflection_ratio(2.1):.6f}")
print("  The paper puts it qualitatively: 'The constant c, as a rule, will be large")
print("  compared to unity, and therefore the isotherm will consist of two regions.'")
print("  Every c on this page, reconstructed or printed, is far above 2.\n")

cs = np.logspace(np.log10(2.05), 4, 300)
rat = np.array([inflection_ratio(cv) for cv in cs])
fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
ax[0].semilogx(cs, rat, lw=1.6, label="BET inflection point")
ax[0].axhline(CEIL_EXACT, color="C3", ls="--", lw=1.2,
              label=r"ceiling $2/\sqrt{3}$ = %.4f" % CEIL_EXACT)
ax[0].plot([C_STAR_NUM], [CEIL_NUM], "o", color="C3", ms=6)
ax[0].scatter(cvals[keep], obs[keep], s=18, color="k", zorder=5,
              label="point B / $v_m$, Tables I and III")
ax[0].axhline(1.0, color="0.6", lw=0.8)
ax[0].set(xlabel="c", ylabel="$v/v_m$ at the landmark", ylim=(0.8, 1.25),
          title="The inflection point is capped; point B sits below it")
ax[0].legend(fontsize=7.5, loc="lower right")

xs = np.linspace(0.001, 0.65, 500)
for cv, col in [(25, "C0"), (84, "C1"), (350, "C2")]:
    ax[1].plot(xs, v_over_vm_A(xs, cv), col, lw=1.4, label=f"c = {cv}")
    xi = x_inflection(cv)
    ax[1].plot([xi], [v_over_vm_A(xi, cv)], "o", color=col, ms=5)
ax[1].axhline(CEIL_EXACT, color="C3", ls="--", lw=1.2)
ax[1].axhline(1.0, color="0.6", lw=0.8)
ax[1].set(xlabel="$p/p_0$", ylabel="$v/v_m$", ylim=(0, 1.8),
          title="Inflection points (dots) of eq. (A) at three c")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()
print("The c on the left axis are RECONSTRUCTED by inverting footnote 16 -- which is how")
print("the E1 - EL column was made in the first place (journal p. 313), so they are the")
C_FLAT = (25.0, 350.0)          # the decade of c the reconstructed values span
_flat = [inflection_ratio(cv) for cv in np.logspace(*np.log10(C_FLAT), 200)]
print(f"authors' own fitted c to the {0.5 / (R_CAL * T_TABLE1) * 100:.2f} % the integer rounding of E "
      "allows. The ceiling")
print(f"above them holds for every c whatever, and the curve varies by {(1 - min(_flat) / max(_flat)) * 100:.1f} % "
      "of its maximum over")
print(f"{C_FLAT[0]:g} < c < {C_FLAT[1]:g}, so this panel is insensitive to c either way --")
print("which is a property of the ceiling, NOT evidence that anything else on the page is.")


# %% [markdown]
# ### 7. The single-point rule, and the exact $c$ it needs
#
# Journal page 315 offers a shortcut: join the 760 mm nitrogen point to the origin
# on an eq. (A) plot and the slope is $1/v_\mathrm{m}$ "with an error of no more
# than 5 % on all solid adsorbents with which we have worked". The paper adds that
# this amounts to $v_\mathrm{m} = v(1-p/p_0)$, and that the 760 mm point sits near
# $p/p_0 = 1/3$.
#
# That statement has an exact content. On the BET isotherm,
#
# $$\frac{v\,(1-x)}{v_\mathrm{m}} \;=\; \frac{c\,x}{1+(c-1)x},$$
#
# which is **always below 1**, the shortcut can only ever *underestimate*
# $v_\mathrm{m}$, a one-sidedness the paper does not mention. Requiring the
# underestimate to be no worse than a fraction $\varepsilon$ at relative pressure
# $x$ gives, in closed form,
#
# $$c \;\ge\; \frac{(1-\varepsilon)(1-x)}{\varepsilon\,x},$$
#
# which at $\varepsilon = 0.05$ and $x = 1/3$ is exactly $c \ge 38$.

# %%
if abs(P["single_point_x"] - 1.0 / 3.0) > 1e-9:
    raise ValueError("the CSV should carry the printed 1/3")
X_760 = 1.0 / 3.0                    # the paper prints the FRACTION 1/3, not a decimal
EPS = P["single_point_error"] / 100.0
C_LO, C_HI = 2.0, 1e9                # c scanned (log grid): from the S-shape limit upwards


def c_threshold(eps=EPS, xv=X_760):
    """Closed form of the smallest c for which the single-point rule is within eps at xv."""
    return (1 - eps) * (1 - xv) / (eps * xv)


def single_point(E_col=t1.E1_minus_EL_cal_per_mol, R=R_CAL, xv=X_760, c_scale=1.0):
    """Worst single-point underestimate %, smallest c and rows past the promised 5 %, Table I."""
    cc = np.exp(np.asarray(E_col, float) / (R * T_TABLE1)) * c_scale
    err = (1.0 - shortcut_ratio(cc, xv)) * 100.0
    return {"worst": float(err.max()), "c_min": float(cc.min()), "err": err,
            "rows_past": int((err > P["single_point_error"]).sum())}


C_THRESH_CLOSED = c_threshold()
C_THRESH_ROOT = unique_root(lambda cv: shortcut_ratio(cv, X_760) - (1 - EPS), C_LO, C_HI,
                            log=True, xtol=1e-12)
C_THRESH_TWO_ROUTE = abs(C_THRESH_ROOT / C_THRESH_CLOSED - 1)
print(f"closed form  c >= {C_THRESH_CLOSED:.10f}")
print(f"root-found   c  = {C_THRESH_ROOT:.10f}   (relative deviation {C_THRESH_TWO_ROUTE:.1e})")

SPT = single_point()
short_err = SPT["err"]
SHORT_WORST = SPT["worst"]
C_T1_MIN = SPT["c_min"]
sp_tbl = t1[["substance"]].assign(**{"c (reconstructed)": c_t1,
                                     "single-point underestimate %": short_err})
display(sp_tbl.round(4))
print(f"smallest reconstructed c in Table I: {C_T1_MIN:.4f}, comfortably above the "
      f"threshold {C_THRESH_CLOSED:.0f}.")
print(f"worst underestimate over the twelve: {SHORT_WORST:.4f} % against the promised "
      f"{P['single_point_error']:.0f} %. The claim holds, and holds one-sidedly.")
n_below = int((c_t3 < C_THRESH_CLOSED).sum())
print("\nAND IT IS EXACTLY AS NARROW AS THE AUTHORS MADE IT. They restrict the rule to")
print(f"nitrogen at -183 C. {n_below} of the seven Table III isotherms have a reconstructed c")
print(f"BELOW the threshold {C_THRESH_CLOSED:.0f} --")
for g, T, cv in zip(t3.gas, t3.temperature_C, c_t3):
    if cv < C_THRESH_CLOSED:
        print(f"   {g:6s} at {T:7.1f} C : c = {cv:7.3f}  ->  the rule would be off by "
              f"{(1-shortcut_ratio(cv, X_760))*100:.2f} %")
print("-- so the sentence is scoped to precisely the gas for which it survives.")

# %% [markdown]
# ### 8. Footnote 16's prefactor is not measurable here, but the paper's two routes to one $c$ differ by 2.24
#
# It is tempting to read Fig. 3's fitted $c = 156.7$ against Fig. 4's
# $E_1-E_\mathrm{L} = 900$ cal, same gas, same sample, same temperature, as a
# measurement of $a_1b_2/b_1a_2$. **It is not one, and the paper says why in three
# places.**
#
# 1. Journal page 313: *"From $c$ one can obtain an approximate value for
#    $E_1-E_\mathrm{L}$."* The energies are made **from** $c$, through footnote 16,
#    with the prefactor already set to 1. There is no independently determined
#    $E_1-E_\mathrm{L}$ anywhere in the paper against which a prefactor could be
#    measured, and inverting footnote 16 on Table I returns the authors' own
#    fitted $c$, to $\pm0.28$ %, the integer rounding of the printed column.
# 2. Fig. 4's energies are Tables I and II's, rounded: nitrogen on catalyst 954 is
#    894 cal in Table I and appears as 900 in Fig. 4; argon on the same sample is
#    704 in Table II and 700 in Fig. 4.
# 3. Journal page 316 states that Fig. 4's −195.8 °C curves *"have been calculated
#    from the −183 ° isotherms"*, on the explicit assumption that
#    $E_1-E_\mathrm{L}$ *"changes only slightly with temperature"*. So the 900 in
#    the −195.8 °C block is the **90.1 K** value carried over, not a 77.3 K fit,
#    both temperature blocks print the same number for each gas, which is the
#    fingerprint of the carry-over.
#
# **What the comparison does measure is an internal inconsistency**, and it is the
# one that matters for every $c$ on this page. The 77.3 K nitrogen isotherm on
# catalyst 954 appears twice: fitted directly in Fig. 3, and *calculated* in
# Fig. 4. The two carry $c$ values a factor 2.24 apart. Below, and then applied as a sensitivity to sections 6 and 7, where
# **section 7 does not survive it**.

# %%
C_FIG4_AT_773 = float(np.exp(P["fig4_N2_1958_E"] / (R_CAL * P["fig3_T"])))
C_TWO_ROUTE = P["fig3_c"] / C_FIG4_AT_773          # 1 / 2.24
E_from_c = R_CAL * P["fig3_T"] * np.log(P["fig3_c"])
E_T1_954 = float(t1.loc[t1.substance == "Fe-Al2O3 catalyst 954",
                        "E1_minus_EL_cal_per_mol"].item())
C_T1_954 = float(np.exp(E_T1_954 / (R_CAL * T_TABLE1)))
print("THE SAME ISOTHERM, N2 on catalyst 954 at 77.3 K, from the paper's two routes:")
print(f"  Fig. 3, FITTED directly                       : c = {P['fig3_c']:g}")
print(f"  Fig. 4, CALCULATED with E1-EL = {P['fig4_N2_1958_E']:g} at "
      f"{P['fig3_T']:g} K : c = {C_FIG4_AT_773:.4f}")
print(f"  ratio {C_TWO_ROUTE:.6f}, i.e. a factor {1/C_TWO_ROUTE:.4f} apart.")
print(f"Read as energies: the fitted c implies E1 - EL = {E_from_c:.4f} cal at 77.3 K, against")
print(f"  Table I's {E_T1_954:.0f} cal for the same adsorbent at {T_TABLE1:g} K "
      f"({(E_from_c/E_T1_954-1)*100:+.3f} %) and against")
print(f"  the {P['fig4_N2_1958_E']:g} cal Fig. 4 reuses ({(E_from_c/P['fig4_N2_1958_E']-1)*100:+.3f} %).")
print("This is NOT footnote 16's prefactor. E1 - EL was obtained FROM c through footnote 16")
print("with the prefactor set to 1 (journal p. 313), so inverting it returns the authors'")
print(f"  own fitted c: Table I's 954 row, {E_T1_954:.0f} cal at {T_TABLE1:g} K, gives "
      f"c = {C_T1_954:.4f},")
print(f"  and the integer rounding of that column costs only "
      f"{0.5/(R_CAL*T_TABLE1)*100:.4f} % in c.")
print(f"The {1 / C_TWO_ROUTE:.2f} is the paper's own inconsistency between a direct fit and a temperature")
print("extrapolation, and it is the honest size of the doubt on any c on this page.\n")

# --- APPLY IT. Which sections survive a c wrong by that factor, and which do not?
SPT_SCALED = single_point(c_scale=C_TWO_ROUTE)
SHORT_WORST_SCALED = SPT_SCALED["worst"]
C_T1_MIN_SCALED = SPT_SCALED["c_min"]
T1_ROWS_BREACHING_SCALED = SPT_SCALED["rows_past"]
PRED_OVER_OBS_SCALED = pointB_stats(c_scale=C_TWO_ROUTE)["inflection_over_pointB_mean"]


def ceiling_with_c_scale(scale, h=H_FD):
    """The inflection ceiling maximised over c when every c is multiplied by scale."""
    f = lambda cv: inflection_ratio(cv * scale)
    cs = unique_root(lambda cv: (f(cv * (1 + h)) - f(cv * (1 - h))) / (2 * h * cv),
                     20.0 / scale, 200.0 / scale, log=True, xtol=1e-12)
    return f(cs)


CEIL_SCALED = ceiling_with_c_scale(C_TWO_ROUTE)
sens = pd.DataFrame([
    ("section 6, inflection ceiling (max over ALL c)", CEIL_NUM, CEIL_SCALED),
    ("section 6, mean inflection / point B", PRED_OVER_OBS, PRED_OVER_OBS_SCALED),
    ("section 7, smallest Table I c (threshold 38)", C_T1_MIN, C_T1_MIN_SCALED),
    ("section 7, worst single-point underestimate %", SHORT_WORST, SHORT_WORST_SCALED),
], columns=["quantity", "as reconstructed", f"with c x {C_TWO_ROUTE:.4f}"])
display(sens.round(6))
print(f"SECTION 6 SURVIVES, and trivially: the ceiling is a MAXIMUM OVER ALL c, so no")
print("  rescaling of c can move it. That is a structural insensitivity, not evidence.")
print(f"SECTION 7 DOES NOT. At c x {C_TWO_ROUTE:.4f} the smallest Table I c falls to "
      f"{C_T1_MIN_SCALED:.4f},")
print(f"  a ratio {C_T1_MIN_SCALED/C_THRESH_CLOSED:.4f} of the threshold {C_THRESH_CLOSED:.0f} "
      f"-- BELOW it -- and")
print(f"  {T1_ROWS_BREACHING_SCALED} of 12 rows breach the paper's promised "
      f"{P['single_point_error']:.0f} %, the worst by {SHORT_WORST_SCALED:.4f} %.")
print(f"So the single-point rule holds on the paper's own numbers ({SHORT_WORST:.2f} % against "
      f"{P['single_point_error']:.0f} %) but")
print("is NOT robust to the size of doubt the paper itself exhibits. Both are stated; the")
print("page claims only the first. Nothing here reports a reconstructed c as a measurement.")


# %% [markdown]
# ### 9. Printed statements that do not survive
#
# Five, ranked, plus two counted elsewhere on the page: the cross-table $\sigma$
# disagreement that exceeds what the printed rounding can produce (§5) and the
# factor-2.24 clash between the paper's two routes to one $c$ (§8). Each is proved
# from the paper's own numbers and none is repaired.
#
# #### 9a. "840 ± 70 cal. […] for nitrogen on all twelve adsorbents"
#
# Verbatim on the crop, journal page 315, and the ellipsis spans a **sentence
# boundary**, which the page shows rather than hides:
#
# > For nitrogen $E_1-E_\mathrm{L}$ is uniformly 840 ± 70 cal. Since
# > $E_\mathrm{L}$ is about 1330 cal., $E_1$ is therefore 2170 ± 70 cal. for
# > nitrogen on all twelve adsorbents.
#
# The first sentence is about "the last column of Table I", which has twelve rows,
# so the substance is unchanged, but the page is publicly reporting an authors'
# error and the quotation is shown at full length for that reason.
#
# The paper never defines its $\pm$. **So pin the convention first, on the
# statement that is not in dispute.** Journal page 315 quotes four bands for the
# two substances of Table II, 840 ± 50, 650 ± 55, 1460 ± 120, 1900 ± 30, and
# Table II prints the eight numbers behind them. Decode there, then apply to
# Table I.

# %%
T2_GASES = [("N2", "t2_N2_centre", "t2_N2_half"), ("A", "t2_A_centre", "t2_A_half"),
            ("CO2", "t2_CO2_centre", "t2_CO2_half"), ("C4H10", "t2_C4H10_centre", "t2_C4H10_half")]


def t2_bands(t2_=t2):
    """Midrange and half-range of each Table II gas against the printed band."""
    out = []
    for g, kc, kh in T2_GASES:
        v = t2_.loc[t2_.gas == g, "E1_minus_EL_cal_per_mol"].to_numpy(float)
        m, h = band(v)
        out.append((g, v.min(), v.max(), m, h, P[kc], P[kh], m - P[kc], h - P[kh]))
    return pd.DataFrame(out, columns=["gas", "lower", "upper", "midrange", "half-range",
                                      "printed centre", "printed half",
                                      "centre error", "half error"])


conv = t2_bands()
rows = list(conv.itertuples(index=False, name=None))
display(conv)
T2_CENTRE_DEV = float(conv["centre error"].abs().max())
T2_HALF_DEV = float(conv["half error"].abs().max())
T2_CO2_HALF_RANGE = float(conv.loc[conv.gas == "CO2", "half-range"].item())

CENTRE_HITS = sum(1 for r in rows if to10(r[3]) == r[5])
HALF_MULT5 = sum(1 for r in rows if abs(r[6] % 5.0) < 1e-9 and abs(r[6] - r[4]) <= 2.5 + 1e-9)
HALF_HITS_10 = sum(1 for r in rows if to10(r[4]) == r[6])
A_HALF_RANGE = float(conv.loc[conv.gas == "A", "half-range"].item())
print("The paper's 'X +- Y' is the MIDRANGE and the HALF-RANGE -- but the two halves are")
print("NOT rounded the same way, and argon is what shows it.")
print(f"  centre, nearest 10 cal   : {CENTRE_HITS} of 4 exact hits "
      f"(worst residual {T2_CENTRE_DEV:.1f} cal)")
print(f"  half-width, a multiple of 5 within 2.5 cal of the half-range: {HALF_MULT5} of 4")
print(f"    -- argon's half-range is {A_HALF_RANGE:.1f} and is PRINTED as {P['t2_A_half']:.0f}, "
      f"which is not a multiple of 10;")
print(f"    -- CO2's is {float(conv.loc[conv.gas=='CO2','half-range'].item()):.1f}, an exact tie, "
      f"printed as {P['t2_CO2_half']:.0f} (resolved downward).")
print(f"  a nearest-10 rule for the half-width would hit only {HALF_HITS_10} of 4. Worst")
print(f"    half-width residual on the correct rule: {T2_HALF_DEV:.1f} cal.")
print("Nothing was free to make that come out. The Table I conclusion below is unchanged")
print("under EITHER half-width rule, which is checked there rather than assumed.")

# %%
E_T1 = t1.E1_minus_EL_cal_per_mol.to_numpy(float)
MID12, HALF12 = band(E_T1)
keep11 = t1.substance != "Cr2O3 gel"
MID11, HALF11 = band(E_T1[keep11.to_numpy()])
outside = np.abs(E_T1 - P["t1_N2_band_centre"]) > P["t1_N2_band_half"]
T1_ROWS_OUTSIDE = int(outside.sum())

res = pd.DataFrame([
    ("all twelve rows, as the sentence says", MID12, HALF12,
     MID12 - P["t1_N2_band_centre"], HALF12 - P["t1_N2_band_half"]),
    ("eleven rows, Cr2O3 gel (738) omitted", MID11, HALF11,
     MID11 - P["t1_N2_band_centre"], HALF11 - P["t1_N2_band_half"]),
], columns=["basis", "midrange", "half-range", "centre error vs printed 840",
            "half error vs printed 70"])
display(res)
print(f"Rows lying OUTSIDE the printed band 840 +- 70: {T1_ROWS_OUTSIDE} of 12 --")
for s, e in zip(t1.substance[outside], E_T1[outside]):
    print(f"   {s:30s} {e:.0f} cal")

RULES = {"half-width to the nearest 5 (the rule Table II pins)": to5,
         "half-width to the nearest 10 (the coarser alternative)": to10}


def omission_hits(E):
    """For each half-width rule, the single-row omissions whose band rounds to the printed one."""
    hits_ = {}
    for name, hrule in RULES.items():
        h_ = []
        for dr in itertools.combinations(range(len(E)), 1):
            m, h = band(np.delete(E, list(dr)))
            if (to10(m), hrule(h)) == (P["t1_N2_band_centre"], P["t1_N2_band_half"]):
                h_.append(tuple(t1.substance[i] for i in dr))
        hits_[name] = h_
    return hits_


def omissions_under_both(E):
    """Number of single-row omissions that reproduce the printed band under BOTH rules."""
    h = omission_hits(E)
    return len(set.intersection(*(set(v) for v in h.values())))


hits = omission_hits(E_T1)
print(f"\nRounded on the convention Table II pins:")
print(f"   twelve rows -> {to10(MID12):.0f} +- {to5(HALF12):.0f}   "
      f"(+- {to10(HALF12):.0f} on the coarser half-width rule)")
print(f"   eleven rows -> {to10(MID11):.0f} +- {to5(HALF11):.0f}   "
      f"(+- {to10(HALF11):.0f} on the coarser rule) = the printed "
      f"{P['t1_N2_band_centre']:.0f} +- {P['t1_N2_band_half']:.0f}")
print("Of the twelve possible single-row omissions, the number that reproduce the printed")
print("band -- and the identity of the row -- is the same under both rules:")
for name, h_ in hits.items():
    print(f"   {name}: {len(h_)} -> {h_[0][0] if h_ else 'none'}")
print(f"The eleven-row pair is {abs(MID11-P['t1_N2_band_centre']):.1f} and "
      f"{abs(HALF11-P['t1_N2_band_half']):.1f} cal from the printed one --")
print(f"well inside the {T2_CENTRE_DEV:.0f} cal rounding Table II already showed.")
T1_OMISSIONS_BOTH = omissions_under_both(E_T1)
print(f"omissions reproducing the printed band under both rules: {T1_OMISSIONS_BOTH} "
      f"({', '.join(sorted({h_[0][0] for h_ in hits.values() if h_}))})")
print("E1 = 2170 +- 70 inherits the same eleven-row basis, being 840 + 1330.")
print("\nREPORTED, NOT REPAIRED. The tables stay as printed. What the page claims is the")
print("narrow claim -- the number belongs to eleven of the twelve rows, and the sentence")
print("says twelve -- not that anybody knows which of the two the authors intended.")

# %% [markdown]
# #### 9b. Footnote 10's lattice sum is mis-grouped
#
# `(8.357 - 8.823/r^3)` cannot be what was meant: it is dimensionally inhomogeneous
# and does not give the 0.466 used two lines below, while $(8.357-8.823)/r^3$ gives
# it exactly. Shown in section 1.
#
# #### 9c. Table IV's $E_1$ column does not close on one row
#
# $E_1$ should be $(E_1-E_\mathrm{L}) + E_\mathrm{L}$. Two rows close to 1 cal or
# better; the third is 5 cal out, on a value that lands exactly on a rounding
# half-way point.
#
# #### 9d. Catalyst 954's $v_\mathrm{m}$ appears twice and disagrees by 15.6 %
#
# Table I gives 2.86 cm³/g for nitrogen at 90.1 K on "Fe–Al₂O₃ catalyst 954";
# Fig. 4 gives 124.7 cm³ for nitrogen at −183 °C on catalyst 954, on a sample it
# states as 50.4 g. The paper prints **no sample mass for Table I's row**, so the
# page cannot decide whether these are two samples or one number is wrong. Both
# readings are printed and neither is adopted.
#
# #### 9e. Table V's two argon rows repeat a $v_\mathrm{m}$ across a 13 K interval
#
# Nitrogen, the only other gas Table V measures at both −195.8 °C and −183 °C, has
# two different $v_\mathrm{m}$ there. Argon has **one** number, 215.5, at both,
# and its solid-packing surface repeats with it, as it must. (Nitrogen's
# *liquid*-packing surface also repeats, at 795, but from two different
# $v_\mathrm{m}$, so that one is a cross-section changing with temperature and not
# an oddity; the cell prints both so the distinction is visible rather than
# asserted.) The CSV's sidecar records the argon repeat, and it is printed here so
# the claim lives on the page and not only in the sidecar. Nothing on this page
# reads the −195.8 °C argon row, §5 merges Tables III and V on gas *and*
# temperature, and Table III has no argon at −195.8 °C, so nothing is affected,
# which is worth stating rather than leaving to be inferred.

# %%
VM4 = t4.v_m_cc_per_g.to_numpy(float); DL4 = t4.d_L_g_per_cc.to_numpy(float)
T4_VM_RATIO = float(VM4[-1] / VM4[0])          # the -80 C row (last) over the 40 C row (first)
T4_DL_RATIO = float((DL4[-1] / DL4[0])**(2.0 / 3.0))
closure = (t4.E1_cal_per_mol - t4.E1_minus_EL_cal_per_mol - t4.EL_cal_per_mol).dropna()
T4_CLOSURE_MAX = float(closure.abs().max())
print(f"Table IV, the two printed ratios (journal page 317):")
print(f"  v_m(-80 C)/v_m(40 C)          = {T4_VM_RATIO:.6f}   printed {P['t4_vm_ratio']:g}")
print(f"  (d_L(-80)/d_L(40))^(2/3)      = {T4_DL_RATIO:.6f}   printed {P['t4_dL_ratio_power']:g}")
print(f"  both exact to the printed precision; the d_L^(2/3) rule is off by a factor "
      f"{T4_VM_RATIO/T4_DL_RATIO:.4f} here,")
print("  which is the point the authors are making with this table.")
print(f"\nTable IV, E1 against (E1 - EL) + EL: residuals {list(closure.astype(int))} cal, "
      f"worst {T4_CLOSURE_MAX:.0f}.")
print(f"  the -5 is on 1705 + 5840 = 7545 printed as 7540: a round-half ambiguity, not a slip.")

vm954 = float(t1.loc[t1.substance == "Fe-Al2O3 catalyst 954", "v_m_cc_per_g"].item())
T1_954_DEV = (vm954 * P["fig3_mass"] / P["fig4_N2_183_vm"] - 1) * 100
print(f"\nCatalyst 954, nitrogen at -183 C / 90.1 K, from two places in the same paper:")
print(f"  Table I : {vm954:g} cc/g  x  {P['fig3_mass']:g} g  = {vm954*P['fig3_mass']:.3f} cc")
print(f"  Fig. 4  : {P['fig4_N2_183_vm']:g} cc               = "
      f"{P['fig4_N2_183_vm']/P['fig3_mass']:.5f} cc/g")
print(f"  Table I is {T1_954_DEV:+.4f} % HIGHER than Fig. 4 on the same nominal sample.")
print("  No mass is printed for Table I's row, so the page cannot say which is wrong.")

both_T = t5.groupby("gas").filter(lambda g: len(g) > 1)
rep = (both_T.groupby("gas")[["v_m_cc_per_g", "surface_solid_m2_per_g",
                              "surface_liquid_m2_per_g"]].nunique() == 1)
print(f"\nTable V, gases measured at BOTH -195.8 and -183 C -- columns that repeat:")
for g in rep.index:
    same = [c for c in rep.columns if bool(rep.loc[g, c])]
    vals = both_T.loc[both_T.gas == g, "v_m_cc_per_g"].tolist()
    print(f"   {g:4s} v_m {vals} -> repeated columns: {same if same else 'none'}")
T5_ARGON_REPEATS = int(rep.loc["A"].sum())
print(f"  Argon repeats {T5_ARGON_REPEATS} of its 3 columns across a 13 K interval; "
      f"nitrogen repeats {int(rep.loc['N2'].sum())}.")
print("  Printed as read, reported, NOT repaired. Nothing on this page reads the -195.8 C")
print("  argon row: section 5 merges Tables III and V on (gas, temperature) and Table III")
print("  has no argon at -195.8 C, so the row is unused rather than silently averaged in.")

# %% [markdown]
# ## Validation
#
# ### The two headlines, each computed a second and independent way
#
# A break table perturbs an input and watches a number move. That establishes
# sensitivity, never correctness, a baseline that is wrong by accident passes
# every row it has. So both headlines are computed twice by routes that share no
# algebra.
#
# | headline | route 1 | route 2 | agreement |
# |---|---|---|---|
# | the inflection ceiling $2/\sqrt3$ | Lagrange condition on the inflection locus, solved symbolically | bracketing root-find on a numerically differentiated root-find | see below |
# | the polarization decay ratio $C$ | eq. (1b), the closed form the paper quotes | the recursion (1a) assembled and solved as a tridiagonal system, decay ratio measured off the solution | see below |
# | eq. (B) on every branch | the closed form, from two geometric series | the layer equilibria solved for $\{s_i\}$, eq. (15)'s sums formed term by term | see below |
# | the single-point $c$ threshold | closed form $c \ge (1-\varepsilon)(1-x)/(\varepsilon x)$ | bracketing root-find | see below |
#
# All four agree to round-off, which is what an identity between two correct
# routes looks like; they are registered as bounds with the basis that sets them,
# never as digits. The location $c^\ast$ of the ceiling agrees only to what the
# derivative step costs, which section 6 measures on a ladder of steps.
#
# ### What the break table can and cannot reach
#
# The physics here is algebra on transcribed constants, so the defects that matter
# are of three kinds and the table is built around them:
#
# 1. **transcription defects**, a digit changed in a table, the mid-dot decimal
#    trap, a lattice sum swapped;
# 2. **model defects**, a sign in eq. (26), the wrong branch of eq. (B), (1c) used
#    where (1b) belongs, the geometric factor $d$ left out;
# 3. **basis defects**, the wrong rows averaged, the butane row put back into a
#    claim the authors exclude it from.
#
# Three rows are in the table *because they cannot move anything*, and an unstated
# blind spot is an implicit claim. **Two are blind spots; the third is a
# structural identity and is labelled as one, because a check that cannot fail is
# not evidence:**
#
# - *Blind spot.* **Adding layers to the BET ladder does not move the branch
#   agreement**, the ladder is exact at every $n$, not converging towards
#   eq. (B), and the page says so rather than presenting an unmoving number as
#   convergence. What it cannot detect is a convergence failure, because there is
#   no convergence.
# - *Blind spot.* **Perturbing the cross-sections does not move the cross-table
#   $\sigma$ identity** if both tables are perturbed together, the identity tests
#   that one set was used, not what the set was.
# - *Structural, and NOT cited as support anywhere.* **Rescaling every
#   reconstructed $c$ does not move the inflection ceiling**, and it never could:
#   the ceiling is a maximum over *all* $c$. Pointing a $c$-perturbation at it is
#   a powerless check. The rows that actually test what a wrong $c$ would do are
#   the four aimed at `single_point_worst_underestimate_pct`,
#   `t1_min_reconstructed_c` and `inflection_over_pointB_mean`, and section 8
#   reports that the first two **fail** under the factor the paper's own two
#   routes to $c$ exhibit.
#
# Every row is declared before it runs; a row whose study raises records `nan`
# for the claims it targets and stays in the table. The cell ends by listing the
# rows that fail no claim and the claims no row fails.

# %%
# ---- every number the page states, registered with how it is checked
RT = 1e-6     # value claims below are arithmetic on transcribed numbers or root-found to
              # far below this; a misread digit moves them by orders more


def val(name, x, rtol=RT, **kw):
    claims.value(name, x, rtol=rtol, **kw)


for _k, _v in [
        ("polar_d_recomputed", D_CONST), ("polar_C_recomputed", C_POL), ("polar_K1_recomputed", K1_POL),
        ("polar_lattice_net", LATTICE_NET), ("bradley_C_from_eq1c", C_BRAD),
        ("C_factor_bradley_over_computed", C_FACTOR), ("energy_factor_bradley_over_computed", E_FACTOR),
        ("bradley_k_over_eq1b_limit", K_LIMIT_RATIO),
        ("eqB_n6_over_n5_at_058_pct", ERR_058), ("eqB_n6_over_n7_at_072_pct", ERR_072),
        ("eqB_n_from_058", N_FROM_058), ("eqB_n_from_072", N_FROM_072),
        ("bet_inflection_ceiling", CEIL_NUM), ("c_at_inflection_ceiling", C_STAR_NUM),
        ("ceiling_excess_pct", (CEIL_NUM - 1) * 100),
        ("pointB_over_vm_mean", PB_MEAN), ("pointB_over_vm_sd", PB_SD),
        ("inflection_over_pointB_mean", PRED_OVER_OBS), ("inflection_above_pointB_pct", INFL_ABOVE_B),
        ("pointB_below_inflection_pct", PB_BELOW_INFL), ("lowest_pointB_over_vm", PBS["lowest_pointB_over_vm"]),
        ("t1_worst_vm_vs_B_pct", T1_WORST_VM_B), ("t3_butane_vm_over_B", T3_BUTANE_RATIO),
        ("single_point_c_threshold", C_THRESH_CLOSED), ("single_point_worst_underestimate_pct", SHORT_WORST),
        ("t1_min_reconstructed_c", C_T1_MIN),
        ("t2_band_centre_max_dev_cal", T2_CENTRE_DEV), ("t2_band_half_max_dev_cal", T2_HALF_DEV),
        ("t2_CO2_half_range_cal", T2_CO2_HALF_RANGE),
        ("t1_band_twelve_row_centre_cal", MID12), ("t1_band_twelve_row_half_cal", HALF12),
        ("t1_band_eleven_row_centre_cal", MID11), ("t1_band_eleven_row_half_cal", HALF11),
        ("t3_solid_mean_m2g", T3_SOLID_MEAN), ("t3_liquid_maxdev_pct", T3_LIQ_DEV),
        ("t3_vm_maxdev_pct", T3_VM_DEV), ("t5_solid_mean_m2g", T5_SOLID_MEAN),
        ("t5_liquid_maxdev_pct", T5_LIQ_DEV), ("t5_vm_maxdev_pct", T5_VM_DEV),
        ("sigma_cross_table_max_pct", SIGMA_CROSS_MAX),
        ("sigma_cross_table_max_over_envelope", SIGMA_OVER_ENVELOPE),
        ("sigma_worst_envelope_pct", float(SIG_WORST_ROW["env_liquid_%"])),
        ("discrimination_factor_silica", DISC_SILICA), ("discrimination_factor_charcoal", DISC_CHARCOAL),
        ("t4_vm_ratio", T4_VM_RATIO), ("t4_dL_power_ratio", T4_DL_RATIO),
        ("t4_E1_closure_max_cal", T4_CLOSURE_MAX),
        ("vm_temperature_extrapolation_pct", VM_TEMP_DEV), ("c_two_routes_ratio_77K", C_TWO_ROUTE),
        ("c_two_routes_factor_77K", 1 / C_TWO_ROUTE), ("c_fig4_at_77K", C_FIG4_AT_773),
        ("E_from_fig3_c_cal", E_from_c), ("c_t1_954", C_T1_954),
        ("c_rounding_cost_pct", 0.5 / (R_CAL * T_TABLE1) * 100),
        ("single_point_worst_underestimate_scaled_pct", SHORT_WORST_SCALED),
        ("t1_min_reconstructed_c_scaled", C_T1_MIN_SCALED),
        ("sigma_tuned_charcoal_spread_pct", SIGMA_TUNED_CHARCOAL),
        ("sigma_printed_charcoal_6_pct", SIGMA_PRINTED_CHARCOAL_6), ("sigma_tuned_power", SIGMA_TUNED_POWER),
        ("sigma_tuned_charcoal_7_pct", SIGMA_TUNED_CHARCOAL_7),
        ("sigma_printed_charcoal_7_pct", SIGMA_PRINTED_CHARCOAL_7),
        ("palmer_clark_ratio_liquid", PC_LIQUID), ("palmer_clark_ratio_adam", PC_ADAM),
        ("t1_954_two_places_pct", T1_954_DEV)]:
    val(_k, _v)
# the truncation spread at Bradley's k comes from near-singular solves (the recursion has
# no decaying solution there), so its digits follow the linear algebra; 5 % still separates
# it from the zero at the physical k by orders of magnitude
val("ladder_truncation_spread_bradley", LADDER_SPREAD_BRAD, rtol=0.05)

C_GRID_N = 40                        # nodes of the sampled-maximum study (log grid, c = 10 to 1000)


def ceiling_sampled(n_grid=C_GRID_N):
    """The ceiling and its c read off a log grid of n_grid nodes instead of root-found."""
    grid = np.logspace(1, 3, n_grid)
    vals = np.array([inflection_ratio(cv) for cv in grid])
    return float(vals.max()), float(grid[int(np.argmax(vals))])   # argmax IS the defect studied


CEIL_SAMPLED, C_STAR_SAMPLED = ceiling_sampled()
CEIL_SAMPLED_RELDEV = abs(CEIL_SAMPLED / CEIL_NUM - 1.0)
C_STAR_SAMPLED_PCT = (C_STAR_NUM / C_STAR_SAMPLED - 1.0) * 100.0
for _k, _v in [("ceiling_sampled_grid", CEIL_SAMPLED), ("ceiling_sampled_grid_reldev", CEIL_SAMPLED_RELDEV),
               ("c_star_sampled_grid", C_STAR_SAMPLED), ("c_star_sampled_grid_pct", C_STAR_SAMPLED_PCT)]:
    val(_k, _v)
print(f"SAMPLED vs ROOT-FOUND, on a {C_GRID_N}-node log grid:")
print(f"  ceiling : sampled {CEIL_SAMPLED:.7f}  root-found {CEIL_NUM:.7f}  -> "
      f"{CEIL_SAMPLED_RELDEV:.3e} relative, in the sixth significant figure.")
print(f"  c*      : sampled {C_STAR_SAMPLED:.6f}  root-found {C_STAR_NUM:.6f}  -> "
      f"{C_STAR_SAMPLED_PCT:.4f} % off; the grid has a node at {C_STAR_SAMPLED:.4f} and the")
print("            argmax lands on it. The maximum is quadratic, so its LOCATION degrades far")
print("            faster than its VALUE.\n")

# bounds: identities, two-route gaps and orders
claims.bound("ladder_vs_eq1b_reldev", LADDER_VS_1B, max=1e-12, basis="identity",
             note="eq. (1b) solves (1a) exactly for k <= 1/2; the gap is round-off of the dense "
                  "solve (C^400 truncation is far below it)")
claims.bound("ladder_truncation_spread_physical", LADDER_SPREAD_PHYS, max=1e-12, basis="identity",
             note="at the physical k the decay ratio does not depend on where the ladder is cut")
claims.bound("symbolic_identity_max_resid", SYMBOLIC_MAX_RESID, max=0.0, basis="identity",
             note="sympy simplifies each of the five identities to exactly zero")
claims.bound("ladder_vs_eqB_max_reldev", LADDER_VS_EQB, max=1e-12, basis="tolerance",
             note="newton step tolerance 1e-12 on occupancies scaled by c x, and round-off")
claims.bound("bet_inflection_ceiling_two_route_reldev", CEIL_TWO_ROUTE, max=1e-12, basis="tolerance",
             note="a quadratic maximum carries the square of the c* error, so round-off sets it")
claims.bound("c_star_two_route_reldev", C_STAR_TWO_ROUTE, max=1e-9, basis="tolerance",
             note="central difference at h = 1e-5: truncation 0.71 h^2 = 7e-11 measured on the "
                  "step ladder, plus round-off of the difference quotient of the same size")
claims.bound("c_star_fd_order", C_STAR_FD_ORDER, min=1.9, max=2.1, basis="theory",
             note="central difference: c* error of order h^2, steps 1e-2, 1e-3, 1e-4")
claims.bound("c_star_ladder_over_floor", C_STAR_LADDER_OVER_FLOOR, min=100.0, basis="chosen",
             note="the finest successive difference of the step ladder must sit far above the "
                  "gap left at the page's step, or the order measures round-off")
claims.bound("single_point_threshold_two_route_reldev", C_THRESH_TWO_ROUTE, max=1e-12,
             basis="tolerance", note="Brent xtol 1e-12 on c = 38 and round-off")
claims.bound("transcription_max_dev_cal", TRANSCRIPTION_MAX_DEV, max=0.0, basis="identity",
             note="five integer cells printed twice in the paper must be equal")

# counts and exact results (integers from integer data: no noise to keep a margin against)
claims.count("t1_rows_outside_printed_band", np.abs(E_T1 - P["t1_N2_band_centre"]),
             P["t1_N2_band_half"], min_margin=0.0, relative=False,
             note="Table I rows outside the printed 840 +- 70")
claims.exact("t1_omissions_both_rules", T1_OMISSIONS_BOTH,
             note="single-row omissions of Table I reproducing 840 +- 70 under both half-width rules")
claims.exact("t2_centre_hits", int(CENTRE_HITS), note="Table II centres equal to the midrange to 10 cal")
claims.exact("t2_half_hits", int(HALF_MULT5),
             note="Table II half-widths a multiple of 5 within 2.5 cal of the half-range")
claims.exact("sigma_pairs_inside_envelope", SIG["sigma_pairs_inside_envelope"],
             note="recovered-sigma pairs inside their own rounding envelope, of 14")
claims.exact("t1_rows_breaching_scaled", T1_ROWS_BREACHING_SCALED,
             note="Table I rows past 5 % with every c scaled by the two-route ratio")
claims.exact("t1_rows_vm_vs_B_over_tolerance", T1_ROWS_OVER_10, note="Table I rows with |v_m - B|/B above 10 %")
claims.exact("t3_rows_below_threshold", n_below, note="Table III isotherms with c below the threshold 38")
print(f"{len(claims.metrics)} claims registered")

# %%
# ---- break rows: each injects one defect through an input or a call argument and
# records what the claims it targets become; a study that raises records nan


def polar(**kw):
    """Section 1's chain from the printed constants P, with any of them overridden."""
    q = dict(P, **kw)
    net = q["polar_lattice_adj_quad"] - q["polar_lattice_adj_cube"]
    d = net / (1.0 + q["polar_lattice_self"] * q["polar_alpha_over_r3"])
    cp = d * q["polar_alpha_over_r3"]
    cb = C_eq1c(q["bradley_k_Al2SO4"])
    return {"polar_lattice_net": net, "polar_d_recomputed": d, "polar_C_recomputed": cp,
            "polar_K1_recomputed": cp**2, "bradley_C_from_eq1c": cb,
            "C_factor_bradley_over_computed": cb / abs(cp),
            "energy_factor_bradley_over_computed": (cb / abs(cp))**2,
            "bradley_k_over_eq1b_limit": q["bradley_k_Al2SO4"] / q["bradley_k_upper_limit"]}


def eqB_pct(xx, ref_n, cc=None, closed=v_eqB):
    """Error % of the n = 6 curve (from `closed`) against the ref_n branch of eq. (B)."""
    cc = C3 if cc is None else cc
    return (float(closed(xx, cc, 6, VM3)) / float(v_eqB(xx, cc, ref_n, VM3)) - 1) * 100


def v_eqB_wrong(xx, cc, nn, vm=1.0):
    """Eq. (B) with the (n+1) of the numerator mistyped as n."""
    return vm * cc * xx / (1 - xx) * (1 - nn * xx**nn + nn * xx**(nn + 1)) / (1 + (cc - 1) * xx - cc * xx**(nn + 1))


C_SCAN = (2.05, 1e4)   # c scanned for the maximum of an inflection value (log grid)


def ceiling_of(expr):
    """Inflection ceiling of an isotherm form expr(x, c) (sympy): the inflection is the root
    of the numerator of its second derivative on (0, 1), maximised over c on C_SCAN by a
    root of the central difference with step H_FD. On eq. (26) it returns 2/sqrt(3)."""
    f = sp.lambdify((x, c), expr)
    d2 = sp.lambdify((x, c), sp.numer(sp.together(sp.diff(expr, x, 2))))
    xi = lambda cc: unique_root(lambda t: d2(t, cc), X_EDGE, 1 - X_EDGE)
    ir = lambda cc: f(xi(cc), cc)
    cs = unique_root(lambda cc: (ir(cc * (1 + H_FD)) - ir(cc * (1 - H_FD))) / (2 * H_FD * cc),
                     *C_SCAN, log=True, xtol=1e-12)
    return float(ir(cs))


def fig_stats(**kw):
    """The figure and Palmer-and-Clark comparisons from P, with any constant overridden."""
    q = dict(P, **kw)
    c4 = float(np.exp(q["fig4_N2_1958_E"] / (R_CAL * q["fig3_T"])))
    return {"c_fig4_at_77K": c4, "c_two_routes_ratio_77K": q["fig3_c"] / c4,
            "c_two_routes_factor_77K": c4 / q["fig3_c"],
            "E_from_fig3_c_cal": R_CAL * q["fig3_T"] * np.log(q["fig3_c"]),
            "vm_temperature_extrapolation_pct": (q["fig4_N2_1958_vm"] / q["fig3_vm"] - 1) * 100,
            "palmer_clark_ratio_liquid": q["pc_area_from_liquid"] / q["pc_area_HF"],
            "palmer_clark_ratio_adam": q["pc_area_from_adam"] / q["pc_area_HF"],
            "t1_954_two_places_pct": (vm954 * q["fig3_mass"] / q["fig4_N2_183_vm"] - 1) * 100}


def with_cell(t, row_i, col, v):
    """A copy of table t with one cell changed (a transcription defect)."""
    t_ = t.copy()
    t_.loc[t_.index[row_i], col] = v
    return t_


def transcription_dev(t1_=t1, t3_=t3):
    """Worst disagreement of the five cells printed twice (Tables I/II and II/III)."""
    d = [t2.loc[(t2.substance == "Silica gel") & (t2.gas == g), "E1_minus_EL_cal_per_mol"].item()
         - t3_.loc[(t3_.gas == g) & (t3_.temperature_C == T), "E1_minus_EL_cal_per_mol"].item()
         for g, T in [("N2", -183), ("A", -183), ("CO2", -78), ("C4H10", 0)]]
    d.append(t2.loc[(t2.substance.str.contains("954")) & (t2.gas == "N2"), "E1_minus_EL_cal_per_mol"].item()
             - t1_.loc[t1_.substance == "Fe-Al2O3 catalyst 954", "E1_minus_EL_cal_per_mol"].item())
    return float(np.max(np.abs(d)))


def t4_stats(t4_=t4):
    v, dl = t4_.v_m_cc_per_g.to_numpy(float), t4_.d_L_g_per_cc.to_numpy(float)
    cl = (t4_.E1_cal_per_mol - t4_.E1_minus_EL_cal_per_mol - t4_.EL_cal_per_mol).dropna()
    return {"t4_vm_ratio": float(v[-1] / v[0]), "t4_dL_power_ratio": float((dl[-1] / dl[0])**(2 / 3)),
            "t4_E1_closure_max_cal": float(cl.abs().max())}


def band_stats(E):
    m12, h12 = band(E)
    return {"t1_band_twelve_row_centre_cal": m12, "t1_band_twelve_row_half_cal": h12,
            "t1_rows_outside_printed_band": int((np.abs(E - P["t1_N2_band_centre"]) > P["t1_N2_band_half"]).sum()),
            "t1_omissions_both_rules": omissions_under_both(E)}


def t2_stats(t2_):
    cv = t2_bands(t2_)
    return {"t2_band_centre_max_dev_cal": float(cv["centre error"].abs().max()),
            "t2_band_half_max_dev_cal": float(cv["half error"].abs().max()),
            "t2_CO2_half_range_cal": float(cv.loc[cv.gas == "CO2", "half-range"].item()),
            "t2_centre_hits": int((cv.midrange.map(to10) == cv["printed centre"]).sum()),
            "t2_half_hits": int(((cv["printed half"] % 5.0 == 0)
                                 & ((cv["printed half"] - cv["half-range"]).abs() <= 2.5)).sum())}


def surface_stats(t3_=t3, t5_=t5, k5_=k5):
    """Area and v_m spreads, the discrimination factors and the recovered-sigma comparison."""
    out = {"t3_solid_mean_m2g": spread(t3_.surface_solid_m2_per_g)[0],
           "t3_liquid_maxdev_pct": spread(t3_.surface_liquid_m2_per_g)[1],
           "t3_vm_maxdev_pct": spread(t3_.v_m_cc_per_g)[1],
           "t5_solid_mean_m2g": spread(t5_.surface_solid_m2_per_g[k5_])[0],
           "t5_liquid_maxdev_pct": spread(t5_.surface_liquid_m2_per_g[k5_])[1],
           "t5_vm_maxdev_pct": spread(t5_.v_m_cc_per_g[k5_])[1]}
    out["discrimination_factor_silica"] = out["t3_vm_maxdev_pct"] / out["t3_liquid_maxdev_pct"]
    out["discrimination_factor_charcoal"] = out["t5_vm_maxdev_pct"] / out["t5_liquid_maxdev_pct"]
    m = sigma_cross(t3_, t5_)
    out.update(sigma_stats(m))
    out["sigma_worst_envelope_pct"] = float(m.loc[m["over_liquid"].idxmax(), "env_liquid_%"])
    sh = t3_.merge(t5_, on=["gas", "temperature_C"], suffixes=("_s", "_c"))
    out["sigma_printed_charcoal_6_pct"] = spread(sh[sh.gas != "C4H10"].surface_liquid_m2_per_g_c)[1]
    out["sigma_printed_charcoal_7_pct"] = spread(sh.surface_liquid_m2_per_g_c)[1]
    return out


def c_t3_below(R=R_CAL):
    cc = np.exp(t3.E1_minus_EL_cal_per_mol.to_numpy(float) / (R * (t3.temperature_C.to_numpy(float) + KELVIN)))
    return int((cc < C_THRESH_CLOSED).sum())


def pick(fn, *keys):
    """Several targets from one study: {key: thunk}, the study run once per row."""
    memo = {}
    def get(k):
        if "r" not in memo:
            memo["r"] = fn()
        return memo["r"][k]
    return {k: (lambda k=k: get(k)) for k in keys}


T1_ROW7_VM = with_cell(t1, 6, "v_m_cc_per_g", 1.4)                          # 0.14 -> 1.4
SCALE_LABEL = f"every reconstructed c scaled by the Fig.3/Fig.4 ratio {C_TWO_ROUTE:.4f}"
ROWS = [
    ("lattice sums swapped: 8.823 - 8.357", pick(
        lambda: polar(polar_lattice_adj_quad=P["polar_lattice_adj_cube"],
                      polar_lattice_adj_cube=P["polar_lattice_adj_quad"]),
        "polar_lattice_net", "polar_d_recomputed", "polar_C_recomputed")),
    ("alpha/r^3 read as 0.29 (lost zero)", pick(
        lambda: polar(polar_alpha_over_r3=0.29),
        "polar_C_recomputed", "polar_K1_recomputed", "C_factor_bradley_over_computed")),
    ("self-term 11.1 dropped from the denominator", pick(
        lambda: polar(polar_lattice_self=0.0), "polar_d_recomputed")),
    ("Bradley's k read as 0.0615 (misplaced point)", pick(
        lambda: polar(bradley_k_Al2SO4=0.0615), "bradley_C_from_eq1c", "C_factor_bradley_over_computed",
        "energy_factor_bradley_over_computed", "bradley_k_over_eq1b_limit")),
    ("eq. (1c) used where (1b) belongs", {
        "ladder_vs_eq1b_reldev": lambda: abs(LADDER_C / C_eq1c(K_POL) - 1)}),
    ("ladder built with the mu_{i+1} coupling dropped", {
        "ladder_vs_eq1b_reldev": lambda: abs(np.linalg.solve(
            np.eye(60) - np.diag(np.full(59, K_POL), -1), np.r_[K_POL, np.zeros(59)])[0] / C_eq1b(K_POL) - 1)}),
    ("truncation spread read at the physical k instead of Bradley's", {
        "ladder_truncation_spread_bradley": lambda: trunc_spread(K_POL)}),
    ("ladder run at Bradley's k instead of the physical one", {
        "ladder_truncation_spread_physical": lambda: trunc_spread(k_brad)}),
    ("n = 4 used as the reference branch at x = 0.58", {
        "eqB_n6_over_n5_at_058_pct": lambda: eqB_pct(P["eqB_x_n5"], 4)}),
    ("n = 8 used as the reference branch at x = 0.72", {
        "eqB_n6_over_n7_at_072_pct": lambda: eqB_pct(P["eqB_x_n7"], 8)}),
    ("eq. (A) used where eq. (B) at n = 6 belongs", {
        "eqB_n6_over_n5_at_058_pct": lambda: eqB_pct(
            P["eqB_x_n5"], 5, closed=lambda xx, cc, nn, vm: v_over_vm_A(xx, cc) * vm)}),
    ("Fig. 3's c misread as 15.67", {
        "eqB_n6_over_n5_at_058_pct": lambda: eqB_pct(P["eqB_x_n5"], 5, cc=15.67),
        "eqB_n_from_058": lambda: n_from_pct(P["eqB_x_n5"], 5, P["eqB_err_at_058"], cc=15.67),
        **pick(lambda: fig_stats(fig3_c=15.67), "c_two_routes_ratio_77K", "c_two_routes_factor_77K",
               "E_from_fig3_c_cal")}),
    ("target set to +7 % instead of +5 % at x = 0.58", {
        "eqB_n_from_058": lambda: n_from_pct(P["eqB_x_n5"], 5, 7.0)}),
    ("target set to -5 % instead of -7 % at x = 0.72", {
        "eqB_n_from_072": lambda: n_from_pct(P["eqB_x_n7"], 7, -5.0)}),
    ("eq. (B) numerator (n+1) mistyped as n", {
        "ladder_vs_eqB_max_reldev": lambda: ladder_worst(closed=v_eqB_wrong, ns=(2, 5, 6, 7)),
        "symbolic_identity_max_resid": lambda: float(abs(sp.N(sp.simplify(
            v_eqB_wrong(x, c, n).subs(n, 1) - c * x / (1 + c * x)).subs({x: sp.Rational(1, 3), c: 100}))))}),
    ("ladder compared against the n+1 branch", {
        "ladder_vs_eqB_max_reldev": lambda: ladder_worst(n_shift=1, ns=(2, 5, 6, 7))}),
    ("ladder run at c 1 % high", {
        "ladder_vs_eqB_max_reldev": lambda: ladder_worst(c_ladder=C3 * 1.01, ns=(2, 5, 6, 7))}),
    ("BLIND SPOT: ladder layer count 5 -> 60 (exact at every n, not converging)", {
        "ladder_vs_eqB_max_reldev": lambda: ladder_worst(ns=(5, 6, 60), xs=(0.05, 0.58))}),
    ("eq. (26) denominator (1-x+cx) -> (1+x+cx)", {
        "bet_inflection_ceiling": lambda: ceiling_of(c * x / ((1 - x) * (1 + x + c * x)))}),
    ("ceiling read at c = 25 instead of maximised over c", {
        "bet_inflection_ceiling": lambda: inflection_ratio(25.0),
        "ceiling_excess_pct": lambda: (inflection_ratio(25.0) - 1) * 100}),
    ("ceiling maximum SAMPLED on a 40-point log grid, not root-found", {
        "bet_inflection_ceiling": lambda: ceiling_sampled(C_GRID_N)[0]}),
    ("ceiling maximum sampled on a 40-point log grid", {
        "c_at_inflection_ceiling": lambda: ceiling_sampled(C_GRID_N)[1]}),
    ("the sampling grid refined from 40 to 400 points", {
        "ceiling_sampled_grid": lambda: ceiling_sampled(400)[0],
        "ceiling_sampled_grid_reldev": lambda: abs(ceiling_sampled(400)[0] / CEIL_NUM - 1),
        "c_star_sampled_grid": lambda: ceiling_sampled(400)[1],
        "c_star_sampled_grid_pct": lambda: (C_STAR_NUM / ceiling_sampled(400)[1] - 1) * 100}),
    ("ceiling compared against 2/sqrt(2) instead of 2/sqrt(3)", {
        "bet_inflection_ceiling_two_route_reldev": lambda: abs(CEIL_NUM / (2 / np.sqrt(2)) - 1)}),
    ("c* derivative step 1e-5 -> 1e-3", {
        "c_star_two_route_reldev": lambda: abs(c_star_root(1e-3) / C_STAR_EXACT - 1)}),
    ("c* derivative taken as a forward difference", {
        "c_star_fd_order": lambda: c_star_study(central=False)["orders"][-1]}),
    ("c* step ladder moved into round-off: 1e-5, 1e-6, 1e-7", {
        "c_star_ladder_over_floor": lambda: (lambda st: abs(st["values"][-2] - st["values"][-1])
                                             / abs(C_STAR_NUM - C_STAR_EXACT))(
            c_star_study(h_ladder=(1e-5, 1e-6, 1e-7)))}),
    ("STRUCTURAL, cannot move: every reconstructed c scaled by the Fig.3/Fig.4 ratio "
     "-- the ceiling is a maximum over ALL c", {
        "bet_inflection_ceiling": lambda: ceiling_with_c_scale(C_TWO_ROUTE)}),
    (SCALE_LABEL, {
        **pick(lambda: single_point(c_scale=C_TWO_ROUTE), "worst", "c_min"),
        "inflection_over_pointB_mean": lambda: pointB_stats(c_scale=C_TWO_ROUTE)["inflection_over_pointB_mean"]}),
    ("the scaled-c sensitivity run at prefactor 1 (i.e. not run at all)", pick(
        lambda: single_point(c_scale=1.0), "worst", "c_min", "rows_past")),
    ("Table I row 7 v_m 0.14 read as 1.4", pick(
        lambda: pointB_stats(t1_=T1_ROW7_VM), "pointB_over_vm_mean", "pointB_over_vm_sd",
        "inflection_over_pointB_mean", "t1_worst_vm_vs_B_pct", "inflection_above_pointB_pct",
        "pointB_below_inflection_pct")),
    ("Table I Commercial Cu v_m 0.09 read as 0.08", pick(
        lambda: pointB_stats(t1_=with_cell(t1, 8, "v_m_cc_per_g", 0.08)), "t1_rows_vm_vs_B_over_tolerance")),
    ("butane row put back in (the authors exclude it)", pick(
        lambda: pointB_stats(with_butane=True), "pointB_over_vm_mean", "inflection_over_pointB_mean",
        "inflection_above_pointB_pct", "pointB_below_inflection_pct")),
    ("Table III butane point B 28.1 read as 58.1", pick(
        lambda: pointB_stats(t3_=with_cell(t3, 6, "point_B_cc_per_g", 58.1)),
        "t3_butane_vm_over_B", "lowest_pointB_over_vm")),
    ("window taken at x = 0.05 instead of 1/3", {
        "single_point_c_threshold": lambda: c_threshold(xv=0.05),
        "single_point_worst_underestimate_pct": lambda: single_point(xv=0.05)["worst"]}),
    ("tolerance taken as 1 % instead of 5 %", {
        "single_point_c_threshold": lambda: c_threshold(eps=0.01)}),
    ("root-found route given x = 0.35 instead of 1/3", {
        "single_point_threshold_two_route_reldev": lambda: abs(unique_root(
            lambda cv: shortcut_ratio(cv, 0.35) - (1 - EPS), C_LO, C_HI, log=True, xtol=1e-12)
            / C_THRESH_CLOSED - 1)}),
    ("R taken as 8.314 (J) instead of 1.9872 (cal)", {
        **pick(lambda: single_point(R=8.314), "c_min", "worst"),
        "c_rounding_cost_pct": lambda: 0.5 / (8.314 * T_TABLE1) * 100,
        "c_t1_954": lambda: float(np.exp(E_T1_954 / (8.314 * T_TABLE1))),
        "t3_rows_below_threshold": lambda: c_t3_below(R=8.314)}),
    ("Table II CO2 silica 1335 -> 1385", pick(
        lambda: t2_stats(with_cell(t2, int(np.flatnonzero((t2.gas == "CO2") & (t2.substance == "Silica gel"))[0]),
                                   "E1_minus_EL_cal_per_mol", 1385)),
        "t2_band_centre_max_dev_cal", "t2_band_half_max_dev_cal", "t2_CO2_half_range_cal",
        "t2_centre_hits", "t2_half_hits")),
    ("Table I Cr2O3 gel 738 -> 838", pick(
        lambda: band_stats(np.where(t1.substance == "Cr2O3 gel", 838.0, E_T1)),
        "t1_band_twelve_row_centre_cal", "t1_band_twelve_row_half_cal", "t1_rows_outside_printed_band",
        "t1_omissions_both_rules")),
    ("catalyst 931 omitted instead of Cr2O3 gel", {
        "t1_band_eleven_row_centre_cal": lambda: band(np.delete(E_T1, 4))[0],
        "t1_band_eleven_row_half_cal": lambda: band(np.delete(E_T1, 4))[1]}),
    ("Table III N2 -183 E1-EL 794 -> 749 (transposition)", {
        "transcription_max_dev_cal": lambda: transcription_dev(
            t3_=with_cell(t3, 1, "E1_minus_EL_cal_per_mol", 749))}),
    ("Table III butane solid surface 504 -> 405", pick(
        lambda: surface_stats(t3_=with_cell(t3, 6, "surface_solid_m2_per_g", 405)), "t3_solid_mean_m2g")),
    ("Table III argon liquid surface 464 -> 646", pick(
        lambda: surface_stats(t3_=with_cell(t3, 2, "surface_liquid_m2_per_g", 646)),
        "t3_liquid_maxdev_pct", "discrimination_factor_silica")),
    ("Table III butane v_m 58.2 -> 85.2 (transposition)", pick(
        lambda: surface_stats(t3_=with_cell(t3, 6, "v_m_cc_per_g", 85.2)),
        "t3_vm_maxdev_pct", "discrimination_factor_silica")),
    ("butane row kept in the charcoal means (the paper drops it)", pick(
        lambda: surface_stats(k5_=np.ones(len(t5), bool)), "t5_solid_mean_m2g", "t5_liquid_maxdev_pct",
        "discrimination_factor_charcoal", "t5_vm_maxdev_pct")),
    ("Table V N2 -183 v_m 173.0 -> 137.0 (transposition)", pick(
        lambda: surface_stats(t5_=with_cell(t5, 1, "v_m_cc_per_g", 137.0)),
        "sigma_cross_table_max_pct", "sigma_cross_table_max_over_envelope", "sigma_pairs_inside_envelope")),
    ("Table V CO2 -78 liquid surface 853 -> 835 (transposition)", pick(
        lambda: surface_stats(t5_=with_cell(t5, 6, "surface_liquid_m2_per_g", 835)),
        "sigma_printed_charcoal_6_pct", "sigma_printed_charcoal_7_pct", "t5_liquid_maxdev_pct")),
    ("BLIND SPOT: V_STP 22414 -> 24450 (a common factor, cancels in the ratio)", {
        "sigma_cross_table_max_pct": lambda: sigma_stats(sigma_cross(t3, t5, v_stp=24450.0))[
            "sigma_cross_table_max_pct"]}),
    ("surface columns assumed to carry one decimal, not to be integers", {
        "sigma_cross_table_max_over_envelope": lambda: sigma_stats(sigma_cross(t3, t5, du_S=0.05))[
            "sigma_cross_table_max_over_envelope"],
        "sigma_worst_envelope_pct": lambda: (lambda m: float(m.loc[m["over_liquid"].idxmax(), "env_liquid_%"]))(
            sigma_cross(t3, t5, du_S=0.05))}),
    ("tuned-sigma bound taken on the butane-INCLUSIVE seven shared pairs", {
        "sigma_tuned_charcoal_spread_pct": lambda: tuned_spread(sh7),
        "sigma_tuned_power": lambda: tuned_spread(sh7) / SIGMA_PRINTED_CHARCOAL_7}),
    ("tuned sigma taken as a constant instead of proportional to 1/v_m", {
        "sigma_tuned_charcoal_spread_pct": lambda: tuned_spread(sh6, power=0.0),
        "sigma_tuned_charcoal_7_pct": lambda: tuned_spread(sh7, power=0.0)}),
    ("Table IV v_m(-80) 131.0 -> 113.0", pick(
        lambda: t4_stats(with_cell(t4, 5, "v_m_cc_per_g", 113.0)), "t4_vm_ratio")),
    ("Table IV d_L(-80) 1.642 -> 1.462", pick(
        lambda: t4_stats(with_cell(t4, 5, "d_L_g_per_cc", 1.462)), "t4_dL_power_ratio")),
    ("Table IV EL(0 C) 5840 -> 5480", pick(
        lambda: t4_stats(with_cell(t4, 2, "EL_cal_per_mol", 5480.0)), "t4_E1_closure_max_cal")),
    ("Fig. 3 v_m 133.0 read as 138.0", pick(
        lambda: fig_stats(fig3_vm=138.0), "vm_temperature_extrapolation_pct")),
    ("Fig. 4 E1-EL 900 read as 700", pick(
        lambda: fig_stats(fig4_N2_1958_E=700.0), "c_two_routes_ratio_77K", "c_two_routes_factor_77K",
        "c_fig4_at_77K")),
    ("Palmer & Clark HF area 4690 -> 4960", pick(
        lambda: fig_stats(pc_area_HF=4960.0), "palmer_clark_ratio_liquid", "palmer_clark_ratio_adam")),
    ("Fig. 3 sample mass 50.4 -> 54.0 g", pick(
        lambda: fig_stats(fig3_mass=54.0), "t1_954_two_places_pct")),
    ("eq. (28) denominator (p0 - p) mistyped as (p0 + p)", {
        "symbolic_identity_max_resid": lambda: float(abs(sp.N(sp.simplify(
            p / ((vm * c * p / ((p0 + p) * (1 + (c - 1) * (p / p0)))) * (p0 - p))
            - (1 / (vm * c) + (c - 1) / (vm * c) * p / p0)).subs({p: sp.Rational(1, 3), p0: 1, c: 100, vm: 1}))))}),
]
# keys of the single_point() and pointB_stats() results, mapped to the claims they are
SP_KEYS = {SCALE_LABEL: {"worst": "single_point_worst_underestimate_pct", "c_min": "t1_min_reconstructed_c"},
           "the scaled-c sensitivity run at prefactor 1 (i.e. not run at all)": {
               "worst": "single_point_worst_underestimate_scaled_pct",
               "c_min": "t1_min_reconstructed_c_scaled", "rows_past": "t1_rows_breaching_scaled"},
           "R taken as 8.314 (J) instead of 1.9872 (cal)": {
               "c_min": "t1_min_reconstructed_c", "worst": "single_point_worst_underestimate_pct"}}
SENSITIVITY = ["BLIND SPOT: ladder layer count 5 -> 60 (exact at every n, not converging)",
               "BLIND SPOT: V_STP 22414 -> 24450 (a common factor, cancels in the ratio)",
               "STRUCTURAL, cannot move: every reconstructed c scaled by the Fig.3/Fig.4 ratio "
               "-- the ceiling is a maximum over ALL c"]

claims.define_rows([r[0] for r in ROWS])
BREAK = []
for row, studies in ROWS:
    for k_, fn in studies.items():
        name = SP_KEYS.get(row, {}).get(k_, k_)
        try:
            out = fn()
        except Exception as exc:                 # noqa: BLE001 - a failed study is a result
            print(f"  {row}: {name} raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
            out = np.nan
        caught_ = claims.breaks(row, name, out)
        base_ = claims.metrics[name]["value"]
        if claims.metrics[name]["kind"] == "bound":   # round-off digits differ between machines
            base_, out = f"(here {base_:.1e})", f"(here {float(out):.1e})"
        BREAK.append((row, name, base_, out, caught_))

bt = pd.DataFrame(BREAK, columns=["injected defect", "metric", "baseline", "with defect", "caught"])
pd.set_option("display.max_colwidth", 60)
pd.set_option("display.max_rows", 300)
print(bt.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
idle = [r for r in bt["injected defect"].unique()
        if not bt[(bt["injected defect"] == r) & bt.caught].shape[0]]
print(f"\n{len(ROWS)} rows defined, {bt['injected defect'].nunique()} in the table; "
      f"rows that fail no claim: {[r for r in idle if r not in SENSITIVITY] or 'none'}; "
      f"documented non-movers (expected to fail none): {[r for r in idle if r in SENSITIVITY] or 'none'}")
missing = [k for k in claims.metrics if not bt[(bt.metric == k) & bt.caught].shape[0]]
print(f"claims no row fails: {missing or 'none'}")

# %% [markdown]
# ### Claims and the prose
#
# Every number the prose of this notebook and of the README states is checked
# below against the claim it quotes, so a number in the text cannot drift from
# the computation and a cut in its precision shows in CI (agreement.json
# "quoted"). Bounds are quoted as bounds.

# %%
SAYS = {
    "polar_C_recomputed": ["-0.010223", "-0.0102"], "polar_d_recomputed": ["-0.352523"],
    "polar_K1_recomputed": ["1.045e-4"], "polar_lattice_net": ["-0.466"],
    "bradley_C_from_eq1c": ["0.989"], "C_factor_bradley_over_computed": ["96.8"],
    "energy_factor_bradley_over_computed": ["9361"], "bradley_k_over_eq1b_limit": ["1.23"],
    "ladder_vs_eq1b_reldev": ["below 1e-12"], "ladder_truncation_spread_physical": ["below 1e-12"],
    "ladder_truncation_spread_bradley": ["4.48"], "ladder_vs_eqB_max_reldev": ["below 1e-12"],
    "t3_vm_maxdev_pct": ["46.9"], "t3_liquid_maxdev_pct": ["10.6"],
    "discrimination_factor_silica": ["4.42"], "t5_vm_maxdev_pct": ["18.6"],
    "t5_liquid_maxdev_pct": ["7.9"], "bet_inflection_ceiling": ["1.1547", "1.1547005", "1.155"],
    "ceiling_excess_pct": ["15.5"], "lowest_pointB_over_vm": ["0.4828"],
    "inflection_over_pointB_mean": ["1.148812"], "inflection_above_pointB_pct": ["14.881"],
    "pointB_below_inflection_pct": ["12.954"], "t1_band_twelve_row_centre_cal": ["824.5"],
    "t1_band_twelve_row_half_cal": ["86.5"], "t1_band_eleven_row_centre_cal": ["841.5"],
    "t1_band_eleven_row_half_cal": ["69.5"], "t2_band_centre_max_dev_cal": ["4"],
    "t2_CO2_half_range_cal": ["122.5"], "t2_centre_hits": ["4"], "t2_half_hits": ["4"],
    "t1_rows_outside_printed_band": ["2"], "t1_omissions_both_rules": ["1"],
    "t4_E1_closure_max_cal": ["5"], "t1_954_two_places_pct": ["15.6"],
    "sigma_cross_table_max_pct": ["0.6631", "0.66"], "sigma_worst_envelope_pct": ["0.221"],
    "sigma_cross_table_max_over_envelope": ["3.0"], "sigma_pairs_inside_envelope": ["13"],
    "ceiling_sampled_grid": ["1.1546960"], "ceiling_sampled_grid_reldev": ["3.9e-6"],
    "c_star_sampled_grid_pct": ["1.43"], "c_star_sampled_grid": ["52.2335"],
    "c_at_inflection_ceiling": ["52.9808"], "c_star_two_route_reldev": ["below 1e-9"],
    "c_star_fd_order": ["2.0"], "c_rounding_cost_pct": ["0.28"],
    "c_two_routes_factor_77K": ["2.24"], "c_fig4_at_77K": ["350.36"], "c_t1_954": ["147.4"],
    "E_from_fig3_c_cal": ["776.4"], "t1_rows_breaching_scaled": ["7"],
    "t1_min_reconstructed_c_scaled": ["27.58"], "single_point_c_threshold": ["38"],
    "single_point_worst_underestimate_pct": ["3.14"], "single_point_worst_underestimate_scaled_pct": ["6.76"],
    "sigma_tuned_charcoal_spread_pct": ["14.38"], "sigma_printed_charcoal_6_pct": ["7.37"],
    "sigma_tuned_power": ["1.95"], "sigma_tuned_charcoal_7_pct": ["31.28"],
    "sigma_printed_charcoal_7_pct": ["31.04"], "eqB_n6_over_n5_at_058_pct": ["5.62"],
    "eqB_n6_over_n7_at_072_pct": ["-6.85"], "vm_temperature_extrapolation_pct": ["-2.26"],
    "t1_worst_vm_vs_B_pct": ["16.7"], "t3_rows_below_threshold": ["4"],
    "t1_rows_vm_vs_B_over_tolerance": ["1"], "single_point_threshold_two_route_reldev": ["below 1e-12"],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)
print(f"prose literals: {sum(len(v) for v in SAYS.values())} checked against their claims with claims.says")
claims.report();   # the semicolon keeps the notebook from echoing the payload


# %% [markdown]
# ## What pymrm adds
#
# **Nothing to the BET equation, and the page says so where a reader would expect a
# solver.** Equations (26), (28), (A), (B) and (E) are closed forms; the whole of
# sections 1, 3, 5, 6, 7, 8 and 9 would run with pymrm uninstalled, exactly as on
# `A1.6` and `A1.1`. There is no grid, no time step and no transport on this page.
# What `newton` and `NumJac` are doing here is real but narrow, and it is four
# things.
#
# **1. A route to the isotherm that shares no algebra with the isotherm.** Eq. (B)
# is two geometric series summed by hand in 1938. Solving the layer equilibria
# (10)–(12) for the occupancies and forming eq. (15)'s sums term by term is a
# different computation with the same answer, and it *catches* a wrong closed
# form: mistyping the $(n+1)$ in eq. (B)'s numerator as $n$ moves the agreement
# from round-off (below 1e-12) to a percent-level disagreement, which is the break row
# that gives this check its power. Exercised on **every branch the paper uses**,
# $n = 1, 2, 5, 6, 7$ and $n\to\infty$, at four relative pressures.
#
# **2. The same treatment applied to the rival theory, which is the point.**
# The polarization theory's closed form (1b) is quoted by BET from a paper nobody
# here has read. Assembling the recursion (1a) and measuring the decay ratio off
# the solution removes that dependence entirely, and it does something the closed
# form cannot: at $k > 1/2$ it *shows what goes wrong*, the moments oscillating and
# the answer depending on where the ladder is cut. Eq. (1b) merely returns a
# complex number.
#
# **3. Root-finds where a sweep would have been wrong.** The inflection point, the
# $c$ that maximises it, the $c$ at which the single-point rule reaches exactly
# 5 %, and the continuous $n$ behind each printed percentage. The break table
# carries the row that shows this matters, and the honest size of it is printed
# there: **sampling** the ceiling on a 40-point log grid instead of root-finding it
# moves the ceiling only in the **sixth** significant figure (1.1546960 against
# 1.1547005, $3.9\times10^{-6}$ relative) but moves $c^\ast$ by **1.43 %**
# (52.2335 against 52.9808, because the grid has a node at 52.2335). The ceiling is
# quadratic at its maximum, which is exactly why the *location* degrades a
# thousandfold faster than the *value*, and it is the location a reader would
# quote.
#
# **4. And one thing that is not pymrm at all, and is the most useful output on
# the page: the inflection ceiling.** $2/\sqrt3$ is not in the 1938 paper and, as
# far as this page's reading goes, is not implied anywhere in it. It converts the
# authors' strongest-looking corroboration, $v_\mathrm{m}$ agreeing with point B,
# a claim journal page 315 makes about the **twelve** isotherms of Table I ("the
# two seldom differing by as much as 10 %"), and which this page extends to the
# eighteen non-butane pairs of Tables I and III, into a statement about the
# *shape of their equation*, and it needed a Lagrange condition and a root-find,
# not a solver. One-sidedly: it caps how far the knee can sit *above*
# $v_\mathrm{m}$, not below.

# %% [markdown]
# ## Reuse
#
# **Use this page for** the BET equation itself, its finite-layer form (B), the
# Langmuir/(E) limit, and, the part that is not in the textbooks, an honest
# account of what a BET plot does and does not establish.
#
# **The equation, ready to lift.** `v_over_vm_A(x, c)` is eq. (26)/(A);
# `v_eqB(x, c, n, vm)` is eq. (B) and is correct at $n = 1$ (Langmuir), at
# non-integer $n$ (Table IV uses 3.5) and in the limit; `ladder_v_over_vm` is the
# independent route and is the one to copy if you are extending the model, because
# adding a layer-dependent energy means editing one residual line rather than
# re-deriving a series.
#
# **Four things to carry away, in the order they will save you time.**
#
# 1. **A BET plot's agreement inside $0.05<p/p_0<0.35$ is a fit.** Two parameters
#    over a fifth of the range. Report it as a fit. The 1938 paper does not lean on
#    it either, its evidence is the *consistency between gases*, which is a
#    different and much stronger thing.
# 2. **$v_\mathrm{m}$ agreeing with point B is weak corroboration, from above.**
#    The BET isotherm's inflection sits at
#    $\le 2/\sqrt3\,v_\mathrm{m} = 1.155\,v_\mathrm{m}$ for **every** $c$, so an
#    eye-read knee at or below the inflection cannot be *high* by more than 15.5 %.
#    **It can be arbitrarily low**: the bound is one-sided, and on this paper's own
#    butane isotherm (Table III) point B is only 0.4828 of $v_\mathrm{m}$. If you are
#    validating a surface-area method, this is not the check to use, and if you
#    quote the ceiling, quote it one-sidedly.
# 3. **The single-point (one-pressure) BET shortcut has an exact condition:**
#    $c \ge (1-\varepsilon)(1-x)/(\varepsilon x)$, i.e. $c \ge 38$ for 5 % at
#    $p/p_0 = 1/3$, and it can only *under*estimate $v_\mathrm{m}$. Nitrogen at
#    90 K on the adsorbents of Table I passes on the paper's own numbers; argon,
#    oxygen, CO₂ and butane on the silica gel of Table III, four of its seven
#    isotherms, do **not**, which is why the authors scoped the sentence to
#    nitrogen. **Do not read the Table I margin as robust**: section 8 shows that
#    a $c$ scaled by the factor the paper's own two routes disagree over puts seven
#    of the twelve rows past 5 % and the smallest $c$ below the threshold.
# 4. **$E_1-E_\mathrm{L}$ is not an independent datum, it *is* $c$, re-expressed.**
#    Journal page 313 says the energies were obtained *from* $c$ through
#    footnote 16 with the prefactor set to 1, so converting back returns the
#    authors' fit and nothing more, and the prefactor cannot be measured from this
#    paper at all. What the paper *does* show is one isotherm carrying two $c$ a
#    factor 2.24 apart (Fig. 3 fitted against Fig. 4 calculated); treat that, not
#    the four figures a printed $c$ suggests, as the size of the doubt.
#
# **Do not use this page for** anything dynamic. Adsorption *breakthrough* is
# `J1.5`; the linear-driving-force rate law and the column mass balance live there,
# and a BET equilibrium can be dropped into that page's isotherm slot. Nothing here
# is a rate.
#
# **Do not use this page as a source on the polarization theory.** De Boer and
# Zwicker (1929), de Boer (1931, 1932), and Bradley (1936) are all read *through*
# this paper's restatement of them and none was consulted. What the page shows is
# that the theory **as Brunauer, Emmett and Teller state it** cannot carry the
# binding energy Bradley's fitted $k$ assigns it. If you need to cite the
# polarization theory, go to the originals.

# %% [markdown]
# ## References
#
# Brunauer, S., Emmett, P. H. and Teller, E. (1938). Adsorption of gases in
# multimolecular layers. *Journal of the American Chemical Society* **60**(2),
# 309–319. [doi:10.1021/ja01269a023](https://doi.org/10.1021/ja01269a023), **the
# paper, and the only document read.** Received November 19, 1937. Identity
# confirmed from its own title page on a native-resolution render: the running head
# "Feb., 1938 … 309", the contribution line "[Contribution from the Bureau of
# Chemistry and Soils and George Washington University]", the title and the by-line
# "By Stephen Brunauer, P. H. Emmett and Edward Teller". Equations (1)–(28), (A),
# (B), (C), (D), (E), footnotes 3, 10, 16 and 18, Tables I–V and the typeset
# constant blocks inside Figs. 3, 4 and 6 were all transcribed from renders at the
# scan's native 300 ppi, each numeric cropped and re-read at that resolution.
#
# **Cited by the paper, not consulted, and nothing here derives from them.**
#
# de Boer, J. H. and Zwicker, C. (1929). *Z. physik. Chem.* **B3**, 407; de Boer,
# J. H. (1931). *ibid.* **B13**, 134; (1931) **B14**, 149; (1932) **B17**, 161,
# the polarization theory. Read only through BET's restatement (eqs. 1–4a and
# footnote 3).
#
# Bradley, R. S. (1936). *J. Chem. Soc.*, 1467 and 1799, the polarization isotherm
# applied to argon on copper and aluminium sulfate, and the source of the $k$
# values 0.6075 and 0.615 this page tests. Read only through BET's restatement.
#
# Emmett, P. H. and Brunauer, S. (1937). *This Journal* **59**, 1553, the table of
# molecular cross-sectional areas behind the surface columns of Tables III and V.
# Not on disk. The cross-sections used on this page are **recovered by inverting
# those two columns**, and the recovery is checked across the two tables
# independently.
#
# McGavack, J. and Patrick, W. A. (1920). *This Journal* **42**, 946, the SO₂
# isotherms behind Table IV. Palmer, W. G. and Clark, R. E. D. (1935). *Proc. Roy.
# Soc. (London)* **A149**, 360, the HF-dissolution surface areas. Goldman, F. and
# Polanyi, M. (1928). *Z. physik. Chem.* **132**, 321, the ethyl-chloride
# isotherms of Fig. 6, which this page scopes out. McBain, J. W. (1932). *The
# Sorption of Gases by Solids*, Routledge, p. 169, "persorption". *International
# Critical Tables* Vol. I p. 103 and Vol. VII p. 11, the density of solid argon
# and the refractive index of argon gas. **None consulted.** Where this page needs
# one of these quantities it says so and recovers it by inversion.
#
# Langmuir, I., the unimolecular derivation the BET treatment generalises. Cited
# by name throughout the 1938 paper without a reference of its own; not consulted,
# and no Langmuir equation on this page comes from anywhere but eq. (B) at
# $n = 1$.
