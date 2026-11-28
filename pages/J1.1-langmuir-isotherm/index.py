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
# title: "Langmuir 1918: the isotherm that could not have proved its own conclusion"
# description: "Langmuir's isotherm is fitted to his own data by a hand-drawn straight line, so the agreement in Tables II-XVII is a goodness of fit and this page labels it as one everywhere. What can be tested is elsewhere, and there are two of them. Langmuir's own contest against Freundlich is a fair one - two parameters each, the same eleven points - and Freundlich loses by a factor 3.82 in RMS on Table VII and on 13 of the 14 tables that can host the comparison, the single exception being the table Langmuir himself excludes. And the monomolecular conclusion does not rest on the isotherm at all: it rests on beta, the saturation amount against a monolayer count computed from the liquid density, whose worst row leaves only 17 % of headroom. The page's own result is why that had to be so. Langmuir's Case VI, eqs. (29) and (30) of the same paper, is shown by symbolic algebra to be EXACTLY the BET isotherm twenty years early; and every pressure he reached lies below p/p0 = 1.6e-3, where his monolayer and multilayer equations differ by at most 0.16 %, against observational scatter 12 to 1100 times larger. Four printed defects are reported and none repaired, including eq. (15)'s logarithm, which contradicts eq. (16) two lines below it; a fifth pattern, the sigma column running low, is reported separately because it is a truncation convention and not an error."
# categories: [sec:J, struct:S1, tier:T0, data:tier5, phase:gas-solid]
# date: 2026-08-07
# ---
#
# # Langmuir 1918: the isotherm that could not have proved its own conclusion
#
# **Catalog ID:** `J1.1` · **Structures:** `S1` (pointwise algebra) · **Tier:** T0
#
# Everybody meets this paper as one equation,
#
# $$\theta_1 \;=\; \frac{\sigma_1\mu}{1+\sigma_1\mu},
# \qquad\text{equivalently}\qquad
# q \;=\; \frac{ab\,p}{1+ap},
# \tag{9, 31}$$
#
# and as the linear plot that fits it,
#
# $$\frac{p}{q} \;=\; \frac{1}{ab} \;+\; \frac{p}{b}.
# \tag{32}$$
#
# Langmuir obtained $a$ and $b$ by **drawing a straight line as nearly as possible
# through the points** of that plot. So the $q_\mathrm{cal}$ columns of Tables II
# to XVII are a **two-parameter fit compared against the two-parameter fit's own
# training data**. That is a goodness of fit. It is labelled one in every section
# of this page, it is never reported as an agreement metric, and a null baseline
# is printed beside it every time.
#
# The paper's real evidence is somewhere else, and this page is organised around
# finding it.
#
# 1. **Langmuir stages a fair contest and wins it.** Table VII carries two extra
#    columns: Freundlich's $q_F = 8.4\,p^{0.417}$, fitted to the *same* eleven
#    points with the *same* number of parameters. That is a discrimination, not a
#    fit, because it can fail. Refitting **both** forms optimally, Langmuir's
#    drawn line is not the least-squares optimum, and neither is his Freundlich , 
#    the best two-parameter Freundlich is **3.82 times worse in RMS** than the
#    best two-parameter Langmuir on Table VII, and worse on **13 of the 14** tables
#    with enough points to host the comparison. The one exception is Table XIII,
#    which is the table Langmuir himself says does not satisfy eq. (31).
# 2. **The monomolecular conclusion never rested on the isotherm.** It rests on
#    $\beta$: the fitted saturation amount $b$ converted to molecules per square
#    centimetre through eq. (37), divided by the number in a close-packed layer
#    computed from the **liquid density**, a conversion with no adsorption in it
#    and nothing adjusted to make it come out below one. It comes out below one on
#    all twenty entries, and the page reports the honest margin: the tightest row
#    leaves **17 %** of headroom.
# 3. **And it could not have rested on the isotherm.** This is the page's own
#    result and it is not in the paper. Langmuir's **Case VI**, his *multilayer*
#    theory, eqs. (26), (29) and (30), printed in the same paper, is shown here by
#    symbolic algebra to be, with his own stated assumption $\sigma_3=\sigma_4=
#    \dots=\sigma_2$, **exactly the BET isotherm**, twenty years before
#    Brunauer, Emmett and Teller wrote it down. Inverting the paper's own Table XXI
#    through eq. (3) recovers the saturation pressures Langmuir never prints, and
#    **every pressure in the entire paper lies below $p/p_0 = 1.5\times10^{-3}$**
#   , §6.9 computes that over *every* printed pressure, both bulbs and both
#    isotherm forms, not only over the entries that can host a fit , 
#    where his monolayer and multilayer equations differ by at most **0.16 %** , 
#    against observational scatter **12 to 1100 times larger**. The isotherm shape
#    was blind to the question the paper is about.
#
# Four printed defects are reported and none repaired. The sharpest is eq. (15),
# whose logarithm contradicts eq. (16) two lines below it. A fifth pattern, the
# $\sigma$ column running systematically low, is reported *separately*, because
# §6.3 identifies it as a truncation convention rather than an error.
#
# **Source.** Langmuir, I. (1918). *The Adsorption of Gases on Plane Surfaces of
# Glass, Mica and Platinum*. J. Am. Chem. Soc. **40**(9), 1361–1403,
# [doi:10.1021/ja02242a004](https://doi.org/10.1021/ja02242a004). Received June 25,
# 1918. Identity confirmed from its own title page on a native-resolution render , 
# the contribution line "[Contribution from the Research Laboratory of the General
# Electric Co.]", the title, the by-line "By Irving Langmuir", "Received June 25,
# 1918", and the ACS download stamp `jacsat/article-pdf/40/9/1361/` printed down
# the margin of every page. **PDF page 1 opens with the numbered summary of the
# preceding article** (absorption spectra of metals in liquid ammonia); Langmuir's
# title sits below it on the same page, and that is what was read.
#
# `pdfimages -list` reports every page as CCITT-G4 bilevel at **300 ppi native**,
# so pages were rendered at 300 ppi and every numeric cropped and re-read at digit
# scale. **The text layer was not used for any digit**, and that discipline was
# load-bearing here: the text layer renders eq. (37)'s coefficient as
# $25.2\times10^{16}$ where the page prints $25.2\times10^{\mathbf{15}}$, a factor
# ten in every surface coverage on the page, and it mangled half the $b$ column of
# Table X.
#

# %% [markdown]
# ## 1. Background
#
# ### What the paper is, and what this page takes from it
#
# Forty-three journal pages, twenty-five tables. Most of it is an experimental
# campaign: outgassing bulbs, thermal-effusion corrections, and, from journal page
# 1393 on, irreversible chemisorption on platinum, which has no isotherm in it at
# all. The isotherm is a small part.
#
# **Scope decision, taken before transcription.** This page is about **the isotherm
# and the constants derived from it**: the theory of journal pages 1368–1376
# (Cases I–VI, eqs. 1–33) and Tables II–XI, XIII–XXII, which are the mica and glass
# isotherms and everything Langmuir computes from them. Table XXV is carried for
# one reason only, it exercises eq. (37) on a *third* adsorbent, and Tables I,
# XII, XXIII and XXIV (outgassing inventories and the platinum oxygen/carbon
# monoxide bookkeeping) are **scoped out and not transcribed**. They are excellent
# self-checking tables, but they are about the apparatus and about chemisorption,
# not about the isotherm.
#
# That decision was taken the way the brief asks: *prefer whatever lets you test
# rather than restate*. Restating $q_\mathrm{cal}$ would be restating a fit. The
# testable content is the discrimination against Freundlich, the $\beta$ bound, and
# the internal identity web that ties $a$ and $b$ to $\sigma$, $N_0$ and $\beta$
# through five printed equations.
#
# ### The theory, in the order Langmuir builds it
#
# The kinetic picture is one sentence long: molecules striking a bare surface
# condense, molecules on the surface evaporate, and equilibrium is the balance.
# Eq. (1)/(2) give the impingement rate, eq. (4) is the balance, eq. (5) is the
# site inventory, and eqs. (6)–(9) are the isotherm. Then he generalises in five
# directions, each of which he calls a *Case*:
#
# | Case | assumption | result |
# |---|---|---|
# | I | one kind of site, one molecule per site | eq. (9), the Langmuir isotherm |
# | II | several kinds of site, fractions $\beta_1,\beta_2,\dots$ | eq. (19), a sum of Langmuir terms |
# | III | a continuum of site types (amorphous solid) | eq. (21), an integral |
# | IV | $n$ molecules per site | eq. (26), a ratio of two polynomials |
# | V | dissociative adsorption, two sites per molecule | eq. (28), $\eta \propto \sqrt{p}$ |
# | VI | **films more than one molecule thick** | eqs. (26), (29), (30) |
#
# Case VI is the one this page is going to press on. Langmuir writes: *"The problem
# is then identical with that considered under Case IV. Equations 22, 23 and 24 are
# applicable without alteration and the solution of the problem is thus given by
# Equation 26."*, and then, because both series become infinite, he transforms
# eq. (26) into eq. (29), whose coefficients are eq. (30).
#
# He also derives a **rate**: eq. (14) is the ODE for the approach to equilibrium,
# eq. (15) its integral, and eq. (16) the half-period. The paper never measures a
# rate, so this is theory only, but it is checkable against itself, and it does
# not survive the check.
#
# ### The other thing in the paper: what Langmuir is arguing against
#
# Journal page 1375 states the rival explicitly. Freundlich's $q = a p^{1/n}$
# *"agrees very poorly with experiment when the range of pressures is large"*, and
# Langmuir explains why on his own theory: at low pressure eq. (9) is linear
# ($1/n \to 1$), at high pressure it saturates ($1/n \to 0$), so a single power law
# must fail at one end or the other. Then, in Table VII, he *fits Freundlich to his
# own data* and prints the residuals beside his own. That contest is reproduced
# here, with both forms refitted optimally, because neither of his fits is the
# least-squares optimum and a rival should be beaten at its best.
#

# %% [markdown]
# ## 2. The published model
#
# All equation numbers are Langmuir's. Every one below was read off a 300 ppi crop.
#
# **Impingement.** With $p$ in *bars* (dynes per sq. cm.), $R = 83.2\times10^6$ erg
# per degree,
#
# $$\mu \;=\; \frac{p}{\sqrt{2\pi MRT}} \;=\; 43.75\times10^{-6}\,\frac{p}{\sqrt{MT}}.
# \tag{2, 3}$$
#
# **Case I.** With $\alpha$ the condensation coefficient, $\nu_1$ the evaporation
# rate from a full surface, $\theta$ bare and $\theta_1$ covered,
#
# $$\alpha\theta\mu = \nu_1\theta_1, \qquad \theta+\theta_1 = 1,
# \qquad \sigma_1 \equiv \frac{\alpha}{\nu_1},
# \tag{4, 5, 7}$$
#
# $$\frac{N}{N_0}\,\eta \;=\; \theta_1 \;=\; \frac{\sigma_1\mu}{1+\sigma_1\mu},
# \qquad N = 6.06\times10^{23}.
# \tag{9}$$
#
# **Lives, and the approach to equilibrium.**
#
# $$\tau \;=\; \frac{N_0\alpha}{N\nu_1} \;=\; \frac{N_0\sigma_1}{N},
# \qquad
# \eta \;=\; \frac{\tau\mu}{1+\sigma_1\mu},
# \qquad
# \eta_\infty = \frac{\tau}{\sigma_1} = \frac{N_0}{N};
# \tag{10, 11, 13}$$
#
# $$\frac{N_0}{N}\frac{d\theta'}{dt} \;=\; \alpha\mu - (\nu_1+\alpha\mu)\theta',
# \tag{14}$$
#
# $$t \;=\; \frac{N_0}{N\nu_1(1+\sigma_1\mu)}\,\ln\frac{\theta'}{\theta_1-\theta'},
# \qquad
# t_{1/2} \;=\; \frac{\tau\ln 2}{\alpha(1+\sigma_1\mu)}.
# \tag{15, 16}$$
#
# Eq. (15) is quoted **exactly as printed**. §6.7 shows it cannot be right and that
# eq. (16), which Langmuir obtains *from* it, is.
#
# **Case II** (several kinds of site) and the two-stage special case:
#
# $$\frac{N\eta}{N_0} = \frac{\beta_1\sigma_1\mu}{1+\sigma_1\mu}
# + \frac{\beta_2\sigma_2\mu}{1+\sigma_2\mu} + \dots,
# \qquad
# q = b_1 + \frac{a_2b_2\,p}{1+a_2p}.
# \tag{19, 33}$$
#
# **Cases IV and VI** (the multilayer ladder). With $\sigma_n = \alpha_n/\nu_n$,
#
# $$\frac{N\eta}{N_0}
# = \frac{\sigma_1\mu + 2\sigma_1\sigma_2\mu^2 + 3\sigma_1\sigma_2\sigma_3\mu^3 + \dots}
#        {1 + \sigma_1\mu + \sigma_1\sigma_2\mu^2 + \sigma_1\sigma_2\sigma_3\mu^3 + \dots}
# \;=\; \frac{\sigma_1}{1/\mu + a + b\mu + c\mu^2 + \dots},
# \tag{26, 29}$$
#
# $$a = \sigma_1 - 2\sigma_2,\quad
# b = \sigma_2(4\sigma_2 - 3\sigma_3 - \sigma_1),\quad
# c = 2\sigma_2(6\sigma_2\sigma_3 - 2\sigma_3\sigma_4 + \sigma_1\sigma_2
#    - \sigma_1\sigma_3 - 4\sigma_2^2).
# \tag{30}$$
#
# **The fitted form and its linearisation.**
#
# $$q = \frac{ab\,p}{1+ap},
# \qquad
# \frac{p}{q} = \frac{1}{ab} + \frac{p}{b}.
# \tag{31, 32}$$
#
# **What $a$ and $b$ mean.** Dividing (31) by $b$ and comparing with (9), and using
# $1\ \text{cu. mm.} = 4.16\times10^{-8}$ gram molecules at 20 °C and 760 mm,
#
# $$\sigma\mu = ap,\qquad N_0q = bN\eta,\qquad \eta = 4.16\times10^{-8}\,q/s,
# \tag{34, 35, 36}$$
#
# $$\boxed{\;N_0 = 25.2\times10^{15}\,b/s\;}
# \qquad\text{and}\qquad
# \boxed{\;\sigma = 22860\,a\sqrt{MT}\;}
# \tag{37, 38}$$
#
# with $s$ the adsorbent area in sq. cm. **Eq. (37) reads $10^{15}$ on the page and
# $10^{16}$ in the PDF's text layer.** Everything downstream, $N_0$, $\beta$, the
# monomolecular conclusion, is a factor ten different between the two.
#
# **Saturation.** For a liquid in equilibrium with its own vapour, eq. (4) becomes
# $\alpha\mu = \nu$, i.e.
#
# $$\sigma\mu \;=\; 1 \quad\text{at saturation}.
# \tag{39}$$
#
# Eq. (39) is what makes Table XXI invertible, and it is the paper's own statement
# of the relation this page uses to recover $p_0$.
#

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.28.7"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "J1.1-langmuir-isotherm"  # equals the page directory

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
# Imports. sympy ships with Colab; nothing else beyond pymrm is needed.
import warnings
from decimal import Decimal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sp

from gallery_utils import Claims, load_data
from pymrm.checks import observed_orders
from model import (caseVI_series, departure, fit_freundlich, fit_langmuir, freu_basis,
                   lang_basis, march_eq14, profile_fit, profile_rms, q_eq19, q_eq31,
                   q_eq33, q_freundlich, rms, sos_at_p0_factor, t_eq15_corrected,
                   t_eq15_printed, t_half_of_march, unique_root)

warnings.filterwarnings("ignore")
pd.set_option("display.width", 200)
claims = Claims(PAGE.split("-")[0])

RES = {}                   # every number the cells compute and the prose quotes
def keep(k, v):
    RES[k] = float(v)
    return v

# %% [markdown]
# ## 3. Parameters and assumptions
#
# **Units, and the one that trips people.** Langmuir's *bar* is a dyne per square
# centimetre, not the modern $10^5$ Pa. His summary says the work used "pressures
# of 100 bars (approximately 0.1 mm of mercury) or less"; the cell below checks
# that against the definition of a millimetre of mercury and reports the result.
# Amounts $q$ are cubic millimetres of gas at 20 °C and 760 mm.
#
# **Constants used.** Only Langmuir's own, and each is verified against the
# definition it came from rather than assumed:
#
# | symbol | printed | where |
# |---|---|---|
# | $R$ | $83.2\times10^6$ erg deg$^{-1}$ | journal p. 1368 |
# | coefficient of eq. (3) | $43.75\times10^{-6}$ | eq. (3) |
# | $N$ | $6.06\times10^{23}$ | eq. (9) |
# | 1 cu. mm. at 20 °C, 760 mm | $4.16\times10^{-8}$ gram mol | journal p. 1389 |
# | coefficient of eq. (37) | $25.2\times10^{15}$ | eq. (37) |
# | coefficient of eq. (38) | $22860$ | eq. (38) |
# | mica surface $s$ | 5750 sq. cm. | journal p. 1376 |
# | glass surface $s$ | 1966 sq. cm. | journal p. 1386 |
# | platinum surface $s$ | 312 sq. cm. | journal p. 1393 |
#
# Those seven numeric constants are not independent: $43.75\times10^{-6}$ must be
# $1/\sqrt{2\pi R}$, $22860$ must be its reciprocal, $4.16\times10^{-8}$ must be
# $p_0V/RT$ at 20 °C and 760 mm, and $25.2\times10^{15}$ must be $N$ times that.
# **None of the four is free**, and checking them is the cheapest possible test
# that the transcription is right and that eq. (37) really carries $10^{15}$.
#
# **Molecular weights** are taken as Langmuir's contemporaries would have:
# N₂ 28, CH₄ 16, CO 28, Ar 39.9, O₂ 32, CO₂ 44. The one that matters is argon , 
# 39.9, not the modern 39.95, and the difference moves $\sigma$ by 0.06 %, far
# inside the rounding of the printed $\sigma$ column.
#
# **Assumptions this page makes and the paper does not state.**
#
# - *Table XXII's combining rule.* Where a gas has two or three adsorbed $\sigma$
#   at one condition, the paper does not say how they were combined into one ratio.
#   §6.4 measures which mean reproduces the printed cells rather than assuming one.
# - *The mica area applies to every mica table.* Verified rather than assumed:
#   §6.3 recovers $s$ independently from each row.
#
# **What is external.** Exactly one thing: the monolayer counts of journal page
# 1391 ($0.66,\,0.63,\,0.66,\,0.77,\,0.77,\,0.61 \times 10^{15}$ per sq. cm.) come
# from the **liquid densities**, and only nitrogen's molecular volume (35.5 cu. cm.)
# is printed. The other four are **recovered** in §6.5 by inverting the printed
# counts, and are labelled recoveries wherever they appear. That external input is
# precisely what makes $\beta$ a test and not a fit.
#

# %%
# --- the constant chain: none of these four numbers is free ---------------
R = 83.2e6                     # erg per degree, journal p. 1368
NAV = 6.06e23                  # eq. (9)
K3 = 43.75e-6                  # eq. (3)
K37 = 25.2e15                  # eq. (37)   <-- 10^15 on the page, 10^16 in the text layer
K38 = 22860.0                  # eq. (38)
CUMM = 4.16e-8                 # gram molecules per cu. mm. at 20 C and 760 mm, p. 1389
S_CM2 = {"mica": 5750.0, "glass": 1966.0, "platinum": 312.0}
MW = {"N2": 28.0, "CH4": 16.0, "CO": 28.0, "Ar": 39.9, "O2": 32.0, "CO2": 44.0}
MMHG = 13.5951 * 980.665 * 0.1  # dyn/cm2 per MILLIMETRE of Hg (rho g h, h = 0.1 cm)
ATM = 760.0 * MMHG

chain = []
chain.append(("eq. (3) coefficient = 1/sqrt(2 pi R)", 1 / np.sqrt(2 * np.pi * R), K3))
chain.append(("eq. (38) coefficient = 1/(eq. 3 coefficient)", 1 / K3, K38))
chain.append(("1 cu.mm at 20 C, 760 mm = pV/RT", ATM * 1e-3 / (R * 293.15), CUMM))
chain.append(("eq. (37) coefficient = N x 4.16e-8", NAV * CUMM, K37))
ch = pd.DataFrame(chain, columns=["identity", "derived", "printed"])
ch["rel_pct"] = (ch.printed / ch.derived - 1) * 100
print(ch.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
keep("chain_max_rel_pct", ch.rel_pct.abs().max())

# eq. (37) with the TEXT LAYER's exponent instead of the page's
keep("eq37_textlayer_factor", 25.2e16 / K37)
print(f"\nworst residual in the constant chain: {RES['chain_max_rel_pct']:.4f} %")
print(f"the text layer's 25.2e16 would multiply every N0 and beta on this page by "
      f"{RES['eq37_textlayer_factor']:.0f} and push most of them above the monolayer bound; "
      f"HOW MANY is not typed here, it is COUNTED in section 6.5 from the data")

# the bar unit, from Langmuir's own parenthesis
keep("hundred_bars_in_mmHg", 100 / MMHG)
print(f"\n100 bars = {RES['hundred_bars_in_mmHg']:.4f} mm Hg, printed as "
      f'"approximately 0.1 mm of mercury" (journal p. 1402)')

# %% [markdown]
# ## 4. The data
#
# Five CSVs, all transcribed from this paper alone. **No other page's dataset is
# loaded**, so none of the cross-page reconciliation obligations apply. **No curve
# is digitised and no figure is used for anything**, the paper's only figure with
# content is Fig. 1, a schematic checkerboard of surface atoms with no numbers on
# it. **No page image is committed anywhere.**
#
# | file | contents |
# |---|---|
# | `langmuir-1918-isotherms.csv` | every $(p, q_\mathrm{obs}, q_\mathrm{cal})$ of Tables II–XI and XIII–XVII, both bulbs |
# | `langmuir-1918-constants.csv` | the fitted $a$, $b$ per table joined to the derived columns of Tables XVIII–XX |
# | `langmuir-1918-lives.csv` | Table XXI (liquefied-gas lives) and Table XXII (the ratios) |
# | `langmuir-1918-platinum-N0.csv` | Table XXV, carried only to exercise eq. (37) on a third adsorbent |
# | `langmuir-1918-printed-claims.csv` | the 36 scalars that live in prose, equations or table headers |
#
# **Tier 5, not 6, and the distinction matters here.** The $q_\mathrm{obs}$ column
# *is* a measurement, pressures read on a McLeod gauge and gas volumes obtained by
# difference. Everything else in the paper ($a$, $b$, $b'$, $\sigma$, $N_0$,
# $\beta$, Tables XXI, XXII, XXV) is a constant Langmuir *derived* from those
# measurements, so reproducing one of those is reproduction, never validation, and
# the page says so at every such comparison.
#
# **Three things the source page says about its own rows, all of which affect us.**
#
# 1. *The mica results are uncorrected.* Journal page 1385: "The results therefore
#    obtained from bulb A are to be considered as being **10 to 30 % too high**",
#    because the blank bulb A′ adsorbed measurably and no subtraction was made.
#    Every mica $\beta$ below is therefore an over-estimate, which is the direction
#    that *helps* Langmuir's conclusion, noted at §6.5.
# 2. *One point is disowned.* Journal page 1384: the $q$ at $p = 12.8$ in Table IX
#    "is undoubtedly in error, as it cannot be made to fit into a smooth curve".
#    It is kept in the CSV, flagged here, and §6.6 reports the discrimination with
#    and without it.
# 3. *The blank bulb went negative.* The largest negative "adsorption" recorded is
#    $-0.7$ cu. mm. (Table XVII), which is the paper's own scale for its
#    measurement error and is used as such below.
#
# **A completeness check before anything else.** Tables XVIII–XX carry a "Number of
# obs." column. It must equal the number of Bulb A rows this page transcribed for
# that table. Twenty independent integers; if the transcription dropped or
# duplicated a row, this catches it.
#

# %%
iso = load_data("langmuir-1918-isotherms.csv", page=PAGE)
con = load_data("langmuir-1918-constants.csv", page=PAGE)
liv = load_data("langmuir-1918-lives.csv", page=PAGE).set_index("gas")
ptn = load_data("langmuir-1918-platinum-N0.csv", page=PAGE)
clm = load_data("langmuir-1918-printed-claims.csv", page=PAGE).set_index("key")

A = iso[iso.bulb == "A"].copy()          # the adsorbent bulb
B = iso[iso.bulb == "Aprime"].copy()     # the blank

counted = A.groupby(["table", "gas"]).size().rename("transcribed")
stated = con.groupby(["table", "gas"]).n_obs.max().rename("printed_n_obs")
chk = pd.concat([stated, counted], axis=1)
chk["ok"] = chk.printed_n_obs == chk.transcribed
print(chk.to_string())
if not chk.ok.all():
    raise ValueError("transcription does not match the printed 'Number of obs.' column")
keep("n_obs_mismatches", (~chk.ok).sum())
keep("n_isotherm_points_A", len(A))
keep("n_isotherm_points_blank", len(B))
print(f"\n{int(RES['n_obs_mismatches'])} mismatches over {len(chk)} tables; "
      f"{int(RES['n_isotherm_points_A'])} Bulb A points, "
      f"{int(RES['n_isotherm_points_blank'])} blank-bulb points")
keep("printed_claims_rows", len(clm))
keep("blank_largest_negative_glass", B[B.adsorbent == "glass"].q_obs_cu_mm.min())
keep("blank_largest_negative_all", B.q_obs_cu_mm.min())
print(f"largest negative blank reading in the GLASS series {RES['blank_largest_negative_glass']:+.1f} "
      f"cu.mm, which is the {clm.loc['blank_largest_negative','value']} cu. mm. of journal page 1389 "
      f"(that sentence is in the glass section). Over ALL blanks, including mica, it is "
      f"{RES['blank_largest_negative_all']:+.1f} cu.mm (Table V, argon).")


# %% [markdown]
# ## 5. PyMRM implementation
#
# **Said plainly: there is no transport on this page.** No grid, no time step, no
# boundary condition, no divergence operator. Eqs. (9), (26), (29), (31) and (33)
# are closed forms in one variable, and most of this notebook would run with pymrm
# uninstalled, exactly as on `A1.6`, `A1.1` and `J1.3`. Manufacturing a mesh here
# would be dishonest.
#
# `newton` and `NumJac` earn three narrow places, and each is chosen because a
# cheaper route would have been *wrong*, not merely slower.
#
# 1. **A route to the fitted constants that shares no algebra with the paper's.**
#    Langmuir got $a$ and $b$ by drawing a line on the $p/q$-versus-$p$ plot, which
#    is a *linear* least squares on a *transformed* variable, and that transform
#    reweights the residuals, hard, towards the low-pressure points. Fitting on $q$
#    itself is a different computation with a different answer, and the difference
#    is a measurement of what his graphical method cost. In both eq. (31) and
#    Freundlich's law the amplitude enters *linearly*, so it is eliminated exactly
#    and the stationarity condition in the one remaining parameter is **root-found**
#   , with `newton` on a `NumJac((1, 1))` Jacobian, and again with a bracketing
#    Brent solve. **Be precise about what that second solve does and does not
#    prove.** The two solvers are handed the *same* objective, the *same* profiled
#    amplitude and the *same* central-difference derivative; only the iteration
#    differs. `fit_two_route_max_rel` is therefore a **root-finder cross-check** , 
#    it catches a Newton that stopped early, a bad initial guess or a converged
#    non-root, and it says nothing whatever about the model or the algebra. The
#    independence in this item is the one in its heading: fitting on $q$ rather
#    than on Langmuir's $p/q$ linearisation, and §7.3 carries the break row that
#    refits his way and moves the discrimination metric.
# 2. **The rate equation integrated rather than quoted.** Eq. (14) is an ODE.
#    §6.7 marches it with a `newton`-solved implicit step and compares the result
#    against eq. (15) **as printed** and against eq. (15) **corrected**. Only one of
#    the two branches is reachable by the integration, and that is how the defect
#    is proved rather than asserted.
# 3. **Root-finds where a sweep would have been wrong.** The pressure at which
#    $d\ln q/d\ln p$ equals a stated exponent, and the relative pressure at which
#    Case VI departs from Case I by 1 %, are both root-found. §7 carries the break
#    row showing what sampling them instead would have cost.
#
# Everything else is `sympy` (the symbolic identities of §6.8) or arithmetic.
#

# %%
# The isotherm forms, the fits and the march live in model.py. The parameters of
# the numerical routes, each one a choice of this page:
FIT = dict(h_rel=1e-6,          # central-difference step of dS/dt, relative to max(|t|, 1)
           newton_tol=1e-12,    # newton step tolerance on the stationarity condition
           backward_tol=1e-8,   # newton backward-error tolerance
           xtol=1e-13)          # Brent's absolute tolerance on the same condition
FIT_TWO_ROUTE_MAX = 1e-9        # the bound on |t_newton/t_brent - 1|, set by the tolerances above
TWO_ROUTE = []                  # (label, |t_newton/t_brent - 1|) of the fits of 6.2 and 6.6
print("pymrm fits: newton on NumJac((1, 1)) and Brent on the same dS; parameters", FIT)


# %% [markdown]
# ## 6. Results
#
# ### 6.1 The fit, labelled as a fit
#
# $q_\mathrm{cal}$ in Tables II–XVII is eq. (31) (or eq. 33, or the two-term
# eq. 19) evaluated at the constants printed in that table's own header, which were
# obtained *from those same points*. Reproducing the column therefore tests the
# transcription and the arithmetic, nothing about the physics. It is reported here
# as a transcription check and **not** as evidence for the model.
#

# %%
def model_for(t, g):
    """Return a callable p -> q for whatever model the paper uses on that table."""
    k = con[(con.table == t) & (con.gas == g)]
    out = {}
    e31 = k[k.model == "eq31"]
    if len(e31):
        out["eq31"] = (lambda p, r=e31.iloc[0]: q_eq31(p, r.a_header, r.b_cu_mm))
    e33 = k[k.model == "eq33"]
    if len(e33) == 2:
        b1 = e33[e33.term == 1].b_cu_mm.iloc[0]
        r2 = e33[e33.term == 2].iloc[0]
        out["eq33"] = (lambda p, b1=b1, r2=r2: q_eq33(p, b1, r2.a_header, r2.b_cu_mm))
    e19 = k[k.model == "eq19"]
    if len(e19) == 2:
        pr = [(r.a_header, r.b_cu_mm) for _, r in e19.iterrows()]
        out["eq19"] = (lambda p, pr=pr: q_eq19(p, pr))
    return out


rows = []
for (t, g), grp in A.groupby(["table", "gas"]):
    for name, f in model_for(t, g).items():
        col = "q_cal2_cu_mm" if name == "eq19" else "q_cal_cu_mm"
        d = grp.dropna(subset=[col])
        if not len(d):
            continue
        r = f(d.p_bars.values) - d[col].values
        rows.append((t, g, name, len(d), np.abs(r).max()))
qc = pd.DataFrame(rows, columns=["table", "gas", "eq", "cells", "max_abs_cu_mm"])
print(qc.sort_values("max_abs_cu_mm", ascending=False).to_string(index=False))
keep("qcal_cells", qc.cells.sum())
keep("qcal_max_abs_cu_mm", qc.max_abs_cu_mm.max())
keep("qcal_median_abs_cu_mm", qc.max_abs_cu_mm.median())
print(f"\nGOODNESS OF FIT, NOT AGREEMENT: eq. (31)/(33)/(19) at the printed constants "
      f"reproduce all {int(RES['qcal_cells'])} printed q_cal cells to "
      f"{RES['qcal_max_abs_cu_mm']:.4f} cu.mm at worst "
      f"(median table {RES['qcal_median_abs_cu_mm']:.4f}). The constants were fitted to "
      f"these very points, so this tests transcription and arithmetic only.")

# %% [markdown]
# ### 6.2 What Langmuir's drawn line cost, and the two tables where it is over-determined
#
# On Tables IV, V and XVII there are exactly **two** observations, so a straight
# line through the $p/q$ plot is not fitted at all, it is *determined*. Those three
# tables let the printed constants be checked against the data with no freedom left,
# and one of them fails.
#

# %%
two = []
for t in ("IV", "V", "XVII"):
    d = A[A.table == t]
    p, q = d.p_bars.values, d.q_obs_cu_mm.values
    y = p / q
    slope = (y[0] - y[1]) / (p[0] - p[1])
    b_exact = 1.0 / slope
    a_exact = 1.0 / ((y[1] - slope * p[1]) * b_exact)
    k = con[(con.table == t) & (con.model == "eq31")].iloc[0]
    s = S_CM2[k.adsorbent]
    b_from_bprime = k.b_prime_cu_mm_per_m2 * s / 1e4
    b_at_printed_a = [qi * (1 + k.a_header * pi) / (k.a_header * pi) for pi, qi in zip(p, q)]
    two.append((t, k.gas, a_exact, b_exact, b_from_bprime,
                float(np.mean(b_at_printed_a)), k.a_header, k.b_cu_mm))
tw = pd.DataFrame(two, columns=["table", "gas", "a_exact", "b_exact", "b_from_bprime",
                                "b_at_printed_a", "a_printed", "b_printed"])
tw["b_printed_vs_exact_pct"] = (tw.b_printed / tw.b_exact - 1) * 100
tw["b_printed_vs_bprime_pct"] = (tw.b_printed / tw.b_from_bprime - 1) * 100
print(tw.to_string(index=False, float_format=lambda v: f"{v:.5g}"))

keep("t4_b_exact", tw[tw.table == "IV"].b_exact.iloc[0])
keep("t4_b_from_bprime", tw[tw.table == "IV"].b_from_bprime.iloc[0])
keep("t4_b_at_printed_a", tw[tw.table == "IV"].b_at_printed_a.iloc[0])
keep("t4_b_printed_dev_pct", tw[tw.table == "IV"].b_printed_vs_exact_pct.iloc[0])
keep("t4_two_route_b_rel_pct",
     abs(RES["t4_b_exact"] / RES["t4_b_from_bprime"] - 1) * 100)
keep("t5_t17_max_b_dev_pct", tw[tw.table != "IV"].b_printed_vs_exact_pct.abs().max())
print(f"\nTable IV: the two observations give b = {RES['t4_b_exact']:.4f}; the paper's own "
      f"b' = 100.0 requires b = {RES['t4_b_from_bprime']:.4f}. Those two independent routes "
      f"agree to {RES['t4_two_route_b_rel_pct']:.4f} %. The PRINTED b = 58.3 is "
      f"{RES['t4_b_printed_dev_pct']:+.2f} % away from both.")
print(f"Tables V and XVII, the same test: printed b is within "
      f"{RES['t5_t17_max_b_dev_pct']:.2f} % of the two-point solution.")

# how much the graphical fit cost, on the best-populated table
d7 = A[A.table == "VII"]; p7, q7 = d7.p_bars.values, d7.q_obs_cu_mm.values
k7 = con[(con.table == "VII") & (con.model == "eq31")].iloc[0]
rms_printed = rms(q_eq31(p7, k7.a_header, k7.b_cu_mm) - q7)
zLa, zLb, rms_opt = fit_langmuir(p7, q7, "Table VII", log=TWO_ROUTE, **FIT)
keep("t7_rms_printed", rms_printed)
keep("t7_rms_optimal", rms_opt)
keep("t7_graphical_penalty", rms_printed / rms_opt)
keep("t7_a_optimal", zLa); keep("t7_b_optimal", zLb)
print(f"\nTable VII: Langmuir's drawn line (a={k7.a_header}, b={k7.b_cu_mm}) gives RMS "
      f"{RES['t7_rms_printed']:.4f} cu.mm; the least-squares optimum "
      f"(a={RES['t7_a_optimal']:.5f}, b={RES['t7_b_optimal']:.4f}) gives "
      f"{RES['t7_rms_optimal']:.4f}. His graphical method cost "
      f"{(RES['t7_graphical_penalty']-1)*100:.2f} % in RMS.")

# %% [markdown]
# ### 6.3 The identity web: $a$ and $b$ to $b'$, $N_0$, $\sigma$ and $\beta$
#
# Tables XVIII–XX are five columns deep and every column after the second is
# computed from the one before it by a printed equation. Chained together they give
# about a hundred independent arithmetic constraints on the transcription, and they
# also locate the defects, because a defect breaks the chain at exactly one link.
#
# The first link is the sharpest. $b' = b/s$, so **every row implies the adsorbent
# area independently**, and the paper prints that area.
#

# %%
c = con.dropna(subset=["b_prime_cu_mm_per_m2"]).copy()
c["s_implied"] = c.b_cu_mm / c.b_prime_cu_mm_per_m2 * 1e4
c["s_printed"] = c.adsorbent.map(S_CM2)
c["s_dev_pct"] = (c.s_implied / c.s_printed - 1) * 100
print(c[["table", "gas", "adsorbent", "term", "b_cu_mm", "b_prime_cu_mm_per_m2",
         "s_implied", "s_dev_pct"]].to_string(index=False, float_format=lambda v: f"{v:.5g}"))

med = c.groupby("adsorbent").s_implied.median()
keep("mica_s_implied_median", med["mica"]); keep("glass_s_implied_median", med["glass"])
inl = c[c.s_dev_pct.abs() < 0.5]
keep("s_inlier_max_dev_pct", inl.s_dev_pct.abs().max())
keep("s_n_inliers", len(inl)); keep("s_n_rows", len(c))
out = c[c.s_dev_pct.abs() >= 0.5].sort_values("s_dev_pct", ascending=False)
keep("s_worst_dev_pct", c.s_dev_pct.max())
keep("s_second_worst_dev_pct", float(out.s_dev_pct.iloc[1]))
print(f"\n{int(RES['s_n_inliers'])} of {int(RES['s_n_rows'])} rows recover the printed area "
      f"to {RES['s_inlier_max_dev_pct']:.2f} % (median mica {RES['mica_s_implied_median']:.0f} "
      f"against the printed 5750, median glass {RES['glass_s_implied_median']:.0f} against 1966).")
print("The two that do not:")
print(out[["table", "gas", "b_cu_mm", "b_prime_cu_mm_per_m2", "s_implied", "s_dev_pct"]]
      .to_string(index=False, float_format=lambda v: f"{v:.5g}"))

# %%
MONO = {"N2": 0.66, "CH4": 0.63, "CO": 0.66, "Ar": 0.77, "O2": 0.77, "CO2": 0.61}

c["N0_from_eq37"] = K37 * c.b_prime_cu_mm_per_m2 / 1e4 / 1e15
c["beta_from_N0"] = c.N0_e15 / c.gas.map(MONO)
c["sigma_from_eq38"] = K38 * c.a_summary * np.sqrt(c.gas.map(MW) * c.T_K)

n0 = c.dropna(subset=["N0_e15"])
keep("N0_max_rel_pct", ((n0.N0_e15 / n0.N0_from_eq37 - 1).abs() * 100).max())
keep("N0_cells", len(n0))
betdf = c.dropna(subset=["beta"])
bd = (betdf.beta - betdf.beta_from_N0.round(2)).abs()
keep("beta_max_abs_dev", bd.max())
keep("beta_n_exact", int((bd < 5e-3).sum())); keep("beta_cells", len(betdf))
sg = c.dropna(subset=["sigma_s"])
srel = (sg.sigma_s / sg.sigma_from_eq38 - 1) * 100
# each cell's own printed granularity (10 ** the exponent of its last printed digit)
gran = np.array([10.0 ** Decimal(str(v)).normalize().as_tuple().exponent for v in sg.sigma_s])


def sigma_computed(k38=K38):
    """Eq. (38) on the cells with the coefficient ``k38``."""
    return k38 * sg.a_summary.values * np.sqrt(sg.gas.map(MW).values * sg.T_K.values)


def sigma_n_below(k38=K38):
    """Cells whose printed sigma lies below the computed one."""
    return int((sg.sigma_s.values / sigma_computed(k38) - 1 < 0).sum())


def sigma_reproduces(rule, k38=K38):
    """Cells reproduced exactly by cutting (``"trunc"``) or rounding (``"round"``) the
    computed value at the cell's own printed granularity."""
    cut = {"trunc": np.floor, "round": np.round}[rule]
    return int((cut(sigma_computed(k38) / gran) * gran == sg.sigma_s.values).sum())


keep("sigma_max_abs_rel_pct", srel.abs().max()); keep("sigma_mean_rel_pct", srel.mean())
keep("sigma_n_below", sigma_n_below()); keep("sigma_cells", len(sg))

print(f"eq. (37): N0 reproduced on {int(RES['N0_cells'])} rows to "
      f"{RES['N0_max_rel_pct']:.3f} % (that worst row is O2 at 155 K, where N0 is printed "
      f"to two decimals only).")
print(f"beta = N0 / monolayer: EXACT on {int(RES['beta_n_exact'])} of "
      f"{int(RES['beta_cells'])} rows. The exception is printed below.")
print(betdf.loc[bd >= 5e-3, ["table", "gas", "adsorbent", "N0_e15", "beta", "beta_from_N0"]]
      .to_string(index=False))
print(f"\neq. (38): sigma reproduced on {int(RES['sigma_cells'])} rows, worst "
      f"{RES['sigma_max_abs_rel_pct']:.2f} %, mean {RES['sigma_mean_rel_pct']:+.2f} %, and "
      f"{int(RES['sigma_n_below'])} of {int(RES['sigma_cells'])} printed values lie BELOW the "
      f"computed one.")
print(sg[["table", "gas", "T_K", "a_summary", "sigma_s", "sigma_from_eq38"]]
      .assign(rel_pct=srel.values).to_string(index=False, float_format=lambda v: f"{v:.6g}"))

# --- WHAT THAT IS AND IS NOT. It is a REPORTING CONVENTION, not an arithmetic error,
# and this page files it separately from the four genuine defects for that reason.
keep("sigma_trunc_reproduces", sigma_reproduces("trunc"))
keep("sigma_round_reproduces", sigma_reproduces("round"))
worst3 = sg.assign(rel_pct=srel.values).nsmallest(3, "rel_pct")
print(f"\nWHY THIS IS FILED AS A REPORTING CONVENTION AND NOT AS AN ERROR AGAINST LANGMUIR:")
print(f"  (a) the three worst rows are exactly the rows printed to ONE or TWO significant "
      f"figures: " + ", ".join(f"Table {r.table} {r.gas} {int(r.T_K)} K, printed "
                                 f"{r.sigma_s:.0f} ({r.rel_pct:+.2f} %)"
                                 for r in worst3.itertuples()) + ";")
print(f"  (b) truncating the computed value at each cell's OWN printed granularity reproduces "
      f"{int(RES['sigma_trunc_reproduces'])} of {int(RES['sigma_cells'])} cells against "
      f"{int(RES['sigma_round_reproduces'])} for rounding, so the convention is discarding "
      f"figures rather than mis-computing them.")
print(f"  The numeric statement stands: {int(RES['sigma_n_below'])} of "
      f"{int(RES['sigma_cells'])} below is not a coin toss, but it belongs beside "
      f"'he truncated', not beside Table IV's b or Table XX's beta, which are wrong.")

# %% [markdown]
# ### 6.4 Table XXII, and the mean the paper does not name
#
# Table XXII is stated to be the ratio of two lives: the adsorbed $\sigma$ of
# Tables XVIII–XX over the liquefied-gas $\sigma$ of Table XXI. For eleven of the
# fourteen cells that is one division. For three of them, nitrogen and methane on
# mica at 90 K, carbon dioxide on mica at 155 K, the gas has **two or three**
# adsorbed $\sigma$, and **the paper does not say how they were combined**. The
# cell below measures it instead of assuming.
#

# %%
sets = {"mica_90K": ("mica", 90, "sigma_liq_90K_s", "ratio_mica_90K"),
        "mica_155K": ("mica", 155, "sigma_liq_155K_s", "ratio_mica_155K"),
        "glass_90K": ("glass", 90, "sigma_liq_90K_s", "ratio_glass_90K")}
rows = []
for name, (ads, T, lc, rc) in sets.items():
    for g in MONO:
        pr = liv.loc[g, rc]
        if not np.isfinite(pr):
            continue
        sl = liv.loc[g, lc]
        s = c[(c.adsorbent == ads) & (c.T_K == T) & (c.gas == g)].sigma_s.dropna()
        s = s[s > 2e4] if len(s) > 1 else s      # the tiny 2nd-site sigma of Table IX is a site, not a run
        am = float(s.mean()) / sl
        gm = float(np.exp(np.log(s).mean())) / sl
        rows.append((name, g, len(s), pr, am, (pr / am - 1) * 100, gm, (pr / gm - 1) * 100))
t22 = pd.DataFrame(rows, columns=["set", "gas", "n_sigma", "printed",
                                  "arith", "arith_pct", "geom", "geom_pct"])
print(t22.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
keep("t22_cells", len(t22))
keep("t22_single_max_pct", t22[t22.n_sigma == 1].arith_pct.abs().max())
multi = t22[t22.n_sigma > 1]
keep("t22_multi_arith_max_pct", multi.arith_pct.abs().max())
keep("t22_multi_geom_max_pct", multi.geom_pct.abs().max())
print(f"\n{int((t22.n_sigma==1).sum())} unambiguous cells reproduce to "
      f"{RES['t22_single_max_pct']:.2f} %. On the {len(multi)} ambiguous ones the ARITHMETIC "
      f"mean is out by up to {RES['t22_multi_arith_max_pct']:.2f} % and the GEOMETRIC mean by "
      f"{RES['t22_multi_geom_max_pct']:.2f} %; the page states the geometric mean as a "
      f"RECONSTRUCTION of an unprinted basis, not as the paper's rule.")

# eq. (37) on a third adsorbent: Table XXV
ptn["N0_from_eq37"] = K37 * ptn.q_cu_mm / S_CM2["platinum"] / 1e15
keep("pt_N0_max_abs_dev", (ptn.N0_e15 - ptn.N0_from_eq37.round(2)).abs().max())
keep("pt_cells", len(ptn))
print(f"\neq. (37) on platinum (s = 312 sq.cm., against 5750 and 1966): all "
      f"{int(RES['pt_cells'])} cells of Table XXV reproduce exactly "
      f"(max |printed - rounded| = {RES['pt_N0_max_abs_dev']:.4f}).")
keep("pt_foil_area", 2 * float(clm.loc["platinum_foil_length", "value"])
     * float(clm.loc["platinum_foil_width", "value"]))
print(f"and the platinum area itself: 2 x 15.3 x 10.2 = {RES['pt_foil_area']:.2f}, "
      f"printed 312 sq. cm. (both faces of the foil).")

# %% [markdown]
# ### 6.5 The test that is not a fit: $\beta < 1$
#
# This is the paper's actual evidence, and it is a genuine out-of-sample bound.
# Nothing in the chain $b \to b' \to N_0 \to \beta$ is adjustable: $b$ came from the
# isotherm, $s$ from a ruler, eq. (37) from Avogadro's number, and the denominator
# from **liquid densities**, which contain no adsorption at all. Journal page 1391
# states the conclusion:
#
# > "The observed values of $N_0$ are all less than these calculated results."
#
# It could have failed on any of twenty entries and did not. But the honest question
# is *by how much*, and the answer is not comfortable.
#
# Only nitrogen's molecular volume (35.5 cu. cm.) is printed, so the chain
# $35.5 \to 3.88\times10^{-8}\ \mathrm{cm} \to 0.66\times10^{15}$ can be checked
# end to end, and the other four molecular volumes are **recovered** by inverting
# the printed counts. They are recoveries and are labelled as such; nothing on this
# page depends on them being right, because $\beta$ uses the printed counts.
#

# %%
# the printed chain, end to end, for the one gas whose molecular volume is printed
Vm_N2 = float(clm.loc["liquid_N2_molecular_volume", "value"])
d_calc = (Vm_N2 / NAV) ** (1 / 3)
n_calc = 1.0 / d_calc ** 2
keep("mono_d_rel_pct", (d_calc / float(clm.loc["N2_molecular_diameter", "value"]) - 1) * 100)
keep("mono_n_rel_pct", (n_calc / float(clm.loc["monolayer_N2", "value"]) - 1) * 100)
print(f"35.5 cu.cm -> d = {d_calc:.4e} cm (printed 3.88e-8, {RES['mono_d_rel_pct']:+.3f} %) "
      f"-> {n_calc:.4e} per sq.cm (printed 0.66e15, {RES['mono_n_rel_pct']:+.3f} %)")
print("\nRECOVERED (not printed) molecular volumes behind the other four counts:")
for g, n in MONO.items():
    print(f"  {g:<4s} {(1/np.sqrt(n*1e15))**3 * NAV:6.2f} cu.cm per gram molecule")

# the bound itself, summed WITHIN a table only (two-term fits are one surface)
g = (c.groupby(["adsorbent", "T_K", "gas", "table"])[["b_prime_cu_mm_per_m2", "N0_e15", "beta"]]
       .sum().reset_index())
g["beta_recomputed"] = g.N0_e15 / g.gas.map(MONO)
print("\n" + g.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
keep("beta_entries", len(g))
keep("beta_max", g.beta_recomputed.max())
keep("beta_headroom", 1.0 / g.beta_recomputed.max())
keep("beta_n_above_one", int((g.beta_recomputed >= 1).sum()))
tight = g.loc[g.beta_recomputed.idxmax()]
print(f"\nTHE TEST: {int(RES['beta_n_above_one'])} of {int(RES['beta_entries'])} entries exceed "
      f"unity. Worst is {tight.gas} on {tight.adsorbent}, Table {tight.table}, at "
      f"beta = {RES['beta_max']:.4f}, only {(RES['beta_headroom']-1)*100:.1f} % of headroom.")
# the text-layer trap of section 3, COUNTED here rather than asserted there
def textlayer_count(factor):
    """Entries above the monolayer bound when eq. (37)'s coefficient is read ``factor``
    times too large."""
    return int((g.beta_recomputed * factor >= 1).sum())


keep("beta_n_above_one_textlayer", textlayer_count(RES["eq37_textlayer_factor"]))
below_tl = g.loc[g.beta_recomputed * RES["eq37_textlayer_factor"] < 1]
keep("beta_textlayer_min_margin",
     float(1.0 / (below_tl.beta_recomputed * RES["eq37_textlayer_factor"]).max()))
print(f"\nTHE TEXT-LAYER TRAP OF SECTION 3, COUNTED: reading eq. (37) as 25.2e16 multiplies every "
      f"beta by {RES['eq37_textlayer_factor']:.0f}, and {int(RES['beta_n_above_one_textlayer'])} of "
      f"{int(RES['beta_entries'])} entries would then exceed unity, Langmuir's conclusion reversed. "
      f"Only {int(RES['beta_entries']) - int(RES['beta_n_above_one_textlayer'])} would survive: "
      + ", ".join(f"{r.gas} {r.adsorbent} {int(r.T_K)} K (beta = {r.beta_recomputed:.4f})"
                  for r in below_tl.itertuples())
      + f". The nearest of those is still a factor {RES['beta_textlayer_min_margin']:.1f} below unity, "
        f"so no rounding convention can move the count.")
lo, hi = [float(clm.loc[k, "value"]) for k in ("mica_blank_low_pct", "mica_blank_high_pct")]
keep("beta_max_corrected_hi", RES["beta_max"] * (1 - hi / 100))
keep("beta_max_corrected_lo", RES["beta_max"] * (1 - lo / 100))
print(f"AND THE CORRECTION GOES THE RIGHT WAY: the paper says the mica figures are "
      f"{lo:.0f} to {hi:.0f} % too high, which puts that worst beta between "
      f"{RES['beta_max_corrected_hi']:.3f} and {RES['beta_max_corrected_lo']:.3f}. "
      f"The margin is real but it is thin, and it is thin on the paper's OWN numbers.")

# the ordering claim of journal page 1391, and the one column that breaks it
print("\nJournal p. 1391: 'if we arrange the gases in order ... as indicated by b\', N0 or "
      "beta we find that this order is the same in all three sets of experiments'")
order_rows = []
for (ads, T), grp in g.groupby(["adsorbent", "T_K"]):
    for col, lab in (("b_prime_cu_mm_per_m2", "b'"), ("N0_e15", "N0"),
                     ("beta", "beta as printed"), ("beta_recomputed", "beta recomputed")):
        o = " > ".join(grp.groupby("gas")[col].max().sort_values(ascending=False).index)
        order_rows.append((f"{ads} {T}K", lab, o))
od = pd.DataFrame(order_rows, columns=["set", "by", "order"])
print(od.to_string(index=False))
ref = {s: od[(od.set == s) & (od.by == "b'")].order.iloc[0] for s in od.set.unique()}
keep("order_breaks_printed_beta",
     int(sum(od[(od.by == "beta as printed")].apply(lambda r: r.order != ref[r.set], axis=1))))
keep("order_breaks_recomputed_beta",
     int(sum(od[(od.by == "beta recomputed")].apply(lambda r: r.order != ref[r.set], axis=1))))
print(f"\nThe order agrees with b' and N0 in all three sets. Using beta AS PRINTED it breaks in "
      f"{int(RES['order_breaks_printed_beta'])} set; using beta RECOMPUTED from the paper's own "
      f"N0 and its own monolayer counts it breaks in "
      f"{int(RES['order_breaks_recomputed_beta'])}. That is the second, independent line of "
      f"evidence on the Table XX defect of section 6.9.")

# %% [markdown]
# ### 6.6 The discrimination: Langmuir against Freundlich, refitted at their best
#
# This is the only place in the paper where two theories are made to disagree about
# something measurable, and it is a **fair** contest: two parameters each, the same
# eleven points, both fitted to the data they are judged on.
#
# Langmuir's own version is not quite fair, in his favour and against it at once:
# his Langmuir line is not the least-squares optimum (§6.2) and neither is his
# $q_F = 8.4\,p^{0.417}$. Both are refitted here, with `newton` on the nonlinear
# normal equations, so that Freundlich is beaten at its best rather than at
# Langmuir's rendering of it.
#
# Two null baselines are printed beside every ratio: a **zero-parameter** model
# ($q = \bar q$) and a **one-parameter** proportional model ($q = kp$, which is what
# eq. (12) predicts at vanishing pressure). Without them, "the RMS is 0.35 cu. mm."
# means nothing.
#

# %%
# the paper's own Freundlich column, against the formula it is labelled with
qF_formula = q_freundlich(p7, 8.4, 0.417)
keep("qF_column_max_dev_pct",
     float((np.abs(qF_formula - d7.q_F_cu_mm.values) / d7.q_F_cu_mm.values).max() * 100))
print(f"The printed q_F column departs from an exact 8.4 p^0.417 by up to "
      f"{RES['qF_column_max_dev_pct']:.2f} %: it was read off the drawn log-log line, "
      f"not computed from the formula. Both are carried; the contest below uses neither.")

def discrimination(iso_a, langmuir="fit", freundlich="fit", fit=FIT, log=None):
    """Best Langmuir against best Freundlich, with the two nulls, on every table
    with four or more points. langmuir="linearised" fits eq. (31) by the p/q line
    instead; freundlich="printed" uses the paper's own 8.4 p^0.417 everywhere."""
    rows = []
    for (t, gname), grp in iso_a.groupby(["table", "gas"]):
        if len(grp) < 4:
            continue
        p, q = grp.p_bars.values, grp.q_obs_cu_mm.values
        if langmuir == "fit":
            rL = fit_langmuir(p, q, f"{t}/{gname}", log=log, **fit)[2]
        else:
            sl, ic = np.polyfit(p, p / q, 1)
            rL = rms(q_eq31(p, sl / ic, 1.0 / sl) - q)
        if freundlich == "fit":
            rF = fit_freundlich(p, q, f"{t}/{gname}", log=log, **fit)[2]
        else:
            rF = rms(q_freundlich(p, 8.4, 0.417) - q)
        r_mean = rms(q - q.mean())
        r_prop = rms((np.dot(p, q) / np.dot(p, p)) * p - q)
        rows.append((t, gname, len(p), p.max() / p.min(), rL, rF, rF / rL,
                     r_mean / rL, r_prop / rL))
    return pd.DataFrame(rows, columns=["table", "gas", "n", "p_span", "rms_Langmuir",
                                       "rms_Freundlich", "F_over_L", "null_mean", "null_kp"])


dsc = discrimination(A, log=TWO_ROUTE)
print("\n" + dsc.to_string(index=False, float_format=lambda v: f"{v:.5g}"))

keep("fit_two_route_max_rel", max(v for _, v in TWO_ROUTE))
keep("fit_two_route_count", len(TWO_ROUTE))
print(f"\nevery fit above was root-found TWICE: pymrm's newton and a bracketing Brent solve. "
      f"Over {int(RES['fit_two_route_count'])} fits the worst disagreement is below "
      f"{FIT_TWO_ROUTE_MAX:.0e} relative (here {RES['fit_two_route_max_rel']:.1e}). READ THAT NARROWLY: both solvers are given "
      f"the SAME objective S, the SAME profiled amplitude and the SAME central-difference dS, "
      f"and differ only in how they iterate, so this is a ROOT-FINDER CROSS-CHECK and not an "
      f"independent computation. What IS independent here is the objective itself: these fits "
      f"are on q, Langmuir's were on the p/q linearisation, and section 7.3 breaks that.")

keep("disc_tables", len(dsc))
keep("disc_n_langmuir_wins", int((dsc.F_over_L > 1).sum()))
keep("disc_median_ratio", dsc.F_over_L.median())
keep("t7_F_over_L", dsc[dsc.table == "VII"].F_over_L.iloc[0])
keep("t7_null_mean", dsc[dsc.table == "VII"].null_mean.iloc[0])
keep("t7_null_kp", dsc[dsc.table == "VII"].null_kp.iloc[0])
keep("disc_min_null_mean", dsc.null_mean.min())
keep("disc_min_null_kp", dsc.null_kp.min())
worst = dsc.loc[dsc.F_over_L.idxmin()]
keep("disc_worst_ratio", worst.F_over_L)
print(f"\nDISCRIMINATION, NOT FIT: Langmuir beats the best two-parameter Freundlich on "
      f"{int(RES['disc_n_langmuir_wins'])} of {int(RES['disc_tables'])} tables; median ratio "
      f"{RES['disc_median_ratio']:.3f}; on Table VII, the best-populated, "
      f"{RES['t7_F_over_L']:.4f}.")
print(f"NULL BASELINES on Table VII: a constant is {RES['t7_null_mean']:.1f} x worse than "
      f"Langmuir and a proportional law {RES['t7_null_kp']:.1f} x; across all tables the "
      f"WEAKEST nulls are {RES['disc_min_null_mean']:.1f} x and {RES['disc_min_null_kp']:.1f} x.")
print(f"THE ONE EXCEPTION is Table {worst.table} ({worst.gas} on glass), where Freundlich is "
      f"{1/RES['disc_worst_ratio']:.2f} x BETTER, and that is exactly the table journal "
      f"page 1388 says 'do not give a straight line when p/q_obs is plotted against p and "
      f"therefore do not satisfy Equation 31'. Langmuir replaces it with eq. (33):")
d13 = A[A.table == "XIII"]
f33 = model_for("XIII", "CO")["eq33"]
keep("t13_rms_eq33", rms(f33(d13.p_bars.values) - d13.q_obs_cu_mm.values))
keep("t13_rms_freundlich", dsc[dsc.table == "XIII"].rms_Freundlich.iloc[0])
print(f"  eq. (33) at the printed constants: RMS {RES['t13_rms_eq33']:.4f} cu.mm, against "
      f"{RES['t13_rms_freundlich']:.4f} for the best Freundlich, but eq. (33) has THREE "
      f"parameters to Freundlich's two, and the page does not count it as a win.")

# the disowned point
d9 = A[A.table == "IX"]
m = d9.p_bars.values != 12.8
_, _, r9_all = fit_langmuir(d9.p_bars.values, d9.q_obs_cu_mm.values, "IX all", **FIT)
_, _, r9_cut = fit_langmuir(d9.p_bars.values[m], d9.q_obs_cu_mm.values[m], "IX cut", **FIT)
_, _, rF9_all = fit_freundlich(d9.p_bars.values, d9.q_obs_cu_mm.values, "IX all", **FIT)
_, _, rF9_cut = fit_freundlich(d9.p_bars.values[m], d9.q_obs_cu_mm.values[m], "IX cut", **FIT)
keep("t9_ratio_all", rF9_all / r9_all); keep("t9_ratio_cut", rF9_cut / r9_cut)
print(f"\nTable IX with and without the point the paper disowns (p = 12.8): the "
      f"Freundlich/Langmuir ratio moves {RES['t9_ratio_all']:.3f} -> {RES['t9_ratio_cut']:.3f}. "
      f"The conclusion does not depend on it.")

# %%
fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))
pp = np.logspace(np.log10(p7.min() * 0.8), np.log10(p7.max() * 1.25), 300)
aFo, nFo, _ = fit_freundlich(p7, q7, "Table VII figure", **FIT)
ax[0].plot(p7, q7, "o", ms=6, mfc="none", color="k", label="Table VII, $q_{obs}$")
ax[0].plot(pp, q_eq31(pp, zLa, zLb), "-", lw=1.8, color="#1f6feb",
           label=f"Langmuir, best fit (RMS {rms_opt:.3f})")
ax[0].plot(pp, q_freundlich(pp, aFo, nFo), "--", lw=1.8, color="#d1242f",
           label=f"Freundlich, best fit (RMS {dsc[dsc.table=='VII'].rms_Freundlich.iloc[0]:.3f})")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("p  (bars = dyn cm$^{-2}$)"); ax[0].set_ylabel("q  (cu. mm.)")
ax[0].set_title("Nitrogen on mica, 90 K, both models at their best")
ax[0].legend(fontsize=8, frameon=False)

ax[1].axhline(0, color="k", lw=0.7)
ax[1].plot(p7, q_eq31(p7, zLa, zLb) - q7, "o-", ms=5, lw=1.2, color="#1f6feb", label="Langmuir")
ax[1].plot(p7, q_freundlich(p7, aFo, nFo) - q7, "s--", ms=5, lw=1.2, color="#d1242f",
           label="Freundlich")
ax[1].fill_between([p7.min() * 0.8, p7.max() * 1.25], -0.7, 0.7, color="0.85", zorder=0,
                   label="paper's own blank-bulb scale, $\pm$0.7 cu. mm.")
ax[1].set_xscale("log"); ax[1].set_xlabel("p  (bars)"); ax[1].set_ylabel("$q_{cal}-q_{obs}$  (cu. mm.)")
ax[1].set_title("Residuals: the shape of the failure, not its size")
ax[1].legend(fontsize=8, frameon=False)
fig.tight_layout(); plt.show()
print("Freundlich's residuals are not merely larger, they SWEEP: positive at both ends and "
      "negative in the middle, which is the signature of a power law forced through a curve "
      "that saturates. That is Langmuir's own argument on journal page 1375, made visible.")


# %% [markdown]
# ### 6.7 Eq. (15) does not survive contact with eq. (16)
#
# Langmuir integrates eq. (14) and prints
#
# $$t \;=\; \frac{N_0}{N\nu_1(1+\sigma_1\mu)}\,\ln\frac{\theta'}{\theta_1-\theta'}
# \tag{15, as printed}$$
#
# and then, two lines below, obtains the half-period "by placing
# $\theta' = \tfrac12\theta_1$":
#
# $$t_{1/2} \;=\; \frac{\tau\ln 2}{\alpha(1+\sigma_1\mu)}.
# \tag{16}$$
#
# **Put $\theta' = \tfrac12\theta_1$ into eq. (15) as printed and the logarithm is
# $\ln 1 = 0$.** Eq. (15) as printed gives $t_{1/2} = 0$, not eq. (16). It also
# gives $t = -\infty$ at $\theta' = 0$, where the integration starts.
#
# Proved the F2.3 way, pin what is *not* free first. Eq. (14) is an ODE and is not
# free. Eq. (16) is printed and is not free. Eq. (10), $\tau = N_0\sigma_1/N$, and
# eq. (7), $\sigma_1 = \alpha/\nu_1$, are printed and are not free. The only free
# thing left is the numerator inside eq. (15)'s logarithm, and putting $\theta_1$
# there, a one-character change, the loss of a subscript, makes the integration
# correct, makes $t(0) = 0$, and returns eq. (16) exactly.
#
# The cell below does not argue: it **integrates eq. (14)** with `newton` and
# compares the march against both branches. **Reported, not repaired**, and the
# correction is labelled an inference.
#
# **Which of the lines it prints are evidence, and which are not.** The cell prints
# four values of $t_{1/2}$ and they are not of equal standing, so it says which is
# which rather than leaving a reader to assume:
#
# | line | what it is |
# |---|---|
# | eq. (15) as printed | $\ln 1 = 0$, the defect itself |
# | eq. (15) corrected | the inference under test |
# | eq. (14) **marched**, then root-found on the trajectory | **the independent numerical route** |
# | eq. (14) solved in closed form, then inverted | an **algebraic identity**, evidence for nothing |
#
# Only the third is independent. The fourth inverts $\theta_1(1-e^{-kt})$ at half
# coverage, which returns $\ln 2/k$, which eqs. (7) and (10) turn into eq. (16)
# *identically*, so its residual, at most 1e-14, is float round-off and **could not have
# come out otherwise**. It is kept because it closes the algebra of eqs. (7), (10),
# (14) and (16), and it is labelled so that it is never read as a fourth
# confirmation. The march is the confirmation: it reaches eq. (16) with a
# *discretisation* error that halves as $\Delta t$ halves, its Richardson value lies
# within the bound printed below of eq. (16), and it cannot reach eq. (15) as printed
# at any step size.
#

# %%
# a concrete case built only from printed numbers: nitrogen on mica, Table VII
kk = con[(con.table == "VII") & (con.model == "eq31")].iloc[0]
sig_mu = kk.a_header * 34.0                       # sigma_1 mu = a p, eq. (34), at the top pressure
theta1 = sig_mu / (1 + sig_mu)
# Work in units of tau/alpha: eq. (16) is t_half = tau ln2 / (alpha (1 + sigma_1 mu)), and
# eqs. (7) and (10) make N0/(N nu1) = tau/alpha exactly, so setting tau = alpha = 1 loses
# nothing and leaves the ONE dimensionless group the comparison depends on, sigma_1 mu.
k_rate = (1.0 + sig_mu)                           # (N/N0)(nu1 + alpha mu) in those units
t_half_eq16 = np.log(2.0) / (1.0 + sig_mu)        # eq. (16) as printed, alpha = 1, tau = 1

NSTEPS = (200, 400, 800, 1600)    # nested ladder: every coarse step time is a fine one
T_END = 6.0 / k_rate              # six relaxation times
STEP_TOL = 1e-14                  # newton step tolerance of each implicit step
ROUNDOFF_MIN = 1e3                # finest successive difference over the march's floor, at least
RICHARDSON_MAX = 5e-5             # Richardson half-coverage time against eq. (16), at most


def rate_study(k_march=k_rate, scheme="implicit", nsteps=NSTEPS, step_tol=STEP_TOL):
    """March eq. (14) on the ladder `nsteps` and read everything section 6.7 states.

    k_march is the rate constant of eq. (14), used by the march and by eq. (14)'s
    closed-form solution; eq. (16) stays as printed. Orders come from successive
    differences (trajectory at the coarsest step times, and the half-coverage time),
    the round-off companion divides the finest difference by the march's floor
    nsteps[-1] * max(step_tol, eps) * theta1, and the Richardson value of the
    half-coverage time is the extrapolated reference."""
    marches = {n: march_eq14(theta1, k_march, T_END, n, scheme=scheme, newton_tol=step_tol)
               for n in nsteps}
    t, th = marches[nsteps[-1]]
    ok = th < theta1 * (1 - 1e-9)
    err = np.abs(t[ok][1:] - t_eq15_corrected(th[ok][1:], theta1, k_march)).max()
    traj = observed_orders(lambda n: marches[n][1][::n // nsteps[0]], nsteps)
    half = observed_orders(lambda n: t_half_of_march(*marches[n], theta1), nsteps)
    finest = np.max(np.abs(traj["values"][-1] - traj["values"][-2]))
    floor = nsteps[-1] * max(step_tol, np.finfo(float).eps) * theta1
    t_cf = unique_root(lambda tt: theta1 * (1 - np.exp(-k_march * tt)) - 0.5 * theta1,
                       1e-12, 100.0 / k_march, log=True)
    return {"eq14_vs_eq15corrected_max_abs": err,
            "eq14_time_order": min(traj["orders"], key=lambda o: -abs(o - 1.0)),
            "eq14_ladder_over_roundoff": finest / floor,
            "t_half_march_vs_eq16_rel": abs(half["values"][-1] / t_half_eq16 - 1),
            "t_half_march_order": min(half["orders"], key=lambda o: -abs(o - 1.0)),
            "t_half_richardson_vs_eq16_rel": abs(half["extrapolated"] / t_half_eq16 - 1),
            "t_half_richardson_p1_vs_eq16_rel": abs(
                (half["values"][-1] + (half["values"][-1] - half["values"][-2])
                 / (nsteps[-1] / nsteps[-2] - 1)) / t_half_eq16 - 1),
            "t_half_closed_form_vs_eq16_rel": abs(t_cf / t_half_eq16 - 1),
            "_traj": traj, "_half": half, "_t_half_closed_form": t_cf}


RATE = rate_study()
for k_, v_ in RATE.items():
    if not k_.startswith("_"):
        keep(k_, v_)
keep("eq14_march_steps", NSTEPS[-1])
keep("t_half_march_rootfind", RATE["_half"]["values"][-1])
keep("t_half_closed_form_rootfind", RATE["_t_half_closed_form"])
print(f"  implicit Euler on {NSTEPS} steps over six relaxation times:")
print(f"  max |t_march - t_eq15_corrected| at {NSTEPS[-1]} steps = "
      f"{RES['eq14_vs_eq15corrected_max_abs']:.3e}")
print(f"  observed orders of the trajectory from successive differences: "
      + ", ".join(f"{o:.4f}" for o in RATE["_traj"]["orders"])
      + f" (implicit Euler is first order; dt is the one axis that carries error here); "
      f"the gap from 1 halves with dt, as the next-order term of a first-order method does")
print(f"  round-off companion: the finest successive difference is "
      f"{RES['eq14_ladder_over_roundoff']:.2e} times the march's floor "
      f"(at least {ROUNDOFF_MIN:.0e} required)")

# now both branches at theta' = theta1/2
t_half_printed = t_eq15_printed(0.5 * theta1, theta1, k_rate)
t_half_corrected = t_eq15_corrected(0.5 * theta1, theta1, k_rate)
keep("t_half_eq15_printed", t_half_printed)
keep("t_half_eq15_corrected", t_half_corrected)
keep("t_half_eq16", t_half_eq16)
keep("t_half_corrected_vs_eq16_rel", abs(t_half_corrected / t_half_eq16 - 1))
print(f"\n  eq. (15) AS PRINTED at theta'=theta1/2 : t_half = {RES['t_half_eq15_printed']:.6g}")
print(f"  eq. (15) CORRECTED                    : t_half = {RES['t_half_eq15_corrected']:.6g}")
print(f"  eq. (16) as printed                   : t_half = {RES['t_half_eq16']:.6g}")
print(f"  eq. (14) MARCHED, then root-found     : t_half = {RES['t_half_march_rootfind']:.6g}"
      f"   ({int(RES['eq14_march_steps'])} steps)")
print(f"  eq. (14) solved in CLOSED FORM, then inverted: t_half = "
      f"{RES['t_half_closed_form_rootfind']:.6g}")
print(f"\n  corrected  vs eq. (16): {RES['t_half_corrected_vs_eq16_rel']:.3e} relative")
print(f"  MARCH      vs eq. (16): {RES['t_half_march_vs_eq16_rel']:.3e} relative, and it is a "
      f"DISCRETISATION error, not a disagreement; refined:")
for n, tm in zip(NSTEPS, RATE["_half"]["values"]):
    print(f"      n_steps={n:5d}  t_half = {tm:.8f}  rel = {abs(tm / t_half_eq16 - 1):.4e}")
print(f"      observed orders in dt from successive differences: "
      + ", ".join(f"{o:.4f}" for o in RATE["_half"]["orders"])
      + f"; the Richardson value is {RES['t_half_richardson_vs_eq16_rel']:.1e} from eq. (16), "
      f"so the march converges TO eq. (16), not merely converges. That gap is dominated by "
      f"the error of the estimated order: extrapolating with the theory order p = 1 instead "
      f"gives {RES['t_half_richardson_p1_vs_eq16_rel']:.1e}, so the bound covers the order "
      f"estimate. It cannot reach eq. (15) as printed at any dt.")
print(f"  CLOSED FORM vs eq. (16): (here {RES['t_half_closed_form_vs_eq16_rel']:.1e}), and "
      f"this line is an ALGEBRAIC IDENTITY, NOT a fourth route. Inverting theta1(1-exp(-kt)) "
      f"at half coverage gives ln2/k, which eqs. (7) and (10) turn into eq. (16) exactly, so "
      f"that number is float round-off and could not have come out otherwise. It is "
      f"shown because it closes the algebra, and it is cited as evidence for nothing.")
print(f"  printed    vs eq. (16): eq. (15) as printed gives EXACTLY ZERO, which is not a small "
      f"discrepancy but a different function.")
print("\nREPORTED, NOT REPAIRED. The paper's eq. (15) is quoted verbatim in section 2. The "
      "substitution theta' -> theta_1 in the numerator is an INFERENCE, and it is the unique "
      "one-character change that makes eq. (14), eq. (15) and eq. (16) mutually consistent.")
print("Note what is NOT claimed: eq. (16) is right, so nothing the paper concludes is affected. "
      "The paper never measures a rate, and no number anywhere in it depends on eq. (15).")

# %% [markdown]
# ### 6.8 Case VI is the BET isotherm, twenty years early
#
# This is the page's own result and it is not in the paper.
#
# **The collapse $\sigma_3=\sigma_4=\dots=\sigma_2$ is Langmuir's, not this page's,
# and so is its consequence.** He prints both immediately under eq. (30), journal
# page 1375, read here on a 300 ppi crop and quoted verbatim:
#
# > "If $\sigma_1$, $\sigma_2$, $\sigma_3$, etc., are equal, all the coefficients in
# > (29) after *a* are zero. **If $\sigma_1$ and $\sigma_2$ are different, but all
# > subsequent values of $\sigma$ (*i. e.*, $\sigma_3$, $\sigma_4$, etc.) are equal
# > to $\sigma_2$, then all the coefficients in (29) after *b* are zero.** Equation
# > 29 thus takes a very simple form and shows that at very low pressures $\eta$ is
# > proportional to $\mu$, but at pressures close to saturation $\eta$ begins to
# > increase rapidly and becomes infinite when saturation is reached."
#
# That is the whole licence, in his own words: he states the assumption **and** he
# states that under it the series *terminates at b*, which is exactly what makes
# truncating eq. (29) at $b$ legitimate rather than convenient. He also states
# where the resulting isotherm diverges. The page's contribution is not the
# assumption; it is noticing what the terminated series *is*.
#
# (Journal page 1374 carries the physical motivation for the same collapse, *"There
# may be a small difference between $\sigma_2$ and $\sigma_3$, but as the number of
# layers increases still further the values of $\sigma$ should remain practically
# constant"*, verified verbatim on its own crop. Read alone that sentence only
# *permits* $\sigma_3\neq\sigma_2$; the p. 1375 statement is the stronger one and is
# the one this section rests on.)
#
# Set $\sigma_3 = \sigma_4 = \dots = \sigma_2$ in eq. (30),
# put it into eq. (29), and change variables to the two the paper itself supplies , 
# $x = \sigma_2\mu$, which is $p/p_0$ by eq. (39), and $c = \sigma_1/\sigma_2$:
#
# $$\frac{N\eta}{N_0}
# \;=\; \frac{\sigma_1}{1/\mu + a + b\mu}
# \;=\; \frac{c\,x}{(1-x)\,\bigl(1+(c-1)x\bigr)}.$$
#
# That is the BET equation. `J1.3` (published today, `pages/J1.3-bet-multimolecular-adsorption/`)
# carries the 1938 paper and reports the same form as its eqs. (26)/(28); the
# identity below is verified against *Langmuir's* coefficients, symbolically, and
# `J1.3`'s data are not loaded and none of its numbers is retyped.
#
# Five symbolic checks in one cell, each of which can fail:
#
# 1. eq. (26) really does follow from eqs. (22)–(25);
# 2. eq. (30)'s $a$, $b$ **and** $c$ really are the coefficients of eq. (29), three
#    printed expressions, none of which the paper derives;
# 3. both of the paper's collapse claims about eq. (30);
# 4. the BET identity itself;
# 5. **and the same ladder summed to infinity, in closed form, with no truncation
#    anywhere.** Checks 1–4 run through eq. (29) truncated at $b$, which is exactly
#    what Langmuir's p.-1375 statement licenses, but it is still a truncation. Under
#    the collapse the occupancies of eqs. (22)/(23) are $\theta_k = \theta_1 x^{k-1}$
#    for $k\ge1$ with $\theta_0 = \theta_1/(cx)$, so both series can be summed to
#    $\infty$ symbolically; the ratio is $cx/[(1-x)(1+(c-1)x)]$ **exactly**. The
#    identity therefore does not depend on stopping the ladder anywhere, and the
#    4000-layer numerical sum below is its independent arithmetic check.
#
# Nothing is more likely to embarrass a page than an identity that is true by
# construction, so the two sides are built from different objects: the left from
# eq. (30) as transcribed (and, in check 5, from the occupancy ladder itself), the
# right from the BET form as `J1.3` states it.
#

# %%
s1, s2, s3, s4, s5, mu = sp.symbols("s1 s2 s3 s4 s5 mu", positive=True)
S = [None, s1, s2, s3, s4, s5]
prod = lambda k: sp.prod([S[i] for i in range(1, k + 1)])

# (1) eq. (26) from eqs. (22)-(25), built from the theta_i rather than quoted
th1 = sp.Symbol("th1", positive=True)
theta = {0: th1 / (s1 * mu), 1: th1}
for k in range(2, 6):
    theta[k] = sp.prod([S[i] for i in range(2, k + 1)]) * mu ** (k - 1) * th1
inv22 = sum(theta[k] for k in range(0, 6))                 # eq. (22)
tot23 = sum(k * theta[k] for k in range(1, 6))             # eq. (23)
eq26_lhs = sp.simplify(tot23 / inv22)
eq26_rhs = (sum(k * prod(k) * mu ** k for k in range(1, 6))
            / (1 + sum(prod(k) * mu ** k for k in range(1, 6))))
r1 = sp.simplify(sp.together(eq26_lhs - eq26_rhs))
print("(1) eq. (26) - (eqs. 22-25 assembled)  =", r1)

# (2) eq. (30) as the coefficients of eq. (29)
ser = sp.series(sp.simplify(s1 * mu * (1 + sum(prod(k) * mu ** k for k in range(1, 6)))
                            / sum(k * prod(k) * mu ** k for k in range(1, 6))),
                mu, 0, 4).removeO().expand()
poly = sp.Poly(sp.expand(ser), mu)
got = [sp.simplify(poly.coeff_monomial(mu ** k)) for k in range(0, 4)]
a30 = s1 - 2 * s2
b30 = s2 * (4 * s2 - 3 * s3 - s1)
c30 = 2 * s2 * (6 * s2 * s3 - 2 * s3 * s4 + s1 * s2 - s1 * s3 - 4 * s2 ** 2)
r2 = [sp.simplify(sp.expand(got[0] - 1))] + [sp.simplify(sp.expand(g - pr))
                                             for g, pr in zip(got[1:], [a30, b30, c30])]
print("(2) constant, a, b, c residuals       =", r2)

# (3) the two printed collapse claims
cl1 = [sp.simplify(e.subs({s2: s1, s3: s1, s4: s1})) for e in (b30, c30)]
cl2 = sp.simplify(c30.subs({s3: s2, s4: s2}))
print("(3) all sigma equal -> (b, c)         =", cl1)
print("    sigma_3.. = sigma_2 -> c          =", cl2)

# (4) the BET identity
cc, xx = sp.symbols("c x", positive=True)
lang = (s1 / (1 / mu + a30 + b30 * mu)).subs({s3: s2, s4: s2})
lang = sp.simplify(lang.subs({mu: xx / s2, s1: cc * s2}))
bet = cc * xx / ((1 - xx) * (1 + (cc - 1) * xx))
r4 = sp.simplify(sp.expand(lang - bet))
print("(4) Case VI (sigma_3..=sigma_2) - BET =", r4)
print("    Case VI, factored                 =", sp.factor(lang))

# (5) the SAME ladder summed to INFINITY in closed form: no truncation anywhere.
# Checks (1)-(4) run through eq. (29) truncated at b, which is what Langmuir's own
# p.-1375 statement licenses. This one never truncates: with sigma_3.. = sigma_2 the
# occupancies of eqs. (22)/(23) are theta_k = theta_1 x^(k-1) for k >= 1 and
# theta_0 = theta_1/(c x), and both series are summed symbolically to infinity.
nlay = sp.Symbol("k", integer=True, positive=True)


def convergent(expr):
    """The |x| < 1 branch of a sympy Sum. Selected by testing each condition at
    x = 1/2, not by position, so it cannot silently pick the wrong branch."""
    e = sp.piecewise_fold(expr)
    if not isinstance(e, sp.Piecewise):
        return e
    for branch, cond in e.args:
        if cond is sp.true or bool(cond.subs(xx, sp.Rational(1, 2))):
            return branch
    raise AssertionError("no convergent branch")


inv_inf = 1 / (cc * xx) + convergent(sp.Sum(xx ** (nlay - 1), (nlay, 1, sp.oo)).doit())
tot_inf = convergent(sp.Sum(nlay * xx ** (nlay - 1), (nlay, 1, sp.oo)).doit())
ladder_inf = sp.simplify(tot_inf / inv_inf)
r5 = sp.simplify(sp.expand(ladder_inf - bet))
print("(5) INFINITE ladder, closed form      =", sp.factor(ladder_inf))
print("    infinite ladder - BET             =", r5)

for k, v in (("sym_eq26_residual", r1), ("sym_eq30_a_residual", r2[1]),
             ("sym_eq30_b_residual", r2[2]), ("sym_eq30_c_residual", r2[3]),
             ("sym_collapse_b_residual", cl1[0]), ("sym_collapse_c_residual", cl1[1]),
             ("sym_collapse2_c_residual", cl2), ("sym_caseVI_minus_BET", r4),
             ("sym_ladder_infinite_minus_BET", r5)):
    keep(k, abs(complex(sp.N(v.subs({s1: sp.Rational(7, 3), s2: sp.Rational(5, 4),
                                     s3: sp.Rational(9, 8), s4: sp.Rational(11, 10),
                                     s5: sp.Rational(13, 12), mu: sp.Rational(1, 7),
                                     cc: sp.Rational(37, 10), xx: sp.Rational(3, 17),
                                     th1: sp.Rational(1, 5)}))).real))

# a numerical route to the same statement, sharing no algebra with either closed form:
# model.caseVI_series sums eqs. (22)/(23) term by term
xt, ct = 0.20, 37.0
SERIES_MAX = 1e-13     # bound on the 4000-layer ladder against the closed form
keep("caseVI_series_vs_closed_form",
     abs(caseVI_series(xt, ct, 4000) / float(bet.subs({cc: ct, xx: xt})) - 1))
xf, cf, nf = 0.60, 37.0, 20
keep("caseVI_ladder_20layer_rel",
     abs(caseVI_series(xf, cf, nf) / float(bet.subs({cc: cf, xx: xf})) - 1))
print(f"\nnumerical ladder (4000 layers, summed term by term) against the closed form: "
      f"below {SERIES_MAX:.0e} relative at x={xt}, c={ct} (here "
      f"{RES['caseVI_series_vs_closed_form']:.1e}), which is summation round-off; the same "
      f"ladder truncated at {nf} layers at x={xf} gives "
      f"{RES['caseVI_ladder_20layer_rel']:.3e}, which is the truncation error and moves when "
      f"eq. (30) is broken.")
print("\nCONCLUSION: eqs. (29) and (30) of Langmuir 1918, with the paper's own stated "
      "assumption about the upper layers, ARE the BET isotherm. What Langmuir does NOT do, "
      "and this is the whole of the difference, is identify sigma_2 mu with p/p0 as a working "
      "variable, introduce c as a fitted constant, linearise, extract v_m, or fit eq. (29) to "
      "any of his own data. He writes it down, remarks that it goes to infinity at saturation, "
      "and moves on.")

# %% [markdown]
# ### 6.9 And therefore: the isotherm could not have decided the question
#
# If Langmuir's monolayer isotherm and Langmuir's multilayer isotherm are the two
# candidate explanations, then the interesting question is **where they disagree**,
# and by how much, at the pressures he actually reached.
#
# Eq. (39) makes that computable from the paper alone: $\sigma\mu = 1$ at
# saturation, so Table XXI's liquefied-gas lives invert through eq. (3) to
# saturation pressures,
#
# $$p_0 \;=\; \frac{\sqrt{MT}}{43.75\times10^{-6}\;\sigma_\mathrm{liq}}.$$
#
# Two of these recoveries deserve a sceptic's attention before the conclusion is
# drawn, and §6.9's cell gives them one: some of the 155 K entries invert to
# saturation pressures of tens of atmospheres, and Langmuir does not print the
# vapour pressures he used. **The conclusion does not depend on them.** The
# headline maximum comes from carbon dioxide, whose recovered $p_0$ is the
# *smallest* in the set, and dropping every 155 K row outright only *lowers* the
# maximum relative pressure and *raises* the scatter-to-separation ratio.
#
# **The paper prints no saturation pressure and no relative pressure anywhere.**
# That claim was checked, not assumed, and the search is stated so it can be
# repeated: every occurrence of the stem "saturat" in the extracted text layer was
# read in context, 38 of them by `pdftotext`, 37 by `pdftotext -layout`, the one
# difference being a hyphenation broken across a line, the theory section (journal
# pages 1361–1376) was read page by page on 300 ppi renders, and every table was
# cropped and read. (The text layer is untrustworthy for *digits* on this file, as
# §3 shows, but it is adequate for locating a word stem, and the hits were then
# read on the renders.) What the paper *does* say is qualitative and correct , 
# journal page 1365, "by using gases at pressures much below their saturation
# pressures", and journal page **1384**, "and with gases far below saturation".
# Langmuir knew where he was. He never says how far, and Table XXI is the only
# route to it.
#
# The ratio of the two isotherms, at the same $\sigma_1$ and $\sigma_2$, is
#
# $$\frac{\text{Case VI}}{\text{Case I}}
# = \frac{1+cx}{(1-x)\bigl(1+(c-1)x\bigr)},$$
#
# whose derivative in $x$ has numerator $c(c-1)x^2 + 2(c-1)x + 2 > 0$ for $c>1$.
# **It is strictly increasing**, so the maximum over a table is at that table's
# highest pressure, an argument, not a sampled maximum.
#
# **Two different scopes, kept apart.** The separation column needs a *fitted*
# isotherm, so it can only be built where the paper prints a single $a$, that is
# the 18 eq. (31) entries of Bulb A, and it excludes Tables XIII and XVI (the
# eq. 33 tables) and the blank bulb $A'$, which ran to *higher* pressures than
# $A$ on six tables. The relative pressure $p/p_0$ needs no fit and no bulb, only
# a pressure and Table XXI, so the cell computes it over **every pressure printed
# in the paper**, both bulbs, all tables, and the "anywhere" claim is reported
# from that wider set. The two are printed side by side below so the difference is
# visible rather than asserted; widening the scope does not move the maximum,
# because the headline entry (Table XI, CO₂ at 155 K) already sits in the narrower
# one.
#

# %%
rows = []
for (t, gname), grp in A.groupby(["table", "gas"]):
    T = grp.T_K.iloc[0]
    sl = liv.loc[gname, f"sigma_liq_{int(T)}K_s"]
    k = con[(con.table == t) & (con.gas == gname) & (con.model == "eq31")]
    if not np.isfinite(sl) or not len(k):
        continue
    a, b = k.a_header.iloc[0], k.b_cu_mm.iloc[0]
    p0 = np.sqrt(MW[gname] * T) / (K3 * sl)
    cpar = a * p0                                     # sigma_1/sigma_2
    p, q = grp.p_bars.values, grp.q_obs_cu_mm.values
    x = p / p0
    qfit = q_eq31(p, a, b)
    dep = (1 + cpar * x) / ((1 - x) * (1 + (cpar - 1) * x)) - 1
    sep = np.abs(dep * qfit)                          # cu.mm between the two theories
    rows.append((t, gname, T, p.max(), p0, p0 / ATM, x.max(), cpar,
                 dep.max() * 100, sep.max(), rms(q - qfit), rms(q - qfit) / sep.max()))
rp = pd.DataFrame(rows, columns=["table", "gas", "T_K", "p_max_bars", "p0_bars", "p0_atm",
                                 "x_max", "c", "departure_pct", "separation_cu_mm",
                                 "scatter_cu_mm", "scatter_over_separation"])
print(rp.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
keep("rp_entries", len(rp))

# ---- and now ANYWHERE really does mean anywhere -----------------------------
# The table above needs an eq. (31) fit, so it covers 18 of the 20 gas/table entries
# and Bulb A only: Tables XIII and XVI are the eq. (33) tables and have no single a,
# and the blank bulb A' ran to HIGHER pressures than A on six tables. p/p0 needs
# neither a fit nor a bulb, only a pressure and Table XXI, so it is computed over
# EVERY pressure in the paper, both bulbs, all tables, and reported from that set.
allx = []
for (t, gname, T, bulb), grp in iso.groupby(["table", "gas", "T_K", "bulb"]):
    sl = liv.loc[gname, f"sigma_liq_{int(T)}K_s"]
    if not np.isfinite(sl):
        continue
    p0 = np.sqrt(MW[gname] * T) / (K3 * sl)
    allx.append((t, gname, int(T), bulb, grp.p_bars.max(), p0, grp.p_bars.max() / p0))
xall = pd.DataFrame(allx, columns=["table", "gas", "T_K", "bulb", "p_max_bars", "p0_bars", "x_max"])
keep("x_entries_all", len(xall))
keep("x_rows_all", int(iso.p_bars.notna().sum()))
keep("x_max_anywhere", xall.x_max.max())
keep("x_max_eq31_bulbA_only", rp.x_max.max())
keep("x_max_eq33_tables", float(xall[xall.table.isin(["XIII", "XVI"])].x_max.max()))
keep("x_max_blank_bulb", float(xall[xall.bulb == "Aprime"].x_max.max()))
top = xall.loc[xall.x_max.idxmax()]
print(f"\nEVERY PRESSURE IN THE PAPER, both bulbs and both isotherm forms: "
      f"{int(RES['x_entries_all'])} table/gas/bulb entries covering "
      f"{int(RES['x_rows_all'])} printed pressures.")
print(xall.sort_values("x_max", ascending=False).head(6)
        .to_string(index=False, float_format=lambda v: f"{v:.5g}"))
print(f"  the eq. (33) tables (XIII, XVI), which the separation table above cannot host: "
      f"max p/p0 = {RES['x_max_eq33_tables']:.4g}")
print(f"  the blank bulb A', which ran HIGHER than A on six tables:  "
      f"max p/p0 = {RES['x_max_blank_bulb']:.4g}")
print(f"  the separation table's own 18 eq.(31) Bulb-A entries:      "
      f"max p/p0 = {RES['x_max_eq31_bulbA_only']:.4g}")
print(f"  => the maximum over EVERYTHING is {RES['x_max_anywhere']:.4g}, Table {top.table}, "
      f"{top.gas} at {top.T_K} K, bulb {top.bulb}. Widening the scope does not move it: the "
      f"headline entry is already in the separation table.")

keep("departure_max_pct", rp.departure_pct.max())
keep("scatter_over_sep_min", rp.scatter_over_separation.min())
keep("scatter_over_sep_t7", rp[rp.table == "VII"].scatter_over_separation.iloc[0])
keep("p0_N2_90K_bars", rp[rp.table == "VII"].p0_bars.iloc[0])
keep("p0_N2_90K_atm", rp[rp.table == "VII"].p0_atm.iloc[0])
keep("c_t7", rp[rp.table == "VII"].c.iloc[0])
best = rp.loc[rp.scatter_over_separation.idxmin()]
print(f"\nRecovered from Table XXI: p0(N2, 90 K) = {RES['p0_N2_90K_bars']:.4g} bars = "
      f"{RES['p0_N2_90K_atm']:.4g} atm, and c = sigma_1/sigma_2 = {RES['c_t7']:.4g} on Table VII.")
print(f"HIGHEST relative pressure reached ANYWHERE in the paper: p/p0 = "
      f"{RES['x_max_anywhere']:.4g} (that count is over every pressure printed, both bulbs). "
      f"At that point Case VI and Case I differ by {RES['departure_max_pct']:.4g} %.")
print(f"THE POINT: the smallest ratio of observational scatter to theory separation, over the "
      f"{int(RES['rp_entries'])} entries that can host the comparison, is "
      f"{RES['scatter_over_sep_min']:.1f} (Table {best.table}, "
      f"{best.gas}); on Table VII, the paper's best nitrogen set, it is "
      f"{RES['scatter_over_sep_t7']:.0f}. The monolayer and multilayer isotherms are NOT "
      f"distinguishable by these data, anywhere, by one to three orders of magnitude.")

# where they WOULD separate: root-find, two routes
cpar7 = RES["c_t7"]
x1 = unique_root(lambda z: departure(z, cpar7) - 0.01, 1e-12, 0.5, log=True, xtol=1e-17)
keep("x_1pct_rootfind", x1)
keep("x_1pct_closed_form", 1.0 / 101.0)           # the c -> infinity limit, where the ratio is 1/(1-x)
keep("x_1pct_two_route_rel", abs(x1 / (1 / 101) - 1))
keep("p_1pct_over_p_max", x1 * RES["p0_N2_90K_bars"] / 34.0)
grid = np.logspace(-6, np.log10(0.4), 40)
dep_grid = (1 + cpar7 * grid) / ((1 - grid) * (1 + (cpar7 - 1) * grid)) - 1
x_sampled = grid[np.argmin(np.abs(dep_grid - 0.01))]
keep("x_1pct_sampled", x_sampled)
keep("x_1pct_sampling_error_pct", abs(x_sampled / x1 - 1) * 100)
print(f"\n1 % separation is reached at p/p0 = {RES['x_1pct_rootfind']:.10f} (root-find) against "
      f"the closed form 1/101 = {RES['x_1pct_closed_form']:.10f} valid as c -> infinity; the two "
      f"differ by {RES['x_1pct_two_route_rel']:.3e}, which is the finite-c correction and not an "
      f"error. That pressure is {RES['p_1pct_over_p_max']:.0f} TIMES the highest Langmuir "
      f"reached on Table VII.")
print(f"Reading it off a 40-point log grid instead would give "
      f"{RES['x_1pct_sampled']:.6f}, wrong by {RES['x_1pct_sampling_error_pct']:.1f} %: that is "
      f"the break row in section 7.")

# --- and the row of that table that a sceptic should attack first ------------
# Some of the recovered p0 are large (tens of atmospheres), and Langmuir does not print
# the vapour pressures he used, so those entries of Table XXI are the least secure part
# of this argument. The conclusion must not depend on them, and it does not:
warm = rp[rp.T_K == 155]
cool = rp[rp.T_K == 90]
keep("x_max_90K_only", cool.x_max.max())
keep("departure_max_90K_only", cool.departure_pct.max())
keep("scatter_over_sep_min_90K_only", cool.scatter_over_separation.min())
keep("x_max_155K_permanent_gases",
     float(warm[warm.gas.isin(["N2", "O2", "Ar", "CO"])].x_max.max()))
print(f"\nSENSITIVITY, because the 155 K entries of Table XXI are the least secure input here "
      f"-- Langmuir does not print the vapour pressures behind them, and their inversion gives "
      f"p0 of tens of atmospheres for the permanent gases:")
print(f"  the headline maximum p/p0 = {RES['x_max_anywhere']:.4g} comes from Table XI, carbon "
      f"dioxide at 155 K, whose recovered p0 is the SMALLEST in the set ({rp[rp.table=='XI'].p0_atm.iloc[0]:.4g} atm) "
      f"and whose Table XXI life (16.5 s) is the largest.")
print(f"  DISCARDING EVERY 155 K ROW ENTIRELY: max p/p0 falls to {RES['x_max_90K_only']:.4g}, "
      f"the largest departure to {RES['departure_max_90K_only']:.4g} % and the smallest "
      f"scatter-to-separation ratio RISES to {RES['scatter_over_sep_min_90K_only']:.1f}. "
      f"The conclusion is unchanged and in fact strengthened.")
print(f"  and the four 155 K permanent-gas rows, the ones with the largest recovered p0, carry "
      f"the SMALLEST relative pressures in the whole paper (max "
      f"{RES['x_max_155K_permanent_gases']:.3g}), so discarding them can only lower the maximum.")

# --- the direction that WEAKENS the result, which the test above does not probe ---
# Dropping the 155 K rows makes the conclusion safer. The honest question is the other
# way round: how badly wrong would the recovered p0 have to be for the conclusion to
# FAIL? Root-found, not swept: shrink p0 of the tightest entry by a factor f (so x -> x/f
# and c -> c f, both as eq. (39) requires) and solve for scatter/separation = 1.
tight_rp = rp.loc[rp.scatter_over_separation.idxmin()]
_gr = A[(A.table == tight_rp.table) & (A.gas == tight_rp.gas)]
_k = con[(con.table == tight_rp.table) & (con.gas == tight_rp.gas) & (con.model == "eq31")].iloc[0]


SOS = lambda f, gr, kk_, p0_ref: sos_at_p0_factor(f, gr.p_bars.values, gr.q_obs_cu_mm.values,
                                                 kk_.a_header, kk_.b_cu_mm, p0_ref)
f_par = unique_root(lambda z: SOS(z, _gr, _k, tight_rp.p0_bars) - 1.0, 1e-6, 1.0, log=True,
                    xtol=1e-17)
keep("p0_shrink_factor_for_parity", 1.0 / f_par)
print(f"\nAND THE DIRECTION THAT WOULD WEAKEN IT, root-found rather than swept. The 155 K "
      f"entries invert Table XXI at a temperature where 'the liquefied gas' of eq. (39) is a "
      f"SUPERCOOLED liquid for CO2 (its triple point is above 155 K). Eq. (39) and BET's p0 "
      f"both mean the liquid, so the page's reading is the self-consistent one, but a reader "
      f"who thinks it too high should know how much too high it would have to be:")
print(f"  the tightest entry is Table {tight_rp.table}, {tight_rp.gas} at "
      f"{int(tight_rp.T_K)} K, at scatter/separation = {tight_rp.scatter_over_separation:.3f}. "
      f"Parity (= 1) needs the recovered p0 to be a factor "
      f"{RES['p0_shrink_factor_for_parity']:.1f} TOO HIGH.")
for ff in (1.0, 0.1, 1.0 / RES["p0_shrink_factor_for_parity"], 0.01):
    print(f"    p0 x {ff:9.6f}  ->  scatter/separation = "
          f"{SOS(ff, _gr, _k, tight_rp.p0_bars):9.4f}")
print(f"  A two-significant-figure life in Table XXI cannot be wrong by a factor "
      f"{RES['p0_shrink_factor_for_parity']:.0f}, so the conclusion is not close to the edge in "
      f"either direction.")
print(f"  BE CLEAR ABOUT WHAT THIS NUMBER IS. At x this small the separation is very nearly "
      f"proportional to 1/p0, so the factor comes out within "
      f"{abs(RES['p0_shrink_factor_for_parity'] / tight_rp.scatter_over_separation - 1) * 100:.1f} % "
      f"of the scatter/separation ratio itself: it is NOT independent information, it is that "
      f"same ratio restated in the units a sceptic of Table XXI would use. Its value is that it "
      f"is ROOT-FOUND on the exact expression rather than assumed proportional.")


# %% [markdown]
# **The sharper version of the same test: give Case VI both parameters and let
# it try.** Everything above compares the two isotherms *at Langmuir's own
# constants*. A sceptic can answer that the comparison was never fair, that Case VI
# was never allowed to fit. So the cell below refits **both** models to Table VII
# from scratch, each with two free parameters and each with its linear amplitude
# profiled out exactly, and puts the best Case VI (BET) against the best Case I:
#
# - Case I: $q = ab\,p/(1+ap)$, free $a$ and $ab$;
# - Case VI: $q = v_m\,cx/[(1-x)(1+(c-1)x)]$ with $x = p/p_0$ from Table XXI, free
#   $c$ and $v_m$.
#
# Both are found on a 4000-point log grid followed by a Brent refinement, so neither
# is a sampled optimum. The Case I result is also a **second, independent route to
# the fit of §6.6**, a different objective assembly, a different search, and the
# two are printed against each other.
#

# %%
# model.profile_rms: a log-grid scan of the one shape parameter, Brent-refined, so the
# reported optimum is not a grid point; it shares no code with profile_fit (different
# objective assembly, different search, no derivative).
SCAN = dict(n_grid=4000, xtol=1e-13)
SECOND_ROUTE_MAX = 1e-12    # Case I here against the fit of 6.6: both solved to xtol 1e-13
OWN_CURVE_MAX = 1e-9        # Case VI refitted to its own curve: round-off of the fit
p0_7 = RES["p0_N2_90K_bars"]
BET_BASIS = lambda p, cpar, p0v=None: (lambda z: cpar * z / ((1 - z) * (1 + (cpar - 1) * z)))(
    p / (p0v if p0v is not None else p0_7))

rL7, a_L7, ab_L7, gridL = profile_rms(lang_basis, p7, q7, 1e-6, 50.0, **SCAN)
rB7, c_B7, vm_B7, gridB = profile_rms(BET_BASIS, p7, q7, 1e-2, 1e9, **SCAN)
keep("t7_BET_free_over_L", rB7 / rL7)
keep("t7_BET_free_c", c_B7)
keep("t7_BET_free_vm", vm_B7)
keep("t7_L_second_route_rel",
     abs(rL7 / dsc[dsc.table == "VII"].rms_Langmuir.iloc[0] - 1))
print(f"CASE I  refitted, both parameters free: a = {a_L7:.6f}, b = {ab_L7 / a_L7:.5f}, "
      f"RMS = {rL7:.10f}")
print(f"        the SAME optimum by the independent route of section 6.6: relative difference "
      f"below {SECOND_ROUTE_MAX:.0e} (here {RES['t7_L_second_route_rel']:.1e}): two searches, two objectives, "
      f"one answer.")
print(f"CASE VI refitted, both parameters free: c = {c_B7:.6g}, v_m = {vm_B7:.5f}, "
      f"RMS = {rB7:.10f}")
print(f"        (grid minima before refinement: {gridL:.8f} and {gridB:.8f}, so neither number "
      f"above is a sampled optimum)")
print(f"\nRATIO RMS(Case VI, free) / RMS(Case I, free) = {RES['t7_BET_free_over_L']:.10f}")
print(f"GIVEN COMPLETE FREEDOM THE MULTILAYER MODEL IMPROVES THE FIT BY "
      f"{(1 - RES['t7_BET_free_over_L']) * 1e5:.2f} PARTS IN 100000. And look where it goes: "
      f"a p0 = {a_L7 * p0_7:.6g} against the fitted c = {c_B7:.6g}, and b = {ab_L7 / a_L7:.5f} "
      f"against the fitted v_m = {vm_B7:.5f}. Case VI, free to be anything, RETURNS Case I, "
      f"which is what x < 1e-5 forces it to do. This confirms section 6.9's conclusion by a "
      f"route that never mentions p/p0 separations at all.")

# does the test have any power? give it data that ARE multilayer.
p0_syn = p7.max() / 0.3                       # same pressures, but reaching x = 0.3


def own_curve_ratio(p0_fit=p0_syn, scan=SCAN):
    """RMS(Case VI free)/RMS(Case I free) on Case VI's OWN curve (c = 50, v_m = 25)
    generated at p0_syn; Case VI is fitted with p0_fit."""
    q_syn = 25.0 * BET_BASIS(p7, 50.0, p0_syn)
    rL_s = profile_rms(lang_basis, p7, q_syn, 1e-6, 50.0, **scan)[0]
    rB_s = profile_rms(lambda p, cpar: BET_BASIS(p, cpar, p0_fit), p7, q_syn, 1e-2, 1e9, **scan)[0]
    return rB_s / rL_s


keep("t7_BET_free_over_L_multilayer", own_curve_ratio())
print(f"\nPOWER CHECK (this is the section 7 break row): feed the same eleven pressures a "
      f"curve that Case VI generated at p/p0 up to 0.3, and the ratio collapses from "
      f"{RES['t7_BET_free_over_L']:.6f} to below {OWN_CURVE_MAX:.0e} (here "
      f"{RES['t7_BET_free_over_L_multilayer']:.1e}): Case VI fits its own curve to round-off. The "
      f"comparison CAN see multilayer adsorption. It does not see it in Langmuir's data "
      f"because it is not there to be seen at his pressures.")

# %%
fig, ax = plt.subplots(figsize=(7.2, 4.3))
xg = np.logspace(-7, -0.3, 400)
for tab, col in (("VII", "#1f6feb"), ("XI", "#d1242f")):
    r = rp[rp.table == tab].iloc[0]
    dep = ((1 + r.c * xg) / ((1 - xg) * (1 + (r.c - 1) * xg)) - 1) * 100
    ax.plot(xg, dep, "-", lw=1.7, color=col,
            label=f"Table {tab} ({r.gas}, {int(r.T_K)} K), c = {r.c:.3g}")
    ax.plot([r.x_max], [r.departure_pct], "o", ms=7, color=col)
ax.axvspan(1e-7, rp.x_max.max(), color="0.88", zorder=0)
ax.axhline(1.0, color="k", lw=0.8, ls=":")
ax.plot([RES["x_1pct_rootfind"]], [1.0], "k*", ms=11, zorder=5)
ax.annotate("every pressure in the paper\nlies in here", xy=(3e-5, 3e-3), fontsize=8.5)
ax.annotate("1 % separation\n(root-found)", xy=(RES["x_1pct_rootfind"], 1.0),
            xytext=(2.5e-3, 6.0), fontsize=8.5,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("$p/p_0$, recovered from the paper's own Table XXI through eq. (39)")
ax.set_ylabel("Case VI vs Case I  (per cent)")
ax.set_title("Where Langmuir's two isotherms disagree, and where he measured")
ax.legend(fontsize=8, frameon=False, loc="upper left")
fig.tight_layout(); plt.show()
print("The filled circles are the highest pressure each table reached. The separation is "
      "strictly increasing in p/p0 (proved above), so those circles ARE the maxima, not "
      "sampled ones.")

# %% [markdown]
# ### 6.10 The log-log slope claim
#
# Journal page 1384, on the log-log plot of Table VII:
#
# > "At the lowest pressures (3 bars) the slope of this line corresponds to an
# > exponent 0.684 while at the highest pressures (100 bars) the exponent decreases
# > to 0.20."
#
# For eq. (31), $d\ln q/d\ln p = 1/(1+ap)$ exactly. With Table VII's own
# $a = 0.156$ that is a one-line prediction, and the two halves of the sentence do
# not behave the same way.
#

# %%
a7 = con[(con.table == "VII") & (con.model == "eq31")].a_header.iloc[0]
slope = lambda p: 1.0 / (1.0 + a7 * p)
p_lo = float(clm.loc["loglog_pressure_low", "value"])
p_hi = float(clm.loc["loglog_pressure_high", "value"])
s_lo = float(clm.loc["loglog_slope_low", "value"])
s_hi = float(clm.loc["loglog_slope_high", "value"])
keep("slope_at_3bars", slope(p_lo))
keep("slope_at_3bars_rel_pct", (slope(p_lo) / s_lo - 1) * 100)
keep("slope_at_100bars", slope(p_hi))
keep("slope_100_factor", s_hi / slope(p_hi))
keep("p_for_slope_0684", unique_root(lambda p: slope(p) - s_lo, 1e-3, 500.0, log=True))
keep("p_for_slope_020", unique_root(lambda p: slope(p) - s_hi, 1e-3, 500.0, log=True))
keep("t7_p_max", A[A.table == "VII"].p_bars.max())
keep("slope_at_t7_pmax", slope(RES["t7_p_max"]))
print(f"at the printed 3 bars      : 1/(1+ap) = {RES['slope_at_3bars']:.4f}, printed {s_lo} "
      f"({RES['slope_at_3bars_rel_pct']:+.2f} %)  <- reproduced")
print(f"at the printed 100 bars    : 1/(1+ap) = {RES['slope_at_100bars']:.4f}, printed {s_hi} "
      f"-- the printed value is {RES['slope_100_factor']:.2f} x larger  <- NOT reproduced")
print(f"the pressure that gives exactly 0.684 : {RES['p_for_slope_0684']:.3f} bars")
print(f"the pressure that gives exactly 0.20  : {RES['p_for_slope_020']:.3f} bars")
print(f"Table VII's own highest pressure      : {RES['t7_p_max']:.1f} bars, where the slope is "
      f"{RES['slope_at_t7_pmax']:.4f}")
print("\nREPORTED, NOT REPAIRED. Pin what is NOT free: a = 0.156 is printed and reproduces the "
      "whole q_cal column (section 6.1), and 1/(1+ap) is eq. (31) differentiated. Given those, "
      "the two exponents 0.684 and 0.20 correspond to 2.96 and 25.6 bars, and 25.6 is NOT the "
      "100 bars printed beside 0.20. Two readings are possible and the page adopts neither: "
      "either the pressure is wrong for the sentence's own table (whose range stops at 34.0 "
      "bars), or '100 bars' is the study's nominal ceiling: the summary on journal page 1402 "
      "says the work used 'pressures of 100 bars ... or less', carried into a sentence that is "
      "otherwise about Table VII. The exponent 0.684 belongs to Table VII either way.")

# %% [markdown]
# ## 7. Validation
#
# ### 7.1 What each check can and cannot fail
#
# | check | class | can it fail? |
# |---|---|---|
# | $q_\mathrm{cal}$ from the printed $a$, $b$ | **goodness of fit** | only on transcription or arithmetic |
# | the four-constant chain (eqs. 3, 37, 38, and the cu. mm. conversion) | internal identity | yes: it is what caught the text layer's exponent |
# | $b' = b/s$ | internal identity | yes: two rows fail it |
# | $N_0$, $\beta$, $\sigma$, Table XXII, Table XXV | internal identity | yes: one $\beta$ cell fails |
# | **Langmuir against Freundlich, both refitted** | **discrimination** | **yes** |
# | **$\beta < 1$ on every entry** | **out-of-sample bound** | **yes** |
# | eq. (26), eq. (30), the collapse claims, Case VI = BET | symbolic identity | yes: a mistyped coefficient moves them |
# | eq. (14) marched against eq. (15) and eq. (16) | two routes | yes |
# | the separation of Case VI from Case I over the data | **the page's own result** | yes |
#
# ### 7.2 Claims
#
# Every number the page stands behind is registered once, with how it is checked
# (`docs/claims-design.md`): a **value** against the committed baseline within its
# stated relative tolerance, a **bound** against its stated limits on every run and
# every machine, an **exact** result by equality.
#
# The bounds are the quantities whose size is not physics. The nine symbolic
# residuals are zero by construction (sympy simplifies them to 0) and are held at
# exactly 0. The numerical ladder against the BET closed form, the two derivative-free
# refits, and the closed-form inversion of eq. (14) are floating-point round-off and
# are held below bounds set by their tolerances. `fit_two_route_max_rel` is the gap
# between two root-finders handed the same objective, held below the bound its
# tolerances set. The observed orders of the implicit march, from successive
# differences on the step ladder, are held between 0.9 and 1.1 because implicit Euler
# is first order, and the finest difference must stand well clear of the march's
# round-off floor. Where a break row moves a value by less than 5 %, the value's
# tolerance is tightened below that move.
#

# %%
# tolerances tightened below the move of the break row aimed at the value
RTOL = {"t7_BET_free_over_L": 1e-6, "t7_graphical_penalty": 1e-3, "slope_at_3bars": 2e-3,
        "beta_headroom": 5e-3, "slope_100_factor": 5e-3, "p_for_slope_020": 5e-3}
def value(k, **kw):
    claims.value(k, RES[k], rtol=RTOL.get(k, 0.05), **kw)

# transcription and the printed identity web
value("chain_max_rel_pct", unit="%")
claims.exact("n_obs_mismatches", int(RES["n_obs_mismatches"]),
             note="tables whose printed 'Number of obs.' differs from the rows transcribed")
value("qcal_max_abs_cu_mm", unit="cu.mm")
value("s_inlier_max_dev_pct", unit="%")
value("s_worst_dev_pct", unit="%")
value("N0_max_rel_pct", unit="%")
value("beta_max_abs_dev")
value("sigma_max_abs_rel_pct", unit="%")
value("t22_single_max_pct", unit="%")
value("t22_multi_geom_max_pct", unit="%")
claims.value("pt_N0_max_abs_dev", RES["pt_N0_max_abs_dev"], rtol=0.0, atol=5e-3,
             note="Table XXV N0 recomputed and rounded to the printed 0.01 against the "
                  "printed cells; 0 when every cell matches, at least 0.01 when one does not")
value("mono_n_rel_pct", unit="%")
value("t4_two_route_b_rel_pct", unit="%",
      note="two routes through printed, rounded constants; the gap is set by the printed digits")
value("t4_b_printed_dev_pct", unit="%")
# fit against test: the discrimination and its nulls
value("disc_median_ratio")
value("t7_F_over_L")
value("t7_null_mean")
value("t7_null_kp")
value("t7_graphical_penalty")
value("qF_column_max_dev_pct", unit="%")
claims.bound("fit_two_route_max_rel", RES["fit_two_route_max_rel"], max=FIT_TWO_ROUTE_MAX,
             basis="tolerance",
             note=f"relative gap between newton (tol {FIT['newton_tol']:g}) and Brent (xtol "
                  f"{FIT['xtol']:g}) on the same central-difference dS, worst over the fits of "
                  "sections 6.2 and 6.6")
claims.bound("t7_L_second_route_rel", RES["t7_L_second_route_rel"], max=SECOND_ROUTE_MAX,
             basis="tolerance",
             note="Table VII Case I RMS by the derivative-free scan against the fit of 6.6; "
                  "both optima located to xtol 1e-13")
claims.bound("t7_BET_free_over_L_multilayer", RES["t7_BET_free_over_L_multilayer"],
             max=OWN_CURVE_MAX, basis="tolerance",
             note="Case VI refitted to its own curve reaching p/p0 = 0.3: round-off of the fit")
# the out-of-sample bound
value("beta_max")
value("beta_headroom")
claims.exact("beta_n_above_one", int(RES["beta_n_above_one"]),
             margin=float((1.0 - g.beta_recomputed).abs().min()), min_margin=0.05,
             note="entries with recomputed beta >= 1; the margin is the smallest |beta - 1|")
claims.count("beta_n_above_one_textlayer", g.beta_recomputed * RES["eq37_textlayer_factor"], 1.0,
             min_margin=0.05,
             note="entries above the monolayer bound if eq. (37) were read 25.2e16; the margin is "
                  "the smallest relative distance of a scaled beta from 1")
claims.count("disc_n_langmuir_wins", dsc.F_over_L, 1.0, min_margin=0.05,
             note="tables on which the best Freundlich RMS exceeds Langmuir's; the margin is the "
                  "smallest relative distance of an RMS ratio from 1")
claims.count("beta_n_exact", bd, 5e-3, above=False, min_margin=0.05,
             note="rows whose printed beta equals N0 over the monolayer count rounded to the "
                  "printed 0.01")
claims.count("sigma_n_below", srel, 0.0, above=False, min_margin=0.0,
             note="cells whose printed sigma lies below eq. (38); a sign count, so the margin is "
                  "the smallest |deviation| in percent")
claims.exact("sigma_trunc_reproduces", int(RES["sigma_trunc_reproduces"]),
             note="cells reproduced exactly by cutting eq. (38) at the cell's printed granularity")
claims.exact("sigma_round_reproduces", int(RES["sigma_round_reproduces"]),
             note="cells reproduced exactly by rounding eq. (38) at the cell's printed granularity")
claims.exact("order_breaks_recomputed_beta", int(RES["order_breaks_recomputed_beta"]),
             note="sets in which the order by recomputed beta differs from the order by b'")
# the symbolic family: zero by construction
for k in ("sym_eq26_residual", "sym_eq30_a_residual", "sym_eq30_b_residual",
          "sym_eq30_c_residual", "sym_collapse_b_residual", "sym_collapse_c_residual",
          "sym_collapse2_c_residual", "sym_caseVI_minus_BET", "sym_ladder_infinite_minus_BET"):
    claims.bound(k, RES[k], max=0.0, basis="identity",
                 note="sympy simplifies the difference to 0, evaluated at rational arguments")
claims.bound("caseVI_series_vs_closed_form", RES["caseVI_series_vs_closed_form"],
             max=SERIES_MAX, basis="chosen",
             note="4000-layer ladder summed term by term against the BET closed form at x = 0.2; "
                  "truncation is x^4000, so the gap is summation round-off")
value("caseVI_ladder_20layer_rel", note="truncation error of the 20-layer ladder at x = 0.6")
# the rate equation
value("t_half_eq15_corrected")
value("eq14_vs_eq15corrected_max_abs", note="discretisation error, 1600 implicit steps")
claims.bound("eq14_time_order", RES["eq14_time_order"], min=0.9, max=1.1, basis="theory",
             note="observed order of implicit Euler from successive differences of the "
                  "trajectory on 200 to 1600 steps, the one furthest from 1; first order")
claims.bound("eq14_ladder_over_roundoff", RES["eq14_ladder_over_roundoff"], min=ROUNDOFF_MIN,
             basis="chosen",
             note="finest successive difference of the trajectory over the march's floor "
                  "n_steps * max(newton tol, eps) * theta1; the orders mean something only "
                  "far above that floor")
value("t_half_march_vs_eq16_rel", note="discretisation error, 1600 implicit steps")
claims.bound("t_half_march_order", RES["t_half_march_order"], min=0.9, max=1.1, basis="theory",
             note="observed order of the half-coverage time from successive differences on "
                  "200 to 1600 steps, the one furthest from 1; first order")
claims.bound("t_half_richardson_vs_eq16_rel", RES["t_half_richardson_vs_eq16_rel"],
             max=RICHARDSON_MAX, basis="chosen",
             note="Richardson value (observed order) of the marched half-coverage time against "
                  "eq. (16): the march converges to eq. (16), not merely converges; the gap "
                  "is dominated by the error of the order estimate, and the bound covers it")
claims.bound("t_half_closed_form_vs_eq16_rel", RES["t_half_closed_form_vs_eq16_rel"],
             max=1e-14, basis="tolerance",
             note="eq. (14)'s closed form inverted at half coverage against eq. (16), an "
                  "algebraic identity; the gap is the root finder's tolerance")
# the page's own result
for k in ("x_max_anywhere", "departure_max_pct", "scatter_over_sep_min",
          "scatter_over_sep_t7", "p0_N2_90K_atm", "x_1pct_rootfind"):
    value(k)
value("x_1pct_two_route_rel", note="the finite-c correction to the c -> infinity limit 1/101")
for k in ("x_1pct_sampling_error_pct", "p_1pct_over_p_max", "x_max_90K_only",
          "departure_max_90K_only", "scatter_over_sep_min_90K_only",
          "x_max_155K_permanent_gases", "t7_BET_free_over_L", "p0_shrink_factor_for_parity",
          # the log-log slope claim
          "slope_at_3bars", "slope_100_factor", "p_for_slope_020"):
    value(k)
print(f"{len(claims.metrics)} claims registered")

# %% [markdown]
# ### 7.3 The break table
#
# Every claim needs a row that fails its own check. Each row below injects one defect,
# a mistyped constant or coefficient, a misread table cell, a different estimator or a
# coarser march, and recomputes every claim the defect reaches. A row whose
# computation cannot be completed under its defect records `nan` for its claims,
# which fails them, and stays in the table.
#
# One row is a sensitivity row, expected to fail nothing: "1 % separation from the
# c -> infinity closed form 1/101" shows that the finite-$c$ correction is far inside
# the tolerance of `x_1pct_rootfind`. Two invariances are stated rather than run as
# rows, because no defect can move them: the scatter-to-separation ratio is
# scale-free in $q$, so rescaling every $b'$ by a common factor leaves it unchanged,
# and a symbolic identity is invariant under a renaming of its symbols.
#

# %%
def chain_max(r=R, nav=NAV):
    """The worst residual of the four-constant chain of section 3, with R and N
    entering every identity they appear in."""
    derived = [1 / np.sqrt(2 * np.pi * r), 1 / K3, ATM * 1e-3 / (r * 293.15), nav * CUMM]
    printed = [K3, K38, CUMM, K37]
    return max(abs(pv / dv - 1) * 100 for pv, dv in zip(printed, derived))


def t7_contest(q7_):
    """Table VII's F/L ratio and both nulls on a q column."""
    rL = fit_langmuir(p7, q7_, **FIT)[2]
    rF = fit_freundlich(p7, q7_, **FIT)[2]
    return {"t7_F_over_L": rF / rL, "t7_null_mean": rms(q7_ - q7_.mean()) / rL,
            "t7_null_kp": rms((np.dot(p7, q7_) / np.dot(p7, p7)) * p7 - q7_) / rL}


def sym_val(e, sub=None):
    """|e| at rational arguments, after sympy simplification."""
    sub = sub or {s1: sp.Rational(7, 3), s2: sp.Rational(5, 4), s3: sp.Rational(9, 8),
                  s4: sp.Rational(11, 10), s5: sp.Rational(13, 12), mu: sp.Rational(1, 7)}
    return abs(complex(sp.N(sp.simplify(e).subs(sub))).real)


def xxi_inverted():
    """Table XXI's lives inverted, carried to every claim that uses them."""
    liv_b = liv.copy()
    for cc_ in ("sigma_liq_90K_s", "sigma_liq_155K_s"):
        liv_b[cc_] = 1.0 / liv_b[cc_]
    t22_bad = []
    for name, (ads, T, lc, rc) in sets.items():
        for gg in MONO:
            pr = liv.loc[gg, rc]
            if not np.isfinite(pr):
                continue
            ss = c[(c.adsorbent == ads) & (c.T_K == T) & (c.gas == gg)].sigma_s.dropna()
            ss = ss[ss > 2e4] if len(ss) > 1 else ss
            if len(ss) == 1:
                t22_bad.append(abs(pr / (float(ss.mean()) / liv_b.loc[gg, lc]) - 1) * 100)
    p0_b = np.sqrt(28 * 90) / (K3 * liv_b.loc["N2", "sigma_liq_90K_s"])
    dep_b = departure(34.0 / p0_b, k7.a_header * p0_b)
    sos_b = float(RES["t7_rms_printed"] / abs(dep_b * q_eq31(34.0, k7.a_header, k7.b_cu_mm)))
    return {"t22_single_max_pct": float(np.max(t22_bad)), "p0_N2_90K_atm": p0_b / ATM,
            "scatter_over_sep_t7": sos_b, "scatter_over_sep_min": sos_b,
            "x_max_90K_only": float(RES["x_max_90K_only"]
                                    * (K3 * liv.loc["N2", "sigma_liq_90K_s"]) ** 2),
            "x_max_155K_permanent_gases": float(RES["x_max_155K_permanent_gases"]
                                                * (K3 * liv.loc["N2", "sigma_liq_155K_s"]) ** 2)}


def departure_at_mean_pressure():
    dep_all, dep_90 = [], []
    for _, r in rp.iterrows():
        grp = A[(A.table == r.table) & (A.gas == r.gas)]
        d = departure(grp.p_bars.mean() / r.p0_bars, r.c) * 100
        dep_all.append(d)
        if r.T_K == 90:
            dep_90.append(d)
    return {"departure_max_pct": float(np.max(dep_all)),
            "departure_max_90K_only": float(np.max(dep_90))}


def slope_with_a(a_):
    sl = lambda pv: 1.0 / (1.0 + a_ * pv)
    return {"slope_at_3bars": sl(p_lo), "slope_100_factor": s_hi / sl(p_hi),
            "p_for_slope_020": unique_root(lambda pv: sl(pv) - s_hi, 1e-3, 500.0, log=True)}


def p0_margin_90K():
    t90 = rp[rp.T_K == 90].loc[rp[rp.T_K == 90].scatter_over_separation.idxmin()]
    g90 = A[(A.table == t90.table) & (A.gas == t90.gas)]
    k90 = con[(con.table == t90.table) & (con.gas == t90.gas) & (con.model == "eq31")].iloc[0]
    return 1.0 / unique_root(lambda z: SOS(z, g90, k90, t90.p0_bars) - 1.0, 1e-9, 1.0,
                             log=True, xtol=1e-17)


def one_fit_gap(**kw):
    log = []
    fit_langmuir(p7, q7, "Table VII", log=log, **dict(FIT, **kw))
    return log[0][1]


TEXT = 25.2e16                                   # the text layer's eq. (37) coefficient
q7b = q7.copy(); q7b[0] = 30.0                   # Table VII top row 33.0 -> 30.0
bad_mono = dict(MONO); bad_mono["CH4"], bad_mono["CO2"] = MONO["CO2"], MONO["CH4"]
x_1pct = lambda cpar: unique_root(lambda z: departure(z, cpar) - 0.01, 1e-12, 0.5, log=True,
                                  xtol=1e-17)
ROWS = [
    ("eq. (37) coefficient 25.2e15 -> 25.2e16 (the text layer's reading)", {
        "N0_max_rel_pct": lambda: float(((c.N0_e15 / (TEXT * c.b_prime_cu_mm_per_m2 / 1e4 / 1e15)
                                          - 1).abs() * 100).max()),
        "beta_n_above_one": lambda: int((g.beta_recomputed * TEXT / K37 >= 1).sum())}),
    ("same, propagated to the bound", {
        "beta_max": lambda: float((g.beta_recomputed * TEXT / K37).max())}),
    ("N 6.06e23 -> 6.60e23 (digit transposition)", {
        "chain_max_rel_pct": lambda: chain_max(nav=6.60e23)}),
    ("R 83.2e6 -> 8.314e7 (SI-flavoured erg value)", {
        "chain_max_rel_pct": lambda: chain_max(r=8.314e7)}),
    ("mica surface 5750 -> 5075 sq.cm (transposition)", {
        "s_inlier_max_dev_pct": lambda: float(
            (inl.b_cu_mm / inl.b_prime_cu_mm_per_m2 * 1e4
             / inl.adsorbent.map({"mica": 5075.0, "glass": 1966.0}) - 1).abs().max() * 100)}),
    ("methane/CO2 monolayer counts swapped", {
        "beta_max": lambda: float((g.N0_e15 / g.gas.map(bad_mono)).max()),
        "beta_headroom": lambda: float(1.0 / (g.N0_e15 / g.gas.map(bad_mono)).max())}),
    ("eq. (38) sqrt(MT) -> MT", {
        "sigma_max_abs_rel_pct": lambda: float(
            ((sg.sigma_s / (K38 * sg.a_summary * (sg.gas.map(MW) * sg.T_K)) - 1) * 100).abs().max())}),
    ("Table XXI sigma_liq inverted (a plausible misreading of 'relative life')", {
        "x_max_anywhere": lambda: 34.0 / (np.sqrt(28 * 90)
                                         / (K3 / liv.loc["N2", "sigma_liq_90K_s"]))}),
    ("Freundlich NOT refitted (the paper's own 8.4 p^0.417 everywhere)", {
        "disc_median_ratio": lambda: float(discrimination(A, freundlich="printed").F_over_L.median()),
        "disc_n_langmuir_wins": lambda: int(
            (discrimination(A, freundlich="printed").F_over_L > 1).sum())}),
    ("Langmuir fitted by Langmuir's linearisation instead of on q", {
        "disc_median_ratio": lambda: float(discrimination(A, langmuir="linearised").F_over_L.median()),
        "disc_n_langmuir_wins": lambda: int(
            (discrimination(A, langmuir="linearised").F_over_L > 1).sum())}),
    ("eq. (30): sign of sigma_1 in b flipped", {
        "sym_caseVI_minus_BET": lambda: sym_val(
            sp.simplify((s1 / (1 / mu + a30 + s2 * (4 * s2 - 3 * s3 + s1) * mu))
                        .subs({s3: s2, s4: s2}).subs({mu: xx / s2, s1: cc * s2})) - bet,
            {cc: sp.Rational(37, 10), xx: sp.Rational(3, 17)})}),
    ("eq. (30): the 4 sigma_2^2 in c changed to 3 sigma_2^2", {
        "sym_collapse_c_residual": lambda: sym_val(
            (2 * s2 * (6 * s2 * s3 - 2 * s3 * s4 + s1 * s2 - s1 * s3 - 3 * s2 ** 2))
            .subs({s2: s1, s3: s1, s4: s1}))}),
    ("1 % separation sampled on a 40-point log grid, not root-found", {
        "x_1pct_rootfind": lambda: RES["x_1pct_sampled"],
        "p_1pct_over_p_max": lambda: RES["x_1pct_sampled"] * RES["p0_N2_90K_bars"] / 34.0}),
    ("eq. (15) taken as printed", {
        "t_half_eq15_corrected": lambda: t_eq15_printed(0.5 * theta1, theta1, k_rate)}),
    ("eq. (14) marched with 200 steps instead of 1600", {
        "eq14_vs_eq15corrected_max_abs":
            lambda: rate_study(nsteps=(25, 50, 100, 200))["eq14_vs_eq15corrected_max_abs"]}),
    ("half-coverage root-found on a 200-step march instead of 1600", {
        "t_half_march_vs_eq16_rel":
            lambda: rate_study(nsteps=(25, 50, 100, 200))["t_half_march_vs_eq16_rel"]}),
    ("eq. (14) marched with a midpoint rule instead of implicit Euler", {
        "eq14_time_order": lambda: rate_study(scheme="midpoint")["eq14_time_order"],
        "t_half_march_order": lambda: rate_study(scheme="midpoint")["t_half_march_order"]}),
    ("each implicit step solved to a newton step of 1e-6", {
        "eq14_ladder_over_roundoff": lambda: rate_study(step_tol=1e-6)["eq14_ladder_over_roundoff"]}),
    ("eq. (14) without its evaporation term (k = alpha mu only)", {
        k_: (lambda k_=k_: rate_study(k_march=sig_mu)[k_])
        for k_ in ("t_half_march_vs_eq16_rel", "t_half_richardson_vs_eq16_rel",
                   "t_half_closed_form_vs_eq16_rel")}),
    ("Table IV: printed b = 58.3 used instead of the two-point solution", {
        "t4_two_route_b_rel_pct": lambda: abs(58.3 / RES["t4_b_from_bprime"] - 1) * 100}),
    ("Table IV b 58.3 -> 57.5", {
        "t4_b_printed_dev_pct": lambda: (57.5 / RES["t4_b_exact"] - 1) * 100}),
    ("Table VII top row q_obs 33.0 -> 30.0", {
        k_: (lambda k_=k_: t7_contest(q7b)[k_]) for k_ in ("t7_F_over_L", "t7_null_mean", "t7_null_kp")}),
    ("one Bulb A row dropped from Table VIII", {
        "n_obs_mismatches": lambda: int(
            (lambda ch_: (ch_.printed_n_obs != ch_.transcribed).sum())(pd.concat(
                [stated, A.drop(A[A.table == "VIII"].index[:1]).groupby(["table", "gas"]).size()
                 .rename("transcribed")], axis=1)))}),
    ("Table XX beta 0.36 accepted rather than recomputed", {
        "order_breaks_recomputed_beta": lambda: RES["order_breaks_printed_beta"]}),
    ("platinum area 312 -> 156 sq.cm (one face)", {
        "pt_N0_max_abs_dev": lambda: float(
            (ptn.N0_e15 - (K37 * ptn.q_cu_mm / 156.0 / 1e15).round(2)).abs().max())}),
    ("Table XXII combined with the arithmetic mean", {
        "t22_multi_geom_max_pct": lambda: RES["t22_multi_arith_max_pct"]}),
    ("monolayer chain: cube root replaced by square root", {
        "mono_n_rel_pct": lambda: float(
            (1.0 / (Vm_N2 / NAV) ** (2 / 2) / float(clm.loc["monolayer_N2", "value"]) - 1) * 100)}),
    ("separation read at the mean pressure instead of the maximum", {
        k_: (lambda k_=k_: departure_at_mean_pressure()[k_])
        for k_ in ("departure_max_pct", "departure_max_90K_only")}),
    ("eq. (22): the bare fraction theta_0 omitted from the inventory", {
        "sym_eq26_residual": lambda: sym_val(
            tot23 / sum(theta[k] for k in range(1, 6)) - eq26_rhs,
            {s1: sp.Rational(7, 3), s2: sp.Rational(5, 4), s3: sp.Rational(9, 8),
             s4: sp.Rational(11, 10), s5: sp.Rational(13, 12), mu: sp.Rational(1, 7),
             th1: sp.Rational(1, 5)})}),
    ("1 % separation from the c -> infinity closed form 1/101", {
        "x_1pct_rootfind": lambda: 1.0 / 101.0}),
    ("Table VIII q_cal 84.8 -> 48.8 (digit transposition)", {
        "qcal_max_abs_cu_mm": lambda: float(np.abs(
            model_for("VIII", "CH4")["eq31"](A[A.table == "VIII"].p_bars.values)
            - np.r_[48.8, A[A.table == "VIII"].q_cal_cu_mm.values[1:]]).max())}),
    ("Table IV b' repaired to 101.4 instead of reported", {
        "s_worst_dev_pct": lambda: float(
            (c.b_cu_mm / c.b_prime_cu_mm_per_m2.where(c.table != "IV", 101.4) * 1e4
             / c.s_printed - 1).max() * 100)}),
    ("Table XX beta repaired to 0.46 instead of reported", {
        "beta_max_abs_dev": lambda: float(
            (betdf.beta.where(betdf.beta != 0.36, 0.46) - betdf.beta_from_N0.round(2)).abs().max()),
        "beta_n_exact": lambda: int(
            ((betdf.beta.where(betdf.beta != 0.36, 0.46) - betdf.beta_from_N0.round(2)).abs()
             < 5e-3).sum())}),
    ("eq. (37) text-layer coefficient read 25.2e15 (factor 1)", {
        "beta_n_above_one_textlayer": lambda: textlayer_count(1.0)}),
    ("eq. (38) coefficient 5 % low", {
        "sigma_n_below": lambda: sigma_n_below(0.95 * K38),
        "sigma_trunc_reproduces": lambda: sigma_reproduces("trunc", 0.95 * K38),
        "sigma_round_reproduces": lambda: sigma_reproduces("round", 0.95 * K38)}),
    ("sigma reproduction tested by rounding instead of cutting", {
        "sigma_trunc_reproduces": lambda: sigma_reproduces("round")}),
    ("sigma reproduction tested by cutting instead of rounding", {
        "sigma_round_reproduces": lambda: sigma_reproduces("trunc")}),
    ("Table XXI sigma_liq inverted", {
        k_: (lambda k_=k_: xxi_inverted()[k_])
        for k_ in ("t22_single_max_pct", "p0_N2_90K_atm", "scatter_over_sep_t7",
                   "scatter_over_sep_min", "x_max_90K_only", "x_max_155K_permanent_gases")}),
    ("printed q_F column checked against 8.4 p^0.42 instead of p^0.417", {
        "qF_column_max_dev_pct": lambda: float(
            (np.abs(q_freundlich(p7, 8.4, 0.42) - d7.q_F_cu_mm.values)
             / d7.q_F_cu_mm.values).max() * 100)}),
    ("Table VII optimum taken on the p/q linearisation", {
        "t7_graphical_penalty": lambda: (lambda sl, ic: rms_printed
                                         / rms(q_eq31(p7, sl / ic, 1 / sl) - q7))(
            *np.polyfit(p7, p7 / q7, 1))}),
    ("eq. (30): a = sigma_1 - 2 sigma_2 -> sigma_1 - sigma_2", {
        "sym_eq30_a_residual": lambda: sym_val(got[1] - (s1 - s2))}),
    ("eq. (30): b's 4 sigma_2 -> 3 sigma_2", {
        "sym_eq30_b_residual": lambda: sym_val(got[2] - s2 * (3 * s2 - 3 * s3 - s1))}),
    ("eq. (30): c's 6 sigma_2 sigma_3 -> 5 sigma_2 sigma_3", {
        "sym_eq30_c_residual": lambda: sym_val(
            got[3] - 2 * s2 * (5 * s2 * s3 - 2 * s3 * s4 + s1 * s2 - s1 * s3 - 4 * s2 ** 2))}),
    ("collapse claim 1 evaluated on b with the 3 sigma_3 dropped", {
        "sym_collapse_b_residual": lambda: sym_val((s2 * (4 * s2 - s1)).subs({s2: s1, s3: s1, s4: s1}))}),
    ("collapse claim 2 evaluated with sigma_4 left free", {
        "sym_collapse2_c_residual": lambda: sym_val(c30.subs({s3: s2}))}),
    ("infinite ladder built with no distinguished first layer", {
        "sym_ladder_infinite_minus_BET": lambda: sym_val(
            tot_inf / (1 / xx + (inv_inf - 1 / (cc * xx))) - bet,
            {cc: sp.Rational(37, 10), xx: sp.Rational(3, 17)})}),
    ("ladder: the bare fraction dropped from the inventory", {
        "caseVI_series_vs_closed_form": lambda: abs(
            caseVI_series(xt, ct, 4000, bare=False) / float(bet.subs({cc: ct, xx: xt})) - 1)}),
    ("ladder truncated at 5 layers instead of 20", {
        "caseVI_ladder_20layer_rel": lambda: abs(
            caseVI_series(xf, cf, 5) / float(bet.subs({cc: cf, xx: xf})) - 1)}),
    ("Table VII a 0.156 -> 0.15", {
        k_: (lambda k_=k_: slope_with_a(0.15)[k_])
        for k_ in ("slope_at_3bars", "slope_100_factor", "p_for_slope_020")}),
    ("c halved (sigma_1/sigma_2 from a different table)", {
        "x_1pct_two_route_rel": lambda: abs(x_1pct(RES["c_t7"] / 2) / (1 / 101) - 1)}),
    ("the sampling grid refined from 40 to 400 points", {
        "x_1pct_sampling_error_pct": lambda: (lambda gr: abs(
            gr[np.argmin(np.abs(departure(gr, cpar7) - 0.01))] / x1 - 1) * 100)(
            np.logspace(-6, np.log10(0.4), 400))}),
    ("every 155 K row kept rather than dropped", {
        "scatter_over_sep_min_90K_only": lambda: rp.scatter_over_separation.min()}),
    ("Table VII's q replaced by Case VI's own curve reaching p/p0 = 0.3", {
        "t7_BET_free_over_L": own_curve_ratio}),
    ("Case VI fitted to its own curve with Table VII's p0", {
        "t7_BET_free_over_L_multilayer": lambda: own_curve_ratio(p0_fit=p0_7)}),
    ("the derivative-free refit stopped at its grid minimum", {
        "t7_L_second_route_rel": lambda: abs(
            profile_rms(lang_basis, p7, q7, 1e-6, 50.0, **dict(SCAN, refine=False))[0]
            / dsc[dsc.table == "VII"].rms_Langmuir.iloc[0] - 1)}),
    ("p0 margin computed on the 90 K entries only", {
        "p0_shrink_factor_for_parity": p0_margin_90K}),
    ("newton tol and backward_tol loosened to 1e-3 on one fit", {
        "fit_two_route_max_rel": lambda: one_fit_gap(newton_tol=1e-3, backward_tol=1e-3)}),
]
SENSITIVITY = ["1 % separation from the c -> infinity closed form 1/101"]

claims.define_rows([r[0] for r in ROWS])
BREAK = []
for row, studies in ROWS:
    for k_, fn in studies.items():
        try:
            out = fn()
        except Exception as exc:                 # noqa: BLE001 - a failed study is a result
            print(f"  {row}: {k_} raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
            out = np.nan
        caught_ = claims.breaks(row, k_, out)
        base_ = claims.metrics[k_]["value"]
        if claims.metrics[k_]["kind"] == "bound":   # round-off digits differ between machines
            base_, out = f"(here {base_:.1e})", f"(here {out:.1e})"
        BREAK.append((row, k_, base_, out, caught_))

bt = pd.DataFrame(BREAK, columns=["injected defect", "metric", "baseline", "with defect", "caught"])
pd.set_option("display.max_colwidth", 70)
print(bt.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
idle = [r for r in bt["injected defect"].unique()
        if not bt[(bt["injected defect"] == r) & bt.caught].shape[0]]
print(f"\n{len(ROWS)} rows defined, {bt['injected defect'].nunique()} in the table; "
      f"rows that fail no claim: {[r for r in idle if r not in SENSITIVITY] or 'none'}; "
      f"sensitivity rows (expected to fail none): {[r for r in idle if r in SENSITIVITY] or 'none'}")
missing = [n for n in claims.metrics if not bt[(bt.metric == n) & bt.caught].shape[0]]
print(f"claims no row fails: {missing or 'none'}")

# %%
# The prose literals of the README and this notebook, recorded so that a cut in prose
# precision shows in CI's loosening check (agreement.json "quoted").
SAYS = {
    "t7_F_over_L": ["3.82", "3.8235"], "disc_median_ratio": ["2.12", "2.1211"], "beta_headroom": ["1.1667"],
    "beta_max": ["0.857"], "beta_n_above_one": ["0"], "x_max_anywhere": ["below 1.6e-3"],
    "departure_max_pct": ["0.160"], "scatter_over_sep_min": ["12.2"],
    "t7_BET_free_over_L": ["0.99999"], "p0_shrink_factor_for_parity": ["12.14"],
    "scatter_over_sep_t7": ["1104"], "qcal_max_abs_cu_mm": ["0.434"],
    "t4_b_printed_dev_pct": ["1.40"], "t4_two_route_b_rel_pct": ["0.011"],
    "chain_max_rel_pct": ["0.136"], "sigma_max_abs_rel_pct": ["3.77"],
    "slope_at_3bars": ["0.6812"], "p_for_slope_020": ["25.64"], "slope_100_factor": ["3.32"],
    "x_1pct_rootfind": ["0.0098993"], "p_1pct_over_p_max": ["1114"],
    "x_1pct_sampling_error_pct": ["6.26"], "t7_graphical_penalty": ["1.063"],
    "qF_column_max_dev_pct": ["1.33"], "p0_N2_90K_atm": ["3.775"],
    "eq14_vs_eq15corrected_max_abs": ["1.78e-3"], "eq14_time_order": ["0.987"], "t_half_march_order": ["1.005"],
    "t_half_march_vs_eq16_rel": ["1.875e-3"], "s_worst_dev_pct": ["1.391"],
    "mono_n_rel_pct": ["0.451"], "t22_single_max_pct": ["2.27"],
    "fit_two_route_max_rel": ["below 1e-9"], "caseVI_series_vs_closed_form": ["below 1e-13"],
    "t7_L_second_route_rel": ["below 1e-12"], "t7_BET_free_over_L_multilayer": ["below 1e-9"],
    "t_half_closed_form_vs_eq16_rel": ["at most 1e-14"],
    "t_half_richardson_vs_eq16_rel": ["at most 5e-5"],
    "t7_null_mean": ["19.4"], "t7_null_kp": ["26.8"], "x_max_90K_only": ["7.7e-4"],
    "departure_max_90K_only": ["0.082"], "scatter_over_sep_min_90K_only": ["31.1"],
    "x_max_155K_permanent_gases": ["1.2e-6"],
    "beta_n_above_one_textlayer": ["17"], "disc_n_langmuir_wins": ["13"], "beta_n_exact": ["22"],
    "sigma_n_below": ["18"], "sigma_trunc_reproduces": ["15"], "sigma_round_reproduces": ["11"],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)

# Numbers the prose states that are not claims: counts compared exactly, the rest relative.
def close(name, quoted, key, tol=5e-3):
    got_ = RES[key]
    if abs(got_ - quoted) > tol * max(1.0, abs(quoted)):
        raise ValueError(f"prose drift on {name}: the markdown says {quoted}, the cells give {got_}")

audit = [
    ("tables compared", 14, "disc_tables", 0),
    ("beta entries", 20, "beta_entries", 0),
    ("separation-table entries", 18, "rp_entries", 0),
    ("table/gas/bulb entries with a p/p0", 30, "x_entries_all", 0),
    ("q_cal cells", 101, "qcal_cells", 0),
    ("eq. (37) text-layer factor", 10.0, "eq37_textlayer_factor", 1e-6),
    ("sigma cells", 21, "sigma_cells", 0),
    ("beta cells", 23, "beta_cells", 0),
    ("fits done", 29, "fit_two_route_count", 0),
    ("implied mica area", 5747.4, "mica_s_implied_median", 1e-4),
    ("second b-prime outlier", 0.659, "s_second_worst_dev_pct", 1e-2),
    ("printed-claims rows", 36, "printed_claims_rows", 0),
]
for name, quoted, key, tol in audit:
    close(name, quoted, key, tol)
print(f"prose audit: {sum(len(v) for v in SAYS.values())} literals checked against their "
      f"claims with claims.says, {len(audit)} other numbers re-derived and matched.\n")
claims.report();   # the semicolon keeps the notebook from echoing the payload

# %% [markdown]
# ## 8. What pymrm adds
#
# **To Langmuir's isotherm, nothing.** Eqs. (9), (26), (29), (31) and (33) are
# closed forms in one variable. There is no grid on this page, no time step in
# space, no boundary condition and no transport, and most of the notebook would run
# with pymrm uninstalled, the same honest answer `A1.6`, `A1.1` and `J1.3` give.
# `newton` and `NumJac` do three specific jobs, and each was chosen because the
# cheap alternative is *wrong*, not slow.
#
# **(1) A fit that shares no algebra with the paper's fit.** Langmuir drew a line
# through $p/q$ against $p$. That is linear least squares on a transformed
# variable, and the transform is not innocent: it reweights the residuals towards
# the low-pressure points, where $q$ is small and $p/q$ is most sensitive. Profiling
# the amplitude out and root-finding the remaining stationarity condition on $q$
# itself is a different computation, and the break table carries the row that shows
# it matters, fitting
# by his linearisation instead moves the Freundlich-versus-Langmuir ratio, which is
# the page's discrimination metric. It also lets Freundlich be beaten **at its
# best** rather than at Langmuir's rendering of it, which is the only version of
# that contest worth reporting.
#
# **(2) The rate equation integrated instead of quoted.** Eq. (15) is printed
# wrongly. The way to establish that without asserting it is to *integrate eq. (14)*
# , which is what `march_eq14` does, with a `newton`-solved implicit step, and see
# which branch the march reaches. It reaches the corrected one to
# $1.8\times10^{-3}$ at 1600 steps, converging at first order in $\Delta t$ (the one
# axis that carries error here), and it cannot reach the printed one at all,
# because the printed one puts $t = -\infty$ at the start of the integration. That
# is a proof rather than a reading of the glyph.
#
# **(3) Root-finds where a sweep would have been wrong.** The relative pressure at
# which Case VI departs from Case I by 1 % is root-found; reading it off a 40-point
# log grid gives a value **6.3 % out**, and it is quoted in the same sentence as
# the number it is compared with. The pressures corresponding to the printed log-log
# exponents are root-found for the same reason.
#
# **(4) And the thing that is not pymrm at all, and is the most useful output on
# the page.** Langmuir's Case VI, eqs. (29) and (30), is the BET isotherm. That
# needed `sympy` and a change of variables the paper itself supplies through
# eq. (39), and it converts the paper's central conclusion from "the isotherm shows
# the film is one molecule thick" into "the isotherm could not have shown anything
# of the kind, and the film being one molecule thick is what $\beta$ shows". The
# separation between his two theories, over his entire data set, is at most
# **0.16 %**, against scatter **12 to 1100 times** larger.
#
# **What this page does not do that a reader might want.** It does not fit eq. (29)
# to Langmuir's data. That would be meaningless: below $p/p_0 = 1.6\times10^{-3}$, where
# every pressure he reached lies, the parameter
# $c$ is not identifiable, and any fit would return the Case I answer with an
# arbitrary $c$. Saying so is the finding.
#

# %% [markdown]
# ## 9. Reuse
#
# ### What you can take from here
#
# - **The isotherm itself**, in either parameterisation: $q = abp/(1+ap)$ with
#   $a$ in reciprocal bars, or $\theta = Kp/(1+Kp)$. If you need it inside a
#   reactor or a breakthrough model, `J1.5` (LDF breakthrough) is the page with the
#   transport; this one has none.
# - **The identity chain** $b \to b' \to N_0 \to \beta$. It is the cheapest
#   self-check that exists for any adsorption dataset reported as a monolayer
#   capacity: convert to molecules per unit area and divide by a close-packing
#   count from the liquid density. If $\beta > 1$ your surface area is wrong, your
#   isotherm is not Langmuir, or both.
# - **The discrimination protocol.** Refit *both* rivals, count the parameters,
#   print a zero-parameter and a one-parameter null beside the ratio, and check
#   whether the residuals *sweep*. A power law forced through a saturating curve
#   leaves a signature you can see (§6.6's right-hand panel) long before the RMS
#   tells you.
# - **The relative-pressure check, which is the transferable lesson.** Before
#   claiming that data support a monolayer model over a multilayer one, compute
#   $p/p_0$ and evaluate how far apart the two models actually are there. On this
#   paper the answer is 0.16 % at best, against 1–3 % scatter.
#
# ### What you must not take from here
#
# - **The $q_\mathrm{cal}$ agreement is a fit.** It is stated as one in the
#   notebook, in `meta.yaml`, in `models_entry.yaml`, in `README.md` and on the
#   case yaml. It is never reported as an agreement metric. If you see it quoted as
#   evidence for the Langmuir model, that is this page being misread.
# - **Tier 5, and only one column of it is measurement.** $q_\mathrm{obs}$ is
#   measured; $a$, $b$, $b'$, $\sigma$, $N_0$, $\beta$ and Tables XXI, XXII, XXV are
#   all derived by Langmuir from those measurements. Reproducing them is
#   reproduction, never validation.
# - **The mica numbers are uncorrected.** The paper says 10–30 % too high, and no
#   correction is applied anywhere on this page.
# - **Four printed defects, reported and none repaired**, §6.7 eq. (15)'s
#   logarithm; §6.2/§6.3 Table IV's $b = 58.3$; §6.3/§6.5 Table XX's $\beta = 0.36$
#   for methane on glass; §6.10 the pair "(100 bars) … 0.20". Where a correction is
#   stated it is labelled an inference.
# - **And, kept separate from them, one reporting convention.** The $\sigma$ column
#   of Tables XVIII–XX runs low against eq. (38): 18 of 21 printed values lie below
#   the computed one, worst 3.77 %. That is a real and non-random pattern, but it
#   is **not an error against Langmuir**, and §6.3 shows why, the three worst rows
#   are exactly the rows he prints to one or two significant figures, and truncating
#   at each cell's own printed granularity reproduces more cells than rounding does.
#   He discarded figures; he did not miscompute them. Filing it as a fifth defect
#   would overstate it.
# - **Nothing here is about platinum.** Journal pages 1393–1401 are irreversible
#   chemisorption with no isotherm, and Tables XXIII and XXIV were deliberately not
#   transcribed.
#
# ### Origins cited but not consulted
#
# Freundlich, H., *Kapillarchemie*, Leipzig (1909), the rival isotherm, read only
# through Langmuir's own restatement of it on journal page 1375 and through the
# $q_F$ column of Table VII. Eucken (1914); Bakker (1915); Haber (1914); Knudsen,
# *Ann. Physik* **31**, 205 (1910), the thermal-effusion correction; Meyer's
# *Kinetic Theory of Gases* (1899), the source of eq. (1); and Langmuir's own
# Parts I and II (*This Journal* **38**, 2221 (1916) and **39**, 1848 (1917)),
# which carry the monomolecular-film evidence this paper leans on. **None is on
# disk and none was consulted.** Nothing on this page derives from any of them.
#
# Brunauer, Emmett & Teller (1938) is on disk and is the subject of `J1.3`. The
# only thing taken from it here is the *form* of the BET isotherm, which `J1.3`
# documents as its eqs. (26)/(28); **no dataset of `J1.3` is loaded and no number
# of `J1.3` is retyped**, and the identity of §6.8 is verified against Langmuir's
# eq. (30) rather than against anything BET print. That the two theories are
# related is not this page's discovery and the 1938 paper says so itself, its
# Section II is headed *"Generalization of Langmuir's Theory to Multimolecular
# Adsorption"*, read here on a 300 ppi crop of journal page 311. What is new is
# that Langmuir's *own* Case VI closed form already **is** that generalisation,
# and that his data lie where it cannot be told apart from Case I.
#
# ### And one thing this page checked for somebody else
#
# Rawlings & Ekerdt, *Chemical Reactor Analysis and Design Fundamentals* (2nd edn),
# book page 443, state that the Danckwerts boundary conditions "were derived at
# least 45 years prior to Danckwerts in a classic paper by Langmuir [22]". **That
# paper is not this one.** Their reference [22] is *Langmuir, I., "The velocity of
# reactions in gases moving through heated vessels and the effect of convection and
# diffusion", J. Am. Chem. Soc.* **30**(11), 1742–1754 (1908), read from their own
# bibliography, and 1953 − 1908 = 45 exactly. The 1918 paper on disk contains no
# flow, no convection and no boundary condition of any kind: the words "convection"
# and "flow" do not occur in its forty-three pages, "diffusion" occurs three times
# and never in a transport equation, and every apparatus in it is a sealed static
# bulb. The claim bears on `A2.1`/`A2.2` and cannot be settled from this file;
# settling it needs JACS **30**(11) 1742, which is not on disk. Recorded on the
# case yaml as an acquisition target. **`A2.1` and `A2.2` were not touched.**
#
