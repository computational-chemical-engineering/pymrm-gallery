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
# title: "The Robeson upper bound, revisited: where the 2008 line moved, and what \"primarily the front factor\" can mean"
# description: "Robeson's 2008 Table 12 prints k and n for thirteen gas pairs, but only nine support a prior-versus-present comparison. On those nine the shift is decomposed into a front-factor part and a slope part at a stated selectivity, and the split is shown to depend on that selectivity, because log-log lines pivot. The paper's own ranking of which bounds moved most is reproduced by the selectivity gain at fixed permeability and by neither delta-log-k nor the permeability gain at fixed selectivity. Tables 13a/13b reproduce from Table 12, one row aside. The only second route, a covering line re-estimated from the paper's own near-bound points, agrees on the size of the shift at a stated reference and not on its split."
# categories: [sec:H, tier:T0, data:tier2, phase:gas, kind:correlation]
# date: 2026-08-08
# ---
#
# # The Robeson upper bound, revisited
#
# **Catalog ID:** `H1.8` · **Structures:** none (see below) · **Tier:** T0
#
# Every polymer membrane paper that plots selectivity against permeability draws
# one diagonal line on it. That line is
#
# $$
# P_i \;=\; k\,\alpha_{ij}^{\,n},
# \qquad \alpha_{ij}=\frac{P_i}{P_j},\qquad n<0,
# $$
#
# Robeson's empirical **upper bound**: for a given selectivity there is a
# permeability no homogeneous polymer film has beaten. In 1991 he drew it for
# nine gas pairs. In 2008 he redrew it against seventeen more years of
# literature and printed both sets of $(k,n)$ side by side, together with the
# 2008 data points that set the new line, together with a Knudsen-transition
# analysis computed from those same $(k,n)$.
#
# The paper's central testable claim is a **decomposition**: the bound moved,
# and *"the shift observed is primarily due to a change in the front factor,
# $k$, whereas the slope of the resultant upper bound relationship remains
# similar to the prior data correlations."* That is arithmetic on one printed
# table, and this page does it, on the **nine** gas pairs that actually support
# it (not the thirteen the table lists), at a **stated** reference selectivity,
# because a shift decomposed into "front factor" and "slope" parts is
# **meaningless until the reference is named**: two log–log lines with different
# slopes pivot, so moving the reference moves the split.
#
# Four results, all from printed numbers:
#
# 1. **Nine, not thirteen.** H2/O2 and He/O2 have a prior bound and no present
#    one; CO2/N2 and N2/CH4 have a present bound and no prior one. Any statement
#    about "the shift" lives on nine pairs. Robeson's own Conclusions name
#    exactly those nine.
# 2. **The claim holds where the paper says it holds, and the number to quote is
#    not $\Delta\log k$.** Because $|n|$ varies by a factor of more than 7 across
#    Table 12, the same $\Delta\log k$ buys wildly different amounts of separation.
#    Ranked by the **selectivity gain at fixed permeability**, the paper's three
#    *"significant"* pairs (He/CH4, He/CO2, He/H2) come out as the top three
#    with a clean gap. Ranked by $\Delta\log k$, or by the permeability gain at
#    fixed selectivity, they do not: CO2/CH4, which the paper calls *modest*,
#    interleaves with them.
# 3. **The two H2/CO2 bounds cross inside the paper's own data range.** For
#    $\alpha \gtrsim 37$ the 2008 "upper bound" lies *below* the 1994 one, and
#    Table 10 has **one** point out there, which is the whole of the "inside the
#    range" half of this claim and is stated as such wherever it appears. The
#    paper notices the slope change and warns its low-permeability end may be
#    skewed; it does not state the crossing.
# 4. **19 of the 117 points the paper tables as "close to the present empirical
#    upper bound" lie above it**, by up to a factor 3.5 in permeability. That is
#    not an error: the line is drawn *"by eye"*, as the paper says. It is what
#    "by eye" amounts to, measured.
#
# The page is a **correlation audit**, not a reactor model. There is no PDE, no
# grid and no pymrm operator in it, and the *What pymrm adds* section says so
# plainly.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.02"  # replaced by the release tag at publication
PYMRM = ">=3.3.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "H1.8-robeson-upper-bound"  # equals the page directory

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
# Imports. Nothing beyond the Colab defaults; pymrm is installed by the bootstrap
# but this page calls none of it (see "PyMRM implementation").
import textwrap
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.optimize import brentq
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import P_bound, alpha_bound, excursion, envelope_fit

np.set_printoptions(precision=6, suppress=False)
pd.set_option("display.width", 130)
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# **The source, precisely.** Lloyd M. Robeson, *"The upper bound revisited"*,
# Journal of Membrane Science **320**(1–2) 390–400 (2008),
# doi:`10.1016/j.memsci.2008.04.030`, Lehigh University, 1801 Mill Creek Road,
# Macungie, PA 18062, United States. Identified from the article's own running
# header (*"Journal of Membrane Science 320 (2008) 390–400"*), display title,
# by-line and affiliation block, all read on a 300 dpi render of p. 390; the
# PDF's own Title metadata carries the DOI.
#
# The file is **born-digital** (Elsevier / Acrobat Distiller 7.0), so the text
# layer is a faithful extraction rather than OCR. It was still not trusted:
# every numeric cell used on this page was read a second time from cropped
# 300 dpi renders at digit scale, and the two readings compared. They agree on
# all 382 **numeric cells**: 44 in Table 12 (its four `NA` rows aside), 234 in
# Tables 1–11 and 104 in Tables 13a/13b, counted in the loading cell below rather
# than asserted. These are cells, not numerals; the digit count is several times
# larger. The reason for the caution is specific to this paper: its tables
# carry **comma thousands separators inside cells** (the loading cell lists the
# Table 12 cells that carry them), which is exactly the class of glyph a naive
# read drops or mis-parses. A break row below drops one.
#
# **The 1991 and 1994 papers are the ORIGIN of the "prior" numbers and were not
# consulted.** Robeson (1991), *J. Membr. Sci.* **62**, 165, and Robeson,
# Burgoyne, Langsam, Savoca & Tien (1994), *Polymer* **35**, 4970, are not on
# disk. Every prior-bound value on this page was read from the 2008 paper, which
# restates them under its own heading *"Prior upper bound data (1,2)"*. That is
# the legitimate reprint route (`AGENTS.md`), and it has a limit this page
# respects throughout: **nothing here can establish that the 2008 restatement of
# the 1991 bounds is faithful to the 1991 paper.** Where a prior value carries a
# result, the page says it is read from the 2008 restatement.
#
# **What the paper does.** It re-reviews the literature since 1991 (*"the
# number of papers where the data were obtained approached 300"*), redraws the
# bound for eleven gas pairs (adding CO2/N2 and N2/CH4, which had no prior
# correlation), tabulates the points that set each new line (Tables 1–11),
# tabulates old and new $(k,n)$ together (Table 12), and then uses the new
# bounds to estimate where solution–diffusion transport must give way to Knudsen
# diffusion (Tables 13a and 13b). Its physical headline is that the big moves
# are almost all in **He-based pairs** and are driven by a family of
# perfluorinated commercial polymers (Nafion, Hyflon AD, Viton, Cytop, Teflon
# AF) that were barely in the 1991 dataset.
#
# **How the line is drawn matters, and the paper is honest about it.** From
# p. 391: *"It is noted that the 'upper bound' line is determined empirically
# ('by eye') as in the original reference [1] with sufficient data to establish
# a realistic bound over several decades of permeability."* So $(k,n)$ are **not
# a fit** to Tables 1–11 in any least-squares sense. Nothing on this page treats
# them as one, and result 4 above is the quantitative consequence.
#
# **Where it sits in the gallery.** [`H1.7`](../H1.7-solution-diffusion/) is the
# transport model underneath, permeability as solubility times diffusivity,
# same units, same barrer. This page supplies the empirical *ceiling* that model
# lives under; it does not derive it. [`H1.1`](../H1.1-sieverts-permeation/) and
# [`H1.9`](../H1.9-zeolite-membrane-maxwell-stefan-mixture/) are the metal and
# zeolite membranes to which the bound explicitly does **not** apply (p. 391:
# *"Heterogeneous membranes, surface modified membranes and molecular sieve
# membranes are not considered in the same class of polymeric materials"*).
# Methodologically the nearest sibling is
# [`B2.1`](../B2.1-voorhies-coking-law/): the origin of an empirical power law,
# tested only on what the source itself prints.
#
# **What this page deliberately does not do.** Twelve figures carry the actual
# data clouds; none is digitised, and none needs to be, because Tables 1–13
# print everything the arguments here use. Figure 12a/b, the correlation of
# $-1/n$ against gas-diameter difference, **cannot** be checked from this
# paper: it needs Lennard-Jones kinetic diameters, which the 2008 paper cites to
# Breck [65] and never prints. That is a checked claim, not an assumption: every
# occurrence of the string `diameter` in the full text layer is prose, a caption
# or eq. (1), and the only numeric diameters printed anywhere in the paper are
# PTMSP's pore size, *"in the range of 0.9–1.2 nm"*. Eq. (2), Freeman's
# theoretical prediction of $k$, is likewise uncheckable here: it needs
# solubility constants the paper does not tabulate. Both are out of scope and
# stay out.

# %% [markdown]
# ## The published model
#
# Equation numbers are the paper's.
#
# **The bound itself** (abstract, and p. 396 in the form Table 12 heads):
#
# $$P_i = k\,\alpha_{ij}^{\,n}$$
#
# with $P_i$ the permeability of the **fast** gas in barrers
# ($1\ \text{barrer} = 10^{-10}\ \mathrm{cm^3(STP)\,cm\,cm^{-2}\,s^{-1}\,cmHg^{-1}}$),
# $\alpha_{ij}=P_i/P_j$ the separation factor, $k$ the *"front factor"*, and $n$
# the slope of the log–log plot. Since $n<0$, in the plane the field actually
# plots, $\log\alpha$ against $\log P$, the bound is the falling diagonal
# $\log\alpha = (\log P - \log k)/n$, and *"below this line ... virtually all
# the experimental data points exist"*.
#
# **The shift, decomposed.** Write the prior bound $(k_0,n_0)$ and the present
# one $(k_1,n_1)$. At a selectivity $\alpha$ the bound permeability moves by
#
# $$
# \underbrace{\Delta \log_{10} P(\alpha)}_{\text{vertical shift}}
# =\underbrace{\Delta\log_{10}k}_{\text{front-factor part}}
# +\underbrace{\Delta n \,\log_{10}\alpha}_{\text{slope part}} .
# $$
#
# Three things follow immediately, and they organise the whole page.
#
# *First*, **the split depends on where you stand.** The front-factor part is
# the shift evaluated at $\alpha = 1$, where the membrane does not separate at
# all. That selectivity has no special standing in the data: it lies inside the
# tabulated range of 5 of the nine pairs and below every point of the others
# (Result 1 counts them). Reparameterise the same two lines about any other
# reference $\alpha_0$, writing $P = k'(\alpha/\alpha_0)^n$, and the
# "front-factor part" becomes $\Delta\log k + \Delta n \log\alpha_0$. So
# "primarily the front factor" is not a property of the pair of lines; it is a
# property of the pair of lines *and a chosen reference*. This page always
# states the reference.
#
# *Second*, the statement that **is** reference-free is that the two lines are
# nearly **parallel**: $|\Delta n|\log_{10}(\alpha_{\max}/\alpha_{\min})$, the
# amount the vertical shift varies across the selectivity range the data
# actually span, compared with the shift itself. That is the honest reading of
# Robeson's claim, and it is what gets measured below.
#
# *Third*, two non-parallel lines **cross**, at
# $\log_{10}\alpha^\star = -\Delta\log_{10}k/\Delta n$. If $\alpha^\star$ falls
# inside the data range, the "new" bound is *lower* than the old one over part
# of that range.
#
# **A second, orthogonal way to measure the same move.** A membrane engineer
# does not ask how much permeability the bound gained at fixed selectivity; they
# ask how much **selectivity** it gained at fixed permeability:
#
# $$
# \Delta \log_{10}\alpha(P)=\frac{\log_{10}P-\log_{10}k_1}{n_1}
# -\frac{\log_{10}P-\log_{10}k_0}{n_0}.
# $$
#
# These two measures are *not* monotonically related across pairs, because they
# differ by a factor $\sim 1/|n|$ and $|n|$ runs from 0.79 to 5.8 in Table 12.
# Which one reproduces Robeson's own verbal ranking is a result, not a
# convention.
#
# **The Knudsen transition** (eqs. 3–5, Tables 13a/13b). Knudsen diffusion
# selectivity is $D_i/D_j=(M_j/M_i)^{1/2}$ (eq. 3), and the solubility
# selectivity $k_s=S_i/S_j$ is taken invariant with permeability, so the
# separation factor at which solution–diffusion selectivity has fallen to the
# Knudsen value is
#
# $$\alpha_{ij}=k_s\,(M_j/M_i)^{1/2},$$
#
# and the permeability there is read off the bound as $P_i = k\,\alpha_{ij}^{n}$.
# $k_s$ itself comes from one of two correlations Freeman gives, eq. (4) from
# critical temperatures, $\ln(S_i/S_j)=0.016(T_{ci}-T_{cj})$, giving Table 13a;
# eq. (5) from Lennard-Jones temperatures,
# $\ln(S_i/S_j)=0.023[(\varepsilon_i/k)-(\varepsilon_j/k)]$, giving Table 13b.
# Neither temperature set is printed, so **$k_s$ is an input here, not a
# checkable quantity**. The last column, however, is entirely checkable against
# Table 12, and reproducing all 26 of its rows is this page's strongest test of
# its own transcription (V1 below).

# %% [markdown]
# ## Parameters and assumptions
#
# | symbol | meaning | source | status here |
# |---|---|---|---|
# | $k$ | front factor, barrers | Table 12, both blocks | transcribed; prior block via the 2008 restatement |
# | $n$ | log–log slope, dimensionless, $<0$ | Table 12, both blocks | transcribed |
# | $P_i,\alpha_{ij}$ | the 117 near-bound points | Tables 1–11 | transcribed; **compiled by Robeson from ~300 papers, none consulted here** |
# | $k_s$ | solubility selectivity | Tables 13a/13b, col. 1 | **input**: eqs. (4)/(5) need temperatures the paper does not print |
# | $(M_j/M_i)^{1/2}$ | Knudsen selectivity | Tables 13a/13b, col. 2 | checkable, and checked, from integer molecular weights |
# | $\alpha_0$ | reference selectivity for the decomposition | **not printed anywhere** | *chosen here*, and always stated |
#
# Assumptions, each one a decision that could have gone otherwise:
#
# 1. **The reference selectivity $\alpha_0$ is the geometric mean of the
#    printed near-bound $\alpha$ for that pair.** It is data-driven, it is
#    deterministic, and it comes from the paper's own tables, but it is a
#    choice, it is not Robeson's, and the page reports the whole $\alpha$ range
#    beside it so the choice can be audited. Robeson's implicit reference is
#    $\alpha_0=1$, which lies inside the data of only 5 of the nine pairs.
# 2. **Molecular weights are integers** (He 4, H2 2, O2 32, N2 28, CO2 44,
#    CH4 16). This is not a convenience: it is *identified from the paper*.
#    For H2/N2 the table prints **3.742**, which is $\sqrt{28/2}$; IUPAC
#    standard atomic weights do not reproduce it (V1 prints both). All eleven
#    pairs follow.
# 3. **The printed rounding band is taken literally**: half a unit in the last
#    printed digit position, so a $k$ printed to the unit carries $\pm 0.5$ and
#    an $n$ printed to three decimals carries $\pm0.0005$. Where a conclusion
#    depends on the band, the page also reports the **generous** band that
#    treats trailing zeros as placeholders (four significant figures on $k$ and
#    on $P$), and states which reading the conclusion survives.
# 4. **A point lies above the bound when $r=P/(k\alpha^{n})>1$.** With $n<0$
#    this is the correct direction and it is worth deriving once: a point beats
#    the bound if its selectivity exceeds the bound's at its permeability,
#    $\log\alpha>(\log P-\log k)/n$; multiplying by the *negative* $n$ flips the
#    inequality to $\log k+n\log\alpha<\log P$, i.e. $r>1$. Getting this
#    backwards counts the points below the bound instead, which is break row
#    (I) below.
# 5. **Nothing is refit as a substitute for the printed bound.** The line is
#    drawn by eye; a least-squares line through Tables 1–11 is a *different
#    object* and is labelled as one wherever it appears.

# %% [markdown]
# ## The data
#
# Four datasets, all transcription of printed tables and printed sentences.
# **Tier 2** throughout: the paper's own tabulated values, no figure digitised,
# no number inferred. Each has a provenance sidecar, and what the sidecars
# record about the rows is collected here.
#
# - `robeson-2008-table12-upper-bounds`: Table 12, a prior and a present row
#   for each of 13 gas pairs. **The `NA` structure is load-bearing**: H2/O2 and He/O2 have
#   no present bound, CO2/N2 and N2/CH4 no prior one, so only nine pairs support
#   a comparison. The `prior` rows are the 1991 and 1994 bounds as the 2008
#   paper restates them; the originals were not consulted. The prior block
#   prints `He/CO2` without the final subscript, a typographic slip normalised
#   in the CSV. The `*_as_printed` columns keep the printed glyphs (commas,
#   trailing zeros) and are loaded as strings, because the printed digits set
#   the rounding bands of V1.
# - `robeson-2008-near-bound-points`: Tables 1–11, one row per printed row, in
#   printed order (not always monotone in permeability). Two rows are flagged
#   **by the paper itself** and used as printed: an O2/N2 point whose position
#   above the bound the paper questions, and an N2/CH4 point whose CH4
#   permeability the paper corrects in prose while its Table 9 keeps the
#   unrepaired $\alpha$. Result 4 locates both from the numbers in the paper's
#   own sentences and checks which value the table carries. Polymer
#   names are verbatim, including Table 2's `Polypyrrole` where Tables 1, 5 and
#   6 print `Polypyrrolone` for the same material, and Table 9's
#   `Poly(pyrrolone)(6FDA-TAB)` where Tables 5 and 6 print `Polypyrrolone
#   (6FDA-TAB)`. The `ref_bracket` column counts the cited references below.
# - `robeson-2008-table13-knudsen-transition`: Tables 13a/13b, 26 rows. The
#   `dataset` column is this page's encoding of the printed "(prior data)"
#   label. $k_s$ is not checkable from this paper. The table is an independent
#   check on this page's transcription of Table 12, which is what V1 tests; it
#   shares Table 12 as an input, so it is not a second route to Robeson's
#   numbers (V2 says why). Its product column is read a second time as strings,
#   because a float drops printed trailing zeros (`0.270`) and with them the
#   rounding band.
# - `robeson-2008-printed-claims`: the prose claims tested here, verbatim. Its
#   `shift_class` column is this page's one-word encoding of the Conclusions'
#   wording. Two quoted sentences are printed defects, kept as printed: the
#   Conclusions name He/CO2 twice and drop He/CH4, and the author repairs the
#   N2/CH4 permeability in prose and not in Table 9.
#
# No dataset from another gallery page is loaded, so the cross-page
# reconciliation rule has nothing to bite on here.
#
# The provenance depth of the second dataset bounds every result drawn from it:
# **these are not Robeson's measurements.** Each row cites the primary paper
# that reported it; the rows cite 44 distinct references (counted below), none
# of which was consulted. This page can show that a tabulated point sits above
# the line Robeson drew; it cannot show that the point is right.

# %%
# `dtype=str` + `keep_default_na=False` so the as-printed columns keep their trailing
# zeros and their literal "NA" - both are load-bearing below.
T12 = load_data("robeson-2008-table12-upper-bounds.csv", page=PAGE,
                dtype=str, keep_default_na=False)
for c in ("k_barrer", "n"):
    T12[c] = pd.to_numeric(T12[c].replace("", np.nan))
PTS = load_data("robeson-2008-near-bound-points.csv", page=PAGE)
T13 = load_data("robeson-2008-table13-knudsen-transition.csv", page=PAGE)
# the product column as printed, for its rounding band: "0.270" read as a float loses its zero
T13["ks_D_as_printed"] = load_data("robeson-2008-table13-knudsen-transition.csv", page=PAGE,
                                   dtype=str).ks_D_ratio
CLM = load_data("robeson-2008-printed-claims.csv", page=PAGE)

print(cite_data(load_meta("robeson-2008-table12-upper-bounds.csv", page=PAGE)))
print(f"Table 12: {len(T12)} rows | Tables 1-11: {len(PTS)} points | "
      f"Tables 13a/13b: {len(T13)} rows | printed claims: {len(CLM)}")
# The transcription burden, counted rather than asserted: every one of these was read twice,
# once from the text layer and once from a digit-scale crop of a 300 dpi render.
# These are numeric CELLS, not numerals - the digit count is several times larger.
N_NUM_T12 = int(2 * T12.k_barrer.notna().sum())
N_NUM_PTS = 2 * len(PTS)
N_NUM_T13 = 4 * len(T13)
N_NUMERALS = N_NUM_T12 + N_NUM_PTS + N_NUM_T13
print(f"Numeric cells transcribed: {N_NUM_T12} (Table 12, excluding the 4 NA rows) + {N_NUM_PTS} "
      f"(Tables 1-11) + {N_NUM_T13} (Tables 13a/13b) = {N_NUMERALS}")
print("Table 12 cells printed with thousands commas:",
      ", ".join(f"{r.dataset} {r.gas_pair} k = {r.k_as_printed}"
                for r in T12.itertuples() if "," in r.k_as_printed))

prior = T12[T12.dataset == "prior"].set_index("gas_pair")
present = T12[T12.dataset == "present"].set_index("gas_pair")
PAIRS_ALL = list(prior.index)
PAIRS = [g for g in PAIRS_ALL if np.isfinite(prior.k_barrer[g]) and np.isfinite(present.k_barrer[g])]
N_PAIRS_TABLE, N_PAIRS_COMPARABLE = len(PAIRS_ALL), len(PAIRS)

print(f"\nTable 12 lists {N_PAIRS_TABLE} gas pairs; {N_PAIRS_COMPARABLE} carry BOTH a prior and a present bound.")
print("  present bound NA (prior only):",
      ", ".join(g for g in PAIRS_ALL if not np.isfinite(present.k_barrer[g])))
print("  prior bound NA (present only):",
      ", ".join(g for g in PAIRS_ALL if not np.isfinite(prior.k_barrer[g])))
print("  comparable:", ", ".join(PAIRS))
print("\nTables 1-11 exist for exactly the 11 pairs with a PRESENT bound:")
print("  ", ", ".join(sorted(PTS.gas_pair.unique())))
print("  pairs with a present bound but no table:",
      sorted(set(present.index[np.isfinite(present.k_barrer)]) - set(PTS.gas_pair.unique())) or "none")

# %%
# The paper's own three-way classification of how far each bound moved, joined onto Table 12.
CLASS = (CLM[CLM.shift_class.notna() & (CLM.gas_pair.notna())]
         .drop_duplicates("gas_pair").set_index("gas_pair").shift_class.to_dict())
print("Robeson's Conclusions classify these pairs, and only these:")
for g in PAIRS:
    print(f"  {g:9s} {CLASS[g]}")
missing = sorted(set(PAIRS) - set(CLASS))
extra = sorted(set(CLASS) - set(PAIRS))
print(f"\nclassified but not comparable: {extra or 'none'}"
      f"  |  comparable but unclassified: {missing or 'none'}")
assert not missing and not extra, "the Conclusions' pair list must be exactly the nine comparable pairs"
print("\nSo the paper itself confines its shift discussion to the nine pairs the NA structure allows.")
print("A page that compared thirteen would be comparing four pairs against nothing.")

# %%
# H2/O2 and He/O2: prior bound printed, present bound NA, no section, no table, no figure.
# A negative claim, so here is the search. The text layer is born-digital; extraction inserts a
# space before the slash for H2 but not for He, so both spellings are counted.
print("Where do H2/O2 and He/O2 appear in the paper?")
print("  Searched: the complete `pdftotext -layout` extraction of all 11 pages, for both the")
print("  hyphen-free spellings and the space-before-slash form the extractor produces.")
print("  Result: they occur only as cells of Table 12, in its prior block and in its present block")
print("  (as 'NA'). No section heading, no data table, no figure, and no sentence anywhere in the")
print("  paper mentions either pair. The paper does not say why they were dropped.")
print("  Consistency check available here: Tables 1-11 cover 11 pairs and neither of these two.")
assert set(PTS.gas_pair.unique()) == set(present.index[np.isfinite(present.k_barrer)]), \
    "every present bound must have a data table and vice versa"
print("  -> verified: {} data tables, {} present bounds, same set.".format(
    PTS.gas_pair.nunique(), int(np.isfinite(present.k_barrer).sum())))


# %% [markdown]
# ## PyMRM implementation
#
# **There is none, and that is the honest answer.**
#
# The Robeson upper bound is a two-parameter algebraic relation between two
# scalars. It has no spatial coordinate, no time coordinate, no conservation
# law, no boundary condition and no unknown field. Nothing in pymrm's operator
# library (`construct_grad`, `construct_div`, `construct_convflux_upwind`,
# `NumJac`, `newton`) has anything to act on. Building a grid here would be
# decoration, and the gallery's structure taxonomy (`S1`–`S13` in
# [`docs/taxonomy.md`](../../docs/taxonomy.md)) has no code that fits: the page's
# `structures` field is **empty**, deliberately, and its `pymrm_api` list is
# empty too.
#
# What the page uses is numpy for the algebra, one `scipy.optimize.brentq` call
# per threshold (thresholds are root-found, never swept), and an exact
# enumeration of the active sets of one small constrained least-squares
# problem for the envelope estimator in *Validation*. The four functions that
# are reusable outside this page, `P_bound`, `alpha_bound`, `excursion` and
# `envelope_fit`, live in `model.py`, with docstrings and a test; the bootstrap
# cell installs pymrm only because every gallery page pins it.
#
# Where the bound *does* meet pymrm is one layer up, and this page does not
# cross that layer: $(P_i, \alpha_{ij})$ is the closure a module-scale membrane
# model consumes. [`H1.7`](../H1.7-solution-diffusion/) builds the transport side
# of that closure from solubility and diffusivity in the same units; the bound
# tells you which $(P,\alpha)$ pairs that model is entitled to be handed. The
# two pages are complementary and neither computes the other.
#
# The rest of this page is therefore arithmetic on printed tables, which is
# exactly what auditing an empirical correlation is.

# %%
# The bound and its inverse (model.py) are written once and used everywhere below, so
# a sign error cannot hide in one branch. The round trip is an identity, registered as a
# claim; break row (N) writes the inverse with n instead of 1/n.
def roundtrip_max(pre, pairs, inv=alpha_bound):
    """Worst relative error of inv(P_bound(alpha)) over four selectivities and the pairs."""
    a = np.array([0.2, 1.0, 7.5, 900.0])
    return max(float(np.max(np.abs(inv(P_bound(a, pre.k_barrer[g], pre.n[g]),
                                       pre.k_barrer[g], pre.n[g]) / a - 1))) for g in pairs)


ROUNDTRIP_MAX = roundtrip_max(present, PAIRS)
print(f"P_bound / alpha_bound round trip, worst over the {len(PAIRS)} pairs: below 1e-12 (here {ROUNDTRIP_MAX:.0e})")

# a point exactly ON the bound has r = 1, by construction; a point at twice the bound
# permeability has r = 2. Sanity, and the sign convention made explicit.
kk, nn = present.k_barrer["He/H2"], present.n["He/H2"]
print(f"sign check on He/H2 (k={kk:,.0f}, n={nn}):")
for fac in (0.5, 1.0, 2.0):
    print(f"  a point at {fac:>3}x the bound permeability at alpha=2.5 has r = "
          f"{float(excursion(fac * P_bound(2.5, kk, nn), 2.5, kk, nn)):.3f}")

# %% [markdown]
# ## Results
#
# ### Result 1: The shift, decomposed at a stated selectivity
#
# For each of the nine comparable pairs: $\Delta\log_{10}k$, $\Delta n$, the
# range of selectivity its own near-bound points span, the vertical shift at the
# geometric-mean selectivity of those points, and the **front-factor share**
# $\Delta\log k/\Delta\log P(\alpha_0)$, the fraction of the shift the front
# factor accounts for *at that reference*.
#
# The column that carries the reference-free version of Robeson's claim is
# **`spread/shift`**: how much the vertical shift varies across the pair's own
# selectivity range, relative to the shift itself. Small means the two lines are
# parallel and "primarily the front factor" is safe wherever you stand; large
# means the split is an artefact of where you stood. The cell also counts the
# pairs whose own data reach $\alpha = 1$, the reference at which
# $\Delta\log k$ is the shift.

# %%
rows = []
for g in PAIRS:
    k0, n0 = prior.k_barrer[g], prior.n[g]
    k1, n1 = present.k_barrer[g], present.n[g]
    sub = PTS[PTS.gas_pair == g]
    a = sub.alpha.to_numpy(float)
    a_lo, a_hi = a.min(), a.max()
    a_gm = float(np.exp(np.mean(np.log(a))))
    dlk = np.log10(k1) - np.log10(k0)
    dn = n1 - n0
    shift = lambda al: dlk + dn * np.log10(al)          # noqa: E731 - one line, used thrice
    spread = abs(dn) * np.log10(a_hi / a_lo)
    rows.append(dict(pair=g, verdict=CLASS[g], dlogk=dlk, dn=dn, n_pres=n1,
                     a_lo=a_lo, a_hi=a_hi, a_gm=a_gm,
                     shift_gm=shift(a_gm), share=dlk / shift(a_gm),
                     spread=spread, spread_over_shift=spread / abs(shift(a_gm)),
                     shift_lo=shift(a_lo), shift_hi=shift(a_hi)))
SH = pd.DataFrame(rows).set_index("pair")

print("\n".join(
    ["Shift of the upper bound, nine comparable pairs.  All logs base 10; 'decades' = log10 units.",
     f"{'pair':9s} {'paper':12s} {'dlog k':>7s} {'dn':>8s} {'alpha range':>18s} "
     f"{'a_gm':>7s} {'shift@a_gm':>10s} {'k-share':>8s} {'spread/shift':>12s}"]
    + [f"{g:9s} {r.verdict:12s} {r.dlogk:7.4f} {r.dn:+8.4f} "
       f"{r.a_lo:8.3g} - {r.a_hi:<7.4g} {r.a_gm:7.3g} {r.shift_gm:10.4f} "
       f"{r.share:8.3f} {r.spread_over_shift:12.2f}" for g, r in SH.iterrows()]))

SHARE_MIN, SHARE_MAX = float(SH.share.min()), float(SH.share.max())
SHARE_MED = float(SH.share.median())
PARALLEL_N = int((SH.spread_over_shift < 0.20).sum())
print(f"\nFront-factor share at the geometric-mean selectivity: median {SHARE_MED:.2f}, "
      f"range {SHARE_MIN:.2f} to {SHARE_MAX:.2f}.")
print(f"Two lines are near-parallel (spread/shift < 0.20) for {PARALLEL_N} of the "
      f"{len(SH)} pairs. For those, 'primarily the front")
print("factor' is a statement about the lines, not about a chosen reference. For the rest it is not:")
for g, r in SH[SH.spread_over_shift >= 0.20].iterrows():
    print(f"  {g:9s} the vertical shift runs {r.shift_lo:+.3f} -> {r.shift_hi:+.3f} decades across "
          f"alpha = {r.a_lo:.3g} to {r.a_hi:.4g}")
print("\nAnd note what the reference does to the headline pair. He/H2's front factor moves")
print(f"{prior.k_barrer['He/H2']:,.0f} -> {present.k_barrer['He/H2']:,.0f} barrers, a factor "
      f"{present.k_barrer['He/H2']/prior.k_barrer['He/H2']:.1f}: dlog k = {SH.dlogk['He/H2']:.4f}.")
print(f"Its slope moves {prior.n['He/H2']} -> {present.n['He/H2']}, i.e. {100*abs(SH.dn['He/H2']/prior.n['He/H2']):.1f} % of itself.")
print(f"Over the selectivity range its OWN data span ({SH.a_lo['He/H2']:.2f} to {SH.a_hi['He/H2']:.2f}), the "
      f"{SH.shift_gm['He/H2']:.2f}-decade shift")
print(f"varies by {100*SH.spread['He/H2']/abs(SH.shift_gm['He/H2']):.1f} % of itself, and the slope term "
      f"dn log10(alpha) is at most {100*abs(SH.dn['He/H2'])*np.log10(SH.a_hi['He/H2'])/SH.shift_gm['He/H2']:.1f} %"
      " of it. There the claim is airtight.")

# Robeson's own reference, alpha = 1, measured against the data rather than assumed outside it.
IN1 = SH[(SH.a_lo < 1) & (SH.a_hi > 1)]
print(f"\nalpha = 1, where dlog k IS the shift, lies inside the tabulated range of {len(IN1)} of the "
      f"{len(SH)} pairs:")
print("\n".join(f"  {g:9s} {int((PTS[PTS.gas_pair == g].alpha < 1).sum())} of "
                f"{int((PTS.gas_pair == g).sum())} points below alpha = 1 (alpha_min {r.a_lo:.3g})"
                for g, r in IN1.iterrows()))
print("and outside it, below every tabulated point, for " + ", ".join(SH.index[SH.a_lo > 1]) + ".")

# %% [markdown]
# ### Result 2: Where the two bounds cross
#
# $\Delta\log P(\alpha)=0$ at
# $\log_{10}\alpha^\star=-\Delta\log_{10}k/\Delta n$. Computed two ways, the
# closed form, and `brentq` root-finding the shift function itself, and then
# asked whether $\alpha^\star$ falls inside the range of the paper's own
# tabulated points for that pair.
#
# Two separate claims live here and they are of very different strengths. That
# the two H2/CO2 lines cross at $\alpha^\star\approx 37$ is a property of four
# printed numbers in Table 12. That the crossing happens **inside the range
# Table 10 spans** is a property of a *single tabulated point*: the cell below
# names it, prints the gap it sits across, and says what is left without it. A
# break row further down mis-reads that one decimal point and takes the result to
# zero.

# %%
_lines = [f"{'pair':9s} {'alpha* (closed form)':>21s} {'alpha* (brentq)':>16s} "
          f"{'|diff|':>9s} {'inside its own alpha range?':>29s}"]
cross_rows = []
for g in PAIRS:
    k0, n0 = prior.k_barrer[g], prior.n[g]
    k1, n1 = present.k_barrer[g], present.n[g]
    dlk, dn = np.log10(k1) - np.log10(k0), n1 - n0
    a_star_cf = 10.0 ** (-dlk / dn)
    f = lambda la: dlk + dn * la                        # noqa: E731 - shift in log-alpha
    lo, hi = np.log10(SH.a_lo[g]), np.log10(SH.a_hi[g])
    inside = f(lo) * f(hi) < 0
    # the same bracket as the registered claim: [-3, 6] decades holds every printed selectivity
    a_star_rf = 10.0 ** brentq(f, -3.0, 6.0, xtol=1e-14, rtol=1e-15) if inside else np.nan
    diff = abs(a_star_rf / a_star_cf - 1) if inside else np.nan
    cross_rows.append(dict(pair=g, a_star=a_star_cf, a_star_rf=a_star_rf, inside=inside))
    _lines.append(f"{g:9s} {a_star_cf:21.4g} {a_star_rf if inside else float('nan'):16.6f} "
                  f"{'< 1e-12' if inside and diff < 1e-12 else 'n/a':>9s} "
                  f"{('YES - inside ' + f'{SH.a_lo[g]:.3g}-{SH.a_hi[g]:.4g}') if inside else 'no':>29s}")
print("\n".join(_lines))
CROSS = pd.DataFrame(cross_rows).set_index("pair")
N_CROSSING_INSIDE = int(CROSS.inside.sum())
G_CROSS = CROSS.index[CROSS.inside][0]
ALPHA_STAR = float(CROSS.a_star[G_CROSS]); ALPHA_STAR_RF = float(CROSS.a_star_rf[G_CROSS])
ALPHA_STAR_AGREE = abs(ALPHA_STAR_RF / ALPHA_STAR - 1)
print(f"\n{N_CROSSING_INSIDE} of {len(PAIRS)} pairs cross inside their own tabulated selectivity range: "
      f"{G_CROSS}, at alpha* = {ALPHA_STAR:.2f}")
print(f"(closed form and brentq agree below 1e-12 relative (here {ALPHA_STAR_AGREE:.0e}) - the root-find exists to")
print(" catch a sign or reciprocal slip in the closed form, not to add precision.)")

k0, n0 = prior.k_barrer[G_CROSS], prior.n[G_CROSS]
k1, n1 = present.k_barrer[G_CROSS], present.n[G_CROSS]
sub = PTS[PTS.gas_pair == G_CROSS].sort_values("alpha")
above = sub[sub.alpha > ALPHA_STAR]
print(f"\nAbove alpha* = {ALPHA_STAR:.2f} the PRESENT {G_CROSS} bound lies BELOW the prior one.")
print(f"{len(above)} of the {len(sub)} points Table 10 tabulates sit in that region:")
for _, r in above.iterrows():
    Pp, Pq = float(P_bound(r.alpha, k1, n1)), float(P_bound(r.alpha, k0, n0))
    print(f"  alpha = {r.alpha:7.4g}  point P = {r.P_fast_barrer:9.4g} | present bound "
          f"{Pp:9.4g} | prior bound {Pq:9.4g}  (present is {100*(Pp/Pq-1):+.1f} % of prior)")
PRESENT_BELOW_PRIOR_PCT = float(100 * (P_bound(sub.alpha.max(), k1, n1)
                                       / P_bound(sub.alpha.max(), k0, n0) - 1))

# How much does this verdict rest on? Exactly one tabulated point, and the page must say so.
_below = sub[sub.alpha < ALPHA_STAR]
print(f"\nAND HOW THIN IS THAT. 'Inside the paper's own tabulated range' rests on "
      f"{len(above)} point{'s' if len(above) != 1 else ''} of")
print(f"the {len(sub)}: alpha* = {ALPHA_STAR:.2f} sits between the tabulated alpha = "
      f"{_below.alpha.max():.4g} and alpha = {above.alpha.min():.4g},")
print(f"a gap of a factor {above.alpha.min()/_below.alpha.max():.1f} with nothing in it. Delete the "
      f"{above.polymer.iloc[0][:38]}")
print(f"point at alpha = {above.alpha.min():.4g} - or mis-read its decimal point, which is a break "
      f"row below - and")
print(f"alpha* falls outside Table 10's range and Result 2 has nothing to stand on. The crossing")
print(f"itself would survive: it is a property of the two printed lines and moves only with them.")
print(f"What would not survive is the claim that it happens where the paper has data. Any use of")
print(f"this result has to carry that sentence with it.")
print(f"\nThe paper anticipates the cause and says so on p. 396 - the H2/CO2 shift is")
print(f'"primarily a slight slope change" and "the limited number of data points at the lower')
print(f'permeability area of the dataset may have skewed the slope versus the original')
print(f'correlation [2]". Searched the full text for any statement that the two lines CROSS or')
print(f"that the new bound falls below the old one anywhere: none - the words cross/crossing/")
print(f"intersect appear only inside the polymer name 'crosslinked' and in reference titles.")
print(f"So the crossing is a consequence of what the paper says, not a contradiction of it -")
print(f"but it is a consequence the paper does not draw, and it is the one place in the table")
print(f"where 'upper bound' is not monotone in time.")

# %% [markdown]
# ### Result 3: Vertical or horizontal? Which measure reproduces the paper's own ranking
#
# The paper sorts the nine pairs into *modest* (5), *minor* (1) and
# *significant* (3). Nothing forces that verbal sorting to agree with any
# particular arithmetic measure, so it is a real test of which measure Robeson
# was actually reasoning with.

# %%
rank_rows = []
for g in PAIRS:
    k0, n0 = prior.k_barrer[g], prior.n[g]
    k1, n1 = present.k_barrer[g], present.n[g]
    P = PTS[PTS.gas_pair == g].P_fast_barrer.to_numpy(float)
    P_gm = float(np.exp(np.mean(np.log(P))))
    # horizontal shift, closed form ...
    dlogA = ((np.log10(P_gm) - np.log10(k1)) / n1) - ((np.log10(P_gm) - np.log10(k0)) / n0)
    # ... and by root-finding each bound's alpha at P_gm independently (no shared algebra).
    # k*alpha^n - P is monotone in alpha for n < 0, and [1e-6, 1e9] spans every separation
    # factor in the paper by orders of magnitude; brentq raises if the root is outside it.
    a1 = brentq(lambda a: P_bound(a, k1, n1) - P_gm, 1e-6, 1e9, xtol=1e-14, rtol=8.9e-16)
    a0 = brentq(lambda a: P_bound(a, k0, n0) - P_gm, 1e-6, 1e9, xtol=1e-14, rtol=8.9e-16)
    dlogA_rf = np.log10(a1) - np.log10(a0)
    rank_rows.append(dict(pair=g, verdict=CLASS[g], dlogk=SH.dlogk[g], dlogP=SH.shift_gm[g],
                          dlogA=dlogA, dlogA_rf=dlogA_rf, alpha_gain=10 ** dlogA,
                          P_gm=P_gm, n_pres=n1))
RK = pd.DataFrame(rank_rows).set_index("pair").sort_values("dlogA", ascending=False)
DLOGA_ROOTFIND_MAX = float(np.max(np.abs(RK.dlogA - RK.dlogA_rf)))

print("\n".join(
    ["Ranked by SELECTIVITY gain at the pair's own geometric-mean permeability:",
     f"{'pair':9s} {'paper':12s} {'n_pres':>7s} {'dlog k':>7s} {'dlogP@a_gm':>10s} "
     f"{'dlogA@P_gm':>10s} {'alpha x':>8s}"]
    + [f"{g:9s} {r.verdict:12s} {r.n_pres:7.3f} {r.dlogk:7.4f} {r.dlogP:10.4f} "
       f"{r.dlogA:10.4f} {r.alpha_gain:8.3f}" for g, r in RK.iterrows()]))
print(f"\n(closed form vs independent root-find of each bound at P_gm: max |difference| below "
      f"1e-12 decades (here {DLOGA_ROOTFIND_MAX:.0e}))")

SIG = [g for g in PAIRS if CLASS[g] == "significant"]
def top3(col):
    return list(RK.sort_values(col, ascending=False).index[:3])
TOP3_A, TOP3_K, TOP3_P = top3("dlogA"), top3("dlogk"), top3("dlogP")
HITS_A, HITS_K, HITS_P = (len(set(TOP3_A) & set(SIG)), len(set(TOP3_K) & set(SIG)),
                          len(set(TOP3_P) & set(SIG)))
gap_A = float(RK.dlogA.iloc[2] - RK.dlogA.iloc[3])
print(f"\nThe paper's three 'significant' pairs: {', '.join(SIG)}")
print(f"  top 3 by selectivity gain dlogA : {', '.join(TOP3_A):40s} {HITS_A} of 3")
print(f"  top 3 by front factor  dlog k   : {', '.join(TOP3_K):40s} {HITS_K} of 3")
print(f"  top 3 by permeability  dlogP    : {', '.join(TOP3_P):40s} {HITS_P} of 3")
print(f"\nOnly the horizontal measure reproduces it, and it does so with a clean gap: "
      f"{RK.dlogA.iloc[2]:.4f} vs {RK.dlogA.iloc[3]:.4f}")
print(f"decades between third and fourth ({gap_A:.4f} decades, {100*gap_A/RK.dlogA.iloc[3]:.0f} % of the fourth).")
print(f"\nThe reason is arithmetic, not taste. CO2/CH4 and He/CH4 have almost the same front-factor")
print(f"shift ({SH.dlogk['CO2/CH4']:.3f} vs {SH.dlogk['He/CH4']:.3f} decades) but slopes differing by "
      f"{abs(present.n['CO2/CH4']/present.n['He/CH4']):.1f}x")
print(f"({present.n['CO2/CH4']} vs {present.n['He/CH4']}). Divide the vertical shift by |n| and the same")
print(f"front-factor move buys CO2/CH4 a factor {RK.alpha_gain['CO2/CH4']:.2f} in selectivity and He/CH4 a factor "
      f"{RK.alpha_gain['He/CH4']:.2f}.")
print(f"So both of Robeson's statements are true at once: the shift IS mostly in k, AND k is the")
print(f"wrong number to rank pairs by. The abstract's sentence is about the SHAPE of the move;")
print(f"the Conclusions' ranking is about its SIZE, and they are different quantities.")

# ---- the evaluation point is a choice here too, so sweep it ---------------------------
# dlogA(P) varies with P for exactly the reason dlogP(alpha) varies with alpha: the two lines
# are not parallel. The RANKING and the GAP are two different questions and get two answers.
Pmin_all = max(PTS[PTS.gas_pair == g].P_fast_barrer.min() for g in PAIRS)
Pmax_all = min(PTS[PTS.gas_pair == g].P_fast_barrer.max() for g in PAIRS)
P_ALL_GM = float(np.exp(np.mean(np.log(PTS.P_fast_barrer.to_numpy(float)))))
def rank_at(Pfun, full=False):
    d = sorted(((g, ((np.log10(Pfun(g)) - np.log10(present.k_barrer[g])) / present.n[g]
                     - (np.log10(Pfun(g)) - np.log10(prior.k_barrer[g])) / prior.n[g]))
                for g in PAIRS), key=lambda z: -z[1])
    out = len(set(g for g, _ in d[:3]) & set(SIG)), d[2][1] - d[3][1]
    return out + (d[2][0], d[3][0]) if full else out
own = lambda f: (lambda g: float(f(PTS[PTS.gas_pair == g].P_fast_barrer.to_numpy(float))))  # noqa
EVAL = [("each pair's own P_min", own(np.min)),
        ("each pair's own P_gm (used above)", own(lambda P: np.exp(np.mean(np.log(P))))),
        ("each pair's own P_median", own(np.median)),
        ("each pair's own P_max", own(np.max)),
        (f"all 117 points' geometric mean, {P_ALL_GM:.0f} barrer", lambda g: P_ALL_GM)]
EVAL += [(f"common P = {c:g} barrer", (lambda g, c=c: float(c)))
         for c in (100.0, 1000.0) ]
print(f"\nThe evaluation point is a CHOICE here just as alpha_0 was in Result 1 - dlogA varies with")
print(f"P because the lines are not parallel. Every reference below lies inside the permeability")
print(f"range all nine pairs share ({Pmin_all:.0f} to {Pmax_all:.0f} barrers) or is each pair's own:")
_h, _g, _lines = [], [], [f"\n{'evaluation point':44s} {'top-3 hits':>11s} {'gap, decades':>13s}"]
for lab, f in EVAL:
    h, gp = rank_at(f)
    _h.append(h); _g.append(gp)
    _lines.append(f"{lab:44s} {h:>7d} of 3 {gp:13.4f}")
print("\n".join(_lines))
EVAL_HITS_MIN = float(min(_h)); EVAL_GAP_MIN = float(min(_g)); EVAL_GAP_MAX = float(max(_g))
_lab_min = EVAL[int(np.argmin(_g))][0]
_, _, THIRD_MIN, FOURTH_MIN = rank_at(EVAL[int(np.argmin(_g))][1], full=True)
print(f"\nSo BOTH things are true and they must be said together. The RESULT is robust: "
      f"{EVAL_HITS_MIN:.0f} of 3 at")
print(f"every one of these references. The GAP is not: it runs {EVAL_GAP_MIN:.4f} to "
      f"{EVAL_GAP_MAX:.4f} decades across them,")
print(f"a factor {EVAL_GAP_MAX/EVAL_GAP_MIN:.0f}, and it is smallest at {_lab_min}, where "
      f"{FOURTH_MIN} - the pair the")
print(f"paper calls '{CLASS[FOURTH_MIN]}' - comes within {100*(10**EVAL_GAP_MIN-1):.0f} % of overtaking "
      f"{THIRD_MIN}. Quoting '{gap_A:.3f} decades' as")
print(f"the margin without its reference is the same mistake this page's Result 1 is about.")
print(f"Outside that shared range it does eventually fail, which is a statement about")
_lines = ["extrapolating past the data rather than about the result:"]
for c in (1.0, 10.0, 1e5):
    h, gp = rank_at(lambda g, c=c: c)
    n_out = sum(not (PTS[PTS.gas_pair == g].P_fast_barrer.min() <= c
                     <= PTS[PTS.gas_pair == g].P_fast_barrer.max()) for g in PAIRS)
    _lines.append(f"  common P = {c:>8g} barrer lies outside {n_out} of the 9 pairs' own data -> "
                  f"{h} of 3")
print("\n".join(_lines))

# %%
# The picture: the nine comparable pairs, both bounds, and the paper's own near-bound points.
# Re-plotted from the transcribed tables - no source figure is traced or reproduced.
order = list(RK.index)
fig, axes = plt.subplots(3, 3, figsize=(12.6, 10.2))
for ax, g in zip(axes.ravel(), order):
    sub = PTS[PTS.gas_pair == g]
    k0, n0 = prior.k_barrer[g], prior.n[g]
    k1, n1 = present.k_barrer[g], present.n[g]
    Pgrid = np.logspace(np.log10(sub.P_fast_barrer.min()) - 0.4,
                        np.log10(sub.P_fast_barrer.max()) + 0.4, 200)
    ax.loglog(Pgrid, alpha_bound(Pgrid, k0, n0), "--", color="0.45", lw=1.4, label="prior (1991/94)")
    ax.loglog(Pgrid, alpha_bound(Pgrid, k1, n1), "-", color="C3", lw=1.6, label="present (2008)")
    r = excursion(sub.P_fast_barrer, sub.alpha, k1, n1)
    ax.loglog(sub.P_fast_barrer[r <= 1], sub.alpha[r <= 1], "o", ms=5, mfc="C0", mec="k",
              mew=0.4, label="tabulated, on/below")
    ax.loglog(sub.P_fast_barrer[r > 1], sub.alpha[r > 1], "^", ms=7, mfc="gold", mec="k",
              mew=0.6, label="tabulated, ABOVE")
    ax.set_title(f"{g}  ({CLASS[g]}, $\\Delta\\log\\alpha$ = {RK.dlogA[g]:.2f})", fontsize=9)
    ax.set_xlabel("$P$ fast gas, barrer", fontsize=8)
    ax.set_ylabel(r"$\alpha$", fontsize=8)
    ax.tick_params(labelsize=7)
axes.ravel()[0].legend(fontsize=7, loc="lower left")
fig.suptitle("Robeson upper bounds, nine comparable gas pairs, ordered by selectivity gain",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97)); plt.show()

# %% [markdown]
# ### Result 4: How many of the paper's own "close to the bound" points are above it
#
# Table 1–11 rows are, by their captions, *"experimental data points close to
# the present empirical upper bound"*. The line was drawn by eye. So: at each
# tabulated point, $r = P/(k\alpha^{n})$ against the pair's own present bound,
# and the same against the prior bound.

# %%
ex_rows = []
for g in PAIRS_ALL:
    if g not in set(PTS.gas_pair):
        continue
    sub = PTS[PTS.gas_pair == g]
    k1, n1 = present.k_barrer[g], present.n[g]
    r1 = np.asarray(excursion(sub.P_fast_barrer, sub.alpha, k1, n1), float)
    a_ratio = r1 ** (-1.0 / n1)              # the same excursion read as a selectivity ratio
    if np.isfinite(prior.k_barrer[g]):
        r0 = np.asarray(excursion(sub.P_fast_barrer, sub.alpha, prior.k_barrer[g], prior.n[g]), float)
    else:
        r0 = np.full(len(sub), np.nan)
    for (_, row), rr1, rr0, aa in zip(sub.iterrows(), r1, r0, a_ratio):
        ex_rows.append(dict(gas_pair=g, polymer=row.polymer, P=row.P_fast_barrer, alpha=row.alpha,
                            r_present=rr1, alpha_ratio=aa, r_prior=rr0))
EX = pd.DataFrame(ex_rows)
N_POINTS = len(EX)
N_ABOVE = int((EX.r_present > 1).sum())
FRAC_ABOVE = N_ABOVE / N_POINTS
R_MAX = float(EX.r_present.max()); R_MIN = float(EX.r_present.min())
R_GEOMEAN = float(np.exp(np.mean(np.log(EX.r_present))))

print("\n".join([f"{'pair':9s} {'N':>3s} {'above':>5s} {'r min':>7s} {'r max':>7s} {'geo-mean r':>10s}"]
                + [f"{g:9s} {len(sub):3d} {int((sub.r_present>1).sum()):5d} {sub.r_present.min():7.3f} "
                   f"{sub.r_present.max():7.3f} {np.exp(np.mean(np.log(sub.r_present))):10.3f}"
                   for g, sub in EX.groupby("gas_pair", sort=False)]))
print(f"\n{N_ABOVE} of {N_POINTS} tabulated points ({100*FRAC_ABOVE:.1f} %) lie ABOVE the very "
      f"bound they are tabulated against.")
print(f"r ranges {R_MIN:.3f} to {R_MAX:.3f}; geometric mean {R_GEOMEAN:.3f} - so the typical "
      f"tabulated point sits a factor")
print(f"{1/R_GEOMEAN:.1f} BELOW the line in permeability, and the line is a genuine envelope in the "
      f"aggregate.")

print("\n".join(["\nThe five largest excursions, and what the paper says about them:"]
                + [f"  {r.gas_pair:9s} P={r.P:8.4g} alpha={r.alpha:8.4g}  r={r.r_present:6.3f} "
                   f"(= {r.alpha_ratio:5.3f}x the bound selectivity)  {r.polymer[:44]}"
                   for _, r in EX.sort_values("r_present", ascending=False).head(5).iterrows()]))

# The two points the paper itself flags, located from the numbers in its own sentences
# (printed-claims dataset), not from numbers typed here.
import re
def quote_numbers(claim_id):
    """The numbers printed in one quoted sentence, in order."""
    q = CLM[CLM.claim_id == claim_id].quote.iloc[0]
    return [float(x) for x in re.findall(r"(?<![\w.])\d+(?:\.\d+)?", q)]
_qP, _qa = quote_numbers("o2n2_questioned")[-2:]           # "P(O2) = 18 barrers; alpha = 9.0"
_q, = EX.index[(EX.gas_pair == "O2/N2") & (EX.P == _qP) & (EX.alpha == _qa)]
Q = EX.loc[_q]
_rank = int((EX.r_present > Q.r_present).sum()) + 1
print(f"\n  The O2/N2 point at P = {Q.P:g}, alpha = {Q.alpha:g} is r = {Q.r_present:.3f}, excursion number "
      f"{_rank} of {len(EX)} in the whole compilation,")
print(f"  and the paper flags exactly it: 'The position of the one data point above the present upper")
print(f"  bound (P(O2) = 18 barrers; alpha(O2/N2) = 9.0) is questioned as only one significant")
print(f"  figure was noted for nitrogen permeability.'")
_o2 = EX[EX.gas_pair == "O2/N2"]
N_ABOVE_O2N2 = int((_o2.r_present > 1).sum())
_other = _o2[(_o2.r_present > 1) & (_o2.index != _q)].iloc[0]
print(f"  The arithmetic sharpens the paper's own wording: it says 'the ONE data point above")
print(f"  the present upper bound' for O2/N2, and Table 1 in fact has {N_ABOVE_O2N2}; the other is")
print(f"  {_other.polymer} at P = {_other.P:g}, alpha = {_other.alpha:g}, r = {_other.r_present:.3f}, "
      f"at the far low-permeability end")
print(f"  where a hand-drawn line is least constrained. Reported, not repaired.")

_listed, _realistic = quote_numbers("n2ch4_repair")         # "listed as 80.2 ... be 180.2"
_t9 = EX[EX.gas_pair == "N2/CH4"]
_m = _t9[np.round(_t9.P / _listed, 1) == _t9.alpha]     # the row whose alpha is P / 80.2
assert len(_m) == 1, "exactly one Table 9 row must carry the as-listed CH4 permeability"
N2 = _m.iloc[0]
_r9 = int((_t9.r_present > N2.r_present).sum()) + 1
print(f"\n  And the N2/CH4 point the paper doubts on other grounds (P = {N2.P:g}, alpha = {N2.alpha:g}, "
      f"the one whose CH4")
print(f"  permeability it says 'would be more realistic to be {_realistic:g} barrers') is excursion "
      f"number {_r9} of {len(_t9)} in Table 9:")
print(f"  r = {N2.r_present:.3f}. Consistency check on the tabulated alpha: {N2.P:g}/{_listed:g} = "
      f"{N2.P/_listed:.3f}, which rounds to the printed {N2.alpha:g};")
print(f"  {N2.P:g}/{_realistic:g} = {N2.P/_realistic:.3f} does not. So Table 9 carries the "
      f"UNREPAIRED value and this page uses it as printed.")
assert round(N2.P / _listed, 1) == N2.alpha != round(N2.P / _realistic, 1)

# %%
# The same measurement against the PRIOR bound: how far the 2008 points sit above the 1991 line.
# This is the shift measured with data instead of with fitted lines.
prior_rows = []
for g in PAIRS:
    sub = EX[EX.gas_pair == g]
    prior_rows.append(dict(pair=g, n=len(sub), n_above=int((sub.r_prior > 1).sum()),
                           rmax=sub.r_prior.max(),
                           gm=float(np.exp(np.mean(np.log(sub.r_prior))))))
PR = pd.DataFrame(prior_rows).set_index("pair")
print("\n".join([f"{'pair':9s} {'N':>3s} {'above prior':>11s} {'max r_prior':>12s} {'geo-mean r_prior':>17s}"]
                + [f"{g:9s} {int(r.n):3d} {int(r.n_above):11d} {r.rmax:12.3g} {r.gm:17.3f}"
                   for g, r in PR.iterrows()]))
N_ABOVE_PRIOR = int(PR.n_above.sum()); N_PRIOR_POINTS = int(PR.n.sum())
print(f"\n{N_ABOVE_PRIOR} of {N_PRIOR_POINTS} points on the nine comparable pairs exceed the 1991/94 bound "
      f"({100*N_ABOVE_PRIOR/N_PRIOR_POINTS:.0f} %),")
print(f"against {int((EX[EX.gas_pair.isin(PAIRS)].r_present>1).sum())} that exceed the 2008 one. That is the shift, seen from the data side.")
print(f"He/H2 is the extreme: all {PR.n['He/H2']} of its points beat the 1991 line, the furthest by a factor "
      f"{PR.rmax['He/H2']:.0f}.")
print(f"He/N2 is the quiet one: {PR.n_above['He/N2']} of {PR.n['He/N2']}, max {PR.rmax['He/N2']:.2f} - matching the paper's")
print(f"'The upper bound shift since 1991 is very minor in spite of the much larger dataset")
print(f"available presently.'")
print("\nCAUTION, and it is not small: these points are 2008 data. That they beat a 1991 line drawn")
print("on 1991 data is expected, not evidence of anything the paper claims. The number is a")
print("MAGNITUDE for the shift, computed without using k or n of the present bound at all - which")
print("is its only real use, and why it appears here rather than among the validations.")

# %% [markdown]
# ## Validation
#
# Four checks, and the first thing to say is which of them is a **second route**
# and which is a **consistency identity**: an identity dressed as corroboration
# is the error that costs a page its credibility.
#
# **V1: internal identity.** Tables 13a/13b are computed *from* Table 12, so
# reproducing all 26 of their rows is a genuine test of every present $(k,n)$
# pair and of the prior pair for He/H2 and H2/CO2. It shares Table 12 as an
# *input*, so it cannot detect an error Robeson himself made; it detects errors
# *this page* could have made, which is what it is for.
#
# **V2: a consistency identity, NOT a second route.** Tables 13a and 13b each
# evaluate He/H2 and H2/CO2 at one selectivity on **both** bounds, so the ratio
# of the two printed permeabilities looks like an independently computed shift.
# It is not independent. Robeson's transition permeability *is* $k\alpha^n$
# evaluated on the same Table 12 line V1 tests, so for two rows at the same
# $\alpha$
#
# $$
# \log_{10}\frac{P_1}{P_0}
# = \underbrace{(\Delta\log_{10}k + \Delta n\log_{10}\alpha)}_{\text{the decomposition}}
# - \underbrace{\big[\log_{10}(1+\rho_1) - \log_{10}(1+\rho_0)\big]}_{\text{the two V1 residuals}},
# $$
#
# with $\rho = P_{k\alpha^n}/P_{\text{printed}} - 1$ each row's V1 relative
# deviation, $P_{\text{printed}} = k\alpha^n/(1+\rho)$. So **V2's residual is the
# difference of the two corresponding V1 residuals**, and nothing else. The cell
# below measures that rather than asserting it. V2 therefore cannot fail unless
# V1 fails on those same four rows; what it *does* have is power over **this
# page's own decomposition code**. A natural-log-for-$\log_{10}$ slip changes the
# value of many claims, and V2 is the one that says the result is wrong rather
# than merely different: a residual that should sit at the level of the
# printing, below $10^{-3}$ decades, becomes of order one. That is the only
# reason V2 is kept, and it is a real one, but it is a unit test of this page,
# not evidence about Robeson.
#
# **Result 1's decomposition has no arithmetic second route, and cannot have
# one.** Given Table 12, $\Delta\log P(\alpha)=\Delta\log k+\Delta n\log\alpha$
# is a rearrangement, not a claim: any quantity built from the same two lines
# reproduces it identically. The cell below demonstrates that on a second
# candidate, the *data-side* shift (the ratio of the geometric-mean excursions
# of a pair's own points above the prior and present bounds), which also comes
# out equal to the decomposition to machine precision. Two identities are not
# two routes.
#
# **V3: the second route, and it is a real one.** The one check here that does
# *not* take Table 12's present block as an input: a covering line re-estimated
# from the paper's own 117 near-bound points by a deterministic envelope
# estimator. Substituting it for the printed present line gives an
# **envelope-implied shift** that shares no arithmetic and no input with the
# decomposition **except two things, both named**: the prior block, which both
# routes read from the same restatement, and the reference selectivity
# $\alpha_0$, which both routes take from Tables 1–11. The residual between the
# two routes is exactly V3's offset at $\alpha_0$, which is why that quantity is
# the headline of this section rather than a footnote, and, because it carries
# $\log\alpha_0$ explicitly, it is itself reference-dependent for the same reason
# the front-factor share is, so it is swept across references and reported with
# its reference attached, exactly as Result 3's gap is. It is weak where the
# lever arm is short and it is labelled so, with the measurement to show which
# pairs those are. It is also the only check on this page whose failure could
# indict *Robeson* rather than this page.
#
# **V4** is the break table: every claim on the page is recomputed with each
# defect injected, and a claim counts as protected only where its own check
# fails.

# %% [markdown]
# ### V1: Tables 13a/13b reproduce from Table 12
#
# Three nested checks, each catching something the previous cannot:
# the molecular-weight column from molecular weights; the product column from
# its two factors; and the permeability column from Table 12's $(k,n)$.

# %%
MW = {"He": 4.0, "H2": 2.0, "O2": 32.0, "N2": 28.0, "CO2": 44.0, "CH4": 16.0}
def sqrtM(pair):
    i, j = pair.split("/")
    return np.sqrt(MW[j] / MW[i])

def ulp(s):
    """Half a unit in the last printed decimal place of the string `s` (commas stripped)."""
    s = str(s).replace(",", "").lstrip("-")
    return 0.5 * (10.0 ** (-len(s.split(".")[1])) if "." in s else 1.0)

# --- (a) the (Mj/Mi)^(1/2) column, from integer molecular weights -----------------
T13 = T13.assign(sqrtM_calc=[sqrtM(g) for g in T13.gas_pair])
T13 = T13.assign(mw_rel=lambda d: d.sqrtM_calc / d.sqrt_M_ratio - 1)
MW_MAX_REL = float(np.max(np.abs(T13.mw_rel)))
IUPAC = {"He": 4.002602, "H2": 2.016, "O2": 31.998, "N2": 28.014, "CO2": 44.009, "CH4": 16.043}
sqrtM_iupac = np.array([np.sqrt(IUPAC[g.split("/")[1]] / IUPAC[g.split("/")[0]]) for g in T13.gas_pair])
MW_MAX_REL_IUPAC = float(np.max(np.abs(sqrtM_iupac / T13.sqrt_M_ratio - 1)))
print(f"(a) (Mj/Mi)^(1/2) column, {len(T13)} rows")
print(f"    integer molecular weights : max |rel. dev| = {MW_MAX_REL:.2e}")
print(f"    IUPAC standard weights    : max |rel. dev| = {MW_MAX_REL_IUPAC:.2e}  "
      f"({MW_MAX_REL_IUPAC/MW_MAX_REL:.0f}x worse)")
print(f"    -> the convention is identified, not assumed: H2/N2 prints "
      f"{T13[(T13.table == '13a') & (T13.gas_pair == 'H2/N2')].sqrt_M_ratio.iloc[0]} = sqrt(28/2), not "
      f"{np.sqrt(IUPAC['N2']/IUPAC['H2']):.3f}.")

# --- (b) the product column ks*(Mj/Mi)^(1/2) --------------------------------------
T13 = T13.assign(prod_calc=lambda d: d.ks * d.sqrt_M_ratio,
                 prod_rel=lambda d: d.ks * d.sqrt_M_ratio / d.ks_D_ratio - 1)
PROD_MAX_REL = float(np.max(np.abs(T13.prod_rel)))
worst_prod = T13.loc[T13.prod_rel.abs().idxmax()]
print(f"\n(b) product column ks*(Mj/Mi)^(1/2), {len(T13)} rows: max |rel. dev| = {PROD_MAX_REL:.2e}")
print(f"    worst row {worst_prod.table} {worst_prod.gas_pair}: printed {worst_prod.ks_D_ratio}, "
      f"computed {worst_prod.prod_calc:.6g}")
print(f"    (printed to {len(worst_prod.ks_D_as_printed.split('.')[-1])} decimals, so a rounding "
      f"band of +-{ulp(worst_prod.ks_D_as_printed)/worst_prod.ks_D_ratio*100:.3f} % on that cell alone)")

# %%
# --- (c) the permeability column, from Table 12's (k, n) --------------------------
# Worst-case band: each printed input contributes its own half-ulp, added in absolute value.
#   dP/P from k : ulp(k)/k
#   dP/P from n : |ln(alpha)| * ulp(n)
#   dP/P from alpha (the product column) : |n| * ulp(alpha)/alpha
#   plus the half-ulp of the printed P itself.
res = []
for _, r in T13.iterrows():
    src = prior if r.dataset == "prior" else present
    k, n = src.k_barrer[r.gas_pair], src.n[r.gas_pair]
    kp = src.k_as_printed[r.gas_pair]; npr = src.n_as_printed[r.gas_pair]
    P_calc = float(P_bound(r.ks_D_ratio, k, n))
    band_tight = (ulp(kp) / k + abs(np.log(r.ks_D_ratio)) * ulp(npr)
                  + abs(n) * ulp(r.ks_D_as_printed) / r.ks_D_ratio
                  + ulp(r.P_transition_barrer) / r.P_transition_barrer)
    # generous band: trailing zeros are placeholders -> 4 significant figures on k and on P
    sig4 = lambda v: 0.5 * 10.0 ** (np.floor(np.log10(abs(v))) - 3)     # noqa: E731
    band_loose = (sig4(k) / k + abs(np.log(r.ks_D_ratio)) * ulp(npr)
                  + abs(n) * ulp(r.ks_D_as_printed) / r.ks_D_ratio
                  + sig4(r.P_transition_barrer) / r.P_transition_barrer)
    res.append(dict(table=r.table, pair=r.gas_pair, dataset=r.dataset, P_printed=r.P_transition_barrer,
                    P_calc=P_calc, rel=P_calc / r.P_transition_barrer - 1,
                    band_tight=band_tight, band_loose=band_loose))
V1 = pd.DataFrame(res)
V1 = V1.assign(inside_tight=lambda d: d.rel.abs() <= d.band_tight,
               inside_loose=lambda d: d.rel.abs() <= d.band_loose)
T13_MAX_REL = float(V1.rel.abs().max())
T13_MED_REL = float(V1.rel.abs().median())
N_INSIDE_TIGHT = int(V1.inside_tight.sum()); N_INSIDE_LOOSE = int(V1.inside_loose.sum())
worst = V1.loc[V1.rel.abs().idxmax()]
# Computed unconditionally, not inside the diagnostic branch: alpha_bound applied to the
# printed permeability of the row that fails, so it exists whether or not the branch runs.
_r13a = T13[(T13.table == "13a") & (T13.gas_pair == "O2/N2")].iloc[0]
T13A_O2N2_BACKOUT_REL = float(alpha_bound(_r13a.P_transition_barrer, present.k_barrer["O2/N2"],
                                          present.n["O2/N2"]) / _r13a.ks_D_ratio - 1)

print("\n".join([f"{'tbl':4s} {'pair':9s} {'set':8s} {'P printed':>12s} {'P = k*alpha^n':>14s} "
                 f"{'rel dev':>9s} {'tight band':>11s} {'in?':>4s}"]
                + [f"{r.table:4s} {r.pair:9s} {r.dataset:8s} {r.P_printed:12,.0f} {r.P_calc:14,.0f} "
                   f"{100*r.rel:+8.3f}% {100*r.band_tight:10.3f}% {'yes' if r.inside_tight else 'NO':>4s}"
                   for _, r in V1.iterrows()]))
print(f"\n{len(V1)} rows. Median |deviation| {100*T13_MED_REL:.3f} %, max {100*T13_MAX_REL:.3f} % "
      f"({worst.table} {worst.pair}).")
print(f"Inside the LITERAL printing band: {N_INSIDE_TIGHT}/{len(V1)}.  Inside the GENEROUS "
      f"(4-significant-figure) band: {N_INSIDE_LOOSE}/{len(V1)}.")
out = V1[~V1.inside_loose]
if len(out):
    for _, r in out.iterrows():
        print(f"\nOUTSIDE EVEN THE GENEROUS BAND: {r.table} {r.pair}, {100*r.rel:+.3f} % against "
              f"+-{100*r.band_loose:.3f} %.")
        src = present if r.dataset == "present" else prior
        k, n = src.k_barrer[r.pair], src.n[r.pair]
        a_needed = float(alpha_bound(r.P_printed, k, n))
        a_print = float(T13[(T13.table == r.table) & (T13.gas_pair == r.pair)].ks_D_ratio.iloc[0])
        print(f"  Backing the printed P out through the same bound needs alpha = {a_needed:.4f}; "
              f"the table prints {a_print}.")
        print(f"  That is a {100*(a_needed/a_print-1):+.2f} % difference in the last digit of a "
              f"{len(_r13a.ks_D_as_printed.replace('.', ''))}-digit cell - a small internal")
        print(f"  inconsistency in one row of Table 13a, reported and not repaired. It changes no")
        print(f"  conclusion here: {r.pair} is a Knudsen-transition estimate, and Table 12's own "
              f"(k, n) for it are")
        print(f"  confirmed by the other table ({100*float(V1[(V1.pair==r.pair)&(V1.table=='13b')].rel.iloc[0]):+.3f} % in 13b).")
_out_t = V1[~V1.inside_tight]
print("\n".join(f"Outside the LITERAL band: {r.table} {r.pair}, |deviation| {100*abs(r.rel):.3f} % against "
                f"{100*r.band_tight:.3f} %, {100*(abs(r.rel)/r.band_tight - 1):.1f} % over it."
                for _, r in _out_t.iterrows()))
print(f"\nWhat V1 can detect: a misread glyph, a lost digit or a dropped comma in any present (k, n)")
print(f"of Table 12, or in the prior (k, n) of He/H2 and H2/CO2, which moves a row far outside its")
print(f"printing band (break rows C, G9, G8, DD, DP). What it cannot detect on its own: an error in")
print(f"the LAST printed digit, which stays inside the band (break row D1 leaves the two present")
print(f"He/H2 rows inside theirs); such an error is caught only against the committed baseline. Nor anything")
print(f"about Table 12 itself, which is an input on both sides, nor the prior (k, n) of the other")
print(f"seven pairs, which Table 13 never touches.")

# %% [markdown]
# ### V2: a consistency identity on the decomposition code
#
# For He/H2 and H2/CO2, Tables 13a and 13b each print the transition
# permeability **twice**, once on the present bound and once on the prior one,
# at *identical* selectivity. Their ratio is $\Delta\log_{10}P(\alpha)$ at that
# selectivity, set beside $\Delta\log_{10}k + \Delta n\log_{10}\alpha$.
#
# **What this can and cannot catch.** It can catch a natural-log-for-$\log_{10}$
# slip, a sign error, or a prior/present swap in *this page's* decomposition
# formula; the break table below shows it failing under all three. It
# cannot catch anything about Table 12, and it is not a second route to Result 1:
# the printed side is $k\alpha^n$ on the same two lines, so the comparison is an
# algebraic tautology whose residual is fixed entirely by V1. The cell measures
# that: the four V2 residuals are set against the corresponding differences of V1
# residuals, and against a second candidate identity, the shift read off the
# *data* rather than off the lines.
#
# What the four rows are genuinely good for is the **reference-dependence**: the
# same two pairs, the same two bounds, four different selectivities, and a
# front-factor share that moves with them. Those are numbers
# Robeson printed, not numbers this page chose.

# %%
v2 = []
for tab in ("13a", "13b"):
    sub = T13[T13.table == tab]
    for g in sub.gas_pair[sub.dataset == "prior"]:
        row_p = sub[(sub.gas_pair == g) & (sub.dataset == "present")].iloc[0]
        row_q = sub[(sub.gas_pair == g) & (sub.dataset == "prior")].iloc[0]
        assert row_p.ks_D_ratio == row_q.ks_D_ratio, "the two rows must sit at the same alpha"
        alpha = float(row_p.ks_D_ratio)
        printed = np.log10(row_p.P_transition_barrer / row_q.P_transition_barrer)
        decomp = SH.dlogk[g] + SH.dn[g] * np.log10(alpha)
        v2.append(dict(table=tab, pair=g, alpha=alpha, printed=printed, decomposed=decomp,
                       diff=decomp - printed, dlogk=SH.dlogk[g],
                       k_share=SH.dlogk[g] / decomp))
V2 = pd.DataFrame(v2)
V2_MAX_ABS = float(V2["diff"].abs().max())
print("\n".join([f"{'tbl':4s} {'pair':8s} {'alpha':>9s} {'shift, printed':>15s} {'shift, decomposed':>18s} "
                 f"{'diff, decades':>14s} {'dlog k':>8s} {'k-share':>8s}"]
                + [f"{r.table:4s} {r.pair:8s} {r.alpha:9.5f} {r.printed:15.5f} {r.decomposed:18.5f} "
                   f"{r['diff']:+14.2e} {r.dlogk:8.4f} {r.k_share:8.3f}" for _, r in V2.iterrows()]))
print(f"\nWorst |difference| over the four comparisons: {V2_MAX_ABS:.2e} decades.")

# --- the identity, measured rather than asserted --------------------------------------
# V2's residual is log10((1+rho_present)/(1+rho_prior)) where rho is each row's V1 relative
# deviation. Not "correlated with" V1: EQUAL to it, up to floating-point re-association.
print("\nIs this a second route? No - and here is the arithmetic that settles it. Robeson's")
print("transition permeability IS k*alpha^n on the same two lines V1 tests, so the ratio of the")
print("two printed permeabilities carries the decomposition plus the two V1 residuals and nothing")
print("else. Set V2's residual against the difference of the corresponding V1 residuals:")
_v2id, _lines = [], [f"\n{'tbl':4s} {'pair':8s} {'V2 residual':>14s} {'V1(pres) - V1(prior)':>22s}"]
for _, r in V2.iterrows():
    _p = V1[(V1.table == r.table) & (V1.pair == r.pair) & (V1.dataset == "present")].rel.iloc[0]
    _q = V1[(V1.table == r.table) & (V1.pair == r.pair) & (V1.dataset == "prior")].rel.iloc[0]
    _d = float(np.log10((1 + _p) / (1 + _q)))
    _v2id.append(abs(r["diff"] - _d))
    _lines.append(f"{r.table:4s} {r.pair:8s} {r['diff']:+14.6e} {_d:+22.6e}")
print("\n".join(_lines))
V2_IS_V1_MAX_ABS = float(max(_v2id))
print(f"\nThey agree below 1e-12 decades (here {V2_IS_V1_MAX_ABS:.0e}) - floating-point re-association, not")
print(f"agreement between two computations. V2 cannot fail unless V1 fails on these four rows;")
print(f"the headline {V2_MAX_ABS:.1e} decades is the V1 table above, differenced.")

# --- and the other candidate identity, for the same reason ----------------------------
# The "data-side" shift: how much further a pair's own points sit above the prior bound than
# above the present one, in the geometric mean. It uses no (k, n) ratio explicitly - and it is
# STILL the decomposition, because averaging log r over the points returns log alpha_gm.
print("\nThe obvious second candidate fails the same way. Take the shift straight off the DATA:")
print("the geometric-mean excursion of a pair's own points above the prior bound, over the same")
print("quantity for the present bound. No Delta log k, no Delta n anywhere in it:")
_ds, _lines = [], [f"\n{'pair':9s} {'data-side shift':>16s} {'decomposition @ a_gm':>22s}"]
for g in PAIRS:
    _sub = EX[EX.gas_pair == g]
    _d = float(np.log10(np.exp(np.mean(np.log(_sub.r_prior))) / np.exp(np.mean(np.log(_sub.r_present)))))
    _ds.append(abs(_d - SH.shift_gm[g]))
    _lines.append(f"{g:9s} {_d:16.10f} {SH.shift_gm[g]:22.10f}")
print("\n".join(_lines))
DATASIDE_IS_DECOMP_MAX_ABS = float(max(_ds))
print(f"\nEqual below 1e-12 decades (here {DATASIDE_IS_DECOMP_MAX_ABS:.0e}), on all nine pairs. The reason is one line:")
print("averaging log r over a pair's points subtracts n*mean(log alpha), and mean(log alpha) IS")
print("log(alpha_gm), so the data-side shift collapses to dlog k + dn*log(alpha_gm) exactly.")
print("\nSo: given Table 12, Result 1's decomposition is a REARRANGEMENT, and no arithmetic on")
print("Table 12 can corroborate it. The route that can is V3 - which replaces the present line")
print("with one estimated from the 117 points instead of read from Table 12.")

print(f"\nAnd read the last two columns of the table above. The SAME two pairs, the SAME two bounds, four different")
print(f"selectivities spanning {V2.alpha.max()/V2.alpha.min():.0f}x - and the front-factor share of the shift runs")
print(f"{V2.k_share.min():.2f} to {V2.k_share.max():.2f}. For He/H2 at alpha = {V2[V2.pair=='He/H2'].alpha.iloc[0]:.3f}, "
      f"dlog k OVERSTATES the shift by {100*(V2[V2.pair=='He/H2'].k_share.iloc[0]-1):.1f} %.")
print(f"For H2/CO2 at alpha = {V2[V2.pair=='H2/CO2'].alpha.iloc[0]:.5f}, dlog k accounts for only "
      f"{100*V2[V2.pair=='H2/CO2'].k_share.iloc[0]:.0f} % of it - and at alpha = {ALPHA_STAR:.1f}")
print(f"the same two lines cross and the shift is exactly zero, so the share is undefined there.")
print(f"THIS is what 'the split depends on where you stand' means, demonstrated on numbers the")
print(f"paper printed rather than on numbers this page chose.")


# %% [markdown]
# ### V3: an envelope re-estimated from the paper's own points, and the second route to Result 1
#
# The only check here that touches the **empirical** content rather than the
# arithmetic, and the only one that does not take Table 12's present block as an
# input. The estimator is stated so it can be argued with: the line
# $\log P = \log k + n\log\alpha$ minimising
# $\sum_i (\log P_i - \log k - n\log\alpha_i)^2$ **subject to every point lying
# on or below it**. That is a convex quadratic program; its optimum activates
# at most two constraints, so it is solved here by **exact enumeration of the
# active sets** rather than by an iterative optimiser: no start point, no seed,
# one feasibility tolerance at round-off level (`envelope_fit` in `model.py`).
#
# It is *not* what "by eye" means. A covering line is pinned by two or three
# extreme points; an eye fits the body of a cloud. The two estimators answer
# different questions, and this page keeps them apart, **nothing below is
# offered as a correction to the printed $(k,n)$.**
#
# **Why this is nevertheless Result 1's second route.** Substitute
# $(k_{\text{env}}, n_{\text{env}})$ for the printed present line and recompute
# the shift at the same reference:
#
# $$
# \Delta\log_{10}P_{\text{env}}(\alpha_0)
# = \log_{10}\frac{k_{\text{env}}}{k_0} + (n_{\text{env}} - n_0)\log_{10}\alpha_0 .
# $$
#
# Its prior half is shared with the decomposition; its present half comes from
# 117 points and a quadratic program instead of from Table 12. Subtracting, the
# prior line cancels exactly, so the **residual between the two routes is
# V3's offset at $\alpha_0$**, the column already tabulated below. That is why
# the offset is the headline of this section: it is not a curiosity about a
# re-fit, it is the disagreement between the only two ways this page can compute
# how far the bound moved. Its scope is exactly the present line; it says
# nothing about the prior block, which both routes read from the same place, and
# nothing about $\alpha_0$, which both routes take from the same 117 points.
#
# **And it carries $\log\alpha_0$, so it is reference-dependent too.** That is
# not a caveat bolted on: it is this page's own thesis pointed at this page's own
# newest headline. Quoting "the two routes agree to *x* decades" without naming
# $\alpha_0$ is the same error as quoting a front-factor share without naming
# $\alpha_0$. So the reference is swept below, over each pair's own minimum,
# geometric-mean, median and maximum tabulated selectivity, all four inside the
# pair's own data, and the residual is reported with its reference attached
# wherever it appears, on this page and on every surface that quotes it.

# %%
v3 = []
for g in sorted(PTS.gas_pair.unique()):
    sub = PTS[PTS.gas_pair == g]
    la = np.log10(sub.alpha.to_numpy(float))
    lp = np.log10(sub.P_fast_barrer.to_numpy(float))
    k_env, n_env, sse = envelope_fit(la, lp)
    k1, n1 = present.k_barrer[g], present.n[g]
    a_gm = float(np.exp(np.mean(np.log(sub.alpha.to_numpy(float)))))
    off = ((np.log10(k_env) + n_env * np.log10(a_gm))
           - (np.log10(k1) + n1 * np.log10(a_gm)))
    v3.append(dict(pair=g, n_pr=n1, n_env=n_env, dn=n_env - n1, k_pr=k1, k_env=k_env,
                   span=la.max() - la.min(), offset_gm=off, sse=sse))
V3 = pd.DataFrame(v3).set_index("pair").sort_values("span")
V3_DN_MED = float(V3.dn.abs().median()); V3_DN_MAX = float(V3.dn.abs().max())
V3_OFFSET_MED = float(V3.offset_gm.median()); V3_OFFSET_MAX = float(V3.offset_gm.abs().max())
V3_DN_WIDE = float(V3[V3.span > 1.5].dn.abs().max()); V3_SPAN_MIN = float(V3.span.min())
print("\n".join(["Sorted by how many decades of selectivity the pair's own tabulated points span.",
                 f"{'pair':9s} {'alpha span, dec':>15s} {'n printed':>10s} {'n envelope':>11s} {'dn':>8s} "
                 f"{'offset @ a_gm, dec':>19s}"]
                + [f"{g:9s} {r.span:15.3f} {r.n_pr:10.4f} {r.n_env:11.4f} {r.dn:+8.4f} {r.offset_gm:+19.4f}"
                   for g, r in V3.iterrows()]))
N_OFFSET_POS = int((V3.offset_gm > 0).sum())
print(f"\nMedian |dn| = {V3_DN_MED:.3f}, worst {V3_DN_MAX:.3f}. Vertical offset at each pair's own")
print(f"geometric-mean selectivity: median {V3_OFFSET_MED:+.3f} decades, largest magnitude "
      f"{V3_OFFSET_MAX:.3f}.")
print(f"\nTwo things, and only two, come out of this.")
print(f"\n(i) The offsets are SMALL and of BOTH SIGNS ({N_OFFSET_POS} of {len(V3)} positive). Over the "
      f"body of each cloud the printed line")
print(f"and the tightest covering line sit within a factor {10**V3_OFFSET_MAX:.1f} of each other, usually far")
print(f"closer - so the {N_ABOVE} excursions of Result 4 are individual points poking through a line that")
print(f"is in the right PLACE, not a line drawn in the wrong place.")
# The sign count of the offsets depends on where it is read; the slope difference does not.
_cnt_at = {lab: int(sum((np.log10(V3.k_env[g]) + V3.n_env[g] * np.log10(af(g)))
                        > (np.log10(present.k_barrer[g]) + present.n[g] * np.log10(af(g)))
                        for g in V3.index))
           for lab, af in [("alpha_min", lambda g: PTS[PTS.gas_pair == g].alpha.min()),
                           ("alpha_gm", lambda g: np.exp(np.mean(np.log(PTS[PTS.gas_pair == g].alpha)))),
                           ("alpha_median", lambda g: PTS[PTS.gas_pair == g].alpha.median()),
                           ("alpha_max", lambda g: PTS[PTS.gas_pair == g].alpha.max())]}
N_DN_POS = int((V3.dn > 0).sum())
print(f"The count of positive offsets is itself a property of the reference: "
      + ", ".join(f"{v} of {len(V3)} at each pair's {k}" for k, v in _cnt_at.items()) + ".")
print(f"What holds at every reference is the slope: dn = n_env - n_printed > 0 for {N_DN_POS} of "
      f"{len(V3)} pairs (smallest |dn| {V3.dn.abs().min():.3f}),")
print(f"so the covering line is shallower than the printed one, not uniformly higher or lower.")
# Stated as a SET match, computed, because the two orderings do not agree rank by rank.
# The largest m for which they agree AS SETS is the statement the data support.
_by_span = list(V3.span.sort_values().index)
_by_dn = list(V3.dn.abs().sort_values(ascending=False).index)
M_LEVER = max(m for m in range(1, len(V3)) if set(_by_span[:m]) == set(_by_dn[:m]))
_off = next(g for i, g in enumerate(_by_span[:M_LEVER]) if _by_dn[i] != g)
print(f"\n(ii) The slope disagreement follows the LEVER ARM, and it lands exactly")
print(f"where the paper itself says it should. The {M_LEVER} pairs whose points span the fewest decades of")
print(f"selectivity are exactly the {M_LEVER} carrying the largest |dn| - as SETS, not rank by rank "
      f"({_off} has the")
print(f"span rank {_by_span.index(_off) + 1} and the |dn| rank {_by_dn.index(_off) + 1}). Those {M_LEVER} are "
      f"{', '.join(f'{g} ({V3.span[g]:.2f} dec, |dn| {abs(V3.dn[g]):.3f})' for g in _by_span[:M_LEVER])}.")
print(f"The {int((V3.span > 1.5).sum())} pairs spanning more than 1.5 decades")
print(f"all agree to |dn| <= {V3[V3.span > 1.5].dn.abs().max():.3f}. For He/H2 - the paper's headline pair, and the one")
print(f"whose {int((PTS.gas_pair == 'He/H2').sum())} points cover only {V3.span['He/H2']:.2f} decades of "
      f"alpha - the tabulated points simply do")
print(f"not determine a slope, and Robeson says so: 'While the front factor has moved")
print(f"significantly, the slope does not appear to have changed although sufficient data does not")
print(f"exist to clearly confirm that observation.' V3 is the arithmetic behind that sentence.")
print(f"\nWhat V3 is NOT: evidence that any printed (k, n) is wrong. It cannot be - it optimises a")
print(f"different objective on a subset of the data Robeson drew the line through (Tables 1-11 are")
print(f"the near-bound points, not the ~300 papers' worth of cloud), and for O2/N2 it is pinned by")
print(f"the very point the paper says it distrusts.")

# %%
# ---- the second route to Result 1, run against the decomposition -----------------------
# Present line from the 117 points (V3) + prior line from Table 12, versus the decomposition.
# The prior half cancels in the difference, so the residual IS V3's offset at alpha_gm - which
# is asserted below rather than taken on trust.
sr = []
for g in PAIRS:
    k0, n0 = prior.k_barrer[g], prior.n[g]
    a_gm = SH.a_gm[g]
    k_env, n_env = V3.k_env[g], V3.n_env[g]
    sh_env = (np.log10(k_env) - np.log10(k0)) + (n_env - n0) * np.log10(a_gm)
    sr.append(dict(pair=g, span=V3.span[g], sh_dec=SH.shift_gm[g], sh_env=sh_env,
                   resid=sh_env - SH.shift_gm[g], share_dec=SH.share[g],
                   share_env=(np.log10(k_env) - np.log10(k0)) / sh_env))
SR = pd.DataFrame(sr).set_index("pair").sort_values("span")
assert np.allclose(SR.resid, V3.offset_gm[SR.index], rtol=0, atol=1e-12), \
    "the two-route residual must BE V3's offset at alpha_gm - the prior line cancels"
SR_WIDE = SR[SR.span > 1.5]
SECOND_ROUTE_MAX = float(SR.resid.abs().max())
SECOND_ROUTE_MED = float(SR.resid.abs().median())
SECOND_ROUTE_WIDE_MAX = float(SR_WIDE.resid.abs().max())
SECOND_ROUTE_SIGN_AGREE = int((np.sign(SR.sh_env) == np.sign(SR.sh_dec)).sum())
SHARE_ENV_MED = float(SR.share_env.median())

print("\n".join(
    ["The shift, computed twice: from Table 12 (the decomposition) and from Tables 1-11 (the",
     "envelope) - sorted by how many decades of selectivity the pair's own points span, because",
     "that is the lever arm the second route is built on.",
     "EVERY number in this table is evaluated at each pair's own geometric-mean selectivity a_gm,",
     "and none of them means anything without that - see the reference sweep below.",
     f"{'pair':9s} {'span, dec':>10s} {'shift, decomp':>14s} {'shift, envelope':>16s} "
     f"{'resid @ a_gm':>12s} {'k-share, decomp':>16s} {'k-share, env':>13s}"]
    + [f"{g:9s} {r.span:10.2f} {r.sh_dec:14.4f} {r.sh_env:16.4f} {r.resid:+12.4f} "
       f"{r.share_dec:16.3f} {r.share_env:13.3f}" for g, r in SR.iterrows()]))
print(f"\nAT a_gm. Sign of the shift: the two routes agree on {SECOND_ROUTE_SIGN_AGREE} of {len(SR)} "
      f"pairs - every bound that moved up by")
print(f"one route moved up by the other. MAGNITUDE: median |residual| {SECOND_ROUTE_MED:.4f} decades, "
      f"worst {SECOND_ROUTE_MAX:.4f}")
print(f"({SR.resid.abs().idxmax()}). Restricted to the {len(SR_WIDE)} pairs whose points span more than "
      f"1.5 decades of selectivity,")
print(f"where the covering line is actually determined: worst |residual| {SECOND_ROUTE_WIDE_MAX:.4f} "
      f"decades, i.e. the")
print(f"two routes agree to a factor {10**SECOND_ROUTE_WIDE_MAX:.2f} in permeability. The two pairs "
      f"that blow that out are")
print(f"{', '.join(SR.index[:2])}, spanning {SR.span.iloc[0]:.2f} and {SR.span.iloc[1]:.2f} decades - "
      f"the same short-lever-arm pairs")
print(f"that carry the largest |dn| above. This is a weak check where it is weak, and it says so.")
print(f"\nAnd now the sharp part, which is the page's central claim arriving from the other side.")
print(f"The two routes agree about the SIZE of the shift and disagree completely about its SPLIT:")
print(f"the front-factor share has median {SH.share.median():.2f} by the decomposition and "
      f"{SHARE_ENV_MED:.2f} by the envelope,")
print(f"ranging {SR.share_env.min():.2f} to {SR.share_env.max():.2f} instead of "
      f"{SH.share.min():.2f} to {SH.share.max():.2f}. That is not a defect in either route.")
print(f"The share is dlog k / dlog P, and dlog k is the shift at alpha = 1, up to "
      f"{np.log10(SH.a_gm).max():.2f} decades from a pair's")
print(f"geometric-mean selectivity. Where alpha = 1 lies outside a pair's data or among its last few")
print(f"points, it is read along a slope the points barely determine, and the two routes' shares part.")
print(f"Pair by pair, ordered by the fraction of its points below alpha = 1:")
_nb = {g: int((PTS[PTS.gas_pair == g].alpha < 1).sum()) for g in SR.index}
print("\n".join(f"  {g:9s} {_nb[g]:2d} of {int((PTS.gas_pair == g).sum()):2d} points below alpha = 1: "
                f"share {SR.share_dec[g]:.2f} (decomposition) against {SR.share_env[g]:.2f} (envelope)"
                for g in sorted(SR.index, key=lambda g: -_nb[g] / (PTS.gas_pair == g).sum())))
print(f"H2/CO2 and He/CO2 agree to {max(abs(SR.share_dec[g] - SR.share_env[g]) for g in ('H2/CO2', 'He/CO2')):.2f}; "
      f"H2/CH4, with as large a fraction of its points below alpha = 1 as He/CO2,")
print(f"does not ({SR.share_dec['H2/CH4']:.2f} against {SR.share_env['H2/CH4']:.2f}), so the position of "
      f"alpha = 1 is not the whole story.")
print(f"A quantity that survives one estimator and not the other is not a property of the data.")

# ---- what a NEGATIVE front-factor share is, since the range contains one ---------------
_gneg = SR.share_env.idxmin()
_k0neg, _kEneg = prior.k_barrer[_gneg], V3.k_env[_gneg]
_dlkneg = np.log10(_kEneg) - np.log10(_k0neg)
_across = 10.0 ** (-_dlkneg / (V3.n_env[_gneg] - prior.n[_gneg]))
print(f"\nOne of those envelope shares is NEGATIVE - {SR.share_env.min():.4f}, on {_gneg} - and a "
      f"negative share needs saying")
print(f"in words rather than being left inside a range. The envelope's front factor for {_gneg} is "
      f"{_kEneg:,.0f} barrers")
print(f"against the prior block's {_k0neg:,.0f} (and the printed present {present.k_barrer[_gneg]:,.0f}): "
      f"dlog k = {_dlkneg:+.4f}, i.e. the front")
print(f"factor moved DOWN, while the envelope-implied shift at a_gm is {SR.sh_env[_gneg]:+.4f} decades. "
      f"So the ENTIRE")
print(f"shift, and a little more, is carried by the slope term. Geometrically: the envelope present")
print(f"line sits BELOW the prior line at alpha = 1 and above it at a_gm = {SH.a_gm[_gneg]:.1f}; the two "
      f"cross at alpha = {_across:.2f},")
print(f"and every decade of separation the bound gained at a_gm was bought by rotation, not by lift.")
_nb1 = int((PTS[PTS.gas_pair == _gneg].alpha < 1).sum())
print(f"A share outside [0, 1] is not a defect - it is what happens when dlog k is read at a")
print(f"selectivity {np.log10(SH.a_gm[_gneg]):.1f} decades below the pair's geometric mean, where "
      f"{_nb1} of its {int((PTS.gas_pair == _gneg).sum())} points lie, which is exactly Result 1's warning.")

# ---- the residual's own reference is a CHOICE too, so sweep it -------------------------
# resid(a0) = (log k_env - log k1) + (n_env - n1) log a0 carries log a0 EXPLICITLY, so it moves
# with the reference for precisely the reason the front-factor share does: the two present lines
# are not parallel. Swept over four references, every one of them inside each pair's own data,
# exactly as EVAL sweeps P for Result 3 two sections up.
A0 = [("each pair's own alpha_min", lambda g: float(SH.a_lo[g])),
      ("each pair's own alpha_gm (used above)", lambda g: float(SH.a_gm[g])),
      ("each pair's own alpha_median",
       lambda g: float(np.median(PTS[PTS.gas_pair == g].alpha.to_numpy(float)))),
      ("each pair's own alpha_max", lambda g: float(SH.a_hi[g]))]
def route_at(af):
    """(max |resid| over the wide-span pairs, median |resid| over all nine, sign agreement)."""
    wide, allr, sgn, worst = [], [], 0, None
    for g in PAIRS:
        a0 = af(g)
        se = ((np.log10(V3.k_env[g]) - np.log10(prior.k_barrer[g]))
              + (V3.n_env[g] - prior.n[g]) * np.log10(a0))
        sd = SH.dlogk[g] + SH.dn[g] * np.log10(a0)
        allr.append(abs(se - sd))
        if np.sign(se) == np.sign(sd):
            sgn += 1
        elif worst is None:
            worst = (g, se, sd)
        if V3.span[g] > 1.5:
            wide.append(abs(se - sd))
    return max(wide), float(np.median(allr)), sgn, worst
print(f"\nAnd now the same discipline applied to THIS section's own headline. The residual is")
print(f"(log k_env - log k1) + (n_env - n1) log a0: it carries log a0 explicitly, so '{SECOND_ROUTE_WIDE_MAX:.3f} decades'")
print(f"is a property of the two present lines AT ONE SELECTIVITY, not of the two lines. Swept:")
_lines = [f"\n{'reference selectivity':40s} {'max|resid|, wide':>17s} {'median|resid|, all 9':>21s} "
          f"{'sign agreement':>15s}"]
_rw, _rm, _rs = [], [], []
for lab, af in A0:
    w, m, sg, worst = route_at(af)
    _rw.append(w); _rm.append(m); _rs.append(sg)
    _lines.append(f"{lab:40s} {w:17.4f} {m:21.4f} {sg:12d} of {len(PAIRS)}")
    if worst is not None:
        _lines.append(f"{'':40s}   ^ sign disagreement: {worst[0]} envelope {worst[1]:+.4f} vs "
                      f"decomposition {worst[2]:+.4f}")
print("\n".join(_lines))
SR_WIDE_MAX_OVER_A0 = float(max(_rw)); SR_SIGN_MIN_OVER_A0 = float(min(_rs))
print(f"\nSo both halves of this section's headline have to be quoted with their reference. The")
print(f"wide-span worst case runs {min(_rw):.4f} to {SR_WIDE_MAX_OVER_A0:.4f} decades across the four - a factor "
      f"{SR_WIDE_MAX_OVER_A0/min(_rw):.1f} - and the")
print(f"median over all nine runs {min(_rm):.4f} to {max(_rm):.4f}. The SIGN count is {min(_rs)} of {len(PAIRS)} at "
      f"{sum(s == min(_rs) for s in _rs)} of the {len(A0)} references")
print(f"and {max(_rs)} of {len(PAIRS)} at the other {sum(s == max(_rs) for s in _rs)}. The reference used above, "
      f"{A0[int(np.argmin(_rw))][0].split(' (')[0]}, is the")
print(f"KINDEST of the {len(A0)}, which is exactly why it must be named: an unattributed "
      f"'{SECOND_ROUTE_WIDE_MAX:.3f} decades' reads")
print(f"as a property of the estimator and is instead the best of a range that reaches "
      f"{SR_WIDE_MAX_OVER_A0:.2f}.")

# ---- and the same substitution applied to Result 3's ranking ---------------------------
rk_env = []
for g in PAIRS:
    k0, n0 = prior.k_barrer[g], prior.n[g]
    P_gm = RK.P_gm[g]
    dA = ((np.log10(P_gm) - np.log10(V3.k_env[g])) / V3.n_env[g]
          - (np.log10(P_gm) - np.log10(k0)) / n0)
    rk_env.append((g, dA))
rk_env.sort(key=lambda t: -t[1])
TOP3_ENV = [g for g, _ in rk_env[:3]]
HITS_ENV = len(set(TOP3_ENV) & set(SIG))
print(f"\nResult 3 gets a second route out of the same substitution, at the SAME reference:")
print(f"  ranked by selectivity gain with the ENVELOPE present line: {', '.join(g for g, _ in rk_env)}")
print(f"  the paper's three 'significant' pairs in the top three: {HITS_ENV} of 3 - unchanged, on a")
print(f"  present line that never touches Table 12.")
# But "survives it" is a robustness claim, and this page already owns a robustness test for
# exactly this quantity: EVAL, the seven references Result 3's primary route is held to. The
# substituted route gets the identical test rather than a kinder one.
print(f"\nThat is one evaluation point. Result 3's PRIMARY route is held to the {len(EVAL)} references of")
print(f"EVAL above, so the substituted route gets the same {len(EVAL)}, not a kinder test:")
_lines = [f"\n{'evaluation point':44s} {'envelope':>9s} {'primary':>9s}   top 3, envelope route"]
_he, _hp, _dropped = [], [], set()
for lab, f in EVAL:
    d = sorted(((g, ((np.log10(f(g)) - np.log10(V3.k_env[g])) / V3.n_env[g]
                     - (np.log10(f(g)) - np.log10(prior.k_barrer[g])) / prior.n[g]))
                for g in PAIRS), key=lambda z: -z[1])
    h = len(set(g for g, _ in d[:3]) & set(SIG))
    hp = rank_at(f)[0]
    _he.append(h); _hp.append(hp)
    _dropped |= set(SIG) - set(g for g, _ in d[:3])
    _lines.append(f"{lab:44s} {h:>4d} of 3 {hp:>4d} of 3   {', '.join(g for g, _ in d[:3])}")
print("\n".join(_lines))
HITS_ENV_MIN = float(min(_he))
_fail = [lab for (lab, _), h in zip(EVAL, _he) if h < max(_he)]
print(f"\nSo the honest sentence is NOT 'it survives'. The substituted route reproduces the ranking")
print(f"{HITS_ENV:.0f} of 3 at the reference used above and at {sum(h == HITS_ENV for h in _he) - 1} of the "
      f"other {len(EVAL) - 1} references, and {int(min(_he))} of 3 at the")
print(f"remaining {len(_fail)} - {'; '.join(_fail)} - where {', '.join(sorted(_dropped))} drops out of "
      f"the top three. The PRIMARY")
print(f"route holds {int(min(_hp))} of 3 at all {len(EVAL)}. That contrast is the interesting result and it is")
print(f"the reported one: the corroboration is real but weaker than the route it corroborates,")
print(f"which is what you would expect of a present line estimated from {len(PTS)} points rather than")
print(f"read off a table.")


# %% [markdown]
# ### V4: break table
#
# Every number the page states is registered below as a claim, with how it is
# checked: a **value** held to a relative tolerance of $10^{-6}$ (the page is
# arithmetic on printed numbers, so its values reproduce to round-off), a
# **bound** for identities and root-finder agreement, and an **exact** result
# for counts, each count with its distance to the threshold that decides it.
# The rows below are the misreadings and slips this page could actually have
# made: sign and convention slips, glyph and digit errors in single cells of
# Tables 12, 1–11 and 13, a lost row in each table, and two slips in the page's
# own algebra (the inverse of the bound, the closed form of the crossing).
#
# `rerun_all` recomputes **every claim** from the (possibly corrupted) inputs,
# by a code path separate from the narrative cells above; its clean call is
# asserted equal to the narrative values. It calls the same functions of
# `model.py` (`P_bound`, `alpha_bound`, `excursion`, `envelope_fit`) and the
# same `ulp` and `sqrtM`, which are written once, so it is a second
# implementation of the page's assembly and not of those primitives. Each row
# hands every claim's recomputed value to `claims.breaks`, and a claim counts as
# protected only by a row that fails its own check; the last cell of the page
# raises when a claim has none.
#
# Some pairs of row and claim are not caught, and the informative ones are exact
# invariances rather than blind spots. Swapping the prior and present blocks
# leaves $\alpha^\star$ where it is ($\alpha^\star$ is symmetric in the two
# lines) and leaves every front-factor share unchanged (numerator and
# denominator both change sign); the natural-log row leaves the shares
# unchanged too (both scale by $\ln 10$); a decimal slip in Table 10 moves
# where the data are but not where the two printed lines cross.

# %%
# Identities registered as claims: each is recomputed under every break row, and the
# natural-log row (L) is the defect that breaks all three.
def v2_rows(t13, sh, lg=np.log10):
    """The four printed prior/present pairs of Table 13 against the decomposition."""
    out = []
    for tab in ("13a", "13b"):
        sub = t13[t13.table == tab]
        for g in sub.gas_pair[sub.dataset == "prior"]:
            row_p = sub[(sub.gas_pair == g) & (sub.dataset == "present")].iloc[0]
            row_q = sub[(sub.gas_pair == g) & (sub.dataset == "prior")].iloc[0]
            if row_p.ks_D_ratio != row_q.ks_D_ratio:
                raise ValueError(f"{tab} {g}: the prior and present rows must sit at one alpha")
            a = float(row_p.ks_D_ratio)
            d = sh.dlogk[g] + sh.dn[g] * lg(a)
            printed = np.log10(row_p.P_transition_barrer / row_q.P_transition_barrer)
            out.append(dict(table=tab, pair=g, alpha=a, diff=d - printed, k_share=sh.dlogk[g] / d))
    return pd.DataFrame(out)


def rerun_all(T12x=None, PTSx=None, T13x=None, CLMx=None, sign=1.0, logfun=np.log10,
              swap_blocks=False, flip_r=False, inv=alpha_bound, cross_cf=None):
    """Recompute every claim of the page from (possibly corrupted) inputs.

    Returns ``(values, margins)``: the claimed numbers, and for the counts the
    smallest distance of a counted quantity to its threshold.
    """
    t = T12 if T12x is None else T12x
    pts = PTS if PTSx is None else PTSx
    t13 = T13 if T13x is None else T13x
    lg = logfun
    cross_cf = cross_cf or (lambda dlk, dn: 10.0 ** (-dlk / dn))
    pri = t[t.dataset == ("present" if swap_blocks else "prior")].set_index("gas_pair")
    pre = t[t.dataset == ("prior" if swap_blocks else "present")].set_index("gas_pair")
    pri = pri.assign(n=sign * pri.n)
    pre = pre.assign(n=sign * pre.n)
    pairs_all = list(t[t.dataset == "prior"].gas_pair)
    pairs = [g for g in pairs_all
             if np.isfinite(pri.k_barrer[g]) and np.isfinite(pre.k_barrer[g])]
    o, mg = {}, {}
    clm = CLM if CLMx is None else CLMx
    cls = (clm[clm.shift_class.notna() & clm.gas_pair.notna()]
           .drop_duplicates("gas_pair").set_index("gas_pair").shift_class)
    sig = [g for g in pairs if cls.get(g) == "significant"]
    for c in ("modest", "minor", "significant"):
        o[f"pairs_classified_{c}"] = int(sum(cls.get(g) == c for g in pairs))

    # ---- the transcription -----------------------------------------------------------
    o["pairs_in_table12"] = len(pairs_all)
    o["pairs_comparable"] = len(pairs)
    o["cells_table12"] = int(2 * t.k_barrer.notna().sum())
    o["cells_tables1_11"] = 2 * len(pts)
    o["cells_table13"] = 4 * len(t13)
    o["t13_rows"] = len(t13)
    o["cells_total"] = o["cells_table12"] + o["cells_tables1_11"] + o["cells_table13"]
    o["references_cited"] = len({x.strip() for r in pts.ref_bracket.astype(str)
                                 for x in r.split(",")})
    nabs = t.n.abs().dropna()
    o["n_abs_min"], o["n_abs_max"] = float(nabs.min()), float(nabs.max())
    o["n_abs_ratio"] = o["n_abs_max"] / o["n_abs_min"]

    # ---- Result 1 -------------------------------------------------------------------
    rr = []
    for g in pairs:
        a = pts[pts.gas_pair == g].alpha.to_numpy(float)
        a_gm = float(np.exp(np.mean(np.log(a))))
        dlk = lg(pre.k_barrer[g]) - lg(pri.k_barrer[g])
        dn = pre.n[g] - pri.n[g]
        s_gm = dlk + dn * lg(a_gm)
        spread = abs(dn) * lg(a.max() / a.min())
        rr.append(dict(pair=g, dlogk=dlk, dn=dn, a_lo=a.min(), a_hi=a.max(), a_gm=a_gm,
                       shift_gm=s_gm, share=dlk / s_gm, spread_over_shift=spread / abs(s_gm)))
    sh = pd.DataFrame(rr).set_index("pair")
    o["front_factor_share_median"] = float(sh.share.median())
    o["front_factor_share_min"] = float(sh.share.min())
    o["front_factor_share_max"] = float(sh.share.max())
    o["pairs_near_parallel"] = int((sh.spread_over_shift < 0.20).sum())
    mg["pairs_near_parallel"] = float(np.min(np.abs(sh.spread_over_shift - 0.20)) / 0.20)
    o["heh2_shift_decades"] = float(sh.shift_gm["He/H2"])
    o["heh2_dlogk"] = float(sh.dlogk["He/H2"])
    o["heh2_slope_change_frac"] = float(abs(sh.dn["He/H2"] / pri.n["He/H2"]))
    o["heh2_spread_over_shift"] = float(sh.spread_over_shift["He/H2"])
    o["heh2_share_gm"] = float(sh.share["He/H2"])
    o["heh2_slope_term_max_frac"] = float(abs(sh.dn["He/H2"]) * lg(sh.a_hi["He/H2"])
                                          / abs(sh.shift_gm["He/H2"]))
    o["h2n2_spread_over_shift"] = float(sh.spread_over_shift["H2/N2"])
    o["h2ch4_spread_over_shift"] = float(sh.spread_over_shift["H2/CH4"])
    o["h2co2_spread_over_shift"] = float(sh.spread_over_shift["H2/CO2"])
    o["h2co2_share_gm"] = float(sh.share["H2/CO2"])
    o["pairs_alpha1_inside_range"] = int(((sh.a_lo < 1) & (sh.a_hi > 1)).sum())
    mg["pairs_alpha1_inside_range"] = float(np.min(np.abs(np.log10(
        np.concatenate([sh.a_lo, sh.a_hi])))))

    # ---- Result 2 -------------------------------------------------------------------
    ins, dist = 0, []
    for g in pairs:
        dlk, dn = sh.dlogk[g], sh.dn[g]
        lo, hi = np.log10(sh.a_lo[g]), np.log10(sh.a_hi[g])
        ins += int((dlk + dn * lo) * (dlk + dn * hi) < 0)
        la = -dlk / dn
        dist.append(min(abs(la - lo), abs(la - hi)))
    o["pairs_crossing_inside_range"] = ins
    mg["pairs_crossing_inside_range"] = float(min(dist))
    gx = "H2/CO2"
    dlk, dn = sh.dlogk[gx], sh.dn[gx]
    a_cf = float(cross_cf(dlk, dn))
    # the root of the shift in log10(alpha): [-3, 6] decades holds every selectivity the
    # paper prints and more; brentq raises when the shift does not change sign on it
    a_rf = 10.0 ** brentq(lambda la: dlk + dn * la, -3.0, 6.0, xtol=1e-14, rtol=1e-15)
    o["h2co2_crossing_alpha"] = a_cf
    o["h2co2_crossing_alpha_rootfind_rel"] = float(abs(a_rf / a_cf - 1))
    ax = pts[pts.gas_pair == gx].alpha.to_numpy(float)
    o["h2co2_alpha_max"] = float(ax.max())
    o["h2co2_points_above_crossing"] = int((ax > a_cf).sum())
    mg["h2co2_points_above_crossing"] = float(np.min(np.abs(np.log10(ax / a_cf))))
    o["h2co2_present_over_prior_at_alpha_max_pct"] = float(100 * (
        pre.k_barrer[gx] * ax.max() ** pre.n[gx] / (pri.k_barrer[gx] * ax.max() ** pri.n[gx]) - 1))
    o["h2co2_present_below_prior_at_alpha_max_pct"] = -o["h2co2_present_over_prior_at_alpha_max_pct"]

    # ---- Result 3 -------------------------------------------------------------------
    rr = []
    for g in pairs:
        P = pts[pts.gas_pair == g].P_fast_barrer.to_numpy(float)
        P_gm = float(np.exp(np.mean(np.log(P))))
        dA = ((lg(P_gm) - lg(pre.k_barrer[g])) / pre.n[g]
              - (lg(P_gm) - lg(pri.k_barrer[g])) / pri.n[g])
        a1 = brentq(lambda a, g=g: pre.k_barrer[g] * a ** pre.n[g] - P_gm, 1e-6, 1e9,
                    xtol=1e-14, rtol=8.9e-16)
        a0 = brentq(lambda a, g=g: pri.k_barrer[g] * a ** pri.n[g] - P_gm, 1e-6, 1e9,
                    xtol=1e-14, rtol=8.9e-16)
        rr.append(dict(pair=g, dlogk=sh.dlogk[g], dlogP=sh.shift_gm[g], dlogA=dA,
                       dlogA_rf=np.log10(a1) - np.log10(a0), alpha_gain=10 ** dA))
    rk = pd.DataFrame(rr).set_index("pair")

    def top3(col):
        v = rk[col].sort_values(ascending=False)
        return len(set(v.index[:3]) & set(sig)), float(v.iloc[2] - v.iloc[3])
    for col, name in (("dlogA", "dlogA_top3_hits"), ("dlogk", "dlogk_top3_hits"),
                      ("dlogP", "dlogP_top3_hits")):
        o[name], mg[name] = top3(col)
    o["dlogA_gap_third_to_fourth"] = mg["dlogA_top3_hits"]
    o["dlogA_rootfind_max_abs"] = float(np.max(np.abs(rk.dlogA - rk.dlogA_rf)))
    o["alpha_gain_hech4"] = float(rk.alpha_gain["He/CH4"])
    o["alpha_gain_co2ch4"] = float(rk.alpha_gain["CO2/CH4"])
    _gm = float(np.exp(np.mean(np.log(pts.P_fast_barrer.to_numpy(float)))))
    _own = lambda f: (lambda g: float(f(pts[pts.gas_pair == g].P_fast_barrer.to_numpy(float))))  # noqa: E731
    evals = ([_own(np.min), _own(lambda P: np.exp(np.mean(np.log(P)))), _own(np.median),
              _own(np.max), lambda g: _gm] + [(lambda g, c=c: float(c)) for c in (100.0, 1000.0)])

    def hits_at(f, kfun, nfun):
        d = sorted(((g, (lg(f(g)) - lg(kfun(g))) / nfun(g)
                     - (lg(f(g)) - lg(pri.k_barrer[g])) / pri.n[g]) for g in pairs),
                   key=lambda z: -z[1])
        return len(set(g for g, _ in d[:3]) & set(sig)), d[2][1] - d[3][1]
    _h, _gp = zip(*(hits_at(f, lambda g: pre.k_barrer[g], lambda g: pre.n[g]) for f in evals))
    o["dlogA_top3_hits_min_over_eval_points"] = int(min(_h))
    mg["dlogA_top3_hits_min_over_eval_points"] = float(min(_gp))
    o["dlogA_gap_min_over_eval_points"] = float(min(_gp))

    # ---- Result 4 -------------------------------------------------------------------
    er = []
    for g in pairs_all:
        if g not in set(pts.gas_pair):
            continue
        sub = pts[pts.gas_pair == g]
        r1 = np.asarray(excursion(sub.P_fast_barrer, sub.alpha, pre.k_barrer[g], pre.n[g]), float)
        r0 = (np.asarray(excursion(sub.P_fast_barrer, sub.alpha,
                                   pri.k_barrer[g], pri.n[g]), float)
              if np.isfinite(pri.k_barrer[g]) else np.full(len(sub), np.nan))
        if flip_r:                       # one inverted function, used in both places
            r1, r0 = 1.0 / r1, 1.0 / r0
        er += [dict(gas_pair=g, r_present=a, r_prior=b) for a, b in zip(r1, r0)]
    ex = pd.DataFrame(er)
    o["points_total"] = len(ex)
    o["points_above_present_bound"] = int((ex.r_present > 1).sum())
    mg["points_above_present_bound"] = float(np.min(np.abs(ex.r_present - 1)))
    o["frac_above_present_bound"] = o["points_above_present_bound"] / len(ex)
    _o2 = ex[ex.gas_pair == "O2/N2"].r_present
    o["o2n2_points_above_present_bound"] = int((_o2 > 1).sum())
    mg["o2n2_points_above_present_bound"] = float(np.min(np.abs(_o2 - 1)))
    o["excursion_max"] = float(ex.r_present.max())
    o["excursion_min"] = float(ex.r_present.min())
    o["excursion_geomean"] = float(np.exp(np.mean(np.log(ex.r_present))))
    _rp = ex[ex.gas_pair.isin(pairs)].r_prior
    o["points_above_prior_bound"] = int((_rp > 1).sum())
    mg["points_above_prior_bound"] = float(np.min(np.abs(_rp - 1)))

    # ---- V1 -------------------------------------------------------------------------
    sq = np.array([sqrtM(g) for g in t13.gas_pair])
    o["t13_mw_column_max_abs_rel"] = float(np.max(np.abs(sq / t13.sqrt_M_ratio - 1)))
    sq_iu = np.array([np.sqrt(IUPAC[g.split("/")[1]] / IUPAC[g.split("/")[0]]) for g in t13.gas_pair])
    o["t13_mw_column_max_abs_rel_iupac"] = float(np.max(np.abs(sq_iu / t13.sqrt_M_ratio - 1)))
    o["t13_mw_h2n2_printed"] = float(t13[(t13.table == "13a") & (t13.gas_pair == "H2/N2")]
                                     .sqrt_M_ratio.iloc[0])
    o["t13_product_column_max_abs_rel"] = float(np.max(np.abs(
        t13.ks * t13.sqrt_M_ratio / t13.ks_D_ratio - 1)))
    vr = []
    for _, r in t13.iterrows():
        src = pri if r.dataset == "prior" else pre
        k, n = src.k_barrer[r.gas_pair], src.n[r.gas_pair]
        kp, npr = src.k_as_printed[r.gas_pair], src.n_as_printed[r.gas_pair]
        rel = k * r.ks_D_ratio ** n / r.P_transition_barrer - 1
        common = (abs(np.log(r.ks_D_ratio)) * ulp(npr)
                  + abs(n) * ulp(r.ks_D_as_printed) / r.ks_D_ratio)
        s4 = lambda v: 0.5 * 10.0 ** (np.floor(np.log10(abs(v))) - 3)     # noqa: E731
        vr.append(dict(key=(r.table, r.gas_pair, r.dataset), rel=rel,
                       bt=ulp(kp) / k + common + ulp(r.P_transition_barrer) / r.P_transition_barrer,
                       bl=s4(k) / k + common + s4(r.P_transition_barrer) / r.P_transition_barrer))
    v1 = pd.DataFrame(vr)
    o["t13_identity_max_abs_rel"] = float(v1.rel.abs().max())
    o["t13_identity_median_abs_rel"] = float(v1.rel.abs().median())
    o["t13_rows_inside_tight_band"] = int((v1.rel.abs() <= v1.bt).sum())
    mg["t13_rows_inside_tight_band"] = float(np.min(np.abs(v1.rel.abs() / v1.bt - 1)))
    o["t13_rows_inside_loose_band"] = int((v1.rel.abs() <= v1.bl).sum())
    mg["t13_rows_inside_loose_band"] = float(np.min(np.abs(v1.rel.abs() / v1.bl - 1)))
    v1k = v1.set_index("key")
    o["t13a_o2n2_rel"] = float(v1k.rel[("13a", "O2/N2", "present")])
    o["t13b_h2n2_tight_excess"] = float(v1k.rel[("13b", "H2/N2", "present")].__abs__()
                                        / v1k.bt[("13b", "H2/N2", "present")] - 1)
    _rn = t13[(t13.table == "13a") & (t13.gas_pair == "N2/CH4")].iloc[0]
    o["t13a_n2ch4_band_tight"] = float(v1k.bt[("13a", "N2/CH4", "present")])
    o["t13a_n2ch4_band_tight_float_read"] = float(
        v1k.bt[("13a", "N2/CH4", "present")]
        + abs(pre.n["N2/CH4"]) * (ulp(repr(float(_rn.ks_D_ratio))) - ulp(_rn.ks_D_as_printed))
        / _rn.ks_D_ratio)
    _heh2 = v1[[k[1] == "He/H2" for k in v1.key]]
    o["v1_heh2_rows_inside_tight_band"] = int((_heh2.rel.abs() <= _heh2.bt).sum())
    mg["v1_heh2_rows_inside_tight_band"] = float(np.min(np.abs(_heh2.rel.abs() / _heh2.bt - 1)))
    o["t13a_o2n2_band_loose"] = float(v1k.bl[("13a", "O2/N2", "present")])
    o["t13b_o2n2_rel"] = float(v1k.rel[("13b", "O2/N2", "present")])
    _r = t13[(t13.table == "13a") & (t13.gas_pair == "O2/N2")].iloc[0]
    o["t13a_o2n2_alpha_printed"] = float(_r.ks_D_ratio)
    o["t13a_o2n2_alpha_needed"] = float(inv(_r.P_transition_barrer, pre.k_barrer["O2/N2"],
                                            pre.n["O2/N2"]))
    o["t13a_o2n2_alpha_backout_rel"] = o["t13a_o2n2_alpha_needed"] / o["t13a_o2n2_alpha_printed"] - 1
    o["_v1_heh2"] = float(abs(v1k.rel[("13a", "He/H2", "present")]))
    o["_v1_h2co2_prior"] = float(abs(v1k.rel[("13a", "H2/CO2", "prior")]))

    # ---- V2 and the two identities --------------------------------------------------
    v2 = v2_rows(t13, sh, lg)
    o["v2_shift_crosscheck_max_abs_decades"] = float(v2["diff"].abs().max())
    o["v2_kshare_min"] = float(v2.k_share.min())
    o["v2_kshare_max"] = float(v2.k_share.max())
    for tab in ("13a", "13b"):
        for g, nm in (("He/H2", "heh2"), ("H2/CO2", "h2co2")):
            o[f"v2_kshare_{nm}_{tab}"] = float(v2[(v2.table == tab) & (v2.pair == g)].k_share.iloc[0])
    o["v2_is_v1_difference_max_abs"] = float(max(
        abs(r["diff"] - np.log10((1 + v1k.rel[(r.table, r.pair, "present")])
                                 / (1 + v1k.rel[(r.table, r.pair, "prior")])))
        for _, r in v2.iterrows()))
    o["dataside_shift_is_decomposition_max_abs"] = float(max(
        abs(np.mean(np.log10(ex[ex.gas_pair == g].r_prior))
            - np.mean(np.log10(ex[ex.gas_pair == g].r_present)) - sh.shift_gm[g])
        for g in pairs))

    # ---- V3, and the second route it carries ------------------------------------------
    v3 = []
    for g in sorted(pts.gas_pair.unique()):
        sub = pts[pts.gas_pair == g]
        la = np.log10(sub.alpha.to_numpy(float))
        lp = np.log10(sub.P_fast_barrer.to_numpy(float))
        k_env, n_env, _ = envelope_fit(la, lp)
        a_gm = float(np.exp(np.mean(np.log(sub.alpha.to_numpy(float)))))
        off = ((np.log10(k_env) + n_env * np.log10(a_gm))
               - (np.log10(pre.k_barrer[g]) + pre.n[g] * np.log10(a_gm)))
        v3.append(dict(pair=g, n_env=n_env, k_env=k_env, dn=n_env - pre.n[g],
                       span=la.max() - la.min(), offset_gm=off))
    v3 = pd.DataFrame(v3).set_index("pair")
    o["envelope_dn_median_abs"] = float(v3.dn.abs().median())
    o["envelope_dn_max_abs"] = float(v3.dn.abs().max())
    o["envelope_offset_median_decades"] = float(v3.offset_gm.median())
    o["envelope_offset_max_abs_decades"] = float(v3.offset_gm.abs().max())
    o["envelope_offsets_positive"] = int((v3.offset_gm > 0).sum())
    mg["envelope_offsets_positive"] = float(v3.offset_gm.abs().min())
    o["envelope_offsets_negative"] = int((v3.offset_gm < 0).sum())
    mg["envelope_offsets_negative"] = mg["envelope_offsets_positive"]
    o["envelope_dn_max_abs_wide_span"] = float(v3[v3.span > 1.5].dn.abs().max())
    o["envelope_span_min_decades"] = float(v3.span.min())
    o["envelope_dn_positive"] = int((v3.dn > 0).sum())
    mg["envelope_dn_positive"] = float(v3.dn.abs().min())
    o["envelope_dn_negative"] = int((v3.dn < 0).sum())
    mg["envelope_dn_negative"] = mg["envelope_dn_positive"]
    o["pairs_span_over_1p5_decades"] = int((v3.span > 1.5).sum())
    mg["pairs_span_over_1p5_decades"] = float(np.min(np.abs(v3.span - 1.5)))
    by_span = list(v3.span.sort_values().index)
    by_dn = list(v3.dn.abs().sort_values(ascending=False).index)
    o["lever_arm_set_size"] = max(m for m in range(1, len(v3)) if set(by_span[:m]) == set(by_dn[:m]))
    # the second-route residual: sh_env minus the decomposition's own shift, on the
    # comparable pairs; at the clean point it equals v3.offset_gm (an identity claim)
    sr = []
    for g in pairs:
        dlk_env = np.log10(v3.k_env[g]) - np.log10(pri.k_barrer[g])
        sh_env = dlk_env + (v3.n_env[g] - pri.n[g]) * np.log10(sh.a_gm[g])
        sr.append(dict(pair=g, dlk_env=dlk_env, sh_env=sh_env, share_env=dlk_env / sh_env,
                       resid=sh_env - sh.shift_gm[g]))
    sr = pd.DataFrame(sr).set_index("pair")
    o["second_route_is_offset_max_abs"] = float(np.max(np.abs(sr.resid - v3.offset_gm[sr.index])))
    o["envelope_offset_max_abs_wide_span"] = float(sr.resid[v3.span[sr.index] > 1.5].abs().max())
    o["sr_resid_hech4"] = float(abs(sr.resid["He/CH4"]))
    o["envelope_implied_share_median"] = float(sr.share_env.median())
    for g, nm in (("H2/CO2", "h2co2"), ("He/CO2", "heco2"), ("H2/CH4", "h2ch4")):
        o[f"{nm}_share_env"] = float(sr.share_env.get(g, np.nan))
    o["heco2_share_gm"] = float(sh.share.get("He/CO2", np.nan))
    o["h2ch4_share_gm"] = float(sh.share["H2/CH4"])
    o["envelope_implied_share_min"] = float(sr.share_env.min())
    o["envelope_implied_share_max"] = float(sr.share_env.max())
    o["envelope_hen2_dlogk"] = float(sr.dlk_env["He/N2"])
    o["envelope_hen2_shift"] = float(sr.sh_env["He/N2"])
    o["envelope_hen2_cross_alpha"] = float(
        10.0 ** (-sr.dlk_env["He/N2"] / (v3.n_env["He/N2"] - pri.n["He/N2"])))
    a0s = [lambda g: float(pts[pts.gas_pair == g].alpha.min()), lambda g: sh.a_gm[g],
           lambda g: float(np.median(pts[pts.gas_pair == g].alpha.to_numpy(float))),
           lambda g: float(pts[pts.gas_pair == g].alpha.max())]
    _w, _s, _sm = [], [], []
    for af in a0s:
        wide, sgn = [], 0
        for g in pairs:
            a0 = af(g)
            se = ((np.log10(v3.k_env[g]) - np.log10(pri.k_barrer[g]))
                  + (v3.n_env[g] - pri.n[g]) * np.log10(a0))
            sd = sh.dlogk[g] + sh.dn[g] * lg(a0)
            sgn += int(np.sign(se) == np.sign(sd))
            _sm.append(min(abs(se), abs(sd)))
            if v3.span[g] > 1.5:
                wide.append(abs(se - sd))
        _w.append(max(wide)); _s.append(sgn)
    o["envelope_offset_max_abs_wide_span_max_over_alpha0"] = float(max(_w))
    o["envelope_sign_agree_min_over_alpha0"] = int(min(_s))
    mg["envelope_sign_agree_min_over_alpha0"] = float(min(_sm))
    o["envelope_sign_agree_at_alpha_gm"] = int(_s[1])
    mg["envelope_sign_agree_at_alpha_gm"] = float(min(min(abs(sr.sh_env)), min(abs(sh.shift_gm))))
    _he, _hg = zip(*(hits_at(f, lambda g: v3.k_env[g], lambda g: v3.n_env[g]) for f in evals))
    o["envelope_implied_dlogA_top3_hits"] = int(_he[1])
    mg["envelope_implied_dlogA_top3_hits"] = float(abs(_hg[1]))
    o["envelope_implied_dlogA_top3_hits_min_over_eval_points"] = int(min(_he))
    mg["envelope_implied_dlogA_top3_hits_min_over_eval_points"] = float(min(map(abs, _hg)))
    o["envelope_eval_points_with_3_hits"] = int(sum(h == 3 for h in _he))
    mg["envelope_eval_points_with_3_hits"] = float(min(map(abs, _hg)))

    o["bound_roundtrip_max_rel"] = roundtrip_max(pre, pairs, inv)
    return o, mg


# ---- the rows -----------------------------------------------------------------------
ROWS = [("(S) n taken POSITIVE (sign convention flipped)", dict(sign=-1.0)),
        ("(B) prior and present blocks swapped", dict(swap_blocks=True)),
        ("(L) natural log used in the decomposition", dict(logfun=np.log)),
        ("(I) excursion inverted: k*alpha^n / P", dict(flip_r=True)),
        ("(N) inverse of the bound written with n instead of 1/n",
         dict(inv=lambda P, k, n: (np.asarray(P, float) / k) ** n)),
        ("(X) closed form of the crossing with its sign dropped",
         dict(cross_cf=lambda dlk, dn: 10.0 ** (dlk / dn)))]
for tag, label, (ds, g, col, val, printed) in [
        ("C", "thousands comma dropped: k(O2/N2, present) 1,396,000 -> 1396",
         ("present", "O2/N2", "k_barrer", 1396.0, "1396")),
        ("G9", "glyph 9<->0: k(He/H2, present) 59,910 -> 59,010",
         ("present", "He/H2", "k_barrer", 59010.0, "59,010")),
        ("G8", "glyph 8<->3: n(He/CH4, present) -0.809 -> -0.309",
         ("present", "He/CH4", "n", -0.309, "-0.309")),
        ("G0", "glyph 8<->0: n(O2/N2, prior) -5.800 -> -5.000",
         ("prior", "O2/N2", "n", -5.0, "-5.000")),
        ("DD", "digit dropped: k(He/CH4, present) 19,800 -> 1,980",
         ("present", "He/CH4", "k_barrer", 1980.0, "1,980")),
        ("D1", "last digit: n(He/H2, present) -4.864 -> -4.865",
         ("present", "He/H2", "n", -4.865, "-4.865")),
        ("D2", "last digit: k(H2/CO2, prior) 1200 -> 1201",
         ("prior", "H2/CO2", "k_barrer", 1201.0, "1201")),
        ("DP", "digit dropped in the PRIOR block: k(He/H2, prior) 960 -> 96",
         ("prior", "He/H2", "k_barrer", 96.0, "96")),
        ("G7", "glyph 7<->1: n(N2/CH4, present) -4.507 -> -4.501",
         ("present", "N2/CH4", "n", -4.501, "-4.501"))]:
    t = T12.copy()
    m = (t.dataset == ds) & (t.gas_pair == g)
    if int(m.sum()) != 1:
        raise ValueError(f"{tag} must hit exactly one Table 12 row")
    t.loc[m, col] = val
    t.loc[m, "k_as_printed" if col == "k_barrer" else "n_as_printed"] = printed
    ROWS.append((f"({tag}) {label}", dict(T12x=t)))
for tag, g in (("RT", "He/CO2"),):
    ROWS.append((f"({tag}) Table 12 row lost: the prior {g} row",
                 dict(T12x=T12[~((T12.dataset == "prior") & (T12.gas_pair == g))])))
# the Tables 1-11 rows: one alpha that carries a headline, one P, one that sets a span,
# and one row lost (the last row whose reference no other row cites)
for tag, label, (gp, col, old, new) in [
        ("T1", "TABLE 10 decimal slip: alpha(H2/CO2) 100.9 -> 10.09",
         ("H2/CO2", "alpha", 100.9, 10.09)),
        ("T2", "TABLE 1 decimal slip: P(O2/N2) 18.0 -> 1.80",
         ("O2/N2", "P_fast_barrer", 18.0, 1.80)),
        ("T3", "TABLE 7 glyph 3<->8: alpha(He/H2) 4.39 -> 4.89",
         ("He/H2", "alpha", 4.39, 4.89))]:
    p = PTS.copy()
    m = (p.gas_pair == gp) & (p[col] == old)
    if int(m.sum()) != 1:
        raise ValueError(f"{tag} must hit exactly one cell, hit {int(m.sum())}")
    p.loc[m, col] = new
    ROWS.append((f"({tag}) {label}", dict(PTSx=p)))
_i9 = PTS[PTS.gas_pair == "N2/CH4"].alpha.idxmax()
ROWS.append((f"(T4) TABLE 9 decimal slip: the largest alpha(N2/CH4) read ten times too large",
             dict(PTSx=PTS.assign(alpha=PTS.alpha.where(PTS.index != _i9, 10 * PTS.alpha[_i9])))))
ROWS.append(("(CL) Conclusions misread: H2/CO2 classed 'modest' instead of 'minor'",
             dict(CLMx=CLM.assign(shift_class=CLM.shift_class.where(
                 ~((CLM.gas_pair == "H2/CO2") & CLM.shift_class.notna()), "modest")))))
_refs = PTS.ref_bracket.astype(str)
_lost = PTS.index[_refs.map(_refs.value_counts()) == 1][-1]
ROWS.append((f"(RP) Tables 1-11 row lost: {PTS.gas_pair[_lost]} {PTS.polymer[_lost][:30]}",
             dict(PTSx=PTS.drop(index=_lost))))
# the Table 13 rows: two molecular-weight cells, one product cell, one row lost
for tag, label, (tb, gp, col, new, printed) in [
        ("M1", "TABLE 13a glyph 42<->24: sqrt(M) H2/N2 3.742 -> 3.724",
         ("13a", "H2/N2", "sqrt_M_ratio", 3.724, None)),
        ("M2", "TABLE 13a glyph 90<->09: sqrt(M) H2/CO2 4.6904 -> 4.6094",
         ("13a", "H2/CO2", "sqrt_M_ratio", 4.6094, None)),
        ("K", "TABLE 13a glyph 73<->37: ks*sqrt(M) O2/N2 1.473 -> 1.437",
         ("13a", "O2/N2", "ks_D_ratio", 1.437, "1.437"))]:
    q = T13.copy()
    m = (q.table == tb) & (q.gas_pair == gp) & (q.dataset == "present")
    if int(m.sum()) != 1:
        raise ValueError(f"{tag} must hit exactly one Table 13 row")
    q.loc[m, col] = new
    if printed is not None:
        q.loc[m, "ks_D_as_printed"] = printed
    ROWS.append((f"({tag}) {label}", dict(T13x=q)))
ROWS.append(("(R13) Table 13 row lost: 13b N2/CH4",
             dict(T13x=T13[~((T13.table == "13b") & (T13.gas_pair == "N2/CH4"))])))

BASE, MARGIN = rerun_all()
BROKEN = {lab: rerun_all(**kw)[0] for lab, kw in ROWS}
BR = pd.DataFrame([BASE] + [BROKEN[lab] for lab, _ in ROWS],
                  index=["(clean) - as transcribed"] + [lab for lab, _ in ROWS])

SHOW = ["heh2_shift_decades", "t13_identity_max_abs_rel", "_v1_heh2", "_v1_h2co2_prior",
        "v2_shift_crosscheck_max_abs_decades", "h2co2_crossing_alpha",
        "envelope_offset_max_abs_decades", "points_above_present_bound", "dlogA_top3_hits"]
HEAD = ["He/H2 shift", "V1 max", "V1 He/H2", "V1 H2/CO2pr", "V2 max", "alpha*",
        "V3 offset", "#above", "top3"]


def _fmt(v):
    return f"{v:11.4f}" if abs(v) >= 1e-3 or v == 0 else f"{v:11.2e}"


_W = max(len(lab) for lab in BR.index)
print("\n".join([f"{'injected defect':{_W}s} " + " ".join(f"{h:>11s}" for h in HEAD)]
                + [f"{lab:{_W}s} " + " ".join(_fmt(BR[c][lab]) for c in SHOW) for lab in BR.index]))

# The clean rerun must reproduce the narrative cells, by its own code path.
NARRATIVE = dict(
    pairs_in_table12=N_PAIRS_TABLE, pairs_comparable=N_PAIRS_COMPARABLE,
    cells_table12=N_NUM_T12, cells_tables1_11=N_NUM_PTS, cells_table13=N_NUM_T13,
    cells_total=N_NUMERALS,
    front_factor_share_median=SHARE_MED, front_factor_share_min=SHARE_MIN,
    front_factor_share_max=SHARE_MAX, pairs_near_parallel=PARALLEL_N,
    heh2_shift_decades=SH.shift_gm["He/H2"], heh2_dlogk=SH.dlogk["He/H2"],
    heh2_spread_over_shift=SH.spread_over_shift["He/H2"], pairs_alpha1_inside_range=len(IN1),
    pairs_crossing_inside_range=N_CROSSING_INSIDE, h2co2_crossing_alpha_rootfind_rel=ALPHA_STAR_AGREE,
    envelope_dn_positive=N_DN_POS, o2n2_points_above_present_bound=N_ABOVE_O2N2, h2co2_crossing_alpha=ALPHA_STAR,
    h2co2_points_above_crossing=len(above),
    h2co2_present_over_prior_at_alpha_max_pct=PRESENT_BELOW_PRIOR_PCT,
    dlogA_top3_hits=HITS_A, dlogk_top3_hits=HITS_K, dlogP_top3_hits=HITS_P,
    dlogA_gap_third_to_fourth=gap_A, dlogA_top3_hits_min_over_eval_points=EVAL_HITS_MIN,
    dlogA_gap_min_over_eval_points=EVAL_GAP_MIN, alpha_gain_hech4=RK.alpha_gain["He/CH4"],
    alpha_gain_co2ch4=RK.alpha_gain["CO2/CH4"], points_total=N_POINTS,
    points_above_present_bound=N_ABOVE, frac_above_present_bound=FRAC_ABOVE,
    excursion_max=R_MAX, excursion_min=R_MIN, excursion_geomean=R_GEOMEAN,
    points_above_prior_bound=N_ABOVE_PRIOR, t13_mw_column_max_abs_rel=MW_MAX_REL,
    t13_mw_column_max_abs_rel_iupac=MW_MAX_REL_IUPAC, t13_product_column_max_abs_rel=PROD_MAX_REL,
    t13_identity_max_abs_rel=T13_MAX_REL, t13_identity_median_abs_rel=T13_MED_REL,
    t13_rows_inside_tight_band=N_INSIDE_TIGHT, t13_rows_inside_loose_band=N_INSIDE_LOOSE,
    t13a_o2n2_alpha_backout_rel=T13A_O2N2_BACKOUT_REL,
    v2_shift_crosscheck_max_abs_decades=V2_MAX_ABS, v2_kshare_min=V2.k_share.min(),
    v2_kshare_max=V2.k_share.max(), envelope_dn_median_abs=V3_DN_MED,
    envelope_dn_max_abs=V3_DN_MAX, envelope_offset_median_decades=V3_OFFSET_MED,
    envelope_offset_max_abs_decades=V3_OFFSET_MAX, envelope_offsets_positive=N_OFFSET_POS,
    envelope_offset_max_abs_wide_span=SECOND_ROUTE_WIDE_MAX,
    envelope_offset_max_abs_wide_span_max_over_alpha0=SR_WIDE_MAX_OVER_A0,
    envelope_sign_agree_min_over_alpha0=SR_SIGN_MIN_OVER_A0,
    envelope_sign_agree_at_alpha_gm=SECOND_ROUTE_SIGN_AGREE,
    envelope_dn_max_abs_wide_span=V3_DN_WIDE, envelope_span_min_decades=V3_SPAN_MIN,
    lever_arm_set_size=M_LEVER, envelope_implied_share_median=SHARE_ENV_MED,
    envelope_implied_share_min=SR.share_env.min(), envelope_implied_share_max=SR.share_env.max(),
    envelope_implied_dlogA_top3_hits=HITS_ENV,
    envelope_implied_dlogA_top3_hits_min_over_eval_points=HITS_ENV_MIN)
_bad = {k: (v, BASE[k]) for k, v in NARRATIVE.items()
        if not abs(float(BASE[k]) - float(v)) <= 1e-12 * max(1.0, abs(float(v)))}
if _bad:
    raise AssertionError(f"rerun_all disagrees with the narrative cells: {_bad}")
print(f"\nThe clean rerun reproduces all {len(NARRATIVE)} narrative values it shares with the cells "
      "above, to 1e-12 relative.")

# %%
# Reading the table row by row. Looked up by label, never by list index.
def pick(frag):
    """The recomputed claims under the one row whose label contains ``frag``."""
    lab, = [lab for lab in BROKEN if frag in lab]
    return BROKEN[lab]


o_ld, o_pk = pick("n(He/H2, present) -4.864"), pick("k(H2/CO2, prior) 1200")
o_g8, o_dd = pick("n(He/CH4, present) -0.809"), pick("k(He/CH4, present) 19,800")
o_t1, o_t2, o_t3 = pick("alpha(H2/CO2) 100.9"), pick("P(O2/N2) 18.0"), pick("alpha(He/H2) 4.39")
o_b = pick("blocks swapped")
_dd_t12, = [kw["T12x"] for lab, kw in ROWS if "k(He/CH4, present) 19,800" in lab]
o_dd_k = float(_dd_t12[(_dd_t12.dataset == "present") & (_dd_t12.gas_pair == "He/CH4")].k_barrer.iloc[0])
print("\n".join([
    "The block swap leaves the crossing selectivity alpha* where it was "
    f"({BASE['h2co2_crossing_alpha']:.4f} -> {o_b['h2co2_crossing_alpha']:.4f}): alpha* is symmetric in",
    "the two lines, so exchanging them cannot move it. Every Table 12 row corrupts one cell and",
    "moves only what that cell feeds; the Tables 1-11 rows corrupt one point each.",
    "",
    "The two LAST-DIGIT rows are the resolving-power rows, and V1's GLOBAL max is the wrong column",
    "to read them in (it stays pinned to the unrelated 13a O2/N2 row at "
    f"{BASE['t13_identity_max_abs_rel']:.2e}). Per row:",
    f"  n(He/H2) -4.864 -> -4.865 moves the He/H2 identity residual {BASE['_v1_heh2']:.2e} -> "
    f"{o_ld['_v1_heh2']:.2e} ({o_ld['_v1_heh2'] / BASE['_v1_heh2']:.0f}x)",
    "                            and the V2 cross-check "
    f"{BASE['v2_shift_crosscheck_max_abs_decades']:.2e} -> {o_ld['v2_shift_crosscheck_max_abs_decades']:.2e} "
    f"({o_ld['v2_shift_crosscheck_max_abs_decades'] / BASE['v2_shift_crosscheck_max_abs_decades']:.1f}x)",
    f"  k(H2/CO2,prior) 1200 -> 1201 moves that row's identity residual {BASE['_v1_h2co2_prior']:.2e} -> "
    f"{o_pk['_v1_h2co2_prior']:.2e} ({o_pk['_v1_h2co2_prior'] / BASE['_v1_h2co2_prior']:.0f}x)",
    f"                            and the crossing selectivity {BASE['h2co2_crossing_alpha']:.4f} -> "
    f"{o_pk['h2co2_crossing_alpha']:.4f}",
    "Both residuals move, but the He/H2 rows stay inside their literal printing bands "
    f"({BASE['v1_heh2_rows_inside_tight_band']} -> {o_ld['v1_heh2_rows_inside_tight_band']} of 4 inside),",
    "so V1 and V2 do NOT detect a last-digit error in n(He/H2) on their own; the regression check",
    "against the committed values does. The crossing result is sensitive to the final digit of",
    "the prior k it depends on.",
    "",
    "The ranking count is deliberately the coarsest quantity here. It survives both last-digit",
    f"rows, and the 8<->3 glyph in n(He/CH4) too ({BASE['dlogA_top3_hits']} -> {o_g8['dlogA_top3_hits']}), "
    "because a shallower slope only",
    "pushes He/CH4 further UP the selectivity-gain ranking. What breaks it is a lost digit:",
    f"k(He/CH4) 19,800 -> 1,980 takes the count {BASE['dlogA_top3_hits']} -> {o_dd['dlogA_top3_hits']}, "
    "because He/CH4's shift then reverses sign.",
    "So the ranking rests on the DECADE of every present k, which V1 pins hard.",
    "",
    "A single slip in ONE cell of Tables 1-11:",
    f"  (T1) alpha(H2/CO2) 100.9 -> 10.09 takes pairs_crossing_inside_range "
    f"{BASE['pairs_crossing_inside_range']} -> {o_t1['pairs_crossing_inside_range']}: with that point gone",
    f"       from the top of Table 10's range, alpha* = {BASE['h2co2_crossing_alpha']:.2f} is no longer "
    "inside it. No Table 12 check can see this.",
    f"  (T2) P(O2/N2) 18.0 -> 1.80 takes points_above_present_bound {BASE['points_above_present_bound']} -> "
    f"{o_t2['points_above_present_bound']}, and the envelope",
    f"       offset {BASE['envelope_offset_max_abs_decades']:.4f} -> "
    f"{o_t2['envelope_offset_max_abs_decades']:.4f} decades, because that point pins O2/N2's covering line.",
    f"  (T3) alpha(He/H2) 4.39 -> 4.89 moves envelope_span_min_decades "
    f"{BASE['envelope_span_min_decades']:.4f} -> {o_t3['envelope_span_min_decades']:.4f},",
    "       the lever arm behind the whole V3 slope discussion.",
    "",
    "(DD) answers the obvious objection to the second route. Tables 1-11 are the points Robeson",
    "SELECTED as close to his 2008 line, so the baseline residual is an optimistic figure. That",
    "says nothing about the route's DETECTING power: k_env never touches Table 12, so",
    "d(resid)/d(log k_present) = -1 identically. (DD) divides k(He/CH4, present) by ten and moves",
    f"that pair's two-route residual {BASE['sr_resid_hech4']:.5f} -> {o_dd['sr_resid_hech4']:.5f}, by "
    f"{o_dd['sr_resid_hech4'] - BASE['sr_resid_hech4']:+.5f} decades for a",
    f"{np.log10(present.k_barrer['He/CH4'] / o_dd_k):.3f}-decade error in the printed line."]))

# %%
# Every number the page states, registered once (docs/claims-design.md). Values are
# arithmetic on printed numbers and reproduce to round-off, so they are held to 1e-6.
VAL = 1e-6
EXACT_COUNTS = ["pairs_in_table12", "pairs_comparable", "cells_table12", "cells_tables1_11",
                "cells_table13", "cells_total", "t13_rows", "references_cited", "points_total",
                "lever_arm_set_size"]
BOUNDS = {
    "bound_roundtrip_max_rel": ("identity", "P_bound and alpha_bound are algebraic inverses"),
    "v2_is_v1_difference_max_abs": ("identity", "V2's residual is the difference of two V1 residuals"),
    "dataside_shift_is_decomposition_max_abs": ("identity", "mean(log alpha) is log(alpha_gm)"),
    "second_route_is_offset_max_abs": ("identity", "the prior line cancels from the two routes"),
    "h2co2_crossing_alpha_rootfind_rel": ("tolerance", "brentq xtol 1e-14 on log10(alpha)"),
    "dlogA_rootfind_max_abs": ("tolerance", "brentq rtol 8.9e-16 on alpha"),
}
for name, x in BASE.items():
    if name.startswith("_"):
        continue
    if name in BOUNDS:
        basis, note = BOUNDS[name]
        claims.bound(name, x, max=1e-12, basis=basis, note=note)
    elif name in MARGIN:
        claims.exact(name, x, margin=MARGIN[name], min_margin=1e-6,
                     note="a count decided against a threshold; margin is its distance to it")
    elif name in EXACT_COUNTS:
        claims.exact(name, x)
    else:
        claims.value(name, x, rtol=VAL)

# every row, every claim: the claim's own check decides whether the row catches it
for lab, _ in ROWS:
    for name in claims.metrics:
        claims.breaks(lab, name, BROKEN[lab][name])

# %% [markdown]
# ### Blind spots: what this page does *not* establish
#
# 1. **Nothing here validates the upper bound against measurement.** Every
#    number is Robeson's: his compiled points, his hand-drawn lines, his derived
#    tables. The page tests *internal consistency and the meaning of his
#    claims*, and that is all. In the gallery's vocabulary this is
#    **reproduction, not validation**.
# 2. **The 117 points are a compilation.** They cite 44 distinct primary
#    references, none of which was consulted. An error made in compiling them is
#    inherited silently.
# 3. **The prior bounds are a restatement.** The 1991 and 1994 papers are not on
#    disk. If the 2008 paper misquotes its own earlier $k$ or $n$, every shift
#    on this page is wrong by exactly that amount and nothing here would show
#    it. V1 touches the prior block for only two of the nine pairs.
# 4. **$k_s$, and therefore Tables 13a/13b's first column, is unchecked**, the
#    critical and Lennard-Jones temperatures behind eqs. (4) and (5) are not
#    printed. V1 and V2 take that column as given; had Robeson mis-evaluated
#    eq. (4), they would both still pass.
# 5. **The reference selectivity is this page's choice.** The geometric mean of
#    the printed points is defensible and stated, but the front-factor *share*
#    is not a property of the data alone, and the page says so in three places
#    because it is the single easiest thing to quote out of context. The same
#    applies to Result 3's evaluation permeability: the ranking survives every
#    reference tried inside the range all nine pairs share, but the *gap*
#    between third and fourth does not, it is swept and printed above, with
#    both its smallest and its largest value, and neither should be quoted
#    without its reference. **And it applies to the second route's own
#    headline**: the residual between the two routes carries $\log\alpha_0$
#    explicitly, so it is swept over in-data references the same way, and both
#    the residual and the sign agreement are reported at their worst as well as
#    at the reference used. Neither number means anything unattached.
# 6. **Figure 12's slope–diameter correlation and Freeman's eq. (2) are out of
#    scope**, for a stated reason: the inputs are not printed in this paper.
#    The page reports no verdict on either.
# 7. **The envelope re-fit (V3) is not a competitor to the printed $(k,n)$** and
#    must not be quoted as a "corrected" bound. It answers a different question
#    (the tightest covering line), and it is pinned by exactly the extreme
#    points the paper distrusts. It is also fitted to points Robeson *selected*
#    as close to his 2008 line, so the second route's small **baseline**
#    residual is an optimistic figure and should not be read as accuracy. That
#    limits the baseline agreement and not the route's **detecting power**:
#    $k_{\text{env}}$ is untouched by Table 12, so the residual responds to a
#    vertical error in the printed present line exactly one for one, measured
#    under break row (DD) above, not argued.
# 8. **The excursion count is a count against *this* paper's own lines.** It
#    says nothing about whether those 19 polymers really beat the state of the
#    art; it says the hand-drawn line does not enclose every point the author
#    tabulated beside it, which he never claimed it would.
# 9. **"The paper does not state the crossing" is a searched claim, not an
#    assumption.** The full text was searched for cross / crossing /
#    intersect / "below the prior"; the only hits are the polymer name
#    *crosslinked*, reference titles containing *crosslinking*, and nothing
#    else. What the paper *does* print, and this page quotes, is the slope
#    change and the warning that the low-permeability end may have skewed it.
# 10. **Result 1's decomposition has no second arithmetic route, and cannot
#    have one.** Given Table 12 it is a rearrangement, so V2 and the data-side
#    shift both reproduce it to machine precision by construction; the page
#    measures that rather than claiming otherwise. The one genuine second route
#    is V3's envelope, which replaces the *present* line with an estimate from
#    Tables 1–11; at each pair's geometric-mean selectivity it agrees on the
#    sign of the shift for every comparable pair (8 of 9 at the least kind of
#    the four references swept), and on the magnitude to the tolerance printed
#    in V3 wherever the lever arm allows, but it does **not** reproduce the front-factor share at all, which
#    is blind spot 5 arriving from the other side. Its scope is the **present
#    line only**: nothing on this page can give the *prior* block a second
#    route, because both routes read it from the same restatement, and $\alpha_0$
#    is likewise common to both. And its corroboration of Result 3's ranking is
#    **weaker than the primary route's**: swept over the same evaluation points
#    the primary route is held to, the substituted route reproduces the ranking
#    at most of them and not at all of them, while the primary route holds at
#    every one. The counts are printed above; the asymmetry is the honest
#    statement and it is the reported one.
# 11. **The "inside its own tabulated range" half of Result 2 rests on a single
#    tabulated point**, named and quantified in Result 2 and removed by break
#    row (T1) above. The crossing itself does not (it is a property of
#    the printed $(k,n)$), but the claim that it lands where the paper has data
#    does.

# %% [markdown]
# ## What pymrm adds
#
# **Nothing.** This page calls no pymrm function, and the honest statement is
# that the Robeson upper bound is not a pymrm problem: there is no field to
# discretise, no operator to assemble, no Jacobian to build. The bootstrap
# installs pymrm only because every gallery page pins it. The `structures` field
# is empty and `pymrm_api` is empty because that is true, not because they were
# not filled in.
#
# What the *page* adds, over reading Table 12:
#
# - It separates the **nine** comparable pairs from the thirteen listed, and
#   shows that Robeson's own Conclusions do the same, so a reader who compares
#   thirteen is departing from the source, not following it.
# - It makes "primarily the front factor" **quantitative and
#   reference-explicit**. On the four references *Robeson's own Table 13*
#   supplies, the front-factor share of the same two shifts runs from **0.56 to
#   1.03**; across the nine pairs at this page's stated reference (the
#   geometric-mean selectivity of each pair's own points, which the paper does
#   not supply) it runs **0.59 to 1.39**. Both ranges are printed above; the
#   distinction between them is the point, and quoting either as "the paper's"
#   without saying which is exactly the error the page exists to warn against.
#   Robeson's own reference, $\alpha = 1$, is not outside the data either: it
#   lies inside the tabulated range of 5 of the nine pairs.
# - It shows that $\Delta\log k$, the quantity the abstract highlights, does
#   **not** reproduce the paper's own ranking of which bounds moved most, and
#   that the selectivity gain at fixed permeability does, exactly, at every
#   evaluation point inside the range all nine pairs share, while the *margin*
#   between third and fourth is itself reference-dependent and runs from 0.011
#   to 0.111 decades across those points. Both of the paper's statements
#   survive; they are statements about different things.
# - It finds the one pair whose bounds **cross inside its own data range**, and
#   says how thin that is: the crossing is a property of four printed numbers,
#   but "inside its own range" rests on a single tabulated point.
# - It measures what *"by eye"* costs: 19 of 117 points above the line. At
#   each pair's $\alpha_{gm}$ the tightest covering line through the same
#   points sits above the printed one for 7 pairs and below it for 4; that
#   count moves with the reference, while the slope does not: the covering line
#   is shallower than the printed one for 10 pairs and steeper for only 1.
# - It says which of its own checks is a **second route** and which is an
#   identity. V2 looks like independent corroboration and is arithmetically
#   incapable of failing unless V1 does; the real second route is the envelope
#   of *Validation V3*, the only computation here that does not read the present
#   bound out of Table 12. That route corroborates the size of the shift, at a
#   stated reference and by an amount that itself depends on that reference
#   (which is why the residual is swept and not quoted at one point), and
#   refuses to corroborate its front-factor split, which is the page's own
#   central claim arriving from a direction the page did not choose. It also
#   corroborates Result 3's ranking **less strongly than the primary route
#   does**, under the identical evaluation-point sweep, and the page reports
#   that gap rather than the flattering half of it.
# - It leaves every printed defect in place and quoted: the doubled `He/CO2`
#   in the Conclusions, the `Polypyrrole`/`Polypyrrolone` inconsistency between
#   Tables 2 and 5/6, the one Table 13a row outside even a generous printing
#   band, and the paper's own "the one data point above the present upper bound"
#   where its Table 1 has 2.
#
# For the transport model that produces a $(P,\alpha)$ pair in the first place,
# go to [`H1.7`](../H1.7-solution-diffusion/); that page *is* a pymrm page.

# %% [markdown]
# ## Reuse
#
# **Use the bound as a sanity ceiling, and quote it correctly.** Given a claimed
# $(P_i,\alpha_{ij})$ for a homogeneous polymer film, `excursion(P, alpha, k, n)`
# in `model.py` returns $r$; $r>1$ means the claim is above the 2008 bound. Take the
# $(k,n)$ from the `present` block of the Table 12 dataset, and note that only
# eleven pairs have one.
#
# **Four cautions, in order of how often they bite.**
#
# 1. **Nine, not thirteen.** If you are comparing 1991 with 2008, H2/O2 and
#    He/O2 have no 2008 bound and CO2/N2 and N2/CH4 have no 1991 one. Comparing
#    them is comparing against `NA`.
# 2. **Never quote a front-factor share without its reference selectivity.**
#    $\Delta\log k$ is the shift at $\alpha=1$, a selectivity with no special
#    standing in the data. On this page's own numbers the share of the He/H2 shift carried
#    by $k$ is 0.98 at the geometric-mean selectivity of its points and 1.02 and
#    1.03 at the two selectivities Robeson's own Tables 13a/13b use. For H2/CO2
#    it is **0.56 and 0.70 at the two references Table 13 supplies**, **1.39 at
#    this page's geometric-mean reference**, and undefined at $\alpha\approx 37$
#    where the two lines cross. Name the reference or the number means nothing:
#    the same pair, the same two bounds, and a share running from 0.56 to 1.39.
# 3. **Rank pairs by selectivity gain, not by $k$.** $|n|$ runs from 0.79
#    (He/CH4) to 5.8 (O2/N2) in Table 12, so equal moves in $\log k$ buy
#    wildly unequal separation. This is the page's sharpest practical result and
#    the easiest thing to get wrong.
# 4. **The bound applies to homogeneous polymer films only.** The paper excludes
#    heterogeneous, surface-modified, mixed-matrix, carbon molecular-sieve and
#    thermally-rearranged membranes explicitly, and shows TR polymers above the
#    CO2/CH4 bound as a deliberate comparison rather than as a refutation. A
#    Pd membrane ([`H1.1`](../H1.1-sieverts-permeation/)) or a zeolite membrane
#    ([`H1.9`](../H1.9-zeolite-membrane-maxwell-stefan-mixture/)) exceeding it is
#    not news.
#
# **If you extend this page**, the highest-value additions are, in order:
# (i) the 1991 paper itself, which would let the prior block be checked rather
# than restated and would make every shift on this page a validated number
# instead of a reproduced one; (ii) Breck's kinetic diameters, which would open
# Fig. 12's $-1/n$ correlation; (iii) solubility constants, which would open
# Freeman's eq. (2). None of the three is reachable from this file, and the page
# claims nothing about any of them.
#
# **Datasets** are reusable as they stand: `robeson-2008-table12-upper-bounds`
# is the canonical $(k,n)$ table, and `robeson-2008-near-bound-points` is 117
# literature points with the compiler's citation on every row. Both carry the
# caveats above in their sidecars; read them before loading, per the gallery's
# cross-page rule.

# %%
# Every number in the prose above, in README.md and in meta.yaml, checked against its claim.
for lit, name in [
        ("382", "cells_total"), ("44", "cells_table12"), ("234", "cells_tables1_11"),   # Background
        ("104", "cells_table13"), ("26", "t13_rows"),
        ("13", "pairs_in_table12"), ("44", "references_cited"),                         # The data
        ("more than 7", "n_abs_ratio"), ("0.79", "n_abs_min"), ("5.8", "n_abs_max"),   # intro, model, Reuse
        ("37", "h2co2_crossing_alpha"), ("37.47", "h2co2_crossing_alpha"),
        ("19", "points_above_present_bound"), ("117", "points_total"),
        ("3.5", "excursion_max"), ("3.49", "excursion_max"),
        ("5", "pairs_alpha1_inside_range"),        # model, assumptions
        ("3.742", "t13_mw_h2n2_printed"),
        ("0.56", "v2_kshare_min"), ("1.03", "v2_kshare_max"),                            # What pymrm adds
        ("0.5647", "v2_kshare_min"), ("1.0334", "v2_kshare_max"),
        ("0.59", "front_factor_share_min"), ("1.39", "front_factor_share_max"),
        ("0.98", "front_factor_share_median"),
        ("0.011", "dlogA_gap_min_over_eval_points"), ("0.111", "dlogA_gap_third_to_fourth"),
        ("7", "envelope_offsets_positive"), ("4", "envelope_offsets_negative"),
        ("2", "o2n2_points_above_present_bound"),
        ("0.98", "heh2_share_gm"), ("1.02", "v2_kshare_heh2_13a"),                       # Reuse
        ("1.03", "v2_kshare_heh2_13b"), ("0.56", "v2_kshare_h2co2_13a"),
        ("0.70", "v2_kshare_h2co2_13b"), ("1.39", "h2co2_share_gm"),
        ("6", "pairs_near_parallel"), ("9", "pairs_comparable"),                         # README, meta.yaml
        ("0.34", "h2n2_spread_over_shift"), ("0.93", "h2ch4_spread_over_shift"),
        ("2.05", "h2co2_spread_over_shift"), ("3.0 %", "heh2_spread_over_shift"),
        ("1.83", "heh2_shift_decades"), ("3.1 %", "heh2_slope_term_max_frac"), ("3", "dlogA_top3_hits"), ("2", "dlogk_top3_hits"),
        ("2", "dlogP_top3_hits"), ("-30.4", "h2co2_present_over_prior_at_alpha_max_pct"),
        ("30.4", "h2co2_present_below_prior_at_alpha_max_pct"), ("30", "h2co2_present_below_prior_at_alpha_max_pct"),
        ("100.9", "h2co2_alpha_max"), ("16.2 %", "frac_above_present_bound"),
        ("0.478", "excursion_geomean"), ("0.017 %", "t13_identity_median_abs_rel"),
        ("25", "t13_rows_inside_loose_band"), ("24", "t13_rows_inside_tight_band"),
        ("-0.612 %", "t13a_o2n2_rel"), ("0.279 %", "t13a_o2n2_band_loose"),
        ("1.4714", "t13a_o2n2_alpha_needed"), ("1.473", "t13a_o2n2_alpha_printed"),
        ("-0.027 %", "t13b_o2n2_rel"), ("1.5e-4", "t13_mw_column_max_abs_rel"),
        ("3.9e-3", "t13_mw_column_max_abs_rel_iupac"), ("1.6e-3", "t13_product_column_max_abs_rel"),
        ("4.9e-4", "v2_shift_crosscheck_max_abs_decades"),
        ("below 1e-12", "v2_is_v1_difference_max_abs"),
        ("below 1e-12", "dataside_shift_is_decomposition_max_abs"),
        ("0.086", "envelope_offset_max_abs_wide_span"),
        ("0.543", "envelope_offset_max_abs_wide_span_max_over_alpha0"),
        ("9", "envelope_sign_agree_at_alpha_gm"), ("8", "envelope_sign_agree_min_over_alpha0"),
        ("0.466", "envelope_offset_max_abs_decades"), ("0.62", "envelope_span_min_decades"),
        ("0.278", "envelope_dn_max_abs_wide_span"), ("4", "lever_arm_set_size"),
        ("3", "envelope_implied_dlogA_top3_hits"),
        ("2", "envelope_implied_dlogA_top3_hits_min_over_eval_points"),
        ("5", "envelope_eval_points_with_3_hits"),
        ("0.43", "envelope_implied_share_median"), ("-0.03", "envelope_implied_share_min"),
        ("1.55", "envelope_implied_share_max"), ("-0.0092", "envelope_hen2_dlogk"),
        ("0.299", "envelope_hen2_shift"), ("1.11", "envelope_hen2_cross_alpha"),
        ("1", "pairs_crossing_inside_range"), ("1", "h2co2_points_above_crossing"),
        ("5", "pairs_classified_modest"), ("1", "pairs_classified_minor"),            # Result 3
        ("3", "pairs_classified_significant"),
        ("10", "envelope_dn_positive"), ("1", "envelope_dn_negative"),                # What pymrm adds
        ("7", "pairs_span_over_1p5_decades"),                                         # meta.yaml
        ("0.92 %", "t13a_n2ch4_band_tight"), ("8.43 %", "t13a_n2ch4_band_tight_float_read"),  # README
        ("4", "v1_heh2_rows_inside_tight_band")]:
    claims.says(lit, name)

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
