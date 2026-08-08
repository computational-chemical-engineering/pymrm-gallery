#!/usr/bin/env python3
"""Generate index.ipynb for page H1.8 (Robeson upper bound, 2008 revisit).

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "The Robeson upper bound, revisited: where the 2008 line moved, and what \"primarily the front factor\" can mean"
description: "Robeson's 2008 Table 12 prints k and n for thirteen gas pairs, but only NINE support a prior-versus-present comparison. On those nine the shift is decomposed into a front-factor part and a slope part at a stated selectivity - and the decomposition is shown to be reference-dependent, because log-log lines pivot. The paper's own verbal ranking of which bounds moved most is reproduced exactly by the SELECTIVITY gain at fixed permeability and by neither delta-log-k nor the permeability gain at fixed selectivity. Tables 13a/13b reproduce all 26 rows from Table 12 - median deviation 0.017 %, 25 rows within 0.11 %, one at 0.61 % and outside even a generous printing band. The decomposition is definitional given Table 12, and the page says so: its only genuine second route replaces the present bound with a covering line re-estimated from the paper's own 117 near-bound points, which agrees on the sign of every shift and on the magnitude wherever those points span enough selectivity to determine a line - AT A STATED REFERENCE SELECTIVITY, since that residual carries log(alpha_0) too and is swept across references rather than quoted at one - and does NOT reproduce the front-factor share, which is the same warning arriving from the other side. 19 of the 117 points the paper tables as \"close to the present upper bound\" lie above it."
categories: [sec:H, tier:T0, data:tier2, phase:gas, kind:correlation]
date: 2026-08-08
---

# The Robeson upper bound, revisited

**Catalog ID:** `H1.8` · **Structures:** none (see below) · **Tier:** T0

Every polymer membrane paper that plots selectivity against permeability draws
one diagonal line on it. That line is

$$
P_i \;=\; k\,\alpha_{ij}^{\,n},
\qquad \alpha_{ij}=\frac{P_i}{P_j},\qquad n<0,
$$

Robeson's empirical **upper bound**: for a given selectivity there is a
permeability no homogeneous polymer film has beaten. In 1991 he drew it for
nine gas pairs. In 2008 he redrew it against seventeen more years of
literature and printed both sets of $(k,n)$ side by side, together with the
2008 data points that set the new line, together with a Knudsen-transition
analysis computed from those same $(k,n)$.

The paper's central testable claim is a **decomposition**: the bound moved,
and *"the shift observed is primarily due to a change in the front factor,
$k$, whereas the slope of the resultant upper bound relationship remains
similar to the prior data correlations."* That is arithmetic on one printed
table — and this page does it, on the **nine** gas pairs that actually support
it (not the thirteen the table lists), at a **stated** reference selectivity,
because a shift decomposed into "front factor" and "slope" parts is
**meaningless until the reference is named**: two log–log lines with different
slopes pivot, so moving the reference moves the split.

Four results, all from printed numbers:

1. **Nine, not thirteen.** H2/O2 and He/O2 have a prior bound and no present
   one; CO2/N2 and N2/CH4 have a present bound and no prior one. Any statement
   about "the shift" lives on nine pairs. Robeson's own Conclusions name
   exactly those nine.
2. **The claim holds where the paper says it holds, and the number to quote is
   not $\Delta\log k$.** Because $n$ varies by a factor of seven across the
   pairs, the same $\Delta\log k$ buys wildly different amounts of separation.
   Ranked by the **selectivity gain at fixed permeability**, the paper's three
   *"significant"* pairs — He/CH4, He/CO2, He/H2 — come out as the top three
   with a clean gap. Ranked by $\Delta\log k$, or by the permeability gain at
   fixed selectivity, they do not: CO2/CH4, which the paper calls *modest*,
   interleaves with them.
3. **The two H2/CO2 bounds cross inside the paper's own data range.** For
   $\alpha \gtrsim 37$ the 2008 "upper bound" lies *below* the 1994 one — and
   Table 10 has **one** point out there, which is the whole of the "inside the
   range" half of this claim and is stated as such wherever it appears. The
   paper notices the slope change and warns its low-permeability end may be
   skewed; it does not state the crossing.
4. **19 of the 117 points the paper tables as "close to the present empirical
   upper bound" lie above it**, by up to a factor 3.5 in permeability. That is
   not an error: the line is drawn *"by eye"*, as the paper says. It is what
   "by eye" amounts to, measured.

The page is a **correlation audit**, not a reactor model. There is no PDE, no
grid and no pymrm operator in it, and the *What pymrm adds* section says so
plainly."""))

# ----------------------------------------------------------------- background
cells.append(md(r"""## Background

**The source, precisely.** Lloyd M. Robeson, *"The upper bound revisited"*,
Journal of Membrane Science **320**(1–2) 390–400 (2008),
doi:`10.1016/j.memsci.2008.04.030`, Lehigh University, 1801 Mill Creek Road,
Macungie, PA 18062, United States. Identified from the article's own running
header (*"Journal of Membrane Science 320 (2008) 390–400"*), display title,
by-line and affiliation block, all read on a 300 dpi render of p. 390; the
PDF's own Title metadata carries the DOI.

The file is **born-digital** (Elsevier / Acrobat Distiller 7.0), so the text
layer is a faithful extraction rather than OCR. It was still not trusted:
every numeric cell used on this page was read a second time from cropped
300 dpi renders at digit scale, and the two readings compared. They agree on
all 382 **numeric cells** — 44 in Table 12 (its four `NA` rows aside), 234 in
Tables 1–11 and 104 in Tables 13a/13b, counted in the loading cell below rather
than asserted. (Cells, not numerals: the digit count is several times larger,
and an earlier draft of this page said "numerals", which was wrong.) The reason for the caution is specific to this paper — its tables
carry **comma thousands separators inside cells** (`1,073,700`, `30,967,000`,
`5,369,140`), which is exactly the class of glyph a naive read drops or
mis-parses. The catalogue's earlier route to this material, an Elsevier
API text dump of the *1991* paper, returned the 1991 Table 2 of kinetic
diameters as `2 6 2 89 3 3 3 46 3 64 3 8` for 2.6, 2.89, 3.3, 3.46, 3.64,
3.8 — a six-value row read as eleven values, with nothing visibly wrong. This
PDF retires that route.

**The 1991 and 1994 papers are the ORIGIN of the "prior" numbers and were not
consulted.** Robeson (1991), *J. Membr. Sci.* **62**, 165, and Robeson,
Burgoyne, Langsam, Savoca & Tien (1994), *Polymer* **35**, 4970, are not on
disk. Every prior-bound value on this page was read from the 2008 paper, which
restates them under its own heading *"Prior upper bound data (1,2)"*. That is
the legitimate reprint route (`AGENTS.md`), and it has a limit this page
respects throughout: **nothing here can establish that the 2008 restatement of
the 1991 bounds is faithful to the 1991 paper.** Where a prior value carries a
result, the page says it is read from the 2008 restatement.

**What the paper does.** It re-reviews the literature since 1991 — *"the
number of papers where the data were obtained approached 300"* — redraws the
bound for eleven gas pairs (adding CO2/N2 and N2/CH4, which had no prior
correlation), tabulates the points that set each new line (Tables 1–11),
tabulates old and new $(k,n)$ together (Table 12), and then uses the new
bounds to estimate where solution–diffusion transport must give way to Knudsen
diffusion (Tables 13a and 13b). Its physical headline is that the big moves
are almost all in **He-based pairs** and are driven by a family of
perfluorinated commercial polymers — Nafion, Hyflon AD, Viton, Cytop, Teflon
AF — that were barely in the 1991 dataset.

**How the line is drawn matters, and the paper is honest about it.** From
p. 391: *"It is noted that the 'upper bound' line is determined empirically
('by eye') as in the original reference [1] with sufficient data to establish
a realistic bound over several decades of permeability."* So $(k,n)$ are **not
a fit** to Tables 1–11 in any least-squares sense. Nothing on this page treats
them as one, and result 4 above is the quantitative consequence.

**Where it sits in the gallery.** [`H1.7`](../H1.7-solution-diffusion/) is the
transport model underneath — permeability as solubility times diffusivity,
same units, same barrer. This page supplies the empirical *ceiling* that model
lives under; it does not derive it. [`H1.1`](../H1.1-sieverts-permeation/) and
[`H1.9`](../H1.9-zeolite-membrane-maxwell-stefan-mixture/) are the metal and
zeolite membranes to which the bound explicitly does **not** apply (p. 391:
*"Heterogeneous membranes, surface modified membranes and molecular sieve
membranes are not considered in the same class of polymeric materials"*).
Methodologically the nearest sibling is
[`B2.1`](../B2.1-voorhies-coking-law/): the origin of an empirical power law,
tested only on what the source itself prints.

**What this page deliberately does not do.** Twelve figures carry the actual
data clouds; none is digitised, and none needs to be, because Tables 1–13
print everything the arguments here use. Figure 12a/b — the correlation of
$-1/n$ against gas-diameter difference — **cannot** be checked from this
paper: it needs Lennard-Jones kinetic diameters, which the 2008 paper cites to
Breck [65] and never prints. That is a checked claim, not an assumption: the
full text layer was searched for the string `diameter` — **17 occurrences**
(8 singular, 9 plural, counted as occurrences in the `pdftotext -layout`
extraction of all eleven pages, not as matching lines), every one of them
prose, a caption or eq. (1) — and the only numeric diameters printed anywhere
in the paper are PTMSP's pore size, *"in the range of 0.9–1.2 nm"*. Eq. (2), Freeman's theoretical
prediction of $k$, is likewise uncheckable here — it needs solubility
constants the paper does not tabulate. Both are out of scope and stay out."""))

# ----------------------------------------------------------- published model
cells.append(md(r"""## The published model

Equation numbers are the paper's.

**The bound itself** (abstract, and p. 396 in the form Table 12 heads):

$$P_i = k\,\alpha_{ij}^{\,n}$$

with $P_i$ the permeability of the **fast** gas in barrers
($1\ \text{barrer} = 10^{-10}\ \mathrm{cm^3(STP)\,cm\,cm^{-2}\,s^{-1}\,cmHg^{-1}}$),
$\alpha_{ij}=P_i/P_j$ the separation factor, $k$ the *"front factor"*, and $n$
the slope of the log–log plot. Since $n<0$, in the plane the field actually
plots — $\log\alpha$ against $\log P$ — the bound is the falling diagonal
$\log\alpha = (\log P - \log k)/n$, and *"below this line ... virtually all
the experimental data points exist"*.

**The shift, decomposed.** Write the prior bound $(k_0,n_0)$ and the present
one $(k_1,n_1)$. At a selectivity $\alpha$ the bound permeability moves by

$$
\underbrace{\Delta \log_{10} P(\alpha)}_{\text{vertical shift}}
=\underbrace{\Delta\log_{10}k}_{\text{front-factor part}}
+\underbrace{\Delta n \,\log_{10}\alpha}_{\text{slope part}} .
$$

Three things follow immediately, and they organise the whole page.

*First*, **the split depends on where you stand.** The front-factor part is
the shift evaluated at $\alpha = 1$ — no separation at all, outside every
dataset in the paper. Reparameterise the same two lines about any other
reference $\alpha_0$, writing $P = k'(\alpha/\alpha_0)^n$, and the
"front-factor part" becomes $\Delta\log k + \Delta n \log\alpha_0$. So
"primarily the front factor" is not a property of the pair of lines; it is a
property of the pair of lines *and a chosen reference*. This page always
states the reference.

*Second*, the statement that **is** reference-free is that the two lines are
nearly **parallel**: $|\Delta n|\log_{10}(\alpha_{\max}/\alpha_{\min})$, the
amount the vertical shift varies across the selectivity range the data
actually span, compared with the shift itself. That is the honest reading of
Robeson's claim, and it is what gets measured below.

*Third*, two non-parallel lines **cross**, at
$\log_{10}\alpha^\star = -\Delta\log_{10}k/\Delta n$. If $\alpha^\star$ falls
inside the data range, the "new" bound is *lower* than the old one over part
of that range.

**A second, orthogonal way to measure the same move.** A membrane engineer
does not ask how much permeability the bound gained at fixed selectivity; they
ask how much **selectivity** it gained at fixed permeability:

$$
\Delta \log_{10}\alpha(P)=\frac{\log_{10}P-\log_{10}k_1}{n_1}
-\frac{\log_{10}P-\log_{10}k_0}{n_0}.
$$

These two measures are *not* monotonically related across pairs, because they
differ by a factor $\sim 1/|n|$ and $|n|$ runs from 0.79 to 5.8 in this table.
Which one reproduces Robeson's own verbal ranking is a result, not a
convention.

**The Knudsen transition** (eqs. 3–5, Tables 13a/13b). Knudsen diffusion
selectivity is $D_i/D_j=(M_j/M_i)^{1/2}$ (eq. 3), and the solubility
selectivity $k_s=S_i/S_j$ is taken invariant with permeability, so the
separation factor at which solution–diffusion selectivity has fallen to the
Knudsen value is

$$\alpha_{ij}=k_s\,(M_j/M_i)^{1/2},$$

and the permeability there is read off the bound as $P_i = k\,\alpha_{ij}^{n}$.
$k_s$ itself comes from one of two correlations Freeman gives — eq. (4) from
critical temperatures, $\ln(S_i/S_j)=0.016(T_{ci}-T_{cj})$, giving Table 13a;
eq. (5) from Lennard-Jones temperatures,
$\ln(S_i/S_j)=0.023[(\varepsilon_i/k)-(\varepsilon_j/k)]$, giving Table 13b.
Neither temperature set is printed, so **$k_s$ is an input here, not a
checkable quantity**. The last column, however, is entirely checkable against
Table 12, and reproducing all 26 of its rows is this page's strongest test of
its own transcription (V1 below)."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

| symbol | meaning | source | status here |
|---|---|---|---|
| $k$ | front factor, barrers | Table 12, both blocks | transcribed; prior block via the 2008 restatement |
| $n$ | log–log slope, dimensionless, $<0$ | Table 12, both blocks | transcribed |
| $P_i,\alpha_{ij}$ | the 117 near-bound points | Tables 1–11 | transcribed; **compiled by Robeson from ~300 papers, none consulted here** |
| $k_s$ | solubility selectivity | Tables 13a/13b, col. 1 | **input** — eqs. (4)/(5) need temperatures the paper does not print |
| $(M_j/M_i)^{1/2}$ | Knudsen selectivity | Tables 13a/13b, col. 2 | checkable, and checked, from integer molecular weights |
| $\alpha_0$ | reference selectivity for the decomposition | **not printed anywhere** | *chosen here*, and always stated |

Assumptions, each one a decision that could have gone otherwise:

1. **The reference selectivity $\alpha_0$ is the geometric mean of the
   printed near-bound $\alpha$ for that pair.** It is data-driven, it is
   deterministic, and it comes from the paper's own tables — but it is a
   choice, it is not Robeson's, and the page reports the whole $\alpha$ range
   beside it so the choice can be audited. Robeson's implicit reference is
   $\alpha_0=1$, which no pair's data reach.
2. **Molecular weights are integers** (He 4, H2 2, O2 32, N2 28, CO2 44,
   CH4 16). This is not a convenience: it is *identified from the paper*.
   IUPAC standard atomic weights give $(M_{N_2}/M_{H_2})^{1/2}=3.728$; the
   table prints **3.742**, which is $\sqrt{28/2}$. All eleven entries follow.
3. **The printed rounding band is taken literally** — half a unit in the last
   printed digit position, so `1,396,000` carries $\pm 0.5$ and `−5.666`
   carries $\pm0.0005$. Where a conclusion depends on the band, the page also
   reports the **generous** band that treats trailing zeros as placeholders
   (four significant figures, so `1,396,000` becomes $\pm 500$), and states
   which reading the conclusion survives.
4. **A point lies above the bound when $r=P/(k\alpha^{n})>1$.** With $n<0$
   this is the correct direction and it is worth deriving once: a point beats
   the bound if its selectivity exceeds the bound's at its permeability,
   $\log\alpha>(\log P-\log k)/n$; multiplying by the *negative* $n$ flips the
   inequality to $\log k+n\log\alpha<\log P$, i.e. $r>1$. Getting this
   backwards turns 19 excursions into 98, which is a break-table row below.
5. **Nothing is refit as a substitute for the printed bound.** The line is
   drawn by eye; a least-squares line through Tables 1–11 is a *different
   object* and is labelled as one wherever it appears."""))

# -------------------------------------------------------------------- setup
cells.append(code('''# Colab environment cell
try:
    import pymrm  # noqa: F401  (not used for physics on this page - see "What pymrm adds")
except ImportError:  # pragma: no cover
    %pip install -q pymrm

import sys, subprocess
from pathlib import Path
for parent in [Path.cwd(), *Path.cwd().parents]:
    if (parent / "shared" / "gallery_utils.py").is_file():
        sys.path.insert(0, str(parent / "shared")); break
else:  # pragma: no cover - Colab
    subprocess.run(["pip", "install", "-q", "pyyaml"], check=False)
    url = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
           "pymrm-gallery/main/shared/gallery_utils.py")
    import urllib.request; urllib.request.urlretrieve(url, "gallery_utils.py")
    sys.path.insert(0, ".")

import textwrap
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.optimize import brentq
from gallery_utils import load_data, load_meta, cite_data, report_agreement

np.set_printoptions(precision=6, suppress=False)
pd.set_option("display.width", 130)
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
PAGE = "H1.8-robeson-upper-bound"
print("numpy", np.__version__, "| pandas", pd.__version__)'''))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

Four datasets, all transcription of printed tables and printed sentences.
**Tier 2** throughout: the paper's own tabulated values, no figure digitised,
no number inferred.

- `robeson-2008-table12-upper-bounds` — Table 12, 26 rows (13 pairs × prior
  and present), with the `NA` structure preserved.
- `robeson-2008-near-bound-points` — Tables 1–11, 117 rows.
- `robeson-2008-table13-knudsen-transition` — Tables 13a/13b, 26 rows.
- `robeson-2008-printed-claims` — the prose claims tested here, verbatim.

No dataset from another gallery page is loaded, so the cross-page
reconciliation rule has nothing to bite on here.

The provenance depth of the second dataset is worth stating plainly, because
it bounds every result drawn from it: **these are not Robeson's measurements.**
Each row cites the primary paper that reported it, and none of those ~65
primary papers was consulted. This page can show that a tabulated point sits
above the line Robeson drew; it cannot show that the point is right."""))

cells.append(code('''# `dtype=str` + `keep_default_na=False` so the as-printed columns keep their trailing
# zeros and their literal "NA" - both are load-bearing below.
T12 = load_data("robeson-2008-table12-upper-bounds.csv", page=PAGE,
                dtype=str, keep_default_na=False)
for c in ("k_barrer", "n"):
    T12[c] = pd.to_numeric(T12[c].replace("", np.nan))
PTS = load_data("robeson-2008-near-bound-points.csv", page=PAGE)
T13 = load_data("robeson-2008-table13-knudsen-transition.csv", page=PAGE)
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

prior = T12[T12.dataset == "prior"].set_index("gas_pair")
present = T12[T12.dataset == "present"].set_index("gas_pair")
PAIRS_ALL = list(prior.index)
PAIRS = [g for g in PAIRS_ALL if np.isfinite(prior.k_barrer[g]) and np.isfinite(present.k_barrer[g])]
N_PAIRS_TABLE, N_PAIRS_COMPARABLE = len(PAIRS_ALL), len(PAIRS)

print(f"\\nTable 12 lists {N_PAIRS_TABLE} gas pairs; {N_PAIRS_COMPARABLE} carry BOTH a prior and a present bound.")
print("  present bound NA (prior only):",
      ", ".join(g for g in PAIRS_ALL if not np.isfinite(present.k_barrer[g])))
print("  prior bound NA (present only):",
      ", ".join(g for g in PAIRS_ALL if not np.isfinite(prior.k_barrer[g])))
print("  comparable:", ", ".join(PAIRS))
print("\\nTables 1-11 exist for exactly the 11 pairs with a PRESENT bound:")
print("  ", ", ".join(sorted(PTS.gas_pair.unique())))
print("  pairs with a present bound but no table:",
      sorted(set(present.index[np.isfinite(present.k_barrer)]) - set(PTS.gas_pair.unique())) or "none")'''))

cells.append(code('''# The paper's own three-way classification of how far each bound moved, joined onto Table 12.
CLASS = (CLM[CLM.shift_class.notna() & (CLM.gas_pair.notna())]
         .drop_duplicates("gas_pair").set_index("gas_pair").shift_class.to_dict())
print("Robeson's Conclusions classify these pairs, and only these:")
for g in PAIRS:
    print(f"  {g:9s} {CLASS[g]}")
missing = sorted(set(PAIRS) - set(CLASS))
extra = sorted(set(CLASS) - set(PAIRS))
print(f"\\nclassified but not comparable: {extra or 'none'}"
      f"  |  comparable but unclassified: {missing or 'none'}")
assert not missing and not extra, "the Conclusions' pair list must be exactly the nine comparable pairs"
print("\\nSo the paper itself confines its shift discussion to the nine pairs the NA structure allows.")
print("A page that compared thirteen would be comparing four pairs against nothing.")'''))

cells.append(code('''# H2/O2 and He/O2: prior bound printed, present bound NA, no section, no table, no figure.
# A negative claim, so here is the search. The text layer is born-digital; extraction inserts a
# space before the slash for H2 but not for He, so both spellings are counted.
print("Where do H2/O2 and He/O2 appear in the paper?")
print("  Searched: the complete `pdftotext -layout` extraction of all 11 pages, for both the")
print("  hyphen-free spellings and the space-before-slash form the extractor produces.")
print("  Result: 4 occurrences in total - two in Table 12's prior block, two in its present block")
print("  (as 'NA'). No section heading, no data table, no figure, and no sentence anywhere in the")
print("  paper mentions either pair. The paper does not say why they were dropped.")
print("  Consistency check available here: Tables 1-11 cover 11 pairs and neither of these two.")
assert set(PTS.gas_pair.unique()) == set(present.index[np.isfinite(present.k_barrer)]), \\
    "every present bound must have a data table and vice versa"
print("  -> verified: {} data tables, {} present bounds, same set.".format(
    PTS.gas_pair.nunique(), int(np.isfinite(present.k_barrer).sum())))'''))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

**There is none, and that is the honest answer.**

The Robeson upper bound is a two-parameter algebraic relation between two
scalars. It has no spatial coordinate, no time coordinate, no conservation
law, no boundary condition and no unknown field. Nothing in pymrm's operator
library — `construct_grad`, `construct_div`, `construct_convflux_upwind`,
`NumJac`, `newton` — has anything to act on. Building a grid here would be
decoration, and the gallery's structure taxonomy (`S1`–`S13` in
[`docs/taxonomy.md`](../../docs/taxonomy.md)) has no code that fits: the page's
`structures` field is **empty**, deliberately, and its `pymrm_api` list is
empty too.

What the page uses is numpy for the algebra, one `scipy.optimize.brentq` call
per threshold (thresholds are root-found, never swept), and one small
constrained least-squares solve for the envelope estimator in *Validation*.

Where the bound *does* meet pymrm is one layer up, and this page does not
cross that layer: $(P_i, \alpha_{ij})$ is the closure a module-scale membrane
model consumes. [`H1.7`](../H1.7-solution-diffusion/) builds the transport side
of that closure from solubility and diffusivity in the same units; the bound
tells you which $(P,\alpha)$ pairs that model is entitled to be handed. The
two pages are complementary and neither computes the other.

The rest of this page is therefore arithmetic on printed tables — which is
exactly what auditing an empirical correlation is."""))

cells.append(code('''# ---------------------------------------------------------------------------
# The only "model" on this page, and its inverse. Both are one-liners; both are
# written once here and used everywhere below, so a sign error cannot hide in
# one branch.
# ---------------------------------------------------------------------------
def P_bound(alpha, k, n):
    """Upper-bound permeability (barrers) at separation factor alpha."""
    return k * np.asarray(alpha, float) ** n

def alpha_bound(P, k, n):
    """Upper-bound separation factor at permeability P (barrers).  Inverse of P_bound."""
    return (np.asarray(P, float) / k) ** (1.0 / n)

def excursion(P, alpha, k, n):
    """r = P / (k alpha^n).  r > 1 means the point lies ABOVE the bound (see Assumption 4)."""
    return np.asarray(P, float) / P_bound(alpha, k, n)

# round-trip identity: the inverse must invert. Structural, but it is the one thing
# that would let a reciprocal-of-n slip through unnoticed.
_a = np.array([0.2, 1.0, 7.5, 900.0])
for g in PAIRS:
    k, n = present.k_barrer[g], present.n[g]
    assert np.allclose(alpha_bound(P_bound(_a, k, n), k, n), _a, rtol=1e-12)
ROUNDTRIP_MAX = max(float(np.max(np.abs(alpha_bound(P_bound(_a, present.k_barrer[g], present.n[g]),
                                                    present.k_barrer[g], present.n[g]) / _a - 1)))
                    for g in PAIRS)
print(f"P_bound / alpha_bound round trip, worst over the {len(PAIRS)} pairs: {ROUNDTRIP_MAX:.3e}")

# a point exactly ON the bound has r = 1, by construction; a point at twice the bound
# permeability has r = 2. Sanity, and the sign convention made explicit.
kk, nn = present.k_barrer["He/H2"], present.n["He/H2"]
print(f"sign check on He/H2 (k={kk:,.0f}, n={nn}):")
for fac in (0.5, 1.0, 2.0):
    print(f"  a point at {fac:>3}x the bound permeability at alpha=2.5 has r = "
          f"{float(excursion(fac * P_bound(2.5, kk, nn), 2.5, kk, nn)):.3f}")''' ))

# ------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1 — The shift, decomposed at a stated selectivity

For each of the nine comparable pairs: $\Delta\log_{10}k$, $\Delta n$, the
range of selectivity its own near-bound points span, the vertical shift at the
geometric-mean selectivity of those points, and the **front-factor share**
$\Delta\log k/\Delta\log P(\alpha_0)$ — the fraction of the shift the front
factor accounts for *at that reference*.

The column that carries the reference-free version of Robeson's claim is
**`spread/shift`**: how much the vertical shift varies across the pair's own
selectivity range, relative to the shift itself. Small means the two lines are
parallel and "primarily the front factor" is safe wherever you stand; large
means the split is an artefact of where you stood."""))

cells.append(code('''rows = []
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

print("Shift of the upper bound, nine comparable pairs.  All logs base 10; 'decades' = log10 units.")
print(f"{'pair':9s} {'paper':12s} {'dlog k':>7s} {'dn':>8s} {'alpha range':>18s} "
      f"{'a_gm':>7s} {'shift@a_gm':>10s} {'k-share':>8s} {'spread/shift':>12s}")
for g, r in SH.iterrows():
    print(f"{g:9s} {r.verdict:12s} {r.dlogk:7.4f} {r.dn:+8.4f} "
          f"{r.a_lo:8.3g} - {r.a_hi:<7.4g} {r.a_gm:7.3g} {r.shift_gm:10.4f} "
          f"{r.share:8.3f} {r.spread_over_shift:12.2f}")

SHARE_MIN, SHARE_MAX = float(SH.share.min()), float(SH.share.max())
SHARE_MED = float(SH.share.median())
PARALLEL_N = int((SH.spread_over_shift < 0.20).sum())
print(f"\\nFront-factor share at the geometric-mean selectivity: median {SHARE_MED:.2f}, "
      f"range {SHARE_MIN:.2f} to {SHARE_MAX:.2f}.")
print(f"Two lines are near-parallel (spread/shift < 0.20) for {PARALLEL_N} of the "
      f"{len(SH)} pairs. For those, 'primarily the front")
print("factor' is a statement about the lines, not about a chosen reference. For the rest it is not:")
for g, r in SH[SH.spread_over_shift >= 0.20].iterrows():
    print(f"  {g:9s} the vertical shift runs {r.shift_lo:+.3f} -> {r.shift_hi:+.3f} decades across "
          f"alpha = {r.a_lo:.3g} to {r.a_hi:.4g}")
print("\\nAnd note what the reference does to the headline pair. He/H2's front factor moves")
print(f"{prior.k_barrer['He/H2']:,.0f} -> {present.k_barrer['He/H2']:,.0f} barrers, a factor "
      f"{present.k_barrer['He/H2']/prior.k_barrer['He/H2']:.1f}: dlog k = {SH.dlogk['He/H2']:.4f}.")
print(f"Its slope moves {prior.n['He/H2']} -> {present.n['He/H2']}, i.e. {100*abs(SH.dn['He/H2']/prior.n['He/H2']):.1f} % of itself.")
print(f"Over the selectivity range its OWN data span ({SH.a_lo['He/H2']:.2f} to {SH.a_hi['He/H2']:.2f}), the slope")
print(f"term contributes at most {100*SH.spread['He/H2']/abs(SH.shift_gm['He/H2']):.1f} % of a "
      f"{SH.shift_gm['He/H2']:.2f}-decade shift. There the claim is airtight.")'''))

cells.append(md(r"""### 2 — Where the two bounds cross

$\Delta\log P(\alpha)=0$ at
$\log_{10}\alpha^\star=-\Delta\log_{10}k/\Delta n$. Computed two ways — the
closed form, and `brentq` root-finding the shift function itself — and then
asked whether $\alpha^\star$ falls inside the range of the paper's own
tabulated points for that pair.

Two separate claims live here and they are of very different strengths. That
the two H2/CO2 lines cross at $\alpha^\star\approx 37$ is a property of four
printed numbers in Table 12. That the crossing happens **inside the range
Table 10 spans** is a property of a *single tabulated point* — the cell below
names it, prints the gap it sits across, and says what is left without it. A
break row further down mis-reads that one decimal point and takes the result to
zero."""))

cells.append(code('''print(f"{'pair':9s} {'alpha* (closed form)':>21s} {'alpha* (brentq)':>16s} "
      f"{'|diff|':>9s} {'inside its own alpha range?':>29s}")
cross_rows = []
for g in PAIRS:
    k0, n0 = prior.k_barrer[g], prior.n[g]
    k1, n1 = present.k_barrer[g], present.n[g]
    dlk, dn = np.log10(k1) - np.log10(k0), n1 - n0
    a_star_cf = 10.0 ** (-dlk / dn)
    f = lambda la: dlk + dn * la                        # noqa: E731 - shift in log-alpha
    lo, hi = np.log10(SH.a_lo[g]), np.log10(SH.a_hi[g])
    inside = f(lo) * f(hi) < 0
    a_star_rf = 10.0 ** brentq(f, lo, hi, xtol=1e-14, rtol=1e-15) if inside else np.nan
    diff = abs(a_star_rf / a_star_cf - 1) if inside else np.nan
    cross_rows.append(dict(pair=g, a_star=a_star_cf, a_star_rf=a_star_rf, inside=inside))
    print(f"{g:9s} {a_star_cf:21.4g} {a_star_rf if inside else float('nan'):16.6f} "
          f"{diff if inside else float('nan'):9.1e} "
          f"{('YES - inside ' + f'{SH.a_lo[g]:.3g}-{SH.a_hi[g]:.4g}') if inside else 'no':>29s}")
CROSS = pd.DataFrame(cross_rows).set_index("pair")
N_CROSSING_INSIDE = int(CROSS.inside.sum())
G_CROSS = CROSS.index[CROSS.inside][0]
ALPHA_STAR = float(CROSS.a_star[G_CROSS]); ALPHA_STAR_RF = float(CROSS.a_star_rf[G_CROSS])
ALPHA_STAR_AGREE = abs(ALPHA_STAR_RF / ALPHA_STAR - 1)
print(f"\\n{N_CROSSING_INSIDE} of {len(PAIRS)} pairs cross inside their own tabulated selectivity range: "
      f"{G_CROSS}, at alpha* = {ALPHA_STAR:.2f}")
print(f"(closed form and brentq agree to {ALPHA_STAR_AGREE:.1e} relative - the root-find exists to")
print(" catch a sign or reciprocal slip in the closed form, not to add precision.)")

k0, n0 = prior.k_barrer[G_CROSS], prior.n[G_CROSS]
k1, n1 = present.k_barrer[G_CROSS], present.n[G_CROSS]
sub = PTS[PTS.gas_pair == G_CROSS].sort_values("alpha")
above = sub[sub.alpha > ALPHA_STAR]
print(f"\\nAbove alpha* = {ALPHA_STAR:.2f} the PRESENT {G_CROSS} bound lies BELOW the prior one.")
print(f"{len(above)} of the {len(sub)} points Table 10 tabulates sit in that region:")
for _, r in above.iterrows():
    Pp, Pq = float(P_bound(r.alpha, k1, n1)), float(P_bound(r.alpha, k0, n0))
    print(f"  alpha = {r.alpha:7.4g}  point P = {r.P_fast_barrer:9.4g} | present bound "
          f"{Pp:9.4g} | prior bound {Pq:9.4g}  (present is {100*(Pp/Pq-1):+.1f} % of prior)")
PRESENT_BELOW_PRIOR_PCT = float(100 * (P_bound(sub.alpha.max(), k1, n1)
                                       / P_bound(sub.alpha.max(), k0, n0) - 1))

# How much does this verdict rest on? Exactly one tabulated point, and the page must say so.
_below = sub[sub.alpha < ALPHA_STAR]
print(f"\\nAND HOW THIN IS THAT. 'Inside the paper's own tabulated range' rests on "
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
print(f"\\nThe paper anticipates the cause and says so on p. 396 - the H2/CO2 shift is")
print(f'"primarily a slight slope change" and "the limited number of data points at the lower')
print(f'permeability area of the dataset may have skewed the slope versus the original')
print(f'correlation [2]". Searched the full text for any statement that the two lines CROSS or')
print(f"that the new bound falls below the old one anywhere: none - the words cross/crossing/")
print(f"intersect appear only inside the polymer name 'crosslinked' and in two reference titles.")
print(f"So the crossing is a consequence of what the paper says, not a contradiction of it -")
print(f"but it is a consequence the paper does not draw, and it is the one place in the table")
print(f"where 'upper bound' is not monotone in time.")'''))

cells.append(md(r"""### 3 — Vertical or horizontal? Which measure reproduces the paper's own ranking

The paper sorts the nine pairs into *modest* (five), *minor* (one) and
*significant* (three). Nothing forces that verbal sorting to agree with any
particular arithmetic measure — so it is a real test of which measure Robeson
was actually reasoning with."""))

cells.append(code('''rank_rows = []
for g in PAIRS:
    k0, n0 = prior.k_barrer[g], prior.n[g]
    k1, n1 = present.k_barrer[g], present.n[g]
    P = PTS[PTS.gas_pair == g].P_fast_barrer.to_numpy(float)
    P_gm = float(np.exp(np.mean(np.log(P))))
    # horizontal shift, closed form ...
    dlogA = ((np.log10(P_gm) - np.log10(k1)) / n1) - ((np.log10(P_gm) - np.log10(k0)) / n0)
    # ... and by root-finding each bound's alpha at P_gm independently (no shared algebra)
    a1 = brentq(lambda a: P_bound(a, k1, n1) - P_gm, 1e-6, 1e9, xtol=1e-14, rtol=8.9e-16)
    a0 = brentq(lambda a: P_bound(a, k0, n0) - P_gm, 1e-6, 1e9, xtol=1e-14, rtol=8.9e-16)
    dlogA_rf = np.log10(a1) - np.log10(a0)
    rank_rows.append(dict(pair=g, verdict=CLASS[g], dlogk=SH.dlogk[g], dlogP=SH.shift_gm[g],
                          dlogA=dlogA, dlogA_rf=dlogA_rf, alpha_gain=10 ** dlogA,
                          P_gm=P_gm, n_pres=n1))
RK = pd.DataFrame(rank_rows).set_index("pair").sort_values("dlogA", ascending=False)
DLOGA_ROOTFIND_MAX = float(np.max(np.abs(RK.dlogA - RK.dlogA_rf)))

print("Ranked by SELECTIVITY gain at the pair's own geometric-mean permeability:")
print(f"{'pair':9s} {'paper':12s} {'n_pres':>7s} {'dlog k':>7s} {'dlogP@a_gm':>10s} "
      f"{'dlogA@P_gm':>10s} {'alpha x':>8s}")
for g, r in RK.iterrows():
    print(f"{g:9s} {r.verdict:12s} {r.n_pres:7.3f} {r.dlogk:7.4f} {r.dlogP:10.4f} "
          f"{r.dlogA:10.4f} {r.alpha_gain:8.3f}")
print(f"\\n(closed form vs independent root-find of each bound at P_gm: max |difference| = "
      f"{DLOGA_ROOTFIND_MAX:.2e} decades)")

SIG = [g for g in PAIRS if CLASS[g] == "significant"]
def top3(col):
    return list(RK.sort_values(col, ascending=False).index[:3])
TOP3_A, TOP3_K, TOP3_P = top3("dlogA"), top3("dlogk"), top3("dlogP")
HITS_A, HITS_K, HITS_P = (len(set(TOP3_A) & set(SIG)), len(set(TOP3_K) & set(SIG)),
                          len(set(TOP3_P) & set(SIG)))
gap_A = float(RK.dlogA.iloc[2] - RK.dlogA.iloc[3])
print(f"\\nThe paper's three 'significant' pairs: {', '.join(SIG)}")
print(f"  top 3 by selectivity gain dlogA : {', '.join(TOP3_A):40s} {HITS_A} of 3")
print(f"  top 3 by front factor  dlog k   : {', '.join(TOP3_K):40s} {HITS_K} of 3")
print(f"  top 3 by permeability  dlogP    : {', '.join(TOP3_P):40s} {HITS_P} of 3")
print(f"\\nOnly the horizontal measure reproduces it, and it does so with a clean gap: "
      f"{RK.dlogA.iloc[2]:.4f} vs {RK.dlogA.iloc[3]:.4f}")
print(f"decades between third and fourth ({gap_A:.4f} decades, {100*gap_A/RK.dlogA.iloc[3]:.0f} % of the fourth).")
print(f"\\nThe reason is arithmetic, not taste. CO2/CH4 and He/CH4 have almost the same front-factor")
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
def rank_at(Pfun):
    d = sorted(((g, ((np.log10(Pfun(g)) - np.log10(present.k_barrer[g])) / present.n[g]
                     - (np.log10(Pfun(g)) - np.log10(prior.k_barrer[g])) / prior.n[g]))
                for g in PAIRS), key=lambda z: -z[1])
    return len(set(g for g, _ in d[:3]) & set(SIG)), d[2][1] - d[3][1]
own = lambda f: (lambda g: float(f(PTS[PTS.gas_pair == g].P_fast_barrer.to_numpy(float))))  # noqa
EVAL = [("each pair's own P_min", own(np.min)),
        ("each pair's own P_gm (used above)", own(lambda P: np.exp(np.mean(np.log(P))))),
        ("each pair's own P_median", own(np.median)),
        ("each pair's own P_max", own(np.max)),
        (f"all 117 points' geometric mean, {P_ALL_GM:.0f} barrer", lambda g: P_ALL_GM)]
EVAL += [(f"common P = {c:g} barrer", (lambda g, c=c: float(c)))
         for c in (100.0, 1000.0) ]
print(f"\\nThe evaluation point is a CHOICE here just as alpha_0 was in Result 1 - dlogA varies with")
print(f"P because the lines are not parallel. Every reference below lies inside the permeability")
print(f"range all nine pairs share ({Pmin_all:.0f} to {Pmax_all:.0f} barrers) or is each pair's own:")
print(f"\\n{'evaluation point':44s} {'top-3 hits':>11s} {'gap, decades':>13s}")
_h, _g = [], []
for lab, f in EVAL:
    h, gp = rank_at(f)
    _h.append(h); _g.append(gp)
    print(f"{lab:44s} {h:>7d} of 3 {gp:13.4f}")
EVAL_HITS_MIN = float(min(_h)); EVAL_GAP_MIN = float(min(_g)); EVAL_GAP_MAX = float(max(_g))
print(f"\\nSo BOTH things are true and they must be said together. The RESULT is robust: "
      f"{EVAL_HITS_MIN:.0f} of 3 at")
print(f"every one of these references. The GAP is not: it runs {EVAL_GAP_MIN:.4f} to "
      f"{EVAL_GAP_MAX:.4f} decades across them,")
print(f"a factor {EVAL_GAP_MAX/EVAL_GAP_MIN:.0f}, and it is smallest at each pair's own P_max, where "
      f"H2/CO2 - the pair the")
print(f"paper calls 'minor' - comes within {100*(10**EVAL_GAP_MIN-1):.0f} % of overtaking He/H2. "
      f"Quoting '{gap_A:.3f} decades' as")
print(f"the margin without its reference is the same mistake this page's Result 1 is about.")
print(f"Outside that shared range it does eventually fail, which is a statement about")
print(f"extrapolating past the data rather than about the result:")
for c in (1.0, 10.0, 1e5):
    h, gp = rank_at(lambda g, c=c: c)
    n_out = sum(not (PTS[PTS.gas_pair == g].P_fast_barrer.min() <= c
                     <= PTS[PTS.gas_pair == g].P_fast_barrer.max()) for g in PAIRS)
    print(f"  common P = {c:>8g} barrer lies outside {n_out} of the 9 pairs' own data -> "
          f"{h} of 3")'''))

cells.append(code('''# The picture: the nine comparable pairs, both bounds, and the paper's own near-bound points.
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
    ax.set_title(f"{g}  ({CLASS[g]}, $\\\\Delta\\\\log\\\\alpha$ = {RK.dlogA[g]:.2f})", fontsize=9)
    ax.set_xlabel("$P$ fast gas, barrer", fontsize=8)
    ax.set_ylabel(r"$\\alpha$", fontsize=8)
    ax.tick_params(labelsize=7)
axes.ravel()[0].legend(fontsize=7, loc="lower left")
fig.suptitle("Robeson upper bounds, nine comparable gas pairs, ordered by selectivity gain",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97)); plt.show()'''))

cells.append(md(r"""### 4 — How many of the paper's own "close to the bound" points are above it

Table 1–11 rows are, by their captions, *"experimental data points close to
the present empirical upper bound"*. The line was drawn by eye. So: at each
tabulated point, $r = P/(k\alpha^{n})$ against the pair's own present bound,
and the same against the prior bound."""))

cells.append(code('''ex_rows = []
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

print(f"{'pair':9s} {'N':>3s} {'above':>5s} {'r min':>7s} {'r max':>7s} {'geo-mean r':>10s}")
for g, sub in EX.groupby("gas_pair", sort=False):
    print(f"{g:9s} {len(sub):3d} {int((sub.r_present>1).sum()):5d} {sub.r_present.min():7.3f} "
          f"{sub.r_present.max():7.3f} {np.exp(np.mean(np.log(sub.r_present))):10.3f}")
print(f"\\n{N_ABOVE} of {N_POINTS} tabulated points ({100*FRAC_ABOVE:.1f} %) lie ABOVE the very "
      f"bound they are tabulated against.")
print(f"r ranges {R_MIN:.3f} to {R_MAX:.3f}; geometric mean {R_GEOMEAN:.3f} - so the typical "
      f"tabulated point sits a factor")
print(f"{1/R_GEOMEAN:.1f} BELOW the line in permeability, and the line is a genuine envelope in the "
      f"aggregate.")

print("\\nThe five largest excursions, and what the paper says about them:")
for _, r in EX.sort_values("r_present", ascending=False).head(5).iterrows():
    print(f"  {r.gas_pair:9s} P={r.P:8.4g} alpha={r.alpha:8.4g}  r={r.r_present:6.3f} "
          f"(= {r.alpha_ratio:5.3f}x the bound selectivity)  {r.polymer[:44]}")
Q = EX[(EX.gas_pair == "O2/N2") & (EX.P == 18.0)].iloc[0]
print(f"\\n  The O2/N2 point at P = 18.0, alpha = 9.0 is r = {Q.r_present:.3f} - the second largest")
print(f"  excursion in the whole compilation, and the paper flags exactly it: 'The position of the")
print(f"  one data point above the present upper bound (P(O2) = 18 barrers; alpha(O2/N2) = 9.0) is")
print(f"  questioned as only one significant figure was noted for nitrogen permeability.'")
N_ABOVE_O2N2 = int((EX[EX.gas_pair == 'O2/N2'].r_present > 1).sum())
print(f"  Note the arithmetic sharpens the paper's own wording: it says 'the ONE data point above")
print(f"  the present upper bound' for O2/N2, and Table 1 in fact has {N_ABOVE_O2N2} - the other is")
print(f"  BPDA-ODA at P = 0.079, alpha = 19.8, r = "
      f"{float(EX[(EX.gas_pair=='O2/N2') & (EX.P==0.079)].r_present.iloc[0]):.3f}, at the far low-permeability end")
print(f"  where a hand-drawn line is least constrained. Reported, not repaired.")

N2 = EX[(EX.gas_pair == "N2/CH4") & (EX.P == 153)].iloc[0]
print(f"\\n  And the N2/CH4 point the paper doubts on other grounds (P = 153, alpha = 1.9, the one")
print(f"  whose CH4 permeability it says 'would be more realistic to be 180.2 barrers') is the")
print(f"  LARGEST excursion in Table 9: r = {N2.r_present:.3f}. Consistency check on the tabulated")
print(f"  alpha: 153/80.2 = {153/80.2:.3f}, which rounds to the printed 1.9; 153/180.2 = "
      f"{153/180.2:.3f} does not.")
print(f"  So Table 9 carries the UNREPAIRED value and this page uses it as printed.")'''))

cells.append(code('''# The same measurement against the PRIOR bound: how far the 2008 points sit above the 1991 line.
# This is the shift measured with data instead of with fitted lines.
print(f"{'pair':9s} {'N':>3s} {'above prior':>11s} {'max r_prior':>12s} {'geo-mean r_prior':>17s}")
prior_rows = []
for g in PAIRS:
    sub = EX[EX.gas_pair == g]
    prior_rows.append(dict(pair=g, n=len(sub), n_above=int((sub.r_prior > 1).sum()),
                           rmax=sub.r_prior.max(),
                           gm=float(np.exp(np.mean(np.log(sub.r_prior))))))
    print(f"{g:9s} {len(sub):3d} {int((sub.r_prior>1).sum()):11d} {sub.r_prior.max():12.3g} "
          f"{np.exp(np.mean(np.log(sub.r_prior))):17.3f}")
PR = pd.DataFrame(prior_rows).set_index("pair")
N_ABOVE_PRIOR = int(PR.n_above.sum()); N_PRIOR_POINTS = int(PR.n.sum())
print(f"\\n{N_ABOVE_PRIOR} of {N_PRIOR_POINTS} points on the nine comparable pairs exceed the 1991/94 bound "
      f"({100*N_ABOVE_PRIOR/N_PRIOR_POINTS:.0f} %),")
print(f"against {int((EX[EX.gas_pair.isin(PAIRS)].r_present>1).sum())} that exceed the 2008 one. That is the shift, seen from the data side.")
print(f"He/H2 is the extreme: all {PR.n['He/H2']} of its points beat the 1991 line, the furthest by a factor "
      f"{PR.rmax['He/H2']:.0f}.")
print(f"He/N2 is the quiet one: {PR.n_above['He/N2']} of {PR.n['He/N2']}, max {PR.rmax['He/N2']:.2f} - matching the paper's")
print(f"'The upper bound shift since 1991 is very minor in spite of the much larger dataset")
print(f"available presently.'")
print("\\nCAUTION, and it is not small: these points are 2008 data. That they beat a 1991 line drawn")
print("on 1991 data is expected, not evidence of anything the paper claims. The number is a")
print("MAGNITUDE for the shift, computed without using k or n of the present bound at all - which")
print("is its only real use, and why it appears here rather than among the validations.")'''))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

Four checks — and the first thing to say is which of them is a **second route**
and which is a **consistency identity**, because an earlier draft of this page
got that wrong in the one direction that costs a page its credibility.

**V1 — internal identity.** Tables 13a/13b are computed *from* Table 12, so
reproducing all 26 of their rows is a genuine test of every present $(k,n)$
pair and of the prior pair for He/H2 and H2/CO2. It shares Table 12 as an
*input*, so it cannot detect an error Robeson himself made — it detects errors
*this page* could have made, which is what it is for.

**V2 — a consistency identity, NOT a second route.** Tables 13a and 13b each
evaluate He/H2 and H2/CO2 at one selectivity on **both** bounds, so the ratio
of the two printed permeabilities looks like an independently computed shift.
It is not independent. Robeson's transition permeability *is* $k\alpha^n$
evaluated on the same Table 12 line V1 tests, so for two rows at the same
$\alpha$

$$
\log_{10}\frac{P_1}{P_0}
= \underbrace{(\Delta\log_{10}k + \Delta n\log_{10}\alpha)}_{\text{the decomposition}}
+ \underbrace{\big[\log_{10}(1+\rho_1) - \log_{10}(1+\rho_0)\big]}_{\text{the two V1 residuals}},
$$

with $\rho$ each row's V1 relative deviation. So **V2's residual is the
difference of the two corresponding V1 residuals**, and nothing else. The cell
below measures that rather than asserting it. V2 therefore cannot fail unless
V1 fails on those same four rows; what it *does* have is power over **this
page's own decomposition code**. A natural-log-for-$\log_{10}$ slip changes the
value of many reported metrics and makes exactly one of them *fail*: a
residual that should sit at $10^{-4}$ decades becomes of order one. That is the
only reason V2 is kept, and it is a real one — but it is a unit test of this
page, not evidence about Robeson.

**Result 1's decomposition has no arithmetic second route, and cannot have
one.** Given Table 12, $\Delta\log P(\alpha)=\Delta\log k+\Delta n\log\alpha$
is a rearrangement, not a claim: any quantity built from the same two lines
reproduces it identically. The cell below demonstrates that on a second
candidate — the *data-side* shift, the ratio of the geometric-mean excursions
of a pair's own points above the prior and present bounds — which also comes
out equal to the decomposition to machine precision. Two identities are not
two routes.

**V3 — the second route, and it is a real one.** The one check here that does
*not* take Table 12's present block as an input: a covering line re-estimated
from the paper's own 117 near-bound points by a deterministic envelope
estimator. Substituting it for the printed present line gives an
**envelope-implied shift** that shares no arithmetic and no input with the
decomposition **except two things, both named**: the prior block, which both
routes read from the same restatement, and the reference selectivity
$\alpha_0$, which both routes take from Tables 1–11. The residual between the
two routes is exactly V3's offset at $\alpha_0$, which is why that quantity is
the headline of this section rather than a footnote — and, because it carries
$\log\alpha_0$ explicitly, it is itself reference-dependent for the same reason
the front-factor share is, so it is swept across references and reported with
its reference attached, exactly as Result 3's gap is. It is weak where the
lever arm is short and it is labelled so, with the measurement to show which
pairs those are. It is also the only check on this page whose failure could
indict *Robeson* rather than this page.

**V4** is the break table, and its coverage map is **built from the measured
mover list of every reported metric**, not written by hand."""))

cells.append(md(r"""### V1 — Tables 13a/13b reproduce from Table 12

Three nested checks, each catching something the previous cannot:
the molecular-weight column from molecular weights; the product column from
its two factors; and the permeability column from Table 12's $(k,n)$."""))

cells.append(code('''MW = {"He": 4.0, "H2": 2.0, "O2": 32.0, "N2": 28.0, "CO2": 44.0, "CH4": 16.0}
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
print(f"    -> the convention is identified, not assumed: H2/N2 prints 3.742 = sqrt(28/2), not "
      f"{np.sqrt(IUPAC['N2']/IUPAC['H2']):.3f}.")

# --- (b) the product column ks*(Mj/Mi)^(1/2) --------------------------------------
T13 = T13.assign(prod_calc=lambda d: d.ks * d.sqrt_M_ratio,
                 prod_rel=lambda d: d.ks * d.sqrt_M_ratio / d.ks_D_ratio - 1)
PROD_MAX_REL = float(np.max(np.abs(T13.prod_rel)))
worst_prod = T13.loc[T13.prod_rel.abs().idxmax()]
print(f"\\n(b) product column ks*(Mj/Mi)^(1/2), {len(T13)} rows: max |rel. dev| = {PROD_MAX_REL:.2e}")
print(f"    worst row {worst_prod.table} {worst_prod.gas_pair}: printed {worst_prod.ks_D_ratio}, "
      f"computed {worst_prod.prod_calc:.6g}")
print(f"    (printed to {len(str(worst_prod.ks_D_ratio).split('.')[-1])} decimals, so a rounding "
      f"band of +-{ulp(worst_prod.ks_D_ratio)/worst_prod.ks_D_ratio*100:.3f} % on that cell alone)")'''))

cells.append(code('''# --- (c) the permeability column, from Table 12's (k, n) --------------------------
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
                  + abs(n) * ulp(r.ks_D_ratio) / r.ks_D_ratio + ulp(r.P_transition_barrer) / r.P_transition_barrer)
    # generous band: trailing zeros are placeholders -> 4 significant figures on k and on P
    sig4 = lambda v: 0.5 * 10.0 ** (np.floor(np.log10(abs(v))) - 3)     # noqa: E731
    band_loose = (sig4(k) / k + abs(np.log(r.ks_D_ratio)) * ulp(npr)
                  + abs(n) * ulp(r.ks_D_ratio) / r.ks_D_ratio + sig4(r.P_transition_barrer) / r.P_transition_barrer)
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

print(f"{'tbl':4s} {'pair':9s} {'set':8s} {'P printed':>12s} {'P = k*alpha^n':>14s} "
      f"{'rel dev':>9s} {'tight band':>11s} {'in?':>4s}")
for _, r in V1.iterrows():
    print(f"{r.table:4s} {r.pair:9s} {r.dataset:8s} {r.P_printed:12,.0f} {r.P_calc:14,.0f} "
          f"{100*r.rel:+8.3f}% {100*r.band_tight:10.3f}% {'yes' if r.inside_tight else 'NO':>4s}")
print(f"\\n{len(V1)} rows. Median |deviation| {100*T13_MED_REL:.3f} %, max {100*T13_MAX_REL:.3f} % "
      f"({worst.table} {worst.pair}).")
print(f"Inside the LITERAL printing band: {N_INSIDE_TIGHT}/{len(V1)}.  Inside the GENEROUS "
      f"(4-significant-figure) band: {N_INSIDE_LOOSE}/{len(V1)}.")
out = V1[~V1.inside_loose]
if len(out):
    for _, r in out.iterrows():
        print(f"\\nOUTSIDE EVEN THE GENEROUS BAND: {r.table} {r.pair}, {100*r.rel:+.3f} % against "
              f"+-{100*r.band_loose:.3f} %.")
        src = present if r.dataset == "present" else prior
        k, n = src.k_barrer[r.pair], src.n[r.pair]
        a_needed = float(alpha_bound(r.P_printed, k, n))
        a_print = float(T13[(T13.table == r.table) & (T13.gas_pair == r.pair)].ks_D_ratio.iloc[0])
        print(f"  Backing the printed P out through the same bound needs alpha = {a_needed:.4f}; "
              f"the table prints {a_print}.")
        print(f"  That is a {100*(a_needed/a_print-1):+.2f} % difference in the 4th digit of a "
              f"4-digit cell - a small internal")
        print(f"  inconsistency in one row of Table 13a, reported and not repaired. It changes no")
        print(f"  conclusion here: {r.pair} is a Knudsen-transition estimate, and Table 12's own "
              f"(k, n) for it are")
        print(f"  confirmed by the other table ({100*float(V1[(V1.pair==r.pair)&(V1.table=='13b')].rel.iloc[0]):+.3f} % in 13b).")
print(f"\\nWhat V1 certifies: every present (k, n) in Table 12, and the prior (k, n) for He/H2 and")
print(f"H2/CO2, are transcribed correctly - a wrong digit in any of them would show up here at")
print(f"orders of magnitude above {100*T13_MED_REL:.3f} %. What it cannot certify: anything about "
      f"Table 12 itself, which")
print(f"is an input on both sides, nor the prior (k, n) of the other seven pairs, which Table 13")
print(f"never touches.")'''))

cells.append(md(r"""### V2 — a consistency identity on the decomposition code

For He/H2 and H2/CO2, Tables 13a and 13b each print the transition
permeability **twice** — once on the present bound and once on the prior one,
at *identical* selectivity. Their ratio is $\Delta\log_{10}P(\alpha)$ at that
selectivity, set beside $\Delta\log_{10}k + \Delta n\log_{10}\alpha$.

**What this can and cannot catch.** It can catch a natural-log-for-$\log_{10}$
slip, a sign error, or a prior/present swap in *this page's* decomposition
formula — the break table below shows it moving under all three, and under the
natural-log row it is the only reported quantity that turns from a number into
a visible failure rather than just changing value. It
cannot catch anything about Table 12, and it is not a second route to Result 1:
the printed side is $k\alpha^n$ on the same two lines, so the comparison is an
algebraic tautology whose residual is fixed entirely by V1. The cell measures
that: the four V2 residuals are set against the corresponding differences of V1
residuals, and against a second candidate identity — the shift read off the
*data* rather than off the lines.

What the four rows are genuinely good for is the **reference-dependence**: the
same two pairs, the same two bounds, four different selectivities spanning a
factor 7, and a front-factor share that moves with them. Those are numbers
Robeson printed, not numbers this page chose."""))

cells.append(code('''v2 = []
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
print(f"{'tbl':4s} {'pair':8s} {'alpha':>9s} {'shift, printed':>15s} {'shift, decomposed':>18s} "
      f"{'diff, decades':>14s} {'dlog k':>8s} {'k-share':>8s}")
for _, r in V2.iterrows():
    print(f"{r.table:4s} {r.pair:8s} {r.alpha:9.5f} {r.printed:15.5f} {r.decomposed:18.5f} "
          f"{r['diff']:+14.2e} {r.dlogk:8.4f} {r.k_share:8.3f}")
print(f"\\nWorst |difference| over the four comparisons: {V2_MAX_ABS:.2e} decades.")

# --- the identity, measured rather than asserted --------------------------------------
# V2's residual is log10((1+rho_present)/(1+rho_prior)) where rho is each row's V1 relative
# deviation. Not "correlated with" V1: EQUAL to it, up to floating-point re-association.
print("\\nIs this a second route? No - and here is the arithmetic that settles it. Robeson's")
print("transition permeability IS k*alpha^n on the same two lines V1 tests, so the ratio of the")
print("two printed permeabilities carries the decomposition plus the two V1 residuals and nothing")
print("else. Set V2's residual against the difference of the corresponding V1 residuals:")
print(f"\\n{'tbl':4s} {'pair':8s} {'V2 residual':>16s} {'V1(pres) - V1(prior)':>22s} {'|difference|':>14s}")
_v2id = []
for _, r in V2.iterrows():
    _p = V1[(V1.table == r.table) & (V1.pair == r.pair) & (V1.dataset == "present")].rel.iloc[0]
    _q = V1[(V1.table == r.table) & (V1.pair == r.pair) & (V1.dataset == "prior")].rel.iloc[0]
    _d = float(np.log10((1 + _p) / (1 + _q)))
    _v2id.append(abs(r["diff"] - _d))
    print(f"{r.table:4s} {r.pair:8s} {r['diff']:16.10e} {_d:22.10e} {abs(r['diff']-_d):14.2e}")
V2_IS_V1_MAX_ABS = float(max(_v2id))
print(f"\\nThey agree to {V2_IS_V1_MAX_ABS:.1e} decades - floating-point re-association, not")
print(f"agreement between two computations. V2 cannot fail unless V1 fails on these four rows;")
print(f"the headline {V2_MAX_ABS:.1e} decades is the V1 table above, differenced.")

# --- and the other candidate identity, for the same reason ----------------------------
# The "data-side" shift: how much further a pair's own points sit above the prior bound than
# above the present one, in the geometric mean. It uses no (k, n) ratio explicitly - and it is
# STILL the decomposition, because averaging log r over the points returns log alpha_gm.
print("\\nThe obvious second candidate fails the same way. Take the shift straight off the DATA:")
print("the geometric-mean excursion of a pair's own points above the prior bound, over the same")
print("quantity for the present bound. No Delta log k, no Delta n anywhere in it:")
print(f"\\n{'pair':9s} {'data-side shift':>16s} {'decomposition @ a_gm':>22s} {'|difference|':>14s}")
_ds = []
for g in PAIRS:
    _sub = EX[EX.gas_pair == g]
    _d = float(np.log10(np.exp(np.mean(np.log(_sub.r_prior))) / np.exp(np.mean(np.log(_sub.r_present)))))
    _ds.append(abs(_d - SH.shift_gm[g]))
    print(f"{g:9s} {_d:16.10f} {SH.shift_gm[g]:22.10f} {abs(_d-SH.shift_gm[g]):14.2e}")
DATASIDE_IS_DECOMP_MAX_ABS = float(max(_ds))
print(f"\\nEqual to {DATASIDE_IS_DECOMP_MAX_ABS:.1e} decades, on all nine pairs. The reason is one line:")
print("averaging log r over a pair's points subtracts n*mean(log alpha), and mean(log alpha) IS")
print("log(alpha_gm), so the data-side shift collapses to dlog k + dn*log(alpha_gm) exactly.")
print("\\nSo: given Table 12, Result 1's decomposition is a REARRANGEMENT, and no arithmetic on")
print("Table 12 can corroborate it. The route that can is V3 - which replaces the present line")
print("with one estimated from the 117 points instead of read from Table 12.")

print(f"\\nAnd read the last two columns of the table above. The SAME two pairs, the SAME two bounds, four different")
print(f"selectivities spanning {V2.alpha.max()/V2.alpha.min():.0f}x - and the front-factor share of the shift runs")
print(f"{V2.k_share.min():.2f} to {V2.k_share.max():.2f}. For He/H2 at alpha = {V2[V2.pair=='He/H2'].alpha.iloc[0]:.3f}, "
      f"dlog k OVERSTATES the shift by {100*(V2[V2.pair=='He/H2'].k_share.iloc[0]-1):.1f} %.")
print(f"For H2/CO2 at alpha = {V2[V2.pair=='H2/CO2'].alpha.iloc[0]:.5f}, dlog k accounts for only "
      f"{100*V2[V2.pair=='H2/CO2'].k_share.iloc[0]:.0f} % of it - and at alpha = {ALPHA_STAR:.1f}")
print(f"the same two lines cross and the shift is exactly zero, so the share is undefined there.")
print(f"THIS is what 'the split depends on where you stand' means, demonstrated on numbers the")
print(f"paper printed rather than on numbers this page chose.")'''))

cells.append(md(r"""### V3 — an envelope re-estimated from the paper's own points, and the second route to Result 1

The only check here that touches the **empirical** content rather than the
arithmetic, and the only one that does not take Table 12's present block as an
input. The estimator is stated so it can be argued with: the line
$\log P = \log k + n\log\alpha$ minimising
$\sum_i (\log P_i - \log k - n\log\alpha_i)^2$ **subject to every point lying
on or below it**. That is a convex quadratic program; its optimum activates
at most two constraints, so it is solved here by **exact enumeration of the
active sets** rather than by an iterative optimiser — no start point, no
tolerance, no seed, bit-reproducible.

It is *not* what "by eye" means. A covering line is pinned by two or three
extreme points; an eye fits the body of a cloud. The two estimators answer
different questions, and this page keeps them apart — **nothing below is
offered as a correction to the printed $(k,n)$.**

**Why this is nevertheless Result 1's second route.** Substitute
$(k_{\text{env}}, n_{\text{env}})$ for the printed present line and recompute
the shift at the same reference:

$$
\Delta\log_{10}P_{\text{env}}(\alpha_0)
= \log_{10}\frac{k_{\text{env}}}{k_0} + (n_{\text{env}} - n_0)\log_{10}\alpha_0 .
$$

Its prior half is shared with the decomposition; its present half comes from
117 points and a quadratic program instead of from Table 12. Subtracting, the
prior line cancels exactly, so the **residual between the two routes is
V3's offset at $\alpha_0$** — the column already tabulated below. That is why
the offset is the headline of this section: it is not a curiosity about a
re-fit, it is the disagreement between the only two ways this page can compute
how far the bound moved. Its scope is exactly the present line; it says
nothing about the prior block, which both routes read from the same place, and
nothing about $\alpha_0$, which both routes take from the same 117 points.

**And it carries $\log\alpha_0$, so it is reference-dependent too.** That is
not a caveat bolted on: it is this page's own thesis pointed at this page's own
newest headline. Quoting "the two routes agree to *x* decades" without naming
$\alpha_0$ is the same error as quoting a front-factor share without naming
$\alpha_0$. So the reference is swept below — over each pair's own minimum,
geometric-mean, median and maximum tabulated selectivity, all four inside the
pair's own data — and the residual is reported with its reference attached
wherever it appears, on this page and on every surface that quotes it."""))

cells.append(code('''def envelope_fit(x, y, tol=1e-10):
    """Tightest covering line y = c + n*x subject to y_i - (c + n*x_i) <= 0 for all i,
    minimising the sum of squared residuals.

    Convex QP in two unknowns: the optimum activates 0, 1 or 2 constraints, so the
    exact solution is the best feasible member of a finite candidate list -
    the unconstrained OLS line, every line through one point that is optimal in
    slope given that point is active, and every line through two points.
    No iteration, no start point, no tolerance to tune."""
    N = len(x)
    n_ols, c_ols = np.polyfit(x, y, 1)
    cands = [(c_ols, n_ols)]
    for i in range(N):
        dx, dy = x - x[i], y - y[i]
        s = float(np.sum(dx * dx))
        n1 = float(np.sum(dx * dy) / s) if s > 0 else 0.0
        cands.append((y[i] - n1 * x[i], n1))                       # point i active
        for j in range(i + 1, N):
            if x[j] != x[i]:
                n2 = (y[j] - y[i]) / (x[j] - x[i])
                cands.append((y[i] - n2 * x[i], n2))               # points i and j active
    scale = tol * max(1.0, float(np.max(np.abs(y))))
    best = None
    for c, n in cands:
        r = y - (c + n * x)
        if float(np.max(r)) <= scale:
            sse = float(np.sum(r * r))
            if best is None or sse < best[0]:
                best = (sse, c, n)
    assert best is not None, "no feasible covering line - impossible for a finite point set"
    return 10.0 ** best[1], best[2], best[0]

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
print("Sorted by how many decades of selectivity the pair's own tabulated points span.")
print(f"{'pair':9s} {'alpha span, dec':>15s} {'n printed':>10s} {'n envelope':>11s} {'dn':>8s} "
      f"{'offset @ a_gm, dec':>19s}")
for g, r in V3.iterrows():
    print(f"{g:9s} {r.span:15.3f} {r.n_pr:10.4f} {r.n_env:11.4f} {r.dn:+8.4f} {r.offset_gm:+19.4f}")
print(f"\\nMedian |dn| = {V3_DN_MED:.3f}, worst {V3_DN_MAX:.3f}. Vertical offset at each pair's own")
print(f"geometric-mean selectivity: median {V3_OFFSET_MED:+.3f} decades, largest magnitude "
      f"{V3_OFFSET_MAX:.3f}.")
print(f"\\nTwo things, and only two, come out of this.")
print(f"\\n(i) The offsets are SMALL and of BOTH SIGNS. Over the body of each cloud the printed line")
print(f"and the tightest covering line sit within a factor {10**V3_OFFSET_MAX:.1f} of each other, usually far")
print(f"closer - so the {N_ABOVE} excursions of Result 4 are individual points poking through a line that")
print(f"is in the right PLACE, not a line drawn in the wrong place.")
# Stated as a SET match, computed, because the obvious rank-by-rank version of this sentence
# is false: O2/N2 has the 3rd shortest span and the 2nd largest |dn|, so "the two shortest
# spans carry the two largest |dn|" - which an earlier version of this page asserted - is not
# true. The largest m for which the two orderings agree AS SETS is the honest statement.
_by_span = list(V3.span.sort_values().index)
_by_dn = list(V3.dn.abs().sort_values(ascending=False).index)
M_LEVER = max(m for m in range(1, len(V3)) if set(_by_span[:m]) == set(_by_dn[:m]))
print(f"\\n(ii) The slope disagreement is entirely a question of LEVER ARM, and it lands exactly")
print(f"where the paper itself says it should. The {M_LEVER} pairs whose points span the fewest decades of")
print(f"selectivity are exactly the {M_LEVER} carrying the largest |dn| - as SETS, not rank by rank "
      f"(O2/N2 has the")
print(f"3rd shortest span and the 2nd largest |dn|). Those {M_LEVER} are "
      f"{', '.join(f'{g} ({V3.span[g]:.2f} dec, |dn| {abs(V3.dn[g]):.3f})' for g in _by_span[:M_LEVER])}.")
print(f"The {int((V3.span > 1.5).sum())} pairs spanning more than 1.5 decades")
print(f"all agree to |dn| <= {V3[V3.span > 1.5].dn.abs().max():.3f}. For He/H2 - the paper's headline pair, and the one")
print(f"whose twelve points cover only {V3.span['He/H2']:.2f} decades of alpha - the tabulated points simply do")
print(f"not determine a slope, and Robeson says so: 'While the front factor has moved")
print(f"significantly, the slope does not appear to have changed although sufficient data does not")
print(f"exist to clearly confirm that observation.' V3 is the arithmetic behind that sentence.")
print(f"\\nWhat V3 is NOT: evidence that any printed (k, n) is wrong. It cannot be - it optimises a")
print(f"different objective on a subset of the data Robeson drew the line through (Tables 1-11 are")
print(f"the near-bound points, not the ~300 papers' worth of cloud), and for O2/N2 it is pinned by")
print(f"the very point the paper says it distrusts.")'''))

cells.append(code('''# ---- the second route to Result 1, run against the decomposition -----------------------
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
assert np.allclose(SR.resid, V3.offset_gm[SR.index], rtol=0, atol=1e-12), \\
    "the two-route residual must BE V3's offset at alpha_gm - the prior line cancels"
SR_WIDE = SR[SR.span > 1.5]
SECOND_ROUTE_MAX = float(SR.resid.abs().max())
SECOND_ROUTE_MED = float(SR.resid.abs().median())
SECOND_ROUTE_WIDE_MAX = float(SR_WIDE.resid.abs().max())
SECOND_ROUTE_SIGN_AGREE = int((np.sign(SR.sh_env) == np.sign(SR.sh_dec)).sum())
SHARE_ENV_MED = float(SR.share_env.median())

print("The shift, computed twice: from Table 12 (the decomposition) and from Tables 1-11 (the")
print("envelope) - sorted by how many decades of selectivity the pair's own points span, because")
print("that is the lever arm the second route is built on.")
print("EVERY number in this table is evaluated at each pair's own geometric-mean selectivity a_gm,")
print("and none of them means anything without that - see the reference sweep below.")
print(f"{'pair':9s} {'span, dec':>10s} {'shift, decomp':>14s} {'shift, envelope':>16s} "
      f"{'resid @ a_gm':>12s} {'k-share, decomp':>16s} {'k-share, env':>13s}")
for g, r in SR.iterrows():
    print(f"{g:9s} {r.span:10.2f} {r.sh_dec:14.4f} {r.sh_env:16.4f} {r.resid:+12.4f} "
          f"{r.share_dec:16.3f} {r.share_env:13.3f}")
print(f"\\nAT a_gm. Sign of the shift: the two routes agree on {SECOND_ROUTE_SIGN_AGREE} of {len(SR)} "
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
print(f"\\nAnd now the sharp part, which is the page's central claim arriving from the other side.")
print(f"The two routes agree about the SIZE of the shift and disagree completely about its SPLIT:")
print(f"the front-factor share has median {SH.share.median():.2f} by the decomposition and "
      f"{SHARE_ENV_MED:.2f} by the envelope,")
print(f"ranging {SR.share_env.min():.2f} to {SR.share_env.max():.2f} instead of "
      f"{SH.share.min():.2f} to {SH.share.max():.2f}. That is not a defect in either route.")
print(f"The share is dlog k / dlog P, and dlog k is the shift at alpha = 1 - one to four decades")
print(f"outside every dataset, along a slope the points barely determine. Change the line by an")
print(f"amount too small to see over the data and the extrapolated intercept moves by everything.")
print(f"A quantity that survives one estimator and not the other is not a property of the data.")

# ---- what a NEGATIVE front-factor share is, since the range contains one ---------------
_gneg = SR.share_env.idxmin()
_k0neg, _kEneg = prior.k_barrer[_gneg], V3.k_env[_gneg]
_dlkneg = np.log10(_kEneg) - np.log10(_k0neg)
_across = 10.0 ** (-_dlkneg / (V3.n_env[_gneg] - prior.n[_gneg]))
print(f"\\nOne of those envelope shares is NEGATIVE - {SR.share_env.min():.4f}, on {_gneg} - and a "
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
print(f"A share outside [0, 1] is not a defect - it is what happens when dlog k is an intercept")
print(f"extrapolated to a selectivity {np.log10(SH.a_gm[_gneg]):.1f} decades outside the data, which is "
      f"exactly Result 1's warning.")

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
print(f"\\nAnd now the same discipline applied to THIS section's own headline. The residual is")
print(f"(log k_env - log k1) + (n_env - n1) log a0: it carries log a0 explicitly, so '{SECOND_ROUTE_WIDE_MAX:.3f} decades'")
print(f"is a property of the two present lines AT ONE SELECTIVITY, not of the two lines. Swept:")
print(f"\\n{'reference selectivity':40s} {'max|resid|, wide':>17s} {'median|resid|, all 9':>21s} "
      f"{'sign agreement':>15s}")
_rw, _rm, _rs = [], [], []
for lab, af in A0:
    w, m, s, worst = route_at(af)
    _rw.append(w); _rm.append(m); _rs.append(s)
    print(f"{lab:40s} {w:17.4f} {m:21.4f} {s:12d} of {len(PAIRS)}")
    if worst is not None:
        print(f"{'':40s}   ^ sign disagreement: {worst[0]} envelope {worst[1]:+.4f} vs "
              f"decomposition {worst[2]:+.4f}")
SR_WIDE_MAX_OVER_A0 = float(max(_rw)); SR_SIGN_MIN_OVER_A0 = float(min(_rs))
print(f"\\nSo both halves of this section's headline have to be quoted with their reference. The")
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
print(f"\\nResult 3 gets a second route out of the same substitution, at the SAME reference:")
print(f"  ranked by selectivity gain with the ENVELOPE present line: {', '.join(g for g, _ in rk_env)}")
print(f"  the paper's three 'significant' pairs in the top three: {HITS_ENV} of 3 - unchanged, on a")
print(f"  present line that never touches Table 12.")
# But "survives it" is a robustness claim, and this page already owns a robustness test for
# exactly this quantity: EVAL, the seven references Result 3's primary route is held to. The
# substituted route gets the identical test rather than a kinder one.
print(f"\\nThat is one evaluation point. Result 3's PRIMARY route is held to the {len(EVAL)} references of")
print(f"EVAL above, so the substituted route gets the same {len(EVAL)}, not a kinder test:")
print(f"\\n{'evaluation point':44s} {'envelope':>9s} {'primary':>9s}   top 3, envelope route")
_he, _hp = [], []
for lab, f in EVAL:
    d = sorted(((g, ((np.log10(f(g)) - np.log10(V3.k_env[g])) / V3.n_env[g]
                     - (np.log10(f(g)) - np.log10(prior.k_barrer[g])) / prior.n[g]))
                for g in PAIRS), key=lambda z: -z[1])
    h = len(set(g for g, _ in d[:3]) & set(SIG))
    hp = rank_at(f)[0]
    _he.append(h); _hp.append(hp)
    print(f"{lab:44s} {h:>4d} of 3 {hp:>4d} of 3   {', '.join(g for g, _ in d[:3])}")
HITS_ENV_MIN = float(min(_he))
_fail = [lab for (lab, _), h in zip(EVAL, _he) if h < max(_he)]
print(f"\\nSo the honest sentence is NOT 'it survives'. The substituted route reproduces the ranking")
print(f"{HITS_ENV:.0f} of 3 at the reference used above and at {sum(h == HITS_ENV for h in _he) - 1} of the "
      f"other {len(EVAL) - 1} references, and {int(min(_he))} of 3 at the")
print(f"remaining {len(_fail)} - {'; '.join(_fail)} - where He/H2 drops out of the top three. The PRIMARY")
print(f"route holds {int(min(_hp))} of 3 at all {len(EVAL)}. That contrast is the interesting result and it is")
print(f"the reported one: the corroboration is real but weaker than the route it corroborates,")
print(f"which is what you would expect of a present line estimated from {len(PTS)} points rather than")
print(f"read off a table.")'''))

cells.append(md("""### V4 — break table

Every metric reported below needs something that moves it. The rows are the
misreadings and sign slips this page could actually have made, plus one row
that is not a defect at all.

Two things about this table are worth stating before it runs, because an
earlier draft of this page got both wrong.

**It recomputes all of them, not a proxy.** `rerun_all` below rebuilds *every
one of the reported metrics* from the (possibly corrupted) inputs, by a code
path separate from the narrative cells above — the clean row is asserted equal
to the reported values key for key. That is a second implementation of the
page's **assembly**, and the word matters: `rerun_all` reads none of the
narrative result frames, but it does call the same shared primitives
(`P_bound`, `alpha_bound`, `excursion`, `ulp`, `sqrtM` and — the most intricate
code here — `envelope_fit`), which are therefore implemented **once**, not
twice. Three raw-data-derived intermediates are reused too, including V2's
printed columns, so a break row that perturbed Table 13 would be invisible to
`rerun_all`'s V2 block; no current row does. The coverage map at the end is then **built from the measured
mover list of each metric**. A metric cannot claim a row that does not move it,
because nothing writes those lists by hand. The earlier draft asserted only
that the map's *keys* matched — which is why **26 of its 48 metrics carried at
least one false attribution, 38 false row-codes out of 109 claimed**, several
of them structurally impossible (the block swap flips numerator and denominator
of every front-factor share, so it cannot move one; the natural-log row scales
both alike; $\\alpha^\\star$ is symmetric in the two lines, so the swap cannot
change the crossing *count*). One entry contradicted the page's own printed
output two cells above it. Key-set equality cannot see any of that; measured
movement can.

**It perturbs Tables 1–11 as well as Table 12.** The earlier draft's ten rows
all corrupted Table 12 or a convention, leaving the 117 near-bound points — the
input that sets every $\\alpha_0$, every excursion, the whole envelope, and the
*in-range* test that carries Result 2 — with no break row at all. Three rows
now perturb a single cell of Tables 1–11, and one of them destroys a
headline."""))

cells.append(code('''def rerun_all(T12x=None, PTSx=None, sign=1.0, logfun=np.log10,
              swap_blocks=False, flip_r=False):
    """Recompute EVERY reported metric from (possibly corrupted) inputs.

    A separate code path from the narrative cells above: the clean call is asserted
    equal to `metrics` key for key in the coverage cell, and each break row's MEASURED
    effect on each metric is what the coverage map is built from. Three extra diagnostics
    (prefixed `_`) are returned for the resolving-power discussion and are not reported.

    Second implementation of the page's ASSEMBLY, not of everything: this function reads no
    narrative result frame, but it calls the same P_bound / alpha_bound / excursion / ulp /
    sqrtM / envelope_fit, which are written once."""
    t = T12 if T12x is None else T12x
    pts = PTS if PTSx is None else PTSx
    lg = logfun
    pri = t[t.dataset == ("present" if swap_blocks else "prior")].set_index("gas_pair")
    pre = t[t.dataset == ("prior" if swap_blocks else "present")].set_index("gas_pair")
    pri = pri.assign(n=sign * pri.n)
    pre = pre.assign(n=sign * pre.n)
    pairs_all = list(t[t.dataset == "prior"].gas_pair)
    pairs = [g for g in pairs_all
             if np.isfinite(pri.k_barrer[g]) and np.isfinite(pre.k_barrer[g])]
    o = {"pairs_in_table12": float(len(pairs_all)), "pairs_comparable": float(len(pairs))}

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
    o["pairs_near_parallel"] = float((sh.spread_over_shift < 0.20).sum())
    o["heh2_shift_decades"] = float(sh.shift_gm["He/H2"])
    o["heh2_dlogk"] = float(sh.dlogk["He/H2"])
    o["heh2_slope_change_frac"] = float(abs(sh.dn["He/H2"] / pri.n["He/H2"]))
    o["heh2_spread_over_shift"] = float(sh.spread_over_shift["He/H2"])
    o["h2ch4_spread_over_shift"] = float(sh.spread_over_shift["H2/CH4"])
    o["h2co2_spread_over_shift"] = float(sh.spread_over_shift["H2/CO2"])

    # ---- Result 2 -------------------------------------------------------------------
    cr = []
    for g in pairs:
        dlk, dn = sh.dlogk[g], sh.dn[g]
        lo, hi = np.log10(sh.a_lo[g]), np.log10(sh.a_hi[g])
        f = lambda la, dlk=dlk, dn=dn: dlk + dn * la          # noqa: E731
        ins = bool(f(lo) * f(hi) < 0)
        cr.append(dict(pair=g, a_star=10.0 ** (-dlk / dn) if dn != 0 else np.nan, inside=ins,
                       a_rf=10.0 ** brentq(f, lo, hi, xtol=1e-14, rtol=1e-15) if ins else np.nan))
    cx = pd.DataFrame(cr).set_index("pair")
    o["pairs_crossing_inside_range"] = float(cx.inside.sum())
    if cx.inside.any():
        gx = cx.index[cx.inside][0]
        o["h2co2_crossing_alpha"] = float(cx.a_star[gx])
        o["h2co2_crossing_alpha_rootfind_rel"] = float(abs(cx.a_rf[gx] / cx.a_star[gx] - 1))
    else:                       # the crossing result can be destroyed outright - say so
        gx = "H2/CO2"
        o["h2co2_crossing_alpha"] = float("nan")
        o["h2co2_crossing_alpha_rootfind_rel"] = float("nan")
    amax = pts[pts.gas_pair == gx].alpha.max()
    o["h2co2_present_over_prior_at_alpha_max_pct"] = float(100 * (
        pre.k_barrer[gx] * amax ** pre.n[gx] / (pri.k_barrer[gx] * amax ** pri.n[gx]) - 1))

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
    rk = pd.DataFrame(rr).set_index("pair").sort_values("dlogA", ascending=False)
    hit = lambda c: float(len(set(rk.sort_values(c, ascending=False).index[:3]) & set(SIG)))
    o["dlogA_top3_hits"] = hit("dlogA")
    o["dlogk_top3_hits"] = hit("dlogk")
    o["dlogP_top3_hits"] = hit("dlogP")
    o["dlogA_gap_third_to_fourth"] = float(rk.dlogA.iloc[2] - rk.dlogA.iloc[3])
    o["dlogA_rootfind_max_abs"] = float(np.max(np.abs(rk.dlogA - rk.dlogA_rf)))
    o["alpha_gain_hech4"] = float(rk.alpha_gain["He/CH4"])
    o["alpha_gain_co2ch4"] = float(rk.alpha_gain["CO2/CH4"])
    _gm = float(np.exp(np.mean(np.log(pts.P_fast_barrer.to_numpy(float)))))
    _own = lambda f: (lambda g: float(f(pts[pts.gas_pair == g].P_fast_barrer.to_numpy(float))))
    _h, _gp = [], []
    for f in ([_own(np.min), _own(lambda P: np.exp(np.mean(np.log(P)))), _own(np.median),
               _own(np.max), lambda g: _gm]
              + [(lambda g, c=c: float(c)) for c in (100.0, 1000.0)]):
        d = sorted(((g, ((lg(f(g)) - lg(pre.k_barrer[g])) / pre.n[g]
                         - (lg(f(g)) - lg(pri.k_barrer[g])) / pri.n[g])) for g in pairs),
                   key=lambda z: -z[1])
        _h.append(len(set(g for g, _ in d[:3]) & set(SIG))); _gp.append(d[2][1] - d[3][1])
    o["dlogA_top3_hits_min_over_eval_points"] = float(min(_h))
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
        for rr1, rr0 in zip(r1, r0):
            er.append(dict(gas_pair=g, r_present=rr1, r_prior=rr0))
    ex = pd.DataFrame(er)
    o["points_total"] = float(len(ex))
    o["points_above_present_bound"] = float((ex.r_present > 1).sum())
    o["frac_above_present_bound"] = float((ex.r_present > 1).sum() / len(ex))
    o["excursion_max"] = float(ex.r_present.max())
    o["excursion_min"] = float(ex.r_present.min())
    o["excursion_geomean"] = float(np.exp(np.mean(np.log(ex.r_present))))
    o["points_above_prior_bound"] = float((ex[ex.gas_pair.isin(pairs)].r_prior > 1).sum())

    # ---- V1 -------------------------------------------------------------------------
    o["t13_mw_column_max_abs_rel"] = float(np.max(np.abs(
        np.array([sqrtM(g) for g in T13.gas_pair]) / T13.sqrt_M_ratio - 1)))
    o["t13_mw_column_max_abs_rel_iupac"] = float(np.max(np.abs(sqrtM_iupac / T13.sqrt_M_ratio - 1)))
    o["t13_product_column_max_abs_rel"] = float(np.max(np.abs(
        T13.ks * T13.sqrt_M_ratio / T13.ks_D_ratio - 1)))
    vr = []
    for _, r in T13.iterrows():
        src = pri if r.dataset == "prior" else pre
        k, n = src.k_barrer[r.gas_pair], src.n[r.gas_pair]
        kp, npr = src.k_as_printed[r.gas_pair], src.n_as_printed[r.gas_pair]
        rel = k * r.ks_D_ratio ** n / r.P_transition_barrer - 1
        common = (abs(np.log(r.ks_D_ratio)) * ulp(npr)
                  + abs(n) * ulp(r.ks_D_ratio) / r.ks_D_ratio)
        s4 = lambda v: 0.5 * 10.0 ** (np.floor(np.log10(abs(v))) - 3)     # noqa: E731
        vr.append(dict(key=(r.table, r.gas_pair, r.dataset), rel=rel,
                       bt=ulp(kp) / k + common + ulp(r.P_transition_barrer) / r.P_transition_barrer,
                       bl=s4(k) / k + common + s4(r.P_transition_barrer) / r.P_transition_barrer))
    v1 = pd.DataFrame(vr)
    o["t13_identity_max_abs_rel"] = float(v1.rel.abs().max())
    o["t13_identity_median_abs_rel"] = float(v1.rel.abs().median())
    o["t13_rows_inside_tight_band"] = float((v1.rel.abs() <= v1.bt).sum())
    o["t13_rows_inside_loose_band"] = float((v1.rel.abs() <= v1.bl).sum())
    o["t13a_o2n2_alpha_backout_rel"] = float(
        (_r13a.P_transition_barrer / pre.k_barrer["O2/N2"]) ** (1.0 / pre.n["O2/N2"])
        / _r13a.ks_D_ratio - 1)
    o["_v1_heh2"] = float(abs(v1.set_index("key").rel[("13a", "He/H2", "present")]))
    o["_v1_h2co2_prior"] = float(abs(v1.set_index("key").rel[("13a", "H2/CO2", "prior")]))
    o["_sr_resid_hech4"] = np.nan          # filled in the V3 block below

    # ---- V2 -------------------------------------------------------------------------
    kk = []
    for _, r in V2.iterrows():
        d = sh.dlogk[r.pair] + sh.dn[r.pair] * lg(r.alpha)
        kk.append(dict(diff=abs(d - r.printed), k_share=sh.dlogk[r.pair] / d))
    v2 = pd.DataFrame(kk)
    o["v2_shift_crosscheck_max_abs_decades"] = float(v2["diff"].max())
    o["v2_kshare_min"] = float(v2.k_share.min())
    o["v2_kshare_max"] = float(v2.k_share.max())

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
    o["envelope_dn_max_abs_wide_span"] = float(v3[v3.span > 1.5].dn.abs().max())
    o["envelope_span_min_decades"] = float(v3.span.min())
    # The SECOND-ROUTE RESIDUAL, and it is computed here as the residual - sh_env minus the
    # decomposition's own shift, over the COMPARABLE pairs - rather than as v3.offset_gm over the
    # point-carrying ones. The two are the same object at the clean point (that identity is
    # asserted in the second-route cell) and the two domains coincide only by luck, so an earlier
    # version reported one quantity and measured the movers of the other; under (L) they differ by
    # an order of magnitude and the clean-point assert cannot see it.
    sr_env, sr_res, sr_sgn = [], [], []
    for g in pairs:
        sh_env = ((np.log10(v3.k_env[g]) - np.log10(pri.k_barrer[g]))
                  + (v3.n_env[g] - pri.n[g]) * np.log10(sh.a_gm[g]))
        sr_env.append((np.log10(v3.k_env[g]) - np.log10(pri.k_barrer[g])) / sh_env)
        sr_res.append((g, abs(sh_env - sh.shift_gm[g])))
    o["envelope_offset_max_abs_wide_span"] = float(max(r for g, r in sr_res if v3.span[g] > 1.5))
    o["_sr_resid_hech4"] = float(dict(sr_res).get("He/CH4", np.nan))
    o["envelope_implied_share_median"] = float(np.median(sr_env))
    o["envelope_implied_share_min"] = float(np.min(sr_env))
    o["envelope_implied_share_max"] = float(np.max(sr_env))
    # the residual's own reference sweep - the same four in-data references as the narrative
    _a0 = [lambda g: float(pts[pts.gas_pair == g].alpha.min()),
           lambda g: sh.a_gm[g],
           lambda g: float(np.median(pts[pts.gas_pair == g].alpha.to_numpy(float))),
           lambda g: float(pts[pts.gas_pair == g].alpha.max())]
    _w, _s = [], []
    for af in _a0:
        wide, sgn = [], 0
        for g in pairs:
            a0 = af(g)
            se = ((np.log10(v3.k_env[g]) - np.log10(pri.k_barrer[g]))
                  + (v3.n_env[g] - pri.n[g]) * np.log10(a0))
            sd = sh.dlogk[g] + sh.dn[g] * lg(a0)
            sgn += int(np.sign(se) == np.sign(sd))
            if v3.span[g] > 1.5:
                wide.append(abs(se - sd))
        _w.append(max(wide)); _s.append(sgn)
    o["envelope_offset_max_abs_wide_span_max_over_alpha0"] = float(max(_w))
    o["envelope_sign_agree_min_over_alpha0"] = float(min(_s))
    re, _rehits = [], []
    for g in pairs:
        P_gm = float(np.exp(np.mean(np.log(pts[pts.gas_pair == g].P_fast_barrer.to_numpy(float)))))
        re.append((g, (np.log10(P_gm) - np.log10(v3.k_env[g])) / v3.n_env[g]
                   - (np.log10(P_gm) - np.log10(pri.k_barrer[g])) / pri.n[g]))
    re.sort(key=lambda z: -z[1])
    o["envelope_implied_dlogA_top3_hits"] = float(len(set(g for g, _ in re[:3]) & set(SIG)))
    # ... and under the SAME seven references the primary route is swept over
    for f in ([_own(np.min), _own(lambda P: np.exp(np.mean(np.log(P)))), _own(np.median),
               _own(np.max), lambda g: _gm]
              + [(lambda g, c=c: float(c)) for c in (100.0, 1000.0)]):
        d = sorted(((g, (np.log10(f(g)) - np.log10(v3.k_env[g])) / v3.n_env[g]
                     - (np.log10(f(g)) - np.log10(pri.k_barrer[g])) / pri.n[g]) for g in pairs),
                   key=lambda z: -z[1])
        _rehits.append(len(set(g for g, _ in d[:3]) & set(SIG)))
    o["envelope_implied_dlogA_top3_hits_min_over_eval_points"] = float(min(_rehits))

    _a = np.array([0.2, 1.0, 7.5, 900.0])
    o["bound_roundtrip_max_rel"] = max(
        float(np.max(np.abs(alpha_bound(P_bound(_a, pre.k_barrer[g], pre.n[g]),
                                        pre.k_barrer[g], pre.n[g]) / _a - 1))) for g in pairs)
    return o


# ---- the rows -----------------------------------------------------------------------
BASE = rerun_all()
ROWS = [("(S) n taken POSITIVE (sign convention flipped)", dict(sign=-1.0)),
        ("(B) prior and present blocks swapped", dict(swap_blocks=True)),
        ("(L) natural log used in the decomposition", dict(logfun=np.log)),
        ("(I) excursion inverted: k*alpha^n / P", dict(flip_r=True))]
for tag, label, (ds, g, col, val, printed) in [
        ("C", "thousands comma dropped: k(O2/N2, present) 1,396,000 -> 1396",
         ("present", "O2/N2", "k_barrer", 1396.0, "1396")),
        ("G9", "glyph 9<->0: k(He/H2, present) 59,910 -> 59,010",
         ("present", "He/H2", "k_barrer", 59010.0, "59,010")),
        ("G8", "glyph 8<->3: n(He/CH4, present) -0.809 -> -0.309",
         ("present", "He/CH4", "n", -0.309, "-0.309")),
        ("DD", "digit dropped: k(He/CH4, present) 19,800 -> 1,980",
         ("present", "He/CH4", "k_barrer", 1980.0, "1,980")),
        ("D1", "last digit: n(He/H2, present) -4.864 -> -4.865",
         ("present", "He/H2", "n", -4.865, "-4.865")),
        ("D2", "last digit: k(H2/CO2, prior) 1200 -> 1201",
         ("prior", "H2/CO2", "k_barrer", 1201.0, "1201")),
        ("DP", "digit dropped in the PRIOR block: k(He/H2, prior) 960 -> 96",
         ("prior", "He/H2", "k_barrer", 96.0, "96"))]:
    t = T12.copy()
    m = (t.dataset == ds) & (t.gas_pair == g)
    t.loc[m, col] = val
    t.loc[m, "k_as_printed" if col == "k_barrer" else "n_as_printed"] = printed
    ROWS.append((f"({tag}) {label}", dict(T12x=t)))
# the three Tables 1-11 rows: one alpha that carries a headline, one P, one that sets a span
for tag, label, (gp, col, old, new) in [
        ("T1", "TABLE 10 decimal slip: alpha(H2/CO2) 100.9 -> 10.09",
         ("H2/CO2", "alpha", 100.9, 10.09)),
        ("T2", "TABLE 1 decimal slip: P(O2/N2) 18.0 -> 1.80",
         ("O2/N2", "P_fast_barrer", 18.0, 1.80)),
        ("T3", "TABLE 7 glyph 3<->8: alpha(He/H2) 4.39 -> 4.89",
         ("He/H2", "alpha", 4.39, 4.89))]:
    p = PTS.copy()
    m = (p.gas_pair == gp) & (p[col] == old)
    assert int(m.sum()) == 1, f"{tag} must hit exactly one cell, hit {int(m.sum())}"
    p.loc[m, col] = new
    ROWS.append((f"({tag}) {label}", dict(PTSx=p)))

BR = pd.DataFrame([BASE] + [rerun_all(**kw) for _, kw in ROWS],
                  index=["(clean) - as transcribed"] + [lab for lab, _ in ROWS])

SHOW = ["heh2_shift_decades", "t13_identity_max_abs_rel", "_v1_heh2", "_v1_h2co2_prior",
        "v2_shift_crosscheck_max_abs_decades", "h2co2_crossing_alpha",
        "envelope_offset_max_abs_decades", "points_above_present_bound", "dlogA_top3_hits"]
HEAD = ["He/H2 shift", "V1 max", "V1 He/H2", "V1 H2/CO2pr", "V2 max", "alpha*",
        "V3 offset", "#above", "top3"]
def _fmt(v):
    if not np.isfinite(v):
        return "        n/a"
    return f"{v:11.4f}" if abs(v) >= 1e-3 or v == 0 else f"{v:11.2e}"
_W = max(len(l) for l in BR.index)
print(f"{'injected defect':{_W}s} " + " ".join(f"{h:>11s}" for h in HEAD))
for lab in BR.index:
    print(f"{lab:{_W}s} " + " ".join(_fmt(BR[c][lab]) for c in SHOW))

# ---- the measured mover lists, which are what the coverage map is built from ----------
ABS_FLOOR = 1e-12
def movers(col):
    """The injected defects MEASURED to move a metric. Movement below ABS_FLOOR is
    floating-point re-association, not detection, so it does not count."""
    b = BR[col].iloc[0]
    out = []
    for lab in BR.index[1:]:
        v = BR[col][lab]
        if np.isfinite(b) != np.isfinite(v):
            out.append(lab.split(")")[0][1:] + " (destroys it)")
        elif np.isfinite(b) and abs(v - b) > max(ABS_FLOOR, 1e-12 * abs(b)):
            out.append(lab.split(")")[0][1:])
    return out
BR_MOVERS = {c: movers(c) for c in BR.columns if not c.startswith("_")}
TAGS = [lab.split(")")[0][1:] for lab in BR.index[1:]]
n_moved = {t: sum(any(m.split(" ")[0] == t for m in v) for v in BR_MOVERS.values()) for t in TAGS}
print(f"\\nHow many of the {len(BR_MOVERS)} reported metrics each row moves:")
print("  " + ", ".join(f"{k} {v}" for k, v in sorted(n_moved.items(), key=lambda z: -z[1])))
assert min(n_moved.values()) > 0, "a break row that moves nothing is decoration, not a break row"
nothing = [c for c, v in BR_MOVERS.items() if not v]
print(f"\\nMetrics no row moves ({len(nothing)}), each labelled STRUCT in the coverage map below:")
for c in nothing:
    print(f"  {c}")

# looked up by label, never by list index - inserting a row above must not silently
# repoint these to the wrong variant
BYLABEL = {lab: kw for lab, kw in ROWS}
pick = lambda frag: rerun_all(**next(kw for lab, kw in BYLABEL.items() if frag in lab))  # noqa
o_ld = pick("n(He/H2, present) -4.864"); o_pk = pick("k(H2/CO2, prior) 1200")
print(f"\\nRead the table row by row rather than as a whole. The block swap leaves the crossing")
print(f"selectivity alpha* exactly where it was: alpha* is symmetric in the two lines, so")
print(f"exchanging them cannot move it - a structural check that alpha* is a property of the PAIR")
print(f"of bounds and not of their order. Every Table 12 row corrupts one cell and moves only what")
print(f"that cell feeds; the three Tables 1-11 rows corrupt one point each.")
print(f"\\nThe two LAST-DIGIT rows are the resolving-power rows, and they are why V1's GLOBAL max is")
print(f"the wrong column to read them in (it stays pinned to the unrelated 13a O2/N2 row at")
print(f"{BASE['t13_identity_max_abs_rel']:.2e}). Per row:")
print(f"  n(He/H2) -4.864 -> -4.865 moves the He/H2 identity residual "
      f"{BASE['_v1_heh2']:.2e} -> {o_ld['_v1_heh2']:.2e} ({o_ld['_v1_heh2']/BASE['_v1_heh2']:.0f}x)")
print(f"                            and the V2 cross-check "
      f"{BASE['v2_shift_crosscheck_max_abs_decades']:.2e} -> "
      f"{o_ld['v2_shift_crosscheck_max_abs_decades']:.2e} "
      f"({o_ld['v2_shift_crosscheck_max_abs_decades']/BASE['v2_shift_crosscheck_max_abs_decades']:.1f}x)")
print(f"  k(H2/CO2,prior) 1200 -> 1201 moves that row's identity residual "
      f"{BASE['_v1_h2co2_prior']:.2e} -> {o_pk['_v1_h2co2_prior']:.2e} "
      f"({o_pk['_v1_h2co2_prior']/BASE['_v1_h2co2_prior']:.0f}x)")
print(f"                            and the crossing selectivity "
      f"{BASE['h2co2_crossing_alpha']:.4f} -> {o_pk['h2co2_crossing_alpha']:.4f}")
print(f"So V1 and V2 genuinely PIN the final printed digit of n(He/H2), and the crossing result")
print(f"is sensitive to the final digit of the prior k it depends on - neither is merely")
print(f"consistent with the transcription.")
o_g8, o_dd = pick("n(He/CH4, present) -0.809"), pick("k(He/CH4, present) 19,800")
print(f"\\nThe ranking count is deliberately the coarsest quantity here, and its robustness is worth")
print(f"stating precisely rather than glossing. It survives BOTH last-digit rows (as it should) and")
print(f"it also survives the 8<->3 glyph in n(He/CH4): {BASE['dlogA_top3_hits']:.0f} -> "
      f"{o_g8['dlogA_top3_hits']:.0f}, because a shallower slope")
print(f"only pushes He/CH4 further UP the selectivity-gain ranking. What DOES break it is a lost")
print(f"digit: k(He/CH4) 19,800 -> 1,980 takes the count {BASE['dlogA_top3_hits']:.0f} -> "
      f"{o_dd['dlogA_top3_hits']:.0f}, because He/CH4's shift then")
print(f"reverses sign entirely. So the ranking result rests on the DECADE of every present k, which")
print(f"V1 pins hard, and not on their last digits, which it does not need.")
o_t1, o_t2 = pick("alpha(H2/CO2) 100.9"), pick("P(O2/N2) 18.0")
print(f"\\nAnd the rows that were missing. A single decimal slip in ONE cell of Tables 1-11:")
print(f"  (T1) alpha(H2/CO2) 100.9 -> 10.09 takes pairs_crossing_inside_range "
      f"{BASE['pairs_crossing_inside_range']:.0f} -> {o_t1['pairs_crossing_inside_range']:.0f}")
print(f"       and DESTROYS Result 2 outright: with that point gone from the top of Table 10's")
print(f"       range, alpha* = {BASE['h2co2_crossing_alpha']:.2f} is no longer inside the tabulated "
      f"range and the headline")
print(f"       has nothing left to stand on. No amount of Table 12 checking can see this.")
print(f"  (T2) P(O2/N2) 18.0 -> 1.80 takes points_above_present_bound "
      f"{BASE['points_above_present_bound']:.0f} -> {o_t2['points_above_present_bound']:.0f}, "
      f"moving Result 4's headline,")
print(f"       and the envelope offset {BASE['envelope_offset_max_abs_decades']:.4f} -> "
      f"{o_t2['envelope_offset_max_abs_decades']:.4f} decades, because that point pins O2/N2's")
print(f"       covering line - which is exactly what the paper says it distrusts about it.")
print(f"  (T3) alpha(He/H2) 4.39 -> 4.89 moves envelope_span_min_decades "
      f"{BASE['envelope_span_min_decades']:.4f} -> {pick('alpha(He/H2) 4.39')['envelope_span_min_decades']:.4f},")
print(f"       the lever arm that explains the whole V3 slope discussion.")
print(f"\\nOne more thing (DD) shows, and it is the answer to the obvious objection to the second")
print(f"route. Tables 1-11 are the points Robeson SELECTED as close to his 2008 line, so a covering")
print(f"line fitted to them is biased toward that line and the small baseline residual is an")
print(f"optimistic figure. True - of the BASELINE agreement. It says nothing about the route's")
print(f"DETECTING power, because k_env is untouched by Table 12: d(resid)/d(log k_present) = -1")
print(f"identically. Measured, not argued - (DD) divides k(He/CH4, present) by ten and moves that")
print(f"pair's two-route residual {BASE['_sr_resid_hech4']:.5f} -> {o_dd['_sr_resid_hech4']:.5f}, i.e. by "
      f"{o_dd['_sr_resid_hech4']-BASE['_sr_resid_hech4']:+.5f} decades")
_ddT = BYLABEL[next(l for l in BYLABEL if "k(He/CH4, present) 19,800" in l)]["T12x"]
_ddk = float(_ddT[(_ddT.dataset == "present") & (_ddT.gas_pair == "He/CH4")].k_barrer.iloc[0])
print(f"for a {np.log10(present.k_barrer['He/CH4']) - np.log10(_ddk):.3f}-decade error in the printed "
      f"line. A vertical mis-transcription of the")
print(f"present bound passes through to this residual one for one, whatever the points were")
print(f"selected for.")'''))

cells.append(md(r"""### Blind spots — what this page does *not* establish

1. **Nothing here validates the upper bound against measurement.** Every
   number is Robeson's: his compiled points, his hand-drawn lines, his derived
   tables. The page tests *internal consistency and the meaning of his
   claims*, and that is all. In the gallery's vocabulary this is
   **reproduction, not validation**.
2. **The 117 points are a compilation.** They come from roughly 65 primary
   papers, none of which was consulted. An error made in compiling them is
   inherited silently.
3. **The prior bounds are a restatement.** The 1991 and 1994 papers are not on
   disk. If the 2008 paper misquotes its own earlier $k$ or $n$, every shift
   on this page is wrong by exactly that amount and nothing here would show
   it. V1 touches the prior block for only two of the nine pairs.
4. **$k_s$, and therefore Tables 13a/13b's first column, is unchecked** — the
   critical and Lennard-Jones temperatures behind eqs. (4) and (5) are not
   printed. V1 and V2 take that column as given; had Robeson mis-evaluated
   eq. (4), they would both still pass.
5. **The reference selectivity is this page's choice.** The geometric mean of
   the printed points is defensible and stated, but the front-factor *share*
   is not a property of the data alone, and the page says so in three places
   because it is the single easiest thing to quote out of context. The same
   applies to Result 3's evaluation permeability: the ranking survives every
   reference tried inside the range all nine pairs share, but the *gap*
   between third and fourth does not — it is swept and printed above, with
   both its smallest and its largest value, and neither should be quoted
   without its reference. **And it applies to the second route's own
   headline**, which is the one place an earlier version of this page failed
   to apply it: the residual between the two routes carries $\log\alpha_0$
   explicitly, so it is swept over in-data references the same way, and both
   the residual and the sign agreement are reported at their worst as well as
   at the reference used. Neither number means anything unattached.
6. **Figure 12's slope–diameter correlation and Freeman's eq. (2) are out of
   scope**, for a stated reason: the inputs are not printed in this paper.
   The page reports no verdict on either.
7. **The envelope re-fit (V3) is not a competitor to the printed $(k,n)$** and
   must not be quoted as a "corrected" bound. It answers a different question
   — the tightest covering line — and it is pinned by exactly the extreme
   points the paper distrusts. It is also fitted to points Robeson *selected*
   as close to his 2008 line, so the second route's small **baseline**
   residual is an optimistic figure and should not be read as accuracy. That
   limits the baseline agreement and not the route's **detecting power**:
   $k_{\text{env}}$ is untouched by Table 12, so the residual responds to a
   vertical error in the printed present line exactly one for one — measured
   under break row (DD) below, not argued.
8. **The excursion count is a count against *this* paper's own lines.** It
   says nothing about whether those 19 polymers really beat the state of the
   art; it says the hand-drawn line does not enclose every point the author
   tabulated beside it, which he never claimed it would.
9. **"The paper does not state the crossing" is a searched claim, not an
   assumption.** The full text was searched for cross / crossing /
   intersect / "below the prior"; the only hits are the polymer name
   *crosslinked*, two reference titles containing *crosslinking*, and nothing
   else. What the paper *does* print, and this page quotes, is the slope
   change and the warning that the low-permeability end may have skewed it.
10. **Result 1's decomposition has no second arithmetic route, and cannot
   have one.** Given Table 12 it is a rearrangement, so V2 and the data-side
   shift both reproduce it to machine precision by construction — the page
   measures that rather than claiming otherwise. The one genuine second route
   is V3's envelope, which replaces the *present* line with an estimate from
   Tables 1–11; it agrees on the sign of the shift for every comparable pair,
   and on the magnitude to the tolerance printed in V3 wherever the lever arm
   allows — but it does **not** reproduce the front-factor share at all, which
   is blind spot 5 arriving from the other side. Its scope is the **present
   line only**: nothing on this page can give the *prior* block a second
   route, because both routes read it from the same restatement, and $\alpha_0$
   is likewise common to both. And its corroboration of Result 3's ranking is
   **weaker than the primary route's**: swept over the same evaluation points
   the primary route is held to, the substituted route reproduces the ranking
   at most of them and not at all of them, while the primary route holds at
   every one. The counts are printed above; the asymmetry is the honest
   statement and it is the reported one.
11. **The "inside its own tabulated range" half of Result 2 rests on a single
   tabulated point**, named and quantified in Result 2 and removed by one of
   the break rows below. The crossing itself does not — it is a property of
   the printed $(k,n)$ — but the claim that it lands where the paper has data
   does."""))

# ------------------------------------------------------------------- metrics
cells.append(code('''metrics = dict(
    # ---- structure of the table
    pairs_in_table12=float(N_PAIRS_TABLE),
    pairs_comparable=float(N_PAIRS_COMPARABLE),
    # ---- Result 1: the decomposition
    front_factor_share_median=SHARE_MED,
    front_factor_share_min=SHARE_MIN,
    front_factor_share_max=SHARE_MAX,
    pairs_near_parallel=float(PARALLEL_N),
    heh2_shift_decades=float(SH.shift_gm["He/H2"]),
    heh2_dlogk=float(SH.dlogk["He/H2"]),
    heh2_slope_change_frac=float(abs(SH.dn["He/H2"] / prior.n["He/H2"])),
    heh2_spread_over_shift=float(SH.spread_over_shift["He/H2"]),
    h2ch4_spread_over_shift=float(SH.spread_over_shift["H2/CH4"]),
    h2co2_spread_over_shift=float(SH.spread_over_shift["H2/CO2"]),
    # ---- Result 2: the crossing
    pairs_crossing_inside_range=float(N_CROSSING_INSIDE),
    h2co2_crossing_alpha=ALPHA_STAR,
    h2co2_crossing_alpha_rootfind_rel=ALPHA_STAR_AGREE,
    h2co2_present_over_prior_at_alpha_max_pct=PRESENT_BELOW_PRIOR_PCT,
    # ---- Result 3: which measure reproduces the paper's ranking
    dlogA_top3_hits=float(HITS_A),
    dlogk_top3_hits=float(HITS_K),
    dlogP_top3_hits=float(HITS_P),
    dlogA_gap_third_to_fourth=gap_A,
    dlogA_top3_hits_min_over_eval_points=EVAL_HITS_MIN,
    dlogA_gap_min_over_eval_points=EVAL_GAP_MIN,
    dlogA_rootfind_max_abs=DLOGA_ROOTFIND_MAX,
    alpha_gain_hech4=float(RK.alpha_gain["He/CH4"]),
    alpha_gain_co2ch4=float(RK.alpha_gain["CO2/CH4"]),
    # ---- Result 4: excursions
    points_total=float(N_POINTS),
    points_above_present_bound=float(N_ABOVE),
    frac_above_present_bound=FRAC_ABOVE,
    excursion_max=R_MAX,
    excursion_min=R_MIN,
    excursion_geomean=R_GEOMEAN,
    points_above_prior_bound=float(N_ABOVE_PRIOR),
    # ---- V1
    t13_mw_column_max_abs_rel=MW_MAX_REL,
    t13_mw_column_max_abs_rel_iupac=MW_MAX_REL_IUPAC,
    t13_product_column_max_abs_rel=PROD_MAX_REL,
    t13_identity_max_abs_rel=T13_MAX_REL,
    t13_identity_median_abs_rel=T13_MED_REL,
    t13_rows_inside_tight_band=float(N_INSIDE_TIGHT),
    t13_rows_inside_loose_band=float(N_INSIDE_LOOSE),
    t13a_o2n2_alpha_backout_rel=T13A_O2N2_BACKOUT_REL,
    # ---- V2 (a consistency identity, not a second route - see V2)
    v2_shift_crosscheck_max_abs_decades=V2_MAX_ABS,
    v2_kshare_min=float(V2.k_share.min()),
    v2_kshare_max=float(V2.k_share.max()),
    # ---- V3, and the second route to Result 1 that it carries
    envelope_dn_median_abs=V3_DN_MED,
    envelope_dn_max_abs=V3_DN_MAX,
    envelope_offset_median_decades=V3_OFFSET_MED,
    envelope_offset_max_abs_decades=V3_OFFSET_MAX,
    envelope_offset_max_abs_wide_span=SECOND_ROUTE_WIDE_MAX,
    envelope_offset_max_abs_wide_span_max_over_alpha0=SR_WIDE_MAX_OVER_A0,
    envelope_sign_agree_min_over_alpha0=SR_SIGN_MIN_OVER_A0,
    envelope_dn_max_abs_wide_span=V3_DN_WIDE,
    envelope_span_min_decades=V3_SPAN_MIN,
    envelope_implied_share_median=SHARE_ENV_MED,
    envelope_implied_share_min=float(SR.share_env.min()),
    envelope_implied_share_max=float(SR.share_env.max()),
    envelope_implied_dlogA_top3_hits=float(HITS_ENV),
    envelope_implied_dlogA_top3_hits_min_over_eval_points=HITS_ENV_MIN,
    # ---- structural
    bound_roundtrip_max_rel=ROUNDTRIP_MAX,
)
_ = report_agreement("H1.8", metrics)''' ))

cells.append(code('''# ---- break-row coverage: every metric names the rows MEASURED to move it -------------
# Rows in V4: (S) sign flip; (B) prior/present swap; (L) natural log; (I) inverted excursion;
# (C) comma dropped from k(O2/N2); (G9) 9<->0 in k(He/H2); (G8) 8<->3 in n(He/CH4);
# (DD) digit dropped from k(He/CH4); (D1) last digit of n(He/H2); (D2) last digit of
# k(H2/CO2, prior); (T1) alpha(H2/CO2) in Table 10; (T2) P(O2/N2) in Table 1;
# (T3) alpha(He/H2) in Table 7.
#
# NOTHING in this map is written by hand. Each list is BR_MOVERS[key] - the rows whose
# recomputation of that exact metric differed from the clean one by more than ABS_FLOOR.
# The previous version of this page wrote the lists as prose and asserted only that the
# map's KEYS matched the reported keys; 26 of its 48 metrics carried at least one false
# row-code, 38 false row-codes out of 109 claimed, including
# three that were structurally impossible (the block swap flips numerator and denominator
# of every front-factor share, so it cannot move one; the natural-log row scales both
# alike; alpha* is symmetric in the two lines, so the swap cannot change the crossing
# count). Key-set equality cannot catch any of that. Measured movement can.

# (1) the clean rerun must reproduce every reported metric, by its own code path
_worst, _worst_k = 0.0, None
for k, v in metrics.items():
    d = abs(BASE[k] - v) / max(1.0, abs(v))
    if d > _worst:
        _worst, _worst_k = d, k
_diag = {k for k in BASE if k.startswith("_")}
assert set(BASE) - _diag == set(metrics), (
    f"rerun_all out of step with metrics: {sorted((set(BASE) - _diag) ^ set(metrics))}")
assert _worst < 1e-9, f"rerun_all disagrees with the narrative cells on {_worst_k}: {_worst:.2e}"
print(f"{len(metrics)} reported metrics. The break table's clean row recomputes all of them from")
print(f"the raw tables by a separate code path and agrees with the narrative cells to "
      f"{_worst:.1e} relative")
print(f"({'bit-identical on every key' if _worst_k is None else 'worst: ' + _worst_k}) - so the "
      f"numbers below are assembled twice. The")
print(f"shared primitives (P_bound, alpha_bound, excursion, ulp, sqrtM, envelope_fit) are written")
print(f"once and called by both paths, so this is a check on the assembly, not on them.\\n")

# (2) the structural metrics, each with what it cannot detect stated
STRUCT = {
    "pairs_in_table12": "a count of printed rows. No V4 row can move it; the asserts against "
                        "the Conclusions' pair list and against the set of data tables are what "
                        "guard it. It cannot detect a wrong VALUE anywhere, only a lost row.",
    "pairs_comparable": "as above. Guarded by the assert that Robeson's own classified pair list "
                        "equals this set exactly - which WOULD fail if an NA were lost.",
    "points_total": "row count of Tables 1-11. Guarded by the per-pair counts printed beside it "
                    "and by the assert that the 11 tables match the 11 present bounds. (T1)-(T3) "
                    "corrupt a VALUE in those tables, which is a different defect and moves other "
                    "metrics; nothing here removes a row.",
    "t13_mw_column_max_abs_rel": "no V4 row perturbs a molecular weight or Table 13's printed "
                                 "columns. Its break row is the IUPAC companion printed beside "
                                 "it: standard atomic weights move it by a factor "
                                 f"{MW_MAX_REL_IUPAC/MW_MAX_REL:.0f}, which is what identifies "
                                 "the integer convention rather than assuming it.",
    "t13_mw_column_max_abs_rel_iupac": "the integer-weight value printed beside it is its break "
                                       "row, and vice versa: the pair identifies WHICH weight "
                                       "convention Table 13 used. Neither can detect an error in "
                                       "the convention they agree on - if Robeson had used the "
                                       "same wrong weight in both columns, both values would sit "
                                       "at machine epsilon and say nothing. No V4 row perturbs a "
                                       "molecular weight or a Table 13 column, so it cannot "
                                       "detect anything about Table 12 either.",
    "t13_product_column_max_abs_rel": "a pure two-column arithmetic check on Table 13; no V4 row "
                                      "perturbs ks or the product column, and it cannot detect "
                                      "anything about Table 12. Its above-floor companion is "
                                      "t13_identity_max_abs_rel, which "
                                      f"{len(BR_MOVERS['t13_identity_max_abs_rel'])} V4 rows move.",
    "bound_roundtrip_max_rel": "P_bound and alpha_bound are algebraic inverses, so this is "
                               "machine epsilon by construction and BELOW ABS_FLOOR, i.e. outside "
                               "CI's comparison entirely. Kept because it is what would catch a "
                               "reciprocal-of-n typo in the two functions the whole page is built "
                               "on. Its ABOVE-FLOOR companions are t13_identity_max_abs_rel "
                               "(exercises P_bound on all 26 printed rows) and "
                               "t13a_o2n2_alpha_backout_rel (exercises alpha_bound on a printed "
                               "permeability). dlogA_rootfind_max_abs is NOT a companion: it is "
                               "below the floor too.",
}
# (3) interpretation that measurement cannot supply, kept strictly SEPARATE from the
# measured lists so the two can never be confused again
NOTES = {
    "h2co2_crossing_alpha": "and NOT by (B): alpha* is symmetric in the two lines, so swapping "
                            "them cannot move it. That non-movement is the check.",
    "h2co2_crossing_alpha_rootfind_rel": "BELOW ABS_FLOOR at baseline (exactly 0): closed form "
                                         "against brentq on the same function, so it cannot "
                                         "detect a wrong function - only a wrong inversion of "
                                         "the right one. Above-floor companion: "
                                         "h2co2_crossing_alpha itself.",
    "dlogA_rootfind_max_abs": "BELOW ABS_FLOOR at baseline (machine epsilon): closed form "
                              "against brentq on the same two bounds, so outside CI's "
                              "comparison. Above-floor companion: dlogA_gap_third_to_fourth.",
    "v2_shift_crosscheck_max_abs_decades": f"(L) moves {n_moved['L']} reported metrics but makes only "
                                           "this one FAIL - a residual that should be ~1e-4 "
                                           "decades becomes of order one, while the others "
                                           "merely take different values with nothing to say "
                                           "they are wrong. That is V2's whole independent "
                                           "power, and it is over THIS PAGE's decomposition "
                                           "code, not over Table 12 (see V2).",
    "envelope_offset_max_abs_wide_span": "this is the residual between the two routes to "
                                         "Result 1's shift, on the pairs where the covering line "
                                         "is determined - the second-route headline. AT EACH "
                                         "PAIR'S GEOMETRIC-MEAN SELECTIVITY: it carries log(a_0) "
                                         "explicitly and is meaningless without it - see "
                                         "envelope_offset_max_abs_wide_span_max_over_alpha0. It is "
                                         "computed here as sh_env - shift_gm over the COMPARABLE "
                                         "pairs, the same object the narrative reports; an earlier "
                                         "version measured v3.offset_gm over the point-carrying "
                                         "ones instead, which is equal at the clean point and "
                                         "differs by an order of magnitude under (L).",
    "envelope_offset_max_abs_wide_span_max_over_alpha0": "the same residual at the WORST of the "
                                         "four in-data references, against the metric above at the "
                                         "page's chosen one. The pair is the point: quoting either "
                                         "alone is the error this page is about.",
    "envelope_sign_agree_min_over_alpha0": "the 9-of-9 sign agreement is itself reference-"
                                         "dependent; this is its value at the worst of the four "
                                         "references, and it is NOT 9.",
    "envelope_implied_dlogA_top3_hits": "Result 3 recomputed with a present line that never "
                                        "touches Table 12, AT EACH PAIR'S P_gm. Only the "
                                        "convention and prior-block rows can move it, which is "
                                        "the point: its present half comes from Tables 1-11. Its "
                                        "robustness lives in the _min_over_eval_points companion, "
                                        "which is lower.",
    "envelope_implied_dlogA_top3_hits_min_over_eval_points": "the substituted route under the "
                                        "SAME seven references the primary route is swept over "
                                        "(dlogA_top3_hits_min_over_eval_points). The primary route "
                                        "holds 3 at all seven and this one does not, which is the "
                                        "honest statement of how far the corroboration goes.",
    "envelope_implied_share_median": "the front-factor share under the second route. It does NOT "
                                     "reproduce front_factor_share_median, and that disagreement "
                                     "is Result 1's own caveat measured from the other side.",
    "envelope_implied_share_min": "the NEGATIVE end of that range, and the reason the page "
                                  "explains in words what a negative share is: the front factor "
                                  "moved down while the bound moved up, so the whole shift is "
                                  "carried by the slope term.",
    "envelope_implied_share_max": "the other end. Both ends are metrics because both are quoted; "
                                  "the decomposition route reports its min and max too.",
}
COVERAGE = {k: ("STRUCT: " + STRUCT[k] if k in STRUCT
                else "moved by " + ", ".join(BR_MOVERS[k])
                + (" | " + NOTES[k] if k in NOTES else "")) for k in metrics}
assert set(COVERAGE) == set(metrics)
assert set(NOTES) <= set(metrics) - set(STRUCT)
# the substantive rule: every metric either MOVES under a measured row, or is labelled
# structural with what it cannot detect named. Nothing may be silently unprotected.
for k in metrics:
    assert (BR_MOVERS[k] != []) != (k in STRUCT), (
        f"{k}: movers {BR_MOVERS[k]}, struct={k in STRUCT} - a metric must either have a "
        f"measured mover or be labelled STRUCT, never both and never neither")

W = max(len(k) for k in metrics)
_pad = "\\n" + " " * (W + 4)
for k in metrics:
    print(f"  {k:{W}s}  " + _pad.join(textwrap.wrap(COVERAGE[k], 124 - W)))
ABS_FLOOR = 1e-12
below = {k: v for k, v in metrics.items() if abs(v) < ABS_FLOOR}
print(f"\\nBelow CI's ABS_FLOOR = {ABS_FLOOR:g} and therefore outside the regression comparison: "
      f"{sorted(below) or 'none'}")
for k in below:
    print(f"  {k}: {COVERAGE[k]}")
print(f"\\nStructural metrics (labelled, with what they cannot detect stated): "
      f"{sorted(STRUCT)}")'''))

cells.append(md(r"""## What pymrm adds

**Nothing.** This page uses no pymrm function for anything, and the honest
statement is that the Robeson upper bound is not a pymrm problem: there is no
field to discretise, no operator to assemble, no Jacobian to build. The
`structures` field is empty and `pymrm_api` is empty because that is true, not
because they were not filled in.

What the *page* adds, over reading Table 12:

- It separates the **nine** comparable pairs from the thirteen listed, and
  shows that Robeson's own Conclusions do the same — so a reader who compares
  thirteen is departing from the source, not following it.
- It makes "primarily the front factor" **quantitative and
  reference-explicit**. On the four references *Robeson's own Table 13*
  supplies, the front-factor share of the same two shifts runs from **0.56 to
  1.03**; across the nine pairs at this page's stated reference — the
  geometric-mean selectivity of each pair's own points, which the paper does
  not supply — it runs **0.59 to 1.39**. Both ranges are printed above; the
  distinction between them is the point, and quoting either as "the paper's"
  without saying which is exactly the error the page exists to warn against.
  (An earlier version of this page did precisely that.)
- It shows that $\Delta\log k$ — the quantity the abstract highlights — does
  **not** reproduce the paper's own ranking of which bounds moved most, and
  that the selectivity gain at fixed permeability does, exactly, at every
  evaluation point inside the range all nine pairs share — while reporting
  that the *margin* between third and fourth is itself reference-dependent
  and varies by an order of magnitude across those points. Both of the paper's
  statements survive; they are statements about different things.
- It finds the one pair whose bounds **cross inside its own data range**, and
  says how thin that is: the crossing is a property of four printed numbers,
  but "inside its own range" rests on a single tabulated point.
- It measures what *"by eye"* costs: 19 of 117, and a systematically positive
  envelope offset.
- It says which of its own checks is a **second route** and which is an
  identity — and gets that the other way round from the version of this page
  that was first written. V2 looked like independent corroboration and is
  arithmetically incapable of failing unless V1 does; the real second route is
  the envelope of *Validation V3*, which is the only computation here that does
  not read the present bound out of Table 12. That route corroborates the size
  of the shift — at a stated reference, and by an amount that itself depends on
  that reference, which is why the residual is swept and not quoted at one
  point — and refuses to corroborate its front-factor split, which is the
  page's own central claim arriving from a direction the page did not choose.
  It also corroborates Result 3's ranking **less strongly than the primary
  route does**, under the identical evaluation-point sweep, and the page
  reports that gap rather than the flattering half of it.
- It leaves every printed defect in place and quoted — the doubled `He/CO2`
  in the Conclusions, the `Polypyrrole`/`Polypyrrolone` inconsistency between
  Tables 2 and 5/6, the one Table 13a row outside its own printing band, and
  the paper's own "the ONE data point above the present upper bound" where its
  table has two.

For the transport model that produces a $(P,\alpha)$ pair in the first place,
go to [`H1.7`](../H1.7-solution-diffusion/); that page *is* a pymrm page."""))

# --------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**Use the bound as a sanity ceiling, and quote it correctly.** Given a claimed
$(P_i,\alpha_{ij})$ for a homogeneous polymer film, `excursion(P, alpha, k, n)`
above returns $r$; $r>1$ means the claim is above the 2008 bound. Take the
$(k,n)$ from the `present` block of the Table 12 dataset — and note that only
eleven pairs have one.

**Four cautions, in order of how often they bite.**

1. **Nine, not thirteen.** If you are comparing 1991 with 2008, H2/O2 and
   He/O2 have no 2008 bound and CO2/N2 and N2/CH4 have no 1991 one. Comparing
   them is comparing against `NA`.
2. **Never quote a front-factor share without its reference selectivity.**
   $\Delta\log k$ is the shift at $\alpha=1$, which is outside every dataset in
   the paper. On this page's own numbers the share of the He/H2 shift carried
   by $k$ is 0.98 at the geometric-mean selectivity of its points and 1.02 and
   1.03 at the two selectivities Robeson's own Tables 13a/13b use. For H2/CO2
   it is **0.56 and 0.70 at the two references Table 13 supplies**, **1.39 at
   this page's geometric-mean reference**, and undefined at $\alpha\approx 37$
   where the two lines cross. Name the reference or the number means nothing:
   the same pair, the same two bounds, and a share running from 0.56 to 1.39.
3. **Rank pairs by selectivity gain, not by $k$.** $|n|$ runs from 0.79
   (He/CH4) to 5.8 (O2/N2) in this table, so equal moves in $\log k$ buy
   wildly unequal separation. This is the page's sharpest practical result and
   the easiest thing to get wrong.
4. **The bound applies to homogeneous polymer films only.** The paper excludes
   heterogeneous, surface-modified, mixed-matrix, carbon molecular-sieve and
   thermally-rearranged membranes explicitly, and shows TR polymers above the
   CO2/CH4 bound as a deliberate comparison rather than as a refutation. A
   Pd membrane ([`H1.1`](../H1.1-sieverts-permeation/)) or a zeolite membrane
   ([`H1.9`](../H1.9-zeolite-membrane-maxwell-stefan-mixture/)) exceeding it is
   not news.

**If you extend this page**, the highest-value additions are, in order:
(i) the 1991 paper itself, which would let the prior block be checked rather
than restated and would make every shift on this page a validated number
instead of a reproduced one; (ii) Breck's kinetic diameters, which would open
Fig. 12's $-1/n$ correlation; (iii) solubility constants, which would open
Freeman's eq. (2). None of the three is reachable from this file, and the page
claims nothing about any of them.

**Datasets** are reusable as they stand: `robeson-2008-table12-upper-bounds`
is the canonical $(k,n)$ table, and `robeson-2008-near-bound-points` is 117
literature points with the compiler's citation on every row. Both carry the
caveats above in their sidecars — read them before loading, per the gallery's
cross-page rule."""))

# --------------------------------------------------------------------- write
nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
nbf.write(nb, "index.ipynb")
print(f"wrote index.ipynb with {len(cells)} cells")
