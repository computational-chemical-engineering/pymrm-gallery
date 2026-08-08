#!/usr/bin/env python3
"""Generate index.ipynb for page B2.1 (Voorhies coking law). Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Voorhies' coking law: carbon on catalyst keeps time, not throughput"
description: "The 1945 origin of C_c = A·theta^n, its two printed tables quantified: over a fourfold feed-rate range carbon-on-catalyst moves with feed rate at exponent ~0 while carbon-on-feed moves at ~-1; the printed algebra chain (2)+(3)+(4)->(5) reproduced exactly by symbolic and numeric routes; the 190-200 F doubling claim measured per catalyst; and the diffusion-hypothesis nomenclature K = A^2 shown inconsistent with the paper's own two lines by exactly a factor 2."
categories: [sec:B, struct:S1, tier:T0, data:tier2, phase:gas-solid]
date: 2026-08-08
---

# Voorhies' coking law: carbon on catalyst keeps time, not throughput

**Catalog ID:** `B2.1` · **Structures:** `S1` · **Tier:** T0

Every deactivation model that writes coke as a function of time-on-stream —
including the Froment–Bischoff reactor model of
[`B2.2`](../B2.2-froment-bischoff-coking/), which was written to replace it —
descends from five pages of *Industrial & Engineering Chemistry*, April 1945:

$$
C_c \;=\; A\,\theta^{\,n},
$$

carbon on catalyst as a power of the time since the catalyst was last
regenerated. The paper's abstract stakes a sharper claim than the power law
itself: *"For all the data examined, the amount of carbon deposited on the
catalyst at given conditions is, within limits, independent of the hydrocarbon
feed rate."* Push twice the oil over the same catalyst for the same two hours
and you get the same carbon on the catalyst — the deposit keeps time, not
throughput.

That is the non-obvious claim, and the paper prints the tables that carry it.
This page quantifies what the printed word "independent" amounts to (a
feed-rate exponent, with the competing hypotheses at exponents $+1$ and $0$ an
order of magnitude apart in resolving power), shows that the paper's two
carbon columns are one measurement seen twice — so the claim gets exactly one
test per table, not two — reproduces the printed conversion–feed-rate–time
algebra exactly, measures the temperature-doubling claim per catalyst, and
settles a factor-2 defect in the printed nomenclature of the diffusion
hypothesis from the paper's own two lines."""))

# ----------------------------------------------------------------- background
cells.append(md(r"""## Background

**The source, precisely.** Alexis Voorhies, Jr., "Carbon Formation in
Catalytic Cracking", *Industrial & Engineering Chemistry* **37**(4) 318–322
(April 1945), doi:10.1021/ie50424a010, Esso Laboratories, Standard Oil Company
of New Jersey (Louisiana Division), Baton Rouge, La. — identified from the
scan's own display title, by-line and abstract on page 318, with the ACS
download stamp confirming volume 37, issue 4, first page 318. The by-line
prints **"ALEXIS VOORHIES, JR."**; the catalogue's "Voorhies, A." lacks the
Jr., and this page carries it. The scan is CCITT-G4 bilevel at 300 ppi native;
every numeral on this page was read from cropped native-resolution renders at
digit scale, never from the text layer. (The neighbouring `B1.7` Mears file on
disk opens with the tail of another paper discussing "the Voorhies (1945)
relationship" — a standing reminder in this repository that a first-page
glance is not an identification.)

**What the paper does.** It observes that for fixed-bed and fluid catalytic
cracking alike, carbon-on-catalyst follows $C_c = A\theta^n$ with $n$ in a
narrow band and $A$ specific to catalyst, feed and temperature; that at a
fixed residence time the deposit is nearly independent of feed rate; and it
then chains its two fitted correlations into a working relation between
conversion, feed rate and cracking-period length (its eq. 5) — the practical
payoff, "applied to the commercial-scale units". A closing hypothesis reads
the low temperature coefficient as diffusion control through the deposit
itself, predicting $n = 0.5$.

**Where it sits.** Downstream, Froment & Bischoff (1961, page
[`B2.2`](../B2.2-froment-bischoff-coking/)) argue that once coke lays down as
a *profile* along a bed, an activity correlated only with time is ill-posed —
and their Section 7(a) checks their mechanisms against the very exponent range
this paper prints. `B2.2` established which of their mechanism/activity pairs
can reach the Voorhies exponents (all of them by the consecutive mechanism,
none by parallel-exponential); that finding is `B2.2`'s and is not restated
here — this page owns the original law, its tables, and its own printed
claims. `B2.3` (Levenspiel's deactivation orders) continues the ladder.

**What this page deliberately does not do.** The fixed-bed carbon-vs-time
points (Fig. 1), the conversion–yield points (Fig. 3), the conversion–time
curves (Fig. 4) and the fluid-cracking data (Fig. 5) exist **only as
figures**. Nothing on this page digitises them: the coefficients and exponents
of eqs. (1), (2), (3), (6), (7) are transcribed as printed and are *not*
refit, and every claim tested here is tested against the paper's **printed
tables and printed equations** alone. In particular, Tables I and II hold
$\theta$ fixed at 120 min, so **nothing on this page tests the time exponent
$n$ against data** — the page says so rather than pretending otherwise."""))

# ----------------------------------------------------------- published model
cells.append(md(r"""## The published model

All equation numbers are the paper's; every constant below was read from a
cropped native-resolution render.

**The law** (p. 319, unnumbered, then fitted twice):

$$C_c = A\,\theta^{\,n}$$

with, from the Nomenclature (p. 322): $C_c$ = "carbon, weight per cent on
catalyst", $\theta$ = "catalyst residence time, minutes", $A$ = "constant,
depending on catalyst, feed stock, and temperature", and $n$ = "constant,
depending **only slightly** on catalyst, feed stock, temperature" (emphasis
added — that asymmetry is measured below). For fixed-bed cracking $\theta$ is
the length of the process period since the last regeneration.

**The four fitted instances**, as printed:

| eq. | correlation | system | mode |
|---|---|---|---|
| (1) | $C_c = 0.86\,\theta^{0.41}$ | West Texas gas oil, natural (activated clay) catalyst, 850 °F | fixed bed |
| (2) | $C_c = 0.65\,\theta^{0.44}$ | East Texas gas oil, synthetic (silica–alumina) catalyst, 850 °F | fixed bed |
| (6) | $C_c = 0.52\,\theta^{0.38}$ | East Texas gas oil, natural catalyst, 900–930 °F | fluid |
| (7) | $C_c = 0.24\,\theta^{0.53}$ | Tinsley gas oil, synthetic catalyst, 950 °F | fluid |

(Eq. 2 prints its subscript with a filled counter — the ink reads "$C_.$" —
but its restatements on pp. 321 and 322 print $C_c$ unambiguously.)

**The carbon-yield correlation** (p. 321), for the same runs as eq. (2):

$$C_f = (3.55)(10^{-5})\,V^{2.93} \tag{3}$$

with $C_f$ = carbon, weight per cent **on feed**, and $V$ = volume per cent
conversion of feed stock, defined on p. 318 as 100 minus volume % cycle gas
oil (400 °F i.b.p.) — the paper's own caveat: a correlating convention, "not
an exact measure of feed stock destruction".

**The mass-balance bridge and the derived design equation** (p. 321):

$$C_c = C_f \div (\text{catalyst-to-oil weight ratio})
      = C_f\,U\,\frac{\theta}{60}\,\frac{D_o}{D_c}, \qquad
\frac{D_o}{D_c} = \frac{1}{0.58} \tag{4}$$

$$\therefore\ 0.65\,\theta^{0.44} = (3.55)(10^{-5})\,V^{2.93}\,
   U\,\frac{\theta}{60}\cdot\frac{1}{0.58}
\quad\Longrightarrow\quad
V = 96\cdot\frac{1}{U^{0.34}}\cdot\frac{1}{\theta^{0.19}} \tag{5}$$

where $U$ is feed rate in v/v/hr and $W = U\,D_o/D_c$ (p. 319) converts it to
a weight basis. The chain (2)+(3)+(4)→(5) is pure printed algebra and is
reproduced exactly below, twice — by a symbolic route and a numeric route with
no solver code in common (they do share the one hand isolation step; the
printed 96 is the external check on that step).

**The diffusion hypothesis** (p. 322): if the rate of diffusion through the
deposit is inversely proportional to the deposit,

$$\frac{dC_c}{d\theta} = \frac{K}{C_c}
\quad\text{or}\quad
C_c = A\,\theta^{0.5} \tag{8}$$

and the Nomenclature prints "$K$ = constant $(= A^2)$". Integrating the
paper's own left-hand line gives $K = A^2/2$ — a factor-2 defect settled
below, from these two lines alone, and **reported, not repaired**.

**Two studies, kept apart.** The paper's fixed-bed material comes from two
distinct experimental campaigns, and conflating them manufactures a
contradiction: (i) the **Fig.-1 study** behind eqs. (2), (3), (5) — small
fixed-bed units, *one* repeatedly regenerated batch of synthetic catalyst of
"such inherent stability that its activity was completely restored after each
successive regeneration", carbon determined from the regeneration gases; and
(ii) the **Table I/II study** — a small-scale fixed-bed unit, a *fresh batch
of catalyst every run*, carbon by combustion of the discharged catalyst.
Same oil, same nominal 850 °F, different units and batches; the measured gap
between their $A$ values is one of the results below."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

The law is two constants per system; nothing else enters. Everything this page
computes uses only the printed values:

| quantity | value | provenance |
|---|---|---|
| $(A, n)$ eq. (1) | 0.86, 0.41 | p. 319, transcribed (fit to figure-only data; not refit here) |
| $(A, n)$ eq. (2) | 0.65, 0.44 | p. 319, transcribed (fit to figure-only data; not refit here) |
| eq. (3) | $3.55\times10^{-5}$, 2.93 | p. 321, transcribed (fit to figure-only data; not refit here) |
| $D_o/D_c$ | $1/0.58$ | p. 321 — printed once, inside the eq.-(5) derivation only |
| $(A, n)$ eqs. (6), (7) | 0.52, 0.38 / 0.24, 0.53 | p. 322, transcribed (fluid; fits to figure-only data) |
| eq. (5) triple | 96, 0.34, 0.19 | p. 321 — printed algebra, re-derived exactly below |
| $\theta$ for Tables I–II | 120 min | "uniformly maintained at 2 hours", p. 320 |

**Conventions.** Feed-rate exponents are $m \equiv d\ln C/d\ln U$ at fixed
$\theta$ (and nominally fixed $T$), by least squares on the printed rows;
deviations are (value − printed)/printed. Temperatures stay in °F as printed.

**Assumptions carried from the source.** Carbon means *catalytic* carbon —
strippable carbon is purged before determination (p. 319) and the paper calls
its contribution negligible for these operations. Conversion is the cycle-gas-oil
convention above. Substantially atmospheric pressure ("pressure was not a
variable in the data presented here", p. 319). The feed-rate-independence
claim is bounded by the paper itself: *"It is not proved or claimed that this
observation will hold for all feed stocks, all catalysts, or an extreme
variation in feed rates"* (p. 320) — the tables span 2× (Table I) and 4×
(Table II) in feed rate, and every statement below lives inside that span."""))

# ------------------------------------------------------------ environment
cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code('''import sys, urllib.request
from pathlib import Path

if not any("shared" in p for p in sys.path):
    local = Path.cwd()
    for _ in range(5):
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
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "B2.1-voorhies-coking-law"
PAGE_B22 = "B2.2-froment-bischoff-coking"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

THETA_TAB = 120.0        # min; Tables I and II, "uniformly maintained at 2 hours"'''))

# --------------------------------------------------------------------- data
cells.append(md(r"""## The data

Three datasets, all transcriptions from this paper — read from cropped
native-resolution renders, never the text layer.

### 1. Table I (p. 320) — East Texas gas oil, 32 runs, the feed-rate and temperature evidence

Four feed rates × four nominal temperature levels × two catalysts, every run a
fresh 2-hour period on a fresh catalyst batch, carbon by combustion of the
discharged catalyst. This is the paper's own supporting evidence for the
feed-rate claim and the data behind its Fig.-2 temperature claim — so every
test against it below is a **consistency check of the paper's claims with the
paper's own data**, and is labelled so; nothing here is held-out validation.

Two transcription notes, both settled at digit scale: the natural-catalyst
793 °F row and the synthetic 911 °F row print feed rate **1.1** (not 1.0; the
synthetic one was verified twice at 4× magnification after a whole-page read
mis-gave it). The mass-balance identity of Check 1 mildly disfavours the
synthetic 1.1 — its implied catalyst-to-oil ratio sits ~6 % below its block
mates, where 1.0 would leave it ~3 % high — but the glyph is unambiguous, so
1.1 is transcribed and the tension is reported here rather than repaired.

### 2. Table II (p. 320) — cetane and Decalin, the fourfold feed-rate range

Pure-compound cracking, natural catalyst, 850 °F, 2-hour periods, feed rates
0.6–2.4 v/v/hr — the widest feed-rate span in the paper, and therefore the
sharpest printed test of the claim.

### 3. The printed constants and prose claims

Every coefficient, exponent and quantitative sentence this page uses, with
locations, in `voorhies-1945-printed-constants.csv`. The five fitted $(A,n)$
pairs and eq. (3) are the author's fits to **figure-only** data; they enter
this page as transcriptions and are never refit.

**Cross-page reconciliation.** `B2.2`'s dataset
`froment-bischoff-1961-printed-claims.csv` carries the Voorhies exponent range
as Froment & Bischoff quote it (their Sec. 7(a)); since this page states the
same range from the original, the two are printed side by side below rather
than retyped. `B2.2`'s findings about those rows — that they are
transcriptions whose OCR mangles the range ("0.88 to 058" for 0.38 to 0.53),
and that the quoted range is other papers' measurements seen through F&B —
are inherited and none affects their use here as a cross-reference."""))

cells.append(code('''t1 = load_data("voorhies-1945-tableI.csv", page=PAGE)
t2 = load_data("voorhies-1945-tableII.csv", page=PAGE)
pc = load_data("voorhies-1945-printed-constants.csv", page=PAGE)
t1_meta = load_meta("voorhies-1945-tableI.csv", page=PAGE)
pcv = dict(zip(pc.item_id, pc.printed_value.astype(float)))

print("Table I (p. 320): fixed-bed, East Texas gas oil, 2-hour periods")
print(t1.to_string(index=False))
print(f"\\n{cite_data(t1_meta)}")
print("\\nTable II (p. 320): cetane and Decalin, natural catalyst, 850 F, 2-hour periods")
print(t2.to_string(index=False))
print("\\nPrinted constants and claims used by this page:")
print(pc[["item_id", "printed_value", "where"]].to_string(index=False))'''))

cells.append(code('''# --- cross-page reconciliation: the exponent range B2.2 quotes vs the original -
fb = load_data("froment-bischoff-1961-printed-claims.csv", page=PAGE_B22)
fb_lo = float(fb.loc[fb.claim_id == "voorhies_n_low",  "printed_value"].iloc[0])
fb_hi = float(fb.loc[fb.claim_id == "voorhies_n_high", "printed_value"].iloc[0])
n_printed = {k: pcv[k] for k in ("eq1_n", "eq2_n", "eq6_n", "eq7_n")}
print("The four gas-oil exponents printed by Voorhies (1945):",
      ", ".join(f"{k[:3]} n = {v:g}" for k, v in n_printed.items()))
print(f"min/max of the four: {min(n_printed.values()):g} / {max(n_printed.values()):g}")
print(f"Froment & Bischoff (1961) Sec. 7(a), as transcribed on B2.2's page: "
      f"'{fb_lo:g} to {fb_hi:g} for gas-oil cracking (Voorhies)'")
same = (min(n_printed.values()) == fb_lo) and (max(n_printed.values()) == fb_hi)
print(f"-> F&B's quoted range is exactly the min/max of the four printed equations: {same}")
assert same, "reconciliation with B2.2's transcription failed - investigate before publishing"'''))

# ------------------------------------------------------ pymrm implementation
cells.append(md(r"""## PyMRM implementation

**There is no field to discretise on this page and it would be dishonest to
invent one.** The law is a two-parameter power function; the feed-rate test is
least squares on log-transformed printed rows; the eq.-(5) chain is algebra;
the one differential equation in the paper (eq. 8) is a scalar separable ODE,
integrated symbolically and — as its independent second route — numerically
with `solve_ivp`. No `construct_grad`, no `construct_div`, no Newton solve,
and nothing here would run differently with pymrm uninstalled. This page
follows [`A1.6`](../A1.6-wen-yu-minimum-fluidisation/),
[`A1.1`](../A1.1-ergun-pressure-drop/) and
[`F1.4`](../F1.4-krishna-ellenberger-holdup/) in saying so in the section
where a reader would otherwise expect a solver. (The bed that *does* need
operators — coke laying down as a profile along a reactor — is the sequel,
[`B2.2`](../B2.2-froment-bischoff-coking/), which builds exactly that on this
law.)

Below: the law, the two slope estimators, and the eq.-(5) closures."""))

cells.append(code('''def voorhies(theta, A, n):
    """C_c = A theta^n, theta in minutes (eq. 1/2/6/7/8 with the printed pairs)."""
    return A * np.asarray(theta, float) ** n


def slope_ols(U, C):
    """m = d ln C / d ln U by least squares - the page's feed-rate exponent."""
    return float(np.polyfit(np.log(np.asarray(U, float)),
                            np.log(np.asarray(C, float)), 1)[0])


def slope_endpoints(U, C):
    """The same exponent from the two extreme feed rates only - a second,
    regression-free estimator (used to show the verdict is not an OLS artifact)."""
    U, C = np.asarray(U, float), np.asarray(C, float)
    i, j = int(np.argmin(U)), int(np.argmax(U))
    return float(np.log(C[j] / C[i]) / np.log(U[j] / U[i]))


def slope_rounding_floor(U, C, half_ulp=0.05):
    """Worst-case |shift| of the OLS slope if every printed C is off by up to
    half a last printed digit (0.05 wt%): linearised, delta_m = sum w_i dlnC_i
    with the OLS weights w_i, maximised by the adversarial sign pattern."""
    x = np.log(np.asarray(U, float))
    w = (x - x.mean()) / np.sum((x - x.mean()) ** 2)
    return float(np.sum(np.abs(w) * (half_ulp / np.asarray(C, float))))


def eq5_V(U, theta, coeff=96.0, eU=0.34, eT=0.19):
    """Eq. (5) as printed: V = coeff / (U^eU theta^eT)."""
    return coeff / (np.asarray(U, float) ** eU * np.asarray(theta, float) ** eT)


def eq5_V_implicit(U, theta, A2=0.65, n2=0.44, c3=3.55e-5, e3=2.93, dr=0.58):
    """V from the chain (2)+(3)+(4), solved numerically: bisection on
       V^e3 = A2 theta^n2 * 60 dr / (c3 U theta), then no further algebra.
    Shares NO SOLVER with the symbolic route below, but it does rest on the
    same hand isolation of V^e3 - an error in that step would fool both routes
    alike, and is instead caught by the comparison against the printed 96."""
    rhs = A2 * theta ** n2 * 60.0 * dr / (c3 * U * theta)
    return brentq(lambda V: V ** e3 - rhs, 1e-3, 1e4, xtol=1e-12, rtol=1e-14)'''))

# ------------------------------------------------------------------ results
cells.append(md(r"""## Results

### 1. What "independent of feed rate" measures as

For each temperature block of Table I and each pure compound of Table II, the
feed-rate exponent $m = d\ln C/d\ln U$ of **both** carbon columns. The two
competing bookkeepings sit an exponent apart:

- if carbon tracked **throughput** — every unit of feed depositing its share —
  $C_c$ at fixed $\theta$ would scale $\propto U$: $m_{cat} = +1$;
- if carbon keeps **time**, $m_{cat} = 0$, and the mass balance then forces
  $m_{feed} = -1$ (fewer per cent of a larger throughput).

The rounding floor beside each slope is the worst-case shift a half-digit
(±0.05 wt %) of printed rounding could produce — the resolution limit of a
2-significant-figure table."""))

cells.append(code('''rows = []
for (cat, blk), g in t1.groupby(["catalyst", "block"]):
    rows.append(dict(dataset="Table I", group=f"{cat} {g.temp_F.mean():.0f}F",
                     m_feed=slope_ols(g.feed_rate_vvhr, g.C_feed_wtpct),
                     m_cat=slope_ols(g.feed_rate_vvhr, g.C_cat_wtpct),
                     floor=slope_rounding_floor(g.feed_rate_vvhr, g.C_cat_wtpct),
                     span=g.feed_rate_vvhr.max() / g.feed_rate_vvhr.min()))
for feed, g in t2.groupby("feed"):
    rows.append(dict(dataset="Table II", group=feed,
                     m_feed=slope_ols(g.feed_rate_vvhr, g.C_feed_wtpct),
                     m_cat=slope_ols(g.feed_rate_vvhr, g.C_cat_wtpct),
                     floor=slope_rounding_floor(g.feed_rate_vvhr, g.C_cat_wtpct),
                     span=g.feed_rate_vvhr.max() / g.feed_rate_vvhr.min()))
sl = pd.DataFrame(rows)
print("feed-rate exponents m = dlnC/dlnU at fixed theta = 120 min")
print(f"{'group':22s} {'span':>5s} {'m_feed':>8s} {'m_cat':>8s} {'floor':>7s}")
for r in sl.itertuples():
    print(f"{r.group:22s} {r.span:4.0f}x {r.m_feed:+8.3f} {r.m_cat:+8.3f} ±{r.floor:6.3f}")

M_FEED_T1 = float(sl[sl.dataset == "Table I"].m_feed.mean())
M_CAT_T1 = float(sl[sl.dataset == "Table I"].m_cat.mean())
M_CAT_T1_WORST = float(sl[sl.dataset == "Table I"].m_cat.abs().max()
                       * np.sign(sl[sl.dataset == "Table I"].m_cat[
                           sl[sl.dataset == "Table I"].m_cat.abs().idxmax()]))
FLOOR_T1 = float(sl[sl.dataset == "Table I"].floor.max())
M_CAT_CET = float(sl.loc[sl.group == "cetane", "m_cat"].iloc[0])
M_FEED_CET = float(sl.loc[sl.group == "cetane", "m_feed"].iloc[0])
M_CAT_DEC = float(sl.loc[sl.group == "decalin", "m_cat"].iloc[0])
M_FEED_DEC = float(sl.loc[sl.group == "decalin", "m_feed"].iloc[0])
FLOOR_T2 = float(sl[sl.dataset == "Table II"].floor.max())
M_CAT_CET_EP = slope_endpoints(t2[t2.feed == "cetane"].feed_rate_vvhr,
                               t2[t2.feed == "cetane"].C_cat_wtpct)

print(f"\\nTable I  (2x span): m_feed mean {M_FEED_T1:+.3f}, m_cat mean {M_CAT_T1:+.3f}, "
      f"worst single block {M_CAT_T1_WORST:+.3f}, rounding floor up to ±{FLOOR_T1:.3f}")
print(f"Table II (4x span): cetane {M_FEED_CET:+.3f} / {M_CAT_CET:+.3f}, "
      f"decalin {M_FEED_DEC:+.3f} / {M_CAT_DEC:+.3f}, floor ±{FLOOR_T2:.3f}")
print(f"cetane m_cat by the endpoint estimator (no regression): {M_CAT_CET_EP:+.3f}")
CLOSER = float(((1 - sl.m_cat).abs() / sl.m_cat.abs()).min())
FEED_DEV_LO = float((sl.m_feed + 1).abs().min())
FEED_DEV_HI = float((sl.m_feed + 1).abs().max())
print(f"\\nReading: against the two hypotheses m_cat = +1 (throughput) and 0 (time),")
print(f"every group sits at |m_cat| <= {abs(M_CAT_T1_WORST):.2f} - at least {CLOSER:.1f}x closer to 0")
print(f"than to +1 in every group - while the")
print(f"SAME rows put m_feed within {100*FEED_DEV_LO:.0f}-{100*FEED_DEV_HI:.0f} % of -1. "
      f"Over Table II's fourfold range,")
print(f"C_cat moves {t2[t2.feed=='cetane'].C_cat_wtpct.min():.1f}-"
      f"{t2[t2.feed=='cetane'].C_cat_wtpct.max():.1f} wt% (cetane) while C_feed falls "
      f"{t2[t2.feed=='cetane'].C_feed_wtpct.max():.1f} -> "
      f"{t2[t2.feed=='cetane'].C_feed_wtpct.min():.1f} - a factor "
      f"{t2[t2.feed=='cetane'].C_feed_wtpct.max()/t2[t2.feed=='cetane'].C_feed_wtpct.min():.1f}.")'''))

cells.append(code('''# --- the temperature confound inside Table I blocks, bounded with the paper's own rule
# Within a block the temperature is NOT constant (spreads up to ~22 F) and is
# correlated with feed rate by accident of operation; at the paper's doubling
# rate that alone can push |m_cat| by up to ~0.1. Correct each C_cat to its
# block-mean temperature with the doubling rule fitted in section 4 below
# (2^(dT/DT2), catalyst's own DT2), and re-read the slopes.
def doubling_dT(temp_F, C):
    s = np.polyfit(np.asarray(temp_F, float), np.log(np.asarray(C, float)), 1)[0]
    return float(np.log(2.0) / s)

DT2 = {cat: doubling_dT(g.temp_F, g.C_cat_wtpct) for cat, g in t1.groupby("catalyst")}

def tcorr_slopes(sign=+1.0):
    out = []
    for (cat, blk), g in t1.groupby(["catalyst", "block"]):
        f = 2.0 ** (sign * (g.temp_F.mean() - g.temp_F) / DT2[cat])
        out.append(slope_ols(g.feed_rate_vvhr, g.C_cat_wtpct * f))
    return np.array(out)

m_raw = sl[sl.dataset == "Table I"].m_cat.to_numpy()
m_cor = tcorr_slopes(+1.0)
M_CAT_T1_TCORR = float(m_cor.mean())
print("Table I m_cat, raw vs corrected to block-mean T with the paper's own doubling rule:")
for r, a, b in zip(sl[sl.dataset == "Table I"].itertuples(), m_raw, m_cor):
    print(f"  {r.group:18s} raw {a:+.3f} -> corrected {b:+.3f}")
print(f"  mean {m_raw.mean():+.3f} -> {M_CAT_T1_TCORR:+.3f}")
print("The correction moves individual blocks by up to "
      f"{np.abs(m_cor-m_raw).max():.3f} but the verdict (|m_cat| << 1) survives it;")
print("the sign-flipped correction is a break-table row below.")'''))

cells.append(code('''# --- figure: both carbon columns against feed rate, log-log ------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=False)
blues = plt.cm.Blues(np.linspace(0.45, 0.95, 4))
for (cat, ls, panel) in [("natural", "-", 0), ("synthetic", "--", 0)]:
    for b in range(1, 5):
        g = t1[(t1.catalyst == cat) & (t1.block == b)]
        lbl = f"{cat} {g.temp_F.mean():.0f}F" if True else None
        ax[0].loglog(g.feed_rate_vvhr, g.C_cat_wtpct, "o" + ls, color=blues[b - 1],
                     mfc=("white" if cat == "synthetic" else None), ms=6, lw=1.2)
        ax[0].loglog(g.feed_rate_vvhr, g.C_feed_wtpct, "s" + ls, color=blues[b - 1],
                     mfc=("white" if cat == "synthetic" else None), ms=4, lw=0.8, alpha=0.55)
ax[0].set_xlabel("feed rate U, v/v/hr"); ax[0].set_ylabel("carbon, wt %")
ax[0].set_title("Table I: $C_{cat}$ (circles, flat) vs $C_{feed}$ (squares, ~1/U)\\n"
                "blues light→dark = 805→958 °F; filled natural, open synthetic")
for feed, c in [("cetane", "tab:blue"), ("decalin", "tab:orange")]:
    g = t2[t2.feed == feed]
    ax[1].loglog(g.feed_rate_vvhr, g.C_cat_wtpct, "o-", color=c, ms=7, label=f"{feed}, on catalyst")
    ax[1].loglog(g.feed_rate_vvhr, g.C_feed_wtpct, "s--", color=c, ms=5, alpha=0.6,
                 label=f"{feed}, on feed")
uu = np.array([0.6, 2.4])
ax[1].loglog(uu, 5.1 * (uu / 0.6) ** 1.0, ":", color="0.4", lw=1)
ax[1].text(1.5, 11.0, "slope +1 (throughput\\nhypothesis)", fontsize=8, color="0.35")
ax[1].loglog(uu, 4.8 * 0.6 / uu, ":", color="0.4", lw=1)
ax[1].text(1.6, 1.05, "slope −1", fontsize=8, color="0.35")
ax[1].set_xlabel("feed rate U, v/v/hr"); ax[1].set_ylabel("carbon, wt %")
ax[1].set_title("Table II (4× span): the claim at its sharpest")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""### 2. The two carbon columns are one measurement — so each table carries *one* test

The paper states how Table I was measured: *"After each run the catalyst was
discharged from the unit and analyzed for carbon by combustion"* — one carbon
determination per run. The two printed columns are then related by the paper's
own eq.-(4) bookkeeping: $C_c = C_f \cdot W\theta/60$ with $W = U\,D_o/D_c$,
so the ratio

$$\frac{W}{U}\Big|_{implied} = \frac{C_c}{C_f\cdot U\cdot(\theta/60)}
   \;\left(= \frac{D_o}{D_c}\right)$$

should be a **constant of the unit and charge, row by row** — and if it is
constant, the two columns carry the same information and $m_{cat}\approx0$
*is* $m_{feed}\approx-1$ restated. Checking its constancy simultaneously (a)
tests that reading, (b) cross-checks all four transcribed columns at once (a
mis-read digit breaks its own row), and (c) recovers the number the paper
never prints for these units — the density ratio the eq.-(5) derivation would
have needed *here* (searched: the only density ratio printed anywhere in the
five pages is the 1/0.58 on p. 321, inside the fluid-oriented eq.-(5)
derivation; the Nomenclature defines $D_o$, $D_c$ without values)."""))

cells.append(code('''def implied_WU(df):
    return df.C_cat_wtpct / (df.C_feed_wtpct * df.feed_rate_vvhr * (THETA_TAB / 60.0))

print(f"{'group':22s} {'mean':>7s} {'CV %':>6s}  rows")
WU = {}
for (cat, blk), g in t1.groupby(["catalyst", "block"]):
    r = implied_WU(g); WU[f"{cat}{blk}"] = r
    print(f"Table I {cat} {g.temp_F.mean():4.0f}F   {r.mean():7.3f} {100*r.std(ddof=0)/r.mean():6.1f}  "
          + " ".join(f"{v:.3f}" for v in r))
for feed, g in t2.groupby("feed"):
    r = implied_WU(g); WU[feed] = r
    print(f"Table II {feed:13s}  {r.mean():7.3f} {100*r.std(ddof=0)/r.mean():6.1f}  "
          + " ".join(f"{v:.3f}" for v in r))

WU_NAT = float(implied_WU(t1[t1.catalyst == "natural"]).mean())
WU_SYN = float(implied_WU(t1[t1.catalyst == "synthetic"]).mean())
WU_CET = float(WU["cetane"].mean())
WU_DEC = float(WU["decalin"].mean())
CV_MAX = max(float(100 * v.std(ddof=0) / v.mean()) for v in WU.values())
WU_RATIO = WU_DEC / WU_CET

# what half-a-last-printed-digit rounding can do to this ratio, per row. Each
# cell's allowance comes from its OWN printed decimals, the sidecars' stated
# convention: +-0.05 on a 1-decimal cell like 4.7, but +-0.005 on Table II's
# 2-decimal 0.75. (Re-reading the CSVs as strings recovers those decimals.)
t1_str = load_data("voorhies-1945-tableI.csv", page=PAGE, dtype=str)
t2_str = load_data("voorhies-1945-tableII.csv", page=PAGE, dtype=str)

def half_ulp(printed):
    return np.array([0.5 * 10.0 ** -(len(s.split(".")[1]) if "." in s else 0)
                     for s in printed])

def row_allowance_pct(df, dfs):
    uC, uF = half_ulp(dfs.C_cat_wtpct), half_ulp(dfs.C_feed_wtpct)
    return ((df.C_cat_wtpct + uC) / (df.C_cat_wtpct - uC)
            * (df.C_feed_wtpct + uF) / (df.C_feed_wtpct - uF) - 1).to_numpy() * 100 / 2

allow_rows = np.concatenate([row_allowance_pct(t1, t1_str), row_allowance_pct(t2, t2_str)])
allow_labels = pd.concat([t1.catalyst + " block " + t1.block.astype(str)
                          + ", U = " + t1.feed_rate_vvhr.astype(str),
                          t2.feed + ", U = " + t2.feed_rate_vvhr.astype(str)],
                         ignore_index=True)
RMAX = float(allow_rows.max())
RMAX_ROW = allow_labels[int(np.argmax(allow_rows))]
print(f"\\nWithin-group scatter: worst CV {CV_MAX:.1f} %, against a per-row half-digit")
print(f"allowance of up to {RMAX:.1f} % ({RMAX_ROW}, whose small 1-decimal carbon")
print(f"cells dominate; the 2-decimal 0.75 cell gets only +-0.005) -")
print("the constancy is at the rounding floor. Consistent with (and at this precision")
print("indistinguishable from) both columns being computed from ONE carbon")
print("determination through the run's fixed feed-to-catalyst mass ratio; two")
print("independent balances agreeing this well on all 40 rows is not excluded, but")
print("either way the columns are redundant AT THE PRINTED PRECISION, so the")
print("feed-rate claim gets ONE test per group, not two.")
print(f"\\nImplied W/U (= D_o/D_c of eq. 4): natural {WU_NAT:.3f}, synthetic {WU_SYN:.3f}, "
      f"cetane {WU_CET:.3f}, decalin {WU_DEC:.3f}")
print(f"  decalin/cetane ratio {WU_RATIO:.3f} - same unit, same catalyst, so W = U D_o/D_c")
print("  attributes the whole difference to the two liquids' density ratio; the paper")
print("  prints no densities, so this is a consistency observation, not a check")
print("  against known values.")
print(f"  None of the four is the 1/0.58 = {1/0.58:.3f} printed for the eq.-(5) unit -")
print("  that ratio belongs to that study and would be wrong here.")

# the one systematic residual, reported not explained:
bm = {(cat, b): float(implied_WU(g).mean()) for (cat, b), g in t1.groupby(["catalyst", "block"])}
DRIFT_NAT = 100 * (bm[("natural", 1)] / bm[("natural", 4)] - 1)
DRIFT_SYN = 100 * (bm[("synthetic", 1)] / bm[("synthetic", 4)] - 1)
print(f"\\nUnexplained residual: the implied ratio drifts DOWN across the temperature")
print(f"blocks - natural {bm[('natural',1)]:.3f} -> {bm[('natural',4)]:.3f}, synthetic "
      f"{bm[('synthetic',1)]:.3f} -> {bm[('synthetic',4)]:.3f} from ~805 to ~958 F")
print(f"(block-1/block-4 ratio - 1: {DRIFT_NAT:+.1f} % and {DRIFT_SYN:+.1f} %).")
print("A charging mass ratio has no business varying with temperature; the drift is")
print("larger than the within-block scatter, and nothing in the paper explains it.")
print("It does not touch the feed-rate verdict, which lives within blocks.")'''))

cells.append(md(r"""### 3. The printed algebra chain, (2)+(3)+(4) → (5), exactly — by two routes

The paper substitutes eq. (2) and eq. (3) into the mass balance (4) with
$D_o/D_c = 1/0.58$ and prints the result as eq. (5): $V = 96/(U^{0.34}\theta^{0.19})$.
The chain is closed-form, so it can be checked to machine precision — and it
is one of this page's two headlines computed twice, with **independent
solvers**: symbolically (exact rationals, sympy), and numerically (the chain is
*root-found* for $V$ by bisection, point by point on a $(U,\theta)$ grid, and
the printed power-law form is recovered by regression on that grid — no
computer algebra involved). One step *is* shared: both routes start from the
same hand isolation of $V^{2.93}$ from (2)+(3)+(4), so their agreement checks
the solvers, not that step. The check on the shared step is external — an
error in it would move the derived coefficient away from the *printed* 96,
which the comparison against the printed triple catches."""))

cells.append(code('''# --- route A: symbolic ------------------------------------------------------
# (exact rationals throughout: 3.55e-5 = 355/10^7, 2.93 = 293/100, so the
# isolation of V is exact and the exponents come out as closed-form rationals)
U_, th_ = sp.symbols("U theta", positive=True)
A2, n2 = sp.Rational(65, 100), sp.Rational(44, 100)
c3, e3, dr = sp.Rational(355, 10 ** 7), sp.Rational(293, 100), sp.Rational(58, 100)
# eq. (2) = eq. (3)*(4):  A2 th^n2 = c3 V^e3 U (th/60) / dr, isolated for V:
Vsol = (A2 * th_ ** n2 * 60 * dr / (c3 * U_ * th_)) ** (1 / e3)
# substituting back is an identity by construction here; route B repeats the
# same isolation but solves it with a different tool chain (bisection +
# regression, no computer algebra), so A-vs-B agreement checks the solvers.
# The isolation step itself is checked against the printed 96 (Check 2).
COEFF_SYM = float(Vsol.subs({U_: 1, th_: 1}))
EXP_U_SYM = -float(sp.simplify(sp.diff(sp.log(Vsol), U_) * U_))
EXP_T_SYM = -float(sp.simplify(sp.diff(sp.log(Vsol), th_) * th_))

# --- route B: numeric root-find + regression (no sympy anywhere) ------------
Ug, Tg = np.meshgrid(np.linspace(0.3, 1.2, 7), np.linspace(30.0, 240.0, 8))
Vg = np.array([eq5_V_implicit(u, t) for u, t in zip(Ug.ravel(), Tg.ravel())])
X = np.column_stack([np.ones_like(Vg), np.log(Ug.ravel()), np.log(Tg.ravel())])
beta, *_ = np.linalg.lstsq(X, np.log(Vg), rcond=None)
COEFF_NUM, EXP_U_NUM, EXP_T_NUM = float(np.exp(beta[0])), -float(beta[1]), -float(beta[2])
resid = float(np.abs(X @ beta - np.log(Vg)).max())

print("eq. (5) re-derived from eqs. (2)+(3)+(4):")
print(f"  route A (symbolic):            V = {COEFF_SYM:.4f} / (U^{EXP_U_SYM:.5f} theta^{EXP_T_SYM:.5f})")
print(f"  route B (root-find + lstsq):   V = {COEFF_NUM:.4f} / (U^{EXP_U_NUM:.5f} theta^{EXP_T_NUM:.5f})")
print(f"    (max |log-residual| of the power-law fit on the grid: {resid:.2e} -")
print(f"     the chain IS exactly a power law, which route B verifies rather than assumes)")
print(f"  printed:                       V = {pcv['eq5_coeff']:g} / (U^{pcv['eq5_expU']:g} theta^{pcv['eq5_expTheta']:g})")
print(f"  -> every printed digit is the rounding of the exact chain: "
      f"{COEFF_SYM:.4f} -> {COEFF_SYM:.0f}, {EXP_U_SYM:.5f} -> {EXP_U_SYM:.2f}, "
      f"{EXP_T_SYM:.5f} -> {EXP_T_SYM:.2f}")

# exact exponent identities, stated:
print(f"\\n  identities: 1/2.93 = {1/2.93:.5f} (= U exponent), "
      f"(1-0.44)/2.93 = {(1-0.44)/2.93:.5f} (= theta exponent)")

# what the printed rounding costs across the paper's own stated ranges:
VV_print = eq5_V(Ug, Tg)
VV_exact = eq5_V(Ug, Tg, COEFF_SYM, EXP_U_SYM, EXP_T_SYM)
ROUND_COST = float(np.abs(VV_print / VV_exact - 1).max()) * 100
print(f"\\n  printed (96, 0.34, 0.19) vs exact chain over U = 0.3-1.2, theta = 30-240:")
print(f"  worst deviation {ROUND_COST:.2f} % - the printed rounding is harmless there.")

# a small exact object the paper does not print:
INST_OVER_AVG = 1.0 - EXP_T_SYM
print(f"\\n  Since V(theta) is the PERIOD-average conversion, the instantaneous")
print(f"  conversion the chain implies is d(theta V)/dtheta = (1 - {EXP_T_SYM:.4f}) V")
print(f"  = {INST_OVER_AVG:.4f} V: at every moment the instantaneous conversion sits at")
print(f"  {100*INST_OVER_AVG:.1f} % of the period average to date - declining, never negative,")
print(f"  so the derived correlation is self-consistent as an averaged quantity.")'''))

cells.append(code('''# --- the "faulty point" remark, made quantitative (labelled inference) -------
V36, V41 = pcv["faulty_measured_pct"], pcv["faulty_predicted_pct"]
CF_IMPLIED = float(pcv["eq3_coeff"] * V41 ** pcv["eq3_exp"])
CF_ON_LINE = float(pcv["eq3_coeff"] * V36 ** pcv["eq3_exp"])
EXCESS = CF_IMPLIED / CF_ON_LINE
print(f"p. 321 twice flags 'the point at about {V36:g}% conversion' as a possibly")
print(f"faulty determination, 'in which case the predicted value would be {V41:g}%'.")
print(f"The point itself lives only in Fig. 3 (not digitised), but eq. (3) makes the")
print(f"remark quantitative without it: a predicted {V41:g} % implies the run's measured")
print(f"carbon-on-feed was C_f = 3.55e-5 * {V41:g}^2.93 = {CF_IMPLIED:.3f} wt%, where the")
print(f"correlation at the measured {V36:g} % gives {CF_ON_LINE:.3f} wt% - i.e. the flagged")
print(f"point sits a factor {EXCESS:.2f} above the line in carbon yield. INFERENCE:")
print(f"exact given eq. (3) and the two printed percentages; not checkable further")
print(f"without the figure, and not checked further.")'''))

cells.append(md(r"""### 4. "Doubles for approximately 190–200 °F" — the paper's fit, measured per catalyst

Fig. 2 is drawn *from Table I* (p. 320: "The data in Table I permit a study of
the effect of temperature"), so testing the printed doubling interval against
Table I is a **goodness-of-fit check of the author's own correlation on the
author's own data** — labelled as such here, in `meta.yaml`, and in the README.
The doubling interval is $\Delta T_2 = \ln 2 / (d\ln C_c/dT)$, fitted per
catalyst over all 16 rows (feed-rate scatter enters as noise; that is what the
claim itself does by drawing one line per catalyst through all feed rates), and
re-estimated by a regression-free two-point route on the extreme temperature
blocks."""))

cells.append(code('''DT2_NAT, DT2_SYN = DT2["natural"], DT2["synthetic"]

def doubling_twopoint(g):
    b1 = g[g.block == 1]; b4 = g[g.block == 4]
    dT = b4.temp_F.mean() - b1.temp_F.mean()
    return float(dT * np.log(2.0) /
                 np.log(np.exp(np.log(b4.C_cat_wtpct).mean() - np.log(b1.C_cat_wtpct).mean())))

DT2_NAT_2P = doubling_twopoint(t1[t1.catalyst == "natural"])
DT2_SYN_2P = doubling_twopoint(t1[t1.catalyst == "synthetic"])

print("doubling interval of C_cat, Table I (theta = 120 min):")
print(f"  natural:   OLS {DT2_NAT:6.1f} F   two-point (block means) {DT2_NAT_2P:6.1f} F")
print(f"  synthetic: OLS {DT2_SYN:6.1f} F   two-point (block means) {DT2_SYN_2P:6.1f} F")
print(f"  printed claim: 'approximately {pcv['doubling_low_F']:g}-{pcv['doubling_high_F']:g} F.'")
print(f"\\nReading: the printed band sits BETWEEN the two catalysts - the natural")
print(f"catalyst doubles every ~{DT2_NAT:.0f} F (slower than the band), the synthetic every")
print(f"~{DT2_SYN:.0f} F (faster). As a pooled one-number summary, 190-200 F is a fair")
print(f"reading of Fig. 2; as a statement about EITHER catalyst it is outside what its")
print(f"own table supports by {100*(DT2_NAT/pcv['doubling_high_F']-1):+.0f} % (natural, vs the 200 end) and "
      f"{100*(DT2_SYN/pcv['doubling_low_F']-1):+.0f} % (synthetic, vs 190).")
print("Goodness of fit, not validation: claim and data are the same paper's.")

fig, ax = plt.subplots(figsize=(6.4, 4.0))
for cat, c, mfc in [("natural", "tab:blue", None), ("synthetic", "tab:orange", "white")]:
    g = t1[t1.catalyst == cat]
    ax.semilogy(g.temp_F, g.C_cat_wtpct, "o", color=c, mfc=mfc, label=f"{cat} (Table I)")
    tt = np.linspace(790, 970, 50)
    lnC = np.polyval(np.polyfit(g.temp_F, np.log(g.C_cat_wtpct), 1), tt)
    ax.semilogy(tt, np.exp(lnC), "-", color=c, lw=1.2,
                label=f"{cat}: doubles / {DT2[cat]:.0f} °F")
ax.set_xlabel("cracking temperature, °F"); ax.set_ylabel("C on catalyst, wt %")
ax.set_title("Table I vs the printed 'doubles per 190–200 °F'")
ax.legend(fontsize=8)
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""### 5. The diffusion hypothesis: the printed exponents against $n = 0.5$, and the $K = A^2$ defect

**The factor-2 defect, from the paper's own two lines.** Eq. (8) prints
$dC_c/d\theta = K/C_c$ "or" $C_c = A\theta^{0.5}$, and the Nomenclature prints
verbatim: **"$K$ = constant $(= A^2)$"**. The two printed lines of eq. (8) fix
the relation between $K$ and $A$ with no freedom: separate and integrate from
a clean catalyst ($C_c(0)=0$) and $C_c = \sqrt{2K\theta}$, so $K = A^2/2$ —
or differentiate the right-hand line: $dC_c/d\theta = 0.5\,A^2/C_c$, same
answer. The printed parenthesis is off by exactly 2 (equivalently $\sqrt2$ in
$A$). $K$ appears nowhere else in the paper, so nothing downstream inherits
the defect. Reported, not repaired; both derivations below, the second by an
ODE integrator that shares nothing with the symbolic route.

**How close is "close enough to 0.5"?** The paper: the four exponents "are all
close enough to 0.5 to suggest a good concordance with the theory advanced
above." The hypothesis is exercised on every printed $(A, n)$ pair by asking
what each correlation says about the quantity the hypothesis holds constant —
$C_c\,dC_c/d\theta = n A^2 \theta^{2n-1}$, constant iff $n = 0.5$ — across a
1.5-decade time window ($\theta$ = 10–300 min, the span of the fixed-bed
discussion)."""))

cells.append(code('''# --- the K = A^2 defect, two routes -----------------------------------------
th = sp.symbols("theta", positive=True)
Ksym, Asym = sp.symbols("K A_", positive=True)
Cfun = sp.Function("C")
Csol = sp.dsolve(sp.Eq(Cfun(th).diff(th), Ksym / Cfun(th)),
                 Cfun(th), ics={Cfun(sp.S(0)): 0})
Cexpr = [s.rhs for s in (Csol if isinstance(Csol, list) else [Csol])
         if (s.rhs.subs({Ksym: 1, th: 1}) > 0)][0]
K_OVER_A2 = float(sp.solve(sp.Eq(Cexpr, Asym * sp.sqrt(th)), Ksym)[0] / Asym ** 2)
print(f"route A (symbolic): dC/dtheta = K/C, C(0)=0  ->  C = {Cexpr}")
print(f"  matching C = A theta^0.5 requires K/A^2 = {K_OVER_A2:g}   (printed: 1)")

A_ex = pcv["eq2_A"]                      # any printed A works; eq. (2)'s is used

def ode_ratio(K, rtol):
    s = solve_ivp(lambda t, c: K / c[0], [1e-12, THETA_TAB],
                  [np.sqrt(2 * K * 1e-12)], rtol=rtol, atol=1e-14)
    return float(s.y[0, -1] / (A_ex * np.sqrt(THETA_TAB)))

K_NUM_RATIO = ode_ratio(A_ex ** 2, 1e-10)
K_NUM_RATIO_HALF = ode_ratio(A_ex ** 2 / 2, 1e-10)
TOL_SHIFT = max(abs(ode_ratio(A_ex ** 2, 1e-12) - K_NUM_RATIO),
                abs(ode_ratio(A_ex ** 2 / 2, 1e-12) - K_NUM_RATIO_HALF))
print(f"route B (solve_ivp, A = {A_ex:g}): integrate with the PRINTED K = A^2 and the")
print(f"  deposit overshoots A theta^0.5 by {K_NUM_RATIO:.6f} (= sqrt(2) = {np.sqrt(2):.6f});")
print(f"  with K = A^2/2 the ratio is {K_NUM_RATIO_HALF:.6f}. Tolerance-refined: tightening")
print(f"  rtol from 1e-10 to 1e-12 moves the ratios by at most {TOL_SHIFT:.1e}.")

# --- every printed (A, n) pair against the n = 0.5 invariant ----------------
TH_LO, TH_HI = 10.0, 300.0
print(f"\\nC dC/dtheta = n A^2 theta^(2n-1): variation over theta = {TH_LO:g}-{TH_HI:g} min")
print(f"{'eq.':>4s} {'A':>6s} {'n':>6s}   {'(theta_hi/lo)^(2n-1)':>22s}   {'drift from constant':>20s}")
drifts = {}
for eq, (Ak, nk) in {"(1)": ("eq1_A", "eq1_n"), "(2)": ("eq2_A", "eq2_n"),
                     "(6)": ("eq6_A", "eq6_n"), "(7)": ("eq7_A", "eq7_n")}.items():
    Av, nv = pcv[Ak], pcv[nk]
    ratio = (TH_HI / TH_LO) ** (2 * nv - 1)
    drifts[eq] = 100 * abs(ratio - 1)
    print(f"{eq:>4s} {Av:6.2f} {nv:6.2f}   {ratio:22.3f}   {drifts[eq]:19.1f} %")
DRIFT_WORST = max(drifts.values())
print(f"\\nReading: under the hypothesis's own invariant, the printed exponents are not")
print(f"'close to 0.5' uniformly - eq. (6)'s n = 0.38 drifts {drifts['(6)']:.0f} % and eq. (7)'s")
print(f"n = 0.53 only {drifts['(7)']:.0f} % across the same window. The paper's qualitative")
print(f"suggestion survives; a quantitative n = 0.5 does not, and the paper itself")
print(f"never claims it does ('n = constant, depending only slightly on...').")

# --- A is local, n is portable: the paper's own two studies -----------------
EQ2_AT_120 = float(voorhies(THETA_TAB, pcv["eq2_A"], pcv["eq2_n"]))
blk = t1[(t1.catalyst == "synthetic") & (t1.block == 2)]
GAP_A = float(blk.C_cat_wtpct.mean() / EQ2_AT_120)
print(f"\\nThe same nominal system in the paper's two campaigns: eq. (2) (Fig.-1 study,")
print(f"one repeatedly regenerated batch) predicts C_c({THETA_TAB:.0f} min) = {EQ2_AT_120:.2f} wt%;")
print(f"Table I's synthetic ~{blk.temp_F.mean():.0f} F block (fresh batch each run) prints "
      f"{blk.C_cat_wtpct.min():.1f}-{blk.C_cat_wtpct.max():.1f} wt%")
print(f"- a factor {GAP_A:.2f}. East Texas gas oil + synthetic catalyst + 850 F + 120 min in")
print(f"both cases; only the unit and the catalyst batch differ. That is the")
print(f"Nomenclature's asymmetry made concrete: A moved ~{100*(GAP_A-1):.0f} % between campaigns")
print(f"of the same laboratory, n stays in a 0.38-0.53 band across two modes, three")
print(f"oils and two catalyst types. A is local; n is the portable part of the law.")
print(f"(The four A values span {min(pcv[k] for k in ('eq1_A','eq2_A','eq6_A','eq7_A')):g}-"
      f"{max(pcv[k] for k in ('eq1_A','eq2_A','eq6_A','eq7_A')):g} - a factor "
      f"{max(pcv[k] for k in ('eq1_A','eq2_A','eq6_A','eq7_A'))/min(pcv[k] for k in ('eq1_A','eq2_A','eq6_A','eq7_A')):.1f} -")
print(f" while the exponents span {min(n_printed.values()):g}-{max(n_printed.values()):g}, "
      f"±{100*(max(n_printed.values())-min(n_printed.values()))/(max(n_printed.values())+min(n_printed.values())):.0f} % about their midpoint.)")'''))

# --------------------------------------------------------------- validation
cells.append(md(r"""## Validation

No measurement outside the paper exists on this page, so validation means:
(a) internal identities the printed tables must satisfy, (b) exact
reproduction of printed algebra by independent routes, and (c) a break table
for **every** reported metric — each row injects a defect a reader of a 1945
bilevel scan could actually commit (a 3↔8 or 5↔6 glyph swap, a dropped
superscript, a swapped column, a unit misreading) and shows which numbers
move. The coverage map at the end is asserted key-for-key against
`agreement.json`.

### Check 1 — the mass-balance identity, and what breaks it

The implied-$W/U$ constancy of section 2 is the one identity the tables
impose on themselves. It is **near-structural for the slope comparison** (if
the columns are linked, $m_{feed}+1 \approx m_{cat}$ identically — that
consequence is labelled, not sold as evidence) but it has real teeth as a
transcription check: a single mis-read carbon digit shifts its own row's ratio
well clear of its block's scatter — the table below quantifies by how much
(the as-transcribed baseline is not at zero either: its worst row is the
reported 1.1-glyph tension row)."""))

cells.append(code('''base = t1.copy()
def wu_flag(df, cat, blk):
    g = df[(df.catalyst == cat) & (df.block == blk)]
    r = implied_WU(g)
    return float(100 * r.std(ddof=0) / r.mean()), float(100 * np.max(np.abs(r / r.mean() - 1)))

b3 = (t1.catalyst == "synthetic") & (t1.block == 3)
print(f"break table 1 - the identity as transcription guard "
      f"(synthetic {t1[b3].temp_F.mean():.0f} F block):")
print(f"{'injected defect':52s} {'CV %':>6s} {'worst row dev %':>15s}")
cases = [("none - as transcribed", None)]
pert = t1.copy(); pert.loc[b3 & (t1.temp_F == 906), "C_feed_wtpct"] = 6.3   # 8 -> 3
cases.append(("C_feed 6.8 mis-read as 6.3 (8->3 glyph)", pert))
pert = t1.copy(); pert.loc[b3 & (t1.temp_F == 911), "C_cat_wtpct"] = 11.3   # 8 -> 3
cases.append(("C_cat 11.8 mis-read as 11.3 (8->3 glyph)", pert))
pert = t1.copy(); pert.loc[b3 & (t1.temp_F == 911), "feed_rate_vvhr"] = 1.0
cases.append(("feed rate 1.1 read as 1.0 (the settled glyph)", pert))
pert = t1.copy()
i = pert.index[b3]
pert.loc[i, ["C_feed_wtpct", "C_cat_wtpct"]] = pert.loc[i, ["C_cat_wtpct", "C_feed_wtpct"]].values
cases.append(("the two carbon columns swapped in the block", pert))
CV0, W0 = wu_flag(t1, "synthetic", 3)
worsts = {}
for name, df in cases:
    cv, worst = wu_flag(df if df is not None else t1, "synthetic", 3)
    worsts[name] = worst
    print(f"{name:52s} {cv:6.1f} {worst:15.1f}")
GLYPH_RAISE = min(worsts[k] / W0 for k in list(worsts)[1:3])
print(f"\\nEach carbon-digit defect raises the block's worst-row deviation by a factor")
print(f">= {GLYPH_RAISE:.1f} over the baseline {W0:.1f} %; the column swap is catastrophic "
      f"(ratio ~ 1/(U theta/60)^2 off).")
print("Note the 1.1->1.0 row: it IMPROVES the identity (see the data section) - which")
print("is exactly why the glyph was re-read at 4x and the printed 1.1 kept, with the")
print("tension reported. The identity cannot adjudicate a glyph this small; the")
print("pixels can, and did.")
print("\\nLABELLED STRUCTURAL for the slope comparison: given the linkage,")
mm = [(slope_ols(g.feed_rate_vvhr, g.C_feed_wtpct) + 1 - slope_ols(g.feed_rate_vvhr, g.C_cat_wtpct))
      for _, g in t1.groupby(["catalyst", "block"])]
print(f"m_feed + 1 - m_cat should be ~0 by construction: worst block "
      f"{np.abs(mm).max():.3f} (it is not independent evidence for the claim).")'''))

cells.append(md(r"""### Check 2 — the eq.-(5) chain: two solvers, and which misreadings it can exclude

The symbolic and numeric routes above agree to 5 decimals on all three
constants — a check on the two solvers; the isolation step they share is
checked against the printed triple itself. Break table 2a injects the large
misreadings a bilevel 1945 scan invites; table 2b sweeps **every** last-digit
±1 alternative of the five inputs and shows which of them the chain does and
does not pin — each row against the *printed* (96, 0.34, 0.19):"""))

cells.append(code('''print(f"routes A and B: coeff {COEFF_SYM:.5f} vs {COEFF_NUM:.5f}, "
      f"expU {EXP_U_SYM:.5f} vs {EXP_U_NUM:.5f}, expT {EXP_T_SYM:.5f} vs {EXP_T_NUM:.5f}")
print(f"largest route disagreement: "
      f"{max(abs(COEFF_SYM-COEFF_NUM)/COEFF_SYM, abs(EXP_U_SYM-EXP_U_NUM), abs(EXP_T_SYM-EXP_T_NUM)):.2e}")

def chain(A2v=0.65, n2v=0.44, c3v=3.55e-5, e3v=2.93, drv=0.58):
    coeff = (A2v * 60 * drv / c3v) ** (1 / e3v)
    return coeff, 1 / e3v, (1 - n2v) / e3v

def hits_printed(c, eu, et):
    return (round(c) == 96) and (round(eu, 2) == 0.34) and (round(et, 2) == 0.19)

print(f"\\nbreak table 2a - leading-digit, transposed-glyph and scale misreadings:")
print(f"{'injected defect':46s} {'coeff':>9s} {'expU':>7s} {'expT':>7s} {'match?':>7s}")
variants = [
    ("none - constants as read at digit scale", {}),
    ("2.93 mis-read as 2.98 (3<->8 glyph)", dict(e3v=2.98)),
    ("2.93 mis-read as 2.63 (9<->6 glyph)", dict(e3v=2.63)),
    ("0.58 mis-read as 0.53 (8<->3 glyph)", dict(drv=0.53)),
    ("10^-5 mis-read as 10^-6 (superscript dot)", dict(c3v=3.55e-6)),
    ("0.44 mis-read as 0.41 (eq. 1's exponent)", dict(n2v=0.41)),
    ("0.65 mis-read as 0.85", dict(A2v=0.85)),
    ("(not a defect) n = 0.5, the hypothesis value", dict(n2v=0.5)),
]
coef_devs = []
for name, kw in variants:
    c, eu, et = chain(**kw)
    ok = hits_printed(c, eu, et)
    if kw and "n2v" not in kw:
        coef_devs.append(abs(c / 96 - 1))
    print(f"{name:46s} {c:9.2f} {eu:7.4f} {et:7.4f} {'YES' if ok else 'no':>7s}")

print(f"\\nbreak table 2b - EVERY last-digit +-1 alternative of the five inputs:")
print(f"{'alternative reading':46s} {'coeff':>9s} {'expU':>7s} {'expT':>7s} {'match?':>7s}")
survivors, n_pass, n_tested = [], 0, 0
for label, key, base, step in [("0.65", "A2v", 0.65, 0.01), ("0.44", "n2v", 0.44, 0.01),
                               ("3.55e-5", "c3v", 3.55e-5, 1e-7),
                               ("2.93", "e3v", 2.93, 0.01), ("0.58", "drv", 0.58, 0.01)]:
    for s in (+step, -step):
        c, eu, et = chain(**{key: base + s})
        ok = hits_printed(c, eu, et)
        n_tested += 1
        if ok:
            n_pass += 1
            if label not in survivors:
                survivors.append(label)
        print(f"{label + ' read as ' + format(base + s, 'g'):46s} "
              f"{c:9.2f} {eu:7.4f} {et:7.4f} {'YES' if ok else 'no':>7s}")

print(f"\\nWhat the chain certifies - and what it cannot. Every table-2a misreading")
print(f"misses the printed triple: the perturbed coefficients miss 96 by "
      f"{100*min(coef_devs):.0f} % to a")
print(f"factor {1+max(coef_devs):.1f}, the perturbed exponents miss a printed second decimal, and")
print(f"2.93's last digit is pinned as well (both its +-0.01 rows in 2b miss 96). So")
print(f"the re-derivation excludes every leading-digit, transposed-glyph and")
print(f"exponent-scale misreading of the five inputs. It does NOT pin the last digit")
print(f"of the other four: {n_pass} of the {n_tested} last-digit alternatives in 2b "
      f"({', '.join(survivors)}")
print(f"each have at least one) still round to the printed (96, 0.34, 0.19). Those")
print(f"final digits rest on the digit-scale crops alone - read independently twice,")
print(f"at transcription and again at verification, in agreement. (n = 0.5 appears in")
print(f"2a to show eq. (5) is sensitive to the 0.44-vs-0.5 distinction the diffusion")
print(f"section blurs.)")'''))

cells.append(md(r"""### Check 3 — can the slope test actually fail?

Three things could fake the feed-rate verdict: shared code between the
compared numbers (excluded — the two columns get the same estimator, but the
*hypotheses* differ by a full unit of exponent), a confound (temperature,
bounded in section 1 with the paper's own rule), and resolving power (the
rounding floor). The break rows:"""))

cells.append(code('''cet = t2[t2.feed == "cetane"]
print(f"{'injected defect':52s} {'m_cat':>8s} {'m_feed':>8s}")
print(f"{'none - as transcribed (cetane, 4x span)':52s} {M_CAT_CET:+8.3f} {M_FEED_CET:+8.3f}")
sw = slope_ols(cet.feed_rate_vvhr, cet.C_feed_wtpct), slope_ols(cet.feed_rate_vvhr, cet.C_cat_wtpct)
print(f"{'carbon columns swapped':52s} {sw[0]:+8.3f} {sw[1]:+8.3f}")
p = cet.copy(); p.loc[p.feed_rate_vvhr == 2.4, "C_cat_wtpct"] = 4.2   # 7 -> 2 glyph
print(f"{'C_cat 4.7 mis-read as 4.2 at U = 2.4':52s} "
      f"{slope_ols(p.feed_rate_vvhr, p.C_cat_wtpct):+8.3f} {M_FEED_CET:+8.3f}")
p = cet.copy(); p.loc[p.feed_rate_vvhr == 2.4, "feed_rate_vvhr"] = 1.4  # 2 -> 1 glyph
print(f"{'feed rate 2.4 mis-read as 1.4 (2->1 glyph)':52s} "
      f"{slope_ols(p.feed_rate_vvhr, p.C_cat_wtpct):+8.3f} "
      f"{slope_ols(p.feed_rate_vvhr, p.C_feed_wtpct):+8.3f}")
print(f"{'throughput hypothesis injected (C_cat ~ U)':52s} "
      f"{slope_ols(cet.feed_rate_vvhr, cet.C_cat_wtpct * cet.feed_rate_vvhr / 0.6):+8.3f}"
      f" {'':>8s}")
m2 = np.array([slope_ols(g.feed_rate_vvhr, g.C_cat_wtpct)
               for _, g in t1.groupby(["catalyst", "block"])])
print(f"{'T-correction sign FLIPPED (Table I mean m_cat)':52s} "
      f"{tcorr_slopes(-1.0).mean():+8.3f}  (raw {m2.mean():+.3f}, corrected {M_CAT_T1_TCORR:+.3f})")

floors_blk = np.array([slope_rounding_floor(g.feed_rate_vvhr, g.C_cat_wtpct)
                       for _, g in t1.groupby(["catalyst", "block"])])
RAW_FLOORS_WORST = float(np.max(np.abs(m_raw) / floors_blk))
CORR_FLOORS_WORST = float(np.max(np.abs(tcorr_slopes(+1.0)) / floors_blk))
print(f"\\nresolving power: the floors ARE the verdict's error bar. Table I's 2x span")
print(f"resolves m to ±{FLOOR_T1:.3f} at worst. The raw block slopes reach "
      f"{RAW_FLOORS_WORST:.1f} floors from")
print(f"zero - too large for rounding alone - but section 1 showed most of that is the")
print(f"within-block temperature confound: corrected with the paper's own doubling")
print(f"rule, every block sits within {CORR_FLOORS_WORST:.1f} floors of zero, with both signs "
      f"occurring.")
FLOOR_CET = slope_rounding_floor(cet.feed_rate_vvhr, cet.C_cat_wtpct)
print(f"Table II's 4x span tightens the floor to ±{FLOOR_CET:.3f} (cetane; decalin "
      f"±{FLOOR_T2:.3f});")
print(f"cetane's m_cat = {M_CAT_CET:+.3f} and the endpoint estimator's {M_CAT_CET_EP:+.3f} "
      f"sit within")
print(f"{max(abs(M_CAT_CET), abs(M_CAT_CET_EP)) / FLOOR_CET:.1f} floors of zero on "
      f"cetane's own floor.")
print(f"So the page reports 'consistent with 0, incompatible with +1' - it does NOT")
print(f"report 'proved exactly zero': a true residual |m_cat| of a few hundredths")
print(f"cannot be excluded at this printing precision, which is presumably what the")
print(f"paper's own 'within limits' hedge is for.")
print(f"The claim's own falsifier (m_cat -> +1 if carbon tracked throughput)")
print(f"moves the metric by ~{abs(slope_ols(cet.feed_rate_vvhr, cet.C_cat_wtpct*cet.feed_rate_vvhr/0.6) - M_CAT_CET):.2f} - "
      f"~{abs(slope_ols(cet.feed_rate_vvhr, cet.C_cat_wtpct*cet.feed_rate_vvhr/0.6) - M_CAT_CET)/FLOOR_CET:.0f} floors of cetane's own. The test can fail; it does not.")'''))

cells.append(md(r"""### Check 4 — the doubling claim's breaks

The claim is a fit to Table I; the break rows establish that the measured
doubling interval is sensitive to exactly the misreadings that would corrupt
it, and that the two estimator routes agree."""))

cells.append(code('''nat = t1[t1.catalyst == "natural"]
print(f"{'injected defect':52s} {'DT2 natural, F':>14s}")
print(f"{'none - as transcribed':52s} {DT2_NAT:14.1f}")
p = nat.copy(); p.loc[p.temp_F == 962, "temp_F"] = 902.0        # 6 -> 0 glyph
print(f"{'962 F mis-read as 902 F (6->0 glyph)':52s} {doubling_dT(p.temp_F, p.C_cat_wtpct):14.1f}")
p = nat.copy()
print(f"{'C_feed column used instead of C_cat':52s} {doubling_dT(p.temp_F, p.C_feed_wtpct):14.1f}")
DT2_C = doubling_dT(nat.temp_F * 5 / 9, nat.C_cat_wtpct)  # column taken as C then fit in C
print(f"{'temperatures taken as Celsius (fit in C units)':52s} {DT2_C:14.1f}")
print(f"{'two-point route (not a defect)':52s} {DT2_NAT_2P:14.1f}")
DT2_NAT_FEED = doubling_dT(nat.temp_F, nat.C_feed_wtpct)
DT2_NAT_LINKED = doubling_dT(nat.temp_F, nat.C_cat_wtpct / nat.feed_rate_vvhr)
print(f"\\nEvery row moves the metric, including the wrong-column one: C_feed doubles")
print(f"every {DT2_NAT_FEED:.0f} F, not {DT2_NAT:.0f} F. Even an EXACTLY constant Check-1 ratio would")
print(f"not make the two answers coincide: it forces C_feed = C_cat/(k U theta/60),")
print(f"and because feed rate happens to correlate with temperature across the pooled")
print(f"rows, the -ln U term alone moves the fit from {DT2_NAT:.0f} to {DT2_NAT_LINKED:.0f} F "
      f"(OLS of")
print(f"ln(C_cat/U) on T). The remaining {DT2_NAT_LINKED - DT2_NAT_FEED:.0f} F down to the "
      f"measured {DT2_NAT_FEED:.0f} F is the")
print(f"cross-block drift of Check 1's implied ratio. So a column mix-up would not go")
print(f"unnoticed - and the printed doubling claim is specifically a claim about")
print(f"carbon ON CATALYST.")'''))

cells.append(md(r"""### Blind spots — what this page does *not* establish

1. **The time exponent $n$ is never tested against data here.** Tables I and
   II sit at a single $\theta = 120$ min; the carbon-vs-time evidence is
   figure-only (Figs. 1 and 5) and was deliberately not digitised. Every $n$
   on this page is a transcription of the author's fit. A page that *tests*
   $n$ needs those figures through the review gate, or another paper's
   tables.
2. **Feed-rate independence is shown inside 0.6–2.4 v/v/hr at 120 min, on one
   gas oil and two pure compounds, at the paper's own 2-figure printing.**
   The paper's "within limits" hedge is honest and this page inherits it: the
   raw Table-I block slopes reach −0.12, most of that is the within-block
   temperature confound (bounded with the paper's own doubling rule), and
   what remains after correction — a few hundredths, of both signs — is
   within ~3 rounding floors and cannot be adjudicated at this precision.
3. **The claim's evidence is one column, not two.** Since the two carbon
   columns are linked by the run's mass ratio, $m_{feed}\approx-1$ is not
   corroboration of $m_{cat}\approx0$ — it is the same fact. This page's
   contribution is establishing the linkage, not doubling the evidence.
4. **The doubling test is goodness of fit** (claim and data share a source),
   and its per-catalyst split (206 vs 176 °F) is itself at the edge of what
   16 rows of 2-figure data determine — the two-point routes land at 208/174,
   so the split is robust in direction, approximate in size.
5. **Nothing validates eq. (3) or eq. (5) against measurement** — the runs
   behind them are figure-only. What is established is that eq. (5) is
   *exactly* what eqs. (2)+(3)+(4) imply, i.e. the paper's algebra is clean,
   and that the printed triple is its correct rounding. The re-derivation
   excludes every leading-digit, transposed-glyph and exponent-scale
   misreading of its five inputs (and pins 2.93's last digit) — but *not* the
   remaining last digits, most of whose ±1 alternatives round to the same
   printed triple (break table 2b); those digits rest on the digit-scale
   crops, read independently twice.
6. **The $K = A^2$ defect is a nomenclature error with no consequences in the
   paper** — $K$ occurs nowhere else. It is reported because the case's whole
   content is the printed constants, and this is the one printed constant
   relation that is wrong.
7. **The implied-$W/U$ temperature drift (~5–8 % down across blocks) is
   unexplained.** It is larger than the within-block scatter, it is not
   feed-rate structure, and no reading offered here accounts for it; it is
   recorded, not resolved."""))

# ------------------------------------------------------------------ metrics
cells.append(code('''metrics = dict(
    # section 1 - feed-rate exponents
    m_feed_tableI_mean=M_FEED_T1, m_cat_tableI_mean=M_CAT_T1,
    m_cat_tableI_worst_block=M_CAT_T1_WORST,
    m_cat_tableI_tempcorrected_mean=M_CAT_T1_TCORR,
    m_feed_tableII_cetane=M_FEED_CET, m_cat_tableII_cetane=M_CAT_CET,
    m_feed_tableII_decalin=M_FEED_DEC, m_cat_tableII_decalin=M_CAT_DEC,
    m_cat_tableII_cetane_endpoint=M_CAT_CET_EP,
    slope_rounding_floor_tableI=FLOOR_T1, slope_rounding_floor_tableII=FLOOR_T2,
    # section 2 - the identity
    WU_implied_natural=WU_NAT, WU_implied_synthetic=WU_SYN,
    WU_implied_cetane=WU_CET, WU_implied_decalin=WU_DEC,
    WU_within_group_cv_max_pct=CV_MAX,
    WU_decalin_over_cetane=WU_RATIO,
    WU_block_drift_natural_pct=DRIFT_NAT, WU_block_drift_synthetic_pct=DRIFT_SYN,
    # section 3 - the chain
    eq5_coeff_symbolic=COEFF_SYM, eq5_expU_symbolic=EXP_U_SYM, eq5_expTheta_symbolic=EXP_T_SYM,
    eq5_coeff_numeric=COEFF_NUM, eq5_expU_numeric=EXP_U_NUM, eq5_expTheta_numeric=EXP_T_NUM,
    eq5_printed_rounding_worst_pct=ROUND_COST,
    inst_over_avg_conversion=INST_OVER_AVG,
    faulty_point_implied_Cf_wtpct=CF_IMPLIED,
    faulty_point_online_Cf_wtpct=CF_ON_LINE,
    faulty_point_excess_ratio=EXCESS,
    # section 4 - doubling
    doubling_dT_natural_F=DT2_NAT, doubling_dT_synthetic_F=DT2_SYN,
    doubling_dT_natural_twopoint_F=DT2_NAT_2P, doubling_dT_synthetic_twopoint_F=DT2_SYN_2P,
    # section 5 - hypothesis and portability
    K_over_A2_required=K_OVER_A2,
    K_printed_numeric_overshoot=K_NUM_RATIO,
    n_printed_min=min(n_printed.values()), n_printed_max=max(n_printed.values()),
    hypothesis_drift_worst_pct=DRIFT_WORST,
    eq2_prediction_at_120min_wtpct=EQ2_AT_120,
    A_campaign_gap_ratio=GAP_A,
)
_ = report_agreement("B2.1", metrics)'''))

cells.append(code('''# ---- break-row coverage: every metric names the row that moves it -----------
# rows: (a) column swap [Chk 1+3]  (b) single-digit glyph perturbations [Chk 1/3/4]
# (c) 1.1->1.0 feed-rate glyph [Chk 1]  (d) T-correction sign flip [Chk 3]
# (e) chain constant misreadings [Chk 2]  (f) K = A^2 vs A^2/2 in the ODE [sec 5]
# (g) 962->902 / Celsius / wrong-column rows [Chk 4]
# (h) throughput-hypothesis injection [Chk 3]  STRUCT = labelled structural.
CHK1, CHK2, CHK3, CHK4 = "Chk 1 rows", "Chk 2 rows", "Chk 3 rows", "Chk 4 rows"
COVERAGE = {
    "m_feed_tableI_mean": f"{CHK3} class: column swap / glyph rows move the sibling "
                          "cetane metrics; same estimator, same columns",
    "m_cat_tableI_mean": f"{CHK3} class (glyph/swap rows); the correction rows bracket "
                         f"it ({tcorr_slopes(-1.0).mean():+.3f} .. {M_CAT_T1_TCORR:+.3f})",
    "m_cat_tableI_worst_block": f"{CHK3} class (glyph row on any block moves it)",
    "m_cat_tableI_tempcorrected_mean": f"(d) directly (sign flip: {M_CAT_T1_TCORR:+.3f} -> "
                                       f"{tcorr_slopes(-1.0).mean():+.3f})",
    "m_feed_tableII_cetane": f"{CHK3} (a) directly (swap: {M_FEED_CET:+.3f} <-> {M_CAT_CET:+.3f})",
    "m_cat_tableII_cetane": f"{CHK3} (a),(b),(h) directly (throughput injection moves it "
                            "to ~+0.95)",
    "m_feed_tableII_decalin": f"{CHK3} class (same estimator/columns as cetane rows)",
    "m_cat_tableII_decalin": f"{CHK3} class",
    "m_cat_tableII_cetane_endpoint": f"{CHK3} (b) directly (the U=2.4 glyph rows move "
                                     "exactly the endpoints it uses)",
    "slope_rounding_floor_tableI": "STRUCT: the resolving-power bound itself; moves only "
                                   "with the printed precision (0.05) and the U design",
    "slope_rounding_floor_tableII": "STRUCT: as above",
    "WU_implied_natural": f"{CHK1} (b),(c) class (any digit in any of 16 rows moves it)",
    "WU_implied_synthetic": f"{CHK1} (b),(c) directly (rows shown)",
    "WU_implied_cetane": f"{CHK1} class via Table II columns; Chk 3 (b) rows share cells",
    "WU_implied_decalin": f"{CHK1} class",
    "WU_within_group_cv_max_pct": f"{CHK1} (b) directly ({CV0:.1f} -> at least doubled "
                                  "by every glyph row)",
    "WU_decalin_over_cetane": f"{CHK1} class (any Table II carbon digit moves it)",
    "WU_block_drift_natural_pct": f"{CHK1} class; reported-unexplained (blind spot 7)",
    "WU_block_drift_synthetic_pct": f"{CHK1} class; reported-unexplained (blind spot 7)",
    "eq5_coeff_symbolic": f"{CHK2} (2a) directly (every leading-digit/scale row misses "
                          "the printed 96; 2b names the last digits it cannot pin)",
    "eq5_expU_symbolic": f"{CHK2} directly (2.93 rows move it in the 2nd decimal)",
    "eq5_expTheta_symbolic": f"{CHK2} directly (0.44 and n=0.5 rows move it)",
    "eq5_coeff_numeric": f"{CHK2} class (same inputs; solver-independent route - its "
                         "agreement with the symbolic value to 5 decimals checks the "
                         "solvers; the printed 96 checks the shared isolation)",
    "eq5_expU_numeric": f"{CHK2} class (as above)",
    "eq5_expTheta_numeric": f"{CHK2} class (as above)",
    "eq5_printed_rounding_worst_pct": f"{CHK2} class (any row that moves the exact chain "
                                      "moves the printed-vs-exact gap it measures)",
    "inst_over_avg_conversion": f"{CHK2} (0.44/n=0.5 rows move expTheta, hence 1-expTheta)",
    "faulty_point_implied_Cf_wtpct": f"{CHK2} class (3.55e-5 and 2.93 rows move it directly)",
    "faulty_point_online_Cf_wtpct": f"{CHK2} class (as above)",
    "faulty_point_excess_ratio": f"{CHK2} class (2.93 rows move it; the coefficient cancels)",
    "doubling_dT_natural_F": f"{CHK4} directly (962->902: {DT2_NAT:.0f} -> "
                             f"{doubling_dT(nat.assign(temp_F=nat.temp_F.where(nat.temp_F != 962, 902.0)).temp_F, nat.C_cat_wtpct):.0f} F; "
                             f"Celsius misreading -> {doubling_dT(nat.temp_F * 5 / 9, nat.C_cat_wtpct):.0f}; "
                             f"wrong column -> {DT2_NAT_FEED:.0f})",
    "doubling_dT_synthetic_F": f"{CHK4} class (same estimator on the sibling catalyst)",
    "doubling_dT_natural_twopoint_F": f"{CHK4} class (independent 2nd estimator; agreement "
                                      "with OLS to ~2 % is the route check)",
    "doubling_dT_synthetic_twopoint_F": f"{CHK4} class",
    "K_over_A2_required": "sec. 5 route pair: STRUCT as an exact symbolic result (0.5 "
                          "identically); its numeric twin below is the row that moves",
    "K_printed_numeric_overshoot": f"sec. 5 (f) directly (K = A^2: {K_NUM_RATIO:.4f}; "
                                   f"K = A^2/2: {K_NUM_RATIO_HALF:.4f})",
    "n_printed_min": "STRUCT: transcription constant; pinned by the digit-scale crops AND "
                     "reconciled against B2.2's independent transcription of F&B's quote",
    "n_printed_max": "STRUCT: as above",
    "hypothesis_drift_worst_pct": f"{CHK2} class (any n misreading moves 2n-1); n=0.5 row "
                                  "sends it to 0 by construction",
    "eq2_prediction_at_120min_wtpct": f"{CHK2} class (0.65/0.44 glyph rows move it)",
    "A_campaign_gap_ratio": f"{CHK2} class + {CHK1} (b) class (either side's digits move it)",
}
missing, extra = set(metrics) - set(COVERAGE), set(COVERAGE) - set(metrics)
assert not missing and not extra, (missing, extra)
print("break-row coverage (which injected defect moves each reported metric):")
for k in metrics:
    print(f"  {k:<34} {COVERAGE[k]}")
ABS_FLOOR = 1e-12
below = sorted(k for k, v in metrics.items() if abs(v) < ABS_FLOOR)
print(f"\\nmetrics below check_agreement's ABS_FLOOR = 1e-12 (outside the CI suite): "
      f"{below if below else 'none'}")'''))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**To the law itself, nothing, and this page does not pretend otherwise.**
$C_c = A\theta^n$ is a two-constant power function; no operator, grid or
solver appears above, and every section would run with pymrm uninstalled —
exactly as on `A1.6`, `A1.1` and `F1.4`, which set the precedent of saying so
where a reader expects a solver. What the reimplementation adds is the
quantification a 1945 paper asserts in words:

- **"Independent of feed rate" gets a number and an error bar.** The printed
  word becomes $m_{cat}$, measured per group against the two hypotheses an
  exponent apart, with the resolving power of a 2-significant-figure table
  computed rather than assumed — the verdict is "within a few rounding floors
  of 0, tens of floors from +1", which is both stronger and more honest than
  "substantially constant".
- **The tables' internal structure is made explicit.** The two carbon columns
  are linked by the run's mass ratio at the rounding floor — so the claim has
  one witness per table, the linkage constant is recovered (a number the paper
  never prints for these units), and one systematic residual (the ~5–8 %
  temperature drift of that constant) is surfaced that a reader of the printed
  table would not see.
- **The printed algebra is verified, and the printed prose made exact.** The
  eq.-(5) triple is re-derived to 5 decimals by symbolic and numeric routes
  with independent solvers — which also pins the five input constants against
  every leading-digit, transposed-glyph and exponent-scale misreading (their
  last digits are within the printed rounding and rest on the digit-scale
  crops, read independently twice — break table 2b); the "faulty point"
  remark and the doubling claim become numbers with stated provenance; and
  the derived correlation gains one small exact object the paper does not
  print (instantaneous conversion = 0.809 × the period average, everywhere).
- **One printed defect is settled from two printed lines**: the Nomenclature's
  "$K$ = constant $(= A^2)$" is inconsistent with eq. (8) by exactly a factor
  2, shown symbolically and by an independent ODE integration, and shown to be
  consequence-free within the paper.

The reactor-scale story this law seeds — coke as a *profile*, activity as a
function of local carbon — is where the operators live, and that is
[`B2.2`](../B2.2-froment-bischoff-coking/)."""))

# ------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**Use $C_c = A\theta^n$ as a first-cut coke inventory** for regenerative
cracking-type systems: carbon on catalyst as a function of catalyst residence
time alone, with feed rate entering only through the carbon-*yield*
bookkeeping ($C_f = C_c/(W\theta/60)$). Two warnings, both measured above:

- **$n$ travels; $A$ does not.** Across this paper's own two fixed-bed
  campaigns on the *same oil, catalyst type and temperature*, $A$ moved by a
  factor ~1.8; across modes and feeds the printed $A$ span a factor 3.6 while
  $n$ stays in 0.38–0.53. Take an exponent from the literature if you must;
  take $A$ only from your own unit.
- **The feed-rate independence is bounded**: demonstrated over 0.6–2.4 v/v/hr
  at $\theta = 120$ min, and the paper itself declines to claim it "for all
  feed stocks, all catalysts, or an extreme variation in feed rates".

**If your deactivation model consumes this law downstream** (activity vs
time-on-stream), read [`B2.2`](../B2.2-froment-bischoff-coking/) first: once
coke forms a profile along a bed, a time-only activity is ill-posed, and
`B2.2` maps which coking mechanisms can and cannot reproduce the exponent band
printed here. `B2.3` (Levenspiel) and `B2.4` (Beeckman–Froment) continue the
ladder; `B1.7` (Mears) is the disk-neighbour whose first page quotes this
paper — and mis-identifies itself to a careless reader.

**Transcription reuse:** the three CSVs under `data/` carry Tables I and II
and every printed constant with per-cell provenance; the mass-balance identity
of Check 1 re-verifies the transcription automatically if you edit them.

**Cite the source, not this page:** Voorhies, A. Jr., *Carbon Formation in
Catalytic Cracking*, Industrial & Engineering Chemistry **37**(4) 318–322
(1945), [doi:10.1021/ie50424a010](https://doi.org/10.1021/ie50424a010).

**Fit/test, one line:** nothing on this page is fitted by this page; every
comparison is either exact reproduction of printed algebra or a consistency
check of the paper's printed claims against the paper's own printed tables,
and is labelled so wherever it appears."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python",
                   "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py", "mimetype": "text/x-python", "name": "python",
        "nbconvert_exporter": "python", "pygments_lexer": "ipython3",
        "version": "3.13.5"},
}
out = Path(__file__).with_name("index.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
