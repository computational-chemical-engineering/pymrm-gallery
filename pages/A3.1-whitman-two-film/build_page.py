#!/usr/bin/env python3
"""Generate index.ipynb for page A3.1 (Whitman's two-film theory). Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "Whitman's two-film theory, and what his three runs can actually test"
description: "The resistance-in-series picture that every absorption calculation still uses, rebuilt as the diffusion problem it claims to be, and then measured against the three runs Whitman published with it. His equilibrium curve is reconstructed from four equilibrium pairs hidden in his own table rather than digitised from his figure; his worked example is reproduced end to end with his printed numbers held out; and the one prediction he holds out is shown to be 99.4 % gas-film controlled, so deleting the liquid film entirely moves it by 0.6 % and moves it closer to the measurement."
categories: [sec:A, struct:S3, tier:T0, data:tier3, phase:gas-liquid]
date: 2026-08-05
---

# Whitman's two-film theory, and what his three runs can actually test

**Catalog ID:** `A3.1` · **Structures:** `S3` (1D steady BVP) · **Tier:** T0

Two films in series, one equilibrium at the interface between them, and the
absorption rate is whatever both films can pass at once:

$$\frac{\mathrm{d}W}{\mathrm{d}\theta} \;=\; k_p\,(p_1-p_2)\;=\;k_c\,(c_2-c_3),
\qquad p_2 = f(c_2).$$

That is the whole of Whitman (1923). It is one of the most-cited results in
chemical engineering, and it is an *algebraic* statement, not a differential
one — which makes it dangerously easy to "validate" by rearranging it.

So this page does three separate things and keeps them apart.

1. **Reproduces Whitman's worked example** — the interface pair he had to find
   by hand on a chart, the liquid-film coefficient he derived from it, and the
   rate he predicted for a run he held out. All of it from Table 1's printed
   numbers, with **his equilibrium curve reconstructed rather than digitised**.
2. **Rebuilds the two films as an actual steady-diffusion problem in pymrm**,
   with the equilibrium jump at the interface as a coupling condition, and says
   plainly that agreement with the algebra is an identity.
3. **Measures what his data can test.** Three runs, two coefficients. The one
   held-out prediction turns out to be 99.4 % controlled by the gas film, so
   **deleting the liquid film entirely** moves it by 0.6 % — and moves it
   *closer* to the observed rate. Sherwood and Pigford wrote in 1951 that the
   theory "has never been adequately checked experimentally". This page puts a
   number on why."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

By 1923 everyone agreed the absorption rate had the form

$$\frac{\mathrm{d}W}{\mathrm{d}\theta} = \text{coefficient} \times \text{driving potential},
\tag{1}$$

and nobody agreed on the potential. One camp (Lewis 1916; Whitman & Keats 1922)
wrote it as a **pressure** difference, the other (Donnan & Masson 1920; Van
Arsdel 1920) as a **concentration** difference:

$$\text{Driving potential} = p_g-p_l = p, \tag{2}\qquad\qquad
\text{Driving potential} = c_g-c_l = c. \tag{3}$$

with, in his words, $p_l$ = "the partial pressure of solute exerted by the
liquid" and $c_g$ = "the concentration of liquid which would be in equilibrium
with the gas" (reprint page 430). **Those two definitions are the whole basis of
this page's equilibrium reconstruction**, so they are quoted here rather than
paraphrased. Two pages later Whitman re-writes the same two potentials in the
station subscripts of his Fig. 1, as $(p_1 - p_3)$ and $(c_1 - c_3)$ in his
eqs. (5) and (6) — i.e. $p_l = p_3$ and $c_g = c_1$. The rest of this page uses
the numbered subscripts throughout; the substitution is his, not ours, and it is
flagged here because it is the one place where a subscript convention carries an
extraction claim.

Whitman's point is that this is a false choice. Both potentials are real and
they act **in series**, on either side of an interface where the two phases are
in equilibrium. Quoting the paper (reprint page 430):

> Conditions at the outside of the gas film (1) are the same as in the main body
> of gas, while those at the inside of the liquid film (3) are the same as in
> the main body of liquid. The gas and liquid at the boundary between the two
> films (2) are in equilibrium. Absorption occurs therefore through two films in
> series.

The two formulations coincide only when $p = kc$; otherwise, he writes, "the
overall coefficients $K_c$ or $K_p$ would have no significance."

**Where the equations on this page come from.** The 1923 original in *Chemical
and Metallurgical Engineering* is pre-DOI and unreachable. It is, however,
**reprinted verbatim** as item 5 of the series *"Pioneer papers in convective
mass transfer"* in *Int. J. Heat Mass Transfer* **5** (1962) 429–433, which
carries its own DOI and is on disk. The reprint's header states its origin:
*"W. G. WHITMAN: The two-film theory of gas absorption, Chemical and
Metallurgical Engineering 29, 146–148 (1923). Reprinted with permission from
Chemical Engineering, Copyright 1923, McGraw-Hill Publishing Co."*

Every equation, every table cell and every constant used here was read off a
**300 ppi render of the 1962 printing, which is that scan's native resolution**
— `pdfimages -list` reports CCITT-G4 bilevel images at 300 ppi, so rendering
larger only interpolates. Each numeric was then **cropped and re-read at that
resolution**, not read at page scale, because the 1962 typesetting sets the
decimal point as a mid-dot (`41·0`, `0·067`) which survives digit-scale
inspection and nothing else. The PDF's text layer is an OCR of the same scan and
is unusable for numbers: it returns eq. (7) as `2,$~=PL~~2Lz` and run 3's $K_c$
as `om7`.

**All page references on this page are to the 1962 reprint, 429–433.** The 1923
range 146–148 appears only in that header line; it is inherited, not verified,
and no individual 1923 page is cited here. The header prints **no issue
number**, so none is given anywhere on this page.

**The Editor's Foreword is not Whitman.** The reprint opens with a 1962
foreword signed "D.B.S." which quotes Sherwood and Pigford's 1951 preface to
*Absorption and Extraction*:

> Methods of applying the Whitman "two-film" theory to various design problems
> have been extended and refined, but it is curious that after 28 years the
> theory itself has never been adequately checked experimentally.

That sentence is editorial commentary about Whitman, not a claim of his, and
this page never attributes it to him. It is quoted because it is the question
this page can answer quantitatively. The Foreword also notes — correctly — that
Whitman knew the film was a fiction: *"the gas and liquid films at the boundary
can be indicated as having a definite thickness, although actually no such sharp
demarcation exists"* is Whitman's own sentence, on reprint page 430.

**The competing pictures are other pages, and the comparison is nobody's yet.**
Higbie's penetration theory (`A3.2`) and Danckwerts' surface renewal (`A3.3`)
replace the steady film with an unsteady one and predict a different exponent on
the diffusivity. That comparison needs all three sources in hand and is not made
here — nor on `A3.3`, which declines it for the same reason, nor on `A3.2`,
which has no source yet. It is open."""))

# ------------------------------------------------------------ published model
cells.append(md(r"""## The published model

Whitman's subscripts are positions in his Fig. 1: **1** the outside of the gas
film (the bulk gas), **2** the interface, **3** the inside of the liquid film
(the bulk liquid). The four equations that matter, transcribed from reprint
pages 431 and 432:

$$\frac{\mathrm{d}W}{\mathrm{d}\theta} \;=\; k_p\,(p_1 - p_2)\;=\;k_c\,(c_2 - c_3)
\tag{4}$$

with $k_p$ "the coefficient of diffusion through the gas film" and $k_c$ that
through the liquid film. The two overall coefficients previously proposed are

$$\frac{\mathrm{d}W}{\mathrm{d}\theta} = K_p\,(p_1-p_3),\tag{5}\qquad\qquad
\frac{\mathrm{d}W}{\mathrm{d}\theta} = K_c\,(c_1-c_3).\tag{6}$$

Dividing eq. (4) by itself and by eq. (5) gives the two working forms:

$$\frac{k_c}{k_p} = \frac{p_1-p_2}{c_2-c_3},\tag{7}\qquad\qquad
\frac{K_p}{k_p} = \frac{p_1-p_2}{p_1-p_3}.\tag{8}$$

**The system closes only with the equilibrium relation.** Eq. (4) is two
equations in three unknowns ($\mathrm{d}W/\mathrm{d}\theta$, $p_2$, $c_2$); what
closes it is $p_2 = f(c_2)$, the equilibrium curve, which Whitman supplies as
Figs 2 and 3. For hydrogen chloride over aqueous HCl at 30 °C that curve is
strongly non-linear — which is the entire point of the paper, because it is
exactly where eqs. (5) and (6) stop meaning anything.

**Read carefully, eqs. (4)–(8) are two different kinds of statement.**
Eqs. (5)–(8) are *definitions and rearrangements*: given $p_2$, the split
$1/K_p = 1/k_p + (p_2-p_3)/(\mathrm{d}W/\mathrm{d}\theta)$ is arithmetic and
cannot be wrong. The **model** is the pair of assertions that $k_p$ and $k_c$
are constants of the apparatus and that the interface is at equilibrium. Those
are what can fail, and this page keeps the two apart everywhere."""))

# --------------------------------------------------------- params/assumptions
cells.append(md(r"""## Parameters and assumptions

Whitman reports no apparatus dimensions, no interfacial area, no flow rates and
no run duration. So:

- $k_p$, $k_c$, $K_p$, $K_c$ are **per-apparatus conductances**, in
  g h⁻¹ mmHg⁻¹ and g h⁻¹ (g/l)⁻¹. They are not per-unit-area film coefficients
  and no film thickness or diffusivity can be recovered from them. The pymrm
  model below therefore uses films of unit thickness so that the conductance
  $D/\delta$ *is* the printed coefficient; the thickness is a gauge choice and
  nothing on this page depends on it. Reconstructing a $\delta$ would be
  fabrication.
- All three runs are at 30 °C, gas and liquid.
- $k_p$ and $k_c$ are assumed the same in all three runs. Whitman states the
  opposite in general — "the values of $k_p$ and $k_c$ will, of course, be
  dependent on experimental conditions" — and then assumes constancy for the
  worked example. Nothing in three runs can check it, and this page does not.
- The interface is at equilibrium, and the equilibrium relation is the one
  plotted in Figs 2 and 3.

**The equilibrium curve is reconstructed, not digitised.** This is the one
reconstruction on the page and it is worth being explicit about, because
digitising Fig. 3 would be the obvious move and it is unnecessary. Whitman's
subscript convention makes $(c_1, p_1)$ and $(c_3, p_3)$ *equilibrium pairs*:
$c_1$ is "the concentration of liquid which would be in equilibrium with the
gas" at $p_1$, and $p_3$ is "the partial pressure of solute exerted by the
liquid" at $c_3$. Table 1 therefore prints **four distinct equilibrium points**,
spanning 2.88 decades of pressure. Fig. 3 plots a straight dashed line on
$\log_{10} p$ against $c$ axes; fitting that form to the four printed points
gives the curve, with every input traceable to a printed number and no pixel
measured anywhere. Run 3's $(9, 0)$ is not usable — the printed zero is a
rounding.

**Table 1 corroborates that reading on its own**, before any fit: runs 2 and 3
share $p_1 = 41$ mmHg and the table prints $c_1 = 368$ g/l for **both**, while
their bulk concentrations, their $K_c$ values and their rates all differ. That
is what an equilibrium concentration must do and what a bulk or inlet
concentration has no reason to do. It is shown below.

The four points are over-determined by two, so the fit **has residuals and can
fail**: a single mis-read digit anywhere in those eight table cells would show
up immediately. That is the page's transcription check.""" ))

# -------------------------------------------------------------------- cell 1
cells.append(code("""# Colab environment cell.
try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm pandas pyyaml matplotlib

import sys, pathlib
_shared = pathlib.Path.cwd()
for _ in range(4):
    if (_shared / "shared" / "gallery_utils.py").is_file():
        sys.path.insert(0, str(_shared / "shared")); break
    _shared = _shared.parent
else:
    if "google.colab" in sys.modules:
        !wget -q https://raw.githubusercontent.com/computational-chemical-engineering/pymrm-gallery/main/shared/gallery_utils.py
        sys.path.insert(0, ".")
import gallery_utils as gu"""))

# ---------------------------------------------------------------- the data
cells.append(md(r"""## The data

Two CSVs, both transcriptions of printed characters read off renders of the 1962
reprint at its native 300 ppi, with every numeric cropped and re-read at that
resolution. Neither is digitised and neither contains a pixel measurement.

- `whitman-1923-table1.csv` — the ten printed columns of Table 1 (reprint
  page 431): three runs, each with the observed rate, the two pressures, the two
  concentrations, the two differences and the two overall coefficients.
- `whitman-1923-printed-results.csv` — everything Whitman prints outside
  Table 1 (reprint pages 431–432): the five numbers of his worked example, three
  ratios he states in words, his "negligible up to approximately 250 g/l", the
  30 °C, and the four printed **axis tick labels** of Fig. 3. The tick labels are
  characters, not positions; no curve or marker is extracted.

**Tier.** Table 1 is a measurement, so it is tier 3 — but a very small one:
three runs, no replicates, no error estimate, no apparatus description, and the
paper's own subtitle is *"A Preliminary Experimental Confirmation"*. Whitman
writes that the data "are insufficient to prove definitely the truth of the
two-film theory". Everything derived from the worked example is tier 6, being
the author's own arithmetic and his own graphical readings.

**What is held out.** The reconstruction below uses only Table 1. Whitman's five
worked-example numbers — $p_2 = 156$, $c_2 = 412$, $k_c = 1.2$, the pair
$(0.6, 224)$, and the predicted 24 g/h — are **never inputs**. They are targets.

This page loads no other page's dataset."""))

cells.append(code("""import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq, least_squares
from pymrm import construct_grad, construct_div, newton

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
np.set_printoptions(linewidth=120)

PAGE = "A3.1-whitman-two-film"
t1 = gu.load_data("whitman-1923-table1.csv", page=PAGE)
pr = gu.load_data("whitman-1923-printed-results.csv", page=PAGE).set_index("quantity")
print(gu.cite_data(gu.load_meta("whitman-1923-table1.csv", page=PAGE)))
display(t1)

def printed(q):
    \"\"\"A number Whitman printed, by name. Nothing in this notebook types one.\"\"\"
    return float(pr.loc[q, "value"])"""))

cells.append(md(r"""### Table 1 is internally consistent — and that check has no power

The two difference columns and the two coefficient columns are recomputed from
the printed rate and the printed pressures and concentrations. This is a
**transcription check only**: $K_p$ and $K_c$ are *defined* as the rate over the
difference, so agreement here is arithmetic, not physics. It is kept because it
is the cheapest way to catch a dropped mid-dot decimal separator — the 1962
typesetting writes 41.0 as `41·0`, and the PDF's own text layer loses it."""))

cells.append(code("""rows = []
for _, r in t1.iterrows():
    rows.append(dict(
        run=int(r["run"]),
        dp_recomputed=r["p1_mmHg"] - r["p3_mmHg"], dp_printed=r["delta_p_mmHg"],
        dc_recomputed=r["c1_g_per_l"] - r["c3_g_per_l"], dc_printed=r["delta_c_g_per_l"],
        Kp_recomputed=r["dW_dtheta_g_per_h"] / r["delta_p_mmHg"], Kp_printed=r["K_p"],
        Kc_recomputed=r["dW_dtheta_g_per_h"] / r["delta_c_g_per_l"], Kc_printed=r["K_c"]))
chk = pd.DataFrame(rows)
chk["dp_dev"] = chk.dp_recomputed - chk.dp_printed
chk["dc_dev"] = chk.dc_recomputed - chk.dc_printed
chk["Kp_reldev_%"] = 100 * (chk.Kp_recomputed / chk.Kp_printed - 1)
chk["Kc_reldev_%"] = 100 * (chk.Kc_recomputed / chk.Kc_printed - 1)
display(chk.round(5))

TABLE1_DIFF_MAX = float(np.abs(np.r_[chk.dp_dev, chk.dc_dev]).max())
TABLE1_COEF_MAX = float(np.abs(np.r_[chk["Kp_reldev_%"], chk["Kc_reldev_%"]]).max())
print(f"difference columns close exactly: worst deviation {TABLE1_DIFF_MAX:.3g}")
print(f"coefficient columns, worst relative deviation: {TABLE1_COEF_MAX:.3f} % "
      "(the printed rounding of K_p and K_c)")"""))

cells.append(md(r"""### Whitman's three stated ratios, recomputed

Three claims are made in words on reprint page 432. Two are right; the third is
loose."""))

cells.append(code("""K = {int(r["run"]): (r["K_p"], r["K_c"]) for _, r in t1.iterrows()}
claims = [
    ("K_p(run 2)/K_p(run 1)", K[2][0]/K[1][0], printed("ratio_Kp_run2_over_run1"),
     "'increases nearly two and one-half fold'"),
    ("K_c(run 1)/K_c(run 3)", K[1][1]/K[3][1], printed("ratio_Kc_run1_over_run3"),
     "'the first run being twelvefold that for run 3'"),
    ("K_c(run 2)/K_c(run 3)", K[2][1]/K[3][1], printed("ratio_Kc_run2_over_run3"),
     "'differ more than twofold' - a lower bound, satisfied"),
]
for name, got, said, quote in claims:
    print(f"  {name} = {got:6.3f}   Whitman states {said:<5g}  {quote}")
RATIO_12FOLD = K[1][1]/K[3][1]
print(f"\\n  The 'twelvefold' is the loose one: the printed columns give "
      f"{RATIO_12FOLD:.2f}, i.e. thirteen, not twelve.")"""))

cells.append(md(r"""### First, the evidence that $c_1$ really is an equilibrium concentration

The reading above is Whitman's own definition, but Table 1 also *tests* it, with
no fit and no figure. If $c_1$ is "the concentration of liquid which would be in
equilibrium with the gas", it is a function of $p_1$ **and of nothing else** —
so two runs at the same $p_1$ must print the same $c_1$, however much else
differs between them. Runs 2 and 3 are exactly that pair."""))

cells.append(code("""shared = t1[t1.p1_mmHg == 41.0]
display(shared[["run", "p1_mmHg", "c1_g_per_l", "c3_g_per_l",
                "dW_dtheta_g_per_h", "K_c"]])
C1_AT_SHARED_P1 = sorted(set(shared["c1_g_per_l"]))
print(f"runs {list(shared.run)} share p1 = 41 mmHg and print c1 = "
      f"{C1_AT_SHARED_P1} g/l - one value, printed twice.")
print(f"Everything else about them differs: c3 = {list(shared.c3_g_per_l)} g/l, "
      f"K_c = {list(shared.K_c)},")
print(f"rate = {list(shared.dW_dtheta_g_per_h)} g/h. A bulk, inlet or mean liquid")
print("concentration would have no reason to coincide; an equilibrium concentration")
print("read off one pressure has no choice.")
C1_SHARED_OK = (len(C1_AT_SHARED_P1) == 1)
print(f"\\nequilibrium-pair reading corroborated by Table 1 alone: {C1_SHARED_OK}")"""))

cells.append(md(r"""### Reconstructing the equilibrium curve from Table 1

Four printed equilibrium pairs, one two-parameter form, two degrees of freedom
left over. The residuals below are what the page's transcription rests on."""))

cells.append(code("""# (c1, p1) and (c3, p3) are equilibrium pairs, by Whitman's own definitions
# (reprint page 430). Run 3's (9, 0) is dropped: the printed 0 is a rounding.
eq_pairs = np.array(
    [[r["c1_g_per_l"], r["p1_mmHg"]] for _, r in t1.iterrows()] +
    [[r["c3_g_per_l"], r["p3_mmHg"]] for _, r in t1.iterrows()])
eq_pairs = np.unique(eq_pairs, axis=0)
eq_pairs = eq_pairs[eq_pairs[:, 1] > 0]            # drop (9, 0)
C_EQ, P_EQ = eq_pairs[:, 0], eq_pairs[:, 1]

def fit_equilibrium(c, p):
    \"\"\"log10 p = A + B c, the straight line Fig. 3 plots.\"\"\"
    B_, A_ = np.polyfit(c, np.log10(p), 1)
    return A_, B_

A_EQ, B_EQ = fit_equilibrium(C_EQ, P_EQ)
peq  = lambda c: 10.0 ** (A_EQ + B_EQ * np.asarray(c, float))
ceq  = lambda p: (np.log10(np.asarray(p, float)) - A_EQ) / B_EQ
dpeq = lambda c: 10.0 ** (A_EQ + B_EQ * np.asarray(c, float)) * np.log(10.0) * B_EQ

res = np.log10(P_EQ) - (A_EQ + B_EQ * C_EQ)
EQ_RESID_MAX_PCT = float(np.abs(10 ** res - 1).max() * 100)
EQ_RESID_MAX_GL  = float(np.abs(res / B_EQ).max())
DECADES = float(np.log10(P_EQ.max() / P_EQ.min()))

print(f"log10(p/mmHg) = {A_EQ:.6f} + {B_EQ:.6f} (c / (g/l))")
for (c_, p_), r_ in zip(eq_pairs, res):
    print(f"   c = {c_:6.1f} g/l   p = {p_:8.3f} mmHg    residual "
          f"{100*(10**r_-1):+7.3f} % in p   {r_/B_EQ:+6.3f} g/l in c")
print(f"\\n  4 points, 2 parameters, 2 degrees of freedom")
print(f"  worst residual {EQ_RESID_MAX_PCT:.3f} % in p  ({EQ_RESID_MAX_GL:.3f} g/l in c) "
      f"over {DECADES:.2f} decades of pressure")"""))

cells.append(md(r"""**One corroboration from Fig. 3's printed axis box, and one from his prose.**

Fig. 3's printed ticks bound a box 280 g/l wide and 3.8 log units tall, and the
dashed line drawn in it **exits at the top-right corner** — visible at a glance,
no measurement involved. Extrapolating the reconstruction to the printed
abscissa ticks should therefore reproduce the printed *top* ordinate tick.

**It cannot do the same at the bottom, and the arithmetic says why before the
number is looked at.** A straight line of the reconstructed slope spans
$B \times 280$ log units across that abscissa, which is *less* than the box is
tall. A line that leaves through the top-right corner must therefore start
somewhere **above** the bottom-left one, by exactly the difference — so the
bottom tick is not a target and a deviation there is not a failure. Both ends
are printed below; only the top one is a check, and the page counts it as one
corroboration, not two.

The second corroboration is Whitman's own sentence that the back pressure is
"negligible up to concentrations of approximately 250 g/l" — the reconstruction
should make it small there, and only there."""))

cells.append(code("""c_lo, c_hi = printed("fig3_abscissa_min"), printed("fig3_abscissa_max")
o_lo, o_hi = printed("fig3_ordinate_min"), printed("fig3_ordinate_max")
lo, hi = A_EQ + B_EQ * c_lo, A_EQ + B_EQ * c_hi
BOX_H = o_hi - o_lo
LINE_SPAN = B_EQ * (c_hi - c_lo)
print("Fig. 3's printed axis box, against the reconstruction extrapolated to it:")
print(f"   at c = {c_hi:.0f} g/l  reconstruction gives log10 p = {hi:+.3f}   "
      f"printed TOP tick    {o_hi:+.1f}   <-- the check")
print(f"   at c = {c_lo:.0f} g/l  reconstruction gives log10 p = {lo:+.3f}   "
      f"printed bottom tick {o_lo:+.1f}   <-- not a target, see below")
AXIS_TOP_DEV = float(abs(hi - o_hi))
AXIS_BOT_DEV = float(abs(lo - o_lo))
print(f"\\n   top-corner deviation    {AXIS_TOP_DEV:.4f} log10 units")
print(f"   bottom-corner deviation {AXIS_BOT_DEV:.4f} log10 units, "
      f"{AXIS_BOT_DEV/AXIS_TOP_DEV:.0f}x larger - and expected:")
print(f"      printed box height       {BOX_H:.3f} log units")
print(f"      line span over 180->460  {LINE_SPAN:.3f} log units")
print(f"      shortfall                {BOX_H-LINE_SPAN:.4f}, against the "
      f"{AXIS_BOT_DEV:.4f} measured at the bottom")
print("   A line of this slope cannot touch both corners of this box. It touches the")
print("   top one, which is where Whitman drew it; the bottom deviation is the")
print("   shortfall plus the top deviation, and carries no independent information.")
print("   (printed tick labels only - no pixel of either figure is measured)")

c_neg = printed("c_negligible_backpressure")
p_neg = float(peq(c_neg))
p1_r1 = float(t1.loc[t1.run == 1, "p1_mmHg"].iloc[0])
p1_r2 = float(t1.loc[t1.run == 2, "p1_mmHg"].iloc[0])
print(f"\\nWhitman's 'negligible up to approximately {c_neg:.0f} g/l':")
print(f"   reconstructed back pressure there = {p_neg:.2f} mmHg")
print(f"   = {100*p_neg/p1_r2:.1f} % of run 2's driving pressure, {100*p_neg/p1_r1:.1f} % of run 1's")"""))

cells.append(code("""fig, ax = plt.subplots(figsize=(7.4, 4.6))
cc = np.linspace(150, 470, 400)
ax.plot(cc, np.log10(peq(cc)), "-", lw=1.6, color="0.35",
        label="reconstruction from Table 1 (4 points, 2 parameters)")
ax.plot(C_EQ, np.log10(P_EQ), "o", ms=8, mfc="none", mew=1.8, color="C0",
        label="equilibrium pairs printed in Table 1")
ax.plot([printed("c2_run1"), printed("c2_run2")],
        np.log10([printed("p2_run1"), printed("p2_run2")]), "s", ms=7, color="C3",
        label="Whitman's own Fig. 3 readings (held out)")
ax.axhline(o_lo, ls=":", lw=1, color="0.6"); ax.axhline(o_hi, ls=":", lw=1, color="0.6")
ax.axvline(c_lo, ls=":", lw=1, color="0.6"); ax.axvline(c_hi, ls=":", lw=1, color="0.6")
ax.set_xlabel("c  /  g HCl per litre"); ax.set_ylabel(r"$\\log_{10}(p\\,/\\,\\mathrm{mmHg})$")
ax.set_title("HCl over aqueous HCl at 30 °C, reconstructed from printed numbers\\n"
             "dotted box = the printed axis ticks of Whitman's Fig. 3")
ax.legend(fontsize=8, loc="upper left"); fig.tight_layout(); plt.show()"""))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

Whitman's "films" are stagnant layers in which the solute moves by diffusion
alone, so the honest pymrm rendering is not the algebra — it is the **two-domain
steady diffusion problem** the algebra summarises:

$$\frac{\mathrm{d}}{\mathrm{d}x}\!\left(-D_g \frac{\mathrm{d}p}{\mathrm{d}x}\right) = 0
\quad\text{on the gas film},\qquad
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(-D_c \frac{\mathrm{d}c}{\mathrm{d}x}\right) = 0
\quad\text{on the liquid film},$$

with $p = p_1$ at the outside of the gas film, $c = c_3$ at the inside of the
liquid film, and at the interface the two conditions that make it a *two-film*
problem: **flux continuity** and **local equilibrium** $p_2 = f(c_2)$.

The state vector is $[\,p_{1..n},\; c_{1..n},\; p_2,\; c_2\,]$. The two interface
values are carried as **external unknowns through `shapes_d`**, exactly as
`AGENTS.md` prescribes: the boundary value changes every Newton step, so it
multiplies a *constant* boundary matrix and the operators are assembled once, in
`__init__`. `construct_div` uses `nu=0` — a plane film, not a curved shell. Both
boundary dictionaries are Dirichlet on the outward normal
($a=0$, $b=1$: $\;b\,\phi = d$, so $\phi = d$ at that end), and the sign
convention never bites because no gradient boundary condition is used.

The system is linear except for the single scalar row $p_2 - f(c_2) = 0$, so the
Jacobian is assembled analytically and `newton` converges in a handful of steps.
`NumJac` is deliberately not used here: with an exact analytic Jacobian available
for a system this small, a numerical one would only add error and cost.

**Film thickness is a gauge.** Both films are given unit thickness so that the
conductance $D/\delta$ equals Whitman's printed coefficient. He publishes no
area, no thickness and no diffusivity, so nothing else is recoverable and
nothing on this page depends on the choice."""))

cells.append(code(r'''def two_film_algebraic(p1, c3, k_p, k_c, pq=None):
    """Whitman's eq. (4) closed with p2 = f(c2), solved as a scalar root.

    No grid, no operators - this is the algebra the paper does by hand, and it
    is what every fit and sweep on this page calls. The bracket is derived from
    the equilibrium relation passed in, so a deliberately perturbed curve (break
    table) still brackets correctly.

    The rate is returned from the GAS side, N = k_p (p1 - p2), never from the
    liquid side N = k_c (c2 - c3). Eq. (4) says the two are equal, but they are
    not equally well conditioned: on this apparatus the gas film carries almost
    the whole drop, so (c2 - c3) is a difference of two nearly equal floats
    while (p1 - p2) is not. At c3 = 204 g/l the spacing of representable
    doubles is 2.8e-14, so at large k_c the liquid-side form quantises the rate
    in steps of k_c * 2.8e-14 - which is why the liquid-side form is not used
    anywhere on this page, and why no large-but-finite k_c has to be chosen.
    The cell below the sweep measures what that choice would have cost.
    """
    pq = peq if pq is None else pq
    f = lambda x: k_c * (x - c3) - k_p * (p1 - pq(x))
    lo, hi, step = c3 + 1e-12, c3 + 1e-12, max(1.0, 0.05 * abs(c3) + 1.0)
    while f(hi) < 0.0 and hi < c3 + 1e7:
        hi += step; step *= 1.6
    c2 = brentq(f, lo, hi)
    return dict(N=k_p * (p1 - float(pq(c2))), p2=float(pq(c2)), c2=c2,
                N_liquid_side=k_c * (c2 - c3))


def no_liquid_film_rate(p1, c3, k_p, pq=None):
    """The k_c -> infinity limit of eq. (4), taken analytically.

    With no liquid-film resistance the interface concentration is the bulk
    concentration, so p2 = f(c3) exactly and the rate is k_p (p1 - f(c3)).
    Nothing is passed a large number and no limit is approached numerically:
    this is the deleted-liquid-film model in closed form, and it is a
    one-parameter model in k_p.
    """
    pq = peq if pq is None else pq
    return dict(N=k_p * (p1 - float(pq(c3))), p2=float(pq(c3)), c2=float(c3))


class TwoFilm:
    """Steady diffusion through a gas film and a liquid film in series.

    Unknowns: p in n gas-film cells, c in n liquid-film cells, and the two
    interface values p2, c2. The interface values enter the films through
    `shapes_d`, so every operator is constant and is built once, here.

    Parameters
    ----------
    k_p, k_c : float
        Film conductances, in Whitman's units g/(h mmHg) and g/(h (g/l)).
        With unit film thickness these are D/delta directly.
    ka : float
        First-order rate constant in the liquid film. ka = 0 is Whitman.
    nu : int
        Geometry index of construct_div. 0 = Cartesian, which is what a plane
        film is. Exposed only so the break table can set it wrongly.
    jump : bool
        True applies the equilibrium condition p2 = f(c2) at the interface.
        False replaces it by p2 = c2, i.e. no equilibrium jump at all - again
        only so the break table can switch it off.
    """

    def __init__(self, k_p, k_c, n=40, ka=0.0, nu=0, jump=True, peq_=None, dpeq_=None):
        self.k_p, self.k_c, self.n, self.ka, self.jump = k_p, k_c, n, ka, jump
        self.peq = peq_ if peq_ is not None else peq
        self.dpeq = dpeq_ if dpeq_ is not None else dpeq
        self.x_f = np.linspace(0.0, 1.0, n + 1)     # unit thickness: D/delta = k
        # Outward normal, a dphi/dn + b phi = d, with a = 0 and b = 1 at both
        # ends: phi = d, i.e. Dirichlet. d = 1 with shapes_d makes the boundary
        # value an external unknown multiplying a CONSTANT matrix.
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
        G, Gl, Gr = construct_grad((n, 1), self.x_f, None, bc, axis=0,
                                   shapes_d=((1, 1), (1, 1)))
        D = construct_div((n, 1), self.x_f, nu=nu, axis=0)   # nu=0: plane film
        self.G, self.Gl, self.Gr, self.D = G, sp.csc_array(Gl), sp.csc_array(Gr), D
        self._assemble()

    def _assemble(self):
        n, G, Gl, Gr, D = self.n, self.G, self.Gl, self.Gr, self.D
        Z, z1 = sp.csc_array((n, n)), sp.csc_array((n, 1))
        zr, z11 = sp.csc_array((1, n)), sp.csc_array((1, 1))
        one = sp.csc_array(np.ones((1, 1)))
        # film rows: div(-k grad phi) = 0, with the interface value as a source
        Ag, Ag_p2, self.bg = -self.k_p * (D @ G), -self.k_p * (D @ Gr), self.k_p * (D @ Gl)
        Al = -self.k_c * (D @ G) + self.ka * sp.eye_array(n)
        Al_c2, self.bl = -self.k_c * (D @ Gl), self.k_c * (D @ Gr)
        # the two interface faces: last face of the gas film, first of the liquid
        e_last = sp.csc_array((np.ones(1), (np.zeros(1), np.array([n]))), shape=(1, n + 1))
        e_frst = sp.csc_array((np.ones(1), (np.zeros(1), np.array([0]))), shape=(1, n + 1))
        Fg, Fg_p2, Fg_p1 = (-self.k_p * (e_last @ G), -self.k_p * (e_last @ Gr),
                            -self.k_p * (e_last @ Gl))
        Fl, Fl_c2, Fl_c3 = (-self.k_c * (e_frst @ G), -self.k_c * (e_frst @ Gl),
                            -self.k_c * (e_frst @ Gr))
        self.M = sp.csc_array(sp.bmat([
            [Ag,   Z,    Ag_p2,  z1   ],     # gas film
            [Z,    Al,   z1,     Al_c2],     # liquid film
            [zr,   zr,   one,    z11  ],     # p2 - f(c2) = 0  (nonlinear part below)
            [Fg,  -Fl,   Fg_p2, -Fl_c2],     # flux continuity across the interface
        ]))
        self.b_p1 = np.r_[self.bg.toarray().ravel(), np.zeros(n), 0.0,
                          Fg_p1.toarray().ravel()]
        self.b_c3 = np.r_[np.zeros(n), self.bl.toarray().ravel(), 0.0,
                          -Fl_c3.toarray().ravel()]

    def _residual(self, y, p1, c3):
        y = np.asarray(y).ravel()
        r = self.M @ y - self.b_p1 * p1 - self.b_c3 * c3
        r[2 * self.n] -= self.peq(y[-1]) if self.jump else y[-1]
        return r

    def _jacobian(self, y):
        J = self.M.tolil()
        J[2 * self.n, 2 * self.n + 1] = -(self.dpeq(np.asarray(y).ravel()[-1])
                                          if self.jump else 1.0)
        return sp.csc_array(J)

    def algebraic(self, p1, c3):
        """Whitman's eq. (4) solved directly, with no grid anywhere."""
        return two_film_algebraic(p1, c3, self.k_p, self.k_c, self.peq)

    def solve(self, p1, c3, tol_rel=1e-11):
        g = self.algebraic(p1, c3) if self.jump else dict(c2=0.5 * (p1 + c3))
        c2g = g["c2"]; p2g = self.peq(c2g) if self.jump else c2g
        y0 = np.r_[np.linspace(p1, p2g, self.n), np.linspace(c2g, c3, self.n), p2g, c2g]
        fun = lambda y: (self._residual(y, p1, c3), self._jacobian(y))
        # `newton` tests ||dy||_inf, which is dimensional: scale the tolerance by
        # the size of the fields so the same tol_rel means the same thing at any
        # grid and any operating point.
        scale = max(abs(p1), abs(c3), 1.0)
        sol = newton(fun, y0.reshape(-1, 1), tol=tol_rel * scale, maxfev=100)
        assert sol.success, "Newton did not converge"        # asserted, never inferred
        y = np.asarray(sol.x).ravel()
        n = self.n
        # one further Newton update, purely as a scale-free convergence measure
        dy = spla.spsolve(self._jacobian(y), self._residual(y, p1, c3))
        upd = float(np.abs(dy).max() / np.abs(y).max())
        assert upd < 1e-10, f"Newton update {upd:.2e} too large"
        return dict(N=self.k_p * (p1 - y[-2]), p2=y[-2], c2=y[-1],
                    p=y[:n], c=y[n:2 * n], nit=sol.nit, update=upd,
                    res=float(np.abs(self._residual(y, p1, c3)).max()))'''))

cells.append(md(r"""### The pymrm solve against the algebra — an identity, and labelled as one

A film with constant $D$, Dirichlet ends and no source has a **linear** profile,
and a conservative finite-volume discretisation of it is exact on any grid. So
the agreement between `TwoFilm.solve` and `TwoFilm.algebraic` below is
**algebraically guaranteed**, and it is reported as a port check: it confirms the
operator assembly, the `shapes_d` plumbing, the flux-continuity row and the
interface equilibrium row are wired together correctly, and nothing else.

Two things follow, and both are measured in the break table further down. It is
**blind to the grid** — three cells and six hundred agree to a few times
$10^{-11}$ relative — because there is no discretisation error to see. And it is
*not* blind to a wrong geometry index or a missing equilibrium jump, which move
it by 144 % and 247 %."""))

cells.append(code("""KP_R3 = float(t1.loc[t1.run == 3, "K_p"].iloc[0])          # Whitman's k_p
KC_W  = printed("k_c_run1")                                # Whitman's k_c = 1.2
CONDS = [(int(r["run"]), r["p1_mmHg"], r["c3_g_per_l"], r["dW_dtheta_g_per_h"])
         for _, r in t1.iterrows()]

tf = TwoFilm(KP_R3, KC_W, n=40)
rows = []
for run, p1, c3, obs in CONDS:
    s, a = tf.solve(p1, c3), tf.algebraic(p1, c3)
    rows.append(dict(run=run, N_pymrm=s["N"], N_algebraic=a["N"],
                     reldev=abs(s["N"] / a["N"] - 1), p2=s["p2"], c2=s["c2"],
                     newton_it=s["nit"], newton_rel_update=s["update"]))
port = pd.DataFrame(rows)
display(port)
PYMRM_VS_ALG = float(port.reldev.max())
NEWTON_UPD_MAX = float(port.newton_rel_update.max())
print(f"worst pymrm-vs-algebra relative deviation : {PYMRM_VS_ALG:.3e}   <-- STRUCTURAL")
print(f"worst converged Newton relative update (asserted, not inferred): {NEWTON_UPD_MAX:.3e}")
print(f"Newton iterations: {list(port.newton_it)}")

grid = [(n, TwoFilm(KP_R3, KC_W, n=n).solve(225.0, 378.0)["N"]) for n in (3, 5, 40, 600)]
print("\\ngrid independence - exact on any grid, because the profile is linear:")
for n, N in grid:
    print(f"   n = {n:4d}   N = {N:.12f}")
GRID_SPREAD = float(max(N for _, N in grid) - min(N for _, N in grid))
print(f"   spread over n = 3 to 600: {GRID_SPREAD:.3e} g/h")"""))

# ------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. Whitman's worked example, reproduced from Table 1 alone

Reprint page 432, in his order:

1. Take $k_p = K_p(\text{run 3})$, "where the liquid is so dilute that the back
   pressure $p_2$ is practically zero".
2. Use run 1 and eq. (4) to get $p_2$ — he prints **156 mm**.
3. Read the matching $c_2$ off Fig. 3 — he prints **412 g/l**.
4. Hence $k_c = 41.0/(412-378)$ — he prints **1.2**.
5. Apply eq. (7) to run 2, solving it together with the equilibrium curve — he
   prints **$p_2 = 0.6$ mm, $c_2 = 224$ g/l**.
6. Hence $\mathrm{d}W/\mathrm{d}\theta = 1.2 \times 20 = $ **24 g/h**, against
   the observed 24.0.

Steps 2 and 4 are arithmetic on Table 1, and reproducing them tests only the
transcription and eq. (4). **Steps 3, 5 and 6 are the held-out ones** — each
needs the equilibrium curve, which Whitman got by eye off Fig. 3 and which here
comes from the reconstruction. All six printed numbers are compared below, but
only the last four are predictions."""))

cells.append(code("""p1_1, c3_1 = 225.0, 378.0
p2_1 = p1_1 - t1.loc[t1.run == 1, "dW_dtheta_g_per_h"].iloc[0] / KP_R3   # eq. (4)
c2_1 = float(ceq(p2_1))
kc_1 = float(t1.loc[t1.run == 1, "dW_dtheta_g_per_h"].iloc[0] / (c2_1 - c3_1))

tf_r = TwoFilm(KP_R3, kc_1, n=40)
r2 = tf_r.solve(41.0, 204.0)                       # run 2, held out entirely

obs2 = float(t1.loc[t1.run == 2, "dW_dtheta_g_per_h"].iloc[0])
comp = pd.DataFrame([
    ("p2, run 1  (mmHg)",  p2_1,      printed("p2_run1"),            "derived by Whitman"),
    ("c2, run 1  (g/l)",   c2_1,      printed("c2_run1"),            "read by Whitman off Fig. 3"),
    ("k_c        (g/h/(g/l))", kc_1,  printed("k_c_run1"),           "derived by Whitman"),
    ("p2, run 2  (mmHg)",  r2["p2"],  printed("p2_run2"),            "read by Whitman off Fig. 3"),
    ("c2, run 2  (g/l)",   r2["c2"],  printed("c2_run2"),            "read by Whitman off Fig. 3"),
    ("dW/dth run 2 (g/h)", r2["N"],   printed("rate_run2_predicted"),
     "Whitman's prediction; equals the OBSERVED 24.0 of Table 1"),
], columns=["quantity", "this page", "Whitman", "what it is"])
comp["deviation"] = comp["this page"] - comp["Whitman"]
comp["reldev_%"] = 100 * (comp["this page"] / comp["Whitman"] - 1)
display(comp.round(4))

WORKED_P2_DEV = float(abs(p2_1 - printed("p2_run1")))
WORKED_C2_DEV = float(abs(c2_1 - printed("c2_run1")))
WORKED_KC_RELDEV = float(abs(kc_1 / printed("k_c_run1") - 1))
RUN2_C2_DEV = float(abs(r2["c2"] - printed("c2_run2")))
RUN2_P2_DEV = float(abs(r2["p2"] - printed("p2_run2")))
RUN2_RATE_RELDEV = float(abs(r2["N"] / obs2 - 1))
print(f"worst deviation on the two Fig. 3 concentrations : "
      f"{max(WORKED_C2_DEV, RUN2_C2_DEV):.2f} g/l")
print(f"held-out run-2 rate: {r2['N']:.3f} g/h against the observed {obs2} "
      f"({100*(r2['N']/obs2-1):+.3f} %)")
print(f"\\nrun 2's interface pressure is the one large relative deviation "
      f"({100*(r2['p2']/printed('p2_run2')-1):+.1f} %): Whitman read 0.6 mmHg where this")
print(f"gives {r2['p2']:.3f}. On his Fig. 3 that is "
      f"{abs(np.log10(r2['p2']/printed('p2_run2'))):.3f} log units on an ordinate spanning "
      f"{o_hi-o_lo:.1f} - a chart-reading difference, and it costs "
      f"{abs(r2['c2']-printed('c2_run2')):.2f} g/l in c2.")"""))

cells.append(md(r"""### 2. How much of that prediction is the liquid film?

The run-2 check is the only prediction in the paper: $k_p$ comes from run 3,
$k_c$ from run 1, and run 2's rate is never used to obtain either. That makes it
the right thing to look at — and the wrong thing to be reassured by.

Run 2 is a **dilute** acid: $c_3 = 204$ g/l, where the back pressure is 0.3 mmHg
against a bulk gas pressure of 41. The interface barely backs up, so almost the
whole resistance sits in the gas film. The table below fixes $k_p$ and sweeps
$k_c$, ending with $k_c \to \infty$ — **deleting the liquid film entirely**.

That last row is taken **analytically**, not by putting a large number in
$k_c$. With no liquid resistance the interface concentration *is* the bulk
concentration, so the limit is the closed form $k_p\,(p_1 - f(c_3))$ with no
tuning constant in it. The cell after the sweep shows what the numerical route
would have cost, because it is the more interesting number of the two."""))

cells.append(code("""def run2_rate(kc, kp=KP_R3):
    return TwoFilm(kp, kc, n=12).algebraic(41.0, 204.0)

base = run2_rate(kc_1)
share_gas_2 = (41.0 - base["p2"]) / (41.0 - 0.3)
r1 = TwoFilm(KP_R3, kc_1, n=12).algebraic(225.0, 378.0)
share_gas_1 = (225.0 - r1["p2"]) / (225.0 - 55.0)
print(f"gas-film share of the total resistance:  run 2 {100*share_gas_2:.2f} %"
      f"   run 1 {100*share_gas_1:.2f} %")

NO_LF = no_liquid_film_rate(41.0, 204.0, KP_R3)      # the k_c -> infinity limit
rows = []
for lab, kc in [("k_c / 8", kc_1/8), ("k_c / 4", kc_1/4), ("k_c / 2", kc_1/2),
                ("k_c as fitted", kc_1), ("k_c x 2", 2*kc_1), ("k_c x 10", 10*kc_1),
                ("k_c x 1000", 1000*kc_1)]:
    s = run2_rate(kc)
    rows.append(dict(case=lab, k_c=kc, c2=s["c2"], p2=s["p2"], rate=s["N"],
                     dev_vs_observed_pct=100*(s["N"]/obs2 - 1)))
rows.append(dict(case="k_c -> infinity, ANALYTIC (NO liquid film)", k_c=np.inf,
                 c2=NO_LF["c2"], p2=NO_LF["p2"], rate=NO_LF["N"],
                 dev_vs_observed_pct=100*(NO_LF["N"]/obs2 - 1)))
power = pd.DataFrame(rows)
display(power.round(6))

NO_LIQUID_RELDEV = float(abs(NO_LF["N"]/obs2 - 1))
NO_LF_DEV_PCT_ABS = 100*NO_LIQUID_RELDEV
print(f"\\nWith the liquid film DELETED, run 2's predicted rate is "
      f"{NO_LF['N']:.6f} g/h, {100*(NO_LF['N']/obs2-1):+.4f} % from the observed {obs2}.")
print(f"With Whitman's k_c = {kc_1:.3f} it is {base['N']:.6f} g/h, "
      f"{100*(base['N']/obs2-1):+.4f} %.")
print(f"Table 1 prints the rate to 0.1 g/h, i.e. +-{100*0.05/obs2:.2f} %.")
print("\\nSo the paper's one held-out prediction cannot distinguish a liquid film")
print("from no liquid film, and of the two the NO-liquid-film answer is closer.")"""))

cells.append(md(r"""#### The same limit, computed a second and independent way

Every other check on this page perturbs an input and watches a number move.
That establishes *sensitivity*; it cannot establish that a baseline is right,
and the number above is the page's headline. So it is computed twice, by two
routes that share no arithmetic.

- **Route A**, used above: the closed form $k_p\,(p_1 - f(c_3))$. One
  multiplication; the root solve is never called.
- **Route B**: solve the *full* two-film problem at two finite, well-conditioned
  values of $k_c$ and extrapolate. The leading correction is
  $N(k_c) = N_\infty - a/k_c + O(k_c^{-2})$, because raising the interface by
  $\delta c = N/k_c$ raises the back pressure by $f'(c_3)\,\delta c$, so two
  points a decade apart eliminate it:
  $N_\infty \approx \bigl(10\,N(10k) - N(k)\bigr)/9$. Route B calls `brentq`,
  never evaluates the closed form, and never uses a number larger than
  $10^5 k_c$.

And the third cell shows what happens if the limit is taken the obvious way
instead — by setting $k_c$ to a large constant and reading the rate off the
**liquid** side of eq. (4)."""))

cells.append(code("""kA = kc_1 * 1e4
NB = (10*run2_rate(10*kA)["N"] - run2_rate(kA)["N"]) / 9.0      # route B
NO_LF_TWO_ROUTES = float(abs(NB/NO_LF["N"] - 1))
print("run 2 with the liquid film deleted, two independent routes:")
print(f"   A  closed form  k_p (p1 - f(c3))          {NO_LF['N']:.12f} g/h")
print(f"   B  Richardson on finite k_c (1e4, 1e5 x)  {NB:.12f} g/h")
print(f"   relative difference                       {NO_LF_TWO_ROUTES:.3e}")
print(f"   deviation from the observed {obs2} g/h:   "
      f"{100*(NO_LF['N']/obs2-1):+.4f} % (A)   {100*(NB/obs2-1):+.4f} % (B)")"""))

cells.append(md(r"""#### And why the obvious route is wrong, which is worth printing

Eq. (4) says $k_p(p_1-p_2) = k_c(c_2-c_3)$, so it looks as though either side
may be evaluated. They are not equally conditioned. As $k_c$ grows, $c_2-c_3$
shrinks like $1/k_c$ while $c_2$ itself stays near 204 g/l, where consecutive
representable doubles are $\mathrm{ulp}(204)$ apart. The liquid-side product
$k_c(c_2-c_3)$ therefore comes in steps of $k_c\,\mathrm{ulp}(204)$ — a
quantisation that grows with $k_c$ while the quantity it is quantising does not.
The table prints the step beside the answer.

The signature is visible without knowing any of that: the liquid-side value
**exceeds** the exact $k_c\to\infty$ limit, which a rate increasing in $k_c$
cannot do."""))

cells.append(code("""ULP = float(np.spacing(204.0))
rows = []
for kc in [1e9, 1e10, 1e11, 1e12, 1e13]:
    s = TwoFilm(KP_R3, kc, n=12).algebraic(41.0, 204.0)
    rows.append(dict(k_c=kc, c2_minus_c3=s["c2"] - 204.0,
                     quantisation_g_per_h=kc*ULP,
                     rate_gas_side=s["N"], rate_liquid_side=s["N_liquid_side"],
                     dev_gas_pct=100*(s["N"]/obs2 - 1),
                     dev_liquid_pct=100*(s["N_liquid_side"]/obs2 - 1),
                     exceeds_the_limit=bool(s["N_liquid_side"] > NO_LF["N"])))
cond = pd.DataFrame(rows)
display(cond)
print(f"exact k_c -> infinity limit: {NO_LF['N']:.6f} g/h;  ulp(204) = {ULP:.3e} g/l")
print(f"the quantity being reported is {abs(NO_LF['N']-obs2):.4f} g/h "
      f"({NO_LF_DEV_PCT_ABS:.4f} % of {obs2}),")
print(f"against a liquid-side quantisation of {1e12*ULP:.4f} g/h at k_c = 1e12 - "
      f"a factor {1e12*ULP/abs(NO_LF['N']-obs2):.1f} larger.")
print("\\nThe gas-side column is flat to twelve figures over four decades of k_c.")
print("The liquid-side column drifts, crosses its own limit, and at 1e13 reverses")
print("the sign of the page's headline. No value of k_c is load-bearing above,")
print("because no value of k_c is used above.")"""))

cells.append(code("""# The classical two-film diagram: ONE potential across both films. In the gas
# film that is the partial pressure itself; in the liquid film it is the
# pressure the local concentration would exert, p_eq(c). The two meet at the
# interface by construction, so the drop across each film is directly readable.
fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
for ax, (run, p1, c3, obs), ttl in zip(
        axes, [CONDS[0], CONDS[1]],
        ["run 1: strong acid, $c_3$ = 378 g/l", "run 2: dilute acid, $c_3$ = 204 g/l"]):
    s = TwoFilm(KP_R3, kc_1, n=60).solve(p1, c3)
    xg = np.linspace(0, 1, 60, endpoint=False) + 1 / 120
    ax.plot(xg, s["p"] / p1, "-", lw=2.2, color="C0", label="gas film: $p$")
    ax.plot(1 + xg, peq(s["c"]) / p1, "-", lw=2.2, color="C3",
            label=r"liquid film: $p_{\\mathrm{eq}}(c)$")
    ax.axvline(1.0, color="0.4", lw=1.2)
    ax.axhline(float(peq(c3)) / p1, ls=":", lw=1, color="0.5")
    ax.set_xlabel("gas film                    |                    liquid film")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{ttl}\\ngas film carries "
                 f"{100*(p1-s['p2'])/(p1-float(peq(c3))):.1f} % of the resistance", fontsize=10)
    ax.set_xticks([]); ax.legend(fontsize=8, loc="lower left")
axes[0].set_ylabel(r"potential  $p\\,/\\,p_1$")
fig.suptitle("Where the resistance sits, at Whitman's own two extremes", y=1.02)
fig.tight_layout(); plt.show()"""))

cells.append(md(r"""### 3. The three runs against three models

This is the comparison with power. Three published rates; three candidate
models, each fitted to all three runs by least squares on the *relative* residual
in $\mathrm{d}W/\mathrm{d}\theta$:

- **single gas film**, eq. (5) with a constant $K_p = k_p$ — one parameter;
- **single liquid film**, eq. (6) with a constant $K_c = k_c$ — one parameter;
- **two films**, eq. (4) closed with the reconstructed equilibrium curve — two
  parameters, so **one degree of freedom** remains.

**Two of these three numbers were settled before any physics entered, and the
page should say which.** The two-film family *nests* the single-gas-film model:
$k_c \to \infty$ is exactly eq. (5) with the back pressure taken from the
equilibrium curve instead of from Table 1's $p_3$ column. Its second parameter
is then free to absorb run 1 — the only run with a real liquid resistance — and
it does, to a residual of a few thousandths of a per cent. So the three-run rms
is set entirely by runs 2 and 3, and their mutual consistency under *one* gas
coefficient is already printed in Table 1: $K_p = 0.59$ for both. The two-film
rms could not have come out at tens of per cent.

What could fail is the one degree of freedom that is left over, and section 4
shows that it does — and that the liquid film makes it *worse*."""))

cells.append(code("""RUNS = t1.to_dict("records")

def rel_resid_two(theta, kc=None, subset=(1, 2, 3), runs=None, pq=None):
    kp = np.exp(theta[0]); kcv = np.exp(theta[1]) if kc is None else kc
    out = []
    for r in (RUNS if runs is None else runs):
        if int(r["run"]) not in subset: continue
        m = two_film_algebraic(r["p1_mmHg"], r["c3_g_per_l"], kp, kcv, pq)
        out.append(m["N"] / r["dW_dtheta_g_per_h"] - 1)
    return out

def rel_resid_gas(theta, runs=None, subset=(1, 2, 3)):
    kp = np.exp(theta[0])
    return [kp*(r["p1_mmHg"]-r["p3_mmHg"])/r["dW_dtheta_g_per_h"] - 1
            for r in (RUNS if runs is None else runs) if int(r["run"]) in subset]

def rel_resid_nolf(theta, runs=None, subset=(1, 2, 3), pq=None):
    \"\"\"The two-film family's k_c -> infinity member: one parameter, k_p.\"\"\"
    kp = np.exp(theta[0])
    return [no_liquid_film_rate(r["p1_mmHg"], r["c3_g_per_l"], kp, pq)["N"]
            / r["dW_dtheta_g_per_h"] - 1
            for r in (RUNS if runs is None else runs) if int(r["run"]) in subset]

def rel_resid_liq(theta, runs=None):
    kc = np.exp(theta[0])
    return [kc*(r["c1_g_per_l"]-r["c3_g_per_l"])/r["dW_dtheta_g_per_h"]
            - 1 for r in (RUNS if runs is None else runs)]

def fit(f, x0):
    s = least_squares(f, x0); r = np.asarray(f(s.x))
    return np.exp(s.x), r, float(np.sqrt((r**2).mean()))

par_g, res_g, rms_g = fit(rel_resid_gas, [np.log(0.4)])
par_l, res_l, rms_l = fit(rel_resid_liq, [np.log(0.3)])
par_2, res_2, rms_2 = fit(lambda th: rel_resid_two(th), np.log([0.59, 1.2]))

mc = pd.DataFrame([
    dict(model="single gas film, eq. (5)", n_par=1,
         fitted=f"k_p = {par_g[0]:.4f}", rms_pct=100*rms_g,
         worst_pct=100*np.abs(res_g).max()),
    dict(model="single liquid film, eq. (6)", n_par=1,
         fitted=f"k_c = {par_l[0]:.4f}", rms_pct=100*rms_l,
         worst_pct=100*np.abs(res_l).max()),
    dict(model="two films, eq. (4)", n_par=2,
         fitted=f"k_p = {par_2[0]:.4f}, k_c = {par_2[1]:.4f}", rms_pct=100*rms_2,
         worst_pct=100*np.abs(res_2).max()),
])
display(mc.round(3))
print("per-run relative residuals, %")
print(f"   single gas film    {np.round(100*res_g, 2)}")
print(f"   single liquid film {np.round(100*res_l, 2)}")
print(f"   two films          {np.round(100*res_2, 3)}")
RMS_GAS, RMS_LIQ, RMS_TWO = 100*rms_g, 100*rms_l, 100*rms_2
print(f"\\nTable 1's own precision on a rate is +-{100*0.05/24.0:.2f} %.")
print("A single film - either one - misses by tens of per cent. Two films do not.")"""))

cells.append(md(r"""#### How much of the 0.38 % was available to fail

Run 1's two-film residual is the smallest of the three, because $k_c$ is the
parameter that fits it. Strip run 1 out and what remains is a **two-point,
one-parameter** problem, whose residual is fixed by how well runs 2 and 3 agree
under a single gas coefficient — a number Table 1 prints directly, as $K_p$
twice."""))

cells.append(code("""Kp2 = obs2/float(t1.loc[t1.run == 2, "delta_p_mmHg"].iloc[0])
Kp3 = float(t1.loc[t1.run == 3, "dW_dtheta_g_per_h"].iloc[0]) / \\
      float(t1.loc[t1.run == 3, "delta_p_mmHg"].iloc[0])
RUNS23_INCONSISTENCY_PCT = 100*abs(Kp2/Kp3 - 1)
print(f"Table 1's K_p column, recomputed: run 2 {Kp2:.5f}, run 3 {Kp3:.5f} "
      f"- both printed as 0.59")
print(f"   they disagree by {RUNS23_INCONSISTENCY_PCT:.3f} %, and that is the whole")
print("   budget any one-parameter gas model has to work with on these two runs.\\n")

par_g23, res_g23, rms_g23 = fit(lambda th: rel_resid_gas(th, subset=(2, 3)), [np.log(0.59)])
par_n23, res_n23, rms_n23 = fit(lambda th: rel_resid_nolf(th, subset=(2, 3)), [np.log(0.59)])
res_2_23 = np.asarray(res_2)[[1, 2]]
GASFILM_RUNS23_RMS = 100*rms_g23
TWOFILM_RUNS23_RMS = float(np.sqrt((res_2_23**2).mean())*100)
cmp23 = pd.DataFrame([
    dict(model="single gas film, eq. (5), printed p3", n_par=1,
         fitted=f"k_p = {par_g23[0]:.5f}", rms_pct=100*rms_g23,
         resid_run2_pct=100*res_g23[0], resid_run3_pct=100*res_g23[1]),
    dict(model="liquid film DELETED (k_c -> inf), f(c3)", n_par=1,
         fitted=f"k_p = {par_n23[0]:.5f}", rms_pct=100*rms_n23,
         resid_run2_pct=100*res_n23[0], resid_run3_pct=100*res_n23[1]),
    dict(model="two films (k_p, k_c fitted to all three)", n_par=2,
         fitted=f"k_p = {par_2[0]:.5f}, k_c = {par_2[1]:.5f}",
         rms_pct=TWOFILM_RUNS23_RMS,
         resid_run2_pct=100*res_2_23[0], resid_run3_pct=100*res_2_23[1]),
], columns=["model", "n_par", "fitted", "rms_pct", "resid_run2_pct", "resid_run3_pct"])
display(cmp23.round(4))
print(f"On the only two runs that constrain k_p, a ONE-parameter gas film fits to "
      f"{100*rms_g23:.3f} % rms")
print(f"and the two-parameter two-film model to {TWOFILM_RUNS23_RMS:.3f} % - "
      f"a factor {TWOFILM_RUNS23_RMS/(100*rms_g23):.2f} WORSE.")
print("Adding the liquid film does not improve these two runs; it degrades them,")
print("because fitting run 1 pushes run 2's interface from 0.30 up to "
      f"{r2['p2']:.3f} mmHg.")"""))

cells.append(md(r"""### 4. Where the one degree of freedom lands, and it does not land cleanly

Three observations minus two parameters leaves one testable combination, and it
is essentially the **ratio of run 3's rate to run 2's**. Both runs share
$p_1 = 41$ mmHg; they differ only in bulk liquid concentration, 9 against
204 g/l. Two-film theory says run 3 must be the faster of the two, because run
2's interface backs up to a finite partial pressure and run 3's does not. The
question is by how much."""))

cells.append(code("""kp_f, kc_f = par_2
m2 = TwoFilm(kp_f, kc_f, n=6).algebraic(41.0, 204.0)
m3 = TwoFilm(kp_f, kc_f, n=6).algebraic(41.0, 9.0)
obs3 = float(t1.loc[t1.run == 3, "dW_dtheta_g_per_h"].iloc[0])
ratio_model, ratio_obs = m3["N"]/m2["N"], obs3/obs2
DOF_OVERPRED_PCT = 100*(ratio_model/ratio_obs - 1)
prec_ratio = 100*np.hypot(0.05/obs2, 0.05/obs3)
print(f"  run 2 interface backs up to p2 = {m2['p2']:.4f} mmHg; run 3 to {m3['p2']:.5f} mmHg")
print(f"  two-film model  run3/run2 = {ratio_model:.5f}   ({100*(ratio_model-1):+.2f} %)")
print(f"  Table 1         run3/run2 = {ratio_obs:.5f}   ({100*(ratio_obs-1):+.2f} %)")
print(f"  the model OVER-predicts the gap by {DOF_OVERPRED_PCT:+.3f} %")
print(f"  Table 1's printed precision on that ratio is +-{prec_ratio:.2f} %")
print("\\n  and the ratio barely depends on k_c, so this test does not probe it either:")
for lab, k in [(f"{kc_f:9.4g}", kc_f), (f"{5*kc_f:9.4g}", 5*kc_f), ("infinity ", None)]:
    if k is None:
        rr = (no_liquid_film_rate(41.0, 9.0, kp_f)["N"] /
              no_liquid_film_rate(41.0, 204.0, kp_f)["N"])
    else:
        rr = (TwoFilm(kp_f, k, n=6).algebraic(41.0, 9.0)["N"] /
              TwoFilm(kp_f, k, n=6).algebraic(41.0, 204.0)["N"])
    print(f"     k_c = {lab}   run3/run2 = {rr:.5f}")"""))

cells.append(md(r"""#### Four models on the one testable quantity

$k_p$ cancels out of a ratio of two runs at the same $p_1$, so three of these
four predictions contain **no fitted parameter at all** — they are properties of
the model, not of the fit. That makes this the cleanest comparison the data
allow, and it puts the two-film model last."""))

cells.append(code("""p3_2 = float(t1.loc[t1.run == 2, "p3_mmHg"].iloc[0])
p3_3 = float(t1.loc[t1.run == 3, "p3_mmHg"].iloc[0])
ratio_rows = [
    ("no back pressure at all", 1.0, "no parameters"),
    ("single gas film, eq. (5), Table 1's p3",
     (41.0 - p3_3)/(41.0 - p3_2), "no parameters"),
    ("liquid film deleted, back pressure from f(c)",
     no_liquid_film_rate(41.0, 9.0, kp_f)["N"]/no_liquid_film_rate(41.0, 204.0, kp_f)["N"],
     "no parameters (k_p cancels)"),
    ("two films", ratio_model, "k_c enters; k_p cancels"),
]
rt = pd.DataFrame([dict(model=m, run3_over_run2=v,
                        dev_vs_table1_pct=100*(v/ratio_obs - 1), free=f)
                   for m, v, f in ratio_rows])
display(rt.round(5))
GASFILM_DOF_PCT = float(rt.dev_vs_table1_pct.iloc[1])
print(f"Table 1:  run3/run2 = {ratio_obs:.5f}   printed precision "
      f"+-{prec_ratio:.2f} %")
print(f"The observed gap sits BETWEEN no back pressure ({rt.dev_vs_table1_pct.iloc[0]:+.3f} %)")
print(f"and a gas film alone ({GASFILM_DOF_PCT:+.3f} %). Adding the liquid film moves the")
print(f"prediction away from it, to {DOF_OVERPRED_PCT:+.3f} %.")
print(f"\\nNote {GASFILM_DOF_PCT:+.3f} % is the same {RUNS23_INCONSISTENCY_PCT:.3f} % as the")
print("disagreement of Table 1's two printed K_p = 0.59 entries - it has to be, since")
print("K_p(run 2)/K_p(run 3) = (N2/N3)(dp3/dp2). One number, seen twice.")"""))

cells.append(md(r"""So the single quantitative test the data offer is **not passed cleanly**: the
observed gap between runs 2 and 3 is about a third of the predicted one, a
discrepancy of order 1 % against a printed precision of order 0.3 %. It is far
too small to refute anything — Whitman quotes no experimental error, and 1 % on a
1923 absorption rate is nothing — but it is worth being exact about what it says.
On **both** quantities the three runs can actually test — the held-out run-2 rate
of section 2, and this ratio — the two-film model is beaten by the same model
with the liquid film taken out, which has one parameter fewer. That is the same
finding twice, and it points the same way as everything else on this page:
**the data show less liquid-film effect than the theory needs.**

### 5. What the three runs can say about $k_c$, and what they cannot

Profile the fit: fix $k_c$, re-fit $k_p$, and plot the residual. Do it twice —
once on all three runs, once with run 1 removed."""))

cells.append(code("""KC_SWEEP_TOP = 3e3
kcs = np.unique(np.r_[np.geomspace(0.3, KC_SWEEP_TOP, 220), kc_f * np.linspace(0.9, 1.1, 21)])
SUBSETS = [((1, 2, 3), "all three runs"), ((1, 2), "runs 1 and 2"),
           ((1, 3), "runs 1 and 3"), ((2, 3), "runs 2 and 3 - run 1 REMOVED")]
prof, bands = {}, {}
THRESH = 0.02      # 2 % rms, about ten times Table 1's own precision on a rate
for sub, lab in SUBSETS:
    v = []
    for kc in kcs:
        s = least_squares(lambda th: rel_resid_two(th, kc=kc, subset=sub), [np.log(0.59)])
        v.append(np.sqrt((np.asarray(rel_resid_two(s.x, kc=kc, subset=sub))**2).mean()))
    prof[lab] = v = np.asarray(v)
    ok = kcs[v <= THRESH]
    bands[lab] = (ok.min(), ok.max()) if len(ok) else (np.nan, np.nan)
    censored = len(ok) and ok.max() >= 0.99 * KC_SWEEP_TOP
    print(f"  {lab:30s}: min rms {100*v.min():7.3f} % at k_c = {kcs[v.argmin()]:8.3f};  "
          f"k_c admitted at rms <= {100*THRESH:.0f} %: "
          f"[{bands[lab][0]:.3g}, {bands[lab][1]:.4g}]"
          f"{'  UNBOUNDED ABOVE (censored at the sweep top)' if censored else ''}")
KC_BAND_ALL = bands["all three runs"]
KC_BAND_NO1 = bands["runs 2 and 3 - run 1 REMOVED"]
KC_BAND_WIDTH_ALL = float(KC_BAND_ALL[1] / KC_BAND_ALL[0])
KC_BAND_WIDTH_NO1 = float(KC_BAND_NO1[1] / KC_BAND_NO1[0])
print(f"\\n  band width with all three runs: x {KC_BAND_WIDTH_ALL:.3f}")
print(f"  band width with run 1 removed : x {KC_BAND_WIDTH_NO1:.4g} and still rising "
      "at the top of the sweep")

fig, ax = plt.subplots(figsize=(7.6, 4.4))
for (sub, lab), st in zip(SUBSETS, ["-", "--", "-.", "-"]):
    ax.loglog(kcs, 100*prof[lab], st, lw=2.2 if len(sub) != 2 else 1.4, label=lab)
ax.axhline(100*THRESH, ls=":", lw=1.2, color="0.4", label=f"{100*THRESH:.0f} % rms threshold")
ax.axvline(kc_f, ls=(0, (6, 3)), lw=1.4, color="0.15",
           label=f"Whitman's $k_c$ = {kc_f:.2f}")
ax.set_xlabel(r"$k_c$  /  g h$^{-1}$ (g/l)$^{-1}$   (held fixed; $k_p$ re-fitted at each point)")
ax.set_ylabel(r"rms relative residual in $dW/d\\theta$  /  %")
ax.set_title("Every subset containing run 1 pins $k_c$; the one without it does not")
ax.legend(fontsize=8); fig.tight_layout(); plt.show()"""))

cells.append(md(r"""**Run 1 is the only run that carries information about $k_c$.** Every subset
that contains it pins the coefficient to within a few per cent; the subset that
does not runs to the top of the sweep, so the remaining data are consistent with
an infinitely fast liquid film.

That distinction matters and it is easy to blur. $k_c$ is tightly *identified* —
one equation, one unknown, and a narrow band. It is not *tested*: the narrowness
is entirely run 1's doing, and run 1 is precisely the run Whitman fits $k_c$ to.
So the liquid-film coefficient of the two-film theory is, in this paper,
determined by one observation and corroborated by none.

That is a quantitative form of the remark the 1962 Foreword quotes from Sherwood
and Pigford. It is not a criticism of Whitman, who says the same thing in
plainer words across reprint pages 432–433 — the data "are insufficient to prove
definitely the truth of the two-film theory" — and calls his own paper "A
Preliminary Experimental Confirmation".

What the data *do* establish, and it is not nothing: **a single film of either
kind cannot describe these three runs**, by tens of per cent, and a two-film
description can. The qualitative claim survives; the quantitative liquid-film
coefficient is untested."""))

cells.append(md(r"""### 6. Where resistance additivity breaks: the textbook formula

The form everyone actually uses is not eq. (4) but its linearised child,

$$\frac{1}{K_p} = \frac{1}{k_p} + \frac{1}{m\,k_c},\qquad m = \frac{\mathrm{d}p}{\mathrm{d}c},$$

which requires a **single** $m$ — Henry's law. Whitman says this is exactly what
fails: "if the deviation from Henry's law … is considerable, a simplification
based on direct proportionality between $p$ and $c$ would be unjustified and the
overall coefficients $K_c$ or $K_p$ would have no significance."

Below, that formula is evaluated with the two most defensible choices of $m$ —
the tangent slope at the bulk liquid concentration, and the secant across the
run's own concentration span — and compared with the exact non-linear interface
solve at the same $k_p$, $k_c$."""))

cells.append(code("""rows = []
for r in RUNS:
    p1, p3, c1, c3 = r["p1_mmHg"], r["p3_mmHg"], r["c1_g_per_l"], r["c3_g_per_l"]
    ex = TwoFilm(KP_R3, kc_1, n=6).algebraic(p1, c3)
    K_exact = ex["N"] / (p1 - float(peq(c3)))
    m_tan, m_sec = float(dpeq(c3)), (p1 - p3) / (c1 - c3)
    K_tan = 1.0 / (1.0/KP_R3 + 1.0/(m_tan*kc_1))
    K_sec = 1.0 / (1.0/KP_R3 + 1.0/(m_sec*kc_1))
    rows.append(dict(run=int(r["run"]), K_p_exact=K_exact,
                     m_tangent=m_tan, K_p_tangent=K_tan, err_tangent_pct=100*(K_tan/K_exact-1),
                     m_secant=m_sec, K_p_secant=K_sec, err_secant_pct=100*(K_sec/K_exact-1)))
add = pd.DataFrame(rows)
display(add.round(4))
CONSTH_R1_TAN = float(abs(add.loc[add.run == 1, "err_tangent_pct"].iloc[0]))
CONSTH_R1_SEC = float(abs(add.loc[add.run == 1, "err_secant_pct"].iloc[0]))
CONSTH_SEC_MAX = float(np.abs(add.err_secant_pct).max())
print(f"run 1 - the only run with a real liquid resistance - the formula is off by "
      f"{CONSTH_R1_TAN:.1f} % (tangent m) and {CONSTH_R1_SEC:.1f} % (secant m).")
print(f"worst over all three runs, secant m: {CONSTH_SEC_MAX:.1f} %")
print("Runs 2 and 3 are gas-controlled and their tangent errors saturate at -100 %:")
print("  the tangent slope at c_3 = 204 and 9 g/l is 0.009 and 0.000 mmHg/(g/l), so")
print("  the formula predicts essentially zero transfer. That number cannot get worse,")
print("  which is why the two run-1 errors are the ones reported as metrics.")
print(f"\\nThe equilibrium curve rises by a factor 10 every "
      f"{1/B_EQ:.1f} g/l, so no single m survives a driving force of 47 g/l.")"""))

cells.append(md(r"""### 7. The overall coefficient is not a property of the apparatus

Same $k_p$, same $k_c$, same equipment. Only the operating point moves."""))

cells.append(code("""c3s = np.linspace(5.0, 415.0, 90)
p1s = [225.0, 120.0, 41.0]
fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
Kp_all, Kc_all = [], []
for p1 in p1s:
    Kp, Kc = [], []
    for c3 in c3s:
        if float(peq(c3)) >= p1 - 1e-9:
            Kp.append(np.nan); Kc.append(np.nan); continue
        s = TwoFilm(KP_R3, kc_1, n=6).algebraic(p1, c3)
        Kp.append(s["N"] / (p1 - float(peq(c3)))); Kc.append(s["N"] / (float(ceq(p1)) - c3))
    axes[0].plot(c3s, Kp, lw=1.8, label=f"$p_1$ = {p1:.0f} mmHg")
    axes[1].plot(c3s, Kc, lw=1.8, label=f"$p_1$ = {p1:.0f} mmHg")
    Kp_all += [v for v in Kp if np.isfinite(v)]
    Kc_all += [v for v in Kc if np.isfinite(v)]
for _, r in t1.iterrows():
    axes[0].plot(r["c3_g_per_l"], r["K_p"], "ks", ms=7, zorder=5)
    axes[1].plot(r["c3_g_per_l"], r["K_c"], "ks", ms=7, zorder=5)
axes[0].set_ylabel(r"$K_p$  /  g h$^{-1}$ mmHg$^{-1}$")
axes[1].set_ylabel(r"$K_c$  /  g h$^{-1}$ (g/l)$^{-1}$")
for ax in axes:
    ax.set_xlabel(r"bulk liquid concentration $c_3$  /  g l$^{-1}$"); ax.legend(fontsize=8)
axes[0].set_title("black squares: Whitman's three runs")
fig.suptitle(r"One apparatus, one $k_p$, one $k_c$ — the 'overall coefficients' still move", y=1.02)
fig.tight_layout(); plt.show()

KP_RANGE_FACTOR = float(max(Kp_all) / min(Kp_all))
KC_RANGE_FACTOR = float(max(Kc_all) / min(Kc_all))
print(f"At fixed k_p = {KP_R3} and k_c = {kc_1:.3f}, across the operating plane "
      "swept above:")
print(f"   K_p ranges over a factor {KP_RANGE_FACTOR:.2f} "
      f"({min(Kp_all):.4f} to {max(Kp_all):.4f})")
print(f"   K_c ranges over a factor {KC_RANGE_FACTOR:.2f} "
      f"({min(Kc_all):.5f} to {max(Kc_all):.5f})")"""))

# ---------------------------------------------------------------- validation
cells.append(md(r"""### 8. One thing the film model can do that the algebra cannot

The film is a *differential* object, and writing it as one buys a limit the
algebra has no way to express. Put a first-order reaction into the liquid film
and let the rate constant grow: the liquid-side resistance must vanish, and the
rate must rise towards $k_p\,p_1$ — the gas film alone — from below, whatever the
equilibrium curve does. That limit is not built into the code anywhere.

**The residual is a grid floor, not a convergence rate, and the page says so.**
The remaining gap is exactly $p_2/p_1$, and $p_2$ stops falling once the reaction
layer — of thickness $\sqrt{D/k a}$, i.e. about a thousandth of the film at
$ka = 10^6$ — is finer than the mesh. So the last two rows of the table below
**saturate** rather than converge, at a value set by $n$. The check is real:
removing the equilibrium jump moves it by two and a half orders of magnitude,
and the sweep in $n$ shows what the floor is made of. But it is a floor.

This is an illustration of the mechanism, not a result about HCl: Whitman
publishes no rate constant, and the enhancement regimes belong to
[`F3.1`](../F3.1-hatta-regimes/), which owns Hatta's paper."""))

cells.append(code("""N_gasfilm_only = KP_R3 * 225.0
REACTION_N = 200
rows = []
for ka in [0.0, 0.1, 1.0, 10.0, 1e3, 1e6]:
    s = TwoFilm(KP_R3, kc_1, n=REACTION_N, ka=ka).solve(225.0, 378.0)
    rows.append(dict(ka=ka, rate=s["N"], p2=s["p2"], c2=s["c2"],
                     frac_of_gasfilm_limit=s["N"]/N_gasfilm_only))
rx = pd.DataFrame(rows)
display(rx.round(9))
REACTION_LIMIT_RELDEV = float(abs(rx.frac_of_gasfilm_limit.iloc[-1] - 1.0))
print(f"gas-film-only rate k_p p_1 = {N_gasfilm_only:.5f} g/h")
print(f"at ka = 1e6 the rate reaches it to {REACTION_LIMIT_RELDEV:.3e} relative, "
      "from below,")
print(f"the residual being exactly p2/p1 = "
      f"{rx.p2.iloc[-1]:.6g}/225 = {rx.p2.iloc[-1]/225.0:.3e}.")
print(f"\\nThat number is a GRID FLOOR at n = {REACTION_N}, not a converged limit "
      "- it moves with n:")
for n in (25, 50, 100, 200, 400):
    v = TwoFilm(KP_R3, kc_1, n=n, ka=1e6).solve(225.0, 378.0)
    print(f"   n = {n:4d}   1 - N/(k_p p_1) = {abs(v['N']/N_gasfilm_only - 1):.3e}")
print("   The reaction layer at ka = 1e6 is ~1/1000 of the film, so it is")
print("   unresolved on every one of these grids and p2 stalls where the mesh")
print("   leaves it. Reported as a floor, not as convergence.")"""))

cells.append(md(r"""## Validation

Six checks, ranked by what they can detect, followed by the table that says so.

1. **Equilibrium consistency of Table 1** — four printed $(c,p)$ pairs, one
   two-parameter form, two degrees of freedom, spanning 2.88 decades. *Can fail:*
   a single mis-read digit in any of eight table cells moves it. Corroborated
   independently by Table 1's own duplication of $c_1 = 368$ g/l at the shared
   $p_1 = 41$ mmHg, which involves no fit at all.
2. **Whitman's worked example, held out** — five printed numbers reproduced from
   Table 1 alone. *Can fail:* it depends on the reconstruction, on eqs. (4)
   and (7), and on the $K_p(\text{run 3}) = k_p$ assumption.
3. **The three-run model comparison** — one degree of freedom on the two-film
   fit. *Can fail,* and it does not pass: the one testable combination is off by
   ~1 %, three times what a one-parameter model achieves on it. The page also
   states which part of the 0.38 % rms was **not** free to fail — the nesting
   argument in section 3 — rather than presenting the whole of it as a result.
4. **The headline computed twice, by routes that share no arithmetic** — the
   deleted-liquid-film rate from the closed form $k_p(p_1-f(c_3))$, and from a
   Richardson extrapolation of the full root-solve at finite $k_c$. They agree
   to round-off. This is the only check on the page that is not a perturbation,
   and it is here for the reason given below.
5. **pymrm against the algebra** — twelve significant figures. *Structural.* A
   linear profile is exact in a conservative finite-volume scheme, so this is a
   port check and is labelled as one. Its value is below `check_agreement.py`'s
   `ABS_FLOOR = 1e-12`, so **CI does not compare it at all**; it is stated here
   rather than relied on. It does move for a wrong geometry index and for a
   missing equilibrium jump, so it is not vacuous — it just cannot see the grid.
6. **Grid independence** — also structural, and for the same reason. Three cells
   and six hundred agree to $10^{-11}$ relative. That is not convergence; it is
   the absence of discretisation error. The fast-reaction limit of section 8 is
   *also* not convergence: it is a floor set by the unresolved reaction layer,
   it moves with $n$, and it is reported as a floor.

### What the break table below cannot do

Every row in it perturbs an input and checks that a number moves. That
establishes **sensitivity**, and sensitivity is not correctness: a row proves a
metric would notice *this* defect, and says nothing about whether the undefected
value is right. A baseline that is wrong rather than insensitive passes every
row it has.

This page was built with one such baseline in it. The deleted-liquid-film rate
was originally obtained by setting $k_c = 10^{12}$ and reading eq. (4) from the
liquid side, where the answer is quantised in steps larger than the quantity
being reported — the conditioning table in section 2 prints both — and it sat
**above** `ABS_FLOOR`, so CI would have compared a machine-dependent number.
The break table below still carries the row that guards that metric, and it
moves it by more than two orders of magnitude; it passed against the wrong
baseline exactly as it passes against the right one. What caught the defect was
not a perturbation at all: the printed value exceeded its own $k_c\to\infty$
limit, which a rate increasing in $k_c$ cannot do.

So one metric on this page — the headline — is now computed a **second,
independent way** rather than only perturbed, and the sweep no longer contains a
tuning constant. Every other reported metric here still rests on perturbation
alone; that is a real limit of this table and it is stated rather than papered
over.

### The break table

Rows are produced by re-running the affected computation with one thing
deliberately wrong. Three kinds appear: **transcription defects** (the mid-dot
decimal trap that this file's own OCR falls into, and a digit transposition),
**model defects** (a wrong `nu`, the equilibrium jump removed, the equilibrium
curvature halved, an unconverged Newton solve), and **data defects** (a rate
replaced by one that a single-film model would fit).

Three rows are there because they barely move anything, and a blind spot left
unstated is an implicit claim.

- Refining the grid does not move the pymrm check — there is no discretisation
  error for it to see.
- A tenfold error in $k_c$ does not move Whitman's held-out prediction. That is
  the page's **result**, not its defect.
- Coarsening the grid from $n = 200$ to 25 moves the fast-reaction residual by
  only a few per cent — below the 10 % threshold this table calls a move. That
  is also the point: a *converging* residual would fall towards zero with $n$,
  and this one does not move much in either direction, because the reaction
  layer is unresolved on every grid tried. It is a floor, and the page calls it
  one rather than calling it convergence.

The cell ends by listing any metric no row moves; the metric block below then
prints the coverage explicitly, naming every reported metric the table does not
reach and why."""))

cells.append(code("""def eq_fit(c, p):
    A_, B_ = fit_equilibrium(c, p)
    worst = float(np.abs(10 ** (np.log10(p) - (A_ + B_*c)) - 1).max() * 100)
    return A_, B_, worst, float(np.abs(np.log10(p) - (A_ + B_*c)).max() / abs(B_))

def worked_chain(A_, B_):
    # Whitman's steps 2-6, from Table 1 alone, on a given equilibrium relation.
    pq  = lambda c: 10.0 ** (A_ + B_*np.asarray(c, float))
    dpq = lambda c: 10.0 ** (A_ + B_*np.asarray(c, float)) * np.log(10.0) * B_
    p2 = 225.0 - 41.0/KP_R3
    c2 = (np.log10(p2) - A_)/B_
    out = dict(worked_p2_run1_dev_mmHg=abs(p2 - printed("p2_run1")),
               worked_c2_run1_dev_gl=abs(c2 - printed("c2_run1")),
               fig3_axis_top_dev_log10=abs(A_ + B_*c_hi - o_hi))
    if c2 <= 378.0:                       # k_c would come out negative: no solution
        out.update(worked_kc_reldev=np.nan, run2_c2_dev_gl=np.nan,
                   run2_p2_dev_mmHg=np.nan, run2_rate_reldev=np.nan,
                   run2_gasfilm_resistance_share=np.nan)
        return out
    kc = 41.0/(c2 - 378.0)
    s = TwoFilm(KP_R3, kc, n=12, peq_=pq, dpeq_=dpq).algebraic(41.0, 204.0)
    out.update(worked_kc_reldev=abs(kc/printed("k_c_run1") - 1),
               run2_c2_dev_gl=abs(s["c2"] - printed("c2_run2")),
               run2_p2_dev_mmHg=abs(s["p2"] - printed("p2_run2")),
               run2_rate_reldev=abs(s["N"]/obs2 - 1),
               run2_gasfilm_resistance_share=(41.0 - s["p2"])/(41.0 - float(pq(204.0))))
    return out

UND = dict(
    eq_curve_worst_resid_pct=EQ_RESID_MAX_PCT, eq_curve_worst_resid_gl=EQ_RESID_MAX_GL,
    fig3_axis_top_dev_log10=AXIS_TOP_DEV, worked_p2_run1_dev_mmHg=WORKED_P2_DEV,
    worked_c2_run1_dev_gl=WORKED_C2_DEV, worked_kc_reldev=WORKED_KC_RELDEV,
    run2_c2_dev_gl=RUN2_C2_DEV, run2_p2_dev_mmHg=RUN2_P2_DEV,
    run2_rate_reldev=RUN2_RATE_RELDEV, twofilm_3run_rms_pct=RMS_TWO,
    gasfilm_3run_rms_pct=RMS_GAS, liquidfilm_3run_rms_pct=RMS_LIQ,
    dof_ratio_overpredict_pct=DOF_OVERPRED_PCT,
    gasfilm_dof_ratio_dev_pct=GASFILM_DOF_PCT,
    gasfilm_runs23_rms_pct=GASFILM_RUNS23_RMS,
    twofilm_runs23_rms_pct=TWOFILM_RUNS23_RMS,
    run2_rate_reldev_no_liquid_film=NO_LIQUID_RELDEV,
    run2_gasfilm_resistance_share=share_gas_2, run1_gasfilm_resistance_share=share_gas_1,
    kc_band_width_all_runs=KC_BAND_WIDTH_ALL, constH_secant_worst_err_pct=CONSTH_SEC_MAX,
    constH_run1_tangent_err_pct=CONSTH_R1_TAN, constH_run1_secant_err_pct=CONSTH_R1_SEC,
    Kp_range_factor=KP_RANGE_FACTOR, Kc_range_factor=KC_RANGE_FACTOR,
    table1_coef_columns_max_pct=TABLE1_COEF_MAX,
    prose_twelvefold_recomputed=RATIO_12FOLD, pymrm_vs_algebra_reldev=PYMRM_VS_ALG,
    newton_worst_rel_update=NEWTON_UPD_MAX, reaction_limit_reldev=REACTION_LIMIT_RELDEV)

breaks = []
def row(label, metric, defected):
    breaks.append((label, metric, UND[metric], float(defected)))

# ---- 1. transcription defects: the pre-1980 mid-dot trap, and a transposition
c_bad, p_bad = C_EQ.copy(), P_EQ.copy()
p_bad[np.argmin(np.abs(c_bad - 378.0))] = 5.5                  # 55 set as 5,5 -> read 5.5
A_b, B_b, w_pct, w_gl = eq_fit(c_bad, p_bad)
row("p3(run 1) 55 -> 5.5 (lost mid-dot)", "eq_curve_worst_resid_pct", w_pct)
row("p3(run 1) 55 -> 5.5 (lost mid-dot)", "eq_curve_worst_resid_gl", w_gl)
for k, v in worked_chain(A_b, B_b).items():
    row("p3(run 1) 55 -> 5.5, curve refitted", k, v)
c_b2 = C_EQ.copy(); c_b2[np.argmin(np.abs(c_b2 - 368.0))] = 386.0
row("c1(runs 2,3) 368 -> 386 (transposed)", "eq_curve_worst_resid_pct",
    eq_fit(c_b2, P_EQ)[2])

# ---- 2. the reconstruction's own leverage: refit it without its extreme point
keep = C_EQ < 400.0
A_3, B_3, _, _ = eq_fit(C_EQ[keep], P_EQ[keep])
for k, v in worked_chain(A_3, B_3).items():
    row("curve refitted without the (425, 225) pair", k, v)

# ---- 3. the printed rounding of k_p
p2_unrounded = 225.0 - 41.0/(24.1/41.0)          # K_p(run 3) unrounded, not 0.59
row("k_p taken unrounded (0.5878, not 0.59)", "worked_p2_run1_dev_mmHg",
    abs(p2_unrounded - printed("p2_run1")))

# ---- 4. pymrm model defects
for lab, kw in [("construct_div nu = 1 instead of 0", dict(nu=1)),
                ("interface equilibrium removed (p2 = c2)", dict(jump=False))]:
    t = TwoFilm(KP_R3, kc_1, n=40, **kw)
    row(lab, "pymrm_vs_algebra_reldev",
        abs(t.solve(225.0, 378.0)["N"]/t.algebraic(225.0, 378.0)["N"] - 1))
t_nj = TwoFilm(KP_R3, kc_1, n=200, ka=1e6, jump=False)
row("interface equilibrium removed (p2 = c2)", "reaction_limit_reldev",
    abs(t_nj.solve(225.0, 378.0)["N"]/N_gasfilm_only - 1))
_bad = TwoFilm(KP_R3, kc_1, n=40)
_y0 = np.r_[np.full(40, 225.0), np.full(40, 378.0), 100.0, 400.0]
_s1 = newton(lambda y: (_bad._residual(y, 225.0, 378.0), _bad._jacobian(y)),
             _y0.reshape(-1, 1), tol=1e-8, maxfev=1)
_dy = spla.spsolve(_bad._jacobian(_s1.x), _bad._residual(_s1.x, 225.0, 378.0))
row("newton maxfev = 1 (unconverged)", "newton_worst_rel_update",
    np.abs(_dy).max()/np.abs(np.asarray(_s1.x)).max())

# ---- 5. the equilibrium curvature, which is what kills additivity
A_f, B_f = A_EQ + 0.5*B_EQ*np.mean(C_EQ), 0.5*B_EQ     # half the slope, same centroid
pq_f = lambda c: 10.0 ** (A_f + B_f*np.asarray(c, float))
sec_f, tan_f, Kp_f_list = [], [], []
for r in RUNS:
    p1_, p3_, c1_, c3_ = (r["p1_mmHg"], r["p3_mmHg"], r["c1_g_per_l"], r["c3_g_per_l"])
    ex = two_film_algebraic(p1_, c3_, KP_R3, kc_1, pq_f)
    Kx = ex["N"]/(p1_ - float(pq_f(c3_)))
    m_s = (float(pq_f(c1_)) - float(pq_f(c3_)))/(c1_ - c3_)
    m_t = float(pq_f(c3_))*np.log(10.0)*B_f
    sec_f.append(abs(1.0/(1.0/KP_R3 + 1.0/(m_s*kc_1))/Kx - 1)*100)
    tan_f.append(abs(1.0/(1.0/KP_R3 + 1.0/(m_t*kc_1))/Kx - 1)*100)
for c3_ in np.linspace(5.0, 415.0, 40):
    if float(pq_f(c3_)) < 225.0 - 1e-9:
        e = two_film_algebraic(225.0, c3_, KP_R3, kc_1, pq_f)
        Kp_f_list.append(e["N"]/(225.0 - float(pq_f(c3_))))
row("equilibrium slope halved (flatter curve)", "constH_secant_worst_err_pct", max(sec_f))
row("equilibrium slope halved (flatter curve)", "constH_run1_secant_err_pct", sec_f[0])
row("equilibrium slope halved (flatter curve)", "constH_run1_tangent_err_pct", tan_f[0])
row("equilibrium slope halved (flatter curve)", "Kp_range_factor",
    max(Kp_f_list)/min(Kp_f_list))
ceq_f = lambda p: (np.log10(np.asarray(p, float)) - A_f)/B_f
Kc_f_list = []
for p1_ in (225.0, 120.0, 41.0):
    for c3_ in np.linspace(5.0, 415.0, 40):
        if float(pq_f(c3_)) < p1_ - 1e-9:
            e = two_film_algebraic(p1_, c3_, KP_R3, kc_1, pq_f)
            Kc_f_list.append(e["N"]/(float(ceq_f(p1_)) - c3_))
row("equilibrium slope halved (flatter curve)", "Kc_range_factor",
    max(Kc_f_list)/min(Kc_f_list))

# ---- 6. defects in the data that the model comparison must see
runs_gasfit = [dict(r) for r in RUNS]
runs_gasfit[0]["dW_dtheta_g_per_h"] = par_2[0]*170.0     # run 1 made gas-film-consistent
row("run 1 rate replaced by k_p*(p1-p3)", "gasfilm_3run_rms_pct",
    100*fit(lambda th: rel_resid_gas(th, runs=runs_gasfit), [np.log(0.4)])[2])
row("run 1 rate replaced by k_p*(p1-p3)", "twofilm_3run_rms_pct",
    100*fit(lambda th: rel_resid_two(th, runs=runs_gasfit), np.log([0.59, 1.2]))[2])
runs_r3 = [dict(r) for r in RUNS]
runs_r3[2]["dW_dtheta_g_per_h"] = 24.0                   # run 3's 24.1 read as 24.0
row("run 3 rate 24.1 -> 24.0 (last digit lost)", "gasfilm_runs23_rms_pct",
    100*fit(lambda th: rel_resid_gas(th, runs=runs_r3, subset=(2, 3)), [np.log(0.59)])[2])
row("p3(run 2) 0.3 -> 3.0 (lost mid-dot)", "gasfilm_dof_ratio_dev_pct",
    100*((41.0 - p3_3)/(41.0 - 3.0)/ratio_obs - 1))

runs_liqfit = [dict(r) for r in RUNS]
for rr in runs_liqfit:
    rr["dW_dtheta_g_per_h"] = par_l[0]*(rr["c1_g_per_l"] - rr["c3_g_per_l"])
row("all rates replaced by k_c*(c1-c3)", "liquidfilm_3run_rms_pct",
    100*fit(lambda th: rel_resid_liq(th, runs=runs_liqfit), [np.log(0.3)])[2])

# ---- 7. the liquid film deleted: which metric notices and which does not
row("k_c -> infinity (liquid film deleted)", "run2_rate_reldev", NO_LIQUID_RELDEV)
row("k_c -> infinity (liquid film deleted)", "twofilm_3run_rms_pct",
    100*fit(lambda th: rel_resid_nolf(th), [np.log(0.59)])[2])
row("k_c -> infinity (liquid film deleted)", "twofilm_runs23_rms_pct",
    100*fit(lambda th: rel_resid_nolf(th, subset=(2, 3)), [np.log(0.59)])[2])
row("k_c -> infinity (liquid film deleted)", "dof_ratio_overpredict_pct",
    100*((no_liquid_film_rate(41.0, 9.0, kp_f)["N"] /
          no_liquid_film_rate(41.0, 204.0, kp_f)["N"])/ratio_obs - 1))
row("k_c -> infinity (liquid film deleted)", "run2_gasfilm_resistance_share", 1.0)
s18 = run2_rate(kc_1/8); s18a = TwoFilm(KP_R3, kc_1/8, n=6).algebraic(225.0, 378.0)
row("k_c / 8", "run2_gasfilm_resistance_share", (41.0 - s18["p2"])/(41.0 - 0.3))
row("k_c / 8", "run1_gasfilm_resistance_share", (225.0 - s18a["p2"])/(225.0 - 55.0))
row("k_c / 8", "run2_rate_reldev_no_liquid_film", abs(s18["N"]/obs2 - 1))

# ---- 8. the back pressure removed altogether: the one degree of freedom
pq_0 = lambda c: np.zeros_like(np.asarray(c, float)) + 1e-30
tf0 = TwoFilm(KP_R3, kc_1, n=6, peq_=pq_0, dpeq_=lambda c: np.zeros_like(np.asarray(c, float)))
ratio0 = tf0.algebraic(41.0, 9.0)["N"]/tf0.algebraic(41.0, 204.0)["N"]
row("equilibrium back pressure set to zero", "dof_ratio_overpredict_pct",
    100*(ratio0/ratio_obs - 1))
_r0 = np.asarray(fit(lambda th: rel_resid_two(th, pq=pq_0), np.log([0.59, 1.2]))[1])
row("equilibrium back pressure set to zero", "twofilm_runs23_rms_pct",
    100*float(np.sqrt((_r0[[1, 2]]**2).mean())))

# ---- 9. transcription of the coefficient columns themselves
t1_bad = t1.copy(); t1_bad.loc[t1_bad.run == 3, "K_c"] = 0.67   # 0.067 misread
row("K_c(run 3) 0.067 -> 0.67", "table1_coef_columns_max_pct",
    100*abs(t1_bad.loc[t1_bad.run == 3, "dW_dtheta_g_per_h"].iloc[0] /
            t1_bad.loc[t1_bad.run == 3, "delta_c_g_per_l"].iloc[0] /
            t1_bad.loc[t1_bad.run == 3, "K_c"].iloc[0] - 1))
row("K_c(run 3) 0.067 -> 0.67", "prose_twelvefold_recomputed",
    t1.loc[t1.run == 1, "K_c"].iloc[0] / 0.67)

# ---- 10. k_c band: removing run 1 is the defect, and it is the result
row("run 1 removed from the fit", "kc_band_width_all_runs", KC_BAND_WIDTH_NO1)

# ---- 11. rows that DO NOT move, kept because a blind spot is a claim
row("grid n = 40 -> 3", "pymrm_vs_algebra_reldev",
    abs(TwoFilm(KP_R3, kc_1, n=3).solve(225.0, 378.0)["N"] /
        TwoFilm(KP_R3, kc_1, n=3).algebraic(225.0, 378.0)["N"] - 1))
row("k_c x 10", "run2_rate_reldev", abs(run2_rate(10*kc_1)["N"]/obs2 - 1))

# ---- 12. and one row that moves for a reason worth stating: the reaction limit
# is a GRID floor, so refining the grid is the defect that exposes it.
row(f"grid n = {REACTION_N} -> 25 (reaction layer unresolved)", "reaction_limit_reldev",
    abs(TwoFilm(KP_R3, kc_1, n=25, ka=1e6).solve(225.0, 378.0)["N"]/N_gasfilm_only - 1))

bt = pd.DataFrame(breaks, columns=["injected defect", "metric", "undefected", "defected"])
bt["moves_by"] = [("nan" if not np.isfinite(d) else
                   (f"x {d/u:.3g}" if u > 0 else f"0 -> {d:.3g}"))
                  for u, d in zip(bt.undefected, bt.defected)]
pd.set_option("display.max_rows", 100, "display.width", 200)
display(bt)

# "moves" = the metric changes by more than 10 % of itself, changes sign, or
# becomes undefined. Anything less is a blind row.
moved = {m for m, u, d in zip(bt.metric, bt.undefected, bt.defected)
         if (not np.isfinite(d)) or u == 0 or abs(d - u) > 0.1 * abs(u)}
print(f"\\nrows: {len(bt)}   distinct metrics in the break table: {len(set(bt.metric))}, "
      f"of which moved by at least one row: {len(moved)}")
unmoved = sorted(set(UND) - moved)
print("metrics in the table with NO moving row:", unmoved if unmoved else "none")
BREAK_METRICS = set(bt.metric)"""))

cells.append(code("""metrics = dict(
    # --- checks that can fail -------------------------------------------------
    eq_curve_worst_resid_pct   = EQ_RESID_MAX_PCT,
    eq_curve_worst_resid_gl    = EQ_RESID_MAX_GL,
    fig3_axis_top_dev_log10    = AXIS_TOP_DEV,
    worked_p2_run1_dev_mmHg    = WORKED_P2_DEV,
    worked_c2_run1_dev_gl      = WORKED_C2_DEV,
    worked_kc_reldev           = WORKED_KC_RELDEV,
    run2_c2_dev_gl             = RUN2_C2_DEV,
    run2_p2_dev_mmHg           = RUN2_P2_DEV,
    run2_rate_reldev           = RUN2_RATE_RELDEV,
    twofilm_3run_rms_pct       = RMS_TWO,
    gasfilm_3run_rms_pct       = RMS_GAS,
    liquidfilm_3run_rms_pct    = RMS_LIQ,
    dof_ratio_overpredict_pct  = DOF_OVERPRED_PCT,
    gasfilm_dof_ratio_dev_pct  = GASFILM_DOF_PCT,
    gasfilm_runs23_rms_pct     = GASFILM_RUNS23_RMS,
    twofilm_runs23_rms_pct     = TWOFILM_RUNS23_RMS,
    # --- what the paper's own check is worth ---------------------------------
    run2_rate_reldev_no_liquid_film = NO_LIQUID_RELDEV,
    run2_gasfilm_resistance_share   = share_gas_2,
    run1_gasfilm_resistance_share   = share_gas_1,
    kc_band_width_all_runs          = KC_BAND_WIDTH_ALL,
    kc_band_width_without_run1      = KC_BAND_WIDTH_NO1,
    # --- consequences of the non-linear equilibrium --------------------------
    constH_run1_tangent_err_pct  = CONSTH_R1_TAN,
    constH_run1_secant_err_pct   = CONSTH_R1_SEC,
    constH_secant_worst_err_pct  = CONSTH_SEC_MAX,
    Kp_range_factor              = KP_RANGE_FACTOR,
    Kc_range_factor              = KC_RANGE_FACTOR,
    # --- structural / port checks (labelled, not evidence) -------------------
    table1_diff_columns_max_dev  = TABLE1_DIFF_MAX,
    table1_coef_columns_max_pct  = TABLE1_COEF_MAX,
    prose_twelvefold_recomputed  = RATIO_12FOLD,
    pymrm_vs_algebra_reldev      = PYMRM_VS_ALG,     # < ABS_FLOOR: not CI-compared
    newton_worst_rel_update      = NEWTON_UPD_MAX,   # < ABS_FLOOR: not CI-compared
    reaction_limit_reldev        = REACTION_LIMIT_RELDEV,
)
gu.report_agreement("A3.1", metrics)

FLOOR = 1e-12
below = sorted(k for k, v in metrics.items() if abs(v) < FLOOR)
print(f"\\nNOTE - metrics below check_agreement.py's ABS_FLOOR = {FLOOR:g}, which CI")
print("therefore does NOT compare. They are structural identities, reported not relied on:")
for k in below:
    print(f"   {k} = {metrics[k]:.4g}")
print("\\nkc_band_width_without_run1 is CENSORED at the top of the k_c sweep "
      f"({KC_SWEEP_TOP:g}):")
print("the true band has no upper end, which is the point of section 5.")
print("\\nTHREE quantities printed on this page are deliberately NOT reported as")
print("metrics. Each is machine-dependent or an identity, and CI comparing it would")
print("manufacture regressions that are not regressions:")
print(f"   grid_spread_n3_to_n600      = {GRID_SPREAD:.3g} g/h - accumulated round-off in")
print("       a solve that has no discretisation error.")
print(f"   no_liquid_film_two_routes   = {NO_LF_TWO_ROUTES:.3g} - the closed form against")
print("       the Richardson extrapolation. It is round-off, and near ABS_FLOOR: the")
print("       agreement is the evidence, its last digits are not.")
print(f"   fig3_axis_bottom_dev_log10  = {AXIS_BOT_DEV:.4f} - an identity, being the top")
print(f"       deviation plus the fixed geometric shortfall {BOX_H-LINE_SPAN:.4f}, so it")
print("       carries no information the top deviation and the slope do not.")

uncovered = sorted(set(metrics) - BREAK_METRICS)
print(f"\\nBREAK-TABLE COVERAGE: {len(metrics)} metrics reported, "
      f"{len(metrics)-len(uncovered)} have a row that moves them.")
print("The remainder, and why no row exists for them:")
for k in uncovered:
    why = {"table1_diff_columns_max_dev":
           "exactly 0 by construction - the printed difference columns are the "
           "subtraction of the printed p and c columns, so there is no scale to move",
           "kc_band_width_without_run1":
           "it IS the defected value of the 'run 1 removed from the fit' row above, "
           "so it cannot also be that row's baseline"}.get(k, "UNEXPLAINED - fix this")
    print(f"   {k}: {why}")"""))

# --------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Not the two-film answer itself.** For Whitman's problem the pymrm solve
reproduces the algebra to twelve digits *by construction*, and the page says so
rather than presenting it as validation. A resistance-in-series identity does
not need a solver.

What it adds is three things the 1923 paper could not do.

1. **The interface state, solved instead of read off a chart.** Steps 3 and 5 of
   the worked example are Whitman finding a point on Fig. 3 by eye. Closing
   eq. (4) with $p_2 = f(c_2)$ makes that a one-line Newton solve, which is why
   this page can put his three hand-worked points inside a continuous sweep of
   the operating plane. That sweep is the result of section 7: at fixed $k_p$ and
   $k_c$ the "overall coefficients" move by a factor 3.80 ($K_p$) and 13.59
   ($K_c$) — both printed by the cell that sweeps them. They are properties of
   the *operating point*, not of the apparatus — Whitman's thesis, which he
   demonstrates at three points and could not map.

2. **A number on how badly the linearised form fails.** The formula everyone
   uses, $1/K_p = 1/k_p + 1/(m k_c)$, is not in Whitman's paper; what is in the
   paper is the warning that it must not be used here. Section 6 turns that
   warning into a measurement against the exact non-linear interface.

3. **The films as a differential object, so an assumption can be removed.** The
   last cell puts a reaction in the liquid film and recovers the gas-film-only
   limit that the algebra cannot express. That is a mechanism demonstration, not
   a result about HCl.

**And one thing that is not pymrm at all.** The most useful output of this page
is the power analysis of section 2 and section 5 — that Whitman's single
held-out prediction is 99.4 % gas-film controlled, and that removing run 1
leaves $k_c$ unbounded above. That needed a careful reading of Table 1 and no
solver whatsoever. It is included because it is the honest answer to the
question the 1962 Foreword raises, and because a page about a
resistance-in-series identity that reported only identities would be worthless."""))

# -------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**The pattern to lift** is `TwoFilm`: two domains, each with its own
`construct_grad`/`construct_div` pair, coupled by (i) flux continuity and
(ii) a **non-linear jump condition** between the two field variables at the
shared face, with both interface values carried as external unknowns through
`shapes_d` so the operators stay constant. Any partition-equilibrium interface
has this shape — gas/liquid, membrane/feed, solid/fluid at an adsorption
isotherm.

**Numbers you may take from this page, and the ones you may not.**

- $k_p = 0.59$ g h⁻¹ mmHg⁻¹ and $k_c \approx 1.2$ g h⁻¹ (g/l)⁻¹ are
  **per-apparatus conductances for one 1923 rig**, not film coefficients. They
  carry no area and no thickness. They are not transferable to anything.
- The reconstructed equilibrium relation
  $\log_{10}(p/\mathrm{mmHg}) = A + B\,c$ is fitted to *Whitman's own four
  printed points* at 30 °C and is a reconstruction of the curve he drew — not a
  thermodynamic correlation for HCl. Use it inside this page's argument and
  nowhere else; modern HCl vapour-pressure data are better and are not here.
- $k_c$ should be treated as **fitted, not measured**. Section 5 shows it rests
  entirely on run 1: with all three runs the admitted band is narrow, but remove
  run 1 and the remaining data admit every value above about 0.4. A tight band
  around a value that one observation sets by itself is identification, not
  corroboration. Do not quote 1.2 as a measured liquid-film coefficient.

**If you are choosing between film, penetration and surface-renewal theory,**
this page cannot help: it establishes only that a *single* film fails on these
data and that two films fit them with one degree of freedom to spare. The
diffusivity exponent that distinguishes the three pictures is not something
Whitman's data can see. **That comparison is still open, and this page does not
send you anywhere for it.** It belongs to a page holding all three sources at
once; `A3.2` (Higbie) and `A3.3` (Danckwerts) are the other two, and `A3.3`
declines it explicitly for the same reason, so do not go there expecting to find
it.

**If you need enhancement by reaction**, go to [`F3.1`](../F3.1-hatta-regimes/);
the reaction cell here is a limit check, not a model.

**Related:** [`F3.1`](../F3.1-hatta-regimes/) (reaction in the liquid film),
[`F3.5`](../F3.5-co2-amine-absorption/) (a modern two-film absorption model with
speciation), [`A2.1`](../A2.1-danckwerts-boundary-conditions/) (outward-normal
boundary conditions), [`A4.4`](../A4.4-knudsen-bosanquet/) (resistances in series
in a different guise, and where that addition is exact)."""))

# ----------------------------------------------------------------- references
cells.append(md(r"""## References

Whitman, W. G. (1923). The two-film theory of gas absorption. *Chemical and
Metallurgical Engineering* **29**, 146–148. — **the origin of the result,
cited and not consulted.** The 1923 printing is pre-DOI and unreachable; the
volume and page range are the ones printed in the 1962 reprint's own header and
are inherited, not verified. **No issue number is given**, because the reprint's
header does not print one and the 1923 printing was not seen.

Whitman, W. G. (1962). The two-film theory of gas absorption. In *Pioneer papers
in convective mass transfer, 5*. *International Journal of Heat and Mass
Transfer* **5**(5), 429–433.
[doi:10.1016/0017-9310(62)90032-7](https://doi.org/10.1016/0017-9310(62)90032-7)
— **the text actually read.** A verbatim reprint of the 1923 article, with an
Editor's Foreword. Eqs. (1)–(8), Table 1, the worked example and the printed axis
labels of Figs 2 and 3 were all transcribed from renders of this printing at its
native 300 ppi, each numeric cropped and re-read at that resolution, and every
page reference on this page is to it.

Sherwood, T. K. and Pigford, R. L. (1951). *Absorption and Extraction*.
McGraw-Hill. — quoted **through** the 1962 Editor's Foreword, which gives it as
*"their book Absorption and Extraction (McGraw-Hill, 1951)"* and reproduces the
sentence used here. No edition is stated there and none is asserted here; the
book itself is not on disk.

Whitman, W. G. and Keats, J. L. (1922). *J. Ind. Eng. Chem.* **14**, 185 —
Whitman's own earlier humidification/dehumidification work, cited by him as the
source of the two-film picture. Not consulted.

Dolezalek, F. (1898). *Z. phys. Chem.* **26**, 334, and Bates, S. J. and
Kirschman, H. D. (1919). *J. Am. Chem. Soc.* **41**, 1897 — the HCl
vapour-pressure data plotted as crosses and circles on Whitman's Fig. 3. Not
consulted, and not digitised: the equilibrium relation used here is reconstructed
from Whitman's own Table 1."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py", "mimetype": "text/x-python", "name": "python",
        "nbconvert_exporter": "python", "pygments_lexer": "ipython3",
        "version": "3.13.5"},
}
out = Path(__file__).with_name("index.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
