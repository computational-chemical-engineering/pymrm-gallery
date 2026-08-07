# A3.5 — Ranz–Marshall / Frössling

`Nu = 2.0 + 0.60 Re^(1/2) Pr^(1/3)` and its mass-transfer twin, from Ranz, W. E.
& Marshall, W. R., Jr., **"Evaporation from Drops", Parts I and II**, *Chemical
Engineering Progress* **48**(3) 141–146 and **48**(4) 173–180 (1952). Pre-DOI, no
DOI exists; both parts read from 400 dpi renders of the scans on disk.

## The point of the page

The correlation contains four numbers with four different statuses, and the page
never mixes them:

| constant | status |
|---|---|
| `2.0` | **theoretical** — Part I eq. (7), derived from the stagnant spherical field. Fitted to nothing. |
| `0.60` | **fitted**, to Figures 6, 7 and 9 — which are Tables 1, 2, 3 and 4. |
| `1/2` on `Re` | **assumed**, from Frössling. |
| `1/3` on `Pr`/`Sc` | **assumed**, from Frössling. |

There is **no held-out set anywhere in either part**, so every rms quoted against
those tables is a goodness of fit and is labelled `IN-SAMPLE` where it is printed,
with null baselines beside it. The intercept is the exception: it was never
fitted, so extrapolating the runs to `Re = 0` is a genuine test of it — and the
paper does that extrapolation three times, in words.

**Headline.** Free two-parameter extrapolations to `Re = 0`:

| set | runs | intercept | vs the theoretical 2.0 |
|---|---|---|---|
| heat, water (Table 1) | 18 | 2.0582 | +2.91 % |
| mass, benzene (Table 4) | 13 | 2.0523 | +2.61 % |
| mass, water (Table 1) | 15 | **1.7322** | **−13.39 %** |

The two non-water sets land on the same answer to within the noise; only water
misses. Run-to-run repeatability, from the paper's own three replicate pairs, is
0.68 % mean and 1.87 % worst, so the last row is about 20 σ. The still-air row says
the same thing more bluntly: it prints `N_Nu = 2.23` and `N_Nu' = 1.79` where eqs.
(21) and (22) both say 2.0 — a ratio of 1.2458 where the correlation demands 1.

**The authors saw this and named the cause.** Part II folio 173: "Data for mass
transfer show a steeper slope and a lower intercept, but the disagreement is always
less than 10 per cent at a given `N_Re`", followed by a list of contributing factors
that opens with "(1) inaccurate values of diffusivity". The page does **not** claim
the gap went unnoticed or unexplained. What it adds is the arithmetic: their 10 %
bound is about the *curves*, where the fitted 0.60 absorbs most of the offset, not
about the intercept, where nothing was fitted; the intercept deficit is 13.4 %; and
it is species-specific — benzene, whose diffusivity came from Hirschfelder, Bird and
Spotz rather than from the authors' own Figure 5, lands where the *heat* data land.

`N_Nu' ∝ 1/D_v`, so closing the water gap needs a diffusivity smaller by 1.117;
Part II prints exactly that back-out (0.204 sq.cm/sec at 290 K and 741 mm Hg) and
calls it "a low value compared with other methods of determination", and Part I says
its water-vapour diffusivities are "approximately 10 per cent lower" than the
International Critical Tables value.

**Where the page does contradict the papers is on the conclusions.** Part I's
abstract claims the study "confirmed the analogy between heat and mass transfer at
low Reynolds numbers, and verified the simple expression for the Nusselt number at
zero Reynolds number" — the still-air analogy ratio is 1.2458. Folio 174 concludes
that the correlation's success implies their calculated `D_v` "may be more accurate
than any reported in the literature" — the intercept says it is 1.117 too large.
Both cannot stand; the page does not say which is right, only that the tables decide
against the stronger claim.

## What pymrm does here

The algebra needs no solver. Two things do.

1. **The `2.0`, computed instead of quoted.** Part I asserts eq. (7); the notebook
   solves the spherical BVP it comes from (`construct_grad`, `construct_div(nu=2)`,
   both ends Dirichlet, geometric radial grid, surface flux read from the pymrm
   face gradient), observed grid order 2.01, Richardson in `h` and then a linear
   extrapolation in `a/b`, giving 2.000000904 — **4.52e-07 from the exact answer**,
   which adaptive quadrature of the series-resistance integral supplies as
   2.000000000000. The page states plainly that this quadrature is *identically* 2
   for any working scheme, so what is reported is a finite-volume solve measured
   against the analytic limit, not two independent estimates agreeing. Load-bearing,
   not decorative — stopping the domain at ten drop radii and calling the answer 2.0
   is 11.1 % wrong, and no grid study would reveal it.
2. **Ranz's own 1993 caveat, made quantitative.** He wrote that the crude
   accounting for "the radial convection at a spherical boundary caused by
   diffusion at the same boundary" was "unfortunately not emphasized in the
   original papers". Adding `construct_convflux_upwind` to the same solve, with
   the total molar flow closed by a scalar `brentq`, settles which average `p_f`
   is (Part II's notation does not say): the **logarithmic** mean, which makes the
   paper's mass-transfer accounting *exact* for the isothermal spherical stagnant
   film; the arithmetic mean is 0.41 % out by `y = 0.20`. The correction reaches
   +1 % on `Sh` at a surface mole fraction of 1.97 % and +10 % at 17.61 %, both
   root-found.

## A convention the papers never state for the tables

The Reynolds column is over-determined by the columns printed beside it, so it
serves as the transcription check. Asking it *at what temperature the air
properties were evaluated* — leaving `T_ref = T_a + w(T_d − T_a)` with `w` free and
minimising by bounded Brent — returns **w = 0.4996**, the film temperature, with
the fitted temperature exponent landing on kinetic theory's −1.7 (−1.765 ± 0.027
against −1.323 ± 0.033 at the free stream) and Table 1's residual collapsing to
**0.126 % rms, 0.288 % worst** over eighteen four-figure Reynolds numbers.

**Two independent corroborations, both new in this revision.** Part II does state a
reference temperature once, for the still-air run: `D_v = 0.204` sq.cm/sec "at an
average temperature of 290° K. and a pressure of 741 mm. Hg." Table 1's still-air row
prints air 24.9 °C and drop 9.1 °C at 741 mm Hg — a film temperature of **290.15 K**,
matching to 0.15 K with the pressure matching exactly. And **Table 3 is a genuine
held-out set**: it prints no Reynolds column at all, so none of its rows can enter
the fit, yet the recovered law predicts its printed abscissa to **0.823 % rms**
against 1.769 % at the free stream and 7.201 % at the drop surface, at reference
temperatures 330–410 K against the 281–329 K the law was fitted on. Table 3 on its
own prefers `w = 0.514`. Scope is stated on the page: out of sample for the
*convention*, not for the correlation, since the 0.60 was fitted to Table 3 through
Figure 9.

That is a reuse instruction: evaluating `Re` at the free stream instead moves it by
up to 15.5 % on the paper's own rows and the predicted `Nu` by up to 5.4 %, and the
error grows with the air-to-drop temperature difference — worst exactly where
drying calculations live. The identification has a stated width (the band of `w`
within 5 % of the best residual is 0.178 wide, a width and not a confidence
interval), so it is "the film temperature, not the free stream and not the drop
surface", and not "exactly one half".

## Files

```
index.ipynb                                  nine sections, executed clean, 2.8 s
meta.yaml                                    page metadata
agreement.json                               62 metrics
data/ranz-marshall-1952-table1.csv           Part II Table 1, water in dry air (18 + still air)
data/ranz-marshall-1952-table2.csv           Part II Table 2, water in 66-90 C air (9)
data/ranz-marshall-1952-table3.csv           Part II Table 3, water in 85-221 C air (9)
data/ranz-marshall-1952-table4.csv           Part II Table 4, benzene (13)
data/ranz-marshall-1952-printed-constants.csv  24 constants, each with its role
```

Every CSV has a `.meta.yaml` sidecar with a `columns:` block. All five are table
transcriptions — **tier 3, no digitisation anywhere**, and no page image is
reproduced.

## Break table and coverage

84 rows against 62 metrics. Every metric has a break row, **and** every metric has
at least one row that moves it past `check_agreement.py`'s own 5 % tolerance — a
row that moves a number by less than CI would notice is not coverage, and the
notebook asserts the distinction (78 of 84 rows clear it). No metric falls below
`ABS_FLOOR = 1e-12`, so all 62 are inside the regression suite; the notebook
asserts that too. Six rows do not clear the tolerance; **exactly one** of them is a
deliberate "barely moves" exhibit and carries that label — the discarded reading of
Table 4 run 11's degraded digit. The other five are ordinary rows whose metric is
insensitive to that particular perturbation and which are covered elsewhere, and the
page says so rather than presenting all six as intentional.

The notebook also states what the table cannot do, **with a worked example from this
page's own history**: every row perturbs an input and watches a number move, which
establishes sensitivity and never correctness. Reading Table 4 run 11 as 10.6 instead
of 10.0 moved every reported metric by under 1.5 % — inside CI's 5 % tolerance — so
the break table, the coverage assertion and `check_agreement.py` all passed on a
wrong number. What catches it is a *constraint*, not a perturbation: that row's
`N_Nu'/N_Nu` ratio must match its neighbours', and it did not, 6.90 % against 1.03 %.
That ratio is now a reported metric with its own break row. The four defences used
instead of the table are the exact-answer check on the 2.0, the root-finding of every
threshold and limit, the over-determination checks on the transcriptions (the
Reynolds column, and Table 4's ratio), and the out-of-sample test of the recovered
reference temperature on Table 3.

## Reading the sources

**No text layer at all** — `pdftotext` returns one byte per page on both parts — so
every digit came off a 400 dpi render (the files' native resolution; `pdfimages
-list` reports JPEG RGB tiles at 400 × 400 ppi, so rendering larger only
interpolates). Each numeric column was cropped on its own and magnified 2–8×, and
ambiguous glyphs were compared against unambiguous instances of both candidate
digits in the same column.

**Part I's page range, recorded both ways.** Part I's scan frames carry two folio
lines — one above the black edge of the leaf, one inside it at the foot — differing
by exactly two. Four of the six frames show both legibly (139/141, 140/142, 141/143,
144/146); on the fourth the upper line is cut off at the leaf edge and on the fifth
there is none, so those two pairs are inferred from the sequence rather than read.
Part II
settles which belongs to the page: its first page has a blank top margin and prints
`Vol. 48, No. 4 | Chemical Engineering Progress | Page 173` at the **foot** of the
text block, so in this issue the folio is a footer and the upper line on a Part I
frame is the preceding leaf's footer caught in the same frame. Read that way Part I
is **141–146**, which is the universal citation and the one Ranz gives himself in
his 1993 ISI *Citation Classic*. All folios cited on the page are the ones printed
at the foot of the page the value appears on. Part II is 173–180 either way.

**Three unreadable cells and one digit settled by arithmetic**, both handled in the
open: an ink blot covers the `N_Nu'` entries of Table 1 runs 12–14, which are left
empty rather than guessed and for which no specific values are asserted. Table 4 run
11's final glyph is degraded and the glyph does **not** decide it — an earlier draft
read it 10.6 from a three-pixel counter and was wrong. It is settled instead by the
table's own arithmetic: both of Table 4's transfer numbers come from the same
measured rate, so their ratio depends only on the printed temperatures and pressure,
and runs 10 and 12 (within 0.2 K, same pressure) bracket run 11's `N_Nu'` in
[9.73, 10.06]. With the first two glyphs reading "10." and the last as wide as this
column's round digits, the value is **10.0**.

**Two files on disk are not this paper.**
`MISC-Ranz1993-citation-classic-commentary-CurrContents22.pdf` is a one-page 1993
reminiscence carrying none of the correlations — quoted once, for the caveat above
— and `Charlesworth1960-evaporation-drops-dissolved-solids-AIChEJ6-9.pdf` is a
companion study on a different problem.

## Scope

Figures 5 (transport properties), 12 (the five-decade master correlation), 13 (free
convection) and 14 (vapour pressure) carry content that exists only as curves. No
maintainer figure review is available for this case, so nothing is digitised and
nothing on the page depends on them. Two consequences are stated rather than worked
around: **the free-convection exponent 1/4 of eqs. (10) and (11) is never exercised
here**, and Table 6's calculated evaporation rates cannot be reproduced. Table 5 is
not used — it carries no transfer coefficient. Packed beds are `A3.4`.
