# J4.4 — Luedeking & Piret, lactic acid fermentation at controlled pH

**Six fitted pairs, one printed fermentation, and what that lets you check.**

Luedeking & Piret print their relation in the Summary itself:

> *"It was found that the instantaneous rate of acid formation dP/dt, could be
> related to the instantaneous rate of bacterial growth dN/dt, and to the
> bacterial density N, throughout a fermentation at a given pH, by the
> expression* dP/dt = α dN/dt + βN *where the constants α and β are determined
> by the pH of the fermentation."*

**The paper uses Greek α and β.** The catalogue entry for this case, and
`queue_cases/J4.4.yaml`, both write the relation with *a* and *b*. The Summary,
eq. (2), eq. (3), Table III, Fig. 10 and the Nomenclature on book p. 411 all
print α and β. This page keeps the paper's symbols and says so rather than
renaming anything silently.

**Three printed tables and nothing digitised.** Table I (book p. 401, 45 rows of
bacterial density), Table II (p. 402, 27 rows of *N*, *P*, d*N*/d*t*, d*P*/d*t*,
*k* and (1/*N*)(d*P*/d*t*)) and Table III (p. 410, the six (pH, α, β) triples).
**Every decimal point is a British mid-dot**, and this scan's text layer
relocates and drops them — it renders book p. 406's `pH 5.4` as `pH 3.4` and,
one line later, as `pH 3 4`, and the U.O.D. definition `N = 0.125r` as
`N = 0.12%`. **No number here came off the text layer**; every cell was read on
a 300 ppi crop enlarged to digit scale, and Table III's β column again at 7×.
Book page = PDF page + 392, read off the printed running heads (PDF 8 → `400`,
PDF 14 → `406`, PDF 15 → `407`, PDF 17 → `409`) on crops, not off the text
layer.

## The scope decision, and the reason for it

**Tables I and II tabulate the pH 6·0 run and no other.** The running text on
book p. 401 says so — *"As an example, the rate information for the
fermentation at pH 6·0 is tabulated in Tables I and II"* — and both tables carry
the same footnote. The other five fermentations appear only as plotted curves in
Figs. 4–11 and as fitted constants in Table III.

The check behind that negative claim, stated so that it reproduces: `pdftotext
-layout` over all twenty pages, then `grep -c "Table"` returns **five** and
`grep -ci "table"` returns **six**. The five case-sensitive hits are every
mention of a table in the paper — the three captions (the text layer renders the
last two as `Table 11` and `Table 111`), the running-text sentence on book
p. 401 quoted above, and the cross-reference on book p. 409, *"can be seen from
Table I11 and from Fig. 10"*. The sixth case-insensitive hit is the word
*"unsuitable"* on book p. 398 and is not a table at all. All twenty rendered
pages were then looked at for ruled blocks; the only three are on PDF pages 9,
10 and 18.

**So one of the six rows can be refitted and five cannot.** No curve is
digitised to manufacture the other five. **Six fitted pairs are six outputs of
the test of the form, not the test:** the form is tested by the linearity of
(1/*N*)(d*P*/d*t*) against *k*, which is Fig. 11, and those points are printed
for exactly one level.

## What it finds

**The authors' straight edge is statistically indistinguishable from least
squares — and that is the headline, and it is a consistency result about their
*fitting*, not about their *model*.** Book p. 409: *"The constants α and β are
determined from the plots of Fig. 11"*, and Fig. 11 plots (1/*N*)(d*P*/d*t*)
against *k*, which **are** Table II's last two columns. The data on both sides
of the comparison are the same measurements from the same run.

| | α (mg lactic/U.O.D.) | β (mg lactic/U.O.D./h) | rms |
|---|---|---|---|
| **printed**, Table III pH 6·0 | 2·2 | 0·55 | |
| A — eq. (3) OLS, the authors' own axes | **2.266562** | **0.528684** | 0.100792 |
| B — eq. (2) OLS, i.e. eq. (3) weighted by *N*² | 2.289913 | 0.510117 | 0.323004 |
| C — the integral form, **no differentiated column** | 2.272660 | 0.510150 | 0.153769 |

Route A is +3.03 % and −3.88 % from the printed pair, and **both are inside one
standard error** — 0.508390 and 0.427714 of one. The three routes span 1.03 % in
α and 3.64 % in β.

**Quote the pair, not the constants.** α and β come out anticorrelated at
**−0.891902**, so those two errors largely cancel. Hold α at the authors' 2·2
and β refits to **0.551282**, 0.23 % from their 0·55; hold β at 0·55 and α
refits to **2.216616**, 0.76 % from their 2·2. Their pair costs **0.65 %** in
rms against the least-squares optimum (0.101450 against 0.100792). That is the
honest measure of the disagreement, and it is small.

**The two-term form beats both of its one-term degenerate cases, and both nulls
are the paper's own.** Book p. 410: *"one can state with equal validity that
during the logarithmic phase the rate of acid production is proportional to the
growth rate of the bacteria, or that during this phase the rate of acid
production is proportional to the quantity of bacteria present. **Neither
statement was found to hold true outside of the period of logarithmic
growth**"*. Each null is refitted with its own best constant:

| null | rms, rate form | × two-term | rms, integral form | × two-term |
|---|---|---|---|---|
| dP/dt = α dN/dt (β = 0) | 1.7927 mg/ml h | **5.5500** | 2.1899 mg/ml | **14.2416** |
| dP/dt = βN (α = 0) | 1.5878 mg/ml h | **4.9158** | 2.7362 mg/ml | **17.7940** |

**The *direction* of every one of those ratios is algebra, not evidence.** Each
null is an exactly nested submodel of the two-term fit on the same rows in the
same variables, so the ratio cannot be below 1. Only the magnitude is a finding,
and the page says so.

**Which is why both halves of the paper's sentence are tested, not just the one
that flatters the model.** That sentence says the two one-term forms *are both
valid* while *k* is constant. The printed *k* column has a gap at 0·45 — it
jumps from 0·40 to 0·46 — and every row above the gap lies on a plateau
0·46–0·48, a spread of 4.3478 %. That plateau **is** the paper's *k*<sub>c</sub>,
and split there:

| rows | *n* | two-term rms | × growth-only | × density-only |
|---|---|---|---|---|
| all | 22 | 0.323004 | 5.549971 | 4.915776 |
| the logarithmic plateau, *k* ≥ 0·45 | 12 | 0.102155 | **1.036042** | **1.033197** |
| outside it, *k* < 0·45 | 10 | 0.447426 | **5.011873** | **3.708246** |

**Inside the window the paper says the one-term forms hold, they cost 3.6 % and
3.3 %; outside it they cost 5.0× and 3.7×.** That is Luedeking & Piret's own
qualifier made quantitative in both directions, and it is a stronger result than
the single all-rows ratios, which average two regimes and take their size
entirely from the second. On the plateau d*N*/d*t* = *k*<sub>c</sub>*N* makes the
two regressors nearly proportional, so a ratio near 1 there is **degeneracy, not
agreement** — the page says that too, and reports the plateau-only fit
(α = 1.614396, β = 0.792957) as the evidence.

*R*² flatters every model in the integral form, because *P*(*t*) is dominated by
its monotone rise — which is why the rms ratios are what is reported.

**A third route to the constants uses no differentiated column at all.**
Integrating eq. (2) and fitting the measured *P* column gives α = 2.272660,
β = 0.510150 — within 0.27 % and 3.6 % of Route A. And with the **authors'**
printed constants, the integrated model reproduces the whole measured acid
curve to **1.70 % rms**, worst 3.84 % at *t* = 6·50 h, and the final
43·3 mg/ml to **+3.12 %** (44.649542 predicted). **It is still not an
independent dataset**: Table II's d*P*/d*t* column *is* the *P* column
graphically differentiated (book p. 400), so this is one run reduced two ways —
a coherence check, labelled as one. For scale, the paper's own dilution and
sampling corrections *"if neglected, can be as high as 15–20 per cent"*.

**The crossover the paper describes in words, given a number.** *"In the early
phases … the first term of equation 3 is the important one, while towards the
end of the fermentation the second term becomes more important."* The two terms
are equal at *k*\* = β/α = **0.25 h⁻¹** at pH 6·0, reached at
*t* ≈ **10.1 h**, which is **72.8 %** of the way through the tabulated window.
**Root-found twice on column pairs that share nothing** — on the printed *k*
column and on the separate printed d*N*/d*t* and *N* columns — and the two
routes differ by 2.450e-05. **That residual is a transcription-consistency
check, not the precision of *t*\***: it says a mis-transcribed *k* cell would
move one route and not the other, and a break row that runs the second route on
the first route's column collapses it to zero. The actual band on *t*\* is two
orders of magnitude wider — swapping the monotone-cubic interpolant for a linear
one moves it 0.16 %, and a single mis-transcribed Table III cell moves it 6.5 %
— so it is quoted to one decimal and never to the digits `brentq` returns. Over
the whole run the two mechanisms contribute **46.7 %** and
**53.3 %** of the acid: neither is a correction to the other, even though each
dominates one half.

**Across pH, β collapses and α barely moves.**

| pH | 6·0 | 5·6 | 5·4 | 5·2 | 4·8 | 4·5 |
|---|---|---|---|---|---|---|
| α | 2·2 | 2·2 | 2·2 | 2·45 | 3·0 | 3·55 |
| β | 0·55 | 0·49 | 0·32 | 0·26 | 0·14 | 0·11 |
| *k*\* = β/α, 1/h | 0.25 | 0.222727 | 0.145455 | 0.106122 | 0.046667 | 0.030986 |

β is strictly decreasing over all six levels, by a factor of **5.0**; α is
non-decreasing and **flat at 2·2 for the top three levels**, rising by only
**1.613636** overall. So *k*\* falls **8.068182**-fold: as the pH drops the
fermentation becomes almost purely growth-associated. **What that does not
establish** is stated on the page — not that the form holds at those five
levels, not a functional α(pH) or β(pH), and not the inflection the paper reads
off its own Fig. 10 (three of the six α values are printed as the same number,
so Table III alone can put an inflection nowhere).

**Table II checks itself, and the strength of that check is measured.** Book
p. 400 defines *k* and (1/*N*)(d*P*/d*t*) as quotients of the other printed
columns, so the identities are the paper's own. They hold to **2.7356 %** across
25 rows (worst at *t* = 3·50 h) and **1.841 %** across 22 — round-off from three
independently rounded two-figure columns, not error. **Of the 2019 single-digit
substitutions available in the (N, dN/dt, k) triple, 1764 — 87.3700 % — break
the identity by more than the band the true table already shows.** The 255 it
cannot see are almost all last-digit ±1 and ±2. And Table II's *N* is footnoted
*"Interpolated"*: interpolating Table I onto its grid reproduces it to
**1.2495 % rms**, worst **3.3774 %**.

## Two printed features, reported and not repaired

**The stray mark in Table III's pH 5·6 β cell.** At 150 ppi it renders as
`-0-49`. **The value is 0·49, and it is settled three ways — arithmetic first,
pixels last.**

1. **The minus reading breaks a relation the notebook asserts.** β falls
   monotonically as the pH falls across all six levels, and the notebook
   executes `assert BETA_MONOTONE`. Read the mark as a sign and the sequence
   becomes 0·55, −0·49, 0·32, … — not monotone, and the notebook refuses to
   execute. The cell prints both readings side by side, and the break row prices
   the violation at **2800 %** on `beta_monotone_margin`.
2. **The paper plots this very column on a zero-based axis.** Fig. 10, book
   p. 407, plots β against pH with the axis labelled 0, 0·2, 0·4, 0·6 and six
   open circles all clearly above zero; calibrated on the panel's own gridlines
   the pH 5·6 marker sits at β ≈ 0·50 and all six land within *about* 0·02 of
   the printed column — "about", because the gridline centres are themselves
   ambiguous by a pixel or two and the largest deviation moves with them, and
   because what the reading has to separate is +0·49 from −0·49.
   **That reading corroborates a sign; it is not a data source**
   — nothing off Fig. 10 enters a CSV, a fit or a metric, and no curve on this
   paper is digitised.
3. **And only then the pixels.** The mark is one connected component of **8
   pixels, 4 × 3**, at (x 1023–1026, y 545–547), 3 px left of a leading zero
   that is complete and undamaged. The page's one true hyphen
   (`Rouy-Photrometer`, same page) is 11 × 4 px with area 38 — four to five
   times this mark — while the five genuine mid-dots in the same column are 3–4
   px wide, area 8–11, exactly this mark. A connected-component census of the
   whole page returns **zero** components of area 4–20 with no other ink within
   12 px: this scan produces no free-floating specks of that size anywhere.

An earlier version of this page rested the reading on the pixels alone, said
Table III *"offers no arithmetic at all"*, and quoted a looser census as *"38
hits, 37 in the Wiley watermark"*. The first two were wrong about the source and
are corrected above; the census does not reproduce — it swings between 43 and
407 depending on how "round blob" is pinned — and has been withdrawn in favour
of the criterion that does. What the mark *is* — broken type, ink speck, or
platen dirt — cannot be settled from one copy and is not guessed at.

**Table III's two footnotes are the only place the paper gives α and β units.**
The Nomenclature on book p. 411 lists *"α, β  Constants"* with no units at all.

**Table III has no *internal* identity** — six rows, three independent columns,
no derived quantity — so nothing inside it constrains anything else inside it,
and each cell rests on its own crop read. It is not short of arithmetic *around*
it, though: the pH column is printed a second time in running text on book
p. 400 (*"a different pH level (6·0, 5·6, 5·4, 5·2, 4·8 and 4·5)"*, asserted by
the notebook against the CSV), β(pH) is plotted zero-based in Fig. 10, β's
monotonicity is asserted here, and the pH 6·0 row can be refitted from Tables I
and II. Table II, by contrast, checks *itself*, and the page measures how far
that goes.

## What pymrm adds

**Not much to the fit, and the page says so.** Refitting a two-parameter linear
model to 22 rows is `numpy.linalg.lstsq`. What pymrm adds is downstream:

- **Composability.** Eq. (2) written as a source term drops unchanged into a
  batch marcher (`newton` + `NumJac((1,2))`) and into a plug-flow fermenter
  (`construct_convflux_upwind` + `construct_div`, `nu=0`, outlet via
  `compute_boundary_values`) — the form the Summary's own *"rational design of
  continuous fermentation processes"* needs and the paper never writes.
- **A closed form to be checked against.** With *k*(*t*) prescribed from the
  paper's own printed *k* column — **no growth law is imported** — eq. (2)
  integrates exactly to *P* = *P*₀ + α(*N* − *N*₀) + β∫*N* d*t*. Both
  discretisations come back first order: **1.008072** for backward Euler,
  **1.004313** for donor-cell upwind, with the Richardson outlet
  **2.224e-05** from the closed form.

**The outlet-read question, re-measured rather than inherited.** With a
zero-gradient outflow condition both reads are first order and the **last-cell
read is the closer**, by **1.016655** at ncell 3200 — the same direction J4.2
measured on different physics. `compute_boundary_values` is used anyway, because
it returns the value the flux operator transports, so a balance written on it
closes: **consistency, not accuracy**.

**No new physics, no extension to unmeasured conditions, and no resolution of
the five pH levels.** The plug-flow result illustrates composability; nothing in
this paper tests it.

## What the page cannot conclude

- **That the relation is correct.** One run, and the constants were fitted on
  it.
- **That the form holds across pH.** That needs the five other runs' raw data.
- **Anything about α(pH) or β(pH) as functions.**
- **That the plug-flow extension describes a real fermenter.**

## Cross-page

`pages/J4.1-monod/` and `pages/J4.2-andrews-substrate-inhibition/` were
published immediately before this page and both carry growth kinetics. **This
page borrows no CSV, and the reason is a finding rather than an omission**: the
notebook searches J4.1's `printed-growth-laws.csv` cell by cell for
`luedeking`, `piret`, `lactic`, `product`, `alpha` and `beta` and prints the
result — **zero hits across all six terms**. That file is about the growth law
μ(*S*); Luedeking & Piret propose no growth law at all, and eq. (2) contains no
substrate concentration. There is no number in it this page could restate, so
there is nothing to reconcile. The boundary runs the other way: J4.1 and J4.2
model growth, this page models product formation *given* growth, and the two
compose.

## Break table and blind spots

**28 defect injections; the coverage map is generated from their measured
moves**, not written by hand. All 28 move a reported metric, the map holds
**235** row-metric links, and all **62** metrics are covered. The rows split
**14 / 14** between mis-transcribing one printed cell and changing a method or a
resolution, and the notebook prints that split rather than asserting it. Among
the cell defects are five of Table III — including reading the pH 5·6 mark as a
minus, and reading the pH 6·0 β cell as 0·35, a 5-for-3 of the kind this scan's
own text layer makes on book p. 406. Among the method changes are a
log-phase-only refit, the unrounded-quotient refit, both interpolants, both
quadratures, both solver resolutions, the last-cell outlet read, and both grid
orders taken by self-convergence against the finest grid instead of against the
closed form.

**No blind spot is declared: every one of the 62 metrics is moved past the 5 %
at which `check_agreement.py` compares.** The weakest cover on the page is
**5.62 %**, on `pfr_boundary_over_last_cell` — which was the one metric under
that tolerance until a review broke it. The page used to defend that metric as
*nearly structural*: it is a ratio of two reads of the same solve, they agree to
1.7 % on every grid across a 16× refinement, and so "any defect that changes the
solution moves both reads and leaves the ratio where it was". **That is false as
a universal, and the counterexample is now a row.** It holds for defects acting
*away* from the outlet — the same single-digit substitution one printed row
earlier moves the ratio 0.15 %. But k at t = 13·50 h is the last printed k cell,
so on a grid whose axial coordinate is residence time it sits in the *outlet*
cell, and a defect there changes the local gradient, which is the one quantity
the two reads differ by. Injected into the solve and into the closed-form
reference together, it moves the ratio **5.62 %**. The ratio was never
structurally protected; it was protected by where this page happened to inject.

An earlier version of this page declared six weak metrics and added that each
was *"moved by more than 5 % by a row that changes the method"*. **That
mitigation was false by construction** — weak is *defined* as the maximum move
over all rows — and it is retracted on the page in the block that replaces it.
The notebook now recomputes the blind spot with the review's rows taken back out
of the coverage map and prints the ladder: six weak metrics before either review
pass, one after the first pass added two rows, none after the second added two
more. Five of the six left on the first pass and the last on the second, and for
all six the strongest mover is one of those four rows — not one left on a row
that was already in the table. A previous version of this paragraph said *three*
of the six: wrong on the page's own coverage data, and spelled as a word, which
no sweep on this page can see. Every count in it is now checked against the live
computation.

No metric sits below `ABS_FLOOR = 1e-12`.

**82 prose and metadata values are checked against the live computation** —
every count this file states about the break table among them — and a
mechanical sweep then reads `meta.yaml`, this file, `models_entry.yaml`, the
three data sidecars **and the notebook's own markdown** for two token classes:
numbers written to five or more decimals, which must match a live value to half
an ulp of their own printed digits, and **integers of two or more digits, which
must match a count the computation produced or a pinned source-derived constant
(book page, ppi, pixel bound, the year)**. The integer half is new: the decimal
sweep could not see a single count on the page, and both of the factual errors a
review found here were counts — *"167 row-metric links"* where the notebook
printed 176, and *"29 rows"* for a CSV whose next cell printed 30. The sweep
reports its
**achieved detection rate** rather than claiming completeness: small counts and
source constants are dense, so a corrupted digit sometimes lands on another
allowed value. Both token counts are shape-dependent and both shapes are pinned.

**And it reports what it cannot see at all.** Its integer class is *digit
strings of two or more digits*, so a count spelled as a word and a count written
as a single digit are outside it entirely — the notebook counts both classes in
the swept prose and prints the totals. That hole is not hypothetical: the third
error a review found here was *"three of those six"*, word-spelled, on three
surfaces. Word numbers were not folded into the matcher, because the small
integers this page legitimately produces are dense enough that almost every such
token would match something and the sweep would then advertise teeth it does not
have. The counts are protected by the pinned-claims list instead.

## Reuse

- **Fit α and β in the variables you will report them in.** Eq. (3) and eq. (2)
  are the same model; least squares on them are two different fits, worth 1.03 %
  on α and 3.64 % on β here.
- **Report *k*\* = β/α.** It is the well-determined combination, it has physical
  meaning, and the paper's own qualitative account of the fermentation is a
  statement about it.
- **Check a printed rate column against the columns it was divided by.** That
  test caught 87.4 % of single-digit substitutions in this transcription.
- **Look for the arithmetic before you reach for the pixels.** Where a table has
  no internal identity — Table III has none — the arithmetic is usually still
  there, one step out: a column reprinted in running text, the same quantity
  plotted on a zero-based axis, a monotonicity your own notebook asserts. This
  page settled an ambiguous glyph on component geometry and *said* the source
  offered nothing better; it did, twice. A pixel argument is the last resort,
  and it should be made to corroborate rather than to carry.
- **Test both halves of the sentence you are quoting.** The paper says the
  one-term forms are valid inside the logarithmic phase and fail outside it.
  Reporting only the 5.5× and 4.9× over the whole run makes the model look
  better than the authors claimed, and hides that the two constants are nearly
  unidentifiable exactly where the one-term forms work.
- **Fitted parameters at *n* conditions and raw data at one means you can check
  one condition.** Check it properly and scope the rest out.
