# A3.1 — Whitman's two-film theory, and what his three runs can actually test

Two films in series, equilibrium at the interface between them, and the
absorption rate is whatever both films can pass at once. It is one of the
most-cited results in chemical engineering, and it is **algebraic** — which
makes it dangerously easy to "validate" by rearranging it.

**Source, and the reprint route.** The 1923 original in *Chemical and
Metallurgical Engineering* **29**, 146–148 is pre-DOI, unreachable, and was
**not consulted**. No issue number is given anywhere here: the reprint's header
prints volume, pages and year only. It is reprinted **verbatim** as item 5 of *"Pioneer papers in
convective mass transfer"* in *Int. J. Heat Mass Transfer* **5**(5) 429–433
(1962), [doi:10.1016/0017-9310(62)90032-7](https://doi.org/10.1016/0017-9310(62)90032-7),
which is on disk and states its own origin in its header. Every equation
(eqs. 1–8), every cell of Table 1, the whole worked example and the printed axis
tick labels of Figs 2 and 3 were read off renders of the 1962 printing at
**that scan's native 300 ppi** on 2026-08-05 — `pdfimages -list` reports
CCITT-G4 bilevel images at 300 ppi, so a larger render only interpolates — with
**every numeric cropped and re-read at that resolution** rather than read at
page scale. **All page references on the page are to the 1962
reprint (429–433).** The 1923 page range appears only in the reprint's header;
it is inherited, not verified, and no individual 1923 page is cited. The PDF's
text layer is an Acrobat Capture OCR of the scan and is not used for any
number — it returns eq. (7) as `2,$~=PL~~2Lz` and run 3's `K_c` as `om7`, and
the 1962 typesetting's mid-dot decimal separator is discarded wholesale.

**The Editor's Foreword is not Whitman.** The reprint opens with a 1962 foreword
signed "D.B.S." quoting Sherwood and Pigford's 1951 preface — that after 28 years
the theory "has never been adequately checked experimentally". That is editorial
commentary *about* Whitman and the page never attributes it to him. It is quoted
because it is the question the page answers quantitatively.

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page: eqs. (1)–(8) transcribed; the equilibrium relation
  reconstructed from Table 1's own printed pairs; Whitman's worked example
  reproduced with his five printed numbers held out; a two-domain steady
  diffusion solve in pymrm with the equilibrium jump as a coupling condition; a
  three-model comparison on the three runs; and the power analysis.
- `data/whitman-1923-table1.csv` — the ten printed columns of Table 1 (tier 3).
- `data/whitman-1923-printed-results.csv` — the worked example, three stated
  ratios, and Fig. 3's printed axis tick labels (tier 6).

No other page's dataset is used. No figure is digitised and no page image is
reproduced.

## The equilibrium curve is reconstructed, not digitised

The obvious move is to digitise Fig. 3, and it is unnecessary. Whitman's
subscript convention makes `(c1, p1)` and `(c3, p3)` **equilibrium pairs** —
`c1` is "the concentration of liquid which would be in equilibrium with the gas"
at `p1`, and `p3` is "the partial pressure of solute exerted by the liquid" at
`c3`. Table 1 therefore prints **four distinct equilibrium points** spanning 2.88
decades of pressure, and Fig. 3 plots a straight dashed line on `log10 p` against
`c`. Fitting that two-parameter form to the four printed points leaves **two
degrees of freedom**, so the reconstruction has residuals and can fail: the worst
is **0.383 % in p (0.128 g/l in c)**. Corrupting `p3(run 1)` from 55 to 5.5 — the
mid-dot decimal trap this PDF's own OCR falls into — moves that by a factor 322.

**Table 1 corroborates that reading with no fit at all.** Runs 2 and 3 share
`p1` = 41 mmHg and the table prints `c1` = **368 g/l for both**, while their bulk
concentrations (204 and 9 g/l), their `K_c` (0.146 and 0.067) and their rates
(24.0 and 24.1 g/h) all differ. A concentration in equilibrium with the gas has
no choice; a bulk or inlet concentration would have no reason to coincide.

Two further checks use printed characters only, and the page counts **one** of
them. Extrapolated to Fig. 3's printed abscissa ticks the reconstruction gives
+2.808 at c = 460 against the printed top tick +2.8 — 0.0081 log units, and the
drawn line does exit at the top-right corner. At c = 180 it gives −0.835 against
the printed −1.0, which is *not* a failure and *not* a second confirmation: a
line of this slope spans 3.643 log units across the 280 g/l abscissa while the
printed box is 3.800 tall, so it cannot touch both corners, and the 0.165
measured at the bottom is exactly that 0.157 shortfall plus the 0.0081 at the
top. The genuine second corroboration is his prose — the back pressure at
"negligible up to approximately 250 g/l" comes out 1.19 mmHg.

## The result

**Whitman's worked example reproduced from Table 1 alone**, with his five printed
numbers held out: `p2` = 155.51 mmHg against his 156, `c2` = 412.63 g/l against
the 412 he read off Fig. 3, `k_c` = 1.184 against his 1.2, the run-2 interface
pair (0.549 mmHg, 224.16 g/l) against his (0.6, 224), and the held-out run-2 rate
**23.87 g/h against the observed 24.0, −0.558 %**.

**The three runs do discriminate.** Fitted to all three rates: a single gas film
(one parameter) misses by 41.7 % rms, a single liquid film by 59.5 %, two films
(two parameters) by **0.379 %** against Table 1's own printed precision of
±0.21 %. The page states which half of that was **not** free to fail: the
two-film family nests the gas-film model, and its second parameter absorbs run 1
to 0.005 %, so the rms is set by runs 2 and 3 alone — whose mutual consistency
Table 1 already prints as `K_p` = 0.59 twice.

**On both quantities the data can actually test, the liquid film makes things
worse.** On runs 2 and 3, a *one-parameter* gas film fits to **0.159 % rms**
where the two-parameter two-film model manages 0.465 %. And on the one remaining
degree of freedom — the run3/run2 rate ratio — Table 1 gives 1.00417 (±0.29 %)
against 1.00000 for no back pressure, 1.00737 for a gas film alone (**+0.319 %**)
and 1.01354 for two films (**+0.934 %**). Three of those four predictions contain
no fitted parameter at all, `k_p` cancelling out of a ratio at fixed `p1`.

**And the headline.** Whitman's only held-out prediction is **99.4 % gas-film
controlled**. Deleting the liquid film entirely moves the predicted rate from
−0.558 % to **+0.054 %** of the observed value — *closer*, not further. That
limit is taken **analytically**, as `k_p (p1 − f(c3))`, so no large-`k_c`
constant is load-bearing; it is confirmed a second and independent way by
Richardson extrapolation of the finite-`k_c` root solve, the two agreeing to
1.3e-12. Profiling the fit shows every subset containing run 1 admits `k_c` only
in [1.14, 1.24], while runs 2 and 3 alone admit [0.40, 3000] and are still rising
at the top of the sweep. So the liquid-film coefficient is tightly **identified**
by one run and **corroborated** by none. That is Sherwood and Pigford's 1951
remark, with a number on it.

## What pymrm adds, honestly

Not the two-film answer: the pymrm two-domain solve reproduces the algebra to
5.2e-13 **by construction** (a linear profile is exact in a conservative
finite-volume scheme on any grid), and the page labels that a port check. It sits
below `check_agreement.py`'s `ABS_FLOOR = 1e-12` and is not CI-compared.

What it adds is the interface state solved instead of read off a chart — so
Whitman's three hand-worked points become a continuous sweep, over which `K_p`
moves by a factor **3.80** and `K_c` by a factor **13.59** at fixed `k_p` and
`k_c`; a number on how badly the
textbook `1/K_p = 1/k_p + 1/(m k_c)` fails on a curve this steep (**87.8 %** with
a tangent `m`, **114.9 %** with a secant, and the error tracks the curvature);
and the films as a differential object, so a reaction can be put in the liquid
film and the gas-film-only limit reached to 3.0e-6 with that limit nowhere in
the code — reported as a **grid floor**, not as convergence, because the reaction
layer is unresolved on every grid tried and the residual moves only from 3.11e-6
at n = 25 to 2.97e-6 at n = 400.

A **50-row defect-injection table** moves 30 of the 32 reported metrics; the page
names the two it does not reach and why. It keeps three rows precisely because
they barely move: the grid does not move the pymrm check, a tenfold error in
`k_c` does not move Whitman's held-out prediction (the result, not the defect),
and coarsening the grid moves the fast-reaction residual by 4 %, which is what
makes it a floor.

**And the page says what that table cannot do.** Every row perturbs an input and
watches a number move, which establishes sensitivity, not correctness. This page
carried exactly that defect: the deleted-liquid-film rate was first obtained by
setting `k_c` = 1e12 and reading eq. (4) from the *liquid* side, where `c2 − c3`
is 2.4e-11 against an `ulp(204)` of 2.8e-14 — quantised at twice the size of the
quantity reported, above `ABS_FLOOR` so CI would have compared it, and
sign-reversing one decade higher. The break row guarding it moved it by a factor
484 and passed throughout. The headline is now computed a **second, independent
way** rather than only perturbed.

## What the page cannot conclude

It cannot choose between film, penetration and surface renewal. The diffusivity
exponent that separates those pictures is not something three runs at one
temperature can see. That comparison needs all three sources at once and is
**still open** — `A3.2` (Higbie) has no source and `A3.3` (Danckwerts) declines
it explicitly, so this page does not send anyone there for it. The reaction cell is a limit check,
not an enhancement model — `F3.1` owns Hatta. And `k_p` and `k_c` are assumed
equal across all three runs, which Whitman explicitly says need not hold and
which three runs cannot check.

Runtime: about 7 s.
