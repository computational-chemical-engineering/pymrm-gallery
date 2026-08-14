# J4.1 — Monod

**One saturating growth law, four printed alternatives, and what a growth curve
can actually decide between them.**

Section 1.5.2 of Froment, De Wilde & Bischoff prints the Monod equation as
eq. (1.5.2-1), defines its two constants by exactly two properties — the
asymptote `r_m` and the half-rate point `K_S` — and couples it to the biomass
balance eq. (1.5.2-2) and the constant-yield relation eq. (1.5.2-3). Its
Chapter 1 reference list carries `Monod, J., Annu. Rev. Microbiol., 3, 371
(1949)` — the exact citation this case catalogues, confirmed here on a 300 ppi
crop. Rawlings & Ekerdt p. 596 prints the same equation beside four
alternatives — Blackman, Tessier, Moser, Contois — and calls it *"the simplest
form of the Langmuir adsorption isotherm"*.

**The 1949 paper was not consulted.** It is not on disk and could not be
obtained. Nothing on this page is attributed to Monod beyond what Froment prints
about him; there is no sentence of the form "Monod showed/assumed/found"
anywhere in this page, its `meta.yaml` or this file.

## What it finds

**Two of the four alternatives cannot be separated from Monod by *any*
constant-yield batch run, and this is exact.** Moser contains Monod at `n = 1`.
And **Contois *is* Monod, reparameterised**: constant yield makes the biomass an
affine function of the substrate, `X = φ(a − S)`, so

```
μ_m S/(K_sx X + S)  =  (μ_m/β) S / (K_sx φ a/β + S),     β = 1 − K_sx φ
```

which is Monod with `μ_m' = μ_m/β`, `K_S' = K_sx φ a/β`. The sympy residual is
exactly zero, and two optimisers started from different places on the two
different parameterisations land on the same curve — agreeing to better than
**1×10⁻⁴** in both constants, with residual RMS differing by less than 1×10⁻¹¹
h⁻¹. Those three residuals are quoted as bounds, here and everywhere else in this
page's metadata, and are kept out of `agreement.json`: the identity makes them
exactly zero, so what they measure is `least_squares`' stopping rule. The exact
values are printed in the notebook, where nothing pins them.
Of the two alternatives that *are* separable at this design, **Blackman needs
1.30 repeats of the run and Tessier 10.3** — and `N*` is the threshold at which
the *expected* excess sum of squares first reaches the critical value, so it is a
~50 %-power point: one clean run is *within 30 % of* enough for Blackman, not
enough.

**Monod and Blackman agree exactly on both of the properties Froment defines the
constants by.** Blackman reaches `μ_m`, and at `S = K_S` — which sits on its
*linear* branch — it gives exactly `μ_m/2`. The two part company at the branch
switch, where the gap is **exactly `μ_m/3`**; on the lower branch the largest gap
is `μ_m(3−2√2)/2 = 0.0858 μ_m` at `S = (√2−1)K_S`. Both closed forms are checked
against a root-found extremum.

**On the one real dataset available, the ranking of the five laws is decided by
the data reduction, not by the data.** Levenspiel's Problem 29.18 prints seven
intervals of a batch run. Fit them the way the exercise invites — `μ_i =
ln(C_C^end/C_C^start)/Δt` against `μ(C̄_A)` — and Monod is **last** of the four
two-parameter laws (RMS 0.05469 h⁻¹ against Blackman's 0.02677) and Moser's
exponent comes out **n = 1.93**, clearing the 95 % F test at **F = 12.27 against
F(1,4) = 7.71** — a verdict **co-caused by the one printed cell this page flags**:
put that cell at its balance-implied 60.18 and route A gives **F = 5.45** against
the same 7.71, so the test is not evidence about the reduction on its own. The
ranking claim does not depend on the cell (Blackman still beats Monod with it
imputed, 0.04065 against 0.04687); what moves is which law comes last of the
four. Do the interval average *exactly* instead — the same seven growth rates,
the same five models, no new parameter, with the substrate trajectory taken from
the table's own material balance rather than from the printed averages — and
Monod's residual falls
**42 %** to 0.03147, Blackman becomes the worst of the four, `n` drops to
**1.51**, and **F falls to 1.44**.

**A seeded bootstrap prices that, with a control, and every draw is refitted
the way route A fits the real data.** On data the Monod model itself generates,
route A's nominal-5 % F test rejects Monod **36.5 % ± 1.1 %** of the time and
puts the median exponent at **1.61 ± 0.0074**. Run the identical test on data drawn from route A's
*own* model — where it is correctly specified — and it sits **at its nominal
size, 4.95 % ± 0.5 % against 5 %**, returning a median exponent of **1.01 ± 0.0049**. So the
small-sample F approximation at n = 7 is essentially exact and is not the
culprit; the whole of the distortion is the reduction. The observed F = 12.27 has
**p = 0.23 ± 0.0094** against a properly generated Monod null. Those six numbers
are Monte-Carlo estimates over 2000 seeded draws — **four frequencies and two
medians**, not six frequencies — and the page prints the sampling error beside
each: binomial for the frequencies, from the order statistics for the medians. At
that draw count the four frequencies are seed-independent in their **leading
figure only** and the two medians to at least **two** figures; the four-figure values on
the page are exact for the seed and reproducible, not resolved.

The first version of this page got those six numbers wrong, and the error is now
a **break row**: it read each draw's fit off the scan grid instead of refining
it, while comparing the result against a *refined* observed F. An unrefined Moser
fit cannot reach its own optimum, so it deflates F on every draw — grid-only, the
control reads a **conservative-looking 2.1 %** and a median exponent of 0.95. A
grid-limited extremum is invisible to every perturbation test, so the page prints
the refinement's convergence instead — as the count of draws whose two-dimensional
fit still ends above the one-dimensional fit nested inside it, taken from the grid
start alone, because the shipped fit's second start sits exactly on the nesting
point and makes that count zero at every refined depth by construction (the grid-only row has no pair of starts at all). The page prints
both counts and says which of them is a measurement. The conclusion is unchanged and stronger: a
control **at** nominal size exonerates the F approximation better than a
conservative one does. One convention worth stating with it — the forward size is
a bias-to-noise ratio and reads 0.3655 ± 0.0108 at route B's residual scatter
against 0.2355 ± 0.0095 at route A's, while the control's size is scale-free
(0.0495 ± 0.0049 and 0.0470 ± 0.0047),
which is exactly what a correctly specified test should do.

**The reduction bias is larger than the signal it is judged against.** At
`C̄_A = 9` the shortcut misstates the rate by **0.0663 h⁻¹**, which is **1.21×**
the residual scatter that route A ranks the five laws on. It has a sign: `μ` is
concave, so `μ(mean) > mean of μ`, and the shortcut therefore makes every
saturating law look too gradual near exhaustion — which is exactly the direction
that inflates Moser's exponent. The last interval consumes **1.89 times** its own
printed average substrate concentration.

**The experiment that does separate Contois from Monod is named and priced.** A
chemostat, sweeping the *feed*. Monod washes out at `D_c = μ_m S_f/(K_S+S_f)`;
Contois at `D_c = μ_m` whatever the feed. The two cross at **exactly `S_f = a`**,
the substrate intercept of the batch run's own material balance (root-found
197.116 against the balance's 197.120), and separate on either side: the gap
exceeds 10 % below `S_f = 95.5` and reaches a factor **1.83** at `S_f = 20`.

## One printed cell, reported and not repaired

All three books assume a constant yield, which forces `C̄_A` to be an affine
function of the interval-average biomass — a check needing no kinetics at all.
**Six of the seven rows sit on that line at R² = 0.99996; all seven give
0.98505.** Route A reads all seven printed `C̄_A` directly; route B never reads
that column at all, and the `φ` and `a` it does use come from the **six-row**
line, the one that leaves row 4 out — so the datum route A leans on hardest is
the one route B's balance line drops, and a break row prices that by refitting
route B on the all-seven line.

A leave-one-out sweep says *which* row: dropping row 4 gives 0.99996,
dropping any other leaves R² in [0.9757, 0.9848]. The six-row line predicts
`C̄_A = 60.18` for row 4 against the printed **43** — a residual of −17.18 where
the largest of the other six is 0.65.

A candidate explanation is recorded **as an inference and then refused as a
systematic**: 43 is within 3.1 % of that interval's *end*-of-interval substrate
concentration (44.33), but applying the same substitution to rows 3 and 5 misses
their printed values by 16.8 % and 53.6 % (74.89 against a printed 90, and 13.46
against 29). **Nothing is edited in the CSV**, and break rows show exactly what
the cell costs: put it back into the balance line, or propagate its
balance-implied value through every route-A fit **and through route A's F test**,
where it takes F from 12.27 to 5.45 and flips the verdict.

## A defect that is not one

Rawlings' Tessier row prints `μ = μ_m(1 − e^{−K_s S})`, in which `K_s` must be a
*reciprocal* concentration — unlike the Monod and Blackman rows. That looks like
a misprint for the commoner `exp(−S/K_s)` until the Moser row is read, where
`K_s` adds to `S^n` and so carries units of concentration^n. The book is reusing
`K_s` as a model-local symbol, and the Moser row settles it. All three
`pdftotext` modes and the 300 ppi render agree on the printed exponent. Not
repaired, and not reported as a defect.

## Sources

**Read from (the E1.1 target)** — Froment, G. F., De Wilde, J. & Bischoff, K. B.,
*Chemical Reactor Analysis and Design*, 3rd edn, Wiley (2011), ISBN
978-0-470-56541-4, section 1.5.2, book pp. 26-29 (PDF 68-71), eqs. (1.5.2-1) to
(1.5.2-3), with the chapter reference list on book p. 58. The filename says
"Froment_Bischoff"; this is the **third** edition and De Wilde is a full author.
The clean born-digital text layer **loses every operator in an equation, and not
by dropping it**: each Symbol-font operator extracts as an unmappable
Private-Use-Area glyph — `U+F02D` for the leading minus of eq. (1.5.2-3),
`U+F03D` for the `=` beside it, `U+F02B` for the `+` of (1.5.2-1) — which renders
as nothing. All three `pdftotext` modes emit the character, so what a reader
copies is an equation with no operators at all rather than one missing a sign.
This corrects the diagnosis the just-published `J4.6` gives for eq. (1.5.1-17) of
this same book; the lesson (read the render) and the transcriptions are right in
both places.

**Origin, cited and NOT consulted** — Monod, J., *Annu. Rev. Microbiol.* **3**,
371-394 (1949). Not on disk, not read.

**Model catalogue** — Rawlings, J. B. & Ekerdt, J. G., *Chemical Reactor Analysis
and Design Fundamentals*, 2nd edn (2025 printing), book p. 596 for the five
growth laws and the Langmuir remark, p. 595 for the one sentence that cites all
five, pp. 596-597 for the chemostat balances
eq. (10.18) and steady states eq. (10.19). Its reference [12] is the same 1949
paper. **Its five growth laws are cited to Bailey & Ollis and Shuler & Kargi**,
not to Blackman, Tessier, Moser or Contois personally, and nothing here is
attributed to those four.

**Data** — Levenspiel, O., *Chemical Reaction Engineering*, 3rd edn (1999),
Problem 29.18, book p. 644, read on a 600 ppi native render. Levenspiel
attributes the table to **"Monod, 1958, p. 74"** — the second edition of
*Recherches sur la Croissance des Cultures Bacteriennes* — **not** to the 1949
paper. That book is not on disk either and was not consulted. Which document
prints this table is not established here. **The table carries no units** for
either concentration column, and its `C̄_A` header carries an **overbar** that
every `pdftotext` mode drops.

## Files

| file | what it is |
|---|---|
| `index.ipynb` | the page; nine sections, executed clean, 110 s |
| `build_page.py` | regenerates `index.ipynb` |
| `data/levenspiel-p644-monod-batch.csv` | the seven printed intervals, with a sidecar that leads on the provenance chain and the flagged row |
| `data/printed-growth-laws.csv` | 30 equations and prose claims from the three books, verbatim, flagged and not repaired |
| `agreement.json` | 74 metrics; 12 sit below CI's `ABS_FLOOR` and are named structural with above-floor companions |

## Validation summary

- **five symbolic identities**, all with sympy residual exactly zero: Monod ≡
  Langmuir; Monod's asymptote and half-rate point; Blackman sharing *both*;
  Moser ≡ Monod at n = 1; Contois ≡ Monod under a constant yield. Plus a
  *positive* statement that Tessier admits no such reduction;
- **a stated numerical result reproduced**: Rawlings' eq. (10.19) at their
  printed parameter set, to 0 (double precision) at three dilution rates, with
  their **printed formula** `D_c = μ_m S_f/(K_s+S_f)` — which evaluates to 5/6 at
  their own printed parameter set, an evaluation of theirs by this page and not a
  number they print — recovered to 2.1×10⁻⁴ from a marched model — and break rows that
  change their printed `K_s` and `y` and watch the reproduction fail;
- **a second, independent computation of the headline fit**: route B's Monod
  constants re-derived through a pymrm backward-Euler marcher
  (`NumJac((1,2))` + `newton`, order 0.9993, Richardson extrapolated) that
  carries `S` as an unknown rather than substituting the yield — agreeing with
  LSODA to **1.7×10⁻⁶** and **9.4×10⁻⁶**, with `S + X/φ` closing at round-off,
  below 1×10⁻¹⁴ (a tautology of the scheme, so it is quoted as a bound);
- **an exact identity for the spatial solve**: the pymrm plug-flow outlet equals
  the batch at `t = τ` (Levenspiel's Fig. 29.1 caption), observed order 0.9981,
  Richardson 3.2×10⁻⁶;
- **79 defect injections, all 79 of which move their metric** by more than
  1×10⁻⁶ relative, covering 60 of the 74 metrics; the other 14 are named
  individually as structural with an above-floor companion each, and both
  generic not-coverage labels are asserted never to be used. **Six** quantities
  the page *proves* are exactly zero are printed but deliberately **not
  reported**, because CI compares everything above `ABS_FLOOR` at 5 % and they
  are optimiser and root-find exit tolerances — and the rule is asserted over
  the whole enumerated class of 18 such quantities, not over the list of
  exclusions, which is how three of the six were found still being compared;
- **96 prose/metadata values plus 10 break-table lookups** checked against the
  live computation by an assertion cell, **plus a mechanical sweep** of this file,
  `meta.yaml`, `models_entry.yaml` and both data sidecars that requires every
  number written in scientific notation — in the plain `e`-notation the YAML
  files use *and* in the typeset-superscript spelling this file uses, the second
  of which it could not read at all until this pass — to match a live value to
  **half an ulp of its own printed digits**. The page fails to execute if any drifts — and the sweep
  measures its own teeth by corrupting the last printed digit of every token it
  finds and reporting how many of those corruptions it rejects, rather than
  claiming completeness;

## Reuse

Integrate the growth law across each interval; do not evaluate it at the
interval-average substrate concentration. Before comparing growth laws, ask
whether your experiment can separate them — on a constant-yield batch run,
Contois and Moser cannot be separated from Monod at all. To separate Contois,
sweep a chemostat's feed concentration. And check any printed batch table against
its own constant-yield material balance first; it is free, and it is sharp.
