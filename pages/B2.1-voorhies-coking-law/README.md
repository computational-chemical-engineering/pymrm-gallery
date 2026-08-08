# B2.1 — Voorhies' coking law

**Catalog ID:** `B2.1` · **Section:** B · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1` · **Data tier:** 2 (tables printed in the paper — the
author's own laboratory measurements)

Reproduces Alexis Voorhies, Jr., "Carbon Formation in Catalytic Cracking",
*Industrial & Engineering Chemistry* **37**(4) 318–322 (1945),
doi:10.1021/ie50424a010 — the origin of $C_c = A\,\theta^n$, the empirical law
every time-on-stream deactivation model descends from. The paper's sharper,
non-obvious claim is the one this page tests: carbon on catalyst at a given
residence time is, "within limits, independent of the hydrocarbon feed rate".

Built from the 300 ppi CCITT-G4 scan; every numeral read from cropped
native-resolution renders at digit scale (the exponents $n$ and coefficients
$A$ *are* the case). The by-line prints "ALEXIS VOORHIES, JR." — the Jr. the
catalogue omits is carried here. **No figure was digitised**: Figs. 1, 3, 4, 5
are figure-only and out of scope, so the time exponent $n$ enters only as the
author's printed fit and is *never tested against data on this page* (Tables I
and II sit at a single $\theta$ = 120 min).

## What the page shows

- **"Independent of feed rate" quantified.** Across ten groups (eight Table-I
  temperature blocks, 2× feed-rate span; cetane and Decalin in Table II, 4×
  span) the feed-rate exponent of carbon-on-catalyst is $m_{cat}$ = −0.116 to
  +0.034 — at least 9.6× closer to 0 (the paper's claim) than to +1 (the
  throughput hypothesis) in every group — while the same rows put
  carbon-on-feed at $m_{feed}$ ≈ −1 (−0.95 to −1.18). Resolving power is
  computed, not assumed: the worst raw block (−0.116) exceeds its rounding
  floor but is mostly the within-block temperature confound; corrected with
  the paper's own doubling rule, every block sits within 3.1 floors of zero,
  both signs occurring. Verdict: consistent with 0, incompatible with +1, and
  *not* "proved exactly zero".
- **The two carbon columns are one measurement, so each table carries one
  test, not two.** The implied catalyst-to-oil ratio
  $C_c/(C_f\,U\,\theta/60)$ is constant within every group to CV ≤ 2.9 % — at
  the 2-significant-figure rounding floor — consistent with both columns
  being computed from the single printed carbon determination (combustion of
  the discharged catalyst). The identity doubles as a per-row transcription
  guard (Check 1 breaks it with single-glyph perturbations), recovers a ratio
  the paper never prints for these units (per-system means 0.86–1.09, *not*
  the 1/0.58 printed for the fluid-derivation unit), and surfaces one unexplained
  ~8 % drift of that ratio across temperature blocks — reported, not
  resolved.
- **The printed algebra chain (2)+(3)+(4)→(5), exact, by two routes with
  independent solvers** (symbolic isolation in exact rationals; the same
  isolation root-found numerically + regression — their 5-decimal agreement
  checks the solvers, the printed 96 checks the shared step):
  $V = 95.7145/(U^{0.34130}\theta^{0.19113})$,
  printed as 96 / 0.34 / 0.19 — every printed digit the correct rounding.
  Every leading-digit, transposed-glyph or exponent-scale misreading of the
  five inputs (3↔8, 9↔6, dropped superscript…) misses the printed triple;
  the last digits are within the printed rounding (six of ten ±1 alternatives
  still round to it, break table 2b) and rest on the digit-scale crops, which
  were read independently twice — at transcription and at verification — in
  agreement. Bonus exact object: the chain implies instantaneous
  conversion = 0.809 × the period average, everywhere.
- **"Doubles per 190–200 °F" measured per catalyst** (goodness of fit — the
  claim was derived from Table I): natural ~206 °F, synthetic ~176 °F, two
  estimator routes each. The printed band is a fair pooled summary and sits
  outside its own table's support for either catalyst individually.
- **A printed factor-2 defect settled from the paper's own two lines**: the
  Nomenclature's "$K$ = constant $(= A^2)$" is inconsistent with eq. (8) —
  integration gives $K = A^2/2$ (symbolically, and by ODE integration which
  overshoots $A\theta^{0.5}$ by exactly $\sqrt2$ with the printed $K$).
  Consequence-free in the paper ($K$ appears nowhere else); reported, never
  repaired.
- **$A$ is local, $n$ is portable — measured from the paper's own two
  campaigns**: eq. (2) vs Table I at the same nominal conditions differ by a
  factor 1.82 in $A$ (different unit, fresh vs regenerated batches), while
  the four printed exponents stay in 0.38–0.53 — exactly the range Froment &
  Bischoff quote, reconciled against `B2.2`'s transcription in the notebook.
- **The "faulty point" remark made quantitative** (labelled inference): the
  printed 36 %→41 % pair implies via eq. (3) that the flagged run's carbon
  yield sat a factor 1.46 above the correlation line.

## Fit vs test

Nothing on this page is fitted by this page, and nothing is validation
against held-out measurement. The slope and doubling tests are consistency
checks of the paper's printed claims against the paper's own printed tables
(the doubling claim is a goodness-of-fit check — Fig. 2 was drawn from
Table I); the eq.-(5) and eq.-(8) work is exact reproduction of printed
algebra; eqs. (1), (2), (3), (6), (7) are the author's fits to figure-only
data, transcribed and never refit. The same labels appear in the notebook,
`meta.yaml`, and the data sidecars.

## Files

- `build_page.py` — generates `index.ipynb` (run `python build_page.py`, then
  execute the notebook).
- `data/voorhies-1945-tableI.csv`, `-tableII.csv` — the two printed tables,
  with per-cell provenance and the glyph notes (the two 1.1 v/v/hr readings)
  in the sidecars.
- `data/voorhies-1945-printed-constants.csv` — every printed constant and
  quantitative prose claim, with locations.
- `agreement.json` — 41 metrics for CI regression checking; break-row
  coverage asserted key-for-key in the notebook, none below `ABS_FLOOR`.

Cross-page: loads `froment-bischoff-1961-printed-claims.csv` from `B2.2` only
to reconcile the quoted exponent range (asserted equal to this page's min/max;
no borrowed number reaches any reported result).

Runtime: about 5 seconds end to end.
