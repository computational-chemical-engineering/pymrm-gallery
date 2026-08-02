# G1.8 — adversarial verification

Verifier pass, 2026-08-02. Target: `queue_cases/G1.8/page/`.
Source re-read independently from `~/papers/pymrm-gallery/AIChE Journal - 1983 - Herskowitz.pdf`
at 600 dpi (journal pages 4 and 8 = PDF pages 4 and 8).

**Verdict: safe to publish after the five fixes listed under "Required fixes".**
The accusation against the figure is sound and I could not break it. What needs
work is one overclaim about evidence independence, one wrong statement about
calibration, an over-precise slope number, an understated limit on check 4, and
a non-sequitur in the rejection of ALT B.

---

## What I did, so the next reader knows what is and is not second-hand

Nothing below is taken from the builder, the staged sidecar, or the case YAML.

1. Re-rendered journal pages 4 and 8 at 600 dpi and read Table 1's sphere row,
   Eqs. 19, 20, 21, the Figure 6 caption, the "Table 2" sentence and the
   notation list off the images.
2. **Re-derived Table 1's sphere row symbolically** (sympy) from Eq. 3 plus a
   Robin condition written on the *V*/*S* length, without looking at the printed
   row. Difference to the printed row simplifies to **exactly 0**.
3. **Re-digitised Figure 6 from scratch** — my own gridline detection, gridline
   erasure by position, diagonal closing, connected components, 620–640 traced
   columns per curve, my own axis calibration.
4. Re-ran the whole reconstruction, the gap statistic, both alternative
   hypotheses and the break table from my own code.
5. Re-executed the notebook (8.4 s, no errors) and diffed `agreement.json`
   against the staged one — byte-identical, so the page is deterministic.
6. Ran 19 extra break tests on check 4 and 4 extra break tests on check 2 that
   the page does not run.
7. Tested nine alternative mis-implementation hypotheses the page does not test.

---

## Confirmed correct (the things the finding rests on)

**Transcription.** Table 1 sphere row, Eqs. 19/20/21 and the caption conditions
are exactly as the page states. The "Table 2" mis-citation is real: the notation
list on page 17 reads "β, γ = constants given in Table 2", and page 12's text
introduces Table 2 as the flow-regime constants for the pressure-drop
correlation. Printed defect 1 is confirmed.

**The sphere row is derivable, not just transcribable.** Solving
∇²C = (3φ)²C on r ∈ [0,R] with (R/3)·dC/dr|_R = α(C* − C) and averaging gives

    eta = C* alpha (3 phi cosh 3 phi - sinh 3 phi)
          / ( phi^2 [ 3 alpha sinh 3 phi + 3 phi cosh 3 phi - sinh 3 phi ] )

which is the builder's `9 alpha C* G / (lam^2 [G + 3 alpha sinh lam])` and
simplifies to Table 1's row **identically**. This is what pins φ to the *V*/*S*
length; it is not an assertion.

**The digitisation.** My independent trace against the CSV:

| curve | my χ(φ=10) | CSV | diff | my slope | CSV slope |
|---|---|---|---|---|---|
| 1 | 3.7851 | 3.7515 | +0.89 % | 1.0991 | 1.1051 |
| 2 | 1.4406 | 1.4412 | −0.04 % | 1.0879 | 1.0976 |
| 3 | 0.3561 | 0.3543 | +0.51 % | 1.1010 | 1.1122 |
| 4 | 0.1432 | 0.1429 | +0.20 % | 1.1233 | 1.1311 |

My fit rms 0.021–0.027 decades against the CSV's 0.024–0.029. My ordinate
calibration is 415.90 px/decade with rms **0.0023 decades** over 3.3 decades
(sidecar says 415.5).

**No missing curve.** The F3.5 failure mode is ruled out: over 650 usable
columns the run-count histogram is 4 in 550 columns and 3 in 90 (gridline
overlaps), never a systematic 5th. No curve leaves the frame through the top or
the bottom — curve 1 tops out near χ = 20 against a frame at χ ≈ 40, curve 4
enters at the left frame at χ ≈ 0.0105 against a frame at 0.01.

**The reconstruction, from my own numbers rather than the CSV:**
L = **0.5007, 2.047, 7.041, 10.170** — i.e. 0.14 %, 2.35 % and 0.59 % from the
printed 0.50, 2.0, 7.0. The builder's 0.5082 / 2.0460 / 7.0609 / 10.1755
reproduce exactly from the CSV. Deviations at φ = 10 from my digitisation:
+0.1 / +69.3 / +312.2 / +151.3 % as printed and +0.09 / +1.90 / +1.04 / +6.20 %
shifted, against the page's +1.0 / +69.2 / +314.3 / +151.8 and
+0.98 / +1.85 / +1.56 / +6.42.

**Tolerance propagation** δL/L = (δχ/χ)/|d log χ/d log L| is correct;
sens −0.60/−0.81/−1.79/−3.68 and tolerances 11.4/8.3/3.8/1.5 % reproduce. It is
**conservative, not flattering**: it charges the whole fit rms as the
uncertainty on the fitted line's value at the *centre* of its own window, which
over-states it by up to √N if residuals were independent. The actual residuals
in χ are 0.98 %, 1.87 %, 1.55 % (0.004–0.008 decades) — far below the fit rms
and comparable to the difference between two independent digitisations. The
agreement is real, not manufactured by a generous band.

**Factorisation.** The fitted model slope is 1.076670 for L_m = 0.5, 1, 2, 7 and
10 to six decimals. Reassignment moves it by exactly zero. Confirmed.

**Powerless-check honesty.** The collapse identity holds at 1e-15 under a sphere
term corrupted ×1.01/×1.5/×3 while check 2 moves to 6.1e-3 / 3.0e-1 / 1.2 —
reproduced. `maxfev=1` at m = 1 changes η by 2.6e-9 — reproduced. Both are
measured, both are labelled, and the page names what they cannot test.

**Rivals.** ALT A 1057 / 16.2 / 4.6 / 9.2 (spread 230×) and ALT B
+0.3 / +9.4 / +13.3 / +3.9 pp — both reproduced exactly. The builder's
correction of the staged sidecar is right: the *correction* is non-monotone
(the required f_e values themselves are monotone, which is what the old note
was describing).

**Routine.** Nine required sections in order; no Quarto-only markdown; sidecar
carries all four `SIDECAR_REQUIRED` fields; deviation is (model − figure)/figure
everywhere; tier 6 with `method: digitised`; nothing calls the curves
experimental or "validated"; `NumJac((n_u, 1))` with no `axes_diagonals`;
`nu=2` commented; outward-normal BCs with the physical equation in a comment;
operators cached outside Newton; grid order reported as a range 2.13 → 2.84;
all four named blind spots present. **G1.8 is at `models.yaml:883` with
`status: planned`, slug `trickle-bed-partial-wetting`** — the builder's
correction of the dispatch brief is confirmed, and `models_entry.yaml` is an
in-place upgrade with matching slug and title.

---

## Findings, by severity

### 1. "Two independent signatures" is false — CONFIRMED (moderate)

`build_page.py` L1082–1085 ("Two independent signatures agree with it — the
curve positions … and the scale-invariant curve spacings"), `meta.yaml`
("reconstructed independently twice"), `README.md` ("Independent support").

`gaps_fig = np.log10(CHI10[:-1] / CHI10[1:])` is computed from **the same four
`CHI10` numbers** that the root-find inverts. There is no second measurement.
The gap statistic is the position statistic with one degree of freedom removed —
a common multiplicative offset in χ.

What it genuinely establishes, and should be labelled as: *no uniform rescaling
of the figure's χ axis can rescue the printed reading.* That is worth saying and
worth keeping. It is not corroboration by an independent witness.

Failure scenario: this sentence is exactly the kind of background-prose claim
that travelled from `F1.4` into `F1.3`. If the four fitted positions were
systematically wrong, both "witnesses" fail together and the page says two
independent things agreed.

### 2. The gap witness is not calibration-free, and it is φ-dependent — CONFIRMED (moderate)

`build_page.py` L531–534: "needs no axis calibration at all … they survive any
rescaling of either axis."

A gap measured in decades scales **linearly** with the assumed
decades-per-pixel on the ordinate. It is invariant to the ordinate *origin* and
to the abscissa entirely — which is what the page's own later sentence ("visible
without knowing where either axis starts") correctly says. Two sentences in the
same section say different things and the stronger one is wrong.

In substance this is harmless: I calibrated the ordinate independently to
415.90 px/decade with rms 0.0023 decades, so the gaps are good to ±0.001
decades and the discrimination (0.415/0.609/0.394 against 0.412/0.611/0.374
shifted and 0.191/0.221/0.611 printed) is untouched.

Second part: the *figure's* gaps are not φ-independent, because the four fitted
slopes differ by 3 %. From the CSV lines they run 0.410/0.620/0.408 at φ = 2 and
0.420/0.600/0.382 at φ = 45.5. The model's are φ-independent. The page quotes
the φ = 10 values without saying that is where they are evaluated.

### 3. The slope residual's magnitude is over-precise, and its stated defence is powerless against the dominant error — CONFIRMED (moderate)

My independent trace gives slopes 1.0991/1.0879/1.1010/1.1233, mean **1.1028**,
against the CSV's mean 1.1115. The ratio 1.1115/1.1028 = 1.0079 is *exactly* the
ratio of the two abscissa calibrations (421.9 vs my 418.6 px/decade), so the
entire difference is calibration, not tracing.

The abscissa is the weak axis. The four labelled verticals (φ = 1, 5, 10, 50) do
not sit on a common log ruler: a least-squares fit leaves **+0.013 decades** at
φ = 10. Drop φ = 10 and the other three are consistent at 417.8 px/decade with
rms 0.002 — as clean as the ordinate. So the φ = 10 gridline is misplaced in the
artwork by ~4 px, and the sidecar's 421.9 is ~1 % above the best-supported value.

| abscissa calibration | mean measured slope | gap vs 1.0767 | end-to-end |
|---|---|---|---|
| 417.8 px/dec (3-point, φ=10 dropped) | 1.1049 | +0.028 | +11 % |
| 418.6 px/dec (4-point LSQ, mine) | 1.1028 | +0.026 | +10 % |
| 421.9 px/dec (the CSV's) | 1.1115 | +0.035 | +14 % |
| 426.0 px/dec (φ=1→10 only) | 1.0837 | +0.007 | +3 % |

The page's defence — "the four measured slopes are systematically *above* 1.077,
not scattered around it" — has **no force against this error class**, because
all four slopes share one abscissa calibration and a common-mode error moves all
four together. This is structurally the `F1.4` confound: four numbers that look
like four samples are one measurement.

The sign survives every wide-baseline calibration, so "reported as unexplained"
remains the right verdict. The *number* should be a range, not +0.035 / +14 %.
Note `agreement.json` pins `slope_gap = 0.03483` as a CI baseline.

### 4. Check 4's declared limit is understated — CONFIRMED (low–moderate)

Baseline 2.30 %, tolerance band 3.8–11.4 %. The page names one non-rejected
defect (Eq. 21 6.91 → 6.31 at 4.14 %). I found four more, and one that makes the
metric *better* than baseline:

| injected defect | check 4 | rejected? |
|---|---|---|
| Eq. 21 6.91 → **7.91** | **1.72 %** | no — *better than baseline* |
| Eq. 21 6.91 → 6.71 | 2.71 % | no |
| Eq. 21 exponent 0.6 → 0.5 | 3.33 % | no |
| Eq. 21 exponent 0.6 → 0.7 | 3.43 % | no |
| Eq. 21 6.91 → 6.31 (the page's) | 4.14 % | marginal |
| Eq. 21 1.05 → 1.03 | 4.48 % | marginal |
| Eq. 21 1.05 → 1.25 | 20.1 % | yes |
| Eq. 20 0.77 → 0.75 / 0.80 | 25.0 / 26.5 % | yes |
| Eq. 20 exponent 0.1 → 0.08 / 0.12 | 47.2 / 23.5 % | yes |

The reason is structural and belongs on the page: at L_m ≈ 1 the
1/(1.05 L^0.3) term carries **87 %** of 1/α_gLs, so the 6.91 L^0.6 term barely
moves χ. So the page's "check 4 catches gross transcription errors in
Eqs. 20–21" is too generous — it catches errors in Eq. 20 and in Eq. 21's *first*
term, and is essentially blind to the second. And the summary line "every
injected defect moved it, by factors 2× to 1302×" would need qualifying if
6.91 → 7.91 were in the table.

Related: check 4 is also weak against an ordinate-scale error — corrupting the
whole χ calibration by ±5 % leaves it at 6.3 % / 10.6 %, inside the band. That
is precisely the error class the gap statistic *does* cover, which is the honest
way to justify keeping the gap test (see finding 1).

### 5. The ALT B rejection is a non-sequitur, though the conclusion is right — CONFIRMED (low)

"No fixed offset, no fixed factor and no rewriting of Eq. 20's coefficient or
exponent can produce that shape, since 0.77 L_m^0.1 is monotone in L_m for any
coefficient and any exponent."

A non-monotone *correction* does not preclude a monotone power law: c·L^n −
0.77·L^0.1 can rise then fall for perfectly ordinary c and n. And the required
f_e values themselves (0.7212, 0.8641, 0.9578, 0.9743) **are** monotone.

What actually kills a power law is curvature. Their log–log slopes over the
three intervals are 0.260, 0.149, 0.014 — strongly concave. The best single
power law is 0.780 L^0.114, which misses the middle two required values by
−9.7 % and −11.8 %. The page's *second* reason ("the corrected values would
contradict the f_e column printed inside the figure") is sufficient on its own,
so the verdict stands. Only the stated reason needs replacing.

### 6. Small internal inconsistencies — CONFIRMED (low)

- α_gLs range quoted twice with different values: the α_gs sweep prints
  "α_gLs itself is only **0.72 to 1.95** over the four curves" (shifted reading,
  L = 0.5…10) and the ALT A paragraph prints "**0.72–1.74**" (printed reading,
  L = 0.5…7). Both correct in context; the first is labelled "over the four
  curves" and refers to a different four.
- `README.md`'s break table gives check 2's baseline as 7.7e-06 (Robin at
  L_m = 2 only) while `meta.yaml`'s headline is 3.1e-4 (the worst case,
  Dirichlet at φ = 45.5 — which is the case actually used for the dry path).
  The break factors "1.2e3×–2.6e5×" are computed against the smaller of the two;
  against 3.1e-4 the n_u = 50 break is 30×.
- "Roughly a third of its excess is the rounding of 10.18 to 10" — the printed
  numbers are 14.6 % → 11.4 %, a 22 % reduction. This is hardcoded prose, not
  interpolated, and it is the only prose/output mismatch I found on the page.
- "Both factors of 3 are asserted, not printed, so both are broken on purpose in
  Validation" — there is one knob (`n_scale`) and it moves both together. I
  broke them independently against a 5.5e-5 baseline: interior only 3 → 1 gives
  **7.6e0**, boundary only 3 → 1 gives **6.6e-1**. So the claim is true in
  substance; the code does not do what the sentence says.

### 7. Stale, now-wrong G1.8 numbers survive in two tracked docs — CONFIRMED (integration, informational)

Not the builder's to fix (`docs/` was out of scope), but they must not survive
integration or the next agent re-derives the wrong thing.

`docs/handoff.md` §"G1.8 is blocked" and `docs/staged-data/G1.8/…meta.yaml` say:

- α_gs = **2462**/16.2/4.7/10.3 — correct is 1057/16.2/4.6/9.2. The old run used
  the figure's *rounded* f_e = 0.72 instead of Eq. 20's 0.7184.
- "a correction that **grows with L_m**" — it is non-monotone.
- "log-log slope **1.03** against 1.10–1.13" — that compared a *local* slope near
  φ = 10 to *global* fits. The like-for-like global fit is 1.0767. The builder's
  correction of the dispatch brief on this point is confirmed.
- "3.767 against 3.752 … **0.4 %**" — with Eq. 20 it is 3.7884, i.e. +1.0 %.
- "the legend rows align with the curve endpoints top to bottom, so the
  curve-to-L_m assignment **is as printed**" — now the opposite of the published
  finding.
- the parallelism note lists slopes "1.098, **1.098**, 1.112, 1.131"; the CSV's
  first two are 1.1051 and 1.0976.

### 8. "The legend block is offset by one row" describes something the artwork does not literally show — CONFIRMED (presentational)

On the 600 dpi render the (L_m, f_e) table is a boxed inset whose four rows are
each positioned at the height where the corresponding curve *terminates* at the
box's left edge (φ ≈ 49). There is no free-floating legend block that could have
slipped; the draughtsman put each label at the end of a curve. The page is
careful to say it cannot distinguish legend-slip from plotting error, so it is
not wrong — but a reader who goes looking for a misplaced block will not find
one. The neutral statement is: *each drawn curve carries a label; from the second
curve down each label is the previous curve's L_m; the value 1.0 labels no curve;
and the fourth curve's value (≈10) appears nowhere in the legend.*

Related: "the whole pattern is displaced by one position" (the gap narrative) is
loose. Only **one** of the three gaps (2 → 7 = 0.611) is shared between the two
readings; the figure's other two gaps correspond to 0.5 → 2 and 7 → 10, which do
not appear in the printed model's gap sequence at all. The evidence is real; the
"displacement" story about the gap *sequence* is rhetorical.

---

## Attacks that failed — the finding survived all of them

**Nine alternative mis-implementations**, each asked to rescue the printed
reading. χ(φ=10) deviations against the four figure curves:

| hypothesis | dev % | figure gaps 0.415/0.609/0.394 reproduced? |
|---|---|---|
| as printed (baseline) | +1 / +69 / +314 / +152 | 0.191/0.221/0.611 — no |
| f_e(L_i) with α_gLs(L_{i+1}) | −20 / +34 / +175 / +124 | 0.189/0.298/0.482 — no |
| f_e(L_{i+1}) with α_gLs(L_i) | −18 / +29 / +53 / +19 | 0.223/0.533/0.503 — no |
| Eq. 21 gas–liquid term only | −15 / +47 / +269 / +132 | no |
| Eq. 21 liquid–solid term only | −84 / −78 / −55 / −80 | no |
| Eq. 21 reciprocal mis-read | −48 / +41 / +449 / +658 | no |
| φ on R instead of R/3 | +182 / +373 / +1057 / +603 | no |
| χ normalised by η_o | −79 / −51 / +68 / +85 | no |
| L_m rescaled ×0.1, ×0.5, ×2, ×10 | all ≥ 35 % somewhere | no |

**A stronger, more general rejection than the page's, worth adopting.** Every
variant of Eqs. 20/21 separates as χ = g(L_m)·h(φ). The figure's own gaps then
force

    d log g / d log L_m  =  -1.380, -2.024, -0.725   over [0.5,1], [1,2], [2,7]

— non-monotone, and steeper in the middle than at either end. No product of
powers of L_m can do that. Under the shift the *same* figure forces

    d log g / d log L_m  =  -0.690, -1.120, -2.546   over [0.5,2], [2,7], [7,10]

against the model's −0.684, −1.122, −2.414. This is one line of arithmetic, it
needs no per-curve fitting, and it kills ALT A, ALT B and every
unstated-parameter hypothesis simultaneously. It would be a stronger closing
argument than the two fitted rivals.

**A fifth curve.** Ruled out (see above).

**Circularity.** Check 2 compares a finite-volume solve to Table 1's closed form
and shares no code with it; I verified Table 1 independently, so check 2 tests
the implementation. Check 1 compares Eq. 20 to the paper's own printed f_e
column — a transcription check, correctly labelled. Check 4 *is* the finding
presented as a check, which is legitimate but means "check 4 passes" is not
independent evidence for the shift; the page does not claim it is.

**Determinism.** Re-executed `agreement.json` is byte-identical to the staged
one.

---

## Required fixes before publishing

1. **Drop "independent"** from the description of the gap statistic (page prose,
   `meta.yaml`, `README.md`). Replace with what it actually shows: it removes a
   common multiplicative offset in the χ axis, so it rules out the hypothesis
   that the printed reading is right and the χ calibration is wrong.
2. **Fix the calibration claim.** Delete "needs no axis calibration at all" and
   "survive any rescaling of either axis"; keep "visible without knowing where
   either axis starts". Say the gaps need only the ordinate's decades-per-pixel,
   which the sidecar reports with residuals under 0.006 decades. Say the figure's
   gaps are evaluated at φ = 10.
3. **Quote the slope gap as a range and name its dominant uncertainty.** The
   abscissa calibration, not the fit rms, dominates: an independent
   re-digitisation gives +0.026 (+10 %) against the page's +0.035 (+14 %), and
   plausible calibrations span +0.007 to +0.060. Remove or qualify "the four
   measured slopes are systematically above … not scattered around it" — the four
   share one calibration, so that is not a scatter argument.
4. **Extend check 4's declared limit.** It is blind to the 6.91 L^0.6 term of
   Eq. 21 (that term carries only ~13 % of the resistance at L_m ≈ 1), and one
   corruption of it (6.91 → 7.91) *improves* the metric. Also note the ±5 %
   ordinate-scale blind spot, which is the honest reason to keep the gap test.
5. **Replace the ALT B argument.** Non-monotonicity of the correction does not
   preclude a power law. Use the curvature (log–log slopes 0.260 → 0.149 → 0.014;
   best power law 0.780 L^0.114 misses by −9.7 % and −11.8 %) or lead with the
   argument that already works — the corrected f_e values contradict the column
   printed inside the figure.

## Recommended (not blocking)

6. Fix the four small inconsistencies in finding 6 (the two α_gLs ranges, the two
   check-2 baselines, "roughly a third" → 22 %, and the "both factors of 3 are
   broken" sentence — or break them independently, which I measured moves check 2
   by 1.4e5× and 1.2e4×).
7. Re-word "the legend block is offset by one row" per finding 8, and soften the
   "whole pattern displaced by one position" gap narrative.
8. Consider adopting the separability argument above as the closing paragraph of
   "The alternative explanations, and why they lose".

## Integration must-do

9. `docs/handoff.md` §"G1.8 is blocked" and `docs/staged-data/G1.8/…meta.yaml`
   carry six now-wrong numbers and one now-reversed conclusion (finding 7). Fix
   or delete both when the page is integrated.

---

## Explicit verdicts

1. **Is the shift reconstruction sound enough to publish as an accusation
   against a published figure?** **Yes.** I verified it three ways the builder
   did not supply: an independent symbolic derivation of Table 1's sphere row
   (exact zero), an independent 600 dpi re-digitisation of all four curves
   (χ at φ = 10 within 0.9 %, reconstructed L_m within 2.4 % of the printed
   0.5/2.0/7.0), and my own reading of Eqs. 19–21 and the Table 2 slip. Nine
   alternative mis-implementations and a general separability argument all fail
   to rescue the printed reading. The accusation is earned.

2. **Is the scale-invariant witness genuinely calibration-independent?** **No,
   as stated** — it requires the ordinate's decades-per-pixel, and it is not
   independent of the position reconstruction (same four numbers). In substance
   it survives, because I verified the ordinate to 0.1 % and the discrimination
   between the two readings is a factor of 2–3 in every gap. Keep the test,
   relabel it.

3. **The residual-slope framing.** The builder's correction of the dispatch
   brief (1.0767 like-for-like, not 1.03 local against 1.10–1.13 global) is
   right, and the factorisation claim is exactly right — every L_m gives
   1.076670 to six decimals, so reassignment moves the slope by zero. The
   magnitude is over-precise: the abscissa calibration, which the page never
   considers, is the dominant uncertainty and takes the gap from +0.035 down to
   +0.026. The "systematically above, not scattered" counter-argument is
   powerless against it. Direction and "unexplained" stand; the number needs a
   range.

4. **Which checks survived my own break tests.** Check 1 (Eq. 20) — survived,
   25–436×. Check 2 (pymrm vs Table 1) — survived, and survived two breaks the
   page does not run (the interior and boundary factors of 3 independently,
   1.4e5× and 1.2e4×). Check 3 (grid) — survived, order 2.13 → 2.84. Check 4
   (reconstruction) — survived only for Eq. 20 and Eq. 21's 1.05 term; blind to
   the 6.91 term (four non-rejected defects, one that improves the metric) and
   weak against a ±5 % ordinate-scale error. Check 5 (collapse) — confirmed
   powerless, correctly labelled and correctly measured. Check 6 (α_gs sweep) —
   reproduced, α_gs = 1000 costs 2.99 % at φ = 30.

**Safe to publish after fixes 1–5.**
