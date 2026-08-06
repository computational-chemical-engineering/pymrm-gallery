# A3.4 — adversarial verification, 2026-08-02

Verifier notes on the staged page at `queue_cases/A3.4/page/`. Source read
independently: `~/papers/pymrm-gallery/Wakao1978-particle-to-fluid-transfer-CES33-1375.pdf`, pages
rendered fresh at 600 dpi, plus the Elsevier PII full text for prose.

**Verdict: send back.** The page's headline claim — that a free fit to Figure 3
rejects the printed α and β by five to six standard errors — does not survive
two analysis choices the page never mentions, one of which the extraction
already knows about and does not correct. Everything else on the page is sound
and several parts are excellent; the α/β section and the break table's
off-diagonal are the problems.

The notebook is bit-for-bit reproducible from `build_page.py` and its stored
outputs match a fresh execution exactly. Every prose number I checked against
printed output matched except the four small slips in §10.

---

## 1. The α/β refit is not robust to two unremarked choices — CONFIRMED

The page reports, on all 182 Figure 3 markers,

    log-space free fit:  alpha = 1.371 +/- 0.053   beta = 0.5459 +/- 0.0089

and calls this "five to six standard errors from the printed pair" and "sharp
enough to REJECT the printed pair over the full range". I reproduced those
numbers exactly. They are not robust.

**(a) The fitting metric.** The page fits a straight line in log–log, i.e.
equal *relative* weight. An unweighted least-squares fit of `y = a Re^b` in
**linear** y, on exactly the same 182 markers, gives

    linear-space OLS:    alpha = 1.181   beta = 0.5799

That is inside ~1.5 of the page's own standard errors of the printed 1.1 and
0.6, and it is not an exotic choice: the paper's Figure 4 — the figure whose
caption is followed by "It is seen that the data are well correlated by
eqn (12)" — plots Sh against (Sc^1/3 Re^0.6)² on **linear axes**, which is a
linear-space fit dominated by the largest values. A Theil–Sen robust fit in log
space gives 1.213 / 0.552, also intermediate. The page does not say that its
headline disagreement is a property of the estimator rather than of the data.

**(b) The axis skew the extraction already measured but did not correct.** See
§2. Correcting it moves the log-space fit to 1.339 / 0.5506.

**Both together:** `alpha = 1.157, beta = 0.5840` — about one of the page's own
standard errors from 1.1 and 0.6.

Failure scenario: a reader takes away "Wakao and Funazkri's own liquid data
reject their constants at 5–6σ". Someone re-runs it in linear space, as the
authors plainly did, and gets 1.18 / 0.58. The gallery has asserted an error in
a 1978 paper that is an artefact of a regression choice.

**What survives, and what the page should lead with.** The direction is real —
after every correction α still runs a little high and β a little low, and the
low-Re excess is visible on the printed page (I looked; the cloud at Re ≈ 4–7,
y ≈ 4–6 sits well above the drawn line, and it is not an extraction artefact).
And the *most robust number on the page is already computed and buried*:

    alpha with beta held at 0.6, whole sample:  1.1044

I could not move that number. It is 1.1040–1.1237 across every erasure-band
setting I tried, and 1.098 after the skew correction. That is the number that
should carry this section.

Also confirmed on the way: the paper's own quotes are accurate. Eq. (11) is
`eps_b J_d = 0.357 Re^-0.359 for 3 < Re < 900` (read on the 600 dpi render of
page 1383); "In liquid-phase system Sh values are large and good to be used for
the determination of α and β values"; and "In order to avoid the possible
natural convection effect the liquid-phase data for Re < 3 are not included in
our data correlation". The digitised set starts at Re = 3.03, so all 182
markers are inside the range the paper accepted — the page states this
correctly.

## 2. The standard errors are not credible as a rejection statistic — CONFIRMED

±0.053 on α assumes iid residuals. They are not, and it is measurable.

| error model | SE(α) | SE(β) | z(α) | z(β) |
|---|---|---|---|---|
| OLS, as printed on the page | 0.053 | 0.0089 | 5.1 | 6.1 |
| Newey–West HAC, L = 5 / 10 / 20 in log Re | 0.079–0.102 | 0.0123–0.0165 | 2.7–3.4 | 3.3–4.4 |
| cluster bootstrap over the page's own Re bins | 0.146 | 0.0245 | 1.9 | 2.2 |

Residual autocorrelation in log-Re order is +0.16 to +0.19 at lags 2–3 and +0.16
at lag 10 — the markers come from eleven laboratories, each contributing a
contiguous Re window with its own systematic offset, which is exactly the
structure a cluster SE is for.

Worse, the model is misspecified over the range: adding a quadratic term in
log Re gives **+0.0553 ± 0.0113, i.e. 4.9σ**. An OLS standard error on a
mis-specified mean function is not a test of "is the printed pair the truth".

And the *dominant* uncertainty is not statistical at all: the extraction-band
sensitivity in §4 spans α = 1.328–1.472, which is wider than ±0.053.

So the sentence "the liquid figure ... is sharp enough to REJECT the printed
pair over the full range" (validation check 5) must go, and the 3.8× "sharper
instrument" comparison should say that both figures' SEs are iid-optimistic.

## 3. The drawn-line control belongs on the page, and it is diagnosable — CONFIRMED

The Figure 3 sidecar reports that fitting the drawn line's own ink returns
α = 1.143, β = 0.5935 against the printed 1.1 / 0.6, and concludes "The axes
are therefore calibrated by the paper's own printed function". **Neither
number appears anywhere on the page, in `meta.yaml`, or in `A3.4.yaml`** — I
grepped. It is the page's own control experiment on an object known to be
exactly 1.1 Re^0.6, it returns +3.9 % in α and −1.1 % in β, and those are the
same signs as the accusation.

I traced the line independently and got 1.1439 / 0.5934 — the sidecar's number
reproduces. Then I found the cause. The offset of the drawn line from the
computed one is not random: −5.45 px at Re 3–10, −2.01 px at 10–100, +1.55 px at
100–1000, i.e. a **+6.8 px tilt across the plot width**, which is precisely the
~5 px vertical skew between the left- and right-axis decade ticks
(left 1144.5/807.5/471.5/135.1, right 1149.7/812.5/476.5/137.5). The sidecar
acknowledges a "~2.5 px skew" as "~0.5 % of systematic axis error" but never
propagates it into β, where it is worth Δβ = −0.0066 — the whole of the
drawn-line residual.

Re-calibrating with a column-dependent row origin from the two tick columns:

    drawn line, tilt-corrected:  alpha = 1.1162   beta = 0.5980   <- the printed pair
    markers,    tilt-corrected:  alpha = 1.3387   beta = 0.5506

So: **the line was drawn from the printed constants, the residual is
calibration, and correcting it takes about 20 % of the α gap and 9 % of the β
gap.** The drawn-line refit therefore *undermines the presentation* — it is a
measured bias in the same direction as the accusation, presented as
confirmation — while leaving the substance of a residual disagreement intact.
It must be on the page, with the correction applied or at least quoted.

## 4. The erasure band moves α by more than the quoted SE — CONFIRMED

I re-ran the Figure 3 pursuit with the computed-position band varied and
everything else fixed (`long_len = 61`, `short = 25`, same gates):

| band half-width | n markers | α | β | α (β = 0.6) |
|---|---|---|---|---|
| 0 (band off) | 210 | 1.328 | 0.5542 | 1.1040 |
| 3 | 201 | 1.348 | 0.5498 | 1.1052 |
| **5 (shipped)** | **182** | **1.371** | **0.5459** | **1.1044** |
| 7 | 162 | 1.404 | 0.5394 | 1.1068 |
| 10 | 146 | 1.457 | 0.5307 | 1.1120 |
| 14 | 133 | 1.472 | 0.5291 | 1.1237 |

The declared direction is right — erasing more of the line pushes α away from
1.1 and β away from 0.6 — but three things follow that the page does not say:

- the shipped setting is **mid-range, not conservative**. Turning the band off
  entirely still leaves the 61 px straight-ink removal doing the work, recovers
  28 more markers, and gives α = 1.328. The disclosure "relaxing the removal
  moved the count 140 → 182 and alpha 1.47 → 1.37" stops short of the setting
  that is actually available;
- the α spread across defensible settings (1.328–1.472) **exceeds the ±0.053
  standard error** the page quotes as its yardstick;
- `A3.4.yaml`'s `answer_changes.fig3_false_positives` says *"The fit is
  insensitive — alpha moved only 1.371 to 1.376 across a detection-threshold
  sweep spanning 132 to 182 markers."* That is true of the *threshold* knob and
  false of the *band* knob at the same marker count: 133 markers gives 1.472.
  A sensitivity claim stated without its knob is the handoff.md "comment
  claiming a sensitivity the check does not have" in another dress.

α with β = 0.6 is stable at 1.104–1.124 across the whole sweep. Again: that is
the number to lead with.

Two things I confirmed on the mechanism itself:
- **Don't-care masking, never ink subtraction.** `setup2.prepare` returns
  `(ink > 0.5)` unmodified with `care0 * (1 - removed)`; `setup3.prepare` does
  the same. The subtracted image `clean2` produces is discarded. Correct, and
  the trap from earlier work is avoided.
- **21 of the 182 shipped centres lie inside the ±5 px band**, so the pursuit
  does recover markers on the line rather than being blind there.

## 5. The break table's off-diagonal is not measured — CONFIRMED

`A3.4.yaml` says "nine injected defects, each moving the metric it should **and
only that one**". The table cannot show that:

```python
v = {k: (ov[k]() if k in ov else BASE[k]) for k in BASE}
```

A column not named in the defect dict is **copied from `BASE`**, not recomputed
under the defect. Every unchanged entry in that table is a copy, and the
impression of orthogonality is manufactured. This is the repository's signature
defect appearing inside the very table built to guard against it.

It happens to be true — the three metrics are genuinely decoupled — but the
table asserts it rather than demonstrating it. Fix: recompute all three metrics
for every defect, or drop the two unaffected columns and print one number per
defect.

Second point on the same table: the metric it breaks for result 5 is
`liquid_alpha`, which is the **β-fixed** α — the most defect-insensitive of the
page's three fits (1.104–1.124 across the whole band sweep). The headline
α = 1.371 has no break test at all.

Everything the table *does* claim, I re-ran and reproduce: Richardson
2.14e-03 → 3.70e-01 / −1.00e+00 / −8.74e-01; re-analysis 1.83 % → 5.81 % /
25.81 % / 3.76 %; α(liq) 1.1044 → 1.0672 / 1.0316 / 1.6286.

## 6. Citation names an author who is not on the paper — CONFIRMED

Title page, verified on the 600 dpi render and in the publisher text:
**"N WAKAO† and T FUNAZKRI"**. Two authors.

The page's closing line says *"Cite the source, not this page: Wakao, N.,
**Kaguei, S.** and Funazkri, T."*, and `page/README.md` repeats it. Kaguei is a
co-author of the 1979 heat-transfer companion, not of CES 33(10) 1375. Every
other metadata block (`meta.yaml`, `models_entry.yaml`, both sidecars) has the
two-author list correctly, so it is an isolated slip in the two places a reader
will actually copy the citation from.

## 7. The re-centring measurement is diluted, and "random" is too strong for the row shift — CONFIRMED

`FIG2_FIT_OK = 0.50` makes a low-scoring fit fall back to the 2026-07-30
crosshair. **16 of the 81 paired rows** (not twelve, as `extraction/README.md`
states) have `col_px == col_px_prev` exactly, so their displacement is zero by
construction, not by measurement. They are inside the mean and the standard
deviation.

| | page (all 81) | genuinely re-centred (65) |
|---|---|---|
| along Re | −0.01 ± 2.59 px (0.0 SEM) | −0.01 ± 2.90 px (0.03 SEM) |
| along Sh | +0.62 ± 2.19 px (2.6 SEM) | **+0.78 ± 2.43 px (2.58 SEM)** |
| median distance | 2.28 px | **3.29 px** |
| 90th pct / max | 5.36 / 9.88 px | 5.71 / 9.88 px |

The significance is unchanged, but the displacement is 25–44 % larger than
reported, and the numbers in `review/README.md` are the ones earmarked for
promotion into `docs/handoff.md`, where they will be quoted as the measured
size of the density-peak offset. Report the re-centred subset, or state that
16 zeros are included.

On the verdict itself: the column shift is genuinely consistent with zero, but
**the row shift at 2.6 SEM is a systematic, not noise** — "random, not
systematic" is the wrong label for a component three standard errors from zero.
The right statement is the one the page already has two sentences later: the
systematic component is −0.41 % in Sh, an order of magnitude below the scatter,
so nothing on the page turns on it. Keep the conclusion, change the word.

Minor, non-blocking: `fit_glyph(search=5)` plus a polish clamped at ±2.5 px caps
the measurable displacement at 7.5 px per axis. Observed max is 7.21 px in rows,
so the truncation is essentially non-binding, but it is a truncated
distribution and one point is within 0.3 px of the wall.

## 8. The 11 double counts are real; the claim about *why* they existed is not — CONFIRMED

Deduplication verified: 81 recorded positions collapse to 70 distinct
`marker_id`s (nine groups of 2, one of 3), plus 9 new − 2 audit rejects = **79
distinct glyphs**. Every merged group's fitted centres agree to under 1 px. The
merge is correctly applied everywhere the point set is used
(`obs = fig2.drop_duplicates("marker_id")`; results 4, 5, 6 and the agreement
metrics all read `obs`). Check 7 deliberately uses all 81 rows, which is correct
for a paired comparison and is labelled.

But `review/README.md` — the text going into `handoff.md` — explains the double
counts as *"Non-maximum suppression at 8 px does not [prevent it]"*. Five of the
eleven merged crosshair pairs were **1.4, 2.0, 2.2, 3.2 and 1.4 px apart** in
the 2026-07-30 CSV. Nothing 1.4 px apart survives 8 px non-maximum suppression.
So either that pass did not apply NMS as documented, or these crosshairs came
from the hand audit. Either way the causal story in the write-up is wrong and
should not be carried into the repository's institutional memory.

The maintainer's actual question — are these double counts or the F2.3 case in
reverse — is properly asked in `A3.4.yaml` and both sidecars, and I cannot
answer it from the data. The evidence for "double count" (two seeds 1–9 px apart
converging to the same sub-pixel optimum) is reasonable but is not independent
of the ±5 px search radius that lets them converge.

## 9. Two smaller calibration points

**"Three known functions, three printed curves" is two.** Validation check 4
says the Figure 2 axis calibration rests on three printed curves recomputed onto
the drawn ones. Eq. (11) contains ε_b, which the page fits to that very line
(`eps_b_figure2 ... fitted, see the page`). It constrains the calibration in
slope but not in offset, so as calibration evidence it is partly circular. The
page is honest about ε_b elsewhere ("inferred from the drawn Eq. 11 line, not
stated"); the calibration sentence just needs to say "two known functions plus
a slope".

**Figure 2's three curves all sit above their computed positions.** I traced
them: eq. (12) −3.77 ± 1.82 px, eq. (11) −4.93 ± 2.94 px, eq. (9) −2.32 ± 4.58
px (negative = drawn ink above computed). Within a line width, so the page's
claim holds — but all three have the *same* sign, which means the Figure 2
ordinate origin is ~4 px off and the extracted Sh values are ~2.6 % high. That
inflates the +9.7 % eq.-12 bias to about +12.6 %. Small, one-directional, and
worth a sentence given the page reports the bias to one decimal.

**Figure 2's dashed-curve band, unquantified.** Zero of the 79 shipped markers
fall inside the ±9 px eq. (9) band, and only 2 inside ±20 px. The band covers
2.2 % of the plot area between Re 1 and 30, where the data cloud crosses the
Ranz–Marshall curve. Zero is consistent both with chance and with the band
having erased one to three markers. The Figure 2 fit is explicitly the weak
instrument so nothing load-bearing rests on it, but the sidecar identifies this
as the one mechanism that can destroy a marker and does not report the count.

## 10. Prose against output

Four slips, all small, all real:

- Validation intro: *"Check 8 injects a defect into each of the load-bearing
  ones — checks 1 and **3** and results 3 and 5."* Check 3 gets no injected
  defect; it is in the blind-spot list as unable to fail. The table covers check
  1 and results 3 and 5.
- `meta.yaml` agreement: *"the markers below Re ~ 10, which sit **~20-30 %**
  above the line"*. The notebook prints Re 3–6 → −21.3 % and Re 6–12 → −10.2 %,
  i.e. 11–27 % above. `A3.4.yaml`'s "10-30 %" is the correct range; `meta.yaml`
  overstates the 6–12 band.
- `extraction/README.md`: *"Twelve of eighty-one did."* Sixteen did.
- `review/README.md`: *"Median displacement here was 2.3 px"* — 3.29 px on the
  markers actually re-centred (§7).

## 11. What survived my break tests

These are genuine and I could not weaken them:

- **Check 1, pymrm vs eq. (7).** Two independent routes — a discretisation of
  eqs. (4)+(5) and a closed form. Error ratios 2.00 over four doublings,
  Richardson 2.0e-06. Injecting a Dirichlet inlet, `nu = 1`, or a flipped
  dispersion sign moves it to 3.7e-01, −1.0e+00, −8.7e-01. It also tests the
  page-image transcription of eq. (7) against that of eqs. (4)+(5); I checked
  the algebra of `theta_eq7` against the printed form and it is right, including
  `4N/Pe_L = 4 a k_f D_ax/(eps_b U²)`.
- **Result 3, the paper's own recalculation.** Petrovic–Thodos eq. (11) pushed
  through eq. (7) at Pe = 2 and inverted with eq. (2) returns eq. (12) to 1.83 %
  mean / 4.48 % worst over 100 ≤ Re ≤ 900, moving to 3.8 %, 5.8 % and 25.8 %
  under three injected defects, and by under half a per cent when L/d_p is
  quadrupled. Nothing fitted, every input printed. **This is the page's best
  result** and the paper never prints the comparison.
- **Check 2**, eq. (7) → eq. (8), and the 1.33 % maximum liquid correction.
  Correct, and it is conclusion 1 of the paper made quantitative.
- **Check 3's own declaration.** "Cannot fail" is accurate: still 2.000000 at
  α = 11. Properly labelled as an identity of the chosen form.
- **The blind-spot list is honest.** I tested each: check 1 is blind to the
  dispersion coefficient, voidage and Sc because both routes are written in
  (N, Pe_L); result 5 leaves α at exactly 1.1044 whatever the bed model does;
  results 1–4 are indeed unchanged by α and β. All three verified.
- **`Sh = j_D Re Sc^(1/3)`** from `j_D = (k_c/U_sup) Sc^(2/3)` — algebra correct.
- Deviation convention `(model − measured)/measured` is stated once and used
  everywhere, on both figures.
- Housekeeping: all five overlays git-ignored (`.gitignore:42`); no `data:`
  payload in the tracked dashboards; `_site/` ignored; `extraction/` contains no
  images; data tier 4 / digitised / **provisional** matches the sibling
  experimental pages (A4.9, C2.1, F1.3, F1.4, F2.3 all tier 4) and the T1 page
  tier is right; `follow_up.blocking: false`; nothing calls reproduction
  validation; no `NumJac` is used anywhere, so the `(n, 1)` rule does not apply
  (the `(n,)` shape passed to `construct_grad`/`construct_div` is harmless);
  `agreement.json` matches the fresh run exactly; sidecars record the
  2026-08-02 review answer and the pass history correctly.

---

## Explicit verdicts

1. **Is the α/β refit stated within what the data support?** **No.** The
   direction of the disagreement is real and visible on the printed page, but
   "five to six standard errors" and "sharp enough to REJECT the printed pair"
   are the most favourable of several defensible analyses. A linear-space fit
   gives 1.181 / 0.5799; correcting the axis skew gives 1.339 / 0.5506; both
   together give 1.157 / 0.5840. The page must present the disagreement as
   estimator-dependent, and should lead with α = 1.104 at β = 0.6, which is the
   only fit I could not move.
2. **Are the quoted standard errors credible?** **No.** Residuals carry a 4.9σ
   quadratic and lag-2/3 autocorrelation of +0.16/+0.19; HAC SEs are 1.5–1.9×
   larger and a cluster bootstrap over Re bins is 2.8× larger, dropping z from
   5.1 / 6.1 to 1.9 / 2.2. And the extraction-band sensitivity (α 1.328–1.472)
   is larger than the quoted ±0.053, so the statistical SE is not the dominant
   uncertainty at all.
3. **Does the drawn-line refit undermine or support the marker refit?** It
   **undermines the presentation and partly bounds the substance**. The sidecar
   frames α = 1.143 / β = 0.5935 on a line known to be exactly 1.1 Re^0.6 as
   confirmation; it is a measured calibration bias in the same direction as the
   accusation. I diagnosed it as the 5 px left/right tick skew: correcting it
   returns the drawn line to 1.116 / 0.598 — so the line *was* drawn from the
   printed constants — and moves the markers to 1.339 / 0.5506, absorbing ~20 %
   of the α gap and ~9 % of the β gap. The residual disagreement survives, but
   the control number must be on the page.
4. **Which checks survived my break tests?** All of §11. In particular check 1,
   result 3 and the whole blind-spot list are genuine and correctly
   characterised. The break table's *diagonal* survived; its *off-diagonal* did
   not (§5).

## Fix list, if the page is to come back

Blocking:

1. Rewrite result 5 and the `meta.yaml` agreement text: add the linear-space
   fit, add the skew-corrected fit, drop "REJECT", replace "five to six
   standard errors" with a stated-error-model comparison (OLS vs cluster), and
   promote α = 1.1044 at β = 0.6 to the headline as the number that does not
   move.
2. Put the drawn-line control (1.143 / 0.5935, and 1.116 / 0.598 skew-corrected)
   on the page, not only in the sidecar, and say what fraction of the marker
   disagreement it accounts for.
3. Correct the Figure 3 y-calibration for the left/right tick skew, or state
   the uncorrected bias in α and β explicitly in the sidecar's
   `estimated_error`.
4. Recompute every column of the break table under every defect, or stop
   claiming "and only that one".
5. Print the erasure-band sweep (band 0 → 14, α 1.328 → 1.472, n 210 → 133) as
   a systematic-uncertainty table, and fix the "the fit is insensitive" line in
   `A3.4.yaml`.
6. Remove "Kaguei, S." from the citation line on the page and in
   `page/README.md`.

Non-blocking but should ride along:

7. Report the re-centring statistics on the 65 genuinely re-centred markers, or
   say that 16 structural zeros are included; correct "twelve" → sixteen in
   `extraction/README.md`; retire "random, not systematic" for the row
   component in favour of the −0.41 %-in-Sh framing already on the page.
8. Drop or rephrase the "NMS at 8 px does not [prevent double counting]" claim
   in `review/README.md` before it reaches `docs/handoff.md`.
9. "Three known functions" → two plus a slope, in validation check 4.
10. Note the ~4 px common offset of Figure 2's three drawn curves and its
    ~2.6 % effect on the extracted Sh.
11. Report how many Figure 2 markers the eq. (9) dashed band could have erased.
12. Fix "checks 1 and 3" → "check 1" in the validation intro, and "~20-30 %" →
    "10-30 %" in `meta.yaml`.
