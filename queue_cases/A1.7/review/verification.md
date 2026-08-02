# A1.7 — adversarial verification

Verifier pass, 2026-08-02. Source read independently at 600 dpi
(`pdftoppm -r 600`, all 8 pages). Notebook re-executed from the staged state
(4.2 s, outputs bit-identical to what is committed); `build_page.py` regenerates
`index.ipynb` sources exactly.

Nothing was written outside this file.

---

## Verdicts asked for

1. **Table 1 is genuinely measurements and the tier is right.** Confirmed.
   Journal page 288, "Experimental results on Group A powders". §4.2 (page 287)
   states the method for the two columns the page's headline rests on: "The
   minimum fluidization velocity was measured in the usual way using the
   pressure drop–gas velocity curve" and U_MB by raising the air velocity until
   the first recognizable bubble broke the surface, repeated several times per
   fraction and averaged. `d_sv` is from a microscope count of ≥650 particles.
   All 22 rows × 6 fields transcribe **exactly** against the 600 dpi render,
   including the em-dash `U_0`, the one-significant-figure `0.4`, the anomalous
   `75–90 → 100 µm` row, and the `ρ_s =` / `ρ_s ≃` distinction. Tier 2 in
   `docs/data-strategy.md` is "tables printed in papers", which this is
   independently of measurement status. **Caveat:** the paper says nothing about
   how `ε_MB` was obtained; it is almost certainly derived from bed height, the
   200 g charge and ρ_s (approximate for both catalysts). Calling it "measured"
   is unsupported — see F5.

2. **The 20/21 agreement is column-independent but nearly powerless, and rests
   on an undisclosed tie-break.** The two routes genuinely share no input column
   (criterion: `d_sv_um`, `drho`; measurement: `U_MB`, `U_0`). But see F1 and F2:
   a constant "everything is group A" predictor scores **19/21** against the
   criterion's 20/21, and the two rows that produce the margin have
   `U_MB/U_0 = 1.000` exactly and are scored group B only because the code uses
   a strict `>` where the page prints eq. (4) as `≥`.

3. **The B/D = Davidson cloud identity is real algebra, not coincidence.**
   Re-derived from the printed eq. (7) without reference to the builder's
   rearrangement: LHS `(g d_B/2)^½ = 0.70711 √(g d_B)`; RHS is eq. (3) divided
   by ε₀, i.e. `U_0/ε_0 = u_f`. So eq. (7) is `u_br ≤ u_f` at d_B = 25 cm, which
   is E1.2's no-cloud pole. `0.711/0.70711 − 1 = 0.5506 %`, against the 1.5928 %
   Geldart rounded away — both reproduced. The identity is correct. Its novelty
   is overstated slightly (F7).

4. **Break tests that survived my own attack:** every number in the three break
   tables reproduces; the eq. (6)/(8) recomputations, the eq. (3) span, the ε₀
   range, the eq. (9)/(10) constants, the worked example, the eq. (2) fit
   residual and the cloud/closed-form agreement all reproduce to the digit. The
   declared blind spots are all real (`ρ_f`, the eq.-(3) spread, the +0.10 ε_MB
   error, the structural H/H₀ = 1 identity). The one break table that is
   *missing* a row is break table 2 — see F1.

---

## Findings, ranked

### F1 — CONFIRMED, highest severity. The headline 20/21 hangs on an undisclosed strict-inequality tie-break, and flipping it to the form the page itself prints inverts the result.

Two Diakon rows print `U_0 = U_MB` exactly (263 µm: 3.11/3.11; 318 µm:
4.11/4.11), so `U_MB/U_0 = 1.000`. Cell 13 uses `ratio_meas > 1.0`. The page's
own "The published model" section prints eq. (4) as `U_MB/U_0 ≥ 1`, which is
what journal page 289 prints.

Measured with the printed `≥`:

| tie-break | criterion, all rows | criterion, Diakon | all-A null model |
|---|---|---|---|
| `>` (what the code does) | 20/21 | 7/8 | 19/21 |
| `≥` (what eq. (4) prints) | 18/21 | 5/8 | 21/21 |

The `>` reading is the *correct* one — §4.1 defines group A as
`U_MB/U_0 > 1`, and §4.3 says those two fractions "bubbled at the incipient
fluidization velocity" and "should be classified as belonging to group B" — so
no number needs changing. But the choice is load-bearing, it is not stated
anywhere, and break table 2, which enumerates nine other plausible mis-readings,
does not contain it.

*Failure scenario:* a reader reimplements from the printed eq. (4) with `>=`,
gets 18/21, and concludes the page fabricated its headline. Or a maintainer
"fixes" the code to match the printed `≥` and silently destroys the result.

**Fix:** state in the results section that the two rows printing `U_MB = U_0`
are scored group B on §4.3's authority; add a `criterion applied with ≥ instead
of >` row to break table 2 with the 18/21 and 5/8 it produces.

### F2 — CONFIRMED, high severity. No null baseline is reported, and the criterion beats one by exactly one row.

19 of the 21 comparable rows are measured group A and 18 of them are predicted A
by any plausible boundary. A constant "group A" predictor scores **19/21**
overall and **6/8** on Diakon. The criterion scores 20/21 and 7/8. The entire
discriminating set is three rows — Diakon 220, 263, 318 µm — on which the
criterion is right **2 of 3** (and 0 of 3 under F1's alternative tie-break).

The page gestures at this ("14 of the 22 fractions sit so far inside group A
that the whole-table score is insensitive") but never quantifies it, and
`agreement.json` exports only `AB_rows_agreeing = 20`, `AB_rows_compared = 21`.
Break table 2 shows the symptom — a 10× mis-read of the boundary (225 → 2250)
still scores 19/21 — without naming the cause.

*Failure scenario:* the 20/21 is quoted downstream as strong validation of the
A/B boundary. It is worth exactly one row of information over a trivial
baseline, and this is precisely the "agreement that is almost guaranteed"
pattern `handoff.md` warns about.

**Fix:** print the all-A baseline (19/21, Diakon 6/8) beside the score, state
"2 of the 3 rows that can discriminate", and add `AB_rows_null_baseline` and
`AB_rows_agreeing_diakon` to the exported metrics.

### F3 — CONFIRMED, medium severity. "The criterion and its author disagree in exactly the same place" is backwards.

On the 220 µm Diakon row: eq. (6) says **B** (product 259); the strict ratio
test on the measured velocities says **A** (1.029); and §4.3 says that fraction
"should be classified as belonging to group B" — i.e. **Geldart's prose sides
with eq. (6)**, against the mechanical ratio test.

The page (cell 36, "What pymrm adds") writes: "the criterion and its author
disagree with each other exactly once, in the same place". `meta.yaml`
validation bullet 2 and `A1.7.yaml` `key_numbers` say the same. On the natural
reading — "the criterion" = eq. (6) — this is false, and it throws away a
stronger and true result.

Cell 14 gets the substance right ("He is reading his own criterion loosely, and
correctly"), so this is a sentence-level defect in the summary, not in the
analysis.

**Fix:** rewrite as "on the one row where eq. (6) and the measured ratio part
company, Geldart's own prose sides with eq. (6)".

### F4 — CONFIRMED, medium-low severity. The eq. (9) placement check does not reproduce the statement it claims to.

Cell 28 reports that eq. (9) crosses line XY at d' = 16310 µm, "a factor 16
beyond the largest labelled tick on his Figure 3 abscissa", and concludes "that
reproduces both of Geldart's statements", the first being that eq. (9) "lies
much too far to the right and is not shown on Fig. 3".

The crossing point sits at Δρ = 0.014 g cm⁻³ — fourteen times below the bottom
of Fig. 3's ordinate (0.2) — so it says nothing about what is on the diagram.
Where the eq. (9) line actually falls inside Fig. 3's window
(20–2000 µm, 0.2–7 g cm⁻³):

| Δρ | eq. (9) | line XY |
|---|---|---|
| 7 | 256 µm | 32 µm |
| 1 | 938 µm | 225 µm |
| 0.2 | 2743 µm | 1125 µm |

Eq. (9) would be plainly visible over most of the plotted region. "Much too far
to the right" is relative to XY and to the data (a factor ≈ 4.2 in size at
Δρ = 1), not off-scale.

I confirmed on the 600 dpi render of Fig. 3 that only four curves are drawn —
band PQ, XY, "EQUATION 8", and O–O — so eq. (9) is indeed absent; and I checked
XY and eq. 8 against the printed axes (Y sits at ≈1000 µm, 0.22 g cm⁻³ →
product 221 ≈ 225; eq. 8's ends give 9.4–9.8 × 10⁵ ≈ 10⁶), and O–O has slope
−1.5 as the d'^1.5 form requires. So the transcriptions are right; only the
justifying sentence is a non-sequitur.

*Failure scenario:* a reader concludes eq. (9) is trivially off-scale rather
than that it misclassifies group B powders as group A.

**Fix:** replace the crossing-point argument with the in-window comparison
(e.g. at Δρ = 1, eq. (9) sits at 938 µm against XY's 225 µm), or drop the claim
to "eq. (9) lies to the right of XY throughout the plotted region, consistent
with Geldart's remark".

### F5 — CONFIRMED, low severity. `ε_MB` and `H_MB/H_0` are labelled "measured"; the paper does not say so.

§4.2 describes the measurement of `U_0` and `U_MB` only. `ε_MB` is not a directly
measurable quantity — it follows from bed height, the 200 g charge and ρ_s,
which is printed as `≃1` and `≃1.5` for the two catalysts. The CSV sidecar and
cell 6 both call all four columns measured.

Consequence is contained: the page uses these columns only for a transcription
check it explicitly labels as such. But the inferred-ε₀ range (0.44–0.57) that
the page calls "physically admissible" inherits the approximate ρ_s for 14 of
22 rows, and the page does not say so.

**Fix:** in the CSV sidecar and cell 6, describe `ε_MB`/`H_MB/H_0` as
"reported"; note that the paper does not state how `ε_MB` was obtained.

### F6 — CONFIRMED, low severity. Blockquote is not verbatim.

Cell 14 quotes "(U_MB/U_0 > 1 for group A)". The printed text is "(U_MB/U_0 > 1
for group A **powders**)". One dropped word inside quotation marks. Everything
else I checked verbatim — the §4.1/3.4/5 quotations, "lies much too far to the
right and is not shown on Fig. 3", "agrees with published results about as well
as eqn. (5)", "drawn so as to separate empirically open and half-closed
points", the "arbitrary (though reasonable) choice of d_B" and the group C/D
descriptions — all match the render.

### F7 — PLAUSIBLE, low severity. The cross-page identity is real but its novelty is oversold, and the page does not say its cloud numbers move with E1.2.

Two small overstatements:

- "Neither paper says this." Geldart's §3.4 and page 289 both say in words that
  the criterion is a bubble rising more slowly than the interstitial gas, and
  the slow-bubble/no-cloud equivalence is standard Kunii–Levenspiel material —
  E1.2 itself uses "slow bubble: no cloud exists" language. What is genuinely
  not in either document is the *arithmetic* identification of eq. (7) with the
  pole. Cell 21's "the whole Geldart diagram is a contour map of the smallest
  cloud-bearing bubble" and cell 36's "something it is not usually read as"
  should be softened by one clause.
- The E1.2 reuse is credited in four places (cell 8, cell 9 docstring, cell 33
  item 5, Reuse), which is thorough. But the page never states the dependency
  the case file records: **if E1.2's operator or closed form changes, this
  page's cloud numbers move with it.** One sentence.

### F8 — CONFIRMED, cosmetic.

- `meta.yaml` has `status: published`; the convention for a staged, unverified
  page in `queue_cases/*/page/meta.yaml` is `in-progress` with the integrator
  flipping it. (Two other queue cases also use `published`, so this is a
  convention drift, not a rule breach.)
- `meta.yaml`/`models_entry.yaml` claim "no row of Table 1 falls in the gap
  **either** creates". True for eq. (8) (max Δρ d'² in Table 1 is 1.19 × 10⁵,
  an order of magnitude below 10⁶) but computed nowhere; only the eq. (6) gap is
  in `agreement.json`.
- Cell 33 item 2 and the README say eq. (8) "starts around 1000 µm at that
  density"; at 1.5 g cm⁻³ it starts at 816 µm.
- Two unrelated quantities both come out at 1.94 % (the eq. (6) rounding and the
  worked example's worst velocity); the README table lists them adjacently.

---

## Checks that passed, with what I did

- **Every Table 1 cell** against three high-resolution crops of the printed
  table. Exact, all 22 rows.
- **Constants and nomenclature** against the list of symbols, journal page 291:
  `g = 981 cm s⁻²`, `d' = cm × 10⁴`, ε₀ = bed voidage at minimum fluidization,
  ρ_s including internal porosity. µ = 1.8 × 10⁻⁴ on pages 289 and 290. All as
  the page's parameter table says.
- **References 10, 11, 36** against journal page 292: Verloop & Heertjes,
  *Chem. Eng. Sci.* 25 (1970) 825; L. Davies & J. F. Richardson, *Trans. Inst.
  Chem. Engrs.* 44 (1966) T293; R. D. Oltrogge, Ph.D. thesis, Univ. of Michigan,
  1972. `models_entry.yaml`'s `origin_not_consulted` block is exactly right, and
  the eq. (3) framing ("used exactly as Geldart prints it", "every statement
  about the correlation is a statement about the form printed on journal page
  289") is careful enough that the factor-6.4 result reads as a test of
  Geldart's printed eq. (3) against Geldart's own data, not as a verdict on
  Davies & Richardson. **This priority is clean.**
- **Eq. (6) recomputation, independently:** `100 × 1.8e-4 / (8e-4 × 981 ×
  1e-4) = 229.3578`, +1.9368 % on 225. Rows in the 225–229.36 gap: **0**
  (nearest are 223.97 and 259.34).
- **Eq. (8) recomputation, independently:** `√(981 × 25/2) = 110.7361`;
  `110.7361 × 1.8e-4 × 0.4 / (8e-4 × 981 × 1e-8) = 1.015928e6`, +1.5928 %.
- **Worked example, journal page 289:** 0.4360 vs 0.43; 1.0 vs 1; 2.2936 vs 2.33
  (= 1/0.43); 1.1767 vs 1.2. Worst 1.94 %. The "about 250 µm" remark is −8.3 %
  from 229 and is correctly demoted to a loose sentence.
- **Eq. (3) inversion:** min 0.4628 (spent, 115 µm), max 2.9393 (fresh, 25 µm),
  median 0.9809, spread 6.351; Diakon 0.7908–1.6630, and 0.7908–1.1581 above
  100 µm (±21 %). All reproduce. The confound — both loose ends are the two
  powders with approximate ρ_s — is named on the page and the clean Diakon
  subset is reported beside it.
- **ε₀ identity:** 0.4414–0.5705, structural at the two H/H₀ = 1.000 rows and
  labelled structural.
- **Eq. (2) fit residual:** 9.645 % mean over the 19 rows measured group A,
  36.905 % worst at the 115 µm spent-catalyst fraction; labelled a fit residual
  in all four places it appears and never counted as agreement. Note it is
  computed over "rows measured group A", which is not identical to Geldart's own
  Fig. 1 set (he excluded the three coarsest Diakon fractions and included the
  as-received catalyst) — differs by two rows, harmless, and the page states its
  own definition.
- **Deliberate breaks, all re-run:** every entry of break tables 1, 2 and 3
  reproduces. `ρ_f = 0` moves nothing (both constants are in the density
  difference) — real blind spot, correctly named. The eq. (3) spread moves only
  6.351 → 5.688 under a `d¹` exponent while the median moves 0.981 → 0.010 —
  real, and the right lesson is drawn. A +0.10 error in one `ε_MB` is invisible
  to the range test — real. `u_br` from radius rather than diameter moves the
  eq. (8) constant −29 %, so the 0.55 % coefficient comparison is not vacuous.
- **Honest non-claims all present:** no C/A boundary (band PQ is hand-drawn — I
  confirmed on the render that Fig. 3 carries a shaded band and no expression);
  no B/D test against data; nothing digitised; the pymrm cloud number labelled
  discretisation error against a formula E1.2 already validated.
- **pymrm conventions:** no `NumJac` on the page (the solve is linear, `spsolve`
  on a constant operator), so the shape rules do not apply; `construct_grad`
  and `construct_div` are called with `(n, 1)`, never a bare 1-D shape; `nu=2`
  is commented as spherical; no diffusivity jump, so the harmonic-mean rule is
  not engaged. The outward-normal comment next to `bc` is correct, and I
  verified the far-field condition `dh/dr + (ν/r)h = (ν+1)A` is satisfied
  exactly by both radial modes, as the docstring claims.
- **Reproducibility:** `build_page.py` regenerates all 38 cell sources
  byte-identically; a fresh nbclient run reproduces every printed output
  byte-identically in 4.2 s. No warm-start, no continuation chain, no
  non-determinism.
- **Prose against output:** every interpolated number matches; every retyped
  number in the markdown cells (229.4, 1.016e6, 1.9 %, 1.6 %, 20 of 21, 0.5 %
  inside, factor 6.4, 0.46–2.94, median 0.98, ±21 %, 0.55 %) matches the
  executed output. No drift found.
- **Housekeeping:** `review/` is empty apart from this file; no image of any
  kind anywhere under `queue_cases/A1.7/`; no Quarto callouts or shortcodes in
  the notebook; Colab cell first; `load_data(..., page=...)` used;
  `report_agreement` called; section order matches the AGENTS.md contract;
  `slug` and `title` identical between `meta.yaml` and `models_entry.yaml`.

---

## Recommendation

**Safe to publish after the fixes in F1, F2, F3 and F4.** F1 and F2 are the
serious ones: they do not make any printed number wrong, but as the page stands
its central claim is presented with more power than it has, and the tie-break
that produces the margin is invisible. F3 and F4 are sentences that misdescribe
the paper. F5–F8 are worth doing while the file is open but would not on their
own hold the page.
