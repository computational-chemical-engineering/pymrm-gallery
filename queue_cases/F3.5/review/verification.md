# F3.5 verification — Bosch, Versteeg & van Swaaij (1989)

Verifier pass, 2026-07-31. Everything below was re-derived from the PDF at
600 dpi and from re-executions of the staged notebook; nothing is taken from
the builder's account.

**Verdict: the reconstruction is LEGITIMATE BUT UNDER-CAVEATED, and one of its
supporting checks points the other way.** The scalar itself survives every
attack I could mount — it is independently corroborated three times over by
figure data that contain no enhancement factor at all. But the page's headline
metric is flattered by the definition of *E*, its statement that the scalar was
"never" fitted to the enhancement factors is false as written, and the
desorption comparison it presents as corroboration in fact *prefers the printed
constants*. Two figure readings on the page are wrong and one contradicts the
page's loading convention.

**Recommendation: safe after the fixes listed in §9. Do not publish as-is.**

---

## 1. What survived the attack (verified, no finding)

- **Table 1 transcription: exact.** All 15 rows re-read on a 600 dpi render of
  page 2738. Every mantissa and every exponent matches the CSV, including the
  three the OCR mangles (`K_c2` 7.29e-8, `K_p` 3.17e-5 / 3.17e-7, `K_c`
  7.06e-4 / 1.41e-4 — the text layer renders the last as `1.41 X lo-*`).
- **Table 2 transcription: exact.** All 32 numbers re-read on page 2739. Note
  the OCR gives `F_HDA = 3.74` in row 3; the printed value is **3.76**, as the
  CSV has it. The builder read the image, not the text layer.
- **The F = J-ratio identity re-run independently:** max deviation 0.18 %,
  matching `table2_F_vs_Jratio_maxdev`.
- **Determinism:** `agreement.json` is **byte-identical** (md5 2e86ceab…) after
  a clean re-execution. No warm-start-dependent quantity appears in it. Runtime
  27 s.
- **Prose vs. code:** all computed numbers quoted in the markdown match the
  executed outputs *except* the two in §4 and §5 below.
- **pymrm conventions:** the source sign matches `A = div(-grad)` (a consumption
  term enters positive — verified against the residual assembly); `nu=0` is
  correct for the slab; `NumJac(shape)` last-axis stencil is right for a
  pointwise source; `J = -D dc/dx|_0` and `E = J/(k_L(c_i-c_b))` are extracted
  with the right signs. The outward-normal trap is *vacuous here* — both ends
  are Dirichlet or homogeneous Neumann, so no sign can be wrong. The page's
  claim to have written the BCs on the outward normal is true but tests nothing.
- **Alternatives correctly rejected (I re-ran them):**
  - `K_c1 / 1.4585`: fixes c_b (0.628 vs implied 0.638) but leaves `E_u` **+19.8 %**
    high. The page's argument holds exactly as stated.
  - scaling `c_carb`: algebraically identical to `K_c1` (OH⁻ = K₃·CO₃/HCO₃ is
    invariant), so it fails the same way.
  - scaling `k_OH` alone: fixes E, cannot touch c_b.
  - **a constant offset in loading** (a hypothesis the page does *not* consider,
    and a live one given the paper's admission that the concentrations "had to
    be guessed"): the best single Δα is +0.052, but it fits the four implied
    c_b only to **1.185 / 0.943 / 0.984 / 0.910** (a ±19 % spread, against
    ±2.5 % for the K₃ scalar), and leaves `E_u` up to 9.2 % and F up to 15.3 %
    off. Decisively rejected by the paper's own numbers.

So among single-constant changes, K₃ = K_w/K_c2 really is selected. That part
of the F2.3 standard is met.

---

## 2. FINDING 1 (HIGH, CONFIRMED) — "E_u within 0.6 %" is a flattered number, and the circularity disclaimer is false as written

The page states:

> the scalar's *magnitude* comes from the four unpromoted implied c_b values
> only — **never from the enhancement factors**

The implied c_b **is computed from the enhancement factors**:
`c_b = c_i − J/(E·k_L)`, printed two cells earlier on the same page. The scalar
is fitted to a specific combination of J and E. The sentence as written is not
true, and it is the single most load-bearing honesty claim in the section.

Worse, the quantity then reported as the validation is *less* sensitive to the
fit residual than the flux is. From the executed run (reconstructed pass):

| condition | `E_u` dev | `J_u` dev | implied s* |
|---|---|---|---|
| α 0.2, 30 kPa | +0.2 % | −0.1 % | 1.420 |
| α 0.4, 30 kPa | +0.6 % | **+3.8 %** | 1.481 |
| α 0.4, 120 kPa | +0.6 % | +0.1 % | 1.419 |
| α 0.6, 120 kPa | +0.6 % | **+3.4 %** | 1.490 |

The four per-condition scalars the identity actually demands span 1.419–1.490
(±2.5 % about the mean). That 5 % spread lands almost entirely in `J_u`, and is
suppressed in `E_u` because the fitted c_b sits in the denominator of E and
partially cancels the numerator error. The four `E_u` deviations are
+0.2/+0.6/+0.6/+0.6 % — essentially **one constant offset repeated four times**,
not four independent agreements.

I confirmed the mechanism directly: `E ∝ s^−0.49` over the fit region
(s = 1.40 → mean `E_u` dev +2.50 %; s = 1.50 → −0.86 %). Refitting the scalar to
`E_u` **alone** gives **s = 1.4735**, 1.0 % from the c_b-fitted 1.4585 — and
half of that 1.0 % is exactly the residual +0.5 % mean `E_u` offset. So "0.6 %"
is, to first order, a restatement of "the two fitting routes agree to 1 %".

**Failure scenario.** A reader takes "unpromoted E reproduced to 0.6 %" as
evidence the film model is right to sub-percent accuracy, and reuses the `Film`
class believing it is validated at that level. The defensible number is
**3.8 % on the flux**.

*This is not fatal circularity* — 8 unpromoted numbers (4 J, 4 E) are matched
with 1 parameter, and the two fitting routes agreeing to 1 % is real content.
But the page must say so in those terms.

---

## 3. FINDING 2 (HIGH, CONFIRMED) — the desorption checks are not independent of the scalar, and they favour the printed constants

The page presents these as out-of-sample corroboration:

> the stated desorption promotion factor 3.81 to 3.35 (−12 %) … the
> driving-force reduction ("about one third": 39 % here)

Both are computed with `k3_scale=s_fit`. I re-ran them with the printed
constants:

| | F_desorption (paper: 3.81) | driving-force reduction (paper: "about one third") |
|---|---|---|
| printed constants | **3.721 (−2.3 %)** | **34.7 %** |
| reconstructed (s = 1.4585) | 3.355 (−12.0 %) | 39.1 % |

The **only two numbers in the paper that are not in Table 2** — the ones a
reconstruction fitted to Table 2 cannot have absorbed — are reproduced *better,
on both counts, by the constants exactly as printed*. The reconstruction makes
the desorption promotion factor five times worse and pushes "about one third"
to 39 %.

The page reports the −12 % as a residual "unexplained" deviation and attributes
it to the paper's HDA uncertainty. That framing is not available once you know
the printed constants give −2.3 %: the deviation is *caused by the
reconstruction*, and the page must say so.

**Failure scenario.** A reader concludes the scalar is corroborated
out-of-sample, when the only genuinely out-of-sample statement in the paper
argues against it. This is exactly the "plausible, confident agreement" the
verifier brief exists to catch.

*Caveat in the reconstruction's favour:* Fig. 4's own bulk CO₂ (§4) strongly
supports the scalar for the same α = 0.6 desorption case, so the paper is
internally inconsistent here, not the page. But the page must surface the
tension rather than report only the half that fits.

---

## 4. FINDING 3 (HIGH, CONFIRMED) — Figures 2 and 3 are the α = 0.2 case, they contradict the page's loading convention, and they are the strongest evidence the page never uses

I calibrated all three figures from their tick marks and verified every axis
assignment by electroneutrality and by the amine balance.

**Figure 4 (α = 0.6, 30 Pa)** — unpromoted HCO₃⁻ = 2382, CO₃²⁻ = 802
(charge 3986 ≈ [K⁺] = 4000 ✓). This *is* the page's convention (2αc_carb,
(1−α)c_carb).

**Figures 2 and 3 (both labelled α = 0.4, 30 kPa)** — unpromoted
**CO₃²⁻ = 1597–1600, HCO₃⁻ = 806–810** (charge 4010 ≈ 4000 ✓; the labels are
unambiguous at 600 dpi and the reverse assignment is not electroneutral). That
is the composition the page's convention assigns to **α = 0.2**, not α = 0.4.
The page's model at α = 0.2 gives HCO₃⁻ = 806, CO₃²⁻ = 1594.

The paper is internally inconsistent about α between Figs 2/3 and Fig 4/Table 2.
Table 2 follows Fig. 4 (the implied-c_b consistency test collapses under the
other convention — the required scalar becomes 6.6–7.9 with a 20 % spread).
The most economical reading is that the "α = 0.4" in the Fig. 2/3 insets is a
misprint for α = 0.2.

**And that makes Figs 2/3 the single best evidence in the paper for the
reconstruction — evidence that touches no enhancement factor at all.** Read
against the page's own α = 0.2 solve:

| quantity | Fig. 2 / Fig. 3 (measured here) | printed constants | reconstructed (s = 1.4585) |
|---|---|---|---|
| bulk CO₂ (unpromoted) | **0.119 / 0.117** | 0.0826 | **0.1194** |
| bulk CO₃²⁻ | 1600 / 1597 | 1591 | 1594 |
| bulk HCO₃⁻ | 810 / 806 | 809 | 806 |
| HDA at interface → bulk (Fig. 2) | **372 → 484** | — | **374 → 484** |

The unpromoted bulk CO₂ alone implies **s = 1.43**, from the figure's own
speciation and the printed K_c1, K_w, K_c2 — no Table 2, no E, no J. Fig. 4
independently implies **s = 1.50** (see Finding 4). Both bracket the fitted
1.4585.

The HDA profile match (372→484 vs 374→484, <1 %) also corroborates the
*promoted* chemistry, which nothing else on the page does independently.

**Consequences for the page as staged:**

1. The claim "the loading convention is confirmed against the paper's own
   Figure 4" is selective: two of the paper's three profile figures contradict
   it. The convention is right *for Table 2*, and the page must say that and
   flag the figure inconsistency.
2. The profile plot is titled **"the paper's Fig. 2 conditions"** and is run at
   α = 0.4. By the paper's own figure content it is not Fig. 2's condition.
   Either re-run the profile figure at α = 0.2 (where it reproduces Fig. 2
   quantitatively) or drop the claim of correspondence.
3. The hard-coded prose under that plot is wrong for the plotted case: I get an
   interface CO₂ gradient ratio of **7.55×**, not the stated "~5×" (the ~5×
   figure is the α = 0.2 value, 4.77×). HDA depletion is 10.3 %, stated as
   "~12 %" — close but also hard-coded rather than computed.

---

## 5. FINDING 4 (MEDIUM, CONFIRMED) — two stated numbers do not match their sources

**(a) Figure 4's bulk CO₂.** The page says:

> its dashed (unpromoted) CO2 profile reaches a bulk level of ≈ 2.3–2.4 mol/m³
> at α = 0.6, where the printed constants give 1.44 and the implied value is 2.15

Calibrated left axis (ticks at 2.0/1.0/0.0, y = 703.5/1433.5/2168.5), the dashed
CO₂ curve reaches **2.16** at x/δ = 1, not 2.3–2.4. The direction of the argument
is unaffected — and the corroboration is in fact *better* than claimed (2.16 vs
implied 2.152, 0.4 %; vs printed 1.445, 50 % off). But a number attributed to a
source figure must be the number the figure shows.

**(b) The printed-constants offset.** The markdown says:

> Every enhancement factor comes out **15–20 %** high — all twelve cases

The notebook's own output two cells earlier gives **+8.4 % to +21.4 %**
(`E_HDA` at α = 0.2 is +8.4 %, `E_DEA` +9.7 %). `meta.yaml` and `README.md`
correctly say 8–21 %. This is a stale hard-coded range contradicted by the
executed cell directly above it — the third instance of this class of defect
this session.

---

## 6. FINDING 5 (MEDIUM, CONFIRMED) — the location choice is three-way, not two-way

The page frames the location as a discrete choice between K₃ and K_c1. There is
a third alternative it does not name: **K_w**. Because K₃ = K_w/K_c2, scaling
K_w down by 1.4585 is *indistinguishable from scaling K_c2 up* for every
unpromoted quantity, but it additionally scales K_B = K_w/K_p and
K_eq = K_c K_c1 K_w/K_p, which only the promoted cases see. The code implements
the K_c2 branch (K_B and K_eq untouched), while the prose, `meta.yaml`, README
and `models_entry.yaml` all say "K_w/K_c2", which is ambiguous between them.

I ran both:

| variant | `E_u` max | promoted E max / rms | F max / rms |
|---|---|---|---|
| K_c2 (as coded) | 0.63 % | 7.58 % / 4.94 % | **5.63 % / 3.27 %** |
| K_w | 0.63 % | 7.52 % / 4.81 % | 8.60 % / 4.40 % |

The unpromoted data cannot separate them at all; the promoted E's cannot
either; only the promotion factors mildly prefer the K_c2 branch, at the level
of the residual scatter. The page should state that the *magnitude* is well
determined, the *K₃-versus-K_c1* choice is decided by the E offset, and the
*K_c2-versus-K_w* split inside K₃ is essentially undetermined by these data.

---

## 7. FINDING 6 (MEDIUM, mixed) — the ionic-strength story is attached, not derived

The page and `meta.yaml` call 1.458 "the size of an ionic-strength correction at
the ~5 M ionic strength of this solution".

- **The paper says the opposite.** Section 3: *"Non-idealities in the liquid
  phase equilibria were not taken into account."* (Compatible with their having
  used an apparent constant, but the page cites ionic strength without
  acknowledging the sentence.) — CONFIRMED
- **The magnitude is not derived.** s = 1.458 is 0.164 pK units on K_c2
  (pK_a2 10.14 → 9.97). Apparent pK_a2 shifts at I ≈ 5–6 M in K⁺ media are
  typically several times larger. Right sign, unverified magnitude — "the size
  of" claims a quantitative match that is not demonstrated. — CONFIRMED
- **A competing story fits at least as well.** Putting the factor on K_w
  instead requires K_w,eff = 2.28e-7 mol²/m⁻⁶, i.e. pK_w = 12.64 against the
  printed 12.48 — about a 10 K difference in the reference temperature for the
  water ion product on the standard pK_w(T) curve. That is as plausible a
  clerical explanation as the ionic-strength one and points at a different
  constant. — PLAUSIBLE

The honest formulation is: the paper's own results fix the *magnitude* of the
correction and localise it to K₃; *why* K₃ differs from its printed parts is
not determined by anything on this page.

---

## 8. FINDING 7 (MEDIUM, CONFIRMED) — what the validation battery can and cannot catch

I asked of each check what error class it can detect. Three of the eight are
weaker than the page implies.

- **Electroneutrality (2.8e-11) and charge flux (1.1e-9) are algebraically
  inevitable.** Every charged species carries the *same* `D_ion` (Am and CO₂,
  the two with different D, are neutral), every reaction conserves charge, and
  both boundaries are Dirichlet at electroneutral states. The charge-weighted
  combination therefore satisfies the same linear operator with zero source and
  zero boundary data — it is identically zero for *any* converged solution,
  correct or not. This is the A4.2 failure mode. The code comment half-admits
  it ("equal ion diffusivities keep it so"); `meta.yaml` lists it as validation.
  The carbon-flux and amine-flux closures are *not* inevitable (mixed
  diffusivities) and are genuine, if weak.
- **"Independent dimensionless second-order assembly", 1.0e-8.** The two solves
  share the grid (`stretched_faces(400, s=6)` in both), the same
  `construct_grad`/`construct_div`/`NumJac`/`newton` calls, and the same
  discretisation. What is genuinely written twice is the *non-dimensionalisation
  and bookkeeping* (δ, D scaling, the ν = 2 stoichiometry, E_inf). So 1e-8
  proves the dimensional bookkeeping, and cannot detect any discretisation
  error common to both. The README's "an independent dimensionless second-order
  assembly to 1e-8" overstates it. (It is *not* the pure A4.2 case — a real
  error class is covered — but the label must be narrowed.)
- **Reversible → irreversible collapse (9.7e-6)** uses the same `Film` class
  with `keq_scale` vs `irreversible2`. It checks that the reverse term is
  written so it vanishes as K_eq → ∞. It is not an independent identity, and
  calling it "the internal identity the reversible model owes the gallery's own
  F3.1 page" is decoration.
- **Working as advertised:** the pseudo-first-order sweep against the closed
  form (3.5e-5 over 0.1 ≤ Ha ≤ 300) is the one check that really tests the
  discretisation on a boundary layer, and grid doubling (1.3e-5) backs it. The
  physical-absorption limit (5.4e-11) is exact for any conservative scheme on a
  linear profile, so it tests unit bookkeeping, not accuracy. VKH at 0.1 % is a
  genuine external reference. Higbie at 10.5 % is a real model-form sensitivity,
  though it rests on dividing out a +1.42 % quadrature bias calibrated on the
  physical-absorption case (disclosed in the output, not in `meta.yaml`).

---

## 9. Lower-severity items

- `build_page.py:254` — "tests 24 of the 32 numbers". The F = J_p/J_u identity
  tests 4 J_u + 4 J_DEA + 4 J_HDA + 4 F_DEA + 4 F_HDA = **20 of 32**. The 12
  numbers it does *not* test are exactly the twelve **E values**, i.e. precisely
  the numbers the reconstruction and the headline both lean on. Worth saying
  explicitly rather than mis-stating the count. (CONFIRMED)
- `bosch1989-table2-results.meta.yaml` — "all eight agree to print precision
  (max 0.6 %)". It is **0.18 %**, as `agreement.json` records. (CONFIRMED)
- Gas-side resistance: the page writes the fraction as `E k_L m / k_g`, which
  evaluates to 2.4e-7, not the quoted bound. The correct grouping is
  `E k_L / (m k_g)` = 2.4e-5 < 3e-5. The bound is right; the formula printed
  next to it is not. (CONFIRMED)
- `vkh_dev` is computed but omitted from `metrics`, so `check_agreement.py`
  will not track a number `meta.yaml` advertises. (CONFIRMED)
- Table 3's experimental entries are printed as "±4" and "±6"; the page renders
  them both ways ("±4"/"about 4"). Fine, but pick one and note the glyph is
  almost certainly "≈".
- The page does **not** claim experimental validation anywhere, states tier 6 in
  four places, and the deviation convention (model − paper)/paper is defined
  once and used consistently in every comparison I checked. No finding.

---

## 10. Required fixes before publication

1. **Rewrite the circularity paragraph.** Delete "never from the enhancement
   factors" — the implied c_b is built from E. State: the scalar is fitted to
   the (J, E) pairs through the definition of E; the four per-condition scalars
   it must reconcile span 1.419–1.490; refitting to `E_u` alone gives 1.4735,
   1.0 % away.
2. **Demote the headline.** Report **`J_u` within 3.8 %** as the reproduction
   metric alongside `E_u` 0.6 %, and say that `E_u` is the less sensitive of the
   two because the fitted c_b sits in its denominator. Keep both in
   `agreement.json` (`J_u` max dev is not currently reported at all).
3. **Report the adverse desorption result.** With the constants exactly as
   printed, F_desorption = 3.72 (−2.3 %) and the driving-force reduction is
   34.7 %; the reconstruction gives 3.35 (−12 %) and 39 %. Present this as a
   tension in the paper, not as corroboration.
4. **Fix the Fig. 4 reading** (2.16, not 2.3–2.4) and the printed-constants
   range (**8–21 %**, not 15–20 %).
5. **Add Figs 2 and 3, and fix the convention claim.** They independently imply
   s = 1.43 with no enhancement factor involved, and the α = 0.2 solve
   reproduces Fig. 2's HDA profile (374→484 vs 372→484) and bulk CO₂ (0.1194 vs
   0.119). This is the page's best evidence and it is currently unused. At the
   same time, say plainly that Figs 2/3 label that composition α = 0.4 while
   Table 2 and Fig. 4 label it α = 0.2 — the paper is inconsistent, and the
   page's convention is confirmed *for Table 2 and Fig. 4 only*.
6. **Fix the profile figure.** Either re-run it at α = 0.2 and keep the "paper's
   Fig. 2 conditions" caption, or keep α = 0.4 and drop the caption. Replace the
   hard-coded "~5×" (the plotted case gives 7.55×) with a computed value.
7. **Name K_w as the third location** and state that the unpromoted data cannot
   separate it from K_c2.
8. **Soften the ionic-strength attribution** to "sign consistent with an
   ionic-strength effect; magnitude not derived; a ~10 K error in the K_w
   reference temperature fits equally well", and note the paper's own
   "non-idealities … were not taken into account".
9. **Re-label the weak checks:** electroneutrality/charge flux as structurally
   guaranteed by equal ion diffusivities (keep them as regression guards, not
   evidence); the dimensionless cross-check as a test of the dimensional
   bookkeeping on a shared grid; the K_eq collapse as a self-consistency check.
10. Small: 20-of-32 not 24-of-32 (and say the 12 untested are the E's); sidecar
    0.18 % not 0.6 %; gas-resistance formula; add `vkh_dev` to `metrics`.

None of these requires re-deriving anything: the model, the transcriptions, the
scalar and the numerics all stand. The work is in the prose, two figure
readings, one plot condition, and three metrics.
