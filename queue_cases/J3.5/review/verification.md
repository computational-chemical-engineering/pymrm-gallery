# J3.5 verification — adversarial review

Verifier pass on `queue_cases/J3.5/page/`, 2026-07-31. Everything below was
re-derived here; nothing is taken from the builder's account. Paper values read
off 600 dpi `pdftoppm` renders of
`~/papers/pymrm-gallery/Marquis_2019_J._Electrochem._Soc._166_A3693.pdf` and
`~/papers/pymrm-gallery/Doyle_1993_J._Electrochem._Soc._140_1526.pdf`.

**Verdict: safe to publish after the fixes in §1–§5.** No finding shows a wrong
physical conclusion. The defects are (a) one demonstrably overstated claim about
what the validation resolves, (b) a headline number fitted on a point the page
itself flags as unconverged, and (c) three untraceable or stale prose numbers.

---

## What was reproduced

Re-executed the staged notebook without overwriting it (nbclient, same kernel,
34.7 s). **Every stream output is byte-identical to the committed
`index.ipynb`**, and every metric in `agreement.json` matches. No warm-start
continuation anywhere: every `Cell.march` starts from `initial(u0)` and the SPMe
from `c = 0`. Determinism requirement satisfied on this machine.

---

## 1. CONFIRMED — the slope test does not resolve two of the four terms the page says it resolves

The Validation section claims:

> a wrong coefficient anywhere in the correction — ohmic, concentration,
> exchange-current or anode — leaves the SPMe stuck at slope 1

I tested this directly by mutating one correction term at a time and refitting,
against a DFN reference converged to `dt_scale = 0.0625` (RMS in mV, currents
0.25 / 0.5 / 1 / 2 / 4 A m⁻²):

| variant | RMS at I = 1 | slope (3 lowest) | local slopes |
|---|---|---|---|
| baseline | 0.0451 | 2.16 | 2.14 2.18 2.36 2.70 |
| ohmic `δ_c/3` → `δ_c/2.7` (10 %) | 0.444 | **0.97** | 0.99 0.96 0.83 0.42 |
| concentration factor 2 → 2.2 (10 %) | 0.675 | **0.98** | 0.99 0.97 0.89 0.50 |
| exchange-current average removed (`pre_bar` → `PRE0`) | 0.0451 | **2.16** | 2.14 2.18 2.36 2.70 |
| anode overpotential at `c_0` instead of `c(0)` | 0.0494 | **1.89** | 1.94 1.85 1.83 2.55 |

The ohmic and concentration coefficients are resolved superbly — a 10 % error
multiplies the I = 1 error by 10–15× and collapses the slope to 1, exactly as
claimed. But:

* **Deleting the electrode-averaged exchange-current correction entirely
  (their Eqs. 48n–o) changes the error by 0.02 %** — below the discretisation
  noise. That term is unresolvable on this cell, and the page's implementation
  of it is untested by any check on the page.
* **Deleting the anode concentration correction changes the I = 1 error by 9 %
  and leaves the slope at 1.89**, not 1.

Failure scenario: a reader copies `SPMe.voltage` for another cell believing all
four terms were verified here. A sign error in the `pre_bar` averaging would
have passed every check on this page.

**Fix:** restrict the claim to the two terms the sweep actually resolves, and
say that the exchange-current average is below resolution on this cell (with the
0.02 % number, which is itself a useful result — it says the term is negligible
for a 1:1 salt at `c_0 = 1000` against `c_max = 3920`).

## 2. CONFIRMED — the headline slope 2.07 is fitted on a point the page's own check 5c disqualifies

`slope_spme_low` is fitted over `low = [0.25, 0.5, 1.0]`. Check 5c then prints
that quartering the DFN time step changes the apparent SPMe error at I = 0.25 by
66 %, and the narrative says the point "is at the edge of the DFN integrator's
resolution even with the tightened step" — and the fit uses it anyway. That is
diagnosing the contamination and then reporting the contaminated number.

Refining further (`dt_scale` 0.25 → 0.0625, 2312 DFN steps, and SPMe `n_t`
800 → 3200 at the two lowest currents, both converged to ≈1 %):

| I | shipped (dt 0.25, n_t 800) | converged | change |
|---|---|---|---|
| 0.25 | 0.00257 | 0.002415 | −6 % |
| 0.50 | 0.01028 | 0.010117 | −1.6 % |
| 1.00 | 0.04543 | 0.0451 | −0.7 % |

The two discretisation errors act in opposite directions and partly cancel:

* local slopes go from **2.00 / 2.14 / 2.36 / 2.70** (shipped) to
  **2.07 / 2.16 / 2.36 / 2.70** (converged) — the pleasing 2.00 at the bottom of
  the sweep is a discretisation artefact;
* the 3-point fit goes from **2.07 → 2.11**.

Window dependence, on the converged data: `[0.25,0.5,1] → 2.11`,
`[0.5,1,2] → 2.27`, `[1,2,4] → 2.53`. The quoted 2.07 is the *smallest* of the
plausible windows *and* uses the least converged point.

The conclusion "the SPMe error scales as C_e², the SPM as C_e" survives — the
converged low-end local slope is 2.07 and the trend extrapolates cleanly to 2.
But the specific number 2.07 is not a converged measurement, it is 2 % below the
converged value, and the README/`models_entry.yaml` carry it with no qualifier.

Related: **the sweep cannot be extended downward.** The DFN fails at I = 0.125
at every `dt_scale` tried (1.0, 0.25, 0.0625 → 36, 91, 7 steps, `u_max` stuck at
u₀ = 0.010). I = 0.25 is not "the edge of the integrator's resolution", it is the
last current the integrator runs at all. The page should say that.

Failure scenario: CI re-executes on another machine, a BLAS difference changes
one Newton retry, the I = 0.25 DFN trajectory shifts, `slope_spme_low` moves by
more than `check_agreement.py`'s 5 % and the page reports a regression that is
not a regression. (This is the `B1.1` lesson in a new dress: a reported number
that depends on an adaptive integrator's step history.)

**Fix, cheapest first:** either drop I = 0.25 from the fit and report the slope
over `[0.5, 1, 2]`, or run the sweep at `dt_scale = 0.0625` and quote 2.11. Say
that the local slope is what converges to 2, not the window fit, and that the
DFN floor is I = 0.25.

## 3. CONFIRMED — "J3.4: utilisation stops at 0.31" is in no source

Asserted twice (`build_page.py:748` printed output, `:1023` in *What pymrm
adds*) and again in `meta.yaml`. Checked every candidate:

* J3.4's stated-results CSV: `u_at_sharp_drop = 0.30`, sourced to Doyle's text
  ("about 30 % of the cathode material");
* J3.4's own notebook output: the published I = 20 curve is drawn to
  **u = 0.299** digitised; the model's 1.9 V crossing is **0.264**;
* **this page's own I = 20 DFN run: `u_max = 0.4006`** (V = 1.511 at the stop).

0.31 matches none of them, and the page attributes it to J3.4, which does not
contain it. Failure scenario: a reader follows the cross-reference to J3.4 to
check the number and cannot find it, which is exactly the trust cost the brief
is written to avoid.

**Fix:** quote Doyle's own "about 30 %" and/or J3.4's 0.299 digitised, with the
source named; or print this page's own I = 20 `u_max`.

## 4. CONFIRMED — stale count: the prose says 20, the code prints 21

* `build_page.py:303` "Recomputing all **twenty** derived values"
* `build_page.py:865` "**20/20** to printed precision"
* both `data/*.meta.yaml` sidecars: "20 of 20", "20 quantities"

against `agreement.json` `marquis_tables_checked = 21`, `meta.yaml`/`README.md`
"21/21", and the notebook's own printed line *immediately above the second
occurrence*: `worst of 21`. 21 is correct (6 Table II timescales + 15 Table III
derived entries; γ_n = 1 is definitional and correctly skipped). The page
contradicts its own output in the same scroll. Same defect class as the A4.2
catch.

## 5. CONFIRMED — "slope 1.0 … 2.0 … over the mid sweep" (build_page.py:1014)

Two errors in one clause. The fit is 2.07, and it is over the three **lowest**
currents, not the mid sweep; the mid-sweep local slopes are 2.14 and 2.36. Also
`build_page.py:1091` (*Reuse*) says the linearised correction should not be
reused "above I ≈ 8 A/m²" — the page computes 10.2 A/m² and 8 appears nowhere.

---

## Lower-severity findings

6. **PLAUSIBLE — the page never states that a shared error cancels.** The
   positive form is there ("same reviewed parameter set … so the measured error
   is the *reduction's* error and nothing else's"), which is the correct and
   honest statement of what the slope measures. But the converse is never
   written: because `U_ocp`, `k_2`, `σ`, `i_0_1` and above all `KAPPA` enter the
   DFN and the reductions identically, **an error in any of them is invisible to
   every comparison on this page**. `KAPPA` matters most — it is J3.4's
   *reconstruction*, with two published routes disagreeing by 19 %, and it
   appears in the DFN's `kf` and in the SPMe's `-I(δ_s/κ + δ_c/3κε^{3/2})` where
   it cancels to first order. One sentence would close this. It is not
   circularity in the harmful sense — nothing is fitted to anything — it is
   unstated blindness.

7. **CONFIRMED — latent silent failure in the error sweep.** `Cell.march`
   returns whatever trajectory it has when the step collapses, and every caller
   then does `np.interp(UGRID, tr[:,0], tr[:,1])`, which *clamps* past the end
   of the data. At I = 0.125 this produces an SPMe "RMS error" of **401 mV**
   with no warning of any kind. All five shipped currents reach u = 0.870, so no
   shipped number is affected — but there is no guard. One assertion
   (`tr[:,0].max() >= UGRID.max()`) removes the whole failure mode.

8. **CONFIRMED — Table VI is evaluated four-fifths of the way.** The page says
   "Evaluating Table VI for the Doyle cell". The printed Table VI has five
   conditions; the page evaluates C_e ≪ 1, σ̄_k ≫ 1, κ̂_e ≫ 1 and the solid
   diffusion condition (via J3.4's S_c), and never evaluates the fifth
   (reaction timescale, `I/(m a c_e^{1/2} c_n,max L) ≪ 1/C_e`).

9. **CONFIRMED — "their rounding" is the wrong word.** The paper *truncates*:
   4.1952e-3 → prints 4.19e-3 (3 s.f. rounding would give 4.20e-3); 0.044251 →
   0.0442; 0.11346 → 0.1134; 0.44444 → 0.4444. The pattern is consistent across
   all five deviations above 5e-4, which actually strengthens the check. Say
   "truncation to four digits", or a careful reader will read 4.20 ≠ 4.19 as a
   real disagreement.

10. **INFO — the SPM's "slope 1.00" carries almost no information.** On this
    cell κ and D are constants and the first-order electrolyte perturbation is
    linear in I, so *every* term the SPM omits is exactly proportional to I. The
    measured SPM errors (3.0393, 6.0742, 12.1387, 24.2999, 49.236) are linear to
    4 digits because they cannot be anything else. Only the SPMe slope is a
    test. Worth one clause so the two headline slopes are not read as equally
    strong evidence.

11. **PLAUSIBLE — SPMe time-step independence is demonstrated only at I = 1.**
    Check 5a gives 0.0002 mV for n_t 800 → 3200 at I = 1, which is 0.4 % of the
    signal there but **6 % of the signal at I = 0.25** — the point that carries
    the headline fit. The check should be run at the sweep's low end, not its
    middle.

---

## What I attacked and could not break

**Both claimed defects in the published paper are real.** Verified independently
on 600 dpi renders, not from the text layer:

* *p. A3701*: "As provided in **Table III**, we have 𝒞ₑ = 5.1 × 10⁻³𝒞 where 𝒞 is
  the C-rate." Table III (p. A3696) prints **4.19 × 10⁻³𝒞**. The sentence names
  the very table it contradicts, so it cannot be a different quantity, cell or
  definition. Recomputation from Table I gives 4.1952e-3𝒞 (τ_e = 94.803 s,
  τ_d = 22598 s), confirming the table. The claim stands.
* *Table II, p. A3695*: the τ*ₙ row prints `(R*ₙ)²/D*_{s,p}` — the subscript is
  crisply `s,p`, identical to the τ*ₚ row beneath it, in a born-digital PDF, so
  it is not a scan artefact. Its printed value 2.5641 × 10³ s = (10 µm)²/
  3.9e-14 = `D*_{s,n}`; `D*_{s,p}` would give 1 × 10³ s, which is the τ*ₚ row.
  A subscript typo, proved from the paper's own value. The claim stands.

**The 21/21 recomputation is not circular.** Every input the check consumes is
printed in Table I (verified against the rendered table row by row); nothing is
back-solved from Tables II/III. Table III's own inputs (c_k,0, cutoff) are
carried as inputs, not recomputed. The check genuinely tests the transcription.

**𝒞ₑ ≡ S_s is an exact algebraic identity, not a numerical coincidence.** Doyle
Eq. 27, read at 600 dpi, is

  S_s = (δ_s + δ_c)² · I / (D F (1 − ε) c_T δ_c)

and `Ce_doyle(I) = (δ_s+δ_c)²/D · I / [(1−ε) c_T δ_c F]` is the same expression
term for term. It holds at any parameter values. One nuance the page could state:
Marquis's *printed* form is 𝒞ₑ = I L/(D F c_n,max) (Table VI), which uses a
nominal capacity c_n,max·L; the equality with S_s follows once τ_d is built from
the *actual* cathode capacity (1−ε)c_T δ_c, which is the natural half-cell
choice and is Doyle's. So: same ratio, same physics, but the identity is with
the page's half-cell form of 𝒞ₑ, not with Marquis's literal formula.

**The applicability finding is computed from the paper's own criterion.** Table
VI's conditions were read on the render and match `Ce_doyle` and
`kappa_hat_doyle` exactly. Doyle's Figure 2 carries I = 5, 10, 13 and 20 A m⁻²
(plus a −10 charge), giving (𝒞ₑ, κ̂ₑ) = (0.077, 2.1), (0.153, 1.04),
(0.199, 0.80), (0.306, 0.52). "Mostly outside the regime" is supported. The
negative-concentration threshold at 10.2 A m⁻² is correctly attributed to the
page's own first-order profile, not to Table VI. *Minor:* the page prints κ̂ₑ
only at I = 1, 10 and 20, so the reader cannot check the I = 5 and I = 13 cases
that the "mostly" rests on.

**Conventions and physics.** Boundary conditions use the outward normal with
homogeneous Neumann at both ends, and the physical zero-*total*-flux condition is
imposed by zeroing the migration flux on the boundary faces — so the discrete
boundary flux is exactly the physical zero and salt conservation is structural
(4.5e-13 confirmed on re-execution). The wall extrapolation
`c(0) = c[0] + x_c[0](1−t⁺)I/(FD)` has the correct sign (salt accumulates at the
foil on discharge, matching Doyle's Fig. 4 text) and is written identically in
`Cell.cell_potential` and `SPMe._record`. `bv_invert` reduces exactly to
`U_ocp` as I → 0 (checked algebraically: z → √((c_T−c_s)/c_s), α = 0.5). The
deviation convention (reduced − full) is stated once and used everywhere. `nu=0`
Cartesian is stated in a comment; operators are assembled once.

**Honesty of the tier-6 framing.** "Nothing on this page touches a measurement",
"tier 6 twice over", "Do not describe anything here as validated against
experiment" — stated in the Background, in `README.md`, in `meta.yaml`'s caveats
and in `models_entry.yaml`. Nothing anywhere calls this experimental validation.
The constant-κ/D caveat is present, correctly says the SPMe's absolute accuracy
is flattered here, and correctly says the exponents rather than the mV values
are the transferable result. (Finding 10 above is the one place that caveat
should reach further: with constant κ the SPM's slope 1 is structural.)

**Housekeeping.** No Quarto-only markdown in the notebook; nine required
sections in order; `report_agreement` called; no new dependencies; `slug` and
`title` identical between `meta.yaml` and `models_entry.yaml`; `review/`
contains no page images (it was empty on arrival, contrary to the case notes,
which is the safe direction).
