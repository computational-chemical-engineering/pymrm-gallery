# H1.4 — verifier report

Adversarial verification of the staged page `queue_cases/H1.4/page/` against
Itoh, *A membrane reactor using palladium*, AIChE J **33**(9) 1576–1578 (1987).

Everything below was re-read from the paper independently of the builder's
crops and `NOTES.md`, and independently of the PDF text layer. **The scan is a
300 dpi bilevel CCITT image** (`pdfimages -list`), so a 600 dpi render is a 2×
upsample and carries no extra information; the marginal glyphs were finally
read off the *native* 300 dpi bitmap (`pdfimages -png`), printed as ASCII maps.

**Verdict: safe to publish after the fixes listed in F1–F4.** The two
substantive claims of error in the paper (Eq. 5, and the K_P expression) are
upheld and are, if anything, understated. The third "misprint" is not one. The
headline agreement is real but has almost no resolving power, and the page does
not say so plainly enough.

---

## Verdicts on the three misprint claims

| Claim | Verdict |
|---|---|
| (a) Eq. 5 is sign-flipped | **UPHELD** — and provable more strongly than the page does |
| (b) `v_A0` exponent is 10⁻⁵ | **Modelling choice UPHELD and robustly proved; the "misprint" framing is OVERTURNED** — the printed text already says 10⁻⁵ and reads unambiguously |
| (c) the printed `K_P` expression cannot be right | **UPHELD** — not explicable by units or by a sign flip |

### (a) Eq. 5 — upheld

Read at native resolution, Eq. 5 prints exactly `v_H = u_H − 3(u_C⁰ − u_C)`,
and Eq. 4 prints `u_B = u_C⁰ − u_C`.

The page argues the printed form is negative. It is worse than that: the
paper's own hydrogen balance (integrate Eq. 2 and Eq. 3 with `dv_H/dL = Q_H`)
gives `u_H + v_H = 3(u_C⁰ − u_C)`. Substituting that into the printed Eq. 5
yields `v_H = −v_H`, i.e. **the printed equation forces `v_H ≡ 0`**, which
contradicts Eq. 3's permeation term outright. That is a cleaner argument than
"it would be negative" and is worth using.

Sign-convention escape routes were checked and closed:

* the flow direction is not a convention mismatch. Eq. 8 gives *all* initial
  conditions at `L = 0` (`v_H = 0`, `v_A = v_A⁰` there), Figure 3 draws `v_A⁰`
  entering the separation side at the left alongside `u_C⁰`/`u_A⁰`, and
  Figure 1 shows the purge entering at the feed end and leaving at the product
  end. The reactor is co-current in the apparatus, in the flow model and in the
  initial conditions;
* the Notation list defines `v_i` as "flow rate of gas *i* in separation-side
  stream", so `v_H ≥ 0` in the `+L` direction.

The page does not use Eq. 5 to close the model (it integrates `dv_H/dL = +Q_H`
and uses the corrected form only as a check), so nothing computed depends on
this call.

### (b) `v_A⁰ = 11.8 × 10⁻⁵` — right number, wrong framing

**The text is not ambiguous.** On p. 1578 the superscript reads `-5` cleanly at
native resolution, in the same line as an unambiguous `-7`; a 5 (flat top bar +
open bowl) is not confusable with the 6 in the Figure 4 legend (top hook +
closed loop), which I extracted for comparison. The page, `README.md` and the
data sidecar all describe this glyph as "marginal"/"ambiguous". It is not.

So there is **no misprint in the text**, and the page asserts an error the
paper does not contain. The only genuinely ambiguous glyph is the
*hand-lettered* abscissa exponent of Figure 4, whose native bitmap is a blotted
blob with a small closed counter — I could not resolve it either.

That figure ambiguity is settled without any glyph reading and without any
model: the text's headline run (`u_C⁰ = 2.90e-7` ⇒ the `0.29` curve;
`X = 0.997`) must be the open circle at abscissa ≈ 11.7 sitting at `X ≈ 1.0`,
so abscissa 11.8 ↔ `v_A⁰ = 11.8e-5` and the scale is × 10⁵. Under the ×10⁶
reading the `0.29` curve could not reach 1.0 at all: its ceiling at
`v_A⁰ = 15e-6` is 0.783.

**The ceiling argument itself survives every attack I could mount.**

* The ceiling formula is correct. With α_H → ∞ the two sides carry equal H₂
  mole fraction at every `L`, which forces `u_H/v_H = F₀/v_A⁰` pointwise and
  therefore also at the outlet, so the outlet split used is exact; with
  instantaneous kinetics the outlet is at equilibrium. Both conditions are
  necessary upper bounds, so `X ≤ ceiling`.
* It depends on the co-current arrangement, which the paper states three
  separate ways (above). Under counter-current the same purge would allow
  `X = 0.93`, so this dependency is load-bearing — but it is pinned.
* It depends on `K_p`, which is a reconstruction. It is robust anyway: the
  ceiling at `11.8e-6` is 0.716 with the reconstructed `K_p` and 0.664 with the
  independent van 't Hoff value. To reach the measured 0.997 at `11.8e-6` would
  need `K_p = 8.0e13 Pa³`, **338× the reconstruction**, implying a
  membrane-free equilibrium conversion of 0.734 against the paper's 0.187. The
  argument therefore needs only that standard thermodynamics is not wrong by
  two and a half orders of magnitude.

**The builder's unquoted supporting claim is true — and quantitatively so.** I
digitised Figure 4's three *calculated* curves myself from the native scan
(frame calibrated on the y-axis 0/0.5 ticks and the x-axis 0/5/10/15 ticks;
23 points) and ran the page's model against them:

| `K_p(473)` | mean dev | rms | max abs |
|---|---|---|---|
| 2.357e11 (the reconstruction) | +0.0057 | **0.0105** | 0.0216 |
| 1.49e11 (van 't Hoff, ÷1.58) | −0.0354 | 0.0370 | 0.0620 |
| 3.73e11 (×1.58) | +0.0447 | 0.0484 | 0.0705 |

So Itoh's own calculated curves reproduce to ~0.01 in conversion **and they
discriminate `K_p` at the factor-1.6 level, selecting the reconstruction over
the van 't Hoff estimate.** This is the strongest independent support the
reconstruction has, and the page does not use it. (It is a reproduction of the
author's calculation, not a validation — but it is evidence about `K_p` that
did not enter the reconstruction.)

### (c) the `K_P` expression — upheld

Read at native resolution: `K_P = 4.89 × 10³⁵ exp (3,190/T)  (Pa³)`, with no
minus sign inside the exponential, on the same line group as
`k = 0.221 exp(−4,270/T)` — whose minus sign is plainly present — and
`K_B = 2.03 × 10⁻¹⁰ exp(6,270/T)`.

* At 473 K that is `4.15e38 Pa³`, and the corresponding equilibrium conversion
  is 1.000, against the paper's own printed 18.7 % and Figure 4's printed
  "Equilibrium conversion (0.187)" intercept (which I confirmed sits at
  X ≈ 0.185 on the traced frame).
* **A units error cannot explain it**, in either direction. Dimensional
  analysis of Eq. 7 forces `[K_p] = Pa³` (`K_p p_C/p_H³` is compared to `p_B`,
  and `K_B K_p p_C/p_H³` must be dimensionless, giving the printed `[K_B] =
  Pa⁻¹`). Re-reading the printed value as atm³ divides it by 1.04e15 and still
  leaves 4e23 atm³ against the 2.3e-4 atm³ required — 27 orders out.
* **The sign of the exponent is independently wrong.** Cyclohexane
  dehydrogenation is endothermic (ΔH° ≈ +206 kJ/mol), so `K_p` must rise with
  temperature, i.e. the argument must be negative. The printed `+3,190/T`
  implies an exothermic reaction with ΔH ≈ −27 kJ/mol.
* The reconstruction uses only printed inputs and I reproduced it by hand:
  `K_p = 27 x⁴ P₀³ / ((1−x)(1 + n_Ar + 3x)³)` with x = 0.187, y_C⁰ = 0.197,
  P₀ = 1.013e5 gives 2.3563e11 Pa³ (page: 2.3566e11). Δn = +3 is handled
  correctly.

*Optional strengthening (not a defect).* The decision table lists only three
candidates and omits the most likely typographic failure of all: a corrupted
**exponent magnitude**. The printed prefactor is thermodynamically sensible —
`exp(ΔS°/R) ≈ 9e33 Pa³` against the printed 4.89e35 — so the natural diagnosis
is that the prefactor survived and the argument did not. The slope that makes
the printed prefactor reproduce the paper's own x_eq is −26,480 K, within 7 %
of the standard van 't Hoff slope −24,810 K. Saying this converts the
reconstruction from "we replaced the whole expression" into "the corruption is
confined to one number, and repairing only that number reproduces both the
paper's x_eq and standard thermodynamics".

---

## Findings, ranked

### F1 — MAJOR, CONFIRMED. The headline agreement has almost no resolving power, and the page's framing implies that it does

`agreement.json` already records `stated_run_fv_vs_ceiling_absdiff = 0.0`: the
reported model conversion **is** the algebraic fast-permeation/equilibrium
ceiling, exactly. It is therefore a closed-form function of
`(K_p, u_C⁰, v_A⁰, y_C⁰, P₀)` alone. Perturbation test on the page's own model:

| variant | X_model | dev vs 0.997 |
|---|---|---|
| baseline | 0.99834 | +0.13 % |
| α_H × 0.3 … × 10 | 0.99834 (unchanged) | +0.13 % |
| k × 0.1 … × 10 | 0.99834 (unchanged) | +0.13 % |
| K_p ÷ 1.58 (the van 't Hoff value) | 0.99739 | **+0.04 %** |
| K_p × 1.58 | 0.99895 | +0.20 % |
| K_p ÷ 3 | 0.99510 | −0.19 % |
| K_p × 10 | 0.99983 | +0.28 % |
| K_p ÷ 10 | 0.98435 | −1.27 % |

So: the permeation constant may be wrong by a factor 3 either way, the kinetic
constants by an order of magnitude either way, `V_r` and the tube-radius
reading by anything at all, and the discretisation is irrelevant (the number is
identical at every `n_z`) — **the reported +0.13 % does not move**. The
measurement constrains `K_p` only to within roughly an order of magnitude, and
constrains nothing else. Note in particular that the van 't Hoff value the page
offers as its uncertainty band gives a *better* agreement (+0.04 %) than the
reconstruction: the comparison cannot tell them apart.

Two further honesty points:

* On the quantity that actually varies the agreement is not 0.13 %:
  `1 − X` is 0.0030 measured against 0.00166 modelled — the model
  under-predicts the unconverted cyclohexane by **45 %**.
* "99.7 %" is quoted to three digits with no uncertainty, so X ∈ [0.9965,
  0.9975] and the deviation is only known to lie in [+0.08 %, +0.19 %]. A
  two-decimal-place deviation on a conversion pinned at 1 is not a meaningful
  precision.

The page's caveats are directionally right ("rides its thermodynamic ceiling …
tests the kinetics only weakly") but too soft, and every summary line —
notebook title block, `meta.yaml` `agreement:`, `README.md`, and
`models_entry.yaml` `description:` — leads with "reproduced to 0.13 % with
nothing fitted".

*Failure scenario:* a later page, or a reader, cites `H1.4` as an experimental
validation of Itoh's Langmuir–Hinshelwood kinetics or of the Sieverts
permeance. Neither was tested at all. Or: someone re-derives `K_p` differently,
gets a value 3× away, sees the same "excellent agreement", and concludes their
value is confirmed.

**Fix.** Print the sensitivity table (or two rows of it) in Validation check 1,
add the `1 − X` comparison, and reword the summary lines to something like:
"the measured 0.997 confirms that the reactor reached its co-current
fast-permeation asymptote; the model puts that asymptote at 0.9983. The
comparison bounds `K_p` to within about an order of magnitude and does not test
the kinetics or the permeance."

### F2 — MODERATE, CONFIRMED. The counter-current prediction is quoted from a solve that is not grid-converged, and the conservation check does not cover it

The counter-current numbers are computed at the class default `n_z = 400` with
no grid study and no conservation check. Measured:

| `v_A⁰` | n=200 | n=400 (**quoted**) | n=800 | n=1600 | n=3200 | n=6400 |
|---|---|---|---|---|---|---|
| 0.5e-6 | 0.20376 | **0.19906** | 0.19539 | 0.19299 | 0.19157 | 0.19078 |
| 8e-6 | 0.93046 | **0.93245** | 0.93338 | 0.93381 | 0.93400 | 0.93408 |

and the hydrogen identity `u_H(1) + v_H(0) = 3(1 − u_C(1))` closes to
**+1.4e-2** (0.5e-6) and **+3.3e-3** (8e-6) relative at `n_z = 400`, decaying
at first order (1.9e-2 → 2.0e-3 from n=200 to n=6400). The assembly is
therefore consistent — I verified it is not a bug, see below — but the printed
third decimals are discretisation, not physics.

Meanwhile the co-current numbers it is compared against (0.228 at 0.5e-6, 0.607
at 8e-6) are *exactly* the algebraic ceiling and identical at every `n_z`, so
the comparison sets a converged number against an unconverged one. The
qualitative conclusions survive and the reversal is in fact *larger* than
reported (0.228 vs ~0.190, not 0.199).

Related: Validation check 5 reports hydrogen conservation "to machine
precision" — but it is evaluated only on the stated co-current run, whose
outlet sits on the source's fixed point, where the identity is exact for any
discretisation. It proves nothing about the scheme; the counter-current case,
where it would have proved something, is not checked.

*Failure scenario:* CI or a reader re-runs at a different `n_z` and gets 0.190
instead of 0.199; or the "machine precision" conservation claim is taken as
covering the counter-current solve, which it does not.

**Fix.** Run the counter-current comparison at `n_z ≥ 3200`, or show its grid
ladder; quote two decimals; and report the counter-current hydrogen balance
next to the co-current one.

### F3 — MODERATE, CONFIRMED. "Three misprints" is an over-count

See verdict (b). Two of the three are misprints; the third is a reading that
the printed text already gives correctly. The page, `README.md`,
`meta.yaml`, the data sidecar and `models_entry.yaml` all say three, and
`meta.yaml` calls the misprint resolution "the archival contribution" of the
page. Asserting an error a published paper does not contain is exactly the
thing `docs/handoff.md` warns about (`G1.8`).

**Fix.** Recast item 2 as "a transcription confirmed, not a misprint": the text
prints 10⁻⁵ and reads cleanly at native resolution; Figure 4's hand-lettered
abscissa exponent is genuinely ambiguous but is fixed by the figure's own data
(the headline run is the marker at abscissa 11.7). Keep the ceiling argument —
it is excellent and it is what makes the ×10⁶ reading impossible — but present
it as corroboration, not as the resolution of a printing error. Change "three
misprints" to two everywhere.

### F4 — MINOR, CONFIRMED. Two stale numbers

* Notebook prose, PyMRM implementation section: *"6 to 14 iterations everywhere
  on this page"*. Instrumenting the page's own solver over the page's own 45
  solves gives **6 to 27** (values 6,7,8,9,10,16,19,21,27). This is not printed
  by any cell, so nothing catches it.
* `data/itoh-1987-stated-values.meta.yaml`: *"An independent van 't Hoff
  estimate … gives 1.6e11 Pa³"*. The notebook computes **1.49e11** (ratio
  1.58); 1.6e11 appears to be the ratio mis-transcribed as the value.

### F5 — MINOR, CONFIRMED. The α_H identity is sound but slightly over-claimed, and it is not load-bearing

Recomputed independently: `2π l₀ D C₀ / ln(r_o/r_i)` = 4.4683e-5 for
r_i = 8.5/r_o = 8.7 mm (−0.03 % vs the printed 4.47e-5) and 4.3646e-5 for the
OD reading (−2.36 %). With α_H printed to three digits, 2.4 % is decisive. The
identity is real.

New corroboration: in Figure 1 the "17" dimension arrows terminate on the
*inner* faces of the two membrane lines and the "28" arrows on the inner faces
of the shell wall, i.e. the drawing dimensions the tube ID. (Soft — the wall is
drawn schematically thick — but it is a second, independent witness.)

Two caveats the page should absorb: strictly, the identity shows what *Itoh
substituted into Eq. 1*, not what the tube physically was ("17.0 mm OD" may
simply have been mis-stated, or the α_H computed with the wrong radii); and the
choice is inconsequential, because the model uses the printed α_H and the only
consumer of `r_i` is `V_r` (4.9 % between the two readings), which affects no
reported number (see F1). Phrase it as "the printed α_H was computed with
r_i = 8.5 mm" and note it does not propagate.

---

## What survived the attack (no finding)

* **Determinism.** Re-executing `index.ipynb` reproduces every stream output
  **byte-identically**, `agreement.json` is unchanged, runtime 6.3 s. No
  warm-start continuation anywhere: the Newton initial guess is a deterministic
  explicit upwind march on the same grid, and the counter-current loop
  constructs each case from scratch. Complies with the `B1.1` lesson.
* **Every other number in the prose traces.** k(473) = 2.653e-5,
  K_B(473) = 1.160e-4, V_r = 3.18e-5 m³, K_p = 2.357e11, α_H recomputation,
  4.2e38 at 473 K, ceiling 0.72 / 0.9983, ×5.3 uplift, refinement ratio ~2,
  |ΔX| = 4.6e-5 at n_z = 800, 0.93 vs 0.61 at 8e-6, reversal at 0.5e-6 — all
  match the executed outputs. Only the two in F4 do not.
* **The pymrm assembly is correct.** I wrote an independent first-order upwind
  FV residual from scratch and evaluated it at the page's converged solutions:
  co-current relative residual 6.8e-15 at the stated run and 5.1e-9 at
  `v_A⁰ = 8e-6` — the two discretisations are identical. The counter-current
  interior matches to 1e-13; the only disagreement is at the two outflow faces
  and it is an O(h) difference in the zero-gradient boundary treatment, which
  is what F2 measures. The per-field velocity `(+1,+1,−1)` and per-field `bc`
  arrays do what the page says they do.
* **The regularised rate is algebraically identical** to Eq. 7 multiplied
  through by `p_H³`, and its `p_H → 0` limit `−k/K_B` is right.
* **The `NumJac` claim is true.** `pymrm/src/pymrm/numjac.py`:
  `dc = −eps_jac·|c|; dc[dc > −eps_jac] = eps_jac` with `eps_jac = 1e-6`, i.e.
  an absolute floor. For unknowns of order 1e-7 the Jacobian is indeed noise.
* **The counter-current convergence limit is genuine, and honestly reported.**
  The solve fails at `v_A⁰ ≥ 1e-5` with three different deterministic initial
  guesses (the inherited march, a loaded-shell guess, and a flat guess), so it
  is not an artefact of the co-current guess. No number outside the stated
  0.5–8e-6 range is quoted anywhere. (The likely cause is not the "deep-removal
  plateau" the page names but the `√v_H` derivative singularity at the
  counter-current shell inlet, where `v_H = 0` exactly.)
* **The model is bounded by its ceiling everywhere on the map** (verified with
  an independent Radau integration, not the FV solve). The bound is tight for
  the two lower feeds and loosens to 0.065 for `u_C⁰ = 1.64e-6` at the highest
  purge — which is the only part of the reported map where kinetics and
  permeance influence the answer at all.
* **Reproduction is not called validation.** The Figure 4 sweep and the
  counter-current comparison are both labelled correctly, and the
  `alpha_H = 0` equilibrium check explicitly says it is a consistency check.
  The decision not to digitise Figure 4's markers is stated and defensible.

---

## Recommendation

**Safe to publish after these fixes:**

1. **F1** — state the resolving power of the headline comparison quantitatively
   (α_H and k make no difference at all; `K_p` is bounded only to ~an order of
   magnitude), add the `1 − X` comparison, and reword the four summary lines
   that currently lead with "reproduced to 0.13 % with nothing fitted".
2. **F3** — "two misprints", not three; recast the `v_A⁰` item as a confirmed
   transcription rather than an error in the paper, keeping the ceiling
   argument as corroboration.
3. **F2** — recompute the counter-current comparison on a converged grid (or
   show its ladder), quote two decimals, and report its hydrogen balance;
   scope the "machine precision" conservation claim to the co-current run.
4. **F4** — fix "6 to 14 iterations" (actual 6–27) and the sidecar's "1.6e11"
   (actual 1.49e11).
5. **F5** — soften the OD/ID claim to what the identity actually shows, and
   note it does not propagate to any reported number.

Optional, both cheap and both strengthening: add the corrupted-exponent
candidate to the `K_P` decision table (the printed prefactor is
thermodynamically sensible; only the exponential argument need be wrong), and
record that the model reproduces Itoh's own calculated curves to rms 0.01 in
conversion — a check that, unlike the headline, does discriminate `K_p`.

*Verifier: Claude Fable 5, 2026-07-31. Nothing outside this file was modified;
no git command was run.*
