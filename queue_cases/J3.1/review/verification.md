# J3.1 — adversarial verification

Verifier pass on the staged page at `queue_cases/J3.1/page/`, against
`~/papers/pymrm-gallery/Doyle_1993_J._Electrochem._Soc._140_1526.pdf` read at
600 dpi, and against the published `pages/J3.4-doyle-fuller-newman/` and
`pages/J3.5-single-particle-model/`.

**Verdict: safe to publish after the fixes listed under F1–F4.** F1 and F2 are
factual errors in prose and must be fixed. F3 and F4 are overclaims about what
two of the checks can detect. Nothing found invalidates a result.

---

## 1. Scope decision — the page is NOT covered by J3.4. Build it.

Each of the builder's four claims was checked against the actual files, not
accepted:

| claim | verdict |
|---|---|
| J3.4 embeds Eq. 17 and never examines the law | **true** — Eq. 17 and Eq. 16 appear in J3.4's markdown and in `Cell.i_transfer`; there is no discussion of the law's structure anywhere |
| J3.4 uses Eq. 30 alongside Eq. 17 without connecting them | **true** — `i0_cathode` (Eq. 30) is called *only* from `kappa_from_nu`, never in the kinetics, and nothing relates it to Eq. 17 |
| neither page asks what happens when α_a+α_c ≠ 1 | **true** — no occurrence in J3.4 or J3.5 |
| J3.4 touches neither limit | **true** — `2·RTF·arcsinh(I/2i01)` is the *exact* inversion of the symmetric law, not a limit; no Tafel, no linear form |

Two qualifications the builder did not state, neither of which changes the
decision:

* `J3.5` already ships `bv_invert` (its line 592), a closed-form inversion of
  Eq. 17 whose algebra silently *contains* the reduction J3.1 derives: at zero
  current its root reproduces Eq. 16 exactly, and its exchange prefactor is
  exactly Eq. 30's structure. J3.5 never says so, never states the α condition,
  and never uses it as a check — but an integrator should know the algebra is
  already in the repository.
* `J3.5`'s `C_r_doyle` is `I/(a·i_0,typ·δ_c)` — the same interfacial-versus-
  superficial ratio J3.1 uses to classify the cathode. J3.5 reports ~1e−4 there
  against J3.1's 1.0e−3, because it scales on `PRE0·c_T` rather than Eq. 30 at
  the initial state. Not a contradiction, but two numbers in the gallery for the
  same idea.

What J3.1 adds is not thin: the derivation, the α_a+α_c = 1 condition, both
closed-form thresholds, the asymptotes and the non-monotonicity, and the η_s
distribution of the J3.4 cell. A section on J3.4 would not carry it. Keep
`J3.1` as its own page.

---

## 2. Transcription and the Eq. 17 → Eq. 30 collapse — genuine, non-circular

I read the equations myself off 600 dpi renders (PDF pages 4 and 6 = journal
1528 and 1530; the article begins on PDF page 2).

* **Eq. 16** (p. 1528, left column):
  `U₂ = U₂^θ − U_ref^θ + (RT/F)(ln((c_T−c_s)/c_s) + βc_s + ζ)` — matches the page.
* **Eq. 17** (p. 1528, right column, breaking across the column):
  `i = Fk₂(c_max−c)^{α_c} c^{α_A} [c_s exp(α_aF/RT (η−U')) − (c_T−c_s) exp(−α_cF/RT (η−U'))]`
  — matches the page exactly, including the typeset `α_A`.
* **Eq. 30** (p. 1530, right column):
  `i_o2 = F(k₂)(c_max−c)^{α_c2}(c)^{α_a2}(c_T−c_s)^{α_a2}(c_s)^{α_c2}`
  — matches the page's derived form **factor for factor**, including which α sits
  on which concentration.
* Eqs. 18, 19, 20, 28, 29 and Table II also match; Table II matches J3.4's
  parameter CSV row for row (R_s = 1.0 µm, ε = 0.3, c_T = 29 000, k₂ = 1e−10,
  σ = 1e4, i_o1 = 12.6, α = 0.5, δ_s = 50 µm, δ_c = 100 µm, c⁰ = 1000, T = 100 °C).

**Non-circular: confirmed.** Eq. 30 is printed two journal pages after Eq. 17,
introduced by "The exchange current density in the cathode can be determined from
the reaction rate parameter k₂ through [30]". The paper does **not** derive it
from Eq. 17. The two transcriptions are independent, and the identity is a real
mutual check.

**"No other reading of α_A reproduces Eq. 30": confirmed.** Break test — coding
Eq. 17's prefactor as `c^{α_c}` while Eq. 30 keeps the printed `(c)^{α_a2}` moves
the residual from 4.5e−12 to **8.2e−1**. The claim holds.

**α_a+α_c = 1 dependence: confirmed by deliberate break.** Off the condition the
residual is 1.0–1.8 (relative) and the Eq.-16 offset is 15–168 mV, matching its
own closed form to 7.3e−16 V.

---

## 3. Deliberate-break results on the identity battery

Break tests run at `scratchpad/J31/attack.py`. "Baseline" is the shipped value.

| check | baseline | breaks under | verdict |
|---|---|---|---|
| `eq17_vs_eq30_max_rel` | 4.6e−12 | Eq.17 prefactor α swap → 8.2e−1; Eq.30 prefactor α swap → 8.2e−1; Eq.30 tail swap → 8.8e−1; Eq.17 bracket prefactors swapped → 2.0 | **strong evidence** |
| `zero_current_at_zero_eta_rel` | 7.8e−15 | bracket prefactors swapped → 1.9e+2; closed form scaled ±10 % → 5.9e−1 | **powered** (blind to the prefactor — correctly, it cancels) |
| `alpha_offset_closed_form_V` | 7.3e−16 | a wrong closed form breaks it | powered, but an exact algebraic consequence |
| `partial_currents_balance_rel` | 1.5e−14 | *nothing tried moved it* | **structural; cannot fail** — see F4 |
| `i0_from_partials_vs_eq30_rel` | 4.4e−16 | Eq.30 tail exponents swapped → **still 4.4e−16** | **decoration; cannot fail** — see F4 |
| `i0_general_vs_eq30_rel` | 8.9e−16 | same swap at (0.3, 0.7) → **7.13** | **powered** — this is the version that works |
| `Rct_identity_max_rel` | 1.8e−9 | this is the central-difference truncation floor, not an identity residual | weakly powered (tests the closed-form equilibrium potential) |
| `arcsinh_inversion_max_rel` | 9.7e−11 | prefactor RTF/(2α) → 7.9e−1; a 1 % perturbation → 3.1e−2 | **powered**, and not circular (Eq. 17 evaluated forward) |
| `tafel_threshold_closed_vs_numeric` | 3.1e−15 | using n → α_a in the closed form breaks it | powered-but-structural |
| `linear_sign_convention_rel` | 0.0 | — | structural, and **already labelled as such on the page** |
| `salt_conservation` | 8.9e−16 | — | structural, and **already labelled as such on the page** |

One point worth a sentence on the page: **at Doyle's own α_a = α_c = 0.5 every
exponent swap I tried is a no-op** (all break tests read 4.6e−12 on the (0.5,0.5)
row). The entire discriminating power of Check 1 comes from the hypothetical
asymmetric rows the sweep adds. That is legitimate and clever — but it should be
said, because a reader will assume the check bites at the published parameters.

---

## 4. The −15.7 % — the number stands, the stated reason does not

* **The ratio really does cancel the unprinted group.** Confirmed algebraically
  against the printed Eq. 28 and Eq. 29 (both carry `(1/κ + 1/σ)` identically),
  so `ν²/δ = (n/α_a)·a·i_o2·δ_c/I`. σ contributes 1e−4 against 1/κ ≈ 20 (bulk)
  or 125 (effective) — negligible either way. **Independent of J3.4's κ
  reconstruction: confirmed.**
* **δ = 1.95 and ν = 68 read correctly at 600 dpi.** Both appear in running text
  immediately below Eq. 30. Confirmed.
* **Rounding analysis reproduced.** ν ± 0.5, δ ± 0.005 → 2331–2412 against 2371,
  i.e. −1.7 % to +1.7 %; closing the gap needs ν = 62.4, 8.2 % off. Confirmed.
* **I could not explain the residual either.** Closing it needs i_o2 = 56.46
  against Eq. 30's 47.57 — which would require u₀ = 0.0141 (Fig. 2's caption
  prints "the initial concentration in the solid was 1 % of maximum") or
  R_s = 0.843 µm (Table II prints 1.0 µm). No printed value accounts for it.
  Leaving it unexplained is right.
* **It is arithmetically J3.4's 19 %, exactly.** I recomputed:
  κ_eff(Eq. 28)/κ_eff(Eq. 29) − 1 = **18.7 %**, and 1/(1 − 0.1574) − 1 = **18.7 %**.
  The page discloses this in a parenthetical; that is honest, but a reader
  skimming will take it as a *second, independent* test when it is one
  disagreement viewed from two sides. Promote it out of the parenthetical.

**F2 (moderate, must fix).** The cell prints: *"the paper does not say which
concentration state or which effective kappa the groups were evaluated at"*.
Both halves are wrong:

1. Doyle p. 1530, in the sentence immediately after Eq. 30: **"For these
   calculations, the concentrations are taken to be at their initial values."**
   Fig. 2's caption pins the initial solid state at 1 % of maximum. The page's
   own sidecar already records `conditions: initial concentrations` for both the
   `delta` and `nu` rows.
2. κ **cancels identically** from ν²/δ, as the page argues two paragraphs
   earlier, so an ambiguity in κ cannot be where the residual sits.

*Failure scenario:* a reader checks the paper, finds the state printed in plain
English, and concludes the page did not read its source carefully — on the one
comparison that uses the paper's numbers rather than its algebra. Replace with:
the paper *does* print the state, and Eq. 30 evaluated there does not reproduce
the printed groups; the residual sits in i_o2, in the specific area a (which
Table II does not print and which must be reconstructed as 3(1−ε)/R_s), or in
the paper's own arithmetic. The same sentence in `meta.yaml` (~line 93) needs
the same fix.

---

## 5. Validity thresholds — derivations and numbers check out

* Tafel deviation `ε_T = 1/(e^{(α_a+α_c)fη}−1)` verified by hand; the inversion
  `η = RT ln(1+1/ε)/((α_a+α_c)F)` is exact, and the sum-only dependence is a
  property of the algebra, not an asymptotic claim. 77.1/97.9/148.4 mV at
  373.15 K and 61.6/78.2/118.6 at 298.15 K all reproduce.
* Linear expansion `ε_L ≈ −bη̂ + (b²−c)η̂²` with `b = (α_a−α_c)/2`,
  `c = (α_a³+α_c³)/6` verified by series; both asymptotes follow.
  51.7/35.9/15.8 mV and the 9.0× narrowing reproduce.
* **The non-monotonicity is NOT a root-finding artefact — confirmed.** I re-ran
  the threshold with a 100× denser scan (2 000 000 points, floor 1e−11 instead
  of 1e−9). Agreement with the shipped algorithm is 1e−15 relative at **every**
  α_a tested, including 0.415 and 0.410 where the jump happens:

  | α_a | 0.50 | 0.45 | 0.42 | 0.415 | 0.410 | 0.405 | 0.30 |
  |---|---|---|---|---|---|---|---|
  | 5 % threshold / mV | 35.86 | 61.15 | **81.33** | 85.03 | **32.65** | 24.99 | 8.29 |

  The mechanism is confirmed independently: the first-lobe overshoot is
  **1.516 %** at α_a = 0.45 (page prose: "only 1.5 %") and first exceeds 5 % at
  α_a ≈ 0.4103 — exactly where the threshold collapses from 85.0 to 32.7 mV.
  The page is right to print it and then recommend the asymptote.
* **F5 (minor).** The page says the 1 % jump "sits near α_a = 0.455". Recomputed
  on a 501-point grid it sits between **0.4592 and 0.4594**; at 0.455 the
  overshoot already reaches 1.225 %, so the window has not jumped yet. This is
  one of only two prose numbers not interpolated from output.
* **F6 (minor, framing).** "The two electrodes of the same cell are in different
  regimes" is tolerance-dependent. At the stated I = 10 both are inside the 5 %
  linear window: cathode 0.032 mV (linear error 4e−6 %), foil 24.9 mV (2.45 %).
  It is unambiguous only at the 1 % tolerance, or at I = 20 where the foil is
  8.3 % out. All the numbers are printed correctly and the 519 mV vs 1.5e−10 mV
  contrast is real — but the section heading is stronger than the I = 10 data.

---

## 6. The three self-caught defects — all genuinely fixed

1. **Bracketing.** Confirmed by reproducing the bug: `brentq(dev, 1e-9, 4.0)` at
   α_a = 0.3 returns **203.4 mV**, exactly the reported wrong number; the shipped
   scan-then-bracket returns 8.29 mV. There are 3 sign changes on (0, 4]. The
   shipped `eta_linear_threshold` scans for the first crossing and raises if
   there is none. **Fixed.**
2. **Time stepping.** `march(dt_seq=seq)` replays the recorded sequence and
   explicitly refuses to shrink a prescribed step (`if dt_seq is not None:
   break`). `dv()` compares on `min(len)` with no interpolation. The `bv`
   control run returns 2.2e−12 mV, which is the right null. **Fixed.**
3. **Both-electrode switch.** `march` calls
   `cell_potential(y, I, anode_law="bv")` unconditionally, so the anode always
   uses the exact law and `kinetics` acts on the porous cathode alone. **Fixed.**
   (Note for the record: column 6 of `out` still uses `self.kinetics` for the
   anode overpotential, but that column is never used in the substitution
   comparison — only in the `env` table, which is built from `eq17` runs only.
   Harmless, but a reader of the code could trip on it.)

---

## 7. Routine

* **Reproducibility: perfect.** `build_page.py` re-run in a temp directory
  regenerates all 42 cells with **zero** source differences from the staged
  `index.ipynb`. Re-executing the staged notebook (50 s) gives **zero** cells
  with differing stream output — byte-identical. No drift.
* **Prose against output.** I checked roughly forty numbers: 77.1/97.9/148.4,
  61.6/78.2/118.6, 51.7/35.9/15.8, 41.3/28.7/12.6, 9.0×, 3.6 %, 6.8 %, 81/36 mV,
  0.032/0.74 mV, 24.8/46.6 mV, 3.705 mV at u = 0.40, 1.5e−10 mV, 519 mV,
  +446.2/+444.2 mV, 0.46 %, 0.6/4.2 mV, aδ_c = 210, 48 A/m², 0.05 A/m², 15.7×,
  1.9×, 0.4 %, 49.4 mV, 1.85 %, ±1.7 %, ν = 62.4, 8.2 %, −15.7 %. **All match.**
  The only two not interpolated are the 1.5 % overshoot (correct) and α_a = 0.455
  (F5).
* **Tier 6.** Confirmed on the render: Fig. 2 carries only smooth computed curves
  plus a dashed open-circuit line, no markers, and its caption describes only
  simulation settings. Nothing on the page is called experimental validation, and
  the reproduced/validated distinction is stated three times.
* **Peak η_s flagged as non-converged.** Correctly flagged and correctly not used
  alone; the 15.7× margin against 1.9× of movement carries the conclusion.
  **F7 (very minor):** the four grids give 0.530 → 0.738 → 0.904 → 1.008 mV,
  increments 0.208, 0.166, 0.104 — a *converging* sequence (Richardson ≈ 1.15–1.2
  mV), so "does NOT converge under refinement" overstates. "Not converged at
  these grids" is accurate. The error is in the conservative direction.
* **ν = 2 open question.** No side taken; run both ways; 49.4 mV / 1.85 %;
  classification unchanged. Correct.
* **Conventions.** `η = Φ₁ − Φ₂` with Φ₁ gauged to zero, `η_s = η − U₂` —
  consistent with Eq. 19 and Eq. 16. `bc` is homogeneous Neumann on the outward
  normal at both ends with the physical fluxes written onto the face arrays and
  commented; `construct_div(..., nu=0)` is commented as Cartesian. Ported from
  J3.4 unchanged. No issues.
* **Metadata.** `slug` and `title` identical between `meta.yaml` and
  `models_entry.yaml`; the notebook front-matter title differing from
  `meta.yaml` matches the house pattern (J3.4 and J3.5 both do it). `S10` is
  justified by the taxonomy (DAE/electroneutrality). `models_entry.yaml` carries
  `status: published` while `meta.yaml` carries `in-progress` — flagged in the
  case notes as an on-merge flip; the integrator must not miss it.
* **F8 (very minor, citation precision).** Verified in the Marquis PDF: the
  hyperbolic-sine kinetics is **Eq. 1g** (`j = j₀ sinh(Fη/2RT)`); **Eq. 1h** is
  the exchange current density `j₀ = m(c_s)^{1/2}(c_s,max−c_s)^{1/2}(c_e)^{1/2}`.
  The page and `meta.yaml` cite "its Eq. 1h" for the kinetics — should be
  "Eqs. 1g–1h". Eq. 18 for the sinh⁻¹ inversion is correct, and the Reuse claim
  that Marquis prints the same i₀ structure with the (c_max−c) factor dropped is
  correct.

---

## Findings, ranked

**F1 — moderate, must fix. The Reuse section's transfer claim is wrong by ~21×.**
CONFIRMED. The page says: *"This cell's factor between the two is aδ_c = 210,
which moves it from the edge of the Tafel regime to three decades inside the
linear one."* With the **superficial** current density, i/i₀ = 0.210,
η_s = 6.75 mV, linear error **0.18 %**, Tafel error 428 %. That is inside even
the **1 %** linear window (15.81 mV). The 5 % Tafel boundary needs i/i₀ = 4.4,
i.e. 21× more current. aδ_c = 210 moves the cathode from *well inside* the linear
regime to *further* inside it; it never approaches Tafel.
*Failure scenario:* this is the sentence a reader lifts into their own work. A
reader with a porous electrode is told that ignoring aδ_c would put them at the
edge of Tafel; they will mis-diagnose their cell in exactly the direction the
page is trying to warn against. The equivalent sentence in *Results* ("three
decades inside the linear regime … because aδ_c = 210 divides the superficial
current") is correct — only the Reuse framing is wrong.

**F2 — moderate, must fix.** The ν²/δ caveat gives a reason the paper
contradicts, and a second reason its own algebra excludes. See §4. Fix in the
notebook cell and in `meta.yaml`.

**F3 — minor, fix by rewording.** "Which is Doyle's Eq. 30 *exactly*, factor by
factor and exponent by exponent" overstates Check 1's coverage. CONFIRMED by
break test: `i_eq17` and `i0_eq30` write the `(c_max−c)^{α_c}c^{α_a}` prefactor as
the *same code expression*, so it cancels in the ratio — swapping α_a↔α_c there
in **both** functions leaves the residual at 4.5e−12, unmoved. Only the
`(c_T−c_s)^{α_a}c_s^{α_c}` half of Eq. 30 is actually tested against Eq. 17's
bracket. (Swapping it in *either one alone* does break, which is why the α_A
claim survives — see §2.) Say what the check covers, and add the sentence that
its power comes entirely from the asymmetric rows because every swap is a no-op
at Doyle's α = 0.5.

**F4 — minor, fix by labelling.** Two of the six Check-2 metrics cannot fail, and
they sit under a clause saying they "break immediately" on a mis-transcribed
coefficient. CONFIRMED: `i0_from_partials_vs_eq30_rel` reads 4.44e−16 both before
and after swapping Eq. 30's tail exponents, because it is evaluated only at
α_a = α_c = 0.5 where the swap is a no-op, and both sides share the `pre`
expression. `partial_currents_balance_rel` is the same algebraic statement as
(a) at n = 1 and is exact for any concentration. The heading is honest
("structural, not numerical"); the clause after it is not. The powered version is
already two lines below (`i0_general_vs_eq30_rel`, 8.9e−16 → 7.13 under the same
swap) — point at it.

**F5 — minor.** "at the 1 % level the jump sits near α_a = 0.455"; computed
0.4593. §5.

**F6 — minor, framing.** "Different regimes" holds at the 1 % tolerance and at
I = 20, not at the 5 % tolerance at the stated I = 10. §5.

**F7 — very minor.** "does NOT converge under refinement" for the peak η_s; the
sequence is converging, just not converged. §7.

**F8 — very minor.** Marquis Eq. 1g vs 1h. §7.

---

## Things I attacked and could not break

* The transcription of Eqs. 16, 17, 18, 19, 20, 28, 29, 30 and Table II —
  re-read independently at 600 dpi, all correct, including which α sits on which
  concentration in Eq. 30.
* Circularity of the Eq. 17 → Eq. 30 collapse — Eq. 30 is printed separately two
  pages later and is not derived from Eq. 17 in the paper.
* The α_a+α_c = 1 condition and its closed-form failure.
* The non-monotonic linear threshold as a root-finding artefact — survives a
  100× denser independent scan to 1e−15.
* All three self-caught defects — reproduced the old wrong answers, confirmed the
  fixes.
* Reproducibility — zero drift, in either direction, on a full re-run.
* Independence of ν²/δ from J3.4's κ reconstruction.
* The tier-6 statement and the absence of any experimental claim.

---

## Recommendation

**Safe to publish after fixing F1 and F2 (factual) and rewording F3 and F4
(overclaim scope).** F5–F8 are cheap and worth taking in the same pass. The
scope decision is sound, the Eq. 17 → Eq. 30 collapse is a genuine and strong
non-circular transcription check, the −15.7 % is correctly left unexplained, and
the page's headline contribution — the closed-form validity thresholds and the
operating envelope — survived every attack I could construct.
