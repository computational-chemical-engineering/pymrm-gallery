# A4.4 — adversarial verification

Verifier pass, 2026-08-01. Source re-read independently at 600 dpi; every
headline number recomputed with code that shares nothing with the page.

## What was done independently

- **Equations re-read off my own 600 dpi renders** of PDF pp. 26, 27, 32, 33
  (= journal 886, 887, 892, 893): eqs. (82), (85), (86), (87), (107), (109),
  (110), the p. 887 scaling sentence, the p. 892 Jackson sentence, the three
  conditions preceding eq. (110), and the Fig. 44 condition box. The builder's
  transcription was **not** consulted first.
- **Closed form re-derived from scratch** in sympy, starting from eq. (109) as
  *printed* (with the ν_j/ν_1 flux-ratio identity), not from the builder's
  rearrangement.
- **Two independent solvers**: a `scipy.integrate.solve_bvp` collocation
  formulation of the binary dusty-gas pellet (no pymrm, no page code), and an
  RK4 shooting solution of the isobaric slab.
- **Break tests re-run**, plus five the builder did not try.
- **Notebook re-executed**: bit-identical output, 6.7 s.

## Transcription — no discrepancy

| eq. | as printed (my render) | page |
|---|---|---|
| (82) | −(1/RT)∇p_i = Σ_{j=1}^n (x_jN_i − x_iN_j)/Ð^e_ij + N_i/Ð^e_iM | same |
| (85) | Ð^e_iM = (ε/τ)(d_0/3)√(8RT/πM_i) | same |
| (86) | (N) = −(1/RT)[B^e]^−1(∇p) — **1/RT**, not c_t/RT | same |
| (87) | B^e_ii = 1/Ð^e_iM + Σ_{k≠i} x_k/Ð^e_ik ; B^e_ij = −x_i/Ð^e_ij | same |
| (107) | Σ N_i√M_i = 0 (∇p = 0; gaseous mixtures) | same |
| (109) | 1/D_1 = 1/Ð^e_1M + Σ_{j≥2}(x_j/Ð^e_1j)(1 − x_1ν_j/(x_jν_{1i})) − α'_1x_1(Ð_visc/Ð_1M)·[Σ(ν_i/ν_1)/Ð^e_iM]/(1 + Ð_visc Σα'_ix_i/Ð^e_iM) | same |
| (110) | 1/D_1 = 1/Ð^e_1M + 1/Ð^e_12 | same |

The `ν_{1i}` is a print artefact *on the page itself* (there is also a stray
`i` under the `j=2` sum limit), not an OCR failure — so the builder's two-way
symbolic confirmation was the right response and it is correct.

## Priority 1 — the contribution: CORRECT

From eq. (109) as printed, with ν_2/ν_1 = N_2/N_1 = −ρ:

    1/D_1 = 1/D_1M + (x_2 + ρ x_1)/D_12

which equals eq. (82) rearranged (residual identically 0, my own derivation),
collapses to eq. (110) at ρ = 1 (residual 0), and gives, with Kn ≡ D_12/D_1M,

    D_Bos/D_1 = 1 + (ρ − 1) x_1 / (1 + Kn)      [residual 0]

    lim_{Kn→∞} = 1        lim_{Kn→0} = 1 + (ρ−1)x_1
    d/dKn = x_1(1−ρ)/(1+Kn)²  → strictly monotone, never worst mid-transition

**Not an artefact of the normalisation.** The coefficient ratio is
normalisation-free. Re-running the pellet with D_AB held fixed and φ defined on
D_AB instead of D_Bos gives 26.96 / 26.74 / 24.75 / 14.2 / 2.7 / 0.30 % at
Kn = 1e-3 … 1e2 — the same monotone decay.

Independent collocation reproduces every headline number:

| quantity | page | independent |
|---|---|---|
| p_0/p at Kn=1e6, ν_B = 1,2,3,4 | dev 4.23e-7 / 7.69e-7 / 1.08e-6 | 4.22e-7 / 7.69e-7 / 1.08e-6 |
| η error, Kn = 1e-3 / 1 / 1e6 | +27.03 / +14.14 / +3.0e-7 | +26.96 / +14.12 / +0.0000 |
| p_0/p bulk, DGM / Bosanquet | 1.00029 / 1.99941 | 1.00029 / 1.99941 |
| isobaric N_A (dgm), Kn=1e-3 | 1.36355 | 1.363554 (RK4 shooting) |
| Bosanquet flux error, Kn = 1e-3 / 1 / 1e3 | −41.33 / −19.17 / −0.04 % | −41.330 / −19.175 / −0.037 % |
| Graham ratio −N_B/N_A | 0.268256 (dev 2.1e-11) | 0.268256 exactly |

## Priority 2 — the Jackson check: agreement, honestly labelled

The page tests the **formula** √ν_B and says so: cell 21 prints "The printed
40 % is the authors' rounding of their own expression, not a different
number". `meta.yaml` and the sidecar say the same. Not a spun discrepancy.
The p. 892 sentence reads exactly `p_0 = √(ν_B) p` … "Thus, for ν_B = 2, we
have a 40% increase in pressure as we proceed towards the centre of the pellet
(Jackson, 1977)". √2 − 1 = 41.42 %; 40 % is a rounding.

## Priority 3 — reprint test: PASSES, all four claims verbatim

1. Named — "a relation usually referred to as the Bosanquet formula. As noted
   above this formula is very restricted in its applicability." ✓
2. Derived — "For the special case of a (i) binary mixture, (ii) with no net
   change in the number of moles, and (iii) satisfying eq. (107), eq. (109)
   simplifies to" ✓
3. Three conditions verbatim ✓
4. Used as a foil — Kaza & Jackson's uphill diffusion "impossible to explain
   with say the Bosanquet formula (110)" ✓ (the Elnashaie & Abasher / Reddy &
   Murty sentence is about "approximate forms of eq. (109)" generally, so the
   *second* use is adjacent rather than explicit; the case-YAML "twice" is a
   mild overstatement, the page itself does not claim it)

Scope split is clean. A4.3 keeps the n-component matrix form, the viscous
term, Fig. 45 and Fig. 44 — the only experimental comparison in the section.
A4.4 leans on none of it. No overlap with A4.2 (S9, bulk, no pore).

## Priority 4/5 — break table and the self-caught defect

Every row reproduced. The M_A/M_B fix is complete: `isobaric_slab` passes
`m_a_over_m_b = MOLAR["H2"]/MOLAR["N2"]`, `dbm = dam·(M_A/M_B)^{1/2}`, which is
what eq. (85) requires; `RHO_GRAHAM = (M_A/M_B)^{1/2} = 0.268256` matches
eq. (107) exactly. No quoted number retains the inverted ratio.

## Findings

Severity order. All CONFIRMED unless noted.

1. **The √ν_B check is blind to a sign flip in the off-diagonals of [B^e], but
   cell 22 claims it tests "the sign structure of eq. (82)".** Measured:
   flipping B_01/B_10 leaves p_0/p at 1.4142130, dev 4.23e-07 — unchanged to
   every printed digit. The same defect moves the dgm-vs-eq109 comparison to
   8.9e-1 and Graham's residual to 2.6e-2, so the page's other checks *do*
   catch it — but the sentence attached to the headline asserts a sensitivity
   it does not have, which is the `B1.6`/`E1.2` failure mode. Also drop "the
   ideal-gas step ∇p_i = RT∇c_i" from the same list: the solver works in
   concentrations throughout, so no such step exists in the code and nothing
   could break it. Fix: move the friction-term signs into V3's "cannot see",
   with the measured 4.23e-07, and point at V1/Graham.

2. **"`dgm` and the other two share no line of code" (cell 14) is false.** All
   three closures use `Cell.faces`, the same `construct_grad`/`construct_div`,
   the same residual assembly and the same `NumJac`; only the constitutive
   closure differs. The constitutive independence is real (break 1 above), but
   the comparison cannot detect a discretisation error — both routes move
   together at n=5. Cell 30 item 5 states this correctly; cell 14 overstates.
   Fix: "share no *constitutive* code".

3. **The headline 27.0 % carries two hidden parameters.** Grid: 27.03 % at the
   sweep's n = 300, converging to 26.96 % at n = 4800; the third significant
   figure printed in the sweep table and stored in `agreement.json`
   (0.2702532880477573) is discretisation. Thiele modulus: the same error is
   25.5 % at φ=10, 26.6 % at φ=20, 27.0 % at φ=30, 27.5 % at φ=100 — and the
   grid study that licenses n = 300 is run at φ = 10, Kn = 1, not at the
   sweep's φ = 30, Kn = 1e-3, where the boundary layer is three times thinner.
   `README.md`, `meta.yaml adds`, and `models_entry.yaml description` all quote
   "27.0 %" with no φ. Fix: name φ = 30 wherever the number appears, and run
   the sweep at n ≥ 600 (total runtime is 6.7 s, so there is room) or quote two
   significant figures.

4. **The Jackson deviation is a finite-Kn model residue, not solver accuracy,
   and the page does not say so.** It scales exactly as 1/Kn — 4.225e-06 at
   Kn=1e5, 4.225e-07 at 1e6, 4.225e-08 at 1e7, 4.227e-09 at 1e8 — and my
   independent collocation returns the same 4.22e-07 at Kn=1e6, so it is not
   discretisation. Cell 21's "the solve reproduces the expression to 4.2e-07
   relative" reads as numerical fidelity. Saying it is the residual bulk term,
   vanishing as 1/Kn, *strengthens* the check.

5. **The "unconverged solve" blindness is under-demonstrated.** With
   `maxfev=1` the answer is already right to seven digits, and `maxfev=2`
   converges (residual 4.58e-15) — the Knudsen-limit pellet is nearly linear,
   so one Newton step nearly solves it. This is not the `B1.6` class (residual
   1e-11 on a physically impossible solution). The conclusion drawn ("never
   infer convergence from an identity"; assert residuals everywhere) is right
   and the practice is right; only the demonstration is weak, and it errs
   toward claiming *less*.

6. **"35.2 °C and 1 atm" is not in Krishna & Wesselingh.** The string "35.2"
   appears nowhere in the paper. P. 872 prints only "D12 = 8.33 x 10⁻⁵,
   D13 = 6.8 x 10⁻⁵, D23 = 1.68 x 10⁻⁵ m² s⁻¹", said to be "estimated from the
   kinetic gas theory", with no temperature or pressure. The conditions are
   Duncan & Toor's (1962), recorded in `A4.9`'s dataset sidecar. Cells 3 and 6
   attribute them to the review's worked example, and the References section
   does not cite Duncan & Toor at all. `T_REV = 308.35` is used quantitatively
   (v̄, the Kn = 1 pore sizes 138.9/113.4/104.4 nm, the transition-curve
   figure, c_t in the scaling study), though no validated number depends on it.
   One clause fixes it.
   *Adjacent, on the published `A4.2`*: its sidecar says the review's values are
   "identical to the values Duncan & Toor **measured**", where K&W say
   "estimated from the kinetic gas theory". Worth checking there.

7. **`EPS_TAU` reads the printed value of τ and uses it as ε/τ.**
   `EPS_TAU = P[("eq84","tau_cylindrical")] = 1.0`, then passed as `eps_tau`.
   Numerically right only because ε = 1 is separately assumed; the prose says
   so, the provenance chain does not. Cosmetic.

8. **No executed check confirms the geometry index `nu`.** Cell 37 says the
   page's "geometry claim rests on the grid study and the effectiveness
   factors" — but neither is compared to a geometry-dependent external
   reference, so nothing on the page would catch `nu=1` where `nu=2` was meant.
   (Independently confirmed correct: my collocation sphere gives η = 0.076140
   against the page's 0.076139 at n = 1200.) Either delete the sentence or add
   η against the classical (3/φ)(coth φ − 1/φ) in the ρ = 1 limit, which is
   free.

9. **Trivial.** The dataset sidecar says Fig. 44's abscissa runs "10 to
   10^4 Pa"; the printed figure runs 10 to 10^5 Pa. All eight transcribed
   Fig. 44 conditions are correct on the render.

10. **For the integrator, not a page defect.** `queue_cases/A4.3.yaml` is
    `status: needs-paper`, blocked on Mason & Malinauskas (1983). A4.4 hands
    A4.3 the n-component matrix form, the viscous term and Fig. 44 — all of
    which K&W print in full on the pages A4.4 already read. The reprint route
    that carries A4.4 unblocks A4.3; re-triage rather than leave it waiting.

11. **Strengthening, not a defect.** Eq. (109) for a binary is
    `1/D_1 = 1/D_1M + (x_2 + ρx_1)/D_12` — Bosanquet's *series-resistance
    structure is exact*; the only thing eq. (110) gets wrong is the bulk
    resistance, which carries a flux-ratio factor. That is *why* the error
    vanishes as the bulk resistance becomes negligible, and it means the
    "interpolation is worst mid-transition" frame is not merely reversed, it is
    the wrong frame: the interpolation is never the problem. Cell 13 gives the
    right physical explanation and stops one sentence short of this.

## What survived unchanged

- Every equation transcription, checked against my own renders.
- The closed form, its two limits and its monotonicity — symbolically and
  under two normalisations.
- All five reported observable numbers, against two independent solvers.
- Graham's law as a prediction (2.1e-11) and its break (0.155).
- The two p. 887 scaling slopes, including the *unshown* claim that the
  residues shrink deeper into each limit: 0.9927 → 0.9993 → 0.99999 (Knudsen)
  and 0.0333 → 0.00036 (bulk). True.
- Every row of the deliberate-break table.
- Prose/output drift: **none**. No computed number is typed into a markdown
  cell; re-execution is bit-identical.
- pymrm conventions: `NumJac((n,2), axes_diagonals=[0], axes_blocks=[-1])` is
  correct here — the divergence is inside the differentiated function, so the
  residual genuinely reads neighbours, and `ndims = 2` makes
  `axes_diagonals` meaningful. Outward-normal BCs commented; `nu` commented.
  Independent collocation agreeing to 6 digits is the strongest available
  confirmation that the stencil is not silently wrong.
- Tier 6 stated repeatedly; Fig. 44 not leaned on; cross-page load honest and
  Colab-safe. `check_metadata.py` and `check_agreement.py` both clean.

## Verdict

**Safe to publish after the fixes in findings 1, 2, 3 and 4** (1, 2, 4 are
prose; 3 needs the sweep re-run at n ≥ 600 and φ named in the three metadata
files). 6, 7, 8, 9 are one-line corrections worth folding in. Nothing found
that changes a conclusion.
