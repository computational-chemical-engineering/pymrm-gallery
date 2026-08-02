# A4.3 — adversarial verification

Verifier pass, 2026-08-02. Source re-read independently off my own 600 dpi
renders of PDF pages 24–33 (= journal 884–893). Notebook re-executed (12 s,
`agreement.json` reproduced **bit-identically**). An independent solver of
eq. (82) was written from scratch (fixed-step RK4 marching + Newton on the three
constant fluxes; no pymrm, no `solve_bvp`, no shared helper, constants read off
my own renders).

**Verdict: send back.** Two of the four headline new results are artefacts of a
parameter held fixed that cannot be held fixed, and one stated validation
coverage claim is false. Nothing is fabricated; the transcription is clean and
the numerics are sound.

---

## What survives the attack

- **Every equation transcribed is correct.** I read eqs. (82), (83), (84), (85),
  (86), (87), (91), (99)–(101), (103), (105), (106), (107), (108)–(110) off my own
  renders. All match the page character for character, including eq. (86)'s
  printed `1/RT` (the OCR's `c_t/RT` is indeed wrong) and the sentence below
  eq. (85), "where d₀ is the pore diameter and the square-root term represents
  the velocity of motion". The p. 887 exponent sentences match the sidecar
  verbatim. The Jackson √ν_B quote (p. 892) and the Haynes and Kaza & Jackson
  sentences (p. 893) are quoted accurately.
- **V1 is real, and stronger than the page claims.** My independent RK4 solver,
  using my own reading of eqs. (83)/(85) and my own pair matrix, reproduces
  pymrm's fluxes at every pore size:

  | d₀ | indep. vs pymrm n=100 | indep. vs pymrm n=1600 |
  |---|---|---|
  | 1 mm | 8.694e-06 | 1.05e-08 |
  | 10 µm | 8.542e-06 | 3.29e-08 |
  | 100 nm | 3.965e-07 | 1.53e-09 |
  | 10 nm | 8.983e-09 | 4.31e-11 |
  | 1 nm | 1.034e-10 | 1.95e-11 |

  The first column reproduces the page's V1 numbers exactly, and the second shows
  the deviation is discretisation error at *every* pore size, not only at the
  100 nm where the page ran its order study. Because my route shares no helper,
  this also independently confirms `knudsen_D` and `pair_matrix` — i.e. the
  declared blind spot is empty in fact, though the page cannot know that.
- **The declared V1 blind spot is real.** Breaking eq. (85)'s mass exponent
  inside the *shared* helper (both routes) gives V1 = 8.4e-09, i.e. completely
  blind. Breaking the prefactor d₀/3 → d₀/2 in both routes gives 2.75e-06, which
  the page's own >100× criterion would print as "BLIND". Correctly declared.
- **The singularity claim is CONFIRMED, and provable.** With the wall term
  deleted, `Σᵢ Bᵢⱼ = Σ_{k≠j} x_k/Ð_jk − Σ_{i≠j} x_i/Ð_ij = 0` because Ð is
  symmetric: **(1,…,1) is an exact left null vector for any composition and any
  n ≥ 2**, so [Bᵉ] is singular identically, not accidentally. Numerically:
  |det| ~ 1e-19, cond ~ 2e16 on random n = 4 draws. With the wall term the column
  sums become 1/Ð^e_jM > 0, making [Bᵉ] **strictly column-diagonally dominant**,
  hence nonsingular by Levy–Desplanques — so the page's Background claim
  ("the wall term is what makes the matrix invertible") is a theorem, not a
  numerical observation. `np.linalg.solve` raised `Singular matrix` at every
  grid I tried (5, 17, 50, 100, 200 cells), so the "raises rather than returning
  a wrong answer" behaviour is robust. The physical reading is also right: summing
  eq. (82) over i gives ∇p/RT = Σᵢ Nᵢ/Ð^e_iM, which without the wall term forces
  ∇p = 0 — Graham's law — so the system is consistent only on that subspace.
- **Every break-table row reproduces.** 19.416 / singular / 0.99674 / 2.9132 /
  1.698e-04 on V1; −1.04e-13, −2.45e-02, +7.34e-02, −2.45e-02 on V2; 5.72e-14,
  9.09e-01, 6.31e-01, 3.76e-02, singular, 1.23e-01, **1.19e-15** on V3. The V3
  5-cell blindness is declared on the page.
- **V5 does have the power it claims.** It carries no break row, so I injected
  one: `b_offdiag = +1` moves it 1.192e-13 → 1.200e+01, `scalar=True` → 8.57e-01.
  The claim "it tests eq. (87) and the inversion" is true. It genuinely *cites*
  A4.4 (retypes eq. 109's binary form, which I checked against the printed
  eq. 109 and against A4.4's `Cell.d_scalar` — identical) rather than re-deriving
  A4.4's Bosanquet error. The local-vs-surface Ð^e_ij bug the builder reports is
  real and correctly diagnosed (Ð^e_ij ∝ 1/p, and c_t is not uniform in a
  reacting pellet).
- **The self-criticism checks out.** A4.4 *does* build a binary 2×2 [Bᵉ]
  (`Cell.b_matrix`, lines 595–604) — the corrected claim is the true one. A4.4's
  own published scope note explicitly hands A4.3 "eqs. (86)–(87), its ternary
  He/Ne/Ar comparison ... in Fig. 44, and the viscous-flow term". A4.2 contains
  no pore, no wall, no Knudsen term (only a forward pointer to A4.3). **Scope
  split confirmed**, with one unflagged repeat: A4.4 already runs the Jackson
  ν_B = 1,2,3,4 pellet sweep, and does the 1/Kn analysis of its residual more
  carefully than A4.3 does. A4.3 does say "Page A4.4 reproduced that number".
- **pymrm conventions.** `NumJac((n_cell, n_species), axes_diagonals=[0],
  axes_blocks=[-1])` — never a bare 1-D shape; the n = 1 case in V2 correctly uses
  `(200, 1)`. I compared the NumJac Jacobian against a dense central-difference
  Jacobian: max relative difference 5.4e-09, identical sparsity, **zero** entries
  the FD Jacobian has that NumJac zeroes. `nu` in `construct_div` commented as
  geometry; both BCs commented on the outward normal. Constant operators built
  once in `__init__`.
- **Prose vs output: no drift.** Every interpolated number in every
  `display(Markdown(...))` cell matches the cell output above it. Nothing is
  retyped.
- **Reference block.** `reference:` (K&W, read) + `origin_not_consulted:`
  (Mason & Malinauskas, Jackson, Mason & del Castillo, Knudsen, plus Kaza &
  Jackson, Haynes, Remick & Geankoplis), **no** `reference_read_from` anywhere.
  `models_entry.yaml` key set is identical to the published A4.4 entry; slug and
  title match `meta.yaml`. Tier 6 stated in the title block, "The data",
  "Parameters and assumptions", the blind-spot table and "What pymrm adds".
- **Pellet continuation is not laundering anything.** Re-solving from a cold
  start at the final φ reproduces the continuation answer to 4e-16. The uphill
  threshold, the one number the page claims determinism for, uses `capillary`,
  which has no continuation at all — that claim is true.
- **Figure quarantine: clean.** `git check-ignore -v` puts
  `A4.3-fig44-overlay.png`, `A4.3-fig44-points.csv`,
  `krishna-wesselingh-1997-fig44.csv`, its `.meta.yaml` and `p887-27.png` all on
  `.gitignore:42 queue_cases/*/review/*`; only `README.md` (line 43) and
  `make_overlay.py` (line 44) are negated. **The `.md` → `.csv` change worked**:
  the marker table is the one artefact that would have been tracked had it stayed
  `.md`. Neither tracked file carries an extracted coordinate or image bytes —
  `make_overlay.py` holds only the tick-pixel calibration and the three mask
  rectangles; `README.md` holds only the questions. `page/data/` contains no
  Fig. 44 file and nothing in the notebook loads one (the only mentions are
  prose saying it is *not* used).

---

## Findings, ranked

### F1 — CONFIRMED, send-back. `d₀* = 212.21 nm` is one arbitrary composition's threshold, and the screening rule it supports is false

The page reports, to five significant figures and with a grid-independence of
7.7e-06, that uphill diffusion "survives only above d₀* = 212.21 nm (Kn* = 2.44)
**for H₂/N₂/CO₂**", and concludes: *"In a pellet whose pores are a few nanometres
across — which is most of them — the dusty gas model says it cannot happen."*
The same claim is in `README.md`, `meta.yaml`'s `adds:` and the `description:`
of `models_entry.yaml`.

The threshold is a function of the near-end composition, which the page fixes at
`X_UP = (0, 0.49, 0.51)` without saying that the answer depends on it. Bisecting
the same deterministic function at other compositions **inside the page's own
stated uphill window** `x_N2 ∈ [0.4823, 0.5002]`:

| x_N2 (near end) | d₀* |
|---|---|
| 0.4830 | 4245.7 nm |
| 0.4850 | 923.8 nm |
| 0.4900 | **212.21 nm  ← the page's number** |
| 0.4950 | 63.7 nm |
| 0.4980 | 20.5 nm |
| 0.4995 | 4.67 nm |
| 0.4999 | 0.913 nm |
| 0.49999 | no sign change down to 0.1 nm |

At x_N2 = 0.5000 exactly (far end also 0.5), N_N2 and ∇c_N2 have the **same sign
at every pore diameter from 1 nm to 1 mm** — uphill by the page's own
`N·∇c > 0` criterion, at 1 nm. That is not a numerical artefact: it is
grid-independent to seven digits over 50→800 cells, Newton residual 1e-14, and
|N_N2/N_H2| = 5.9e-05, four orders above roundoff.

*Failure scenario.* A reader designing a nanoporous catalyst reads "a scalar
effective diffusivity is not missing anything on that count", drops the matrix
closure, and misses an uphill N₂ flux that the dusty gas model does predict at
1 nm for compositions the page's own window admits. The number is also now in the
catalogue description, so it propagates to anyone searching `models.yaml`.

This is the B1.4 pattern the handoff warns about: a bisection that converges
beautifully on a function whose *definition* contains an unstated free parameter.
The grid-independence number gives the result an authority it has not earned —
it measures the solver, not the threshold.

*Fix.* Either (a) report the threshold as a function of composition and take the
supremum over the window as the physical statement (in which case uphill survives
below 1 nm and the screening rule inverts), or (b) state plainly that 212.21 nm is
the threshold *at x_N2 = 0.49* and delete the "a few nanometres ⇒ cannot happen"
rule. (a) is the honest version and is still an interesting result — the pore size
below which uphill needs an increasingly fine-tuned composition.

### F2 — CONFIRMED, send-back. The 6.9 % viscous relief is measured at a Knudsen number a 1 µm pore cannot have

The page derives D_visc/Ð^e_iM = 3pd₀/(32ηv̄ᵢ) — **I verified this assembly
against the printed eqs. (85), (91), (106); it is exactly right**, and the
2.2e-16 agreement is a transcription check (an algebraic identity, correctly not
presented as more). It then evaluates the group at d₀ = 1 nm (2.639e-04) and
1 µm (0.2639) and feeds both into a pellet sweep run at **Kn = 1e4**, reporting
that viscous flow changes the pellet-centre pressure by 0.01 % "in the nanoporous
catalyst" and **6.9 %** "in the macroporous support".

But Kn and D_visc/Ð^e_iM are both functions of d₀ and move in *opposite*
directions — the page's own Reuse section says so ("Both are linear or
inverse-linear in d₀, so pore size decides both"). For H₂/N₂ at 1 atm:

| d₀ | D_visc/Ð^e_AM | **actual** Kn | p₀/p no visc | p₀/p visc | change |
|---|---|---|---|---|---|
| 1 nm | 2.639e-04 | 138.9 | 1.409187 | 1.409057 | **−0.009 %** |
| 10 nm | 2.639e-03 | 13.89 | 1.369428 | 1.368308 | −0.082 % |
| 100 nm | 2.639e-02 | 1.389 | 1.193678 | 1.188841 | −0.405 % |
| **1 µm** | 2.639e-01 | **0.1389** | 1.034560 | 1.028110 | **−0.62 %** |
| 10 µm | 2.639e+00 | 0.01389 | 1.003262 | 1.000955 | −0.230 % |
| 100 µm | 2.639e+01 | 0.001389 | 1.000294 | 1.000011 | −0.028 % |

The 1 nm figure is fine (0.009 % either way, because 1 nm really is near-Knudsen).
The 1 µm figure is not: the physically consistent answer is **0.62 %, not 6.9 %**,
an eleven-fold overstatement, and the whole physically reachable range peaks at
0.62 %. The reason is that a macroporous pellet has almost no
reaction-generated pressure rise to bleed off in the first place (p₀/p = 1.035,
not 1.414).

Consequently the page's conclusion is inverted. It says:

> **Anyone modelling a macroporous membrane reactor and dropping the viscous
> term is making a tens-of-percent error in the internal pressure, and the
> review's sentence does not license that.**

Measured consistently, the error is sub-percent at every pore size, i.e. Haynes's
"not very serious" holds *everywhere*, not just for catalysts — the opposite of
the page's headline. The 15.4 % and 29.0 % relief figures quoted at
D_visc/Ð^e = 1 and 100 correspond to Kn ≈ 0.037 and 0.0004, where the pressure
rise is 1.003 and 1.0003 — those percentages are percentages of nothing.

This is the verifier brief's confounding check: two "groups" (1 nm, 1 µm) that do
not overlap in the variable held fixed.

*Caveat, in fairness.* The two groups are only locked together **because the page
uses eq. (91)'s cylindrical-pore B₀ = d₀²/32 with the same d₀ that fixes Ð^e_iM**
and labels the sweep points "1 nm pore" and "1 µm pore". A model taking B₀ from
eqs. (92)–(94) (packed bed) could decouple them. But that is not what this page
does, and it names the pellets.

*Fix.* Replace the fixed-Kn sweep with the consistent one above (it is six pellet
pairs and runs in seconds), keep the screening group — which is correct and
genuinely useful — and rewrite the conclusion to what the consistent sweep says.
The result is still new and still worth printing; it just says the reverse.

### F3 — CONFIRMED, fix required. "eq. (85) is pinned separately by V3 and V2" is false for the prefactor, and nothing on the page pins it

The claim appears four times: page V1 preamble ("eqs. (83) and (85) are pinned
separately, by the symbolic derivation of Graham's law from eq. (85)'s exponent
and **by eq. (103) for the pore-size factor**"), `README.md`, `meta.yaml`'s
`validation:` block, and `A4.3.yaml`'s notes.

Measured, by replacing eq. (85)'s `d₀/3` with `d₀/2` and `d₀/30`:

| route | intact | d₀/2 | d₀/30 |
|---|---|---|---|
| V2 (eq. 103) | −1.04e-13 | **−2.81e-13** | **−1.54e-11** |
| V3 (Graham) | 5.72e-14 | **8.18e-14** | **0.00e+00** |
| V1 (both routes) | 3.97e-07 | 2.75e-06 (page prints "BLIND") | — |

V2 is blind because its "printed closed form" reference is
`N_EXACT = (knudsen_D(...)·Δc + RTB₀/η·Δ(c²)/2)/L` — it calls the same
`knudsen_D` the solver calls, so a wrong prefactor moves both sides together.
V3 is blind because the prefactor cancels out of Σ Nᵢ√Mᵢ = 0 identically; Graham
pins the **mass exponent only** (9.09e-01 when inverted), which is what the
on-page blind-spot *table* correctly says. V4A and V5 both declare themselves
blind to eq. (85). V6 tests linearity in d₀, which d₀/2 preserves.

So: **nothing on this page pins eq. (85)'s d₀/3.** Nothing pins eq. (83) either —
the page's own closing paragraph says "Nothing tests ε/τ, which is assumed",
directly contradicting the V1 preamble's "eqs. (83) and (85) are pinned
separately" two screens above.

*Failure scenario.* Had the builder made the classic Knudsen error — reading d₀
as a radius, or transcribing d₀/2 — V1 would have printed 2.75e-06 and been
labelled "(no defect)", V2 −2.8e-13, V3 8.2e-14, and every headline number
(d₀* = 212.21 nm, Kn* = 2.44, the viscous group, every Knudsen number on the
page) would be wrong by a factor of 1.5 with the validation section reporting
machine precision throughout. I confirmed by independent transcription that this
did **not** happen — eq. (85) really is (ε/τ)(d₀/3)√(8RT/πMᵢ) and d₀ really is the
diameter — but the page's statement about its own coverage is wrong, and it is
wrong in the direction that reads as candour.

*Fix.* Replace the clause with: "eq. (85)'s **mass exponent** is pinned separately
by V3; its d₀/3 prefactor and eq. (83)'s ε/τ are pinned by nothing here and are
transcription-only." Four files.

### F4 — CONFIRMED, fix. "worst relative Newton residual anywhere on this page" is not

`RES_WORST` is a max over V1, V2, the inert runs, the viscous sweep, V5, `jack`
and `surf` — and reports 1.2e-10. It omits `mass`, `knsw`, the uphill capillary
solves, V4, V6 and V7. V4B's own printed output, two cells earlier, shows a
residual of **1.8e-08** at d₀ = 100 mm. The summary line and README both say
"anywhere". 1.8e-08 is still small and the claim's spirit survives, but the word
is false against the page's own output.

*Fix.* Either fold V4 into `RES_WORST` and print 1.8e-08, or say "anywhere in the
runs listed above".

### F5 — CONFIRMED, fix. "it shrinks when either is pushed" — false for φ

The page attributes the −5.0e-05 Jackson residual to "the finite Kn = 10000 and
finite φ = 60 of the run, and it shrinks when either is pushed". Measured:

| Kn | φ | dev from √2 |
|---|---|---|
| 1e4 | 60 | −5.004e-05 |
| 1e5 | 60 | −5.005e-06 |
| 1e6 | 60 | −5.005e-07 |
| 1e4 | 120 | −5.028e-05 |
| 1e4 | 200 | **−5.093e-05** |

Kn: exact 1/Kn, as A4.4 established. φ: the deviation *grows* slightly. The
claim is half wrong and neither half is demonstrated on this page (A4.4
demonstrates the Kn half).

*Fix.* "...it is the finite Kn = 10⁴ of the run: it falls as 1/Kn (5.0e-06 at
Kn = 10⁵, 5.0e-07 at 10⁶), as `A4.4` established."

### F6 — CONFIRMED, minor. "every scalar closure" is true only for finite D > 0

The Results prose states the premise correctly ("N_C = −D_C∇c_C with D_C > 0, so
N_C = 0 forces ∇c_C = 0"), and I verified the diagonal closure gives an exactly
flat inert profile (variation 0.00e+00 at Kn = 0.1, 1, 10 — the page asserts this
without running it; it is true). But the headline bullet, `README.md`,
`meta.yaml`'s `adds:` and "What pymrm adds" all drop the premise and say "every
scalar closure". The review's own scalar closure, eq. (109), does *not* satisfy
it for an inert: with ν_C = 0 the term x_C ν_j/(x_j ν_C) diverges and eq. (109)
gives D_C → 0, so it predicts nothing about the profile rather than predicting a
flat one. *Fix:* add "with a finite positive diffusivity" to the three universal
restatements.

### F7 — PLAUSIBLE, minor. A conclusion is attributed to a work not consulted

"**Haynes's conclusion is a statement about pore size**, not about the dusty gas
model." Haynes (1978) was not consulted, and the page's own caveats say "no
result here is attributed to any of them". The review's sentence is also narrower
than the page implies: it is made about the third RHS term of eq. (109) in
Schnitzlein & Hofmann's (1988) catalytic-reforming calculation (Fig. 55), which
the page does not mention. The "What pymrm adds" version ("*the review's* remark
is a statement about pore size") is fine. *Fix:* use the "What pymrm adds"
wording in Results too — and it needs rewriting anyway under F2.

### F8 — noted, not a defect. V5 is an algebraic identity given the flux ratio

1.192e-13 is machine precision because the steady state enforces N_B = −ν_B N_A
exactly, and eq. (109) *is* the closed form of the 2×2 inversion under that
ratio. The page does not oversell it — it is listed fifth, labelled "Shared: the
diffusivity values", and the injected defects above show it does discriminate. No
change needed; noting it so the next reader is not surprised.

---

## Explicit verdicts requested

1. **V1's blind spot: correctly bounded on the page's blind-spot table (it names
   `knudsen_D` and `pair_matrix`), but the accompanying claim that eq. (85) is
   pinned elsewhere is FALSE for the d₀/3 prefactor.** V2 is measurably blind to
   it (−1.0e-13 → −2.8e-13), V3 is measurably blind to it (5.7e-14 → 8.2e-14) and
   pins only the mass exponent. Nothing on the page pins the prefactor. See F3.
   The transcription is nonetheless correct — I checked it independently.
2. **Singularity claim: CONFIRMED, and provable.** (1,…,1) is an exact left null
   vector of the wall-free [Bᵉ] for any composition and any n ≥ 2; adding the wall
   term makes the matrix strictly column-diagonally dominant, hence nonsingular.
   The `LinAlgError` is robust across grids.
3. **The four new results.** (a) The screening group 3pd₀/(32ηv̄ᵢ) is correctly
   assembled from the printed eqs. (85)/(91)/(106) — but the 6.9 % it is used to
   produce is a fixed-Kn artefact and the conclusion drawn from it inverts under a
   consistent sweep (F2). (b) The uphill threshold is deterministic,
   grid-independent and bisected on a properly bracketed function, but it is
   composition-specific and the screening rule built on it is false (F1).
   (c) Jackson's generalisation p₀/p = 1 + (ν_B√(M_B/M_A) − 1)x_A,s is correct,
   the three conditions are each shown necessary by a sweep that moves the number,
   and it does reduce to A4.4's verified √ν_B; only the "shrinks when φ is
   pushed" rider is wrong (F5). (d) The inert-species result holds, and the
   scalar closure really does give an exactly flat profile; only the word "every"
   needs a qualifier (F6).
4. **Scope split: CONFIRMED** against both published `build_page.py` files —
   A4.2 has no pore, A4.4's own scope note hands this page eqs. (86)–(87), the
   viscous term and Fig. 44, and A4.4's binary [Bᵉ] is acknowledged. One
   undeclared repeat: A4.4 already runs the ν_B = 1–4 Jackson pellet sweep.
   **Figure quarantine: CONFIRMED clean** by `git check-ignore -v` — all five
   extraction artefacts ignored, only `README.md` and `make_overlay.py` tracked,
   neither carrying a coordinate or image bytes, `page/data/` free of Fig. 44,
   and the `.md` → `.csv` change is exactly what keeps the marker table out.

---

## Recommendation

**Send back.** F1 and F2 are the page's two headline new results and both are
wrong as stated — not by rounding, but in their conclusions. F3 is a false
statement about the page's own validation coverage, repeated in four files, and
is the kind that reads as candour. F4–F7 are one-line edits.

Everything else — the transcription, V1, the singularity proof, the break tables,
the numerics, the scope split, the figure quarantine, the reference block, the
tier-6 honesty — is sound, and F1/F2 are both fixable in-place without touching
the solver: the consistent viscous sweep is six pellet pairs, and the uphill
threshold needs a composition axis rather than a single number.

---

*Verifier artefacts (scratch, not committed): independent RK4 solver, defect
injections into `knudsen_D`, dense-FD Jacobian comparison, consistent viscous
sweep, composition sweep of the uphill threshold.*
