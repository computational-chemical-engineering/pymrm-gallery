# E1.2 — adversarial verification

Verifier pass, 2026-08-01. Staged page: `queue_cases/E1.2/page/`.
Source re-read independently on 600 dpi `pdftoppm` renders of
`~/papers/pymrm-gallery/Kunii1968-bubbling-bed-model-IECFund7-481.pdf`, journal pages 481, 482, 490, 491, 492.
Notebook re-executed from the staged file (8.4 s); **every stream output is
byte-identical to the staged outputs**, so there is no execution drift.

Verdict up front: **safe after the listed fixes.** The page's central claim — the
through-flow coefficient 4.5 rebuilt from a flow field that never sees it — is
real and survives every deliberate break I could construct. Two of the nine
reported checks are decoration presented as evidence, one traceability citation
points at the wrong journal page, and the scope decision's third argument is
factually wrong. None of these makes the physics wrong.

---

## Verdict 1 — the attribution

**Adequately sourced *as the page is currently written*, but the sidecar,
the slug and one prose sentence claim more than the source supports.**

What I confirmed myself on 600 dpi renders:

| claim | verdict |
|---|---|
| K&L 1968 does not list Davidson & Harrison in its literature cited | **CONFIRMED.** Journal p. 492, right column: one alphabetical list, `Chu, J. C.` → `Zabrodsky, S. S.`, no entry between Chu and Donnadieu. `Davidson, J. F.` appears once, as second author of `Orcutt, Davidson & Pigford, Chem. Eng. Progr. Symp. Ser. 58 (38), 1 (1958)`. No 1963 book, no Cambridge University Press entry anywhere. |
| K&L defer the derivations to K&L 1968b | **CONFIRMED.** p. 481: "The bubbling bed model, recently proposed by Kunii and Levenspiel (1968b)"; p. 482 left column repeats it. Reference list carries `Kunii, D., Levenspiel, O., Ind. Eng. Chem. Fundamentals, 7, 446 (1968b)`. Not on disk. |
| the attribution rests on the phrase "the Davidson bubble" | **CONFIRMED, and it is thinner than the sidecar says.** `pdftotext` over all 12 pages returns exactly two Davidson hits: the Orcutt reference, and one parenthesis. That parenthesis reads: *"Finally, many simplifying assumptions were made in developing this model (the Davidson bubble, bed with single size of bubbles, emulsion at minimum fluidizing conditions at all gas velocities, etc.)."* It is a list of simplifying assumptions in the closing discussion — **not** a section about the Davidson–Harrison model. |
| it is on journal page 489 | **FALSE.** It is on journal **page 490**, bottom of the left column (PDF page 10). See finding F4. |
| eqs. 9, 10, 11 are printed and unattributed | **CONFIRMED** on p. 482 right column; transcriptions on the page are term-for-term correct, including `4.5`, `5.85`, the `1/2`/`1/4`/`5/4` exponents and the `3u_mf/ε_mf` numerator. |

**Applying the E1.1 test.** E1.1 was parked because K&L "state the relation but
never attribute it, never name it and never test it." E1.2 clears two and a half
of those three:

* *name it* — yes, once, "the Davidson bubble";
* *test it* — yes, decisively. Eqs. 9, 10, 11 and 17 are **evaluated with the
  arithmetic shown and the answer printed** in three worked appendices. This is
  the qualitative difference from E1.1 and it is what makes the reprint carry the
  case: what the page reproduces is the equations K&L themselves print and use.
* *attribute it* — half. They name Davidson; they never name Harrison, never
  cite the 1963 book, and give no citation for eq. 9, eq. 10 or the `0.711`.

So: **build, do not park.** But the page must not be published as "the
Davidson–Harrison bubble model", and — to its credit — the notebook does not do
that. Its title is *"The Davidson bubble: where 0.711's companions come from"*,
and its opening states plainly that "the equations are Kunii & Levenspiel's as
printed; the attribution to Davidson & Harrison is theirs, not this page's".
That is the right framing already. Three things still need correcting:

1. **The slug is `davidson-harrison-bubble`.** The slug is the URL, and it
   asserts a two-author attribution that appears nowhere in any document that was
   read. Recommend `davidson-bubble`, matching the title and matching the only
   phrase K&L actually print. The catalogue key `E1.2` is stable, so the
   catalogue title may keep "Davidson–Harrison bubble" without the page doing so.
2. **"the attribution to Davidson & Harrison is theirs, not this page's"** is not
   quite true. K&L's attribution is to *Davidson*, singular, in one parenthesis
   inside a list of simplifying assumptions on journal page 490. Say exactly that
   — one sentence — rather than letting "theirs" imply a citation.
3. **The `also:` bibliographic detail traces to nothing that was read.** "*Fluidised
   Particles*, Cambridge University Press, Cambridge, 1963" — the year, the title,
   the publisher and the city are all absent from K&L, who give only the surname.
   The equations are not from memory, but this citation is. Mark it as supplied
   from general bibliographic knowledge rather than from a consulted source, or
   drop the publisher and city.

---

## Verdict 2 — the scope decision

`build` is right, but **argument 3 is false and must be corrected.**

| scope argument | verdict |
|---|---|
| 1. catalogue assigns S3 vs S7 | **CONFIRMED.** `docs/catalog-B-reactors.md` line 76 gives E1.2 as S3; line 86 gives E2.1 as S7. |
| 2. E2.1 never solves a flow field | **CONFIRMED.** `pages/E2.1-.../build_page.py` line 314 `ubr = 0.711*np.sqrt(G*db)` and line 319 `fc = 3.0*(umf/emf)/(ubr - umf/emf) + alpha` are one-line inputs to `hydrodynamics()`. E2.1 has no radial grid and no potential solve anywhere. |
| 3. "NO HEADLINE NUMBER IS DUPLICATED … the only overlap is u_br = 42.8 and u_b = 53.9" | **FALSE — CONFIRMED.** |
| 4. E1.2 produces a result E2.1 cannot (the pole) | **CONFIRMED.** Nothing in E2.1 evaluates the cloud outside appendix C. |

On argument 3: **all eight** of E1.2's route-1 reproductions are already published
on E2.1, to every printed digit, from the same CSV:

| quantity | E2.1 (published) | E1.2 (staged) |
|---|---|---|
| A `K_bc` | 46.5552 · 0.12 % | 46.555 · 0.12 % |
| A `ub_minus_u0` | 11.9579 · 0.35 % | 11.958 · 0.35 % |
| B `H_bc` | 0.0358451 · 0.43 % | 0.035845 · 0.43 % |
| B `ub_minus_u0` | 5.73866 · 1.06 % | 5.7387 · 1.06 % |
| C `u_br` | 42.8138 · 0.03 % | 42.814 · 0.03 % |
| C `u_b` | 53.9138 · 0.03 % | 53.914 · 0.03 % |
| C `K_bc` | 5.43494 · 0.09 % | 5.4349 · 0.09 % |
| C `gamma_c` | 0.398154 · 0.46 % | 0.39815 · 0.46 % |

The overlap is 8/8, not 2/8, and both of E1.2's route-1 headline figures — "worst
1.06 %" and "K_bc/H_bc worst 0.43 %" — are E2.1's numbers. The *route* differs for
`K_bc`/`H_bc` (E1.2's leading term comes from the solve, E2.1's is the literal
`4.5`), but since the solve returns 4.499995 the totals are indistinguishable, so
that difference carries no information beyond check 3's direct 4.5 comparison.

This does not change the verdict — the novel content (the flow field, the pole,
the 2-D bubble, the one-field-read-twice result) is genuinely absent from E2.1 —
but the page should say once that route 1 is E2.1's already-published
reproduction reached by a second route, and `scope_decision.answer` argument 3
should be rewritten.

I re-read every CSV row used against the 600 dpi renders myself. All 25 rows of
`kunii_levenspiel_1968_appendix_values.csv` are correct transcriptions.

---

## Verdict 3 — which checks survived the deliberate-break test

I rebuilt `solve_mode`/`davidson_bubble` outside the notebook and injected
defects. Baseline is appendix C's bed, n = 1600.

| injected defect | 4.5 coefficient | R_c/R error | field error | **integral identity** |
|---|---|---|---|---|
| *(baseline)* | 4.499995 | 7.2e-08 | 2.0e-06 | 1.96e-12 |
| ν = 1 everywhere (2-D operator, 3-D truth) | **1.032** (−77 %) | 4.7e-03 | 2.5e-01 | **4.7e-13** |
| ν = 1 in `construct_div` only | **1.078** | 2.7e-02 | 3.5e+00 | 4.2e-01 |
| ν = 1 in the sink term only | **12.07** | 6.1e-02 | 4.9e+00 | 1.0e+00 |
| sink coefficient 6 (ℓ = 2 harmonic) | **0.281** (−94 %) | 4.9e-02 | 6.7e+00 | 6.7e-01 |
| far-field amplitude ν instead of ν+1 | **3.000** (−33 %) | 7.2e-08 | 6.7e+00 | 1.96e-12 |
| outer BC `a = −1` (outward-normal flipped) | **13.507** | 7.5e-05 | 4.0e+01 | **2.0e-12** |
| outer BC `d` sign flipped | **−4.500** | 7.2e-08 | 4.0e+01 | **1.96e-12** |
| bubble-surface solids BC `a = −1` | 4.499995 | 7.2e-08 | 2.0e-06 | 1.96e-12 |
| bubble-surface solids BC → Dirichlet | 4.499995 | *no cloud* | 2.0e-06 | **1.0e-12** |

**Survived — real evidence**

* **The through-flow → 4.5 (1.1e-6).** The page's strongest claim and it holds.
  `4.5` enters the code nowhere: `grep` finds it only in the comparison
  (`abs(coeff - 4.5)/4.5`), in the printed column and in prose. The coefficient
  comes out of `−f′(R) = 3u_f` times the upper-hemisphere integral `πR²`, and it
  moves under **every** operator defect above — 77 %, 68 %, 94 %, 33 %, 200 %, and a
  sign flip. It is correctly independent of the solids field (the two
  bubble-surface breaks leave it untouched), which is the asymmetry the page
  claims. I verified the algebra independently: `f = −u_f(r − R³/r²)`,
  `f′(R) = −3u_f`, `q = ε_mf·3u_f·πR² = 3u_mf πR²`, `q/V_b = 4.5 u_mf/d_b`. This is
  a genuine reproduction of a printed constant from an independent derivation,
  in the `B3.1` mould. **Keep, unqualified.**
* **Grid convergence and the field error (2.0e-6, observed order 2.000).** Moves
  under every defect, and is the *only* metric that catches the outer-boundary
  sign flip (4.0e+01 vs 2.0e-06). The page labels it correctly as testing the
  discretisation and not the physics.
* **The streamline-bisection route (1.2e-14).** The caveat is **on the page and
  correctly worded**, in two places — in the Validation markdown ("it should [agree
  at machine precision]: the separatrix *is* the surface W = 0 … not a second
  estimate carrying its own error") and printed by cell 16 itself. Confirmed.
* **The "no cloud below u_br = u_f" branch** of check 6 (root finder returns `nan`
  at 0.90) is a genuine topology test that can fail.
* **The far-field "no truncation" claim.** Verified: at n = 1600, r_out/R = 2, 5,
  20, 200 give R_c/R = 1.0987060, 1.0987060, 1.0987060, 1.0987058, and the 4.5
  error goes 3.3e-09 → 2.1e-08 → 7.2e-08 → 2.3e-07, i.e. it grows only because
  the same 1600 cells stretch over a larger range. The claim holds exactly as
  stated.
* **The 2-D bubble (3.9e-8)** plus its through-flow, 1.445298 from the ν = 1 field
  against an analytic 1.445299 = 8u_mf/(πd_b). Two routes, genuinely different.
  I re-derived both closed forms independently and they are right.

**Did not survive — decoration presented as evidence**

* **The discrete integral identity (2.0e-12).** See F1. It is structural.
* **"on all three beds"** for the 4.5. See F5. It is one result, not three.
* **R_c/R (7.2e-8) is blind to a common-mode far-field error** — flipping the outer
  `d` sign destroys both potentials (field error 4.0e+01) and leaves R_c/R at
  1.0987060, unchanged to eight digits, because both potentials scale together
  and the root of `W = g′ − f′` is unmoved. Not a defect in the page — the 4.5
  catches it — but worth one clause so nobody treats R_c/R as the guard.

---

## Verdict 4 — the Appendix B pole claim

**CONFIRMED in every part.** Re-read on a 600 dpi render of journal page 490,
right column, bottom:

> **Appendix B** … Bed. *u*_mf = 10 cm./sec., ε_mf = 0.5
> γ_b = 0.001 (estimated)
> *d*_b = 0.5 cm. (estimated)

and corroborated twice on page 491, left column, by the paper's own substituted
arithmetic: `H_bc = 4.5(10)(0.24)(1.18×10⁻³)/0.50 + …` and
`u_b = u_0 − 10 + 0.711(980 × 0.50)^{1/2} = u_0 + 5.8`. Both carry `10` and `0.50`
explicitly. There is no mis-read constant here.

Arithmetic: u_br = 0.711√(980·0.5) = 15.739; u_f = 10/0.5 = 20; ratio 0.787;
eq. 9 gives −14.08; the threshold is d_b = (20/0.711)²/980 = 0.807 cm. All match
the page.

**Appendix B never evaluates eq. 9 — CONFIRMED.** Its worked solution computes
(Re)_t = 20.3, (Nu*)_t = 4.39, H_bc = 3.60×10⁻², δ, (1−ε_f) = 8.70/(u_0+5.8), and
then (Nu)_apparent = (u_0+5.8)/8.70 · [(0.001)(4.39) + (0.806)(0.036)²/(6·6.25×10⁻⁵)(0.0360)].
Only γ_b and H_bc appear; γ_c and V_c/V_b are absent. The paper is not in error,
exactly as the page says. (Appendix A likewise never evaluates eq. 9 — only
appendix C does.)

**The asymmetry claim is sound.** The gas potential's data are φ(R) = 0 and
φ → −u_f r cosθ; u_br appears nowhere in it, and the break test confirms the
coefficient is unmoved when the solids boundary condition is changed outright.
One nuance the page glides over: what is speed-independent is the *volumetric
through-flow* q = 3u_mf πR². Whether q/V_b still means "bubble → cloud"
interchange when there is no cloud is an interpretive step, since the two-phase
series structure bubble → cloud → emulsion has no middle term below the
threshold. One clause would settle it.

**Robustness to ε_mf.** E2.1 already established (and prints) that appendix B's own
printed (1−ε_f)u_b = 8.70 is inconsistent with its stated ε_mf = 0.50, implying
ε_mf ≈ 0.447. That *strengthens* E1.2's claim — u_br/u_f falls from 0.787 to
0.703 — and the conclusion survives for any ε_mf below 0.635. Worth one clause,
because a reader who knows E2.1 will ask.

---

## Findings, ranked

### F1 — CONFIRMED, severity 1. The discrete integral identity claims power it does not have

`build_page.py` lines 688–692 (notebook cell 7) state:

> *"Evaluating both sides on the DISCRETE solution tests the finite-volume
> weights, the geometry index nu and the boundary gradients together: a wrong nu
> or a mis-signed boundary flux breaks it immediately."*

and `meta.yaml` validation bullet 8 repeats it. **Both halves of that sentence are
false**, verified by deliberate break:

* solving the 3-D bubble with ν = 1 everywhere makes the through-flow coefficient
  1.032 instead of 4.500 and the field error 2.5e-01 instead of 2.0e-06 — and the
  identity reads **4.7e-13, better than baseline**;
* flipping the outward-normal sign at the outer boundary (`{"a": -1}`) makes the
  field error 4.0e+01 and the coefficient 13.507 — and the identity reads
  **2.0e-12, unchanged**;
* flipping the sign of the outer `d` makes the coefficient −4.500 — identity
  **1.96e-12, unchanged**;
* replacing the bubble-surface Neumann condition with Dirichlet destroys the
  cloud entirely — identity **1.0e-12**.

The reason is structural: `lhs − rhs` is the volume-weighted sum of the discrete
residuals, which telescopes for *any* consistent (`construct_div`, cell-volume,
ν) triple and *any* boundary data. Its one genuine power is detecting a
**mismatch** between the ν inside `construct_div` and the ν inside the sink term
(0.42 and 1.00 in the two mixed cases above).

*Failure scenario.* This page advertises itself in Reuse as the `S3` skeleton to
copy. A future agent copies it, flips a boundary sign — the single most common
pymrm error, called out at the top of `AGENTS.md` — sees `2e-12` on the
conservation line, and ships a field that is wrong by a factor of three.

*Fix.* Relabel, as `A4.2` now does: *"an algebraic identity — the volume-weighted
sum of the discrete residuals, which closes for any consistent operator and any
boundary data. It detects a mismatch between the geometry index in
`construct_div` and the one in the sink term. It cannot detect a wrong-but-
consistent ν, and it cannot detect a mis-signed boundary condition; the
through-flow check below catches both."* Update `meta.yaml` bullet 8 to match.

### F2 — CONFIRMED, severity 2. Scope argument 3 is wrong; route 1 duplicates E2.1

Table above. 8/8, not 2/8, with identical digits. Rewrite
`queue_cases/E1.2.yaml` `scope_decision.answer` argument 3, and add one sentence
to the page's Validation section saying route 1 is E2.1's published reproduction
reached by a second route for `K_bc`/`H_bc`.

*Failure scenario.* A reader counts "1.06 % over the appendices" as independent
corroboration of two pages when it is one transcription checked once.

### F3 — CONFIRMED, severity 2. Wrong page citation inside the traceability table

The Data section's table says

> `| `C u_br`, `C u_b` | u_br = 0.711(980×3.7)^{1/2} = 42.8; u_b = 13.2−2.1+42.8 = 53.9 | 490, right column |`

Those two lines are on journal page **491, left column** (bottom, under "from
Equation 3"). Journal page 490's right column is **appendix A**, whose
corresponding line is `u_b = u_0 − 1.21 + 0.711(980×0.35)^{1/2} = u_0 + 12.0`. The
other five rows of that table are correct as cited. This is a traceability error
inside the table that exists to provide traceability.

### F4 — CONFIRMED, severity 2. "Journal page 489" is wrong; it is 490

The phrase "the Davidson bubble" is on journal page **490** (PDF page 10), bottom
of the left column. `pdftotext` per page confirms zero Davidson hits on page 489.
The wrong number appears in `queue_cases/E1.2.yaml` `notes:` and in the header
comment of `queue_cases/E1.2/models_entry.yaml`. It is **not** on the rendered
page, so this is metadata-only — but it is the single citation on which the whole
attribution rests, and it would send the next agent to the wrong page.

### F5 — CONFIRMED, severity 2. "On all three beds" is one result, not three

Cell 17 prints, for beds A, B and C, coefficients 4.499995 / 4.499995 / 4.499995
and deviations 1.11e-06 / 1.11e-06 / 1.11e-06 — identical to every digit. They
must be: ε_mf, u_mf and d_b all cancel algebraically
(coeff = ε_mf · 3(u_mf/ε_mf) · πR² / (4/3 πR³) · d_b/u_mf ≡ 4.5), and the grid is
geometric in r/R so the discretisation error is scale-invariant too. The prose —
"it comes out of the solve to 1.1e-06 relative on all three beds", repeated in
`meta.yaml` and `models_entry.yaml` — reads as triple confirmation.

*Fix.* Say what it is: *"identical on all three beds, because ε_mf, u_mf and d_b
cancel — the through-flow coefficient is a pure number, which is itself the
point."* That is a stronger sentence than the current one, and true.

Related, same cell: the reported `1.1e-6` is **discretisation error**, not the
strength of the agreement with the printed constant. The check's real content is
binary — 4.5 versus anything else — and it would fail loudly, by tens of per
cent, if the transcription or the derivation were wrong. Worth one clause so the
`1e-6` is not read as evidence strength.

### F6 — CONFIRMED, severity 2. The opening paragraph overclaims by one result, and is unsourced

> *"Davidson and Harrison obtained all three in 1963 … This page solves that flow
> problem numerically and recovers the three results."*

The page recovers **two**. Three paragraphs later it says so itself: "0.711 —
**taken as given.** It is the Davies–Taylor spherical-cap rise velocity and does
not follow from the percolation problem at all." And the historical claim that
Davidson and Harrison obtained the rise velocity is sourced to nothing that was
read: K&L state eq. 3 bare on p. 482 with no citation for `0.711`, and credit the
model to "Kunii and Levenspiel (1968b)". Change to "recovers two of the three,
and takes the rise velocity as given", and drop or hedge the 1963 provenance of
`0.711`.

### F7 — severity 3. Attribution wording, slug and `also:` detail

Three items, listed under Verdict 1. The slug change is the one that matters,
because it is the URL.

### F8 — severity 3. "with no warning" is slightly unfair to the paper

"What pymrm adds" says a reader carrying eq. 9 into such a bed "gets a nonsense
number **with no warning**". K&L do preface eq. 9 with *"In beds with fast rising
bubbles, the volume of cloud surrounding each of the bubbles is given by"* —
verified on the 600 dpi render of p. 482. The page quotes that qualifier
correctly in "The published model", so it is not hidden, but the two sentences
sit oddly together. Suggest "with no *quantitative* warning — the paper says
'in beds with fast rising bubbles' and never says how fast is fast". The finding
survives intact; it is the phrasing that overreaches.

### F9 — severity 3. Appendix B's ε_mf inconsistency, already published on E2.1, goes unmentioned

See Verdict 4. Adding it strengthens the claim rather than weakening it, and it
pre-empts the obvious objection.

### F10 — severity 3. Reuse and `related:` name pages that do not exist

`B1.6`, `E1.3` and `E1.4` are in neither `models.yaml` nor `pages/`. The Reuse
section names `B1.6` (Prater relation) alongside the published `B1.1` and `J1.5`
as though it were a sibling page, and `related:` lists all three. No CI break —
`scripts/check_metadata.py` does not validate `related` — but the prose promises
pages a reader will not find. `E1.1` exists in `models.yaml` but is parked with
no page; naming it as "the earlier picture this refined" is fine as prose.

### F11 — severity 3. `meta.yaml` bullet 1 misdescribes what the 7.2e-8 measures

> *"the cloud volume computed from the pymrm flow field reproduces Kunii and
> Levenspiel's equation 9 to 8e-8 in R_c/R at n = 1600"*

The 7.2e-08 is measured against `Rc_exact`, the closed form derived on the page,
and is pure discretisation error. What actually shows eq. 9 was read correctly is
(a) the derived `(R_c/R)³ = (u_br + 2u_f)/(u_br − u_f)` matching the printed eq. 9
term for term, and (b) γ_c = 0.39815 against the paper's printed 0.40 (0.46 %).
The Validation markdown gets this right; the meta bullet does not. Also, "8e-8"
here versus "7.2e-8" in `agreement:` two blocks below — pick one.

---

## Things I attacked and found clean

* **Prose against output.** Every number in every markdown cell checked against
  the executed outputs: 1.099, 1.103, 0.787, −14.1, 0.81 cm, 2.55, 10 % thick,
  1.1e-6, 7.2e-8, 2.000, 5.3e-7, 2.0e-12, 3.9e-8, 1.2e-14, 1.06 %, 0.43 %, 0.46 %.
  All match. No drift, and the notebook re-executes to byte-identical output.
* **Is `4.5` an input?** No. `grep` over `build_page.py`: it appears in the
  comparison, the printed column, the LaTeX of eqs. 10/11, and prose. Never in
  the solve.
* **Is `q/V_b` an input?** No. It is computed from `-df[0]`, which comes from
  `Grad @ h + grad_bc`.
* **Is the closed form circular?** `Rc_exact` is typed in, but the numerical solve
  that reproduces it to 7e-8 contains no trace of it, and I re-derived it
  independently: W = 0 gives r³/R³ = (u_br + 2u_f)/(u_br − u_f), hence
  V_c/V_b = 3u_f/(u_br − u_f), which is eq. 9 term for term. The 2-D counterpart
  likewise: r²/R² = (u_br + u_f)/(u_br − u_f). Both correct.
* **ν and boundary conditions.** `construct_div(..., nu=nu)` carries the comment
  "nu: 2 spherical, 1 cylindrical"; both bubble-surface conditions carry the
  physical equation and an explicit note that the outward normal points along
  −r̂ there. Correct, and correctly noted that the solids condition is
  sign-blind because d = 0.
* **Deviation convention.** `100·|computed − printed| / printed` everywhere,
  unsigned, consistent across cells 18 and 19 and consistent with E2.1. Never
  stated on the page — worth one line, but no direction ambiguity arises.
* **Tier 6 honesty.** Stated in three places, and the Validation section calls the
  appendix numbers "the authors' own arithmetic". Nothing is called experimental.
  "What pymrm adds" ends with an explicit negative ("This page cannot say whether
  any of it is *true*"). This is well done.
* **Cross-page data load.** `load_data(..., page="E2.1-kunii-levenspiel-bubbling-bed")`
  resolves to `<repo>/pages/E2.1-…/data/` locally and to the raw GitHub URL on
  Colab. Works, and the page states the provenance and the dependency openly.
  `datasets: []` is correct.
* **Structure.** All nine required sections present and in order; no Quarto
  callouts or shortcodes; Colab cell is cell 1; `report_agreement` called.
* **The figures.** Both re-rendered and inspected. The three-panel streamplot is
  physically right — downward far field with closed recirculation for the two
  fast bubbles, uniformly upward with gas passing straight through for
  u_br/u_f = 0.55. The pole figure's axis label and legend are correct.

---

## Recommendation

**Safe to publish after these fixes:**

1. **F1** — relabel the discrete integral identity in the cell-7 comment/output
   and in `meta.yaml` bullet 8. *Blocking.*
2. **F3** — correct `C u_br`/`C u_b` to "491, left column" in the Data table.
3. **F4** — correct "journal page 489" → "490" in `E1.2.yaml` and
   `models_entry.yaml`.
4. **F5** — rewrite the "all three beds" sentence to say why they are identical.
5. **F2** — rewrite `scope_decision` argument 3; add one sentence to Validation
   noting route 1 is E2.1's published reproduction by a second route.
6. **F6** — "recovers two of the three"; hedge the 1963 provenance of 0.711.
7. **F7** — slug `davidson-harrison-bubble` → `davidson-bubble`; correct the
   "attribution … is theirs" sentence; flag the `also:` bibliographic detail as
   unsourced.

**Recommended but not blocking:** F8, F9, F10, F11, and one clause noting that
R_c/R is blind to a common-mode far-field error while the 4.5 is not.

Nothing found requires sending the case back. The physics is right, the
transcriptions are right, the headline result is real and survives breaking, and
the page's own caveats are unusually honest — the defects are in how two of the
nine checks are described and in four citations.
