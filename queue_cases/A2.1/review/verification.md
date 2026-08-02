# A2.1 — adversarial verification, 2026-08-02

Scope: the staged page `queue_cases/A2.1/page/`, its `meta.yaml`,
`models_entry.yaml`, `queue_cases/A2.1.yaml` and `queue_cases/A2.2.yaml`, against
both source PDFs read at 600 dpi.

**Verdict: safe to publish after the five fixes listed under "Required fixes".**
Nothing on the page is fabricated, no transcription is wrong, the boundary
conditions are encoded correctly, and the headline claims survive. The findings
are all of the "a check is weaker or more grid-conditional than the page says"
class, plus one source sentence the page does not report.

---

## Verdicts on the four questions asked

1. **Covering `A2.2` is right.** Confirmed on the page image.
2. **The boundary-condition encoding is correct.** Confirmed three independent ways.
3. **The undetectable-outlet blind spot is real, and stronger than the page says.**
4. **The harmonic-mean rule is right; the factor 347 is a grid artefact.** This
   needs a wording fix in `AGENTS.md` as well as on the page.

---

## Required fixes

### F1 — MEDIUM. "347× worse" is a property of the grid, not of the arithmetic mean

Re-running `three_section(..., face_average=...)` with only `n_b` varied
(`stretch=False`, `P=20`, everything else as the page):

| `n_b` | arithmetic, fore vs eq. 19 | harmonic, fore vs eq. 19 | ratio |
|---|---|---|---|
| 100 | 1.241e-03 | 1.471e-05 | 84 |
| 200 | 6.248e-04 | 3.632e-06 | 172 |
| 400 | 3.134e-04 | 9.026e-07 | **347** (the page's grid) |
| 800 | 1.569e-04 | 2.254e-07 | 696 |

The arithmetic mean converges at **first** order at the jump, the harmonic mean
at **second**. The ratio therefore grows in proportion to `n` without bound, and
347 is simply its value at `n_b = 400`.

`AGENTS.md` now states, as a measured constant, "the arithmetic mean was **347×**
worse in the fore section". *Failure scenario:* an agent follows that rule on a
50-cell multi-domain model, measures ~40× rather than 347×, and concludes either
that the number was wrong or that the rule does not apply to their problem. The
grid-invariant statement is the order.

Everything else about the harmonic-mean finding **is** confirmed:

- the arithmetic mean **fails silently** — with it, `f(0+)` in the first bed cell
  is 0.667011 against 0.667012, the exit value moves by 4e-09, and the profile is
  smooth and plausible;
- the damage is localised to the fore section (bed error 1.6e-06 against
  8.5e-07, exit unchanged), which the page already scopes correctly.

*Fix:* on the page and in `AGENTS.md`, say "the arithmetic mean drops the scheme
from second order to first at the jump — 347× worse in the fore section at this
page's `n_b = 400`, and growing in proportion to `n`."

### F2 — MEDIUM. The "after section flat" residual is structural and is not labelled

The page prints, in a block of five deviations:

```
after section  flat?         : peak-to-peak    2.64e-11
```

and `meta.yaml` lists "the after section flat to 2.6e-11" as a validation bullet.
The line immediately below it correctly labels the *fore-section flux* as
structural, but this one is left to read as a measurement.

It is structural. With no reaction, constant `D` and the zero-gradient outlet,
flux constancy propagates backwards cell by cell and forces `c_i = c_{i+1}`
exactly, for **any** `Pe_c` and **any** truncation length. Measured, truncating
the after section from 30 decay lengths to 0.2:

| after section | cells | peak-to-peak |
|---|---|---|
| `P = 30` (as the page) | 102 | 2.25e-15 |
| `P = 2` | 56 | 2.00e-15 |
| `P = 0.2` (absurdly short) | 22 | 2.61e-15 |

*Failure scenario:* a reader, or a CI baseline, treats 2.6e-11 as evidence that
the after-section physics or the truncation length is right. It would read the
same with `Pe_c` wrong by six orders of magnitude or the domain 150× too short.
This is the exact defect class of `handoff.md#the-check-that-cannot-fail`.

*Fix:* one clause — "structural, like the fore-section flux below: the discrete
equations force it for any `Pe_c` and any truncation."

### F3 — MEDIUM. The recovered `df/dz(1)` is mostly the diagnostic's own fit bias, and the identity behind it is exact

Two measurements.

(a) Apply the page's own 3-point quadratic fit to the **exact** eq. 33 samples on
the same cell centres — where the true gradient at `z = 1` is exactly zero:

| `n_b` | page's recovered `df/dz(1)` | same fit on exact eq. 33 | genuine residual |
|---|---|---|---|
| 100 | −3.640e-04 | −3.044e-04 | −5.96e-05 |
| 400 | −2.340e-05 | −1.973e-05 | −3.67e-06 |

84 % of the printed number is the extrapolation error of the diagnostic, not the
solution's gradient. (Both parts do converge at second order, so the *conclusion*
stands; the magnitude does not mean what it appears to.)

(b) In the discrete system the last bed cell equals the after-section value to
1.1e-16 (measured at `n_b` = 100 and 400). So the chain "after section flat ⇒
zero diffusive flux at `z = 1`" holds at machine precision, and the effective
discrete condition at `z = 1` in the three-section solve is **identical** to
`bc_out = {a:1, b:0, d:0}` in the closed vessel.

The page's conceptual caveat is already honest ("not an independent proof of its
premise"). What is missing is that the identity is exact and the reported
second-order decay is a post-processing artefact.

*Fix:* print the same fit applied to the exact closed form next to it, or state
that the residual is the fit's truncation error.

### F4 — LOW/MEDIUM (honesty). Wehner & Wilhelm qualify the inlet discontinuity and the page does not report it

Two sentences on the page image:

- p. 90, of Danckwerts' eq. (2): *"There is a discontinuity at this boundary as
  reported which will be shown later to be **not necessarily correct**."*
- p. 91: *"The boundary condition as written by DANCKWERTS is the equivalent of
  the step function mentioned above since he did not include the term
  (1/Pe_a)·df(z<0)/dz."*

The page's opening sentence — "a **flux** balance at the inlet, which makes the
concentration drop discontinuously across the bed entrance" — and its first
figure (a dotted vertical drop from 1.0 to 0.668 at `z = 0`) present the step as
the physics. In Wehner & Wilhelm's own three-section model, `f` is **continuous**
at `z = 0`; the fall from 1 to `f(0)` happens over the fore section with decay
length `1/Pe_a`. Measured, adjacent-cell difference across `z = 0` at `n_b = 200`:

| `Pe_a` | 0.3 | 1 | 3 | 30 | 1000 |
|---|---|---|---|---|---|
| `f(0−) − f(0+)` | 2.4e-03 | 3.0e-03 | 4.6e-03 | 2.6e-02 | 2.7e-01 |

i.e. a kink in slope resolved by the mesh, not a jump — the discontinuity is the
`Pe_a → ∞` limit. The page's three-section figure shows this correctly, but no
prose reconciles the two pictures, and the paper's own qualifying sentence is one
the page's framing argues against. Verifier check 9.

*Fix:* one sentence in the Background or in check 3 — the step at `z = 0` is the
closed-vessel idealisation, exact only as `Pe_a → ∞`, and Wehner & Wilhelm say so.

### F5 — LOW. The page description says "eight decades" where the naive-inlet sweep is 4.2

Notebook front matter: *"…the naive inlet's error measured across **eight
decades** of Peclet number."* The naive-inlet sweep is
`np.logspace(-1.5, 2.7, 43)`, i.e. Pe = 0.0316–501, **4.2 decades**. Eight decades
belongs to check 2 (pymrm against eq. 34). `models_entry.yaml` and `meta.yaml`
both already say "four decades" correctly; only the notebook `description:` — the
gallery card subtitle — is wrong.

---

## Optional, not blocking

- **F6.** "This is the check with the most power, because it obtains the same
  profile from a **completely different problem statement**" oversells the
  three-section route. It shares `construct_grad`, `construct_div`,
  `construct_convflux_upwind`, the van Leer deferred correction and the solver
  with `solve_closed`, and by F3(b) its `z = 1` condition is discretely identical.
  What it genuinely establishes is that the Danckwerts *closure* is the correct
  reduction of the three-section problem — it separates Danckwerts from Hulburt by
  50 % at the reference point — plus invariance to `Pe_a` and `Pe_c`. Suggest
  "different problem statement, same operator library".
- **F7.** `queue_cases/A2.2.yaml` contradicts itself: "Nothing in this paper is
  left unused" and, two lines later, "the only material not already on A2.1 is the
  non-steady-state discussion in their Appendix". Both also miss W&W eqs. (23)–(26)
  and the p. 92 paragraph reconciling `df/dz(1) = 0` with the plug-flow gradient
  `−R·exp(−R)` (*"The key term in this development is {1 − exp(a·Pe_b(z−1))}"*) —
  a boundary-layer / non-uniform-limit argument. Neither omission changes the
  `covered` verdict; both would be natural additions to the page.
- **F8.** The TVD deferred-correction loop returns silently after `max_it = 80`.
  Every solve on the page converges (worst 26 iterations, at Pe = 1e4) and the
  notebook re-executes bit-identically, so there is no live defect — but
  `handoff.md`'s `B1.6` lesson is to assert it.
- **F9 (cosmetic).** The eq. 37 decay exponent is printed as 1.95 from a fit over
  Pe = 10–10^4.5; the asymptotic exponent is 2.00 (local slope at the top of the
  range 1.9993). The page's conclusion is right and the number understates it.
  Also `f_exact_ref :=` in the check-5 cell is an unused assignment.

---

## What I attacked and could not break

### Transcriptions — every symbol, on my own 600 dpi renders

`pdftoppm -r 600` of Danckwerts p. 10 and Wehner & Wilhelm pp. 89–93, read as
images. All of the following match the page exactly:

- Danckwerts eq. (30) `d²c/dy² − (u/D)dc/dy − kc/D = 0`; eq. (31)
  `uc* = uc − D·dc/dy, y = 0`; the unnumbered outlet balance `ufc* = uc − D·dc/dy`;
  eq. (32) `dc/dy = 0, y = L`.
- Danckwerts eq. (33) — including the **leading 2 on both** numerator terms and
  `a = √(1 + 4kD/u²)`.
- Danckwerts eq. (34) — including `−uL/2D` on **both** exponentials.
- Danckwerts eqs. (35) `1 − exp(−kL/u)`, (36) `kL/(u + kL)`, (37)
  `1 − (1 + k²DL/u³)exp(−kL/u)`, and the `k²DL/u³ ≪ 1` / `[ln f]²D/Lu ≪ 1` sentence.
- W&W eq. (1), eqs. (2)/(3), eqs. (8)–(15), eqs. (16)–(18), eqs. (19), (20), (21),
  (22), `a = √(1 + 4R/Pe_b)`, `g₀`, and the Notation page (`Pe = Lu/D`, `R = kL/u`).

All quoted prose is verbatim, including the two long Danckwerts quotations and
the "Intuition suggests…" sentence.

### The decisive scope facts

- **The p. 91 sentence exists**, read on the render with the folio "91" in the
  same crop: *"The solution for the reaction section (equation 20) is identical
  with that of DANCKWERTS."*
- **The two equations really are one formula.** Transcribed independently, eq. (33)
  is eq. (20) with the factor 2 moved out of `g₀`; eq. (34) is eq. (21) with
  `exp(Pe/2)` distributed. The 1.1e-16 / 5.6e-17 agreements are correct and the
  page labels them for what they are ("a check on the reading, not on the
  mathematics"), with a measured sensitivity (2.9e-02 for one wrong digit).
- **The "(13)" reading is right.** The text layer gives "(18)" — the 3→8 trap.
  Eq. (13) is the only reading that is a *boundary condition* and the only one
  from which "using also equation (14)" yields `df(1−)/dz = 0`.
- **89–93 confirmed.** The header glyph is unreadable even at 600 dpi, but the PDF
  holds exactly five pages and the printed folios "91" and "93" are legible. The
  catalogue correction is right; `docs/papers-on-disk.yaml` already carries it.
- **Fig. 4 parameters.** "Pe_b = 2·667", "Pe_b = ∞", "Pe_b = 0", "R = 2" read as
  printed labels; `Pe_a = 1` is stated in the **body text** on p. 92, not only in
  the figure, so the page is if anything conservative. No curve digitised anywhere.

### The closed form actually solves the stated problem

The page says its transcription check cannot test this. An independent
`scipy.solve_bvp` on `f'' = Pe(f' + Rf)` with `ya[0] − D·ya[1] − 1 = 0` and
`yb[1] = 0` — no shared code with either the closed form or pymrm:

| Pe | bvp `f(1)` vs eq. 34 | bvp vs eq. 33 over the bed |
|---|---|---|
| 0.1 | 2.2e-16 | 1.1e-16 |
| 1 | 2.2e-16 | 1.1e-16 |
| 2.667 | 2.2e-16 | 2.2e-16 |
| 10 | 2.2e-16 | 5.6e-16 |
| 200 | 7.3e-15 | 6.0e-15 |

### Boundary-condition encoding

Inlet `{a: D/u, b: 1, d: c*}` with `n = −z`, outlet `{a: 1, b: 0, d: 0}` with
`n = +z`. Correct — verified by the order-2.00 convergence to eq. 34, by the 197 %
error when the inlet sign is flipped, and by the independent BVP above, which
writes the same physical condition in a completely different form.

### The outlet blind spot — stronger than the page claims

Not merely "the deviation does not move". With `b = 0` and `d = 0`, the
**assembled matrix and right-hand side are bit-identical**:

| outlet `a` | max &#124;ΔA&#124; vs `a = 1` | max &#124;Δb&#124; |
|---|---|---|
| −1 | 0.0 | 0.0 |
| 1e7 | 0.0 | 0.0 |

(The inlet, for contrast: max &#124;ΔA&#124; = 1.8e+03 for a flipped sign.) The
claim is correctly *scoped* by its because-clause: with `d = 0.05` the outlet sign
does matter (`f(exit)` 0.24700 vs 0.22202). Prominence is adequate — it appears in
the BC section, the break table, "What pymrm adds" point 2, `meta.yaml` caveats
and the README table, and it is stated as a limitation, not buried.

### The quantified Dirichlet error

- Hulburt = Danckwerts / `f(0)` proved analytically: `f_D(z)/f_D(0)` satisfies the
  same linear ODE, equals 1 at `z = 0`, and has zero gradient at `z = 1`. It is
  also W&W's own printed sentence (p. 92).
- Table reproduced by hand: Pe = 0.1 → +179.1 %, Pe = 2.667 → +49.7 %,
  Pe = 200 → +1.0 %.
- Pe → 0: `f(0) → 1/(1+R)` and `f_exit → 1/(1+R)`, so the naive ratio → **1**
  exactly, against the stirred tank's 0.3333. Confirmed analytically, not just
  numerically.
- The pymrm Dirichlet solve is a genuinely independent route to the same number
  (7.1e-05 worst over Pe = 0.0316–501 at n = 1200) — this check *can* fail.

### The numerical claims

- Order study reproduced: upwind 1.00, van Leer 2.00, 2.2e-07 at n = 1600.
- **The numerical-diffusion claim holds and is worth keeping.** Comparing bare
  upwind against eq. 34 at `Pe_eff = 1/(1/Pe + 1/2n)` leaves a residual falling at
  order 1.97 → 2.00 over five refinements. `Pe_eff = 133` at Pe = 200, n = 200 ✓.
- Pe sweep: worst 1.877e-04, located at Pe = 1e4 as `meta.yaml` says, where the
  cell Peclet number is 6.25. Every TVD iteration converges (max 26).
- Fore-section flux value 1.0 as a truncation check: correct. Analytically the
  fore-section flux equals `N₁ = 1 − N₂·exp(−Pe_a·Z)`, and the measured 7.0e-10
  tracks `exp(−20) = 2e-09`. The "structural in form, informative in value" label
  is honest.
- Break table: all nine defects behave as printed; each moves its own target; only
  the outlet sign does not move, as advertised; the baseline row is recovered.

### Routine

- Notebook re-executes in 4.4 s with **byte-identical** outputs (`runtime_seconds: 5` ✓).
- Every number in the markdown matches the cells (50 %, 179 %, 197 %, 133, 0.6 %,
  1 %), and every number in `README.md` matches too. Only F5 is wrong.
- Nine required sections present, in the AGENTS.md order.
- Tier 6 stated plainly; nothing called experimental; the Raschig-ring figure is
  correctly excluded with Danckwerts' own sentence as the reason.
- Deviation convention `(model − reference)/reference` used throughout.
- Shape `(n, 1)` everywhere; no `NumJac` (the problem is linear, and the page says
  so); `nu = 0` commented at both `construct_div` calls.
- `covers:` follows the existing `A1.1` / `B1.1` precedent; `A2.1` is absent from
  `models.yaml` as the integrator note states; `A2.3.related` already names `A2.1`.
- `check_agreement.py` skips a page with no committed baseline, so the new
  `agreement.json` integrates cleanly.
