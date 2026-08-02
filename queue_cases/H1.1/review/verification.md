# H1.1 — verifier report

Adversarial verification of the staged page `queue_cases/H1.1/page/` against
Itoh, *A membrane reactor using palladium*, AIChE J **33**(9) 1576–1578 (1987).

Everything below was re-read from the paper independently of the builder's
`data/*.meta.yaml`, independently of `queue_cases/H1.4/review/`, and
independently of the PDF text layer. **The scan is a 300 dpi bilevel CCITT
image** (`pdfimages -list`: `2400×3239 gray 1 1 ccitt`), so a `pdftoppm -r 600`
render is a 2× upsample carrying no extra information; every glyph and every
drawing feature quoted below was read on the *native* bitmap (`pdfimages -png`),
printed as ASCII ink maps.

Every computation reported here was done with code I wrote, not the page's,
except where the page's own solver is explicitly driven as an object under test.

**Verdict: safe to publish after the fixes listed in F1–F3 (F4–F7 are cheap and
should ride along).** The page is unusually honest — it breaks eleven checks on
purpose and publishes the sensitivity tables — and its arithmetic is right
everywhere I could test it. But it makes one false mathematical claim, it argues
its headline in a way that *is* circular when it does not have to be, and two of
its three new "Itoh's own tube" numbers describe a different problem from the one
the summary lines attach them to.

---

## Findings, ranked

### F1 — MAJOR, CONFIRMED. The free-boundary "limitation" is a Newton globalisation failure described as a mathematical impossibility

The page states, in five places (notebook §4, "What pymrm adds", "Reuse" trap 3,
`meta.yaml` `caveats:`, `README.md`, and again in `queue_cases/H1.1.yaml`
`notes:`):

> the fixed-grid Newton solve has no non-negative solution at the transition cell

and

> the first-order upwind system has no non-negative solution at the transition
> cell … a free-boundary problem needs a free-boundary method

**Both claims are false.** The discrete system the page assembles,

```
(w_i − w_{i−1})/h + Π (max(w_i,0)/(max(w_i,0)+β))^{1/2} = 0 ,
```

has a unique non-negative solution at *every* cell and for *every* Π. The
cell residual `g(w)` is continuous and strictly increasing on `w ≥ 0`, with
`g(0) = −w_{i−1}/h < 0` and `g(w_{i−1}) = Π√y > 0`, so a root always exists in
`(0, w_{i−1})`. Downstream of the front the iteration `w_i ≈ w_{i−1}²/(hΠ/√β)²`
decays *quadratically per cell*, underflows to exactly 0 within a few cells, and
`w = 0` then satisfies the cell equation identically.

Solved cell-by-cell with a bracketed scalar root find (same equations, same grid,
no pymrm):

| Π/Π_crit | n_z | all w ≥ 0 | cells exactly 0 | max residual | max\|w − w_exact\| |
|---|---|---|---|---|---|
| 1.05 | 800 | yes | 28 | 4.4e-14 | 1.13e-03 |
| 1.50 | 800 | yes | 256 | 4.5e-14 | 1.62e-03 |
| 3.00 | 800 | yes | 523 | 1.3e-13 | 3.23e-03 |
| 10.0 | 800 | yes | 711 | 4.4e-14 | 1.07e-02 |
| **18.1** (Itoh's own point) | 800 | yes | 747 | 4.6e-14 | 1.91e-02 |
| 18.1 | 6400 | yes | 6036 | 6.3e-13 | 2.44e-03 |

First order, as expected, and the front is placed to O(h) — the fixed grid does
*not* fail to place it.

**And the page's own solver finds these solutions.** Seeding
`Permeator.w` with the marched profile and calling the page's unmodified
`solve()`:

| Π/Π_crit | flat guess (what the page does) | seeded with the march |
|---|---|---|
| 0.95 | converged | converged |
| 1.05 | **FAILED**, min w = −2.05e-06 | converged, min w = 0, ‖r‖∞ = 1.3e-13 |
| 1.50 | **FAILED**, min w = −1.56e-07 | converged |
| 3.00 | **FAILED**, min w = −1.31e-06 | converged, min w = 0, ‖r‖∞ = 2.0e-13 |
| 18.1 | **FAILED**, min w = −1.70e-07 | converged, min w = 0, ‖r‖∞ = 7.8e-14 |

So the page's `tol = 1e-12` Newton converges over the entire sweep range,
including Itoh's operating point, and returns a certified non-negative solution.
What actually fails is the **flat initial guess plus a `max(w,0)` clip**: the
sink's Jacobian diverges like `w^{−1/2}` as `w → 0`, `NumJac`'s perturbation is
an *absolute* floor (`eps_jac = 1e-6`, see `pymrm/src/pymrm/numjac.py`) so the
differenced derivative is meaningless for the sub-1e-6 values the solution
takes, and the line search parks the iterate at a slightly negative `w` where
the clipped residual is non-smooth. That is a globalisation problem with a
standard fix (march, continuation in Π, or a complementarity formulation), not a
property of the half-power law.

*Failure scenario.* The "Reuse" section tells the next builder that *"a sink with
an exponent below 1 drives its variable to zero in finite space, and a fixed-grid
Newton solve has no non-negative solution past that point"*. A builder who
believes that will reach for front tracking (`S12`) on a problem a marched
initial guess solves in one line. This is precisely the class `handoff.md` flags
— *"a code comment that asserted a sensitivity the check did not have"*, here in
the stronger form of an asserted impossibility.

*Mitigating.* **No reported number comes from the unconverged region.** I checked
every one: `out_fv` is confined to `Π ≤ 0.95 Π_crit`, `SWEEP_MAXDEV` is computed
over that range only, and every quantity quoted at Itoh's Π (`L_ext_over_tube`,
`y_H_required`, `Pi_itoh/Pi_crit`) comes from the closed-form `G`, not from a
solve. So this is a wrong statement, not a wrong number.

**Fix.** Either (a) restate it accurately — "the flat-start damped Newton stalls
once the extinction point enters the domain, because the sink's Jacobian is
singular at `w = 0` and `NumJac`'s absolute perturbation floor cannot resolve the
sub-1e-6 tail; a marched or continued initial guess converges to residual 1e-13
at every Π tested, including Itoh's" — or (b) better, add the march as the
initial guess (≈ 6 lines) and run the FV sweep over the whole range, which turns
a stated defeat into a demonstration. Either way the sentence about no
non-negative solution must go from all six locations.

---

### F2 — MAJOR, CONFIRMED. The headline α_H check is circular as the page argues it — and it need not be, because Figure 1 independently dimensions the bore and the page discards that witness

**(a) The circularity is real.** `r_i = 8.5 mm / r_o = 8.7 mm` are chosen *by
requiring the α_H identity to close* (H1.4 established the choice this way; the
page re-derives it the same way and says so). The page then reports the residual
of that same identity, −0.032 %, as its headline validation, and check 2 closes
with:

> The two defects it DOES catch decisively are a wrong constant and a
> **wrong pair of radii** — which is what check 1 is actually for.

That sentence is false. Check 1 cannot catch a wrong pair of radii: it selected
them. What it can do is *discriminate between the two candidate readings of one
printed sentence*, which is the derivation of the radii, not a test of them. The
same wording propagates verbatim to `README.md`, `meta.yaml` (`validation:`
bullet 2) and `queue_cases/H1.1.yaml` (`notes:`) as "decisive about the constants
**and the radii**".

**Quantified, as asked.** α_H = 2π l₀ D C₀ / ln(r_o/r_i) constrains only the
*ratio*, i.e. (to first order) the wall-thickness-to-mean-radius ratio. Given the
printed 200 µm wall:

| r_i | α_H | dev vs printed |
|---|---|---|
| 8.30 mm (the OD reading) | 4.3646e-05 | −2.357 % |
| 8.40 mm | 4.4166e-05 | −1.195 % |
| **8.50 mm (used)** | 4.4686e-05 | **−0.032 %** |
| 8.60 mm | 4.5205e-05 | +1.131 % |
| 8.70 mm | 4.5725e-05 | +2.293 % |

Agreement degrades to 1 % at **r_i = 8.417 mm and 8.589 mm** — a window of
−83/+89 µm, i.e. **±1.0 % in r_i**. So the answer to "how far do the radii have
to move" is *barely at all* — but the flip side is that the identity does pin the
mean radius to ±1 % **given** the printed thickness. Drop the thickness and
nothing is pinned: the locus α_H = 4.47e-5 is `r_o/r_i = 1.023522` exactly, which
(r_i, wall) = (4.25 mm, 100 µm), (8.5 mm, 200 µm), (17 mm, 400 µm) all satisfy.
Thickness sensitivity at fixed r_i is much sharper: ±10 µm on the wall moves
α_H by ∓5 %.

**(b) The fix is free, and the page is leaving evidence on the table.**
**Figure 1 prints the dimensions**: `0.2`, `17φ`, `28φ` and `140`, with "[mm]".
I read the arrowheads on the native bitmap of journal page 1577 (page-image rows
quoted; upper Pd membrane line = rows 262–266, lower = rows 343–348; upper shell
band = rows 231–238, lower = rows 369–377):

| dimension | arrow column | upper apex | lower apex | terminates on |
|---|---|---|---|---|
| `0.2` | 705 | row 259 (pointing **down**) | row 267 (pointing **up**) | *straddles* the membrane band from outside — it is the wall thickness |
| `17φ` | 754 | row 268 (pointing **up**) | row 342 (pointing **down**) | the **inner** faces of the two membrane lines |
| `28φ` | 808 | row 239 (pointing **up**) | row 368 (pointing **down**) | the **inner** faces of the shell wall |

The convention is consistent across all three, and it is self-consistent: the
wall is dimensioned *separately* as 0.2, so 17 cannot include it. **Figure 1
dimensions the palladium tube's inner diameter as 17 mm** — i.e. r_i = 8.5 mm,
r_o = 8.7 mm — independently of α_H, and in contradiction to the text's "17.0 mm
OD". (H1.4's verifier reported this as a "soft" second witness; re-read from the
bitmap it is not soft. The drawing is schematic — the 0.2 mm wall is drawn ~5 px
where 17 mm is 75 px — but the arrowhead *placement* is unambiguous at native
resolution, and the draughtsman's choice of which faces to dimension is exactly
the information needed.)

The page instead states flatly, in the parameters table and in `README.md`,
*"They are not printed"*, and rests the whole reconstruction on the α_H match.
The sidecar `data/itoh-1987-permeation.meta.yaml` does mention the Figure 1
observation, but only as something H1.4's verifier noticed — a borrowed claim,
not re-derived, and it never reaches the page.

*Failure scenario.* A reader who notices that the radii were fitted to α_H reads
the −0.032 % as arithmetic self-agreement and discounts the page's one external
comparison entirely; or a later page cites "α_H reproduced to 0.03 %" as
evidence about the tube geometry, which as written it is not.

**Fix (three sentences, no new computation).**
1. Say plainly in check 1 what the residual tests: *the transcription of D, C₀,
   l₀ and the wall thickness, the structure of Eq. 1, and the radial assembly —
   given the geometry*. Delete "and the radii" from check 2's closing bullet,
   `README.md`, `meta.yaml` and the case file, or restate it as "it discriminates
   the two candidate readings of '17.0 mm' by 2.4 % against a 0.03 % residual".
2. Add Figure 1 as the **independent** witness, with the arrowhead reading above.
   Then the two readings become mutually validating in the `J3.4` pattern
   `handoff.md` records ("two independent readings of the same paper can validate
   each other"), and the circularity disappears rather than being confessed.
3. Note that α_H fixes only `r_o/r_i`, so the 200 µm wall is load-bearing: it is
   what converts the ratio into a radius.

---

### F3 — MAJOR, CONFIRMED. Two of the three "Itoh's own tube" headline numbers describe a clean-permeate separator, not Itoh's reactor, and every summary line drops the qualifier

The notebook is careful — the cell prints *"Itoh's membrane, **used as a
separator** for the duty his reactor sets"*. `README.md`, `meta.yaml` (`adds:`),
`models_entry.yaml` (`description:`) and `queue_cases/H1.1.yaml`
(`scope_decision:`) all drop that clause and assert the numbers of Itoh's
reactor:

> Itoh's own tube: `Pi = 51.4`, so it needs **5.5 %** of its length to strip the
> hydrogen his reaction makes and **holds the driving mole fraction at 3.8e-4**.
> That *explains* — rather than restates — `H1.4`'s verified finding …

Both italicised claims are wrong for Itoh's reactor, because his permeate side is
**loaded** (a 11.8e-5 mol/s argon purge that accumulates H₂), while Π_crit,
L_ext and `y = Π^{−2}` are all derived for `y' = 0`.

At H1.4's stated run, the co-current fast-permeation/equilibrium ceiling forces
**equal H₂ mole fraction on both sides**. Using H1.4's own printed inputs
(u_C⁰ = 2.90e-7, y_C⁰ = 0.197, v_A⁰ = 11.8e-5, X = 0.9983):

```
q = 3 X u_C0 / (F0 + v_A0) = 8.6852e-7 / 1.19473e-4 = 7.2696e-3
y = q/(1+q) = 7.217e-3
```

(identical to `ceiling()` in `pages/H1.4-…/build_page.py`, lines 564–574.) So the
reaction-side hydrogen mole fraction in Itoh's reactor is **7.2 × 10⁻³, a factor
19 above the 3.8 × 10⁻⁴ the page quotes** — and the *driving* mole fraction
difference tends to **zero**, not to 3.8e-4, because that is what the ceiling
means. Likewise the hydrogen is never "stripped" in 5.5 % of the length: it is
stripped until the two sides equalise, which the page's own check 5b establishes
is an asymptotic approach to `β y'/(1−y')` with **no finite extinction length at
all**.

The *conclusion* survives — Π = 51 ≫ Π_crit does quantify "permeation is not the
limiting resistance", which is one of the two conditions (with equilibrium
kinetics) that put H1.4 on its algebraic ceiling, and H1.4's verifier measured
that ceiling as exact (`stated_run_fv_vs_ceiling_absdiff = 0.0`). But the
*mechanism sentence* names the wrong mechanism, and the case file leans on it as
"the decisive argument for a separate page".

*Failure scenario.* A reader takes 3.8e-4 as Itoh's hydrogen partial pressure and
uses it to size a purge, or to argue the membrane held the reaction 19× further
from equilibrium than it did. Or a later page cites "the tube strips its hydrogen
in 5.5 % of its length" as a fact about the published reactor.

**Fix.** Carry the notebook's own qualifier into all four sidecars: *"the same
membrane, run as a separator on a clean permeate against the duty his reaction
sets, would strip it in 5.5 % of the tube at a driving mole fraction of 3.8e-4"*.
Then state the H1.4 link at the level it actually holds: *"Π = 51 against
Π_crit = 2.8 says permeation is nowhere near limiting, which is the condition
that puts H1.4's model on its co-current fast-permeation ceiling — and once there
the conversion is a function of K_p and the purge split alone, so α_H has no
observable left to move."* Optionally add the loaded-permeate number (7.2e-3) —
it costs one line and it makes the connection quantitative *and* correct.

---

### F4 — MODERATE, CONFIRMED. Π_crit is conditional on a clean permeate, and §3/§4 present it unconditionally

`Π_crit = G(1,β)` exists only because `∫₀ dw/√w` converges, which requires
`y' = 0`. With any permeate loading the finite extinction length disappears
entirely. The page *knows* this — check 5b says so explicitly and measures it —
but the statement lives 3 000 words downstream of §3's

> A ½-law strips its feed completely in a finite length … `Π_crit = G(1,β)` — a
> number, not a fit

and of §4's map, which is drawn on the `y' = 0` axis with no caveat, and of the
"What pymrm adds" summary. Given F3, this is the caveat that matters most.

**Fix.** One clause in §3(a) and one in §4: "for a clean permeate side; with
`y' > 0` the tube strips only to `β y'/(1−y')` and the extinction length becomes
infinite — see check 5b." Π_crit is *not* an artefact of the closure otherwise:
I reproduced `G(1,β) = 2.8388708811` by direct quadrature of
`∫₀¹ √((w+β)/w) dw = 2.8388708810` (1.6e-11 relative), and the inverted-`G`
solution against an independent RK45 integration to 5.2e-12.

---

### F5 — MINOR/JUDGEMENT. The reference block should follow `J3.1`, not `F1.3`

The builder set `reference` = Bohmholdt & Wicke (1967), `reference_read_from` =
Itoh. That is a literal reading of the `AGENTS.md` reprint convention, and the
attribution behind it is **verified**: see "What survived" below.

I recommend changing it anyway, to the `J3.1` form:

```yaml
reference:                       # the paper actually read, and the source of every equation
  authors: ["Itoh, N."]
  year: 1987
  container: "AIChE Journal 33(9) 1576-1578"
  doi: 10.1002/aic.690330921
origin_not_consulted: >
  Bohmholdt, G. and Wicke, E. (1967), Z. physik. Chem. Neue Folge 56, 133 — the
  half-power pressure law, as Itoh attributes it; Sieverts, A. and Danz, W.
  (1936), Z. physik. Chem. (B) 34, 158 — the solubility constant C_0 = 1280
  mol/m3; Nagamoto, H. and Inoue, H. (1985), Chem. Eng. Commun. 34, 315 — the
  diffusivity D. All three are recorded without DOIs because no copy was
  available and no citation could be verified against an original. Not
  consulted; nothing on the page derives from them beyond the two constants
  Itoh quotes.
```

Reasons:

1. **The reprint test is not met for B&W.** `AGENTS.md` requires that the paper
   on disk *print the result in full with attribution* — Krishna & Ellenberger
   reprint Wilkinson's Eqs. 1–4 verbatim, which is why `F1.3` may put Wilkinson
   in `reference`. Itoh does not reprint anything of Bohmholdt & Wicke's: he
   cites them for the *exponent*, then writes his own Eq. 1 with his own
   cylindrical permeance and his own constants. Everything this page implements
   and checks is Itoh's.
2. **The citation is unverifiable.** Container and year come from Itoh's
   reference list only. Putting it in `reference` makes the canonical index's
   primary source for H1.1 a paper that cannot be checked and has no DOI —
   exactly the situation `J3.1`'s `origin_not_consulted` was written for, and
   exactly the "a DOI resolved from a terse citation is usually wrong" hazard.
3. **It is not the origin of the thing the case is named for.** The case is
   "Sieverts law H₂ permeation"; Sieverts' law is Sieverts', and B&W is Itoh's
   attribution for the permeation form. `origin_not_consulted` can record both
   without asserting either as *the* reference.
4. **Consistency with H1.4**, which carries `reference: Itoh` for the same paper.
   `B1.6` — also a reprint-route page — likewise keeps the read paper in
   `reference` and puts Prater in `also`.

Keep the excellent `note:` text; move it under `reference`. If the integrator
prefers to keep `reference_read_from`, at minimum add `doi: null` visibility and
say in the note that the B&W citation was transcribed from Itoh's reference list
and not verified.

---

### F6 — MINOR, CONFIRMED. The "217×" is partly computed from non-solutions

`channel_break_min_ratio = 217.29` is `2.113e-01 / 9.722e-04`, and the numerator
comes from the "mole fraction replaced by the flow w itself" row, which the table
itself annotates *"(Newton did not converge)"*. The `nu = 1` row is also
unconverged. The claim "every injected defect moves the error by at least 217×"
is therefore in part a statement about a failed solve rather than about a wrong
solution. Harmless (the defects that *did* converge move it 340× and 2900×), but
say so, or quote the minimum over the converged breaks.

### F7 — MINOR. "1.16 %" and "−1.19 %" are the same defect against different references

The break table prints the `nu = 0` row as −1.19 % (vs the *printed* α_H) while
`README.md`, `meta.yaml`, the case file and the "Reuse" trap all quote 1.16 %
(vs the *baseline solve*). Both are correct and the page's stated convention
covers it, but a reader comparing the sidecar with the table will see two
numbers. Pick one and name its reference.

### F8 — MINOR. `meta.yaml` claims `runtime_seconds: 7`; measured 6.7 s. Fine. `pymrm_min_version: "2.3"` is satisfied (installed 2.3.1.dev3).

---

## What survived the attack (no finding)

* **The attribution is exactly as the builder reports it, verified on the native
  bitmap.** Page 1576, right column: *"assumed to obey the half-power pressure
  law (Bohmholdt and Wicke, 1967), i.e., it is proportional to the difference
  between the roots of the hydrogen partial pressures in the reaction and
  separation sides, p_H and p′_H:"*, followed by Eq. 1 printed exactly as the
  page quotes it, `Q_H = α_H(√(p_H/P_o) − √(p′_H/P_o))`, `α_H = 2πl_o/ln(r_o/r_i)
  · D C_o`. Page 1577, left column: *"9.23 × 10⁻¹⁰ m²·s⁻¹ for D (Nagamoto and
  Inoue, 1985) and 1,280 mol·m⁻³ for C_o (**Sieverts and Danz, 1936**) at 473 K,
  were used, α_H calculated was 4.47 × 10⁻⁵ mol·s⁻¹."* All three superscripts
  (−10, −3, −5) read cleanly at native resolution, as the page claims. **Sieverts
  & Danz is cited for C₀ and for nothing else**, and the page asserts nothing
  about Sieverts' own paper anywhere — the Background derives `c = C₀√(p/P₀)`
  from dissociative equilibrium rather than citing it, and the Sources block
  lists all three cited classics under "not consulted here". No misprint anywhere
  in the permeation law, confirmed independently.
* **Determinism.** Re-executing `index.ipynb` reproduces every stream output
  byte-identically apart from an `ipykernel` PID inside a `MatrixRankWarning`
  path; `agreement.json` is unchanged; runtime 6.7 s. No warm-start continuation
  anywhere — every `Permeator` starts from a flat guess and the wall solve is
  linear. Complies with the `B1.1` determinism lesson.
* **The wall assembly is right, verified against three independent analytic
  geometries I derived myself** (not the page's code):
  cylindrical `2πl₀DC₀/ln(r_o/r_i)` = 4.4686e-05; the `nu = 0` slab
  `2πr_i l₀DC₀/(r_o−r_i)` = 4.4167e-05; the `nu = 2` shell
  `2πl₀DC₀·r_o/(r_o−r_i)` = 4.5206e-05. Every row of the break table matches to
  4–5 digits, including −89.00 % for the 10× wall, +899.68 % for the slipped `D`
  digit, and 45.1 % for the `nu` defect at `r_o/r_i = 3`. **So the two admitted
  blind spots are real and correctly attributed**: the 1.16 % `nu` insensitivity
  is a property of a wall 2.35 % of the tube radius, not of the code, and the
  2.2e-9 at `n_r = 3` is a property of a log profile whose flux a 3-cell FV
  scheme integrates almost exactly. Publishing both is strictly better than
  publishing the residual alone, exactly as `B1.6` now does.
* **The channel assembly is right.** I wrote a first-order upwind FV residual
  and a bidiagonal Newton from scratch, sharing no code with the page, and got
  `max|w − w_exact|` = 7.719e-03 / 9.722e-04 / 2.432e-04 at n_z = 100 / 800 /
  3200 — the page's numbers to four significant figures.
* **The new results reproduce independently.** β = 1.692047, Π = 51.3793,
  Π_crit = 2.838871 (quadrature 2.8388708810 vs the closed form 2.8388708811),
  L_ext/L = 0.055253, y = Π⁻² = 3.7881e-04, Π/Π_crit = 18.1. The inverted-`G`
  reference agrees with an independent RK45 integration to 5.2e-12. The `y ~
  Π^{−1/n}` scaling is right by construction (`q = α_n y^n`), and the finite-vs-
  infinite extinction argument is right (`∫₀ dw/√w` converges, `∫₀ dw/w` does
  not).
* **The 18× oversizing survives a redefinition of the duty**, which the builder's
  `NOTES.md` asked the verifier to test. β is correctly assembled from printed
  numbers (Ar + benzene over 3 u_C⁰). Alternatives: counting argon only
  (β = 1.359) gives Π/Π_crit = 19.8; using the measured X = 0.997 instead of
  complete conversion gives 18.1; taking the duty at the *membrane-free
  equilibrium* X = 0.187 gives 48.8. Every defensible variant lands in
  18–50×, so the "far above Π_crit" conclusion is robust even though the
  headline 18.1 is not a sharp number.
* **No quoted number comes from the unconverged region** — checked line by line
  (see F1).
* **Deviation convention.** `(model − reference)/reference` is stated once and
  used everywhere I checked, including the break tables (where the reference is
  named as the baseline solve rather than the printed value — see F7).
* **pymrm conventions.** `NumJac(self.shape)` with `shape = (n_z, 1)` — the
  corrected `AGENTS.md` guidance, correct; no bare 1-D shape and no
  `axes_diagonals` on `ndims = 1` anywhere. `construct_div(..., nu=1)  # nu=1:
  cylindrical` and `nu=0: axial` are both commented. Boundary conditions are
  written on the outward normal with the physical equation in a comment at every
  `bc`, and the outlet-`None` singularity trap is named. Operators are assembled
  once in `__init__`.
* **Tier 6 is stated plainly and nothing is called experimental.** "Nothing on
  this page is compared with a permeation measurement, because the source does
  not print one" appears in the title block, the Data section, Validation's
  preamble, `meta.yaml` `caveats:`, `README.md` and the CSV sidecar. Itoh indeed
  prints no flux, no permeance, no pressure sweep and no temperature dependence —
  confirmed by reading the whole paper. The Caravella non-use is disclosed in
  four places and the follow-up is recorded.
* **The series-resistance limits are labelled as the algebraic identity they
  are** ("Both limits come out of the same solver … it cannot fail on physics").
  Correct and correctly labelled.
* **Check 4 (graded `D` vs a quadrature of the resistance integral) is the one
  genuinely independent numerical check on the page**, it is break-tested
  (`nu = 0` moves it to −1.03 %), and the "what pymrm adds" claim built on it is
  honest — including the opening admission that for Itoh's actual uniform wall
  pymrm adds nothing.
* **Scope vs H1.4: genuinely distinct.** H1.4 quotes the printed α_H and consumes
  `r_i` only through `V_r`; H1.1 derives α_H from a radial BVP, adds the
  single-domain permeator with a prescribed permeate, and adds Π/Π_crit, the
  extinction-length structure and the surface-resistance apparent exponent. None
  of that is on H1.4. The `structures` change `[S7] → [S3, S5]` is correct — the
  permeate side is a prescribed boundary value, so there is no second solved
  domain. No prose, number or dataset is reused from H1.4 except the radii
  reconstruction, which is re-derived here (and see F2 for what re-derivation
  should have included).

---

## Explicit verdicts

**(1) Is the α_H check circular, and how should the page state it?**
**Partly circular as argued — CONFIRMED — and the page is having it both ways.**
It concedes blindness to geometry and grid, then claims the identity is
"decisive about the constants **and the radii**", which it cannot be, because the
radii are its own output. The residual tests the transcription of D, C₀, l₀ and
the 200 µm wall, the structure of Eq. 1, and the radial assembly — *given* a
geometry, and given that the true geometry is one of the two candidate readings.
It should say exactly that. **But the page can do much better than confess:
Figure 1 dimensions the bore as 17 mm on the inner faces (F2b), independently of
α_H. With that added, the −0.032 % becomes two independent readings agreeing,
and the circularity is gone rather than admitted.**

**(2) The attribution and the reference block.** The attribution is **verified on
the native bitmap and is handled honestly on the page** — Bohmholdt & Wicke for
the exponent, Sieverts & Danz for C₀ alone, no claim about Sieverts' own paper,
Caravella unused and said so. The **reference block should change**: put Itoh in
`reference` (with the DOI) and record Bohmholdt & Wicke, Sieverts & Danz and
Nagamoto & Inoue in `origin_not_consulted`, the `J3.1` form. Itoh does not
reprint B&W's result, the B&W citation is unverifiable, and B&W is not the origin
of the law the case is named for. See F5.

**(3) Which checks survived my own break tests.** All of them, numerically. I
reproduced every wall break row from independent analytics, reproduced the
channel FV errors with a solver written from scratch, and reproduced Π, Π_crit,
L_ext, y and G by quadrature and by RK45. The two admitted blind spots (1.16 %
for a wrong `nu`, 2.2e-9 at `n_r = 3`) are real, correctly diagnosed as
consequences of a thin wall and a log profile rather than of the code, and
honestly framed *except* for the "and the radii" clause in F2. The 217× channel
figure is real but partly drawn from unconverged runs (F6). Nothing on this page
is a check that cannot fail *and is presented as though it could*.

**(4) Scope vs H1.4.** Distinct, and the page earns its own entry — on the
derivation of α_H, the Π/Π_crit structure, the finite-extinction result and §5's
apparent exponent. It does **not** earn it on the sentence the case file calls
"the decisive argument": the claim that Π = 51.4 explains H1.4's insensitivity
*via* a 5.5 % stripping length and a 3.8e-4 driving mole fraction is the wrong
mechanism for a loaded permeate, and the mole fraction is out by 19× (F3). Fix
the framing and the scope argument still stands on the rest.

---

## Recommendation

**Safe to publish after these fixes:**

1. **F1** — remove the "no non-negative solution" claim from all six locations.
   Either restate it as the Newton globalisation failure it is, or (better) seed
   the solve with a march and extend the FV sweep over the whole Π range.
2. **F2** — delete "and the radii" from check 2, `README.md`, `meta.yaml` and the
   case file; add Figure 1's `0.2 / 17φ / 28φ / 140` dimensions and the
   inner-face arrowhead reading as the independent witness; note that α_H fixes
   only `r_o/r_i`, so the printed 200 µm wall is load-bearing.
3. **F3** — restore the notebook's "used as a separator" qualifier to
   `README.md`, `meta.yaml`, `models_entry.yaml` and `queue_cases/H1.1.yaml`, and
   restate the H1.4 link as "Π ≫ Π_crit ⇒ the fast-permeation ceiling", not as
   stripping in 5.5 % of the length. Quote 7.2e-3 as the reactor's actual
   both-sides H₂ mole fraction if the link is to stay quantitative.
4. **F4** — add the `y' = 0` condition to §3(a) and §4 where Π_crit is
   introduced.
5. **F5** — reference block to the `J3.1` form.
6. **F6, F7** — one clause each.

Optional and cheap: the marched initial guess of F1 would let the Π sweep show
finite-volume points at Itoh's own operating point, which is where the page's
most interesting claim lives; and adding the Figure 1 reading of F2 converts the
page's weakest-looking result into its `J3.4`-pattern strongest.

*Verifier: Claude Opus 5, 2026-08-01. Nothing outside this file was modified; no
git command was run.*
