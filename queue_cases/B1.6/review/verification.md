# B1.6 — verification report

Verifier pass on `queue_cases/B1.6/page/`, 2026-07-31. Source consulted directly:
`~/papers/pymrm-gallery/1-s2.0-0009250962850052-main.pdf` (Weisz & Hicks 1962),
journal page 266 rendered at 600 dpi. Attack scripts in the session scratchpad
(`VB16/attack_a.py`, `attack_b.py`, `attack_f.py`, `attack_t.py`).

**Verdict: safe after fixes.** The physics, the derivations and the transcription
are all correct — I could not break any of them. What is wrong is the page's
account of *what its numbers prove*. Three specific claims about the power of the
324-solve residual are demonstrably false, and the page's own honest caveat is
weaker than the truth it is trying to state.

---

## Verdict on the two questions asked

**(1) Does the 324-solve residual have any power?** Partially — much less than the
page says, and less than the page's own caveat concedes.

`ε` is **blind** to: the discretisation, the geometry index, and whether the
nonlinear solve converged at all. It **does** move for: a sign or scale error in
either source term, a boundary-condition *type* mismatch between the two fields,
and a rate evaluated at inconsistent states. Demonstrations below.

**(2) Is the claimed sign contradiction in Weisz & Hicks real?** **Yes —
CONFIRMED, and it is not a convention error by the builder.** Read off the 600 dpi
render of page 266, verbatim:

- eq. (4) `D∇²c − dn/dt = 0`
- eq. (5) `K∇²T − H dn/dt = 0`
- eq. (7) `ΔT = T − T₀ = − (HD/K)(c₀ − c)`  ← the leading minus is unambiguous
- eq. (8) `γ = Q/RT₀,  β = c₀HD/(KT₀) = (ΔT/T₀)_max,  y = c/c₀`
- eq. (8) rate `dn/dt = k₀c₀ y exp{γβ(1−y)/[1 + β(1−y)]}`
- eq. (6) `dn/dt = k₀c exp[(Q/RT₀)(ΔT/T₀)·1/(1 + ΔT/T₀)]`

The paper's own prose (right column, page 266) says β "is the maximum temperature
variation (ΔT)max ... relative to the boundary temperature; **this follows from
equation (7), by c → 0**." Substituting c = 0 into eq. (7) as printed gives
(ΔT/T₀)_max = **−**c₀HD/(KT₀), the negative of the middle member of eq. (8). The
two printed members of eq. (8) therefore cannot both be right, and the paper
itself names eq. (7) as the route between them. This is a statement about printed
symbols and is independent of what sign convention `H` carries — it holds whether
`H` is the enthalpy (negative exothermic) or the heat released (positive
exothermic). **Not an artefact of assuming exothermic ⇒ H > 0.**

All three of the builder's selecting arguments check out independently:

1. **Check C** — verified against the images. eq. (6)'s exponent is `γx/(1+x)`
   with `x = ΔT/T₀`; eq. (8)'s is `γβ(1−y)/[1+β(1−y)]`. `x = +β(1−y)` closes,
   `x = −β(1−y)` does not. **Not circular**: eq. (6) and eq. (8) are two
   separately printed equations, and the contradiction in finding (2) is
   established from (7) vs (8), not from (C).
2. **The sweep range** — page 268 as printed: "β in the range from 0 to +0·8
   (exothermic reaction) and 0 to −0.8 (endothermic reaction)". Verified.
3. **η > 1** — "For the case of exothermic reactions it is seen that η can become
   larger than unity for large enough values for β". Verified.

Reference [11] on the reference page reads `PRATER C. D., Chem. Engng. Sci. 1958
8 284`, and [10] `DAMKÖHLER G., Z. Phys. Chem. 1943 A193 16` — both cited
correctly on the page. The Prater-not-consulted framing is honest and stated in
three places.

---

## Findings, ranked

### 1. CONFIRMED — the page claims ε "measures the nonlinear solve". It does not.

`ε` is at roundoff even when Newton has not converged and the solution is
physically impossible.

| `maxfev` | Newton residual | y at centre | ε |
|---|---|---|---|
| 100 | 3.3e-11 | 1.0e-06 | 2.3e-13 |
| 3 | 7.5e+00 | 0.730 | 1.1e-11 |
| 2 | 2.2e+01 | 0.573 | 6.7e-13 |
| 1 | 1.8e+01 | **2.12** (y > 1, impossible) | **8.5e-12** |

The `maxfev=1` iterate is garbage — the concentration exceeds the surface value —
and its ε is *inside* the page's headline worst of 1.1e-11.

**Failure scenario.** A reader — or a future builder copying this page as the
recommended "free unit test" — concludes that a small ε means the pellet solve
converged. It does not. A silently non-converged Newton passes the test.

**Fix.** ε is a linear-algebra roundoff measurement. Say that: the discrete
identity holds at any *root* of the discrete system and is inherited from the
shared operator; ε measures conditioning and floating-point error, not
convergence. Report the Newton residual alongside it (the code already computes
`rn` and asserts `rn < 1e-8` — surface that number on the page).

### 2. CONFIRMED — "a mismatched geometry index" is not something ε can catch.

Two separate problems with this claim:

- A *per-field* mismatch is impossible by construction: `construct_div` is called
  once for shape `(n_u, 2)`, so one `nu` serves both columns. The page cannot
  produce the defect it says the test would catch.
- A *wrong but shared* `nu` is exactly what the identity is blind to — Prater
  holds for any geometry, which is the page's own headline. Solving a sphere with
  `nu = 0`:

| run | η | ε |
|---|---|---|
| `nu = 2` (correct) | 3.856 | 2.1e-13 |
| `nu = 1` (wrong) | 2.951 | 6.7e-13 |
| `nu = 0` (wrong) | **1.670** (57 % error) | **4.7e-12** |

ε stays inside the sweep's reported band while η is wrong by more than half.

The same claim is repeated in **"What pymrm adds"** ("...means a sign, a boundary
condition or a **geometry index** is wrong") and in `meta.yaml`'s `adds:` field.
Three places, all wrong.

**Fix.** Delete "geometry index" from all three. The correct list is: a sign or
scale error in either source term, a boundary-condition *type* mismatch between
the fields, and a rate evaluated at inconsistent states — the three I confirmed
do move it (below).

### 3. CONFIRMED — ε is blind to discretisation error, so the sweep's severity is decoration.

The page argues the sweep is "not a soft test of the physics: it reaches η = 162
and depletion to y = 1.6e-195". Severity has no bearing on ε:

| `n_u` | η | ε |
|---|---|---|
| 200 | 3.856 | 2.1e-13 |
| 40 | 3.859 | 1.6e-14 |
| 12 | 3.892 | 6.8e-15 |
| 6 | 3.734 | 1.5e-15 |
| **3** | **5.276** (37 % error) | **2.4e-15** |

ε gets *better* as the grid gets worse. The 324 solves and their extreme states
add no evidential weight to the identity beyond what one solve gives.

### 4. CONFIRMED — Section 5 does **not** remove the shared-operator objection.

`bvp_pellet` genuinely shares no code with the pymrm path (self-contained `rhs`,
`bcs`, scipy only) — that part of the claim is true. But the collocation inherits
the invariant for exactly the same structural reason: `w = θ + βy` and
`w' = θ' + βy'` form a *closed linear subsystem* of the first-order form
(`dw' = −ν dw/x`, the rate cancels), and scipy applies the same mesh and the same
collocation nodes to every component. So the discrete `w` is forced to `1+β`
regardless of mesh quality:

| `tol` | mesh nodes used | ε | y_min |
|---|---|---|---|
| 1e-8 | 19775 | 1.3e-15 | 1.0e-06 |
| 1e-3 | 508 | 4.8e-15 | 1.0e-06 |
| 1e-1 | 29 | 4.2e-13 | 1.0e-06 |
| 1.0 | 7 | 1.7e-11 | 6.0e-06 |
| **1.0** | **5** | **2.0e-11** | 2.6e-05 |

A 5-node mesh — a useless solution — gives ε = 2.0e-11, worse than nothing but
still in the same band as the headline. The sentence "Section 5 removes even the
shared-operator objection" is false and should go.

**What would have removed it**, and costs almost nothing: the page never compares
the pymrm profile with the `solve_bvp` profile. It computes ε separately inside
each. Comparing `y(u)` between the two solvers *would* be a genuine independent
check of the pymrm discretisation. As written, Section 5 adds nothing that
Section 4 does not already have.

### 5. CONFIRMED — the film "closed form matched to 1.9e-11" is tautological too.

The **derivation is correct** — I re-derived it independently and it is right:
`w = θ + βy` is harmonic ⇒ constant ⇒ `dw/dn|_s = 0` ⇒
`Bi_h(1−θ_s) + βBi_m(1−y_s) = 0` ⇒ `θ_s − 1 = β(Bi_m/Bi_h)(1−y_s)`, and
subtracting gives `β(1−y_s)(Bi_m/Bi_h − 1)` uniformly in position. Both steps are
forced by the *discrete* harmonicity, so the numerical match cannot fail:

| `n_u` | Bi_m/Bi_h | ε_bulk | closed form | relative difference |
|---|---|---|---|---|
| 400 | 100/10 | 2.5795e-01 | 2.5795e-01 | 4.3e-12 |
| 400 | 10/100 | 1.7269e-01 | 1.7269e-01 | 1.9e-11 |
| 40 | 100/10 | 2.5795e-01 | 2.5795e-01 | 9.3e-13 |
| **8** | 100/10 | 2.5795e-01 | 2.5795e-01 | **5.6e-14** |
| **8** | 10/100 | 1.7251e-01 | 1.7251e-01 | **1.4e-14** |

At `n_u = 8` — eight cells for φ = 3 — the match is *better* than at 400. The
surface balance itself is exact to 3e-14…5e-12 at every Biot pair I tried.

The closed form is also implicit, not predictive: it contains `y_s`, which is a
solve output. This does not make it wrong, but "the failure is not merely observed
— it is predicted exactly" (page) overstates it. What the section genuinely
contributes is the **derivation** and the **magnitudes** (5.9 %, 26 %, 2881 % of
β; 8.8 %, 24 %, 99.6 % of the rise sitting in the film) — those come from the
solved `y_s` and are real. Reword the "reproduced to 1e-11" line accordingly.

### 6. CONFIRMED — the cross-page β recovery is the identity restated, not a cross-check.

`BETA_RECOVERY = max |(θ−1)/(1−y) − 0.6|` is algebraically the same quantity as
the sweep's ε, divided by `(1−y)`. It is not an independent measurement of
anything. Its only genuine content is that B1.1/B1.5 use the same *numerical
value* 0.6 for the same symbol — a convention check, worth having, but not the
validation item it is listed as.

Likewise the two-field ↔ reduced-equation match (6.1e-12): I checked B1.1's
notebook, and its formulation **does** presuppose Prater — its assumption list
says so explicitly ("Prater relation exact"), and `wh_rate` folds the relation
into the rate. Given the exact discrete identity, `θ = 1 + β(1−y)` exactly, and
`γ(1 − 1/θ) ≡ γβ(1−y)/(1+β(1−y))` identically, so the two-field `y` is an exact
root of the reduced discrete system on the same grid. The agreement is guaranteed.

It does still check one real thing — that the two pages' `γ` and `β` conventions
coincide — and the page should say that is what it checks.

**The one genuinely independent numerical check on the page is the shooting
comparison** (3.9e-05 at `n_u = 800`, observed order 2.11). That tests pymrm's
discretisation against a method that discretises nothing. It should be given more
weight than the 1e-11s, which currently dominate the headline.

### 7. NOT A FINDING — the transient breakdown is real. Independently confirmed.

This is the page's strongest genuine contribution and it survives attack. I
re-ran the transient on a **completely independent node-based finite-difference
discretisation** (nodes at `u_i = i·h`, my own sphere stencil with the `6(f₁−f₀)/h²`
centre form, no pymrm anywhere):

| Le | page (N=100) | independent FD, N=100 | independent FD, N=200 |
|---|---|---|---|
| 10 | 0.251926 | 0.251938 | 0.251918 |
| 1 | 3.5e-13 | 1.1e-15 | 1.0e-15 |
| 0.1 | 0.123960 | 0.123966 | 0.123968 |
| 0.01 | 0.145900 | 0.145902 | 0.145903 |

**The Le = 10 number the brief asked about specifically: 0.2519 — confirmed to
four digits on two grids and two discretisations, peak at τ ≈ 0.05 in both.**
The pymrm run at N = 200 gives 0.2519078 against 0.2519255 at N = 100 — converged
to four digits. The `t_eval` geomspace samples at 4.5 % spacing, far finer than
the peak's width.

**The violation genuinely vanishes at steady state and is not masked by the
convergence criterion.** In the independent run the final value is 1.1e-15
(Le=10) and 4.6e-15 (Le=0.1) — floating-point zero, on a scheme with no shared
solver machinery. Le = 0.01 ends at 1.03e-05 in *both* implementations, and the
page correctly states that this case "is still relaxing at the end of the window"
rather than claiming zero. Honest.

The `Le` scaling in the markdown is dimensionally correct: I checked
`Le·β·φ²·R` against `εL²(−ΔH)r/(ρ_p c_p T_s D_e)` and they are the same group.

### 8. CONFIRMED (minor) — assumption-table row 3 is wrong as written.

> | 3 | constant D_e and λ | fails unless each coefficient depends only on its own variable |

If `D = D(c)` and `λ = λ(T)`, define Kirchhoff potentials `Φ_c = ∫D dc`,
`Φ_T = ∫λ dT`. Then `∇²Φ_c = r` and `∇²Φ_T = −(−ΔH)r`, so what is harmonic is
`Φ_T + (−ΔH)Φ_c` — a relation between the *potentials*. The linear relation
between `T` and `c` that the whole page is about does **not** survive; it becomes
nonlinear. The row currently reads as "the relation holds if each coefficient
depends only on its own variable", which is false.

**Fix.** "holds only after a Kirchhoff transformation — the relation between T and
c itself becomes nonlinear (not tested here)". `meta.yaml` already says this is
untested; the claim just needs to be correct.

### 9. CONFIRMED (minor) — "any particle geometry ... tested rather than assumed" was not tested.

"What pymrm adds" says `construct_div` "also takes a callable `nu` for an
arbitrary area profile, **which is how Prater's 'any particle geometry' is tested
rather than assumed**." No callable `nu` appears anywhere in the notebook; only
`nu = 0, 1, 2`. And per finding 2, ε could not detect a geometry error even if it
were varied. Reword to "could be tested with", or drop.

### 10. CONFIRMED (minor) — a printed symbolic result is wrong for the cylinder.

`sp.dsolve` prints `general solution of laplacian(w) = 0 : C1 + C2*u**(1 - nu)`.
At `nu = 1` this degenerates (the second solution is `ln u`, not `u⁰`), so the
displayed general solution is invalid for exactly one of the three geometries the
page sweeps. The *conclusion* is unaffected — the regular symmetric solution is
the constant in every case, and the markdown's maximum-principle argument is
correct and geometry-free. But the printed line is a false statement on the page.
Either note the `ν = 1` exception or drop the `dsolve` line and keep the
maximum-principle argument.

---

## What I could not break

- **The transcription.** Eqs. (4), (5), (6), (7), (8) are exactly as the page
  states them. I read all five off the 600 dpi render myself.
- **The sign analysis.** Real contradiction, correctly diagnosed, correctly
  resolved, alternatives printed, nothing repaired by inference.
- **The symbolic proof.** `∇²(θ + βy) = 0` for arbitrary `R` and arbitrary `ν`;
  correct.
- **The film derivation.** Re-derived independently; identical.
- **The transient theory and numbers.** Confirmed against an independent
  discretisation (finding 7).
- **What ε *does* catch.** All three of the mechanisms the corrected claim should
  list do move it, sharply:

  | injected defect | ε |
  |---|---|
  | baseline | 2.1e-13 |
  | heat source uses `−β` | 4.2e-01 |
  | heat source uses `1.01β` (1 % scale error) | 1.0e-02 |
  | y Dirichlet, θ Robin (BC type mismatch) | 1.6e-01 |
  | rate at `θ = 1` in the energy equation only | 5.4e-01 |

  The 1 % scale error landing at exactly ε = 1.0e-02 is a nice property — the
  test is linear in the defect. Worth saying on the page; it is a stronger
  statement than the ones currently made.

- **Determinism.** Re-executed the staged notebook end to end. Every printed
  number is byte-identical to the staged outputs. `solve_ramp`'s fixed
  `(8, 24, 72)` ladder is genuinely deterministic — no warm-start chain along a
  swept parameter (the B1.1 lesson is respected).
- **"All 324 converged" is asserted, not assumed.** `assert not failed` with
  `ok = rn < 1e-8 and isfinite and y.min() >= 0`. Real.
- **Prose numbers match the code.** Every hardcoded number in `README.md`,
  `meta.yaml` and `agreement.json` matches the executed outputs — including the
  whole five-row film table (2.1e-12 / 5.854e-02 / 2.580e-01 / 2.881e+01 /
  1.727e-01 and the 2.8 / 8.8 / 24.2 / 99.6 / 2.7 % column) and the transient
  0.146 / 0.124 / 0.252. No stale markdown numbers.
- **Boundary conditions.** Outward normal, physical equation in a comment beside
  every `bc`. I checked the film pair by hand: `D_e dc/dr|_R = k_g(c_b − c_s)` →
  `dy/dn + Bi_m y = Bi_m`, and `λ dT/dr|_R = h(T_b − T_s)` →
  `dθ/dn + Bi_h θ = Bi_h`. Same sign pattern, correct, and the comment about the
  inward-normal trap is right.
- **`clip_approach` with a per-field lower-bound array.** Real and necessary. The
  bounds differ per field (`0.0` on y, `0.2` on θ to keep `1/θ` finite); a scalar
  cannot express that. `clip_approach` broadcasts, and the `(2·n_u, 1)` tile
  matches the flattened unknown vector. Correct.
- **`y_min > 0` as a dead-core test.** The reasoning in the docstring is right:
  zero order has no `y` in the rate, so a dead core means the idealised model has
  no non-negative steady state, whereas first-order/LH may legitimately underflow
  (the sweep reaches `y = 1.6e-195`). *Small inconsistency:* the sweep passes the
  default `no_dead_core=True` for **all three** rate laws, not just zero order, so
  the code is stricter than the docstring's rationale. Harmless here — nothing hit
  it — but a first-order case that underflowed to exactly 0 would be recorded as
  a convergence failure and trip the assert. Consider passing
  `no_dead_core=(kin == "zero order")`.
- **Tier 6 / not-experimental.** Stated in the categories, "The data", the
  Validation preamble, `meta.yaml` caveats and `README.md`. Nothing is described
  as experimental validation. `models_entry.yaml` carries all the required fields
  and its `slug`/`title` match `meta.yaml` exactly.
- **No Quarto-specific markdown**, section order correct, no data directory
  needed, no figure digitised, no maintainer review gate.

---

## Required fixes before publishing

1. **Section 4, "Two honest remarks".** Replace with an accurate account:
   ε is inherited exactly by any discretisation applying the same operator to
   both fields, so it measures floating-point error — **not** the discretisation,
   **not** the geometry index, and **not** whether Newton converged (demonstrate
   or simply state). What it does catch, sharply and linearly in the defect: a
   sign or scale error in either source term, a BC *type* mismatch between the
   fields, and a rate evaluated at inconsistent states.
2. **Delete the sentence "Section 5 removes even the shared-operator objection."**
   The collocation inherits the invariant for the same reason (finding 4). Either
   drop the claim, or add the check that would earn it: compare the pymrm and
   `solve_bvp` profiles to each other.
3. **Remove "geometry index"** from the "free unit test" list in *What pymrm
   adds*, from `meta.yaml`'s `adds:`, and from Section 4's list. Three places.
4. **Reword the film claim** from "predicted exactly / reproduced to 1e-11" to
   note that the match is forced by the discrete harmonicity (it holds at 8 cells
   too), and that the contribution is the derivation plus the magnitudes.
5. **Label the β-recovery and reduced-equation cross-checks** as consistency
   checks on conventions and implementation, not corroboration of the identity —
   B1.1's formulation presupposes Prater. Promote the shooting comparison
   (3.9e-05, order 2.11) as the page's real discretisation test.
6. **Fix assumption-table row 3** (Kirchhoff, finding 8).
7. **Fix or drop the "callable `nu` ... tested rather than assumed"** sentence
   (finding 9) and the `dsolve` line's `ν = 1` degeneracy (finding 10).
8. Optional: soften the title-cell description and README's "the sweep is not a
   soft test" paragraph, which currently invite the reader to weight the 1e-11s
   as evidence of correctness.

None of these touch a number, a derivation, a transcription or a figure. They are
all statements about what the numbers mean. The underlying work is sound.
