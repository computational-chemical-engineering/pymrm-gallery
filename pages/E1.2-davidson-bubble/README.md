# E1.2: The Davidson bubble

A single bubble in an incipiently fluidised bed is a Laplace problem with a
uniform-pressure hole in it. Solve it once and you get the cloud, the
through-flow, and, because both come from the same field, the reason two
constants that look unrelated in the reactor literature are not.

- **Structures:** `S3` (1D steady BVP)
- **Runtime:** about 10 s
- **Tier:** 6. Nothing here is compared with a measurement.

## Sources, and the distinction that matters

**Kunii, D. & Levenspiel, O.**, *Bubbling bed model for kinetic processes in
fluidized beds*, Ind. Eng. Chem. Process Des. Dev. **7**(4) 481–492 (1968),
[doi:10.1021/i260028a001](https://doi.org/10.1021/i260028a001), **the only text
read.** It prints the results as its equations 3, 9, 10, 11 and 17 and, what makes
the reprint carry the case, *evaluates* each of them with the arithmetic shown and
the answer printed in three worked appendices.

**The attribution, at its actual strength.** The name Davidson appears in the body
of that paper exactly **once**: on journal page 490, inside a parenthesis listing
the model's simplifying assumptions: *"many simplifying assumptions were made in
developing this model (the Davidson bubble, bed with single size of bubbles,
…)"*. That is the whole of it. There is **no** Davidson & Harrison entry in their
literature cited (journal page 492), where `Davidson, J. F.` appears only as a
co-author of Orcutt, Davidson & Pigford (1958); equations 3 and 9 and the `0.711`
carry no citation at all; and every derivation is deferred to a companion paper
(Kunii & Levenspiel, *I&EC Fundamentals* **7** 446, 1968b) that is not on disk.

The monograph usually named as the origin of these results was **not consulted**
and is the source of nothing here, and even the knowledge that it is the usual
citation is general background rather than something traced to a document read for
this page, so its bibliographic detail is deliberately absent from the metadata.
**The derivation on the page is done here from scratch.**

The slug is `davidson-bubble`, not `davidson-harrison-bubble`: nothing read here
attributes the model to two authors.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | `DavidsonBubble`, the pymrm solve of both potentials, and the trajectory route to the cloud radius |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata |

**No `data/` directory.** The printed appendix numbers this page checks against
are already transcribed, sidecar and all, in
`pages/E2.1-kunii-levenspiel-bubbling-bed/data/`, and are loaded **cross-page**
with `load_data(..., page="E2.1-kunii-levenspiel-bubbling-bed")`. One
transcription in the repository, no chance of two copies drifting. Every row used
was re-read here on a 600 dpi render before use. **This page must be integrated
while `E2.1` remains published.**

## The model in one paragraph

In the frame of the bubble, the solids stream down past it at `u_br` in inviscid
potential flow, and the interstitial gas percolates through those solids at
`u_f = u_mf/eps_mf` by Darcy's law. The bubble contains gas of negligible
viscosity and weight, so its surface is at **uniform pressure**. Both the solids
potential and the gas potential are then harmonic with the same `cos(theta)`
angular dependence, so both reduce to the **same 1D radial ODE**

```
(1/r^nu) d/dr( r^nu dh/dr ) - nu h / r^2 = 0
```

with `nu = 2` for a spherical bubble and `nu = 1` for the 2-D one. One operator,
two boundary conditions at `r = R`: `h = 0` for the gas (uniform pressure),
`dh/dr = 0` for the solids (nothing crosses). The problem is linear: one sparse
solve per potential.

## Validation

**Route 1: worked examples with printed intermediates.** All three of Kunii and
Levenspiel's appendices evaluate one of these quantities with the arithmetic
shown. Appendix C's `gamma_c` line is a fully substituted evaluation of
**equation 9** with a printed answer. **But these values are already published on
`E2.1`, from the one shared transcription this page loads cross-page.** They are
the same reading checked twice, not independent corroboration.

**Route 2: an internal identity.** The through-flow computed from the numerical
field **rebuilds the `4.5`** that opens equations 10 and 11. That constant is never
an input; it falls out of `-f'(R) = 3 u_f` and a hemisphere integral. The closed
form derived on the page gives the same number; the solve shows that the discrete
operator reproduces it, to its discretisation error.

| check | result |
| --- | --- |
| through-flow coefficient vs 4.5, n = 1600 | 1.1e-6 <!-- claim: throughflow_coefficient_rel --> |
| cloud radius vs the closed form, n = 1600, fraction of R | 7.2e-8 <!-- claim: Rc_rel_error_n1600 --> |
| observed order of the field error | within [1.95, 2.05] <!-- claim: grid_order_observed --> |
| far field: R_c moved by extending the domain from 20 R to 400 R | below 1e-9 <!-- claim: far_field_robin_shift_Rc --> |
| cloud radius, zero of W vs trajectory bisection | below 1e-10 <!-- claim: Rc_two_routes_rel --> |
| R_c/R limits, u_br/u_f from 1.001 to 1000 | 5.3e-7 <!-- claim: limit_worst_rel --> |
| 2-D cloud vs its own closed form | 3.9e-8 <!-- claim: cloud_2d_rel_error --> |
| `K_bc`, `H_bc` rebuilt vs printed (E2.1's numbers), worst | 0.43 <!-- claim: Kbc_Hbc_worst_pct --> per cent |
| printed appendix values (E2.1's numbers), worst | 1.06 <!-- claim: appendix_worst_pct --> per cent |

**Two reported numbers are weaker than they look, and the break table says how
much.** The discrete integral identity is an *algebraic* identity, the
volume-weighted sum of the discrete residuals, and it telescopes for any
consistent operator with any boundary data. It catches a geometry index in
`construct_div` or a sink coefficient that differs from the model's `nu`, and
nothing else: it stays at round-off under the four boundary-condition defects,
including a mis-signed outer condition that ruins the field. `R_c/R` is blind to
the two far-field defects that rescale both potentials together. **If you copy
this skeleton, watch the through-flow line, not the conservation line.**

## The finding worth carrying forward

`V_c/V_b = 3/(u_br/u_f - 1)` has a **pole**, not a zero. Below `u_br = u_f` there
is no cloud at all and equation 9 returns a negative number.

**Appendix B of the source paper is on the wrong side of that threshold.** Its
printed `d_b`, `u_mf` and `eps_mf` give `u_br/u_f` = 0.787 <!-- claim: ubr_over_uf_B -->
and `V_c/V_b` = -14.1 <!-- claim: Vc_over_Vb_eq9_B -->. Nothing in the paper is
wrong (appendix B is a heat-transfer example that needs `gamma_b` and `H_bc` and
never evaluates equation 9), but the bubble there would have to exceed
0.81 <!-- claim: db_threshold_B --> cm for a cloud to exist, and the paper says
equation 9 is for "beds with fast rising bubbles" without saying how fast is fast.
The verdict is robust to the one appendix-B number in doubt: `E2.1` records that
its printed `(1-eps_f)(u_0 + 5.8)` implies `eps_mf` = 0.447 <!-- claim: emf_implied_B -->
rather than the stated value, which moves `u_br/u_f` down to
0.704 <!-- claim: ubr_over_uf_B_implied -->, and nothing below `eps_mf` =
0.635 <!-- claim: emf_needed_B --> gives that bed a cloud.

The through-flow has no such threshold: `u_br` never enters the gas potential, so
the `4.5` holds on both sides of it while the `3` does not.

## The reusable trick

Only two radial modes exist, `h = A r` and `h = B r^-nu`, and the combination

```
dh/dr + (nu/r) h = (nu+1) A
```

is satisfied **exactly by both**. Imposing it as a Robin condition at a finite
outer radius pins the far field without truncating the decaying mode. What moving
the outer radius still changes is the discretisation of the outer cells, and the
page measures it.

## What this page is not

It is not `E2.1`. That page takes `u_br`, `V_c/V_b` and `K_bc` as **inputs** and
builds the bubbling-bed reactor model (`S7`) on them; it never solves a flow
field. This page is the `S3` hydrodynamic result those inputs come from.

## Corrections

- 2026-11-02: said the `R_c/R` limits hold "over seven decades" of `u_br/u_f`.
  The seven values run from 1.001 to 1000: three decades of `u_br/u_f`, six of
  `u_br/u_f - 1`.
- 2026-11-02: gave the through-flow coefficient as "1.1e-6". That number is its
  relative error against 4.5, not the coefficient, which is 4.5 to that error
  <!-- claim: throughflow_coefficient_rel -->.
- 2026-11-02: called the 1.2e-14 agreement between the cloud radius from the zero
  of `W` and from trajectory bisection "machine precision". Both routes
  interpolated `W` with the same cubic spline, and 45 halvings of the bracket
  resolve the radius to about 3e-14 of it, so the number measured the bisection.
  With the zero now read from the face values and the trajectories run through
  splines, the gap is below 1e-10 <!-- claim: Rc_two_routes_rel -->.
- 2026-11-02: said the integral identity "reads better than baseline under a
  wrong-but-consistent `nu` that makes the coefficient 77 % wrong". That row
  solved the correct two-dimensional problem and divided its through-flow per unit
  length by a sphere's volume, a dimensionally inconsistent ratio, and evaluated
  the identity with the same `nu`. A two-dimensional operator inside the
  three-dimensional model moves the coefficient, and the identity, evaluated with
  the model's `nu`, catches it.
- 2026-11-02: said the through-flow coefficient is "the only one of these numbers
  that catches every defect". The field error against the closed form also moved
  on every defect of the table; the coefficient is the one that needs no closed
  form.
