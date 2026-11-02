# B1.6: The Prater relation

Inside a catalyst particle the temperature is not an independent field. Reactant
is consumed and heat is released at the same points in fixed proportion, so

    lambda (T - T_s) = (-dH) D_e (c_s - c)

exactly: no rate constant, no activation energy, no particle shape. This page
proves it, confirms it numerically on two-field solves, establishes what that
confirmation is and is not sensitive to by injecting defects on purpose, and
then breaks the relation itself to find out what its assumptions are worth.

- **Structures:** `S3` (1D steady BVP)
- **Sources:** Weisz & Hicks (1962) Chem Eng Sci 17(4) 265-275,
  doi:10.1016/0009-2509(62)85005-2, **the paper actually read**: eqs. (4), (5),
  (7) and (8) transcribed from a 600 dpi render of journal page 266.
  Prater (1958) Chem Eng Sci 8(3-4) 284-286,
  doi:10.1016/0009-2509(58)85035-6, **the origin of the result, not consulted**;
  no copy is on disk and none is open access, so it is cited as Weisz & Hicks
  cite it.
- **Runtime:** about 95 s

## Results

The identity is proved symbolically for an arbitrary rate function and an
arbitrary geometry index, then measured as
`eps = max|(theta-1) - beta(1-y)| / |beta|` on solutions of the **two-field**
system, in which the relation is never substituted.

| check | result | can it fail? |
|---|---|---|
| 324 pymrm solves (3 geometries x 3 rate laws x gamma 10-30 x beta -0.4 to +0.6) | worst `eps` below 1e-10 <!-- claim: sweep_worst_identity_violation --> | only for one class of defect, see below |
| backward error of those solves | below 1e-12 <!-- claim: sweep_worst_backward_error --> | yes; this is what says they converged |
| pymrm vs an independent `scipy.solve_bvp` collocation, y(u) | 1.0e-6 <!-- claim: pymrm_vs_collocation_profile_n800 --> at 800 cells, second order | **yes** |
| pymrm eta vs the shooting reference of `B1.1`, which discretises nothing, at the same phi | 2.0e-6 <!-- claim: vs_shooting_reference_n800 --> at 800 cells, second order | **yes** |
| beta recovered from the solution vs `B1.5`'s beta = 0.6 | below 1e-8 <!-- claim: beta_recovered_vs_0p6 --> | no: the identity restated; a convention check |
| two-field vs the reduced equation `B1.1`/`B1.5` solves | below 1e-10 <!-- claim: two_field_vs_reduced_max_dy --> in y | no: guaranteed by the identity; an implementation check |

**What `eps` is worth.** The identity is linear and one `construct_div` serves
both fields, so `w = theta + beta*y` is a closed linear subsystem: `eps` is a
floating-point measurement, and its size carries no information about the
discretisation or the solve. The page establishes that by injecting defects on
one case of the sweep (sphere, phi = 2, beta = 0.3, gamma = 20) and reading
three numbers each time: `eps`, the effectiveness factor and the backward
error. A sign error, a scale error and a rate evaluated at the wrong state in
the heat source, and a film on the temperature field only, lift `eps` from
round-off to the size of the defect; a heat source scaled by 1 + delta gives
`eps = delta (1 - y_min)` exactly. A Newton iterate stopped after one step (with
more reactant at the centre than at the surface), the wrong geometry index and a
3-cell grid leave `eps` at round-off (below 1e-10 <!-- claim: blind_rows_max_identity_violation -->
for all three); the effectiveness factor and the backward error catch those. The break table in the notebook prints every row.

**Which steady state.** Where several steady states exist, the solve lands on
the least-ignited one, the state a start-up from surface conditions reaches.
27 <!-- claim: sweep_cases_several_states --> of the sweep cases with beta > 0
have a second steady state, and the solve is on the less ignited one in
27 <!-- claim: sweep_several_states_page_lowest --> of them. At phi = 0.4 (sphere,
beta = 0.6, gamma = 20) that is eta
1.159 <!-- claim: phi04_eta_page --> against
44.5 <!-- claim: phi04_eta_upper --> on the upper branch, and the state equals a
resolved Le = 1 start-up to below 1e-8 <!-- claim: phi04_state_vs_startup -->.

## Where it breaks

**A film with unequal Biot numbers.** The interior relation survives, referenced
to the *surface*. The bulk-referenced form, the one usually written down, is
wrong by `beta (1-y_s) (Bi_m/Bi_h - 1)`, uniformly in position. That closed form
is derived on the page. The numerics match it to round-off, but that match is
forced (both steps of the derivation hold in the discrete system, and an 8-cell
grid matches it as well), so it confirms the algebra, not the solve. The content
is the derivation and the magnitudes: at Bi_m = 100, Bi_h = 10 the bulk form is
off by 0.258 <!-- claim: film_bulk_eps_bim100_bih10 --> of beta and
24.2 % <!-- claim: film_external_share_bih10 --> of the temperature rise sits in
the film; at Bi_h = 3 the pellet ignites and
99.6 % <!-- claim: film_external_share_bih3 --> of it does.

**A transient.** The relation is a linear invariant of the time-dependent system
exactly when `Le = 1`, and only then. The control stays at round-off over the
whole transient; away from it the violation peaks at
0.252 <!-- claim: transient_worst_violation --> of beta at Le = 10,
0.124 <!-- claim: transient_peak_Le0p1 --> at Le = 0.1 and
0.146 <!-- claim: transient_peak_Le0p01 --> at Le = 0.01, early in the
start-up, and decays to zero at steady state. The peaks come from variable-step
BDF2 taken to its limit in the step, checked against `solve_ivp` on the same
semi-discrete system.

## Data

**Provenance tier 6: not experimental, and this is a proof plus numerical
confirmation, not validation against measurement.** Neither source reports
measurements: Weisz & Hicks (1962) is computational and Prater (1958) is a
theoretical note. No dataset is shipped and nothing is digitised.

## A note on the source

Weisz & Hicks print two expressions for beta in their eq. (8) that carry
**opposite signs**: `c_0 H D / (K T_0)` and `(dT/T_0)_max`, which eq. (7) at
c -> 0 makes negatives of one another. The page prints both alternatives and shows
which one the paper's own rate expression, its stated sign convention (beta from
0 to +0.8 exothermic) and its eta > 1 result all select, rather than repairing it
by inference.

The scan's text layer is unusable for equations: eq. (7) OCRs as
`AT = T - To = - F` and `[grad ~1,~~`. Render journal page 266 at 600 dpi.

## Corrections

- 2026-11-01: said pymrm's effectiveness factor agrees with `B1.1`'s shooting
  reference to 3.9e-05 at n_u = 800, at observed order 2.11, compared "at equal
  centre concentration, ... so no interpolation enters it". The pymrm value used as
  the centre concentration was the first cell's, half a cell from the centre. On
  the ignited state (phi = 0.8, a nearly depleted centre) eta is so sensitive to
  it that this offset set the number, and the reference's absolute tolerance
  (1e-12) was the size of the centre value itself. Compared at the same phi, the
  agreement is 2.0e-6 <!-- claim: vs_shooting_reference_n800 --> at 800 cells,
  at second order.
- 2026-11-01: said the sweep reaches "depletion to y = 1.6e-195" (its notebook
  printed 4.27e-230 for the same quantity). A centre concentration that far
  below the surface value is below what the solve resolves: Newton judges each
  field against its largest value, which is 1. The page no longer quotes it.

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        the two-field pellet, the transient march, the reduced equation,
                the collocation and shooting routes
test_model.py   fast checks of model.py and one stored claim
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
```
