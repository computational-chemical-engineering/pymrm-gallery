# B1.6 — The Prater relation

Inside a catalyst particle the temperature is not an independent field. Reactant
is consumed and heat is released at the same points in fixed proportion, so

    lambda (T - T_s) = (-dH) D_e (c_s - c)

exactly — no rate constant, no activation energy, no particle shape. This page
proves it, confirms it numerically, establishes what that confirmation is and is
not sensitive to by injecting defects on purpose, and then breaks the relation
itself to find out what its assumptions are worth.

- **Structures:** `S3` (1D steady BVP)
- **Sources:** Weisz & Hicks (1962) Chem Eng Sci 17(4) 265–275,
  doi:10.1016/0009-2509(62)85005-2 — **the paper actually read**, eqs. (4), (5),
  (7) and (8) transcribed from a 600 dpi render of journal page 266.
  Prater (1958) Chem Eng Sci 8(3–4) 284–286,
  doi:10.1016/0009-2509(58)85035-6 — **the origin of the result, not consulted**;
  no copy is on disk and none is open access, so it is cited as Weisz & Hicks
  cite it.
- **Runtime:** ~1 min

## Results

The identity is proved symbolically for an arbitrary rate function and an
arbitrary geometry index, then measured as
`eps = max|(theta-1) - beta(1-y)| / |beta|` on solutions of the **two-field**
system, in which the relation is never substituted.

| check | result | can it fail? |
|---|---|---|
| 324 pymrm solves (3 geometries × 3 rate laws × gamma 10–30 × beta −0.4…+0.6) | worst `eps` = 1.1e-11, median 8.9e-13 | only for one class of defect — see below |
| Newton residual over those 324 solves | worst 4.6e-11 | yes; this is what says they converged |
| independent `scipy.solve_bvp` collocation, no pymrm | worst `eps` = 2.6e-13 | no — it inherits the invariant structurally |
| pymrm vs the collocation **profiles**, y(u) | 1.0e-06 at n_u = 800, second order | **yes** |
| vs `B1.1`'s shooting reference, which discretises nothing | 3.9e-05 at n_u = 800, order 2.11 | **yes — the page's real discretisation test** |
| beta recovered from the solution vs `B1.5`'s beta = 0.6 | 3.3e-10 | no — the identity restated; a convention check |
| two-field vs the reduced equation `B1.1`/`B1.5` solves | 6.1e-12 in y | no — guaranteed by the identity; an implementation check |

**What `eps` is worth.** The identity is linear and one `construct_div` serves
both fields, so `w = theta + beta*y` is a closed linear subsystem: `eps` is a
floating-point measurement, and its size carries no information about the
discretisation or the solve. The page establishes that by injecting defects
rather than by asserting it:

| injected defect | `eps` |
|---|---|
| none (baseline, sphere, phi=2, beta=0.3, gamma=20) | 2.1e-13 |
| heat source uses `-beta` | 4.2e-01 |
| heat source uses `1.01*beta` (1 % scale error) | 1.0e-02 |
| y Dirichlet but theta Robin (boundary-condition type mismatch) | 1.6e-01 |
| rate evaluated at theta = 1 in the energy equation only | 5.4e-01 |

| what it is blind to | consequence | `eps` |
|---|---|---|
| Newton stopped at `maxfev=1` | y = 2.12 at the centre — impossible | 8.5e-12 |
| the wrong geometry index | eta 57 % wrong | 4.7e-12 |
| `n_u = 3` | eta 37 % wrong | 2.4e-15 |

So `eps` detects an inconsistency between the two source terms or between the
two boundary conditions — linearly in the size of the defect — and nothing else.
The severity of the sweep (it reaches eta = 162 and depletion to y = 1.6e-195)
makes the *physics* varied but adds no weight to `eps`, which one solve would
have given just as well.

## Where it breaks

**A film with unequal Biot numbers.** The interior relation survives, referenced
to the *surface*. The bulk-referenced form — the one usually written down — is
wrong by `beta (1-y_s) (Bi_m/Bi_h - 1)`, uniformly in position. That closed form
is derived on the page. The numerics match it to 1.9e-11 relative, but that
match is forced — both steps of the derivation hold in the discrete system, and
an 8-cell grid matches it to 5.6e-14 — so it confirms the algebra, not the
solve. The content is the derivation and the magnitudes:

| Bi_m | Bi_h | surface-referenced eps | bulk-referenced eps | closed form | share of dT outside the pellet |
|---|---|---|---|---|---|
| 100 | 100 | 2.1e-12 | 2.0e-12 | 0 | 2.8% |
| 100 | 30 | 1.3e-12 | 5.854e-02 | 5.854e-02 | 8.8% |
| 100 | 10 | 3.9e-12 | 2.580e-01 | 2.580e-01 | 24.2% |
| 100 | 3 | 5.5e-12 | 2.881e+01 | 2.881e+01 | 99.6% |
| 10 | 100 | 3.3e-12 | 1.727e-01 | 1.727e-01 | 2.7% |

**A transient.** The relation is a linear invariant of the time-dependent system
exactly when `Le = 1`, and only then. The control holds to 3.5e-13 over the whole
transient; away from it the violation reaches 0.124 of beta at Le = 0.1 and 0.252
at Le = 10, peaking early and decaying to zero at steady state.

## Data

**Provenance tier 6 — not experimental, and this is a proof plus numerical
confirmation, not validation against measurement.** Neither source reports
measurements: Weisz & Hicks (1962) is computational and Prater (1958) is a
theoretical note. No dataset is shipped and nothing is digitised.

## A note on the source

Weisz & Hicks print two expressions for beta in their eq. (8) that carry
**opposite signs**: `c_0 H D / (K T_0)` and `(dT/T_0)_max`, which eq. (7) at
c → 0 makes negatives of one another. The page prints both alternatives and shows
which one the paper's own rate expression, its stated sign convention (beta from
0 to +0.8 exothermic) and its eta > 1 result all select, rather than repairing it
by inference.

The scan's text layer is unusable for equations — eq. (7) OCRs as
`AT = T - To = - F` and `[grad ~1,~~`. Render journal page 266 at 600 dpi.
