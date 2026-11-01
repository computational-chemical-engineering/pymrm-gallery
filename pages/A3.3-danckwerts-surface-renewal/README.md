# A3.3: Danckwerts surface renewal

Replace the stagnant interfacial film with a surface that is continually torn up
and replaced, assume the chance of an element being replaced does not depend on
its age, and the age distribution is forced: `phi(theta) = s exp(-s theta)`.
Averaging the penetration flux over it gives `k_L = sqrt(Ds)`, against film
theory's `k_L ~ D`.

**Source.** P. V. Danckwerts, *Significance of liquid-film coefficients in gas
absorption*, Ind. Eng. Chem. **43**(6) 1460-1467 (1951),
[doi:10.1021/ie50498a055](https://doi.org/10.1021/ie50498a055). Read from the ACS
PDF at its **native 300 ppi** (`pdfimages -list`: CCITT-G4 bilevel, 300 ppi;
rendering higher is interpolation), with the two constants the text layer
mangles re-read on nearest-neighbour crops at that resolution. The PDF is not in
this repository.

## How it is solved

`Element` (`model.py`) is one surface element with an optional first-order
reaction: the diffusion operator from `construct_grad` and `construct_div`,
assembled once and kept in banded form, backward Euler on age steps uniform in
`sqrt(theta)` (one `solve_banded` call per step), Richardson-extrapolated in the
step. The surface flux comes from the discrete mass balance, with boundary fluxes
read by `compute_boundary_values`. `TwoSpeciesElement` solves the second-order
pair with one `newton` call per step and the analytic Jacobian of the source,
assembled with `construct_coefficient_matrix`; its results are extrapolated in
the step and in the grid. Refinement studies use `pymrm.checks.observed_orders`,
and the source Jacobian is checked with `pymrm.checks.check_jacobian`.

## What the page establishes

- **Every number the paper prints** reproduces, and three of them are loose
  enough to be worth stating rather than rounding away (the depth at 1 hour is
  6.83 <!-- claim: depth_1h_mm --> mm against his "about 6 mm").
- **The element against the closed forms**: eqs. 8, 10/30 and 12/34 with
  relative deviations -4.0e-5 <!-- claim: R_vs_eq8_rel -->,
  -1.6e-5 <!-- claim: R_vs_eq30_rel --> and -3.9e-5 <!-- claim: R_vs_eq34_rel -->,
  each being the grid error of 400 cells
  (3.9e-5 <!-- claim: space_err_n400 --> for eq. 8) plus the extrapolated step
  error. Orders from successive differences: 1.00 <!-- claim: time_order --> in
  the step, 2.00 <!-- claim: space_order --> in the grid.
- **A printed defect in eq. 25**, proved from the paper's own eqs. 9 and 29: the
  exponent needs `H^2`. With it, eq. 25 reproduces eq. 9 at every `(H, k_G/H)`
  swept to round-off; as printed it reproduces it nowhere (closest
  4.2 % <!-- claim: eq25_printed_best_dev_over_Hk --> off).
- **The case he says he could not solve**: the second-order reaction, solved
  directly, turns his 10 % bound into a measurement. At his illustrative point
  the substitution is in error by
  +0.210 <!-- claim: so_err_at_printed_50_pct --> %; it reaches 10 % at
  c0'/c* = 0.856 <!-- claim: so_10pct_point_c0p_over_cstar -->, a factor
  62 <!-- claim: eq36_conservative_factor --> below the 53 that eq. 36 demands.
- **A second route to eqs. 15 and 16**: the same solver at r'c0'/s = 3000 lands
  within 0.48 <!-- claim: inst_vs_eq15_worst_pct --> % of eq. 15, a deviation
  that shrinks as the rate grows.
- **Film theory and surface renewal as two ends of one calculation**:
  `k_L = sqrt(Ds) coth(d sqrt(s/D))`, whose exponent `n = 1/2 + z/sinh 2z` runs
  from 1 to 1/2 as the layer deepens.
- **Why the sqrt(D) exponent cannot be tested from this paper**, stated as the
  page's main negative result.

## What it does not do

The paper has no data. `s` comes from inverting the very relation one would want
to test (an inference, flagged as one, since the paper says only "calculated
from published values of `k_L`"), so film theory and surface renewal fit its one
physical number equally well, with zero residual each. Nothing here is validated
against a measurement. Whitman (`A3.1`) and Higbie (`A3.2`) are separate cases;
`A3.2` takes the three-way comparison this page declines. That does not touch
this paper's own finding: the `sqrt(D)` exponent is not testable from anything
Danckwerts prints.

## Corrections

- 2026-10-31: gave the error of the pseudo-first-order substitution at
  c0'/c* = 50 as +0.22 %. That value was extrapolated in the step only, on 250
  cells, and carried the grid error (the known-answer case's
  +0.0104 <!-- claim: so_control_extrap_pct --> % after step
  extrapolation is that grid error, not step error). Extrapolated in both axes it
  is +0.210 <!-- claim: so_err_at_printed_50_pct --> %.
- 2026-10-31: said eq. 36 is conservative by a factor of 53 in the concentration
  ratio, with the exact 10 % point near 1. The 10 % point is
  c0'/c* = 0.856 <!-- claim: so_10pct_point_c0p_over_cstar -->, root-found, and
  the factor is 62 <!-- claim: eq36_conservative_factor --> against the threshold
  of 53 <!-- claim: eq36_threshold_c0p_over_cstar -->.
- 2026-10-31: said the substitution is accurate to better than 1 % over
  essentially the whole range where diffusion and reaction play roughly equal
  parts. At r'c0'/s = 1 it is 4.8 <!-- claim: so_err_at_c0p_2_pct --> % at c0'/c* = 2 and
  1.04 <!-- claim: so_err_at_c0p_10_pct --> % at 10; it is below 1 % from 25
  upward.
- 2026-10-31: said a Dirichlet base on the reacting element moves the rate by
  +3.8 % because the deep liquid must decay as c0 exp(-r theta). The +3.8 % came
  from reading the base outflow from the last cell centre, which is inconsistent
  with the operator's face gradient; read at the face, a Dirichlet base 45
  diffusion lengths down changes the rate by
  less than 1e-9 <!-- claim: dirichlet_base_shift_depth45 -->. It does matter on a shallow
  element: 3 diffusion lengths deep, the Dirichlet base changes the rate by
  -2.3e-3 <!-- claim: dirichlet_base_shift_depth3 --> against the no-flux base on
  the same grid.
- 2026-10-31: listed the agreement of the two flux readings (mass balance and
  boundary gradient) as a (conditionally) structural row that "cannot be broken",
  because backward Euler conserves mass. It is the conservation identity of the finite-volume scheme, its gap is
  the residual of the banded solve, and a cell-centre read or a balance without
  the reaction breaks it. The same inconsistent base read made the finite-layer
  comparison with coth z 4.7e-6; read at the face it is
  7.4e-6 <!-- claim: kL_coth_worst_rel -->, the grid and step error.
- 2026-10-31: listed "arithmetic instead of harmonic face mean" as a structural
  break row. With a uniform diffusivity it injects nothing; it is gone.

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        Element, TwoSpeciesElement, age_average, rate, rate_second_order,
                eq15, lam_eq16
test_model.py   fast checks of model.py and two stored claims
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
data/           the paper's printed numbers, with a provenance sidecar
```

Runtime about 200 s. No other page's dataset is loaded.
