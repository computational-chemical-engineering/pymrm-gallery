# A4.2: Maxwell–Stefan vs Fick for multicomponent mixtures

Reproduces the bulk-fluid core of Krishna & Wesselingh (1997), *The
Maxwell–Stefan approach to mass transfer*, Chem. Eng. Sci. 52(6) 861–911
([doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7)):
the demonstration that the Maxwell–Stefan formulation and the generalized Fick
matrix `[D] = [B]^-1 [Gamma]` are two bookkeepings of the same physics, and
where the scalar effective-diffusivity simplification breaks (osmotic
diffusion, reverse diffusion, diffusion barrier).

## Headline

- **The printed worked examples, recomputed.** The ideal-gas `[D]` to
  1.03 % <!-- claim: ex1_fick_matrix_maxdev --> per element and `[B]` to
  3.53 % <!-- claim: ex1_b_maxdev --> (its most heavily rounded element); the
  non-ideal triple closes on `[D] = [B]^-1 [Gamma]` to
  1.65 % <!-- claim: ex2_identity_maxdev --> once the sign of the printed `B12`
  is corrected. As printed it misses by 75.1 % <!-- claim: ex2_as_printed_maxdev -->,
  and the rival reading (the typo in `B21`) by
  301.6 % <!-- claim: ex2_rival_b21_maxdev -->. A back-solve from `B11`, `B21`,
  `B22` alone returns the authors' unrounded diffusivities
  3.375e-9 <!-- claim: d12_backsolved -->, 2.507e-9 <!-- claim: d13_backsolved -->
  and 1.661e-9 <!-- claim: d23_backsolved --> m²/s and predicts
  `B12` = +0.0359e9 <!-- claim: b12_backsolved -->: the printed magnitude to
  0.21 % <!-- claim: b12_backsolve_magnitude_dev -->, the opposite sign.
- **One physics, two bookkeepings, tested by a route that could fail.** The raw
  friction balances with the `sum(J) = 0` bootstrap, sharing no code with the
  `[B]` route, reproduce `[D](x)` over a fixed random sample of compositions to at most
  1e-9 <!-- claim: fick_matrix_independent_maxdev --> and the whole two-bulb
  transient to at most 1e-10 <!-- claim: ms_independent_equivalence --> in mole
  fraction; Maxwell–Stefan is invariant to the eliminated species to at most
  1e-10 <!-- claim: ms_elimination_spread -->.
- **Against the Duncan–Toor measurements** (the `A4.9` dataset, 28
  <!-- claim: n_exp_points --> points): mean absolute deviation
  0.59 % <!-- claim: mean_abs_dev_ms --> (mole %) for Maxwell–Stefan and
  3.14 % <!-- claim: mean_abs_dev_wilke --> for the Wilke scalar law with CO₂
  eliminated. Over the 22 <!-- claim: n_paired_points --> points after t = 0 (the
  t = 0 points are the initial condition, paired difference exactly zero) the
  paired difference is 3.25 % <!-- claim: mean_abs_dev_wilke_minus_ms -->
  with a standard error of 0.58 % <!-- claim: mean_abs_dev_wilke_minus_ms_se -->
  (points taken as independent). On the nitrogen bulb difference (5
  <!-- claim: n_n2_points --> paired times): 1.38 % <!-- claim: n2_diff_err_ms -->
  against 9.07 % <!-- claim: n2_diff_err_wilke -->; over the 4
  <!-- claim: n_paired_n2_points --> times after t = 0 the paired difference is
  9.61 % <!-- claim: n2_diff_err_wilke_minus_ms --> with standard error
  1.50 % <!-- claim: n2_diff_err_wilke_minus_ms_se -->. Nothing fitted.
- **The scalar law is a bookkeeping choice.** Its mean deviation is
  1.44 % <!-- claim: mean_abs_dev_wilke_elim_n2 --> with N₂ eliminated and
  4.31 % <!-- claim: mean_abs_dev_wilke_elim_h2 --> with H₂ eliminated. Applied
  to all three species its fluxes miss `sum(J) = 0` by
  104.8 % <!-- claim: wilke_closure_defect --> of the largest flux.
- **The anomalies, root-found in time.** Point O (zero N₂ difference) at
  0.959 <!-- claim: n2_point_o_time_h --> h, the diffusion barrier (zero N₂
  flux) at 7.504 <!-- claim: n2_barrier_time_h --> h with a bulb difference of
  0.1231 <!-- claim: n2_barrier_dx -->, extrapolated to zero step. The review's frozen-matrix estimate of
  the initial flux ratio is off by 6.7 % <!-- claim: frozen_matrix_ratio_error -->.

The mean deviation against the measurements is the number page `A4.9` already
publishes (same solver, grid and dataset), so it confirms a faithful port, not
an independent corroboration.

## Numerics

The capillary is quasi-steady; the bulbs are marched by backward Euler. Grid
and time errors are measured apart: the grid order on the steady capillary
alone, 2.03 <!-- claim: grid_order --> at worst on five nested grids, with
`J2/J1` on the production grid
3.3e-5 <!-- claim: grid_error_ratio_nz40 --> (relative) from its Richardson
limit; the time order on bulb-1 N₂ at step times common to a five-level step
ladder, 0.92 <!-- claim: time_order --> at worst, approaching 1 from below;
halving the production step moves bulb-1 N₂ by
1.48e-4 <!-- claim: timestep_sensitivity -->. Both are far inside the
digitisation error of the data.

## What pymrm adds

The capillary solve: `construct_grad` with `shapes_d` for the moving bulb
compositions, `construct_boundary_value_matrices` for the end-face compositions,
`construct_div`, `NumJac`, `newton`, and `pymrm.checks` for the refinement
ladders and the roots in time. The solver is the `A4.9` page's; the page's own
contribution sits at the formulation level (the friction-system route, the
elimination test, the structure of the scalar law's failure, the B12 sign).

## Files

- `index.py` / `index.ipynb`: the page (paired by jupytext).
- `model.py`: the flux closures, the friction-system route, the operators and
  the capillary and bulb solvers.
- `test_model.py`: one Maxwell–Stefan transient against the stored deviation.
- `data/krishna-wesselingh-1997-worked-examples.csv`: the printed values of the
  two worked examples, transcribed from high-resolution renders (tier 6, the authors'
  own computations; provenance in the sidecar).
- The experimental comparison loads the digitised Duncan & Toor dataset
  published with page `A4.9` cross-page; nothing was re-digitised.

## Corrections

- 2026-11-10: the N₂ difference "peaks +0.1227 at t = 7.52 h" was the largest step
  value of the march; extrapolated to zero step the N₂ flux vanishes at
  7.504 <!-- claim: n2_barrier_time_h --> h, where the difference is
  0.1231 <!-- claim: n2_barrier_dx -->.
- 2026-11-10: "D1_eff (H2) range        : 6.43e-5 ... 7.55e-5 m2/s - positive and finite throughout, straying at most 5.4%"
  divided each step's flux by the previous step's bulb gradient. With both at the
  same time and extrapolated to zero step, hydrogen's effective diffusivity is
  4.9 % <!-- claim: h2_eff_excursion --> below its lower pair value at the end of the
  run, where the excursion is still growing.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```
