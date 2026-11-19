# A4.9: Duncan-Toor ternary gas diffusion

Two bulbs of gas joined by a capillary. Nitrogen starts slightly apart,
closes that gap, then keeps going *against* its own gradient until it is twice as
far apart the other way, and finally stops dead at the moment its gradient is
largest. Fick's law can produce none of that.

- **Structures:** `S3` (1D steady BVP), `S9` (implicit multicomponent flux)
- **Reference:** Duncan & Toor (1962), AIChE J 8(1) 38-41, doi:10.1002/aic.690080112
- **Runtime:** about 100 s (`meta.yaml` `runtime_seconds`)

## Agreement

0.59 % <!-- claim: mean_abs_deviation --> mean absolute deviation in mole
fraction on the 28 <!-- claim: n_points --> digitised points, largest 1.79 % <!-- claim: max_abs_deviation -->,
nothing fitted.

Point O (osmotic diffusion) at 0.96 <!-- claim: point_o_time_h --> h and the
diffusion barrier at 7.50 <!-- claim: barrier_time_h --> h: roots of the N2
difference and of the N2 flux on the time-exact `solve_ivp` route. At the barrier
the N2 difference is 2.0 <!-- claim: n2_extreme_ratio --> times its initial value,
with the opposite sign.

## The comparison with the data is the only guard

Flux uniformity along the capillary (below 1e-12 <!-- claim: flux_nonuniformity -->)
and the closed-cell mole balance (below 1e-13 <!-- claim: mole_balance_error -->)
are one algebraic identity: the capillary solve drives `div(N)` to the Newton
tolerance with `construct_div(nu=0)`, which on a Cartesian grid is the flux
difference between faces, and a uniform flux makes the two bulb updates cancel.
Neither can see a wrong constant. The notebook's break table perturbs each of
the five transcribed constants by one digit: every row moves the data
comparison, the regime times or both, and none lifts the flux non-uniformity
above 1e-10 (5 <!-- claim: n_injections_flux_within_bound --> of 5). They are kept, labelled structural, because they catch what
nothing else here catches: an unconverged capillary solve and a wrong `nu`.

## Numerical error

Grid: refined on the time-exact route over `n_z` 10 to 160, observed orders
at least 1.8 <!-- claim: grid_observed_order --> and at most 2.2 <!-- claim: grid_observed_order -->, spatial error at the
production `n_z = 40` of 1.6e-6 <!-- claim: grid_error_nz40 --> in mole
fraction. Time: backward Euler over `n_steps` 100 to 800, observed orders
at least 0.9 <!-- claim: timestep_observed_order --> and at most 1.1 <!-- claim: timestep_observed_order -->, time error at the
production 400 steps of 3.0e-4 <!-- claim: timestep_error_ns400 --> against
`solve_ivp` on the same grid. Time dominates, 183 <!-- claim: time_over_grid_error_ratio -->
times the spatial error; the mean
model-vs-data deviation is 3.6e3 <!-- claim: dev_over_grid_error --> times the
spatial error.

## What pymrm adds

The discretisation: constant `construct_grad` and `construct_div` operators with
the bulb compositions entering through `shapes_d`, a `NumJac` Jacobian for the
implicit Maxwell-Stefan flux and a checked `newton` solve. Against Fick with
Wilke effective diffusivities, which cannot cross zero, the N2 difference is off
by 9.1 % <!-- claim: fick_n2_dev --> against 1.4 % <!-- claim: ms_n2_dev --> for
Maxwell-Stefan; against Toor's linearisation, resolving the composition
dependence of [B] changes the N2 difference by up to 5.0 % <!-- claim: toor_dev_frac -->
of its value at the barrier. That dependence is the whole gap: with every pair
diffusivity equal, so that [B] is constant, Toor's closed form and the
time-exact route agree to below 1e-10 <!-- claim: toor_gap_equal_d -->.

## Data

Digitised from Figure 2 by programmatic marker extraction, validated by closure,
monotonicity and the closed-cell exchange balance. See
`data/duncan-toor-1962-run1.meta.yaml` for the method and its limits.

All model parameters were checked against the paper text. The paper uses an
*effective* (L/A) measured from binary runs, not the geometric value; the
notebook prints both.

## Corrections

- 2026-11-19: "Point O (osmotic diffusion) at 0.96 h, diffusion barrier at 7.52 h."
  The barrier is at 7.50 <!-- claim: barrier_time_h --> h, the root of the N2
  flux on the time-exact route; the old value was the backward-Euler step with
  the largest N2 difference.
- 2026-11-19: "composition dependence of [B] along the capillary (5.5% of the peak N2 separation)"
  The gap is 5.0 % <!-- claim: toor_dev_frac --> of the N2 difference at the
  barrier, measured against the time-exact route; the old value compared Toor's
  closed form with the backward-Euler trajectory and so included its time error.
- 2026-11-19: "building a difference three times larger than it started with"
  The N2 difference at the barrier is 2.0 <!-- claim: n2_extreme_ratio --> times
  the initial one in size; the largest measured is at least 2.2 <!-- claim: exp_n2_diff_ratio -->
  and at most 2.3 <!-- claim: exp_n2_diff_ratio --> times.
