# J1.5: The linear driving force

Almost every adsorption, ion-exchange and chromatography model replaces
diffusion inside the particle with a single first-order relaxation at
`k = 15 D/r²`. This page solves the sphere that approximation stands in for.

- **Structures:** `S3` (1-D BVP), `S4` (1-D transient PDE)
- **Reference:** Glueckauf (1955), Trans. Faraday Soc. 51 1540–1551,
  doi:10.1039/TF9555101540
- **Runtime:** ~12 s

## Agreement

pymrm's sphere matches the Crank series (summed to a measured relative tail
below 1e-15 <!-- claim: series_tail_bound -->) to **2.6e-3** <!-- claim: solver_worst_rel_err -->
at worst across `τ = 0.001` to `0.5`, and to 6.5e-5 <!-- claim: solver_rel_err_n400 -->
at `τ = 0.05` on 400 cells. At the worst time 97 % <!-- claim: time_error_share_tau1e3 -->
of the error is the backward-Euler step, not the grid (the spatial part is 7.5e-5 <!-- claim: spatial_err_tau1e3 -->);
refinement of the grid alone converges at second order (the observed order departs from 2 by
below 0.2 <!-- claim: grid_order_worst_dev -->) and of the step at first order (it departs from 1
by below 0.1 <!-- claim: time_order_worst_dev -->).

## Where 15 comes from

Assume a parabolic profile `q = a + b(r'/r)²`. Then `q* − q̄ = (2/5)b` while the
surface flux gives `dq̄/dt = 6Db/r²`; eliminating `b` leaves exactly `15 D/r²`.
In a geometry with index `ν` the same argument gives `(ν+1)(ν+3)`: 15, 8 and 3.

## What the page shows

The coefficient the LDF *would need* at each instant,

```
k_eff(τ) = (dq̄/dτ) / (q* − q̄)
```

computed in closed form from the series. It is not flat:

| regime | behaviour |
|---|---|
| `τ → 0` | diverges as `3/√(πτ)`: uptake grows as a square root, and no exponential can follow that |
| `τ = 0.0217` <!-- claim: k_eff_equals_15_tau --> | passes through **15**, particle 43 % <!-- claim: uptake_at_k_eff_15 --> loaded |
| `τ → ∞` | settles on **π² = 9.87** <!-- claim: long_time_decay_constant --> |

15 is exact at one instant in the sense that the surface gradient equals that of
the parabola with the same mean. It does not mean the profile is parabolic there:
that parabola's gap to the real profile, largest at the centre, is 0.42 <!-- claim: profile_gap_at_k_eff_15 -->
while the real centre is still empty. The LDF with 15 under-predicts uptake from
an empty particle by 0.170 <!-- claim: ldf_worst_abs_err --> at worst, at `τ = 0.0118` <!-- claim: ldf_worst_abs_err_tau -->,
crosses the exact curve at `τ = 0.094` <!-- claim: ldf_crossing_tau -->, and
over-predicts by 0.036 <!-- claim: ldf_over_prediction_max --> at worst.

## What pymrm adds

Nothing to the answer above, which comes from the closed form. The pymrm solve is
the inner problem of a coupled column, which has no closed form once the isotherm is
nonlinear or the diffusivity varies; the page measures its
discretisation (second order in space, first order in time, an exact flux balance
below 1e-12 <!-- claim: flux_balance_gap -->, a geometry index `nu`).

## Honest limits

Provenance tier 6: the reference is the analytic series, not measurement.
Glueckauf's Tables 3–5 are not transcribed (the 1955 scan drops decimal points
wholesale, rendering 0.8647 as `8647`), and nothing is lost, since the exact
series is a stronger reference than his tabulated approximation to it.

The page treats a step in surface concentration with a linear isotherm. The
coupled case, where the surface value tracks a bulk the particle is itself
depleting, is the `S8` problem of `D1.4` and is not solved here.

## Corrections

- 2026-11-03: said the volume weights "sum to one". They were the midpoint rule
  `3 r² Δr`, which does not (the page's break row for `weights_sum_gap` prints its
  gap). The weights are now the exact cell volumes and sum to one to round-off.
  With the exact volumes the 400-cell error at `τ = 0.05` is 6.5e-5 <!-- claim: solver_rel_err_n400 -->, not 6.6e-5.
- 2026-11-03: printed the LDF crossing as `τ = 0.0912`, `k_eff = 15` at `τ = 0.0219`
  and the worst under-prediction at `τ = 0.0119`, and said the particle is 44 % loaded
  there. Those were nodes of log-spaced grids. The roots are `τ = 0.0940`
  <!-- claim: ldf_crossing_tau -->, `0.0217` <!-- claim: k_eff_equals_15_tau --> (43 %
  <!-- claim: uptake_at_k_eff_15 --> loaded) and `0.0118` <!-- claim: ldf_worst_abs_err_tau -->.
- 2026-11-03: said 15 is "the best a constant can do near there because that is where
  the profile really is parabolic". At the instant `k_eff = 15` the profile is far
  from parabolic: the largest gap between it and the parabola with the same mean is
  0.42 <!-- claim: profile_gap_at_k_eff_15 --> at the centre, which is empty in the real
  profile. `k_eff = 15` means the surface gradient equals that parabola's.
- 2026-11-03: gave "max relative error 95.5 % for `τ < 0.01`" as an error range. It was
  the value at `τ = 10⁻⁴`, the first point of the scan; the shortfall tends to 100 %
  as `τ → 0`.
