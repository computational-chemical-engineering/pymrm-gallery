# A2.5 — Where the axial dispersion coefficient comes from

`A2.1` takes the axial-dispersion reactor as given. `A2.3` derives a dispersion
coefficient for a tube. This page is the experimental one: a pulse of argon in
air, two response curves, and the two-point moment method that turns them into
`D_L`.

- **Structures:** `S4` (1-D transient advection-dispersion)
- **Reference:** Edwards & Richardson (1968), Chem. Eng. Sci. 23(2) 109-123,
  doi:10.1016/0009-2509(68)87056-3
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`TracerBed` and `TwoZoneBed`, the simulated two-detector experiment,
  the moment inversion and the two-zone shape diagnostic), tested by
  `test_model.py`. The correlations, the fits and the bootstrap stay in the
  notebook.
- **Runtime:** about 940 s.

## Agreement

**Experimental**, against the markers digitised from Figure 9.

| Check | Result |
|---|---|
| Eq. (18) vs the Figure 9 markers | 11.75 % <!-- claim: eq18_vs_fig9_mean_abs_dev --> mean absolute deviation, bias 0.003 % <!-- claim: eq18_vs_fig9_bias --> |
| null: Eq. (15), no radial-mixing term | 29.91 % <!-- claim: eq15_vs_fig9_mean_abs_dev --> with a -29.83 % <!-- claim: eq15_vs_fig9_bias --> bias; 30 <!-- claim: n_points_above_pe2 --> of the 55 <!-- claim: n_fig9_markers --> points lie above its ceiling of Pe = 2 |
| refit β to the digitised points (γ fixed) | 8.98 <!-- claim: beta_refit_gamma_fixed --> against the printed 9.7; block bootstrap 95 % interval 7.50 <!-- claim: beta_ci_gamma_fixed_lo --> to 9.94 <!-- claim: beta_ci_gamma_fixed_hi --> |
| refit both parameters, a *different* estimator | γ = 0.702 <!-- claim: gamma_refit_free -->, β = 8.58 <!-- claim: beta_refit_free --> |
| Eq. (18)'s constants from Eq. (17)'s | 1.28 % <!-- claim: eq18_coeff1_from_eq17_rel_dev --> and 0.31 % <!-- claim: eq18_coeff2_from_eq17_rel_dev --> |
| gas crossover Re, printed as 1.8 | 1.830 <!-- claim: gas_crossover_Re --> from printed values only, +1.7 % <!-- claim: gas_crossover_Re_rel_dev -->, not exact |
| Table 1 voidages vs 1 − bulk/material density | worst 0.0018 <!-- claim: voidage_density_resid_worst --> absolute over all eight rows |
| two-point moment inversion of a simulated experiment | 1.38e-5 <!-- claim: moment_recovery_rel_err_n3200 --> relative at n = 3200 (van Leer, limiter correction implicit, Crank-Nicolson otherwise), third order on successive grids |
| numerical dispersion predicted before the runs | `u h/2 + u²Δt/2` accounts for the excess to 0.42 % <!-- claim: numerical_dispersion_prediction_worst --> |
| two-zone model vs an independent closed form | 3.07e-5 <!-- claim: two_zone_vs_closed_form_rel_err --> at n = 2000, third order |
| deferred correction vs `newton` with `ConvFluxTVD`, same step equations | below 1e-7 <!-- claim: tvd_two_route_D_gap --> in the recovered `D_L` |

Deviations are `(model − measured)/measured` everywhere on the page. Every claim
has a row in the defect-injection table that fails its own check.

## The 11.75 % is a goodness of fit, not a prediction

γ and β were fitted by Edwards and Richardson to the very markers this page
digitises. What *is* independent:

- the **shape**. Eq. (15) has the right molecular branch and the right high-Re
  limit and still cannot exceed Pe = 2, where over half the data do.
- the **refit**. Recovering β from a fresh digitisation tests the extraction and
  the page-image transcription of Eq. (17) against each other. No difference from
  the printed 9.7 is claimed: each estimator's bootstrap interval is quoted
  beside its own estimate, and whether 9.7 falls inside the γ-fixed one is
  decided in the second decimal by the block size.
- the **location of the maximum**: Pe = 3.186 <!-- claim: eq18_pe_max --> at
  Re = 3.195 <!-- claim: eq18_re_at_pe_max -->, against the paper's own bound
  ("greater than 3"). Gunn (1993) puts a maximum at Re ≈ 4, but for his Figure 1's
  *combined* gas-phase set, not for E&R's alone.

## What the page adds

The measurement is an inverse problem, so the page simulates it. A pulse
propagated on a pymrm bed, sampled at two detectors and inverted with Eqs. (6)
and (7), returns the input `D_L`, and the numerical dispersion of the cheaper
schemes is **predicted before the runs**: `u h/2` for upwind plus `u²Δt/2` for
implicit Euler. At n = 400 the bare upwind scheme reports a dispersion
coefficient 4.5 <!-- claim: upwind_n400_over_true --> times the truth from
response curves that look perfectly smooth.

The fine-particle section turns the authors' own shape diagnostic round. They
report that for fine particles the computed second response curve peaks *later*
than the measured one and attribute it to channels forming. A two-zone bed
constrained to their measured `D_L` reproduces that sign when the zone holding
the smaller share of the void is the **slower** one, in all 8
<!-- claim: split_scan_sign_ok --> distinct beds of the volume-split scan on
1000 cells and in all 8 <!-- claim: flux_scan_sign_ok --> of the flow-split
scan. Every bed of the volume-split scan was also run on 2000 cells: the sign
rule holds in 8 <!-- claim: split_scan_sign_ok_n2000 --> of 8 there too, the
largest change of a lag where the minority zone is slow is 7.9 %
<!-- claim: split_lag_change_slow_minority -->, where it is fast 48 % <!-- claim: split_lag_change_fast_minority -->: the
fast-minority lags are grid values good to about 50 %, not converged numbers.

The *size* of the displacement bounds nothing. On matched grids of 3200 cells the
largest lag is 0.49 % <!-- claim: flux_scan_lag_max_n3200 --> of the transit time
at the published slice and 1.34 % <!-- claim: split_scan_lag_max_n3200 --> once
the volume fraction moves, a ratio of 2.73
<!-- claim: lag_ratio_matched_n3200 --> (2.93
<!-- claim: lag_ratio_matched_n1600 --> on 1600 cells); the shape residual
reaches 18.0 % <!-- claim: split_scan_resid_max --> on 1000 cells and moves by
1.2 % <!-- claim: split_resid_change_fine --> on 2000.

The time step is Crank-Nicolson on the upwind and diffusive operator with the van
Leer limiter correction implicit. Its step equations are solved by deferred
correction on one factorisation; `newton` with `ConvFluxTVD` solves the same
equations, agrees to the bound above and is slower here, so it is the check, not
the production route.

## Honest limits

- No particle-size test. Figure 9's five marker shapes were deliberately not
  extracted, because Eq. (18) contains no particle-size term.
- One gas pair only, so the ν in Re is untested. The authors say so themselves.
- The extraction was audited and corrected: four merged clusters re-fitted glyph
  by glyph, one row moved off a bare stretch of the dashed curve, one missed
  marker added, three glyphs still fused and listed in the sidecar with the
  measured effect of including the most nearly resolvable of them. The figure
  review is **pending and non-blocking**; the overlay is
  `queue_cases/A2.5/review/fig9-overlay.png` (git-ignored).
- The solver diagnostics are complementary but **not jointly sufficient**:
  reading the voidage off the adjacent row of Table 1 moves `u` by 6.1 %
  <!-- claim: parameter_chain_u_shift --> and `D_L` by 6.9 %
  <!-- claim: parameter_chain_D_L_shift_unseen --> and moves no diagnostic on
  the page beyond its fifth decimal.
- One printed check does not come out. The liquid crossover Reynolds number
  needs ν ≈ 0.0180 <!-- claim: liquid_crossover_nu_implied --> cm²/s where water at
  20 °C is 0.0100. Stated, not repaired.
- Figure 15's detector separation is **not printed**. The reconstruction assumes
  Figure 13's; the page reports what the paper's other test section gives.
- The outlet boundary condition is a **measured** blind spot: with the bed as
  published, an absorbing outlet and an outflow face return identical numbers.

## Reproduce

```bash
cd pages/A2.5-edwards-richardson-dispersion
jupytext --sync index.py      # index.py is the source; index.ipynb is paired with it
jupyter nbconvert --to notebook --execute --inplace index.ipynb
GALLERY_OFFLINE=1 python -m pytest -q test_model.py
```

No PDF needed: the three datasets are committed CSVs. `review/extract-fig9.py`
regenerates the Figure 9 dataset from a 600 dpi render of the page, which the
repository deliberately does not carry.

## Corrections

- 2026-10-31: stated that the two-zone sign rule holds "in 10 of the 12 volume
  splits scanned, and in all 8 of the 8 where the peak moves by more than
  0.05 %", with two exceptions at θ₁ = 0.35 and 0.65, and a largest lag of 1.14 %
  with an 18.7 % shape residual. The 12 rows were 8 distinct beds: (θ₁, u₁/u₂)
  and (1 − θ₁, u₂/u₁) are the same bed with the zones renamed, and the two
  exceptions were one bed counted twice. That exception was the 500-cell grid:
  on 1000 cells the rule holds in all 8 beds, 7 <!-- claim: split_scan_above_floor -->
  of them above 0.05 %, and the same on 2000 cells (7
  <!-- claim: split_scan_above_floor_n2000 --> above 0.05 %); the largest lag
  and residual are 1.33 % <!-- claim: split_scan_lag_max --> and 18.0 % on
  1000 cells.
- 2026-10-31: stated that the peak displacement "more than quadruples" when the
  volume fraction moves; the numbers then published, 0.36 % and 1.14 %, give a
  factor 3.2, and on matched refined grids it is 2.93 on 1600 cells and 2.73
  <!-- claim: lag_ratio_matched_n3200 --> on 3200.
- 2026-10-31: stated that the recovered velocity is exact because
  `Δt̄ = Δx/u` "holds in the discrete system too for a uniform velocity". It
  holds to round-off for a linear scheme (bare upwind, with implicit Euler or
  Crank-Nicolson); the nonlinear van Leer limiter breaks it, and with the limiter
  it is a discretisation error, 3.3e-6 <!-- claim: v5_u_ratio_dev --> at n = 800
  and 4.9e-7 <!-- claim: fig13_u_rel_err --> at n = 1600, order 2.7
  <!-- claim: u_ratio_order -->.
- 2026-10-31: stated that seven merged-cluster centres were re-fitted. The
  sidecar replaces ten rows by eleven re-fitted centres; seven of the ten
  replaced rows moved by 5 px or more.
