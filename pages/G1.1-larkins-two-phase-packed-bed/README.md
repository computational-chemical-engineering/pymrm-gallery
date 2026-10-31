# G1.1 — Larkins two-phase concurrent flow in packed beds

**Source.** Larkins, R. P., White, R. R. and Jeffrey, D. W. (1961),
*Two-Phase Concurrent Flow in Packed Beds*, AIChE Journal **7**(2) 231–239,
[doi:10.1002/aic.690070213](https://doi.org/10.1002/aic.690070213). Read on
300 ppi renders, the scan's native resolution. Book page = PDF page + 230.
**PDF page 1 opens with the previous article's nomenclature**; Larkins et al.
start below the rule on the same page.

- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`Larkins` for the printed constants and Eqs. (17)–(23),
  `recover_mu_air`, `Column`, the pymrm column `ColumnFV`, the quadrature route
  `quad_outlet`, the paper's sectional recipe `recipe_outlet`, the hand-written
  `backward_euler_march` and the root `threshold_length`), tested by
  `test_model.py`.
- **Runtime:** about 35 seconds.
- `data/`: Table 3 (measured), Table 4 (the authors' arithmetic) and every
  printed constant with its book page, each with a provenance sidecar.

## Why this paper got a page

It prints a **worked example**. Table 3 (book p. 234) is eight runs of measured
two-phase data on ⅜-in. Raschig rings; Table 4 (book p. 235) is the authors'
own arithmetic on the identically numbered rows. Every constant that turns one
into the other is printed somewhere in the nine pages, so the chain can be
redone cell by cell.

## What the page establishes

| | |
|---|---|
| liquid Reynolds column reproduced | 0.174 <!-- claim: re_liquid_max_rel_error_pct --> % worst, 0.147 <!-- claim: re_liquid_rms_rel_error_pct --> % r.m.s. |
| δ_l (Ergun, α = 266, β = 2.33) | 0.218 <!-- claim: delta_l_max_rel_error_pct --> % worst, 0.176 <!-- claim: delta_l_rms_rel_error_pct --> % r.m.s. |
| δ_g | 0.416 <!-- claim: delta_g_max_rel_error_pct --> % worst |
| δ_lg (a measurement, corrected by two densities) | 0.639 <!-- claim: delta_lg_max_rel_error_pct --> % worst, 0.260 <!-- claim: delta_lg_rms_rel_error_pct --> % r.m.s. |
| **run 5 does not close** | **-5.99 <!-- claim: run5_delta_lg_rel_error_pct --> %**, an inconsistency between the paper's two tables |
| air viscosity recovered, linear in T | 0.0259 <!-- claim: mu_air_linear_fit_max_resid_pct --> % worst residual over seven rows |
| that line extrapolated to Table 2's own 80 °F | 0.819 <!-- claim: mu_air_at_80F_vs_table2_pct --> % from the printed 0.0192 cP (0.672 <!-- claim: mu_air_at_80F_offset_corrected_vs_table2_pct --> % once the offset measured on the liquid rows is divided out) |
| average-pressure claim, one phase, horizontal | exact: proved, and at round-off on the page (below 1e-11 <!-- claim: gas_horizontal_arith_mean_max_rel_error -->) |
| the same with the geometric mean of the same terminals | 1.12e-05 <!-- claim: gas_horizontal_geom_mean_max_rel_error --> |
| one-section recipe over a 4 ft bed, two phases | 0.332 % <!-- claim: two_phase_recipe_rel_error_4ft --> on the pressure drop |
| bed length at which it reaches 1 %, **root-found**, on Table 3's largest-friction-loss row | 5.6012 <!-- claim: two_phase_recipe_1pct_length_ft --> ft, inside the paper's own 7 ft column |
| the same root-find over the other two-phase rows | up to 466.9537 <!-- claim: two_phase_recipe_1pct_length_max_ft --> ft; 1 <!-- claim: two_phase_recipe_1pct_runs_inside_experimental_column --> of the 6 inside the 7 ft column, 4 <!-- claim: two_phase_recipe_1pct_runs_inside_commercial_reactor --> inside the 40 ft reactor |
| sectional march, observed order | 1.9999 <!-- claim: sectioning_observed_order --> (the midpoint rule, second order) |
| pymrm finite volumes with the `{"outflow": True}` outlet, error at 1600 cells | 3.78e-05 <!-- claim: fv_rel_error_at_n_fine --> (first order, observed 1.0003 <!-- claim: fv_observed_order -->) |
| its Richardson pair (800, 1600) against the adaptive quadrature | 1.47e-08 <!-- claim: fv_extrapolated_vs_quadrature_rel --> relative, 5.35 <!-- claim: fv_extrap_orders_below_recipe --> orders below the recipe error it measures; the remainder is O(h²), observed order 2.00 <!-- claim: fv_richardson_remainder_order --> (the smaller zero-gradient remainders come from a cancellation of errors of opposite sign) |
| the outflow solve against a hand-written backward-Euler march, and its own discrete balance | both below 1e-10 <!-- claim: fv_outflow_vs_hand_march_rel --> |
| a zero-gradient outlet instead | folded into the flux as (9 P_N − P_{N−1})/8, last step shortened to 8/9 of h; its face value lands on the outflow solve to O(h²); its last cell is 3.84 <!-- claim: cbv_over_upwind_face_outlet_error_ratio -->× closer by a cancellation, and a balance written on it fails, by 1.2e-3 <!-- claim: zg_lastcell_balance_gap_n100 --> of the pressure drop at 100 cells (pitfall P8) |

## The run-5 discrepancy, stated and not repaired

Run 5 has no gas, so R_l = 100 %, so Eq. (21)'s mixture density equals the
liquid density and the manometer correction of book p. 236 cancels exactly:
Eq. (20) reduces to δ_lg = ΔP/ΔL. **Table 3 prints `2.321`; Table 4 prints
`2.469`.** The run's liquid rate is corroborated twice inside Table 4, and the
same procedure closes the other seven rows, run 160 included, at
+0.134 <!-- claim: run160_delta_lg_rel_error_pct --> %.

**Run 160 is the mirror single-phase row, not the same reduction with the
phases swapped.** With R_l = 0 nothing cancels: the manometer correction still
applies in full and Eq. (21) still contributes ρ_g, and dropping the correction
alone would move that row by +22.29 <!-- claim: run160_manometer_correction_worth_pct --> %.

**Both readings are weighed on the paper's own scale.** Taking Table 3's cell
makes the two-phase parameter 0.9408 <!-- claim: run5_two_phase_parameter_from_table3 -->,
-5.92 <!-- claim: run5_parameter_gap_from_unity_pct --> % off the unity the
paper says the parameter approaches as χ → ∞; the paper's own scatter for its
correlation is a 13 % standard deviation, and run 160 at χ → 0 is itself
+0.90 <!-- claim: run160_parameter_off_unity_pct --> % off unity. Against it,
run 5's 2.469 sits 0.081 <!-- claim: run5_vs_ergun_pct --> % from Ergun's δ_l,
11.5 <!-- claim: run5_closer_factor --> times closer than run 160's
0.93 <!-- claim: run160_vs_ergun_pct --> %. That second arm is the same
statistic as the first (with no gas, |δ_lg/δ_l − 1| *is* the parameter minus
one), and against a 13 % scatter alone a 0.081 % residual is the more probable
of the two; it is evidence only under the hypothesis that Table 4's cell was
back-computed from Eq. (17). **The page does not choose.**

## A second printed discrepancy: Eq. (23) is printed twice

Book p. 237 prints the liquid-saturation fit with constant term −0.774 as the
typeset equation and −0.744 inside Figure 7 on the facing column. The difference
is worth 7.15 <!-- claim: eq23_fig7_constant_worth_pct --> % in R_l. Against the
six Table 4 rows inside the declared range the Figure 7 reading is better in
r.m.s. (18.33 <!-- claim: eq23_fig7_rms_dev_pct --> % against
19.55 <!-- claim: eq23_rms_dev_pct --> %) and worse on the paper's own "within
20 %" count, so the page adopts neither. Both readings put the ceiling below a
full column (0.7730 <!-- claim: eq23_fig7_ceiling_liquid_saturation --> and
0.7214 <!-- claim: eq23_ceiling_liquid_saturation -->).

## Fit versus test

Eqs. (22) and (23) were fitted to the body of data these eight rows are a
sample of, so comparing them here is a **goodness of fit, not a test**. Eq. (22)
beats the best single constant by 4.05 <!-- claim: eq22_gain_over_bestconst -->×
in log units over eight rows; two of those rows are the χ limits where it
returns 1 for any constants, and on the six that discriminate the gain is
2.48 <!-- claim: eq22_free_rows_gain_over_bestconst -->× and the r.m.s.
14.79 <!-- claim: eq22_free_rows_rms_dev_pct --> %, the figure to set beside the
paper's own 13 %. Eq. (23)'s best-constant null is
117.9 <!-- claim: eq23_null_bestconst_rms_dev_pct --> % r.m.s.

## Scope

Nothing is digitised; the 600+ points behind the correlations are not
tabulated. Reference 4, Larkins's 1959 thesis, is not on disk. The column is a
design calculation on run 63's rates at the printed 90 psig equipment limit,
validated against three independent routes rather than a measurement; run 63
is the extreme row, so the 1 % length is a worst case.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/G1.1-larkins-two-phase-packed-bed/index.py   # after editing index.py
python scripts/run_kernels.py pages/G1.1-larkins-two-phase-packed-bed
GALLERY_OFFLINE=1 PYTHONPATH=shared python -m pytest -q pages/G1.1-larkins-two-phase-packed-bed/test_model.py
```

## Corrections

- 2026-10-30: stated that neither extremum of the two correlations was sampled
  off a grid. Eq. (22)'s maximum was the scan point itself: the scan landed
  exactly on χ = 1, where the difference is zero, and the code fell back to the
  sampled maximum. The scan now brackets the root and `brentq` finds it; the
  value is unchanged (χ = 1, Eq. (22) maximum 4.2133 <!-- claim: eq22_max_parameter -->).
- 2026-10-30: stated that z(P) by adaptive quadrature is inverted by bisection.
  It is inverted by `brentq`.
- 2026-10-30: the sectional recipe did not check that its fixed point
  converged, and at the upper end of run 86's root bracket (0.97 of the choking
  length) it had not. That evaluation only fixed the sign at the bracket end,
  so no published number changed (466.9537 <!-- claim: two_phase_recipe_1pct_length_max_ft --> ft);
  the recipe now raises when it does not converge, and the bracket is widened
  from short beds through converged evaluations only.
- 2026-10-30: nine numbers were labelled structural, that is not movable by
  any defect, four of them because "nothing transcribed on this page enters
  it". Eight of those labels were false: each of the eight moves past its check
  under an injected defect (a misread table cell or constant, a mis-paired
  grid). The ninth, Eq. (22)'s argmax, was held at exactly 1 by the scan
  fallback of the first correction, so no row could move it then; with the
  root found by `brentq` a looser root tolerance moves it. Every claim now has
  a break row that fails its check.
