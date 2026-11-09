# A3.4: Wakao and Funazkri, particle-to-fluid mass transfer in packed beds

`Sh = 2 + 1.1 Sc^(1/3) Re^0.6`, from Wakao, N. & Funazkri, T., *Chemical
Engineering Science* **33**(10) 1375-1384 (1978),
[doi:10.1016/0009-2509(78)85120-3](https://doi.org/10.1016/0009-2509(78)85120-3).

The paper measures nothing. It takes thirty-five years of published packed-bed
mass-transfer data, observes that the coefficients were extracted with the wrong
bed model, and extracts them again with the right one. The page rebuilds that
re-analysis, and tests the two fitted constants against the digitised figure they
were fitted on.

> **Status of the data.** Both figure extractions are pending maintainer review;
> each sidecar's `review` block records what was asked. Figure 2's marker
> identifications were confirmed on 2026-08-02 (*"mostly symbols are identified
> well, but often the marker is off centre"*); the shape-fit re-centring that
> followed, the markers it added, and the whole of Figure 3 have not been reviewed.
> Results 1 to 3 use no figure data.

- **Structures:** `S3` (steady 1-D convection-diffusion-reaction, Danckwerts)
- **Source of the page:** `index.py` (paired with `index.ipynb`); the bed model,
  its closed form, the dispersion laws, the inversion and the three correlations
  in `model.py`, tested by `test_model.py`. The figure analysis and the break
  table stay in the notebook.
- **Runtime:** about 75 s, most of it the bin bootstrap rerun under the break
  table's defects.

## Contents

| file | what it is |
| --- | --- |
| `data/wakao-funazkri-1978-fig2.csv` | Figure 2, gas phase: one row per marker position ever recorded (crosshair pass and shape-fit pass), so the re-centring is a paired sample; group by `marker_id` |
| `data/wakao-funazkri-1978-fig3.csv` | Figure 3, liquid phase, `(Sh-2)/Sc^(1/3)` against `Re`: the data alpha and beta were fitted on. One marker set per erasure-band half-width; `band_hw == 5` is the shipped set, the rest the systematic uncertainty |
| `data/wakao-funazkri-1978-fig3-line.csv` | the correlation Figure 3 draws inside itself, traced from its own ink: a calibration control whose answer is known |
| `data/wakao-funazkri-1978-parameters.csv` | every constant read off 600 dpi page renders, plus both figures' axis calibrations, so the page converts pixels itself |

Each CSV has a `.meta.yaml` provenance sidecar with the source page, the figure,
the method and the reading precision. Read the figure sidecars before using
those files: they record what each extraction is known to miss, and for Figure 3
which direction the recall loss biases a refit.

## What the page shows

1. **The whole paper is one ratio.** Eq. (2) and the `Pe = 2` assumption it
   replaces share the convective term, so `D_ax(eq. 2) / D_ax(Pe = 2) = 1 + 40/(Sc Re)`,
   which doubles the dispersion coefficient at `Sc Re` = 40 <!-- claim: dax_double_ScRe -->:
   `Re` = 66.7 <!-- claim: dax_double_Re_gas --> for a gas and
   0.04 <!-- claim: dax_double_Re_liquid --> for a liquid. That predicts both
   halves of the paper's first conclusion.
2. **The liquid correction is small everywhere in the range.** For a liquid
   (`Sc = 1000`) the largest `Sh/Sh'` correction over `3 <= Re <= 10^4` is
   1.33 % <!-- claim: liquid_correction_max --> at `L/d_p = 10` and
   1.39 % <!-- claim: liquid_correction_max_all_Ldp --> at most for `L/d_p` up
   to `100`, both at the low end of the range. For a gas at `Re = 3` and
   `L/d_p = 5` it is a factor 6.92 <!-- claim: gas_correction_Re3_Ldp5 -->.
3. **The paper's own recalculation, redone from printed algebra.** Pushing
   Petrovic and Thodos' eq. (11) through eq. (7) twice (at their `Pe = 2`, then
   with eq. (2)) returns eq. (12) to 1.67 % <!-- claim: eq11_reanalysed_vs_eq12_mean_abs_dev_logRe -->
   on average over `log Re` and 4.48 % <!-- claim: eq11_reanalysed_vs_eq12_worst --> at worst
   over `100 <= Re <= 900`, nothing fitted. Changing the bed height from
   `L/d_p = 10` to `40` moves the mean deviation by
   0.24 % <!-- claim: eq11_reanalysis_Ldp_effect_mean --> of eq. (12).
4. **The constants against the data they were fitted on.** Hold beta at the
   printed value and the 182 <!-- claim: n_fig3_markers --> Figure 3 markers
   return alpha = 1.1000 <!-- claim: refit_liquid_alpha_beta_fixed -->; across
   every erasure band and both axis calibrations it stays between
   1.098 <!-- claim: refit_liquid_alpha_beta_fixed_logmetric_lo --> and
   1.118 <!-- claim: refit_liquid_alpha_beta_fixed_logmetric_hi -->, and between
   1.007 <!-- claim: refit_liquid_alpha_beta_fixed_grid_lo --> and
   1.118 <!-- claim: refit_liquid_alpha_beta_fixed_grid_hi --> once the estimator
   is varied too. The free alpha over the same choices runs from
   1.05 <!-- claim: refit_liquid_alpha_free_grid_lo --> to
   1.46 <!-- claim: refit_liquid_alpha_free_grid_hi -->. That contrast is the
   result; the page does not claim the printed constants are wrong.

## Agreement

| Check | Result |
|---|---|
| pymrm solve of eq. (4) against eq. (7), ladder of 200 to 3200 cells | observed orders 1.00 <!-- claim: pymrm_order_min --> to 1.03 <!-- claim: pymrm_order_max -->, first-order upwind; 0.00214 <!-- claim: pymrm_vs_eq7_n1600 --> at 1600 cells; the extrapolated value below 1e-5 <!-- claim: pymrm_vs_eq7_richardson --> from eq. (7) |
| the re-analysis on the pymrm solve instead of eq. (7) | below 1e-3 <!-- claim: reanalysis_route_gap --> from the eq. (7) route |
| eq. (7) to eq. (8) as `Pe_L` grows | `ln(theta_7/theta_8) Pe_L / N^2` = 0.99995 <!-- claim: eq7_eq8_rate --> at `Pe_L = 1e5`, limit `1` |
| eq. (11) against the stagnant floor `Sh = 2` | crosses it at `Re` = 4.59 <!-- claim: eq11_crosses_sh2_Re --> |
| eq. (12) against the 79 <!-- claim: n_fig2_markers --> Figure 2 gas markers | mean absolute deviation 27.2 % <!-- claim: eq12_vs_fig2_mean_abs_dev -->, bias +9.7 % <!-- claim: eq12_vs_fig2_bias --> with a standard error of 3.6 % <!-- claim: eq12_vs_fig2_bias_se --> |
| eq. (12) against the Figure 3 liquid markers | mean absolute deviation 17.6 % <!-- claim: eq12_vs_fig3_mean_abs_dev -->, bias +2.3 % <!-- claim: eq12_vs_fig3_bias --> with a standard error of 1.6 % <!-- claim: eq12_vs_fig3_bias_se --> |
| the trend the paper claims for Figure 2 | eq. (11)'s binned bias changes 3.5 <!-- claim: fig2_span_ratio_eq11_over_eq12 --> times as much across the range as eq. (12)'s |
| Figure 3's drawn line, refitted from its own ink | 1.107 <!-- claim: fig3_drawn_line_alpha_skew --> / 0.598 <!-- claim: fig3_drawn_line_beta_skew --> with the skew-corrected calibration, 1.126 <!-- claim: fig3_drawn_line_alpha_flat --> / 0.595 <!-- claim: fig3_drawn_line_beta_flat --> without |
| free log-log fit of Figure 3 | alpha 1.344 <!-- claim: refit_liquid_alpha -->, beta 0.5496 <!-- claim: refit_liquid_beta -->; bin-bootstrap `95 %` interval of alpha 1.04 <!-- claim: refit_liquid_alpha_cluster_ci_lo --> to 1.53 <!-- claim: refit_liquid_alpha_cluster_ci_hi -->, which contains the printed value |
| its error model | quadratic term at ordinary t = 5.0 <!-- claim: refit_liquid_quadratic_term_z --> on 179 <!-- claim: refit_liquid_quadratic_dof --> degrees of freedom; bin-bootstrap standard error 2.5 <!-- claim: refit_liquid_cluster_se_inflation_alpha --> times the ordinary one |

Deviations are `(model - measured)/measured`.

## What pymrm adds

The bed model as a solved boundary-value problem: pymrm discretises eq. (4) with
Danckwerts conditions on the outward normal, reads the outlet through the same
condition, and on a grid ladder shows that the paper's eq. (7) is its solution.
The page's inversions run on eq. (7), which is exact; result 3 runs the
re-analysis on the pymrm solve as well, with the same root finder, and the two
agree to the discretisation error. That route is what a bed model with no closed
form would need. The correlation and the figure statistics owe pymrm nothing.

## Reusing this

**The reusable object is the bed model, not the correlation.** `theta_pymrm`
solves eq. (4) in `(N, Pe_L)` only, so it applies unchanged to a first-order
catalytic bed, an adsorber, or any closed vessel with a linear sink.
`theta_eq7` is its analytical solution and a free regression test.

**`Sh_of_theta` converts between dispersion assumptions.** Every packed-bed
transfer coefficient carries one; pre-1978 correlations and modern ones are
usually not the same quantity.

**Use eq. (12) with eq. (2), never on its own.** Conclusion 3 of the paper: the
correlation was extracted under a particular dispersion coefficient and is only
self-consistent when paired with it.

## Corrections

- 2026-11-09: said "the liquid-phase correction never exceeds 1.33 % over `3 ≤ Re ≤ 10⁴`".
  That maximum was computed at one bed length, `L/d_p = 10`, where it is
  1.33 % <!-- claim: liquid_correction_max -->; the bed length is swept, not
  known, and the maximum grows with it, to 1.39 % <!-- claim: liquid_correction_max_all_Ldp -->
  at `L/d_p = 100`.
- 2026-11-09: said "Below Re ~ 50 for a gas the inversion is ill conditioned".
  Where a gas bed saturates depends on its length. With an outlet resolution of
  `1e-3`, a bed obeying eq. (12) saturates below `Re` = 85.0 <!-- claim: gas_saturation_Re_Ldp10 -->
  at `L/d_p = 10`, and nowhere in `3 <= Re <= 10^4` at `L/d_p = 5`.
- 2026-11-09: described the break table's first column as "The metric of check 1: pymrm's outlet against eq. (7), after one Richardson step."
  It was the error of the 1600-cell solve with no extrapolation,
  0.00214 <!-- claim: pymrm_vs_eq7_n1600 -->; the extrapolated value is below
  1e-5 <!-- claim: pymrm_vs_eq7_richardson --> from eq. (7).
- 2026-11-09: said that with eq. (4) discretised "the closed form becomes a check on the solver rather than the method".
  Every inversion, correction factor and re-analysis on the page ran on eq. (7),
  the closed form. pymrm validates eq. (7), and result 3 now also runs on the
  pymrm solve, below 1e-3 <!-- claim: reanalysis_route_gap --> from the eq. (7)
  route.
- 2026-11-09: said "The free α over the same choices runs 1.05 to 1.47".
  The largest value is 1.4648 <!-- claim: refit_liquid_alpha_free_grid_hi -->, so the range is 1.05 <!-- claim: refit_liquid_alpha_free_grid_lo -->
  to 1.46 <!-- claim: refit_liquid_alpha_free_grid_hi -->.
- 2026-11-09: the Figure 2 sidecar said "12 of the 81 scored below that". In the file 16 <!-- claim: n_crosshair_kept --> of the 81 <!-- claim: n_paired --> paired markers score below the cut and keep the crosshair; 65 <!-- claim: n_recentred --> were re-centred.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: every dataset is a committed CSV in `data/`.
