# A3.8: The Onda correlations, and which of their fourteen constants the paper can see

Onda, K., Takeuchi, H. & Okumoto, Y., *"Mass transfer coefficients between gas
and liquid phases in packed columns"*, **Journal of Chemical Engineering of
Japan** **1**(1) 56-62 (1968),
[doi:10.1252/jcej.1.56](https://doi.org/10.1252/jcej.1.56).

Three coupled correlations (wetted area, liquid-side coefficient, gas-side
coefficient) with the area feeding the liquid side and not the gas side.

- **Structures:** `S3` (1-D steady BVP)
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`OndaSet`, the counter-current `Absorber` and its closed form),
  tested by `test_model.py`. The paper's algebra, the studies and the break table
  stay in the notebook.
- **Data:** three CSVs in `data/`, each with a `.meta.yaml` sidecar: the fourteen
  printed constants of eqs. (1) to (3) and the four of eqs. (4) and (5), Table 1
  (page 57), and the scalar results printed in the text. Every value is a printed
  numeral read from 300 ppi renders, each equation and table row cropped to one
  line; eqs. (1) to (3) were read twice, from the article and from their reprint
  in the companion article that starts on page 62. Nothing is digitised.
- **Runtime:** about 8 s.

The article was identified by reading it: journal page 56 carries the title, the
byline and *"Received on July 10, 1967"*, the running feet read *"VOL.1 NO.1
1968"*, and CrossRef returns the same title, authors, volume, issue and pages for
`10.1252/jcej.1.56`.

## What the page shows

1. **Eq. (4) is eq. (3) rearranged, and recovering it pins all four gas-side
   constants.** A `sympy` reduction and a numerical evaluate-and-regress over
   random property combinations both give a prefactor of
   0.774443 <!-- claim: jd_prefactor_symbolic --> against the printed `0.771`,
   +0.447 <!-- claim: jd_prefactor_dev_pct --> %, and the exponent `-0.30` with a
   gap below 1e-12 <!-- claim: jd_exponent_dev -->. The Schmidt number cancels
   exactly, which is what tests eq. (3)'s `1/3`. The gap is inside the printed
   precision of `a_t D_p = 3.4`: 3.4076 <!-- claim: jd_atDp_for_printed_0771 -->
   reproduces `0.771` exactly.
2. **A ceiling the correlation cannot exceed, which Onda's own table breaks.**
   Eqs. (1) and (2) force `dln(k_L a)/dlnL` into (2/3, 0.8 <!-- claim: window_hi -->]
   for any packing, liquid and flow rate; the bound holds no prefactor, so no refit
   moves it. Of the six exponents of Table 1, 2 <!-- claim: table1_n_above_ceiling -->
   (both spheres) sit above the ceiling, the worst by
   0.06 <!-- claim: table1_worst_excess -->, each at least six half-units of the
   printed `n` from an edge.
3. **Which constants a design can see.** Placed in a pymrm counter-current
   absorber, the packed height's elasticity to every constant is computed three
   ways. Eq. (1)'s `0.75` is exactly inert wherever the packing and the liquid share
   a surface tension (change below 1e-14 <!-- claim: inert_075_reldev --> for a
   tenfold exponent), and moves the height by 119.6 % <!-- claim: live_075_reldev -->
   for the same change at `sigma_c/sigma = 0.85`. Eq. (3)'s `1/3` is worth
   0.73 % <!-- claim: e_ScG_doubled_reldev --> of packed height when doubled for an
   air system, 0.024 <!-- claim: e_ScG_doubled_over_height_band --> of the paper's
   own `±30 %` band on height: no design calculation can see it.

## Agreement

| Check | Result |
|---|---|
| eq. (4) prefactor from eq. (3), two routes | 0.774443 <!-- claim: jd_prefactor_numeric --> numerically; routes agree below 1e-12 <!-- claim: jd_two_routes_spread --> |
| eq. (4) against Shulman's eq. (5) over `100 < X < 10,000` | -14.95 <!-- claim: eq45_dev_at_X100_pct --> % to +12.12 <!-- claim: eq45_dev_at_X10000_pct --> %, crossing at X = 1486 <!-- claim: eq45_crossover_X -->, mean over ln X 6.96 <!-- claim: eq45_mean_dev_pct --> % |
| window, analytic against numerical differentiation | below 1e-8 <!-- claim: n_window_two_routes --> |
| Table 1 inside the window | 4 <!-- claim: table1_n_inside --> of six rows |
| mean error in `n`: ceiling, no-coupling null, fitted constant | 0.0667 <!-- claim: mae_window_ceiling -->, 0.0967 <!-- claim: mae_null_no_area_coupling -->, 0.0533 <!-- claim: mae_null_fitted_constant -->; paired t 0.80 <!-- claim: t_ceiling_vs_no_coupling --> and 0.60 <!-- claim: t_fitted_vs_ceiling -->: six rows separate none of these point predictors, while every exponent exceeds 2/3 (t = 3.50 <!-- claim: t_n_above_no_coupling -->) |
| elasticities: chain rule against differences on `Z`, and through the column | below 1e-6 <!-- claim: elasticity_two_routes --> and below 1e-5 <!-- claim: elasticity_column_route --> |
| pymrm column against its closed form, 800 cells | 1.90e-6 <!-- claim: col_raw_reldev -->; orders 1.930 <!-- claim: col_order_face_min --> to 1.990 <!-- claim: col_order_face_max --> on 25 to 800 cells; Richardson limit below 1e-7 <!-- claim: col_richardson_reldev --> |
| the same outlet read at the last cell centre | 7.29e-4 <!-- claim: col_cell_centre_reldev -->, order 1.020 <!-- claim: col_order_cell_finest --> on the finest triple, falling toward 1 from above, solute balance off by 6.1e-4 <!-- claim: col_cell_imbalance --> (face values: below 1e-13 <!-- claim: col_imbalance -->) |

Deviations are `(computed - printed)/printed`.

## What pymrm adds

Eq. (6) as a boundary value problem instead of a quadrature: two counter-current
convective fluxes with velocity `[[G_M, -L_M]]`, a per-cell linear coupling and
outlet values read at the faces. With upwind faces the problem is linear and takes
one sparse solve; van Leer faces make it nonlinear, and `newton` on the exact
Jacobian of `ConvFluxTVD` solves it. On this page that buys a third route to the
elasticities of packed height, through the solved column. It adds no accuracy:
the closed form of the same problem is exact. The correlations, the eq. (4)
identity and the exponent window are `sympy` and `numpy` and owe pymrm nothing.

## What this page does NOT do

- **No out-of-sample test of eqs. (2) and (3).** They were fitted to the data any
  comparison here would use, Table 1's runs included (page 57 says so).
- **No test of `a = a_w`.** The area never appears except multiplied by a
  coefficient, so nothing in the article separates them.
- **No figure digitisation.** Figs. 1, 2, 3 and 7 are scatter plots of the fits,
  Fig. 8 is `Z_cal` against `Z_act`, Fig. 9 the resistance split for distillation.
- **No comparison with A3.9 (Billet-Schultes) or A3.10 (Rocha-Bravo-Fair).**
  Neither source is on disk.
- **The illustrative absorber is not Onda's.** Air and water at 25 °C are
  standard-table properties; `sigma_c/sigma = 0.85`, `m = 1` and the `95 %` duty are
  page choices. The geometry, `a_t D_p = 6(1-eps) = 3.4` for spheres (page 58), is
  the only one the article supplies.

## Corrections

- 2026-11-09: said "The *outflow face* stays first order". The outlet face value
  from `compute_boundary_values` converges at orders
  1.930 <!-- claim: col_order_face_min --> to 1.990 <!-- claim: col_order_face_max -->
  on 25 to 800 cells; the first-order reading is the last cell centre.
- 2026-11-09: said "finds $L_M/(mG_M) = 1$ as the boundary where a counter-current column
  stops existing; eq. (6)'s integrand diverges there, which is harder to read
  than a solve that stops converging". The column solve converges for every `m` (at
  `m = 1.1 L_M/G_M` and the design height, `y_out/y_in` =
  0.4449 <!-- claim: y_out_beyond_pinch_over_y_in -->). What stops existing is a height
  that meets the duty, beyond `m` = 6.769 <!-- claim: m_duty_limit -->, where
  `L_M/(mG_M)` equals the removal; at `m = 1.1 L_M/G_M` the root search finds
  0 <!-- claim: duty_heights_beyond_pinch --> heights.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: the three datasets are committed CSVs in `data/`.
