# A3.12: Yagi–Kunii effective thermal conductivity

A packed bed conducts heat by several routes at once. Yagi & Kunii's Equation (15)
adds them (solid conduction, the thin gas film at the contact points, radiation
surface-to-surface, radiation void-to-void), and Equation (6) puts Ranz's lateral
mixing in parallel once gas flows. It is still the skeleton of every stagnant
conductivity correlation in use.

- **Structures:** `S6` (cylindrical `nu=1` radial conduction, for the apparatus)
- **Reference:** Yagi & Kunii (1957), *AIChE J* **3**(3) 373–381,
  doi:10.1002/aic.690030317
- **Runtime:** ~25 s

## Three kinds of result, kept apart

The paper contains theory, *other people's* data re-analysed, and the authors' own
measurements. The page labels each:

| | where | label |
|---|---|---|
| Equations (1)–(19) | pages 374–376 | reproduced, cross-checked against their own printed limits |
| Table 1 (twenty flow studies) and Table 3 (twenty stagnant rows) | pages 373, 377 | third-party measurement + the authors' reduction of it |
| iron spheres, porcelain, cement clinker, firebrick, Raschig rings, 0–1000 °C | **Figures 13–17, the plotted points** | **scoped out**: figures, and no maintainer is available for a digitisation review |
| the smooth **theoretical** lines in those same figures | Figures 13–17 | *read*: the authors' computed output, 17 curve-end values, used as a reproduction check |
| the endpoints of Figure 12's five temperature profiles | page 378 | *read*: 10 values off a labelled ordinate, to fix the Δt the apparatus ran |

## The finding

**One of Equation (15)'s five ingredients is untestable on the only table of numbers
in the paper, and it is the radiation term the paper's central claim rests on.**
Equation (15) needs `phi`, the effective contact-film thickness, which the paper
supplies only as Figures 9 and 10. With `phi` free per row, the solid-to-solid
coefficient `h_rs` enters only through the sum `1/phi + Dp·h_rs/kg` and is
algebraically degenerate with it. Deleting one term at a time and re-solving for
`phi` on the Table 3 rows that carry a particle diameter:

| term deleted | rows `phi` can no longer match | worst shortfall | largest `phi` it demands |
|---|---|---|---|
| solid-solid radiation, Eq. (7) | 0 <!-- claim: delete_hrs_rows_unmatched --> | none | 0.054 <!-- claim: delete_hrs_max_phi_needed -->, inside Table 4's printed band |
| void-void radiation, Eq. (8) | 2 <!-- claim: delete_hrv_rows_unmatched --> | 14.2 <!-- claim: delete_hrv_worst_shortfall_pct --> % | |
| solid conduction, `gamma(kg/ks)` | 0 <!-- claim: delete_gamma_rows_unmatched --> | none | 0.42 <!-- claim: delete_gamma_max_phi_needed -->: 12 <!-- claim: delete_gamma_rows_phi_above_printed --> rows above every printed value, 2 <!-- claim: delete_gamma_rows_off_fig9_axis --> off Figure 9's axis |
| *control:* numerator halved | 6 <!-- claim: delete_control_beta_half_rows_unmatched --> | 46.7 <!-- claim: delete_control_beta_half_worst_shortfall_pct --> % | |

Deleting Eq. (7) entirely and re-solving for `phi` reproduces every printed
calculated value to round-off (below 1e-13 <!-- claim: delete_hrs_exact_reproduction_resid -->),
with a `phi` that stays inside Table 4's printed band. The control shows the test is
not blind. **`gamma`, by contrast, is identifiable**: it looks degenerate only if the
search may run to `phi = 1`, a "film" as thick as the particle. And the paper's
central conclusion, radiation dominant at temperatures "higher than 400 °C",
concerns a regime no row of Table 3 reaches.

What breaks the degeneracy is temperature: `phi` is a geometric ratio, `h_rs` scales
as `T³`. That experiment exists, it is Figures 13 to 17, and its measured points are
the part of the paper this page cannot read.

## Agreement

| check | result |
|---|---|
| Eq. (16) = Eq. (15) with `h_r = 0`, `gamma = 1` | below 1e-13 <!-- claim: V1_identity_eq16_from_eq15 --> over random draws |
| Eq. (14) = the `kg → 0` limit of Eq. (13) | residual vanishes as `O(kg)`, order 1.0000 <!-- claim: V2_eq14_limit_order --> |
| Eq. (13) rebuilt from the Eqs. (9) to (12) network | below 1e-13 <!-- claim: V3_eq13_from_network --> |
| Eq. (19)'s constant `0.743`, from the closed form with `2*pi*l` and the paper's `2.3` for `ln 10` | +0.048 <!-- claim: eq19_constant_2pi_vs_printed_pct --> %: the printed value to three digits |
| the same from a pymrm `nu=1` annulus solve, which carries `ln` exactly | +0.16 <!-- claim: V4_pymrm_constant_vs_printed_pct --> %: the paper's own rounding of `ln 10` to `2.3` |
| Table 4's `ks/kg = 4.48` from its own `ks` and the printed `kg(400 °C)` | -0.13 <!-- claim: table4_ks_over_kg_roundtrip_pct --> % |
| **Eq. (15) vs the ends of the paper's own printed theoretical curves, Figs. 13 to 17** | **mean 2.4 <!-- claim: fig13to17_curve_ends_mean_abs_pct --> %, worst 4.2 <!-- claim: fig13to17_curve_ends_worst_abs_pct --> %**, all 14 <!-- claim: fig13to17_curve_ends_positive_count --> deviations positive |
| Table 3, the measured rows vs the paper's calculated column | mean \|dev\| 15.3 <!-- claim: paper_calcd_vs_exp_mean_abs_pct --> %, mean dev -3.5 <!-- claim: paper_calcd_vs_exp_bias_pct --> % with a standard error of 5.5 <!-- claim: paper_calcd_vs_exp_bias_se_pct --> %, worst 52.0 <!-- claim: paper_calcd_vs_exp_worst_abs_pct --> % |
| Table 1's intercept vs Table 3's `Exp.`, five <!-- claim: table1_vs_table3_daggered_rows --> shared references | 8.8 <!-- claim: table1_vs_table3_exp_max_abs_pct --> % worst, but identical on 3 <!-- claim: table1_vs_table3_identical_rows -->, see below |

Three stated conclusions reproduced from printed parameters: radiation reaches half
of `ke0` at 461 °C <!-- claim: radiation_half_temp_iron_sphere_C --> (iron
spheres) and 265 °C <!-- claim: radiation_half_temp_firebrick_sphere_C -->
(insulating firebrick; 233 °C <!-- claim: radiation_half_temp_firebrick_sphere_C_ks_over_kg_const -->
if Table 4's printed `ks/kg` is held instead of its `ks`) against the paper's "higher
than 400 °C"; `0.18 mm` cement clinker varies 9.5 <!-- claim: fine_clinker_spread_pct_200_1000C --> %
over `200 to 1000 °C` against 154 <!-- claim: coarse_clinker_spread_pct_200_1000C --> %
for the 5 mm grade ("nearly constant for particles smaller than 0.35 mm"); randomly
packed Raschig rings run 12.9 <!-- claim: raschig_random_over_regular_min_pct --> %
to 24.4 <!-- claim: raschig_random_over_regular_max_pct --> % above regularly
packed with no fitted parameter ("about 20 % greater").

**The one comparison that is not what it looks like.** Table 3 daggers all five of the
shared references as "Extrapolated values, `N_ReM → 0`", so its `Exp.` column and
Table 1's intercept are *the same quantity*: the authors' own zero-flow extrapolation
of the same third-party flow experiments. Three of the five agree exactly because
they are one number printed twice. The page keeps the comparison as a transcription
and internal-consistency check and labels it as one.

## Two printed defects, both proved from the paper's own results

- **Equation (17) is a factor of two low.** It sets `pi*l` where the steady
  cylindrical solution requires `2*pi*l`. The paper's own Eq. (19) constant, `0.743`,
  proves it: the printed form gives 2.00095 <!-- claim: eq17_as_printed_over_eq19_ratio -->
  times that. Anyone re-deriving the reduction from Eq. (17) would report every `ke0`
  twice too large. Eq. (19) itself is right, so no measurement in the paper is
  affected.
- **`gamma = 1.2 mm/9.0 mm` is printed as `0.0134`**, ten times too small for that
  ratio. Table 4 carries the ratio, so the paper's own calculations are unaffected;
  using the decimal shifts `ke0/kg` by +9.2 <!-- claim: gamma_printed_decimal_cost_at_400C_pct --> %
  at 400 °C.

## What pymrm adds

The correlation is algebra and needs no solver. pymrm's operators earn their place in
the *apparatus*: the cylindrical annulus with `nu=1`, solved with a temperature-
dependent conductivity by `pymrm.newton`, wall values read with
`compute_boundary_values`, and the heat rate put on a cell ladder.

- **The reduction formula is biased against its own model, by about ten per cent.** For
  steady radial conduction the conductivity Eq. (19) returns is the average of `ke0`
  uniform in **temperature**, while the abscissa is the **volume**-weighted mean
  temperature, which sits 0.703 <!-- claim: tm_fraction_of_dT_constant_k --> of the
  drop below the hot face. At the five temperature drops Figure 12 shows, 195 <!-- claim: fig12_dT_min_C -->
  to 260 °C <!-- claim: fig12_dT_max_C -->, a bed obeying Eq. (15) exactly plots
  +9.8 <!-- claim: reduction_bias_min_pct --> % to +11.5 <!-- claim: reduction_bias_max_pct --> %
  above the Eq. (15) curve. The 400-cell value is within 1e-3 <!-- claim: reduction_bias_grid_error_n400 -->
  percentage points of its Richardson limit. Swept to a 800 °C drop it reaches
  +21.7 <!-- claim: reduction_bias_swept_max_pct --> %, but that is about three times
  the largest drop the paper ran and the page labels it an extrapolation.
- **"Nearly straight" is a measurement, and it cuts one way harder than the other.**
  A constant-conductivity annulus at this bore ratio departs from a straight line by
  0.137 <!-- claim: nonstraightness_constant_k --> of the temperature drop; with
  `ke0(T)` from Eq. (15) it still departs by up to 0.083 <!-- claim: nonstraightness_varK_at_fig12_max -->
  at the drops Figure 12 shows, straighter by a factor of only 1.66 <!-- claim: straightness_factor_at_fig12_min -->
  at the worst of them. So the remark is strong evidence *against* a constant
  conductivity and only partial evidence *for* Eq. (15).

## Caveats

- **The measurement comparison is weak by construction**, see the finding above.
  Treat the radiation terms as *supported by evidence this page has not examined*,
  not as validated here.
- **The emissivity `p` is never printed in the paper.** `p = 1` is assumed. The
  curve-end comparison is the only quantitative evidence about it. At `p = 1` every
  deviation is positive, the mean 2.4 <!-- claim: fig13to17_curve_ends_mean_abs_pct --> %
  has a standard error of 0.35 <!-- claim: fig13to17_curve_ends_se_pct --> % (the 14 ends treated as independent), and 4 <!-- claim: fig13to17_within_reading_error_p10 -->
  of the readings lie within their own reading error (the sidecar's `+-0.15` in
  `ke0/kg`, on average 1.30 <!-- claim: fig13to17_reading_error_mean_pct --> % of a
  reading). The mean deviation is zero at `p` = 0.935 <!-- claim: fig13to17_p_zero_mean_dev -->,
  and between 0.899 <!-- claim: fig13to17_p_zero_offset_lo --> and 0.970 <!-- claim: fig13to17_p_zero_offset_hi -->
  if every reading carries one offset of up to `+-0.15`; at `p = 0.9` it is -1.28 <!-- claim: fig13to17_mean_dev_pct_at_p09 --> % with 9 <!-- claim: fig13to17_within_reading_error_p09 -->
  readings within their reading error, and at `p = 0.8` it is -4.92 <!-- claim: fig13to17_mean_dev_pct_at_p08 --> %.
  Only a large emissivity deficit is excluded. At `p = 1` both
  radiation coefficients collapse to the same `0.1952(T/100)³`, so every other number
  on the page is blind to the emissivity structure of Eqs. (7) and (8); three break
  rows exist for exactly that.
- **Table 3's calculated column cannot be reproduced from printed material.** At
  `eps = 0.40` it implies `phi` spanning a factor of 7.56 <!-- claim: phi_implied_spread_at_eps040 -->,
  or 2.30 <!-- claim: phi_implied_spread_at_eps040_no_celite --> without the one
  pathological row, while Table 4 assigns that voidage the single value `0.034`.
- **For insulating firebrick, Table 4's printed `ks/kg = 4.48` is the input the
  authors used, not its `ks = 0.20`.** Holding `ks` misses their own Figure 16 curve
  ends by 37.2 <!-- claim: firebrick_curve_ends_mean_abs_pct_ks_const --> % on
  average; holding the ratio, by 3.5 <!-- claim: firebrick_curve_ends_mean_abs_pct_ratio_const --> %.
- **Ranz (1952) is cited, not consulted.** It is not on disk; nothing about it is
  written from memory.

## What the break table cannot detect

Every row of the break table perturbs an input and watches a claim fail, so **none of
them can catch a baseline that is wrong rather than insensitive**. The three exact
agreements above are the example. The defence is the one number computed a second,
independent way: Eq. (15) against the ends of the paper's own printed theoretical
curves.

## Corrections

- 2026-11-09: said of the curve-end comparison that it was "consistent with `p = 1` (+2.4 % mean) and equally with `p = 0.9` (−1.3 %), inconsistent with `p ≲ 0.8` (−4.9 %)",
  and the notebook that the deviations were "well inside it once the ±0.15 reading
  error is allowed for". Measured against the sidecar's reading error, only 4 <!-- claim: fig13to17_within_reading_error_p10 -->
  of the readings at `p = 1` lie within it, all 14 <!-- claim: fig13to17_curve_ends_positive_count -->
  deviations are positive, and the mean of 2.4 <!-- claim: fig13to17_curve_ends_mean_abs_pct --> %
  has a standard error of 0.35 <!-- claim: fig13to17_curve_ends_se_pct --> % (the 14 ends treated as independent). The mean
  deviation is zero at `p` = 0.935 <!-- claim: fig13to17_p_zero_mean_dev -->, and between
  0.899 <!-- claim: fig13to17_p_zero_offset_lo --> and 0.970 <!-- claim: fig13to17_p_zero_offset_hi -->
  if every reading carries one offset of up to `+-0.15`; only a large emissivity
  deficit is excluded.
- 2026-11-09: said "The variable-`k` annulus solve is n-limited." and that "the fixed
  point does not reach tolerance at `n = 1200` without relaxation". The fixed point
  contracts by a factor of at most 0.174 <!-- claim: picard_n1200_contraction_max -->
  per sweep at 1200 cells, as fast as the 0.175 <!-- claim: picard_contraction_max -->
  at 400 cells. What stops it is the stopping rule: its step settles at about 7e-10 <!-- claim: picard_n1200_stall_step_C -->
  degC, the round-off of the linear solve on that grid, and never falls to the
  `1e-10` degC it asks for; with a `1e-8` degC rule it converges at 1200 cells in 10 <!-- claim: picard_n1200_sweeps --> sweeps. The page now solves the
  nonlinear annulus with `pymrm.newton` and keeps the fixed point as a second route,
  which agrees with Newton to 1e-11 <!-- claim: picard_vs_newton_max_dT --> of the
  temperature drop.

## Reuse

`phi` is the parameter that will bite you. If you need a number, Table 4's printed
pairs are the defensible source: `0.034` at `eps = 0.40`, `0.040` at `0.43`, `0.050`
at `0.50`, `0.060` at `0.54`, which the exponential fit on the page passes through
within the residual it prints. If you re-use Figures 13 to 17 as data, subtract about
ten per cent for the reduction bias.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: all six datasets are committed CSVs in `data/`.
