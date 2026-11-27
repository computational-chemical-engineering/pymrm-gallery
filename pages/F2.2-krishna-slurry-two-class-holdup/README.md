# F2.2: Krishna et al. (1997), the two-class holdup split in slurry bubble columns

From Krishna, R., de Swart, J. W. A., Ellenberger, J., Martina, G. B. & Maretto, C.,
*Gas holdup in slurry bubble columns: effect of column diameter and slurry
concentrations*, *AIChE Journal* **43**(2) 311-316 (1997),
[doi:10.1002/aic.690430204](https://doi.org/10.1002/aic.690430204).

- **Structures:** `S7` (multiphase coupling)
- **Runtime:** a few seconds.

The page tests the paper's split of the gas holdup, a scale-independent dense-phase
voidage (eq. 4) and a slurry-independent, diameter-dependent large-bubble holdup
(eq. 6), against the paper's own figures, and works out what the split implies for a
first-order gas-consuming reaction (rate constant illustrative). The reactor numbers
are consequences of the correlations, never measurements.

## Findings

- **Dense phase.** At `eps_s = 0.16`, the only concentration run in both columns, the lowest visible `0.10 m` centre in Fig. 7 sits 0.0198 <!-- claim: dense_016_gap_D010_minus_fig6 --> in eps_df from the `0.38 m` Fig. 6 mean, 0.0149 <!-- claim: dense_016_margin --> clear of the read error and two standard errors: the paper's "virtually independent of the column diameter" is not borne out at this concentration, unless `0.10 m` circles are hidden in the dark `0.38 m` diamond stack. For scale, the `0.10 m` column's own `eps_s = 0.25` group misses eq. 4 by 0.0207 <!-- claim: dense_025_D010_eq4_miss -->, so a group can sit this far off within one column, and the figures do not show that the diameter causes the gap. The `0.38 m` Fig. 6 means sit below eq. 4 in 4 <!-- claim: fig6_series_below_eq4 --> of 4 series.
- **Dilute phase.** Relative to the `0.10 m` markers (`eps_s` 0, `0.05` and `0.10` pooled, 2 <!-- claim: fig8_n_D010_es0 -->, 7 <!-- claim: fig8_n_D010_es005 --> and 8 <!-- claim: fig8_n_D010_es010 --> markers), the large-bubble amplitude is 1.104 <!-- claim: dilute_ratio_D019_low --> at `0.19 m` and 0.880 <!-- claim: dilute_ratio_D038_low --> at `0.38 m`, where `N = 0.18` predicts 0.891 <!-- claim: N018_factor_D019 --> and 0.786 <!-- claim: N018_factor_D038 -->; with only the eps_s 0 circles as reference the ratios are 1.108 <!-- claim: dilute_ratio_D019_es0ref --> and 0.884 <!-- claim: dilute_ratio_D038_es0ref -->, so the direction stands. No single power of `D_T` fits the three columns (chi-square 59 <!-- claim: dilute_free_N_chi2 --> on one degree of freedom; the formal weighted fit N = 0.11 <!-- claim: dilute_free_N --> is not a fit): pair by pair the exponent is -0.154 <!-- claim: dilute_pair_N_010_019 -->, 0.096 <!-- claim: dilute_pair_N_010_038 --> and 0.327 <!-- claim: dilute_pair_N_019_038 -->, and on the abscissa range the three series share the formal fit gives 0.30 <!-- claim: dilute_free_N_shared_range -->, so a single N depends on the range chosen. The standard errors are of the marker scatter: the `eps_s` `0.05` and `0.10` series of the `0.10 m` reference differ by 11.9 % <!-- claim: dilute_D010_subseries_gap --> (4.7 <!-- claim: dilute_D010_subseries_z --> marker standard errors). With the series as the unit of scatter the `0.19 m` ratio is 1.4 <!-- claim: dilute_series_z_vs_1_D019_low --> from 1 and 3.0 <!-- claim: dilute_series_z_vs_N018_D019_low --> from 0.891 <!-- claim: N018_factor_D019 --> series-level standard errors, the `0.38 m` ratio -1.8 <!-- claim: dilute_series_z_vs_1_D038_low --> from 1. At `0.38 m` the high-slurry amplitude is 1.152 <!-- claim: slurry_ratio_D038 --> times the low-slurry one: against 1 that is +6.8 <!-- claim: slurry_ratio_z_vs_1_D038 --> marker-level and 1.9 <!-- claim: slurry_series_z_vs_1_D038 --> series-level standard errors; against the paper's own alpha2 step of 1.119 <!-- claim: alpha2_ratio --> it is +1.4 <!-- claim: slurry_ratio_z_vs_alpha2_D038 --> marker-level and 0.4 <!-- claim: slurry_series_z_vs_alpha2_D038 --> series-level. "Virtually independent of slurry concentration" is the paper's wording while its own correlation carries that step. The isolated Fig. 8b markers sit 4.4 % <!-- claim: slurry_selected_over_band_D038 --> higher than all 99 <!-- claim: fig8_n_band_b --> Fig. 8b markers, the size of the selection effect.
- **Printed constants.** The Fig. 5 lines for `eps_s` 0 and `0.05` are -16.6 % <!-- claim: fig5_gap_eq2_es0 --> and 17.5 % <!-- claim: fig5_gap_eq2_es005 --> off eq. 2; the Fig. 7, 8 and 9 lines are eqs. 4, 6 and 7 as printed.
- **Reactor consequence.** At `D_T = 0.38 m` and `U = 0.30 m/s` the slurry halves the holdup and cuts the low-conversion capacity 2.61 <!-- claim: capacity_ratio_k0 -->-fold, but raises the k H0/U needed for `X = 0.9` only 1.236 <!-- claim: k_ratio_T_T0_stirred -->-fold (rate constant illustrative).

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | the correlations (eqs. 1 to 7), the closed-form two-class reactor (`split`, `conversion`, `d_required`, `crossing`) and `TwoClassReactor` (the reactor on pymrm operators, one `newton` step) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/krishna-1997-printed-numbers.csv` | the printed constants and Table 1, as strings |
| `data/krishna-1997-fig*.csv` | pixel reads of Figs. 5 to 9 with their values, each with a provenance sidecar holding the tick calibration and read error |
| `data/check_values.csv` | independent check values of the page's functions |
