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

- **Dense phase.** At `eps_s = 0.16` the `0.10 m` column of Fig. 7 sits 0.0198 <!-- claim: dense_016_gap_D010_minus_fig6 --> in eps_df from the `0.38 m` Fig. 6 mean, a margin of 0.0149 <!-- claim: dense_016_margin --> after the read error and two standard errors: the paper's "virtually independent of the column diameter" does not hold at the one concentration where the figures separate the diameters. The `0.38 m` Fig. 6 means sit below eq. 4 in 4 <!-- claim: fig6_series_below_eq4 --> of 4 series.
- **Dilute phase.** Relative to `0.10 m`, the large-bubble amplitude is 1.104 <!-- claim: dilute_ratio_D019_low --> at `0.19 m` and 0.880 <!-- claim: dilute_ratio_D038_low --> at `0.38 m`, where `N = 0.18` predicts `0.891` and `0.786`; the three columns give N = 0.11 <!-- claim: dilute_free_N -->. At `0.38 m` the high-slurry amplitude is 1.152 <!-- claim: slurry_ratio_D038 --> times the low-slurry one, against 1 for "virtually independent" and the paper's own alpha2 ratio.
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
