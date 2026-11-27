# F1.2: Hikita et al., gas holdup with liquid and gas properties

The gas-holdup correlation of Hikita, Asai, Tanigawa, Segawa and Kitao, "Gas hold-up in
bubble columns", *Chem. Eng. J.* **20** 59-67 (1980),
[doi:10.1016/0300-9467(80)85006-4](https://doi.org/10.1016/0300-9467(80)85006-4): eq. 8,
four dimensionless groups including the gas-to-liquid density and viscosity ratios, and
eq. 10 for electrolyte solutions. The page sets eq. 8 beside Akita-Yoshida on the
paper's own Table 2 systems, checks the printed window against those systems and the
plotted points, and reads Figures 2 and 3. A correlation page; it calls no pymrm
operator.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/hikita-1980-printed-numbers.csv` | every constant, window bound, range and printed deviation the page uses, as printed |
| `data/hikita-1980-table2-properties.csv` | Table 2, the properties of the systems, as printed (CGS) |
| `data/hikita-1980-fig2.csv` | holdup against gas velocity for four gases in water, digitised from Figure 2 |
| `data/hikita-1980-fig3.csv` | observed against calculated holdup of the fitted systems, digitised from Figure 3 |
| `data/digitise.py` | the readings of both figures and their calibration; writes the two figure CSVs |
| `data/check_values.csv` | the page's functions at one test condition, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it. Akita and Yoshida's
constants and water at 20 C are loaded from
[F1.1](../F1.1-akita-yoshida-bubble-column/), and one check value from
[F1.8](../F1.8-shah-bubble-column-kla/).

## Things to know before reusing this

- Eq. 8 and Akita-Yoshida differ by more than the maximum deviation printed for eq. 8
  somewhere in the paper's velocity range for 12 <!-- claim: n_fitted_maxdev_gt_10 -->
  of the fitted systems (9 <!-- claim: n_fitted_maxdev_gt_10_beyond_rounding --> by more
  than the rounding of the printed $Mo$ exponent moves the ratio). The ratio is a
  velocity-shape factor common to all systems, peaking at
  1.0545 <!-- claim: shape_factor -->, times a property factor: 1.2671 <!-- claim: P_H2_water -->
   for H2 in water, 0.8285 <!-- claim: P_air_nbutanol --> for
  n-butanol. Neither formula has a column diameter.
- The printed window excludes air-methanol, a fitted system: its gas-to-liquid viscosity
  ratio is 1.511 <!-- claim: mr_methanol_over_bound --> times the upper bound. It
  excludes water systems below a gas velocity of 0.0793 <!-- claim: u_enter_Ca_window_air_water -->
   m/s (air-water, strict bound), although
  the experiments start lower, and 3 <!-- claim: f3_n_outside_window --> of the
  assigned Figure 3 markers lie outside it.
- In Figure 2, at matched velocity, H2 lies -0.158 <!-- claim: f2_H2_air --> ±
  0.011 <!-- claim: f2_H2_air_se --> in ln(holdup) from air: Akita-Yoshida (no gas
  effect) is far off, and eq. 8's gas factor is off by t = 7.3 <!-- claim: f2_H2_air_t8 -->
  . In Figure 3, eq. 8 is closer than Akita-Yoshida for H2 by
  0.120 <!-- claim: f3_pair_H2_water --> ± 0.025 <!-- claim: f3_pair_se_H2_water --> in
  absolute log ratio. Both figures are fitting data of eq. 8: this is how eq. 8
  represents its own data, not validation.
- Figures 4 (held-out alcohols), 12 and 13 (electrolytes) are not digitised; the page
  says nothing about the held-out data.
