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
  somewhere in the paper's velocity range for all 12 <!-- claim: n_fitted_maxdev_gt_10 --> of the
  12 <!-- claim: n_fitted --> fitted systems (9 <!-- claim: n_fitted_maxdev_gt_10_beyond_rounding --> by more
  than the rounding of the printed $Mo$ exponent moves the ratio); air-water clears
  `10 %` by only 0.0032 <!-- claim: maxdev_air_water_excess -->. The ratio of Akita-Yoshida to eq. 8 is a
  velocity-shape factor common to all systems, peaking at
  1.0545 <!-- claim: shape_factor -->, times a property factor: 1.2671 <!-- claim: P_H2_water --> for H2 in
  water, 0.8285 <!-- claim: P_air_nbutanol --> for n-butanol. Neither formula has a
  column diameter.
- The printed window (bounds taken to half a unit of their last digit) excludes
  air-methanol, a fitted system: its gas-to-liquid viscosity ratio is 1.511 <!-- claim: mr_methanol_over_bound -->
  times the upper bound. H2-water is inside only through that half unit: its
  $\rho_G/\rho_L$ is -0.2645 <!-- claim: halfunits_rr_lo --> half-units from the printed lower bound
  `8.4e-5` (positive is inside). The window excludes water systems below a gas velocity of
  0.0793 <!-- claim: u_enter_Ca_window_air_water --> m/s (air-water, strict bound), although the experiments start
  lower. It is narrower than the plotted fitting data at both ends: 1 <!-- claim: f3_n_below_Ca --> of
  the 33 <!-- claim: f3_n_assigned --> assigned Figure 3 markers has a $Ca$ below the printed bound, and
  2 <!-- claim: f3_n_above_umax --> H2 markers plot where eq. 8 needs a velocity beyond the highest of
  the experiments, which the paper does not explain.
- In Figure 2, at matched velocity, H2 lies -0.158 <!-- claim: f2_H2_air --> ±
  0.011 <!-- claim: f2_H2_air_se --> in ln(holdup) from air, where Akita-Yoshida (no gas
  effect) gives 0. Eq. 8 puts H2 at 0.7892 <!-- claim: gasfactor_H2_water --> of air, so the markers show
  a smaller drop: a misfit of 8.222 % <!-- claim: f2_H2_air_misfit --> in holdup, inside eq. 8's printed `10 %`
  maximum deviation; the ± is about one reading error of a marker, the misfit several. H2 and the 5:1 mixture miss eq. 8's gas factor
  (t = +7.3 <!-- claim: f2_H2_air_t8 --> and -3.8 <!-- claim: f2_51_air_t8 -->); the 1:1 mixture follows it
  (t = -0.4 <!-- claim: f2_11_air_t8 -->). The miss is in the ratio between gases.
- In Figure 3, eq. 8 is closer than Akita-Yoshida for H2 by 0.1315 <!-- claim: f3_pair_H2_water --> ±
  0.0231 <!-- claim: f3_pair_se_H2_water --> in absolute log ratio, and by 0.1592 <!-- claim: f3_pair_H2_water_in_range --> ±
  0.0181 <!-- claim: f3_pair_se_H2_water_in_range --> when the 2 H2 markers beyond that velocity are left out; the ± are
  scatter only (a 10 px shift of every H2 marker moves the mean by 0.03402 <!-- claim: f3_pair_H2_water_shift10px --> at the worst of four directions).
  Both figures are fitting data of eq. 8: this is how eq. 8 represents its own data, not
  validation.
- Figures 4 (held-out alcohols), 12 and 13 (electrolytes) are not digitised; the page
  says nothing about the held-out data.
