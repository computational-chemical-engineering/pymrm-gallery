# F1.9: Zahradnik et al., the homogeneous to heterogeneous regime transition

The bubble-flow models of Zahradnik, Fialova, Ruzicka, Drahos, Kastanek and Thomas,
"Duality of the gas-liquid flow regimes in bubble column reactors", *Chem. Eng. Sci.*
**52** 3811-3826 (1997),
[doi:10.1016/S0009-2509(97)00226-1](https://doi.org/10.1016/S0009-2509(97)00226-1): the
slip-velocity model of the homogeneous branch with the six hindrance functions of Table
2 (eqs. 3 to 6), Zuber-Findlay for the heterogeneous branch (eqs. 7 to 12) and the drift
velocities of eqs. 13 and 14. The page asks whether either model, or the two together,
fixes the printed transition (`0.04` m/s, `0.210`), and checks the printed numbers
against each other and against Figs. 3(a), 11(a) and 12(b). A correlation page; it calls
no pymrm operator.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/zahradnik-1997-printed-numbers.csv` | every constant, limit and printed result the page uses, as printed, with a short quote |
| `data/zahradnik-1997-fig11a.csv` | slip-velocity ratio against holdup, plate A, and the drawn Table 2 lines, digitised from Fig. 11(a) |
| `data/zahradnik-1997-fig3a.csv` | holdup against gas velocity, plates A and B, digitised from Fig. 3(a) |
| `data/zahradnik-1997-fig12b.csv` | the drift-flux plot behind eq. 10, plate B, digitised from Fig. 12(b) |
| `data/digitise.py` | the readings of the three figures and their calibration; writes the figure CSVs |
| `data/check_values.csv` | the page's functions at one test condition, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it.

## Things to know before reusing this

- Neither fitted model fixes the transition. The homogeneous branch's flux maximum is
  the only point it fixes by itself; for Griffith-Wallis it lies at a gas velocity of
  0.05875 <!-- claim: r2_uGmax_gw --> m/s, and the printed velocity is 0.6809 <!-- claim: r2_ratio_u_gw -->
  of it (0.8641 <!-- claim: r2_ratio_u_rz --> for Richardson-Zaki; the other four forms have
  no maximum). The two branches never meet at the printed constants: 0 <!-- claim: r3_n_forms_crossing -->
  of the 6 forms cross eq. 10, because the homogeneous slip velocity never reaches the
  eq. 10 intercept. With the other single-bubble velocities the paper prints, the
  crossings sit at holdups of 0.0354 <!-- claim: r3_cross_eps_max_267 --> or less.
- The change of regime is a jump in holdup at fixed gas velocity: at `0.04` m/s the
  Griffith-Wallis holdup is 1.837 <!-- claim: r3_jump_ratio --> times that of eq. 10.
- The printed pair is a property of the data: it is the last Fig. 11(a) marker within
  three reading errors of the Griffith-Wallis line (marker 7 <!-- claim: r4_last_on -->,
  0.6 <!-- claim: r4_hu_eps_last --> half-units from `0.210`), and the next marker is
  7.7 <!-- claim: r4_z_first_off --> reading errors off. On the line itself, `0.210` sits at
  0.03899 <!-- claim: r1_uG_gw --> m/s, which rounds to `0.04`; the same line at exactly
  `0.04` m/s gives 0.2175 <!-- claim: r1_eps_gw -->. The printed digits do not select the
  form: 4 <!-- claim: r1_n_forms_round_to_ucrit --> of the 6 forms put `0.210` at a velocity
  that rounds to `0.04`. The fit to the markers does: Griffith-Wallis has an rms residual
  of 0.0078 <!-- claim: r4_rms_gw --> on 7 <!-- claim: r4_n_homog --> homogeneous markers, against
  0.0343 <!-- claim: r4_rms_marrucci --> for the runner-up.
- The pair holds for air-water, a column of `0.15` m at an aspect ratio of `10` and plate A only.
  The plate A data rejoin eq. 10 (fitted on plate B) from `0.125` m/s, within
  2.78 <!-- claim: r4_z_rejoin_max --> standard deviations of the plate B scatter.
- In the airlift riser, eq. 6 inserted in eq. 3 is Zuber-Findlay with a distribution
  parameter of one, an identity; the first printed slope break then holds a holdup of
  0.1593 <!-- claim: r5_eps6_1 -->, beyond the printed limit `0.15`; the verdict would flip
  only at a riser terminal velocity of 0.2679 <!-- claim: r5_ub_flip --> m/s, higher than
  every single-bubble velocity the paper prints.
- Eq. 14 cannot give the printed `0.267` m/s for water at any bubble size: its supremum is
  0.2529 <!-- claim: r3_ub_sup --> m/s with corrected groups; the groups as printed in the
  eq. 14 block are not dimensionless.
- No temperature is printed; eqs. 13 and 14 use water at 20 C as an input of this page.
