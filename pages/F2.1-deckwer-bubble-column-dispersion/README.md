# F2.1: Axial dispersion and kLa in tall bubble columns

From Deckwer, W.-D., Burckhart, R. and Zoll, G., *Mixing and mass transfer in tall
bubble columns*, *Chemical Engineering Science* **29** 2177-2188 (1974),
[doi:10.1016/0009-2509(74)80025-4](https://doi.org/10.1016/0009-2509(74)80025-4).

- **Structures:** `S7` (gas-liquid coupling), `S4` (1-D convection, dispersion and exchange)
- **Runtime:** see `meta.yaml`.

## What the page does

The paper reads kLa from measured liquid oxygen profiles by fitting a liquid-only axial
dispersion model with a prescribed gas composition, and defends that evaluation with
three statements: S1, a dispersion coefficient off by "up to 100 per cent" leaves
Stanton numbers that "do not differ more than 3 per cent"; S2, the plug-flow NTU evaluation is low by
an amount that "does not exceed more than 20 per cent"; S3, a gas composition varying
along the column raises kLa by "up to 20 per cent". The page tests each at the paper's
own operating points (Figs. 9 and 10) against the scatter of the paper's own Table 1
correlation, and solves with pymrm the coupled liquid-dispersion and gas plug-flow
model that the paper named and set aside.

## Findings

- **S1 holds at the measuring ports in 17 <!-- claim: s1_D1_count_within_3pct --> of
  20 <!-- claim: s1_D1_case_total --> cases, not in all.** With D halved or doubled and St
  refitted on the ports of Fig. 10, three cases lie outside the paper's bound: the Fig. 9
  BC I water point with D doubled, 15.3 % <!-- claim: s1_D1_dev_out1 -->, and the Fig. 9
  BC I Na2SO4 point with D doubled, 12.7 % <!-- claim: s1_D1_dev_out2 -->, whose profiles
  are flat after the first tenth of the column, and the Fig. 10 BC I point at the lower
  liquid flow with D halved, 5.3 % <!-- claim: s1_D1_dev_out3 -->. The count is fragile to the input
  reads: the nearest case inside the bound reads 2.99 % <!-- claim: s1_D1_nearest_dev -->,
  and 3.03 % <!-- claim: s1_D1_nearest_dev_epsp01 --> with every gas holdup raised by the step the page
  prints,
  which leaves 16 <!-- claim: s1_D1_count_within_3pct_epsp01 --> within. On an evenly spaced
  design the largest change is 24.2 % <!-- claim: s1_D2_worst_dev -->: S1 depends on where
  the ports are.
- **S2 holds for Fig. 11 as drawn; it is not a general bound.** Fig. 11's NTU line lies at
  0.840 <!-- claim: s2_fig11_ratio_min --> to 0.878 <!-- claim: s2_fig11_ratio_max --> of the
  BC II line, inside "20 per cent". The dispersion model gives an NTU St of
  0.812 <!-- claim: s2_ntu_ratio_VL2p5_min --> to 0.908 <!-- claim: s2_ntu_ratio_VL2p5_max -->
  times the dispersion St at the BC II water flow of Fig. 10, but
  0.464 <!-- claim: s2_ntu_ratio_VL0p8_min --> to 0.622 <!-- claim: s2_ntu_ratio_VL0p8_max -->
  at the flow of Fig. 9, beyond "20 per cent" and Table 1's
  14.6 % <!-- claim: t1_scatter_BCII_water --> scatter; the NTU St falls 20 per cent below the dispersion St at
  1.46 <!-- claim: s2_VL_ratio0p8_uG0p6 --> to 2.38 <!-- claim: s2_VL_ratio0p8_uG1p5 --> m3/h.
  At the printed BC II operating points the ratio is
  0.611 <!-- claim: s2_ntu_ratio_printed_mol_min --> to
  0.686 <!-- claim: s2_ntu_ratio_printed_mol_max --> for the molasses runs of Fig. 9 and
  0.918 <!-- claim: s2_ntu_ratio_printed_H2O --> for the water run of Fig. 10. Matching each
  read ratio of Fig. 11 needs 1.73 <!-- claim: s2_implied_VL_uG0p6 --> to
  3.38 <!-- claim: s2_implied_VL_uG1p5 --> m3/h, rising with u_G while Fig. 11's ratio rises
  too: no single liquid flow reproduces the line, and the paper does not print the flow of
  its NTU runs. A stirred-tank evaluation gives
  2.63 <!-- claim: s2_cstr_ratio_min --> times the dispersion St or more: "too high", as
  p. 2180 says.
- **S3 holds.** The gas mole fraction changes by 0.0150 <!-- claim: s3_dyt_F10BCII --> at the most affected point, inside the paper's "value of 0.04". Eq. 7 raises St
  over constant y by at most 13.4 % <!-- claim: s3_eq7_gain_max --> at the four points,
  inside "up to 20 per cent". Against the coupled model, constant y gives
  0.831 <!-- claim: s3_fit_cy_F10BCII --> of the true St at the BC II water point of Fig. 10,
  beyond the 14.6 % <!-- claim: t1_scatter_BCII_water --> scatter; eq. 7, used for Fig. 10,
  gives 0.943 <!-- claim: s3_fit_e7_F10BCII -->, beyond it at
  0 <!-- claim: s3_count_e7_beyond_scatter --> of 4 points. The constant-y margin is smallest
  at 6 m3/h, 15.2 % <!-- claim: s3_BCII_cy_bias_VL6 -->, and
  14.5 % <!-- claim: s3_BCII_cy_bias_VL6_csat0p95 --> with c*_sat lowered; the bias is beyond
  the scatter at all 13 <!-- claim: s3_BCII_scan_cy_beyond_scatter --> scan points of V_L, and at
  12 <!-- claim: s3_BCII_scan_cy_beyond_scatter_csat0p95 --> of them with c*_sat lowered, so "beyond
  the scatter at every liquid flow" depends on the saturation value. With n_G0 at the top
  pressure instead of the entry pressure the constant-y ratio is
  0.775 <!-- claim: s3_sens_ptop_cy_F10BCII --> and eq. 7's gain
  18.3 % <!-- claim: s3_sens_ptop_gain_max -->, still inside "up to 20 per cent".

These are self-consistency tests of the evaluation, not measurements: the profiles of
Figs. 9 and 10 are not digitised, and Table 1 is the paper's regression of its own fits.

## Assumptions

- The feed of every BC I absorption point is an assumed value (the page's Parameters
  section): Fig. 10's two open-circle feed marks cannot be assigned to one curve by marker
  shape. Reading the larger mark instead moves the largest St change of S1 at the ports of
  Fig. 10 to 16.1 % <!-- claim: s1_D1_worst_dev_feed1p05 --> (15.3 % <!-- claim: s1_D1_worst_dev -->
  as read) and on the evenly spaced design to 25.0 % <!-- claim: s1_D2_worst_dev_feed1p05 -->
  (24.2 % <!-- claim: s1_D2_worst_dev -->), and leaves the counts inside the paper's bound at
  17 <!-- claim: s1_D1_count_within_3pct_feed1p05 --> and
  13 <!-- claim: s1_D2_count_within_3pct_feed1p05 -->.
- The feed of the BC II water desorption point is read at the square mark on Fig. 10's
  ordinate.
- H is not printed: it is set from the oxygen saturation of water under air (Benson and
  Krause 1984); the page measures the effect of changing that saturation value.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | the groups of eqs. 1, 4, 5, 8; `adm_profile` (the closed forms A6 to A20, misprints corrected, printed forms as options); `fit_stanton`; `CoupledColumn` (liquid dispersion and gas plug flow on pymrm operators, one Newton system); `coupled_bvp` (the same model by `solve_bvp`) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/deckwer-1974-*.csv` | printed numbers, Table 1, the operating points, holdup reads, port positions, feed marks and Fig. 11 reads, each with its provenance sidecar |
| `data/check_values.csv` | values of the page's functions at one test condition from a second route, asserted before any claim |
