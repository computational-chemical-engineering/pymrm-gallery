# F2.1: Axial dispersion and kLa in tall bubble columns

From Deckwer, W.-D., Burckhart, R. and Zoll, G., *Mixing and mass transfer in tall
bubble columns*, *Chemical Engineering Science* **29** 2177-2188 (1974),
[doi:10.1016/0009-2509(74)80025-4](https://doi.org/10.1016/0009-2509(74)80025-4).

- **Structures:** `S7` (gas-liquid coupling), `S4` (1-D convection, dispersion and exchange)
- **Runtime:** see `meta.yaml`.

## What the page does

The paper reads kLa from measured liquid oxygen profiles by fitting a liquid-only axial
dispersion model with a prescribed gas composition, and defends that evaluation with
three statements: S1, a dispersion coefficient off by "up to 100 per cent" moves the
Stanton number by "not more than 3 per cent"; S2, the plug-flow NTU evaluation is low by
an amount that "does not exceed more than 20 per cent"; S3, a gas composition varying
along the column raises kLa by "up to 20 per cent". The page tests each at the paper's
own operating points (Figs. 9 and 10) against the scatter of the paper's own Table 1
correlation, and solves with pymrm the coupled liquid-dispersion and gas plug-flow
model that the paper named and set aside.

## Findings

- **S1 holds at the measuring ports except where the profile saturates early.** With D
  halved or doubled and St refitted on the ports of Fig. 10, St stays within the paper's
  bound in 17 <!-- claim: s1_D1_count_within_3pct --> of the cases. The largest change
  is 15.3 % <!-- claim: s1_D1_worst_dev -->, at the Fig. 9 BC I points with D doubled,
  whose profiles are flat after the first tenth of the column; one case sits at the
  bound itself, so the count is fragile to the input reads. On an evenly spaced design
  the largest change is 24.2 % <!-- claim: s1_D2_worst_dev -->: S1 depends on where the
  ports are.
- **S2 does not hold at the liquid flow of Fig. 9.** For desorption in BC II the NTU St
  is 0.464 <!-- claim: s2_ntu_ratio_VL0p8_min --> to 0.622 <!-- claim: s2_ntu_ratio_VL0p8_max -->
  times the dispersion St at the liquid flow of Fig. 9, and
  0.812 <!-- claim: s2_ntu_ratio_VL2p5_min --> to 0.908 <!-- claim: s2_ntu_ratio_VL2p5_max -->
  at the BC II water flow of Fig. 10. The ratios of Fig. 11's dashed NTU line imply
  liquid flows of 1.73 <!-- claim: s2_implied_VL_uG0p6 --> to 3.38 <!-- claim: s2_implied_VL_uG1p5 -->
  m3/h for the NTU runs, which the paper does not print. A stirred-tank evaluation gives
  2.63 <!-- claim: s2_cstr_ratio_min --> times the dispersion St or more: "too high", as
  p. 2180 says.
- **S3: the gas phase is the larger effect for BC II desorption.** Against the coupled
  model, constant-y evaluation gives 0.831 <!-- claim: s3_fit_cy_F10BCII --> of the true
  St at the BC II water point of Fig. 10, beyond Table 1's scatter, and stays beyond it
  over the whole liquid-flow range at that gas velocity. Eq. 7 recovers it to
  0.943 <!-- claim: s3_fit_e7_F10BCII -->, a gain of 13.4 % <!-- claim: s3_eq7_gain_max -->.
  The largest change of the gas mole fraction, 0.0150 <!-- claim: s3_dyt_F10BCII -->,
  is inside the paper's bound.

These are self-consistency tests of the evaluation, not measurements: the profiles of
Figs. 9 and 10 are not digitised, and Table 1 is the paper's regression of its own fits.

## Assumptions

- The feed of every BC I absorption point is an assumed value (the page's Parameters
  section): Fig. 10's two open-circle feed marks cannot be assigned to one curve by
  marker shape.
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
