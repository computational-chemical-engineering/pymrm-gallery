# C2.1: Steam methane reforming, the Xu-Froment intrinsic kinetics

Three rate equations from 1989 that almost every published reformer simulation
still uses. This page rebuilds them from the printed tables and puts them back
against the measurements they were fitted to.

- **Structures:** `S1` (0D/1D reaction network), `S2` (plug flow with reaction)
- **Reference:** Xu & Froment (1989), AIChE J 35(1) 88-96, doi:10.1002/aic.690350109
- **Runtime:** about 10 s
- **Files:** `index.py` (the page, paired with `index.ipynb`), `model.py`
  (`XuFroment` kinetics, `SteamReformer` plug flow, the Radau reference and the
  round trip), `test_model.py`

## Agreement

0.0017 <!-- claim: mad_conversion --> mean absolute deviation in conversion over
61 digitised measurements at four temperatures (773-848 K), worst point
0.0058 <!-- claim: max_dev_conversion -->, which is
2.7 <!-- claim: mad_rel_mean_pct --> % of the mean measured conversion. The
digitisation error is about 0.0006 <!-- claim: digitisation_error_conversion -->. **Nothing was fitted**: every kinetic
parameter is Table 6 as printed.

The model is biased high, and the bias in methane conversion grows with space
time. At 823 and 848 K a 5 % shift in the equilibrium constants moves the
predictions at the long-space-time measurements by
0.0025 <!-- claim: k_halfspan_long_823K --> and
0.0035 <!-- claim: k_halfspan_long_848K -->, the size of the bias there
(0.0035 <!-- claim: bias_long_823K --> and 0.0023 <!-- claim: bias_long_848K -->);
the equilibrium constants are the one ingredient that does not come from this
paper. At 773 and 798 K the same shift moves them by only
0.0005 <!-- claim: k_halfspan_long_773K --> and
0.0014 <!-- claim: k_halfspan_long_798K --> against biases of
0.0028 <!-- claim: bias_long_773K --> and 0.0039 <!-- claim: bias_long_798K -->,
so there the residual is not explained by them.

The discretisation does not enter: first-order upwind on 1000 cells meets an
independent Radau integration of the same ODEs to
4.1e-5 <!-- claim: disc_error_max_n1000 --> at the measured points,
2.4 <!-- claim: disc_error_share_of_mad_pct --> % of the mean deviation.

## Data

Two datasets, both with provenance sidecars.

`xu-froment-1989-parameters.csv`: Tables 5 and 6. **The PDF text layer mangles
every exponent** (`8.664 lo-'` for 8.664e-7), so these were read from a 600 dpi
render of journal page 94 and none was repaired by inference. The page reads the
CSV a second time as text, so the round trip scores every entry against the
precision it was printed with.

`xu-froment-1989-conversion.csv`: the experimental markers of Figures 2 and 3.
The plotted curves are the authors' own model and were deliberately not
extracted. The sidecar records the cross-figure pairing and carbon-closure
checks of the two independent digitisations.

## A summary that cannot see the error it is meant to catch

The worst per-row round-trip deviation is
1.4432 <!-- claim: param_round_trip_worst_pct --> %. It is a `max()` over seven
rows, and its baseline scatter comes from the *printed precision of the
activation energies*: E is given to 0.1 kJ/mol, which at T_ref = 648 K is worth
0.93 <!-- claim: e_precision_share_pct --> % on A by itself. A single-digit slip
landing under that ceiling leaves it unchanged; a transposition in `k3`, the row
that sets the maximum, moves it down. Dividing each row by its own printed
precision and summing the seven ratios
(10.438 <!-- claim: param_round_trip_ratio_sum -->) gives a statistic every row
reaches: it fails its check for every injected slip except a one-unit change in
the fourth figure of `k1`'s A, which no round trip built from these digits can
resolve. Compare the sum, not the max. The notebook's break table has every row.

## Three traps in this paper

1. **Split reference temperature.** T_ref = 648 K for k1, k2, k3, K_CO, K_H2 but
   **823 K** for K_CH4 and K_H2O. Forcing 648 K on all seven puts K_CH4 out by a
   factor of 0.221 <!-- claim: tref_slip_factor_K_CH4 --> and K_H2O by
   33 <!-- claim: tref_slip_factor_K_H2O -->.
2. **Three activity levels.** Tables 5 and 6 are the *steam-reforming reference*
   level (partially deactivated catalyst), which is what the Figure 2/3 data are,
   so no correction. Multiply the rate coefficients by **1.225** to reproduce the
   paper's Figures 4 and 5, and by **2.246** for fresh catalyst. **Neither factor
   is validated here**: the runs they belong to are chemistry this page does not
   model. The page executes them and says so.
3. **No K_CO2 term.** It was never statistically significant and correctly has
   no term in the denominator. Adding one silently changes the model.

Also: the paper does not tabulate the equilibrium constants K1, K2, K3 at all.
They must come from elsewhere; the notebook checks the correlation it uses
against the paper's own Table 7 reaction enthalpies (worst
1.93 <!-- claim: equilibrium_vs_table7_worst_pct --> %). **That check compares
slopes only.** A +1 slip in either intercept multiplies that K by e at every
temperature and leaves it unchanged. What catches it is the comparison against
the 61 measurements, where both slips fail every deviation claim.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/C2.1-xu-froment-smr/index.py   # after editing index.py
python scripts/run_kernels.py pages/C2.1-xu-froment-smr
GALLERY_OFFLINE=1 python -m pytest -q pages/C2.1-xu-froment-smr/test_model.py
```

No PDF is needed to rebuild: both datasets are committed CSVs. To regenerate the
digitised dataset you need your own lawfully obtained copy of the PDF (it is not
in the repository, and must not be):

```bash
python pages/C2.1-xu-froment-smr/extract_figures.py ~/papers/pymrm-gallery/"Xu1989-methane-steam-reforming-kinetics-AIChEJ35-88.pdf"
```

That script is deterministic (the candidate indices that survived the visual
audit are constants in it) and re-runs the monotonicity, cross-figure pairing
and carbon-closure checks, exiting non-zero if any of them fails.

## Corrections

- 2026-11-02: said the discretisation error at 1000 cells is 2.1e-5, 1 % of the
  mean deviation from the data. That was one point (x_CH4 at 848 K, W/F = 0.40),
  estimated from the grid ladder; against an independent Radau integration at
  every measured point it is 4.1e-5 <!-- claim: disc_error_max_n1000 -->, 2.4
  <!-- claim: disc_error_share_of_mad_pct --> % of the mean deviation.
- 2026-11-02: said a 5 % shift in the equilibrium constants moves the
  long-space-time predictions by 0.0026 to 0.0039, the size of the bias, so the
  residual there cannot be charged to the kinetics. That was computed at 823 and
  848 K only, at W/F = 0.35, on a coarser grid than the model it was compared
  with. At the measured points it holds at 823 and 848 K, but at 773 and 798 K
  the shift is 0.0005 <!-- claim: k_halfspan_long_773K --> and 0.0014
  <!-- claim: k_halfspan_long_798K --> against biases of 0.0028
  <!-- claim: bias_long_773K --> and 0.0039 <!-- claim: bias_long_798K -->.
- 2026-11-02: said reactions I and III are 85-97 % of the way to equilibrium at
  the longest space times. Those are the 848 K values of reaction I between W/F = 0.28 and
  0.45; at the measured points above W/F = 0.28, beta is 0.57
  <!-- claim: beta_I_long_min --> to 0.90 <!-- claim: beta_I_long_max --> for
  reaction I and 0.45 <!-- claim: beta_III_long_min --> to 0.88
  <!-- claim: beta_III_long_max --> for reaction III.
- 2026-11-02: printed the space times at which beta reaches 0.90 as the first
  point above 0.90 on a 200-point grid, late by more than the printed digit: 0.345, 0.432, 0.242 and
  0.327 are 0.343 <!-- claim: tau90_II_798K -->, 0.430
  <!-- claim: tau90_I_823K -->, 0.240 <!-- claim: tau90_II_823K --> and 0.326
  <!-- claim: tau90_I_848K --> when root-found.
- 2026-11-02: printed beta at W/F = 0.45 as "the end of the measured range"; the
  last measurement is at 0.405 <!-- claim: W_F_measured_max -->, and 0.45 is the outlet of the computed bed.
