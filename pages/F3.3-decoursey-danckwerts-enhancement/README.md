# F3.3: DeCoursey (1974), an explicit enhancement factor for the Danckwerts model

W. J. DeCoursey, "Absorption with chemical reaction: development of a new relation for the
Danckwerts model", *Chem. Eng. Sci.* **29**, 1867-1872 (1974),
doi:10.1016/0009-2509(74)85003-7. DeCoursey derived an explicit enhancement factor, Eq. (21),
for a second-order reaction under Danckwerts surface renewal, but "No numerical solution ...
for the Danckwerts model appears to be available" (p. 1870), so he compared it only with
Higbie-model numbers and predicted it "should be within about 12 per cent of the exact results
for the Danckwerts model at all practical values of the parameters" (p. 1870). This page solves
the exact Danckwerts model over the paper's own grid and decides that prediction.

- **Structure:** `S4` (1-D transient PDE; the catalogue line says `S3`), tier T1, data tier 6
  (nothing measured).

## Findings

- **The prediction holds, with room to spare.** Over the 12 cases of Table 2 the worst error
  of Eq. (21) against the exact Danckwerts model is 6.2705 <!-- claim: r1_worst_abs_pct --> %
  in magnitude (at E_a = 11, D_B/D_A = `0.02`, M = 46.78 <!-- claim: r1_worst_M -->), so
  12 <!-- claim: r1_n_within --> of 12 cases are within the paper's "about 12 per cent", a margin of
  1.9137 <!-- claim: r1_margin -->. For D_B/D_A of 1 and above the worst error is
  3.2053 <!-- claim: r1_max_r_ge1_pct --> %. Every worst point is root-found inside the
  window and Richardson-extrapolated (largest error estimate
  0.0039 <!-- claim: r1_max_errest --> points).
- **The film comparison of F3.1 is mostly a model difference.** At E_a = 100, D_B/D_A = 1 the
  worst error against the Danckwerts model is 2.8551 <!-- claim: r2_ea100_worst_pct --> %,
  0.32704 <!-- claim: r2_share_of_f31 --> of F3.1's film figure (loaded from F3.1, not
  retyped); at F3.1's worst point the Danckwerts error is
  -0.0011552 <!-- claim: r2_at_f31_M_pct --> %.
- **Compensation.** At equal diffusivities the bridging relation (15) is exact for the
  Danckwerts model (extrapolated gap at most 1e-6 <!-- claim: r2_idgap_ext -->), so the error is
  the "B_i constant in time" step. That step alone (Eq. 19 fed the exact mean) overestimates
  E by up to 29.988 <!-- claim: r3_e19_max_pct --> % on the window, still rising at its end,
  while Eq. (21) stays within 2.4208 <!-- claim: r3_e21_max_r1_pct --> %: the self-consistent
  mean compensates. A measured gap on this grid, not a general result.
- **Table 1 recomputed** with the exact Higbie average: 2 <!-- claim: r4_n_round --> of 7
  printed corrections round from the recomputed value, 7 <!-- claim: r4_n_allow --> of 7 are
  within the accuracy Brian et al. state for their charts plus half a unit; the largest gap is
  1.6888 <!-- claim: r4_max_gap --> points, at M = 24 (printed 6, recomputed
  4.3112 <!-- claim: r4_corr_4 -->). **Table 2**: 24 <!-- claim: r5_n_inside --> of 24
  printed range ends are within that allowance, and the exact maxima are
  4.8421 <!-- claim: r5_max_equal_pct --> % (equal diffusivities) and
  12.907 <!-- claim: r5_max_002_pct --> % (D_B/D_A = `0.02`), under the paper's "no more than 6"
  and "no more than 14".
- **The premise** that the two models "agree within 3 per cent" (p. 1870) holds: the Higbie
  pseudo-first-order form is below the Danckwerts one by -2.5494 <!-- claim: r6_pfo_worst_pct --> %
  at worst (M =
  1.2659 <!-- claim: r6_pfo_worst_M -->).
- **Misprint:** Eq. (23) prints erfc where erf is meant (as printed it tends to infinity as
  M goes to 0).

## Things to know before reusing this

- Eq. (21) is a Danckwerts-model estimate; it is not a film-model (F3.1) or Higbie estimate.
- At unequal diffusivities q = B0/(z A_i) comes from the exact penetration relation for E_a,
  not from Eq. (14), which the paper says is exact only for equal diffusivities.
- The fast-reaction approach falls as 1/M: d ln(E_a - E_D)/d ln M is
  -0.99164 <!-- claim: v1_fast_slope --> between M = `1e4` and `1e5`.
- Not covered: reversible reactions, A_0 not zero, other orders, experiments.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | Eqs. (20) to (23), the pseudo-first-order forms, the penetration relation for q, `RenewalMarch` (the exact model, pymrm finite volumes and newton per BDF2 step) and `History` (Danckwerts and Higbie averages and their derivatives) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/decoursey-1974-printed-numbers.csv` | Tables 1 and 2 and the stated accuracies, as printed strings, with their provenance sidecar |
| `data/check_values.csv` | values of the page's functions at the test condition from a second route, asserted before any claim |
