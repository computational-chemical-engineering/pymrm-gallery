# F1.8: Shah, Kelkar, Godbole and Deckwer, the kLa recommendation

The kLa statements of the bubble-column review by Shah, Kelkar, Godbole and Deckwer,
*AIChE J.* **28**(3) 353-379 (1982),
[doi:10.1002/aic.690280302](https://doi.org/10.1002/aic.690280302): Akita-Yoshida as
the conservative estimate, kLa = b uG^n with one exponent n, and that "an influence of
$u_L$ on $k_La$ can be neglected". The page tests them against the correlations, table and figures the review
itself prints. A correlation page; it calls no pymrm operator.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/shah-1982-printed-numbers.csv` | every constant, stated range and figure condition the page uses, as printed |
| `data/shah-1982-fig6a1.csv` | kLa against uG at three liquid velocities, digitised from Figure 6-a.1 |
| `data/shah-1982-fig6a2.csv` | measured kLa against the Akita-Yoshida group, digitised from Figure 6-a.2 |
| `data/check_values.csv` | the correlation functions at one test condition, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it. The water properties and
the Akita-Yoshida constants of the original papers are loaded from
[F1.1](../F1.1-akita-yoshida-bubble-column/).

## Things to know before reusing this

- Akita-Yoshida is below the nozzle row of Table 4-a.1 over the part of the row's window
  it covers, at the row's diameter: at most 0.7336 <!-- claim: noz_max_ratio --> of it,
  1.51 <!-- claim: noz_gap_over_err --> printed mean errors or more. Against Hikita it is
  lower up to D = 0.2752 <!-- claim: hik_cross_D --> m, beyond Hikita's own diameters;
  at Akita-Yoshida's largest stated diameter it would be up to
  1.1417 <!-- claim: hik_max_ratio_ay_Dmax --> times Hikita. Choosing it is a choice of
  a lower kLa by these factors.
- The local exponent of Akita-Yoshida is not the n of the fitted rows: it falls with the holdup, and its
  secant exponent over the nozzle window is 0.904 <!-- claim: secant_exp_noz_window -->.
- The Nakanoh-Yoshida velocity group is printed without its square root; the page
  restores it (cm-g-s over SI is 0.0100 <!-- claim: v2_ny_printed_cgs_ratio --> as
  printed). Water is outside that correlation's viscosity range.
- The Akita-Yoshida liquid-velocity range is printed in two tables of the review a
  factor of ten apart; the page does not decide which is right.
- Figure 6-a.1, the review's evidence that "an influence of $u_L$ on $k_La$ can be
  neglected", shows at matched gas velocity the highest liquid velocity lying
  0.222 <!-- claim: ul_pair_diff --> ± 0.033 <!-- claim: ul_pair_se --> above the lowest
  in log ratio to the drawn line (t = 6.7 <!-- claim: ul_pair_t -->,
  5 <!-- claim: ul_pair_n --> pairs). The unmatched comparison of all markers gives
  0.300 <!-- claim: ul_all_diff --> ± 0.074 <!-- claim: ul_all_diff_se --> (Welch
  t = 4.06 <!-- claim: ul_all_diff_t -->).
