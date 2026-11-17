# C1.1: Langmuir-Hinshelwood-Hougen-Watson kinetics

**Catalog ID:** `C1.1` · **Section:** C · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1`, `S2` · **Data tier:** 2 (printed experimental tables;
everything is *fit* data, see below)

Reproduces the founding worked example of the LHHW formalism: O. A. Hougen &
K. M. Watson, *Chemical Process Principles, Part Three: Kinetics and
Catalysis*, John Wiley & Sons (1947), Chapter XIX, Illustration 2 (book
pp. 943-958): the vapor-phase hydrogenation of codimer over supported nickel,
with the complete 18-mechanism parameter estimation on the printed data of
Tschernitz, Bornstein, Beckmann & Hougen (*Trans. AIChE* **42**, 883, 1946,
not consulted; the book reprints the data in full).

The source is a book with no DOI, held as the combined volume (PDF page = book
page + 16), identified from its own title and imprint pages. Its scan is the
lowest native resolution in the gallery, so every numeral was read from
cropped native-resolution enlargements and the text layer was never trusted
for a digit.

## What the page shows

- **The transcription is proven, not asserted.** The book prints its own
  least-squares sums for the 200 °C block; the transcribed columns reproduce
  the three linear sums to at most 1e-9 <!-- claim: tableA_linear_psum_max_abs --> (digit for digit) and the product sums
  within at most 1 <!-- claim: tableA_product_psum_max_half_units --> half-unit of the printed last digit.
- **The book's estimation reproduces.** Mechanism d's constants at 200 °C by
  two independent routes (transcribed data, and the printed normal equations)
  agree to 0.0063 <!-- claim: mechd_200C_two_routes_max_rel --> relative; the book's average deviation recomputes to
  8.42 <!-- claim: tableC_avg_abs_pct --> % against its printed "8.44".
- **A 1947 worksheet slip, diagnosed.** The c and g/m families of Table D
  descend from one corrupted Σ(R·pS) entry: two independent printed rows imply
  the same wrong value (685 <!-- claim: DEF_implied_sum_RpS_c_row --> against the data's 278 <!-- claim: DEF_data_sum_RpS -->) to 0.011 <!-- claim: DEF_corrupted_sum_two_witness_gap_rel -->
  relative. The joint acceptance verdict survives it.
- **A misprint in the final equation, adjudicated.** Eq. (q) prints
  ΔS_S = -30.96 where Table G prints -30.46; the book's own K_S row decides
  for Table G (0.025 <!-- claim: KS_row_rms_rel_tableG_dSS --> rms against 0.236 <!-- claim: KS_row_rms_rel_eqq_dSS -->). Typing the boxed equation
  costs 22 <!-- claim: eqq_KS_low_pct --> % of K_S.
- **The discrimination question, answered with standard errors.** The joint
  sign-test decision {d, h} holds at every threshold swept, but fit quality
  separates little: d and h are 0.36 <!-- claim: discrim_d_vs_h_gap_pct_points --> points apart with a paired standard
  error of 0.58 <!-- claim: discrim_d_vs_h_gap_se_pts -->, the best rejected mechanism trails by 2.7 <!-- claim: discrim_third_best_gap_pct_points --> points with a
  standard error of 2.2 <!-- claim: discrim_third_best_gap_se_pts -->, and 3 <!-- claim: nls_unresolved_from_best --> mechanisms lie within two standard errors
  of the best fit. Even noiseless h-rates are fitted by mechanism d to
  3.6 <!-- claim: discrim_power_d_fits_exact_h_pct --> %. The book's fit-quality ground for choosing d does not survive;
  its chemical argument does.
- **The fitted law in a pymrm bed.** Isothermal plug flow at 200 °C and the top pressure of the data
  (a demonstration: the book poses the problem but prints no answer):
  W/F = 106.62 <!-- claim: reactor_WF_extrapolated --> for the conversion of Chapter XXI's Problem 1, extrapolated from a grid ladder with
  observed first order; adaptive quadrature agrees to at most 1e-4 <!-- claim: reactor_WF_richardson_vs_quad_rel -->. Ignoring
  the mole contraction (a fraction 0.42 <!-- claim: reactor_total_flow_contraction --> of the inlet flow) lengthens the bed by a
  relative 0.51 <!-- claim: reactor_naive_fixed_total_bias_rel -->.

## What pymrm adds

`pymrm.checks.find_roots` for the rate maximum of eq. (q) (a root of the
analytic derivative of ln r) and for the bed length that reaches the target
conversion; the plug-flow bed of Section 5 (upwind convection assembled once,
`NumJac((n, 3))`, `newton` with its success checked, an outflow outlet read by
`compute_boundary_values`) and `observed_orders` for its grid ladder. The
estimation itself (Sections 1 to 3) is NumPy and SciPy, and the page says so.

## Fit vs test

Everything in the source is fitting data: every printed constant was fitted
by the book to the same 40 runs, and no held-out measurement exists. All
agreements on the page are labelled reproduction or internal-consistency
adjudication, never validation. The one synthetic computation (noiseless
h-rates refitted by d, to measure the design's discriminating power) is
labelled as a numerical experiment in the cell that runs it.

## Corrections

- 2026-11-16: stated "So the printed data do discriminate the family shape - dual-site surface control with all three species adsorbed - while genuinely failing to resolve molecular vs atomic hydrogen." The best rejected mechanism, desorption-controlled, trails the best fit by 2.7 <!-- claim: discrim_third_best_gap_pct_points --> points with a paired standard error of 2.2 <!-- claim: discrim_third_best_gap_se_pts -->, and 3 <!-- claim: nls_unresolved_from_best --> mechanisms lie within two standard errors of the best fit. Fit quality resolves the product term (F = 8.4 <!-- claim: nls_noKS_F_200C -->, p = 0.02 <!-- claim: nls_noKS_F_pvalue_200C --> at 200 °C) and the forms 1.5 <!-- claim: nls_resolved_ratio_min --> to 3.4 <!-- claim: nls_resolved_ratio_max --> times worse than d, not the control step, which the sign test selects.
- 2026-11-16: stated "In each case the recomputation sides with the siblings." By the implied normal-equation sums, only the 275 °C c column is decided, for the siblings' 4.45 (14 <!-- claim: DEF_dup_c275_gap_over_halfunit --> times the shift of half a printed unit); the 325 °C b/q a and c-family a fall below that shift (0.95 <!-- claim: DEF_dup_bq325_gap_over_halfunit -->, 0.42 <!-- claim: DEF_dup_c325_gap_over_halfunit -->).

## Files

- `index.py` and `index.ipynb`: the page, paired with jupytext (edit `index.py`,
  then `jupytext --sync`).
- `model.py`: the mechanism library, the rate-space fit, eq. (q) and its rate
  maximum, the plug-flow bed and the quadrature; `test_model.py` checks it.
- `data/hougen-watson-1947-*.csv` + `.meta.yaml`: the six transcribed
  tables with provenance, including both printings of every constant the book
  prints twice and the printed intermediates of the estimation.
- `agreement.json`: the claims, written by the last cell, checked by CI.

Runtime: see `meta.yaml` (`runtime_seconds`).
