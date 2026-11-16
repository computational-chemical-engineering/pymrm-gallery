# C1.2: Eley-Rideal kinetics: the impact rate law, and the case that the name is wrong

**Catalog ID:** `C1.2` · **Section:** C · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1` · **Data tier:** 6 (analytical; the source prints no
data, see below)

Built from Prins, R., "Eley-Rideal, the Other Mechanism", *Topics in
Catalysis* **61** (2018) 714-721, doi:10.1007/s11244-018-0948-8, a
historical and terminological review whose thesis is that the mechanism
universally called **Eley-Rideal** (a gas molecule colliding with a
chemisorbed one) was proposed by **Langmuir in 1922** and should be called
**Langmuir-Rideal**, while what Eley and Rideal actually studied (1940-41)
is a chemisorbed atom reacting with a *physisorbed* molecule. The page
reports that finding plainly (same class as `B1.4` and `H1.1`), builds the
impact rate law, and contrasts it with the Langmuir-Hinshelwood law
**imported from the published page `C1.1`** rather than re-derived.

- **pymrm:** `pymrm.checks.find_roots` locates every extremum and threshold on
  its physical range (exactly one root required), and pymrm's `newton` on a
  `NumJac((1, 1))` Jacobian is the independent cross-check for the two maxima.
  The page says plainly that pymrm contributes nothing else: the rate laws, the
  fits and the comparisons are NumPy, SciPy and sympy.

## What the page shows

- **Every printed limit of Eqs. (3a)-(3c), derived.** All 11 <!-- claim: printed_limits_agree --> limiting-order
  claims of p. 719 hold; derived symbolically and measured numerically on
  independently coded laws.
- **One printed defect, reported and not repaired.** "The rate of the L-H
  reaction should even reach a maximum for c_A = K_A." [sic] is dimensionally
  inconsistent with the paper's own Eq. (3a); the derived condition
  K_A·c_A = 1 + K_B·c_B is unit-invariant where the printed one moves by λ²
  under a unit rescale. The maximum is the single root of the numeric law on
  its physical range, and agrees with the closed form to at most 1e-8 <!-- claim: maxloc_closed_minus_root_rel -->,
  and with `newton` to at most 1e-6 <!-- claim: maxloc_newton_minus_root_rel -->.
- **The p. 719 paragraph turned into a number.** The identifiability window:
  the chemisorbing species must span a factor w* = 11.8 <!-- claim: wstar_tol844 --> (root-found,
  two optimizer routes to at most 1e-6 <!-- claim: wstar_tol844_two_routes_rel -->) centred on the L-H maximum
  before the best impact-law fit misses noiseless L-H rates by more than the
  1947 fit's own 8.44 <!-- claim: codimer_tableC_avg_abs_pct --> % average deviation, and centring the window where
  both laws are first order collapses the misfit at w* to 0.000053 <!-- claim: dev_at_wstar_center001 -->:
  width buys nothing without curvature.
- **The impact law on the founding LHHW dataset (numerical experiment).**
  Hougen & Watson's own mechanism (o) *is* Prins's Eq. (3b) (with the
  product's adsorption term retained). Refit to noiseless mechanism-(d)
  rates over the 1947 design it misses by 27.4 <!-- claim: lr_o_on_exact_d_mean_pct --> % mean (48 <!-- claim: lr_o_on_exact_d_max_pct --> % worst) at the
  least-squares-optimal fit, the declared protocol; directly minimising the mean
  |Δ| instead gives 22.8 <!-- claim: lr_o_direct_mean_pct --> %, conclusion unchanged. The noise floor is 16.8 <!-- claim: codimer_d_scatter_pct --> %. Both exceed it, so on noiseless rates the design separates the two laws. But allow a hydrogen term in the
  first-power denominator (3.3 <!-- claim: lr_oKH_on_exact_d_mean_pct --> %, nonzero) or collapse the design's decade-wide
  hydrogen span (3.2 <!-- claim: lr_o_collapsedpH_mean_pct --> %) and it drops below the noise floor, and the
  noiseless misfit falls short of the scatter: the design separates them only because it spans a decade in p_H (the design also
  straddles the L-H turnover, p_U* = 2.50 <!-- claim: codimer_pUstar_atm --> atm, inside the design range: an
  observation, not a second necessary condition, since collapsing the hydrogen
  span alone already suffices). The denominator exponent (1 vs 2) is
  indistinguishable there at the design's noise floor even when the rates carry
  no noise. Against the measured rates the d-versus-o gap is 9.8 <!-- claim: lr_gap_measured_pct_points --> points with a
  standard error of 6.3 <!-- claim: lr_gap_measured_se_pts --> over the 12 <!-- claim: runs_at_200C --> runs, so it is not resolved at two standard errors.
- **Cross-page reconciliation that caught what it should.** The imported
  constants reproduce the book's printed Table C to within one unit of the last
  printed digit on 11 <!-- claim: codimer_Rcalc_rows_within_one_unit --> of the 12 <!-- claim: runs_at_200C --> runs; run 25a is off by 0.503 <!-- claim: codimer_25a_abs_gap -->,
  independently reproducing `C1.1`'s documented digit slip.

## No data, and what that means

The source prints **no measurement, no table and no figure** (checked by
full-text read, page-by-page render inspection, and the image inventory). The
page therefore claims **no experimental validation**: everything is analytical
or a labelled numerical experiment on `C1.1`'s fitted constants, which are
themselves fit data, as `C1.1`'s own metadata records.

## Corrections

- 2026-11-10: stated that "fit quality alone could have rejected it" (Eq. (3b) on the 1947 design). Fitted to the 12 <!-- claim: runs_at_200C --> measured runs at 200 °C, the impact law's mean residual exceeds mechanism (d)'s by 9.8 <!-- claim: lr_gap_measured_pct_points --> points with a paired standard error of 6.3 <!-- claim: lr_gap_measured_se_pts -->, so fit quality alone does not separate the two laws at two standard errors on these runs.

## Files

- `index.py` and `index.ipynb`: the page, paired with jupytext (edit `index.py`,
  then `jupytext --sync`).
- `data/prins-2018-printed-statements.csv` + `.meta.yaml`: every quantitative
  statement the paper prints, with page provenance (all numerics verified on
  render crops).
- Three `C1.1` datasets are loaded cross-page; that page's findings on the exact
  rows used are listed in the notebook's data section.
- `agreement.json`: every number the page states, with how it is checked and the
  break-table row that moves it.

Runtime: about 30 s end to end.
