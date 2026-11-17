# F3.2: Van Krevelen–Hoftijzer, the chart replaced by a solve

The enhancement factor for a second-order gas-liquid reaction is an implicit
relation that textbooks hand over as a chart. This page solves the two film
equations with pymrm, prices the Van Krevelen-Hoftijzer approximation against
them over the whole printed chart, tests the assumption it is said to rest on in
both of the wordings the book prints, and checks Levenspiel's three printed
expansions of the same result.

- **Structure:** `S3` (1-D steady BVP), data tier 6 (nothing measured).
- **Origin, cited but not consulted:** Van Krevelen, D. W. & Hoftijzer, P. J.
  (1948), *Kinetics of gas-liquid reactions. Part I. General theory*, Recueil
  des Travaux Chimiques des Pays-Bas **67**(7), 563–586,
  doi:10.1002/recl.19480670708. Not on disk. Both books on disk print only
  "67, 563 (1948)": the issue number, the end page and the DOI are a Crossref
  record lookup on that DOI, not a reading of the article.
- **Read, and the source of every equation:** Froment, De Wilde & Bischoff
  (2011), *Chemical Reactor Analysis and Design*, 3rd edn, Ch. 6, sections
  6.3.1 to 6.3.5 (book pp. 326–341), eqs. (6.3.2-11) and (6.3.5-1).
- **Read, as an independent check on the constants:** Levenspiel (1999),
  *Chemical Reaction Engineering*, 3rd edn, Ch. 23, Fig. 23.4 (book p. 530).

## Headline

The reference is a Richardson extrapolation of the same graded film on
n_x = 400 and 800 cells, because a single graded solve is grid-limited exactly
where the headline sits. Against it, over the chart's own printed domain:

| | value |
|---|---|
| worst error of Van Krevelen-Hoftijzer | -2.7955 <!-- claim: vkh_worst_rel_error_pct --> % at γ = 4.728 <!-- claim: vkh_worst_gamma -->, q = 1.522 <!-- claim: vkh_worst_q --> |
| worst over the ten printed curve labels only | -2.7294 <!-- claim: vkh_worst_on_printed_labels_pct --> % |
| margin against "valid to within 10 percent" (book p. 340) | 3.5772 <!-- claim: froment_10pct_claim_margin --> |

Both extrema are roots of a central difference of the error in ln γ and ln q,
bracketed by a scan of the whole window. That the approximation never lies above
the film solve is also a theorem: β′ = (α′ + F_A)/q ≥ 0, so β(ξ) ≥ β(0)
everywhere, and (6.3.5-1) is the fixed point of the decreasing map that freezes
β at β(0).

**The stated assumption, in both of its wordings, does not order the error.**
Book p. 335 says the curves assume "B is only weakly depleted near the
interface"; book p. 340 says "the concentration of B remains approximately
constant close to the interface". Along the worst curve:

| | γ | β(0) | rate-weighted β / β(0) | VKH error |
|---|---|---|---|---|
| B half gone (root-found) | 2.204 <!-- claim: gamma_at_half_depletion --> | 0.5 | | -1.8304 <!-- claim: vkh_error_at_half_depletion_pct --> % |
| the worst point | | 0.1647 <!-- claim: beta_interface_at_worst --> | 2.115 <!-- claim: beta_constancy_at_worst --> | |
| B 99 % gone (root-found) | 13.038 <!-- claim: gamma_at_99pct_depletion --> | 0.01 | 19.20 <!-- claim: beta_constancy_at_99pct_depletion --> | -1.5035 <!-- claim: vkh_error_at_99pct_depletion_pct --> % |
| B gone (γ = 1000, axis end) | | at most 1e-12 <!-- claim: beta_interface_at_axis_end --> | | -0.00037 <!-- claim: vkh_error_at_full_depletion_pct --> % |

On a line of 21 points from γ = 0.5 to 100, β(0) falls on
20 <!-- claim: depletion_falls_steps --> steps of 20 and the constancy ratio rises
on 20 <!-- claim: constancy_rises_steps -->, while the error peaks in between and
vanishes at both ends, because the ceiling F_A → 1 + q is a limit the
approximation reproduces by construction. Froment's own sentence is qualified to
moderately fast reactions, and measured, that qualified claim stands.

**Levenspiel's expansions.** 1 + M_H²/3 holds to 1 % up to a root-found
M_H = 0.7231 <!-- claim: lev_small_M_1pct_threshold -->. At the boundary
E_i = 5 M_H his box declares, his pseudo-first-order expansion errs by
-0.5366 <!-- claim: lev_pfo_expansion_err_pct --> % against
9.3004 <!-- claim: uncorrected_E_eq_M_err_pct --> % for E = M_H, and the same
expansion of (6.3.5-1) by -0.7374 <!-- claim: vkh_form_pfo_expansion_err_pct --> %.
His instantaneous correction E_i²(E_i−1)/M_H² is the large-γ asymptote of
(6.3.5-1) only when E_i is large too (on the worst curve the VKH deficit at the
axis end is 0.9727 <!-- claim: vkh_deficit_over_lev_correction_axis_end --> of it),
and never of the film equations: it tracks the VKH deficit to
0.9953 <!-- claim: lev_inst_asymptote_ratio_vkh --> of itself at M_H = 400 and the
film's deficit to only 0.2286 <!-- claim: lev_inst_asymptote_ratio_num -->. The
film's deficit is q β(0) exactly (the first integral), exponentially small: its
local exponent steepens to -5.3385 <!-- claim: inst_deficit_exponent_num_steepest -->
while the VKH deficit's converges on -2 (-1.9949 <!-- claim: inst_deficit_exponent_vkh_finest -->
on the finest octave), so at M_H = 800 the printed correction is
44.24 <!-- claim: lev_inst_overstatement_factor --> times the true deficit. Two
independent methods agree on that deficit to at most
1e-3 <!-- claim: inst_deficit_routes_max_rel_diff --> relative.

**The two charts' labels.** Froment labels his curves with q, Levenspiel the same
family with E_i, and both print the same ten numerals, so the same label implies
plateaus that differ by a factor up to 2.0 <!-- claim: chart_label_max_plateau_ratio -->.

## Validation

| check | result |
|---|---|
| one graded n_x = 800 solve vs. scipy collocation, 20 chart points | at most 1e-5 <!-- claim: fv_vs_bvp_max_rel_diff --> in F_A, at most 5e-6 <!-- claim: fv_vs_bvp_max_abs_diff_beta --> in β(0) |
| the reference vs. collocation | at most 1e-7 <!-- claim: ref_vs_bvp_max_rel_diff --> |
| the reference vs. a finer extrapolation (800, 1600) | at most 1e-7 <!-- claim: ref_self_refinement_max_rel_diff --> |
| first integral q β(0) = 1 + q − F_A | at most 1e-8 <!-- claim: first_integral_max_abs_resid --> |
| grid order at (γ, q) = (10, 4), successive differences, 100 to 800 cells | 2.00 <!-- claim: fv_convergence_order --> |
| finest ladder difference over the Newton tolerance | at least 1e3 <!-- claim: fv_ladder_over_newton_tol --> |
| graded n_x = 800 error at (10, 4) | 1.90e-6 <!-- claim: fv_rel_error_n800 --> |
| uniform n_x = 800 error at the same point | at most 1e-8 <!-- claim: fv_rel_error_n800_uniform --> |
| the implicit relation, for F_A and for γ′ | at most 1e-12 <!-- claim: vkh_two_routes_max_rel_diff --> |
| (6.3.5-1) at q = 1e10 against γ/tanh γ | at most 2e-8 <!-- claim: vkh_pfo_limit_max_rel_error --> |
| (6.3.5-1) at γ = 1e6 against 1 + q | at most 1e-9 <!-- claim: vkh_instantaneous_limit_rel_error --> |

**What grading costs and buys.** At the grid-study point (γ, q) = (10, 4) the graded
grid is at least 100 <!-- claim: grading_cost_ratio --> times less accurate than
a uniform one of the same size, which is why the reference extrapolates. On
`F3.1`'s own pseudo-first-order check to Ha = 200 it is the other way round: this
page's graded grid is 893.9 <!-- claim: f31_first_order_error_ratio --> times
tighter than `F3.1`'s stored number, of which grading at the same n_x is a factor
124.8 <!-- claim: grading_worth_ratio --> and the refinement from 400 to 800
cells a factor 7.163 <!-- claim: refinement_worth_ratio -->.

**Reconciliation with `F3.1`.** `F3.1`'s stored maximum VKH error is loaded, not
retyped. Its definition recomputed here on the reference gives
2.0769 <!-- claim: f31_definition_vkh_max_error_pct --> %, within at most
2e-5 <!-- claim: f31_stored_vs_recomputed_rel --> relative of the stored value,
and on `F3.1`'s own uniform 400-cell grid it reproduces the stored value to at
most 1e-9 <!-- claim: f31_cold_own_grid_rel -->, so the residual is that grid's
discretisation. Root-finding γ at the same three E_i gives
2.0785 <!-- claim: f31_definition_rootfound_pct --> %: almost all of the gap to the
headline comes from letting q move. The first-order check is reproduced on
`F3.1`'s grid to at most 1e-12 <!-- claim: f31_first_order_reproduced_rel -->;
graded n_x = 800 gives 7.088e-6 <!-- claim: pfo_max_rel_error_graded -->, uniform
n_x = 800 8.847e-4 <!-- claim: pfo_max_rel_error_uniform_same_n --> and uniform
n_x = 400 6.34e-3 <!-- claim: pfo_max_rel_error_uniform_half_n -->.

A break table injects one defect per row (a sign, offset or group of (6.3.5-1),
the depletion of B, the grid, the reference, the bulk boundary condition, the
extremum search, the label set, Levenspiel's four constants, three solver
tolerances) and every claim fails at least one row.

## What pymrm adds

The approximation itself needs nothing from pymrm: it is one scalar equation.
What pymrm supplies is the reference: the two coupled, nonlinear film equations
on a graded grid (`construct_grad` and `construct_div` take arbitrary face
positions), per-species boundary conditions, `newton` with `NumJac` on the
pointwise rate and a lower bound of zero, `compute_boundary_values` for β(0), and
`pymrm.checks` for the root scans and the refinement ladder. With it, the
accuracy claim, the stated assumption and one of Levenspiel's corrections become
measurements.

## Data

Provenance tier 6, nothing measured. The two CSVs are transcriptions: the
numeric constants, thresholds, axis ends, the two printed wordings of the
assumption, and the ten curve labels both charts share. No curve is digitised.
Froment's equations lose every Symbol-font operator under text extraction, so
everything was read on renders; his Fig. 6.3.2-1 is an embedded 148 ppi raster
and its labels were read at that native resolution. Levenspiel spells the first
author "van Krevelens" in both places he names him; Froment cites the same paper
(volume 67, page 563) as Van Krevelen. The two books' citations do not
conflict: Froment cites the 1948 paper, Levenspiel's caption a 1954 one, and his
reference list carries both.

## Corrections

- 2026-11-17: stated "Putting F_A = 1 + q − d into (6.3.5-1) gives d → q(1+q)²/γ² for large γ, which in Levenspiel's symbols is E_i²(E_i−1)/M_H²: his correction is the algebraic asymptote of the approximation." The large-γ limit at fixed q is d → q x²/γ² with x = (1 + q) tanh x, which is Levenspiel's correction only as E_i grows; on the worst curve the ratio is 0.9727 <!-- claim: vkh_deficit_over_lev_correction_axis_end -->.

## What this page cannot conclude

Nothing about absorbers. Both sides of every comparison are equations, and the
1948 origin was not read, so nothing here is evidence about what Van Krevelen and
Hoftijzer wrote, only about the equation two monographs attribute to them.

## Files

| file | contents |
|---|---|
| `index.py`, `index.ipynb` | the page (paired by jupytext) |
| `model.py` | the film by finite volumes, its Richardson extrapolation, the collocation route, the implicit relation by two routes, the root helpers |
| `test_model.py` | fast offline checks of `model.py` |
| `data/vkh-printed-claims.csv` | transcribed constants and sentences of both books |
| `data/enhancement-chart-labels.csv` | the ten printed curve labels |
| `agreement.json` | every claim, written by the notebook |

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```
