# F3.1: Hatta number and the gas-liquid reaction regimes

Whether a dissolving gas reacts in the film or in the bulk decides which
absorber you should build. One number settles it.

- **Structure:** `S3` (1D steady BVP), data tier 6 (nothing measured)
- **References:** Hatta (1932); Van Krevelen & Hoftijzer (1948)
  doi:10.1002/recl.19480670708; DeCoursey (1974) doi:10.1016/0009-2509(74)85003-7

## Results

The stagnant-film equations are solved with pymrm on a uniform 400-cell film.
Pseudo-first-order matches the exact `Ha/tanh(Ha)` to
6.3e-3 <!-- claim: first_order_max_rel_error --> (largest relative error over
the sweep), concentrated at high Ha: at the top of the sweep the reaction layer
spans 2 <!-- claim: reaction_layer_cells_at_ha_max --> cells. A grid study at Ha = 5 gives an observed order between 2.9 and 3.1 <!-- claim: convergence_order --> on 100 to 800 cells, a pre-asymptotic range: at Ha = 2 the order is at least 1.9 <!-- claim: convergence_order_ha2 --> and at most 2.1 <!-- claim: convergence_order_ha2 -->, the second order of the scheme. The two-component solve reduces to the first-order
limit to 5.4e-6 <!-- claim: excess_reduces_to_first_order --> when the
co-reactant is in large excess, and approaches the instantaneous ceiling
E -> E_i at large Ha.

Accuracy of the standard approximations against the film solution, at
E_i = 5, 20 and 100:

| approximation | largest error on the sweep | worst point, refined | where |
|---|---|---|---|
| Van Krevelen-Hoftijzer | 2.1 <!-- claim: vkh_max_error_pct --> % | -2.08 <!-- claim: vkh_worst_error_pct --> % at Ha = 8.23 <!-- claim: vkh_worst_ha -->, E_i = 5 | Ha = 1.65 <!-- claim: vkh_worst_ha_over_ei --> E_i: fast to instantaneous |
| DeCoursey | 8.7 <!-- claim: decoursey_max_error_pct --> % | 8.73 <!-- claim: decoursey_worst_error_pct --> % at Ha = 1.43 <!-- claim: decoursey_worst_ha -->, E_i = 100 | Ha = 0.0143 <!-- claim: decoursey_worst_ha_over_ei --> E_i: slow to fast |

The worst points are refined from a scan of the whole Ha window, and refining
the film fourfold moves neither by more than its claim bound.

**Why DeCoursey is worst where it is.** DeCoursey's relation was derived for
Danckwerts' surface-renewal model. With the co-reactant in excess it becomes
Danckwerts' `sqrt(1 + Ha^2)`, which differs from the film's `Ha/tanh(Ha)` by up
to 8.81 <!-- claim: decoursey_pfo_limit_max_error_pct --> %. Its worst
second-order error is 0.991 <!-- claim: decoursey_worst_share_of_pfo_limit -->
of that: nearly all of it is the model difference, not a second-order effect.

**Reconciliation with `F3.2`.** `F3.2` loads this page's sampled VKH maximum
and root-finds the same three E_i on its own graded, extrapolated reference;
this page loads that value and agrees with it within the bound of
`vkh_worst_vs_f32_rootfound_rel`.

## Data

**Provenance tier 6, not experimental.** The sources are analytical: Hatta's
regimes and the first-order enhancement are exact, and VKH and DeCoursey are
published approximations. Validation is against exact solutions and limiting
cases.

## Corrections

- 2026-11-19. The page said "DeCoursey is worst in the transition between the
  fast and instantaneous regimes". It is worst at
  Ha = 1.43 <!-- claim: decoursey_worst_ha --> (Ha = 0.0143 <!-- claim: decoursey_worst_ha_over_ei --> E_i),
  the transition from the slow to the fast regime, and the error there is
  0.991 <!-- claim: decoursey_worst_share_of_pfo_limit --> of the difference
  that remains with the co-reactant in excess, the film against the
  surface-renewal model.
- 2026-11-19. The page said "where the reaction layer is thinner than a cell".
  At the top of the pseudo-first-order sweep the layer spans
  2 <!-- claim: reaction_layer_cells_at_ha_max --> cells.
- 2026-11-19. The page said "reproducing it to $10^{-4}$". The largest relative error over the sweep is 6.3e-3 <!-- claim: first_order_max_rel_error -->.
