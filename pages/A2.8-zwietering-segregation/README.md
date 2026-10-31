# A2.8: Zwietering's micromixing bounds

**Catalog ID:** `A2.8` · **Section:** A · **Tier:** T0 · **Priority:** P2 ·
**Structures:** `S1`, `S3`, `S5` · **Data tier:** 6 (no measurement exists)

Reproduces Th. N. Zwietering, *The degree of mixing in continuous flow systems*,
**Chemical Engineering Science 11(1) 1-15 (1959)**,
[doi:10.1016/0009-2509(59)80068-3](https://doi.org/10.1016/0009-2509(59)80068-3).

> **The page range is not 1-11.** The article is cited that way essentially
> everywhere. Its last page carries the printed folio **15**, holds eqs. (II, 11)
> and (II, 12) and the three-item reference list, and the printed folios run 1 to
> 15 without a gap. Read off a 300 dpi render, not from metadata.

## What the page shows

For a reaction of any order other than one, knowing the residence-time
distribution is **not enough** to fix the conversion. Zwietering's two extremes,
Danckwerts' complete segregation (mixing as late as possible) and his own state of
maximum mixedness (as early as possible), bracket every reactor with that RTD.

The page reproduces all 44 <!-- claim: n_printed_conversions --> conversions and
all 8 <!-- claim: n_printed_J --> degrees of segregation he printed, with **zero
fitted parameters**, and then measures the thing the bounds exist for: how wide the
bracket is. It is wide, and it widens with reaction rate: a factor
**2.12** <!-- claim: bracket_max_ratio_hi_over_lo --> between the two bounds at
$kc_0\tau = 50$ for a two-tank RTD. Both of Zwietering's *intermediate* reactors
sit well inside the band, between 0.23 <!-- claim: bracket_mu_min --> and
0.85 <!-- claim: bracket_mu_max --> of its width, so the bracketing check is a check
that can fail.

## Two findings about the paper

1. **The conversion bodies of Tables 1 and 2 are interchanged with respect to
   their captions.** Proved from three closed forms Zwietering printed himself:
   eq. (42), derived from eq. (41)'s $2/\tau$ and therefore unambiguously the
   two-vessel chain, reproduces Table 2's row (captioned *three* vessels) to
   4.1e-4 <!-- claim: eq42_vs_table2_max_dev --> and Table 1's only to
   3.3e-2 <!-- claim: eq42_vs_table1_max_dev -->; eqs. (38) and (39)+(40) do the
   same for the other two algebraic rows. Scored over all 44 values the
   interchanged reading gives a worst deviation of
   4.5e-3 <!-- claim: tables_max_dev_interchanged --> against
   3.5e-2 <!-- claim: tables_max_dev_as_captioned -->. The $J$ column does **not**
   move with the body, which is what makes it an interchange rather than a
   relabelling.
2. **One printed value is not reproduced, and is not repaired.** The three-tank
   degree of segregation under maximum mixedness computes to
   0.080296 <!-- claim: J_mm_3tank_computed --> against the printed 0.0831. The
   two-tank companion meets its closed form and Appendix I's variance identity
   closes. Reported as unresolved.

Also measured: Zwietering's own hand-integration error. Every *algebraic* row in
both tables is inside one unit of his last decimal (0.001), though
4 <!-- claim: tables_algebraic_cells_outside_half_unit --> algebraic cells miss by
more than half a unit; the only row outside 0.001 is the
maximum-mixedness row he had to integrate by hand for the two-tank RTD, which
carries a **+1.10 %** <!-- claim: mm_2tank_mean_bias --> bias with all four printed
deviations of one sign.

## What pymrm contributes

Eq. (31) is normally solved as an initial-value problem integrated backwards from
infinite life expectancy. Written conservatively,

```
div(v gamma) = E(x) - K w(x) gamma^p ,    v = -w(x) ,    x = lambda/tau
```

it becomes an ordinary convection-reaction boundary-value problem with a decaying
velocity field and a distributed side feed: `construct_convflux_upwind` and
`construct_div`, van Leer face fluxes with their exact Jacobian from `ConvFluxTVD`,
the source derivative as a diagonal from `construct_coefficient_matrix`, and
`newton(damping="monotone", stagnation=True)`. The life-expectation coordinate
behaves like any other spatial axis: second order on a graded grid, the side feed
conserved to round-off with the exit read on the face (`compute_boundary_values`),
and a handful of Newton iterations. Van Leer deferred correction is kept as a
second route to the same discrete solution. For this scalar equation the grid is
not needed: the gridless DOP853 integration at a loose tolerance is as fast and
more accurate than the 800-cell solve. What the conservative form adds is the
closed discrete balance, the whole profile as one Newton system, and the operator
stack a model with side feeds, a second phase or a measured RTD on a grid reuses.

The dashes in Zwietering's maximum-mixedness rows are filled in. He integrated
eq. (31) only for $K = 5, 10, 20, 30$.

## Validation summary

| Check | Result |
|---|---|
| printed conversions, nothing fitted | worst 4.5e-3 <!-- claim: tables_max_dev_interchanged -->, rms 8.3e-4 <!-- claim: tables_rms_dev_interchanged -->, 40 <!-- claim: tables_within_0p001_interchanged --> of 44 inside 0.001 |
| the same, read as captioned | worst 3.5e-2 <!-- claim: tables_max_dev_as_captioned -->, 0 <!-- claim: tables_within_0p001_as_captioned --> of 44 inside 0.001 |
| pymrm finite volume vs gridless DOP853 integration | 8.6e-7 <!-- claim: mm_pymrm_vs_dop853_max --> over 12 cases (the 800-cell grid error); Richardson limit below 1e-8 <!-- claim: mm_richardson_vs_dop853_max --> |
| deferred correction vs Newton, same discrete equations | below 1e-11 <!-- claim: mm_deferred_vs_newton_max --> in every cell |
| discrete balance with the exit read on the face | below 1e-12 <!-- claim: mm_balance_gap_face --> |
| grid refinement | order 2.00 <!-- claim: mm_grid_order_vanleer --> (van Leer), 1.00 <!-- claim: mm_grid_order_upwind --> (bare upwind), every rung of one sign |
| domain truncation | decay 8.44 <!-- claim: mm_truncation_decay_rate --> per unit X at K = 10; below 1e-13 <!-- claim: mm_truncation_error_at_X6_all_cases --> at X = 6 in every printed case (the DOP853 tolerance floor) |
| exit at the first cell centre, uniform grids | first order (0.90 to 0.99 from 200 to 3200 cells), face read second order |
| exponential RTD to stirred-tank root (his eq. II, 12) | below 1e-12 <!-- claim: limit_exponential_rtd_max_dev --> |
| first order: both bounds collapse onto (1+K/n)^-n | below 1e-12 <!-- claim: limit_first_order_vs_exact --> (the routes share the RTD functions) |
| plug-flow limit (300 tanks) | bracket 4.0 % <!-- claim: limit_plugflow_rel_width_n300 --> wide, down from 83 % <!-- claim: bracket_rel_width_K30_n2 --> at 2 tanks |
| bracket reverses through order one | +8.3e-3 <!-- claim: order_margin_p_1p1 --> at p = 1.1, -7.4e-3 <!-- claim: order_margin_p_0p9 --> at p = 0.9 |
| printed degrees of segregation | 7 <!-- claim: J_reproduced_within_1e-4 --> to 4.3e-5 <!-- claim: J_max_dev_seven_of_eight --> (four exact rationals); 1 unresolved |
| Appendix I variance identity | below 1e-12 <!-- claim: J_variance_identity_residual --> |
| defect injection | every claim failed by at least one row |

**Labelled untested:** the far-field boundary condition. Replacing the root by
`gamma = 1` moves the answer by below 1e-12 <!-- claim: break_farfield_root_shift -->
at X = 6. Shorten the domain to X = 0.5 and the same defect is worth
1.7e-2 <!-- claim: break_farfield_root_shift_X0p5 -->, which is what shows the
insensitivity belongs to the domain length rather than to the check.

## Corrections

- 2026-10-31: the break table said the flipped reaction sign "sends the answer to
  28984", with the interior reactor still inside the bracket, and that the solver
  refusing to converge catches it. 28984 was the iterate at which the
  deferred-correction loop was stopped unconverged (the same table printed
  "solver converged: False"). With the sign flipped there is no steady solution.
  Eq. (31) is first order with one inflow condition, so the backward integration
  from x = 6 is the only candidate, and with the flip it blows up at
  x = 5.8 <!-- claim: flip_blowup_x_dop853 --> (DOP853; the Riccati equation with
  E/w frozen gives the same place), and the flipped far-field equation has no real
  root (discriminant -76 <!-- claim: flip_farfield_discriminant -->). Newton with
  the exact Jacobian does not converge either.
- 2026-10-31: stated that replacing the far-field root by `gamma = 1` moves the
  answer by 5.9e-7, and by 2.1e-2 at X = 0.5. Both were distances from the
  X = 6 reference and carried the grid error (5.9e-7) and the X = 0.5 truncation
  error (4.1e-3). Against the same solve with the root, the defect moves the answer
  by round-off at X = 6 and by 1.7e-2 <!-- claim: break_farfield_root_shift_X0p5 -->
  at X = 0.5.
- 2026-10-31: stated that reading the exit at the first cell centre costs 5.0e-4 on
  a uniform 800-cell grid (585 times the reported precision) and 9.7e-7 on the
  graded grid. Both were distances from the reference and included each grid's
  own face error. The cell read itself moves the exit by
  4.7e-4 <!-- claim: break_exit_cell_uniform_grid --> (uniform) and
  3.9e-7 <!-- claim: exit_cell_shift_graded --> (graded, below the grid error).
- 2026-10-31: stated the domain-truncation error at the production X = 6 as
  4.1e-15, a factor 2e8 below the grid error. 4.1e-15 is the floor set by the
  DOP853 tolerance, reached at X = 4 for K = 10; over all printed cases the
  difference is 3.5e-10 <!-- claim: mm_truncation_error_at_X4_all_cases --> at X = 4
  (K = 3, two tanks) and below 1e-13 at X = 6.
- 2026-10-31: said the segregation quadrature and the backward integration of the
  first-order check "share no code". They share the RTD functions.
- 2026-10-31: the table of filled dashes divided the computed maximum-mixedness
  value by the *printed* three-decimal segregated value, giving 1.151 and 1.044 at
  K = 50; on the same (computed) basis the widths are
  1.119 <!-- claim: filled_dash_width_K50_2tank --> and
  1.010 <!-- claim: filled_dash_width_K50_3tank -->. The prose compared 0.0753 with
  the printed 0.035; the computed segregated value is 0.0355, whose three-decimal
  rounding is 0.036, not the printed 0.035.
- 2026-10-31: said the break table measures what Appendix I's variance identity
  catches; no row touched it. A within-point age distribution taken from the
  three-tank RTD now breaks it.
- 2026-10-31: said every algebraic row is inside Zwietering's three-decimal
  rounding and that the two-tank maximum-mixedness row is the only place outside
  it. Every algebraic row is inside one unit of his last decimal (0.001), but
  4 <!-- claim: tables_algebraic_cells_outside_half_unit --> algebraic cells miss
  by more than half a unit (two-tank segregated chain at K = 3 and 10, two-tank
  segregation at K = 50, three-tank segregated chain at K = 3); of the cells that do
  not round to the printed digits, half are in the maximum-mixedness row.
- 2026-10-31: said `A2.4` (tanks in series) is not built; it was published on
  2026-08-07.

## Files

```
index.py           the page source (jupytext percent format), paired with index.ipynb
index.ipynb        the executed notebook
model.py           RTD functions, the four mixing states, MaxMixedness, mm_ode, J
test_model.py      fast checks of model.py and one stored claim
meta.yaml          page metadata
agreement.json     the claims, format 2, written by the last cell
data/zwietering-1959-conversions.csv          printed values of ce/c0
data/zwietering-1959-segregation-degree.csv   printed values of J
```

Runtime about 40 s.

## Related pages

[`A2.4`](../A2.4-tanks-in-series/) (tanks in series),
[`A2.1`](../A2.1-danckwerts-boundary-conditions/) (Danckwerts boundary
conditions), [`A2.3`](../A2.3-taylor-aris-dispersion/) (Taylor-Aris),
[`A2.5`](../A2.5-edwards-richardson-dispersion/) (measured dispersion),
[`A2.6`](../A2.6-gunn-dispersion-correlations/) (dispersion correlations). These
end with an RTD; this page is about what an RTD still leaves undetermined. No
sibling page's dataset is loaded.
