# A2.6: Gunn's dispersion correlations, and what the number 2 depends on

`A2.1` takes the axial-dispersion reactor as given. `A2.3` derives a dispersion
coefficient for an empty tube. `A2.5` measures one in a packed bed and fits two
constants to the measurement. This page is the fourth question: a stochastic model
with **one** empirical function in it, covering axial and radial dispersion,
spheres and cylinders, gases and liquids, and what its constants say about the
Péclet number 2 that everyone quotes.

- **Structures:** `S6` (2-D axisymmetric dispersion-reaction)
- **Reference:** Gunn, D. J. (1987), Chem. Eng. Sci. 42(2) 363-373,
  doi:10.1016/0009-2509(87)85066-2
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (the correlations, each taking the dictionary of printed constants,
  and `Bed`, the 2-D bed on the pymrm operators), tested by `test_model.py`. The
  closed forms, the collocation route and the studies stay in the notebook.
- **Runtime:** about 200 s.

## The headline

**Pe → 2 is a property of spheres, not of packed beds.** The high-Reynolds-number
limit of Gunn's eq. (1) is `Ud/D → 2p∞/(1−p∞)`, and eq. (3) gives
`p∞ = 0.17 + 0.33 = 0.50` exactly, so the limit is exactly 2. Eqs (4) and (5),
fitted to cylinder measurements that never see a sphere, give
1.7037 <!-- claim: pe_limit_solid_cylinders --> and
1.1746 <!-- claim: pe_limit_hollow_cylinders -->. The habit of writing
`D = u dp / 2` at high flow is a sphere habit.

Two printed digits fix the sphere limit only to about
4.1 % <!-- claim: pe_limit_rounding_width -->. The spread between shapes,
a factor of 1.703 <!-- claim: pe_limit_spread_spheres_over_hollow --> from spheres to
hollow cylinders, is outside rounding: with every constant pushed to the opposite
end of its rounding interval, spheres still sit above solid cylinders and solid
cylinders above hollow ones: the ratios at those worst corners are
1.083 <!-- claim: pe_limit_spheres_over_solid_worst_rounding --> and
1.335 <!-- claim: pe_limit_solid_over_hollow_worst_rounding -->.

## Agreement

| Check | Result |
|---|---|
| high-Re limit, closed form vs eq. (1) at `Re = 1e8` | below 1e-7 <!-- claim: pe_limit_vs_eq1_worst_rel_dev --> on all three shapes |
| eq. (8), a worked example with every input printed | `1.65 × 0.37 / 0.6` against the printed `1.02`: -0.245 % <!-- claim: wall_velocity_ratio_rel_dev --> |
| seven numbers the paper states in words | 5 <!-- claim: stated_numbers_within_half_unit --> of the seven within half a unit of their last printed digit; worst deviation 11.43 % <!-- claim: stated_numbers_worst_rel_dev --> |
| eq. (1) [1987] vs eq. (42) [1993] | below 1e-9 <!-- claim: eq1_vs_eq42_identity -->, an identity of two readings |
| eq. (1)+(3) vs the markers `A2.5` digitised | 20.30 % <!-- claim: eq1_vs_fig9_mean_abs_dev --> mean absolute deviation, bias 15.77 % <!-- claim: eq1_vs_fig9_bias -->, 4.75 <!-- claim: eq1_vs_fig9_bias_over_se --> standard errors from zero |
| E&R's own eq. (18) on the same rows, recomputed | 11.75 % <!-- claim: eq18_vs_fig9_mean_abs_dev_recomputed -->, equal to `A2.5`'s committed claim |
| null: `p` a free constant | 22.96 % <!-- claim: null_free_constant_p_mean_abs_dev --> at `p` = 0.0793 <!-- claim: null_free_constant_p_value --> |
| `p(Re)` against that null, ε = `0.37` | 0.89 <!-- claim: p_cost_vs_free_constant_over_se --> standard errors of the paired difference: not resolved |
| the same, each model at its preferred porosity | 2.59 <!-- claim: p_cost_at_preferred_eps_over_se --> standard errors |
| the porosity, ε = `0.34` vs `0.40` | 26.62 % <!-- claim: eq1_vs_fig9_mad_eps_034 --> vs 17.55 % <!-- claim: eq1_vs_fig9_mad_eps_040 --> |
| 2-D bed vs Gunn's eq. (17), flat inlet | 7.454e-7 <!-- claim: eq17_flat_sup_error --> at `n_z = 800`, second order |
| 2-D bed vs a Bessel-mode closed form derived here | 2.896e-5 <!-- claim: bessel_mode_sup_error --> at `(800, 80)`, second order on each axis |
| wall channel: `1.65` read as interstitial | 5.461 % <!-- claim: wall_channel_cost_without_radial --> more unconverted with no radial mixing, 1.142 % <!-- claim: wall_channel_cost_with_radial --> with eq. (9)'s |
| zero-gradient outlet instead of eq. (18) | 1.035 % <!-- claim: zerograd_outlet_error_long_bed --> at `L/d = 200`, the same at `L/d = 5` |
| the same two, recomputed by `solve_bvp` collocation | outlet 1.0353 % <!-- claim: outlet_error_collocation -->, finite-volume gap -3.5e-6 <!-- claim: outlet_error_fv_vs_collocation --> relative; wall channel gap 1.5e-6 <!-- claim: wall_cost_fv_vs_collocation --> |

Deviations are `(model − reference)/reference` everywhere on the page.

## Are these markers inside Gunn's fit set? Probably not

- Gunn 1987 p. 365: `p` was estimated from experiments "in which the effect of
  molecular diffusion was very small".
- Gunn 1993 p. 335 says what that meant, `ReSc ≫ 1`, i.e. **liquid-phase**
  measurements, with gas-phase results used "to compare".
- The markers are gas phase at `Sc = 0.77`: `ReSc` runs
  0.0066 <!-- claim: resc_min --> to 37.8 <!-- claim: resc_max -->, median
  1.40 <!-- claim: resc_median -->, with 26 <!-- claim: n_resc_below_1 --> below
  one and 47 <!-- claim: n_resc_below_10 --> below ten. They fail that criterion.
- The 1987 paper attributes eqs (4) and (5) to England & Gunn (1970) and eq. (3) to
  nothing, and says of *these* points that they "are well supported by eq. (1)",
  the language of comparison.

So the comparison is **probably out of sample, though not provably held out**, and
it is not called a validation. What the markers resolve, measured:

- **`p(Re)` against a free constant `p`.** At the assumed porosity `p(Re)` is
  2.65 % <!-- claim: p_cost_vs_free_constant --> (percentage points) better, which is
  0.89 <!-- claim: p_cost_vs_free_constant_over_se --> standard errors of the
  paired difference: the markers do not distinguish the two. Given each the
  porosity it prefers, the difference is 4.62 % <!-- claim: p_cost_at_preferred_eps -->
  (percentage points), 2.59 <!-- claim: p_cost_at_preferred_eps_over_se --> standard errors.
  Which comparison applies is set by a porosity the paper does not print.
- **The porosity moves the headline more than either**, over ε = `0.34` to `0.40`, and
  the residual has an interior minimum at ε = 0.4429 <!-- claim: eq1_vs_fig9_best_eps -->
  (16.38 % <!-- claim: eq1_vs_fig9_mad_at_best_eps -->).

## What the page adds

**A number on Gunn's wall argument.** Eqs (7) and (8) conclude the interstitial
velocity in the wall region is close to the bulk ("we find that U_w/U_b is 1.02"), "suggesting that there may
be no significant differences". A 2-D bed says how much: reading Price's
superficial `1.65` as interstitial raises unconverted reactant by
5.461 % <!-- claim: wall_channel_cost_without_radial --> with no radial mixing and
1.142 % <!-- claim: wall_channel_cost_with_radial --> with eq. (9)'s coefficient
(160 radial cells, 0.088 % <!-- claim: wall_cost_grid_error --> above its
Richardson estimate), so radial dispersion erodes
79.1 % <!-- claim: radial_erosion_fraction --> of the channel; Gunn's corrected
ratio moves it by 0.0009 % <!-- claim: wall_cost_gunn_ratio_with_radial -->.

**A closed form for the 2-D problem.** Feeding the bed a zero-flux radial
eigenmode makes eq. (17) exact again with `k → k + D_R(β₁/R)²`. It is the only
check on the page that reads the radial operator: the flat-inlet check against
eq. (17) leaves `nu=0` at its baseline value, while the 2-D one moves by a factor
2646 <!-- claim: bessel_break_nu0_factor -->.

**The cost of the wrong outlet condition.** The zero-gradient outlet's error at
the *outlet* is a length-independent multiplicative bias of size `λD/U`:
1.035 % <!-- claim: zerograd_outlet_error_short_bed --> at L/d = 5 and the same, to within
1e-5 <!-- claim: zerograd_outlet_error_length_spread --> relative, at L/d = 200. What length buys is the *interior*: 4 diameters from the inlet the same
error is 0.150 % <!-- claim: zerograd_error_at_probe_5d_bed --> in a 5-diameter
bed, 9.65e-8 <!-- claim: zerograd_error_interior_10d_bed --> in a 10-diameter bed
and round-off beyond. Over a 400-fold sweep in `kL/U` the outlet error is
0.912 <!-- claim: zerograd_error_over_lambda_D_U_min --> to
0.9997 <!-- claim: zerograd_error_over_lambda_D_U_max --> times `λD/U`.

**A second, independent computation of two of those numbers** by
`scipy.integrate.solve_bvp`, which also shows what reading the outlet at the last
cell centre instead of the face costs: -11.42 % <!-- claim: outlet_error_centre_vs_collocation_h8 -->
at `h = d/8`.

**The solver.** The bed is linear apart from the van Leer limiter; it is solved by
`ConvFluxTVD` with `newton`, and checked against van Leer deferred correction on
the same discrete equations (gap below 1e-12 <!-- claim: newton_vs_deferred_gap -->).
Outlet values are read with `compute_boundary_values` under the outlet's own bc.

## Honest limits

- **The porosity is an assumption**, from a sentence on journal page 368 about
  the bulk of a bed of spheres. It is the largest lever on every comparison with
  the markers.
- **No liquid-phase data.** The Schmidt-group dependence, which Gunn's Figure 1
  exists to show, is untested. Jacques & Vermeulen (1957) or Miller & King (1966)
  would supply it.
- **The borrowed dataset carries `A2.5`'s open questions.** Three glyphs remain
  fused; including the nearly resolvable one moves the headline to
  20.65 % <!-- claim: eq1_vs_fig9_mad_fused_included -->; the ink-centroid
  systematic moves it to 19.27 % <!-- claim: eq1_vs_fig9_mad_shift_plus5px --> or
  21.59 % <!-- claim: eq1_vs_fig9_mad_shift_minus5px -->. No new digitisation was
  done.
- **The two closed-form solver checks cannot see a coefficient that is wrong
  everywhere**: an error of a tenth in `D_R` made consistently in operator and closed form
  moves the 2-D check by 6.3 % <!-- claim: bessel_consistent_DR_move -->; the same
  error in the operator alone by a factor 445 <!-- claim: bessel_break_DR_operator_only_factor -->.
- **One stated range does not come out**: the upper end of "0.2 < Re < 300" is not
  a crossing. The liquid/gas ratio has 1 <!-- claim: liquid_gas_crossings -->
  crossing and is still 1.084 <!-- claim: liquid_gas_ratio_at_Re_300 --> at
  Re = 300.
- **Three printing inconsistencies in the source are recorded, not repaired**: the
  Figure 5 caption names "eq. (4)" where eqs (9) and (10) are required; the
  tortuosity for spheres is `1.4` in eq. (3) and `1.2` in eq. (10); and eq. (2) prints
  `16α₁⁴(1−ε)` where eq. (1) prints `16α₁⁴(1−ε)²`, which Gunn 1993's eq. (42)
  confirms.
- **The reactor operating point is chosen, not printed.** Its Danckwerts inlet is
  written with the area-mean velocity; the local `u(r)` changes the wall-channel
  answer by 1.05 % <!-- claim: local_inlet_effect --> of itself.

## Corrections

- 2026-11-09: said the outlet metric read at the last cell centre "The metric came out 13 % low". Against the collocation reference it was -11.42 % <!-- claim: outlet_error_centre_vs_collocation_h8 --> at `h = d/8` and -2.97 % <!-- claim: outlet_error_centre_vs_collocation_h32 --> at `h = d/32`.
- 2026-11-09: said "Deleting `p(Re)` costs 2.7 percentage points against the strongest null", a difference of 0.89 <!-- claim: p_cost_vs_free_constant_over_se --> standard errors of the paired difference: at the assumed porosity the markers do not distinguish `p(Re)` from a free constant.
- 2026-11-09: said the seven stated numbers reproduce "each inside the precision of the phrase it tests". Counted in half-units of the last printed digit, 5 <!-- claim: stated_numbers_within_half_unit --> are; `p(1000)` deviates by -1.6 % <!-- claim: p_at_Re_1000_rel_dev --> from its plateau, and the liquid/gas maximum is 10.89 <!-- claim: liquid_gas_max_ratio --> against "a factor of 10".

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed. The two Gunn datasets are committed CSVs; the measured points are
loaded from `pages/A2.5-edwards-richardson-dispersion/data/` through
`load_data(..., page=...)`.
