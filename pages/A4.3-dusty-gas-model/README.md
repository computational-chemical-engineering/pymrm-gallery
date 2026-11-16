# A4.3: The dusty gas model

Pore diffusion in its *n*-component form: the matrix that one effective
diffusivity per species is an approximation of, plus the viscous term that
page `A4.4` sets to zero.

**Source.** Every equation and constant on the page is taken from **Krishna, R. &
Wesselingh, J. A. (1997)**, "The Maxwell–Stefan approach to mass transfer",
*Chemical Engineering Science* **52**(6) 861–911,
[doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7),
read from renders of the printed page (the PDF text layer is an OCR that mangles
eq. (86)'s prefactor and eq. (109)'s stoichiometric subscript). The origins of
the model (Mason & Malinauskas 1983, Jackson 1977, Knudsen 1909) and Kaza &
Jackson (1980), Haynes (1978) and Remick & Geankoplis (1974) are cited through
the review.

## Headline

Three things the matrix form says that no scalar effective diffusivity can
express, one that quantifies a sentence the review leaves qualitative, and one
theorem.

1. An **inert species** in a pellet running A → 2B carries zero flux
   (conservation, not a prediction) with a concentration that varies between the
   centre and the surface by 24.4 % <!-- claim: inert_concentration_variation_kn1 -->
   at Kn = 1 and 45.1 % <!-- claim: inert_concentration_variation_kn0p1 --> at
   Kn = 10⁻¹. A scalar closure **with a finite positive diffusivity**, run on the
   same pellets, gives a variation of at most 1e-12 <!-- claim: inert_scalar_closure_variation -->. The review's own eq. (109)
   degenerates to D → 0 for an inert and predicts nothing there.
2. **Uphill diffusion has no critical pore size.** The pore diameter at which N₂'s
   flux reverses depends on the near-end composition and moves by a factor
   4650 <!-- claim: uphill_threshold_range_over_window -->, from
   4246 <!-- claim: uphill_threshold_nm_at_x_0p483 --> nm to
   0.913 <!-- claim: uphill_threshold_nm_at_x_0p4999 --> nm, across compositions
   that all lie inside the uphill window at 1 mm (x_N2 from
   0.4823 <!-- claim: uphill_window_lower_1mm --> to
   0.5002 <!-- claim: uphill_window_upper_1mm -->). What pore size controls is the
   *width of the composition window*: its half-width falls from
   1.77e-2 <!-- claim: uphill_window_halfwidth_1mm --> at 1 mm to
   1.09e-4 <!-- claim: uphill_window_halfwidth_1nm --> at 1 nm, a log–log slope of
   0.976 <!-- claim: uphill_window_width_slope_in_d0 --> in d₀. Every edge and
   threshold is a root found by a scan of the whole physical range.
3. **Jackson's** *p*₀ = √ν_B *p* is the corner case of
   *p*₀/*p* = 1 + (ν_B √(M_B/M_A) − 1)·x_A,s, needing pure A at the surface,
   mass-conserving stoichiometry **and** Knudsen control. Under bulk control the
   rise nearly vanishes (1.0023 <!-- claim: jackson_p0_bulk_kn0p01 --> at
   Kn = 10⁻²). At finite Kn the centre falls short of √2 by 1/(2 Kn): Kn times
   the departure, extrapolated to zero cell size, is
   -0.4999 <!-- claim: jackson_dev_times_kn -->.
4. **"Neglect of the viscous flow contribution is not very serious (Haynes,
   1978)"** becomes a screening group from the review's eqs. (85), (91) and (106):
   D_visc/Ðᵉ_iM = 3 *p* d₀/(32 η v̄_i), 2.64e-4 <!-- claim: visc_group_1nm --> for
   H₂ in a 1 nm pore and 0.264 <!-- claim: visc_group_1um --> in a 1 µm one. It
   cannot be swept alone: Kn depends on the same d₀, running the other way
   (138.9 <!-- claim: kn_1nm --> at 1 nm, 0.139 <!-- claim: kn_1um --> at 1 µm).
   Tying both to d₀, viscous flow changes the pellet-centre pressure by
   0.0092 % <!-- claim: viscous_relief_1nm_pore --> at 1 nm and
   0.634 % <!-- claim: viscous_relief_1um_pore --> at 1 µm, and is largest,
   0.640 % <!-- claim: viscous_relief_worst_over_pore_sizes -->, at
   744 <!-- claim: viscous_relief_worst_d0_nm --> nm where
   Kn = 0.187 <!-- claim: viscous_relief_worst_at_kn --> (all extrapolated to zero
   cell size). Freezing Kn at 10⁴ instead gives
   6.90 % <!-- claim: viscous_relief_fixed_kn_artefact --> at the 1 µm viscous
   group, 10.9 <!-- claim: viscous_fixed_kn_overstatement --> times what a 1 µm
   pore shows.
5. **[Bᵉ] is invertible for a reason.** Σᵢ Bᵢⱼ = 0 for the wall-free matrix because
   Ð is symmetric, so (1,…,1) is a left null vector for any composition and any
   n ≥ 2; the wall term makes every column sum 1/Ðᵉ_jM > 0 and the matrix strictly
   column-diagonally dominant. Checked on the assembled `b_matrix` for n = 2 to 6:
   null vector at most 1e-14 <!-- claim: bmatrix_wall_free_null_vector -->, column
   sums 1/Ðᵉ_jM to at most 1e-14 <!-- claim: bmatrix_colsum_knudsen_dev -->, and
   deleting the wall term in the solver leaves no solution.

**Data tier 6**: nothing is compared with a measurement.

## Validation

| check | result |
|---|---|
| independent collocation solution (shares no transport code), 100 cells | at most 1e-5 <!-- claim: collocation_worst_dev --> over six decades of pore size |
| the same, finite-volume flux extrapolated to zero cell size | at most 1e-9 <!-- claim: collocation_extrapolated_dev --> |
| grid order (successive differences, 25 to 800 cells) / 100-cell error | 2.01 <!-- claim: grid_order --> at worst / at most 1e-6 <!-- claim: grid_error_n100 --> |
| eq. (103) recovered, viscous term included | at most 1e-10 <!-- claim: eq103_viscous_dev --> |
| collapse to `A4.4`'s eq. (109) | at most 1e-10 <!-- claim: eq109_collapse_worst --> |
| Graham's law, predicted not imposed (an algebraic identity) | at most 1e-9 <!-- claim: graham_residual_worst --> |
| free-molecule limit at d₀ = 10⁻² nm, and its rate in d₀ | 1.76e-4 <!-- claim: free_molecule_limit_dev -->, 1.00 <!-- claim: free_molecule_rate --> |
| Maxwell–Stefan limit at d₀ = 1 mm, and its rate in d₀ | at most 5.8e-4 <!-- claim: maxwell_stefan_limit_dev -->, -1.00 <!-- claim: maxwell_stefan_rate --> |
| the four p. 887 exponents | worst departure 0.0205 <!-- claim: p887_exponent_worst_dev --> |
| eqs. (105), (107), (103) derived from eq. (101) symbolically | nonzero residuals: 0 <!-- claim: symbolic_nonzero_residuals --> |

The collocation route shares `knudsen_D` and the pair products with the solver,
so it cannot catch a mistake inside those. Eq. (85)'s mass exponent is pinned by
Graham's law; its d₀/3 prefactor and eq. (83)'s ε/τ are pinned by **nothing on the
page**. That is measured: breaking d₀/3 → d₀/2 in every route at once moves the
collocation gap by a factor of at most 20 <!-- claim: prefactor_v1_ratio --> (the change of
the 100-cell discretisation error with Kn) and leaves eq. (103) within 1e-10 <!-- claim: break_eq85_prefactor_v2 -->
and Graham's law within 1e-9 <!-- claim: break_eq85_prefactor_v3 -->.

Every claim has a defect row in the page's table that fails it: eq. (87)'s
off-diagonals, the wall term, eq. (85)'s mass exponent and prefactor, a pair
diffusivity, eq. (91)'s 32, α′, the face-composition rule, the Newton tolerance,
the scalar closure, the stoichiometric mass ratio, the grids of a study.

**Scope.** The scalar closure and the Bosanquet relation are `A4.4`; bulk
multicomponent diffusion with no wall is `A4.2`; uphill diffusion measured in
free space is `A4.9`. Not implemented here, though the review discusses all of
them: surface diffusion, pore-size distributions, viscous selectivity (α′ᵢ ≠ 1),
thermal diffusion, non-ideal mixtures. The review's Fig. 44 (the dusty gas model
against Remick & Geankoplis's He/Ne/Ar capillary measurements) is not used.

## What pymrm adds

The transport skeleton of the n-component solver: `construct_grad` for the face
gradients, `construct_div` with `nu` for slab, cylinder or sphere,
`construct_boundary_value_matrices` for the boundary face compositions, `NumJac`
with a cell-and-species stencil, `newton` with bounds; and `pymrm.checks` for the
root scans and the refinement ladders. The flux law (a dense solve of eq. 87 at
every face) and the studies are the page's own.

## Files

- `index.py` / `index.ipynb`: the page (paired by jupytext).
- `model.py`: `DustyGas`, the capillary and pellet set-ups, and the two
  collocation routes (the dusty gas equations and Maxwell–Stefan with Graham's
  law).
- `test_model.py`: one capillary solve against the collocation route and the
  stored deviation, and the Knudsen-limit pellet.
- `data/krishna-wesselingh-1997-dgm-constants.csv`: the structural constants of
  eqs. (83) to (85), (87), (91), (101), (103), (107), plus the four exponents the
  review states in words on its p. 887, with a provenance sidecar. Molar masses
  are IUPAC and flagged as not from the review.
- The pair diffusivities are loaded cross-page from the `A4.2` dataset (the
  review's p. 872 values).

## Corrections

- 2026-11-16: "a concentration that varies by 22.5 % at Kn = 1 and 36.6 % under bulk control"
  read the inert's concentration in the cells next to the centre and the surface,
  half a cell short of each boundary, where the profile is steep; those numbers
  moved with the grid. On the boundary faces the variation is
  24.4 % <!-- claim: inert_concentration_variation_kn1 --> at Kn = 1 and
  45.1 % <!-- claim: inert_concentration_variation_kn0p1 --> at Kn = 10⁻¹.
- 2026-11-16: "residual at Kn = 10⁴ is exactly 1/Kn"
  equates the departure from √2 with 1/Kn. It falls as 1/Kn, but Kn times it,
  extrapolated to zero cell size, is -0.4999 <!-- claim: jackson_dev_times_kn -->, so
  the departure is −1/(2 Kn).
- 2026-11-16: "pushing φ does not shrink it (−5.09 × 10⁻⁵ at φ = 200)"
  quoted a value whose difference from φ = 60 was discretisation error of the
  200-cell grid. Extrapolated to zero cell size, φ = 200 changes the departure by
  at most 1e-3 <!-- claim: jackson_phi_effect --> of itself.
- 2026-11-16: "0.62 % at 1 µm, and never by more than 0.63 % anywhere from 1 nm to 100 µm"
  and "the maximum sits at 562 nm, Kn = 0.25" were read on a 200-cell grid and at a
  sampled pore diameter. Extrapolated to zero cell size, the change is
  0.634 % <!-- claim: viscous_relief_1um_pore --> at 1 µm and largest,
  0.640 % <!-- claim: viscous_relief_worst_over_pore_sizes -->, at
  744 <!-- claim: viscous_relief_worst_d0_nm --> nm, Kn =
  0.187 <!-- claim: viscous_relief_worst_at_kn -->, found by a bounded maximisation of
  the extrapolated change.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```
