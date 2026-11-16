# A4.4: Knudsen diffusion and the Bosanquet relation

What happens when the pore is small enough that molecule–wall collisions
compete with molecule–molecule collisions, and where the Bosanquet
interpolation `1/D_eff = 1/D_K + 1/D_AB` is right.

**Source.** Every equation on the page is read from **Krishna, R. &
Wesselingh, J. A. (1997)**, "The Maxwell–Stefan approach to mass transfer",
*Chemical Engineering Science* **52**(6) 861–911,
[doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7),
which prints the dusty gas model (eqs. 82–87), derives the Bosanquet formula as
its eq. (110), names it, and states the three conditions it needs. The equations
were read off 600 dpi renders of the printed page (the PDF text layer is an OCR
that renders eq. (86)'s prefactor as `c_t/RT` where the page prints `1/RT`, and
eq. (109)'s `nu_1` as `nu_{1i}`). The origins, Knudsen (1909) and Bosanquet's
1944 British report BR-507, were **not consulted**.

## Headline

1. **A closed form for the error the review does not quantify.**
   `D_Bos/D_1 = 1 + (rho - 1) x_1 / (1 + Kn)`, derived symbolically from the
   review's eqs. (82), (109) and (110) with 0 <!-- claim: symbolic_nonzero_residuals -->
   nonzero residuals, where `rho` is the flux ratio the application imposes. The
   relation is exact under Knudsen control and worst under bulk control, the
   opposite of the instinct that an interpolation is worst in the middle.
2. **The effectiveness factor of a reacting pellet.** A sphere running `A -> 2B`
   at a Thiele modulus of 30: Bosanquet overestimates the effectiveness factor
   by 26.96 % <!-- claim: bosanquet_eta_error_bulk --> at Kn = 1e-3,
   14.12 % <!-- claim: bosanquet_eta_error_kn1 --> at Kn = 1 and
   2.97e-7 <!-- claim: bosanquet_eta_error_knudsen --> (relative) at Kn = 1e6,
   each extrapolated to zero cell size from a ladder of 300 to 4800 cells. The
   error falls at every step of the sweep in Kn
   (18 <!-- claim: bosanquet_eta_error_decreasing_in_kn --> steps on the 600-cell sweep). The numbers
   are conditions-specific: the bulk-limit error is
   25.49 % <!-- claim: bosanquet_eta_error_bulk_phi10 --> at phi = 10 and
   27.43 % <!-- claim: bosanquet_eta_error_bulk_phi100 --> at phi = 100, rising
   with phi at both steps (2 <!-- claim: bosanquet_eta_error_grows_with_phi -->).
3. **The isobaric H₂/N₂ counter-diffusion cell.** Graham's law fixes the flux
   ratio at `sqrt(M_A/M_B)`, so Bosanquet's flux error is
   -41.33 % <!-- claim: bosanquet_flux_error_isobaric_bulk --> at Kn = 1e-3 and
   -0.037 % <!-- claim: bosanquet_flux_error_isobaric_knudsen --> at Kn = 1e3
   (eq. 109 with Graham's ratio integrated across the slab in closed form).
4. **The internal pressure.** At Kn = 1e-3 the dusty gas model puts the pellet
   centre at 1.000287 <!-- claim: dgm_centre_pressure_bulk --> times the surface
   pressure and the Bosanquet closure at
   1.99941 <!-- claim: bosanquet_centre_pressure_bulk -->; under Knudsen control
   both sit on Jackson's `sqrt(nu_B)`.

## Validation

Tier 6: no measurement anywhere on this page.

| check | result |
|---|---|
| Jackson's `p_0 = sqrt(nu_B) p` (K&W p. 892, "a 40 % increase" for `nu_B = 2`), at Kn = 1e6 | 41.421 % <!-- claim: jackson_nu2_solved_increase --> increase; departure 1.08e-6 <!-- claim: jackson_centre_pressure_dev --> at worst over `nu_B` = 1 to 4 |
| that departure is the finite-Kn term of the model | departure × Kn varies by at most 1e-2 <!-- claim: jackson_dev_x_kn_spread_1e5_to_1e8 --> (relative) from Kn = 1e5 to 1e8 on a fixed grid |
| invariance over geometry, Thiele modulus and grid | spread at most 1e-6 <!-- claim: jackson_invariance_spread --> |
| defects the Knudsen-limit check cannot see: `[B^e]` off-diagonals flipped, slab for sphere, five cells | shifts at most 1e-5 <!-- claim: jackson_shift_b_offdiag_flipped -->, 1e-9 <!-- claim: jackson_shift_slab --> and 1e-6 <!-- claim: jackson_shift_five_cells -->; the break table shows each at Kn = 1 or phi = 3 |
| Graham's law, eq. (107), predicted by the isobaric slab, not imposed (an identity of eqs. 82 and 85) | residual at most 1e-10 <!-- claim: graham_law_residual -->; flux ratio at most 1e-10 <!-- claim: graham_flux_ratio_dev --> from `sqrt(M_A/M_B)` |
| finite-volume slab (200 cells) against eq. (109) integrated in closed form | at most 1e-5 <!-- claim: isobaric_fv_vs_closed_form --> |
| p. 887 exponents as local slopes, Knudsen (d0 = 20 nm, 1 kPa) / bulk (20 µm, 1 MPa) | departures 9.00e-4 <!-- claim: scaling_slope_dev_knudsen --> / 1.30e-3 <!-- claim: scaling_slope_dev_bulk -->, falling with Kn at rates -1.00 <!-- claim: scaling_rate_knudsen --> / 1.00 <!-- claim: scaling_rate_bulk --> |
| the same at T = 400 K / with the borrowed pair diffusivities × 1.3 | Knudsen 1.03e-3 <!-- claim: scaling_slope_dev_knudsen_T400 --> / 6.93e-4 <!-- claim: scaling_slope_dev_knudsen_D1p3 -->; bulk 1.14e-3 <!-- claim: scaling_slope_dev_bulk_T400 --> / 1.69e-3 <!-- claim: scaling_slope_dev_bulk_D1p3 --> |
| finite-volume slopes against closed-form slopes | at most 1e-6 <!-- claim: scaling_slope_vs_closed_form --> |
| `[B^e]` route against eq. (109) with the stoichiometric ratio (an identity of the scheme) / against eq. (110) | at most 1e-12 <!-- claim: dgm_vs_eq109_correct_rho --> / at least 0.1 <!-- claim: dgm_vs_eq110_wrong_rho --> |
| grid order on the effectiveness factor, 50 to 6400 cells | 2.04 <!-- claim: grid_convergence_order --> (finest triple); finest difference at least 1e3 <!-- claim: grid_ladder_over_roundoff --> times the Newton tolerance |
| quoted errors: finest grid against the extrapolated value / what the 600-cell sweep misses at Kn = 1e-3 | at most 1e-4 <!-- claim: eta_extrapolation_worst --> / 1.29e-4 <!-- claim: bosanquet_eta_error_bulk_grid_error --> |
| geometry index against the classical effectiveness factor of slab, cylinder and sphere | at most 1e-5 <!-- claim: geometry_eta_vs_closed_form_dev --> |
| Bosanquet centre pressure against its closed form `nu_B D_A/D_B` | at most 1e-10 <!-- claim: bosanquet_p0_vs_closed_form --> |
| identities: eq. (85) scalings / closed-form error map against the solver's coefficients | at most 1e-14 <!-- claim: eq85_scaling_dev --> / 1e-14 <!-- claim: closed_form_error_map_dev --> |

Every claim fails at least one row of the defect-injection table on the page.

**Scope.** Binary and scalar (`S3`) by design. The n-component matrix form of the
dusty gas model, the viscous-flow term, and the review's Fig. 44 comparison with
the Remick & Geankoplis capillary measurements are page `A4.3`. Bulk
multicomponent diffusion with no pore is `A4.2`, built from the same review.

## What pymrm adds

The transport skeleton of the binary solver: `construct_grad` for the face
gradients, `construct_div` with `nu` for slab, cylinder or sphere,
`construct_boundary_value_matrices` for the boundary face compositions (the
pellet-centre pressure is read there), `NumJac` with a cell-and-species stencil,
`newton` with bounds; and `pymrm.checks.observed_orders` for the refinement
ladders. The three flux closures, the closed forms and the studies are the
page's own.

## Files

- `index.py` / `index.ipynb`: the page (paired by jupytext).
- `model.py`: `Cell` with the three flux closures, the pellet and isobaric-slab
  set-ups, the closed-form slab flux and the classical effectiveness factor.
- `test_model.py`: one pellet solve against Jackson's limit and the stored value,
  and the slab against its closed form.
- `data/krishna-wesselingh-1997-porous-media.csv`: the printed constants and
  stated results of the review's dusty-gas sections (tier 6; provenance in the
  sidecar). The review's pair diffusivities are loaded cross-page from the `A4.2`
  dataset.

## Corrections

- 2026-11-16: "the same error is 25.5 % at phi = 10 and 27.7 % at phi = 100"
  quoted values read on the 600-cell sweep grid. Extrapolated to zero cell size the error is
  25.49 % <!-- claim: bosanquet_eta_error_bulk_phi10 --> at phi = 10 and
  27.43 % <!-- claim: bosanquet_eta_error_bulk_phi100 --> at phi = 100.
- 2026-11-16: "the sign flip instead caught by the closure comparison (1.0e0)"
  rested on a solve that had not converged: with the off-diagonals of `[B^e]`
  flipped, the pellet at Kn = 0.01 stops on the non-negativity bound and Newton
  does not report success. Every solve on the page raises unless Newton reports
  success. A smaller defect in the off-diagonals that does converge (a row of the
  break table) moves the closure comparison outside its bound of
  1e-12 <!-- claim: dgm_vs_eq109_correct_rho -->.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```
