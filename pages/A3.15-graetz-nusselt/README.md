# A3.15: Graetz–Nusselt problem

The developing thermal boundary layer in a tube, solved as a 2-D PDE, and the
two amplitudes Graetz declared out of reach and calibrated out of his own
apparatus instead.

- **Structure:** `S6` (2-D PDE)
- **Reference:** Graetz (1882), *Ueber die Wärmeleitungsfähigkeit von
  Flüssigkeiten. Erste Abhandlung*, Annalen der Physik **254**(1) 79–94,
  doi:10.1002/andp.18822540106
- **Data:** tier 1, 51 measured runs on six liquids, transcribed from the
  300 dpi page bitmaps

## What the paper does, and where it stops

Graetz's problem was not heat transfer; it was measuring the conductivity of
liquids with a flow calorimeter. To read his instrument he had to solve

$$2\alpha(1-r^2/R^2)\,\partial u/\partial z = a^2\,\big(u_{rr} + u_r/r\big),$$

which is the problem now named after him. He gets the eigenvalues, and then
stops:

> "Es wäre natürlich von wesentlichem Vortheil, wenn sich die $p_i$ bestimmt
> durch $R$ und $\mu_i$ ausdrücken liessen. Indess gehört dazu eine eingehende
> Untersuchung der Function $V(r,\beta,R)$, die wahrscheinlich nicht einfach
> ist … Für physikalische Zwecke lassen sich dagegen die $p_i$ sehr einfach als
> Constanten des Apparates experimentell bestimmen."

He fits $p_1$, $p_2$ **and** $k/c$ to three runs at three pressure heads.

## How it is solved

`GraetzTube` (`model.py`) puts the parabolic problem on the pymrm operators:
`construct_div(nu=1)` radially, van Leer face values axially with the exact
Jacobian from `ConvFluxTVD`, and `newton(damping="monotone", stagnation=True)`.
Van Leer deferred correction on one LU of the upwind operator solves the same
equations as a second route; it lands on the Newton field to below
1e-10 <!-- claim: newton_vs_deferred_field --> and is faster on these grids, but
only under-relaxed: at $\omega = 1$ it settles into a period-two cycle. The
eigenfunction route (shooting plus projection) shares no code with the solve.

## Agreement

Against what the paper prints:

| Check | Result |
|---|---|
| his power series for $V(\mu)=0$, vs exact rationals | 3 <!-- claim: series_coeffs_right_to_printed_digits --> of his 5 non-trivial coefficients right to every printed digit; $\mu^8$ off by 0.58 % <!-- claim: series_coeff_mu8_printed_err -->, $\mu^{10}$ by a factor of 2.65 |
| $\mu_1 = 2.7043$ | 2.4e-5 <!-- claim: mu1_printed_rel_err --> |
| the four $J_0$ zeros of his plug-flow case | worst 4.3e-5 <!-- claim: bessel_zeros_printed_max_err --> |
| $\mu_2 = 6.50$ | **2.7 % <!-- claim: mu2_printed_rel_err --> low**; the 10-term series cannot place it |
| his calibration re-solved, $\log p_1$ | 2.3e-5 <!-- claim: calibration_logp1_abs_err --> |
| his working formula vs his own $k/c$ column, 51 runs | 0.71 % <!-- claim: formula_reproduces_kc_mean_dev --> mean |

Against limits he did not have, neither of which the series supplies:

| Check | Result |
|---|---|
| Lévêque constant $2/(\Gamma(4/3)\,9^{1/3}) = 1.076732$, derived | 2.9e-4 <!-- claim: leveque_const_rel_err --> at 400 x 160 cells, order 2.0 |
| $\mathrm{Nu}_\infty = \mu_1^2/2$ | 3.65679 <!-- claim: nu_infinity --> |
| 2-D solve vs eigenfunction projection (no shared code) | $\mu_1$ 5.5e-6 <!-- claim: mu1_rel_err_2d -->, $p_1$ 4.2e-4 <!-- claim: p1_rel_err_2d -->, $p_{1,\text{area}}$ 4.3e-4 <!-- claim: p1_area_rel_err_2d --> |
| observed order in $\mu_1$, axial / radial (successive differences) | 1.94 / 2.00, both bounds checked on every run |

**The 5.5e-6 is not the accuracy of the solve.** At $n_z = 600$, $n_r = 128$ the
axial error is +5.9e-5 <!-- claim: mu1_err_axial_component_n600 --> and the radial
-5.3e-5 <!-- claim: mu1_err_radial_component_n128 -->, each measured on its own
ladder without the published grid; they cancel. Their sum reproduces the error of
the two grids not used to build them to below 2e-6. The solver's accuracy there is
that of the larger component.

## Two printed defects, both settled from the paper's own numbers

- **The exponent in the working formula on page 90 reads 2.4013**, where page 87
  derives 2.7043. Inverting the formula for all 51 runs decides it: 0.71 % mean
  deviation from his own $k/c$ column with 2.7043, 26.3 % <!-- claim: formula_misprint_kc_mean_dev -->
  with 2.4013. (2.401 is the brass-wall-corrected Bessel root four pages
  earlier, a plausible compositor's substitution, but the data settle it, not the
  guess.)
- **Two series coefficients are arithmetically wrong**, found against an
  independent exact-rational solution of his own recurrence. One caution,
  recorded on the page and in the data sidecar: the third digit-group of the
  $\mu^{10}$ coefficient is a **broken glyph** on the 300 ppi bitmap. At 16× it is
  a **3**, which gives the factor of 2.65; had it been a 9 the printed value would
  sit 2.4 % from exact and the claim would collapse to an ordinary rounding slip.
  The $\mu^{8}$ slip does not depend on it.

## The result, and it is not the expected one

His calibration has three constants and three observations, so it fits exactly
and proves nothing. Substituting the computed amplitudes removes the fit, and
those three runs then **no longer agree on a single $k/c$: they spread by
15.6 % <!-- claim: kc_spread_calibration_unfitted -->, monotonically with flow
rate**, against 0.77 % <!-- claim: kc_spread_thermometry_sigma --> propagated from
the half-division reading error of the 1/10-degree thermometers he describes:
20.4 <!-- claim: kc_spread_significance --> sigma on his *stated* thermometry. That
budget carries E₁ and E₂ only; adding the bath thermometer, a 1 % weighing error
and dropping the √5 averaging brings it down to
8.0 <!-- claim: kc_spread_significance_wide --> sigma. **The honest headline is
that range.**

**And it holds out of sample, which needs no error model at all.** Graetz sets
this test himself on page 89 (whether it matters which pressure head you use) and
asserts on page 90 that it passes. Run on every liquid he measured at more than
one head (48 <!-- claim: runs_never_fitted --> of the 51 runs were never fitted
to anything):

| liquid | runs | his fitted $p_i$ | computed $p_i$ |
|---|---|---|---|
| water | 13 | 3.6 % <!-- claim: kc_spread_out_of_sample_water_fitted --> | **17.0 % <!-- claim: kc_spread_out_of_sample_water -->** |
| copper sulfate | 11 | 4.8 % <!-- claim: kc_spread_out_of_sample_cuso4_fitted --> | **20.9 % <!-- claim: kc_spread_out_of_sample_cuso4 -->** |
| alcohol | 7 | 7.2 % <!-- claim: kc_spread_out_of_sample_alcohol_fitted --> | **17.8 % <!-- claim: kc_spread_out_of_sample_alcohol -->** |

In all 5 <!-- claim: liquids_lowest_head_highest_kc --> liquids measured at more
than one head the lowest head gives the highest $k/c$. 20 modes give the same
three $k/c$ as 10 to below 1e-12 <!-- claim: kc_10_vs_20_modes -->, and inverting
through the 2-D solver's own $\theta_{\rm cup}(x)$ with *no eigenfunction
anywhere* gives 15.6 % <!-- claim: kc_spread_2d_route -->.

Removing the fit moves every conductivity down, by
16.6 % <!-- claim: k_drop_unfitted_least --> to
24.5 % <!-- claim: k_drop_unfitted_most -->, and destroys the agreement he had
with Lorberg's recomputation of Weber (rows `k_water_graetz_30C` and
`k_water_lorberg_carried_30C` of the constants file). So
the empirical calibration was absorbing a real systematic effect in the apparatus
that his own method could not see, and computing the constants makes his numbers
*worse*. The page does not diagnose the effect; see *What this page cannot
conclude*.

## Corrections

- 2026-10-31: gave the radial observed order of the 2-D solve as 2.12. It was
  read against the exact eigenvalue on a ladder at $n_z = 800$, whose own axial
  error (+3.3e-5, printed in §3b) sits in every rung. From successive
  differences, which cancel it, the order is 2.00 (2.01 and 2.00 on the
  triples); the axial order stays 1.94.
- 2026-10-31: stated that the solver's honest accuracy at the published grid is
  about 4.5e-5 in $\mu_1$, from refining one axis alone (600/256: 4.5e-5
  <!-- claim: mu1_rel_err_600x256 -->, 1200/128: 3.8e-5
  <!-- claim: mu1_rel_err_1200x128 -->). Those grids still cancel partly. The
  components at the published grid are +5.9e-5 (axial) and -5.3e-5 (radial).
- 2026-10-31: gave the order of the Lévêque ladder as 1.9. It was fitted against
  the exact constant, which the ladder does not reach: the Richardson limits of
  its last pairs remain 3e-5 to 4e-5 <!-- claim: leveque_richardson_rel_err -->
  from it. From successive differences the order is 2.0.
- 2026-10-31: said the axial-conduction shift in the 2-D solve gives C = -4.6091,
  0.7 % from the perturbation result. -4.6091 was C at Pe = 80 on 300 x 64 cells,
  fitted up to the outlet, where the zero-gradient condition (a real one once
  axial conduction is in) acts on the last cells. Fitted before the outlet, C
  agrees with the perturbation result to below 1e-3
  <!-- claim: axial_pert_rel_gap --> once extrapolated to large Pe.
- 2026-10-31: said first-order upwind needs 16 times the cells of the van Leer
  scheme for the same error. With orders 1 and 2 the factor grows with the
  accuracy asked: about 15 <!-- claim: upwind_cells_factor_n75 --> at the axial
  error of 75 van Leer cells, about 50 at 300.
- 2026-10-31: said the lowest head gives the highest $k/c$ on four liquids; it
  does on all five liquids measured at more than one head.
- 2026-10-31: said the amplitudes were computed two independent ways. Only
  $p_1$ (cup and area mean) has a second route, the far-field decay of the 2-D
  solve; $p_2$ to $p_{10}$ come from the projection alone, and the 2-D solve
  confirms the summed series at the calibration runs.

## Files

```
index.py           the page source (jupytext percent format), paired with index.ipynb
index.ipynb        the executed notebook
model.py           GraetzTube (2-D solve, two solver routes), the eigenfunction route,
                   the exact series, the k/c inversion
test_model.py      fast checks of model.py and two stored claims
meta.yaml          page metadata
agreement.json     the claims, format 2, written by the last cell
data/graetz-1882-observations.csv       all 51 runs of his Section 4
data/graetz-1882-printed-constants.csv  every constant printed in Sections 2-4
```

Runtime about 190 s.

## Cite the source, not this page

Graetz, L., *Ueber die Wärmeleitungsfähigkeit von Flüssigkeiten. Erste
Abhandlung*, Annalen der Physik **254**(1) 79–94 (1882).
The catalogue records this as "Graetz (1883)"; the paper is signed *München,
15. Oct. 1882*. The *zweite Abhandlung*, Ann. Phys. **261** (1885), is a
different paper and was not read here.
