# A2.3 — Taylor–Aris dispersion

Shear plus molecular diffusion produces spreading that looks exactly like
diffusion but is thousands of times faster and, counter-intuitively, gets
*slower* as the solute diffuses faster.

- **Structures:** `S3` (1-D steady BVP), `S6` (2-D PDE)
- **Reference:** Taylor (1953), Proc. R. Soc. Lond. A 219(1137) 186–203,
  doi:10.1098/rspa.1953.0139
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`closure`, the homogenisation closure; `Slug` and `run_history`,
  the 2-D transient with per-step moments and settling roots), tested by
  `test_model.py`. Taylor's closed forms and the modal series stay in the
  notebook.
- **Runtime:** about 210 s.

## Agreement

Two independent routes to Taylor's Eq. 25, `k = a²u₀²/(192 D)`:

| Route | Result |
|---|---|
| homogenisation closure (1-D, no convection), grid mean subtracted | 1.9e-9 <!-- claim: closure_rel_err_n200 --> relative at `n_r = 200`, observed order 4 |
| the same closure over Pe = 2 to 398 | an identity: `k ∝ Pe²` exactly, so the relative error is the grid error, `-3/n_r⁴`, at every Pe |
| direct 2-D transient, window τ = 0.6 to 1.2 | 1.112 % <!-- claim: direct_late_rel_err_n3200 --> at `n_z = 3200`, first order in Δz, extrapolated limit -2.4e-5 <!-- claim: direct_late_richardson_rel_err --> |
| the same, window τ = 0.2 to 1.0 | extrapolated limit 0.416 % <!-- claim: direct_richardson_rel_err --> below Eq. 25: the start-up deficit, 0.407 % <!-- claim: modal_window_deficit_early --> by the modal series |

Against the paper: Taylor's worked capillary run recomputes to 0.04 %
<!-- claim: taylor_k_rel_err --> in `k` and 0.38 % <!-- claim: taylor_D_rel_err -->
in the diffusivity it implies, and that value, 5.98e-6 cm²/s, lies inside the
4.35e-6 to 1.50e-5 range Furth & Ullmann measured independently.

## Watch the velocity convention

`u₀` is the **centreline** velocity, as in Taylor's paper. The mean is `u₀/2`,
so the equivalent form is `k = a²ū²/(48 D)`. Using the mean velocity with 192
gives an answer four times too small, the commonest way to get this wrong.

## What the page adds

Taylor bounded when his result applies (wait "much longer than" the radial
equilibration time `a²/(3.83²D)`) but left it as an inequality. The modal series
of the closure, derived on the page, gives the time after which the variance
growth rate stays within a band of Taylor's value:

| band | τ = tD/a² | equilibration times |
|---|---|---|
| ±10 % | 0.1498 <!-- claim: tau_settle_10pct_modal --> | 2.20 |
| ±5 % | 0.1970 <!-- claim: tau_settle_5pct_modal --> | 2.89 |
| ±3 % | 0.2318 <!-- claim: tau_settle_3pct_modal --> | 3.40 |
| ±2 % | 0.2594 <!-- claim: tau_settle_2pct_modal --> | 3.80 |
| ±1 % | 0.3066 <!-- claim: tau_settle_1pct_modal --> | 4.50 |

close to `ln(A/b)/3.83²` with `A` = 0.9014 <!-- claim: modal_amplitude --> at
Pe = 50: each factor e of accuracy costs one equilibration time. The 2-D pymrm
simulation reproduces these once its upwind numerical-dispersion floor (1.114 %
<!-- claim: numerical_dispersion_floor_n3200 --> of the rate at `n_z = 3200`,
equal to its closed form `Pe Δz/8 / (1 + Pe²/192)` to the radial error of the
24-cell grid) is removed: the corrected times lie 0.58 %
<!-- claim: tau_settle_3pct_corrected_vs_modal --> above the series in every band,
the bias of backward Euler and the radial grid in the fitted decay rate. Read off
a single mesh the same times are biased low (0.2115
<!-- claim: tau_settle_3pct_n3200 --> for ±3 % at `n_z = 3200`), and a band
narrower than the floor is never entered.

The series is scipy, not pymrm; pymrm supplies the operators of the closure and
the 2-D transient (`construct_grad`, `construct_div`, `construct_convflux_upwind`).

## Honest limits

Data tier 2, and the comparison is indirect: Taylor's raw concentration profiles
are not transcribed, only his worked inference and the independently measured
range it is checked against. His second inferred value (Eq. 42) is left out
because the tube radius for those runs is not stated beside the calculation, and
back-solving it implies a different tube.

The direct simulation uses first-order upwind; its numerical dispersion is
measured, checked against a closed form and removed from the settling times by a
correction, not by the scheme. The settling times hold for a slug released
uniformly across the section at Pe = 50. Every claim on the page has a row in
the defect-injection table that fails its own check, except one labelled
structural: the spread of the closure's relative error across Pe, zero by
construction because the closure is linear in the velocity.

The 2-D domain is closed: no diffusive flux and a face velocity of zero on the
two end faces, so the volume-weighted operator's column sums are zero to
round-off and mass is conserved for any state.

## Reproduce

```bash
cd pages/A2.3-taylor-aris-dispersion
jupytext --sync index.py      # index.py is the source; index.ipynb is paired with it
jupyter nbconvert --to notebook --execute --inplace index.ipynb
GALLERY_OFFLINE=1 python -m pytest -q test_model.py
```

No PDF needed: the dataset is a committed CSV.

## Corrections

- 2026-08-05: the settling times were stated from one unrefined mesh
  (`n_z = 1600`), with the ±3 % threshold at τ = 0.1950 and the ±2 % band
  "never settling". Both were the mesh's numerical-dispersion floor: the ±2 %
  band is entered, and the thresholds are higher (see the table above).
- 2026-08-05: the Reuse section said `τ ≳ 1`; the page's own result was `τ ≳ 0.2`.
- 2026-10-31: stated that the closure reproduces Eq. 25 to 1.0e-4 at
  `n_r = 200`, converging at second order, and that its singular system can be
  pinned at any cell. The exact mean `u₀/2` was subtracted, whose average on the
  grid is not zero, so the system was inconsistent: moving the pinned cell from
  the axis to the wall changes `k` by 1.5e-4 <!-- claim: closure_u0half_pin_gap_n200 -->,
  and the second-order error is that inconsistency. With the grid's own mean
  subtracted the pin changes `k` by round-off and the error is 1.9e-9
  <!-- claim: closure_rel_err_n200 -->, at fourth order.
- 2026-10-31: stated that the direct simulation converges to Eq. 25 at first
  order in Δz and that its error is upwind numerical dispersion. Over the window
  τ = 0.2 to 1.0 it measured, it converges to 0.416 %
  <!-- claim: direct_richardson_rel_err --> below Eq. 25, the start-up deficit
  inside that window; over τ = 0.6 to 1.2 the limit is -2.4e-5
  <!-- claim: direct_late_richardson_rel_err -->.
- 2026-10-31: stated that a wrong geometry (`nu = 0`) and a dropped mean
  subtraction move the closure error by 3 and 2 orders of magnitude; from the
  1.0e-4 then published they move it to 0.47 and 4.0, 3.7 and 4.6 orders.
- 2026-10-31: stated a "worst relative deviation over Pe = 2 to 398" of
  0.016 %, as if the closure's agreement varied with the Peclet number. The
  closure is linear in the velocity, so `k ∝ Pe²` exactly: 0.016 % was the
  error at `n_r = 160`, the same at every Pe to round-off.
- 2026-10-31: stated that every boundary of the 2-D simulation is no-flux and
  that mass conservation is exact by construction. The zero-gradient end
  conditions let upwind carry `u' C` through the two end faces; mass was
  conserved because the slug never reached them. The ends now carry no flux,
  which moves no published number.
