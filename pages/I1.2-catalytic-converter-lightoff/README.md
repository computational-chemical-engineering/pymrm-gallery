# I1.2 — Light-off in a monolithic catalytic converter

An ignition front sweeping down a monolith channel, and the wrong-way
temperature rise that appears when the same converter is suddenly cooled.

- **Structures:** `S4` (transient 1-D convection and reaction), `S7` (two-phase
  gas/solid with interfacial transfer)
- **Reference:** Oh, S. H.; Cavendish, J. C. (1982), *Ind. Eng. Chem. Prod. Res.
  Dev.* 21(1) 29–37, doi:10.1021/i300005a006
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`OhCavendish` for Table I and the kinetics, `Feed` for a case, and
  `Monolith`, the twelve-field channel on one grid with constant operators,
  `NumJac` and one bounded `newton` call per backward-Euler step; `march` with
  root-found events, and `peak`, which removes the sampling error of a maximum
  in time), tested by
  `test_model.py`.
- **Runtime:** about 7 minutes, half of it the defect-injection table, which
  re-runs every transient result with one defect in place.
- `data/`: Table I (parameters), Table II (cell density) and the printed case-2
  results, each with its `.meta.yaml` provenance sidecar.

## Agreement

**All six numbers of Table III, with nothing fitted.** 0.63 <!-- claim: table3_T_mean_abs_dev_pct --> %
mean absolute deviation in the peak wall temperature (990.4 <!-- claim: table3_T_model_1.675e-03 -->,
986.1 <!-- claim: table3_T_model_1.675e-02 --> and 966.9 <!-- claim: table3_T_model_1.675e-01 --> K
against 996 <!-- claim: paper_T_s_max_1.675e-03 -->, 991 <!-- claim: paper_T_s_max_1.675e-02 -->
and 975 <!-- claim: paper_T_s_max_1.675e-01 --> K), and at worst
1.25 <!-- claim: table3_x_max_dev_cells --> cells of the paper's own grid in its
location. Both trends reproduced: raising the solid conductivity lowers the peak
and moves it downstream. The transient rise is 236 <!-- claim: temperature_rise_K --> °C
against their 241 <!-- claim: paper_temperature_rise_K -->, and the peak is
reached at 9.73 <!-- claim: table3_t_peak_1.675e-02 --> s with a 0.05 s step
and at 9.76 <!-- claim: t_peak_dt_richardson_s --> s extrapolated in the step
(first-order Richardson), 3.4 <!-- claim: t_peak_dt_richardson_dev_pct --> %
before their 10.1 <!-- claim: paper_t_peak_s --> s. The wall at the inlet ends
97.9 <!-- claim: overshoot_K --> K above the adiabatic temperature on 80 cells;
that value is first order in the grid, and its Richardson estimate is about
108 <!-- claim: overshoot_richardson_K --> K.

Free checks that test different things: Table II reproduces from the cell
density and the stated wall thickness alone to 0.09 <!-- claim: table2_geometry_max_dev_pct --> %;
the Lewis number of H₂ comes out 3.89 <!-- claim: lewis_H2 --> against the
about 4 <!-- claim: paper_lewis_H2 --> the paper states, which tests the
diffusivities, the gas conductivity correlation and the assumed molecular weight
together. The species balance closes below 1e-7 <!-- claim: species_balance_closure -->
of the smallest species flux fed at every time level, and the energy balance below
1e-10 <!-- claim: energy_balance_net_closure --> of the heat released.

Light-off (a prediction, not a comparison; see below): 47.12 <!-- claim: lightoff_CO_s --> s
to 50 % CO conversion, root-found on the step that crosses it, moving by
0.45 <!-- claim: lightoff_grid_dt_spread_s --> s over a fourfold change in grid
and a fourfold change in time step.

**The Table III residual is a bias, not scatter, and it is not the grid.** All
three temperatures are low by roughly the same fraction; the value changes by
1.5 <!-- claim: case2_change_80_320_K --> K over a fourfold grid refinement while
the gap to the paper is 4.9 <!-- claim: table3_gap_to_paper_std_K --> K. Something
small and systematic separates this reimplementation from theirs, most plausibly
one of the four unprinted quantities below. The page does not claim to know which.

## Provenance: tier 6, not experimental

The paper contains **no experimental data of any kind**. Its two experimental
references, Hegedus (1975) for the steady-state overshoot above adiabatic and
Mondt (1981) for the absence of melting during a deceleration, are qualitative
appeals to other people's work with no numbers attached. Everything checked here
is the authors' own computed result, from a Galerkin discretisation integrated
with GEAR. Do not describe this page as validated against experiment.

The asymmetry between the two cases is the thing an integrator should know:
**case 2 (step decrease) has the published numbers; case 1 (light-off) has
none.** The light-off time on the page is a prediction, and the only claims made
about it are that it is converged and that the qualitative sequence the paper
describes in words comes out.

## Four things the paper needs and does not print

Each is taken from a source the paper itself names, and each is checked on the
page rather than assumed:

1. **`S`**, geometric area per unit volume, follows from the channel geometry,
   `S = 2*eps/R_h`. Checked by reproducing all six numbers of Table II from the
   cell density and the stated wall thickness.
2. **`D_i`**, the Slattery and Bird formula the paper cites (Bird et al. 1960,
   p. 505). Checked against the paper's own statement that H₂ has a Lewis number
   near 4.
3. **`Nu_inf`, `Sh_inf`**: the paper says only "constant wall heat flux, Shah
   and London (1971)", which for a square duct is three different numbers (H1,
   H2 and T). Table III is used as the arbiter; the page prints what each
   alternative gives.
4. **`(-dH)_i`**, standard formation enthalpies, water as vapour.

## Two readings of the printed equations, and how they were settled

Equation 12 is printed as `d(C_ps*T_s)/dt`, not `C_ps*dT_s/dt`, and the
Nomenclature confirms it with `Psi = (1-eps)*rho_s*[C_ps + T_s*dC_ps/dT_s]`.
With the tabulated `C_ps` the two readings differ by a factor of
2.1 <!-- claim: inertia_ratio_300K --> in effective thermal inertia at the start
of case 1. The page runs both and lets Table III choose; the time of the peak
discriminates much more sharply than the peak temperature.

`R_h` is defined in the paper as `2*(area/perimeter)`, that is **half** the usual
hydraulic diameter, which is why eqs 17 and 18 divide by `2*R_h`. Table II
confirms the reading.

## pymrm points recorded here

1. **Bounds keep Newton on the physical root.** The inhibition denominator of the
   kinetics gives the residual roots with negative mole fractions; without bounds,
   full Newton steps take the first step of case 2 to one of them and report
   convergence. `bounds=` with every mole fraction in `[0, y_in]` and the
   temperatures in a range where the rates are defined is the only guard the
   solve needs, with `damping="monotone"`.
2. **Check every step.** With the old-time accumulation term given the wrong
   sign, `newton` stops on the temperature bound with status 8 ("stopped on a
   bound"); a march that does not check `success` reads that as a result. Here a
   step that does not converge raises `StepFailed`.
3. **A first-order outlet is `{"outflow": True}`.** The gas equations have no
   outlet condition; a zero-gradient dictionary imposes one and extrapolates the
   face. With the outflow marker the face is the last cell and both balances close
   to the Newton tolerance.
4. **A refined extreme is still a fixed-step value.** Refining the peak in time
   over the length of one step removes the sampling error of up to half a step,
   but backward Euler's own first-order error has the same sign; extrapolate in
   the step before comparing with a reference. The light-off root differs from
   linear interpolation between levels by milliseconds, far inside the
   refinement spread.
5. **Nothing that touches a DataFrame belongs inside the residual.** Resolve the
   parameter table to plain floats before the time loop (`OhCavendish.from_table`).

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/I1.2-catalytic-converter-lightoff/index.py   # after editing index.py
python scripts/run_kernels.py pages/I1.2-catalytic-converter-lightoff
GALLERY_OFFLINE=1 PYTHONPATH=shared python -m pytest -q pages/I1.2-catalytic-converter-lightoff/test_model.py
```

No PDF is needed to rebuild: all three datasets are committed CSVs.

## Corrections

- 2026-10-30: stated that the species balance closes to 3.7e-8 and the energy
  balance to 1.3e-11 of the heat released; the run gives
  5.3e-9 <!-- claim: species_balance_closure --> and
  3.0e-12 <!-- claim: energy_balance_abs_closure -->, both set by the Newton
  tolerance.
- 2026-10-30: stated that a species balance written on the last cell instead of
  the outlet face closes to 1 % instead of 1e-12. With the zero-gradient outlet
  the old page used, the last-cell balance misses by
  0.8 <!-- claim: zero_gradient_species_closure_cell --> of the smallest species
  flux fed (face and last cell differ by up to
  3e-3 <!-- claim: zero_gradient_face_cell_gap --> of the feed), and the face
  balance closes below 1e-7 <!-- claim: species_balance_closure -->, not to 1e-12.
- 2026-10-30: stated that without the two-sided scaling the first sparse solve
  returns "matrix is exactly singular", and that an undamped Newton step from the
  initial guess overshoots into negative temperature. With the bounds, unscaled
  and undamped steps both converge; what the bounds prevent is a converged root
  with negative mole fractions, not a negative temperature.
- 2026-10-30: stated that with the zero-gradient outlet the face value differs
  from the last cell centre by 3 parts in 10^4. Over the first minute of case 1
  the largest difference is 4.2e-5 <!-- claim: zero_gradient_face_cell_gap_abs -->
  in a mole fraction, which is 3e-3 <!-- claim: zero_gradient_face_cell_gap -->
  of that species' feed.
- 2026-10-30: stated that the Jacobian is not approximated, the exact pointwise
  block Jacobian; `NumJac` builds the pointwise blocks by finite differences (the
  constant operators are analytic), so the sparse structure is complete but the
  pointwise entries are differenced.
- 2026-10-30: stated that the twelve fields span eleven orders of magnitude; they
  span 7 <!-- claim: state_span_decades --> decades.
- 2026-10-30: `meta.yaml` called the refinement study fifteen runs; it is eleven.
- 2026-10-30: this README said the page was in progress and staged in
  `queue_cases/`; it is published.
