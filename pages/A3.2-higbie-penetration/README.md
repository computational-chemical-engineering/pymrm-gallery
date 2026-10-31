# A3.2: Higbie penetration theory

A liquid surface is exposed to a gas for a definite length of time, absorbs by
unsteady diffusion the whole while, and is then swept away. Nothing reaches a
steady state and no film thickness appears anywhere, so the coefficient is
`k_L = 2 sqrt(D / pi t_exp)` and `k_L ~ D^(1/2)` against film theory's `k_L ~ D`.

## Where the model comes from

R. Higbie, *Trans. AIChE* **31**, 365-389 (1935) is the origin. It is pre-DOI,
it is not on disk, and **it was not consulted**. Everything here is read from
**Bird, Stewart & Lightfoot, *Transport Phenomena*, 2nd edn (Wiley, 2002),
section 18.5** (book pp. 558-561), which derives the result, names the
"penetration model", attributes it to Higbie in the same paragraph and prints the
full citation in a footnote that also credits him with having "provided the basis
for the 'penetration model' of mass transfer". Every equation number on the
page is BSL's, and every digit was read on 300 dpi renders at the scan's native
resolution; its ABBYY text layer destroys equations and was used only to find
pages.

That one sentence and that footnote are **all** BSL prints about the 1935 paper,
so nothing on the page is attributed to it beyond what BSL prints and, where the
Froment cross-check is described, what Froment prints about it (his "Higbie's
uniform age"). The two geometries the contact time is read off, the falling film
(`t_exp = L/v_max`) and the rising bubble (`t_exp = D/v_t`), are BSL's section
18.5 and BSL's Example 18.5-1; the page does not claim they are Higbie's.

**Froment, De Wilde & Bischoff, 3rd edn (Wiley, 2011), section 6.4** is a
cross-check only: an independent statement of the same constant, the closed-form
enhancement factor for Higbie's uniform surface age, and Table 6.4.2.1.
Levenspiel names and cites Higbie but neither derives nor tests him, and is not
used.

## How it is solved

`Film` (`model.py`) marches BSL's Eq. 18.5-7 (parabolic velocity profile, no-flux
wall) or Eq. 18.5-11 (plug flow, deep domain) in the flow direction. The
transverse operator comes from `construct_grad` and `construct_div`, assembled
once; the interfacial flux is read with `compute_boundary_values`. Each step is
one banded solve (scipy) on steps uniform in `sqrt(zeta)`, by backward Euler with
Richardson extrapolation in the step. The page measures that choice against BDF2
and Crank-Nicolson at equal cost, and locates the 1 % threshold with
`pymrm.checks.find_roots`; refinement studies use `pymrm.checks.observed_orders`.

## What the page establishes

- **Problem 18A.4 reproduced two ways, each scoped to what it tests.**
  Eq. 18.5-18 gives 0.27351 <!-- claim: wa_closed_gmol_per_hr --> g-mol/hr against the printed 0.273
  (+0.188 % <!-- claim: eq18518_vs_printed_18A4_rel -->); a pymrm march of
  Eq. 18.5-11, read as the outlet mixing cup, gives the same to the march's grid
  error. The two routes are *not* independent end to end: the wetted width, the
  3/2, the solubility conversion and the units enter both and cancel out of the
  ratio, so the cross-route gap tests the dimensionless `Phi(Lambda)` against
  `sqrt(4 Lambda / pi)` and nothing else. The +0.188 % against the printed
  answer is what tests the prefactor. It is a hair *outside* three-figure
  rounding; two unprinted inputs (the density of water, the molecular weight of
  Cl2) are each worth about that much, and so is half a unit in the last printed
  digit of the diffusivity. Neither unprinted value is adopted.
- **A printed defect in the pinned source.** Problem 18A.7(a) says "Use
  Eq. 18.5-20" and prints Eq. 18.5-19's answer. The two differ by exactly
  sqrt(3). Reported, not repaired.
- **The one empirical test on disk.** For the same bubble the penetration
  prediction carries no adjustable constant, `t_exp = D/v_t`, and comes out
  **-12.0 % <!-- claim: higbie_kL_vs_hg_measured_rel -->** from Hammerton &
  Garner's measured `k_c`. Creeping flow would be
  -49.2 % <!-- claim: creeping_kL_vs_hg_measured_rel --> from it.
- **The three-way comparison, taken.** Froment's Table 6.4.2.1 is the only place
  in either book with film, surface renewal and penetration side by side. Its
  surface-renewal column rounds correctly everywhere; its film column is high by
  one unit in the last printed place at 2 <!-- claim: froment_t6421_film_cells_off_rounding -->
  of six; **its penetration column disagrees with the formula printed at its own
  head at all 6 <!-- claim: froment_t6421_penetration_cells_off_rounding -->**,
  and two cells are impossible on the printed digits alone (0.94 where an
  enhancement factor cannot be below 1, and 10.39 above the same row's
  surface-renewal 10.05).
- **The penetration approximation, decomposed.** BSL justifies Eq. 18.5-11 with
  two arguments as though they were one. Solving Eq. 18.5-7, which BSL sets up
  and then declines to solve, shows they are worth a factor
  4810 <!-- claim: profile_over_wall_effect_at_lambda_0p1 --> apart: at
  Lambda = 0.1 the velocity profile costs
  1.92 % <!-- claim: profile_effect_at_lambda_0p1 --> and the finite film
  4.0e-6 <!-- claim: wall_effect_series_at_lambda_0p1 -->, from the exact
  series and confirmed by the solver refined to 4800 cells
  (4.0e-6 <!-- claim: wall_effect_solver_refined_at_lambda_0p1 -->). The 1 %
  threshold is Lambda* = 0.0564 <!-- claim: lambda_star_1pct -->, the
  penetration depth at which it is reached is 0.475 <!-- claim: pen_depth_over_delta_at_1pct -->
  of delta.

## What it does not do

It does **not** test the sqrt(D) exponent. One measured coefficient at one
diffusivity cannot resolve an exponent, and the page carries the search behind
that statement rather than asserting it. The exponent question stands where
[`A3.3`](../A3.3-danckwerts-surface-renewal/) left it: film theory and surface
renewal each fit a single datum exactly, one free constant each, zero residual
each. What is new is that on this datum penetration theory has no free constant at
all, so it is the only one of the three that can be wrong, and it is, by 12 %.

The measurement is quoted from BSL, not consulted; BSL prints no error bar for
it and none is invented. The film thickness of Problem 18A.4's column is not
printed and is not inferred.

## Corrections

- 2026-10-31: said the solver "cannot see the wall effect at all", because at its
  settings its discretisation error was of the same size (its estimate 1.04e-5
  against the exact 4.0e-6), and called that the strongest statement that the
  effect is negligible. That holds for the 600-cell grid only, whose grid error
  (6.4e-6 at Lambda = 0.1) exceeds the effect. Refined to 4800 cells (grid order
  2.00) with steps 6400/12800, the solver gives
  3.988e-6 <!-- claim: wall_effect_solver_refined_at_lambda_0p1 -->, within
  1.0e-4 <!-- claim: wall_effect_refined_vs_series_rel --> of the series; the
  remainder is the extrapolated step error, which falls about 16 times for four
  times the steps.
- 2026-10-31: called the worst deviation of the march from the exact wall series,
  7.2e-6 <!-- claim: wall_series_vs_pymrm_worst_rel -->, the size of the solver's
  marching error. It is the grid error: doubling the cells divides it by 3.9 at
  every Lambda; the extrapolated step error at the same steps is about 3e-8.
- 2026-10-31: said the marching step is the larger error, and that this holds
  a fortiori at the production settings. It holds for the raw march; every
  production value is extrapolated in the step, and then the grid error is the
  larger. On the deep plug-flow case (800 cells) the grid error is
  2.8e-6 <!-- claim: space_err_n800 --> against an extrapolated step error of
  3.3e-8 at steps 800/1600; on the production film (600 cells, steps 800/1600)
  the smallest ratio of grid error to step error over Lambda from 1e-4 to 0.2 is
  46 <!-- claim: production_grid_over_step_err_min -->.
- 2026-10-31: gave the grid order as 2.07, from errors against the finest grid
  taken as exact. From successive differences of three levels it is 2.00 (1.997
  to 1.999); the grid error at 400 cells is 1.12e-5
  <!-- claim: space_err_n400 -->, not 1.05e-5.
- 2026-10-31: called the 1.3e-12 gap between Froment's Eq. 6.4.2-8 and the
  quadrature of his Eq. 6.4.2-4 pure round-off. It was the default tolerance of
  `scipy.integrate.quad`; at relative tolerances of 1e-10 and 1e-13 the gap is
  about 2e-14 both times, which is round-off.
- 2026-10-31: labelled the mixing-cup identity (outlet inventory equals the
  integrated interfacial flux) structural. It is an identity of backward Euler on
  the conservative scheme, but a defect in the cup weight or in the flux reading
  breaks it, so it is checked as a bound with break rows. A break row "arithmetic
  instead of harmonic face mean" injected no defect (the diffusivity is uniform)
  and is gone.

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        Film (the finite-volume march), phi_extrapolated, phi_penetration,
                phi_wall_series
test_model.py   fast checks of model.py and two stored claims
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
data/           three CSVs of printed numbers, each with a provenance sidecar
```

Runtime about 35 s.
