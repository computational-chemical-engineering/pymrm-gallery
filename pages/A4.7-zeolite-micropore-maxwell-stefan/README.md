# A4.7: Maxwell–Stefan diffusion in zeolite micropores

Reproduces the single-component core of Krishna, R. & Baur, R. (2003),
*Modelling issues in zeolite based separation processes*, Separation and
Purification Technology **33**(3) 213–254
([doi:10.1016/S1383-5866(03)00008-X](https://doi.org/10.1016/S1383-5866(03)00008-X)):
the Maxwell–Stefan flux law for a sorbate against a stationary framework, the
thermodynamic correction factor `D = eth * Gamma`, the two loading scenarios the
review distinguishes, and the transient uptake problem it solves in a spherical
crystallite.

## Files

- `index.py` / `index.ipynb`: the page (paired by jupytext).
- `model.py`: `Isotherm` (the dual-site Langmuir of eq. 32, its exact inverse,
  eq. 33 and the roots that locate its extrema), `karger` (eqs. 23 and 24), the
  series eq. (30) with its tail bound, eq. (31), and `Crystal`, the spherical
  crystallite with a loading-dependent diffusivity solved with pymrm.
- `test_model.py`: the Sherwood number at three-quarters uptake against the stored value, the
  isotherm inverse and the single-site collapse, a constant-D crystal against
  eq. (30), and one loading-dependent Newton run.
- `data/krishna-baur-2003-sorption-parameters.csv`: the review's printed
  *inputs*: the DSL isotherm parameters of Tables 1 and 2, the benzene
  parameters printed inside the caption of Fig. 9, the jump frequencies and MFI
  cell dimensions behind eq. (23), and the two diffusivities of Section 4.1.
- `data/krishna-baur-2003-printed-results.csv`: the review's printed
  *outputs*, the targets the page is scored against, including **both** of the
  two values it prints for one orientation-averaged diffusivity, each with the
  string as printed.

## The result

**Tier 6: reproduction, not validation.** Every target on this page is a number
Krishna & Baur themselves computed; none is a measurement. Every experimental
comparison inside this page's scope lives only in a figure (Garg & Ruthven's
ethane uptake, Shah's benzene diffusivities, Millot's Arrhenius data), and
nothing here rests on one.

The printed targets reproduce to 0.78 % <!-- claim: targets_worst_rel_dev --> at
worst (the review's LDF constant "corresponding to a 75% approach to
equilibrium", against a computed 75.59 % <!-- claim: F_at_Sh10 -->). The DSL isotherm gives the review's two worked
3MP loadings within half a unit of their printed last digit (at worst
0.41 <!-- claim: dsl_3MP_halfunits_max --> of it); the permeation flux ratio lies
-0.17 % <!-- claim: flux_ratio_rel_dev --> from the printed fluxes' ratio, and
dropping the thermodynamic factor moves it by a factor of
3.73 <!-- claim: flux_ratio_gamma_lever -->. One further printed number, the
running text's orientation-averaged diffusivity, is out by
4.2 % <!-- claim: karger_eq23_vs_text_rel --> and is reported as a printed
inconsistency rather than scored.

**A factor of two in eq. (23).** The review prints two values of one
orientation-averaged zero-loading diffusivity, 4.1 % <!-- claim: karger_table3_vs_text_rel -->
apart. Recomputing Kärger's eq. (23), x and y reproduce to
0.29 % <!-- claim: karger_xy_max_rel_dev -->, and the printed D_z is
0.50005 <!-- claim: karger_z_printed_over_computed --> of the computed one: half the
computed value lies within half a unit of the printed last digit. Table 3 is
the value consistent with eq. (23) as printed
(0.070 % <!-- claim: karger_table3_rel_dev -->), the running text the one
consistent with its own printed D_z (-0.013 % <!-- claim: karger_text_rel_dev -->).
The review applies eq. (23) a second time, to CH4, and prints only the average:
the halved relation reproduces it to -0.25 % <!-- claim: karger_CH4_halved_rel_dev -->
and the relation as printed misses by +5.0 % <!-- claim: karger_CH4_unhalved_rel_dev -->.
So the review's practice is the halved relation in both worked examples and
Table 3 is the odd one out. Where the factor originates, in Kärger's 1973
original or in the review's arithmetic, the page does not say: Kärger is not on
disk.

**The claim under test.** "The M–S diffusivity is the loading-independent one"
holds by a factor of 10 <!-- claim: swing_weak_advantage --> in one of the
review's two scenarios and fails by the same factor in the other. The factor is
algebraically forced (Γ = 1/(1−θ) and eth = 1−θ are reciprocals) and its size is
set by the loading range chosen. The direction is the finding; what survives
both scenarios is that the pole in D at saturation lives entirely in Γ.

**The uptake problem.** With the non-linear eq. (18) the time-averaged
Sherwood number at three-quarters uptake runs from 12.82 <!-- claim: sh75_weak_03_ext --> at
θ_s = 0.3 to 35.98 <!-- claim: sh75_weak_09_ext --> at θ_s = 0.9, extrapolated in step and cell
size, against the series' 10.078 <!-- claim: sh75_series_exact -->; the production curves read
them low by at most 0.5 <!-- claim: sh75_production_rel_err --> %.
In the review's 3MP example adsorption is nearly complete at the Fourier number
the review names while desorption is not, under both loading scenarios (smallest gap in F
0.114 <!-- claim: mp3_asymmetry_gap_min -->).

## Validation

- Eq. (33) against a central difference of ln P in ln Θ: worst gap
  2.5e-6 <!-- claim: gamma_eq33_vs_numerical -->, inside the difference's own
  truncation and round-off budget; it fails when a site factor is dropped.
- The crystal against the exact series eq. (30), cell size and step refined on
  separate ladders: observed order 0.99 <!-- claim: dt_observed_order --> in the
  step and 1.95 <!-- claim: space_observed_order --> in cell size; the spatial error
  of the finest grid 8.6e-6 <!-- claim: space_refined_rel_err --> relative, a figure that
  depends on the order used to extrapolate; the
  production curves carry 1.54e-3 <!-- claim: production_rel_err --> of
  time-integration error, the resolution limit used for eq. (31).
- The same ladders on the 3MP desorption: the production value is
  1.6e-3 <!-- claim: nl_production_err --> from the value extrapolated in both,
  against the asymmetry gap above.
- Every solve is checked: `newton` reports success on every implicit step.
- A defect-injection table: every claim that is not a named identity fails at
  least one row (cell dimensions, the printed eq. 23c, transposed printed digits,
  isotherm rows, factors of eq. 33, the difference step, the loading range, the
  loading dependence of D, the geometry index, the crystal resolution, the
  step schedule, the ladders, the prefactor and length of the series, the
  membrane stoichiometry).

## What pymrm adds

pymrm supplies the crystallite: `construct_grad` and `construct_div` for the
sphere, `construct_boundary_value_matrices` for the boundary-face loadings,
`newton` per backward-Euler step for the loading-dependent diffusivity, and from
`pymrm.checks` `find_roots` (extrema of Γ, roots of eqs. 30 and 31, crossings on
the solved model) and `observed_orders` (the refinement ladders). The algebra of
eqs. (23), (32) and (33) and the two-extrema survey are the page's own.

**One headline has no printed target behind it.** The "two extrema in Γ" survey
rests on hand-transcribed Table 1 rows; only C1 and nC4 are certified by
anything, and the eq. (33) identity check is blind to a wrong parameter. The
page measures the exposure instead: of the injected single-character slips,
1 <!-- claim: extrema_survey_verdicts_moved --> rewrites the survey's count, and
another moves the 2MB minimum from 5.71 <!-- claim: extrema_2MB_thmin --> to
4.19 <!-- claim: extrema_2MB_thmin_slipped --> with every verdict unchanged.

**What the page does not reproduce.** The review's Fig. 10(a) reports an
inflection in the 3MP desorption kinetics; under neither stated loading
dependence does the desorption rate rise at any step
(0 <!-- claim: mp3_desorption_rate_rises --> steps). The claim lives in an
undigitised figure and the page does not resolve it.

**Three printed defects are recorded, none silently fixed:** the sign of
eq. (25); `nu_zz + nu_zz` in the denominator of eq. (23c); and the symbol
`nu_str` printed twice in the sentence supplying the CH4 jump frequencies, where
the second must be `nu_zz`. The third is a repair by inference and is flagged as
one in the CSV row's note, in both sidecars and on the page.

## Scope

Single component only. Mixture diffusion, the exchange coefficients, IAST and
the membrane selectivities are case `H1.9`; the argument for cutting there rather
than at micropore/membrane is on the page.

**Source text.** The Elsevier full-text API, PII `S138358660300008X`, with the
authorised institutional key. The 2003 review is born-digital, so the
"API text is unusable for numbers" warning in `docs/pdf-findings.md`, which is
about pre-1980 scans, does not apply. No PDF or page image is stored.

**The catalogue's other reference is not this paper.** `A4.7` was previously
pointed at `10.1016/0009-2509(95)00102-B`, which is van den Broeke & Krishna
(note the author order; the catalogue reverses it), a packed-bed breakthrough
study on activated carbon and a carbon molecular sieve with no zeolite
experiment in it. See `docs/source-sweep-2026-08-02.md`.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```
