# A4.5: Fuller-Schettler-Giddings binary gas diffusivities

Fuller, E. N., Schettler, P. D. & Giddings, J. C., *"A New Method for Prediction
of Binary Gas-Phase Diffusion Coefficients"*, **Ind. Eng. Chem.** **58**(5)
18-27 (1966), [doi:10.1021/ie50677a007](https://doi.org/10.1021/ie50677a007).

The most-used empirical estimate of binary gas-phase diffusivity, rebuilt from
its own printed tables, and the difference between the number it is quoted for
and the number it earns.

- **Structures:** `S3` (scalar diffusion)
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (eq. (4) with every published constant in one dictionary, the
  least-squares problem of eq. (3), and the heated tube solved with pymrm and
  checked against a quadrature), tested by `test_model.py`.
- **Data:** two CSVs in `data/`, each with a `.meta.yaml` sidecar: Table I, the
  fitted diffusion volumes (tier 6, fitted constants), and Table III in full,
  340 <!-- claim: n_rows --> measured binary diffusivities from
  153 <!-- claim: n_systems --> systems, each with its literature reference and
  the FSG error the authors printed beside it (tier 5, measurements, every one
  inside the correlation's own fitting set). Read from the page images at their
  native 300 ppi; 23 <!-- claim: n_rows_image_corrected --> rows were corrected
  against the image where the text layer had altered a digit.
- **The 1966 erratum, `10.1021/ie50680a601`, is not on disk and was not
  consulted.** The 1966 table is used as printed.

## What the page shows

**The quoted average error is a goodness of fit, not a validation.** The
diffusion volumes were fitted by least squares to the same measurements the
average error is reported against, and the authors say so: *"Despite the large
sample, this is not a totally objective test since the parameters were obtained
from these particular data."*

**A held-out test for half the method.** The atomic-increment scheme is
additive, so an organic vapour can be deleted from the fitting set and still be
predicted from the increments the other compounds fix. Leave-one-compound-out
for 42 <!-- claim: loco_n_compounds --> organic vapours
(153 <!-- claim: loco_n_points --> points), refitting every volume each time:

| on the organic-vapour points | mean \|error\| |
|---|---|
| published volumes, in sample (Table II prints `5.0`) | 5.03 <!-- claim: loco_in_sample_published_pct --> % |
| own full refit, in sample | 5.04 <!-- claim: loco_in_sample_refit_pct --> % |
| **leave-one-compound-out, held out** | **6.61 <!-- claim: loco_held_out_pct --> %** |

The excess is 1.56 <!-- claim: loco_gap_pp --> percentage points with a standard
error of 0.48 <!-- claim: loco_gap_se_pp --> over the compounds, clustered by
compound (n = 42 <!-- claim: loco_n_compounds --> compounds; t above
2 <!-- claim: loco_gap_t -->). In 2 <!-- claim: loco_fits_at_bound --> of the
refits (aniline or benzene held out) the nitrogen increment sits on its lower
bound. The simple molecules each carry their own fitted volume, so their points
admit no held-out test, and the held-out figure must not be transferred to them.

## Agreement

Every check below is on the transcription and on the implementation of eq. (4),
not on the physics: the printed errors are the authors' own arithmetic on the
same equation.

| Check | Result |
|---|---|
| printed FSG errors recomputed | 338 <!-- claim: n_rows_reproduced --> of 340 <!-- claim: n_rows --> rows to 0.043 <!-- claim: max_abs_delta_pp --> pp, median 0.0038 <!-- claim: median_abs_delta_pp --> pp; 2 <!-- claim: n_rows_irreproducible --> rows are not reproducible (He-H₂O at `352.5` K, CO₂-ethylene oxide at `298.0` K), both confirmed as printed on the page image |
| route B: the `T^1.75` identity within a series, no volumes and no eq. (4) | 142 <!-- claim: routeB_n_ratios --> ratios in 51 <!-- claim: routeB_n_series --> series, every deviation below 1e-3 <!-- claim: routeB_max_dev_excl_anomaly --> (the printed rounding) except the 1 <!-- claim: routeB_n_beyond_rounding --> ratio through the He-H₂O row; exponent recovered from the printed table 1.7492 <!-- claim: routeB_recovered_exponent --> with every ratio, or 1.7500 <!-- claim: routeB_recovered_exponent_excl --> with the ratios through the He-H₂O row left out |
| printed aggregates | `AVE. ERROR` 4.3211 <!-- claim: ave_error_pct --> against `4.32`; the tail of errors beyond ten per cent 7.35 <!-- claim: pct_over_10 --> against `7.4`; `STANDARD DEVIATION` 6.7101 <!-- claim: residual_sd_pct --> against `6.71`, as the residual standard deviation with the 28 printed parameters (27 rows of Table I plus b); 10 <!-- claim: aggregates_within_print --> of the 11 within half a unit of the printed digit, the Table II N₂ mean at 1.023 <!-- claim: aggregate_worst_halfunits --> half-units, past half a unit by less than the 0.1 <!-- claim: aggregate_input_rounding_halfunits --> half-units by which rounding of the printed errors can move a mean |
| eq. (3), `phi'` | printed error column 1.4048 <!-- claim: phi_prime_printed_column -->, printed volumes 1.3996 <!-- claim: phi_prime_published_volumes -->, free refit 1.3701 <!-- claim: phi_prime_free_refit -->: all below the paper's constrained minimum `1.4762` |
| free refit | simple-molecule volumes within 1.8 % <!-- claim: refit_simple_volume_max_rel_dev --> of the printed ones; no volume on a bound; random starts agree below 1e-10 <!-- claim: refit_start_spread --> |
| pymrm tube against quadrature | 1.79e-7 <!-- claim: flux_pymrm_vs_quadrature --> at 400 cells; observed orders 1.84 <!-- claim: flux_grid_order_min --> to 1.98 <!-- claim: flux_grid_order_max --> on 25 to 800 cells; Richardson limit below 1e-8 <!-- claim: flux_richardson_vs_quadrature --> |

## Findings about the paper

**The printed `phi'` minima do not belong to the published table.** Eq. (4) is a
special case of eq. (2), so a minimum of eq. (2) cannot lie above a value
reachable inside eq. (4). The printed error column gives a `phi'` that sits
0.071 <!-- claim: phi_prime_gap_printed_min --> below the paper's constrained
minimum, far more than the rounding of the printed errors can move it; the
printed volumes give one 0.077 <!-- claim: phi_prime_gap_published --> below it.
Why is not resolvable from what is printed.

**The temperature exponent.** The paper says observed exponents *"usually lie in
the range 1.6 to 1.8"*. Of 26 <!-- claim: exponent_n_series --> series with four
or more temperatures, 20 <!-- claim: exponent_in_band --> slopes fall in the band,
median 1.704 <!-- claim: exponent_median -->, but only
6 <!-- claim: exponent_certain_in --> lie inside by more than two standard errors of the least-squares slope
and 4 <!-- claim: exponent_certain_out --> outside; the rest are too short to
place. The slopes run from 1.582 <!-- claim: exponent_min --> (H₂-NH₃) to
1.951 <!-- claim: exponent_max --> (H₂-SF₆), 7.9 <!-- claim: exponent_extremes_z -->
combined standard errors apart and 14.2 <!-- claim: exponent_spread_over_ci -->
times the full width of the paper's interval on b: that interval is the
precision of one global exponent, not the uncertainty for a gas pair.

**Two inferences from the printed error column.** The ring term applies to
aromatic rings only (the named CO₂-ethylene oxide row aside; on the non-aromatic heterocycles it leaves rows
8.7 <!-- claim: ring_delta_min_pp --> to 21.9 <!-- claim: ring_delta_max_pp -->
percentage points off), and M(D₂) is `4.00`. The
22 <!-- claim: n_rows_inferred --> rows they rest on are not independent tests; of
the remaining rows 1 <!-- claim: n_free_rows_failing --> is irreproducible.

**A printed inconsistency.** The `AVE. ERROR` and `STANDARD DEVIATION` rows at the
foot of Table III follow Table II's method order, not Table III's own column
order. Only the FSG column is transcribed here; anyone using the other eight
should settle the ordering first.

## What pymrm adds

Nothing to the correlation: eq. (4), the reproduction, the cross-validation and
the exponents are numpy and scipy. pymrm supplies the step a correlation cannot
take, turning an error in D into an error in an answer: a linear
variable-coefficient diffusion solve in a tube from 300 K to 600 K, one sparse
solve with the flux read at the boundary. Moving b across the measured range
changes the benzene flux by -5.85 <!-- claim: flux_change_b_min_anchor_cold --> %
to +7.33 <!-- claim: flux_change_b_max_anchor_cold --> % when D is pinned at the
cold end, and by -0.22 <!-- claim: flux_change_b_min_anchor_mid --> % to
+0.12 <!-- claim: flux_change_b_max_anchor_mid --> % when it is pinned at the
geometric mean temperature. On this tube the quadrature gives the same changes
to within 1e-4 <!-- claim: flux_change_routes_gap_pp --> pp; the solve is what
carries over to a geometry or source term with no quadrature. A scale error in D needs no solve: the problem is
linear in it.

## Scope

`A4.6` (Chapman-Enskog) is the first-principles route to the same quantity and
should own the comparison of the two on one measured axis; no number here comes
from Chapman & Cowling. Molar masses are modern IUPAC values; the paper prints
none.

## Corrections

- 2026-11-16: said "with $M = 4.00$ all twelve reproduce". Table III holds
  14 <!-- claim: n_rows_d2 --> D₂ rows, and all of them reproduce within the
  0.043 <!-- claim: max_abs_delta_pp --> pp of every reproducible row.
- 2026-11-16: said "every one of the fourteen D₂ rows sits about 0.3 percentage points below its printed error"
  (with `M(D₂) = 4.028`). The D₂ rows then sit
  -0.34 <!-- claim: d2_isotopic_delta_min_pp --> to
  -0.11 <!-- claim: d2_isotopic_delta_max_pp --> percentage points from their
  printed errors.
- 2026-11-16: said "rows are free tests, and 2 of them" fail. The CO₂-ethylene
  oxide row is among the 22 <!-- claim: n_rows_inferred --> rows the inferences
  rest on, so 1 <!-- claim: n_free_rows_failing --> free row fails.
- 2026-11-16: said withholding the ring term from the non-aromatic heterocycles "reproduces every one of them". The CO₂-ethylene oxide row at `298.0` K is irreproducible either way; the other 7 <!-- claim: n_hetero_reproduced --> reproduce.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: both datasets are committed CSVs in `data/`.
