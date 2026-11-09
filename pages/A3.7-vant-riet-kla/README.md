# A3.7: van 't Riet's k_L a correlation for stirred gas-liquid vessels

`k_L A = 2.6e-2 (P/V)^0.4 v_s^0.5` for water and
`k_L A = 2.0e-3 (P/V)^0.7 v_s^0.2` for strong ionic solutions, from
van 't Riet, K., *Ind. Eng. Chem. Process Des. Dev.* **18**(3) 357-364 (1979),
[doi:10.1021/i260071a001](https://doi.org/10.1021/i260071a001).

**It is a review, and the page is built accordingly.** Roughly ninety papers are
cited and five figures compile other people's measurements. Two things in it are
van 't Riet's own, eqs. 8 and 9, and they are what the page is about. Every
number the page reads carries an `origin` field saying whose result it is.

- **Structures:** `S3` (1-D steady convection with a linear sink), `S4` (1-D
  transient PDE coupled to a lumped phase)
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (the two correlations and their closed-form arithmetic, the gas
  column, the gassing-out experiment, eq. 2 and the Lagrangian route), tested by
  `test_model.py`. The studies, the claims and the break table stay in the
  notebook.
- **Runtime:** about two minutes, most of it the defect-injection table.

## Data

| file | what it is |
| --- | --- |
| `data/vant-riet-1979-correlations.csv` | eqs. 8 and 9: prefactors, exponents, validity windows, accuracy claim |
| `data/vant-riet-1979-printed-numbers.csv` | every other number the review prints as text, each tagged with **whose result it is** |
| `data/vant-riet-1979-figure-legends.csv` | the legend tables printed inside Figures 1-5: **printed text, not extracted curve data** |

Each CSV has a provenance sidecar with a `columns:` block. Read
`vant-riet-1979-figure-legends.meta.yaml` before using that file: it says
exactly what the legends are and are not.

## What the page shows

1. **Figure 5's topology, from the printed constants.** Pure-water curves cannot
   cross pure-water curves and ionic cannot cross ionic, so all nine candidate
   crossings are pure-vs-ionic, each a closed-form root. Exactly 1 <!-- claim: fig5_crossings_in_window --> falls inside the plotted window: the
   `4 cm/s` water curve meeting the `0.5 cm/s` ionic curve at P/V = 826.5 <!-- claim: fig5_crossing_pv --> W/m³, which is what the printed figure shows.
   No figure coordinate is measured. The window's left edge, 400 W/m³, is
   Figure 5's printed axis label `4.10²`, transcribed text; the count is the same
   for any left edge from 328.0 <!-- claim: fig5_next_crossing_pv --> W/m³ up to
   the crossing, and at the correlations' printed validity bound of 500 W/m³.
2. **Which constants the data exercise.** With the prefactor re-fitted at the
   window centre in every case, deleting a term costs:

   | term | cost of deleting it | inside the review's own 20-40 %? |
   | --- | --- | --- |
   | eq. 8 `(P/V)^0.4` | 82.1 % <!-- claim: null_pure_alpha0_cost --> | no, load-bearing |
   | eq. 8 `v_s^0.5` | 73.6 % <!-- claim: null_pure_beta0_cost --> | no, load-bearing |
   | eq. 9 `(P/V)^0.7` | 185.3 % <!-- claim: null_ion_alpha0_cost --> | no, load-bearing |
   | **eq. 9 `v_s^0.2`** | **26.4 % <!-- claim: null_ion_beta0_cost --> over Figure 4's drawn lines; 35.2 % <!-- claim: null_ion_beta0_cost_dataspan --> over its data span** | **yes, not resolved** |
   | the branching itself | 83.3 % <!-- claim: null_single_branch_worst --> under-prediction | no, the most load-bearing feature |

   The whole pure-water literature spread of exponents the review quotes costs
   at most 16.2 % <!-- claim: purewater_alpha_range_cost --> (alpha) and
   24.7 % <!-- claim: purewater_beta_range_cost --> (beta).
3. **An internal inconsistency in the review.** It expects ions to raise `k_L A`
   by 2-10×, argued from bubble mechanism and Zlokarnik and **not** from eqs. 8
   and 9. The ratio of the two fitted correlations spans 1.30 <!-- claim: ion_ratio_min --> to 5.98 <!-- claim: ion_ratio_max --> times over
   the printed box, its extremes at two corners fixed by the exponent signs, and
   lies inside 2-10× on 83.7 % <!-- claim: ion_ratio_frac_in_2_10 --> of the box
   (log-area share, in closed form). In the high-`v_s`, low-`P/V` corner it falls
   short of the 2× floor. Reported unresolved; what would settle it is in
   Figures 3 and 4, which are not digitised.
4. **The electrode ceiling closes.** The printed `k_L A = 0.1` s⁻¹ is Linek's
   `1/(5 τ_p)` at `τ_p = 2` s, 0.1000 <!-- claim: electrode_ceiling_from_linek -->
   s⁻¹. Over P/V 500-10 000 W/m³ and Figure 4's C1..C2 span of `v_s`, eq. 8 lies
   above it on 23.2 % <!-- claim: frac_window_above_ceiling_pure --> of the window
   and eq. 9 on 80.5 % <!-- claim: frac_window_above_ceiling_ionic --> (log-area shares, in closed form); the
   ionic branch has passed it by 1215 <!-- claim: pv_where_ionic_passes_ceiling -->
   W/m³ at every velocity Figure 5 draws.
5. **What the well-mixed-gas assumption costs, with pymrm.** Eq. 2, the method
   most of the data come from, assumes a uniform gas composition. The review
   warns about this four times and never computes anything. At the centre of the
   printed `P/V` window, at the lowest velocity Figure 5 draws, in a 2600 L
   vessel, the page's transient gassing-out experiment (two-axis limit of a step
   ladder and a cell ladder) shows eq. 2 under-reporting `k_L A` by 19.1 % <!-- claim: bias_pure_transient --> on the pure-water branch (the quasi-steady
   closed form gives 16.0 % <!-- claim: bias_pure_quasisteady -->, optimistic
   because it assumes `τ_G ≪ 1/k_L A`) and by 45.4 % <!-- claim: bias_ionic_transient --> on the ionic one at the `(0.2, 0.8)` fit
   window. On the ionic branch `τ_G` = 29.8 <!-- claim: tau_G_worked --> s
   against `1/(kφ)` = 12.0 <!-- claim: liquid_time_ionic --> s, the response is
   not a single exponential, and what eq. 2 returns depends on the window:
   34.1 % <!-- claim: bias_ionic_window_min --> to 56.5 % <!-- claim: bias_ionic_window_max --> across the three windows swept, with the
   quasi-steady 45.6 % <!-- claim: bias_ionic_quasisteady --> inside that range.
   Either way the bias is the same size as the entire accuracy the review claims.
   Over the vessel volumes its own legends print, the quasi-steady ionic bias runs
   from 6.5 % <!-- claim: bias_ionic_2p5L --> at `2.5 L` to 52.2 % <!-- claim: bias_ionic_5100L --> at 5100 L; across the printed aspect-ratio
   window the tall vessel is the worst case, 59.3 % <!-- claim: bias_ionic_TH05 -->
   at `T/H = 0.5` against 38.0 % <!-- claim: bias_ionic_TH15 --> at `T/H = 1.5`.

**This is a self-consistency statement, not a correction.** If the underlying
measurements were biased low by φ, eq. 9 as fitted already contains that bias.
What is established is that eq. 9 taken at face value predicts conditions under
which the method that produced most of its data does not measure what it is
assumed to measure.

## Agreement

| Check | Result |
| --- | --- |
| Figure 5's crossing count and pair, from the six printed constants | one crossing, pair `c x d`, as printed |
| electrode ceiling from Linek's criterion | 0.1000 <!-- claim: electrode_ceiling_from_linek --> s⁻¹ against the printed `0.1` |
| eq. 3's two printed forms | equal below 1e-12 <!-- claim: eq3_forms_gap --> for the vessel built |
| pymrm column against `(1 - exp(-N_G))/N_G` | 9.33e-5 <!-- claim: pymrm_vs_closed_n1600_pure --> (pure) and 2.00e-4 <!-- claim: pymrm_vs_closed_n1600_ionic --> (ionic) at 1600 cells, first order; the Richardson limit below 1e-6 <!-- claim: pymrm_vs_closed_richardson_ionic --> |
| transient, step ladder and cell ladder | backward Euler and upwind, both first order; the ionic `(0.2, 0.8)` fit does not depend on the cells, below 1e-8 <!-- claim: kapp_space_effect_ionic -->, because it closes before `τ_G` |
| transient against the Lagrangian route along the gas characteristics | the two limits agree below 1e-4 <!-- claim: lagrangian_limit_gap_pure --> on both branches |
| two-phase oxygen balance, faces read with `compute_boundary_values` | below 1e-8 <!-- claim: o2_balance_ionic --> |

Nothing is compared against a measurement: the data are in Figures 1-5, which
are not digitised, so no fit quality is computed for eqs. 8 and 9.

## What pymrm adds

To the correlation, nothing: eqs. 8 and 9 are five lines of arithmetic, and
Results 1 to 4 are `numpy`. To the measurement underneath it, a number. The
upwind convection operator with an outflow outlet, the divergence and the face
reads discretise the gas column exchanging with a well-mixed liquid; with them the
page gives the governing group `N_G = (1-ε) k H_t/(v_s H)`, which is not the
review's `τ_G · k_L A` (they differ by `(1-ε)/(ε H)`, 0.300 <!-- claim: criterion_factor --> here), and the bias in a measured `k_L A`
listed above. Starting the gas phase inert instead of pre-equilibrated lowers the
recovered ionic `k_L A` from 0.546 <!-- claim: kapp_over_k_ionic --> to 0.284 <!-- claim: kapp_over_k_ionic_inert_start --> of the true value, a number that
moves only from 0.280 <!-- claim: kapp_over_k_ionic_inert_start_w19 --> to 0.292 <!-- claim: kapp_over_k_ionic_inert_start_w37 --> across the fit windows, while the
pre-equilibrated one moves from 0.435 <!-- claim: kapp_over_k_ionic_w19 --> to
0.659 <!-- claim: kapp_over_k_ionic_w37 -->.

## What this page does NOT do

- **No figure is digitised.** Figures 1-5 hold the only measurements in the
  review; extracting them needs a maintainer overlay review and none is
  available. So the page is a **reproduction, not a validation**, and the data
  tier is 6.
- **Three inputs are not printed in the review** and all three are marked: the
  Henry constant `H` (the load-bearing one: φ_ionic runs from 0.424 <!-- claim: phi_ionic_H20 --> at `H = 20` to 0.625 <!-- claim: phi_ionic_H40 -->
  at `H = 40`), the gas holdup `ε` (which cancels except through `1-ε`) and the
  aspect ratio, taken at the centre of the review's printed `0.5 < T/H < 1.5`
  with both extremes swept at constant vessel volume, the implied volume printed
  beside the sweep.
- **It does not repair the two printed defects.** The sentence under eq. 7 says
  *"where N = stirrer diameter"*; the Nomenclature gives `N` = stirrer speed and
  `D` = stirrer diameter. And four of the five figure legends cite the drawn
  correlation by an equation number one lower than the text does (Figures 1-3
  print `CORRELATION EQUATION (7)` for eq. 8's line, Figure 4 prints
  `CORRELATION EQUATION (8)` for eq. 9's); Figure 5's legend cites none. The
  offset is a column of the legend dataset.
- **Eq. 8's validity statement prints no lower volume bound**, so the CSV cell is
  empty; the running text one sentence above eq. 8 says the compiled volumes
  "may range between 2 and 2600 L", and the sidecar records both facts.

## Corrections

- 2026-11-09: said that "every defect the docstring names moves the gap by more
  than two orders of" magnitude, where the docstring of the Lagrangian route
  claimed it "can catch a mis-assembled operator, a wrong velocity, a wrong sign or a wrong"
  (1-eps). No velocity defect was run, so the statement covered a case it never
  measured. The page now claims only what its rows measure: a flipped transfer
  sign and `ε` written where `1-ε` belongs each move the gap between the pymrm and
  Lagrangian limits past its bound (below 1e-4 <!-- claim: lagrangian_limit_gap_ionic -->), and under `nu = 1` the trace cannot be fitted, so the row records nan.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: the three datasets are committed CSVs in `data/`.
