# J3.5: Single particle models with electrolyte (SPM / SPMe)

Marquis et al.'s (2019) asymptotic reduction of the Doyle–Fuller–Newman model,
the SPM at leading order and the SPMe at first order, applied to the gallery's
own 1993 half-cell (`J3.4`), with the error scaling the asymptotics promise
measured against the full model.

- **Structures:** `S8` + `S10`, reduced towards `S3`
- **Reference:** Marquis, Sulzer, Timms, Please & Chapman (2019),
  *J. Electrochem. Soc.* **166**(15) A3693–A3706, doi:10.1149/2.0341915jes
  (open access)
- **Foundation:** the `J3.4` page: its `model.py` is copied unchanged into this
  page's `model.py` (the page checks that the copy reproduces J3.4's stored
  claims), and its maintainer-reviewed Doyle 1993 parameter set is loaded from it
- **Runtime:** about 3 minutes (most of it the DFN reference ladder)

## Provenance: tier 6, twice over

Every comparison on this page is reduced model against full model, and the
full model (`J3.4`) is itself a reproduction of a published *simulation*. No
measurement exists anywhere in this chain. Do **not** describe anything here
as validated against experiment.

## Agreement

| what | pymrm | paper / theory |
|---|---|---|
| Tables II + III recomputed from Table I (21 values) | worst 0.0012 <!-- claim: marquis_tables_max_rel_dev --> relative | printed digits, truncated: 3 <!-- claim: marquis_above_half_unit_count --> values more than half a unit below the recomputed one, none above |
| in-text C_e on page A3701 | 0.004195 <!-- claim: marquis_C_e_recomputed --> C recomputed | Table III's `4.19e-3 C` confirmed; the sentence's `5.1e-3 C` is wrong |
| SPM voltage error against the DFN, order along the current ladder | 1.00 <!-- claim: spm_order_min --> to 1.04 <!-- claim: spm_order_max --> | O(C_e): order 1 |
| SPMe voltage error, order of the lowest triple (`I = 0.5, 0.25, 0.125`) | 2.08 <!-- claim: spme_order_I05_I0125 --> | O(C_e²): order 2 as C_e → 0 |
| SPMe order, higher triples | 2.18 <!-- claim: spme_order_I1_I025 -->, 2.41 <!-- claim: spme_order_I2_I05 -->, 2.78 <!-- claim: spme_order_I4_I1 --> | pre-asymptotic |
| SPM / SPMe error at I = 1 A/m² | 12.14 <!-- claim: rms_spm_I1_mV --> mV / 0.0451 <!-- claim: rms_spme_I1_mV --> mV (factor 269 <!-- claim: ratio_spm_over_spme_I1 -->) | |
| steady SPMe operator against the closed-form quadratic, 24 + 48 cells | 0.000035 <!-- claim: spme_steady_op_rel_err -->, order 2.00 <!-- claim: spme_steady_op_order_min --> | exact solution |
| SPMe electrolyte profile against the DFN, mid-discharge | 0.057 % <!-- claim: profile_dev_I1_frac --> of c_0 (I = 1), 1.4 % <!-- claim: profile_dev_I5_frac --> (I = 5) | |
| salt conservation (SPMe) / coulomb identity (DFN) | below 1e-11 <!-- claim: salt_conservation_spme --> / below 1e-12 <!-- claim: coulomb_identity_dev --> | identities; their Eq. 32 |
| C_e for this cell at I = 10 | 0.153 <!-- claim: Ce_doyle_at_I10 --> | Doyle's S_s, printed `0.15` in 1993 |

The DFN reference on every rung of the current ladder (4 down to 1/8 A/m²) is
extrapolated in time (nested steps) and in space (12 + 24 to 96 + 192 cells);
what is left of its discretisation error is below 0.01 <!-- claim: dfn_ref_error_over_spme_error_max -->
of the SPMe error it is compared with on every rung, and the same holds for the
SPMe's own space error (below 0.01 <!-- claim: spme_error_space_uncertainty_max -->).
The SPMe is solved exactly in time by the eigendecomposition of its constant
operator. The DFN runs down to the bottom of the ladder with
0 <!-- claim: dfn_failed_solves_sweep --> failed Newton solves.

## What the scaling does *not* establish

The page ships a mutation table for this, because a scaling test is only as
strong as the terms it can see:

- a `10 %` error in the ohmic (δ_c/3) or concentration (factor 2) coefficient
  multiplies the I = 1 error by 9.8 <!-- claim: mut_ohm_error_ratio_I1 --> and
  14.9 <!-- claim: mut_conc_error_ratio_I1 --> and takes the lowest order to
  0.98 <!-- claim: mut_ohm_order_lowest --> and 0.99 <!-- claim: mut_conc_order_lowest -->:
  those two are resolved sharply;
- deleting the electrode-averaged exchange-current correction multiplies the
  error by 1.0001 <!-- claim: mut_exchange_avg_error_ratio_I1 -->: that term is
  *not* tested by anything on this page;
- evaluating the anode overpotential at c_0 multiplies it by
  1.10 <!-- claim: mut_anode_error_ratio_I1 --> (lowest order
  1.90 <!-- claim: mut_anode_order_lowest -->): visible, not sharp;
- raising the reconstructed κ by the gap between its two published routes moves
  the DFN curve by 0.82 <!-- claim: kappa_shift_mV --> mV and multiplies the
  measured error by 0.91 <!-- claim: kappa_plus19pc_error_ratio_I1 -->: shared
  parameters cancel, so **no comparison here tests the parameter set**.

## What the page shows that neither paper does

- All five a-priori applicability conditions (Marquis Table VI) evaluated for
  the 1993 cell: valid at a few A/m², marginal at I = 10 (C_e = 0.15 <!-- claim: Ce_doyle_at_I10 -->,
  κ̂_e = 1.04 <!-- claim: kappa_hat_at_I10 -->), failed at I = 20
  (0.31 <!-- claim: Ce_at_I20 -->, 0.52 <!-- claim: kappa_hat_at_I20 -->), where
  the first-order profile itself predicts zero salt at the collector above
  I = 10.2 <!-- claim: I_depletion_Am2 --> A/m², the transport-limited collapse
  Doyle et al. describe at I = 20 ("about 30 % of the cathode material", their
  text, p. 1529).
- Two of those five conditions are groups Doyle, Fuller & Newman printed
  themselves in 1993: C_e ≡ S_s (Eq. 27) and the solid-diffusion condition
  ≡ S_c (Eq. 26).
- A closed-form (piecewise quadratic) SPMe(S) for the half-cell, and an
  exact-in-time SPMe.

## What is deliberately not attempted

Their Table IV (RMS errors on their graphite/LCO cell) needs OCP and
electrolyte property *functions* that live in their cited codes (DUALFOIL,
fastDFN), not in the printed paper. Those values are carried as context only.

## Files

- `index.py` / `index.ipynb`: the page (paired with jupytext; edit `index.py`)
- `model.py`: J3.4's DFN (copied unchanged) and the reductions (`v_spm`, `SPMe`,
  `bv_invert`); `test_model.py` its fast checks
- `data/marquis-2019-table1-parameters.csv`: their Table I (+ sidecar)
- `data/marquis-2019-stated-results.csv`: their Tables II to IV and one in-text
  value (+ sidecar)
- physical inputs load from `J3.4`'s parameter and stated-results tables

## Corrections

- 2026-11-08: said "The bottom of the sweep (I = 0.25 A/m²) is set by the DFN
  integrator, which does not run at I = 0.125 at any time step tried". The Newton
  solves at `I = 0.125` converged; the page's own acceptance test rejected them, an
  absolute threshold on a residual scaled by factors that grow as
  1/I. J3.4's DFN, which checks Newton's own convergence, runs the whole ladder
  down to `I = 0.125` with 0 <!-- claim: dfn_failed_solves_sweep --> failed solves,
  and the bottom triple gives the SPMe order 2.08 <!-- claim: spme_order_I05_I0125 -->.
- 2026-11-08: said the SPM's slope 1 is structural because "with constant kappa
  and D every term the SPM omits is exactly proportional to I". Only the
  ohmic and the quasi-steady concentration corrections are; the SPM also omits
  terms that grow faster, the SPMe's own error among them, and the SPM order rises to 1.04 <!-- claim: spm_order_max --> at the
  top of the ladder.
