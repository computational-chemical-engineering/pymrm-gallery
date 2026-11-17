# F3.5: CO2 absorption into amine-promoted carbonate solutions

Why adding a little amine to hot potash multiplies the CO2 absorption flux,
answered by solving the film with every reversible reaction in it: the model
that ended the shuttle-versus-homogeneous-catalysis argument.

- **Structure:** `S3` (1D steady BVP) + `S10` (instantaneous-equilibrium
  constraints)
- **Reference:** Bosch, Versteeg & van Swaaij (1989), Chem. Eng. Sci. 44(11)
  2735-2743, doi:10.1016/0009-2509(89)85216-9. (The catalogue's original
  citation, "Versteeg & van Swaaij (1988)", is the group's kinetics paper,
  the source of this paper's correlation (11), not the model paper.)

## Results

With Table 1 exactly as printed, every one of the twelve enhancement factors
comes out high, by up to 21.4% <!-- claim: E_maxdev_printed -->, systematically.
The paper's own identity J = E kL (ci − cb) implies a bulk dissolved CO2 above
what the printed constants give, and the model's structure localises the
discrepancy to the single constant K3 = Kw/Kc2 (whose rescaling fixes the bulk
CO2 and, through the bulk OH−, the Hatta number at once). One scalar, the root
of the mean log residual of that identity, s = 1.4585 <!-- claim: k3_scale_fitted -->,
gives:

| comparison | largest deviation (model − paper)/paper |
|---|---|
| **unpromoted flux J_u (4 values), headline** | **3.8% <!-- claim: Ju_maxdev_reconstructed -->** |
| unpromoted E (4 values) | 0.6% <!-- claim: Eu_maxdev_reconstructed -->, see below |
| promoted E, DEA + HDA (8 values) | 7.6% <!-- claim: E_promoted_maxdev_reconstructed --> |
| promotion factors F (8 values) | 5.6% <!-- claim: F_maxdev_reconstructed --> |

**Read the E_u number with care.** The scalar is fitted to
cb = ci − J/(E kL), so it consumes the printed fluxes *and* the printed
enhancement factors. The identity does not admit one scalar exactly: the four
per-condition scalars spread by 5.0% <!-- claim: k3_scale_spread_per_condition -->
(n = 4 conditions), and that residual lands in the flux, because the fitted cb
sits in E's denominator and makes E less sensitive to s than J is. The four
E_u deviations are one offset repeated, not four independent agreements;
refitting to E_u alone gives 1.4735 <!-- claim: k3_scale_refitted_to_Eu -->. The
flux is the metric to read.

**Where the scalar sits.** On K_c1 instead it fixes the bulk CO2 but not the
OH−, and E_u stays 19.8% <!-- claim: Eu_max_branch_Kc1 --> high: rejected. On
K_w instead of K_c2 every unpromoted quantity is identical, and the promotion
factors miss by 8.6% <!-- claim: Fmax_branch_Kw --> against
5.6% <!-- claim: Fmax_branch_Kc2 --> on K_c2: the split is not resolved.

**Evidence for, using no enhancement factor at all:** the unpromoted bulk CO2
read off Figures 2, 3 and 4 implies a mean scalar of
1.445 <!-- claim: s_from_figures_mean -->. The figure readings are a
cross-check *pending maintainer visual review*; no dataset is shipped. They
pass two charge balances, each to at most 5e-3 <!-- claim: fig_charge_dev_unpromoted_max -->
(unpromoted) and at most 5e-3 <!-- claim: fig_charge_dev_promoted_max -->
(promoted).

**Evidence against, stated plainly on the page:** the paper's only two
out-of-sample numbers both prefer the constants exactly as printed. The
desorption promotion factor misses by -2.3% <!-- claim: desorption_F_dev_printed_constants -->
with the printed constants against -12.0% <!-- claim: desorption_F_dev -->
reconstructed, and the amine cuts the driving force by
34.7% <!-- claim: desorption_drive_cut_printed_constants --> against
39.1% <!-- claim: desorption_drive_cut --> for the paper's "about one third".
The paper's own record is internally inconsistent here, since its Figure 4
supports the reconstruction at the same condition.

## Validation

Ranked by what each check can catch. Every one has a row in the notebook's
break table that fails it.

- **Pseudo-first-order film** against its closed form over a Hatta sweep
  from a reaction layer ten times the film to a small fraction of it, at most 1e-4 <!-- claim: pfo_max_rel_err --> relative; the check that
  resolves the reaction layer. Grid doubling moves E of the stiffest case by at
  most 1e-4 <!-- claim: grid_delta -->, a further decade of k_inst by at most
  1e-4 <!-- claim: kinst_delta -->.
- **Van Krevelen-Hoftijzer**, an external approximation, within at most
  1e-3 <!-- claim: vkh_dev --> at the DEA Hatta number, nothing tuned to it.
- **Penetration instead of film.** The same network as a Higbie penetration
  model (backward Euler on a geometric time grid, `newton` per step) gives E
  larger by up to 11.3% <!-- claim: pen_vs_film_max --> at n_t = 640; the time
  orders lie between 0.9 and 1.1 <!-- claim: pen_time_order --> on the ladder
  160 to 1280 (finest successive difference at least
  1e3 <!-- claim: pen_ladder_over_newton_tol --> times the Newton tolerance), and
  the Richardson limit from the finest triple is
  11.57% <!-- claim: pen_vs_film_max_richardson -->, and the finest level, n_t = 1280,
  falls short of it by 1.3e-3 <!-- claim: pen_richardson_error_estimate -->.
- **Consistency checks, true by construction:** physical absorption E = 1 to at
  most 1e-9 <!-- claim: physical_E1 -->; a dimensionless two-species assembly
  sharing grid and operators, at most 1e-6 <!-- claim: second_order_cross -->;
  the reversible-to-irreversible collapse, at most
  1e-4 <!-- claim: collapse_residual -->. Carbon-flux uniformity (at most
  1e-8 <!-- claim: carbon_flux_closure -->), the amine flux (at most
  1e-8 <!-- claim: amine_flux_closure -->) and the charge flux (at most
  1e-8 <!-- claim: charge_flux -->) test the
  stoichiometry matrix, not the diffusivities; electroneutrality, at most
  1e-9 <!-- claim: electroneutrality -->, tests the equal ion diffusivities.
- **Physical bounds:** all 12 <!-- claim: cases_below_ceiling --> cases lie below
  their instantaneous-reaction ceiling; the gas-side share of the resistance is
  at most 1e-4 <!-- claim: gas_side_fraction -->, which licenses the Dirichlet
  interface. The printed F agree with the ratios of printed fluxes to at most
  5e-3 <!-- claim: table2_F_vs_Jratio_maxdev -->, print rounding.

## What pymrm adds

The 1989 calculation was a dedicated research code; here it is `model.py` on
pymrm operators: `construct_grad` with per-species boundary conditions on a
graded grid, `construct_div`, `NumJac` on the pointwise reaction network and
`newton` with bounds, continued in k_inst, with the flux read from the
face-flux operator. The same classes run absorption, desorption and the
penetration variant (`newton` per backward-Euler step), and reduce to the
classical limits on demand; `pymrm.checks` supplies the roots of both fits and
the time-refinement orders. The model itself is the paper's: no improvement
is claimed. What the reimplementation adds is the forensic reading of the
parameter set above, with the evidence against it reported alongside.

## Data

**Provenance tier 6: reproduction, not validation.** The paper tabulates no
measurements; its Table 2 is the authors' own numerics and its only
experimental numbers print as "±4" and "±6" (almost certainly "≈4" and "≈6").
Both tables are transcribed from 600 dpi page renders (the scan's OCR drops
exponents) with the transcription checks described in the sidecars.

## Corrections

- 2026-08-05: the film-vs-penetration difference was published as "~10.5%" from a single run at n_t = 160. Backward Euler is first order in n_t; the difference is 11.3% <!-- claim: pen_vs_film_max --> at n_t = 640 and extrapolates to 11.57% <!-- claim: pen_vs_film_max_richardson -->.
- 2026-08-05: the carbon- and amine-flux closures were described as "genuine" because "they mix the three different diffusivities". Every reaction conserves carbon and the amine count, so the closures hold for any set of diffusivities; they test the stoichiometry matrix.
- 2026-11-17: stated that physical_E1, second_order_cross, collapse_residual, charge_flux and electroneutrality "have no break row and none is constructible". Each now has one: the film made thicker than D_CO2/k_L, one amine per CO2 in reaction (2), the K_eq scan stopped early, the stoichiometry rows and a halved D(HCO3−).
- 2026-11-17: stated "Electroneutrality (2.8e-11) and charge flux (1.1e-9) are algebraically guaranteed by the equal ion diffusivities". The charge flux vanishes for any diffusivities, because every reaction conserves charge and only neutral CO2 crosses the interface; only electroneutrality needs the equal ion diffusivities.

## Files

| file | contents |
|---|---|
| `index.py`, `index.ipynb` | the page (paired by jupytext) |
| `model.py` | bulk speciation, the film and penetration models, closures and limits |
| `test_model.py` | fast offline checks of `model.py` |
| `data/bosch1989-table1-parameters.csv` | Table 1, every constant |
| `data/bosch1989-table2-results.csv` | Table 2, the paper's numerical results |
| `agreement.json` | every claim, written by the notebook |

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```
