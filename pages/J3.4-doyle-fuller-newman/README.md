# J3.4: Doyle–Fuller–Newman lithium/polymer/insertion cell (P2D)

The 1993 pseudo-two-dimensional cell model that every physics-based lithium-ion
battery model descends from, rebuilt in pymrm from the paper, including the one
transport property the paper never prints.

- **Structures:** `S8` (coupled multi-phase transport with electroneutrality),
  `S10` (porous electrode with a distributed interfacial reaction)
- **Reference:** Doyle, Fuller & Newman (1993), *J. Electrochem. Soc.* **140**(6)
  1526–1533, doi:10.1149/1.2221597
- **Runtime:** about 4 minutes

## Provenance: tier 6, not experimental

**Figure 2 contains no measured points.** Every curve on it (charge,
open-circuit and four discharges) is the authors' own simulation; the text cites
measured discharge curves by reference only and plots none. This was confirmed
by the gallery maintainer on review of the original figure. The page is a
reference-solution reproduction in the sense of `D2.2`. Do **not** describe it as
validated against experiment.

## Agreement

| what | pymrm | paper |
|---|---|---|
| Eq. 16 against the plotted open-circuit curve | 3.25 <!-- claim: ocp_eq16_mad_mV --> mV mean abs. dev. | digitisation good to about `3 mV` |
| four discharge curves, every digitised point in range | 25.5 <!-- claim: discharge_pooled_mad_mV --> mV pooled, bias -23.8 <!-- claim: discharge_pooled_bias_mV --> mV | |
| *u* at the `1.9 V` crossing, all four currents | short by 0.023 <!-- claim: u19_shortfall_mean --> on average, 0.014 <!-- claim: u19_shortfall_min --> to 0.031 <!-- claim: u19_shortfall_max --> | |
| charge curve | 65.1 <!-- claim: mad_charge_mV --> mV mean abs. dev., bias 64.9 <!-- claim: charge_bias_mV --> mV | |
| *S*_c (Eq. 26) at *I* = 10 | 1.02e-4 <!-- claim: S_c_at_I10 --> | `1.0e-4` |
| *S*_s (Eq. 27) at *I* = 10 | 0.153 <!-- claim: S_s_at_I10 --> | `0.15` |
| σ/κ | 2.1e5 <!-- claim: sigma_over_kappa --> | of order `1e5` |
| *u* at the `1.7 V` cutoff, *I* = 10 | 0.831 <!-- claim: u_at_cutoff_I10 --> | `0.84` stated; 0.821 <!-- claim: u17_I10_digitised --> digitised |
| optimum cathode porosity | 0.591 <!-- claim: porosity_optimum --> | `0.60` |
| utilisation at that optimum | 0.928 <!-- claim: u_at_optimum --> | `0.97` |
| short-time salt at the lithium face against the semi-infinite solution | 0.83 <!-- claim: short_time_face_dev_mean --> mol/m³ mean | |

Numerically: the spatial error of the production grid in the cell potential at
`u = 0.5` is 0.028 <!-- claim: grid_error_mV --> mV (second order), the step
error 0.27 <!-- claim: dt_spread_mV --> mV (first order), each from its own nested
ladder with the other axis fixed. On the headline itself the step dominates: the
production steps leave the pooled deviation 0.8 <!-- claim: headline_dt_error_pooled_mV -->
mV below its zero-step limit, and the `I = 20` deviation 1.9 <!-- claim: headline_dt_error_I20_mV --> mV. Both are small against the model-figure
deviation.

## What the identities test

Three numbers on the page are identities of the scheme, each exact for the
interior physics and each testing one thing:

- **Total salt** is constant to below 1e-12 <!-- claim: salt_conservation_drift -->.
  With the salt flux zero at both ends, `sum(dx * div N)` telescopes to zero for
  any interior flux, so the migration term, the transference number and the
  reaction coupling cannot move it (structural, measured under all three). What
  it tests is the two salt-flux boundary conditions.
- **Faraday's law**: the utilisation computed from the state equals the charge
  passed, to below 1e-12 <!-- claim: faraday_gap -->. It is the only check that
  sees the pore-wall current entering the solid balance differently from the
  charge balance.
- **The mean pore-wall flux** equals `-I/(F a δ_c)` at any time. It sees the
  current in the charge balance differing from the one reported.

The check that sees the electrolyte physics is the short-time comparison at the
lithium face: no conductivity in it and nothing fitted. A control run with the
transference number held constant, where the semi-infinite solution is exact,
agrees to below 0.2 <!-- claim: short_time_control_dev_max_10s --> mol/m³ up to
`10 s`; what remains in the page's own run is the concentration dependence of
`t0+`.

## The gap in the paper, and how it is closed

**Appendix A never gives the conductivity.** It says the conductivity of
PEO-LiCF₃SO₃ "was fit to a third-order polynomial" and cites its source; the
coefficients appear nowhere in the article and there is no figure of κ against
*c*. The salt diffusivity, the solubility limit and the transference-number fit
are all printed. The conductivity is not.

It is recovered by inverting the paper's own Eqs. 28 and 29, which contain
(1/κ + 1/σ) and whose values the paper prints. The two routes give
κ_eff = 7.97e-3 <!-- claim: kappa_eff_from_eq28 --> and 6.72e-3 <!-- claim: kappa_eff_from_eq29 --> S/m, 18.7 <!-- claim: kappa_route_spread_pct -->
% apart, and both put σ/κ at the paper's order of magnitude. What cannot be
recovered is the *concentration dependence*, and the opposite biases on charge and
discharge are where that most likely shows.

## Two ambiguities in the printed model

1. **The Bruggeman exponent in the cathode salt flux.** Table I writes the flux
   as `∇·(εD∇c)` while the text defines `D_eff = Dε^0.5`. Read together they give
   `Dε^1.5`, which exhausts the electrolyte at `u` = 0.193 <!-- claim: brug_eps15_u_salt_out_I10 --> at `I = 10`, contradicting the
   paper's own Figure 3, and puts the `1.9 V` crossing at 0.38 <!-- claim: brug_eps15_ratio_min --> to 0.62 <!-- claim: brug_eps15_ratio_max -->
   of the published utilisation (mean |Δu| 0.294 <!-- claim: brug_mean_du_eps15 -->
   against 0.023 <!-- claim: u19_shortfall_mean -->). Dropping the correction
   entirely overshoots at `I = 10`, 13 and 20 (mean |Δu| 0.086 <!-- claim: brug_mean_du_eps0 -->). The page uses `Dε^0.5`.
2. **A missing factor ν = 2 in Eq. 5.** As printed, the coefficient of ∇ln c is
   half what a 1:1 salt requires: at zero current the equation must reproduce the
   emf of a lithium concentration cell, which is (2RT/F)∫(1−t⁰₊)dln c. The page
   uses the corrected form. **Figure 2 does not settle this one**: the printed
   version gives a smaller mean deviation (18.2 <!-- claim: nu1_mad_mV --> against
   25.5 <!-- claim: discharge_pooled_mad_mV --> mV) but puts every
   transport-limited crossing late, by 0.026 <!-- claim: nu1_u19_bias --> on
   average, where the corrected form puts them early. The two bracket the
   published curves at every current.

## What the page adds

The porosity optimisation the paper reports in one sentence, swept over the
window it discusses and refined by a bounded maximisation: the maximum lies at a
porosity of 0.591 <!-- claim: porosity_optimum -->, just below the paper's, and
it is flat, with the loss above the optimum the larger at equal distance,
0.0216 <!-- claim: loss_above_015 --> against 0.0163 <!-- claim: loss_below_015 --> at `0.15`. No new physics: the model is theirs, unchanged.

## Corrections

- 2026-11-04: said of the short-time check "Mean 7.7 mol/m3 over five times, no
  conductivity in it and nothing fitted". The model value was read at the centre
  of the first cell, which sits 5.8 <!-- claim: short_time_first_cell_offset_mean -->
  mol/m³ below the lithium face at the page's grid, and the steps were the
  production ones. Read at the face under the boundary condition, with the steps
  refined 16-fold, the mean deviation is 0.83 <!-- claim: short_time_face_dev_mean --> mol/m³, nearly all of it the
  concentration dependence of `t0+`.
- 2026-11-04: listed the paper's side of the short-time comparison as "Fig. 4 to ±5",
  a hand read's precision, as if model and figure agreed to it. The model's face
  value differs from that read by up to 15.7 <!-- claim: fig4_dev_max --> mol/m³,
  beyond the reading at 4 <!-- claim: fig4_beyond_reading_count --> of the five times.
- 2026-11-04: said "The porosity optimum is found at 0.60, exactly where the paper
  puts it". That was the largest node of a scan in steps of `0.05`; the maximum,
  refined by a bounded maximisation, is at 0.591 <!-- claim: porosity_optimum -->.
- 2026-11-04: said "the penalty for being below the optimum is far steeper than for
  being above it". At equal distance from the optimum the loss is larger above it,
  0.0216 <!-- claim: loss_above_015 --> against 0.0163 <!-- claim: loss_below_015 --> at `0.15`.
- 2026-11-04: said "every transport-limited crossing arrives 2–3 percentage points
  of utilisation early, at every current". The shortfalls in utilisation are 0.014 <!-- claim: u19_shortfall_min --> to 0.031 <!-- claim: u19_shortfall_max -->.
- 2026-11-04: said "the 1.9 V crossing comes out at a third to a half of the
  published utilisation at every current" for `Dε^1.5`. It is 0.38 <!-- claim: brug_eps15_ratio_min --> to 0.62 <!-- claim: brug_eps15_ratio_max -->
  of it.
- 2026-11-04: said "the charge curve is over-polarised by 64 mV in the opposite
  direction". With the face values and the march ending at the exhausted-solid
  state, the mean absolute deviation is 65.1 <!-- claim: mad_charge_mV --> mV and
  the bias 64.9 <!-- claim: charge_bias_mV --> mV.
- 2026-11-04: called total salt conservation "the cleanest check on the page". It is an
  identity of the scheme for any interior flux and tests only the two salt-flux
  boundary conditions.
