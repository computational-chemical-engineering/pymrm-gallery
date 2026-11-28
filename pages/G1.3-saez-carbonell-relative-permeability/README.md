# G1.3: Sáez and Carbonell (1985), trickle-flow holdup and pressure drop from relative permeabilities

**Source.** Sáez, A. E. and Carbonell, R. G. (1985), *Hydrodynamic parameters for
gas-liquid cocurrent flow in packed beds*, AIChE Journal **31**(1) 52-62,
[doi:10.1002/aic.690310105](https://doi.org/10.1002/aic.690310105). Read on the
native 300 ppi page images; journal page = PDF page + 51 (PDF page 1 opens on the end
of the previous article).

- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (groups, the stagnant-gas closed form, the coupled pair with its analytic
  Jacobian solved by `pymrm.newton`, the scalar elimination solved by `brentq`, the
  Table 5 ratios and the printed-arithmetic helpers), tested by `test_model.py`.
- `data/`: the printed numbers (strings), Figure 8 digitised by `data/digitise.py`
  (it takes the PDF path; no image is stored), and independent check values, each
  with a provenance sidecar.

## The question

Do the paper's printed results follow from its printed relations? The model is the
pair of two-fluid momentum balances with saturation-only relative permeabilities,
$k_L = \delta^{a}$ and $k_G = S_G^{c}$ with the two printed exponents, which gives holdup (Eq. 51) and pressure
drop (Eq. 50) together.

## What the page establishes

| | |
|---|---|
| derivation chain | re-derived; Eq. 51 is Eq. 21 without the gas weight, which moves the holdup by -2.0155e-4 <!-- claim: T_eps_L_shift_eq51_vs_eq21 --> relative at the test condition; Eq. 34's rounded exponent sits 0.305 <!-- claim: exp34_halfunits_from_0p41 --> half-units of its last digit from 1/a and moves the stagnant-gas δ by a relative 0.0036895 <!-- claim: T_delta_shift_041_vs_inv_a --> |
| coupled model at the test condition of the notebook (air-water, spheres) | holdup 0.15705 <!-- claim: T_eps_L -->, ψ_G 292.88 <!-- claim: T_psi_G -->, -dP/dz 3436.1 <!-- claim: T_minus_dPdz_Pa_per_m --> Pa/m, coupling term 0.35217 <!-- claim: T_coupling --> |
| sensitivity of the pressure drop to the holdup | d ln ψ / d ln ε_L = 3.3813 <!-- claim: T_dlnpsi_dlneps -->; a holdup high or low by the paper's mean holdup deviation moves ψ_G by a relative 0.5681 <!-- claim: T_psi_factor_plus_12p7 --> or -0.33721 <!-- claim: T_psi_factor_minus_12p7 --> |
| Table 5, largest ratio the corrected model can give (densest packing) | 1.0361 <!-- claim: t5_corrected_max_sol1_dense --> and 1.0510 <!-- claim: t5_corrected_max_sol2_dense -->, against the printed Calc. values |
| measured ratios above that bound | 8 <!-- claim: t5_n_exp_above_bound --> of 8; smallest margin 1.0328 <!-- claim: t5_exp_min_over_bound_sol1 --> times the bound |
| where the printed Eq. 35 returns both Calc. values | one root, porosity 0.70506 <!-- claim: t5_printed_eps -->, Eo*_o 19.006 <!-- claim: t5_printed_Eo_o -->; correct ratios there 1.0121 <!-- claim: t5_corrected_at_printed_sol1 --> and 1.0168 <!-- claim: t5_corrected_at_printed_sol2 --> |
| printed statistics recombined | 4 <!-- claim: stats_n_within_allowance --> of 4 inside their rounding allowance; Table 7 holdup at 0.763 <!-- claim: t7_eps_gap_over_allowance --> of it |
| Figure 8, readable markers | 13 <!-- claim: fig8_n_readable --> of 59; 7 <!-- claim: fig8_n_above_plus --> above the upper dashed line, 2 <!-- claim: fig8_n_below_minus --> below the lower one; every blob of touching markers inside the band |
| newton against brentq over the sweep | below 1e-10 <!-- claim: sweep_newton_vs_brentq --> |

## The finding

The paper writes "From these results we can see that the ratio of dynamic holdups
tends to be overestimated by our correlation" (p. 58). The printed Eq. 35 drops a
factor Eo*_o; with it restored, the model's ratio is bounded for every packed bed and
lies below all eight measured ratios, so the corrected model under-predicts them. The
paper writes that "the liquid holdup data lie between the −20% and the +20% lines
(Figure 8)" (p. 59); of the readable markers, 9 <!-- claim: fig8_n_outside --> lie
outside the drawn band, each by more than its reading error. Fitted constants are
not refitted, and nothing on the page says whether saturation-only relative
permeability is right for trickle beds.

## Not done

Figure 7 (pressure drop, 49 points) and Figure 5 (the refit of the gas exponent) are not
digitised on this page.
