# J2.7: Coulaloglou and Tavlarides (1977), drop breakage and coalescence in a stirred vessel

C. A. Coulaloglou and L. L. Tavlarides, "Description of interaction processes in agitated
liquid-liquid dispersions", *Chem. Eng. Sci.* **32**, 1289-1297 (1977),
doi:10.1016/0009-2509(77)85023-9. The page rebuilds the paper's breakage, daughter,
collision and film-drainage kernels in a steady population balance for a continuous-flow
stirred vessel (fixed pivots, `pymrm.newton`), fixes what the paper does not print from its
own printed model numbers, and asks whether the scaling it derives for d_max from breakage
alone, N^(-6/5) where the exponential of eq. 14 dominates (p. 1291), holds for d32 in its
model, curves and data over 190 to 310 rpm; the paper does not claim it for d32.

## Findings

- The paper prints no physical properties. With S = sigma/rho_d calibrated to the printed
  calculated d32 of Fig. 2 (S = 48.6847 <!-- claim: r3_S --> cm3/s2) and X = mu_c rho_c/sigma^2 from the printed
  efficiency "6 %", the model's omega_v is 1.0228 <!-- claim: r3_omega_v_ratio --> times the printed calculated
  value in 1/s, inside the band of the rounding of "6 %" (0.9914 <!-- claim: r3_omega_v_band_lo --> to
  1.0528 <!-- claim: r3_omega_v_band_hi -->); read in 1/min it is 61.4 <!-- claim: r3_omega_v_ratio_if_per_min --> times off,
  and counted per cm3 0.0048 <!-- claim: r3_pv_ratio_if_per_min --> times even in 1/min. Counts per vessel and
  frequencies in 1/s are the reading that reproduces the paper.
- The exponent of d32 in N* over the paper's range is not -6/5. The full model gives
  -1.075 <!-- claim: r5_m_full_0.05 -->, -1.042 <!-- claim: r5_m_full_0.10 --> and -1.020 <!-- claim: r5_m_full_0.15 --> (phi `0.05`, `0.10`, `0.15`);
  breakage alone with the residence time -1.309 <!-- claim: r5_m_brk_0.05 --> to -1.311 <!-- claim: r5_m_brk_0.15 -->; the
  closed form at d_max -1.301 <!-- claim: r2_slope_0.05 --> to -1.303 <!-- claim: r2_slope_0.15 -->, since y at d_max is
  10.46 <!-- claim: r2_y_min --> to 11.55 <!-- claim: r2_y_max --> and -6/5 is the limit y -> infinity. Even breakage alone
  misses -6/5 at the paper's y, and d32 is flatter still. At the finite y that g theta = 1 sets,
  the prefactor of eq. 14 steepens the exponent; coalescence flattens it by
  0.234 <!-- claim: r5_shift_coal_0.05 --> to 0.292 <!-- claim: r5_shift_coal_0.15 -->.
- The paper's own data lie between: flatter than -6/5, steeper than the reconstruction (phi `0.10`:
  -1.092 <!-- claim: r5_data_m_0.10 --> against -1.042 <!-- claim: r5_m_full_0.10 -->). The paper's model lines give
  -1.114 <!-- claim: r5_line_m_0.05 --> to -1.074 <!-- claim: r5_line_m_0.10 -->. Pooled over the three holdups the
  slope is -1.0951 <!-- claim: r5_pooled_m --> ± 0.0106 <!-- claim: r5_pooled_se -->, t = 9.9 <!-- claim: r5_pooled_t --> against -1.2, and
  3 <!-- claim: r5_n_data_excl --> of the 3 holdups exclude -6/5 on their own at the two-sided `95 %` point of Student's t. Each
  data slope carries the larger of its regression error and the error that the figure reading
  passes on.
- The reconstruction lies within the reading error of the printed model lines of Fig. 3 at
  1 <!-- claim: r4_fig3_n_inside --> of 15 <!-- claim: r4_fig3_n_points --> conditions (largest deviation 0.0241 <!-- claim: r4_fig3_max_dev_mm --> mm) and of
  Fig. 4 at 1 <!-- claim: r4_fig4_n_inside --> of 9 <!-- claim: r4_fig4_n_points --> (mean ratio 1.025 <!-- claim: r4_fig4_mean_ratio -->). Its Fig. 2 density
  is narrower than the printed one: 0.00146 <!-- claim: r4_fig2_frac_above_04_model --> of the volume above `0.4 mm`
  against 0.054 <!-- claim: r4_fig2_frac_04_0495_read --> between `0.4` and `0.495 mm` alone in the printed curve, a lower
  bound for its tail. The page names no cause.
- The paper's marching iteration with coalescence does not converge at the Fig. 2 point within
  its sweep limit. With the paper's stop rule (0.001 <!-- claim: v4_march_inc --> per increment, p. 1293)
  its marching iteration stops after 357 <!-- claim: v4_march_stop_sweeps --> sweeps with d32 lower than the steady state by a
  relative 0.0735 <!-- claim: v4_march_stop_gap -->. Newton from the breakage-only state converges at every condition of the page.

## Things to know before reusing this

- C1 to C4 were fitted to the drop size distributions and coalescence frequencies at phi `0.05`
  (p. 1294), and S here is calibrated
  to one printed model result: agreement with those figures is not a validation.
- C3 carries the 12 l vessel volume (counts per vessel); scale-up to another vessel is not tested.
- Eq. 25 is used with the exponent 4 (eqs. 23, 24, 27), eq. 32 with its sign and variance
  corrected, and the daughter distribution normalised on [0, v'].
- The fixed pivot converges at second order (V3); the 15 conditions use q = 8 and 16 with one
  Richardson step.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | the kernels (eqs. 36, 37), d_max, the paper's estimates, `SteadyPBE` (fixed-pivot steady balance with analytic Jacobian, solved by `pymrm.newton`) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/coulaloglou-1977-*.csv` | printed numbers, Figs. 3 and 4 symbols and model lines, the Fig. 2 theory curve, each with its provenance sidecar |
| `data/check_values.csv` | values of the page's functions at one test condition from a second route, asserted before any claim |
