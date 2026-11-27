# J3.6: Springer, Zawodzinski and Gottesfeld (1991), water balance and resistance of a PEM fuel cell membrane

T. E. Springer, T. A. Zawodzinski and S. Gottesfeld, "Polymer electrolyte fuel cell model",
*J. Electrochem. Soc.* **138**, 2334-2342 (1991), doi:10.1149/1.2085971. The page asks whether
the membrane results of Table I, `alpha = 0.4` (`beta = 0.2`) and `R_mem = 0.285` Ohm cm2,
follow from the paper's printed membrane model (Eqs. 16 to 18, 22 to 26) with the printed
interface water fractions `x_w2 = 0.1013` and `x_w3 = 0.2327` as boundary data. The electrodes
and the polarisation curves cannot be rebuilt from the paper (porosity, critical constants and
the `lambda <= 4` diffusion table are not printed) and are not judged.

## Findings

- Read literally, Eqs. 23 and 24 count the factor `lambda/22` twice (Eq. 18 already holds it;
  p. 2340 gives `n_drag = 1.59` at `lambda = 14`). Then a rising profile needs alpha below the
  drag at the anode face, 0.2262 <!-- claim: r5_literal_bound -->, for any `rho_dry/M_m`: a
  scan of `c_f` from `1e-4` to `1e-2` mol/cm3 finds 0 <!-- claim: r5_literal_n_roots --> roots of `alpha = 0.4`.
- With the drag counted once, `alpha = 0.4` needs `rho_dry/M_m` =
  0.0016942 <!-- claim: r4_cf_star --> mol/cm3 (not printed; a dry density of
  1.8636 <!-- claim: r4_rhodry_star --> g/cm3 at an assumed equivalent weight of `1100`).
  There `R_mem` = 0.27875 <!-- claim: r4_R_star --> Ohm cm2 from Eq. 26 on the dry coordinate,
  -2.19 <!-- claim: r4_R_gap_pct --> % against the printed `0.285`, which is
  -12.5 <!-- claim: r4_R_hu --> half-units of its last digit. Over the printed half-unit band
  of alpha (`0.35` to `0.45`) it stays between 0.27718 <!-- claim: r4_R_band_lo --> and
  0.28048 <!-- claim: r4_R_band_hi -->; rounding the interface fractions and using the
  printed inlet fraction `0.1558` keep it below `0.285` in 7 <!-- claim: r4_sens_n_below -->
  cases out of `7`.
- Read with the swollen thickness, `R_mem` = 0.30470 <!-- claim: r6_R_swollen_star -->, so
  the printed value lies between the two readings, at 0.241 <!-- claim: r6_bracket_pos --> of
  the way. The paper does not say which it used; the page names no cause for the gap.
- The gas side of Table I closes on alpha: Eqs. 7 and 8 back-solved from `x_w1`, `x_w4` and
  `x_O4` give 0.3993 <!-- claim: r2_alpha_xw1 -->, 0.3980 <!-- claim: r2_alpha_xw4 --> and
  0.3968 <!-- claim: r2_alpha_xO4 -->. Eq. 27 gives `V_cell` = 0.66844 <!-- claim: r7_vcell_xliq077 -->
  V with the printed `x_liq,4 = 0.077`, 0.88 <!-- claim: r7_vcell_hu --> half-units from `0.668`.
- Table I's inlet fraction `0.1558` is not Eq. 15's `P_sat/P` = 0.155642 <!-- claim: r1_xsat_80C -->,
  and Fig. 3 confirms Eq. 21, not Eq. 20 (Eq. 21's 2.4252 <!-- claim: r1_Dratio3 --> at
  `lambda = 3` is 0.10 <!-- claim: r1_fig3_dev --> reading errors from the read peak).

## Things to know before reusing this

- `rho_dry/M_m` is not printed: the list of symbols (p. 2341) gives `M_m` and `rho` without values, Table I and Eqs. 19 to 26 give none, and the paper's one equivalent weight, `1155` (p. 2337), is that of Yeo and Eisenberg's Nafion. Any value used with this model is an assumption or a fit.
- Eq. 22 holds only for `lambda > 4`; the paper's table below that is not printed, so
  high-current cases with a drying anode cannot be run from the paper.
- The Fig. 4 comparison is illustrative only: its caption holds gas flows constant at
  `nu_H = 1` at `1` A/cm2, so its `0.5` A/cm2 curve is at `nu_H = 2` by the flow rule of p. 2338, against Table I's `4`.
- The gas constant of Eq. 27 is not printed; the page uses `8.314` J/(mol K).
