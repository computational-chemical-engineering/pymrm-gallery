# F2.6 Van 't Riet (1979): the sulfite method as a stirred-tank gas-liquid reactor

van 't Riet, K. (1979), "Review of measuring methods and results in nonviscous gas-liquid mass
transfer in stirred vessels", *Ind. Eng. Chem. Process Des. Dev.* 18(3), 357-364,
doi:10.1021/i260071a001.

The review judges the sulfite method "even more restricted" than the other ways of measuring
k_L A, because the apparent k_L A must be converted with a rate constant and an order that are
hard to know. This page builds the reactor the method uses (not in the review): a stagnant film
coupled to a well-mixed tank, with an m-th order reaction of the dissolved oxygen and k_L A from
the review's eq. 9. One group decides everything, Al = k_L^2/(D k_L A), from printed numbers plus
an assumed diffusivity. The page asks, inside the review's own window (P/V and v_s of eq. 9 and
Figure 4), whether there is a band of rate constants where the method reads eq. 9's k_L A within
its printed accuracy, "accurate for 20-40%", without any conversion, how wide it is, and whether the printed 25-fold shift of k_R fits:

- at the test condition (centre of eq. 9's P/V range, Figure 4's Reith sulphite line,
  `k_L = 4e-4 m/s`) Al = 333.43 <!-- claim: Al_T -->; the band width W = k_hi/k_lo is
  54.548 <!-- claim: m1_tol20_W_T --> (first order) and 20.429 <!-- claim: m2_tol20_W_T -->
  (second order) at `20 %`, 292.18 <!-- claim: m1_tol40_W_T --> and
  175.26 <!-- claim: m2_tol40_W_T --> at `40 %`: the 25-fold shift fits in
  3 <!-- claim: fits_at_T_count --> of the four cases;
- over the box, Al runs from 116.85 <!-- claim: box_kL4e-4_Al_min --> to
  1521 <!-- claim: box_kL4e-4_Al_max -->, W at `20 %` from 20.684 <!-- claim: box_m1_tol20_kL4e-4_W_min -->
  to 240 <!-- claim: box_m1_tol20_kL4e-4_W_max --> (first order) and from
  9.7118 <!-- claim: box_m2_tol20_kL4e-4_W_min --> to 74.822 <!-- claim: box_m2_tol20_kL4e-4_W_max -->
  (second order); the shift fits on the whole box in 2 <!-- claim: fits_whole_box_count --> of
  the four cases (both at `40 %`), and at `20 %` on 0.97722 <!-- claim: box_m1_tol20_kL4e-4_frac_W_ge_25 -->
  and 0.49087 <!-- claim: box_m2_tol20_kL4e-4_frac_W_ge_25 --> of it (log-area shares);
- the review's "This implies that the reaction largely takes place in the boundary layer"
  depends on the order: where the bulk is nil to `1 %`, the film
  carries 0.031989 <!-- claim: box_m1_nil01_film_share_at_Al_max --> to
  0.29101 <!-- claim: box_m1_nil01_film_share_at_Al_min --> of the reaction for the first order
  (the implication fails) and 0.62338 <!-- claim: box_m2_nil01_film_share_at_Al_max --> to
  0.92285 <!-- claim: box_m2_nil01_film_share_at_Al_min --> for the second (it holds); where
  the film carries half, the method reads 1.5181 <!-- claim: m1_film50_phi_T --> and
  1.5699 <!-- claim: m2_film50_phi_T --> times the true k_L A, beyond the `40 %` band;
- what can flip the verdicts: at the review's ionic k_L floor (`1e-4 m/s`) the shift fits nowhere
  at `20 %` (W reaches only 17.266 <!-- claim: box_m1_tol20_kL1e-4_W_max --> and
  8.5519 <!-- claim: box_m2_tol20_kL1e-4_W_max -->); the whole box has W of at least 25 at `20 %`
  only for a diffusivity below 1.6184e-09 <!-- claim: box_m1_tol20_kL4e-4_D_all_W_ge_25 --> m2/s
  (first order) and 5.4384e-10 <!-- claim: box_m2_tol20_kL4e-4_D_all_W_ge_25 --> m2/s (second);
  surface renewal moves the upper edge to 0.70451 <!-- claim: lim_m1_tol20_upper_edge_ratio_sr_over_film -->
  of the film-theory edge; the nil-bulk verdict flips with Al: for the first order the film
  carries most of the reaction at the nil point only below Al = 45.416 <!-- claim: nil_flip_Al_m1 -->,
  which the box reaches at the review's k_L floor at 1 <!-- claim: nil_film_majority_m1_kL_floor_count -->
  of its two corners, and for the second order only up to Al = 2764.6 <!-- claim: nil_flip_Al_m2 -->;
- the printed sulfite numbers: 100 g/L is 0.7934 <!-- claim: printed_na2so3_molar_100gL --> M,
  consistent with Figure 4's `0.79` M; "1 to 2 M" matches 10-100 g/L only read as an ionic
  strength; the pseudo-m-th-order rate holds at the band edges (E_inf about
  17.529 <!-- claim: E_inf_10gL --> at `10 g/L`, under assumed C* and diffusivity ratio).

Only the kinetic part of "more restricted" is tested. The page does not say that any published
sulfite experiment sat inside or outside the band: the review prints no k_R, catalyst level,
partial pressure or temperature. The 25 is the effect of "5% glycerin" (Onken and Schalk), not the
uncertainty of k_R in a clean experiment. Ion effects on A and k_L, gas-phase depletion and
scale-up are not modelled.

## Files

- `index.py` / `index.ipynb`: the page (jupytext pair).
- `model.py`: eq. 9 and Al, the closed form for m = 1, `FilmTank` (the film coupled to the tank
  as one bordered pymrm system), the root and Richardson helpers, and `SulfiteTank` (the studies
  over the box from one dictionary of constants).
- `test_model.py`: fast checks: pymrm against the closed form at m = 1, the oxygen balance, one
  check value, one stored claim.
- `data/vant-riet-1979-sulfite-statements.csv`: the review's sulfite sentences, verbatim, with
  pages (sidecar `.meta.yaml`).
- `data/check_values.csv`: independent values from a second route (no pymrm), asserted before
  the first claim.
- Eq. 9, k_L, the accuracy and the Figure 4 legend are loaded from page `A3.7`'s CSVs.

## Rebuilding

```
jupytext --sync index.py
python scripts/run_kernels.py pages/F2.6-stirred-tank-sulfite-reactor
```
