# H1.5: Tsotsis et al. (1992), packed-bed catalytic membrane reactor

T. T. Tsotsis, A. M. Champagnie, S. P. Vasileiadis, Z. D. Ziaka and R. G. Minet, "Packed bed
catalytic membrane reactors", *Chem. Eng. Sci.* **47**(9-11), 2903-2908 (1992),
doi:10.1016/0009-2509(92)87149-K. The paper prints an isothermal packed-bed catalytic membrane
reactor model (eqs. 1, 3, 5: Knudsen permeation of every species through a porous alumina wall),
fits it to ethane dehydrogenation and plots model curves beside its data, but prints no kinetics
and no permeances. The page closes the model with named assumptions and asks whether any
permeation number and any Damkohler number produce the plotted model curves, above all the rise
of conversion with residence time in Fig. 3b. It recomputes every equilibrium line first.

## Findings

- Fig. 3b is not reproduced. Over eight closure variants (pressure read as gauge or absolute,
  tabulated or the paper's own equilibrium constant, reaction in the tube only or also in the
  shell), with both scales proportional to residence time and the reaction at its fast limit,
  the best variant misses the printed model line by 1.08 <!-- claim: r4_best_maxdev_pt --> points
  at its worst node, 16 <!-- claim: r4_best_maxdev_re --> reading errors. The model rises with
  falling slope (second difference over 2, 6 and 10 s: -0.83 <!-- claim: r4_best_d2_wide_pt -->
  points) where the printed line is convex (0.81 <!-- claim: r4_printed_d2_wide_pt --> points,
  5.6 <!-- claim: r4_printed_d2_wide_sig --> reading errors).
- The cap is drainage: at a fixed pressure difference the Knudsen wall passes ethane too, and the
  tube empties into the shell at a permeation number of 2.82 <!-- claim: r3_Pi_d_DaT --> at the
  test rate constant. The best variant's supremum (the tube at drainage, fast reaction) is
  10.07 <!-- claim: r4_a_a7_ts_sup_pct --> %, 1.07 <!-- claim: r4_gap_sup_vs_3b_10s_pt --> points
  below the printed conversion at 10 s.
- The equilibrium lines follow gauge pressures: on the dilute lines of Fig. 2a the rms deviation is
  0.65 <!-- claim: r1_rms_2a_gauge --> points against 5.66 <!-- claim: r1_rms_2a_abs --> on the
  absolute reading. The dense-feed lines lie above what the printed dense feed gives (mean residual,
  computed minus read, -1.14 <!-- claim: r1_mean_dense_gauge --> points); an H2 co-feed of
  4.47 <!-- claim: r1_fit_yH_pct --> % or an equilibrium constant times
  1.37 <!-- claim: r1_fit_Kmult --> closes the gap. Two printed copies of one shell-side line
  disagree by 0.61 <!-- claim: r1_m3_gap_pt --> points.
- All 29 <!-- claim: r2_n_symbols --> digitised measured points lie above the shell-side
  equilibrium computed for the printed feeds by more than two reading errors
  (29 <!-- claim: r2_n_above_shell -->), and
  27 <!-- claim: r2_n_above_mix --> above the sweep-diluted equilibrium, the reference a wall that
  passes every gas must beat.

## Things to know before reusing this

- The paper defines neither conversion, sweep ratio, residence time nor its two equilibrium lines;
  each reading is a named assumption on the page, and the missing ingredient may be in the thesis
  the paper cites. The page does not say the paper's curves are wrong.
- No rate constant, permeance, layer property or porosity on the page is the paper's; the page
  predicts no physical reactor.
- Not covered: the energy balances, reaction inside the membrane, any steam reforming model; the
  Fig. 4 measured points are not digitised.
