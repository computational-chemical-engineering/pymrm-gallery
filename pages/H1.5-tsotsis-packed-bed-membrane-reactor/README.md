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

The paper measured its kinetics and the permeability of all species involved in the reactions
(p. 2906) but prints neither, and its eq. 1 carries a reaction in the wall of a Pt-impregnated
tube; the page drops that term (A6) and leaves Pi and Da free.

- Fig. 3b is not reproduced by the printed model without its wall reaction. One bound holds for
  all eight closure variants (pressure read as gauge or absolute, tabulated or the paper's own
  equilibrium constant, reaction in the tube only or also in the shell; both scales proportional
  to residence time, the reaction at its fast limit): the model's highest supremum (the tube at
  drainage) is 10.07 <!-- claim: r4_max_sup_over_variants_pct --> %, 1.07 <!-- claim: r4_min_gap_10s_pt -->
  points short of the printed conversion at 10 s, 12.7 <!-- claim: r4_min_gap_10s_re --> times that
  node's reading error. The best variant (absolute pressures, tabulated equilibrium constant,
  reaction in the shell too) misses the printed model line by up to
  1.08 <!-- claim: r4_best_maxdev_pt --> points (at 10 <!-- claim: r4_best_maxdev_pt_node_s --> s,
  12.9 <!-- claim: r4_best_maxdev_pt_re --> of that node's reading errors) and by up to
  16 <!-- claim: r4_best_maxdev_re --> reading errors (at 9 <!-- claim: r4_best_maxdev_re_node_s --> s).
  On the gauge reading the best variant (the paper's own equilibrium constant, reaction in the shell
  too) misses by 3.07 <!-- claim: r4_g_pk_ts_maxdev_pt --> points; its supremum is
  8.07 <!-- claim: r4_g_pk_ts_sup_pct --> %.
- The model rises with falling slope in every variant: 0 <!-- claim: r4_n_convex_windows --> of
  1600 <!-- claim: r4_n_windows --> windows (second difference at the 2, 6 and 10 s nodes, 200
  values of Pi per variant) are convex. In the best variant the second difference over 2, 6 and
  10 s is -0.83 <!-- claim: r4_best_d2_wide_pt --> points, where the printed line is convex
  (0.81 <!-- claim: r4_printed_d2_wide_pt --> points,
  5.6 <!-- claim: r4_printed_d2_wide_sig --> reading errors). That conversion rises with Da, which
  makes the fast-reaction limit the envelope, was checked on the gauge tube-only variant
  (72 <!-- claim: r3_map_n_rise_with_Da --> of 72 <!-- claim: r3_map_n_pairs --> map pairs) and on the
  variants a_a7_ts, a_a7_t and g_pk_ts (135 <!-- claim: r4_var_n_rise_with_Da --> of
  135 <!-- claim: r4_var_n_pairs --> pairs).
- The cap is drainage: at a fixed pressure difference the Knudsen wall passes ethane too, and the
  tube is drained (its outflow at the page's drainage threshold) at a permeation number of 2.82 <!-- claim: r3_Pi_d_DaT -->
  (gauge reading, tube-only reaction, test rate constant) and 1.91 <!-- claim: r4_best_Pi_d --> in
  the best variant (absolute reading, shell reaction, fast reaction). With the tube drained, its
  whole content reaches the shell and equilibrates there, so the supremum of the best variant,
  10.07 <!-- claim: r4_a_a7_ts_sup_pct --> %, is the sweep-diluted equilibrium at the shell pressure
  (the two differ by 0.00018 <!-- claim: r4_sup_minus_xeq_mix_pt --> points). That equilibrium
  reaches the printed conversion at 10 s only at a sweep ratio of
  1.36 <!-- claim: r4_Pr_needed_best -->, against the printed sweep ratio. The finite-volume solve (pymrm)
  confirms the best variant at its 6 s node: observed orders 1.01 <!-- claim: r4_fv_best_order_min -->
  to 1.02 <!-- claim: r4_fv_best_order_max -->, extrapolated conversion
  9.3475 <!-- claim: r4_fv_best_extrap_pct --> %, within 0.008 <!-- claim: r4_fv_best_gap_over_err -->
  of its extrapolation-error estimate of the Radau value; the test and the suprema themselves run
  on the local-equilibrium integration and Radau.
- The equilibrium lines favour the gauge reading; the page carries both. On the dilute lines of
  Fig. 2a the rms deviation is 0.65 <!-- claim: r1_rms_2a_gauge --> points against
  5.66 <!-- claim: r1_rms_2a_abs --> on the absolute reading; the dense lines do not separate the
  readings (rms 1.19 <!-- claim: r1_rms_dense_gauge --> against 1.41 <!-- claim: r1_rms_dense_abs -->).
  The dense-feed lines lie above what the printed dense feed gives (mean residual, computed minus
  read, -1.14 <!-- claim: r1_mean_dense_gauge --> points); an H2 fraction of
  4.47 <!-- claim: r1_fit_yH_pct --> % in place of the printed one, or K_p times
  1.37 <!-- claim: r1_fit_Kmult -->, lowers the dense rms from 1.19 <!-- claim: r1_rms_dense_gauge -->
  to 0.28 <!-- claim: r1_fit_yH_rms --> or 0.35 <!-- claim: r1_fit_Kmult_rms --> points; both are
  fits, not the authors' inputs. Two printed copies of one shell-side line disagree by
  0.61 <!-- claim: r1_m3_gap_pt --> points.
- On the gauge reading all 29 <!-- claim: r2_n_symbols --> digitised measured points lie above the
  shell-side equilibrium computed for the printed feeds by more than two reading errors
  (29 <!-- claim: r2_n_above_shell -->), 29 <!-- claim: r2_n_above_tube --> above the tube-side one, and
  27 <!-- claim: r2_n_above_mix --> above the sweep-diluted equilibrium, the reference a wall that
  passes every gas must beat; the smallest margin of that count is
  0.0015 <!-- claim: r2_mix_min_margin_pt --> points from its two-reading-error threshold. On the
  absolute reading the counts are 25 <!-- claim: r2_n_above_shell_abs -->,
  29 <!-- claim: r2_n_above_tube_abs --> and 2 <!-- claim: r2_n_above_mix_abs -->.

## Things to know before reusing this

- The paper defines neither conversion, sweep ratio, residence time nor its two equilibrium lines;
  each reading is a named assumption on the page, and the missing ingredient may be in the thesis
  the paper cites. The page does not say the paper's curves are wrong.
- No rate constant, permeance, layer property or porosity on the page is the paper's; the page
  predicts no physical reactor.
- Not covered: the energy balances, the pressure drop (the paper adds Ergun-type equations, p. 2905;
  A8 neglects it), reaction inside the membrane, any steam reforming model; the Fig. 4 measured
  points are not digitised.
