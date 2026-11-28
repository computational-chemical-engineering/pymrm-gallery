# H1.6: Gallucci et al. (2013), membrane or bed-to-wall film, which limits the hydrogen flux

F. Gallucci, E. Fernandez, P. Corengia and M. van Sint Annaland, "Recent advances on membranes
and membrane reactors for hydrogen production", *Chem. Eng. Sci.* **92**, 40-66 (2013),
doi:10.1016/j.ces.2013.01.008. A review: it prints no model. On pp. 54 to 55 it states that thin
Pd membranes made the bed-to-wall film the limiting step, backed by Fig. 15 (from Tiemersma et
al., 2006) and Fig. 16 (from Peters et al., 2008). The page also uses the effectiveness factor
the review reports for fluidized beds (p. 57) and Table 2's permeances (p. 53). This page puts those statements on one scale, the flux elasticity of a
Sieverts membrane behind a film (not the share of the pressure drop), decides each from the
review's own numbers where they suffice, and solves a membrane module with pymrm for the
statement on membrane area.

- **Structure:** `S7` (retentate and permeate coupled by the membrane flux), tier T1, data tier
  4 (figure reads and printed numbers; nothing measured is fitted).

## Findings

- **Fig. 16 decides "the limiting factor" (S2, p. 55) without any unprinted number.** From the
  ratio of the film and membrane shares alone, the bound on membrane control is
  0.22335 <!-- claim: r3_bound_10 --> (10 bar), 0.22936 <!-- claim: r3_bound_15 --> (15 bar),
  0.22847 <!-- claim: r3_bound_20 --> (20 bar) and 0.20565 <!-- claim: r3_bound_26 --> (26 bar),
  for any Pd/support interface pressure and support share. Below one half at
  4 <!-- claim: r3_n_below_half --> of the four pressures, the closest
  23.8 <!-- claim: r3_margin_u_20 --> propagated read uncertainties away. A bound on membrane control does not by itself make the
  film limiting, since the support could take the control; the film itself keeps at least
  0.7826 <!-- claim: r3_epsk_min_15 --> of it at 15 bar (the lowest of the four) when the permeate
  H2 pressure is not negative. The device is a
  permeation cell with sweep gas (Caravella et al., 2009, p. 623): "bed-to-wall" is the
  review's label.
- **Fig. 15 cannot decide S1 on its own.** The radial reads are parabolas
  (RMS residual 0.19563 <!-- claim: r1_rms_u --> read units, all
  27 <!-- claim: r1_n_within_u --> points within one). Under vacuum permeate the membrane keeps
  0.81176 <!-- claim: r2_epsQ_vac_2D_4x --> of the control or more even at four times the
  permeability, so the film limits only through the permeate pressure. The thresholds for "a limiting step":
  at the experimental permeability the film is the limiting step only if the permeate H2 pressure
  exceeds 0.84119 <!-- claim: r2_theta_star_2D_1x --> of the bed's (centre reference; read
  uncertainty ±0.010658 <!-- claim: r2_theta_star_u_2D_1x -->, 0.83944 <!-- claim: r2_theta_star_n1_2D_1x -->
  at Sieverts exponent 1) or 0.91682 <!-- claim: r2_theta_star_1D_1x --> (cup-mix;
  ±0.0058526 <!-- claim: r2_theta_star_u_1D_1x -->, 0.91636 <!-- claim: r2_theta_star_n1_1D_1x --> at exponent 1);
  at four times, if it exceeds 0.40305 <!-- claim: r2_theta_star_2D_4x --> (±0.0093406 <!-- claim: r2_theta_star_u_2D_4x -->, 0.36631 <!-- claim: r2_theta_star_n1_2D_4x --> at exponent 1) or
  0.63443 <!-- claim: r2_theta_star_1D_4x --> (±0.0071952 <!-- claim: r2_theta_star_u_1D_4x -->,
  0.62351 <!-- claim: r2_theta_star_n1_1D_4x --> at exponent 1); each ± is one read uncertainty on
  the last marker. The review prints no permeate pressure; at a fixed centre value the three
  curves give a control of 0.99033 <!-- claim: r2_fd_control_1x_4x --> ± 0.050544 <!-- claim: r2_fd_control_u --> (one read uncertainty) between the experimental and the quadrupled permeability.
- **Area (S3), pymrm module.** At the review's WGS feed (p. 46) with a permeate of pure H2 at
  one bar (assumed), the area elasticity falls to one half, for a recovery of 0.8 <!-- claim: tc_r_star -->, at a conductance ratio of
  1.1403 <!-- claim: r4_beta_half_pymrm --> (pymrm, grid error of the finest-grid value
  0.00056223 <!-- claim: r4_beta_half_pymrm_est -->; quadrature
  1.1403 <!-- claim: r4_beta_half_quad -->); 0.95425 <!-- claim: r4_beta_half_R80_vac --> with a
  vacuum permeate, and across the recoveries of the page's table 1.2546 <!-- claim: r4_beta_half_R50_tc -->
  at the lowest to 1.0522 <!-- claim: r4_beta_half_R90_tc --> at the highest: above it a better membrane buys less area than a
  better film. At the feed, with equal conductances, the film takes
  0.5033 <!-- claim: r4_share_in --> of the drop while the membrane keeps
  0.59292 <!-- claim: r4_epsQ_in --> of the control: half the drop is not half the control. At equal conductances the area is
  1.7095 <!-- claim: r4_N80_over_mem --> times the membrane-limit area, so neglecting the film does
  not follow from membrane control.
- **Table 2.** Permeance times thickness equals permeability in
  22 <!-- claim: r5_n_consistent_two --> of 23 <!-- claim: r5_n_rows --> rows when every entry carries at least two
  significant figures; ZIF-8/Al2O3 gives 0.8 <!-- claim: r5_zif8_ratio --> of its printed
  permeability. With every digit significant 20 <!-- claim: r5_n_consistent_strict --> of 23 <!-- claim: r5_n_rows -->;
  with every trailing zero a placeholder 23 <!-- claim: r5_n_consistent_one --> of 23 <!-- claim: r5_n_rows -->. The
  SINTEF permeance in the text (p. 46) is 3.078 <!-- claim: r5_sintef_ratio_secant --> or
  1.539 <!-- claim: r5_sintef_ratio_tangent --> times Table 2's at Sieverts exponent one half
  (two readings of the 1 bar footnote), outside the printed digits under both; which number is
  wrong is not decided.
- **p. 57.** An effectiveness factor "higher than 0.9" read as a series film gives a film elasticity of
  0.10497 <!-- claim: r6_epsk_vac_lo --> (vacuum permeate) to 0.1 <!-- claim: r6_epsk_y1_lo -->;
  "as high as 0.98" gives 0.0202 <!-- claim: r6_epsk_vac_hi -->.

## Things to know before reusing this

- The review prints no model; the equations are the page's standard resistances-in-series
  relations. The figure verdicts are closed forms and reads; pymrm solves only the module.
- Fig. 15's ordinate is a weight-fraction ratio, used as a partial-pressure ratio; the bounds and
  thresholds of Fig. 15 are worded so that a smaller depletion only strengthens them. Markers 1 to 4 of the two triangle
  curves are hidden behind the squares and not read.
- Fig. 15, Fig. 16 and p. 57 are different membranes and devices; the one-scale table places
  them, it does not compare packed and fluidized beds.
- Not covered: reaction, pressure drop and particle size, temperature control, sweep gas,
  Tables 1, 3, 4 and 5. The originals (Tiemersma 2006, Peters 2008 and 2011a, Chen 2007) were
  not available; nothing here comes from them.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | the local series problem and its elasticities, the figure conversions, the module by quadrature and closed forms, and `MembraneModule` (pymrm finite volumes, area as a bordering unknown, discrete dN/dbeta) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/gallucci-2013-fig15-reads.csv` | figure reads of Fig. 15, with a sidecar stating the method and read uncertainty |
| `data/gallucci-2013-fig16-reads.csv` | figure reads of Fig. 16's segment boundaries, with a sidecar |
| `data/gallucci-2013-table2.csv` | Table 2 as printed strings, with a sidecar |
| `data/gallucci-2013-printed-numbers.csv` | the printed numbers of p. 46, 53 and 57, as strings, with a sidecar |
| `data/check_values.csv` | values of the page's functions from a second route, asserted before any claim |
