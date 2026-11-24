# J5.6: Smith and Ewart (1948), emulsion polymerisation kinetics and particle number

W. V. Smith and R. H. Ewart, "Kinetics of emulsion polymerization", *J. Chem. Phys.* **16**,
592-599 (1948), doi:10.1063/1.1746951. The page asks whether the particle-number law
`N = k (rho/mu)^(2/5) (a_s S)^(3/5)`, its exponents and its bracket `0.37 < k < 0.53` follow exactly
from the paper's two capture assumptions and its Case 2 growth premise, and over what part of
the paper's own (alpha, m) space of eq. 1, solved without the limiting approximations, that
premise (`1/2` radical per particle, independent of particle size) holds. The paper prints no
physical values; the page works in the paper's dimensionless groups.

## Findings

- Every printed Stage I constant follows from the premises, within 1.6244 <!-- claim: r34_hu_max -->
  half-units of its last printed digit; 3 <!-- claim: r34_n_round --> of the
  5 <!-- claim: r34_n_printed --> round correctly. The upper-bound constant is exactly
  `5^(3/5)/(3 (4 pi)^(1/5))` = 0.52774 <!-- claim: r3_k_upper --> (printed `0.53`), the lower one
  0.37040 <!-- claim: r4_k_lower --> (printed `0.370`). The constant of eq. 30 is
  1.1309 <!-- claim: r3_hu_c_y --> half-units from the printed `4.36` and that of t_1
  1.6244 <!-- claim: r4_hu_t1_const --> from the printed `0.650`: both equal the exact value cut, not rounded.
  The root y_1 = -2.1337 <!-- claim: r4_y1 --> holds as printed, but the five printed terms of eq. 29
  give a root of -2.1354 <!-- claim: r4_y1_n4 -->, 1.0851 <!-- claim: r4_hu_y1_n4 --> half-units off.
- The exponents `2/5` and `3/5` follow from the area law alone (a sphere growing at constant volume
  rate, `a ~ (t - tau)^(2/3)`); dimensional analysis leaves one exponent free and the capture rule
  sets only k. The area bound ends Stage I 1.2332 <!-- claim: r4_ratio_t1 --> times later with
  0.70187 <!-- claim: r4_ratio_k --> times the particles.
- A population balance in particle volume, which does not use the paper's integral equation,
  reproduces both bounds: the extrapolated stage end is 1.3587 <!-- claim: r5_upper_t1_ext -->
  (upper, error estimate 0.000052 <!-- claim: r5_upper_t1_err_est -->) and
  1.6755 <!-- claim: r5_lower_t1_ext --> (lower, error estimate 0.000030 <!-- claim: r5_lower_t1_err_est -->)
  in units of t*, against the closed forms 1.3587 <!-- claim: r3_ts1_upper --> and
  1.6755 <!-- claim: r4_ts1_lower -->; the gaps are far inside those error estimates.
- The paper reads `beta > 1` from eq. 10 (p. 595) and expects nbar near one half "if beta is
  sufficiently large". The exact solution of eq. 1 at `m = 0` puts that at beta of
  9.6723 <!-- claim: r1_beta_10pct_m0 --> for `10 %` and 99.667 <!-- claim: r1_beta_1pct_m0 --> for
  `1 %`; at `beta = 1` it is 0.88782 <!-- claim: r2_exact_nbar_beta1 -->. The series (14) is exact to
  second order (the next coefficient, 0.16667 <!-- claim: r2_exact_c3 --> (`1/6`), is not printed), and
  the Appendix construction is better: relative error
  -0.00089038 <!-- claim: r2_appendix_err_beta1 --> at `beta = 1`, against
  -0.061369 <!-- claim: r2_series_err_beta1 --> for the series.
- The Case 2 plateau needs small desorption: the smallest slope d ln nbar/d ln alpha reaches `0.1`
  only for m = k_0 a/k_t up to 0.0057763 <!-- claim: r1_m_plateau_01 -->, so by that threshold Case 2 is
  already gone above it (at `m = 0.3` the smallest slope is 0.45584 <!-- claim: r1_m03_min_slope -->).
  Any dip below `1/2` is gone from `m = 1/2` on: the closed form at `m = 1/2` has none, and none of
  the 8 <!-- claim: r1_scan_n_at_or_above_half --> scanned m at or above `1/2` (of
  25 <!-- claim: r1_scan_n_m --> up to `10`) has a minimum (minima found: 0 <!-- claim: r1_scan_n_min_at_or_above_half -->);
  there the slope falls from `1` to `1/2`, the paper's "transition from case 1 to case 3 without ever
  realizing case 2" (p. 596). The boundary sits at `m = 1/2` because the large-alpha slope is
  `s = 1/2 + (2m - 1)/(4a) + O(a^-2)` with `a = sqrt(8 alpha)`; at `m = 0.3` and `alpha = 1e10` the
  page finds s - 1/2 = -3.5363e-7 <!-- claim: r1_asym_dev -->, within 0.00022 <!-- claim: r1_asym_gap -->
  (relative) of the asymptote.

## Things to know before reusing this

- The closed form of eq. 1 in Bessel functions is later work (Stockmayer 1957, O'Toole 1965), not
  Smith and Ewart's.
- Nothing here says where a real system lies on the map, or that real particle numbers fall inside
  the bracket: the paper's bracket argument (p. 597) is heuristic and untested here.
- Eq. 8 is read with re-entry of desorbed radicals, which the paper does not state.
- Not covered: molecular weight, Case 1 water-phase kinetics (eqs. 3 to 6 and 9), experiments.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | closed form and truncated recursion of eq. 1, the Appendix I construction, exact series, the Stage I constants (Mittag-Leffler series), `StageIPBE` (population balance in particle volume, solved with `pymrm.newton`) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/smith-ewart-1948-printed-numbers.csv` | the printed constants, exponents, coefficients and expressions, as strings, with their provenance sidecar |
| `data/check_values.csv` | values of the page's functions at test conditions from a second route, asserted before any claim |
