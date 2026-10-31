# A3.1: Whitman's two-film theory, and what his three runs can actually test

Two films in series, equilibrium at the interface between them, and the
absorption rate is whatever both films can pass at once. The result is
algebraic, which makes it easy to "validate" by rearranging it. This page
reproduces Whitman's worked example, rebuilds the films as a pymrm diffusion
problem, and measures what his three runs can and cannot test.

- **Structure:** `S3` (1-D steady BVP, two domains coupled by a nonlinear interface)
- **Reference (origin, not consulted):** Whitman (1923), *Chemical and
  Metallurgical Engineering* **29**, 146–148. Pre-DOI and unreachable; volume
  and pages inherited from the reprint's header, no issue number given.
- **Read from:** the verbatim reprint in *Int. J. Heat Mass Transfer* **5**(5)
  429–433 (1962), [doi:10.1016/0017-9310(62)90032-7](https://doi.org/10.1016/0017-9310(62)90032-7),
  rendered at the scan's native 300 ppi, every numeric cropped and re-read; the
  OCR text layer is not used for any number. All page references are to the
  reprint. The 1962 Editor's Foreword (signed "D.B.S.", quoting Sherwood and
  Pigford 1951) is editorial commentary and is never attributed to Whitman.
- **Data:** Table 1 (tier 3, three runs) and every number printed outside it
  (tier 6). No figure is digitised; no other page's dataset is used.

## How it is solved

`TwoFilm` (`model.py`) puts the gas film `p`, the interface values `p2`, `c2`
and the liquid film `c` on one `pymrm.BlockLayout`. The interface values enter
the films through `shapes_d`, so every operator is constant; the equilibrium
row `p2 = f(c2)` is the one nonlinear row, with its derivative placed
analytically; the interface faces are named subdomains selected with `pick`.
Newton starts from the deleted-liquid-film state with `scale=` per block, and
each solve is judged by its backward error. The worked example, the
model fits and the operating-plane sweeps need only the scalar root of eq. (4)
(`two_film_algebraic`); pymrm's part is the coupled finite-volume film, checked
against the exact reacting-film profile, and the BlockLayout coupling pattern. The equilibrium curve is
reconstructed from the four equilibrium pairs that Whitman's own subscript
convention puts in Table 1, not digitised from Fig. 3.

## Results

| Check | Result |
|---|---|
| four equilibrium pairs on one semi-log line | worst residual 0.383<!-- claim: eq_curve_worst_resid_pct --> % in p |
| Whitman's held-out run-2 rate, from Table 1 alone | 0.558 % <!-- claim: run2_rate_reldev --> low |
| same, liquid film deleted (closed form) | 0.054 % <!-- claim: run2_rate_reldev_no_liquid_film --> high: closer |
| gas-film share of run 2's resistance | 99.4 % <!-- claim: run2_gasfilm_resistance_share --> |
| three runs, rms: gas film / liquid film / two films | 41.7<!-- claim: gasfilm_3run_rms_pct --> % / 59.5<!-- claim: liquidfilm_3run_rms_pct --> % / 0.379<!-- claim: twofilm_3run_rms_pct --> % |
| the one degree of freedom (run3/run2): two films / gas film alone | +0.934<!-- claim: dof_ratio_overpredict_pct --> % / +0.319<!-- claim: gasfilm_dof_ratio_dev_pct --> % |
| k_c admitted at 2 % rms by every subset with run 1 | 1.126 <!-- claim: kc_band_low_run1_subsets --> to 1.244 <!-- claim: kc_band_high_run1_subsets --> |
| k_c admitted with run 1 removed | from 0.396 <!-- claim: kc_band_low_without_run1 -->, no upper end |
| 1/K_p = 1/k_p + m/k_c on run 1, tangent / secant m | 34.1<!-- claim: constH_run1_tangent_err_pct --> % / 12.8<!-- claim: constH_run1_secant_err_pct --> % |
| K_p, K_c over the operating window at fixed k_p, k_c | factors 3.80 <!-- claim: Kp_range_factor --> and 13.59 <!-- claim: Kc_range_factor --> |
| reacting liquid film vs the exact film profile, 200 cells | 9.9e-10 <!-- claim: reaction_vs_exact_film_n200 -->, second order |

The pymrm solve reproduces the algebra to round-off by construction (a linear
profile is exact in the finite-volume scheme); the page labels it a port check.
With a fast reaction in the liquid film the rate approaches
`k_p (p1 - f(0))`, the gap to `k_p p1` being the intercept of the reconstructed
line, 2.96e-6 <!-- claim: intercept_over_p1 --> of it.

The page cannot choose between film, penetration and surface renewal; that
comparison is on [`A3.2`](../A3.2-higbie-penetration/). The reaction cell is a
limit check; [`F3.1`](../F3.1-hatta-regimes/) owns enhancement.

## Corrections

- 2026-10-31: wrote the linearised overall coefficient as
  1/K_p = 1/k_p + 1/(m k_c), which is dimensionally inconsistent, and published
  its errors on run 1 as 87.8 % (tangent m) and 114.9 % (secant m), 29.1 % and
  15.9 % with the slope halved, and runs 2 and 3 "saturating at -100 %". For
  p = m c the formula is 1/K_p = 1/k_p + m/k_c: run 1 is off by 34.1 % and
  12.8 %, 14.3 % and 4.9 % with the slope halved, and the gas-controlled runs
  2 and 3 by at most 0.17<!-- claim: constH_runs23_tangent_worst_err_pct --> %
  with the tangent. With a linear equilibrium both choices of m are exact.
- 2026-10-31: called the 3.0e-6 gap of the fast-reaction rate to k_p p1 a grid
  floor of an unresolved reaction layer, not convergence, and said the rate
  rises towards k_p p1. The gap converges with the grid to the exact
  reacting-film value 2.969e-6 <!-- claim: reaction_fast_exact_film_resid -->
  at ka = 1e6, and as ka grows tends to f(0)/p1, the intercept of the
  reconstructed equilibrium line. Removing the intercept moves the exact-film value
  to about 1e-8;
  refining from 25 to 200 cells moves it by 4 %.
- 2026-10-31: gave the k_c band of every subset containing run 1 as
  [1.14, 1.24], and runs 2 and 3 alone as [0.40, 3000]. The edges were sweep
  samples; root-found they are [1.126, 1.244], and without run 1 the band starts
  at 0.396 and has no upper end (3000 was the top of the sweep).
- 2026-10-31: said 3 and 600 cells give "the same twelve figures". They agree
  to below 1e-9 relative (round-off), about ten figures.
- 2026-10-31: said the two routes to the deleted-liquid-film rate share no
  arithmetic. They share the equilibrium function and no solver.

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        Equilibrium, two_film_algebraic, no_liquid_film_rate,
                reacting_film_rate, TwoFilm (BlockLayout, shapes_d, newton)
test_model.py   fast checks of model.py and two stored claims
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
data/whitman-1923-table1.csv              the ten printed columns of Table 1
data/whitman-1923-printed-results.csv     the worked example, stated ratios, Fig. 3 ticks
```

Runtime about 9 s.
