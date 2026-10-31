# D2.2: Parametric sensitivity and runaway in a fixed bed

Two 1970 criteria for when a cooled tubular reactor stops being controllable,
rebuilt and then swept over the whole operating plane the original could only
sample.

- **Structures:** `S2` (plug flow with reaction)
- **Reference:** Van Welsenaere & Froment (1970), Chem Eng Sci 25(10) 1503–1516,
  doi:10.1016/0009-2509(70)85073-4
- **Runtime:** about 25 s
- **Source:** `index.py` (paired with `index.ipynb`); model code in `model.py`
  (`VWFCase`, the criteria `t_max_peak`, `t_tangency`, `inlet_limits`, the
  phase-plane `back_integrate`, the pymrm reactor `RunawayTube` with
  `trace_inlet` for continuation, and `critical_p0`), tested by `test_model.py`

## Agreement

**0.054 <!-- claim: examples_mean_dev_pct --> % mean absolute deviation over all
30 <!-- claim: n_examples --> numbers printed in their Section 6**, worst
0.40 <!-- claim: examples_worst_dev_pct --> %, nothing fitted. Split by what they
test: 0.051 <!-- claim: criterion_mean_dev_pct --> % over the 25 closed-form
criterion values, 0.069 <!-- claim: integration_mean_dev_pct --> % over the 5
numerical-integration values.

Separately, the pymrm reactor solve and an independent phase-plane
back-integration agree on the critical inlet pressure to
**0.18 <!-- claim: cross_method_worst_pct --> %** over 600–700 K, the first-order
upwind error at 1500 cells. Those two numbers are the only pair on the page that
does not involve the paper.

## Provenance: tier 6, not experimental

The paper contains **no experimental data of any kind.** Every value here is the
authors' own computed result: partly closed-form extrapolation, partly their
fourth-order Runge–Kutta. The page is validated against a published reference
solution in the same sense as `B1.1` and `F3.1`. Do not describe it as validated
against experiment.

## Two traps

1. **`c_p` is printed as kcal/m³·°C: it is already volumetric.** Multiplying by
   `rho_g` "to fix the units" changes `B` and `C` by that factor and moves the
   runaway boundary with no other symptom. `rho_g` appears only in `A`. The
   break table runs this trap first.
2. **The Elsevier full-text API is useless for the numbers in this paper.** It
   returns the publisher's OCR of the 1970 scan, which discards the mid-dot
   decimal separator: `R = 00125 m` for 0.0125 m, `b = 19837` for 19.837,
   `001353` for 0.01353 atm. Excellent for prose, unusable for parameters.
   Everything numeric here was read from a 600 dpi page render.

The check that the reading is right costs nothing and is in the notebook:
quantities the paper prints (ln K, t_w, T_M and the Example 1 bracket) each
depend on the whole parameter table, and none was used to obtain it.

## What the page adds

Their Figs. 6 and 7 are graphical solutions of one implicit equation (Eq. 20),
included explicitly so no computer was needed. Solving it directly turns both
into checks. Their Fig. 8, recomputed by the same back-integration at 21 wall
temperatures, shows that the second criterion has the narrower bracket and the
more conservative critical value, as claimed, but that the **first**
criterion's midpoint is the better estimate at every wall temperature. The paper
compares them only at 625 K, where the midpoints are already
0.8 <!-- claim: mean_error_c1_625K_pct --> % and
3.2 <!-- claim: mean_error_c2_625K_pct --> % above their exact values.

The pymrm reactor confirms the criterion-1 values, and `pymrm.continuation`
draws the hot spot through the runaway region: at 625 K, 1 % more inlet pressure
from 0.0182 atm raises it by 304 <!-- claim: dT_hotspot_1pct_converged_K --> K
(converged in the grid). The steepest relative sensitivity lies about
11 <!-- claim: steepest_sensitivity_above_crit_pct --> % above criterion 1's critical
inlet pressure: the criterion marks the start of the sensitive region, not the
cliff. On the default 1500 cells the steep part is under-resolved by a few kelvin
(3.1 <!-- claim: hotspot_grid_err_0182_K --> K at 0.0182 atm).

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/D2.2-van-welsenaere-froment-runaway/index.py   # after editing index.py
python scripts/run_kernels.py pages/D2.2-van-welsenaere-froment-runaway
GALLERY_OFFLINE=1 python -m pytest -q pages/D2.2-van-welsenaere-froment-runaway/test_model.py
```

No PDF is needed to rebuild: both datasets are committed CSVs.

## Corrections

- 2026-10-30: stated that past the critical inlet pressure the reactor has no
  bounded hot spot and that a failed Newton solve is the signal of runaway; the
  solution stays bounded (1165 <!-- claim: hotspot_p0190_K --> K at 0.019 atm and
  625 K on 1500 cells, 1167 <!-- claim: hotspot_0190_converged_K --> K converged,
  below the adiabatic 1417 K) and Newton converges to it from a neighbouring
  solution.
- 2026-10-30: stated that the continuous Fig. 8 boundary came from bisection on
  reactor solves; it is the phase-plane back-integration, and the pymrm reactor
  checks criterion 1 at six wall temperatures.
- 2026-10-30: stated that the two criteria are indistinguishable at 625 K; there
  the first criterion's midpoint is 0.8 <!-- claim: mean_error_c1_625K_pct --> % and
  the second's 3.2 <!-- claim: mean_error_c2_625K_pct --> % above their exact values,
  and the exact critical pressures differ by
  0.9 <!-- claim: c2_below_c1_625K_pct --> %.
- 2026-10-30: stated that following the inlet pressure by continuation makes the
  critical value depend on the path; continuation agrees with independent solves
  to below 1e-6 <!-- claim: trace_vs_independent_K --> K at every traced point.
- 2026-10-30: printed the hot-spot grid error "at n_z = 1500" from the
  n_z = 1600 solve (3.8e-2 K); at n_z = 1500 it is
  4.0e-2 <!-- claim: hotspot_error_n1500_K --> K.
