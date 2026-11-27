# D3.5 Ethylene oxide multitubular: is the kinetic window parametrically sensitive?

Borman, P. C. and Westerterp, K. R. (1995), *Ind. Eng. Chem. Res.* 34, 49-58,
doi:10.1021/ie00040a002. Criterion: Morbidelli, M. and Varma, A. (1988), *Chem. Eng. Sci.*
43, 91-102.

The paper measures the kinetics of ethene epoxidation and combustion over silver for work on
wall-cooled packed beds, but prints no tube and no reactor data. This page puts its model III
kinetics (rebuilt on `C2.11`) into a 1-D wall-cooled plug-flow tube, with the tube diameter,
wall coefficient and heat capacity of `C2.10` as named page choices, and tests the paper's own
window (air, ethene up to the lower explosion limit, coolant over the kinetic temperatures,
the standard pressure of Table 4) with Morbidelli and Varma's generalized criterion:

- with Table 6 as printed (combustion factor 1) the window holds no sensitive point: the
  largest normalised sensitivity of the hot spot to the cooling is
  0.0634 <!-- claim: max_abs_Sh_M_F1 --> at 500 K and the explosion limit, and
  0.0770 <!-- claim: max_abs_Sh_window_F1 --> at the worst corner of the window;
- with the factor `F_fig6` that reproduces the authors' Figure 6 (`C2.11`) the boundary runs
  through the window: at the explosion limit and 500 K the critical cooling is
  14.449 <!-- claim: hcrit_MV_M_F6 --> W/(kg K), 0.4326 <!-- claim: hcrit_over_hSem_M_F6 -->
  of the Semenov tangency, with a sensitivity peak of
  -414.8 <!-- claim: smin_MV_M_F6 -->; over the corners of the window the ratio to the
  tangency runs from 0.327 <!-- claim: map_ratio_min --> to 0.493 <!-- claim: map_ratio_max -->;
  the `C2.10` tube is critical for coolant temperatures near
  507.86 <!-- claim: Tc_crit_C210_03_F6 --> K at the explosion limit and
  526.1 <!-- claim: Tc_crit_C210_018_F6 --> K at the standard ethene fraction of Table 4,
  inside the top kelvin of the measured range;
- all 7 <!-- claim: Tstar_crit_above_top --> critical hot spots computed lie above the highest
  kinetic temperature (733.8 <!-- claim: Tstar_crit_M_F6 --> K at the explosion limit and
  500 K), and the hot spot reaches that temperature only at
  1.2675 <!-- claim: h527_over_hcrit_M_F6 --> times the critical cooling: the boundary is
  located by extrapolated kinetics;
- the criterion is intrinsic here (the ethene-fraction criterion at the critical cooling
  returns 0.0300 <!-- claim: ycrit_MV_at_hcrit_M_F6 -->); normalising by the rise instead of
  the hot-spot temperature moves the `C2.10` tube's critical coolant temperature at the
  standard ethene fraction by
  -0.25 <!-- claim: Tc_crit_rise_shift_018_K --> K; with `C2.11`'s Table 4 row 1 factor the
  critical cooling is 7.796 <!-- claim: hcrit_MV_M_Frow1 -->, proportional to the factor
  (logarithmic slope 1.031 <!-- claim: dlnh_dlnF_M -->);
- the hot spot costs selectivity: at the base tube 0.0149 <!-- claim: sel_cost_hotspot_B -->
  at half conversion, near the boundary 0.0663 <!-- claim: sel_cost_hotspot_nc -->.

This supports `C2.11`'s finding that the printed combustion constant is too small from the
reactor side; it does not prove it. Nothing here is compared with reactor or plant data.

## Files

- `index.py` / `index.ipynb`: the page (jupytext pair).
- `model.py`: model III rates and heat release, the inlet, `CooledTube` (the steady FV tube
  with pymrm, its sensitivity and the criterion), and the critical-point march.
- `test_model.py`: fast checks: one check value, the element balances, the frozen-limit root,
  one stored claim.
- `data/`: the printed reactor numbers and statements of the paper and the check values of a
  second route, each with a `.meta.yaml` provenance sidecar. Table 6, Table 3 and the factors
  are loaded from `C2.11`; U, d_t and c_p from `C2.10`.
