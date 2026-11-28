# D3.1: Steam methane reformer tube (Xu and Froment, 1989, Part II)

Is the conversion of an industrial steam reformer limited by the heat flux? Part II's own
model (Part I kinetics, slab pellet with the fitted tortuosity, Table 2's tube) rebuilt with pymrm and
probed: the exit-conversion elasticities to catalyst activity and to heat input under an
imposed flux and an imposed wall temperature, a printed-material audit, the pellet at the
printed 4 m state and the tube against the digitised Figs. 3 to 7.

- **Structures:** `S2` (axial march), `S3` (pellet boundary-value problem), `S8` (nested coupling)
- **Reference:** Xu & Froment (1989), AIChE J 35(1) 97-103, doi:10.1002/aic.690350110;
  kinetics from Part I through page `C2.1`
- **Runtime:** about 5 min (tube runs in a process pool)
- **Files:** `index.py` (the page, paired with `index.ipynb`), `model.py` (`Case`, `Pellet`,
  `Tube`, `run_tube`), `test_model.py`

## Agreement

There are no plant data: every comparison is with the authors' simulation, digitised. The
ratio of the elasticities e_a/e_Q at the end of heating is 0.0082 <!-- claim: R4_flux_ratio -->
under imposed flux and 0.0237 <!-- claim: R4_Tw_ratio --> under imposed wall temperature;
over the robustness sweeps the largest is 0.0127 <!-- claim: R4_sweep_max_ratio -->.
