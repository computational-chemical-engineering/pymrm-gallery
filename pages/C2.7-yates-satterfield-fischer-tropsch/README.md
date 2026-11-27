# C2.7 Fischer-Tropsch on cobalt: the Yates-Satterfield rate law

Yates and Satterfield (1991), *Energy & Fuels* 5, 168-173,
doi:10.1021/ef00025a029.

The cobalt Fischer-Tropsch rate law -R = a P_H2 P_CO/(1 + b P_CO)^2 (eq. 16),
fitted in a slurry CSTR at 220 and 240 C (Table V), which slurry-reactor models
take as given (`F2.3`). The page

- digitises the experimental design (Figures 1 and 2): the partial pressures are
  correlated (r = 0.797 <!-- claim: r1_r_220 --> and 0.773 <!-- claim: r1_r_240 -->),
  where the paper says they are not;
- audits Table V: the t value of b at 220 C is not estimate over standard error
  (32.9 <!-- claim: r2_allow_b220 --> times the rounding allowance), and the
  design that reproduces the printed 240 C standard errors implies a standard
  error of a at 220 C of 14.7 <!-- claim: r2_sd_a220_design -->, a factor
  10.6 <!-- claim: r2_factor_220 --> above the printed one;
- computes the apparent activation energy Table V carries, at fixed partial
  pressures as the authors did: 133.0 <!-- claim: r3_band_lo --> to
  159.4 <!-- claim: r3_band_hi --> kJ/mol over the range both temperatures
  cover, against "The apparent activation energies were 92.7 and 94.5 kJ/mol"
  (p. 172), values Table V reaches only at
  P_CO 0.0420 <!-- claim: r3_root_927 --> and 0.0444 <!-- claim: r3_root_945 -->
  MPa, below every run. With the a-b covariance from the design the smallest
  one-sample z is 9.56 <!-- claim: r3_z_945 -->; without it,
  0.90 <!-- claim: r3_z_945_rho0 -->, so the covariance decides the verdict. Which of
  Table V and the quoted values is wrong the paper does not let one say;
- checks Table V against Figure 8, Table II and Table IV;
- puts eq. 16 in the paper's CSTR and an ideal plug-flow reactor with pymrm
  (CO conversion 0.2048 <!-- claim: r6_cstr_X --> and
  0.2093 <!-- claim: r6_pfr_X --> at the Table II feed).

The product (chain-length) distribution of the catalogue entry is not in this
paper and is not covered.

## Files

- `index.py` / `index.ipynb`: the page (jupytext pair).
- `model.py`: eq. 16, the apparent activation energy and its closed forms, the
  design covariance, the CSTR (pymrm `newton`) and `PlugFlow` (pymrm upwind).
- `test_model.py`: fast checks, including one stored claim.
- `extract_figures.py`: regenerates the digitised CSVs from a local copy of the
  PDF (not in the repository).
- `data/`: the printed numbers, the digitised design of Figures 1 and 2, the
  parity plot of Figure 8 and the check values, each with a `.meta.yaml`
  provenance sidecar.
