# G1.7: Wammes–Westerterp high-pressure trickle-bed hydrodynamics

Wammes, W. J. A. and Westerterp, K. R., *Hydrodynamics in a pressurized
cocurrent gas–liquid trickle-bed reactor*, Chemical Engineering & Technology
**14**(6) 406–413 (1991), doi:10.1002/ceat.270140608 (open access via the
University of Twente repository). The atmospheric companion study is Wammes,
Mechielsen and Westerterp, Chem. Eng. Sci. 46 (1991) 409–417; the CET paper is
the one that derives the high-pressure correlations and is the one built here.

The page closes the paper's two mutually dependent correlations, pressure
gradient (Eq. 8, containing the hold-up) and dynamic hold-up (Eq. 9, containing
the pressure gradient), as a coupled system, with nothing fitted, against 16
digitised hold-up points (Fig. 7) and 20 digitised pressure-gradient points
(Fig. 6). It also solves the paper's plug-flow CO₂-absorption balance (the
model behind every interfacial-area measurement) with pymrm operators and
verifies it against the paper's own closed form, Eq. (3).

- **Source:** `index.py` (paired with `index.ipynb`); model code in `model.py`
  (`OperatingPoints`, `HoldupPressureClosure` with `solve`, `solve_scalar` and
  `sign_changes`, `CO2PlugFlow`), tested by `test_model.py`.
- **Runtime:** a few seconds.
- `data/`: digitised points and transcribed parameters, each with a
  provenance sidecar. **Digitisation reviewed and confirmed** by the maintainer
  on 2026-08-02 against the numbered overlays in `../review/` (git-ignored,
  they contain the page image): every cross on a real marker, one per marker,
  none missed, series correctly told apart by shape, and the four reconstructed
  centres in the touching pairs (fig7 at v_l = 2 and 4 mm/s; fig6 at
  ρ_g ≈ 0.45 and 1.6 kg/m³) checked individually. No CSV row changed.

## Agreement

Deviation convention: (model − measured)/measured. Hold-up mean |dev|
2.9 % <!-- claim: holdup_mean_abs_dev --> (paper's own error: 8 %); pressure
gradient 10.4 % <!-- claim: pgrad_mean_abs_dev --> (paper: 12 %), with a bias of
+10.3 % <!-- claim: pgrad_bias --> that a porosity of
0.3987 <!-- claim: eps_zero_bias --> instead of the printed 0.39 would remove, as
would a sphere diameter of 3.18 <!-- claim: dp_zero_bias_mm --> mm inside the
printed 3 ± 0.5 mm; neither is fitted.
The matched-pair gas-density test: a 5× density change moves the measured
hold-up by 9.2 <!-- claim: gas_density_shift_meas_pct --> %, Eq. (9) predicts
13.1 <!-- claim: gas_density_shift_model_pct --> %.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/G1.7-wammes-high-pressure-trickle-bed/index.py   # after editing index.py
python scripts/run_kernels.py pages/G1.7-wammes-high-pressure-trickle-bed
GALLERY_OFFLINE=1 PYTHONPATH=shared python -m pytest -q pages/G1.7-wammes-high-pressure-trickle-bed/test_model.py
```

## Corrections

- 2026-10-30: read the absorber's outlet CO₂ from the last cell centre and
  stated its conversion error at n = 200 as 1.19e-3, with an Eq. (3) round trip
  of −0.32 %; on the exit face, the value the operator transports, the error is
  9.8e-4 <!-- claim: xi_upwind_err_n200 --> and the round trip
  -0.27 <!-- claim: a_gl_roundtrip_n200_pct --> %.
- 2026-10-30: stated that dP/L scales as eps^-2.6 at fixed beta_t; at fixed
  beta_t the exponent is -3.9 <!-- claim: dpl_eps_exponent_fixed_beta -->, and in
  the coupled solve -4.1 <!-- claim: dpl_eps_exponent_coupled_max --> to
  -5.6 <!-- claim: dpl_eps_exponent_coupled_min -->.
- 2026-10-30: stated that reading a_s per particle volume would raise every
  hold-up prediction by 44 %; 1.44 <!-- claim: as_alt_factor --> is the factor on
  Eq. (9) at a given pressure gradient, and in the coupled solve the hold-ups
  rise by 21 <!-- claim: as_alt_coupled_rise_min_pct --> % to
  41 <!-- claim: as_alt_coupled_rise_max_pct --> %.
- 2026-10-30: stated the deviations of the v_g = 0.24 m/s hold-up series as
  −7 to +1 %; they are -6.8 % <!-- claim: ser024_dev_min --> to
  +1.6 % <!-- claim: ser024_dev_max -->.
