# C2.5 Methanol synthesis from CO2 (Vanden Bussche and Froment)

Vanden Bussche and Froment (1996), *J. Catal.* 161, 1-10,
doi:10.1006/jcat.1996.0156.

Methanol from CO2 only, CO entering through the reverse water gas shift: eq. 3
with the parameters of Table 2 and two printed equilibrium constants. The
experiments behind it ("276 experiments", p. 5) are not printed, so the page does not refit them. It asks
whether the printed model reproduces the paper's own simulations at the
conditions the paper states.

- **Structures:** `S2` (plug flow with reaction)
- **Files:** `index.py` (the page, paired with `index.ipynb`), `model.py`
  (`Kinetics`, `Thermo`, `PlugFlowBed`, the Radau reference, the equilibria),
  `test_model.py`, `extract_figures.py` (regenerates the digitised data from
  your copy of the PDF)

## What the page finds

- **Printed constants.** Table 1 converts to Table 2 within the printed digits
  only for a gas constant between 8.31405 <!-- claim: r1_R_lo --> and
  8.31471 <!-- claim: r1_R_hi --> J/(mol K); with R rounded to three decimals the water
  constant misses by 1.81 <!-- claim: r1_hu_KH2O_R8314 --> half-units (1.24 <!-- claim: r1_hu_KH2O_R8314_min --> at the smallest,
  allowing for the rounding of the printed A*, B and T_av). Table 1's intervals use
  the multiplier 1.96 <!-- claim: r1_multiplier --> (a normal quantile, while
  the two-sided t quantile for 543 <!-- claim: r1_dof --> degrees of freedom is
  1.9643 <!-- claim: r1_t975 -->). Back-transformed, the
  frequency factor of the k5 lump has a t ratio of at most
  0.701 <!-- claim: r1_tA_max_k5 --> for any correlation: besides the water
  adsorption constant the paper names, it is a second factor that is not
  significant.
- **Fig. 5 (bench reactor).** With 298 K heats, argon and the rates as printed,
  the bed ends at 550.92 <!-- claim: r3_T_plateau --> K,
  0.33 <!-- claim: r3_plateau_minus_fig --> K from Fig. 5b's plateau; the van 't Hoff
  heats end at 553.21 <!-- claim: r3_T_plateau_HvH_Ar --> K. The RWGS reverses at
  2.89 <!-- claim: r3_zrev_mm --> mm (the text: "3 mm", p. 7). The rate factor that fits
  Fig. 5b's shape best is 1.096 <!-- claim: r3_F -->, but the heat-capacity
  perturbation set before the run (the heat capacities are not printed) moves it by 1.05 <!-- claim: r3_shift_over_half_distance -->
  times half its distance from 1: by the rule set before the run, the factor is
  undecided, though every reading tested gives a factor above 1 (lowest
  1.046 <!-- claim: r3_F_lowest -->).
- **Figs. 6 and 7**, predicted without refitting: methanol plateaus within
  0.033 <!-- claim: r4_fig7b_dyM_maxabs --> mol% of Fig. 7b; Fig. 7a's
  temperatures sit 2.8 <!-- claim: r4_fig7a_offset_min --> to
  3.7 <!-- claim: r4_fig7a_offset_max --> K above the model, and its base case
  3.3 <!-- claim: r4_fig7a_minus_fig5b --> K above Fig. 5b's: an offset axis.
- **Fig. 8 (Klier's maximum).** At the caption's condition the curve
  exceeds the equilibrium allowed by the paper's own constants (carbon basis,
  ideal gas) by up to 8.3 <!-- claim: r5_excess_max --> points, at
  57 <!-- claim: r5_points_above_eq --> of 114 <!-- claim: r5_points_total -->
  points by more than three pixels. Eq. 3 is written in partial pressures, so
  it cannot reach that curve there. At 225 <!-- claim: r5_t_alt_C --> C the curve lies below equilibrium
  (0 <!-- claim: r5_points_above_eq_alt --> of the points above it); on a
  CO-only basis equilibrium reaches it (0 <!-- claim: r5_points_above_eq_co -->
  of the points above by more than three pixels), but eq. 3 then fits with an RMS of
  at best 23.7 <!-- claim: r5_co_rms_caption --> points at the caption's temperature. Eq. 3
  reproduces the whole traced line at 504.2 <!-- claim: r5_T_fit --> K and
  W/F0 = 10.5 <!-- claim: r5_tau_fit --> kg s/mol (RMS
  0.70 <!-- claim: r5_rms_fit --> points, largest residual
  1.58 <!-- claim: r5_resid_max --> points); at the caption
  temperature the best RMS is 5.04 <!-- claim: r5_rms_caption --> points. With
  the water terms of eq. 3 removed, the conversion at the curve's highest CO2 rises from
  21.2 <!-- claim: r5_model_X21 --> (grid error -0.08 <!-- claim: r5_X21_pymrm_minus_ivp --> points against a Radau solve) to 37.1 <!-- claim: r5_removed_X21 --> %,
  the equilibrium value.

Every figure compared is the authors' own simulation, so agreement tests the
transcription and the unprinted inputs, not the chemistry.

## Data

Tables 1 to 3 and the printed numbers are vector text in the PDF; each entry was
checked against the page image and is read as text. Figs. 5 to 8 are 1-bit
scans, digitised by `extract_figures.py` at native resolution and calibrated on
their tick marks; each CSV has a provenance sidecar. Heats and heat capacities
come from the NIST Chemistry WebBook (sidecars). `check_values.csv` holds the
values of an independent check script.
