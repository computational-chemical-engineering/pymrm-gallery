# H1.3: Caravella, Barbieri and Drioli (2009), concentration polarisation along a laminar Pd-membrane tube

A. Caravella, G. Barbieri and E. Drioli, "Concentration polarization analysis in
self-supported Pd-based membranes", *Sep. Purif. Technol.* **66**, 613-624 (2009),
doi:10.1016/j.seppur.2009.01.008. The paper computes a concentration polarisation coefficient
(CPC, Eq. 13) from a film model at "a generic point of the tube-in-tube membrane purifier"
(p. 615), with the length-averaged Sherwood correlation of Eq. (9), and states that "the
results are valid in every position of the separation equipment" (p. 618). This page solves
what the paper does not, a resolved 2-D laminar tube with a Sieverts wall, and tests the film
along it; it also checks the paper's printed numbers. The catalogue also cites an IJHE (2015)
paper, which was not available: the page is built from the 2009 paper only.

- **Structure:** `S6` (2-D steady, marched along the tube; the catalogue line says S6+S7, but
  the permeate is pure H2 at fixed pressure, so there is no second domain), tier T1, data
  tier 6 (nothing measured is fitted).

## Findings

- **Module total right, local CPC wrong.** At the test condition (2-D depletion
  5.677 <!-- claim: r1_depl_pct --> %, inside the paper's "(4–7%)", p. 623), the film module with
  Eq. (9)'s averaged coefficient recovers 0.10802 <!-- claim: r1_rec_1d --> of the inlet H2
  against 0.10744 <!-- claim: r1_rec_2d --> from the 2-D solve
  (0.53437 <!-- claim: r1_rec_gap_pct --> %), but its outlet CPC is
  0.72028 <!-- claim: r1_cpc_avg_film --> against 0.80081 <!-- claim: r1_cpc_2d -->. A
  local Lévêque coefficient gives 0.80287 <!-- claim: r1_cpc_local_film -->. The paper's
  recipe, the mean of inlet and outlet film fluxes, is
  0.55193 <!-- claim: r1_recipe_gap_pct --> % above the 2-D module mean.
- **Over the Graetz number** (fixed film beta), the module-total error of the averaged
  coefficient is not monotone: most negative -5.3664 <!-- claim: r2_err_ext_pct --> % near
  Gz = 10.777 <!-- claim: r2_err_ext_gz -->, past -2 % between
  Gz = 4.2044 <!-- claim: r2_gz_err2_lo --> and 31.903 <!-- claim: r2_gz_err2_hi -->, and back
  towards zero where the outlet nears equilibrium. There is no single Graetz number below
  which the module total fails.
- **The profile matters.** The wall Reynolds number is 7.4612 <!-- claim: r1_rewall_lo --> to
  13.43 <!-- claim: r1_rewall_hi --> (assumed Schmidt range), far from small suction. With a
  plug profile the averaged film misses the module total by
  -43.595 <!-- claim: v5_plug_rec_gap_pct --> %: the module-total verdict holds for the
  parabolic profile only.
- **Printed numbers.** Fig. 9's reference fluxes give a membrane permeance times thickness
  that spreads by 1.0638 <!-- claim: r3_spread_pct --> %; 1 <!-- claim: r3_n_exceed --> of 3
  pairs against the 100 micrometre panel differ by more than the half-units of the printed
  fluxes. Eq. (1) with Eq. (6) as printed gives CPC 0.80493 <!-- claim: r5_cpc_printed -->
  at the test condition against 0.70788 <!-- claim: r5_cpc_exact --> for the exact film (M5:
  which form the authors coded is not decidable).
- **Fig. 10 (figure reads).** The membrane permeance implied by Eq. (10) drops by
  20.463 <!-- claim: r6_drop_pct --> % from 10 to 26 bar, against a reading error of
  4.1761 <!-- claim: r6_read_err_pct --> %. With the printed Re and an assumed viscosity the
  printed depletion needs a tube 1.5071 <!-- claim: r6_L_mm_lo --> to
  2.5635 <!-- claim: r6_L_mm_hi --> mm long: a question the printed material raises, not
  decided here.
- **Misprints:** Eq. (1) lacks the minus sign (M1), Eq. (5) the subscript of x (M2), Eq. (6)
  prints "=" for "-" (M3); q and r of Eqs. (3) and (8) are not given (M4).

## Things to know before reusing this

- The 2-D tube is the page's extension; the paper prints no module, and not the inputs of its
  maps. The maps are turbulent (Table 4) and nothing here transfers to them.
- The Schmidt number and the viscosity are assumptions, not the paper's.
- Not covered: the six-component model, surface steps, sweep gas, counter-current permeate.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | the binary film relations (exact, Eq. 1 with Eq. 6 or 7, low-flux), CPC and its zero-flux limit, the 1-D film module, and `TubeMarch`, the 2-D tube (pymrm finite volumes, bordered newton per BDF2 step) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/caravella-2009-printed-numbers.csv` | Eq. (9)'s constants, Tables 4 and 5, the stated ranges and Fig. 9's printed fluxes, as strings, with a provenance sidecar |
| `data/caravella-2009-fig10-model-reads.csv` | figure reads of Fig. 10's model flux and CPC, with a sidecar stating the method and resolution |
| `data/check_values.csv` | values of the page's functions at the test condition from a second route, asserted before any claim |
