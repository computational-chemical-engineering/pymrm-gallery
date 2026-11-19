# C2.13 Three-way catalyst: Voltz CO and propylene oxidation kinetics

Voltz, Morgan, Liederman and Jacob (1973), *Ind. Eng. Chem. Prod. Res. Dev.*
12, 294-301, doi:10.1021/i360048a006.

The Langmuir-Hinshelwood rate equations for CO and propylene oxidation on
platinum, with inhibition by CO, propylene and NO, that most catalytic-converter
models use. The page

- transcribes Table I and evaluates the rate law, locating the rate maximum in
  CO by root finding (0.330 <!-- claim: co_max_550F --> mol % CO at the temperature and propylene of curve G);
- rebuilds the paper's isothermal integral reactor (eq 1) with pymrm
  (`IntegralReactor` in `model.py`) and checks it against a `solve_ivp` march;
- tests the reconstruction against the authors' own calculated lines of
  Figures 19 and 20, which shows that the rate equations take propylene in ppm
  as C6 (half the propylene ppm), as the paper's text defines: the spread of the
  inferred space-time factor is 0.047 <!-- claim: a_spread_dex --> dex on that
  basis and 0.143 <!-- claim: a_spread_other_basis_dex --> dex with plain ppm;
- compares with the measured markers: root-mean-square misses of 6.3 <!-- claim: rms_fig19_pp -->
  percentage points for CO and 14.7 <!-- claim: rms_fig20_pp --> for propylene.
  Table I was fitted to these runs, so this is reproduction, not validation;
- checks that the adsorption constants and activation energies used by `I1.2`
  (Oh and Cavendish 1982) are Voltz's, unit-converted at face value; on the C6
  basis found here, K2 and K3 for propylene mole fractions are half and a
  quarter of theirs.

The space velocity of the figures relates to V/F by a factor the paper does not
state; the page infers it from the calculated lines, and it is 0.59 <!-- claim: a_fit_over_a_0C -->
times the value for gas flow at 550 F, space velocity in litres at the ice point
and the stated bulk density.

## Files

- `index.py` / `index.ipynb`: the page (jupytext pair).
- `model.py`: rate constants, rates, the CO rate maximum, `IntegralReactor`
  (pymrm) and `march` (`solve_ivp`).
- `test_model.py`: fast checks, including one stored claim.
- `extract_figures.py`: regenerates the digitised CSVs from a local copy of the
  PDF (not in the repository).
- `data/`: Table I, the legend compositions, the measured markers and the
  calculated-line crossings, each with a `.meta.yaml` provenance sidecar.
