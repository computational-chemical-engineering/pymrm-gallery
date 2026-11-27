# I1.9: Hufton, Mayorga and Sircar (1999), the equilibrium limit of sorption-enhanced reforming

J. R. Hufton, S. Mayorga and S. Sircar, "Sorption-enhanced reaction process for hydrogen
production", *AIChE J.* **45**, 248-256 (1999), doi:10.1002/aic.690450205. The paper's case
for sorption-enhanced reforming rests on one batch equilibrium of steam reforming and the
water-gas shift with a fraction f of the CO2 removed (Eqs. 3 to 10, p. 250), printed as
numbers on p. 251 and in Table 1 (p. 255). This page re-solves that equilibrium and tests
the printed numbers and the sentences built on them.

- **Structure:** `S10` (algebraic equilibrium system), tier T1, data tier 6 (nothing
  measured). The catalogue line gives S4 + S12; the paper prints no column model, only the
  batch equilibrium.
- **Equilibrium constants:** not printed in the paper; taken from page `C2.1`'s
  correlations, with every verdict re-evaluated over an assumed band in both constants.

## Findings

- **The printed equilibrium states are reproduced within the K uncertainty.** The f = 0
  state at 450 C and 50 psig: 5 <!-- claim: r1_n_consistent --> of 5 printed cells are
  consistent within the K band; Table 1's 645 C row: 4 <!-- claim: r5_n_consistent --> of 4.
- **Table 1's "CH4 Conversion" is the carbon-based conversion**, not the chi/4 the text
  defines: the carbon-based conversion reaches the printed value at 641.12 <!-- claim: r3_Tstar_X --> C
  and chi/4 at 667.26 <!-- claim: r3_Tstar_chi4 --> C (55 psig, f = 0), so "~650°C" matches
  only the former.
- **The 450 C row labelled 55 psig carries the 50 psig numbers:** 4 <!-- claim: r4_n_p50 -->
  of 4 cells match p. 251's 50 psig values, 1 <!-- claim: r4_n_p55 --> of 4 match them moved
  to 55 psig.
- **"> 95%" H2 and "< 100 ppm" carbon oxides "when 99.9% or more of the CO2 produced" is removed do not follow.** At 1 - f of
  one in a thousand the dry carbon oxides are 257.41 <!-- claim: r2_999_COx_ppm --> ppm and
  the dry H2 is 94.129 <!-- claim: r2_999_H2 --> %; the thresholds are 1 - f =
  0.00037425 <!-- claim: r2_om_star_COx_ppm --> and 0.00073198 <!-- claim: r2_om_star_H2 -->.
  Both verdicts hold across the K band; the "> 80%" conversion is undecided there.

## Files

- `index.py` / `index.ipynb`: the page (paired by jupytext).
- `model.py`: the equilibrium (`Case`, `solve_extents`, `State`), the balance route
  (`solve_species`), `unique_root`, `Branch`; `test_model.py` tests it offline.
- `data/hufton-1999-printed-numbers.csv`: the printed numbers, as strings.
- `data/check_values.csv`: independent values from a second route.
