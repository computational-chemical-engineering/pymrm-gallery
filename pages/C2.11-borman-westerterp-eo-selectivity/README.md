# C2.11 Ethene epoxidation selectivity over silver (Borman and Westerterp, 1995)

Borman, P. C. and Westerterp, K. R. (1995), *Ind. Eng. Chem. Res.* 34, 49-58,
doi:10.1021/ie00040a002.

The paper fits three Langmuir-Hinshelwood rate laws for ethene epoxidation and combustion
over a silver catalyst, prefers model III, and draws selectivity conclusions from it. This
page tests whether those conclusions follow from model III with Table 6 as printed:

- Table 6 is consistent with its own rounding (every adsorption constant at the reference temperature and every
  lumped constant inside the interval its rounded inputs allow), so the defect below cannot
  be seen inside the table;
- as printed, model III gives 99.62 <!-- claim: sigma_T_F1_pct --> % selectivity at the
  Figure 6 condition, where the authors' own model III curves lie far lower, and
  99.60 <!-- claim: sigma_T4_F1_pct --> % at the standard condition of Table 4, where the
  measured selectivities are 80.0 <!-- claim: sigma_meas_row1_pct -->,
  70.95 <!-- claim: sigma_meas_row2_pct --> and 68.52 <!-- claim: sigma_meas_row3_pct --> %;
- one factor of 112.8 <!-- claim: F_fig6 --> on the combustion rate constant reproduces all
  three Figure 6 curves (root mean square 0.031 <!-- claim: fig6_rms_Ffit --> selectivity
  points); a round factor of a hundred leaves an offset of 2.46 <!-- claim: fig6_mean_F100 --> points;
- the trend signs (CO2 raises, ethene oxide lowers, water raises the selectivity) and the
  fall of selectivity with temperature hold for any factor; at the Figure 6 composition the
  apparent activation energies are 78.9 <!-- claim: Eapp1_T_kJmol --> and
  103.3 <!-- claim: Eapp2_T_kJmol --> kJ/mol, a larger gap than the lumped values the paper
  prints;
- eq. 8 gives r3/r1 = 0.727 <!-- claim: r3_over_r1_T --> at the Figure 6 condition, outside
  the range the paper prints.

The page says what factor the paper's own curves imply, never why. `D3.5` (planned) should
load `F_fig6` from this page's `agreement.json` with `load_claims`.

## Files

- `index.py` / `index.ipynb`: the page (jupytext pair).
- `model.py`: model III rates, selectivity, analytic derivatives, Table 8 and A-4/A-5, eq. 8
  and `BertyReactor` (the forward solve with pymrm).
- `test_model.py`: fast checks, including one check value and one stored claim.
- `extract_figures.py`: regenerates the digitised Figure 6 from a local copy of the PDF (not
  in the repository).
- `data/`: Tables 3, 4, 6 and 7, the printed numbers, the digitised Figure 6 and the check
  values of a second route, each with a `.meta.yaml` provenance sidecar.
