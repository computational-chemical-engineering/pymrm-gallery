# F1.5: Vandu and Krishna, kLa per unit gas holdup in the churn-turbulent regime

C.O. Vandu and R. Krishna, "Volumetric mass transfer coefficients in slurry bubble
columns operating in the churn-turbulent flow regime", *Chem. Eng. Process.* **43**
987-995 (2004), [doi:10.1016/j.cep.2003.09.007](https://doi.org/10.1016/j.cep.2003.09.007).
The paper's estimation rule is kLa = eps_G P, with P = kLa/eps_G independent of the gas
velocity in the churn-turbulent regime and proportional to a power of the Schmidt
number. The page tests it against the paper's own exact vector plateaus (Figures 5, 7,
9, 11), its digitised points (Figures 4c, 8c, 10f) and its printed probe constants. A
correlation-and-data page; it calls no pymrm operator.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/vandu-krishna-2004-printed-numbers.csv` | Tables 1 to 3, the Fig. 2 and 3 caption values, as printed |
| `data/vandu-krishna-2004-vector.csv` | plotted values of Figs. 2b, 5, 7, 9 and 11, exact from the PDF vector paths |
| `data/vandu-krishna-2004-fig2a-curves.csv` | the five drawn sensor-response curves of Fig. 2a |
| `data/vandu-krishna-2004-fig4c.csv`, `-fig8c.csv`, `-fig10f.csv` | kLa/eps_G against U, read from the 600 ppi bitmaps |
| `data/digitise.py`, `data/digitise_raster.py` | the extraction; take a local PDF path (the PDF is not stored) |
| `data/check_values.csv` | the page's functions at one test condition, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it. Akita and Yoshida's constants
and Table II ranges are loaded from [F1.1](../F1.1-akita-yoshida-bubble-column/).

## Things to know before reusing this

- The rule's exponent is not pinned by the paper's data: the four exact plateaus give
  -0.409 <!-- claim: exp_free_slope --> ± 0.061 <!-- claim: exp_free_se --> on two
  degrees of freedom, the free slope differs from the paper's -1/3 by t = -1.24 <!-- claim: exp_t_vs_minus_third -->
  standard errors, and one liquid (Tellus oil, leverage 0.989 <!-- claim: exp_leverage_tellus -->)
  sets the slope. Schmidt number and viscosity correlate at 0.9925 <!-- claim: corr_lnSc_lnmu -->
  across the four liquids, so the data cannot tell them apart.
- No constant is printed. Read from the drawn line of Figure 5 it is
  3.2045 <!-- claim: fig5_line_C --> 1/s, and it misses the plateaus by
  -27.4 % <!-- claim: rule_miss_tetradecane --> (tetradecane) to
  +29.6 % <!-- claim: rule_miss_tellus --> (Tellus oil): use the plateau of the liquid
  when it is one of the four.
- "Practically independent of U": of the 6 <!-- claim: r3_n_series --> digitised
  series, 4 <!-- claim: r3_n_t_beyond_2 --> have a slope beyond two standard errors,
  and water is marginal (5 <!-- claim: r3_n_t_beyond_2_nolow -->
  without each series' lowest marker).
- The probe correction is not small: an uncorrected evaluation would return
  0.562 <!-- claim: probe_rate_min --> to 0.930 <!-- claim: probe_rate_max --> of the
  corrected kLa_L in the captioned cases.
- C carries 1/s: the rule holds for these columns, spargers and ambient air only.
