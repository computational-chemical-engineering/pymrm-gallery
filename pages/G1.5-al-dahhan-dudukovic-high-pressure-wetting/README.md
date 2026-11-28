# G1.5: Al-Dahhan and Dudukovic (1995), catalyst wetting efficiency at high pressure

M. H. Al-Dahhan and M. P. Dudukovic, "Catalyst wetting efficiency in trickle-bed
reactors at high pressure", *Chem. Eng. Sci.* **50**(15) 2377-2389 (1995),
[doi:10.1016/0009-2509(95)00092-J](https://doi.org/10.1016/0009-2509(95)00092-J). The
page evaluates eq. 9 with its printed constants at the paper's own measured pressure
drops (Figs 4b, 5b) and sets it against the wetting efficiencies it was fitted to (Figs
4a, 5a), the paper's parity plot (Fig. 6) and trend figure (Fig. 7); it tests El-Hisnawi's
eq. 1 at the lowest pressure, and solves with pymrm the pellet moment problem behind the tracer
method of Tables 1 and 2.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `model.py` | eq. 9, eq. 1, the metrics and refits, the input checks, the tracer relations and the pymrm pellet moment BVP |
| `test_model.py` | fast checks of `model.py` |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/al-dahhan-1995-printed-numbers.csv` | every constant, range, legend value and stated error the page uses, as printed |
| `data/al-dahhan-1995-fig4.csv`, `-fig5.csv` | wetting efficiency and pressure drop against L, spheres and extrudates, digitised |
| `data/al-dahhan-1995-fig6.csv`, `-fig7.csv` | the parity plot and the trend figure, digitised |
| `data/digitise.py` | the readings of Figs 4 to 7 and their calibration; writes the four figure CSVs |
| `data/nist-n-hexane-298k.csv` | n-hexane density and viscosity at 298 K (NIST WebBook), the one external input |
| `data/check_values.csv` | the page's functions at one test condition, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it.

## Things to know before reusing this

- Eq. 9 reduces to $K L^{1/3}\varepsilon_B^{-1/3}\psi_L^{1/9}$: the particle diameter
  cancels, so it has no particle-size dependence, and the bed enters only through
  $\varepsilon_B$. Do not use it to compare particle sizes.
- Route A (eq. 9 recomputed) gives a mean relative error of 3.32 % <!-- claim: r2_mre_capped_exp -->
  with predictions capped at 1 and 4.01 % <!-- claim: r2_mre_uncapped_exp --> uncapped; the paper
  prints "The mean relative error is 3.7%.", which its own Fig. 6 reproduces (3.67 % <!-- claim: r3_mre_B -->).
  Whether the authors capped, and which liquid properties they used, is not printed.
- The residuals depend on the bed: spheres -4.46 % <!-- claim: r4_bias_sphere -->, extrudates
  +2.11 % <!-- claim: r4_bias_extrudate -->, a split of 5.1 <!-- claim: r4_bias_z_cases --> standard
  errors with the nine case means as units. Refit the bed term, and the exponents if possible, for a
  new bed: within the paper's beds the pressure-drop exponent is 0.0790 <!-- claim: r4_m --> (case
  bootstrap 95 % 0.0704 <!-- claim: r4_m_lo --> to 0.0889 <!-- claim: r4_m_hi -->), below the printed 1/9.
- El-Hisnawi's eq. 1 lies above the low-pressure data (cases 2 and 3) at 16 <!-- claim: r5_above --> of
  17 <!-- claim: r5_n --> points, and "within the error margin of the correlation (± 10%)" at only
  9 <!-- claim: r5_within10 -->.
- Table 2 and the liquid-filled lines of Table 1 of the paper write the tracer variance without the factor 2 of the second cumulant:
  an effective diffusivity recovered with them is 0.5 <!-- claim: r6_de_recovered_over_true --> of the
  true value, which cancels in eq. 8 only when both variances use the same form.
