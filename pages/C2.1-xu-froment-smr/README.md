# C2.1 — Steam methane reforming: the Xu–Froment intrinsic kinetics

Three rate equations from 1989 that almost every published reformer simulation
still uses. This page rebuilds them from the printed tables and puts them back
against the measurements they were fitted to.

- **Structures:** `S1` (0D/1D reaction network), `S2` (plug flow with reaction)
- **Reference:** Xu & Froment (1989), AIChE J 35(1) 88–96, doi:10.1002/aic.690350109
- **Runtime:** ~5 s

## Agreement

0.0017 mean absolute deviation in conversion over 61 digitised measurements at
four temperatures (773–848 K), worst point 0.0058, which is 2.7 % of the mean
measured conversion. The digitisation error is about 0.0006. **Nothing was
fitted** — every kinetic parameter is Table 6 as printed.

The model is biased high, and the bias grows with space time: +0.0012 for
W/F < 0.20, +0.0032 for W/F > 0.28. That is the equilibrium end of the range,
and a 5 % shift in the equilibrium constants moves the predictions there by
0.0026–0.0039 — the same size. The equilibrium constants are the one ingredient
that does *not* come from this paper.

## Data

Two datasets, both with provenance sidecars.

`xu-froment-1989-parameters.csv` — Tables 5 and 6. **The PDF text layer mangles
every exponent** (`8.664 lo-'` for 8.664e-7), so these were read from a 600 dpi
render of journal page 94 and none was repaired by inference. The Table 5 ⇄
Table 6 round trip in the notebook reproduces all seven preexponential factors
to within 1.44 %, which is the check on that reading.

`xu-froment-1989-conversion.csv` — the experimental markers of Figures 2 and 3.
The plotted curves are the authors' own model and were deliberately not
extracted. All 30 Figure-2 points pair with a Figure-3 point to within 0.0027 in
W/F, from two independently calibrated digitisations.

## Three traps in this paper

1. **Split reference temperature.** T_ref = 648 K for k1, k2, k3, K_CO, K_H2 but
   **823 K** for K_CH4 and K_H2O. Forcing 648 K on all seven puts K_CH4 out by a
   factor of 0.22 and K_H2O by a factor of 33.
2. **Three activity levels.** Tables 5 and 6 are the *steam-reforming reference*
   level (partially deactivated catalyst) — what the Figure 2/3 data are, so no
   correction. Multiply the rate coefficients by **1.225** to reproduce the
   paper's Figures 4 and 5, and by **2.246** for fresh catalyst.
3. **No K_CO2 term.** It was never statistically significant and correctly has
   no term in the denominator. Adding one silently changes the model.

Also: the paper does not tabulate the equilibrium constants K1, K2, K3 at all.
They must come from elsewhere; the notebook checks the correlation it uses
against the paper's own Table 7 reaction enthalpies (worst 1.9 %).

## Rebuilding the page

```bash
python build_page.py                    # regenerate index.ipynb from the builder
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

To regenerate the digitised dataset you need your own lawfully obtained copy of
the PDF (it is not in the repository, and must not be):

```bash
python extract_figures.py ~/papers/pymrm-gallery/"AIChE Journal - January 1989 - Xu.pdf"
```

That script is deterministic — the candidate indices that survived the visual
audit are constants in it — and re-runs the monotonicity, cross-figure pairing
and carbon-closure checks, exiting non-zero if any of them fails.
