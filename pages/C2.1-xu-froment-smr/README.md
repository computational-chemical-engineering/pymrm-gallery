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
Table 6 round trip in the notebook is the check on that reading — and section 6a
measures how sharp it actually is.

## The headline that could not see the error it exists to catch

The round trip used to be reported as `rt.rel_pct.max()` = **1.4432 %**, "worst
deviation". That is a `max()` over seven rows, and its baseline scatter comes
from the *printed precision of the activation energies*: E is given to
0.1 kJ/mol, which at T_ref = 648 K is worth 0.93 % on A by itself, so `k1` and
`k3` sit at 1.28 and 1.44 % with nothing wrong. **Any single-digit slip landing
under that ceiling leaves the headline unchanged.** Measured, injecting the
defect that actually happens when a table is read off a page render:

| injected defect | worst % (old headline) | sum of per-row precision ratios |
|---|---|---|
| as published | 1.4432 | 10.438 |
| `K_CO` A 8.23e-5 → 8.32e-5 | **1.4432** | 13.670 |
| `k2` A 1.955e6 → 1.965e6 | **1.4432** | 14.504 |
| `k3` A 1.020e15 → 1.002e15 | **1.2844** (moves *down*) | 9.321 |
| `k1` E 240.1 → 240.2 kJ/mol | **1.4432** | 9.694 |
| `k1` A 4.225e15 → 4.226e15 (+0.024 %) | 1.4432 | 10.462 — **not detectable** |
| `K_CH4` value 0.1791 → 0.1719 | 3.8854 | 31.550 |
| `k1` A lost decade | 887 | 926 |

**A worst-case summary is only as sharp as its own baseline scatter.** The fix
is to divide each row by *its own* printed precision and to summarise with a
statistic every row reaches — the sum of the seven ratios, which moves for every
injection except one, and that one is a 0.024 % change in a four-significant-
figure number that no round trip built from these digits could resolve. Both
numbers are published; compare the sum, not the max.

This is a fourth shape of "the check that cannot fail", distinct from the three
in `docs/guard-structure-audit-2026-08-05.md`: the two sides are genuinely
independent and the check does have power — it just has less resolution than the
error class it is aimed at.

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
   paper's Figures 4 and 5, and by **2.246** for fresh catalyst. **Neither of
   those two factors is validated here** — the runs they belong to are chemistry
   this page does not model. Section 6d executes them and says so.
3. **No K_CO2 term.** It was never statistically significant and correctly has
   no term in the denominator. Adding one silently changes the model.

Also: the paper does not tabulate the equilibrium constants K1, K2, K3 at all.
They must come from elsewhere; the notebook checks the correlation it uses
against the paper's own Table 7 reaction enthalpies (worst 1.9 %). **That check
compares slopes only.** A +1 slip in either hard-coded intercept multiplies that
K by e = 2.718 at every temperature and leaves it at 1.93 %. What catches it is
the comparison against the 61 measurements, where the mean deviation goes
0.0017 → 0.0160 and 0.0148. Both are in section 6b/6c.

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
