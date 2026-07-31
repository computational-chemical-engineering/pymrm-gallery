# J3.5 — Single particle models with electrolyte (SPM / SPMe)

Marquis et al.'s (2019) asymptotic reduction of the Doyle–Fuller–Newman model —
the SPM at leading order, the SPMe at first order — applied to the gallery's
own 1993 half-cell (`J3.4`), with the error scaling the asymptotics promise
measured against the full model.

- **Structures:** `S8` + `S10`, reduced towards `S3`
- **Reference:** Marquis, Sulzer, Timms, Please & Chapman (2019),
  *J. Electrochem. Soc.* **166**(15) A3693–A3706, doi:10.1149/2.0341915jes
  (open access)
- **Foundation:** the `J3.4` page — its DFN implementation and its
  maintainer-reviewed Doyle 1993 parameter set are reused unchanged
- **Runtime:** ~2.5 min (most of it the converged DFN reference sweep)

## Provenance: tier 6, twice over

Every comparison on this page is reduced-model against full-model, and the
full model (`J3.4`) is itself a reproduction of a published *simulation*. No
measurement exists anywhere in this chain. Do **not** describe anything here
as validated against experiment.

## Agreement

| what | pymrm | paper / theory |
|---|---|---|
| Tables II + III recomputed from Table I (21 values) | ≤ 1.2e−3 rel | printed digits (truncated, not rounded) |
| in-text C_e = 5.1e−3 C (p. A3701) | 4.195e−3 C recomputed | Table III's 4.19e−3 C confirmed; sentence wrong |
| SPM voltage error vs DFN, slope in I | 1.00 (structural here) | O(C_e): slope 1 |
| SPMe voltage error vs DFN, local exponent | 2.07 → 2.70 across the sweep; window fits 2.11 / 2.26 / 2.53 | O(C_e²): slope 2 as C_e → 0 |
| SPM / SPMe error at I = 1 A/m² | 12.14 mV / 0.0452 mV (×268) | "order of magnitude", their Fig. 4 |
| steady SPMe operator vs closed-form quadratic | 3.5e−5, O(h²) ratio 4.02 | exact solution |
| SPMe electrolyte profile vs DFN, mid-discharge | 0.06 % of c_0 (I=1), 1.4 % (I=5) | — |
| salt conservation (SPMe) / coulomb identity (DFN) | 4e−13 / 7e−16 | structural / their Eq. 32 |
| C_e for this cell at I = 10 | 0.153 | ≡ Doyle's S_s, printed 0.15 in 1993 |

## What the slope test does *not* establish

The page ships a mutation table for this, because a slope test is only as
strong as the terms it can see:

- a 10 % error in the ohmic (δ_c/3) or concentration (factor 2) coefficient
  multiplies the I = 1 error by 9.8× / 14.9× and collapses the exponent to
  ~1 — those two are resolved sharply;
- deleting the electrode-averaged exchange-current correction entirely changes
  the error by **0.01 %** — that term is *not* tested by anything on this page;
- evaluating the anode overpotential at c_0 changes it by 10 % (exponent
  1.83) — visible, not sharp;
- moving the reconstructed κ by its published 19 % uncertainty changes the
  measured error by 9 %: shared parameters cancel, so **no comparison here
  tests the parameter set**.

The bottom of the sweep (I = 0.25 A/m²) is set by the DFN integrator, which
does not run at I = 0.125 at any time step tried — not by the SPMe.

## What the page shows that neither paper does

- All five a-priori applicability conditions (Marquis Table VI) evaluated for
  the 1993 cell: valid at a few A/m², marginal at I = 10 (C_e = 0.15,
  κ̂_e = 1.04), failed at I = 20 (0.31, 0.52) — where the first-order profile
  itself predicts negative electrolyte concentration above I = 10.2 A/m²,
  which is exactly the transport-limited collapse Doyle et al. describe at
  I = 20 ("about 30 % of the cathode material", their text, p. 1529).
- Two of those five conditions are groups Doyle, Fuller & Newman printed
  themselves in 1993: C_e ≡ S_s (Eq. 27) and the solid-diffusion condition
  ≡ S_c (Eq. 26).
- A closed-form (piecewise quadratic) SPMe(S) for the half-cell.

## What is deliberately not attempted

Their Table IV (RMS errors on their graphite/LCO cell) needs OCP and
electrolyte property *functions* that live in their cited codes (DUALFOIL,
fastDFN), not in the printed paper. Those values are carried as context only.

## Files

- `build_page.py` — generates `index.ipynb` (run from this directory)
- `data/marquis-2019-table1-parameters.csv` — their Table I (+ sidecar)
- `data/marquis-2019-stated-results.csv` — their Tables II–IV and one in-text
  value (+ sidecar)
- physical inputs load cross-page from
  `pages/J3.4-doyle-fuller-newman/data/doyle-fuller-newman-1993-parameters.csv`
