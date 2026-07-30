# J3.4 — Doyle–Fuller–Newman lithium/polymer/insertion cell (P2D)

The 1993 pseudo-two-dimensional cell model that every physics-based lithium-ion
battery model descends from, rebuilt in pymrm from the paper — including the one
transport property the paper never prints.

- **Structures:** `S8` (coupled multi-phase transport with electroneutrality),
  `S10` (porous electrode with a distributed interfacial reaction)
- **Reference:** Doyle, Fuller & Newman (1993), *J. Electrochem. Soc.* **140**(6)
  1526–1533, doi:10.1149/1.2221597
- **Runtime:** ~1 min

## Provenance: tier 6, not experimental

**Figure 2 contains no measured points.** Every curve on it — charge,
open-circuit and four discharges — is the authors' own simulation; the text
cites measured discharge curves by reference only and plots none. This was
confirmed by the gallery maintainer on review of the original figure. The page
is a reference-solution reproduction in the sense of `D2.2`. Do **not** describe
it as validated against experiment.

## Agreement

| what | pymrm | paper |
|---|---|---|
| Eq. 16 vs the plotted open-circuit curve, 237 points | 3.3 mV mean abs. dev. | digitisation good to ~3 mV |
| four discharge curves, 884 points | 25.5 mV pooled, bias −23.8 mV | — |
| *u* at the 1.9 V crossing, all four currents | short by 0.023 in *u*, same sign each time | — |
| *S*_c (Eq. 26) at *I* = 10 | 1.02e−4 | 1.0e−4 |
| *S*_s (Eq. 27) at *I* = 10 | 0.153 | 0.15 |
| σ/κ | 2.1e5 | O(1e5) |
| *u* at the 1.7 V cutoff, *I* = 10 | 0.831 | 0.84 stated; 0.821 digitised |
| optimum cathode porosity | 0.60 | 0.60 |
| utilisation at that optimum | 0.928 | 0.97 |
| total salt conservation | 1e−15 | — |
| grid spread over a 4× refinement | 0.14 mV | — |

## The gap in the paper, and how it is closed

**Appendix A never gives the conductivity.** It says the conductivity of
PEO-LiCF₃SO₃ "was fit to a third-order polynomial" and cites its source; the
coefficients appear nowhere in the article and there is no figure of κ against
*c*. The salt diffusivity, the solubility limit and the transference-number fit
are all printed. The conductivity is not.

It is recovered by inverting the paper's own Eqs. 28 and 29, which contain
(1/κ + 1/σ) and whose values the paper prints (δ = 1.95, ν = 68). The two routes
give κ_eff = 7.97e−3 and 6.72e−3 S/m — 19 % apart — and both put σ/κ at O(1e5)
as the paper states. What cannot be recovered is the *concentration dependence*,
and the residual +64 mV bias on charge against −24 mV on discharge is where that
shows.

## Two ambiguities in the printed model

1. **The Bruggeman exponent in the cathode salt flux.** Table I writes the flux
   as ∇·(εD∇c) while the text defines D_eff = Dε^0.5. Read together they give
   Dε^1.5, which drives the electrolyte concentration to zero before half of the
   discharge — contradicting the paper's own Figure 3 — and the 1.9 V crossing
   comes out at a third to a half of the published utilisation at every current
   (mean |Δu| 0.294 against 0.023). Dropping the correction entirely overshoots
   *I* = 13 and 20 (mean |Δu| 0.086). The page uses Dε^0.5 and prints the table.
2. **A missing factor ν = 2 in Eq. 5.** As printed, the coefficient of ∇ln c is
   half what a 1:1 salt requires: at zero current the equation must reproduce the
   emf of a lithium concentration cell, which is (2RT/F)∫(1−t⁰₊)dln c. The page
   uses the corrected form. **Figure 2 does not settle this one** — the printed
   version gives a smaller mean deviation (18.2 vs 25.5 mV) but overshoots every
   transport-limited crossing by +0.026 in *u* where the corrected form
   undershoots by −0.023. The two bracket the published curves. The page says so
   instead of choosing by fit.

## What the page adds

The porosity optimisation the paper reports in one sentence, swept
continuously. Holding the theoretical capacity fixed, utilisation at the 1.7 V
cutoff peaks at a porosity of 0.60 — exactly where the paper puts it — and the
curve shows that the penalty for being below the optimum is far steeper than for
being above it. No new physics: the model is theirs, unchanged.

## Rebuilding the page

```bash
python build_page.py                    # regenerate index.ipynb from the builder
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF is needed to rebuild: all three datasets are committed CSVs.
