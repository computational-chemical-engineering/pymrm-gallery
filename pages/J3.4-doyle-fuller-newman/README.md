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
| short-time anode concentration vs the semi-infinite solution | 7.7 mol/m³ | Fig. 4 to ±5 |
| grid spread over a 4× refinement | 0.14 mV | — |
| time-step spread (added 2026-08-05) | 0.27 mV, observed order 0.95 | — |
| total salt conservation — **structural, see below** | 3e−16 | — |

## The check that was ranked first and catches nothing it claimed

`salt_conservation_drift` was **first of six** in the Validation section, under
the claim that *"any sign error in the migration term, the transference number or
the reaction coupling breaks it"*. It breaks on none of the three.

It is a telescoping identity. With `r_c = eps*(c-c_o)/dt + div @ n_salt` and
`n_salt[0] = n_salt[-1] = 0`, `sum(dx * div(n)) = n[N] - n[0] = 0` exactly, so
`sum(eps*dx*c)` is constant for **any** interior flux expression and **any**
parameters. Measured on the page's own residual:

| injected defect | salt drift | short-time check | pore-wall check | *V*(*u* = 0.5) |
| --- | --- | --- | --- | --- |
| as published | 3.33e−16 | 7.7 mol/m³ | 1.1e−12 | 1.9850 V |
| migration sign flipped in `n_salt` | **4.44e−16** | **169.8** | 4.4e−12 | **2.1643 V** (+179 mV) |
| `(1−t⁰₊)` written as `t⁰₊` | **3.33e−16** | **72.7** | 4.8e−12 | **2.0911 V** (+106 mV) |
| reaction coupling sign wrong in the solid | **3.33e−16** | 6.7 | — | **cell does not discharge** |
| reaction coupling ×1.1 in the charge eq. | **2.22e−16** | 7.7 | **9.1e−02** | 1.9849 V |
| no-anion BC on `n_salt` dropped at the Li | **8.67e−01** | 88.2 | — | cell dies |

It is now ranked **last** and labelled `STRUCTURAL`. It is kept, because the last
row is a real error class — the salt-flux boundary conditions — and nothing else
on the page sees it. It is also 3.3e−16, below `check_agreement.py`'s
`ABS_FLOOR = 1e-12`, so CI never compared it either.

**The coverage the old ranking claimed does exist — it is just elsewhere.** The
short-time semi-infinite comparison moves 22× on the migration sign and 9.4× on
the transference number, so it is now ranked first; the Figure 2 comparison
catches the reaction-coupling sign error, because the cell stops discharging.
One gap is left stated rather than papered over: a coupling error in the *solid*
balance alone (`ai` scaled, not sign-flipped) moves nothing here by more than a
millivolt, and no check is constructible from the paper's published results.

Note the pore-wall check is partly structural too: `r_p = div@i2 - ai` with
`i2[0] = I` and `i2[-1] = 0` forces `sum(dx*ai) = -I` identically, so the *mean*
flux is guaranteed. What it does catch is the same `ai` failing to appear in both
balances.

## The grid study was refining the wrong knob alone

`Cell.march` fixes `dt, dt_max = 2e-5*t_ref, 6e-3*t_ref` growing ×1.3 — a schedule
that is a function of the **current alone** and completely decoupled from `n_s`
and `n_c`. So the grid-independence study held the temporal error *constant by
construction*, and the page had no bound on it anywhere.

`J3.5` had already implemented the correct fix — a `dt_scale` that multiplies
`dt0` **and** `dt_max`, since on a growing schedule the step reached at a given
time is set by the ceiling, not by the starting value. It was never inherited
back. It is ported here. Measured: observed order 0.95, and the temporal error at
the production schedule is **0.273 mV, 2.0× the 0.138 mV grid spread**. Both are
far below the 25.5 mV model-vs-figure gap, so no conclusion moves — but the
larger of the two numbers was the one nobody had.

**`J3.1` shares this `march` and still does not have `dt_scale`.** Its
`grid_V_spread_mV`, `grid_mean_eta_spread_frac` and `grid_peak_eta_ratio` are all
in `agreement.json` and all hold `dt` constant. Not fixed here.

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
