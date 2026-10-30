# H1.4 — Itoh Pd membrane reactor: shifting a dehydrogenation equilibrium

Reproduces Itoh (1987), *A membrane reactor using palladium*, AIChE Journal
**33**(9) 1576–1578, [doi:10.1002/aic.690330921](https://doi.org/10.1002/aic.690330921).

- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`ItohCase` for the printed constants, the equilibrium closure
  `kp_of_xeq` and its inverse, the pointwise `source` of Eqs. 1 and 7, the
  kinetics-free `ceiling`, the pymrm `MembraneReactor` with its tube and shell
  channels in a `BlockLayout`, co-current or counter-current, and the
  independent `conversion_ivp`), tested by `test_model.py`.
- **Runtime:** about 15 seconds.
- `data/itoh-1987-stated-values.csv` and its `.meta.yaml`: the stated values
  (tier 2: numbers printed in the paper) with their provenance.

Cyclohexane dehydrogenates over Pt/Al₂O₃ inside a 200 µm palladium tube;
hydrogen, and only hydrogen, leaves through the wall into an argon purge, so a
reaction whose equilibrium conversion is 18.7 % <!-- claim: equilibrium_conversion_stated -->
runs to a measured 99.7 % <!-- claim: conversion_measured -->. The page solves the
paper's two-channel plug-flow model (Eqs. 1 to 8) as a pymrm finite-volume
Newton system, lands on the measured conversion with nothing fitted, and then
measures how much that agreement can actually resolve, which is the part worth
reading.

## Validation route

Stated numerical results plus internal identities, no figure digitisation:

1. the stated measured conversion at the stated flows, the one experimental
   datum, modelled as 0.9983 <!-- claim: conversion_model --> (deviation
   0.13 % <!-- claim: conversion_dev_frac -->). **Read the resolving power with
   it:** at this operating point the model sits exactly on its algebraic
   fast-permeation and equilibrium ceiling, so varying α_H over a factor 33 or
   `k·V_r` over a factor 100 moves the number by below 1e-8 <!-- claim: sens_alphaH_x0.3_to_x10_maxabs_dX -->,
   and K_p may be divided by 3 or multiplied by 10 with the deviation between
   -0.19 % <!-- claim: sens_Kp_div3_dev_frac --> and +0.28 % <!-- claim: sens_Kp_x10_dev_frac -->.
   On the quantity that does vary, `1 − X`, the model deviates by
   -45 % <!-- claim: unconverted_dev_frac -->. The comparison establishes that
   the reactor reached its co-current fast-permeation asymptote and that the
   model puts that asymptote in the right place; it is not evidence about the
   kinetics or the permeance;
2. the printed permeation constant α_H recomputed from Eq. 1's own inputs
   (-0.03 % <!-- claim: alphaH_recomputed_dev_frac -->), a genuine internal
   identity that propagates to nothing;
3. K_p(473 K) reconstructed from the paper's stated equilibrium conversion,
   cross-checked against an independent van 't Hoff estimate (ratio
   1.58 <!-- claim: Kp_473_vanthoff_ratio -->);
4. the finite-volume solve against an independent stiff-IVP integration, outlet
   and whole profiles on a grid ladder (first order, mid-map outlet error
   4.4e-5 <!-- claim: fv_vs_ivp_absdiff_midmap_nz800 --> at 800 cells), hydrogen
   conservation below 1e-12 <!-- claim: h_balance_rel_error_countercurrent_nz3200 -->
   in both flow arrangements with outlet values read from the faces, and a
   kinetics-free fast-permeation ceiling that bounds the whole operating map.

Every number on the page is a registered claim, and a defect-injection table
re-runs each computation with one defect in place and checks that the claim's
own test catches it.

## Two misprints in the paper, and one printed value confirmed

The scan is a 300 dpi bilevel CCITT image, so a 600 dpi render is a 2×
upsample: every marginal glyph was read on the native bitmap.

- **Eq. 5** prints `v_H = u_H - 3(u_C0 - u_C)`. Substituting the paper's own
  hydrogen balance `u_H + v_H = 3(u_C0 - u_C)` into it gives `v_H = -v_H`, i.e.
  the printed equation forces `v_H ≡ 0` and contradicts Eq. 3's permeation term.
  Sign flip. No convention escape: Eq. 8, Figure 3 and Figure 1 all establish
  co-current flow, and the sentence below Eq. 6 makes `v_i` the separation-side
  flow rate, positive in `+L`.
- **K_P = 4.89e35 exp(+3190/T) Pa³** contradicts the paper's own equilibrium
  conversion: it leaves an unconverted fraction of 1e-25 <!-- claim: Kp_printed_unconverted_eq -->.
  A units re-reading leaves it 27 <!-- claim: Kp_printed_vs_needed_decades -->
  orders out, and the exponent's sign is independently wrong for an endothermic
  reaction. Decision table in the notebook. `K_p(473) = 2.357e11 <!-- claim: Kp_473_reconstructed_Pa3 --> Pa³`
  is reconstructed from the stated equilibrium conversion, every input printed.
- **`v_A0 = 11.8e-5` mol/s is NOT a misprint.** The superscript reads `-5`
  cleanly on the native bitmap. The page confirms it rather than correcting it,
  and the confirmation is kinetics-free: under a `10⁻⁶` reading the
  thermodynamic ceiling, infinite kinetics *and* infinite permeance, is
  X = 0.716 <!-- claim: ceiling_at_1e-6_reading -->, so the measured value would be
  impossible. Figure 4's hand-lettered abscissa exponent *is* illegible on the
  scan; nothing here depends on reading it.

Native-resolution verification crops are in `queue_cases/H1.4/review/`
(git-ignored: they are page images).

## Honesty notes

The single measured point sits exactly on the model's thermodynamic ceiling.
Check 1 of the notebook quantifies what that costs, and the page's summary lines
are written to match it. Figure 4's experimental markers at other purge rates
were deliberately not digitised (ranked validation policy). The counter-current
comparison in "What pymrm adds" is a pure prediction, computed on a fine grid
(n_z = 3200, ladder shown to 6400) with its hydrogen closure reported; its
observed order is still rising toward 1 on that ladder.

A non-blocking follow-up is recorded in `queue_cases/H1.4.yaml`: a trial
digitisation of Figure 4's *calculated* curves, which would add a check that
does discriminate K_p, is staged unreviewed in `queue_cases/H1.4/review/` and
is deliberately not used anywhere on this page.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/H1.4-itoh-membrane-dehydrogenation/index.py   # after editing index.py
python scripts/run_kernels.py pages/H1.4-itoh-membrane-dehydrogenation
GALLERY_OFFLINE=1 python -m pytest -q pages/H1.4-itoh-membrane-dehydrogenation/test_model.py
```

## Corrections

- 2026-10-30: stated that the counter-current hydrogen closure is a genuine
  O(h) discretisation error and the one conservation
  check that is not automatic. With the outlet values read from the faces the
  operators transported, the closure holds below 1e-12 <!-- claim: h_balance_rel_error_countercurrent_nz3200 -->
  in both arrangements; the O(h) came from reading the last cells next to a
  zero-gradient outlet (pitfall P8), which a break row reproduces.
- 2026-10-30: stated that the counter-current solve stalls above 8e-6 mol/s of
  purge, most likely at the shell inlet where `sqrt(v_H)` has an infinite
  derivative. On 400 cells it converges up to about
  9.4e-6 <!-- claim: cc_continuation_last_vA0 --> mol/s (at 9.0e-6 mol/s the
  outlet unconverted fraction is 9.8e-3 <!-- claim: cc_unconverted_9e-6_nz400 -->,
  n_z = 400). The stop is a solver limit, not a reactor limit: past it the
  Newton iterates convert the cyclohexane completely inside the reactor
  (unconverted below 1e-6 <!-- claim: cc_past_edge_iterate_unconverted --> at the
  outlet) with a residual of order 100, and with the square root replaced by a
  smooth `y/sqrt(y + 1e-8)` continuation stops at the same purge
  (9.4e-6 <!-- claim: cc_continuation_last_vA0_smoothed_sqrt --> mol/s), so the
  square root is not the cause.
- 2026-10-30: stated that the counter-current solution converges first order
  in h on a converged grid; the observed order on the ladder to 6400 cells is
  still rising toward 1 (0.87 <!-- claim: cc_observed_order --> from the three
  finest grids), so the grid is fine rather than converged, which the
  two-decimal quotes absorb.
