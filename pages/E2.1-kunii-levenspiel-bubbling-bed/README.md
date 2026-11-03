# E2.1: Kunii–Levenspiel bubbling bed model

From Kunii, D. & Levenspiel, O., *Bubbling bed model for kinetic processes in
fluidized beds. Gas–solid mass and heat transfer and catalytic reactions*,
Ind. Eng. Chem. Process Des. Dev. **7**(4) 481–492 (1968).
[doi:10.1021/i260028a001](https://doi.org/10.1021/i260028a001)

- **Structures:** `S7` (multiphase coupling)
- **Runtime:** about 5 s
- **Tier:** 6. Nothing here is compared with a measurement.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | the closed form (equations 2 to 57), `BubblingBed` (the three-region distributed model on pymrm operators) and `extended_exact` (its eigenmode solution) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata |
| `data/kunii_levenspiel_1968_appendix_values.csv` | the numbers the authors printed in their own three worked appendices, plus its provenance sidecar |

The CSV is also loaded cross-page by `E1.2`.

## The check the paper pays for

The paper prints **three worked appendices**, A for gas–solid mass transfer, B
for gas–solid heat transfer, C for catalytic conversion, each ending in a stated
result. Between them they exercise every equation the model has, and none of them
requires digitising anything:

* appendix C's eight intermediates recomputed to 0.46 % <!-- claim: appendixC_max_dev --> worst case;
* appendix A's chain ending in `(Sh)_over-all = 0.045 Re − 0.0101` to 0.57 % <!-- claim: appendixA_max_dev -->;
* appendix B's `(Nu)_apparent = 0.070 + 0.051 Re` to 1.8 % <!-- claim: appendixB_max_dev -->.

## Identities, labelled as such

Three internal checks hold for any values of the four coefficients the paper
prints (`0.711`, `4.5`, `5.85`, `6.78`): the `K_r → 0` limit of equation 49, the
reactant balance of the distributed model, and the ceiling's independence of
`u_0`. They catch a wrong *formula*, never a wrong *number*, and the notebook's
break table shows both halves: they stay inside their bounds under every
coefficient injection and fail on the formula defect each is built for (a
`gamma_e` off equation 15's closure, an exchange block that does not conserve the
reactant, a spurious `u_0` dependence in `K_bc`). The limit is blind to `gamma_c`
and `gamma_e` being swapped, since only their sum survives it.

## The distributed model, and its exact solution

Equations 48b and 48c freeze the cloud and the emulsion. `BubblingBed` puts their
convection back: bubble gas up at `u_0 − u_mf`, cloud and wake gas up with the
bubble, emulsion gas taking the rest, which is downward for this bed, and the two
ends coupled by the recirculation loop of the paper's figure 6. The problem is
linear, so it is one sparse solve per `K_m`:

* one `construct_convflux_upwind` operator per region on shape `(n, 1)`, embedded
  in the `(n, 3)` state with `kron`, because the regions enter at opposite ends;
* inlets `{"a": 0, "b": 1, "d": 1}` with the inlet values as two extra unknowns
  through `shapes_d`, and the two mixing equations as two extra rows;
* outlets `{"outflow": True}`, so the exit is the face value the operator
  transports, read with `compute_boundary_values`.

The same three linear ODEs are solved exactly by `extended_exact` through their
eigenmodes, each normalised at the end where it is largest. The pymrm exit
converges to it at first order; at 1600 cells its worst relative error in the
unconverted fraction, 8.4e-4 <!-- claim: extension_pymrm_vs_exact_n1600 -->, is at the
largest `K_m`, 8. The extension raises the unconverted fraction by
7.2 % <!-- claim: extension_min_change --> to 14.6 % <!-- claim: extension_max_change -->
over figure 9's range of `K_m`, and equation 54 can match it with one fitted
bubble size only at one `K_m`: the `d_b` it needs runs from
5.00 <!-- claim: equivalent_db_max --> cm at `K_m = 0.5` to
3.85 <!-- claim: equivalent_db_min --> cm at `K_m = 8`.

## Two printed slips, recorded not repaired

**Appendix C** prints `u_0 = (6.6 + 9.9 + 13.2 + 20)/5 = 13.2 cm/sec`, four terms
divided by five. The result 13.2 is used in every later step and is confirmed by
the paper's own `u_b = 13.2 − 2.1 + 42.8 = 53.9`, so the value is not in doubt; a
term was dropped in typesetting. It is stored as printed.

**Appendix B** prints `(1 − ε_f)(u_0 + 5.8) = 8.70` where its own stated
`ε_mf = 0.50` and `u_br` give a smaller number. The printed final answer follows
from 8.70, so 8.70 is what the authors used; it corresponds to `ε_mf` =
0.447 <!-- claim: appendixB_emf_implied -->, and the stated value would move the
line by 10.6 % <!-- claim: appendixB_own_emf_shift -->.

## Read the equations off renders, not the text layer

This scan's text layer renders equations 49, 50 and 56 as loose fragments and
drops decimal points: appendix C's `13.2 𝒦_m` comes out as `1.32 𝒦_m`. Every
equation and every constant on this page was read from a 600 dpi `pdftoppm`
render.

## No experimental data, deliberately

The paper's measurements are in figures 3, 5, 7, 8 and 9, all scatter plots. On
figure 9, three velocity series against a common set of model curves, the markers
at low `𝒦_m` overlap into a solid mass in which individual centres cannot be
located, and the three glyphs are not reliably separable there either. The
extraction is staged for maintainer review under `queue_cases/E2.1/review/` and is
**not** on this page.

## What the model turns out to say

As `K_r → ∞` equation 49's bracket collapses to `K_bc/K_r`, so
`𝒦_f → K_bc L_f/u_b`, a ceiling with no rate constant in it. Substituting
equation 4 and `(1 − δ)u_b = u_br` gives

```
K_f_max = K_bc(d_b) (1 − ε_m) L_m / [(1 − ε_mf) u_br(d_b)]
```

which at fixed bubble size contains **no u_0 at all**. For the appendix C bed that
is 99.13 % <!-- claim: ceiling_conversion_db37 --> conversion and no more,
whatever catalyst is loaded.

## Corrections

- 2026-11-03: said the global mass balance of the distributed model closes to
  "8e-6" and that this residual "is first-order upwind truncation (observed order
  0.91, halving on every refinement)", that "check 6 measures grid error, not
  physics", and that the 0.711 transposition moved it by enlarging that
  truncation. The residual came from the zero-gradient outlets: the upwind flux
  carried an extrapolated outlet face value while the exit and the recirculated
  emulsion gas were read from the last cells. With `{"outflow": True}` outlets the
  balance closes to round-off on every grid and `K_m`, below 1e-10
  <!-- claim: mass_balance_dev -->; the zero-gradient construction is kept as a
  break row and leaves a gap of -8.1e-6 <!-- claim: mass_balance_zero_gradient_n2000 -->
  at 2000 cells and `K_m = 2`.
- 2026-11-03: said the extension's effect on the unconverted fraction is
  "comparable to moving d_b from 3.7 to 4.2 cm" and that "the fitted bubble size
  absorbs it". The bubble size at which equation 54 reproduces the extended exit is
  4.20 <!-- claim: equivalent_db_Km2 --> cm only at `K_m = 2`; it runs from
  5.00 <!-- claim: equivalent_db_max --> cm at `K_m = 0.5` to
  3.85 <!-- claim: equivalent_db_min --> cm at `K_m = 8`, so one fitted bubble size
  does not absorb it across figure 9. The 5.00 cm at `K_m = 0.5` reflects how weakly
  the unconverted fraction depends on `d_b` there.
- 2026-11-03: printed the extension at 1600 cells from the finite-volume solve:
  14.7 % <!-- claim: extension_pymrm_change_Km4_n1600 --> at `K_m = 4` and
  10.7 % <!-- claim: extension_pymrm_change_Km8_n1600 --> at `K_m = 8`, and 8.2 %
  <!-- claim: extension_paper_pymrm_change_Km8_n1600 --> at `K_m = 8` with the
  emulsion velocity of equations 5 and 7. The exact eigenmode route gives
  14.6 % <!-- claim: extension_exact_change_Km4 -->,
  10.6 % <!-- claim: extension_exact_change_Km8 --> and
  8.1 % <!-- claim: extension_paper_exact_change_Km8 -->; the earlier values carried the
  first-order error of the 1600-cell solve.
