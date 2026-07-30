# E2.1 — Kunii–Levenspiel bubbling bed model

From Kunii, D. & Levenspiel, O., *Bubbling bed model for kinetic processes in
fluidized beds. Gas–solid mass and heat transfer and catalytic reactions*,
Ind. Eng. Chem. Process Des. Dev. **7**(4) 481–492 (1968).
[doi:10.1021/i260028a001](https://doi.org/10.1021/i260028a001)

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 4 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata |
| `data/kunii_levenspiel_1968_appendix_values.csv` | the numbers the authors printed in their own three worked appendices, plus its provenance sidecar |

## The check the paper pays for

The paper prints **three worked appendices** — A for gas–solid mass transfer, B
for gas–solid heat transfer, C for catalytic conversion — each ending in a stated
result. Between them they exercise every equation the model has, and none of them
requires digitising anything. That is the whole validation strategy of this page:

* appendix C's eight intermediates recomputed to **0.46 %** worst case;
* appendix A's chain ending in `(Sh)_over-all = 0.045 Re − 0.0101` to **0.57 %**;
* appendix B's `(Nu)_apparent = 0.070 + 0.051 Re` to **1.8 %**.

## Two printed slips, recorded not repaired

**Appendix C** prints `u_0 = (6.6 + 9.9 + 13.2 + 20)/5 = 13.2 cm/sec` — four
terms divided by five. The result 13.2 is used in every later step and is
confirmed independently by the paper's own `u_b = 13.2 − 2.1 + 42.8 = 53.9`, so
the value is not in doubt; a term was dropped in typesetting. It is stored as
printed.

**Appendix B** prints `(1 − ε_f)(u_0 + 5.8) = 8.70` where its own stated
`ε_mf = 0.50` and `u_br = 15.7 cm/s` give **7.87**. The printed final answer
follows from 8.70, so 8.70 is what the authors used — it corresponds to
`ε_mf = 0.447`. The page prints both lines and lets the reader see the 10.6 %
difference. Neither slip is fixed by inference.

## Read the equations off renders, not the text layer

This scan has the best text layer of any paper in the gallery — 12.4k characters
per page — and it is **still** unusable for equations. It renders equations 49,
50 and 56 as loose fragments, and it drops decimal points: appendix C's
`13.2 𝒦_m` comes out as `1.32 𝒦_m`. Every equation and every constant on this
page was read from a 600 dpi `pdftoppm` render.

## No experimental data, deliberately

The paper's measurements are in figures 3, 5, 7, 8 and 9, all scatter plots. On
figure 9 — the sharpest of them, three velocity series against a common set of
model curves — the markers below `𝒦_m ≈ 1.8` overlap into a solid mass in which
individual centres cannot be located, and the three glyphs (`⊙`, `●`, `○`) are
not reliably separable there either. The extraction is staged for maintainer
review under `queue_cases/E2.1/review/` and is **not** on this page.

So this is a **tier 6** page: everything is checked against the authors' own
arithmetic, and the page says so instead of implying more.

## What the model turns out to say

As `K_r → ∞` equation 49's bracket collapses to `K_bc/K_r`, so
`𝒦_f → K_bc L_f/u_b` — a ceiling with no rate constant in it. Substituting
equation 4 and `(1 − δ)u_b = u_br` gives

```
K_f_max = K_bc(d_b) (1 − ε_m) L_m / [(1 − ε_mf) u_br(d_b)]
```

which contains **no u_0 at all**. For the appendix C bed that is 99.13 %
conversion and no more, whatever catalyst is loaded.

## The pymrm trap on this page

Equations 48b and 48c are algebraic — the cloud and emulsion carry no gas. The
distributed model puts the convection back, and then the emulsion flows
**downward**, which makes the problem two-point: the emulsion enters at the top
and rejoins the feed at the bottom. Two things matter:

* an outlet left as `bc=None` makes the matrix singular (same trap as `F2.3`);
  every outlet here gets an explicit zero-gradient condition;
* the emulsion superficial velocity is a **stated closure** — the three
  velocities are required to sum to `u_0`. The paper's own equations 5 and 7 give
  a different value that overshoots `u_0` by 9.7 %, and the notebook runs both.
