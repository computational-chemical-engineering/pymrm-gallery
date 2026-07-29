# B3.1 — Shrinking-core model for non-catalytic gas–solid reaction

From Yagi, S. & Kunii, D., *Studies on combustion of carbon particles in flames
and fluidized beds*, 5th Symposium (International) on Combustion, 231–244 (1955).

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 4 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata, validated against `models.yaml` |

No `data/` directory: this is a tier 6 analytic page and has no experimental
content. That is stated on the page rather than papered over.

## Why this paper and not a textbook

Every textbook gives three limiting conversion–time laws for the shrinking core.
This paper gives the **single equation all three come from** — Eq. 6 — and it is
rarely quoted. Reproducing it is the point of the page.

The scan's text layer is unusable for equations (it renders `θ_B` as `0B` and
drops the exponents in Eq. 6's numerator), so both equations were read off 600 dpi
renders.

## The check that establishes the reading

Eq. 6 was read from a page image, so it needs verifying. Two independent routes:

1. It is exactly 0 at *r*/*R* = 1 and exactly 1 at *r*/*R* = 0 for **any**
   (ω, γ) — which only happens if the numerator collapses to the denominator term
   by term, i.e. only for the printed coefficients.
2. Integrating the moving boundary from Eq. 5's three resistances, without ever
   looking at Eq. 6, gives the same function to **6.9 × 10⁻¹⁶** over six decades
   of resistance ratio. That route independently recovers the factor 3 on the film
   term and the 12 in *k*_d1 = 12𝔻/*D*_p.

It also collapses to all three classical laws in the right limits, to 2.4 × 10⁻⁸.

## What the page adds

A regime map. For every (ω, γ), how far is the closest single-resistance law from
the full equation? Over six decades in each group, some single law is within 5 %
on 71 % of the plane and 2 % on 52 %; the worst case anywhere is 20 %, at ω ≈ 11,
γ ≈ 0.8 — where no resistance dominates and all three must be kept.

## The sign trap

Both shell boundary conditions use the **outward** normal. At the core face that
normal points toward the centre, so the reacting condition is
`{"a": D, "b": kc, "d": 0.0}` — structurally the same as the film condition, not
opposite to it. Write it with the inward normal and the core appears to emit gas.
Same trap as `H1.7`.

## Not computed, deliberately

`k_c1` cannot be recovered from the paper: Parker and Hottel's correlation is
given for the specific combustion rate `K_c`, and the unit conversion to `k_c1`
is not printed. So no absolute burnout times appear, and the paper's Figure 6 is
not reproduced. Everything is in ω, γ and θ/θ_B, which is what Eq. 6 needs.
