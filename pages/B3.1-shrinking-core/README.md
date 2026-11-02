# B3.1 Shrinking-core model for non-catalytic gas–solid reaction

From Yagi, S. & Kunii, D., *Studies on combustion of carbon particles in flames
and fluidized beds*, 5th Symposium (International) on Combustion, 231–244 (1955).

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb` |
| `index.ipynb` | the executed notebook; runs in about 27 s |
| `model.py` | eq. (6), the integrated resistances, the regime map and its crest, and the pymrm ash shell |
| `test_model.py` | fast checks of `model.py` and one stored claim |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `agreement.json` | the claims, format 2, written by the last cell |

No `data/` directory: this is a tier 6 analytic page and has no experimental
content. That is stated on the page rather than papered over.

## Why this paper and not a textbook

Every textbook gives three limiting conversion–time laws for the shrinking core.
This paper gives the **single equation all three come from**, Eq. 6, and it is
rarely quoted. Reproducing it is the point of the page.

The scan's text layer is unusable for equations (it renders `θ_B` as `0B` and
drops the exponents in Eq. 6's numerator), so both equations were read off 600 dpi
renders.

## The check that establishes the reading

Eq. 6 was read from a page image, so it needs verifying. Two independent routes:

1. It is exactly 0 at *r*/*R* = 1 and exactly 1 at *r*/*R* = 0 for **any**
   (ω, γ), which only happens if the numerator collapses to the denominator term
   by term, i.e. only for the printed coefficients.
2. Integrating the moving boundary from Eq. 5's three resistances, without ever
   looking at Eq. 6, gives the same function to round-off, 6.9e-16
   <!-- claim: eq6_vs_integrated_worst -->, over twelve decades of *k*_c/*k*_f. That
   route independently recovers the factor 3 on the film term and the 12 in
   *k*_d1 = 12𝔻/*D*_p; the break table drops each in turn and the claim fails.

Eq. 6 also collapses to all three classical laws in the right limits: the gap
is proportional to the small parameter, 8.2e-10 <!-- claim: eq6_to_film_limit -->
for γ = 1e-9 (film), 2.3e-8 <!-- claim: eq6_to_ash_limit --> for ω = 1e7 (ash) and
3.8e-10 <!-- claim: eq6_to_reaction_limit --> for ω = 1e-9 (reaction).

## The pymrm shell solve

The quasi-steady gas profile in the ash shell is solved by pymrm (spherical
`construct_div`, Robin condition on each face) and the flow to the core is read as
the discrete gradient on the core face with `compute_boundary_values`. It meets the
closed-form series resistance to 2.1e-6 <!-- claim: shell_solve_worst_rel_err -->
at 400 cells on six cases (thin to thick shell), with observed orders 1.90
<!-- claim: shell_order_finest_min --> to 1.97 <!-- claim: shell_order_finest_max --> on
the finest triple of a 25 to 400 cell ladder (the thick shell approaches 2 from
below). The solve conserves the flow: the flow to the core equals the flow delivered
by the film to round-off, 3.3e-12 <!-- claim: shell_conservation_gap -->, and the
finest ladder difference is at least 1e3 times that round-off.

## What the page adds

A regime map. For every (ω, γ), how far is the closest single-resistance law from
the full equation? Over six decades in each group, some single law is within 5 %
on 70.8 % <!-- claim: frac_within_5pct --> of the plane and within 2 % on
51.3 % <!-- claim: frac_within_2pct -->; the worst error on the grid is
20.5 % <!-- claim: worst_single_law_error -->. The worst case itself lies on a crest
in γ near 0.855 <!-- claim: crest_gamma_omega1e6 --> where the film law and the
reaction law are equally far from Eq. 6, a kink that a grid in γ undersamples. The
crest error rises with ω: 20.49 % <!-- claim: crest_error_omega1e3 --> at the
map edge ω = 1e3 and 20.49 % <!-- claim: crest_error_omega1e6 --> at ω = 1e6, the
limit ω → ∞ being the worst case. There the reaction resistance is under 0.1 %
<!-- claim: crest_reaction_share_omega1e3 --> of the total, so along the crest only
the film and the ash layer matter. The map is a grid count and is computed on a
ladder of grids.

## The sign trap

Both shell boundary conditions use the **outward** normal. At the core face that
normal points toward the centre, so the reacting condition is
`{"a": D, "b": kc, "d": 0.0}`, structurally the same as the film condition, not
opposite to it. Write it with the inward normal and the core appears to emit gas.

## Not computed, deliberately

`k_c1` cannot be recovered from the paper: Parker and Hottel's correlation is
given for the specific combustion rate `K_c`, and the unit conversion to `k_c1`
is not printed. So no absolute burnout times appear, and the paper's Figure 6 is
not reproduced. Everything is in ω, γ and θ/θ_B, which is what Eq. 6 needs.

## Corrections

- 2026-11-02: said some single law is within 2 % on 52 % of the six-decade
  plane. The count was taken on a 90 × 90 grid (51.56 %); on 180, 360 and 720
  points per side it is 51.42, 51.35 and 51.30 % <!-- claim: frac_within_2pct -->,
  so the fraction is 51 %. The 5 % fraction (71 %) moves from 71.2 to
  70.8 % <!-- claim: frac_within_5pct --> and stays 71 % to two digits.
- 2026-11-02: said the worst case of the map is at ω = 11, γ = 0.8. That is the
  maximum of the 90 × 90 grid, taken on an undersampled kinked crest. The grid
  argmax moves with the grid (ω = 11, 25, 52 and 128 on 90, 180, 360 and 720
  points per side) and is not the worst case. Maximising over γ directly, the
  crest error rises with ω to 20.49 % <!-- claim: crest_error_omega1e6 --> at
  ω = 1e6, γ = 0.855 <!-- claim: crest_gamma_omega1e6 -->; the supremum is the
  limit ω → ∞. The level stays 20 % to two digits, 20.5 % on the finest grid.
