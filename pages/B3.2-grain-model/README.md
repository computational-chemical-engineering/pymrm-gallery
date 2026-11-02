# B3.2 Szekely–Evans grain model for gas–solid reactions with a moving boundary

From Szekely, J. & Evans, J. W., *A structural model for gas–solid reactions
with a moving boundary*, Chemical Engineering Science **25**(6), 1091–1107
(1970), doi:10.1016/0009-2509(70)85053-9. Read from the 17-page scan at its
native 300 ppi; the text layer of this 1970 Pergamon scan is unusable for
mathematics and captions (mid-dot decimals vanish), so every equation and
every constant was transcribed from cropped native renders.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb` |
| `index.ipynb` | the executed notebook; runs in about 95 s |
| `model.py` | the exact route (grain law, penetration law, profiles, zone width, departure times) and the pymrm march `GrainModel` |
| `test_model.py` | fast checks of `model.py` and one stored claim |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `agreement.json` | the claims, format 2, written by the last cell |

No `data/` directory: tier 6 by necessity. The paper prints parameter sets
(figure captions, Table 1, F_T = 2.75) and closed forms, but no measurements;
every computed result (Figs. 5–19) exists only as figures, and C0, which sets
the absolute time axis, is printed nowhere. Nothing is digitised and nothing
is fitted.

## What the page establishes

- **An exact closed form.** The quasi-steady, first-order grain model is
  integrable: with the local exposure ζ = ∫ψ dθ, the whole model collapses to
  E.P./ℓ = √(6X(θ)) with X a closed-form polynomial. Both classical limits fall
  out, and beyond the surface burnout the specimen-scale law
  E² = 6(θ − Δ) holds exactly, with Δ = 1/4 + g/30.
- **The pymrm march reproduces it** to 3.4e-5 <!-- claim: ep_vs_exact_g05_relmax -->
  (g = 0.5) and 6.2e-5 <!-- claim: ep_vs_exact_g5_relmax --> (g = 5) over every
  step of the history from θ = 0.1, and its profiles at the cell centres to
  2.9e-5 <!-- claim: profile_psi_dev_n600 --> (gas) and 1.7e-5
  <!-- claim: profile_xi_dev_n600 --> (grain) on 600 cells. Two fields share one
  grid through `BlockLayout.from_fields`, with an analytic Jacobian that agrees
  with `NumJac` to 9.0e-7 <!-- claim: jacobian_gap_vs_numjac --> (relative,
  mid-march, a measurement that a Jacobian missing a term fails); the march with
  `NumJac` reaches the same state. Observed orders are 2.00 in
  the grid and in the time step, measured *through* grain-exhaustion events,
  which the scheme handles as a complementarity condition (Eq. 21's printed
  0 ≤ R′).
- **The model reduces to the published `B3.1`** at the exposed face: the
  simulated surface-grain history converges (order 1.94
  <!-- claim: surface_b31_order_n600 -->) onto Yagi & Kunii's eq. 6 in its no-film
  limit with ω = g/6, a closed form from a different paper, transcribed
  independently.
- **The diffuse zone, quantified.** The paper's central qualitative claim
  (reaction in a zone "of the same order as the size of the reacting
  specimen") becomes a number: a travelling zone of exactly constant width
  w(g)·ℓ, 0.295 cm <!-- claim: zone_width_cm_fig7base --> for the printed Fig. 7
  base case, from printed constants alone (no C0 enters). The 5 % window where
  neither textbook limit works spans a factor 3.88
  <!-- claim: window_ratio_g05 --> in time at g = 0.5 and 5.27e4
  <!-- claim: window_ratio_g500 --> at the paper's base case g = 500.

## The checks that can fail, and the one that cannot

Every claim has a break-table row that injects one defect through an argument
of `model.py`, except the spatial flux-telescoping identity, which is labelled
structural. The temporal mass balance is the check with power against the
paper-specific hazard: Eqs. 21 and 27 print the same series resistance in two
different algebraic arrangements, and dropping the (1 − ξ) from one of them
moves that metric from 1.0e-6 <!-- claim: mass_balance_rel_dt5e3 --> up by orders of magnitude.
The dimensional closure on the printed Fig. 7 constants (2.3e-4
<!-- claim: dimensional_closure_relmax -->) fails when the domain is cut
to 6 ℓ, the far-field hazard that no interior residual sees; the far boundary
condition itself is invisible to every check once the domain is long enough.

## Scope, honestly

Part I's geometry is a **semi-infinite slab** (`nu=0`; the grains' sphericity
lives in the closure). The finite spherical pellet is Part II (*CES* **26**
(1971) 1901), which has no readable scan in this collection; its Elsevier
api-text drops decimal points and must not be used for constants. The paper's
own validation claim is reproducing "general trends" with "reasonable
values"; this page inherits exactly that epistemic status plus internal
mathematical consistency, and contains no experimental comparison.

## For `B3.4` (Sohn's additive reaction times)

Sohn's paper cites this one and bridges the shrinking-core and grain
pictures. Import `E_exact`, `G_of_xi`, `X_of_theta` from this page's
`model.py` as the exact reference for the grain-model side.

## Corrections

- 2026-11-02: said the march meets the exact law to 4.3e-5 at g = 5 over
  0.1 ≤ θ ≤ 20. The window dropped θ = 0.1, where the error is largest, because
  the times were accumulated in floating point; over every step from θ = 0.1 it
  is 6.2e-5 <!-- claim: ep_vs_exact_g5_relmax -->. The g = 0.5 value is
  unchanged.
- 2026-11-02: said the time step converges at observed order 2.02. The errors
  were taken against a finer run of the same code, which biases the order; from
  successive differences, which need no reference, it is 2.00
  <!-- claim: dt_order_observed -->.
- 2026-11-02: said the gas and grain profiles agree with the exact solution to
  2.4e-2 and 1.9e-2, limited by interpolation across the front corner. Both
  maxima were at η = 0, ahead of the first cell centre, an artefact of the
  interpolation. At the cell centres the profiles agree to 2.9e-5
  <!-- claim: profile_psi_dev_n600 --> and 1.7e-5
  <!-- claim: profile_xi_dev_n600 --> on 600 cells.
- 2026-11-02: said the first pair of the dimensional step ladder (order 0.82)
  was pre-asymptotic. The old ladder compared maxima taken at different times;
  from successive differences at a common time, which every rung lands on,
  the observed order is 2.00 <!-- claim: dimensional_order_observed -->.
