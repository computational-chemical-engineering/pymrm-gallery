# B3.2 — Szekely–Evans grain model for gas–solid reactions with a moving boundary

From Szekely, J. & Evans, J. W., *A structural model for gas–solid reactions
with a moving boundary*, Chemical Engineering Science **25**(6), 1091–1107
(1970), doi:10.1016/0009-2509(70)85053-9. Read from the 17-page scan at its
native 300 ppi; the text layer of this 1970 Pergamon scan is unusable for
mathematics and captions (mid-dot decimals vanish), so every equation and
every constant was transcribed from cropped native renders.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 37 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata |

No `data/` directory: tier 6 by necessity. The paper prints parameter sets
(figure captions, Table 1, F_T = 2.75) and closed forms, but no measurements;
every computed result (Figs. 5–19) exists only as figures, and C0 — which sets
the absolute time axis — is printed nowhere. Nothing is digitised and nothing
is fitted ("fit vs test" reduces to: everything is test; the page's only
least-squares operation is a labelled diagnostic line on its own computed
asymptote).

## What the page establishes

- **An exact closed form.** The quasi-steady, first-order grain model is
  integrable: with the local exposure ζ = ∫ψ dθ, the whole model collapses to
  E.P./ℓ = √(6X(θ)) with X a closed-form polynomial. Szekely & Evans solved
  this numerically in 1969; for their exact formulation the computer was
  avoidable. Both classical limits fall out, with the lag Δ = 1/4 + g/30.
- **The pymrm march reproduces it** to 3.4e-5 (g = 0.5) and 4.3e-5 (g = 5)
  over the whole history, at observed order 2.00 in the grid and 2.02 in the
  time step — orders measured *through* grain-exhaustion events, which the
  scheme handles as a complementarity condition (Eq. 21's printed 0 ≤ R′).
- **The model reduces to the published `B3.1`** at the exposed face: the
  simulated surface-grain history converges (order ~1.9) onto Yagi & Kunii's
  eq. 6 in its no-film limit with ω = g/6 — a closed form from a different
  paper, transcribed independently. The underlying algebraic identity between
  the two printed functions is recorded separately and labelled structural.
- **The diffuse zone, quantified.** The paper's central qualitative claim
  (reaction in a zone "of the same order as the size of the reacting
  specimen") becomes a number: a travelling zone of constant width w(g)·ℓ —
  0.295 cm for the printed Fig. 7 base case, from printed constants alone (no
  C0 enters). The 5 % window where neither textbook limit works spans a factor
  ~4 in time at g = 0.5 and more than four orders of magnitude at the paper's
  base case g = 500.

## The checks that can fail, and the ones that cannot

Every metric has a break-table row; the table recomputes **every** metric
under **every** defect, and the no-move cells are stated as blindness
findings. Labelled as structural (kept, not cited as evidence): the spatial
flux-telescoping identity (4.6e-15) and the G ↔ eq. 6 algebraic identity
(4.4e-13, below the CI floor). The temporal mass balance is the check with
power against the paper-specific hazard: Eqs. 21 and 27 print the same series
resistance in two different algebraic arrangements, and transcribing them
inconsistently moves that metric from 4e-6 to 2.5e-1.

The dimensional-closure check **fired during the build**: a first draft
truncated the domain at 6 ℓ and the check returned 48 % — the far-field
hazard that the break table separately demonstrates is invisible to every
interior residual once L is adequate.

## Scope, honestly

Part I's geometry is a **semi-infinite slab** (`nu=0`; the grains' sphericity
lives in the closure). The finite spherical pellet is Part II (*CES* **26**
(1971) 1901), which has no readable scan in this collection — its Elsevier
api-text drops decimal points and must not be used for constants. The paper's
own validation claim is reproducing "general trends" with "reasonable
values"; this page inherits exactly that epistemic status plus internal
mathematical consistency, and contains no experimental comparison.

## For `B3.4` (Sohn's additive reaction times)

Sohn's paper cites this one and bridges the shrinking-core and grain
pictures. Import `E_exact`, `G_of_xi`, `X_of_theta` from this page as the
exact reference for the grain-model side.
