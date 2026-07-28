# H1.7 — Wijmans–Baker solution–diffusion model

From Wijmans, J. G. & Baker, R. W., *The solution-diffusion model: a review*,
*J. Membr. Sci.* **107**(1–2) 1–21 (1995),
[doi:10.1016/0376-7388(95)00102-I](https://doi.org/10.1016/0376-7388(95)00102-I).

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 4 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `data/wijmans-baker-1995-fig5.csv` | 12 markers digitised from Figure 5, reviewed |
| `meta.yaml` | page metadata, validated against `models.yaml` |

## The argument

`A` is fitted to the water flux alone (Eq. 37) and `B` to the salt flux alone
(Eq. 40). Rejection is then not free — the permeate concentration is fixed by
the ratio of the two fluxes — so the third panel of Figure 5 is a **prediction**
from the first two, and that is what makes the figure worth reproducing.

Two results:

**The rejection panel cannot test the model it appears to test.** The prediction
is a rise of 0.31 percentage points across the measured range: about 4 px on a
figure whose printed curve is 6 px thick. The markers agree, and would have
agreed with any prediction above 99 %.

**Concentration polarisation, assumed away throughout the derivation, is
quantified and then ruled out.** The pymrm film solve gives a wall enrichment of
23 % for a 50 µm boundary layer at the highest measured flux. But a film thick
enough to explain the gap between the fitted osmotic pressure (322 psi) and the
stated one (~350) would also halve the apparent permeability, which the data
does not show.

## The one trap

The film's wall boundary condition uses the **outward** normal, which at the
downstream face points in +x, so the rejected-salt term enters with a minus
sign. Get it wrong and the wall concentration comes out *below* bulk — salt
depleting against a membrane it cannot cross.

## Open

The fitted osmotic pressure sits 28 psi below the caption's stated value, which
is four standard errors of the fit. Polarisation is ruled out; nothing is ruled
in. The page says so rather than explaining it away.
