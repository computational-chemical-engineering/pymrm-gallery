# B1.2 — Aris generalised shape modulus

From Aris, R., *On shape factors for irregular particles — I. The steady state
problem. Diffusion and reaction*, Chemical Engineering Science **6**(6) 262–268
(1957), doi:10.1016/0009-2509(57)85028-3.

**Staged, not published.** `B1.2` has no entry in `models.yaml`. Promoting this
page means moving the directory to `pages/B1.2-aris-shape-factor/` and adding the
matching `models.yaml` entry. Nothing outside `queue_cases/B1.2/` has been
touched.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 34 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata (`status: staged`) |
| `data/aris1957-table1.csv` | his Table 1, 21 entries, + provenance sidecar |
| `data/aris1957-spherical-shell.csv` | his section 3 shell table, 5 entries, + sidecar |

Both datasets are the author's own computed values, not measurements. Aris (1957)
reports no experiments; this is provenance tier 6 and the page says so.

## Getting the paper — read the 1995 reprint, not the 1957 scan

Aris (1957) is closed access. It was retrieved through the Elsevier article
retrieval API (PII `0009250957850283`) under an institutional subscription. That
scan's text layer is the publisher's OCR of a mid-dot-decimal typesetting and is
useless for numbers: Table 1's 0.698 comes out as `898`, 0.432 as `482`.

**Chemical Engineering Science reprinted the paper verbatim in 1995** for its
Golden Jubilee — 50(24) 3899–3903, doi:`10.1016/0009-2509(96)81819-7`, PII
`0009250996818197` — and that re-typesetting OCRs cleanly. Both were retrieved
and every transcribed number is one on which the two agree. **This trick is worth
trying for any pre-1970 CES classic**; the Golden Jubilee issue reprinted a dozen
of them, and it is a cheaper route than page-image reading where it applies.

## What could NOT be obtained

The entitlement returns a **one-page preview PDF only**, for both printings
(journal pages 262 and 3899). Table 1 is on page 265 / 3901. So the gallery's
usual discipline of reading numbers off a 600 dpi render was not available.

What stands in for it: the paper's own algebra reproduces 19 of the 21 Table 1
entries to ≤ 0.0008, and all five spherical-shell entries to ≤ 0.0004. No
mis-transcription survives that. **Two entries are flagged rather than
corrected** — see `meta.yaml`'s `blocker` block.

## The check the paper pays for

The spherical-shell expression has six distinct terms and had to be read from a
mangled scan. It is over-determined seven ways:

1. it must collapse to eq. 17 (sphere) at *p* = 0 — the paper says so — and does,
   to 1.1 × 10⁻¹⁶;
2. it must collapse to eq. 15 (flat plate) at *p* = 1 — likewise — to
   5.8 × 10⁻¹⁰;
3. it must reproduce his five printed values at Λ = 1, and does, to 0.00041.

A single mis-read coefficient breaks all seven at once. As a by-product it
decides which of the two printings is right where they disagree: at *p* = 0.75
the 1957 OCR gives `728`, the 1995 reprint 0.723, and the expression gives
0.7227.

## The ambiguity only a solver can settle

Aris's case (iv), the finite cylinder with porous ends, is given as a double
Fourier–Bessel series that **he never evaluates**. Its `p` is defined two
incompatible ways in the same paper: the nomenclature says radius:length, the
running text says "length = 2*a*/*p*". They differ by exactly the factor of 2 to
which the axial eigenvalue is most sensitive, and since no number is printed,
nothing on the page breaks the tie.

A 2-D axisymmetric pymrm solve does: with `p` read as the nomenclature defines
it, the printed series matches to 4.2 × 10⁻⁴; with the running text's reading it
is out by 47 %. **The series is correct as printed and the running text is the
error.**

## What the page adds

Aris wrote that estimating the divergence in the general case "would be an
excessively difficult task, equivalent to solving the problem completely". The
page solves it over fourteen shapes — his three, four spherical shells, five
finite cylinders from long rod to flat disc, a square rod and a cube:

- the band is **0.116** wide at Λ ≈ 1.1, against his plate-to-sphere **0.092** —
  27 % wider;
- **the sphere is not the floor.** At the same *v*_p/*s*_x = *a*/3 the cube lies
  0.026 (3.7 %) below it, with the cylinder of length 2*a* between them. A corner
  admits reactant from two directions at once.

Both of Aris's own statements survive and are not in conflict: at *equal volume*
the sphere is the worst case (his Leva shape-factor argument), while at *equal Λ*
it is the best of the three shapes sharing *v*_p/*s*_x = *a*/3.

## The sign trap

The low end of every axis carries `{"a": 1, "b": 0, "d": 0}` — with the
**outward** normal that is ∂χ/∂n = 0, i.e. symmetry for a solid shape and a
*sealed* face for the shell. Mark the shell's inner face exposed by mistake and
you get a perfectly convergent solution to a different problem, with no warning.
The *p* → 0 and *p* → 1 limits are what catch it.

## Numerical limit worth knowing

Everything except the cube is cheap. In 3-D a direct sparse solve runs out of
resolution before it runs out of memory: the reacting layer is *O*(1/λ) thick, so
at Λ = 20 the 24³ and 32³ solutions differ by 11 %. The sweep therefore takes the
cube from its eigenfunction series beyond Λ ≈ 5, validated against the 3-D solve
at Λ = 1 (7.5 × 10⁻⁵, converging second order). Grading the grid toward the
exposed faces buys about a decade in Λ for free and is used for the eq. 18
asymptote.

## Relationship to B1.1

`B1.1` already plots the Aris collapse for Thiele's three geometries and reports
the residual spread. This page deliberately does not repeat that. It is about the
shapes that have no closed form, which is what needs a solver — and it
cross-links to `B1.1` rather than restating it.
