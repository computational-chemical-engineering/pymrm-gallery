# A3.8 — The Onda correlations

Onda, K., Takeuchi, H. & Okumoto, Y., *"Mass transfer coefficients between gas
and liquid phases in packed columns"*, **Journal of Chemical Engineering of
Japan** **1**(1) 56–62 (1968),
[doi:10.1252/jcej.1.56](https://doi.org/10.1252/jcej.1.56).

Three coupled correlations — wetted area, liquid-side coefficient, gas-side
coefficient — with the area feeding the liquid side and not the gas side. Still
the default packed-column correlation set in process simulators.

> **Status: staged, ready.** Nothing is digitised, so nothing waits on a
> maintainer figure review. Every number comes from a printed numeral.

## How the paper was identified

The file on disk is `~/papers/pymrm-gallery/1_56.pdf`. A bare numeric filename
carries no metadata, so **the article was identified by reading it**: journal
page 56 carries the title, the byline *"KAKUSABURO ONDA, HIROSHI TAKEUCHI\*\* AND
YOSHIO OKUMOTO / Dept. of Chem. Eng., University of Nagoya, Nagoya"* and the
footnote *"Received on July 10, 1967"*; the running feet read *"VOL.1 NO.1
1968"* and the printed folios run 56–62, so the filename encodes volume 1,
page 56.

That reading was then **verified** against CrossRef for `10.1252/jcej.1.56`,
which returns the same title, the same three surnames, container *Journal of
Chemical Engineering of Japan*, volume 1, issue 1, pages 56–62, 1968. This is a
verified DOI, not an auto-resolution from a terse citation.

**Two other articles share the same seven-page PDF.** Page 56 opens with the
closing column of a preceding bubble-column paper (porous-plate distributors);
page 62 begins *"Gas absorption with chemical reaction in packed columns"* by
Onda, Sada and Takeuchi. Neither is used for its own content. The companion
article earns exactly one job: **it reprints eqs. (1)–(3) verbatim** in its
introduction, giving a second, independently typeset printing against which every
constant was checked.

**Resolution.** `pdfimages -list` reports all seven pages as CCITT-G4 bilevel,
2456 × 3330 px at **300 ppi native**. Every render is at 300 dpi; rendering
larger would only interpolate. Each equation and table row was **cropped to a
single line and re-read at that resolution** — necessary because the 1968
typesetting sets exponents at about a quarter of body size, and at page scale
`-0.05` / `0.05` and `-2.0` / `-2.6` are not separable.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 8 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, never the notebook |
| `data/onda-1968-correlation-constants.csv` | all fourteen printed constants and exponents of eqs. (1)–(3), plus the four of eqs. (4) and (5). **Every constant used in the notebook is read from here; none is typed** |
| `data/onda-1968-table1.csv` | Table 1 (page 57) — the authors' six CO₂-into-organic-solvent runs reduced to `k_L a = alpha L^n` |
| `data/onda-1968-stated-results.csv` | the scalar claims printed in the text: four error bands, `a_t D_p = 6(1-eps) = 3.4`, the eq. (4)/(5) validity range, the liquid loadings printed in Figs. 4 and 6, and the "0.61 of Norman" remark |
| `meta.yaml` | page metadata |

Each CSV has a `.meta.yaml` provenance sidecar with a `columns:` block.

## What the page establishes

**The eq. (4) identity — the headline.** Eq. (4) is eq. (3) rearranged into a
`j_D` factor for spheres, printed by the authors with no working shown. Recovering
it tests 5.23, 0.7, 1/3 and −2.0 against a fifth number printed independently of
them. Done two independent ways — a `sympy` reduction and a purely numerical
evaluate-and-regress over 4000 random property combinations — both giving
**0.774443** and **−0.300000** against the printed **0.771** and **−0.30**:
+0.447 % and 1.7e-16, with the two routes agreeing to 2.2e-16. The numerical
route also pins the 1/3 *structurally*: it samples `Sc_G` over 0.02–83 and the
one-variable fit has a worst residual of 3.6e-15, which is only possible if the
Schmidt number cancels exactly.

**A ceiling the correlation cannot exceed, which Onda's own table breaks.**
Eqs. (1) and (2) together force `dln(k_L a)/dlnL` into **(2/3, 0.800)** for any
packing, any liquid and any flow rate — the bound contains no prefactor, so no
refit can move it. Table 1 prints six measured exponents. **Two, both spheres,
sit above the ceiling**, the worst by +0.060, and the table's spread (0.16) is
wider than the whole window (0.1333).

**Two printed exponents are essentially untestable, and the page says so.**
Eq. (1)'s **0.75** is *exactly* inert whenever `sigma_c = sigma` — replacing it
by 7.5 leaves all twelve printed digits of the packed height — and this article
prints `sigma_c` for no material at all. Eq. (3)'s **1/3** has semi-elasticity
`−f_G ln(Sc_G)`, and `Sc_G ≈ 0.97` for air, so replacing it by 2/3 moves the
packed height by 0.73 %.

**Onda's eq. (6) as a pymrm BVP, not a quadrature.** A counter-current absorber:
two phases as two fields of one array with opposite velocities
(`v = [[G_M, -L_M]]`), inlets at opposite ends, a pointwise interphase source.
The elasticity of the packed height to every printed constant is then computed
three ways — analytic chain rule, finite differences on `Z`, and the same
differences pushed through the column solve — agreeing to 1.8e-7 and 1.4e-7.

## What it does not do, deliberately

- **No out-of-sample test of eqs. (2) and (3).** They were fitted to the data any
  comparison here would use, *including* Table 1's runs (page 57 says so). The
  exponent window is the only quantity a refit could not absorb, and even there a
  single fitted constant scores better on six rows with no error bars — reported
  beside it as the null baseline.
- **No test of `a = a_w`.** Onda's central assumption. The area never appears
  except multiplied by a coefficient, so nothing in the article separates them
  and nothing here can either. No break row exists for it.
- **No figure digitisation.** Figs. 1, 2, 3 and 7 are scatter plots of the fits;
  Fig. 8 is `Z_cal` against `Z_act`; Fig. 9 is the resistance split for
  distillation. All are out of scope and none is needed.
- **No comparison with A3.9 (Billet–Schultes) or A3.10 (Rocha–Bravo–Fair).**
  Neither source is on disk. That three-way comparison belongs to whichever of
  those two is built second. **A3.6 and A3.7** are stirred-tank correlations with
  a different geometry and a different `a`; nothing here transfers and none of
  their content is duplicated.
- **The illustrative absorber is not Onda's.** Air and water at 25 °C are
  standard-table properties; `sigma_c/sigma = 0.85`, `m = 1` and the 95 % duty
  are page choices, all labelled as such. The one thing taken from the paper is
  the geometry — `a_t D_p = 6(1-eps) = 3.4` for spheres, page 58, the only
  packing geometry the article supplies, which is why every calculation uses
  spheres.

## The A2.6 trap, measured rather than avoided

The gas outlet is read through `compute_boundary_values`, the same reconstruction
the flux operator used: order **1.990**, 1.90e-6 at n = 800, and the solute
balance closes to 4.4e-16. Reading it at the **centre of the last cell** instead
gives 7.29e-4 at the same n — **384× worse** — at order 1.007, and breaks the
balance. Both are printed side by side, because the defective reading looks
perfectly monotone and converged.

## Break table

65 rows over 22 metrics; **every metric has at least one row that moves it**.
Seven of the 46 reported metrics fall below `check_agreement.py`'s
`ABS_FLOOR = 1e-12` and are therefore not compared by CI at all — the page lists
them explicitly as pinned rather than proven.
