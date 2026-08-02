# A4.7 — Maxwell–Stefan diffusion in zeolite micropores

Reproduces the single-component core of Krishna, R. & Baur, R. (2003),
*Modelling issues in zeolite based separation processes*, Separation and
Purification Technology **33**(3) 213–254
([doi:10.1016/S1383-5866(03)00008-X](https://doi.org/10.1016/S1383-5866(03)00008-X)):
the Maxwell–Stefan flux law for a sorbate against a stationary framework, the
thermodynamic correction factor `D = eth * Gamma`, the two loading scenarios the
review distinguishes, and the transient uptake problem it solves in a spherical
crystallite.

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page.
- `data/krishna-baur-2003-sorption-parameters.csv` — the review's printed
  *inputs*: the DSL isotherm parameters of Tables 1 and 2, the benzene
  parameters printed inside the caption of Fig. 9, the jump frequencies and MFI
  cell dimensions behind eq. (23), and the two diffusivities of Section 4.1.
- `data/krishna-baur-2003-printed-results.csv` — the review's printed
  *outputs*, which are the targets the page is scored against, including **both**
  of the two values it prints for one orientation-averaged diffusivity.

**Tier 6 — reproduction, not validation.** Every target on this page is a number
Krishna & Baur themselves computed; none is a measurement. The review does
contain experimental comparisons, but every one inside this page's scope lives
only in a figure (Garg & Ruthven's ethane uptake, Shah's benzene diffusivities,
Millot's Arrhenius data), digitising a figure needs a maintainer review that is
not available, and nothing here rests on one. Reproducing these numbers
establishes that the page implements the review's equations and reads its
parameters correctly. It says nothing about whether the Maxwell–Stefan
description of zeolite diffusion is right. The distinction is kept in those words
throughout the notebook and in `meta.yaml`.

**What the page does not reproduce.** The review's Fig. 10(a) reports an
inflection in the 3MP desorption kinetics; neither of the review's two stated
loading dependences of `eth` produces a turning point in the desorption rate
here, the review does not say which one 3MP follows, and the claim lives in an
undigitised figure. The page reports the non-reproduction and does not resolve
it.

**One headline has no printed target behind it.** The "two extrema in Gamma"
survey rests on ten hand-transcribed Table 1 rows; only the C1 and nC4 rows are
certified by anything (the permeation flux ratio, 0.17 %), and the one diagnostic
that touches those rows, `GAMMA_IDENTITY`, differentiates the same isotherm
object it tests and so is blind to a wrong parameter by construction. The page
says so and measures the exposure with a perturbation table instead of implying a
validation it cannot have.

**Three printed defects are recorded, none silently fixed:** the sign of
eq. (25); `nu_zz + nu_zz` in the denominator of eq. (23c); and the symbol
`nu_str` printed twice in the sentence supplying the CH4 jump frequencies, where
the second must be `nu_zz`. The third is a repair by inference and is flagged as
one in the CSV row's note, in both sidecars and on the page.

**Source text.** The Elsevier full-text API, PII `S138358660300008X`, with the
authorised institutional key. The 2003 review is born-digital, so the
"API text is unusable for numbers" warning in `docs/pdf-findings.md` — which is
about pre-1980 scans — does not apply; three independent internal checks confirm
the transcription and are on the page. No PDF or page image is stored.

**The catalogue's other reference is not this paper.** `A4.7` was previously
pointed at `10.1016/0009-2509(95)00102-B`, which is van den Broeke & Krishna
(note the author order; the catalogue reverses it), a packed-bed breakthrough
study on activated carbon and a carbon molecular sieve with no zeolite
experiment in it. See `docs/source-sweep-2026-08-02.md`.

Key numbers: eight printed targets reproduced, worst deviation 0.78 %; the DSL
isotherm gives the review's own worked loadings 5.5958 against 5.596 and 0.09282
against 0.093; the permeation flux ratio 10.949 against 34/3.1 = 10.968, and
dropping the thermodynamic factor moves that by 3.73×; the pymrm sphere agrees
with the exact series to 1.21e-5 after the time error is extrapolated out. A
ninth printed number, the running text's orientation-averaged diffusivity, is out
by 4.2 % and is reported as a printed inconsistency rather than scored.

The page's own finding is the 4.1 % disagreement between the review's two
printed values of one quantity, resolved as a factor of exactly two in its
printed `eth_z`. The review applies eq. (23) a second time, to CH4, and that
number — 2.75e-7 m²/s, printed with no `eth_z` beside it — is reproduced to
0.25 % by the halved relation and missed by 5.0 % by the relation as printed. So
the review's *practice* is the halved relation in both of its worked examples and
Table 3 is the odd one out. Where the factor originates, in Kärger's 1973
original or in the review's arithmetic, the page does not say: Kärger is not on
disk.

The second finding is that "the M–S diffusivity is the loading-independent one"
holds in one of the review's two scenarios and fails by the same factor in the
other. Both are reproduced; the review states both in consecutive paragraphs and
contrasts them explicitly twice, so what is new is only that they are quantified
against each other on one axis. The factor of ten is algebraically forced —
Γ = 1/(1−θ) and eth = (1−θ) are reciprocals — and its size is set entirely by the
choice θ_max = 0.9. The direction is the finding; the magnitude is arithmetic.

**Scope.** Single component only. Mixture diffusion, the exchange coefficients,
IAST and the membrane selectivities are case `H1.9`; the argument for cutting
there rather than at micropore/membrane is on the page and in both case files.
Runtime ~105 s.
