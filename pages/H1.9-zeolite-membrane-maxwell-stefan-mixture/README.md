# H1.9 — Maxwell–Stefan mixture diffusion in a zeolite membrane

Mixture diffusion in an MFI (silicalite) membrane, from Krishna, R. & Baur, R.,
*Modelling issues in zeolite based separation processes*, **Separation and
Purification Technology 33**(3) 213–254 (2003),
[doi:10.1016/S1383-5866(03)00008-X](https://doi.org/10.1016/S1383-5866(03)00008-X).

A 95–5 methane/n-butane mixture is fed to the membrane at 100 kPa. The
pure-component fluxes predict a butane-over-methane permeation selectivity of
1.73. The measured and computed value is two to three orders of magnitude
larger. This page implements the two reasons the review gives — a mixture
isotherm that needs IAST rather than the multicomponent Langmuir, and **finite**
exchange coefficients between the two sorbates — and scores both against the
review's numbers and against the measurements it quotes.

## What is here

| | |
|---|---|
| Sections of the review | 3.1–3.4 (mixture theory) and 4.2–4.3 (MFI membrane separations) |
| Equations implemented | (38)–(47) and (53)–(65) |
| Data tier | **2** — four *measured* selectivities, printed as numbers in the review's running text |
| Validation | Funke et al.'s measured 24 reproduced as 25.0 (+4.3 %), nothing fitted — inside a 1.6× band of three undisclosed-in-the-review readings, all swept on the page |
| Reproduction | nine printed steady-state targets, worst 0.74 %; the transient peak **time** misses by ~10 % and is reported as an open tension |
| Runtime | ~10 minutes — the transient runs at `n_z` = 120 with a grid study out to 240 and a three-level step study |

The single-component half of the same review — the flux law, the scalar
thermodynamic factor, Kärger's zero-loading relations, transient uptake in a
crystallite — is [`A4.7`](../A4.7-zeolite-micropore-maxwell-stefan/) and is not
repeated here. **The cut between the two cases is single-component / mixture,
not micropore / membrane**; the reasoning is on both pages under *Reuse*.

## Data

| file | tier | what it is |
|---|---|---|
| `bakker-funke-mfi-membrane-selectivities.csv` | **2** | Four **measured** permeation selectivities — Bakker's 380 and 60 for CH₄/n-C₄H₁₀, Funke et al.'s 24 and 1.3 for n-C₆/3MP. Stated as numbers in the review's running text, not plotted, so there is no digitisation and no digitisation error. Both origins recorded under `origin_not_consulted`: neither is on disk and neither was opened. |
| `krishna-baur-2003-mixture-printed-results.csv` | 6 | Every mixture result the review prints as a number. The authors' own computed output — comparing against it is reproduction. |
| `krishna-baur-2003-mixture-parameters.csv` | 6 | The parameters that appear only in the mixture sections (3.3, 3.4, 4.3). |

Tables 1 and 2, the MFI density and the Section 4.1 single-component results are
**not** duplicated here: they are loaded cross-page from `A4.7`. Every finding
`A4.7` states about the rows used is listed in *The data*, and **no number that
is a row in a borrowed CSV is *stated as a result* without this page's own value
being printed beside it and the two reconciled.** (The narrative markdown — the
title and the framing cells — does retype 1.73, 34, 3.1, 95 and 5 kPa, because a
markdown cell cannot interpolate. Each is loaded, printed beside this page's
value and reconciled in the cells that use it; the obligation is on results, not
on prose that names the numbers the reader is about to see.)

## Headline results

**Against measurement** (the only validation on the page):

| | measured | this page | model / measured |
|---|---|---|---|
| Funke, nC₆/3MP at 15 kPa | 24 | 25.0 | **1.04** |
| Bakker, 95–5 CH₄/nC₄ | 380 | 488 | 1.28 |
| Bakker, 50–50 CH₄/nC₄ | 60 | 203 | 3.38 |

With the exchange coefficients removed the same model gives 12.4 and 5.4 for
Bakker's two feeds — 31× and 11× *too low*. So the measurements discriminate
decisively between the two implementations while agreeing quantitatively with
neither. The review calls the first comparison "quite close" and quotes the
second without comment.

**Reproducing the review's own computations:** all five printed permeation
selectivities (487, 12.3, 202, 5.4, 43.2) to better than 0.75 %, and all four
printed mixture fluxes to 0.41 % — but see the caveat below on what the flux
comparison is *not* a test of.

## The transient peak time, and the membrane thickness

The review never prints the membrane thickness. `A4.7` reconstructed
δ ≈ 40 µm by inverting the printed single-component methane flux. **The transient
methane peak time is the only quantity on this page that carries δ at all** — the
four "absolute" mixture fluxes look like a second, independent handle on it and
are not one, because δ, ρ_MFI, the unit-cell mass and N_A cancel out of them
identically (demonstrated on the page by swinging ρ_MFI over a factor of four and
getting bit-identical fluxes).

Refined **under the transient** — which the steady grid study never did, because
the transient carries a self-sharpening butane front and the steady profile does
not — the peak time lands about 10 % below the printed 0.73 s, and both further
grid and step refinement move it further away. Closing that needs δ about 5 %
larger, which nothing else here forbids. Reported as an open tension, with the
candidates named: ρ_MFI = 1800 read off a figure caption, the hard-coded
Si₉₆O₁₉₂ cell mass, the two-figure printed 34, or the review's transient not
being the calculation reproduced here.

The two peak **heights** move with the grid too: at `n_z` = 30 they sat 1.5 % and
3.5 % *below* the printed values, converged they sit about 1 % and 2 % *above*.

## Two printed numbers that do not reproduce

Both are investigated rather than explained away.

- **The Henry-limit sorption selectivity, printed as 2200.** Table 1's dual-site
  Henry coefficients give **2649.9**, and the IAST solve reaches the same limit
  from an entirely different code path — so it is not an arithmetic slip. The
  obvious explanation, a transcription error in one of the four C1/nC4 Table 1
  parameters, is *excluded*: every single-parameter change that would make the
  ratio 2200 moves `A4.7`'s certifying flux-ratio check by 4.2 % to 24 %,
  against the 0.17 % actually observed. Where 2200 comes from is left open — the
  review's Fig. 24 is drawn against CBMC fugacities, so a Henry coefficient
  taken from the simulation rather than from the fit remains possible.
- **The multicomponent-Langmuir counterfactual, printed as 800.** Its
  *structural* claim reproduces exactly — with equal saturation capacities
  $S_P$ is composition-independent to 1.02× over feeds from 95 % down to 5 %
  methane, where IAST spans more than a factor of three over the same five — and its *value* does not (about 1680). The review does not say which
  affinity it used once the dual-site fit is collapsed to a single site, and
  nothing here was tuned to reach 800.

## Notes for the next agent

- The review's conclusion for zeolite 4A points the **other way**: there
  $\eth_{ij}\to\infty$ fits better. The page carries that; a page showing only
  finite exchange coefficients winning would misreport it.
- Sections 4.4–4.7 (4A uptake, O₂/N₂ chromatography, C₅/C₆ breakthroughs) and
  Figs. 17–21 are scoped out: every one is validated only by a figure and the
  review prints no numerical result for any of them. Their parameter tables are
  clean, so they are viable future cases — a packed-bed chromatography page in
  section J is the natural home.
- Two printed defects are recorded on the mixture side and neither is repaired
  silently: eq. (37) prints `L_ij = L_ij` where the Onsager reciprocal relations
  require `L_ij = L_ji`, and the caption of Fig. 21 prints `Theta_1,sat` twice.
  **Both readings come from the Elsevier full text; there is no page image of
  this paper on disk and nothing here was read off a render.** Both involve only
  subscripts, which the text dump reproduces reliably.
- **There is no third defect.** An earlier version of this page claimed eq. (60)
  is *printed* inverted to the way it is *used*. It is not: the text dump
  flattens eqs. (60) and (61) identically, as compound fractions, the page reads
  eq. (61) that way, and the review's own worked value `3.1/34.0/(5/95) = 1.73`
  agrees with the same reading of eq. (60). The break-table row stays — it
  settles *which reading is meant*, by showing the alternative lands 361× away —
  but the claim about the printed document is withdrawn.
- **The four-flux comparison cannot fail on δ.** It is the δ-free flux-ratio
  check rescaled by the stored 34; both deviations are printed side by side on
  the page so the coincidence is visible. It is kept, and labelled, because it
  does test the mixture physics.

## Reproduce

```bash
cd pages/H1.9-zeolite-membrane-maxwell-stefan-mixture
python build_page.py     # regenerate index.ipynb from the builder script
jupyter nbconvert --execute --inplace index.ipynb
```

Nothing on the page is stochastic; the one random draw (200 states for the
friction-system check) is seeded. Two consecutive executions give identical
content, figures and `agreement.json`.
