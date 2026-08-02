# A1.1 — Ergun, Kozeny–Carman, Darcy–Forchheimer and the Eisfeld wall correction

One page, four correlations, one dataset — Ergun's own.

From Ergun, S., *Fluid flow through packed columns*, **Chem. Eng. Prog.**
48(2) 89–94 (1952), with the wall correction from Eisfeld, B. and Schnitzlein,
K., **Chem. Eng. Sci.** 56(14) 4321–4329 (2001),
[doi:10.1016/S0009-2509(00)00533-9](https://doi.org/10.1016/S0009-2509(00)00533-9).

This directory covers catalog entries **`A1.1`, `A1.2`, `A1.3` and `A1.4`**, which
is what `docs/catalog-A-foundations.md` asks for: *"One page overlaying Ergun,
Kozeny–Carman, Forchheimer, and the Eisfeld wall correction against a single Δp
dataset, with the wall-effect regime highlighted, is far more useful than four
pages."*

## Status: ready — the figure digitisation was reviewed and approved

The maintainer reviewed the numbered overlays on 2026-08-02: the ringed centres
are on real markers, including in the dense chain where markers merge into the
drawn eq. (14b) line, and nothing was picked up off the line, the graph paper or
the panel labels. A few markers per panel were missed. That is recorded as a
**stated limitation of the dataset** — `limitations.recall` in
`data/ergun-1952-fig7-markers.meta.yaml` — because it is a recall limitation:
it costs precision and distorts counts, and does not move the fitted constants
so long as the misses are scattered, which is what was reported.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 6 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, never the notebook |
| `data/ergun-1952-fig7-markers.csv` | **244** markers from Ergun's Figure 7, three panels |
| `data/ergun-1952-fig5-present-markers.csv` | 74 markers from Figure 5's top panel — the *same* runs, other ordinate, digitised independently as a cross-check |
| `data/ergun-1952-parameters.csv` | his constants, and the measured level of the line he labels "Kozeny–Carman" |
| `data/eisfeld-2001-wall-correction.csv` | Eisfeld & Schnitzlein Tables 2 and 3 |
| `data/eisfeld-2001-table1-sources.csv` | their Table 1 rows for the sources Ergun plotted — this is where D/dₚ comes from |
| `meta.yaml` | page metadata |

Every CSV has a `.meta.yaml` provenance sidecar. Read
`ergun-1952-fig7-markers.meta.yaml` before reusing the marker file.

## Three things to know before reusing this

**The paper has no text layer at all.** `pdftotext` returns four bytes for four
pages, and there are no tables anywhere in it. Every number on this page was read
off a 600 dpi render, and the dataset was digitised from figures. The catalog's
note that "Ergun's original data are tabulated in the 1952 paper" is wrong — it
is not.

**The round trip is the point, and it is not circular.** Ergun's constants were
least-squares fitted to these very points, so recovering them is the check the
paper pays for — but only if the extraction cannot be reading his drawn line
instead of his markers. It cannot: detection keys on the enclosed *white
interior* of an open circle, a line has none, and the half of the top panel that
contains the printed line and no markers returns zero points. That test is
repeated on the page from the CSV itself.

**"Least squares" needs a weight here.** f_v spans two and a half decades, so an
unweighted fit in Ergun's linear coordinates is dominated by a handful of
high-Reynolds points and returns k₁ = 168.5 (+12 %). The page quotes the fit with
equal *relative* weight, k₁ = 151.9, and prints both so the round trip is not
resting on a hidden choice.

## Headline results

| | |
| --- | --- |
| refit of k₁, k₂ from Ergun's own figure | **151.9 and 1.697** vs his printed 150 and 1.75 |
| per-panel refit (3 sources, 3 calibrations) | k₁ = 146.6 / 154.6 / 151.0 |
| Ergun's eq. vs the 244 points | bias +0.35 %, mean abs 5.16 %, rms 7.13 % |
| Fig. 5 ⇄ Fig. 7 cross-check | 2.9 % on binned median f_v, no trend |
| the line he calls "Kozeny–Carman" | f_v = 149.2 ± 0.5 — his 150, not Carman's 180 |
| Kozeny–Carman vs Ergun | +7.5 % at x = 10, −45 % at 100, −90 % at 1000 |
| Eisfeld (infinite bed) vs Ergun | +2.7 % viscous, **−25 % inertial** |
| Eisfeld at D/dₚ = 5, vs Ergun | +46 % on f_v at x = 10, −11 % at x = 1000 |
| the wall correction itself, D/dₚ = 1.6 | +178 % at x = 3, −24.5 % at x = 30000 — it changes sign |

## Regenerating

```bash
python build_page.py                       # rewrite index.ipynb
python - <<'EOF'
import nbformat; from nbclient import NotebookClient
nb=nbformat.read("index.ipynb",as_version=4)
NotebookClient(nb,timeout=1800,kernel_name="python3",resources={"metadata":{"path":"."}}).execute()
nbformat.write(nb,"index.ipynb"); print("OK")
EOF
```

The extraction itself is not in this directory — it needs the source scan, which
must never enter the repository. The method is described in full in the dataset
sidecars, in enough detail to reproduce from the PDF.
