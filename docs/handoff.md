# Handoff — state as of 2026-07-26

Start here if you are picking this up fresh. Read this, then
[`AGENTS.md`](../AGENTS.md), then [`pdf-findings.md`](pdf-findings.md).

## Built and live

**https://computational-chemical-engineering.github.io/pymrm-gallery/**

4 pages, 5 published catalog entries, both CI workflows green.

| Page | What it shows | Validation |
|---|---|---|
| `A4.9` Duncan–Toor ternary diffusion | Osmotic → uphill → diffusion barrier | **experimental** — 0.59 mole % vs the paper's own 0.45%, 28 digitised points |
| `C2.1` Xu–Froment steam reforming | The most-used SMR kinetics, against the runs they were fitted to | **experimental** — 0.0017 in conversion (2.7%) over 61 digitised points, nothing fitted |
| `B1.1`+`B1.5` Thiele + Weisz–Hicks | η(φ), and 3 steady states at one φ | 2.2e-4 vs exact; both methods agree on all three branches |
| `F3.1` Hatta regimes | Enhancement factor, 3 regimes | 6.3e-3 vs exact; VKH good to 2.1%, DeCoursey to 8.7% |

Status counts live in `models.yaml`: 5 published, 24 planned, 2 deferred
(`H1.12`, `B1.12` — unpublished manuscripts, see the published-work-only policy
in [`blueprint.md §9`](blueprint.md#published-work-only-policy)).

**Two of four pages are now validated against experiment**, up from one of
three. Keeping that ratio moving is still the most valuable thing to do next.

## Papers available

Nine PDFs at **`~/papers/pymrm-gallery/`** (priorities 1 and 2 complete).
Inventory with per-file text-layer quality is in
[`pdf-findings.md`](pdf-findings.md).

## Elsevier full text — use this instead of OCR

An institutional API key is stored at **`~/.config/elsevier/apikey`** (mode 600,
outside the repo; `.gitignore` blocks credential filenames — never commit it).

The **DOI endpoint 404s** on older Chemical Engineering Science articles because
their DOIs contain parentheses. **Use the PII endpoint**, and note the PII is
recoverable straight from the ScienceDirect filename:
`1-s2.0-`**`0009250970850734`**`-main.pdf`.

```bash
K=$(cat ~/.config/elsevier/apikey)
curl -sS -H "X-ELS-APIKey: $K" -H "Accept: text/plain" \
  "https://api.elsevier.com/content/article/pii/0009250970850734"
```

Verified working on Van Welsenaere & Froment (47 kB) and Wakao & Funazkri
(45 kB) — both of which had the *worst* local OCR of the set.

**Correction, 2026-07-26 — do not trust the API text for numbers.** An earlier
version of this file claimed the API "returns properly encoded text, so
exponents and subscripts survive". That is wrong for pre-1980 scans. What comes
back is the *publisher's* OCR of the same scan, and on Van Welsenaere & Froment
it drops decimal points wholesale:

| API text | Truth (600 dpi page image) |
|---|---|
| `R = 00125 m` | *R* = 0.0125 m |
| `M = 2948 kg/kmole` | *M* = 29.48 kg/kmole |
| `(p°),,, = 001353 atm` | (*p*⁰)ₗ = 0.01353 atm |
| `b = 19837` | *b* = 19.837 |
| `t,„ = 21 -818` | *t*ₘ = 21.818 |

The 1970 typesetting uses a mid-dot decimal separator, which the OCR discards.
So the API is excellent for **prose** — section structure, what each figure
shows, which assumptions are stated — and useless for **parameters**. Read
those from a 600 dpi page render, the same as for a Wiley scan. Treat any
integer with an implausible magnitude as a lost decimal point and go to the
image; never "fix" it by inference.

**Other limits:** Elsevier only. Xu & Froment and Krishna & Ellenberger are
Wiley and have no API route at all. And full text never yields figure data —
plots still need digitising. Be polite with request volume; this is not a
bulk-download tool.

## Extraction cost, revised

| Page | Paper | Extraction | Value |
|---|---|---|---|
| `D2.2` | Van Welsenaere & Froment | API for prose, **page image for every number** | Runaway criteria; sweep-based figures are striking |
| `A3.4` | Wakao & Funazkri | API for prose, page image for numbers | Sh–Re dataset, but likely a scatter *figure* → digitise |
| `E2.1` | Kunii & Levenspiel | good OCR (12.4k chars/page) | The canonical fluidised-bed model |
| `I1.2` | Oh & Cavendish | good OCR (9.0k) | Converter light-off, `S4`+`S7` |
| `F1.4` | Krishna & Ellenberger | tables clean, holdup in figures | 2,787 experiments |
| ~~`C2.1`~~ | ~~Xu & Froment~~ | ~~hard~~ | **done 2026-07-26** |

## Recommended next moves

1. **`D2.2` Van Welsenaere & Froment** — the runaway boundary is a sweep over
   many solves, which makes a strong figure the original could only sketch. The
   API gives the prose cheaply, but every number has to come off a page render
   (see the correction above). The parameter set has already been read and is
   recorded in [`pdf-findings.md`](pdf-findings.md#4-van-welsenaere--froment-1970--page-d22);
   the paper has no tables, only ten figures, so validation means digitising
   Figures 1–3 (profiles) or Fig. 8 (the Barkelew comparison).
2. **`F1.4` Krishna & Ellenberger** — transcribe Tables 1–2, digitise the holdup
   figures. Budget real time for the figures: see the marker-extraction note
   below.
3. **`E2.1` Kunii & Levenspiel** — best text layer of the set, and the canonical
   fluidised-bed model.

## Hard-won lessons — read before building

The first two are also stored as user memories and will load automatically.

- **Never fabricate data.** If a dataset cannot be obtained, mark it
  `status: placeholder`, keep the page `status: planned`, and say so on the page.
  Re-simulating the model and presenting it as data is circular and would
  discredit the gallery.
- **Reported numbers must be deterministic.** On `B1.1` the fold was first
  located on a warm-start continuation curve; CI re-executed on another machine
  and got η_ignited 44.45 against 36.15 locally. Locate features on a smooth
  deterministic reference (analytical, or a shooting solve) and rank turning
  points by prominence.
- **Validation catches what inspection does not.** Three defects in `B1.1` were
  found by the monotonicity and closure checks, not by reading the code: a sign
  error, a sweep that missed a solution branch entirely, and a comparison that
  reported 89% deviation where the truth was 1e-5.
- **`check_agreement.py` failing is a signal, not a nuisance.** Update the
  baseline only when the model genuinely changed, and say why in the commit.
- **Look for a check the *paper* pays for.** `C2.1` gained three that cost
  nothing: Tables 5 and 6 are related by the Arrhenius form, so recomputing one
  from the other tests the page-image reading *and* the split reference
  temperature at once; the van 't Hoff slope of an equilibrium correlation is a
  reaction enthalpy, so it can be checked against the paper's own table; and two
  figures plotting the same runs must pair up point-for-point across
  independently fitted axes. Look for these before writing the notebook — they
  are worth more than any amount of code review.
- **Marker extraction is per-figure, not a recipe.** Morphological opening
  worked on Duncan & Toor because its markers are solid and much larger than the
  curves. Xu & Froment's are ~20 px glyphs on ~10 px curves and needed a
  different method: local ink density minus the largest value explainable by a
  locally straight structure (grey-scale opening with a long *line* element,
  maximised over six orientations). Take the maximum over orientations, or every
  steep near-origin curve section is flagged as a marker. Then audit every
  candidate visually at 600 dpi — on `C2.1` roughly a quarter of the automatic
  candidates were curvature artefacts with no glyph at the crosshair.

## Working commands

```bash
cd ~/Code/pymrm_suite/pymrm-gallery
source ~/Code/pymrm_suite/.venv/bin/activate

python scripts/check_metadata.py            # metadata + provenance validation
python scripts/check_metadata.py --report-missing   # catalog IDs not yet in models.yaml
python scripts/run_pages.py                 # execute every page notebook (slow: ~minutes)
python scripts/check_agreement.py           # metric regression check

export QUARTO_PYTHON=$(which python)
quarto render                               # build the site into _site/
```

Pages are generated from builder scripts rather than hand-edited JSON.
`pages/C2.1-xu-froment-smr/build_page.py` is the reference example: a list of
`nbformat` markdown/code cells written to `index.ipynb`, so prose and code stay
in one reviewable file. Nine required sections, listed in
[`AGENTS.md`](../AGENTS.md).
