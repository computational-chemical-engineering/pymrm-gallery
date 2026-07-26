# Handoff — state as of 2026-07-26

Start here if you are picking this up fresh. Read this, then
[`AGENTS.md`](../AGENTS.md), then [`pdf-findings.md`](pdf-findings.md).

## Built and live

**https://computational-chemical-engineering.github.io/pymrm-gallery/**

3 pages, 4 published catalog entries, both CI workflows green.

| Page | What it shows | Validation |
|---|---|---|
| `A4.9` Duncan–Toor ternary diffusion | Osmotic → uphill → diffusion barrier | **0.59 mole %** vs the paper's own 0.45%, 28 digitised points |
| `B1.1`+`B1.5` Thiele + Weisz–Hicks | η(φ), and 3 steady states at one φ | 2.2e-4 vs exact; both methods agree on all three branches |
| `F3.1` Hatta regimes | Enhancement factor, 3 regimes | 6.3e-3 vs exact; VKH good to 2.1%, DeCoursey to 8.7% |

Status counts live in `models.yaml`: 4 published, 25 planned, 2 deferred
(`H1.12`, `B1.12` — unpublished manuscripts, see the published-work-only policy
in [`blueprint.md §9`](blueprint.md#published-work-only-policy)).

**Only `A4.9` is validated against experiment.** The other two are provenance
tier 6 (exact/reference solutions) and say so on the page. Correcting that
balance is the single most valuable thing to do next.

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
(45 kB) — both of which had the *worst* OCR of the set. The API returns properly
encoded text, so exponents and subscripts survive.

**Limits:** Elsevier only. Xu & Froment and Krishna & Ellenberger are Wiley and
remain OCR-only. And full text never yields figure data — plots still need
digitising. Be polite with request volume; this is not a bulk-download tool.

## Extraction cost, revised

The API changes the ranking. Cheapest first:

| Page | Paper | Extraction | Value |
|---|---|---|---|
| `D2.2` | Van Welsenaere & Froment | **API, clean text** | Runaway criteria; sweep-based figures are striking |
| `A3.4` | Wakao & Funazkri | **API, clean text** | Sh–Re dataset, but likely a scatter *figure* → digitise |
| `E2.1` | Kunii & Levenspiel | good OCR (12.4k chars/page) | The canonical fluidised-bed model |
| `I1.2` | Oh & Cavendish | good OCR (9.0k) | Converter light-off, `S4`+`S7` |
| `C2.1` | Xu & Froment | **hard** — Table 6 needs page-image reading | Highest: most-used kinetics in the catalog |
| `F1.4` | Krishna & Ellenberger | tables clean, holdup in figures | 2,787 experiments |

## Recommended next moves

1. **`C2.1` Xu & Froment** despite the extraction cost — it is the most-used
   kinetics in the whole catalog and would be the second experimentally-validated
   page. Read Table 6 from a 600 dpi page image (the route that worked for Duncan
   & Toor Figure 2). **Do not repair the OCR by inference** — a mis-read exponent
   is a silently wrong rate constant. Three traps in that paper are recorded in
   [`pdf-findings.md`](pdf-findings.md#1-xu--froment-1989--page-c21).
2. **`D2.2` Van Welsenaere & Froment** — now cheap via the API, and the runaway
   boundary is a sweep over many solves, which makes a strong figure the original
   could only sketch.
3. **`F1.4` Krishna & Ellenberger** — transcribe Tables 1–2, digitise the holdup
   figures.

## Hard-won lessons — read before building

Both are also stored as user memories and will load automatically.

- **Never fabricate data.** If a dataset cannot be obtained, mark it
  `status: placeholder`, keep the page `status: planned`, and say so on the page.
  Re-simulating the model and presenting it as data is circular and would
  discredit the gallery.
- **Reported numbers must be deterministic.** On `B1.1` the fold was first
  located on a warm-start continuation curve; CI re-executed on another machine
  and got η_ignited 44.45 against 36.15 locally. Locate features on a smooth
  deterministic reference (analytical, or a shooting solve) and rank turning
  points by prominence. See
  [`pdf-findings.md`](pdf-findings.md) and the `B1.1` page README.
- **Validation catches what inspection does not.** Three defects in `B1.1` were
  found by the monotonicity and closure checks, not by reading the code: a sign
  error, a sweep that missed a solution branch entirely, and a comparison that
  reported 89% deviation where the truth was 1e-5.
- **`check_agreement.py` failing is a signal, not a nuisance.** Update the
  baseline only when the model genuinely changed, and say why in the commit.

## Working commands

```bash
cd ~/Code/pymrm_suite/pymrm-gallery
source ~/Code/pymrm_suite/.venv/bin/activate

python scripts/check_metadata.py            # metadata + provenance validation
python scripts/check_metadata.py --report-missing   # catalog IDs not yet in models.yaml
python scripts/run_pages.py                 # execute every page notebook
python scripts/check_agreement.py           # metric regression check

export QUARTO_PYTHON=$(which python)
quarto render                               # build the site into _site/
```

Pages are generated from builder scripts rather than hand-edited JSON; the
pattern is in the session scratchpad but is easy to recreate — build a list of
`nbformat` markdown/code cells and write the notebook. Nine required sections,
listed in [`AGENTS.md`](../AGENTS.md).
