# Gallery blueprint — structure, publishing, and agent-friendliness

Proposed shape of the `pymrm-gallery` repository. Written so that the first page
can be built without re-litigating any of these decisions.

---

## 0. Decision record: Quarto as the platform

**Decided (2026-07-25): Quarto.** Confirmed by the project owner. The reasoning
is recorded here so it does not get re-argued.

The two obvious reasons — `.ipynb` renders natively, and the output looks good —
are correct but are not the decisive ones. Three others matter more over the life
of a 266-entry gallery:

1. **`freeze`.** Quarto caches executed notebook output, so rebuilding the site
   after a typo fix does not re-run forty coupled reactor models. Without this,
   site build time grows linearly with the gallery and editing prose becomes
   painful. This is the single most important feature for a gallery of
   *expensive* notebooks.
2. **Citations and cross-references.** BibTeX support, `@xu1989` style citations,
   and numbered cross-refs are load-bearing for a literature gallery. A
   documentation-oriented static-site generator would need this bolted on.
3. **Listings as a first-class concept.** The gallery view, category filtering,
   and per-page YAML metadata come from the same mechanism, which is what lets
   the catalog metadata drive the site rather than being duplicated in it.

A fourth, smaller benefit: multi-format output means a citable PDF "handbook"
of the whole gallery is close to free once the site exists.

### Cost accepted

The rest of the suite (`pymrm-book`, `pymrm-book-teacher`) uses **MyST /
Jupyter Book** (`myst.yml`). Choosing Quarto means the suite carries two
publishing toolchains, with different CI, different contributor instructions,
and no shared config. That is a real cost and it was accepted deliberately: a
book is linear and TOC-driven, a gallery is non-linear and facet-driven, and
Quarto's listings are materially better for the latter. Revisit only if the
gallery turns out to be read like a book.

### Consequence 1 — keep notebook Markdown portable

The pymrm model style guide (§3.1) requires notebooks to render correctly in
Jupyter, GitHub preview, and VS Code, and therefore bans MyST-only directives.
**The same discipline applies to Quarto syntax.** Quarto's callouts
(`::: {.callout-note}`), shortcodes, and div attributes are exactly as
non-portable as the MyST directives the style guide already prohibits, and the
gallery's whole promise is that the `.ipynb` is downloadable and usable
standalone.

Rule: Quarto-specific syntax is allowed in `.qmd` wrapper and index pages, and
in YAML front matter. It is **not** allowed inside a page's `index.ipynb`
Markdown cells. Page notebooks use plain Markdown, `$$…$$` maths, and standard
images only.

### Consequence 2 — faceted filtering needs prefixed category tags

The taxonomy is deliberately two-dimensional (physical domain × structure code),
plus tier, priority, and data tier. Quarto listings give one category filter
plus text search, not true multi-axis faceting.

Cheap solution, adopted: flatten every facet into prefixed tags in the page's
`categories` field, so a single category cloud filters on all axes.

```yaml
categories: [sec:D, struct:S8, struct:S3, tier:T0, data:tier2, phase:gas-solid]
```

If that proves too coarse once ~50 pages exist, replace the listing index with a
custom page that renders a filterable table from `models.json`. Do not build
that until the need is demonstrated.

### Consequence 3 — one source of truth for metadata

`models.yaml` (§4) and per-page YAML front matter describe the same facts and
*will* drift if both are hand-maintained.

Adopted: **`models.yaml` is canonical.** It must be, because `status: planned`
entries have no page and therefore no front matter. CI generates each page's
`categories` from it and fails if a page's front matter contradicts its record.
`models.json` is emitted at build time for agent consumption.

### Consequence 4 — freeze hides breakage, so CI must run notebooks separately

With `freeze` on, a green site build says nothing about whether the notebooks
still execute. The `execute.yml` workflow (§6) is therefore not optional and
must run on a schedule as well as on PRs — otherwise a `pymrm` release can break
every page while the site keeps rendering cached output that looks fine.

### Operational notes

- Pin the Quarto version in CI (`quarto-dev/quarto-actions/setup` with an
  explicit version). Quarto is a single binary, so this is cheap and it keeps
  renders reproducible.
- Commit `_freeze/` so contributors and CI share the cache and PRs only re-run
  what actually changed.
- `execute-dir: project` so notebook-relative data paths behave identically
  locally, in CI, and on Colab.

---

## 1. Repository layout

```
pymrm-gallery/
├── README.md
├── AGENTS.md                  # how a coding agent should use this repo
├── _quarto.yml                # site config, listings, theme
├── index.qmd                  # landing page
├── models.yaml                # MACHINE-READABLE INDEX — the agent entry point
├── docs/                      # this planning material
│   ├── taxonomy.md
│   ├── catalog-A-foundations.md
│   ├── catalog-B-reactors.md
│   ├── catalog-C-adjacent.md
│   ├── data-strategy.md
│   └── blueprint.md
├── templates/
│   ├── page-template.ipynb    # the canonical page skeleton
│   └── meta-template.yaml
├── shared/
│   └── gallery_utils.py       # data loading, parity plots, Colab detection
├── pages/
│   ├── C2.1-xu-froment-smr/
│   │   ├── index.ipynb        # the page (Quarto renders this)
│   │   ├── meta.yaml          # structured metadata
│   │   ├── data/
│   │   │   ├── fig4-rates.csv
│   │   │   └── fig4-rates.meta.yaml
│   │   └── README.md          # short, for GitHub browsing
│   ├── A4.9-duncan-toor/
│   └── ...
└── .github/workflows/
    ├── execute.yml            # run every notebook on PR
    └── publish.yml            # render + deploy to GitHub Pages
```

One directory per page, named `<catalog-id>-<slug>`. The catalog ID is the join
key between `models.yaml`, the docs, and the page — never rename it.

## 2. The page template

Sections in fixed order, matching the page contract in
[taxonomy.md](taxonomy.md):

1. **Header** — title, catalog ID, structure codes, tier, one-sentence summary.
2. **Background** — what the model is for, why it mattered, who uses it.
3. **The published model** — original equations in original notation, plus a
   symbol → code-variable table.
4. **Parameters and assumptions** — one cell, all constants named and united.
5. **The data** — load the CSV, print provenance, plot the raw points.
6. **PyMRM implementation** — class-based per the pymrm style guide.
7. **Results** — pymrm overlaid on data, with a quantitative agreement metric.
8. **Validation** — analytical limit, grid independence, or conservation check.
9. **What pymrm adds** — the improvement, or an explicit statement that this is
   a faithful reproduction with no extension.
10. **Reuse** — "to adapt this to your system, change X, Y, Z" + links to sibling
    pages sharing structure codes.

Sections 1–9 mirror the pymrm model style guide's notebook section order, so
existing tutorials and teacher solutions can be lifted in with light editing.
Section 10 is new and exists specifically for the "adapt it to my problem" user.

## 3. Colab compatibility

Every notebook must run top-to-bottom on a fresh Colab VM with no local files.
Two requirements:

**Cell 1 is always the environment cell:**

```python
# Runs everywhere; no-op outside Colab
try:
    import pymrm
except ImportError:
    %pip install -q pymrm
```

**Data loads by URL, not relative path.** `shared/gallery_utils.py` provides:

```python
from gallery_utils import load_data   # falls back to raw.githubusercontent.com
df = load_data("fig4-rates.csv")      # works locally AND on Colab
```

An "Open in Colab" badge on every page, generated from the path — Quarto can
inject it from `meta.yaml` so it is never hand-maintained.

Pin `pymrm` to a minimum version in each page's `meta.yaml` and have CI test
against both the pinned floor and current `main`, so pages fail loudly when the
library moves rather than silently producing wrong figures.

## 4. `models.yaml` — the machine-readable index

This is the single most important file for agent-friendliness. One record per
catalogued model, whether or not a page exists yet:

```yaml
- id: C2.1
  slug: xu-froment-smr
  title: Steam methane reforming kinetics (Xu & Froment)
  section: C
  status: published        # planned | in-progress | published
  tier: T1
  priority: P1
  structures: [S1, S2]
  phases: [gas, solid]
  reference:
    authors: [Xu J., Froment G.F.]
    year: 1989
    doi: 10.1002/aic.690350109
  predicts: [reaction rate, conversion, selectivity]
  pymrm_api: [NumJac, newton, construct_convflux_upwind, construct_div]
  data:
    tier: 2
    method: table
  page: pages/C2.1-xu-froment-smr/
  related: [C2.2, D3.1, B1.11]
```

Publish it as JSON too. An agent asked to "build a fluidised bed model with
interphase mass transfer" greps `structures: S7` and `section: E`, finds
E2.1, and copies a working page rather than inventing an assembly from the API
docs.

## 5. `AGENTS.md` — instructions for coding agents

A short, imperative file at the repo root. Contents:

- Read `models.yaml` first; do not crawl the site HTML.
- The catalog ID is the stable identifier; slugs and titles may change.
- To build a new model: find the nearest entry by `structures`, copy that page
  directory, then substitute physics. Do not start from a blank notebook.
- The pymrm conventions that must be followed (link to the style guide; the
  outward-normal BC convention is the most common agent error).
- Never fabricate data. If no dataset exists, mark the page `status: planned`
  rather than inventing numbers.
- How to add a page: required files, required metadata fields, CI checks.

Also add `llms.txt` at the site root pointing to `models.yaml` and `AGENTS.md`,
since that is becoming the convention for agent-readable sites.

## 6. Continuous integration

- **Execute** every notebook on every PR. The `pymrm` repo already has
  `scripts/run_notebook.py`, `run_examples.py`, `run_tutorials.py` — reuse them.
- **Validate metadata**: every page has `meta.yaml`; every CSV has a
  `.meta.yaml`; every `meta.yaml` ID exists in `models.yaml`.
- **Check the agreement metric** hasn't regressed: each page writes its RMSE (or
  equivalent) to a small JSON; CI compares against a stored baseline. This is
  what stops a pymrm change silently degrading forty pages.
- **Link check** and **DOI resolution check**.
- **Publish** to GitHub Pages on merge to `main`.

## 7. Licensing

Three-way split, stated in the README and in `_quarto.yml`:

- **Code** (notebooks, `shared/`): MIT or Apache-2.0 — match whatever `pymrm`
  itself uses, for frictionless reuse.
- **Prose and figures**: CC-BY-4.0.
- **Data we generate** (simulation outputs, our own measurements): CC0 or CC-BY.
- **Data digitised from third-party papers**: not relicensed. The `.meta.yaml`
  records the source and the basis for extraction (see
  [data-strategy.md](data-strategy.md)); the CSV is distributed as factual data
  with attribution, and the source figure is never reproduced.

## 8. Contribution model

The catalog is ~260 entries; no single group will build it. Design for outside
contribution from the first commit:

- A `page-template.ipynb` that already passes CI when filled in.
- A GitHub issue template per proposed model, pre-populated from `models.yaml`
  (`status: planned` entries become "good first contribution" issues).
- Explicit authorship on every page, with ORCID, so contributors get citable
  credit.
- A `CITATION.cff` for the gallery as a whole, and a Zenodo DOI per release, so
  contributing is academically worthwhile.
- Consider a lightweight review checklist (physics correct, data provenance
  complete, notebook executes, validation present) rather than full peer review.

## 9. Suggested first six pages

Chosen for data availability, structural diversity, and in-house ownership —
they collectively exercise `S1`–`S9` and prove every part of the workflow.

| Order | Page | Why first |
|---|---|---|
| 1 | **A4.9 Duncan–Toor ternary diffusion** | Small tabulated dataset; Fick fails, Maxwell–Stefan succeeds; the clearest "why pymrm" story |
| 2 | **B1.1/B1.5 Thiele + Weisz–Hicks** | Analytical validation; exercises `S3` and nonlinear multiplicity; teacher solution exists |
| 3 | **C2.1 Xu–Froment SMR kinetics** | Tier-2 data; the most-used kinetics in the catalog |
| 4 | **H1.12 Ammonia membrane reactor** | In-house, published, data owned — proves the gallery handles research-scale models |
| 5 | **F1.4 Krishna–Ellenberger holdup** | Abundant tabulated data; a correlation page, testing that format |
| 6 | **D1.4 Fixed bed + particle coupling** | The `S8` flagship; shows what pymrm does that spreadsheet-level tools cannot |

After these six, the workflow (data provenance, template, CI, Colab) is proven
and the catalog can be worked through in priority order by more than one person.
