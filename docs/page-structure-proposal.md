# Proposal: pages that are easy to download, run and reuse

Status: PROPOSAL, 2026-09-26, not implemented. The gallery queue is paused; this
note is for the session that resumes it. Written alongside the pymrm agent
plugin (pymrm repo, branch `feature/claude-plugin`), whose output-format rule
this proposal follows.

## 1. What a page is today (measured on 2026-09-26, 82 pages)

- One notebook per page, `pages/<id>/index.ipynb`, generated in 79 of 82 pages
  by `build_page.py`, which writes the cells as Python strings. The generator
  is the real source; the notebook is its output.
- Everything is inline. 50 of 82 notebooks define their model classes in cells.
  Code per notebook: median about 980 lines, largest about 3250. The model sits
  between the reproduction of the source, digitised data, break tables,
  agreement metrics and figures.
- Download works through the bootstrap cells: `%pip install -q pymrm` when
  pymrm is missing, and `shared/gallery_utils.py` and the page's data fetched
  from `raw.githubusercontent.com/.../main`. A single downloaded `index.ipynb`
  runs in Colab, as the README promises, given network access.
- Pinning: `requirements.txt` has `pymrm>=2.3`; CI runs every notebook against
  the pinned floor and the latest pymrm (good). The per-page minimum version in
  `meta.yaml` and the generated "Open in Colab" badge that the blueprint (section
  3) asks for do not exist yet.

## 2. What falls short

1. **The model is not reusable.** The gallery is for researchers "who want to
   adapt rather than rebuild", but reusing a model means excavating a class from
   a thousand-line notebook. There is nothing to import.
2. **Downloads drift.** The bootstrap installs the newest pymrm and fetches
   helpers and data from `main`. A later change can alter an old page's numbers
   without anyone looking at that page. The pymrm API guards now on
   `feature/api-guards` are an example: a new `NumJac` warning and a new stencil
   default for 1-D shapes with `axes_diagonals`.
3. **The model is buried in its validation.** A reader looking for the model
   meets the reproduction machinery first.
4. **The source is a generator.** Editing a page means editing Python strings
   that become cells, which is hard for people and error-prone for agents.

## 3. Proposal

### 3.1 A page is a folder with an importable model

```
pages/<id>/
  model.py        the reusable model: class(es) and functions, docstrings,
                  no plotting, no data loading, no page-specific studies
  index.ipynb     the page: narrative, drives model.py, data, validation,
                  figures, agreement metrics (as today, minus the model code)
  data/           as today
  meta.yaml       as today, plus `pymrm: ">=X.Y"` (minimum version)
  README.md       as today
  agreement.json  as today
```

This is the pymrm agent plugin's format rule for a full reactor model: `.py`
module(s) plus a notebook that drives and reports. Section 10 of the page
template ("Reuse") then becomes concrete: `from model import TanksInSeries`.

Pages whose model is a closed form or a few lines (many of sections A1, A4, J)
keep everything in the notebook; the rule applies where a model class exists.

### 3.2 Downloads that keep working

- The bootstrap cell resolves `model.py` the way it already resolves
  `gallery_utils.py`: the page folder when it exists, otherwise a download.
- Helpers, data and `model.py` are fetched from a **tag**, not from `main`
  (for example `gallery-2026.10`), written into the notebook when the page is
  published. A downloaded notebook then fetches the files it was published with.
- The bootstrap installs `pymrm` with the page's pinned specifier from
  `meta.yaml`, not the newest release.
- The publish workflow also builds a per-page zip (`index.ipynb`, `model.py`,
  `data/`, a copy of `gallery_utils.py`) that runs offline, and links it next to
  the "Open in Colab" badge; both links generated from the page path.
- Quarto renders with `execute-dir: project`, so the notebook must put its page
  folder on `sys.path` before `import model`; `gallery_utils` already locates
  the page folder (`_page_dir`) and can do this.

### 3.3 What stays

`agreement.json`, `check_agreement.py` with the pinned and latest pymrm matrix,
the published-work-only policy, provenance tiers, the break-row and
second-route rules. The restructuring must not change any agreement metric.

## 4. Migration

- **New pages** follow the structure from the next queue item on.
- **Existing pages** move when they are next touched for another reason, not in
  a mechanical sweep. Priority when choosing: pages with a model class that
  others will want to reuse (reactor models in D, H, G, I; nested and
  multi-domain structures S7, S8).
- **Acceptance per migrated page:** `check_agreement.py` passes with every metric
  unchanged within its stored tolerance; the notebook runs from a clean
  download (only `index.ipynb`) and from the zip; `model.py` imports without
  the notebook.
- **`build_page.py`:** decide per page. Where it only concatenates cells, retire
  it and edit the notebook directly (with `jupytext` pairing if a text source
  is wanted); where it computes content, keep it but have it import `model.py`.

## 5. Decisions for the maintainer

1. Tag scheme for published assets, and whether old pages are re-tagged when
   migrated.
2. Keep, pair (`jupytext`) or retire `build_page.py`.
3. Whether each `model.py` gets a tiny test file run in CI, or the notebook's
   own checks are enough.
4. Whether the zip artefact is worth the publish-workflow complexity, given that
   Colab and tagged single-notebook downloads already cover most users; the zip
   is for offline use.
