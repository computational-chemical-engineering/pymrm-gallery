# Instructions for coding agents

This repository is designed to be used by coding agents assisting researchers.
Read this file before doing anything else in the repo.

## Start here, not with the website

Read [`models.yaml`](models.yaml). It is the canonical machine-readable index of
every model in the gallery, built and planned. Do **not** scrape the rendered
HTML — it is generated from this file and lags behind it.

`models.json` is emitted at build time if you prefer JSON.

## The catalog ID is the stable key

`A4.9`, `D1.4`, `J3.4` — these never change. Slugs, titles, and directory names
may. Reference models by catalog ID in anything you write.

## To build a model for a user

1. Identify the **mathematical structure** of the user's problem, not just its
   physical domain. The `structures` field uses codes `S1`–`S13`, defined in
   [`docs/taxonomy.md`](docs/taxonomy.md).
2. Search `models.yaml` for entries sharing those structure codes, preferring
   `status: published`.
3. **Copy that page's directory and substitute the physics.** Do not start from
   a blank notebook and assemble operators from the API documentation — that is
   how the common pymrm mistakes below get made.
4. If nothing shares the structure, say so, and pick the closest `status:
   planned` entry as a specification to implement.

Structure matching crosses fields deliberately. A diffusion-limited enzyme
pellet (`J4.7`), a catalyst pellet (`B1.1`), and the solid phase of a lithium-ion
electrode (`J3.4`) are the same problem. If the user's field has no page, another
field's probably does.

## pymrm conventions you must follow

The authority is `pymrm/docs/pymrm-model-style-guide.md` in the sibling `pymrm`
repository. The errors agents actually make:

- **Boundary conditions use the OUTWARD normal.** `bc` is a 2-tuple of dicts
  `{"a": ..., "b": ..., "d": ...}` meaning `a·∂c/∂n + b·c = d` with `n` pointing
  *out of* the domain. The sign of `a` therefore flips meaning between the two
  ends: `{a:1, b:0, d:q}` at the left boundary gives `dc/dx = −q`, and the same
  dict at the right gives `dc/dx = +q`. Write the physical equation in a comment
  next to every `bc`.
- **Assemble constant operators once**, in `__init__` or `_build_operators` —
  never inside a Newton iteration or time step. If a boundary value changes
  every step, that is what `shapes_d` is for: it makes the boundary value an
  external unknown multiplying a *constant* matrix.
- **Spatial axes first, fields last**, and one layout for the whole model:
  `(n_x, n_c)`, `(n_z, n_phase, n_c)`, `(n_z, n_r, n_c)`. Flatten only at the
  residual/Jacobian interface.
- **Use the operators.** `construct_grad`, `construct_div`,
  `construct_convflux_upwind`, `construct_coefficient_matrix`. Do not hand-build
  finite-difference stencils with `scipy.sparse.diags`.
- **`nu` in `construct_div`** is geometry: `0` Cartesian, `1` cylindrical, `2`
  spherical, or a callable for an arbitrary area profile. State it in a comment.
- **`NumJac` stencil**: `NumJac(shape)` couples only the last axis (correct for
  a pointwise `reaction(c)`); add `axes_diagonals=[0]` when the source term
  depends on neighbouring cells; `axes_blocks=[-2,-1]` for phase-and-species
  coupling.

## Data rules — the ones that matter most

- **Never fabricate data.** If no dataset exists for a model, set
  `data: {status: placeholder}` and mark the page `status: planned`. A page with
  invented numbers is worse than no page.
- Every `data/*.csv` requires a `data/*.meta.yaml` provenance sidecar. The
  schema and the legal basis are in [`docs/data-strategy.md`](docs/data-strategy.md).
- Never commit source PDFs. DOIs and extracted CSVs only.
- Never reproduce a source figure image. Extract the points and re-plot.

## Adding a page

```
pages/<catalog-id>-<slug>/
├── index.ipynb          # from templates/page-template.ipynb
├── meta.yaml            # must not contradict models.yaml
├── data/
│   ├── <name>.csv
│   └── <name>.meta.yaml
└── README.md
```

Required in the notebook, in this order: Background · The published model ·
Parameters and assumptions · The data · PyMRM implementation · Results ·
Validation · What pymrm adds · Reuse.

Constraints:

- Cell 1 is the Colab environment cell (`try: import pymrm / except ImportError:
  %pip install -q pymrm`).
- Load data with `gallery_utils.load_data(name, page=...)`, never a bare
  relative path — it must work on a fresh Colab VM with no checkout.
- **No Quarto-specific Markdown inside `index.ipynb`.** No `::: {.callout-*}`,
  no shortcodes. The notebook must render correctly in Jupyter, GitHub preview,
  and VS Code. Quarto syntax belongs in `.qmd` pages and YAML front matter only.
- Include at least one validation: analytical limit, conservation check, grid or
  time-step independence, or a physical bound.
- Call `report_agreement(...)` so CI can detect regressions.

Then run `python scripts/check_metadata.py` before opening a PR.

## Honesty requirements

The "What pymrm adds" section must be truthful. If the page merely reproduces
the original with no extension, write that. Do not manufacture an improvement,
and do not describe agreement as good when the metric says otherwise — the
agreement number is printed on the page and checked by CI.
