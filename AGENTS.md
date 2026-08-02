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
- **At a jump in diffusivity, the face value is the HARMONIC mean.** Two cells
  with different D are two resistances in series, so the face conductivity is
  `2 D_L D_R / (D_L + D_R)`, not `(D_L + D_R)/2`. The difference is an **order**,
  not a factor: at the jump the arithmetic mean converges at *first* order and the
  harmonic at *second*, so the error ratio grows ∝ n without bound. Measured on
  `A2.1`'s three-section vessel, refining only the bed: **84× at n = 100, 172× at
  200, 347× at 400, 696× at 800.** Quote the order, never a single factor — on a
  coarse grid the gap looks small enough to dismiss the rule. And it **fails
  silently**: the wrong mean still gives a smooth, plausible profile (f(0+)
  0.667011 against 0.667012, exit unchanged to 4e-9). Only matters where D
  actually jumps; harmless where it is smooth.
- **`NumJac` stencil**: `NumJac(shape)` couples the **last axis in full** —
  correct for a pointwise `reaction(c)` when the last axis is the field index.
  `axes_blocks=[-2,-1]` for phase-and-species coupling. Two traps, both measured
  on published pages (2026-08-01):

  **Never pass a bare 1-D shape.** For a single field write `(n, 1)`, not `(n,)`.
  With `(n,)` the last axis *is space*, so the default stencil declares every cell
  coupled to every other and builds a **dense n×n Jacobian** — 3.2 s to construct
  at n=400 and 70 s at n=1600, against 0.5 ms for the right shape. The answers are
  bit-identical; only the cost changes. This was live on `B1.1`, `B1.6` and
  `F3.1`, and cost `B1.1` 6.3× its runtime.

  **`axes_diagonals=[0]` is only meaningful once `ndims ≥ 2`.** An earlier version
  of this file said to add it "when the source term depends on neighbouring
  cells", with no dimension caveat. On a 1-D shape that is not merely over-broad,
  it is *wrong*: `axes_blocks` still defaults to `[-1]`, which is the same axis 0,
  so `stencil_block_diagonals` treats axis 0 as a fixed axis and reinterprets the
  `[-1, 0, 1]` offsets as **absolute indices** n−1, 0, 1. The Jacobian comes out
  with no diagonal at all and the solve converges to a different answer. Write
  `NumJac((n, 1), axes_diagonals=[0])` instead — and only when the *source term
  itself* reads neighbouring cells, which is rarer than it sounds: the coupling
  from the Laplacian normally arrives analytically through the divergence
  operator, not through the function handed to `NumJac`.

## When the model is read from a paper that reprints it

Many classics are unreachable — pre-DOI, no open-access route — while a paper
that **is** on disk prints the model in full with attribution. Building from that
reprint is legitimate and is how `F1.3` (Wilkinson, via Krishna & Ellenberger
1996) and `B1.6` (Prater, via Weisz & Hicks 1962) were built. Building from
memory or from a textbook you have not opened is not: that is the "textbook
restatement passed off as the source's" the builder brief forbids, and it is the
fabrication route for this class of case.

Record both, and keep them distinct:

```yaml
reference:                    # the ORIGIN of the result — cited, not consulted
  authors: ["Wilkinson, P. M.", ...]
  year: 1992
reference_read_from:          # the text actually read, and where each equation came from
  authors: ["Krishna, R.", "Ellenberger, J."]
  year: 1996
  note: "Reprints the correlation as its Eqs. 1-4; verified on a 600 dpi render."
```

Say the same thing in prose on the page. A reader must never have to guess which
document a transcription came from, and the next agent must not go looking for a
paper nobody has.

**Check the reprint actually carries the case before using it.** `E1.1` was
parked, not built, because Kunii & Levenspiel state the two-phase relation but
never attribute it to Toomey & Johnstone, never name it and never test it — so
the reprint could not source a page *about* that postulate. Ask whether the
reprint supports the claim the case is for, not merely whether the equation
appears in it.

## Published work only

The gallery reproduces **published** models. If a model's source is a
manuscript in preparation, under review, or in revision, it does not get a page
yet — regardless of how good the data situation is, and regardless of whether
the authors are the gallery maintainers.

Such entries carry `status: deferred` plus a `blocked_by` field saying what
lifts the block. Do not build a page for a `deferred` entry, and do not promote
one to `planned` without being told to. `scripts/check_metadata.py` enforces
that every deferred entry explains itself and has no page directory.

## Data rules — the ones that matter most

- **Never fabricate data.** If no dataset exists for a model, set
  `data: {status: placeholder}` and mark the page `status: planned`. A page with
  invented numbers is worse than no page.
- Every `data/*.csv` requires a `data/*.meta.yaml` provenance sidecar. The
  schema and the legal basis are in [`docs/data-strategy.md`](docs/data-strategy.md).
- **Loading another page's dataset means reading that page.** A CSV borrowed
  through `load_data(..., page=...)` is not a bare table of numbers — it belongs
  to a page that has already established findings about those very rows, and
  those findings do not travel with the file. Two obligations:

  1. In your own *The data* section, list every finding the source page states
     about the rows you use — flagged rows, derived quantities it already
     computes, known bad cells, discrepancies against the paper — and say whether
     each affects you.
  2. **If a number you are about to state also exists in a dataset you loaded,
     print it beside yours and reconcile the two. Never retype a value that is a
     row in a CSV you already read.**

  Rule 2 is the cheap one and it catches most of this class. Measured on the
  2026-08-02 sweep of the eight pages that borrow data, it alone would have
  caught five of the nine findings, plus the defect that prompted the sweep:
  `A1.6` built its central argument on an inverted voidage of 0.399 while
  `A1.7` — the page supplying the CSV — already computed 0.444 from Geldart's own
  columns and already flagged the two rows that give it with no inference at all.
  At that voidage the reference balance `A1.6` held up as the standard came out
  **+58 % biased, worse than the correlation it was judging**, and the error
  reached the Reuse advice. `J3.1` did the same thing in the other direction: it
  loaded a `V_cutoff` and a stated result, used neither, and reported an
  overpotential taken past the cutoff — 3.7 mV where its own data give 1.7.

  Reuse advice is where this class does its damage, so check it there explicitly.
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
