# Staged data

Datasets that are extracted, reviewed and provenance-documented, but whose page
is not built yet.

They live here rather than under `pages/<id>/data/` because a directory under
`pages/` declares a page, and `scripts/check_metadata.py` rightly errors on one
with no `meta.yaml`. Staging them keeps the extraction and its review in git
without pretending the page exists.

Move the whole directory to `pages/<id>-<slug>/data/` when the page is built.

| Directory | Catalog | State |
|---|---|---|
| `F1.4/` | `F1.4` Krishna & Ellenberger | Figure 11 digitised, 63 positions, reviewed by the maintainer 2026-07-27. Gas labels deliberately incomplete — see the sidecar's `review` block. |
