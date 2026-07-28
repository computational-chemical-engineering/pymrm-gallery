# Staged data

Datasets that are extracted, reviewed and provenance-documented, but whose page
is not built yet.

They live here rather than under `pages/<id>/data/` because a directory under
`pages/` declares a page, and `scripts/check_metadata.py` rightly errors on one
with no `meta.yaml`. Staging them keeps the extraction and its review in git
without pretending the page exists.

Move the whole directory to `pages/<id>-<slug>/data/` when the page is built.

All four below were reviewed by the maintainer on 2026-07-29 against the
original figures; each sidecar's `review:` block records the verdict verbatim and
what changed because of it.

| Directory | Catalog | State |
|---|---|---|
| `H1.7/` | Wijmans & Baker Fig. 5 | 13 markers, **all correct on review**. Self-checks: flux linear to r²=0.9998, intercept 322 psi vs the stated ~350 psi osmotic pressure. |
| `F2.3/` | Maretto & Krishna Fig. 2 | 95 resolved + 20 flagged unresolved. Review corrected the series assignment; identity now comes from curve proximity, not matched shape. |
| `J3.4/` | Doyle–Fuller–Newman Fig. 2 | 6 traced curves. **Contains no experimental data** — all simulation, confirmed on review. Reference-solution page only. |
| `G1.8/` | Herskowitz & Smith Fig. 6 | 4 fitted lines. Replaces the abandoned Fig. 2 attempt on the maintainer's decision; Fig. 2 is cited, not digitised. |

Earlier: `F1.4` moved to `pages/F1.4-krishna-ellenberger-holdup/data/` when that
page was built, 2026-07-28.
