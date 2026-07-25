# Page skeleton

Copy this whole directory to `pages/<catalog-id>-<slug>/` and replace the
physics. It is a working page (the Duncan-Toor model, `A4.9`) rather than an
empty shell, so it passes CI before you start editing and you can run it to see
what "done" looks like.

Replace, in order:

1. `meta.yaml`      - id, slug, title, section, tier, structures, authors
2. `data/`          - your CSV plus its `.meta.yaml` provenance sidecar
3. `index.ipynb`    - the nine required sections

Then:

    python scripts/check_metadata.py
    python scripts/run_pages.py

See `AGENTS.md` and `contributing.qmd` for the full rules.
