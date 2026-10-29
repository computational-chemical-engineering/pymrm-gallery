# Design: infrastructure for refreshed pages (refresh plan, step 1)

Status: REVISED 2026-10-29 after review (`think`, "adopt with changes"; changes folded
in). Implements decisions 1 to 3 of
[`page-structure-proposal.md`](page-structure-proposal.md) section 5 and step 1 of
[`refresh-plan.md`](refresh-plan.md). Nothing here changes an unmigrated page.

## 1. Release tags

Decision 1 says bootstraps never reference `main` and each published notebook fetches
the files it was published with. A notebook cannot contain the name of a tag that points
at itself, so the release commit lives only on the tag. Tags pushed with `GITHUB_TOKEN`
trigger no workflow, so there is no loop.

`publish.yml` gets a `release` job before `build`, only on `refs/heads/main`, with
job-level `permissions: contents: write`:

1. Full checkout with tags; a git identity for the bot.
2. If `HEAD` already carries a `gallery-*` tag (a re-run), reuse it. Otherwise the tag
   is `gallery-YYYY.MM.DD`, with `.2`, `.3` for further deploys that day, found with
   `git ls-remote --tags`.
3. On a detached `HEAD`, rewrite `GALLERY_REF = "main"` to the tag in every migrated
   page's `index.ipynb` (through `nbformat`, editing the bootstrap cell's source; a text
   match fails on the escaped quotes) and `index.py`, asserting exactly one replacement
   per file. Commit "release gallery-YYYY.MM.DD"; push the tag only.
4. `build` checks out **the tag** and renders from it, so the site shows the tag in every
   bootstrap. Quarto does not execute the notebooks (`_quarto.yml` has no
   `execute: enabled`) and `_freeze/` is not committed, so rendering the tag tree is
   safe. The stale comment in `publish.yml` about a committed `_freeze/` is corrected.
5. The "Open in Colab" badge and download link are generated from the page path and the
   tag through a Quarto filter (`_extensions/` or a `_metadata.yml` the job writes),
   never through Quarto syntax in a notebook. They exist only once the tag is pushed,
   because `build` needs `release`.

An organisation tag ruleset could refuse the bot's push; check before the first deploy.

## 2. The bootstrap

The first code cell of every migrated page, identical apart from its literals:

```python
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "main"           # replaced by the release tag at publication
PYMRM = ">=3.1.0,<4"           # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A2.4-tanks-in-series"  # equals the page directory

import sys, urllib.request
from pathlib import Path
_shared = next((p / "shared" for p in (Path.cwd(), *Path.cwd().parents)
                if (p / "shared" / "gallery_bootstrap.py").is_file()), None)
if _shared is None:
    _shared = Path(urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/computational-chemical-engineering/"
        f"pymrm-gallery/{GALLERY_REF}/shared/gallery_bootstrap.py",
        Path.home() / f"gallery_bootstrap_{GALLERY_REF}.py")[0]).parent
sys.path.insert(0, str(_shared))
from gallery_bootstrap import setup
setup(PAGE, PYMRM, GALLERY_REF)
```

No `%pip` in the cell: `jupytext` comments magics out in `index.py`, which would leave an
`if` with no body. `shared/gallery_bootstrap.py` uses the **standard library only** (a
test enforces it), because it runs before numpy or pymrm may exist. `setup`:

1. **pymrm.** Missing: install with `subprocess.check_call([sys.executable, "-m",
   "pip", "install", "-q", f"pymrm{PYMRM}"])`. Present and inside the specifier
   (checked with `prereleases=True`, so CI's `latest` job's release candidates count):
   nothing. Present and newer: warn, never reinstall (CI's `latest` job must keep its
   version). Present and below the floor: reinstall if pymrm has not been imported in
   this kernel, otherwise raise "restart the runtime".
2. **Ref.** Sets `os.environ["GALLERY_REF"]`; `gallery_utils` reads it when it builds a
   URL, not at import.
3. **Files.** A checkout puts `shared/` and the page folder on `sys.path`. Without one,
   `gallery_utils.py` and `pages/<PAGE>/model.py` are downloaded to
   `tempfile.gettempdir()/pymrm-gallery/<ref>/<PAGE>/`, never to the working directory
   (a user's own `model.py` must survive), that folder goes first on `sys.path`, and a
   `model` module already in `sys.modules` from another page is dropped.
4. **Offline mode.** With `GALLERY_OFFLINE=1` (set by `run_pages.py` and the model
   tests), nothing is fetched: `_resolve`, `load_claims` and `setup` raise when a file
   is missing locally. Otherwise a file deleted in a PR would pass CI by loading
   `main`'s copy.

## 3. The pymrm pin

`meta.yaml` `pymrm_min_version` is the single source. The bootstrap literal is
`">={floor},<{major+1}"`. `check_metadata.py` checks, for every migrated page: the pin
matches `meta.yaml`; `PAGE` equals the directory name; `GALLERY_REF == "main"` on `main`;
the cell equals the template apart from its three literals. CI's `locked` job tests the
lock, which must satisfy every page's pin.

## 4. `jupytext`

`index.py` (percent format) is the edited source of a migrated page and replaces
`build_page.py`. Pairing is per notebook (`jupytext.formats: "ipynb,py:percent"` in its
metadata). `jupytext` joins `requirements-lock.txt`. A CI step fails when an `index.py`
and its `index.ipynb` differ in cell sources (outputs ignored). Content a `build_page.py`
computed at build time moves into code cells or `model.py`, never into a new generator.

## 5. `model.py` and its test

Where a page defines a model class: `pages/<id>/model.py` (classes and functions,
NumPy-style docstrings, no plotting, no data loading) and `pages/<id>/test_model.py`
(import, one fast solve, one stored value read with `load_claims` and checked by that
claim's own rule; seconds). CI runs each `test_model.py` in its **own** pytest process
(same file name in folders without `__init__.py` collides, and `import model` would reuse
the first page's module), with `GALLERY_OFFLINE=1`.

## 6. The pymrm plugin's checkers

`check_model.py` and `check_notebook_math.py` come from the pymrm repo at the tag that
matches the lock (`v3.1.0`), never a branch. On migrated pages: `check_notebook_math.py`
fails the job; `check_model.py` warns (`continue-on-error`) until the pilot measures its
false-positive rate.

## 7. Smoke test of the download path

After `release`, a job takes one migrated page and runs its `index.ipynb` from the tag in
an empty directory and a fresh venv without pymrm. Nothing else tests the fetch path,
and "runs from a clean download" is an acceptance criterion of lane A.

## 8. Order of implementation (one PR)

1. `shared/gallery_bootstrap.py` (stdlib only), `GALLERY_REF` and `GALLERY_OFFLINE` in
   `gallery_utils`, unit tests (including a fresh-venv run of `setup` and the
   stdlib-only check).
2. `check_metadata.py` checks (section 3) and the `jupytext` sync check.
3. CI: per-page `test_model.py`, the two checkers, offline mode in `run_pages.py`.
4. `publish.yml`: `release` job, render from the tag, badge filter, smoke job.
5. `templates/page-template.ipynb` and `templates/index.py`; `AGENTS.md` "Adding a page".

The pilot pages are its first users and come after this PR.

## 9. Known risk outside this repo

pymrm 3.1.0 declares `ipython`, `ipykernel`, `nbclient`, `pytest`, `pytest-cov`, `sphinx`,
`sphinx_rtd_theme`, `myst-parser` and `pip` as runtime dependencies (its
`requirements.txt` feeds `pyproject.toml`). Installing it on Colab can upgrade IPython
or ipykernel under the live kernel. The fix belongs in pymrm (move them to extras); until
it is released, the smoke job (section 7) and a manual Colab run of the first pilot page
are the check.
