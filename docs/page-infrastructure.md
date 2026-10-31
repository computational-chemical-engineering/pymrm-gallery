# Design: infrastructure for refreshed pages (refresh plan, step 1)

Status: REVISED 2026-10-29 after review (`think`, "adopt with changes"; changes folded
in). Implements decisions 1 to 3 of
[`page-structure-proposal.md`](page-structure-proposal.md) section 5 and step 1 of
[`refresh-plan.md`](refresh-plan.md). Nothing here changes an unmigrated page.

Built 2026-10-30 on branch `refresh/step1-infrastructure`; the sections below
describe the implementation, with its deviations from the reviewed design stated
where they occur (sections 1 to 7).

## 1. Release tags

Decision 1 says bootstraps never reference `main` and each published notebook fetches
the files it was published with. A notebook cannot contain the name of a tag that points
at itself, so the release commit lives only on the tag. Tags pushed with `GITHUB_TOKEN`
trigger no workflow, so there is no loop.

`publish.yml` gets a `release` job before `build`, only on `refs/heads/main`, with
job-level `permissions: contents: write`:

1. Full checkout with tags; a git identity for the bot.
2. If a `gallery-*` tag already belongs to `HEAD` (a re-run), reuse it. The tag sits
   on the release commit, not on `HEAD`, so "belongs" means the tag's commit is
   `HEAD` or its parent is `HEAD`. Otherwise the tag
   is `gallery-YYYY.MM.DD`, with `.2`, `.3` for further deploys that day, found with
   `git ls-remote --tags`.
3. On a detached `HEAD`, rewrite `GALLERY_REF = "main"` to the tag in every migrated
   page's `index.ipynb` (through `nbformat`, editing the bootstrap cell's source; a text
   match fails on the escaped quotes) and `index.py`, asserting exactly one replacement
   per file (`scripts/release_ref.py`). Commit "release gallery-YYYY.MM.DD"
   (`--allow-empty`, so there is a release commit even with no migrated page), tag
   it, run `check_metadata.py --gallery-ref <tag>` on that tree, and only then push
   the tag (nothing else). A broken tag therefore never leaves the job; a pushed one is reused
   by every re-run.
4. `build` checks out **the tag** and renders from it. Since `build` needs `release`,
   which runs on `main` only, a `workflow_dispatch` from another branch no longer
   deploys. Rendering the tag makes the site show the tag in every bootstrap. Quarto does not execute the notebooks (`_quarto.yml` has no
   `execute: enabled`) and `_freeze/` is not committed, so rendering the tag tree is
   safe. The stale comment in `publish.yml` about a committed `_freeze/` is corrected.
5. The "Open in Colab" badge and download link are generated from the page path and the
   tag through a Quarto filter (`_extensions/` or a `_metadata.yml` the job writes),
   never through Quarto syntax in a notebook. They exist only once the tag is pushed,
   because `build` needs `release`. Not built yet (a TODO in `publish.yml`).

An organisation tag ruleset could refuse the bot's push; check before the first deploy.

## 2. The bootstrap

The first code cell of every migrated page, identical apart from its literals. It
replaces the old Colab environment cell (`try: import pymrm / except ImportError:
%pip install`), so it is the first code cell, not the second. The canonical text,
with placeholders, is `templates/bootstrap_cell.py`; as built:

```python
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "main"  # replaced by the release tag at publication
PYMRM = ">=3.3.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A2.4-tanks-in-series"  # equals the page directory

import os, sys, urllib.request
from pathlib import Path
_shared = next((p / "shared" for p in (Path.cwd(), *Path.cwd().parents)
                if (p / "shared" / "gallery_bootstrap.py").is_file()), None)
if _shared is None:
    _shared = (Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
               / "pymrm-gallery" / GALLERY_REF)
    _shared.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/computational-chemical-engineering/"
        f"pymrm-gallery/{GALLERY_REF}/shared/gallery_bootstrap.py",
        _shared / "gallery_bootstrap.py")
sys.path.insert(0, str(_shared))
sys.modules.pop("gallery_bootstrap", None)
from gallery_bootstrap import setup
setup(PAGE, PYMRM, GALLERY_REF);  # the semicolon stops the cell echoing the path
```

No `%pip` in the cell: `jupytext` comments magics out in `index.py`, which would leave an
`if` with no body. `shared/gallery_bootstrap.py` uses the **standard library only** (a
test enforces it), because it runs before numpy or pymrm may exist. `setup`:

1. **pymrm.** Missing: install with `subprocess.check_call([sys.executable, "-m",
   "pip", "install", "-q", f"pymrm{PYMRM}"])`. Present and inside the specifier
   (checked with `prereleases=True`, so CI's `latest` job's release candidates count):
   nothing. Present and newer: warn, never reinstall (CI's `latest` job must keep its
   version). Present and below the floor: in a checkout, raise (never pip-install
   into a development environment); otherwise reinstall if pymrm has not been
   imported in this kernel, and raise "restart the runtime" if it has. The version
   comparison is the module's own (PEP 440 order, pre-releases admitted as with
   `prereleases=True`), since `packaging` may be absent.
2. **Ref.** Sets `os.environ["GALLERY_REF"]`; `gallery_utils` reads it when it builds a
   URL, not at import.
3. **Files.** A checkout puts `shared/` and the page folder on `sys.path`. Without one,
   `gallery_utils.py` and `pages/<PAGE>/model.py` (a 404 is fine: not every page has
   one) are downloaded to `$XDG_CACHE_HOME/pymrm-gallery/<ref>/<PAGE>/`, falling back
   to `~/.cache/pymrm-gallery/<ref>/<PAGE>/`, never to the working directory (a
   user's own `model.py` must survive) and never to a predictable path in a shared
   `/tmp` (it goes first on `sys.path`). The cell downloads `gallery_bootstrap.py` to
   `.../pymrm-gallery/<ref>/` the same way. A `model` module already in
   `sys.modules` from another page is dropped.
4. **Offline mode.** With `GALLERY_OFFLINE=1` (set by `run_pages.py` and the model
   tests), nothing is fetched: `_resolve`, `load_claims` and `setup` raise when a file
   is missing locally. Otherwise a file deleted in a PR would pass CI by loading
   `main`'s copy.

## 3. The pymrm pin

`meta.yaml` `pymrm_min_version` is the single source. The bootstrap literal is
`">={floor},<{major+1}"`. Migration sets `pymrm_min_version` to the lock's pymrm (all
unmigrated pages say `"2.3"`, which would pin `>=2.3,<3` and install 2.5 on Colab).
A page is migrated when it has an `index.py` or a code cell imports
`gallery_bootstrap` or assigns `GALLERY_REF`. `check_metadata.py` checks, for every
migrated page: the bootstrap is the first code cell and has the exact
`GALLERY_REF = ` line; the pin matches `meta.yaml`; the `pymrm==` of
`requirements-lock.txt` satisfies the pin (through `gallery_bootstrap`'s comparison);
`PAGE` equals the directory name; `GALLERY_REF` equals `--gallery-ref` (default
`main`; the release job passes the tag); the cell equals the template apart from its
three literals; `index.py` exists and `build_page.py` does not. CI's `locked` job tests the
lock, which must satisfy every page's pin.

## 4. `jupytext`

`index.py` (percent format) is the edited source of a migrated page and replaces
`build_page.py`. `jupytext` is in `requirements.txt` as well as the lock, because the
publish build runs `check_metadata.py` from `requirements.txt`. Pairing is per notebook (`jupytext.formats: "ipynb,py:percent"` in its
metadata). `jupytext` joins `requirements-lock.txt`. A CI step fails when an `index.py`
and its `index.ipynb` differ in cell sources (outputs ignored). Content a `build_page.py`
computed at build time moves into code cells or `model.py`, never into a new generator.

## 5. `model.py` and its test

Where a page defines a model class: `pages/<id>/model.py` (classes and functions,
NumPy-style docstrings, no plotting, no data loading) and `pages/<id>/test_model.py`
(import, one fast solve, one stored value read with `load_claims` and checked by that
claim's own rule; seconds). CI runs each `test_model.py` in its **own** pytest process
(same file name in folders without `__init__.py` collides, and `import model` would reuse
the first page's module), with `GALLERY_OFFLINE=1` and `PYTHONPATH=shared`, so a test
imports `gallery_utils` without running the bootstrap.

## 6. The pymrm plugin's checkers

`check_model.py` and `check_notebook_math.py` come from the pymrm repo at the tag that
matches the lock (`v3.3.0`), never a branch. On migrated pages: `check_notebook_math.py`
fails the job, with KaTeX installed in CI (`npm install --prefix ~/.cache/pymrm-katex
katex@0.16`) so its parse step runs; `check_model.py` warns (`continue-on-error`) until
the pilot measures its false-positive rate. `check_model.py` runs a script, and
`model.py` only defines classes, so it runs `index.py` in the page directory. That
rewrites the tracked `agreement.json` without the BLAS pin, so the checker steps come
after every step that reads those files, and the last one ends with
`git checkout -- pages`.

## 7. Smoke test of the download path

After `release`, a job takes one migrated page and runs its `index.ipynb` from the tag in
an empty directory and a fresh venv without pymrm. Nothing else tests the fetch path,
and "runs from a clean download" is an acceptance criterion of lane A.

As built (`smoke` in `publish.yml`): the first page listed by `check_metadata.py
--list-migrated`; the job succeeds without running anything while there is none. The
venv gets `nbclient`, `ipykernel` and the Colab defaults the template assumes
(`pandas`, `matplotlib`, `pyyaml`), never pymrm, and is first on `PATH` so the kernel
runs in it. The notebook comes from the tag's raw URL and is executed with nbclient;
an error fails the job, and `deploy` needs both `build` and `smoke`.

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

pymrm 3.1.0, 3.2.0 and 3.3.0 declare `ipython`, `ipykernel`, `nbclient`, `pytest`, `pytest-cov`, `sphinx`,
`sphinx_rtd_theme`, `myst-parser` and `pip` as runtime dependencies (its
`requirements.txt` feeds `pyproject.toml`). Installing it on Colab can upgrade IPython
or ipykernel under the live kernel. The fix belongs in pymrm (move them to extras); until
it is released, the smoke job (section 7) and a manual Colab run of the first pilot page
are the check.
