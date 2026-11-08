#!/usr/bin/env python3
"""Prepare a disposable copy of the gallery for ``quarto render``.

    python scripts/prepare_site.py --copy                  # in the build checkout
    python scripts/prepare_site.py --copy --repo ~/pymrm-gallery   # in a git-archive copy

Run by the ``build`` job of publish.yml before rendering. Two steps:

1. Remove every ``pages/*/index.py``. Quarto silently skips an ``index.ipynb``
   that has a jupytext-paired ``index.py`` beside it, and an exclusion in
   ``_quarto.yml`` does not help. The scripts are part of the page source, so
   this must never run in a working checkout: ``--copy`` is required.
2. Write ``pages/<page>/_metadata.yml`` with ``date-modified`` (the date of the
   last commit that changed the page notebook) and ``catalog-id``. Quarto merges
   a directory's ``_metadata.yml`` into its documents, so the title block shows
   both the published and the modified date, and the listings can show and sort
   by them, without editing a notebook.

The dates come from the history of ``--repo`` (default: the current directory).
When its HEAD is a release commit (``release gallery-...``, which rewrites the
bootstrap of every page), the dates are read from its parent, so a release
does not count as an update.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--copy", action="store_true",
                        help="confirm this directory is a disposable copy")
    parser.add_argument("--repo", type=Path, default=Path.cwd(),
                        help="git checkout to read the page history from")
    args = parser.parse_args()
    if not args.copy:
        parser.error("deletes pages/*/index.py: run only in a disposable copy, with --copy")

    rev = "HEAD"
    if _git(args.repo, "log", "-1", "--format=%s", "HEAD").startswith("release gallery-"):
        rev = "HEAD^"

    pages = sorted(p.parent for p in Path("pages").glob("*/index.ipynb"))
    missing = []
    for page in pages:
        (page / "index.py").unlink(missing_ok=True)
        modified = _git(args.repo, "log", "-1", "--format=%cs", rev, "--",
                        f"{page.as_posix()}/index.ipynb")
        if not modified:
            missing.append(page.name)
            continue
        catalog_id = page.name.split("-", 1)[0]
        (page / "_metadata.yml").write_text(
            f"date-modified: {modified}\ncatalog-id: \"{catalog_id}\"\n")
    print(f"prepared {len(pages) - len(missing)} of {len(pages)} pages (history at {rev})")
    if missing:
        print("no history for: " + ", ".join(missing), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
