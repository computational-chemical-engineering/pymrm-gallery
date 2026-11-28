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
   last commit that changed the page notebook), ``catalog-id``, ``paper-year``
   (``reference.year`` of the page's record in models.yaml, the year the original
   paper appeared), a ``decade:<NNNN>s`` category, and ``paper-citation`` and
   ``paper-doi`` for the "Original paper" line of the title block
   (``_title-metadata.html``). Quarto merges a
   directory's ``_metadata.yml`` into its documents, categories included, so the
   title block shows both the published and the modified date, and the listings
   can show, sort and filter by them, without editing or executing a notebook.

The dates come from the history of ``--repo`` (default: the current directory).
When its HEAD is a release commit (``release gallery-...``, which rewrites the
bootstrap of every page), the dates are read from its parent, so a release
does not count as an update.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import yaml


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


def references(repo_root: Path) -> dict[str, dict]:
    """Map catalog ID to the ``reference`` of every record in models.yaml."""
    index = yaml.safe_load((repo_root / "models.yaml").read_text(encoding="utf-8"))
    return {m["id"]: m.get("reference") or {} for m in index["models"]}


def short_authors(authors: list[str]) -> str:
    """``Ergun``, ``Ranz and Marshall`` or ``Henze et al.`` from "Surname, I." entries."""
    names = [a.split(",")[0].strip() for a in authors]
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return f"{names[0]} and {names[1]}"
    return f"{names[0]} et al."


def own_categories(nb_path: Path) -> list[str]:
    """The ``categories`` of a notebook's front-matter cell."""
    cells = json.loads(nb_path.read_text(encoding="utf-8")).get("cells", [])
    src = "".join(cells[0].get("source", "")) if cells else ""
    if not src.startswith("---"):
        return []
    return list((yaml.safe_load(src.split("---")[1]) or {}).get("categories") or [])


def page_metadata(catalog_id: str, modified: str, ref: dict,
                  categories: list[str] = ()) -> dict:
    """The ``_metadata.yml`` content of one page; ``ref`` must carry a year.

    ``listing-order`` is the page's own categories, joined: the gallery's
    default order. Sorting on ``categories`` itself would put the merged
    ``decade:`` tag first and reorder the listing by decade.
    """
    year = ref["year"]
    meta = {
        "listing-order": ",".join(categories),
        "date-modified": modified,
        "catalog-id": catalog_id,
        "paper-year": str(year),
        "categories": [f"decade:{year // 10 * 10}s"],
        "paper-citation": f"{short_authors(ref.get('authors') or ['?'])} ({year})",
    }
    if ref.get("doi"):
        meta["paper-doi"] = str(ref["doi"])
    return meta


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
    refs = references(Path.cwd())
    missing = []
    no_year = []
    for page in pages:
        (page / "index.py").unlink(missing_ok=True)
        modified = _git(args.repo, "log", "-1", "--format=%cs", rev, "--",
                        f"{page.as_posix()}/index.ipynb")
        if not modified:
            missing.append(page.name)
            continue
        catalog_id = page.name.split("-", 1)[0]
        ref = refs.get(catalog_id, {})
        if not isinstance(ref.get("year"), int):
            no_year.append(page.name)
            continue
        (page / "_metadata.yml").write_text(yaml.safe_dump(
            page_metadata(catalog_id, modified, ref,
                          own_categories(page / "index.ipynb")), allow_unicode=True,
            sort_keys=False), encoding="utf-8")
    ok = len(pages) - len(missing) - len(no_year)
    print(f"prepared {ok} of {len(pages)} pages (history at {rev})")
    if missing:
        print("no history for: " + ", ".join(missing), file=sys.stderr)
    if no_year:
        print("no reference.year in models.yaml for: " + ", ".join(no_year),
              file=sys.stderr)
    return 1 if missing or no_year else 0


if __name__ == "__main__":
    sys.exit(main())
