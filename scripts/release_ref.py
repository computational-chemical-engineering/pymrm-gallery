#!/usr/bin/env python3
"""Rewrite GALLERY_REF = "main" to a release tag in every migrated page.

    python scripts/release_ref.py gallery-2026.10.30

Used by the ``release`` job of publish.yml on a detached HEAD, before the
release commit is tagged (docs/page-infrastructure.md, section 1). Edits the
bootstrap cell of each migrated ``index.ipynb`` through nbformat (a text match
on the JSON would fail on the escaped quotes) and ``index.py`` as text. Each file
must contain exactly one ``GALLERY_REF = "main"`` line; otherwise that file is
not written and the script fails, so the job stops before the release commit.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from page_layout import ROOT, migrated_pages

PATTERN = re.compile(r'^(GALLERY_REF = )"main"', re.MULTILINE)


def _replace(text: str, tag: str, where: str) -> str:
    new, n = PATTERN.subn(lambda m: f'{m[1]}"{tag}"', text)
    if n != 1:
        raise ValueError(f'{where}: expected exactly one GALLERY_REF = "main", found {n}')
    return new


def rewrite_notebook(path: Path, tag: str) -> None:
    """Set GALLERY_REF to ``tag`` in the notebook's code cells, exactly once."""
    import nbformat

    nb = nbformat.read(path, as_version=4)
    hits = [c for c in nb.cells if c.cell_type == "code" and PATTERN.search(c.source)]
    total = sum(len(PATTERN.findall(c.source)) for c in hits)
    if total != 1:
        raise ValueError(f'{path}: expected exactly one GALLERY_REF = "main", found {total}')
    hits[0].source = _replace(hits[0].source, tag, str(path))
    nbformat.write(nb, path)


def rewrite_py(path: Path, tag: str) -> None:
    """Set GALLERY_REF to ``tag`` in a percent-format ``index.py``, exactly once."""
    path = Path(path)
    path.write_text(_replace(path.read_text(encoding="utf-8"), tag, str(path)),
                    encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) != 1 or not re.fullmatch(r"[A-Za-z0-9._-]+", argv[0]):
        print(__doc__)
        return 2
    tag = argv[0]
    pages = migrated_pages(ROOT)
    for page in pages:
        rewrite_notebook(page / "index.ipynb", tag)
        if (page / "index.py").is_file():
            rewrite_py(page / "index.py", tag)
        print(f"{page.name}: GALLERY_REF = {tag!r}")
    print(f"{len(pages)} migrated page(s) set to {tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
