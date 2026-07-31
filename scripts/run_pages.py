#!/usr/bin/env python3
"""Execute gallery page notebooks in place.

    python scripts/run_pages.py             # all pages
    python scripts/run_pages.py --changed   # only pages touched since origin/main
    python scripts/run_pages.py 900         # all pages, 900 s timeout

Each notebook runs with its own directory as the working directory so that
page-relative data paths behave the same as they do locally.
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parents[1]
ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
TIMEOUT = int(ARGS[0]) if ARGS else 1800


def run(path: Path) -> tuple[bool, float, str]:
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(
        nb, timeout=TIMEOUT, kernel_name="python3",
        resources={"metadata": {"path": str(path.parent)}},
    )
    t0 = time.perf_counter()
    try:
        client.execute()
    except CellExecutionError as exc:
        return False, time.perf_counter() - t0, str(exc).splitlines()[-1][:200]
    nbformat.write(nb, path)
    return True, time.perf_counter() - t0, ""


def changed_pages() -> set[str]:
    """Page directories touched since origin/main, from git."""
    try:
        out = subprocess.run(["git", "diff", "--name-only", "origin/main...HEAD"],
                             cwd=ROOT, capture_output=True, text=True, check=True).stdout
        out += subprocess.run(["git", "status", "--porcelain"],
                              cwd=ROOT, capture_output=True, text=True, check=True).stdout
    except Exception:
        return set()
    hits = set()
    for line in out.splitlines():
        parts = line.split()
        path = parts[-1] if parts else ""
        bits = Path(path).parts
        if len(bits) > 1 and bits[0] == "pages":
            hits.add(bits[1])
    return hits


def main() -> int:
    notebooks = sorted((ROOT / "pages").glob("*/index.ipynb"))
    if not notebooks:
        print("no page notebooks found")
        return 0

    # Executing every page on every merge costs ~13 min and grows with the
    # gallery, while a merge usually touches one page. Run the changed pages by
    # default; the nightly workflow passes --all for the full sweep that catches
    # cross-page breakage (a shared helper, a pymrm release).
    if "--changed" in sys.argv:
        only = changed_pages()
        if only:
            notebooks = [n for n in notebooks if n.parent.name in only]
            print(f"--changed: {len(notebooks)} of {len(list((ROOT/'pages').glob('*/index.ipynb')))} "
                  f"page(s) touched since origin/main\n")
        else:
            print("--changed: nothing under pages/ has changed; running nothing\n")
            return 0

    failures = []
    for nb_path in notebooks:
        page = nb_path.parent.name
        ok, dt, msg = run(nb_path)
        status = "OK  " if ok else "FAIL"
        print(f"{status} {page:<40} {dt:6.1f} s" + (f"  {msg}" if msg else ""))
        if not ok:
            failures.append(page)

    print(f"\n{len(notebooks) - len(failures)}/{len(notebooks)} pages executed")
    if failures:
        print("failed: " + ", ".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
