#!/usr/bin/env python3
"""Execute every gallery page notebook in place.

Used by CI. Each notebook runs with its own directory as the working directory
so that page-relative data paths behave the same as they do locally.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parents[1]
TIMEOUT = int(sys.argv[1]) if len(sys.argv) > 1 else 1800


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


def main() -> int:
    notebooks = sorted((ROOT / "pages").glob("*/index.ipynb"))
    if not notebooks:
        print("no page notebooks found")
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
