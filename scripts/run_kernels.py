#!/usr/bin/env python3
"""Run one page under three OpenBLAS kernels and compare what it publishes.

    python scripts/run_kernels.py pages/A2.4-tanks-in-series [--out DIR]

The page is first brought up to date from its source (``jupytext --to ipynb
index.py`` for a refreshed page, ``build_page.py`` for an old-layout one), then
executed with ``OPENBLAS_CORETYPE`` Sandybridge, Nehalem and Haswell, one thread,
``GALLERY_OFFLINE=1``, NumPy's SIMD dispatch pinned as in ``run_pages.py``. Haswell runs last, so the files left in the page folder
are the pinned baseline. Each run's notebook and ``agreement.json`` are kept in
``DIR/<page>/`` (default: a folder under the system temp directory).

Then two checks, which is what a refreshed page must pass
(docs/refresh-plan.md, lane A):

1. every claim of the Sandybridge and Nehalem runs passes its own check against
   the Haswell value (format 2 pages);
2. the printed text of the three runs is identical once numbers inside
   ``(here ...)`` are masked: a published sentence must not depend on the machine.
   The mask works line by line, so keep each ``(here ...)`` on one printed line.

Exit status 1 when a run fails or either check finds a difference.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "shared"))
from gallery_utils import _passes  # noqa: E402

KERNELS = ("Sandybridge", "Nehalem", "Haswell")     # the pinned one last
# NumPy's SIMD dispatch, pinned as in scripts/run_pages.py (an AVX-512 CPU would
# otherwise compute exp, log and pow differently in every run)
NPY_PIN = "X86_V4 AVX512_ICL AVX512_SPR"
HERE = re.compile(r"\(here [^)]*\)")


def build(page: Path) -> None:
    if (page / "index.py").is_file():
        subprocess.run([sys.executable, "-m", "jupytext", "--to", "ipynb", "index.py"],
                       cwd=page, check=True, capture_output=True)
    elif (page / "build_page.py").is_file():
        subprocess.run([sys.executable, "build_page.py"], cwd=page, check=True,
                       capture_output=True)


def execute(page: Path, kernel: str, out: Path) -> tuple[bool, float, str]:
    env = {"OPENBLAS_CORETYPE": kernel, "OPENBLAS_NUM_THREADS": "1",
           "OMP_NUM_THREADS": "1", "GALLERY_OFFLINE": "1",
           "NPY_DISABLE_CPU_FEATURES": NPY_PIN}
    saved = {k: os.environ.get(k) for k in env}
    os.environ.update(env)          # the kernel subprocess inherits this
    nb = nbformat.read(page / "index.ipynb", as_version=4)
    t0 = time.perf_counter()
    try:
        NotebookClient(nb, timeout=7200, kernel_name="python3",
                       resources={"metadata": {"path": str(page)}}).execute()
        ok, msg = True, ""
    except CellExecutionError as exc:
        ok, msg = False, str(exc).strip().splitlines()[-1][:300]
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    nbformat.write(nb, out / f"{kernel}.ipynb")
    if (page / "agreement.json").is_file():
        (out / f"{kernel}.json").write_text((page / "agreement.json").read_text())
    return ok, time.perf_counter() - t0, msg


def printed(path: Path) -> list[str]:
    nb = nbformat.read(path, as_version=4)
    lines = []
    for i, cell in enumerate(nb.cells):
        for o in cell.get("outputs", []):
            text = o.get("text") or o.get("data", {}).get("text/plain", "")
            text = "".join(text) if isinstance(text, list) else text
            lines += [f"cell {i}: {HERE.sub('(here #)', ln)}" for ln in text.splitlines()]
    return lines


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    page = Path(args[0]).resolve()
    base = (Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv
            else Path(tempfile.gettempdir()) / "pymrm-gallery-kernels")
    out = base / page.name
    out.mkdir(parents=True, exist_ok=True)

    build(page)
    status = 0
    for k in KERNELS:
        ok, dt, msg = execute(page, k, out)
        print(f"{'ok  ' if ok else 'FAIL'} {k:<12} {dt:7.1f} s  {msg}")
        status |= not ok
    nbformat.write(nbformat.read(out / "Haswell.ipynb", as_version=4), page / "index.ipynb")

    ref_path = out / "Haswell.json"
    if ref_path.is_file() and json.loads(ref_path.read_text()).get("format") == 2:
        ref = json.loads(ref_path.read_text())["metrics"]
        for k in KERNELS[:-1]:
            path = out / f"{k}.json"
            if not path.is_file():
                continue
            run = json.loads(path.read_text())["metrics"]
            bad = [n for n, c in ref.items()
                   if n not in run or not _passes(c, c["value"], run[n]["value"])]
            print(f"claims {k}: {len(ref) - len(bad)}/{len(ref)} pass against Haswell"
                  + (f"; failing: {', '.join(bad)}" if bad else ""))
            status |= bool(bad)

    ref_text = printed(out / "Haswell.ipynb")
    for k in KERNELS[:-1]:
        diff = sorted(set(printed(out / f"{k}.ipynb")) ^ set(ref_text))
        print(f"printed text {k}: {'identical' if not diff else f'{len(diff)} differing line(s)'}")
        for line in diff[:20]:
            print(f"    {line[:160]}")
        status |= bool(diff)
    print(f"runs kept in {out}")
    return status


if __name__ == "__main__":
    sys.exit(main())
