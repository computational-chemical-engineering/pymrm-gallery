#!/usr/bin/env python3
"""Execute gallery page notebooks in place.

    python scripts/run_pages.py             # all pages
    python scripts/run_pages.py --changed   # only pages touched since origin/main
    python scripts/run_pages.py 900         # all pages, 900 s timeout

Each notebook runs with its own directory as the working directory so that
page-relative data paths behave the same as they do locally.
"""

from __future__ import annotations

import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError, DeadKernelError

ROOT = Path(__file__).resolve().parents[1]
ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
TIMEOUT = int(ARGS[0]) if ARGS else 1800

# The arithmetic is pinned here and not only in CI, so that a local run, the
# baselines and CI compute the same way (docs/claims-design.md, section 7). The
# kernels inherit this environment. A value already set is kept: the portability
# job sets OPENBLAS_CORETYPE=Sandybridge to find claims that hold on one machine only.
# NumPy picks its own SIMD code per CPU, independent of OpenBLAS: on an AVX-512 CPU
# exp, log and pow come from other code and differ in the last bits. Disabling the
# AVX-512 targets gives every x86-64 machine the AVX2 (X86_V3) path the baselines
# were made with (docs/claims-design.md, section 7). Names a CPU lacks are ignored.
# OPENBLAS_L2_SIZE: OpenBLAS 0.3.34 enlarges the GEMM blocks on AMD CPUs that expose
# AVX-512 and report a 1024 KB L2 (Zen 4, Zen 5), for every core table, including
# the forced Haswell one, whose kernel then overflows its stack buffer: H1.7 and J1.3
# segfaulted in numpy.linalg.solve on EPYC 9V74 and 9V45 runners (OpenBLAS #6013).
# The Haswell, Sandybridge and Nehalem tables read L2 nowhere else, so 512 only
# switches that override off and keeps the blocking of the baselines everywhere.
NPY_PIN = "X86_V4 AVX512_ICL AVX512_SPR"
PINNED_ENV = {"OPENBLAS_CORETYPE": "Haswell", "OPENBLAS_NUM_THREADS": "1",
              "OMP_NUM_THREADS": "1", "OPENBLAS_L2_SIZE": "512"}
if platform.machine().lower() in ("x86_64", "amd64"):
    PINNED_ENV["NPY_DISABLE_CPU_FEATURES"] = NPY_PIN
for _k, _v in PINNED_ENV.items():
    os.environ.setdefault(_k, _v)

# Every page runs from this checkout, so nothing may be fetched from GitHub: a
# file deleted in a PR must fail here, not load main's copy (gallery_utils and
# shared/gallery_bootstrap.py honour this; docs/page-infrastructure.md, section 2).
os.environ.setdefault("GALLERY_OFFLINE", "1")
# A kernel that dies in C code leaves only "Kernel died"; with the fault handler
# on it prints the Python stack of every thread to the log first.
os.environ.setdefault("PYTHONFAULTHANDLER", "1")


def run(path: Path, retry: bool = True) -> tuple[bool, float, str]:
    nb = nbformat.read(path, as_version=4)
    started = []                    # cells begun, to say where a kernel died
    client = NotebookClient(
        nb, timeout=TIMEOUT, kernel_name="python3",
        resources={"metadata": {"path": str(path.parent)}},
        on_cell_start=lambda cell, cell_index: started.append(cell_index),
        # leave the kernel's stderr on ours, so a fault-handler stack reaches the log
        extra_arguments=["--IPKernelApp.capture_fd_output=False"],
    )
    t0 = time.perf_counter()
    try:
        client.execute()
    except CellExecutionError as exc:
        return False, time.perf_counter() - t0, str(exc).splitlines()[-1][:200]
    except (DeadKernelError, RuntimeError, TimeoutError) as exc:
        # the kernel, not the page: one retry, then report it as this page's
        # failure instead of aborting every page after it
        where = f" in cell {started[-1]}" if started else ""
        if retry:
            print(f"     {path.parent.name}: {type(exc).__name__}{where}, retrying once")
            return run(path, retry=False)
        return False, time.perf_counter() - t0, f"{type(exc).__name__}{where}: {exc}"[:200]
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

    print("arithmetic: " + ", ".join(f"{k}={os.environ[k]}" for k in PINNED_ENV)
          + f"; GALLERY_OFFLINE={os.environ['GALLERY_OFFLINE']}\n")
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
