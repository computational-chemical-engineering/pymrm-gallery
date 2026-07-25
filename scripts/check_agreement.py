#!/usr/bin/env python3
"""Compare each page's agreement metrics against its committed baseline.

Pages call ``gallery_utils.report_agreement`` which writes ``agreement.json``
next to the notebook. That file is committed and acts as the baseline. This
script re-reads it after execution and fails if any metric has drifted beyond
tolerance, which is how a pymrm change that quietly degrades a page gets caught.

Update a baseline deliberately (and say why in the PR) rather than loosening the
tolerance.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REL_TOL = 0.05      # 5% relative drift allowed
ABS_FLOOR = 1e-12   # metrics below this are numerical noise; ignore drift


def committed(path: Path) -> dict | None:
    """The baseline as committed in git HEAD, not the freshly written file."""
    rel = path.relative_to(ROOT)
    try:
        blob = subprocess.run(
            ["git", "show", f"HEAD:{rel.as_posix()}"],
            cwd=ROOT, capture_output=True, check=True, text=True).stdout
    except subprocess.CalledProcessError:
        return None
    return json.loads(blob)


def main() -> int:
    failures = []
    checked = 0

    for path in sorted((ROOT / "pages").glob("*/agreement.json")):
        page = path.parent.name
        base = committed(path)
        if base is None:
            print(f"NEW  {page}: no committed baseline, skipping")
            continue
        current = json.loads(path.read_text(encoding="utf-8"))
        checked += 1

        for key, base_val in base.get("metrics", {}).items():
            if key not in current.get("metrics", {}):
                failures.append(f"{page}: metric {key!r} disappeared")
                continue
            new_val = current["metrics"][key]
            if abs(base_val) < ABS_FLOOR and abs(new_val) < ABS_FLOOR:
                continue
            denom = max(abs(base_val), ABS_FLOOR)
            drift = abs(new_val - base_val) / denom
            if drift > REL_TOL:
                failures.append(
                    f"{page}: {key} {base_val:.6g} -> {new_val:.6g} "
                    f"({drift:.1%} drift, tolerance {REL_TOL:.0%})")

        for key in set(current.get("metrics", {})) - set(base.get("metrics", {})):
            print(f"NEW  {page}: metric {key!r} added")

    for f in failures:
        print(f"FAIL {f}")
    print(f"\n{checked} page(s) checked, {len(failures)} regression(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
