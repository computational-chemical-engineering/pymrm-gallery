#!/usr/bin/env python3
"""Measure how far each agreement metric moves between runs on different arithmetic.

    python scripts/measure_spread.py LABEL=DIR [LABEL=DIR ...]

Each DIR is a checkout whose pages have been executed (``pages/*/agreement.json``)
under one arithmetic: another BLAS kernel (``OPENBLAS_CORETYPE``), another CPU,
another runner. Writes ``docs/reproducibility/spread.json`` with, per page and
metric, the value in the first sample, the absolute spread ``max - min`` over
all samples and the spread relative to the largest magnitude.

``scripts/check_claims.py`` reads this file: a ``value`` claim must satisfy
``rtol*|x0| + atol >= 100 * spread`` and a quoted literal ``unit >= 10 * spread``
(docs/claims-design.md). The spread lives here, not in ``agreement.json``,
because it changes with every portability run.
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "reproducibility" / "spread.json"


def numbers(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for k, v in data["metrics"].items():
        v = v["value"] if isinstance(v, dict) else v
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            out[k] = float(v)
    return out


def main() -> int:
    samples = [a.split("=", 1) for a in sys.argv[1:]]
    if not samples or any(len(s) != 2 for s in samples):
        print(__doc__)
        return 2
    runs = {label: {p.parent.name: numbers(p)
                    for p in sorted(Path(d).glob("pages/*/agreement.json"))}
            for label, d in samples}
    first = samples[0][0]
    pages = {}
    for page, metrics in runs[first].items():
        rows = {}
        for k, v in metrics.items():
            vs = [runs[label].get(page, {}).get(k) for label, _ in samples]
            vs = [x for x in vs if x is not None]
            if len(vs) < 2:
                continue
            s = max(vs) - min(vs)
            scale = max(abs(x) for x in vs)
            rows[k] = {"value": v, "abs": s, "rel": s / scale if scale else 0.0, "n": len(vs)}
        pages[page] = rows
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"measured": date.today().isoformat(),
                               "samples": [label for label, _ in samples],
                               "pages": pages}, indent=1) + "\n", encoding="utf-8")
    n = sum(len(r) for r in pages.values())
    moved = sum(1 for r in pages.values() for m in r.values() if m["rel"] > 1e-6)
    print(f"{n} metrics on {len(pages)} pages, {moved} with relative spread above 1e-6 -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
