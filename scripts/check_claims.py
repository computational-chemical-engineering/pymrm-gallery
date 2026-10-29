#!/usr/bin/env python3
"""Check the claims of pages migrated to agreement.json format 2.

    python scripts/check_claims.py

For every page whose committed ``agreement.json`` is format 2
(docs/claims-design.md):

1. **Admissible values.** A ``value`` claim must exceed its measured noise:
   ``rtol*|x0| + atol >= 100 * spread``, with the spread from
   ``docs/reproducibility/spread.json``. A value that fails this is a bound or
   needs better numerics.
2. **Quoted precision.** Every prose literal (``quoted``, and README numbers
   tagged ``<!-- claim: name -->`` right after the literal) must have
   ``unit >= 10 * spread``.
3. **Tagged README numbers** pass the same rule as ``Claims.says``.
4. **No ABS_FLOOR prose.** A migrated page protects every claim by its kind,
   so text saying a metric is below ``ABS_FLOOR`` or outside CI is stale.

Pages still in format 1 are skipped.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "shared"))
from gallery_utils import _NUM, _says_ok, parse_literal  # noqa: E402

SPREAD = ROOT / "docs" / "reproducibility" / "spread.json"
PHRASE = r"(?:(?:below|under|less than|at most|above|over|more than|at least)\s+)?"
TAG = re.compile(rf"({PHRASE}{_NUM}\s*%?)\s*<!--\s*claim:\s*([\w.\-]+)\s*-->", re.I)
STALE = re.compile(r"ABS_FLOOR|outside (?:the )?CI|below (?:the )?(?:CI'?s? )?floor", re.I)


def spreads() -> dict:
    if not SPREAD.is_file():
        return {}
    return json.loads(SPREAD.read_text(encoding="utf-8"))["pages"]


def check_page(page_dir: Path, spread: dict) -> list[str]:
    path = page_dir / "agreement.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("format") != 2:
        return []
    page, problems = page_dir.name, []
    metrics, sp = data["metrics"], spread.get(page_dir.name, {})
    if data.get("breaks_required") is False:
        print(f"WARN {page}: no break table yet (breaks_required: false)")

    for name, c in metrics.items():
        s = sp.get(name, {}).get("abs")
        if c["kind"] == "value" and s is not None:
            allowed = c["rtol"] * abs(c["value"]) + c.get("atol", 0.0)
            if allowed < 100 * s:
                problems.append(f"{page}: value {name} tolerance {allowed:.3g} is not 100x its "
                                f"measured spread {s:.3g}; make it a bound or fix the numerics")

    quoted = [(name, lit, "notebook") for name, lits in data.get("quoted", {}).items()
              for lit in lits]
    readme = page_dir / "README.md"
    if readme.is_file():
        for m in TAG.finditer(readme.read_text(encoding="utf-8")):
            lit, name = m.group(1), m.group(2)
            if name not in metrics:
                problems.append(f"{page}: README tags unknown claim {name!r}")
                continue
            if not _says_ok(parse_literal(lit), metrics[name]):
                problems.append(f"{page}: README says {lit!r} for {name}, "
                                f"computed {metrics[name]['value']!r}")
            quoted.append((name, lit, "README"))
    for name, lit, where in quoted:
        s = sp.get(name, {}).get("abs")
        if s is None:
            continue
        unit = parse_literal(lit)["unit"]
        if parse_literal(lit)["op"] == "eq" and unit < 10 * s:
            problems.append(f"{page}: {where} quotes {lit!r} for {name}; its last place "
                            f"{unit:.1g} is below 10x the measured spread {s:.3g}")

    for f in (readme, page_dir / "meta.yaml"):
        if f.is_file():
            for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
                if STALE.search(line):
                    problems.append(f"{page}: {f.name}:{i} mentions ABS_FLOOR or 'outside CI', "
                                    "stale on a format 2 page")
    return problems


def main() -> int:
    spread = spreads()
    problems, migrated = [], 0
    for page_dir in sorted((ROOT / "pages").glob("*/")):
        if not (page_dir / "agreement.json").is_file():
            continue
        problems += check_page(page_dir, spread)
        migrated += json.loads((page_dir / "agreement.json").read_text()).get("format") == 2
    for p in problems:
        print(f"FAIL {p}")
    print(f"{migrated} page(s) in format 2 checked, {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
