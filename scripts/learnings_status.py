"""Status of docs/refresh-learnings.md against the pymrm plugin.

Prints one summary line: how many rows are open, how many of those target pymrm,
and which open rows the pymrm plugin already cites as ``(gl N)`` (merged there,
not marked here). With ``--pymrm`` pointing at a pymrm checkout, the citations
are read from its ``origin/main`` (fetched first unless ``--no-fetch``); without
it, the cross-check is skipped.

    python scripts/learnings_status.py --pymrm ~/Code/pymrm_suite/pymrm

Exit status 1 when an open row is cited by pymrm, or when 10 or more rows that
target pymrm are open (the flush trigger in docs/refresh-plan.md section 5).
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEARNINGS = ROOT / "docs" / "refresh-learnings.md"
FLUSH_AT = 10
ROW = re.compile(r"^\| (\d+) \|")


def rows(text: str) -> dict:
    """Row number -> (target, status) from the learnings table."""
    out = {}
    for line in text.splitlines():
        m = ROW.match(line)
        if not m:
            continue
        cells = line.replace("\\|", "\x00").rstrip().rstrip("|").split("|")
        out[int(m.group(1))] = (cells[-2].strip(), cells[-1].strip())
    return out


def is_open(status: str) -> bool:
    return status == "open"


def targets_pymrm(target: str) -> bool:
    """Whether a row's target is in the pymrm repo (anything not marked gallery)."""
    t = target.lower().strip("`[( ")
    return not t.startswith(("gallery", "build briefs", "checklist"))


def cited(pymrm: Path, fetch: bool) -> set:
    """Row numbers the pymrm plugin cites as ``(gl N)`` on origin/main."""
    if fetch:
        subprocess.run(["git", "-C", str(pymrm), "fetch", "-q", "origin", "main"], check=True)
    grep = subprocess.run(["git", "-C", str(pymrm), "grep", "-ohE", r"gl [0-9]+(, ?[0-9]+)*",
                           "origin/main", "--", "plugins"], capture_output=True, text=True)
    return {int(n) for hit in grep.stdout.splitlines() for n in re.findall(r"\d+", hit.split(":")[-1])}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pymrm", type=Path, help="a pymrm checkout to cross-check citations against")
    ap.add_argument("--no-fetch", action="store_true", help="do not fetch the pymrm checkout first")
    a = ap.parse_args(argv)
    table = rows(LEARNINGS.read_text(encoding="utf-8"))
    open_rows = sorted(n for n, (_, s) in table.items() if is_open(s))
    open_pymrm = [n for n in open_rows if targets_pymrm(table[n][0])]
    stale = sorted(set(open_rows) & cited(a.pymrm, not a.no_fetch)) if a.pymrm else []
    flush = len(open_pymrm) >= FLUSH_AT
    print(f"learnings: {len(table)} rows, {len(open_rows)} open, {len(open_pymrm)} open for pymrm"
          + (f" (flush: {FLUSH_AT} or more)" if flush else "")
          + (f"; cited by pymrm but open: {stale}" if stale else "")
          + ("" if a.pymrm else "; pymrm citations not checked"))
    if open_pymrm:
        print("open for pymrm:", ", ".join(map(str, open_pymrm)))
    return 1 if (stale or flush) else 0


if __name__ == "__main__":
    sys.exit(main())
