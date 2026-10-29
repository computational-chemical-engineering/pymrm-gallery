#!/usr/bin/env python3
"""Compare each page's agreement metrics against its committed baseline.

    python scripts/check_agreement.py                    # run against git HEAD
    python scripts/check_agreement.py --base origin/main # also refuse loosening

Pages write ``agreement.json`` next to the notebook, in one of two formats:

format 1 (``report_agreement``)
    ``{"metrics": {name: float}}``. Every metric must stay within ``REL_TOL``
    of the committed value; pairs both below ``ABS_FLOOR`` are not compared.
format 2 (``Claims``, docs/claims-design.md)
    ``{"format": 2, "metrics": {name: {"kind": ..., "value": ...}}}``. Each
    metric is checked by its committed kind: a ``value`` within its ``rtol`` and
    ``atol``, a ``bound`` inside its bounds, an ``exact`` equal. The set of
    names must match the committed set in both directions.

With ``--base REF`` the committed files are also compared with the merge-base of
``REF`` and ``HEAD``, and a loosening fails unless ``docs/loosened.yaml`` has an
entry for it: a widened tolerance or bound, a changed kind, a removed claim, a
removed break row, or a prose literal made coarser or dropped. A page moving
from format 1 to format 2 is listed for review, not failed.

Update a baseline deliberately (and say why in the PR) rather than loosening the
tolerance.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "shared"))
from gallery_utils import _passes, parse_literal  # noqa: E402

REL_TOL = 0.05      # format 1: 5% relative drift allowed
ABS_FLOOR = 1e-12   # format 1: metrics below this are not compared


def git_json(ref: str, path: Path) -> dict | None:
    """``path`` as committed at ``ref``, or None when it does not exist there."""
    rel = path.relative_to(ROOT)
    try:
        blob = subprocess.run(
            ["git", "show", f"{ref}:{rel.as_posix()}"],
            cwd=ROOT, capture_output=True, check=True, text=True).stdout
    except subprocess.CalledProcessError:
        return None
    return json.loads(blob)


def is_v2(d: dict | None) -> bool:
    return bool(d) and d.get("format") == 2


def check_v1(page: str, base: dict, current: dict, failures: list) -> None:
    cur = {k: (v["value"] if isinstance(v, dict) else v)
           for k, v in current.get("metrics", {}).items()}
    for key, base_val in base.get("metrics", {}).items():
        if key not in cur:
            failures.append(f"{page}: metric {key!r} disappeared")
            continue
        new_val = cur[key]
        if isinstance(new_val, dict) or isinstance(new_val, str):
            failures.append(f"{page}: {key} is not a number ({new_val!r})")
            continue
        if abs(base_val) < ABS_FLOOR and abs(new_val) < ABS_FLOOR:
            continue
        denom = max(abs(base_val), ABS_FLOOR)
        drift = abs(new_val - base_val) / denom
        if drift > REL_TOL:
            failures.append(
                f"{page}: {key} {base_val:.6g} -> {new_val:.6g} "
                f"({drift:.1%} drift, tolerance {REL_TOL:.0%})")
    for key in set(cur) - set(base.get("metrics", {})):
        print(f"NEW  {page}: metric {key!r} added")


def check_v2(page: str, base: dict, current: dict, failures: list) -> None:
    bm, cm = base["metrics"], current.get("metrics", {})
    if not is_v2(current):
        failures.append(f"{page}: committed in format 2, ran in format 1")
        return
    for key in sorted(set(bm) - set(cm)):
        failures.append(f"{page}: claim {key!r} disappeared")
    for key in sorted(set(cm) - set(bm)):
        failures.append(f"{page}: claim {key!r} is not in the committed file")
    for key in sorted(set(bm) & set(cm)):
        c, x = bm[key], cm[key]["value"]
        if cm[key]["kind"] != c["kind"]:
            failures.append(f"{page}: {key} kind {c['kind']} -> {cm[key]['kind']}")
        elif not _passes(c, c["value"], x):
            failures.append(f"{page}: {key} {c['value']!r} -> {x!r} fails its "
                            f"{c['kind']} check")
    benv, cenv = base.get("env", {}), current.get("env", {})
    diff = {k: (benv.get(k), cenv.get(k)) for k in set(benv) | set(cenv)
            if benv.get(k) != cenv.get(k)}
    if diff:
        print(f"NOTE {page}: environment differs from the baseline's: {diff}")


def coarser(old: str, new_list: list[str]) -> bool:
    """Whether every remaining literal is coarser than ``old`` or ``old`` is gone."""
    if old in new_list:
        return False
    try:
        u_old = parse_literal(old)["unit"]
    except ValueError:
        return True
    units = []
    for n in new_list:
        try:
            units.append(parse_literal(n)["unit"])
        except ValueError:
            pass
    return not units or min(units) > u_old


def loosenings(page: str, old: dict, new: dict) -> list[tuple[str, str]]:
    """``(metric, what)`` for every way ``new`` checks less than ``old``."""
    out = []
    om, nm = old["metrics"], new["metrics"]
    for key, o in om.items():
        if key not in nm:
            out.append((key, "claim removed"))
            continue
        n = nm[key]
        if n["kind"] != o["kind"]:
            out.append((key, f"kind {o['kind']} -> {n['kind']}"))
            continue
        if o["kind"] == "value":
            if n.get("rtol", 0) > o.get("rtol", 0) or n.get("atol", 0) > o.get("atol", 0):
                out.append((key, "tolerance widened"))
        elif o["kind"] == "bound":
            if o.get("max") is not None and (n.get("max") is None or n["max"] > o["max"]):
                out.append((key, "max widened"))
            if o.get("min") is not None and (n.get("min") is None or n["min"] < o["min"]):
                out.append((key, "min widened"))
            if n.get("structural") and not o.get("structural"):
                out.append((key, "made structural"))
        elif o["kind"] == "exact":
            if o.get("min_margin") is not None and (
                    n.get("min_margin") is None or n["min_margin"] < o["min_margin"]):
                out.append((key, "min_margin lowered"))
        if set(old.get("breaks", {}).get(key, [])) - set(new.get("breaks", {}).get(key, [])):
            out.append((key, "break row removed"))
        for lit in old.get("quoted", {}).get(key, []):
            if coarser(lit, new.get("quoted", {}).get(key, [])):
                out.append((key, f"prose literal {lit!r} dropped or made coarser"))
    return out


def allowed_loosenings() -> set[tuple[str, str]]:
    path = ROOT / "docs" / "loosened.yaml"
    if not path.is_file():
        return set()
    import yaml

    entries = yaml.safe_load(path.read_text(encoding="utf-8")) or []
    for e in entries:
        if not e.get("reason"):
            raise SystemExit(f"docs/loosened.yaml: entry {e} has no reason")
    return {(str(e["page"]), str(e["metric"])) for e in entries}


def main() -> int:
    base_ref = None
    if "--base" in sys.argv:
        ref = sys.argv[sys.argv.index("--base") + 1]
        base_ref = subprocess.run(["git", "merge-base", ref, "HEAD"], cwd=ROOT,
                                  capture_output=True, check=True, text=True).stdout.strip()
    allowed = allowed_loosenings() if base_ref else set()
    failures = []
    checked = 0

    for path in sorted((ROOT / "pages").glob("*/agreement.json")):
        page = path.parent.name
        pid = page.split("-")[0]
        base = git_json("HEAD", path)
        if base is None:
            print(f"NEW  {page}: no committed baseline, skipping")
            continue
        current = json.loads(path.read_text(encoding="utf-8"))
        checked += 1
        if is_v2(base):
            check_v2(page, base, current, failures)
        else:
            check_v1(page, base, current, failures)

        if base_ref:
            old = git_json(base_ref, path)
            if old is None:
                continue
            if is_v2(old) and is_v2(base):
                for key, what in loosenings(page, old, base):
                    if (pid, key) in allowed or (page, key) in allowed:
                        print(f"LOOSENED {page}: {key}: {what} (in docs/loosened.yaml)")
                    else:
                        failures.append(f"{page}: {key}: {what}, with no entry in "
                                        "docs/loosened.yaml")
            elif is_v2(base):
                kinds = {}
                for c in base["metrics"].values():
                    kinds[c["kind"]] = kinds.get(c["kind"], 0) + 1
                print(f"MIGRATED {page}: format 1 -> 2, {kinds}; review every bound "
                      "and every value whose rtol is not 0.05")

    for f in failures:
        print(f"FAIL {f}")
    print(f"\n{checked} page(s) checked, {len(failures)} regression(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
