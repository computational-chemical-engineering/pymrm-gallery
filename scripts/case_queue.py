#!/usr/bin/env python3
"""Track every catalogued case through to a published page.

One YAML per case under docs/queue/. Agents working a case write their own
entry and nothing else shared, so several can run at once without racing on
models.yaml or the git index.

    python scripts/case_queue.py seed      # create entries for catalogued IDs
    python scripts/case_queue.py status    # counts by status
    python scripts/case_queue.py export    # docs/queue/_index.json for the dashboards

Statuses
    unclaimed    nobody has looked at it yet
    in-progress  an agent is working it now
    needs-paper  cannot proceed: the source is not on disk and not open access
    needs-input  cannot proceed: a decision or visual check is owed by the maintainer
    ready        built and green, waiting to be committed
    published    live on the site
    deferred     deliberately not built; `blocker.detail` says why
"""
from __future__ import annotations

import json
import re
import sys
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
QUEUE = ROOT / "queue_cases"
STATUSES = ["unclaimed", "in-progress", "needs-paper", "needs-input",
            "ready", "published", "deferred"]

ROW = re.compile(r"^\|\s*\**`?([A-J]\d+\.\d+)`?\**\s*\|(.+)$")


def parse_catalog():
    """{id: {title, reference, priority, tier, section}} from the catalog tables."""
    out = {}
    for f in sorted(ROOT.glob("docs/catalog-*.md")):
        for line in f.read_text(encoding="utf-8").splitlines():
            m = ROW.match(line.strip())
            if not m:
                continue
            cid, rest = m.group(1), m.group(2)
            cells = [c.strip() for c in rest.split("|")]
            clean = lambda s: re.sub(r"[*`]", "", s).strip()
            rec = dict(id=cid, section=cid[0], title=clean(cells[0]) if cells else "",
                       reference=clean(cells[1]) if len(cells) > 1 else "",
                       catalog=f.name)
            for c in cells:
                cc = clean(c)
                if re.fullmatch(r"P[123]", cc):
                    rec["priority"] = cc
                elif re.fullmatch(r"T[0-6]", cc):
                    rec["tier"] = cc
                elif re.fullmatch(r"S\d+(\+\w+)?(\+S?\d+)*", cc):
                    rec.setdefault("structures", cc)
            out[cid] = rec
    return out


def published_ids():
    idx = yaml.safe_load((ROOT / "models.yaml").read_text(encoding="utf-8"))
    return {m["id"]: m for m in idx["models"]}


def load(cid):
    p = QUEUE / f"{cid}.yaml"
    return yaml.safe_load(p.read_text(encoding="utf-8")) if p.is_file() else None


def save(rec):
    QUEUE.mkdir(exist_ok=True)
    rec["updated"] = str(date.today())
    (QUEUE / f"{rec['id']}.yaml").write_text(
        yaml.safe_dump(rec, sort_keys=False, allow_unicode=True), encoding="utf-8")


def cmd_seed():
    cat, models = parse_catalog(), published_ids()
    made = skipped = 0
    for cid, rec in sorted(cat.items()):
        if load(cid):
            skipped += 1
            continue
        m = models.get(cid, {})
        status = "unclaimed"
        if m.get("status") == "published":
            status = "published"
        elif m.get("status") == "deferred":
            status = "deferred"
        entry = dict(
            id=cid, section=rec["section"], title=rec["title"],
            catalog_reference=rec["reference"],
            priority=rec.get("priority", "P3"), tier=rec.get("tier", ""),
            structures=rec.get("structures", ""),
            status=status,
            reference=m.get("reference", {}) or {},
            paper=dict(on_disk=None, open_access=None, checked=None),
            blocker=None, notes=None,
        )
        if m.get("status") == "deferred":
            entry["blocker"] = dict(kind="policy", question="",
                                    detail=m.get("blocked_by", ""), artifacts=[])
        save(entry)
        made += 1
    print(f"seeded {made} new, {skipped} already present")


def all_entries():
    return [yaml.safe_load(p.read_text(encoding="utf-8"))
            for p in sorted(QUEUE.glob("*.yaml"))]


def cmd_status():
    from collections import Counter
    e = all_entries()
    c = Counter(x["status"] for x in e)
    print(f"{len(e)} cases")
    for s in STATUSES:
        if c.get(s):
            print(f"   {s:14s} {c[s]:4d}")
    for s in sorted(set(c) - set(STATUSES)):
        print(f"   {s:14s} {c[s]:4d}   <- unknown status")


def cmd_export():
    e = all_entries()
    (QUEUE / "_index.json").write_text(json.dumps(e, indent=1), encoding="utf-8")
    print(f"wrote {QUEUE/'_index.json'} ({len(e)} cases)")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    {"seed": cmd_seed, "status": cmd_status, "export": cmd_export}[cmd]()
