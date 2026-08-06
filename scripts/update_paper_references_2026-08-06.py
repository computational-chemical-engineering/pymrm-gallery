#!/usr/bin/env python3
"""Rewrite every reference to an old PDF filename after the 2026-08-06 rename.

Companion to `scripts/rename_papers_2026-08-06.sh`. Run it AFTER the rename.

The old→new map lives in `docs/papers-inventory.yaml` (`file` + `old_names`),
which is the authoritative catalogue of the library — this script reads it, so
there is exactly one place where a name is decided.

By default it leaves `models.yaml`, `AGENTS.md`, `pages/**` and
`docs/handoff.md` alone, because those are updated centrally. Pass
`--include-pages` to sweep them too; pass `--dry-run` to see what would change.

    python scripts/update_paper_references_2026-08-06.py --dry-run
    python scripts/update_paper_references_2026-08-06.py
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "docs" / "papers-inventory.yaml"

SKIP_DIRS = {".git", "_site", "__pycache__", ".quarto", ".ipynb_checkpoints"}
# The inventory IS the map: it must keep the old names in `old_names`, so it can
# never be a target of the rewrite. Rewriting it would erase the provenance
# chain this whole exercise exists to preserve.
NEVER = {"docs/papers-inventory.yaml",
         "scripts/update_paper_references_2026-08-06.py",
         "scripts/rename_papers_2026-08-06.sh"}
CENTRAL = {"models.yaml", "AGENTS.md", "models.json", "docs/handoff.md"}
CENTRAL_DIRS = {"pages"}
TEXT_SUFFIXES = {".yaml", ".yml", ".md", ".py", ".ipynb", ".html", ".json",
                 ".qmd", ".bib", ".txt", ".csv", ".sh"}
MAX_BYTES = 8_000_000


def load_map() -> dict[str, str]:
    entries = yaml.safe_load(INVENTORY.read_text(encoding="utf-8")) or []
    m: dict[str, str] = {}
    for e in entries:
        for old in e.get("old_names") or []:
            m[old] = e["file"]
    # longest first, so a name that is a substring of another cannot shadow it
    return dict(sorted(m.items(), key=lambda kv: -len(kv[0])))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--include-pages", action="store_true",
                    help="also rewrite models.yaml, pages/** and docs/handoff.md")
    args = ap.parse_args()

    mapping = load_map()
    if not mapping:
        print("no old_names in the inventory — nothing to do", file=sys.stderr)
        return 1

    changed = skipped = 0
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            p = Path(dirpath) / fn
            rel = p.relative_to(ROOT).as_posix()
            if p.suffix.lower() not in TEXT_SUFFIXES:
                continue
            if p.stat().st_size > MAX_BYTES or rel in NEVER:
                continue
            central = rel in CENTRAL or rel.split("/")[0] in CENTRAL_DIRS
            if central and not args.include_pages:
                try:
                    t = p.read_text(encoding="utf-8", errors="ignore")
                except OSError:
                    continue
                if any(old in t for old in mapping):
                    skipped += 1
                    print(f"SKIPPED (central): {rel}")
                continue
            try:
                t = p.read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError):
                continue
            new = t
            for old, repl in mapping.items():
                if old in new:
                    new = new.replace(old, repl)
            if new != t:
                changed += 1
                print(("would rewrite " if args.dry_run else "rewrote ") + rel)
                if not args.dry_run:
                    p.write_text(new, encoding="utf-8")

    print(f"\n{changed} file(s) {'would be ' if args.dry_run else ''}rewritten; "
          f"{skipped} left for the central pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
