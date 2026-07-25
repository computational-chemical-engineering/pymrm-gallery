#!/usr/bin/env python3
"""Validate gallery metadata consistency.

Checks, in order:
  1. every page directory has meta.yaml and index.ipynb;
  2. every page's meta.yaml has a matching record in models.yaml (canonical);
  3. page metadata does not contradict the canonical record;
  4. every data/*.csv has a data/*.meta.yaml sidecar with required fields;
  5. no source PDFs are committed;
  6. page notebooks contain no Quarto-only Markdown syntax.

Exit code 0 on success, 1 on any error. Warnings do not fail the build.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CANONICAL_FIELDS = ["id", "slug", "title", "section", "status", "tier", "priority", "structures"]
SIDECAR_REQUIRED = ["dataset_id", "source", "acquisition", "columns"]
QUARTO_ONLY = re.compile(r"^:::|\{\{<\s*\w+|\{\.callout-", re.MULTILINE)

errors: list[str] = []
warnings: list[str] = []


def load_yaml(path: Path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-missing", action="store_true",
                    help="list catalog IDs in docs/ that are not yet in models.yaml")
    ap.add_argument("--emit-json", action="store_true",
                    help="write models.json next to models.yaml")
    args = ap.parse_args()

    index = load_yaml(ROOT / "models.yaml")
    records = {m["id"]: m for m in index["models"]}

    dup = len(index["models"]) - len(records)
    if dup:
        errors.append(f"models.yaml contains {dup} duplicate id(s)")

    # Deferred entries must say why, so the block is auditable and can be lifted.
    for mid, rec in records.items():
        status = rec.get("status")
        if status not in {"planned", "in-progress", "published", "deferred"}:
            errors.append(f"models.yaml {mid}: unknown status {status!r}")
        if status == "deferred" and not rec.get("blocked_by"):
            errors.append(
                f"models.yaml {mid}: status 'deferred' requires a `blocked_by` "
                f"field explaining the block and what lifts it")
        if status == "deferred" and (ROOT / "pages").joinpath(
                f"{mid}-{rec.get('slug', '')}").is_dir():
            errors.append(f"models.yaml {mid}: deferred but a page directory exists")

    # ---- pages -------------------------------------------------------------
    pages_dir = ROOT / "pages"
    for page in sorted(p for p in pages_dir.iterdir() if p.is_dir()) if pages_dir.is_dir() else []:
        rel = page.relative_to(ROOT)
        meta_path = page / "meta.yaml"
        if not meta_path.is_file():
            errors.append(f"{rel}: missing meta.yaml")
            continue
        if not (page / "index.ipynb").is_file():
            errors.append(f"{rel}: missing index.ipynb")
            continue

        meta = load_yaml(meta_path)
        pid = meta.get("id")
        if pid not in records:
            errors.append(f"{rel}: id {pid!r} has no record in models.yaml (canonical)")
            continue

        canon = records[pid]
        for field in CANONICAL_FIELDS:
            if field in meta and field in canon and meta[field] != canon[field]:
                errors.append(
                    f"{rel}: {field} = {meta[field]!r} contradicts models.yaml "
                    f"({canon[field]!r}). models.yaml is canonical."
                )
        if canon.get("status") == "published" and canon.get("page") != f"{rel}/":
            warnings.append(f"{rel}: models.yaml `page:` is {canon.get('page')!r}")

        # ---- notebook portability -----------------------------------------
        nb = json.loads((page / "index.ipynb").read_text(encoding="utf-8"))
        for i, cell in enumerate(nb.get("cells", [])):
            if cell.get("cell_type") != "markdown":
                continue
            src = "".join(cell.get("source", []))
            # the YAML front-matter cell legitimately starts with ---
            if i == 0:
                continue
            if QUARTO_ONLY.search(src):
                errors.append(
                    f"{rel}: markdown cell {i} uses Quarto-only syntax; page "
                    f"notebooks must render in Jupyter/GitHub/VS Code (AGENTS.md)"
                )

        # ---- datasets -------------------------------------------------------
        for csv in sorted((page / "data").glob("*.csv")) if (page / "data").is_dir() else []:
            sidecar = csv.with_suffix(".meta.yaml")
            if not sidecar.is_file():
                errors.append(f"{csv.relative_to(ROOT)}: missing provenance sidecar")
                continue
            side = load_yaml(sidecar)
            for field in SIDECAR_REQUIRED:
                if field not in side:
                    errors.append(f"{sidecar.relative_to(ROOT)}: missing {field!r}")
            if side.get("status") == "placeholder":
                warnings.append(
                    f"{csv.relative_to(ROOT)}: placeholder dataset — page shows "
                    f"simulation only")

    # ---- no PDFs anywhere ---------------------------------------------------
    for pdf in ROOT.rglob("*.pdf"):
        if ".git" in pdf.parts or "_site" in pdf.parts:
            continue
        errors.append(f"{pdf.relative_to(ROOT)}: source PDFs must not be committed")

    # ---- optional reports ---------------------------------------------------
    if args.report_missing:
        catalog_ids = set()
        for doc in (ROOT / "docs").glob("catalog-*.md"):
            catalog_ids |= set(re.findall(r"^\| ([A-J]\d+\.\d+)", doc.read_text(encoding="utf-8"),
                                          re.MULTILINE))
        missing = sorted(catalog_ids - set(records))
        print(f"\n{len(records)} of {len(catalog_ids)} catalogued models are in models.yaml")
        print(f"{len(missing)} not yet migrated:")
        print("  " + " ".join(missing))

    if args.emit_json:
        (ROOT / "models.json").write_text(json.dumps(index, indent=2), encoding="utf-8")
        print("wrote models.json")

    # ---- report -------------------------------------------------------------
    for w in warnings:
        print(f"WARNING  {w}")
    for e in errors:
        print(f"ERROR    {e}")
    if errors:
        print(f"\n{len(errors)} error(s), {len(warnings)} warning(s)")
        return 1
    print(f"\nmetadata OK ({len(warnings)} warning(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
