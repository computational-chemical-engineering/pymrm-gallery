#!/usr/bin/env python3
"""Validate gallery metadata consistency.

Checks, in order:
  1. every page directory has meta.yaml and index.ipynb;
  2. every page's meta.yaml has a matching record in models.yaml (canonical);
  3. page metadata does not contradict the canonical record;
  4. every data/*.csv has a data/*.meta.yaml sidecar with required fields;
  5. no source PDFs are committed;
  6. page notebooks contain no Quarto-only Markdown syntax;
  7. migrated pages (an index.py or a bootstrap cell, docs/page-infrastructure.md)
     have a bootstrap equal to templates/bootstrap_cell.py apart from its
     literals, a pymrm pin from meta.yaml that the locked pymrm satisfies, PAGE
     equal to the directory, GALLERY_REF equal to --gallery-ref (default main),
     an index.py in step with index.ipynb, and no build_page.py.

Exit code 0 on success, 1 on any error. Warnings do not fail the build.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import yaml

from page_layout import (MARKER, bootstrap_position, code_sources, is_migrated, literals,
                         locked_pymrm, migrated_pages, normalise, pin_for, read_notebook,
                         template_cell)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
from gallery_bootstrap import satisfies  # noqa: E402  (the bootstrap's own comparison)

ROOT = Path(__file__).resolve().parents[1]
CANONICAL_FIELDS = ["id", "slug", "title", "section", "status", "tier", "priority", "structures"]
SIDECAR_REQUIRED = ["dataset_id", "source", "acquisition", "columns"]
QUARTO_ONLY = re.compile(r"^:::|\{\{<\s*\w+|\{\.callout-", re.MULTILINE)

errors: list[str] = []
warnings: list[str] = []


def load_yaml(path: Path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def check_migrated(page: Path, meta: dict, gallery_ref: str = "main",
                   lock_version: str | None = None) -> list[str]:
    """Errors of a migrated page: bootstrap cell, pin, index.py pairing, layout.

    Called only for pages that ``page_layout.is_migrated`` detects; an
    unmigrated page is never passed here.
    """
    rel = page.name
    out: list[str] = []
    if not (page / "index.py").is_file():
        out.append(f"{rel}: migrated page without index.py (its edited source)")
    if (page / "build_page.py").is_file():
        out.append(f"{rel}: migrated page still has build_page.py; index.py replaces it")
    nb = read_notebook(page / "index.ipynb")
    pos = bootstrap_position(nb)
    if pos is None:
        return out + [f"{rel}: no bootstrap cell (no code cell imports gallery_bootstrap)"]
    if pos != 0:
        out.append(f"{rel}: the bootstrap cell is code cell {pos}, must be the first")
    cell = code_sources(nb)[pos]
    if not MARKER.search(cell):
        out.append(f"{rel}: the bootstrap cell has no line starting `GALLERY_REF = `, "
                   f"which the release rewrite needs")

    found = literals(cell)
    for name, values in found.items():
        if len(values) != 1:
            out.append(f"{rel}: bootstrap assigns {name} {len(values)} times "
                       f"as a string literal, expected once")
    value = {n: v[0] for n, v in found.items() if len(v) == 1}

    floor = meta.get("pymrm_min_version")
    if floor is None:
        out.append(f"{rel}: meta.yaml has no pymrm_min_version, the bootstrap pin's source")
    elif "PYMRM" in value and value["PYMRM"] != pin_for(floor):
        out.append(f"{rel}: bootstrap PYMRM = {value['PYMRM']!r}, but meta.yaml "
                   f"pymrm_min_version {floor!r} requires {pin_for(floor)!r}")
    if floor is not None:
        locked = lock_version if lock_version is not None else locked_pymrm()
        if locked is None:
            out.append(f"{rel}: requirements-lock.txt pins no pymrm==")
        elif not satisfies(locked, pin_for(floor)):
            out.append(f"{rel}: the locked pymrm {locked} does not satisfy the pin "
                       f"{pin_for(floor)!r}; set meta.yaml pymrm_min_version to the "
                       f"version the page was checked with ({locked})")
    if "PAGE" in value and value["PAGE"] != page.name:
        out.append(f"{rel}: bootstrap PAGE = {value['PAGE']!r}, not the directory name")
    if "GALLERY_REF" in value and value["GALLERY_REF"] != gallery_ref:
        out.append(f"{rel}: bootstrap GALLERY_REF = {value['GALLERY_REF']!r}, "
                   f"expected {gallery_ref!r}")
    if normalise(cell) != normalise(template_cell()):
        out.append(f"{rel}: bootstrap cell differs from templates/bootstrap_cell.py "
                   f"apart from its three literals")

    if (page / "index.py").is_file():
        out.extend(check_paired(page))
    return out


def check_paired(page: Path) -> list[str]:
    """Errors when index.py and index.ipynb differ in cells (outputs ignored)."""
    import jupytext

    rel = page.name
    nb = jupytext.read(page / "index.ipynb")
    py = jupytext.read(page / "index.py", fmt="py:percent")
    out = []
    formats = nb.metadata.get("jupytext", {}).get("formats")
    if formats != "ipynb,py:percent":
        out.append(f"{rel}: index.ipynb metadata jupytext.formats is {formats!r}, "
                   f"expected 'ipynb,py:percent'")

    def cells(n):
        return [(c.cell_type, c.source.strip()) for c in n.cells]

    a, b = cells(nb), cells(py)
    if a != b:
        diff = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
        out.append(f"{rel}: index.py and index.ipynb differ ({len(b)} vs {len(a)} cells, "
                   f"first difference at cell {diff}); run `jupytext --sync`")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-missing", action="store_true",
                    help="list catalog IDs in docs/ that are not yet in models.yaml")
    ap.add_argument("--emit-json", action="store_true",
                    help="write models.json next to models.yaml")
    ap.add_argument("--gallery-ref", default="main",
                    help="GALLERY_REF expected in migrated pages (the release tag "
                         "when checking a release tree)")
    ap.add_argument("--list-migrated", action="store_true",
                    help="print the migrated page directories, one per line, and exit")
    args = ap.parse_args()

    if args.list_migrated:
        for page in migrated_pages(ROOT):
            print(page.relative_to(ROOT))
        return 0

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

        # ---- migrated pages: bootstrap and pairing ----------------------------
        if is_migrated(page):
            errors.extend(f"pages/{e}" for e in check_migrated(page, meta, args.gallery_ref))

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
        if any(part.startswith(".venv") for part in pdf.parts):
            continue          # a local virtual environment, ignored by git
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
