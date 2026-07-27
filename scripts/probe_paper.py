#!/usr/bin/env python3
"""Triage a downloaded paper: how expensive will its page be to build?

The cost of a gallery page is dominated by where its numbers live. Numbers in a
text layer transcribe in minutes; numbers in a scanned table need page images;
numbers that exist only as markers in a figure need extraction and a human
sign-off, which is the slow path. This reports which case each paper is, so a
batch can be ordered cheapest-first instead of discovered one at a time.

    python scripts/probe_paper.py ~/papers/pymrm-gallery/*.pdf
    python scripts/probe_paper.py --render 7 paper.pdf   # 600 dpi page renders

Rule of thumb from the pages built so far: below about 5,000 characters per page
the text layer will mangle sub- and superscripts, so every number has to be read
off a render regardless of what the character count promises.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

THIN_TEXT = 5000        # chars/page below which the text layer is untrustworthy


def pdftotext(path: Path, first=None, last=None) -> str:
    cmd = ["pdftotext", "-layout"]
    if first:
        cmd += ["-f", str(first), "-l", str(last or first)]
    cmd += [str(path), "-"]
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=120).stdout
    except Exception:
        return ""


def n_pages(path: Path) -> int:
    try:
        out = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True,
                             timeout=60).stdout
        m = re.search(r"^Pages:\s+(\d+)", out, re.M)
        return int(m.group(1)) if m else 0
    except Exception:
        return 0


def probe(path: Path) -> dict:
    pages = n_pages(path)
    text = pdftotext(path)
    chars = len(text)
    per_page = chars / pages if pages else 0
    tables = sorted(set(int(m) for m in re.findall(r"\bTable\s+(\d{1,2})\b", text)))
    figures = sorted(set(int(m) for m in re.findall(r"\bFig(?:ure)?\.?\s+(\d{1,2})\b", text)))
    # a table caption followed by digit-dense lines suggests the table survived
    numeric_lines = sum(1 for ln in text.splitlines()
                        if len(re.findall(r"\d", ln)) >= 8 and len(ln.split()) >= 4)
    # the tell-tale of a mangled scan: "lo-'" style exponents, stray decimal loss
    mangled = len(re.findall(r"\b\d+\s+lo[-‘'`]", text)) + len(re.findall(r"\bl0[-–]", text))
    if per_page == 0:
        verdict, cost = "no text layer at all", "page images for everything"
    elif per_page < THIN_TEXT:
        verdict, cost = "thin text layer", "read numbers off 600 dpi renders"
    elif mangled:
        verdict, cost = "text layer mangles exponents", "read numbers off renders"
    else:
        verdict, cost = "usable text layer", "tables may transcribe directly"
    return dict(name=path.name, pages=pages, per_page=int(per_page),
                tables=tables, figures=figures, numeric_lines=numeric_lines,
                mangled=mangled, verdict=verdict, cost=cost)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pdfs", nargs="+", type=Path)
    ap.add_argument("--render", type=int, metavar="PAGE",
                    help="also render this page at 600 dpi next to the PDF")
    args = ap.parse_args()

    rows = []
    for p in args.pdfs:
        if p.suffix.lower() != ".pdf" or not p.is_file():
            continue
        r = probe(p)
        rows.append(r)
        print(f"\n{p.name}")
        print(f"   {r['pages']:>3} pages · {r['per_page']:>6} chars/page · "
              f"{r['numeric_lines']:>3} digit-dense lines")
        print(f"   tables mentioned : {r['tables'] or '—'}")
        print(f"   figures mentioned: {r['figures'] or '—'}")
        print(f"   -> {r['verdict']}; {r['cost']}")
        if args.render:
            out = p.with_suffix("")
            subprocess.run(["pdftoppm", "-r", "600", "-f", str(args.render),
                            "-l", str(args.render), "-png", str(p), str(out)],
                           check=False)
            print(f"   rendered page {args.render} at 600 dpi next to the PDF")

    if len(rows) > 1:
        print("\n\ncheapest first:")
        order = {"usable text layer": 0, "text layer mangles exponents": 1,
                 "thin text layer": 2, "no text layer at all": 3}
        for r in sorted(rows, key=lambda d: (order[d["verdict"]], -d["per_page"])):
            print(f"  {r['verdict']:30s} {len(r['figures']):>2} figs  {r['name'][:46]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
