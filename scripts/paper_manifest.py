#!/usr/bin/env python3
"""What papers does the gallery still need, and where does each one live?

Turns the fuzzy question "which paper do you need next" into a fixed checklist:
for every catalog entry that has no page yet, resolve a DOI (via CrossRef when
models.yaml does not carry one), say whether the PDF is already on disk, and
emit a click-through page of publisher links.

    python scripts/paper_manifest.py                 # checklist to stdout
    python scripts/paper_manifest.py --html out.html # click-through page
    python scripts/paper_manifest.py --update-yaml   # write resolved DOIs back

The intent is that acquiring papers becomes one batch errand rather than
seventeen interruptions. Nothing here downloads anything: fetching is left to a
human with a browser, or to a sanctioned TDM API once tokens exist. Driving an
authenticated publisher session automatically is what gets an institution's
whole IP range blocked, so this script deliberately stops at the link.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from html import escape
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
PAPERS = Path(os.environ.get("PYMRM_PAPERS", Path.home() / "papers" / "pymrm-gallery"))


def _mailto():
    """CrossRef's polite pool wants a contact address; never hardcode one here.

    This repository is public, so the address comes from the environment or the
    local git config at run time and is not committed.
    """
    env = os.environ.get("CROSSREF_MAILTO")
    if env:
        return env
    try:
        return subprocess.run(["git", "config", "user.email"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""


MAILTO = _mailto()

PUBLISHER = [
    (r"chemical engineering science|catalysis today|journal of membrane|"
     r"chemical engineering journal|powder technology|fuel\b|applied catalysis|"
     r"international journal of hydrogen|computers & chemical", "Elsevier", "insttoken"),
    (r"aiche|canadian journal of chemical", "Wiley", "TDM token"),
    (r"i&ec|ind\.? eng\.? chem|industrial & engineering|energy & fuels", "ACS", "no API"),
    (r"catalysis reviews|chemical engineering communications", "Taylor & Francis", "no API"),
    (r"faraday|journal of the chemical society", "RSC", "no API"),
    (r"electrochemical society", "IOP/ECS", "no API"),
    (r"proc\.? r\.? soc|proceedings of the royal", "Royal Society", "often free"),
]


def classify(container: str):
    for pat, pub, route in PUBLISHER:
        if container and re.search(pat, container, re.I):
            return pub, route
    return "other / not a journal", "manual"


def crossref_doi(authors, year, container, title_hint=""):
    """Best-effort DOI lookup. Returns (doi, matched_title) or (None, reason)."""
    q = " ".join(filter(None, [title_hint, " ".join(authors[:2]), str(year or ""), container]))
    url = ("https://api.crossref.org/works?rows=3"
           + (f"&mailto={urllib.parse.quote(MAILTO)}" if MAILTO else "")
           + "&query.bibliographic=" + urllib.parse.quote(q))
    try:
        with urllib.request.urlopen(url, timeout=25) as fh:
            items = json.load(fh)["message"]["items"]
    except Exception as exc:                                  # network, rate limit
        return None, f"lookup failed ({exc.__class__.__name__})"
    if not items:
        return None, "no CrossRef match"
    best = items[0]
    got_year = (best.get("issued", {}).get("date-parts") or [[None]])[0][0]
    if year and got_year and abs(int(got_year) - int(year)) > 1:
        return None, f"best match is {got_year}, wanted {year}"
    return best["DOI"], (best.get("title") or [""])[0]


def _norm(s):
    s = unicodedata.normalize("NFKD", s or "").encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def on_disk(entry, files, doi=None):
    """Match a catalog entry to a downloaded PDF by surname + year, or by PII.

    `doi` may be one just resolved from CrossRef: many of these files are named
    after the publisher's PII (Kunii1968-bubbling-bed-model-IECFund7-481.pdf, Wakao1978-particle-to-fluid-transfer-CES33-1375.pdf)
    and carry neither author nor year, so without it the check misses papers that
    are already on disk and sends you to fetch them again.
    """
    r = entry.get("reference") or {}
    surnames = [_norm(a.split(",")[0]) for a in r.get("authors", [])]
    year = str(r.get("year", ""))
    doi = doi or r.get("doi") or ""
    pii = re.sub(r"[^0-9A-Za-z]", "", doi.split("/")[-1]) if doi else ""
    for f in files:
        n = _norm(f.name)
        if surnames and surnames[0] and surnames[0] in n and year and year in f.name:
            return f.name
        if pii and pii[-10:] and pii[-10:] in n:
            return f.name
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--html", type=Path, help="write a click-through link page")
    ap.add_argument("--update-yaml", action="store_true",
                    help="write DOIs resolved via CrossRef back into models.yaml")
    ap.add_argument("--no-lookup", action="store_true", help="skip CrossRef")
    ap.add_argument("--json", type=Path, help="write the resolved rows as JSON")
    args = ap.parse_args()

    models = yaml.safe_load((ROOT / "models.yaml").read_text(encoding="utf-8"))
    files = sorted(PAPERS.glob("*.pdf")) if PAPERS.is_dir() else []

    rows, resolved = [], {}
    for m in models["models"]:
        if m["status"] != "planned":
            continue
        r = m.get("reference") or {}
        container = r.get("container", "")
        doi = r.get("doi")
        note = ""
        if not doi and not args.no_lookup:
            doi, note = crossref_doi(r.get("authors", []), r.get("year"),
                                     container, m.get("title", ""))
            if doi:
                resolved[m["id"]] = doi
                note = "DOI from CrossRef"
            time.sleep(0.4)                       # polite pool
        pub, route = classify(container)
        rows.append(dict(id=m["id"], pri=m.get("priority", "P?"),
                         title=m.get("title", ""), container=container,
                         authors=", ".join(r.get("authors", [])),
                         year=r.get("year", ""), doi=doi or "",
                         pub=pub, route=route, have=on_disk(m, files, doi),
                         note=note))

    rows.sort(key=lambda d: (d["pri"], d["id"]))
    need = [d for d in rows if not d["have"]]
    print(f"{len(rows)} planned entries · {len(rows) - len(need)} already on disk · "
          f"{len(need)} to fetch\n")
    print(f"{'':2s}{'pri':4s} {'id':6s} {'publisher':16s} {'route':11s} {'doi':30s} title")
    for d in rows:
        mark = "OK" if d["have"] else "  "
        print(f"{mark:2s}{d['pri']:4s} {d['id']:6s} {d['pub']:16s} {d['route']:11s} "
              f"{d['doi'] or '-':30s} {d['title'][:42]}")
        if d["note"]:
            print(f"{'':41s}{d['note']}")

    by_pub = {}
    for d in need:
        by_pub.setdefault(d["pub"], []).append(d)
    print("\nstill to fetch, by publisher:")
    for pub, ds in sorted(by_pub.items(), key=lambda kv: -len(kv[1])):
        print(f"  {len(ds):2d}  {pub:18s} ({ds[0]['route']})")

    if args.update_yaml and resolved:
        text = (ROOT / "models.yaml").read_text(encoding="utf-8")
        for cid, doi in resolved.items():
            # insert a doi line after the container line of that entry only
            pat = re.compile(rf"(- id: {re.escape(cid)}\n(?:.*\n)*?      container: .*\n)")
            text = pat.sub(lambda mo: mo.group(1) + f"      doi: {doi}\n", text, count=1)
        (ROOT / "models.yaml").write_text(text, encoding="utf-8")
        print(f"\nwrote {len(resolved)} resolved DOI(s) into models.yaml")

    if args.json:
        args.json.write_text(json.dumps(rows, indent=1, default=str), encoding="utf-8")
        print(f"\nwrote {args.json}")

    if args.html:
        parts = ["<meta charset='utf-8'><title>pymrm-gallery: papers to fetch</title>",
                 "<style>body{font:15px/1.5 system-ui;margin:2rem;max-width:60rem}"
                 "li{margin:.45rem 0}code{background:#f2f2f2;padding:0 .25rem}"
                 "h2{margin-top:1.6rem;font-size:1.05rem;color:#444}</style>",
                 "<h1>Papers still needed</h1>",
                 f"<p>Save into <code>{escape(str(PAPERS))}</code>. "
                 "Open each link, download the PDF, keep any filename.</p>"]
        for pub, ds in sorted(by_pub.items(), key=lambda kv: -len(kv[1])):
            parts.append(f"<h2>{escape(pub)} — {len(ds)}</h2><ol>")
            for d in ds:
                href = (f"https://doi.org/{d['doi']}" if d["doi"] else
                        "https://scholar.google.com/scholar?q=" +
                        urllib.parse.quote(f"{d['authors']} {d['title']} {d['year']}"))
                parts.append(
                    f"<li><a href='{escape(href)}' target='_blank'>"
                    f"<b>{escape(d['id'])}</b> {escape(d['title'][:70])}</a><br>"
                    f"<small>{escape(d['authors'][:60])} ({d['year']}) — "
                    f"{escape(d['container'][:52])}"
                    + ("" if d["doi"] else " — <i>no DOI, search link</i>") + "</small></li>")
            parts.append("</ol>")
        args.html.write_text("\n".join(parts), encoding="utf-8")
        print(f"\nwrote {args.html}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
