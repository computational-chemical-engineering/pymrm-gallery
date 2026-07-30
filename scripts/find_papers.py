#!/usr/bin/env python3
"""Resolve DOIs for queued cases and check whether the paper is open access.

    python scripts/find_papers.py --limit 40          # resolve + check
    python scripts/find_papers.py --limit 40 --fetch  # also download OA PDFs

A case should only ever reach the maintainer's papers dashboard after this has
tried and failed, so that the list they are asked to upload is genuinely the
irreducible remainder.

Order of attack per case:
  1. already on disk in ~/papers/pymrm-gallery/  (matched on DOI, then author+year)
  2. CrossRef  -> DOI, from the catalogue citation string
  3. Unpaywall -> a legitimate open-access PDF link for that DOI
  4. otherwise -> needs-paper, with the citation and DOI recorded for the human

Contact address for the CrossRef and Unpaywall polite pools is read from
$CROSSREF_MAILTO or git config user.email — never hardcoded, the repo is public.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
QUEUE = ROOT / "queue_cases"
PAPERS = Path.home() / "papers" / "pymrm-gallery"
UA = "pymrm-gallery/1.0 (https://github.com/computational-chemical-engineering/pymrm-gallery)"


def contact() -> str:
    import os
    m = os.environ.get("CROSSREF_MAILTO")
    if m:
        return m
    try:
        return subprocess.run(["git", "config", "user.email"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""


def get_json(url: str, timeout: int = 25):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:  # noqa: S310
        return json.loads(fh.read().decode("utf-8", "replace"))


def crossref_doi(citation: str, title: str, mail: str):
    """Best-effort DOI. Returns (doi, score, matched_title, year).

    CrossRef always returns *something*, so the caller must filter. The single
    most effective filter is the publication year: a surname match against a
    modern paper that merely cites the classic is the characteristic failure
    ("Carman (1937)" once resolved to a 2025 permeability paper).
    """
    q = f"{title} {citation}".strip()
    if not q:
        return None, 0.0, "", None
    url = ("https://api.crossref.org/works?rows=3&query.bibliographic="
           + urllib.parse.quote(q) + (f"&mailto={urllib.parse.quote(mail)}" if mail else ""))
    try:
        items = get_json(url).get("message", {}).get("items", [])
    except Exception:
        return None, 0.0, "", None
    if not items:
        return None, 0.0, "", None
    it = items[0]
    yr = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
    return it.get("DOI"), float(it.get("score", 0)), (it.get("title") or [""])[0], yr


def unpaywall(doi: str, mail: str):
    if not doi or not mail:
        return None
    try:
        d = get_json(f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}?email="
                     + urllib.parse.quote(mail))
    except Exception:
        return None
    if not d.get("is_oa"):
        return None
    loc = d.get("best_oa_location") or {}
    return loc.get("url_for_pdf") or loc.get("url")


MANUAL = ROOT / "docs" / "papers-on-disk.yaml"


def on_disk(entry) -> str | None:
    """Match a PDF already in the papers directory.

    The hand-maintained map is consulted first and deliberately wins: many of
    these files are named by publisher PII with no author or year in the name,
    so nothing can be inferred from the filename. Reporting needs-paper for a
    paper the maintainer already supplied is the one failure mode worth extra
    machinery to avoid.
    """
    if not PAPERS.is_dir():
        return None
    if MANUAL.is_file():
        try:
            m = yaml.safe_load(MANUAL.read_text(encoding="utf-8")) or {}
            name = m.get(entry["id"])
            if name and (PAPERS / name).is_file():
                return str(PAPERS / name)
        except Exception:
            pass
    ref = entry.get("reference") or {}
    doi = str(ref.get("doi") or "").lower()
    pdfs = sorted(PAPERS.glob("*.pdf"))
    if doi:
        tail = re.sub(r"[^a-z0-9]", "", doi.split("/")[-1])
        for p in pdfs:
            if tail and tail in re.sub(r"[^a-z0-9]", "", p.name.lower()):
                return str(p)
    authors = ref.get("authors") or []
    year = str(ref.get("year") or "")
    if authors and year:
        surname = re.sub(r"[^a-z]", "", authors[0].split(",")[0].lower())
        for p in pdfs:
            n = p.name.lower()
            if surname and surname in n and year in n:
                return str(p)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=30)
    ap.add_argument("--fetch", action="store_true", help="download open-access PDFs")
    ap.add_argument("--status", default="unclaimed",
                    help="only touch cases in this status (comma separated)")
    args = ap.parse_args()
    mail = contact()
    if not mail:
        print("no contact address (set CROSSREF_MAILTO or git config user.email); "
              "Unpaywall needs one")
    want = set(args.status.split(","))

    todo = []
    for p in sorted(QUEUE.glob("*.yaml")):
        e = yaml.safe_load(p.read_text(encoding="utf-8"))
        if not e or e.get("status") not in want:
            continue
        if (e.get("paper") or {}).get("checked"):
            continue
        todo.append((p, e))
    todo = todo[:args.limit]
    print(f"checking {len(todo)} cases\n")

    n_disk = n_oa = n_doi = n_none = 0
    for p, e in todo:
        ref = e.get("reference") or {}
        doi = ref.get("doi")
        hit = on_disk(e)
        oa = None
        if hit:
            n_disk += 1
            note = f"on disk: {Path(hit).name}"
        else:
            if not doi:
                doi, score, matched, cyear = crossref_doi(e.get("catalog_reference", ""),
                                                          e.get("title", ""), mail)
                # CrossRef will confidently return something for anything; only keep
                # a hit that also looks like the right title.
                # CrossRef returns something confident-looking for any query, so a
                # score alone is not a filter. Require the returned title to share
                # real words with ours before believing it.
                if doi:
                    ours = set(re.findall(r"[a-z]{4,}", (e.get("title") or "").lower()))
                    theirs = set(re.findall(r"[a-z]{4,}", (matched or "").lower()))
                    want_year = re.search(r"(1[89]\d\d|20\d\d)",
                                          e.get("catalog_reference", "") or "")
                    year_ok = (not want_year or not cyear
                               or abs(cyear - int(want_year.group(1))) <= 1)
                    if score < 40 or not (ours & theirs) or not year_ok:
                        doi = None
                    else:
                        e["doi_source"] = "crossref-auto"   # unverified; flagged on the dashboard
                if doi:
                    n_doi += 1
                    ref = dict(ref); ref["doi"] = doi
                    e["reference"] = ref
                time.sleep(0.4)
            oa = unpaywall(doi, mail) if doi else None
            if oa:
                n_oa += 1
                note = f"open access: {oa[:70]}"
                if args.fetch:
                    try:
                        PAPERS.mkdir(parents=True, exist_ok=True)
                        dest = PAPERS / f"{e['id'].replace('.', '_')}.pdf"
                        req = urllib.request.Request(oa, headers={"User-Agent": UA})
                        with urllib.request.urlopen(req, timeout=60) as fh:  # noqa: S310
                            blob = fh.read()
                        if blob[:4] == b"%PDF":
                            dest.write_bytes(blob)
                            hit = str(dest)
                            note += f"  -> {dest.name}"
                        else:
                            note += "  (not a PDF, skipped)"
                    except Exception as exc:
                        note += f"  (download failed: {type(exc).__name__})"
            else:
                n_none += 1
                note = "no OA route" + ("" if doi else ", no DOI either")
                e["status"] = "needs-paper"
                e["blocker"] = dict(
                    kind="paper",
                    question="",
                    detail=("Not in ~/papers/pymrm-gallery/ and no open-access copy found via "
                            "Unpaywall." + ("" if doi else " CrossRef could not resolve a DOI "
                                            "confidently from the catalogue citation either, so "
                                            "the citation below is as the catalogue records it.")),
                    artifacts=[])
        e["paper"] = dict(on_disk=hit, open_access=oa, checked=str(date.today()))
        e["updated"] = str(date.today())
        p.write_text(yaml.safe_dump(e, sort_keys=False, allow_unicode=True), encoding="utf-8")
        print(f"  {e['id']:7s} {note}")

    print(f"\non disk {n_disk}   open access {n_oa}   new DOIs {n_doi}   "
          f"-> needs-paper {n_none}")


if __name__ == "__main__":
    main()
