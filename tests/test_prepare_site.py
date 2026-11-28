"""Tests for the page metadata written by scripts/prepare_site.py."""

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_site  # noqa: E402


def test_page_metadata_year_decade_and_doi():
    ref = {"authors": ["Ranz, W. E.", "Marshall, W. R."], "year": 1952,
           "doi": "10.1000/x"}
    meta = prepare_site.page_metadata("A3.5", "2026-11-01", ref)
    assert meta["paper-year"] == "1952"
    assert meta["categories"] == ["decade:1950s"]
    assert meta["paper-citation"] == "Ranz and Marshall (1952)"
    assert meta["paper-doi"] == "10.1000/x"
    # the written file must parse back to the same values
    assert yaml.safe_load(yaml.safe_dump(meta, sort_keys=False)) == meta


def test_page_metadata_without_doi_and_many_authors():
    ref = {"authors": ["Henze, M.", "Grady, C.", "Gujer, W."], "year": 1987}
    meta = prepare_site.page_metadata("J4.8", "2026-11-01", ref)
    assert meta["paper-citation"] == "Henze et al. (1987)"
    assert "paper-doi" not in meta


def test_every_page_has_a_paper_year():
    refs = prepare_site.references(ROOT)
    pages = sorted(p.parent.name for p in (ROOT / "pages").glob("*/index.ipynb"))
    missing = [p for p in pages
               if not isinstance(refs.get(p.split("-", 1)[0], {}).get("year"), int)]
    assert not missing
