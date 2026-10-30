"""Tests for the Claims helper in shared/gallery_utils.py (docs/claims-design.md)."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
from gallery_utils import Claims, ClaimsFailed, load_claims, parse_literal  # noqa: E402


@pytest.mark.parametrize("text, op, nums, unit", [
    ("0.7186", "eq", [0.7186], 1e-4),
    ("3.9e-5", "eq", [3.9e-5], 1e-6),
    ("40", "eq", [40.0], 1.0),
    ("2.5 %", "eq", [0.025], 1e-3),
    ("below 1e-8", "max", [1e-8], 1e-8),
    ("less than 0.1%", "max", [0.001], 1e-3),
    ("at least 1.8", "min", [1.8], 0.1),
    ("between 1.8 and 2.2", "between", [1.8, 2.2], 0.1),
])
def test_parse_literal(text, op, nums, unit):
    lit = parse_literal(text)
    assert lit["op"] == op
    assert lit["nums"] == pytest.approx(nums)
    assert lit["unit"] == pytest.approx(unit)


def test_parse_literal_rejects_words():
    with pytest.raises(ValueError):
        parse_literal("an order of magnitude")


def _one(kind="value", x=0.7186, **kw):
    c = Claims("T0.0")
    getattr(c, kind)("m", x, **kw)
    return c


@pytest.mark.parametrize("literal, x, ok", [
    ("0.7186", 0.7186, True),
    ("0.7186", 0.71865, True),    # rounding boundary: half a unit would flip here
    ("0.7186", 0.71869, True),
    ("0.7186", 0.7188, False),
    ("0.7186", 0.69, False),      # inside 5 % CI tolerance, but not what the prose says
    ("0.72", 0.7186, True),
])
def test_says_value(literal, x, ok):
    assert _one(x=x).says(literal, "m") is ok


def test_says_bound_phrase_must_be_as_loose_as_the_bound():
    c = _one("bound", 4e-11, max=1e-8, basis="tolerance")
    assert c.says("below 1e-8", "m")
    assert c.says("below 1e-6", "m")
    assert not c.says("below 1e-10", "m")   # true on this run, not guaranteed


def test_says_value_phrase_needs_whole_band():
    c = _one(x=0.098)
    assert not c.says("below 0.1", "m")     # 0.098 * 1.05 > 0.1
    assert c.says("below 0.2", "m")


def test_bound_violation_raises_at_report():
    c = _one("bound", 2e-8, max=1e-8, basis="tolerance", structural=True)
    with pytest.raises(ClaimsFailed, match="violates"):
        c.report(write=False)


def test_missing_break_row_raises():
    c = _one()
    with pytest.raises(ClaimsFailed, match="no break row"):
        c.report(write=False)
    assert c.breaks("off by 1 %", "m", 0.7186 * 1.01) is False   # inside 5 %: not caught
    with pytest.raises(ClaimsFailed, match="no break row"):
        c.report(write=False)
    assert c.breaks("sign flip", "m", -0.7186) is True
    c.report(write=False)


def test_identity_needs_its_own_bound_broken():
    c = _one("bound", 1e-16, max=1e-12, basis="identity")
    assert c.breaks("moves it a little", "m", 1e-13) is False
    assert c.breaks("breaks the identity", "m", 1e-6) is True
    c.report(write=False)


def test_count_margin():
    c = Claims("T0.0")
    n = c.count("caught", [1e-3, 0.5, 2.0], threshold=0.1, min_margin=0.05)
    assert n == 2
    c.breaks("all caught", "caught", 3)
    c.report(write=False)
    c = Claims("T0.0")
    c.count("caught", [0.1001, 0.5], threshold=0.1, min_margin=0.05)
    with pytest.raises(ClaimsFailed, match="decided by noise"):
        c.report(write=False, require_breaks=False)


def test_twice_and_zero_value():
    c = _one()
    with pytest.raises(ValueError):
        c.value("m", 1.0)
    with pytest.raises(ValueError):
        Claims("T").value("z", 0.0)
    with pytest.raises(ValueError):
        Claims("T").bound("b", 1.0, max=2.0, basis="chosen")


def test_report_writes_format_2_and_load_claims_reads_both(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    c = _one()
    c.breaks("sign flip", "m", -1.0)
    c.says("0.7186", "m")
    c.report(write=True)
    data = json.loads((tmp_path / "agreement.json").read_text())
    assert data["format"] == 2 and data["metrics"]["m"]["kind"] == "value"
    assert data["quoted"] == {"m": ["0.7186"]}
    assert data["breaks"] == {"m": ["sign flip"]}
    for page, payload in [("P2", data), ("P1", {"page": "P1", "metrics": {"m": 0.5}})]:
        (tmp_path / "pages" / page).mkdir(parents=True)
        (tmp_path / "pages" / page / "agreement.json").write_text(json.dumps(payload))
    assert load_claims("P2") == {"m": 0.7186}
    assert load_claims("P1") == {"m": 0.5}


def test_says_exact_integer():
    c = Claims("T0.0")
    c.exact("n", 10)
    assert c.says("10", "n")
    assert not c.says("9", "n")
    assert c.says("between 9 and 11", "n")


def test_unbroken_recorded_when_breaks_not_required():
    c = _one()
    payload = c.report(write=False, require_breaks=False)
    assert payload["breaks_required"] is False and payload["unbroken"] == ["m"]


def test_report_writes_only_in_a_page_folder(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    c = _one()
    c.breaks("sign flip", "m", -1.0)
    c.report()
    assert not (tmp_path / "agreement.json").exists()      # a reader's folder
    (tmp_path / "meta.yaml").write_text("id: X\n")
    (tmp_path / "index.py").write_text("# %%\n")
    c.report()
    assert (tmp_path / "agreement.json").exists()          # a page folder
