"""Tests for scripts/check_claims.py on a fixture page."""

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_claims", ROOT / "scripts" / "check_claims.py")
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)


def page(tmp_path, readme=""):
    d = tmp_path / "X1.1-fixture"
    d.mkdir()
    (d / "agreement.json").write_text(json.dumps({"format": 2, "metrics": {
        "err": {"kind": "value", "value": 0.7186, "rtol": 0.05},
        "res": {"kind": "bound", "value": 1e-11, "max": 1e-8, "basis": "tolerance"}},
        "quoted": {"err": ["0.7186"]}}))
    (d / "README.md").write_text(readme)
    return d


def test_clean(tmp_path):
    d = page(tmp_path, "error 0.72<!-- claim: err --> and residual below 1e-8 <!-- claim: res -->")
    assert cc.check_page(d, {"X1.1-fixture": {"err": {"abs": 1e-7}}}) == []


def test_problems(tmp_path):
    d = page(tmp_path, "error 0.69<!-- claim: err -->, residual below 1e-10<!-- claim: res -->,"
                       " and one metric is below ABS_FLOOR")
    p = cc.check_page(d, {"X1.1-fixture": {"err": {"abs": 1e-3}}})
    text = "\n".join(p)
    assert "not 100x" in text                  # 5 % of 0.7186 < 100 * 1e-3
    assert "README says '0.69'" in text
    assert "README says 'below 1e-10'" in text
    assert "last place" in text                # 0.7186 quoted, spread 1e-3
    assert "ABS_FLOOR" in text
