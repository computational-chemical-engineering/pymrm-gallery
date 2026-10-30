"""Tests for format 2 and the loosening check in scripts/check_agreement.py."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_agreement", ROOT / "scripts" / "check_agreement.py")
ca = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ca)


def v2(**metrics):
    return {"format": 2, "metrics": metrics, "quoted": {}, "breaks": {}}


VAL = {"kind": "value", "value": 1.0, "rtol": 0.05}
BND = {"kind": "bound", "value": 1e-11, "max": 1e-8, "basis": "tolerance"}
EXA = {"kind": "exact", "value": 9}


def run_v2(base, cur):
    failures = []
    ca.check_v2("P", base, cur, failures)
    return failures


def test_v2_kinds():
    base = v2(a=VAL, b=BND, c=EXA)
    assert run_v2(base, v2(a={**VAL, "value": 1.04}, b={**BND, "value": 9e-9}, c=EXA)) == []
    f = run_v2(base, v2(a={**VAL, "value": 1.06}, b={**BND, "value": 2e-8}, c={**EXA, "value": 8}))
    assert len(f) == 3


def test_v2_name_set_is_fixed():
    f = run_v2(v2(a=VAL), v2(a=VAL, extra=VAL))
    assert any("not in the committed file" in s for s in f)
    f = run_v2(v2(a=VAL, b=BND), v2(a=VAL))
    assert any("disappeared" in s for s in f)


def test_v1_unchanged():
    failures = []
    ca.check_v1("P", {"metrics": {"a": 1.0, "z": 1e-14}}, {"metrics": {"a": 1.04, "z": 5e-13}}, failures)
    assert failures == []
    ca.check_v1("P", {"metrics": {"a": 1.0}}, {"metrics": {"a": 1.06}}, failures)
    assert len(failures) == 1


def test_loosenings():
    old = v2(a=VAL, b=BND, c={**EXA, "margin": 0.3, "min_margin": 0.1})
    old["quoted"] = {"a": ["1.000"]}
    old["breaks"] = {"a": ["sign flip"]}
    assert ca.loosenings("P", old, old) == []
    new = v2(a={**VAL, "rtol": 0.1}, b={**BND, "max": 1e-6}, c={**EXA, "min_margin": 0.01})
    new["quoted"] = {"a": ["1.0"]}
    whats = {w for _, w in ca.loosenings("P", old, new)}
    assert {"tolerance widened", "max widened", "min_margin lowered", "break row removed"} <= whats
    assert any("coarser" in w for w in whats)
    assert ca.loosenings("P", old, v2(a=VAL)) [-1] == ("c", "claim removed")
