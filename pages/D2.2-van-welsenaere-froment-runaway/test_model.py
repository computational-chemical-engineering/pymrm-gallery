"""Fast checks of model.py for page D2.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (RunawayTube, VWFCase, critical_p0, critical_points,  # noqa: E402
                   inlet_limits, t_max_peak)

PAGE = HERE.name
TW = 625.0


def _case():
    par = load_data("van-welsenaere-froment-1970-parameters.csv", page=PAGE)
    return VWFCase.from_table(zip(par.symbol, par.value))


def test_criterion1_against_example_1a():
    """Eq. 8 and the Eq. 27/28 bracket against the paper's printed Example 1(a)."""
    case = _case()
    ex = load_data("van-welsenaere-froment-1970-examples.csv", page=PAGE)
    printed = dict(zip(ex[ex.example == "1a"].quantity, ex[ex.example == "1a"].value))
    TM = t_max_peak(case, TW)
    lo, up, _ = inlet_limits(case, TM, TW)
    assert abs(TM - printed["T_cr"]) < 0.05
    assert abs(lo - printed["p0_lower"]) / printed["p0_lower"] < 1e-3
    assert abs(up - printed["p0_upper"]) / printed["p0_upper"] < 1e-3


def test_continuation_matches_independent_solve():
    """trace_inlet over a short sweep against a solve from the inlet state."""
    case = _case()
    p0s = np.linspace(0.0160, 0.0170, 6)
    res = RunawayTube(case, p0s[0], TW, n_z=400).trace_inlet(p0s)
    assert res.success, res.message
    direct = RunawayTube(case, p0s[-1], TW, n_z=400).solve().hot_spot()
    assert abs(res.xs[-1][:, 1].max() - direct) < 1e-6


def test_stored_claim_reproduced():
    """The committed reactor critical p0 at 625 K, held to its claim's rtol (1e-6)."""
    case = _case()
    TM = t_max_peak(case, TW)
    lo, _, mean = inlet_limits(case, TM, TW)
    p_crit = critical_p0(case, TW, lo, mean)
    x0 = load_claims(PAGE)["p0_critical_reactor_625K"]
    assert abs(p_crit - x0) <= 1e-6 * abs(x0)
    # and it sits above the exact (back-integrated) value by the upwind error only
    exact = critical_points(case, TW)["p0_1"]
    assert 0.0 < (p_crit - exact) / exact < 1e-3


def test_critical_p0_raises_the_upper_end():
    """Eq. 30 is an estimate: on a coarse grid at 670 K the root lies above it."""
    case = _case()
    TM = t_max_peak(case, 670.0)
    lo, _, mean = inlet_limits(case, TM, 670.0)
    p_crit = critical_p0(case, 670.0, lo, mean, n_z=300)
    assert p_crit > mean
