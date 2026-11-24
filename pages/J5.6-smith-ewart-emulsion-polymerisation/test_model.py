"""Fast checks of model.py for page J5.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from pymrm.checks import find_roots  # noqa: E402
from model import nbar_closed, nmax_for, p01_closed, recursion, stage_closed, stats  # noqa: E402

PAGE = HERE.name
ALPHA_T, M_T = 0.1, 0.01


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_radical_balance_and_normalisation():
    """The truncated recursion conserves particles and radicals."""
    _, p = nmax_for(ALPHA_T, M_T)
    s = stats(p)
    assert abs(s["sum"] - 1) < 1e-14
    assert abs(ALPHA_T * (1 - p[-1]) - M_T * s["nbar"] - 2 * s["S2"]) < 1e-15


def test_closed_form_against_recursion():
    """At T the closed form meets the converged recursion in nbar, p0 and p1."""
    _, p = nmax_for(ALPHA_T, M_T)
    s = stats(p)
    p0, p1 = p01_closed(ALPHA_T, M_T)
    assert abs(s["nbar"] / nbar_closed(ALPHA_T, M_T) - 1) < 1e-12
    assert abs(s["p0"] / p0 - 1) < 1e-12 and abs(s["p1"] / p1 - 1) < 1e-12
    ref, rtol = check_value("T_nbar")
    assert abs(nbar_closed(ALPHA_T, M_T) / ref - 1) < rtol


def test_stage_check_value():
    """The lower-bound constant k of eq. 26a meets the second route's value."""
    c = stage_closed(find_roots=find_roots)
    ref, rtol = check_value("k_lower")
    assert abs(c["k_lower"] / ref - 1) < rtol


def test_stored_claim():
    """A stored claim of the page is reproduced: the truncation error at n_max = 3."""
    stored = load_claims(PAGE)
    err = stats(recursion(ALPHA_T, M_T, 3))["nbar"] - nbar_closed(ALPHA_T, M_T)
    assert np.isclose(err, stored["v1_trunc_err_3"], rtol=1e-6, atol=0)
