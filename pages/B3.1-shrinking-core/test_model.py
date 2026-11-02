"""Fast checks of model.py for page B3.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import Shell, eq6, flow_analytic, groups, theta_from_resistances  # noqa: E402

PAGE = HERE.name
R, D = 0.05, 1e-4


def test_eq6_is_the_integral_of_the_resistances():
    """Eq. 6 and the integrated moving boundary agree to round-off, endpoints exact."""
    x = np.linspace(0, 1, 101)
    om, ga = groups(1e-2, 1e-3, D, R)
    assert np.abs(eq6(x, om, ga) - theta_from_resistances(x, 1e-2, 1e-3, D, R)).max() < 1e-13
    assert abs(eq6(1.0, om, ga)) < 1e-14 and abs(eq6(0.0, om, ga) - 1.0) < 1e-14


def test_shell_solve_against_series_resistance():
    """One shell solve on 200 cells: flow within 1e-5 of the series resistance, conserved."""
    rc = 0.6 * R
    s = Shell(rc, R, D, 1e-2, 1e-2, 200)
    C = s.solve()
    W = s.flow_core(C)
    assert abs(W / flow_analytic(rc, R, D, 1e-2, 1e-2) - 1) < 1e-5
    assert abs(W - s.flow_outer(C)) < 1e-9 * W
    assert abs(W - s.flow_core(C, read="cells")) < 1e-9 * W


def test_stored_claims_reproduced():
    """One committed claim, held to its own rtol."""
    stored = load_claims(PAGE)
    rc = 0.3 * R
    s = Shell(rc, R, D, 1.0, 1e-3, 400)
    err = abs(s.flow_core(s.solve()) / flow_analytic(rc, R, D, 1.0, 1e-3) - 1)
    assert abs(err / stored["shell_solve_worst_rel_err"] - 1) <= 0.05
