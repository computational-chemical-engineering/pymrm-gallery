"""Fast checks of model.py for page B1.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import FirstOrder, TwoFieldPellet, ReducedPellet, violation  # noqa: E402

PAGE = HERE.name


def test_identity_holds_on_a_two_field_solve():
    """The relation is never used by the solver, yet holds to round-off."""
    pel = TwoFieldPellet("sphere", 100)
    v, res = pel.solve(1.0, 0.3, FirstOrder(20.0))
    assert res.steady.success
    assert violation(v, 0.3) < 1e-10
    assert pel.backward_error(v, 1.0, 0.3, FirstOrder(20.0)) < 1e-12


def test_heat_scale_defect_response():
    """A heat source scaled by 1 + delta gives eps = delta (1 - y_min)."""
    pel = TwoFieldPellet("slab", 100)
    v, _ = pel.solve(1.0, 0.3, FirstOrder(20.0), heat_factor=1.01)
    assert abs(violation(v, 0.3) / (0.01 * (1.0 - v[:, 0].min())) - 1.0) < 1e-9


def test_reduced_equation_matches_two_field():
    """The B1.1 single equation has the two-field y as its root."""
    pel = TwoFieldPellet("cylinder", 100)
    v, _ = pel.solve(0.5, 0.6, FirstOrder(20.0))
    y_red = ReducedPellet("cylinder", 100).solve(0.5, 0.6, 20.0, v[:, :1])
    assert np.max(np.abs(y_red[:, 0] - v[:, 0])) < 1e-10


def test_stored_claim_reproduced():
    """The committed baseline_eta, held to its claim's rtol."""
    stored = load_claims(PAGE)["baseline_eta"]
    pel = TwoFieldPellet("sphere", 200)
    v, _ = pel.solve(2.0, 0.3, FirstOrder(20.0))
    assert abs(pel.eta(v, 2.0) - stored) <= 1e-6 * abs(stored)
