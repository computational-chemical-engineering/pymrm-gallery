"""Fast checks of model.py for page J4.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_marchers_reproduce_the_closed_form():
    """BatchReduced follows the Lambert W solution to O(dt); BatchMassAction keeps both balances."""
    t = np.array([0.0, 3.0, 9.0, 18.0])
    exact = model.s_closed(t, 1.0, 0.16, 0.5)
    be = model.BatchReduced(0.16, 0.5).run(1.0, 18.0, 360, t)
    assert np.max(np.abs(be - exact)) < 2e-3
    p = model.elem(1.0, 1e-2, 0.16, 0.5)
    out = model.BatchMassAction(p).run([1.0, 1e-2, 0.0, 0.0], 18.0, 180, 18)
    assert np.max(np.abs(out[:, 1] + out[:, 2] - 1e-2)) < 1e-13
    assert np.max(np.abs(out[:, 0] + out[:, 2] + out[:, 3] - 1.0)) < 1e-13


def test_max_gap_and_stored_value():
    """max_gap finds the same maximum as a fine scan; the stored rho* is a root of route 1."""
    p = model.elem(1.0, 1e-3, 0.16, 0.5)
    full = model.mass_action(p, 1e-3, np.array([0.0, 18.0]), 1.0, dense=True)
    red = model.ReducedCurve(1.0, 0.16, 0.5)
    v = model.max_gap(full, red)
    ts = np.linspace(0.0, 18.0, 20001)
    assert v >= np.max(np.abs(full(ts) - red(ts))) * (1 - 1e-9)
    stored = load_claims(PAGE)["KM_ratio_at_rho_1e4"]
    assert stored == 1.0 + 1e4
