"""Fast checks of model.py for page A2.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import closure, run_history  # noqa: E402

PAGE = HERE.name
PE = 50.0


def test_closure_pin_independent():
    """With the grid's own mean subtracted, the pinned cell does not change k."""
    k_axis = closure(PE, n_r=100).k_over_D
    k_wall = closure(PE, n_r=100, pin=99).k_over_D
    assert abs(k_wall - k_axis) <= 1e-11 * k_axis


def test_slug_conserves_mass_and_rate_matches_variance():
    """A short coarse run conserves mass, and its per-step rate is the variance slope."""
    h = run_history(PE, 200, dt=2e-4, n_r=8, tau_end=0.05)
    assert abs(h.mass[-1] / h.mass[0] - 1) < 1e-12
    slope = np.diff(h.variance) / np.diff(h.t) / 2
    mid = 0.5 * (h.rate[1:] + h.rate[:-1])
    assert np.max(np.abs(slope[50:] / mid[50:] - 1)) < 5e-3


def test_stored_claim_reproduced():
    """The committed closure_rel_err_n200, held to the rtol of its claim."""
    stored = load_claims(PAGE)["closure_rel_err_n200"]
    err = abs(closure(PE, n_r=200).k_over_D / (PE ** 2 / 192) - 1)
    assert abs(err - stored) <= 0.05 * abs(stored)
