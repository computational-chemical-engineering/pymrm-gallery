"""Fast checks of model.py for page A2.8 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import MaxMixedness, gamma_far, gamma_mm, mm_ode  # noqa: E402

PAGE = HERE.name


def test_newton_and_deferred_reach_the_same_solution():
    """On a coarse grid the two solver routes give the same discrete solution."""
    a = MaxMixedness(10, 2, n_cells=100).solve()
    b = MaxMixedness(10, 2, n_cells=100).solve_deferred()
    assert np.max(np.abs(a.g - b.g)) < 1e-11
    assert abs(a.balance_gap()) < 1e-12


def test_first_order_collapses_onto_closed_form():
    """For a first-order reaction maximum mixedness is (1 + K/n)^-n."""
    exact = (1 + 5 / 2) ** -2
    assert abs(mm_ode(5, 2, X=25.0, p=1.0) - exact) < 1e-12
    assert abs(gamma_mm(5, 2, p=1.0, n_cells=400) - exact) < 1e-5


def test_exponential_rtd_is_the_stirred_tank_root():
    """One tank: maximum mixedness equals the stirred-tank root."""
    assert abs(mm_ode(5, 1, X=25.0) - gamma_far(5, 1)) < 1e-12


def test_stored_claim_reproduced():
    """The committed grid error at K = 10, two tanks, within mm_pymrm_vs_dop853_max."""
    stored = load_claims(PAGE)["mm_pymrm_vs_dop853_max"]
    err = abs(gamma_mm(10, 2) - mm_ode(10, 2, rtol=1e-13, atol=1e-15))
    # K = 10 is one of the twelve cases whose worst gap is the claim
    assert 0.5 * stored <= err <= stored
