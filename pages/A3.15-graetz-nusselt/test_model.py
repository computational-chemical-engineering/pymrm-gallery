"""Fast checks of model.py for page A3.15 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import GraetzTube, amplitudes, decay, eigen  # noqa: E402

PAGE = HERE.name


def test_newton_and_deferred_reach_the_same_field():
    """On a coarse grid the two solver routes give the same discrete solution."""
    tube = GraetzTube(n_z=60, n_r=16)
    a = tube.solve()
    b = tube.solve_deferred()
    assert a["nit"] <= 10
    assert np.max(np.abs(a["theta"] - b["theta"])) < 1e-10


def test_coarse_solve_close_to_series():
    """A coarse 2-D solve already gives mu_1 to about 1e-2."""
    s = GraetzTube(n_z=60, n_r=16).solve()
    mu, _ = decay(s["x"], s["cup"])
    assert abs(mu / 2.7043644 - 1) < 2e-2


def test_stored_claims_reproduced():
    """nu_infinity and p_1 from the series route, by each claim's own tolerance."""
    stored = load_claims(PAGE)
    mu = eigen(n_modes=1)
    p_cup, _, _ = amplitudes(mu)
    assert abs(mu[0] ** 2 / 2 / stored["nu_infinity"] - 1) <= 1e-6
    assert abs(p_cup[0] / stored["p1_cup"] - 1) <= 1e-6
