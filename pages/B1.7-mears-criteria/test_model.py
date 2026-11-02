"""Fast checks of model.py for page B1.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (AxialVessel, Pellet, Power, cross, eta_iso1, film_heat_dev,  # noqa: E402
                   frac_unconverted)

PAGE = HERE.name


def test_pellet_against_closed_form():
    """First-order sphere at phi = 1: converged, second order in the grid."""
    err = []
    for n in (100, 200):
        p = Pellet(n)
        y, sol = p.solve(1.0, Power(1.0))
        assert sol.success
        err.append(abs(p.eta(y, 1.0) / eta_iso1(1.0) - 1))
    assert err[1] < 2e-6
    assert 1.9 <= np.log2(err[0] / err[1]) <= 2.1


def test_axial_vessel_against_closed_form():
    """First order at Pe_L = 20: van Leer TVD by newton, exit within 1e-4 of the closed form."""
    v = AxialVessel(20.0, 200)
    c, sol = v.solve(1.05)
    assert sol.success
    assert abs(v.outlet(c) / frac_unconverted(1.05, 20.0) - 1) < 1e-4


def test_stored_claims_reproduced():
    """Two committed claims, each held to its own rtol."""
    stored = load_claims(PAGE)
    s1 = cross(Power(1.0), lambda p, e, P, o, ys: P - 1.0, lo=1e-3)
    assert abs((1 - s1[2]) / stored["eq1_dev_at_bound_n1"] - 1) <= 1e-5
    assert abs(film_heat_dev(0.15, 30.0) / stored["eq14_dev_exo_gamma30"] - 1) <= 1e-9
