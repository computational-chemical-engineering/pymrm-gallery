"""Fast checks of model.py for page B1.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (Pellet, Power, WeiszHicks, shoot, eta_exact, threshold_exact,  # noqa: E402
                   phi0_eq10)

PAGE = HERE.name


def test_pymrm_against_closed_form():
    """First-order sphere at phi_0 = 1: converged, second order in the grid."""
    err = []
    for n in (100, 200):
        p = Pellet("sphere", n)
        y, sol = p.solve(1.0, Power(1))
        assert sol.success
        err.append(abs(p.eta(y, 1.0) / eta_exact(2, 1.0) - 1))
    assert err[1] < 2e-6
    assert 1.9 <= np.log2(err[0] / err[1]) <= 2.1


def test_shooting_against_pymrm_nonisothermal():
    """One ignited state of eq. (9): shooting and a seeded pymrm solve agree."""
    phi, eta, _, sol = shoot(2, WeiszHicks(0.6, 20.0), 2.0, rtol=1e-11, dense=True)
    p = Pellet("sphere", 400)
    y0 = np.clip(sol.sol(np.clip(p.u_c * phi, sol.t[0], phi))[0], 0.0, 1.0)
    y, _ = p.solve(phi, WeiszHicks(0.6, 20.0), y_init=y0)
    assert abs(p.eta(y, phi) / eta - 1) < 1e-3


def test_stored_claims_reproduced():
    """Two committed claims, each held to its own rtol."""
    stored = load_claims(PAGE)
    assert abs(phi0_eq10(0.3, 20.0) / stored["eq10_worked_case_phi0"] - 1) <= 1e-6
    eta_sphere = threshold_exact(2)[1]
    assert abs(eta_sphere / stored["eta_at_Phi1_sphere_first_order"] - 1) <= 1e-6
