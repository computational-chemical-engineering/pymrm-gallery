"""Fast checks of model.py for page B1.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import CentrePinned, Pellet, trace_branch, wh_rate, wh_reference  # noqa: E402

PAGE = HERE.name


def test_isothermal_sphere_against_thiele():
    """Linear solve at phi = 5 against the exact sphere result, second order."""
    exact = (3 / 25.0) * (5.0 / np.tanh(5.0) - 1)
    err = [abs(Pellet("sphere", n_r=n).solve_linear(5.0)[1] - exact) / exact for n in (100, 200)]
    assert err[1] < 2e-5
    assert 1.9 <= np.log2(err[0] / err[1]) <= 2.2


def test_nonlinear_solve_matches_reference():
    """Newton at fixed phi on the Weisz-Hicks rate, against the shooting reference."""
    phi, eta_ref = wh_reference(0.2, 20.0, 0.5)
    y, eta, res = Pellet("sphere", n_r=200).solve_nonlinear(phi, wh_rate(0.2, 20.0))
    assert res.success
    assert abs(eta - eta_ref) / eta_ref < 1e-4


def test_stored_claim_reproduced():
    """The committed eta_ignited_pymrm, held to the rtol of its claim (0.05).

    The library continuation in the centre value reaches the ignited branch; a
    Newton solve at the stored probe modulus from that point reproduces it.
    """
    stored = load_claims(PAGE)
    pellet = Pellet("sphere", n_r=400)
    rate = wh_rate(0.6, 20.0)
    yc = np.concatenate([1.0 - np.geomspace(1e-6, 0.9, 30), np.geomspace(0.1, 1e-5, 40)[1:]])
    branch, xs, res = trace_branch(CentrePinned(pellet, rate), yc)
    assert res.success, res.message
    _, eta, sol = pellet.solve_nonlinear(stored["phi_probe"], rate, y_init=xs[-1][:-1])
    assert sol.success
    x0 = stored["eta_ignited_pymrm"]
    assert abs(eta - x0) <= 0.05 * abs(x0)
