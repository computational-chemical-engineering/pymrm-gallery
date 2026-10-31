"""Fast checks of model.py for page A2.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import ClosedVessel, ThreeSection  # noqa: E402

PAGE = HERE.name
PE_B, R = 2.667, 2.0


def f_exit_d34(pe, r):
    """Danckwerts eq. (34), written independently of the page's transcription."""
    a = np.sqrt(1.0 + 4.0 * r / pe)
    return 4 * a / ((1 + a) ** 2 * np.exp(-pe / 2 * (1 - a))
                    - (1 - a) ** 2 * np.exp(-pe / 2 * (1 + a)))


def test_outlet_sign_absent_from_system():
    """With b = d = 0 at the outlet, the sign of a does not enter the assembled system."""
    a0, g0 = ClosedVessel(PE_B, R, n=50).system()
    a1, g1 = ClosedVessel(PE_B, R, n=50, outlet_a=-1.0).system()
    assert (a1 - a0).nnz == 0 or abs(a1 - a0).max() == 0.0
    assert np.max(np.abs(g1 - g0)) == 0.0


def test_newton_equals_deferred():
    """Newton on the TVD correction and deferred correction solve the same equations."""
    m = ClosedVessel(PE_B, R, n=200)
    assert abs(m.solve().c_exit - m.solve_deferred().c_exit) < 1e-13


def test_three_section_bed_matches_closed_vessel():
    """The three-section bed reproduces the closed vessel's exit without a Danckwerts bc."""
    ts = ThreeSection(1.0, PE_B, 1.0, R, n_b=100)
    _, c, _, _ = ts.solve()
    assert abs(ts.exit_value(c) / f_exit_d34(PE_B, R) - 1) < 1e-4


def test_stored_claim_reproduced():
    """The committed closed_form_rel_err_n800, held to the rtol of its claim."""
    stored = load_claims(PAGE)["closed_form_rel_err_n800"]
    dev = abs(ClosedVessel(PE_B, R, n=800).solve().c_exit / f_exit_d34(PE_B, R) - 1)
    assert abs(dev - stored) <= 0.05 * abs(stored)
