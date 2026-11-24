"""Fast checks of model.py for page J2.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import PRINTED, SteadyPBE, x_from_efficiency  # noqa: E402

PAGE = HERE.name
N_T, PHI_T, S_T = 31 / 6, 0.10, 48.0


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_volume_balance():
    """Volume is exact per breakage column and at the steady state."""
    m = SteadyPBE(8, N_T, PHI_T, S_T, x_from_efficiency())
    assert np.max(np.abs(m.x @ m.B[:, 1:] / m.x[1:] - 1)) < 1e-13
    y, info = m.solve()
    assert info["newton"].success
    assert abs(m.stats(y)["vol"] - 1) < 1e-12


def test_residence_only_identity():
    """Without breakage and coalescence the steady state is the feed times theta."""
    m = SteadyPBE(8, N_T, PHI_T, S_T, x_from_efficiency(), p=replace(PRINTED, C1=0.0, C3=0.0))
    assert np.allclose(m.breakage_only(), m.feed * PRINTED.theta, rtol=1e-14, atol=0)


def test_check_value_d32():
    """The steady d32 at the test condition meets the second route's value."""
    m = SteadyPBE(8, N_T, PHI_T, S_T, x_from_efficiency())
    y, _ = m.solve()
    ref, rtol = check_value("T_q8_d32_cm")
    assert abs(m.stats(y)["d32"] / ref - 1) < rtol


def test_stored_claim():
    """A stored claim of the page is reproduced: breakage alone at the calibrated S."""
    stored = load_claims(PAGE)
    m = SteadyPBE(8, N_T, PHI_T, stored["r3_S"], x_from_efficiency())
    d32b = m.stats(m.breakage_only())["d32"] * 10
    assert abs(d32b / stored["r3_d32_brk_mm"] - 1) < 1e-6
