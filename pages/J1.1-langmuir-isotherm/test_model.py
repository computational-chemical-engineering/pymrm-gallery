"""Fast checks of model.py for page J1.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_fit_recovers_exact_isotherm():
    """Eq. (31) fitted to its own exact curve returns its constants, and newton and
    Brent agree on the stationarity condition."""
    p = np.array([2.0, 4.0, 8.0, 15.0, 25.0, 34.0])
    q = model.q_eq31(p, 0.156, 38.9)
    log = []
    a, b, r = model.fit_langmuir(p, q, "exact", log=log)
    assert abs(a / 0.156 - 1) < 1e-7 and abs(b / 38.9 - 1) < 1e-7 and r < 1e-6
    assert log[0][1] < 1e-9


def test_march_reaches_eq16_at_first_order():
    """Implicit Euler on eq. (14) approaches ln2/k, the error halving with the step."""
    theta1, k = 0.84, 6.3
    gaps = []
    for n in (200, 400):
        t, th = model.march_eq14(theta1, k, 6.0 / k, n)
        gaps.append(model.t_half_of_march(t, th, theta1) / (np.log(2.0) / k) - 1)
    assert 0 < gaps[1] < gaps[0] and abs(gaps[0] / gaps[1] - 2.0) < 0.05


def test_ladder_truncation_matches_stored_value():
    """The 20-layer ladder at x = 0.6, c = 37 against the BET closed form gives the
    stored caseVI_ladder_20layer_rel (its rtol, 5 %)."""
    x, c = 0.6, 37.0
    bet = c * x / ((1 - x) * (1 + (c - 1) * x))
    err = abs(model.caseVI_series(x, c, 20) / bet - 1)
    stored = load_claims(PAGE)["caseVI_ladder_20layer_rel"]
    assert abs(err - stored) <= 0.05 * stored
    assert abs(model.theta_caseVI(x, c) / bet - 1) < 1e-14
