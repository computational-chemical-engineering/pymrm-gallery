"""Fast checks of model.py for page J1.10 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (Column, cumulants, delta_of, plate_height, width_red,  # noqa: E402
                   one_root)

PAGE = HERE.name


def test_eq21_limits_and_eq22_root():
    """The corrected eq (21) meets eq (24) and eq (25); as printed it tends to 2."""
    assert width_red(1e-8, "a") == pytest.approx(4.0, abs=1e-6)
    assert width_red(5.0, "a") == pytest.approx(5.0 + np.sqrt(2 * np.pi), rel=1e-5)
    assert width_red(1e-8, "delta") == pytest.approx(2.0, abs=1e-6)
    a = 0.7
    d = delta_of(a)
    assert d * np.exp(-d * d / 2) == pytest.approx((a + d) * np.exp(-(a + d) ** 2 / 2), rel=1e-12)


def test_one_root_raises_without_a_root():
    with pytest.raises(RuntimeError):
        one_root(lambda x: x * x + 1.0, -1.0, 1.0, "no root")


def test_short_column_plate_height_and_newton_route():
    """A coarse run: the two-station H is close to eq (38), and the deferred and newton
    routes give the same H."""
    u, D, alpha, K, F_I, F_II = 0.2275, 0.1512, 0.04825, 1.3, 0.4, 0.4   # the page's column
    xK = K * F_I / F_II
    H38 = 2 * D / u + 2 * (u * F_I / alpha) / (1 + xK) ** 2
    col = Column(u, D, alpha, K, F_I, F_II, L=100.0, n_z=100)
    ts, rec, zs = col.run(1500.0, 4.0, 12.0, [20.0, 80.0])
    Kc = [cumulants(ts, rec[:, j]) for j in range(2)]
    H = plate_height(Kc, zs)
    assert abs(H / H38 - 1) < 0.01
    ts2, rec2, zs2 = col.run_newton(1500.0, 4.0, 12.0, [20.0, 80.0])
    Kc2 = [cumulants(ts2, rec2[:, j]) for j in range(2)]
    assert abs(plate_height(Kc2, zs2) / H - 1) < 1e-3


def test_stored_claim_reproduced():
    """The committed eq (21) a -> 0 limit as printed, held to its own rtol."""
    stored = load_claims(PAGE)
    assert abs(width_red(1e-8, "delta") / stored["eq21_printed_a0_limit"] - 1) <= 0.05
