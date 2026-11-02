"""Fast checks of model.py for page B2.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import Bed, ce_exact, ce_patch, pe_exact  # noqa: E402

PAGE = HERE.name


def test_fresh_bed_second_order():
    """Fresh bed against y = exp(-x): newton and the deferred sweeps agree, second order."""
    err = []
    for n in (100, 200):
        bed = Bed("PE", n=n)
        y = bed.solve_y(np.zeros(n))
        y_d = Bed("PE", n=n, route="deferred", tol=1e-13).solve_y(np.zeros(n))
        assert np.abs(y - y_d).max() < 1e-11
        err.append(np.abs(y - np.exp(-bed.x_c)).max())
    assert 1.9 <= np.log2(err[0] / err[1]) <= 2.1


def test_short_march_against_closed_form():
    """Parallel/exponential bed marched to T = 1 on 100 cells, 60 steps."""
    bed = Bed("PE", n=100)
    rec = bed.march(1.0, 60, record_T=[1.0])
    y, C = rec[1.0]
    ye, Ce = pe_exact(bed.x_c, 1.0)
    assert np.abs(y - ye).max() < 1e-4 and np.abs(C - Ce).max() < 1e-4


def test_stored_claims_reproduced():
    """Two committed claims, each held to its own rtol (eq. 49's az reading at a z = 2, T = 5)."""
    stored = load_claims(PAGE)
    bed = Bed("PE", n=200)
    y = bed.solve_y(np.zeros(200))
    assert abs(np.abs(y - np.exp(-bed.x_c)).max() / stored["fresh_bed_maxerr_n200"] - 1) <= 0.05
    y_exact = float(ce_exact([2.0], [5.0])[0][0, 0])
    dev = abs(ce_patch(2.0, 5.0)[0] - y_exact) / y_exact
    assert abs(dev / stored["eq49_relerr_az_reading"] - 1) <= 1e-6
