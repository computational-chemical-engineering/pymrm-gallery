"""Fast checks of model.py for page B2.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import Bed, PlugExact, _cascade, run_mixed  # noqa: E402

PAGE = HERE.name


def test_frozen_bed_second_order():
    """Frozen activity a = exp(-3 s): newton and the deferred sweeps agree, outlet second order."""
    ref = float(np.exp(-1.5 * (1.0 - np.exp(-3.0)) / 3.0))
    err = []
    for m in (100, 200):
        bed = Bed(m)
        c = bed.steady_c(1.5, np.exp(-3.0 * bed.s_c))
        c_d = Bed(m, route="deferred").steady_c(1.5, np.exp(-3.0 * bed.s_c))
        assert np.abs(c - c_d).max() < 1e-11
        err.append(abs(bed.outlet(c) - ref))
    assert 1.9 <= np.log2(err[0] / err[1]) <= 2.1


def test_short_march_against_quadrature():
    """Parallel bed, 100 cells, 40 steps, against the quadrature route."""
    X, _, _ = Bed(100).march("parallel", 1.5, 3.0, 1.0, 1.0, nt=40)
    X_ref = PlugExact(41).run("parallel", 1.5, 3.0, 1.0, 1.0, rtol=1e-10)[0]
    assert np.abs(X - X_ref).max() < 1e-4


def test_stored_claims_reproduced():
    """Two committed claims, each held to its own check."""
    stored = load_claims(PAGE)
    # the outlet activity of the quadrature route against the bed's doubly
    # extrapolated value, within the claim's bound (1e-7)
    a_quad = float(PlugExact(81).run("parallel", 1.5, 3.0, 1.0, 1.0, rtol=1e-12)[1][-1, -1])
    assert abs(a_quad - stored["bed_a_outlet_end_par"]) < 1e-7
    # degeneracy 2 at one order combination, within its bound (1e-8)
    ys = run_mixed("series", 1.5, 3.0, 1.0, 1.0, n=1.0)[0]
    yp = run_mixed("parallel", 1.5, 3.0 * 1.5, 2.0, 1.0, n=1.0)[0]
    assert np.abs(ys - yp).max() <= 1e-8
    assert stored["degen_mixed_series_is_parallel"] <= 1e-8
    assert _cascade(np.array([0.0]), 1.0, 1.0)[0] == 1.0
