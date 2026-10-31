"""Fast checks of model.py for page A2.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (D_two_zone, TwoZoneBed, invert_two_point, simulate_experiment,  # noqa: E402
                   simulate_two_zone)

# Figure 13 conditions, as the page computes them from the printed values
NU_AIR = 0.72 * 0.205
U13 = 24.6 * NU_AIR / (0.300 * 0.383)
D13 = 0.73 * 0.205 + 0.5 * U13 * 0.300 / (1.0 + 9.7 * 0.205 / (U13 * 0.300))
PAGE = HERE.name


def test_upwind_numerical_dispersion_is_predicted():
    """Bare upwind with implicit Euler reports D + u h/2 + u^2 dt/2 on a coarse grid."""
    bed, t, c1, c2, dx = simulate_experiment(U13, D13, n=200, theta=1.0, limiter=None)
    D_hat = invert_two_point(t, c1, c2, dx)[1]
    pred = D13 + U13 * bed.h / 2 + U13 ** 2 * (0.5 * bed.h / U13) / 2
    assert abs(D_hat / pred - 1) < 1e-2


def test_exchange_block_conserves_tracer():
    """The per-cell exchange block has zero column sums: it moves tracer between zones."""
    bed = TwoZoneBed(50.0, 40, 0.15, 5.0, 2.0, 0.15, 0.4, kappa_scale=1.0)
    exch = bed.A - TwoZoneBed(50.0, 40, 0.15, 5.0, 2.0, 0.15, 0.0).A
    assert np.max(np.abs(exch.sum(axis=0))) < 1e-14
    assert exch.nnz == 4 * 40


def test_stored_claim_reproduced():
    """The committed two_zone_rel_err_n500, held to the rtol of its claim."""
    stored = load_claims(PAGE)["two_zone_rel_err_n500"]
    u = 0.0710 * NU_AIR / (0.0097 * 0.41)
    D0 = 0.73 * 0.205 + 0.5 * u * 0.0097 / (1.0 + 9.7 * 0.205 / (u * 0.0097))
    u1 = 2.0 * u
    u2 = (u - 0.15 * u1) / 0.85
    kappa = D_two_zone(0.15, u1, u2, 1.0) / (0.503 - D0)
    _, t, a, b, dx, Dg = simulate_two_zone(0.15, u1, u2, D0, kappa, n=500)
    err = abs(invert_two_point(t, a, b, dx)[1] / Dg - 1)
    assert abs(err - stored) <= 0.05 * abs(stored)
