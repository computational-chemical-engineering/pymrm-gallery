"""Fast checks of model.py for page B1.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import EXPOSED, spherical_shell, sphere, square_rod  # noqa: E402

PAGE = HERE.name


def test_sphere_against_eq17_second_order():
    """Linear solve at Lambda = 1 against Aris eq. 17, second order in the grid."""
    exact = (3.0 / np.tanh(3.0) - 1) / 3.0
    err = []
    for n in (100, 200):
        pel, va = sphere(n)
        err.append(abs(pel.eta(1.0 / va) / exact - 1))
    assert err[1] < 1e-5
    assert 1.9 <= np.log2(err[0] / err[1]) <= 2.1


def test_flux_route_equals_volume_average():
    """Surface flux and volume average agree to the solve residual on a 2-D grid."""
    pel, va = square_rod(n=40)
    chi = pel.solve(1.0 / va)
    assert abs(pel.eta_flux(1.0 / va, chi) - pel.eta(1.0 / va, chi)) < 1e-10
    assert abs(pel.eta_flux(1.0 / va, chi, read="cell") - pel.eta(1.0 / va, chi)) > 1e-4


def test_shell_inner_face_sealed():
    """An exposed inner face changes the shell's eta (the silent trap of the page)."""
    pel, va = spherical_shell(0.5, n=100)
    bad, _ = spherical_shell(0.5, n=100, bc_low=EXPOSED)
    assert abs(bad.eta(1.0 / va) - pel.eta(1.0 / va)) > 1e-2


def test_stored_claim_reproduced():
    """The committed square_rod_vs_series, held to the rtol of its claim (0.05)."""
    stored = load_claims(PAGE)["square_rod_vs_series"]
    m = 2 * np.arange(800) + 1
    w = 8.0 / (m**2 * np.pi**2)
    b = (m * np.pi / 2) ** 2
    pel, va = square_rod(n=200)
    worst = 0.0
    for L in (0.2, 1.0, 5.0):
        lam = L / va
        ref = 1.0 - float(np.sum(np.outer(w, w) * lam**2 / (lam**2 + b[:, None] + b[None, :])))
        worst = max(worst, abs(pel.eta(lam) / ref - 1))
    assert abs(worst - stored) <= 0.05 * abs(stored)
