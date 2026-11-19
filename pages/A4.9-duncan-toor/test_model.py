"""Fast checks of model.py for page A4.9 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name

D_MS = np.zeros((3, 3))
D_MS[0, 1] = D_MS[1, 0] = 83.3e-6
D_MS[0, 2] = D_MS[2, 0] = 68.0e-6
D_MS[1, 2] = D_MS[2, 1] = 16.8e-6
CELL = dict(d_ms=D_MS, length=85.9e-3, la_eff=258.1e2, vol_1=77.99e-6, vol_2=78.63e-6,
            c_t=101.325e3 / (8.314462 * 308.35))
X1 = np.array([0.6, 0.4, 0.0])
X2 = np.array([0.0, 0.5, 0.5])


def test_steady_capillary_flux_is_uniform():
    """The converged capillary solve carries the same flux through every face."""
    cell = model.TwoBulbCell(**CELL, n_z=10)
    assert cell.flux_nonuniformity(X1, X2) < 1e-12


def test_binary_limit_is_fick():
    """With one species absent on both sides the flux is Fick's with the pair diffusivity."""
    cell = model.TwoBulbCell(**CELL, n_z=10)
    xb1, xb2 = np.array([0.7, 0.3, 0.0]), np.array([0.2, 0.8, 0.0])
    n = cell.flux(cell.solve_capillary(xb1, xb2), xb1, xb2)
    fick = cell.c_t * D_MS[0, 1] * (xb1[0] - xb2[0]) / cell.length
    assert abs(n[0, 0] / fick - 1.0) < 1e-10


def test_capillary_volume_fraction_matches_stored_claim():
    """The capillary-to-bulb volume ratio equals the page's stored claim."""
    stored = load_claims(PAGE)["capillary_volume_fraction"]
    value = model.TwoBulbCell(**CELL, n_z=10).capillary_volume_fraction
    assert abs(value - stored) <= 1e-9 * stored
