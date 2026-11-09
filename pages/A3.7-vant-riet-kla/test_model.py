"""Fast checks of model.py for page A3.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (crossing_pv, kla, phi_closed, phi_pymrm,  # noqa: E402
                   stripping_number)

PAGE = HERE.name


def test_branches_cross_once():
    """Figure 5's crossing: the pure-water curve at 4 cm/s meets the ionic curve at
    0.5 cm/s where the two correlations are equal."""
    x = crossing_pv(0.040, 0.005)
    assert abs(kla(x, 0.040, "pure_water") / kla(x, 0.005, "ionic") - 1) < 1e-12


def test_pymrm_column_against_closed_form():
    """The pymrm gas column agrees with (1 - exp(-N_G))/N_G to first order in the grid,
    and its face-flux balance closes to round-off."""
    k = float(kla(np.sqrt(500.0 * 1e4), 0.005, "ionic"))
    H_t = (4.0 * 2.6 / np.pi) ** (1 / 3)
    ng = stripping_number(k, H_t, 0.005, 0.10, 30.0)
    exact = float(phi_closed(ng))
    phi, gap = phi_pymrm(200, H_t, 0.005, 0.10, k, 30.0, return_balance=True)
    assert abs(phi / exact - 1) < 2e-3
    assert gap < 1e-12


def test_stored_metric_reproduced():
    """The committed phi_ionic claim, to rtol 1e-9."""
    stored = load_claims(PAGE)["phi_ionic"]
    k = float(kla(np.sqrt(500.0 * 1e4), 0.005, "ionic"))
    H_t = (4.0 * 2.6 / np.pi) ** (1 / 3)
    ng = float(stripping_number(k, H_t, 0.005, 0.10, 30.0))
    assert abs(float(phi_closed(ng)) - stored) <= 1e-9 * abs(stored)
