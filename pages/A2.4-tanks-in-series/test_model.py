"""Fast checks of model.py for page A2.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import Cascade, dispersion_exit  # noqa: E402

PAGE = HERE.name


def test_cascade_is_eq9():
    """N stirred tanks with the library outflow bc give Levenspiel's Eq. 9."""
    c_exit, _ = Cascade(4, tbar=1.0, c_in=1.0).first_order(2.0)
    assert abs(c_exit - (1.0 + 2.0 / 4) ** -4) < 1e-12


def test_outflow_equals_sink_construction():
    """{"outflow": True} is the page's own sink construction, entry for entry."""
    lib, own = Cascade(6, c_in=1.0), Cascade(6, c_in=1.0, outlet="stirred")
    assert abs(lib.jac_const - own.jac_const).max() <= 1e-14
    assert np.max(np.abs(lib.g_const - own.g_const)) <= 1e-14


def test_stored_claim_reproduced():
    """The committed disp_reldiff_N16_kt2, held to the rtol of its claim (0.05)."""
    stored = load_claims(PAGE)["disp_reldiff_N16_kt2"]
    tis = (1.0 + 2.0 / 16) ** -16
    gap = abs(dispersion_exit(32.0, 2.0) - tis) / tis
    assert abs(gap - stored) <= 0.05 * abs(stored)
