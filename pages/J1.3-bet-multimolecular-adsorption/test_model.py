"""Fast checks of model.py for page J1.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_ladder_matches_eqB():
    """The layer equilibria solved by newton and summed term by term give eq. (B)."""
    for n in (1, 6, 20):
        v, sol = model.ladder_v_over_vm(0.58, 156.7, n)
        assert sol.success
        assert abs(v / float(model.v_eqB(0.58, 156.7, n)) - 1) < 1e-12


def test_dipole_ladder_matches_eq1b():
    """At k below 1/2 the solved recursion decays at the ratio of eq. (1b)."""
    k = 0.0102
    assert abs(model.dipole_ladder(k, 200)[0] / model.C_eq1b(k) - 1) < 1e-12


def test_ceiling_matches_stored_value():
    """The root-found inflection ceiling is 2/sqrt(3) and equals the stored claim."""
    c_star = model.c_star_root(1e-5)
    ceiling = model.inflection_ratio(c_star)
    assert abs(ceiling / (2 / np.sqrt(3)) - 1) < 1e-12
    assert abs(c_star / (27 + 15 * np.sqrt(3)) - 1) < 1e-9
    stored = load_claims(PAGE)["bet_inflection_ceiling"]
    assert abs(ceiling - stored) <= 1e-6 * stored
