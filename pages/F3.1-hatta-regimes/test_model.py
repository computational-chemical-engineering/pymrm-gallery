"""Fast checks of model.py for page F3.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_first_order_film_matches_closed_form():
    """The linear film reproduces Ha/tanh(Ha) on 400 cells at Ha = 5."""
    E = model.Film(400, 1).solve_first_order(5.0)[1]
    assert abs(E / (5.0 / np.tanh(5.0)) - 1.0) < 1e-6


def test_second_order_excess_and_stored_limit():
    """With B in excess the second-order film is the first-order one; the small-Ha
    end of the sweep matches the stored claim."""
    f2 = model.Film(400, 2)
    E2 = f2.solve_second_order(2.0, 1e6)[1]
    E1 = model.Film(400, 1).solve_first_order(2.0)[1]
    assert abs(E2 / E1 - 1.0) < 1e-5
    lo = f2.solve_second_order(0.05, 20.0)[1]
    stored = load_claims(PAGE)["limit_E_small_Ha"]
    assert abs(lo - stored) <= 1e-9 * stored


def test_vkh_reduces_to_first_order():
    """Van Krevelen-Hoftijzer is Ha/tanh(Ha) when E_i is very large."""
    assert abs(model.vkh(3.0, 1e9) / (3.0 / np.tanh(3.0)) - 1.0) < 1e-8
