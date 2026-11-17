"""Fast checks of model.py for page F3.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_film_error_against_collocation_and_stored_value():
    """At (gamma, q) = (10, 4) one graded n_x = 800 solve differs from collocation
    by the stored fv_rel_error_n800 (its rtol, 5 %)."""
    f_fv = model.Film(n_x=800, stretch=8.0).f_a(10.0, 4.0)[0]   # raises if Newton fails
    f_col = model.f_a_bvp(10.0, 4.0, tol=1e-11)[0]
    err = abs(f_fv / f_col - 1.0)
    stored = load_claims(PAGE)["fv_rel_error_n800"]
    assert abs(err - stored) <= 0.05 * stored


def test_first_integral_and_vkh_routes():
    """q beta(0) = 1 + q - F_A on the film, and the two VKH routes agree."""
    f_a, b0 = model.Film(n_x=400, stretch=8.0).f_a(4.7, 1.5)
    assert abs(1.5 * b0 - (2.5 - f_a)) / 2.5 < 1e-8
    assert abs(model.f_a_vkh(4.7, 1.5) / model.f_a_vkh_route_b(4.7, 1.5) - 1.0) < 1e-12
    assert np.isnan(model.unique_root_log(lambda g: 1.0 + g, 0.5, 1000.0))
