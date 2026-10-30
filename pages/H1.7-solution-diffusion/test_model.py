"""Fast checks of model.py for page H1.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import film_analytic, film_wall_ratio, polarised_fluxes  # noqa: E402

PAGE = HERE.name
D, R_INT, LMH = 1.5e-9, 0.994, 1.0 / 3.6e6


def _j_max():
    obs = load_data("wijmans-baker-1995-fig5.csv", page=PAGE)
    return obs[obs.series == "water_flux"].value.max() * LMH


def test_stored_claim_reproduced():
    """The committed film error at n = 400, held to its claim's rtol (1e-3)."""
    Jv = _j_max()
    ex = film_analytic(Jv, 50e-6, D, R_INT)
    err = abs(film_wall_ratio(Jv, 50e-6, D, R_INT, n=400) - ex) / ex
    x0 = load_claims(PAGE)["film_err_n400"]
    assert abs(err - x0) <= 1e-3 * abs(x0)


def test_first_order_and_face_beats_cell():
    """Upwind is first order; the face value is closer than the last cell centre."""
    Jv = _j_max()
    ex = film_analytic(Jv, 100e-6, D, R_INT)
    e = [abs(film_wall_ratio(Jv, 100e-6, D, R_INT, n=n) - ex) for n in (200, 400)]
    assert 0.95 < np.log2(e[0] / e[1]) < 1.05
    assert e[1] < abs(film_wall_ratio(Jv, 100e-6, D, R_INT, n=400, wall="cell") - ex)


def test_polarised_fluxes_solve_their_equation():
    """Each flux satisfies J = A (p - pi w(J)), and no film gives the linear law."""
    p, A, pi = np.array([600.0, 900.0]), 0.05, 350.0
    w = lambda J: film_analytic(J * LMH, 50e-6, D, R_INT)
    J = polarised_fluxes(p, A, pi, w)
    assert np.allclose(J, A * (p - pi * np.array([w(j) for j in J])), rtol=1e-10)
    assert np.allclose(polarised_fluxes(p, A, pi, lambda J: 1.0), A * (p - pi))
