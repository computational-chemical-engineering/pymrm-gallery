"""Fast checks of model.py for page H1.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import LEVEQUE_AVG, TubeMarch, cpc_zero_flux, film_cpc, graetz_lambda0  # noqa: E402

PAGE = HERE.name
THETA, GZ, DA = 1.01 / 20.0, 500.0, 25.0
BETA = 2.0 * DA / (LEVEQUE_AVG * GZ ** (1.0 / 3.0))


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_film_and_constants():
    """The exact film, its zero-flux limit and the Graetz eigenvalue meet their check values."""
    for name, x in (("TC_film_exact_CPC", film_cpc(0.5, BETA, THETA)),
                    ("TC_film_printed_CPC", film_cpc(0.5, BETA, THETA, "printed")),
                    ("TC_CPC_zero_flux_limit", cpc_zero_flux(BETA, THETA)),
                    ("graetz_lambda0", graetz_lambda0())):
        ref, rtol = check_value(name)
        assert abs(x / ref - 1) < rtol
    assert abs(film_cpc(0.5, BETA, THETA, "endform") / film_cpc(0.5, BETA, THETA) - 1) < 1e-10


def test_one_march():
    """A coarse 2-D march at TC: N2 conserved to round-off, the outlet CPC within its
    coarse-grid error of the check value and of the stored claim."""
    h = TubeMarch(0.5, THETA, GZ, DA, n=50, k=100).run()
    assert h.n2_balance() < 1e-12
    ref, _ = check_value("TC_2D_CPC_out")
    assert abs(h.cpc_local()[-1] / ref - 1) < 2e-4
    stored = load_claims(PAGE)
    assert abs(h.cpc_local()[-1] / stored["r1_cpc_2d"] - 1) < 2e-4
