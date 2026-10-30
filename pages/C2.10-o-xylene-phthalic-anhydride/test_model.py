"""Fast checks of model.py for page C2.10 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (FromentCase, Tube1D, Tube1DSteady, Tube2D, exchanger_mean,  # noqa: E402
                   froment_source)

PAGE = HERE.name


def _case():
    par = load_data("froment-1967-parameters.csv", page=PAGE)
    return FromentCase.from_table(zip(par.symbol, par.value))


def test_radial_operator_against_series():
    """Reactionless exchanger: pymrm radial operator against the paper's Bessel series."""
    from scipy.integrate import solve_ivp
    case = _case()
    m = Tube2D(lambda u: np.zeros_like(u), (case.a2,), (case.alpha,), case.R, n_r=40)
    sol = solve_ivp(m.rhs, (0.0, 50.0), np.ones(40), method="BDF", jac=m.jac_disp,
                    rtol=1e-11, atol=1e-13)
    ref = exchanger_mean(50.0, case.a2, case.alpha, case.R)[0]
    assert abs(m.weights @ sol.y[:, -1] - ref) / ref < 5e-4


def test_steady_1d_matches_march():
    """The steady upwind route against the march, at the first-order upwind gap."""
    case = _case()
    src = froment_source(case, 357.0)
    wall = (0.0, 0.0, case.wall_sink(82.7))
    marched = Tube1D(src, wall).march(case.z_end, cap=300.0,
                                      max_step=case.z_end / 200).hot_spot()[1]
    b = Tube1DSteady(src, wall, case.z_end, n_z=2000)
    b.solve(np.zeros(b.shape), scale=np.array([1.0, 1.0, 100.0]))
    assert abs(b.hot_spot()[1] - marched) / marched < 5e-3


def test_stored_claim_reproduced():
    """The committed 2-D hot spot at 357 degC, held to the rtol of its claim (0.05)."""
    case = _case()
    m = Tube2D(froment_source(case, 357.0), (case.a1, case.a1, case.a2),
               (0.0, 0.0, case.alpha), case.R, n_r=80).march(case.z_end, cap=150.0)
    x0 = load_claims(PAGE)["tm2d_hotspot_357_C"]
    assert abs(m.hot_spot()[1] - x0) <= 0.05 * abs(x0)
