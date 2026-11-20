"""Fast checks of model.py for page C2.13 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (IntegralReactor, co_rate_maximum, co_rate_maximum_closed_form,  # noqa: E402
                   march)

PAGE = HERE.name


def _par():
    t = load_data("voltz-1973-table1.csv", page=PAGE)
    return dict(zip(t.symbol, t.value.astype(float)))


def test_reactor_against_march():
    """Curve A of Figure 19: the pymrm reactor meets the solve_ivp march at first order."""
    par = _par()
    inlet = (4.0, 800.0, 4.5, 0.0)
    gaps = [abs(IntegralReactor(par, 550.0, inlet, n=n).solve(0.08)[0]
                - march(par, 550.0, inlet, 0.08)[0]) for n in (200, 400)]
    assert gaps[1] < 0.05
    assert 1.8 < gaps[0] / gaps[1] < 2.2


def test_closed_form_maximum():
    """Without k_a3 the root of d ln r/dCO is (1 + k_a2 HC)/k_a1."""
    par = dict(_par(), k0_a3=0.0)
    assert abs(co_rate_maximum(par, 550.0, 50.0) / co_rate_maximum_closed_form(par, 550.0, 50.0)
               - 1.0) < 1e-9


def test_stored_claim_reproduced():
    """The committed CO of the rate maximum at 550 degF, held to its claim's rtol (0.05)."""
    x0 = load_claims(PAGE)["co_max_550F"]
    assert np.isclose(co_rate_maximum(_par(), 550.0, 50.0), x0, rtol=0.05)
