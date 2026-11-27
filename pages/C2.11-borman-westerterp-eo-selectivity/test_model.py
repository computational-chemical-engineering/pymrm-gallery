"""Fast checks of model.py for page C2.11 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (BertyReactor, air_feed, element_balances, invert_berty,  # noqa: E402
                   outlet_fractions, rates, selectivity, table6_params)

PAGE = HERE.name
P_FIG = np.array([0.09, 1.06, 0.07, 0.06, 0.07])     # Figure 6 caption, bar


def _par():
    return table6_params(load_data("borman-1995-table6.csv", page=PAGE, dtype=str), "III")


def _check(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_check_value_sigma_printed():
    """Selectivity at the Figure 6 condition with Table 6 as printed, against the second route."""
    v, rtol = _check("sigma_T_F1")
    assert abs(selectivity(*rates(_par(), 476.0, P_FIG)) / v - 1.0) <= rtol


def test_element_balances_close():
    """Table 8 conserves C, H and O at a conversion of the size of Table 4."""
    xI, _ = air_feed(0.018)
    gaps = element_balances(xI, outlet_fractions(xI, 0.016, 0.005))
    assert max(abs(g) for g in gaps.values()) < 1e-13


def test_forward_round_trip():
    """The pymrm Berty solve returns the conversions of the explicit inversion (row 2)."""
    par = _par()
    xI, M = air_feed(0.018)
    phi, W, R1, R2 = 0.74 / M, 32.5e-3, 160e-6, 131e-6
    zE, zC = invert_berty(R1, R2, W, phi, xI)
    P = np.array([outlet_fractions(xI, zE, zC)[j] * 5.0 for j in ("E", "O", "C", "EO", "W")])
    m1, m2 = rates(par, 500.0, P)
    rx = BertyReactor(par, 500.0, 5.0, xI, phi, W, 1.0, R1 / m1, R2 / m2, scale=(zE, zC))
    z1, z2, result, check = rx.solve()
    assert result.success and check["ok"]
    assert abs(z1 / zE - 1.0) < 1e-9 and abs(z2 / zC - 1.0) < 1e-9


def test_stored_claim_reproduced():
    """The committed selectivity at the Figure 6 condition (F = 1), to its claim's rtol."""
    x0 = load_claims(PAGE)["sigma_T_F1_pct"]
    assert np.isclose(100.0 * selectivity(*rates(_par(), 476.0, P_FIG)), x0, rtol=0.005)
