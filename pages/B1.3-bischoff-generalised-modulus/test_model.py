"""Fast checks of model.py for page B1.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (Power, VolChange, Curve, Slab, eq33, m_eta_n1, m_eta_n2,  # noqa: E402
                   mb_eta)

PAGE = HERE.name


def test_route_b_against_closed_forms():
    """Quadrature of eqs. (14)-(15) against eqs. (21) and (23)-(24)."""
    for closed, n in ((m_eta_n1, 1.0), (m_eta_n2, 2.0)):
        for u in (0.05, 0.5, 0.9):
            m_b, e_b = mb_eta(Power(n), u)
            m_c, e_c = closed(u)
            assert abs(m_b / m_c - 1) < 1e-10 and abs(e_b / e_c - 1) < 1e-10


def test_slab_second_order_against_tanh():
    """First-order slab at m = 1 against tanh(m)/m, second order in the grid."""
    err = [abs(Slab(n).eta(Power(1), 1.0) / np.tanh(1.0) - 1) for n in (100, 200)]
    assert err[1] < 1e-6
    assert 1.9 <= np.log2(err[0] / err[1]) <= 2.1


def test_variable_diffusivity_newton():
    """One Newton solve with d(c) on the faces converges and matches Route B."""
    rf = VolChange(1.0)
    assert abs(Slab(200).eta(rf, 1.0) / Curve(rf).eta(1.0) - 1) < 1e-5


def test_deadzone_regularised_converges():
    """n = 1/2 past m* = 3: the regularised solve converges near eta = 1/m."""
    rf, m = Power(0.5), 5.0
    slab = Slab(100)
    phi = m * np.sqrt(2.0 * rf.Rint(0.0, 1.0))
    c, sol = slab.solve(rf, phi, eps=1e-12)
    assert sol.success
    assert abs(slab.eta_of(c, rf, phi) * m - 1) < 1e-4


def test_stored_claim_reproduced():
    """The committed table1_eq33_three_cells_worst_rel, held to its claim's rtol."""
    stored = load_claims(PAGE)["table1_eq33_three_cells_worst_rel"]
    tab1 = load_data("bischoff-table1.csv", page=PAGE)
    rows = tab1[tab1["omega_co"] != -0.25]
    worst = max(abs(p / eq33(w) - 1) for w, p in zip(rows["omega_co"], rows["eq33_printed"]))
    assert abs(worst - stored) <= 1e-3 * abs(stored)
