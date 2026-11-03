"""Fast checks of model.py for page E1.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import DavidsonBubble, cloud_radius_exact, rise_velocity  # noqa: E402

PAGE = HERE.name

# Appendix C's bed of Kunii and Levenspiel (1968), CGS.
D_B, U_MF, EPS_MF = 3.7, 2.1, 0.50
R, U_F, U_BR = D_B / 2, U_MF / EPS_MF, rise_velocity(D_B)


def bubble(n=200, **kw):
    return DavidsonBubble(R, U_F, U_BR, n=n, **kw).solve()


def test_second_order_and_face_reading():
    """R_c and the through-flow converge at second order; the bubble-face gradient
    from compute_boundary_values equals the operator's face gradient."""
    exact = R * cloud_radius_exact(U_F, U_BR)
    err_rc, err_q = [], []
    for n in (100, 200):
        s = bubble(n)
        assert s.dfdr_R == pytest.approx(s.df[0], rel=1e-12)
        err_rc.append(abs(s.cloud_radius() - exact))
        err_q.append(abs(s.throughflow_per_volume(EPS_MF) * D_B / U_MF - 4.5))
    assert 3.6 < err_rc[0] / err_rc[1] < 4.4
    assert 3.6 < err_q[0] / err_q[1] < 4.4


def test_no_cloud_and_cloud_outside_domain():
    """No zero of W below the threshold gives nan; beyond r_out it raises."""
    assert np.isnan(DavidsonBubble(1.0, 1.0, 0.9, n=200).solve().cloud_radius())
    with pytest.raises(RuntimeError):
        DavidsonBubble(1.0, 1.0, 1.0001, n=200).solve().cloud_radius()


def test_stored_claim_reproduced():
    """The committed through-flow error at n = 1600, held to its own rtol."""
    stored = load_claims(PAGE)
    s = bubble(1600)
    got = abs(s.throughflow_per_volume(EPS_MF) * D_B / U_MF - 4.5) / 4.5
    assert abs(got / stored["throughflow_coefficient_rel"] - 1) <= 0.05
