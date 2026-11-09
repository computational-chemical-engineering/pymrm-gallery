"""Fast checks of model.py for page A3.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (Sh_of_theta, eq12, theta_eq7, theta_eq8, theta_of_Sh,  # noqa: E402
                   theta_pymrm)

PAGE = HERE.name
C = dict(alpha=1.1, beta=0.6, E0=20.0, slope=0.5, pe_pt=2.0, eps_b=0.4, alpha_ranz=0.6,
         beta_ranz=0.5, jd_pt=0.357, jd_exp_pt=-0.359)


def test_solve_converges_to_eq7():
    """First order: the error halves when the cells double, and the outlet of a
    bed without a sink is exactly the inlet value."""
    ex = theta_eq7(5.0, 10.0)
    e1, e2 = theta_pymrm(5.0, 10.0, 200) / ex - 1, theta_pymrm(5.0, 10.0, 400) / ex - 1
    assert 1.9 < e1 / e2 < 2.1
    assert abs(theta_pymrm(0.0, 10.0, 50) - 1.0) < 1e-12


def test_inversion_round_trip():
    """Sh -> outlet -> Sh under eq. (2) and plug flow; eq. (7) tends to eq. (8)."""
    for law in ("eq2", "none"):
        th = theta_of_Sh(12.0, 50.0, 0.6, 10.0, law, C)
        assert abs(Sh_of_theta(th, 50.0, 0.6, 10.0, law, C) / 12.0 - 1) < 1e-9
    assert abs(theta_eq7(2.0, 1e7) / theta_eq8(2.0) - 1) < 1e-6
    assert abs(eq12(1e-30, 0.6, C) - 2.0) < 1e-12


def test_stored_claim_reproduced():
    """The committed pymrm_vs_eq7_n1600, held to the rtol of its claim."""
    stored = load_claims(PAGE)["pymrm_vs_eq7_n1600"]
    err = max(abs(theta_pymrm(N, Pe, 1600) / theta_eq7(N, Pe) - 1)
              for N, Pe in ((1.0, 40.0), (5.0, 10.0), (0.5, 4.0), (20.0, 2.0)))
    assert abs(err - stored) <= 1e-6 * abs(stored)
