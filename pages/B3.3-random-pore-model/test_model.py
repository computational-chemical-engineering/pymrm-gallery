"""Fast checks of model.py for page B3.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (X_of_tau, X_pet_via_eps, X_petersen, XM_eq33, dSS0_dX, eps_o_opt_block,  # noqa: E402
                   m_star, march37, stat35, tau_at_X_block, unique_root)

PAGE = HERE.name


def test_closed_forms():
    """Eq. (33) is the root of dS/dX, and Petersen's eq. (6) is the porosity route."""
    for pv in (3.0, 10.0):
        assert abs(unique_root(lambda X: dSS0_dX(X, pv), 1e-12, 1 - 1e-12) - XM_eq33(pv)) < 1e-14
    s = np.linspace(0.0, 1.5, 11)
    assert np.abs(X_petersen(s, 0.3) - X_pet_via_eps(s, 0.3)).max() < 1e-14


def test_pymrm_blocks():
    """The tau block inverts eq. (32), the eq. (35) block is a root, the march converges."""
    psis = np.array([0.0, 2.0, 100.0])
    tau = tau_at_X_block(0.5, psis, 10.0)
    assert np.abs(X_of_tau(tau, psis, 10.0) - 0.5).max() < 1e-12
    e0 = eps_o_opt_block(0.3, 0.67)
    assert abs(stat35(e0, 0.3, 0.67)[0]) < 1e-12
    X = march37([0.0, 5.0], 1.0, 100)
    assert np.abs(X - X_of_tau(1.0, np.array([0.0, 5.0]))).max() < 1e-4


def test_stored_claims_reproduced():
    """One committed claim, held to its own rtol."""
    stored = load_claims(PAGE)
    assert abs(m_star(1.0) / stored["m_star_psi1"] - 1) <= 1e-6
