"""Fast checks of model.py for page B3.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import E_exact, G_of_xi, GrainModel, X_of_theta, lag, xi_of_zeta, yagi_kunii_eq6  # noqa: E402

PAGE = HERE.name


def test_exact_route():
    """G inverts, E^2 = 6 (theta - Delta) beyond theta*, and G/theta* is eq. 6 without film."""
    g = 0.5
    z = np.linspace(0.0, 1.0 + g / 6, 41)
    assert np.abs(G_of_xi(xi_of_zeta(z, g), g) - z).max() < 1e-14
    th = np.array([1.2, 5.0, 20.0])
    assert np.abs(E_exact(th, g) ** 2 / 6 - (th - lag(g))).max() < 1e-12
    assert abs(lag(g) - (0.25 + g / 30)) < 1e-15
    assert abs(float(X_of_theta(1e-4, g)) / (1.5e-8) - 1) < 1e-3          # X ~ 3 theta^2 / 2
    xi = np.linspace(0, 1, 101)
    assert np.abs(G_of_xi(xi, g) / (1 + g / 6) - yagi_kunii_eq6(xi, g / 6, 1e12)).max() < 1e-11


def test_march_exact_and_numjac():
    """A short march agrees with the closed form, both Jacobians give the same state."""
    a = GrainModel(n=150, L=30.0, g=0.5).march(1.5, 0.01)
    b = GrainModel(n=150, L=30.0, g=0.5, jacobian="numjac").march(1.5, 0.01)
    assert abs(a.EP() / float(E_exact(1.5, 0.5)) - 1) < 3e-4
    assert np.abs(a.xi - b.xi).max() < 1e-12 and np.abs(a.psi - b.psi).max() < 1e-12
    assert a.balance_gap() < 1e-12


def test_stored_claims_reproduced():
    """One committed claim, held to its own rtol."""
    stored = load_claims(PAGE)
    mm = GrainModel(n=300, L=30.0, g=0.5).march(2.0, 0.005)
    gap = abs(mm.I_gas - mm.solid_consumed()) / mm.solid_consumed()
    assert abs(gap / stored["mass_balance_rel_dt5e3"] - 1) <= 0.05


def test_jacobian_gap_detects_defect():
    """The analytic Jacobian agrees with NumJac; dropping d(sink)/d(xi) breaks that."""
    m = GrainModel(n=60, L=20.0, g=0.5).march(1.05, 0.01)
    assert m.jacobian_gap(0.01) < 1e-5
    bad = GrainModel(n=60, L=20.0, g=0.5, jac_defect="drop_dsink_dxi").march(1.05, 0.01)
    assert bad.jacobian_gap(0.01) > 1e-3
