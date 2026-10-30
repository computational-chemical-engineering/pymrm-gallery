"""Fast checks of model.py for page G1.8 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import PAPER, WettedSphere, chi_collapsed, chi_table1, eta_sphere, reconstruct_L  # noqa: E402

PAGE = HERE.name


def test_fv_pellet_matches_table1_and_stored_value():
    """One Robin solve against Table 1's sphere row; the stored n_u = 800 claim to its rtol."""
    s = WettedSphere(n_u=800)
    a = float(PAPER.alpha_gLs(2.0))
    phis = np.logspace(np.log10(1.07), np.log10(45.5), 60)
    num = np.array([s.eta(s.solve(p, a)) for p in phis])
    dev = np.max(np.abs(num * (phis ** 2 / a + 1.0 / eta_sphere(phis)) - 1.0))
    x0 = load_claims(PAGE)["pymrm_vs_table1_L2"]
    assert abs(dev - x0) <= 0.05 * abs(x0)


def test_surface_flux_identity():
    """The volume-average eta equals the one from the boundary-face gradient."""
    s = WettedSphere(n_u=200)
    for alpha in (float(PAPER.alpha_gLs(0.5)), np.inf):
        y = s.solve(10.0, alpha)
        assert abs(s.eta_surface(y, 10.0, alpha) / s.eta(y) - 1.0) < 1e-10


def test_reconstruction_stored_values():
    """The reconstructed L_m of the four drawn curves, to the claims' rtol (1e-6)."""
    chi10 = load_data("herskowitz-smith-1983-fig6.csv", page=PAGE)["chi_at_phi_10"].to_numpy()
    stored = load_claims(PAGE)
    for i, c in enumerate(chi10):
        x0 = stored[f"L_rec_curve{i+1}"]
        assert abs(reconstruct_L(c) - x0) <= 1e-6 * abs(x0)


def test_collapse_identity():
    """chi_collapsed is chi_table1 at alpha_gs -> inf, C*_L = 1."""
    phi = np.logspace(0, 1.7, 9)
    assert np.max(np.abs(chi_collapsed(phi, 2.0) / chi_table1(phi, 2.0) - 1.0)) < 1e-12
