"""Fast checks of model.py for page A4.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_sherwood_at_75_percent_is_the_stored_value():
    """Sh at F = 0.75 of eq. (30), a root of the series, equals the stored sh75_series_exact."""
    fo = model.one_root(lambda f: float(model.series_one_minus_F(f)[0]) - 0.25, 1e-3, 1.0,
                        "F = 0.75", xtol=1e-15)
    sh = float(model.sherwood(fo, 0.25))
    stored = load_claims(PAGE)["sh75_series_exact"]
    assert abs(sh - stored) <= 1e-6 * stored


def test_isotherm_inverse_and_single_site_collapse():
    """P_of inverts eq. (32); with one site eq. (33) is 1/(1 - theta)."""
    iso = model.Isotherm(0.4, 4.0, 3e-4, 4.2)
    P = np.logspace(-4, 8, 13)
    assert np.max(np.abs(iso.P_of(iso.theta(P)) / P - 1)) < 1e-10
    one = model.Isotherm(4.86e-6, 11.0)
    th = np.linspace(0.1, 10.0, 7)
    assert np.max(np.abs(one.Gamma(th) * (1 - th / 11.0) - 1)) < 1e-12


def test_constant_d_crystal_against_series():
    """A constant-D sphere (one run, 100 cells) against eq. (30) within its step error."""
    c = model.Crystal(1.0, n_r=100)
    F = c.run(0.0, 1.0, [0.05], s=0.25)[0][0]
    assert abs(F / model.series_F(0.05)[0] - 1) < 1e-3


def test_loading_dependent_crystal_one_newton_run():
    """A loading-dependent sphere solves (Newton reports success on every step) and stays bounded."""
    c = model.Crystal(lambda th: 1.0 / (1.0 - th), n_r=40)
    m, prof = c.run(0.0, 0.6, [0.02])
    assert 0.0 < m[0] < 0.6
    assert np.all((prof[0.02] > -1e-12) & (prof[0.02] < 0.6 + 1e-12))
