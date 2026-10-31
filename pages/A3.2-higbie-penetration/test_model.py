"""Fast checks of model.py for page A3.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import Film, phi_extrapolated, phi_penetration, phi_wall_series  # noqa: E402

PAGE = HERE.name


def test_deep_plug_flow_matches_penetration():
    """A coarse deep plug-flow march reproduces sqrt(4 Lam/pi) after extrapolation."""
    lam = 0.01
    ph = phi_extrapolated(lam, 100, XI=12 * np.sqrt(4 * lam), n_x=100, stretch=4.0, profile="plug")
    assert abs(ph["raw"] / phi_penetration(lam) - 1) < 2e-3
    assert abs(ph["extrapolated"] / phi_penetration(lam) - 1) < 2e-4


def test_mixing_cup_equals_flux_integral():
    """Backward Euler on the conservative scheme: outlet inventory = integrated flux."""
    out = Film(0.05, n_x=60).march(50, record_flux=True)
    assert abs(out["flux_integral"] / out["cup"] - 1) < 1e-11


def test_stored_claims_reproduced():
    """The series wall effect and the profile effect at Lambda = 0.1, by their tolerance."""
    stored = load_claims(PAGE)
    wall = 1 - phi_wall_series(0.1) / phi_penetration(0.1)
    assert abs(wall / stored["wall_effect_series_at_lambda_0p1"] - 1) <= 1e-6
    full = phi_extrapolated(0.1, 200, n_x=150, profile="parabolic")["extrapolated"]
    plug = phi_extrapolated(0.1, 200, n_x=150, profile="plug")["extrapolated"]
    assert abs((1 - full / plug) / stored["profile_effect_at_lambda_0p1"] - 1) <= 0.05
