"""Fast checks of model.py for page J1.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (Particle, crank_profile, crank_series, k_effective,  # noqa: E402
                   k_effective_crossing, ldf_constant, ldf_error_extremes)

PAGE = HERE.name


def test_series_tail_and_profile_mean():
    """The truncation meets its tail bound and the profile's volume mean is 1 - deficit."""
    tau = 0.02
    s = crank_series(tau, rtol=1e-15)
    assert s["tail"] <= 1e-15
    rho = np.linspace(0.0, 1.0, 40001)
    q = crank_profile(rho, tau, s["n_terms"])
    assert np.trapezoid(3 * rho**2 * q, rho) == pytest.approx(1.0 - s["deficit"][0], abs=1e-7)
    assert ldf_constant(2) == 15 and ldf_constant(1) == 8 and ldf_constant(0) == 3


def test_march_matches_series_and_balances_flux():
    """One backward-Euler march: close to the series, and the flux balance closes."""
    p = Particle(n_r=100)
    taus = np.array([0.02, 0.05])
    res = p.march(taus, 1e-4, keep=True)
    exact = 1.0 - crank_series(taus)["deficit"]
    assert np.all(np.abs(res["mean"] / exact - 1.0) < 1e-3)
    means = p.mean(res["states"])
    cum = (p.nu + 1) * 1e-4 * np.cumsum(p.surface_flux(res["states"]))
    assert np.max(np.abs(means - cum)) < 1e-12
    assert abs(p.w.sum() - 1.0) < 1e-14


def test_ldf_roots():
    """k_eff = 15 at one time, the LDF error has two extrema, and a double root raises."""
    t = k_effective_crossing(15.0)
    assert k_effective(t)[0] == pytest.approx(15.0, rel=1e-9)
    (_, _), (lo, hi) = ldf_error_extremes(15.0)
    assert lo < 0.0 < hi
    with pytest.raises(RuntimeError):
        k_effective_crossing(5.0)      # k_eff never falls below pi^2: no root


def test_stored_claim_reproduced():
    """The committed LDF crossing time, held to its own rtol."""
    stored = load_claims(PAGE)
    got = k_effective_crossing(15.0)
    assert abs(got / stored["k_eff_equals_15_tau"] - 1) <= 0.05
