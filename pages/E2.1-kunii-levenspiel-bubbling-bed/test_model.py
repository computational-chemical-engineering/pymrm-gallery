"""Fast checks of model.py for page E2.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (BubblingBed, Kf_closed, Kf_over_Km, backmix_crossing,  # noqa: E402
                   ceiling_Kf, extended_exact, hydrodynamics)

PAGE = HERE.name
# Appendix C's bed of Kunii and Levenspiel (1968), CGS.
C = dict(u0=13.2, umf=2.1, emf=0.50, em=0.45, Lm=34.0, db=3.7, alpha=0.47, D=0.204)
H = hydrodynamics(**C)


def test_frozen_is_backward_euler_and_limit():
    """Frozen solve equals (1 + K_f/n)^-n; the K_r -> 0 limit is u_0/(u_0 - u_mf)."""
    n, km = 200, 2.0
    s = BubblingBed(H, n, convect=False).solve(km)
    assert s.c_exit == pytest.approx((1 + Kf_closed(km, H) / n) ** -n, rel=1e-12)
    assert Kf_over_Km(0.0, H) == pytest.approx(C["u0"] / (C["u0"] - C["umf"]), rel=1e-14)


def test_balance_and_first_order_to_exact():
    """Outflow outlets close the reactant balance; the error halves with n."""
    ex = extended_exact(2.0, H)["c_exit"]
    errs = []
    for n in (200, 400):
        s = BubblingBed(H, n).solve(2.0)
        assert abs(s.mass_balance_dev()) < 1e-10
        errs.append(abs(s.c_exit - ex))
    assert 1.9 < errs[0] / errs[1] < 2.1
    assert abs(BubblingBed(H, 200, outlet="zero_gradient").solve(2.0).mass_balance_dev()) > 1e-6


def test_crossing_needs_the_range():
    """d_b = 3 cm crosses backmix flow, but not below K_m = 50."""
    h3 = hydrodynamics(**{**C, "db": 3.0})
    assert backmix_crossing(h3, hi=50.0) is None
    assert backmix_crossing(h3) > 50.0


def test_stored_claims_reproduced():
    """The committed ceiling and extension maximum, each held to its own rtol."""
    stored = load_claims(PAGE)
    assert 1 - np.exp(-ceiling_Kf(H)) == pytest.approx(stored["ceiling_conversion_db37"], rel=1e-6)
    base = [np.exp(-Kf_closed(k, H)) for k in (0.5, 1.0, 2.0, 4.0, 8.0)]
    chg = max(extended_exact(k, H)["c_exit"] / b - 1 for k, b in zip((0.5, 1.0, 2.0, 4.0, 8.0), base))
    assert chg == pytest.approx(stored["extension_max_change"], rel=1e-6)
