"""Fast checks of model.py for page H1.8 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import P_bound, alpha_bound, envelope_fit, excursion  # noqa: E402

PAGE = HERE.name


def _present(pair):
    t = load_data("robeson-2008-table12-upper-bounds.csv", page=PAGE)
    row = t[(t.dataset == "present") & (t.gas_pair == pair)].iloc[0]
    return float(row.k_barrer), float(row.n)


def test_inverse_and_sign_convention():
    """alpha_bound inverts P_bound, and a point above the bound has r > 1."""
    k, n = _present("He/H2")
    a = np.array([0.2, 1.0, 7.5, 900.0])
    assert np.allclose(alpha_bound(P_bound(a, k, n), k, n), a, rtol=1e-12)
    assert np.allclose(excursion(np.array([0.5, 1.0, 2.0]) * P_bound(2.5, k, n), 2.5, k, n),
                       [0.5, 1.0, 2.0], rtol=1e-12)


def test_envelope_covers_and_is_tight():
    """The covering line lies on or above every point and touches at least one."""
    rng = np.random.default_rng(0)
    x = rng.uniform(0.0, 3.0, 15)
    y = 4.0 - 1.5 * x - rng.exponential(0.3, 15)
    k, s, _ = envelope_fit(x, y)
    r = y - (np.log10(k) + s * x)
    assert r.max() <= 1e-9 and r.max() > -1e-9


def test_stored_claim_reproduced():
    """The committed He/H2 envelope span and the count of points above the bound."""
    pts = load_data("robeson-2008-near-bound-points.csv", page=PAGE)
    t = load_data("robeson-2008-table12-upper-bounds.csv", page=PAGE)
    pre = t[t.dataset == "present"].set_index("gas_pair")
    above = sum(int(np.sum(excursion(sub.P_fast_barrer, sub.alpha,
                                     pre.k_barrer[g], pre.n[g]) > 1))
                for g, sub in pts.groupby("gas_pair"))
    stored = load_claims(PAGE)
    assert above == stored["points_above_present_bound"]
    sub = pts[pts.gas_pair == "He/H2"]
    la = np.log10(sub.alpha.to_numpy(float))
    k, s, _ = envelope_fit(la, np.log10(sub.P_fast_barrer.to_numpy(float)))
    x0 = stored["envelope_dn_max_abs"]
    assert abs(abs(s - pre.n["He/H2"]) - x0) <= 1e-6 * abs(x0)
