"""Fast checks of model.py for page H1.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (G_half, Permeator, alpha_from_wall, exact_half,  # noqa: E402
                   q_series)

PAGE = HERE.name


def _stated():
    vals = load_data("itoh-1987-permeation.csv", page=PAGE)
    return dict(zip(vals.quantity, vals.value.astype(float)))


def test_alpha_h_stored_claim():
    """The wall solve reproduces the committed deviation from the printed alpha_H."""
    s = _stated()
    a = alpha_from_wall(s["r_i"], s["r_o"], s["D_H_in_Pd"], s["C_0"], s["l_0"])
    dev = (a - s["alpha_H"]) / s["alpha_H"] * 100
    x0 = load_claims(PAGE)["alpha_H_dev_vs_printed_pct"]
    assert abs(dev - x0) <= 1e-3 * abs(x0)
    closed = 2 * np.pi * s["l_0"] * s["D_H_in_Pd"] * s["C_0"] / np.log(s["r_o"] / s["r_i"])
    assert abs(a - closed) / closed < 1e-8


def _beta(s):
    n_ar = (1 - s["feed_cyclohexane_fraction"]) / s["feed_cyclohexane_fraction"]
    return (n_ar + 1) / 3.0


def test_channel_page_route():
    """The page's solve in s from the flat state, below and above Pi_crit, against
    the exact solution and the march."""
    s = _stated()
    beta = _beta(s)
    pi_crit = G_half(1.0, beta)
    for ratio in (0.9, 10.0):
        p = Permeator(ratio * pi_crit, beta, n_z=200)
        p.solve()
        assert p.w.min() >= 0.0
        assert np.max(np.abs(p.w - p.march_guess())) < 1e-9
        assert np.max(np.abs(p.w.ravel() - exact_half(p.z_c, beta, ratio * pi_crit))) < 0.05
    assert p.outlet_value() < 1e-20


def test_series_limits():
    """The surface-step series model recovers Sieverts at lambda -> 0."""
    assert abs(q_series(0.5, 1e-12) - np.sqrt(0.5)) < 1e-5
