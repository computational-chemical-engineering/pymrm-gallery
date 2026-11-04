"""Fast checks of model.py for page C1.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (PAR, OptimumCurve, OptimumProfileConverter, alpha_eq, conversion,  # noqa: E402
                   f14, split_quad)

PAGE = HERE.name


def test_closed_form_inversions():
    """The Lambert W and first-order inversions reproduce eq. (14) and its first-order form."""
    aeq, y = 0.93, np.array([0.05, 0.5, 2.0, 8.0])
    al, _ = conversion(y, aeq, "mvk")
    assert np.allclose(f14(al, aeq), y, rtol=1e-12)
    al, _ = conversion(y, aeq, "first_order")
    assert np.allclose(f14(al, aeq, retard=False), y, rtol=1e-12)


def test_alpha_eq_satisfies_eq8():
    """The one root on the physical range satisfies the equilibrium condition."""
    T, s, w = 823.15, 0.08, 0.106
    x = alpha_eq(T, s, w)
    K = 10.0 ** (PAR["slope"] / T - PAR["intercept"])
    assert abs((x / (1 - x)) ** 2 * (1 - s * x / 2) / (w - s * x / 2) / K - 1.0) < 1e-10


def test_split_routes_and_stored_claim():
    """Newton converges, the outflow outlet is the last cell, and the stored split holds."""
    t3 = load_data("mars-vankrevelen-1954-table3-converter.csv", page=PAGE).sort_values(
        "temperature_C")
    group = tuple((T + 273.15, A) for T, A in zip(t3.temperature_C, t3.kstar_group_atm_per_hr))
    curve = OptimumCurve(group, 0.08, 0.106, n_T=300, n_alpha=200)
    split, tau = split_quad(curve)
    conv = OptimumProfileConverter(curve.R, 1.05 * tau, 1000).solve()
    assert conv.result.success and conv.backward_error() < 1e-8
    assert conv.outlet() == conv.alpha[-1, 0]
    pym = 100 * conv.crossing(0.81) / conv.crossing(0.97)
    assert abs(pym - split) < 0.3
    x0 = load_claims(PAGE)["split_converter1_pct_quad"]
    assert abs(split - x0) <= 0.05 * abs(x0)
