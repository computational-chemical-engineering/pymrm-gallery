"""Fast checks of model.py for page I1.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (Channel, Graetz, Voltz, cross, parse_flux, reconstruct_groups,  # noqa: E402
                   wall_flux)

PAGE = HERE.name


def _channel():
    pv = load_data("heck-1976-printed-values.csv", page=PAGE)
    P = {r.quantity: r.value for r in pv.itertuples()}
    f = lambda k: float(P[k])                                   # noqa: E731
    dTad, _, b0 = reconstruct_groups(f("delta_H") * 1e6, f("CG_CO_inlet"),
                                     f("channel_diameter") * 1e-3, f("Nu_used_1D"),
                                     f("TG_inlet_fig1") + 273.0)
    kin = Voltz(f("k_r0"), f("k_a0"), f("Er_over_R"), f("Ea_over_R"), f("CG_CO_inlet"),
                f("CG_O2_inlet"))
    return Channel(kin, dTad, b0), f, P


def test_upper_threshold_reproduces_stored_claim():
    """The 1-D inlet light-off threshold, held to its claim's rtol (1e-5)."""
    ch, f, _ = _channel()
    T = brentq(lambda T: ch.fold_CG(T, ch.b0 * f("Nu_used_1D")) - ch.cg0, 280.0, 420.0,
               xtol=1e-8)
    T0 = load_claims(PAGE)["T_up_1d_square"]
    assert abs(T - T0) <= 1e-5 * abs(T0)


def test_lightoff_root_inside_its_step():
    """The refined 2-D light-off lies inside the step the unrefined march reports."""
    ch, _, _ = _channel()
    g0 = ch.march_2d(316.0, 1e-3, 0.4, n_r=16, refine=False)
    g1 = ch.march_2d(316.0, 1e-3, 0.4, n_r=16)
    assert g0 < g1 <= g0 + 1e-3


def test_newton_and_deferred_correction_agree():
    """ConvFluxTVD with newton solves the same van Leer equations as the deferred correction."""
    _, _, P = _channel()
    Q = wall_flux(*parse_flux(P["Qw_fig9"]))
    g = Graetz(X=0.05, n_x=200, n_r=10, r_grade=1.0)
    a, b = g.solve(("flux", Q)), g.solve_newton(("flux", Q))
    assert np.max(np.abs(a["T"] - b["T"])) < 1e-9
    assert abs(cross(a["x"], a["Tw"] - a["cup"]) - cross(b["x"], b["Tw"] - b["cup"])) < 1e-9
