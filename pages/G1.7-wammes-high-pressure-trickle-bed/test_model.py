"""Fast checks of model.py for page G1.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (CO2PlugFlow, HoldupPressureClosure, OperatingPoints,  # noqa: E402
                   conversion_exact, gas_viscosity)

PAGE = HERE.name


def _params():
    par = load_data("wammes-1991-parameters.csv", page=PAGE)
    return dict(zip(par.quantity, par.value.astype(float)))


def _fig7():
    P = _params()
    h = load_data("wammes-1991-fig7.csv", page=PAGE)
    eta = gas_viscosity(h.rho_g_kg_m3.values, P["n2_density_per_bar"],
                        P["n2_viscosity_1bar"], P["n2_viscosity_100bar"])
    op = OperatingPoints(h.v_g_m_s.values, h.v_l_m_s.values, h.rho_g_kg_m3.values, eta,
                         P["water_density"], P["water_viscosity"],
                         P["cylinder_bed_porosity"], P["cylinder_dp"], P["eps_beta_residual"])
    return h, op


def test_stored_claim_reproduced():
    """The committed Fig. 7 mean absolute deviation, held to its claim's rtol (0.05)."""
    h, op = _fig7()
    s = HoldupPressureClosure(op).solve()
    dev = np.abs(s["beta_dyn"] / h.beta_dyn.values - 1.0).mean()
    x0 = load_claims(PAGE)["holdup_mean_abs_dev"]
    assert abs(dev - x0) <= 0.05 * abs(x0)


def test_newton_matches_brentq():
    """The vector Newton solve and the scalar root-finder agree point by point."""
    _, op = _fig7()
    cl = HoldupPressureClosure(op)
    assert np.abs(cl.solve()["beta_dyn"] - cl.solve_scalar()["beta_dyn"]).max() < 1e-9
    assert np.all(cl.sign_changes(200) == 1)


def test_absorber_first_order_and_stored_value():
    """Exit-face conversion converges at first order onto Eq. 2 integrated."""
    P = _params()
    a_gl = 0.5 * 6.0 * (1.0 - P["cylinder_bed_porosity"]) / P["cylinder_dp"]
    args = (a_gl, P["K_2M_DEA"], P["bed_height_absorption"], 0.25)
    exact = conversion_exact(*args)
    errs = []
    for n in (100, 200):
        m = CO2PlugFlow(*args, n_x=n)
        m.solve()
        errs.append(abs(m.conversion() - exact))
    assert 0.95 < np.log2(errs[0] / errs[1]) < 1.05
    x0 = load_claims(PAGE)["xi_upwind_err_n200"]
    assert abs(errs[1] - x0) <= 0.05 * abs(x0)
