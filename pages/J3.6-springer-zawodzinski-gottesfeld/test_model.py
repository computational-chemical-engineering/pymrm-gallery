"""Fast checks of model.py for page J3.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import Membrane, constants, lam_iso, psat  # noqa: E402

PAGE = HERE.name


def table():
    pn = load_data("springer-1991-printed-numbers.csv", page=PAGE, dtype=str)
    C = constants(pn.set_index("item")["printed"].to_dict())
    ps = psat(C["T_cell"], C)
    l2 = lam_iso(C["x_w2"] * C["P_A"] / ps, C)
    l3 = lam_iso(C["x_w3"] * C["P_A"] / ps, C)
    return C, l2, l3


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_interface_water_contents():
    """Eqs. 15 to 17 at the Table I interface fractions meet the second route."""
    _, l2, l3 = table()
    for name, v in (("lambda2_anode", l2), ("lambda3_cathode", l3)):
        ref, rtol = check_value(name)
        assert abs(v / ref - 1) < rtol


def test_one_solve_and_flux_balance():
    """One membrane solve: one flux on every face, alpha near the second route's value."""
    C, l2, l3 = table()
    m = Membrane(400, l2, l3, C["J"], C["T_cell"], C["t_mem"], C).solve(1.8)
    f = m.face_flux(m.lam, m.u)
    assert np.ptp(f) / abs(f[0]) < 1e-9
    ref, _ = check_value("alpha_cfref")
    assert abs(m.alpha() / ref - 1) < 1e-4          # 400 cells, second order, not extrapolated


def test_stored_claim():
    """A stored claim of the page is reproduced: the literal-drag bound 5 lambda_2^2/484."""
    C, l2, _ = table()
    stored = load_claims(PAGE)
    bound = 2 * C["eq18_ndrag_sat"] * l2**2 / C["eq18_lambda_ref"] ** 2
    assert np.isclose(bound, stored["r5_literal_bound"], rtol=1e-9, atol=0)
