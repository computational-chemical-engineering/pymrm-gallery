"""Fast checks of model.py for page D3.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import CooledTube, Kinetics, inlet, table6_model3  # noqa: E402
from pymrm.checks import find_roots  # noqa: E402

PAGE = HERE.name
C211 = "C2.11-borman-westerterp-eo-selectivity"
CP = 0.25 * 4184.0                       # J/(kg K), C2.10's c_p (page choice)


def _kin(F=None):
    par = table6_model3(load_data("borman-1995-table6.csv", page=C211, dtype=str))
    F = load_claims(C211)["F_fig6"] if F is None else F
    return Kinetics(par, F=F, P_total=5.0)


def _check(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_check_value_inlet_heat():
    """Inlet heat release at B with F_fig6 against the second route."""
    v, rtol = _check("q_in_B_F6")
    assert abs(_kin().heat(500.0, inlet(0.018)) / v - 1.0) <= rtol


def test_element_balances_on_a_tube():
    """C, H and O are kept cell by cell by the first-order upwind tube."""
    tb = CooledTube(_kin(), 0.018, 500.0, 19.56, CP, 600.0, 200)
    tb.trace("h", np.geomspace(80.0, 19.56, 4))
    f = tb.faces
    for inv in (2 * f[:, 0] + f[:, 2] + 2 * f[:, 3], 4 * f[:, 0] + 4 * f[:, 3] + 2 * f[:, 4],
                2 * f[:, 1] + 2 * f[:, 2] + f[:, 3] + f[:, 4]):
        assert np.max(np.abs(inv - inv[0])) / abs(inv[0]) < 1e-12


def test_frozen_limit_root():
    """A long frozen-composition tube ends at the lower root of q(T) = h (T - T_c)."""
    k, h, Tc = _kin(), 19.56, 476.0
    g = lambda T: k.heat(T, inlet(0.018)) - h * (T - Tc)
    root = find_roots(g, Tc, 527.0, n_scan=64, xtol=1e-12)["roots"]
    assert len(root) == 1
    tb = CooledTube(k, 0.018, Tc, h, CP, 60.0 * CP / h, 100, mode="frozen").solve()
    assert abs(tb.outlet()[6] / root[0] - 1.0) < 1e-12


def test_stored_claim():
    """The committed critical cooling at M is the check row's within its rtol."""
    v, rtol = _check("hcrit_MV_M_F6")
    assert abs(load_claims(PAGE)["hcrit_MV_M_F6"] / v - 1.0) <= rtol
