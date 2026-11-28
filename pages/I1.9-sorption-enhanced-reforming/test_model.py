"""Fast checks of model.py for page I1.9 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_data  # noqa: E402
from model import Case, psig_to_atm, solve_extents, solve_species  # noqa: E402

PAGE = HERE.name
C21 = "C2.1-xu-froment-smr"


def kfun():
    """The K(T) of page C2.1, from its model.py in the checkout, K_SMR in atm^2."""
    spec = importlib.util.spec_from_file_location("c2_1_model", HERE.parent / C21 / "model.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    kin = mod.XuFroment.from_table(load_data("xu-froment-1989-parameters.csv", page=C21))

    def k(T):
        k1, k2, _ = kin.equilibrium(T)
        return k1 / (101.325 / 100.0) ** 2, k2
    return k


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_f0_state_and_routes():
    """The f = 0 state meets its check values, and route B agrees with route A."""
    case = Case(450.0, psig_to_atm(50.0), 1.0, kfun())
    st = solve_extents(case)
    for name, x in (("f0_XCH4", st.X), ("f0_ydry_H2", st.ydry["H2"]),
                    ("P_atm_50psig", case.P_atm)):
        ref, rtol = check_value(name)
        assert abs(x / ref - 1) < rtol
    assert np.max(np.abs(st.element_gaps())) < 1e-14
    yb = solve_species(case, [0.1] * 5)
    assert max(abs(yb[k] / v - 1) for k, v in st.y.items()) < 1e-10


def test_deep_removal():
    """At 1 - f = 1e-3 the dry carbon oxides meet their check value."""
    st = solve_extents(Case(450.0, psig_to_atm(50.0), 1e-3, kfun()))
    ref, rtol = check_value("f999_COx_dry_ppm")
    assert abs(st.cox_ppm / ref - 1) < rtol
