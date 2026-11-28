"""Fast checks of model.py for page H1.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

import numpy as np  # noqa: E402
from gallery_utils import load_claims, load_data  # noqa: E402
from model import MembraneModule, elasticities, local_flux, module_number  # noqa: E402

PAGE = HERE.name
THETA, PHI = 1.0 / 14.95, 0.425 / 0.575


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_local_forms_and_check_values():
    """The three flux forms agree; flux and elasticity meet their check values."""
    pi = np.linspace(0.1, 1.0, 7)
    forms = [local_flux(pi, THETA, 0.7, f) for f in ("ratio", "membrane", "film")]
    assert np.max(np.abs(forms[1] - forms[0])) < 1e-14 and np.max(np.abs(forms[2] - forms[0])) < 1e-14
    for name, x in (("TC_j_inlet", local_flux(1.0, THETA, 1.0)),
                    ("TC_epsQ_inlet", elasticities(1.0, THETA, 1.0)[0])):
        ref, rtol = check_value(name)
        assert abs(x / ref - 1) < rtol


def test_one_module_solve():
    """A coarse module solve: H2 conserved to round-off, the area within its first-order
    grid error of the quadrature and of the stored claim."""
    m = MembraneModule(THETA, PHI, 1.0, 100).solve_area(0.8)
    assert m.balance_error() < 1e-12
    ref, _ = check_value("TC_N_R80_beta1")
    assert abs(m.N / ref - 1) < 1e-2
    assert abs(module_number(0.8, THETA, PHI, 1.0) / ref - 1) < 1e-10
    stored = load_claims(PAGE)
    assert abs(m.N / stored["r4_N80_pymrm"] - 1) < 1e-2
