"""Fast checks of model.py for page J1.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_iast_equal_capacity_is_extended_langmuir():
    """For equal capacities the IAST solve reproduces extended Langmuir, which it never uses."""
    P = np.array([0.01, 0.3, 2.0, 10.0])
    y = np.array([0.1, 0.4, 0.7, 0.95])
    s = model.iast(P, y, model.langmuir_pure(1.0, 10.0), model.langmuir_pure(1.0, 1.0))
    n1, n2 = model.ext_langmuir(P, y, 1.0, 10.0, 1.0, 1.0)
    assert np.max(np.abs(s["n1"] - n1) / (n1 + n2)) < 1e-12
    assert np.max(np.abs(s["n2"] - n2) / (n1 + n2)) < 1e-12


def test_dilute_limit_matches_stored_value():
    """The closed-form dilute-in-1 selectivity of pair U is 5(1 + sqrt 3) and the stored claim."""
    s1, s2 = model.dilute_selectivities(1.0, model.langmuir_pure(1.0, 10.0),
                                        model.langmuir_pure(0.5, 2.0))
    assert abs(s1 / (5 * (1 + np.sqrt(3))) - 1) < 1e-12
    assert abs(s2 / 60.0 - 1) < 1e-12
    stored = load_claims(PAGE)["selectivity_dilute1_P1"]
    assert abs(s1 - stored) <= 1e-9 * stored
