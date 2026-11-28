"""Fast checks of model.py for page G1.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_data  # noqa: E402
from model import (MOMENT_EXACT, eta9, eta9_reduced, metrics, pellet_moment_coef)  # noqa: E402

PAGE = HERE.name
K = {"c9": 1.104, "n9": 1 / 3, "m9": 1 / 9, "c1": 1.617, "n1": 0.1461, "m1": -0.0711}
CV = load_data("check_values.csv", page=PAGE).set_index("name")["value"]


def test_eq9_full_groups_equal_reduced_form():
    """Eq. 9 with its full Re and Ga groups against K L^(1/3) eps^(-1/3) psi^(1/9), and the check value."""
    L = np.geomspace(0.4, 2.8, 25)
    full = eta9(L, 0.5, 1.52e-3, 0.412, 654.8, 3.0e-4, K, cap=False)
    red = eta9_reduced(L, 0.5, 0.412, 654.8, 3.0e-4, K)
    assert np.max(np.abs(full / red - 1)) < 1e-13
    assert abs(eta9(1.0, 0.5, 1.52e-3, 0.412, 654.8, 3.0e-4, K) / CV["eta9_sphere"] - 1) < 1e-11


def test_sphere_moment_coefficient_converges_to_one_fifteenth():
    """The pymrm pellet BVP: second-order approach to 1/15, Richardson value within 1e-8."""
    q = [pellet_moment_coef(n, 2) for n in (80, 160, 320)]
    order = np.log2((q[0] - q[1]) / (q[1] - q[2]))
    assert 1.95 < order < 2.05
    rich = q[2] + (q[2] - q[1]) / 3.0
    assert abs(rich - MOMENT_EXACT[2]) < 1e-8
    assert abs(MOMENT_EXACT[2] / CV["tau1_coef_sphere"] - 1) < 1e-11


def test_metric_row():
    """The synthetic metric set against its check values, capped, (p - e)/e."""
    synth = [(1.0, 0.5, 0.70, 1.52e-3, 0.412), (2.0, 4.0, 0.98, 1.52e-3, 0.412),
             (0.7, 3.0, 0.66, 1.9914454e-3, 0.355), (1.5, 0.1, 0.85, 1.9914454e-3, 0.355)]
    p = np.array([eta9(a, b, d, e, 654.8, 3.0e-4, K) for a, b, _, d, e in synth])
    m = metrics(p, [s[2] for s in synth], "exp")
    assert m["inside"] == CV["synth_inside10_capped_exp"]
    assert abs(m["maxdev"] / CV["synth_maxdev_capped_exp"] - 1) < 1e-6
