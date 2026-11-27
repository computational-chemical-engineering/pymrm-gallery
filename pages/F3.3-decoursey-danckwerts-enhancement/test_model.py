"""Fast checks of model.py for page F3.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import RenewalMarch, e20_residual, e21, e23_hikita_asai, q_from_ea  # noqa: E402

PAGE = HERE.name


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_closed_forms():
    """Eq. (21) solves Eq. (20) and meets its check value; Eq. (23) and q likewise."""
    E = float(e21(24.0, 6.0))
    assert abs(e20_residual(E, 24.0, 6.0)) < 1e-14
    for name, x in (("E21_Ea6_M24", E), ("E23_HA_erf_Ea6_M24", e23_hikita_asai(24.0, 6.0)),
                    ("q_Bo_over_zAi_Ea6_r0.02", q_from_ea(6.0, 0.02))):
        ref, rtol = check_value(name)
        assert abs(x / ref - 1) < rtol
    assert abs(q_from_ea(6.0, 1.0) - 5.0) < 1e-12


def test_one_march():
    """A coarse march at the test condition: frozen B gives sqrt(1 + M); the full
    model lands within its coarse-grid error of the extrapolated check value and of
    the stored claim."""
    frozen = RenewalMarch(np.inf, 1.0, n=100, dsigma=0.08, theta_end=1e3).run()
    assert abs(frozen.danckwerts(24.0) / 5.0 - 1) < 2e-3
    h = RenewalMarch(5.0, 1.0, n=100, dsigma=0.08, theta_end=1e3).run()
    ref, _ = check_value("E_danckwerts_exact_Ea6_r1_M24")
    assert abs(h.danckwerts(24.0) / ref - 1) < 5e-4
    stored = load_claims(PAGE)
    assert abs(h.danckwerts(24.0) / stored["t_ED"] - 1) < 5e-4
    assert math.isclose(stored["t_ED"], ref, rel_tol=1e-6)
