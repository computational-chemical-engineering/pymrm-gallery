"""Fast checks of model.py for page E1.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (Penetration, int_kd_s_over_v, int_s_over_v, k_diffusion, k_throughflow,  # noqa: E402
                   renewal_factor, tracer_march)

PAGE = HERE.name


def check_value(name):
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    return float(cv.loc[name, "value"]), float(cv.loc[name, "rtol"])


def test_penetration_against_higbie():
    """The FV solve, exact in time, meets 2/sqrt(pi); the amount equals the integrated face
    fluxes; the face flux from compute_boundary_values is the operator's own face gradient."""
    p = Penetration(200)
    s = p.exact_time()
    assert abs(s["amount"] - 2 / math.sqrt(math.pi)) < 1e-7
    assert abs(s["amount"] - (s["q_in"] - s["q_out"])) < 1e-12
    face = -(p.grad @ s["c"] + p.grad_bc)[0]
    assert abs(s["flux"] / face - 1) < 1e-12


def test_renewal_factors():
    """Surface renewal: the sphere gives the printed term exactly, the circle the check value."""
    assert abs(renewal_factor("sphere") - 1) < 1e-10
    ref, rtol = check_value("renewal_factor_2d")
    assert abs(renewal_factor("circle") / ref - 1) < rtol


def test_tracer_march_against_closed_form():
    """The implicit march converges to the closed form of E9, with eqn (1) along the bubble."""
    I = int_s_over_v(0.3, 0.6, 160, 280, 0.9)
    ref, rtol = check_value("T1_int_SoverV_front_s_per_cm")
    assert abs(I / ref - 1) < rtol
    kt = float(k_throughflow(11.2, 2))
    K = lambda d: kt + float(k_diffusion(0.2, 0.4, 90.0, d))
    closed = 2.44 * math.exp(-kt * I - int_kd_s_over_v(0.3, 0.6, 160, 280, 0.9, 0.2, 0.4, 90.0))
    ref, rtol = check_value("T1_VC_front_0.6_eqn1_growing")
    assert abs(closed / ref - 1) < rtol
    assert abs(tracer_march(K, 0.3, 0.6, 160, 280, 2.44, 640, 0.9) / closed - 1) < 1e-3


def test_stored_claim():
    """A stored claim is reproduced: the 100-cell amount error of the penetration solve."""
    stored = load_claims(PAGE)
    err = Penetration(100).exact_time()["amount"] - 2 / math.sqrt(math.pi)
    assert np.isclose(err, stored["s4_amount_err_n100"], rtol=1e-3, atol=0)
