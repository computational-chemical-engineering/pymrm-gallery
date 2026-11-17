"""Fast checks of model.py for page C1.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import IsothermalBed, lin_fit, rate_eqq, rate_max_temperature, wf_quadrature  # noqa: E402

PAGE = HERE.name


def _thermo():
    return load_data("hougen-watson-1947-eqq-thermo.csv", page=PAGE).set_index("name")["value"]


def test_mechanism_d_reproduces_printed_constants():
    """The book's linearised fit of mechanism d at 200 C lands within 1 % of its printed solution."""
    A = load_data("hougen-watson-1947-tableA-rates.csv", page=PAGE)
    d = A[A.temp_C == 200]
    coef, _, _ = lin_fit("d", *(d[k].to_numpy() for k in ("pH_atm", "pU_atm", "pS_atm", "r_lbmol_per_lb_hr")))
    assert np.max(np.abs(coef - [2.764, 1.526, 1.010, 1.129]) / [2.764, 1.526, 1.010, 1.129]) < 0.01


def test_bed_converges_and_matches_quadrature():
    """Newton succeeds, the outlet face is the last cell under outflow, and a coarse bed
    is within its first-order error of the quadrature."""
    th = _thermo()
    rate = lambda pH, pU, pS: rate_eqq(473.16, pH, pU, pS, th)
    W, _ = wf_quadrature(rate, 0.998)
    bed = IsothermalBed(rate, 200, W).solve()
    assert bed.result.success
    assert np.allclose(bed.outlet(), bed.N[-1])
    assert abs(bed.conversion_out() - 0.998) < 0.01


def test_stored_claims():
    """The quadrature W/F and the rate maximum reproduce the stored claims."""
    th = _thermo()
    stored = load_claims(PAGE)
    W, _ = wf_quadrature(lambda pH, pU, pS: rate_eqq(473.16, pH, pU, pS, th), 0.998)
    assert abs(W / stored["reactor_WF_998_quadrature"] - 1) < 1e-6
    T = rate_max_temperature(th, 2.459, 0.527, 0.515, 380.0, 750.0)
    assert abs(T / stored["rate_max_T_K_run3c_composition"] - 1) < 1e-6
