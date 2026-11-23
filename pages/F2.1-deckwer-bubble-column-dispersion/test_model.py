"""Fast checks of model.py for page F2.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (CoupledColumn, adm_printed_A6_A11, adm_profile,  # noqa: E402
                   coupled_bvp, cstar_coefficients, danckwerts_outlet_ratio)

PAGE = HERE.name
PE, ST, LAM, H = 12.4805248618, 1.60170774114, 0.0844048147856, 0.683204540882
TH0 = 0.6 / 9.86


def test_closed_vessel_limit():
    """With c* = 0 the closed form equals Danckwerts' closed-vessel outlet."""
    assert abs(adm_profile(1.0, PE, ST, 0, 0, 0, 1.0) / danckwerts_outlet_ratio(PE, ST) - 1) < 1e-12


def test_printed_A6_A11_solve_the_danckwerts_inlet():
    """Eqs. A6 to A11 as printed agree with the corrected inlet (misprint M1 is in A4 only)."""
    a, b, _ = cstar_coefficients(0.21, 0.0, 0.21, H)
    assert abs(adm_printed_A6_A11(1.0, PE, ST, a, b, TH0)
               / adm_profile(1.0, PE, ST, a, b, 0.0, TH0) - 1) < 1e-12


def test_coupled_fv_against_bvp():
    """The pymrm solve at n = 640 is within its first-order error of solve_bvp."""
    m = CoupledColumn(PE, ST, LAM, TH0, 1.0, 0.21, 0.21, H, n=640).solve()
    sol = coupled_bvp(PE, ST, LAM, TH0, 1.0, 0.21, 0.21, H, tol=1e-8)
    th0f, th1f, psi1 = m.boundary()
    assert abs(th1f / sol.sol(1.0)[0] - 1) < 1e-3
    assert abs(psi1 / sol.sol(1.0)[2] - 1) < 1e-4
    # the overall O2 balance holds on the faces
    assert abs((th1f - TH0) * LAM / (1.0 - psi1) - 1) < 1e-9


def test_stored_claim():
    """A stored claim of the page is reproduced by the closed form."""
    stored = load_claims(PAGE)
    a, b, g = cstar_coefficients(0.21, 0.0, 0.21, H)
    c1 = adm_profile(1.0, PE, ST, a, b, g, TH0)
    c1_a4 = adm_profile(1.0, PE, ST, a, b, g, TH0, inlet_sign=-1.0)
    assert abs((c1_a4 / c1 - 1) / stored["m1_A4_printed_sign_c_out_shift"] - 1) < 1e-6
