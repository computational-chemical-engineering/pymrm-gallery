"""Fast checks of model.py: one check value, the single-phase Ergun limit, two routes at T."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

import model as M  # noqa: E402
from gallery_utils import load_data  # noqa: E402

CV = load_data("check_values.csv", page=HERE.name).set_index("name")

T = dict(d_e=3.0e-3, eps=0.38, L=5.0, G=0.1, rho_L=998.0, mu_L=1.0e-3, sigma=0.0728,
         rho_G=1.2, mu_G=1.8e-5, g=9.81)


def test_check_value_psi():
    psi = float(np.squeeze(M.state(T, M.params())["psi"]))
    row = CV.loc["T_psi_G"]
    assert abs(psi / row.value - 1) < row.rtol


def test_ergun_limit():
    P = M.params()
    g = M.groups(T, P)
    dim = M.ergun_gradient(T["rho_G"], T["mu_G"], T["G"], T, P["A"], P["B"])
    assert abs(T["rho_G"] * T["g"] * g["X_G"] / dim - 1) < 1e-13


def test_newton_against_brentq():
    P = M.params()
    d_n, psi_n, res = M.solve_pair(T, P)
    d_b = M.solve_scalar(T, P)
    assert res.success and abs(float(d_n[0]) / float(d_b[0]) - 1) < 1e-10
