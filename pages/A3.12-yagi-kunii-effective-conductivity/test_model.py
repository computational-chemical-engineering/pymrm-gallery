"""Fast checks of model.py for page A3.12 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import R_I, R_O, eq15, eq16, ke0_iron, h_rs, h_rv, set_air_table, set_iron_sphere_row, solve_annulus  # noqa: E402,E501


def test_eq16_is_eq15_without_radiation():
    """Eq. (16) is Eq. (15) with both radiation coefficients zero and gamma = 1."""
    a = eq15(0.4, 1.0, 1.0, 0.034, 11.0e-3, 45.0, 0.0447, 0.0, 0.0)
    b = eq16(0.4, 1.0, 0.034, 45.0, 0.0447)
    assert abs(a / b - 1) < 1e-12
    assert h_rs(400.0) > 0 and h_rv(400.0, 0.4) > 0


def test_annulus_profile_stored_metric():
    """The constant-k annulus on 800 cells against the exact logarithmic profile,
    held to the page's stored V4_annulus_profile_err_n800 at rtol 1e-9."""
    stored = load_claims(HERE.name)["V4_annulus_profile_err_n800"]
    Tin, Tout = 500.0, 100.0
    _, r_c, T, _ = solve_annulus(800, Tin, Tout, k_const=1.0)
    exact = Tout + (Tin - Tout) * np.log(R_O / r_c) / np.log(R_O / R_I)
    err = float(np.max(np.abs(T - exact)) / (Tin - Tout))
    assert abs(err - stored) <= 1e-9 * abs(stored)


def test_newton_and_fixed_point_agree():
    """The k(T) annulus: Newton and the fixed point give the same profile, and the
    heat rate is the same through every face."""
    import pandas as pd
    data = HERE / "data"
    set_air_table(pd.read_csv(data / "yagi-kunii-1957-air-conductivity.csv", comment="#"))
    t4 = pd.read_csv(data / "yagi-kunii-1957-table4.csv", comment="#")
    set_iron_sphere_row(t4[t4.solid == "Iron sphere"].iloc[0])
    _, _, Tn, Qn = solve_annulus(200, 858.0, 598.0, kfun=ke0_iron)
    _, _, Tp, _ = solve_annulus(200, 858.0, 598.0, kfun=ke0_iron, method="picard")
    assert np.max(np.abs(Tn - Tp)) < 1e-8
    assert np.ptp(Qn) / abs(Qn[0]) < 1e-10
