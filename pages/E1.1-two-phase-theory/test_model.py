"""Fast checks of model.py for page E1.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import PATM, CompressibleBed, mu_air, rho_air  # noqa: E402

PAGE = HERE.name

# The reference run of the page: bead 3 at G = 66 lb/(hr ft2), read from the page's CSVs.
INH2O, IN, FT = 1000.0 * 9.80665 * 0.0254, 0.0254, 0.3048
G_SI = 0.45359237 / (3600 * 0.09290304)
T = 366.483                                            # 200 F
_t1 = load_data("toomey_johnstone_1952_table1_beads.csv", page=PAGE)
_t2 = load_data("toomey_johnstone_1952_table2_pressure_losses.csv", page=PAGE)
_pr = load_data("toomey_johnstone_1952_printed_results.csv", page=PAGE)
_P = dict(zip(_pr.symbol, _pr.value))
_ref = _t2[(_t2.bead_no == 3) & (_t2.G_lb_hr_ft2 == 66.0)].iloc[0]
DP = float(_t1[_t1.bead_no == 3].dp_micron.iloc[0]) * 1e-6
W, L_MF = _ref.dPmf_inH2O * INH2O / (_ref.Lmf_in * IN), _ref.Lmf_in * IN
A = _P["eq1_coefficient"] / np.sqrt(DP / FT) / FT
U_MF = _P["Re_mf_bead3"] * mu_air(T) / (DP * rho_air(T))
R_MEAN = np.mean([_P[f"umf_over_ue_{i}"] for i in (1, 2, 3, 4)])
U_TOP = _ref.G_lb_hr_ft2 * G_SI / rho_air(T)


def bed(n=200, **kw):
    return CompressibleBed(W, A, U_MF / R_MEAN, U_TOP, L_MF, n=n, **kw)


def test_validation_routes_agree():
    """The closed form and the Runge-Kutta route give the same bed height."""
    b = bed()
    assert abs(b.Lf_ivp() / b.Lf_analytic() - 1) < 1e-10


def test_pymrm_solve_first_order_and_conservative():
    """The finite-volume height converges at first order, and the distributor
    face pressure equals P_top + w L_mf at the root."""
    err = []
    for n in (100, 200):
        b = bed(n)
        Lf = b.Lf()
        s = b.solve(Lf)
        assert abs(s["P_bottom"] - (PATM + W * L_MF)) < 1e-9 * W * L_MF
        err.append(abs(Lf / b.Lf_analytic() - 1))
    assert 1.9 < err[0] / err[1] < 2.1


def test_stored_claim_reproduced():
    """The committed L_f/L_mf at n = 1600, held to its own rtol."""
    stored = load_claims(PAGE)
    b = bed(1600)
    assert abs(b.Lf() / L_MF / stored["Lf_over_Lmf"] - 1) <= 0.05
