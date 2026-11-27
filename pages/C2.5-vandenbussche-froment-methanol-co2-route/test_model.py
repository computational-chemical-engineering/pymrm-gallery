"""Fast checks of model.py for page C2.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import ATOMS, Kinetics, PlugFlowBed, Thermo, reference_solution, table2  # noqa: E402

PAGE = HERE.name
Y_FEED = np.array([0.04, 0.03, 0.82, 0.0, 0.0, 0.11])     # Table 3


def _kin():
    return Kinetics(table2(load_data("vandenbussche-1996-table2.csv", page=PAGE, dtype=str)))


def _thermo(**kw):
    return Thermo(load_data("nist-webbook-thermo.csv", page=PAGE, dtype=str),
                  load_data("nist-webbook-methanol-cp.csv", page=PAGE, dtype=str), **kw)


def test_check_value_rate():
    """One row of check_values.csv: the methanol rate at the test condition."""
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    p = 50.0 * np.array([0.030, 0.025, 0.800, 0.010, 0.005, 0.130])
    r = _kin().rates(p, 523.15)
    assert abs(r[0] / cv.value["TC_r_MeOH"] - 1.0) < cv.rtol["TC_r_MeOH"]


def test_bed_balances_and_outlet_against_radau():
    """Adiabatic bench bed: C, H, O conserved at every cell; outlet T meets solve_ivp at first order."""
    th = _thermo()
    w_f0 = 0.0348 / (2.8e-5 / th.m_feed_of(Y_FEED))
    beds = [PlugFlowBed(_kin(), Y_FEED, 493.2, 50.0, w_f0, thermo=th, n=n).solve()
            for n in (400, 800)]
    assert np.abs((beds[1].u[:, :6] - Y_FEED) @ ATOMS.T).max() < 1e-12
    T_ref = 493.2 * reference_solution(beds[0]).y[6, -1]
    gaps = [abs(493.2 * b.outlet()[6] - T_ref) for b in beds]
    assert gaps[1] < 0.05 and beds[1].backward_error() < 1e-8


def test_stored_claim_reproduced():
    """The committed plateau temperature of Fig. 5 (298 K heats, argon), to its claim's rtol."""
    th = _thermo()
    w_f0 = 0.0348 / (2.8e-5 / th.m_feed_of(Y_FEED))
    bed = PlugFlowBed(_kin(), Y_FEED, 493.2, 50.0, w_f0, thermo=th, n=600).solve()
    x0 = load_claims(PAGE)["r3_T_plateau"]
    assert abs(493.2 * bed.outlet()[6] - x0) <= 1e-3 * x0
