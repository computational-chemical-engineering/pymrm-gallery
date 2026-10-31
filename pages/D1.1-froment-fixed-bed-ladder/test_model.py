"""Fast checks of model.py for page D1.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import json
import sys
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (Ladder, Pellet, Reactor, eta_sphere, ode_reference,  # noqa: E402
                   p0_critical, p0_critical_ode, phi_prime)

PAGE = HERE.name
D22 = "D2.2-van-welsenaere-froment-runaway"


def _reactor():
    """The design case, built from the CSVs as the notebook builds it."""
    bk = load_data("froment-2011-ch11-printed.csv", page=PAGE).set_index("key")["value"]
    bk = {k: float(v) for k, v in bk.items() if _is_number(v)}
    par = load_data("van-welsenaere-froment-1970-parameters.csv", page=D22)
    v = dict(zip(par.symbol, par.value))
    d_t, d_p = bk["tube_id"], bk["particle_diameter"]
    eps = bk["voidage_a"] + bk["voidage_b"] * (1.0 + (d_t / d_p - 2.0) ** 2 / (d_t / d_p) ** 2)
    return Reactor(M=v["M"], p_t=v["P"], rho_b=v["rho_b"], rho_g=v["rho_g"],
                   minus_dH=v["minus_dH"], cp_vol=v["c_p"], u=v["u"], U=v["U"], R=v["R"],
                   p_B0=v["p_B0"], a=v["a"], b=v["b"], length=bk["tube_length"], d_p=d_p,
                   T_r=bk["coolant_temperature"], T_perm=bk["temp_range_high"], eps=eps,
                   lam_er=bk["lam_er"], alpha_w=bk["alpha_w"], pe_mr=bk["pe_mr"],
                   pe_hr=bk["pe_hr"], re_film=bk["Re_1173"]), bk


def _is_number(x):
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def test_rung1_threshold_against_grid_free_route():
    """Rung 1's bisected threshold on 300 cells (which p0_critical's crossing
    control has already accepted) lies within first-order upwind's error of the
    grid-free Radau threshold."""
    rx, _ = _reactor()
    p_grid = p0_critical(Ladder(rx, rung=1, n_z=300))
    p_ode = p0_critical_ode(rx)
    assert abs(p_grid / p_ode - 1.0) < 1e-2


def test_pellet_against_closed_form():
    """The spherical pellet's eta against eq. (11.9.1-11), second order in n_r."""
    rx, _ = _reactor()
    T_s, p_s = 680.0, 0.01
    phi = phi_prime(rx, T_s, 1.0e-6)
    e = Pellet(rx, n_r=64).solve(rx.theta * p_s, T_s)["eta"]
    assert abs(e / eta_sphere(phi) - 1.0) < 1e-3


def test_stored_claim_reproduced():
    """T_hot_0018_bar_reading_C: rung 1 in the printed bar, held to its claim's rtol."""
    rx, bk = _reactor()
    rx_bar = replace(rx, p_B0=bk["oxygen_pressure"], p_t=bk["total_pressure"])
    T = ode_reference(rx_bar, bk["runaway_p0_safe"])[0] - 273.15
    claim = json.loads((HERE / "agreement.json").read_text())["metrics"][
        "T_hot_0018_bar_reading_C"]
    x0 = load_claims(PAGE)["T_hot_0018_bar_reading_C"]
    assert abs(T - x0) <= claim["rtol"] * abs(x0)
