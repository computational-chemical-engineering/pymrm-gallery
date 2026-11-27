"""Fast checks of model.py for page F2.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_data  # noqa: E402
from model import (TwoClassReactor, conversion, d_required, eps_b, eps_total,  # noqa: E402
                   split)

PAGE = HERE.name
pn = load_data("krishna-1997-printed-numbers.csv", page=PAGE, dtype=str, keep_default_na=False)
PR = {r.quantity: float(r.value) for r in pn.itertuples() if r.quantity != "H0"}
P = dict(V0=PR["V_small0"], cV=PR["eq2_slope"], edf0=PR["eps_df0"], cE=PR["eq4_slope"],
         a2_low=PR["alpha2_low"], a2_high=PR["alpha2_high"], es_low=PR["eps_s_alpha2_low"],
         es_high=PR["eps_s_alpha2_high"], N=PR["N"], m=PR["eq6_exponent"], D_ref=PR["D_ref"])
CHECK = {r.name: float(r.value) for r in load_data("check_values.csv", page=PAGE).itertuples()}


def test_correlation_check_value():
    """Eq. 6 and eq. 7 at the test condition reproduce the independent check values."""
    assert abs(eps_b(0.30, 0.38, 0.36, P) / CHECK["eps_b_T"] - 1) < 1e-11
    assert abs(eps_total(np.array([0.30]), 0.38, 0.36, P)[0] / CHECK["eps_T"] - 1) < 1e-11


def test_fv_against_closed_form():
    """The pymrm solve at 800 cells is within its first-order error of the closed form."""
    s = split(0.30, 0.38, 0.0, P, H0=1.75)
    D = d_required(s, 0.9, "dispersed", 5.0)
    m = TwoClassReactor(800, 5.0, D * s["tb"], D * s["tdf"], s["beta"]).solve()
    assert abs(m.conversion() - 0.9) < 1e-3
    assert abs(conversion(D, s, "dispersed", 5.0) - 0.9) < 1e-12


def test_mass_balance_closure():
    """Large-bubble field: inlet flux minus exit face flux equals the volume integral of the sink."""
    s = split(0.30, 0.38, 0.36, P, H0=1.92)
    n, D = 400, 2.4
    m = TwoClassReactor(n, 5.0, D * s["tb"], D * s["tdf"], s["beta"]).solve()
    cb = m.u[:n]
    b_out, _ = m.exits()
    sink = D * s["tb"] * np.sum(cb) / n
    assert abs((1.0 - b_out) - sink) < 1e-12
