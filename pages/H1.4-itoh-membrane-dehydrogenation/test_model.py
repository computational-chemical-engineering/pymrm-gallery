"""Fast checks of model.py for page H1.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (ItohCase, MembraneReactor, ceiling, conversion_ivp,  # noqa: E402
                   kp_of_xeq, xeq_of_kp)

PAGE = HERE.name


def _case():
    vals = load_data("itoh-1987-stated-values.csv", page=PAGE)
    s = dict(zip(vals.quantity, vals.value.astype(float)))
    K_p = kp_of_xeq(s["equilibrium_conversion"], s["feed_cyclohexane_fraction"], s["pressure"])
    return ItohCase.printed(s["temperature"], s["pressure"], s["feed_cyclohexane_fraction"], K_p), s


def test_equilibrium_inverse():
    """kp_of_xeq and xeq_of_kp are inverses at the stated equilibrium conversion."""
    case, s = _case()
    assert abs(xeq_of_kp(case.K_p, case.y_C0, case.P0) - s["equilibrium_conversion"]) < 1e-12


def test_stored_claim_reproduced():
    """The committed headline conversion, held to its claim's rtol (1e-6), on 200 cells."""
    case, s = _case()
    m = MembraneReactor(case, s["u_C0"], s["v_A0"], n_z=200).solve()
    x0 = load_claims(PAGE)["conversion_model"]
    assert abs(m.conversion() - x0) <= 1e-6 * abs(x0)
    assert abs(m.conversion() - ceiling(case, s["u_C0"], s["v_A0"])) < 1e-10


def test_closure_and_ivp_midmap():
    """Mid-map: hydrogen closes on the faces in both arrangements; FV is near the IVP."""
    case, s = _case()
    co = MembraneReactor(case, 0.8e-6, 3e-5, n_z=200).solve()
    cc = MembraneReactor(case, s["u_C0"], 2e-6, n_z=200, flow="counter").solve()
    assert abs(co.h_balance_error()) < 1e-12 and abs(cc.h_balance_error()) < 1e-12
    assert abs(co.conversion() - conversion_ivp(case, 0.8e-6, 3e-5)) < 1e-3


def test_permeation_off_gives_equilibrium():
    """alpha_H = 0 and a slow feed: the tube equilibrates at the stated x_eq."""
    case, s = _case()
    m = MembraneReactor(replace(case, alpha_H=0.0), s["u_C0"] / 10, 1e-9, n_z=100).solve()
    assert abs(m.conversion() - s["equilibrium_conversion"]) < 1e-8
    assert np.all(m.w[:, 2] == 0.0)
