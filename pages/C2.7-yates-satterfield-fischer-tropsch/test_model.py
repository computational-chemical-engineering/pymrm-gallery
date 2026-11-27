"""Fast checks of model.py for page C2.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (PlugFlow, cstr, cstr_wf, e_app, e_app_root, feed_flow,  # noqa: E402
                   pfr_wf)

PAR = dict(a240=75.76, b240=11.61, a220=53.11, b220=22.26, T1=493.15, T2=513.15)


def test_e_app_root_closed_form():
    """E_app at its closed-form root equals the quoted value."""
    p = e_app_root(94.5, PAR)
    assert abs(float(e_app(p, PAR)) - 94.5) < 1e-10


def test_cstr_against_closed_form():
    """The newton CSTR meets its closed-form W/F at the Table II feed."""
    tau = 3.0 / feed_flow(0.067)
    c = cstr(PAR["a240"], PAR["b240"], 0.79, tau)
    assert abs(cstr_wf(c["X"], PAR["a240"], PAR["b240"], 0.79) / tau - 1.0) < 1e-11


def test_pfr_first_order():
    """The pymrm PFR approaches the closed form at first order."""
    tau = 3.0 / feed_flow(0.067)
    errs = []
    for n in (100, 200):
        x = PlugFlow(PAR["a240"], PAR["b240"], 0.79, n=n).solve(tau)
        errs.append(abs(pfr_wf(x, PAR["a240"], PAR["b240"], 0.79) / tau - 1.0))
    assert 1.8 < errs[0] / errs[1] < 2.2


def test_stored_claim_reproduced():
    """The committed lower end of the E_app band, held to its claim's rtol."""
    x0 = load_claims(HERE.name)["r3_band_lo"]
    assert np.isclose(float(e_app(0.145, PAR)), x0, rtol=0.0112)
