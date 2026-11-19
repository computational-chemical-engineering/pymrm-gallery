"""Fast checks of model.py for page J4.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name
P = dict(muh=1.0, Ks=0.03, Ki=2.0, Y=0.5, theta=3.0, S0=5.0, X0=0.0, Si=10.0, Xi=0.005)


def test_closed_form_batch_time_matches_lsoda_and_stored_value():
    """The partial-fraction batch time equals an LSODA event and the stored claim."""
    t = model.t_batch(P["Xi"], 5.0, P["Si"], P)
    assert abs(t / model.batch_ivp(P["Xi"], 5.0, P["Si"], P) - 1) < 1e-9
    stored = load_claims(PAGE)["batch_time_direct_h"]
    assert abs(t - stored) <= 1e-9 * stored


def test_steady_states_and_plug_flow_solve():
    """Two nontrivial steady states with mu = 1/theta; a coarse plug-flow solve keeps S + X/Y."""
    lo, hi = model.steady_states(P)
    assert lo < np.sqrt(P["Ks"] * P["Ki"]) < hi < P["S0"]
    for s in (lo, hi):
        assert abs(float(model.mu(s, P)) * P["theta"] - 1) < 1e-12
    p = model.PlugFlowFermenter(5.0, 0.1, 3.0, P, ncell=50).solve()
    a = 5.0 + 0.1 / P["Y"]
    assert np.max(np.abs((p.c[:, 0] + p.c[:, 1] / P["Y"]) / a - 1)) < 1e-9
    assert np.allclose(p.outlet(), p.c[-1])      # outflow: the face carries the last cell
