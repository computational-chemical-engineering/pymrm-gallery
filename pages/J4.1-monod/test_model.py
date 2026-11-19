"""Fast checks of model.py for page J4.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def monod(S, X):
    """Rawlings' Fig. 10.12 set: mu_m = 1, K_s = 1."""
    return S/(1.0 + S)


def test_chemostat_reproduces_rawlings_10_19():
    """The marched chemostat lands on eq. (10.19): S = D K_s/(mu_m - D), X = y (S_f - S)."""
    S, X = model.Chemostat(monod, 1.0, 5.0, 0.5).steady()
    assert abs(S - 1.0) < 1e-9 and abs(X - 4.0) < 1e-9


def test_plug_flow_outlet_is_the_batch_state():
    """Upwind plug flow converges at first order to the batch state at t = tau."""
    ref = model.batch_reference(monod, 0.3, 50.0, 1.0, 2.0)
    e = [abs(model.PlugFlowFermenter(monod, 0.3, 50.0, 1.0, 2.0, ncell=n).solve().outlet()[0]
             - ref[0])/ref[0] for n in (100, 200)]
    assert e[1] < 1e-2 and abs(np.log2(e[0]/e[1]) - 1.0) < 0.1


def test_washout_matches_stored_value():
    """The washout root of the marched chemostat against Rawlings' D_c = 5/6, as stored."""
    Dc = model.one_root(lambda D: model.Chemostat(monod, 1.0, 5.0, D).steady()[1] - 1e-6*5.0,
                        0.05, 2.0, n_scan=16, xtol=1e-9)
    stored = load_claims(PAGE)["rawlings_Dc_relerr"]
    assert abs(abs(Dc - 5/6)/(5/6) - stored) <= 0.05*stored
