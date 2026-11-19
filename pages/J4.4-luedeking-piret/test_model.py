"""Fast checks of model.py for page J4.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def _cf():
    """Closed form on a smooth synthetic k(t), independent of the page's data."""
    t = np.linspace(0.0, 10.0, 11)
    return model.ClosedForm(t, 0.5 - 0.04 * t, 0.1)


def test_plug_flow_matches_closed_form_at_first_order():
    """Outflow outlet: boundary read equals the last cell, and the error halves with h."""
    cf = _cf()
    p_ref = cf.P(cf.t_hi, 2.0, 0.5)
    err = []
    for n in (100, 200):
        f = model.PlugFlowFermenter(n, 2.0, 0.5, cf).solve()
        assert f.outlet()[1] == f.outlet_last_cell()[1]
        err.append(abs(f.outlet()[1] / p_ref - 1))
    assert 1.8 < err[0] / err[1] < 2.2


def test_batch_march_first_order():
    cf = _cf()
    p_ref = cf.P(cf.t_hi, 2.0, 0.5)
    e = [abs(model.batch_march(n, 2.0, 0.5, cf)[1] / p_ref - 1) for n in (50, 100)]
    assert 1.8 < e[0] / e[1] < 2.2


def test_refit_matches_stored_value():
    """The eq. (3) refit of the printed tables reproduces the stored claim."""
    d = HERE / "data"
    t1 = pd.read_csv(d / "luedeking1959-table1-bacterial-density.csv", comment="#")
    t2 = pd.read_csv(d / "luedeking1959-table2-growth-acid.csv", comment="#")
    t3 = pd.read_csv(d / "luedeking1959-table3-alpha-beta-vs-ph.csv", comment="#")
    r = model.refit_all(t1, t2, t3)
    stored = load_claims(PAGE)["alpha_refit_eq3"]
    assert abs(r["alpha_refit_eq3"] - stored) <= 1e-9 * stored
