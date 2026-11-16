"""Fast checks of model.py for page A4.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def _configure():
    """The page's parameters (Duncan & Toor cell) and operators, as on the page."""
    T, p, R_GAS = 308.35, 101.325e3, 8.314462
    length, vol_1, vol_2, LA_EFF = 85.9e-3, 77.99e-6, 78.63e-6, 258.1e2
    d_ms = np.zeros((3, 3))
    d_ms[0, 1] = d_ms[1, 0] = 8.33e-5
    d_ms[0, 2] = d_ms[2, 0] = 6.8e-5
    d_ms[1, 2] = d_ms[2, 1] = 1.68e-5
    model.configure(c_t=p / (R_GAS * T), d_ms=d_ms, n_c=3, n_i=2, area=length / LA_EFF,
                    vol_1=vol_1, vol_2=vol_2, t_end=16 * 3600.0,
                    **model.build_operators(40, length))


def test_ms_transient_against_stored_deviation():
    """One Maxwell-Stefan two-bulb solve reproduces the stored mean absolute
    deviation against the A4.9 measurements (mean_abs_dev_ms), rtol 1e-9."""
    _configure()
    exp = load_data("duncan-toor-1962-run1.csv", page="A4.9-duncan-toor")
    t0 = exp[exp.time_h == 0]
    x1_0, x2_0 = [np.array([float(t0[(t0.bulb == b) & (t0.species == s)].x.iloc[0])
                            for s in ["H2", "N2", "CO2"]]) for b in (1, 2)]
    x1_0, x2_0 = x1_0 / x1_0.sum(), x2_0 / x2_0.sum()
    t, h1, h2, _ = model.simulate(model.flux_ms, x1_0, x2_0)   # raises if a solve fails
    idx = {"H2": 0, "N2": 1, "CO2": 2}
    dev = np.abs([np.interp(r.time_h, t, (h1 if r.bulb == 1 else h2)[:, idx[r.species]]) - r.x
                  for r in exp.itertuples()])
    stored = load_claims(PAGE)["mean_abs_dev_ms"]
    assert abs(float(dev.mean()) - stored) <= 1e-9 * abs(stored)
