"""Fast checks of model.py for page F3.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

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
T1 = pd.read_csv(HERE / "data" / "bosch1989-table1-parameters.csv", comment="#")
P = model.params(T1)


def test_physical_absorption_and_pfo_limit():
    """All reactions off gives E = 1; a pseudo-first-order film at Ha = 10 matches
    its closed form to the stored pfo_max_rel_err bound."""
    cb = model.bulk(P, 0.2)
    c_gas = 30e3 / (P["R_GAS"] * P["T"])
    f = model.Film(P, cb, c_gas, None, enable_r1=False, enable_r34=False)
    f.solve(k_ladder=(0.0,))
    assert abs(f.enhancement() - 1.0) < 1e-9
    ha = 10.0
    fp = model.Film(P, cb, c_gas, None, enable_r34=False,
                    pfo_k1=(ha * P["KL"]) ** 2 / P["DCO2"])
    fp.solve(k_ladder=(1.0,), tol=1e-8, bounds=None)
    exact = model.pfo_enhancement(ha, cb["CO2"] / model.c_interface(P, 30.0))
    assert abs(fp.enhancement() / exact - 1.0) <= 1e-4


def test_unpromoted_fit_matches_stored_value():
    """The scalar that reconciles the unpromoted bulk CO2 with Table 2 at one
    condition lies within the stored per-condition spread of k3_scale_fitted."""
    stored = load_claims(PAGE)
    s, spread = stored["k3_scale_fitted"], stored["k3_scale_spread_per_condition"]
    t2 = pd.read_csv(HERE / "data" / "bosch1989-table2-results.csv", comment="#")
    r = t2.iloc[0]
    implied = model.c_interface(P, r.p_co2_kpa) - r.J_u / (r.E_u * P["KL"])
    per = implied / model.bulk(P, r.alpha)["CO2"]
    assert abs(per / s - 1.0) <= spread
    assert np.isfinite(model.vkh(10.0, 50.0))
