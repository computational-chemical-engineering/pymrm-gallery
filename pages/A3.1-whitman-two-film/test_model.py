"""Fast checks of model.py for page A3.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import Equilibrium, TwoFilm, reacting_film_rate, two_film_algebraic  # noqa: E402

PAGE = HERE.name


def _equilibrium():
    t1 = load_data("whitman-1923-table1.csv", page=PAGE)
    pairs = np.unique(np.r_[t1[["c1_g_per_l", "p1_mmHg"]].to_numpy(float),
                            t1[["c3_g_per_l", "p3_mmHg"]].to_numpy(float)], axis=0)
    pairs = pairs[pairs[:, 1] > 0]
    return t1, Equilibrium.fit(*pairs.T)


def test_pymrm_solve_matches_algebra():
    """The two-film solve reproduces Whitman's eq. 4 (a linear profile is exact)."""
    t1, eq = _equilibrium()
    r = t1.set_index("run")
    tf = TwoFilm(float(r.loc[3, "K_p"]), 1.2, eq.p, eq.dpdc, n=8)
    s = tf.solve(float(r.loc[1, "p1_mmHg"]), float(r.loc[1, "c3_g_per_l"]))
    a = tf.algebraic(float(r.loc[1, "p1_mmHg"]), float(r.loc[1, "c3_g_per_l"]))
    assert s["success"] and s["nit"] <= 10
    assert abs(s["N"] / a["N"] - 1) < 1e-11
    assert abs(s["N_liquid_face"] / s["N_gas_face"] - 1) < 1e-10


def test_reacting_film_second_order():
    """With a reaction in the liquid film the solve converges to the exact film at order 2."""
    t1, eq = _equilibrium()
    r = t1.set_index("run")
    p1, c3, kp = float(r.loc[1, "p1_mmHg"]), float(r.loc[1, "c3_g_per_l"]), float(r.loc[3, "K_p"])
    ex = reacting_film_rate(p1, c3, kp, 1.2, 10.0, eq.p, eq.c)["N"]
    e = [abs(TwoFilm(kp, 1.2, eq.p, eq.dpdc, n=n, ka=10.0).solve(p1, c3)["N"] / ex - 1)
         for n in (20, 40)]
    assert 1.9 < np.log2(e[0] / e[1]) < 2.1


def test_stored_claims_reproduced():
    """The equilibrium residual and the held-out run-2 rate, by each claim's tolerance."""
    stored = load_claims(PAGE)
    t1, eq = _equilibrium()
    r = t1.set_index("run")
    c = np.r_[t1.c1_g_per_l, t1.c3_g_per_l].astype(float)
    p = np.r_[t1.p1_mmHg, t1.p3_mmHg].astype(float)
    keep = p > 0
    assert abs(eq.residuals(c[keep], p[keep])[0] / stored["eq_curve_worst_resid_pct"] - 1) <= 0.01
    kp = float(r.loc[3, "K_p"])
    p2 = float(r.loc[1, "p1_mmHg"]) - float(r.loc[1, "dW_dtheta_g_per_h"]) / kp
    kc = float(r.loc[1, "dW_dtheta_g_per_h"]) / (float(eq.c(p2)) - float(r.loc[1, "c3_g_per_l"]))
    n2 = two_film_algebraic(float(r.loc[2, "p1_mmHg"]), float(r.loc[2, "c3_g_per_l"]),
                            kp, kc, eq.p)["N"]
    dev = abs(n2 / float(r.loc[2, "dW_dtheta_g_per_h"]) - 1)
    assert abs(dev / stored["run2_rate_reldev"] - 1) <= 0.01
