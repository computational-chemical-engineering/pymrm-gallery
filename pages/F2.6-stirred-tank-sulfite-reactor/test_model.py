"""Fast checks of model.py for page F2.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import FilmTank, closed_m1, hinterland, kla_eq9, richardson  # noqa: E402

PAGE = HERE.name
AL_T = hinterland(4e-4, kla_eq9(np.sqrt(500.0 * 1e4), 0.047, 2.0e-3, 0.7, 0.2), 2e-9)


def test_pymrm_against_closed_form_m1():
    """At m = 1 the bordered pymrm system meets the closed form: second order on the
    grid, and the Richardson limit of three levels to 1e-9."""
    exact = closed_m1(1.0, AL_T)
    ns = (200, 400, 800)
    sols = [FilmTank(n).solve(1.0, AL_T, 1) for n in ns]
    assert abs(sols[-1]["phi"] / exact["phi"] - 1) < 1e-6
    for k in ("phi", "ub", "g"):
        lim, _ = richardson([s[k] for s in sols], ns)
        assert abs(lim / exact[k] - 1) < 1e-9, k


def test_oxygen_balance_m2():
    """M7: absorbed = reacted in the film + passed to the bulk, to round-off, and the
    bulk balance g = Ha^2 (Al - 1) u_b^2 holds at the solution."""
    s = FilmTank(400).solve(1.0, AL_T, 2)
    assert abs(s["balance"]) < 1e-11
    assert abs(s["g"] / ((AL_T - 1.0) * s["ub"] ** 2) - 1) < 1e-9


def test_check_value_m2_phi():
    """One check value from the second route: phi at T, Ha^2 = 1, m = 2."""
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    ns = (200, 400, 800)
    lim, _ = richardson([FilmTank(n).solve(1.0, AL_T, 2)["phi"] for n in ns], ns)
    assert abs(lim / cv.loc["m2_Ha2_1_phi", "value"] - 1) < cv.loc["m2_Ha2_1_phi", "rtol"]


def test_stored_metric_reproduced():
    """The committed m = 1 phi claim, to rtol 1e-9."""
    stored = load_claims(PAGE)["m1_Ha2_1_phi"]
    assert abs(closed_m1(1.0, AL_T)["phi"] / stored - 1) < 1e-9
