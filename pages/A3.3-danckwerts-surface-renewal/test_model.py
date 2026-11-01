"""Fast checks of model.py for page A3.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import Element, TwoSpeciesElement, age_average, eq15, rate, rate_second_order  # noqa: E402

PAGE = HERE.name


def test_coarse_element_matches_eq8():
    """A coarse deep element, extrapolated in the step, gives sqrt(Ds) to its grid error."""
    r = rate(1.0, 200, n_x=100)
    assert abs(r["raw"] - 1) < 6e-3
    assert abs(r["extrapolated"] - 1) < 1e-3


def test_mass_balance_equals_face_flux():
    """The discrete mass balance and the face flux agree to round-off, also at a Dirichlet base."""
    for el in (Element(n_x=60, r=2.0), Element(n_x=60, depth=1.0, base="mixed", stretch=3.0)):
        out = el.march(c0=0.3, n_t=100)
        assert np.max(np.abs(out["psi"][1:] / out["psi_face"][1:] - 1)) < 1e-9


def test_second_order_frozen_reagent_is_pseudo_first_order():
    """With the reagent frozen the pair is eq. 12 with r = r' c_o'; a coarse solve gives sqrt(2)."""
    r = rate_second_order(50.0, 1.0, n_x=80, n_t=60, deplete=False)
    assert abs(r["extrapolated"] / np.sqrt(2) - 1) < 2e-3
    out = TwoSpeciesElement(rp=1 / 50.0, n_x=40).march(n_t=10)
    assert max(out["nit"]) <= 4


def test_stored_claims_reproduced():
    """Eq. 8 at production settings and the eq. 16 transcription check, by their tolerance."""
    stored = load_claims(PAGE)
    r = rate(1.0, 1600, c0=0.3)["extrapolated"]
    assert abs((r / 0.7 - 1) / stored["R_vs_eq8_rel"] - 1) <= 0.05
    swapped = abs(eq15(1.0, 3.0, 1.0, 1.0, 1.0, swap=True)[0] / 4.0 - 1)
    assert abs(swapped / stored["eq17_swapped_erf_dev"] - 1) <= 0.05
    o = Element(n_x=40).march(n_t=20)
    assert np.isfinite(age_average(o["theta"], o["Q"], 1.0))
