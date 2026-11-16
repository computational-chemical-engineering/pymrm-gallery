"""Fast checks of model.py for page A4.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_jackson_limit_and_stored_value():
    """A -> 2B under Knudsen control: the centre pressure matches the stored
    jackson_nu2_solved_increase (rtol 1e-6) on the page's 600-cell sphere."""
    p, c = model.pellet(1e6, 2.0, 30.0, n=600)              # raises if Newton fails
    inc = p.centre_total(c) - 1.0
    stored = load_claims(PAGE)["jackson_nu2_solved_increase"]
    assert abs(inc - stored) <= 1e-6 * stored
    assert abs(inc + 1.0 - np.sqrt(2.0)) < 1e-5


def test_isobaric_slab_against_closed_form():
    """H2/N2 slab at Kn = 1: the finite-volume flux matches eq. (109) integrated."""
    mr = 28.0134 / 2.01588
    s, c = model.isobaric_slab(1.0, mr, n=100)
    dam, dab = model.d_groups(1.0)
    ref = model.isobaric_flux_closed_form(dam, dab, mr ** -0.5, 0.9, 0.1)
    assert abs(s.flux(c)[50, 0] / ref - 1.0) < 1e-5
