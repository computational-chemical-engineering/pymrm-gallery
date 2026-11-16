"""Fast checks of model.py for page A4.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name
PAGE_A42 = "A4.2-maxwell-stefan-vs-fick"


def _gas():
    """The page's capillary constants: H2/N2/CO2 at 35.2 C and 1 atm, 8.6 cm."""
    K = load_data("krishna-wesselingh-1997-dgm-constants.csv", page=PAGE)
    P = {(r.group, r.symbol): r.value for r in K.itertuples()}
    wx = load_data("krishna-wesselingh-1997-worked-examples.csv", page=PAGE_A42)
    W = {(r.example, r.symbol): r.value for r in wx.itertuples()}
    T, p = 308.35, 101325.0
    ct = p / (model.RGAS * T)
    pairs = {(0, 1): W[("ideal", "D12")], (0, 2): W[("ideal", "D13")],
             (1, 2): W[("ideal", "D23")]}
    molar = np.array([P[("molar", f"M_{s}")] for s in ["H2", "N2", "CO2"]])
    return dict(molar=molar, dij_ct=model.pair_matrix(pairs, ct), T=T, p=p, L=0.086)


def test_capillary_against_collocation_and_stored_value():
    """At 1 mm the finite-volume flux (100 cells) matches the collocation route to
    the stored collocation_worst_dev (rtol 1e-3), which sits at this pore size."""
    gas = _gas()
    xa, xb = [0.0, 0.49, 0.51], [0.5, 0.5, 0.0]
    m, c = model.capillary(1e-3, xa, xb, **gas)            # raises if Newton fails
    nb = model.capillary_bvp(1e-3, xa, xb, **gas)
    dev = float(np.max(np.abs(m.flux(c)[50] / nb - 1)))
    stored = load_claims(PAGE)["collocation_worst_dev"]
    assert abs(dev - stored) <= 1e-3 * stored


def test_pellet_knudsen_limit():
    """A -> 2B under Knudsen control approaches sqrt(2) at the centre."""
    m, c = model.pellet(2.0, 1e4, 60.0)
    assert abs(m.faces(c)[0][0].sum() / np.sqrt(2.0) - 1.0) < 1e-4
