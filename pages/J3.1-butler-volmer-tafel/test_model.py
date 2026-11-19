"""Fast checks of model.py for page J3.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name
PAR = {r.symbol: float(r.value) for r in
       load_data("doyle-fuller-newman-1993-parameters.csv", page="J3.4-doyle-fuller-newman").itertuples()}


def test_eq17_reduces_to_eq30_when_alphas_sum_to_one():
    """Eq. 17 with U_2 from Eq. 16 is Eq. 30 in the standard form, for alpha_a + alpha_c = 1."""
    kin = model.Kinetics(PAR)
    c, cs = np.array([300.0, 1000.0, 3000.0]), np.array([500.0, 14500.0, 28000.0])
    eta_s = np.array([-0.1, 0.02, 0.2])
    for aa in (0.5, 0.3):
        lhs = kin.i_eq17(kin.U_ocp(cs / kin.c_T) + eta_s, c, cs, aa, 1 - aa)
        rhs = kin.i_bv(eta_s, kin.i0_eq30(c, cs, aa, 1 - aa), aa, 1 - aa)
        assert np.max(np.abs(lhs / rhs - 1)) < 1e-10


def test_tafel_threshold_matches_stored_claim():
    """The closed-form Tafel threshold equals its root on the full law and the stored claim."""
    kin = model.Kinetics(PAR)
    closed = kin.eta_tafel_threshold(0.10)
    assert abs(kin.eta_tafel_threshold_numeric(0.10, 0.3, 0.7) / closed - 1) < 1e-12
    stored = load_claims(PAGE)["eta_tafel_10pct_mV"]
    assert abs(1e3 * closed - stored) <= 1e-9 * stored


def test_first_crossing_is_not_the_last():
    """At alpha_a = 0.3 the 5 % linear threshold is the first of several crossings."""
    kin = model.Kinetics(PAR)
    first = kin.eta_linear_threshold(0.05, 0.3, 0.7)
    last = kin.eta_linear_threshold(0.05, 0.3, 0.7, crossing="last")
    assert last > 10 * first


def test_cell_first_steps_conserve_salt():
    """Ten implicit-Euler steps of the coarse cell converge and conserve salt."""
    kin = model.Kinetics(PAR)
    stated = {r.quantity: float(r.value) for r in load_data(
        "doyle-fuller-newman-1993-stated-results.csv", page="J3.4-doyle-fuller-newman").itertuples()}
    kappa = model.kappa_from_delta(PAR, stated["delta"]) / PAR["epsilon"] ** 1.5
    out, _ = model.Cell(kin, kappa, n_s=6, n_c=12).march(10.0, max_steps=10)
    assert len(out) == 10
    assert abs(out[:, 3].max() / out[:, 3].min() - 1) < 1e-12
