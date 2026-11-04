"""Fast checks of model.py for page J3.4 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import Cell, crossing, kappa_from_delta, ocp  # noqa: E402

PAGE = HERE.name


@pytest.fixture(scope="module")
def params():
    par = load_data("doyle-fuller-newman-1993-parameters.csv", page=PAGE)
    P = {r.symbol: float(r.value) for r in par.itertuples()}
    kappa = kappa_from_delta(1.95, P) / P["epsilon"] ** 1.5
    return P, kappa


def test_boundary_fluxes_and_ocp(params):
    """The operators carry exactly I and 0 on the end faces, and Eq. 16 is symmetric."""
    P, kappa = params
    cell = Cell(P, kappa, n_s=6, n_c=12)
    c, phi2, _ = cell.unpack(cell.initial(P["u_0"]))
    i2, n_salt = cell.fluxes(c, phi2, 10.0)
    assert i2[0] == pytest.approx(10.0, rel=1e-13) and abs(i2[-1]) < 1e-12
    assert abs(n_salt[0]) < 1e-15 and abs(n_salt[-1]) < 1e-15
    beta_free = dict(P, beta=0.0, zeta=0.0)
    assert float(ocp(0.5, beta_free)) == pytest.approx(P["U_theta"], abs=1e-12)


def test_short_march_identities(params):
    """One coarse march to 500 s: every solve converges, Faraday's law and salt hold."""
    P, kappa = params
    cell = Cell(P, kappa, n_s=6, n_c=12)
    r = cell.march(10.0, t_stop=500.0)
    assert r["halvings"] == 0 and not cell.failures
    u_faraday = P["u_0"] + 10.0 * r["t"] / cell.t_full
    assert np.max(np.abs(r["u"] - u_faraday)) < 1e-12
    assert np.max(np.abs(r["salt"] / r["salt"][0] - 1)) < 1e-12
    with pytest.raises(RuntimeError):
        crossing(r["t"], r["V"], 5.0)          # no crossing: it raises


def test_stored_claim_reproduced(params):
    """The committed 1.7 V utilisation at I = 10, held to its own rtol."""
    P, kappa = params
    r = Cell(P, kappa).march(10.0, stop_on_failure=True)
    got = crossing(r["u"], r["V"], P["V_cutoff"])
    stored = load_claims(PAGE)["u_at_cutoff_I10"]
    assert abs(got / stored - 1) <= 0.05
