"""Fast checks of model.py for page J3.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import SPMe, bv_invert, kappa_from_delta, pre0, v_spm  # noqa: E402

PAGE = HERE.name
PAGE_J34 = "J3.4-doyle-fuller-newman"


@pytest.fixture(scope="module")
def params():
    par = load_data("doyle-fuller-newman-1993-parameters.csv", page=PAGE_J34)
    P = {r.symbol: float(r.value) for r in par.itertuples()}
    stated = load_data("doyle-fuller-newman-1993-stated-results.csv", page=PAGE_J34)
    delta = float(stated.loc[(stated.quantity == "delta")
                             & (stated.conditions == "initial concentrations"), "value"].iloc[0])
    return P, kappa_from_delta(delta, P) / P["epsilon"] ** 1.5


def test_dfn_part_is_j34_unchanged():
    """The first part of model.py is J3.4's model.py, character for character."""
    j34 = (HERE.parent / PAGE_J34 / "model.py")
    if not j34.is_file():
        pytest.skip("J3.4 not in this checkout")
    body = j34.read_text().split('"""', 2)[2]
    mine = (HERE / "model.py").read_text().split('"""', 2)[2].split("\n\n\n# ----")[0]
    assert mine == body.rstrip("\n")


def test_spme_operator_and_salt(params):
    """Zero total flux at both ends, salt conserved, steady profile near the closed form."""
    P, kappa = params
    sp = SPMe(P, kappa, 12, 24)
    assert sp.flux_I[0] == 0.0 and sp.flux_I[-1] == 0.0
    t = np.linspace(10.0, 2.0e4, 5)
    salt = (sp.profile(t, 1.0) * (sp.eps * sp.dx)).sum(axis=1)
    assert np.max(np.abs(salt)) < 1e-12 * P["c_0"] * np.sum(sp.eps * sp.dx)
    exact = sp.steady_closed_form(1.0)
    assert np.abs(sp.steady(1.0) - exact).max() / np.ptp(exact) < 2e-4
    with pytest.raises(ValueError):
        bv_invert(0.3 * P["c_T"], 0.5, pre0(P), dict(P, alpha_a=0.3, alpha_c=0.7))


def test_stored_claim_reproduced(params):
    """The committed steady-operator error on 24 + 48 cells, held to its own rtol."""
    P, kappa = params
    sp = SPMe(P, kappa, 24, 48)
    exact = sp.steady_closed_form(1.0)
    got = np.abs(sp.steady(1.0) - exact).max() / np.ptp(exact)
    stored = load_claims(PAGE)["spme_steady_op_rel_err"]
    assert abs(got / stored - 1) <= 0.05
    assert np.all(np.isfinite(v_spm(np.linspace(0.05, 0.85, 5), 1.0, P)))
