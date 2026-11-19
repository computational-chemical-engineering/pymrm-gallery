"""Fast checks of model.py for page J4.8 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

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


def _tables():
    read = lambda name, **kw: pd.read_csv(HERE / "data" / name, comment="#", **kw)  # noqa: E731
    comp = read("henze-1987-table1-components.csv")
    stoi = read("henze-1987-table1-stoichiometry.csv", dtype=str)
    pars = read("henze-1987-table4-parameters.csv")
    return list(comp.symbol), stoi, pars


def test_charge_and_nitrogen_close_exactly():
    """Charge and nitrogen continuity are identities of the matrix with S_N2."""
    c13, stoi, _ = _tables()
    nu = model.petersen_matrix(stoi, c13)
    gamma = model.composition_vectors(c13 + ["S_N2"])
    for law in ("charge", "N"):
        assert all(r == 0 for r in model.continuity(nu, gamma[law]).values())
    stored = load_claims(PAGE)["cod_residual_p3_symbolic"]
    assert str(model.continuity(nu, gamma["COD"])[3]) == stored


def test_lane_conserves_nitrogen_and_matches_stored_leak():
    """One lane solve: nitrogen closes to round-off, the COD leak is the stored value."""
    c13, stoi, pars = _tables()
    comps = c13 + ["S_N2"]
    ix = {c: k for k, c in enumerate(comps)}
    p = {r.symbol: float(np.sqrt(r.value_min * r.value_max)) for r in pars.itertuples()}
    p["k_a"] = 0.016
    nu = model.petersen_matrix(stoi, c13)
    nu_m = model.numeric_matrix(nu, p)
    ww = dict(S_I=30., S_S=70., X_I=50., X_S=200., X_BH=0., X_BA=0., X_P=0.,
              S_O=0., S_NO=0., S_NH=30., S_ND=5., X_ND=10., S_ALK=6., S_N2=0.)
    ras = dict(S_I=30., S_S=3., X_I=1500., X_S=120., X_BH=3000., X_BA=160., X_P=1400.,
               S_O=0., S_NO=8., S_NH=1., S_ND=1., X_ND=12., S_ALK=4., S_N2=18.)
    feed = np.array([0.5 * (ww[c] + ras[c]) for c in comps])
    s = model.Lane(nu_m, p, ix, feed, n=150).solve()
    import sympy as sp
    gamma = model.composition_vectors(comps)
    sub = {model.SYMS[k]: sp.Rational(str(p[k])) for k in model.STOICH_SYMBOLS}
    w = {law: np.array([float(gamma[law][c].subs(sub)) for c in comps]) for law in gamma}
    assert abs(model.balance(s, w["N"])) < 1e-8
    leak = abs(model.balance(s, w["COD"], cod=True))
    stored = load_claims(PAGE)["reactor_cod_leak"]
    assert abs(leak - stored) <= 1e-6 * stored
