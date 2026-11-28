"""Fast checks of model.py for page H1.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_data  # noqa: E402
from model import (ATM, PSI, Case, MembraneReactor, Thermo, conversion,  # noqa: E402
                   exchange_tube_ethane, ln_K, solve_ivp_reactor)

PAGE = HERE.name


def _case(**kw):
    th = load_data("svna-thermo.csv", page=PAGE)
    table = {r.species: (r.A, r.B, r.C, r.D, r.dHf298_J_mol, r.dGf298_J_mol) for r in th.itertuples()}
    K = float(np.exp(ln_K(823.15, Thermo(table))))
    return Case(np.array([0.834, 0.0, 0.083, 0.083]), 1.0, 30 * PSI + ATM, 15 * PSI + ATM, K, **kw)


def test_check_value_X_T():
    """The IVP route reproduces the check value X_T to its rtol."""
    cv = load_data("check_values.csv", page=PAGE).set_index("name")["value"]
    case = _case()
    x = conversion(solve_ivp_reactor(case, 0.5, 0.1, rtol=1e-12).y[:, -1], case)
    assert abs(x / cv["X_T"] - 1) < 1e-8


def test_atom_balances_fv():
    """The finite-volume solve closes C, H and Ar on the faces to round-off."""
    m = MembraneReactor(_case(), 0.5, 0.1, n_z=100).solve()
    assert m.atom_balances().max() < 1e-12


def test_exchange_closed_form():
    """q = 1, Da = 0, beta = 1: the IVP matches the closed form of E17 ii."""
    case = _case(q=np.ones(4))
    from dataclasses import replace
    case = replace(case, PP=case.PF)
    u = solve_ivp_reactor(case, 0.5, 0.0, rtol=1e-12).y[:, -1]
    assert abs(u[0] - exchange_tube_ethane(0.834, 0.5, 1.0)) < 1e-9
