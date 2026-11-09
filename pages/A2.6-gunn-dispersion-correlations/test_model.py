"""Fast checks of model.py for page A2.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import Bed, D_over_Ud, D_over_Ud_1993, Pe_axial, Pe_radial, dD_over_Ud_dRe  # noqa: E402

PAGE = HERE.name
G = dict(zip(*load_data("gunn-1987-stated.csv", page=PAGE)[["quantity", "value"]].T.values))
EPS, SC = G["voidage_bulk_spheres"], G["schmidt_argon_air"]


def test_sphere_limit_is_two():
    """Eq. (1) at Re = 1e8 is within 1e-7 of 2 p_inf / (1 - p_inf) = 2 for spheres."""
    assert abs(float(Pe_axial(1e8, SC, G, EPS)) / 2.0 - 1) < 1e-7


def test_eq1_matches_eq42_and_its_derivative():
    """The 1987 and 1993 forms agree; the analytic derivative matches a central difference."""
    Re = np.logspace(-2, 3, 41)
    assert np.max(np.abs(D_over_Ud(Re, SC, G, EPS) / D_over_Ud_1993(Re, SC, G, EPS) - 1)) < 1e-9
    h = 1e-6 * Re
    fd = (D_over_Ud(Re + h, SC, G, EPS) - D_over_Ud(Re - h, SC, G, EPS)) / (2 * h)
    scale = D_over_Ud(Re, SC, G, EPS) / Re
    assert np.max(np.abs(dD_over_Ud_dRe(Re, SC, G, EPS) - fd) / scale) < 1e-6


def test_newton_and_deferred_agree():
    """The two solver routes reach the same limited equations on a small bed."""
    pe_ax = float(Pe_axial(50.0, SC, G, EPS))
    pe_r = float(Pe_radial(50.0, SC, G, EPS))
    b = Bed(1 / pe_ax, 1 / pe_r, 0.02, L=100.0, n_z=200, n_r=8, inlet_mode=1).solve()
    c_newton = b.c.copy()
    b.solve(route="deferred")
    assert np.max(np.abs(c_newton - b.c)) < 1e-12


def test_stored_claim_reproduced():
    """The committed eq17_flat_sup_error, held to the rtol of its claim."""
    stored = load_claims(PAGE)["eq17_flat_sup_error"]
    pe_ax = float(Pe_axial(50.0, SC, G, EPS))
    pe_r = float(Pe_radial(50.0, SC, G, EPS))
    b = Bed(1 / pe_ax, 1 / pe_r, 0.02, L=100.0, d_t=20.0, n_z=800, n_r=20).solve()
    ex = b.closed_form()
    err = float(np.max(np.abs(b.c - ex)) / np.max(np.abs(ex)))
    assert abs(err - stored) <= 1e-3 * abs(stored)
