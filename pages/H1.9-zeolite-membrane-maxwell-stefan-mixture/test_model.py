"""Fast checks of model.py for page H1.9 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import DSL, IAST, B_matrix, Membrane, friction_solve, selectivity  # noqa: E402

PAGE = HERE.name
A47 = "A4.7-zeolite-micropore-maxwell-stefan"


def _c1_nc4():
    par = load_data("krishna-baur-2003-sorption-parameters.csv", page=A47)

    def row(component, table):
        s = par[(par.component == component) & (par.source_table == table)]
        return dict(zip(s.symbol, s.value))

    isos = [DSL(g["b_A"], g["Theta_sat_A"], g["b_B"], g["Theta_sat_B"])
            for g in (row("C1", "Table 1"), row("nC4", "Table 1"))]
    D = [row("C1", "Section 4.1")["D_MS"], row("nC4", "Section 4.1")["D_MS"]]
    return isos, D


def test_iast_roundtrip():
    """Forward then inverse IAST returns the pressures it started from."""
    isos, _ = _c1_nc4()
    mix = IAST(isos)
    p = np.array([[95e3, 5e3], [1e3, 1e3]])
    assert np.max(np.abs(mix.inverse(mix.forward(p)) / p - 1)) < 1e-10


def test_friction_two_routes():
    """eq. (44) inverted against the matrix-free eq. (42)."""
    th, D, d = np.array([0.2, 0.1]), np.array([1e-9, 1e-11]), np.array([1.0, -0.5])
    u1 = friction_solve(th, D, d)
    u2 = np.linalg.solve(B_matrix(th[None, :], D[None, :])[0], d)
    assert np.max(np.abs(u1 / u2 - 1)) < 1e-10


def test_selectivity_stored_claim():
    """The 95-5 membrane with finite exchange, direct solve, against the committed S_P."""
    isos, D = _c1_nc4()
    m = Membrane(isos, [i.Th_sat for i in isos], D, n_z=40)
    p0 = np.array([95e3, 5e3])
    _, J = m.steady(p0)
    assert not m.continued
    assert m.flux_spread < 1e-9
    x0 = load_claims(PAGE)["SP_955_finite"]
    assert abs(selectivity(J, p0) - x0) <= 1e-4 * abs(x0)
