"""Fast checks of model.py for page A3.8 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from pymrm.checks import find_roots  # noqa: E402
from model import (HR, Absorber, OndaSet, column_closed_form,  # noqa: E402
                   set_printed_constants)

PAGE = HERE.name

# The page's illustrative operating point (index.py), air and water at 25 C on spheres.
K = load_data("onda-1968-correlation-constants.csv", page=PAGE)
SR = load_data("onda-1968-stated-results.csv", page=PAGE).set_index("quantity")
PROPS = dict(rho_L=997.0, mu_L=8.90e-4 * HR, D_L=1.92e-9 * HR, sigma=0.0720 * HR**2,
             rho_G=1.184, mu_G=1.84e-5 * HR, D_G=1.60e-5 * HR, T=298.15, P=1.0)
D_P, G_REF, M, Y_IN, X_IN, REMOVAL = 0.0254, 2000.0, 1.0, 0.02, 0.0, 0.95
L_REF = float(SR.loc["L_fig6_hi_kg_m2_hr", "value"])
A_T = float(SR.loc["atDp_spheres", "value"]) / D_P
C_AV, G_M, L_M = PROPS["rho_L"] / 18.015, G_REF / 28.96, L_REF / 18.015
NAMES = ("C_aw", "e_sigma", "e_Re_L", "e_Fr_L", "e_We_L", "C_kL", "e_Re_Lw", "e_Sc_L",
         "e_atDp_L", "C_kG", "e_Re_G", "e_Sc_G", "e_atDp_G")


def _onda():
    set_printed_constants(**{n: float(K.loc[K.symbol == n, "value"].iloc[0]) for n in NAMES})
    return OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=0.85, **PROPS)


def _stored(name):
    return load_claims(PAGE)[name]


def _height(o):
    KGa = o.KGa(L_REF, G_REF, M, C_AV)
    f = lambda Z: column_closed_form(Z, KGa, M, G_M, L_M, Y_IN, X_IN)[0] - (1 - REMOVAL) * Y_IN
    roots = find_roots(f, 1e-4, 500.0, log=True, xtol=1e-14)["roots"]
    assert len(roots) == 1
    return roots[0], KGa


def test_stored_metrics_reproduced():
    """The gas-side resistance share and the packed height equal the committed
    metrics f_G_ref and Z_star_m to 1e-9 relative."""
    o = _onda()
    f_G = o.f_G(L_REF, G_REF, M, C_AV)
    Z, _ = _height(o)
    assert abs(f_G / _stored("f_G_ref") - 1) < 1e-9
    assert abs(Z / _stored("Z_star_m") - 1) < 1e-9


def test_column_matches_closed_form():
    """The pymrm column agrees with the closed form of the same BVP to 1e-5 at 800
    cells, the solute balance closes to round-off, the deferred-correction route
    reaches the same solution, and the upwind solve is first order."""
    o = _onda()
    Z, KGa = _height(o)
    col = Absorber(Z, KGa, M, G_M, L_M, Y_IN, X_IN, n=800).solve()
    y_ex, x_ex = column_closed_form(Z, KGa, M, G_M, L_M, Y_IN, X_IN)
    assert abs(col.y_out / y_ex - 1) < 1e-5
    assert abs(col.x_out / x_ex - 1) < 1e-5
    assert col.imbalance() < 1e-12
    up = Absorber(Z, KGa, M, G_M, L_M, Y_IN, X_IN, n=800).solve(scheme="upwind")
    assert abs(up.y_out / y_ex - 1) > 1e-3        # first order: far coarser than van Leer
    dc = Absorber(Z, KGa, M, G_M, L_M, Y_IN, X_IN, n=800).solve(route="deferred")
    assert abs(dc.y_out / col.y_out - 1) < 1e-10
