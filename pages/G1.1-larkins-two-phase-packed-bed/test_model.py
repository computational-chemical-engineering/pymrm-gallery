"""Fast checks of model.py for page G1.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (Column, ColumnFV, Larkins, backward_euler_march, quad_outlet,  # noqa: E402
                   recipe_outlet, recover_mu_air, threshold_length)

PAGE = HERE.name
T3 = load_data("larkins1961-table3-measured.csv", page=PAGE)
T4 = load_data("larkins1961-table4-calculated.csv", page=PAGE)
K = load_data("larkins1961-printed-constants.csv", page=PAGE)
KV = dict(zip(K.quantity, K.value.astype(float)))
STORED = load_claims(PAGE)


def _model():
    m = Larkins.from_constants(KV)
    mg = T3.air_rate_lb_sqft_min.values > 0
    s, i, _, _ = recover_mu_air(m, T3.air_rate_lb_sqft_min.values[mg],
                                T3.column_temperature_F.values[mg], T4.reynolds_air.values[mg])
    return replace(m, mu_air_line=(s, i))


def _column(m, run):
    r = T3[T3.run == run].iloc[0]
    return Column(m, float(r.liquid_rate_lb_sqft_min), float(r.air_rate_lb_sqft_min),
                  float(r.column_temperature_F), float(r.liquid_viscosity_cP))


P_IN = KV["equipment_pressure_limit"] + Larkins.from_constants(KV).p_atm


def test_fv_outlet_error_matches_stored_claim():
    """The outflow solve on 1600 cells against the quadrature, to the claim's rtol."""
    col = _column(_model(), 63)
    fv = ColumnFV(col, P_IN, 4.0, 1600)
    P = fv.solve()
    err = abs(fv.outlet_value(P) / quad_outlet(col, P_IN, 4.0) - 1.0)
    x0 = STORED["fv_rel_error_at_n_fine"]
    assert abs(err - x0) <= 1e-3 * abs(x0)
    assert abs(fv.balance_gap(P)) < 1e-10


def test_outflow_is_the_backward_euler_march():
    """pymrm with {"outflow": True} equals the hand-written march to round-off."""
    col = _column(_model(), 63)
    fv = ColumnFV(col, P_IN, 4.0, 100)
    P = fv.solve()
    assert abs(fv.outlet_value(P) / backward_euler_march(col, P_IN, 4.0, 100) - 1.0) < 1e-10


def test_one_phase_horizontal_identity():
    """Arithmetic-mean recipe is exact for one phase without gravity."""
    m = replace(_model(), gravity=False, two_phase=False)
    r = T3[T3.run == 160].iloc[0]
    cg = Column(m, 0.0, float(r.air_rate_lb_sqft_min), float(r.column_temperature_F),
                float(r.liquid_viscosity_cP))
    Pe = quad_outlet(cg, P_IN, 1.0)
    assert abs((recipe_outlet(cg, P_IN, 1.0) - Pe) / (P_IN - Pe)) < 1e-11


def test_one_percent_length_matches_stored_claim():
    """The root-found 1 % bed length of run 63."""
    L1, _ = threshold_length(_column(_model(), 63), P_IN, 0.01)
    x0 = STORED["two_phase_recipe_1pct_length_ft"]
    assert abs(L1 - x0) <= 1e-3 * abs(x0)
