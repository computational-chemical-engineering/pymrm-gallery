"""Fast checks of model.py for page D2.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import (BarkelewTube, Case, blowup_z, linearisations, p0_crit_vwf,  # noqa: E402
                   t_max_peak, tangency, tau_max_phase, tau_max_pymrm)

PAGE = HERE.name


def _case():
    par = load_data("van-welsenaere-froment-1970-parameters.csv", page=PAGE)
    return Case.from_table(zip(par.symbol, par.value))


def test_parameter_reading():
    """T_M (Eq. 8) and the criterion-1 back-integration against the printed Example 1."""
    case = _case()
    assert abs(t_max_peak(case, 625.0) - 656.6) < 0.05
    assert abs(p0_crit_vwf(case, 625.0) - 0.01651) / 0.01651 < 3.1e-4


def test_stored_claim_reproduced():
    """The committed tangency hot spot at S = 4, held to its claim's rtol (1e-5)."""
    x0 = load_claims(PAGE)["tau_m_tangency_S4"]
    assert abs(tangency(4.0)[1] - x0) <= 1e-5 * abs(x0)


def test_pymrm_tube_against_quadrature():
    """The finite-volume tube converges at first order onto the quadrature hot spot."""
    N, S = 0.944 * 4.0, 4.0
    ref = tau_max_phase(N, S)
    e1 = abs(tau_max_pymrm(N, S, n_z=200) - ref)
    e2 = abs(tau_max_pymrm(N, S, n_z=400) - ref)
    assert e2 < 0.01 * ref
    assert 1.8 < e1 / e2 < 2.2


def test_failed_solve_raises():
    """A solve that does not converge raises instead of returning an iterate."""
    tube = BarkelewTube(4.0, 4.0, 1.0, 200, cool=-1.0)
    with pytest.raises((RuntimeError, ValueError)):
        tube.solve()


def test_steep_exponent_ignites_at_inlet():
    """The square-dropped linearisation runs away within the first millimetre."""
    case = _case()
    kw, slope = linearisations(case, 625.0)["square dropped, a(T-T_w)/T_w"]
    z, Tmax = blowup_z(case, 625.0, kw, slope, 0.017)
    assert z is not None and z < 1e-3 and Tmax == pytest.approx(800.0)
    with pytest.raises(RuntimeError):
        blowup_z(case, 625.0, kw, slope, 0.017, route="z")
    kw, slope = linearisations(case, 625.0)["tangent at T_w  (this page)"]
    z_phase, _ = blowup_z(case, 625.0, kw, slope, 0.017)
    z_march, _ = blowup_z(case, 625.0, kw, slope, 0.017, route="z")
    assert np.isclose(z_phase, z_march, atol=1e-6)
