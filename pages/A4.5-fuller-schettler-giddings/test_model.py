"""Fast checks of model.py for page A4.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import D_fsg, Tube, flux_quadrature, fsg_params  # noqa: E402

PAGE = HERE.name


def test_tube_against_quadrature():
    """The pymrm tube at 400 cells reproduces the stored flux deviation from the
    quadrature (flux_pymrm_vs_quadrature), rtol 1e-6, and conserves the flux."""
    par = fsg_params(load_data("fuller-1966-table1-diffusion-volumes.csv", page=PAGE))
    T_mid = float(np.sqrt(300.0 * 600.0))
    case = dict(L=10.0, T0=300.0, T1=600.0, D_ref=D_fsg("AIR", "BENZENE", T_mid, par),
                T_ref=T_mid, b=1.75)
    sol = Tube(400, **case).solve()
    dev = abs(sol["flux"] / flux_quadrature(**case)[0] - 1)
    stored = load_claims(PAGE)["flux_pymrm_vs_quadrature"]
    assert abs(dev - stored) <= 1e-6 * stored
    assert sol["spread"] <= 1e-10 * abs(sol["flux"])


def test_reproduction_of_printed_errors():
    """Eq. (4) reproduces the stored largest deviation from the printed errors
    (max_abs_delta_pp) over the rows other than the two named ones."""
    par = fsg_params(load_data("fuller-1966-table1-diffusion-volumes.csv", page=PAGE))
    d = load_data("fuller-1966-table3-diffusion-data.csv", page=PAGE)
    calc = np.array([D_fsg(a, b, T, par) for a, b, T in
                     zip(d.component_a, d.component_b, d.temperature_K)])
    delta = np.abs(100 * (calc - d.D_obs_cm2_s) / d.D_obs_cm2_s - d.err_fsg_pct).values
    named = (((d.component_a == "HE") & (d.component_b == "WATER") & (d.temperature_K == 352.5))
             | ((d.component_a == "CO2") & (d.component_b == "ETHYLENEOXIDE")
                & (d.temperature_K == 298.0))).values
    stored = load_claims(PAGE)["max_abs_delta_pp"]
    assert abs(delta[~named].max() - stored) <= 1e-9 * stored
