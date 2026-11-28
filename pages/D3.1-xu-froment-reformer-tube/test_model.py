"""Fast checks of model.py for page D3.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
import model as M  # noqa: E402

PAGE = HERE.name
C21 = "C2.1-xu-froment-smr"
C25 = "C2.5-vandenbussche-froment-methanol-co2-route"
A45 = "A4.5-fuller-schettler-giddings"


def _module(page):
    spec = importlib.util.spec_from_file_location(page.split("-")[0].replace(".", "_") + "_model",
                                                  HERE.parent / page / "model.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _case():
    c21, a45 = _module(C21), _module(A45)
    kin = c21.XuFroment.from_table(load_data("xu-froment-1989-parameters.csv", page=C21))
    sh = load_data("nist-webbook-thermo.csv", page=C25)
    sh = sh[~((sh.species == "N2") & (sh.T_min_K < 400))].set_index("species")
    ch4 = load_data("nist-webbook-ch4-shomate.csv", page=PAGE).set_index("species")
    shom = np.array([ch4.loc["CH4", list("ABCDEF")].to_numpy(float)]
                    + [sh.loc[s, list("ABCDEF")].to_numpy(float) for s in ("H2O", "H2", "CO", "CO2", "N2")])
    fv = load_data("fuller-1966-table1-diffusion-volumes.csv", page=A45).set_index("entry").diffusion_volume
    vol = np.array([fv["C"] + 4 * fv["H"], fv["H2O"], fv["H2"], fv["CO"], fv["CO2"], fv["N2"]])
    aw = a45.AW
    mw = np.array([aw["C"] + 4 * aw["H"], 2 * aw["H"] + aw["O"], 2 * aw["H"], aw["C"] + aw["O"],
                   aw["C"] + 2 * aw["O"], 2 * aw["N"]])
    pn = load_data("xu-froment-1989-ii-printed-numbers.csv", page=PAGE, dtype=str).set_index("quantity")
    dv = tuple(float(pn.loc[k, "value"]) for k in pn.index if k.startswith("pore_d_"))
    return M.Case(kin=kin, shomate=shom, mw=mw, vol=vol,
                  pore_d=(30.0, 40.0, 65.0, 100.0, 160.0, 350.0, 750.0, 5500.0, 55000.0), pore_dv=dv), c21


def test_reuse_identity():
    """V5: the page's partial-pressure rates equal C2.1's XuFroment.rates."""
    case, c21 = _case()
    p = c21.partial_pressures(0.25, 0.2, 28.0, 3.358, 0.122)
    assert np.allclose(case.rates(p, 900.0), case.kin.rates(0.25, 0.2, 900.0, 28.0, 3.358, 0.122),
                       rtol=1e-13, atol=0.0)


def test_tanh_closed_form():
    """The graded slab against tanh(phi)/phi at second order."""
    e = [M.thiele_slab_eta(10.0, n)[0] for n in (160, 320)]
    assert abs((e[1] + (e[1] - e[0]) / 3) / (math.tanh(10.0) / 10.0) - 1) < 1e-8


def test_pellet_check_row_and_claim():
    """One C check row (eta_CH4, ring basis) and the stored cut ratio claim."""
    case, _ = _case()
    pb, yf = case.state(0.25, 0.20, 28.0)
    De = case.d_eff(yf, 900.0, 28.0)
    r = case.rates(pb, 900.0)
    ext = M.pellet_rates(case, pb, 900.0, De, n=80, levels=2)
    cv = load_data("check_values.csv", page=PAGE).set_index("name")
    ref = cv.loc["C_eta_CH4_ring"]
    assert abs((ext[0] + ext[2]) / (r[0] + r[2]) / case.ring_factor / ref.value - 1) <= ref.rtol
    cut = M.pellet_rates(case, pb, 900.0, De, n=80, levels=2, cut=0.08e-3)
    x0 = load_claims(PAGE)["C_robs_CH4_cut_over_full"]
    assert abs((cut[0] + cut[2]) / (ext[0] + ext[2]) / x0 - 1) < 1e-6


def test_energy_invariant_short_tube():
    """V1 on a 0.5 m tube: the enthalpy flow gains exactly the march's heat sum."""
    case, _ = _case()
    tb = M.Tube(case, lambda z: 29.0, q=lambda z: 70.0e3, n_r=40)
    Y, Qs = tb.march(np.linspace(0.0, 0.5, 11))
    gap = case.h_flow(*Y[-1]) - case.h_flow(*Y[0]) - Qs[-1]
    assert abs(gap) / Qs[-1] < 1e-12
    assert 0.0 < Y[-1, 0] < 1.0
