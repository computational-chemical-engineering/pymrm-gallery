"""Fast checks of model.py for page C2.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import ATOMS, STOICH, SteamReformer, XuFroment, reference_profile  # noqa: E402

PAGE = HERE.name


def _kin():
    return XuFroment.from_table(load_data("xu-froment-1989-parameters.csv", page=PAGE))


def test_solve_against_radau():
    """One reactor at 848 K: Newton converges and meets the Radau route at first order."""
    kin = _kin()
    m = SteamReformer(kin, 848.0).solve()
    taus = np.array([0.15, 0.30, 0.40])
    ref = reference_profile(kin, 848.0, taus)
    assert m.backward_error() < 1e-8
    assert np.abs(np.stack([m.at(taus, 0), m.at(taus, 1)], axis=-1) - ref).max() < 1e-4
    assert np.abs(ATOMS @ (STOICH @ m.u.T)).max() < 1e-12
    assert np.allclose(m.outlet(), m.u[-1])


def test_stored_claim_reproduced():
    """The committed mean absolute deviation, held to the rtol of its claim (0.05)."""
    kin = _kin()
    obs = load_data("xu-froment-1989-conversion.csv", page=PAGE)
    ms = {T: SteamReformer(kin, T).solve() for T in sorted(obs.temperature_K.unique())}
    dev = [ms[r.temperature_K].at(r.W_F, 0 if r.quantity == "x_CH4" else 1) - r.value
           for r in obs.itertuples()]
    x0 = load_claims(PAGE)["mad_conversion"]
    assert abs(np.mean(np.abs(dev)) - x0) <= 0.05 * abs(x0)
