"""Fast checks of model.py for page A1.7 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import DavidsonBubble, cloud_radius_exact, rise_velocity  # noqa: E402

PAGE = HERE.name
PAGE_E12 = "E1.2-davidson-bubble"

# Geldart's own constants, CGS, and a 25 cm bubble in 2.6 g/cm3 sand
G, MU, C3, EPS_0, RHO_F, D_B = 981.0, 1.8e-4, 8.0e-4, 0.4, 1.2e-3, 25.0


def cloud(d_um, n=400):
    """R_c/R of a D_B bubble in sand of size d_um: u_f = U_0 / eps_0, U_0 from eq. (3)."""
    u_f = C3 * G * (d_um * 1e-4) ** 2 * (2.6 - RHO_F) / MU / EPS_0
    u_br = rise_velocity(D_B, g=G)
    sol = DavidsonBubble(0.5 * D_B, u_f, u_br, n=n).solve()
    return sol.cloud_radius() / sol.R, cloud_radius_exact(u_f, u_br)


def test_model_is_e12_unchanged():
    """The model is E1.2's model.py, character for character."""
    e12 = HERE.parent / PAGE_E12 / "model.py"
    if not e12.is_file():
        pytest.skip("E1.2 not in this checkout")
    assert (HERE / "model.py").read_bytes() == e12.read_bytes()


def test_cloud_matches_closed_form_and_vanishes_in_group_d():
    """A cloud in the group B sand agrees with E1.2's closed form; the group D sand has none."""
    rc, exact = cloud(600.0)
    assert rc == pytest.approx(exact, rel=1e-4)
    assert np.isnan(cloud(700.0)[0]) and np.isnan(cloud(700.0)[1])


def test_stored_claim_reproduced():
    """The committed R_c/R of the 600 um sand at n = 1600, held to its own rtol."""
    stored = load_claims(PAGE)["cloud_RcR_sand_600um"]
    assert abs(cloud(600.0, n=1600)[0] / stored - 1) <= 1e-6
