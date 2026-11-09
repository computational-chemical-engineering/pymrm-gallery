"""Fast checks of model.py for page A3.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import BoundaryLayer, higbie, leveque, sh_sphere  # noqa: E402

PAGE = HERE.name


def test_sphere_floor_and_flux_identity():
    """The stagnant sphere is within 1e-3 of 2/(1 - ratio) on 200 cells, and the
    surface and shell fluxes agree to the solve's residual."""
    sh = sh_sphere(1e-4, 200)
    assert abs(sh * (1 - 1e-4) / 2.0 - 1) < 1e-3
    assert abs(sh / sh_sphere(1e-4, 200, at="shell") - 1) < 1e-9


def test_boundary_layer_limits():
    """Both limits on a small grid: within 1 % of their closed forms, and the flux
    sum equal to the exit-plane balance."""
    for flow, exact in ((dict(u_s=0.0, gdot=1.0), leveque(1e-9, 1.0, 1.0)),
                        (dict(u_s=1e-3, gdot=0.0), higbie(1e-9, 1e-3, 1.0))):
        r = BoundaryLayer(D=1e-9, L=1.0, n_y=100, **flow).run(100)
        assert abs(r["kbar"] / exact - 1) < 1e-2
        assert abs(r["kbar"] / r["kbar_balance"] - 1) < 1e-9


def test_stored_claim_reproduced():
    """The committed leveque_rel_err_raw, held to the rtol of its claim."""
    stored = load_claims(PAGE)["leveque_rel_err_raw"]
    k = BoundaryLayer(D=1e-9, u_s=0.0, gdot=1.0, L=1.0, n_y=800).run(800)["kbar"]
    err = k / leveque(1e-9, 1.0, 1.0) - 1
    assert abs(err - stored) <= 1e-4 * abs(stored)
