"""Fast checks of model.py for page A3.5 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import (nu_shell, sh_closed, sh_stefan, stefan_threshold,  # noqa: E402
                   stefan_threshold_lambertw)

PAGE = HERE.name


def test_stagnant_sphere_floor():
    """The stagnant shell gives 2/(1 - a/b) to within 1e-3 on 400 cells (error 1.3e-4,
    falling by 4 per doubling), and the surface and shell fluxes agree to 1e-9."""
    ratio = 1.0e4
    assert abs(nu_shell(400, ratio) * (1 - 1 / ratio) / 2.0 - 1) < 1e-3
    assert abs(nu_shell(400, ratio) / nu_shell(400, ratio, at="shell") - 1) < 1e-9


def test_threshold_two_routes():
    """The root and the Lambert-W closed form of the 5 % threshold agree to 1e-12."""
    assert abs(stefan_threshold(0.05) - stefan_threshold_lambertw(0.05)) < 1e-12


def test_stored_claim_reproduced():
    """The committed stefan_pymrm_vs_closed_rel, to rtol 1e-9."""
    stored = load_claims(PAGE)["stefan_pymrm_vs_closed_rel"]
    ys = [0.01, 0.05, 0.10, 0.20]
    worst = max(abs(sh_stefan(600, 1.0e5, ya)[0] / sh_closed(1.0e5, ya) - 1) for ya in ys)
    assert abs(worst - stored) <= 1e-9 * abs(stored)
