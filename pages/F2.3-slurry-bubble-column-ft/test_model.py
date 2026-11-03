"""Fast checks of model.py for page F2.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
from model import Operating, SlurryColumn, reactor_ode  # noqa: E402

PAGE = HERE.name


def test_balance_closes_and_species_agree():
    """Feed minus face exit equals the liquid uptake; CO and H2 convert alike."""
    sol = SlurryColumn(Operating(0.40, 0.30), 200).solve()
    assert np.max(np.abs(sol.conversion("face") - sol.conversion("reaction"))) < 1e-12
    assert abs(np.diff(sol.conversion())[0]) < 1e-12
    assert sol.resid < 1e-10


def test_jacobian_is_exact():
    """The analytic Jacobian against central differences on a small grid."""
    col = SlurryColumn(Operating(0.40, 0.30), 10)
    u = np.concatenate([np.linspace(1.0, 0.5, 20), [0.1, 0.2]])
    jac = col.residual(u)[1].toarray()
    fd = np.empty_like(jac)
    for k in range(u.size):
        e = np.zeros_like(u)
        e[k] = 1e-6
        fd[:, k] = (col.residual(u + e)[0] - col.residual(u - e)[0]) / 2e-6
    assert np.max(np.abs(jac - fd)) < 1e-7 * np.max(np.abs(fd))


def test_first_order_to_the_ode_route():
    """The finite-volume conversion halves its gap to the ODE route per refinement."""
    op = Operating(0.40, 0.30)
    ref = reactor_ode(op)["syngas_conversion"]
    e1, e2 = (abs(SlurryColumn(op, n).solve().syngas_conversion() - ref) for n in (200, 400))
    assert 1.8 < e1 / e2 < 2.2


def test_stored_conversion_reproduced():
    """The committed U = 0.40 conversion, held to its own rtol."""
    stored = load_claims(PAGE)
    x = SlurryColumn(Operating(0.40, 0.30), 1600).solve().syngas_conversion()
    assert x == pytest.approx(stored["conversion_U040"], rel=1e-4)
