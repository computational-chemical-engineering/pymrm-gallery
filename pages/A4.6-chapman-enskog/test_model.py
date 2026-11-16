"""Fast checks of model.py for page A4.6 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def test_printed_limits_and_stored_value():
    """Section 14.3's Lorentz limit: [D12]_2/[D12]_1 = 13/12 for rigid spheres,
    equal to the stored lorentz_D2_over_D1."""
    d1 = model.delta_generic(1.0, 1.0, 0.0, 1.0, 1.0, 0.4, 0.6, 0.2)
    assert abs(1.0 / (1.0 - d1) - 13.0 / 12.0) < 1e-14
    stored = load_claims(PAGE)["lorentz_D2_over_D1"]
    assert abs(1.0 / (1.0 - d1) - stored) <= 1e-6 * stored


def test_turning_point_is_a_root():
    """The outermost turning point satisfies the radial equation to round-off."""
    e = np.array([0.1, 1.0, 10.0])
    b = np.array([1.5, 1.0, 0.5])
    r0, n_pos = model.turning_point(e, b)
    f = 1.0 - (b / r0) ** 2 - 4.0 * (r0 ** -12 - r0 ** -6) / e
    assert np.all(np.abs(f) < 1e-10)
    assert set(n_pos.tolist()) <= {1, 3}


def test_constant_d_cell_against_series():
    """A constant-D cell (one solve, 100 cells, 200 steps) against the cosine series."""
    cell = model.DiffusionCell(50.0, 100, 0.6)
    s = cell.split(cell.solve(300.0, 200))
    assert abs(s / model.analytic_split(0.6, 50.0, 300.0) - 1.0) < 5e-3


def test_variable_d_cell_one_newton_run():
    """A composition-dependent cell solves (Newton reports success on every step)."""
    cell = model.DiffusionCell(50.0, 40, lambda x: 0.6 * (1.0 + 0.03 * x), tol=1e-12)
    x = cell.solve(300.0, 20)
    assert np.all((x > -1e-12) & (x < 1.0 + 1e-12))
