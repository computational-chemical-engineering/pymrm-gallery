"""Fast checks of model.py for page A4.1 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import SMFilm, exact_solve, set_table  # noqa: E402

PAGE = HERE.name


def test_pymrm_film_against_closed_form():
    """The pymrm film at n_z = 320 reproduces the stored worst flux error against
    the closed-form solution (pymrm_vs_closed_form_reldiff_n320), rtol 1e-9."""
    base = load_data("wilke-1950-table1.csv", page=PAGE).drop_duplicates("example")
    set_table(base.set_index("example"))
    err = 0.0
    for e in base.example:
        exact = exact_solve(e)[0]
        N = SMFilm(e, n_z=320).solve()
        err = max(err, max(abs(N[i] / exact[i] - 1) for i in (0, 1)))
    stored = load_claims(PAGE)["pymrm_vs_closed_form_reldiff_n320"]
    assert abs(err - stored) <= 1e-9 * abs(stored)
