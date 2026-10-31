"""Fast checks of model.py for page I1.2 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
from model import I_TS, Feed, Monolith, OhCavendish  # noqa: E402

PAGE = HERE.name


def _inputs():
    par = load_data("oh-cavendish-1982-parameters.csv", page=PAGE)
    vals = {(r.symbol, r.case): float(r.value) for r in par.itertuples()}

    def value(sym, case="both"):
        return vals.get((sym, case), vals.get((sym, "both")))
    return (OhCavendish.from_table(value), Feed.from_table(value, "case 1"),
            Feed.from_table(value, "case 2"))


def test_case2_peak_reproduces_stored_claim():
    """The refined peak wall temperature on 40 cells, held to its claim's rtol (1e-4)."""
    prm, _, feed2 = _inputs()
    m = Monolith(prm, feed2, n_x=40)
    mr = m.march(14.0, 0.05, record=lambda s: s.u.copy())
    k = int(np.argmax([u[:, I_TS].max() for u in mr.records]))
    assert 0 < k < len(mr.records) - 1
    _, T_max, _ = m.peak(mr.records[k - 1], 0.05)
    T0 = load_claims(PAGE)["case2_T_max_n40_dt0.05"]
    assert abs(T_max - T0) <= 1e-4 * abs(T0)
    assert np.all(mr.nit >= 1)


def test_outflow_face_is_last_cell():
    """With the outflow outlet the transported face value is the last cell's."""
    prm, feed1, _ = _inputs()
    m = Monolith(prm, feed1, n_x=20)
    m.step(0.1)
    assert np.array_equal(m.outlet(), m.u[-1])


def test_lightoff_root_inside_its_step():
    """The light-off event is root-found inside the step that crosses 50 % conversion."""
    prm, feed1, _ = _inputs()
    m = Monolith(prm, feed1, n_x=20)
    ev = lambda u: 1.0 - m.outlet(u)[2] / m.y_in[0] - 0.5   # noqa: E731
    mr = m.march(60.0, 0.5, record=lambda s: 1.0 - s.outlet()[2] / s.y_in[0], events=[ev])
    t_lo = mr.t_event[0]
    X = np.array(mr.records)
    i = int(np.flatnonzero(X >= 0.5)[0])
    assert mr.t[i - 1] < t_lo <= mr.t[i]
