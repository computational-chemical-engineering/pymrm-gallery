"""Fast checks of model.py for page J5.3 (run with GALLERY_OFFLINE=1, PYTHONPATH=shared)."""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "shared"))

from gallery_utils import load_claims, load_data  # noqa: E402
import model  # noqa: E402

PAGE = HERE.name


def _condition(TC=50, I0=0.0258):
    t1 = load_data("chiu1983-table1-rate-constants.csv", page=PAGE)
    t2 = load_data("chiu1983-table2-model-parameters.csv", page=PAGE)
    r = t2[(t2.polym_temp_C == TC) & np.isclose(t2.initiator_loading_mol_per_L, I0)].iloc[0]
    return model.Chiu(TC, I0, r.theta_t_min, r.theta_p_min, r.A, r.B,
                      model.table1_constants(t1))


def test_qssa_root_solves_eq12():
    """The quasi-steady lambda_0 zeroes eq. 12's right-hand side."""
    m = _condition()
    lam = m.lam0_qssa(0.3, 0.02)
    y = np.array([0.02, 0.3, lam, 0.0, 0.0, 0.0, 0.0, 0.0])
    assert abs(m.rhs(0.0, y)[2]) < 1e-12 * 2.0 * m.f * m.kd * 0.02


def test_shell_bvp_matches_stored_order_and_flux():
    """The spherical shell flux converges on K_exact; the stored rms error is reproduced."""
    b = model.shell_bvp(1280, 0.05)
    assert abs(b["K_mid"] / b["K_exact"] - 1.0) < 1e-4
    stored = load_claims(PAGE)["bvp_rms_err_finest"]
    assert abs(b["err"] - stored) <= 1e-6 * stored


def test_moment_identity_holds():
    """Eq. 12-17 reproduce the summed chain-length equations to round-off."""
    assert model.moment_identity(ktc=4.0e8)[0] < 1e-14
    assert abs(model.moment_identity(ktc=4.0e8, break_term="m2")[0] - 2.0) < 1e-12
