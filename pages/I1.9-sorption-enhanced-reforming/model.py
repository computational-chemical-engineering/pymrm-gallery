"""Batch equilibrium of steam reforming and the water-gas shift with a fraction f of
the CO2 removed (Hufton, Mayorga and Sircar 1999, Eqs. 3 to 10, p. 250).

Model code of gallery page I1.9. The equilibrium constants are not part of the
paper; every function takes them as ``kfun(T_K) -> (K_SMR in atm^2, K_WGS)``.

Contents
--------
psig_to_atm
    The paper's pressure conversion, kPa = (psig + 14.7) x 6.89 (Table 1 footnote).
Case
    One condition: temperature, pressure, feed ratio, 1 - f and the K functions.
extent_residual, solve_extents
    Route A: the paper's extents (beta, delta) with Eqs. (8) and (9), solved by
    ``pymrm.newton`` in the log unknowns (ln beta, ln t), t = (beta - delta)/beta.
State
    Mole fractions (Eq. 7), the two conversions (Eq. 10 and the carbon-based one),
    dry fractions, element balances and the Eq. (5)/(6) identities.
species_residual, solve_species
    Route B: five gas amounts with the C, H and O balances (sorbed CO2 included)
    and Eqs. (3) and (4), also by ``pymrm.newton``.
unique_root
    ``pymrm.checks.find_roots`` that raises unless exactly one root is found.
langmuir_theta, serp_closure
    Closed forms of printed numbers (p. 253, p. 254 to 255).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from pymrm import newton
from pymrm.checks import find_roots

__all__ = ["KPA_PER_ATM", "psig_to_atm", "Case", "State", "extent_residual", "solve_extents",
           "species_residual", "solve_species", "unique_root", "langmuir_theta",
           "serp_closure", "half_unit_interval", "with_om", "Branch"]

KPA_PER_ATM = 101.325


def psig_to_atm(psig, kpa_per_psi=6.89, p_atm_psi=14.7):
    """Absolute pressure (atm) from gauge psig by the paper's footnote conversion."""
    return (psig + p_atm_psi) * kpa_per_psi / KPA_PER_ATM


@dataclass(frozen=True)
class Case:
    """One batch-equilibrium condition.

    Parameters
    ----------
    T_C : float
        Temperature (degC).
    P_atm : float
        Total gas pressure (atm, absolute).
    om : float
        1 - f, the fraction of the CO2 made that stays in the gas (1 means no sorbent).
    kfun : callable
        ``kfun(T_K) -> (K_SMR, K_WGS)`` with K_SMR in atm^2.
    alpha : float, optional
        Feed H2O/CH4 (6 in the paper).
    """

    T_C: float
    P_atm: float
    om: float
    kfun: object
    alpha: float = 6.0

    def K(self):
        """``(K_SMR, K_WGS)`` at this temperature."""
        return self.kfun(self.T_C + 273.15)


def _terms(case, beta, delta, s, t):
    """Linear factors of Eqs. (8) and (9): value and d/d(beta), d/d(delta).

    The two factors that vanish as f -> 1 are taken from the unknowns themselves,
    1 - beta (1 + alpha) = s and beta - delta = beta t, so no cancellation enters."""
    a, f = case.alpha, 1.0 - case.om
    return {"bmd": (beta * t, 1.0, -1.0), "h2": (3 * beta + delta, 3.0, 1.0),
            "tot": (1 + 2 * beta - f * delta, 2.0, -f),
            "ch4": (s, -(1 + a), 0.0),
            "h2o": (a - (1 + a) * (beta + delta), -(1 + a), -(1 + a)), "d": (delta, 0.0, 1.0)}


# exponent of each factor in ln(RHS of Eq. 8) and ln(RHS of Eq. 9)
_EQ8 = {"bmd": 1, "h2": 3, "tot": -2, "ch4": -1, "h2o": -1}
_EQ9 = {"d": 1, "h2": 1, "bmd": -1, "h2o": -1}


def extent_residual(u, case, jac_slip=0.0):
    """Residual and Jacobian of Eqs. (8) and (9) in ``u = (ln s, ln t)``.

    ``s = 1 - beta (1 + alpha)`` is the unconverted CH4 fraction and
    ``t = (beta - delta)/beta`` the unshifted CO fraction; both go to zero as f -> 1,
    so their logs keep the iteration inside the feasible set.

    Each residual is ``ln K - ln(right-hand side)``, so it is dimensionless and of
    order one. ``jac_slip`` adds a wrong term to the Jacobian (a break-row knob; 0).
    """
    a, P = case.alpha, case.P_atm
    s, t = np.exp(u[0]), np.exp(u[1])
    beta = (1.0 - s) / (1 + a)
    delta = beta * (1.0 - t)
    tm = _terms(case, beta, delta, s, t)
    k_smr, k_wgs = case.K()
    g = np.empty(2)
    jac = np.zeros((2, 2))
    dbeta = np.array([-s / (1 + a), 0.0])                 # d beta / d u
    ddelta = np.array([-s / (1 + a) * (1 - t), -beta * t])  # d delta / d u
    for i, (expo, const) in enumerate(((_EQ8, 2 * np.log(P) + 2 * np.log(1 + a)),
                                       (_EQ9, np.log(case.om) + np.log(1 + a)))):
        lhs = np.log(k_smr) if i == 0 else np.log(k_wgs)
        g[i] = lhs - const - sum(e * np.log(tm[k][0]) for k, e in expo.items())
        for k, e in expo.items():
            v, db, dd = tm[k]
            jac[i] -= e * (db * dbeta + dd * ddelta) / v
    jac[0, 0] += jac_slip
    return g, jac


def solve_extents(case, guess=None, tol=1e-10, maxfev=60, jac_slip=0.0, backward_tol=1e-10):
    """Route A: ``(beta, delta)`` at ``case`` by ``pymrm.newton``.

    ``tol`` bounds the Newton step in the log unknowns and ``backward_tol`` the
    backward error; both must hold.

    Raises
    ------
    RuntimeError
        When newton does not report success.
    """
    if guess is None:
        guess = np.log([0.5, 0.5])
        if case.om < 0.3:   # follow the branch down from om = 1, three steps a decade
            for om in np.geomspace(1.0, case.om, int(np.ceil(-3 * np.log10(case.om))) + 1)[:-1]:
                guess = solve_extents(replace(case, om=om), guess, tol, maxfev, jac_slip,
                                      backward_tol).u
    with np.errstate(invalid="ignore", divide="ignore"):   # a trial step outside the
        res = newton(extent_residual, np.asarray(guess, float),   # feasible set is nan,
                     args=(case, jac_slip), tol=tol, maxfev=maxfev,  # which newton damps
                     backward_tol=backward_tol)
    if not res.success:
        raise RuntimeError(f"newton failed at {case}: {res.message}")
    beta = (1.0 - np.exp(res.x[0])) / (1 + case.alpha)
    return State(case, float(beta), float(beta * (1 - np.exp(res.x[1]))), res.x.copy(),
                 int(getattr(res, "nit", -1)))


@dataclass(frozen=True)
class State:
    """An equilibrium state on the paper's extents.

    Attributes
    ----------
    case : Case
    beta, delta : float
        Moles of CO made by SMR and shifted by WGS, per mole of feed.
    u : numpy.ndarray
        The solved log unknowns, a warm start for a neighbouring case.
    nit : int
        Newton iterations.
    """

    case: Case
    beta: float
    delta: float
    u: np.ndarray
    nit: int = -1

    @property
    def y(self):
        """Wet mole fractions of Eq. (7): CH4, H2O, CO, CO2, H2."""
        a, f, b, d = self.case.alpha, 1 - self.case.om, self.beta, self.delta
        tot = 1 + 2 * b - f * d
        s_, t_ = np.exp(self.u)     # 1 - beta (1 + alpha) and (beta - delta)/beta, exact
        return {"CH4": s_ / ((1 + a) * tot),
                "H2O": (a - (b + d) * (1 + a)) / ((1 + a) * tot),
                "CO": b * t_ / tot, "CO2": self.case.om * d / tot, "H2": (3 * b + d) / tot}

    @property
    def ydry(self):
        """Dry mole fractions (water removed)."""
        y = self.y
        return {k: v / (1 - y["H2O"]) for k, v in y.items() if k != "H2O"}

    @property
    def X(self):
        """Carbon-based CH4 conversion, beta (1 + alpha)."""
        return self.beta * (1 + self.case.alpha)

    @property
    def chi4(self):
        """The paper's 'CH4 to H2 conversion', chi/4 = (3 beta + delta)(1 + alpha)/4 (Eq. 10)."""
        return (3 * self.beta + self.delta) * (1 + self.case.alpha) / 4

    @property
    def cox_ppm(self):
        """Dry CO + CO2 in ppm."""
        d = self.ydry
        return (d["CO"] + d["CO2"]) * 1e6

    def element_gaps(self):
        """C, H, O balances per mole of feed (sorbed CO2 included), minus the feed."""
        a, f, b, d = self.case.alpha, 1 - self.case.om, self.beta, self.delta
        n = {"CH4": 1 / (1 + a) - b, "H2O": a / (1 + a) - b - d, "CO": b - d,
             "CO2": (1 - f) * d, "H2": 3 * b + d}
        s = f * d
        return np.array([n["CH4"] + n["CO"] + n["CO2"] + s - 1 / (1 + a),
                         4 * n["CH4"] + 2 * n["H2O"] + 2 * n["H2"] - (4 + 2 * a) / (1 + a),
                         n["H2O"] + n["CO"] + 2 * (n["CO2"] + s) - a / (1 + a)])

    def eq5_eq6(self):
        """Relative gaps of Eq. (5), Eq. (6) corrected (y_H2O) and printed Eq. (6) (y_H2)
        against the state's own fractions; the printed one is a ratio, not a gap."""
        y, P = self.y, self.case.P_atm
        ks, kw = self.case.K()
        e5 = (ks * kw) ** 0.25 / P ** 0.5 * y["CH4"] ** 0.25 * y["H2O"] ** 0.5 / y["CO2"] ** 0.25
        common = ks ** 0.25 / (kw ** 0.75 * P ** 0.5) * y["CH4"] ** 0.25 * y["CO2"] ** 0.75
        return (e5 / y["H2"] - 1, common / y["H2O"] ** 0.5 / y["CO"] - 1,
                common / y["H2"] ** 0.5 / y["CO"])


_SP = ("CH4", "H2O", "CO", "CO2", "H2")
# element content of each gas species: C, H, O
_ATOMS = np.array([[1, 4, 0], [0, 2, 1], [1, 0, 1], [1, 0, 2], [0, 2, 0]], float)


def species_residual(v, case, o_sorbed=2.0):
    """Route B: residual and Jacobian in ``v = ln n`` (gas moles per mole feed).

    Rows: C, H, O balances (sorbed CO2 = (1 - om)/om x gas CO2, with ``o_sorbed``
    oxygens each; 2 is right, the knob is a break row), each relative to the feed's
    content, then ln K_SMR and ln K_WGS minus Eqs. (3) and (4).
    """
    n = np.exp(v)
    a, om = case.alpha, case.om
    feed = np.array([1.0, 4 + 2 * a, a]) / (1 + a)
    srb = (1 - om) / om
    sorbed = np.array([srb, 0.0, o_sorbed * srb])
    atoms = _ATOMS.copy()
    atoms[3] += sorbed
    g = np.empty(5)
    jac = np.zeros((5, 5))
    g[:3] = (atoms.T @ n) / feed - 1
    jac[:3] = atoms.T * n / feed[:, None]
    N = n.sum()
    ks, kw = case.K()
    c_smr = np.array([-1, -1, 1, 0, 3.0])
    c_wgs = np.array([0, -1, -1, 1, 1.0])
    g[3] = np.log(ks) - 2 * np.log(case.P_atm) - c_smr @ v + 2 * np.log(N)
    jac[3] = -c_smr + 2 * n / N
    g[4] = np.log(kw) - c_wgs @ v
    jac[4] = -c_wgs
    return g, jac


def solve_species(case, guess, tol=1e-12, o_sorbed=2.0):
    """Route B solve from ``guess`` (gas moles); returns the wet mole fractions."""
    res = newton(species_residual, np.log(np.asarray(guess, float)), args=(case, o_sorbed),
                 tol=tol, maxfev=80, backward_tol=1e-10)
    if not res.success:
        raise RuntimeError(f"route B newton failed at {case}: {res.message}")
    n = np.exp(res.x)
    return dict(zip(_SP, n / n.sum()))


def unique_root(fun, lo, hi, log=False, n_scan=48, xtol=1e-12):
    """The single root of ``fun`` on ``[lo, hi]`` by ``pymrm.checks.find_roots``.

    Raises
    ------
    ValueError
        Unless exactly one root is found and every scan point evaluated.
    """
    out = find_roots(fun, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(out["roots"]) != 1 or out["failed"]:
        raise ValueError(f"find_roots on [{lo:g}, {hi:g}]: {out['message']}")
    return out["roots"][0]


def langmuir_theta(p, b):
    """Langmuir coverage b p / (1 + b p), p in atm, b in atm^-1 (p. 253)."""
    return b * p / (1 + b * p)


def serp_closure(X, h_net=None, h_press=0.0):
    """Dry H2 fraction of a product holding all unconverted CH4 and no carbon oxides.

    Per mole CH4 fed: 4X H2 and (1 - X) CH4, so y_H2 = 4X/(1 + 3X). With ``h_net``
    (net H2, mol/kg) and ``h_press`` (pressurisation H2 counted into the product,
    mol/kg), the CH4 fed is h_net/(4X) per kg.
    """
    if h_net is None:
        return 4 * X / (1 + 3 * X)
    fed = h_net / (4 * X)
    h2 = h_net + h_press
    return h2 / (h2 + (1 - X) * fed)


def half_unit_interval(printed):
    """``(lo, hi)`` of the values that round to the printed string."""
    s = printed.strip()
    dec = len(s.split(".")[1]) if "." in s else 0
    h = 0.5 * 10.0 ** (-dec)
    v = float(s)
    return v - h, v + h


def with_om(case, om):
    """A copy of ``case`` at another 1 - f."""
    return replace(case, om=om)


class Branch:
    """Route A solves along one branch, each warm-started from the last state.

    Root scans and sweeps call ``solve`` at neighbouring conditions; the last
    solution is the guess, and the continuation from 1 - f = 1 of
    :func:`solve_extents` is the fallback when that guess does not converge.
    """

    def __init__(self, **solve_kwargs):
        self.kw = solve_kwargs
        self.u = None

    def solve(self, case):
        """The state at ``case``."""
        if self.u is not None:
            try:
                st = solve_extents(case, self.u, **self.kw)
                self.u = st.u
                return st
            except RuntimeError:
                pass
        st = solve_extents(case, None, **self.kw)
        self.u = st.u
        return st
