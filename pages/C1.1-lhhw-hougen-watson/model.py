"""Hougen and Watson (1947): the codimer LHHW estimation and the fitted law in a bed.

Model code of gallery page C1.1. The tables are passed in as arrays or data
frames; nothing here loads data or plots.

Contents
--------
MECH_LIN, TRANSFORMS, MECH_ZEROS
    The book's linearisation of each of the 18 mechanisms of Table B
    (transformed observable, hydrogen basis) and which constants the theory
    requires to vanish.
MECH_RATE
    The same 18 mechanisms as rate laws with positive parameters.
lin_fit
    The book's method: unweighted linear least squares on the transformed
    variable.
fit_rate
    A rate-space fit with relative residuals, positive parameters through a log
    reparameterisation, Levenberg-Marquardt from fixed starts, with the solver
    status, the gradient at the solution and per-parameter standard errors.
eqq_constants, rate_eqq, rate_max_temperature
    Eq. (q) and the temperature of its rate maximum at a fixed composition, a
    root of the analytic temperature derivative of ``ln r``.
IsothermalBed, bed_wf, wf_quadrature
    Steady isothermal plug flow in the catalyst-mass coordinate ``w = W/F_U0``
    with pymrm (upwind convection of the three molar flows, ``newton``), the
    ``w`` at which the outlet conversion reaches a target (a root), and the
    same quantity by adaptive quadrature of ``dX / r``.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

from pymrm import NumJac, compute_boundary_values, construct_convflux_upwind, construct_div, newton
from pymrm.checks import find_roots

SQ = np.sqrt
RGAS = 1.987  # cal/(g-mol K), the book's value

# mechanism -> (R-transform key, H-basis key); Table B, book pp. 947-949
MECH_LIN = {
    "a": ("pH/r", "pH"), "b": ("pU/r", "pH"),
    "c": ("pHpU/r", "pH"), "d": ("sqrt(pHpU/r)", "pH"),
    "e": ("sqrt(pH/r)", "sqrt"), "f": ("pU/r", "sqrt"),
    "g": ("pHpU/r", "sqrt"), "h": ("cbrt(pHpU/r)", "sqrt"),
    "i": ("pH/r", "pH"), "j": ("pHpU/r", "pH"),
    "k": ("pHpU/r", "pH"), "l": ("sqrt(pH/r)", "sqrt"),
    "m": ("pHpU/r", "sqrt"), "n": ("sqrt(pHpU/r)", "sqrt"),
    "o": ("pHpU/r", "pH"), "p": ("pHpU/r", "pH"),
    "q": ("pU/r", "pH"), "r": ("pHpU/r", "pH"),
}
TRANSFORMS = {
    "pH/r": lambda H, U, r: H / r,
    "pU/r": lambda H, U, r: U / r,
    "pHpU/r": lambda H, U, r: H * U / r,
    "sqrt(pHpU/r)": lambda H, U, r: np.sqrt(H * U / r),
    "sqrt(pH/r)": lambda H, U, r: np.sqrt(H / r),
    "cbrt(pHpU/r)": lambda H, U, r: np.cbrt(H * U / r),
}
# the constants each mechanism requires to vanish (the rest must be >= 0),
# from the structure of Table B's equations
MECH_ZEROS = {"a": "f", "b": "b", "c": "c", "d": "", "e": "f", "f": "b",
              "g": "c", "h": "", "i": "bf", "j": "bc", "k": "b", "l": "bf",
              "m": "bc", "n": "b", "o": "f", "p": "cf", "q": "bf", "r": "bcf"}

# rate laws with positive parameters p (lumped prefactor first): (n_par, law)
MECH_RATE = {
    "a": (3, lambda p, H, U, S: p[0] * H / (1 + p[1] * U + p[2] * S)),
    "b": (3, lambda p, H, U, S: p[0] * U / (1 + p[1] * H + p[2] * S)),
    "c": (3, lambda p, H, U, S: p[0] * H * U / (1 + p[1] * H + p[2] * U)),
    "d": (4, lambda p, H, U, S: p[0] * H * U / (1 + p[1] * H + p[2] * U + p[3] * S) ** 2),
    "e": (3, lambda p, H, U, S: p[0] * H / (1 + p[1] * U + p[2] * S) ** 2),
    "f": (3, lambda p, H, U, S: p[0] * U / (1 + SQ(p[1] * H) + p[2] * S)),
    "g": (3, lambda p, H, U, S: p[0] * H * U / (1 + SQ(p[1] * H) + p[2] * U)),
    "h": (4, lambda p, H, U, S: p[0] * H * U / (1 + SQ(p[1] * H) + p[2] * U + p[3] * S) ** 3),
    "i": (2, lambda p, H, U, S: p[0] * H / (1 + p[1] * S)),
    "j": (2, lambda p, H, U, S: p[0] * H * U / (1 + p[1] * H)),
    "k": (3, lambda p, H, U, S: p[0] * H * U / (1 + p[1] * H + p[2] * S)),
    "l": (2, lambda p, H, U, S: p[0] * H / (1 + p[1] * S) ** 2),
    "m": (2, lambda p, H, U, S: p[0] * H * U / (1 + SQ(p[1] * H))),
    "n": (3, lambda p, H, U, S: p[0] * H * U / (1 + SQ(p[1] * H) + p[2] * S) ** 2),
    "o": (3, lambda p, H, U, S: p[0] * H * U / (1 + p[1] * U + p[2] * S)),
    "p": (2, lambda p, H, U, S: p[0] * H * U / (1 + p[1] * U)),
    "q": (2, lambda p, H, U, S: p[0] * U / (1 + p[1] * S)),
    "r": (1, lambda p, H, U, S: p[0] * H * U),
}
# mechanism (d) without the product's adsorption term K_S p_S
RATE_D_NO_KS = (3, lambda p, H, U, S: p[0] * H * U / (1 + p[1] * H + p[2] * U) ** 2)


def lin_fit(mech, H, U, S, r, basis=None):
    """Unweighted least squares on the transformed variable (the book's method).

    Parameters
    ----------
    mech : str
        Mechanism letter, ``"a"`` to ``"r"``.
    H, U, S, r : ndarray
        Partial pressures of hydrogen, codimer and saturate (atm) and the rate.
    basis : {None, "pH", "sqrt"}
        Hydrogen basis of the linearisation; ``None`` takes the mechanism's own.

    Returns
    -------
    coef : ndarray
        The constants ``a, b, c, f`` of ``R = a + b pU + c pS + f h(pH)``.
    X : ndarray
        The design matrix, columns ``1, pU, pS, h(pH)``.
    R : ndarray
        The transformed observable.
    """
    tk, bk = MECH_LIN[mech]
    bk = bk if basis is None else basis
    R = TRANSFORMS[tk](H, U, r)
    h = H if bk == "pH" else SQ(H)
    X = np.column_stack([np.ones_like(H), U, S, h])
    coef = np.linalg.lstsq(X, R, rcond=None)[0]
    return coef, X, R


def fit_rate(law, H, U, S, r, max_nfev=40000, starts=(0.0, -1.0, 1.0)):
    """Fit a rate law to rates by relative residuals, parameters kept positive.

    Parameters
    ----------
    law : tuple
        ``(n_par, f(p, H, U, S))``, an entry of ``MECH_RATE``.
    H, U, S, r : ndarray
        Partial pressures (atm) and the rates the law is fitted to.
    max_nfev : int
        Function-evaluation cap of each Levenberg-Marquardt run.
    starts : sequence of float
        Each start sets every log-parameter to this value.

    Returns
    -------
    dict
        ``p`` (parameters), ``pct`` (per-run ``100 |r_fit - r| / r``),
        ``status`` (least_squares status of the best start, positive when
        converged), ``grad`` (infinity norm of the gradient of the cost in the
        log-parameters at the solution), ``spread`` (largest relative cost
        difference between the starts), ``se_log`` (standard error of each
        ``ln p``, from ``s^2 (J^T J)^-1`` with ``s^2 = sum(res^2) / dof``;
        ``inf`` where ``J^T J`` is singular), ``n`` and ``dof``.
    """
    npar, fun = law

    def res(lp):
        return (fun(np.exp(lp), H, U, S) - r) / r

    sols = [least_squares(res, np.full(npar, x0), method="lm", max_nfev=max_nfev) for x0 in starts]
    best = min(sols, key=lambda s: s.cost)
    costs = np.array([s.cost for s in sols])
    n = r.size
    dof = n - npar
    jtj = best.jac.T @ best.jac
    s2 = 2.0 * best.cost / dof if dof > 0 else np.nan
    try:
        cond = np.linalg.cond(jtj)
        se = np.sqrt(np.diag(np.linalg.inv(jtj)) * s2) if cond < 1e14 else np.full(npar, np.inf)
    except np.linalg.LinAlgError:
        se = np.full(npar, np.inf)
    return {"p": np.exp(best.x), "pct": 100 * np.abs(best.fun), "status": int(best.status),
            "grad": float(np.max(np.abs(best.jac.T @ best.fun))),
            "spread": float(np.max(costs - best.cost) / max(best.cost, 1e-300)),
            "se_log": se, "n": n, "dof": dof}


def eqq_constants(T_K, thermo, dS_S=None):
    """``alpha, K_H, K_U, K_S`` of eq. (q) at ``T_K``.

    Parameters
    ----------
    T_K : float or ndarray
        Temperature, K.
    thermo : mapping
        ``dH_H, dH_U, dH_S, dS_H, dS_U, dS_S, A, B`` as printed (cal/g-mol and
        cal/(g-mol K)); ``ln K = -dH/(R T) + dS/R``, ``ln alpha = -A/(R T) + B/R``.
    dS_S : float, optional
        Replaces ``thermo["dS_S"]``.
    """
    dss = thermo["dS_S"] if dS_S is None else dS_S

    def K(dH, dS):
        return np.exp(-dH / (RGAS * T_K) + dS / RGAS)

    alpha = np.exp(-thermo["A"] / (RGAS * T_K) + thermo["B"] / RGAS)
    return alpha, K(thermo["dH_H"], thermo["dS_H"]), K(thermo["dH_U"], thermo["dS_U"]), K(thermo["dH_S"], dss)


def rate_eqq(T_K, pH, pU, pS, thermo, dS_S=None, expo=2):
    """Forward rate of eq. (q), lb-mol/(lb hr); ``expo`` is the adsorption-group exponent."""
    al, KH, KU, KS = eqq_constants(T_K, thermo, dS_S)
    return al * KH * KU * pH * pU / (1 + KH * pH + KU * pU + KS * pS) ** expo


def rate_max_temperature(thermo, pH, pU, pS, lo, hi, expo=2, xtol=1e-9):
    """Temperature of the rate maximum of eq. (q) at a fixed composition.

    The maximum is the root of ``R T^2 d ln r / dT = (A + dH_H + dH_U)
    - expo * sum_i K_i p_i dH_i / (1 + sum_i K_i p_i)``, found by
    ``find_roots`` on ``[lo, hi]``.

    Returns
    -------
    float
        The root, K.

    Raises
    ------
    ValueError
        Unless exactly one root lies on ``[lo, hi]``, or when it is not a
        maximum (the derivative must fall through zero).
    """
    num = thermo["A"] + thermo["dH_H"] + thermo["dH_U"]

    def g(T):
        _, KH, KU, KS = eqq_constants(T, thermo)
        terms = (KH * pH, KU * pU, KS * pS)
        dH = (thermo["dH_H"], thermo["dH_U"], thermo["dH_S"])
        return num - expo * sum(t * h for t, h in zip(terms, dH)) / (1 + sum(terms))

    found = find_roots(g, lo, hi, n_scan=64, xtol=xtol)
    if len(found["roots"]) != 1:
        raise ValueError(f"expected one root of d ln r/dT on [{lo}, {hi}] K: {found['message']}")
    T = found["roots"][0]
    if not g(lo) > 0 > g(hi):
        raise ValueError("the root is not a maximum of the rate")
    return T


def _partials(N, P):
    """Partial pressures (pH, pU, pS) from molar flows ordered (U, H2, S)."""
    tot = N.sum(axis=-1, keepdims=True)
    p = P * N / tot
    return p[..., 1], p[..., 0], p[..., 2]


class IsothermalBed:
    """Steady isothermal plug flow of U + H2 -> S in the catalyst-mass coordinate.

    ``dN_i/dw = nu_i r(p)`` on ``0 < w < w_end``, ``w = W/F_U0``, for the molar
    flows ``N = (N_U, N_H2, N_S)`` per mole of codimer fed, shape ``(n, 3)``
    (space first, species last). The flows are convected at unit velocity in
    ``w`` by first-order upwind fluxes; the rate couples the species pointwise,
    so its Jacobian comes from ``NumJac((n, 3))``.

    Parameters
    ----------
    rate : callable
        ``rate(pH, pU, pS)``, lb-mol/(lb hr).
    n : int
        Number of cells.
    w_end : float
        Bed length in ``w``, lb catalyst per (lb-mol U/hr).
    feed : sequence of float
        Inlet molar flows (U, H2, S) per mole of U.
    nu : sequence of float
        Stoichiometric coefficients (U, H2, S).
    P : float
        Total pressure, atm.
    moving_total : bool
        Partial pressures from the local total flow (True) or from the inlet
        total (False).
    nu_geom : int
        Geometry factor of ``construct_div`` (0 for the mass coordinate).
    outlet_bc : dict, optional
        Outlet boundary condition; default ``{"outflow": True}``.
    tol : float
        Newton step tolerance.
    """

    def __init__(self, rate, n, w_end, feed=(1.0, 1.4, 0.0), nu=(-1.0, -1.0, 1.0), P=3.5,
                 moving_total=True, nu_geom=0, outlet_bc=None, tol=1e-11):
        self.rate, self.n, self.P = rate, n, P
        self.feed = np.asarray(feed, dtype=float)
        self.nu = np.asarray(nu, dtype=float)
        self.moving_total = moving_total
        self.tol = tol
        self.w_f = np.linspace(0.0, w_end, n + 1)
        self.w_c = 0.5 * (self.w_f[:-1] + self.w_f[1:])
        # outward normal n, a dN/dn + b N = d:
        #   w = 0 (inlet, n = -w): N_i = N_i,feed, the face carries the feed flow -> a = 0, b = 1
        #   w = w_end (outlet): pure convection, the face carries the adjacent cell's flow
        self.bc = ({"a": 0.0, "b": 1.0, "d": self.feed[None, :]},
                   {"outflow": True} if outlet_bc is None else outlet_bc)
        shape = (n, 3)
        conv, conv_bc = construct_convflux_upwind(shape, self.w_f, self.w_c, self.bc, v=1.0)
        div = construct_div(shape, self.w_f, nu=nu_geom)
        self.L = div @ conv                       # assembled once
        self.b = np.asarray((div @ conv_bc).todense()).reshape(shape)
        self.numjac = NumJac(shape)
        self.N = np.tile(self.feed, (n, 1))
        self.result = None

    def source(self, N):
        """``nu_i r`` per cell, shape ``(n, 3)``."""
        if self.moving_total:
            pH, pU, pS = _partials(N, self.P)
        else:
            p = self.P * N / self.feed.sum()
            pH, pU, pS = p[:, 1], p[:, 0], p[:, 2]
        return self.nu[None, :] * self.rate(pH, pU, pS)[:, None]

    def residual(self, x):
        N = x.reshape(self.N.shape)
        s, ds = self.numjac(self.source, N)
        g = (self.L @ x.ravel() + self.b.ravel()) - s.ravel()
        return g, self.L - ds

    def solve(self):
        """Solve with ``newton``; raises ``RuntimeError`` unless it reports success."""
        res = newton(self.residual, self.N.ravel(), tol=self.tol, maxfev=50, stagnation=True)
        if not res.success:
            raise RuntimeError(f"bed solve failed: {res.message}")
        self.result = res
        self.N = res.x.reshape(self.N.shape)
        return self

    def outlet(self):
        """Molar flows on the outlet face, read with ``compute_boundary_values``."""
        _, _, out, _ = compute_boundary_values(self.N, self.w_f, self.w_c, self.bc)
        return np.ravel(out)

    def conversion_out(self):
        """Codimer conversion on the outlet face."""
        return 1.0 - self.outlet()[0] / self.feed[0]


def bed_wf(rate, n, X_target, lo, hi, n_scan=8, xtol=1e-9, **kw):
    """Bed length ``w`` at which the outlet conversion of an ``n``-cell bed reaches ``X_target``.

    The root of ``X_out(w_end) - X_target`` on ``[lo, hi]`` by ``find_roots``;
    raises ``ValueError`` unless exactly one root is found.
    """
    def f(w_end):
        return IsothermalBed(rate, n, w_end, **kw).solve().conversion_out() - X_target

    found = find_roots(f, lo, hi, n_scan=n_scan, xtol=xtol)
    if len(found["roots"]) != 1:
        raise ValueError(f"expected one root of X_out - {X_target} on [{lo}, {hi}]: {found['message']}")
    return found["roots"][0]


def wf_quadrature(rate, X_target, feed=(1.0, 1.4, 0.0), nu=(-1.0, -1.0, 1.0), P=3.5,
                  moving_total=True, epsrel=1e-12):
    """``w = int_0^X dX / r`` by adaptive quadrature, with the extent ``X`` per mole of U fed.

    Returns
    -------
    (float, float)
        The integral and quad's error estimate.
    """
    feed = np.asarray(feed, dtype=float)
    nu = np.asarray(nu, dtype=float)

    def inv_rate(X):
        N = feed + nu * X
        tot = N.sum() if moving_total else feed.sum()
        pU, pH, pS = P * N / tot
        return 1.0 / rate(pH, pU, pS)

    val, err = quad(inv_rate, 0.0, X_target, epsabs=0.0, epsrel=epsrel, limit=200)
    return val, err
