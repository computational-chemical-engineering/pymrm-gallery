"""Krishna et al. (1997) slurry holdup correlations and a two-class gas reactor.

Correlations (paper, eqs. 1 to 7) and the page's reactor consequence (closed forms
and a pymrm finite-volume solve). Every function takes the parameter dictionary
``p`` built from the paper's printed constants, so a defect in one constant is one
changed entry (keys listed in :data:`KEYS`). No plotting, no data loading.
"""

import numpy as np
from scipy import sparse
from scipy.special import lambertw
from pymrm import (compute_boundary_values, construct_convflux_upwind, construct_div,
                   construct_grad, newton)
from pymrm.checks import find_roots

#: Entries of the parameter dictionary and their meaning (all SI).
KEYS = {
    "V0": "V_small,0 of eq. 2 [m/s]",
    "cV": "slope of V_small in eps_s, eq. 2 [m/s]",
    "edf0": "eps_df,0 of eq. 4 [-]",
    "cE": "slope of eps_df in eps_s, eq. 4 [-]",
    "a2_low": "alpha2 for eps_s <= es_low [m^0.18 (m/s)^-0.58]",
    "a2_high": "alpha2 for eps_s >= es_high",
    "es_low": "upper end of the low-slurry range of alpha2",
    "es_high": "lower end of the high-slurry range of alpha2",
    "N": "diameter exponent of eq. 6",
    "m": "velocity exponent of eq. 6",
    "D_ref": "reference diameter of Fig. 8 [m]",
}


def v_small(eps_s, p):
    """Small-bubble rise velocity, eq. 2: V_small,0 + 0.8 eps_s [m/s]."""
    return p["V0"] + p["cV"] * np.asarray(eps_s, float)


def eps_df(eps_s, p):
    """Dense-phase voidage, eq. 4: eps_df,0 (1 - (0.7/eps_df,0) eps_s) [-]."""
    return p["edf0"] - p["cE"] * np.asarray(eps_s, float)


def u_df(eps_s, p, form="paper"):
    """Dense-phase superficial velocity at the regime transition.

    ``form="paper"`` is eq. 5, U_df = V_small eps_df. ``form="reilly"`` is the
    U = V_small eps (1 - eps) convention of Reilly and of Krishna and Ellenberger
    (1996), a break row (misprints and ambiguities M4).
    """
    e = eps_df(eps_s, p)
    if form == "paper":
        return v_small(eps_s, p) * e
    if form == "reilly":
        return v_small(eps_s, p) * e * (1.0 - e)
    raise ValueError(form)


def alpha2(eps_s, p):
    """alpha2 of eq. 6: 0.268 up to eps_s 0.10, 0.3 from 0.16; raises in between (M1)."""
    if eps_s <= p["es_low"] + 1e-12:
        return p["a2_low"]
    if eps_s >= p["es_high"] - 1e-12:
        return p["a2_high"]
    raise ValueError(f"alpha2 is not defined for {p['es_low']} < eps_s = {eps_s} < {p['es_high']}")


def eps_b(U, D_T, eps_s, p, udf_form="paper"):
    """Dilute-phase (large-bubble) holdup, eq. 6, for U > U_df [-]."""
    x = np.asarray(U, float) - u_df(eps_s, p, udf_form)
    if np.any(x <= 0):
        raise ValueError("eq. 6 needs U > U_df")
    return alpha2(eps_s, p) * D_T ** (-p["N"]) * x ** p["m"]


def eps_total(U, D_T, eps_s, p, overlap=True, udf_form="paper"):
    """Total gas holdup: eq. 3 (U <= U_df) or eq. 7 (U > U_df) [-].

    ``overlap=False`` drops the (1 - eps_b) factor of eq. 7 (a break row).
    """
    U = np.asarray(U, float)
    udf = u_df(eps_s, p, udf_form)
    out = U / v_small(eps_s, p)
    het = U > udf
    if np.any(het):
        eb = eps_b(U[het], D_T, eps_s, p, udf_form)
        ed = eps_df(eps_s, p)
        out = np.array(out, float, copy=True)
        out[het] = eb + ed * ((1.0 - eb) if overlap else 1.0)
    return out


def eps_b_from_total(eps, e_df):
    """Eq. 1, the dynamic-gas-disengagement inversion: (eps - eps_df)/(1 - eps_df)."""
    return (eps - e_df) / (1.0 - e_df)


def fig8_amplitude(eps_s, p, D_ref=None):
    """Amplitude of the Fig. 8 line, alpha2 D_ref^-N (E8)."""
    D_ref = p["D_ref"] if D_ref is None else D_ref
    return alpha2(eps_s, p) * D_ref ** (-p["N"])


def udf_maximum(p):
    """eps_s and U_df at the maximum of U_df(eps_s) = (V0 + cV eps_s)(edf0 - cE eps_s)."""
    es = (p["cV"] * p["edf0"] - p["cE"] * p["V0"]) / (2.0 * p["cV"] * p["cE"])
    return es, float(u_df(es, p))


# ---------------------------------------------------------------- reactor (page model)

def split(U, D_T, eps_s, p, H0=None, overlap_in_tau=True, udf_form="paper", edf=None):
    """Flow and holdup shares of the two gas classes (page equations R1 to R3).

    Returns a dict with ``eps``, ``eps_b``, ``eps_df``, ``U_df``, ``beta`` (flow share of
    the large bubbles), ``phi`` (dense-phase share of the gas holdup), ``tb`` and ``tdf``
    (residence times over the mean), ``r = tdf/tb`` and, when ``H0`` is given,
    ``H = H0/(1 - eps)`` and ``tau_mean = eps H/U`` [s].
    ``overlap_in_tau=False`` takes the dense gas per column volume as eps_df instead of
    eps_df (1 - eps_b) in tau_df (a break row; tau_b keeps eps_b/eps). ``edf`` replaces eq. 4 by a measured dense-phase
    voidage (eqs. 5 and 6 then use it too).
    """
    if edf is None:
        udf = float(u_df(eps_s, p, udf_form))
        eb = float(eps_b(U, D_T, eps_s, p, udf_form))
        ed = float(eps_df(eps_s, p))
    else:
        ed = float(edf)
        udf = float(v_small(eps_s, p)) * ed * ((1.0 - ed) if udf_form == "reilly" else 1.0)
        eb = float(alpha2(eps_s, p) * D_T ** (-p["N"]) * (U - udf) ** p["m"])
    eps = eb + ed * (1.0 - eb)
    beta = 1.0 - udf / U
    phi = ed * ((1.0 - eb) if overlap_in_tau else 1.0) / eps
    out = dict(eps=eps, eps_b=eb, eps_df=ed, U_df=udf, beta=beta, phi=phi,
               tb=(eb / eps) / beta, tdf=phi / (1.0 - beta))     # eb/eps = 1 - phi in eq. 7
    out["r"] = out["tdf"] / out["tb"]
    if H0 is not None:
        out["H"] = H0 / (1.0 - eps)
        out["tau_mean"] = eps * out["H"] / U
    return out


def x_dense(Da, mode, Pe=None):
    """Conversion of the dense-phase gas: stirred, plug or axially dispersed (Danckwerts)."""
    Da = np.asarray(Da, float)
    if mode == "stirred":
        return Da / (1.0 + Da)
    if mode == "plug":
        return -np.expm1(-Da)
    if mode == "dispersed":
        q = np.sqrt(1.0 + 4.0 * Da / Pe)
        return 1.0 - 4.0 * q * np.exp(Pe * (1.0 - q) / 2.0) / ((1.0 + q) ** 2 - (1.0 - q) ** 2 * np.exp(-q * Pe))
    raise ValueError(mode)


def conversion(D, s, mode="stirred", Pe=None, k_ratio=1.0):
    """Overall gas conversion at D = k tau_mean (R4); large bubbles in plug flow.

    ``k_ratio`` is k_df/k_b (1 in the base case).
    """
    xb = -np.expm1(-D * s["tb"])
    return s["beta"] * xb + (1.0 - s["beta"]) * x_dense(k_ratio * D * s["tdf"], mode, Pe)


def d_required(s, X=0.9, mode="stirred", Pe=None, k_ratio=1.0):
    """The one D = k tau_mean giving conversion ``X`` (R6), by ``find_roots``."""
    res = find_roots(lambda D: conversion(D, s, mode, Pe, k_ratio) - X, 1e-3, 1e3, n_scan=200, log=True)
    if len(res["roots"]) != 1:
        raise RuntimeError(f"expected one design root, found {res['roots']}: {res['message']}")
    return res["roots"][0]


def kH0_over_U(D, s):
    """k H0/U for a given D = k tau_mean: D (1 - eps)/eps."""
    return D * (1.0 - s["eps"]) / s["eps"]


def crossing(s, route="lambertw"):
    """Per-gas crossing a = k tau_b where 1 - exp(-a) = r a/(1 + r a) (R7), r = tdf/tb > 1."""
    r = s["r"]
    if route == "lambertw":
        return float(np.real(-lambertw(-np.exp(-1.0 / r) / r, -1) - 1.0 / r))
    res = find_roots(lambda a: -np.expm1(-a) - r * a / (1.0 + r * a), 1e-3, 1e3, n_scan=400, log=True)
    if len(res["roots"]) != 1:
        raise RuntimeError(f"expected one crossing, found {res['roots']}")
    return res["roots"][0]


class TwoClassReactor:
    """Two gas classes in parallel on x = z/H, assembled from pymrm operators.

    Large bubbles: plug flow, dc_b/dx + D tb c_b = 0, first-order upwind with an
    outflow outlet. Dense-phase gas: c' - c''/Pe + D tdf c = 0 with a Danckwerts
    inlet and a zero-gradient outlet. The system is linear and solved by one
    ``newton`` step. Each field has its own ``(n, 1)`` operators.

    Parameters
    ----------
    n : int
        Cells.
    Pe : float
        Dense-phase axial Peclet number.
    Da_b, Da_df : float
        D tb and D tdf.
    beta : float
        Flow share of the large bubbles.
    b_outlet : {"outflow", "zero_gradient"}
        Large-bubble outlet (zero gradient is a break row, P6).
    inlet_sign : float
        Sign of ``a`` in the dense-phase inlet dictionary (-1 is a break row, P1).
    """

    def __init__(self, n, Pe, Da_b, Da_df, beta, b_outlet="outflow", inlet_sign=1.0):
        self.n, self.Pe, self.Da_b, self.Da_df, self.beta = n, Pe, Da_b, Da_df, beta
        self.x_f = np.linspace(0.0, 1.0, n + 1)
        self.x_c = 0.5 * (self.x_f[1:] + self.x_f[:-1])
        # large bubbles, x = 0: c_b = 1; x = 1: convective outflow (no condition imposed)
        out = {"outflow": True} if b_outlet == "outflow" else {"a": 1.0, "b": 0.0, "d": 0.0}
        self.bc_b = ({"a": 0.0, "b": 1.0, "d": 1.0}, out)
        # dense phase, x = 0 (outward normal -x): c - c'/Pe = 1 (Danckwerts);
        # x = 1: c' = 0
        self.bc_d = ({"a": inlet_sign / Pe, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        shape = (n, 1)
        div = construct_div(shape, self.x_f, nu=0)          # nu = 0: Cartesian
        cb, cb_bc = construct_convflux_upwind(shape, self.x_f, self.x_c, self.bc_b, v=1.0)
        cd, cd_bc = construct_convflux_upwind(shape, self.x_f, self.x_c, self.bc_d, v=1.0)
        gd, gd_bc = construct_grad(shape, self.x_f, self.x_c, self.bc_d)
        eye = sparse.eye_array(n, format="csc")
        self.jac = sparse.block_array([[div @ cb + Da_b * eye, None],
                                       [None, div @ (cd - gd / Pe) + Da_df * eye]], format="csc")
        self.g = np.concatenate([np.ravel(_dense(div @ cb_bc)),
                                 np.ravel(_dense(div @ (cd_bc - gd_bc / Pe)))])

    def residual(self, u):
        """Residual and (constant) Jacobian, for ``newton``."""
        return self.jac @ u + self.g, self.jac

    def solve(self, lin_solver_kwargs=None):
        """One newton step from zero (the system is linear); raises unless ``success``."""
        sol = newton(self.residual, np.zeros(2 * self.n), tol=1e-9, maxfev=3,
                     lin_solver_kwargs=lin_solver_kwargs)
        if not sol.success:
            raise RuntimeError(f"two-class solve failed: {sol.message}")
        self.u = np.asarray(sol.x)
        return self

    def exits(self):
        """Exit face values (c_b, c_df) by ``compute_boundary_values`` under the same bc."""
        cb = self.u[:self.n].reshape(self.n, 1)
        cd = self.u[self.n:].reshape(self.n, 1)
        _, _, b1, _ = compute_boundary_values(cb, self.x_f, self.x_c, self.bc_b)
        _, _, d1, _ = compute_boundary_values(cd, self.x_f, self.x_c, self.bc_d)
        return float(np.ravel(b1)[0]), float(np.ravel(d1)[0])

    def conversion(self):
        """beta (1 - c_b,out) + (1 - beta)(1 - c_df,out)."""
        b1, d1 = self.exits()
        return self.beta * (1.0 - b1) + (1.0 - self.beta) * (1.0 - d1)


def _dense(m):
    """A boundary vector from a pymrm constructor as a dense array."""
    return np.asarray(m.todense()) if sparse.issparse(m) else np.asarray(m)
