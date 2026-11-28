"""Saez and Carbonell (1985): trickle-flow holdup and pressure drop from relative permeabilities.

Closed forms and one coupled algebraic pair, no discretisation. Every published
constant sits in one dictionary ``P`` that every function takes, so a misread
constant is one changed entry (refresh checklist rule 34). Operating conditions
are a dictionary ``op`` of SI scalars or arrays (``d_e`` m, ``eps`` -, ``L`` and
``G`` kg/(m2 s), ``rho_L``, ``rho_G`` kg/m3, ``mu_L``, ``mu_G`` Pa s, ``sigma``
N/m, ``g`` m/s2).

Symbols follow the paper with L for its beta (liquid) and G for its gamma (gas).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq

from pymrm import BlockLayout, newton
from pymrm.checks import find_roots

#: The constants as the paper prints them (journal pages in comments).
PRINTED = dict(
    A=180.0,       # Eq. 50 companion, p. 60: viscous Ergun constant (Macdonald et al.)
    B=1.8,         # p. 60: inertial Ergun constant
    a=2.43,        # Eq. 32, p. 56: k_L = delta^a
    c=4.8,         # Eq. 48, p. 58: k_G = S_G^c
    e34=0.41,      # Eq. 34, p. 57: delta = X_L^0.41, printed for 1/a
    eo_a=20.0,     # Eq. 8, p. 54: eps0 = 1/(20 + 0.9 Eo*)
    eo_b=0.9,      # Eq. 8
)


def params(**changes):
    """The printed constants with ``changes`` applied, plus the reading switches.

    Switches: ``re_mu`` (True: Re* carries the viscosity, as Ergun requires;
    False: the nomenclature's printed form without it), ``exact`` (True: the
    gas weight kept, Eq. 21; False: Eq. 51 as printed), ``stagnant_exp``
    (None: 1/a; a number: that exponent in the stagnant-gas closed form).
    """
    P = dict(PRINTED, re_mu=True, exact=False, stagnant_exp=None)
    unknown = set(changes) - set(P)
    if unknown:
        raise KeyError(f"unknown parameter(s) {sorted(unknown)}")
    P.update(changes)
    return P


# ---------------------------------------------------------------- groups

def reynolds(flux, mu, op, P):
    """Re* = rho <v> d_e / (mu (1 - eps)) with rho <v> the mass flux (nomenclature, p. 61)."""
    re = flux * op["d_e"] / (1.0 - op["eps"])
    return re / mu if P["re_mu"] else re


def galileo(rho, mu, op):
    """Ga* = rho^2 g d_e^3 eps^3 / (mu^2 (1 - eps)^3) (nomenclature, p. 61)."""
    e = op["eps"]
    return rho**2 * op["g"] * op["d_e"]**3 * e**3 / (mu**2 * (1.0 - e)**3)


def eotvos(op, sigma=None):
    """Eo* = rho_L g d_e^2 eps^2 / (sigma (1 - eps)^2) (Eq. 7, p. 54)."""
    s = op["sigma"] if sigma is None else sigma
    e = op["eps"]
    return op["rho_L"] * op["g"] * op["d_e"]**2 * e**2 / (s * (1.0 - e)**2)


def static_holdup(eo, P):
    """eps0 = 1/(20 + 0.9 Eo*) (Eq. 8, p. 54)."""
    return 1.0 / (P["eo_a"] + P["eo_b"] * eo)


def drag(re, ga, P):
    """X = A Re*/Ga* + B Re*^2/Ga*: single-phase Ergun drag over rho g (Eqs. 11, 12)."""
    return P["A"] * re / ga + P["B"] * re**2 / ga


def groups(op, P):
    """Every dimensionless group of one operating condition (or an array of them)."""
    re_l = reynolds(op["L"], op["mu_L"], op, P)
    re_g = reynolds(op["G"], op["mu_G"], op, P)
    ga_l = galileo(op["rho_L"], op["mu_L"], op)
    ga_g = galileo(op["rho_G"], op["mu_G"], op)
    eo = eotvos(op)
    return dict(Re_L=re_l, Re_G=re_g, Ga_L=ga_l, Ga_G=ga_g, Eo=eo,
                eps0=static_holdup(eo, P), X_L=drag(re_l, ga_l, P), X_G=drag(re_g, ga_g, P))


# codehealth: allow params -- physical inputs and the solver settings the break table turns
def ergun_gradient(rho, mu, flux, op, A, B):
    """Dimensional Ergun pressure gradient (Pa/m) for one phase filling the bed."""
    e, d = op["eps"], op["d_e"]
    u = flux / rho
    return A * mu * u * (1 - e)**2 / (d**2 * e**3) + B * rho * u**2 * (1 - e) / (d * e**3)


# ---------------------------------------------------------------- stagnant gas

def stagnant_delta(x_l, P):
    """Reduced saturation with stagnant gas: delta = X_L^(1/a), or X_L^e34 (Eqs. 33, 34)."""
    e = 1.0 / P["a"] if P["stagnant_exp"] is None else P["stagnant_exp"]
    return x_l**e


def holdup(delta, eps, eps0):
    """eps_L from the reduced saturation (Eq. 24 inverted)."""
    return eps0 + delta * (eps - eps0)


def gas_saturation(delta, eps, eps0):
    """S_G = 1 - eps_L/eps (Eq. 52)."""
    return (eps - eps0) * (1.0 - delta) / eps


def trickle_end_re(ga_l, P):
    """Re_L where X_L = 1, so delta = 1 with stagnant gas: the root of B Re^2 + A Re - Ga = 0."""
    A, B = P["A"], P["B"]
    if B == 0:
        return ga_l / A
    return (-A + np.sqrt(A * A + 4.0 * B * ga_l)) / (2.0 * B)


# ---------------------------------------------------------------- the coupled pair

def _pieces(op, P):
    g = groups(op, P)
    k = op["rho_G"] / op["rho_L"]
    one = 1.0 - k if P["exact"] else 1.0
    return g, k, one


def pair_residual(x, op, P, drop_coupling_jac=False):
    """Residual and analytic sparse Jacobian of the pair (delta, psi) per point.

    r1 = X_L delta^-a - (rho_G/rho_L) psi - 1   (Eq. 51; 1 - rho_G/rho_L with ``exact``)
    r2 = psi S_G^c - X_G                         (Eq. 50 with the pole multiplied out)

    ``x`` has shape (n, 2): column 0 delta, column 1 psi. The Jacobian is block
    diagonal with one 2x2 block per point, placed by ``BlockLayout``. ``drop_coupling_jac``
    leaves out dr1/dpsi (a defect for the break table, not a model option).
    """
    g, k, one = _pieces(op, P)
    a, c, eps, eps0 = P["a"], P["c"], op["eps"], g["eps0"]
    x = np.asarray(x, dtype=float).reshape(-1, 2)
    n = x.shape[0]
    d, psi = x[:, 0], x[:, 1]
    s = gas_saturation(d, eps, eps0) * np.ones(n)
    ds = -(eps - eps0) / eps * np.ones(n)
    r1 = g["X_L"] * d**(-a) - k * psi - one
    r2 = psi * s**c - g["X_G"]
    j11 = -a * g["X_L"] * d**(-a - 1) * np.ones(n)
    j12 = (0.0 if drop_coupling_jac else -k) * np.ones(n)
    j21 = psi * c * s**(c - 1) * ds
    j22 = s**c
    lay = BlockLayout.from_fields(["delta", "psi"], (n,))   # the (n, 2) state, fields interleaved
    dg = lambda v: sp.diags_array(v, format="csc")  # noqa: E731
    jac = (lay.place(dg(j11), "delta", "delta") + lay.place(dg(j12), "delta", "psi")
           + lay.place(dg(j21), "psi", "delta") + lay.place(dg(j22), "psi", "psi"))
    return np.column_stack([r1, r2]).ravel(), jac


# codehealth: allow params -- physical inputs and the solver settings the break table turns
def solve_pair(op, P, maxfev=100, require_success=True, drop_coupling_jac=False):
    """Solve the pair for every operating point at once with ``pymrm.newton``.

    Started from the stagnant-gas closed form; ``bounds`` keep delta inside
    (0, 1) where the residual is defined; the unknowns are scaled (delta ~ 1,
    psi ~ 300) and ``rtol`` 1e-12 with ``tol`` 0 makes the step test relative
    (pitfall P7). Returns (delta, psi, result).
    """
    g = groups(op, P)
    n = np.broadcast(*[np.asarray(v) for v in op.values()]).size
    d0 = np.clip(stagnant_delta(g["X_L"] * np.ones(n), P), 1e-12, 1 - 1e-12)
    psi0 = g["X_G"] / gas_saturation(d0, op["eps"], g["eps0"])**P["c"] * np.ones(n)
    x0 = np.column_stack([d0, np.maximum(psi0, 1e-300)])
    lo = np.column_stack([np.zeros(n), np.zeros(n)])
    hi = np.column_stack([np.ones(n), np.full(n, np.inf)])
    res = newton(lambda x: pair_residual(x, op, P, drop_coupling_jac), x0, tol=0.0, rtol=1e-12,
                 bounds=(lo, hi), scale=np.array([1.0, max(float(np.max(psi0)), 1.0)]),
                 maxfev=maxfev)
    if require_success:
        _require_converged(res)
    return res.x[:, 0], res.x[:, 1], res


def _require_converged(res):
    """Raise unless newton converged with every delta inside (0, 1)."""
    if not res.success:
        raise RuntimeError(f"newton failed: {res.message}")
    if np.any(res.x[:, 0] <= 0) or np.any(res.x[:, 0] >= 1):
        raise RuntimeError("newton returned delta outside (0, 1)")


def eliminated(delta, g, op, P):
    """The scalar equation f(delta) = 0 with psi eliminated (Eq. 51 with Eq. 50)."""
    k = op["rho_G"] / op["rho_L"]
    one = 1.0 - k if P["exact"] else 1.0
    s = gas_saturation(delta, op["eps"], g["eps0"])
    return g["X_L"] * delta**(-P["a"]) - one - k * g["X_G"] / s**P["c"]


def solve_scalar(op, P):
    """Second route: brentq on the eliminated equation, one point at a time.

    f falls from +inf at delta -> 0 (L > 0) to -inf at delta -> 1 (G > 0), and is
    monotone between, so the bracket (0, 1) holds exactly one root.
    """
    keys = list(op)
    arrs = np.broadcast_arrays(*[np.asarray(op[kk], dtype=float) for kk in keys])
    flat = [a.ravel() for a in arrs]
    out = np.empty(flat[0].size)
    for i in range(out.size):
        o = {kk: flat[j][i] for j, kk in enumerate(keys)}
        g = groups(o, P)
        out[i] = brentq(eliminated, 1e-15, 1 - 1e-15, args=(g, o, P), xtol=1e-16, rtol=1e-15)
    return out


def state(op, P, delta=None, psi=None):
    """Every derived quantity at the solved state (or solved here when not given)."""
    g = groups(op, P)
    if delta is None:
        delta, psi, _ = solve_pair(op, P)
    eps = op["eps"]
    eps_l = holdup(delta, eps, g["eps0"])
    s = gas_saturation(delta, eps, g["eps0"])
    rg = op["rho_G"] * op["g"]
    return dict(g, delta=delta, psi=psi, eps_L=eps_l, S_G=s,
                gas_mod_grad=rg * psi, minus_dpdz=rg * (psi - 1.0),
                coupling=op["rho_G"] / op["rho_L"] * psi,
                dlnpsi_dlneps=P["c"] * eps_l / (eps - eps_l))


def psi_factor(eps_l, eps, factor, P):
    """Relative change of psi at fixed flows (Eq. 50) when eps_L is multiplied by ``factor``."""
    return ((eps - factor * eps_l) / (eps - eps_l))**(-P["c"]) - 1.0


def momentum_residuals(st, op, P):
    """Eqs. 13 and 14 over rho g of their phase, at a state with -dP/dz from the gas balance.

    Returns (liquid, gas). With Eq. 51 the gas balance closes exactly and the
    liquid one is off by -rho_G/rho_L (the gas weight Eq. 51 drops).
    """
    a, c = P["a"], P["c"]
    mdp = st["minus_dpdz"]
    liq = (mdp + op["rho_L"] * op["g"]) / (op["rho_L"] * op["g"]) - st["X_L"] / st["delta"]**a
    gas = (mdp + op["rho_G"] * op["g"]) / (op["rho_G"] * op["g"]) - st["X_G"] / st["S_G"]**c
    return liq, gas


# ---------------------------------------------------------------- Table 5

def ratio_direct(eps, eo_o, s, P):
    """eps_D,i / eps_D,o = (eps - eps0_i)/(eps - eps0_o) at equal delta; s = sigma_o/sigma_i."""
    return (eps - static_holdup(eo_o * s, P)) / (eps - static_holdup(eo_o, P))


def ratio_rederived(eps, eo_o, s, P):
    """The same ratio in the paper's zeta notation, with Eo*_o kept (Eqs. 35 to 37 re-derived)."""
    za = P["eo_b"] + P["eo_a"] / eo_o
    zb = P["eo_b"] * s + P["eo_a"] / eo_o
    return za / zb * (zb * eo_o * eps - 1.0) / (za * eo_o * eps - 1.0)


def ratio_printed(eps, eo_o, s, P):
    """Eq. 35 as printed: (zeta_A/zeta_B)(zeta_B eps - 1)/(zeta_A eps - 1), Eo*_o missing."""
    za = P["eo_b"] + P["eo_a"] / eo_o
    zb = P["eo_b"] * s + P["eo_a"] / eo_o
    return za / zb * (zb * eps - 1.0) / (za * eps - 1.0)


def ratio_slope_numerator(eo_o, eps, s, P):
    """Numerator of d(ratio_direct)/dEo*_o (pole multiplied out): s u^2 (eps - v) - v^2 (eps - u)."""
    u = static_holdup(eo_o * s, P)
    v = static_holdup(eo_o, P)
    return s * u * u * (eps - v) - v * v * (eps - u)


# codehealth: allow params -- physical inputs and the solver settings the break table turns
def ratio_max(eps, s, P, eo_lo=1e-3, eo_hi=1e4):
    """Largest ratio over Eo*_o at porosity ``eps``: the root of the analytic slope.

    Returns (ratio, Eo*_o). Raises unless ``find_roots`` finds exactly one root.
    """
    r = find_roots(lambda e: ratio_slope_numerator(e, eps, s, P), eo_lo, eo_hi, n_scan=200, log=True,
                   xtol=1e-14)
    if len(r["roots"]) != 1:
        raise RuntimeError(f"ratio maximum: {len(r['roots'])} roots ({r['message']})")
    e = r["roots"][0]
    return float(ratio_direct(eps, e, s, P)), float(e)


def ratio_sup(eps, P):
    """Bound over every Eo*: eps0 lies in (0, 1/eo_a], so the ratio is below 1 + eps0max/(eps - eps0max)."""
    m = 1.0 / P["eo_a"]
    return 1.0 + m / (eps - m)


def eps_for_ratio(target, P):
    """Porosity below which the bound ratio_sup reaches ``target``."""
    m = 1.0 / P["eo_a"]
    return m + m / (target - 1.0)


# codehealth: allow params -- physical inputs and the solver settings the break table turns
def printed_eo_at(eps, target, s, P, eo_lo=1e-2, eo_hi=1e4):
    """Eo*_o that makes the printed Eq. 35 equal ``target`` at ``eps`` (nan unless exactly one)."""
    f = lambda e: ratio_printed(eps, e, s, P) - target  # noqa: E731
    # the printed form has a pole where zeta_A eps = 1: scan each side of it
    spans = [(eo_lo, eo_hi)]
    if 1.0 / eps > P["eo_b"]:
        pole = P["eo_a"] / (1.0 / eps - P["eo_b"])
        spans = [(lo, hi) for lo, hi in ((eo_lo, min(eo_hi, pole * (1 - 1e-9))),
                                         (max(eo_lo, pole * (1 + 1e-9)), eo_hi)) if hi > lo]
    roots = []
    for lo, hi in spans:
        r = find_roots(f, lo, hi, n_scan=120, log=True, xtol=1e-14)
        roots += [e for e in r["roots"] if abs(f(e)) < 1e-9]
    return roots[0] if len(roots) == 1 else np.nan


# codehealth: allow params -- physical inputs and the solver settings the break table turns
def printed_inversion(t1, t2, s1, s2, P, eps_lo=0.27, eps_hi=0.95):
    """Every (eps, Eo*_o) at which the printed Eq. 35 returns both Calc values.

    Nested: Eo*_o from the first target at each eps, then the roots in eps of the
    second, all by ``find_roots``. Returns a list of (eps, Eo*_o).
    """
    def g(eps):
        e = printed_eo_at(eps, t1, s1, P)
        return ratio_printed(eps, e, s2, P) - t2 if np.isfinite(e) else np.nan

    r = find_roots(g, eps_lo, eps_hi, n_scan=137, xtol=1e-14)
    return [(x, printed_eo_at(x, t1, s1, P)) for x in r["roots"] if abs(g(x)) < 1e-9]


# ---------------------------------------------------------------- printed arithmetic

def half_unit(printed: str) -> float:
    """Half a unit in the last printed digit of a numeric string."""
    p = printed.strip().replace(",", "")
    dec = len(p.split(".")[1]) if "." in p else 0
    return 0.5 * 10.0**(-dec)


def recombined(ns, means_printed, total_printed):
    """Count-weighted mean of printed subset means against the printed total.

    Returns (mean, |mean - total| over the rounding allowance). The allowance is
    the shift the rounding of the inputs can cause plus half a unit of the total.
    """
    ns = np.asarray(ns, dtype=float)
    vals = np.array([float(m) for m in means_printed])
    w = ns / ns.sum()
    mean = float(w @ vals)
    allowance = float(w @ np.array([half_unit(m) for m in means_printed])) + half_unit(total_printed)
    return mean, abs(mean - float(total_printed)) / allowance
