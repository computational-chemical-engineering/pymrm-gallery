"""Yates and Satterfield (1991): eq. 16 on cobalt, its apparent activation energy, its
parameter covariance from the experimental design, and the slurry CSTR and an ideal PFR.

Units follow the paper: partial pressures in MPa, the rate -R_(H2+CO) (H2 plus CO consumed)
in mmol/(min g of unreduced catalyst), a in mmol/(min g MPa^2), b in 1/MPa, feed in
L(STP)/(min g). Energies are returned in kJ/mol.

Every published constant sits in one dictionary ``par`` (keys ``a240, b240, a220, b220``
and the temperatures ``T1`` = 220 C, ``T2`` = 240 C in kelvin), which every function takes,
so a break row per constant is one call. No plotting and no data loading here.
"""
from __future__ import annotations

import numpy as np
from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)
from pymrm.checks import find_roots

R_GAS = 8.314462618          # J/(mol K)


def rate(a, b, p_h2, p_co):
    """Eq. 16: -R_(H2+CO) = a P_H2 P_CO / (1 + b P_CO)^2 (mmol/(min g))."""
    return a * p_h2 * p_co / (1.0 + b * p_co) ** 2


def rate_dpco(a, b, p_h2, p_co):
    """Analytic d(-R)/dP_CO at fixed P_H2: a P_H2 (1 - b P_CO)/(1 + b P_CO)^3."""
    return a * p_h2 * (1.0 - b * p_co) / (1.0 + b * p_co) ** 3


def inv_dT(par):
    """D = 1/T1 - 1/T2 (1/K), T1 the 220 C and T2 the 240 C temperature."""
    return 1.0 / par["T1"] - 1.0 / par["T2"]


def e_app(p_co, par, log=np.log):
    """Apparent activation energy at fixed P (kJ/mol); P_H2 cancels.

    E = (R/D) [ln(a240/a220) + 2 ln((1 + b220 P)/(1 + b240 P))]. ``log`` is the
    logarithm used (``np.log``; a break row passes ``np.log10``).
    """
    p = np.asarray(p_co, dtype=float)
    s = log(par["a240"] / par["a220"]) + 2.0 * log((1.0 + par["b220"] * p)
                                                     / (1.0 + par["b240"] * p))
    return R_GAS * s / inv_dT(par) / 1e3


def arrhenius(par, log=np.log):
    """Two-point temperature dependence of a and b.

    Returns ``(E_a, dH_b, T_a, T_b)``: E_a of a and the adsorption enthalpy of b
    (kJ/mol), and the activation temperatures ln(a240/a220)/D and ln(b220/b240)/D (K).
    """
    d = inv_dT(par)
    ta = log(par["a240"] / par["a220"]) / d
    tb = log(par["b220"] / par["b240"]) / d
    return R_GAS * ta / 1e3, -R_GAS * tb / 1e3, ta, tb


def e_app_limits(par, log=np.log):
    """Closed-form limits of E_app: P -> 0 gives E_a, P -> inf gives E_a - 2 dH_b."""
    ea, dhb, _, _ = arrhenius(par, log)
    return ea, ea - 2.0 * dhb


def e_app_root(e_q, par):
    """Closed-form P_CO (MPa) where E_app = ``e_q`` (kJ/mol)."""
    q = np.exp((e_q * 1e3 * inv_dT(par) / R_GAS - np.log(par["a240"] / par["a220"])) / 2.0)
    return (q - 1.0) / (par["b220"] - q * par["b240"])


# codehealth: allow params -- root-finder settings, not physics
def one_root(f, lo, hi, what, n_scan=64, log=False):
    """The single root of ``f`` on [lo, hi] by find_roots; raises unless exactly one."""
    res = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=1e-14)
    if len(res["roots"]) != 1 or res["failed"]:
        raise ValueError(f"{what}: {len(res['roots'])} roots on [{lo}, {hi}]; {res['message']}")
    return res["roots"][0]


# codehealth: allow params -- eq. 16 parameters and the design arrays
def design_jacobian(a, b, p_h2, p_co, b_col=1.0):
    """Analytic sensitivities of eq. 16 to (a, b) at each run, shape (n, 2).

    ``b_col`` multiplies the b column (1 is correct; a break row passes 0.5).
    """
    p_h2, p_co = np.asarray(p_h2, float), np.asarray(p_co, float)
    return np.column_stack([p_h2 * p_co / (1.0 + b * p_co) ** 2,
                            -2.0 * b_col * a * p_h2 * p_co ** 2 / (1.0 + b * p_co) ** 3])


# codehealth: allow params -- eq. 16 parameters, the design arrays and the step
def fd_jacobian(a, b, p_h2, p_co, rel=1e-6):
    """Central-difference sensitivities of eq. 16 to (a, b), the check of the analytic ones."""
    ha, hb = rel * a, rel * b
    return np.column_stack([(rate(a + ha, b, p_h2, p_co) - rate(a - ha, b, p_h2, p_co)) / (2 * ha),
                            (rate(a, b + hb, p_h2, p_co) - rate(a, b - hb, p_h2, p_co)) / (2 * hb)])


# codehealth: allow params -- eq. 16 parameters, the design arrays and the scale
def design_covariance(a, b, p_h2, p_co, sd_b=None, b_col=1.0):
    """Covariance of (a, b) from the design: C = s^2 (J^T J)^-1 at the estimate.

    With ``sd_b`` given, s^2 is chosen so that sqrt(C[1, 1]) = sd_b; otherwise s^2 = 1.
    The sd ratio sqrt(C00/C11) and the correlation do not depend on s^2.
    """
    j = design_jacobian(a, b, p_h2, p_co, b_col)
    c = np.linalg.inv(j.T @ j)
    if sd_b is not None:
        c = c * sd_b ** 2 / c[1, 1]
    return c


def sd_ratio_rho(c):
    """sd_a/sd_b and the correlation of a and b from a covariance matrix."""
    return np.sqrt(c[0, 0] / c[1, 1]), c[0, 1] / np.sqrt(c[0, 0] * c[1, 1])


# codehealth: allow params -- pressure, parameters and one covariance per temperature
def sd_e_app(p_co, par, c240, c220, rho_zero=False):
    """Standard error of E_app(P) (kJ/mol) from the two covariances, temperatures independent.

    ``rho_zero`` drops the a-b covariance (the variances are kept).
    """
    p = np.asarray(p_co, dtype=float)
    if rho_zero:
        c240, c220 = np.diag(np.diag(c240)), np.diag(np.diag(c220))
    g240 = np.stack([np.full_like(p, 1.0 / par["a240"]), -2.0 * p / (1.0 + par["b240"] * p)])
    g220 = np.stack([np.full_like(p, -1.0 / par["a220"]), 2.0 * p / (1.0 + par["b220"] * p)])
    var = (np.einsum("i...,ij,j...->...", g240, c240, g240)
           + np.einsum("i...,ij,j...->...", g220, c220, g220))
    return R_GAS * np.sqrt(var) / inv_dT(par) / 1e3


# ---------------------------------------------------------------- reactors
def feed_flow(q_l, t_std=273.15, p_std=101325.0):
    """Syngas feed (mmol/(min g)) from L(STP)/(min g) at the stated standard state."""
    return q_l / (R_GAS * t_std / p_std * 1e3) * 1e3


# codehealth: allow params -- the Table II row: rate, feed, ratio, pressure, usage ratio, standard state
def closure(r_meas, q_l, m, p_tot, u_ratio=2.0, t_std=273.15):
    """Outlet P_H2, P_CO (MPa) and X_CO implied by a measured rate in a CSTR.

    Gas H2, CO and H2O only (hydrocarbons not in the gas), CO used = rate/(1 + U).
    """
    f0 = feed_flow(q_l, t_std)
    f_co0, f_h20 = f0 / (1.0 + m), m * f0 / (1.0 + m)
    used = r_meas / (1.0 + u_ratio)
    f = np.array([f_h20 - u_ratio * used, f_co0 - used, used])
    return p_tot * f[0] / f.sum(), p_tot * f[1] / f.sum(), used / f_co0


# codehealth: allow params -- state, rate-law parameters and the reading of the rate
def _flows_rate(y, a, b, p_tot, rate_is_co, u_ratio):
    """Partial pressures from scaled flows (CO, H2, H2O) and the CO consumption rate."""
    tot = y[..., 0] + y[..., 1] + y[..., 2]
    p_co, p_h2 = p_tot * y[..., 0] / tot, p_tot * y[..., 1] / tot
    r = rate(a, b, p_h2, p_co)
    return r if rate_is_co else r / (1.0 + u_ratio)


# codehealth: allow params -- rate-law parameters, operating point and solver settings
def cstr(a, b, p_tot, tau, m=2.0, u_ratio=2.0, rate_is_co=False, tol=1e-13):
    """Isothermal isobaric CSTR by pymrm ``newton`` on the three species balances.

    Unknowns y = (F_CO, F_H2, F_H2O)/F_CO0; tau = W/F_CO0 (g min/mmol). Balances
    y_in - y - nu tau r_CO = 0 with nu = (1, U, -1) consumed. Returns a dict with X,
    the H2 + CO rate, P_CO, P_H2, the newton result and its final residual.
    """
    y_in = np.array([[1.0, m, 0.0]])
    nu = np.array([[1.0, u_ratio, -1.0]])
    jac = NumJac((1, 3))

    def residual(y):
        def g(v):
            return y_in - v - nu * tau * _flows_rate(v, a, b, p_tot, rate_is_co, u_ratio)[..., None]
        return jac(g, y)

    res = newton(residual, y_in.copy(), tol=tol, bounds=(0.0, None), stagnation=True, maxfev=60)
    if not res.success:
        raise RuntimeError(f"cstr: {res.message}")
    y = res.x.reshape(1, 3)
    tot = y.sum()
    p_co, p_h2 = p_tot * y[0, 0] / tot, p_tot * y[0, 1] / tot
    return {"X": 1.0 - y[0, 0], "rate": rate(a, b, p_h2, p_co), "P_CO": p_co, "P_H2": p_h2,
            "residual": float(np.max(np.abs(residual(y)[0]))), "result": res}


def cstr_wf(x, a, b, p_tot):
    """Closed-form W/F_CO0 of the CSTR at CO conversion x (H2/CO fed = U = 2)."""
    u = 1.0 - x
    c = 2.0 + b * p_tot
    return 3.0 * x * (1.0 + c * u) ** 2 / (2.0 * a * p_tot ** 2 * u ** 2)


def pfr_wf(x, a, b, p_tot):
    """Closed-form W/F_CO0 of the ideal isothermal PFR at CO conversion x (H2/CO = U = 2)."""
    u = 1.0 - x
    c = 2.0 + b * p_tot
    return 3.0 / (2.0 * a * p_tot ** 2) * (1.0 / u - 1.0 - 2.0 * c * np.log(u) + c ** 2 * (1.0 - u))


def outlet_pco(x, p_tot):
    """Outlet P_CO (MPa) at CO conversion x for H2/CO fed = U = 2: P u/(1 + 2u)."""
    u = 1.0 - x
    return p_tot * u / (1.0 + 2.0 * u)


class PlugFlow:
    """Ideal isothermal isobaric PFR of the slurry kinetics, solved with pymrm.

    Unknowns y = (F_CO, F_H2, F_H2O)/F_CO0 on the scaled catalyst mass xi = W/W_tot in
    [0, 1], shape ``(n, 3)``. Steady convection of the molar flows at unit velocity
    (``construct_convflux_upwind``, first order; ``construct_div`` with ``nu=0``,
    Cartesian) with the source -tau nu r_CO, tau = W_tot/F_CO0. The varying total flow
    needs no velocity: the flows themselves are convected (checklist rule 22).

    Parameters
    ----------
    a, b, p_tot : float
        Eq. 16 parameters and total pressure (MPa).
    n : int
        Number of cells.
    m, u_ratio : float
        H2/CO fed and the usage ratio U.
    """

    # codehealth: allow params -- rate-law parameters, operating point and grid
    def __init__(self, a, b, p_tot, n=200, m=2.0, u_ratio=2.0):
        self.a, self.b, self.p_tot, self.m, self.u_ratio = a, b, p_tot, m, u_ratio
        self.shape = (int(n), 3)
        self.x_f = np.linspace(0.0, 1.0, self.shape[0] + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.nu = np.array([1.0, u_ratio, -1.0])
        # xi = 0, outward normal -xi: y = y_in per field (Dirichlet) -> a = 0, b = 1, d = y_in
        # xi = 1, outward normal +xi: convective outflow, the face carries the cell value
        self.bc = ({"a": 0.0, "b": 1.0, "d": np.array([[1.0, m, 0.0]])}, {"outflow": True})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c, self.bc,
                                                  v=1.0, axis=0)
        div = construct_div(self.shape, self.x_f, nu=0, axis=0)   # nu=0: Cartesian
        self.jac_const = div @ conv
        self.g_const = (div @ conv_bc).toarray().reshape(-1, 1)
        self.numjac = NumJac(self.shape)          # pointwise source, the 3 flows coupled

    def residual(self, y, tau):
        def src(v):
            return -tau * self.nu * _flows_rate(v, self.a, self.b, self.p_tot, False,
                                                self.u_ratio)[..., None]
        g_src, jac_src = self.numjac(src, y)
        g = self.g_const + self.jac_const @ y.reshape(-1, 1) - g_src.reshape(-1, 1)
        return g, self.jac_const - jac_src

    def solve(self, tau, permc_spec=None):
        """Outlet CO conversion at tau = W/F_CO0 (g min/mmol), read from the outlet face."""
        y0 = np.tile([1.0, self.m, 0.0], (self.shape[0], 1))
        kw = {} if permc_spec is None else {"lin_solver_kwargs": {"permc_spec": permc_spec}}
        res = newton(lambda y: self.residual(y, tau), y0, tol=1e-12, bounds=(0.0, None),
                     stagnation=True, maxfev=60, **kw)
        if not res.success:
            raise RuntimeError(f"PlugFlow.solve: {res.message}")
        self.y = res.x.reshape(self.shape)
        _, _, y_out, _ = compute_boundary_values(self.y, self.x_f, self.x_c, self.bc)
        return 1.0 - float(np.asarray(y_out).reshape(-1)[0])
