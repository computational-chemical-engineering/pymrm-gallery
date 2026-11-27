"""Membrane water balance of Springer, Zawodzinski and Gottesfeld (1991), page J3.6.

The property correlations (Eqs. 15 to 18, 21, 22, 25), the gas-side relations (Eqs. 7,
8, 27) and the membrane water balance (Eqs. 23, 24, 26) as finite volumes with pymrm.
Every printed constant arrives through one dictionary ``C`` (``constants``), so a
mistyped constant is one changed entry. Units are the paper's: cm, mol, s, atm,
A/cm2, Ohm cm2, temperatures in degC unless the name says K.
"""

import math
from dataclasses import dataclass

import numpy as np
from numpy.polynomial import Polynomial
from scipy.sparse import block_array, csc_array

from pymrm import (NumJac, compute_boundary_values, construct_div, construct_grad,
                   interp_cntr_to_stagg, newton)
from pymrm.checks import find_roots

C_REF = 1e-3          # mol/cm3, scale of the fixed-charge concentration rho_dry/M_m
D_UNIT = 1e-6         # cm2/s, the unit of the Eq. 22 polynomial


def constants(printed):
    """Floats from the printed strings ``{item: str}`` of the printed-numbers CSV."""
    return {k: float(v) for k, v in printed.items()}


def psat(T, C):
    """Eq. 15: saturation pressure of water, atm, at T in degC."""
    lg = C["eq15_c0"] + C["eq15_c1"] * T + C["eq15_c2"] * T**2 + C["eq15_c3"] * T**3
    return 10.0**lg


def lam_eq16(a, C):
    """Eq. 16: membrane water content at 30 C for activity 0 < a <= 1."""
    return C["eq16_c0"] + C["eq16_c1"] * a + C["eq16_c2"] * a**2 + C["eq16_c3"] * a**3


def lam_iso(a, C):
    """Water content at a membrane face: Eq. 16 for a <= 1, Eq. 17 for 1 < a <= 3."""
    if 0.0 < a <= 1.0:
        return lam_eq16(a, C)
    if 1.0 < a <= 3.0:
        return C["eq17_c0"] + C["eq17_c1"] * (a - 1.0)
    raise ValueError(f"activity {a:g} outside (0, 3], where Eqs. 16 and 17 hold")


def ndrag(lam, C):
    """Eq. 18: water molecules dragged per proton, 2.5 lambda/22."""
    return C["eq18_ndrag_sat"] * lam / C["eq18_lambda_ref"]


def drag(lam, C, reading="corrected"):
    """Drag water flux per H2 flux I at water content lam.

    ``corrected``: N_w,drag = n_drag (2I) (Eq. 18), so 2 n_drag. ``literal``: Eqs. 23
    and 24 as printed, n_drag (2I) lambda/22, which counts the lambda/22 twice.
    """
    d = 2.0 * ndrag(lam, C)
    if reading == "literal":
        d = d * lam / C["eq18_lambda_ref"]
    return d


def arrhenius(E, T, C):
    """exp[E (1/303 - 1/(273 + T))], the temperature factor of Eqs. 22 and 25."""
    return math.exp(E * (1.0 / C["T_ref_K"] - 1.0 / (273.0 + T)))


def d_poly(C):
    """The Eq. 22 polynomial in lambda (units of 1e-6 cm2/s)."""
    return Polynomial([C["eq22_c0"], C["eq22_c1"], C["eq22_c2"], C["eq22_c3"]])


def d_lambda(lam, T, C, arrh=True):
    """Eq. 22: diffusion coefficient D_lambda, cm2/s, for lambda > 4."""
    f = arrhenius(C["eq22_E"], T, C) if arrh else 1.0
    return D_UNIT * f * d_poly(C)(lam)


def d_kirchhoff(l2, l3, T, C, arrh=True):
    """Integral of D_lambda over lambda from l2 to l3, cm2/s, in closed form."""
    f = arrhenius(C["eq22_E"], T, C) if arrh else 1.0
    P = d_poly(C).integ()
    return D_UNIT * f * (P(l3) - P(l2))


def sigma(lam, T, C):
    """Eq. 25: membrane conductivity, (Ohm cm)^-1, for lambda > 1."""
    return arrhenius(C["eq25_E"], T, C) * (C["eq25_c1"] * lam + C["eq25_c0"])


def d_ratio(lam, C, eq20=False):
    """Eq. 21 correction factor D_lambda/D' at 30 C; ``eq20`` keeps Eq. 20's extra lambda.

    The activity is the single root of Eq. 16 on (0, 1] (``find_roots``).
    """
    r = find_roots(lambda a: lam_eq16(a, C) - lam, 1e-9, 1.0, n_scan=64, xtol=1e-15)
    if len(r["roots"]) != 1:
        raise ValueError(f"Eq. 16 at lambda {lam:g}: {r['message']}")
    a = r["roots"][0]
    dlda = C["eq16_c1"] + 2 * C["eq16_c2"] * a + 3 * C["eq16_c3"] * a**2
    s = C["eq21_s"]
    ratio = lam / ((1.0 + s * lam) ** 2 * a * dlda)
    return ratio * lam if eq20 else ratio


def alpha_eq8(xw1, xI, nuH):
    """Eq. 8 solved for alpha (it is linear-fractional in alpha)."""
    return (nuH * xI - xw1 * (xI + nuH - 1.0)) / ((1.0 - xI) * (1.0 - xw1))


def alpha_eq7_xw4(xw4, xI, nuO, xON):
    """Eq. 7 (x_w4) solved for alpha."""
    q = (1.0 - xI) * xON
    return (xw4 * (nuO + q) - xI * nuO - 2.0 * q) / (2.0 * q * (1.0 - xw4))


def alpha_eq7_xO4(xO4, xI, nuO, xON):
    """Eq. 7 (x_O4) solved for alpha."""
    q = (1.0 - xI) * xON
    return ((nuO - 1.0) * q / xO4 - nuO - q) / (2.0 * q)


def vcell_eq27(J, R, xO3, xliq, C, T_K, Rgas):
    """Eq. 27 solved for V_cell, with eta = V_oc - V_cell - J R_mem."""
    F = C["F"]
    eta = Rgas * T_K / (C["eq27_tafel_factor"] * F) * math.log(
        J * (1.0 - xliq) / (C["j_o"] * C["P_C"] * xO3))
    return C["V_oc"] - J * R - eta


def closed_const_d(l2, l3, kd, pe):
    """Constant D and linear drag (V1): alpha and lambda(zeta) in closed form.

    alpha = kd lambda - (c_f D/(I t)) dlambda/dzeta with lambda(0) = l2, lambda(1) = l3
    gives lambda = alpha/kd + K exp(Pe zeta), K = (l3 - l2)/(exp(Pe) - 1),
    alpha = kd (l2 - K), with Pe = kd I t/(c_f D).
    """
    K = (l3 - l2) / math.expm1(pe)
    alpha = kd * (l2 - K)
    return alpha, (lambda z: alpha / kd + K * np.exp(pe * np.asarray(z)))


@dataclass(frozen=True)
class Scheme:
    """Discretisation choices (break rows change them)."""
    face: str = "mean"          # "mean": arithmetic mean of the two cells; "left": one-sided
    tol: float = 1e-11


class Membrane:
    """Steady membrane water balance d/dz [I drag(lambda) - c_f D(lambda) dlambda/dz] = 0.

    Parameters
    ----------
    n : int
        Number of cells on zeta = z/t_mem in [0, 1] (z from the anode face).
    lam2, lam3 : float
        Water content at the anode face (interface 2) and the cathode face (3).
    J, T, t : float
        Current density A/cm2, cell temperature degC, membrane thickness cm.
    C : dict
        The printed constants (``constants``).
    reading : str
        ``corrected`` or ``literal`` drag (see ``drag``).
    arrh : bool
        Include the Arrhenius factor of Eq. 22.
    d_const : float or None
        Constant D in cm2/s in place of Eq. 22 (V1).
    scheme : Scheme
        Discretisation choices.
    """

    def __init__(self, n, lam2, lam3, J, T, t, C, reading="corrected", arrh=True,
                 d_const=None, scheme=Scheme()):
        self.n, self.lam2, self.lam3, self.T, self.t, self.C = n, lam2, lam3, T, t, C
        self.reading, self.arrh, self.d_const, self.scheme = reading, arrh, d_const, scheme
        self.I = J / (2.0 * C["F"])                    # Eq. 1, mol H2/(cm2 s)
        self.q = C_REF * D_UNIT / t                    # flux unit, mol/(cm2 s)
        self.x_f = np.linspace(0.0, 1.0, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # zeta = 0 (anode face): lambda = lambda_2 (Eq. 16 at x_w2), a Dirichlet value;
        # zeta = 1 (cathode face): lambda = lambda_3 (Eq. 17 at x_w3), a Dirichlet value.
        # Dirichlet rows have a = 0, so the outward-normal sign does not enter.
        self.bc = ({"a": 0.0, "b": 1.0, "d": lam2}, {"a": 0.0, "b": 1.0, "d": lam3})
        shape = (n, 1)
        self.grad, gbc = construct_grad(shape, self.x_f, self.x_c, self.bc, axis=0)
        self.gbc = np.asarray(gbc.todense() if hasattr(gbc, "todense") else gbc).ravel()
        self.div = construct_div(shape, self.x_f, nu=0, axis=0)   # nu = 0: Cartesian slab
        self.numjac = NumJac(shape, axes_diagonals=[0])          # fluxes couple neighbours
        self.lam = np.linspace(lam2, lam3, n).reshape(shape)
        self.u = 1.8                                          # c_f / C_REF

    def d_hat(self, lam):
        """D in units of D_UNIT."""
        if self.d_const is not None:
            return np.full_like(np.asarray(lam, dtype=float), self.d_const / D_UNIT)
        return d_lambda(lam, self.T, self.C, self.arrh) / D_UNIT

    def face_lambda(self, lam):
        """lambda on the faces: the Dirichlet values at both ends, inside the arithmetic
        mean of the two cells. D is a smooth cubic in lambda (no jump), so D at the mean
        face lambda is second order like the harmonic mean of D (P2 applies at jumps)."""
        lam = lam.reshape(-1)
        if self.scheme.face == "left":
            lf = np.concatenate(([lam[0]], lam))
        else:
            lf = interp_cntr_to_stagg(lam.reshape(-1, 1), self.x_f, self.x_c).reshape(-1)
        lf[0], lf[-1] = self.lam2, self.lam3
        return lf

    def face_flux(self, lam, u):
        """Water flux on every face in units of q: drag minus back-diffusion."""
        lf = self.face_lambda(lam)
        g = self.grad @ lam.reshape(-1) + self.gbc
        return (self.I / self.q) * drag(lf, self.C, self.reading) - u * self.d_hat(lf) * g

    def _res(self, lam, u):
        return (self.div @ self.face_flux(lam, u)).reshape(-1, 1)

    def _fun(self, lam):
        lam = lam.reshape(self.n, 1)
        return self.numjac(lambda x: self._res(x, self.u), lam)

    def solve(self, u=None):
        """Solve for lambda at c_f = u C_REF (u kept when None)."""
        if u is not None:
            self.u = u
        r = newton(self._fun, self.lam, tol=self.scheme.tol, stagnation=True)
        if not r.success:
            raise RuntimeError(f"membrane solve failed: {r.message}")
        self.lam = r.x.reshape(self.n, 1)
        return self

    def _fun_b(self, x, alpha_target):
        n = self.n
        lam, u = x[:n].reshape(n, 1), x[n]
        g_l, j_l = self.numjac(lambda y: self._res(y, u), lam)
        lf = self.face_lambda(lam)
        g = self.grad @ lam.reshape(-1) + self.gbc
        col = csc_array((self.div @ (-self.d_hat(lf) * g)).reshape(-1, 1))
        k = self.q / self.I
        f0 = self.face_flux(lam, u)[0]
        row = csc_array(-k * u * self.d_hat(lf[0]) * self.grad[[0], :])
        corner = csc_array(np.array([[-k * self.d_hat(lf[0]) * g[0]]]))
        jac = block_array([[j_l, col], [row, corner]], format="csc")
        res = np.concatenate((g_l.reshape(-1), [k * f0 - alpha_target]))
        return res, jac

    def solve_alpha(self, alpha_target, u0=None):
        """Bordered Newton: lambda and c_f together, with the row alpha - target = 0."""
        x0 = np.concatenate((self.lam.reshape(-1), [u0 if u0 is not None else self.u]))
        r = newton(self._fun_b, x0, args=(alpha_target,), tol=self.scheme.tol,
                   stagnation=True)
        if not r.success:
            raise RuntimeError(f"bordered solve failed: {r.message}")
        self.lam, self.u = r.x[: self.n].reshape(self.n, 1), float(r.x[self.n])
        return self

    @property
    def c_f(self):
        return self.u * C_REF

    def boundary_flux(self):
        """Face fluxes (mol/(cm2 s)) at the anode and cathode faces from
        compute_boundary_values under the same bc (gradients along +zeta)."""
        v0, g0, v1, g1 = compute_boundary_values(self.lam, self.x_f, self.x_c, self.bc)
        v0, g0, v1, g1 = (float(np.ravel(v)[0]) for v in (v0, g0, v1, g1))
        f = [self.I * drag(v, self.C, self.reading) - self.q * self.u * self.d_hat(v) * g
             for v, g in ((v0, g0), (v1, g1))]
        return f[0], f[1]

    def alpha(self):
        """Net water flux over I, read at the anode face (compute_boundary_values)."""
        return self.boundary_flux()[0] / self.I

    def resistance(self, swollen=False):
        """Eq. 26, R = t int dzeta/sigma (midpoint rule on the cells); ``swollen`` adds
        the factor (1 + s lambda) of the wet thickness."""
        lam = self.lam.reshape(-1)
        w = (1.0 + self.C["eq21_s"] * lam) if swollen else 1.0
        return float(self.t * np.sum(np.diff(self.x_f) * w / sigma(lam, self.T, self.C)))

    def lam_mid(self):
        """lambda at zeta = 1/2: a face for even n, a cell centre for odd n."""
        if self.n % 2:
            return float(self.lam[self.n // 2, 0])
        return float(self.face_lambda(self.lam)[self.n // 2])
