"""Effectiveness factors for general rate forms (Bischoff 1965) and their references.

Model code of gallery page B1.3. Three routes to the effectiveness factor of a
slab, sharing no assembly, plus the finite-cylinder machinery of the page's
out-of-sample audit and the closed forms of the two papers it reads.

Contents
--------
Power, LHHW, VolChange
    Rate forms in normalised variables, ``r(1) = 1`` and ``d = D(C)/D``, each
    with the closed-form rate integral ``Rint(a, b) = int_a^b d r dc``.
mb_eta, Curve, ZeroOrder, mstar_power
    Route B: quadrature of Bischoff's parametric solution, eqs. (14)-(15).
m_eta_n1, m_eta_n2, m_eta_n3, m_eta_nhalf, m_eq20
    Route C: the paper's closed forms, eqs. (20)-(28).
Slab
    Route A: ``d/dx[d(c) dc/dx] = Phi^2 r(c)`` by pymrm finite volume.
Cyl2D
    A finite cylinder, radial times axial, every face exposed.
mu1_series, mu2_eq30, eta_series_1st, pz_coeffs, eta_pz
    Pan and Zhu (1998): their shape coefficients, the exact first-order
    series, and their polynomial for arbitrary kinetics.
eq33, eq34
    Bischoff's Table 1 correction factors for volume change.

Several functions take an argument whose default is the correct reading and
whose alternative is a defect the page's break table injects (a printed
misprint, a wrong modulus, a wrong geometry factor).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ellipeinc, ellipkinc, expit, jn_zeros, logit

from pymrm import (NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_div, construct_grad, interp_cntr_to_stagg, newton)

__all__ = ["Power", "LHHW", "VolChange", "mb_eta", "Curve", "ZeroOrder", "mstar_power",
           "K15SQ", "K45SQ", "phi_from_cos", "m_eta_n1", "m_eta_n2", "m_eta_n3",
           "m_eta_nhalf", "m_eq20", "SYMMETRY", "SURFACE", "Slab", "Cyl2D",
           "mu1_series", "mu2_eq30", "eta_series_1st", "pz_coeffs", "eta_pz",
           "eq33", "eq34"]


# ---------------------------------------------------------------- rate forms
def _t_minus_log1p(t):
    """``t - ln(1 + t)`` without cancellation for small ``|t|``.

    For ``|t| <= 0.05`` the series ``t^2 sum_j (-t)^j/(j + 2)``, 16 terms
    (truncation below ``1e-19`` relative).
    """
    t = float(t)
    if abs(t) > 0.05:
        return t - math.log1p(t)
    acc = 0.0
    for j in range(15, -1, -1):
        acc = acc * (-t) + 1.0 / (j + 2)
    return t * t * acc


class Power:
    """Power law ``r = c^n`` with constant diffusivity (Bischoff eq. 16).

    Parameters
    ----------
    n : float
        Reaction order.
    """

    variable_d = False

    def __init__(self, n):
        self.n = float(n)
        self.label = f"n = {n:g}"

    def r(self, c):
        """Normalised rate; zero for ``c <= 0``."""
        return np.maximum(c, 0.0) ** self.n

    def dr(self, c):
        """Derivative of :meth:`r` for ``c > 0``."""
        return self.n * np.asarray(c, float) ** (self.n - 1.0)

    def d(self, c):
        """Normalised diffusivity, 1."""
        return np.ones_like(np.asarray(c, float))

    def Rint(self, a, b):
        """``int_a^b r dc``."""
        return (b ** (self.n + 1) - a ** (self.n + 1)) / (self.n + 1)


class LHHW:
    """Adsorption kinetics ``r = c (z + 1)/(z + c)``, ``z = zeta/C_o`` (Bischoff eq. 29).

    Parameters
    ----------
    z : float
        Reduced adsorption constant ``zeta* = zeta/C_o``.
    stable : bool, optional
        ``False`` writes :meth:`Rint` as ``(b - a) - z ln(1 + (b - a)/(z + a))``,
        which cancels as ``b -> a``; it exists for the page's break table.
    """

    variable_d = False

    def __init__(self, z, stable=True):
        self.z, self.stable = float(z), stable
        self.label = f"LHHW z* = {z:g}"

    def r(self, c):
        """Normalised rate; zero for ``c <= 0``."""
        cp = np.maximum(c, 0.0)
        return cp * (self.z + 1.0) / (self.z + cp)

    def d(self, c):
        """Normalised diffusivity, 1."""
        return np.ones_like(np.asarray(c, float))

    def Rint(self, a, b):
        """``int_a^b r dc`` in closed form, without cancellation as ``b -> a``."""
        z = self.z
        t = (b - a) / (z + a)
        if not self.stable:
            return (z + 1.0) * ((b - a) - z * np.log1p(t))
        return (z + 1.0) * (a * t + z * _t_minus_log1p(t))


class VolChange:
    """First order with ``D(C) = D/(1 + w c)``, ``w = omega C_o`` (Bischoff eq. 32).

    Parameters
    ----------
    w : float
        Volume-change group ``omega C_o``.
    stable : bool, optional
        ``False`` writes :meth:`Rint` as ``((b - a) - ln(1 + t)/w)/w``, which
        cancels as ``b -> a``; it exists for the page's break table.
    """

    variable_d = True

    def __init__(self, w, stable=True):
        self.w, self.stable = float(w), stable
        self.label = f"volume change w = {w:+g}"

    def r(self, c):
        """Normalised rate; zero for ``c <= 0``."""
        return np.maximum(c, 0.0)

    def d(self, c):
        """Normalised diffusivity ``1/(1 + w c)``."""
        return 1.0 / (1.0 + self.w * np.asarray(c, float))

    def Rint(self, a, b):
        """``int_a^b d r dc`` in closed form, without cancellation as ``b -> a``."""
        w = self.w
        t = w * (b - a) / (1.0 + w * a)
        if not self.stable:
            return ((b - a) - np.log1p(t) / w) / w
        return (w * a * t + _t_minus_log1p(t)) / w ** 2


# ------------------------------------------------- Route B: eqs. (14)-(15)
U_MIN = 1e-10


def mb_eta(rf, u, denominator_exponent=0.5, epsrel=1e-10):
    """Generalised modulus and effectiveness factor at centre concentration ``u``.

    The outer integral of eq. (14) is taken on ``c = u + s^2``, which removes
    its square-root singularity; the inner rate integral is ``rf.Rint``.

    Parameters
    ----------
    rf : Power, LHHW or VolChange
        Rate form.
    u : float
        Centre concentration ``C_L/C_o`` in ``(0, 1)``.
    denominator_exponent : float, optional
        Exponent of the bracket ``[int_0^Co D r dC]`` in the denominator of
        eq. (14). Eqs. (8) and (13) force ``+1/2``; the paper prints ``-1/2``.
    epsrel : float, optional
        Relative tolerance of the outer quadrature.

    Returns
    -------
    m, eta : float
        Eq. (14) and eq. (15).
    """
    i0 = rf.Rint(0.0, 1.0)
    d = (lambda c: float(rf.d(c))) if rf.variable_d else (lambda c: 1.0)

    def outer(s):
        b = u + s * s
        v = rf.Rint(u, b)
        return d(b) * 2.0 * s / math.sqrt(v) if v > 0 else 0.0

    i_out = quad(outer, 0.0, np.sqrt(1.0 - u), epsabs=1e-12, epsrel=epsrel, limit=200)[0]
    m = 0.5 * i_out / i0 ** denominator_exponent
    return m, (1.0 / m) * np.sqrt(rf.Rint(u, 1.0) / i0)


class Curve:
    """Exact ``eta(m)`` of one rate form (Route B), with the dead-zone branch.

    A table of ``(u, m)`` built once brackets the centre concentration, which
    ``brentq`` then finds for each query, in ``logit(u)``.

    Parameters
    ----------
    rf : Power, LHHW or VolChange
        Rate form.
    nu : int, optional
        Size of the bracketing table.
    root_in : {"logit", "u"}, optional
        Variable of the root search. ``"u"`` resolves ``u`` only to an
        absolute tolerance; it exists for the page's break table.

    Attributes
    ----------
    m_lim : float
        Largest modulus in the table, at ``u = U_MIN``. Above it
        :meth:`eta` returns ``1/m``: exact on a dead-zone branch, and within
        ``1e-9`` of ``eta`` elsewhere.
    """

    def __init__(self, rf, nu=44, root_in="logit"):
        self.rf, self.root_in = rf, root_in
        us = np.unique(np.concatenate([np.geomspace(U_MIN, 0.5, nu // 2),
                                       1.0 - np.geomspace(1e-9, 0.5, nu // 2)]))
        ms = np.array([mb_eta(rf, u)[0] for u in us])
        order = np.argsort(ms)
        self._m_tab, self._u_tab = ms[order], us[order]
        self.m_lim = float(ms[0])

    def eta(self, m):
        """Effectiveness factor at generalised modulus ``m``."""
        if m >= self.m_lim * (1.0 - 1e-12):
            return 1.0 / m
        j = int(np.searchsorted(self._m_tab, m))
        lo = self._u_tab[j] if j < len(self._u_tab) else U_MIN
        hi = self._u_tab[j - 1] if j > 0 else 1.0 - 1e-12

        # root-find in v = logit(u), so that u -> 0 (large m) and u -> 1 (small
        # m) are both resolved to a relative tolerance
        to_u, to_v = (expit, logit) if self.root_in == "logit" else (float, float)

        def f(v):
            return mb_eta(self.rf, to_u(v))[0] - m

        lo, hi = to_v(lo), to_v(hi)
        if not (f(lo) > 0 > f(hi) or f(lo) < 0 < f(hi)):
            lo, hi = to_v(U_MIN), to_v(1.0 - 1e-12)
        v = brentq(f, lo, hi, xtol=1e-12, rtol=8.9e-16)
        return mb_eta(self.rf, to_u(v))[1]


class ZeroOrder:
    """Bischoff eq. (22), attributed by him to Wheeler: the exact zero-order curve."""

    rf = Power(0.0)
    m_lim = 1.0

    def eta(self, m):
        """``1`` below ``m = 1``, ``1/m`` above."""
        return 1.0 if m < 1.0 else 1.0 / m


def mstar_power(n):
    """Dead-zone onset ``m* = (n+1)/(1-n)`` of ``r = c^n``, ``n < 1``.

    Derived from eq. (14) with ``C_L = 0``; not printed in the paper.
    """
    return (n + 1.0) / (1.0 - n)


# --------------------------------------------- Route C: eqs. (20)-(28)
SQ3 = np.sqrt(3.0)
K15SQ = (2.0 - SQ3) / 4.0     # k^2 with k = sin 15 deg
K45SQ = 0.5                   # k^2 with k = sin 45 deg


def phi_from_cos(t):
    """Amplitude ``phi`` of eqs. (23)-(24) from ``cos phi = (sqrt3+1-t)/(sqrt3-1+t)``.

    Equal to the printed ``sin phi = sqrt(4 sqrt3 (t-1))/(sqrt3-1+t)``, but it
    carries ``phi`` past ``pi/2`` (at ``t = 1 + sqrt3``).
    """
    return np.arccos((SQ3 + 1.0 - t) / (SQ3 - 1.0 + t))


def _phi_from_sin(t):
    """The printed ``sin phi`` form through ``arcsin``: stops at ``pi/2``."""
    s = np.sqrt(4.0 * SQ3 * (t - 1.0)) / (SQ3 - 1.0 + t)
    return np.arcsin(np.minimum(s, 1.0))


def m_eta_n1(u):
    """Eq. (21) region: ``m = arccosh(1/u)``, ``eta = tanh(m)/m``."""
    m = np.arccosh(1.0 / u)
    return m, np.tanh(m) / m


def m_eta_n2(u, k2=K15SQ, branch="cos"):
    """Eqs. (23)-(24), second order.

    Parameters
    ----------
    u : float
        Centre concentration ``C_L/C_o``.
    k2 : float, optional
        Squared elliptic modulus; the paper's is ``sin^2 15 deg``.
    branch : {"cos", "sin"}, optional
        ``"cos"`` continues ``phi`` past ``pi/2``; ``"sin"`` is the printed
        ``sin phi`` through ``arcsin``, which does not.

    Returns
    -------
    m, eta : float
    """
    t = 1.0 / u
    phi = phi_from_cos(t) if branch == "cos" else _phi_from_sin(t)
    m = 1.5 * np.sqrt(t) * ellipkinc(phi, k2) / 3.0 ** 0.25
    return m, (1.0 / m) * np.sqrt(1.0 - u ** 3)


def m_eta_n3(u):
    """Eqs. (25)-(26), third order."""
    phi = np.arcsin(np.sqrt(1.0 - u ** 2))
    m = np.sqrt(2.0) / u * ellipkinc(phi, K45SQ)
    return m, (1.0 / m) * np.sqrt(1.0 - u ** 4)


def _bracket27(phi):
    s, c = np.sin(phi), np.cos(phi)
    return ((1.0 + SQ3) * ellipkinc(phi, K15SQ) - 2.0 * SQ3 * ellipeinc(phi, K15SQ)
            + 2.0 * SQ3 * s * np.sqrt(1.0 - K15SQ * s * s) / (1.0 + c))


def m_eta_nhalf(u, printed=False):
    """Eqs. (27)-(28), order one-half.

    Parameters
    ----------
    u : float
        Centre concentration ``C_L/C_o``.
    printed : bool, optional
        ``False``: the reduction of eq. (20) through ``y = sqrt(t)``, with
        prefactor ``(C_L/C_o)^(1/4)`` and ``phi`` built on ``(C_o/C_L)^(1/2)``
        (the page's amendment). ``True``: exactly as printed, prefactor
        ``sqrt(C_L/C_o)`` and ``phi`` from the ``n = 2`` where-block, built
        on ``C_o/C_L``.

    Returns
    -------
    m, eta : float
    """
    if printed:
        m = 1.5 * np.sqrt(u) / 3.0 ** 0.25 * _bracket27(phi_from_cos(1.0 / u))
    else:
        m = 1.5 * u ** 0.25 / 3.0 ** 0.25 * _bracket27(phi_from_cos(np.sqrt(1.0 / u)))
    return m, (1.0 / m) * np.sqrt(1.0 - u ** 1.5)


def m_eq20(n, u):
    """Eq. (20): ``m`` of ``r = c^n`` at ``u`` by one quadrature."""
    t = 1.0 / u
    val = quad(lambda s: 2.0 * s / np.sqrt(np.expm1((n + 1.0) * np.log1p(s * s))),
               0.0, np.sqrt(t - 1.0), epsabs=1e-13, epsrel=1e-11, limit=400)[0]
    return (n + 1.0) / 2.0 * t ** ((n - 1.0) / 2.0) * val


# ------------------------------------------------ Route A: pymrm slab
# a dc/dn + b c = d with the OUTWARD normal n (pitfalls P1)
SYMMETRY = {"a": 1.0, "b": 0.0, "d": 0.0}   # dc/dn = 0: the pellet centre (x = 0)
SURFACE = {"a": 0.0, "b": 1.0, "d": 1.0}    # c = 1: the surface (x = 1)


class Slab:
    """``d/dx[d(c) dc/dx] = Phi^2 r(c)`` on ``(0, 1)``, ``c'(0) = 0``, ``c(1) = 1``.

    Parameters
    ----------
    n : int, optional
        Number of cells.
    nu : int, optional
        Geometry factor of ``construct_div`` (0 slab). Other values exist for
        the page's break table.

    Attributes
    ----------
    jac_diff : scipy.sparse.csc_array
        Constant operator ``div @ (-grad)`` for ``d = 1``.
    g_const : numpy.ndarray
        Its constant boundary contribution, shape ``(n, 1)``.
    weights : numpy.ndarray
        Cell volumes divided by the total volume.
    """

    def __init__(self, n=400, nu=0):
        self.n, self.nu = n, nu
        self.shape = (n, 1)                        # one field: (n, 1), never (n,)
        self.x_f = np.linspace(0.0, 1.0, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.bc = (SYMMETRY, SURFACE)              # x = 0: dc/dx = 0; x = 1: c = 1
        self.grad_mat, self.grad_bc = construct_grad(self.shape, self.x_f, self.x_c,
                                                     self.bc, axis=0)
        # nu = 0: slab (Cartesian)
        self.div_mat = construct_div(self.shape, self.x_f, nu=nu, axis=0)
        self.jac_diff = (self.div_mat @ (-self.grad_mat)).tocsc()
        self.g_const = (self.div_mat @ (-self.grad_bc)).toarray().reshape(-1, 1)
        self.numjac = NumJac(self.shape)
        # the full residual with d(c) on the faces reads both neighbours of a cell
        self.numjac_full = NumJac(self.shape, axes_diagonals=[0])
        volume = np.diff(self.x_f ** (nu + 1)) / (nu + 1)
        self.weights = volume / volume.sum()

    def face_concentration(self, c):
        """``c`` on the faces: linear between cells, the boundary values from the bc."""
        c_face = interp_cntr_to_stagg(c, self.x_f, self.x_c, axis=0)
        c_lo, _, c_hi, _ = compute_boundary_values(c, self.x_f, self.x_c, self.bc, axis=0)
        c_face[0], c_face[-1] = c_lo.ravel()[0], c_hi.ravel()[0]
        return c_face

    def full_residual(self, c, rf, Phi):
        """``div(-d(c_face) grad c) + Phi^2 r(c)``, shape ``(n, 1)``."""
        flux = -rf.d(self.face_concentration(c)) * (self.grad_mat @ c + self.grad_bc)
        return self.div_mat @ flux + Phi ** 2 * rf.r(c)

    def residual(self, c_flat, rf, Phi, eps=0.0, exact_jacobian=True):
        """Residual and Jacobian for :func:`pymrm.newton`.

        Constant diffusivity: the constant operator plus the reaction, whose
        Jacobian is ``NumJac``. Variable diffusivity: one ``NumJac`` over the
        full residual, ``d(c)`` on the faces included. ``eps > 0``: the
        regularised rate ``r(c + eps) - r(eps)`` with its exact Jacobian, or,
        with ``exact_jacobian=False``, with ``NumJac`` (constant diffusivity
        only).
        """
        c = c_flat.reshape(self.shape)
        if eps > 0.0:
            assert not rf.variable_d, "the regularised rate assumes a constant diffusivity"
        if eps > 0.0 and not exact_jacobian:
            g_r, jac_r = self.numjac(
                lambda cc: Phi ** 2 * (rf.r(np.maximum(cc, 0.0) + eps) - rf.r(eps)), c)
            return self.jac_diff @ c + self.g_const + g_r.reshape(-1, 1), self.jac_diff + jac_r
        if eps > 0.0:
            cp = np.maximum(c, 0.0)
            g = (self.jac_diff @ c + self.g_const
                 + Phi ** 2 * (rf.r(cp + eps) - rf.r(eps)))
            jac = self.jac_diff + construct_coefficient_matrix(Phi ** 2 * rf.dr(cp + eps),
                                                               shape=self.shape)
            return g, jac
        if rf.variable_d:
            g, jac = self.numjac_full(lambda cc: self.full_residual(cc, rf, Phi), c)
            return g.reshape(-1, 1), jac
        g_r, jac_r = self.numjac(lambda cc: Phi ** 2 * rf.r(cc), c)
        return self.jac_diff @ c + self.g_const + g_r.reshape(-1, 1), self.jac_diff + jac_r

    def solve(self, rf, Phi, eps=0.0, check=True, exact_jacobian=True, **newton_kw):
        """Steady profile from ``c = 1``.

        Parameters
        ----------
        rf : Power, LHHW or VolChange
            Rate form.
        Phi : float
            Thiele modulus of the normalised problem.
        eps : float, optional
            Regularisation of the rate (see :meth:`residual`).
        check : bool, optional
            Raise ``RuntimeError`` unless newton reports success.
        exact_jacobian : bool, optional
            See :meth:`residual`.
        **newton_kw
            Options of :func:`pymrm.newton`; default ``tol=1e-12``,
            ``maxfev=60``, ``bounds=(0, inf)`` (``c^n`` with ``n < 1`` is
            undefined below zero) and ``stagnation=True`` (on fine grids the
            step test meets the round-off floor near ``1e-12``).

        Returns
        -------
        c : numpy.ndarray
            Shape ``(n, 1)``.
        sol : scipy.optimize.OptimizeResult
        """
        kw = {"tol": 1e-12, "maxfev": 60, "bounds": (0.0, np.inf), "stagnation": True,
              **newton_kw}
        sol = newton(lambda x: self.residual(x, rf, Phi, eps, exact_jacobian),
                     np.ones(self.shape), **kw)
        if check and not sol.success:
            raise RuntimeError(f"slab solve failed ({rf.label}, Phi = {Phi:g}): {sol.message}")
        return sol.x.reshape(self.shape), sol

    def eta_of(self, c, rf, Phi, read="volume"):
        """Effectiveness factor of a profile.

        ``"volume"``: the volume average of the rate. ``"face"``: the surface
        flux from :func:`pymrm.compute_boundary_values`. ``"cell"``: the
        surface flux from the last cell centre and the surface value, a
        two-point gradient.
        """
        if read == "volume":
            return float(self.weights @ rf.r(c).ravel())
        if read == "face":
            _, _, _, grad_hi = compute_boundary_values(c, self.x_f, self.x_c, self.bc, axis=0)
            grad = float(np.ravel(grad_hi)[0])
        else:
            grad = (1.0 - float(c[-1, 0])) / (self.x_f[-1] - self.x_c[-1])
        return (self.nu + 1) * float(rf.d(1.0)) * grad / Phi ** 2

    def eta(self, rf, m, eps=0.0, read="volume"):
        """Effectiveness factor at generalised modulus ``m``.

        ``Phi = m sqrt(2 I0)`` inverts eq. (13), ``I0 = int_0^1 d r dc``.
        """
        Phi = m * np.sqrt(2.0 * rf.Rint(0.0, 1.0))
        c, _ = self.solve(rf, Phi, eps)
        return self.eta_of(c, rf, Phi, read)


# ------------------------------------------- the finite-cylinder audit
class Cyl2D:
    """Finite cylinder, radial (``nu = 1``) times axial (Cartesian), ``c = 1`` on
    the mantle and the end faces, symmetry on the axis and the midplane.

    Parameters
    ----------
    nr : int, optional
        Radial cells; the axial count is ``max(round(nr Z), 20)``.
    Z : float, optional
        Aspect ratio ``H/(2 R_p)``: the half-length in units of the radius.
    nu_radial : int, optional
        Geometry factor of the radial axis (1). Other values exist for the
        page's break table.
    """

    def __init__(self, nr=100, Z=1.0, nu_radial=1):
        nz = max(int(round(nr * Z)), 20)
        self.shape = (nr, nz, 1)                   # spatial axes first, one field
        jac, g, measures = 0, 0, []
        for ax, (npts, ext, nu) in enumerate(((nr, 1.0, nu_radial), (nz, Z, 0))):
            x_f = np.linspace(0.0, ext, npts + 1)
            x_c = 0.5 * (x_f[:-1] + x_f[1:])
            # low end: dc/dn = 0 (axis, midplane); high end: c = 1 (exposed face)
            grad_mat, grad_bc = construct_grad(self.shape, x_f, x_c, (SYMMETRY, SURFACE),
                                               axis=ax)
            # radial axis nu = 1 (cylindrical), axial nu = 0
            div_mat = construct_div(self.shape, x_f, nu=nu, axis=ax)
            jac = jac + div_mat @ (-grad_mat)
            g = g + div_mat @ (-grad_bc)
            measures.append(np.diff(x_f ** (nu + 1)) / (nu + 1))
        self.jac_diff = jac.tocsc()
        self.g_const = g.toarray().reshape(-1, 1)
        self.numjac = NumJac(self.shape)
        w = np.multiply.outer(measures[0], measures[1])
        self.weights = (w / w.sum()).ravel()

    def eta(self, rf, phi):
        """Effectiveness factor at Thiele modulus ``phi = R_p sqrt(k_v C_o^(n-1)/D)``."""
        def residual(c_flat):
            c = c_flat.reshape(self.shape)
            g_r, jac_r = self.numjac(lambda cc: phi ** 2 * rf.r(cc), c)
            return (self.jac_diff @ c_flat.reshape(-1, 1) + self.g_const
                    + g_r.reshape(-1, 1), self.jac_diff + jac_r)

        sol = newton(residual, np.ones(self.shape), tol=1e-12, maxfev=40,
                     bounds=(0.0, np.inf), stagnation=True)
        if not sol.success:
            raise RuntimeError(f"cylinder solve failed ({rf.label}, phi = {phi:g}): "
                               f"{sol.message}")
        return float(self.weights @ rf.r(sol.x).ravel())


# ------------------------------------------------------- Pan and Zhu (1998)
def mu1_series(Z, M=200, N=400, bessel_order=0):
    """Their eq. (26), with ``alpha_m`` the zeros of ``J0`` and
    ``theta_n = (2n+1) pi/(2Z)`` (eq. 24).

    ``M``, ``N`` are the numbers of terms; ``bessel_order`` exists for the
    page's break table.
    """
    am = jn_zeros(bessel_order, M)[:, None]
    k = 2 * np.arange(N)[None, :] + 1
    th = k * np.pi / (2.0 * Z)
    return 32.0 / np.pi ** 2 * float(np.sum(1.0 / (am ** 2 * k ** 2 * (am ** 2 + th ** 2))))


def mu2_eq30(Z):
    """Their eq. (30), ``mu2 = 2 + 1/Z``."""
    return 2.0 + 1.0 / Z


def eta_series_1st(phi, Z=1.0, M=200, N=400, bessel_order=0):
    """Exact first-order effectiveness factor of the finite cylinder, their eq. (37).

    The same double series as Aris's case (iv), which page B1.2 checks against
    a 2-D solve.
    """
    am = jn_zeros(bessel_order, M)[:, None]
    k = 2 * np.arange(N)[None, :] + 1
    th = k * np.pi / (2.0 * Z)
    return 1.0 - 32.0 * phi ** 2 / np.pi ** 2 * float(
        np.sum(1.0 / (k ** 2 * am ** 2 * (am ** 2 + th ** 2 + phi ** 2))))


def pz_coeffs(Z, Q1, Q2, sign405=+1, mu1=None):
    """``q1, q2, q3`` of their eq. (32) from the two matching conditions they state.

    ``phi -> 0`` to their eq. (25) and ``phi -> inf`` to their eq. (29). With
    ``sign405 = +1`` the coefficients satisfy both conditions (asserted);
    ``-1`` is their eq. (36) as printed.
    """
    a = (mu1_series(Z) if mu1 is None else mu1) * Q1
    b = mu2_eq30(Z) * Q2
    q1 = b
    q2 = 27.0 - 135.0 * a - 6.0 * b
    q3 = -54.0 + sign405 * 405.0 * a + 9.0 * b
    if sign405 == +1:
        assert abs(q1 / 3 + q2 / 9 + q3 / 27 - 1.0) < 1e-12      # eta -> 1 at phi = 0
        assert abs(3 * q1 + 2 * q2 + q3 - 135.0 * a) < 1e-9      # slope of eq. (25)
    return q1, q2, q3


def eta_pz(phi, Z, Q1, Q2, sign405=+1, mu1=None):
    """Their eq. (32): ``q1 e + q2 e^2 + q3 e^3`` with ``e = (coth phi - 1/phi)/phi``."""
    q1, q2, q3 = pz_coeffs(Z, Q1, Q2, sign405, mu1)
    e = (1.0 / np.tanh(phi) - 1.0 / phi) / phi
    return q1 * e + q2 * e ** 2 + q3 * e ** 3


# --------------------------------------------------------- Table 1
def eq33(w, root=0.5):
    """Bischoff eq. (33): ``m/M`` for volume change, ``|w|/sqrt(2 (w - ln(1 + w)))``.

    ``root`` is the exponent of the bracket; it exists for the page's break table.
    """
    return abs(w) / np.sqrt(2.0) / (w - np.log(w + 1.0)) ** root


def eq34(w, exponent=0.7):
    """Hawthorn's factor, Bischoff eq. (34): ``(w/ln(1 + w))^0.7``."""
    return (w / np.log(1.0 + w)) ** exponent
