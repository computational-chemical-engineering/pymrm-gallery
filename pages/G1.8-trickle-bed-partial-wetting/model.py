"""Partially wetted catalyst pellet with a volatile reactant (Herskowitz and Smith 1983).

Model code for gallery page G1.8. Two routes to the same quantity:

- the paper's closed form: Table 1's sphere row (page 4) with the closures of
  Eqs. 20 and 21 (page 8), in :func:`eta_o_table1`, :func:`chi_table1` and the
  collapsed form :func:`chi_collapsed` for ``alpha_gs -> inf`` and ``C*_L = 1``;
- a pymrm finite-volume solve of the particle problem, Eq. 3 with the Robin
  conditions 4a/4b, combined by the paper's Eq. 6 superposition, in
  :class:`WettedSphere`.

The root finders :func:`reconstruct_L` and :func:`alpha_gs_for` invert the
closed form on brackets that are justified by monotonicity and widened until
they bracket. No plotting and no data loading here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from pymrm import NumJac, compute_boundary_values, construct_div, construct_grad, newton


# ---------------------------------------------------------------------------
# The paper's closures and closed form
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Closures:
    """The two closures Figure 6 uses (page 8), with their printed constants.

    Parameters
    ----------
    fe_coef, fe_expo : float
        Eq. 20, ``f_e = fe_coef * L_m**fe_expo`` (0.77 and 0.1 as printed).
    c21 : tuple of float
        Eq. 21, ``1/alpha_gLs = 1/(c21[0] L_m**c21[1]) + 1/(c21[2] L_m**c21[3])``
        (1.05, 0.3, 6.91, 0.6 as printed): gas-liquid and liquid-solid
        resistances in series, so the reciprocals add.

    The defect-injection rows of the page pass a ``Closures`` with one constant
    mis-read, so that every computation that uses the closures sees the defect.
    """

    fe_coef: float = 0.77
    fe_expo: float = 0.1
    c21: tuple = (1.05, 0.3, 6.91, 0.6)

    def f_e(self, L_m):
        """Wetting efficiency, Eq. 20."""
        return self.fe_coef * np.asarray(L_m, float) ** self.fe_expo

    def alpha_gLs(self, L_m):
        """Gas-liquid plus liquid-solid Biot number on the V/S length, Eq. 21."""
        L = np.asarray(L_m, float)
        c = self.c21
        return 1.0 / (1.0 / (c[0] * L ** c[1]) + 1.0 / (c[2] * L ** c[3]))

    def L_full_wetting(self):
        """The ``L_m`` at which Eq. 20 reaches ``f_e = 1`` (and ``chi = 0``)."""
        return (1.0 / self.fe_coef) ** (1.0 / self.fe_expo)


PAPER = Closures()


def eta_sphere(phi, length_ratio=3.0):
    """Sphere effectiveness factor with ``phi`` defined on ``R / length_ratio``.

    Table 1's sphere row contains ``3 phi^2 / (3 phi coth 3 phi - 1)``, which is
    ``1/eta_sphere``: with the paper's length ``V_p/S_x = R/3`` the Thiele
    modulus on the radius is ``lam = 3 phi`` and
    ``eta = 3 (lam coth lam - 1) / lam^2``.

    Parameters
    ----------
    phi : float or array_like
        Thiele modulus on the paper's length scale.
    length_ratio : float, optional
        ``R`` over the length ``phi`` is defined on; 3 for ``V/S`` (the paper).
        Any other value is a defect the page injects on purpose.
    """
    lam = length_ratio * np.asarray(phi, float)
    return 3.0 * (lam / np.tanh(lam) - 1.0) / lam ** 2


def eta_o_table1(phi, L_m, alpha_gs=np.inf, C_star_L=1.0, closures=PAPER, length_ratio=3.0):
    """Overall effectiveness factor, Table 1 sphere row (page 4), m = 1.

    The paper's page 8 cites this as "Table 2"; it is the third row of Table 1.
    A fraction ``f_e`` of the surface is wetted (Robin with ``alpha_gLs`` and
    bulk value ``C*_L``), the rest is dry (Robin with ``alpha_gs`` and bulk
    value 1); the two whole-surface particles are superposed (Eq. 6).
    """
    inv_eta = 1.0 / eta_sphere(phi, length_ratio)
    a, fe = closures.alpha_gLs(L_m), closures.f_e(L_m)
    wet = fe * C_star_L / (phi ** 2 / a + inv_eta)
    dry = (1.0 - fe) / (phi ** 2 / alpha_gs + inv_eta)
    return wet + dry


def chi_table1(phi, L_m, alpha_gs=np.inf, C_star_L=1.0, closures=PAPER, length_ratio=3.0):
    """Eq. 19, ``chi = eta_o / eta_o(f_e = 1) - 1``, on Table 1's sphere row."""
    # denominator: the same row at f_e = 1 and C*_L = 1, the wetted term alone
    eta_full = 1.0 / (phi ** 2 / closures.alpha_gLs(L_m) + 1.0 / eta_sphere(phi, length_ratio))
    return eta_o_table1(phi, L_m, alpha_gs, C_star_L, closures, length_ratio) / eta_full - 1.0


def chi_collapsed(phi, L_m, closures=PAPER, length_ratio=3.0):
    """``chi`` for ``alpha_gs -> inf`` and ``C*_L = 1``: ``(1 - f_e) eta_s phi^2 / alpha_gLs``.

    An algebraic identity with ``chi_table1(..., alpha_gs=inf)``: agreement of
    the two checks the collapse step, never the physics.
    """
    return ((1.0 - closures.f_e(L_m)) * eta_sphere(phi, length_ratio) * phi ** 2
            / closures.alpha_gLs(L_m))


# ---------------------------------------------------------------------------
# Inverting the closed form
# ---------------------------------------------------------------------------

def _widen(fun, x, sign, factor, limit):
    """Multiply ``x`` by ``factor`` until ``sign * fun(x) > 0``; raise past ``limit``."""
    while not sign * fun(x) > 0:
        x *= factor
        if not min(limit, 1.0 / limit) <= x <= max(limit, 1.0 / limit):
            raise ValueError(f"no sign change found up to {limit:g}")
    return x


def reconstruct_L(chi_target, phi=10.0, closures=PAPER, length_ratio=3.0):
    """The ``L_m`` at which :func:`chi_collapsed` equals ``chi_target``.

    On ``(0, L_full_wetting)`` both ``1 - f_e`` and ``1/alpha_gLs`` are positive
    and decrease with ``L_m`` (positive exponents), so ``chi`` falls
    monotonically from ``+inf`` to 0 and a positive target has exactly one root.
    The upper end of the bracket is ``L_full_wetting`` (``chi = 0``); the lower
    end starts three decades below it and is widened tenfold until it brackets.

    Parameters
    ----------
    chi_target : float
        A positive ``chi``.
    phi : float, optional
        Thiele modulus at which to match.
    closures : Closures, optional
        Eqs. 20 and 21.
    length_ratio : float, optional
        Passed to :func:`eta_sphere`.

    Returns
    -------
    float
        The root, to a relative accuracy of about 1e-14.
    """
    if not chi_target > 0:
        raise ValueError(f"chi_target must be positive, got {chi_target!r}")
    hi = closures.L_full_wetting()

    def f(L):
        return float(chi_collapsed(phi, L, closures, length_ratio)) - chi_target

    lo = _widen(f, 1e-3 * hi, +1, 0.1, 1e-30)     # f(hi) = -chi_target < 0
    # xtol relative to the lower end, so a root far below 1 keeps its relative accuracy
    return brentq(f, lo, hi, xtol=1e-14 * lo, rtol=1e-14)


def alpha_gs_for(chi_target, phi, L_m, closures=PAPER):
    """The finite gas-solid Biot number at which :func:`chi_table1` equals ``chi_target``.

    ``chi_table1`` increases monotonically with ``alpha_gs`` (only the dry term
    depends on it, and that term increases), from ``f_e - 1`` at 0 to
    :func:`chi_collapsed` at infinity, so a root exists if and only if the
    target lies between the two, and it is unique. The bracket starts at
    ``[1, 10]`` and each end is widened tenfold until it brackets.

    The argument needs ``f_e < 1``: at ``f_e >= 1`` there is no dry surface and
    the dependence on ``alpha_gs`` vanishes or reverses.

    Raises
    ------
    ValueError
        When ``f_e(L_m) >= 1``, or when the target is at or above the
        ``alpha_gs -> inf`` limit (no root).
    """
    fe = float(closures.f_e(L_m))
    if not fe < 1.0:
        raise ValueError(f"f_e = {fe:.4g} >= 1 at L_m = {L_m:g}: no dry surface, so no "
                         "alpha_gs to solve for")
    if not chi_target < float(chi_collapsed(phi, L_m, closures)):
        raise ValueError("chi_target is not below the alpha_gs -> inf limit: no root")

    def f(a):
        return float(chi_table1(phi, L_m, alpha_gs=a, closures=closures)) - chi_target

    hi = _widen(f, 10.0, +1, 10.0, 1e30)       # f(hi) > 0
    lo = _widen(f, 1.0, -1, 0.1, 1e-30)        # f(lo) < 0
    return brentq(f, lo, hi, xtol=1e-10, rtol=1e-14)


# ---------------------------------------------------------------------------
# The pymrm route: Eq. 3 with a Robin surface condition
# ---------------------------------------------------------------------------

class WettedSphere:
    """Isothermal spherical pellet, Eq. 3 with one Robin surface condition.

    Grid ``u = r/R`` on ``[0, 1]``. The paper's length scale is
    ``V_p/S_x = R/3``, so ``d/dn|paper = (1/3) d/du`` and
    ``lap|paper = (1/9) lap_u``: the interior equation is
    ``lap_u C = (3 phi)^2 C^m`` and the Robin condition
    ``(1/3) dC/du + alpha C = alpha C*``.

    Parameters
    ----------
    n_u : int, optional
        Number of cells.
    nu : int, optional
        Geometry of ``construct_div``: 0 slab, 1 cylinder, 2 sphere (the model).
    n_scale : float, optional
        ``R`` over the length of the paper's normal derivative and Thiele
        modulus, 3 for ``V/S``. Any other value is a defect the page injects.
    tol : float, optional
        ``newton`` step tolerance, above the round-off floor of the unknowns.
    """

    def __init__(self, n_u=800, nu=2, n_scale=3.0, tol=1e-10):
        self.n_u, self.shape, self.nu, self.n_scale, self.tol = n_u, (n_u, 1), nu, n_scale, tol
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)                  # faces
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])            # centres
        self.div = construct_div(self.shape, self.u_f, nu=nu)      # nu: 0 slab, 1 cylinder, 2 sphere
        vf = self.u_f ** (nu + 1) / (nu + 1.0)                     # cell volumes for this nu
        self.dv, self.v_tot = np.diff(vf), vf[-1]
        self.numjac = NumJac(self.shape)   # (n_u, 1): the last axis is the field, pointwise source
        self._ops = {}
        self.last_result = None

    def bc(self, alpha, C_star):
        """Boundary dictionaries for a surface Biot number ``alpha`` (OUTWARD normal)."""
        s = self.n_scale
        # u = 0 (centre): symmetry, dC/du = 0; the outward normal is -u, so
        #   a dC/dn + b C = d with a = 1, b = 0, d = 0
        left = {"a": 1.0, "b": 0.0, "d": 0.0}
        if np.isinf(alpha):
            # u = 1 (surface), alpha -> inf: C = C*
            right = {"a": 0.0, "b": 1.0, "d": C_star}
        else:
            # u = 1 (surface): dC/dn|paper = alpha (C* - C), with dC/dn|paper = (1/s) dC/du
            #   and the outward normal +u:  dC/du + s alpha C = s alpha C*
            right = {"a": 1.0, "b": s * alpha, "d": s * alpha * C_star}
        return (left, right)

    def operators(self, alpha, C_star):
        """Laplacian and its boundary vector, assembled once per surface condition."""
        key = (float(alpha), float(C_star))
        if key not in self._ops:
            g, gb = construct_grad(self.shape, self.u_f, self.u_c, self.bc(alpha, C_star))
            lap = (self.div @ g).tocsc()
            lap_bc = (self.div @ gb).toarray().reshape(-1, 1)
            self._ops[key] = (lap, lap_bc)
        return self._ops[key]

    def residual(self, y, phi, alpha, C_star=1.0, m=1.0):
        """Residual ``lap_u C - (s phi)^2 C^m`` and its Jacobian."""
        lap, lap_bc = self.operators(alpha, C_star)
        lam2 = (self.n_scale * phi) ** 2

        def source(c):
            return -lam2 * c ** m          # c >= 0 is kept by newton's bounds

        g_s, j_s = self.numjac(source, y.reshape(self.shape))
        g = lap @ y.reshape((-1, 1)) + lap_bc + np.asarray(g_s).reshape((-1, 1))
        return g, lap + j_s

    def solve(self, phi, alpha, C_star=1.0, m=1.0, maxfev=60):
        """Steady concentration profile, shape ``(n_u, 1)``.

        Raises
        ------
        RuntimeError
            When ``newton`` returns without convergence (a breakdown raises
            ``pymrm.NewtonBreakdown`` itself).
        """
        guess = np.full(self.shape, min(1.0, C_star))
        # bounds: C**m is undefined below 0 for non-integer m
        r = newton(lambda y: self.residual(y, phi, alpha, C_star, m), guess,
                   tol=self.tol, maxfev=maxfev, bounds=(0.0, None))
        self.last_result = r
        if not r.success:
            raise RuntimeError(f"WettedSphere.solve at phi = {phi:g}, alpha = {alpha:g}, "
                               f"m = {m:g}: status {r.status}, {r.message}")
        return r.x.reshape(self.shape)

    def eta(self, y, m=1.0):
        """``eta = V^-1 int C^m dV`` (Eq. 2a), the volume average of the local rate."""
        return float(np.sum(np.ravel(y) ** m * self.dv) / self.v_tot)

    def eta_surface(self, y, phi, alpha, C_star=1.0, surface="face"):
        """``eta`` from the surface flux, ``eta = (dC/dn|paper) / phi^2``.

        Integrating Eq. 3 over the pellet turns the volume average of the rate
        into the flux through the surface; in grid units
        ``eta = (dC/du at u = 1) / (s phi^2)`` with ``s = n_scale``. With
        ``surface="face"`` the gradient is the face value the operator
        transports, from ``compute_boundary_values`` (pitfall P8); ``"cell"``
        takes it from the last two cell centres, the defect the page injects.
        """
        if surface == "face":
            _, grad = compute_boundary_values(y, self.u_f, self.u_c, self.bc(alpha, C_star)[1],
                                              bound_id=1)
            grad = float(np.ravel(grad)[0])                        # along +u, the outward normal
        else:
            yv = np.ravel(y)
            grad = float((yv[-1] - yv[-2]) / (self.u_c[-1] - self.u_c[-2]))
        # the volume integral of u^2 du over the unit sphere is v_tot = 1/3
        return grad / (self.n_scale ** 2 * phi ** 2 * self.v_tot) * (self.u_f[-1] ** self.nu)

    def chi(self, phi, L_m, m=1.0, alpha_gs=np.inf, C_star_L=1.0, closures=PAPER):
        """Eq. 19 from three finite-volume solves and the Eq. 6 superposition."""
        a, fe = float(closures.alpha_gLs(L_m)), float(closures.f_e(L_m))
        eta_wet = self.eta(self.solve(phi, a, C_star_L, m), m)    # wetted path
        eta_dry = self.eta(self.solve(phi, alpha_gs, 1.0, m), m)  # dry path
        eta_ref = self.eta(self.solve(phi, a, 1.0, m), m)         # f_e = 1, C*_L = 1
        eta_o = fe * eta_wet + (1.0 - fe) * eta_dry               # Eq. 6
        return eta_o / eta_ref - 1.0
