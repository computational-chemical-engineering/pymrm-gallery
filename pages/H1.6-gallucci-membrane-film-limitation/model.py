"""A Sieverts membrane behind a bed-to-wall film: local elasticities and a membrane module.

Model code of gallery page H1.6. The review the page tests (Gallucci et al., 2013) prints no
model; these are the page's standard resistances-in-series relations, scaled by the feed H2
partial pressure ``p_in``:

- local problem (E4, E5): ``pi = p_b/p_in`` (bulk), ``theta = p_p/p_in`` (permeate),
  ``beta = Q_s/(k_p sqrt(p_in))`` (membrane over film conductance),
  ``u = sqrt(p_w/p_in)`` (membrane surface); ``u^2 + beta u = pi + beta sqrt(theta)``;
  flux ``j = J/(Q_s sqrt(p_in)) = u - sqrt(theta)``;
- flux elasticities (E6): ``eps_Q = 2u/(2u + beta)``, ``eps_k = beta/(2u + beta)``;
- figure conversions (E7 to E10, E15);
- the module (E11 to E14): plug-flow retentate with H2 fraction ``f`` and pure-H2 permeate
  ``g`` on ``xi = A/A_tot``, ``df/dxi = -N j``, ``dg/dxi = +N j``.

Contents
--------
wall_root, local_flux, film_share, elasticities
    The local series problem.
eps_q_vacuum, theta_star, lambda_star, depletion_vacuum
    Fig. 15 conversions (vacuum permeate; the permeate ratio that makes the film limiting).
share_elasticities, eps_q_bound
    Fig. 16: elasticities from the shares of the pressure drop, and the bound.
eps_k_from_eta
    p. 57: an effectiveness factor read as a series film.
wall_factor, centre_depletion, cupmix_depletion, fit_profile
    Fig. 15: the radial profile normalised at the first cell.
pi_of_f, pinch, module_number, area_elasticity, n_membrane_vacuum, a_film, recovery_at
    The module by quadrature and closed forms (the second route).
MembraneModule
    The module as a pymrm finite-volume system solved by ``pymrm.newton``.
"""

from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import quad

from pymrm import (BlockLayout, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, newton)
from pymrm.checks import find_roots

__all__ = ["wall_root", "local_flux", "film_share", "elasticities", "eps_q_vacuum", "theta_star",
           "lambda_star", "depletion_vacuum", "share_elasticities", "eps_q_bound", "eps_k_from_eta",
           "wall_factor", "centre_depletion", "cupmix_depletion", "fit_profile", "pi_of_f", "pinch",
           "module_number", "area_elasticity", "n_membrane_vacuum", "a_film", "recovery_at",
           "MembraneModule"]


# ---------------------------------------------------------------- local series problem
def wall_root(pi, theta, beta, n=0.5):
    """Membrane-surface variable ``u = sqrt(p_w/p_in)`` (E4), cancellation-free.

    For ``n = 1`` (a linear membrane, the sensitivity of the exponent) the returned value is
    still ``sqrt(p_w/p_in)``, from ``p_w (1 + beta') = pi + beta' theta`` with
    ``beta' = beta`` read as the linear conductance ratio.
    """
    if n == 1.0:
        return np.sqrt((pi + beta * theta) / (1.0 + beta))
    c = pi + beta * np.sqrt(theta)
    return 2.0 * c / (beta + np.sqrt(beta * beta + 4.0 * c))


# codehealth: allow params -- the scaled groups of E4 and E5 plus the form selector
def local_flux(pi, theta, beta, form="ratio", n=0.5):
    """Scaled flux ``j`` (E5) in one of three algebraically equal forms.

    ``"ratio"`` ``(pi - theta)/(beta + u + sqrt(theta))`` (cancellation-free near the pinch),
    ``"membrane"`` ``u - sqrt(theta)``, ``"film"`` ``(pi - u^2)/beta``. ``n = 1``: a linear
    membrane, ``j = u^2 - theta``.
    """
    u = wall_root(pi, theta, beta, n)
    if n == 1.0:
        return {"ratio": (pi - theta) / (1.0 + beta), "membrane": u * u - theta,
                "film": (pi - u * u) / beta}[form]
    st = np.sqrt(theta)
    if form == "ratio":
        return (pi - theta) / (beta + u + st)
    if form == "membrane":
        return u - st
    if form == "film":
        return (pi - u * u) / beta
    raise ValueError(f"form {form!r} not recognised")


def film_share(pi, theta, beta):
    """Film share of the partial-pressure drop, ``(p_b - p_w)/(p_b - p_p)`` (E5)."""
    u = wall_root(pi, theta, beta)
    return beta / (beta + u + np.sqrt(theta))


def elasticities(pi, theta, beta):
    """Flux elasticities ``(eps_Q, eps_k)`` at fixed bulk and permeate (E6); they sum to 1."""
    u = wall_root(pi, theta, beta)
    return 2.0 * u / (2.0 * u + beta), beta / (2.0 * u + beta)


# ---------------------------------------------------------------- Fig. 15 conversions
def eps_q_vacuum(x2):
    """Membrane elasticity under vacuum permeate from ``x^2 = p_w/p_b`` (E7)."""
    return 2.0 * x2 / (1.0 + x2)


def theta_star(x2):
    """Permeate-to-bed H2 ratio at which ``eps_Q = 1/2`` (E7); the film limits above it.

    Raises
    ------
    ValueError
        When ``x2 <= 1/3``: the film limits at any permeate pressure.
    """
    if x2 <= 1.0 / 3.0:
        raise ValueError(f"x^2 = {x2:.4g} <= 1/3: eps_Q < 1/2 already under vacuum")
    return (3.0 * x2 - 1.0) ** 2 / (4.0 * x2)


def lambda_star(d1):
    """Permeability multiplier at which the vacuum film takes 2/3 of the drop (E8)."""
    return (2.0 / math.sqrt(3.0)) / (d1 / math.sqrt(1.0 - d1))


def depletion_vacuum(beta_b):
    """Wall depletion ``1 - p_w/p_b`` of the vacuum series problem at bulk ratio ``beta_b`` (E8)."""
    u = 2.0 / (beta_b + math.sqrt(beta_b * beta_b + 4.0))
    return 1.0 - u * u


# ---------------------------------------------------------------- Fig. 16 conversions
# codehealth: allow params -- the three shares of Fig. 16, the interface ratio and the exponent of E9
def share_elasticities(c, m, sigma, s, n=0.5):
    """``(eps_Q, eps_k, eps_support)`` from the shares of the total drop (E9).

    ``c`` film, ``m`` membrane, ``sigma`` support shares; ``s = p_s/p_w`` at the Pd/support
    interface; ``n`` the membrane exponent.
    """
    g = n * (1.0 - s) / (1.0 - s ** n)
    a = (c / m) * g
    b = (sigma / m) * g * s ** (n - 1.0)
    return 1.0 / (1.0 + a + b), a / (1.0 + a + b), b / (1.0 + a + b)


def eps_q_bound(c, m, n=0.5):
    """Upper bound of ``eps_Q`` for any interface ratio and support share (E9)."""
    return 1.0 / (1.0 + n * c / m)


# ---------------------------------------------------------------- p. 57 conversion
def eps_k_from_eta(eta, y):
    """Film elasticity of an effectiveness factor ``eta`` read as a series film (E10).

    ``y = sqrt(p_p/p_b)``.
    """
    x = y + eta * (1.0 - y)
    beta_b = (1.0 - x * x) / (x - y)
    return beta_b / (2.0 * x + beta_b)


# ---------------------------------------------------------------- Fig. 15 profile (E15)
def wall_factor(r0, rl, k=2.0):
    """Wall over last-cell depletion of a profile ``1 - D r^k`` normalised at ``r0``."""
    return (1.0 - r0 ** k) / (rl ** k - r0 ** k)


def centre_depletion(d_last, r0, rl, k=2.0):
    """Centre-to-wall depletion ``D`` from the last-cell depletion of the normalised profile."""
    d_w = d_last * wall_factor(r0, rl, k)
    return d_w / (1.0 - r0 ** k * (1.0 - d_w))


def cupmix_depletion(D, k=2.0):
    """Wall depletion relative to the area (cup-mix, plug flow) mean of ``1 - D r^k``."""
    w = 2.0 / (k + 2.0)
    return D * (1.0 - w) / (1.0 - D * w)


def fit_profile(r, v, r0, k=2.0):
    """Least-squares ``D`` of ``1 - v = D (r^k - r0^k)/(1 - D r0^k)`` and the residuals.

    The model is linear in ``q = D/(1 - D r0^k)``, so the fit is exact linear least squares.
    """
    x = np.asarray(r) ** k - r0 ** k
    y = 1.0 - np.asarray(v)
    q = float((x * y).sum() / (x * x).sum())
    return q / (1.0 + q * r0 ** k), (1.0 - q * x) - np.asarray(v)


# ---------------------------------------------------------------- module, second route
def pi_of_f(f, phi):
    """Bulk H2 ratio ``p_b/p_in`` at retained H2 fraction ``f`` (inert per H2 ``phi``)."""
    return f * (1.0 + phi) / (f + phi)


def pinch(theta, phi):
    """``(f_p, R_max)``: the retained fraction where ``p_b = p_p``, and the largest recovery."""
    f_p = theta * phi / (1.0 + phi - theta)
    return f_p, 1.0 - f_p


def _quad(fun, a, b):
    val, _ = quad(fun, a, b, epsabs=0.0, epsrel=1e-12, limit=200)
    return val


# codehealth: allow params -- the groups of E11 and the target recovery
def module_number(r_star, theta, phi, beta, n=0.5):
    """Membrane number ``N`` for recovery ``r_star`` (E12), by quadrature."""
    return _quad(lambda f: 1.0 / local_flux(pi_of_f(f, phi), theta, beta, n=n), 1.0 - r_star, 1.0)


# codehealth: allow params -- the groups of E11, the target recovery and the weight of E14
def area_elasticity(r_star, theta, phi, beta, weight="elasticity"):
    """``e_A = -dln a/dln beta`` at fixed recovery (E14): the area-weighted mean of ``eps_Q``.

    ``weight="share"`` averages the membrane's share of the drop instead of ``eps_Q``: not
    ``e_A``, kept only as a defect for the page's break table.
    """
    def w(f):
        p = pi_of_f(f, phi)
        if weight == "share":
            return (1.0 - film_share(p, theta, beta)) / local_flux(p, theta, beta)
        return elasticities(p, theta, beta)[0] / local_flux(p, theta, beta)
    return (_quad(w, 1.0 - r_star, 1.0)
            / _quad(lambda f: 1.0 / local_flux(pi_of_f(f, phi), theta, beta), 1.0 - r_star, 1.0))


def n_membrane_vacuum(r_star, phi):
    """Closed-form ``N`` of the membrane-limited module under vacuum permeate (E13)."""
    def F(f):
        return math.sqrt(f * (f + phi)) + phi * math.log(math.sqrt(f) + math.sqrt(f + phi))
    return (F(1.0) - F(1.0 - r_star)) / math.sqrt(1.0 + phi)


def a_film(r_star, theta, phi):
    """Closed-form film units ``a = k_p p_in A/F_H2,in`` of the film-limited module (E13)."""
    f1 = 1.0 - r_star
    c1, c0 = 1.0 + phi - theta, theta * phi
    return (1.0 - f1) / c1 + ((phi + c0 / c1) / c1) * math.log((c1 - c0) / (c1 * f1 - c0))


def recovery_at(N, theta, phi, beta):
    """Recovery of a module with membrane number ``N`` (E12 inverted; one root)."""
    r_max = pinch(theta, phi)[1]
    # N(R) diverges logarithmically at the pinch; 1e-6 short of it, N is already about
    # 11 at TC, and closer the quadrature loses its tolerance
    hit = find_roots(lambda r: module_number(r, theta, phi, beta) - N, 1e-9, r_max - 1e-6, n_scan=16)
    if len(hit["roots"]) != 1:
        raise ValueError(f"recovery_at: {hit['message']}")
    return hit["roots"][0]


# ---------------------------------------------------------------- module, pymrm
class MembraneModule:
    """Retentate ``f`` and permeate ``g`` on ``n_z`` cells, one pymrm finite-volume system.

    The state is ``(n_z, 2)``: domains ``"ret"`` and ``"perm"`` (``(n_z, 1)`` each) side by
    side on the field axis (``pymrm.BlockLayout``). Both streams flow along ``+xi`` with unit
    velocity (flow-rate form), first-order upwind, so each cell balance is backward Euler in
    ``xi``. The permeation source ``N j(f)`` leaves the retentate and enters the permeate in
    the same cell; its exact derivative enters through ``construct_coefficient_matrix``.

    Parameters
    ----------
    theta, phi, beta : float
        Permeate ratio, inert per H2 in the feed, conductance ratio.
    n_z : int
        Number of cells.
    source : {"cell", "trapezoid"}, optional
        ``"cell"``: the flux at the cell value (the scheme of the page). ``"trapezoid"``: the
        mean of the flux at the two faces of the cell, a second-order variant kept only as a
        defect for the order check.
    """

    # codehealth: allow params -- the three groups of E11, the grid and the source variant
    def __init__(self, theta, phi, beta, n_z=200, source="cell"):
        if source not in ("cell", "trapezoid"):
            raise ValueError(f"source {source!r} not recognised")
        self.theta, self.phi, self.beta, self.source = float(theta), float(phi), float(beta), source
        self.n_z = int(n_z)
        self.x_f = np.linspace(0.0, 1.0, self.n_z + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.layout = BlockLayout({"ret": (self.n_z, 1), "perm": (self.n_z, 1)}, axis=-1)
        self.shape = self.layout.shape
        self._build_operators()

    def _build_operators(self):
        # Outward normal, a dw/dn + b w = d.
        #   xi = 0 (inlet):  f = 1, g = 0        -> a = 0, b = 1, d = 1 (ret) and d = 0 (perm)
        #   xi = 1 (outlet): pure outflow, nothing imposed -> {"outflow": True}
        self.bc_ret = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"outflow": True})
        self.bc_perm = ({"a": 0.0, "b": 1.0, "d": 0.0}, {"outflow": True})
        sh = (self.n_z, 1)
        div = construct_div(sh, self.x_f, nu=0, axis=0)           # nu = 0: xi is Cartesian
        conv_r, conv_r_bc = construct_convflux_upwind(sh, self.x_f, self.x_c, self.bc_ret, v=1.0)
        conv_p, conv_p_bc = construct_convflux_upwind(sh, self.x_f, self.x_c, self.bc_perm, v=1.0)
        self.face_ret = (conv_r, conv_r_bc.toarray())               # face values of f (v = 1)
        place = self.layout.place
        self.jac_const = place(div @ conv_r, "ret", "ret") + place(div @ conv_p, "perm", "perm")
        g = np.zeros(self.shape)
        self.layout.put(g, "ret", (div @ conv_r_bc).toarray())
        self.layout.put(g, "perm", (div @ conv_p_bc).toarray())
        self.g_const = g.reshape((-1, 1))
        self.div = div

    # local flux and its derivatives at retentate values f
    def _flux(self, f):
        p = pi_of_f(f, self.phi)
        u = wall_root(p, self.theta, self.beta)
        j = (p - self.theta) / (self.beta + u + math.sqrt(self.theta))
        dj_dpi = 1.0 / (2.0 * u + self.beta)
        dpi_df = (1.0 + self.phi) * self.phi / (f + self.phi) ** 2
        dj_dbeta = -j / (2.0 * u + self.beta)
        return j, dj_dpi * dpi_df, dj_dbeta

    def _source(self, f):
        """Cell source ``s = j`` per unit ``N`` with ``ds/df`` (matrix) and ``ds/dbeta``."""
        if self.source == "cell":
            j, dj, jb = self._flux(f)
            return j, construct_coefficient_matrix(dj), jb
        conv, conv_bc = self.face_ret
        ff = conv @ f.reshape((-1, 1)) + conv_bc                  # face values, n_z + 1
        j, dj, jb = self._flux(ff.ravel())
        avg = 0.5 * (j[:-1] + j[1:])
        dmat = construct_coefficient_matrix(dj) @ conv
        dmat = 0.5 * (dmat[:-1] + dmat[1:])
        return avg, dmat, 0.5 * (jb[:-1] + jb[1:])

    def _residual(self, w, N):
        f = self.layout.view(w, "ret").ravel()
        s, ds, sb = self._source(f)
        src = np.zeros(self.shape)
        self.layout.put(src, "ret", (N * s).reshape((-1, 1)))
        self.layout.put(src, "perm", (-N * s).reshape((-1, 1)))
        res = self.g_const + self.jac_const @ w.reshape((-1, 1)) + src.reshape((-1, 1))
        place = self.layout.place
        jac = self.jac_const + place(N * ds, "ret", "ret") - place(N * ds, "perm", "ret")
        col_n = np.zeros(self.shape)
        self.layout.put(col_n, "ret", s.reshape((-1, 1)))
        self.layout.put(col_n, "perm", -s.reshape((-1, 1)))
        col_b = np.zeros(self.shape)
        self.layout.put(col_b, "ret", (N * sb).reshape((-1, 1)))
        self.layout.put(col_b, "perm", (-N * sb).reshape((-1, 1)))
        return res, jac, col_n.reshape((-1, 1)), col_b.reshape((-1, 1))

    def _guess(self, r):
        w = np.zeros(self.shape)
        self.layout.put(w, "ret", (1.0 - r * self.x_c).reshape((-1, 1)))
        self.layout.put(w, "perm", (r * self.x_c).reshape((-1, 1)))
        return w

    def solve_recovery(self, N, tol=1e-12):
        """Steady state at membrane number ``N``; returns ``self`` with ``w`` set.

        Raises
        ------
        RuntimeError
            When Newton does not converge.
        """
        self.N = float(N)

        def fun(x):
            res, jac, _, _ = self._residual(x.reshape(self.shape), self.N)
            return res, jac
        out = newton(fun, self._guess(0.5), tol=tol, maxfev=60, damping="monotone",
                     bounds=(0.0, None))
        if not out.success:
            raise RuntimeError(f"MembraneModule N = {N:g}, n_z = {self.n_z}: {out.message}")
        self.w, self.result = out.x.reshape(self.shape), out
        return self

    def solve_area(self, r_star, N0=None, tol=1e-12):
        """Membrane number ``N`` for recovery ``r_star``: ``N`` is a bordering unknown whose
        row is ``g(xi = 1) = r_star`` (the outlet face of the outflow boundary is the last
        cell). Also returns the exact discrete ``dN/dbeta`` from the same Jacobian.

        Raises
        ------
        RuntimeError
            When Newton does not converge.
        """
        m = self.n_z * 2
        last = m - 1                      # flat index of the last permeate cell: (n_z, 2), fields last
        e_last = np.zeros((1, m))
        e_last[0, last] = 1.0
        e_last = sp.csc_array(e_last)

        def fun(x):
            w, N = x[:m].reshape(self.shape), x[m]
            res, jac, col_n, _ = self._residual(w, N)
            r_full = np.vstack([res, [[w.reshape(-1)[last] - r_star]]])
            j_full = sp.block_array([[jac, sp.csc_array(col_n)], [e_last, None]], format="csc")
            return r_full, j_full
        x0 = np.concatenate([self._guess(r_star).ravel(), [N0 if N0 else 2.0 * (1.0 + self.beta)]])
        out = newton(fun, x0, tol=tol, maxfev=60, damping="monotone", bounds=(0.0, None))
        if not out.success:
            raise RuntimeError(f"MembraneModule area at R = {r_star:g}, beta = {self.beta:g}, "
                               f"n_z = {self.n_z}: {out.message}")
        x = out.x.ravel()
        self.w, self.N, self.result = x[:m].reshape(self.shape), float(x[m]), out
        # sensitivity: J dx/dbeta = -dR/dbeta (the bordering row has no beta in it)
        res, jac, col_n, col_b = self._residual(self.w, self.N)
        j_full = sp.block_array([[jac, sp.csc_array(col_n)], [e_last, None]], format="csc")
        rhs = -np.vstack([col_b, [[0.0]]])
        dx = spla.spsolve(j_full, rhs)
        self.dN_dbeta = float(dx[-1])
        return self

    def area_elasticity(self):
        """Discrete ``e_A = 1 - (beta/N) dN/dbeta`` at the last ``solve_area``."""
        return 1.0 - self.beta / self.N * self.dN_dbeta

    def outlet(self):
        """Outlet face values ``(f, g)`` read by ``compute_boundary_values`` under the bc."""
        f = self.layout.view(self.w, "ret")
        g = self.layout.view(self.w, "perm")
        _, _, f1, _ = compute_boundary_values(f, self.x_f, self.x_c, self.bc_ret)
        _, _, g1, _ = compute_boundary_values(g, self.x_f, self.x_c, self.bc_perm)
        return float(np.ravel(f1)[0]), float(np.ravel(g1)[0])

    def inlet_permeate(self):
        """Permeate face value at the inlet, read by ``compute_boundary_values``: zero by the bc."""
        g = self.layout.view(self.w, "perm")
        g0, _, _, _ = compute_boundary_values(g, self.x_f, self.x_c, self.bc_perm)
        return float(np.ravel(g0)[0])

    def recovery(self):
        """Permeated H2 over fed H2, the permeate outlet face value."""
        return self.outlet()[1]

    def balance_error(self):
        """Largest ``|f + g - 1|`` over the cells (H2 conservation, an identity of the scheme)."""
        f = self.layout.view(self.w, "ret").ravel()
        g = self.layout.view(self.w, "perm").ravel()
        return float(np.max(np.abs(f + g - 1.0)))

    def cell_flux_forms(self):
        """The three forms of E5 at every cell value of ``f``, shape ``(3, n_z)``."""
        p = pi_of_f(self.layout.view(self.w, "ret").ravel(), self.phi)
        return np.array([local_flux(p, self.theta, self.beta, form)
                         for form in ("ratio", "membrane", "film")])
