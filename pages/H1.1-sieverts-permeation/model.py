"""Sieverts-law hydrogen permeation (Itoh 1987, Eq. 1): the wall and the channel.

Model code of gallery page H1.1. Two separate single-field models, because they
answer different questions:

* the palladium WALL, steady radial diffusion of dissolved atomic hydrogen across
  a cylindrical shell with Sieverts equilibrium at both faces; its throughput at
  unit driving force is the permeance ``alpha_H`` of Eq. 1;
* a membrane CHANNEL, isothermal plug flow of hydrogen in an inert, hydrogen
  removed through the wall by an ``n``-power law against a permeate side at a
  fixed hydrogen mole fraction.

Contents
--------
sieverts_c
    Sieverts' law, the dissolved concentration in equilibrium with the gas.
wall_solve, alpha_from_wall
    The wall as a pymrm finite-volume solve, and the permeance read off it.
Permeator
    The channel as a pymrm finite-volume system solved by ``pymrm.newton`` in
    ``s = sqrt(w)``, and a cell-by-cell march of the same discrete equations as
    a second route.
G_half, exact_half, exact_quadrature, exact_lin
    Exact channel solutions that share no code with the finite-volume assembly.
q_series, n_apparent
    A surface step in series with the wall, and the apparent pressure exponent
    it produces.
root_log
    A root of an increasing function of a positive parameter, bracket widened
    until it brackets.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.optimize import brentq
from scipy.sparse import diags_array
from scipy.sparse.linalg import spsolve

from pymrm import (NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, construct_grad, newton)

__all__ = ["sieverts_c", "wall_solve", "alpha_from_wall", "Permeator", "G_half",
           "exact_half", "exact_quadrature", "exact_lin", "q_series", "n_apparent", "root_log"]


def sieverts_c(p_over_P0, C_0):
    """Sieverts' law, ``c = C_0 sqrt(p/P_0)``.

    Parameters
    ----------
    p_over_P0 : float or numpy.ndarray
        Hydrogen partial pressure over the reference pressure ``P_0``.
    C_0 : float
        Dissolved concentration in equilibrium with ``P_0`` (mol/m3).
    """
    return C_0 * np.sqrt(p_over_P0)


def wall_solve(c_in, c_out, r_i, r_o, D, l_0, n_r=200, nu=1, D_profile=None):
    """Steady radial diffusion of dissolved hydrogen through a cylindrical wall.

    ``(1/r) d/dr (r D dc/dr) = 0`` on ``r_i <= r <= r_o`` with the face values
    fixed by Sieverts equilibrium. The problem is linear, so one sparse solve
    finishes it. ``D`` is given at the faces: uniform, or ``D_profile(r_f)``
    evaluated at the face positions, so a smoothly graded wall needs no face
    averaging (a jump in ``D`` between cells would need the harmonic mean).

    Parameters
    ----------
    c_in, c_out : float
        Dissolved concentration at ``r_i`` and ``r_o`` (mol/m3).
    r_i, r_o : float
        Inner and outer radius (m).
    D : float
        Diffusivity of dissolved hydrogen (m2/s), used when ``D_profile`` is None.
    l_0 : float
        Tube length (m).
    n_r : int, optional
        Number of cells across the wall.
    nu : int, optional
        Geometry of ``construct_div``: 1 cylindrical (the wall). 0 and 2 exist
        only so that a defect can be injected on purpose.
    D_profile : callable, optional
        ``D_profile(r)`` (m2/s) for a wall whose diffusivity varies with radius.

    Returns
    -------
    r_c : numpy.ndarray
        Cell centres (m).
    c : numpy.ndarray
        Dissolved concentration in the cells (mol/m3).
    r_f : numpy.ndarray
        Face positions (m).
    Q_f : numpy.ndarray
        Molar rate through each cylindrical face, ``2 pi r l_0 (-D dc/dr)``
        (mol/s); constant in ``r`` for the cylindrical assembly.
    """
    shape = (n_r, 1)                           # one field: (n, 1), never (n,)
    r_f = np.linspace(r_i, r_o, n_r + 1)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    # outward normal, a dc/dn + b c = d; Sieverts equilibrium at both faces
    bc = ({"a": 0.0, "b": 1.0, "d": c_in},    # r = r_i: c = C_0 sqrt(p_H/P_0)
          {"a": 0.0, "b": 1.0, "d": c_out})   # r = r_o: c = C_0 sqrt(p'_H/P_0)
    grad_mat, grad_bc = construct_grad(shape, r_f, r_c, bc, axis=0)
    div_mat = construct_div(shape, r_f, nu=nu, axis=0)     # nu = 1: cylindrical
    D_f = np.full(n_r + 1, float(D)) if D_profile is None else np.asarray(D_profile(r_f), float)
    diff_mat = construct_coefficient_matrix(D_f.reshape((-1, 1)), shape=shape, axis=0)
    flux_mat = -diff_mat @ grad_mat
    flux_bc = -diff_mat @ grad_bc
    jac = (div_mat @ flux_mat).tocsc()
    g_const = (div_mat @ flux_bc).toarray().ravel()
    c = spsolve(jac, -g_const)
    J_f = (flux_mat @ c + flux_bc.toarray().ravel())        # -D dc/dr at every face
    Q_f = 2.0 * np.pi * r_f * l_0 * J_f
    return r_c, c, r_f, Q_f


def alpha_from_wall(r_i, r_o, D, C_0, l_0, n_r=200, nu=1):
    """Permeance ``alpha_H`` read off the wall solve (mol/s).

    The rate through the inner face at unit square-root driving force,
    ``p_H = P_0`` against ``p'_H = 0``. Arguments as in ``wall_solve``.
    """
    _, _, _, Q_f = wall_solve(sieverts_c(1.0, C_0), 0.0, r_i, r_o, D, l_0, n_r=n_r, nu=nu)
    return Q_f[0]


class Permeator:
    """Isothermal plug-flow membrane channel, hydrogen removed by an ``n``-power law.

    ``dw/dL = -Pi [ (w/(w + beta))^n - y_perm^n ]``, ``w(0) = 1``, on
    ``0 <= L <= length``, with ``w = u_H/u_H0`` the scaled hydrogen flow,
    ``beta`` the inert flow over ``u_H0`` and ``Pi = alpha_H sqrt(P_T/P_0)/u_H0``
    the capacity ratio. First-order upwind convection at unit velocity (the
    balance is in flow-rate form), assembled once; the pointwise sink is
    linearised by ``NumJac`` on the ``(n_z, 1)`` state. ``solve`` iterates on
    ``s = sqrt(w)``, in which the half-power sink is smooth through zero.

    Parameters
    ----------
    Pi : float
        Capacity ratio.
    beta : float
        Inert flow over the inlet hydrogen flow.
    n : float, optional
        Pressure exponent, 1/2 for Sieverts' law.
    y_perm : float, optional
        Hydrogen mole fraction on the permeate side.
    n_z : int, optional
        Number of cells.
    length : float, optional
        Dimensionless channel length.
    eps_jac : float, optional
        Perturbation of ``NumJac`` (the pymrm default, 1e-6, is an absolute
        step for unknowns below 1).
    sink_sign, use_mole_fraction, nu, outlet
        Defects for the page's break table, not modelling options:
        ``sink_sign=+1`` adds hydrogen, ``use_mole_fraction=False`` drives the
        flux with ``w`` itself, ``nu=1`` treats ``L`` as a radius, and
        ``outlet`` is ``"outflow"`` (the model), ``"zero_gradient"`` or
        ``"dirichlet"`` (``w = 0`` imposed at the exit).

    Attributes
    ----------
    w : numpy.ndarray
        The converged state, ``(n_z, 1)``, set by ``solve``.
    result : scipy.optimize.OptimizeResult
        The last ``pymrm.newton`` result.
    """

    #: lower bound 0 on the unknown (s, or w in the comparison form)
    BOUNDS = (0.0, None)

    def __init__(self, Pi, beta, n=0.5, y_perm=0.0, n_z=800, length=1.0, eps_jac=1e-6,
                 sink_sign=-1.0, use_mole_fraction=True, nu=0, outlet="outflow"):
        if outlet not in ("outflow", "zero_gradient", "dirichlet"):
            raise ValueError(f"outlet {outlet!r} not recognised")
        self.Pi, self.beta, self.n, self.y_perm = float(Pi), float(beta), float(n), float(y_perm)
        self.sink_sign, self.use_mole_fraction = float(sink_sign), use_mole_fraction
        self.nu, self.outlet, self.eps_jac = nu, outlet, eps_jac
        self.z_f = np.linspace(0.0, length, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.shape = (n_z, 1)                   # one field: (n, 1), never (n,)
        self._build_operators()
        self.w, self.result = None, None

    def _build_operators(self):
        # outward normal, a dw/dn + b w = d
        #   inlet  (L = 0):      w = 1                     -> a = 0, b = 1, d = 1
        #   outlet (L = length): nothing imposed, convection only -> {"outflow": True}
        #   (comparisons: dw/dL = 0, a = 1, b = 0; or w = 0, a = 0, b = 1, d = 0)
        exit_bc = {"outflow": {"outflow": True},
                   "zero_gradient": {"a": 1.0, "b": 0.0, "d": 0.0},
                   "dirichlet": {"a": 0.0, "b": 1.0, "d": 0.0}}[self.outlet]
        self.bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, exit_bc)
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c, self.bc,
                                                      v=1.0, axis=0)
        div_mat = construct_div(self.shape, self.z_f, nu=self.nu, axis=0)   # nu = 0: L is axial
        self.jac_const = (div_mat @ conv_mat).tocsc()
        self.g_const = (div_mat @ conv_bc).toarray().reshape((-1, 1))
        self.numjac = NumJac(self.shape, eps_jac=self.eps_jac)   # pointwise sink

    def driving(self, w):
        """``(w/(w + beta))^n - y_perm^n``, the driving force of the sink."""
        y = w / (w + self.beta) if self.use_mole_fraction else w
        return y ** self.n - self.y_perm ** self.n

    def source(self, w):
        """The sink ``-Pi * driving(w)`` (with ``sink_sign = -1``)."""
        return self.sink_sign * self.Pi * self.driving(w)

    def residual(self, w):
        """Steady residual ``dw/dL - source(w)`` and its Jacobian."""
        g_src, jac_src = self.numjac(self.source, w)
        g = self.g_const + self.jac_const @ w.reshape((-1, 1)) - g_src.reshape((-1, 1))
        return g, self.jac_const - jac_src

    def residual_s(self, s):
        """Steady residual in ``s = sqrt(w)`` and its Jacobian with respect to ``s``.

        The same discrete equations as ``residual``, with ``w = s^2`` substituted:
        ``jac_const @ s^2 - source(s^2)``. For ``n = 1/2`` the sink becomes
        ``-Pi s / sqrt(s^2 + beta)``, smooth and linear through ``s = 0``. The
        Jacobian is ``jac_const @ diag(2 s)`` minus ``NumJac`` of the composed sink.
        """
        g_src, jac_src = self.numjac(lambda x: self.source(x * x), s)
        w = s * s
        g = self.g_const + self.jac_const @ w.reshape((-1, 1)) - g_src.reshape((-1, 1))
        jac = self.jac_const @ diags_array(2.0 * s.ravel()) - jac_src
        return g, jac.tocsc()

    def flat_guess(self):
        """The inlet value in every cell."""
        return np.ones(self.shape)

    def march_guess(self):
        """The discrete equations solved one cell at a time, as an initial guess.

        First-order upwind on a uniform grid with the outflow outlet makes the
        system bidiagonal, so cell ``i`` couples only to cell ``i - 1``:
        ``(w_i - w_(i-1))/h + Pi [ (w_i/(w_i + beta))^n - y_perm^n ] = 0``.
        That residual is continuous and strictly increasing in ``w_i``, at most
        0 at the flux-free floor ``beta y_perm/(1 - y_perm)`` and positive at
        ``w_(i-1)``, so each cell has exactly one root in between, found with
        ``brentq``. Past the extinction point the root decays quadratically per
        cell and underflows to 0, which then satisfies the cell equation. The
        march solves the assembled cell equations by scalar root finding, an
        independent route to the same discrete solution as ``solve``.
        """
        h = np.diff(self.z_f)
        w_stop = self.beta * self.y_perm / (1.0 - self.y_perm) if self.y_perm > 0 else 0.0
        w, prev = np.empty(len(h)), 1.0
        for i in range(len(h)):
            def cell(x, prev=prev, h_i=h[i]):
                return (x - prev) / h_i + self.Pi * self.driving(x)
            if prev <= w_stop or cell(w_stop) >= 0.0:
                prev = w_stop                  # already at the floor
            else:
                prev = brentq(cell, w_stop, prev, xtol=1e-300, rtol=8.9e-16, maxiter=500)
            w[i] = prev
        return w.reshape(self.shape)

    def solve(self, w0=None, tol=1e-10, maxfev=100, damping="monotone", variable="s",
              bounds=BOUNDS, require=True):
        """Solve with ``pymrm.newton`` from ``w0`` (default ``flat_guess``).

        ``variable="s"`` (the model) solves for ``s = sqrt(w)`` with the lower
        bound ``s >= 0``; ``w = s^2`` is then non-negative by construction.
        ``variable="w"`` solves for ``w`` itself with the bound ``w >= 0``, kept
        as a comparison: in ``w`` every tail cell's Newton step crosses the bound,
        the bound rule then caps the step length, and the iteration stalls above
        ``Pi_crit``. ``tol`` is on the order-one unknown, above its round-off floor.

        Parameters
        ----------
        w0 : numpy.ndarray, optional
            Initial guess for ``w``, ``(n_z, 1)``.
        tol, maxfev, damping
            Passed to ``pymrm.newton``.
        variable : {"s", "w"}, optional
            The unknown Newton iterates on.
        bounds : tuple, optional
            Passed to ``pymrm.newton``; ``None`` removes the lower bound.
        require : bool, optional
            True: raise ``RuntimeError`` when Newton does not converge.
            False: return the result for the caller to inspect; ``w`` is then
            set only on success.

        Returns
        -------
        scipy.optimize.OptimizeResult
            The ``pymrm.newton`` result (``x`` is ``s`` or ``w``).
            ``NewtonBreakdown`` (a singular Jacobian, a non-finite step) propagates.
        """
        if variable not in ("s", "w"):
            raise ValueError(f"variable {variable!r} not recognised")
        w0 = self.flat_guess() if w0 is None else np.asarray(w0, dtype=float)
        if variable == "s":
            fun, x0 = self.residual_s, np.sqrt(w0.reshape(self.shape))
        else:
            fun, x0 = self.residual, w0.reshape(self.shape)
        # a damping trial beyond the bound evaluates a power of a negative number
        # (nan); newton rejects such a trial and shortens the step
        with np.errstate(invalid="ignore"):
            result = newton(fun, x0, tol=tol, maxfev=maxfev, damping=damping, bounds=bounds)
        self.result = result
        if result.success:
            x = result.x.reshape(self.shape)
            self.w = x * x if variable == "s" else x
        elif require:
            raise RuntimeError(f"Permeator at Pi = {self.Pi:g}, n_z = {self.shape[0]}: "
                               f"{result.message}")
        return result

    def outlet_value(self):
        """``w`` on the exit face, from ``compute_boundary_values``."""
        _, _, w_out, _ = compute_boundary_values(self.w, self.z_f, self.z_c, self.bc)
        return float(np.ravel(w_out)[0])


def G_half(w, beta):
    """Antiderivative for ``n = 1/2``, ``y_perm = 0``: ``dG/dw = sqrt((w + beta)/w)``.

    Then ``G(w0) - G(w(L)) = Pi L`` exactly, and ``G(0) = 0``, so the hydrogen
    is exhausted at the finite length ``L_ext = G(w0)/Pi``.
    """
    return np.sqrt(w * (w + beta)) + beta * np.arcsinh(np.sqrt(w / beta))


def exact_half(L, beta, Pi, w0=1.0):
    """``w(L)`` for ``n = 1/2``, ``y_perm = 0``, by inverting ``G``."""
    out = []
    for t in np.atleast_1d(G_half(w0, beta) - Pi * np.asarray(L, float)):
        out.append(0.0 if t <= 0.0 else
                   brentq(lambda w: G_half(w, beta) - t, 0.0, w0, xtol=1e-16, rtol=8.9e-16))
    return np.array(out)


def exact_quadrature(L, beta, Pi, n=0.5, y_perm=0.0, w0=1.0):
    """``w(L)`` for any ``n`` and permeate loading, by separation and quadrature.

    ``L(w) = (1/Pi) int_w^w0 dw' / [ (w'/(w' + beta))^n - y_perm^n ]``,
    inverted by ``brentq``. Shares nothing with the finite-volume assembly.
    """
    def f(x):
        return (x / (x + beta)) ** n - y_perm ** n
    w_stop = beta * y_perm / (1.0 - y_perm) if y_perm > 0 else 0.0   # flux-free floor

    def L_of(w):
        # the integrand is singular at w_stop, approached only asymptotically;
        # every w evaluated here lies above it, so quad's endpoint warning is noise
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            return quad(lambda x: 1.0 / f(x), w, w0, limit=400)[0] / Pi
    out = []
    for Li in np.atleast_1d(np.asarray(L, float)):
        lo = w_stop * (1.0 + 1e-12) + 1e-300
        if L_of(lo) <= Li:
            out.append(lo)
            continue
        out.append(brentq(lambda w: L_of(w) - Li, lo, w0, xtol=1e-15, rtol=8.9e-16))
    return np.array(out)


def exact_lin(L, beta, Pi, w0=1.0):
    """``w(L)`` for ``n = 1``, ``y_perm = 0``: ``w + beta ln w = w0 + beta ln w0 - Pi L``.

    ``w = 0`` needs ``L`` infinite: a linear pressure law never empties the tube.
    """
    rhs = w0 + beta * np.log(w0) - Pi * np.asarray(L, float)
    out = []
    for r in np.atleast_1d(rhs):
        out.append(brentq(lambda w: w + beta * np.log(w) - r, 1e-300, w0,
                          xtol=1e-300, rtol=8.9e-16))
    return np.array(out)


def q_series(pi1, lam, pi2=0.0, surface_power=2.0):
    """Flux through (surface | wall | surface) in series, scaled on ``alpha_H``.

    ``q = a (pi1 - theta1^2) = alpha_H (theta1 - theta2) = a (theta2^2 - pi2)``,
    ``lam = alpha_H/a``. Eliminating ``theta1`` and ``theta2`` leaves
    ``(q + sqrt(lam q + pi2))^2 = pi1 - lam q``, solved for ``q`` with
    ``brentq``; the bracket is doubled until it brackets. Reduces to
    ``sqrt(pi1) - sqrt(pi2)`` at ``lam -> 0`` and ``(pi1 - pi2)/(2 lam)`` at
    ``lam -> inf``. ``surface_power`` is the stoichiometric power of the
    surface step (2 for dissociation); another value is a break-table defect
    and is only implemented for ``pi2 = 0``.
    """
    if lam <= 0.0:
        return np.sqrt(pi1) - np.sqrt(pi2)
    if surface_power == 2.0:
        def f(q):
            return (q + np.sqrt(lam * q + pi2)) ** 2 - (pi1 - lam * q)
    else:
        if pi2 != 0.0:
            raise ValueError("surface_power != 2 is implemented for pi2 = 0 only")
        p = surface_power

        def f(q):
            theta2 = (lam * q) ** (1.0 / p)
            return (q + theta2) ** p - (pi1 - lam * q)
    hi = np.sqrt(pi1)
    while f(hi) < 0:
        hi *= 2.0
    return brentq(f, 1e-300, hi, xtol=1e-18, rtol=8.9e-16)


def n_apparent(pi1, lam, d=1e-4, surface_power=2.0):
    """Local exponent ``d ln q / d ln pi1`` at ``pi1``, permeate side at zero."""
    lo = np.log(q_series(pi1 * np.exp(-d), lam, surface_power=surface_power))
    hi = np.log(q_series(pi1 * np.exp(+d), lam, surface_power=surface_power))
    return (hi - lo) / (2 * d)


def root_log(f, lo, hi, grow=10.0, max_grow=20, xtol=1e-12):
    """Root of an increasing ``f`` in ``ln p`` for a positive parameter ``p``.

    ``lo`` and ``hi`` are first guesses: while ``f(lo) > 0`` ``lo`` is divided
    by ``grow``, while ``f(hi) < 0`` ``hi`` is multiplied by it, at most
    ``max_grow`` times each. Raises ``ValueError`` when no bracket is found.
    """
    for _ in range(max_grow):
        if f(lo) <= 0.0:
            break
        lo /= grow
    else:
        raise ValueError(f"root_log: f > 0 down to {lo:g}")
    for _ in range(max_grow):
        if f(hi) >= 0.0:
            break
        hi *= grow
    else:
        raise ValueError(f"root_log: f < 0 up to {hi:g}")
    return float(np.exp(brentq(lambda s: f(np.exp(s)), np.log(lo), np.log(hi), xtol=xtol)))
