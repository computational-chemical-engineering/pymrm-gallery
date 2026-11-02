"""Levenspiel's four deactivation classes in four batch-solids contactings.

Model code of gallery page B2.3 (Levenspiel 1972). Everything is
dimensionless: time ``theta = t/t_obs``, reactant ``c = C_A/C_A0``, poison
``p = C_P/C_P0``, activity ``a``, reaction group ``kappa`` and decay group
``lam``. The rate of decay is ``-da/dtheta = lam phi^n' a^d`` with the class
factor ``phi = c`` (parallel), ``1 - c`` (series), ``p`` (side-by-side) or 1
(independent).

The activity is carried in the paper's own linearising variable
``w = (a^(1-d) - 1)/(1-d)`` (``ln a`` at ``d = 1``), in which every class reads
``dw/dtheta = -lam phi^n'``.

Contents
--------
a_of_w, w_of_a
    The activity and its linearising variable.
run_batch, run_mixed, run_const_c
    The three lumped contactings (Fig. 1a, 1c at constant flow, 1c with
    ``C_A`` held), ODEs in ``theta``.
PlugExact
    Plug flow over a batch of solids, pseudo-steady, by quadrature: the exact
    first-order cascade in the weight coordinate.
respond
    The observable of one contacting.
Bed
    The pymrm finite-volume bed: van Leer TVD convection (``ConvFluxTVD``),
    each fluid profile by ``newton``, the activity by Heun's method; and the
    full transient with the fluid accumulation term.

Several methods take an argument whose default is the correct model and whose
alternative is a defect the page's break table injects.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.optimize import brentq
from scipy.sparse import block_array, csr_array, vstack
from scipy.sparse.linalg import splu

from pymrm import (ConvFluxTVD, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, interp_cntr_to_stagg,
                   interp_cntr_to_stagg_tvd, newton, upwind, vanleer)

__all__ = ["TH", "SolveFailed", "a_of_w", "w_of_a", "run_batch", "run_mixed",
           "run_const_c", "PlugExact", "PLUG", "respond", "Bed"]

TH = np.linspace(0.0, 1.0, 41)
"""Observation times ``theta = t/t_obs`` of every observable."""

TINY_D = 1e-8


class SolveFailed(RuntimeError):
    """A fluid solve that did not converge (newton or the deferred sweeps)."""


def _pow(x, e):
    """``x**e``, finite at ``x = 0`` for a fractional ``e``."""
    return np.exp(e * np.log(np.clip(x, 1e-300, None)))


def _cascade(I, g, order):
    """Exact solution of ``dy/ds = -g y^order``, ``y(0) = 1``.

    Parameters
    ----------
    I : numpy.ndarray
        ``int_0^s`` of the coefficient profile.
    g : float
        Rate group.
    order : float
        Reaction order.

    Returns
    -------
    numpy.ndarray
        ``y`` at the points of ``I``.
    """
    if g == 0.0:
        return np.ones_like(I)
    if order == 1.0:
        return np.exp(-g * I)
    return _pow(np.clip(1.0 - (1.0 - order) * g * I, 0.0, None), 1.0 / (1.0 - order))


# -------------------------------------------- the paper's linearising variable
def a_of_w(w, d):
    """Activity from ``w = (a^(1-d) - 1)/(1-d)``.

    Parameters
    ----------
    w : array_like
        Linearising variable.
    d : float
        Deactivation order.

    Returns
    -------
    numpy.ndarray
        ``a``; 0 once the activity is used up (possible only for ``d < 1``).
    """
    w = np.asarray(w, float)
    if abs(1.0 - d) < TINY_D:
        return np.exp(w)
    x = (1.0 - d) * w
    ok = x > -1.0
    return np.where(ok, np.exp(np.log1p(np.where(ok, x, 0.0)) / (1.0 - d)), 0.0)


def w_of_a(a, d):
    """Linearising variable ``w`` of an activity ``a`` (inverse of :func:`a_of_w`)."""
    a = np.asarray(a, float)
    if abs(1.0 - d) < TINY_D:
        return np.log(a)
    return (_pow(a, 1.0 - d) - 1.0) / (1.0 - d)


def _extinction_event(d, idx=0):
    """Terminal event at the finite extinction time of ``d < 1`` (``w = -1/(1-d)``)."""
    def ev(_, y):
        return 1.0 + (1.0 - d) * y[idx] if d < 1.0 - TINY_D else 1.0
    ev.terminal, ev.direction = True, -1
    return ev


# ------------------------------------------------------------- Fig. 1a, batch
def run_batch(cls, kappa, lam, d, nprime, n=1.0, sig=0.0, theta=TH, rtol=1e-9):
    """Batch of solids and batch of fluid.

    Parameters
    ----------
    cls : {"parallel", "series", "side-by-side", "independent"}
        Deactivation class.
    kappa, lam : float
        Reaction and decay groups.
    d, nprime, n : float
        Deactivation order, concentration order of the decay, reaction order.
    sig : float, optional
        Poison consumption group (side-by-side).
    theta : numpy.ndarray, optional
        Observation times.
    rtol : float, optional
        LSODA relative tolerance (``atol`` 1e-13).

    Returns
    -------
    X, a : numpy.ndarray
        Conversion and activity at ``theta``.
    """
    def rhs(_, y):
        c, a, p = max(y[0], 1e-300), float(a_of_w(y[1], d)), max(y[2], 1e-300)
        f = {"parallel": c, "series": 1.0 - c, "side-by-side": p}.get(cls, 1.0)
        dw = -lam * _pow(max(f, 0.0), nprime)
        return [-kappa * _pow(c, n) * a, dw, sig * dw * _pow(a, d)]
    s = solve_ivp(rhs, (0.0, theta[-1]), [1.0, 0.0, 1.0], t_eval=theta,
                  rtol=rtol, atol=1e-13, method="LSODA",
                  events=_extinction_event(d, idx=1))
    c, a = np.empty(theta.size), np.zeros(theta.size)   # dead catalyst: a = 0,
    k = s.y[0].size                                     # and c frozen
    c[:k], a[:k] = s.y[0], a_of_w(s.y[1], d)
    c[k:] = s.y[0][-1] if k else 1.0
    return 1.0 - np.clip(c, 0.0, 1.0), a


# ------------------------------------------------ Fig. 1c, constant flow rate
def _mixed_c(kappa, a, n):
    """Root of ``1 - c = kappa c^n a``, the mixed-flow performance relation (eq. 16)."""
    if a <= 0.0:
        return 1.0
    if n == 1.0:
        return 1.0 / (1.0 + kappa * a)
    return brentq(lambda c: 1.0 - c - kappa * _pow(c, n) * a, 1e-14, 1.0,
                  xtol=1e-15, rtol=8.9e-16)


def _mixed_p(g, nprime):
    """Root of ``1 - p = g p^n'``, the poison balance over a well-mixed reactor."""
    if g <= 0.0:
        return 1.0
    if nprime == 1.0:
        return 1.0 / (1.0 + g)
    if nprime == 0.0:
        return max(1.0 - g, 0.0)
    return brentq(lambda p: 1.0 - p - g * _pow(p, nprime), 1e-14, 1.0,
                  xtol=1e-15, rtol=8.9e-16)


def run_mixed(cls, kappa, lam, d, nprime, n=1.0, sig=0.0, theta=TH, rtol=1e-9):
    """Batch of solids, mixed flow of fluid at constant flow rate.

    Parameters and returns as :func:`run_batch`.
    """
    def rhs(_, y):
        a = float(a_of_w(y[0], d))
        c = _mixed_c(kappa, a, n)
        f = (_mixed_p(sig * lam * _pow(a, d), nprime) if cls == "side-by-side"
             else {"parallel": c, "series": 1.0 - c}.get(cls, 1.0))
        return [-lam * _pow(max(f, 0.0), nprime)]
    s = solve_ivp(rhs, (0.0, theta[-1]), [0.0], t_eval=theta, rtol=rtol,
                  atol=1e-13, method="LSODA", events=_extinction_event(d))
    a = np.zeros(theta.size)                            # dead catalyst: a = 0
    a[:s.y[0].size] = a_of_w(s.y[0], d)
    return 1.0 - np.array([_mixed_c(kappa, ai, n) for ai in a]), a


# ------------------------------------------------- Fig. 1c, C_A held constant
def run_const_c(cls, kappa, lam, d, nprime, c_set, n=1.0, sig=0.0, theta=TH,
                rtol=1e-11):
    """Mixed flow with the flow rate lowered so that ``C_A`` stays at ``c_set``.

    The contacting the paper recommends.

    Parameters
    ----------
    c_set : float
        Held reactant level ``C_A/C_A0``.
    cls, kappa, lam, d, nprime, n, sig, theta, rtol
        As :func:`run_batch` (``atol`` 1e-14).

    Returns
    -------
    tau, a : numpy.ndarray
        Space time relative to its value on fresh catalyst, and activity.
    """
    r_set = 1.0 - c_set

    def tau_rel(a):
        return r_set / (kappa * _pow(c_set, n) * max(a, 1e-300))

    def rhs(_, y):
        a = float(a_of_w(y[0], d))
        f = (_mixed_p(sig * lam * _pow(a, d) * tau_rel(a), nprime)
             if cls == "side-by-side"
             else {"parallel": c_set, "series": r_set}.get(cls, 1.0))
        return [-lam * _pow(max(f, 0.0), nprime)]
    s = solve_ivp(rhs, (0.0, theta[-1]), [0.0], t_eval=theta, rtol=rtol,
                  atol=1e-14, method="LSODA", events=_extinction_event(d))
    a = np.zeros(theta.size)          # past extinction the flow rate needed to
    a[:s.y[0].size] = a_of_w(s.y[0], d)   # hold C_A is zero, i.e. tau' -> inf
    return np.array([tau_rel(ai) for ai in a]), a


# ----------------------------------------------- Fig. 1b, plug flow, exactly
class PlugExact:
    """Plug flow of fluid over a batch of solids, pseudo-steady, by quadrature.

    At frozen activity the fluid is a first-order cascade in the weight
    coordinate ``s``, ``c = [1 - (1 - n) kappa int_0^s a]^(1/(1-n))``, so the
    profile is a cumulative Simpson quadrature: no operator, no grid stencil
    and no linear solve in common with :class:`Bed`. Cumulative Simpson is
    linear in its integrand, so its matrix is built once.

    Parameters
    ----------
    ns : int, optional
        Odd number of nodes on ``s`` in [0, 1].
    """

    def __init__(self, ns=21):
        assert ns % 2 == 1, "Simpson needs an odd node count"
        self.s = np.linspace(0.0, 1.0, ns)
        self.ns = ns
        self.Q = np.column_stack([cumulative_simpson(np.eye(ns)[i], x=self.s, initial=0.0)
                                  for i in range(ns)])

    def profiles(self, a, kappa, lam, n, nprime, d, sig, cls):
        """Reactant and poison profiles at the nodes, for activity ``a``."""
        c = _cascade(self.Q @ a, kappa, n)
        p = (_cascade(self.Q @ _pow(a, d), sig * lam, nprime)
             if (cls == "side-by-side" and sig > 0.0) else np.ones_like(a))
        return c, p

    def run(self, cls, kappa, lam, d, nprime, n=1.0, sig=0.0, theta=TH, rtol=1e-8):
        """Conversion history and activity at the nodes.

        Parameters as :func:`run_batch` (LSODA, ``atol`` 1e-12).

        Returns
        -------
        X : numpy.ndarray
            Outlet conversion at ``theta``.
        A : numpy.ndarray
            Activity, nodes by times.
        """
        def rhs(_, w):
            a = a_of_w(w, d)
            c, p = self.profiles(a, kappa, lam, n, nprime, d, sig, cls)
            f = {"parallel": c, "series": 1.0 - c,
                 "side-by-side": p}.get(cls, np.ones_like(a))
            return -lam * _pow(np.clip(f, 0.0, None), nprime)
        sol = solve_ivp(rhs, (0.0, theta[-1]), np.zeros(self.ns), t_eval=theta,
                        rtol=rtol, atol=1e-12, method="LSODA")
        A = a_of_w(sol.y, d)
        X = np.array([1.0 - self.profiles(A[:, j], kappa, lam, n, nprime, d, sig, cls)[0][-1]
                      for j in range(theta.size)])
        return X, A


PLUG = PlugExact(21)
"""The quadrature bed of the separability sweep (21 nodes)."""


def respond(config, cls, kappa, lam, d, nprime, sig, n=1.0, c_set=0.4, theta=TH):
    """The observable of one contacting.

    Parameters
    ----------
    config : {"batch", "plug", "mixed", "constC"}
        Contacting: batch fluid, plug flow (by :data:`PLUG`), mixed flow at
        constant rate, mixed flow with ``C_A`` held at ``c_set``.
    cls, kappa, lam, d, nprime, sig, n, theta
        As :func:`run_batch`.
    c_set : float, optional
        Held level of the ``constC`` contacting.

    Returns
    -------
    numpy.ndarray
        Conversion at ``theta``, or ``ln tau'`` for ``constC``.
    """
    if config == "batch":
        return run_batch(cls, kappa, lam, d, nprime, n=n, sig=sig, theta=theta)[0]
    if config == "mixed":
        return run_mixed(cls, kappa, lam, d, nprime, n=n, sig=sig, theta=theta)[0]
    if config == "plug":
        return PLUG.run(cls, kappa, lam, d, nprime, n=n, sig=sig, theta=theta)[0]
    if config == "constC":
        return np.log(run_const_c(cls, kappa, lam, d, nprime, c_set, n=n,
                                  sig=sig, theta=theta)[0])
    raise ValueError(config)


# ------------------------------------------------------------ the pymrm bed
class Bed:
    """Plug flow over a batch of solids on a finite-volume grid.

    ``eps dc/dtheta + dc/ds = -kappa c^n a`` on ``s = w/W`` in [0, 1],
    ``c(0) = 1``; the activity obeys ``dw/dtheta = -lam phi^n'`` cell by cell.
    The pseudo-steady bed (``eps = 0``, :meth:`march`) solves the fluid at
    every instant; :meth:`transient` keeps the accumulation term.

    Parameters
    ----------
    m : int, optional
        Number of cells.
    limiter : callable or None, optional
        TVD limiter; None is first-order upwind.
    outlet : {"extrapolated", "outflow", "zero_gradient"}, optional
        Value carried by the outlet face: the linear extrapolation of the last
        two cells (default); the last cell value (``{"outflow": True}``); or
        the zero-gradient reconstruction of ``{"a": 1, "b": 0, "d": 0}``.
    route : {"newton", "deferred"}, optional
        ``"newton"``: ``ConvFluxTVD`` fluxes and Jacobian, ``newton``.
        ``"deferred"``: the same discrete equations by van Leer deferred
        correction, swept to a fixed point (a validation route).
    nu : int, optional
        Geometry of ``construct_div``; 0 (Cartesian) is the model.
    feed : bool, optional
        True imposes ``c = 1`` at the inlet; False a zero gradient there.
    tol : float, optional
        Newton step tolerance, or the largest change of a deferred sweep.
    """

    def __init__(self, m=200, limiter=vanleer, outlet="extrapolated", route="newton",
                 nu=0, feed=True, tol=1e-12):
        self.m, self.limiter, self.outlet_kind, self.route = int(m), limiter, outlet, route
        self.tol = float(tol)
        self.shape = (self.m, 1)
        self.s_f = np.linspace(0.0, 1.0, self.m + 1)
        self.s_c = 0.5 * (self.s_f[:-1] + self.s_f[1:])
        # bc: a dc/dn + b c = d with n the OUTWARD normal
        if feed:
            # inlet s = 0: c = 1 (the feed)                               -> a = 0, b = 1, d = 1
            bc_in = {"a": 0.0, "b": 1.0, "d": 1.0}
        else:
            # inlet s = 0: dc/ds = 0, no feed imposed                     -> a = 1, b = 0, d = 0
            bc_in = {"a": 1.0, "b": 0.0, "d": 0.0}
        if outlet == "zero_gradient":
            # outlet s = 1: dc/ds = 0, dc/dn = dc/ds                      -> a = 1, b = 0, d = 0
            bc_out = {"a": 1.0, "b": 0.0, "d": 0.0}
        else:
            # outlet s = 1: pure outflow, the face carries the last cell value
            bc_out = {"outflow": True}
        self.bc = (bc_in, bc_out)
        self.div = construct_div(self.shape, self.s_f, nu=nu)   # nu = 0: s is a weight fraction
        self.tvd = ConvFluxTVD(self.shape, self.s_f, self.s_c, self.bc,
                               tvd_limiter=upwind if limiter is None else limiter)
        # The outlet face value. With {"outflow": True} the face carries the last
        # cell value, half a cell upstream of the face: an O(h) error where c has
        # a gradient there. "extrapolated" adds the constant correction from that
        # value to the linear extrapolation of the last two cells (the last row
        # of interp_cntr_to_stagg without a bc), assembled once.
        K = csr_array((self.m + 1, self.m))
        if outlet == "extrapolated":
            w = np.zeros(self.m)
            w[-2:] = interp_cntr_to_stagg(np.eye(self.m, 2, k=2 - self.m), self.s_f,
                                          self.s_c, axis=0)[-1]
            w[-1] -= 1.0
            K = vstack([csr_array((self.m, self.m)), csr_array(w[np.newaxis, :])], format="csr")
        self.divK = (self.div @ K).tocsc()
        if route == "deferred":
            conv, conv_bc = construct_convflux_upwind(self.shape, self.s_f, self.s_c,
                                                      self.bc, v=1.0)
            self.L_up = (self.div @ conv).tocsc() + self.divK
            self.b_up = (self.div @ conv_bc).toarray().ravel()
        self.last_nit = 0

    # -- the fluid at frozen activity -------------------------------------------
    def residual(self, y, coef, order, mass=0.0, y_old=None):
        """Residual and Jacobian of ``mass (y - y_old) + d(F(y))/ds + coef y^order``.

        ``F`` are the van Leer face fluxes of ``ConvFluxTVD`` (velocity 1) with
        the outlet face set by ``outlet``.
        """
        y = np.asarray(y, float).reshape(self.shape)
        flux, jac = self.tvd(y, 1.0)
        yr = y.ravel()
        g = self.div @ np.asarray(flux).ravel() + self.divK @ yr + coef * _pow(yr, order)
        diag = coef * order * _pow(yr, order - 1.0) + mass
        if mass:
            g = g + mass * (yr - y_old)
        J = self.div @ jac + self.divK + construct_coefficient_matrix(
            np.broadcast_to(diag, (self.m,)).reshape(self.shape))
        return g, J.tocsc()

    def cascade(self, coef, order, y0=None):
        """Solve ``dy/ds = -coef(s) y^order``, ``y(0) = 1``, on the cells.

        Parameters
        ----------
        coef : numpy.ndarray or float
            Rate coefficient on the cells.
        order : float
            Reaction order.
        y0 : numpy.ndarray, optional
            Initial guess (default: the exact cascade of the cell coefficients).

        Returns
        -------
        numpy.ndarray
            ``y`` on the cells.

        Raises
        ------
        SolveFailed
            When ``newton`` reports no success, or the sweeps do not settle.
        """
        coef = np.broadcast_to(np.asarray(coef, float), (self.m,))
        if y0 is None:
            h = np.diff(self.s_f)
            y0 = _cascade(np.cumsum(coef * h) - 0.5 * coef * h, 1.0, order)
        if self.route == "deferred":
            return self._deferred(coef, order, y0)
        sol = newton(lambda y: self.residual(y, coef, order), np.reshape(y0, self.shape),
                     tol=self.tol, damping="monotone", stagnation=True)
        if not sol.success:
            raise SolveFailed(f"fluid profile, m = {self.m}: {sol.message}")
        self.last_nit = sol.nit
        return np.asarray(sol.x).ravel()

    def _deferred(self, coef, order, y0, max_it=300):
        """Upwind operator and reaction implicit (Newton in ``y^order``), van Leer
        correction of the previous sweep explicit."""
        y = np.asarray(y0, float).ravel()
        for sweep in range(1, max_it + 1):
            _, dface = interp_cntr_to_stagg_tvd(y.reshape(self.shape), self.s_f, self.s_c,
                                                self.bc, 1.0, tvd_limiter=self.limiter, axis=0)
            dcorr = self.div @ np.asarray(dface).ravel()
            diag = coef * order * _pow(y, order - 1.0)
            J = (self.L_up + construct_coefficient_matrix(diag.reshape(self.shape))).tocsc()
            r = self.L_up @ y + self.b_up + coef * _pow(y, order) + dcorr
            y_new = y + splu(J).solve(-r)
            change = np.max(np.abs(y_new - y))
            y = y_new
            if change < self.tol:
                self.last_nit = sweep
                return y
        raise SolveFailed(f"deferred correction: change {change:.1e} after {max_it} sweeps")

    def steady_c(self, kappa, a, n=1.0, c0=None):
        """Reactant profile at frozen activity ``a``."""
        return self.cascade(kappa * np.asarray(a, float), n, c0)

    def steady_p(self, g, a, d, nprime, p0=None):
        """Poison profile at frozen activity ``a`` (``g = sig lam``)."""
        if g == 0.0:
            return np.ones(self.m)
        return self.cascade(g * _pow(a, d), nprime, p0)

    def outlet(self, c):
        """The outlet concentration the operators transported (the outlet face value)."""
        c = np.asarray(c, float).reshape(self.shape)
        if self.outlet_kind == "extrapolated":
            return float(interp_cntr_to_stagg(c, self.s_f, self.s_c, axis=0)[-1, 0])
        return float(np.ravel(compute_boundary_values(c, self.s_f, self.s_c, self.bc)[2])[0])

    @staticmethod
    def phi(cls, c, p):
        """Class factor of the decay rate."""
        return {"parallel": c, "series": 1.0 - c, "side-by-side": p}.get(cls, np.ones_like(c))

    # -- pseudo-steady bed ------------------------------------------------------
    def march(self, cls, kappa, lam, d, nprime, n=1.0, sig=0.0, nt=200, theta=TH,
              scheme="heun", on_grid=True):
        """Pseudo-steady bed (``eps = 0``): fluid solved at each instant, ``w`` marched.

        Parameters
        ----------
        cls, kappa, lam, d, nprime, n, sig
            As :func:`run_batch`.
        nt : int, optional
            Number of equal steps to ``theta[-1]``.
        theta : numpy.ndarray, optional
            Times the conversion is returned at; each must be a step time.
        scheme : {"heun", "euler"}, optional
            Heun's method (second order) or forward Euler.
        on_grid : bool, optional
            True asserts that every time in ``theta`` is a step time. False
            interpolates the conversion linearly between steps (a defect the
            page injects: an O(dt^2) read-off error of its own).

        Returns
        -------
        X : numpy.ndarray
            Outlet conversion at ``theta``.
        A : numpy.ndarray
            Activity on the cells, cells by steps (``nt + 1`` columns).
        ts : numpy.ndarray
            Step times.
        """
        dt = float(theta[-1]) / nt
        if on_grid:
            k = np.asarray(theta, float) / dt
            assert np.all(np.abs(k - np.round(k)) < 1e-9), "observation time off the step grid"
        guess = {}

        def rate(w):
            a = a_of_w(w, d)
            c = self.steady_c(kappa, a, n, guess.get("c"))
            guess["c"] = c
            if cls == "side-by-side":
                p = self.steady_p(sig * lam, a, d, nprime, guess.get("p"))
                guess["p"] = p
            else:
                p = np.ones(self.m)
            f = np.clip(self.phi(cls, c, p), 0.0, None)
            return -lam * _pow(f, nprime), c, a

        w = np.zeros(self.m)
        f1, c, a = rate(w)
        ts, Xs, Aall = [0.0], [1.0 - self.outlet(c)], [a.copy()]
        for j in range(nt):
            if scheme == "euler":
                w = w + dt * f1
            else:
                f2, _, _ = rate(w + dt * f1)
                w = w + 0.5 * dt * (f1 + f2)
            f1, c, a = rate(w)
            ts.append((j + 1) * dt)
            Xs.append(1.0 - self.outlet(c))
            Aall.append(a.copy())
        return np.interp(theta, ts, Xs), np.array(Aall).T, np.array(ts)

    # -- full transient ---------------------------------------------------------
    def transient(self, cls, kappa, lam, d, nprime, eps, n=1.0, theta=TH, rtol=1e-8):
        """The bed with the fluid accumulation term ``eps dc/dtheta`` retained.

        Method of lines on the state ``[c, w]`` (cells first), BDF with the
        exact Jacobian: ``ConvFluxTVD`` for the convection, coefficient
        matrices for the pointwise couplings. Classes parallel, series and
        independent.

        Parameters
        ----------
        cls, kappa, lam, d, nprime, n, theta
            As :func:`run_batch`.
        eps : float
            Fluid residence time over the observation window.
        rtol : float, optional
            BDF relative tolerance (``atol`` 1e-11).

        Returns
        -------
        X, A : numpy.ndarray
            Outlet conversion at ``theta``; activity, cells by times.
        """
        if cls == "side-by-side":
            raise ValueError("transient: side-by-side is not implemented")
        dphi = {"parallel": 1.0, "series": -1.0}.get(cls, 0.0)
        m = self.m

        def ccm(v):
            return construct_coefficient_matrix(np.asarray(v, float).reshape(self.shape))

        def parts(y):
            c, w = y[:m], y[m:]
            a = a_of_w(w, d)
            f = np.clip(self.phi(cls, c, np.ones(m)), 0.0, None)
            return c, a, f

        def rhs(_, y):
            c, a, f = parts(y)
            flux, _ = self.tvd(c.reshape(self.shape), 1.0)
            rc = -(self.div @ np.asarray(flux).ravel() + self.divK @ c
                   + kappa * a * _pow(c, n)) / eps
            return np.concatenate([rc, -lam * _pow(f, nprime)])

        def jac(_, y):
            c, a, f = parts(y)
            _, jconv = self.tvd(c.reshape(self.shape), 1.0)
            dcc = -(self.div @ jconv + self.divK
                    + ccm(kappa * a * n * _pow(c, n - 1.0))) / eps
            dcw = ccm(-kappa * _pow(c, n) * _pow(a, d) / eps)        # da/dw = a^d
            dwc = ccm(-lam * nprime * _pow(f, nprime - 1.0) * dphi * (f > 0))
            return block_array([[dcc, dcw], [dwc, None]], format="csc")

        c0 = self.steady_c(kappa, np.ones(m), n)
        sol = solve_ivp(rhs, (0.0, float(theta[-1])), np.concatenate([c0, np.zeros(m)]),
                        t_eval=theta, rtol=rtol, atol=1e-11, method="BDF", jac=jac)
        if not sol.success:
            raise SolveFailed(f"transient bed: {sol.message}")
        X = np.array([1.0 - self.outlet(sol.y[:m, j]) for j in range(theta.size)])
        return X, a_of_w(sol.y[m:, :], d)
