"""Zwietering's four states of mixing for an Erlang residence-time distribution.

Model code of gallery page A2.8 (Zwietering 1959). The page imports these; the
three closed forms Zwietering printed for two vessels (his eqs. 38, 39 with 40,
and 42) stay in the notebook, so that the transcribed algebra and the model
never share a line of code.

Everything is dimensionless: time and life expectation in units of the mean
residence time ``tau`` (``x = t/tau`` or ``x = lambda/tau``), concentration in
units of the feed (``gamma = c/c0``), and the rate ``R = k c^p`` enters through
``K = k c0^(p-1) tau``.

Contents
--------
E_rtd, W_rtd
    Frequency function ``f`` and the fraction still inside ``w = 1 - F`` of
    ``n`` equal stirred tanks in series (Erlang distribution).
batch
    A batch element of order ``p``.
gamma_seg, gamma_real_seg, gamma_real_mixed
    Complete segregation (his eqs. 36 and 37), ``n`` real internally
    segregated vessels, ``n`` real ideally mixed vessels (his eq. 41).
gamma_far
    The far-field root ``K g^p + n (g - 1) = 0`` at infinite life expectation.
MaxMixedness
    Maximum mixedness, his eq. (31), in conservative form as a steady 1-D
    convection-reaction problem on the pymrm operators, van Leer TVD solved by
    Newton with the exact Jacobian (``ConvFluxTVD``, pitfalls P10), or by
    deferred correction as a second solver route.
mm_ode
    The same state from eq. (31) in its non-conservative form, integrated as
    an initial-value problem by DOP853: no grid, no flux, no pymrm.
var_alpha, J_real_mixed, J_real_seg, alpha_P, J_mm, var_within_mm
    Danckwerts' degree of segregation ``J = var(alpha_P)/var(alpha)`` for the
    four states (Zwietering's eqs. 11, 27 and 28 and Appendix I).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import gammaincc, gammaln

from pymrm import (ConvFluxTVD, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, interp_cntr_to_stagg_tvd,
                   newton, vanleer)

__all__ = ["E_rtd", "W_rtd", "batch", "gamma_seg", "gamma_real_seg", "gamma_real_mixed",
           "gamma_far", "MaxMixedness", "mm_ode", "var_alpha", "J_real_mixed",
           "J_real_seg", "alpha_P", "J_mm", "var_within_mm"]


# ---- residence-time distribution of n equal stirred tanks, mean 1 ----------------
def E_rtd(x, n):
    """Frequency function of ``n`` equal tanks in series.

    Parameters
    ----------
    x : array_like
        Age or life expectation over the mean residence time.
    n : int
        Number of tanks.

    Returns
    -------
    numpy.ndarray
        ``tau f(tau x)``, the Erlang density of shape ``n`` and scale ``1/n``,
        evaluated through logarithms so that it does not overflow at large ``n``.
    """
    x = np.asarray(x, float)
    if n == 1:
        return np.exp(-x)
    out = np.zeros_like(x)
    m = x > 0
    out[m] = np.exp(n * np.log(n) + (n - 1) * np.log(x[m]) - n * x[m] - gammaln(n))
    return out


def W_rtd(x, n):
    """Fraction of the feed still inside, ``w = 1 - F``.

    Parameters
    ----------
    x : array_like
        Age or life expectation over the mean residence time.
    n : int
        Number of tanks.

    Returns
    -------
    numpy.ndarray
        The regularised upper incomplete gamma function ``Q(n, n x)``.
    """
    return gammaincc(n, n * np.asarray(x, float))


def batch(x, K, p=2.0):
    """Batch element ``d gamma/d theta = -K gamma^p``, ``gamma(0) = 1``.

    Parameters
    ----------
    x : array_like
        Time over the mean residence time.
    K : float
        ``k c0^(p-1) tau``.
    p : float, optional
        Reaction order.

    Returns
    -------
    numpy.ndarray
        ``gamma(x)``; below order one the element is exhausted in finite time
        and stays at zero.
    """
    x = np.asarray(x, float)
    if p == 1.0:
        return np.exp(-K * x)
    if p > 1.0:
        return (1.0 + (p - 1.0) * K * x) ** (-1.0 / (p - 1.0))
    return np.clip(1.0 - (1.0 - p) * K * x, 0.0, None) ** (1.0 / (1.0 - p))


# ---- the three states that need no PDE -------------------------------------------
def gamma_seg(K, n, p=2.0):
    """Complete segregation, Zwietering eqs. (36) and (37).

    Parameters
    ----------
    K : float
        ``k c0^(p-1) tau``.
    n : int
        Number of tanks of the RTD.
    p : float, optional
        Reaction order.

    Returns
    -------
    float
        ``c_e/c0 = int f(t) gamma_batch(t) dt`` by adaptive quadrature, split at
        multiples of ``1/n`` and, below order one, at the exhaustion time.
    """
    f = lambda x: float(E_rtd(np.array([x]), n)[0] * batch(np.array([x]), K, p)[0])
    knots = [0.0, .5 / n, 1. / n, 2. / n, 4. / n, 10. / n, 30. / n]
    if p < 1.0:
        t_ex = 1.0 / ((1.0 - p) * K)
        if not np.any(np.isclose(knots, t_ex, rtol=1e-12, atol=0.0)):
            knots = sorted(knots + [t_ex])
    v = sum(quad(f, a, b, limit=400)[0] for a, b in zip(knots[:-1], knots[1:]))
    return v + quad(f, knots[-1], np.inf, limit=400)[0]


def gamma_real_seg(K, n):
    """``n`` real vessels in series, each internally segregated, second order.

    Batch elements inside each tank, the feed to tank ``i + 1`` being tank
    ``i``'s cup-mixed outflow: the chain Zwietering's eqs. (39) and (40) write
    out for ``n = 2``.

    Parameters
    ----------
    K : float
        ``k c0 tau``.
    n : int
        Number of tanks.

    Returns
    -------
    float
        ``c_e/c0``.
    """
    g = 1.0
    for _ in range(n):
        Kv = K * g
        g *= quad(lambda u: n * np.exp(-n * u) / (1.0 + Kv * u), 0, np.inf, limit=500)[0]
    return g


def gamma_real_mixed(K, n):
    """``n`` real ideally mixed vessels in series, second order (his eq. 41).

    Parameters
    ----------
    K : float
        ``k c0 tau``.
    n : int
        Number of tanks.

    Returns
    -------
    float
        ``c_e/c0`` from the chain ``K g_i^2 = n (g_(i-1) - g_i)``.
    """
    g = 1.0
    for _ in range(n):
        g = (-n + np.sqrt(n * n + 4 * K * n * g)) / (2 * K)
    return g


def gamma_far(K, n, p=2.0):
    """Far-field root ``K g^p + n (g - 1) = 0``, Zwietering's ``dc/dlambda = 0``.

    Parameters
    ----------
    K : float
        ``k c0^(p-1) tau``.
    n : int
        Number of tanks; ``f/(1 - F)`` tends to ``n`` at large life expectation.
    p : float, optional
        Reaction order.

    Returns
    -------
    float
        The root in ``(0, 1]``.
    """
    return brentq(lambda g: K * max(g, 0.0) ** p + n * (g - 1.0), 1e-15, 1.0)


# ---- maximum mixedness on the pymrm operators -------------------------------------
class MaxMixedness:
    """Zwietering eq. (31) in conservative form, a steady 1-D convection-reaction problem.

    ``div(v gamma) = E(x) - K w(x) gamma^p`` with ``v = -w(x)`` on the faces,
    ``x = lambda/tau`` the life expectation, ``gamma = c/c0`` and ``w = 1 - F``.
    The flow runs toward decreasing ``x`` and leaves the domain at the exit
    ``x = 0``; at ``x = X`` it enters with the far-field root. The operators are
    assembled once in the constructor; :meth:`solve` only solves.

    Parameters
    ----------
    K : float
        ``k c0^(p-1) tau``.
    n : int
        Number of tanks of the RTD.
    X : float, optional
        Truncated domain length in life expectation.
    n_cells : int, optional
        Number of cells, graded as ``x = X t^grade`` so that they pack at the exit.
    grade : float, optional
        Grid grading exponent; 1 gives a uniform grid.
    p : float, optional
        Reaction order. The rate is written ``gamma |gamma|^(p-1)`` so that a
        negative iterate meets a restoring rate (pitfalls P9).
    limiter : callable or None, optional
        TVD limiter of the convective face values; ``None`` is first-order upwind.
    rxn_sign : float, optional
        Sign of the reaction term; ``-1`` exists for the page's break table.
    feed : bool, optional
        Include the side feed ``E(x)``; ``False`` exists for the break table.
    v_const : bool, optional
        Face velocity ``-1`` instead of ``-w``; for the break table.
    far_field : {"root", "one"}, optional
        Inflow value at ``x = X``: the far-field root, or 1 for the break table.
    """

    def __init__(self, K, n, X=6.0, n_cells=800, grade=2.0, p=2.0, limiter=vanleer,
                 rxn_sign=1.0, feed=True, v_const=False, far_field="root"):
        self.K, self.n, self.p, self.X = K, n, p, X
        self.limiter, self.rxn_sign, self.feed = limiter, rxn_sign, feed
        self.shape = (n_cells, 1)                     # one field: never a bare (n,)
        t = np.linspace(0.0, 1.0, n_cells + 1)
        self.x_f = X * t ** grade                     # cells packed at the exit x = 0
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.g_inf = gamma_far(K, n, p) if far_field == "root" else 1.0
        # Outward normal: a dgamma/dn + b gamma = d at each end.
        #   x = 0, the exit (outflow, v < 0): no physical condition belongs here.
        #     dgamma/dx = 0 is a numerical closure: it sets the outflow face value
        #     by extrapolation from the cells (the true exit slope is K gamma(0)^p,
        #     not 0), and the discrete balance closes with that face value.
        #                                                -> a = 1, b = 0, d = 0
        #   x = X, the entrance (lambda -> infinity, inflow): gamma = gamma_inf
        #                                                -> a = 0, b = 1, d = gamma_inf
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0},
                   {"a": 0.0, "b": 1.0, "d": self.g_inf})
        w_f = np.ones_like(self.x_f) if v_const else W_rtd(self.x_f, n)
        self.v = -w_f.reshape(-1, 1)                  # negative: flow toward x = 0
        self.w_c = W_rtd(self.x_c, n).reshape(self.shape)
        self.E_c = E_rtd(self.x_c, n).reshape(self.shape)
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                      self.bc, v=self.v, axis=0)
        self.div = construct_div(self.shape, self.x_f, nu=0, axis=0)  # nu=0: Cartesian
        self.conv_mat = conv_mat.tocsc()
        self.conv_bc = np.asarray(conv_bc.todense()).reshape(-1, 1)
        self.jac_up = (self.div @ self.conv_mat).tocsc()
        self.g_up = np.asarray(self.div @ self.conv_bc).reshape(-1)
        self.tvd = (None if limiter is None else
                    ConvFluxTVD(self.shape, self.x_f, self.x_c, self.bc,
                                tvd_limiter=limiter, axis=0))

    # -- the physics -------------------------------------------------------------
    def source(self, g):
        """Source ``S = E(x) - K w gamma |gamma|^(p-1)`` and its diagonal derivative.

        Parameters
        ----------
        g : numpy.ndarray
            ``gamma`` at the cell centres, shape ``(n_cells, 1)``.

        Returns
        -------
        tuple of numpy.ndarray
            ``(S, dS/dgamma)``, both of shape ``(n_cells, 1)``; exact, no
            finite differences.
        """
        a = np.abs(g)
        rate = self.rxn_sign * self.K * self.w_c * g * a ** (self.p - 1.0)
        drate = self.rxn_sign * self.K * self.w_c * self.p * a ** (self.p - 1.0)
        s = (self.E_c - rate) if self.feed else -rate
        return s, -drate

    def residual(self, x):
        """Steady residual ``div(flux) - S`` and its exact Jacobian.

        Parameters
        ----------
        x : numpy.ndarray
            State, any array with ``n_cells`` elements.

        Returns
        -------
        tuple
            ``(residual, jacobian)``: a vector of ``n_cells`` entries and a sparse
            CSC matrix.
        """
        g = np.asarray(x).reshape(self.shape)
        s, ds = self.source(g)
        jac_s = construct_coefficient_matrix(ds, shape=self.shape)
        if self.tvd is None:
            res = self.jac_up @ g.reshape(-1) + self.g_up - s.reshape(-1)
            return res, (self.jac_up - jac_s).tocsc()
        flux, jac_f = self.tvd(g, self.v)
        res = self.div @ flux.reshape(-1) - s.reshape(-1)
        return res, (self.div @ jac_f - jac_s).tocsc()

    # -- the two solver routes ---------------------------------------------------
    def solve(self, tol=1e-11, maxfev=60):
        """Solve by Newton with the exact Jacobian of the limited flux.

        Starts from the far-field root, with ``damping="monotone"`` and
        ``stagnation=True``.

        Parameters
        ----------
        tol : float, optional
            Newton step tolerance (absolute; ``gamma`` is of order one).
        maxfev : int, optional
            Iteration limit.

        Returns
        -------
        MaxMixedness
            ``self``, with ``g`` (shape ``(n_cells, 1)``) and ``nit`` set.

        Raises
        ------
        RuntimeError
            When Newton does not succeed.
        """
        g0 = np.full(self.shape, self.g_inf)
        sol = newton(self.residual, g0, tol=tol, maxfev=maxfev, damping="monotone",
                     stagnation=True)
        if not sol.success:
            raise RuntimeError(f"maximum-mixedness Newton solve failed: {sol.message}")
        self.g, self.nit = np.asarray(sol.x).reshape(self.shape), int(sol.nit)
        return self

    def solve_deferred(self, tol=1e-13, max_sweeps=200, maxfev=60):
        """Solve by van Leer deferred correction: the limiter correction lagged.

        Each sweep solves the upwind equations with the previous sweep's TVD
        correction ``div(v dgamma_f)`` in the residual, by Newton with the exact
        upwind-plus-source Jacobian. A second solver route for the same discrete
        equations as :meth:`solve`.

        Parameters
        ----------
        tol : float, optional
            Largest change of a cell value between sweeps at which to stop.
        max_sweeps : int, optional
            Sweep limit.
        maxfev : int, optional
            Newton iteration limit per sweep.

        Returns
        -------
        MaxMixedness
            ``self``, with ``g`` and ``nit`` (the number of sweeps) set.

        Raises
        ------
        RuntimeError
            When a sweep's Newton solve fails or the sweeps do not converge.
        """
        if self.limiter is None:
            raise ValueError("deferred correction needs a limiter")
        corr = np.zeros(self.shape[0])

        def residual(x):
            g = np.asarray(x).reshape(self.shape)
            s, ds = self.source(g)
            jac_s = construct_coefficient_matrix(ds, shape=self.shape)
            res = self.jac_up @ g.reshape(-1) + self.g_up - s.reshape(-1) + corr
            return res, (self.jac_up - jac_s).tocsc()

        g = np.full(self.shape, self.g_inf)
        change = np.inf
        for sweep in range(1, max_sweeps + 1):
            sol = newton(residual, g, tol=tol, maxfev=maxfev, damping="monotone",
                         stagnation=True)
            if not sol.success:
                raise RuntimeError(f"deferred-correction sweep {sweep}: {sol.message}")
            g_new = np.asarray(sol.x).reshape(self.shape)
            change = float(np.max(np.abs(g_new - g)))
            g = g_new
            if change < tol:
                self.g, self.nit = g, sweep
                return self
            _, dg = interp_cntr_to_stagg_tvd(g, self.x_f, self.x_c, self.bc, self.v,
                                             tvd_limiter=self.limiter, axis=0)
            corr = np.asarray(self.div @ (self.v * dg.reshape(-1, 1))).reshape(-1)
        raise RuntimeError(f"deferred correction did not reach {tol:g} in {max_sweeps} "
                           f"sweeps (last change {change:.2e})")

    # -- reading the solution ----------------------------------------------------
    def exit_value(self):
        """``gamma`` on the exit face ``x = 0``, from ``compute_boundary_values`` (P8)."""
        return float(np.asarray(compute_boundary_values(self.g, self.x_f, self.x_c,
                                                        self.bc)[0]).ravel()[0])

    def exit_cell(self):
        """``gamma`` at the first cell centre, half a cell short of the exit."""
        return float(self.g[0, 0])

    def balance_gap(self, exit_value=None):
        """Discrete balance: outflow at the exit minus inflow plus net source.

        ``w(0) gamma(0) - w(X) gamma_inf - int (E - K w gamma^p) dx`` over the
        cells, with the exit value read on the face unless one is given.

        Parameters
        ----------
        exit_value : float, optional
            The exit value to test; the face value by default.

        Returns
        -------
        float
            Zero to round-off for a conservative scheme read on the face.
        """
        ge = self.exit_value() if exit_value is None else exit_value
        s, _ = self.source(self.g)
        dx = np.diff(self.x_f)
        out = -self.v[0, 0] * ge
        inflow = -self.v[-1, 0] * self.g_inf
        return float(out - inflow - np.sum(dx * s.ravel()))


def gamma_mm(K, n, tol=1e-11, **kw):
    """Maximum-mixedness exit value by :class:`MaxMixedness` (Newton, face read).

    Parameters
    ----------
    K : float
        ``k c0^(p-1) tau``.
    n : int
        Number of tanks.
    tol : float, optional
        Newton step tolerance, passed to :meth:`MaxMixedness.solve`.
    **kw
        Passed to :class:`MaxMixedness`.

    Returns
    -------
    float
        ``gamma(0)``.
    """
    return MaxMixedness(K, n, **kw).solve(tol=tol).exit_value()


__all__.append("gamma_mm")


def mm_ode(K, n, X=6.0, p=2.0, rtol=1e-12, atol=1e-14, g_start=None, backward=True,
           ratio=None):
    """Maximum mixedness from eq. (31) as an initial-value problem (DOP853).

    ``d gamma/dx = K gamma^p + (E/w)(gamma - 1)``, integrated from ``x = X``
    down to the exit. No grid, no flux, no pymrm operator and no Newton solve;
    explicit, which is legitimate because the equation contracts in the
    backward direction and is not stiff there.

    Parameters
    ----------
    K : float
        ``k c0^(p-1) tau``.
    n : int
        Number of tanks.
    X : float, optional
        Starting life expectation.
    p : float, optional
        Reaction order.
    rtol, atol : float, optional
        Tolerances of ``solve_ivp``.
    g_start : float, optional
        Starting value; the far-field root by default.
    backward : bool, optional
        Integrate from ``X`` to 0 (the model); ``False`` integrates forward from
        0 to ``X`` and exists for the break table.
    ratio : callable, optional
        Replaces ``E/w``; for the break table.
    Returns
    -------
    float
        ``gamma`` at the end of the integration, or ``nan`` when ``solve_ivp``
        does not succeed.
    """
    g0 = gamma_far(K, n, p) if g_start is None else g_start
    rr = ratio if ratio is not None else (
        lambda x: float(E_rtd(np.atleast_1d(x), n)[0] / W_rtd(x, n)))
    rhs = lambda x, y: K * max(y[0], 0.0) ** p + rr(x) * (y[0] - 1.0)
    s = solve_ivp(rhs, (X, 0.0) if backward else (0.0, X), [g0], method="DOP853",
                  rtol=rtol, atol=atol)
    return float(s.y[0, -1]) if s.success else np.nan


# ---- the degree of segregation ----------------------------------------------------
def var_alpha(n):
    """Variance of the age distribution of ``n`` tanks, his eq. (11), over ``tau^2``."""
    return (n + 1) * (n + 5) / (12 * n * n)


def J_real_mixed(n):
    """``J`` of ``n`` real ideally mixed vessels, ``(n - 1)/(n + 5)``."""
    return (n - 1) / (n + 5)


def J_real_seg(n):
    """``J`` of ``n`` real internally segregated vessels, ``1 - 6(n-1)/((n+1)(n+5))``."""
    return 1 - 6 * (n - 1) / ((n + 1) * (n + 5))


def alpha_P(lam, n, U=60.0):
    """Mean age of the points with life expectation ``lam``, his eq. (27).

    Parameters
    ----------
    lam : float
        Life expectation over ``tau``.
    n : int
        Number of tanks.
    U : float, optional
        Upper limit replacing infinity in the quadrature.

    Returns
    -------
    float
        ``int_lam^U w ds / w(lam)``.
    """
    return quad(lambda s: W_rtd(s, n), lam, U, limit=500)[0] / W_rtd(lam, n)


def J_mm(n, abar=None, U=60.0):
    """``J`` under maximum mixedness, his eq. (28).

    Parameters
    ----------
    n : int
        Number of tanks.
    abar : float, optional
        Mean age; ``(n + 1)/(2n)`` by default.
    U : float, optional
        Upper limit of the inner quadrature.

    Returns
    -------
    float
        ``int (alpha_P - abar)^2 w dlam / var(alpha)``.
    """
    ab = (n + 1) / (2 * n) if abar is None else abar
    v = quad(lambda l: (alpha_P(l, n, U) - ab) ** 2 * W_rtd(l, n), 0, 30, limit=600)[0]
    return v / var_alpha(n)


def var_within_mm(n, U=60.0, density=None):
    """Variance of the ages within a point under maximum mixedness (Appendix I).

    Parameters
    ----------
    n : int
        Number of tanks.
    U : float, optional
        Upper limit of the inner quadratures.
    density : callable, optional
        Replaces the frequency function ``E_rtd(., n)`` in the within-point age
        distribution ``phi_P(a) = f(lam + a)/(1 - F(lam))``; for the break table.

    Returns
    -------
    float
        ``int (<a^2>_P - alpha_P^2) w dlam``, over ``tau^2``.
    """
    f = (lambda y: E_rtd(np.array([y]), n)[0]) if density is None else density

    def inner(lam):
        m2 = quad(lambda a: a * a * f(lam + a) / W_rtd(lam, n), 0, U, limit=500)[0]
        return (m2 - alpha_P(lam, n, U) ** 2) * W_rtd(lam, n)
    return quad(inner, 0, 30, limit=400)[0]
