"""The Graetz problem: a 2-D finite-volume solve and the eigenfunction series.

Model code of gallery page A3.15 (Graetz 1882). Two routes to the same
constants that share no code: :class:`GraetzTube` solves the parabolic PDE on
the pymrm operators and never forms an eigenfunction; :func:`eigen` and
:func:`amplitudes` shoot the radial ODE and project the inlet profile onto its
eigenfunctions, with no grid and no pymrm.

Everything is dimensionless: ``s = r/R``, ``x = a^2 z/(2 alpha R^2)`` with
``alpha`` the mean velocity, and ``theta = (T - T_0)/(T_1 - T_0)``. In these
variables the problem has no parameter left:

    (1 - s^2) dtheta/dx = (1/s) d/ds(s dtheta/ds) + (1/Pe^2) d2theta/dx2,

``theta(0, s) = 1``, ``theta(x, 1) = 0``; ``1/Pe^2 = 0`` is Graetz's equation.

Contents
--------
GraetzTube
    The 2-D solve: van Leer TVD convection along ``x`` by Newton with the exact
    Jacobian (``ConvFluxTVD``), and the same discrete equations by deferred
    correction as a second solver route.
decay
    ``(mu_1, p_1)`` from the far-field slope and intercept of a mean of theta.
V_profile, V_end, eigen, amplitudes
    The eigenfunction route: radial ODE, eigenvalues, cup and area amplitudes.
series_coeffs
    The power series of ``V(mu) = 0`` in exact rational arithmetic, with the
    Graetz weight or with it dropped (Bessel).
theta_series, invert
    The cup-mean outlet temperature from a set of ``(mu_i, p_i)``, and its
    inversion for ``k/c`` at one run.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from pymrm import (ConvFluxTVD, compute_boundary_values, construct_convflux_upwind,
                   construct_div, construct_grad, interp_cntr_to_stagg_tvd, newton, vanleer)

__all__ = ["GraetzTube", "decay", "V_profile", "V_end", "eigen", "amplitudes",
           "series_coeffs", "theta_series", "invert"]


class GraetzTube:
    """Steady 2-D Graetz problem on ``x`` in ``[0, X]``, ``s`` in ``[0, 1]``.

    The state ``theta`` has shape ``(n_z, n_r)``: the axial (marching) axis
    first, the radial axis second, one field. All operators are assembled once
    in the constructor; :meth:`solve` and :meth:`solve_deferred` only solve.

    Parameters
    ----------
    X : float, optional
        Length of the domain in ``x``.
    n_z, n_r : int, optional
        Number of axial and radial cells.
    nu_r : int, optional
        Geometry of the radial divergence: 1 cylindrical (the tube), 0 a slab.
    inv_pe2 : float, optional
        ``1/Pe^2``, the axial conduction term; 0 is Graetz's equation.
    plug : bool, optional
        Uniform velocity instead of Poiseuille (Graetz's first model).
    z_grade, r_grade : float, optional
        Grading of the faces: ``x = X t^z_grade`` packs cells at the inlet,
        ``s = 1 - (1 - t)^r_grade`` packs them at the wall.
    limiter : callable or None, optional
        TVD limiter of the axial convective face values; ``None`` is
        first-order upwind (one linear solve).
    """

    def __init__(self, X=1.0, n_z=600, n_r=128, nu_r=1, inv_pe2=0.0, plug=False,
                 z_grade=3.0, r_grade=2.0, limiter=vanleer):
        self.X, self.shape, self.limiter = X, (n_z, n_r), limiter
        t = np.linspace(0.0, 1.0, n_z + 1)
        self.x_f = X * t ** z_grade
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        t = np.linspace(0.0, 1.0, n_r + 1)
        self.s_f = 1.0 - (1.0 - t) ** r_grade
        self.s_c = 0.5 * (self.s_f[:-1] + self.s_f[1:])
        # u(s) over twice the mean velocity; the 2 is absorbed into x
        self.w = np.ones_like(self.s_c) if plug else 1.0 - self.s_c ** 2
        self.v = np.broadcast_to(self.w, (n_z + 1, n_r))     # on the axial faces

        # a dtheta/dn + b theta = d, with n the OUTWARD normal of each end.
        #   inlet  x = 0: theta = 1 (the feed temperature)   -> a = 0, b = 1, d = 1
        #   outlet x = X: dtheta/dx = 0 (n = +x there)        -> a = 1, b = 0, d = 0
        #     With 1/Pe^2 = 0 the outlet needs no condition; this one only closes
        #     the outflow face value and reaches no cell upstream.
        self.bc_z = ({"a": 0.0, "b": 1.0, "d": 1.0},
                     {"a": 1.0, "b": 0.0, "d": 0.0})
        #   axis   s = 0: dtheta/ds = 0, symmetry (n = -s; zero flux either way)
        #                                                     -> a = 1, b = 0, d = 0
        #   wall   s = 1: theta = 0, the wall at T_0          -> a = 0, b = 1, d = 0
        self.bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0},
                     {"a": 0.0, "b": 1.0, "d": 0.0})

        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c, self.bc_z,
                                                  v=self.v, axis=0)
        grad_z, grad_z_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc_z, axis=0)
        self.div_z = construct_div(self.shape, self.x_f, nu=0, axis=0)    # nu=0: Cartesian
        grad_r, grad_r_bc = construct_grad(self.shape, self.s_f, self.s_c, self.bc_r, axis=1)
        div_r = construct_div(self.shape, self.s_f, nu=nu_r, axis=1)     # nu=1: cylindrical
        # conduction (radial, and axial when 1/Pe^2 > 0): linear, constant
        self.jac_cond = (-(div_r @ grad_r) - inv_pe2 * (self.div_z @ grad_z)).tocsc()
        self.g_cond = np.asarray((-(div_r @ grad_r_bc) - inv_pe2 * (self.div_z @ grad_z_bc))
                                 .todense()).ravel()
        # first-order upwind convection: the whole operator is linear
        self.jac_up = (self.div_z @ conv + self.jac_cond).tocsc()
        self.g_up = np.asarray((self.div_z @ conv_bc).todense()).ravel() + self.g_cond
        self.tvd = (None if limiter is None else
                    ConvFluxTVD(self.shape, self.x_f, self.x_c, self.bc_z,
                                tvd_limiter=limiter, axis=0))
        area = 2.0 * self.s_c * np.diff(self.s_f)                # area weights, sum 1
        self.area_w = area
        self.cup_w = area * self.w / np.sum(area * self.w)       # flow weights, sum 1

    def residual(self, x):
        """Steady residual ``div(v theta_f) - conduction`` and its exact Jacobian.

        Parameters
        ----------
        x : numpy.ndarray
            State with ``n_z * n_r`` elements.

        Returns
        -------
        tuple
            ``(residual, jacobian)``: a vector and a sparse CSC matrix.
        """
        th = np.asarray(x).reshape(self.shape)
        flux, jac_f = self.tvd(th, self.v)
        res = self.div_z @ flux.reshape(-1) + self.jac_cond @ th.reshape(-1) + self.g_cond
        return res, (self.div_z @ jac_f + self.jac_cond).tocsc()

    def solve(self, tol=1e-11, maxfev=30):
        """Solve by Newton with the exact Jacobian of the limited flux.

        Starts from the inlet value ``theta = 1`` everywhere, with
        ``damping="monotone"`` and ``stagnation=True``. Without a limiter the
        equations are linear and one sparse LU solve gives the answer.

        Parameters
        ----------
        tol : float, optional
            Newton step tolerance (absolute; theta is of order one).
        maxfev : int, optional
            Iteration limit.

        Returns
        -------
        dict
            ``theta`` (shape ``(n_z, n_r)``), ``nit`` (Newton iterations, 0 for
            the linear upwind solve), ``cup`` and ``area`` (the flow-weighted
            and area-weighted cross-section means along ``x``).

        Raises
        ------
        RuntimeError
            When Newton does not succeed.
        """
        if self.tvd is None:
            th = splu(self.jac_up).solve(-self.g_up).reshape(self.shape)
            return self._result(th, 0)
        sol = newton(self.residual, np.ones(self.shape), tol=tol, maxfev=maxfev,
                     damping="monotone", stagnation=True)
        if not sol.success:
            raise RuntimeError(f"Graetz Newton solve failed: {sol.message}")
        return self._result(np.asarray(sol.x).reshape(self.shape), int(sol.nit))

    def solve_deferred(self, omega=0.6, tol=1e-11, maxit=200, raise_on_stall=True):
        """Solve the same equations by van Leer deferred correction.

        One LU of the first-order upwind operator; each sweep puts the
        previous sweep's limiter correction ``div(v dtheta_f)`` in the
        right-hand side and relaxes the update by ``omega``. A second solver
        route for the discrete equations of :meth:`solve`.

        Parameters
        ----------
        omega : float, optional
            Relaxation of the update; 1 is plain deferred correction.
        tol : float, optional
            Largest change of a cell value per unrelaxed sweep at which to stop.
        maxit : int, optional
            Sweep limit.
        raise_on_stall : bool, optional
            Raise when ``maxit`` is reached; ``False`` returns the last sweep.

        Returns
        -------
        dict
            As :meth:`solve`, with ``nit`` the number of sweeps and ``incr``
            the last unrelaxed change.

        Raises
        ------
        RuntimeError
            When the sweeps do not converge and ``raise_on_stall`` is set.
        """
        if self.tvd is None:
            raise ValueError("deferred correction needs a limiter")
        lu = splu(self.jac_up)
        th = lu.solve(-self.g_up).reshape(self.shape)
        incr = np.inf
        for nit in range(1, maxit + 1):
            _, dth = interp_cntr_to_stagg_tvd(th, self.x_f, self.x_c, self.bc_z, self.v,
                                              self.limiter, axis=0)
            corr = np.asarray(self.div_z @ (self.v * dth).reshape(-1, 1)).ravel()
            new = th + omega * (lu.solve(-self.g_up - corr).reshape(self.shape) - th)
            incr = float(np.max(np.abs(new - th))) / omega
            th = new
            if incr < tol:
                break
        else:
            if raise_on_stall:
                raise RuntimeError(f"deferred correction did not converge: {incr:.2e}")
        out = self._result(th, nit)
        out["incr"] = incr
        return out

    def _result(self, th, nit):
        return dict(theta=th, nit=nit, x=self.x_c, cup=th @ self.cup_w, area=th @ self.area_w)

    def nusselt(self, theta):
        """Local Nusselt number ``-2 (dtheta/ds)|_wall / theta_cup`` along ``x``.

        The wall gradient is the one the radial operator uses, read with
        ``compute_boundary_values`` (pitfalls P8).

        Parameters
        ----------
        theta : numpy.ndarray
            Solution of shape ``(n_z, n_r)``.

        Returns
        -------
        tuple of numpy.ndarray
            ``(x, Nu)`` at the axial cell centres.
        """
        _, grad_wall = compute_boundary_values(theta, self.s_f, self.s_c, self.bc_r[1],
                                               axis=1, bound_id=1)
        return self.x_c, -2.0 * np.ravel(grad_wall) / (theta @ self.cup_w)


def decay(x, mean, xmin=0.35, xmax=None):
    """``(mu_1, p_1)`` from a straight-line fit of ``ln(mean)`` against ``x``.

    Parameters
    ----------
    x : numpy.ndarray
        Axial positions.
    mean : numpy.ndarray
        A cross-section mean of theta at ``x``.
    xmin : float, optional
        Start of the fit window; higher modes must have decayed there.
    xmax : float, optional
        End of the fit window; ``None`` fits to the end of the domain.

    Returns
    -------
    tuple of float
        ``(sqrt(-slope), exp(intercept))``.
    """
    m = (x > xmin) if xmax is None else (x > xmin) & (x < xmax)
    slope, icept = np.polyfit(x[m], np.log(mean[m]), 1)
    return float(np.sqrt(-slope)), float(np.exp(icept))


# ---- the eigenfunction route: no grid, no pymrm --------------------------------------
def _rhs(mu, weight):
    def rhs(t, y):
        return [y[1], -(y[1] / t if t > 0 else 0.0) - mu ** 2 * weight(t) * y[0]]
    return rhs


def V_profile(mu, n=40001, weight=lambda s: 1.0 - s * s):
    """``V(s)`` and ``V'(s)`` for ``V'' + V'/s + mu^2 w(s) V = 0``, ``V(0) = 1``.

    Parameters
    ----------
    mu : float
        Eigenvalue candidate.
    n : int, optional
        Number of equally spaced output points on ``[0, 1]``.
    weight : callable, optional
        ``w(s)``; ``1 - s^2`` is Poiseuille flow.

    Returns
    -------
    tuple of numpy.ndarray
        ``(s, V, V')``.
    """
    s0 = 1e-10
    w0 = weight(0.0)
    sol = solve_ivp(_rhs(mu, weight), [s0, 1.0],
                    [1 - w0 * mu ** 2 * s0 ** 2 / 4, -w0 * mu ** 2 * s0 / 2],
                    t_eval=np.linspace(s0, 1.0, n), rtol=1e-13, atol=1e-15)
    return sol.t, sol.y[0], sol.y[1]


def V_end(mu, rtol=1e-12, weight=lambda s: 1.0 - s * s):
    """``V(1)`` alone, without output sampling.

    Parameters
    ----------
    mu : float
        Eigenvalue candidate.
    rtol : float, optional
        Relative tolerance of ``solve_ivp``.
    weight : callable, optional
        ``w(s)``.

    Returns
    -------
    float
        ``V(1; mu)``; its roots are the eigenvalues.
    """
    s0 = 1e-10
    w0 = weight(0.0)
    return solve_ivp(_rhs(mu, weight), [s0, 1.0],
                     [1 - w0 * mu ** 2 * s0 ** 2 / 4, -w0 * mu ** 2 * s0 / 2],
                     rtol=rtol, atol=1e-15).y[0, -1]


def eigen(n_modes=10, hi=40.0, weight=lambda s: 1.0 - s * s):
    """Roots of ``V(1; mu) = 0``: a coarse scan to bracket, then ``brentq``.

    Parameters
    ----------
    n_modes : int, optional
        Number of eigenvalues returned.
    hi : float, optional
        Upper end of the scan; the roots are about 4 apart.
    weight : callable, optional
        ``w(s)``.

    Returns
    -------
    numpy.ndarray
        The first ``n_modes`` eigenvalues ``mu_i``.
    """
    grid = np.linspace(0.5, hi, 300)        # spacing at most 0.24; roots ~4 apart
    val = np.array([V_end(m, rtol=1e-9, weight=weight) for m in grid])
    roots = [brentq(V_end, grid[i], grid[i + 1], args=(1e-12, weight), xtol=1e-11)
             for i in np.nonzero(val[:-1] * val[1:] < 0)[0]]
    return np.array(roots[:n_modes])


def amplitudes(mus, n=40001):
    """Cup-mean and area-mean expansion coefficients of a uniform inlet.

    Parameters
    ----------
    mus : array_like
        Eigenvalues of the Poiseuille problem.
    n : int, optional
        Number of quadrature points (trapezoidal) on ``[0, 1]``.

    Returns
    -------
    tuple of numpy.ndarray
        ``(p_cup, p_area, identity)``: the amplitudes of the flow-weighted and
        the area-weighted mean, and ``|N_i / (-V'(1)/mu_i^2) - 1|``, which the
        ODE makes zero, so it measures the quadrature.
    """
    p_cup, p_area, ident = [], [], []
    for mu in mus:
        s, V, dV = V_profile(mu, n)
        N = np.trapezoid((1 - s ** 2) * V * s, s)       # = -V'(1)/mu^2, by the ODE
        D = np.trapezoid((1 - s ** 2) * V ** 2 * s, s)
        M = np.trapezoid(V * s, s)
        p_cup.append(4 * N * N / D)          # 4 = 1 / int (1 - s^2) s ds
        p_area.append(2 * (N / D) * M)       # 2 = 1 / int s ds
        ident.append(abs(N / (-dV[-1] / mu ** 2) - 1))
    return np.array(p_cup), np.array(p_area), np.array(ident)


def series_coeffs(N=60, weight_term=True):
    """Coefficients of ``V(mu; s = 1)`` as a polynomial in ``mu^2``, exactly.

    ``V = sum a_k (mu s)^(2k)`` with ``a_k = sum_j c[k][j] mu^(-2j)``, from the
    recurrence of the ODE in rational arithmetic; setting ``s = 1`` collapses
    it to the polynomial in ``mu^2`` that Graetz prints.

    Parameters
    ----------
    N : int, optional
        Order of truncation of the recurrence.
    weight_term : bool, optional
        ``False`` drops the ``-s^2`` of the weight ``1 - s^2``, which turns the
        equation into Bessel's.

    Returns
    -------
    list of float
        The coefficients of ``mu^0`` to ``mu^10``.
    """
    c = [[Fraction(1)]]
    for m in range(N):
        row = [Fraction(0)] * (m + 2)
        for j in range(m + 2):
            t1 = c[m][j] if j < len(c[m]) else Fraction(0)
            t2 = (c[m - 1][j - 1] if (weight_term and m >= 1 and j >= 1
                                      and j - 1 < len(c[m - 1])) else Fraction(0))
            row[j] = -(t1 - t2) / Fraction((2 * m + 2) ** 2)
        c.append(row)
    return [float(sum(c[k][k - p] for k in range(p, N + 1) if k - p < len(c[k])))
            for p in range(6)]


def theta_series(kc, l, G, mus, ps, half=2.0):
    """Cup-mean outlet temperature ``sum p_i exp(-(k/c) pi l mu_i^2 / (half G))``.

    Parameters
    ----------
    kc : float
        ``k/c`` in the paper's units (cm^2 g / min / cm^3, i.e. g/(cm min)).
    l : float
        Tube length in cm.
    G : float
        Mass flow in g/min.
    mus, ps : numpy.ndarray
        Eigenvalues and amplitudes.
    half : float, optional
        2 for Poiseuille (``x = (k/c) pi z / 2G``), 1 for Graetz's plug-flow
        variables.

    Returns
    -------
    float
        ``theta`` at the outlet.
    """
    return float(np.sum(ps * np.exp(-kc * np.pi * l * mus ** 2 / (half * G))))


def invert(th, l, G, mus, ps, half=2.0):
    """``k/c`` at which :func:`theta_series` equals the measured ``th``.

    Parameters
    ----------
    th : float
        Measured ``(U - T_0)/(T_1 - T_0)``.
    l, G : float
        Tube length (cm) and mass flow (g/min).
    mus, ps : numpy.ndarray
        Eigenvalues and amplitudes.
    half : float, optional
        As in :func:`theta_series`.

    Returns
    -------
    float
        The root, by ``brentq``.
    """
    return brentq(lambda kc: theta_series(kc, l, G, mus, ps, half) - th,
                  1e-8, 20.0, xtol=1e-15, rtol=1e-15)
