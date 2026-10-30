"""Tanks in series and the closed dispersion vessel, assembled from pymrm operators.

Model code of gallery page A2.4 (Levenspiel, *Chemical Reaction Engineering*,
3rd ed., 1999, Chapter 14). The page imports these; the closed forms it audits
(Levenspiel's Eqs. 3, 8, 9, 10 and 11) stay in the notebook, so that the pymrm
solve and the book's algebra never share a line of code.

Contents
--------
Cascade
    ``N`` stirred tanks in series as the ``N``-cell first-order upwind scheme,
    with a pure-outflow exit (``{"outflow": True}``, pitfalls P6).
dispersion_exit
    Exit concentration of a closed axial-dispersion vessel with a first-order
    reaction, TVD convection solved by Newton with the exact Jacobian
    (``ConvFluxTVD``, pitfalls P10).
dispersion_exit_deferred
    The same vessel by deferred correction, the page's second solver route.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import eye_array, identity
from scipy.sparse.linalg import splu

from pymrm import (ConvFluxTVD, NumJac, compute_boundary_values,
                   construct_coefficient_matrix, construct_convflux_upwind,
                   construct_div, construct_grad, interp_cntr_to_stagg_tvd,
                   newton, vanleer)

__all__ = ["Cascade", "dispersion_exit", "dispersion_exit_deferred"]


class Cascade:
    """``N`` stirred tanks in series, assembled from the pymrm operators.

    A tracer balance on tank ``i``, ``dC_i/dt = (C_{i-1} - C_i)/t_i``, is
    exactly the first-order upwind finite-volume scheme for plug flow on ``N``
    cells, so the cell count IS Levenspiel's ``N`` and nothing here is a
    discretisation of anything. The line is ``0 <= z <= 1`` with velocity
    ``u = 1/tbar``, so the transit time of the whole line is ``tbar``.

    Parameters
    ----------
    N : int
        Number of tanks (cells).
    tbar : float, optional
        Mean residence time of the whole chain.
    c_in : float, optional
        Feed concentration, imposed on the inlet face.
    nu : int, optional
        Geometry of ``construct_div``: 0 Cartesian (the model); other values
        exist only for the page's break table.
    outlet : {"outflow", "stirred", "zerograd"}, optional
        Exit face. ``"outflow"`` (default) is pymrm's pure-outflow boundary:
        the face carries the last tank's own value, which is what a stirred
        tank discharges. ``"stirred"`` is the equivalent construction the page
        used before pymrm had that boundary: no operator flux on the exit face
        and the outflow ``u C_N / dz_N`` as a sink on the last tank; it is
        kept to show that the two are the same operator. ``"zerograd"`` is the
        usual PDE outflow ``dc/dz = 0``, which reconstructs the exit face to
        second order and is WRONG for this model (pitfalls P6); it exists for
        the break table.
    faces : array_like, optional
        Face positions on ``[0, 1]``; cell ``i`` then holds a tank of residence
        time ``tbar * (faces[i+1] - faces[i])``. Default: ``N`` equal tanks.
    v_sign : float, optional
        Sign of the velocity; ``-1`` exists only for the break table.

    Attributes
    ----------
    jac_const : scipy.sparse.csc_array
        Constant operator ``div @ conv`` (plus the sink for ``"stirred"``).
    g_const : numpy.ndarray
        Constant boundary contribution, shape ``(N,)``.
    t_i : numpy.ndarray
        Residence time of each tank.
    """

    def __init__(self, N, tbar=1.0, c_in=0.0, nu=0, outlet="outflow",
                 faces=None, v_sign=1.0):
        self.N, self.tbar, self.c_in, self.outlet = N, tbar, c_in, outlet
        self.shape = (N, 1)                      # (cells, fields), never (N,)
        self.z_f = (np.linspace(0.0, 1.0, N + 1) if faces is None
                    else np.asarray(faces, float))
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.dz = np.diff(self.z_f)
        self.u = v_sign / tbar                   # transit time of the whole line = tbar
        self.t_i = self.dz / abs(self.u)         # residence time of each tank
        self._build_operators(nu)

    def _build_operators(self, nu):
        # inlet at z = 0: the feed concentration on the face, c = c_in -> a=0, b=1, d=c_in
        bc_in = {"a": 0.0, "b": 1.0, "d": self.c_in}
        if self.outlet == "outflow":
            # exit at z = 1: a stirred tank discharges its own contents, so the
            # face value is the last cell's value (pitfalls P6)
            bc_out = {"outflow": True}
        elif self.outlet == "stirred":
            # exit at z = 1: c = 0 on the face, so the operator puts no flux there;
            # the physical outflow u C_N is added below as a sink on the last tank
            bc_out = {"a": 0.0, "b": 1.0, "d": 0.0}
        elif self.outlet == "zerograd":
            # exit at z = 1: dc/dz = 0, n = +z -> a=1, b=0, d=0. Right for a
            # discretised PDE, wrong for this model.
            bc_out = {"a": 1.0, "b": 0.0, "d": 0.0}
        else:
            raise ValueError(f"unknown outlet {self.outlet!r}")
        self.bc = (bc_in, bc_out)
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                      self.bc, v=self.u)
        div_mat = construct_div(self.shape, self.z_f, nu=nu)    # nu=0: Cartesian
        jac_conv = div_mat @ conv_mat
        if self.outlet == "stirred":
            coef = np.zeros(self.shape)
            coef[-1, 0] = abs(self.u) / self.dz[-1]
            jac_conv = jac_conv + construct_coefficient_matrix(coef, shape=self.shape)
        self.jac_const = jac_conv.tocsc()
        self.g_const = np.asarray((div_mat @ conv_bc).todense()).ravel()

    # ---- steady state ---------------------------------------------------
    def first_order(self, k):
        """Steady exit concentration for the rate ``r = k C``.

        Parameters
        ----------
        k : float
            First-order rate constant, in units of 1/time.

        Returns
        -------
        c_exit : float
            Concentration in the last tank, which is the exit stream.
        c : numpy.ndarray
            Concentration in every tank, shape ``(N,)``.
        """
        jac = (self.jac_const + k * eye_array(self.N, format="csc")).tocsc()
        c = splu(jac).solve(-self.g_const)
        return c[-1], c

    def reaction(self, c, k, order):
        """Consumption rate ``k C**order`` per tank (negative concentrations clipped)."""
        return k * np.maximum(c, 0.0) ** order

    def order_n(self, k, order=2.0, c_guess=None):
        """Steady exit concentration for ``r = k C**order``, by Newton.

        Parameters
        ----------
        k : float
            Rate constant.
        order : float, optional
            Reaction order.
        c_guess : float, optional
            Uniform initial guess; default the feed concentration.

        Returns
        -------
        c_exit : float
            Concentration in the last tank.
        c : numpy.ndarray
            Concentration in every tank, shape ``(N,)``.
        sol : pymrm newton result
            The solver result (``nit``, ``success``).
        """
        numjac = NumJac(self.shape)      # (N, 1): last axis is the field, block-diagonal

        def residual(c):
            c = np.asarray(c).reshape(self.shape)
            g_rxn, jac_rxn = numjac(lambda x: self.reaction(x, k, order), c)
            g = self.jac_const @ c.ravel() + self.g_const + np.asarray(g_rxn).ravel()
            return g, (self.jac_const + jac_rxn).tocsc()

        guess = np.full(self.shape, self.c_in if c_guess is None else c_guess)
        sol = newton(residual, guess, tol=1e-14)
        c = np.asarray(sol.x).ravel()
        return c[-1], c, sol

    # ---- transient ------------------------------------------------------
    def pulse_response(self, tmax, nt, theta=0.5):
        """Exit-age distribution ``E(t)`` from a unit pulse into tank 1.

        Parameters
        ----------
        tmax : float
            End of the time window.
        nt : int
            Number of fixed time steps.
        theta : float, optional
            Theta scheme: 0.5 Crank-Nicolson (second order), 1 implicit Euler
            (first order).

        Returns
        -------
        ts, es : numpy.ndarray
            Times and ``E = u C_N``. One mole is injected into tank 1, so the
            curve carries no numerical normalisation and its area is a check.
        """
        c = np.zeros(self.N)
        c[0] = 1.0 / self.dz[0]
        dt = tmax / nt
        lu = splu((identity(self.N, format="csc") + theta * dt * self.jac_const).tocsc())
        ts = np.empty(nt + 1)
        es = np.empty(nt + 1)
        ts[0], es[0] = 0.0, abs(self.u) * c[-1]
        for j in range(nt):
            c = lu.solve(c - (1 - theta) * dt * (self.jac_const @ c) - dt * self.g_const)
            ts[j + 1], es[j + 1] = (j + 1) * dt, abs(self.u) * c[-1]
        return ts, es


def _dispersion_operators(Pe, ktbar, n):
    """Grid, bc and the linear (diffusion, reaction, upwind) parts of the vessel."""
    shape = (n, 1)
    z_f = np.linspace(0.0, 1.0, n + 1)
    z_c = 0.5 * (z_f[:-1] + z_f[1:])
    d_ax = 1.0 / Pe
    # inlet at z = 0 (Danckwerts): u c - D dc/dz = u c_in with u = 1, c_in = 1.
    #   n = -z, so dc/dn = -dc/dz and the condition is D dc/dn + c = 1 -> a=D, b=1, d=1
    # outlet at z = 1 (closed vessel): dc/dz = 0. n = +z -> a=1, b=0, d=0
    bc = ({"a": d_ax, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
    conv_mat, conv_bc = construct_convflux_upwind(shape, z_f, z_c, bc, v=1.0)
    grad_mat, grad_bc = construct_grad(shape, z_f, z_c, bc)
    div_mat = construct_div(shape, z_f, nu=0)                  # nu=0: Cartesian
    return shape, z_f, z_c, bc, d_ax, conv_mat, conv_bc, grad_mat, grad_bc, div_mat


def dispersion_exit(Pe, ktbar, n=1500, limiter=vanleer, tol=1e-10, return_info=False):
    """Exit concentration of a closed dispersion vessel with a first-order reaction.

    ``u = 1``, ``L = 1``, ``tbar = 1``; Danckwerts inlet, zero-gradient outlet.
    With a limiter, the TVD face fluxes and their exact Jacobian come from
    ``ConvFluxTVD`` and the steady equations are solved by ``newton`` with
    ``damping="monotone"`` from the upwind solution (pitfalls P10).

    Parameters
    ----------
    Pe : float
        Peclet number ``uL/D``.
    ktbar : float
        Damkoehler number ``k tbar``.
    n : int, optional
        Number of cells.
    limiter : callable or None, optional
        TVD limiter (default van Leer); None gives first-order upwind, solved
        as one linear system.
    tol : float, optional
        Newton step tolerance; the concentrations are of order one. It must
        stay above the round-off floor of the discrete solution: at small
        ``Pe`` on a fine grid the system is ill-conditioned, each iteration
        refactorises, and the steps jitter at about 1e-12 after the residual
        has reached round-off, so ``tol=1e-12`` never stops (measured at
        ``Pe = 0.1``, ``k tbar = 0.1``, 1500 cells). Newton converges
        quadratically, so the last step taken below 1e-10 leaves the iterate
        at round-off.
    return_info : bool, optional
        Also return the number of Newton iterations (0 for upwind).

    Returns
    -------
    float or tuple
        Exit concentration ``c(1)/c_in`` read from the face with
        ``compute_boundary_values`` (pitfalls P8); with ``return_info``,
        ``(c_exit, nit)``.
    """
    (shape, z_f, z_c, bc, d_ax, conv_mat, conv_bc, grad_mat, grad_bc,
     div_mat) = _dispersion_operators(Pe, ktbar, n)
    # first-order upwind: one linear solve (and the initial guess for TVD)
    jac_up = (div_mat @ (conv_mat - d_ax * grad_mat)
              + ktbar * eye_array(n, format="csc")).tocsc()
    g_up = np.asarray((div_mat @ (conv_bc - d_ax * grad_bc)).todense()).ravel()
    c = splu(jac_up).solve(-g_up)
    nit = 0
    if limiter is not None:
        jac_diff_rxn = (-d_ax * (div_mat @ grad_mat)
                        + ktbar * eye_array(n, format="csc")).tocsc()
        g_diff = -d_ax * np.asarray((div_mat @ grad_bc).todense()).ravel()
        tvd = ConvFluxTVD(shape, z_f, z_c, bc, tvd_limiter=limiter)

        def residual(c):
            flux, jac_conv = tvd(c, 1.0)
            g = div_mat @ flux.reshape(-1) + jac_diff_rxn @ c.reshape(-1) + g_diff
            return g, (div_mat @ jac_conv + jac_diff_rxn).tocsc()

        sol = newton(residual, c.reshape(shape), tol=tol, damping="monotone")
        if not sol.success:
            raise RuntimeError(f"TVD Newton solve failed: {sol.message}")
        c, nit = np.asarray(sol.x).ravel(), sol.nit
    c_exit = float(compute_boundary_values(c.reshape(shape), z_f, z_c, bc)[2].ravel()[0])
    return (c_exit, nit) if return_info else c_exit


def dispersion_exit_deferred(Pe, ktbar, n=1500, limiter=vanleer, max_it=200, tol=1e-13,
                             return_info=False):
    """The same vessel by deferred correction (a second solver route).

    The limiter correction is lagged in the right-hand side and the upwind
    operator is factorised once; the iteration is a fixed point, not Newton,
    and for steady problems it can converge slowly or stall (pitfalls P10).

    Parameters
    ----------
    Pe, ktbar, n, limiter
        As for :func:`dispersion_exit`.
    max_it : int, optional
        Cap on the number of correction sweeps.
    tol : float, optional
        Stop when the largest change of a cell value in one sweep is below it.
    return_info : bool, optional
        Also return the number of sweeps.

    Returns
    -------
    float or tuple
        Exit concentration; with ``return_info``, ``(c_exit, sweeps)``.

    Raises
    ------
    RuntimeError
        When ``max_it`` sweeps do not reach ``tol``.
    """
    (shape, z_f, z_c, bc, d_ax, conv_mat, conv_bc, grad_mat, grad_bc,
     div_mat) = _dispersion_operators(Pe, ktbar, n)
    jac = (div_mat @ (conv_mat - d_ax * grad_mat) + ktbar * eye_array(n, format="csc")).tocsc()
    g_const = np.asarray((div_mat @ (conv_bc - d_ax * grad_bc)).todense()).ravel()
    lu = splu(jac)
    c = lu.solve(-g_const)
    sweeps = 0
    if limiter is not None:
        for sweeps in range(1, max_it + 1):
            _, dc = interp_cntr_to_stagg_tvd(c.reshape(shape), z_f, z_c, bc, 1.0,
                                             tvd_limiter=limiter, axis=0)
            c_new = lu.solve(-g_const - np.asarray(div_mat @ dc.reshape(-1, 1)).ravel())
            done = np.max(np.abs(c_new - c)) < tol
            c = c_new
            if done:
                break
        else:
            raise RuntimeError("deferred correction hit its iteration cap")
    c_exit = float(compute_boundary_values(c.reshape(shape), z_f, z_c, bc)[2].ravel()[0])
    return (c_exit, sweeps) if return_info else c_exit
