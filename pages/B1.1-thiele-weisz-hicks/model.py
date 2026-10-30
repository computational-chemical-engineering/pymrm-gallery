"""Catalyst pellet with diffusion and reaction, isothermal (Thiele) and Weisz-Hicks.

Model code of gallery page B1.1 (Thiele 1939; Weisz and Hicks 1962). The page
imports these; the exact isothermal solutions it compares with are tabulated in
its data file and written out in the notebook, so the pymrm solve and the closed
forms never share a line of code.

Contents
--------
Pellet
    Dimensionless pellet ``div(grad y) = phi^2 R(y)``, ``y = 1`` at the surface,
    slab, cylinder or sphere, assembled from the pymrm operators.
wh_rate
    The Weisz-Hicks rate with the Prater relation folded in.
CentrePinned
    The pellet with ``psi = phi^2`` as an unknown and the centre value pinned:
    the residual that makes the multivalued ``eta(phi)`` curve single-valued.
trace_branch
    The whole S-curve by :func:`pymrm.continuation` in the centre value.
trace_branch_own
    The same sweep by a plain warm-start loop, the page's second route.
wh_reference
    Independent reference: one shooting integration from the centre.
locate_folds
    Turning points of the reference curve, bracketed on a scan and refined.
reference_states
    The steady states of the reference at one Thiele modulus, root-found.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse import bmat, csc_array, eye_array
from scipy.sparse.linalg import spsolve

from pymrm import NumJac, construct_div, construct_grad, continuation, newton

__all__ = ["NU", "Pellet", "wh_rate", "CentrePinned", "trace_branch",
           "trace_branch_own", "wh_reference", "locate_folds", "reference_states"]

NU = {"slab": 0, "cylinder": 1, "sphere": 2}


class Pellet:
    """Dimensionless pellet: ``div(grad y) = phi^2 R(y)`` with ``y(1) = 1``.

    The coordinate ``u = r/R`` runs from the centre (0) to the surface (1) on
    ``n_r`` equal cells; ``phi`` is built on the half-thickness (slab) or the
    radius (cylinder, sphere).

    Parameters
    ----------
    geom : {"slab", "cylinder", "sphere"}, optional
        Geometry; sets ``nu`` of ``construct_div``.
    n_r : int, optional
        Number of cells.
    nu : int, optional
        Overrides the geometry index of ``construct_div`` while ``eta_from``
        keeps the factor of ``geom``. Exists only for the page's break table
        (a pellet assembled in the wrong geometry).

    Attributes
    ----------
    jac_diff : scipy.sparse.csc_array
        Constant operator ``div @ (-grad)``.
    g_diff : numpy.ndarray
        Its constant boundary contribution, shape ``(n_r,)``.
    numjac : pymrm.NumJac
        Jacobian of a pointwise rate on shape ``(n_r, 1)``.
    """

    def __init__(self, geom="sphere", n_r=200, nu=None):
        self.geom, self.n_r = geom, n_r
        self.nu = NU[geom] if nu is None else nu
        self._build_grid()
        self._build_operators()

    def _build_grid(self):
        self.shape = (self.n_r, 1)                          # one field: (n_r, 1), P3
        self.r_f = np.linspace(0.0, 1.0, self.n_r + 1)      # faces
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])     # centres

    def _build_operators(self):
        # outward normal: -r at the centre, +r at the surface (pitfalls P1)
        bc = ({"a": 1.0, "b": 0.0, "d": 0.0},   # u = 0: -dy/du = 0, symmetry
              {"a": 0.0, "b": 1.0, "d": 1.0})   # u = 1: y = 1, surface value
        grad_mat, grad_bc = construct_grad(self.shape, self.r_f, self.r_c, bc)
        div_mat = construct_div(self.shape, self.r_f, nu=self.nu)  # 0 slab, 1 cylinder, 2 sphere
        self.flux_mat, self.flux_bc = -grad_mat, -grad_bc
        self.jac_diff = (div_mat @ self.flux_mat).tocsc()
        self.g_diff = (div_mat @ self.flux_bc).toarray().ravel()
        # the surface-flux row: eta is read from the flux the solver balanced
        self.surf = np.asarray(self.flux_mat[[-1], :].todense()).ravel()
        self.surf_bc = float(self.flux_bc.toarray().ravel()[-1])
        # The rate is pointwise in y, so its Jacobian is diagonal; the
        # tridiagonal coupling is in jac_diff, added analytically.
        self.numjac = NumJac(self.shape)

    def eta_from(self, y, phi):
        """Effectiveness factor ``(nu + 1)/phi^2 dy/du`` at the surface.

        Parameters
        ----------
        y : numpy.ndarray
            Cell values, ``n_r`` of them in any shape.
        phi : float
            Thiele modulus.

        Returns
        -------
        float
        """
        dydu = -(self.surf @ np.ravel(y) + self.surf_bc)
        return (NU[self.geom] + 1.0) * dydu / phi**2

    def solve_linear(self, phi):
        """First-order isothermal reaction: one sparse solve.

        Parameters
        ----------
        phi : float
            Thiele modulus.

        Returns
        -------
        y : numpy.ndarray
            Cell values, shape ``(n_r,)``.
        eta : float
            Effectiveness factor.
        """
        a = self.jac_diff + phi**2 * eye_array(self.n_r, format="csc")
        y = spsolve(a.tocsc(), -self.g_diff)
        return y, self.eta_from(y, phi)

    def residual(self, y, phi, rate):
        """Residual and Jacobian of ``div(grad y) - phi^2 R(y) = 0`` (flattened)."""
        g_r, jac_r = self.numjac(rate, np.reshape(y, self.shape))
        g = (self.jac_diff @ np.reshape(y, (-1, 1)) + self.g_diff.reshape((-1, 1))
             + phi**2 * g_r.reshape((-1, 1)))
        return g, self.jac_diff + phi**2 * jac_r

    def solve_nonlinear(self, phi, rate, y_init=None, **newton_kwargs):
        """Steady state at a fixed Thiele modulus by Newton.

        Parameters
        ----------
        phi : float
            Thiele modulus.
        rate : callable
            Pointwise dimensionless rate ``R(y)``.
        y_init : numpy.ndarray, optional
            Initial guess (default: the uniform pellet, ``y = 1``). Where
            several steady states exist, the guess decides which one is found.
        **newton_kwargs
            Options for :func:`pymrm.newton` (default ``maxfev=60``,
            ``tol=1e-11``).

        Returns
        -------
        y : numpy.ndarray
            Cell values, shape ``(n_r,)``.
        eta : float
            Effectiveness factor.
        result : scipy.optimize.OptimizeResult
            The Newton result; check ``result.success``.

        Raises
        ------
        pymrm.NewtonBreakdown
            When no Newton step can be formed.
        """
        y0 = np.ones(self.n_r) if y_init is None else np.array(y_init, float).ravel()
        opts = {"maxfev": 60, "tol": 1e-11} | newton_kwargs
        res = newton(lambda y: self.residual(y, phi, rate), y0, **opts)
        return res.x, self.eta_from(res.x, phi), res


def wh_rate(beta, gamma):
    """Dimensionless Weisz-Hicks rate with the Prater relation folded in.

    ``R(y) = y exp[gamma beta (1 - y) / (1 + beta (1 - y))]``. It is smooth and
    first order through ``y = 0``, so a Newton iterate that overshoots below
    zero is pulled back (pitfalls P9); it needs no clipping and no bound.

    Parameters
    ----------
    beta : float
        Prater number, the largest relative temperature rise.
    gamma : float
        Arrhenius number ``E / (R_g T_s)``.

    Returns
    -------
    callable
        ``rate(y)``, elementwise.
    """
    def rate(y):
        return y * np.exp(gamma * beta * (1 - y) / (1 + beta * (1 - y)))
    return rate


class CentrePinned:
    """The pellet with ``psi = phi^2`` as an unknown and the centre cell pinned.

    Unknowns ``x = [y_1 .. y_n, psi]``; the extra equation is
    ``y_1 - y_c = 0``. Along ``y_c`` the Weisz-Hicks curve is single-valued,
    so continuation in ``y_c`` passes the folds of ``eta(phi)`` without
    turning.

    Parameters
    ----------
    pellet : Pellet
        The discretised pellet (its operators are reused).
    rate : callable
        Pointwise dimensionless rate.
    """

    def __init__(self, pellet, rate):
        self.pellet, self.rate = pellet, rate
        self._pin = csc_array(eye_array(1, pellet.n_r, format="csc"))

    def __call__(self, x, y_c):
        """Residual and Jacobian at centre value ``y_c``, for :func:`pymrm.continuation`."""
        p = self.pellet
        y, psi = x[:-1], x[-1]
        g_r, jac_r = p.numjac(self.rate, y.reshape(p.shape))
        g = np.vstack([
            p.jac_diff @ y.reshape((-1, 1)) + p.g_diff.reshape((-1, 1))
            + psi * g_r.reshape((-1, 1)),
            np.array([[y[0] - y_c]]),
        ])
        jac = bmat([[p.jac_diff + psi * jac_r, csc_array(g_r.reshape((-1, 1)))],
                    [self._pin, None]], format="csc")
        return g, jac

    def point(self, x):
        """``(phi, eta)`` of a solution vector ``x``."""
        phi = float(np.sqrt(x[-1]))
        return phi, self.pellet.eta_from(x[:-1], phi)


def trace_branch(problem, y_centres, psi_init=0.005, **kwargs):
    """Trace the solution branch with :func:`pymrm.continuation` in ``y_c``.

    Parameters
    ----------
    problem : CentrePinned
        Residual in ``(x, y_c)``.
    y_centres : array_like
        Centre values in the order to visit them. Consecutive values must be
        distinct (see the Notes).
    psi_init : float, optional
        Guess of ``phi^2`` at the first centre value.
    **kwargs
        Passed to :func:`pymrm.continuation`; ``newton_kwargs`` defaults to
        ``maxfev=40, tol=1e-11``.

    Returns
    -------
    branch : numpy.ndarray
        Rows ``(phi, eta, y_c)`` of every point reached.
    xs : list of numpy.ndarray
        The solution vectors ``[y, psi]``.
    result : scipy.optimize.OptimizeResult
        The continuation result; ``success`` is False when it stopped early.
    """
    kwargs.setdefault("newton_kwargs", {"maxfev": 40, "tol": 1e-11})
    n = problem.pellet.n_r
    res = continuation(problem, np.concatenate([np.ones(n), [psi_init]]),
                       np.asarray(y_centres, float), **kwargs)
    branch = np.array([(*problem.point(x), p) for p, x in zip(res.params, res.xs)])
    return branch.reshape(-1, 3), res.xs, res


def trace_branch_own(problem, y_centres, psi_init=0.005, **newton_kwargs):
    """The page's own warm-start sweep in ``y_c``, a second route to the branch.

    Each centre value starts from the previous solution (no predictor); a
    failed solve is skipped and the march goes on. The page compares it with
    :func:`trace_branch` point by point.

    Parameters
    ----------
    problem : CentrePinned
        Residual in ``(x, y_c)``.
    y_centres : array_like
        Centre values in the order to visit them.
    psi_init : float, optional
        Guess of ``phi^2`` at the first centre value.
    **newton_kwargs
        Options for :func:`pymrm.newton` (default ``maxfev=40, tol=1e-11``).

    Returns
    -------
    branch : numpy.ndarray
        Rows ``(phi, eta, y_c)`` of the converged points.
    n_failed : int
        Number of centre values whose solve did not converge.
    """
    opts = {"maxfev": 40, "tol": 1e-11} | newton_kwargs
    y_state, psi_state = np.ones(problem.pellet.n_r), psi_init
    branch, n_failed = [], 0
    for y_c in y_centres:
        res = newton(lambda v, _c=y_c: problem(v, _c),
                     np.concatenate([y_state, [psi_state]]), **opts)
        y_state = res.x[:-1]
        if res.success and res.x[-1] > 0:
            psi_state = res.x[-1]
            branch.append((*problem.point(res.x), y_c))
        else:
            n_failed += 1
    return np.array(branch).reshape(-1, 3), n_failed


def wh_reference(beta, gamma, y_c, nu=2):
    """Independent reference: shoot outward from the centre.

    With ``s = phi u`` the equation no longer contains ``phi``, so the value
    of ``s`` at which ``y`` reaches 1 IS ``phi``. That makes ``eta`` a
    single-valued function of the centre value ``y_c``, without continuation.

    Parameters
    ----------
    beta, gamma : float
        Prater and Arrhenius numbers.
    y_c : float
        Centre value, ``0 < y_c < 1``.
    nu : int, optional
        Geometry index (2: sphere).

    Returns
    -------
    tuple of float or None
        ``(phi, eta)``, or None when ``y`` does not reach 1 by ``s = 400``.
    """
    rate = wh_rate(beta, gamma)

    def rhs(s, v):
        y, dy = v
        return [dy, rate(y) - (nu / s) * dy]

    def hit_surface(s, v):
        return v[0] - 1.0
    hit_surface.terminal, hit_surface.direction = True, 1

    s0 = 1e-6                     # series start: y = y_c + R s^2 / (2 (nu + 1))
    v0 = [y_c + rate(y_c) * s0**2 / (2 * (nu + 1)), rate(y_c) * s0 / (nu + 1)]
    sol = solve_ivp(rhs, (s0, 400.0), v0, events=hit_surface,
                    rtol=1e-10, atol=1e-12)
    if not sol.t_events[0].size:
        return None
    s_star = float(sol.t_events[0][0])
    return s_star, (nu + 1.0) * float(sol.y_events[0][0][1]) / s_star


def locate_folds(beta, gamma, y_scan, nu=2, n_folds=2):
    """Turning points of ``phi(y_c)`` on the reference curve.

    The scan brackets every turning point; the ``n_folds`` most prominent are
    refined by bounded scalar minimisation of ``-phi`` (a maximum) or ``phi``
    (a minimum) over the two neighbouring scan intervals, so the result does
    not depend on the scan spacing.

    Parameters
    ----------
    beta, gamma : float
        Prater and Arrhenius numbers.
    y_scan : array_like
        Centre values to scan, increasing.
    nu : int, optional
        Geometry index.
    n_folds : int, optional
        Number of folds to return.

    Returns
    -------
    folds : list of tuple
        ``(phi_fold, y_c_fold, phi_sampled)`` sorted by ``phi_fold``;
        ``phi_sampled`` is the scan point the fold was bracketed at.
    """
    y_scan = np.asarray(y_scan, float)
    pts = [(wh_reference(beta, gamma, v, nu), v) for v in y_scan]
    curve = np.array([(q[0], v) for q, v in pts if q])
    phi, yc = curve[:, 0], curve[:, 1]
    turn = np.where(np.diff(np.sign(np.diff(phi))) != 0)[0] + 1
    # prominence, so that ripple cannot masquerade as a fold
    prom = [min(abs(phi[i] - phi[:i].min()), abs(phi[i] - phi[i:].max())) for i in turn]
    turn = [t for _, t in sorted(zip(prom, turn), reverse=True)[:n_folds]]
    folds = []
    for i in turn:
        sign = -1.0 if phi[i] > phi[i - 1] else 1.0         # -1: a maximum
        opt = minimize_scalar(lambda v: sign * wh_reference(beta, gamma, v, nu)[0],
                              bounds=(yc[i - 1], yc[i + 1]), method="bounded",
                              options={"xatol": 1e-12})
        if not (opt.success and yc[i - 1] < opt.x < yc[i + 1]):
            raise RuntimeError(f"fold near y_c = {yc[i]:.4g} not refined: {opt.message}")
        folds.append((sign * float(opt.fun), float(opt.x), float(phi[i])))
    return sorted(folds)


def reference_states(beta, gamma, phi, y_scan, nu=2):
    """Every steady state of the reference at Thiele modulus ``phi``.

    The scan brackets each crossing of ``phi(y_c) = phi``; ``brentq`` finds
    the centre value to ``xtol = 1e-14``.

    Parameters
    ----------
    beta, gamma, phi : float
        Prater number, Arrhenius number, Thiele modulus.
    y_scan : array_like
        Centre values to scan, increasing.
    nu : int, optional
        Geometry index.

    Returns
    -------
    list of tuple
        ``(eta, y_c)`` per steady state, sorted by ``eta``.
    """
    y_scan = np.asarray(y_scan, float)
    curve = np.array([(q[0], v) for v in y_scan
                      for q in [wh_reference(beta, gamma, v, nu)] if q])
    f = curve[:, 0] - phi
    out = []
    for i in np.where(np.diff(np.sign(f)) != 0)[0]:
        y_c = brentq(lambda v: wh_reference(beta, gamma, v, nu)[0] - phi,
                     curve[i, 1], curve[i + 1, 1], xtol=1e-14, rtol=1e-13)
        out.append((wh_reference(beta, gamma, y_c, nu)[1], y_c))
    return sorted(out)
