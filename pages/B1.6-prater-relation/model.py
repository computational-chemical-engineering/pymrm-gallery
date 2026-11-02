"""The Prater relation: species and energy in a catalyst particle as two fields.

Model code of gallery page B1.6 (Prater 1958, as printed by Weisz and Hicks
1962). The pellet is solved as two coupled fields, concentration and
temperature, and the relation ``theta - 1 = beta (1 - y)`` is never used by the
solver: the page measures it on the solutions.

Contents
--------
FirstOrder, LangmuirHinshelwood, ZeroOrder
    Rate laws pointwise in ``(y, theta)``, each returning the rate and its two
    partial derivatives.
TwoFieldPellet
    ``div(grad y) = phi^2 R``, ``div(grad theta) = -beta phi^2 R`` on a slab,
    cylinder or sphere, with Dirichlet or film (Robin) surfaces, assembled from
    the pymrm operators with the exact Jacobian.
violation
    The identity violation ``max |(theta - 1) - beta (1 - y)| / |beta|``.
march, march_ivp, peak_of
    The transient with Lewis number ``Le``: variable-step BDF2 with a Newton
    solve per step, and the same semi-discrete system by ``solve_ivp`` (BDF).
ReducedPellet
    The single equation with the Prater relation substituted (``B1.1``).
collocation
    The two-field problem by ``scipy.integrate.solve_bvp``, no pymrm.
shooting, shooting_at_phi
    The reduced sphere equation integrated outward from the centre.
prater_beta
    The Prater number from dimensional transport properties.

Several functions take an argument whose default is the correct model and
whose alternative is a defect the page's break table injects.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse import diags_array

from pymrm import (BlockLayout, compute_boundary_values, construct_coefficient_matrix,
                   construct_div, construct_grad, newton, pseudo_transient)

__all__ = ["NU", "FirstOrder", "LangmuirHinshelwood", "ZeroOrder", "TwoFieldPellet",
           "violation", "march", "march_ivp", "peak_of", "ReducedPellet", "collocation",
           "shooting", "shooting_at_phi", "prater_beta"]

NU = {"slab": 0, "cylinder": 1, "sphere": 2}


# ------------------------------------------------------------------ rate laws
class FirstOrder:
    """Weisz-Hicks kinetics written in ``theta``: ``R = y exp[gamma (1 - 1/theta)]``.

    Parameters
    ----------
    gamma : float
        Arrhenius number ``E/(R_g T_s)``.
    """

    def __init__(self, gamma):
        self.gamma = float(gamma)

    def __call__(self, y, theta):
        """Return ``(R, dR/dy, dR/dtheta)``."""
        e = np.exp(self.gamma * (1.0 - 1.0 / theta))
        r = y * e
        return r, e, r * self.gamma / theta**2


class LangmuirHinshelwood:
    """``R = y exp[gamma (1 - 1/theta)] / (1 + K y)^2``, ``K = K_s exp[gamma_a (1/theta - 1)]``.

    A second temperature dependence sits in the adsorption constant.

    Parameters
    ----------
    gamma : float
        Arrhenius number of the rate constant.
    gamma_a : float, optional
        Arrhenius number of the adsorption constant.
    k_ads : float, optional
        Adsorption constant at the surface temperature.
    """

    def __init__(self, gamma, gamma_a=8.0, k_ads=3.0):
        self.gamma, self.gamma_a, self.k_ads = float(gamma), float(gamma_a), float(k_ads)

    def __call__(self, y, theta):
        """Return ``(R, dR/dy, dR/dtheta)``."""
        e = np.exp(self.gamma * (1.0 - 1.0 / theta))
        k = self.k_ads * np.exp(self.gamma_a * (1.0 / theta - 1.0))
        q = 1.0 + k * y
        r = y * e / q**2
        dr_dth = r * self.gamma / theta**2 + 2.0 * r * y * k * self.gamma_a / (theta**2 * q)
        return r, e * (1.0 - k * y) / q**3, dr_dth


class ZeroOrder:
    """Zero order in the reactant: ``R = exp[gamma (1 - 1/theta)]``.

    Parameters
    ----------
    gamma : float
        Arrhenius number.
    """

    def __init__(self, gamma):
        self.gamma = float(gamma)

    def __call__(self, y, theta):
        """Return ``(R, dR/dy, dR/dtheta)``."""
        r = np.exp(self.gamma * (1.0 - 1.0 / theta)) * np.ones_like(y)
        return r, np.zeros_like(y), r * self.gamma / theta**2


def violation(v, beta, y_ref=1.0, theta_ref=1.0):
    """``max |(theta - theta_ref) - beta (y_ref - y)| / |beta|`` over the cells.

    Parameters
    ----------
    v : numpy.ndarray
        State, shape ``(n, 2)``: ``y`` then ``theta``.
    beta : float
        Prater number.
    y_ref, theta_ref : float, optional
        Reference values; 1 and 1 is the bulk-referenced relation.

    Returns
    -------
    float
    """
    v = np.asarray(v).reshape(-1, 2)
    return float(np.max(np.abs((v[:, 1] - theta_ref) - beta * (y_ref - v[:, 0]))) / abs(beta))


# ------------------------------------------------------------------ the pellet
class TwoFieldPellet:
    """Species and energy as two coupled fields; the Prater relation is not used.

    ``div(grad y) = phi^2 R(y, theta)`` and ``div(grad theta) = -beta phi^2
    R(y, theta)`` on ``u`` in (0, 1), symmetric at the centre, with
    ``y = theta = 1`` at the surface or a film on either field.

    Parameters
    ----------
    geom : {"slab", "cylinder", "sphere"}, optional
        Geometry; sets the geometry index of ``construct_div`` and of
        :meth:`eta`.
    n : int, optional
        Number of cells.
    bi_m, bi_h : float or None, optional
        Biot numbers of a film on the species and the energy field; None is a
        Dirichlet surface. One of them alone is the page's boundary-type defect.
    nu : int, optional
        Overrides the geometry index of ``construct_div`` only (a defect the
        page injects); :meth:`eta` keeps the one of ``geom``.
    theta_floor : float, optional
        Lower bound on ``theta`` during Newton, a numerical safeguard: the
        Arrhenius factors overflow below about ``theta = 0.012``. By the Prater
        identity every root has ``theta >= 1 + beta`` for ``beta < 0``, so the
        floor cannot bind at a root of the page's sweep (``beta >= -0.4``); the
        solve checks that no converged state comes near it.

    Attributes
    ----------
    layout : pymrm.BlockLayout
        The fields ``"y"`` and ``"theta"`` interleaved in a state ``(n, 2)``.
    lap, lap_bc : scipy.sparse.csc_array, numpy.ndarray
        The constant operator ``div grad`` on both fields and its boundary
        contribution.
    """

    def __init__(self, geom="sphere", n=200, bi_m=None, bi_h=None, nu=None,
                 theta_floor=0.2):
        self.geom, self.n = geom, n
        self.nu_geom = NU[geom]
        nu_div = self.nu_geom if nu is None else nu
        self.layout = BlockLayout.from_fields(["y", "theta"], (n,))
        self.shape = self.layout.shape                     # (n, 2), fields last
        self.x_f = np.linspace(0.0, 1.0, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # bc: a dc/dn + b c = d with n the OUTWARD normal, one column per field.
        # centre u = 0: symmetry, dy/du = 0 and dtheta/du = 0
        bc_centre = {"a": [[1.0, 1.0]], "b": [[0.0, 0.0]], "d": [[0.0, 0.0]]}
        # surface u = 1, outward normal +u. Dirichlet: y = 1, theta = 1. A film
        # supplies the flux into the pellet:
        #   D_e dc/du = k_g (c_b - c_s)  ->  dy/du + Bi_m y = Bi_m
        #   lam dT/du = h (T_b - T_s)    ->  dtheta/du + Bi_h theta = Bi_h
        # with the same sign pattern for both fields (a = 1, b = d = Bi).
        a_s = [0.0 if bi_m is None else 1.0, 0.0 if bi_h is None else 1.0]
        b_s = [1.0 if bi_m is None else bi_m, 1.0 if bi_h is None else bi_h]
        bc_surface = {"a": [a_s], "b": [b_s], "d": [b_s]}
        self.bc = (bc_centre, bc_surface)
        grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc, axis=0)
        # nu: 0 slab (Cartesian), 1 cylinder, 2 sphere
        div_mat = construct_div(self.shape, self.x_f, nu=nu_div, axis=0)
        self.lap = (div_mat @ grad_mat).tocsc()
        self.lap_bc = (div_mat @ grad_bc).toarray().reshape(-1, 1)
        self.bounds = (np.array([-np.inf, theta_floor]), np.inf)   # per field

    def source(self, v, phi, beta, rate, heat_factor=1.0, energy_rate_isothermal=False):
        """Source terms and their pointwise Jacobian.

        Parameters
        ----------
        v : numpy.ndarray
            State ``(n, 2)``.
        phi, beta : float
            Thiele and Prater numbers.
        rate : callable
            ``rate(y, theta) -> (R, dR/dy, dR/dtheta)``.
        heat_factor : float, optional
            Multiplies ``beta`` in the energy source only; 1 is the model,
            ``-1`` and ``1.01`` are the page's sign and scale defects.
        energy_rate_isothermal : bool, optional
            Evaluate the energy equation's rate at ``theta = 1`` (a defect).

        Returns
        -------
        src : numpy.ndarray
            Shape ``(n, 2)``: ``-phi^2 R`` and ``beta phi^2 R``.
        dsrc : numpy.ndarray
            Shape ``(n, 2, 2)``, ``dsrc[i, a, b] = d src_a / d v_b``.
        """
        y, th = v[:, 0], v[:, 1]
        p2, hb = phi * phi, heat_factor * beta
        r, r_y, r_th = rate(y, th)
        if energy_rate_isothermal:
            r_e, r_e_y, _ = rate(y, np.ones_like(th))
            r_e_th = np.zeros_like(th)
        else:
            r_e, r_e_y, r_e_th = r, r_y, r_th
        src = np.stack([-p2 * r, hb * p2 * r_e], axis=-1)
        dsrc = np.empty((self.n, 2, 2))
        dsrc[:, 0, 0], dsrc[:, 0, 1] = -p2 * r_y, -p2 * r_th
        dsrc[:, 1, 0], dsrc[:, 1, 1] = hb * p2 * r_e_y, hb * p2 * r_e_th
        return src, dsrc

    def residual(self, v, phi, beta, rate, **defects):
        """``div(grad v) + src`` and its exact Jacobian, for :func:`pymrm.newton`."""
        v = np.asarray(v).reshape(self.shape)
        src, dsrc = self.source(v, phi, beta, rate, **defects)
        # the 2 x 2 block of each cell couples the fields at the same point
        jac_src = construct_coefficient_matrix(dsrc, shape=((self.n, 2, 1), (self.n, 1, 2)))
        g = self.lap @ v.reshape(-1, 1) + self.lap_bc + src.reshape(-1, 1)
        return g, self.lap + jac_src

    def solve(self, phi, beta, rate, check=True, **defects):
        """Steady state reached from the surface state ``y = theta = 1``.

        A pseudo-time march (:func:`pymrm.pseudo_transient`, equal pseudo
        capacities for both fields), finished by a steady Newton solve with
        ``tol=1e-12`` and ``stagnation=True``.

        Parameters
        ----------
        phi, beta : float
            Thiele and Prater numbers.
        rate : callable
            Rate law.
        check : bool, optional
            Raise ``RuntimeError`` unless the steady solve converged.
        **defects
            ``heat_factor``, ``energy_rate_isothermal`` (see :meth:`source`).

        Returns
        -------
        v : numpy.ndarray
            Shape ``(n, 2)``.
        result : scipy.optimize.OptimizeResult
            The :func:`pymrm.pseudo_transient` result.
        """
        res = pseudo_transient(
            lambda x: self.residual(x, phi, beta, rate, **defects), np.ones(self.shape),
            newton_kwargs={"maxfev": 10, "bounds": self.bounds},
            steady_kwargs={"tol": 1e-12, "stagnation": True, "maxfev": 100,
                           "bounds": self.bounds})
        ok = bool(res.success and res.steady is not None and res.steady.success)
        if check and not ok:
            raise RuntimeError(f"pellet solve failed (phi={phi}, beta={beta}): {res.message}")
        v = res.x.reshape(self.shape)
        if check and v[:, 1].min() < 1.5 * self.bounds[0][1]:
            raise RuntimeError("theta reached the neighbourhood of its floor")
        return v, res

    def newton_iterate(self, phi, beta, rate, maxfev):
        """The plain Newton iterate after ``maxfev`` steps from ``y = theta = 1``,
        converged or not (the page's non-converged defect)."""
        res = newton(lambda x: self.residual(x, phi, beta, rate), np.ones(self.shape),
                     tol=1e-12, maxfev=maxfev, bounds=self.bounds)
        return res.x.reshape(self.shape)

    def backward_error(self, v, phi, beta, rate, **defects):
        """Componentwise backward error of ``v``, as :func:`pymrm.newton` forms it.

        ``max_i |g_i| / (|J| (|v| + v_ref) + |J v - g|)_i`` with ``v_ref`` the
        largest ``|v|`` of each field, so that a cell far below its field's
        magnitude (a depleted centre) is judged against the field.
        """
        v = np.asarray(v, float).reshape(self.shape)
        g, jac = self.residual(v, phi, beta, rate, **defects)
        x = v.reshape(-1, 1)
        x_abs = (np.abs(v) + np.abs(v).max(axis=0)).reshape(-1, 1)
        denom = abs(jac) @ x_abs + np.abs(jac @ x - g)
        return float(np.max(np.abs(g) / denom))

    def surface(self, v):
        """Surface values ``(y_s, theta_s)`` and gradients ``(dy/du, dtheta/du)`` at u = 1."""
        val, grad = compute_boundary_values(np.asarray(v).reshape(self.shape), self.x_f,
                                            self.x_c, self.bc[1], axis=0, bound_id=1)
        return np.ravel(val), np.ravel(grad)

    def eta(self, v, phi, read="face"):
        """Effectiveness factor from the surface flux, ``(nu + 1) dy/du(1) / phi^2``.

        ``read="face"`` takes the surface gradient from
        :func:`pymrm.compute_boundary_values`; ``"cell"`` from the last cell
        centre and the surface value, a two-point gradient (a defect the page
        injects).
        """
        val, grad = self.surface(v)
        if read == "cell":
            v = np.asarray(v).reshape(self.shape)
            grad = [(val[0] - v[-1, 0]) / (self.x_f[-1] - self.x_c[-1])]
        return (self.nu_geom + 1.0) * float(grad[0]) / phi**2


# ------------------------------------------------------------------ transient
def _rhs(pellet, v, phi, beta, rate, le, heat_factor):
    """``dv/dtau`` and its Jacobian: ``L (div grad v + src)``, ``L = diag(1, Le)``."""
    g, jac = pellet.residual(v, phi, beta, rate, heat_factor=heat_factor)
    lewis = np.tile([1.0, le], pellet.n)
    return lewis * g.ravel(), diags_array(lewis) @ jac


def march(pellet, phi, beta, rate, le, t_end, per_decade, t_first=1e-5, heat_factor=1.0,
          order=2):
    """Start-up from ``y = theta = 1`` by variable-step BDF2, one Newton solve per step.

    ``dy/dtau = div grad y - phi^2 R`` and ``dtheta/dtau = Le (div grad theta +
    beta phi^2 R)``, time on the mass-diffusion scale. The first step is
    backward Euler to ``t_first``; then the steps grow geometrically with
    ``per_decade`` steps per decade up to ``t_end``.

    Parameters
    ----------
    pellet : TwoFieldPellet
        Operators and grid.
    phi, beta : float
        Thiele and Prater numbers.
    rate : callable
        Rate law.
    le : float
        Lewis number.
    t_end : float
        End of the window.
    per_decade : int
        Steps per decade of time.
    t_first : float, optional
        Length of the first step.
    heat_factor : float, optional
        See :meth:`TwoFieldPellet.source`.
    order : {2, 1}, optional
        2 is BDF2; 1 is backward Euler throughout (a defect the page injects).

    Returns
    -------
    t : numpy.ndarray
        Time levels, starting at 0.
    dev : numpy.ndarray
        :func:`violation` at each level.
    v : numpy.ndarray
        The state at ``t_end``.
    """
    n_steps = int(round(per_decade * np.log10(t_end / t_first)))
    t = np.concatenate([[0.0], t_first * (t_end / t_first) ** (np.arange(n_steps + 1) / n_steps)])
    eye = diags_array(np.ones(2 * pellet.n))
    v_old, v = None, np.ones(2 * pellet.n)
    dev = [violation(v, beta)]
    for k in range(1, len(t)):
        h = t[k] - t[k - 1]
        if k == 1 or order == 1:            # backward Euler
            a0, a1, a2 = 1.0, -1.0, 0.0
        else:                               # BDF2 for a step ratio w
            w = h / (t[k - 1] - t[k - 2])
            a0, a1, a2 = (1 + 2 * w) / (1 + w), -(1 + w), w * w / (1 + w)
        v_n, v_p = v.copy(), (v if v_old is None else v_old).copy()

        def step(x):
            f, jac = _rhs(pellet, x, phi, beta, rate, le, heat_factor)
            return (a0 * x.ravel() + a1 * v_n + a2 * v_p) / h - f, (a0 / h) * eye - jac

        res = newton(step, v, tol=1e-12, stagnation=True, maxfev=30)
        if not res.success:
            raise RuntimeError(f"time step {k} failed: {res.message}")
        v_old, v = v, res.x.ravel()
        dev.append(violation(v, beta))
    return t, np.array(dev), v.reshape(pellet.shape)


def peak_of(t, dev):
    """Maximum of ``dev(t)``: a parabola in ``log t`` through the largest level
    and its two neighbours. Returns ``(peak, t_peak)``."""
    i = int(np.argmax(dev))
    if not 1 < i < len(t) - 1:
        raise ValueError("maximum at the window edge")
    x = np.log(t[i - 1:i + 2])
    c = np.polyfit(x, dev[i - 1:i + 2], 2)
    x_m = -c[1] / (2.0 * c[0])
    return float(np.polyval(c, x_m)), float(np.exp(x_m))


def march_ivp(pellet, phi, beta, rate, le, t_end, rtol=1e-9):
    """The same semi-discrete start-up by ``solve_ivp`` (BDF, error-controlled),
    with the peak of :func:`violation` root-found on the dense output.

    Returns
    -------
    peak, t_peak : float
    """
    def f(t, w):
        return _rhs(pellet, w, phi, beta, rate, le, 1.0)[0]

    def jac(t, w):
        return _rhs(pellet, w, phi, beta, rate, le, 1.0)[1]

    sol = solve_ivp(f, (0.0, t_end), np.ones(2 * pellet.n), method="BDF", jac=jac,
                    rtol=rtol, atol=1e-2 * rtol, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    tt = np.geomspace(1e-5, t_end, 400)
    d = np.array([violation(sol.sol(x), beta) for x in tt])
    i = int(np.argmax(d))
    res = minimize_scalar(lambda lx: -violation(sol.sol(np.exp(lx)), beta),
                          bracket=(np.log(tt[i - 1]), np.log(tt[i]), np.log(tt[i + 1])),
                          method="brent", options={"xtol": 1e-10})
    return float(-res.fun), float(np.exp(res.x))


# ------------------------------------------------------------- reduced equation
class ReducedPellet:
    """``div(grad y) = phi^2 y exp[gamma beta (1 - y)/(1 + beta (1 - y))]``, ``y(1) = 1``.

    The single equation of ``B1.1``: the Prater relation already substituted.

    Parameters
    ----------
    geom : {"slab", "cylinder", "sphere"}, optional
        Geometry.
    n : int, optional
        Number of cells.
    """

    def __init__(self, geom="sphere", n=400):
        self.nu, self.n = NU[geom], n
        self.shape = (n, 1)                                 # one field: (n, 1), never (n,)
        self.x_f = np.linspace(0.0, 1.0, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0},           # centre: dy/du = 0
                   {"a": 0.0, "b": 1.0, "d": 1.0})           # surface: y = 1
        grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc, axis=0)
        div_mat = construct_div(self.shape, self.x_f, nu=self.nu, axis=0)  # 0/1/2 geometry
        self.lap = (div_mat @ grad_mat).tocsc()
        self.lap_bc = (div_mat @ grad_bc).toarray().reshape(-1, 1)

    def residual(self, y, phi, beta, gamma):
        """Residual and exact Jacobian."""
        y = np.asarray(y).reshape(self.shape)
        s = beta * (1.0 - y)
        e = np.exp(gamma * s / (1.0 + s))
        r = y * e
        dr = e * (1.0 - y * gamma * beta / (1.0 + s) ** 2)
        g = self.lap @ y + self.lap_bc - phi**2 * r
        return g, self.lap - construct_coefficient_matrix(phi**2 * dr, shape=self.shape)

    def solve(self, phi, beta, gamma, y_init):
        """Newton from ``y_init`` (``tol=1e-12``, ``stagnation=True``); raises on failure."""
        res = newton(lambda y: self.residual(y, phi, beta, gamma),
                     np.asarray(y_init, float).reshape(self.shape), tol=1e-12,
                     stagnation=True, maxfev=60)
        if not res.success:
            raise RuntimeError(f"reduced solve failed: {res.message}")
        return res.x.reshape(self.shape)


# -------------------------------------------------------- independent routes
def collocation(nu, phi, beta, gamma, tol=1e-8, u0=1e-4, n_ramp=8, max_nodes=200000,
                heat_factor=1.0):
    """The two-field first-order problem by ``solve_bvp``, no pymrm.

    The domain starts at ``u0`` with the regular-singular series condition
    ``y'(u0) = phi^2 R u0/(nu + 1)`` (and ``theta' = -beta`` times it), and
    ``beta`` is ramped in ``n_ramp`` equal steps inside this solver.
    ``heat_factor`` multiplies ``beta`` in the energy equation only (a defect
    the page injects).

    Returns
    -------
    scipy.integrate._bvp.BVPResult
    """
    k = u0 / (nu + 1.0)
    u_g = np.linspace(u0, 1.0, 300)
    w = np.vstack([np.ones_like(u_g), np.zeros_like(u_g), np.ones_like(u_g),
                   np.zeros_like(u_g)])
    sol = None
    for b in beta * np.linspace(0.0, 1.0, n_ramp + 1)[1:]:
        def rhs(x, ww, b=b):
            yy, dy, tt, dt = ww
            r = yy * np.exp(gamma * (1.0 - 1.0 / tt))
            return np.vstack([dy, phi**2 * r - nu * dy / x,
                              dt, -heat_factor * b * phi**2 * r - nu * dt / x])

        def bcs(wa, wb, b=b):
            r = wa[0] * np.exp(gamma * (1.0 - 1.0 / wa[2]))
            return np.array([wa[1] - k * phi**2 * r, wa[3] + k * heat_factor * b * phi**2 * r,
                             wb[0] - 1.0, wb[2] - 1.0])

        sol = solve_bvp(rhs, bcs, u_g, w, tol=tol, max_nodes=max_nodes)
        u_g, w = sol.x, sol.y
    return sol


def shooting(beta, gamma, y_c, nu=2, rtol=1e-11):
    """The reduced equation integrated outward from the centre value ``y_c``.

    With ``s = phi u`` the equation does not contain ``phi``: the ``s`` at
    which ``y`` reaches 1 is ``phi``.

    Returns
    -------
    tuple of float or None
        ``(phi, eta)``, or None when ``y`` does not reach 1 by ``s = 400``.
    """
    def rate(y):
        s = beta * (1.0 - y)
        return y * np.exp(gamma * s / (1.0 + s))

    def rhs(s, v):
        return [v[1], rate(v[0]) - (nu / s) * v[1]]

    def hit(s, v):
        return v[0] - 1.0
    hit.terminal, hit.direction = True, 1
    s0 = 1e-6                     # series start: y = y_c + R s^2 / (2 (nu + 1))
    v0 = [y_c + rate(y_c) * s0**2 / (2 * (nu + 1)), rate(y_c) * s0 / (nu + 1)]
    sol = solve_ivp(rhs, (s0, 400.0), v0, events=hit, rtol=rtol, atol=1e-3 * rtol)
    if not sol.t_events[0].size:
        return None
    s_star = float(sol.t_events[0][0])
    return s_star, (nu + 1.0) * float(sol.y_events[0][0][1]) / s_star


def shooting_at_phi(beta, gamma, phi, y_guess, nu=2, rtol=1e-11):
    """``eta`` at Thiele modulus ``phi`` on the branch through the centre value
    ``y_guess``, by root-finding the centre value of :func:`shooting`.

    The bracket ``(y_guess/f, y_guess f)`` starts at ``f = 1.1`` and widens
    (``f -> f^2``) until :func:`shooting` brackets ``phi``; on a branch with
    several steady states this keeps the root on the branch of ``y_guess``.

    Returns
    -------
    eta, y_c : float
    """
    def miss(y_c):
        return shooting(beta, gamma, y_c, nu, rtol)[0] - phi

    f = 1.1
    for _ in range(6):
        lo, hi = y_guess / f, min(y_guess * f, 1.0 - 1e-12)
        if miss(lo) * miss(hi) < 0.0:
            break
        f *= f
    else:
        raise RuntimeError("no bracket around the centre value")
    y_c = brentq(miss, lo, hi, xtol=1e-300, rtol=1e-14)
    return shooting(beta, gamma, y_c, nu, rtol)[1], y_c


def prater_beta(dh_kj_per_mol, d_e, c_s, lam, t_s):
    """``beta = (-dH) D_e c_s/(lam T_s)``; ``dH`` in kJ/mol, negative when exothermic.

    Parameters
    ----------
    dh_kj_per_mol : float
        Reaction enthalpy [kJ/mol].
    d_e : float
        Effective diffusivity [m2/s].
    c_s : float
        Surface concentration [mol/m3].
    lam : float
        Effective conductivity [W/m/K].
    t_s : float
        Surface temperature [K].

    Returns
    -------
    float
    """
    return (-dh_kj_per_mol * 1e3) * d_e * c_s / (lam * t_s)
