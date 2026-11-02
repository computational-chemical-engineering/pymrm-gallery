"""The Szekely-Evans grain model: exact penetration law and the pymrm march.

Model code of gallery page B3.2 (Szekely and Evans 1970, Part I: a
semi-infinite porous solid of spherical shrinking-core grains, quasi-steady gas
phase, first-order kinetics). Dimensionless variables: core radius
``xi = R'/R_g``, gas concentration ``psi = C'_p/C_0``, depth ``eta = y/ell``,
time ``theta = t/tau`` and the one group ``g = k R_g/D``.

Contents
--------
G_of_xi, xi_of_zeta, X_of_theta, E_exact
    The grain law ``G(xi) = zeta`` (zeta the local exposure, the time
    integral of psi), its inverse, ``X(theta)`` and the exact equivalent
    penetration ``E = sqrt(6 X(theta))``.
lag, profile_exact, zone_width, departure_thetas
    The late-time lag ``1/4 + g/30``, the exact profiles at any set of depths,
    the width of the reaction zone, and where the kinetic and the
    diffusion-limited laws stop being within a given fraction of E.
yagi_kunii_eq6
    Eq. 6 of Yagi and Kunii (page B3.1), the shrinking core with a film.
scales
    ``ell``, ``tau`` and ``g`` from dimensional parameters.
GrainModel
    The gas and grain fields marched in time with pymrm: a finite-volume gas
    equation, the grain equation by the trapezoidal rule (or another theta
    method), the grain exhaustion as a complementarity condition, and one
    ``newton`` solve per step with the exact Jacobian.

Several functions take an argument whose default is the correct model and
whose alternative is a defect the page's break table injects.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.sparse import diags_array
from scipy.sparse.linalg import spsolve

from pymrm import (BlockLayout, NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_div, construct_grad, newton)

__all__ = ["G_of_xi", "xi_of_zeta", "X_of_theta", "E_exact", "lag", "profile_exact",
           "zone_width", "departure_thetas", "yagi_kunii_eq6", "scales", "GrainModel"]


# ------------------------------------------------------------------ exact route
# The grain law is written in u = 1 - xi, the converted fraction of the grain
# radius, so that small exposures carry no cancellation:
#   G = u + g u^2/2 - g u^3/3,   x = 1 - xi^3 = u (3 - 3u + u^2),
#   X = Q(u) = integral_0^u x(u') G'(u') du'.
def _G_u(u, g):
    return u + g * u**2 / 2 - g * u**3 / 3


def _dG_u(u, g):
    return 1.0 + g * u * (1.0 - u)


def _Q(u, g):
    return (1.5 * u**2 + (g - 1.0) * u**3 + (1.0 - 6.0 * g) / 4.0 * u**4 + 0.8 * g * u**5
            - g / 6.0 * u**6)


def G_of_xi(xi, g):
    """Grain law as the paper's equations give it, ``G(xi) = zeta``.

    Parameters
    ----------
    xi : array_like
        Core radius over grain radius, in ``[0, 1]``.
    g : float
        ``k R_g / D``.

    Returns
    -------
    numpy.ndarray
        ``(1 - xi) + (g/6)(1 - 3 xi^2 + 2 xi^3)``, the exposure at which a
        grain has shrunk to ``xi``.
    """
    xi = np.asarray(xi, float)
    return (1.0 - xi) + (g / 6.0) * (1.0 - 3.0 * xi**2 + 2.0 * xi**3)


def _u_of_zeta(zeta, g):
    """Invert ``G`` in ``u`` (monotone on [0, 1]): bisection, then Newton."""
    z = np.clip(np.asarray(zeta, float), 0.0, 1.0 + g / 6.0)
    lo, hi = np.zeros_like(z), np.ones_like(z)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        up = _G_u(mid, g) < z
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    u = 0.5 * (lo + hi)
    for _ in range(2):
        u = np.clip(u - (_G_u(u, g) - z) / _dG_u(u, g), 0.0, 1.0)
    return u


def xi_of_zeta(zeta, g):
    """Core radius of a grain that has seen the exposure ``zeta``.

    Parameters
    ----------
    zeta : array_like
        Local exposure, the time integral of ``psi``.
    g : float
        ``k R_g / D``.

    Returns
    -------
    numpy.ndarray
        ``xi`` with ``G(xi) = zeta``; 1 for ``zeta <= 0`` and 0 for
        ``zeta >= 1 + g/6`` (the grain is consumed).
    """
    return 1.0 - _u_of_zeta(zeta, g)


def X_of_theta(theta, g):
    """``X(theta)``, the integral of the local conversion ``x = 1 - xi^3`` over exposure.

    Parameters
    ----------
    theta : array_like
        Exposure (at the surface: the time).
    g : float
        ``k R_g / D``.

    Returns
    -------
    numpy.ndarray
        Closed form; linear in ``theta`` beyond ``theta* = 1 + g/6``.
    """
    th = np.asarray(theta, float)
    ts = 1.0 + g / 6.0
    return np.where(th >= ts, _Q(1.0, g) + (th - ts), _Q(_u_of_zeta(np.minimum(th, ts), g), g))


def E_exact(theta, g):
    """Exact equivalent penetration over ``ell``, ``E = sqrt(6 X(theta))``.

    Parameters
    ----------
    theta : array_like
        Time ``k C_0 t / (rho R_g)``.
    g : float
        ``k R_g / D``.

    Returns
    -------
    numpy.ndarray
    """
    return np.sqrt(6.0 * X_of_theta(theta, g))


def lag(g):
    """Late-time lag ``Delta``: ``E^2 = 6 (theta - Delta)`` exactly for ``theta >= 1 + g/6``.

    Parameters
    ----------
    g : float

    Returns
    -------
    float
        ``theta* - X(theta*) = 1/4 + g/30``.
    """
    return 1.0 + g / 6.0 - _Q(1.0, g)


def profile_exact(eta, theta, g, rtol=1e-12):
    """Exact gas and grain profiles at depths ``eta`` and time ``theta``.

    The exposure obeys ``zeta' = -sqrt(2 X(zeta)/3)`` from ``zeta(0) = theta``.
    In the consumed layer (``zeta >= theta*``) that is integrated in closed form,
    beyond it as an ODE for ``u = 1 - xi`` evaluated at the requested depths.

    Parameters
    ----------
    eta : array_like
        Increasing depths, for example the cell centres of a grid.
    theta, g : float
        Time and ``k R_g / D``.
    rtol : float, optional
        Relative tolerance of the ODE integration (DOP853).

    Returns
    -------
    psi, xi : numpy.ndarray
        ``psi = sqrt(X(zeta)/X(theta))`` and ``xi`` at ``eta``.
    """
    eta = np.asarray(eta, float)
    ts, Q1 = 1.0 + g / 6.0, _Q(1.0, g)
    Xt = float(X_of_theta(theta, g))
    a = Q1 - ts                                       # X = a + zeta for zeta >= ts
    eta_f = max(np.sqrt(6.0 * (a + theta)) - np.sqrt(6.0 * Q1), 0.0) if theta > ts else 0.0
    psi, xi = np.empty_like(eta), np.empty_like(eta)
    burnt = (eta <= eta_f) if theta > ts else np.zeros(eta.shape, bool)
    if burnt.any():
        s = np.sqrt(a + theta) - eta[burnt] / np.sqrt(6.0)
        psi[burnt], xi[burnt] = np.sqrt(s**2 / Xt), 0.0
    rest = ~burnt
    if rest.any():
        u0 = 1.0 if theta > ts else float(_u_of_zeta(theta, g))
        rhs = lambda e, u: -np.sqrt(2.0 * np.maximum(_Q(u, g), 0.0) / 3.0) / _dG_u(u, g)
        sol = solve_ivp(rhs, (eta_f, eta[rest][-1]), [u0], method="DOP853", t_eval=eta[rest],
                        rtol=rtol, atol=1e-300)
        u = sol.y[0]
        psi[rest], xi[rest] = np.sqrt(_Q(u, g) / Xt), 1.0 - u
    return psi, xi


def zone_width(g, x_lo=0.05, x_hi=0.95, xi_exponent=3):
    """Width of the reaction zone between local conversions ``x_lo`` and ``x_hi``.

    The exposure profile is a travelling wave, so the width is exactly constant
    once the conversion ``x_hi`` has appeared at the surface:
    ``w = integral dzeta / sqrt(2 X(zeta)/3)`` between the two exposures.

    Parameters
    ----------
    g : float
        ``k R_g / D``.
    x_lo, x_hi : float, optional
        Conversion levels ``x = 1 - xi^3``.
    xi_exponent : int, optional
        3 is the local conversion; 1 reads the levels on ``1 - xi`` instead (a
        defect: an observable definition slip).

    Returns
    -------
    float
        Width in units of ``ell``.
    """
    u_lo, u_hi = (1.0 - (1.0 - x) ** (1.0 / xi_exponent) for x in (x_lo, x_hi))
    f = lambda u: _dG_u(u, g) / np.sqrt(2.0 * _Q(u, g) / 3.0)
    return quad(f, u_lo, u_hi, epsabs=0.0, epsrel=1e-13, limit=200)[0]


def departure_thetas(g, frac=0.95):
    """Times at which the two classical laws leave a band around the exact E.

    Parameters
    ----------
    g : float
        ``k R_g / D``.
    frac : float, optional
        ``E_kinetic`` within ``E/frac`` and ``E_diffusion >= frac E``.

    Returns
    -------
    theta_early, theta_late : float
        The kinetic law ``3 theta`` holds to ``frac`` for ``theta < theta_early``,
        the diffusion law ``sqrt(6 (theta - Delta))`` for ``theta > theta_late``.
    """
    D = lag(g)
    t_e = brentq(lambda t: float(E_exact(t, g)) / (3.0 * t) - frac, 1e-6, 1e5, xtol=1e-14, rtol=1e-14)
    t_l = brentq(lambda t: np.sqrt(6.0 * (t - D)) / float(E_exact(t, g)) - frac,
                 D * (1.0 + 1e-12), 1.0 + g / 6.0, xtol=1e-14, rtol=1e-14)
    return t_e, t_l


def yagi_kunii_eq6(x, omega, gamma, cubic="printed"):
    """Eq. 6 of Yagi and Kunii as page B3.1 prints it: ``theta/theta_B`` at ``x = r/R``.

    Parameters
    ----------
    x : array_like
        Core radius over particle radius.
    omega, gamma : float
        ``k_c1/k_d1`` and ``3 k_f1/k_d1``; a large ``gamma`` is a negligible film.
    cubic : {"printed", "misread"}, optional
        ``omega (1 - 2 gamma)`` as printed, or ``omega (1 - gamma)`` (a defect).

    Returns
    -------
    numpy.ndarray
    """
    x = np.asarray(x, float)
    c3 = omega * (1.0 - 2.0 * gamma) if cubic == "printed" else omega * (1.0 - gamma)
    num = gamma * x + 3.0 * omega * gamma * x**2 + c3 * x**3
    return 1.0 - num / (omega + omega * gamma + gamma)


def scales(Dp, Rg, D, k, Sg, FT, rho, C0, factor=3.0):
    """Effective diffusivity and the reduction scales of the dimensionless model.

    Parameters
    ----------
    Dp, Rg, D, k, Sg, FT, rho, C0 : float
        Pore diffusivity, grain radius, product-layer diffusivity, rate
        constant, porosity, tortuosity, molar density of the solid and bulk
        gas concentration (cgs).
    factor : float, optional
        The 3 in ``ell``; another value is a defect.

    Returns
    -------
    dict
        ``Deff = Sg Dp / FT``, ``ell = sqrt(Deff Rg / (3 (1 - Sg) k))``,
        ``tau = rho Rg / (k C0)``, ``g = k Rg / D``.
    """
    Deff = Sg * Dp / FT
    return {"Deff": Deff, "ell": np.sqrt(Deff * Rg / (factor * (1.0 - Sg) * k)),
            "tau": rho * Rg / (k * C0), "g": k * Rg / D}


# ------------------------------------------------------------------ pymrm march
class GrainModel:
    """Szekely-Evans grain model on a truncated semi-infinite slab.

    ``diff psi'' = sink_pref xi^2 psi / (1 + g xi (1 - xi))`` (Eq. 27) and
    ``dxi/dt = -rate_pref psi / (1 + g xi (1 - xi))`` with ``xi >= 0`` (Eq. 21).
    All prefactors are 1 in dimensionless form; in cgs ``diff = D'``,
    ``sink_pref = 3 (1 - S_g) k / R_g``, ``rate_pref = k / (rho R_g)`` and
    ``c_surf = C_0``.

    Parameters
    ----------
    n : int, optional
        Number of cells on ``[0, L]``.
    L : float, optional
        Domain length; the far face stands for depth infinity.
    g : float, optional
        ``k R_g / D``.
    diff, sink_pref, rate_pref, c_surf : float, optional
        As above.
    nu : int, optional
        Geometry index of ``construct_div``: 0, the slab, is the model.
    far_bc : {"dirichlet", "neumann"}, optional
        ``psi = 0`` (the model) or zero flux at ``L``.
    area_exponent : int, optional
        Exponent of ``xi`` in the sink, the reacting area of a grain: 2.
    shield_sink, shield_rate : bool, optional
        Keep the ``(1 - xi)`` of the shell resistance in Eq. 27 and in Eq. 21.
    theta_scheme : float, optional
        Weight of the new time level in the grain equation: 0.5 is the
        trapezoidal rule, 1 backward Euler.
    jacobian : {"exact", "numjac"}, optional
        The analytic pointwise Jacobian, or ``NumJac((n, 2))`` of the source
        terms (the validation route).
    newton_tol : float, optional
        Tolerance of each ``newton`` solve.
    jac_defect : {None, "drop_dsink_dxi"}, optional
        A deliberate error in the analytic Jacobian: drop the derivative of the
        sink with respect to ``xi``. The residual is untouched, so the converged
        state is the same and only the convergence of ``newton`` suffers.

    Attributes
    ----------
    layout : pymrm.BlockLayout
        Fields ``"psi"`` and ``"xi"`` interleaved in a state ``(n, 2)``.
    lin, b : scipy.sparse.csc_array, numpy.ndarray
        The constant operator ``-div(diff grad)`` placed on the psi rows and its
        boundary contribution.
    """

    def __init__(self, n=600, L=35.0, g=0.5, diff=1.0, sink_pref=1.0, rate_pref=1.0, c_surf=1.0,
                 nu=0, far_bc="dirichlet", area_exponent=2, shield_sink=True, shield_rate=True,
                 theta_scheme=0.5, jacobian="exact", newton_tol=1e-12, jac_defect=None):
        self.n, self.g = n, g
        self.diff, self.sink_pref, self.rate_pref = diff, sink_pref, rate_pref
        self.area_exponent, self.shield_sink, self.shield_rate = area_exponent, shield_sink, shield_rate
        self.theta_scheme, self.jacobian, self.newton_tol = theta_scheme, jacobian, newton_tol
        self.jac_defect = jac_defect
        self.layout = BlockLayout.from_fields(["psi", "xi"], (n,))
        self.shape = self.layout.shape                               # (n, 2), fields last
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.dx = np.diff(self.x_f)
        # bc: a dpsi/dn + b psi = d with n the OUTWARD normal.
        #   eta = 0: psi = c_surf (Eq. 28)
        #   eta = L: psi = 0 (Eq. 29, psi -> 0 at depth infinity, truncated);
        #            the defect "neumann" is dpsi/dn = 0
        far = {"a": 1.0, "b": 0.0, "d": 0.0} if far_bc == "neumann" else {"a": 0.0, "b": 1.0, "d": 0.0}
        self.bc = ({"a": 0.0, "b": 1.0, "d": c_surf}, far)
        grad, grad_bc = construct_grad((n, 1), self.x_f, self.x_c, self.bc)
        # nu = 0: Cartesian slab, the semi-infinite sample; the grains' sphericity
        # enters through the sink and rate closures, not through this operator
        div = construct_div((n, 1), self.x_f, nu=nu)
        self.grad, self.grad_bc = grad, grad_bc
        self.lop = (div @ (-diff * grad)).tocsc()                   # -(diff psi')' on (n, 1)
        self.lop_bc = (div @ (-diff * grad_bc)).toarray().ravel()
        self.lin = self.layout.place(self.lop, "psi", "psi").tocsc()
        self.b = np.zeros(self.layout.size)
        self.layout.put(self.b, "psi", self.lop_bc)
        self.b = self.b.reshape(-1, 1)
        self.numjac = NumJac(self.shape) if jacobian == "numjac" else None
        self.t = 0.0
        self.xi = np.ones(n)
        self.psi = spsolve((self.lop + diags_array(self.sink_coeff(self.xi)[0])).tocsc(), -self.lop_bc)
        self.I_gas = 0.0
        self.q_prev = self.influx()
        self.hist = {"t": [], "EP": [], "xi_s": [], "xi_0": [], "front": []}
        self.snapshots = {}
        self.newton_iters = []

    # ---------------------------------------------------------------- closures
    @staticmethod
    def _shell(xi, g, shield):
        """Shell-plus-reaction resistance ``1 + g xi (1 - xi)`` and its xi-derivative."""
        xp = np.clip(xi, 0.0, 1.0)
        inside = (xi > 0.0) & (xi < 1.0)
        if shield:
            return 1.0 + g * xp * (1.0 - xp), np.where(inside, g * (1.0 - 2.0 * xp), 0.0)
        return 1.0 + g * xp, np.where(inside, g, 0.0)

    def sink_coeff(self, xi):
        """Eq. 27 sink per unit psi and its xi-derivative."""
        xp = np.clip(xi, 0.0, 1.0)
        inside = (xi > 0.0) & (xi < 1.0)
        p = self.area_exponent
        den, dden = self._shell(xi, self.g, self.shield_sink)
        area = xp**p
        darea = np.where(inside, p * xp ** (p - 1), 0.0)
        return self.sink_pref * area / den, self.sink_pref * (darea * den - area * dden) / den**2

    def rate(self, psi, xi):
        """``-dxi/dt`` (Eq. 21) and its derivatives in psi and xi, unclamped."""
        den, dden = self._shell(xi, self.g, self.shield_rate)
        r = self.rate_pref * psi / den
        return r, self.rate_pref / den, -r * dden / den

    # ---------------------------------------------------------------- residual
    def source(self, v):
        """Pointwise terms of the residual and their Jacobian blocks.

        Parameters
        ----------
        v : numpy.ndarray
            State ``(n, 2)``.

        Returns
        -------
        src : numpy.ndarray
            ``(n, 2)``: the sink ``s psi`` and the grain equation
            ``min(xi - xi_old + dt [w r + (1 - w) r_old], xi)``, the complementarity
            form of Eq. 21's ``0 <= R'``.
        dsrc : numpy.ndarray
            ``(n, 2, 2)``, ``dsrc[i, a, b] = d src_a / d v_b``.
        """
        psi, xi = v[:, 0], v[:, 1]
        w, dt = self.theta_scheme, self.dt
        s, ds = self.sink_coeff(xi)
        r, r_psi, r_xi = self.rate(psi, xi)
        step = xi - self.xi_old + dt * (w * r + (1.0 - w) * self.r_old)
        on = step <= xi                                  # the grain equation holds
        src = np.stack([s * psi, np.where(on, step, xi)], axis=-1)
        dsrc = np.empty((self.n, 2, 2))
        dsrc[:, 0, 0], dsrc[:, 0, 1] = s, (0.0 if self.jac_defect == "drop_dsink_dxi" else ds * psi)
        dsrc[:, 1, 0] = np.where(on, dt * w * r_psi, 0.0)
        dsrc[:, 1, 1] = np.where(on, 1.0 + dt * w * r_xi, 1.0)
        return src, dsrc

    def residual(self, v):
        """``-div(diff grad psi) + sink`` and the grain equation, with the Jacobian."""
        v = np.asarray(v).reshape(self.shape)
        if self.numjac is not None:
            src, jac_src = self.numjac(lambda u: self.source(u)[0], v)
        else:
            src, dsrc = self.source(v)
            # the 2 x 2 block of each cell couples the fields at the same point
            jac_src = construct_coefficient_matrix(dsrc, shape=((self.n, 2, 1), (self.n, 1, 2)))
        g = self.lin @ v.reshape(-1, 1) + self.b + src.reshape(-1, 1)
        return g, self.lin + jac_src

    def jacobian_gap(self, dt):
        """Relative gap between the analytic and the NumJac Jacobian at the current state.

        The state is perturbed off the converged values so that every entry of the
        pointwise block is non-trivial; the step ``dt`` and the old level are those
        of a step starting from the current state.

        Returns
        -------
        float
            ``max|J_exact - J_numjac| / max|J_exact|`` of the pointwise source terms.
        """
        self.dt = dt
        self.xi_old = self.xi.copy()
        self.r_old = np.where(self.xi_old > 0.0, self.rate(self.psi, self.xi_old)[0], 0.0)
        v = np.column_stack([self.psi, self.xi * 0.999 + 1e-4])
        _, dsrc = self.source(v)
        j_exact = construct_coefficient_matrix(dsrc, shape=((self.n, 2, 1), (self.n, 1, 2))).toarray()
        _, j_num = NumJac(self.shape)(lambda u: self.source(u)[0], v)
        return float(np.abs(j_exact - j_num.toarray()).max() / np.abs(j_exact).max())

    # ---------------------------------------------------------------- march
    def influx(self):
        """Gas flux into the solid at ``eta = 0``, from the gradient the solve uses."""
        _, grad0 = compute_boundary_values(self.psi.reshape(-1, 1), self.x_f, self.x_c, self.bc[0],
                                           bound_id=0)
        return float(-self.diff * np.ravel(grad0)[0])

    def march(self, t_end, dt, snap_at=()):
        """March to ``t_end`` in steps ``dt``, recording the history at every step.

        Parameters
        ----------
        t_end, dt : float
            End time and step; ``t_end/dt`` is rounded to an integer step count.
        snap_at : sequence of float, optional
            Step times at which to store ``(psi, xi)``.

        Returns
        -------
        GrainModel
            ``self``; ``hist`` holds ``t``, ``EP``, ``xi_s`` (surface
            ``xi``, extrapolated), ``xi_0`` (first cell) and ``front`` per step.

        Raises
        ------
        RuntimeError
            If a ``newton`` solve does not converge.
        """
        self.dt = dt
        t0 = self.t
        snap = sorted(snap_at)
        v = np.column_stack([self.psi, self.xi])
        for i in range(int(round(t_end / dt))):
            self.xi_old = self.xi.copy()
            r_old = self.rate(self.psi, self.xi_old)[0]
            self.r_old = np.where(self.xi_old > 0.0, r_old, 0.0)
            res = newton(self.residual, v, tol=self.newton_tol, maxfev=30, damping="monotone",
                         stagnation=True)
            if not res.success:
                raise RuntimeError(f"newton failed at step {i} (t = {self.t:.4g}): {res.message}")
            self.newton_iters.append(int(res.nit))
            v = np.asarray(res.x).reshape(self.shape).copy()
            v[:, 1] = np.maximum(v[:, 1], 0.0)
            self.psi, self.xi = v[:, 0].copy(), v[:, 1].copy()
            self.t = t0 + (i + 1) * dt
            q = self.influx()
            self.I_gas += 0.5 * dt * (self.q_prev + q)                # trapezoid in time
            self.q_prev = q
            for key, val in (("t", self.t), ("EP", self.EP()), ("xi_s", self.xi_surface()),
                             ("xi_0", self.xi_surface(read="cell")), ("front", self.front())):
                self.hist[key].append(val)
            while snap and self.t >= snap[0] - 0.5 * dt:
                self.snapshots[snap.pop(0)] = (self.psi.copy(), self.xi.copy())
        return self

    def history(self):
        """The recorded history as arrays, keyed as ``hist``."""
        return {k: np.asarray(v) for k, v in self.hist.items()}

    # ---------------------------------------------------------------- observables
    def EP(self):
        """Eq. 34, the equivalent penetration ``sum (1 - xi^3) dx``."""
        return float(np.sum((1.0 - np.clip(self.xi, 0.0, 1.0) ** 3) * self.dx))

    def xi_surface(self, read="extrapolated"):
        """``xi`` at ``eta = 0``.

        Parameters
        ----------
        read : {"extrapolated", "cell"}, optional
            Linear extrapolation from the first two centres (second order), or
            the first cell's value (first order, a defect).
        """
        if read == "cell":
            return float(self.xi[0])
        x0, x1 = self.x_c[0], self.x_c[1]
        return float(self.xi[0] + (self.xi[0] - self.xi[1]) * x0 / (x1 - x0))

    def front(self):
        """Depth of the fully consumed layer: the face before the first cell with ``xi > 0``."""
        alive = np.nonzero(self.xi > 0.0)[0]
        return float(self.x_f[alive[0]] if alive.size else self.x_f[-1])

    def solid_consumed(self):
        """Solid consumed per unit face area, in units of the gas uptake."""
        return self.EP() * self.sink_pref / (3.0 * self.rate_pref)

    def balance_gap(self):
        """``(influx - outflux) - sum(sink psi dx)`` at the current state."""
        f = -self.diff * (self.grad @ self.psi.reshape(-1, 1) + self.grad_bc.toarray()).ravel()
        return abs((f[0] - f[-1]) - float(np.sum(self.sink_coeff(self.xi)[0] * self.psi * self.dx)))
