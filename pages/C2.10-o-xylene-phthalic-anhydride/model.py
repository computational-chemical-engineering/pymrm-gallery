"""Fixed-bed tube models of Froment (1967): 2-D pseudo-homogeneous and 1-D.

Model code of gallery page C2.10 (o-xylene to phthalic anhydride). The tube
classes know nothing about the chemistry: they take a pointwise source function
of the state, so another fixed-bed problem (other kinetics, other fields) uses
them unchanged. The o-xylene chemistry and the printed constants are separate.

Contents
--------
FromentCase
    The printed constants of the o-xylene case and the dimensionless groups of
    the paper's Eq. 1 derived from them.
froment_source
    The triangular scheme A -> B -> C, A -> C as a pointwise source of
    ``(x, w, t)``: conversions to B and to C, temperature rise.
u_overall
    Overall heat-transfer coefficient of the 1-D model from the paper's
    moment-matching formula.
Tube2D
    Two-dimensional pseudo-homogeneous tube: radial dispersion of every field,
    a Robin wall, marched along the axis with a stiff integrator.
Tube1D
    The radially lumped tube, marched along the axis.
Tube1DSteady
    The same 1-D tube as a steady boundary-value problem: first-order upwind
    convection and ``pymrm.newton``.
critical_inlet
    Bisection of a runaway predicate in the inlet temperature.
exchanger_mean
    The paper's Bessel series for the reactionless packed-bed heat exchanger
    (radial-mean temperature), an analytic check of the radial operator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import j0, j1

from pymrm import (NumJac, construct_coefficient_matrix, construct_convflux_upwind,
                   construct_div, construct_grad, newton)

__all__ = ["FromentCase", "froment_source", "u_overall", "Tube2D", "Tube1D",
           "Tube1DSteady", "critical_inlet", "exchanger_mean"]


@dataclass(frozen=True)
class FromentCase:
    """Printed constants of Froment's o-xylene case (units kg, m, hr, kcal, degC).

    Parameters
    ----------
    NA0, N0 : float
        Inlet mole fractions of o-xylene and oxygen.
    dH1, dH3 : float
        Heats of the partial oxidation and of the complete combustion,
        kcal/kmol, positive for exothermic (``-Delta H``).
    d_t, d_p : float
        Tube and particle diameters (m).
    rho_b : float
        Catalyst bulk density (kg/m3).
    G : float
        Superficial mass flow velocity (kg/(m2 hr)).
    c_p : float
        Gas specific heat, mass basis (kcal/(kg degC)).
    M_m : float
        Mean molecular weight of the feed (kg/kmol).
    lam_R, alpha_w : float
        Effective radial conductivity (kcal/(m hr degC)) and wall coefficient
        (kcal/(m2 hr degC)).
    Pe_mR : float
        Radial mass Peclet number.
    L_bed : float
        Bed length (m).
    E, lnA : tuple of float
        Activation energies (cal/mol) and natural-log frequency factors of
        k1, k2, k3, as printed: ``ln k_i = -E_i / (R_gas T) + lnA_i``.
    R_gas : float
        Gas constant as printed in the Arrhenius lines (1.98 cal/(mol K)).
    kelvin : float, optional
        Offset from degC to K. Not printed; 273.15 by default.
    """

    NA0: float
    N0: float
    dH1: float
    dH3: float
    d_t: float
    d_p: float
    rho_b: float
    G: float
    c_p: float
    M_m: float
    lam_R: float
    alpha_w: float
    Pe_mR: float
    L_bed: float
    E: tuple
    lnA: tuple
    R_gas: float
    kelvin: float = 273.15

    @classmethod
    def from_table(cls, table, **changes):
        """Build the case from ``{symbol: value}`` as in the page's parameter CSV."""
        P = dict(table)
        kw = dict(NA0=P["N_A0"], N0=P["N_0"], dH1=P["minus_dH1"], dH3=P["minus_dH3"],
                  d_t=P["d_t"], d_p=P["d_p"], rho_b=P["rho_b"], G=P["G"], c_p=P["c_p"],
                  M_m=P["M_m"], lam_R=P["lambda_R"], alpha_w=P["alpha_w"],
                  Pe_mR=P["Pe_mR"], L_bed=P["L_bed"],
                  E=(P["E1"], P["E2"], P["E3"]), lnA=(P["lnA1"], P["lnA2"], P["lnA3"]),
                  R_gas=P["R_gas"])
        kw.update(changes)
        return cls(**kw)

    # --- the dimensionless groups of Eq. 1 (p. 24); z = z'/d_p, r = r'/d_p ---
    @property
    def R_tube(self):
        """Tube radius R' (m)."""
        return self.d_t / 2

    @property
    def R(self):
        """Dimensionless tube radius R'/d_p."""
        return self.R_tube / self.d_p

    @property
    def z_end(self):
        """Dimensionless bed length L/d_p."""
        return self.L_bed / self.d_p

    @property
    def a1(self):
        """Radial mass dispersion, 1/Pe_mR."""
        return 1.0 / self.Pe_mR

    @property
    def a2(self):
        """Radial heat dispersion, lam_R/(G c_p d_p) = 1/Pe_hR."""
        return self.lam_R / (self.G * self.c_p * self.d_p)

    @property
    def alpha(self):
        """Dimensionless wall film coefficient alpha_w d_p / lam_R."""
        return self.alpha_w * self.d_p / self.lam_R

    @property
    def b1(self):
        """Conversion scale rho_b d_p M_m / (G N_A0)."""
        return self.rho_b * self.d_p * self.M_m / (self.G * self.NA0)

    @property
    def b2(self):
        """Temperature scale of A -> B, rho_b d_p (-dH1) / (G c_p) (degC)."""
        return self.rho_b * self.d_p * self.dH1 / (self.G * self.c_p)

    @property
    def b3(self):
        """Temperature scale of the combustion, rho_b d_p (-dH3) / (G c_p) (degC)."""
        return self.rho_b * self.d_p * self.dH3 / (self.G * self.c_p)

    def wall_sink(self, U):
        """Coefficient of the 1-D wall term, 4 U d_p / (G c_p d_t), per unit z."""
        return 4.0 * U * self.d_p / (self.G * self.c_p * self.d_t)


def froment_source(case, T0C, *, b=None, heat_scale=1.0, k_order=(0, 1, 2)):
    """Pointwise source of the triangular scheme, for state ``(..., 3)``.

    Fields: ``x`` conversion to phthalic anhydride, ``w`` conversion to
    CO + CO2, ``t`` temperature rise above the inlet (= coolant) temperature.
    With ``y = x + w``,
    ``r_B = N_A0 N_0 [k1 (1 - y) - k2 x]`` and ``r_C = N_A0 N_0 [k2 x + k3 (1 - y)]``;
    the source is ``(b1 r_B, b1 r_C, b2 r_B + b3 r_C)``.

    Parameters
    ----------
    case : FromentCase
        Constants and groups.
    T0C : float
        Inlet and coolant temperature (degC).
    b : tuple of float, optional
        ``(b1, b2, b3)`` instead of the case's groups.
    heat_scale : float, optional
        Factor on ``b2`` and ``b3``. ``N_A0`` cancels from the conversion
        equations, so an inlet concentration ``s`` times the base one is exactly
        ``heat_scale = s``.
    k_order : tuple of int, optional
        Which Arrhenius line serves as k1, k2, k3 (``(0, 2, 1)`` swaps k2 and k3;
        used for a break row).

    Returns
    -------
    callable
        ``source(u) -> ndarray`` of the shape of ``u``.
    """
    b1, b2, b3 = (case.b1, case.b2, case.b3) if b is None else b
    b2, b3 = b2 * heat_scale, b3 * heat_scale
    lnA = np.asarray(case.lnA, dtype=float)[list(k_order)]
    E = np.asarray(case.E, dtype=float)[list(k_order)]
    Rg, kelvin = case.R_gas, case.kelvin
    scale = case.NA0 * case.N0

    def source(u):
        x, w, t = u[..., 0], u[..., 1], u[..., 2]
        T = t + T0C + kelvin
        k1 = np.exp(lnA[0] - E[0] / (Rg * T))
        k2 = np.exp(lnA[1] - E[1] / (Rg * T))
        k3 = np.exp(lnA[2] - E[2] / (Rg * T))
        y = x + w
        rB = scale * (k1 * (1 - y) - k2 * x)
        rC = scale * (k2 * x + k3 * (1 - y))
        return np.stack([b1 * rB, b1 * rC, b2 * rB + b3 * rC], axis=-1)

    return source


def u_overall(lam_R, alpha_w, R_tube):
    """Overall coefficient of the 1-D model, ``1/U = 1/alpha_w + R'/(4 lam_R)``.

    The paper's zero-order-moment match to the 2-D exchanger (pp. 25-26).
    """
    return 1.0 / (1.0 / alpha_w + R_tube / (4.0 * lam_R))


class Tube2D:
    """Two-dimensional pseudo-homogeneous fixed-bed tube, marched in ``z``.

    Every field obeys ``du_k/dz = a_k (1/r) d/dr (r du_k/dr) + s_k(u)`` on
    ``0 < r < R``, with symmetry at the axis and
    ``du_k/dr = -h_k (u_k - u_wall_k)`` at the wall (``h_k = 0``: no flux).
    State ``(n_r, n_f)``, spatial axis first; flattened only at the integrator
    interface.

    Assumptions, all built in: plug flow with the same unit convective
    coefficient for every field (the source and the coefficients are scaled
    to it, as in Froment's dimensionless form); no axial dispersion, so ``z``
    is a marching variable and the discretisation is radial only; ``a_k`` and
    ``h_k`` constant in ``r`` and ``z``; a wall value constant in ``z`` (it is
    folded into the constant ``g_disp``); a uniform radial grid. A model with
    a velocity profile, property variation or a coolant that changes along the
    tube needs more than a new source function.

    Parameters
    ----------
    source : callable
        ``source(u)`` for ``u`` of shape ``(n_r, n_f)``, returning the same
        shape. Pointwise: its Jacobian is taken by ``NumJac`` per cell.
    dispersion : array_like
        ``a_k`` per field, shape ``(n_f,)``.
    wall : array_like
        ``h_k`` per field, shape ``(n_f,)``.
    R : float
        Dimensionless tube radius.
    u_wall : array_like, optional
        Wall (coolant) values per field; only fields with ``h_k != 0`` use it.
    n_r : int, optional
        Number of radial cells (uniform).
    nu : int, optional
        Geometry index of ``construct_div``: 1, cylindrical, for a tube. Exists
        for the page's break table (the radial mean keeps the tube's area
        weights whatever ``nu`` is).

    Attributes
    ----------
    jac_disp : scipy.sparse.csc_array
        Constant dispersion operator ``div @ (a grad)``.
    g_disp : numpy.ndarray
        Its constant boundary contribution (zero for a zero wall value).
    weights : numpy.ndarray
        Annular area fractions for the radial mean, shape ``(n_r,)``.
    """

    def __init__(self, source, dispersion, wall, R, u_wall=0.0, n_r=80, nu=1):
        self.source, self.R, self.n_r, self.nu = source, float(R), int(n_r), nu
        self.dispersion = np.asarray(dispersion, dtype=float)
        self.wall = np.asarray(wall, dtype=float)
        self.n_f = self.dispersion.size
        self.u_wall = np.broadcast_to(np.asarray(u_wall, dtype=float), (self.n_f,))
        self.shape = (self.n_r, self.n_f)
        self._build_grid()
        self._build_operators()

    def _build_grid(self):
        self.r_f = np.linspace(0.0, self.R, self.n_r + 1)
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        # the radial mean zeta_m = 2 int_0^1 zeta (r/R) d(r/R): annular areas
        self.weights = (self.r_f[1:] ** 2 - self.r_f[:-1] ** 2) / self.R ** 2

    def _build_operators(self):
        # r = 0, outward normal -r: symmetry, du/dr = 0       -> a = 1, b = 0, d = 0
        # r = R, outward normal +r: du/dr = -h (u - u_wall)
        #                           -> du/dn + h u = h u_wall  -> a = 1, b = h, d = h u_wall
        bc = ({"a": 1.0, "b": 0.0, "d": 0.0},
              {"a": 1.0, "b": self.wall[None, :], "d": (self.wall * self.u_wall)[None, :]})
        grad_mat, grad_bc = construct_grad(self.shape, self.r_f, self.r_c, bc, axis=0)
        div_mat = construct_div(self.shape, self.r_f, nu=self.nu, axis=0)  # nu=1: cylinder
        disp_mat = construct_coefficient_matrix(self.dispersion[None, :], self.shape, axis=0)
        self.jac_disp = (div_mat @ disp_mat @ grad_mat).tocsc()
        self.g_disp = np.asarray((div_mat @ disp_mat @ grad_bc).todense()).ravel()
        self.numjac = NumJac(self.shape)          # pointwise source: last axis in full

    def rhs(self, z, y):
        """``dy/dz`` for the flat state ``y``."""
        return self.jac_disp @ y + self.g_disp + self.source(y.reshape(self.shape)).ravel()

    def jac(self, z, y):
        """Jacobian of :meth:`rhs`: the constant operator plus ``NumJac`` of the source."""
        _, jac_src = self.numjac(self.source, y.reshape(self.shape))
        return self.jac_disp + jac_src

    def march(self, z_end, u0=None, cap=None, field=-1, rtol=1e-8, atol=None):
        """Integrate from ``z = 0`` to ``z_end`` (BDF with the exact sparse Jacobian).

        Parameters
        ----------
        z_end : float
            Dimensionless length.
        u0 : array_like, optional
            Inlet state, shape ``(n_r, n_f)`` or broadcastable; zero by default.
        cap : float, optional
            Stop when the radial mean of ``field`` rises through ``cap``
            (``self.ran_away`` is then True).
        rtol, atol : float, optional
            Integrator tolerances; ``atol`` defaults to ``rtol / 100``.

        Returns
        -------
        Tube2D
            ``self``, with ``sol`` (the ``solve_ivp`` result) and ``ran_away``.
        """
        y0 = np.broadcast_to(0.0 if u0 is None else np.asarray(u0, float),
                             self.shape).ravel().copy()
        events = None
        if cap is not None:
            n_f, f = self.n_f, field % self.n_f

            def hit_cap(z, y):
                return self.weights @ y[f::n_f] - cap
            hit_cap.terminal, hit_cap.direction = True, 1
            events = hit_cap
        self.sol = solve_ivp(self.rhs, (0.0, z_end), y0, method="BDF", jac=self.jac,
                             rtol=rtol, atol=rtol * 1e-2 if atol is None else atol,
                             dense_output=True, events=events)
        if not self.sol.success:
            raise RuntimeError(f"Tube2D.march failed: {self.sol.message}")
        self.ran_away = events is not None and bool(self.sol.t_events[0].size)
        return self

    def state(self, z):
        """State at ``z`` (scalar or 1-D array): shape ``(n_r, n_f)`` or ``(n_z, n_r, n_f)``."""
        y = self.sol.sol(np.asarray(z, dtype=float))
        if y.ndim == 1:
            return y.reshape(self.shape)
        return y.T.reshape((-1,) + self.shape)

    def mean(self, z):
        """Radial (area) means of every field at ``z``, shape ``(..., n_f)``."""
        return np.einsum("...rf,r->...f", self.state(z), self.weights)

    def profiles(self, n_z=1500):
        """``(z, means, axis values)`` on ``n_z`` points up to where the march ended.

        The axis value is the first cell centre, ``r = h/2``.
        """
        z = np.linspace(0.0, self.sol.t[-1], n_z)
        u = self.state(z)
        return z, np.einsum("zrf,r->zf", u, self.weights), u[:, 0, :]

    def hot_spot(self, field=-1, n_z=4000):
        """``(z, value)`` at the maximum of the radial mean of ``field``.

        Returns ``(nan, inf)`` when the march stopped at its cap.
        """
        if self.ran_away:
            return np.nan, np.inf
        z, means, _ = self.profiles(n_z)
        i = int(np.argmax(means[:, field]))
        return z[i], means[i, field]


class Tube1D:
    """Radially lumped tube, ``du/dz = s(u) - k_w (u - u_wall)``, marched in ``z``.

    Parameters
    ----------
    source : callable
        ``source(u)`` for ``u`` of shape ``(n_f,)`` (or ``(..., n_f)``).
    wall : array_like
        Wall exchange coefficient ``k_w`` per field (for Froment's model
        ``4 U d_p / (G c_p d_t)`` on the temperature, 0 elsewhere).
    u_wall : array_like, optional
        Wall values per field.
    """

    def __init__(self, source, wall, u_wall=0.0):
        self.source = source
        self.wall = np.asarray(wall, dtype=float)
        self.u_wall = np.broadcast_to(np.asarray(u_wall, dtype=float), self.wall.shape)

    def rhs(self, z, u):
        """``du/dz``."""
        return self.source(u) - self.wall * (u - self.u_wall)

    def march(self, z_end, u0=None, cap=None, field=-1, rtol=1e-10, max_step=None):
        """Integrate with LSODA; stop when ``field`` rises through ``cap``."""
        y0 = np.zeros(self.wall.shape) if u0 is None else np.asarray(u0, float)
        events = None
        if cap is not None:
            def hit_cap(z, u):
                return u[field] - cap
            hit_cap.terminal, hit_cap.direction = True, 1
            events = hit_cap
        self.sol = solve_ivp(self.rhs, (0.0, z_end), y0, method="LSODA", rtol=rtol,
                             atol=rtol * 1e-2, dense_output=True, events=events,
                             max_step=np.inf if max_step is None else max_step)
        if not self.sol.success:
            raise RuntimeError(f"Tube1D.march failed: {self.sol.message}")
        self.ran_away = events is not None and bool(self.sol.t_events[0].size)
        return self

    def hot_spot(self, field=-1, n_z=4000):
        """``(z, value)`` at the maximum of ``field``; ``(nan, inf)`` after a runaway."""
        if self.ran_away:
            return np.nan, np.inf
        z = np.linspace(0.0, self.sol.t[-1], n_z)
        u = self.sol.sol(z)[field]
        i = int(np.argmax(u))
        return z[i], u[i]


class Tube1DSteady:
    """The lumped tube as a steady boundary-value problem on ``n_z`` cells.

    ``d(u)/dz = s(u) - k_w (u - u_wall)``, first-order upwind convection with
    unit velocity, solved by ``pymrm.newton``. A second discretisation of the
    same equations as :class:`Tube1D` (steady finite volumes against an
    adaptive march).

    Parameters
    ----------
    source, wall, u_wall
        As for :class:`Tube1D`.
    z_end : float
        Dimensionless length.
    u_in : array_like, optional
        Inlet values per field (zero by default).
    n_z : int, optional
        Number of cells.
    """

    def __init__(self, source, wall, z_end, u_wall=0.0, u_in=0.0, n_z=4000):
        self.source = source
        self.wall = np.asarray(wall, dtype=float)
        self.n_f = self.wall.size
        self.u_wall = np.broadcast_to(np.asarray(u_wall, dtype=float), (self.n_f,))
        self.u_in = np.broadcast_to(np.asarray(u_in, dtype=float), (self.n_f,))
        self.shape = (int(n_z), self.n_f)
        self.z_f = np.linspace(0.0, z_end, self.shape[0] + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self._build_operators()

    def _build_operators(self):
        # z = 0, outward normal -z: u = u_in (all fields)   -> a = 0, b = 1, d = u_in
        # z = L, outward normal +z: pure outflow, du/dz = 0 -> a = 1, b = 0, d = 0
        bc = ({"a": 0.0, "b": 1.0, "d": self.u_in[None, :]}, {"a": 1.0, "b": 0.0, "d": 0.0})
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c, bc,
                                                      v=1.0, axis=0)
        div_mat = construct_div(self.shape, self.z_f, nu=0, axis=0)   # nu=0: axial
        sink = construct_coefficient_matrix(np.broadcast_to(self.wall, self.shape))
        self.jac_const = (div_mat @ conv_mat + sink).tocsc()
        self.g_const = (div_mat @ conv_bc).toarray().reshape(-1, 1) \
            - np.broadcast_to(self.wall * self.u_wall, self.shape).reshape(-1, 1)
        self.numjac = NumJac(self.shape)           # pointwise source: last axis in full

    def residual(self, u):
        """Steady residual and Jacobian: convection + wall exchange - source."""
        g_src, jac_src = self.numjac(self.source, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_src.reshape((-1, 1))
        return g, self.jac_const - jac_src

    def solve(self, u0, scale=None, tol=1e-10, maxfev=50):
        """Newton from ``u0`` (shape ``(n_z, n_f)``); raises if it does not converge.

        ``scale`` (per field) makes the step test relative to each field's size;
        ``tol`` applies to the scaled step.
        """
        result = newton(self.residual, np.asarray(u0, dtype=float).reshape(self.shape),
                        tol=tol, maxfev=maxfev, scale=scale)
        if not result.success:
            raise RuntimeError(f"Tube1DSteady.solve: {result.message}")
        self.u, self.result = result.x.reshape(self.shape), result
        return self

    def hot_spot(self, field=-1):
        """``(z, value)`` at the largest cell value of ``field``."""
        i = int(np.argmax(self.u[:, field]))
        return self.z_c[i], float(self.u[i, field])


def critical_inlet(runs_away, lo, hi, tol):
    """Bisect ``runs_away(T0) -> bool`` on ``[lo, hi]`` to a bracket below ``tol``.

    Both ends are evaluated first: ``runs_away(lo)`` must be False and
    ``runs_away(hi)`` True, else ``ValueError``. Every probe is an independent
    solve from the inlet, so the result is a deterministic function of the
    bracket and the tolerance.

    Returns
    -------
    tuple
        ``(midpoint, lo, hi)`` of the final bracket.
    """
    if runs_away(lo) or not runs_away(hi):
        raise ValueError(f"critical_inlet: [{lo}, {hi}] does not bracket the threshold")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if runs_away(mid):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi), lo, hi


def exchanger_mean(z, a2, alpha, R, n_max=150.0):
    """Radial-mean temperature of the reactionless exchanger, the paper's series (p. 25).

    ``t_m/t_0 = 4 (alpha R)^2 sum exp(-(a2/R^2) l_n^2 z) / (l_n^2 (l_n^2 + (alpha R)^2))``
    over the roots ``l_n`` of ``l J1(l) = alpha R J0(l)`` below ``n_max``, for a
    uniform inlet ``t_0`` and wall (coolant) value 0.

    Returns
    -------
    tuple
        ``(t_m / t_0, roots)``.
    """
    aR = alpha * R
    f = lambda s: s * j1(s) - aR * j0(s)       # noqa: E731
    xs = np.linspace(1e-6, n_max, 30001)
    fs = f(xs)
    roots = np.array([brentq(f, xs[i], xs[i + 1])
                      for i in range(len(xs) - 1) if fs[i] * fs[i + 1] < 0])
    tm = 4 * aR ** 2 * np.sum(np.exp(-(a2 / R ** 2) * roots ** 2 * z)
                              / (roots ** 2 * (roots ** 2 + aR ** 2)))
    return tm, roots
