"""Five fixed-bed models of Froment, De Wilde and Bischoff (2011, ch. 11) on one reactor.

Model code of gallery page D1.1. One assembly builds all five rungs of the
chapter's model classification on section 11.5.2's hydrocarbon-oxidation design
case, with the pseudo-first-order rate held fixed:

====  =======  =========  ==============================================
rung  section  structure  what it adds
====  =======  =========  ==============================================
1     11.5     S2         basic 1-D pseudohomogeneous plug flow
2     11.6     S4         + axial dispersion, Danckwerts inlet
3     11.8     S7         + a solid phase behind a film
4     11.9     S8         + intraparticle profiles, via eq. (11.9.1-11)
5     11.7     S6         2-D, radial dispersion and a wall coefficient
====  =======  =========  ==============================================

Units are Van Welsenaere and Froment's (1970): kmol, kg, m, h, kcal, K and
pressures in atm. Every axial coordinate is a length in m; the velocity enters
the source term, so each operator acts on ``(n_z, n_f)`` (rung 5:
``(n_z, n_r, 2)``) with the spatial axes first and the fields last.

Contents
--------
KCAL
    kJ per thermochemical kcal.
Reactor
    The constants of the design case and the groups derived from them.
film
    Wakao and Funazkri's Sherwood correlation and the analogous Nusselt number.
eta_sphere, phi_prime
    Eqs. (11.9.1-11) and (11.9.1-12).
Ladder
    The five rungs, one finite-volume assembly, solved with ``pymrm.newton``.
Pellet
    Eqs. (11.9.1-3) and (11.9.1-4) on a spherical grid, one node of the bed.
hotspot
    A hot spot root-found on a local cubic, never a grid maximum.
guess_from_rung1, solve_robust, T_max_safe
    Fixed-inlet solves from a deterministic guess.
NotACrossing, crossing_control, p0_critical
    The inlet partial pressure at which the hot spot reaches a ceiling, by
    bisection along the start-up branch, checked on every call.
cross_scan
    The parameter value at which a threshold shift crosses a target, between
    two points the crossing control validates.
ode_reference, p0_critical_ode
    Rung 1 by an adaptive Radau solve with no grid and no Newton.
richardson
    Observed order and Richardson limit from the last three of a chain.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from pymrm import (NumJac, construct_coefficient_matrix, construct_convflux_upwind,
                   construct_div, construct_grad, newton)

__all__ = ["KCAL", "Reactor", "film", "eta_sphere", "phi_prime", "Ladder", "Pellet",
           "hotspot", "guess_from_rung1", "solve_robust", "T_max_safe", "RAMP",
           "CTL_DELTA", "CTL_JUMP", "CTL_STATS", "NotACrossing", "crossing_control",
           "p0_critical", "cross_scan", "ode_reference", "p0_critical_ode",
           "richardson"]

KCAL = 4.1868                      # kJ per thermochemical kcal


@dataclass(frozen=True)
class Reactor:
    """The design case of section 11.5.2, in Van Welsenaere and Froment's units.

    Parameters
    ----------
    M : float
        Mean molar mass of the gas (kg/kmol).
    p_t : float
        Total pressure (atm).
    rho_b, rho_g : float
        Catalyst bulk density and gas density (kg/m3).
    minus_dH : float
        Heat of reaction, positive for exothermic (kcal/kmol).
    cp_vol : float
        Volumetric heat capacity of the gas (kcal/(m3 K)), as printed.
    u : float
        Superficial velocity (m/h).
    U : float
        Overall wall heat-transfer coefficient (kcal/(m2 h K)).
    R : float
        Tube radius (m).
    p_B0 : float
        Oxygen partial pressure, in excess (atm).
    a, b : float
        Arrhenius constants of ``k = exp(b - a/T)``, ``a = E/R_gas`` (K).
    length : float
        Tube length (m).
    d_p : float
        Particle diameter (m).
    T_r : float
        Coolant temperature, equal to the inlet temperature (K).
    T_perm : float
        Top of the stated operating range, the ceiling of every threshold (K).
    eps : float
        Bed voidage.
    lam_er, alpha_w : float
        Effective radial conductivity (kJ/(m s K)) and wall coefficient
        (kJ/(m2 s K)) of rung 5.
    pe_mr, pe_hr : float
        Radial Peclet numbers for mass and heat.
    re_film : float
        The printed particle Reynolds number that fixes the gas viscosity of
        the film correlation.
    """

    M: float
    p_t: float
    rho_b: float
    rho_g: float
    minus_dH: float
    cp_vol: float
    u: float
    U: float
    R: float
    p_B0: float
    a: float
    b: float
    length: float
    d_p: float
    T_r: float
    T_perm: float
    eps: float
    lam_er: float
    alpha_w: float
    pe_mr: float
    pe_hr: float
    re_film: float

    @property
    def A(self):
        """Group A of eq. (11.5.2-1), ``M p_t rho_b p_B0 / rho_g`` (1/h per atm)."""
        return self.M * self.p_t * self.rho_b / self.rho_g * self.p_B0

    @property
    def B(self):
        """Group B, ``(-dH) rho_b p_B0 / (rho_g c_p)`` (K/h per atm)."""
        return self.minus_dH * self.rho_b / self.cp_vol * self.p_B0

    @property
    def theta(self):
        """Concentration per partial pressure, ``rho_g / (M p_t)`` (kmol/(m3 atm))."""
        return self.rho_g / (self.M * self.p_t)

    @property
    def cp_mass(self):
        """Mass heat capacity (kJ/(kg K))."""
        return self.cp_vol * KCAL / self.rho_g

    @property
    def rho_s(self):
        """Particle density (kg/m3)."""
        return self.rho_b / (1.0 - self.eps)

    @property
    def mass_flux(self):
        """Superficial mass flux ``G`` (kg/(m2 s))."""
        return self.rho_g * self.u / 3600.0

    def C(self, R=None):
        """Wall group ``C = 2 U / (c_p R)`` (1/h) at tube radius ``R``."""
        return 2.0 * self.U / (self.cp_vol * (self.R if R is None else R))

    def kexp(self, T):
        """Rate constant ``exp(b - a/T)``, clipped against overflow of an iterate."""
        return np.exp(np.clip(self.b - self.a / np.maximum(T, 1.0), -700, 50))


def film(rx, scale=1.0, lam=4.75e-5, Dm=1.0e-5):
    """Film coefficients from Wakao and Funazkri's Sherwood correlation (page A3.4).

    The gas viscosity is fixed by the printed particle Reynolds number
    ``rx.re_film``; ``lam`` and ``Dm`` are ordinary air values near 360 C that
    the chapter does not print.

    Parameters
    ----------
    rx : Reactor
        The design case.
    scale : float, optional
        Factor on both film coefficients.
    lam : float, optional
        Gas thermal conductivity (kJ/(m s K)).
    Dm : float, optional
        Molecular diffusivity (m2/s).

    Returns
    -------
    dict
        ``Re, Sc, Pr, Sh, Nu``; ``k_g`` (m/s), ``h_f`` (kJ/(m2 s K)), ``a_v``
        (m2/m3); ``kga`` (1/h) and ``hfa`` (kcal/(m3 h K)).
    """
    G = rx.mass_flux
    mu = rx.d_p * G / rx.re_film
    Re = rx.d_p * G / mu
    Sc = mu / (rx.rho_g * Dm)
    Pr = rx.cp_mass * mu / lam
    Sh = 2.0 + 1.1 * Sc ** (1 / 3) * Re ** 0.6
    Nu = 2.0 + 1.1 * Pr ** (1 / 3) * Re ** 0.6
    a_v = 6.0 * (1.0 - rx.eps) / rx.d_p
    k_g, h_f = scale * Sh * Dm / rx.d_p, scale * Nu * lam / rx.d_p
    return dict(Re=Re, Sc=Sc, Pr=Pr, Sh=Sh, Nu=Nu, k_g=k_g, h_f=h_f, a_v=a_v,
                kga=k_g * a_v * 3600.0,                   # 1/h
                hfa=h_f * a_v * 3600.0 / KCAL)            # kcal/(m3 h K)


def eta_sphere(phi):
    """Eq. (11.9.1-11): the isothermal effectiveness factor of a sphere, first order.

    Parameters
    ----------
    phi : float or array_like
        Thiele modulus on the particle radius.

    Returns
    -------
    float or ndarray
        ``3 (phi coth phi - 1) / phi^2``, by its series below ``phi = 1e-6``.
    """
    phi = np.atleast_1d(np.asarray(phi, float))
    out = np.empty_like(phi)
    s = phi < 1e-6
    out[s] = 1.0 - phi[s] ** 2 / 15.0
    q = phi[~s]
    out[~s] = 3.0 / q ** 2 * (q / np.tanh(q) - 1.0)
    return out if out.size > 1 else float(out[0])


def phi_prime(rx, T_s, D_e):
    """Eq. (11.9.1-12): Thiele modulus at surface temperature ``T_s`` (K).

    The rate constant in atm-kmol-kg-h is converted to a first-order 1/s, and
    ``D_e`` is in m2/s.
    """
    k_v = rx.kexp(T_s) * rx.p_B0 * rx.rho_s / (rx.theta * 3600.0)
    return (rx.d_p / 2.0) * np.sqrt(k_v / D_e)


class Ladder:
    """Chapter 11's five model configurations, one assembly.

    Parameters
    ----------
    rx : Reactor
        The design case.
    rung : {1, 2, 3, 4, 5}
        The model (see the module docstring).
    n_z, n_r : int
        Axial cells, and radial cells for rung 5.
    R, T_r, length : float, optional
        Tube radius (m), coolant and inlet temperature (K) and tube length (m);
        default the design case's.
    pe_ma, pe_ha : float
        Axial Peclet numbers for mass and heat on the particle diameter (rung 2).
    film_scale : float
        Factor on both film coefficients (rungs 3 and 4).
    D_e : float
        Effective diffusivity in the particle (m2/s, rung 4).
    pe_mr, pe_hr, alpha_w, lam_er : float, optional
        Radial Peclet numbers, wall coefficient and radial conductivity of
        rung 5; default the design case's.
    eta_on : bool
        Rung 4 with eq. (11.9.1-11); False forces ``eta = 1``.

    Notes
    -----
    Fields: ``(p, T)`` for rungs 1, 2 and 5, ``(p, T, p_s, T_s)`` for rungs 3
    and 4, where ``p_s, T_s`` are the algebraic surface values.
    """

    def __init__(self, rx, rung=1, n_z=1200, n_r=20, R=None, T_r=None, length=None,
                 pe_ma=1.5, pe_ha=1.5, film_scale=1.0, D_e=1.0e-6,
                 pe_mr=None, pe_hr=None, alpha_w=None, lam_er=None, eta_on=True):
        self.rx = rx
        R = rx.R if R is None else R
        T_r = rx.T_r if T_r is None else T_r
        length = rx.length if length is None else length
        pe_mr = rx.pe_mr if pe_mr is None else pe_mr
        pe_hr = rx.pe_hr if pe_hr is None else pe_hr
        alpha_w = rx.alpha_w if alpha_w is None else alpha_w
        lam_er = rx.lam_er if lam_er is None else lam_er
        self.rung, self.n_z, self.n_r = rung, n_z, n_r
        self.R, self.T_r, self.L = R, T_r, length
        self.C = rx.C(R)
        self.pe_ma, self.pe_ha, self.D_e, self.eta_on = pe_ma, pe_ha, D_e, eta_on
        self.f = film(rx, film_scale)
        self.pe_mr, self.pe_hr = pe_mr, pe_hr
        self.wall_ratio = alpha_w / lam_er                        # 1/m
        self.z_f = np.linspace(0.0, self.L, n_z + 1)
        self.z_c = 0.5 * (self.z_f[1:] + self.z_f[:-1])
        self.nf = 4 if rung in (3, 4) else 2
        if rung == 5:
            self.r_f = np.linspace(0.0, self.R, n_r + 1)
            self.r_c = 0.5 * (self.r_f[1:] + self.r_f[:-1])
            self.shape = (n_z, n_r, 2)
        else:
            self.shape = (n_z, self.nf)

    def _ops(self, p0, T0):
        """The constant linear operator and its boundary vector, for inlet ``(p0, T0)``.

        Assembled once per solve, never inside a Newton iteration.
        """
        d_p = self.rx.d_p
        nf = self.shape[-1]
        din = np.zeros(nf); din[0], din[1] = p0, T0
        bin_ = np.zeros(nf); bin_[0] = bin_[1] = 1.0
        if self.rung == 2:
            # eps*D_ea/u_s = d_p/Pe_ma and lambda_ea/(rho_g c_p u_s) = d_p/Pe_ha
            # are LENGTHS.  Danckwerts inlet, A2.1's published dict:
            #   u_s(c_in - c) = -D dc/dz at z=0; n = -z so dc/dn = -dc/dz,
            #   giving (D/u_s) dc/dn + c = c_in  ->  a = d_p/Pe, b = 1, d = c_in.
            # Outlet dc/dz = 0 with n = +z  ->  a = 1, b = 0, d = 0.
            Dv = np.array([d_p / self.pe_ma, d_p / self.pe_ha])
            bc = ({"a": Dv, "b": bin_, "d": din},
                  {"a": np.ones(nf), "b": np.zeros(nf), "d": np.zeros(nf)})
            grad, dgrad = construct_grad(self.shape, self.z_f, self.z_c, bc, axis=0)
            conv, dconv = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                    bc, v=1.0, axis=0)
            Dm = construct_coefficient_matrix(Dv, shape=(self.n_z + 1, 2))
            div = construct_div(self.shape, self.z_f, nu=0, axis=0)     # Cartesian
            return div @ (conv - Dm @ grad), div @ (dconv - Dm @ dgrad)
        if self.rung == 5:
            # z = 0: p = p0, T = T0 (a = 0, b = 1, d = inlet); z = L: dc/dz = 0, n = +z.
            bcz = ({"a": np.zeros(2), "b": np.ones(2), "d": np.array([p0, T0])},
                   {"a": np.ones(2), "b": np.zeros(2), "d": np.zeros(2)})
            conv, dconv = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                    bcz, v=1.0, axis=0)
            divz = construct_div(self.shape, self.z_f, nu=0, axis=0)     # Cartesian
            Dr = np.array([d_p / self.pe_mr, d_p / self.pe_hr])          # m
            # r = 0: symmetry, dc/dn = 0 for both fields (n = -r there).
            # r = R_t: p is impermeable; for T, eq. (11.7.1-1) with n = +r,
            #   lambda_er dT/dn + alpha_w T = alpha_w T_w, scaled by 1/(rho_g c_p u_s)
            #   so that a keeps the same units as the operator's diffusivity.
            bcr = ({"a": np.ones(2), "b": np.zeros(2), "d": np.zeros(2)},
                   {"a": np.array([1.0, Dr[1]]),
                    "b": np.array([0.0, Dr[1] * self.wall_ratio]),
                    "d": np.array([0.0, Dr[1] * self.wall_ratio * self.T_r])})
            gradr, dgradr = construct_grad(self.shape, self.r_f, self.r_c, bcr, axis=1)
            divr = construct_div(self.shape, self.r_f, nu=1, axis=1)     # cylindrical
            Dm = construct_coefficient_matrix(Dr, shape=(self.n_z, self.n_r + 1, 2))
            return (divz @ conv - divr @ (Dm @ gradr),
                    divz @ dconv - divr @ (Dm @ dgradr))
        # z = 0: fixed inlet values (a = 0, b = 1); z = L: dc/dz = 0, n = +z.
        bc = ({"a": np.zeros(nf), "b": bin_, "d": din},
              {"a": np.ones(nf), "b": np.zeros(nf), "d": np.zeros(nf)})
        conv, dconv = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                bc, v=1.0, axis=0)
        div = construct_div(self.shape, self.z_f, nu=0, axis=0)          # Cartesian
        if nf == 2:
            return div @ conv, div @ dconv
        # rungs 3 and 4: fields 2 and 3 are ALGEBRAIC, so the divergence must not
        # reach them.  A coefficient matrix is the mask.
        mask = construct_coefficient_matrix(np.array([1.0, 1.0, 0.0, 0.0]),
                                            shape=self.shape)
        return mask @ (div @ conv), mask @ (div @ dconv)

    def source(self, c):
        """The pointwise source term, per unit length, at state ``c``.

        A diverging Newton iterate far above a threshold drives ``p k(T)`` past
        the float range; that iterate is rejected by the residual check in
        ``solve``, so the floating-point warnings are suppressed here only.
        """
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            return self._source(c)

    def _source(self, c):
        rx = self.rx
        s = np.empty_like(c)
        if self.rung in (1, 2, 5):
            p, T = c[..., 0], c[..., 1]
            q = rx.A / rx.u * p * rx.kexp(T)                  # atm/m
            s[..., 0] = -q
            s[..., 1] = (rx.B / rx.A) * q
            if self.rung != 5:                 # rung 5 cools through the wall bc
                s[..., 1] -= self.C / rx.u * (T - self.T_r)
            return s
        p, T, ps, Ts = c[..., 0], c[..., 1], c[..., 2], c[..., 3]
        eta = (eta_sphere(phi_prime(rx, Ts, self.D_e))
               if (self.rung == 4 and self.eta_on) else 1.0)
        rate = eta * rx.A * ps * rx.kexp(Ts)                  # atm/h at SURFACE cond.
        s[..., 0] = -rate / rx.u
        s[..., 1] = (rx.B / rx.A) * rate / rx.u - self.C / rx.u * (T - self.T_r)
        # (11.8.1-3) and (11.8.1-4), each divided by its own transfer coefficient
        # so that all four residual components are the same size.
        s[..., 2] = rate / self.f["kga"] - (p - ps)
        s[..., 3] = rx.cp_vol * (rx.B / rx.A) * rate / self.f["hfa"] - (Ts - T)
        return s

    def solve(self, p0, T0=None, tol=1e-10, guess=None):
        """Steady state at inlet partial pressure ``p0`` (atm) and temperature ``T0``.

        Parameters
        ----------
        p0 : float
            Inlet partial pressure (atm).
        T0 : float, optional
            Inlet temperature (K); default the coolant temperature.
        tol : float
            Newton step tolerance.
        guess : ndarray, optional
            Initial state of shape ``self.shape``; default the flat inlet state.

        Returns
        -------
        ndarray
            The state, shape ``self.shape``.

        Raises
        ------
        RuntimeError
            When ``newton`` does not succeed, or the final residual is not
            finite or exceeds 1e-6.
        """
        T0 = self.T_r if T0 is None else T0
        Amat, g = self._ops(p0, T0)
        gv = (np.asarray(g.todense()).ravel() if hasattr(g, "todense")
              else np.asarray(g).ravel())
        jac = NumJac(self.shape)     # last axis in full: correct for a pointwise source

        def fun(x):
            cc = x.reshape(self.shape)
            return (Amat @ x + gv - self.source(cc).reshape(-1),
                    Amat - jac(self.source, cc)[1])
        if guess is None:
            x0 = np.empty(self.shape)
            x0[..., 0] = p0; x0[..., 1] = T0
            if self.nf == 4:
                x0[..., 2] = p0; x0[..., 3] = T0
        else:
            x0 = guess
        # stagnation=True: a state at the round-off floor, where the steps stop
        # shrinking above tol, ends the iteration with success (backward error
        # permitting) instead of running to maxfev; anything else that does not
        # succeed is raised, never used.
        sol = newton(fun, np.asarray(x0).reshape(-1), tol=tol, maxfev=100,
                     stagnation=True)
        if not sol.success:
            raise RuntimeError(f"rung {self.rung}: newton did not converge ({sol.message})")
        x = sol.x
        res = float(np.max(np.abs(fun(x)[0])))
        if not np.isfinite(res) or res > 1e-6:
            raise RuntimeError(f"rung {self.rung}: Newton residual {res:.3e}")
        self.last_residual = res
        return x.reshape(self.shape)

    def bulk_T(self, c):
        """Gas temperature along z; for rung 5 the radial mean of page 574."""
        if self.rung != 5:
            return c[..., 1]
        w = 2.0 * self.r_c * (self.r_f[1:] - self.r_f[:-1])
        return (c[..., 1] * w).sum(axis=1) / w.sum()

    def axis_T(self, c):
        """Temperature on the tube axis (rung 5) or the gas temperature."""
        return c[:, 0, 1] if self.rung == 5 else c[..., 1]


class Pellet:
    """Eqs. (11.9.1-3) and (11.9.1-4) on a spherical pymrm grid, one node of the bed.

    Parameters
    ----------
    rx : Reactor
        The design case.
    n_r : int
        Radial cells.
    D_e : float
        Effective diffusivity (m2/s).
    lam_e : float
        Effective conductivity (kJ/(m s K)).
    isothermal : bool
        Solve the mass balance only, at the surface temperature.
    nu : int
        Geometry of the divergence: 2 is the sphere; another value is a defect
        injected by the page's break table.
    """

    def __init__(self, rx, n_r=64, D_e=1.0e-6, lam_e=1.0e-3, isothermal=True, nu=2):
        self.rx = rx
        self.n_r, self.D_e, self.lam_e, self.iso = n_r, D_e, lam_e, isothermal
        self.nu = nu
        self.nf = 1 if isothermal else 2
        self.R_p = rx.d_p / 2.0
        self.xi_f = np.linspace(0.0, self.R_p, n_r + 1)
        self.xi_c = 0.5 * (self.xi_f[1:] + self.xi_f[:-1])
        # (n_r, nf), never a bare (n_r,): NumJac couples the last axis in full.
        self.shape = (n_r, self.nf)

    def solve(self, C_surf, T_surf, tol=1e-10):
        """Profiles at surface concentration ``C_surf`` (kmol/m3) and ``T_surf`` (K).

        ``tol`` is ``newton``'s absolute step tolerance; with the temperature as an
        unknown near 680 K a step of 1e-12 is below its round-off and the solve
        would never stop.

        Raises
        ------
        RuntimeError
            When ``newton`` does not succeed.

        Returns
        -------
        dict
            ``eta`` from the surface flux of the gradient operator, the state
            ``c``, the centre temperature ``T_centre`` and ``C_ratio``, the
            centre over the surface concentration.
        """
        rx, nf = self.rx, self.nf
        aD = np.array([self.D_e] + ([self.lam_e] if nf == 2 else []))
        # xi = 0 symmetry: dc/dn = 0.  xi = R_p: surface conditions imposed, so the
        # film sits outside this solve and eta is a pure intraparticle quantity.
        bc = ({"a": np.ones(nf), "b": np.zeros(nf), "d": np.zeros(nf)},
              {"a": np.zeros(nf), "b": np.ones(nf),
               "d": np.array([C_surf] + ([T_surf] if nf == 2 else []))})
        grad, dgrad = construct_grad(self.shape, self.xi_f, self.xi_c, bc, axis=0)
        div = construct_div(self.shape, self.xi_f, nu=self.nu, axis=0)    # 2: spherical
        Dm = construct_coefficient_matrix(aD, shape=(self.n_r + 1, nf))
        Lin = -(div @ (Dm @ grad))
        rhs = -(div @ (Dm @ dgrad))
        rhs = (np.asarray(rhs.todense()).ravel() if hasattr(rhs, "todense")
               else np.asarray(rhs).ravel())
        Ts0 = np.full(self.n_r, T_surf)

        def src(c):
            Cs = c[..., 0]
            Ts = c[..., 1] if nf == 2 else Ts0
            r = rx.kexp(Ts) * rx.p_B0 * (Cs / rx.theta) / 3600.0       # kmol/(kg s)
            s = np.empty_like(c)
            s[..., 0] = -rx.rho_s * r
            if nf == 2:
                s[..., 1] = rx.rho_s * r * rx.minus_dH * KCAL
            return s
        jac = NumJac(self.shape)

        def fun(x):
            cc = x.reshape(self.shape)
            return Lin @ x + rhs - src(cc).reshape(-1), Lin - jac(src, cc)[1]
        x0 = np.empty(self.shape); x0[..., 0] = C_surf
        if nf == 2:
            x0[..., 1] = T_surf
        sol = newton(fun, x0.reshape(-1), tol=tol, maxfev=60)
        if not sol.success:
            raise RuntimeError(f"pellet: newton did not converge ({sol.message})")
        c = sol.x.reshape(self.shape)
        gf = (grad @ c.reshape(-1) + (np.asarray(dgrad.todense()).ravel()
              if hasattr(dgrad, "todense") else np.asarray(dgrad).ravel())
              ).reshape(self.n_r + 1, nf)
        # eta from the SURFACE FLUX row of the gradient operator, never a volume mean
        flux = -self.D_e * gf[-1, 0]
        r_obs = -flux * (3.0 / self.R_p) / rx.rho_s
        r_surf = rx.kexp(T_surf) * rx.p_B0 * (C_surf / rx.theta) / 3600.0
        return dict(eta=float(r_obs / r_surf), c=c,
                    T_centre=float(c[0, 1]) if nf == 2 else T_surf,
                    C_ratio=float(c[0, 0] / C_surf))


def hotspot(z_c, T):
    """The maximum of ``T(z)``, root-found as ``dT/dz = 0`` on a local cubic.

    Returns
    -------
    tuple
        ``(T_max, z_max, at_edge)``; at an edge, or with no interior root, the
        grid value is returned and ``at_edge`` says which happened.
    """
    i = int(np.argmax(T))
    if i in (0, len(T) - 1):
        return float(T[i]), float(z_c[i]), True
    sl = slice(max(i - 3, 0), min(i + 4, len(T)))
    co = np.polyfit(z_c[sl], T[sl], 3)
    rt = [r.real for r in np.roots(np.polyder(co))
          if abs(r.imag) < 1e-12 and z_c[i - 1] <= r.real <= z_c[i + 1]]
    if not rt:
        return float(T[i]), float(z_c[i]), False
    zm = max(rt, key=lambda r: np.polyval(co, r))
    return float(np.polyval(co, zm)), float(zm), False


def guess_from_rung1(lad, p0):
    """A deterministic initial guess for ``lad`` at ``p0``: the rung-1 profile.

    A fixed function of ``(p0, grid)``, so no history of a sweep enters it.
    """
    r1 = Ladder(lad.rx, rung=1, n_z=lad.n_z, R=lad.R, T_r=lad.T_r, length=lad.L)
    try:
        c1 = r1.solve(p0)
    except RuntimeError:                       # rung 1 has itself run away
        c1 = np.column_stack([np.full(lad.n_z, p0), np.full(lad.n_z, lad.T_r)])
    g = np.empty(lad.shape)
    if lad.rung == 5:
        g[:, :, 0] = c1[:, [0]]; g[:, :, 1] = c1[:, [1]]
        return g
    g[..., 0] = c1[:, 0]; g[..., 1] = c1[:, 1]
    if lad.nf == 4:
        g[..., 2] = c1[:, 0]; g[..., 3] = c1[:, 1]
    return g


RAMP = (0.4, 0.6, 0.8, 1.0)      # a FIXED ladder of fractions, not a sweep


def solve_robust(lad, p0):
    """A fixed-inlet solve: from the rung-1 guess, else a FIXED four-step ramp in p0.

    The ramp is the same four fractions every time, so the result is a function
    of ``(lad, p0)`` alone.
    """
    try:
        return lad.solve(p0, guess=guess_from_rung1(lad, p0))
    except RuntimeError:
        pass
    g = guess_from_rung1(lad, RAMP[0] * p0)
    for frac in RAMP:
        g = lad.solve(frac * p0, guess=g)     # a failure here propagates, by design
    return g


def T_max_safe(lad, p0):
    """Hot spot of ``solve_robust(lad, p0)``, or NaN when no solve converges."""
    try:
        return hotspot(lad.z_c, lad.bulk_T(solve_robust(lad, p0)))[0]
    except Exception:
        return np.nan


CTL_DELTA = 1.0e-5     # the fractional step either side of a returned threshold
CTL_JUMP = 1.0         # K: the widest T(p0*(1+d)) - T(p0*(1-d)) a crossing may show
CTL_STATS = {"checked": 0, "rejected": 0, "margin_K": np.inf}
#   margin_K: the smallest distance, over every control run, of either control
#   solve from the ceiling or of the jump from its allowance.


class NotACrossing(RuntimeError):
    """A bisected threshold that is not a crossing of the ceiling.

    Raised rather than returned, so that no sweep, bisection or table can
    report the number by accident.
    """


def _T_from(lad, q, guess):
    """Hot spot of the steady state Newton reaches from ``guess``, and that state.

    ``(nan, None)`` if the solve does not converge.
    """
    try:
        c = lad.solve(q, guess=guess)
        return hotspot(lad.z_c, lad.bulk_T(c))[0], c
    except Exception:
        return np.nan, None


def crossing_control(lad, p0, ceiling=None, delta=CTL_DELTA, guess=None, jump=None):
    """Whether the hot spot crosses ``ceiling`` at ``p0``, and crosses it continuously.

    Solves at ``p0 (1 - delta)`` and ``p0 (1 + delta)`` from ``guess``, the
    converged state at the lower end of the bisection's final bracket (by
    default ``lad.branch``, set by ``p0_critical``).

    Returns
    -------
    tuple
        ``(T_below, T_above, ok)``: ``ok`` is True when both solves converge,
        lie on either side of the ceiling and differ by less than ``jump``
        (default ``CTL_JUMP``).
    """
    ceiling = lad.rx.T_perm if ceiling is None else ceiling
    if guess is None:
        pb, guess = lad.branch                      # set by p0_critical on this lad
        assert abs(pb / p0 - 1.0) < 1e-8, "crossing_control: no branch state at p0"
    lo, _ = _T_from(lad, p0 * (1.0 - delta), guess)
    hi, _ = _T_from(lad, p0 * (1.0 + delta), guess)
    jump = CTL_JUMP if jump is None else jump
    ok = bool(np.isfinite(lo) and np.isfinite(hi) and lo < ceiling < hi
              and abs(hi - lo) < jump)
    return lo, hi, ok


def p0_critical(lad, ceiling=None, lo=0.004, hi=0.0200, n_scan=81, iters=90,
                jump=None):
    """The inlet partial pressure at which the hot spot reaches ``ceiling``.

    Bisection along the start-up branch: one cold solve at ``lo``
    (``solve_robust``), then every scan point and every midpoint started from
    the converged state at the current lower end of the bracket, which moves up
    only on a sub-ceiling result. A non-converged solve counts as above the
    ceiling, so the result is checked by ``crossing_control`` on every call.

    Parameters
    ----------
    lad : Ladder
        The model.
    ceiling : float, optional
        Hot-spot ceiling (K); default ``lad.rx.T_perm``.
    lo, hi : float
        Bracket (atm).
    n_scan : int
        Points of the bracketing scan.
    iters : int
        Bisection halvings at most; the bracket stops at 1e-12.
    jump : float, optional
        The control's jump allowance (K); default ``CTL_JUMP``.

    Returns
    -------
    float
        The threshold (atm). The branch state is stored as ``lad.branch``.

    Raises
    ------
    NotACrossing
        When the control rejects the threshold.
    RuntimeError
        When the bracket holds no sub-ceiling point or no crossing.
    """
    ceiling = lad.rx.T_perm if ceiling is None else ceiling

    def side(q, guess):
        v, c = _T_from(lad, q, guess)
        return bool(np.isfinite(v) and v < ceiling), c
    qs = np.linspace(lo, hi, n_scan)
    try:
        c_a = solve_robust(lad, qs[0])              # the anchor: one cold state at lo
    except Exception:
        raise RuntimeError("no sub-ceiling point on this bracket")
    if not hotspot(lad.z_c, lad.bulk_T(c_a))[0] < ceiling:
        raise RuntimeError("no sub-ceiling point on this bracket")
    a = qs[0]
    for q in qs[1:]:
        ok, c = side(q, c_a)
        if ok:
            a, c_a = q, c
        else:
            b = q
            break
    else:
        raise RuntimeError("the hot spot never reaches the ceiling on this bracket")
    for _ in range(iters):
        m = 0.5 * (a + b)
        ok, c = side(m, c_a)
        if ok:
            a, c_a = m, c
        else:
            b = m
        if b - a < 1e-12 * max(1.0, b):
            break
    p = 0.5 * (a + b)
    lad.branch = (p, c_a)
    CTL_STATS["checked"] += 1
    lo_, hi_, ok = crossing_control(lad, p, ceiling, guess=c_a, jump=jump)
    if np.isfinite(lo_) and np.isfinite(hi_):
        CTL_STATS["margin_K"] = min(CTL_STATS["margin_K"], abs(lo_ - ceiling),
                                    abs(hi_ - ceiling),
                                    abs((CTL_JUMP if jump is None else jump)
                                        - abs(hi_ - lo_)))
    if not ok:
        CTL_STATS["rejected"] += 1
        raise NotACrossing(
            f"rung {lad.rung} at n_z = {lad.n_z}: p0* = {p:.8g} is not a "
            f"ceiling crossing - T(p0*(1-{CTL_DELTA:g})) = {lo_-273.15:.2f} C, "
            f"T(p0*(1+{CTL_DELTA:g})) = {hi_-273.15:.2f} C")
    return p


def _shift_or_none(quantity, x):
    """``(quantity(x), True)``, or ``(nan, False)`` when its threshold fails."""
    try:
        return float(quantity(x)), True
    except NotACrossing:
        return np.nan, False
    except RuntimeError:                       # no threshold inside the bracket
        return np.nan, False


def cross_scan(quantity, xs, target, levels=8, probes=3):
    """The parameter value at which ``|quantity|`` crosses ``target``.

    ``xs`` is swept outward from the chapter's own value until two consecutive
    points that pass the crossing control straddle the target; ``|quantity|``
    must be monotone over the validated points up to that pair. The bracket is
    then cut into ``probes + 1`` parts per level, keeping the interior points
    that validate (one finer retry per level), for at most ``levels`` levels.

    Returns
    -------
    dict
        ``value`` (the endpoint of the final bracket nearest the target, a
        validated point), ``shift_pct`` there, the sweep table and the counts
        of the search.

    Raises
    ------
    AssertionError
        When the sweep does not straddle the target between validated points.
    """
    rows, good, pair = [], [], None
    for x in xs:
        v, ok = _shift_or_none(quantity, float(x))
        rows.append(dict(x=float(x), shift_pct=v, controlled=ok))
        if not ok:
            continue
        good.append((float(x), float(v)))
        if len(good) >= 2 and (abs(good[-2][1]) - target) * (abs(v) - target) <= 0.0:
            mono = np.abs([g[1] for g in good])
            assert np.all(np.diff(mono) >= 0) or np.all(np.diff(mono) <= 0), (
                "cross_scan: |shift| is not monotone from the chapter's own value "
                f"to the straddling pair: {mono}")
            pair = (good[-2], good[-1])
            break
    tab = pd.DataFrame(rows)
    assert len(good) >= 2, "cross_scan: fewer than two validated sweep points"
    assert pair is not None, (
        "cross_scan: the sweep ends without straddling the target between two "
        f"VALIDATED points - it reaches |shift| = {max(abs(g[1]) for g in good):.4f}"
        f" against a target of {target}. Extend the sweep; a bound is not an answer.")

    (a, va), (b, vb) = pair
    bracket0, done, tried, refused, mid_fail = abs(b - a), 0, 0, 0, None
    for _ in range(levels):
        pts = [(a, va), (b, vb)]
        for k in (probes, 2 * probes + 1):     # ONE finer retry before giving up
            pts = [(a, va), (b, vb)]
            for j, xi in enumerate(a + (b - a) * np.arange(1, k + 1) / (k + 1)):
                vi, oki = _shift_or_none(quantity, float(xi))
                tried += 1
                refused += 0 if oki else 1
                if k % 2 and j == (k - 1) // 2 and not oki and mid_fail is None:
                    mid_fail = done + 1
                if oki:
                    pts.append((float(xi), float(vi)))
            if len(pts) > 2:
                break
        if len(pts) == 2:                      # nothing inside validates: stop here
            break
        pts.sort(key=lambda t: t[0])
        for (x0, v0), (x1, v1) in zip(pts, pts[1:]):
            if (abs(v0) - target) * (abs(v1) - target) <= 0.0:
                (a, va), (b, vb) = (x0, v0), (x1, v1)
                break
        done += 1
    x_star, v_star = min(((a, va), (b, vb)), key=lambda t: abs(abs(t[1]) - target))
    return dict(kind="root", value=x_star, shift_pct=v_star, sweep=tab,
                n_controlled=len(good), n_evaluated=len(rows), levels=done,
                probed=tried, refused=refused, mid_fail=mid_fail,
                bracket=abs(b - a), bracket0=bracket0,
                pair=(pair[0][0], pair[1][0]), pair_shift=(pair[0][1], pair[1][1]))


def ode_reference(rx, p0, T_r=None, R=None, length=None, rtol=1e-12, atol=1e-14):
    """Rung 1 by an adaptive Radau solve: no grid, no pymrm operator, no Newton.

    Returns
    -------
    tuple
        ``(T_max, z_max, solution)``, the maximum root-found as ``dT/dz = 0``.
    """
    T_r = rx.T_r if T_r is None else T_r
    R = rx.R if R is None else R
    length = rx.length if length is None else length
    C = rx.C(R)
    A, B = rx.A, rx.B

    def rhs(z, y):
        p, T = y
        q = A / rx.u * p * rx.kexp(T)
        return [-q, (B / A) * q - C / rx.u * (T - T_r)]
    s = solve_ivp(rhs, (0.0, length), [p0, T_r], method="Radau", rtol=rtol,
                  atol=atol, dense_output=True)
    zg = np.linspace(0.0, length, 8001)
    Tg = s.sol(zg)[1]
    i = int(np.argmax(Tg))
    if 0 < i < len(zg) - 1:
        zm = brentq(lambda z: rhs(z, s.sol(z))[1], zg[i - 1], zg[i + 1],
                    xtol=1e-14, rtol=8.9e-16)
        return float(s.sol(zm)[1]), float(zm), s
    return float(Tg[i]), float(zg[i]), s


def p0_critical_ode(rx, T_r=None, R=None, ceiling=None, lo=0.010, hi=0.0195):
    """Rung 1's threshold from ``ode_reference``, root-found with brentq (atm)."""
    ceiling = rx.T_perm if ceiling is None else ceiling
    return brentq(lambda q: ode_reference(rx, q, T_r=T_r, R=R)[0] - ceiling,
                  lo, hi, xtol=1e-12, rtol=8.9e-16)


def richardson(vals, ns):
    """Observed order and the ``h -> 0`` limit from the last three of a chain.

    Parameters
    ----------
    vals : array_like
        The quantity on each grid.
    ns : array_like
        Cells of each grid.

    Returns
    -------
    tuple
        ``(order, limit)``.
    """
    v, n = np.asarray(vals, float), np.asarray(ns, float)
    order = float(np.log(abs((v[-3] - v[-2]) / (v[-2] - v[-1]))) / np.log(n[-1] / n[-2]))
    lim = float(v[-1] + (v[-1] - v[-2]) / ((n[-1] / n[-2]) ** order - 1.0))
    return order, lim
