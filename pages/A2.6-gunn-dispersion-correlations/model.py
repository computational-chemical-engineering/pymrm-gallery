"""Gunn's (1987) dispersion correlations and a 2-D axisymmetric fixed bed.

Model code of gallery page A2.6 (Gunn, Chem. Eng. Sci. 42 (1987) 363-373). The
page imports these; the closed forms the bed is checked against (Gunn's eq. 17,
its Bessel-mode generalisation) and the collocation route stay in the notebook,
so that the pymrm solve and the reference never share a line of code.

Every printed constant enters through one dictionary ``G``, the rows of the
page's ``gunn-1987-stated.csv`` as ``{quantity: value}``, which every
correlation function takes. A break-table row that mis-reads a constant is then
one call with ``{**G, key: wrong_value}``.

Contents
--------
p_of_Re
    Eqs. (3) to (5), the probability ``p(Re)`` for the three particle shapes.
D_over_Ud, dD_over_Ud_dRe, Pe_axial
    Eq. (1), its analytic derivative in Re, and the axial Peclet group.
D_over_Ud_1993
    The same model as Gunn (1993) prints it, eq. (42) with
    ``Gamma = 4 (1 - eps) alpha_1**2 / eps``.
Pe_f, Pe_radial
    Eqs. (10) to (12) and eq. (9), the radial Peclet group.
Bed
    Gunn's eq. (16) in two dimensions, axial and radial dispersion with a
    first-order reaction, on the pymrm operators.

Notes
-----
The bed equations are linear: a first-order reaction, constant coefficients
and a fixed velocity profile. With van Leer convection the discrete equations
are nonlinear only through the limiter. ``Bed.solve`` starts from the
first-order upwind solution and solves the limited equations by Newton with
the exact ``ConvFluxTVD`` Jacobian (pitfalls P10); the residual is written as
the correction to the upwind solution, as on page A2.1. The same equations
solved by van Leer deferred correction (the limiter correction lagged, one LU
of the upwind operator reused) are the second solver route.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu
from scipy.special import j0, jn_zeros

from pymrm import (ConvFluxTVD, compute_boundary_values, construct_convflux_upwind,
                   construct_div, construct_grad, interp_cntr_to_stagg_tvd, newton,
                   vanleer)

__all__ = ["ALPHA1", "BETA1", "SHAPES", "p_of_Re", "D_over_Ud", "dD_over_Ud_dRe",
           "D_over_Ud_1993", "Pe_axial", "Pe_f", "Pe_radial", "Bed"]

#: Eq. (1): "alpha_1 is the first root of J_0(u) = 0".
ALPHA1 = float(jn_zeros(0, 1)[0])
#: First root of J_1: dJ_0/dr = 0 at the wall, the first zero-flux radial mode.
BETA1 = float(jn_zeros(1, 1)[0])

#: Per shape, the keys of ``G``: amplitude of p, axial tortuosity, Pe_f
#: asymptote, Pe_f amplitude, radial tortuosity.
SHAPES = {
    "spheres": ("p_amplitude_spheres", "tau_axial_spheres", "pef_asymptote_spheres",
                "pef_amplitude_spheres", "tau_radial_spheres"),
    "solid cylinders": ("p_amplitude_solid_cylinders", "tau_axial_solid_cylinders",
                        "pef_asymptote_solid_cylinders", "pef_amplitude_solid_cylinders",
                        "tau_radial_solid_cylinders"),
    "hollow cylinders": ("p_amplitude_hollow_cylinders", "tau_axial_hollow_cylinders",
                         "pef_asymptote_hollow_cylinders", "pef_amplitude_hollow_cylinders",
                         "tau_radial_hollow_cylinders"),
}


def p_of_Re(Re, G, shape="spheres"):
    """Eqs. (3) to (5): ``p = p_offset + amplitude exp(-scale / Re)``.

    Parameters
    ----------
    Re : float or numpy.ndarray
        Reynolds number on the superficial velocity, ``d U_0 rho / mu``.
    G : dict
        Printed constants, ``{quantity: value}`` of ``gunn-1987-stated.csv``.
    shape : {"spheres", "solid cylinders", "hollow cylinders"}, optional
        Particle shape, selecting the amplitude.

    Returns
    -------
    numpy.ndarray
        The probability ``p``.
    """
    Re = np.asarray(Re, float)
    return G["p_offset"] + G[SHAPES[shape][0]] * np.exp(-G["p_reynolds_scale"] / Re)


def _dp_dRe(Re, G, shape):
    """Derivative of eqs. (3) to (5) in Re."""
    Re = np.asarray(Re, float)
    s = G["p_reynolds_scale"]
    return G[SHAPES[shape][0]] * s / Re ** 2 * np.exp(-s / Re)


def D_over_Ud(Re, Sc, G, eps, shape="spheres", p=None, alpha1=ALPHA1,
              keep_exp_term=True, stable=True):
    """Eq. (1), the axial dispersion group ``D / (U d)``.

    Parameters
    ----------
    Re, Sc : float or numpy.ndarray
        Reynolds number (superficial) and Schmidt number.
    G : dict
        Printed constants (see :func:`p_of_Re`); the axial tortuosity of the
        shape is read from it.
    eps : float
        Bed porosity (not printed in the paper; the page states its choice).
    shape : str, optional
        Particle shape.
    p : float or numpy.ndarray, optional
        Overrides ``p(Re)``: a frozen or free constant ``p`` for the nulls.
    alpha1 : float, optional
        First root of ``J_0``; another value exists for the break table.
    keep_exp_term : bool, optional
        Drop the second (exponential) term when False; for the break table.
    stable : bool, optional
        Evaluate ``exp(x) - 1`` with ``expm1`` (True) or literally.

    Returns
    -------
    numpy.ndarray
        ``D / (U d)``; the axial Peclet group is its reciprocal.
    """
    Re, Sc = np.asarray(Re, float), np.asarray(Sc, float)
    tau = G[SHAPES[shape][1]]
    p = p_of_Re(Re, G, shape) if p is None else np.asarray(p, float)
    A = Re * Sc / (4.0 * alpha1 ** 2 * (1.0 - eps))
    arg = -1.0 / (A * p * (1.0 - p))
    em1 = np.expm1(arg) if stable else np.exp(arg) - 1.0
    term2 = A ** 2 * p * (1.0 - p) ** 3 * em1
    return (A * (1.0 - p) ** 2 + (term2 if keep_exp_term else 0.0)
            + eps / (tau * Re * Sc))


def dD_over_Ud_dRe(Re, Sc, G, eps, shape="spheres", alpha1=ALPHA1, keep_exp_term=True,
                   **_):
    """Analytic derivative of eq. (1) in Re at fixed Sc, ``p = p(Re)``.

    Parameters are those of :func:`D_over_Ud`; ``stable`` is accepted and
    ignored (the derivative always uses ``expm1``).

    Returns
    -------
    numpy.ndarray
        ``d(D/(U d))/dRe``.
    """
    Re, Sc = np.asarray(Re, float), np.asarray(Sc, float)
    tau = G[SHAPES[shape][1]]
    p, dp = p_of_Re(Re, G, shape), _dp_dRe(Re, G, shape)
    A = Re * Sc / (4.0 * alpha1 ** 2 * (1.0 - eps))
    dA = A / Re
    q, dq = p * (1.0 - p), (1.0 - 2.0 * p) * dp
    y = -1.0 / (A * q)
    dy = (dA * q + A * dq) / (A * q) ** 2
    pref = A ** 2 * p * (1.0 - p) ** 3
    dpref = (2.0 * A * dA * p * (1.0 - p) ** 3
             + A ** 2 * dp * ((1.0 - p) ** 3 - 3.0 * p * (1.0 - p) ** 2))
    d1 = dA * (1.0 - p) ** 2 - 2.0 * A * (1.0 - p) * dp
    d2 = (dpref * np.expm1(y) + pref * np.exp(y) * dy) if keep_exp_term else 0.0
    d3 = -eps / (tau * Re ** 2 * Sc)
    return d1 + d2 + d3


def D_over_Ud_1993(Re, Sc, G, eps, shape="spheres", alpha1=ALPHA1):
    """The same model as Gunn (1993) prints it: eq. (42), ``Gamma = 4(1-eps)alpha_1^2/eps``.

    Parameters are those of :func:`D_over_Ud`.

    Returns
    -------
    numpy.ndarray
        ``D / (U d)``.
    """
    Re, Sc = np.asarray(Re, float), np.asarray(Sc, float)
    tau = G[SHAPES[shape][1]]
    p = p_of_Re(Re, G, shape)
    gam = 4.0 * (1.0 - eps) * alpha1 ** 2 / eps
    return ((Re * Sc / (eps * gam)) * (1.0 - p) ** 2
            + (Re ** 2 * Sc ** 2 / (eps ** 2 * gam ** 2)) * p * (1.0 - p) ** 3
            * np.expm1(-eps * gam / (p * (1.0 - p) * Re * Sc))
            + eps / (tau * Re * Sc))


def Pe_axial(Re, Sc, G, eps, **kw):
    """Axial Peclet group ``U d / D`` from eq. (1); keywords as :func:`D_over_Ud`."""
    return 1.0 / D_over_Ud(Re, Sc, G, eps, **kw)


def Pe_f(Re, G, shape="spheres"):
    """Eqs. (10) to (12): ``Pe_f = asymptote - amplitude exp(-scale / Re)``.

    Parameters
    ----------
    Re : float or numpy.ndarray
        Reynolds number.
    G : dict
        Printed constants.
    shape : str, optional
        Particle shape.

    Returns
    -------
    numpy.ndarray
        The fluid-mechanical radial Peclet group.
    """
    keys = SHAPES[shape]
    return G[keys[2]] - G[keys[3]] * np.exp(-G["pef_reynolds_scale"] / np.asarray(Re, float))


def Pe_radial(Re, Sc, G, eps, shape="spheres"):
    """Eq. (9): ``1/Pe = 1/Pe_f + eps / (tau Re Sc)``, with the radial tortuosity.

    Parameters
    ----------
    Re, Sc : float or numpy.ndarray
        Reynolds and Schmidt numbers.
    G : dict
        Printed constants.
    eps : float
        Bed porosity.
    shape : str, optional
        Particle shape.

    Returns
    -------
    numpy.ndarray
        Radial Peclet group ``U d / D_R``.
    """
    tau = G[SHAPES[shape][4]]
    Re, Sc = np.asarray(Re, float), np.asarray(Sc, float)
    return 1.0 / (1.0 / Pe_f(Re, G, shape) + eps / (tau * Re * Sc))


class Bed:
    """Steady 2-D axisymmetric fixed bed with a first-order reaction.

    ``d(u(r) c)/dz = D d2c/dz2 + (1/r) d/dr (r D_R dc/dr) - k c`` on
    ``0 <= z <= L``, ``0 <= r <= R``, lengths in particle diameters. The state
    has the shape ``(n_z, n_r)``, spatial axes first. All operators are
    assembled once in the constructor; :meth:`solve` only solves.

    Parameters
    ----------
    D, D_R : float
        Axial and radial dispersion coefficients.
    k : float
        First-order rate constant.
    U : float, optional
        Area-mean interstitial velocity.
    d : float, optional
        Particle diameter (the length unit).
    L : float, optional
        Bed length.
    d_t : float, optional
        Tube diameter in particle diameters.
    n_z, n_r : int, optional
        Equal cells along ``z`` and ``r``.
    wall_ratio : float, optional
        Interstitial velocity ratio ``U_w / U_b`` in an annulus ``wall_extent``
        diameters thick. The profile is renormalised so that its area mean is
        ``U`` whatever the ratio: the wall comparison redistributes the feed.
    wall_extent : float, optional
        Thickness of the wall annulus in particle diameters.
    outlet : {"gunn", "zero-gradient"}, optional
        Gunn's eq. (18), ``dc/dz = -lam c``, or Danckwerts' ``dc/dz = 0``.
    inlet : {"danckwerts", "dirichlet"}, optional
        The Danckwerts flux balance or ``c = c_0`` at ``z = 0``.
    nu_r : int, optional
        Geometry of the radial divergence: 1 cylindrical (the model); 0 exists
        for the break table.
    inlet_mode : {0, 1}, optional
        Feed profile: 0 flat, 1 the first zero-flux Bessel mode
        ``J_0(beta_1 r / R)``.
    D_disc, DR_disc : float, optional
        Coefficients used in the discretised operator only, while ``lam``,
        ``k_eff`` and :meth:`closed_form` keep ``D`` and ``D_R``: an
        operator/reference mismatch, for the break table.
    local_inlet : bool, optional
        Write the Danckwerts inlet with the local velocity ``u(r)`` instead of
        the area mean ``U``.
    lam_scale : float, optional
        Multiplies ``lambda`` in the eq. (18) dict only; values other than 1
        exist for the break table.

    Attributes
    ----------
    lam : float
        Gunn's ``lambda = sqrt(U^2/(4 D^2) + k_eff/D) - U/(2 D)``.
    k_eff : float
        ``k`` for a flat feed, ``k + D_R (beta_1/R)^2`` for the Bessel mode.
    c : numpy.ndarray
        Solution at the cell centres after :meth:`solve`, shape ``(n_z, n_r)``.
    nit : int
        Newton iterations or deferred-correction sweeps of the last solve.
    """

    def __init__(self, D, D_R, k, U=1.0, d=1.0, L=100.0, d_t=20.0, n_z=400, n_r=20,
                 wall_ratio=1.0, wall_extent=0.5, outlet="gunn", inlet="danckwerts",
                 nu_r=1, inlet_mode=0, D_disc=None, DR_disc=None, local_inlet=False,
                 lam_scale=1.0):
        self.D, self.DR, self.k, self.U, self.L, self.d = D, D_R, k, U, L, d
        D_d = D if D_disc is None else D_disc            # the operator's coefficients
        DR_d = D_R if DR_disc is None else DR_disc
        self.R = R = 0.5 * d_t * d
        self.z_f = np.linspace(0.0, L, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.r_f = np.linspace(0.0, R, n_r + 1)
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        self.shape = (n_z, n_r)

        self.w = (self.r_f[1:] ** 2 - self.r_f[:-1] ** 2) / R ** 2   # exact annulus areas
        prof = np.where(self.r_c > R - wall_extent * d, wall_ratio, 1.0)
        self.prof = prof / np.sum(self.w * prof)
        self.v = np.broadcast_to(U * self.prof, (n_z + 1, n_r)).copy()   # face velocity

        self.shape_in = np.ones(n_r) if inlet_mode == 0 else j0(BETA1 * self.r_c / R)
        self.k_eff = k + (D_R * (BETA1 / R) ** 2 if inlet_mode == 1 else 0.0)
        self.lam = (np.sqrt(U ** 2 / (4 * D ** 2) + self.k_eff / D) - U / (2 * D))

        c0 = 1.0
        # z = 0, outward normal -z.  Danckwerts: u c0 = u c - D dc/dz, i.e.
        # D dc/dn + u c = u c0, with u the area mean U (or the local u(r)).
        u_in = U * self.prof if local_inlet else np.full(n_r, U)
        bc_in = ({"a": D_d, "b": u_in, "d": u_in * c0 * self.shape_in}
                 if inlet == "danckwerts"
                 else {"a": 0.0, "b": 1.0, "d": c0 * self.shape_in})   # c = c0 (break row)
        # z = L, outward normal +z.  Gunn eq. (18): dc/dz = -lam c, i.e. dc/dn + lam c = 0.
        # Danckwerts: dc/dz = 0, i.e. dc/dn = 0.
        bc_out = ({"a": 1.0, "b": lam_scale * self.lam, "d": 0.0} if outlet == "gunn"
                  else {"a": 1.0, "b": 0.0, "d": 0.0})
        self.bc_z = (bc_in, bc_out)
        # r = 0: symmetry, dc/dn = 0.  r = R: impermeable wall, dc/dn = 0.
        self.bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 0.0})

        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c, self.bc_z,
                                                  v=self.v, axis=0)
        gz, gz_bc = construct_grad(self.shape, self.z_f, self.z_c, self.bc_z, axis=0)
        self.div_z = construct_div(self.shape, self.z_f, nu=0, axis=0)   # nu=0: axial, Cartesian
        gr, gr_bc = construct_grad(self.shape, self.r_f, self.r_c, self.bc_r, axis=1)
        div_r = construct_div(self.shape, self.r_f, nu=nu_r, axis=1)     # nu=1: radial, cylindrical
        n_tot = n_z * n_r
        self.conv = conv.tocsc()
        self.conv_bc = np.asarray(conv_bc.todense()).ravel()
        self.jac_up = (self.div_z @ (conv - D_d * gz) - DR_d * (div_r @ gr)
                       + k * eye_array(n_tot, format="csc")).tocsc()
        self.g_up = np.asarray((self.div_z @ (conv_bc - D_d * gz_bc)
                                - DR_d * (div_r @ gr_bc)).todense()).ravel()
        self.jac_no_conv = (self.jac_up - self.div_z @ self.conv).tocsc()
        self.lu_up = splu(self.jac_up)            # constant operator: factorised once
        self.tvd = ConvFluxTVD(self.shape, self.z_f, self.z_c, self.bc_z,
                               tvd_limiter=vanleer, axis=0)
        self.c = None
        self.nit = 0

    def solve(self, limiter=vanleer, route="newton", tol=1e-12, maxfev=60, max_sweeps=300,
              jacobian="exact"):
        """Solve the steady bed.

        Parameters
        ----------
        limiter : callable or None, optional
            TVD limiter (van Leer by default); None gives bare first-order
            upwind, one linear solve.
        route : {"newton", "deferred"}, optional
            Newton with the exact ``ConvFluxTVD`` Jacobian, or van Leer
            deferred correction swept to a fixed point.
        tol : float, optional
            Step tolerance relative to the largest value of each axial row
            (the decaying envelope), for either route.
        maxfev, max_sweeps : int, optional
            Iteration caps.
        jacobian : {"exact", "upwind"}, optional
            Newton's Jacobian: the exact ``ConvFluxTVD`` one, or the upwind
            operator alone (a break-table row, pitfalls P10).

        Returns
        -------
        Bed
            ``self``, with ``c`` and ``nit`` set.

        Raises
        ------
        RuntimeError
            When Newton does not report success or the sweeps do not reach
            ``tol``.
        """
        c_up = self.lu_up.solve(-self.g_up).reshape(self.shape)
        if limiter is None:
            self.c, self.nit = c_up, 0
            return self
        scale = np.broadcast_to(np.max(np.abs(c_up), axis=1, keepdims=True), self.shape)
        if limiter is not vanleer:
            tvd = ConvFluxTVD(self.shape, self.z_f, self.z_c, self.bc_z,
                              tvd_limiter=limiter, axis=0)
        else:
            tvd = self.tvd
        if route == "newton":
            c_up_flat = c_up.ravel()

            def residual(x):
                c = x.reshape(-1)
                flux, jac = tvd(x.reshape(self.shape), self.v)
                flux_up = self.conv @ c + self.conv_bc
                g = self.jac_up @ (c - c_up_flat) + self.div_z @ (flux.reshape(-1) - flux_up)
                if jacobian == "upwind":
                    return g, self.jac_up
                return g, (self.jac_no_conv + self.div_z @ jac).tocsc()

            sol = newton(residual, c_up, tol=tol, maxfev=maxfev, damping="monotone",
                         stagnation=True, scale=scale)
            if not sol.success:
                raise RuntimeError(f"TVD Newton solve failed: {sol.message}")
            self.c, self.nit = np.asarray(sol.x).reshape(self.shape), int(sol.nit)
            return self
        c = c_up
        for sweep in range(1, max_sweeps + 1):
            _, dc_f = interp_cntr_to_stagg_tvd(c, self.z_f, self.z_c, self.bc_z, self.v,
                                               tvd_limiter=limiter, axis=0)
            c_new = self.lu_up.solve(-self.g_up - self.div_z @ (self.v * dc_f).ravel())
            c_new = c_new.reshape(self.shape)
            change = float(np.max(np.abs(c_new - c) / scale))
            c = c_new
            if change < tol:
                self.c, self.nit = c, sweep
                return self
        raise RuntimeError(f"deferred correction did not reach {tol:g} in {max_sweeps} "
                           f"sweeps (last change {change:.2e})")

    def closed_form(self):
        """Gunn's eq. (17) at the cell centres, with ``k -> k_eff``.

        For a flat feed this is eq. (17) unchanged; for the Bessel-mode feed
        the mode is an eigenfunction of the radial operator and the same
        expression holds with ``k + D_R (beta_1/R)^2``.
        """
        return (self.U / (self.U + self.D * self.lam)
                * np.exp(-self.lam * self.z_c)[:, None] * self.shape_in[None, :])

    def outlet_values(self, at="face"):
        """Values at the outlet ``z = L``.

        ``at="face"`` reads the outlet face with ``compute_boundary_values``
        under the bed's own outlet dict (pitfalls P8); ``at="centre"`` returns
        the last cell centres, half a cell short of the outlet, and exists for
        the break table.
        """
        if at == "centre":
            return self.c[-1].copy()
        return np.asarray(compute_boundary_values(self.c, self.z_f, self.z_c, self.bc_z,
                                                  axis=0)[2]).reshape(-1)

    def outlet_amplitude(self, at="face"):
        """Projection of the outlet values on the feed profile.

        ``sum(w c_L s) / sum(w s^2)`` with ``s`` the feed profile at the cell
        centres and ``w`` the exact annulus areas: the area mean for a flat
        feed, the amplitude of the Bessel mode otherwise. The closed form gives
        ``U / (U + D lam) exp(-lam L)``.
        """
        s = self.shape_in
        return float(np.sum(self.w * self.outlet_values(at) * s) / np.sum(self.w * s * s))

    def cup(self, at="face"):
        """Mixing-cup concentration at the outlet, weighted by the local flux."""
        flux = self.w * self.prof
        return float(np.sum(self.outlet_values(at) * flux) / np.sum(flux))
