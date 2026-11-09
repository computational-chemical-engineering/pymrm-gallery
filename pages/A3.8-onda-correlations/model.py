"""Models of page A3.8 (Onda, Takeuchi and Okumoto 1968).

* :class:`OndaSet`, Onda's eqs. (1) to (3) and the ``K_G a`` they imply, with every
  printed constant a named argument;
* :class:`Absorber`, the counter-current packed absorber as a pymrm boundary value
  problem (convection of two fields with opposite velocities and a pointwise
  linear interphase source; van Leer by Newton on the exact Jacobian, or by
  deferred correction);
* :func:`column_closed_form`, the exact solution of the same boundary value problem.

Onda's units throughout: m, kg, hr, kg-mole, atm, K. The printed constants are
read from the page's CSV, not typed here: the page passes them to
:func:`set_printed_constants` once, and :class:`OndaSet` falls back to them for
every constant it is not given.
"""

import numpy as np
from scipy.sparse.linalg import splu
from pymrm import (ConvFluxTVD, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, generate_grid,
                   interp_cntr_to_stagg_tvd, newton, vanleer)

PRINTED = {}   # the fourteen printed constants, filled by set_printed_constants

_NAMES = ("C_aw", "e_sigma", "e_Re_L", "e_Fr_L", "e_We_L", "C_kL", "e_Re_Lw", "e_Sc_L",
          "e_atDp_L", "C_kG", "e_Re_G", "e_Sc_G", "e_atDp_G")


def set_printed_constants(**constants):
    """Set the printed constants that :class:`OndaSet` uses as defaults.

    Parameters
    ----------
    **constants : float
        Any of ``C_aw, e_sigma, e_Re_L, e_Fr_L, e_We_L, C_kL, e_Re_Lw, e_Sc_L,
        e_atDp_L, C_kG, e_Re_G, e_Sc_G, e_atDp_G``, as read from the CSV.
    """
    unknown = set(constants) - set(_NAMES)
    if unknown:
        raise KeyError(f"not a printed constant: {sorted(unknown)}")
    PRINTED.update(constants)


HR = 3600.0                     # s per hr, for converting tabulated SI properties
R_GAS = 0.0820574               # m3 atm / (kg-mole K)   -- Onda's R
G_ACC = 9.80665 * HR**2         # m/hr2                  -- Onda's g


class OndaSet:
    """Onda's eqs. (1), (2) and (3), and the K_G a they imply.

    Every printed constant is a named argument defaulting to the CSV value, so a
    single one can be perturbed without touching anything downstream. That is
    what the break table and the elasticity study both use.

    Parameters
    ----------
    a_t : float
        Total packing area per bed volume, m2/m3.
    D_p : float
        Packing size, m.
    rho_L, mu_L, D_L, sigma : float
        Liquid density (kg/m3), viscosity (kg/(m hr)), diffusivity (m2/hr) and
        surface tension (kg/hr2).
    sigma_ratio : float
        Critical surface tension of the packing over ``sigma``.
    rho_G, mu_G, D_G, T : float
        Gas density (kg/m3), viscosity (kg/(m hr)), diffusivity (m2/hr) and
        temperature (K).
    P : float, optional
        Total pressure, atm.
    C_aw, ..., e_atDp_G : float, optional
        The thirteen continuous printed constants of eqs. (1) to (3); each defaults
        to the value given to :func:`set_printed_constants`.
    """

    def __init__(self, *, a_t, D_p, rho_L, mu_L, D_L, sigma, sigma_ratio,
                 rho_G, mu_G, D_G, T, P=1.0,
                 C_aw=None, e_sigma=None, e_Re_L=None, e_Fr_L=None, e_We_L=None,
                 C_kL=None, e_Re_Lw=None, e_Sc_L=None, e_atDp_L=None,
                 C_kG=None, e_Re_G=None, e_Sc_G=None, e_atDp_G=None):
        given = dict(C_aw=C_aw, e_sigma=e_sigma, e_Re_L=e_Re_L, e_Fr_L=e_Fr_L,
                     e_We_L=e_We_L, C_kL=C_kL, e_Re_Lw=e_Re_Lw, e_Sc_L=e_Sc_L,
                     e_atDp_L=e_atDp_L, C_kG=C_kG, e_Re_G=e_Re_G, e_Sc_G=e_Sc_G,
                     e_atDp_G=e_atDp_G)
        self.__dict__.update(a_t=a_t, D_p=D_p, rho_L=rho_L, mu_L=mu_L, D_L=D_L,
                             sigma=sigma, sigma_ratio=sigma_ratio, rho_G=rho_G,
                             mu_G=mu_G, D_G=D_G, T=T, P=P)
        # a constant not given is the printed one, from set_printed_constants
        self.__dict__.update({k: PRINTED[k] if v is None else v for k, v in given.items()})
        self.Sc_L = mu_L / (rho_L * D_L)
        self.Sc_G = mu_G / (rho_G * D_G)
        self.atDp = a_t * D_p

    # ---- eq. (1) -----------------------------------------------------------
    def psi(self, L):
        """The exponent argument of eq. (1); a_w/a_t = 1 - exp(-psi)."""
        Re = L / (self.a_t * self.mu_L)
        Fr = L**2 * self.a_t / (self.rho_L**2 * G_ACC)
        We = L**2 / (self.rho_L * self.sigma * self.a_t)
        return (self.C_aw * self.sigma_ratio**self.e_sigma * Re**self.e_Re_L
                * Fr**self.e_Fr_L * We**self.e_We_L)

    def a_w(self, L):
        """Wetted area per bed volume, ``a_t (1 - exp(-psi))``, m2/m3."""
        return self.a_t * (1.0 - np.exp(-self.psi(L)))

    def phi(self, L):
        """dln(a_w)/dln(psi) = psi e^-psi / (1 - e^-psi); -> 1 as psi -> 0, 0 as psi -> inf."""
        p = self.psi(L)
        return p * np.exp(-p) / (1.0 - np.exp(-p))

    # ---- eq. (2) -----------------------------------------------------------
    def k_L(self, L):
        """Liquid-side coefficient of eq. (2) at liquid mass velocity ``L``, m/hr."""
        aw = self.a_w(L)
        return (self.C_kL * (L / (aw * self.mu_L))**self.e_Re_Lw
                * self.Sc_L**self.e_Sc_L * self.atDp**self.e_atDp_L
                / (self.rho_L / (self.mu_L * G_ACC))**(1.0 / 3.0))

    # ---- eq. (3) -----------------------------------------------------------
    def k_G(self, G):
        """Gas-side coefficient of eq. (3) at gas mass velocity ``G``, kg-mole/(m2 hr atm)."""
        return (self.C_kG * (G / (self.a_t * self.mu_G))**self.e_Re_G
                * self.Sc_G**self.e_Sc_G * self.atDp**self.e_atDp_G
                * self.a_t * self.D_G / (R_GAS * self.T))

    # ---- resistances in series (A3.1), with a = a_w (Onda's assumption) -----
    def resistances(self, L, G, m, c_av):
        """Gas-side and liquid-side resistances ``(R_G, R_L)`` per unit mole fraction.

        ``m`` is the equilibrium slope ``y* = m x`` and ``c_av`` the liquid molar
        density, kg-mole/m3.
        """
        aw = self.a_w(L)
        R_G = 1.0 / (self.P * self.k_G(G) * aw)      # gas-side, per unit mole fraction
        R_L = m / (self.k_L(L) * aw * c_av)          # liquid-side, same units
        return R_G, R_L

    def KGa(self, L, G, m, c_av):
        """Overall gas-side coefficient ``K_G a``, kg-mole/(m3 hr) per unit mole fraction."""
        R_G, R_L = self.resistances(L, G, m, c_av)
        return 1.0 / (R_G + R_L)

    def f_G(self, L, G, m, c_av):
        """Gas-side share ``R_G / (R_G + R_L)`` of the total resistance."""
        R_G, R_L = self.resistances(L, G, m, c_av)
        return R_G / (R_G + R_L)


class Absorber:
    """Counter-current packed absorber: Onda's eq. (6) as a boundary value problem.

    Fields ``[y, x]`` on one axial grid, ``z`` from the bottom of the packing. Gas
    rises (``+G_M``), liquid falls (``-L_M``), so the two phases are two components of
    one array with opposite velocities and inlets at opposite ends:

        d(G_M y)/dz + N = 0,   d(-L_M x)/dz - N = 0,   N = K_G a (y - m x).

    The source is linear in ``(y, x)``, so with upwind convection the discrete problem
    is one constant sparse matrix and one solve. The van Leer limiter makes the face
    values depend on the field, so the second-order problem is nonlinear: it is
    solved by ``newton`` on :class:`pymrm.ConvFluxTVD`, which supplies the exact
    Jacobian of the limited fluxes, or by deferred correction on the factorised
    upwind matrix (the second route).

    Parameters
    ----------
    Z : float
        Packed height, m.
    KGa : float
        Overall coefficient, kg-mole/(m3 hr) per unit mole fraction.
    m : float
        Equilibrium slope, ``y* = m x``.
    G_M, L_M : float
        Molar mass velocities of gas and liquid, kg-mole/(m2 hr).
    y_in, x_in : float
        Inlet mole fractions of gas (at ``z = 0``) and liquid (at ``z = Z``).
    n : int, optional
        Number of axial cells.
    limiter : callable, optional
        TVD limiter of the second-order faces, a pymrm limiter function.
    cocurrent : bool, optional
        Liquid fed at the bottom and flowing up with the gas instead.
    """

    def __init__(self, Z, KGa, m, G_M, L_M, y_in, x_in, n=200, limiter=vanleer,
                 cocurrent=False):
        self.n, self.shape = n, (n, 2)
        self.z_f, self.z_c = generate_grid(n, [0.0, Z], generate_x_c=True)
        self.Z, self.KGa, self.m, self.G_M, self.L_M = Z, KGa, m, G_M, L_M
        self.y_in, self.x_in = y_in, x_in
        self.v = np.array([[G_M, -L_M]])          # gas up, liquid down
        # Outward normal n, a dphi/dn + b phi = d, one entry per field [y, x]:
        #   z = 0 (n = -z): gas inlet, y = y_in                 -> a = 0, b = 1, d = y_in
        #                   liquid outlet, dx/dz = 0            -> a = 1, b = 0, d = 0
        #   z = Z (n = +z): gas outlet, dy/dz = 0               -> a = 1, b = 0, d = 0
        #                   liquid inlet, x = x_in              -> a = 0, b = 1, d = x_in
        # {"outflow": True} marks a whole face; each face here is the inlet of one
        # field and the outlet of the other, so the outlets are the per-field
        # zero-gradient entries, and outlet values are read with
        # compute_boundary_values under this same bc.
        self.bc = ({"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[y_in, 0.0]]},
                   {"a": [[1.0, 0.0]], "b": [[0.0, 1.0]], "d": [[0.0, x_in]]})
        self.cocurrent = cocurrent
        if cocurrent:
            # both phases enter at z = 0 (y = y_in, x = x_in) and leave at z = Z,
            # a pure-outflow face for both fields
            self.v = np.array([[G_M, L_M]])
            self.bc = ({"a": 0.0, "b": 1.0, "d": [[y_in, x_in]]}, {"outflow": True})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                 self.bc, v=self.v, axis=0)
        self.div = construct_div(self.shape, self.z_f, nu=0, axis=0)   # nu=0: Cartesian axis
        # per-cell coupling of the residual rows [y, x] to the fields [y, x]: +N, -N
        self.source = construct_coefficient_matrix(
            np.array([[KGa, -KGa * m], [-KGa, KGa * m]]), shape=((n, 2, 1), (n, 1, 2)))
        self.upwind = (self.div @ conv + self.source).tocsc()
        self.g_upwind = (self.div @ conv_bc).toarray().ravel()
        self.tvd = ConvFluxTVD(self.shape, self.z_f, self.z_c, bc=self.bc,
                               tvd_limiter=limiter, axis=0)
        self.limiter = limiter
        self.u, self.its = None, 0

    def residual(self, u):
        """Residual and exact Jacobian of the van Leer balance at the field ``u``."""
        flux, jac = self.tvd(u, self.v)
        return (self.div @ flux.ravel() + self.source @ np.ravel(u),
                self.div @ jac + self.source)

    def solve(self, scheme="vanleer", route="newton", rtol=1e-12, max_it=200):
        """Solve the boundary value problem.

        Parameters
        ----------
        scheme : {"vanleer", "upwind"}, optional
            Face values: van Leer limited (second order), or first-order upwind,
            which is linear and takes one sparse solve.
        route : {"newton", "deferred"}, optional
            For ``"vanleer"``: ``newton`` on the exact Jacobian, or deferred
            correction (the limiter correction lagged, the upwind matrix factorised
            once).
        rtol : float, optional
            Relative step that ends either iteration.
        max_it : int, optional
            Cap on deferred-correction iterations; reaching it raises.

        Returns
        -------
        Absorber
            ``self``, with ``u`` the solved field and ``its`` the iteration count.
        """
        lu = splu(self.upwind)
        u = lu.solve(-self.g_upwind).reshape(self.shape)    # the upwind solution
        self.its = 0
        if scheme == "upwind":
            self.u = u
            return self
        if route == "newton":
            floor = 1e3 * np.finfo(float).eps * np.abs(u).max()      # above the round-off of u
            res = newton(self.residual, u, damping="monotone", tol=floor, rtol=rtol,
                         stagnation=True)
            if not res.success:
                raise RuntimeError(f"newton failed on the van Leer column: {res.message}")
            self.u, self.its = res.x.reshape(self.shape), int(res.nit)
            return self
        for it in range(1, max_it + 1):
            _, dc_f = interp_cntr_to_stagg_tvd(u, self.z_f, self.z_c, self.bc, self.v,
                                              tvd_limiter=self.limiter, axis=0)
            corr = self.div @ (self.v * dc_f).ravel()
            u_new = lu.solve(-self.g_upwind - corr).reshape(self.shape)
            step = np.max(np.abs(u_new - u))
            u = u_new
            if step <= rtol * np.abs(u).max():
                self.u, self.its = u, it
                return self
        raise RuntimeError(f"deferred correction did not converge in {max_it} iterations")

    def faces(self):
        """Outlet values at the faces, from ``compute_boundary_values`` under ``self.bc``.

        The gas leaves through the face at ``z = Z`` and the liquid through the face
        at ``z = 0``; the last cell centres sit half a cell short of them.
        """
        v_lo, _, v_hi, _ = compute_boundary_values(self.u, self.z_f, self.z_c,
                                                   self.bc, axis=0)
        x_out = v_hi[0, 1] if self.cocurrent else v_lo[0, 1]
        return float(v_hi[0, 0]), float(x_out)           # y at z=Z, x at its outlet

    @property
    def y_out(self):
        """Gas leaving the top, at the face."""
        return self.faces()[0]

    @property
    def x_out(self):
        """Liquid leaving the bottom, at the face."""
        return self.faces()[1]

    @property
    def y_out_cell(self):
        """Gas at the centre of the top cell, half a cell short of the outlet face."""
        return float(self.u[-1, 0])

    def imbalance(self, cell=False):
        """Solute lost by the gas against solute gained by the liquid, relative.

        With the face values this is a telescoping sum of the discrete fluxes and
        holds to round-off for any ``K_G a``; it tests the transport bookkeeping and
        nothing about the physics. ``cell=True`` reads both outlets at the last cell
        centres instead.
        """
        x_cell = float(self.u[-1 if self.cocurrent else 0, 1])
        y_out, x_out = ((self.y_out_cell, x_cell) if cell
                        else (self.y_out, self.x_out))
        gas = self.G_M * (self.y_in - y_out)
        liq = self.L_M * (x_out - self.x_in)
        return abs(gas / liq - 1.0)


def column_closed_form(Z, KGa, m, G_M, L_M, y_in, x_in):
    """Exact solution of the same BVP, by hand. It shares no operator with :class:`Absorber`.

    D = y - m x obeys dD/dz = -lambda D with lambda = KGa (1/G_M - m/L_M), so the
    whole profile is one exponential and the two-point problem closes analytically.

    Parameters
    ----------
    Z, KGa, m, G_M, L_M, y_in, x_in : float
        As for :class:`Absorber`.

    Returns
    -------
    tuple of float
        Gas mole fraction at ``z = Z`` and liquid mole fraction at ``z = 0``.
    """
    lam = KGa * (1.0 / G_M - m / L_M)
    E = np.exp(-lam * Z)
    beta = (m * KGa / L_M) / lam
    D0 = (y_in - m * x_in) / (1.0 + beta * (1.0 - E))
    y_out = y_in - (KGa / G_M) * D0 * (1.0 - E) / lam
    x_bot = x_in + (KGa / L_M) * D0 * (1.0 - E) / lam
    return float(y_out), float(x_bot)
