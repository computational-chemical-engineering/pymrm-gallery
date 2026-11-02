"""The Mears criteria: the pellet, film and axial-dispersion problems behind them.

Model code of gallery page B1.7 (Mears 1971). Every criterion of the paper is
an inequality in observables; the page measures what each admits by solving
the underlying problem exactly. Three routes to the pellet and two to the
reactor:

``shoot`` and ``cross``
    the pellet integrated outward from a prescribed centre concentration, so a
    whole solution branch is traversed deterministically and the state where
    a criterion binds is root-found on it;
``Pellet``
    the same pellet by pymrm finite volume, Dirichlet or film surface;
``eta_iso1``
    the closed-form isothermal first-order sphere;
``AxialVessel`` and ``frac_unconverted``
    the closed dispersion vessel by pymrm (van Leer TVD, exact Jacobian) and
    the classical first-order closed form.

Contents
--------
Power, WeiszHicks, InhibitedLH
    Rate laws ``R(y)`` normalised to ``R(1) = 1``, with their derivative.
shoot, cross
    The shooting reference and the root-found states where a condition binds.
eta_iso1
    Closed-form isothermal first-order effectiveness factor.
Pellet
    ``nabla^2 y = phi^2 R(y)`` by pymrm, exact Jacobian.
film_mass_dev, film_heat_dev
    The interphase deviations, closed functions of the observables.
AxialVessel, frac_unconverted, axial_dev, axial_dev_covary, petersen_dev
    The axial-dispersion criteria's reactor and its required-length deficit.
joint_worst
    Worst first-order deviation among states that pass eqs. 1, 17 and 21.

Several functions take an argument whose default is the correct model and
whose alternative is a defect the page's break table injects.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu

from pymrm import (ConvFluxTVD, NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, construct_grad,
                   interp_cntr_to_stagg_tvd, newton, vanleer)
from pymrm.checks import find_roots

__all__ = ["Power", "WeiszHicks", "InhibitedLH", "shoot", "cross", "eta_iso1", "SolveFailed",
           "Pellet", "film_mass_dev", "film_heat_dev", "AxialVessel", "frac_unconverted",
           "axial_dev", "axial_dev_covary", "petersen_dev", "joint_worst"]


# ------------------------------------------------------------------ rate laws
class Power:
    """Power law ``R = y^n``; ``n = 0`` is zero order, ``n < 0`` a negative order.

    For ``n < 0`` the concentration is floored at 1e-9 so that the rate stays
    finite in a depleted cell.

    Parameters
    ----------
    n : float
        Reaction order.
    """

    def __init__(self, n):
        self.n = float(n)
        self.floor = 1e-9 if self.n < 0 else 0.0

    def r(self, y):
        """Normalised rate.

        Parameters
        ----------
        y : array_like
            Dimensionless concentration ``C/C_s``.

        Returns
        -------
        numpy.ndarray
            Rate ``R(y)``, same shape as ``y``.
        """
        y = np.clip(np.asarray(y, float), self.floor, 1.0)
        return np.ones_like(y) if self.n == 0.0 else y ** self.n

    def dr(self, y):
        """Derivative ``dR/dy``.

        Parameters
        ----------
        y : array_like
            Dimensionless concentration ``C/C_s``.

        Returns
        -------
        numpy.ndarray
            ``dR/dy``, same shape as ``y``.
        """
        y = np.clip(np.asarray(y, float), self.floor, 1.0)
        return np.zeros_like(y) if self.n == 0.0 else self.n * y ** (self.n - 1.0)


class WeiszHicks:
    """Power-law kinetics with the Prater-reduced Arrhenius factor.

    ``R = y^n exp[gb (1 - y)/(1 + beta (1 - y))]`` with ``beta = gb/gamma``,
    parametrised by the product ``gb = gamma beta``.

    Parameters
    ----------
    n : float
        Reaction order.
    gb : float
        ``gamma beta``; positive for an exothermic reaction.
    gamma : float
        Arrhenius number ``E/(R_g T_s)``.
    """

    def __init__(self, n, gb, gamma):
        self.n, self.gb, self.gamma = float(n), float(gb), float(gamma)
        self.beta = self.gb / self.gamma

    def r(self, y):
        """Normalised rate.

        Parameters
        ----------
        y : array_like
            Dimensionless concentration ``C/C_s``.

        Returns
        -------
        numpy.ndarray
            Rate ``R(y)``, same shape as ``y``.
        """
        y = np.asarray(y, float)
        yc = np.clip(y, 0.0, 1.0)
        return yc ** self.n * np.exp(self.gb * (1 - y) / (1 + self.beta * (1 - y)))

    def dr(self, y):
        """Derivative ``dR/dy``.

        Parameters
        ----------
        y : array_like
            Dimensionless concentration ``C/C_s``.

        Returns
        -------
        numpy.ndarray
            ``dR/dy``, same shape as ``y``.
        """
        y = np.asarray(y, float)
        yc = np.clip(y, 0.0, 1.0)
        w = 1 + self.beta * (1 - y)
        e = np.exp(self.gb * (1 - y) / w)
        return e * (self.n * yc ** (self.n - 1.0) - yc ** self.n * self.gb / w ** 2)


class InhibitedLH:
    """Self-inhibited Langmuir-Hinshelwood rate ``R = (1 + K)^2 y/(1 + K y)^2``.

    Parameters
    ----------
    K : float, optional
        Adsorption constant times surface concentration.
    """

    def __init__(self, K=5.0):
        self.K = float(K)

    def r(self, y):
        """Normalised rate.

        Parameters
        ----------
        y : array_like
            Dimensionless concentration ``C/C_s``.

        Returns
        -------
        numpy.ndarray
            Rate ``R(y)``, same shape as ``y``.
        """
        y = np.clip(np.asarray(y, float), 0.0, 1.0)
        return (1.0 + self.K) ** 2 * y / (1.0 + self.K * y) ** 2

    def dr(self, y):
        """Derivative ``dR/dy``.

        Parameters
        ----------
        y : array_like
            Dimensionless concentration ``C/C_s``.

        Returns
        -------
        numpy.ndarray
            ``dR/dy``, same shape as ``y``.
        """
        y = np.clip(np.asarray(y, float), 0.0, 1.0)
        return (1.0 + self.K) ** 2 * (1.0 - self.K * y) / (1.0 + self.K * y) ** 3


# ------------------------------------------------------------- shooting route
def shoot(rate, t, nu=2, rtol=1e-10, bim=None, dense=False):
    """One branch state from its centre concentration ``y_c = 1 - exp(-t)``.

    Integrates ``y'' + (nu/s) y' = R(y)`` in ``s = phi u`` outward from the
    centre. With a Dirichlet surface (``bim=None``) it stops where ``y = 1``,
    and that radius is ``phi``; with a film surface it stops where
    ``s y'(s) = Bi_m (1 - y)`` and ``y`` is scaled on the bulk value. No grid,
    no operator, no Newton iteration: it shares only the rate law with
    :class:`Pellet`.

    Parameters
    ----------
    rate : Power, WeiszHicks or InhibitedLH
        Rate law.
    t : float
        Branch parameter, ``t > 0``; small ``t`` is the deep-diffusion end.
    nu : int, optional
        Geometry index: 0 slab, 1 cylinder, 2 sphere.
    rtol : float, optional
        Relative tolerance of the DOP853 integration.
    bim : float or None, optional
        Mass Biot number ``k_c r_p/D_e`` of a film surface; None is Dirichlet.
    dense : bool, optional
        Also return the ``solve_ivp`` result with dense output.

    Returns
    -------
    tuple or None
        ``(phi, eta, Phi, omega, y_s)``: Thiele modulus, effectiveness factor
        and observable ``Phi = phi^2 eta`` (on the surface or, with a film, the
        bulk normalisation), ``omega = Phi/Bi_m`` (None without a film) and
        the surface value; plus the integration result when ``dense``. None
        when the surface condition is never met.
    """
    y_c = -np.expm1(-t)
    r0 = float(rate.r(y_c))
    if not (0.0 < y_c < 1.0) or r0 <= 0.0:
        return None
    # start on the series y = y_c + R s^2/(2 (nu + 1)), offset 1e-10 of the
    # distance from y_c to the nearer of 0 and 1
    s0 = min(1e-3, np.sqrt(2 * (nu + 1) * 1e-10 * min(y_c, np.exp(-t)) / abs(r0)))

    def rhs(s, v):
        return [v[1], float(rate.r(min(max(v[0], 0.0), 1.0))) - (nu / s) * v[1]]

    if bim is None:
        def hit(s, v):
            return v[0] - 1.0
    else:
        def hit(s, v):
            return s * v[1] - bim * (1.0 - v[0])
    hit.terminal, hit.direction = True, 1
    v0 = [y_c + r0 * s0 ** 2 / (2 * (nu + 1)), r0 * s0 / (nu + 1)]
    sol = solve_ivp(rhs, (s0, 1e5), v0, events=hit, rtol=rtol, atol=1e-16,
                    method="DOP853", dense_output=dense)
    if not sol.t_events[0].size:
        return None
    s = float(sol.t_events[0][0])
    y_s, yp = (float(x) for x in sol.y_events[0][0])
    out = (s, (nu + 1) * yp / s, (nu + 1) * s * yp,
           None if bim is None else (nu + 1) * s * yp / bim, y_s)
    return out + (sol,) if dense else out


def cross(rate, val, lo=1e-6, hi=20.0, n_scan=90, nu=2, bim=None, rtol=1e-10, which="deep"):
    """States on the branch where ``val(phi, eta, Phi, omega, y_s)`` changes sign.

    A :func:`pymrm.checks.find_roots` scan on a fixed logarithmic ladder of the
    branch parameter ``t`` (no continuation), each sign change refined by
    Brent's method on :func:`shoot` itself to 1e-13 in ``t``.

    Parameters
    ----------
    rate : Power, WeiszHicks or InhibitedLH
        Rate law.
    val : callable
        Condition; its zeros are returned.
    lo, hi : float, optional
        Range of ``t``.
    n_scan : int, optional
        Scan points.
    nu, bim, rtol
        As for :func:`shoot`.
    which : {"deep", "kinetic", "all"}, optional
        ``"deep"`` returns the root with the smallest ``t``, the largest-``Phi``
        state (the deep-side entry into the condition's band);
        ``"kinetic"`` the one with the largest ``t``; ``"all"`` every root,
        deep end first.

    Returns
    -------
    tuple, list of tuple or None
        ``(t, phi, eta, Phi, omega, y_s)`` per root; None when there is none.
    """
    def f(t):
        o = shoot(rate, t, nu=nu, bim=bim, rtol=rtol)
        return np.nan if o is None else val(*o)

    roots = find_roots(f, lo, hi, n_scan=n_scan, log=True, xtol=1e-13)["roots"]
    out = [(t,) + shoot(rate, t, nu=nu, bim=bim, rtol=rtol) for t in roots]
    if which == "all":
        return out
    if not out:
        return None
    return out[0] if which == "deep" else out[-1]


# ----------------------------------------------------------------- closed form
def eta_iso1(phi, nu=2):
    """Closed-form isothermal first-order effectiveness factor.

    Slab ``tanh(phi)/phi``; sphere ``3 (coth(phi) - 1/phi)/phi``, with its
    Taylor series below ``phi = 0.1`` where the closed form cancels.

    Parameters
    ----------
    phi : float
        Thiele modulus.
    nu : {0, 2}, optional
        Geometry index: 0 slab, 2 sphere.

    Returns
    -------
    float
        Effectiveness factor.
    """
    phi = float(phi)
    if nu == 0:
        return np.tanh(phi) / phi
    if phi < 0.1:
        x = phi * phi
        return (1.0 - x / 15 + 2 * x ** 2 / 315 - x ** 3 / 1575 + 2 * x ** 4 / 31185
                - 1382 * x ** 5 / 212837625)
    return (3.0 / phi) * (1.0 / np.tanh(phi) - 1.0 / phi)


# ------------------------------------------------------------------ pymrm pellet
class SolveFailed(RuntimeError):
    """Raised when ``newton`` reports no success."""


class Pellet:
    """``nabla^2 y = phi^2 R(y)`` on ``0 < u < 1`` by pymrm finite volume.

    Symmetry at the centre; ``y = 1`` at the surface, or a film. Layout
    ``(n_u, 1)``: the spatial axis first, the single field last. The operators
    are assembled once.

    Parameters
    ----------
    n_u : int, optional
        Number of cells (uniform).
    nu : int, optional
        Geometry index for ``construct_div`` and the observable: 0 slab,
        1 cylinder, 2 sphere.
    bim : float or None, optional
        Mass Biot number of a film surface; None is Dirichlet.
    """

    def __init__(self, n_u=400, nu=2, bim=None):
        self.n_u, self.nu, self.bim, self.shape = n_u, nu, bim, (n_u, 1)
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])
        # bc: a dy/dn + b y = d with n the OUTWARD normal (-u at u = 0, +u at u = 1)
        #   centre  u = 0: symmetry  dy/du = 0                     -> a = 1, b = 0, d = 0
        if bim is None:
            #   surface u = 1: Dirichlet  y = 1                    -> a = 0, b = 1, d = 1
            bc_s = {"a": 0.0, "b": 1.0, "d": 1.0}
        else:
            #   surface u = 1: film flux  dy/du = Bi_m (1 - y), i.e.
            #   dy/dn + Bi_m y = Bi_m                              -> a = 1, b = d = Bi_m
            bc_s = {"a": 1.0, "b": float(bim), "d": float(bim)}
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, bc_s)
        grad, grad_bc = construct_grad(self.shape, self.u_f, self.u_c, self.bc)
        div = construct_div(self.shape, self.u_f, nu=nu)        # 0 slab, 1 cylinder, 2 sphere
        self.lap = (div @ grad).tocsc()
        self.lap_bc = (div @ grad_bc).toarray().reshape(-1, 1)
        self.numjac = NumJac(self.shape)                        # only for jacobian="numjac"

    def residual(self, y, phi, rate, jacobian="exact"):
        """Residual ``lap y - phi^2 R(y)`` and its Jacobian.

        Parameters
        ----------
        y : numpy.ndarray
            State, any shape with ``n_u`` entries.
        phi : float
            Thiele modulus.
        rate : Power, WeiszHicks or InhibitedLH
            Rate law.
        jacobian : {"exact", "numjac"}, optional
            ``"exact"`` puts ``phi^2 dR/dy`` on the diagonal with
            ``construct_coefficient_matrix``; ``"numjac"`` differences the
            source with ``NumJac``.

        Returns
        -------
        g : numpy.ndarray
            Residual, shape ``(n_u, 1)``.
        jac : scipy.sparse matrix
            Jacobian.
        """
        y = np.asarray(y, float).reshape(self.shape)
        p2 = phi * phi
        if jacobian == "exact":
            g_src = p2 * rate.r(y)
            jac_src = construct_coefficient_matrix(p2 * rate.dr(y))
        else:
            g_src, jac_src = self.numjac(lambda v: p2 * rate.r(v), y)
        g = self.lap @ y.reshape((-1, 1)) + self.lap_bc - np.asarray(g_src).reshape((-1, 1))
        return g, self.lap - jac_src

    def solve(self, phi, rate, y_init=None, tol=1e-12, maxfev=60, jacobian="exact",
              check=True, lin_solver_kwargs=None):
        """Steady profile at Thiele modulus ``phi``.

        ``newton`` with ``bounds=(0, 1)`` (the maximum principle keeps ``y`` in
        ``[0, 1]`` for every rate law here) and ``stagnation=True``.

        Parameters
        ----------
        phi : float
            Thiele modulus.
        rate : Power, WeiszHicks or InhibitedLH
            Rate law.
        y_init : numpy.ndarray, optional
            Initial guess; default ``y = 1``.
        tol : float, optional
            Step tolerance of ``newton``.
        maxfev : int, optional
            Iteration limit.
        jacobian : {"exact", "numjac"}, optional
            See :meth:`residual`.
        check : bool, optional
            Raise :class:`SolveFailed` when ``newton`` reports no success.
        lin_solver_kwargs : dict, optional
            Passed to ``newton`` (``splu`` options such as ``permc_spec``).

        Returns
        -------
        y : numpy.ndarray
            Cell values, shape ``(n_u,)``.
        sol : OptimizeResult
            The ``newton`` result.
        """
        y0 = np.ones(self.shape) if y_init is None else np.asarray(y_init, float)
        sol = newton(lambda y: self.residual(y, phi, rate, jacobian), y0.reshape((-1, 1)),
                     tol=tol, maxfev=maxfev, bounds=(0.0, 1.0), stagnation=True,
                     lin_solver_kwargs=lin_solver_kwargs)
        if check and not sol.success:
            raise SolveFailed(f"phi = {phi:g}: {sol.message}")
        return np.asarray(sol.x).ravel(), sol

    def Phi(self, y, read="face"):
        """Observable ``(nu + 1) dy/du`` at the surface.

        ``read="face"`` takes the surface gradient the operator transports,
        from ``compute_boundary_values``; ``"cell"`` the two-point gradient
        from the last cell centre to the surface value (a defect the page
        injects).

        Parameters
        ----------
        y : numpy.ndarray
            Concentration profile on the cells.
        read : {"face", "cell"}, optional
            How the surface gradient is read.

        Returns
        -------
        float
            The observable ``Phi``.
        """
        y = np.asarray(y, float).reshape(self.shape)
        y_s, g_s = compute_boundary_values(y, self.u_f, self.u_c, self.bc[1], bound_id=1)
        g_s = float(np.ravel(g_s)[0])
        if read == "cell":
            g_s = (float(np.ravel(y_s)[0]) - y[-1, 0]) / (self.u_f[-1] - self.u_c[-1])
        return (self.nu + 1.0) * g_s

    def eta(self, y, phi, read="face"):
        """Effectiveness factor ``Phi/phi^2``.

        Parameters
        ----------
        y : numpy.ndarray
            Concentration profile on the cells.
        phi : float
            Thiele modulus.
        read : {"face", "cell"}, optional
            See :meth:`Phi`.

        Returns
        -------
        float
            Effectiveness factor.
        """
        return self.Phi(y, read) / phi ** 2


# ------------------------------------------------------------------- the film
def film_mass_dev(omega, n, area=3.0):
    """Rate deviation at bulk conditions across a lumped mass film.

    The film balance with the observed rate gives ``1 - C_s/C_b = omega/area``
    (``area = a r_p``, 3 for a sphere), so ``dev = 1 - (1 - omega/area)^n``.
    ``area=2`` is the cylinder's factor (a defect the page injects).

    Parameters
    ----------
    omega : float
        Mears group of eq. 17, observed rate over film transport.
    n : float
        Reaction order.
    area : float, optional
        Film area factor ``a r_p``.

    Returns
    -------
    float
        Relative rate deviation, bulk against surface conditions.
    """
    return 1.0 - (1.0 - omega / area) ** n


def film_heat_dev(gamma_chi, gamma_b, area=3.0):
    """Rate deviation at bulk conditions across a lumped heat film.

    ``(T_s - T_b)/T_b = chi/area`` and the Arrhenius factor give
    ``dev = exp[gamma_b u/(1 + u)] - 1`` with ``u = gamma_chi/(gamma_b area)``;
    ``gamma_chi`` is signed, positive for an exothermic reaction.

    Parameters
    ----------
    gamma_chi : float
        Arrhenius number times the film temperature group of eq. 14.
    gamma_b : float
        Arrhenius number at bulk temperature.
    area : float, optional
        Film area factor ``a r_p``.

    Returns
    -------
    float
        Relative rate deviation, bulk against surface conditions.
    """
    u = gamma_chi / (gamma_b * area)
    return float(np.expm1(gamma_b * u / (1.0 + u)))


# ----------------------------------------------------------- axial dispersion
def frac_unconverted(ktau, Pe_L):
    """First order, closed vessel: exit ``c/c*``, the classical closed form.

    Parameters
    ----------
    ktau : float
        Damkoehler number ``k tau``.
    Pe_L : float
        Axial Peclet number.

    Returns
    -------
    float
        Exit concentration over feed concentration.
    """
    q = np.sqrt(1.0 + 4.0 * ktau / Pe_L)
    den = (1 + q) ** 2 * np.exp(q * Pe_L / 2) - (1 - q) ** 2 * np.exp(-q * Pe_L / 2)
    return 4 * q * np.exp(Pe_L / 2) / den


class AxialVessel:
    """Closed vessel with axial dispersion, ``0 <= z <= 1``, ``v = 1``, ``D = 1/Pe_L``.

    ``d(c)/dz - (1/Pe_L) d2c/dz2 + ktau c^n = 0`` with Danckwerts conditions.
    Convection by van Leer TVD face values. The operators are assembled once.

    Parameters
    ----------
    Pe_L : float
        Vessel Peclet number ``v L/D_a``.
    n_z : int, optional
        Number of cells.
    inlet : {"danckwerts", "dirichlet"}, optional
        ``"dirichlet"`` fixes ``c = c*`` at the inlet (a defect the page
        injects).
    limiter : callable or None, optional
        TVD limiter; None is first-order upwind convection (a defect the page
        injects).
    """

    def __init__(self, Pe_L, n_z=800, inlet="danckwerts", limiter=vanleer):
        self.Pe_L, self.n_z, self.shape = float(Pe_L), n_z, (n_z, 1)
        self.z_f = np.linspace(0.0, 1.0, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        d_ax = 1.0 / self.Pe_L
        # bc: a dc/dn + b c = d with n the OUTWARD normal (-z at the inlet, +z at the outlet)
        if inlet == "danckwerts":
            # inlet z = 0: v c* = v c - D dc/dz, and dc/dn = -dc/dz, so
            #   (D/v) dc/dn + c = c*                                -> a = D, b = 1, d = 1
            bc_in = {"a": d_ax, "b": 1.0, "d": 1.0}
        else:
            # inlet z = 0: c = c*                                   -> a = 0, b = 1, d = 1
            bc_in = {"a": 0.0, "b": 1.0, "d": 1.0}
        # outlet z = 1: dc/dz = 0, dc/dn = dc/dz                    -> a = 1, b = 0, d = 0
        self.bc = (bc_in, {"a": 1.0, "b": 0.0, "d": 0.0})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c, self.bc, v=1.0)
        grad, grad_bc = construct_grad(self.shape, self.z_f, self.z_c, self.bc)
        self.div = construct_div(self.shape, self.z_f, nu=0)   # nu = 0: Cartesian
        # diffusion and its boundary part, and upwind convection plus diffusion
        self.jac_diff = (-d_ax * (self.div @ grad)).tocsc()
        self.g_diff = -d_ax * (self.div @ grad_bc).toarray().ravel()
        self.jac_up = (self.div @ (conv - d_ax * grad)).tocsc()
        self.g_up = (self.div @ (conv_bc - d_ax * grad_bc)).toarray().ravel()
        self.conv, self.conv_bc, self.limiter = conv, conv_bc, limiter
        if limiter is not None:
            self.tvd = ConvFluxTVD(self.shape, self.z_f, self.z_c, self.bc, tvd_limiter=limiter)
        self.eye = eye_array(n_z, format="csc")

    def _upwind_guess(self, ktau):
        """First-order upwind solution of the first-order problem: one linear solve."""
        return splu((self.jac_up + ktau * self.eye).tocsc()).solve(-self.g_up)

    def solve(self, ktau, n_react=1, tol=1e-12, maxfev=50):
        """Steady profile by ``newton`` with the exact TVD Jacobian.

        ``ConvFluxTVD`` gives the van Leer face fluxes and their Jacobian;
        ``newton(damping="monotone", stagnation=True)`` from the upwind
        solution of the first-order problem.

        Parameters
        ----------
        ktau : float
            Damkoehler number ``k tau``.
        n_react : float, optional
            Reaction order.
        tol : float, optional
            Newton step tolerance (concentrations are of order one).
        maxfev : int, optional
            Iteration limit.

        Returns
        -------
        c : numpy.ndarray
            Cell values, shape ``(n_z,)``.
        sol : OptimizeResult
            The ``newton`` result.

        Raises
        ------
        SolveFailed
            When ``newton`` reports no success.
        """
        def residual(c):
            c = np.asarray(c).reshape(-1)
            if self.limiter is None:
                flux, jac_conv = self.conv @ c + self.conv_bc.toarray().ravel(), self.conv
            else:
                flux, jac_conv = self.tvd(c.reshape(self.shape), 1.0)
            cp = np.clip(c, 0.0, None)
            src = ktau * cp ** n_react
            dsrc = ktau * n_react * cp ** (n_react - 1.0)
            g = self.div @ np.asarray(flux).reshape(-1) + self.jac_diff @ c + self.g_diff + src
            jac = self.div @ jac_conv + self.jac_diff + construct_coefficient_matrix(
                dsrc.reshape(self.shape))
            return g, jac.tocsc()

        sol = newton(residual, self._upwind_guess(ktau).reshape(self.shape), tol=tol,
                     maxfev=maxfev, damping="monotone", stagnation=True)
        if not sol.success:
            raise SolveFailed(f"axial vessel, Pe_L = {self.Pe_L:g}, ktau = {ktau:g}: "
                              f"{sol.message}")
        return np.asarray(sol.x).ravel(), sol

    def solve_deferred(self, ktau, n_react=1, tol=1e-13, max_it=300):
        """The same discrete equations by a fixed-point iteration (a validation route).

        Each sweep solves the upwind operator plus the source linearised as
        ``ktau c_old^(n-1) c`` with the van Leer correction of the previous
        iterate on the right-hand side (deferred correction).

        Returns
        -------
        c : numpy.ndarray
            Cell values.
        sweeps : int
            Number of sweeps.

        Raises
        ------
        SolveFailed
            When ``max_it`` sweeps do not bring the largest change below ``tol``.
        """
        c = self._upwind_guess(ktau)
        if self.limiter is None and n_react == 1:
            return c, 0
        for sweep in range(1, max_it + 1):
            w = ktau * np.clip(c, 0.0, None) ** (n_react - 1.0)
            lu = splu((self.jac_up + construct_coefficient_matrix(w.reshape(self.shape))).tocsc())
            if self.limiter is None:
                corr = 0.0
            else:
                _, dc_f = interp_cntr_to_stagg_tvd(c.reshape(self.shape), self.z_f, self.z_c,
                                                   self.bc, 1.0, tvd_limiter=self.limiter,
                                                   axis=0)
                corr = (self.div @ dc_f.reshape(-1, 1)).ravel()
            c_new = lu.solve(-self.g_up - corr)
            done = np.max(np.abs(c_new - c)) < tol
            c = c_new
            if done:
                return c, sweep
        raise SolveFailed(f"deferred correction: no convergence in {max_it} sweeps")

    def outlet(self, c, read="face"):
        """Outlet ``c/c*``: the face value from ``compute_boundary_values``, or
        ``read="cell"`` the last cell centre (a defect the page injects).

        Parameters
        ----------
        c : numpy.ndarray
            Concentration on the cells, scaled by the feed.
        read : {"face", "cell"}, optional
            How the outlet value is read.

        Returns
        -------
        float
            Outlet ``c/c*``.
        """
        c = np.asarray(c, float).reshape(self.shape)
        if read == "cell":
            return float(c[-1, 0])
        return float(np.ravel(compute_boundary_values(c, self.z_f, self.z_c, self.bc)[2])[0])


def axial_dev(X, n_react=1, n_z=None, route="newton", inlet="danckwerts", read="face",
              limiter=vanleer, tol=1e-12, sweep_tol=1e-13):
    """Deviation in required ``k tau`` at eq. 30/32's own boundary, conversion ``X``.

    At the boundary ``Pe_L = 20 n ln(1/(1 - X))``. The ``k tau`` the dispersed
    vessel needs for conversion ``X`` is root-found by Brent's method between
    the plug-flow value and twice it, and compared with plug flow.

    Parameters
    ----------
    X : float
        Conversion.
    n_react : float, optional
        Reaction order.
    n_z : int or None, optional
        None uses the closed form (first order only); otherwise the pymrm
        vessel on ``n_z`` cells.
    route : {"newton", "deferred"}, optional
        Solver of the pymrm vessel.
    inlet : {"danckwerts", "dirichlet"}, optional
        See :class:`AxialVessel`.
    read : {"face", "cell"}, optional
        See :meth:`AxialVessel.outlet`.
    limiter : callable or None, optional
        See :class:`AxialVessel`.
    tol, sweep_tol : float, optional
        Newton step tolerance; largest change per sweep of the deferred route.

    Returns
    -------
    float
        ``ktau/ktau_plug - 1``.
    """
    lnr = np.log(1.0 / (1.0 - X))
    Pe_L = 20.0 * n_react * lnr
    ktau_plug = lnr if n_react == 1 else X / (1.0 - X)        # plug flow, c* = 1
    if n_z is None:
        if n_react != 1:
            raise ValueError("the closed form is first order only")

        def f(k):
            return frac_unconverted(k, Pe_L) - (1.0 - X)
    else:
        vessel = AxialVessel(Pe_L, n_z, inlet=inlet, limiter=limiter)

        def f(k):
            c = (vessel.solve(k, n_react, tol=tol)[0] if route == "newton"
                 else vessel.solve_deferred(k, n_react, tol=sweep_tol)[0])
            return vessel.outlet(c, read) - (1.0 - X)
    ktau = brentq(f, ktau_plug, 2.5 * ktau_plug, xtol=1e-14, rtol=8.9e-16)
    return ktau / ktau_plug - 1.0


def _exit(ktau, Pe_L, n_z, inlet):
    """First-order exit ``c/c*``: the closed form, or the pymrm vessel on ``n_z`` cells.

    Parameters
    ----------
    ktau, Pe_L : float
        Damkoehler and axial Peclet numbers.
    n_z : int or None
        Cells of the pymrm vessel; None uses the closed form.
    inlet : {"danckwerts", "dirichlet"}
        Inlet condition of the vessel.

    Returns
    -------
    float
        Exit concentration over feed concentration.
    """
    if n_z is None:
        return frac_unconverted(ktau, Pe_L)
    vessel = AxialVessel(Pe_L, n_z, inlet=inlet)
    return vessel.outlet(vessel.solve(ktau)[0])


def axial_dev_covary(X, n_z=None, inlet="danckwerts"):
    """Extra required length at fixed ``k``, ``v`` and ``D_a``, first order.

    The length ratio ``l`` scales ``k tau`` and ``Pe_L`` together; root-found
    on the closed form (``n_z=None``) or on the pymrm vessel. Returns ``l - 1``.

    Parameters
    ----------
    X : float
        Conversion.
    n_z : int or None, optional
        Cells of the pymrm vessel; None uses the closed form.
    inlet : {"danckwerts", "dirichlet"}, optional
        Inlet condition of the vessel.

    Returns
    -------
    float
        Relative extra length ``l - 1``.
    """
    lnr = np.log(1.0 / (1.0 - X))
    return brentq(lambda l: _exit(l * lnr, l * 20.0 * lnr, n_z, inlet) - (1.0 - X),
                  1.0, 2.0, xtol=1e-14, rtol=8.9e-16) - 1.0


def petersen_dev(X, n_z=None, inlet="danckwerts"):
    """Deviation in required ``k tau`` at Petersen's boundary ``alpha = 1``, first order.

    ``alpha^2 = k tau/Pe_L = 1``, so the vessel with ``k tau = Pe_L`` that
    reaches conversion ``X`` is root-found and compared with plug flow.

    Parameters
    ----------
    X : float
        Conversion.
    n_z : int or None, optional
        Cells of the pymrm vessel; None uses the closed form.
    inlet : {"danckwerts", "dirichlet"}, optional
        Inlet condition of the vessel.

    Returns
    -------
    float
        Relative deviation of the required ``k tau`` from plug flow.
    """
    lnr = np.log(1.0 / (1.0 - X))
    pe = brentq(lambda p: _exit(p, p, n_z, inlet) - (1.0 - X), lnr, 3.0 * lnr, xtol=1e-14,
                rtol=8.9e-16)
    return pe / lnr - 1.0


# ------------------------------------------------------- stacked tolerances
def joint_worst(Bim, n_scan=110, detail=False, refine=True, c21=0.33, c17=0.15):
    """Worst first-order deviation among film-pellet states passing eqs. 1, 17, 21.

    States on the branch at mass Biot number ``Bim`` (sphere, :func:`shoot`
    with a film) pass when ``omega < c17`` (eq. 17), ``Phi_b < 1/(1 + c21
    omega)`` (eq. 21, isothermal part) and ``Phi_b/y_s < 1`` (eq. 1 on the
    surface observable). The supremum sits on the edge of the passing set, so
    every edge crossing on the scan ladder is root-found.

    Parameters
    ----------
    Bim : float
        Mass Biot number.
    n_scan : int, optional
        Scan points in ``t`` on ``[1e-4, 10^1.25]``.
    detail : bool, optional
        Also return the state at the worst point.
    refine : bool, optional
        ``False`` keeps only the passing scan points (a defect the page
        injects: a grid maximum under-reads the supremum).
    c21, c17 : float, optional
        The constants of eqs. 21 and 17.

    Returns
    -------
    float or tuple
        ``1 - eta_b`` at the worst passing state; with ``detail``, also
        ``(phi, eta_b, Phi_b, omega, y_s)`` there.
    """
    rate = Power(1.0)

    def state(t):
        return shoot(rate, t, bim=Bim)

    def margin(o):                                     # > 0: some criterion fails
        phi, eta_b, Phi_b, om, ys = o
        return max(om - c17, Phi_b - 1.0 / (1.0 + c21 * om), Phi_b / ys - 1.0)

    best, arg, prev = 0.0, None, None
    for t in np.geomspace(1e-4, 10 ** 1.25, n_scan):
        o = state(t)
        if o is None:
            prev = None
            continue
        g = margin(o)
        if g < 0 and 1 - o[1] > best:
            best, arg = 1 - o[1], o
        if refine and prev is not None and prev[1] * g < 0:
            oo = state(brentq(lambda x: margin(state(x)), prev[0], t, xtol=1e-13))
            if 1 - oo[1] > best:
                best, arg = 1 - oo[1], oo
        prev = (t, g)
    return (best, arg) if detail else best
