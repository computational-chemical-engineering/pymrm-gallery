"""Heck, Wei and Katzer (1976): duct Nusselt numbers, the Graetz problem and a reacting channel.

Three groups of solvers, all built from pymrm operators assembled once:

* :func:`duct`: fully developed flow and heat transfer in a duct cross-section
  (Table 1 of the paper): the Poiseuille velocity, the constant-flux (H1)
  temperature and the u-weighted eigenproblem of the constant-wall-temperature
  case, on a Cartesian quarter square or a cylindrical radius.
* :class:`Graetz`: the developing temperature field in a round tube (eqs. 15-16
  of the paper) on the whole (X, R) domain, van Leer fluxes in X, with an
  imposed wall flux or wall temperature.
* :class:`Channel`: the reacting monolith channel, the one-dimensional model
  (eqs. 8-14) and the two-dimensional one (eqs. 15-18 with Le = 1), on the
  Voltz kinetics of :class:`Voltz`. Light-off positions are roots inside the
  axial step that crosses the fold, not the step itself.

No plotting and no data loading: the notebook reads the printed values and
passes them in.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.sparse import diags_array
from scipy.sparse.linalg import eigsh, splu

from pymrm import (ConvFluxTVD, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, construct_grad,
                   interp_cntr_to_stagg_tvd, newton, vanleer)


# ----------------------------------------------------------------- helpers

def richardson(vals, ratio=2.0, order=2.0):
    """Richardson-extrapolate the last two rungs of a refinement ladder.

    Parameters
    ----------
    vals : sequence of float
        Values on successively refined grids.
    ratio : float, optional
        Refinement ratio between the last two rungs.
    order : float, optional
        Assumed order of the leading error term.

    Returns
    -------
    float
        The extrapolated value.
    """
    return vals[-1] + (vals[-1] - vals[-2]) / (ratio ** order - 1.0)


def observed_order(errs, ratio=2.0):
    """Rung-to-rung observed orders ``log(e_i / e_{i+1}) / log(ratio)``.

    Parameters
    ----------
    errs : sequence of float
        Errors on successively refined grids.
    ratio : float, optional
        Refinement ratio.

    Returns
    -------
    list of float
    """
    return [float(np.log(errs[i] / errs[i + 1]) / np.log(ratio))
            for i in range(len(errs) - 1)]


def fitted_order(ns, errs):
    """Least-squares slope of ``log(error)`` against ``log(n)``.

    Parameters
    ----------
    ns : sequence of int
        Cell counts.
    errs : sequence of float
        Errors at those counts.

    Returns
    -------
    float
        Minus the slope, so a second-order ladder gives 2.
    """
    return float(-np.polyfit(np.log(np.asarray(ns, float)),
                             np.log(np.asarray(errs, float)), 1)[0])


def cross(x, y):
    """Root of ``y(x)`` at its last sign change, on the piecewise-linear interpolant.

    Parameters
    ----------
    x, y : numpy.ndarray
        Samples.

    Returns
    -------
    float
        The root, or NaN when ``y`` does not change sign.
    """
    i = np.where(np.diff(np.sign(y)) != 0)[0]
    if len(i) == 0:
        return np.nan
    j = i[-1]
    return brentq(lambda t: np.interp(t, x, y), x[j], x[j + 1])


_FLUX = re.compile(r"^\s*([+-]?[\d.]+)\s*\+\s*\(\s*([\d.]+)\s*X\s*\)\s*\*\*\s*([\d.]+)\s*$")


def parse_flux(expr):
    """Coefficients ``(a, b, p)`` of a printed wall flux ``a + (b X)**p``.

    Parameters
    ----------
    expr : str
        The flux as transcribed, for example ``"-3+(40 X)**8"``.

    Returns
    -------
    tuple of float
    """
    m = _FLUX.match(expr)
    if m is None:
        raise ValueError(f"not a wall flux of the form a+(b X)**p: {expr!r}")
    return tuple(float(g) for g in m.groups())


def wall_flux(a, b, p):
    """The wall flux ``Q_w(X) = a + (b X)**p`` as a function of ``X``."""
    return lambda X: a + (b * np.asarray(X, float)) ** p


def correlation(expr):
    """A printed correlation in ``X`` (Python syntax, ``exp`` allowed) as a function.

    Parameters
    ----------
    expr : str
        For example ``"4.364+8.68*X**-0.506*exp(-0.041*X)"``.

    Returns
    -------
    callable
        ``f(X)``, evaluated with numpy.
    """
    code = compile(expr, "<correlation>", "eval")
    if set(code.co_names) - {"X", "exp"}:
        raise ValueError(f"unexpected names in {expr!r}")
    return lambda X: eval(code, {"__builtins__": {}, "exp": np.exp}, {"X": X})


# ----------------------------------------------------------- Table 1 ducts

def closed_form_nusselt():
    """Exact constant-flux Nusselt numbers of the circular duct and the equilateral triangle.

    Solves ``lap u = -1``, ``u = 0`` at the wall, and ``lap th = u``, ``th = 0``
    at the wall, by symbolic polynomials (sympy), then
    ``Nu = -4 ub A^2 / (P^2 thb)`` with ``thb`` the mixing-cup mean. For the
    triangle the velocity is the product of the three wall distances and the
    temperature that product times a quadratic.

    Returns
    -------
    tuple of sympy.Rational
        ``(Nu_circle, Nu_triangle)``, 48/11 and 28/9.
    """
    import sympy as sp

    r_, R0_ = sp.symbols("r R0", positive=True)
    u_c = (R0_ ** 2 - r_ ** 2) / 4
    th_c = sp.integrate(sp.integrate(u_c * r_, r_) / r_, r_)
    th_c = sp.simplify(th_c - th_c.subs(r_, R0_))
    I_c = lambda f: sp.integrate(f * 2 * sp.pi * r_, (r_, 0, R0_))   # noqa: E731
    A_c, P_c = sp.pi * R0_ ** 2, 2 * sp.pi * R0_
    nu_c = sp.nsimplify(sp.simplify(
        -4 * (I_c(u_c) / A_c) * A_c ** 2 / (P_c ** 2 * (I_c(th_c * u_c) / I_c(u_c)))))

    x_, y_ = sp.symbols("x y", real=True)
    s3 = sp.sqrt(3)
    l1, l2, l3 = y_ + 1, (-y_ + s3 * x_) / 2 + 1, (-y_ - s3 * x_) / 2 + 1   # incircle radius 1
    lap = lambda f: sp.diff(f, x_, 2) + sp.diff(f, y_, 2)               # noqa: E731
    cc = sp.symbols("cc")
    cc_val = sp.solve(sp.Poly(sp.expand(lap(cc * l1 * l2 * l3) + 1), x_, y_).coeffs(), cc)[cc]
    u_t = sp.expand(cc_val * l1 * l2 * l3)
    aa, bb = sp.symbols("aa bb")
    th_t = sp.expand(l1 * l2 * l3 * (aa + bb * (x_ ** 2 + y_ ** 2)))
    sol_t = sp.solve(sp.Poly(sp.expand(lap(th_t) - u_t), x_, y_).coeffs(), [aa, bb])
    th_t = sp.expand(th_t.subs(sol_t))
    assert sp.simplify(lap(th_t) - u_t) == 0 and sp.simplify(th_t.subs(y_, -1)) == 0
    xl, xr = sp.solve(sp.Eq(l3, 0), x_)[0], sp.solve(sp.Eq(l2, 0), x_)[0]
    I_t = lambda f: sp.integrate(sp.integrate(f, (x_, xr, xl)), (y_, -1, 2))   # noqa: E731
    A_t, P_t = I_t(1), 6 * sp.sqrt(3)
    nu_t = sp.nsimplify(sp.simplify(
        -4 * (I_t(u_t) / A_t) * A_t ** 2 / (P_t ** 2 * (I_t(th_t * u_t) / I_t(u_t)))))
    return nu_c, nu_t


def duct(nx, ny=None, geom="square", nu_r=1, plug=False, wall_neumann=False,
         weight_eigen=True):
    """Fully developed duct flow: Nu for constant wall flux (H1) and wall temperature (T).

    Parameters
    ----------
    nx : int
        Cells along x (square) or r (circle).
    ny : int, optional
        Cells along y for the square; ``nx`` when omitted.
    geom : {"square", "circle"}
        Quarter of a square of half-width 1 with symmetry planes at x = y = 0,
        or a circle of radius 1.
    nu_r : int, optional
        Geometry factor of the radial divergence for the circle (1 cylindrical;
        0 is a break-table defect).
    plug, wall_neumann, weight_eigen : bool, optional
        Break-table defects: plug velocity, no wall condition, plain Laplacian
        eigenvalue.

    Returns
    -------
    dict
        ``Nu_H1``, ``Nu_T``, ``fRe`` and the mean velocity ``ub``.
    """
    if geom == "circle":
        f = np.linspace(0.0, 1.0, nx + 1)
        c = 0.5 * (f[:-1] + f[1:])
        shape, grids = (nx, 1), [(f, c, nu_r)]
        dv = (f[1:] ** 2 - f[:-1] ** 2) / 2.0                 # int r dr, 2*pi dropped
        A_full, P_full = np.pi, 2 * np.pi
    else:
        ny = ny or nx
        fx = np.linspace(0.0, 1.0, nx + 1)
        cx = 0.5 * (fx[:-1] + fx[1:])
        fy = np.linspace(0.0, 1.0, ny + 1)
        cy = 0.5 * (fy[:-1] + fy[1:])
        shape, grids = (nx, ny), [(fx, cx, 0), (fy, cy, 0)]
        dv = np.outer(np.diff(fx), np.diff(fy)).ravel()
        A_full, P_full = 4.0, 8.0                             # full square, half-width 1
    # a dphi/dn + b phi = d, n the OUTWARD normal.
    #   symmetry plane (lower): dphi/dn = 0 -> a=1, b=0, d=0
    #   wall (upper):          phi = 0      -> a=0, b=1, d=0 (Dirichlet)
    wall = {"a": 1.0, "b": 0.0, "d": 0.0} if wall_neumann else {"a": 0.0, "b": 1.0, "d": 0.0}
    bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, wall)
    lap, lap_bc = 0, 0
    for ax, (f, c, nu) in enumerate(grids):
        g, gbc = construct_grad(shape, f, c, bc, axis=ax)
        d = construct_div(shape, f, nu=nu, axis=ax)            # nu=0 Cartesian, 1 cylindrical
        lap = lap + d @ g
        lap_bc = lap_bc + d @ gbc
    lap_bc = np.asarray(lap_bc.todense()).ravel()
    lu = splu(lap.tocsc())
    u = np.ones(int(np.prod(shape))) if plug else lu.solve(-1.0 - lap_bc)
    th = lu.solve(u - lap_bc)
    ub = (u * dv).sum() / dv.sum()
    thb = (th * u * dv).sum() / (u * dv).sum()
    nu_h1 = -4.0 * ub * A_full ** 2 / (P_full ** 2 * thb)
    w = dv * u if weight_eigen else dv
    # v0 fixed: ARPACK's default start vector is random
    lam = eigsh((-diags_array(dv) @ lap).tocsc(), k=1, M=diags_array(w).tocsc(),
                sigma=0.0, which="LM", v0=np.ones(int(np.prod(shape))))[0][0]
    nu_t = 4.0 * lam * ub * A_full ** 2 / P_full ** 2
    fre = (4.0 * A_full / P_full) ** 2 / (2.0 * ub)
    return dict(Nu_H1=nu_h1, Nu_T=nu_t, fRe=fre, ub=ub)


def square_spectral(N, nq=160, weight_eigen=True):
    """Square duct by a sine-series Galerkin method on (0, 2)^2, sharing no code with :func:`duct`.

    Parameters
    ----------
    N : int
        Sine modes per direction.
    nq : int, optional
        Gauss-Legendre points per direction for the weighted mass matrix.
    weight_eigen : bool, optional
        False drops the velocity weight (a break-table defect).

    Returns
    -------
    tuple of float
        ``(Nu_H1, Nu_T)``. Neither is a bound: the eigenvalue converges from
        above (Rayleigh-Ritz) and the truncated mean velocity from below.
    """
    import scipy.linalg as sl

    m = np.arange(1, N + 1)
    k = m * np.pi / 2.0
    b1 = 2.0 * (1.0 - np.cos(m * np.pi)) / (m * np.pi)
    K2 = k[:, None] ** 2 + k[None, :] ** 2
    U = np.outer(b1, b1) / K2                      # -lap u = 1
    ints = (1.0 - np.cos(m * np.pi)) / k
    IU = (U * np.outer(ints, ints)).sum()
    ub = IU / 4.0
    TH = -U / K2                                   # lap th = u
    thb = (TH * U).sum() / IU
    A, Pm = 4.0, 8.0
    nu_h1 = -4 * ub * A ** 2 / (Pm ** 2 * thb)
    xg, wg = np.polynomial.legendre.leggauss(nq)
    xg = xg + 1.0
    S = np.sin(np.outer(k, xg))
    ug = np.einsum("ab,ai,bj->ij", U, S, S, optimize=True)
    Wq = np.outer(wg, wg) * (ug if weight_eigen else 1.0)
    Amp = np.einsum("mi,pi->mpi", S, S).reshape(N * N, nq)
    # (Amp Wq Amp^T) is indexed ((m,p),(n,q)); the eigenproblem needs ((m,n),(p,q))
    Mm = np.transpose((Amp @ Wq @ Amp.T).reshape(N, N, N, N), (0, 2, 1, 3)).reshape(N * N, N * N)
    lam0 = sl.eigh(np.diag(K2.ravel()), Mm, eigvals_only=True, subset_by_index=[0, 0])[0]
    return nu_h1, 4 * lam0 * ub * A ** 2 / Pm ** 2


# ------------------------------------------------------------ Graetz tube

class Graetz:
    """Developing temperature field in a round tube, eqs. (15)-(16) of the paper.

    ``half * (1 - R^2) dT/dX = (1/R) d/dR (R dT/dR)`` on ``0 < X < X_end``,
    ``0 < R < 1``, state shape ``(n_x, n_r)``: axis 0 carries convection, axis 1
    diffusion. The operators are assembled once per wall condition.

    Parameters
    ----------
    X : float
        Length of the domain in X.
    n_x, n_r : int
        Cells along X and R.
    half : float, optional
        Coefficient of the axial term, 1/2 as printed in eq. (15); 2 is the
        coefficient the Notation's ``X* = 4x/(D Re Pr)`` would require.
    x_grade, r_grade : float, optional
        Grid grading: ``x_f = X t**x_grade``, ``r_f = 1 - (1 - t)**r_grade``.
    """

    def __init__(self, X, n_x, n_r, half=0.5, x_grade=1.0, r_grade=2.0):
        t = np.linspace(0.0, 1.0, n_x + 1)
        self.x_f = X * t ** x_grade
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        t = np.linspace(0.0, 1.0, n_r + 1)
        self.r_f = 1.0 - (1.0 - t) ** r_grade
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        self.shape = (n_x, n_r)
        self.v = np.broadcast_to(half * (1.0 - self.r_c ** 2), (n_x + 1, n_r))
        aw = 2.0 * self.r_c * np.diff(self.r_f)
        self.area_w = aw / aw.sum()
        self.cup_w = aw * (1.0 - self.r_c ** 2) / np.sum(aw * (1.0 - self.r_c ** 2))

    def _bcs(self, wall, inlet, outflow):
        n_x = self.shape[0]
        # a dT/dn + b T = d, n the OUTWARD normal.
        #   inlet X=0: T = inlet -> a=0, b=1, d=inlet (Dirichlet)
        #   outlet X=X_end: convection only in X, so no condition: {"outflow": True}
        #   (outflow=False imposes dT/dX = 0, a=1, b=0, d=0, a condition the
        #   equation does not have)
        bc_x = ({"a": 0.0, "b": 1.0, "d": inlet},
                {"outflow": True} if outflow else {"a": 1.0, "b": 0.0, "d": 0.0})
        if wall[0] == "flux":
            q = np.asarray(wall[1](self.x_c), float).reshape(n_x, 1)
            #   axis R=0: dT/dR = 0      -> a=1, b=0, d=0
            #   wall R=1: dT/dR = Qw(X)  -> a=1, b=0, d=Qw (n = +R there)
            bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": q})
        else:
            q = None
            #   wall R=1: T = value -> a=0, b=1, d=value (Dirichlet)
            bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0},
                    {"a": 0.0, "b": 1.0, "d": float(wall[1])})
        return bc_x, bc_r, q

    def _diffusion(self, bc_r):
        gr, gr_bc = construct_grad(self.shape, self.r_f, self.r_c, bc_r, axis=1)
        dr = construct_div(self.shape, self.r_f, nu=1, axis=1)      # nu=1: cylindrical, radial
        return (dr @ gr).tocsc(), np.asarray((dr @ gr_bc).todense()).ravel()

    def _result(self, T, bc_r, q, nit, incr, wall_read="face"):
        # the wall value from pymrm's boundary reconstruction on the face;
        # wall_read="cell" takes the last cell centre instead (a break-table defect)
        if wall_read == "cell":
            Tw = T[:, -1]
        else:
            Tw, _ = compute_boundary_values(T, self.r_f, self.r_c, bc_r[1], axis=1, bound_id=1)
        return dict(x=self.x_c, r=self.r_c, T=T, Tw=np.asarray(Tw).ravel(),
                    cup=T @ self.cup_w, area=T @ self.area_w,
                    q=(None if q is None else q.ravel()), nit=nit, incr=incr)

    def solve(self, wall, inlet=0.0, tvd=True, outflow=True, tol=1e-10, maxit=600,
              omega=0.6, accept=1e-8, wall_read="face"):
        """Solve by van Leer deferred correction on one LU of the upwind operator.

        Parameters
        ----------
        wall : tuple
            ``("flux", Qw)`` with ``Qw(X)`` the wall gradient, or
            ``("temp", value)``.
        inlet : float, optional
            Inlet temperature.
        tvd : bool, optional
            False keeps the first-order upwind solution.
        outflow : bool, optional
            Outflow outlet; False imposes a zero gradient.
        wall_read : {"face", "cell"}, optional
            Wall value from ``compute_boundary_values`` (the face) or from the
            last cell centre (a break-table defect).
        tol, maxit, omega : optional
            Stop when the largest change per unit relaxation is below ``tol``;
            relaxation factor ``omega``.
        accept : float, optional
            Raise if the last change after ``maxit`` sweeps exceeds this.

        Returns
        -------
        dict
            ``x``, ``r``, ``T``, wall value ``Tw``, mixing-cup ``cup`` and area
            ``area`` means per axial cell, the wall gradient ``q`` (flux walls),
            sweeps ``nit`` and last change ``incr``.
        """
        bc_x, bc_r, q = self._bcs(wall, inlet, outflow)
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                  bc_x, v=self.v, axis=0)
        dx = construct_div(self.shape, self.x_f, nu=0, axis=0)      # nu=0: Cartesian, axial
        lap, lap_bc = self._diffusion(bc_r)
        A = (dx @ conv - lap).tocsc()
        b = np.asarray((dx @ conv_bc).todense()).ravel() - lap_bc
        lu = splu(A)
        T = lu.solve(-b).reshape(self.shape)
        nit, incr = 0, 0.0
        if tvd:
            for nit in range(1, maxit + 1):
                _, dT = interp_cntr_to_stagg_tvd(T, self.x_f, self.x_c, bc_x,
                                                 self.v, vanleer, axis=0)
                corr = np.asarray(dx @ (self.v * dT).reshape(-1, 1)).ravel()
                new = T + omega * (lu.solve(-b - corr).reshape(self.shape) - T)
                incr = np.max(np.abs(new - T)) / omega
                T = new
                if incr < tol:
                    break
            else:
                if incr > accept:
                    raise RuntimeError(f"deferred correction did not converge: {incr:.2e}")
        return self._result(T, bc_r, q, nit, incr, wall_read)

    def solve_newton(self, wall, inlet=0.0, outflow=True, tol=1e-10, guess="upwind"):
        """Solve the same discrete equations by Newton with ``ConvFluxTVD``.

        The van Leer face fluxes and their Jacobian come from
        :class:`pymrm.ConvFluxTVD`; ``newton(damping="monotone")`` starts from
        the upwind solution (``guess="upwind"``) or from the inlet value
        everywhere (``guess="inlet"``). Each iteration factorises the full
        Jacobian.

        Returns
        -------
        dict
            As :meth:`solve`; ``nit`` is the number of Newton iterations.
        """
        bc_x, bc_r, q = self._bcs(wall, inlet, outflow)
        tvd = ConvFluxTVD(self.shape, self.x_f, self.x_c, bc_x, tvd_limiter=vanleer, axis=0)
        dx = construct_div(self.shape, self.x_f, nu=0, axis=0)
        lap, lap_bc = self._diffusion(bc_r)

        def residual(T):
            flux, jac = tvd(T, self.v)
            g = dx @ flux.reshape(-1) - lap @ T.reshape(-1) - lap_bc
            return g, (dx @ jac - lap).tocsc()

        T0 = (self.solve(wall, inlet, tvd=False, outflow=outflow)["T"] if guess == "upwind"
              else np.full(self.shape, float(inlet)))
        res = newton(residual, T0, damping="monotone", tol=tol)
        if not res.success:
            raise RuntimeError(f"newton did not converge: {res.message}")
        return self._result(res.x.reshape(self.shape), bc_r, q, res.nit, 0.0)


def local_nusselt(sol):
    """Eq. (19): ``Nu = 2 (dT/dR)_w / (T_w - Tbar)``, with the imposed wall gradient.

    Parameters
    ----------
    sol : dict
        A flux-wall result of :meth:`Graetz.solve`.

    Returns
    -------
    numpy.ndarray
        Local Nusselt number per axial cell, mixing-cup mean.
    """
    return 2.0 * sol["q"] / (sol["Tw"] - sol["cup"])


# -------------------------------------------------------- reacting channel

def reconstruct_groups(delta_H, CG0, D_ch, Sh, T_props, M_G=29.0, C_PG=1093.0,
                       P=101325.0, R_g=8314.0, D_AB300=0.208e-4, n_T=1.75, molar=True):
    """The two dimensional groups the paper does not print.

    Parameters
    ----------
    delta_H : float
        Heat of reaction, J/kmol.
    CG0 : float
        Inlet CO mole fraction.
    D_ch : float
        Channel diameter, m.
    Sh : float
        Sherwood number (= Nu).
    T_props : float
        Temperature for the gas properties, K.
    M_G, C_PG : float, optional
        Molar mass (kg/kmol) and heat capacity (J/kg K) of the exhaust gas.
    P, R_g : float, optional
        Pressure (Pa) and gas constant (J/kmol K).
    D_AB300 : float, optional
        CO diffusivity in air at 300 K (m^2/s).
    n_T : float, optional
        Temperature exponent of the diffusivity, ``D_AB300 (T/300)**n_T``.
    molar : bool, optional
        False reads ``C_PG`` as a molar heat capacity (a units defect).

    Returns
    -------
    tuple of float
        ``(dTad, beta, b0)``: adiabatic rise (K), ``beta = rho_G k_m / M_G``
        (kmol/m^2 s per unit mole fraction) and ``b0 = beta / Sh``.
    """
    D_AB = D_AB300 * (T_props / 300.0) ** n_T
    dTad = delta_H * CG0 / ((M_G if molar else 1.0) * C_PG)
    b0 = P / (R_g * T_props) * D_AB / D_ch
    return dTad, b0 * Sh, b0


@dataclass(frozen=True)
class Voltz:
    """Voltz et al. (1973) CO oxidation kinetics as printed by Heck et al., eqs. (5)-(7).

    ``r = k_r C C_O2 / (1 + k_a C)^2`` with Arrhenius ``k_r`` and ``k_a`` in the
    wall temperature (degC + 273). Concentrations are mole fractions.

    Parameters
    ----------
    kr0, ka0 : float
        Pre-exponential factors.
    err, ear : float
        ``E_r/R`` and ``E_a/R``, K.
    cg0, co20 : float
        Inlet CO and O2 mole fractions.
    o2 : {"const", "depleted"}
        Oxygen at its inlet value, or depleted stoichiometrically (1/2 per CO).
    """

    kr0: float
    ka0: float
    err: float
    ear: float
    cg0: float
    co20: float
    o2: str = "const"

    def rate(self, Cw, Tw):
        """Rate at wall mole fraction ``Cw`` and wall temperature ``Tw`` (degC); complex-safe."""
        TK = Tw + 273.0
        kr = self.kr0 * np.exp(-self.err / TK)
        ka = self.ka0 * np.exp(-self.ear / TK)
        if self.o2 == "const":
            CO2 = self.co20
        else:
            CO2 = self.co20 - 0.5 * (self.cg0 - Cw)
            CO2 = np.where(np.real(CO2) > 0.0, CO2, 0.0)
        return kr * Cw * CO2 / (1.0 + ka * Cw) ** 2


@dataclass
class Channel:
    """One monolith channel: the 1-D and the 2-D light-off models on the same groups.

    The wall temperature follows the adiabatic locus of eq. (14) as derived
    from eq. (13) with Le = 1 and Sh = Nu: ``T_w = T_G0 + dTad (1 - C_w/C_G0)``.
    Axial coordinate ``G = x/(D Re Pr)``.

    Parameters
    ----------
    kin : Voltz
        Kinetics.
    dTad : float
        Adiabatic temperature rise, K.
    b0 : float
        ``beta / Sh``, a property group; ``beta = b0 Nu``.
    """

    kin: Voltz
    dTad: float
    b0: float
    _folds: dict = field(init=False, default_factory=dict, repr=False, compare=False)

    @property
    def cg0(self):
        return self.kin.cg0

    def r_of(self, Cw, TG0):
        """Rate on the adiabatic locus."""
        return self.kin.rate(Cw, TG0 + self.dTad * (1.0 - Cw / self.cg0))

    def drdc(self, Cw, TG0):
        """``dr/dCw`` on the locus by complex step (no step-size choice)."""
        return np.imag(self.r_of(np.asarray(Cw, dtype=complex) + 1e-30j, TG0)) / 1e-30

    def fold(self, TG0, beta, c_max=None):
        """Fold of the cold branch for a mass-transfer coefficient ``beta``.

        The total sink ``F(c) = r(c) + beta c`` meets ``beta C_G``; the cold root
        merges with the middle one at the largest local minimum of F, the
        largest root of ``r'(c) = -beta`` below ``c_max``. NaN when there is
        none: then the branch is single valued (the paper's p. 480 condition).

        Returns
        -------
        tuple of float
            ``(c_fold, CG_fold)``: wall and bulk mole fractions at the fold.
        """
        key = (float(TG0), float(beta))
        if key not in self._folds:
            cs = np.geomspace(1e-8, self.cg0, 600)
            h = self.drdc(cs, TG0) + beta
            roots = [brentq(lambda c: self.drdc(c, TG0) + beta, cs[i], cs[i + 1],
                            xtol=1e-18, rtol=8.9e-16)
                     for i in np.where(np.diff(np.sign(h)) != 0)[0]]
            self._folds[key] = [(c2, (self.r_of(c2, TG0) + beta * c2) / beta) for c2 in roots]
        below = [rc for rc in self._folds[key] if c_max is None or rc[0] < c_max]
        return below[-1] if below else (np.nan, np.nan)

    def fold_CG(self, TG0, beta):
        """Bulk mole fraction at which the cold branch folds (light-off); NaN if none."""
        return self.fold(TG0, beta)[1]

    def cw_cold(self, CG, TG0, beta):
        """Wall mole fraction on the cold branch at bulk ``CG``."""
        return brentq(lambda c: self.r_of(c, TG0) - beta * (CG - c), 1e-14, CG,
                      xtol=1e-18, rtol=8.9e-16)

    def lightoff_quad(self, TG0, Sh, npts=120):
        """Light-off position of the constant-Nu 1-D model by quadrature in ``C_w``.

        ``G* = int_{c_fold}^{c_in} [beta + r'(c)] / [4 Sh r(c)] dc`` on the cold
        branch, Gauss-Legendre with ``npts`` points. The integrand vanishes at
        the fold, so the form has no endpoint singularity.

        Returns
        -------
        float
            ``G*``; ``inf`` without a fold, 0 when lit at the inlet.
        """
        beta = self.b0 * Sh
        cf, CGs = self.fold(TG0, beta)
        if not np.isfinite(cf):
            return np.inf
        if CGs >= self.cg0:
            return 0.0
        cin = self.cw_cold(self.cg0, TG0, beta)
        x, w = np.polynomial.legendre.leggauss(npts)
        c = 0.5 * (cin - cf) * x + 0.5 * (cin + cf)
        f = (beta + self.drdc(c, TG0)) / (4.0 * Sh * self.r_of(c, TG0))
        return float(0.5 * (cin - cf) * np.sum(w * f))

    def lightoff_quad_CG(self, TG0, Sh, limit=600):
        """The same integral parametrised by ``C_G``, adaptive ``quad``: a square-root endpoint."""
        beta = self.b0 * Sh
        CGs = self.fold_CG(TG0, beta)
        if not np.isfinite(CGs):
            return np.inf
        if CGs >= self.cg0:
            return 0.0
        f = lambda CG: 1.0 / (4.0 * Sh * (CG - self.cw_cold(CG, TG0, beta)))   # noqa: E731
        return quad(f, CGs * (1.0 + 1e-12), self.cg0, limit=limit)[0]

    def _cold_root(self, w0, TG0, beta):
        """Cold-branch wall value (scaled by C_G0) for scaled bulk ``w0``; NaN when folded."""
        cf, CGs = self.fold(TG0, beta, c_max=w0 * self.cg0)
        if np.isfinite(cf) and w0 * self.cg0 <= CGs:
            return np.nan
        lo = cf / self.cg0 if np.isfinite(cf) else 1e-12
        phi = lambda cw: w0 - self.r_of(cw * self.cg0, TG0) / (beta * self.cg0) - cw   # noqa: E731
        return brentq(phi, lo, w0, xtol=1e-16, rtol=8.9e-16)

    def _margin(self, w0, TG0, beta):
        """Scaled distance of the bulk from the fold; positive while the cold branch exists."""
        cf, CGs = self.fold(TG0, beta, c_max=w0 * self.cg0)
        return w0 - CGs / self.cg0 if np.isfinite(cf) else 1.0

    def march_1d(self, TG0, nu_fun, dG, Gmax, refine=True):
        """1-D channel ``dc/dG = -4 Nu(G) (c - c_w)``, explicit in G, ``beta(G) = b0 Nu(G)``.

        Parameters
        ----------
        TG0 : float
            Inlet temperature, degC.
        nu_fun : callable
            ``Nu(G)``.
        dG, Gmax : float
            Step and end of the channel.
        refine : bool, optional
            Root-find light-off inside the crossing step (the linear Euler
            update meets the fold); False returns the first point past the fold.

        Returns
        -------
        float
            Light-off position, ``inf`` if none before ``Gmax``.
        """
        c, G, prev = 1.0, 0.0, None
        while G < Gmax - 1e-15:
            nu = nu_fun(G)
            beta = self.b0 * nu
            cw = self._cold_root(c, TG0, beta)
            if not np.isfinite(cw):
                if prev is None or not refine:
                    return G
                G0, c0, s0 = prev
                m = lambda h: self._margin(c0 - h * s0, TG0, self.b0 * nu_fun(G0 + h))  # noqa: E731
                return G0 + brentq(m, 0.0, dG, xtol=1e-15, rtol=8.9e-16)
            s = 4.0 * nu * (c - cw)
            prev = (G, c, s)
            c = c - dG * s
            G += dG
        return np.inf

    def march_2d(self, TG0, dG, Gmax, n_r=48, grade=2.0, nu_r=1, refine=True):
        """2-D channel, eqs. (15)-(18) with Le = 1, backward Euler in G.

        With Le = 1 the scaled temperature equals ``1 - c`` everywhere, so one
        field remains. The interior is linear in the wall flux q, so the
        response to a unit flux is computed once per step length and each step
        is one linear solve (at q = 0) and a scalar root for the wall value:
        the 1-D algebra with an effective ``beta = gam / (S C_G0)`` and the
        wall value at zero flux as the bulk.

        Parameters
        ----------
        TG0 : float
            Inlet temperature, degC.
        dG, Gmax : float
            Step and end of the channel.
        n_r, grade : int, float, optional
            Radial cells and wall grading.
        nu_r : int, optional
            Geometry factor of the radial divergence (0 is a break-table defect).
        refine : bool, optional
            Root-find light-off on the length of the crossing step.

        Returns
        -------
        float
            Light-off position, ``inf`` if none before ``Gmax``.
        """
        rad = Radial(n_r, grade=grade, nu_r=nu_r)
        gam = 2.0 * self.b0 * self.cg0                 # eq. (18) group
        step = rad.stepper(dG)
        c, G = np.ones((n_r, 1)), 0.0
        while G < Gmax - 1e-15:
            c0, c1, w0, S = step(c)
            beta = gam / (S * self.cg0)
            cw = self._cold_root(w0, TG0, beta)
            if not np.isfinite(cw):
                if not refine:
                    return G
                def m(h, c=c):
                    _, _, w0h, Sh = rad.stepper(h)(c)
                    return self._margin(w0h, TG0, gam / (Sh * self.cg0))
                if m(1e-9 * dG) <= 0.0:              # folded within any step
                    return G
                return G + brentq(m, 1e-9 * dG, dG, xtol=1e-15, rtol=8.9e-16)
            q = -self.r_of(cw * self.cg0, TG0) / gam
            c = c0 + (c1 - c0) * q
            G += dG
        return np.inf

    def first_step_margin(self, TG0, dG, n_r=48, grade=2.0, nu_r=1):
        """Fold margin of the first 2-D step from the inlet; zero where light-off enters it."""
        rad = Radial(n_r, grade=grade, nu_r=nu_r)
        _, _, w0, S = rad.stepper(dG)(np.ones((n_r, 1)))
        return self._margin(w0, TG0, 2.0 * self.b0 / S)


class Radial:
    """Radial operators of the 2-D reacting channel, assembled once.

    Parameters
    ----------
    n_r : int
        Radial cells.
    grade : float, optional
        Wall grading, ``r_f = 1 - (1 - t)**grade``.
    nu_r : int, optional
        Geometry factor of the divergence, 1 cylindrical.
    """

    def __init__(self, n_r=48, grade=2.0, nu_r=1):
        t = np.linspace(0.0, 1.0, n_r + 1)
        self.r_f = 1.0 - (1.0 - t) ** grade
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        self.shape = (n_r, 1)
        # a dc/dn + b c = d, n the OUTWARD normal.
        #   R=0 axis: dc/dR = 0, symmetry -> a=1, b=0, d=0
        #   R=1 wall: dc/dR = q, imposed  -> a=1, b=0, d=1 times the external q
        #   (shapes_d returns the boundary vector per unit q; by linearity it equals
        #   the plain boundary vector of the dict with d = 1)
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 1.0})
        g, _, g_q = construct_grad(self.shape, self.r_f, self.r_c, self.bc, axis=0,
                                   shapes_d=(None, (1, 1)))
        d = construct_div(self.shape, self.r_f, nu=nu_r, axis=0)   # nu=1: cylindrical
        self.lap = (d @ g).tocsc()
        self.lap_q = np.asarray((d @ g_q).todense()).reshape(self.shape)   # per unit q
        self.w = 0.5 * (1.0 - self.r_c ** 2).reshape(self.shape)          # eq. (15) coefficient

    def wall_value(self, c, q):
        """Wall value of ``c`` at wall gradient ``q``, by ``compute_boundary_values``."""
        bc = {"a": 1.0, "b": 0.0, "d": q}
        return float(compute_boundary_values(c.reshape(self.shape), self.r_f, self.r_c, bc,
                                             axis=0, bound_id=1)[0].ravel()[0])

    def stepper(self, h):
        """One backward-Euler step of length ``h``: ``c -> (c0, c1, w0, S)``.

        ``c0`` and ``c1`` are the new fields at wall gradient 0 and 1, ``w0`` the
        wall value of ``c0`` and ``S`` the wall value gained per unit gradient.
        """
        accum = construct_coefficient_matrix(self.w / h, shape=self.shape)
        lu = splu((accum - self.lap).tocsc())
        dc = lu.solve(self.lap_q.ravel()).reshape(self.shape)   # response to unit q
        S = self.wall_value(dc, 1.0)       # the wall value is affine in (c, q)

        def step(c):
            c0 = lu.solve((self.w * c / h).ravel()).reshape(self.shape)
            return c0, c0 + dc, self.wall_value(c0, 0.0), S
        return step
