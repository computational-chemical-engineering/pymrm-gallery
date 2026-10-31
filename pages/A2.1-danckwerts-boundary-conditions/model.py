"""The Danckwerts closed vessel and Wehner & Wilhelm's three-section problem.

Model code of gallery page A2.1 (Danckwerts 1953; Wehner & Wilhelm 1956). The
page imports these; the closed forms it checks them against (Danckwerts' Eqs. 33
to 37, Wehner & Wilhelm's Eqs. 19 to 22) stay in the notebook, so that the pymrm
solve and the transcribed algebra never share a line of code.

Everything is dimensionless in Wehner & Wilhelm's groups: the bed is
``0 <= z <= 1``, the velocity is ``u = 1``, the axial diffusivity is ``1/Pe`` and
the first-order rate constant is ``R = kL/u``.

Contents
--------
ClosedVessel
    Danckwerts' closed vessel with his inlet and outlet conditions, one field
    on a uniform grid, van Leer TVD convection solved by Newton with the exact
    Jacobian (``ConvFluxTVD``, pitfalls P10), or by deferred correction as a
    second solver route.
ThreeSection
    Wehner & Wilhelm's fore | reaction | after problem as one pymrm domain with
    a piecewise diffusivity, harmonic face means at the two jumps (pitfalls P2)
    and no condition imposed at the bed ends.
Solution
    A converged closed-vessel solve.
geom_faces
    Faces of a geometrically stretched grid.

Notes
-----
Both models are linear apart from the limiter. Newton starts from the
first-order upwind solution ``c_up``, and the residual is written as the
correction to it, ``jac_up @ (c - c_up) + div @ (F_tvd(c) - F_up(c))``.
That is the full TVD residual minus the upwind residual of ``c_up``, which is
zero to round-off. Written this way the residual is a small quantity computed
from small ones: the full form, ``div @ F_tvd(c) + (diffusion + reaction) @ c +
g``, evaluates diffusion entries of order ``D/dz^2`` on order-one fractions,
and its round-off floor is set by them. Newton then stops on a state that
depends on the BLAS kernel. The page measures both forms against deferred
correction (``residual="full"`` reproduces the first). The Jacobian is the same
either way.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu

from pymrm import (ConvFluxTVD, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, construct_grad,
                   interp_cntr_to_stagg_tvd, newton, vanleer)

__all__ = ["ClosedVessel", "ThreeSection", "Solution", "geom_faces"]


@dataclass
class Solution:
    """A converged solve.

    Attributes
    ----------
    z_c : numpy.ndarray
        Cell centres, shape ``(n,)``.
    c : numpy.ndarray
        Fraction of reactant remaining at the cell centres, shape ``(n,)``.
    c_exit : float
        Value on the outlet face, from ``compute_boundary_values`` (pitfalls P8).
    nit : int
        Newton iterations or deferred-correction sweeps (0 for bare upwind).
    """

    z_c: np.ndarray
    c: np.ndarray
    c_exit: float
    nit: int


class _Linear1D:
    """Shared assembly and solvers of a steady, linear 1-D convection-dispersion model.

    Subclasses set ``shape``, ``z_f``, ``z_c``, ``bc``, ``u`` and call
    :meth:`_assemble` with the face-diffusivity matrix and the reaction matrix.
    """

    def _assemble(self, d_mat, k_mat, nu=0, permc_spec=None):
        shape, z_f, z_c, bc = self.shape, self.z_f, self.z_c, self.bc
        conv_mat, conv_bc = construct_convflux_upwind(shape, z_f, z_c, bc, v=self.u)
        grad_mat, grad_bc = construct_grad(shape, z_f, z_c, bc)
        self.div_mat = construct_div(shape, z_f, nu=nu)
        self.conv_mat = conv_mat.tocsc()
        self.conv_bc = np.asarray(conv_bc.todense()).ravel()
        self.diff_flux_mat = (-d_mat @ grad_mat).tocsc()
        self.diff_flux_bc = np.asarray((-d_mat @ grad_bc).todense()).ravel()
        self.jac_diff_rxn = (self.div_mat @ self.diff_flux_mat + k_mat).tocsc()
        self.jac_up = (self.div_mat @ self.conv_mat + self.jac_diff_rxn).tocsc()
        self.g_up = np.asarray(self.div_mat @ (self.conv_bc + self.diff_flux_bc)).ravel()
        self.lu_up = (splu(self.jac_up) if permc_spec is None
                      else splu(self.jac_up, permc_spec=permc_spec))
        self._tvd = {}

    def _tvd_flux(self, limiter):
        if limiter not in self._tvd:
            self._tvd[limiter] = ConvFluxTVD(self.shape, self.z_f, self.z_c, self.bc,
                                             tvd_limiter=limiter)
        return self._tvd[limiter]

    def _solve(self, limiter, tol, jacobian="exact", maxfev=50, residual_form="correction"):
        """Upwind solution, then the TVD correction to it by Newton."""
        c_up = self.lu_up.solve(-self.g_up)
        if limiter is None:
            return c_up, 0
        tvd = self._tvd_flux(limiter)
        g_diff = np.asarray(self.div_mat @ self.diff_flux_bc).ravel()

        def residual(x):
            c = x.reshape(-1)
            flux, jac_conv = tvd(x.reshape(self.shape), self.u)
            if residual_form == "correction":
                # upwind residual of c_up (zero to round-off) taken out exactly
                flux_up = self.conv_mat @ c + self.conv_bc
                g = (self.jac_up @ (c - c_up)
                     + self.div_mat @ (flux.reshape(-1) - flux_up))
            else:                                   # the full residual of c
                g = self.div_mat @ flux.reshape(-1) + self.jac_diff_rxn @ c + g_diff
            if jacobian == "upwind":
                return g, self.jac_up
            return g, (self.div_mat @ jac_conv + self.jac_diff_rxn).tocsc()

        sol = newton(residual, c_up.reshape(self.shape), tol=tol, damping="monotone",
                     stagnation=True, maxfev=maxfev)
        if not sol.success:
            raise RuntimeError(f"TVD Newton solve failed: {sol.message}")
        return np.asarray(sol.x).ravel(), int(sol.nit)

    def _solve_deferred(self, limiter, tol, max_it):
        """Deferred correction: the limiter correction lagged, one LU reused.

        Raises RuntimeError when ``max_it`` sweeps do not bring the largest change
        of a cell value below ``tol``.
        """
        c = self.lu_up.solve(-self.g_up)
        change = np.inf
        for sweep in range(1, max_it + 1):
            _, dc_f = interp_cntr_to_stagg_tvd(c.reshape(self.shape), self.z_f, self.z_c,
                                               self.bc, self.u, tvd_limiter=limiter, axis=0)
            c_new = self.lu_up.solve(-self.g_up - np.asarray(
                self.div_mat @ (self.u * dc_f.reshape(-1, 1))).ravel())
            change = np.max(np.abs(c_new - c))
            c = c_new
            if change < tol:
                return c, sweep
        raise RuntimeError(f"deferred correction did not reach {tol:g} in {max_it} sweeps "
                           f"(last change {change:.2e})")

    def exit_value(self, c):
        """Value on the downstream face, from ``compute_boundary_values`` (pitfalls P8)."""
        return float(compute_boundary_values(c.reshape(self.shape), self.z_f, self.z_c,
                                             self.bc)[2].ravel()[0])

    def system(self):
        """The assembled upwind system ``(jac_up, g_up)``; its residual is ``jac_up c + g_up``."""
        return self.jac_up, self.g_up


class ClosedVessel(_Linear1D):
    """Danckwerts' closed vessel, ``0 <= z <= 1``, on the pymrm operators.

    Steady ``d(u c)/dz - (1/Pe) d2c/dz2 + R c = 0`` with the Danckwerts inlet
    ``u c - D dc/dz = u c*`` (``c* = 1``) and the outlet ``dc/dz = 0``. The
    operators are assembled once in the constructor; :meth:`solve` only solves.

    Parameters
    ----------
    Pe : float
        Peclet number ``uL/D`` of the bed.
    R : float
        Damkoehler number ``kL/u``.
    n : int, optional
        Number of equal cells.
    nu : int, optional
        Geometry of ``construct_div``: 0 Cartesian (the model); other values
        exist for the page's break table.
    inlet : {"danckwerts", "dirichlet"}, optional
        Danckwerts' flux condition, or the naive ``c(0) = c*``.
    inlet_a : float, optional
        Overrides ``a`` of the Danckwerts inlet dictionary (default ``D/u``);
        for the break table.
    outlet : {"zero-gradient", "dirichlet-zero"}, optional
        Danckwerts' ``dc/dz = 0``, or ``c(1) = 0`` for the break table.
    outlet_a, outlet_d : float, optional
        ``a`` and ``d`` of the zero-gradient dictionary (default 1 and 0); the
        page uses them to show that the sign of ``a`` is absent from the
        problem when ``b = d = 0``, and present once ``d != 0``.
    inlet_d : float, optional
        ``d`` of the inlet dictionary, the feed fraction ``c*`` (default 1);
        other values exist for the break table.
    r_sign : float, optional
        Sign of the reaction term; ``-1`` exists for the break table.
    permc_spec : str, optional
        Column ordering handed to ``splu`` for the upwind operator (default
        SuperLU's); the page uses a second ordering to measure how far the
        solution depends on the factorisation.

    Attributes
    ----------
    bc : tuple of dict
        The two boundary dictionaries.
    jac_up : scipy.sparse.csc_array
        Upwind convection, diffusion and reaction, ``div @ (conv - D grad) + R I``.
    g_up : numpy.ndarray
        Boundary contribution to the upwind system, shape ``(n,)``.
    """

    def __init__(self, Pe, R, n=800, nu=0, inlet="danckwerts", inlet_a=None, inlet_d=1.0,
                 outlet="zero-gradient", outlet_a=None, outlet_d=0.0, r_sign=1.0,
                 permc_spec=None):
        self.Pe, self.R, self.n = Pe, R, n
        self.u = 1.0
        self.d_ax = 1.0 / Pe
        self.shape = (n, 1)                      # (cells, fields), never a bare (n,)
        self.z_f = np.linspace(0.0, 1.0, n + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        if inlet == "danckwerts":
            # inlet at z = 0, Danckwerts eq. 31: u c - D dc/dz = u c*, c* = 1.
            # n = -z, so dc/dn = -dc/dz and (D/u) dc/dn + c = c*  ->  a=D/u, b=1, d=1
            bc_in = {"a": self.d_ax / self.u if inlet_a is None else inlet_a,
                     "b": 1.0, "d": inlet_d}
        elif inlet == "dirichlet":
            # the naive alternative, c(0) = c*  ->  a=0, b=1, d=1
            bc_in = {"a": 0.0, "b": 1.0, "d": inlet_d}
        else:
            raise ValueError(f"unknown inlet {inlet!r}")
        if outlet == "zero-gradient":
            # outlet at z = 1, Danckwerts eq. 32: dc/dz = 0. n = +z  ->  a=1, b=0, d=0
            bc_out = {"a": 1.0 if outlet_a is None else outlet_a, "b": 0.0, "d": outlet_d}
        elif outlet == "dirichlet-zero":
            # c(1) = 0  ->  a=0, b=1, d=0 (a defect)
            bc_out = {"a": 0.0, "b": 1.0, "d": 0.0}
        else:
            raise ValueError(f"unknown outlet {outlet!r}")
        self.bc = (bc_in, bc_out)
        d_mat = self.d_ax * eye_array(n + 1, format="csc")       # constant D on every face
        k_mat = r_sign * R * eye_array(n, format="csc")
        self._assemble(d_mat, k_mat, nu=nu, permc_spec=permc_spec)  # nu=0: Cartesian slab

    def solve(self, limiter=vanleer, tol=1e-12, jacobian="exact", maxfev=50,
              residual="correction"):
        """Steady profile; TVD by Newton with the exact Jacobian.

        Parameters
        ----------
        limiter : callable or None, optional
            TVD limiter (default van Leer). None gives first-order upwind,
            solved as one linear system.
        tol : float, optional
            Newton step tolerance on order-one fractions.
            ``stagnation=True`` ends the solve with success should the steps
            reach the round-off floor first.
        jacobian : {"exact", "upwind"}, optional
            ``"exact"`` uses the limited-flux Jacobian of ``ConvFluxTVD``;
            ``"upwind"`` hands Newton the upwind operator instead, which is
            what deferred correction iterates with (pitfalls P10).
        maxfev : int, optional
            Iteration cap of ``newton``.
        residual : {"correction", "full"}, optional
            How the residual is evaluated (module Notes): as the correction to
            the upwind solution (default), or in full, with its higher
            round-off floor.

        Returns
        -------
        Solution

        Raises
        ------
        RuntimeError
            When Newton does not converge.
        """
        c, nit = self._solve(limiter, tol, jacobian, maxfev, residual)
        return Solution(self.z_c, c, self.exit_value(c), nit)

    def solve_deferred(self, limiter=vanleer, tol=1e-14, max_it=80):
        """The same discrete equations by deferred correction (a second route).

        The limiter correction is lagged in the right-hand side and the upwind
        operator is factorised once; the iteration is a fixed point, not
        Newton (pitfalls P10).

        Parameters
        ----------
        limiter : callable, optional
            TVD limiter.
        tol : float, optional
            Stop when the largest change of a cell value in one sweep is below it.
        max_it : int, optional
            Cap on the number of sweeps.

        Returns
        -------
        Solution
            ``nit`` is the number of sweeps.

        Raises
        ------
        RuntimeError
            When ``max_it`` sweeps do not reach ``tol``.
        """
        c, sweeps = self._solve_deferred(limiter, tol, max_it)
        return Solution(self.z_c, c, self.exit_value(c), sweeps)


def geom_faces(length, dz0, ratio=1.06):
    """Faces ``0 .. length`` with cells growing geometrically from ``dz0``.

    Parameters
    ----------
    length : float
        Length to cover.
    dz0 : float
        First cell width (before the final rescaling to ``length``).
    ratio : float, optional
        Growth ratio of consecutive cells.

    Returns
    -------
    numpy.ndarray
        Face positions, starting at 0 and ending at ``length``.
    """
    n = int(np.ceil(np.log1p(length * (ratio - 1) / dz0) / np.log(ratio)))
    f = np.concatenate([[0.0], np.cumsum(dz0 * ratio ** np.arange(n))])
    return f * (length / f[-1])


class ThreeSection(_Linear1D):
    """Wehner & Wilhelm's fore | reaction | after problem as one pymrm domain.

    The grid runs from far upstream to far downstream, with diffusivity
    ``1/Pe_a`` for ``z < 0``, ``1/Pe_b`` in the bed and ``1/Pe_c`` for
    ``z > 1``, and the reaction only in the bed. The only boundary conditions
    are Wehner & Wilhelm's Eqs. 8 and 15 at the two truncations; the bed ends
    are ordinary interior faces.

    Parameters
    ----------
    Pe_a, Pe_b, Pe_c : float
        Peclet groups of the fore, reaction and after sections.
    R : float
        Damkoehler number of the bed.
    n_b : int, optional
        Number of equal cells in the bed.
    P : float, optional
        Length of each end section, in its own decay lengths ``1/Pe``.
    P_after : float, optional
        Overrides ``P`` for the after section alone.
    stretch : bool, optional
        Geometrically stretched end sections (default) or uniform ones with the
        bed's cell size.
    face_average : {"harmonic", "arithmetic"}, optional
        Face diffusivity at the two jumps. Harmonic is the model (pitfalls P2);
        arithmetic exists to measure what it costs.
    react_past_bed : bool, optional
        Let the reaction run on into the after section (a defect for the
        break table: the after section is then no longer flat).

    Attributes
    ----------
    z_f, z_c : numpy.ndarray
        Faces and cell centres.
    d_face : numpy.ndarray
        Face diffusivities, shape ``(n + 1,)``.
    """

    def __init__(self, Pe_a, Pe_b, Pe_c, R, n_b=400, P=30.0, P_after=None, stretch=True,
                 face_average="harmonic", react_past_bed=False):
        self.u = 1.0
        dz0 = 1.0 / n_b
        p_after = P if P_after is None else P_after
        if stretch:
            z_up = -geom_faces(P / Pe_a, dz0)[::-1]
            z_dn = geom_faces(p_after / Pe_c, dz0) + 1.0
        else:
            z_up = np.linspace(-P / Pe_a, 0.0, int(P / Pe_a / dz0) + 1)
            z_dn = np.linspace(1.0, 1.0 + p_after / Pe_c, int(p_after / Pe_c / dz0) + 1)
        self.z_f = np.concatenate([z_up[:-1], np.linspace(0.0, 1.0, n_b + 1)[:-1], z_dn])
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.n = len(self.z_c)
        self.shape = (self.n, 1)
        bed = (self.z_c > 0) & (self.z_c < 1.0)
        d_cell = np.where(bed, 1 / Pe_b, np.where(self.z_c < 0, 1 / Pe_a, 1 / Pe_c))
        d_l, d_r = d_cell[:-1], d_cell[1:]
        self.d_face = np.empty(self.n + 1)
        if face_average == "harmonic":
            # two cells of different D are two resistances in series; the cells on
            # either side of each jump have the same width, so the unweighted
            # harmonic mean is the series resistance (pitfalls P2)
            self.d_face[1:-1] = 2 * d_l * d_r / (d_l + d_r)
        elif face_average == "arithmetic":
            self.d_face[1:-1] = 0.5 * (d_l + d_r)
        else:
            raise ValueError(f"unknown face_average {face_average!r}")
        self.d_face[0], self.d_face[-1] = d_cell[0], d_cell[-1]
        self.k_cell = np.where(self.z_c > 0 if react_past_bed else bed, R, 0.0)
        # W&W eq. 8, f(-inf) = 1: Dirichlet at the upstream truncation -> a=0, b=1, d=1
        # W&W eq. 15, f(+inf) finite: no growing mode; with no reaction downstream
        #   this is df/dz = 0 at any truncation. n = +z -> a=1, b=0, d=0
        self.bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        d_mat = construct_coefficient_matrix(self.d_face.reshape(-1, 1), shape=self.shape,
                                             axis=0)
        k_mat = construct_coefficient_matrix(self.k_cell.reshape(-1, 1), shape=self.shape)
        self._assemble(d_mat, k_mat, nu=0)                      # nu=0: Cartesian

    def solve(self, tol=1e-12, limiter=vanleer, residual="correction"):
        """Steady profile, van Leer TVD by Newton with the exact Jacobian.

        Parameters
        ----------
        tol : float, optional
            Newton step tolerance on order-one fractions.
        limiter : callable or None, optional
            van Leer (the model); None gives first-order upwind, one linear
            solve, which exists for the page's break table.
        residual : {"correction", "full"}, optional
            As for :meth:`ClosedVessel.solve`.

        Returns
        -------
        z_c : numpy.ndarray
            Cell centres.
        c : numpy.ndarray
            Fraction remaining at the cell centres.
        z_f : numpy.ndarray
            Faces.
        flux : numpy.ndarray
            Total (convective plus diffusive) flux on every face.

        Raises
        ------
        RuntimeError
            When Newton does not converge.
        """
        c, self.nit = self._solve(limiter, tol, residual_form=residual)
        if limiter is None:
            conv_flux = self.conv_mat @ c + self.conv_bc
        else:
            conv_flux = np.asarray(self._tvd_flux(limiter)(c.reshape(self.shape), self.u)[0])
        flux = np.ravel(conv_flux) + self.diff_flux_mat @ c + self.diff_flux_bc
        return self.z_c, c, self.z_f, flux

    def solve_deferred(self, tol=1e-11, max_it=80):
        """The same discrete equations by deferred correction (a second route).

        Parameters
        ----------
        tol : float, optional
            Stop when the largest change of a cell value in one sweep is below it.
        max_it : int, optional
            Cap on the number of sweeps.

        Returns
        -------
        z_c, c : numpy.ndarray
            Cell centres and the fraction remaining there.

        Raises
        ------
        RuntimeError
            When ``max_it`` sweeps do not reach ``tol``.
        """
        c, self.sweeps = self._solve_deferred(vanleer, tol, max_it)
        return self.z_c, c
