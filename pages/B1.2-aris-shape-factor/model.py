"""Linear first-order reaction and diffusion in a pellet of arbitrary shape.

Model code of gallery page B1.2 (Aris 1957, the generalised shape modulus).
The page imports these. The closed forms and eigenfunction series it compares
with are written out in the notebook, so the pymrm solve and its references
never share a line of code.

Contents
--------
SEALED, EXPOSED
    The two boundary dictionaries every pellet uses.
graded
    Face coordinates on one axis, optionally clustered toward the exposed end.
Pellet
    ``div(-grad chi) + lambda^2 chi = 0`` on 1 to 3 axes, ``chi = 1`` on the
    exposed faces, with a geometry factor ``nu`` on axis 0.
slab, long_cylinder, sphere, spherical_shell, finite_cylinder, square_rod, cube
    The shapes of the page, each returned with its ``(v_p/s_x)/a``.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import spsolve

from pymrm import (compute_boundary_values, construct_coefficient_matrix,
                   construct_div, construct_grad)

__all__ = ["SEALED", "EXPOSED", "graded", "Pellet", "slab", "long_cylinder", "sphere",
           "spherical_shell", "finite_cylinder", "square_rod", "cube"]

# a dchi/dn + b chi = d with the OUTWARD normal n (pitfalls P1)
SEALED = {"a": 1.0, "b": 0.0, "d": 0.0}    # dchi/dn = 0: symmetry plane, or a sealed face
EXPOSED = {"a": 0.0, "b": 1.0, "d": 1.0}   # chi = 1: face exposed to the bulk


def graded(lo, hi, n, q=1.0):
    """Face coordinates on ``[lo, hi]``.

    Parameters
    ----------
    lo, hi : float
        Ends of the axis; ``hi`` is the exposed face.
    n : int
        Number of cells.
    q : float, optional
        Grading exponent: ``x = lo + (hi - lo) (1 - (1 - u)^q)`` on uniform
        ``u``. ``q = 1`` is uniform; ``q > 1`` clusters cells at ``hi``. The
        map is nested: doubling ``n`` keeps every face.

    Returns
    -------
    numpy.ndarray
        ``n + 1`` increasing face coordinates.
    """
    u = np.linspace(0.0, 1.0, n + 1)
    return lo + (hi - lo) * (1.0 - (1.0 - u) ** q)


class Pellet:
    """Normalised concentration ``chi`` in a pellet: ``div(-grad chi) + lambda^2 chi = 0``.

    Each axis runs from a sealed or symmetry face at its low end to a face
    exposed to the bulk (``chi = 1``) at its high end. Axis 0 carries the
    geometry factor ``nu`` (0 Cartesian, 1 cylindrical radial, 2 spherical
    radial); every other axis is Cartesian.

    Parameters
    ----------
    n : tuple of int
        Cells per axis.
    extent : tuple of float
        High end of each axis, in units of the characteristic length ``a``.
    nu : tuple of int
        Geometry factor per axis, passed to ``construct_div``.
    origin : float, optional
        Low end of axis 0 (``p`` for a spherical shell between ``pa`` and ``a``).
    q : float, optional
        Grading exponent of every axis, see :func:`graded`.
    bc_low : dict, optional
        Boundary dictionary of the low end of every axis. Default
        :data:`SEALED`. Exists for the page's break table (a shell whose inner
        face is exposed by mistake).

    Attributes
    ----------
    x_f, x_c : list of numpy.ndarray
        Faces and cell centres per axis.
    jac_diff : scipy.sparse.csc_array
        Constant operator ``div @ (-grad)``.
    g_diff : numpy.ndarray
        Its constant boundary contribution, flattened.
    jac_react : scipy.sparse.csc_array
        Unit reaction operator; the system is ``jac_diff + lambda^2 jac_react``.
    volume : numpy.ndarray
        Cell volumes (per unit angle and unit length of the omitted axes).
    """

    def __init__(self, n, extent, nu, origin=0.0, q=1.0, bc_low=SEALED):
        self.shape, self.nu, self.bc = tuple(n), tuple(nu), (bc_low, EXPOSED)
        self.x_f, self.x_c = [], []
        jac_diff, g_diff, measures = 0, 0, []
        for ax, (n_ax, ext_ax, nu_ax) in enumerate(zip(n, extent, nu)):
            x_f = graded(origin if ax == 0 else 0.0, ext_ax, n_ax, q)
            x_c = 0.5 * (x_f[:-1] + x_f[1:])
            self.x_f.append(x_f)
            self.x_c.append(x_c)
            # low end: dchi/dn = 0 (SEALED); high end: chi = 1 (EXPOSED)
            grad_mat, grad_bc = construct_grad(self.shape, x_f, x_c, self.bc, axis=ax)
            # nu_ax = 0 Cartesian, 1 cylindrical radial, 2 spherical radial
            div_mat = construct_div(self.shape, x_f, nu=nu_ax, axis=ax)
            jac_diff = jac_diff + div_mat @ (-grad_mat)        # flux = -grad(chi)
            g_diff = g_diff + div_mat @ (-grad_bc)
            measures.append(np.diff(x_f ** (nu_ax + 1)) / (nu_ax + 1))
        self.jac_diff = jac_diff.tocsc()
        self.g_diff = g_diff.toarray().ravel()
        self.jac_react = construct_coefficient_matrix(np.ones(self.shape), shape=self.shape)
        self._measures = measures
        volume = measures[0]
        for m in measures[1:]:
            volume = np.multiply.outer(volume, m)
        self.volume = volume

    def solve(self, lam):
        """Cell values of ``chi`` at ``lambda = lam`` (one sparse direct solve).

        Parameters
        ----------
        lam : float
            ``lambda a``, the reaction-diffusion group on the length ``a``.

        Returns
        -------
        numpy.ndarray
            ``chi``, shape :attr:`shape`.
        """
        mat = (self.jac_diff + lam**2 * self.jac_react).tocsc()
        return spsolve(mat, -self.g_diff).reshape(self.shape)

    def eta(self, lam, chi=None):
        """Effectiveness factor as the volume average of ``chi`` (Aris eq. 13).

        Parameters
        ----------
        lam : float
            ``lambda a``.
        chi : numpy.ndarray, optional
            A solution from :meth:`solve`; solved here when omitted.

        Returns
        -------
        float
        """
        chi = self.solve(lam) if chi is None else chi
        return float(np.sum(self.volume * chi) / np.sum(self.volume))

    def eta_flux(self, lam, chi=None, read="face"):
        """Effectiveness factor from the flux through the exposed faces.

        ``eta = (1 / (lambda^2 v_p)) * surface integral of dchi/dn``. The face
        gradient is read with ``compute_boundary_values``, the gradient the
        operator transports (pitfalls P8).

        Parameters
        ----------
        lam : float
            ``lambda a``.
        chi : numpy.ndarray, optional
            A solution from :meth:`solve`; solved here when omitted.
        read : {"face", "cell"}, optional
            ``"cell"`` reads a two-point gradient ``(1 - chi_N) / (x_f - x_N)``
            from the last cell centre instead. Exists for the page's break table.

        Returns
        -------
        float
        """
        chi = self.solve(lam) if chi is None else chi
        total = 0.0
        for ax, (x_f, x_c, nu_ax) in enumerate(zip(self.x_f, self.x_c, self.nu)):
            if read == "face":
                _, grad = compute_boundary_values(chi, x_f, x_c, self.bc[1], axis=ax,
                                                  bound_id=1)
            else:
                last = np.take(chi, [-1], axis=ax)
                grad = (1.0 - last) / (x_f[-1] - x_c[-1])
            area = x_f[-1] ** nu_ax                  # face measure on this axis
            others = [m for k, m in enumerate(self._measures) if k != ax]
            weight = np.ones(())
            for m in others:
                weight = np.multiply.outer(weight, m)
            total += area * float(np.sum(np.reshape(grad, weight.shape) * weight))
        return total / (lam**2 * float(np.sum(self.volume)))


# The shape catalogue. Each returns (pellet, va) with va = (v_p/s_x)/a, which
# converts Aris's modulus Lambda = lambda v_p/s_x to lambda a = Lambda / va.

def slab(n=400, q=1.0):
    """Flat plate of half-thickness ``a``, both faces exposed.

    Parameters
    ----------
    n : int, optional
        Cells across the half-thickness.
    q : float, optional
        Grading exponent, see :func:`graded`.

    Returns
    -------
    pellet : Pellet
        The assembled pellet.
    va : float
        ``(v_p/s_x)/a``; ``lambda a = Lambda / va``.
    """
    return Pellet((n,), (1.0,), (0,), q=q), 1.0


def long_cylinder(n=400, q=1.0):
    """Long circular cylinder of radius ``a``, ends sealed.

    Parameters
    ----------
    n : int, optional
        Radial cells.
    q : float, optional
        Grading exponent, see :func:`graded`.

    Returns
    -------
    pellet : Pellet
        The assembled pellet.
    va : float
        ``(v_p/s_x)/a``; ``lambda a = Lambda / va``.
    """
    return Pellet((n,), (1.0,), (1,), q=q), 0.5


def sphere(n=400, q=1.0):
    """Sphere of radius ``a``.

    Parameters
    ----------
    n : int, optional
        Radial cells.
    q : float, optional
        Grading exponent, see :func:`graded`.

    Returns
    -------
    pellet : Pellet
        The assembled pellet.
    va : float
        ``(v_p/s_x)/a``; ``lambda a = Lambda / va``.
    """
    return Pellet((n,), (1.0,), (2,), q=q), 1.0 / 3.0


def spherical_shell(p, n=400, q=1.0, bc_low=SEALED):
    """Spherical shell between ``pa`` and ``a``; only the outer surface is exposed.

    Parameters
    ----------
    p : float
        Inner to outer radius ratio, ``0 < p < 1``.
    n : int, optional
        Radial cells.
    q : float, optional
        Grading exponent, see :func:`graded`.
    bc_low : dict, optional
        Boundary dictionary of the inner face (default sealed); for the break
        table only.

    Returns
    -------
    pellet : Pellet
        The assembled pellet.
    va : float
        ``(v_p/s_x)/a``; ``lambda a = Lambda / va``.
    """
    return Pellet((n,), (1.0,), (2,), origin=p, q=q, bc_low=bc_low), (1.0 - p**3) / 3.0


def finite_cylinder(tau, n_r=120, n_z=None, q=1.0):
    """Circular cylinder of radius ``a`` and half-length ``tau a``, every face exposed.

    Parameters
    ----------
    tau : float
        Half-length over radius.
    n_r : int, optional
        Radial cells.
    n_z : int, optional
        Axial cells over the half-length; default ``n_r tau`` clipped to
        ``[40, 2 n_r]``.
    q : float, optional
        Grading exponent of both axes, see :func:`graded`.

    Returns
    -------
    pellet : Pellet
        The assembled pellet.
    va : float
        ``(v_p/s_x)/a``; ``lambda a = Lambda / va``.
    """
    if n_z is None:
        n_z = int(np.clip(round(n_r * tau), 40, 2 * n_r))
    return Pellet((n_r, n_z), (1.0, tau), (1, 0), q=q), tau / (2.0 * tau + 1.0)


def square_rod(n=200, q=1.0):
    """Infinitely long rod of square section, half-width ``a``.

    Parameters
    ----------
    n : int, optional
        Cells per axis over the half-width.
    q : float, optional
        Grading exponent, see :func:`graded`.

    Returns
    -------
    pellet : Pellet
        The assembled pellet.
    va : float
        ``(v_p/s_x)/a``; ``lambda a = Lambda / va``.
    """
    return Pellet((n, n), (1.0, 1.0), (0, 0), q=q), 0.5


def cube(n=24, q=1.0):
    """Cube of half-width ``a``: the same ``v_p/s_x`` as a sphere of radius ``a``.

    Parameters
    ----------
    n : int, optional
        Cells per axis over the half-width.
    q : float, optional
        Grading exponent, see :func:`graded`.

    Returns
    -------
    pellet : Pellet
        The assembled pellet.
    va : float
        ``(v_p/s_x)/a``; ``lambda a = Lambda / va``.
    """
    return Pellet((n, n, n), (1.0, 1.0, 1.0), (0, 0, 0), q=q), 1.0 / 3.0
