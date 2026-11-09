"""Models of page A3.6 (Calderbank and Moo-Young 1961).

Two linear pymrm problems and the closed forms they are checked against:

* :func:`sh_sphere`, steady diffusion from a sphere into a stagnant shell, the
  stagnant floor ``Sh -> 2`` of the paper's page-44 form;
* :class:`BoundaryLayer`, the concentration boundary layer
  ``u(y) dc/dx = D d2c/dy2`` with ``u = u_s + gdot y``, marched in ``x`` with
  backward Euler on one factorised operator;
* :func:`leveque` and :func:`higbie`, the mean coefficients of its two pure
  limits.

Both problems are linear, so each is one sparse solve (the sphere) or one LU
factorisation reused at every step (the march). Units are those of the caller;
only ratios matter.
"""

import numpy as np
from scipy.sparse.linalg import splu, spsolve
from scipy.special import gamma as gammafn
from pymrm import (compute_boundary_values, construct_coefficient_matrix, construct_div,
                   construct_grad)


def stretched_faces(lo, hi, n, stretch):
    """Faces ``lo .. hi`` refined toward ``lo`` by a fixed exponential map.

    Parameters
    ----------
    lo, hi : float
        Ends of the domain.
    n : int
        Number of cells.
    stretch : float
        Map parameter, nonzero.

    Returns
    -------
    numpy.ndarray
        ``n + 1`` face positions. The map is a fixed function of ``i / n``, so
        doubling ``n`` nests the grids (every old face is a new face).
    """
    w = np.linspace(0.0, 1.0, n + 1)
    return lo + (hi - lo) * np.expm1(stretch * w) / np.expm1(stretch)


def sh_sphere(ratio=1e-4, n_r=800, nu=2, outer="dirichlet", at="surface"):
    """Sherwood number of a sphere in a stagnant shell, ``Sh = k d / D``.

    Solves ``(1/r^2) d/dr (r^2 dc/dr) = 0`` (for ``nu = 2``) on ``R <= r <= R/ratio``
    with one sparse solve: the problem is linear.

    Parameters
    ----------
    ratio : float
        ``R / R_inf``, the sphere radius over the shell radius.
    n_r : int
        Number of cells, spaced geometrically (``r_f = R (R_inf/R)^(i/n_r)``),
        so the grids of a doubling ladder nest.
    nu : int
        Geometry of ``construct_div``: 0 slab, 1 cylinder, 2 sphere (correct).
    outer : {"dirichlet", "noflux"}
        Clean bulk ``c = 0`` at ``R_inf``, or a closed shell (a break row).
    at : {"surface", "shell"}
        Where the flux is read: at the sphere, or at the shell face scaled by
        the area ratio ``(R_inf/R)^2``. For ``nu = 2`` the two agree to the
        residual of the solve (the discrete divergence telescopes).

    Returns
    -------
    float
        ``Sh``, with the surface gradient read by ``compute_boundary_values``
        under the solve's own boundary conditions. Exact for ``nu = 2``:
        ``2 / (1 - ratio)``.
    """
    R, R_inf = 1.0, 1.0 / ratio
    r_f = R * (R_inf / R)**np.linspace(0.0, 1.0, n_r + 1)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n_r, 1)
    # Outward normal n, a dc/dn + b c = d:
    #   r = R     : c = c* = 1 (saturated surface)       -> a=0, b=1, d=1
    #   r = R_inf : c = 0 (clean bulk)                    -> a=0, b=1, d=0
    #               or dc/dn = 0 (closed shell, a defect) -> a=1, b=0, d=0
    bc = ({"a": 0.0, "b": 1.0, "d": 1.0},
          {"a": 0.0, "b": 1.0, "d": 0.0} if outer == "dirichlet"
          else {"a": 1.0, "b": 0.0, "d": 0.0})
    div = construct_div(shape, r_f, nu=nu)
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc)
    lap = (div @ grad).tocsc()
    rhs = -np.asarray((div @ grad_bc).todense()).ravel()
    c = spsolve(lap, rhs).reshape(shape)
    if at == "surface":
        _, dcdr = compute_boundary_values(c, r_f, r_c, bc[0], bound_id=0)
        area = 1.0
    else:
        _, dcdr = compute_boundary_values(c, r_f, r_c, bc[1], bound_id=1)
        area = (R_inf / R)**2
    return float(-2.0 * R * area * np.ravel(dcdr)[0])    # d = 2R, D = 1, c* = 1


class BoundaryLayer:
    """Concentration boundary layer over an interface of contact length ``L``.

    ``u(y) dc/dx = D d2c/dy2`` with ``u = u_s + gdot y`` and ``v = 0`` (an exact
    solution of continuity), ``c(x, 0) = 1``, ``c(0, y) = 0`` and no flux at
    ``y_max``. ``u_s = 0`` is the rigid interface (Leveque, ``k ~ D^(2/3)``),
    ``gdot = 0`` the mobile one (Higbie, ``k ~ D^(1/2)``).

    The problem is linear and marched in ``x`` with backward Euler: one
    constant operator, factorised once per step length.

    Parameters
    ----------
    D : float
        Diffusivity.
    u_s, gdot : float
        Interface velocity and shear rate.
    L : float
        Contact length.
    n_y : int
        Cells across the layer, stretched toward the interface.
    stretch : float
        Grid stretching parameter.
    y_span : float
        Domain height in units of the smaller of the two layer thicknesses
        ``(9 D L / gdot)^(1/3)`` and ``sqrt(D L / u_s)``.
    """

    def __init__(self, D=1e-9, u_s=0.0, gdot=1.0, L=1.0, n_y=400, stretch=6.0,
                 y_span=16.0):
        self.D, self.u_s, self.gdot, self.L = D, u_s, gdot, L
        d_lev = (9 * D * L / gdot)**(1 / 3) if gdot > 0 else np.inf
        d_hig = np.sqrt(D * L / u_s) if u_s > 0 else np.inf
        self.y_max = y_span * min(d_lev, d_hig)
        self.n_y = n_y
        self.y_f = stretched_faces(0.0, self.y_max, n_y, stretch)
        self.y_c = 0.5 * (self.y_f[:-1] + self.y_f[1:])
        self.shape = (n_y, 1)
        self.V = np.diff(self.y_f)
        self.u = u_s + gdot * self.y_c

    def bc(self, interface="dirichlet"):
        """Boundary conditions, outward normal, ``a dc/dn + b c = d``.

        y = 0 (normal -y): saturated interface c = 1 -> a=0, b=1, d=1;
        or, as a defect, zero flux dc/dn = 0 -> a=1, b=0, d=0.
        y = y_max (normal +y): zero flux dc/dn = 0 -> a=1, b=0, d=0.
        """
        return ({"a": 0.0, "b": 1.0, "d": 1.0} if interface == "dirichlet"
                else {"a": 1.0, "b": 0.0, "d": 0.0},
                {"a": 1.0, "b": 0.0, "d": 0.0})

    def run(self, n_x=400, interface="dirichlet", read="face"):
        """March from ``x = 0`` to ``x = L`` in ``n_x`` backward-Euler steps.

        Parameters
        ----------
        n_x : int
            Number of steps.
        interface : {"dirichlet", "noflux"}
            The interface condition (``"noflux"`` is a defect).
        read : {"face", "centre"}
            How the interface flux is read: the face gradient of
            ``compute_boundary_values`` under the solve's own boundary
            conditions, or (a defect) the two-point difference
            ``(1 - c_0) / y_0`` to the first cell centre.

        Returns
        -------
        dict
            ``kbar`` the mean coefficient from the step fluxes (the sum of
            ``dx`` times the face flux, read with ``compute_boundary_values``,
            over ``L``); ``kbar_balance`` the same from the exit-plane balance
            ``sum(u V c) / L``; ``k_local`` the coefficient of the last step;
            ``c`` the exit profile, shape ``(n_y, 1)``.
        """
        bc = self.bc(interface)
        div = construct_div(self.shape, self.y_f, nu=0)
        grad, grad_bc = construct_grad(self.shape, self.y_f, self.y_c, bc)
        A = (self.D * (div @ grad)).tocsc()
        b = self.D * np.asarray((div @ grad_bc).todense()).ravel()
        U = construct_coefficient_matrix(self.u.reshape(-1, 1), self.shape).tocsc()
        dx = self.L / n_x
        lu = splu((U / dx - A).tocsc())        # constant operator, factorised once
        c = np.zeros(self.n_y)
        kx = np.empty(n_x)
        for j in range(n_x):
            c = lu.solve(self.u * c / dx + b)
            if read == "face":
                _, g0 = compute_boundary_values(c.reshape(self.shape), self.y_f, self.y_c,
                                                bc[0], bound_id=0)
                kx[j] = -self.D * float(np.ravel(g0)[0])  # flux into the fluid, c* = 1
            else:
                kx[j] = self.D * (1.0 - c[0]) / self.y_c[0]
        return {"kbar": float(kx.sum() * dx / self.L),
                "kbar_balance": float((self.u * self.V * c).sum() / self.L),
                "k_local": float(kx[-1]), "c": c.reshape(self.shape)}


def leveque(D, gdot, L):
    """Mean coefficient over ``0..L`` for ``u = gdot y`` (Leveque similarity).

    ``f'' + 3 eta^2 f' = 0`` gives ``k(x) = D / (Gamma(4/3) (9 D x / gdot)^(1/3))``
    and ``kbar = (3/2) k(L)``.
    """
    return 1.5 * D**(2 / 3) * (gdot / (9 * L))**(1 / 3) / gammafn(4 / 3)


def higbie(D, u_s, L):
    """Mean coefficient over ``0..L`` for ``u = u_s`` (penetration, contact time ``L/u_s``)."""
    return 2.0 * np.sqrt(D * u_s / (np.pi * L))
