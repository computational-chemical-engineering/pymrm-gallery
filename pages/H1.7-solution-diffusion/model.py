"""Wijmans and Baker (1995): the feed-side concentration-polarisation film.

Model code of gallery page H1.7. The membrane itself is algebra (Eqs. 37 and
40 of the paper are closed forms) and is fitted on the page. What needs a
transport solve is the feed-side film the derivation assumes away: salt is
carried to the membrane at the permeate velocity, a fraction ``1 - R_int``
passes, and the rest diffuses back against the flow.

Contents
--------
film_wall_ratio
    ``c_wall / c_bulk`` from a pymrm steady convection-diffusion solve.
film_analytic
    The closed form of the same balance.
polarised_fluxes
    Water fluxes ``J = A (p - pi_bulk * c_wall/c_bulk)`` with the film in the
    loop, one bracketed root per pressure.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from pymrm import (compute_boundary_values, construct_convflux_upwind,
                   construct_div, construct_grad)

__all__ = ["film_wall_ratio", "film_analytic", "polarised_fluxes"]


def film_wall_ratio(Jv, delta, D, R_int, n=400, profile=False, wall="face", solver=None):
    """Wall-to-bulk salt ratio in the feed-side film, from a pymrm solve.

    Solves ``d/dx (v c - D dc/dx) = 0`` on ``0 < x < delta`` with the bulk
    value at ``x = 0`` and, at the membrane ``x = delta``, the salt balance
    ``v c - D dc/dx = v (1 - R_int) c``: what arrives by convection minus
    what permeates leaves by diffusion. Upwind convection, so first order.

    Parameters
    ----------
    Jv : float
        Permeate (water) velocity toward the membrane (m/s), >= 0.
    delta : float
        Film thickness (m), >= 0.
    D : float
        Salt diffusivity in the film (m2/s).
    R_int : float
        Intrinsic rejection of the membrane, ``1 - c_permeate/c_wall``.
    n : int, optional
        Number of cells.
    profile : bool, optional
        Return ``(x_c, c)`` instead of the wall value.
    wall : {"face", "cell"}, optional
        Where the wall value is read: on the membrane face with
        ``compute_boundary_values`` (the value the operators transport), or
        from the last cell centre, half a cell short (pitfall P8; kept for
        the page's break table).
    solver : callable, optional
        ``solver(matrix, rhs) -> c`` for the assembled sparse system;
        ``scipy.sparse.linalg.spsolve`` when None. A second solver gives the
        round-off floor of the solve.

    Returns
    -------
    float or tuple of numpy.ndarray
        ``c_wall / c_bulk``, or the cell centres and normalised profile.
    """
    if delta <= 0.0 or Jv <= 0.0:
        return (np.array([0.0]), np.array([1.0])) if profile else 1.0
    x_f = np.linspace(0.0, delta, n + 1)
    x_c = 0.5 * (x_f[:-1] + x_f[1:])
    shape = (n, 1)
    # bulk, x = 0: c = 1 (normalised)
    # membrane, x = delta, outward normal +x: D dc/dx - Jv R_int c = 0
    bc = ({"a": 0.0, "b": 1.0, "d": 1.0},
          {"a": D, "b": -Jv * R_int, "d": 0.0})
    grad, grad_bc = construct_grad(shape, x_f, x_c, bc)
    conv, conv_bc = construct_convflux_upwind(shape, x_f, x_c, bc, v=Jv)
    div = construct_div(shape, x_f, nu=0)          # nu = 0: Cartesian film
    lhs = div @ (conv - D * grad)
    rhs = -(div @ (conv_bc - D * grad_bc)).toarray().ravel()
    c = (spsolve if solver is None else solver)(lhs.tocsc(), rhs).reshape(shape)
    if profile:
        return x_c, c[:, 0]
    if wall == "cell":
        return float(c[-1, 0])
    _, _, c_wall, _ = compute_boundary_values(c, x_f, x_c, bc)
    return float(np.ravel(c_wall)[0])


def film_analytic(Jv, delta, D, R_int):
    """Closed form of the film balance of :func:`film_wall_ratio`.

    ``c_wall/c_bulk = exp(Pe) / (R_int + (1 - R_int) exp(Pe))`` with
    ``Pe = Jv delta / D``.
    """
    Pe = Jv * delta / D
    return np.exp(Pe) / (R_int + (1.0 - R_int) * np.exp(Pe))


def polarised_fluxes(p, A, pi_bulk, wall_ratio):
    """Water fluxes with the wall osmotic pressure raised by the film.

    Solves ``J = A (p - pi_bulk * wall_ratio(J))`` at each pressure. The
    residual ``J - A (p - pi_bulk w(J))`` is increasing in ``J`` when ``w`` is
    non-decreasing and at least 1, negative at ``J = 0`` when ``p > pi_bulk``
    and positive at ``J = A p``, so ``[0, A p]`` brackets exactly one root.

    Parameters
    ----------
    p : array_like
        Applied pressures (same unit as ``pi_bulk``), each above ``pi_bulk``.
    A : float
        Water permeability (flux per unit pressure).
    pi_bulk : float
        Osmotic pressure of the bulk feed.
    wall_ratio : callable
        ``wall_ratio(J) -> c_wall/c_bulk`` for a flux ``J`` in the unit of
        ``A * p``.

    Returns
    -------
    numpy.ndarray
        The fluxes.

    Raises
    ------
    ValueError
        If a pressure does not exceed ``pi_bulk`` (no positive flux), or if
        the root-finder does not converge.
    """
    p = np.asarray(p, dtype=float)
    if np.any(p <= pi_bulk):
        raise ValueError("every pressure must exceed the bulk osmotic pressure")
    out = []
    for pv in p:
        J, r = brentq(lambda J: J - A * (pv - pi_bulk * wall_ratio(J)), 0.0, A * pv,
                      xtol=1e-12, full_output=True)
        if not r.converged:
            raise ValueError(f"flux root not converged at p = {pv}")
        out.append(J)
    return np.array(out)
