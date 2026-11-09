"""Models of page A3.4 (Wakao and Funazkri 1978).

The bed model of the paper, its closed form and the inversion that turns a
measured outlet into a Sherwood number:

* :func:`theta_pymrm`, eq. (4) with Danckwerts boundary conditions,
  ``theta' - Pe_L^-1 theta'' + N theta = 0`` on ``0 < z < 1``, discretised with
  pymrm. The problem is linear, so it is one sparse solve;
* :func:`theta_eq7` and :func:`theta_eq8`, the paper's closed form and its
  plug-flow limit;
* :func:`pe_particle`, :func:`theta_of_Sh` and :func:`Sh_of_theta`, the
  dispersion laws and the forward and inverse maps between ``Sh`` and the outlet;
* :func:`eq12`, :func:`eq9_ranz` and :func:`eq11_petrovic`, the three
  correlations the paper compares.

Every function that uses a printed constant takes the dictionary ``C`` of them
(keys in :data:`CONSTANT_KEYS`), so a mis-read constant is one changed entry.
Everything is dimensionless.
"""

import numpy as np
from scipy.sparse import eye_array
from scipy.sparse.linalg import spsolve
from pymrm import (compute_boundary_values, construct_convflux_upwind, construct_div,
                   construct_grad)
from pymrm.checks import find_roots

CONSTANT_KEYS = ("alpha", "beta", "E0", "slope", "pe_pt", "eps_b",
                 "alpha_ranz", "beta_ranz", "jd_pt", "jd_exp_pt")


def theta_pymrm(N, Pe_L, n=1600, bc_in=None, nu=0, disp=1.0, profile=False):
    """Outlet of eq. (4) with Danckwerts boundary conditions, solved with pymrm.

    Parameters
    ----------
    N : float
        Number of transfer units, ``a k_f L / (eps_b U)``.
    Pe_L : float
        Bed Peclet number ``U L / D_ax``.
    n : int, optional
        Number of uniform cells on ``0 < z < 1``.
    bc_in : dict, optional
        Inlet condition; the default is the Danckwerts inlet. Given only to
        inject a defect.
    nu : int, optional
        Geometry code of the divergence; ``0`` is the slab of the bed model.
    disp : float, optional
        Factor on the dispersion term; ``1`` is the model.
    profile : bool, optional
        Also return the cell values, the cell centres and the residual.

    Returns
    -------
    float or tuple
        The outlet ``theta(1) = (C_s - C_exit)/(C_s - C_in)``, read with
        ``compute_boundary_values``; with ``profile`` also the cell values,
        the cell centres and the relative residual of the linear solve.
    """
    shape = (n, 1)
    z_f = np.linspace(0.0, 1.0, n + 1)
    z_c = 0.5 * (z_f[:-1] + z_f[1:])
    # a dtheta/dn + b theta = d, n the OUTWARD normal (dtheta/dn = -dtheta/dz at z = 0).
    #   z = 0: U (C - C_in) = D_ax dC/dx, i.e. theta - Pe_L^-1 dtheta/dz = 1   (eq. 5)
    #   z = 1: dC/dx = 0, i.e. dtheta/dz = 0                                    (eq. 5)
    bc = (bc_in if bc_in is not None else {"a": 1.0 / Pe_L, "b": 1.0, "d": 1.0},
          {"a": 1.0, "b": 0.0, "d": 0.0})
    conv, conv_bc = construct_convflux_upwind(shape, z_f, z_c, bc, v=1.0, axis=0)
    grad, grad_bc = construct_grad(shape, z_f, z_c, bc, axis=0)
    div = construct_div(shape, z_f, nu=nu, axis=0)
    jac = (div @ (conv - disp * grad / Pe_L) + N * eye_array(n, format="csr")).tocsc()
    rhs = -(div @ (conv_bc - disp * grad_bc / Pe_L)).toarray().ravel()
    th = spsolve(jac, rhs)
    out = float(np.ravel(compute_boundary_values(th.reshape(shape), z_f, z_c,
                                                 bc=bc[1], axis=0, bound_id=1)[0])[0])
    if not profile:
        return out
    resid = float(np.abs(jac @ th - rhs).max() / np.abs(rhs).max())
    return out, th, z_c, resid


def theta_eq7(N, Pe_L):
    """Eq. (7), the closed-form outlet of eq. (4) with Danckwerts conditions.

    Written with ``exp(A Pe_L / 2)`` factored out so it stays finite for a long bed.

    Parameters
    ----------
    N, Pe_L : float or numpy.ndarray
        Transfer units and bed Peclet number.

    Returns
    -------
    float or numpy.ndarray
        ``theta`` at the outlet.
    """
    A = np.sqrt(1.0 + 4.0 * N / Pe_L)
    B = 0.5 * Pe_L
    return 4 * A * np.exp(B - A * B) / ((1 + A) ** 2 - (1 - A) ** 2 * np.exp(-2 * A * B))


def theta_eq8(N):
    """Eq. (8), the plug-flow limit ``D_ax -> 0``: ``theta = exp(-N)``."""
    return np.exp(-N)


def pe_particle(Sc_Re, law, C):
    """Particle Peclet number ``d_p U / D_ax`` under one dispersion assumption.

    Parameters
    ----------
    Sc_Re : float
        ``Sc Re``.
    law : {"none", "pe2", "eq2"}
        No dispersion, the constant Peclet number of Petrovic and Thodos, or
        eq. (2), ``eps_b D_ax / D_v = E0 + slope Sc Re``.
    C : dict
        Printed constants.

    Returns
    -------
    float
    """
    if law == "none":
        return np.inf
    if law == "pe2":
        return C["pe_pt"]
    if law == "eq2":
        return Sc_Re / (C["E0"] + C["slope"] * Sc_Re)
    raise ValueError(law)


def n_units(Sh, Re, Sc, L_dp, C):
    """``N = a k_f L / (eps_b U)`` with ``a = 6 (1 - eps_b) / d_p`` for spheres."""
    return 6.0 * (1.0 - C["eps_b"]) * Sh * L_dp / (Sc * Re)


def theta_of_Sh(Sh, Re, Sc, L_dp, law, C):
    """The outlet a bed with Sherwood number ``Sh`` gives under ``law`` (eq. 7 or 8)."""
    N = n_units(Sh, Re, Sc, L_dp, C)
    pe = pe_particle(Sc * Re, law, C)
    return theta_eq8(N) if not np.isfinite(pe) else theta_eq7(N, pe * L_dp)


def Sh_of_theta(theta, Re, Sc, L_dp, law, C, lo=1e-6, hi=1e6):
    """The Sherwood number that reproduces the outlet ``theta`` under ``law``.

    The plug-flow inverse is closed form. Otherwise the root of
    ``theta_eq7(Sh) - theta`` over ``lo <= Sh <= hi`` from
    ``pymrm.checks.find_roots``; it raises unless there is exactly one.

    Returns
    -------
    float
    """
    k = 6.0 * (1.0 - C["eps_b"]) * L_dp / (Sc * Re)
    pe = pe_particle(Sc * Re, law, C)
    if not np.isfinite(pe):
        return -np.log(theta) / k
    r = find_roots(lambda s: theta_eq7(s * k, pe * L_dp) - theta, lo, hi,
                   n_scan=48, log=True, xtol=1e-14)
    if len(r["roots"]) != 1:
        raise ValueError(f"Sh_of_theta: {len(r['roots'])} roots ({r['message']})")
    return float(r["roots"][0])


def eq12(Re, Sc, C):
    """Eq. (12): ``Sh = 2 + alpha Sc^(1/3) Re^beta``."""
    return 2.0 + C["alpha"] * Sc ** (1 / 3) * np.asarray(Re, float) ** C["beta"]


def eq9_ranz(Re, Sc, C):
    """Eq. (9), Ranz and Marshall: ``Sh = 2 + 0.6 Sc^(1/3) Re^(1/2)``."""
    return 2.0 + C["alpha_ranz"] * Sc ** (1 / 3) * np.asarray(Re, float) ** C["beta_ranz"]


def eq11_petrovic(Re, Sc, C):
    """Eq. (11), ``eps_b j_D = 0.357 Re^-0.359``, as a Sherwood number.

    With ``j_D = (k_c / U_sup) Sc^(2/3)``, ``Sh = j_D Re Sc^(1/3)``.
    """
    Re = np.asarray(Re, float)
    return C["jd_pt"] / C["eps_b"] * Re ** C["jd_exp_pt"] * Re * Sc ** (1 / 3)
