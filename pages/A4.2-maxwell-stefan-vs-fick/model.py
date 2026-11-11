"""Maxwell-Stefan and Fick fluxes for a ternary gas in a two-bulb cell (page A4.2).

Three flux closures (Maxwell-Stefan via [B], generalized Fick [D](x), scalar
Wilke effective diffusivity), an independent route through the raw friction
system, and the quasi-steady capillary plus bulb march (Duncan-Toor cell).

The page's parameters and pymrm operators are module state, handed over once by
`configure`: c_t, d_ms, n_c, n_i, area, vol_1, vol_2, t_end, and the capillary
operators from `build_operators` (n_z, shape, grad_mat, grad_bc_1, grad_bc_2,
bv_1, bv_bc_1, bv_2, bv_bc_2, div_mat, numjac, face). `operators` swaps the
operators for one block of code (grid ladder, defect rows) and restores them.
"""
from contextlib import contextmanager

import numpy as np
from pymrm import (NumJac, construct_boundary_value_matrices, construct_div,
                   construct_grad, newton)

c_t = d_ms = n_c = n_i = area = vol_1 = vol_2 = t_end = None
n_z = shape = grad_mat = grad_bc_1 = grad_bc_2 = div_mat = numjac = None
bv_1 = bv_bc_1 = bv_2 = bv_bc_2 = face = None
OPERATOR_KEYS = ("n_z", "shape", "grad_mat", "grad_bc_1", "grad_bc_2", "bv_1",
                 "bv_bc_1", "bv_2", "bv_bc_2", "div_mat", "numjac", "face")


def configure(**params):
    """Set the module state that the flux functions read.

    Parameters
    ----------
    **params
        Any of c_t (mol/m3), d_ms ((3, 3) pair diffusivities, m2/s), n_c, n_i,
        area (m2), vol_1, vol_2 (m3), t_end (s), and the keys returned by
        `build_operators`.
    """
    globals().update(params)


def build_operators(n_cells, length, n_ind=2, nu=0, face_rule="mean"):
    """Capillary grid and constant pymrm operators, assembled once.

    Both ends carry the bulb composition as a Dirichlet condition. Outward
    normal, a dx/dn + b x = d with d = 1 and the bulb composition passed as the
    external vector of `shapes_d`:
    z = 0: x_i = x_i,bulb1  ->  a = 0, b = 1, d = 1 (times x_bulb1);
    z = L: x_i = x_i,bulb2  ->  a = 0, b = 1, d = 1 (times x_bulb2).

    Parameters
    ----------
    n_cells : int
        Number of capillary cells.
    length : float
        Capillary length, m.
    n_ind : int, optional
        Independent mole fractions (fields per cell).
    nu : int, optional
        Geometry of `construct_div`: 0 for the straight capillary.
    face_rule : {"mean", "left"}, optional
        Interior face composition: mean of the two neighbouring cells, or the
        left cell alone (a first-order rule, used as a defect).

    Returns
    -------
    dict
        Keys of OPERATOR_KEYS, ready for `configure`.
    """
    shp = (n_cells, n_ind)
    z_f = np.linspace(0.0, length, n_cells + 1)
    z_c = 0.5 * (z_f[:-1] + z_f[1:])
    bc = ({"a": 0.0, "b": 1.0, "d": 1.0},          # z = 0: x = x_bulb1
          {"a": 0.0, "b": 1.0, "d": 1.0})          # z = L: x = x_bulb2
    g, g1, g2 = construct_grad(shp, z_f, z_c, bc, axis=0,
                               shapes_d=((1, n_ind), (1, n_ind)))
    b1, b1_bc = construct_boundary_value_matrices(shp, z_f, z_c, bc[0], axis=0,
                                                  bound_id=0, shape_d=(1, n_ind))
    b2, b2_bc = construct_boundary_value_matrices(shp, z_f, z_c, bc[1], axis=0,
                                                  bound_id=1, shape_d=(1, n_ind))
    return dict(n_z=n_cells, shape=shp, grad_mat=g, grad_bc_1=g1, grad_bc_2=g2,
                bv_1=b1, bv_bc_1=b1_bc, bv_2=b2, bv_bc_2=b2_bc,
                div_mat=construct_div(shp, z_f, nu=nu, axis=0),
                numjac=NumJac(shp, axes_diagonals=[0], axes_blocks=[-1]),
                face=face_rule)


@contextmanager
def operators(n_cells, length, **kw):
    """Use the operators of `build_operators(n_cells, length, **kw)` inside a
    `with` block, and restore the configured ones afterwards."""
    old = {k: globals()[k] for k in OPERATOR_KEYS}
    configure(**build_operators(n_cells, length, n_ind=n_i, **kw))
    try:
        yield
    finally:
        configure(**old)



def full_x(x_ind):
    """Append the dependent mole fraction: x_n = 1 - sum(x_1..x_{n-1}).

    Parameters
    ----------
    x_ind : numpy.ndarray
        Independent mole fractions, shape (..., n_c - 1).

    Returns
    -------
    numpy.ndarray
        All mole fractions, shape (..., n_c).
    """
    return np.concatenate([x_ind, 1.0 - x_ind.sum(axis=-1, keepdims=True)], axis=-1)


def build_b(x_ind, dms=None):
    """[B] of review eq. (26), shape (..., 2, 2). Units s/m2.

    Species n_c-1 is the one eliminated; `dms` is passed explicitly so the
    same code can be run with the species reordered (see the permutation
    table on the page).

    Parameters
    ----------
    x_ind : numpy.ndarray
        Independent mole fractions, shape (..., n_i).
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        [B], shape (..., n_i, n_i), s/m2.
    """
    dms = d_ms if dms is None else dms
    x = full_x(x_ind)
    b = np.zeros(x_ind.shape[:-1] + (n_i, n_i))
    for i in range(n_i):
        acc = x[..., i] / dms[i, n_c - 1]
        for k in range(n_c):
            if k != i:
                acc = acc + x[..., k] / dms[i, k]
        b[..., i, i] = acc
        for j in range(n_i):
            if j != i:
                b[..., i, j] = -x[..., i] * (1.0 / dms[i, j] - 1.0 / dms[i, n_c - 1])
    return b


def fick_matrix(x_ind, dms=None):
    """[D](x) = [B]^-1 for the ideal gas, as explicit closed-form elements.

    Writing the four composition-dependent Fick coefficients out (adjugate over
    determinant) is the point: this is the 'generalized Fick' description as a
    modeller would code it, with no linear solve anywhere.

    Parameters
    ----------
    x_ind : numpy.ndarray
        Independent mole fractions, shape (..., n_i).
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        [D], shape (..., n_i, n_i), m2/s.
    """
    dms = d_ms if dms is None else dms
    b = build_b(x_ind, dms)
    det = b[..., 0, 0] * b[..., 1, 1] - b[..., 0, 1] * b[..., 1, 0]
    d = np.empty_like(b)
    d[..., 0, 0] =  b[..., 1, 1] / det
    d[..., 0, 1] = -b[..., 0, 1] / det
    d[..., 1, 0] = -b[..., 1, 0] / det
    d[..., 1, 1] =  b[..., 0, 0] / det
    return d


def friction_solve(x_full, grad_ind, dms=None):
    """The independent route: eqs. (16)-(17) solved as they stand.

    Assembles the n friction balances -grad(x_i) = sum_j (x_j J_i - x_i J_j)
    / (c_t D_ij) directly, replaces the last (linearly dependent) one by the
    closed-cell constraint sum_k J_k = 0, and returns all n fluxes. Nothing is
    eliminated, [B] is never formed and `build_b` is never called, so this
    path can disagree with `flux_ms`/`flux_fick` if eq. (26) or the [B]^-1
    transformation is coded wrongly.

    Parameters
    ----------
    x_full : numpy.ndarray
        All mole fractions at the faces, shape (n_faces, n_c).
    grad_ind : numpy.ndarray
        Gradients of the independent species, shape (n_faces, n_i), 1/m.
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        Fluxes of all species, shape (n_faces, n_c), mol/(m2 s).
    """
    dms = d_ms if dms is None else dms
    n_f = x_full.shape[0]
    a = np.zeros((n_f, n_c, n_c))
    for i in range(n_c - 1):
        for j in range(n_c):
            if j != i:
                a[:, i, i] += x_full[:, j] / dms[i, j]
                a[:, i, j] -= x_full[:, i] / dms[i, j]
    a[:, n_c - 1, :] = 1.0                        # bootstrap: sum_k J_k = 0
    rhs = np.zeros((n_f, n_c))
    rhs[:, :n_c - 1] = -c_t * grad_ind            # the n-1 independent forces
    return np.linalg.solve(a, rhs[..., None])[..., 0]


def friction_residual(x_full, grad_ind, j_full, dms=None):
    """Residual of the friction balance that `friction_solve` did not use.

    It must be satisfied automatically, because sum_i grad(x_i) = 0 makes the
    n balances linearly dependent. Scaled by the largest driving force.

    Parameters
    ----------
    x_full : numpy.ndarray
        All mole fractions at the faces, shape (n_faces, n_c).
    grad_ind : numpy.ndarray
        Gradients of the independent species, shape (n_faces, n_i), 1/m.
    j_full : numpy.ndarray
        Fluxes of all species from `friction_solve`, shape (n_faces, n_c).
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    float
        Largest scaled residual of the unused balance.
    """
    dms = d_ms if dms is None else dms
    i = n_c - 1
    g_i = -grad_ind.sum(axis=-1)
    lhs = sum((x_full[:, j] * j_full[:, i] - x_full[:, i] * j_full[:, j])
              / (c_t * dms[i, j]) for j in range(n_c) if j != i)
    return np.abs(lhs + g_i).max() / np.abs(grad_ind).max()


def fick_matrix_from_friction(x_ind, dms=None):
    """[D] recovered column by column from `friction_solve`.

    Apply a unit gradient in one independent species (species n taking -1 by
    closure) and read the fluxes. Never touches build_b.

    Parameters
    ----------
    x_ind : numpy.ndarray
        Independent mole fractions, shape (..., n_i).
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        [D], shape (..., n_i, n_i), m2/s.
    """
    dms = d_ms if dms is None else dms
    x_full = full_x(x_ind).reshape(-1, n_c)
    cols = []
    for j in range(n_i):
        g = np.zeros((x_full.shape[0], n_i))
        g[:, j] = 1.0
        cols.append(-friction_solve(x_full, g, dms)[:, :n_i] / c_t)
    return np.stack(cols, axis=-1).reshape(x_ind.shape[:-1] + (n_i, n_i))


def wilke_d(x_ind, dms=None, which=None):
    """Wilke effective diffusivity, review eq. (32).

    Parameters
    ----------
    x_ind : numpy.ndarray
        Independent mole fractions, shape (..., n_i).
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.
    which : iterable of int, optional
        Species to evaluate; the default is the n-1 species the scalar law is
        imposed on, i.e. all but the eliminated one.

    Returns
    -------
    numpy.ndarray
        Effective diffusivities, shape (..., len(which)), m2/s.
    """
    dms = d_ms if dms is None else dms
    x = full_x(x_ind)
    which = range(n_i) if which is None else which
    return np.stack([(1.0 - x[..., i]) /
                     sum(x[..., k] / dms[i, k] for k in range(n_c) if k != i)
                     for i in which], axis=-1)


def grad_of(x, xb1, xb2):
    """Mole-fraction gradients at all faces (constant operators + bc terms).

    Parameters
    ----------
    x : numpy.ndarray
        Cell mole fractions, shape (n_z, n_i).
    xb1, xb2 : numpy.ndarray
        Bulb compositions at z = 0 and z = L (first n_i entries used).

    Returns
    -------
    numpy.ndarray
        Gradients, shape (n_z + 1, n_i), 1/m.
    """
    return (grad_mat @ x.reshape((-1, 1))
            + grad_bc_1 @ xb1[:n_i].reshape((-1, 1))
            + grad_bc_2 @ xb2[:n_i].reshape((-1, 1))).reshape(n_z + 1, n_i)


def face_of(x, xb1, xb2):
    """Mole fractions at all faces.

    Interior faces: the mean of the two neighbouring cells (uniform grid,
    second order), or the left cell when the operators were built with
    face_rule="left". End faces: the boundary values of the Dirichlet
    conditions, from the boundary-value operators (they return the bulb
    compositions).

    Parameters
    ----------
    x : numpy.ndarray
        Cell mole fractions, shape (n_z, n_i).
    xb1, xb2 : numpy.ndarray
        Bulb compositions at z = 0 and z = L (first n_i entries used).

    Returns
    -------
    numpy.ndarray
        Face mole fractions, shape (n_z + 1, n_i).
    """
    xf = np.empty((n_z + 1, n_i))
    xf[1:-1] = 0.5 * (x[:-1] + x[1:]) if face == "mean" else x[:-1]
    xc = x.reshape((-1, 1))
    xf[0] = (bv_1 @ xc + bv_bc_1 @ xb1[:n_i].reshape((-1, 1))).ravel()
    xf[-1] = (bv_2 @ xc + bv_bc_2 @ xb2[:n_i].reshape((-1, 1))).ravel()
    return xf


def flux_ms(x, xb1, xb2, dms=None):
    """Maxwell-Stefan: solve [B] J = -c_t grad(x) at every face.

    Parameters
    ----------
    x : numpy.ndarray
        Cell mole fractions, shape (n_z, n_i).
    xb1, xb2 : numpy.ndarray
        Bulb compositions at z = 0 and z = L.
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        Face fluxes, shape (n_z + 1, n_i), mol/(m2 s).
    """
    dms = d_ms if dms is None else dms
    g, xf = grad_of(x, xb1, xb2), face_of(x, xb1, xb2)
    return -c_t * np.linalg.solve(build_b(xf, dms), g[..., None])[..., 0]


def flux_fick(x, xb1, xb2, dms=None):
    """Generalized Fick: J = -c_t [D](x) grad(x), [D] explicit.

    Parameters
    ----------
    x : numpy.ndarray
        Cell mole fractions, shape (n_z, n_i).
    xb1, xb2 : numpy.ndarray
        Bulb compositions at z = 0 and z = L.
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        Face fluxes, shape (n_z + 1, n_i), mol/(m2 s).
    """
    dms = d_ms if dms is None else dms
    g, xf = grad_of(x, xb1, xb2), face_of(x, xb1, xb2)
    return -c_t * np.einsum('fij,fj->fi', fick_matrix(xf, dms), g)


def flux_ms_friction(x, xb1, xb2, dms=None):
    """Maxwell-Stefan from the full friction system, no [B], nothing eliminated.

    Parameters
    ----------
    x : numpy.ndarray
        Cell mole fractions, shape (n_z, n_i).
    xb1, xb2 : numpy.ndarray
        Bulb compositions at z = 0 and z = L.
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        Face fluxes, shape (n_z + 1, n_i), mol/(m2 s).
    """
    dms = d_ms if dms is None else dms
    g, xf = grad_of(x, xb1, xb2), face_of(x, xb1, xb2)
    return friction_solve(full_x(xf), g, dms)[:, :n_i]


def flux_wilke(x, xb1, xb2, dms=None):
    """Scalar effective diffusivity: J_i = -c_t D_i^W(x) grad(x_i).

    Parameters
    ----------
    x : numpy.ndarray
        Cell mole fractions, shape (n_z, n_i).
    xb1, xb2 : numpy.ndarray
        Bulb compositions at z = 0 and z = L.
    dms : numpy.ndarray, optional
        Pair diffusivities (n_c, n_c), m2/s; default the configured d_ms.

    Returns
    -------
    numpy.ndarray
        Face fluxes, shape (n_z + 1, n_i), mol/(m2 s).
    """
    dms = d_ms if dms is None else dms
    g, xf = grad_of(x, xb1, xb2), face_of(x, xb1, xb2)
    return -c_t * wilke_d(xf, dms) * g


def solve_capillary(flux, x_init, xb1, xb2, tol=1e-12):
    """Quasi-steady capillary: div(J_i) = 0 with bulb-composition Dirichlet bc.

    The unknowns are mole fractions of order one, so the absolute Newton step
    tolerance is a relative one; 1e-12 sits above the round-off floor of the
    flux evaluation and `stagnation=True` stops a solve that reaches that
    floor first.

    Parameters
    ----------
    flux : callable
        Flux closure `flux(x, xb1, xb2)` returning face fluxes.
    x_init : numpy.ndarray
        Initial guess, shape (n_z, n_i).
    xb1, xb2 : numpy.ndarray
        Bulb compositions at z = 0 and z = L.
    tol : float, optional
        Newton step tolerance (mole fraction).

    Returns
    -------
    numpy.ndarray
        Converged cell mole fractions, shape (n_z, n_i).

    Raises
    ------
    RuntimeError
        If Newton does not converge.
    """
    def residual(x):
        g, jac = numjac(
            lambda xx: (div_mat @ flux(xx, xb1, xb2).reshape((-1, 1))).reshape(shape), x)
        return g.reshape((-1, 1)), jac
    result = newton(residual, x_init, maxfev=30, tol=tol, stagnation=True)
    if not result.success:
        raise RuntimeError(f"capillary solve failed: {result.message}")
    return result.x.reshape(shape)


def simulate(flux, x1_init, x2_init, n_steps=400, t_end=None, theta=1.0, tol=1e-12):
    """March the two bulbs in time; the capillary is quasi-steady.

    Bulb balances V_1 c_t dx_1/dt = -A J(0), V_2 c_t dx_2/dt = +A J(L), stepped
    with the theta rule: theta = 1 is backward Euler (the page's scheme),
    theta = 0.5 the trapezoidal rule. The new bulb compositions and the
    capillary profile are coupled by three fixed-point iterations per step.

    Parameters
    ----------
    flux : callable
        Flux closure `flux(x, xb1, xb2)`.
    x1_init, x2_init : numpy.ndarray
        Initial bulb compositions (all n_c species).
    n_steps : int, optional
        Number of time steps.
    t_end : float, optional
        End time in s; default the configured t_end.
    theta : float, optional
        Weight of the end-of-step flux in the bulb update.
    tol : float, optional
        Newton step tolerance of every capillary solve.

    Returns
    -------
    tuple
        (t_h, h1, h2, hJ): time in h, bulb histories (n_steps + 1, n_c), and
        the quasi-steady capillary flux (face mean) at the same times,
        (n_steps + 1, n_i).
    """
    t_end = globals()["t_end"] if t_end is None else t_end
    dt = t_end / n_steps
    xb1, xb2 = x1_init.copy(), x2_init.copy()
    x = np.linspace(0, 1, n_z)[:, None] * (xb2[:n_i] - xb1[:n_i]) + xb1[:n_i]
    x = solve_capillary(flux, x, xb1, xb2, tol)
    j_old = flux(x, xb1, xb2)
    t_h, h1, h2, hJ = [0.0], [xb1.copy()], [xb2.copy()], [j_old.mean(axis=0)]

    for _ in range(n_steps):
        xb1n, xb2n = xb1.copy(), xb2.copy()
        for _ in range(3):                     # fixed-point iterations on the coupling
            x = solve_capillary(flux, x, xb1n, xb2n, tol)
            j_face = flux(x, xb1n, xb2n)
            j_0 = theta * j_face[0] + (1.0 - theta) * j_old[0]       # face z = 0
            j_l = theta * j_face[-1] + (1.0 - theta) * j_old[-1]     # face z = L
            xb1n[:n_i] = xb1[:n_i] - area * j_0 / (vol_1 * c_t) * dt
            xb2n[:n_i] = xb2[:n_i] + area * j_l / (vol_2 * c_t) * dt
            xb1n[n_i] = 1.0 - xb1n[:n_i].sum()
            xb2n[n_i] = 1.0 - xb2n[:n_i].sum()
        xb1, xb2, j_old = xb1n, xb2n, j_face
        t_h.append(t_h[-1] + dt); h1.append(xb1.copy()); h2.append(xb2.copy())
        hJ.append(j_face.mean(axis=0))
    return (np.asarray(t_h) / 3600.0, np.asarray(h1), np.asarray(h2), np.asarray(hJ))
