"""Models of page A3.5 (Ranz and Marshall 1952).

Spherical boundary-value problems that give the correlation's constant 2.0 and
the radial-convection correction Ranz raised in 1993, plus the closed forms they
are checked against:

* :func:`nu_shell`, steady conduction (or diffusion) from a sphere into a
  stagnant shell, ``Nu = 2/(1 - a/b)``; linear, one sparse solve;
* :func:`sh_stefan`, diffusion out of the sphere through stagnant B with the
  radial convection that diffusion drives; :func:`sh_closed`, its closed form;
* :func:`pf_log`, :func:`pf_arith` and :func:`corr_excess`, the two readings of
  the pressure factor and the Sherwood excess they correct;
  :func:`stefan_threshold` and :func:`stefan_threshold_lambertw`, the surface
  mole fraction at which that excess reaches a given size, by a root and in
  closed form.

Units are those of the caller; only ratios matter. The sphere has unit radius.
"""

import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.special import lambertw
from pymrm import (compute_boundary_values, construct_convflux_upwind, construct_div,
                   construct_grad)
from pymrm.checks import find_roots


def _faces(n, ratio, geometric=True):
    """Faces from ``a = 1`` to ``b = ratio``; geometric faces nest when ``n`` doubles."""
    w = np.linspace(0.0, 1.0, n + 1)
    return float(ratio) ** w if geometric else 1.0 + (float(ratio) - 1.0) * w


def nu_shell(n, ratio, D=1.0, Ca=1.0, Cb=0.0, nu=2, geometric=True, outer="dirichlet",
             at="surface"):
    """Nusselt number of a unit sphere in a stagnant shell, by finite volumes.

    Solves ``(1/r^2) d/dr (r^2 D dC/dr) = 0`` (for ``nu = 2``) on ``a <= r <= b``.
    The problem is linear: one sparse solve.

    Parameters
    ----------
    n : int
        Number of cells in the radial direction.
    ratio : float
        Outer-to-inner radius ``b/a``.
    D : float, optional
        Diffusivity (or conductivity).
    Ca, Cb : float, optional
        Dirichlet values at the surface ``r = a`` and the far field ``r = b``.
    nu : int, optional
        Geometry exponent of ``construct_div``: 0 slab, 1 cylinder, 2 sphere.
    geometric : bool, optional
        Geometric faces (True) or uniform faces (False).
    outer : {"dirichlet", "noflux"}, optional
        Far field held at ``Cb``, or a closed shell (a defect).
    at : {"surface", "shell"}, optional
        Where the flux is read: at ``r = a``, or at ``r = b`` times the area
        ratio ``(b/a)^2``. For ``nu = 2`` the two agree to the residual of the
        solve, because the discrete divergence telescopes.

    Returns
    -------
    float
        ``2 a (-D dC/dr) / (D (Ca - Cb))``, the gradient read with
        ``compute_boundary_values`` under the solve's own boundary conditions.
        Exact for ``nu = 2``: ``2/(1 - a/b)``.
    """
    a, b = 1.0, float(ratio)
    r_f = _faces(n, b, geometric)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n, 1)
    # alpha dC/dn + beta C = d, n the outward normal:
    #   r = a : C = Ca (saturated drop surface)          -> a=0, b=1, d=Ca
    #   r = b : C = Cb (far field)                       -> a=0, b=1, d=Cb
    #           or dC/dn = 0 (closed shell, a defect)    -> a=1, b=0, d=0
    bc = ({"a": 0.0, "b": 1.0, "d": Ca},
          {"a": 0.0, "b": 1.0, "d": Cb} if outer == "dirichlet"
          else {"a": 1.0, "b": 0.0, "d": 0.0})
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc)
    div = construct_div(shape, r_f, nu=nu)
    lhs = (div @ (-D * grad)).tocsc()
    rhs = -np.asarray((div @ (-D * grad_bc)).todense()).ravel()
    C = spsolve(lhs, rhs).reshape(shape)
    if at == "surface":
        _, dcdr = compute_boundary_values(C, r_f, r_c, bc[0], bound_id=0)
        area = 1.0
    else:
        _, dcdr = compute_boundary_values(C, r_f, r_c, bc[1], bound_id=1)
        area = (b / a) ** 2
    # gradients come back along +r, so -D dC/dr is the flux out of the sphere
    return float(2.0 * a * area * (-D * np.ravel(dcdr)[0]) / (D * (Ca - Cb)))


def sh_stefan(n, ratio, ya, yb=0.0, D=1.0, c=1.0, convection=True, n_scan=24,
              balance=False):
    """Sherwood number of a unit sphere with radial (Stefan) convection.

    Species A diffuses out of the sphere through stagnant B (``N_B = 0``), so
    the molar-average velocity is ``v = N_A/c = W/(4 pi c r^2)`` with ``W`` the
    total molar flow. For a given ``W`` the balance
    ``(1/r^2) d/dr (r^2 (c v y - c D dy/dr)) = 0`` is linear in ``y``; ``W``
    is then the one root of "the flow the solve produces equals the flow it
    was given".

    The convection is first-order upwind: ``v > 0`` everywhere, so the upwind
    operator is linear in ``W`` and is assembled once for ``W = 1``.

    Parameters
    ----------
    n : int
        Number of cells (geometric radial grid).
    ratio : float
        Outer-to-inner radius ``b/a``.
    ya, yb : float
        Mole fraction of A at the surface and in the far field.
    D, c : float, optional
        Diffusivity and total molar concentration.
    convection : bool, optional
        False drops the radial convection from the solve (a defect).
    n_scan : int, optional
        Scan points of the root search.
    balance : bool, optional
        Also return the balance gap.

    Returns
    -------
    sh : float
        ``W / (2 pi a c D (ya - yb))``.
    W : float
        Total molar flow out of the sphere.
    gap : float
        Only with ``balance``: ``|4 pi b^2 N_A(b) - W| / W``, the flow through the
        outer face against the flow through the surface. The discrete divergence
        telescopes, so this is the residual of the solve and the root.

    Raises
    ------
    RuntimeError
        Unless exactly one root lies in ``[W_d/2, 2 W_d/(1 - ya)]``, with
        ``W_d = 4 pi a c D (ya - yb)/(1 - a/b)`` the flow of pure diffusion;
        the exact flow lies between ``W_d`` and ``W_d/(1 - ya)``.
    """
    a, b = 1.0, float(ratio)
    r_f = _faces(n, b)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n, 1)
    # alpha dy/dn + beta y = d, n the outward normal; both ends fixed mole fractions:
    #   r = a : y = ya (saturated surface)   -> a=0, b=1, d=ya
    #   r = b : y = yb (far field)           -> a=0, b=1, d=yb
    bc = ({"a": 0.0, "b": 1.0, "d": ya}, {"a": 0.0, "b": 1.0, "d": yb})
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc)
    div = construct_div(shape, r_f, nu=2)
    v1 = (1.0 / (4.0 * np.pi * c)) / r_f ** 2            # face velocity per unit W
    conv1, conv1_bc = construct_convflux_upwind(shape, r_f, r_c, bc, v=v1.reshape(-1, 1))
    on = 1.0 if convection else 0.0
    dif, dif_bc = div @ (-c * D * grad), div @ (-c * D * grad_bc)
    cv, cv_bc = div @ (c * conv1), div @ (c * conv1_bc)

    def produced(W, end=0):
        lhs = (dif + on * W * cv).tocsc()
        rhs = -np.asarray((dif_bc + on * W * cv_bc).todense()).ravel()
        y = spsolve(lhs, rhs).reshape(shape)
        y_f, dydr = compute_boundary_values(y, r_f, r_c, bc[end], bound_id=end)
        # radial flux of A at the face: c v y - c D dy/dr (gradient along +r)
        N = (on * W * v1[-end] * c * float(np.ravel(y_f)[0])
             - c * D * float(np.ravel(dydr)[0]))
        return 4.0 * np.pi * r_f[-end] ** 2 * N

    W_d = 4.0 * np.pi * a * c * D * (ya - yb) / (1.0 - a / b)
    found = find_roots(lambda W: produced(W) - W, 0.5 * W_d, 2.0 * W_d / (1.0 - ya),
                       n_scan=n_scan, xtol=1e-14 * W_d)
    if len(found["roots"]) != 1:
        raise RuntimeError(f"sh_stefan: {len(found['roots'])} roots for W ({found['message']})")
    W = float(found["roots"][0])
    sh = float(W / (2.0 * np.pi * a * c * D * (ya - yb)))
    if balance:
        return sh, W, abs(produced(W, end=1) - W) / W
    return sh, W


def sh_closed(ratio, ya, yb=0.0):
    """Exact Sherwood number for the spherical stagnant film.

    Parameters
    ----------
    ratio : float
        Outer-to-inner radius ``b/a``.
    ya, yb : float
        Mole fraction of A at the surface and in the far field.

    Returns
    -------
    float
        ``2 ln((1 - yb)/(1 - ya)) / ((1 - a/b) (ya - yb))``.
    """
    return 2.0 * np.log((1 - yb) / (1 - ya)) / ((1 - 1 / ratio) * (ya - yb))


def pf_log(y):
    """Pressure factor ``p_f/pi`` if ``p_f`` is the logarithmic mean.

    Parameters
    ----------
    y : float or ndarray
        Surface mole fraction of A.

    Returns
    -------
    float or ndarray
        ``y / ln(1/(1 - y))``; removes the radial-convection enhancement exactly.
    """
    return y / np.log(1.0 / (1.0 - y))        # /pi


def pf_arith(y):
    """Pressure factor ``p_f/pi`` if ``p_f`` is the arithmetic mean.

    Parameters
    ----------
    y : float or ndarray
        Surface mole fraction of A.

    Returns
    -------
    float or ndarray
        ``1 - y/2``.
    """
    return 1.0 - 0.5 * y                       # /pi


def corr_excess(y):
    """Relative Sherwood enhancement from radial convection, ``Sh/2 - 1``.

    Parameters
    ----------
    y : float or ndarray
        Surface mole fraction of A, 0 < y < 1.

    Returns
    -------
    float or ndarray
        ``-ln(1 - y)/y - 1``.
    """
    return -np.log1p(-y) / y - 1.0       # Sh/2 - 1


def stefan_threshold(excess, lo=1e-6, hi=0.999):
    """Surface mole fraction at which radial convection raises ``Sh`` by ``excess``.

    Parameters
    ----------
    excess : float
        Relative enhancement, ``Sh/2 - 1``.
    lo, hi : float, optional
        Range of the root search; ``corr_excess`` rises from 0 at ``y = 0`` to
        infinity at ``y = 1``.

    Returns
    -------
    float
        The one root of ``corr_excess(y) = excess`` on ``[lo, hi]``.

    Raises
    ------
    RuntimeError
        Unless exactly one root is found.
    """
    found = find_roots(lambda y: corr_excess(y) - excess, lo, hi, n_scan=64, xtol=1e-15)
    if len(found["roots"]) != 1:
        raise RuntimeError(f"stefan_threshold: {len(found['roots'])} roots ({found['message']})")
    return float(found["roots"][0])


def stefan_threshold_lambertw(excess, branch=0):
    """The same surface mole fraction in closed form.

    With ``s = 1 + excess`` and ``u = 1 - y``, ``-ln(1 - y) = s y`` reads
    ``-s u exp(-s u) = -s exp(-s)``, whose nontrivial solution is the principal
    branch: ``y = 1 + W_0(-s exp(-s))/s``.

    Parameters
    ----------
    excess : float
        Relative enhancement, ``Sh/2 - 1``.
    branch : int, optional
        Branch of Lambert's W; ``-1`` returns the trivial solution ``y = 0``
        (a defect).

    Returns
    -------
    float
        ``1 + W_0(-s e^-s)/s``.
    """
    s = 1.0 + excess
    return float(1.0 + lambertw(-s * np.exp(-s), branch).real / s)
