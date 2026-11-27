"""Sit and Grace (1981): model functions for page E1.4.

The single-bubble interphase coefficient of eqns (1), (2) as a throughflow term plus a
penetration term, its enhanced forms (eqns 3, 4 and the page's alternatives), the
potential-flow surface-renewal mean of the penetration term, and two small pymrm solves:
the scaled penetration problem behind the diffusive term (finite volumes, Dirichlet face,
exact in time through the matrix exponential of the augmented constant operator, or
backward Euler) and the tracer balance of a growing two-dimensional bubble (an implicit
march in time). Every quantity also has a closed form here. cgs units throughout
(cm, s). No plotting, no data loading.
"""

import math

import numpy as np
from scipy import integrate, linalg, sparse
from scipy.sparse.linalg import expm_multiply, splu
from pymrm import compute_boundary_values, construct_div, construct_grad

# Throughflow q/S_b in units of U_mf, by dimension of the bubble. Murray's is the one
# eqns (1), (2) use (U_mf/pi, U_mf/4); Davidson's (q = 4 U_mf R per unit depth in 2-D,
# 3 pi R^2 U_mf in 3-D) is the one page E1.2 solves for.
THROUGHFLOW = {"murray": {2: 1.0 / math.pi, 3: 0.25},
               "davidson": {2: 2.0 / math.pi, 3: 0.75}}


def k_throughflow(umf, dim, rule="murray", table=THROUGHFLOW):
    """Throughflow term of the coefficient, q/S_b.

    Parameters
    ----------
    umf : float or ndarray
        Minimum fluidization velocity (cm/s).
    dim : int
        2 for a circular two-dimensional bubble, 3 for a sphere.
    rule : str
        ``"murray"`` (eqns 1, 2) or ``"davidson"`` (page E1.2).
    table : dict
        Coefficient table; a break row passes a changed one.

    Returns
    -------
    float or ndarray
        k_T in cm/s.
    """
    return table[rule][dim] * np.asarray(umf, dtype=float)


# codehealth: allow params -- the printed formula's four inputs and the two readings a break row changes
def k_diffusion(DG, eps, UB, DB, convention="capacity", size_factor=1.0):
    """Penetration term of eqns (1) to (4), (4 D_G eps_mf U_B/(pi D_B))^(1/2).

    Parameters
    ----------
    DG, eps, UB, DB : float or ndarray
        Molecular diffusivity (cm2/s), voidage at minimum fluidization, bubble velocity
        (cm/s) and bubble diameter (cm).
    convention : str
        ``"capacity"``: capacity eps_mf with diffusivity D_G, the printed form.
        ``"effective"``: effective diffusivity eps_mf D_G as well, which puts eps_mf
        outside the root.
    size_factor : float
        Multiplies D_B (1 as printed; 0.5 puts the radius in its place).

    Returns
    -------
    float or ndarray
        k_D in cm/s.
    """
    DG, eps, UB, DB = (np.asarray(a, dtype=float) for a in (DG, eps, UB, DB))
    D_eff = DG * eps if convention == "effective" else DG
    return np.sqrt(4.0 * D_eff * eps * UB / (math.pi * DB * size_factor))


def bubble_velocity(DB, c, g=981.0):
    """Rise velocity c (g D_B)^(1/2) (cm/s), the page's bracket for unprinted U_B."""
    return c * np.sqrt(g * np.asarray(DB, dtype=float))


def renewal_factor(geometry, us_coeff=None, epsrel=1e-13):
    """Potential-flow surface-renewal mean coefficient over (4 D U/(pi d))^(1/2).

    The emulsion moves along the bubble surface at u_s = a U sin(theta) (a = 2 for a
    circle, 1.5 for a sphere); an element arriving at theta has been exposed for the
    time it took to travel from the nose, and its local flux is
    C (D/pi)^(1/2) u_s r^j / (int_0^theta u_s r^(2j) R dtheta')^(1/2), with r = R
    sin(theta) and j = 0 (circle) or 1 (sphere). Both integrals by adaptive
    quadrature, no closed form used.

    Parameters
    ----------
    geometry : str
        ``"circle"`` or ``"sphere"``.
    us_coeff : float, optional
        Surface-velocity coefficient a; default 2 (circle) or 1.5 (sphere).
    epsrel : float
        Relative tolerance of both quadratures.

    Returns
    -------
    float
        Surface mean of the local coefficient divided by (4 D U/(pi d))^(1/2).
    """
    j = 0 if geometry == "circle" else 1
    a = us_coeff if us_coeff is not None else (2.0 if j == 0 else 1.5)
    R, U, D = 1.0, 1.0, 1.0
    us = lambda th: a * U * math.sin(th)
    r = lambda th: R * math.sin(th)

    def inner(th):
        return integrate.quad(lambda s: us(s) * r(s) ** (2 * j) * R, 0.0, th,
                              epsabs=0.0, epsrel=epsrel, limit=200)[0]

    local = lambda th: math.sqrt(D / math.pi) * us(th) * r(th) ** j / math.sqrt(inner(th))
    weight = lambda th: r(th) ** j * R            # surface element per dtheta
    num = integrate.quad(lambda th: local(th) * weight(th), 0.0, math.pi,
                         epsabs=0.0, epsrel=epsrel, limit=200)[0]
    den = integrate.quad(weight, 0.0, math.pi, epsabs=0.0, epsrel=epsrel, limit=200)[0]
    return (num / den) / math.sqrt(4.0 * D * U / (math.pi * 2.0 * R))


class Penetration:
    """The scaled penetration problem c_t = c_xx on 0 < x < L with pymrm.

    Dirichlet faces, c(0, t) = 1 and c(L, t) = 0, and c(x, 0) = 0. Finite volumes on a
    uniform grid: ``construct_grad`` with the two Dirichlet dictionaries and
    ``construct_div`` (Cartesian, nu = 0), assembled once, give dc/dt = A c + b.

    Parameters
    ----------
    n : int
        Number of cells.
    L : float
        Domain length (12 puts the far face where erfc(L/2) = 2e-17).
    d_left : float
        Datum of the left Dirichlet dictionary (1; a mis-signed dictionary gives -1).
    first_cell : float
        Initial value of the first cell (0; a defect puts the boundary value there).
    """

    def __init__(self, n, L=12.0, d_left=1.0, first_cell=0.0):
        self.n, self.L = n, L
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[1:] + self.x_f[:-1])
        self.vol = np.diff(self.x_f)                      # exact cell volumes
        # outward normal (P1): a dc/dn + b c = d. Left: c(0) = d_left; right: c(L) = 0.
        self.bc = ({"a": 0.0, "b": 1.0, "d": d_left}, {"a": 0.0, "b": 1.0, "d": 0.0})
        self.grad, grad_bc = construct_grad((n, 1), self.x_f, self.x_c, self.bc)
        self.div = construct_div((n, 1), self.x_f, nu=0)   # Cartesian slab
        self.grad_bc = grad_bc.toarray().ravel()
        self.A = self.div @ self.grad
        self.b = self.div @ self.grad_bc
        self.c0 = np.zeros(n)
        self.c0[0] = first_cell

    def augmented(self, bc_scale=1.0):
        """Dense constant operator on z = (c, Q_in, Q_out, 1).

        dQ_in/dt is the flux entering at x = 0 and dQ_out/dt the flux leaving at
        x = L, both from the operator's own face gradients, so the cell amount equals
        Q_in - Q_out up to round-off (discrete conservation). ``bc_scale`` scales the
        boundary column (1; the round-off companion's break row perturbs it)."""
        n = self.n
        G = self.grad.toarray()
        M = np.zeros((n + 3, n + 3))
        M[:n, :n] = self.A.toarray()
        M[:n, -1] = self.b * bc_scale
        M[n, :n], M[n, -1] = -G[0], -self.grad_bc[0] * bc_scale
        M[n + 1, :n], M[n + 1, -1] = -G[-1], -self.grad_bc[-1] * bc_scale
        return M

    def exact_time(self, t=1.0, method="expm", bc_scale=1.0):
        """State at time t, exact in time: expm of the augmented operator.

        Returns
        -------
        dict
            ``c`` (cells), ``amount`` (sum of exact cell volumes times c), ``q_in``,
            ``q_out`` (time-integrated face fluxes) and ``flux`` (the face flux into
            the domain at x = 0 at time t, read with ``compute_boundary_values``).
        """
        M = self.augmented(bc_scale) * t
        z0 = np.concatenate([self.c0, [0.0, 0.0, 1.0]])
        if method == "expm":
            z = linalg.expm(M) @ z0
        else:
            z = expm_multiply(sparse.csc_array(M), z0)
        c = z[:self.n]
        return {"c": c, "amount": float(self.vol @ c), "q_in": float(z[self.n]),
                "q_out": float(z[self.n + 1]), "flux": self.face_flux(c)}

    def face_flux(self, c):
        """Flux -dc/dx into the domain at x = 0, from the face gradient under self.bc."""
        _, grad_left, _, _ = compute_boundary_values(c.reshape(self.n, 1), self.x_f,
                                                     self.x_c, self.bc)
        return -float(np.ravel(grad_left)[0])   # gradient along +x (P1)

    # codehealth: allow params -- time, scheme weight and the two solver switches of the ladder
    def stepped(self, nt, t=1.0, theta=1.0, route="splu", bc_scale=1.0):
        """Amount at time t after nt steps of the theta scheme (1: backward Euler).

        One sparse LU per call (the operator is constant). ``route="banded"`` solves the same
        tridiagonal system with LAPACK's banded solver instead, the round-off companion."""
        dt = t / nt
        I = sparse.eye_array(self.n, format="csc")
        L_ = I - theta * dt * self.A
        if route == "banded":
            ab = np.zeros((3, self.n))
            for k in (-1, 0, 1):
                d = L_.diagonal(k)
                ab[1 - k, max(k, 0):self.n + min(k, 0)] = d
            assert np.isclose(abs(L_).sum(), np.abs(ab).sum(), rtol=1e-12), "the operator is not tridiagonal"
            lhs = lambda r: linalg.solve_banded((1, 1), ab, r)
        else:
            lu = splu(L_, permc_spec="COLAMD")
            lhs = lu.solve
        rhs_op = I + (1.0 - theta) * dt * self.A
        src = dt * self.b * bc_scale
        c = self.c0.copy()
        for _ in range(nt):
            c = lhs(rhs_op @ c + src)
        return float(self.vol @ c)


# ----------------------------------------------------------- tracer balance of a 2-D bubble

def equivalent_diameter(V, w):
    """Diameter of the circle with the frontal area V/w of a bubble of volume V (cm)."""
    return np.sqrt(4.0 * np.asarray(V, dtype=float) / (math.pi * w))


# codehealth: allow params -- the interval, the two volumes, the thickness and the perimeter factor
def int_s_over_v(t0, t1, V0, V1, w, phi=1.0):
    """Closed form of int S/V dt for V linear in t, S = phi w pi D_eq (s/cm).

    S/V = 4 phi/D_eq; with V0 = V1 the volume is constant."""
    if V1 == V0:
        return 4.0 * phi / float(equivalent_diameter(V0, w)) * (t1 - t0)
    return 4.0 * phi * math.sqrt(math.pi * w) * (math.sqrt(V1) - math.sqrt(V0)) * (t1 - t0) / (V1 - V0)


# codehealth: allow params -- the closed form's interval, volumes, thickness and the eqn (1) inputs
def int_kd_s_over_v(t0, t1, V0, V1, w, DG, eps, UB, phi=1.0, convention="capacity",
                    size_factor=1.0):
    """Closed form of int k_D(D_eq) S/V dt, the penetration part of eqn (1) along a
    bubble growing linearly in volume (dimensionless)."""
    pref = float(k_diffusion(DG, eps, UB, 1.0, convention, size_factor))   # k_D at D = 1
    cc = math.sqrt(4.0 / (math.pi * w))
    if V1 == V0:
        return pref * 4.0 * phi * float(equivalent_diameter(V0, w)) ** -1.5 * (t1 - t0)
    return pref * 4.0 * phi * cc ** -1.5 * 4.0 * (V1 ** 0.25 - V0 ** 0.25) * (t1 - t0) / (V1 - V0)


# codehealth: allow params -- the march's interval, volumes, thickness and the scheme switches
def tracer_march(K_of_D, t0, t1, V0, V1, y0, nsteps, w, phi=1.0, theta=1.0, dilution=True,
                 accumulate="product"):
    """Implicit march of d(V C)/dt = -K S C for y = V C, V linear in time.

    Parameters
    ----------
    K_of_D : callable
        Coefficient K (cm/s) as a function of the equivalent diameter D_eq (cm).
    t0, t1, V0, V1, y0 : float
        Interval, volumes at its ends (cm3) and V C at t0 (cm3).
    nsteps : int
        Number of steps.
    w : float
        Bubble thickness (cm).
    phi : float
        Perimeter factor on S.
    theta : float
        1 for backward Euler, 0.5 for the trapezoidal rule.
    dilution : bool
        True marches V C; False marches C alone (the growth dilution dropped), a defect.
    accumulate : str
        ``"product"`` multiplies the step factors; ``"log"`` sums their logarithms (the
        round-off companion).

    Returns
    -------
    float
        V C at t1.
    """
    t = np.linspace(t0, t1, nsteps + 1)
    V = V0 + (V1 - V0) * (t - t0) / (t1 - t0)
    D = equivalent_diameter(V, w)
    rate = np.array([K_of_D(d) for d in D]) * 4.0 * phi / D       # K S/V at each time
    dt = (t1 - t0) / nsteps
    g = (1.0 - (1.0 - theta) * dt * rate[:-1]) / (1.0 + theta * dt * rate[1:])
    f = np.prod(g) if accumulate == "product" else np.exp(np.sum(np.log(g)))
    return float(y0 * f) if dilution else float(y0 / V0 * f * V1)
