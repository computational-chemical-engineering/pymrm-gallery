"""The dusty gas model on a 1-D pore (page A4.3).

Krishna and Wesselingh (1997), eqs. (82)-(87) with the viscous term of
eq. (101): the species fluxes follow from ``[B^e] (N) = -(1/RT) grad p_i``
minus the viscous contribution, solved at every face, while the transport
skeleton is pymrm's finite-volume operators. Two configurations: a capillary
(slab, composition imposed at both ends, no reaction) and a dimensionless
pellet (slab, cylinder or sphere by ``geom``, ``A -> nu_B B`` first order in A).

Two independent routes that use no pymrm operator: ``capillary_bvp`` (eq. 82
as an ODE system with the fluxes as unknown parameters) and ``ms_bvp``
(Maxwell-Stefan without the wall term, closed by Graham's law).
"""
import numpy as np
import scipy.sparse as sps
from scipy.integrate import solve_bvp
from pymrm import (NumJac, construct_boundary_value_matrices, construct_div,
                   construct_grad, newton)

RGAS = 8.314462618  # J/(mol K)

# r = 0 of a pellet (outward normal -r): dc/dr = 0  ->  1*dc/dn + 0*c = 0
SYMMETRY = {"a": 1.0, "b": 0.0, "d": 0.0}


def dirichlet(cvec):
    """Imposed composition at a pore end: ``0*dc/dn + 1*c = c_bnd`` (outward normal).

    Parameters
    ----------
    cvec : array_like
        Concentration of each species at the boundary.

    Returns
    -------
    dict
        Boundary-condition dictionary for pymrm.
    """
    return {"a": 0.0, "b": 1.0, "d": np.asarray(cvec, float).reshape(1, -1)}


def knudsen_D(d0, T, M, eps_tau=1.0, den=3.0, exp=-0.5):
    """Effective Knudsen diffusivity, eq. (85).

    ``D^e_iM = (eps/tau) (d0/den) sqrt(8 R T / pi) M_i**exp`` with the printed
    ``den = 3`` and ``exp = -1/2``; ``d0`` is the pore diameter.

    Parameters
    ----------
    d0 : float
        Pore diameter, m.
    T : float
        Temperature, K.
    M : array_like
        Molar masses, kg/mol.
    eps_tau, den, exp : float
        Porosity over tortuosity, the denominator of ``d0`` and the mass
        exponent (arguments so that each can be broken on purpose).

    Returns
    -------
    numpy.ndarray
        ``D^e_iM`` per species, m2/s.
    """
    return eps_tau * (d0 / den) * np.sqrt(8.0 * RGAS * T / np.pi) * np.asarray(M, float) ** exp


def pair_matrix(pairs, ct, n=3):
    """Symmetric matrix of the pressure-independent products ``D_ij c_t``.

    Parameters
    ----------
    pairs : dict
        ``{(i, j): D_ij}`` in m2/s at total concentration ``ct``.
    ct : float
        Total concentration the diffusivities belong to, mol/m3.
    n : int
        Number of species.

    Returns
    -------
    numpy.ndarray
        ``(n, n)`` array of ``D_ij c_t``, mol/(m s), zero on the diagonal.
    """
    m = np.zeros((n, n))
    for (i, j), v in pairs.items():
        if max(i, j) < n:
            m[i, j] = m[j, i] = v * ct
    return m


class DustyGas:
    """n-component dusty gas model on a 1-D grid; fields of shape ``(n_cell, n)``.

    Parameters
    ----------
    n_cell : int
        Number of cells.
    L : float
        Domain length, m (or 1 for the dimensionless pellet).
    dij_ct : array_like
        ``(n, n)`` products ``D_ij c_t``, pressure independent.
    diM : array_like
        Knudsen diffusivities ``D^e_iM``, eq. (85).
    T : float
        Temperature (``RT`` enters only the viscous term).
    geom : int
        ``nu`` of ``construct_div``: 0 slab, 1 cylinder, 2 sphere.
    bc : tuple of dict
        Boundary conditions on the outward normal, left and right.
    eps_tau : float
        Porosity over tortuosity applied to the pair diffusivities, eq. (83).
    alpha : array_like, optional
        Viscous selectivities ``alpha'_i`` of eq. (101); 1 by default.
    B0, eta : float
        Viscous permeability (eq. 91) and viscosity; ``B0 = 0`` drops the term.
    sink, kvol : array_like, float
        Stoichiometric sink per species and rate constant of a reaction first
        order in species 0.
    b_offdiag : float
        Sign of the off-diagonal elements of eq. (87); the printed value is -1.
    knudsen : bool
        Keep the wall term ``1/D^e_iM``.
    scalar : bool
        Replace ``[B^e]^-1`` by the inverse of its diagonal (a scalar closure).
    face_rule : {"mean", "left"}
        Interior face composition: arithmetic mean of the neighbours, or the
        left cell (first order).
    """

    def __init__(self, n_cell, L, dij_ct, diM, T, geom=0, bc=(None, None), eps_tau=1.0,
                 alpha=None, B0=0.0, eta=1.0, sink=None, kvol=0.0, b_offdiag=-1.0,
                 knudsen=True, scalar=False, face_rule="mean"):
        self.n = len(np.atleast_1d(diM))
        self.n_cell, self.shape = n_cell, (n_cell, self.n)
        self.diM = np.atleast_1d(np.array(diM, float))
        self.T = T
        self.alpha = np.ones(self.n) if alpha is None else np.atleast_1d(np.array(alpha, float))
        self.B0, self.eta = B0, eta
        self.sink = None if sink is None else np.array(sink, float)
        self.kvol = kvol
        self.b_offdiag = b_offdiag
        self.scalar = scalar
        self.x_f = np.linspace(0.0, L, n_cell + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # constant operators, assembled once
        self.grad, gb = construct_grad(self.shape, self.x_f, self.x_c, bc)
        self.grad_bc = gb.toarray().reshape((-1, 1))
        self.div = construct_div(self.shape, self.x_f, nu=geom)
        bl, blb = construct_boundary_value_matrices(self.shape, self.x_f, self.x_c, bc[0], bound_id=0)
        br, brb = construct_boundary_value_matrices(self.shape, self.x_f, self.x_c, bc[1], bound_id=1)
        if face_rule == "mean":
            w = 0.5 * (sps.eye(n_cell - 1, n_cell) + sps.eye(n_cell - 1, n_cell, k=1))
        elif face_rule == "left":
            w = sps.eye(n_cell - 1, n_cell)
        else:
            raise ValueError(f"face_rule {face_rule!r}")
        # face compositions: boundary faces from the bc, interior faces from w
        self.face = sps.vstack([bl, sps.kron(w, sps.eye(self.n)), br], format="csc")
        self.face_bc = np.concatenate([blb.toarray().ravel(), np.zeros((n_cell - 1) * self.n),
                                       brb.toarray().ravel()]).reshape((-1, 1))
        self.numjac = NumJac(self.shape, axes_diagonals=[0], axes_blocks=[-1])
        # eq. (83): eps/tau on the pair diffusivities; eq. (85) carries it in diM
        dij_ct = np.array(dij_ct, float).reshape(self.n, self.n)
        nz = dij_ct > 0
        self.inv_ij = np.zeros_like(dij_ct)
        self.inv_ij[nz] = 1.0 / (eps_tau * dij_ct[nz])
        self.inv_iM = (1.0 / self.diM) if knudsen else np.zeros(self.n)

    def faces(self, c):
        """Face compositions and face gradients, each ``(n_cell + 1, n)``."""
        cc = np.reshape(c, (-1, 1))
        cf = (self.face @ cc + self.face_bc).reshape(self.n_cell + 1, self.n)
        g = (self.grad @ cc + self.grad_bc).reshape(self.n_cell + 1, self.n)
        return cf, g

    def b_matrix(self, cf):
        """Eq. (87) at each face; ``x_k / D^e_ik = c_k / (eps_tau D_ik c_t)``."""
        B = self.b_offdiag * cf[:, :, None] * self.inv_ij[None, :, :]
        d = self.inv_iM[None, :] + (cf[:, None, :] * self.inv_ij[None, :, :]).sum(axis=2)
        i = np.arange(self.n)
        B[:, i, i] = d
        return B

    def flux(self, c):
        """Face fluxes from eq. (86), with eq. (101)'s viscous term when ``B0 > 0``."""
        cf, g = self.faces(c)
        rhs = g.copy()                                    # (1/RT) grad p_i = grad c_i
        if self.B0:
            ct = cf.sum(axis=1)
            d_visc = ct * self.B0 * RGAS * self.T / self.eta          # eq. (106)
            rhs = rhs + (self.alpha[None, :] * (cf / ct[:, None]) * d_visc[:, None]
                         * self.inv_iM[None, :] * g.sum(axis=1)[:, None])
        B = self.b_matrix(cf)
        if self.scalar:
            return -rhs / np.diagonal(B, axis1=1, axis2=2)
        return -np.linalg.solve(B, rhs[..., None])[..., 0]

    def _f(self, c):
        out = self.div @ self.flux(c).reshape((-1, 1))
        if self.kvol:          # sink_i times a rate first order in species 0
            out = out + (self.sink[None, :] * (self.kvol * c[:, 0:1])).reshape((-1, 1))
        return out.reshape(self.shape)

    def residual(self, c):
        """Residual and Jacobian for ``newton``."""
        g, jac = self.numjac(self._f, np.reshape(c, self.shape))
        return g.reshape((-1, 1)), jac

    def solve(self, c0, tol=1e-12):
        """Steady state by Newton from ``c0``; raises unless it converges.

        Parameters
        ----------
        c0 : array_like
            Initial guess, shape ``(n_cell, n)``.
        tol : float
            Newton step tolerance relative to the largest initial concentration
            (round-off is near 1e-15).

        Returns
        -------
        numpy.ndarray
            Concentrations, shape ``(n_cell, n)``.
        """
        c0 = np.asarray(c0, float).reshape(self.shape)
        scale = float(np.max(np.abs(c0)))
        r = newton(self.residual, c0, tol=tol, scale=scale, bounds=(0.0, np.inf),
                   stagnation=True, maxfev=100)
        if not r.success:
            raise RuntimeError(f"dusty gas solve failed: {r.message}")
        return r.x.reshape(self.shape)


def capillary(d0, xa, xb, *, molar, dij_ct, T, p, L, n_cell=100, eps_tau=1.0, kn_den=3.0,
              mass_exp=-0.5, tol=1e-12, **kw):
    """Isothermal slab, composition imposed at both ends, no reaction.

    Parameters
    ----------
    d0 : float
        Pore diameter, m.
    xa, xb : array_like
        Mole fractions at ``z = 0`` and ``z = L``.
    molar, dij_ct, T, p, L : see DustyGas
        Molar masses, pair products ``D_ij c_t``, temperature, pressure, length.
    n_cell, eps_tau, kn_den, mass_exp, tol : optional
        Grid, eq. (83)/(85) constants and the Newton tolerance.
    **kw
        Passed to ``DustyGas``.

    Returns
    -------
    tuple
        ``(model, c)``.
    """
    ct = p / (RGAS * T)
    diM = knudsen_D(d0, T, molar, eps_tau=eps_tau, den=kn_den, exp=mass_exp)
    ca, cb = ct * np.asarray(xa, float), ct * np.asarray(xb, float)
    m = DustyGas(n_cell, L, dij_ct, diM, T, eps_tau=eps_tau,
                 bc=(dirichlet(ca), dirichlet(cb)), **kw)
    c0 = np.outer(np.linspace(0.0, 1.0, n_cell), cb - ca) + ca
    return m, m.solve(np.maximum(c0, 1e-12 * ct), tol=tol)


def pellet(nu_b, kn, phi, x_as=1.0, mass_ratio=None, n_cell=200, dvisc=0.0, n_extra=0,
           x_extra=0.0, geom=2, mass_exp=-0.5, sink_extra=0.0, **kw):
    """Dimensionless pellet, ``A -> nu_B B`` first order in A.

    ``D^e_AM = 1``, ``D^e_ij = kn`` (so Kn = kn), surface total concentration 1,
    ``phi**2 = k R**2 / D^e_AM`` and ``RT = 1``; the viscous group
    ``D_visc / D^e_AM`` is ``dvisc``. ``D^e_BM / D^e_AM = (M_B/M_A)**mass_exp``
    by eq. (85), with ``M_B/M_A = 1/nu_B`` (mass conserved) by default. Inert
    species (``n_extra``, mole fraction ``x_extra`` at the surface) have the
    molar mass of A; ``sink_extra`` consumes them with A, for a defect row.

    Returns
    -------
    tuple
        ``(model, c)``; the solve continues in ``phi`` from 0.5.
    """
    mr = (1.0 / nu_b) if mass_ratio is None else mass_ratio
    n = 2 + n_extra
    diM = np.array([1.0, mr ** mass_exp] + [1.0] * n_extra)
    dij = np.full((n, n), float(kn))
    np.fill_diagonal(dij, 0.0)
    cs = np.maximum(np.array([x_as, 1.0 - x_as - n_extra * x_extra] + [x_extra] * n_extra), 1e-12)
    sink = np.array([1.0, -nu_b] + [sink_extra] * n_extra)
    m = DustyGas(n_cell, 1.0, dij, diM, 1.0 / RGAS, geom=geom, eps_tau=1.0,
                 bc=(SYMMETRY, dirichlet(cs)), B0=dvisc, eta=1.0, sink=sink,
                 kvol=phi ** 2, **kw)
    c = np.outer(np.ones(n_cell), cs)
    for f in np.geomspace(0.5, max(phi, 0.5), 14):   # continuation in phi
        m.kvol = f ** 2
        c = m.solve(c)
    return m, c


def capillary_bvp(d0, xa, xb, *, molar, dij_ct, T, p, L, eps_tau=1.0, kn_den=3.0,
                  mass_exp=-0.5, nz=41, tol=1e-10):
    """Eq. (82) as an ODE system in ``c_i`` with the fluxes as unknown parameters.

    Uses no pymrm operator, no ``[B^e]`` and no finite-volume grid; shares only
    ``knudsen_D`` and the pair products with the finite-volume route.

    Returns
    -------
    numpy.ndarray
        The three fluxes, mol/(m2 s); raises unless ``solve_bvp`` converges.
    """
    ct = p / (RGAS * T)
    dij = np.array(dij_ct, float) / ct * eps_tau        # D^e_ij at this pressure
    diM = knudsen_D(d0, T, molar, eps_tau=eps_tau, den=kn_den, exp=mass_exp)
    ca, cb = ct * np.asarray(xa, float), ct * np.asarray(xb, float)
    n = len(ca)

    def rhs(z, c, N):
        x = c / c.sum(axis=0)
        out = np.empty_like(c)
        for i in range(n):
            s = N[i] / diM[i] + 0.0 * z
            for j in range(n):
                if j != i:
                    s = s + (x[j] * N[i] - x[i] * N[j]) / dij[i, j]
            out[i] = -s                                   # = (1/RT) dp_i/dz, eq. (82)
        return out

    z = np.linspace(0.0, L, nz)
    y0 = ca[:, None] + (cb - ca)[:, None] * (z / L)[None, :]
    sol = solve_bvp(rhs, lambda c0, c1, N: np.concatenate([c0 - ca, c1 - cb]),
                    z, y0, p=-(cb - ca) * diM / L, tol=tol, max_nodes=50000)
    if sol.status != 0:
        raise RuntimeError(f"collocation route failed: {sol.message}")
    return sol.p


def ms_bvp(xa, xb, *, molar, dij_ct, T, p, L, eps_tau=1.0, gr_exp=0.5, tol=1e-11):
    """Maxwell-Stefan without the wall term, in mole fractions, closed by Graham's law.

    ``N_3`` follows from ``sum_i N_i M_i**gr_exp = 0`` (eq. 107). Independent of
    ``DustyGas`` and of ``capillary_bvp``.

    Returns
    -------
    numpy.ndarray
        The three fluxes; raises unless ``solve_bvp`` converges.
    """
    ct = p / (RGAS * T)
    dij = np.array(dij_ct, float) / ct * eps_tau
    wt = np.asarray(molar, float) ** gr_exp
    xa, xb = np.asarray(xa, float), np.asarray(xb, float)

    def fluxes(q):
        return np.array([q[0], q[1], -(q[0] * wt[0] + q[1] * wt[1]) / wt[2]])

    def rhs(z, x, q):
        N = fluxes(q)
        xs = [x[0], x[1], 1.0 - x[0] - x[1]]
        out = np.empty_like(x)
        for i in range(2):
            s = 0.0 * z
            for j in range(3):
                if j != i:
                    s = s + (xs[j] * N[i] - xs[i] * N[j]) / dij[i, j]
            out[i] = -s / ct
        return out

    z = np.linspace(0.0, L, 41)
    y0 = xa[:2, None] + (xb - xa)[:2, None] * (z / L)[None, :]
    sol = solve_bvp(rhs, lambda a, b, q: np.array([a[0] - xa[0], a[1] - xa[1],
                                                   b[0] - xb[0], b[1] - xb[1]]),
                    z, y0, p=np.array([1e-3, -1e-4]), tol=tol, max_nodes=50000)
    if sol.status != 0:
        raise RuntimeError(f"Maxwell-Stefan collocation failed: {sol.message}")
    return fluxes(sol.p)
