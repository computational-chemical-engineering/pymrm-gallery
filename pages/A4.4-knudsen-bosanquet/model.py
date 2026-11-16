"""The binary dusty gas model and the Bosanquet relation on a 1-D grid (page A4.4).

Krishna and Wesselingh (1997), eqs. (82)-(87) for two species without the
viscous term, with three flux closures on one finite-volume skeleton: the
2 x 2 dusty-gas matrix ``[B^e]`` solved at every face (no flux ratio assumed),
the scalar effective diffusivity of eq. (109) with a prescribed flux ratio, and
the Bosanquet formula, eq. (110). Two configurations: a slab with the
composition imposed at both ends (no reaction) and a dimensionless pellet
(slab, cylinder or sphere by ``geom``) running ``A -> nu_B B`` first order in A.

Two routes that use no pymrm operator: ``isobaric_flux_closed_form`` (eq. 109
with Graham's flux ratio integrated across the slab) and ``eta_classical`` (the
first-order effectiveness factor at constant diffusivity).
"""
import numpy as np
import scipy.sparse as sps
from scipy.special import iv
from pymrm import (NumJac, construct_boundary_value_matrices, construct_div,
                   construct_grad, newton)

RGAS = 8.314462618  # J/(mol K)

# r = 0 of a pellet (outward normal -r): dc/dr = 0  ->  1*dc/dn + 0*c = 0
SYMMETRY = {"a": 1.0, "b": 0.0, "d": 0.0}


def dirichlet(ca, cb):
    """Imposed composition at a boundary: ``0*dc/dn + 1*c = (c_A, c_B)``.

    Parameters
    ----------
    ca, cb : float
        Concentrations of A and B on the boundary face.

    Returns
    -------
    dict
        Boundary-condition dictionary for pymrm (outward normal).
    """
    return {"a": 0.0, "b": 1.0, "d": np.array([[ca, cb]], float)}


def knudsen_D(d0, T, M, eps_tau=1.0, den=3.0, exp=-0.5):
    """Effective Knudsen diffusivity, eq. (85).

    ``D^e_iM = (eps/tau) (d0/den) sqrt(8 R T / pi) M_i**exp`` with the printed
    ``den = 3`` and ``exp = -1/2``; ``d0`` is the pore diameter.

    Parameters
    ----------
    d0 : float or array_like
        Pore diameter, m.
    T : float
        Temperature, K.
    M : float or array_like
        Molar mass, kg/mol.
    eps_tau, den, exp : float
        Porosity over tortuosity, the denominator of ``d0`` and the mass
        exponent (arguments so that each can be broken on purpose).

    Returns
    -------
    float or numpy.ndarray
        ``D^e_iM``, m2/s.
    """
    return eps_tau * (np.asarray(d0, float) / den) * np.sqrt(8.0 * RGAS * T / np.pi) * \
        np.asarray(M, float) ** exp


def d_groups(kn):
    """``(D^e_AM, D^e_AB)`` scaled by ``D_Bos``, given ``Kn = D^e_AB / D^e_AM``.

    With this scaling ``1/D^e_AM + 1/D^e_AB = 1``, so the Bosanquet diffusivity
    of A is 1 at every Kn.
    """
    return 1.0 + 1.0 / kn, 1.0 + kn


class Cell:
    """Binary porous medium on a 1-D grid; fields of shape ``(n, 2) = (c_A, c_B)``.

    Parameters
    ----------
    n : int
        Number of cells on ``[0, 1]`` (or ``[0, L]``).
    geom : int
        ``nu`` of ``construct_div``: 0 slab, 1 cylinder, 2 sphere.
    bc : tuple of dict
        Boundary conditions on the outward normal, left and right.
    dam, dbm, dab : float
        ``D^e_AM``, ``D^e_BM`` (eq. 85) and ``D^e_AB`` (eq. 83).
    phi : float
        Thiele modulus of ``A -> nu_B B`` first order in A; 0 switches it off.
    nu_b : float
        Stoichiometric coefficient of B in the source term.
    model : {"dgm", "eq109", "bosanquet"}
        Flux closure: the 2 x 2 matrix of eq. (87), eq. (109) with flux ratio
        ``rho``, or eq. (110).
    rho : float, optional
        Flux ratio ``-N_B/N_A`` handed to eq. (109); ``nu_b`` by default.
    b_offdiag : float
        Factor on the off-diagonals of eq. (87); the printed value is 1.
    nu_swap : bool
        Read eq. (109)'s stoichiometric ratio upside down (``x_A/rho`` in place
        of ``rho x_A``), the transcription error the OCR invites.
    L : float
        Domain length.
    face_rule : {"mean", "left"}
        Interior face composition: arithmetic mean of the neighbours, or the
        left cell (first order).
    """

    def __init__(self, n=200, geom=0, bc=(None, None), dam=1.0, dbm=1.0, dab=1.0,
                 phi=0.0, nu_b=1.0, model="dgm", rho=None, b_offdiag=1.0, nu_swap=False,
                 L=1.0, face_rule="mean"):
        if model not in ("dgm", "eq109", "bosanquet"):
            raise ValueError(f"model {model!r}")
        self.n, self.geom, self.shape = n, geom, (n, 2)
        self.dam, self.dbm, self.dab = dam, dbm, dab
        self.phi, self.nu_b, self.model = phi, nu_b, model
        self.rho = nu_b if rho is None else rho
        self.b_offdiag, self.nu_swap = b_offdiag, nu_swap
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # constant operators, assembled once
        self.grad, gb = construct_grad(self.shape, self.x_f, self.x_c, bc)
        self.grad_bc = gb.toarray().reshape((-1, 1))
        self.div = construct_div(self.shape, self.x_f, nu=geom)
        bl, blb = construct_boundary_value_matrices(self.shape, self.x_f, self.x_c, bc[0],
                                                    bound_id=0)
        br, brb = construct_boundary_value_matrices(self.shape, self.x_f, self.x_c, bc[1],
                                                    bound_id=1)
        # face compositions: boundary faces from the bc, interior faces from w
        if face_rule == "mean":
            w = 0.5 * (sps.eye(n - 1, n) + sps.eye(n - 1, n, k=1))
        elif face_rule == "left":
            w = sps.eye(n - 1, n)
        else:
            raise ValueError(f"face_rule {face_rule!r}")
        self.face = sps.vstack([bl, sps.kron(w, sps.eye(2)), br]).tocsr()
        self.face_bc = np.concatenate([blb.toarray().ravel(), np.zeros(2 * (n - 1)),
                                       brb.toarray().ravel()]).reshape((-1, 1))
        self.numjac = NumJac(self.shape, axes_diagonals=[0], axes_blocks=[-1])

    def faces(self, c):
        """Face compositions and face gradients, each ``(n + 1, 2)``."""
        cc = np.reshape(c, (-1, 1))
        cf = (self.face @ cc + self.face_bc).reshape(self.n + 1, 2)
        g = (self.grad @ cc + self.grad_bc).reshape(self.n + 1, 2)
        return cf, g

    def b_matrix(self, xf):
        """Eq. (87) for two species: ``B_ii = 1/D_iM + x_k/D_ik``, ``B_ij = -x_i/D_ij``."""
        xa, xb = xf[:, 0], xf[:, 1]
        b = np.empty(xf.shape[:-1] + (2, 2))
        b[..., 0, 0] = 1.0 / self.dam + xb / self.dab
        b[..., 0, 1] = -self.b_offdiag * xa / self.dab
        b[..., 1, 0] = -self.b_offdiag * xb / self.dab
        b[..., 1, 1] = 1.0 / self.dbm + xa / self.dab
        return b

    def d_scalar(self, xf, bosanquet=False):
        """Eq. (109) without the viscous term, or eq. (110) if ``bosanquet``.

        With ``N_B/N_A = -rho``: ``1/D_A = 1/D_AM + (x_B + rho x_A)/D_AB`` and
        ``1/D_B = 1/D_BM + (x_A + x_B/rho)/D_AB``.

        Returns
        -------
        numpy.ndarray
            ``(n_faces, 2)`` diffusivities of A and B.
        """
        xa, xb, r = xf[:, 0], xf[:, 1], self.rho
        if bosanquet:
            ia = np.full_like(xa, 1.0 / self.dam + 1.0 / self.dab)
            ib = np.full_like(xb, 1.0 / self.dbm + 1.0 / self.dab)
        elif self.nu_swap:
            ia = 1.0 / self.dam + (xb + xa / r) / self.dab
            ib = 1.0 / self.dbm + (xa + xb * r) / self.dab
        else:
            ia = 1.0 / self.dam + (xb + r * xa) / self.dab
            ib = 1.0 / self.dbm + (xa + xb / r) / self.dab
        return np.stack([1.0 / ia, 1.0 / ib], axis=1)

    def flux(self, c):
        """Face fluxes ``(n + 1, 2)`` from the chosen closure (ideal gas, grad p_i = RT grad c_i)."""
        cf, g = self.faces(c)
        xf = cf / cf.sum(axis=1, keepdims=True)
        if self.model == "dgm":       # eq. (86): (N) = -[B^e]^-1 (grad c)
            return -np.linalg.solve(self.b_matrix(xf), g[..., None])[..., 0]
        return -self.d_scalar(xf, self.model == "bosanquet") * g

    def _f(self, c):
        out = self.div @ self.flux(c).reshape((-1, 1))
        if self.phi:          # A -> nu_B B, first order in A
            r = self.phi ** 2 * c[:, 0:1]
            out = out + np.concatenate([r, -self.nu_b * r], axis=1).reshape((-1, 1))
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
            Initial guess, shape ``(n, 2)``.
        tol : float
            Newton step tolerance relative to the largest initial concentration
            (round-off is near 1e-15).

        Returns
        -------
        numpy.ndarray
            Concentrations, shape ``(n, 2)``.
        """
        c0 = np.asarray(c0, float).reshape(self.shape)
        scale = float(np.max(np.abs(c0)))
        r = newton(self.residual, c0, tol=tol, scale=scale, bounds=(0.0, np.inf),
                   stagnation=True, maxfev=100)
        if not r.success:
            raise RuntimeError(f"binary dusty gas solve failed: {r.message}")
        return r.x.reshape(self.shape)

    def centre_total(self, c):
        """Total concentration on the left boundary face (the pellet centre)."""
        return float(self.faces(c)[0][0].sum())

    def eta(self, c):
        """Effectiveness factor ``(nu + 1) (-N_A at the surface face) / phi**2``."""
        return float((self.geom + 1.0) * (-self.flux(c)[-1, 0]) / self.phi ** 2)


def pellet(kn, nu_b, phi, n=600, geom=2, model="dgm", mass_exp=-0.5, nu_src=None,
           b_offdiag=1.0, nu_swap=False, face_rule="mean", tol=1e-12):
    """Dimensionless pellet, pure A outside, ``A -> nu_B B`` first order in A.

    Mass conservation fixes ``M_B = M_A/nu_B``, so eq. (85) gives
    ``D^e_BM/D^e_AM = (M_B/M_A)**mass_exp``. Concentrations are scaled by the
    surface value of A; diffusivities by ``D_Bos`` of A (``d_groups``).

    Parameters
    ----------
    kn, nu_b, phi : float
        Knudsen number ``D^e_AB/D^e_AM``, stoichiometric coefficient, Thiele modulus.
    n, geom, model : optional
        Cells, geometry index and closure (see ``Cell``).
    mass_exp : float
        Eq. (85)'s mass exponent (printed -1/2).
    nu_src : float, optional
        Coefficient of B in the source term; ``nu_b`` by default.
    b_offdiag, nu_swap, face_rule, tol : optional
        Passed to ``Cell`` and ``Cell.solve``.

    Returns
    -------
    tuple
        ``(cell, c)``.
    """
    dam, dab = d_groups(kn)
    p = Cell(n=n, geom=geom, dam=dam, dbm=dam * (1.0 / nu_b) ** mass_exp, dab=dab,
             phi=phi, nu_b=nu_b if nu_src is None else nu_src, model=model, rho=nu_b,
             b_offdiag=b_offdiag, nu_swap=nu_swap, face_rule=face_rule,
             bc=(SYMMETRY,                  # centre: dc/dr = 0
                 dirichlet(1.0, 0.0)))       # surface: c_A = 1, c_B = 0
    c0 = np.zeros(p.shape)
    c0[:, 0] = 1.0
    return p, p.solve(c0, tol=tol)


def isobaric_slab(kn, mass_ratio, xa=(0.9, 0.1), model="dgm", n=200, mass_exp=-0.5,
                  rho=None, b_offdiag=1.0, ct=1.0, dab=None, dam=None, L=1.0):
    """Binary slab, equal total concentration at both ends, no reaction.

    Parameters
    ----------
    kn : float
        Knudsen number ``D^e_AB/D^e_AM`` (diffusivities scaled by ``D_Bos``),
        unless ``dam`` and ``dab`` are given.
    mass_ratio : float
        ``M_B/M_A``; eq. (85) gives ``D^e_BM/D^e_AM = mass_ratio**mass_exp``.
    xa : tuple of float
        Mole fraction of A at ``z = 0`` and ``z = L``.
    model, n, rho, b_offdiag : optional
        Closure, cells, flux ratio for eq. (109), eq. (87) off-diagonal factor.
    mass_exp : float
        Eq. (85)'s mass exponent.
    ct, dab, dam, L : optional
        Total concentration, dimensional diffusivities and length.

    Returns
    -------
    tuple
        ``(cell, c)``.
    """
    if dam is None:
        dam, dab = d_groups(kn)
    x1, x2 = xa
    s = Cell(n=n, geom=0, dam=dam, dbm=dam * mass_ratio ** mass_exp, dab=dab, model=model,
             rho=rho, b_offdiag=b_offdiag, L=L,
             bc=(dirichlet(x1 * ct, (1 - x1) * ct),        # z = 0: imposed composition
                 dirichlet(x2 * ct, (1 - x2) * ct)))       # z = L: imposed composition
    xs = np.linspace(x1, x2, n)
    c0 = ct * np.stack([xs, 1.0 - xs], axis=1)
    return s, s.solve(c0)


def isobaric_flux_closed_form(dam, dab, rho, x1, x2, ct=1.0, L=1.0):
    """Flux of A from eq. (109) with a constant flux ratio, integrated across the slab.

    At uniform total concentration ``N_A L = ct * int_{x2}^{x1} D_A(x) dx`` with
    ``1/D_A = a + b x``, ``a = 1/D_AM + 1/D_AB``, ``b = (rho - 1)/D_AB``. Uses no
    pymrm operator and no grid.

    Returns
    -------
    float
        ``N_A``.
    """
    a, b = 1.0 / dam + 1.0 / dab, (rho - 1.0) / dab
    if b == 0.0:
        return ct * (x1 - x2) / (a * L)
    return ct * np.log((a + b * x1) / (a + b * x2)) / (b * L)


def eta_classical(geom, phi):
    """First-order effectiveness factor at constant diffusivity, per geometry index."""
    if geom == 0:
        return np.tanh(phi) / phi
    if geom == 1:
        return 2.0 * iv(1, phi) / (phi * iv(0, phi))
    return (3.0 / phi) * (1.0 / np.tanh(phi) - 1.0 / phi)
