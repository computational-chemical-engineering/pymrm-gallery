"""Wilke mixture rule and exact ternary Stefan-Maxwell film (page A4.1).

Methods 1 and 2 of Wilke (1950) as root-found film equations, three exact
routes to the same flux pair (closed form, matrix exponential, pymrm finite
volume), all for a ternary gas film with species C stagnant.
"""
import numpy as np
import scipy.sparse as sp
from pymrm import (construct_grad, construct_div, compute_boundary_values, NumJac,
                   newton)
from pymrm.checks import find_roots

P_MM, XLEN, T_K, R_GAS = 150.0, 0.1, 313.15, 82.06      # mm Hg, cm, K, cc atm/mol/K
P_ATM = P_MM / 760.0
C_TOT = P_ATM / (R_GAS * T_K)                            # mol/cm3
SCALE = 1e5                                              # Table 1's multiplier
N_SCAN = 1500      # log-grid points per sign and root scan (see roots_both_signs)

BASE = None   # Table-1 base rows indexed by example; set by set_table()


def set_table(base):
    """Give the model the Table-1 base rows that `unpack` reads.

    Parameters
    ----------
    base : pandas.DataFrame
        One row per example (index: example label) with the columns gas_A,
        gas_B, gas_C, p1_A..p2_C (mm Hg) and D_AB, D_AC, D_BC (cm2/s).
    """
    global BASE
    BASE = base


# =========================================================== problem definition
def unpack(e):
    """(names, y1, y2, D_AB, D_AC, D_BC) for a Table-1 example, mole fractions.

    Parameters
    ----------
    e : str
        Example label of Table 1.

    Returns
    -------
    tuple
        (names, y1, y2, D_AB, D_AC, D_BC): gas names, terminal mole
        fractions at the two film ends, and the three binary diffusivities (cm2/s).
    """
    r = BASE.loc[e]
    y1 = np.array([r.p1_A, r.p1_B, r.p1_C]) / P_MM
    y2 = np.array([r.p2_A, r.p2_B, r.p2_C]) / P_MM
    return (r.gas_A, r.gas_B, r.gas_C), y1, y2, r.D_AB, r.D_AC, r.D_BC


def Dmat(DAB, DAC, DBC):
    """Symmetric matrix of the three binary diffusivities.

    Parameters
    ----------
    DAB, DAC, DBC : float
        Binary diffusivities (cm2/s).

    Returns
    -------
    numpy.ndarray
        Symmetric 3 x 3 matrix of binary diffusivities, NaN on the diagonal.
    """
    return np.array([[np.nan, DAB, DAC], [DAB, np.nan, DBC], [DAC, DBC, np.nan]])


# ================================================== the two mixture rules
def wilke_rule(y, i, D):
    """Eq. (6)/(25): weighted HARMONIC mean of the binaries.

    Parameters
    ----------
    y : array_like
        Mole fractions (3,).
    i : int
        Index of the diffusing species.
    D : numpy.ndarray
        Binary diffusivity matrix from `Dmat` (cm2/s).

    Returns
    -------
    float
        Effective diffusivity of species i (cm2/s).
    """
    return (1.0 - y[i]) / sum(y[j] / D[i, j] for j in range(3) if j != i)


def arith_rule(y, i, D):
    """Null baseline: weighted ARITHMETIC mean.  The paper's own printed
    Discussion (J. H. Arnold, journal p. 104) says Wilke's expression 'is simply
    a weighted harmonic mean, as opposed to the weighted arithmetic mean
    suggested by Hougen and Watson; the latter is badly in error in some cases'.
    Same weights y_j, arithmetic instead of harmonic.  Hougen and Watson's own
    text was NOT consulted; this is the form Arnold's sentence names.

    Parameters
    ----------
    y : array_like
        Mole fractions (3,).
    i : int
        Index of the diffusing species.
    D : numpy.ndarray
        Binary diffusivity matrix from `Dmat` (cm2/s).

    Returns
    -------
    float
        Effective diffusivity of species i (cm2/s).
    """
    return sum(y[j] * D[i, j] for j in range(3) if j != i) / (1.0 - y[i])


# ================================================== exact route 1: closed form
def exact_closed(NA, NB, y1, DAB, DAC, DBC, x=XLEN, c=C_TOT):
    """(y_A, y_B) at x from the analytic solution of the linear SM system.

    Parameters
    ----------
    NA, NB : array_like
        Fluxes of A and B (mol/cm2/s).
    y1 : array_like
        Mole fractions at x = 0 (3,).
    DAB, DAC, DBC : float
        Binary diffusivities (cm2/s).
    x : float, optional
        Position (cm); the film length by default.
    c : float, optional
        Total molar concentration P/RT (mol/cm3).

    Returns
    -------
    tuple of numpy.ndarray
        (y_A, y_B) at x; NaN where the closed form is singular (N_t = 0).
    """
    NA, NB = np.asarray(NA, float), np.asarray(NB, float)
    yA1, yB1, yC1 = y1
    Nt = NA + NB
    lamC = (NA / DAC + NB / DBC) / c               # the y_C eigenmode
    lam2 = Nt / (c * DAB)                          # the A-B eigenmode
    Q = NA * NB * (1.0 / DAC - 1.0 / DBC)
    yC = yC1 * np.exp(lamC * x)
    d = lamC - lam2
    fp = np.where(d != 0, -Q * yC1 / (c * np.where(d == 0, 1.0, d)), np.nan)
    f = (NB * yA1 - NA * yB1 - fp) * np.exp(lam2 * x) + fp * np.exp(lamC * x)
    yA = np.where(Nt != 0, (f + NA * (1.0 - yC)) / np.where(Nt == 0, 1.0, Nt), np.nan)
    return yA, 1.0 - yC - yA


# ============================================ exact route 2: matrix exponential
def _matfun2(M, t, f, fp):
    """f(M t) for a batch of 2x2 matrices, from the two eigenvalues.

    Any analytic f satisfies f(A) = f0 I + f1 A with
    f1 = (f(l1)-f(l2))/(l1-l2) and f0 = (l1 f(l2) - l2 f(l1))/(l1-l2),
    and the confluent limit uses f'.  No inversion of M anywhere, which is what
    lets the equimolal case N_t = 0 -- where M is exactly singular -- be
    evaluated instead of returning NaN.

    Parameters
    ----------
    M : numpy.ndarray
        Batch of 2 x 2 matrices (..., 2, 2).
    t : float
        Scalar multiplying M.
    f, fp : callable
        The function and its derivative, vectorised over complex arrays.

    Returns
    -------
    numpy.ndarray
        Real part of f(M t), shape (..., 2, 2).
    """
    A = M * t
    a, b, c, d = A[..., 0, 0], A[..., 0, 1], A[..., 1, 0], A[..., 1, 1]
    tr, det = a + d, a * d - b * c
    disc = np.sqrt((tr / 2) ** 2 - det + 0j)
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    dl = l1 - l2
    near = np.abs(dl) < 1e-12 * (np.abs(l1) + np.abs(l2) + 1e-300)
    dls = np.where(near, 1.0, dl)
    F1, F2 = f(l1), f(l2)
    f1 = np.where(near, fp(l1), (F1 - F2) / dls)
    f0 = np.where(near, F1 - l1 * fp(l1), (l1 * F2 - l2 * F1) / dls)
    I = np.zeros(A.shape, dtype=complex)
    I[..., 0, 0] = I[..., 1, 1] = 1.0
    return (f0[..., None, None] * I + f1[..., None, None] * A).real


def _phi1(z):
    """(e^z - 1)/z, analytic at 0.

    Parameters
    ----------
    z : numpy.ndarray
        Argument (complex allowed).

    Returns
    -------
    numpy.ndarray
        phi_1(z).
    """
    small = np.abs(z) < 1e-8
    zz = np.where(small, 1.0, z)
    return np.where(small, 1.0 + z / 2 + z ** 2 / 6, np.expm1(zz) / zz)


def _phi1p(z):
    """Derivative of `_phi1`.

    Parameters
    ----------
    z : numpy.ndarray
        Argument (complex allowed).

    Returns
    -------
    numpy.ndarray
        d phi_1 / dz.
    """
    small = np.abs(z) < 1e-6
    zz = np.where(small, 1.0, z)
    return np.where(small, 0.5 + z / 3, (zz * np.exp(zz) - np.expm1(zz)) / zz ** 2)


def exact_expm(NA, NB, y1, DAB, DAC, DBC, x=XLEN, c=C_TOT):
    """Exact film composition from the matrix exponential.

    Parameters
    ----------
    NA, NB : array_like
        Fluxes of A and B (mol/cm2/s).
    y1 : array_like
        Mole fractions at x = 0 (3,).
    DAB, DAC, DBC : float
        Binary diffusivities (cm2/s).
    x : float, optional
        Position (cm); the film length by default.
    c : float, optional
        Total molar concentration P/RT (mol/cm3).

    Returns
    -------
    tuple of numpy.ndarray
        (y_A, y_B) at x. Finite at N_t = 0, where `exact_closed` is not.
    """
    NA = np.atleast_1d(np.asarray(NA, float))
    NB = np.atleast_1d(np.asarray(NB, float))
    M = np.empty(NA.shape + (2, 2))
    M[..., 0, 0] = NA / DAC + NB / DAB
    M[..., 0, 1] = NA / DAC - NA / DAB
    M[..., 1, 0] = NB / DBC - NB / DAB
    M[..., 1, 1] = NB / DBC + NA / DAB
    M /= c
    b = np.stack([-NA / DAC, -NB / DBC], -1) / c
    # variation of constants without inverting M:
    #     y(x) = exp(Mx) y_1 + x phi_1(Mx) b,   phi_1(z) = (e^z - 1)/z
    E = _matfun2(M, x, np.exp, np.exp)
    Ph = _matfun2(M, x, _phi1, _phi1p)
    ye = (np.einsum("...ij,...j->...i", E, y1[:2])
          + x * np.einsum("...ij,...j->...i", Ph, b))
    return ye[..., 0], ye[..., 1]


# ================================================================ root machinery
def roots_both_signs(f, lo, hi, n_scan=N_SCAN, xtol=1e-14):
    """Every root of a scalar f with lo <= |z| <= hi, of either sign.

    Each sign is scanned by pymrm.checks.find_roots on its own logarithmic grid,
    so roots spread over many decades are bracketed; z = 0 is left out, where the
    film equations are singular. Callers decide how many roots they expect.

    Parameters
    ----------
    f : callable
        Scalar function of one variable; nan where it cannot be evaluated.
    lo, hi : float
        Smallest and largest |z| scanned (both positive).
    n_scan : int, optional
        Scan points per sign.
    xtol : float, optional
        Absolute tolerance of each refined root.

    Returns
    -------
    list of float
        All roots, sorted.
    """
    neg = find_roots(lambda z: f(-z), lo, hi, n_scan=n_scan, log=True, xtol=xtol)
    pos = find_roots(f, lo, hi, n_scan=n_scan, log=True, xtol=xtol)
    return sorted([-z for z in neg["roots"]] + pos["roots"])


def single(roots, what):
    """The one element of a root list that must hold exactly one.

    Parameters
    ----------
    roots : list
        Roots returned by a scan.
    what : str
        Name of the problem, for the error message.

    Returns
    -------
    object
        The single root.
    """
    if len(roots) != 1:
        raise RuntimeError(f"{what}: {len(roots)} roots in the scanned range, expected one")
    return roots[0]


def eq23_NB(NA, y1, y2, DAC, DBC, x=XLEN, c=C_TOT):
    """Wilke eq. (23) solved for N_B - an exact integral, not an approximation.

    Parameters
    ----------
    NA : array_like
        Flux of A (mol/cm2/s).
    y1, y2 : array_like
        Mole fractions at the two film ends (3,).
    DAC, DBC : float
        Binary diffusivities (cm2/s).
    x : float, optional
        Film length (cm).
    c : float, optional
        Total molar concentration P/RT (mol/cm3).

    Returns
    -------
    numpy.ndarray
        Flux of B (mol/cm2/s).
    """
    rhs = c / x * np.log(y2[2] / y1[2])
    return (rhs - np.asarray(NA, float) / DAC) * DBC


def exact_solve(e, route=exact_closed, lim=2e-2, x=XLEN, binaries=None, c=C_TOT,
                n_scan=N_SCAN):
    """All exact (N_A, N_B) for a Table-1 example.  Eq. (23) removes one unknown,
    so this is a 1-D scan over N_A: deterministic, no seed, no continuation, and
    it reports how many roots exist rather than assuming one.

    Parameters
    ----------
    e : str
        Example label of Table 1.
    route : callable, optional
        Exact composition route (`exact_closed` or `exact_expm`).
    lim : float, optional
        Largest |N_A| scanned (mol/cm2/s); the smallest is 1e-12.
    x : float, optional
        Film length (cm).
    binaries : tuple of float, optional
        Replacement (D_AB, D_AC, D_BC) in cm2/s.
    c : float, optional
        Total molar concentration P/RT (mol/cm3).
    n_scan : int, optional
        Scan points per sign of N_A.

    Returns
    -------
    list of tuple
        One (N_A, N_B) pair per root.
    """
    _, y1, y2, DAB, DAC, DBC = unpack(e)
    if binaries is not None:
        DAB, DAC, DBC = binaries

    def F(u):                      # u = N_A * SCALE, the flux in Table 1's units
        NA = np.atleast_1d(u / SCALE)
        NB = eq23_NB(NA, y1, y2, DAC, DBC, x, c)
        return float(np.ravel(route(NA, NB, y1, DAB, DAC, DBC, x, c)[0])[0] - y2[0])

    return [(u / SCALE, float(eq23_NB(u / SCALE, y1, y2, DAC, DBC, x, c)))
            for u in roots_both_signs(F, 1e-12 * SCALE, lim * SCALE, n_scan)]


# ===================================== exact route 3: the pymrm finite volume
class SMFilm:
    """Exact ternary Maxwell-Stefan film with species C stagnant, on a pymrm grid.

    State is y = (y_A, y_B) with shape (n_z, 2): SPATIAL AXIS FIRST, FIELD LAST.
    The face closure is N = -c [B]^-1 grad(y), with

        [B] = [[ y_B/D_AB + y_C/D_AC ,  -y_A/D_AB          ],
               [ -y_B/D_AB           ,  y_A/D_AB + y_C/D_BC]]

    which is the same 2x2 that appears in the analytic route, so the two are NOT
    independent in the matrix; they are independent in everything downstream of
    it (integration versus discretisation).  The residual is nonlinear in y
    through [B]^-1, so it is solved by pymrm.newton.  At steady state div N = 0,
    so N is the same at every face: that uniformity is the solve's residual
    re-expressed, and the iteration count never is a convergence assertion.
    """

    def __init__(s, e, n_z=80, binaries=None, c=C_TOT, nu=0, face="mean"):
        """Set up the grid, operators and boundary conditions of one example.

        Parameters
        ----------
        e : str
            Example label of Table 1.
        n_z : int, optional
            Number of cells across the film.
        binaries : tuple of float, optional
            Replacement (D_AB, D_AC, D_BC) in cm2/s.
        c : float, optional
            Total molar concentration P/RT (mol/cm3).
        nu : int, optional
            Geometry of the divergence: 0 is the Cartesian film.
        face : {'mean', 'left'}, optional
            Interior face composition: the mean of the two cells (second order)
            or the left cell's value (first order; a break-table defect).
        """
        _, y1, y2, DAB, DAC, DBC = unpack(e)
        if binaries is not None:
            DAB, DAC, DBC = binaries
        s.DAB, s.DAC, s.DBC, s.c = DAB, DAC, DBC, c
        s.y1, s.y2, s.n_z, s.face = y1, y2, n_z, face
        s.x_f = np.linspace(0.0, XLEN, n_z + 1)
        s.x_c = 0.5 * (s.x_f[:-1] + s.x_f[1:])
        s.shape = (n_z, 2)                              # spatial axis first, never (n,)
        s.div = construct_div(s.shape, s.x_f, nu=nu)    # nu = 0: Cartesian film
        nf, nc = n_z + 1, 2
        s._cols = (np.arange(nf)[:, None] * nc
                   + np.arange(nc)[None, :]).repeat(nc, 0).ravel()
        s._ptr = np.arange(nf * nc + 1) * nc
        s._shape2 = (nf * nc, nf * nc)
        # OUTWARD normal n: a dy/dn + b y = d.  Both ends are Dirichlet, the film's
        # terminal compositions:
        #   x = 0    : y = y_1  ->  a = 0, b = 1, d = y_1
        #   x = XLEN : y = y_2  ->  a = 0, b = 1, d = y_2
        s.bc = ({"a": 0.0, "b": 1.0, "d": y1[None, :2]},
                {"a": 0.0, "b": 1.0, "d": y2[None, :2]})
        s.grad, gb = construct_grad(s.shape, s.x_f, s.x_c, s.bc)
        s.gbc = np.asarray(gb.todense()).reshape(-1)
        # the residual reads NEIGHBOURING cells through the face mobility c [B]^-1,
        # so the space axis is tridiagonal; the two species couple in full
        s.numjac = NumJac(s.shape, axes_diagonals=[0])

    def _faces(s, y):
        """Face compositions: interior faces the mean of the two cells (uniform
        grid; or the left cell with face='left'), the two boundary faces from compute_boundary_values under s.bc.

        Parameters
        ----------
        y : numpy.ndarray
            Cell compositions, shape (n_z, 2).

        Returns
        -------
        numpy.ndarray
            Compositions at the n_z + 1 faces, shape (n_z + 1, 2).
        """
        yf = np.empty((s.n_z + 1, 2))
        yf[1:-1] = 0.5 * (y[:-1] + y[1:]) if s.face == "mean" else y[:-1]
        y_left, _, y_right, _ = compute_boundary_values(y, s.x_f, s.x_c, s.bc)
        yf[0], yf[-1] = np.ravel(y_left), np.ravel(y_right)
        return yf

    def _M(s, y):
        """Face-wise mobility matrix c [B]^-1.

        Parameters
        ----------
        y : numpy.ndarray
            Cell compositions, shape (n_z, 2).

        Returns
        -------
        scipy.sparse.csr_array
            Block-diagonal matrix of the 2 x 2 face blocks.
        """
        yf = s._faces(y)
        yA, yB = yf[:, 0], yf[:, 1]
        yC = 1.0 - yA - yB
        B = np.empty((s.n_z + 1, 2, 2))
        B[:, 0, 0] = yB / s.DAB + yC / s.DAC
        B[:, 0, 1] = -yA / s.DAB
        B[:, 1, 0] = -yB / s.DAB
        B[:, 1, 1] = yA / s.DAB + yC / s.DBC
        M = s.c * np.linalg.inv(B)
        return sp.csr_array((M.ravel(), s._cols, s._ptr), shape=s._shape2)

    def flux(s, y):
        """Face fluxes N = -c [B]^-1 grad(y).

        Parameters
        ----------
        y : numpy.ndarray
            Cell compositions, shape (n_z, 2).

        Returns
        -------
        numpy.ndarray
            Fluxes of A and B at the n_z + 1 faces (mol/cm2/s).
        """
        return -(s._M(y) @ (s.grad @ y.reshape(-1) + s.gbc)).reshape(-1, 2)

    def residual(s, y):
        """Steady-state residual -div N.

        Parameters
        ----------
        y : numpy.ndarray
            Cell compositions, flat or shape (n_z, 2).

        Returns
        -------
        numpy.ndarray
            Residual, shape (n_z, 2).
        """
        y = y.reshape(s.n_z, 2)
        return -(s.div @ (s._M(y) @ (s.grad @ y.reshape(-1) + s.gbc))).reshape(s.n_z, 2)

    def _newton_fun(s, y):
        """Residual and its NumJac Jacobian, for pymrm.newton."""
        return s.numjac(s.residual, y.reshape(s.shape))

    def solve(s, tol=1e-12):
        """Newton solve of the steady film; sets its, y, res, flux_spread, N.

        Parameters
        ----------
        tol : float, optional
            Newton step tolerance on the mole fractions (absolute; they are of
            order one, and 1e-12 sits above their round-off).

        Returns
        -------
        numpy.ndarray
            Fluxes (N_A, N_B), mean over the faces (mol/cm2/s).
        """
        w = (s.x_c / XLEN)[:, None]
        y0 = (1 - w) * s.y1[:2] + w * s.y2[:2]          # linear, deterministic start
        res = newton(s._newton_fun, y0, tol=tol, stagnation=True)
        if not res.success:
            raise RuntimeError(f"SMFilm Newton solve failed: {res.message}")
        y = res.x.reshape(s.shape)
        s.its, s.y = res.nit, y
        s.res = float(np.max(np.abs(s.residual(y))))
        N = s.flux(y)
        s.flux_spread = float(np.max(np.abs(N / N.mean(0) - 1.0)))
        s.N = N.mean(0)
        return s.N


# ================================ Methods 1 and 2: the film equations, root-found
def g(s, y1, y2):
    """phi * ln[(phi - y2)/(phi - y1)] written with s = 1/phi, so that the
    equimolal limit N_t -> 0 (phi -> infinity, s -> 0) is finite, not a 0/0.

    Parameters
    ----------
    s : array_like
        1/phi.
    y1, y2 : float
        Mole fraction of the species at the two film ends.

    Returns
    -------
    numpy.ndarray
        The value of phi * ln[(phi - y2)/(phi - y1)].
    """
    s = np.asarray(s, float)
    small = np.abs(s) < 1e-7
    ss = np.where(small, 1.0, s)
    raw = ss * (y1 - y2) / (1.0 - ss * y1)
    arg = np.where(small | ~np.isfinite(raw) | (raw <= -1.0), 0.0, raw)
    big = np.where(np.isfinite(raw) & (raw > -1.0), np.log1p(arg) / ss, np.nan)
    tay = (y1 - y2) * (1.0 + np.where(small, s, 0.0) * (y1 + y2) / 2.0)
    return np.where(small, tay, big)


R_RANGE = (1e-7, 1e5)   # |N_B/N_A| scanned by film_roots, either sign


def film_roots(Kfun, y1, y2, n_scan=N_SCAN):
    """Solve eq. (14) for A and for B simultaneously.

    Both equations say N_t = K_i * ln[(phi_i - y_i2)/(phi_i - y_i1)], and
    phi_A + phi_B = 1 because N_C = 0, so everything is a function of the single
    ratio r = N_B/N_A.  Every sign change on 1e-7 <= |r| <= 1e5 (R_RANGE) is
    bracketed by find_roots, so multiple solutions are found rather than missed.

    Parameters
    ----------
    Kfun : callable
        r -> (K_A, K_B), the mass-transfer coefficients c D / x per species.
    y1, y2 : array_like
        Mole fractions at the two film ends (3,).
    n_scan : int, optional
        Scan points per sign of r.

    Returns
    -------
    list of tuple
        One (N_A, N_B) pair per root.
    """
    yA1, yB1, yA2, yB2 = y1[0], y1[1], y2[0], y2[1]

    def H(r):
        r = np.atleast_1d(np.asarray(r, float))
        KA, KB = Kfun(r)
        with np.errstate(divide="ignore", invalid="ignore"):
            h = r * KA * g(1.0 + r, yA1, yA2) - KB * g(1.0 + 1.0 / r, yB1, yB2)
        return float(np.ravel(h)[0])

    out = []
    for r in roots_both_signs(H, *R_RANGE, n_scan=n_scan):
        KA, _ = Kfun(np.array([r]))
        NA = float(np.atleast_1d(KA)[0] * np.atleast_1d(g(np.array([1.0 + r]), yA1, yA2))[0])
        out.append((NA, r * NA))
    return out


def method1(e, rule=wilke_rule, comp="mean", x=XLEN, binaries=None, c=C_TOT):
    """Wilke's rule frozen at one composition, then eq. (14) for both species.

    Parameters
    ----------
    e : str
        Example label of Table 1.
    rule : callable, optional
        Mixture rule (`wilke_rule` or `arith_rule`).
    comp : {'mean', 'point1', 'point2'}, optional
        Composition at which the rule is frozen.
    x : float, optional
        Film length (cm).
    binaries : tuple of float, optional
        Replacement (D_AB, D_AC, D_BC) in cm2/s.
    c : float, optional
        Total molar concentration P/RT (mol/cm3).

    Returns
    -------
    tuple
        (roots, (D_A, D_B)): the (N_A, N_B) pairs and the effective
        diffusivities used.
    """
    _, y1, y2, DAB, DAC, DBC = unpack(e)
    if binaries is not None:
        DAB, DAC, DBC = binaries
    D = Dmat(DAB, DAC, DBC)
    yref = {"mean": 0.5 * (y1 + y2), "point1": y1, "point2": y2}[comp]
    DA, DB = rule(yref, 0, D), rule(yref, 1, D)
    KA, KB = c * DA / x, c * DB / x

    def K(r):
        r = np.asarray(r, float)
        return np.full(r.shape, KA), np.full(r.shape, KB)
    return film_roots(K, y1, y2), (DA, DB)


def m2_pieces(yref, phiA, D, clip=True):
    """Eqs. (27)-(38): effective film composition and D'' for A and for B.

    Parameters
    ----------
    yref : array_like
        Reference mole fractions (3,).
    phiA : float
        N_A / (N_A + N_B).
    D : numpy.ndarray
        Binary diffusivity matrix from `Dmat` (cm2/s).
    clip : bool, optional
        Clip negative partial pressures to zero, as the paper states.

    Returns
    -------
    dict
        Effective diffusivities under keys 0 and 1, and the partial
        pressures, mole fractions and negativity flags under
        ('pf', i), ('ypp', i), ('neg', i).
    """
    phi = (phiA, 1.0 - phiA)
    out = {}
    for i in (0, 1):
        others = [1, 2] if i == 0 else [0, 2]
        pf, neg = {}, {}
        for m in others:
            phi_im = 1.0 if m == 2 else phi[i]        # phi_iC = N_i/(N_i+N_C) = 1
            ypr = yref[i] / (yref[i] + yref[m])
            pf[m] = P_MM / phi_im * (yref[i] + yref[m]) * (phi_im - ypr)
            neg[m] = bool(pf[m] < 0)
        if clip:                                      # the paper's stated branch
            pf = {m: max(pf[m], 0.0) for m in pf}
        tot = sum(pf.values())
        ypp = {m: (pf[m] / tot if tot != 0 else np.nan) for m in pf}
        out[i] = 1.0 / sum(ypp[m] / D[i, m] for m in others)
        out[("pf", i)], out[("ypp", i)], out[("neg", i)] = pf, ypp, neg
    return out


def method2(e, clip=True, x=XLEN, binaries=None, c=C_TOT):
    """Wilke's Method 2: the rule at a flux-weighted film composition.

    Parameters
    ----------
    e : str
        Example label of Table 1.
    clip : bool, optional
        Clip negative partial pressures to zero, as the paper states.
    x : float, optional
        Film length (cm).
    binaries : tuple of float, optional
        Replacement (D_AB, D_AC, D_BC) in cm2/s.
    c : float, optional
        Total molar concentration P/RT (mol/cm3).

    Returns
    -------
    tuple
        (roots, info): the (N_A, N_B) pairs and `m2_pieces` for each.
    """
    _, y1, y2, DAB, DAC, DBC = unpack(e)
    if binaries is not None:
        DAB, DAC, DBC = binaries
    D = Dmat(DAB, DAC, DBC)
    yref = 0.5 * (y1 + y2)

    def K(r):
        r = np.atleast_1d(np.asarray(r, float))
        KA, KB = np.empty_like(r), np.empty_like(r)
        for k, rk in enumerate(r):
            p = m2_pieces(yref, 1.0 / (1.0 + rk), D, clip)
            KA[k], KB[k] = c * p[0] / x, c * p[1] / x
        return KA, KB
    roots = film_roots(K, y1, y2)
    info = [m2_pieces(yref, NA / (NA + NB), D, clip) for NA, NB in roots]
    return roots, info
