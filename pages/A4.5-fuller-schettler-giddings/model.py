"""Fuller-Schettler-Giddings binary gas diffusivity and a tube it is used in (page A4.5).

Eq. (4) of Fuller, Schettler and Giddings (1966) with every published constant in
one parameter dictionary, the least-squares problem of their eq. (3) for refitting
the Table I volumes, and steady diffusion of a dilute species through a stagnant
gas in a tube with a linear wall temperature, solved with pymrm and checked
against an adaptive quadrature of 1/Gamma.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy.sparse.linalg import spsolve
from pymrm import (construct_grad, construct_div, construct_coefficient_matrix,
                   compute_boundary_values)

R_GAS = 82.0574          # cm3 atm / (mol K)

# Formulae and AROMATIC ring counts of the 42 organic vapours of Table III.
ORG = {
    "METHANE": (dict(C=1, H=4), 0), "ETHANE": (dict(C=2, H=6), 0),
    "PROPANE": (dict(C=3, H=8), 0), "N-BUTANE": (dict(C=4, H=10), 0),
    "ISOBUTANE": (dict(C=4, H=10), 0), "N-HEXANE": (dict(C=6, H=14), 0),
    "2,3-DIMETHYLBUTANE": (dict(C=6, H=14), 0),
    "N-HEPTANE": (dict(C=7, H=16), 0), "2,4-DIMETHYLPENTANE": (dict(C=7, H=16), 0),
    "N-OCTANE": (dict(C=8, H=18), 0), "2,2,4-TRIMETHYLPENTANE": (dict(C=8, H=18), 0),
    "N-DECANE": (dict(C=10, H=22), 0), "2,3,3-TRIMETHYLHEPTANE": (dict(C=10, H=22), 0),
    "N-DODECANE": (dict(C=12, H=26), 0), "CYCLOHEXANE": (dict(C=6, H=12), 0),
    "METHYLCYCLOPENTANE": (dict(C=6, H=12), 0), "ETHYLENE": (dict(C=2, H=4), 0),
    "BENZENE": (dict(C=6, H=6), 1), "TOLUENE": (dict(C=7, H=8), 1),
    "CHLOROBENZENE": (dict(C=6, H=5, Cl=1), 1),
    "NITROBENZENE": (dict(C=6, H=5, N=1, O=2), 1),
    "ANILINE": (dict(C=6, H=7, N=1), 1), "DIPHENYL": (dict(C=12, H=10), 2),
    "PYRIDINE": (dict(C=5, H=5, N=1), 1), "PIPERIDINE": (dict(C=5, H=11, N=1), 0),
    "THIOPHENE": (dict(C=4, H=4, S=1), 1), "TETRAHYDROTHIOPHENE": (dict(C=4, H=8, S=1), 0),
    "METHANOL": (dict(C=1, H=4, O=1), 0), "ETHANOL": (dict(C=2, H=6, O=1), 0),
    "PROPANOL": (dict(C=3, H=8, O=1), 0), "2-PROPANOL": (dict(C=3, H=8, O=1), 0),
    "BUTANOL": (dict(C=4, H=10, O=1), 0), "2-BUTANOL": (dict(C=4, H=10, O=1), 0),
    "PENTANOL": (dict(C=5, H=12, O=1), 0), "2-PENTANOL": (dict(C=5, H=12, O=1), 0),
    "HEXANOL": (dict(C=6, H=14, O=1), 0), "ETHYLACETATE": (dict(C=4, H=8, O=2), 0),
    "ETHYLENEOXIDE": (dict(C=2, H=4, O=1), 0), "HCN": (dict(C=1, H=1, N=1), 0),
    "CYANOGENCHLORIDE": (dict(C=1, N=1, Cl=1), 0),
    "PHOSGENE": (dict(C=1, O=1, Cl=2), 0),
    "CHLOROPICRIN": (dict(C=1, N=1, O=2, Cl=3), 0),
}
# Atomic weights (IUPAC); the paper prints none.
AW = {"H": 1.008, "D": 2.014, "He": 4.0026, "C": 12.011, "N": 14.007,
      "O": 15.999, "F": 18.998, "S": 32.06, "Cl": 35.45, "Ne": 20.180,
      "Ar": 39.948, "Kr": 83.798, "Xe": 131.29, "Br": 79.904}
MW_SIMPLE = {
    "H2": 2 * AW["H"], "D2": 4.00, "HE": AW["He"], "N2": 2 * AW["N"], "O2": 2 * AW["O"],
    "AIR": 28.96, "NE": AW["Ne"], "AR": AW["Ar"], "KR": AW["Kr"], "XE": AW["Xe"],
    "CO": AW["C"] + AW["O"], "CO2": AW["C"] + 2 * AW["O"], "N2O": 2 * AW["N"] + AW["O"],
    "NH3": AW["N"] + 3 * AW["H"], "WATER": 2 * AW["H"] + AW["O"],
    "FREON12": AW["C"] + 2 * AW["Cl"] + 2 * AW["F"], "SF6": AW["S"] + 6 * AW["F"],
    "CL2": 2 * AW["Cl"], "BR2": 2 * AW["Br"], "SO2": AW["S"] + 2 * AW["O"]}
ATOMS = ["C", "H", "O", "N", "Cl", "S"]
_RENAME = {"H2O": "WATER", "CCL2F2": "FREON12"}   # Table I names to Table III keys


def fsg_params(volumes):
    """Every published constant of eq. (4) in one dictionary.

    Parameters
    ----------
    volumes : pandas.DataFrame
        Table I with columns entry, diffusion_volume and kind.

    Returns
    -------
    dict
        ``C`` (1.00e-3), ``b`` (1.75), ``v_atom`` (atomic increments), ``v_ring``
        (aromatic ring increment), ``v_simple`` (simple-molecule volumes, Table III
        keys), ``mw_simple`` (molar masses, g/mol) and ``ring_extra`` (species given
        one ring term beyond their aromatic rings; empty as published).
    """
    v_atom = {r.entry: float(r.diffusion_volume) for r in volumes.itertuples()
              if r.kind == "atomic increment"}
    v_ring = float(volumes.loc[volumes.kind == "structural increment",
                               "diffusion_volume"].iloc[0])
    v_simple = {_RENAME.get(r.entry.upper(), r.entry.upper()): float(r.diffusion_volume)
                for r in volumes.itertuples() if r.kind == "simple molecule"}
    return dict(C=1.00e-3, b=1.75, v_atom=v_atom, v_ring=v_ring, v_simple=v_simple,
                mw_simple=dict(MW_SIMPLE), ring_extra=())


def sigma_v(name, par):
    """Diffusion volume of a species: Table I value or the sum of its increments."""
    if name in par["v_simple"]:
        return par["v_simple"][name]
    f, nr = ORG[name]
    nr = nr + (name in par["ring_extra"])
    return sum(par["v_atom"][a] * n for a, n in f.items()) + nr * par["v_ring"]


def mw(name, par):
    """Molar mass of a species in g/mol."""
    if name in par["mw_simple"]:
        return par["mw_simple"][name]
    return sum(AW[a] * n for a, n in ORG[name][0].items())


def D_fsg(A, B, T, par, p=1.0):
    """Eq. (4): binary diffusivity in cm2/s at T in K and p in atm.

    Parameters
    ----------
    A, B : str
        Species keys as in Table III.
    T : float or numpy.ndarray
        Temperature, K.
    par : dict
        Constants from :func:`fsg_params`.
    p : float, optional
        Pressure, atm.
    """
    vA, vB = sigma_v(A, par), sigma_v(B, par)
    return (par["C"] * T**par["b"] * (1.0 / mw(A, par) + 1.0 / mw(B, par))**0.5
            / (p * (vA**(1 / 3) + vB**(1 / 3))**2))


class VolumeFit:
    """The least-squares problem of eq. (3) for the 27 Table I volumes at fixed b.

    Parameters
    ----------
    comp_a, comp_b : sequence of str
        Species of each data row.
    T, D_obs : numpy.ndarray
        Temperature (K) and observed diffusivity (cm2/s) of each row.
    par : dict
        Constants from :func:`fsg_params`; its volumes are the published start.
    """

    def __init__(self, comp_a, comp_b, T, D_obs, par):
        self.par = par
        self.names = ATOMS + ["ring"] + list(par["v_simple"])
        idx = {p: i for i, p in enumerate(self.names)}

        def incidence(sp):
            v = np.zeros(len(self.names))
            if sp in par["v_simple"]:
                v[idx[sp]] = 1.0
                return v
            f, nr = ORG[sp]
            for a, n in f.items():
                v[idx[a]] += n
            v[idx["ring"]] += nr
            return v
        self.MA = np.array([incidence(s) for s in comp_a])
        self.MB = np.array([incidence(s) for s in comp_b])
        self.massfac = np.array([(1 / mw(a, par) + 1 / mw(b, par))**0.5
                                 for a, b in zip(comp_a, comp_b)])
        self.T, self.D_obs = np.asarray(T, float), np.asarray(D_obs, float)
        self.published = np.array([par["v_atom"][k] for k in ATOMS] + [par["v_ring"]]
                                  + list(par["v_simple"].values()))
        # lower bounds: volumes positive, the ring increment free to be negative
        self.lower = np.r_[np.full(len(ATOMS), 1e-3), -1e3,
                           np.full(len(par["v_simple"]), 1e-3)]

    def D(self, theta, b):
        """Eq. (4) on every row for volumes ``theta``."""
        va, vb = self.MA @ theta, self.MB @ theta
        return self.par["C"] * self.T**b * self.massfac / (va**(1 / 3) + vb**(1 / 3))**2

    def resid(self, theta, mask=None, b=None, objective="relative"):
        """Residuals of eq. (3), (E - Y)/Y, or E - Y for ``objective='absolute'``."""
        d = self.D(theta, self.par["b"] if b is None else b) - self.D_obs
        d = d / self.D_obs if objective == "relative" else d
        return d if mask is None else d[mask]

    def fit(self, start=None, mask=None, b=None, objective="relative", tol=1e-14,
            lower=None):
        """Bounded least squares from ``start`` (default: the published volumes).

        ``lower`` replaces the lower bounds (default: 1e-3 for every volume, -1e3
        for the ring increment).

        Returns
        -------
        scipy.optimize.OptimizeResult
            The ``least_squares`` result; ``x`` holds the fitted volumes.
        """
        x0 = self.published if start is None else start
        return least_squares(self.resid, x0, args=(mask, b, objective),
                             bounds=(self.lower if lower is None else lower, np.inf),
                             xtol=tol, ftol=tol, gtol=tol)


def gamma(T, D_ref, T_ref, b, p=1.0):
    """c_t D = (p/RT) D_ref (T/T_ref)^b in mol/(cm s)."""
    return (p / (R_GAS * T)) * D_ref * (T / T_ref)**b


class Tube:
    """Steady diffusion of a dilute species through stagnant gas in a heated tube.

    d/dz(-Gamma(z) dx/dz) = 0 on 0 < z < L, x(0) = x0, x(L) = xL, with
    T(z) = T0 + (T1 - T0) z/L and Gamma = gamma(T). The problem is linear in x:
    one sparse solve, no Newton iteration.

    Parameters
    ----------
    n : int
        Number of cells.
    L, T0, T1 : float
        Length (cm) and the end temperatures (K).
    D_ref, T_ref, b : float
        D = D_ref (T/T_ref)^b, cm2/s.
    p, x0, xL : float, optional
        Pressure (atm) and the end mole fractions.
    nu : int, optional
        Geometry index of construct_div; 0 is the tube's constant cross-section.
    face_rule : {"face", "left_cell"}, optional
        Gamma at a face: ``"face"`` evaluates gamma at the face temperature, the
        one rule for every face, end faces included; ``"left_cell"`` takes the
        value of the cell on the low-z side (the first cell at z = 0).
    bc_right : dict, optional
        Replaces the Dirichlet condition at z = L.
    """

    def __init__(self, n, L, T0, T1, D_ref, T_ref, b, p=1.0, x0=0.10, xL=0.0,
                 nu=0, face_rule="face", bc_right=None):
        self.z_f = np.linspace(0.0, L, n + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        T_f = T0 + (T1 - T0) * self.z_f / L
        if face_rule == "face":
            self.G_f = gamma(T_f, D_ref, T_ref, b, p)
        elif face_rule == "left_cell":
            G_c = gamma(T0 + (T1 - T0) * self.z_c / L, D_ref, T_ref, b, p)
            self.G_f = np.r_[G_c[0], G_c]
        else:
            raise ValueError(f"unknown face_rule {face_rule!r}")
        # a dx/dn + b x = d on the OUTWARD normal n:
        #   z = 0 (n = -z): x = x0  ->  a = 0, b = 1, d = x0
        #   z = L (n = +z): x = xL  ->  a = 0, b = 1, d = xL
        # with a = 0 the sign of n does not enter; a flux condition would need it.
        self.bc = ({"a": 0.0, "b": 1.0, "d": x0},
                   {"a": 0.0, "b": 1.0, "d": xL} if bc_right is None else bc_right)
        grad, grad_bc = construct_grad((n, 1), self.z_f, self.z_c, self.bc)
        div = construct_div((n, 1), self.z_f, nu=nu)
        G = construct_coefficient_matrix(self.G_f.reshape(-1, 1))
        self.flux_op, self.flux_bc = -(G @ grad), -(G @ grad_bc)   # face flux = op x + bc
        self.A = div @ self.flux_op                       # div(flux) = 0
        self.rhs = -np.asarray((div @ self.flux_bc).todense()).ravel()

    def solve(self, permc_spec="COLAMD"):
        """Solve the linear system once.

        Returns
        -------
        dict
            ``x`` (cell mole fractions), ``flux`` (molar flux into the tube at
            z = 0 from the boundary gradient, mol/(cm2 s)) and ``spread`` (max
            minus min of the flux over all faces).
        """
        x = spsolve(self.A, self.rhs, permc_spec=permc_spec).reshape(-1, 1)
        _, dxdz0, _, _ = compute_boundary_values(x, self.z_f, self.z_c, self.bc)
        face = (self.flux_op @ x).ravel() + np.asarray(self.flux_bc.todense()).ravel()
        flux = -self.G_f[0] * float(np.ravel(dxdz0)[0])       # gradient along +z
        return dict(x=x.ravel(), flux=flux, spread=float(face.max() - face.min()))


def flux_quadrature(L, T0, T1, D_ref, T_ref, b, p=1.0, x0=0.10, xL=0.0, epsrel=1e-12):
    """Flux (x0 - xL) / integral of dz/Gamma by adaptive quadrature; no grid, no pymrm.

    Returns
    -------
    tuple of float
        The flux and quad's error estimate relative to the integral.
    """
    f = lambda z: 1.0 / gamma(T0 + (T1 - T0) * z / L, D_ref, T_ref, b, p)
    integral, err = quad(f, 0.0, L, epsabs=0.0, epsrel=epsrel, limit=400)
    return (x0 - xL) / integral, err / integral
