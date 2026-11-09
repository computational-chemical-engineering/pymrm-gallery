"""Yagi and Kunii (1957) effective conductivity of a packed bed, and the annulus
apparatus it was measured in.

Equations are numbered as in the paper. Temperatures are in degC, lengths in
metres except where a docstring says Dp in mm, conductivities in
kcal/(m hr degC). No plotting and no data loading: the page passes the data
tables in (``set_air_table``, ``set_iron_sphere_row``) and the Table rows as
objects with the attribute names of the data files.
"""
import numpy as np
from scipy.sparse.linalg import spsolve

from pymrm import (NumJac, compute_boundary_values, construct_coefficient_matrix,
                   construct_div, construct_grad, newton)
from pymrm.checks import find_roots, residual_check

SIG4 = 0.1952   # = 4*sigma in the paper's units: kcal/(m2 hr degC) per (T/100 K)^3

# Eq. (19)'s printed apparatus geometry, metres
L_BED, D0_M, DI_M = 0.25, 0.146, 0.0375
R_I, R_O = DI_M / 2, D0_M / 2

_KGT = None            # air conductivity table, set by set_air_table
_IRON_SPHERE_ROW = None  # Table 4 iron-sphere packing, set by set_iron_sphere_row


def set_air_table(kgt):
    """Register the air conductivity table used by `kg_air`.

    Parameters
    ----------
    kgt : pandas.DataFrame
        Columns ``temp_C`` and ``kg_kcal_m_hr_C`` (the paper's printed table).
    """
    global _KGT
    _KGT = kgt


def set_iron_sphere_row(row):
    """Register the Table 4 iron-sphere packing used by `ke0_iron`.

    Parameters
    ----------
    row : pandas.Series
        Table 4 row with ``eps``, ``beta``, ``gamma``, ``phi``, ``Dp_mm`` and
        ``ks_kcal_m_hr_C``.
    """
    global _IRON_SPHERE_ROW
    _IRON_SPHERE_ROW = row


def h_rs(t_C, p=1.0):
    """Eq. (7): radiation coefficient, solid surface to solid surface.

    Parameters
    ----------
    t_C : float or ndarray
        Temperature, degC.
    p : float, optional
        Emissivity-like factor p of the paper; 1 is a black body.

    Returns
    -------
    float or ndarray
        h_rs in kcal/(m2 hr degC).
    """
    return SIG4 * (p / (2.0 - p)) * ((t_C + 273.0) / 100.0) ** 3


def h_rv(t_C, eps, p=1.0):
    """Eq. (8): radiation coefficient, void to void.

    Parameters
    ----------
    t_C : float or ndarray
        Temperature, degC.
    eps : float
        Bed voidage.
    p : float, optional
        Emissivity-like factor p of the paper; 1 is a black body.

    Returns
    -------
    float or ndarray
        h_rv in kcal/(m2 hr degC).
    """
    return SIG4 / (1.0 + eps / (2.0 * (1.0 - eps)) * (1.0 - p) / p) \
        * ((t_C + 273.0) / 100.0) ** 3


def eq15(eps, beta, gamma, phi, Dp_m, ks, kg, hrs, hrv):
    """Eq. (15): ke0/kg for gas-filled voids. Dp in METRES (h_r is per m2).

    Parameters
    ----------
    eps : float
        Bed voidage.
    beta : float
        Ratio of the effective path length to the particle diameter.
    gamma : float
        Ratio of the solid-path length to the particle diameter.
    phi : float
        Effective gas-film thickness at the contacts, as a fraction of Dp.
    Dp_m : float
        Particle diameter, m.
    ks, kg : float
        Solid and gas conductivities, kcal/(m hr degC).
    hrs, hrv : float
        Radiation coefficients from `h_rs` and `h_rv`, kcal/(m2 hr degC).

    Returns
    -------
    float
        Stagnant effective conductivity over the gas conductivity, ke0/kg.
    """
    contact = 1.0 / ((1.0 / phi) + Dp_m * hrs / kg)     # film + solid-solid radiation
    solid_path = beta * (1.0 - eps) / (gamma * (kg / ks) + contact)
    void_rad = eps * beta * Dp_m * hrv / kg
    return solid_path + void_rad


def eq16(eps, beta, phi, ks, kg):
    """Eq. (16): the radiation-free form, for fine particles or liquid-filled voids.

    Parameters
    ----------
    eps : float
        Bed voidage.
    beta : float
        Ratio of the effective path length to the particle diameter.
    phi : float
        Effective gas-film thickness at the contacts, as a fraction of Dp.
    ks, kg : float
        Solid and gas conductivities, kcal/(m hr degC).

    Returns
    -------
    float
        ke0/kg.
    """
    return beta * (1.0 - eps) / ((kg / ks) + phi)


def eq14(eps, beta, Dp_m, ks, hrs, hrv, delta=0.0):
    """Eq. (14): the vacuum limit, normalised on ks rather than kg.

    Parameters
    ----------
    eps : float
        Bed voidage.
    beta : float
        Ratio of the effective path length to the particle diameter.
    Dp_m : float
        Particle diameter, m.
    ks : float
        Solid conductivity, kcal/(m hr degC).
    hrs, hrv : float
        Radiation coefficients, kcal/(m2 hr degC).
    delta : float, optional
        Fraction of the cross-section occupied by the contact conduction
        path (0 in the limit).

    Returns
    -------
    float
        ke0/ks.
    """
    return delta + beta * (1.0 - eps) / (1.0 + ks / (Dp_m * hrs)) \
        + eps * beta * Dp_m * hrv / ks


def eq6(ke0_over_kg, alpha_beta, N_PeM):
    """Eq. (6): the static term plus Ranz's lateral mixing.

    Parameters
    ----------
    ke0_over_kg : float
        Stagnant conductivity over the gas conductivity.
    alpha_beta : float
        Lateral-mixing coefficient (alpha times beta).
    N_PeM : float
        Modified Peclet number of the flow.

    Returns
    -------
    float
        ke/kg with flow.
    """
    return ke0_over_kg + alpha_beta * N_PeM


def kg_air(t_C):
    """Air conductivity, kcal/(m hr degC), from the paper's own printed table.

    Parameters
    ----------
    t_C : float or ndarray
        Temperature, degC. `set_air_table` must have been called.

    Returns
    -------
    float or ndarray
        Linear interpolation between the printed nodes (0 to 1200 degC every
        200 degC). Outside that range ``np.interp`` holds the end value: above
        1200 degC kg stays at the 1200 degC entry.
    """
    return np.interp(t_C, _KGT.temp_C.values, _KGT.kg_kcal_m_hr_C.values)


def one_root(f, lo, hi, **kw):
    """The single root of ``f`` on ``[lo, hi]``, by `pymrm.checks.find_roots`.

    Parameters
    ----------
    f : callable
        Scalar function of one variable.
    lo, hi : float
        Physically possible range; the scan covers all of it.
    **kw
        Passed to `find_roots` (``log``, ``n_scan``, ``xtol``).

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        When the scan finds no root or more than one.
    """
    found = find_roots(f, lo, hi, **kw)
    if len(found["roots"]) != 1:
        raise ValueError(f"expected one root on [{lo:g}, {hi:g}], found "
                         f"{len(found['roots'])}: {found['message']}")
    return float(found["roots"][0])


def annulus_operators(n, T_in, T_out, nu=1):
    """Grid, boundary conditions and constant operators of the annulus.

    Parameters
    ----------
    n : int
        Number of equal cells between ``R_I`` and ``R_O``.
    T_in, T_out : float
        Inner and outer wall temperatures, degC.
    nu : int, optional
        Geometry: 1 = cylindrical (face area ~ r), 0 = Cartesian.

    Returns
    -------
    dict
        ``r_f``, ``r_c``, ``bc``, ``grad`` (sparse, faces x cells), ``grad_bc``
        (dense, faces x 1), ``div`` (sparse, cells x faces) and ``shape``.
    """
    r_f = np.linspace(R_I, R_O, n + 1)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n, 1)                                   # spatial axis first, one field
    # a dT/dn + b T = d with n the OUTWARD normal; a = 0 makes both Dirichlet:
    #   inner wall r = r_i:  T = T_in
    #   outer wall r = r_o:  T = T_out
    bc = ({"a": 0.0, "b": 1.0, "d": T_in}, {"a": 0.0, "b": 1.0, "d": T_out})
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc, axis=0)
    div = construct_div(shape, r_f, nu=nu, axis=0)   # nu=1: (1/r) d(r q)/dr
    return dict(r_f=r_f, r_c=r_c, bc=bc, grad=grad.tocsr(), shape=shape,
                grad_bc=np.asarray(grad_bc.todense()).reshape(-1, 1), div=div.tocsr())


def face_temperatures(T, ops):
    """Temperatures on every face: interior faces the mean of the two cells
    (the midpoint of equal cells), wall faces from `compute_boundary_values`
    under the solve's own boundary conditions.

    Parameters
    ----------
    T : ndarray
        Cell temperatures, shape (n, 1).
    ops : dict
        From `annulus_operators`.

    Returns
    -------
    ndarray
        Face temperatures, shape (n + 1, 1).
    """
    T = T.reshape(-1, 1)
    t_in, _, t_out, _ = compute_boundary_values(T, ops["r_f"], ops["r_c"], ops["bc"], axis=0)
    Tf = np.empty((T.shape[0] + 1, 1))
    Tf[1:-1] = 0.5 * (T[:-1] + T[1:])
    Tf[0], Tf[-1] = np.ravel(t_in)[0], np.ravel(t_out)[0]
    return Tf


def annulus_flux(T, ops, kfun=None, k_const=None, k_f=None):
    """Radial heat flux on every face, kcal/(m2 hr): q = -k(T_face) dT/dr.

    The gradient is the operator's, wall rows included; the face
    conductivity is ``k_f`` when given, else ``k_const``, else ``kfun`` at
    `face_temperatures`.
    """
    T = T.reshape(-1, 1)
    if k_f is None:
        k_f = (np.full((T.shape[0] + 1, 1), float(k_const)) if k_const is not None
               else np.asarray(kfun(face_temperatures(T, ops)), float).reshape(-1, 1))
    return -k_f * (ops["grad"] @ T + ops["grad_bc"])


def solve_annulus(n, T_in, T_out, kfun=None, k_const=None, nu=1, method="newton",
                  tol=1e-10, maxit=300, info=None):
    """Steady radial conduction in the packed annulus, d(r q)/dr = 0.

    Finite volumes with pymrm operators, Dirichlet at both walls. A constant
    conductivity makes the problem linear: one sparse solve. A conductivity
    ``kfun(T)`` makes it nonlinear; ``method="newton"`` solves it with
    `pymrm.newton` on the scaled unknown ``(T - T_out)/(T_in - T_out)`` and a
    tridiagonal ``NumJac((n, 1))``; ``method="picard"`` is the second route:
    freeze the face conductivities, solve the linear problem, repeat until the
    largest temperature change of a sweep is below ``tol`` degC.

    Parameters
    ----------
    n : int
        Number of cells between ``R_I`` and ``R_O``.
    T_in, T_out : float
        Inner and outer wall temperatures, degC.
    kfun : callable, optional
        Conductivity as a function of temperature; used when ``k_const`` is None.
    k_const : float, optional
        Constant conductivity.
    nu : int, optional
        Geometry: 1 = cylindrical, 0 = Cartesian. An argument so a break row can
        change it.
    method : {"newton", "picard"}, optional
        Route for a temperature-dependent conductivity.
    tol : float, optional
        Newton: step tolerance on the scaled unknown. Picard: largest change of
        a sweep, degC.
    maxit : int, optional
        Picard: maximum number of sweeps.
    info : dict, optional
        Filled with ``iterations`` and, for Picard, ``steps`` (largest change of
        every sweep).

    Returns
    -------
    r_f : ndarray
        Face radii, m.
    r_c : ndarray
        Cell-centre radii, m.
    T : ndarray
        Cell temperatures, degC.
    Qprime : ndarray
        Heat rate per unit bed length through each face, kcal/(m hr):
        ``2 pi r_f q`` (the 2 pi is the angle of a full annulus).
    """
    ops = annulus_operators(n, T_in, T_out, nu)
    info = {} if info is None else info
    if k_const is not None:
        A = (ops["div"] @ (-float(k_const) * ops["grad"])).tocsc()
        b = ops["div"] @ (-float(k_const) * ops["grad_bc"])
        T = spsolve(A, -b.ravel()).reshape(-1, 1)
        info["iterations"] = 1
    elif method == "newton":
        dT = T_in - T_out
        numjac = NumJac(ops["shape"], axes_diagonals=[0])     # 3-point stencil

        def heat(theta):
            return ops["div"] @ annulus_flux(T_out + dT * theta, ops, kfun=kfun)

        def residual(theta):
            g, jac = numjac(heat, theta.reshape(ops["shape"]))
            return g.reshape(-1, 1), jac

        theta0 = np.linspace(1.0, 0.0, n).reshape(ops["shape"])
        res = newton(residual, theta0, tol=tol, stagnation=True)
        check = residual_check(residual, res.x)
        if not (res.success and check["ok"]):
            raise RuntimeError(f"annulus Newton failed: {res.message}, backward error "
                               f"{check['backward_error']:.1e}")
        T = T_out + dT * res.x.reshape(-1, 1)
        info["iterations"] = int(res.nit) if hasattr(res, "nit") else None
    elif method == "picard":
        T = np.linspace(T_in, T_out, n).reshape(-1, 1)
        steps = []
        for _ in range(maxit):
            k_f = np.asarray(kfun(face_temperatures(T, ops)), float).reshape(-1, 1)
            Km = construct_coefficient_matrix(k_f, shape=ops["shape"], axis=0)
            A = (ops["div"] @ (Km @ ops["grad"])).tocsc()
            b0 = ops["div"] @ (k_f * ops["grad_bc"])
            T_new = spsolve(A, -b0.ravel()).reshape(-1, 1)
            steps.append(float(np.max(np.abs(T_new - T))))
            T = T_new
            if steps[-1] < tol:
                break
        info["steps"] = steps
        info["iterations"] = len(steps)
        if steps[-1] >= tol:
            raise RuntimeError(f"annulus fixed point: last step {steps[-1]:.1e} degC is not "
                               f"below {tol:g} after {maxit} sweeps")
    else:
        raise ValueError(f"unknown method {method!r}")
    q = annulus_flux(T, ops, kfun=kfun, k_const=k_const)
    return ops["r_f"], ops["r_c"], T.ravel(), 2.0 * np.pi * ops["r_f"] * q.ravel()


def invert_phi(row, p=1.0, drop_hrs=False, drop_hrv=False, gamma=1.0):
    """phi that makes Eq. (15)/(16) return the paper's printed Calcd. value.

    Parameters
    ----------
    row : pandas.Series
        Table 3 row with ``eps``, ``ks_kcal_m_hr_C``, ``kg``,
        ``calcd_ke0_over_kg``, ``Dp_mm`` and ``temp_C``.
    p : float, optional
        Emissivity-like factor p of the paper.
    drop_hrs, drop_hrv : bool, optional
        Delete the solid-solid or void-void radiation term (the deletion test).
    gamma : float, optional
        Value of gamma used in Eq. (15).

    Returns
    -------
    float
        The film thickness phi, or nan when the row has no printed kg or no
        positive solution.
    """
    eps, ks, kg, tgt = row.eps, row.ks_kcal_m_hr_C, row.kg, row.calcd_ke0_over_kg
    if np.isnan(kg):
        return np.nan                     # He / H2 rows: the paper prints no kg
    if np.isnan(row.Dp_mm):               # no Dp -> Eq. (16), radiation cannot enter
        return (1.0 - eps) / tgt - kg / ks
    Dp = row.Dp_mm * 1e-3
    hs = 0.0 if drop_hrs else h_rs(row.temp_C, p)
    hv = 0.0 if drop_hrv else h_rv(row.temp_C, eps, p)
    resid = tgt - eps * Dp * hv / kg
    if resid <= 0:
        return np.nan
    c = (1.0 - eps) / resid - gamma * kg / ks
    if c <= 0 or (1.0 / c - Dp * hs / kg) <= 0:
        return np.nan
    return 1.0 / (1.0 / c - Dp * hs / kg)


def forward(row, phi, p=1.0, drop_hrs=False, drop_hrv=False, gamma=1.0, beta=1.0):
    """Eq. (15)/(16) evaluated for one Table 3 row at a given phi.

    Parameters
    ----------
    row : pandas.Series
        Table 3 row with ``eps``, ``ks_kcal_m_hr_C``, ``kg``, ``Dp_mm`` and
        ``temp_C``.
    phi : float
        Effective gas-film thickness, as a fraction of Dp.
    p : float, optional
        Emissivity-like factor p of the paper.
    drop_hrs, drop_hrv : bool, optional
        Delete the solid-solid or void-void radiation term.
    gamma, beta : float, optional
        Values of gamma and beta used in Eq. (15).

    Returns
    -------
    float
        ke0/kg; Eq. (16) when the row has no particle diameter.
    """
    eps, ks, kg = row.eps, row.ks_kcal_m_hr_C, row.kg
    if np.isnan(row.Dp_mm):
        return eq16(eps, beta, phi, ks, kg)
    Dp = row.Dp_mm * 1e-3
    hs = 0.0 if drop_hrs else h_rs(row.temp_C, p)
    hv = 0.0 if drop_hrv else h_rv(row.temp_C, eps, p)
    return eq15(eps, beta, gamma, phi, Dp, ks, kg, hs, hv)


def phi_needed(row, phi_lo=1e-8, phi_hi=10.0, **kw):
    """The phi that reproduces the printed Calcd. value, or nan if none in the bracket.

    Eq. (15) is monotone *decreasing* in phi -- a thicker gas film in the void is more
    resistance, not less -- so a scan of the whole bracket finds at most one root;
    more than one raises.

    Parameters
    ----------
    row : pandas.Series
        Table 3 row (see `forward`).
    phi_lo, phi_hi : float, optional
        Bracket for the root find.
    **kw
        Passed to `forward` (``p``, ``drop_hrs``, ``drop_hrv``, ``gamma``, ``beta``).

    Returns
    -------
    float
        The single root found by `pymrm.checks.find_roots`, or nan when the scan finds none.
    """
    f = lambda ph: forward(row, ph, **kw) - row.calcd_ke0_over_kg
    found = find_roots(f, phi_lo, phi_hi, log=True, xtol=1e-14)
    if len(found["roots"]) > 1:
        raise ValueError(f"more than one phi reproduces the row: {found['roots']}")
    return float(found["roots"][0]) if found["roots"] else np.nan


def can_match(row, phi_lo=1e-4, phi_hi=1.0, **kw):
    """Shortfall (%) in reproducing the printed Calcd. with phi anywhere in [phi_lo, phi_hi].

    The shortfall is measured from whichever END OF THE BRACKET is nearer the target.
    Taking max() instead reports the far end whenever
    the target is unreachable from one side, because Eq. (15) decreases with phi.

    Parameters
    ----------
    row : pandas.Series
        Table 3 row (see `forward`).
    phi_lo, phi_hi : float, optional
        Bracket for phi.
    **kw
        Passed to `forward`.

    Returns
    -------
    float
        Zero when the target lies inside the bracket's range, else the distance
        in per cent from the nearer end.
    """
    tgt = row.calcd_ke0_over_kg
    ends = [forward(row, phi_lo, **kw), forward(row, phi_hi, **kw)]
    if (ends[0] - tgt) * (ends[1] - tgt) <= 0:
        return 0.0
    return 100 * min(abs(v / tgt - 1) for v in ends)


def sweep(row, T, p=1.0, beta=None, gamma=None, phi=None, drop_rad=False, hr_scale=1.0):
    """ke0/kg and its radiation-free counterpart for one Table 4 packing.

    `hr_scale` exists only so the break table can attack the 0.1952 prefactor.

    Parameters
    ----------
    row : pandas.Series
        Table 4 row with ``Dp_mm``, ``eps``, ``beta``, ``gamma``, ``phi`` and
        ``ks_kcal_m_hr_C``.
    T : float or ndarray
        Temperature, degC.
    p : float, optional
        Emissivity-like factor p of the paper.
    beta, gamma, phi : float, optional
        Overrides of the row's values.
    drop_rad : bool, optional
        Delete both radiation terms.
    hr_scale : float, optional
        Factor on both radiation coefficients.

    Returns
    -------
    full, norad : float or ndarray
        ke0/kg with and without radiation.
    radiation_pct : float or ndarray
        Share of ke0 carried by radiation, per cent.
    """
    kg = kg_air(T)
    Dp = row.Dp_mm * 1e-3
    b = row.beta if beta is None else beta
    g = row.gamma if gamma is None else gamma
    ph = row.phi if phi is None else phi
    hs = 0.0 if drop_rad else hr_scale * h_rs(T, p)
    hv = 0.0 if drop_rad else hr_scale * h_rv(T, row.eps, p)
    full = eq15(row.eps, b, g, ph, Dp, row.ks_kcal_m_hr_C, kg, hs, hv)
    norad = eq15(row.eps, b, g, ph, Dp, row.ks_kcal_m_hr_C, kg, 0.0, 0.0)
    return full, norad, 100 * (1 - norad / full)


def eq15_row(row, T, ks=None, beta=None, p=1.0):
    """ke0/kg for one Table 4 packing, with optional overrides.

    Parameters
    ----------
    row : pandas.Series
        Table 4 row (see `sweep`).
    T : float or ndarray
        Temperature, degC.
    ks, beta : float, optional
        Overrides of the row's solid conductivity and beta.
    p : float, optional
        Emissivity-like factor p of the paper.

    Returns
    -------
    float or ndarray
        ke0/kg from Eq. (15).
    """
    kg = kg_air(T)
    return eq15(row.eps, row.beta if beta is None else beta, row.gamma, row.phi,
                row.Dp_mm * 1e-3, row.ks_kcal_m_hr_C if ks is None else ks, kg,
                h_rs(T, p), h_rv(T, row.eps, p))


def beta_of(r):
    """The figure's own beta where it labels one (Figure 13), else Table 4's.

    Parameters
    ----------
    r : pandas.Series
        Curve-end row with a ``beta_curve`` entry (nan when the figure does not
        label one).

    Returns
    -------
    float or None
        The labelled beta, or None to fall back on Table 4's value.
    """
    return float(r.beta_curve) if np.isfinite(r.beta_curve) else None


def eq13_from_network(eps, beta, gamma, phi, Dp, ks, kg, hrs, hrv, delta, dx=1.0):
    """ke0/kg rebuilt from Eqs. (9)-(12) rather than read off Eq. (13).

    Parameters
    ----------
    eps, beta, gamma, phi : float or ndarray
        As in `eq15`.
    Dp : float or ndarray
        Particle diameter, m.
    ks, kg : float or ndarray
        Solid and gas conductivities, kcal/(m hr degC).
    hrs, hrv : float or ndarray
        Radiation coefficients, kcal/(m2 hr degC).
    delta : float or ndarray
        Fraction of the cross-section occupied by the contact conduction path.
    dx : float, optional
        Slab thickness; drops out of the result.

    Returns
    -------
    float or ndarray
        ke0/kg.
    """
    NDx = dx / (beta * Dp)                                     # Eq. (3): N = 1/(beta Dp)
    inv_Us = NDx * (gamma * Dp / ks + 1.0 / (kg / (phi * Dp) + hrs))   # Eq. (10)
    inv_Uv = NDx / hrv                                                  # Eq. (11)
    q_over_dt = delta * ks / dx + (1.0 - eps - delta) / inv_Us + eps / inv_Uv   # Eq. (9)
    return q_over_dt * dx / kg                                  # ke0/kg


def ke0_iron(T, p=1.0):
    """Absolute ke0 for the Table 4 iron-sphere bed, kcal/(m hr degC).

    Parameters
    ----------
    T : float or ndarray
        Temperature, degC. `set_air_table` and `set_iron_sphere_row` must have
        been called.
    p : float, optional
        Emissivity-like factor p of the paper.

    Returns
    -------
    float or ndarray
        ke0 = kg(T) * Eq. (15).
    """
    iron_sphere_row = _IRON_SPHERE_ROW
    kg = kg_air(T)
    return kg * eq15(iron_sphere_row.eps, iron_sphere_row.beta, iron_sphere_row.gamma,
                     iron_sphere_row.phi, iron_sphere_row.Dp_mm * 1e-3,
                     iron_sphere_row.ks_kcal_m_hr_C, kg, h_rs(T, p), h_rv(T, iron_sphere_row.eps, p))


def reduction_row(Ti, To, n=400, method="newton", tol=1e-10, nu=1):
    """Eq. (19)'s reduction of the iron-sphere bed run in the annulus.

    Solves the variable-conductivity annulus with ``ke0_iron`` and reports the
    conductivity Eq. (19) would infer, against the true value at the
    volume-averaged temperature.

    Parameters
    ----------
    Ti, To : float
        Inner (hot) and outer (cold) wall temperatures, degC.
    n : int, optional
        Number of cells of the annulus solve.
    method, tol, nu : optional
        Passed to `solve_annulus` for the variable-k solve.

    Returns
    -------
    dict
        Wall temperatures, drop ``dT``, volume-averaged temperatures for variable
        and constant k, apparent conductivity ``k_app`` and its closed form
        ``k_app_closed``, ``k_at_tm``, ``bias_pct`` and the departures from a
        straight profile as a fraction of ``dT``.
    """
    dT = Ti - To
    r_f, r_c, T, Qp = solve_annulus(n, Ti, To, kfun=ke0_iron, method=method, tol=tol, nu=nu)
    vol = 2 * np.pi * r_c * np.diff(r_f); vol /= vol.sum()
    tm = float(np.sum(vol * T))
    k_app = float(Qp[0] / (2 * np.pi) * np.log(R_O / R_I) / dT)     # what Eq. (19) reports
    Tg = np.linspace(To, Ti, 20001)
    k_app_closed = float(np.trapezoid(ke0_iron(Tg), Tg) / dT)       # the identity above
    lin = np.polyval(np.polyfit(r_c, T, 1), r_c)
    r_fc, r_cc, Tc_, _ = solve_annulus(n, Ti, To, k_const=1.0, nu=nu)
    volc = 2 * np.pi * r_cc * np.diff(r_fc); volc /= volc.sum()
    linc = np.polyval(np.polyfit(r_cc, Tc_, 1), r_cc)
    return dict(t_out=To, t_in=Ti, dT=dT, tm_varK=tm, tm_constK=float(np.sum(volc * Tc_)),
                k_app=k_app, k_app_closed=k_app_closed,
                k_at_tm=float(ke0_iron(tm)), face_spread=float(np.ptp(Qp) / abs(Qp[0])), bias_pct=100 * (k_app / ke0_iron(tm) - 1),
                nonstraight_varK=float(np.max(np.abs(T - lin)) / dT),
                nonstraight_constK=float(np.max(np.abs(Tc_ - linc)) / dT))
