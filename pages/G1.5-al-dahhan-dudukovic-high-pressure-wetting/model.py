"""Al-Dahhan and Dudukovic (1995): catalyst wetting efficiency at high pressure.

Closed forms of the paper (eq. 9 and El-Hisnawi's eq. 1), the error metrics and
least-squares refits the page applies to digitised data, the input checks of
Tables 3 and 4, the tracer-moment relations of Tables 1 and 2, and one pymrm
model: the pellet first-moment boundary value problem that gives the moment
coefficients 1/15, 1/8 and 1/3 behind those tables.

Every printed constant enters through a dictionary ``K`` (keys ``c9``, ``n9``,
``m9`` for eq. 9 and ``c1``, ``n1``, ``m1`` for eq. 1), so a mistyped constant is
one argument. SI units throughout: L in kg/(m2 s), lengths in m, rho in kg/m3,
mu in Pa s; ``dpz`` is (dP/Z)/(rho_L g), dimensionless, as plotted in Figs 4b, 5b.
"""
import numpy as np
from scipy.sparse.linalg import spsolve
from pymrm import construct_div, construct_grad

G_GRAV = 9.80665          # m/s2, standard gravity
R_GAS = 8.314462618       # J/(mol K)
M_N2 = 28.0134e-3         # kg/mol
PA_PER_PSI = 6894.757293168   # Pa per psi: the exact 6894.757293168361 rounded to 13 digits, as in the check-value script (relative difference 5e-14)
P_ATM_MPA = 0.101325      # MPa


# ---- eq. 9 (Al-Dahhan and Dudukovic) and eq. 1 (El-Hisnawi) -------------------------

def psi_l(dpz, one_plus=True):
    """psi_L = 1 + (dP/Z)/(rho_L g) (Notation); ``one_plus=False`` drops the 1."""
    return (1.0 if one_plus else 0.0) + np.asarray(dpz, dtype=float)


# codehealth: allow params -- the groups of eq. 9; voidage=False is the defect the page's break table injects
def re9(L, d_p, eps_B, mu, voidage=True):
    """Re_L of eq. 9, U_L rho_L d_p/(mu_L (1 - eps_B)); without (1 - eps_B) if not voidage."""
    return np.asarray(L, dtype=float) * d_p / (mu * ((1.0 - eps_B) if voidage else 1.0))


# codehealth: allow params -- the groups of eq. 9; voidage=False is the defect the page's break table injects
def ga9(d_p, eps_B, rho, mu, g=G_GRAV, voidage=True):
    """Ga_L of eq. 9, d_p^3 rho^2 g eps_B^3/(mu^2 (1 - eps_B)^3); plain d_p^3 rho^2 g/mu^2 if not voidage."""
    v = (eps_B / (1.0 - eps_B)) ** 3 if voidage else 1.0
    return d_p ** 3 * rho ** 2 * g * v / mu ** 2


# codehealth: allow params -- eq. 9 as printed: flow, bed, liquid and the printed constants; cap, voidage and one_plus are the readings the page tests
def eta9(L, dpz, d_p, eps_B, rho, mu, K, cap=True, voidage=True, one_plus=True, g=G_GRAV):
    """Eq. 9: eta_CE = c9 Re_L^n9 (psi_L/Ga_L)^m9, capped at 1 when ``cap``."""
    e = K["c9"] * re9(L, d_p, eps_B, mu, voidage) ** K["n9"] * \
        (psi_l(dpz, one_plus) / ga9(d_p, eps_B, rho, mu, g, voidage)) ** K["m9"]
    return np.minimum(e, 1.0) if cap else e


def k9(rho, mu, K, g=G_GRAV):
    """Prefactor of the reduced form (n9 = 3 m9): c9 (mu rho^2 g)^(-m9), SI."""
    return K["c9"] * (mu * rho ** 2 * g) ** (-K["m9"])


# codehealth: allow params -- the reduced form of eq. 9 with the same inputs as eq. 9, compared with it
def eta9_reduced(L, dpz, eps_B, rho, mu, K, g=G_GRAV):
    """Eq. 9 reduced for n9 = 3 m9: K L^n9 eps_B^(-3 m9) psi^m9 (d_p and 1 - eps_B cancel)."""
    return k9(rho, mu, K, g) * np.asarray(L, dtype=float) ** K["n9"] * eps_B ** (-3 * K["m9"]) * \
        psi_l(dpz) ** K["m9"]


# codehealth: allow params -- eq. 9 solved for L, with the same inputs as eq. 9
def l_full9(dpz, d_p, eps_B, rho, mu, K, g=G_GRAV):
    """The L at which uncapped eq. 9 reaches 1 (closed form)."""
    re = (1.0 / (K["c9"] * (psi_l(dpz) / ga9(d_p, eps_B, rho, mu, g)) ** K["m9"])) ** (1.0 / K["n9"])
    return re * mu * (1.0 - eps_B) / d_p


# codehealth: allow params -- eq. 1 as printed; eps_B switches to the voidage groups, a defect the page tests
def eta1(L, d_p, rho, mu, K, g=G_GRAV, eps_B=None):
    """Eq. 1 (El-Hisnawi): 1.617 Re_L^0.1461 Ga_L^-0.0711 with Re_L = L d_p/mu, Ga_L = d_p^3 g rho^2/mu^2.
    With ``eps_B`` given, the voidage groups of eq. 9 are used instead (a defect)."""
    if eps_B is None:
        re, ga = np.asarray(L, dtype=float) * d_p / mu, d_p ** 3 * g * rho ** 2 / mu ** 2
    else:
        re, ga = re9(L, d_p, eps_B, mu), ga9(d_p, eps_B, rho, mu, g)
    return K["c1"] * re ** K["n1"] * ga ** K["m1"]


def eq1_exponents(K):
    """Exponents of d_p and mu_L in eq. 1 at fixed L, rho_L."""
    return K["n1"] + 3 * K["m1"], -K["n1"] - 2 * K["m1"]


# ---- metrics ------------------------------------------------------------------------

def deviations(pred, exp, convention="exp"):
    """(p - e)/e (convention 'exp') or (p - e)/p (convention 'pred')."""
    p, e = np.asarray(pred, dtype=float), np.asarray(exp, dtype=float)
    return (p - e) / (e if convention == "exp" else p)


def metrics(pred, exp, convention="exp", band=0.10):
    """Mean relative error, largest |deviation|, count with |deviation| <= band and the
    smallest distance of any |deviation| to the band edge."""
    d = np.abs(deviations(pred, exp, convention))
    return {"mre": float(d.mean()), "maxdev": float(d.max()), "inside": int(np.sum(d <= band)),
            "margin": float(np.min(np.abs(d - band)))}


def mean_se(x):
    """Mean and its standard error (sample SD over sqrt(n))."""
    x = np.asarray(x, dtype=float)
    return float(x.mean()), float(x.std(ddof=1) / np.sqrt(x.size))


# ---- refits (E9) --------------------------------------------------------------------

def ols(X, y):
    """Ordinary least squares: coefficients and their standard errors."""
    X, y = np.asarray(X, dtype=float), np.asarray(y, dtype=float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    dof = X.shape[0] - X.shape[1]
    cov = (r @ r / dof) * np.linalg.inv(X.T @ X)
    return beta, np.sqrt(np.diag(cov))


# codehealth: allow params -- the groups of eq. 6 need the bed and the liquid
def design_eq6(L, dpz, d_p, eps_B, rho, mu, g=G_GRAV):
    """Columns of eq. 6 in logs: 1, ln Re_L, ln(psi_L/Ga_L) (eq. 9 groups)."""
    return np.column_stack([np.ones(len(L)), np.log(re9(L, d_p, eps_B, mu)),
                            np.log(psi_l(dpz) / ga9(d_p, eps_B, rho, mu, g))])


def design_bed(L, dpz, bed):
    """Columns of the bed-offset form: one intercept per bed, ln L, ln psi_L."""
    bed = np.asarray(bed)
    return np.column_stack([(bed == "sphere").astype(float), (bed == "extrudate").astype(float),
                            np.log(L), np.log(psi_l(dpz))])


# codehealth: allow params -- design, response, clusters, draw count and generator of one bootstrap
def cluster_bootstrap(X, y, cluster, draws, rng):
    """Refit on resamples of whole clusters; returns the coefficient draws (rows) and
    the number of singular resamples (rank-deficient design, e.g. one bed only)."""
    ids = np.unique(cluster)
    out, singular = [], 0
    for _ in range(draws):
        pick = rng.choice(ids, size=ids.size, replace=True)
        idx = np.concatenate([np.nonzero(cluster == c)[0] for c in pick])
        Xb = X[idx]
        if np.linalg.matrix_rank(Xb) < X.shape[1]:
            singular += 1
            continue
        out.append(np.linalg.lstsq(Xb, y[idx], rcond=None)[0])
    return np.array(out), singular


def leave_one_out(X, y, cluster):
    """Refit with each cluster left out in turn."""
    return np.array([np.linalg.lstsq(X[cluster != c], y[cluster != c], rcond=None)[0]
                     for c in np.unique(cluster)])


def interp_loglog(x, xs, ys):
    """Log-log interpolation of (xs, ys) at x inside [min xs, max xs]; nan outside."""
    lx, lxs, lys = np.log(np.atleast_1d(x)), np.log(xs), np.log(ys)
    o = np.argsort(lxs)
    lxs, lys = lxs[o], lys[o]
    out = np.full(lx.shape, np.nan)
    for i, v in enumerate(lx):
        j = np.searchsorted(lxs, v)
        inside = lxs[0] <= v <= lxs[-1]
        if inside:
            j = min(max(j, 1), lxs.size - 1)
            t = (v - lxs[j - 1]) / (lxs[j] - lxs[j - 1])
            out[i] = np.exp(lys[j - 1] + t * (lys[j] - lys[j - 1]))
    return out


def interp_linear(x, xs, ys):
    """Linear interpolation of (xs, ys) at x inside the range; nan outside."""
    return np.log(interp_loglog(np.exp(np.atleast_1d(x)), np.exp(xs), np.exp(ys)))


# ---- input checks (Tables 3 and 4) --------------------------------------------------

def deq_surface(d, l):
    """Surface-volume diameter 6 V_p/S_ex of a cylinder d x l (end faces included)."""
    return 3.0 * d * l / (2.0 * l + d)


def d_volume(d, l):
    """Diameter of the sphere of equal volume."""
    return (1.5 * d * d * l) ** (1.0 / 3.0)


def shape_factor(d, l):
    """S_ex r_p/(2 V_p) of Table 1 for a cylinder d x l: 1 + r/l."""
    return 1.0 + 0.5 * d / l


def psig_to_mpa(psig):
    """Gauge psi to absolute MPa."""
    return psig * PA_PER_PSI * 1e-6 + P_ATM_MPA


def half_units(value, printed):
    """|value - printed| in half-units of the last printed digit (printed is a string)."""
    s = printed.strip()
    dec = len(s.split(".")[1]) if "." in s else 0
    return abs(value - float(s)) / (0.5 * 10.0 ** (-dec))


def gas_flux_ideal(p_mpa, u_cm_s, T=298.0):
    """N2 mass flux rho_N2 U_G, ideal gas, kg/(m2 s)."""
    return p_mpa * 1e6 * M_N2 / (R_GAS * T) * u_cm_s * 1e-2


def rho_implied(L, q_ml_min, d_r):
    """Liquid density implied by a mass flux L and a volumetric flow Q_L in a tube d_r."""
    return L / (q_ml_min * 1e-6 / 60.0 / (np.pi * d_r ** 2 / 4.0))


# ---- the pellet first moment (E13, E14) with pymrm ----------------------------------

# codehealth: allow params -- grid, geometry, film and the three defect or solver switches of the break table
def pellet_moment_coef(n, nu, bi=None, weight="volume", permc_spec=None, dtype=float):
    """tau_1/(eps_p R^2/De) from the pellet BVP (1/x^nu) d/dx(x^nu du/dx) = 1.

    Finite volumes on n uniform cells of [0, 1]; ``nu``: 0 slab, 1 cylinder, 2 sphere.
    Centre: du/dx = 0. Surface: u = 0 (``bi=None``) or du/dx + Bi u = 0 (film).
    The coefficient is minus the volume average of u; ``weight="cell"`` uses plain
    cell widths instead of the exact volumes and ``weight="drop_outer"`` leaves the
    outermost cell out of the average (defects); ``dtype=np.float32`` solves in
    single precision (a round-off defect).
    """
    x_f = np.linspace(0.0, 1.0, n + 1)
    x_c = 0.5 * (x_f[1:] + x_f[:-1])
    bc = ({"a": 1, "b": 0, "d": 0},               # x = 0, outward normal -x: -du/dx = 0 (symmetry)
          {"a": 0, "b": 1, "d": 0} if bi is None  # x = 1: u = 0 (no film)
          else {"a": 1, "b": bi, "d": 0})         # x = 1, outward normal +x: du/dx + Bi u = 0
    grad, grad_bc = construct_grad((n, 1), x_f, x_c, bc)
    div = construct_div((n, 1), x_f, nu=nu)       # nu: 0 slab, 1 cylinder, 2 sphere
    A = (div @ grad).astype(dtype)
    rhs = (np.ones(n) - div @ grad_bc.toarray().ravel()).astype(dtype)   # grad u = grad @ u + grad_bc
    if permc_spec == "dense":                     # LAPACK LU of the dense matrix: a second solver
        u = np.linalg.solve(A.toarray(), rhs)
    else:
        u = spsolve(A, rhs) if permc_spec is None else spsolve(A, rhs, permc_spec=permc_spec)
    w = np.diff(x_f ** (nu + 1)) if weight != "cell" else np.diff(x_f)
    if weight == "drop_outer":
        w = w.copy()
        w[-1] = 0.0
    return float(-(w @ u) / w.sum())


MOMENT_EXACT = {0: 1.0 / 3.0, 1: 1.0 / 8.0, 2: 1.0 / 15.0}      # tau_1 coefficient, no film
FILM_EXACT = {0: 1.0, 1: 1.0 / 2.0, 2: 1.0 / 3.0}               # film term times Bi


# codehealth: allow params -- the reactor, bed, pellet and flow of the tracer moment relations, plus the factor-2 switch
def tracer_standard(Vr, eps, eps_B, eps_p, R, De, Q, nu=2, bi=None, factor2=True):
    """Plug flow through a bed of pellets (E13): first moment (s), variance (s2) and
    sigma_D^2 = variance/mu_1^2. ``eps`` is the external liquid fraction (eps_B liquid
    filled, eps_L trickle flow). ``factor2=False`` drops the 2 of the second cumulant."""
    t_f = Vr * eps / Q
    d0 = (1.0 - eps_B) * eps_p / eps
    tau1 = eps_p * R ** 2 / De * (MOMENT_EXACT[nu] + (FILM_EXACT[nu] / bi if bi else 0.0))
    mu1 = t_f * (1.0 + d0)
    var = (2.0 if factor2 else 1.0) * t_f * d0 * tau1
    return mu1, var, var / mu1 ** 2


# codehealth: allow params -- the reactor, bed, pellet and flow of Table 2 as printed
def tracer_table2(Vr, eps, eps_B, eps_p, R, De, Q):
    """Table 2 as printed, no film and no axial dispersion: sigma_D^2 and the slope
    of sigma_D^2 against Q_L."""
    d2 = (1.0 - eps_B) * eps_p / eps
    sig = (Q / (Vr * eps)) * d2 * R ** 2 * eps_p / (15.0 * De) / (1.0 + d2) ** 2
    slope = ((1.0 - eps_B) * eps_p ** 2 / eps) * (R ** 2 / (15.0 * De)) / (Vr * eps * (1.0 + d2) ** 2)
    return sig, slope


# codehealth: allow params -- the inputs of the tracer form that is inverted, plus the form and the factor-2 switch
def de_recovered(sig_measured, Vr, eps, eps_B, eps_p, R, Q, form="table2", factor2=True):
    """De that makes a form reproduce a measured sigma_D^2 (both are proportional to 1/De)."""
    if form == "table2":
        s1 = tracer_table2(Vr, eps, eps_B, eps_p, R, 1.0, Q)[0]
    else:
        s1 = tracer_standard(Vr, eps, eps_B, eps_p, R, 1.0, Q, factor2=factor2)[2]
    return s1 / sig_measured
