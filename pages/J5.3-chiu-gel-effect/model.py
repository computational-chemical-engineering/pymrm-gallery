"""Chiu, Carratt and Soong (1983): the gel effect without a switch (page J5.3).

Chiu, Carratt and Soong, Macromolecules 16, 348-357 (1983). The page needs:

- Table I's rate constants as functions of temperature (:func:`table1_constants`,
  :func:`k_d`, :func:`k_t0`, :func:`k_p0`, :func:`d_m`);
- the eight moment equations 10-17 with the constitutive equations 29-32, for
  one condition of Table II (:class:`Chiu`), with a BDF reference solve, the
  quasi-steady radical root of eq. 12 and the quasi-steady two-equation solve;
- the pymrm backward-Euler marcher on the eight fields (:func:`march`);
- eq. 21 as a 1-D spherical-shell boundary-value problem on pymrm operators
  (:func:`shell_bvp`);
- the moment-closure identity on a seeded finite chain population
  (:func:`moment_identity`);
- the paper's graphical onset construction, reconstructed
  (:func:`onset_from`, :func:`onset_window`, :func:`onset_x`), the Arrhenius
  fit of Figure 13 (:func:`arrhenius`) and a checked single-root finder
  (:func:`single_root`).

Every defect the page's break table injects enters through an argument of
these functions or an attribute of a :class:`Chiu` instance set before its
solve.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.sparse.linalg import spsolve
from pymrm import NumJac, construct_div, construct_grad, newton
from pymrm.checks import find_roots

R_CAL = 1.987          # cal/(mol K), the gas constant Table I writes into k_t0
MW_MMA = 100.117       # g/mol, methyl methacrylate; not printed in the paper


# ------------------------------------------------------------- root finding
def single_root(f, lo, hi, n_scan=64, log=False, xtol=1e-12):
    """The one root of ``f`` on ``[lo, hi]``, located by ``find_roots``.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        The physical range scanned.
    n_scan : int, optional
        Number of scan points of ``pymrm.checks.find_roots``.
    log : bool, optional
        Scan on a logarithmic grid.
    xtol : float, optional
        Absolute tolerance of the refined root.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        Unless exactly one root is found on the range.
    """
    res = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(res["roots"]) != 1:
        raise ValueError(f"expected one root on [{lo:g}, {hi:g}]: {res['message']}")
    return float(res["roots"][0])


def roots_with_direction(f, lo, hi, n_scan, xtol=1e-14):
    """All sign changes of ``f`` on ``[lo, hi]``, each with its direction.

    Parameters
    ----------
    f : callable
        Scalar function (vectorised or not).
    lo, hi : float
        The range scanned.
    n_scan : int
        Number of scan points of ``pymrm.checks.find_roots``.
    xtol : float, optional
        Absolute tolerance of each refined root.

    Returns
    -------
    list of tuple
        ``(root, +1)`` where ``f`` goes from negative to positive and
        ``(root, -1)`` where it goes from positive to negative, sorted.
    """
    res = find_roots(f, lo, hi, n_scan=n_scan, xtol=xtol)
    h = 1e-9 * (hi - lo)        # far inside one scan interval, far above xtol
    out = []
    for r in res["roots"]:
        a, b = float(f(max(lo, r - h))), float(f(min(hi, r + h)))
        out.append((r, 1 if b > a else -1))
    return out


# ------------------------------------------------------------------ Table I
def table1_constants(t1):
    """Table I's printed constants as one dictionary.

    Parameters
    ----------
    t1 : pandas.DataFrame
        Table I as read from its CSV, columns ``symbol``, ``a`` and ``b``.

    Returns
    -------
    dict
        ``f``, ``kd_a``, ``kd_b``, ``kt0_a``, ``kt0_e``, ``kp0_a``, ``kp0_e``,
        ``ktc``, ``d_p``, ``dm_a``, ``dm_b`` and ``T_gp``.
    """
    t = t1.set_index("symbol")

    def ab(sym):
        return float(t.loc[sym, "a"]), float(t.loc[sym, "b"])

    return {"f": float(t.loc["f", "a"]), "kd_a": ab("k_d")[0], "kd_b": ab("k_d")[1],
            "kt0_a": ab("k_t0")[0], "kt0_e": ab("k_t0")[1],
            "kp0_a": ab("k_p0")[0], "kp0_e": ab("k_p0")[1],
            "ktc": float(t.loc["k_tc", "a"]), "d_p": float(t.loc["d_p", "a"]),
            "dm_a": ab("d_m")[0], "dm_b": ab("d_m")[1], "T_gp": float(t.loc["T_gp", "a"])}


def k_d(T, c):
    """Initiator decomposition constant, 1/min (Table I prints L/min; eq. 1 is first order)."""
    return c["kd_a"] * np.exp(-c["kd_b"] / T)


def k_t0(T, c):
    """Intrinsic termination constant, L/(mol min)."""
    return c["kt0_a"] * np.exp(-c["kt0_e"] / (R_CAL * T))


def k_p0(T, c):
    """Intrinsic propagation constant, L/(mol min)."""
    return c["kp0_a"] * np.exp(-c["kp0_e"] / (R_CAL * T))


def d_m(T, c):
    """Monomer density, g/cm3, with the printed 273 (not 273.15)."""
    return c["dm_a"] - c["dm_b"] * (T - 273.0)


# ----------------------------------------------------------------- the model
class Chiu:
    """One experimental condition of Table II: eq. 10-17 with eq. 29-32.

    States, in the order the right-hand side uses them:
    ``I, x, lam0, lam1, lam2, mu0, mu1, mu2``. Eq. 18 (the quasi-steady
    primary radical) is already substituted, which turns ``k_i R M`` into
    ``2 f k_d I``.

    Parameters
    ----------
    T_C : float
        Polymerisation temperature, degC.
    I0 : float
        Initiator loading, mol/L.
    theta_t, theta_p, A, B : float
        Table II's parameters of eq. 31 and 32 (theta in min).
    const : dict
        Table I's constants, from :func:`table1_constants`.
    ktc : float, optional
        Termination by combination; Table I's value when omitted.
    ln10 : float, optional
        The factor eq. 31 and 32 print as 2.3.
    switch_xc : float, optional
        ``None`` is the published, branch-free model. A value gives the
        counterfactual with diffusion control switched on at that conversion.
    """

    FIELDS = ("I", "x", "lam0", "lam1", "lam2", "mu0", "mu1", "mu2")

    def __init__(self, T_C, I0, theta_t, theta_p, A, B, const, ktc=None,
                 ln10=2.3, switch_xc=None):
        self.T_C, self.I0 = float(T_C), float(I0)
        self.T = self.T_C + 273.15
        self.theta_t, self.theta_p = float(theta_t), float(theta_p)
        self.A, self.B = float(A), float(B)
        self.f = const["f"]
        self.ktc = const["ktc"] if ktc is None else ktc
        self.ln10 = ln10
        self.switch_xc = switch_xc
        self.kd, self.kt0, self.kp0 = k_d(self.T, const), k_t0(self.T, const), k_p0(self.T, const)
        self.dm = d_m(self.T, const)
        self.eps = (self.dm - const["d_p"]) / const["d_p"]   # volume expansion factor
        self.M0 = self.dm * 1000.0 / MW_MMA                  # mol/L

    def phi_m(self, x):
        """Eq. 30: volume fraction of monomer."""
        return (1.0 - x) / (1.0 + self.eps * x)

    def fujita(self, x):
        """D/D_0 of eq. 29 as eq. 31 and 32 write it, exp[2.3 phi_m/(A + B phi_m)]."""
        pm = self.phi_m(x)
        return np.exp(self.ln10 * pm / (self.A + self.B * pm))

    def kt_kp(self, x, P):
        """Eq. 31 and 32, one expression over the whole range (no branch)."""
        g = self.fujita(x)
        if self.switch_xc is None:
            return (1.0 / (1.0 / self.kt0 + self.theta_t * P / g),
                    1.0 / (1.0 / self.kp0 + self.theta_p * P / g))
        on = np.asarray(x) >= self.switch_xc
        return (np.where(on, 1.0 / (1.0 / self.kt0 + self.theta_t * P / g), self.kt0),
                np.where(on, 1.0 / (1.0 / self.kp0 + self.theta_p * P / g), self.kp0))

    def rhs(self, t, y):
        """Eq. 10-17, the time derivatives of the eight states."""
        I, x, l0, l1, l2, m0, m1, m2 = y
        x = min(max(x, 0.0), 1.0 - 1e-12)
        kt, kp = self.kt_kp(x, l0)
        ktd = kt - self.ktc                   # k_t = k_tc + k_td, eq. 4b
        e, M0 = self.eps, self.M0
        dil = e * (1.0 - x) * kp * l0 / (1.0 + e * x)   # (1/V) dV/dt
        M = M0 * (1.0 - x) / (1.0 + e * x)              # eq. 7
        gen = 2.0 * self.f * self.kd * I                # eq. 18 substituted
        return np.array([
            -self.kd * I - dil * I,                                    # eq. 10
            kp * (1.0 - x) * l0,                                       # eq. 11
            -dil * l0 + gen - kt * l0 * l0,                            # eq. 12
            -dil * l1 + gen - kt * l0 * l1 + kp * l0 * M,              # eq. 13
            -dil * l2 + gen - kt * l0 * l2 + kp * M * (2.0 * l1 + l0),  # eq. 14
            -dil * m0 + ktd * l0 * l0 + 0.5 * self.ktc * l0 * l0,      # eq. 15
            -dil * m1 + ktd * l0 * l1 + self.ktc * l0 * l1,            # eq. 16
            -dil * m2 + ktd * l0 * l2 + self.ktc * (l2 * l0 + l1 * l1),  # eq. 17
        ])

    def y0(self):
        """Book p. 350: I = I_0 and x = lam = mu = 0 at t = 0."""
        return np.array([self.I0, 0., 0., 0., 0., 0., 0., 0.])

    def solve(self, t_end, rtol=1e-11, atol=1e-20, x_stop=None):
        """BDF reference solve with dense output, optionally stopped at ``x = x_stop``.

        Raises
        ------
        RuntimeError
            When ``solve_ivp`` does not succeed.
        """
        ev = None
        if x_stop is not None:
            def ev(t, y):
                return y[1] - x_stop
            ev.terminal, ev.direction = True, 1
        s = solve_ivp(self.rhs, (0.0, t_end), self.y0(), method="BDF",
                      rtol=rtol, atol=atol, dense_output=True, events=ev)
        if not s.success:
            raise RuntimeError(f"BDF failed for {self.T_C} degC, {self.I0} mol/L: {s.message}")
        return s

    def lam0_qssa(self, x, I):
        """Eq. 12 with its left-hand side set to zero, solved for lambda_0.

        The residual ``2 f k_d I - lambda_0^2 (k_t + dilution)`` is strictly
        decreasing in lambda_0, because ``lambda_0^2 k_t(lambda_0)`` and
        ``lambda_0^2 k_p(lambda_0)`` both grow with lambda_0, so the root on
        ``[1e-18, 1e-1]`` mol/L is unique; :func:`single_root` checks it.
        """
        def h(lam):
            kt, kp = self.kt_kp(x, lam)
            return (-self.eps * (1.0 - x) * kp * lam * lam / (1.0 + self.eps * x)
                    + 2.0 * self.f * self.kd * I - kt * lam * lam)
        return single_root(h, 1e-18, 1e-1, n_scan=9, log=True, xtol=1e-28)

    def solve_qssa(self, t_end, rtol=1e-10, atol=1e-16):
        """Eq. 12-14 replaced by algebraic equations; only I and x are marched.

        Raises
        ------
        RuntimeError
            When ``solve_ivp`` does not succeed.
        """
        def rhs2(t, y):
            I, x = y
            x = min(max(x, 0.0), 1.0 - 1e-12)
            l0 = self.lam0_qssa(x, I)
            kt, kp = self.kt_kp(x, l0)
            dil = self.eps * (1.0 - x) * kp * l0 / (1.0 + self.eps * x)
            return [-self.kd * I - dil * I, kp * (1.0 - x) * l0]
        s = solve_ivp(rhs2, (0.0, t_end), [self.I0, 0.0], method="BDF",
                      rtol=rtol, atol=atol, dense_output=True)
        if not s.success:
            raise RuntimeError(f"QSSA solve failed: {s.message}")
        return s


# --------------------------------------------------------- pymrm marcher
def march(model, nt, t_end, y0=None, scale=None, tol=1e-12):
    """Backward Euler on the eight fields with pymrm's ``newton`` and ``NumJac``.

    ``NumJac((1, 8))``: the last axis is the field index, so the default
    stencil couples the eight fields in full and nothing else, the right
    shape for a pointwise source. The fields are scaled to order one before
    marching, because ``newton`` tests the infinity norm of the update and the
    states span twenty orders of magnitude.

    Parameters
    ----------
    model : Chiu
        The condition marched.
    nt : int
        Number of equal steps.
    t_end : float
        Length of the marched interval, min.
    y0 : ndarray, optional
        Start state; ``model.y0()`` when omitted.
    scale : ndarray, optional
        Scale of each field; the start state when omitted.
    tol : float, optional
        ``newton`` step tolerance on the scaled fields.

    Returns
    -------
    ndarray
        The eight states at the end of the interval.

    Raises
    ------
    RuntimeError
        When a step's ``newton`` solve does not succeed.
    """
    dt = t_end / nt
    nf = len(model.FIELDS)
    jac = NumJac((1, nf))
    y = (model.y0() if y0 is None else y0).reshape(1, -1).copy()
    sc = (np.maximum(np.abs(np.asarray(scale, float)), 1e-30) if scale is not None
          else np.maximum(np.abs(y.ravel()), 1e-12))
    u = y / sc

    def src(uu):
        out = np.empty_like(uu)
        for k in range(uu.shape[0]):
            out[k] = model.rhs(0.0, uu[k] * sc) / sc
        return out

    for step in range(nt):
        u_old = u.copy()

        def resid(uu):
            s, js = jac(src, uu.reshape(1, -1))
            r = (uu.reshape(1, -1) - u_old) / dt - s
            return r.reshape(-1, 1), np.eye(nf) / dt - js.toarray()

        sol = newton(resid, u, tol=tol, maxfev=200)
        if not sol.success:
            raise RuntimeError(f"backward-Euler step {step} did not converge: {sol.message}")
        u = sol.x.reshape(1, -1)
    return (u * sc).ravel()


# ------------------------------------------------- eq. 21, spherical shell
def shell_bvp(n, ratio, Cm_=0.3, Cb_=1.0, D_=1.0, rm_=1.0, nu=2):
    """Eq. 21 on the shell ``r_m <= r <= r_D`` with pymrm operators.

    The problem is linear, so one sparse solve gives the profile; no Newton
    iteration is involved.

    Parameters
    ----------
    n : int
        Number of uniform cells.
    ratio : float
        ``r_m / r_D``.
    Cm_, Cb_ : float, optional
        Dirichlet values at ``r_m`` and ``r_D``.
    D_, rm_ : float, optional
        Diffusivity and inner radius.
    nu : int, optional
        Geometry of ``construct_div``: 2 is spherical, the geometry of eq. 21.

    Returns
    -------
    dict
        ``err`` (rms error of the cell values against the exact profile at
        the cell centres), ``K_spread`` (range of the face fluxes),
        ``K_mid`` (flux ``4 pi r^2 D dC/dr`` at the middle face), ``K_exact``,
        ``K_approx`` (eq. 22) and ``r_mid`` (radius of the middle face).
    """
    rD_ = rm_ / ratio
    shape = (n, 1)
    r_f = np.linspace(rm_, rD_, n + 1)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    # a dC/dn + b C = d on the outward normal; both ends are Dirichlet (a = 0):
    #   r = r_m : C = C_m          r = r_D : C = C_b
    bc = ({"a": 0.0, "b": 1.0, "d": Cm_}, {"a": 0.0, "b": 1.0, "d": Cb_})
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc, axis=0)
    div = construct_div(shape, r_f, nu=nu, axis=0)       # nu = 2: spherical
    c = spsolve(div @ (-D_ * grad),
                -(div @ (-D_ * grad_bc)).toarray().ravel())
    dCdr = np.asarray((grad @ c.reshape(-1, 1) + grad_bc)).ravel()
    K_faces = 4 * np.pi * r_f**2 * D_ * dCdr               # eq. 21's left side
    C_exact = Cb_ + (Cm_ - Cb_) * (1 / r_c - 1 / rD_) / (1 / rm_ - 1 / rD_)
    K_exact = 4 * np.pi * D_ * (Cb_ - Cm_) / (1 / rm_ - 1 / rD_)
    K_approx = 4 * np.pi * D_ * rm_ * (Cb_ - Cm_)          # eq. 22
    return dict(err=float(np.sqrt(np.mean((c - C_exact)**2))),
                K_spread=float(K_faces.max() - K_faces.min()),
                K_mid=float(K_faces[n // 2]), K_exact=float(K_exact),
                K_approx=float(K_approx), r_mid=float(r_f[n // 2]))


# ------------------------------------------------- the moment identity
def moment_identity(N=60, ktc=0.0, seed=19830348, break_term=None, support=None):
    """Eq. 12-17 against the chain-length equations 4a, 4b and 5 summed.

    A seeded random population ``P_1..P_N`` with a zero upper half, so the
    telescoping of eq. 4b and the convolution of eq. 5 truncate exactly.

    Parameters
    ----------
    N : int, optional
        Number of chain lengths.
    ktc : float, optional
        Termination by combination, L/(mol min).
    seed : int, optional
        Seed of the population.
    break_term : str, optional
        Moment term (``"l0"`` ... ``"m2"``) whose sign is flipped.
    support : int, optional
        First zero chain length; ``N // 2`` when omitted. A larger value
        violates the convolution truncation.

    Returns
    -------
    tuple
        Largest relative residual, the summed chain-length side and the
        moment side.
    """
    rng = np.random.default_rng(seed)
    n = np.arange(1, N + 1)
    P = rng.random(N) * np.exp(-n / 12.0)
    P[(N // 2 if support is None else support):] = 0.0
    kp, kt, M, gen = 3.1e4, 1.7e9, 8.4, 2.9e-6
    ktd = kt - ktc
    lam = np.array([np.sum(n**k * P) for k in (0, 1, 2)])
    dP = np.empty(N)
    dP[0] = gen - kp * M * P[0] - kt * P[0] * lam[0]
    dP[1:] = kp * M * (P[:-1] - P[1:]) - kt * P[1:] * lam[0]
    conv = np.array([np.sum(P[:i] * P[i - 1::-1][:i]) for i in range(N)])
    dMn = ktd * P * lam[0] + 0.5 * ktc * conv
    lhs = np.array([np.sum(n**k * dP) for k in (0, 1, 2)]
                   + [np.sum(n**k * dMn) for k in (0, 1, 2)])
    c = dict(l0=gen - kt * lam[0]**2,
             l1=gen - kt * lam[0] * lam[1] + kp * M * lam[0],
             l2=gen - kt * lam[0] * lam[2] + kp * M * (2 * lam[1] + lam[0]),
             m0=ktd * lam[0]**2 + 0.5 * ktc * lam[0]**2,
             m1=ktd * lam[0] * lam[1] + ktc * lam[0] * lam[1],
             m2=ktd * lam[0] * lam[2] + ktc * (lam[2] * lam[0] + lam[1]**2))
    if break_term is not None:
        c[break_term] = -c[break_term]
    rhs = np.array([c["l0"], c["l1"], c["l2"], c["m0"], c["m1"], c["m2"]])
    scale = np.maximum(np.abs(lhs), 1e-300)
    return float(np.max(np.abs(rhs - lhs) / scale)), lhs, rhs


# ---------------------------------------------------- Figure 13 and onsets
def arrhenius(theta, T_K):
    """Arrhenius line through ``ln(1/theta)`` against ``1/T``.

    Returns
    -------
    tuple
        Activation energy in kcal/mol, rms and largest ln-residual.
    """
    y = np.log(1.0 / np.asarray(theta, float))
    X = np.vstack([np.ones_like(T_K), 1.0 / T_K]).T
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    return (-c[1] * R_CAL / 1000.0, float(np.sqrt((res**2).mean())),
            float(np.max(np.abs(res))))


N_FEATURE_SCAN = 20000     # scan points for the curvature features of log10 k_t


def onset_from(x, L, L0, mode="flattest"):
    """The paper's graphical onset construction, as a tangent extrapolated to L0.

    The bend is the first local minimum of ``L''`` (an upward sign change of
    ``L'''``), the flattest point of the sloped branch the first local
    maximum of ``L''`` after it (a downward sign change of ``L'''``); both are
    roots from ``find_roots``, the first of their kind, not the largest.

    Parameters
    ----------
    x, L : ndarray
        Conversion grid and log10 k_t on it.
    L0 : float
        The initial value the tangent is extrapolated to.
    mode : str, optional
        ``"flattest"``, ``"decade3"`` (the point three decades below L0) or
        ``"steepest"`` (the steepest point of the window, the edge included).

    Returns
    -------
    dict
        ``x_onset``, ``x_tangent``, ``slope``, ``x_bend``, ``L_tangent`` and
        ``refined`` (``"root"`` or ``"edge"``).

    Raises
    ------
    ValueError
        When the feature the reading needs is not in the window.
    """
    cs = CubicSpline(x, L)
    d1, d2, d3 = cs.derivative(1), cs.derivative(2), cs.derivative(3)
    feats = roots_with_direction(d3, x[0], x[-1], N_FEATURE_SCAN)
    ups = [r for r, s in feats if s > 0]
    if not ups:
        raise ValueError("no local minimum of L'' (the bend) in the window")
    x_bend = ups[0]
    refined = "root"
    if mode == "flattest":
        downs = [r for r, s in feats if s < 0 and r > x_bend]
        if not downs:
            raise ValueError("no local maximum of L'' past the bend in the window")
        xs = downs[0]
    elif mode == "decade3":
        xs = single_root(lambda v: float(cs(v)) - (L0 - 3.0), x_bend, x[-1],
                         n_scan=400, xtol=1e-14)
    elif mode == "steepest":
        cand = [r for r, s in roots_with_direction(d2, x_bend, x[-1], N_FEATURE_SCAN)
                if s > 0]
        cand.append(float(x[-1]))
        xs = min(cand, key=lambda v: float(d1(v)))
        refined = "edge" if xs == float(x[-1]) else "root"
    else:
        raise ValueError(mode)
    return dict(x_onset=float(xs + (L0 - cs(xs)) / d1(xs)), x_tangent=float(xs),
                slope=float(d1(xs)), x_bend=float(x_bend), L_tangent=float(cs(xs)),
                refined=refined)


def onset_window(x, L, L0, lo=2.0, hi=4.0):
    """The primary reading: a least-squares line through log10 k_t over the
    decade window ``L0 - lo .. L0 - hi``, extrapolated to L0.

    Returns
    -------
    tuple
        Onset conversion, number of points in the window and the slope.

    Raises
    ------
    ValueError
        When the window holds fewer than ten points.
    """
    m = (L <= L0 - lo) & (L >= L0 - hi)
    if m.sum() < 10:
        raise ValueError("the decade window holds too few points")
    p = np.polyfit(x[m], L[m], 1)
    return float((L0 - p[1]) / p[0]), int(m.sum()), float(p[0])


def onset_x(x, L, L0, mode="window24"):
    """One name for every reading of 'extrapolating the sloped region'."""
    windows = {"window24": (2.0, 4.0), "window13": (1.0, 3.0),
               "window35": (3.0, 5.0), "window46": (4.0, 6.0)}
    if mode in windows:
        return onset_window(x, L, L0, *windows[mode])[0]
    return onset_from(x, L, L0, mode)["x_onset"]
