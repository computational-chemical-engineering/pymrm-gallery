"""Smith and Ewart (1948): model functions for page J5.6.

Part I of the paper: the steady distribution of radicals over latex particles (eq. 1) in
the page groups ``alpha = rho' v/(N k_t)`` (entry), ``m = k_0 a/k_t`` (desorption), with
the closed form in modified Bessel functions (later work: Stockmayer 1957 for m = 0,
O'Toole 1965 in general, not Smith and Ewart's), the truncated recursion, the slope
``d ln nbar/d ln alpha``, the paper's Case 2 construction of Appendix I and exact power
series. Part II: the Stage I particle number bounds (eqs. 21 to 30, 26a) through
Mittag-Leffler series, and a population balance in particle volume solved with pymrm
(upwind convection, implicit Euler, ``pymrm.newton``). No plotting, no data loading.
"""

import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from scipy.special import gamma, ive, rgamma
from scipy.sparse import csc_array, diags_array, eye_array, vstack
from scipy.sparse.linalg import splu
from pymrm import construct_convflux_upwind, construct_div, newton


# ----------------------------------------------------------------------------- Part I

@dataclass(frozen=True)
class Law:
    """The Part I rate laws; each field is one reading of eq. 1.

    Attributes
    ----------
    term : float
        Factor on the termination events ``k_t n(n-1)/v`` (1: the paper's eq. 1).
    desorb : str
        ``"linear"``: desorption ``k_0 a n/v`` (eq. 1); ``"constant"``: ``k_0 a/v`` for
        every occupied particle (a defect).
    leak : bool
        Truncation that lets entry leave the last class (a defect).
    """
    term: float = 1.0
    desorb: str = "linear"
    leak: bool = False


LAW = Law()


def _ratio(m, a):
    """I_m(a)/I_(m-1)(a), exponentially scaled Bessel functions."""
    return ive(m, a) / ive(m - 1.0, a)


def nbar_closed(alpha, m, law=LAW):
    """Mean radicals per particle from the closed form, (a/4) I_m(a)/I_(m-1)(a), a = sqrt(8 alpha).

    A termination factor c is the same balance with alpha/c and m/c.
    """
    al, mm = alpha / law.term, m / law.term
    a = np.sqrt(8.0 * al)
    return a / 4.0 * _ratio(mm, a)


def p01_closed(alpha, m, law=LAW):
    """Fractions of empty and singly occupied particles, p_0 and p_1, from the closed form."""
    al, mm = alpha / law.term, m / law.term
    a, b = np.sqrt(8.0 * al), 2.0 * np.sqrt(al)
    p0 = 2.0 ** (-(1.0 - mm) / 2.0) * ive(mm - 1.0, b) / ive(mm - 1.0, a) * np.exp(b - a)
    return p0, p0 * np.sqrt(al) * _ratio(mm, b)


def slope_closed(alpha, m, law=LAW):
    """Local slope s = d ln nbar/d ln alpha = 1 - m + (a/2)(1/R - R), R = I_m/I_(m-1)."""
    mm = m / law.term
    a = np.sqrt(8.0 * alpha / law.term)
    r = _ratio(mm, a)
    return 1.0 - mm + a / 2.0 * (1.0 / r - r)


def dslope_da(a, m):
    """Analytic derivative ds/da of the slope at fixed m (a = sqrt(8 alpha)).

    Uses R' = 1 - R^2 - (2m - 1) R/a, from I_m' = I_(m-1) - (m/a) I_m and
    I_(m-1)' = I_m + ((m-1)/a) I_(m-1).
    """
    r = _ratio(m, a)
    dr = 1.0 - r * r - (2.0 * m - 1.0) * r / a
    return 0.5 * (1.0 / r - r) + a / 2.0 * (-1.0 / r**2 - 1.0) * dr


def recursion_matrix(alpha, m, nmax, law=LAW):
    """Eq. 1 for n = 0..nmax in the page groups, equation 0 replaced by sum p = 1.

    Returns a CSC matrix A and the right-hand side e_0 with A p = e_0. No entry leaves
    the last class unless ``law.leak``.
    """
    n = np.arange(nmax + 1, dtype=float)
    d = n if law.desorb == "linear" else (n > 0).astype(float)
    ent = np.full(nmax + 1, float(alpha))
    if not law.leak:
        ent[-1] = 0.0
    # bands of equations 1..nmax (row 0 zeroed: it becomes sum p = 1)
    main = -ent - m * d - law.term * n * (n - 1.0)
    sub = np.full(nmax, float(alpha))                       # p_(n-1) in equation n
    sup1 = m * d[1:]                                        # p_(n+1) in equation n
    sup2 = law.term * n[2:] * (n[2:] - 1.0)                 # p_(n+2) in equation n
    main[0] = 0.0
    sup1[0] = 0.0
    sup2[0] = 0.0
    bands = diags_array([sub, main, sup1, sup2][: 4 if nmax >= 2 else 3],
                        offsets=[-1, 0, 1, 2][: 4 if nmax >= 2 else 3], shape=(nmax + 1, nmax + 1), format="csc")
    norm = vstack([csc_array(np.ones((1, nmax + 1))), csc_array((nmax, nmax + 1))], format="csc")
    A = bands + norm
    rhs = np.zeros(nmax + 1)
    rhs[0] = 1.0
    return A, rhs


def equation0_residual(p, alpha, m, law=LAW):
    """Residual of the replaced equation 0 of eq. 1 (entry, desorption, termination into 0)."""
    p = np.asarray(p)
    d1 = 1.0
    r = -alpha * p[0] + m * d1 * p[1] if p.size > 1 else -alpha * p[0]
    if p.size > 2:
        r += law.term * 2.0 * p[2]
    return r


def recursion(alpha, m, nmax, law=LAW, method="splu", permc_spec="COLAMD"):
    """Distribution p_0..p_nmax from the truncated recursion (a linear solve).

    ``method`` "splu" (sparse LU with ``permc_spec``) or "dense" (LAPACK).
    """
    A, rhs = recursion_matrix(alpha, m, nmax, law)
    if method == "dense":
        return np.linalg.solve(A.toarray(), rhs)
    return splu(A, permc_spec=permc_spec).solve(rhs)


def nmax_for(alpha, m, law=LAW, tol=1e-17, start=None, cap=None):
    """Smallest checked n_max with p_nmax/max p below ``tol`` (doubling search), and p."""
    nmax = start or int(20 + 6.0 * math.sqrt(alpha / law.term))
    while True:
        if cap is not None:
            nmax = min(nmax, cap)
        p = recursion(alpha, m, nmax, law)
        if cap is not None and nmax >= cap:
            return nmax, p
        if abs(p[-1]) / np.max(np.abs(p)) < tol:
            return nmax, p
        nmax = int(nmax * 1.5) + 5


def stats(p):
    """nbar, p_0, p_1, termination sum S2 = sum n(n-1) p_n."""
    n = np.arange(p.size, dtype=float)
    return {"nbar": float(n @ p), "p0": float(p[0]), "p1": float(p[1]),
            "S2": float((n * (n - 1.0)) @ p), "sum": float(np.sum(p))}


def dist(alpha, m, law=LAW, cap=None):
    """(nbar, p_0, p_1) at (alpha, m): closed form for desorption proportional to n,
    otherwise the converged recursion."""
    if law.desorb == "linear" and not law.leak and cap is None:
        p0, p1 = p01_closed(alpha, m, law)
        return float(nbar_closed(alpha, m, law)), float(p0), float(p1)
    _, p = nmax_for(alpha, m, law, cap=cap)
    s = stats(p)
    return s["nbar"], s["p0"], s["p1"]


def appendix_nbar(beta, eq13_from=3, nterms=40):
    """nbar of the paper's Case 2 construction (m = 0, eq. 12 exact at n = 0 and 1,
    eq. 13 from ``eq13_from`` on), N_n continued by eq. 13 to ``nterms``.

    Returns (nbar, N_1/N_0).
    """
    N = np.zeros(nterms + 1)
    N[0] = 1.0
    if eq13_from == 3:
        N[2] = 1.0 / (2.0 * beta)                       # eq. 12, n = 0 (exact)
        N[3] = N[2] / (1.0 + 6.0 * beta)                # eq. 13 at n = 3
        N[1] = 1.0 + 6.0 * beta * N[3]                  # eq. 12, n = 1
    elif eq13_from == 2:                                 # eq. 13 also at n = 2 (a defect)
        g = 6.0 * beta / ((1.0 + 2.0 * beta) * (1.0 + 6.0 * beta))
        N[1] = 1.0 / (1.0 - g)
        N[2] = N[1] / (1.0 + 2.0 * beta)
        N[3] = N[2] / (1.0 + 6.0 * beta)
    else:
        raise ValueError("eq13_from must be 2 or 3")
    for n in range(4, nterms + 1):
        N[n] = N[n - 1] / (1.0 + beta * n * (n - 1.0))
    n = np.arange(nterms + 1)
    return float(n @ N / np.sum(N)), float(N[1] / N[0])


# exact power series in alpha = 1/beta (Fractions), for the coefficients of eq. 14

def _smul(a, b, K):
    return [sum((a[i] * b[k - i] for i in range(k + 1)), Fraction(0)) for k in range(K + 1)]


def _sinv(a, K):
    out = [Fraction(1) / a[0]]
    for k in range(1, K + 1):
        out.append(-sum((a[i] * out[k - i] for i in range(1, k + 1)), Fraction(0)) / a[0])
    return out


def _over(c, K):
    """Series of alpha/(alpha + c)."""
    return [Fraction(0)] + [Fraction((-1) ** (k - 1), 1) / Fraction(c) ** k for k in range(1, K + 1)]


def exact_series(K=4, term=1):
    """Coefficients of 2 nbar = sum c_k alpha^k at m = 0 from the closed form
    (I_0 and I_1 power series), termination factor ``term``."""
    t = Fraction(term)
    A = [Fraction(2) ** k / t**k / (math.factorial(k) ** 2) for k in range(K + 1)]
    B = [Fraction(2) ** k / t**k / (math.factorial(k) * math.factorial(k + 1)) for k in range(K + 1)]
    return _smul(A, _sinv(B, K), K)


def appendix_series(K=4, eq13_from=3):
    """Coefficients of 2 nbar = sum c_k alpha^k for the paper's construction (eqs. 12, 13)."""
    nt = K + 4
    one = [Fraction(1)] + [Fraction(0)] * K
    N = [None] * (nt + 1)
    N[0] = one
    if eq13_from == 3:
        N[2] = [Fraction(0), Fraction(1, 2)] + [Fraction(0)] * (K - 1)
        N[3] = _smul(N[2], _over(6, K), K)
        N[1] = [x + y for x, y in zip(one, [3 * v for v in _over(6, K)])]
    else:
        g = [6 * v for v in _smul(_over(2, K), _over(6, K), K)]     # 6 alpha/((alpha+2)(alpha+6))
        N[1] = _sinv([x - y for x, y in zip(one, g)], K)
        N[2] = _smul(N[1], _over(2, K), K)
        N[3] = _smul(N[2], _over(6, K), K)
    for n in range(4, nt + 1):
        N[n] = _smul(N[n - 1], _over(n * (n - 1), K), K)
    num = [sum((n * N[n][k] for n in range(nt + 1)), Fraction(0)) for k in range(K + 1)]
    den = [sum((N[n][k] for n in range(nt + 1)), Fraction(0)) for k in range(K + 1)]
    return [2 * v for v in _smul(num, _sinv(den, K), K)]


# ---------------------------------------------------------------------------- Part II

@dataclass(frozen=True)
class Stage:
    """The Stage I premises; each field is one reading of eqs. 21 to 30.

    Attributes
    ----------
    q : float
        Exponent of the area law a ~ (t - tau)^q (2/3: constant growth rate, sphere).
    pi_exp : float
        Exponent of 4 pi inside eq. 22 (1/2 as printed).
    ml_a : float or None
        Exponent in Gamma(1 + ml_a n) of eq. 29; None: q + 1.
    n26a : int or None
        Last term of eq. 26a; None: converged.
    capture : str
        ``"micelle"``: dN/dt = rho A_m/A (eq. 21a); ``"particle"``: rho A_p/A (a defect).
    ml_tol : float
        Size of the last term at which every Mittag-Leffler series is stopped.
    """
    q: float = 2.0 / 3.0
    pi_exp: float = 0.5
    ml_a: float = None
    n26a: int = None
    capture: str = "micelle"
    ml_tol: float = 1e-17

    @property
    def a(self):
        return self.q + 1.0 if self.ml_a is None else self.ml_a


STAGE = Stage()


def mittag_leffler(z, a, b=1.0, tol=1e-17, nmax=400, nterms=None):
    """E_(a,b)(z) = sum z^n/Gamma(a n + b) by its series.

    Returns (value, last term) with the series stopped when a term is below ``tol``
    (relative to 1), or after ``nterms`` + 1 terms when given.
    """
    z = np.asarray(z, dtype=float)
    s = np.zeros_like(z)
    zn = np.ones_like(z)
    last = np.inf
    for n in range(nmax + 1):
        t = zn * rgamma(a * n + b)
        s = s + t
        last = float(np.max(np.abs(t)))
        if nterms is not None:
            if n >= nterms:
                return s, last
        elif n > 2 and last < tol:
            return s, last
        zn = zn * z
    raise RuntimeError("Mittag-Leffler series did not converge")


def theta_c(st=STAGE):
    """theta/mu^(2/3) of eq. 23: [(4 pi)^pi_exp 3]^(2/3)."""
    return ((4.0 * math.pi) ** st.pi_exp * 3.0) ** (2.0 / 3.0)


def stage_closed(st=STAGE, find_roots=None, y_lo=-3.0):
    """Every Stage I constant in closed form.

    ``find_roots`` is ``pymrm.checks.find_roots`` (passed in so that model.py does not
    import the checks). Returns a dict; ``y1`` and what depends on it are nan when the
    nucleation flux has no sign change on [y_lo, 0] (then ``noroot`` is True).
    """
    q, a = st.q, st.a
    th = theta_c(st)
    out = {"c_y": th * gamma(q + 1.0), "c_theta": th ** (-1.0 / (q + 1.0)),
           "ts1_upper": (q + 1.0) ** (1.0 / (q + 1.0))}
    out["k_upper"] = out["c_theta"] * out["ts1_upper"]
    out["exp_rho"], out["exp_aS"] = q / (q + 1.0), 1.0 / (q + 1.0)
    out["exp_mu"] = -(2.0 / 3.0) / (q + 1.0)
    f0 = 1.0 if st.capture == "micelle" else 0.0
    flux = lambda y: f0 * float(mittag_leffler(y, a, tol=st.ml_tol)[0])
    out["noroot"] = False
    r = find_roots(flux, y_lo, 0.0, n_scan=64, xtol=1e-14)
    roots = r["roots"]
    if st.capture != "micelle":
        roots = []            # the flux is zero at every scan point: no particle is ever formed
    if len(roots) == 0:
        out["noroot"] = True
        out["flux_scan_max"] = max(abs(flux(y)) for y in np.linspace(y_lo, 0.0, 64))
        for k in ("y1", "bracket", "t1_const", "k_lower", "ts1_lower", "Ns_lower",
                  "ratio_t1", "ratio_k"):
            out[k] = float("nan")
        return out
    if len(roots) != 1:
        raise ValueError(f"nucleation flux has {len(roots)} roots on [{y_lo}, 0]")
    y1 = roots[0]
    out["y1"] = y1
    if st.n26a is None:
        out["bracket"] = float(mittag_leffler(y1, a, 2.0, tol=st.ml_tol)[0])
    else:
        out["bracket"] = float(mittag_leffler(y1, a, 2.0, nterms=st.n26a)[0])
    out["ts1_lower"] = (-y1 / gamma(q + 1.0)) ** (1.0 / (q + 1.0))
    out["t1_const"] = out["c_theta"] * out["ts1_lower"]
    out["k_lower"] = out["t1_const"] * out["bracket"]
    out["Ns_lower"] = out["ts1_lower"] * out["bracket"]
    out["ratio_t1"] = out["ts1_lower"] / out["ts1_upper"]
    out["ratio_k"] = out["k_lower"] / out["k_upper"]
    return out


def _dot(a, b):
    """Dot product without BLAS: a threaded BLAS ddot on long vectors costs more than it
    saves inside the step loop (measured: 20 times the CPU time at 25600 cells)."""
    return float(np.einsum("i,i->", a, b))


class _ReversedLU:
    """LU of P J P (P reverses the cells), natural ordering, no row pivoting: the dense
    inflow row is factorised last and fills nothing. Has the ``solve`` of a SuperLU."""

    def __init__(self, J):
        n = J.shape[0]
        self.rev = np.arange(n)[::-1]
        Jr = J[self.rev][:, self.rev]
        self.lu = splu(csc_array(Jr), permc_spec="NATURAL", diag_pivot_thresh=0.0)

    def solve(self, rhs):
        rhs = np.asarray(rhs, dtype=float).ravel()
        return self.lu.solve(rhs[self.rev])[self.rev]


class _RankOneLU:
    """(M + u w^T)^(-1) by the Sherman-Morrison formula on a sparse LU of M (COLAMD)."""

    def __init__(self, M, u, w):
        self.lu = splu(csc_array(M), permc_spec="COLAMD")
        self.w = w
        self.mu = self.lu.solve(u)
        self.den = 1.0 + w @ self.mu

    def solve(self, rhs):
        y = self.lu.solve(np.asarray(rhs, dtype=float).ravel())
        return y - self.mu * _dot(self.w, y) / self.den


class StageIPBE:
    """Stage I as a population balance in particle volume (scaled, E15 of the plan).

    dn/dt + dn/dv = 0 on v in [0, v_max] (v in units of mu t*, t in t*), inflow
    n(t, 0) = g: g = 1 (upper bound, eq. 21) or g = 1 - A_p (lower bound, eq. 27), with
    A_p = sum w_i n_i the particle area over a_s S and w_i the exact cell integral of
    v^q. First-order upwind (``construct_convflux_upwind``, the inflow an external
    boundary value through ``shapes_d``), ``construct_div`` (Cartesian, nu = 0),
    implicit Euler, each step solved by ``pymrm.newton`` with the analytic Jacobian
    (bidiagonal plus the one dense row of the nonlocal inflow), factorised once.

    Parameters
    ----------
    n_c : int
        Cells on [0, v_max].
    bound : str
        "upper" or "lower".
    q : float
        Area-law exponent (2/3).
    weights : str
        "exact" (cell integral of v^q) or "centre" (v_i^q dv, a defect).
    lag : bool
        Inflow from the previous step (a defect).
    capture : str
        "micelle" (g = 1 - A_p) or "particle" (g = A_p, a defect) for the lower bound.
    inflow_defect : float
        Inflow times (1 + inflow_defect dv^(1/2)) (a grid-dependent defect).
    factor : str
        "reversed": sparse LU of the Jacobian with the cells numbered from the right and
        natural ordering, so the dense inflow row comes last and causes no fill;
        "rank_one": sparse LU of the bidiagonal part (COLAMD) and the Sherman-Morrison
        formula for the dense row (the round-off companion).
    """

    def __init__(self, n_c, bound="lower", q=2.0 / 3.0, v_max=2.0, weights="exact", lag=False,
                 capture="micelle", inflow_defect=0.0, factor="reversed", inflow_scale=1.0):
        self.n_c, self.bound, self.lag, self.capture = n_c, bound, lag, capture
        self.dv = v_max / n_c
        self.dt = self.dv                       # CFL 1
        x_f = np.linspace(0.0, v_max, n_c + 1)
        self.x_f = x_f
        x_c = 0.5 * (x_f[1:] + x_f[:-1])
        shape = (n_c, 1)
        # left face: Dirichlet value n(t, 0) = g, passed as an external vector (d = 1 is
        # a coefficient); outward normal points to -v, b c = d with a = 0 gives n = g.
        # right face (v = v_max): pure outflow, never reached in the stage.
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"outflow": True})
        conv, conv_bc_l, _ = construct_convflux_upwind(shape, x_f, x_c, bc=bc, v=1.0,
                                                       shapes_d=((1, 1), None))
        div = construct_div(shape, x_f, nu=0)          # Cartesian coordinate in volume
        self.A = div @ conv
        self.B = div @ conv_bc_l                       # (n_c, 1): -1/dv in cell 0
        self.b = self.B.toarray().ravel()
        if weights == "exact":
            self.w = (x_f[1:] ** (q + 1.0) - x_f[:-1] ** (q + 1.0)) / (q + 1.0)
        elif weights == "centre":
            self.w = x_c**q * self.dv
        else:
            raise ValueError(weights)
        self.gfac = (1.0 + inflow_defect * math.sqrt(self.dv)) * inflow_scale
        M = eye_array(n_c, format="csc") / self.dt + self.A
        u = np.zeros(n_c)
        if bound == "lower" and not lag:
            u = (-1.0 if capture == "micelle" else 1.0) * self.gfac * self.b
        if factor == "reversed":
            J = M + csc_array(u[:, None]) @ csc_array(self.w[None, :])
            self.lu = _ReversedLU(csc_array(J))
        elif factor == "rank_one":
            self.lu = _RankOneLU(M, u, self.w)
        else:
            raise ValueError(factor)

    def inflow(self, n):
        if self.bound == "upper":
            return self.gfac
        ap = _dot(self.w, n)
        return self.gfac * ((1.0 - ap) if self.capture == "micelle" else ap)

    def march(self, t_end):
        """Implicit Euler to t_end; histories at t_k = k dt, initial state included.

        Returns dict of arrays t, Ap, N and the largest newton iteration count.
        """
        K = int(round(t_end / self.dt))
        n = np.zeros(self.n_c)
        t = np.arange(K + 1) * self.dt
        Ap = np.zeros(K + 1)
        N = np.zeros(K + 1)
        it = 0
        for k in range(1, K + 1):
            n_old = n
            g_lag = self.inflow(n_old)

            def fun(x):
                g = g_lag if (self.lag or self.bound == "upper") else self.inflow(x)
                res = (x - n_old) / self.dt + self.A @ x + self.b * g
                return res, self.lu

            sol = newton(fun, n_old, tol=1e-12, maxfev=10, solver="splu", stagnation=True)
            if not sol.success:
                raise RuntimeError(f"PBE step {k}: {sol.message}")
            it = max(it, sol.nit if hasattr(sol, "nit") else 0)
            n = np.asarray(sol.x).ravel()
            Ap[k] = _dot(self.w, n)
            N[k] = self.dv * np.sum(n)
        return {"t": t, "Ap": Ap, "N": N, "nit": it, "n": n}
