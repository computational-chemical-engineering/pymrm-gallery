"""Drop breakage and coalescence in a continuous-flow stirred vessel.

Coulaloglou and Tavlarides, Chem. Eng. Sci. 32, 1289-1297 (1977). Units as in the
paper (cgs): v in cm3, d and D in cm, N* in rev/s, every frequency in 1/s. The
physical properties are not printed; the model works with the groups
S = sigma/rho_d (cm3/s2) and X = mu_c rho_c/sigma^2 (s3/cm4).

Closed forms
    ``breakage_frequency`` (eqs. 14, 36), ``breakage_group`` (the exponent y),
    ``collision_frequency`` (eqs. 19, 37), ``efficiency`` (eqs. 27, 37),
    ``d_max`` (the root of g theta = 1, breakage alone), ``contour_exponent``,
    ``x_from_efficiency`` (X from the printed 6 % efficiency) and the paper's own
    estimates of pp. 1295-1296 (``estimates``).
Population balance
    ``SteadyPBE``: the steady balance (eqs. 3, 4) on fixed pivots (Kumar and
    Ramkrishna), dense residual and analytic Jacobian, solved by ``pymrm.newton``
    from the breakage-only state; ``breakage_only`` is the exact back substitution of
    the balance without coalescence; ``march`` is the paper's marching iteration.

No plotting and no data loading here.
"""
import math
from functools import lru_cache
from dataclasses import dataclass, replace

import numpy as np

from pymrm import newton
from pymrm.checks import find_roots

SIX_PI = 6.0 / math.pi


@dataclass(frozen=True)
class Params:
    """Printed constants and the page's readings; every field is a break-row knob.

    C1, C2 (eq. 36), C3 in cm^-3 per vessel and C4 in cm^-2 (eq. 37); D, V_T, theta from
    p. 1294. ``damp`` multiplies every power of (1 + phi) of eqs. 36 and 37 (1: as
    printed, 0: no damping). ``shape_power`` is the exponent of
    v^(1/3) v'^(1/3)/(v^(1/3) + v'^(1/3)) in the efficiency (4: eq. 27). ``beta`` is
    'normalised' (on [0, v']) or 'printed' (the printed 2.4). ``assign`` is 'pivot'
    (fixed pivot, two neighbours) or 'nearest' (whole product to the nearest pivot).
    ``per_volume`` reads the counts per cm3 (C3 / V_T).
    """
    C1: float = 0.40
    C2: float = 0.08
    C3: float = 2.8e-6
    C4: float = 1.83e9
    D: float = 10.0
    VT: float = 12000.0
    theta: float = 600.0
    damp: float = 1.0
    shape_power: float = 4.0
    beta: str = "normalised"
    assign: str = "pivot"
    per_volume: bool = False


PRINTED = Params()


def volume(d):
    """Sphere volume of diameter d."""
    return math.pi / 6.0 * np.asarray(d) ** 3


def diameter(v):
    """Sphere diameter of volume v."""
    return np.cbrt(SIX_PI * np.asarray(v))


def breakage_group(v, n, phi, s, p=PRINTED):
    """Exponent y of eq. 36: C2 S (1+phi)^2/(v^(5/9) D^(4/3) N*^2)."""
    return p.C2 * s * (1 + phi) ** (2 * p.damp) / (np.asarray(v) ** (5 / 9) * p.D ** (4 / 3) * n**2)


def breakage_frequency(v, n, phi, s, p=PRINTED):
    """g(v) of eq. 36, 1/s."""
    v = np.asarray(v, dtype=float)
    return (p.C1 * v ** (-2 / 9) * p.D ** (2 / 3) * n / (1 + phi) ** p.damp
            * np.exp(-breakage_group(v, n, phi, s, p)))


def collision_frequency(v, w, n, phi, p=PRINTED):
    """h(v, v') of eq. 37 per pair of drops, 1/s (counts per vessel, or per cm3)."""
    c3 = p.C3 / p.VT if p.per_volume else p.C3
    return (c3 * (v ** (2 / 3) + w ** (2 / 3)) * np.sqrt(v ** (2 / 9) + w ** (2 / 9))
            * p.D ** (2 / 3) * n / (1 + phi) ** p.damp)


def shape_factor(v, w, p=PRINTED):
    """(v^(1/3) w^(1/3)/(v^(1/3) + w^(1/3)))^shape_power."""
    a, b = np.cbrt(v), np.cbrt(w)
    return (a * b / (a + b)) ** p.shape_power


def efficiency(v, w, n, phi, x, p=PRINTED):
    """lambda(v, v') of eq. 37, with X = mu_c rho_c/sigma^2."""
    return np.exp(-p.C4 * x * p.D**2 * n**3 / (1 + phi) ** (3 * p.damp) * shape_factor(v, w, p))


def x_from_efficiency(lam=0.06, v=1.88e-5, w=3.75e-5, n=4.17, phi=0.10, p=PRINTED):
    """X that gives efficiency ``lam`` for the drop pair of p. 1295."""
    return -math.log(lam) / (p.C4 * p.D**2 * n**3 / (1 + phi) ** (3 * p.damp) * shape_factor(v, w, p))


def d_max(n, phi, s, p=PRINTED):
    """Diameter (cm) where g theta = 1 (breakage alone), by find_roots in ln v."""
    def f(lv):
        return math.log(float(breakage_frequency(math.exp(lv), n, phi, s, p)) * p.theta)
    r = find_roots(f, math.log(1e-10), math.log(1e-2), n_scan=64, xtol=1e-13)["roots"]
    if len(r) != 1:
        raise ValueError(f"g theta = 1 has {len(r)} roots")
    return float(diameter(math.exp(r[0])))


def contour_exponent(y):
    """d ln d / d ln N* along g = const (eq. 36): -3 (1 + 2y)/(5y - 2)."""
    return -3.0 * (1 + 2 * y) / (5 * y - 2)


def estimates(rotta_in_c3=False, contact="d_over_u", c2="pi/12", cutter=70.0):
    """The paper's estimates of pp. 1295-1296 under a named reading of the text.

    Returns C4 at h = 1000 and 500 Angstrom (cm^-2), C3 per vessel, C2/c1 and C1 c5.
    ``contact`` 'd_over_u' takes tbar = (d+d')/u with the Rotta velocity; 'unity' takes
    the constant of eq. 26 as one. ``c2`` is the kinetic-energy factor of eq. 9.
    """
    rotta, eps_coef, vt = 1.88, 0.407, 12000.0
    base = 3.0 / (16.0 * math.pi) * eps_coef * SIX_PI ** (4 / 3)
    k4 = rotta**1.5 if contact == "d_over_u" else rotta
    c4 = {h: base * k4 / (h * 1e-8) ** 2 for h in (1000.0, 500.0)}
    kin = (9 * math.pi / 2) ** (1 / 3) * eps_coef ** (1 / 3)
    c3 = kin / vt * (math.sqrt(rotta * SIX_PI ** (2 / 9)) if rotta_in_c3 else 1.0)
    eps_imp = cutter * eps_coef
    c2f = math.pi / 12 if c2 == "pi/12" else 0.5
    c2_over_c1 = 1.0 / (c2f * rotta * eps_imp ** (2 / 3) * SIX_PI ** (5 / 9))
    c1_c5 = eps_imp ** (1 / 3) * SIX_PI ** (-2 / 9)
    return {"C4_h1000": c4[1000.0], "C4_h500": c4[500.0], "C3": c3,
            "C2_over_c1": c2_over_c1, "C1_c5": c1_c5}


def _gauss_hat(a, b, c, w, mu, sd):
    """int_a^b (v - c)/w exp(-(v - mu)^2/(2 sd^2)) dv in closed form."""
    if b <= a:
        return 0.0
    r2 = math.sqrt(2.0) * sd
    i0 = sd * math.sqrt(math.pi / 2) * (math.erf((b - mu) / r2) - math.erf((a - mu) / r2))
    i1 = sd * sd * (math.exp(-(a - mu) ** 2 / (2 * sd * sd)) - math.exp(-(b - mu) ** 2 / (2 * sd * sd)))
    return (i1 + (mu - c) * i0) / w


@lru_cache(maxsize=16)
def _breakage_matrix(q, x1, nx, beta):
    """Fixed-pivot breakage matrix B[i, k]: daughters at pivot i per breakage of pivot k (cached)."""
    xs = x1 * 2.0 ** (np.arange(nx) / q)
    xe = np.concatenate(([0.0], xs))
    B = np.zeros((nx, nx))
    for k in range(1, nx):
        mu, sd = xs[k] / 2, xs[k] / 6
        if beta == "printed":
            pref = 2 * 2.4 / xs[k]                   # nu = 2 times the printed 2.4/v'
        else:
            pref = 2 / (2 * sd * math.sqrt(math.pi / 2) * math.erf(3 / math.sqrt(2)))
        for i in range(k + 1):
            val = _gauss_hat(xe[i], xs[i], xe[i], xs[i] - xe[i], mu, sd)
            if i < k:
                val -= _gauss_hat(xs[i], xs[i + 1], xs[i + 1], xs[i + 1] - xs[i], mu, sd)
            B[i, k] = pref * val
    B.flags.writeable = False
    return B


class SteadyPBE:
    """Steady population balance of eqs. 3 and 4 on fixed pivots.

    Pivots x_k = x1 2^(k/q), k = 0 .. octaves q. Unknowns y_k: drops per vessel at
    pivot k. Residual R = feed + B (g y) + births - (g + 1/theta) y - y (K y), with
    K = lambda h (pairs whose volume exceeds the last pivot excluded from gain and loss).

    Parameters
    ----------
    q : int
        Pivots per doubling of volume.
    n, phi, s, x : float
        N* (rev/s), holdup, S = sigma/rho_d (cm3/s2), X = mu_c rho_c/sigma^2.
    p : Params
        Printed constants and readings.
    x1, octaves : float, int
        Smallest pivot (cm3) and the number of volume doublings of the grid.
    d0, rel0 : float
        Feed: Gaussian in v, mean diameter d0 (cm), standard deviation rel0 times the mean volume.
    coalescence : bool
        False drops K (breakage and throughflow only).
    """

    def __init__(self, q, n, phi, s, x, p=PRINTED, x1=1e-8, octaves=17, d0=0.1, rel0=0.1,
                 coalescence=True):
        self.p, self.phi, self.q = p, phi, q
        nx = octaves * q + 1
        xs = x1 * 2.0 ** (np.arange(nx) / q)
        self.n, self.x = nx, xs
        g = breakage_frequency(xs, n, phi, s, p)
        g[0] = 0.0                       # the smallest pivot has no pivot to break into
        self.g = g
        self.B = self._breakage_matrix()
        q_flow = p.VT / p.theta
        v0 = float(volume(d0))
        w = np.exp(-((xs - v0) ** 2) / (2 * (rel0 * v0) ** 2))
        self.feed = q_flow * phi * w / np.sum(xs * w)
        I, J = np.triu_indices(nx)
        vs = xs[I] + xs[J]
        ok = vs <= xs[-1] * (1 + 1e-13)
        I, J, vs = I[ok], J[ok], vs[ok]
        k = np.clip(np.searchsorted(xs, vs, side="right") - 1, 0, nx - 1)
        kp = np.minimum(k + 1, nx - 1)
        al = np.where(kp > k, (xs[kp] - vs) / np.where(kp > k, xs[kp] - xs[k], 1.0), 1.0)
        if p.assign == "nearest":       # whole product to the nearest pivot (break-row reading)
            near = np.abs(xs[kp] - vs) < np.abs(vs - xs[k])
            al = np.where(near, 0.0, 1.0)
        K = np.zeros((nx, nx))
        if coalescence:
            K[I, J] = (collision_frequency(xs[I], xs[J], n, phi, p)
                       * efficiency(xs[I], xs[J], n, phi, x, p))
            K[J, I] = K[I, J]
        self.K = K
        self.I, self.J, self.k, self.kp, self.al = I, J, k, kp, al
        self.half = np.where(I == J, 0.5, 1.0)

    def _breakage_matrix(self):
        return _breakage_matrix(self.q, float(self.x[0]), self.n, self.p.beta)

    def rates(self, y):
        """Coalescence events per pair (1/s), coalescence births, breakage births, coalescence loss."""
        R = self.half * self.K[self.I, self.J] * y[self.I] * y[self.J]
        cb = (np.bincount(self.k, self.al * R, self.n)
              + np.bincount(self.kp, (1 - self.al) * R, self.n))
        return R, cb, self.B @ (self.g * y), y * (self.K @ y)

    def residual(self, y):
        """Steady residual, drops per second per pivot."""
        R, cb, br, loss = self.rates(y)
        return self.feed + br + cb - (self.g + 1 / self.p.theta) * y - loss

    def jacobian(self, y):
        """Analytic dense Jacobian of ``residual``."""
        jac = self.B * self.g[None, :] - np.diag(self.g + 1 / self.p.theta + self.K @ y) - y[:, None] * self.K
        kij = self.half * self.K[self.I, self.J]
        dI, dJ = kij * y[self.J], kij * y[self.I]
        for tgt, wt in ((self.k, self.al), (self.kp, 1 - self.al)):
            np.add.at(jac, (tgt, self.I), wt * dI)
            np.add.at(jac, (tgt, self.J), wt * dJ)
        return jac

    def fun(self, y):
        """(residual, Jacobian) for pymrm.newton."""
        return self.residual(y), self.jacobian(y)

    def breakage_only(self):
        """Exact steady state without coalescence: back substitution from the largest pivot."""
        y = np.zeros(self.n)
        for k in range(self.n - 1, -1, -1):
            # B[k, k] g_k y_k: daughters of pivot k that land on pivot k itself (fixed pivot)
            y[k] = ((self.feed[k] + self.B[k, k + 1:] @ (self.g[k + 1:] * y[k + 1:]))
                    / (self.g[k] * (1 - self.B[k, k]) + 1 / self.p.theta))
        return y

    def solve(self, y0=None, tol=1e-12, maxfev=50):
        """Steady state by ``pymrm.newton`` (analytic dense Jacobian, bounds at zero,
        scaled by the largest unknown), started from the breakage-only state.

        Returns (y, info) with info = {'newton': OptimizeResult, 'residual'}, the
        residual relative to the largest feed entry. Raises RuntimeError when newton
        does not report success.
        """
        y = self.breakage_only() if y0 is None else np.array(y0, dtype=float)
        res = newton(self.fun, y, tol=tol, maxfev=maxfev, bounds=(0.0, None),
                     scale=max(float(np.max(y)), 1.0), stagnation=True)
        if not res.success:
            raise RuntimeError(f"newton did not converge: {res.message}")
        y = res.x
        return y, {"newton": res,
                   "residual": float(np.max(np.abs(self.residual(y))) / np.max(self.feed))}

    def march(self, rtol=1e-12, sweeps=500, increment=None, floor=1e-6):
        """The paper's marching iteration (p. 1293): sweeps from the largest pivot down,
        coalescence births and losses lagged by one sweep. Stops when the largest change
        of a sweep relative to the largest entry is at most ``rtol``, or, with
        ``increment`` given (the paper's rule, 0.001), when the relative change per
        increment is at most ``increment`` at every pivot holding more than ``floor``
        of the largest entry. Returns (y, sweeps used) or (y, -1) when it has not
        stopped after ``sweeps``."""
        y = self.breakage_only()
        for it in range(1, sweeps + 1):
            _, cb, _, _ = self.rates(y)
            ky = self.K @ y
            new = np.zeros(self.n)
            for k in range(self.n - 1, -1, -1):
                new[k] = ((self.feed[k] + cb[k] + self.B[k, k + 1:] @ (self.g[k + 1:] * new[k + 1:]))
                          / (self.g[k] + 1 / self.p.theta + ky[k]))
            if increment is None:
                done = np.max(np.abs(new - y)) <= rtol * np.max(new)
            else:
                sig = new > floor * np.max(new)
                done = np.max(np.abs(new - y)[sig] / new[sig]) <= increment
            y = new
            if done:
                return y, it
        return y, -1

    def stats(self, y):
        """d32 (cm), omega_v (eq. 34, 1/s), omega_n (eq. 33, 1/s), N (drops) and the volume ratio."""
        d = diameter(self.x)
        R, cb, br, loss = self.rates(y)
        return {"d32": float(np.sum(y * d**3) / np.sum(y * d**2)),
                "omega_v": float(np.sum(self.x * cb) / (self.p.VT * self.phi)),
                "omega_n": float(2 * np.sum(R) / np.sum(y)),
                "N": float(np.sum(y)),
                "vol": float(np.sum(self.x * y) / (self.p.VT * self.phi))}

    def number_deficit(self, y):
        """Daughters lost below the smallest pivot: per daughter formed and per fed drop."""
        lost = np.sum(self.g * y * (2.0 - self.B.sum(axis=0)))
        return float(lost / np.sum(2 * self.g * y)), float(lost / np.sum(self.feed))


def constant_kernel_number(K, n0, theta):
    """E20: steady number in a CSTR with a constant kernel and no breakage."""
    return (-1 / theta + math.sqrt(1 / theta**2 + 2 * K * n0)) / K
