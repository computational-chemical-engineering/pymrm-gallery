"""Single-component Maxwell-Stefan diffusion in zeolite micropores (Krishna and Baur 2003).

The dual-site Langmuir isotherm with its thermodynamic correction factor
(eqs. 32 and 33), Karger's zero-loading diffusivities of MFI (eqs. 23 and 24),
and transient uptake into a spherical crystallite with a loading-dependent
Fick diffusivity (eqs. 25 to 28), solved with pymrm. Classes and functions
only: no plotting and no data loading.
"""

import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import splu

from pymrm import construct_boundary_value_matrices, construct_div, construct_grad, newton
from pymrm.checks import find_roots

GAMMA_FORMS = ("eq33", "no site B factor", "no site A factor")


class Isotherm:
    """Dual-site Langmuir isotherm, eq. (32), and its correction factor, eq. (33).

    ``Theta(P) = tA x/(1+x) + tB y/(1+y)`` with ``x = b_A P`` and ``y = b_B P``.
    With ``Th_sat_B = 0`` it is the single-site Langmuir of eq. (15).

    Parameters
    ----------
    b_A, b_B : float
        Site affinities, 1/Pa.
    Th_sat_A, Th_sat_B : float
        Site capacities, molecules per unit cell.
    form : str
        ``"eq33"``, or a defect for the break table: ``"no site B factor"``
        drops ``(1 - Theta_B/Theta_sat,B)`` from eq. (33), ``"no site A factor"``
        drops ``(1 - Theta_A/Theta_sat,A)``.
    """

    def __init__(self, b_A, Th_sat_A, b_B=0.0, Th_sat_B=0.0, form="eq33"):
        if form not in GAMMA_FORMS:
            raise ValueError(f"form {form!r}")
        self.b_A, self.tA, self.b_B, self.tB, self.form = b_A, Th_sat_A, b_B, Th_sat_B, form
        self.Th_sat = Th_sat_A + Th_sat_B

    def theta(self, P):
        """Loading Theta (molecules per unit cell) at pressure ``P`` (Pa), eq. (32)."""
        x, y = self.b_A * P, self.b_B * P
        return self.tA * x / (1 + x) + self.tB * y / (1 + y)

    def surface_potential(self, P):
        """``int_0^P Theta(p)/p dp = tA ln(1 + b_A P) + tB ln(1 + b_B P)``."""
        return self.tA * np.log1p(self.b_A * P) + self.tB * np.log1p(self.b_B * P)

    def P_of(self, Th):
        """Pressure at loading ``Th``: the positive root of a quadratic in P.

        Clearing the denominators of eq. (32) gives
        ``a2 P^2 + a1 P + Th = 0`` with ``a2 = b_A b_B (Th - Th_sat)`` and
        ``a1 = Th (b_A + b_B) - tA b_A - tB b_B``. For ``0 < Th < Th_sat``
        ``a2 < 0 < Th``, so exactly one root is positive; it is taken in the
        form that avoids cancellation.
        """
        Th = np.asarray(Th, float)
        a2 = self.b_A * self.b_B * (Th - self.Th_sat)
        a1 = Th * (self.b_A + self.b_B) - self.tA * self.b_A - self.tB * self.b_B
        sq = np.sqrt(a1 * a1 - 4.0 * a2 * Th)
        with np.errstate(divide="ignore", invalid="ignore"):
            big = np.where(a2 != 0.0, (a1 + sq) / (-2.0 * np.where(a2 != 0.0, a2, 1.0)), np.inf)
        return np.where(a1 < 0.0, 2.0 * Th / (sq - a1), big)

    def _gamma_P(self, P):
        """Eq. (33) at pressure ``P``, divided through by P so that P = 0 gives 1."""
        x, y = self.b_A * P, self.b_B * P
        num = self.tA * self.b_A / (1 + x) + self.tB * self.b_B / (1 + y)
        fa = 1.0 if self.form == "no site A factor" else 1.0 / (1 + x)
        fb = 1.0 if self.form == "no site B factor" else 1.0 / (1 + y)
        return num / (self.tA * self.b_A / (1 + x) * fa + self.tB * self.b_B / (1 + y) * fb)

    def Gamma(self, Th):
        """Thermodynamic correction factor at loading ``Th``, eq. (33)."""
        return self._gamma_P(self.P_of(Th))

    def extremum_function(self, u):
        """``Theta_u^2 - Theta Theta_uu`` in ``u = ln P``, zero at an extremum of Gamma.

        Gamma = Theta / Theta_u (eq. 33 is ``d ln P / d ln Theta``), so
        ``dGamma/du = (Theta_u^2 - Theta Theta_uu) / Theta_u^2``; the positive
        denominator is multiplied out. Positive where Gamma rises with loading.
        """
        P = np.exp(u)
        th, tu, tuu = 0.0, 0.0, 0.0
        for t, b in ((self.tA, self.b_A), (self.tB, self.b_B)):
            x = b * P
            th += t * x / (1 + x)
            tu += t * x / (1 + x) ** 2
            tuu += t * x * (1 - x) / (1 + x) ** 3
        return tu * tu - th * tuu

    def extrema(self, lo=1e-4, hi=0.9999, n_scan=2000):
        """Interior extrema of Gamma over ``lo Th_sat < Theta < hi Th_sat``.

        Every root of :meth:`extremum_function` on the whole window in
        ``ln P``, from :func:`pymrm.checks.find_roots`; a root where it falls
        through zero is a maximum of Gamma, one where it rises a minimum.

        Returns
        -------
        list of tuple
            ``(Theta, Gamma, kind)`` in order of loading, ``kind`` ``"max"`` or ``"min"``.
        """
        if self.form != "eq33":                 # the defect forms have no closed-form derivative
            return _extrema_by_derivative(self, lo, hi, n_scan)
        u_lo, u_hi = (float(np.log(self.P_of(f * self.Th_sat))) for f in (lo, hi))
        r = find_roots(self.extremum_function, u_lo, u_hi, n_scan=n_scan, xtol=1e-13)
        out = []
        for u in r["roots"]:
            d = 1e-6 * max(1.0, abs(u))
            kind = "max" if self.extremum_function(u - d) > 0 else "min"
            P = np.exp(u)
            out.append((float(self.theta(P)), float(self._gamma_P(P)), kind))
        return out


def _extrema_by_derivative(iso, lo, hi, n_scan):
    """Extrema of Gamma for a defect form: roots of a central difference of ln Gamma in u."""
    u_lo, u_hi = (float(np.log(iso.P_of(f * iso.Th_sat))) for f in (lo, hi))

    def dg(u, h=1e-5):
        return (np.log(iso._gamma_P(np.exp(u + h))) - np.log(iso._gamma_P(np.exp(u - h)))) / (2 * h)

    r = find_roots(dg, u_lo, u_hi, n_scan=n_scan, xtol=1e-11)
    out = []
    for u in r["roots"]:
        kind = "max" if dg(u - 1e-4) > 0 else "min"
        P = np.exp(u)
        out.append((float(iso.theta(P)), float(iso._gamma_P(P)), kind))
    return out


def first_max_min(ext):
    """``((Theta_max, Gamma_max), (Theta_min, Gamma_min))`` of the first of each kind; nan if absent."""
    mx = next(((t, g) for t, g, k in ext if k == "max"), (np.nan, np.nan))
    mn = next(((t, g) for t, g, k in ext if k == "min"), (np.nan, np.nan))
    return mx, mn


def gamma_numerical(iso, P, dlnP=1e-3):
    """Central difference of ln P against ln Theta at pressure ``P``, step ``dlnP`` in ln P.

    The definition eq. (33) encodes, ``Gamma = d ln f / d ln Theta`` with
    ``f = P``; its truncation error is of order ``dlnP^2``.
    """
    lo, hi = iso.theta(P * np.exp(-dlnP)), iso.theta(P * np.exp(dlnP))
    return 2 * dlnP / (np.log(hi) - np.log(lo))


def karger(nu_zz, nu_str, a, b, c, dz_denominator="corrected", dz_factor=1.0):
    """Karger's zero-loading diffusivities of MFI, eq. (23), and their average, eq. (24).

    Parameters
    ----------
    nu_zz, nu_str : float
        Jump frequencies along the zig-zag and straight channels, 1/s.
    a, b, c : float
        Unit-cell dimensions, m.
    dz_denominator : str
        ``"corrected"`` (``nu_str + nu_zz``) or ``"as printed"`` (``nu_zz + nu_zz``).
    dz_factor : float
        Multiplies D_z (``0.5``: the halved relation the review's worked
        examples use).

    Returns
    -------
    tuple of float
        ``(D_x, D_y, D_z, D_avg)`` in m^2/s.
    """
    den = nu_str + nu_zz if dz_denominator == "corrected" else 2.0 * nu_zz
    Dx = 0.25 * nu_zz * a ** 2
    Dy = 0.25 * nu_str * b ** 2
    Dz = dz_factor * 0.25 * (nu_str * nu_zz / den) * c ** 2
    return Dx, Dy, Dz, (Dx + Dy + Dz) / 3.0


def series_one_minus_F(fo, n_terms=800):
    """``1 - F`` of eq. (30), summed directly so that no cancellation occurs."""
    fo = np.atleast_1d(np.asarray(fo, float))
    m = np.arange(1, n_terms + 1)[:, None]
    return 6.0 / np.pi ** 2 * np.sum(np.exp(-m ** 2 * np.pi ** 2 * fo) / m ** 2, axis=0)


def series_F(fo, n_terms=800):
    """Fractional approach to equilibrium at constant D, eq. (30)."""
    return 1.0 - series_one_minus_F(fo, n_terms)


def series_tail(fo, n_terms=800):
    """Bound on the tail of eq. (30) after ``n_terms`` terms, relative to ``1 - F``.

    The terms fall with a ratio that itself falls, so the tail is at most
    ``t_{N+1} / (1 - t_{N+2}/t_{N+1})``.
    """
    k1, k2 = n_terms + 1, n_terms + 2
    t1 = np.exp(-k1 ** 2 * np.pi ** 2 * fo) / k1 ** 2
    t2 = np.exp(-k2 ** 2 * np.pi ** 2 * fo) / k2 ** 2
    tail = t1 / (1.0 - t2 / t1) if t1 > 0.0 else 0.0
    return float(6.0 / np.pi ** 2 * tail / series_one_minus_F(fo, n_terms)[0])


def sherwood(fo, one_minus_F, prefactor=2.0 / 3.0):
    """Time-averaged Sherwood number, eq. (31): ``-prefactor ln(1 - F) / Fo``."""
    return -prefactor * np.log(np.asarray(one_minus_F)) / np.asarray(fo)


def k_instantaneous(fo, n_terms=800):
    """Instantaneous LDF constant ``(dF/dFo) / (1 - F)`` of eq. (30), in units of D/r_c^2."""
    m = np.arange(1, n_terms + 1)
    e = np.exp(-m ** 2 * np.pi ** 2 * fo)
    return float(6.0 * np.sum(e) / (6.0 / np.pi ** 2 * np.sum(e / m ** 2)))


def one_root(f, lo, hi, what, **kw):
    """The single root of ``f`` on ``[lo, hi]`` by :func:`pymrm.checks.find_roots`; raises otherwise."""
    r = find_roots(f, lo, hi, **kw)
    if len(r["roots"]) != 1:
        raise RuntimeError(f"{what}: {len(r['roots'])} roots on [{lo:.3g}, {hi:.3g}]")
    return r["roots"][0]


class Crystal:
    """Transient diffusion in a spherical zeolite crystallite, eqs. (25) to (27).

    Dimensionless: radius 1, time the Fourier number ``Fo = D_MS(0) t / r_c^2``.
    ``dTheta/dFo = -div N`` with ``N = -D(Theta) dTheta/dr``, symmetry at the
    centre and ``Theta = Theta_s`` at the surface, by backward Euler on a step
    schedule ``dt = min(dt_max, max(dt0, dt_frac Fo))`` truncated at each output
    time, all three scaled by ``s``. A constant ``D`` makes the problem linear
    (one sparse LU per step size); a loading-dependent ``D`` is solved by
    :func:`pymrm.newton` per step with the assembled Jacobian.

    Parameters
    ----------
    D : float or callable
        The Fick diffusivity over D_MS(0): a constant, or ``D(Theta)``.
    n_r : int
        Number of cells.
    nu : int
        Geometry index of ``construct_div`` (2 sphere, 1 cylinder, 0 slab); it
        also sets the volume weights of the mean.
    centre : str
        ``"symmetry"``, or ``"dirichlet"`` (a defect: Theta held at its
        initial value at r = 0).
    mean : str
        Interior face diffusivity: ``"harmonic"`` or ``"arithmetic"`` mean of
        the two cell values. Boundary faces take D at the boundary loading from
        ``construct_boundary_value_matrices``.
    tol : float
        Newton step tolerance on Theta.
    upper : float or None
        Loading at which ``D`` diverges (saturation); Newton iterates are kept
        between zero and it through ``bounds``, where ``D`` is defined. A
        safeguard: the page's runs never reach either bound.
    """

    def __init__(self, D, n_r=150, nu=2, centre="symmetry", mean="harmonic", tol=1e-11, upper=None):
        self.D, self.n_r, self.nu, self.centre, self.mean, self.tol = D, n_r, nu, centre, mean, tol
        self.bounds = None if upper is None else (0.0, upper)
        self.shape = (n_r, 1)                               # (space, field)
        self.r_f = np.linspace(0.0, 1.0, n_r + 1)
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        self.div = construct_div(self.shape, self.r_f, nu=nu)
        # exact cell volumes over the domain volume: they sum to one
        self.w = np.diff(self.r_f ** (nu + 1))
        self.eye = sps.eye(n_r, format="csc")
        # left and right cell of each face; boundary faces (rows 0 and n_r) have none
        self.L_face = sps.eye(n_r + 1, n_r, k=-1, format="csr")
        self.R_face = sps.eye(n_r + 1, n_r, format="csr")
        self.L_face[n_r, n_r - 1] = 0.0
        self.R_face[0, 0] = 0.0
        self.L_face.eliminate_zeros()
        self.R_face.eliminate_zeros()

    def _operators(self, Th0, Th_s):
        # outward normal n, a dTheta/dn + b Theta = d:
        #   centre  r = 0 (n = -r): symmetry, dTheta/dr = 0   -> a = 1, b = 0, d = 0
        #           ('dirichlet' defect: Theta = Theta_0      -> a = 0, b = 1, d = Theta_0)
        #   surface r = 1 (n = +r): Theta = Theta_s, eq. (27) -> a = 0, b = 1, d = Theta_s
        left = ({"a": 1.0, "b": 0.0, "d": 0.0} if self.centre == "symmetry"
                else {"a": 0.0, "b": 1.0, "d": float(Th0)})
        bc = (left, {"a": 0.0, "b": 1.0, "d": float(Th_s)})
        grad, gb = construct_grad(self.shape, self.r_f, self.r_c, bc)
        bl, blv = construct_boundary_value_matrices(self.shape, self.r_f, self.r_c, bc[0], bound_id=0)
        br, brv = construct_boundary_value_matrices(self.shape, self.r_f, self.r_c, bc[1], bound_id=1)
        self.grad, self.grad_bc = grad.tocsr(), np.asarray(gb.toarray()).reshape(-1, 1)
        zero = sps.csr_matrix((self.n_r - 1, self.n_r))
        # boundary face loadings: Theta_face = B @ Theta + b
        self.Bface = sps.vstack([bl, zero, br]).tocsr()
        bv = np.zeros((self.n_r + 1, 1))
        bv[0] = np.asarray(blv.toarray() if sps.issparse(blv) else blv).ravel()[0]
        bv[-1] = np.asarray(brv.toarray() if sps.issparse(brv) else brv).ravel()[0]
        self.bface = bv

    def _face_D(self, Th):
        """Face diffusivities, shape ``(n_r + 1, 1)``, and their Jacobian with respect to Theta."""
        Dc = self.D(Th)
        h = 1e-7 * (1.0 + np.abs(Th))
        dDc = (self.D(Th + h) - self.D(Th - h)) / (2 * h)
        Dl, Dr = Dc[:-1], Dc[1:]
        if self.mean == "harmonic":
            Di = 2 * Dl * Dr / (Dl + Dr)
            wl, wr = 2 * Dr ** 2 / (Dl + Dr) ** 2, 2 * Dl ** 2 / (Dl + Dr) ** 2
        else:
            Di, wl, wr = 0.5 * (Dl + Dr), 0.5 * np.ones_like(Dl), 0.5 * np.ones_like(Dr)
        thb = (self.Bface @ Th + self.bface)[[0, -1]]       # the two boundary-face loadings
        hb = 1e-7 * (1.0 + np.abs(thb))
        Db = self.D(thb)
        dDb = (self.D(thb + hb) - self.D(thb - hb)) / (2 * hb)
        z = np.zeros((1, 1))
        Df = np.vstack([Db[:1], Di, Db[1:]])
        dDf = (self.L_face.multiply(np.vstack([z, wl * dDc[:-1], z]))
               + self.R_face.multiply(np.vstack([z, wr * dDc[1:], z]))
               + self.Bface.multiply(np.vstack([dDb[:1], np.zeros((self.n_r - 1, 1)), dDb[1:]])))
        return Df, dDf

    def _residual(self, Th_old, dt):
        def fun(y):
            y = np.reshape(y, self.shape)
            g = self.grad @ y + self.grad_bc
            Df, dDf = self._face_D(y)
            r = (y - Th_old) / dt + self.div @ (-Df * g)
            dflux = -(self.grad.multiply(Df) + dDf.multiply(g))
            return r, (self.eye / dt + self.div @ dflux).tocsc()
        return fun

    def run(self, Th0, Th_s, fo_out, s=1.0, dt0=5e-6, dt_frac=0.02, dt_max=5e-3, keep_rates=False):
        """Mean loading at each Fourier number of ``fo_out``, from ``Theta = Th0`` everywhere.

        Raises unless every Newton solve reports success.

        Returns
        -------
        means : ndarray
            Volume-mean loading at each output time.
        profiles : dict
            Cell loadings at each output time.
        """
        self._operators(Th0, Th_s)
        Th = np.full(self.shape, float(Th0))
        fo, out, profiles = 0.0, [], {}
        self.n_steps, self.max_iter, self.floor = 0, 0, 0.0
        self.rates, self.step_fo = [], []
        linear = not callable(self.D)
        if linear:
            lap = (self.div @ (-self.D * self.grad)).tocsc()
            src = self.div @ (-self.D * self.grad_bc)
            norm_lap = float(abs(lap).sum(axis=1).max())
        eps = np.finfo(float).eps
        for target in np.atleast_1d(np.asarray(fo_out, float)):
            while fo < target - 1e-15:
                dt = min(dt_max * s, max(dt0 * s, dt_frac * s * fo), target - fo)
                if linear:
                    lu = splu((self.eye / dt + lap).tocsc())
                    new = lu.solve(Th / dt - src).reshape(self.shape)
                    # backward error of LU on an M-matrix: at most cond eps, cond <= 1 + dt |L|
                    self.floor += (1.0 + dt * norm_lap) * eps * max(1.0, float(np.abs(new).max()))
                else:
                    res = newton(self._residual(Th, dt), Th, tol=self.tol, stagnation=True, maxfev=50,
                             bounds=self.bounds)
                    if not res.success:
                        raise RuntimeError(f"crystal step at Fo = {fo:.3e} failed: {res.message}")
                    new = np.reshape(res.x, self.shape)
                    self.max_iter = max(self.max_iter, int(res.nit))
                    self.floor += self.tol
                if keep_rates:
                    self.rates.append(float(self.w @ (new - Th).ravel()) / dt)
                    self.step_fo.append(fo + dt)
                Th = new
                fo += dt
                self.n_steps += 1
            out.append(float(self.w @ Th.ravel()))
            profiles[float(target)] = Th.ravel().copy()
        return np.array(out), profiles

    def fo_at_mean(self, Th0, Th_s, target, fo_max, s=1.0, dt0=5e-6, dt_frac=0.02, dt_max=5e-3):
        """Fourier number at which the mean loading reaches ``target``, on the solved model.

        The schedule is marched until a step crosses ``target``; the crossing
        is then the root, in the size of that last (truncated) step, of the
        mean after one backward-Euler step: exactly the value :meth:`run`
        returns for an output time inside that step.
        """
        self._operators(Th0, Th_s)
        Th = np.full(self.shape, float(Th0))
        fo = 0.0
        sign = np.sign(target - Th0)

        def step(Th, dt):
            if not callable(self.D):
                lap = (self.div @ (-self.D * self.grad)).tocsc()
                src = self.div @ (-self.D * self.grad_bc)
                return splu((self.eye / dt + lap).tocsc()).solve(Th / dt - src).reshape(self.shape)
            res = newton(self._residual(Th, dt), Th, tol=self.tol, stagnation=True, maxfev=50,
                             bounds=self.bounds)
            if not res.success:
                raise RuntimeError(f"crystal step failed: {res.message}")
            return np.reshape(res.x, self.shape)

        while fo < fo_max:
            dt = min(dt_max * s, max(dt0 * s, dt_frac * s * fo))
            new = step(Th, dt)
            if sign * (float(self.w @ new.ravel()) - target) >= 0.0:
                tau = one_root(lambda t: float(self.w @ step(Th, t).ravel()) - target,
                               1e-9 * dt, dt, "mean crossing", n_scan=8, xtol=1e-15)
                return fo + tau
            Th, fo = new, fo + dt
        raise RuntimeError(f"mean did not reach {target} by Fo = {fo_max}")
