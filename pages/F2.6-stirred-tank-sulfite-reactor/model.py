"""A stirred gas-liquid tank with absorption followed by liquid reaction.

Film theory coupled to a well-mixed liquid bulk (S1 + S7). The gas is uniform at
the interface value C*; dissolved gas reacts in the film and in the bulk at the
pseudo-m-th-order rate k_m c^m (sulfite and catalyst in excess). In
dimensionless form, xi = x / delta and u = c / C*::

    u'' = Ha^2 u^m,   u(0) = 1,   u(1) = u_b,
    g := -u'(1) = Ha^2 (Al - 1) u_b^m            (bulk balance per liquid volume)

with Ha^2 = k_m C*^(m-1) D / k_L^2 and Al = k_L^2 / (D k_L A), the Hinterland
ratio. phi = -u'(0) is the absorption rate relative to k_L A C*: what the
sulfite method reports as its k_L A when it takes the bulk as nil. None of this
is printed in van 't Riet (1979); it is the page's model of the reactor the
review describes.

Two routes. :func:`closed_m1` is the exact solution for m = 1. :class:`FilmTank`
solves any order with pymrm: a bordered Newton system whose last unknown is the
bulk level u_b, entering the film through the right Dirichlet face of a constant
gradient operator (``shapes_d``). The study functions below take either route as
a callable ``solve(Ha2, Al, m) -> dict``.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from pymrm import construct_div, construct_grad, newton
from pymrm.checks import find_roots

__all__ = ["kla_eq9", "hinterland", "closed_m1", "FilmTank", "PymrmRoute", "edge", "band",
           "point_where", "al_root", "box_share", "richardson", "SulfiteTank",
           "sulfite_numbers", "NS_LIMIT", "mapped_grid"]


# codehealth: allow params -- physical inputs and the defects the page injects
def kla_eq9(pv, vs, K, alpha, beta, vs_unit=1.0):
    """Eq. 9 of van 't Riet (1979): k_L A = K (P/V)^alpha v_s^beta, in 1/s.

    Parameters
    ----------
    pv : float or array
        Power per unit volume, W/m3.
    vs : float or array
        Superficial gas velocity, m/s.
    K, alpha, beta : float
        The printed constants.
    vs_unit : float
        Factor applied to ``vs`` before the power law (1 for m/s; 100 is the
        cm/s misreading the break table injects).
    """
    return K * np.asarray(pv, float) ** alpha * (np.asarray(vs, float) * vs_unit) ** beta


def hinterland(kL, kla, D):
    """Al = k_L^2 / (D k_L A): liquid volume per film volume, k_L / (A D)."""
    return kL ** 2 / (D * kla)


def _power(u, m):
    """u^m and its derivative, sign-safe for an iterate that dips below zero."""
    a = np.abs(u)
    return u * a ** (m - 1.0), m * a ** (m - 1.0)


# codehealth: allow params -- physical inputs and the defects the page injects
def closed_m1(Ha2, Al, m=1, *, film_excluded=True, bulk="coupled", share="flux"):
    """Exact solution for a first-order reaction (m must be 1).

    Returns the same fields as :meth:`FilmTank.solve`: ``phi``, ``ub``, ``g``,
    ``film_share`` = 1 - g/phi and ``E`` = phi/(1 - u_b). ``film_excluded=False``
    writes the bulk balance with Al for Al - 1; ``bulk="dirichlet0"`` sets
    u_b = 0 (the film alone, Al -> infinity); ``share="ub"`` writes the film share
    as 1 - u_b. The last three are defects the page injects.
    """
    if m != 1:
        raise ValueError("the closed form holds for m = 1 only")
    Ha = np.sqrt(Ha2)
    c, s = np.cosh(Ha), np.sinh(Ha)
    vol = Al - 1.0 if film_excluded else Al
    ub = 0.0 if bulk == "dirichlet0" else 1.0 / (c + Ha * vol * s)
    phi = Ha * (c - ub) / s
    g = Ha * (1.0 - ub * c) / s
    fs = 1.0 - g / phi if share == "flux" else 1.0 - ub
    return dict(phi=phi, ub=ub, g=g, film_share=fs, E=phi / (1.0 - ub))


class FilmTank:
    """The film coupled to the tank bulk, on a pymrm finite-volume grid.

    Parameters
    ----------
    n : int
        Number of cells across the film.
    x_f : array, optional
        Face coordinates on [0, 1]; uniform when omitted.
    permc_spec : str, optional
        Column ordering of the sparse LU in ``newton`` (the round-off companion
        uses a second one).

    Notes
    -----
    The operators are assembled once. The left face carries the Dirichlet value
    u(0) = 1 as a constant vector; the right face carries u(1) = u_b through
    ``shapes_d``, with d = 1 as the coefficient of the external unknown, so the
    gradient operator never changes while u_b does.
    """

    def __init__(self, n, x_f=None, permc_spec=None):
        self.n = n
        self.x_f = np.linspace(0.0, 1.0, n + 1) if x_f is None else np.asarray(x_f, float)
        self.x_c = 0.5 * (self.x_f[1:] + self.x_f[:-1])
        self.h = np.diff(self.x_f)
        # left, xi = 0: u = 1 (b u = d with b = 1, d = 1: a fixed interface value)
        # right, xi = 1: u = u_b (b u = d with d = 1 the COEFFICIENT of the
        # external unknown u_b, pitfall P5)
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
        self.G, gl, gbr = construct_grad((n, 1), self.x_f, self.x_c, bc, axis=0,
                                         shapes_d=(None, (1, 1)))
        self.gl = gl.toarray().ravel()            # face gradients from u(0) = 1
        self.gbr = gbr.toarray().ravel()          # face gradients per unit u_b
        self.Dv = construct_div((n, 1), self.x_f, nu=0)   # nu = 0: Cartesian film
        self.L = self.Dv @ self.G
        self.sL = self.Dv @ self.gl
        self.bR = self.Dv @ self.gbr
        self.G_last = self.G[[-1], :]
        self.lin = {} if permc_spec is None else {"permc_spec": permc_spec}

    def face_gradient(self, u, ub):
        """du/dxi on every face, the bc vectors included (pitfall P8)."""
        return self.G @ u + self.gl + self.gbr * ub

    def _guess(self, Ha2, Al, m, film_excluded, bulk):
        if m == 1:
            Ha = np.sqrt(Ha2)
            vol = Al - 1.0 if film_excluded else Al
            ub = 0.0 if bulk == "dirichlet0" else 1.0 / (np.cosh(Ha) + Ha * vol * np.sinh(Ha))
            u = np.cosh(Ha * self.x_c) + (ub - np.cosh(Ha)) / np.sinh(Ha) * np.sinh(Ha * self.x_c)
            return u, ub
        # m != 1: the m = 1 profile at the same Ha^2 and Al
        return self._guess(Ha2, Al, 1, film_excluded, bulk)

    # codehealth: allow params -- physical inputs and the defects the page injects
    def solve(self, Ha2, Al, m, *, film_excluded=True, bulk="coupled", read="face",
              share="flux", tol=1e-10, maxfev=60, guess=None):
        """Solve for one (Ha^2, Al, m).

        Parameters
        ----------
        Ha2, Al, m : float
            Hatta number squared, Hinterland ratio, reaction order.
        film_excluded : bool
            Bulk volume per area Al - 1 (True) or Al (a defect).
        bulk : {"coupled", "dirichlet0"}
            u_b an unknown of the bordered system, or u_b = 0 (film only).
        read : {"face", "cell"}
            phi from the operator's left face gradient, or (1 - u_1)/(h_1/2) from
            the first cell (a defect).
        share : {"flux", "ub"}
            film share 1 - g/phi, or 1 - u_b (a defect).

        Returns
        -------
        dict
            ``phi``, ``ub``, ``g``, ``film_share``, ``E``, ``balance`` (phi - g -
            Ha^2 sum h u^m, relative to phi), ``nit``, ``u``.
        """
        n = self.n
        vol = Al - 1.0 if film_excluded else Al
        coupled = bulk == "coupled"
        u0, ub0 = self._guess(Ha2, Al, m, film_excluded, bulk) if guess is None else guess

        def fun(x):
            u = x[:n]
            ub = x[n] if coupled else 0.0
            r, dr = _power(u, m)
            res_f = self.L @ u + self.sL + self.bR * ub - Ha2 * r
            J_ff = self.L - sp.diags_array(Ha2 * dr, format="csc")
            if not coupled:
                return res_f, J_ff
            rb, drb = _power(np.array([ub]), m)
            g = -(self.G_last @ u)[0] - self.gl[-1] - self.gbr[-1] * ub
            res_b = g - Ha2 * vol * rb[0]
            J = sp.block_array(
                [[J_ff, sp.csc_array(self.bR.reshape(n, 1))],
                 [-self.G_last, sp.csc_array([[-self.gbr[-1] - Ha2 * vol * drb[0]]])]],
                format="csc")
            return np.append(res_f, res_b), J

        x0 = np.append(u0, ub0) if coupled else np.asarray(u0, float)
        sol = newton(fun, x0, tol=tol, maxfev=maxfev, lin_solver_kwargs=self.lin or None)
        if not sol.success:
            raise RuntimeError(f"newton failed at Ha2={Ha2:g}, Al={Al:g}, m={m}: {sol.message}")
        u = sol.x[:n]
        ub = float(sol.x[n]) if coupled else 0.0
        gf = self.face_gradient(u, ub)
        phi = -gf[0] if read == "face" else (1.0 - u[0]) / (0.5 * self.h[0])
        g = -gf[-1]
        r, _ = _power(u, m)
        bal = (phi - g - Ha2 * np.sum(self.h * r)) / phi   # M7, with phi as read
        fs = 1.0 - g / phi if share == "flux" else 1.0 - ub
        return dict(phi=phi, ub=ub, g=g, film_share=fs, E=phi / (1.0 - ub),
                    balance=bal, nit=sol.nit, u=u)


def mapped_grid(n, a=4.0):
    """Faces x = (e^(a xi) - 1)/(e^a - 1) on a uniform xi grid of n cells: fine at
    the interface (xi = 0), and nested under doubling, so a Richardson ladder on
    it stays valid."""
    return np.expm1(a * np.linspace(0.0, 1.0, n + 1)) / np.expm1(a)


class PymrmRoute:
    """``solve(Ha2, Al, m) -> dict`` on a fixed :class:`FilmTank`, warm-started.

    A solve close to the last one (same Al and m, Ha^2 within a factor e^0.5)
    starts from the last profile; any other starts from the m = 1 closed-form
    profile. ``flat_guess`` starts every solve from that uniform profile instead
    (a different Newton path, which the round-off companion uses). Keyword arguments go to :meth:`FilmTank.solve`.
    """

    # codehealth: allow params -- physical inputs and the defects the page injects
    def __init__(self, n, x_f=None, permc_spec=None, mapped=None, flat_guess=None, **opts):
        if mapped is not None:
            x_f = mapped_grid(n, mapped)
        self.tank = FilmTank(n, x_f=x_f, permc_spec=permc_spec)
        self.opts = opts
        self.flat = flat_guess
        self.last = None
        self.solves = 0

    def __call__(self, Ha2, Al, m):
        guess = None
        if self.flat is not None:
            guess = (np.full(self.tank.n, self.flat), 0.5 * self.flat)
        elif self.last is not None:
            h2, al, mm, u, ub = self.last
            if al == Al and mm == m and abs(np.log(Ha2 / h2)) < 0.5:
                guess = (u, ub)
        s = self.tank.solve(Ha2, Al, m, guess=guess, **self.opts)
        self.last = (Ha2, Al, m, s["u"], s["ub"])
        self.solves += 1
        return s


def _one_root(f, lo, hi, n_scan, what):
    """The single root of f on [lo, hi] (log scan), or raise."""
    res = find_roots(lambda s: f(np.exp(s)), np.log(lo), np.log(hi), n_scan=n_scan,
                     xtol=1e-14)
    if len(res["roots"]) != 1 or res["failed"]:
        raise RuntimeError(f"{what}: {len(res['roots'])} roots on [{lo:g}, {hi:g}], "
                           f"{len(res['failed'])} failed points ({res['message']})")
    return float(np.exp(res["roots"][0]))


# codehealth: allow params -- physical inputs and the defects the page injects
def edge(solve, Al, m, level, lo=1e-8, hi=1e3, n_scan=9):
    """Ha^2 at which phi = level (phi rises monotonically with Ha^2)."""
    return _one_root(lambda h2: solve(h2, Al, m)["phi"] - level, lo, hi, n_scan,
                     f"edge phi = {level:g} at Al = {Al:g}, m = {m}")


# codehealth: allow params -- physical inputs and the defects the page injects
def band(solve, Al, m, tol, **kw):
    """Lower edge (phi = 1 - tol), upper edge (phi = 1 + tol) and W = hi/lo."""
    lo = edge(solve, Al, m, 1.0 - tol, **kw)
    hi = edge(solve, Al, m, 1.0 + tol, **kw)
    return lo, hi, hi / lo


# codehealth: allow params -- physical inputs and the defects the page injects
def point_where(solve, Al, m, field, target, lo=1e-8, hi=1e3, n_scan=12):
    """Ha^2 at which ``field`` (``ub`` or ``film_share``) equals ``target``,
    and the solution there."""
    h2 = _one_root(lambda x: solve(x, Al, m)[field] - target, lo, hi, n_scan,
                   f"{field} = {target:g} at Al = {Al:g}, m = {m}")
    return h2, solve(h2, Al, m)


def al_root(Wfun, target, lo, hi, n_scan=6):
    """Al at which the band width W(Al) equals ``target`` (log scan)."""
    return _one_root(lambda al: Wfun(al) - target, lo, hi, n_scan,
                     f"W = {target:g}")


# codehealth: allow params -- physical inputs and the defects the page injects
def box_share(c, pv_box, vs_box, alpha, beta, log=True):
    """Share of the (P/V, v_s) box where alpha ln(P/V) + beta ln v_s <= c.

    ``log=True`` measures the share as a log-area (the box is read on log axes,
    as eq. 9 is a power law); ``log=False`` as a linear area (a defect). The
    region is a half-plane in log coordinates; the share is integrated exactly
    over v_s piece by piece.
    """
    from scipy.integrate import quad
    xl, xh = np.log(pv_box[0]), np.log(pv_box[1])
    yl, yh = np.log(vs_box[0]), np.log(vs_box[1])

    def frac_x(y):                      # y = ln v_s
        xs = (c - beta * y) / alpha     # ln(P/V) on the dividing line
        if log:
            return np.clip((xs - xl) / (xh - xl), 0.0, 1.0)
        return np.clip((np.exp(xs) - pv_box[0]) / (pv_box[1] - pv_box[0]), 0.0, 1.0)

    # the line leaves the box's P/V range at two values of ln v_s: split there
    ys = sorted(y for y in ((c - alpha * xl) / beta, (c - alpha * xh) / beta) if yl < y < yh)
    knots = [yl, *ys, yh]
    if log:
        tot = sum(quad(frac_x, a, b, epsabs=0, epsrel=1e-13)[0] for a, b in zip(knots, knots[1:]))
        return tot / (yh - yl)
    tot = sum(quad(lambda v: frac_x(np.log(v)), np.exp(a), np.exp(b), epsabs=0,
                   epsrel=1e-13)[0] for a, b in zip(knots, knots[1:]))
    return tot / (vs_box[1] - vs_box[0])


def richardson(values, ns, orders=(2.0, 3.0)):
    """Richardson limit of a quantity computed on grids ``ns`` (ratio 2).

    The finite-volume error expands as c2 h^2 + c3 h^3 + ... (the one-sided
    Dirichlet face gradients add the odd term), so three levels remove both
    leading terms. ``values`` holds the quantity on each of the last
    ``len(orders) + 1`` levels of ``ns``.

    Returns ``(limit, error)``: the limit, and the change made by the last
    elimination (the estimate of what remains).
    """
    v = [np.asarray(x, float) for x in values][-(len(orders) + 1):]
    r = ns[-1] / ns[-2]
    prev = None
    for p in orders:
        f = r ** p - 1.0
        prev = v[-1]
        v = [b + (b - a) / f for a, b in zip(v[:-1], v[1:])]
    lim = v[-1]
    return lim, np.abs(lim - prev)


# ----------------------------------------------------------------------------
# The studies: every published constant in one dictionary ``P``; the solver
# defects the break table injects in ``opt`` (keywords of FilmTank.solve).
# ----------------------------------------------------------------------------

NS_LIMIT = (200, 400, 800)   # the Richardson levels of every pymrm value


class SulfiteTank:
    """The reactor over van 't Riet's box, one dictionary of constants.

    Parameters
    ----------
    P : dict
        ``K``, ``alpha``, ``beta`` (eq. 9), ``vs_unit`` (1 for m/s), ``kL``
        (m/s), ``D`` (m2/s, assumed), ``pv_T``, ``vs_T`` (the test condition),
        ``pv_box``, ``vs_box`` (pairs), ``F`` (the printed k_R shift).
    ns : tuple of int
        Richardson levels for the pymrm route.
    m1_route : {"closed", "pymrm"}
        Route for m = 1 (the closed form is exact; pymrm is its check).
    **opt
        Solver options and defects passed to the route (``film_excluded``,
        ``bulk``, ``read``, ``share``; ``read`` acts on pymrm only).
    """

    def __init__(self, P, ns=NS_LIMIT, m1_route="closed", **opt):
        self.P, self.ns, self.m1_route, self.opt = dict(P), tuple(ns), m1_route, opt

    # -- k_L A and Al --------------------------------------------------------
    def kla(self, pv, vs):
        P = self.P
        return kla_eq9(pv, vs, P["K"], P["alpha"], P["beta"], P.get("vs_unit", 1.0))

    def Al(self, pv, vs, kL=None):
        return hinterland(self.P["kL"] if kL is None else kL, self.kla(pv, vs), self.P["D"])

    def conditions(self, kL=None):
        """k_L A and Al at T and at the box corners (Al falls as k_L A rises)."""
        P = self.P
        pv, vs = P["pv_box"], P["vs_box"]
        return dict(kla_T=self.kla(P["pv_T"], P["vs_T"]), Al_T=self.Al(P["pv_T"], P["vs_T"], kL),
                    kla_min=self.kla(pv[0], vs[0]), kla_max=self.kla(pv[1], vs[1]),
                    Al_min=self.Al(pv[1], vs[1], kL), Al_max=self.Al(pv[0], vs[0], kL))

    # -- routes ----------------------------------------------------------------
    def evaluate(self, fun, m, **extra):
        """``fun(solve)`` on the route for order m: exact for m = 1 by the closed
        form, else the Richardson limit of the pymrm levels. Returns
        ``(value, grid_error)``, elementwise for a tuple."""
        o = {**self.opt, **extra}
        if m == 1 and self.m1_route == "closed":
            co = {k: v for k, v in o.items() if k in ("film_excluded", "bulk", "share")}
            val = np.asarray(fun(lambda h2, al, mm: closed_m1(h2, al, mm, **co)), float)
            return val, np.zeros_like(val)
        vals = [fun(PymrmRoute(n, **o)) for n in self.ns]
        return richardson(vals, self.ns)

    # -- results -------------------------------------------------------------
    def point(self, m, Ha2=1.0, Al=None):
        """phi, u_b, film share and E at (Ha^2, Al); Al defaults to Al_T."""
        al = self.conditions()["Al_T"] if Al is None else Al
        keys = ("phi", "ub", "film_share", "E", "g")
        return self.evaluate(lambda s: [s(Ha2, al, m)[k] for k in keys], m)

    def band(self, m, tol, Al=None):
        """Ha2_lo, Ha2_hi and W at Al (default Al_T)."""
        al = self.conditions()["Al_T"] if Al is None else Al
        return self.evaluate(lambda s: band(s, al, m, tol), m)

    def al_F(self, m, tol, F=None, lo=1.5, hi=3000.0):
        """Al at which W = F (F defaults to the printed k_R shift)."""
        F = self.P["F"] if F is None else F
        if m == 1 and self.m1_route == "closed":
            return self.evaluate(lambda s: [al_root(lambda al: band(s, al, m, tol)[2], F, lo, hi)], m)
        # level by level, each finer level scanning a bracket around the coarser root
        vals = []
        for n in self.ns:
            s = PymrmRoute(n, **self.opt)
            a, b, k = (lo, hi, 6) if not vals else (vals[-1] / 1.25, vals[-1] * 1.25, 3)
            vals.append(al_root(lambda al: band(s, al, m, tol)[2], F, a, b, n_scan=k))
        return richardson([[v] for v in vals], self.ns)

    def box(self, m, tol, kL=None, al_F=None):
        """W at the two Al corners, the log-area share with W >= F, the margin
        Al_min/Al_F, and the diffusivities at which the whole box (D_all) or none
        of it (D_none) has W >= F (Al is proportional to 1/D)."""
        c = self.conditions(kL)
        alF = self.al_F(m, tol)[0][0] if al_F is None else al_F
        w_min = self.band(m, tol, Al=c["Al_min"])[0][2]
        w_max = self.band(m, tol, Al=c["Al_max"])[0][2]
        P = self.P
        kl = P["kL"] if kL is None else kL
        # W >= F  <=>  Al >= Al_F  <=>  alpha ln(P/V) + beta ln v_s <= ln(kL^2/(D Al_F K vs_unit^beta))
        cc = np.log(kl ** 2 / (P["D"] * alF * P["K"] * P.get("vs_unit", 1.0) ** P["beta"]))
        share = box_share(cc, P["pv_box"], P["vs_box"], P["alpha"], P["beta"],
                          log=P.get("log_area", True))
        return dict(W_min=w_min, W_max=w_max, share=share, margin_min=c["Al_min"] / alF,
                    margin_max=c["Al_max"] / alF, D_all=P["D"] * c["Al_min"] / alF,
                    D_none=P["D"] * c["Al_max"] / alF)

    def where(self, m, field, target, Al=None):
        """Ha^2 at which ``field`` equals ``target``, with phi, u_b and the film
        share there."""
        al = self.conditions()["Al_T"] if Al is None else Al

        def fun(s):
            h2, r = point_where(s, al, m, field, target)
            return [h2, r["phi"], r["ub"], r["film_share"]]
        return self.evaluate(fun, m)

    def nil_flip(self, m, ub_nil=0.01, share=0.5, lo=1.5, hi=1e6):
        """Al at which the film share at the nil point (u_b = ``ub_nil``) equals
        ``share``: the film share there falls monotonically with Al."""
        def fun(s):
            f = lambda al: point_where(s, al, m, "ub", ub_nil)[1]["film_share"] - share
            return [_one_root(f, lo, hi, 6, f"film share {share:g} at u_b = {ub_nil:g}")]
        return self.evaluate(fun, m)

    def film_limit(self, m, tol, Al=None):
        """The upper edge of the film alone (Al -> infinity: u_b = 0), or at a
        finite Al when one is passed (a defect: the limit read at a point)."""
        if Al is None:
            return self.evaluate(lambda s: [edge(s, 1.0, m, 1.0 + tol)], m, bulk="dirichlet0")
        return self.evaluate(lambda s: [edge(s, Al, m, 1.0 + tol)], m)


def sulfite_numbers(c):
    """Arithmetic on the printed sulfite numbers (R5).

    ``c``: ``M`` (g/mol, Na2SO3), ``g_lo``, ``g_hi`` (g/L), ``ion_factor`` (ionic
    strength per mol of salt, 3 for Na2SO3), ``C_star`` (mol/m3) and ``DB_over_D``
    (both assumed), ``nu`` (mol sulfite per mol O2).
    """
    lo, hi = c["g_lo"] / c["M"], c["g_hi"] / c["M"]
    I_lo, I_hi = c["ion_factor"] * lo, c["ion_factor"] * hi
    E_inf = 1.0 + c["DB_over_D"] * lo * 1e3 / (c["nu"] * c["C_star"])
    return dict(molar_lo=lo, molar_hi=hi, ionic_lo=I_lo, ionic_hi=I_hi, E_inf_10gL=E_inf)
