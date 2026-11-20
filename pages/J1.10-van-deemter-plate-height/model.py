"""The van Deemter (1956) column, its moments, and the printed formulae it is checked with.

The transient model is eqs. (29)/(30) of van Deemter, Zuiderweg and Klinkenberg
(1956): a linear two-phase chromatographic column, mobile phase ``c_I`` with
axial dispersion ``D`` and interstitial velocity ``u``, immobile phase ``c_II``
exchanging with it at the rate ``alpha (c_I - K c_II)`` per unit volume,

    F_I  dc_I/dt  = F_I D d2c_I/dz2 - F_I u dc_I/dz + alpha (K c_II - c_I),
    F_II dc_II/dt = alpha (c_I - K c_II),

with a Danckwerts inlet ``u c_I - D dc_I/dz = u c_in(t)`` at ``z = 0`` and a
zero gradient at ``z = L``. The plate height is read from the first two time
moments at two interior stations, so every entrance and exit effect cancels.

Contents
--------
Column
    Operators from ``construct_convflux_upwind``, ``construct_grad`` and
    ``construct_div`` assembled once; a theta-scheme march with the van Leer
    flux as a deferred correction (one linear solve per sweep), and the same
    theta scheme with the TVD flux of ``ConvFluxTVD`` solved by ``newton``.
cumulants, plate_height, kappas_per_z
    Moments of a recorded curve and the two-station readings built on them.
one_root
    ``find_roots`` on a physical range, raising unless exactly one root.
delta_of, width_red, cmax, n_from_height, n_from_width
    Eqs. (20), (21) and (22) of the paper and the two routes from an elution
    curve to a plate number.
vd, refit, fit_eq49, dii_from_slope
    The three-term plate height and the least-squares fits made with it.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import bmat, eye_array
from scipy.sparse.linalg import splu
from scipy.special import erf, erfinv

import pymrm
from pymrm import (ConvFluxTVD, construct_convflux_upwind, construct_div, construct_grad,
                   interp_cntr_to_stagg_tvd, vanleer)
from pymrm.checks import find_roots

__all__ = ["Column", "cumulants", "plate_height", "kappas_per_z", "one_root", "delta_of",
           "width_red", "cmax", "n_from_height", "n_from_width", "vd", "refit", "fit_eq49",
           "dii_from_slope"]

SQ2 = np.sqrt(2.0)


class Column:
    """Eqs. (29)/(30): a linear two-phase chromatographic column.

    Parameters
    ----------
    u : float
        Interstitial velocity (mm/s).
    D : float
        Axial dispersion coefficient (mm2/s).
    alpha : float
        Mass-transfer coefficient per unit volume of packing (1/s).
    K : float
        Distribution factor ``(c_I / c_II)`` at equilibrium.
    F_I, F_II : float
        Fractional volumes of the mobile and the immobile phase.
    L : float
        Column length (mm).
    n_z : int
        Number of cells.
    limiter : callable or None, optional
        TVD limiter of the convective face value; ``None`` is first-order upwind.
    nu : int, optional
        Geometry index of ``construct_div`` (0 Cartesian, the column; 1 and 2
        exist only so that a break row can apply the wrong geometry).
    bc_left : dict, optional
        Inlet condition for ``c_in = 1``; default Danckwerts.

    Notes
    -----
    Fields have shape ``(n_z, 1)`` per phase and the two phases are stacked into
    one ``2 n_z`` system, ``dc/dt = M c + b c_in(t) + (TVD correction)``. Only
    ``c_I`` carries a spatial operator; ``c_II`` exchanges with it pointwise.
    """

    def __init__(self, u, D, alpha, K, F_I, F_II, L, n_z, limiter=vanleer,
                 nu=0, bc_left=None):
        self.__dict__.update(u=u, D=D, alpha=alpha, K=K, F_I=F_I, F_II=F_II,
                             L=L, n_z=n_z, limiter=limiter)
        self.x_f = np.linspace(0.0, L, n_z + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        shape = (n_z, 1)                       # one field on a grid: never (n_z,)
        # outward normal n, both dictionaries read a dc/dn + b c = d:
        #   z = 0 (n = -z): Danckwerts inflow u c - D dc/dz = u c_in  -> a = D, b = u, d = u
        #   z = L (n = +z): zero gradient dc/dz = 0                   -> a = 1, b = 0, d = 0
        # d is written for c_in = 1; every boundary term is linear in d and scaled by c_in(t).
        self.bc = (bc_left or {"a": D, "b": u, "d": u},
                   {"a": 1.0, "b": 0.0, "d": 0.0})
        conv, conv_bc = construct_convflux_upwind(shape, self.x_f, self.x_c, self.bc, v=u)
        grad, grad_bc = construct_grad(shape, self.x_f, self.x_c, self.bc)
        div = construct_div(shape, self.x_f, nu=nu)        # nu = 0: Cartesian column
        self.div, self.grad, self.grad_bc = div, grad, grad_bc
        T = div @ (conv - D * grad)                        # upwind transport of c_I
        self.bvec = -np.asarray((div @ (conv_bc - D * grad_bc)).todense()).ravel()
        eye = eye_array(n_z, format="csc")
        # d/dt [c_I ; c_II] = M [c_I ; c_II] + b c_in(t)
        self.exchange = bmat([[-(alpha / F_I) * eye, (alpha * K / F_I) * eye],
                              [(alpha / F_II) * eye, -(alpha * K / F_II) * eye]], format="csc")
        self.M = bmat([[-T - (alpha / F_I) * eye, (alpha * K / F_I) * eye],
                       [(alpha / F_II) * eye, -(alpha * K / F_II) * eye]], format="csc")
        self.v_face = np.full((n_z + 1, 1), u)

    def _corr(self, cI, cin=1.0):
        """Deferred van Leer correction of the convective term (zero when ``limiter`` is None).

        The TVD face value minus the upwind one, with the inlet data of the
        current ``c_in`` (the first interior face reads the boundary value).
        """
        if self.limiter is None:
            return np.zeros(self.n_z)
        bc = (dict(self.bc[0], d=self.bc[0]["d"] * cin), self.bc[1])
        _, dg = interp_cntr_to_stagg_tvd(cI.reshape(-1, 1), self.x_f, self.x_c, bc,
                                         self.v_face, tvd_limiter=self.limiter)
        return -np.asarray(self.div @ (self.v_face * dg.reshape(-1, 1))).ravel()

    def _stations(self, stations):
        """Cell index for each station ``z``: the cell whose right face is the first face at or beyond ``z``."""
        idx = [int(np.searchsorted(self.x_f, z) - 1) for z in stations]
        if min(idx) < 0 or max(idx) >= self.n_z:
            raise ValueError(f"station outside the column: {stations}")
        return idx

    @staticmethod
    def inlet(t, t_sig):
        """Gaussian inlet pulse of width ``t_sig`` centred at ``5 t_sig``."""
        return float(np.exp(-0.5 * ((t - 5 * t_sig) / t_sig) ** 2))

    def run(self, t_end, dt, t_sig, stations, sweeps=1, theta=0.5, inlet_data=True):
        """March the theta scheme with the TVD flux as a deferred correction.

        Parameters
        ----------
        t_end, dt : float
            End time and step (s); ``t_end / dt`` is rounded to whole steps.
        t_sig : float
            Width of the Gaussian inlet pulse (s).
        stations : sequence of float
            Positions (mm) at which ``c_I`` is recorded.
        sweeps : int, optional
            Linear solves per step. The correction in sweep ``k`` is
            ``(1 - theta) g(c^n) + theta g(c^(n+1, k-1))`` with
            ``c^(n+1, 0) = c^n``: one sweep lags the correction a whole step
            (first order in time), two sweeps are formally second order with
            ``theta = 0.5``.
        theta : float, optional
            0.5 Crank-Nicolson, 1 backward Euler.
        inlet_data : bool, optional
            Give the correction the inlet data of the current ``c_in``; False
            keeps the data of ``c_in = 1`` at every step.

        Returns
        -------
        ts : ndarray
            Times of the records.
        rec : ndarray
            ``c_I`` at the stations, shape ``(len(ts), len(stations))``.
        zs : ndarray
            Cell centres of the stations.
        """
        n = self.n_z
        lu = splu(eye_array(2 * n, format="csc") / dt - theta * self.M)
        c = np.zeros(2 * n)
        idx = self._stations(stations)
        nt = int(round(t_end / dt))
        rec = np.zeros((nt + 1, len(idx)))
        ts = np.arange(nt + 1) * dt
        zeros = np.zeros(n)
        for k in range(nt):
            cin0, cin1 = self.inlet(k * dt, t_sig), self.inlet((k + 1) * dt, t_sig)
            b0 = np.concatenate([self.bvec * cin0, zeros])
            b1 = np.concatenate([self.bvec * cin1, zeros])
            rhs0 = c / dt + (1 - theta) * (self.M @ c) + ((1 - theta) * b0 + theta * b1)
            g0 = self._corr(c[:n], cin0 if inlet_data else 1.0)
            cn = c
            for _ in range(sweeps):
                g1 = self._corr(cn[:n], cin1 if inlet_data else 1.0)
                cn = lu.solve(rhs0 + np.concatenate([(1 - theta) * g0 + theta * g1, zeros]))
            c = cn
            rec[k + 1] = c[:n][idx]
        return ts, rec, self.x_c[idx]

    def run_newton(self, t_end, dt, t_sig, stations, theta=0.5, tol=1e-12):
        """The same theta scheme with the TVD flux of ``ConvFluxTVD``, solved by ``newton``.

        Each step solves ``(c - c^n)/dt = theta f(c) + (1 - theta) f(c^n)`` with
        the exact Jacobian of the limited flux, ``damping="monotone"`` and
        ``stagnation=True``; a step that does not converge raises.
        Returns ``ts, rec, zs`` as :meth:`run`.
        """
        n = self.n_z
        if self.limiter is None:
            raise ValueError("run_newton needs a limiter; the upwind scheme is linear (use run)")
        tvd = ConvFluxTVD((n, 1), self.x_f, self.x_c, self.bc, tvd_limiter=self.limiter,
                          shapes_d=((1, 1), None))
        dgrad = self.D * self.grad

        def f(c, cin):
            flux, jac, jac_l, _ = tvd(c[:n], self.u, d=(np.array([[cin]]), None))
            diff = dgrad @ c[:n] + self.D * grad_bc * cin
            fI = -(self.div @ (flux.reshape(-1) - diff))
            jI = -(self.div @ (jac - dgrad))
            return np.concatenate([fI, np.zeros(n)]) + self.exchange @ c, jI

        eye = eye_array(2 * n, format="csc")
        zero = eye_array(n, format="csc") * 0.0
        grad_bc = np.asarray(self.grad_bc.todense()).ravel()
        c = np.zeros(2 * n)
        idx = self._stations(stations)
        nt = int(round(t_end / dt))
        rec = np.zeros((nt + 1, len(idx)))
        ts = np.arange(nt + 1) * dt
        for k in range(nt):
            cin0, cin1 = self.inlet(k * dt, t_sig), self.inlet((k + 1) * dt, t_sig)
            f0, _ = f(c, cin0)
            c_old = c

            def residual(x):
                f1, jI = f(x, cin1)
                jac = eye / dt - theta * (bmat([[jI, None], [None, zero]], format="csc") + self.exchange)
                return (x - c_old) / dt - theta * f1 - (1 - theta) * f0, jac

            sol = pymrm.newton(residual, c_old, tol=tol, damping="monotone", stagnation=True)
            if not sol.success:
                raise RuntimeError(f"newton failed at step {k + 1}: {sol.message}")
            c = sol.x
            rec[k + 1] = c[:n][idx]
        return ts, rec, self.x_c[idx]


def cumulants(ts, y):
    """Zeroth moment, mean, and the second to fourth cumulants of a recorded curve.

    Trapezoidal quadrature over the record; the fourth cumulant is
    ``m4 - 3 k2^2`` with ``m4`` the fourth central moment.
    """
    m0 = np.trapezoid(y, ts)
    m1 = np.trapezoid(ts * y, ts) / m0
    d = ts - m1
    k2 = np.trapezoid(d ** 2 * y, ts) / m0
    k3 = np.trapezoid(d ** 3 * y, ts) / m0
    m4 = np.trapezoid(d ** 4 * y, ts) / m0
    return m0, m1, k2, k3, m4 - 3 * k2 ** 2


def plate_height(Kc, zs):
    """Plate height from two stations: ``dz (k2(z2) - k2(z1)) / (m1(z2) - m1(z1))^2``.

    Differences of the cumulants of a linear system cancel every entrance and
    exit effect and the shape of the injected pulse.
    """
    dz = zs[1] - zs[0]
    return dz * (Kc[1][2] - Kc[0][2]) / (Kc[1][1] - Kc[0][1]) ** 2


def kappas_per_z(Kc, zs):
    """Third and fourth cumulants per unit length from two stations."""
    dz = zs[1] - zs[0]
    return ((Kc[1][3] - Kc[0][3]) / dz, (Kc[1][4] - Kc[0][4]) / dz)


def one_root(f, lo, hi, what, log=False, n_scan=64, xtol=1e-12):
    """The single root of ``f`` on the physical range ``[lo, hi]``.

    ``pymrm.checks.find_roots`` scans the range and refines each sign change;
    this raises ``RuntimeError`` unless exactly one root is found.
    """
    res = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(res["roots"]) != 1:
        raise RuntimeError(f"{what}: {len(res['roots'])} roots on [{lo:g}, {hi:g}] "
                           f"({res['message']})")
    return float(res["roots"][0])


def delta_of(a):
    """Eq. (22) for ``delta``: ``delta exp(-delta^2/2) = (a + delta) exp(-(a + delta)^2/2)``.

    Solved in log form, ``h(d) = log d - d^2/2 - log(a + d) + (a + d)^2/2``,
    which increases with ``d`` (``h' = a (1 + 1/(d (a + d))) > 0``), so its root
    in ``(0, 1)`` is unique.
    """
    def h(d):
        return np.log(d) - 0.5 * d * d - np.log(a + d) + 0.5 * (a + d) ** 2
    return one_root(h, 1e-300, 1.0 - 1e-15, f"eq (22) at a = {a:g}", xtol=1e-16)


def width_red(a, denominator="a"):
    """Eq. (21), the reduced width ``Delta s`` at the inflection points.

    ``denominator="a"`` is the form derived from eq. (19); ``"delta"`` is the
    third term exactly as printed.
    """
    d = delta_of(a)
    den = a if denominator == "a" else d
    return (a + 2 * d + np.sqrt(2 * np.pi) * (a + d) / den * np.exp(0.5 * d * d)
            * (erf((a + d) / SQ2) - erf(d / SQ2)))


def cmax(a):
    """Eq. (20): the band maximum ``erf(a / (2 sqrt 2))`` for a slug of reduced width ``a``."""
    return erf(a / (2 * SQ2))


def n_from_height(A, hmax, nv):
    """Plate number from the height of the elution curve (eq. 20); nan when saturated."""
    if hmax >= 1.0:
        return np.nan
    a = 2 * SQ2 * erfinv(hmax)
    return (nv / (A / a)) ** 2


def n_from_width(A, W, nv, denominator="a"):
    """Plate number from the base width ``W`` of the elution curve (eq. 21).

    Solves ``width_red(A / w) w = W`` for ``w = v sqrt(n)`` on ``[A/25, 1e4]``
    ml; raises unless exactly one root is found.
    """
    def f(w):
        return width_red(A / w, denominator) * w - W
    return (nv / one_root(f, A / 25.0, 1e4, f"width route, A = {A:g}, W = {W:g}",
                          log=True, xtol=1e-13)) ** 2


def vd(u, B, A, C):
    """The three-term plate height ``B/u + A + C u``."""
    return B / u + A + C * u


def refit(u, h, terms=("B", "A", "C")):
    """Least squares on ``H = B/u + A + C u`` (or a subset of the terms).

    Returns the coefficients, their standard errors and the rmse, each by term name.
    """
    cols = {"B": 1 / u, "A": np.ones_like(u), "C": u}
    X = np.vstack([cols[t] for t in terms]).T
    beta, *_ = np.linalg.lstsq(X, h, rcond=None)
    r = X @ beta - h
    dof = len(u) - X.shape[1]
    s2 = (r ** 2).sum() / dof
    se = np.sqrt(np.diag(np.linalg.inv(X.T @ X)) * s2)
    return dict(zip(terms, beta)), dict(zip(terms, se)), float(np.sqrt((r ** 2).mean()))


def fit_eq49(X, Y, three_term=False):
    """Least squares of eq. (49), ``H/d_p = 2 lambda + S X`` (with ``B'/X`` when ``three_term``).

    Returns the coefficients (``[1/X,] 1, X`` order) and the rmse.
    """
    cols = [1 / X, np.ones_like(X), X] if three_term else [np.ones_like(X), X]
    A = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(A, Y, rcond=None)
    r = A @ b - Y
    return b, float(np.sqrt((r ** 2).mean()))


def dii_from_slope(S, F_I, F_II, K, D_I=None, visc=75.0):
    """Invert eq. (49)'s ``C`` for ``D_II`` from the slope ``S`` of ``H/d_p`` against ``u_sup d_p``.

    ``D_I=None`` drops the small viscous first term; ``visc`` is the 75 of
    that term. Units follow ``S`` and ``D_I`` (mm2/s on this page).
    """
    a_c = F_I ** 2 / (visc * (1 - F_I) ** 2)
    b_c = F_I * K / (2 * np.pi ** 2 * (1 - F_I))
    xK = K * F_I / F_II
    if D_I is None:
        return b_c / (S * F_I * (1 + xK) ** 2)
    return b_c * D_I / (S * F_I * D_I * (1 + xK) ** 2 - a_c)
