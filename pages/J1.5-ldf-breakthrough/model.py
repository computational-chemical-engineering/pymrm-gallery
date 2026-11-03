"""Diffusion into a sphere, the Crank series and the linear driving force (LDF).

Dimensionless throughout: radius 1, diffusivity 1, so time is ``tau = D t / r^2``,
the surface value is held at 1 and the particle starts empty. ``q`` is the local
loading over its surface value and ``qbar`` its volume mean.

    dq/dtau = (1/rho^nu) d/drho( rho^nu dq/drho ),
    dq/drho = 0 at rho = 0,  q = 1 at rho = 1,  q = 0 at tau = 0,

with ``nu = 2`` for a sphere (``1`` a long cylinder, ``0`` a slab). The problem is
linear with constant coefficients, so every implicit step is one linear solve with
a matrix that is assembled and factorised once per step length.

Contents
--------
crank_series
    The exact series for the sphere, truncated where a measured tail bound is
    below a requested relative tolerance.
crank_uptake, ldf_uptake, k_effective, ldf_constant
    Mean loading from the series, from the LDF, the coefficient the LDF would
    need at each instant, and the parabolic-profile constant of a geometry.
Particle
    The finite-volume particle: operators from ``construct_grad`` and
    ``construct_div`` assembled once, a backward-Euler march, the same problem
    exact in time by an eigendecomposition, and the surface flux read with
    ``compute_boundary_values``.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eig
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu

from pymrm import compute_boundary_values, construct_div, construct_grad
from pymrm.checks import find_roots

__all__ = ["crank_tail", "crank_series", "crank_uptake", "crank_profile", "ldf_uptake",
           "k_effective", "ldf_constant", "ldf_error_extremes", "ldf_crossing",
           "k_effective_crossing", "Particle"]


def crank_tail(tau, n_terms):
    """Bound on the tail of both Crank sums after ``n_terms`` terms, relative to each sum.

    Returns the larger of the two relative bounds at every ``tau``.

    The terms ``t_n = exp(-n^2 pi^2 tau) / n^2`` (deficit) and ``e_n`` without
    the ``1/n^2`` (rate) have ratios ``t_(m+1)/t_m`` that decrease with ``m``, so
    the tail after ``N`` terms is at most ``t_(N+1) / (1 - t_(N+2)/t_(N+1))``,
    whenever that ratio is below one.
    """
    tau = np.asarray(tau, dtype=float)
    out = []
    for power in (2, 0):
        m = np.arange(1, n_terms + 1)[:, None]
        total = np.sum(np.exp(-m**2 * np.pi**2 * tau) / m**power, axis=0)
        a = np.exp(-(n_terms + 1) ** 2 * np.pi**2 * tau) / (n_terms + 1) ** power
        b = np.exp(-(n_terms + 2) ** 2 * np.pi**2 * tau) / (n_terms + 2) ** power
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(a > 0, b / a, 0.0)
        out.append(np.where(r < 1.0, a / (1.0 - r) / total, np.inf))
    return np.maximum(out[0], out[1])


def crank_series(tau, rtol=1e-15, n_terms=None):
    """Deficit and rate of the Crank series, with the number of terms it needed.

    Parameters
    ----------
    tau : float or numpy.ndarray
        Dimensionless time ``D t / r^2``.
    rtol : float, optional
        Largest tolerated relative tail of either sum, from a bound proved for
        decreasing term ratios (see ``crank_tail``); the number of terms is the
        smallest that meets it at the smallest ``tau`` given.
    n_terms : int, optional
        Use this many terms instead (for a scan that calls the function many
        times); raises when the tail bound at some ``tau`` exceeds ``rtol``.

    Returns
    -------
    dict
        ``deficit`` ``1 - qbar`` and ``rate`` ``d qbar / d tau`` at every tau,
        ``n_terms`` and ``tail`` (the measured worst relative tail bound).
    """
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    if np.any(tau <= 0):
        raise ValueError("tau must be positive")
    if n_terms is None:
        lo, hi = 1, 1
        while crank_tail(tau, hi).max() > rtol:
            lo, hi = hi, 2 * hi
            if hi > 10**6:
                raise ValueError("the series needs more than 1e6 terms; tau is too small")
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if crank_tail(tau, mid).max() <= rtol else (mid, hi)
        n_terms = hi
    elif crank_tail(tau, n_terms).max() > rtol:
        raise ValueError(f"{n_terms} terms leave a tail bound above {rtol:g}")
    n = np.arange(1, n_terms + 1)[:, None]
    e = np.exp(-n**2 * np.pi**2 * tau)
    return {"deficit": 6.0 / np.pi**2 * np.sum(e / n**2, axis=0),   # 1 - qbar
            "rate": 6.0 * np.sum(e, axis=0),                         # d qbar / d tau
            "n_terms": n_terms, "tail": float(crank_tail(tau, n_terms).max())}


def crank_uptake(tau, rtol=1e-15, n_terms=None):
    """Exact mean loading ``qbar`` of a sphere, ``1 - deficit``."""
    return 1.0 - crank_series(tau, rtol, n_terms)["deficit"]


def ldf_uptake(tau, k=15.0):
    """Mean loading of the LDF ``d qbar / d tau = k (1 - qbar)`` from an empty particle."""
    return 1.0 - np.exp(-k * np.asarray(tau, dtype=float))


def k_effective(tau, rtol=1e-15, n_terms=None):
    """Coefficient the LDF would need at each instant, ``(d qbar/d tau) / (1 - qbar)``.

    Both the rate and the deficit are summed from the series, so nothing is
    differentiated numerically and nothing is formed as ``1 - (1 - eps)``.
    """
    s = crank_series(tau, rtol, n_terms)
    return s["rate"] / s["deficit"]


def crank_profile(rho, tau, n_terms):
    """Exact loading profile ``q(rho, tau)`` of a sphere.

    ``q = 1 + 2 sum_n (-1)^n sinc(n rho) exp(-n^2 pi^2 tau)``, with
    ``sinc(x) = sin(pi x)/(pi x)`` so that the centre needs no special case.
    Pass the ``n_terms`` that ``crank_series`` used at this ``tau``.
    """
    rho = np.asarray(rho, dtype=float)
    n = np.arange(1, n_terms + 1)[:, None]
    terms = (-1.0) ** n * np.sinc(n * rho[None, :]) * np.exp(-n**2 * np.pi**2 * tau)
    return 1.0 + 2.0 * terms.sum(axis=0)


def ldf_constant(nu=2):
    """LDF constant of the parabolic profile ``q = a + b rho^2`` in geometry ``nu``.

    With ``q(1) = 1`` the mean is ``a + b (nu + 1)/(nu + 3)``, so ``1 - qbar =
    2 b/(nu + 3)``, while the surface gradient ``2 b`` gives ``d qbar/d tau =
    (nu + 1) 2 b``. Their ratio is ``(nu + 1)(nu + 3)``: 15, 8 and 3.
    """
    b = 1.0
    deficit = 2.0 * b / (nu + 3)
    rate = (nu + 1) * 2.0 * b
    return rate / deficit


def _one_root(f, lo, hi, what):
    res = find_roots(f, lo, hi, n_scan=200, log=True, xtol=1e-14)
    if len(res["roots"]) != 1:
        raise RuntimeError(f"{what}: {res['message']}")
    return res["roots"][0]


def ldf_error_extremes(k=15.0, tau_lo=1e-4, tau_hi=3.0):
    """The two extrema of ``qbar_LDF - qbar_exact`` over ``tau``.

    The extrema are the zeros of ``d/d tau`` of the error, ``k exp(-k tau) -
    d qbar/d tau``, found with ``find_roots`` on ``[tau_lo, tau_hi]``.

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        The two times and the error there: first the under-prediction (negative),
        then the over-prediction (positive). Raises unless there are exactly two.
    """
    n = crank_series(tau_lo)["n_terms"]

    def slope(t):
        return k * np.exp(-k * t) - crank_series(t, n_terms=n)["rate"][0]

    res = find_roots(slope, tau_lo, tau_hi, n_scan=200, log=True, xtol=1e-14)
    if len(res["roots"]) != 2:
        raise RuntimeError(f"expected two extrema of the LDF error: {res['message']}")
    t = np.array(res["roots"])
    return t, ldf_uptake(t, k) - crank_uptake(t, n_terms=n)


def ldf_crossing(k=15.0, tau_lo=1e-4, tau_hi=3.0):
    """The one ``tau`` at which the LDF mean loading crosses the exact one."""
    n = crank_series(tau_lo)["n_terms"]
    return _one_root(lambda t: ldf_uptake(t, k) - crank_uptake(t, n_terms=n)[0],
                     tau_lo, tau_hi, "crossing of the LDF and the exact uptake")


def k_effective_crossing(target=15.0, tau_lo=1e-4, tau_hi=3.0):
    """The one ``tau`` at which ``k_effective`` equals ``target``."""
    n = crank_series(tau_lo)["n_terms"]
    return _one_root(lambda t: k_effective(t, n_terms=n)[0] - target, tau_lo, tau_hi,
                     "k_effective equal to the target")


class Particle:
    """Transient diffusion into a sphere (or cylinder, or slab) on ``n_r`` cells.

    The state has shape ``(n_r, 1)``: space first, one field. The operator
    ``(1/rho^nu) d/drho(rho^nu d/drho)`` and its boundary vector are assembled once
    in the constructor, from ``construct_grad`` and ``construct_div``. The
    problem is linear, so a backward-Euler step is one linear solve with the
    constant matrix ``I/dt + J``, factorised once per step length (a Newton
    iteration would converge in one step and add nothing).

    Parameters
    ----------
    n_r : int, optional
        Number of cells, uniform in ``rho`` from 0 to 1.
    nu : int, optional
        Geometry: 2 sphere, 1 cylinder, 0 slab.
    nu_div : int, optional
        Geometry index handed to ``construct_div``; defaults to ``nu``. A
        different value is a deliberate defect for the break table.
    weights : {"volume", "midpoint"}, optional
        Cell weights of the mean: the exact cell volumes ``rho_f^(nu+1)``
        differences, or the midpoint rule ``(nu + 1) rho_c^nu d rho`` (a
        deliberate defect: it does not sum to one).
    flux_from : {"face", "cell"}, optional
        Surface flux from ``compute_boundary_values`` (the face value the
        operators transport) or from the last cell centre (a deliberate defect).
    """

    def __init__(self, n_r=200, nu=2, nu_div=None, weights="volume", flux_from="face"):
        self.n_r, self.nu, self.flux_from = int(n_r), int(nu), flux_from
        self.shape = (self.n_r, 1)                       # (space, field), never (n_r,)
        self.r_f = np.linspace(0.0, 1.0, self.n_r + 1)
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        # Outward normal n, a dq/dn + b q = d:
        #   centre  (n = -rho): symmetry, dq/drho = 0  -> a = 1, b = 0, d = 0
        #   surface (n = +rho): fixed surface value, q = 1 -> a = 0, b = 1, d = 1
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
        grad, grad_bc = construct_grad(self.shape, self.r_f, self.r_c, self.bc)
        div = construct_div(self.shape, self.r_f, nu=self.nu if nu_div is None else nu_div)
        self.jac = (-(div @ grad)).tocsc()                # dq/dtau = -jac q + g
        self.g = (div @ grad_bc).toarray().ravel()
        if weights == "volume":
            w = np.diff(self.r_f ** (self.nu + 1))
        elif weights == "midpoint":
            w = (self.nu + 1) * self.r_c ** self.nu * np.diff(self.r_f)
        else:
            raise ValueError("weights must be 'volume' or 'midpoint'")
        self.w = w
        self._lu = {}

    def mean(self, q):
        """Volume-weighted mean of one profile ``(n_r,)`` or a stack ``(n_r, k)``."""
        return self.w @ np.asarray(q)

    def surface_flux(self, q):
        """``dq/drho`` at the surface for a stack of profiles ``(n_r, k)``.

        Read with ``compute_boundary_values``, which returns gradients along
        ``+rho``, so this is the flux into the particle per unit diffusivity.
        """
        q = np.asarray(q, dtype=float).reshape(self.n_r, -1)
        if self.flux_from == "cell":
            return (1.0 - q[-1]) / (0.5 * (self.r_f[-1] - self.r_f[-2]))
        _, _, _, grad_right = compute_boundary_values(q, self.r_f, self.r_c, self.bc)
        return np.ravel(grad_right)

    def _factor(self, dt):
        if dt not in self._lu:
            self._lu[dt] = splu((eye_array(self.n_r, format="csc") / dt + self.jac).tocsc())
        return self._lu[dt]

    def step(self, q, dt):
        """One backward-Euler step ``(I/dt + J) q_new = q/dt + g``."""
        return self._factor(dt).solve(q / dt + self.g)

    def march(self, taus, dt, keep=False):
        """March from the empty particle, reading the mean at each of ``taus``.

        Parameters
        ----------
        taus : array_like
            Ascending times; each must be a whole number of steps.
        dt : float
            Step length.
        keep : bool, optional
            Also return every state, shape ``(n_r, n_steps)``, for flux balances.

        Returns
        -------
        dict
            ``mean`` at ``taus``, ``profiles`` shape ``(n_r, len(taus))`` and,
            with ``keep``, ``states``.
        """
        taus = np.atleast_1d(np.asarray(taus, dtype=float))
        steps = np.rint(taus / dt).astype(int)
        if np.any(np.abs(steps * dt - taus) > 1e-9 * taus):
            raise ValueError("every tau must be a whole number of steps")
        q = np.zeros(self.n_r)
        out, states, k = [], [], 0
        for target in steps:
            while k < target:
                q = self.step(q, dt)
                if keep:
                    states.append(q)
                k += 1
            out.append(q.copy())
        prof = np.array(out).T
        res = {"mean": self.mean(prof), "profiles": prof}
        if keep:
            res["states"] = np.array(states).T
        return res

    def spectrum(self):
        """Eigenpairs of the semi-discrete operator, ``J = V diag(lam) V^-1``."""
        if not hasattr(self, "_eig"):
            lam, v = eig(self.jac.toarray())
            order = np.argsort(lam.real)
            self._eig = (lam.real[order], v.real[:, order])
        return self._eig

    def exact_in_time(self, taus):
        """Mean loading of the semi-discrete problem, exact in time.

        The steady state of ``dq/dtau = -J q + g`` is ``q = 1`` (the surface
        value), so ``q(tau) = 1 - V exp(-lam tau) V^-1 1``. What is left is the
        spatial error alone.
        """
        taus = np.atleast_1d(np.asarray(taus, dtype=float))
        lam, v = self.spectrum()
        c = np.linalg.solve(v, np.ones(self.n_r))
        prof = 1.0 - v @ (np.exp(-lam[:, None] * taus[None, :]) * c[:, None])
        return {"mean": self.mean(prof), "profiles": prof}
