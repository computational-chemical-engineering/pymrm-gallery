"""The stagnant-film model of a second-order gas-liquid reaction (page F3.2).

Froment, De Wilde and Bischoff (2011), section 6.3: gas A dissolves into a
liquid film that holds a non-volatile reactant B, ``aA + bB -> products``,
first order in each. In ``xi = y/y_L``, ``alpha = C_A/C_Ai``,
``beta = C_B/C_Bb``::

    alpha'' = gamma^2 alpha beta,      beta'' = (gamma^2/q) alpha beta,
    alpha(0) = 1, beta'(0) = 0,        alpha(1) = 0, beta(1) = 1,

and the enhancement factor is ``F_A = -alpha'(0)``.

Three ways to the enhancement factor, kept apart:

- :class:`Film`, the film equations by finite volumes (pymrm operators,
  ``newton`` with ``NumJac``), and :class:`Extrapolated`, the Richardson
  extrapolation of one graded film on two grids;
- :func:`f_a_bvp`, the same two-point problem by adaptive collocation
  (``scipy.integrate.solve_bvp``), sharing no assembly, Jacobian or mesh;
- :func:`f_a_vkh` and :func:`f_a_vkh_route_b`, the Van Krevelen-Hoftijzer
  implicit relation (6.3.5-1) solved for ``F_A`` and for ``gamma'``.

The extremum and root helpers at the end locate the worst point of the
approximation and the depletion thresholds.
"""
import functools

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
from pymrm import (NumJac, compute_boundary_values, construct_div,
                   construct_grad, newton)
from pymrm.checks import find_roots

#: Brent tolerances for the implicit relation and the extremum roots.
XTOL = dict(xtol=1e-13, rtol=8.9e-16)


# ------------------------------------------- Froment (6.3.2-11) + (6.3.5-1)
def f_a_vkh(gamma, q, sign=-1.0, offset=1.0, denom="q", xtol=None):
    """Van Krevelen-Hoftijzer: ``F_A = gamma'/tanh(gamma')`` solved for ``F_A``.

    Parameters
    ----------
    gamma, q : float
        The two groups of the film problem.
    sign, offset, denom : float, float, str
        The printed form is ``gamma' = gamma sqrt(1 - (F_A - 1)/q)``: ``sign``
        -1, ``offset`` 1, ``denom`` ``"q"``. The break table mistranscribes it
        through these (``denom="ei"`` divides by ``1 + q``).
    xtol : float, optional
        Brent ``xtol``; ``None`` uses :data:`XTOL`.

    Returns
    -------
    float
        The enhancement factor of the approximation.
    """
    d = q if denom == "q" else 1.0 + q
    tol = XTOL if xtol is None else dict(xtol=xtol, rtol=8.9e-16)

    def resid(F):
        s = np.sqrt(max(1.0 + sign * (F - offset) / d, 0.0))
        g = gamma * s
        return F - (g / np.tanh(g) if g > 1e-8 else 1.0 + g * g / 3.0)

    hi = 1.0 + q
    while resid(hi) < 0 and hi < 1e12:      # a broken form can exceed the ceiling
        hi *= 2.0
    return brentq(resid, 1.0, hi, **tol, maxiter=300)


def f_a_vkh_route_b(gamma, q):
    """The same two printed equations, solved for ``gamma'`` instead of ``F_A``.

    ``F = u/tanh u`` and ``F = 1 + q (1 - u^2/gamma^2)`` with ``u = gamma'``;
    eliminating ``F`` leaves one scalar equation in ``u``.

    Parameters
    ----------
    gamma, q : float
        The two groups of the film problem.

    Returns
    -------
    float
        The enhancement factor of the approximation.
    """
    def resid(u):
        return u / np.tanh(u) - (1.0 + q * (1.0 - (u / gamma) ** 2))

    u = brentq(resid, 1e-12, gamma, **XTOL, maxiter=300)
    return u / np.tanh(u)


# ------------------------------------------------------- the film, in pymrm
def faces(n, stretch):
    """Face positions on ``[0, 1]``, graded towards ``xi = 0`` (the interface).

    Parameters
    ----------
    n : int
        Number of cells.
    stretch : float
        Grading ``s`` of ``x_f = (exp(s t) - 1)/(exp(s) - 1)``; 0 gives a
        uniform grid.

    Returns
    -------
    numpy.ndarray
        ``n + 1`` face positions.
    """
    t = np.linspace(0.0, 1.0, n + 1)
    return t if stretch <= 1e-12 else np.expm1(stretch * t) / np.expm1(stretch)


class Film:
    """Stagnant film on ``xi`` in ``[0, 1]``, fields ``(alpha, beta)`` = (A, B).

    Parameters
    ----------
    n_x : int
        Number of cells.
    stretch : float
        Grid grading towards the interface (see :func:`faces`).
    bc_beta_bulk : {"dirichlet", "neumann"}
        Condition on ``beta`` at the bulk edge: ``beta(1) = 1`` (the model), or
        ``dbeta/dxi = 0`` there (a break-table defect that cuts the supply of B).
    """

    def __init__(self, n_x=800, stretch=8.0, bc_beta_bulk="dirichlet"):
        self.n_x, self.stretch = n_x, stretch
        self._memo = {}
        self.shape = (n_x, 2)                  # spatial axis first, then (alpha, beta)
        self.x_f = faces(n_x, stretch)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # a dc/dn + b c = d, n the OUTWARD normal (n = -xi at xi = 0, +xi at xi = 1);
        # one column per field, [alpha, beta].
        # xi = 0, the interface: alpha = 1 (saturated, C_A = C_Ai) and
        # -dbeta/dxi = 0 (B non-volatile, no flux of B through the interface).
        left = {"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[1.0, 0.0]]}
        if bc_beta_bulk == "dirichlet":
            # xi = 1, the bulk edge: alpha = 0 (C_Ab = 0, no reaction in the bulk)
            # and beta = 1 (B at its bulk value C_Bb).
            right = {"a": [[0.0, 0.0]], "b": [[1.0, 1.0]], "d": [[0.0, 1.0]]}
        else:
            # xi = 1, the break-table variant: alpha = 0 and dbeta/dxi = 0, a
            # zero-gradient bulk edge, so no B enters the film at all.
            right = {"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[0.0, 0.0]]}
        self.bc = (left, right)
        grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c,
                                           self.bc, axis=0)
        div_mat = construct_div(self.shape, self.x_f, nu=0, axis=0)   # slab
        # the constant diffusion operator, assembled once
        self.flux_mat, self.flux_bc = -grad_mat, -grad_bc
        self.jac_diff = div_mat @ self.flux_mat
        self.g_diff = (div_mat @ self.flux_bc).toarray().ravel()
        # F_A = -dalpha/dxi on the interface face: row 0 of the face-flux operator
        # (face 0, field alpha), boundary term included, so the flux the solve
        # transports through xi = 0 and not a cell-centre difference
        self.iface = np.asarray(self.flux_mat[[0], :].todense()).ravel()
        self.iface_bc = float(self.flux_bc.toarray().ravel()[0])
        # (n_x, 2), never a bare (n_x,): the rate is pointwise in the two fields
        self.numjac = NumJac((n_x, 2))

    def enhancement(self, u):
        """``F_A``, the interfacial flux of A, from the face-flux operator.

        Parameters
        ----------
        u : numpy.ndarray
            Solution, shape ``(n_x, 2)``.

        Returns
        -------
        float
        """
        return float(self.iface @ np.asarray(u).ravel() + self.iface_bc)

    def solve(self, gamma, q, frozen_b=False, tol=1e-10):
        """Steady film by Newton from a closed-form cold start.

        Parameters
        ----------
        gamma, q : float
            The two groups of the film problem.
        frozen_b : bool
            Break-table defect: B held at its bulk value in the rate and its
            consumption dropped, so ``beta = 1`` and the film is
            pseudo-first-order.
        tol : float
            Newton step tolerance on unknowns of order one.

        Returns
        -------
        numpy.ndarray
            ``(alpha, beta)`` at the cell centres, shape ``(n_x, 2)``.

        Raises
        ------
        RuntimeError
            If Newton does not converge.
        """
        def rate(u):
            b = np.ones_like(u[:, 1]) if frozen_b else np.clip(u[:, 1], 0, None)
            r = np.clip(u[:, 0], 0, None) * b
            return np.stack([gamma**2 * r,
                             np.zeros_like(r) if frozen_b else (gamma**2 / q) * r],
                            axis=1)

        # cold start: the pseudo-first-order profile at this same gamma
        g = min(gamma, 700.0)                      # sinh overflows past ~710
        u0 = np.stack([np.sinh(g * (1.0 - self.x_c)) / np.sinh(g),
                       np.ones_like(self.x_c)], axis=1)

        def residual(u):
            g_r, jac_r = self.numjac(rate, u)
            return (self.jac_diff @ u.reshape((-1, 1))
                    + self.g_diff.reshape((-1, 1))
                    + g_r.reshape((-1, 1))), self.jac_diff + jac_r

        res = newton(residual, u0, tol=tol, maxfev=200,
                     bounds=(0.0, np.inf), stagnation=True)
        if not res.success:
            raise RuntimeError(f"film solve failed at gamma={gamma}, q={q}: "
                               f"{res.message}")
        return np.asarray(res.x).reshape(self.shape)

    def beta_interface(self, u):
        """``beta(0)`` on the interface face, by ``compute_boundary_values``.

        Parameters
        ----------
        u : numpy.ndarray
            Solution, shape ``(n_x, 2)``.

        Returns
        -------
        float
        """
        v, _ = compute_boundary_values(u, self.x_f, self.x_c, self.bc[0], axis=0,
                                       bound_id=0)
        return float(np.asarray(v).ravel()[1])

    def beta_constancy(self, gamma, q, **kw):
        """Rate-weighted mean of ``beta`` over the film, divided by ``beta(0)``.

        ``int alpha beta^2 dxi / int alpha beta dxi / beta(0)``: 1 exactly when
        ``beta`` does not vary where the reaction happens (book p. 340's
        wording of the assumption). ``nan`` when there is no reaction to
        average over (no B left, as under the zero-gradient bulk defect).

        Parameters
        ----------
        gamma, q : float
            The two groups of the film problem.
        **kw
            Passed to :meth:`solve`.

        Returns
        -------
        float
        """
        key = ("constancy", gamma, q, tuple(sorted(kw.items())))
        if key not in self._memo:          # a pure function of its arguments
            u = self.solve(gamma, q, **kw)
            a = np.clip(u[:, 0], 0.0, None)
            b = np.clip(u[:, 1], 0.0, None)
            w = np.diff(self.x_f) * a * b
            s, b0 = float(w.sum()), self.beta_interface(u)
            self._memo[key] = (float("nan") if not (s > 0.0) or b0 == 0.0
                               else float((w * b).sum() / s) / b0)
        return self._memo[key]

    def f_a(self, gamma, q, **kw):
        """``(F_A, beta(0))`` at one point, memoised on ``(gamma, q, kw)``.

        Parameters
        ----------
        gamma, q : float
            The two groups of the film problem.
        **kw
            Passed to :meth:`solve`.

        Returns
        -------
        tuple of float
        """
        key = (gamma, q, tuple(sorted(kw.items())))
        if key not in self._memo:
            u = self.solve(gamma, q, **kw)
            self._memo[key] = (self.enhancement(u), self.beta_interface(u))
        return self._memo[key]


class Extrapolated:
    """Richardson extrapolation of the same film on two grids.

    ``F(2n) + [F(2n) - F(n)]/(2^p - 1)`` with the observed order ``p = 2``.

    Parameters
    ----------
    coarse, fine : Film
        The same grading on ``n`` and ``2n`` cells.
    order : float
        Order of the discretisation error removed.
    """

    def __init__(self, coarse, fine, order=2.0):
        self.coarse, self.fine = coarse, fine
        self.n_x, self.stretch = fine.n_x, fine.stretch
        self.w = 1.0 / (2.0**order - 1.0)

    def f_a(self, gamma, q, **kw):
        """Extrapolated ``(F_A, beta(0))``; see :meth:`Film.f_a`."""
        f_c, b_c = self.coarse.f_a(gamma, q, **kw)
        f_f, b_f = self.fine.f_a(gamma, q, **kw)
        return f_f + (f_f - f_c) * self.w, b_f + (b_f - b_c) * self.w

    def solve(self, gamma, q, **kw):
        """The fine grid's profile; only scalars are extrapolated."""
        return self.fine.solve(gamma, q, **kw)

    def beta_constancy(self, gamma, q, **kw):
        """Extrapolated constancy ratio; see :meth:`Film.beta_constancy`."""
        c = self.coarse.beta_constancy(gamma, q, **kw)
        f = self.fine.beta_constancy(gamma, q, **kw)
        return f + (f - c) * self.w


_GRIDS = {}


def graded_film(n_x, stretch):
    """One :class:`Film` per ``(n_x, stretch)``, so its solve memo is shared.

    Parameters
    ----------
    n_x : int
    stretch : float

    Returns
    -------
    Film
    """
    return _GRIDS.setdefault((n_x, float(stretch)), Film(n_x=n_x, stretch=stretch))


# ------------------------------------------------ the film, by collocation
@functools.lru_cache(maxsize=None)
def f_a_bvp(gamma, q, tol=1e-10, n0=200):
    """The same two-point problem by adaptive collocation (``solve_bvp``).

    Parameters
    ----------
    gamma, q : float
        The two groups of the film problem.
    tol : float
        ``solve_bvp`` tolerance.
    n0 : int
        Initial mesh size (graded towards the interface).

    Returns
    -------
    tuple of float
        ``(F_A, beta(0))``.

    Raises
    ------
    RuntimeError
        If ``solve_bvp`` does not converge.
    """
    def fun(x, y):
        r = gamma**2 * y[0] * y[2]
        return np.vstack([y[1], r, y[3], r / q])

    def bc(ya, yb):
        return np.array([ya[0] - 1.0, yb[0], ya[3], yb[2] - 1.0])

    x = np.expm1(6.0 * np.linspace(0, 1, n0)) / np.expm1(6.0)
    y = np.zeros((4, x.size))
    y[0], y[1], y[2] = 1.0 - x, -1.0, 1.0
    s = solve_bvp(fun, bc, x, y, tol=tol, max_nodes=400000)
    if s.status != 0:
        raise RuntimeError(f"solve_bvp failed at gamma={gamma}, q={q}: {s.message}")
    return float(-s.y[1][0]), float(s.sol(0.0)[2])


# -------------------------------------------------------- root-find helpers
def argext_log(f, lo, hi, n_scan=13, h=1e-3, rootfind=True):
    """The extremum of ``|f|`` on ``[lo, hi]``, root-found from ``f' = 0``.

    A scan of ``n_scan`` log-spaced points over the whole window brackets the
    largest ``|f|``; the extremum is then the Brent root of the central
    difference ``f(x e^h) - f(x e^-h)`` between the scan neighbours.

    Parameters
    ----------
    f : callable
    lo, hi : float
        The window, scanned on a logarithmic grid.
    n_scan : int
        Scan points.
    h : float
        Step of the central difference in ``ln x``.
    rootfind : bool
        False returns the largest sampled value (a break-table defect).

    Returns
    -------
    tuple
        ``(x, f(x), root_found)``. ``root_found`` is False when the largest
        sampled value sits on a window end (the extremum is at or beyond it)
        or the bracket holds no sign change.
    """
    xs = np.geomspace(lo, hi, n_scan)
    fs = np.array([f(x) for x in xs])
    i = int(np.argmax(np.abs(fs)))
    if not rootfind or i in (0, n_scan - 1):
        return float(xs[i]), float(fs[i]), False

    def dfd(lx):
        x = np.exp(lx)
        return f(x * np.exp(h)) - f(x * np.exp(-h))

    a, b = np.log(xs[i - 1]), np.log(xs[i + 1])
    if dfd(a) * dfd(b) >= 0:
        return float(xs[i]), float(fs[i]), False
    lx = brentq(dfd, a, b, **XTOL, maxiter=200)
    return float(np.exp(lx)), float(f(np.exp(lx))), True


def root_log(f, lo, hi):
    """The root of ``f`` on ``[lo, hi]``, Brent in ``ln x``; ``f`` must change sign.

    Parameters
    ----------
    f : callable
    lo, hi : float

    Returns
    -------
    float
    """
    return float(np.exp(brentq(lambda lx: f(np.exp(lx)), np.log(lo), np.log(hi),
                               **XTOL, maxiter=200)))


def unique_root_log(f, lo, hi, n_scan=13):
    """The single root of ``f`` on the physical range ``[lo, hi]``.

    ``pymrm.checks.find_roots`` scans the whole range on a log grid; one sign
    change is refined by :func:`root_log`.

    Parameters
    ----------
    f : callable
    lo, hi : float
    n_scan : int

    Returns
    -------
    float
        The root, or ``nan`` when the scan finds no sign change on the range.

    Raises
    ------
    RuntimeError
        If the scan finds more than one root, or cannot evaluate ``f``.
    """
    found = find_roots(f, lo, hi, n_scan=n_scan, log=True)
    if found["failed"]:
        raise RuntimeError(f"root scan could not evaluate f at {found['failed']}")
    roots = found["roots"]
    if not roots:
        return float("nan")
    if len(roots) > 1:
        raise RuntimeError(f"{len(roots)} roots on [{lo}, {hi}]: {roots}")
    xs = np.geomspace(lo, hi, n_scan)
    j = int(np.searchsorted(xs, roots[0]))
    return root_log(f, xs[max(j - 1, 0)], xs[min(j, n_scan - 1)])
