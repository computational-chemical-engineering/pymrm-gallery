"""Froment-Bischoff coke deactivation in a fixed bed: exact references and a pymrm bed.

Model code of gallery page B2.2 (Froment and Bischoff 1961). In the canonical
variables of the page, ``x = a z`` (bed depth) and ``T`` (coking time, ``alpha b
eta`` or ``K b eta``), the gas phase is quasi-steady in ``x`` at every ``T``
and the coke is marched in ``T``:

    dy/dx = -g(C) y,   y(0, T) = 1,         dC/dT = r_c(C, y),   C(x, 0) = 0.

Routes to the same solutions:

``pe_exact``, ``ph_exact``
    the paper's closed forms (32)-(33) and (36)-(37);
``ce_exact``
    an exact scalar-ODE reduction of the consecutive case (not in the paper);
``pe_gamma_ode``
    the paper's eq. (29), an ODE in ``x`` exact for any ``gamma``;
``Bed``
    the pymrm finite-volume bed: van Leer TVD convection (``ConvFluxTVD``),
    each gas profile by ``newton``, coke by Heun's method in ``T``.

Contents
--------
pe_exact, ph_exact, ce_exact, pe_gamma_ode
    Exact references.
ei, ce_patch
    The paper's large-eta approximation (46)+(49) of the consecutive case.
cbar_slope
    Distance-averaged coke and its exact local Voorhies exponent.
half_contrast, contrast_of_profile
    The half-bed contrast of a coke profile (the profile diagnostic).
Bed, err_vs_exact
    The pymrm bed and its maximum error against the exact references.

Several functions take an argument whose default is the correct model and
whose alternative is a defect the page's break table injects.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.optimize import brentq
from scipy.sparse import csr_array, vstack
from scipy.sparse.linalg import splu
from scipy.special import expi, lambertw

from pymrm import (ConvFluxTVD, compute_boundary_values, construct_coefficient_matrix,
                   construct_convflux_upwind, construct_div, interp_cntr_to_stagg,
                   interp_cntr_to_stagg_tvd, newton, upwind, vanleer)

__all__ = ["SolveFailed", "pe_exact", "ph_exact", "ce_exact", "pe_gamma_ode", "ei", "ce_patch",
           "cbar_slope", "contrast_of_profile", "half_contrast", "Bed", "err_vs_exact"]


class SolveFailed(RuntimeError):
    """A solve that did not converge (newton or the deferred-correction sweeps)."""


# ----------------------------------------------------------- exact references
def pe_exact(x, T):
    """Parallel mechanism, exponential activity, ``gamma = 0``: eqs. (32)-(33).

    Parameters
    ----------
    x, T : array_like
        Bed depth ``a z`` and coking time ``alpha b eta`` (broadcast).

    Returns
    -------
    y, C : numpy.ndarray
        Reactant mole fraction and ``alpha c``. ``exp(-C)(x, T)`` equals
        ``y(T, x)``: the pair is symmetric under the swap of its arguments.
    """
    x, T = np.broadcast_arrays(np.asarray(x, float), np.asarray(T, float))
    y = 1.0 / (1.0 + np.exp(-T) * np.expm1(x))
    C = np.log1p(np.exp(-x) * np.expm1(T))
    return y, C


def ph_exact(x, T):
    """Parallel mechanism, hyperbolic activity: eqs. (36)-(37).

    ``s = sqrt(1 + 2 T) - 1``, ``W exp(W) = s exp(s - x)`` (Lambert W), and
    ``y = exp(-x + s - W)``, which equals ``W/s``.

    Parameters
    ----------
    x, T : array_like
        Bed depth ``a z`` and coking time ``K b eta`` (broadcast).

    Returns
    -------
    y, W : numpy.ndarray
        Reactant mole fraction and ``W = K c``.
    """
    x, T = np.broadcast_arrays(np.asarray(x, float), np.asarray(T, float))
    s = np.sqrt(1.0 + 2.0 * T) - 1.0
    W = np.real(lambertw(s * np.exp(s - x)))
    return np.exp(-x + s - W), W


def ce_exact(x, T, rtol=1e-12, atol=1e-14):
    """Consecutive mechanism, exponential activity, by an exact scalar reduction.

    At fixed ``x``, ``Chat = -ln y`` obeys ``dChat/dT = 1 - Chat - exp(-Chat)``
    with ``Chat(0) = x``, and ``dC/dT = 1 - exp(-Chat)`` with ``C(0) = 0``
    (integrate ``d(Chat_T)/dx = (exp(-Chat) - 1) Chat_x`` from the inlet, where
    ``Chat = 0``). Integrated with DOP853 for all ``x`` at once.

    Parameters
    ----------
    x : array_like
        Bed depths ``a z``.
    T : array_like
        Coking times ``alpha b eta``, ascending.
    rtol, atol : float, optional
        Integrator tolerances.

    Returns
    -------
    y, C : numpy.ndarray
        Shape ``(len(x), len(T))``.
    """
    x = np.atleast_1d(np.asarray(x, float))
    T = np.atleast_1d(np.asarray(T, float))
    nx = x.size

    def rhs(_, u):
        e = np.exp(-u[:nx])
        return np.concatenate([1.0 - u[:nx] - e, 1.0 - e])

    sol = solve_ivp(rhs, [0.0, max(T.max(), 1e-12)], np.concatenate([x, np.zeros(nx)]),
                    method="DOP853", t_eval=T, rtol=rtol, atol=atol)
    if not sol.success:
        raise SolveFailed(f"consecutive-case reference: {sol.message}")
    return np.exp(-sol.y[:nx]), sol.y[nx:]


def pe_gamma_ode(T_vals, gamma, x_eval, rtol=1e-12):
    """Parallel/exponential case for any ``gamma`` by the paper's eq. (29).

    In canonical form ``C_x = exp(-C) - 1 - gamma C`` with ``C(0) = T``, and
    ``d ln y/dx = -(exp(-C) + gamma)``. One ODE in ``x`` per ``T``; no time
    marching.

    Parameters
    ----------
    T_vals : array_like
        Coking times.
    gamma : float
        ``k2_0/k1_0``.
    x_eval : array_like
        Ascending depths, starting at 0.
    rtol : float, optional
        Integrator tolerance.

    Returns
    -------
    Y, C : numpy.ndarray
        Shape ``(len(T_vals), len(x_eval))``.
    """
    x_eval = np.asarray(x_eval, float)
    T_vals = np.atleast_1d(np.asarray(T_vals, float))
    Y = np.empty((T_vals.size, x_eval.size))
    C = np.empty_like(Y)
    for i, T in enumerate(T_vals):
        sol = solve_ivp(lambda _, u: [np.exp(-u[0]) - 1.0 - gamma * u[0], -(np.exp(-u[0]) + gamma)],
                        [0.0, x_eval[-1]], [T, 0.0], method="DOP853", t_eval=x_eval,
                        rtol=rtol, atol=1e-14)
        if not sol.success:
            raise SolveFailed(f"eq. (29) route: {sol.message}")
        C[i], Y[i] = sol.y[0], np.exp(sol.y[1])
    return Y, C


# ------------------------------------------------- the paper's large-eta patch
def ei(u):
    """The paper's ``ei(u) = int_u^inf exp(-l)/l dl = -Ei(-u)``."""
    return -expi(-u)


def ce_patch(x, T):
    """The large-eta approximation (46)+(49) of the consecutive case.

    Solves ``ei(T y) = ei(T) + x exp(-T)`` for ``y`` and sets ``C = (1 - y) T``.

    Parameters
    ----------
    x : float
        The depth that enters eq. (49): ``a z``, the reading eq. (48)
        requires (the page passes ``alpha z`` for the reading as printed).
    T : float
        ``alpha b eta``.

    Returns
    -------
    y, C : float
    """
    y = brentq(lambda yy: ei(T * yy) - (ei(T) + x * np.exp(-T)), 1e-300, 1.0, xtol=1e-15,
               rtol=4 * np.finfo(float).eps)
    return y, (1.0 - y) * T


# ------------------------------------------------- averages and the exponent
def _profiles(mech, xg, T):
    """``C`` and ``dC/dT`` on ``xg`` (rows) at the times ``T`` (columns)."""
    T = np.atleast_1d(np.asarray(T, float))
    if mech == "PE":
        y, C = pe_exact(xg[:, None], T[None, :])
        return C, y
    if mech == "PH":
        y, W = ph_exact(xg[:, None], T[None, :])
        return W, y / (1.0 + W)
    y, C = ce_exact(xg, T)
    return C, 1.0 - y


def cbar_slope(mech, aL, T, nx=801):
    """Distance-averaged coke and its exact local Voorhies exponent.

    ``cbar = (1/aL) int_0^aL C dx`` (eqs. 51-52) and the exponent
    ``d ln cbar/d ln T = T cbar'(T)/cbar``, where ``cbar'`` is the average of the
    coke rate ``dC/dT`` (``y`` for PE, ``y/(1+W)`` for PH, ``1 - y`` for CE):
    no difference quotient.

    Parameters
    ----------
    mech : {"PE", "PH", "CE"}
        Case.
    aL : float
        Bed length ``a L``.
    T : array_like
        Ascending coking times.
    nx : int, optional
        Simpson nodes along the bed (odd).

    Returns
    -------
    cbar, slope : numpy.ndarray
    """
    T = np.atleast_1d(np.asarray(T, float))
    xg = np.linspace(0.0, aL, nx)
    C, rate = _profiles(mech, xg, T)
    cb = simpson(C, x=xg, axis=0) / aL
    return cb, T * simpson(rate, x=xg, axis=0) / aL / cb


# ----------------------------------------------------- the profile diagnostic
def contrast_of_profile(C, x, rule="simpson"):
    """Half-bed contrast ``(rear-half mean - front-half mean)/cbar``.

    Parameters
    ----------
    C : numpy.ndarray
        Coke profile.
    x : numpy.ndarray
        Its abscissae: Simpson nodes from 0 to ``aL`` (odd count) for
        ``rule="simpson"``, cell centres of a uniform grid (even count) for
        ``rule="cells"``.
    rule : {"simpson", "cells"}, optional

    Returns
    -------
    float
    """
    if rule == "cells":
        h = C.size // 2
        return float((C[h:].mean() - C[:h].mean()) / C.mean())
    h = C.size // 2
    front = simpson(C[:h + 1], x=x[:h + 1])
    rear = simpson(C[h:], x=x[h:])
    return float((rear - front) / (0.5 * (rear + front)))


def half_contrast(mech, aL, T, nx=801):
    """Half-bed contrast of the exact coke profile at time ``T``.

    Parameters
    ----------
    mech : {"PE", "CE"}
        Case (``"PE"`` in place of ``"CE"`` is the mechanism swap the break
        table injects).
    aL, T : float
        Bed length and coking time.
    nx : int, optional
        Simpson nodes (odd).

    Returns
    -------
    float
    """
    xg = np.linspace(0.0, aL, nx)
    C, _ = _profiles(mech, xg, [T])
    return contrast_of_profile(C[:, 0], xg)


# --------------------------------------------------------------- the pymrm bed
class Bed:
    """Quasi-steady gas phase on a finite-volume grid; coke marched in ``T``.

    Parameters
    ----------
    mech : {"PE", "PH", "CE", "MIX"}
        Parallel/exponential (with ``gamma``), parallel/hyperbolic (``C`` holds
        ``W = K c``), consecutive/exponential, or exponential activity with the
        mixed coke source ``(1 - w) y + w (1 - y)``.
    L : float, optional
        Bed length ``a L``.
    n : int, optional
        Number of cells.
    gamma : float, optional
        ``k2_0/k1_0`` (PE only).
    w_mix : float, optional
        Consecutive weight of the mixed source.
    limiter : callable or None, optional
        TVD limiter; None is first-order upwind.
    outlet : {"extrapolated", "outflow", "zero_gradient"}, optional
        Value carried by the outlet face: the linear extrapolation of the last
        two cells (default); the last cell value (``{"outflow": True}``); or
        the zero-gradient reconstruction of ``{"a": 1, "b": 0, "d": 0}``.
    route : {"newton", "deferred"}, optional
        ``"newton"``: ``ConvFluxTVD`` fluxes and Jacobian, ``newton``.
        ``"deferred"``: the same discrete equations by van Leer deferred
        correction, swept to a fixed point (a validation route).
    y_in : float, optional
        Feed mole fraction.
    activity : callable, optional
        Replaces the activity ``phi(C)`` of the mechanism.
    tol : float, optional
        Default tolerance of :meth:`solve_y`.
    """

    def __init__(self, mech, L=2.0, n=200, gamma=0.0, w_mix=0.0, limiter=vanleer,
                 outlet="extrapolated", route="newton", y_in=1.0, activity=None, tol=1e-12):
        self.mech, self.L, self.n, self.tol = mech, float(L), int(n), float(tol)
        self.gamma, self.w_mix, self.route = float(gamma), float(w_mix), route
        self.shape = (self.n, 1)
        self.x_f = np.linspace(0.0, self.L, self.n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self._activity = activity
        # bc: a dy/dn + b y = d with n the OUTWARD normal
        # inlet x = 0: y = y_in (pure feed)                               -> a = 0, b = 1, d = y_in
        bc_in = {"a": 0.0, "b": 1.0, "d": float(y_in)}
        if outlet == "zero_gradient":
            # outlet x = L: dy/dx = 0, dy/dn = dy/dx                       -> a = 1, b = 0, d = 0
            bc_out = {"a": 1.0, "b": 0.0, "d": 0.0}
        else:
            # outlet x = L: pure outflow, the face carries the last cell value
            bc_out = {"outflow": True}
        self.bc = (bc_in, bc_out)
        self.div = construct_div(self.shape, self.x_f, nu=0)            # nu = 0: Cartesian bed
        self.tvd = ConvFluxTVD(self.shape, self.x_f, self.x_c, self.bc,
                               tvd_limiter=upwind if limiter is None else limiter)
        self.limiter = limiter
        # The outlet face value. With {"outflow": True} the face carries the last
        # cell value, half a cell upstream of the face: an O(h) error where y has
        # a gradient there. "extrapolated" adds the constant correction from that
        # value to the linear extrapolation of the last two cells (the last row
        # of interp_cntr_to_stagg without a bc), assembled once.
        self.K = csr_array((self.n + 1, self.n))
        if outlet == "extrapolated":
            # outlet-face weights of the last two cells (unit columns interpolated)
            w = np.zeros(self.n)
            w[-2:] = interp_cntr_to_stagg(np.eye(self.n, 2, k=2 - self.n), self.x_f, self.x_c, axis=0)[-1]
            w[-1] -= 1.0
            self.K = vstack([csr_array((self.n, self.n)), csr_array(w[np.newaxis, :])], format="csr")
        self.outlet = outlet
        self.divK = (self.div @ self.K).tocsc()
        if route == "deferred":
            conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c, self.bc, v=1.0)
            self.L_up = (self.div @ conv).tocsc() + self.divK
            self.b_up = (self.div @ conv_bc).toarray().ravel()

    # -- physics ---------------------------------------------------------------
    def activity(self, C):
        """Activity ``phi(C)``: ``1/(1 + W)`` for PH, ``exp(-C)`` otherwise."""
        if self._activity is not None:
            return self._activity(C)
        return 1.0 / (1.0 + C) if self.mech == "PH" else np.exp(-C)

    def gcoef(self, C):
        """``g`` in ``dy/dx = -g y``: linear in ``y`` at frozen ``C``."""
        g = self.activity(C)
        return g + self.gamma if self.mech == "PE" else g

    def coke_rate(self, C, y):
        """``dC/dT``: ``y`` (PE), ``y/(1+W)`` (PH), ``1 - y`` (CE), mixed (MIX)."""
        if self.mech == "CE":
            return 1.0 - y
        if self.mech == "MIX":
            return (1.0 - self.w_mix) * y + self.w_mix * (1.0 - y)
        return y * self.activity(C) if self.mech == "PH" else y

    # -- the gas phase ---------------------------------------------------------
    def residual(self, y, C, mass=0.0, y_old=None):
        """Residual and Jacobian of ``mass (y - y_old) + d(F(y))/dx + g(C) y``.

        ``F`` are the van Leer face fluxes of ``ConvFluxTVD`` (velocity 1) with
        the outlet face set by ``outlet``; ``mass`` is ``eps_s/dt`` in a
        backward-Euler step of the full transient, 0 for the quasi-steady gas.
        """
        y = np.asarray(y, float).reshape(self.shape)
        flux, jac = self.tvd(y, 1.0)
        g = self.div @ np.asarray(flux).ravel() + self.divK @ y.ravel() + self.gcoef(C) * y.ravel()
        diag = self.gcoef(C) + mass
        if mass:
            g = g + mass * (y.ravel() - y_old)
        J = self.div @ jac + self.divK + construct_coefficient_matrix(diag.reshape(self.shape))
        return g, J.tocsc()

    def solve_y(self, C, y0=None, mass=0.0, y_old=None, tol=None):
        """Gas profile at frozen coke ``C``.

        Parameters
        ----------
        C : numpy.ndarray
            Coke on the cells.
        y0 : numpy.ndarray, optional
            Initial guess (default: ``exp(-int g dx)`` at the cell centres).
        mass, y_old : float, numpy.ndarray, optional
            Backward-Euler term of the full transient.
        tol : float, optional
            Newton step tolerance, or the largest change of a deferred sweep
            (default: the bed's ``tol``).

        Returns
        -------
        numpy.ndarray
            ``y`` on the cells.

        Raises
        ------
        SolveFailed
            When ``newton`` reports no success, or the sweeps do not settle.
        """
        tol = self.tol if tol is None else tol
        g = self.gcoef(C)
        if y0 is None:
            h = np.diff(self.x_f)
            y0 = np.exp(-(np.cumsum(g * h) - 0.5 * g * h))
        if self.route == "deferred":
            return self._deferred(C, y0, mass, y_old, tol)
        sol = newton(lambda y: self.residual(y, C, mass, y_old), np.reshape(y0, self.shape),
                     tol=tol, damping="monotone", stagnation=True)
        if not sol.success:
            raise SolveFailed(f"gas profile, {self.mech}, n = {self.n}: {sol.message}")
        self.last_nit = sol.nit
        return np.asarray(sol.x).ravel()

    def _deferred(self, C, y0, mass, y_old, tol, max_it=300):
        """Upwind operator implicit, van Leer correction of the previous sweep explicit."""
        diag = self.gcoef(C) + mass
        lu = splu((self.L_up + construct_coefficient_matrix(diag.reshape(self.shape))).tocsc())
        rhs0 = -self.b_up + (mass * y_old if mass else 0.0)
        y = np.asarray(y0, float).ravel()
        for sweep in range(1, max_it + 1):
            _, dface = interp_cntr_to_stagg_tvd(y.reshape(self.shape), self.x_f, self.x_c, self.bc,
                                                1.0, tvd_limiter=self.limiter, axis=0)
            y_new = lu.solve(rhs0 - self.div @ np.asarray(dface).ravel())
            change = np.max(np.abs(y_new - y))
            y = y_new
            if change < tol:
                self.last_nit = sweep
                return y
        raise SolveFailed(f"deferred correction: change {change:.1e} after {max_it} sweeps")

    def outlet_value(self, y):
        """The outlet mole fraction the operators transported (the outlet face value)."""
        if self.outlet == "extrapolated":
            return float(interp_cntr_to_stagg(np.asarray(y, float).reshape(self.shape),
                                              self.x_f, self.x_c, axis=0)[-1, 0])
        return float(np.ravel(compute_boundary_values(np.asarray(y, float).reshape(self.shape),
                                                      self.x_f, self.x_c, self.bc)[2])[0])

    # -- coking time -----------------------------------------------------------
    def march(self, T_end, m, record_T=(), scheme="heun", on_grid=True):
        """March ``C`` from a fresh bed to ``T_end`` in ``m`` equal steps.

        Parameters
        ----------
        T_end : float
            Final coking time.
        m : int
            Number of steps.
        record_T : sequence of float, optional
            Times at which ``(y, C)`` is recorded; each must lie on the step
            grid (an assertion), never at the first step past it.
        scheme : {"heun", "euler"}, optional
            Heun's method (second order) or forward Euler.
        on_grid : bool, optional
            False records each time at the first step at or past it, under
            the target time's key (a defect the page injects).

        Returns
        -------
        dict
            ``{T: (y, C)}`` for the recorded times. ``self.C`` holds the final
            coke profile.
        """
        dt = T_end / m
        C = np.zeros(self.n)
        rec, k, y = {}, 0, None
        record_T = np.asarray(sorted(record_T), float)
        for j in range(m):
            y = self.solve_y(C, y)
            f1 = self.coke_rate(C, y)
            if scheme == "euler":
                C = C + dt * f1
            else:
                y2 = self.solve_y(C + dt * f1, y)
                C = C + 0.5 * dt * (f1 + self.coke_rate(C + dt * f1, y2))
            t_now = (j + 1) * dt
            while k < record_T.size and record_T[k] <= t_now + 1e-9 * dt:
                assert not on_grid or abs(record_T[k] - t_now) < 1e-9, "record time off the step grid"
                rec[float(record_T[k])] = (self.solve_y(C, y), C.copy())
                k += 1
        self.C = C
        return rec

    def march_transient(self, eps, y, C, T_span, m, tol=1e-12):
        """The full transient ``eps dy/dT + dy/dx = -g(C) y``, ``dC/dT = r_c``.

        Backward Euler in ``y`` (one ``newton`` solve per step) and
        ``C <- C + dt r_c(C, (y_old + y_new)/2)``: first order in the step.
        With ``route="deferred"`` and ``tol=numpy.inf`` each step is a single
        sweep with the van Leer correction taken from the previous step's
        profile (the correction lagged by one step).

        Parameters
        ----------
        eps : float
            ``eps_s = alpha b/a``, coking time per unit of gas transit.
        y, C : numpy.ndarray
            Initial state on the cells.
        T_span : float
            Duration.
        m : int
            Number of steps.
        tol : float, optional
            Passed to :meth:`solve_y`.

        Returns
        -------
        y, C : numpy.ndarray
        """
        dt = T_span / m
        y, C = np.array(y, float), np.array(C, float)
        for _ in range(m):
            y_new = self.solve_y(C, y, mass=eps / dt, y_old=y, tol=tol)
            C = C + dt * self.coke_rate(C, 0.5 * (y + y_new))
            y = y_new
        return y, C


def err_vs_exact(mech, T_end, n, m, nrec=5, L=2.0, on_grid=True, return_bed=False, **kw):
    """Maximum error of the pymrm bed against the exact solution.

    The maximum over the cells and over ``nrec`` equally spaced record times
    up to ``T_end``, in ``y`` and in ``C``.

    Parameters
    ----------
    mech : {"PE", "PH", "CE"}
    T_end : float
    n, m : int
        Cells and coking steps (``m`` a multiple of ``nrec``).
    nrec : int, optional
    L : float, optional
    on_grid : bool, optional
        Passed to :meth:`Bed.march`.
    return_bed : bool, optional
        Also return the bed after the march.
    **kw
        Passed to :class:`Bed` (``scheme`` to :meth:`Bed.march`).

    Returns
    -------
    err_y, err_C : float
        And the bed when ``return_bed``.
    """
    scheme = kw.pop("scheme", "heun")
    bed = Bed(mech, L=L, n=n, **kw)
    times = np.linspace(T_end / nrec, T_end, nrec)
    rec = bed.march(T_end, m, record_T=times, scheme=scheme, on_grid=on_grid)
    ey = eC = 0.0
    for T, (y, C) in rec.items():
        if mech == "PE":
            ye, Ce = pe_exact(bed.x_c, T)
        elif mech == "PH":
            ye, Ce = ph_exact(bed.x_c, T)
        else:
            ye, Ce = ce_exact(bed.x_c, [T])
            ye, Ce = ye[:, 0], Ce[:, 0]
        ey, eC = max(ey, np.abs(y - ye).max()), max(eC, np.abs(C - Ce).max())
    return (ey, eC, bed) if return_bed else (ey, eC)
