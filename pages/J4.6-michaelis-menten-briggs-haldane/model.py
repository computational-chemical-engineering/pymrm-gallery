"""Models of page J4.6: the enzyme scheme (1.5.1-1) of Froment, De Wilde & Bischoff.

Three models of one well-mixed batch:

* :class:`BatchMassAction`, the scheme with no approximation, four fields marched
  by backward Euler with ``pymrm.newton`` and a ``NumJac((1, 4))`` Jacobian;
* :class:`BatchReduced`, the rate law (1.5.1-15) marched the same way;
* :func:`s_closed`, the closed-form integral of (1.5.1-15) through Lambert W.

:func:`mass_action` integrates the scheme with ``scipy`` as the reference, and
:func:`max_gap` locates the largest difference of two solutions in time as a
root of the analytic time derivative of their difference.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import identity as speye
from scipy.special import lambertw

from pymrm import NumJac, newton
from pymrm.checks import find_roots


def s_closed(t, s0, rm, KM):
    """Closed-form substrate concentration of the rate law (1.5.1-15) in a batch.

    Parameters
    ----------
    t : float or array_like
        Time, min.
    s0 : float
        Initial substrate concentration, cmol/L.
    rm, KM : float
        Maximum rate r_m (cmol/(L min)) and K_M (cmol/L).

    Returns
    -------
    ndarray
        C_A(t) from K_M ln(s0/s) + (s0 - s) = r_m t, inverted with the principal
        branch of Lambert W.
    """
    z = (s0/KM)*np.exp((s0 - rm*np.asarray(t, float))/KM)
    return KM*np.real(lambertw(z, 0))


class ReducedCurve:
    """C_A(t) of the rate law (1.5.1-15), with its time derivative.

    Parameters
    ----------
    s0, rm, KM : float
        Initial substrate (cmol/L), r_m (cmol/(L min)) and K_M (cmol/L).
    """

    def __init__(self, s0, rm, KM):
        self.s0, self.rm, self.KM = s0, rm, KM

    def __call__(self, t):
        return s_closed(t, self.s0, self.rm, self.KM)

    def deriv(self, t):
        """dC_A/dt = -r_m C_A/(K_M + C_A), the rate law itself."""
        s = self(t)
        return -self.rm*s/(self.KM + s)


def elem(rho, CE0, rm, KM):
    """Elementary constants from the observable ones.

    Parameters
    ----------
    rho : float
        k_2/k_-1.
    CE0 : float
        Total enzyme, cmol/L.
    rm, KM : float
        r_m = k_2 C_E0 (cmol/(L min)) and K_M = (k_-1 + k_2)/k_1 (cmol/L).

    Returns
    -------
    tuple of float
        (k_1, k_-1, k_2); the map is bijective at fixed C_E0.
    """
    kc = rm/CE0
    kr = kc/rho
    kf = (kr + kc)/KM
    return kf, kr, kc


class DenseMassAction:
    """Dense solution of the scheme (1.5.1-1), with the analytic time derivative.

    Parameters
    ----------
    sol : OdeSolution
        Dense output of ``solve_ivp`` for the state (C_A, C_A-E).
    p : tuple of float
        (k_1, k_-1, k_2).
    CE0 : float
        Total enzyme, cmol/L.
    total : bool
        Return C_A + C_A-E instead of the free substrate C_A.
    """

    def __init__(self, sol, p, CE0, total):
        self.sol, self.p, self.CE0, self.total = sol, p, CE0, total

    def __call__(self, t):
        y = self.sol(t)
        return y[0] + y[1] if self.total else y[0]

    def deriv(self, t):
        """Time derivative from the right-hand side of (1.5.1-2) to (1.5.1-4)."""
        kf, kr, kc = self.p
        y = self.sol(t)
        if self.total:
            return -kc*y[1]
        return -kf*y[0]*(self.CE0 - y[1]) + kr*y[1]


def mass_action(p, CE0, te, a0, method="LSODA", rtol=1e-11, atol=1e-14,
                dense=False, total=False):
    """Integrate the scheme (1.5.1-1) with no approximation.

    Parameters
    ----------
    p : tuple of float
        (k_1, k_-1, k_2).
    CE0 : float
        Total enzyme, cmol/L.
    te : array_like
        Output times (min); ``te[0]`` and ``te[-1]`` set the span.
    a0 : float
        Initial substrate, cmol/L. The complex starts at zero.
    method, rtol, atol
        ``solve_ivp`` settings; Radau is tried when the first method fails.
    dense : bool
        Return a :class:`DenseMassAction` instead of samples at ``te``.
    total : bool
        C_A + C_A-E instead of the free substrate C_A. Which of the two a reduced
        model is compared with is a convention; the page reports both where it
        changes a headline.

    Returns
    -------
    ndarray or DenseMassAction

    Raises
    ------
    RuntimeError
        When neither integrator succeeds.
    """
    kf, kr, kc = p

    def f(_, y):
        v1 = kf*y[0]*(CE0 - y[1])
        return [-v1 + kr*y[1], v1 - (kr + kc)*y[1]]

    def jac(_, y):
        return [[-kf*(CE0 - y[1]), kf*y[0] + kr],
                [kf*(CE0 - y[1]), -kf*y[0] - kr - kc]]

    kw = dict(method=method, jac=jac, rtol=rtol, atol=atol)
    span = (te[0], te[-1])
    s = solve_ivp(f, span, [a0, 0.0], t_eval=None if dense else te,
                  dense_output=dense, **kw)
    if s.status != 0:
        kw["method"] = "Radau"
        s = solve_ivp(f, span, [a0, 0.0], t_eval=None if dense else te,
                      dense_output=dense, **kw)
    if s.status != 0:
        raise RuntimeError(f"mass-action integration failed: {s.message}")
    if dense:
        return DenseMassAction(s.sol, p, CE0, total)
    return s.y[0] + s.y[1] if total else s.y[0]


def max_gap(f, g, window=(0.0, 18.0), n_scan=181, with_t=False, interp_tol=1e-9):
    """Largest |f(t) - g(t)| over a window, at a root of the analytic derivative.

    The difference is scanned on a uniform grid of ``n_scan`` points joined with
    a logarithmic one when the window starts at zero (the induction layer sits
    near t = 1e-7 min and a uniform scan misses it). Inside the two scan
    intervals around the largest sampled value, the maximum is the single root
    of ``f.deriv - g.deriv`` (the largest of several where the gap is flat). A
    maximum on a window end is that end.

    Parameters
    ----------
    f, g : callable
        Curves with a ``deriv`` method (:class:`DenseMassAction`,
        :class:`ReducedCurve`).
    window : tuple of float
        (t_start, t_end), min.
    n_scan : int
        Uniform scan points.
    with_t : bool
        Return (value, t) instead of the value.
    interp_tol : float
        Absolute amount (cmol/L) by which the value at the root may fall short
        of the scan: the derivative comes from the right-hand side at the
        interpolated state, which matches the dense interpolant's own slope
        only to the interpolation accuracy (measured 3.6e-10 cmol/L at
        rtol 1e-11; the page passes 1e3 rtol c_A0, at least 1e-8).

    Raises
    ------
    RuntimeError
        When the derivative has no root next to an interior maximum, or
        the root gives less than the scan.
    """
    a, b = window
    ts = np.unique(np.concatenate([
        np.linspace(a, b, n_scan),
        np.geomspace(max(a, 1e-9), b, 61) if a < 1e-6 else np.array([])]))
    d = np.abs(np.asarray(f(ts), float) - np.asarray(g(ts), float))
    i = int(np.argmax(d))
    if i in (0, len(ts) - 1):                 # the maximum is on a window end
        val, tpk = float(d[i]), float(ts[i])
    else:
        lo, hi = float(ts[i - 1]), float(ts[i + 1])
        r = find_roots(lambda t: float(f.deriv(t) - g.deriv(t)), lo, hi, n_scan=9,
                       xtol=1e-10*(hi - lo))
        if not r["roots"] or r["failed"]:
            raise RuntimeError(f"max_gap: no root of the derivative "
                               f"on [{lo:.3e}, {hi:.3e}]: {r['message']}")
        # several stationary points in the bracket (a flat gap): the maximum is
        # the largest of them
        vals = [float(abs(f(t) - g(t))) for t in r["roots"]]
        tpk = float(r["roots"][int(np.argmax(vals))])
        val = max(vals)
        if val < d[i] - interp_tol:
            raise RuntimeError(f"max_gap: root at t = {tpk:.3e} gives {val:.6e}, "
                               f"below the scan's {d[i]:.6e}")
    return (val, tpk) if with_t else val


class BatchMassAction:
    """The scheme (1.5.1-1) in a batch: one cell, four fields, theta time stepping.

    Parameters
    ----------
    p : tuple of float
        (k_1, k_-1, k_2).
    theta : float
        1 is backward Euler (the page's scheme), 0.5 Crank-Nicolson.
    tol : float
        Absolute ``newton`` step tolerance, cmol/L.
    backward_tol : float or None
        ``newton`` backward-error test; None stops on the step alone.
    p_row_factor : float
        Multiplies the product row of the source; 1 is the scheme.
    reverse : bool
        Keep the k_-1 term; False deletes it.

    Notes
    -----
    Fields (C_A, C_E, C_A-E, C_P) on shape (1, 4): ``NumJac((1, 4))`` couples the
    last axis in full, the pointwise block of a reaction; a bare ``(4,)`` would
    make space the last axis.
    """
    NC = 4

    def __init__(self, p, theta=1.0, tol=1e-13, p_row_factor=1.0, reverse=True,
                 backward_tol=1e-8):
        self.p, self.theta, self.tol, self.backward_tol = p, theta, tol, backward_tol
        self.p_row_factor, self.reverse = p_row_factor, reverse
        self.numjac = NumJac((1, self.NC))
        self.eye = speye(self.NC, format="csc")

    def source(self, c):
        """Net production rates of (C_A, C_E, C_A-E, C_P), eqs. (1.5.1-2) to (1.5.1-4)."""
        kf, kr, kc = self.p
        cA, cE, cAE = c[..., 0], c[..., 1], c[..., 2]
        v1, vm, v2 = kf*cA*cE, (kr*cAE if self.reverse else 0.0*cAE), kc*cAE
        return np.stack([-v1 + vm, -v1 + vm + v2, v1 - vm - v2,
                         self.p_row_factor*v2], axis=-1)

    def step(self, c_old, dt):
        """One theta step; raises when ``newton`` does not converge."""
        th = self.theta
        f_old = self.source(c_old)

        def res(w):
            c = w.reshape(1, self.NC)
            f, jf = self.numjac(self.source, c)
            return (((c - c_old)/dt - th*f - (1 - th)*f_old).reshape(-1, 1),
                    (self.eye/dt - th*jf).tocsc())
        r = newton(res, c_old.reshape(-1, 1).copy(), tol=self.tol, maxfev=50,
                   stagnation=True, backward_tol=self.backward_tol)
        if not r.success:
            raise RuntimeError(f"BatchMassAction step: newton failed: {r.message}")
        return r.x.reshape(1, self.NC)

    def run(self, c0, t_end, n_steps, n_out):
        """March to ``t_end`` in ``n_steps`` steps; return ``n_out + 1`` equidistant states."""
        dt, every = t_end/n_steps, n_steps//n_out
        c = np.asarray(c0, float).reshape(1, self.NC)
        out = np.empty((n_out + 1, self.NC))
        out[0], j = c[0], 1
        for i in range(1, n_steps + 1):
            c = self.step(c, dt)
            if i % every == 0:
                out[j] = c[0]
                j += 1
        return out


class BatchReduced:
    """dC_A/dt = -r_m C_A/(K_M + C_A), the rate law (1.5.1-15), by backward Euler.

    Parameters
    ----------
    rm, KM : float
        r_m (cmol/(L min)) and K_M (cmol/L).

    Notes
    -----
    One field on shape (1, 1), so ``NumJac((1, 1))``, never ``NumJac((1,))``.
    """

    def __init__(self, rm, KM):
        self.rm, self.KM = rm, KM
        self.numjac = NumJac((1, 1))
        self.eye = speye(1, format="csc")

    def source(self, c):
        """The rate law (1.5.1-15) as a sink of C_A."""
        return -self.rm*c/(self.KM + c)

    def run(self, c0, t_end, n_steps, t_out):
        """March to ``t_end``; return C_A at the step times nearest ``t_out``.

        Raises
        ------
        RuntimeError
            When ``newton`` does not converge in a step.
        """
        dt = t_end/n_steps
        c = np.array([[c0]], float)
        idx = np.round(np.asarray(t_out)/dt).astype(int)
        out, k = [], 0
        for i in range(n_steps + 1):
            while k < len(idx) and idx[k] == i:
                out.append(c[0, 0])
                k += 1
            if i == n_steps:
                break
            c_old = c

            def res(w, c_old=c_old):
                cc = w.reshape(1, 1)
                f, jf = self.numjac(self.source, cc)
                return ((cc - c_old)/dt - f).reshape(-1, 1), (self.eye/dt - jf).tocsc()
            r = newton(res, c.reshape(-1, 1).copy(), tol=1e-14, maxfev=50,
                       stagnation=True)
            if not r.success:
                raise RuntimeError(f"BatchReduced step: newton failed: {r.message}")
            c = r.x.reshape(1, 1)
        return np.array(out)
