"""A two-detector tracer experiment on a packed bed, and a two-zone bed.

Model code of gallery page A2.5 (Edwards & Richardson 1968). The page imports
these; the correlations (Eqs. 15 to 18), the moment formulae applied to the
digitised data and every fit stay in the notebook.

Units are those of the paper: cm, s, cm2/s. ``u`` is the interstitial velocity.

Contents
--------
TracerBed
    ``dc/dt + u dc/dx = D d2c/dx2`` on a uniform grid, van Leer TVD convection
    and a theta-method in time. Each step is solved by deferred correction on
    one factorisation, or by ``newton`` with the exact TVD Jacobian
    (``ConvFluxTVD``) as a second route.
TwoZoneBed
    Two mobile zones that exchange laterally, layout ``(n, 2)``, the exchange
    assembled by ``construct_coefficient_matrix``.
simulate_experiment, simulate_two_zone
    One pulse, two detectors: the response curves a two-point experiment records.
moments, invert_two_point
    Edwards & Richardson's Eqs. (2) to (4) and (6), (7).
emg, inverse_gaussian_kernel, peak_time, D_two_zone
    The injected pulse, the impulse response of the dispersion equation, a
    refined peak time and the two-zone closed form.
channelling_case, two_zone_split
    A two-zone bed tuned to an apparent ``D_L``, and the shape diagnostic of
    the second response curve computed from the first.
STEP_LOG
    Inner iterations of every time step taken, appended by both beds.

Notes
-----
The time step is the theta-method on the upwind and diffusive operator ``A``,
with the TVD limiter correction (the TVD face flux minus the upwind one) taken
at the new time level:

    (M/dt + theta A) c_new + div(F_tvd(c_new) - F_up(c_new))
        = M c_old/dt - (1 - theta) A c_old - b0.

Both solvers solve exactly these equations. Deferred correction iterates the
correction on a factorisation of ``M/dt + theta A`` until the change of a sweep
is below ``tol``; ``newton`` uses the exact Jacobian and refactorises at every
iteration.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu
from scipy.special import erfc

from pymrm import (ConvFluxTVD, construct_coefficient_matrix, construct_convflux_upwind,
                   construct_div, construct_grad, interp_cntr_to_stagg_tvd, newton,
                   vanleer)

__all__ = ["TracerBed", "TwoZoneBed", "simulate_experiment", "simulate_two_zone",
           "moments", "invert_two_point", "emg", "inverse_gaussian_kernel",
           "peak_time", "D_two_zone", "channelling_case", "two_zone_split", "STEP_LOG",
           "X0", "X1", "SIGMA", "TAIL", "L_TAIL"]

# The bed layout, in cm: pulse released at X0, first detector at X1, and the bed
# continues L_TAIL past the second detector. SIGMA and TAIL shape the pulse.
X0, X1, SIGMA, TAIL, L_TAIL = 30.0, 60.0, 4.0, 8.0, 80.0

STEP_LOG: list[tuple[str, int]] = []
"""``(method, iterations)`` for every time step taken by either bed."""

# inlet, outward normal n = -x:  c = 0 (no tracer enters)  ->  a = 0, b = 1, d = 0
BC_INLET = {"a": 0.0, "b": 1.0, "d": 0.0}
# outlet: the face carries the last cell's value and dc/dx = 0 for diffusion
BC_OUTFLOW = {"outflow": True}
# the defect: an absorbing outlet, c = 0  ->  a = 0, b = 1, d = 0
BC_DIRICHLET = {"a": 0.0, "b": 1.0, "d": 0.0}


class _ThetaTVD:
    """The time step shared by both beds (see the module notes).

    Subclasses set ``shape``, ``x_f``, ``x_c``, ``bc``, ``v`` (face velocity,
    scalar or per field), ``div``, ``conv``, ``conv_bc``, ``A``, ``b0``, ``M``,
    ``limiter``, ``tol``, ``max_corr`` and ``method``, then call ``_init_step``.
    """

    def _init_step(self):
        self._key, self._lu, self._lhs = None, None, None
        self._last_corr = np.zeros(int(np.prod(self.shape)))
        self._conv_bc = np.asarray(self.conv_bc.todense()).ravel()
        self._tvd = (ConvFluxTVD(self.shape, self.x_f, self.x_c, self.bc,
                                 tvd_limiter=self.limiter)
                     if self.limiter is not None and self.method == "newton" else None)

    def _factor(self, dt, theta):
        if self._key != (dt, theta):
            self._lhs = (self.M / dt + theta * self.A).tocsc()
            self._lu = splu(self._lhs)
            self._key = (dt, theta)
        return self._lu

    def correction(self, c):
        """``div(F_tvd(c) - F_up(c))``, the limiter correction to the upwind flux."""
        _, dc_f = interp_cntr_to_stagg_tvd(c.reshape(self.shape), self.x_f, self.x_c,
                                           self.bc, self.v, tvd_limiter=self.limiter,
                                           axis=0)
        return np.asarray(self.div @ (self.v * dc_f).reshape(-1, 1)).ravel()

    def step(self, c, dt, theta=1.0):
        """One theta-method step from the flattened state ``c``.

        Raises
        ------
        RuntimeError
            When deferred correction reaches ``max_corr`` sweeps without meeting
            ``tol``, or ``newton`` does not succeed.
        """
        lu = self._factor(dt, theta)
        rhs = self.M @ c / dt - (1.0 - theta) * (self.A @ c) - self.b0
        if self.limiter is None:
            STEP_LOG.append(("upwind", 0))
            return lu.solve(rhs)
        if self.method == "newton":
            return self._step_newton(lu, rhs)
        c_new = lu.solve(rhs - self._last_corr)          # warm start
        for it in range(1, self.max_corr + 1):
            corr = self.correction(c_new)
            cc = lu.solve(rhs - corr)
            done = np.max(np.abs(cc - c_new)) <= self.tol * max(1.0, np.max(np.abs(cc)))
            c_new = cc
            if done:
                break
        else:
            raise RuntimeError(f"deferred correction did not converge in {self.max_corr} sweeps")
        STEP_LOG.append(("deferred", it))
        self._last_corr = corr
        return c_new

    def _step_newton(self, lu, rhs):
        """The same step by Newton with the exact TVD Jacobian.

        The residual is written as the correction to the upwind step ``c_up``
        (``lhs @ c_up = rhs``): ``lhs @ (c - c_up) + div(F_tvd(c) - F_up(c))``.
        No bounds: a Crank-Nicolson step may undershoot zero by round-off, and
        the residual is defined there.
        """
        c_up = lu.solve(rhs)
        lhs, div, conv, shape = self._lhs, self.div, self.conv, self.shape

        def residual(c):
            flux, jac_tvd = self._tvd(c.reshape(shape), self.v)
            f_up = conv @ c.reshape(-1) + self._conv_bc
            g = lhs @ (c.reshape(-1) - c_up) + div @ (flux.reshape(-1) - f_up)
            return g, (lhs + div @ (jac_tvd - conv)).tocsc()

        sol = newton(residual, c_up.reshape(shape), tol=self.tol, damping="monotone",
                     stagnation=True)
        if not sol.success:
            raise RuntimeError(f"TVD Newton step failed: {sol.message}")
        STEP_LOG.append(("newton", int(sol.nit)))
        return np.asarray(sol.x).ravel()


class TracerBed(_ThetaTVD):
    """``dc/dt + u dc/dx = D d2c/dx2`` on ``0 <= x <= L``, uniform grid.

    The pulse is released in the interior and must never reach either end, so
    the boundary conditions only have to be harmless: no tracer enters at
    ``x = 0``, and the outlet is an outflow face.

    Parameters
    ----------
    L : float
        Bed length, cm.
    n : int
        Number of cells.
    u : float
        Interstitial velocity, cm/s.
    D : float
        Axial dispersion coefficient, cm2/s.
    limiter : callable or None, optional
        TVD limiter (van Leer); None gives first-order upwind.
    nu : int, optional
        Geometry of ``construct_div``: 0 Cartesian (the model), 1 cylindrical
        (a defect for the break table).
    outlet : {"outflow", "dirichlet"}, optional
        ``"dirichlet"`` makes the outlet absorbing, ``c = 0`` (a defect).
    tol : float, optional
        Deferred correction stops when a sweep changes no cell by more than
        ``tol * max(1, max|c|)``; it is also ``newton``'s step tolerance.
    max_corr : int, optional
        Cap on deferred-correction sweeps per step.
    method : {"deferred", "newton"}, optional
        How each TVD step is solved.
    """

    def __init__(self, L, n, u, D, limiter=vanleer, nu=0, outlet="outflow", tol=1e-9,
                 max_corr=40, method="deferred"):
        self.n, self.u, self.D, self.L, self.h = n, u, D, L, L / n
        self.shape = (n, 1)                       # (cells, fields), never (n,)
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.bc = (BC_INLET, BC_OUTFLOW if outlet == "outflow" else BC_DIRICHLET)
        self.v = u
        self.conv, self.conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                            self.bc, v=u)
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc)
        self.div = construct_div(self.shape, self.x_f, nu=nu)   # nu=0: Cartesian slab
        self.A = (self.div @ (self.conv - D * grad)).tocsc()
        self.b0 = np.asarray((self.div @ (self.conv_bc - D * grad_bc)).todense()).ravel()
        self.M = eye_array(n, format="csc")
        self.limiter, self.tol, self.max_corr, self.method = limiter, tol, max_corr, method
        self._init_step()


class TwoZoneBed(_ThetaTVD):
    """Two mobile zones exchanging laterally, layout ``(n, 2)``.

    ``th_i dc_i/dt + d(th_i u_i c_i)/dx = d(th_i D0 dc_i/dx)/dx -/+ kappa (c_1 - c_2)``

    Parameters
    ----------
    L, n : float, int
        Bed length (cm) and number of cells.
    th1 : float
        Volume fraction of zone 1; zone 2 holds ``1 - th1``.
    u1, u2 : float
        Interstitial velocities of the zones, cm/s.
    D0 : float
        Dispersion coefficient inside each zone, cm2/s.
    kappa : float
        Exchange coefficient per unit bed volume, 1/s.
    limiter, tol, max_corr, method
        As for :class:`TracerBed`.
    kappa_scale : float, optional
        Multiplies ``kappa`` in the assembly only (a defect for the break table).
    mass : {"theta", "identity"}, optional
        ``"identity"`` drops the volume fractions from the accumulation (a defect).
    """

    def __init__(self, L, n, th1, u1, u2, D0, kappa, limiter=vanleer, tol=1e-9,
                 max_corr=40, kappa_scale=1.0, mass="theta", method="deferred"):
        kappa = kappa * kappa_scale
        th = np.array([th1, 1.0 - th1])
        uu = np.array([u1, u2])
        self.n, self.th, self.uu, self.h = n, th, uu, L / n
        self.shape = (n, 2)                        # spatial axis first, fields last
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.bc = (BC_INLET, BC_OUTFLOW)           # c = 0 in, outflow face out
        self.v = np.broadcast_to(th * uu, (n + 1, 2))
        self.conv, self.conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                            self.bc, v=self.v, axis=0)
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc, axis=0)
        self.div = construct_div(self.shape, self.x_f, nu=0, axis=0)   # nu=0: Cartesian
        # face dispersion th_i D0, constant per zone, so no face mean is needed
        d_face = construct_coefficient_matrix(th * D0, shape=self.shape, axis=0)
        # per-cell exchange: +kappa (c1 - c2) in zone 1, -kappa (c1 - c2) in zone 2
        exchange = construct_coefficient_matrix(kappa * np.array([[1.0, -1.0], [-1.0, 1.0]]),
                                                shape=((n, 2, 1), (n, 1, 2)))
        self.A = (self.div @ (self.conv - d_face @ grad) + exchange).tocsc()
        self.b0 = np.asarray((self.div @ (self.conv_bc - d_face @ grad_bc)).todense()).ravel()
        # accumulation carries th_i: the volume fraction multiplies dc_i/dt
        self.M = (construct_coefficient_matrix(th, shape=self.shape) if mass == "theta"
                  else eye_array(2 * n, format="csc"))
        self.limiter, self.tol, self.max_corr, self.method = limiter, tol, max_corr, method
        self._init_step()


def emg(x, x0, sigma, tail):
    """The injected pulse: a Gaussian with an exponential tail.

    Edwards & Richardson's Fig. 5 shows "a long tail" and their Eq. (8) fits it
    as an exponential decay. ``z = x0 - x > 0`` lies behind the pulse, so it
    arrives later.
    """
    z = x0 - x
    return (0.5 / tail) * np.exp(0.5 * (sigma / tail) ** 2 - z / tail) * \
        erfc((sigma / tail - z / sigma) / np.sqrt(2))


def moments(t, c):
    """Area, mean and variance of a response curve, E&R Eqs. (2) to (4).

    Returns
    -------
    tuple of float
        ``(A, t_mean, variance)`` by the trapezoidal rule on the samples.
    """
    A = np.trapezoid(c, t)
    m1 = np.trapezoid(c * t, t) / A
    s2 = np.trapezoid(c * (t - m1) ** 2, t) / A
    return A, m1, s2


def invert_two_point(t, c1, c2, dx):
    """E&R Eqs. (6) and (7): velocity and dispersion from two response curves.

    Returns
    -------
    tuple of float
        ``(u_hat, D_hat, A1, A2)``.
    """
    A1, t1, s1 = moments(t, c1)
    A2, t2, s2 = moments(t, c2)
    u_hat = dx / (t2 - t1)
    D_hat = u_hat ** 3 * (s2 - s1) / (2.0 * dx)
    return u_hat, D_hat, A1, A2


def _record(t_end, dt):
    ns = int(np.ceil(t_end / dt))
    return ns, np.arange(ns + 1) * dt


def simulate_experiment(u, D, dx_probe=100.0, n=1200, cfl=0.5, theta=1.0, limiter=vanleer,
                        sigma=SIGMA, tail=TAIL, n_spread=7.0, l_tail=L_TAIL, t_end_mult=1.0,
                        **bed_kw):
    """Run E&R's experiment: one pulse, two detectors a distance ``dx_probe`` apart.

    The detectors read the cells whose centres are nearest ``X1`` and
    ``X1 + dx_probe``; the separation used in the inversion is that of those
    centres. ``t_end_mult`` scales how long the record runs and ``l_tail`` is
    how far the bed continues past the second detector; both are defects
    waiting to be injected. ``bed_kw`` goes to :class:`TracerBed`.

    Returns
    -------
    tuple
        ``(bed, t, c1, c2, dx)``.
    """
    x2 = X1 + dx_probe
    bed = TracerBed(x2 + l_tail, n, u, D, limiter=limiter, **bed_kw)
    c = emg(bed.x_c, X0, sigma, tail)
    i1 = int(np.argmin(np.abs(bed.x_c - X1)))
    i2 = int(np.argmin(np.abs(bed.x_c - x2)))
    dt = cfl * bed.h / u
    spread = np.sqrt(2.0 * D * (bed.x_c[i2] - X0) / u)
    t_end = t_end_mult * (bed.x_c[i2] - X0 + 6.0 * (sigma + tail) + n_spread * spread) / u
    ns, t = _record(t_end, dt)
    c1, c2 = np.empty(ns + 1), np.empty(ns + 1)
    c1[0], c2[0] = c[i1], c[i2]
    for k in range(1, ns + 1):
        c = bed.step(c, dt, theta=theta)
        c1[k], c2[k] = c[i1], c[i2]
    return bed, t, c1, c2, bed.x_c[i2] - bed.x_c[i1]


def D_two_zone(th1, u1, u2, kappa):
    """Added dispersion of the two-zone model, ``th1^2 th2^2 (u1 - u2)^2 / kappa``.

    From expanding ``c_i = c + a_i`` with ``th_1 a_1 + th_2 a_2 = 0``; it shares
    no code with :class:`TwoZoneBed`.
    """
    th2 = 1.0 - th1
    return (th1 * th2) ** 2 * (u1 - u2) ** 2 / kappa


def simulate_two_zone(th1, u1, u2, D0, kappa, dx_probe=100.0, n=1200, cfl=0.5, sigma=SIGMA,
                      tail=TAIL, theta=0.5, detector="flux", **bed_kw):
    """The two-detector experiment on a :class:`TwoZoneBed`.

    Both zones start with the same pulse. Each detector records the flux
    average ``sum(th_i u_i c_i) / u`` (the flowing concentration), or with
    ``detector="volume"`` the void average ``sum(th_i c_i)`` (a defect for the
    break table). ``bed_kw`` goes to :class:`TwoZoneBed`.

    Returns
    -------
    tuple
        ``(bed, t, c1, c2, dx, D_closed_form)``.
    """
    u = th1 * u1 + (1.0 - th1) * u2
    x2 = X1 + dx_probe
    bed = TwoZoneBed(x2 + L_TAIL, n, th1, u1, u2, D0, kappa, **bed_kw)
    c0 = emg(bed.x_c, X0, sigma, tail)
    c = np.column_stack([c0, c0]).ravel()
    i1 = int(np.argmin(np.abs(bed.x_c - X1)))
    i2 = int(np.argmin(np.abs(bed.x_c - x2)))
    dt = cfl * bed.h / max(u1, u2)
    Dg = D0 + D_two_zone(th1, u1, u2, kappa)
    spread = np.sqrt(2.0 * Dg * (bed.x_c[i2] - X0) / u)
    t_end = (bed.x_c[i2] - X0 + 6.0 * (sigma + tail) + 8.0 * spread) / u
    ns, t = _record(t_end, dt)
    # the detector sees the flux average (or, as a defect, the void average)
    w = bed.th * bed.uu / u if detector == "flux" else bed.th
    c1, c2 = np.empty(ns + 1), np.empty(ns + 1)
    cr = c.reshape(n, 2)
    c1[0], c2[0] = cr[i1] @ w, cr[i2] @ w
    for k in range(1, ns + 1):
        c = bed.step(c, dt, theta=theta)
        cr = c.reshape(n, 2)
        c1[k], c2[k] = cr[i1] @ w, cr[i2] @ w
    return bed, t, c1, c2, bed.x_c[i2] - bed.x_c[i1], Dg


def inverse_gaussian_kernel(t, dx, u, D):
    """Impulse response of the dispersion equation over ``dx`` in an infinite bed.

    E&R's Appendix B calculation in closed form: area 1, mean ``dx/u`` and
    variance ``2 D dx/u^3``, exactly their Eqs. (5) to (7).
    """
    g = np.zeros_like(t)
    m = t > 0
    g[m] = dx / np.sqrt(4 * np.pi * D * t[m] ** 3) * np.exp(-(dx - u * t[m]) ** 2
                                                             / (4 * D * t[m]))
    return g


def peak_time(t, c):
    """Time of the maximum of a sampled curve, refined by a parabola through the
    largest sample and its two neighbours."""
    i = int(np.argmax(c))
    y0, y1, y2 = c[i - 1], c[i], c[i + 1]
    return t[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * (t[1] - t[0])


def channelling_case(flux_frac, u, D0, D_apparent, th1=0.15, n=1000, dx_probe=100.0,
                     sigma=SIGMA, tail=TAIL, detector="flux", **bed_kw):
    """A two-zone bed whose two-point moments return a given apparent ``D_L``.

    Zone 1 holds the volume fraction ``th1`` and the share ``flux_frac`` of the
    total flow; ``kappa`` is set from :func:`D_two_zone` so that
    ``D0 + D_add = D_apparent``. The second response curve is then computed
    from the first with the fitted ``u_hat`` and ``D_hat``, as Edwards &
    Richardson's Appendix B does, and compared with the simulated one.

    Parameters
    ----------
    flux_frac : float
        Share of the total flow carried by zone 1.
    u : float
        Mean interstitial velocity, cm/s.
    D0 : float
        Dispersion coefficient inside each zone, cm2/s.
    D_apparent : float
        The apparent dispersion coefficient the bed must reproduce, cm2/s.
    th1, n, dx_probe, sigma, tail, detector
        As for :func:`simulate_two_zone`.
    **bed_kw
        Passed to :class:`TwoZoneBed` (``kappa_scale`` and ``mass`` are defects).

    Returns
    -------
    dict
        ``u1``, ``u2``, ``kappa``, ``D_hat``, ``lag`` (computed peak time minus
        simulated, over the simulated), ``resid`` (peak-to-peak of the scaled
        difference of the two curves) and the curves ``t``, ``sim``, ``comp``,
        ``first``.
    """
    u1 = flux_frac * u / th1
    u2 = (1.0 - flux_frac) * u / (1.0 - th1)
    kappa = D_two_zone(th1, u1, u2, 1.0) / (D_apparent - D0)
    _, t, a, b, dxx, _ = simulate_two_zone(th1, u1, u2, D0, kappa, n=n, sigma=sigma,
                                           tail=tail, dx_probe=dx_probe, detector=detector,
                                           **bed_kw)
    uh, Dh, _, _ = invert_two_point(t, a, b, dxx)
    g = inverse_gaussian_kernel(t, dxx, uh, Dh)
    g = g / np.trapezoid(g, t)
    comp = np.convolve(a, g)[:len(t)] * (t[1] - t[0])
    p_sim, p_comp = peak_time(t, b), peak_time(t, comp)
    resid = float(np.ptp(comp / comp.max() - b / b.max()))
    return dict(flux_frac=flux_frac, th1=th1, u1=u1, u2=u2, kappa=kappa, D_hat=Dh,
                lag=(p_comp - p_sim) / p_sim, resid=resid, dx=dx_probe, t=t, sim=b,
                comp=comp, first=a)


def two_zone_split(th1, ratio, u, D0, D_apparent, n=500, **kw):
    """:func:`channelling_case` at volume fraction ``th1`` and velocity ratio
    ``u1/u2 = ratio``, still reproducing ``D_apparent``."""
    u2 = u / (th1 * ratio + (1.0 - th1))
    return channelling_case(th1 * ratio * u2 / u, u, D0, D_apparent, th1=th1, n=n, **kw)
