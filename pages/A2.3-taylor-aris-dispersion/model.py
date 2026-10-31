"""Taylor-Aris dispersion in a tube, assembled from pymrm operators.

Model code of gallery page A2.3 (Taylor, *Proc. R. Soc. Lond. A* 219, 1953,
186-203). Everything is dimensionless: lengths in tube radii ``a``, time in
``a^2/D``, velocity in ``D/a``, so the only parameter is the Peclet number
``Pe = u0 a / D`` with ``u0`` the CENTRELINE velocity. Taylor's Eq. 25 then
reads ``k/D = Pe^2/192``; the page keeps that closed form, and his worked
capillary run, in the notebook, so that the solves and the paper's algebra
never share a line of code.

Contents
--------
radial_grid
    Uniform cells on ``0 <= r <= 1`` and their area weights.
velocity_deviation
    Poiseuille velocity and its deviation from the section mean.
radial_eigenvalue
    The slowest decay rate of the discrete radial operator.
closure
    The homogenisation closure: a 1-D steady boundary-value problem in ``r``
    whose solution gives ``k/D``. Linear, solved with one sparse solve.
Slug
    The full 2-D axisymmetric transient: a slug released across the section,
    in a frame moving with the mean velocity, by backward Euler on one
    factorisation of a constant operator.
History
    The section-mean moments of a :class:`Slug` run at every time step, with
    the instantaneous variance growth rate and the settling times derived
    from it by root finding.
run_history
    Build a :class:`Slug`, march it, and return its :class:`History`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu, spsolve

from pymrm import construct_convflux_upwind, construct_div, construct_grad

__all__ = ["radial_grid", "velocity_deviation", "closure", "ClosureResult", "radial_eigenvalue",
           "Slug", "History", "run_history"]

# no flux through the axis (symmetry) and the wall: a dC/dn + b C = d with
# a = 1, b = 0, d = 0 at both ends; with b = d = 0 the outward normal's sign
# does not enter
NO_FLUX = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 0.0})


def radial_grid(n_r):
    """Uniform radial cells on ``[0, 1]``.

    Parameters
    ----------
    n_r : int
        Number of cells.

    Returns
    -------
    r_f, r_c, w : numpy.ndarray
        Faces ``(n_r + 1,)``, centres ``(n_r,)`` and area weights
        ``w_i = r_{f,i+1}^2 - r_{f,i}^2 = 2 r_{c,i} dr``, which sum to exactly 1:
        ``sum(w * f)`` is the cell-wise area average that the cylindrical
        divergence conserves.
    """
    r_f = np.linspace(0.0, 1.0, n_r + 1)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    return r_f, r_c, 2.0 * r_c * np.diff(r_f)


def velocity_deviation(r_c, w, Pe, mean="discrete"):
    """Poiseuille velocity and its deviation from the section mean.

    Parameters
    ----------
    r_c, w : numpy.ndarray
        Cell centres and area weights from :func:`radial_grid`.
    Pe : float
        ``u0 a / D``, with ``u0`` the centreline velocity.
    mean : {"discrete", "continuous"}, optional
        ``"discrete"`` (default) subtracts the area average on the grid,
        ``sum(w * u)``, so that ``sum(w * u')`` is zero to round-off and the
        pure-Neumann closure is consistent. ``"continuous"`` subtracts the exact
        mean ``Pe/2``, whose grid average differs from zero by a second-order
        amount; it exists for the page's comparison of the two.

    Returns
    -------
    u, u_dev : numpy.ndarray
        Velocity ``Pe (1 - r^2)`` and its deviation, shape ``(n_r,)``.
    """
    u = Pe * (1.0 - r_c ** 2)
    if mean == "discrete":
        return u, u - np.sum(w * u)
    if mean == "continuous":
        return u, u - 0.5 * Pe       # the mean of Poiseuille is half the centreline
    raise ValueError(f"mean must be 'discrete' or 'continuous', not {mean!r}")


@dataclass
class ClosureResult:
    """Solution of the homogenisation closure.

    Attributes
    ----------
    k_over_D : float
        ``k/D = -<u' B>``, Taylor's dispersion coefficient over ``D``.
    B : numpy.ndarray
        Closure field, shape ``(n_r, 1)``.
    r_c, w, u_dev : numpy.ndarray
        Cell centres, area weights and the velocity deviation used.
    """

    k_over_D: float
    B: np.ndarray
    r_c: np.ndarray
    w: np.ndarray
    u_dev: np.ndarray


def closure(Pe, n_r=200, nu=1, mean="discrete", subtract_mean=True, pin=0):
    """``k/D`` from the homogenisation closure.

    ``B(r)`` solves ``(1/r) d/dr (r dB/dr) = u'(r)`` with ``dB/dr = 0`` at the
    axis and the wall, and ``k/D = -<u' B>``. The system is pure Neumann, so it
    is singular by one constant; one cell's equation is replaced by
    ``B_pin = 0``. Because the cylindrical divergence conserves ``sum(w * .)``,
    the dropped equation is implied by the others exactly when
    ``sum(w * u') = 0``, which the discrete mean (default) guarantees; then the
    pinned cell cannot affect ``k``.

    Parameters
    ----------
    Pe : float
        Peclet number ``u0 a / D``.
    n_r : int, optional
        Number of radial cells.
    nu : int, optional
        Geometry of ``construct_div``: 1 cylindrical (the tube). 0, the plane
        channel, exists for the page's break table.
    mean : {"discrete", "continuous"}, optional
        Passed to :func:`velocity_deviation`.
    subtract_mean : bool, optional
        False puts the full velocity ``u`` on the right-hand side instead of
        ``u'``, an inconsistent system; it exists for the break table.
    pin : int, optional
        Index of the cell whose equation is replaced by ``B = 0``.

    Returns
    -------
    ClosureResult
    """
    r_f, r_c, w = radial_grid(n_r)
    shape = (n_r, 1)
    grad, grad_bc = construct_grad(shape, r_f, r_c, NO_FLUX)
    div = construct_div(shape, r_f, nu=nu)          # nu = 1: cylindrical radial
    lap = (div @ grad).tolil()
    rhs = -(div @ grad_bc).toarray().ravel()
    u, u_dev = velocity_deviation(r_c, w, Pe, mean)
    rhs = rhs + (u_dev if subtract_mean else u)
    lap[pin, :] = 0.0
    lap[pin, pin] = 1.0
    rhs[pin] = 0.0
    B = spsolve(lap.tocsc(), rhs)
    return ClosureResult(-float(np.sum(w * u_dev * B)), B.reshape(shape), r_c, w, u_dev)


def radial_eigenvalue(n_r, nu=1):
    """Smallest non-zero eigenvalue of the discrete no-flux radial operator.

    The operator is ``-(1/r) d/dr (r d/dr)`` on ``n_r`` cells of
    :func:`radial_grid`, assembled as in :class:`Slug`; its continuous
    counterpart is ``j_1^2 = 3.8317^2``.

    Parameters
    ----------
    n_r : int
        Number of radial cells.
    nu : int, optional
        Geometry of ``construct_div``.

    Returns
    -------
    float
    """
    r_f, r_c, _ = radial_grid(n_r)
    shape = (n_r, 1)
    grad, _ = construct_grad(shape, r_f, r_c, NO_FLUX)
    div = construct_div(shape, r_f, nu=nu)
    lam = np.sort(np.linalg.eigvals(-(div @ grad).toarray()).real)
    return float(lam[1])


class Slug:
    """Transient 2-D axisymmetric release, in a frame moving at the mean speed.

    ``dC/dt + u'(r) dC/dz = d2C/dz2 + (1/r) d/dr (r dC/dr)`` on
    ``-L <= z <= L``, ``0 <= r <= 1``, with first-order upwind convection along
    ``z`` (``construct_convflux_upwind``), state shape ``(n_z, n_r, 1)``
    (spatial axes first, one field). The constant operator is assembled once in
    ``__init__`` and factorised once per time step size.

    Parameters
    ----------
    Pe : float
        Peclet number ``u0 a / D``.
    half_len : float, optional
        Half length ``L`` of the axial domain, in tube radii. The slug never
        reaches the ends; the ends are closed (no convective and no diffusive
        flux), so the operator conserves mass exactly whatever the state.
    n_z, n_r : int, optional
        Axial and radial cell counts.
    width : float, optional
        Axial width of the initial slug, centred on ``z = 0``.
    release_radius : float, optional
        The slug fills ``r < release_radius`` at ``t = 0``; 1 (default) is a
        release across the whole section.
    nu_r : int, optional
        Geometry of the radial divergence: 1 cylindrical (the tube). 0 exists
        for the page's break table.
    wall_b : float, optional
        ``b`` of the wall condition ``dC/dn + b C = 0``: 0 (default) is the
        impermeable wall; ``b > 0`` is a first-order wall uptake, which Taylor's
        model excludes; it exists for the break table.
    mean : {"discrete", "continuous"}, optional
        Passed to :func:`velocity_deviation`.
    closed_ends : bool, optional
        True (default) sets the velocity to zero on the two end faces, so no
        solute crosses them. False leaves upwind carrying ``u' C`` through
        them; it exists for the page's break table.

    Attributes
    ----------
    jac : scipy.sparse.csc_array
        Constant operator ``J``, so that ``dC/dt = -(J C + g)``.
    g : numpy.ndarray
        Constant boundary contribution (zero for the default boundaries).
    c : numpy.ndarray
        Current state, flattened.
    t : float
        Current time.
    """

    def __init__(self, Pe, half_len=40.0, n_z=1600, n_r=24, width=1.0,
                 release_radius=1.0, nu_r=1, wall_b=0.0, mean="discrete", closed_ends=True):
        self.Pe, self.shape = Pe, (n_z, n_r, 1)
        self.z_f = np.linspace(-half_len, half_len, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.r_f, self.r_c, self.w = radial_grid(n_r)
        # axial ends: no diffusive flux, dC/dn = 0, and (below) no convective
        # flux, so the total flux u' C - dC/dz through both end faces is zero
        bc_z = NO_FLUX
        # axis: dC/dn = 0 by symmetry; wall: dC/dn + b C = 0 with n = +r, i.e.
        # -dC/dr = b C, an uptake flux b C (b = 0: the impermeable wall)
        bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": wall_b, "d": 0.0})
        _, u_dev = velocity_deviation(self.r_c, self.w, Pe, mean)
        v = np.broadcast_to(u_dev[None, :, None], (n_z + 1, n_r, 1)).copy()   # face velocities
        # closed ends: u' = 0 on the two end faces, u' C = 0 at z = -L and z = +L
        if closed_ends:
            v[0] = 0.0
            v[-1] = 0.0

        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                  bc_z, v=v, axis=0)
        gz, gz_bc = construct_grad(self.shape, self.z_f, self.z_c, bc_z, axis=0)
        div_z = construct_div(self.shape, self.z_f, nu=0, axis=0)        # nu = 0: Cartesian
        gr, gr_bc = construct_grad(self.shape, self.r_f, self.r_c, bc_r, axis=1)
        div_r = construct_div(self.shape, self.r_f, nu=nu_r, axis=1)     # nu = 1: radial
        self.jac = (div_z @ (conv - gz) - div_r @ gr).tocsc()
        g = div_z @ (conv_bc - gz_bc) - div_r @ gr_bc
        self.g = np.asarray(g.todense() if hasattr(g, "todense") else g).ravel()

        dz = np.diff(self.z_f)
        vol = (dz[:, None] * self.w[None, :]).ravel()           # cell volume / (pi a^2)
        z = np.repeat(self.z_c, n_r)
        self._moments = np.vstack([vol, vol * z, vol * z ** 2])

        c = np.zeros(self.shape)
        inside = (np.abs(self.z_c)[:, None] < width / 2) & (self.r_c[None, :] < release_radius)
        c[inside] = 1.0
        c *= 1.0 / float(np.sum(self._moments[0] * c.ravel()))   # unit mass
        self.c, self.t, self._dt, self._lu = c.ravel(), 0.0, None, None

    def step(self, dt, n=1):
        """Advance ``n`` backward-Euler steps of size ``dt``."""
        if self._dt != dt:
            n_tot = self.c.size
            self._lu = splu((eye_array(n_tot, format="csc") / dt + self.jac).tocsc())
            self._dt = dt
        for _ in range(n):
            self.c = self._lu.solve(self.c / dt - self.g)
            self.t += dt

    def column_sums(self):
        """``V @ J``: zero in every column exactly when the operator conserves mass."""
        return self._moments[0] @ self.jac

    def rhs(self):
        """``dC/dt`` of the semi-discrete system at the current state."""
        return -(self.jac @ self.c + self.g)

    def moments(self):
        """Raw moments ``(M0, M1, M2)`` of the state and of ``dC/dt``.

        ``M_k = sum over cells of V C z^k``, the exact discrete integrals.
        After a backward-Euler step, ``dC/dt`` at the new state equals the
        step's difference quotient.
        """
        return self._moments @ self.c, self._moments @ self.rhs()

    def profile(self):
        """Section-mean concentration ``C_m(z) = sum(w C)``, shape ``(n_z,)``."""
        return self.c.reshape(self.shape)[..., 0] @ self.w


def _variance_rate(M, dM):
    """Half the time derivative of the axial variance, from raw moments."""
    M0, M1, M2 = M
    d0, d1, d2 = dM
    m1 = M1 / M0
    dvar = d2 / M0 - M2 * d0 / M0 ** 2 - 2.0 * m1 * (d1 / M0 - M1 * d0 / M0 ** 2)
    return 0.5 * dvar


@dataclass
class History:
    """Moments of a :class:`Slug` run at every backward-Euler step.

    Attributes
    ----------
    t : numpy.ndarray
        Times, shape ``(n + 1,)``, starting at 0.
    M, dM : numpy.ndarray
        Raw moments of the state and of its time derivative, ``(n + 1, 3)``.
    profiles : dict
        Section-mean profiles kept at requested times.
    slug : Slug
        The marched model, at its final state.
    """

    t: np.ndarray
    M: np.ndarray
    dM: np.ndarray
    profiles: dict = field(default_factory=dict)
    slug: Slug | None = None

    @property
    def mass(self):
        """Total mass ``M0`` at every step."""
        return self.M[:, 0]

    @property
    def variance(self):
        """Axial variance of the section-mean profile at every step."""
        M0, M1, M2 = self.M.T
        return M2 / M0 - (M1 / M0) ** 2

    @property
    def rate(self):
        """Instantaneous effective diffusivity ``(1/2) d(variance)/dt`` over ``D``.

        Computed from ``dC/dt`` of the semi-discrete system at each step, not
        by differencing sampled variances.
        """
        return _variance_rate(self.M.T, self.dM.T)

    def index(self, tau):
        """Step index of time ``tau`` (raises unless ``tau`` is a step time)."""
        dt = self.t[1] - self.t[0]
        k = int(round(tau / dt))
        if abs(self.t[k] - tau) > 1e-9 * max(1.0, tau):
            raise ValueError(f"tau = {tau} is not a step time")
        return k

    def window_dispersion(self, tau0, tau1):
        """``(var(tau1) - var(tau0)) / (2 (tau1 - tau0))``, the mean rate over a window."""
        k0, k1 = self.index(tau0), self.index(tau1)
        var = self.variance
        return (var[k1] - var[k0]) / (2.0 * (self.t[k1] - self.t[k0]))

    def rate_on_step(self, k, s):
        """The rate along the step from ``t[k]`` to ``t[k+1]``, at fraction ``s``.

        The state is taken linear in time across the step (backward Euler's own
        interpolant), so its moments and those of ``dC/dt`` are linear in ``s``.
        """
        M = self.M[k] + s * (self.M[k + 1] - self.M[k])
        dM = self.dM[k] + s * (self.dM[k + 1] - self.dM[k])
        return _variance_rate(M, dM)

    def fit_relaxation(self, target, floor, lo=0.2, hi=0.45):
        """Fit ``floor - (rate/target - 1) = A exp(-lam tau)`` on ``lo <= tau <= hi``.

        Parameters
        ----------
        target : float
            Reference rate (Taylor's ``1 + Pe^2/192``).
        floor : float
            The run's own late-time offset ``rate/target - 1``.
        lo, hi : float, optional
            Fit window, after the higher radial modes have decayed.

        Returns
        -------
        lam, A : float
            Decay rate and amplitude, from a least-squares line through the
            logarithm of the deficit.
        """
        sel = (self.t >= lo) & (self.t <= hi)
        deficit = floor - (self.rate[sel] / target - 1.0)
        if np.any(deficit <= 0.0):
            raise ValueError("the deficit is not positive over the fit window")
        slope, intercept = np.polyfit(self.t[sel], np.log(deficit), 1)
        return float(-slope), float(np.exp(intercept))

    def settle(self, target, band):
        """First time after which ``|rate/target - 1|`` stays below ``band``.

        The last step on which the deviation is at least ``band`` is found from
        the per-step rates, and the crossing inside the following step by a
        root on :meth:`rate_on_step`. Returns NaN when the deviation is still
        at least ``band`` at the end of the run.
        """
        dev = np.abs(self.rate / target - 1.0)
        bad = np.nonzero(dev >= band)[0]
        if len(bad) == 0:
            return float(self.t[0])
        k = int(bad[-1])
        if k == len(self.t) - 1:
            return float("nan")
        s = brentq(lambda s: abs(self.rate_on_step(k, s) / target - 1.0) - band,
                   0.0, 1.0, xtol=1e-13)
        return float(self.t[k] + s * (self.t[k + 1] - self.t[k]))


def run_history(Pe, n_z, dt=5e-4, tau_end=1.2, keep=(), **kw):
    """March a :class:`Slug` to ``tau_end`` and record its :class:`History`.

    Parameters
    ----------
    Pe : float
        Peclet number.
    n_z : int
        Axial cell count.
    dt : float, optional
        Backward-Euler step.
    tau_end : float, optional
        End time, a whole number of steps.
    keep : iterable of float, optional
        Step times at which to store the section-mean profile.
    **kw
        Passed to :class:`Slug`.

    Returns
    -------
    History
    """
    s = Slug(Pe, n_z=n_z, **kw)
    n = int(round(tau_end / dt))
    t = np.zeros(n + 1)
    M = np.empty((n + 1, 3))
    dM = np.empty((n + 1, 3))
    M[0], dM[0] = s.moments()
    keep_k = {int(round(tau / dt)): tau for tau in keep}
    profiles = {keep_k[0]: s.profile()} if 0 in keep_k else {}
    for k in range(1, n + 1):
        s.step(dt)
        t[k] = k * dt
        M[k], dM[k] = s.moments()
        if k in keep_k:
            profiles[keep_k[k]] = s.profile()
    s.t = t[-1]
    return History(t, M, dM, profiles, s)
