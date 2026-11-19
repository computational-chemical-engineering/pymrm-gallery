"""Duncan and Toor (1962): ternary gas diffusion in a two-bulb cell (page A4.9).

Duncan and Toor, AIChE J. 8, 38 (1962). Two well-mixed bulbs joined by a
capillary; the capillary holds a quasi-steady Maxwell-Stefan profile at every
instant and its end fluxes change the bulb compositions.

- the Maxwell-Stefan matrix ``[B]`` of the ``n - 1`` independent fractions
  (:func:`ms_b_matrix`) and the Wilke effective diffusivity (:func:`wilke_deff`);
- the cell itself (:class:`TwoBulbCell`): constant operators assembled once,
  the face fluxes, the steady capillary solve by ``newton`` on
  ``NumJac((n_z, n_c - 1))``, a backward-Euler (or theta) march of the bulbs
  and the same bulb equations integrated by ``solve_ivp`` with dense output;
- the two closed-form approximations the page compares with: Fick with Wilke
  effective diffusivities (:func:`fick_difference`) and Toor's linearisation
  with ``[B]`` frozen at the mean composition (:func:`toor_difference`).

Every defect the page's break table injects enters through an argument of
these functions or of the class constructor.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from pymrm import construct_grad, construct_div, NumJac, newton

NEWTON_TOL = 1e-12   # newton step tolerance; the unknowns are mole fractions, order one


def ms_b_matrix(x_ind, d_ms):
    """Maxwell-Stefan matrix ``[B]`` from the independent mole fractions.

    Parameters
    ----------
    x_ind : ndarray, shape (..., n_c - 1)
        Mole fractions of the first ``n_c - 1`` species; the last species
        closes the sum.
    d_ms : ndarray, shape (n_c, n_c)
        Symmetric Maxwell-Stefan pair diffusivities [m2/s].

    Returns
    -------
    ndarray, shape (..., n_c - 1, n_c - 1)
        ``B_ii = x_i/D_in + sum_{k != i} x_k/D_ik`` and
        ``B_ij = -x_i (1/D_ij - 1/D_in)``, with ``n`` the last species.
    """
    n_c = d_ms.shape[0]
    n_i = n_c - 1
    x = np.concatenate([x_ind, 1.0 - x_ind.sum(axis=-1, keepdims=True)], axis=-1)
    b = np.zeros(x_ind.shape[:-1] + (n_i, n_i))
    for i in range(n_i):
        acc = x[..., i] / d_ms[i, n_i]
        for k in range(n_c):
            if k != i:
                acc = acc + x[..., k] / d_ms[i, k]
        b[..., i, i] = acc
        for j in range(n_i):
            if j != i:
                b[..., i, j] = -x[..., i] * (1.0 / d_ms[i, j] - 1.0 / d_ms[i, n_i])
    return b


def wilke_deff(x, d_ms):
    """Wilke effective diffusivity of every species at composition ``x``.

    Parameters
    ----------
    x : ndarray, shape (n_c,)
        Mole fractions.
    d_ms : ndarray, shape (n_c, n_c)
        Maxwell-Stefan pair diffusivities [m2/s].

    Returns
    -------
    ndarray, shape (n_c,)
        ``(1 - x_i) / sum_{j != i} x_j / D_ij`` [m2/s].
    """
    n_c = len(x)
    return np.array([(1.0 - x[i]) / sum(x[j] / d_ms[i, j] for j in range(n_c) if j != i)
                     for i in range(n_c)])


class TwoBulbCell:
    """Two well-mixed bulbs joined by a quasi-steady Maxwell-Stefan capillary.

    Parameters
    ----------
    d_ms : ndarray, shape (n_c, n_c)
        Maxwell-Stefan pair diffusivities [m2/s].
    length : float
        Capillary length [m].
    la_eff : float
        Effective length-to-area ratio of the capillary [1/m]; the area is
        ``length / la_eff``.
    vol_1, vol_2 : float
        Bulb volumes [m3].
    c_t : float
        Total molar concentration [mol/m3].
    n_z : int, optional
        Capillary cells.
    nu : int, optional
        Geometry of the divergence: 0 Cartesian (the capillary), 1 cylindrical
        and 2 spherical exist only as injected defects.
    face : {"mean", "upwind"}, optional
        Composition at which ``[B]`` is evaluated on an interior face: the mean
        of the two cells (second order), or the left cell (an injected defect).
    tol : float, optional
        ``newton`` step tolerance of the capillary solve.

    Notes
    -----
    The bulb compositions are the boundary values of the capillary. They enter
    through the ``shapes_d`` matrices of ``construct_grad``, so the operators
    are assembled once although the boundary values change at every step.
    """

    def __init__(self, d_ms, length, la_eff, vol_1, vol_2, c_t, n_z=40, nu=0,
                 face="mean", tol=NEWTON_TOL):
        self.d_ms = np.asarray(d_ms, dtype=float)
        self.n_c = self.d_ms.shape[0]
        self.n_i = self.n_c - 1
        self.length, self.la_eff = length, la_eff
        self.area = length / la_eff
        self.vol_1, self.vol_2, self.c_t = vol_1, vol_2, c_t
        self.n_z, self.face, self.tol = n_z, face, tol
        self.shape = (n_z, self.n_i)                 # space first, species last
        self.z_f = np.linspace(0.0, length, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        # Outward normal, a dx/dn + b x = d. Both ends are Dirichlet on the bulb
        # composition: x(0) = x_bulb1 and x(L) = x_bulb2. With shapes_d, d = 1 is
        # the coefficient of the external boundary vector (the bulb fractions).
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0},        # z = 0: x = x_bulb1
              {"a": 0.0, "b": 1.0, "d": 1.0})        # z = L: x = x_bulb2
        self.grad, self.grad_bc_1, self.grad_bc_2 = construct_grad(
            self.shape, self.z_f, self.z_c, bc, axis=0,
            shapes_d=((1, self.n_i), (1, self.n_i)))
        self.div = construct_div(self.shape, self.z_f, nu=nu, axis=0)  # nu=0: Cartesian
        # The residual handed to NumJac is div(N(x)) with N built from [B] at the
        # face compositions, so it reads both neighbours: a [-1, 0, 1] stencil on
        # the spatial axis, species coupled in full at each point.
        self.numjac = NumJac(self.shape, axes_diagonals=[0], axes_blocks=[-1])
        self._x_warm = None

    @property
    def capillary_volume_fraction(self):
        """Capillary volume over the volume of bulb 1 [-]."""
        return self.area * self.length / self.vol_1

    def initial_profile(self, xb1, xb2):
        """Linear capillary profile between the two bulb compositions, shape (n_z, n_c - 1)."""
        s = (self.z_c / self.length)[:, None]
        return xb1[:self.n_i] + s * (xb2[:self.n_i] - xb1[:self.n_i])

    def flux(self, x, xb1, xb2):
        """Molar fluxes at all capillary faces [mol/(m2 s)], shape (n_z + 1, n_c - 1).

        Parameters
        ----------
        x : ndarray, shape (n_z, n_c - 1)
            Independent mole fractions in the capillary cells.
        xb1, xb2 : ndarray, shape (n_c,)
            Bulb compositions, the boundary values at ``z = 0`` and ``z = L``.
        """
        n_i, n_z = self.n_i, self.n_z
        grad_x = (self.grad @ x.reshape((-1, 1))
                  + self.grad_bc_1 @ xb1[:n_i].reshape((-1, 1))
                  + self.grad_bc_2 @ xb2[:n_i].reshape((-1, 1))).reshape(n_z + 1, n_i)
        x_face = np.empty((n_z + 1, n_i))
        x_face[1:-1] = 0.5 * (x[:-1] + x[1:]) if self.face == "mean" else x[:-1]
        x_face[0], x_face[-1] = xb1[:n_i], xb2[:n_i]
        return -self.c_t * np.linalg.solve(ms_b_matrix(x_face, self.d_ms),
                                           grad_x[..., None])[..., 0]

    def solve_capillary(self, xb1, xb2, x_init=None, maxfev=30, check=True):
        """Steady capillary profile, ``div N = 0``, for given bulb compositions.

        Parameters
        ----------
        xb1, xb2 : ndarray, shape (n_c,)
            Bulb compositions.
        x_init : ndarray, shape (n_z, n_c - 1), optional
            Initial guess; the linear profile when omitted.
        maxfev : int, optional
            Newton iterations allowed.
        check : bool, optional
            Raise when ``newton`` does not report success. ``False`` exists
            only for the injected unconverged solve of the break table.

        Returns
        -------
        ndarray, shape (n_z, n_c - 1)
        """
        if x_init is None:
            x_init = self.initial_profile(xb1, xb2)

        def residual(x):
            g, jac = self.numjac(
                lambda xx: (self.div @ self.flux(xx, xb1, xb2).reshape((-1, 1))).reshape(self.shape),
                x.reshape(self.shape))
            return g.reshape((-1, 1)), jac

        res = newton(residual, x_init.reshape((-1, 1)), maxfev=maxfev, tol=self.tol,
                     stagnation=True)
        if check and not res.success:
            raise RuntimeError(f"capillary newton failed: {res.message}")
        return res.x.reshape(self.shape)

    def flux_nonuniformity(self, xb1, xb2, maxfev=30, check=True):
        """Largest deviation of the face fluxes from their mean, over the largest flux [-]."""
        nf = self.flux(self.solve_capillary(xb1, xb2, maxfev=maxfev, check=check), xb1, xb2)
        return float(np.max(np.abs(nf - nf.mean(axis=0))) / np.max(np.abs(nf)))

    def end_fluxes(self, xb1, xb2):
        """Fluxes leaving bulb 1 (face ``z = 0``) and entering bulb 2 (face ``z = L``)."""
        self._x_warm = self.solve_capillary(xb1, xb2, x_init=self._x_warm)
        nf = self.flux(self._x_warm, xb1, xb2)
        return nf[0], nf[-1]

    def bulb_rates(self, xb1, xb2):
        """Time derivatives of the independent bulb fractions [1/s], two arrays (n_c - 1,)."""
        n0, nl = self.end_fluxes(xb1, xb2)
        return (-self.area * n0 / (self.vol_1 * self.c_t),
                self.area * nl / (self.vol_2 * self.c_t))

    def _full(self, x_ind):
        return np.append(x_ind, 1.0 - x_ind.sum())

    def march(self, x1_0, x2_0, t_end, n_steps, theta=1.0, coupling_tol=1e-14,
              max_coupling=30):
        """Theta march of the bulbs, backward Euler for ``theta = 1``.

        Each step solves ``xb_new = xb + dt (theta r(xb_new) + (1 - theta) r(xb))``
        by fixed-point iteration on the bulb compositions, a steady capillary
        solve per pass, until the largest change of a bulb fraction is below
        ``coupling_tol``.

        Parameters
        ----------
        x1_0, x2_0 : ndarray, shape (n_c,)
            Initial bulb compositions.
        t_end : float
            Duration [s].
        n_steps : int
            Number of equal steps.
        theta : float, optional
            Implicitness; 0.5 (trapezoidal) is an injected defect of the time
            study.
        coupling_tol, max_coupling : float, int, optional
            Convergence of the fixed-point iteration; raises when not met.

        Returns
        -------
        t : ndarray, shape (n_steps + 1,)
            Times [s].
        h1, h2 : ndarray, shape (n_steps + 1, n_c)
            Bulb compositions.
        passes : int
            Largest number of fixed-point passes in any step.
        """
        n_i = self.n_i
        dt = t_end / n_steps
        b1, b2 = np.asarray(x1_0, float)[:n_i].copy(), np.asarray(x2_0, float)[:n_i].copy()
        self._x_warm = None
        h1, h2, worst = [self._full(b1)], [self._full(b2)], 0
        p1, p2 = b1.copy(), b2.copy()                # last accepted step, for the predictor
        for _ in range(n_steps):
            r1_old, r2_old = (self.bulb_rates(self._full(b1), self._full(b2))
                              if theta < 1.0 else (0.0, 0.0))
            n1, n2 = 2.0 * b1 - p1, 2.0 * b2 - p2      # linear extrapolation as first guess
            p1, p2 = b1.copy(), b2.copy()
            for k in range(1, max_coupling + 1):
                r1, r2 = self.bulb_rates(self._full(n1), self._full(n2))
                m1 = b1 + dt * (theta * r1 + (1.0 - theta) * r1_old)
                m2 = b2 + dt * (theta * r2 + (1.0 - theta) * r2_old)
                change = max(np.max(np.abs(m1 - n1)), np.max(np.abs(m2 - n2)))
                n1, n2 = m1, m2
                if change < coupling_tol:
                    break
            else:
                raise RuntimeError(f"bulb coupling not converged in {max_coupling} passes "
                                   f"(last change {change:.1e})")
            worst = max(worst, k)
            b1, b2 = n1, n2
            h1.append(self._full(b1))
            h2.append(self._full(b2))
        return np.linspace(0.0, t_end, n_steps + 1), np.array(h1), np.array(h2), worst

    def integrate(self, x1_0, x2_0, t_end, rtol=1e-10, atol=1e-13):
        """The same bulb equations integrated by ``solve_ivp`` (DOP853) with dense output.

        Returns
        -------
        scipy OdeResult
            ``sol.sol(t)`` gives ``[x1_1..x1_{n-1}, x2_1..x2_{n-1}]`` at any
            ``t`` in ``[0, t_end]`` [s]. Raises unless the integration succeeds.
        """
        n_i = self.n_i
        self._x_warm = None

        def rhs(_t, y):
            r1, r2 = self.bulb_rates(self._full(y[:n_i]), self._full(y[n_i:]))
            return np.concatenate([r1, r2])

        y0 = np.concatenate([np.asarray(x1_0, float)[:n_i], np.asarray(x2_0, float)[:n_i]])
        sol = solve_ivp(rhs, (0.0, t_end), y0, method="DOP853", rtol=rtol, atol=atol,
                        dense_output=True)
        if sol.status != 0:
            raise RuntimeError(f"solve_ivp failed: {sol.message}")
        return sol


def fick_difference(t, x1_0, x2_0, d_ms, area, length, vol_1, vol_2):
    """Bulb difference ``x1 - x2`` of every species under Fick with Wilke diffusivities.

    Each species relaxes on its own, with the Wilke effective diffusivity at the
    mean initial composition, so every difference decays exponentially.

    Parameters
    ----------
    t : ndarray
        Times [s].
    x1_0, x2_0 : ndarray, shape (n_c,)
        Initial bulb compositions.
    d_ms : ndarray, shape (n_c, n_c)
        Pair diffusivities [m2/s].
    area, length, vol_1, vol_2 : float
        Capillary area [m2] and length [m], bulb volumes [m3].

    Returns
    -------
    ndarray, shape (len(t), n_c)
    """
    d_eff = wilke_deff(0.5 * (x1_0 + x2_0), d_ms)
    beta = (area / length) * (1.0 / vol_1 + 1.0 / vol_2) * d_eff
    return (x1_0 - x2_0)[None, :] * np.exp(-np.outer(np.asarray(t), beta))


def toor_difference(t, x1_0, x2_0, d_ms, area, length, vol_1, vol_2):
    """Bulb difference of the independent species with ``[B]`` frozen at the mean.

    Toor's linearisation: a matrix exponential of the constant
    ``(A/L)(1/V1 + 1/V2)[B]^-1``.

    Returns
    -------
    ndarray, shape (len(t), n_c - 1)
    """
    n_i = d_ms.shape[0] - 1
    b_avg = ms_b_matrix(0.5 * (x1_0 + x2_0)[None, :n_i], d_ms)[0]
    beta = (area / length) * (1.0 / vol_1 + 1.0 / vol_2) * np.linalg.inv(b_avg)
    dx0 = (x1_0 - x2_0)[:n_i]
    return np.array([expm(-beta * tt) @ dx0 for tt in np.atleast_1d(t)])
