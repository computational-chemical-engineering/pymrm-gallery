"""Absorption into a falling liquid film, marched in the axial coordinate.

Model code of gallery page A3.2 (Higbie's penetration model, read from Bird,
Stewart & Lightfoot section 18.5). Dimensionless: ``xi = x/delta`` across the
film from the gas-liquid interface, ``zeta = z/L`` down it, ``u = c_A/c_A0`` and
``Lam = D L/(v_max delta^2)``. BSL Eq. 18.5-7 then reads

    cap(xi) du/dzeta = Lam d2u/dxi2,   u(xi, 0) = 0,   u(0, zeta) = 1,

with ``cap = 1 - xi^2`` (the parabolic velocity profile) and a no-flux wall at
``xi = 1``. BSL's penetration form, Eq. 18.5-11, has ``cap = 1`` and no wall.

Contents
--------
Film
    The finite-volume march: the transverse operator from ``construct_grad`` and
    ``construct_div``, assembled once; backward Euler (default), BDF2 or
    Crank-Nicolson in ``zeta`` on steps uniform in ``sqrt(zeta)``.
phi_extrapolated
    The outlet mixing cup by backward Euler at ``n_t`` and ``2 n_t`` steps, and
    its Richardson extrapolation.
phi_penetration, phi_wall_series
    The two exact answers for plug flow: semi-infinite (BSL Eq. 18.5-16
    integrated) and with a no-flux wall (eigenfunction series).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from pymrm import compute_boundary_values, construct_div, construct_grad

__all__ = ["Film", "phi_extrapolated", "phi_penetration", "phi_wall_series"]


class Film:
    """Falling-film absorption, marched in ``zeta`` from 0 to 1.

    The operators are built on the field shape ``(n_x, 1)``: space first, one
    field. The march carries the state as the flat vector of that shape,
    ``(n_x,)``, which is what ``solve_banded`` takes and returns. The
    transverse operator ``Lam d2/dxi2`` and its boundary vector are assembled
    once, in the constructor, from ``construct_grad`` and ``construct_div``. The
    step matrix ``cap/dzeta - Lam d2/dxi2`` changes only through ``dzeta``, so
    the tridiagonal operator is also stored once in banded form and each step
    adds ``cap/dzeta`` to its main diagonal before one ``solve_banded`` call.

    Parameters
    ----------
    Lam : float
        ``D L/(v_max delta^2)``; ``sqrt(4 Lam)`` is the penetration depth over
        the film thickness at the outlet.
    XI : float, optional
        Domain depth in units of ``delta``: 1 for the film, or several
        penetration depths for the semi-infinite penetration problem.
    n_x : int, optional
        Number of cells across the film.
    stretch : float, optional
        Exponential grid stretch, fine at the interface; the map of the cell
        index is fixed, so doubling ``n_x`` halves every cell. 0 is uniform.
    profile : {"parabolic", "plug"}, optional
        Capacity ``1 - xi^2`` (BSL Eq. 18.5-7) or 1 (Eq. 18.5-11).
    wall : {"noflux", "sink"}, optional
        Condition at ``xi = XI``: no flux (BSL Eq. 18.5-10, or the far field of
        Eq. 18.5-14 on a deep domain), or ``u = 0``.
    """

    def __init__(self, Lam, XI=1.0, n_x=600, stretch=6.0, profile="parabolic",
                 wall="noflux"):
        if profile not in ("parabolic", "plug"):
            raise ValueError("profile must be 'parabolic' or 'plug'")
        if wall not in ("noflux", "sink"):
            raise ValueError("wall must be 'noflux' or 'sink'")
        self.Lam, self.XI, self.n_x = float(Lam), float(XI), int(n_x)
        self.shape = (self.n_x, 1)                 # (space, field), never (n_x,)
        self._build_grid(stretch)
        self.cap = (1.0 - self.x_c ** 2) if profile == "parabolic" else np.ones(self.n_x)
        # Outward normal n, a du/dn + b u = d:
        #   xi = 0  (n = -xi): gas-liquid interface, u = 1 (BSL Eq. 18.5-9)
        #   xi = XI (n = +xi): wall, du/dxi = 0 (BSL Eq. 18.5-10), or u = 0
        self.bc = ({"a": 0.0, "b": 1.0, "d": 1.0},
                   {"a": 1.0, "b": 0.0, "d": 0.0} if wall == "noflux"
                   else {"a": 0.0, "b": 1.0, "d": 0.0})
        self._build_operators()

    def _build_grid(self, stretch):
        s = np.linspace(0.0, 1.0, self.n_x + 1)
        self.x_f = (self.XI * np.expm1(stretch * s) / np.expm1(stretch) if stretch
                    else self.XI * s)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.dx = np.diff(self.x_f)

    def _build_operators(self):
        grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc)
        div_mat = construct_div(self.shape, self.x_f, nu=0)   # nu = 0: Cartesian slab
        self.jac_diff = (self.Lam * (div_mat @ grad_mat)).tocsr()
        self.g_const = self.Lam * (div_mat @ grad_bc).toarray().ravel()
        # the same tridiagonal operator in solve_banded's (3, n_x) layout
        self.ab_diff = np.zeros((3, self.n_x))
        self.ab_diff[0, 1:] = self.jac_diff.diagonal(1)
        self.ab_diff[1] = self.jac_diff.diagonal(0)
        self.ab_diff[2, :-1] = self.jac_diff.diagonal(-1)

    def cup(self, u):
        """Mixing cup ``int cap u dxi`` of a profile (BSL Problem 18C.3)."""
        return float(self.dx @ (self.cap * np.ravel(u)))

    def interface_flux(self, u):
        """``-Lam du/dxi`` at ``xi = 0`` for one profile or a stack ``(n_x, k)``."""
        u = np.asarray(u, dtype=float).reshape(self.n_x, -1)
        _, grad_0, _, _ = compute_boundary_values(u, self.x_f, self.x_c, self.bc)
        return -self.Lam * np.ravel(grad_0)    # the gradient is along +xi

    def _solve(self, a0_over_dz, rhs, theta=1.0):
        ab = -theta * self.ab_diff
        ab[1] += a0_over_dz * self.cap
        return solve_banded((1, 1), ab, rhs, check_finite=False)

    def march(self, n_t=800, scheme="euler", n_start=4, spacing="sqrt",
              record_flux=False):
        """March from ``zeta = 0`` to 1 in ``n_t`` steps.

        Parameters
        ----------
        n_t : int, optional
            Number of steps.
        scheme : {"euler", "bdf2", "cn"}, optional
            Backward Euler; variable-step BDF2 with a backward-Euler first
            step; or Crank-Nicolson after ``n_start`` backward-Euler steps.
        n_start : int, optional
            Backward-Euler start-up steps of ``"cn"``.
        spacing : {"sqrt", "uniform"}, optional
            Steps uniform in ``sqrt(zeta)`` or in ``zeta``.
        record_flux : bool, optional
            Also integrate the interfacial flux over the march, with the flux
            of every step from ``compute_boundary_values`` (backward Euler
            only: it is the quadrature that matches that scheme).

        Returns
        -------
        dict
            ``u`` the outlet profile, ``cup`` its mixing cup and, with
            ``record_flux``, ``flux_integral``, the step lengths ``dz`` and the
            profile after every step, ``steps`` with shape ``(n_x, n_t)``.
        """
        w = np.linspace(0.0, 1.0, n_t + 1)
        z = w ** 2 if spacing == "sqrt" else w
        u = np.zeros(self.n_x)
        u_prev = None
        steps = np.empty((self.n_x, n_t)) if record_flux else None
        if record_flux and scheme != "euler":
            raise ValueError("record_flux is the backward-Euler quadrature")
        for j in range(1, n_t + 1):
            dz = z[j] - z[j - 1]
            if scheme == "euler" or j == 1 or (scheme == "cn" and j <= n_start):
                u_new = self._solve(1.0 / dz, self.cap * u / dz + self.g_const)
            elif scheme == "bdf2":
                om = dz / (z[j - 1] - z[j - 2])
                a0 = (1 + 2 * om) / (1 + om)
                rhs = self.cap * ((1 + om) * u - om ** 2 / (1 + om) * u_prev) / dz
                u_new = self._solve(a0 / dz, rhs + self.g_const)
            elif scheme == "cn":
                rhs = self.cap * u / dz + 0.5 * (self.jac_diff @ u) + self.g_const
                u_new = self._solve(1.0 / dz, rhs, theta=0.5)
            else:
                raise ValueError("scheme must be 'euler', 'bdf2' or 'cn'")
            u_prev, u = u, u_new
            if record_flux:
                steps[:, j - 1] = u
        out = {"u": u, "cup": self.cup(u)}
        if record_flux:
            out["flux_integral"] = float(np.diff(z) @ self.interface_flux(steps))
            out["dz"], out["steps"] = np.diff(z), steps
        return out


def phi_extrapolated(Lam, n_t=800, **film):
    """Mixing cup by backward Euler at ``n_t`` and ``2 n_t`` steps.

    Backward Euler is first order in the step, so ``2 Phi(2 n_t) - Phi(n_t)``
    removes the leading error.

    Parameters
    ----------
    Lam : float
        As in :class:`Film`.
    n_t : int, optional
        The coarser of the two step counts.
    **film
        Further :class:`Film` arguments.

    Returns
    -------
    dict
        ``raw`` (``n_t`` steps), ``fine`` (``2 n_t``) and ``extrapolated``.
    """
    f = Film(Lam, **film)
    raw, fine = f.march(n_t)["cup"], f.march(2 * n_t)["cup"]
    return {"raw": raw, "fine": fine, "extrapolated": 2.0 * fine - raw}


def phi_penetration(Lam):
    """``sqrt(4 Lam/pi)``: BSL Eq. 18.5-16 integrated across a semi-infinite film.

    Parameters
    ----------
    Lam : float or numpy.ndarray

    Returns
    -------
    float or numpy.ndarray
    """
    return np.sqrt(4.0 * Lam / np.pi)


def phi_wall_series(Lam, n_terms=20000):
    """Exact mixing cup for plug flow with a no-flux wall at ``xi = 1``.

    Separation of variables on ``v = 1 - u`` (``v(0) = 0``, ``dv/dxi(1) = 0``,
    ``v(xi, 0) = 1``) gives ``Phi = 1 - sum_n 2/lam_n^2 exp(-Lam lam_n^2)``
    with ``lam_n = (n + 1/2) pi``. It shares no code with :class:`Film`.

    Parameters
    ----------
    Lam : float
    n_terms : int, optional
        Number of terms kept.

    Returns
    -------
    float
    """
    lam = (np.arange(n_terms) + 0.5) * np.pi
    return 1.0 - float(np.sum(2.0 / lam ** 2 * np.exp(-Lam * lam ** 2)))
