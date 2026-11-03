"""The Davidson bubble: the l = 1 potential-flow mode around a uniform-pressure bubble.

Both the solids potential and the gas potential around a bubble in an incipiently
fluidised bed have the angular dependence cos(theta), and their radial factors h(r)
satisfy the same linear ODE

    (1/r^nu) d/dr( r^nu dh/dr ) - nu h / r^2 = 0,

nu = 2 for a spherical bubble, nu = 1 for a circular (two-dimensional) one. The
problem is linear: one sparse solve per potential, no Newton iteration.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.sparse.linalg import spsolve
from pymrm import (compute_boundary_values, construct_coefficient_matrix, construct_div,
                   construct_grad)
from pymrm.checks import find_roots

G = 980.0  # cm/s^2; the source paper works in CGS


def geometric_faces(R, r_out, n):
    """Face positions geometric in r, clustering cells at the bubble surface.

    Parameters
    ----------
    R : float
        Bubble radius.
    r_out : float
        Outer radius of the domain.
    n : int
        Number of cells.

    Returns
    -------
    numpy.ndarray
        ``n + 1`` face positions from ``R`` to ``r_out``. Doubling ``n`` keeps
        every face, so a ladder in ``n`` is nested.
    """
    return R * np.exp(np.linspace(0.0, np.log(r_out / R), n + 1))


def rise_velocity(d_b, c_ubr=0.711, g=G):
    """Single-bubble rise velocity, Kunii and Levenspiel's equation 3.

    Parameters
    ----------
    d_b : float or numpy.ndarray
        Bubble diameter [cm].
    c_ubr : float, optional
        The printed coefficient 0.711, taken as given.
    g : float, optional
        Gravitational acceleration [cm/s^2].

    Returns
    -------
    float or numpy.ndarray
        ``u_br`` [cm/s].
    """
    return c_ubr * np.sqrt(g * d_b)


def cloud_radius_exact(u_f, u_br, nu=2):
    """Closed-form cloud radius ratio R_c/R, or nan when no cloud exists.

    Parameters
    ----------
    u_f : float
        Interstitial gas velocity in the emulsion, ``u_mf / eps_mf``.
    u_br : float
        Bubble rise velocity relative to the solids.
    nu : int, optional
        2 spherical, 1 circular.

    Returns
    -------
    float
        ``[(u_br + nu u_f)/(u_br - u_f)]^(1/(nu+1))`` for ``u_br > u_f``.
    """
    if u_br <= u_f:
        return np.nan
    return ((u_br + nu * u_f) / (u_br - u_f)) ** (1.0 / (nu + 1))


class DavidsonBubble:
    """Gas and solids potentials around one bubble, on a finite-volume grid in r.

    Parameters
    ----------
    R : float
        Bubble radius.
    u_f : float
        Interstitial gas velocity far from the bubble, relative to the solids.
    u_br : float
        Bubble rise velocity, the far-field solids velocity in the bubble frame.
    n : int, optional
        Number of cells, geometric in r.
    r_out_mult : float, optional
        Outer radius as a multiple of ``R``.
    nu : int, optional
        Geometry: 2 spherical bubble, 1 circular bubble.
    far_field : {"robin", "dirichlet"}, optional
        ``"robin"`` imposes ``dh/dr + (nu/r) h = (nu+1) A``, which both radial
        modes ``A r`` and ``B r^-nu`` satisfy, so it does not truncate the
        decaying mode. ``"dirichlet"`` imposes ``h = A r`` and truncates it
        (a deliberate defect).
    nu_div, sink_coeff, far_factor, outer_a, outer_d_sign : float, optional
        Deliberate defects for the break table: the geometry index inside
        ``construct_div``, the coefficient of the ``-c h / r^2`` term, the
        far-field factor ``nu + 1``, the sign of ``a`` and of ``d`` at the
        outer boundary. ``None`` (the default) is the correct model.
    """

    def __init__(self, R, u_f, u_br, n=1600, r_out_mult=20.0, nu=2, far_field="robin",
                 nu_div=None, sink_coeff=None, far_factor=None, outer_a=1.0,
                 outer_d_sign=1.0):
        self.R, self.u_f, self.u_br, self.n, self.nu = R, u_f, u_br, n, nu
        self.far_field = far_field
        self.nu_div = nu if nu_div is None else nu_div
        self.sink_coeff = nu if sink_coeff is None else sink_coeff
        self.far_factor = nu + 1.0 if far_factor is None else far_factor
        self.outer_a, self.outer_d_sign = outer_a, outer_d_sign
        self._build_grid(r_out_mult)
        self._build_operators()

    def _build_grid(self, r_out_mult):
        self.x_f = geometric_faces(self.R, r_out_mult * self.R, self.n)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.r_out = self.x_f[-1]
        self.shape = (self.n, 1)

    def _bc(self, bubble, amplitude):
        """(left, right) bc dictionaries, ``a dh/dn + b h = d``, n the OUTWARD normal."""
        if self.far_field == "dirichlet":
            # r = r_out, defect: h = A r_out, which truncates the decaying mode B r^-nu.
            outer = {"a": 0.0, "b": 1.0, "d": amplitude * self.r_out}
        else:
            # r = r_out, outward normal +r: dh/dr + (nu/r_out) h = (nu + 1) A,
            # satisfied by both h = A r and h = B r^-nu, so no truncation error.
            outer = {"a": self.outer_a, "b": self.nu / self.r_out,
                     "d": self.outer_d_sign * self.far_factor * amplitude}
        return (bubble, outer)

    def _build_operators(self):
        """Assemble the constant operator of each potential once."""
        # GAS, relative to the solids: u_rel = -grad(phi), phi = f(r) cos(theta).
        # Uniform pressure over the bubble surface: f(R) = 0 (a = 0, b = 1, d = 0).
        # Far field: f -> -u_f r.
        # SOLIDS: v = grad(Phi_p), Phi_p = g(r) cos(theta). No solids cross the
        # bubble surface: dg/dr(R) = 0. The outward normal at r = R is -r, so
        # a dg/dn = -dg/dr; with d = 0 the sign does not matter (a = 1, b = 0).
        # Far field: g -> -u_br r.
        self.bc = {"gas": self._bc({"a": 0.0, "b": 1.0, "d": 0.0}, -self.u_f),
                   "solids": self._bc({"a": 1.0, "b": 0.0, "d": 0.0}, -self.u_br)}
        div_mat = construct_div(self.shape, self.x_f, nu=self.nu_div)  # nu: 2 spherical, 1 cylindrical
        sink_mat = construct_coefficient_matrix(-self.sink_coeff / self.x_c.reshape(-1, 1) ** 2)
        self.ops = {}
        for name, bc in self.bc.items():
            grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, bc)
            jac = (div_mat @ grad_mat + sink_mat).tocsc()
            self.ops[name] = (grad_mat, grad_bc, jac, div_mat @ grad_bc)

    def solve(self):
        """Solve both potentials (one sparse direct solve each).

        Returns
        -------
        DavidsonBubble
            ``self``, with ``f``, ``g`` (cell values), ``df``, ``dg`` (face
            gradients), ``W`` (radial velocity amplitude at faces), ``T``
            (tangential amplitude at cells) and ``dfdr_R`` (the gas potential's
            gradient at the bubble surface, from ``compute_boundary_values``).
        """
        out = {}
        for name, (grad_mat, grad_bc, jac, g_const) in self.ops.items():
            h = spsolve(jac, -g_const).reshape(self.shape)
            out[name] = (h, (grad_mat @ h + grad_bc))
        (f, df), (g, dg) = out["gas"], out["solids"]
        self.f, self.g, self.df, self.dg = f.ravel(), g.ravel(), df.ravel(), dg.ravel()
        # In the bubble frame the gas velocity is v + u_rel:
        #   w_r = W(r) cos(theta), w_theta = -T(r) sin(theta)
        self.W = self.dg - self.df
        self.T = (self.g - self.f) / self.x_c
        # Gradient at the bubble face, along +r (compute_boundary_values).
        _, self.dfdr_R = compute_boundary_values(f, self.x_f, self.x_c, self.bc["gas"][0],
                                                 bound_id=0)
        self.dfdr_R = float(np.ravel(self.dfdr_R)[0])
        return self

    def W_at(self, r):
        """Radial velocity amplitude between faces, linear in s = r^-(nu+1).

        For the two radial modes ``W = a + b r^-(nu+1)`` exactly, so this
        interpolation adds no error to the face values.
        """
        s_f = self.x_f[::-1] ** -(self.nu + 1.0)
        return np.interp(np.asarray(r, dtype=float) ** -(self.nu + 1.0), s_f, self.W[::-1])

    def cloud_radius(self):
        """Radius of the cloud, the zero of W on [R, r_out].

        Returns
        -------
        float
            The root, or nan when W has no zero and ``u_br <= u_f`` (no cloud
            exists). Raises when W has more than one zero, and when it has none
            although ``u_br > u_f`` (the cloud lies beyond ``r_out``).
        """
        res = find_roots(self.W_at, self.R, self.r_out, n_scan=256, log=True,
                         xtol=1e-13 * self.R)
        if len(res["roots"]) > 1:
            raise RuntimeError(f"W has {len(res['roots'])} zeros: {res['message']}")
        if not res["roots"]:
            if self.u_br > self.u_f:
                raise RuntimeError("no zero of W inside r_out although u_br > u_f: "
                                   "the cloud lies beyond r_out, enlarge the domain")
            return np.nan
        return res["roots"][0]

    def throughflow_per_volume(self, eps_mf):
        """Gas through-flow of the bubble per unit bubble volume [1/s].

        Gas crosses the bubble surface at ``u_rel`` only (the solids have no
        normal velocity there), so the superficial outflow over the upper half
        is ``eps_mf (-f'(R))`` times ``pi R^2`` (sphere) or ``2 R`` per unit
        length (circle).
        """
        if self.nu == 2:
            return eps_mf * (-self.dfdr_R) * np.pi * self.R**2 / (4.0 / 3.0 * np.pi * self.R**3)
        return eps_mf * (-self.dfdr_R) * 2.0 * self.R / (np.pi * self.R**2)

    def integral_identity(self, field="gas"):
        """Relative gap of the once-integrated radial equation on the discrete field.

        ``sum(vol nu h / r^2) = [r^nu dh/dr]`` over the domain, with the
        model's own ``nu``. It telescopes for any consistent operator and any
        boundary data, so it only detects a mismatch between the geometry
        index in ``construct_div`` and the sink coefficient.
        """
        h, dh = (self.f, self.df) if field == "gas" else (self.g, self.dg)
        x_f, nu = self.x_f, self.nu
        vol = (x_f[1:] ** (nu + 1) - x_f[:-1] ** (nu + 1)) / (nu + 1)
        lhs = np.sum(vol * nu * h / self.x_c**2)
        rhs = x_f[-1] ** nu * dh[-1] - x_f[0] ** nu * dh[0]
        return abs(lhs - rhs) / abs(rhs)

    def f_exact(self, r):
        """Closed-form gas radial factor, ``-u_f (r - R^(nu+1) / r^nu)``."""
        return -self.u_f * (r - self.R ** (self.nu + 1) / r**self.nu)


def cloud_radius_by_trajectories(bubble, th0=0.1, t_max=60.0, n_halvings=45):
    """Cloud radius from the fate of gas trajectories (an independent route).

    Gas elements are released at angle ``th0`` from the upward axis and
    integrated through cubic-spline interpolants of the computed velocity
    amplitudes; the release radius is bisected until a trajectory stops
    returning to the bubble.

    Parameters
    ----------
    bubble : DavidsonBubble
        A solved spherical bubble with a cloud.
    th0 : float, optional
        Release angle [rad].
    t_max : float, optional
        Integration time after which a trajectory counts as not returning.
    n_halvings : int, optional
        Number of bisection steps.

    Returns
    -------
    float
        The release radius separating returning from escaping trajectories.
    """
    R, Rc = bubble.R, bubble.cloud_radius()
    s_df, s_dg = CubicSpline(bubble.x_f, bubble.df), CubicSpline(bubble.x_f, bubble.dg)
    s_f = CubicSpline(bubble.x_c, bubble.f, extrapolate=True)
    s_g = CubicSpline(bubble.x_c, bubble.g, extrapolate=True)

    def rhs(t, y):
        r, th = y
        w_r = (s_dg(r) - s_df(r)) * np.cos(th)
        w_t = -((s_g(r) - s_f(r)) / r) * np.sin(th)
        return [w_r, w_t / r]

    def hit(t, y):
        return y[0] - R
    hit.terminal, hit.direction = True, -1

    def gone(t, y):
        return y[0] - 4 * Rc
    gone.terminal, gone.direction = True, 1

    def returns(r0):
        s = solve_ivp(rhs, (0, t_max), [r0, th0], events=[hit, gone], rtol=1e-10,
                      atol=1e-12, max_step=0.05)
        return len(s.t_events[0]) > 0

    lo, hi = R * 1.001, Rc * 3
    for _ in range(n_halvings):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if returns(mid) else (lo, mid)
    return 0.5 * (lo + hi)
