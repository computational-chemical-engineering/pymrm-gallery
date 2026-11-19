"""Growth laws, a batch marcher, a plug-flow fermenter and a chemostat for J4.1.

The five printed growth laws are written as mu = mu_m g(theta; S, X), with
theta the parameters that do not enter linearly. ``BatchMarch`` and
``PlugFlowFermenter`` are built from pymrm operators; ``Chemostat`` and
``batch_reference`` are the same source term integrated by ``solve_ivp``.
No plotting and no data loading: the page passes every datum in.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import identity as speye

from pymrm import (construct_convflux_upwind, construct_div, NumJac, newton,
                   compute_boundary_values)
from pymrm.checks import find_roots

# mu / mu_m for the five printed laws and the zero-order null; theta is a list.
GROWTH = {
    "Monod": lambda th, S, X: S/(th[0] + S),
    "Blackman": lambda th, S, X: np.minimum(1.0, S/(2*th[0])),
    "Tessier": lambda th, S, X: 1 - np.exp(-th[0]*S),
    "Moser": lambda th, S, X: S**th[1]/(th[0] + S**th[1]),
    "Contois": lambda th, S, X: S/(th[0]*X + S),
    "zero order": lambda th, S, X: np.ones_like(S),
}


def integrate_interval(rhs, t_end, y0, rtol, atol, method="LSODA", dense_output=False):
    """Integrate ``dy/dt = rhs(t, y)`` from 0 to ``t_end`` and check the solve.

    Parameters
    ----------
    rhs : callable
        ``rhs(t, y) -> list``.
    t_end : float
        End time.
    y0 : sequence of float
        Initial state.
    rtol, atol : float
        Integrator tolerances.
    method : str, optional
        ``solve_ivp`` method.
    dense_output : bool, optional
        Return the dense interpolant as ``sol``.

    Returns
    -------
    OdeResult
        The ``solve_ivp`` result.

    Raises
    ------
    RuntimeError
        When ``solve_ivp`` reports a status other than 0.
    """
    s = solve_ivp(rhs, [0.0, t_end], list(y0), rtol=rtol, atol=atol, method=method,
                  dense_output=dense_output)
    if s.status != 0:
        raise RuntimeError(f"solve_ivp failed (status {s.status}): {s.message}")
    return s


def one_root(f, lo, hi, *, n_scan=64, log=False, xtol=1e-12, what="root"):
    """The single root of ``f`` on the physical range ``[lo, hi]``.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        Physical range scanned for sign changes.
    n_scan : int, optional
        Scan points.
    log : bool, optional
        Scan on a logarithmic grid.
    xtol : float, optional
        Absolute tolerance of the refined root.
    what : str, optional
        Name used in the error message.

    Returns
    -------
    float
        The root.

    Raises
    ------
    RuntimeError
        Unless ``pymrm.checks.find_roots`` finds exactly one root and no failed
        scan point.
    """
    r = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(r["roots"]) != 1 or r["failed"]:
        raise RuntimeError(f"{what}: {len(r['roots'])} roots on [{lo:g}, {hi:g}],"
                           f" {len(r['failed'])} failed points ({r['message']})")
    return float(r["roots"][0])


def _source(growth, phi, c, sign=1.0):
    """Source (dS/dt, dX/dt) = (-sign m/phi, m) with m = growth(S, X) X."""
    S = np.clip(c[..., 0:1], 0.0, None)
    m = growth(S, c[..., 1:2])*c[..., 1:2]
    return np.concatenate([-sign*m/phi, m], axis=-1)


class BatchMarch:
    """Backward Euler on the well-mixed batch (S, X) with ``newton`` and ``NumJac``.

    The state has shape (1, 2): no spatial extent and two fields, and the
    growth law is pointwise in those fields, which is the default ``NumJac``
    stencil (last axis coupled in full).

    Parameters
    ----------
    growth : callable
        ``growth(S, X) -> mu``, the specific growth rate in 1/h.
    phi : float
        Yield, biomass formed per substrate consumed.
    nstep : int, optional
        Number of equal time steps.
    newton_kw : dict, optional
        Extra keyword arguments for ``newton`` (tolerances).
    """

    def __init__(self, growth, phi, nstep=200, newton_kw=None):
        self.g, self.phi, self.nstep = growth, phi, nstep
        self.newton_kw = {"maxfev": 100, **(newton_kw or {})}
        self.shape = (1, 2)                      # never a bare (2,)
        self.numjac = NumJac(self.shape)         # last axis coupled in full

    def source(self, c):
        """Source term of the two fields."""
        return _source(self.g, self.phi, c)

    def step(self, c_old, h):
        """One backward-Euler step of length ``h``; raises when newton fails."""
        eye = speye(c_old.size, format="csc")

        def res(c):
            s, js = self.numjac(self.source, c)
            return ((c.reshape(-1, 1) - c_old.reshape(-1, 1))/h - s.reshape(-1, 1),
                    eye/h - js)
        out = newton(res, c_old, **self.newton_kw)
        if not out.success:
            raise RuntimeError(f"batch step did not converge: {out.message}")
        return out.x.reshape(self.shape)

    def march(self, c0, T):
        """March the state ``c0`` (S, X) over ``[0, T]``; returns shape (1, 2)."""
        c = np.asarray(c0, float).reshape(self.shape)
        h = T/self.nstep
        for _ in range(self.nstep):
            c = self.step(c, h)
        return c


class PlugFlowFermenter:
    """Steady 1-D plug flow of (S, X) at constant yield, from pymrm operators.

    ``construct_convflux_upwind`` and ``construct_div`` assemble the
    convection operator once; ``NumJac((n, 2))`` differentiates the pointwise
    growth source.

    Parameters
    ----------
    growth : callable
        ``growth(S, X) -> mu`` in 1/h.
    phi : float
        Yield.
    Sf, Xf : float
        Feed substrate and biomass.
    tau : float
        Residence time; the length is ``v tau``.
    ncell : int, optional
        Number of cells.
    v : float, optional
        Velocity.
    nu : int, optional
        Geometry of ``construct_div``: 0 Cartesian (the plug-flow tube).
    substrate_sign : float, optional
        Sign of the substrate consumption term; 1 is the physical model.
    newton_kw : dict, optional
        Extra keyword arguments for ``newton`` (tolerances).
    """

    def __init__(self, growth, phi, Sf, Xf, tau, ncell=400, v=1.0, nu=0,
                 substrate_sign=1.0, newton_kw=None):
        self.newton_kw = {"maxfev": 200, **(newton_kw or {})}
        self.g, self.phi, self.v, self.tau = growth, phi, v, tau
        self.sign = substrate_sign
        self.shape = (ncell, 2)                       # never a bare (ncell,)
        self.z_f = np.linspace(0.0, v*tau, ncell + 1)
        self.z_c = 0.5*(self.z_f[:-1] + self.z_f[1:])
        d0 = np.array([[Sf, Xf]])
        # OUTWARD normal, a dc/dn + b c = d:
        #   z = 0, inlet: c = (Sf, Xf), a Dirichlet feed  a=0, b=1, d=(Sf, Xf)
        #   z = L, outlet: pure convective outflow, the face carries the last
        #          cell's value  {"outflow": True}
        self.bc = ({"a": 0.0, "b": 1.0, "d": d0}, {"outflow": True})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                  self.bc, v=v, axis=0)
        self.div = construct_div(self.shape, self.z_f, nu=nu, axis=0)  # nu=0 Cartesian
        self.A, self.b = self.div @ conv, self.div @ conv_bc
        self.numjac = NumJac(self.shape)               # pointwise in the 2 fields
        self.c0 = np.tile(d0, (ncell, 1))

    def source(self, c):
        """Source term of the two fields."""
        return _source(self.g, self.phi, c, self.sign)

    def residual(self, c):
        """Steady residual div(v c) - source and its Jacobian."""
        s, js = self.numjac(self.source, c)
        return (self.b + self.A @ c.reshape((-1, 1)) - s.reshape((-1, 1)),
                self.A - js)

    def solve(self):
        """Solve the steady state from the feed; raises when newton fails."""
        r = newton(self.residual, self.c0, **self.newton_kw)
        if not r.success:
            raise RuntimeError(f"plug-flow solve did not converge: {r.message}")
        self.c = r.x.reshape(self.shape)
        return self

    def outlet(self):
        """Outlet (S, X) from ``compute_boundary_values`` under the same bc."""
        return np.asarray(compute_boundary_values(self.c, self.z_f, self.z_c,
                                                  self.bc, axis=0)[2]).reshape(-1)[:2]


def batch_reference(growth, phi, S0, X0, tau, rtol=1e-12):
    """The batch state at ``t = tau`` by LSODA (Levenspiel Fig. 29.1).

    Parameters
    ----------
    growth : callable
        ``growth(S, X) -> mu`` in 1/h.
    phi : float
        Yield.
    S0, X0 : float
        Initial state.
    tau : float
        End time.
    rtol : float, optional
        Relative tolerance; ``atol`` is ``1e-2 rtol``.

    Returns
    -------
    ndarray
        (S, X) at ``tau``.
    """
    def rhs(t, y):
        m = growth(max(y[0], 0.0), y[1])*y[1]
        return [-m/phi, m]
    return integrate_interval(rhs, tau, [S0, X0], rtol, rtol*1e-2).y[:, -1]


class Chemostat:
    """Rawlings eq. (10.18), marched to steady state by LSODA.

    Parameters
    ----------
    growth : callable
        ``growth(S, X) -> mu`` in 1/h.
    y_yield : float
        Yield.
    Sf : float
        Feed substrate.
    D : float
        Dilution rate in 1/h.
    """

    def __init__(self, growth, y_yield, Sf, D):
        self.g, self.y, self.Sf, self.D = growth, y_yield, Sf, D

    def rhs(self, t, c):
        """Right-hand side of the two balances."""
        S, X = max(c[0], 0.0), c[1]
        m = self.g(S, X)
        return [self.D*(self.Sf - S) - m*X/self.y, (m - self.D)*X]

    def steady(self, tend=40000.0, rtol=1e-11):
        """(S, X) at ``tend`` from S = Sf and X = 0.1 Sf y."""
        return integrate_interval(self.rhs, tend, [self.Sf, 0.1*self.Sf*self.y],
                                  rtol, rtol*1e-3).y[:, -1]
