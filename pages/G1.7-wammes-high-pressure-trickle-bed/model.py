"""Wammes and Westerterp (1991): trickle-bed hydrodynamics and the CO2 absorption balance.

Model code of gallery page G1.7. Two pieces:

- the correlation pair of the paper, Eq. 8 (pressure gradient, containing the
  hold-up) and Eq. 9 (dynamic hold-up, containing the pressure gradient),
  closed as one algebraic system per operating point and solved for many
  operating points at once with ``pymrm.NumJac`` and ``pymrm.newton``;
- the plug-flow gas-phase CO2 balance of Eqs. 1 and 2, the model behind the
  paper's interfacial-area measurements, discretised with pymrm operators.

Contents
--------
gas_viscosity
    Table 2 gas viscosity, linear in pressure between its 1 and 100 bar values.
OperatingPoints
    Gas and liquid velocities, gas density and viscosity per point, and the
    liquid and packing constants shared by all points.
closure_residual
    Residual of Eq. 9 with Eq. 8 substituted, and the pressure gradient.
HoldupPressureClosure
    The coupled pair for many points: ``solve`` (``newton``), ``solve_scalar``
    (``brentq`` per point, an independent route) and ``sign_changes`` (a scan
    for further roots).
CO2PlugFlow
    The steady plug-flow balance of Eq. 2 with the linear sink of Eq. 1.
conversion_exact, a_gl_from_conversion
    Eq. 2 integrated, and its inversion Eq. 3.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import eye_array
from scipy.sparse.linalg import spsolve

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, newton)

__all__ = ["G", "gas_viscosity", "OperatingPoints", "closure_residual",
           "HoldupPressureClosure", "CO2PlugFlow", "conversion_exact",
           "a_gl_from_conversion"]

G = 9.81  # m/s2


def gas_viscosity(rho_g, rho_per_bar, eta_1bar, eta_100bar):
    """Gas viscosity from Table 2, linear in pressure between 1 and 100 bar.

    Parameters
    ----------
    rho_g : float or numpy.ndarray
        Gas density (kg/m3).
    rho_per_bar : float
        Table 2's density coefficient, ``rho = rho_per_bar * P`` with ``P`` in bar.
    eta_1bar, eta_100bar : float
        Viscosity at 1 and 100 bar (Pa s).

    Returns
    -------
    numpy.ndarray
        Viscosity (Pa s), same shape as ``rho_g``.
    """
    p_bar = np.asarray(rho_g, float) / rho_per_bar
    return eta_1bar + (eta_100bar - eta_1bar) * (p_bar - 1.0) / 99.0


@dataclass(frozen=True)
class OperatingPoints:
    """Operating points of one gas-liquid-packing system.

    Parameters
    ----------
    v_g, v_l : numpy.ndarray
        Superficial gas and liquid velocities (m/s), one entry per point.
    rho_g, eta_g : numpy.ndarray
        Gas density (kg/m3) and viscosity (Pa s) per point.
    rho_l, eta_l : float
        Liquid density (kg/m3) and viscosity (Pa s).
    eps : float
        Bed porosity.
    d_p : float
        Particle diameter (m).
    eps_beta_r : float
        Residual hold-up per unit bed volume, ``eps * beta_r`` (0.06 in the paper).
    a_s_per_bed : bool
        True: ``a_s = 6 (1 - eps) / d_p`` per unit bed volume, the paper's
        convention. False: per unit particle volume, ``6 / d_p``.
    """

    v_g: np.ndarray
    v_l: np.ndarray
    rho_g: np.ndarray
    eta_g: np.ndarray
    rho_l: float
    eta_l: float
    eps: float
    d_p: float
    eps_beta_r: float = 0.06
    a_s_per_bed: bool = True

    def __post_init__(self):
        n = np.size(self.v_g)
        for name in ("v_g", "v_l", "rho_g", "eta_g"):
            arr = np.broadcast_to(np.asarray(getattr(self, name), float), (n,))
            object.__setattr__(self, name, np.array(arr))

    @property
    def n(self) -> int:
        """Number of operating points."""
        return self.v_g.size

    @property
    def beta_r(self) -> float:
        """Residual hold-up per unit void volume."""
        return self.eps_beta_r / self.eps

    @property
    def a_s(self) -> float:
        """Packing surface area (m2/m3) in the convention of ``a_s_per_bed``."""
        return 6.0 * (1.0 - self.eps) / self.d_p if self.a_s_per_bed else 6.0 / self.d_p

    def with_(self, **changes) -> "OperatingPoints":
        """A copy with some fields changed (``dataclasses.replace``)."""
        return replace(self, **changes)


def closure_residual(b_dyn, op: OperatingPoints):
    """Residual of Eq. 9 with Eq. 8's pressure gradient substituted.

    Parameters
    ----------
    b_dyn : numpy.ndarray
        Dynamic hold-up per unit void volume, shape ``(n, 1)`` (or anything of
        size ``n``).
    op : OperatingPoints
        The operating points.

    Returns
    -------
    residual : numpy.ndarray
        ``beta_dyn - Eq. 9``, shape ``(n, 1)``.
    dPL : numpy.ndarray
        Two-phase pressure gradient of Eqs. 5 and 8 (N/m3), shape ``(n, 1)``.
    """
    b = np.reshape(b_dyn, (-1, 1))
    v_g, v_l = op.v_g[:, None], op.v_l[:, None]
    rho_g, eta_g = op.rho_g[:, None], op.eta_g[:, None]
    eps, d_p = op.eps, op.d_p
    re_l = op.rho_l * v_l * d_p / op.eta_l
    ga_l = d_p ** 3 * op.rho_l ** 2 * G / op.eta_l ** 2
    beta_t = op.beta_r + b
    u_g = v_g / (eps * (1.0 - beta_t))                       # Eq. 6
    re_g = rho_g * u_g * eps * d_p / (eta_g * (1.0 - eps))   # Eq. 8 Reynolds group
    dPL = (155.0 * re_g ** -0.37 * 0.5 * rho_g * u_g ** 2 / d_p
           * (1.0 - eps) / (eps * (1.0 - beta_t)))           # Eqs. 5 and 8
    rhs = (3.8 * re_l ** 0.55
           * (ga_l * (1.0 + dPL / (op.rho_l * G))) ** -0.42
           * (op.a_s * d_p / eps) ** 0.65)                   # Eq. 9
    return b - rhs, dPL


class HoldupPressureClosure:
    """Eqs. 8 and 9 closed as one algebraic system per operating point.

    Each point is an independent scalar equation in its dynamic hold-up, so the
    Jacobian is diagonal; ``NumJac((n, 1))`` builds it, once per instance.

    Parameters
    ----------
    op : OperatingPoints
        The operating points.
    """

    def __init__(self, op: OperatingPoints):
        self.op = op
        self.shape = (op.n, 1)
        self.numjac = NumJac(self.shape)   # pointwise: diagonal Jacobian
        # the residual has u_g = v_g / (eps (1 - beta_t)), undefined for
        # beta_t >= 1: an upper bound keeps every Newton iterate inside
        self.b_max = 1.0 - op.beta_r

    def residual(self, b):
        """Residual and Jacobian for ``newton``."""
        g, jac = self.numjac(lambda x: closure_residual(x, self.op)[0], b.reshape(self.shape))
        return g.reshape(-1, 1), jac

    def solve(self, guess=0.2, tol=1e-10, **newton_kwargs):
        """Solve all points with ``pymrm.newton``.

        Parameters
        ----------
        guess : float
            Initial dynamic hold-up at every point.
        tol : float
            Newton step tolerance; the unknowns are of order 0.1, so 1e-10 is
            far above their round-off and far below every reported digit.
        **newton_kwargs
            Passed to ``newton``.

        Returns
        -------
        dict
            ``beta_dyn``, ``beta_t``, ``dPL`` (N/m3), ``Re_g`` (Eq. 8 group),
            ``Re_l``, ``dP_grav`` (``dPL / (rho_l g)``), and ``result``.

        Raises
        ------
        RuntimeError
            If ``newton`` returns without convergence.
        """
        kwargs = {"tol": tol, "bounds": (None, self.b_max)}
        kwargs.update(newton_kwargs)
        res = newton(self.residual, np.full(self.shape, guess), **kwargs)
        if not res.success:
            raise RuntimeError(f"closure did not converge: status {res.status}, {res.message}")
        return self.diagnostics(res.x, result=res)

    def solve_scalar(self, xtol=1e-14):
        """Solve each point separately with ``brentq``, an independent route.

        The bracket is ``[0, 1 - beta_r)``: the residual is negative at zero
        hold-up (Eq. 9 is positive) and tends to ``+inf`` as ``beta_t -> 1``.
        """
        b = np.empty(self.op.n)
        for i in range(self.op.n):
            one = self._point(i)
            b[i] = brentq(lambda x: float(closure_residual(np.array([x]), one)[0][0, 0]),
                          0.0, self.b_max * (1.0 - 1e-9), xtol=xtol, rtol=4 * np.finfo(float).eps)
        return self.diagnostics(b)

    def sign_changes(self, n_scan=2000):
        """Sign changes of each point's residual on a uniform scan of ``[0, 1 - beta_r)``.

        Returns
        -------
        numpy.ndarray
            Number of sign changes per point (1 for a unique root in the scan).
        """
        b = np.linspace(0.0, self.b_max, n_scan + 1)[:-1]
        counts = np.empty(self.op.n, dtype=int)
        for i in range(self.op.n):
            one = self._point(i)
            r = closure_residual(b, one.with_(v_g=np.full(b.size, one.v_g[0]),
                                              v_l=np.full(b.size, one.v_l[0]),
                                              rho_g=np.full(b.size, one.rho_g[0]),
                                              eta_g=np.full(b.size, one.eta_g[0])))[0].ravel()
            counts[i] = int(np.sum(np.signbit(r[1:]) != np.signbit(r[:-1])))
        return counts

    def _point(self, i) -> OperatingPoints:
        op = self.op
        return op.with_(v_g=op.v_g[i:i + 1], v_l=op.v_l[i:i + 1],
                        rho_g=op.rho_g[i:i + 1], eta_g=op.eta_g[i:i + 1])

    def diagnostics(self, b_dyn, result=None):
        """Hold-up, pressure gradient and the dimensionless groups at ``b_dyn``."""
        op = self.op
        b = np.reshape(b_dyn, -1)
        dPL = closure_residual(b, op)[1].ravel()
        beta_t = op.beta_r + b
        u_g = op.v_g / (op.eps * (1.0 - beta_t))
        return {"beta_dyn": b, "beta_t": beta_t, "dPL": dPL,
                "Re_g": op.rho_g * u_g * op.eps * op.d_p / (op.eta_g * (1.0 - op.eps)),
                "Re_l": op.rho_l * op.v_l * op.d_p / op.eta_l,
                "dP_grav": dPL / (op.rho_l * G), "result": result}


def conversion_exact(a_gl, K, length, v_g):
    """CO2 conversion from Eq. 2 integrated: ``1 - exp(-K a_gl L / v_g)``."""
    return 1.0 - np.exp(-K * a_gl * length / v_g)


def a_gl_from_conversion(xi, K, length, v_g):
    """Eq. 3: interfacial area (m2/m3) from a measured conversion."""
    return -v_g / (K * length) * np.log(1.0 - xi)


class CO2PlugFlow:
    """Steady gas-phase CO2 balance, Eqs. 1 and 2, normalised by the inlet value.

    ``v_g dc/dx = -K a_gl c`` on ``0 < x < L`` with ``c(0) = 1``, first-order
    upwind convection, solved as one sparse linear system.

    Parameters
    ----------
    a_gl : float
        Gas-liquid interfacial area (m2/m3).
    K : float
        ``m sqrt(k_1p c_DEA^p D_CO2)`` of Table 1 (m/s).
    length : float
        Bed length (m).
    v_g : float
        Superficial gas velocity (m/s).
    n_x : int
        Number of cells.
    """

    def __init__(self, a_gl, K, length, v_g, n_x=200):
        self.a_gl, self.K, self.length, self.v_g, self.n_x = a_gl, K, length, v_g, n_x
        self.shape = (n_x, 1)
        self._build_grid()
        self._build_operators()

    def _build_grid(self):
        self.x_f = np.linspace(0.0, self.length, self.n_x + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # bc on the OUTWARD normal, a dc/dn + b c = d:
        #   inlet  x = 0: c = 1           -> a = 0, b = 1, d = 1
        #   outlet x = L: dc/dx = 0       -> a = 1, b = 0, d = 0 (the convective
        #                                    operator needs a value on the exit face)
        self.bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})

    def _build_operators(self):
        conv_mat, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                      self.bc, v=self.v_g)
        div_mat = construct_div(self.shape, self.x_f, nu=0)   # nu = 0: Cartesian axial
        jac_sink = self.K * self.a_gl * eye_array(self.n_x, format="csc")  # J a_gl = K a_gl c
        self.jac_const = (div_mat @ conv_mat + jac_sink).tocsc()
        self.g_const = (div_mat @ conv_bc).toarray().reshape(-1, 1)

    def solve(self):
        """Solve the linear system; returns ``c`` of shape ``(n_x, 1)``."""
        self.c = spsolve(self.jac_const, -self.g_const.ravel()).reshape(self.shape)
        return self.c

    def outlet(self):
        """CO2 on the exit face, the value the operator transports (``compute_boundary_values``)."""
        _, _, c_out, _ = compute_boundary_values(self.c, self.x_f, self.x_c, self.bc)
        return float(np.ravel(c_out)[0])

    def conversion(self):
        """Gas-phase CO2 conversion ``1 - c(L)``, from the exit-face value."""
        return 1.0 - self.outlet()

    def balance_gap(self, outlet="face"):
        """Gap of the CO2 balance: inlet = outlet + absorbed, per unit ``v_g``.

        ``|1 - c_out - K a_gl sum(c dx) / v_g|`` with ``c_out`` on the exit face
        (``outlet="face"``, zero to round-off) or in the last cell (``"cell"``).
        """
        c_out = self.outlet() if outlet == "face" else float(self.c[-1, 0])
        absorbed = self.K * self.a_gl * float(np.sum(self.c[:, 0] * np.diff(self.x_f)))
        return abs(1.0 - c_out - absorbed / self.v_g)
