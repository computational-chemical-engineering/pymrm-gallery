"""Maretto and Krishna's (1999) slurry bubble column for Fischer-Tropsch synthesis.

Hydrodynamics (equations 6 to 11), mass transfer (equations 13 and 14),
Yates-Satterfield kinetics (equations 1 to 3) and the reactor (equations 15 to
21): large bubbles in plug flow with a contracting gas, small bubbles and slurry
well mixed.

Two routes solve the reactor. ``SlurryColumn`` discretises the large-bubble
balances with pymrm finite volumes and solves gas profiles and liquid
concentrations together with ``pymrm.newton`` and an exact Jacobian.
``reactor_ode`` integrates the same balances as ODEs with ``scipy`` at a tight
tolerance and solves the two liquid unknowns with ``scipy.optimize.root``; it
shares the closures with the first route and nothing of its discretisation or
solver.

Dimensionless unknowns: ``y = c_G/c_G^i`` for the large-bubble gas (shape
``(n, 2)``, species CO, H2 on the last axis) and ``x = m c_L/c_G^i`` for the
liquid (shape ``(2,)``), the gas concentration in equilibrium with the liquid
relative to the feed.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.sparse import bmat, bsr_array, csr_array, diags_array, kron
from pymrm import compute_boundary_values, construct_convflux_upwind, construct_div, newton

RGAS = 8.314                      # J/(mol K)
SPECIES = ("CO", "H2")

#: Table 1, section 3 and section 4 of the paper, SI units.
PARAMS = dict(
    T=513.15,                     # K, 240 C
    P=30e5,                       # Pa
    DT=7.0,                       # m, reactor diameter
    H=30.0,                       # m, dispersion height
    eps_df_ref=0.27,              # small-bubble holdup, paraffin oil, eq. 7
    v_small_ref=0.095,            # m/s, small-bubble rise velocity, eq. 8
    rho_g_ref=1.3,                # kg/m3, gas at ambient conditions, eq. 11
    rho_g=7.0,                    # kg/m3, syngas at 30 bar and 240 C
    C9=0.3,                       # eq. 9 coefficient
    dt_exp=-0.18,                 # eq. 9 diameter exponent
    dt_cap=1.0,                   # m, eq. 9 diameter dependence capped here
    m=(2.478, 2.964),             # c_G = m c_L, CO and H2 in C16H34
    D_L=(45.5e-9, 17.2e-9),       # m2/s, CO and H2 at 240 C
    D_L_ref=2e-9,                 # m2/s, eq. 13
    kla_large=0.5,                # (k_L a)_large / eps_b, eq. 12-13 [1/s]
    kla_small=1.0,                # (k_L a)_small / eps_df, eq. 14 [1/s]
    rho_p=647.0,                  # kg/m3, catalyst partial density (Table 1)
    alpha=-0.648,                 # contraction factor at 100 % conversion, eq. 21
    inert=0.05,                   # inert fraction of the feed
    h2_co=2.0,                    # feed ratio H2/CO
    a_pre=8.8533e-3,              # eq. 2 prefactor [mol/(s kg_cat bar^2)], printed 8.8533e3
    a_act=4494.41,                # eq. 2 temperature coefficient [K]
    b_pre=2.226,                  # eq. 3 [1/bar]
    b_act=-8236.0,                # eq. 3 temperature coefficient [K]
    T_kin=493.15,                 # K, reference temperature of eqs. 2 and 3
    nu=(1.0, 2.0),                # moles CO and H2 per mole of the eq. 1 rate
    u_liq=0.0,                    # m/s, liquid outflow of eqs. 19-20 (U_L^o, not in Table 1)
)


def hydrodynamics(U, eps_s, d_t=None, rho_g=None, p=None, **override):
    """Two-bubble-class holdup: equations 6 to 11.

    Parameters
    ----------
    U : float or numpy.ndarray
        Superficial gas velocity [m/s].
    eps_s : float
        Solids volume fraction of the slurry, gas-free basis [-].
    d_t : float, optional
        Column diameter [m]; default ``PARAMS["DT"]``.
    rho_g : float, optional
        Gas density [kg/m3]; default ``PARAMS["rho_g"]``.
    p : dict, optional
        Parameter set; default ``PARAMS``. Keywords override single entries
        (``C9``, ``dt_exp``, ...).

    Returns
    -------
    dict
        ``eps_df`` (eqs. 7 and 11), ``v_small`` (eq. 8), ``U_df`` (eq. 10),
        ``eps_b`` (eq. 9) and ``eps`` (eq. 6), each shaped like ``U``.

    Notes
    -----
    Eq. 9's diameter dependence is capped: ``min(d_t, dt_cap) ** dt_exp``. At
    any column wider than ``dt_cap`` the base is exactly 1, so the exponent has
    no effect there.
    """
    p = {**(PARAMS if p is None else p), **override}
    d_t = p["DT"] if d_t is None else d_t
    rho_g = p["rho_g"] if rho_g is None else rho_g
    eps_df = (p["eps_df_ref"] * (rho_g / p["rho_g_ref"]) ** 0.48
              * (1.0 - 0.7 / p["eps_df_ref"] * eps_s))                    # eqs. 7, 11
    v_small = p["v_small_ref"] + 0.8 * eps_s                              # eq. 8
    U_df = eps_df * v_small                                               # eq. 10
    excess = np.maximum(np.asarray(U, dtype=float) - U_df, 0.0)
    eps_b = (p["C9"] * np.minimum(d_t, p["dt_cap"]) ** p["dt_exp"]
             * excess ** 0.58)                                            # eq. 9
    return dict(eps_df=eps_df, v_small=v_small, U_df=U_df, eps_b=eps_b,
                eps=eps_b + eps_df * (1.0 - eps_b))                       # eq. 6


def kinetics(T, p=None):
    """Yates-Satterfield constants at temperature ``T``: equations 2 and 3.

    Parameters
    ----------
    T : float
        Temperature [K].
    p : dict, optional
        Parameter set; default ``PARAMS``.

    Returns
    -------
    tuple of float
        ``a`` [mol/(s kg_cat bar^2)] and ``b`` [1/bar].
    """
    p = PARAMS if p is None else p
    s = 1.0 / p["T_kin"] - 1.0 / T
    return p["a_pre"] * np.exp(p["a_act"] * s), p["b_pre"] * np.exp(p["b_act"] * s)


def rate(p_co, p_h2, a, b):
    """Equation 1, ``a p_H2 p_CO / (1 + b p_CO)^2`` [mol/(s kg_cat)], and its gradient.

    Parameters
    ----------
    p_co, p_h2 : float
        Partial pressures in equilibrium with the liquid [bar].
    a, b : float
        Constants from ``kinetics``.

    Returns
    -------
    tuple
        The rate and ``(dR/dp_CO, dR/dp_H2)``.
    """
    den = 1.0 + b * p_co
    r = a * p_h2 * p_co / den**2
    return r, (a * p_h2 * (1.0 - b * p_co) / den**3, a * p_co / den**2)


def block_diag2(blocks):
    """Sparse block-diagonal matrix from ``blocks`` of shape ``(n, 2, 2)``.

    Parameters
    ----------
    blocks : numpy.ndarray
        One 2 x 2 block per cell.

    Returns
    -------
    scipy.sparse.csr_array
        Shape ``(2 n, 2 n)``, rows and columns cell-major like the flattened state.
    """
    n = blocks.shape[0]
    # block row k holds one block, in block column k
    return bsr_array((blocks, np.arange(n), np.arange(n + 1)), shape=(2 * n, 2 * n)).tocsr()


def contraction(r, alpha, law="molar"):
    """Large-bubble superficial velocity relative to the inlet, eq. 21.

    Parameters
    ----------
    r : numpy.ndarray
        Syngas concentration relative to the feed, ``sum c_i / sum c_i^i``.
    alpha : float
        Contraction factor [-].
    law : {"molar", "concentration"}, optional
        ``"molar"``: the conversion in eq. 21 is the molar conversion of the
        large-bubble syngas, ``chi = 1 - (U - U_df) r / (U^i - U_df)``, which with
        eq. 21 gives ``(1 + alpha) / (1 + alpha r)``. ``"concentration"``: the
        conversion read as ``1 - r``, which gives ``1 + alpha (1 - r)``.

    Returns
    -------
    tuple of numpy.ndarray
        ``s = (U - U_df)/(U^i - U_df)`` and ``ds/dr``.
    """
    if law == "molar":
        den = 1.0 + alpha * r
        return (1.0 + alpha) / den, -alpha * (1.0 + alpha) / den**2
    if law == "concentration":
        return 1.0 + alpha * (1.0 - r), np.full_like(r, -alpha)
    raise ValueError(f"unknown contraction law {law!r}")


class Operating:
    """Everything that does not depend on the solution, at one operating point.

    Parameters
    ----------
    U : float
        Inlet superficial gas velocity [m/s].
    eps_s : float
        Solids volume fraction of the slurry [-].
    kla_factor : float, optional
        Multiplies both volumetric mass-transfer coefficients (the paper's
        sensitivity test).
    d_t : float, optional
        Column diameter for eq. 9 [m]; default ``PARAMS["DT"]``.
    p : dict, optional
        Parameter set; default ``PARAMS``. Keywords override single entries.

    Attributes
    ----------
    c_in : numpy.ndarray
        Feed concentrations of CO and H2 [mol/m3], ideal gas at ``P``, ``T``.
    Ub_in, U_df : float
        Inlet large-bubble and small-bubble superficial velocities [m/s].
    K, K_s : numpy.ndarray
        ``(k_L a)/m`` for large and small bubbles [1/s], per species.
    cat : float
        Catalyst per reactor volume, ``eps_L eps_s rho_p`` [kg/m3].
    """

    def __init__(self, U, eps_s, kla_factor=1.0, d_t=None, p=None, **override):
        p = {**(PARAMS if p is None else p), **override}
        self.p, self.U, self.eps_s = p, float(U), float(eps_s)
        self.hyd = hyd = hydrodynamics(U, eps_s, d_t=d_t, p=p)
        m, D_L = np.array(p["m"]), np.array(p["D_L"])
        root_d = np.sqrt(D_L / p["D_L_ref"])
        self.K = kla_factor * p["kla_large"] * root_d * hyd["eps_b"] / m       # eq. 13
        self.K_s = kla_factor * p["kla_small"] * root_d * hyd["eps_df"] / m    # eq. 14
        self.U_df = float(hyd["U_df"])
        self.Ub_in = self.U - self.U_df
        c_syn = p["P"] / (RGAS * p["T"]) * (1.0 - p["inert"])
        self.c_in = c_syn * np.array([1.0, p["h2_co"]]) / (1.0 + p["h2_co"])
        self.w = self.c_in / self.c_in.sum()      # r = w . y
        self.cat = (1.0 - hyd["eps"]) * eps_s * p["rho_p"]
        self.a, self.b = kinetics(p["T"], p)
        self.nu = np.array(p["nu"])
        self.to_bar = self.c_in * RGAS * p["T"] / 1e5    # p_i = x_i * to_bar_i

    def small(self, x):
        """Small-bubble concentrations from eqs. 17 and 18, relative to the feed."""
        g = self.K_s * self.p["H"]
        return (self.U_df + g * x) / (self.U_df + g), g / (self.U_df + g)

    def liquid_sink(self, x):
        """Reaction ``H nu cat R / c_in`` per species [m/s] and its 2 x 2 Jacobian in x."""
        R, (dR_co, dR_h2) = rate(x[0] * self.to_bar[0], x[1] * self.to_bar[1],
                                 self.a, self.b)
        f = self.p["H"] * self.cat * self.nu / self.c_in
        jac = np.outer(f, [dR_co * self.to_bar[0], dR_h2 * self.to_bar[1]])
        return f * R, jac, R

    def liquid_out(self, x):
        """Liquid outflow ``U_L c_L / c_in`` per species [m/s] and its diagonal in x."""
        g = self.p["u_liq"] / np.array(self.p["m"])
        return g * x, g


class SlurryColumn:
    """The reactor of eqs. 15 to 21 on pymrm finite volumes, one Newton solve.

    The large-bubble balance (eqs. 15 and 16) is written for the molar flux,
    ``d/dz[(U - U_df) c] = -k_L a (c/m - c_L)``. In the unknowns ``y`` the flux
    is ``q = s(r) y`` with ``s`` from ``contraction``, so first-order upwind with
    unit velocity on ``q`` is the conservative scheme for any velocity profile:
    the face flux is the upstream cell's ``(U - U_df) c``. The convection
    operator is therefore constant and assembled once.

    Parameters
    ----------
    op : Operating
        The operating point.
    n : int
        Cells over the height.
    law : {"molar", "concentration"}, optional
        Reading of the conversion in eq. 21, see ``contraction``.
    form : {"conservative", "advective"}, optional
        ``"advective"`` writes ``(U - U_df) dc/dz`` instead of the divergence of
        the flux (a defect: it drops ``c d(U - U_df)/dz``).
    outlet : {"outflow", "zero_gradient", "none"}, optional
        Outlet condition of the convection operator. ``"outflow"`` carries the
        last cell's flux out; ``"zero_gradient"`` imposes ``dq/dz = 0`` at the
        face; ``"none"`` passes no condition (a defect).
    """

    def __init__(self, op, n, law="molar", form="conservative", outlet="outflow"):
        self.op, self.n, self.law, self.form = op, n, law, form
        self.z_f = np.linspace(0.0, op.p["H"], n + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.dz = np.diff(self.z_f)
        # Inlet z = 0: (U - U_df) c = (U^i - U_df) c^i, i.e. q = 1 (a = 0, b = 1, d = 1).
        inlet = {"a": 0.0, "b": 1.0, "d": 1.0}
        # Outlet z = H: the gas leaves with the flux of the last cell, no condition
        # imposed ("outflow"). "zero_gradient" imposes dq/dz = 0 instead (outward
        # normal +z, so a = 1, b = 0, d = 0); "none" leaves the face unset.
        out = {"outflow": {"outflow": True},
               "zero_gradient": {"a": 1.0, "b": 0.0, "d": 0.0},
               "none": None}[outlet]
        self.bc = (inlet, out)
        shape = (n, 2)
        conv_mat, conv_bc = construct_convflux_upwind(shape, self.z_f, self.z_c, self.bc,
                                                      v=1.0, axis=0)
        div_mat = construct_div(shape, self.z_f, nu=0, axis=0)   # nu = 0: constant area
        self.jac_conv = (div_mat @ conv_mat).tocsr()
        self.g_conv = np.asarray((div_mat @ conv_bc).todense()).ravel()
        self.k_ex = np.tile(op.K / op.Ub_in, n)                  # K/(U^i - U_df) per row
        # x enters every cell's gas rows as -K/(U^i - U_df) x: a (2n, 2) constant block
        self.jac_gx = kron(csr_array(np.ones((n, 1))), csr_array(np.diag(-op.K / op.Ub_in)),
                           format="csr")
        # the large-bubble transfer integral of eqs. 19 and 20, divided by U
        self.w_int = np.repeat(self.dz, 2) * np.tile(op.K, n) / op.U
        self.jac_xg = kron(csr_array(self.dz[None, :]), csr_array(np.diag(op.K / op.U)),
                           format="csr")

    def flux(self, y):
        """``q = s(r) y`` per cell and its block-diagonal Jacobian (2 x 2 per cell)."""
        r = y @ self.op.w
        s, ds = contraction(r, self.op.p["alpha"], self.law)
        q = s[:, None] * y
        blocks = s[:, None, None] * np.eye(2) + (ds[:, None] * y)[:, :, None] * self.op.w
        return q, s, block_diag2(blocks)

    def residual(self, u):
        """Residual of the gas rows (eqs. 15, 16) and the two liquid rows (19, 20).

        Parameters
        ----------
        u : numpy.ndarray
            ``y`` flattened (cell-major) followed by ``x``; length ``2 n + 2``.

        Returns
        -------
        tuple
            Residual vector and its exact sparse Jacobian.
        """
        op, n = self.op, self.n
        y, x = u[:2 * n].reshape(n, 2), u[2 * n:]
        q, s, jac_q = self.flux(y)
        if self.form == "conservative":
            g_gas = self.jac_conv @ q.ravel() + self.g_conv
            jac_gg = self.jac_conv @ jac_q
        else:   # (U - U_df) dy/dz: the operator on y, scaled by s per cell
            dydz = self.jac_conv @ y.ravel() + self.g_conv
            g_gas = np.repeat(s, 2) * dydz
            ds = contraction(y @ op.w, op.p["alpha"], self.law)[1]
            ds_blocks = ds[:, None, None] * np.ones((2, 1)) * op.w   # ds_k w_j in row (k, i)
            jac_gg = (diags_array(np.repeat(s, 2)) @ self.jac_conv
                      + diags_array(dydz) @ block_diag2(ds_blocks))
        g_gas = g_gas + self.k_ex * (y.ravel() - np.tile(x, n))
        jac_gg = jac_gg + diags_array(self.k_ex)
        y_s, dys = op.small(x)
        sink, jac_sink, _ = op.liquid_sink(x)
        out, dout = op.liquid_out(x)
        g_liq = ((self.w_int * (y.ravel() - np.tile(x, n))).reshape(n, 2).sum(axis=0)
                 + (op.p["H"] * op.K_s * (y_s - x) - sink - out) / op.U)
        jac_xx = (np.diag(-self.w_int.reshape(n, 2).sum(axis=0)
                          + (op.p["H"] * op.K_s * (dys - 1.0) - dout) / op.U)
                  - jac_sink / op.U)
        jac = bmat([[jac_gg, self.jac_gx], [self.jac_xg, csr_array(jac_xx)]], format="csc")
        return np.concatenate([g_gas, g_liq]), jac

    def solve(self, guess=None, tol=1e-12, maxfev=50, check=True):
        """Solve with ``pymrm.newton``; returns a ``ColumnSolution``.

        Parameters
        ----------
        guess : numpy.ndarray, optional
            Start vector (``y`` flattened, then ``x``); default ``y = 1``,
            ``x = 0.05``.
        tol : float, optional
            Newton step tolerance on the order-one unknowns.
        maxfev : int, optional
            Newton iteration limit.
        check : bool, optional
            Raise when newton does not report success; False returns the last
            iterate (used to show what an unconverged solve looks like).

        Raises
        ------
        RuntimeError
            When newton does not report success and ``check`` is True.
        """
        n = self.n
        u0 = (np.concatenate([np.ones(2 * n), [0.05, 0.05]]) if guess is None
              else np.asarray(guess, dtype=float))
        res = newton(self.residual, u0, tol=tol, maxfev=maxfev)
        if check and not res.success:
            raise RuntimeError(f"newton: {res.message}")
        return ColumnSolution(self, res)


class ColumnSolution:
    """Result of ``SlurryColumn.solve``.

    Attributes
    ----------
    y : numpy.ndarray
        Large-bubble concentrations relative to the feed, ``(n, 2)``.
    x : numpy.ndarray
        Liquid concentrations as equilibrium gas relative to the feed, ``(2,)``.
    q_out : numpy.ndarray
        Outlet face flux ``(U - U_df) y / (U^i - U_df)`` per species, from
        ``compute_boundary_values``.
    nfev : int
        Newton iterations.
    resid : float
        Largest residual entry at the returned state (rows are order one).
    """

    def __init__(self, col, res):
        self.col, self.res, self.nfev = col, res, int(res.nfev)
        n = col.n
        self.u = np.asarray(res.x).ravel()
        self.y, self.x = self.u[:2 * n].reshape(n, 2), self.u[2 * n:]
        q, self.s, _ = col.flux(self.y)
        self.q = q
        if col.bc[1] is None:      # no outlet condition: nothing leaves through the face
            self.q_out = np.zeros(2)
        else:
            self.q_out = np.ravel(compute_boundary_values(q, col.z_f, col.z_c, col.bc)[2])
        self.q_last = q[-1]
        self.resid = float(np.max(np.abs(col.residual(self.u)[0])))

    def conversion(self, read="face"):
        """Syngas conversion, ``1 - (large + small bubble flow out) / flow in``.

        Parameters
        ----------
        read : {"face", "cell", "reaction"}, optional
            ``"face"``: the large-bubble exit from the outlet face.
            ``"cell"``: from the last cell value (a defect with
            ``outlet="zero_gradient"``). ``"reaction"``: from what the liquid
            takes up, ``(H nu cat R + U_L c_L) / (U c_in)``, which uses no gas
            exit value.

        Returns
        -------
        numpy.ndarray
            Conversion of CO and of H2.
        """
        op = self.col.op
        if read == "reaction":
            return (op.liquid_sink(self.x)[0] + op.liquid_out(self.x)[0]) / op.U
        q_out = self.q_out if read == "face" else self.q_last
        y_s = op.small(self.x)[0]
        return 1.0 - (op.Ub_in * q_out + op.U_df * y_s) / op.U

    def syngas_conversion(self, read="face"):
        """Conversion of CO + H2 together (feed-weighted)."""
        return float(self.col.op.w @ self.conversion(read))

    def large_bubble_conversion(self):
        """Molar conversion of the large-bubble syngas at the outlet face."""
        return float(1.0 - self.col.op.w @ self.q_out)


def reactor_ode(op, law="molar", rtol=1e-12, atol=1e-14, x0=(0.05, 0.05)):
    """The same reactor with the large bubbles integrated as ODEs.

    For given liquid concentrations, ``d q_i/dz = -K_i (y_i - x_i)/(U^i - U_df)``
    with ``q = s y`` is integrated by ``solve_ivp`` (DOP853); the transfer
    integral of eqs. 19 and 20 is then ``(U^i - U_df)(1 - q_i(H))`` exactly, and
    ``scipy.optimize.root`` solves the two liquid balances.

    Parameters
    ----------
    op : Operating
        The operating point.
    law : {"molar", "concentration"}, optional
        See ``contraction``.
    rtol, atol : float, optional
        ODE tolerances.
    x0 : tuple, optional
        Start for the liquid unknowns.

    Returns
    -------
    dict
        ``x``, ``q_out``, ``conversion`` (CO, H2), ``syngas_conversion`` and
        ``residual``, the largest liquid-balance residual (relative to ``U``).

    Raises
    ------
    RuntimeError
        When the liquid balances are not met to 1e-12.
    """
    alpha = op.p["alpha"]

    def y_of_q(q):
        if law == "molar":    # q gives chi directly: chi = 1 - w . q, s = 1 + alpha chi
            return q / (1.0 + alpha * (1.0 - op.w @ q))
        # concentration law: s = 1 + alpha (1 - w . y) and q = s y, so w.q = s (w.y)
        wq = op.w @ q
        a_, b_ = -alpha, 1.0 + alpha     # s = b_ + a_ r with r = w.y and wq = s r
        r = (-b_ + np.sqrt(b_**2 + 4 * a_ * wq)) / (2 * a_)
        return q / (b_ + a_ * r)

    def exit_flux(x):
        sol = solve_ivp(lambda z, q: -op.K * (y_of_q(q) - x) / op.Ub_in,
                        (0.0, op.p["H"]), [1.0, 1.0], method="DOP853",
                        rtol=rtol, atol=atol)
        return sol.y[:, -1]

    def balance(x):
        q_out = exit_flux(x)
        y_s = op.small(x)[0]
        return (op.Ub_in * (1.0 - q_out) + op.p["H"] * op.K_s * (y_s - x)
                - op.liquid_sink(x)[0] - op.liquid_out(x)[0]) / op.U

    sol = root(balance, np.asarray(x0, dtype=float), method="hybr",
               options={"xtol": 1e-13})
    x = sol.x
    resid = float(np.max(np.abs(balance(x))))
    if not resid < 1e-12:
        raise RuntimeError(f"reactor_ode: liquid balance residual {resid:.2e} ({sol.message})")
    q_out = exit_flux(x)
    conv = 1.0 - (op.Ub_in * q_out + op.U_df * op.small(x)[0]) / op.U
    return dict(x=x, q_out=q_out, conversion=conv, syngas_conversion=float(op.w @ conv),
                residual=resid)
