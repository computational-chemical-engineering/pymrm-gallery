"""Henze et al. (1987): the ASM1 Petersen matrix and a plug-flow lane (page J4.8).

Henze, Grady, Gujer, Marais and Matsuo, Water Research 21(5), 505-515 (1987).
Two parts:

- the stoichiometric matrix held symbolically, with the composition vectors
  the component definitions fix (theoretical oxygen demand, nitrogen, ionic
  charge) and the continuity residual of each process
  (:func:`petersen_matrix`, :func:`exact_factor_matrix`,
  :func:`composition_vectors`, :func:`continuity`, :func:`recover_cell`,
  :func:`numeric_matrix`);
- the eight rate expressions of Table 1 (:func:`rates`) on a steady plug-flow
  lane with an unaerated head end and internal mixed-liquor recycle
  (:class:`Lane`), solved as one Newton system with the recycle inside it, and
  the same lane by a fixed-point iteration on the inlet value as an
  independent route (:meth:`Lane.solve_picard`).

Every defect the page's break table injects enters through an argument of
these functions. No data are loaded here: the page passes the tables in.
"""
import numpy as np
import sympy as sp
from pymrm import construct_div, construct_convflux_upwind, NumJac, newton

STOICH_SYMBOLS = ("Y_H", "Y_A", "f_P", "i_XB", "i_XP")
SYMS = {s: sp.Symbol(s) for s in STOICH_SYMBOLS}
ORGANIC = ("S_I", "S_S", "X_I", "X_S", "X_BH", "X_BA", "X_P")
TINY = 1e-10          # guards the X_S / X_BH and X_ND / X_S ratios at zero solids


# --------------------------------------------------------------- the matrix
def cell_expr(x):
    """One matrix cell as an exact sympy expression.

    Parameters
    ----------
    x : str or float
        The transcribed cell; a blank cell (``nan``) is exactly zero.

    Returns
    -------
    sympy.Expr
        ``rational=True`` keeps 2.86 as 143/50, so every balance is an exact
        rational identity and not a floating-point near-miss.
    """
    if isinstance(x, float) and np.isnan(x):
        return sp.Integer(0)
    return sp.sympify(str(x), locals=SYMS, rational=True)


def petersen_matrix(stoi, components, add_n2=True):
    """The Petersen matrix as a dictionary of exact cells.

    Parameters
    ----------
    stoi : pandas.DataFrame
        Table 1 as strings, one row per process with column ``j`` and one
        column per component.
    components : list[str]
        The paper's components in order; ``S_N2`` is appended.
    add_n2 : bool, optional
        Give process 2 the dinitrogen cell ``-nu_{2,NO}``. ``False`` keeps the
        ``S_N2`` column but leaves it blank, as the paper prints it.

    Returns
    -------
    dict
        ``nu[j][component]`` for processes ``j = 1..8``.
    """
    nu = {int(r.j): {c: cell_expr(getattr(r, c)) for c in components}
          for r in stoi.itertuples()}
    for j in nu:
        nu[j]["S_N2"] = -nu[2]["S_NO"] if (j == 2 and add_n2) else sp.Integer(0)
    return nu


def exact_factor_matrix(nu):
    """The matrix with 20/7 and 32/7 in place of the printed 2.86 and 4.57.

    Parameters
    ----------
    nu : dict
        Matrix from :func:`petersen_matrix`.

    Returns
    -------
    dict
        A copy with the three process-2 cells and the process-3 oxygen cell
        rebuilt from the exact fractions.
    """
    Y_H, Y_A, i_XB = SYMS["Y_H"], SYMS["Y_A"], SYMS["i_XB"]
    out = {j: dict(v) for j, v in nu.items()}
    out[2]["S_NO"] = -(1 - Y_H) / (sp.Rational(20, 7) * Y_H)
    out[2]["S_N2"] = -out[2]["S_NO"] if nu[2]["S_N2"] != 0 else sp.Integer(0)
    out[2]["S_ALK"] = (1 - Y_H) / (14 * sp.Rational(20, 7) * Y_H) - i_XB / 14
    out[3]["S_O"] = -(sp.Rational(32, 7) - Y_A) / Y_A
    return out


def composition_vectors(components, cod_no=sp.Rational(-64, 14),
                        charge_nh=sp.Rational(1, 14)):
    """Content of COD, nitrogen and charge per unit of each component.

    COD is theoretical oxygen demand with nitrogen referenced to -III
    (ammonium and organic amine), the state dichromate COD leaves it in.

    Parameters
    ----------
    components : list[str]
        Component names, ``S_N2`` included.
    cod_no : sympy.Expr, optional
        COD of nitrate per g N (NO3- to NH4+ accepts 8 electrons, -64/14).
    charge_nh : sympy.Expr, optional
        Charge of ammonium in mol per g N (NH4+, +1/14).

    Returns
    -------
    dict
        ``gamma[law][component]`` for ``law`` in ``"COD"``, ``"N"``,
        ``"charge"``, exact sympy values (``i_XB`` and ``i_XP`` symbolic).
    """
    zero = {c: sp.Integer(0) for c in components}
    cod, nit, chg = dict(zero), dict(zero), dict(zero)
    for c in ORGANIC:
        cod[c] = sp.Integer(1)                        # declared in COD units
    cod["S_O"] = sp.Integer(-1)                       # "oxygen (negative COD)"
    cod["S_NO"] = cod_no
    cod["S_N2"] = sp.Rational(-24, 14)                # N2 -> NH4+ accepts 6 e- per N2
    for c in ("S_NO", "S_NH", "S_ND", "X_ND", "S_N2"):
        nit[c] = sp.Integer(1)                        # "all given as nitrogen"
    nit["X_BH"] = nit["X_BA"] = SYMS["i_XB"]          # mass N / mass COD in biomass
    nit["X_P"] = SYMS["i_XP"]                         # ... in products from biomass
    chg["S_NH"] = charge_nh                           # NH4+, per g N
    chg["S_NO"] = sp.Rational(-1, 14)                 # NO3-, per g N
    chg["S_ALK"] = sp.Integer(-1)                     # alkalinity, molar units
    return {"COD": cod, "N": nit, "charge": chg}


def continuity(nu, gamma_law):
    """Residual of one conservation law on every process; zero means it closes.

    Parameters
    ----------
    nu : dict
        Matrix from :func:`petersen_matrix`.
    gamma_law : dict
        One composition vector of :func:`composition_vectors`.

    Returns
    -------
    dict
        ``{j: simplified residual}``.
    """
    return {j: sp.simplify(sum(nu[j][c] * gamma_law[c] for c in nu[j]))
            for j in sorted(nu)}


def recover_cell(nu, gamma_law, j, cell):
    """Solve one conservation law of process ``j`` for one of its cells.

    Parameters
    ----------
    nu : dict
        Matrix from :func:`petersen_matrix`.
    gamma_law : dict
        The composition vector of the law used.
    j : int
        Process.
    cell : str
        Component whose cell is deleted and recovered. For the process-2
        nitrate cell the dinitrogen cell follows it, as nitrogen continuity
        ties them.

    Returns
    -------
    sympy.Expr
        The recovered cell minus the printed one (zero when identical).
    """
    x = sp.Symbol("x")
    trial = dict(nu[j])
    trial[cell] = x
    if j == 2 and cell == "S_NO" and nu[2]["S_N2"] != 0:
        trial["S_N2"] = -x
    sol = sp.solve(sum(trial[c] * gamma_law[c] for c in trial), x)[0]
    return sp.simplify(sol - nu[j][cell])


def numeric_matrix(nu, values):
    """The matrix at numerical stoichiometric parameters.

    Parameters
    ----------
    nu : dict
        Matrix from :func:`petersen_matrix`.
    values : dict
        Floats for ``Y_H``, ``Y_A``, ``f_P``, ``i_XB``, ``i_XP``.

    Returns
    -------
    numpy.ndarray
        Shape ``(8, n_components)``.
    """
    sub = {SYMS[k]: sp.Rational(str(values[k])) for k in STOICH_SYMBOLS}
    return np.array([[float(nu[j][c].subs(sub)) for c in nu[j]] for j in sorted(nu)])


# ---------------------------------------------------------------- the rates
def rates(c, p, ix):
    """Process rates of Table 1, right-hand column.

    Parameters
    ----------
    c : numpy.ndarray
        Concentrations, shape ``(..., n_components)``, g m-3 (alkalinity in
        mol m-3); negative values are read as zero.
    p : dict
        Kinetic parameters by the paper's symbols. An extra ``K_ALK`` (mol
        m-3) multiplies the nitrification rate by ``S_ALK / (K_ALK + S_ALK)``;
        ASM1 has no such term.
    ix : dict
        Component name to index on the last axis.

    Returns
    -------
    numpy.ndarray
        ``rho_j``, shape ``(..., 8)``, g COD m-3 d-1 (process 6 in g N).
    """
    def q(name):
        return np.maximum(c[..., ix[name]], 0.0)
    S_S, S_O, S_NO, S_NH, S_ND = q("S_S"), q("S_O"), q("S_NO"), q("S_NH"), q("S_ND")
    X_S, X_BH, X_BA, X_ND = q("X_S"), q("X_BH"), q("X_BA"), q("X_ND")
    ox = S_O / (p["K_OH"] + S_O)                      # oxygen switch, eq. (4)
    anox = p["K_OH"] / (p["K_OH"] + S_O)              # its complement, eq. (5)
    nit = S_NO / (p["K_NO"] + S_NO)
    sub = S_S / (p["K_S"] + S_S)
    r = np.empty(c.shape[:-1] + (8,))
    r[..., 0] = p["mu_hat_H"] * sub * ox * X_BH
    r[..., 1] = p["mu_hat_H"] * sub * anox * nit * p["eta_g"] * X_BH
    r[..., 2] = p["mu_hat_A"] * S_NH / (p["K_NH"] + S_NH) * S_O / (p["K_OA"] + S_O) * X_BA
    if "K_ALK" in p:      # not in ASM1: an alkalinity switch, used only as a defect
        S_ALK = q("S_ALK")
        r[..., 2] *= S_ALK / (p["K_ALK"] + S_ALK)
    r[..., 3] = p["b_H"] * X_BH
    r[..., 4] = p["b_A"] * X_BA
    r[..., 5] = p["k_a"] * S_ND * X_BH
    ratio = X_S / np.maximum(X_BH, TINY)
    r[..., 6] = p["k_h"] * ratio / (p["K_X"] + ratio) * (ox + p["eta_h"] * anox * nit) * X_BH
    r[..., 7] = r[..., 6] * X_ND / np.maximum(X_S, TINY)
    return r


# ----------------------------------------------------------------- the lane
def faces(n, f_anox):
    """Cell faces on the unit lane, with one face on the aeration boundary.

    Letting the anoxic/aerobic step fall inside a cell makes the
    discretisation error jump with ``n``; the operators take any face array.

    Parameters
    ----------
    n : int
        Number of cells.
    f_anox : float
        Unaerated fraction of the lane.

    Returns
    -------
    numpy.ndarray
        The ``n + 1`` faces.
    """
    if not 0.0 < f_anox < 1.0:
        return np.linspace(0.0, 1.0, n + 1)
    n1 = min(max(int(round(f_anox * n)), 1), n - 1)
    return np.concatenate([np.linspace(0.0, f_anox, n1 + 1),
                           np.linspace(f_anox, 1.0, n - n1 + 1)[1:]])


class Lane:
    """Steady plug-flow aeration lane with an unaerated head end and recycle.

    ``d(u c_i)/dz = sum_j nu_ji rho_j(c) + delta_{i,S_O} kLa (S_O* - S_O)`` on
    the unit lane, with the inlet the blend of feed and recycled outlet,
    ``c_in = (c_feed + R c_out) / (1 + R)``.

    Parameters
    ----------
    nu_m : numpy.ndarray
        Numerical matrix, shape ``(8, n_c)``.
    p : dict
        Kinetic parameters.
    ix : dict
        Component name to index.
    c_feed : numpy.ndarray
        Blended feed, shape ``(n_c,)``.
    n : int, optional
        Cells.
    f_anox : float, optional
        Unaerated fraction.
    R : float, optional
        Internal recycle ratio.
    tau : float, optional
        Lane volume over (influent + return sludge) flow, d.
    so_set : float, optional
        Dissolved-oxygen setpoint of the aerated section, g m-3.
    kla : float, optional
        Transfer coefficient, d-1, large enough to hold the setpoint.
    flux_read : {"face", "upstream", "last_cell"}, optional
        Outlet flux of the balances from the operator's exit face, from the
        last-but-one cell (an off-by-one defect for the break table), or from
        the last cell, ``u c[-1]``, which differs from the exit face under a
        zero-gradient outlet.
    kla_per_cell : bool, optional
        Multiply the transfer coefficient by each cell's width, so the oxygen
        transfer is per cell rather than per volume (a defect for the break
        table).
    outlet : {"outflow", "zero_gradient"}, optional
        Exit condition: ``{"outflow": True}`` (the exit face carries the last
        cell) or ``dc/dz = 0`` at z = 1, whose exit face value the upwind
        operator extrapolates from the interior. The recycle reads the exit
        face in both cases.
    """

    def __init__(self, nu_m, p, ix, c_feed, n=150, f_anox=0.30, R=2.0, tau=8.0 / 24.0,
                 so_set=2.0, kla=2000.0, flux_read="face", kla_per_cell=False,
                 outlet="outflow"):
        self.nu_m, self.p, self.ix, self.R, self.so_set = nu_m, p, ix, R, so_set
        self.c_feed = np.asarray(c_feed, dtype=float)
        nc = nu_m.shape[1]
        self.n, self.nc, self.shape = n, nc, (n, nc)
        self.z_f = faces(n, f_anox)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.dz = np.diff(self.z_f)
        self.kla = np.where(self.z_c < f_anox, 0.0, kla)
        if kla_per_cell:               # a defect: transfer scaled by the cell width
            self.kla = self.kla * self.dz
        self.u = (1.0 + R) / tau                      # velocity on the unit lane, d-1
        self.flux_read = flux_read
        if outlet == "outflow":
            exit_bc = {"outflow": True}               # z=1: exit face carries the last cell
        elif outlet == "zero_gradient":
            exit_bc = {"a": 1.0, "b": 0.0, "d": 0.0}  # z=1: dc/dz = 0 (outward normal +z)
        else:
            raise ValueError(f"unknown outlet {outlet!r}")
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0},         # z=0: c = 1 * c_in, c_in external
              exit_bc)
        div = construct_div(self.shape, self.z_f, nu=0)        # nu=0: Cartesian, constant area
        conv, conv_l, conv_r = construct_convflux_upwind(
            self.shape, self.z_f, self.z_c, bc, v=self.u, shapes_d=((1, nc), None))
        self.conv, self.conv_l = conv, conv_l
        self.conv_r = np.asarray(conv_r.todense() if hasattr(conv_r, "todense")
                                 else conv_r, dtype=float).reshape(-1)
        self.A = (div @ conv).tocsc()
        self.B = (div @ conv_l).tocsc()               # inlet value -> cell balances
        self.b_r = np.asarray(div @ self.conv_r).reshape(-1)
        # outlet face value as a matrix on the state: the recycle reads what the
        # operator transports through the exit face
        last = slice(n * nc, (n + 1) * nc)
        self.E = (conv[last, :] / self.u).tocsc()
        self.e_r = self.conv_r[last] / self.u
        self.K = (self.B @ self.E * (R / (1.0 + R))).tocsc()   # recycle coupling, constant
        self.numjac = NumJac(self.shape)              # pointwise source: last axis only

    def source(self, c):
        """Reaction plus oxygen transfer, shape ``(n, n_c)``."""
        s = rates(c, self.p, self.ix) @ self.nu_m
        s[..., self.ix["S_O"]] += self.kla * (self.so_set - c[..., self.ix["S_O"]])
        return s

    def inlet(self, x):
        """Lane inlet value from the state: feed blended with the outlet."""
        c_out = self.E @ x.reshape(-1) + self.e_r
        return (self.c_feed + self.R * c_out) / (1.0 + self.R)

    def _residual(self, x, c_in=None):
        x = x.reshape(self.shape)
        g_s, jac_s = self.numjac(self.source, x)
        if c_in is None:                              # recycle inside the system
            b = self.B @ ((self.c_feed + self.R * self.e_r) / (1.0 + self.R)) + self.b_r
            g = (self.A + self.K) @ x.reshape(-1) + b - g_s.reshape(-1)
            return g, self.A + self.K - jac_s
        g = self.A @ x.reshape(-1) + self.B @ c_in + self.b_r - g_s.reshape(-1)
        return g, self.A - jac_s

    def _newton(self, fun, c0, tol, maxfev):
        # concentrations are bounded below by zero: the clipped rates of Table 1
        # give no restoring term below it (pitfalls P9)
        sol = newton(fun, c0, maxfev=maxfev, tol=tol, stagnation=True, bounds=(0.0, None))
        if not sol.success:
            raise RuntimeError(f"lane newton failed (n={self.n}): {sol.message}")
        return sol

    def solve(self, tol=1e-11, maxfev=80):
        """Steady state with the recycle in the Newton system.

        Parameters
        ----------
        tol : float, optional
            Newton step tolerance, g m-3.
        maxfev : int, optional
            Newton iteration cap.

        Returns
        -------
        dict
            See :meth:`_result`.
        """
        c0 = np.tile(self.c_feed, (self.n, 1))
        sol = self._newton(self._residual, c0, tol, maxfev)
        return self._result(sol.x.reshape(self.shape), nit=sol.nit)

    def solve_picard(self, tol=1e-9, maxouter=300, newton_tol=1e-11):
        """The same steady state by a fixed-point iteration on the inlet value.

        Each outer iterate solves the lane for a fixed inlet, then sets the
        inlet to the blend of feed and the new outlet; a passive component
        contracts at ``R / (1 + R)`` per iterate.

        Parameters
        ----------
        tol : float, optional
            Stop when the inlet moves less than this, g m-3.
        maxouter : int, optional
            Iteration cap; reaching it raises.
        newton_tol : float, optional
            Step tolerance of each inner Newton solve.

        Returns
        -------
        dict
            See :meth:`_result`, plus ``outer`` (iterations) and
            ``contraction`` (geometric-mean ratio of successive inlet updates
            between 1e-6 and 1e-3 g m-3).
        """
        c = np.tile(self.c_feed, (self.n, 1))
        c_in, steps = self.c_feed.copy(), []
        for it in range(maxouter):
            sol = self._newton(lambda x: self._residual(x, c_in), c, newton_tol, 60)
            c = sol.x.reshape(self.shape)
            target = self.inlet(c)
            steps.append(np.max(np.abs(target - c_in)))
            c_in = target
            if steps[-1] < tol:
                out = self._result(c, c_in=c_in)
                # contraction per iterate, from the updates between 1e-6 and 1e-3,
                # well above the inner solve's tolerance
                k = [i for i, st in enumerate(steps) if 1e-6 <= st <= 1e-3]
                if len(k) < 2:
                    raise RuntimeError("too few fixed-point updates to measure the contraction")
                q = (steps[k[-1]] / steps[k[0]]) ** (1.0 / (k[-1] - k[0]))
                out.update(outer=it + 1, contraction=q)
                return out
        raise RuntimeError(f"recycle fixed point not reached in {maxouter} iterations")

    def _result(self, c, c_in=None, nit=None):
        """Profiles, face fluxes and oxygen transfer of a converged state."""
        c_in = self.inlet(c) if c_in is None else c_in
        flux = np.asarray(self.conv @ c.reshape(-1) + self.conv_l @ c_in
                          + self.conv_r).reshape(self.n + 1, self.nc)
        if self.flux_read == "upstream":
            flux_out = self.u * c[-2]
        elif self.flux_read == "last_cell":
            flux_out = self.u * c[-1]
        else:
            flux_out = flux[-1]
        otr = float(np.sum(self.kla * (self.so_set - c[:, self.ix["S_O"]]) * self.dz))
        return dict(c=c, c_in=c_in, c_out=self.E @ c.reshape(-1) + self.e_r, z=self.z_c,
                    dz=self.dz, n=self.n, otr=otr, nit=nit, flux_in=flux[0],
                    flux_out=flux_out, feed=self.c_feed)


def balance(state, w, cod=False):
    """Steady balance of a conserved quantity over the lane.

    Parameters
    ----------
    state : dict
        A :class:`Lane` result.
    w : numpy.ndarray
        Composition vector of the quantity, shape ``(n_c,)``.
    cod : bool, optional
        Subtract the oxygen transferred in (``gamma_COD(S_O) = -1``).

    Returns
    -------
    float
        Inflow plus boundary source minus outflow; zero when it is conserved.
    """
    src = -state["otr"] if cod else 0.0
    return float(state["flux_in"] @ w) + src - float(state["flux_out"] @ w)
