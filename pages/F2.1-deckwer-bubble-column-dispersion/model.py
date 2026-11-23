"""Liquid axial dispersion with interphase oxygen transfer in a tall bubble column.

Deckwer, Burckhart and Zoll, Chem. Eng. Sci. 29, 2177-2188 (1974). Units as in the
paper (cgs): lengths in cm, velocities in cm/s, D in cm2/s, kLa in 1/s, c in mg/l,
p in atm. The models work in x = z/L, theta = c / c_ref and psi = F / (y_ref n_G0).

Closed forms
    ``adm_profile``: the liquid model (eq. 3, Danckwerts inlet, closed outlet) with an
    equilibrium concentration that is a polynomial of degree two in x (eqs. A6 to A20,
    with the misprints of A4 and A19 corrected; ``inlet_sign`` and ``b3_prefactor``
    reproduce the printed forms).
Coupled model (not solved in the paper)
    ``CoupledColumn``: the liquid model with the gas in plug flow, its oxygen flux
    following the transfer, on pymrm finite-volume operators; ``coupled_bvp`` is the
    same model by ``scipy.integrate.solve_bvp``, an independent route.

No plotting and no data loading here.
"""
import numpy as np
from scipy import sparse
from scipy.integrate import solve_bvp

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind,
                   construct_div, construct_grad, newton)
from pymrm.checks import find_roots

ATM_DYN_CM2 = 1.01325e6  # dyn/cm2 per atm


# ---------------------------------------------------------------- groups (eqs. 1, 4, 5, 8)
def liquid_velocity(VL_m3h, d_cm):
    """Superficial liquid velocity u_L0 in cm/s from a flow in m3/h."""
    return VL_m3h * 1e6 / 3600.0 / (np.pi * d_cm**2 / 4.0)


# codehealth: allow params -- eq. 1 with its printed constants as defaults
def dispersion_eq1(d_cm, uG, coef=2.7, exp_d=1.4, exp_uG=0.3):
    """Liquid axial dispersion coefficient D (cm2/s), eq. 1."""
    return coef * d_cm**exp_d * uG**exp_uG


def kla_eq8(uG, b0, b1):
    """Volumetric mass transfer coefficient kLa (1/s), eq. 8 with Table 1 constants."""
    return b0 * uG**b1


def peclet(uL, L, D, eps_G):
    """Pe = u_L0 L / (D (1 - eps_G)), eq. 4."""
    return uL * L / (D * (1.0 - eps_G))


def stanton(kla, L, uL):
    """St = kLa L / u_L0, eq. 5."""
    return kla * L / uL


def head_atm(eps_G, L, rho, g):
    """Hydrostatic head of the dispersion, rho g (1 - eps_G) L, in atm."""
    return rho * g * (1.0 - eps_G) * L / ATM_DYN_CM2


# codehealth: allow params -- the dimensional inputs of one group, each a physical quantity
def capacity_ratio(uL, c_ref, M, y_ref, uG0, p0, R, T):
    """Lambda = u_L0 c_ref / (M y_ref n_G0): liquid oxygen capacity over gas oxygen flow.

    ``c_ref`` in mg/l (1e-6 g/cm3), ``M`` in g/mol, ``n_G0 = uG0 p0 / (R T)`` in
    mol/(cm2 s) with ``R`` in cm3 atm/(mol K).
    """
    n_G0 = uG0 * p0 / (R * T)
    return uL * c_ref * 1e-6 / (M * y_ref * n_G0)


def roots_A7(Pe, St):
    """Roots r1 > 0 > r2 of r**2 / Pe - r - St = 0 (eq. A7)."""
    q = np.sqrt(1.0 + 4.0 * St / Pe)
    return 0.5 * Pe * (1.0 + q), 0.5 * Pe * (1.0 - q)


# codehealth: allow params -- the gas composition, the head and the printed-sign switch
def cstar_coefficients(y0, dy, y_ref, head_ratio, bprime_sign=-1.0):
    """theta* = alpha + beta x + gamma x**2 for y = y0 + dy x and p = p_t (1 + h (1 - x)).

    In units of c_ref = y_ref p_t / H (eqs. 6, 7, A12, A13 with b' = -rho g eps_L L / H).
    ``bprime_sign = +1`` is A13's b' as printed (pressure rising with height).
    """
    s0, s1 = y0 / y_ref, dy / y_ref
    hb = bprime_sign * head_ratio
    return s0 * (1.0 + head_ratio), s1 * (1.0 + head_ratio) + s0 * hb, s1 * hb


# ---------------------------------------------------------------- closed forms
# codehealth: allow params -- the published closed form's groups plus the two printed-form switches
def adm_profile(x, Pe, St, alpha, beta=0.0, gamma=0.0, theta0=0.0, inlet_sign=1.0,
                b3_prefactor="St"):
    """Liquid profile theta(x) of (1/Pe) theta'' - theta' + St (theta* - theta) = 0.

    theta* = alpha + beta x + gamma x**2. Inlet (Danckwerts)
    theta(0) - inlet_sign (1/Pe) theta'(0) = theta0; outlet theta'(1) = 0.
    ``inlet_sign = -1`` is eq. A4 as printed. ``b3_prefactor = "Pe"`` is eq. A19 as
    printed (2 gamma / Pe instead of 2 gamma / St). Accepts complex ``St`` (complex-step
    derivatives). The growing exponential is written as exp(r1 (x - 1)).
    """
    x = np.asarray(x, dtype=float)
    r1, r2 = roots_A7(Pe, St)
    b5 = gamma
    b4 = beta - 2.0 * gamma / St
    pre = 2.0 * gamma / (St if b3_prefactor == "St" else Pe)
    b3 = pre * (1.0 / Pe + 1.0 / St) + alpha - beta / St
    s = inlet_sign / Pe
    # unknowns A (of exp(r1 (x - 1))) and B (of exp(r2 x))
    a11, a12 = np.exp(-r1) * (1.0 - s * r1), 1.0 - s * r2
    a21, a22 = r1, r2 * np.exp(r2)
    f1 = theta0 - (b3 - s * b4)
    f2 = -(b4 + 2.0 * b5)
    det = a11 * a22 - a12 * a21
    A = (f1 * a22 - a12 * f2) / det
    B = (a11 * f2 - a21 * f1) / det
    return A * np.exp(r1 * (x - 1.0)) + B * np.exp(r2 * x) + b3 + b4 * x + b5 * x**2


# codehealth: allow params -- the printed closed form's own arguments
def adm_printed_A6_A11(x, Pe, St, a, b, c0):
    """Eqs. A6 and A8 to A11 as printed, for c* = a + b x (any consistent unit)."""
    r1, r2 = roots_A7(Pe, St)
    Bp = (a - c0 - b / St) * Pe - b
    N = r1**2 * np.exp(r1) - r2**2 * np.exp(r2)
    A1 = (Bp * r2 * np.exp(r2) - b * r1) / N
    A2 = (-Bp * r1 * np.exp(r1) + b * r2) / N
    return A1 * np.exp(r1 * x) + A2 * np.exp(r2 * x) + a + b * x - b / St


def danckwerts_outlet_ratio(Pe, St):
    """c(1)/c0 of the closed vessel with c* = 0 (Danckwerts 1953; page A2.1)."""
    q = np.sqrt(1.0 + 4.0 * St / Pe)
    return 4.0 * q * np.exp(0.5 * Pe * (1.0 - q)) / (
        (1.0 + q)**2 - (1.0 - q)**2 * np.exp(-q * Pe))


# codehealth: allow params -- the closed form's groups
def plug_profile(x, St, alpha, beta, theta0):
    """Plug flow, theta' = St (alpha + beta x - theta), theta(0) = theta0."""
    x = np.asarray(x, dtype=float)
    return alpha + beta * x - beta / St + (theta0 - alpha + beta / St) * np.exp(-St * x)


def cstr_outlet(St, alpha, beta, theta0):
    """Stirred tank with the column-mean equilibrium alpha + beta / 2."""
    return (theta0 + St * (alpha + 0.5 * beta)) / (1.0 + St)


def st_two_point_cstr(theta_out, theta0, alpha, beta):
    """St a stirred-tank evaluation infers from feed and outlet (p. 2180)."""
    return (theta_out - theta0) / (alpha + 0.5 * beta - theta_out)


def st_two_point_ntu(theta_out, theta0):
    """St of the plug-flow NTU evaluation for c* = 0: ln(c0 / c(1)) (p. 2184)."""
    return np.log(theta0 / theta_out)


# ---------------------------------------------------------------- the fit (eq. 15 of the plan)
# codehealth: allow params -- the model, the design, the data and the bracket
def fit_stanton(model, x, data, St_guess, n_scan=24):
    """St minimising sum (model(x, St) - data)**2, Pe and theta0 held in ``model``.

    The minimum is the root of the analytic gradient (complex-step derivative of
    the model), found by ``find_roots`` on [St_guess / 3, 3 St_guess]. Returns
    ``nan`` unless exactly one root with a positive curvature lies in that bracket
    (the fit is then not identifiable there).
    """
    data = np.asarray(data, dtype=float)

    def grad(S):
        h = 1e-20 * S
        m = model(x, S + 1j * h)
        return float(np.sum((m.real - data) * m.imag / h))

    lo, hi = St_guess / 3.0, 3.0 * St_guess
    res = find_roots(grad, lo, hi, n_scan=n_scan, log=True, xtol=1e-14 * St_guess)
    mins = [r for r in res["roots"] if grad(r * (1 + 1e-6)) > grad(r * (1 - 1e-6))]
    return mins[0] if len(mins) == 1 else np.nan


# ---------------------------------------------------------------- gas phase
def y_over_yref(psi, y0, y_ref):
    """y / y_ref from psi = F / (y_ref n_G0), F = G y / (1 - y), G = (1 - y0) n_G0."""
    return psi / ((1.0 - y0) + y_ref * psi)


# codehealth: allow params -- the overall balance's groups and the balance switch
def dy_from_outlet(theta_out, theta0, Lam, psi0, y0, y_ref, balance="inert"):
    """Total change of the gas mole fraction from the overall oxygen balance.

    ``"inert"``: psi(1) = psi0 - Lam (theta(1) - theta0), inert gas conserved (E11).
    ``"naive"``: constant total molar flow, dy = -y_ref Lam (theta(1) - theta0) (E12).
    """
    if balance == "naive":
        return -y_ref * Lam * (theta_out - theta0)
    psi1 = psi0 - Lam * (theta_out - theta0)
    return y_ref * y_over_yref(psi1, y0, y_ref) - y0


class CoupledColumn:
    """Liquid axial dispersion coupled to gas plug flow, on pymrm operators.

    Liquid: (1/Pe) theta'' - theta' + St (theta* - theta) = 0, Danckwerts inlet,
    closed outlet. Gas: psi' = -Lam St (theta* - theta), psi(0) = psi0, with
    theta* = (y / y_ref) P(x), P = 1 + h (1 - x). kLa and eps_G constant along x.

    Parameters
    ----------
    Pe, St, Lam : float
        Peclet number, Stanton number, capacity ratio (``capacity_ratio``).
    theta0, psi0 : float
        Liquid feed and gas oxygen flux at x = 0.
    y0, y_ref : float
        Gas mole fraction at the entry and the reference of c_ref.
    head_ratio : float
        h = rho g eps_L L / p_t.
    n : int
        Number of finite volumes (uniform).
    inlet_sign : float
        +1 the Danckwerts inlet; -1 eq. A4 as printed (a break row).
    gas_outlet : {"outflow", "zero_gradient"}
        Gas outlet condition; ``"zero_gradient"`` is a break row (pitfall P6).
    outlet_read : {"face", "cell"}
        Outlet values from ``compute_boundary_values`` or the last cell (break row, P8).
    bprime_sign : float
        -1 the hydrostatic pressure falling with height; +1 A13's b' as printed (break row).
    """

    # codehealth: allow params -- the model's dimensionless groups and the break-row switches
    def __init__(self, Pe, St, Lam, theta0, psi0, y0, y_ref, head_ratio, n=320,
                 inlet_sign=1.0, gas_outlet="outflow", outlet_read="face", bprime_sign=-1.0):
        self.Pe, self.St, self.Lam = Pe, St, Lam
        self.theta0, self.psi0, self.y0, self.y_ref, self.h = theta0, psi0, y0, y_ref, head_ratio
        self.n, self.outlet_read, self.bprime_sign = n, outlet_read, bprime_sign
        self.x_f = np.linspace(0.0, 1.0, n + 1)
        self.x_c = 0.5 * (self.x_f[1:] + self.x_f[:-1])
        self.P_c = 1.0 + head_ratio * (1.0 + bprime_sign * self.x_c)  # p / p_t
        # Liquid, x = 0, outward normal n = -x: theta - (1/Pe) dtheta/dx = theta0 is
        # (1/Pe) dtheta/dn + theta = theta0. x = 1, n = +x: dtheta/dx = 0.
        self.bc_l = ({"a": inlet_sign / Pe, "b": 1.0, "d": theta0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        # Gas, x = 0: psi = psi0 (Dirichlet). x = 1: convective outflow, the face carries
        # the last cell (P6); the zero-gradient variant extrapolates the face value.
        out = {"outflow": True} if gas_outlet == "outflow" else {"a": 1.0, "b": 0.0, "d": 0.0}
        self.bc_g = ({"a": 0.0, "b": 1.0, "d": psi0}, out)
        self._build_operators()

    def _build_operators(self):
        """Assemble the constant linear part once, in the (n, 2) layout [theta, psi]."""
        shape = (self.n, 1)  # each field on its own operators: outflow marks a whole boundary
        div = construct_div(shape, self.x_f, nu=0)  # nu = 0: Cartesian (constant cross-section)
        grad_l, grad_l_bc = construct_grad(shape, self.x_f, self.x_c, self.bc_l)
        conv_l, conv_l_bc = construct_convflux_upwind(shape, self.x_f, self.x_c, self.bc_l, v=1.0)
        conv_g, conv_g_bc = construct_convflux_upwind(shape, self.x_f, self.x_c, self.bc_g, v=1.0)
        lin_l = div @ (conv_l - grad_l / self.Pe)
        lin_g = div @ conv_g
        e00 = sparse.diags_array([1.0, 0.0], format="csc")  # selects field 0 (liquid)
        e11 = sparse.diags_array([0.0, 1.0], format="csc")  # selects field 1 (gas)
        self.lin = (sparse.kron(lin_l, e00, format="csc") + sparse.kron(lin_g, e11, format="csc"))
        g_l = div @ (_dense(conv_l_bc) - _dense(grad_l_bc) / self.Pe)
        g_g = div @ _dense(conv_g_bc)
        self.g_const = np.column_stack([np.ravel(g_l), np.ravel(g_g)]).ravel()
        self.numjac = NumJac((self.n, 2))  # pointwise source: couples the two fields of a cell

    def theta_star(self, psi):
        """Equilibrium liquid concentration theta* = (y / y_ref) P at the cell centres."""
        return y_over_yref(psi, self.y0, self.y_ref) * self.P_c

    def _source(self, c):
        drive = self.St * (self.theta_star(c[:, 1]) - c[:, 0])
        return np.column_stack([-drive, self.Lam * drive])

    def residual(self, c_flat):
        """Residual and Jacobian of the discrete system, for ``newton``."""
        c = c_flat.reshape(self.n, 2)
        src, jac_src = self.numjac(self._source, c)
        return self.lin @ c_flat.ravel() + self.g_const + src.ravel(), self.lin + jac_src

    def solve(self, tol=1e-11, maxfev=30, lin_solver_kwargs=None):
        """Newton solve from the constant-gas closed form; raises unless ``success``."""
        al, be, ga = cstar_coefficients(self.y0, 0.0, self.y_ref, self.h, self.bprime_sign)
        th = adm_profile(self.x_c, self.Pe, self.St, al, be, ga, self.theta0)
        guess = np.column_stack([th, np.full(self.n, self.psi0)])
        sol = newton(self.residual, guess.ravel(), tol=tol, maxfev=maxfev, stagnation=True,
                     lin_solver_kwargs=lin_solver_kwargs)
        if not sol.success:
            raise RuntimeError(f"coupled Newton solve failed: {sol.message}")
        self.c = np.asarray(sol.x).reshape(self.n, 2)
        self.nit = int(sol.nit)
        return self

    def boundary(self):
        """theta(0), theta(1), psi(1) on the faces (or the last cell, a break row)."""
        th0f, _, th1f, _ = compute_boundary_values(self.c[:, :1], self.x_f, self.x_c, self.bc_l)
        _, _, psi1f, _ = compute_boundary_values(self.c[:, 1:], self.x_f, self.x_c, self.bc_g)
        th0f, th1f, psi1f = (float(np.ravel(v)[0]) for v in (th0f, th1f, psi1f))
        if self.outlet_read == "cell":
            th1f, psi1f = self.c[-1, 0], self.c[-1, 1]
        return th0f, th1f, psi1f

    def at(self, x):
        """theta and psi at interior points ``x``, by linear interpolation of cell values.

        The page reads them on faces of the grid (x = k/20 with n a multiple of 20), where
        this is the mean of the two adjacent cells, second order on a uniform grid.
        """
        return np.interp(x, self.x_c, self.c[:, 0]), np.interp(x, self.x_c, self.c[:, 1])


def _dense(m):
    """A boundary vector from a pymrm constructor as a dense column."""
    return np.asarray(m.todense()) if sparse.issparse(m) else np.asarray(m)


# codehealth: allow params -- the same groups as CoupledColumn and the solver tolerance
def coupled_bvp(Pe, St, Lam, theta0, psi0, y0, y_ref, head_ratio, tol=1e-10, n0=401):
    """The coupled model by ``scipy.integrate.solve_bvp`` (collocation, not finite volumes).

    State [theta, theta', psi]; returns the ``solve_bvp`` result, raising unless
    ``status == 0``.
    """
    def rhs(x, Y):
        ts = y_over_yref(Y[2], y0, y_ref) * (1.0 + head_ratio * (1.0 - x))
        drive = St * (ts - Y[0])
        return np.vstack([Y[1], Pe * (Y[1] - drive), -Lam * drive])

    def bc(Ya, Yb):
        return np.array([Ya[0] - Ya[1] / Pe - theta0, Yb[1], Ya[2] - psi0])

    xm = np.linspace(0.0, 1.0, n0)
    al, be, ga = cstar_coefficients(y0, 0.0, y_ref, head_ratio)
    th = adm_profile(xm, Pe, St, al, be, ga, theta0)
    dth = np.gradient(th, xm)
    sol = solve_bvp(rhs, bc, xm, np.vstack([th, dth, np.full_like(xm, psi0)]), tol=tol,
                    max_nodes=200000)
    if sol.status != 0:
        raise RuntimeError(f"solve_bvp failed: {sol.message}")
    return sol
