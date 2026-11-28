"""Concentration polarisation at a Pd membrane: the paper's film relations and a 2-D tube.

The film functions are the binary (H2/N2) reduction of Caravella, Barbieri and Drioli
(2009), Eqs. (1) to (15), with Sieverts' law at the membrane. ``TubeMarch`` is the
page's own extension: a laminar tube with a permeable wall, resolved in r and z.
Dimensionless groups (see the page): x the H2 mole fraction, theta = p_perm/P,
beta = Pi sqrt(P)/(c k), Da = Pi sqrt(P) R/(c D), Gz = u_mean d^2/(D L).
"""
import math

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.special import gamma, hyp1f1

from pymrm import (NumJac, compute_boundary_values, construct_convflux_upwind, construct_div,
                   construct_grad, newton)
from pymrm.checks import find_roots

LEVEQUE_LOCAL = (8.0 / 9.0) ** (1.0 / 3.0) / gamma(4.0 / 3.0)   # Sh = C Gz_z^(1/3), tube, local
LEVEQUE_AVG = 1.5 * LEVEQUE_LOCAL                               # averaged over 0..L


def unique_root(f, lo, hi, n_scan=64, log=False, xtol=1e-13):
    """The single root of ``f`` on ``[lo, hi]`` by ``pymrm.checks.find_roots``; raises otherwise."""
    r = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if len(r["roots"]) != 1 or r["failed"]:
        raise ValueError(f"expected one root on [{lo}, {hi}]: {r['message']}")
    return r["roots"][0]


def graetz_lambda0():
    """First Graetz eigenvalue of a tube at constant wall concentration.

    The eigenfunction is exp(-l xi^2/2) M(1/2 - l/4, 1, l xi^2) (Kummer), zero at the wall.
    Returns lambda_0; the developed Sherwood number is lambda_0^2 / 2.
    """
    return unique_root(lambda lam: hyp1f1(0.5 - lam / 4.0, 1.0, lam), 1.0, 4.0)


def z_bernoulli(phi, form="eq6"):
    """Eq. (6), phi/(e^phi - 1) (minus sign restored, M3), or its series Eq. (7)."""
    if form == "eq7":
        return 1.0 - phi / 2.0 + phi * phi / 12.0
    return phi / math.expm1(phi) if phi != 0.0 else 1.0


def film_residual(xs, xb, beta, theta, kind="exact"):
    """Film relation minus the Sieverts flux, in units of c k, at surface fraction ``xs``.

    kind : 'exact' (Stefan film, n = ln((1 - xs)/(1 - xb))), 'printed' (Eq. 1 with the
    surface fraction in the convective term and Eq. 6 as printed for Z, M5), 'printed_eq7'
    (the same with Eq. 7), 'endform' (Eq. 1 with Z = phi/(1 - e^-phi), identical to the
    exact film) or 'linear' (n = xb - xs, no Stefan term).
    """
    n = beta * (math.sqrt(xs) - math.sqrt(theta))
    if kind == "exact":
        return math.log((1.0 - xs) / (1.0 - xb)) - n
    if kind == "linear":
        return (xb - xs) - n
    if kind == "endform":      # Eq. (6) in the end-state form phi/(1 - e^-phi), with x_s (M5)
        z = n / -math.expm1(-n) if n != 0.0 else 1.0
        return z * (xb - xs) + n * xs - n
    z = z_bernoulli(n, "eq7" if kind == "printed_eq7" else "eq6")
    return z * (xb - xs) + n * xs - n


def film_xs(xb, beta, theta, kind="exact"):
    """Surface H2 fraction of the film in series with the membrane (one root, checked)."""
    lo, hi = sorted((xb, theta))
    return unique_root(lambda s: film_residual(s, xb, beta, theta, kind), lo + 1e-15, hi - 1e-15)


def film_cpc(xb, beta, theta, kind="exact"):
    """CPC of the film at bulk fraction ``xb`` (Eq. 13)."""
    return cpc(film_xs(xb, beta, theta, kind), xb, theta)


def cpc(xs, xb, theta):
    """Concentration polarisation coefficient, Eq. (13) with Sieverts' law (Eq. 15)."""
    return 1.0 - (math.sqrt(xs) - math.sqrt(theta)) / (math.sqrt(xb) - math.sqrt(theta))


def cpc_zero_flux(beta, theta, stefan_factor=True):
    """Closed-form CPC at x_b -> theta for the exact film: G/(1 + G).

    G = beta (1 - theta)/(2 sqrt(theta)); ``stefan_factor=False`` drops (1 - theta).
    """
    g = beta * ((1.0 - theta) if stefan_factor else 1.0) / (2.0 * math.sqrt(theta))
    return g / (1.0 + g)


def film_module(x_in, theta, gz, da, sh, kind="exact", rtol=1e-12):
    """1-D plug-flow module with a film coefficient ``sh(zeta)`` (Sherwood number k d/D).

    Integrates d x_b/d zeta = -8 nu (1 - x_b)^2/((1 - x_in) Gz), nu = Da (sqrt(x_s) - sqrt(theta)),
    with x_s from the exact film at every point; with a constant ``sh`` the state is x_s
    (x_b explicit), otherwise x_b. Returns a dict with x_out, recovery, CPC_out, nu_in, nu_out.
    """
    const = not callable(sh)
    s_of = (lambda z: sh) if const else sh

    def rate(xb, xs):
        return -8.0 * da * (math.sqrt(xs) - math.sqrt(theta)) * (1.0 - xb) ** 2 / ((1.0 - x_in) * gz)

    def xb_of_xs(xs, beta):
        return 1.0 - (1.0 - xs) * math.exp(-beta * (math.sqrt(xs) - math.sqrt(theta)))

    if const and kind == "exact":
        beta = 2.0 * da / sh

        def rhs(z, y):
            xs = y[0]
            n = beta * (math.sqrt(xs) - math.sqrt(theta))
            dxb_dxs = math.exp(-n) * (1.0 + (1.0 - xs) * beta / (2.0 * math.sqrt(xs)))
            return [rate(xb_of_xs(xs, beta), xs) / dxb_dxs]
        y0 = [film_xs(x_in, beta, theta)]
        sol = solve_ivp(rhs, (0.0, 1.0), y0, method="DOP853", rtol=rtol, atol=1e-14)
        if sol.status != 0:
            raise RuntimeError(sol.message)
        xs_out = sol.y[0, -1]
        x_out = xb_of_xs(xs_out, beta)
    else:
        def rhs(z, y):
            xs = film_xs(y[0], 2.0 * da / s_of(max(z, 1e-300)), theta, kind)
            return [rate(y[0], xs)]
        sol = solve_ivp(rhs, (1e-12, 1.0), [x_in], method="DOP853", rtol=rtol, atol=1e-14)
        if sol.status != 0:
            raise RuntimeError(sol.message)
        x_out = sol.y[0, -1]
        xs_out = film_xs(x_out, 2.0 * da / s_of(1.0), theta, kind)
    f_out = (1.0 - x_in) / (1.0 - x_out)
    return {"x_out": x_out, "recovery": 1.0 - f_out * x_out / x_in,
            "cpc_out": cpc(xs_out, x_out, theta),
            "nu_out": da * (math.sqrt(xs_out) - math.sqrt(theta))}


class TubeMarch:
    """Steady 2-D laminar tube with a Pd-membrane wall, marched in the axial coordinate.

    Unknowns per axial step: x on n radial cells (field shape (n, 1)), the wall value x_w
    and the total molar flow f (both bordering the cell block). Conservative form

        (Gz/4) d[a(xi) f x]/d zeta + (1/xi) d/dxi[xi (nu_w w(xi) x - dx/dxi)] = 0,
        f' = -8 nu_w/Gz,  nu_w = Da (sqrt(x_w) - sqrt(theta)),

    a = 2(1 - xi^2), w = 2 xi - xi^3 (parabolic, low suction) or a = 1, w = xi (plug).

    Parameters
    ----------
    x_in, theta, gz, da : float
        Inlet H2 fraction, permeate/retentate pressure ratio, Graetz and Damkoehler numbers.
    n : int
        Radial cells; faces at xi = s (1 + clus (1 - s)), s uniform, so a refinement is nested.
    k : int
        Axial steps; zeta_j = (j/k)^p, every step refined by the same factor.
    profile : {'parabolic', 'plug'}
    wall : {'sieverts', 'dirichlet'}
        'dirichlet' fixes x_w = 0 with no suction (the Graetz problem, a validation).
    stefan : bool
        False drops the suction flow and the (1 - x_w) factor (low-flux wall, f = 1).
    gz_axial : float or None
        Graetz number used in the axial term only (a defect when it differs from ``gz``).
    conv : {'central', 'upwind'}
        Face value of the suction flux.
    scheme : {'bdf2', 'euler'}
        Axial step (variable-step BDF2 after a backward Euler start, or backward Euler).
    tol : float
        ``newton`` step tolerance per axial step.
    """

    def __init__(self, x_in, theta, gz, da, n=100, k=200, p=3.0, clus=0.8, profile="parabolic",
                 wall="sieverts", stefan=True, gz_axial=None, conv="central", scheme="bdf2", tol=1e-11):
        self.x_in, self.theta, self.gz, self.da = x_in, theta, gz, da
        self.gz_ax = gz if gz_axial is None else gz_axial
        self.n, self.k, self.wall, self.tol, self.conv, self.scheme = n, k, wall, tol, conv, scheme
        self.stefan = stefan and wall == "sieverts"
        s = np.linspace(0.0, 1.0, n + 1)
        self.xf = s * (1.0 + clus * (1.0 - s))
        self.xc = 0.5 * (self.xf[1:] + self.xf[:-1])
        self.zeta = (np.arange(k + 1) / k) ** p
        self.profile = profile
        self._build_operators()

    def _build_operators(self):
        xf, xc, n = self.xf, self.xc, self.n
        shape = (n, 1)
        # Centre xi = 0: symmetry, -dx/dxi = 0 (outward normal -xi): a = 1, b = 0, d = 0.
        # Wall xi = 1: x = x_w, Dirichlet with the wall value an external unknown (shapes_d, d = 1).
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
        self.G, _, self.Gb = construct_grad(shape, xf, xc, self.bc, shapes_d=(None, (1, 1)))
        self.Div = construct_div(shape, xf, nu=1)     # nu = 1: cylindrical, (1/xi) d(xi F)/dxi
        if self.profile == "parabolic":
            w = 2.0 * xf - xf ** 3                     # radial molar flux / nu_w, low suction
            vol_a = xf ** 2 - xf ** 4 / 2.0            # int xi a dxi, a = 2(1 - xi^2)
        else:
            w = xf.copy()                              # plug flow: radial flux nu_w xi
            vol_a = xf ** 2 / 2.0                      # int xi a dxi, a = 1
        vol = np.diff(xf ** 2) / 2.0                   # exact cell volumes / (2 pi), sum 1/2
        self.a = np.diff(vol_a) / vol                  # cell mean of a over the cell volume
        self.wcup = self.a * vol / np.sum(self.a * vol)
        # central face value of the suction flux: mean of the two upwind operators
        cp, _, cpb = construct_convflux_upwind(shape, xf, xc, self.bc, v=w.reshape(-1, 1),
                                               shapes_d=(None, (1, 1)))
        cm, _, cmb = construct_convflux_upwind(shape, xf, xc, self.bc, v=-w.reshape(-1, 1),
                                               shapes_d=(None, (1, 1)))
        if self.conv == "central":
            self.C, self.Cb = 0.5 * (cp - cm), 0.5 * (cpb - cmb)
        else:                                          # first-order upwind (suction is outward)
            self.C, self.Cb = cp, cpb
        self.DC, self.DG = self.Div @ self.C, self.Div @ self.G
        # wall-value columns of the two cell operators, as dense vectors (one column each)
        self.DCb, self.DGb = (self.Div @ self.Cb).toarray(), (self.Div @ self.Gb).toarray()
        self.gw, self.gwb = self.G[[n], :], float(self.Gb[n, 0])   # wall-face gradient row
        self.wall_jac = NumJac((1, 1))                 # d/dx_w of the pointwise wall flux
        self._jacobian_pattern()

    def _jacobian_pattern(self):
        """Fixed CSC pattern of the bordered Jacobian: tridiagonal cell block, wall row, two columns."""
        n = self.n
        band = sparse.diags_array([np.ones(n - 1), np.ones(n), np.ones(n - 1)], offsets=[-1, 0, 1],
                                  format="csc")
        for m in (self.DC, self.DG):                   # the cell block must be tridiagonal
            if abs(m - m.multiply(band)).sum() > 0:
                raise ValueError("cell operator wider than tridiagonal")
        i, j = band.nonzero()
        self._tri = (i, j)
        rows = np.concatenate([i, np.full(n, n), np.arange(n + 2), np.arange(n + 2)])
        cols = np.concatenate([j, np.arange(n), np.full(n + 2, n), np.full(n + 2, n + 1)])
        self._perm = np.lexsort((rows, cols))
        self._indices = rows[self._perm]
        self._indptr = np.concatenate([[0], np.cumsum(np.bincount(cols, minlength=n + 2))])
        dense = lambda m: np.asarray(m[i, j]).ravel()
        self._dc, self._dg, self._eye = dense(self.DC), dense(self.DG), (i == j).astype(float)
        self._gw_row = self.gw.toarray().ravel()

    def nu(self, xw):
        """Dimensionless wall flux, Sieverts: Da (sqrt(x_w) - sqrt(theta))."""
        if self.wall == "dirichlet":
            return np.zeros_like(xw)
        return self.da * (np.sqrt(np.maximum(xw, 0.0)) - math.sqrt(self.theta))

    def _residual(self, u, c0, hist_fx, hist_f, h):
        n, i = self.n, self._tri[0]
        x, xw, f = u[:n], u[n], u[n + 1]
        nu, dnu = self.wall_jac(self.nu, np.array([[xw]]))   # value and d nu/d x_w
        nu, dnu = float(nu.ravel()[0]), float(dnu.toarray()[0, 0])
        nu_c, dnu_c = (nu, dnu) if self.stefan else (0.0, 0.0)
        acc = 0.25 * self.gz_ax * self.a
        conv = self.DC @ x + self.DCb[:, 0] * xw
        r_x = acc * (c0 * f * x - hist_fx) / h + nu_c * conv - (self.DG @ x + self.DGb[:, 0] * xw)
        gw = float(self._gw_row @ x) + self.gwb * xw
        sieverts = self.wall == "sieverts"
        if sieverts:   # wall row: -dx/dxi = nu_w (1 - x_w), or nu_w without the Stefan factor
            stf = (1.0 - xw) if self.stefan else 1.0
            r_w, d_w = -gw - nu * stf, -self.gwb - dnu * stf + (nu if self.stefan else 0.0)
        else:          # Graetz validation: x_w = 0
            r_w, d_w = xw, 1.0
        if self.stefan:
            r_f, d_fw, d_ff = (c0 * f - hist_f) / h + 8.0 * nu / self.gz_ax, 8.0 * dnu / self.gz_ax, c0 / h
        else:
            r_f, d_fw, d_ff = f - 1.0, 0.0, 1.0
        data = np.concatenate([
            self._eye * (acc * c0 * f / h)[i] + nu_c * self._dc - self._dg,
            -self._gw_row if sieverts else np.zeros(n),
            np.concatenate([nu_c * self.DCb[:, 0] + dnu_c * conv - self.DGb[:, 0], [d_w, d_fw]]),
            np.concatenate([acc * c0 * x / h, [0.0, d_ff]])])
        jac = sparse.csc_array((data[self._perm], self._indices, self._indptr), shape=(n + 2, n + 2))
        return np.concatenate([r_x, [r_w, r_f]]), jac

    def run(self):
        """March from the inlet; returns self with per-step arrays zeta, f, x_cup, x_w, nu_w, g_w."""
        n = self.n
        u = np.concatenate([np.full(n, self.x_in), [self.x_in, 1.0]])
        hist = [u.copy()]
        out = {"x_cup": [self.x_in], "x_w": [np.nan], "f": [1.0], "nu": [np.nan], "g_w": [np.nan]}
        for j in range(1, self.k + 1):
            h = self.zeta[j] - self.zeta[j - 1]
            if j == 1 or self.scheme == "euler":         # backward Euler (start)
                c0, hfx, hf = 1.0, hist[-1][n + 1] * hist[-1][:n], hist[-1][n + 1]
            else:                                        # variable-step BDF2
                om = h / (self.zeta[j - 1] - self.zeta[j - 2])
                c0 = (1.0 + 2.0 * om) / (1.0 + om)
                a1, a2 = (1.0 + om), om * om / (1.0 + om)
                fx = [hh[n + 1] * hh[:n] for hh in hist[-2:]]
                hfx = a1 * fx[1] - a2 * fx[0]
                hf = a1 * hist[-1][n + 1] - a2 * hist[-2][n + 1]
            guess = u.copy()
            if j == 1:
                guess[n] = 0.5 * (self.x_in + max(self.theta, 1e-12))
            sol = newton(self._residual, guess, args=(c0, hfx, hf, h), tol=self.tol, maxfev=60,
                         stagnation=True)
            if not sol.success:
                raise RuntimeError(f"newton failed at step {j}: {sol.message}")
            u = sol.x
            hist = [hist[-1], u.copy()]
            self._record(u, out)
        for key, v in out.items():
            setattr(self, key, np.asarray(v))
        return self

    def _record(self, u, out):
        n = self.n
        x = u[:n].reshape(-1, 1)
        bc = (self.bc[0], {"a": 0.0, "b": 1.0, "d": u[n]})   # the same bc, wall value filled in
        _, _, xw, gw = compute_boundary_values(x, self.xf, self.xc, bc)
        out["x_cup"].append(float(self.wcup @ u[:n]))
        out["x_w"].append(float(xw[0, 0]))
        out["f"].append(float(u[n + 1]))
        out["nu"].append(float(self.nu(np.array([[u[n]]]))[0, 0]))
        out["g_w"].append(float(gw[0, 0]))

    # results ------------------------------------------------------------------------------
    def recovery(self):
        """Fraction of the inlet H2 that has permeated at the outlet."""
        return 1.0 - self.f[-1] * self.x_cup[-1] / self.x_in

    def cpc_local(self):
        """Local CPC of the 2-D solve per step (Eq. 13 with the cup-mean bulk fraction)."""
        sq = math.sqrt(self.theta)
        return 1.0 - (np.sqrt(self.x_w[1:]) - sq) / (np.sqrt(self.x_cup[1:]) - sq)

    def sherwood_local(self):
        """Local Sherwood number k d/D = 2 (-dx/dxi)_w/(x_cup - x_w) per step."""
        return 2.0 * (-self.g_w[1:]) / (self.x_cup[1:] - self.x_w[1:])

    def n2_balance(self):
        """Largest |f (1 - x_cup) - (1 - x_in)| over the steps (N2 conservation)."""
        return float(np.max(np.abs(self.f * (1.0 - self.x_cup) - (1.0 - self.x_in))))
