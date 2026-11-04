"""Doyle, Fuller and Newman (1993): lithium foil | polymer separator | TiS2 cathode.

The pseudo-two-dimensional cell on one finite-volume grid with pymrm operators:
salt concentration, electrolyte potential and solid lithium concentration in one
Newton system per backward-Euler step. Parameters come in as a dictionary ``P``
keyed by the symbols of the page's parameter table; nothing here loads data.
"""

from __future__ import annotations

import numpy as np
from pymrm import (NumJac, construct_boundary_value_matrices, construct_div, construct_grad,
                   interp_cntr_to_stagg, newton)
from pymrm.checks import find_roots

R_GAS = 8.3143      # J/(mol K), the value of the paper's era
F_CONST = 96487.0   # C/mol


def rtf(P):
    """Thermal voltage RT/F in volts."""
    return R_GAS * P["T"] / F_CONST


def ocp(u, P):
    """Eq. 16, the open-circuit potential of the insertion material at utilisation ``u``.

    Parameters
    ----------
    u : float or numpy.ndarray
        Utilisation c_s / c_T.
    P : dict
        Parameters (``U_theta``, ``beta``, ``zeta``, ``c_T``, ``T``).

    Returns
    -------
    numpy.ndarray
        Potential in volts against lithium.
    """
    cs = np.asarray(u, float) * P["c_T"]
    return P["U_theta"] + rtf(P) * (np.log((P["c_T"] - cs) / cs) + P["beta"] * cs + P["zeta"])


def s_c(I, P):
    """Eq. 26: particle diffusion time over discharge time at current density ``I``."""
    return P["R_s"] ** 2 * I / (P["D_s"] * F_CONST * (1 - P["epsilon"]) * P["c_T"] * P["delta_c"])


def s_s(I, P):
    """Eq. 27: electrolyte diffusion time over discharge time at current density ``I``."""
    return ((P["delta_s"] + P["delta_c"]) ** 2 * I
            / (P["D"] * F_CONST * (1 - P["epsilon"]) * P["c_T"] * P["delta_c"]))


def i0_cathode(c, cs, P):
    """Eq. 30, exchange current density of the insertion reaction in A/m2."""
    aa, ac = P["alpha_a"], P["alpha_c"]
    return (F_CONST * P["k_2"] * (P["c_max"] - c) ** ac * c ** aa
            * (P["c_T"] - cs) ** aa * cs ** ac)


def kappa_from_delta(delta, P, I=10.0):
    """Invert Eq. 28 for the effective electrolyte conductivity in S/m.

    Parameters
    ----------
    delta : float
        The printed value of the group delta.
    P : dict
        Parameters (``alpha_a``, ``delta_c``, ``T``, ``sigma``).
    I : float, optional
        Current density of the printed value, A/m2.
    """
    pre = P["alpha_a"] * F_CONST * I * P["delta_c"] / (R_GAS * P["T"])
    return 1.0 / (delta / pre - 1.0 / P["sigma"])


def kappa_from_nu(nu, P):
    """Invert Eq. 29 at the initial state for the effective conductivity in S/m."""
    a_spec = 3.0 * (1 - P["epsilon"]) / P["R_s"]
    i0 = i0_cathode(P["c_0"], P["u_0"] * P["c_T"], P)
    pre = ((P["alpha_a"] + P["alpha_c"]) * F_CONST * a_spec * i0 * P["delta_c"] ** 2
           / (R_GAS * P["T"]))
    return 1.0 / (nu ** 2 / pre - 1.0 / P["sigma"])


def crossing(t, v, level, falling=True):
    """The one time at which a piecewise-linear trace passes ``level``.

    The trace is the backward-Euler solution joined linearly between step times
    (the integrator's own dense output, accurate to second order in the step).
    Its root is located with :func:`pymrm.checks.find_roots` on the whole span of
    the trace, scanned finer than the steps, and the function raises unless
    there is exactly one.

    Parameters
    ----------
    t, v : numpy.ndarray
        Abscissa (increasing) and trace values.
    level : float
        The level to cross.
    falling : bool, optional
        Sign convention only: a falling trace is negated, so the root is the
        same either way.

    Returns
    -------
    float
        The abscissa of the crossing.
    """
    t, v = np.asarray(t, float), np.asarray(v, float)
    s = -1.0 if falling else 1.0
    res = find_roots(lambda x: s * (np.interp(x, t, v) - level), t[0], t[-1],
                     n_scan=4 * len(t), xtol=1e-12 * max(1.0, abs(t[-1])))
    if len(res["roots"]) != 1:
        raise RuntimeError(f"{len(res['roots'])} crossings of {level} on "
                           f"[{t[0]:.4g}, {t[-1]:.4g}]: {res['message']}")
    return float(res["roots"][0])


class Cell:
    """The 1993 cell sandwich on one grid with a face at the separator/cathode boundary.

    State ``y`` of shape ``(n_x, 3)``: ``ln(c/c_0)``, ``Phi_2 F/(RT)`` and
    ``c_s/c_T``. The matrix potential is gauged to zero, so the overpotential is
    ``-Phi_2``.

    Parameters
    ----------
    P : dict
        Parameters keyed by the symbols of the parameter table.
    kappa : float
        Bulk electrolyte conductivity in S/m; ``kappa_eff = kappa eps**1.5``.
    n_s, n_c : int, optional
        Cells in the separator and in the cathode.
    nu_factor : float, optional
        Coefficient of ``(1 - t0+) grad ln c`` in Eq. 5 in units of RT/F; 2 for a
        1:1 salt, 1 as Eq. 5 is printed.
    p_D : float, optional
        Bruggeman exponent of the salt diffusivity: ``D_eff = D eps**p_D``.
    eps, delta_c : float, optional
        Cathode porosity and thickness (default: the parameter table).
    mig_sign : float, optional
        Sign of the migration term in the salt flux, ``-1`` as derived (defect
        injection only).
    t0_form : {"1-t0", "t0"}, optional
        Factor multiplying the current in the salt flux and Eq. 5 (defect
        injection only).
    ai_solid, ai_charge : float, optional
        Multipliers of the pore-wall current in the solid and charge balances
        (defect injection only).
    anion_blocking : bool, optional
        Zero salt flux at the lithium face; False imposes a zero concentration
        gradient there instead (defect injection only).
    face_mean : {"harmonic", "arithmetic"}, optional
        Mean of the cell transport coefficients on interior faces; the porosity
        jumps on a face, where only the harmonic mean is consistent (defect
        injection only).
    """

    def __init__(self, P, kappa, n_s=24, n_c=48, nu_factor=2.0, p_D=0.5, eps=None,
                 delta_c=None, mig_sign=-1.0, t0_form="1-t0", ai_solid=1.0, ai_charge=1.0,
                 anion_blocking=True, face_mean="harmonic"):
        self.P = dict(P)
        self.face_mean = face_mean
        self.kappa, self.nu_factor, self.p_D = kappa, nu_factor, p_D
        self.eps_c = P["epsilon"] if eps is None else eps
        self.dc = P["delta_c"] if delta_c is None else delta_c
        self.mig_sign, self.t0_form = mig_sign, t0_form
        self.ai_solid, self.ai_charge = ai_solid, ai_charge
        self.anion_blocking = anion_blocking
        self.n_s, self.n_x = n_s, n_s + n_c
        self.rtf = rtf(P)
        self._build_grid()
        self._build_operators()
        # discharge time of the full theoretical capacity is t_full / |I|
        self.t_full = (1 - self.eps_c) * P["c_T"] * self.dc * F_CONST
        self.failures = []      # why each failed newton solve stopped

    def _build_grid(self):
        P, ds = self.P, self.P["delta_s"]
        self.x_f = np.concatenate([np.linspace(0.0, ds, self.n_s + 1),
                                   np.linspace(ds, ds + self.dc, self.n_x - self.n_s + 1)[1:]])
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.dx = np.diff(self.x_f)
        sep = np.arange(self.n_x) < self.n_s
        self.cathode = ~sep
        self.eps = np.where(sep, 1.0, self.eps_c)          # the separator is pure polymer
        self.a = np.where(sep, 0.0, 3.0 * (1 - self.eps_c) / P["R_s"])

    def _to_face(self, b):
        """Face values of a cell property: harmonic mean inside (the porosity jumps
        on a face), the adjacent cell on the two boundary faces."""
        w = self.dx
        bf = np.empty(self.n_x + 1)
        if self.face_mean == "harmonic":
            bf[1:-1] = (w[:-1] + w[1:]) / (w[:-1] / b[:-1] + w[1:] / b[1:])
        else:
            bf[1:-1] = (w[:-1] * b[:-1] + w[1:] * b[1:]) / (w[:-1] + w[1:])
        bf[0], bf[-1] = b[0], b[-1]
        return bf.reshape(-1, 1)

    def _build_operators(self):
        shape = (self.n_x, 1)
        self.kf = self.kappa * self._to_face(self.eps ** 1.5)     # kappa_eff on faces
        self.Df = self.P["D"] * self._to_face(self.eps ** self.p_D)
        # One gradient operator for c, ln c and Phi_2, Neumann at both ends, with
        # the value at x = 0 an external datum (shapes_d, d = 1): its face
        # gradient is -g for the datum g, because the outward normal there is -x.
        # The data are the physical boundary gradients each step computes from
        # the state (bc_left below); at the collector all three gradients vanish.
        self.bc_unit = ({"a": 1.0, "b": 0.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        self.grad, self.grad_bc_left, self.grad_bc_right = construct_grad(
            shape, self.x_f, self.x_c, self.bc_unit, shapes_d=((1, 1), None))
        # the face value at x = 0 under the same condition, as a matrix on the cells
        # plus one on the datum: value = bval @ f + bval_bc @ [g]
        self.bval, self.bval_bc = construct_boundary_value_matrices(
            shape, self.x_f, self.x_c, bc=self.bc_unit[0], bound_id=0, shape_d=(1, 1))
        self.div = construct_div(shape, self.x_f, nu=0)          # nu=0: Cartesian slab
        # pointwise coupling of the three fields, plus the neighbours along x that
        # the face fluxes read
        self.numjac = NumJac((self.n_x, 3), axes_diagonals=[0])

    # -- physical properties ------------------------------------------------
    def t0_plus(self, c):
        """Appendix A: transference number of Li+, linear in c."""
        return self.P["t0_a"] + self.P["t0_b"] * c

    def _t0_factor(self, t0):
        return (1 - t0) if self.t0_form == "1-t0" else t0

    def unpack(self, y):
        """Salt concentration, electrolyte potential (V) and solid concentration."""
        y = np.asarray(y).reshape(self.n_x, 3)
        return (self.P["c_0"] * np.exp(y[:, 0]), y[:, 1] * self.rtf,
                y[:, 2] * self.P["c_T"])

    def bc_left(self, c, I):
        """Boundary data at x = 0 for c, ln c and Phi_2, from the physics.

        No anion enters the lithium: N = -D_eff dc/dx - (1 - t0) i2 / F = 0 with
        i2 = I, so -dc/dx = (1 - t0) I / (F D_eff). The current is all in the
        electrolyte: i2 = -kappa_eff dPhi_2/dx + nu kappa_eff (RT/F)(1 - t0)
        dln c/dx = I, so -dPhi_2/dx = I/kappa_eff - nu (RT/F)(1 - t0) dln c/dx.
        Written on the outward normal (-x) these are the d of a Neumann
        dictionary {a: 1, b: 0, d}. t0 is taken at the first cell.
        """
        t0b = self._t0_factor(self.t0_plus(c[0]))
        kf0, Df0 = self.kf[0, 0], self.Df[0, 0]
        g_c = -self.mig_sign * t0b * I / (F_CONST * Df0) if self.anion_blocking else 0.0
        c_b = self.face_value(c, g_c)
        g_w = g_c / c_b                                   # -dln c/dx
        g_p = I / kf0 + self.nu_factor * self.rtf * t0b * g_w   # -dPhi_2/dx
        return g_c, g_w, g_p, c_b

    def face_value(self, f, g):
        """Value at x = 0 of a cell field whose datum there is ``g`` (= -df/dx)."""
        return float((self.bval @ f.reshape(-1, 1) + self.bval_bc @ np.array([[g]]))[0, 0])

    def fluxes(self, c, phi2, I):
        """Face current density i2 (A/m2) and salt flux N (mol/(m2 s))."""
        g_c, g_w, g_p, _ = self.bc_left(c, I)
        cc, lc, pp = c.reshape(-1, 1), np.log(c).reshape(-1, 1), phi2.reshape(-1, 1)
        gc = self.grad @ cc + self.grad_bc_left @ np.array([[g_c]])
        glc = self.grad @ lc + self.grad_bc_left @ np.array([[g_w]])
        gp = self.grad @ pp + self.grad_bc_left @ np.array([[g_p]])
        t0f = interp_cntr_to_stagg(self.t0_plus(cc), self.x_f, self.x_c)
        t0f[0], t0f[-1] = self.t0_plus(c[0]), self.t0_plus(c[-1])   # as in bc_left
        w = self._t0_factor(t0f)
        i2 = -self.kf * gp + self.nu_factor * self.kf * self.rtf * w * glc
        n_salt = -self.Df * gc + self.mig_sign * w * i2 / F_CONST
        return i2.ravel(), n_salt.ravel()

    def i_transfer(self, c, phi2, cs):
        """Eq. 17 with Phi_1 = 0, so eta = -Phi_2; A/m2 of pore wall.

        Valid for alpha_a + alpha_c = 1.
        """
        P = self.P
        U_prime = P["U_theta"] + self.rtf * (P["beta"] * cs + P["zeta"])     # Eqs. 18, 20
        z = (-phi2 - U_prime) / self.rtf
        pre = F_CONST * P["k_2"] * (P["c_max"] - c) ** P["alpha_c"] * c ** P["alpha_a"]
        return pre * (cs * np.exp(P["alpha_a"] * z) - (P["c_T"] - cs) * np.exp(-P["alpha_c"] * z))

    # -- residual -----------------------------------------------------------
    def residual(self, y, y_old, dt, I):
        """Scaled backward-Euler residual of the three balances, shape (n_x, 3)."""
        P = self.P
        c, phi2, cs = self.unpack(y)
        c_o, _, cs_o = self.unpack(y_old)
        i2, n_salt = self.fluxes(c, phi2, I)
        ai = self.a * self.i_transfer(c, phi2, cs)
        r_c = self.eps * (c - c_o) / dt + (self.div @ n_salt.reshape(-1, 1)).ravel()
        r_p = (self.div @ i2.reshape(-1, 1)).ravel() - self.ai_charge * ai
        r_s = np.where(self.cathode,
                       (1 - self.eps_c) * (cs - cs_o) / dt + self.ai_solid * ai / F_CONST,
                       cs)                    # no solid in the separator: c_s = 0
        ts = self.t_full / abs(I)
        return np.stack([r_c * ts / P["c_0"], r_p * self.dc / abs(I), r_s * ts / P["c_T"]],
                        axis=1)

    # -- outputs ------------------------------------------------------------
    def boundary(self, y, I):
        """Face values at x = 0: salt concentration and electrolyte potential."""
        c, phi2, _ = self.unpack(y)
        _, _, g_p, c_b = self.bc_left(c, I)
        return c_b, self.face_value(phi2, g_p)

    def cell_potential(self, y, I):
        """Phi_1(collector) - Phi_1(lithium), Eq. 22, from the face values at x = 0."""
        P = self.P
        c_b, phi2_b = self.boundary(y, I)
        i01 = P["i_0_1"] * np.sqrt((P["c_max"] - c_b) * c_b
                                   / ((P["c_max"] - P["c_0"]) * P["c_0"]))   # Eq. 8
        eta_s1 = 2 * self.rtf * np.arcsinh(I / (2 * i01))     # Eq. 6, alpha_a1 = alpha_c1
        return -phi2_b - eta_s1 - I * self.dc / (2 * P["sigma"])

    def utilisation(self, y):
        """Eq. 25: mean solid concentration over c_T."""
        cs = self.unpack(y)[2]
        m = self.cathode
        return float(np.sum(cs[m] * self.dx[m]) / (self.P["c_T"] * self.dc))

    def salt(self, y):
        """Total salt per unit area, mol/m2."""
        return float(np.sum(self.eps * self.dx * self.unpack(y)[0]))

    def wall_flux_mean(self, y):
        """Mean pore-wall flux j_n = i/F over the cathode thickness, mol/(m2 s)."""
        c, phi2, cs = self.unpack(y)
        jn = self.i_transfer(c, phi2, cs) / F_CONST
        m = self.cathode
        return float(np.sum(jn[m] * self.dx[m]) / self.dc), jn[m]

    # -- time marching ------------------------------------------------------
    def initial(self, u0):
        """Uniform salt at c_0, solid at u0, electrolyte at the open-circuit level."""
        y = np.zeros((self.n_x, 3))
        y[:, 1] = -ocp(u0, self.P) / self.rtf
        y[:, 2] = np.where(self.cathode, u0, 0.0)
        return y

    def schedule(self, I, u0, u_end, nodes=(), dt_max=6e-3):
        """Base step times from 0 to the time at which the utilisation reaches
        ``u_end``: dt grows by 1.3 from 2e-5 to ``dt_max`` of t_full/|I|, with
        ``nodes`` inserted as step times."""
        t_ref = self.t_full / abs(I)
        t_end = t_ref * abs(u_end - u0)
        ts, t, dt = [0.0], 0.0, 2e-5 * t_ref
        while t < t_end:
            t = min(t + dt, t_end)
            ts.append(t)
            dt = min(1.3 * dt, dt_max * t_ref)
        extra = [x for x in nodes if 0.0 < x < t_end]
        ts = np.unique(np.concatenate([ts, extra]))
        keep = np.concatenate([[True], np.diff(ts) > 1e-9 * t_ref])
        ts = ts[keep]
        for x in extra:   # a node merged with a nearby step time takes its place exactly
            ts[np.argmin(np.abs(ts - x))] = x
        return ts

    def step(self, y_old, dt, I, tol=1e-9):
        """One backward-Euler step; the newton result, or None when it does not converge."""
        y_flat_old = y_old.ravel()

        def fun(z):
            return self.numjac(lambda w: self.residual(w, y_flat_old, dt, I), z)

        try:
            # a trial iterate may overflow the exponentials or take ln of a negative
            # number; newton's damping rejects such points (non-finite residual),
            # so the floating-point warnings they raise carry no information
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                sol = newton(fun, y_old.copy(), tol=tol, maxfev=30)
        except ValueError as err:       # NewtonBreakdown: singular or non-finite
            self.failures.append(str(err).split(":")[0])
            return None
        if not sol.success:
            self.failures.append(sol.message)
            return None
        return sol

    def march(self, I, u0=None, v_stop=1.55, v_top=3.2, u_stop=(0.9995, 1e-4), nodes=(),
              refine=1, t_stop=None, max_halvings=6, tol=1e-9, dt_max=6e-3,
              stop_on_failure=False):
        """Galvanostatic backward-Euler march at current density ``I`` (A/m2).

        Each base step of :meth:`schedule` is split into ``refine`` equal steps,
        so a refined march contains every step time of a coarser one (``dt_max`` is
        the ceiling of the base step as a fraction of t_full/|I|). A step whose
        newton solve fails is retried as two half steps, recursively, at most
        ``max_halvings`` levels deep; a failure beyond that raises. The march stops
        after the step that takes the cell potential below ``v_stop`` (discharge)
        or above ``v_top`` (charge), at ``t_stop``, or where the utilisation
        reaches ``u_stop[0]`` (discharge) or ``u_stop[1]`` (charge). With
        ``stop_on_failure`` a step that fails at the deepest halving ends the
        march instead of raising (``end`` is then ``"newton"`` and ``failure``
        says where); every accepted step is a converged solve either way.

        Returns
        -------
        dict
            ``t``, ``u``, ``V``, ``salt``, ``c_min`` (arrays from t = 0), ``snaps`` (state at
            each requested node reached), ``halvings`` (number of split steps),
            ``iterations`` (newton iterations, largest per step) and ``end``
            (why it stopped).
        """
        u0 = self.P["u_0"] if u0 is None else u0
        base = self.schedule(I, u0, u_stop[0] if I > 0 else u_stop[1], nodes, dt_max)
        node_set = {float(x) for x in nodes}
        y = self.initial(u0)
        t = 0.0
        out = {"t": [0.0], "u": [self.utilisation(y)], "V": [self.cell_potential(y, I)],
               "salt": [self.salt(y)], "c_min": [self.P["c_0"]], "snaps": {}, "halvings": 0,
               "iterations": 0,
               "end": "utilisation"}

        def advance(y0, dt, depth):
            sol = self.step(y0, dt, I, tol)
            if sol is not None:
                out["iterations"] = max(out["iterations"], sol.nit)
                return sol.x.reshape(self.n_x, 3)
            if depth >= max_halvings:
                raise RuntimeError(f"newton failed after t = {t:.6g} s with dt = {dt:.3g} s "
                                   f"at I = {I} A/m2, u = {self.utilisation(y0):.4f}")
            out["halvings"] += 1
            y1 = advance(y0, 0.5 * dt, depth + 1)
            return advance(y1, 0.5 * dt, depth + 1)

        for k in range(len(base) - 1):
            t0, t1 = base[k], base[k + 1]
            for j in range(refine):
                ta = t0 + (t1 - t0) * j / refine
                tb = t1 if j == refine - 1 else t0 + (t1 - t0) * (j + 1) / refine
                try:
                    y = advance(y, tb - ta, 0)
                except RuntimeError as err:
                    if not stop_on_failure:
                        raise
                    out["end"], out["failure"] = "newton", str(err)
                    break
                t = tb
                out["t"].append(t)
                out["u"].append(self.utilisation(y))
                v = self.cell_potential(y, I)
                out["V"].append(v)
                out["salt"].append(self.salt(y))
                out["c_min"].append(float(self.unpack(y)[0].min()))
            if out["end"] == "newton":
                break
            if t1 in node_set:
                out["snaps"][t1] = y.copy()
            if (I > 0 and v < v_stop) or (I < 0 and v > v_top):
                out["end"] = "potential"
                break
            if t_stop is not None and t >= t_stop:
                out["end"] = "t_stop"
                break
        for key in ("t", "u", "V", "salt", "c_min"):
            out[key] = np.array(out[key])
        return out
