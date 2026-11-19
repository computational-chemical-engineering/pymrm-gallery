"""Butler-Volmer kinetics, its linear and Tafel limits, and the Doyle cell (page J3.1).

Doyle, Fuller and Newman, J. Electrochem. Soc. 140, 1526 (1993). Two parts:

- :class:`Kinetics`, the interfacial law of the paper for one parameter set:
  the open-circuit potential (Eq. 16) and its reduced form (Eqs. 18, 20), the
  insertion kinetics exactly as printed (Eq. 17), the exchange current density
  (Eq. 30) and its general-alpha counterpart, the standard two-exponential
  form (Eq. 6), its linear and Tafel limits, the charge-transfer resistance,
  and the overpotentials at which each limit deviates from the full law by a
  given fraction (:meth:`Kinetics.eta_tafel_threshold`,
  :meth:`Kinetics.eta_linear_threshold`, located with
  ``pymrm.checks.find_roots``);
- :class:`Cell`, the one-dimensional porous-electrode cell of page ``J3.4``
  (lithium foil, separator, porous cathode) with the kinetic law of the
  cathode selectable, marched in time by implicit Euler with ``newton`` on a
  ``NumJac`` Jacobian.

Every defect the page's break table injects enters through an argument of
these functions or through the parameter dictionary.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from pymrm import (construct_grad, construct_div, NumJac, newton,
                   stencil_block_diagonals)
from pymrm.checks import find_roots

R_GAS, F_CONST = 8.3143, 96487.0     # J/(mol K), C/mol: the constants J3.4 uses


def one_root(f, lo, hi, what, **kw):
    """The single root of ``f`` on ``[lo, hi]``; raises unless there is exactly one.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        The physical range scanned.
    what : str
        Name of the quantity, for the error message.
    **kw
        Passed to ``pymrm.checks.find_roots`` (``n_scan``, ``log``, ``xtol``).

    Returns
    -------
    float
        The root.
    """
    res = find_roots(f, lo, hi, **kw)
    if len(res["roots"]) != 1 or res["failed"]:
        raise RuntimeError(f"{what}: {res['message']} on [{lo:g}, {hi:g}]")
    return res["roots"][0]


def kappa_from_delta(par, delta, I=10.0):
    """Effective electrolyte conductivity from Doyle's Eq. 28, as reconstructed on J3.4.

    Parameters
    ----------
    par : dict
        Doyle's parameters (``alpha_a``, ``delta_c``, ``T``, ``sigma``).
    delta : float
        The printed dimensionless group delta.
    I : float
        Current density at which delta is quoted, A/m2.

    Returns
    -------
    float
        Effective conductivity, S/m.
    """
    pre = par["alpha_a"] * F_CONST * I * par["delta_c"] / (R_GAS * par["T"])
    return 1.0 / (delta / pre - 1.0 / par["sigma"])


class Kinetics:
    """The interfacial law of Doyle et al. (1993) for one parameter set.

    Parameters
    ----------
    par : dict
        Doyle's parameters by symbol, as in J3.4's parameter file (``T``,
        ``c_T``, ``c_max``, ``k_2``, ``alpha_a``, ``alpha_c``, ``U_theta``,
        ``beta``, ``zeta``, ...).

    Notes
    -----
    ``aa`` and ``ac`` default to Doyle's ``alpha_a``; every method takes
    them so that the transfer coefficients can be swept off the published
    pair.
    """

    def __init__(self, par):
        self.par = dict(par)
        self.rtf = R_GAS * self.par["T"] / F_CONST
        self.alpha = self.par["alpha_a"]
        self.c_T, self.c_max = self.par["c_T"], self.par["c_max"]

    def _a(self, aa, ac):
        return (self.alpha if aa is None else aa), (self.alpha if ac is None else ac)

    def U_ocp(self, u):
        """Eq. 16, the open-circuit potential at utilisation ``u = c_s/c_T``, V."""
        cs = np.asarray(u, float) * self.c_T
        p = self.par
        return p["U_theta"] + self.rtf * (np.log((self.c_T - cs) / cs) + p["beta"] * cs
                                          + p["zeta"])

    def U_prime(self, cs):
        """Eq. 18 with the Eq. 20 constants, V."""
        p = self.par
        return p["U_theta"] + self.rtf * (p["beta"] * np.asarray(cs, float) + p["zeta"])

    def i_eq17(self, eta, c, cs, aa=None, ac=None, pre=None, brk=None):
        """Eq. 17 exactly as printed, ``eta = Phi_1 - Phi_2`` (Eq. 19), A/m2.

        Parameters
        ----------
        eta : array_like
            Electrode potential against the solution, V.
        c, cs : array_like
            Electrolyte and solid-surface concentrations, mol/m3.
        aa, ac : float, optional
            Transfer coefficients.
        pre : tuple, optional
            Exponents on ``c`` and on ``c_max - c`` in the prefactor; default
            ``(aa, ac)``, as printed.
        brk : tuple, optional
            Exponents in the anodic and the cathodic exponential; default
            ``(aa, ac)``, as printed.
        """
        aa, ac = self._a(aa, ac)
        pa, pc = (aa, ac) if pre is None else pre
        ba, bc = (aa, ac) if brk is None else brk
        c, cs = np.asarray(c, float), np.asarray(cs, float)
        z = (np.asarray(eta, float) - self.U_prime(cs)) / self.rtf
        p = F_CONST * self.par["k_2"] * (self.c_max - c) ** pc * c ** pa
        return p * (cs * np.exp(ba * z) - (self.c_T - cs) * np.exp(-bc * z))

    def i0_eq30(self, c, cs, aa=None, ac=None, pre=None, tail=None):
        """Eq. 30, the exchange current density of the insertion reaction, A/m2.

        ``pre`` as in :meth:`i_eq17`; ``tail`` gives the exponents on
        ``c_T - c_s`` and on ``c_s``, default ``(aa, ac)`` as printed.
        """
        aa, ac = self._a(aa, ac)
        pa, pc = (aa, ac) if pre is None else pre
        ta, tc = (aa, ac) if tail is None else tail
        c, cs = np.asarray(c, float), np.asarray(cs, float)
        return (F_CONST * self.par["k_2"] * (self.c_max - c) ** pc * c ** pa
                * (self.c_T - cs) ** ta * cs ** tc)

    def eta_equilibrium(self, cs, aa, ac, drop_n=False):
        """Zero-current potential of Eq. 17 in closed form, ``U' + (RT/nF) ln((c_T-c_s)/c_s)``.

        ``drop_n=True`` omits the ``1/n``, ``n = aa + ac`` (a defect).
        """
        cs = np.asarray(cs, float)
        n = 1.0 if drop_n else aa + ac
        return self.U_prime(cs) + self.rtf * np.log((self.c_T - cs) / cs) / n

    def i0_general(self, c, cs, aa, ac):
        """Common value of the two partial currents of Eq. 17 at equilibrium, any alpha, A/m2."""
        c, cs = np.asarray(c, float), np.asarray(cs, float)
        L = np.log((self.c_T - cs) / cs)
        return (F_CONST * self.par["k_2"] * (self.c_max - c) ** ac * c ** aa
                * cs * np.exp(aa * L / (aa + ac)))

    def i_bv(self, eta_s, i0, aa=None, ac=None):
        """The standard form, Eq. 6 with a general ``i0``: ``i0 (e^(aa f eta_s) - e^(-ac f eta_s))``."""
        aa, ac = self._a(aa, ac)
        z = np.asarray(eta_s, float) / self.rtf
        return np.asarray(i0, float) * (np.exp(aa * z) - np.exp(-ac * z))

    def i_linear(self, eta_s, i0, aa=None, ac=None):
        """Low-overpotential limit, ``i0 n f eta_s``: a charge-transfer resistance."""
        aa, ac = self._a(aa, ac)
        return np.asarray(i0, float) * (aa + ac) * np.asarray(eta_s, float) / self.rtf

    def i_tafel(self, eta_s, i0, aa=None, ac=None, branch="anodic"):
        """High-overpotential limit: one exponential, the branch chosen beforehand."""
        aa, ac = self._a(aa, ac)
        e = np.asarray(eta_s, float) / self.rtf
        if branch == "anodic":
            return np.asarray(i0, float) * np.exp(aa * e)
        return -np.asarray(i0, float) * np.exp(-ac * e)

    def R_ct(self, i0, aa=None, ac=None):
        """Charge-transfer resistance ``RT/(n F i0)``, ``n = aa + ac``, Ohm m2."""
        aa, ac = self._a(aa, ac)
        return self.rtf / ((aa + ac) * np.asarray(i0, float))

    # ----------------------------------------------------- validity thresholds
    def eta_tafel_threshold(self, eps, n=1.0, rtf=None, form="log1p"):
        """Closed form: the overpotential at which the Tafel current deviates by ``eps``.

        ``RT ln(1 + 1/eps) / (n F)``. ``form="log"`` writes ``ln(1/eps)`` (a
        defect).
        """
        rtf = self.rtf if rtf is None else rtf
        e = np.asarray(eps, float)
        return rtf / n * (np.log1p(1.0 / e) if form == "log1p" else np.log(1.0 / e))

    def eta_tafel_threshold_numeric(self, eps, aa, ac, rtf=None):
        """The same, the single root of the Tafel deviation on the full law, V."""
        f = 1.0 / (self.rtf if rtf is None else rtf)

        def dev(e):
            full = np.exp(aa * f * e) - np.exp(-ac * f * e)
            return np.exp(aa * f * e) / full - 1.0 - eps

        return one_root(dev, 1e-6, 5.0, "Tafel threshold", log=True, n_scan=128,
                        xtol=1e-16)

    def eta_linear_threshold(self, eps, aa, ac, sign=+1, rtf=None, crossing="first"):
        """Overpotential at which the linear form deviates by ``eps`` in current, V.

        For ``aa != ac`` the deviation is not monotonic in eta: the linear form
        overshoots, comes back through the exact law and then falls away, so
        ``|dev| - eps`` has up to four roots on ``[1e-6, 4]`` V. The validity
        boundary is the smallest (``crossing="first"``); ``"last"`` returns the
        largest (a defect). Raises when the scan finds no crossing.
        """
        f = 1.0 / (self.rtf if rtf is None else rtf)
        n = aa + ac

        def dev(e):
            e = sign * e
            full = np.exp(aa * f * e) - np.exp(-ac * f * e)
            return abs(n * f * e / full - 1.0) - eps

        res = find_roots(dev, 1e-6, 4.0, n_scan=1200, log=True, xtol=1e-16)
        if not res["roots"] or res["failed"]:
            raise RuntimeError(f"linear threshold: {res['message']}")
        return res["roots"][0] if crossing == "first" else res["roots"][-1]


def first_lobe_peak(aa, slope=None, x_hi=124.0, n_scan=4000):
    """Peak of the linear form's overshoot over the exact law, as a fraction.

    For ``aa < 1/2`` on the anodic side the linear form ``slope x`` runs above
    the exact ``e^(aa x) - e^(-(1-aa) x)`` before falling below it, with
    ``x = F eta_s / RT``. The maximum of ``slope x / exact - 1`` is refined by
    a bounded maximisation from a scan of ``x`` in ``[1e-4, x_hi]``; a scan
    maximum on the upper end raises. On the lower end the overshoot is below
    ``|aa - 1/2| 1e-4`` and that end value is returned.

    Parameters
    ----------
    aa : float
        Anodic transfer coefficient, ``alpha_c = 1 - aa``.
    slope : float, optional
        Slope of the linear form; default ``n = 1`` (a defect passes another).
    """
    ac, s = 1.0 - aa, (1.0 if slope is None else slope)

    def g(x):
        # e^(aa x) - e^(-ac x) = e^(-ac x) expm1(x) for aa + ac = 1, no cancellation
        return s * x * np.exp(ac * x) / np.expm1(x) - 1.0

    x = np.geomspace(1e-4, x_hi, n_scan)
    v = g(x)
    k = int(np.argmax(v))
    if k == n_scan - 1:
        raise RuntimeError("first-lobe peak at the upper scan end")
    if k == 0:
        return float(v[0])
    r = minimize_scalar(lambda t: -g(t), bounds=(x[k - 1], x[k + 1]), method="bounded",
                        options={"xatol": 1e-12 * x[k]})
    return float(max(-r.fun, v[k]))


class Cell:
    """The J3.4 cell: lithium foil, separator, porous cathode, one-dimensional.

    Electrolyte concentration, electrolyte potential and solid-surface
    concentration on one cell-centred grid with a face at the separator
    boundary, marched by implicit Euler.

    Parameters
    ----------
    kin : Kinetics
        The parameter set and its kinetic law.
    kappa : float
        Bulk electrolyte conductivity, S/m (``kappa_eff = kappa eps^1.5``).
    n_s, n_c : int
        Cells in the separator and in the cathode.
    nu_factor : float
        The factor of the concentration term of Doyle's Eq. 5 (J3.4 leaves 2).
    p_D : float
        Bruggeman exponent of the diffusivity.
    kinetics : {"eq17", "bv", "linear", "tafel"}
        Law of the porous cathode: Eq. 17 as printed, the standard form with
        Eq. 30's ``i0``, the linear limit, or the cathodic Tafel branch.
    face_mean : {"harmonic", "arithmetic"}
        Face value of the effective transport coefficients; harmonic is the
        series-resistance value at the porosity jump (arithmetic is a defect).
    leak : float
        Salt flux through the current collector as a fraction of ``|I|/F``;
        zero is Doyle's Eq. 12 (anything else is a defect).
    newton_tol : float
        Step tolerance of each implicit-Euler ``newton`` solve.
    """

    def __init__(self, kin, kappa, n_s=24, n_c=48, nu_factor=2.0, p_D=0.5,
                 kinetics="eq17", face_mean="harmonic", leak=0.0, newton_tol=1e-9):
        p = kin.par
        self.newton_tol = newton_tol
        self.kin, self.par, self.kinetics, self.leak = kin, p, kinetics, leak
        self.rtf, self.alpha = kin.rtf, kin.alpha
        self.c_T, self.c_max = kin.c_T, kin.c_max
        self.eps_c, self.dc = p["epsilon"], p["delta_c"]
        self.kappa, self.nu_factor, self.sigma = kappa, nu_factor, p["sigma"]
        self.a_spec = 3.0 * (1 - self.eps_c) / p["R_s"]
        N = n_s + n_c
        self.N, self.n_s = N, n_s
        # one grid, with a face exactly at x = delta_s
        x_f = np.concatenate([np.linspace(0.0, p["delta_s"], n_s + 1),
                              np.linspace(p["delta_s"], p["delta_s"] + self.dc, n_c + 1)[1:]])
        self.x_f, self.x_c = x_f, 0.5 * (x_f[:-1] + x_f[1:])
        self.dx = np.diff(x_f)

        sep = np.arange(N) < n_s
        self.eps = np.where(sep, 1.0, self.eps_c)              # separator is pure polymer
        self.a = np.where(sep, 0.0, self.a_spec)
        self.brug_k = self.eps ** 1.5                          # kappa_eff = kappa eps^1.5
        self.brug_D = self.eps ** p_D                          # D_eff = D eps^0.5 (J3.4)
        w = self.dx

        def to_face(b):
            """Face values of a cell coefficient; harmonic: the porosity jumps."""
            bf = np.empty(N + 1)
            if face_mean == "harmonic":
                bf[1:-1] = (w[:-1] + w[1:]) / (w[:-1] / b[:-1] + w[1:] / b[1:])
            else:
                bf[1:-1] = (w[:-1] * b[:-1] + w[1:] * b[1:]) / (w[:-1] + w[1:])
            bf[0], bf[-1] = b[0], b[-1]
            return bf

        self.kf = self.kappa * to_face(self.brug_k)
        self.Df = p["D"] * to_face(self.brug_D)

        # Outward normal, a dy/dn + b y = d. Homogeneous Neumann at both ends,
        # dy/dn = 0, so the gradient operator returns zero on the boundary faces;
        # the physical boundary fluxes are written on those faces in fluxes():
        # i2 = I at the anode (Eq. 24), i2 = 0 and N = 0 at the collector (Eq. 12).
        bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        self.grad, _ = construct_grad((N,), x_f, self.x_c, bc)
        self.div = construct_div((N,), x_f, nu=0)              # nu=0: Cartesian slab
        # fields (c, Phi_2, c_s) on the last axis; neighbouring cells couple only
        # through the operators, the kinetic source is pointwise in x
        self.numjac = NumJac((N, 3), stencil=stencil_block_diagonals(
            ndims=2, axes_diagonals=[0], axes_blocks=[-1]))
        self.t_full = (1 - self.eps_c) * self.c_T * self.dc * F_CONST   # / I = discharge time

    def t0_plus(self, c):
        """Appendix A of Doyle et al.: transference number, linear in c."""
        return self.par["t0_a"] + self.par["t0_b"] * c

    def unpack(self, y):
        """Concentration, electrolyte potential and solid concentration from the scaled state."""
        y = y.reshape(self.N, 3)
        c = self.par["c_0"] * np.exp(np.clip(y[:, 0], -40.0, 5.0))    # positivity is structural
        return c, y[:, 1] * self.rtf, y[:, 2] * self.c_T

    def fluxes(self, c, phi2, I):
        """Face current density in the electrolyte and face salt flux."""
        gc, glnc, gp = self.grad @ c, self.grad @ np.log(c), self.grad @ phi2
        t0f = np.empty(self.N + 1)
        t0f[1:-1] = 0.5 * (self.t0_plus(c[:-1]) + self.t0_plus(c[1:]))
        t0f[0], t0f[-1] = self.t0_plus(c[0]), self.t0_plus(c[-1])
        i2 = -self.kf * gp + self.nu_factor * self.kf * self.rtf * (1 - t0f) * glnc
        i2[0], i2[-1] = I, 0.0                     # Eq. 24 at the anode; Eq. 12 at the collector
        n_salt = -self.Df * gc - (1 - t0f) * i2 / F_CONST
        n_salt[0] = 0.0                            # no anion into the lithium
        n_salt[-1] = self.leak * abs(I) / F_CONST  # Eq. 12: nothing through the collector
        return i2, n_salt

    def eta_surface(self, c, phi2, cs):
        """Surface overpotential ``eta_s = eta - U_2``, with ``Phi_1`` gauged to zero."""
        css = np.clip(cs, 1.0, self.c_T - 1.0)
        return -phi2 - (self.kin.U_prime(css) + self.rtf * np.log((self.c_T - css) / css))

    def i_transfer(self, c, phi2, cs):
        """Transfer current density of the cathode in the selected law, A/m2."""
        cc = np.clip(c, 1.0, self.c_max - 1.0)
        css = np.clip(cs, 1.0, self.c_T - 1.0)
        al, f = self.alpha, 1.0 / self.rtf
        if self.kinetics == "eq17":
            z = np.clip(al * f * (-phi2 - self.kin.U_prime(css)), -200.0, 200.0)
            pre = F_CONST * self.par["k_2"] * (self.c_max - cc) ** al * cc ** al
            return pre * (css * np.exp(z) - (self.c_T - css) * np.exp(-z))
        i0 = self.kin.i0_eq30(cc, css)
        eta_s = self.eta_surface(cc, phi2, css)
        if self.kinetics == "bv":
            z = np.clip(al * f * eta_s, -200.0, 200.0)
            return i0 * (np.exp(z) - np.exp(-z))
        if self.kinetics == "linear":
            return i0 * 2.0 * al * f * eta_s
        if self.kinetics == "tafel":            # cathodic branch, chosen a priori
            return -i0 * np.exp(np.clip(-al * f * eta_s, -200.0, 200.0))
        raise ValueError(self.kinetics)

    def residual(self, y, y_old, dt, I):
        """Implicit-Euler residual, each field scaled to order one."""
        c, phi2, cs = self.unpack(y)
        c_o, _, cs_o = self.unpack(y_old)
        i2, n_salt = self.fluxes(c, phi2, I)
        ai = self.a * self.i_transfer(c, phi2, cs)
        r_c = self.eps * (c - c_o) / dt + self.div @ n_salt
        r_p = self.div @ i2 - ai
        r_s = np.where(self.a > 0,
                       (1 - self.eps_c) * (cs - cs_o) / dt + ai / F_CONST, cs)
        ts = self.t_full / abs(I)
        return np.stack([r_c * ts / self.par["c_0"], r_p * self.dc / abs(I),
                         r_s * ts / self.c_T], axis=1).ravel()

    def eta_anode(self, c0, I, law=None):
        """Doyle Eqs. 6 to 8 at the lithium foil, in the given law; returns ``(eta_s1, i01)``."""
        law = self.kinetics if law is None else law
        p, rtf = self.par, self.rtf
        i01 = p["i_0_1"] * np.sqrt((self.c_max - c0) * c0 / ((self.c_max - p["c_0"]) * p["c_0"]))
        if law in ("eq17", "bv"):
            return 2 * rtf * np.arcsinh(I / (2 * i01)), i01
        if law == "linear":
            return rtf * I / i01, i01
        if law == "tafel":                       # anodic branch, chosen a priori
            return 2 * rtf * np.log(np.maximum(I / i01, 1e-300)), i01
        raise ValueError(law)

    def cell_potential(self, y, I, anode_law=None):
        """``Phi_1(collector) - Phi_1(anode)`` with ``Phi_1(anode) = 0`` (Eq. 22), V."""
        c, phi2, _ = self.unpack(y)
        eta_s1, _ = self.eta_anode(c[0], I, anode_law)
        k0 = self.kappa * self.brug_k[0]
        dlnc = -(1 - self.t0_plus(c[0])) * I / (F_CONST * self.Df[0] * c[0])
        dphi = (-I + self.nu_factor * k0 * self.rtf * (1 - self.t0_plus(c[0])) * dlnc) / k0
        phi2_0 = phi2[0] - self.x_c[0] * dphi
        return -phi2_0 - eta_s1 - I * self.dc / (2 * self.sigma)

    def utilisation(self, y):
        """Cathode utilisation, Eq. 25."""
        _, _, cs = self.unpack(y)
        m = self.a > 0
        return float(np.sum(cs[m] * self.dx[m]) / (self.c_T * self.dc))

    def salt(self, y):
        """Total salt per unit area, mol/m2."""
        return float(np.sum(self.eps * self.dx * self.unpack(y)[0]))

    def initial(self, u0, I):
        """Uniform-current initial guess, consistent with the chosen law."""
        y = np.zeros((self.N, 3))
        i_avg = -I / (self.a_spec * self.dc)
        i0v = float(self.kin.i0_eq30(self.par["c_0"], u0 * self.c_T))
        al, rtf = self.alpha, self.rtf
        if self.kinetics in ("eq17", "bv"):
            e0 = rtf / al * np.arcsinh(i_avg / (2 * i0v))
        elif self.kinetics == "linear":
            e0 = rtf * i_avg / (2 * al * i0v)
        else:
            e0 = -rtf / al * np.log(max(-i_avg / i0v, 1e-300))
        y[:, 1] = -(self.kin.U_ocp(u0) + e0) / rtf
        y[:, 2] = u0
        return y

    def _step(self, y_old, dt, I):
        """One implicit-Euler step; ``None`` when ``newton`` does not succeed."""
        def fun(z):
            return self.numjac(lambda w: self.residual(w, y_old.ravel(), dt, I), z)
        try:
            sol = newton(fun, y_old.copy(), tol=self.newton_tol, maxfev=30, stagnation=True)
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            return None                          # NewtonBreakdown: the caller shortens dt
        if not sol.success or not np.all(np.isfinite(sol.x)):
            return None
        r = self.residual(sol.x.ravel(), y_old.ravel(), dt, I)
        return sol.x if np.max(np.abs(r)) < 1e-6 else None

    def march(self, I, u0=None, v_stop=1.55, u_stop=0.9995, max_steps=2500,
              dt_seq=None, dt_scale=1.0, scale_dt_max=True, u_probe=None,
              anode_law="bv"):
        """Implicit-Euler discharge at constant current; ``dt`` shrinks on a failed step.

        Parameters
        ----------
        I : float
            Current density, A/m2 (positive on discharge).
        u0 : float, optional
            Initial utilisation; default Doyle's ``u_0``.
        v_stop, u_stop : float
            The stepper stops below ``v_stop`` or above ``u_stop``.
        max_steps : int
            Step cap.
        dt_seq : array_like, optional
            A step sequence recorded from another run, replayed exactly, so two
            kinetic laws are compared on an identical time grid.
        dt_scale : float
            Multiplies the initial step and the step ceiling together (1 is the
            production schedule).
        scale_dt_max : bool
            ``False`` scales the initial step alone (a defect: the ceiling then
            sets the step over almost the whole discharge).
        u_probe : float, optional
            Utilisation at which the cell potential is also computed exactly:
            the run's state is advanced from the last step before the time at
            which Faraday's law puts that utilisation, by one implicit-Euler step
            that lands on it. Stored in ``self.v_probe``.
        anode_law : str
            Law at the lithium foil in the recorded cell potential; ``"bv"`` is
            exact there, so ``kinetics`` acts on the cathode alone.

        Returns
        -------
        out : ndarray
            Per step: ``u, V, t, salt, max|eta_s|, mean eta_s, eta_s1, mean|eta_s|``.
        dts : ndarray
            The steps taken.
        """
        u0 = self.par["u_0"] if u0 is None else u0
        y = self.initial(u0, I)
        t_ref = self.t_full / abs(I)
        dt = 2e-5 * t_ref * dt_scale
        dt_max = 6e-3 * t_ref * (dt_scale if scale_dt_max else 1.0)
        t_probe = None if u_probe is None else (u_probe - u0) * self.t_full / I
        self.v_probe = np.nan
        t, out, dts, dt_peak = 0.0, [], [], dt
        m = self.a > 0
        wm = self.dx[m] / self.dx[m].sum()
        n_max = max_steps if dt_seq is None else len(dt_seq)
        for k in range(n_max):
            step_dt = dt if dt_seq is None else dt_seq[k]
            y_new = None
            for _retry in range(8):
                y_new = self._step(y, step_dt, I)
                if y_new is not None:
                    break
                if dt_seq is not None:
                    break            # a prescribed sequence must not be altered
                step_dt *= 0.4
            if y_new is None:
                break
            if t_probe is not None and t < t_probe <= t + step_dt:
                y_p = self._step(y, t_probe - t, I)
                if y_p is None:
                    raise RuntimeError("the landing step on u_probe did not converge")
                self.v_probe = self.cell_potential(y_p, I, anode_law=anode_law)
            y, t = y_new, t + step_dt
            dts.append(step_dt)
            if dt_seq is None:
                dt = min(step_dt * 1.3, dt_max)
                dt_peak = max(dt_peak, dt)
                if dt < 1e-3 * dt_peak:
                    break    # the step size has collapsed: the model is done
            c, phi2, cs = self.unpack(y)
            es = self.eta_surface(c, phi2, cs)[m]
            u, v = self.utilisation(y), self.cell_potential(y, I, anode_law=anode_law)
            out.append((u, v, t, self.salt(y), float(np.max(np.abs(es))),
                        float(np.sum(wm * es)), float(self.eta_anode(c[0], I)[0]),
                        float(np.sum(wm * np.abs(es)))))
            if dt_seq is None and ((I > 0 and v < v_stop) or u > u_stop or u < 1e-4):
                break
        return np.array(out), np.array(dts)


def u_at_level(u, v, level):
    """Utilisation at which a sampled cell potential crosses ``level``.

    The root of the piecewise-linear interpolant of ``v(u)``: its sign changes
    are exactly those of the samples, so the crossing is unique when exactly
    one sample interval changes sign, and the function raises otherwise.

    Parameters
    ----------
    u, v : array_like
        Utilisation and cell potential per step.
    level : float
        Potential, V.

    Returns
    -------
    float
        The utilisation of the crossing.
    """
    u, g = np.asarray(u, float), np.asarray(v, float) - level
    k = np.flatnonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)
    if len(k) != 1:
        raise RuntimeError(f"{len(k)} crossings of {level} V; expected one")
    k = int(k[0])
    return float(u[k] + (u[k + 1] - u[k]) * g[k] / (g[k] - g[k + 1]))
