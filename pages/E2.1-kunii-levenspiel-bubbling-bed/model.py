"""Kunii and Levenspiel's (1968) bubbling-bed model, in closed form and distributed.

The published model (equations 2 to 57 of the paper, CGS units throughout) is
algebra: the bubbling-bed flow model gives the interchange coefficients and the
solids distribution, and for a first-order reaction equations 48a to 48c reduce
to the continued fraction of equation 49 and the exponential of equation 54.

The distributed model puts back the convection that equations 48b and 48c leave
out: bubble, cloud-and-wake and emulsion gas each move along the bed height with
their own superficial velocity and exchange with their neighbours. It is a linear
system of three first-order ODEs in the height, solved two ways here: by finite
volumes with pymrm operators (``BubblingBed``), and exactly by its eigenmodes
(``extended_exact``).
"""

import numpy as np
from scipy.sparse import bmat, csc_array, identity, kron
from scipy.sparse.linalg import spsolve
from pymrm import compute_boundary_values, construct_convflux_upwind, construct_div
from pymrm.checks import find_roots

G = 980.0                  # cm/s^2
RHO_G, MU = 1.18e-3, 1.8e-4  # g/cm^3, g/(cm s): air, as stated in appendices A and B

#: The four numerical coefficients printed inside equations 10, 12 and 51.
KL = dict(c_ubr=0.711, c_kbc_visc=4.5, c_kbc_diff=5.85, c_kce=6.78)


def hydrodynamics(u0, umf, emf, em, Lm, db, alpha, D, gamma_b=0.0, coef=None):
    """The bubbling-bed flow model: equations 2, 3, 4, 10, 12, 15, 17 and 51.

    Parameters
    ----------
    u0, umf : float
        Superficial and minimum-fluidisation gas velocities [cm/s].
    emf, em : float
        Voidage at minimum fluidisation and of the fixed bed [-].
    Lm : float
        Fixed-bed height [cm].
    db : float
        Effective bubble diameter [cm], the model's one fitted parameter.
    alpha : float
        Wake-to-bubble volume ratio [-].
    D : float
        Gas diffusivity, also used for the emulsion [cm^2/s].
    gamma_b : float, optional
        Solids dispersed in the bubbles per bubble volume [-].
    coef : dict, optional
        Replaces the printed coefficients ``KL`` (used to inject transcription
        defects).

    Returns
    -------
    dict
        ``ubr, ub, delta, Kbc, Kce, fc`` (V_c/V_b + alpha), ``gb, gc, ge``,
        ``Lf`` and the inputs.
    """
    k = KL if coef is None else coef
    ubr = k["c_ubr"] * np.sqrt(G * db)                               # eq. 51
    ub = u0 - umf + ubr                                              # eq. 3
    delta = (u0 - umf) / ub                                          # eq. 2
    Kbc = (k["c_kbc_visc"] * (umf / db)
           + k["c_kbc_diff"] * (D**0.5 * G**0.25 / db**1.25))        # eq. 10
    Kce = k["c_kce"] * np.sqrt(emf * D * ub / db**3)                 # eq. 12
    fc = 3.0 * (umf / emf) / (ubr - umf / emf) + alpha               # eq. 9 + wake
    gc = (1 - emf) * fc                                              # eq. 17
    ge = (1 - emf) * (1 - delta) / delta - gc - gamma_b              # eq. 15
    Lf = Lm * (1 - em) / ((1 - delta) * (1 - emf))                   # eq. 4, twice
    return dict(ubr=ubr, ub=ub, delta=delta, Kbc=Kbc, Kce=Kce, fc=fc,
                gb=gamma_b, gc=gc, ge=ge, Lf=Lf, emf=emf, u0=u0, umf=umf,
                em=em, Lm=Lm, db=db, alpha=alpha, D=D)


def Kr_of_Km(Km, h):
    """Equation 45 read backwards: the rate constant K_r [1/s] behind K_m.

    Parameters
    ----------
    Km : float or numpy.ndarray
        Fixed-bed rate group ``(1 - eps_m) K_r L_m / u_0`` [-].
    h : dict
        Output of ``hydrodynamics``.

    Returns
    -------
    float or numpy.ndarray
        ``K_r`` [1/s].
    """
    return Km * h["u0"] / ((1 - h["em"]) * h["Lm"])


def Kf_over_Km(Km, h):
    """Equations 49 and 50 divided by K_m; finite at K_m = 0.

    Parameters
    ----------
    Km : float or numpy.ndarray
        Fixed-bed rate group [-]; 0 gives the slow-reaction limit exactly.
    h : dict
        Output of ``hydrodynamics``.

    Returns
    -------
    float or numpy.ndarray
        ``K_f / K_m``.
    """
    Kr = Kr_of_Km(Km, h)
    inner = Kr / h["Kce"] + 1.0 / h["ge"]
    mid = h["gc"] + 1.0 / inner
    brack = h["gb"] + 1.0 / (Kr / h["Kbc"] + 1.0 / mid)
    return (1.0 / (1 - h["emf"])) * (h["u0"] / h["ubr"]) * brack


def Kf_closed(Km, h):
    """Equations 49 and 50: the fluidised-bed rate group K_f.

    Parameters
    ----------
    Km : float or numpy.ndarray
        Fixed-bed rate group [-].
    h : dict
        Output of ``hydrodynamics``.

    Returns
    -------
    float or numpy.ndarray
        ``K_f``; the fraction unconverted is ``exp(-K_f)`` (equation 54).
    """
    return Km * Kf_over_Km(Km, h)


def ceiling_Kf(h):
    """The limit of K_f as K_r goes to infinity: ``K_bc L_f / u_b``.

    Parameters
    ----------
    h : dict
        Output of ``hydrodynamics``.

    Returns
    -------
    float
        The largest K_f the bed can reach; ``1 - exp(-ceiling)`` is the largest
        conversion.
    """
    return h["Kbc"] * h["Lf"] / h["ub"]


def backmix_crossing(h, lo=1e-6, hi=1e9, n_scan=400):
    """The K_m at which equation 54 meets backmix flow, ``1/(1 + K_m)``.

    Parameters
    ----------
    h : dict
        Output of ``hydrodynamics``.
    lo, hi : float, optional
        Range scanned on a logarithmic grid.
    n_scan : int, optional
        Scan points for ``pymrm.checks.find_roots``.

    Returns
    -------
    float or None
        The crossing, or None when the range holds no sign change.

    Raises
    ------
    RuntimeError
        When the range holds more than one crossing.
    """
    res = find_roots(lambda K: np.log1p(K) - Kf_closed(K, h), lo, hi,
                     n_scan=n_scan, log=True)
    if len(res["roots"]) > 1 or res["failed"]:
        raise RuntimeError(f"backmix crossing: {res['message']}")
    return res["roots"][0] if res["roots"] else None


def velocities(h, convect=True, ve_mode="balance"):
    """Superficial gas velocities (bubble, cloud and wake, emulsion) [cm/s].

    Parameters
    ----------
    h : dict
        Output of ``hydrodynamics``.
    convect : bool, optional
        False freezes the cloud and the emulsion (equations 48b and 48c).
    ve_mode : {"balance", "paper"}, optional
        ``"balance"``: the emulsion takes ``u_0`` minus the other two.
        ``"paper"``: the emulsion velocity of equations 5 and 7.

    Returns
    -------
    numpy.ndarray
        Shape (3,).
    """
    vb = h["u0"] - h["umf"]                                          # eq. 2
    if not convect:
        return np.array([vb, 0.0, 0.0])
    d, emf = h["delta"], h["emf"]
    vc = d * h["fc"] * emf * h["ub"]          # cloud and wake gas rise with the bubble
    if ve_mode == "balance":
        ve = h["umf"] - vc
    elif ve_mode == "paper":
        us = h["alpha"] * d * h["ub"] / (1 - d - h["alpha"] * d)    # eq. 5
        ve = emf * (1 - d * (1 + h["fc"])) * (h["umf"] / emf - us)  # eq. 7
    else:
        raise ValueError(f"unknown ve_mode {ve_mode!r}")
    return np.array([vb, vc, ve])


def exchange_matrix(h, Kr, kbc_asym=1.0):
    """Exchange and reaction per unit bed volume, the 3 x 3 block of eqs. 48a-c.

    Parameters
    ----------
    h : dict
        Output of ``hydrodynamics``.
    Kr : float
        Rate constant [1/s].
    kbc_asym : float, optional
        Factor on K_bc in the cloud row only; 1 is the model, anything else
        injects an exchange that does not conserve the reactant.

    Returns
    -------
    numpy.ndarray
        ``S`` with ``dF/dz = S C`` (F the upward fluxes).
    """
    d = h["delta"]
    S = np.zeros((3, 3))
    S[0, 0] = -d * (h["gb"] * Kr + h["Kbc"])
    S[0, 1] = d * h["Kbc"]
    S[1, 0] = d * h["Kbc"] * kbc_asym
    S[1, 1] = -d * (h["Kbc"] + h["gc"] * Kr + h["Kce"])
    S[1, 2] = d * h["Kce"]
    S[2, 1] = d * h["Kce"]
    S[2, 2] = -d * (h["Kce"] + h["ge"] * Kr)
    return S


class BubblingBed:
    """Three-region axially distributed bed, finite volumes, one linear solve.

    The state is ``C`` of shape ``(n, 3)``: height first, then (bubble, cloud and
    wake, emulsion). Each region is convected by first-order upwind fluxes at its
    own velocity; a region moving down enters at the top. When one does, the
    inlet values become two extra unknowns: ``C_in`` at the bottom (feed mixed
    with the emulsion gas leaving there) and ``C_top`` at the top (the mixing cup
    of the gas leaving upward), which is the recirculation loop of the paper's
    figure 6.

    Parameters
    ----------
    h : dict
        Output of ``hydrodynamics``.
    n : int
        Number of cells over the fluidised height ``L_f``.
    convect : bool, optional
        False freezes the cloud and the emulsion (the published model).
    ve_mode : {"balance", "paper"}, optional
        Emulsion velocity closure, see ``velocities``.
    outlet : {"outflow", "zero_gradient"}, optional
        ``"outflow"``: the face value at an outlet is the adjacent cell value,
        exact for a convection-only equation. ``"zero_gradient"`` imposes
        ``dC/dn = 0`` there, which the upwind flux folds in as an extrapolation;
        kept to show what it does to the balance.
    """

    def __init__(self, h, n, convect=True, ve_mode="balance", outlet="outflow"):
        self.h, self.n, self.outlet = h, n, outlet
        self.v = velocities(h, convect, ve_mode)
        self.x_f = np.linspace(0.0, h["Lf"], n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.up, self.down = self.v > 0, self.v < 0
        self._build_operators(outlet)

    def _build_operators(self, outlet):
        n, v = self.n, self.v
        div = construct_div((n, 1), self.x_f, nu=0, axis=0)    # nu = 0: constant area
        N = 3 * n
        self.transport = csc_array((N, N))
        self.b_in = np.zeros((n, 3))     # boundary vector per unit C_in
        self.b_top = np.zeros((n, 3))    # boundary vector per unit C_top
        self.bc = [None, None, None]
        for k in range(3):
            if v[k] == 0.0:
                continue                 # a frozen region carries nothing
            # Inlet: C = C_inlet (a = 0, b = 1); d = 1 multiplies the external
            # inlet value through shapes_d. Outlet: material leaves, no condition
            # is physical; "outflow" takes the face value from the last cell.
            inlet = {"a": 0.0, "b": 1.0, "d": 1.0}
            out = ({"outflow": True} if outlet == "outflow"
                   else {"a": 1.0, "b": 0.0, "d": 0.0})       # dC/dn = 0, outward normal
            if v[k] > 0:
                bc, shapes_d = (inlet, out), ((1,), None)     # enters at z = 0
            else:
                bc, shapes_d = (out, inlet), (None, (1,))     # enters at z = L_f
            conv, b0, b1 = construct_convflux_upwind((n, 1), self.x_f, self.x_c, bc=bc,
                                                     v=v[k], axis=0, shapes_d=shapes_d)
            E = np.zeros((3, 3))
            E[k, k] = 1.0
            self.transport = self.transport + kron(div @ conv, csc_array(E), format="csc")
            col = (div @ (b0 if v[k] > 0 else b1)).toarray().ravel()
            (self.b_in if v[k] > 0 else self.b_top)[:, k] = col
            self.bc[k] = bc
        # mixing rows, on the outlet cells (an outflow face carries the cell value)
        self.su = float(self.v[self.up].sum())
        self.row_in = np.zeros((n, 3))
        self.row_in[0, self.down] = -np.abs(v[self.down])
        self.row_top = np.zeros((n, 3))
        self.row_top[-1, self.up] = -v[self.up]

    def solve(self, Km, kbc_asym=1.0):
        """Steady profiles at one K_m: a single sparse linear solve.

        Parameters
        ----------
        Km : float
            Fixed-bed rate group [-].
        kbc_asym : float, optional
            Passed to ``exchange_matrix`` (defect injection).

        Returns
        -------
        BedSolution
        """
        h, n, N = self.h, self.n, 3 * self.n
        Kr = Kr_of_Km(Km, h)
        S = exchange_matrix(h, Kr, kbc_asym)
        M = self.transport - kron(identity(n, format="csc"), csc_array(S), format="csc")
        # residual: div(F) - S C = 0 with F = conv C + b C_inlet
        if self.down.any():
            K = bmat([[M, csc_array(self.b_in.reshape(N, 1)),
                       csc_array(self.b_top.reshape(N, 1))],
                      [csc_array(self.row_in.reshape(1, N)), csc_array([[self.su]]), None],
                      [csc_array(self.row_top.reshape(1, N)), None, csc_array([[self.su]])]],
                     format="csc")
            rhs = np.zeros(N + 2)
            rhs[N] = h["u0"]           # su C_in = u_0 * 1 + sum |v_down| C_bottom
            x = spsolve(K, rhs)
            C, c_in, c_top = x[:N].reshape(n, 3), x[N], x[N + 1]
        else:                          # no recirculation: every inlet carries the feed
            c_in = c_top = 1.0
            C = spsolve(M.tocsc(), -(self.b_in + self.b_top).ravel()).reshape(n, 3)
        return BedSolution(self, Km, Kr, C, c_in, c_top)


class BedSolution:
    """Result of ``BubblingBed.solve``.

    Attributes
    ----------
    C : numpy.ndarray
        Cell concentrations ``(n, 3)``, relative to the feed.
    c_in, c_top : float
        Inlet values at the bottom (upward regions) and the top (downward ones).
    exit_faces : numpy.ndarray
        Outlet face value of each upward region, from ``compute_boundary_values``.
    c_exit : float
        Mixing cup of the gas leaving the top.
    """

    def __init__(self, bed, Km, Kr, C, c_in, c_top):
        self.bed, self.Km, self.Kr, self.C = bed, Km, Kr, C
        self.c_in, self.c_top = c_in, c_top
        v = bed.v
        faces = np.zeros(3)
        for k in np.flatnonzero(bed.up):
            if bed.outlet == "outflow":
                inlet, out = bed.bc[k]
                bc = ({**inlet, "d": c_in}, out)
                _, _, top, _ = compute_boundary_values(C[:, k:k + 1], bed.x_f, bed.x_c, bc)
                faces[k] = float(np.ravel(top)[0])
            else:
                faces[k] = C[-1, k]   # zero-gradient outlet: the last cell, as the mixing rows
        self.exit_faces = faces
        self.c_exit = float(v[bed.up] @ faces[bed.up] / bed.su)

    def consumed(self):
        """Reaction integrated over the bed, per unit cross-section [cm/s x C_feed]."""
        h = self.bed.h
        gam = np.array([h["gb"], h["gc"], h["ge"]])
        dz = np.diff(self.bed.x_f)
        return float(h["delta"] * self.Kr * ((self.C * gam).sum(axis=1) * dz).sum())

    def mass_balance_dev(self):
        """(consumed - feed flow x (1 - c_exit)) / feed flow x (1 - c_exit)."""
        bed = self.bed
        feed = bed.h["u0"] if bed.down.any() else bed.su
        lost = feed * (1.0 - self.c_exit)
        return (self.consumed() - lost) / lost


def extended_exact(Km, h, ve_mode="balance", scaled=True):
    """The distributed model solved exactly through its eigenmodes.

    ``V dC/dz = S C`` with ``V = diag(v)`` has the solution
    ``C(z) = sum_j a_j w_j phi_j(z)``; each ``phi_j`` is normalised to 1 at the end
    where it is largest, so no exponential overflows. The three inlet conditions
    and the two mixing equations give five linear equations for ``a_j``,
    ``C_in`` and ``C_top``.

    Parameters
    ----------
    Km : float
        Fixed-bed rate group [-].
    h : dict
        Output of ``hydrodynamics``.
    ve_mode : {"balance", "paper"}, optional
        Emulsion velocity closure; the emulsion must flow downward.
    scaled : bool, optional
        False writes every mode as ``exp(lam z)`` (defect injection: the
        growing modes then dominate the system by ``exp(lam L)``).

    Returns
    -------
    dict
        ``c_exit``, ``c_in``, ``c_top``, the eigenvalues ``lam`` [1/cm], the
        condition number ``cond`` of the 5 x 5 system and ``imag``, the largest
        imaginary part discarded.
    """
    v = velocities(h, True, ve_mode)
    if not (v[0] > 0 and v[1] > 0 and v[2] < 0):
        raise ValueError("extended_exact needs bubble and cloud up, emulsion down")
    L = h["Lf"]
    lam, W = np.linalg.eig(exchange_matrix(h, Kr_of_Km(Km, h)) / v[:, None])

    def modes(z):
        if not scaled:
            return W * np.exp(lam * z)
        return W * np.where(lam.real <= 0, np.exp(lam * z), np.exp(lam * (z - L)))

    P0, PL = modes(0.0), modes(L)
    su = v[0] + v[1]
    A = np.zeros((5, 5), complex)
    r = np.zeros(5, complex)
    A[0, :3], A[0, 3] = P0[0], -1.0          # bubble enters at the bottom at C_in
    A[1, :3], A[1, 3] = P0[1], -1.0          # cloud and wake likewise
    A[2, :3], A[2, 4] = PL[2], -1.0          # emulsion enters at the top at C_top
    A[3, :3], A[3, 3], r[3] = -abs(v[2]) * P0[2], su, h["u0"]   # bottom mixing
    A[4, :3], A[4, 4] = -(v[0] * PL[0] + v[1] * PL[1]), su       # top mixing cup
    x = np.linalg.solve(A, r)
    return dict(c_exit=x[4].real, c_in=x[3].real, c_top=x[4].real, lam=lam,
                cond=np.linalg.cond(A), imag=float(np.abs(x.imag).max()))


def appendix_A(p, coef=None):
    """Appendix A's chain for gas-solid mass transfer: equations 1, 10, 27, 51.

    Parameters
    ----------
    p : dict
        ``dp, phis, ut, D, y, Sc, umf, emf, gb, db`` as printed (CGS).
    coef : dict, optional
        Replaces ``KL``.

    Returns
    -------
    dict
        The nine quantities the appendix prints, keyed by the dataset's symbols.
    """
    k = KL if coef is None else coef
    Re_t = p["dp"] * RHO_G * p["ut"] / MU
    Sh_t = 2 + 0.6 * p["Sc"]**(1 / 3) * Re_t**0.5                           # eq. 1
    Kbc = (k["c_kbc_visc"] * p["umf"] / p["db"]
           + k["c_kbc_diff"] * (p["D"]**0.5 * G**0.25 / p["db"]**1.25))     # eq. 10
    ubr = k["c_ubr"] * np.sqrt(G * p["db"])                                 # eq. 51
    # delta/(1 - eps_f) = (u_0 - u_mf)/[(1 - eps_mf) u_br], as (1 - delta) u_b = u_br
    brack = p["gb"] * Sh_t + p["y"] * p["phis"] * p["dp"]**2 / (6 * p["D"]) * Kbc
    slope_u0 = brack / ((1 - p["emf"]) * ubr)                               # eq. 27
    Re_per_u0 = p["dp"] * RHO_G / MU
    return {"Re_t": Re_t, "Sh_star_t": Sh_t, "K_bc": Kbc, "ub_minus_u0": ubr - p["umf"],
            "one_minus_ef_times_ub": (1 - p["emf"]) * ubr, "Re_per_u0": Re_per_u0,
            "Sh_slope_u0": slope_u0, "Sh_intercept": -slope_u0 * p["umf"],
            "Sh_slope_Re": slope_u0 / Re_per_u0}


def appendix_B(p, coef=None):
    """Appendix B's chain for gas-solid heat transfer: equations 1, 11, 38, 51.

    Parameters
    ----------
    p : dict
        ``dp, phis, ut, kg, Cpg, Pr, umf, emf, gb, db`` as printed (CGS).
    coef : dict, optional
        Replaces ``KL``.

    Returns
    -------
    dict
        ``Re_t, Nu_star_t, H_bc, ubr, ub_minus_u0, Re_per_u0``, the bracket of
        equation 38 and ``one_minus_ef_own`` = (1 - eps_mf) u_br from the stated
        eps_mf. The apparent Nusselt line is ``brack / x * (ub_minus_u0 + Re /
        Re_per_u0)`` for a chosen ``x`` = (1 - eps_f)(u_0 + 5.8).
    """
    k = KL if coef is None else coef
    Re_t = p["dp"] * RHO_G * p["ut"] / MU
    Nu_t = 2 + 0.6 * p["Pr"]**(1 / 3) * Re_t**0.5                           # eq. 1
    Hbc = (k["c_kbc_visc"] * p["umf"] * RHO_G * p["Cpg"] / p["db"]
           + k["c_kbc_diff"] * (p["kg"] * RHO_G * p["Cpg"])**0.5
           * G**0.25 / p["db"]**1.25)                                       # eq. 11
    ubr = k["c_ubr"] * np.sqrt(G * p["db"])                                 # eq. 51
    brack = p["gb"] * Nu_t + p["phis"] * p["dp"]**2 / (6 * p["kg"]) * Hbc
    return dict(Re_t=Re_t, Nu_star_t=Nu_t, H_bc=Hbc, ubr=ubr, ub_minus_u0=ubr - p["umf"],
                Re_per_u0=p["dp"] * RHO_G / MU, brack=brack,
                one_minus_ef_own=(1 - p["emf"]) * ubr)
