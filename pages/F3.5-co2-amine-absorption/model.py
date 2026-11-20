"""Film and penetration models of CO2 absorption into amine-promoted carbonate (page F3.5).

Bosch, Versteeg and van Swaaij (1989): CO2 dissolves into hot potash that holds
a small amount of amine. Seven species, two kinetic reactions and the
instantaneous proton transfers, all reversible:

    (1) CO2 + OH-      <=> HCO3-               (kinetic, k_OH)
    (2) CO2 + 2 Am     <=> AmCOO- + AmH+       (kinetic, k_2)
    (3) CO3-- + H2O    <=> HCO3- + OH-         (instantaneous, K3 = Kw/Kc2)
    (4) AmH+ + OH-     <=> Am + H2O            (instantaneous, KB = Kw/Kp)

The instantaneous equilibria are fast reversible kinetics with a rate constant
``k_inst``. Every function takes the parameter dictionary ``P`` built by
:func:`params` from the paper's Table 1, so a defect in one constant is one
changed entry.

Contents:

- :func:`bulk`, the well-mixed liquid at a loading, solved in log variables;
- :class:`Film`, the steady stagnant film by finite volumes (pymrm operators,
  ``newton`` with ``NumJac``), continued in ``k_inst`` to the instantaneous limit;
- :class:`Penetration`, the same network as a Higbie penetration model,
  backward Euler with ``newton`` per step on a geometric time grid;
- the closures and limits the page checks: :func:`closures`,
  :func:`ceiling`, :func:`pfo_enhancement`, :func:`vkh`,
  :func:`second_order_dimensionless`.
"""
import numpy as np
from scipy.optimize import root
from scipy.sparse import diags
from pymrm import NumJac, construct_div, construct_grad, newton
from pymrm.checks import find_roots

#: Species order of the film; the unpromoted film uses the first four.
SP = ["CO2", "HCO3", "CO3", "OH", "Am", "AmH", "AmCOO"]
#: Charge, carbon and amine count per species.
CHARGE = {"CO2": 0, "HCO3": -1, "CO3": -2, "OH": -1, "Am": 0, "AmH": +1, "AmCOO": -1}
CARBON = {"CO2": 1, "HCO3": 1, "CO3": 1, "OH": 0, "Am": 0, "AmH": 0, "AmCOO": 1}
AMINE = {"CO2": 0, "HCO3": 0, "CO3": 0, "OH": 0, "Am": 1, "AmH": 1, "AmCOO": 1}
#: Reaction columns of the stoichiometry matrix.
R1, R3, R2, R4 = 0, 1, 2, 3
#: The continuation in the instantaneous rate constant, 1/s (see :meth:`Film.solve`).
K_LADDER = (1e2, 1e3, 1e4, 1e5, 1e6, 1e7)


def params(table1, R_GAS=8.314, T=363.0):
    """The parameter dictionary from the paper's Table 1.

    Parameters
    ----------
    table1 : pandas.DataFrame
        Columns ``parameter``, ``value_DEA``, ``value_HDA``, as transcribed.
    R_GAS, T : float
        Gas constant, J/(mol K), and temperature, K (stated in the paper).

    Returns
    -------
    dict
        Scalars ``C_AM, C_CARB`` (mol/m3), ``KL`` (m/s), ``M_DIST`` (-),
        ``KG`` (m/s), ``KC1`` (m3/mol), ``KC2, KW`` (mol/m3, mol2/m6),
        ``KOH`` (m3/(mol s)), ``DCO2, DAM, DION`` (m2/s); per amine ``KP``,
        ``KCC``, ``K2``; ``D``, the diffusivity of each species.
    """
    def p(name, amine="DEA"):
        return float(table1.loc[table1.parameter == name, f"value_{amine}"].iloc[0])

    P = dict(R_GAS=R_GAS, T=T, C_AM=p("c_Am"), C_CARB=p("c_carb"), KL=p("k_L"),
             M_DIST=p("m"), KG=p("k_g"), KC1=p("K_c1"), KC2=p("K_c2"), KW=p("K_w"),
             KOH=p("k_OH"), DCO2=p("D_CO2"), DAM=p("D_Am"), DION=p("D_ion"))
    for key, name in (("KP", "K_p"), ("KCC", "K_c"), ("K2", "k_2")):
        P[key] = {a: p(name, a) for a in ("DEA", "HDA")}
    P["D"] = {"CO2": P["DCO2"], "HCO3": P["DION"], "CO3": P["DION"], "OH": P["DION"],
              "Am": P["DAM"], "AmH": P["DION"], "AmCOO": P["DION"]}
    return P


def c_interface(P, p_kpa):
    """Interfacial dissolved CO2, ``c_i = m p/(RT)``, mol/m3.

    Parameters
    ----------
    P : dict
        Parameters (:func:`params`).
    p_kpa : float
        CO2 partial pressure, kPa.

    Returns
    -------
    float
    """
    return P["M_DIST"] * p_kpa * 1e3 / (P["R_GAS"] * P["T"])


def bulk(P, alpha, amine=None, k3_scale=1.0):
    """Equilibrium bulk composition at loading ``alpha`` (mol CO2 per mol K2CO3).

    The equilibria of (1), (3), (4) and the carbamate, with total carbon
    ``c_carb (1 + alpha)``, total amine and electroneutrality against
    ``[K+] = 2 c_carb``, solved in log-concentrations.

    Parameters
    ----------
    P : dict
        Parameters (:func:`params`).
    alpha : float
        Loading.
    amine : {None, "DEA", "HDA"}
        Promoter, or none.
    k3_scale : float
        ``K3 = Kw/Kc2`` is divided by it (1 = as printed): the reconstruction
        scalar, a multiplication of ``Kc2``.

    Returns
    -------
    dict
        Concentration of every species in :data:`SP`, mol/m3, and ``res``, the
        largest residual of the equations.

    Raises
    ------
    RuntimeError
        If the solve fails or its residual exceeds 1e-9.
    """
    k3 = P["KW"] / P["KC2"] / k3_scale
    c_carb, kc1 = P["C_CARB"], P["KC1"]
    carbon, kplus, xi = c_carb * (1.0 + alpha), 2.0 * c_carb, alpha * c_carb
    if amine is None:
        def eqs(logy):
            co2, hco3, co3, oh = np.exp(np.clip(logy, -80, 80))
            return [np.log(hco3) - np.log(kc1 * co2 * oh),          # equilibrium (1)
                    np.log(hco3 * oh) - np.log(k3 * co3),           # equilibrium (3)
                    (co2 + hco3 + co3 - carbon) / c_carb,           # carbon
                    (kplus - (hco3 + 2 * co3 + oh)) / c_carb]       # charge
        guess = np.log([0.3, max(2 * xi, 1.0), max(c_carb - xi, 1.0), 5.0])
    else:
        kc, kb, c_am = P["KCC"][amine], P["KW"] / P["KP"][amine], P["C_AM"]

        def eqs(logy):
            co2, hco3, co3, oh, am, amh, amcoo = np.exp(np.clip(logy, -80, 80))
            return [np.log(hco3) - np.log(kc1 * co2 * oh),
                    np.log(hco3 * oh) - np.log(k3 * co3),
                    np.log(amh * oh) - np.log(kb * am),             # equilibrium (4)
                    np.log(amcoo) - np.log(kc * am * hco3),         # carbamate
                    (am + amh + amcoo - c_am) / c_am,               # amine total
                    (co2 + hco3 + co3 + amcoo - carbon) / c_carb,
                    (kplus + amh - (hco3 + 2 * co3 + amcoo + oh)) / c_carb]
        guess = np.log([0.3, max(1.5 * xi, 1.0), c_carb - 0.6 * xi, 5.0,
                        0.6 * c_am, 0.05 * c_am, 0.3 * c_am])
    sol = root(eqs, guess, method="hybr", options={"xtol": 1e-13})
    resid = float(np.max(np.abs(eqs(sol.x))))
    if not (sol.success and resid < 1e-9):
        raise RuntimeError(f"bulk speciation failed at alpha={alpha}, {amine}: "
                           f"residual {resid:.1e}")
    y = np.exp(sol.x)
    out = (dict(zip(SP, y)) if amine
           else {**dict(zip(SP[:4], y)), "Am": 0.0, "AmH": 0.0, "AmCOO": 0.0})
    out["res"] = resid
    return out


def stretched_faces(n, s=6.0):
    """Face positions on ``[0, 1]``, clustered at the interface ``xi = 0``.

    Parameters
    ----------
    n : int
        Number of cells.
    s : float
        Grading of ``x_f = (exp(s t) - 1)/(exp(s) - 1)``.

    Returns
    -------
    numpy.ndarray
        ``n + 1`` face positions.
    """
    u = np.linspace(0.0, 1.0, n + 1)
    return (np.exp(s * u) - 1.0) / (np.exp(s) - 1.0)


class Film:
    """Stagnant film with the reactions (1) to (4), steady, on ``xi = x/L``.

    The residual is ``A c + g + s(c) = 0`` with ``A = div(-grad)`` and
    ``s_j = -(L^2/D_j) sum_r nu_jr R_r``, so a consumption enters with a
    positive sign.

    Parameters
    ----------
    P : dict
        Parameters (:func:`params`).
    cb : dict
        Bulk concentrations, mol/m3 (:func:`bulk`).
    c_gas : float
        CO2 concentration in the gas, ``p/(RT)``, mol/m3.
    amine : {None, "DEA", "HDA"}
        Promoter, or none (then four species).
    n_x : int
        Number of cells.
    s_grid : float
        Grading towards the interface (:func:`stretched_faces`).
    k3_scale : float
        Divides ``K3`` (the reconstruction scalar).
    keq_scale : float
        Multiplies ``Keq`` of reaction (2).
    enable_r1, enable_r34 : bool
        Reaction (1); the instantaneous reactions (3) and (4).
    irreversible2 : bool
        Reaction (2) without its reverse term.
    length : float, optional
        Domain length, m; the film thickness ``delta = D_CO2/k_L`` by default.
    pfo_k1 : float, optional
        Pseudo-first-order limit: reaction (1) becomes ``pfo_k1 * [CO2]``,
        irreversible, with OH- buffered (its column of the stoichiometry zero).
    nu_set : dict, optional
        ``{(species, reaction column): coefficient}`` overriding the
        stoichiometry matrix.
    """

    def __init__(self, P, cb, c_gas, amine=None, n_x=400, s_grid=6.0, k3_scale=1.0,
                 keq_scale=1.0, enable_r1=True, enable_r34=True, irreversible2=False,
                 length=None, pfo_k1=None, nu_set=None):
        self.P, self.amine = P, amine
        self.k3 = P["KW"] / P["KC2"] / k3_scale
        self.enable_r1, self.enable_r34, self.irrev2 = enable_r1, enable_r34, irreversible2
        self.pfo_k1 = pfo_k1
        self.species = SP if amine is not None else SP[:4]
        self.idx = {sp: i for i, sp in enumerate(self.species)}
        self.n_c, self.n_x = len(self.species), n_x
        self.diff = np.array([P["D"][sp] for sp in self.species])
        self.L = P["DCO2"] / P["KL"] if length is None else length
        if amine is not None:
            self.k2, self.kc = P["K2"][amine], P["KCC"][amine]
            self.kb = P["KW"] / P["KP"][amine]
            # thermodynamic cycle (2) = (1) + carbamate + (4) reversed:
            # Keq = Kc1 Kc KB = Kc Kc1 Kw/Kp, the paper's definition
            self.keq = P["KC1"] * self.kc * self.kb * keq_scale
        self.cb = np.array([cb[sp] for sp in self.species])
        self.c_gas, self.k_inst = c_gas, K_LADDER[-1]
        self.c_i = P["M_DIST"] * c_gas
        self.shape = (n_x, self.n_c)                 # spatial axis first, then species
        self.x_f = stretched_faces(n_x, s_grid)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # a dc/dn + b c = d with n the OUTWARD normal (n = -xi at xi = 0, +xi at 1).
        # xi = 0, the interface: CO2 at its interfacial value, c = m c_gas (the
        # gas-side resistance is negligible, measured on the page); every other
        # species -dc/dxi = 0, non-volatile, no flux through the interface.
        a0, b0, d0 = np.ones(self.n_c), np.zeros(self.n_c), np.zeros(self.n_c)
        a0[0], b0[0], d0[0] = 0.0, 1.0, self.c_i
        # xi = 1, the bulk edge: c_j = c_bulk,j for every species.
        self.bc = ({"a": [a0.tolist()], "b": [b0.tolist()], "d": [d0.tolist()]},
                   {"a": [[0.0] * self.n_c], "b": [[1.0] * self.n_c],
                    "d": [self.cb.tolist()]})
        grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc, axis=0)
        div_mat = construct_div(self.shape, self.x_f, nu=0, axis=0)   # slab
        # the constant diffusion operator, assembled once; flux_mat c + flux_bc is
        # -dc/dxi on every face, boundary terms included
        self.flux_mat, self.flux_bc = -grad_mat, -grad_bc
        self.jac_diff = div_mat @ self.flux_mat
        self.g_diff = np.asarray((div_mat @ self.flux_bc).todense()).ravel()
        # the rate is pointwise: the default stencil couples the species axis only
        self.numjac = NumJac(self.shape)
        n_r = 4 if amine is not None else 2
        NU = np.zeros((self.n_c, n_r))
        i = self.idx
        NU[i["CO2"], R1], NU[i["OH"], R1], NU[i["HCO3"], R1] = -1, -1, 1
        NU[i["CO3"], R3], NU[i["HCO3"], R3], NU[i["OH"], R3] = -1, 1, 1
        if amine is not None:
            NU[i["CO2"], R2], NU[i["Am"], R2] = -1, -2
            NU[i["AmCOO"], R2], NU[i["AmH"], R2] = 1, 1
            NU[i["AmH"], R4], NU[i["OH"], R4], NU[i["Am"], R4] = -1, -1, 1
        if pfo_k1 is not None:
            NU[i["OH"], R1] = 0.0                    # OH- buffered
        for (sp, r), v in (nu_set or {}).items():
            NU[i[sp], r] = v
        self.NU = NU

    def rates(self, c):
        """Reaction rates, mol/(m3 s), one column per reaction.

        Parameters
        ----------
        c : numpy.ndarray
            Concentrations, shape ``(n_x, n_c)``.

        Returns
        -------
        numpy.ndarray
        """
        i, P = self.idx, self.P
        co2, hco3, co3, oh = (c[:, i["CO2"]], c[:, i["HCO3"]], c[:, i["CO3"]],
                              c[:, i["OH"]])
        r = np.zeros((c.shape[0], self.NU.shape[1]))
        if self.pfo_k1 is not None:
            r[:, R1] = self.pfo_k1 * co2
        elif self.enable_r1:
            r[:, R1] = P["KOH"] * (co2 * oh - hco3 / P["KC1"])
        if self.enable_r34:
            r[:, R3] = self.k_inst * (co3 - hco3 * oh / self.k3)
        if self.amine is not None:
            am, amh, amcoo = c[:, i["Am"]], c[:, i["AmH"]], c[:, i["AmCOO"]]
            rev = 0.0 if self.irrev2 else amcoo * amh / (self.keq * np.maximum(am, 1e-12))
            r[:, R2] = self.k2 * (co2 * am - rev)
            if self.enable_r34:
                r[:, R4] = self.k_inst * (amh * oh - self.kb * am)
        return r

    def source(self, c):
        """``-(L^2/D_j) sum_r nu_jr R_r``, the sign matching ``A = div(-grad)``."""
        return -(self.rates(c) @ self.NU.T) * (self.L**2 / self.diff)

    def residual(self, c):
        """Residual and Jacobian of the steady film."""
        g_r, jac_r = self.numjac(self.source, c.reshape(self.shape))
        g = (self.jac_diff @ c.reshape((-1, 1)) + self.g_diff.reshape((-1, 1))
             + g_r.reshape((-1, 1)))
        return g, self.jac_diff + jac_r

    def solve(self, k_ladder=K_LADDER, tol=None, bounds=(1e-12, np.inf)):
        """Steady film by Newton, continued along a ladder of ``k_inst``.

        The instantaneous equilibria are stiff at ``k_inst = 1e7`` 1/s: from the
        flat bulk start Newton does not reach that solution directly, so each
        rung starts from the converged profile of the one below. The page
        measures that a further decade moves ``E`` far below its tolerance.

        Parameters
        ----------
        k_ladder : sequence of float
            The rungs of ``k_inst``, 1/s; the last one is the model.
        tol : float, optional
            Newton step tolerance, mol/m3. The round-off floor of the linear
            solves grows with the inverse square of the first cell, so with
            ``n_x^2`` on this grid; the default ``1e-8 (n_x/400)^2`` sits above it.
        bounds : tuple or None
            Newton bounds, mol/m3: positivity by default. ``None`` for a linear
            rate whose solution falls below any positive floor (the
            pseudo-first-order limit, ``exp(-Ha/2)`` mid-film).

        Returns
        -------
        numpy.ndarray
            Concentrations, shape ``(n_x, n_c)``.

        Raises
        ------
        RuntimeError
            If Newton does not converge on a rung.
        """
        if tol is None:
            tol = 1e-8 * (self.n_x / 400.0) ** 2
        c = np.tile(np.maximum(self.cb, 1e-12), (self.n_x, 1))
        for k_i in k_ladder:
            self.k_inst = k_i
            res = newton(self.residual, c, tol=tol, maxfev=200, bounds=bounds,
                         stagnation=True)
            if not res.success:
                raise RuntimeError(f"film newton failed at k_inst={k_i}: {res.message}")
            c = np.asarray(res.x).reshape(self.shape)
        self.c = c
        return c

    def face_fluxes(self):
        """Species fluxes along ``+x`` (into the liquid), mol/(m2 s), at every face.

        From the face-flux operator with its boundary terms, so the flux the
        solve transports; shape ``(n_x + 1, n_c)``.
        """
        f = (self.flux_mat @ self.c.reshape((-1, 1))).ravel() \
            + np.asarray(self.flux_bc.todense()).ravel()
        return f.reshape((self.n_x + 1, self.n_c)) * self.diff / self.L

    def J_co2(self):
        """Absorption flux of CO2 through the interface face, mol/(m2 s)."""
        return float(self.face_fluxes()[0, 0])

    def enhancement(self):
        """``E = J/(k_L (c_i - c_b))``."""
        return self.J_co2() / (self.P["KL"] * (self.c_i - self.cb[0]))


def run_case(P, alpha, p_kpa, amine, k3_scale=1.0, solve_kw=None, **kw):
    """Bulk speciation and a solved film at one condition.

    Parameters
    ----------
    P : dict
        Parameters (:func:`params`).
    alpha, p_kpa : float
        Loading and CO2 partial pressure, kPa.
    amine : {None, "DEA", "HDA"}
    k3_scale : float
        The reconstruction scalar.
    solve_kw : dict, optional
        Passed to :meth:`Film.solve`.
    **kw
        Passed to :class:`Film`.

    Returns
    -------
    Film
        Solved, with ``cb_dict`` the bulk composition.
    """
    cb = bulk(P, alpha, amine, k3_scale=k3_scale)
    film = Film(P, cb, p_kpa * 1e3 / (P["R_GAS"] * P["T"]), amine, k3_scale=k3_scale, **kw)
    film.solve(**(solve_kw or {}))
    film.cb_dict = cb
    return film


def contact_time(P):
    """Higbie contact time with the same ``k_L``: ``theta = 4 D/(pi k_L^2)``, s."""
    return 4.0 * P["DCO2"] / (np.pi * P["KL"] ** 2)


class Penetration(Film):
    """The film's network as a Higbie penetration model, transient.

    A semi-infinite liquid, truncated at ``fac sqrt(D_CO2 theta)``, starts at
    the bulk composition and is exposed for the contact time ``theta``.

    Parameters
    ----------
    P, cb, c_gas, amine
        As :class:`Film`.
    fac : float
        Domain depth in penetration lengths.
    theta_factor : float
        Multiplies the contact time (1: ``4 D/(pi k_L^2)``).
    **kw
        Passed to :class:`Film`.
    """

    def __init__(self, P, cb, c_gas, amine=None, fac=12.0, theta_factor=1.0, **kw):
        self.theta = theta_factor * contact_time(P)
        super().__init__(P, cb, c_gas, amine,
                         length=fac * np.sqrt(P["DCO2"] * self.theta), **kw)

    def solve_transient(self, n_t=640, t0_frac=1e-5, tol_rel=1e-9, uniform=False):
        """Mean absorption flux over the contact time.

        Backward Euler on the geometric time grid ``theta * geomspace(t0_frac,
        1, n_t)`` after ``t = 0``, ``newton`` per step from the previous level.
        The flux on ``[0, t_1]`` is the exact ``1/sqrt(t)`` integral of physical
        absorption; the rest is the trapezoid rule over the interface fluxes.

        Parameters
        ----------
        n_t : int
            Number of time levels.
        t0_frac : float
            First level as a fraction of ``theta``.
        tol_rel : float
            Newton step tolerance per step, as a fraction of the largest
            concentration (bulk or interface).
        uniform : bool
            Equal time steps from ``t = 0`` instead of the geometric grid.

        Returns
        -------
        float
            Mean CO2 flux, mol/(m2 s).

        Raises
        ------
        RuntimeError
            If Newton fails on a step.
        """
        theta = self.theta
        self.k_inst = K_LADDER[-1]
        tol = tol_rel * max(1.0, self.cb.max(), self.c_i)
        tgrid = (np.linspace(0.0, theta, n_t + 1) if uniform else
                 np.concatenate([[0.0], theta * np.geomspace(t0_frac, 1.0, n_t)]))
        c = np.tile(np.maximum(self.cb, 1e-12), (self.n_x, 1))
        tau = self.L ** 2 / self.diff
        j_hist = np.empty(n_t + 1)
        j_hist[0] = np.nan
        fbc0 = float(np.asarray(self.flux_bc.todense()).ravel()[0])
        for n in range(1, n_t + 1):
            dt = tgrid[n] - tgrid[n - 1]
            w = np.tile(tau / dt, self.n_x)
            W = diags(w, format="csc")
            c_old = c.ravel().copy()

            def residual(u):
                g, jac = self.residual(u)
                return g + (w * (u.ravel() - c_old)).reshape((-1, 1)), jac + W

            res = newton(residual, c, tol=tol, maxfev=30, bounds=(1e-12, np.inf),
                         stagnation=True)
            if not res.success:
                raise RuntimeError(f"penetration step {n} of {n_t} failed: {res.message}")
            c = np.asarray(res.x).reshape(self.shape)
            # interface flux from row 0 of the face-flux operator (face 0, CO2)
            j_hist[n] = ((self.flux_mat @ c.reshape((-1, 1))).ravel()[0] + fbc0) \
                * self.diff[0] / self.L
        j_int = 2.0 * (self.c_i - self.cb[0]) * np.sqrt(self.diff[0] * tgrid[1] / np.pi)
        j_int += np.trapezoid(j_hist[1:], tgrid[1:])
        self.c = c
        self.J_avg = j_int / theta
        return self.J_avg

    def enhancement_avg(self):
        """Mean enhancement, ``J_avg/(k_L (c_i - c_b))``."""
        return self.J_avg / (self.P["KL"] * (self.c_i - self.cb[0]))


def closures(film):
    """Conservation checks on a solved film, each divided by the carbon flux.

    Returns
    -------
    dict
        ``carbon`` (spread of the carbon flux over the faces), ``amine`` (largest
        amine flux), ``charge`` (largest charge flux) and ``electroneutrality``
        (largest ``|2 c_carb + sum z_j c_j|/(2 c_carb)`` over the cells).
    """
    fl = film.face_fluxes()
    sp = film.species
    carbon = fl @ np.array([CARBON[s] for s in sp], float)
    charge_v = np.array([CHARGE[s] for s in sp], float)
    amine = fl @ np.array([AMINE[s] for s in sp], float)
    j0 = abs(carbon[0])
    k_plus = 2.0 * film.P["C_CARB"]
    return {"carbon": float(np.max(np.abs(carbon - carbon[0])) / j0),
            "amine": float(np.max(np.abs(amine)) / j0),
            "charge": float(np.max(np.abs(fl @ charge_v)) / j0),
            "electroneutrality": float(np.max(np.abs(k_plus + film.c @ charge_v)) / k_plus)}


def ceiling(film, with_carbonate=True):
    """Instantaneous-reaction ceiling of ``E`` from the transportable base capacity.

    ``1 + (D_ion ([OH-] + [CO3--]) + D_Am [Am]/2)/(D_CO2 (c_i - c_b))`` at the
    bulk composition; ``with_carbonate=False`` drops the carbonate's capacity.
    """
    P, cbd = film.P, film.cb_dict
    drive = P["DCO2"] * (film.c_i - cbd["CO2"])
    cap = P["DION"] * (cbd["OH"] + (cbd["CO3"] if with_carbonate else 0.0)) + (
        P["DAM"] * cbd["Am"] / 2.0 if film.amine else 0.0)
    return 1.0 + cap / drive


def pfo_enhancement(ha, cb_ci):
    """Closed-form film enhancement of a pseudo-first-order reaction.

    ``E = Ha (cosh Ha - c_b/c_i)/(sinh Ha (1 - c_b/c_i))``, with a bulk
    concentration ``c_b`` at equilibrium.
    """
    return ha * (np.cosh(ha) - cb_ci) / (np.sinh(ha) * (1.0 - cb_ci))


def vkh(ha, e_inf):
    """Van Krevelen-Hoftijzer enhancement, ``E = b/tanh b``, ``b = Ha sqrt((E_inf - E)/(E_inf - 1))``.

    Raises
    ------
    RuntimeError
        Unless :func:`pymrm.checks.find_roots` returns exactly one root on
        ``[1, E_inf)``.
    """
    def f(E):
        b = ha * np.sqrt(max((e_inf - E) / (e_inf - 1.0), 1e-30))
        return E - b / np.tanh(b)

    r = find_roots(f, 1.0, e_inf * (1.0 - 1e-12), xtol=1e-13)
    if len(r["roots"]) != 1:
        raise RuntimeError(f"VKH: {len(r['roots'])} roots ({r['message']})")
    return float(r["roots"][0])


def second_order_dimensionless(ha, q, n_x=400, s_grid=6.0, tol=1e-11):
    """Irreversible second-order film, written dimensionless (two species).

    ``a'' = Ha^2 a b``, ``b'' = (Ha^2/q) a b`` on ``[0, 1]``, ``a = CO2/c_i``,
    ``b = Am/Am_b``, ``q = E_inf - 1``; an assembly written apart from
    :class:`Film`, on the same grid.

    Returns
    -------
    float
        ``E = -da/dxi`` on the interface face.

    Raises
    ------
    RuntimeError
        If Newton fails.
    """
    x_f = stretched_faces(n_x, s_grid)
    x_c = 0.5 * (x_f[:-1] + x_f[1:])
    # a dc/dn + b c = d, outward normal. xi = 0: a = 1 (saturated), -db/dxi = 0
    # (amine non-volatile). xi = 1: a = 0, b = 1 (bulk).
    bc = ({"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[1.0, 0.0]]},
          {"a": [[0.0, 0.0]], "b": [[1.0, 1.0]], "d": [[0.0, 1.0]]})
    grad, gbc = construct_grad((n_x, 2), x_f, x_c, bc, axis=0)
    div = construct_div((n_x, 2), x_f, nu=0, axis=0)
    A = div @ -grad
    g = np.asarray((div @ -gbc).todense()).reshape((-1, 1))
    nj = NumJac((n_x, 2))

    def rate(u):
        r = np.clip(u[:, 0], 0, None) * np.clip(u[:, 1], 0, None)
        return np.stack([ha**2 * r, (ha**2 / q) * r], axis=1)

    def res(u):
        gr, jr = nj(rate, u.reshape((n_x, 2)))
        return A @ u.reshape((-1, 1)) + g + gr.reshape((-1, 1)), A + jr

    sol = newton(res, np.tile([0.5, 1.0], (n_x, 1)), tol=tol, maxfev=100,
                 bounds=(0.0, np.inf), stagnation=True)
    if not sol.success:
        raise RuntimeError(f"dimensionless second-order film failed: {sol.message}")
    u = np.asarray(sol.x).reshape((-1, 1))
    return float((-grad[[0], :] @ u).ravel()[0] - gbc.toarray().ravel()[0])
