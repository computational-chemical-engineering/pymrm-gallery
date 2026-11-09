"""Models of page A3.7 (van 't Riet 1979).

The two printed k_L a correlations and the gas column behind the gassing-out
method that produced most of their data:

* :func:`kla`, :func:`crossing_pv`, :func:`refit_cost`, :func:`log_box_fraction`,
  the branched power law (eqs. 8 and 9) and the closed-form arithmetic done on it;
* :func:`phi_closed` and :func:`phi_pymrm`, the quasi-steady gas column, in
  closed form and with the pymrm operators;
* :func:`gassing_out`, :func:`crossing_time`, :func:`fit_eq2` and
  :func:`lagrangian`, the transient experiment, eq. 2 applied to it, and an
  operator-free route along the gas characteristics.

``CORR`` holds the printed constants of eqs. 8 and 9 (the page checks them
against its CSV). Units are SI unless a docstring says otherwise.
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve, splu
from pymrm import construct_div, construct_convflux_upwind, compute_boundary_values
from pymrm.checks import find_roots

# --- the two branches, exactly as printed (eqs. 8 and 9) ---------------
CORR = {
    "pure_water": dict(eq=8, K=2.6e-2, alpha=0.4, beta=0.5),
    "ionic":      dict(eq=9, K=2.0e-3, alpha=0.7, beta=0.2),
}

def kla(pv, vs, branch=None, K=None, alpha=None, beta=None):
    """k_L*A [1/s] from P/V [W/m3] and v_s [m/s]. Eq. 6 with the branch's constants.

    Parameters
    ----------
    pv : array_like
        Power per unit volume, W/m3.
    vs : array_like
        Superficial gas velocity, m/s.
    branch : {"pure_water", "ionic"}, optional
        Key of ``CORR``; when given it overrides ``K``, ``alpha`` and ``beta``.
    K, alpha, beta : float, optional
        Prefactor and exponents, used when ``branch`` is None.

    Returns
    -------
    numpy.ndarray
        ``K * pv**alpha * vs**beta``, in 1/s.
    """
    if branch is not None:
        c = CORR[branch]
        K, alpha, beta = c["K"], c["alpha"], c["beta"]
    return K * np.asarray(pv, float) ** alpha * np.asarray(vs, float) ** beta

def crossing_pv(vs_pure, vs_ionic, C=CORR):
    """P/V at which the pure-water curve at vs_pure meets the ionic curve at vs_ionic.

    Setting eq. 8 = eq. 9 and solving: the exponents differ, so there is exactly
    one positive root and it is available in closed form.

    Parameters
    ----------
    vs_pure, vs_ionic : float
        Superficial gas velocities of the pure-water and ionic curves, m/s.
    C : dict, optional
        Constants in the layout of ``CORR``.

    Returns
    -------
    float
        The crossing P/V, W/m3.
    """
    p, i = C["pure_water"], C["ionic"]
    num = p["K"] * vs_pure ** p["beta"]
    den = i["K"] * vs_ionic ** i["beta"]
    return (num / den) ** (1.0 / (i["alpha"] - p["alpha"]))

def refit_cost(exp_alt, exp_ref, lo, hi, which):
    """How much a different exponent costs, once the prefactor is re-fitted.

    The prefactor is re-fitted so the altered correlation is EXACT at the
    geometric centre of the printed window. What is returned is the largest
    fractional departure from the printed correlation anywhere in that window -
    i.e. the most a reader could be misled by using `exp_alt` instead.

    This is the quantity van 't Riet asks for and does not compute: "Relatively
    much store is set by the values of alpha and beta. To what extent these can
    be varied without affecting the outcome significantly ...".

    Parameters
    ----------
    exp_alt, exp_ref : float
        The altered and the printed exponent.
    lo, hi : float
        Ends of the window over which the exponent's variable is compared.
    which : str
        Label of the variable (``"pv"`` or ``"vs"``); not used in the arithmetic.

    Returns
    -------
    cost : float
        Largest fractional departure from the printed correlation in the window.
    factors : numpy.ndarray
        The two end-point ratios ``(end / centre)**(exp_alt - exp_ref)``.
    """
    centre = np.sqrt(lo * hi)
    ends = np.array([lo, hi]) / centre
    factors = ends ** (exp_alt - exp_ref)
    return float(np.max(np.abs(factors - 1.0))), factors

def log_box_fraction(a, b, c, x_box, y_box):
    """Fraction of a log-log box where ``a ln x + b ln y > c``, in closed form.

    The region is a half-plane in ``(ln x, ln y)``, so its share of the box is
    the area of a convex polygon: the box clipped by one straight line. This is
    the area fraction that a uniform geometric grid over the box converges to.

    Parameters
    ----------
    a, b, c : float
        Coefficients of the half-plane in logarithmic coordinates.
    x_box, y_box : tuple of float
        Ends of the box in ``x`` and ``y`` (positive).

    Returns
    -------
    float
        Area fraction in ``(ln x, ln y)``, between 0 and 1.
    """
    X, Y = np.log(x_box), np.log(y_box)
    poly = [(X[0], Y[0]), (X[1], Y[0]), (X[1], Y[1]), (X[0], Y[1])]
    s = lambda p: a * p[0] + b * p[1] - c
    out = []
    for p, q in zip(poly, poly[1:] + poly[:1]):          # Sutherland-Hodgman, one edge
        if s(p) > 0:
            out.append(p)
        if (s(p) > 0) != (s(q) > 0):
            t = s(p) / (s(p) - s(q))
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    if len(out) < 3:
        return 0.0
    xs, ys = np.array(out).T
    area = 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1)))
    return float(area / ((X[1] - X[0]) * (Y[1] - Y[0])))

# ---------------------------------------------------------------------
# Gas plug flow up the vessel, liquid perfectly mixed.
#
#   gas    (per unit DISPERSION volume, holdup eps, superficial velocity v_s):
#          eps dC_G/dt + d(v_s C_G)/dz = -(1-eps) k (C_G/H - C_L)
#   liquid (per unit LIQUID volume, so that it reduces to the review's eq. 2):
#          dC_L/dt = k (<C_G>/H - C_L),        <C_G> = (1/H_t) int_0^{H_t} C_G dz
#
# k is k_L*A on a liquid-volume basis, matching eq. 2. H is the Henry constant
# of the review's eq. 1, C_G/H being the liquid-side equilibrium concentration.
# The gas enters at z = 0, where the sparger is, and leaves at z = H_t.
# The problem is linear: the steady column is one sparse solve, and the
# transient is one LU factorisation reused at every step.
# ---------------------------------------------------------------------

def _operators(n, H_t, v_s, nu=0, inlet_top=False, dirichlet_both=False, C_in=1.0):
    """Upwind convection operator of the gas column, assembled once.

    First-order upwind face values on ``n`` equal cells of ``[0, H_t]``.

    Parameters
    ----------
    n : int
        Number of cells.
    H_t : float
        Column height, m.
    v_s : float
        Superficial gas velocity, m/s.
    nu : int, optional
        Geometry code of ``construct_div``; 0 is a constant cross-section.
    inlet_top : bool, optional
        Defect-injection switch: feed the gas at the top, flowing down.
    dirichlet_both : bool, optional
        Defect-injection switch: impose the inlet value at the outlet too.
    C_in : float, optional
        Inlet gas concentration.

    Returns
    -------
    dict
        ``L`` (``Div @ Conv``, sparse), ``g`` (its boundary vector, flat),
        ``bc``, ``x_f`` and the signed velocity ``v``, so that the face
        values can be read back with ``compute_boundary_values``.
    """
    shape, x_f = (n, 1), np.linspace(0.0, H_t, n + 1)
    # inlet, a dC/dn + b C = d on the outward normal: C_G = C_in (a = 0, b = 1)
    inlet = {"a": 0.0, "b": 1.0, "d": C_in}
    # outlet: pure convection leaves; the face carries the last cell's value
    outlet = {"outflow": True}
    if inlet_top:                       # BREAK ROW: gas fed at z = H_t, flowing down
        bc, v = (outlet, inlet), -v_s
    elif dirichlet_both:                # BREAK ROW: C_G = C_in imposed at the outlet face too
        bc, v = (inlet, dict(inlet)), v_s
    else:
        bc, v = (inlet, outlet), v_s
    Conv, conv_bc = construct_convflux_upwind(shape, x_f, bc=bc, v=v)
    Div = construct_div(shape, x_f, nu=nu)
    return dict(L=(Div @ Conv).tocsc(), g=np.asarray((Div @ conv_bc).todense()).ravel(),
                bc=bc, x_f=x_f, v=v)

def _net_inflow(C, op, read="face"):
    """Convective inflow minus outflow per unit cross-section.

    ``read="face"`` takes the face values the operator transported, from
    ``compute_boundary_values`` under the same bc; ``read="cell"`` takes the
    two end cells instead (a defect-injection switch).
    """
    if read == "face":
        c_lo, _, c_hi, _ = compute_boundary_values(C.reshape(-1, 1), op["x_f"], bc=op["bc"])
        c_lo, c_hi = np.ravel(c_lo)[0], np.ravel(c_hi)[0]
    else:                               # BREAK ROW: end cells read as face values
        c_lo, c_hi = C[0], C[-1]
    return float(op["v"] * (c_lo - c_hi))

def stripping_number(k, H_t, v_s, eps, Henry):
    """N_G: how many e-foldings of gas-side driving force the column strips.

    A gas parcel spends tau_G = eps*H_t/v_s in the column (the review's eq. 3)
    and loses oxygen at (1-eps) k (C_G/H - C_L) per unit dispersion volume, i.e.
    at (1-eps) k / eps per unit gas volume. A driving force u in liquid units is
    H*u in gas units, so

        N_G = (1-eps) k tau_G / (eps H) = (1-eps) k H_t / (v_s H)

    and the gas holdup cancels except through (1-eps). The review's own
    criterion for this problem, tau_G against 1/(k_L A), omits the factor
    (1-eps)/(eps H).

    Parameters
    ----------
    k : float
        k_L*A on a liquid-volume basis, 1/s.
    H_t : float
        Dispersion height, m.
    v_s : float
        Superficial gas velocity, m/s.
    eps : float
        Gas holdup.
    Henry : float
        Henry constant H of the review's eq. 1.

    Returns
    -------
    float
        ``N_G``, dimensionless.
    """
    return (1.0 - eps) * k * H_t / (v_s * Henry)

def phi_closed(N_G):
    """Mean driving force over inlet driving force for a quasi-steady gas column.

    Closed form, derived independently of any discretisation: with C_L frozen
    the gas driving force decays exponentially along z, so its column mean is
    (1 - exp(-N_G))/N_G. The apparent k_L*A that eq. 2 returns is phi times the
    true one.

    Parameters
    ----------
    N_G : array_like
        Stripping number from :func:`stripping_number`.

    Returns
    -------
    numpy.ndarray
        ``phi``, between 0 and 1.
    """
    N_G = np.asarray(N_G, float)
    return np.where(N_G < 1e-8, 1.0 - N_G / 2.0, (1.0 - np.exp(-N_G)) / np.where(N_G == 0, 1.0, N_G))

def phi_pymrm(n, H_t, v_s, eps, k, Henry, C_in=8.5, C_L_frac=0.3,
              nu=0, phase=None, sign=+1.0, return_balance=False, read="face", **bckw):
    """The same quantity from the pymrm operators: quasi-steady gas at frozen C_L.

    Linear, so one sparse solve.

    Parameters
    ----------
    n : int
        Number of cells.
    H_t, v_s, eps, k, Henry : float
        As in :func:`stripping_number`.
    C_in : float, optional
        Inlet gas concentration.
    C_L_frac : float, optional
        Frozen liquid concentration as a fraction of ``C_in / Henry``.
    nu : int, optional
        Geometry code of ``construct_div``. Defect-injection switch.
    phase : float, optional
        Interfacial-area weight; default ``1 - eps``. Defect-injection switch.
    sign : float, optional
        Sign of the transfer term. Defect-injection switch.
    return_balance : bool, optional
        Also return the relative gap of the column balance: convective net
        inflow through the two faces against the transfer summed over the
        cells. Zero by construction for ``nu = 0``; its size is the solve's
        round-off.
    read : {"face", "cell"}, optional
        How the balance reads the inlet and outlet values (:func:`_net_inflow`).
    **bckw
        Boundary switches passed to :func:`_operators` (``inlet_top``,
        ``dirichlet_both``).

    Returns
    -------
    float or tuple of float
        ``phi``, the mean driving force over the inlet driving force, and the
        balance gap when asked for.
    """
    ph = (1.0 - eps) if phase is None else phase
    C_L = C_L_frac * C_in / Henry
    op = _operators(n, H_t, v_s, nu=nu, C_in=C_in, **bckw)
    A = op["L"] + sp.eye(n, format="csc") * (sign * ph * k / Henry)
    b = -op["g"] + sign * ph * k * C_L
    c = spsolve(A.tocsc(), b)
    phi = float((c.mean() / Henry - C_L) / (C_in / Henry - C_L))
    if not return_balance:
        return phi
    net = _net_inflow(c, op, read)
    sink = sign * ph * k * H_t * float(np.mean(c / Henry - C_L))
    return phi, abs(net - sink) / max(abs(net), 1e-300)

def gassing_out(n, nt, t_end, H_t, v_s, eps, k, Henry, C_in=8.5,
                gas_initially_inert=False, nu=0, phase=None, sign=+1.0, read="face", **bckw):
    """The dynamic gassing-out experiment: gas column and liquid, backward Euler.

    Both phases in one linear system per step (monolithic backward Euler, first
    order in the step), one LU factorisation for every step. The liquid starts
    oxygen-free; the gas starts at the inlet composition, or oxygen-free.

    Parameters
    ----------
    n : int
        Number of gas cells.
    nt : int
        Number of time steps; the trace is sampled at ``(i + 1) * t_end / nt``.
    t_end : float
        End time, s.
    H_t, v_s, eps, k, Henry : float
        As in :func:`stripping_number`.
    C_in : float, optional
        Inlet gas concentration.
    gas_initially_inert : bool, optional
        Start with an oxygen-free gas column instead of one at ``C_in``.
    nu, phase, sign, read, **bckw
        Defect-injection switches, as in :func:`phi_pymrm`.

    Returns
    -------
    ts : numpy.ndarray
        Times, ``nt + 1`` values.
    CLs : numpy.ndarray
        Liquid concentration at those times.
    Cs : numpy.ndarray
        Gas concentration profiles, shape ``(nt + 1, n)``.
    balance : float
        Relative gap of the two-phase oxygen balance: inventory change against
        the convective net inflow read at the faces with
        ``compute_boundary_values``. Backward Euler telescopes, so the gap is
        the solves' round-off.
    """
    ph = (1.0 - eps) if phase is None else phase
    op = _operators(n, H_t, v_s, nu=nu, C_in=C_in, **bckw)
    dt, w = t_end / nt, np.full(n, 1.0 / n)
    I = sp.eye(n, format="csc")
    A = sp.bmat([[eps / dt * I + op["L"] + (sign * ph * k / Henry) * I,
                  sp.csc_matrix(-sign * ph * k * np.ones((n, 1)))],
                 [sp.csc_matrix(-(k / Henry) * w.reshape(1, n)),
                  sp.csc_matrix([[1.0 / dt + k]])]], format="csc")
    lu = splu(A)
    C0 = np.zeros(n) if gas_initially_inert else np.full(n, C_in)
    C, C_L, ts, CLs, Cs = C0.copy(), 0.0, [0.0], [0.0], [C0.copy()]
    net_in = 0.0        # cumulative convective input per unit dispersion volume
    for m in range(nt):
        sol = lu.solve(np.concatenate([eps / dt * C - op["g"], [C_L / dt]]))
        C, C_L = sol[:n], sol[n]
        net_in += _net_inflow(C, op, read) / H_t * dt
        ts.append((m + 1) * dt); CLs.append(C_L); Cs.append(C.copy())
    balance = abs(eps * (C - C0).mean() + (1.0 - eps) * C_L - net_in) / max(abs(net_in), 1e-300)
    return np.array(ts), np.array(CLs), np.array(Cs), float(balance)

def crossing_time(ts, CLs, level):
    """Time at which the piecewise-linear interpolant of a trace crosses ``level``.

    A root of ``interp(t) - level`` on the whole trace with
    ``pymrm.checks.find_roots``; raises unless there is exactly one.

    Parameters
    ----------
    ts, CLs : numpy.ndarray
        Time and liquid-concentration trace.
    level : float
        Concentration to cross.

    Returns
    -------
    float
        The crossing time, s.
    """
    res = find_roots(lambda t: float(np.interp(t, ts, CLs)) - level, ts[0], ts[-1],
                     n_scan=400, xtol=1e-12 * ts[-1])
    if len(res["roots"]) != 1:
        raise RuntimeError(f"trace crosses {level:.4g} {len(res['roots'])} times: {res['message']}")
    return float(res["roots"][0])

def fit_eq2(ts, CLs, C_star, f1=0.2, f2=0.8):
    """What an experimenter gets by applying the review's eq. 2 to a C_L trace.

    Eq. 2 is  k_L*A = ln[(C_G/H - C_L(t1)) / (C_G/H - C_L(t2))] / (t2 - t1),
    with C_G taken as the inlet gas composition, which is the well-mixed-gas
    assumption this page tests. t1 and t2 are where C_L crosses f1 and f2 of
    saturation (:func:`crossing_time`), so C_L(t1) = f1 C_star exactly. The
    window (f1, f2) is chosen; the review prescribes none.

    Parameters
    ----------
    ts, CLs : numpy.ndarray
        Time and liquid-concentration trace.
    C_star : float
        Saturation concentration ``C_in / Henry``.
    f1, f2 : float, optional
        Fractions of saturation that place ``t1`` and ``t2``.

    Returns
    -------
    float
        The apparent k_L*A, 1/s.
    """
    t1, t2 = crossing_time(ts, CLs, f1 * C_star), crossing_time(ts, CLs, f2 * C_star)
    return float(np.log((1.0 - f1) / (1.0 - f2)) / (t2 - t1))

def lagrangian(M, nt, t_end, H_t, v_s, eps, k, Henry, C_in=8.5,
               gas_initially_inert=False):
    """Independent route: integrate along the gas characteristics.

    No spatial operator, no pymrm matrix, no finite-volume cell. Every gas
    parcel is labelled by its age a in [0, tau_G]; a sub-step h <= da first
    strips each parcel (explicit midpoint rule), then moves the age labels on
    by h, with linear interpolation between the age nodes and the inlet value
    for ages below h. The column mean is the trapezoidal average over ages.
    It shares the physical model with :func:`gassing_out` and nothing else.

    Parameters
    ----------
    M : int
        Number of age intervals.
    nt : int
        Number of output steps.
    t_end : float
        End time, s.
    H_t, v_s, eps, k, Henry : float
        As in :func:`stripping_number`.
    C_in : float, optional
        Inlet gas concentration.
    gas_initially_inert : bool, optional
        Start with an oxygen-free gas column.

    Returns
    -------
    ts, CLs : numpy.ndarray
        Time and liquid-concentration trace, ``nt + 1`` values each.
    """
    tau_G = eps * H_t / v_s
    a = np.linspace(0.0, tau_G, M + 1)
    parcels = np.zeros(M + 1) if gas_initially_inert else np.full(M + 1, C_in)
    lam, dt = (1.0 - eps) * k / (eps * Henry), t_end / nt
    da = a[1] - a[0]
    C_L, ts, CLs = 0.0, [0.0], [0.0]
    for m in range(nt):
        sub = max(1, int(np.ceil(dt / da)))
        h = dt / sub
        for _ in range(sub):
            k1 = -lam * (parcels - C_L * Henry)
            k2 = -lam * (parcels + 0.5 * h * k1 - C_L * Henry)
            parcels = parcels + h * k2
            parcels = np.interp(a - h, a, parcels, left=C_in)   # move the age labels on by h
            parcels[a - h < 0] = C_in
            mean = np.trapezoid(parcels, a) / tau_G
            l1 = k * (mean / Henry - C_L)
            C_L = C_L + h * k * (mean / Henry - (C_L + 0.5 * h * l1))
        ts.append((m + 1) * dt); CLs.append(C_L)
    return np.array(ts), np.array(CLs)
