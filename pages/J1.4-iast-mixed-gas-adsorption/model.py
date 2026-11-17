"""Myers & Prausnitz (1965): the ideal adsorbed solution theory (page J1.4).

Myers and Prausnitz, AIChE J. 11, 121 (1965). Everything here is algebra plus
one scalar root-find per state:

- pure-component isotherms as dictionaries of the amount ``n(P)``, the reduced
  spreading pressure ``z(P) = pi A/RT`` of eq. (19), its inverse ``Pinv(z)``
  and the Henry constant ``K`` (:func:`langmuir_pure`, :func:`henry_pure`,
  :func:`dualsite_pure`, whose inverse is a vectorised ``newton`` solve);
- the IAST solve of eqs. (36) to (42) for a vector of binary states, one
  scalar equation per state, with ``newton`` on ``NumJac((K, 1))``
  (:func:`iast`), and the extended-Langmuir rule it is compared with
  (:func:`ext_langmuir`, :func:`eqm_dev`, :func:`nt_dev`,
  :func:`nt_dev_sup`, :func:`n1_dilute_dev`, :func:`dilute_selectivities`);
- the spreading-pressure quadrature of eq. (19) for an isotherm known only
  pointwise (:func:`z_quad`) and the Gibbs adsorption isotherm (21) integrated
  along a composition path (:func:`path_rhs_gauss`, :func:`path_rhs_trap`);
- the loop integral of the Gibbs one-form for any mixture model
  (:func:`circulation`), the lever rules (45), (46) of Fig. 1
  (:func:`lever_rule`) and the relative spread of a selectivity
  (:func:`rel_spread`).

Every defect the page's break table injects enters through an argument of
these functions.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar
from pymrm import NumJac, newton

TOL = 1e-13          # newton step tolerance on ln z (or ln P), order-one unknowns


# --------------------------------------------------------- pure isotherms
def langmuir_pure(m, b, t=1.0):
    """Single-site Langmuir isotherm, or its Sips form for ``t != 1``.

    Parameters
    ----------
    m : float
        Capacity.
    b : float
        Affinity.
    t : float
        Heterogeneity exponent of ``n = m (bP)^t / (1 + (bP)^t)``; 1 is
        Langmuir. Below 1 the isotherm has no Henry region.

    Returns
    -------
    dict
        ``n``, ``z`` (eq. 19 in closed form, ``(m/t) ln(1 + (bP)^t)``),
        ``Pinv`` (its inverse) and ``K`` (Henry constant, ``m b``; inf for
        ``t < 1``).
    """
    if t == 1.0:
        return dict(n=lambda Pp: m * b * Pp / (1.0 + b * Pp),
                    z=lambda Pp: m * np.log1p(b * Pp),
                    Pinv=lambda zz: np.expm1(zz / m) / b,
                    K=m * b)
    return dict(n=lambda Pp: m * (b * Pp)**t / (1.0 + (b * Pp)**t),
                z=lambda Pp: m / t * np.log1p((b * Pp)**t),
                Pinv=lambda zz: np.expm1(t * zz / m)**(1.0 / t) / b,
                K=np.inf)


def henry_pure(K):
    """Pure Henry isotherm ``n = K P``, the eq. (30) special case.

    Parameters
    ----------
    K : float
        Henry constant.

    Returns
    -------
    dict
        ``n``, ``z``, ``Pinv`` and ``K`` as in :func:`langmuir_pure`.
    """
    return dict(n=lambda Pp: K * Pp, z=lambda Pp: K * Pp, Pinv=lambda zz: zz / K, K=K)


def dualsite_pure(ma, ba, mb, bb, scale=1.0, tol=TOL):
    """Two-site Langmuir isotherm evaluated at ``scale * P``.

    ``z`` is analytic; the inverse ``Pinv`` is numerical, ``newton`` on
    ``ln P`` with ``NumJac((n, 1))``, the machinery an isotherm without an
    analytic inverse needs.

    Parameters
    ----------
    ma, ba, mb, bb : float
        Capacities and affinities of the two sites.
    scale : float
        Pressure scale: the same curve shifted in ``log P``.
    tol : float
        Newton step tolerance on ``ln P``.

    Returns
    -------
    dict
        ``n``, ``z``, ``Pinv`` and ``K`` as in :func:`langmuir_pure`.

    Raises
    ------
    RuntimeError
        When the inverse does not converge.
    """
    def zf(Pp):
        return ma * np.log1p(ba * scale * Pp) + mb * np.log1p(bb * scale * Pp)

    def nf(Pp):
        return (ma * ba * scale * Pp / (1 + ba * scale * Pp)
                + mb * bb * scale * Pp / (1 + bb * scale * Pp))

    def Pinv(zz):
        zz = np.atleast_1d(np.asarray(zz, float))
        jac = NumJac((zz.size, 1))

        def res(wc):
            return (zf(np.exp(wc[:, 0])) - zz)[:, None]

        w0 = np.log(np.maximum(np.expm1(zz / (ma + mb)) / (ba * scale), 1e-300))
        sol = newton(lambda wc: jac(res, wc), w0.reshape(-1, 1), tol=tol, maxfev=200,
                     stagnation=True)
        if not sol.success:
            raise RuntimeError(f"dual-site isotherm inverse: newton failed ({sol.message})")
        return np.exp(sol.x[:, 0])

    return dict(n=nf, z=zf, Pinv=Pinv, K=(ma * ba + mb * bb) * scale)


# ------------------------------------------------------------- IAST solve
def iast(Pt, y1v, pure1, pure2, tol=TOL, backward_tol=1e-8, zfac2=1.0, nt_mode="harmonic"):
    """Solve eqs. (36) to (42) for a vector of binary states.

    One unknown per state, ``w = ln z`` with ``z = pi A/RT``, and the residual
    ``R = P y1/P1(z) + P y2/P2(z) - 1``: ``K`` independent scalar equations,
    solved together by ``newton`` on ``NumJac((K, 1))``.

    Parameters
    ----------
    Pt : float or ndarray
        Total pressure, broadcast against ``y1v``.
    y1v : ndarray
        Gas mole fraction of component 1.
    pure1, pure2 : dict
        Pure isotherms (:func:`langmuir_pure` and siblings).
    tol, backward_tol : float or None
        Newton step tolerance and backward-error tolerance.
    zfac2 : float
        Defect: component 2 read at ``zfac2 * z``, a broken eq. (40). 1 is IAST.
    nt_mode : {"harmonic", "arith"}
        Eq. (26) as printed (harmonic) or read as an arithmetic mean (defect).

    Returns
    -------
    dict
        ``z``, ``P1``, ``P2``, ``x1``, ``nt``, ``n1``, ``n2`` and ``s12``
        (selectivity, eq. 28; nan where y1 is 0 or 1).

    Raises
    ------
    RuntimeError
        When newton does not converge.
    """
    Pt = np.broadcast_to(np.asarray(Pt, float), np.shape(y1v)).ravel().copy()
    y1v = np.asarray(y1v, float).ravel()
    jac = NumJac((y1v.size, 1))

    def res(wc):
        zz = np.exp(wc[:, 0])
        return (Pt * y1v / pure1["Pinv"](zz) + Pt * (1 - y1v) / pure2["Pinv"](zfac2 * zz)
                - 1.0)[:, None]

    z0 = np.maximum(pure1["z"](Pt), pure2["z"](Pt))     # upper start: R(z0) <= 0, monotone
    sol = newton(lambda wc: jac(res, wc), np.log(z0).reshape(-1, 1), tol=tol, maxfev=200,
                 backward_tol=backward_tol, stagnation=True)
    if not sol.success:
        raise RuntimeError(f"IAST newton failed ({sol.message})")
    zz = np.exp(sol.x[:, 0])
    P1v, P2v = pure1["Pinv"](zz), pure2["Pinv"](zfac2 * zz)
    x1v = Pt * y1v / P1v                                          # eq. (22)
    n1o, n2o = pure1["n"](P1v), pure2["n"](P2v)
    if nt_mode == "harmonic":
        ntv = 1.0 / (x1v / n1o + (1 - x1v) / n2o)                 # eq. (26)
    elif nt_mode == "arith":
        ntv = x1v * n1o + (1 - x1v) * n2o
    else:
        raise ValueError(f"nt_mode must be 'harmonic' or 'arith', not {nt_mode!r}")
    s12 = (x1v / np.where(y1v > 0, y1v, np.nan)) / ((1 - x1v) / np.where(y1v < 1, 1 - y1v, np.nan))
    return dict(z=zz, P1=P1v, P2=P2v, x1=x1v, nt=ntv, n1=x1v * ntv, n2=(1 - x1v) * ntv, s12=s12)


def ext_langmuir(Pt, y1v, m1, b1, m2, b2):
    """Extended Langmuir, the mixture rule the paper reports as inconsistent.

    Parameters
    ----------
    Pt, y1v : float or ndarray
        Total pressure and gas mole fraction of component 1.
    m1, b1, m2, b2 : float
        Capacities and affinities.

    Returns
    -------
    tuple of ndarray
        ``n1, n2``.
    """
    Dd = 1.0 + b1 * Pt * y1v + b2 * Pt * (1 - y1v)
    return m1 * b1 * Pt * y1v / Dd, m2 * b2 * Pt * (1 - y1v) / Dd


def eqm_dev(pure1, pure2, lang, Pv, yv, **iast_kw):
    """Largest deviation of IAST from extended Langmuir over a set of states.

    Parameters
    ----------
    pure1, pure2 : dict
        The isotherms IAST is solved with.
    lang : tuple
        ``(m1, b1, m2, b2)`` of the extended-Langmuir rule.
    Pv, yv : ndarray
        The states.
    **iast_kw
        Passed to :func:`iast`.

    Returns
    -------
    float
        ``max |n_i(IAST) - n_i(EL)| / n_t(EL)`` over states and components.
    """
    s = iast(Pv, yv, pure1, pure2, **iast_kw)
    n1e, n2e = ext_langmuir(Pv, yv, *lang)
    return float(np.max(np.abs(np.array([s["n1"] - n1e, s["n2"] - n2e]))
                        / np.maximum(n1e + n2e, 1e-300)))


def nt_dev(Pt, yv, pure1, pure2, lang):
    """Relative deviation of the extended-Langmuir total amount from IAST.

    Parameters
    ----------
    Pt : float
        Total pressure.
    yv : float or ndarray
        Gas mole fraction of component 1.
    pure1, pure2 : dict
        Pure isotherms.
    lang : tuple
        ``(m1, b1, m2, b2)`` of extended Langmuir.

    Returns
    -------
    ndarray
        ``|n_t(EL) - n_t(IAST)| / n_t(IAST)`` per state.
    """
    yv = np.atleast_1d(np.asarray(yv, float))
    s = iast(Pt, yv, pure1, pure2)
    n1e, n2e = ext_langmuir(Pt, yv, *lang)
    return np.abs((n1e + n2e) - s["nt"]) / s["nt"]


def nt_dev_sup(Pt, pure1, pure2, lang, n_scan, y_edge, refine=True, xatol=1e-13):
    """Maximum over y1 of :func:`nt_dev` at one pressure: scan, then refine.

    Parameters
    ----------
    Pt : float
        Total pressure.
    pure1, pure2 : dict
        Pure isotherms.
    lang : tuple
        ``(m1, b1, m2, b2)`` of extended Langmuir.
    n_scan : int
        Nodes of the uniform scan of ``[y_edge, 1 - y_edge]``.
    y_edge : float
        Distance of the scan ends from 0 and 1.
    refine : bool
        Refine the scan maximum by a bounded maximisation between the scan
        nodes next to it; False returns the scan read-off (a defect).
    xatol : float
        Absolute tolerance on y1 of the refinement.

    Returns
    -------
    tuple of float
        The maximum and the y1 where it sits.

    Raises
    ------
    RuntimeError
        When the scan maximum sits on a scan end (the maximum is then at or
        beyond the end).
    """
    ys = np.linspace(y_edge, 1.0 - y_edge, n_scan)
    dev = nt_dev(Pt, ys, pure1, pure2, lang)
    i = int(np.argmax(dev))
    if i in (0, n_scan - 1):
        raise RuntimeError(f"maximum on the scan end y1 = {ys[i]:g}: at or beyond the end")
    if not refine:
        return float(dev[i]), float(ys[i])
    r = minimize_scalar(lambda y: -float(nt_dev(Pt, y, pure1, pure2, lang)[0]),
                        bounds=(ys[i - 1], ys[i + 1]), method="bounded",
                        options=dict(xatol=xatol))
    if not r.success:
        raise RuntimeError(f"bounded maximisation failed ({r.message})")
    return float(-r.fun), float(r.x)


def n1_dilute_dev(Pt, pure1, pure2, m1, b1, b2):
    """Closed-form y1 -> 0 limit of the component-1 deviation of extended Langmuir.

    Both models are linear in y1 there; the ratio of slopes is
    ``[m1 b1 P/(1 + b2 P)] / [(P/P1(z2(P))) n2(P)]``.

    Parameters
    ----------
    Pt : float
        Total pressure.
    pure1, pure2 : dict
        Pure isotherms of the IAST side.
    m1, b1, b2 : float
        Extended-Langmuir constants (the slope does not contain m2).

    Returns
    -------
    float
        ``|slope_EL / slope_IAST - 1|``.
    """
    slope_el = m1 * b1 * Pt / (1.0 + b2 * Pt)
    slope_iast = (Pt / pure1["Pinv"](pure2["z"](Pt))) * pure2["n"](Pt)
    return float(abs(slope_el / slope_iast - 1.0))


def dilute_selectivities(Pt, pure1, pure2):
    """Closed-form infinite-dilution selectivities of IAST at pressure ``Pt``.

    Parameters
    ----------
    Pt : float
        Total pressure.
    pure1, pure2 : dict
        Pure isotherms.

    Returns
    -------
    tuple of float
        ``s(y1 -> 0) = P / P1(z2(P))`` and ``s(y1 -> 1) = P2(z1(P)) / P``.
    """
    return (float(Pt / pure1["Pinv"](pure2["z"](Pt))),
            float(pure2["Pinv"](pure1["z"](Pt)) / Pt))


# ----------------------------------------------------- quadrature, paths
def _quad(f, x, rule):
    """Trapezoid (``"trapezoid"``) or left rectangle (``"rectangle"``, a defect) sum."""
    if rule == "trapezoid":
        return np.trapezoid(f, x)
    if rule == "rectangle":
        return float(np.sum(f[:-1] * np.diff(x)))
    raise ValueError(f"rule must be 'trapezoid' or 'rectangle', not {rule!r}")


def z_quad(Pp, pure, n_nodes, t0_frac=1e-8, anchor=True, grid="log", rule="trapezoid"):
    """Eq. (19) for an isotherm known only pointwise: ``z = int_0^P n d ln t``.

    The segment ``[0, t0]`` is the Henry segment, added in closed form as
    ``K t0``; the rest is a quadrature on ``n_nodes`` nodes.

    Parameters
    ----------
    Pp : float
        Pressure.
    pure : dict
        Pure isotherm (``n`` and ``K`` are used).
    n_nodes : int
        Quadrature nodes on ``[t0, Pp]``.
    t0_frac : float
        ``t0 / Pp``.
    anchor : bool
        Add the Henry segment ``K t0`` (False is a defect).
    grid : {"log", "linear"}
        Nodes uniform in ``ln t`` or in ``t`` (the latter a defect).
    rule : {"trapezoid", "rectangle"}
        Quadrature rule.

    Returns
    -------
    float
    """
    t0 = t0_frac * Pp
    if grid == "log":
        tg = np.geomspace(t0, Pp, n_nodes)
        q = _quad(pure["n"](tg), np.log(tg), rule)
    elif grid == "linear":
        tg = np.linspace(t0, Pp, n_nodes)
        q = _quad(pure["n"](tg) / tg, tg, rule)
    else:
        raise ValueError(f"grid must be 'log' or 'linear', not {grid!r}")
    return float(pure["K"] * t0 + q) if anchor else float(q)


def path_rhs_gauss(ng, Pt, ya, yb, pure1, pure2, **iast_kw):
    """Right side of eq. (21) along ``y1 = ya .. yb`` at constant ``Pt``, Gauss-Legendre.

    Parameters
    ----------
    ng : int
        Gauss points.
    Pt, ya, yb : float
        Pressure and path ends.
    pure1, pure2 : dict
        Pure isotherms.
    **iast_kw
        Passed to :func:`iast` (defects).

    Returns
    -------
    float
        ``int (n1/y1 - n2/y2) dy1``, to be compared with ``z(yb) - z(ya)``.
    """
    xg, wg = leggauss(ng)
    yv = 0.5 * (yb - ya) * xg + 0.5 * (ya + yb)
    s = iast(Pt, yv, pure1, pure2, **iast_kw)
    return 0.5 * (yb - ya) * float(np.sum(wg * (s["n1"] / yv - s["n2"] / (1 - yv))))


def path_rhs_trap(n_nodes, Pt, ya, yb, pure1, pure2, rule="trapezoid", **iast_kw):
    """As :func:`path_rhs_gauss`, on ``n_nodes`` uniform nodes.

    Parameters
    ----------
    n_nodes : int
        Nodes on ``[ya, yb]``.
    Pt, ya, yb : float
        Pressure and path ends.
    pure1, pure2 : dict
        Pure isotherms.
    rule : {"trapezoid", "rectangle"}
        Quadrature rule.
    **iast_kw
        Passed to :func:`iast` (defects).

    Returns
    -------
    float
    """
    yv = np.linspace(ya, yb, n_nodes)
    s = iast(Pt, yv, pure1, pure2, **iast_kw)
    return float(_quad(s["n1"] / yv - s["n2"] / (1 - yv), yv, rule))


# ------------------------------------------------- loop, lever, spreads
def circulation(model, ngauss=48, skip_edge4=False, Pa=0.2, Pb=1.0, ya=0.2, yb=0.8):
    """Loop integral of ``n_t d ln P + (n1/y1 - n2/y2) d y1`` around a rectangle.

    Zero for any mixture model with a spreading pressure (a state function).

    Parameters
    ----------
    model : callable
        ``model(P, y1) -> (n1, n2)`` on arrays.
    ngauss : int
        Gauss-Legendre points per edge.
    skip_edge4 : bool
        Drop the left edge (a defect).
    Pa, Pb, ya, yb : float
        The rectangle in ``(P, y1)``.

    Returns
    -------
    float
        The circulation, counter-clockwise in ``(ln P, y1)``.
    """
    xg, wg = leggauss(ngauss)
    hl = 0.5 * (np.log(Pb) - np.log(Pa))
    lnv = hl * xg + 0.5 * (np.log(Pa) + np.log(Pb))
    yv = 0.5 * (yb - ya) * xg + 0.5 * (ya + yb)
    tot = 0.0
    n1, n2 = model(np.exp(lnv), np.full_like(lnv, ya))
    tot += hl * np.sum(wg * (n1 + n2))                                  # bottom, ->
    n1, n2 = model(np.full_like(yv, Pb), yv)
    tot += 0.5 * (yb - ya) * np.sum(wg * (n1 / yv - n2 / (1 - yv)))     # right, ^
    n1, n2 = model(np.exp(lnv), np.full_like(lnv, yb))
    tot -= hl * np.sum(wg * (n1 + n2))                                  # top, <-
    if not skip_edge4:
        n1, n2 = model(np.full_like(yv, Pa), yv)
        tot -= 0.5 * (yb - ya) * np.sum(wg * (n1 / yv - n2 / (1 - yv)))  # left, v
    return float(tot)


def lever_rule(pure1, pure2, zstar, x1star, swap=False):
    """The construction of Fig. 1 and its lever rules (45), (46), against (22), (29).

    Parameters
    ----------
    pure1, pure2 : dict
        Pure isotherms.
    zstar : float
        Spreading pressure of the horizontal ABC.
    x1star : float
        Adsorbed composition that places D on BC (eq. 29).
    swap : bool
        Draw the chords OB and OC the wrong way round (a defect).

    Returns
    -------
    dict
        Abscissae ``P1``, ``P2``, ``Pd``, ordinates ``F``, ``E``, the lever
        readings ``y1``, ``x1``, the direct ``y1_direct`` and ``dev``, the
        larger deviation of the two lever readings from (22) and (29).
    """
    P1s, P2s = float(pure1["Pinv"](zstar)), float(pure2["Pinv"](zstar))
    Pd = P1s * x1star + P2s * (1 - x1star)               # eq. (29): D chosen via x1
    Fh, Eh = zstar * Pd / P1s, zstar * Pd / P2s          # chords OB, OC at abscissa Pd
    if swap:
        Fh, Eh = Eh, Fh
    y1_lever = (zstar - Eh) / (Fh - Eh)                  # eq. (45)
    x1_lever = (P2s - Pd) / (P2s - P1s)                  # eq. (46)
    y1_direct = P1s * x1star / Pd                        # eq. (22)
    return dict(P1=P1s, P2=P2s, Pd=Pd, F=Fh, E=Eh, y1=y1_lever, x1=x1_lever,
                y1_direct=y1_direct,
                dev=max(abs(y1_lever - y1_direct), abs(x1_lever - x1star)))


def rel_spread(s):
    """Relative spread ``(max - min) / mean`` of the finite entries of ``s``.

    Parameters
    ----------
    s : ndarray

    Returns
    -------
    float
    """
    s = np.asarray(s, float)
    s = s[np.isfinite(s)]
    return float(np.ptp(s) / np.mean(s))
