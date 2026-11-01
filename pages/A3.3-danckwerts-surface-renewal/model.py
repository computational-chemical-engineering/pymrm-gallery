"""Surface elements of Danckwerts (1951), marched in age and averaged over ages.

Model code of gallery page A3.3 (P. V. Danckwerts, Ind. Eng. Chem. 43, 1460,
1951). An element of liquid brought to the gas surface absorbs by unsteady
diffusion while it stays there,

    dc/dtheta = D d2c/dx2 - r c,      c(x, 0) = c_o,

and the mean rate per unit surface is the age average of its surface flux psi
over the distribution phi(theta) = s exp(-s theta) (the paper's eq. 7). All
quantities are in units where only ratios matter (``D = s = c* = 1`` on the
page unless stated).

Contents
--------
Element
    One element with a first-order reaction: the diffusion operator from
    ``construct_grad`` and ``construct_div`` assembled once, stored in banded
    form; backward Euler on age steps uniform in ``w = sqrt(theta)``.
TwoSpeciesElement
    Gas A absorbed at the surface reacting with dissolved reagent B at rate
    ``r' c c'`` (the paper's case 5): two fields, an analytic Jacobian of the
    pointwise source, one ``newton`` solve per age step.
age_average
    Eq. 7 integrated by parts onto the cumulative uptake, by Simpson's rule in w.
rate, rate_second_order
    The age-averaged rate by backward Euler at ``n_t`` and ``2 n_t`` steps, and
    its Richardson extrapolation in the step.
lam_eq16, eq15
    The instantaneous-reaction closed form (eqs. 15 and 16).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from scipy.sparse import eye_array
from scipy.special import erf, erfc, erfcx

from pymrm import (compute_boundary_values, construct_coefficient_matrix,
                   construct_div, construct_grad, newton)

__all__ = ["Element", "TwoSpeciesElement", "age_average", "rate",
           "rate_second_order", "lam_eq16", "eq15", "stretched_faces"]


def stretched_faces(depth, n_x, stretch):
    """Face positions on ``[0, depth]``, exponentially fine at the surface.

    The map of the cell index is fixed, so doubling ``n_x`` halves every cell.

    Parameters
    ----------
    depth : float
        Depth of the element.
    n_x : int
        Number of cells.
    stretch : float
        Stretch exponent; 0 gives uniform cells.

    Returns
    -------
    numpy.ndarray
        ``n_x + 1`` face positions.
    """
    u = np.linspace(0.0, 1.0, n_x + 1)
    return depth * np.expm1(stretch * u) / np.expm1(stretch) if stretch else depth * u


def age_grid(th_max, n_t):
    """Ages ``theta_j = (j dw)^2``, uniform in ``w = sqrt(theta)``."""
    return np.linspace(0.0, np.sqrt(th_max), n_t + 1) ** 2


class Element:
    """One surface element with an optional first-order reaction.

    The state has the field shape ``(n_x, 1)``: space first, one field. The
    operator ``D d2/dx2`` and its boundary vector are assembled once per set of
    boundary conditions from ``construct_grad`` and ``construct_div``, and the
    tridiagonal matrix is stored in ``solve_banded`` layout; a step adds
    ``1/dtheta + r`` to its diagonal and makes one banded solve.

    Parameters
    ----------
    D : float, optional
        Diffusivity.
    r : float, optional
        First-order rate constant.
    depth : float, optional
        Depth of the element in units of ``sqrt(D/s)``.
    n_x : int, optional
        Number of cells.
    stretch : float, optional
        Grid stretch exponent (see :func:`stretched_faces`).
    base : {"deep", "mixed"}, optional
        Condition at ``x = depth``: no flux (a truncated semi-infinite liquid)
        or ``c = c_o`` (a layer over a well-mixed bulk).
    """

    def __init__(self, D=1.0, r=0.0, depth=45.0, n_x=400, stretch=10.0, base="deep"):
        if base not in ("deep", "mixed"):
            raise ValueError("base must be 'deep' or 'mixed'")
        self.D, self.r, self.depth, self.n_x, self.base = float(D), float(r), depth, int(n_x), base
        self.x_f = stretched_faces(depth, self.n_x, stretch)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.V = np.diff(self.x_f)                       # cell volumes per unit area
        self.shape = (self.n_x, 1)                       # (space, field), never (n_x,)
        self.div = construct_div(self.shape, self.x_f, nu=0)   # nu = 0: Cartesian slab

    def bc(self, cstar=1.0, c0=0.0, k=None, sign=1.0):
        """Boundary conditions, outward normal: ``a dc/dn + b c = d``.

        Parameters
        ----------
        cstar, c0 : float
            Saturation and bulk concentrations.
        k : float or None
            Surface coefficient for eq. 29's condition; None for eq. 27's.
        sign : float
            Multiplies ``b`` of the Robin condition (1 is correct).

        Returns
        -------
        tuple of dict
        """
        # x = 0 (n = -x, dc/dn = -dc/dx):
        #   eq. 27, surface saturated:          c = c*
        #   eq. 29, surface resistance:         -D dc/dx = k (c* - c)
        #                                       <=> D dc/dn + k c = k c*
        left = ({"a": 0.0, "b": 1.0, "d": cstar} if k is None
                else {"a": self.D, "b": sign * k, "d": k * cstar})
        # x = depth (n = +x):
        #   deep liquid, no flux:               dc/dx = 0
        #   well-mixed bulk under a layer:      c = c_o
        right = ({"a": 1.0, "b": 0.0, "d": 0.0} if self.base == "deep"
                 else {"a": 0.0, "b": 1.0, "d": c0})
        return (left, right)

    def operators(self, bc):
        """``D d2/dx2`` for ``bc``: (matrix, boundary vector, banded matrix)."""
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, bc)
        jac = (self.D * (self.div @ grad)).tocsr()
        g = self.D * (self.div @ grad_bc).toarray().ravel()
        ab = np.zeros((3, self.n_x))
        ab[0, 1:], ab[1], ab[2, :-1] = jac.diagonal(1), jac.diagonal(0), jac.diagonal(-1)
        return jac, g, ab

    def face_fluxes(self, c, bc):
        """Surface flux into the element and flux out through the base.

        Parameters
        ----------
        c : numpy.ndarray
            One profile ``(n_x,)`` or a stack ``(n_x, m)``.
        bc : tuple of dict

        Returns
        -------
        tuple of numpy.ndarray
            ``(-D dc/dx at x = 0, -D dc/dx at x = depth)``, each of length m.
        """
        c = np.asarray(c, dtype=float).reshape(self.n_x, -1)
        _, g0, _, g1 = compute_boundary_values(c, self.x_f, self.x_c, bc)
        return -self.D * np.ravel(g0), -self.D * np.ravel(g1)   # gradients along +x

    def march(self, cstar=1.0, c0=0.0, k=None, th_max=30.0, n_t=1600, bc_sign=1.0,
              drop_reaction_in_flux=False, base_flux="face", keep_profiles=False,
              profiles_at=()):
        """March from ``theta = 0`` to ``th_max`` by backward Euler.

        Parameters
        ----------
        cstar, c0, k, bc_sign : float
            See :meth:`bc`.
        th_max : float, optional
            Last age.
        n_t : int, optional
            Number of steps, uniform in ``sqrt(theta)``.
        drop_reaction_in_flux : bool, optional
            Leave the reaction out of the mass balance that gives ``psi``.
        base_flux : {"face", "cell"}, optional
            Base outflow from ``compute_boundary_values`` (the flux the operator
            transports), or ``D (c_N - c_o)/(x_f - x_N)`` from the last cell
            centre (inconsistent with the operator, for the break table).
        keep_profiles : bool, optional
            Return every profile, shape ``(n_x, n_t)``.
        profiles_at : sequence of float, optional
            Ages at which to keep a profile (the first step at or past each).

        Returns
        -------
        dict
            ``theta``; ``psi``, the surface flux of every step from the discrete
            mass balance (accumulation + reaction + base outflow);
            ``Q``, its cumulative integral; ``psi_face``, the same flux read
            with ``compute_boundary_values``; ``prof``; ``profiles``.
        """
        bc = self.bc(cstar, c0, k, bc_sign)
        jac, g, ab_diff = self.operators(bc)
        th = age_grid(th_max, n_t)
        c = np.full(self.n_x, float(c0))
        psi = np.zeros(n_t + 1)
        Q = np.zeros(n_t + 1)
        P = np.empty((self.n_x, n_t))
        prof = {}
        r_in_flux = 0.0 if drop_reaction_in_flux else self.r
        for j in range(1, n_t + 1):
            dt = th[j] - th[j - 1]
            ab = -ab_diff
            ab[1] += 1.0 / dt + self.r
            c_old, c = c, solve_banded((1, 1), ab, c / dt + g, check_finite=False)
            P[:, j - 1] = c
            psi[j] = (self.V @ (c - c_old)) / dt + r_in_flux * (self.V @ c)
            Q[j] = Q[j - 1] + psi[j] * dt
            for t in profiles_at:
                if th[j - 1] < t <= th[j]:
                    prof[t] = (th[j], c.copy())
        f_surf, f_base = self.face_fluxes(P, bc)
        if base_flux == "cell" and self.base == "mixed":
            f_base = self.D * (P[-1] - c0) / (self.x_f[-1] - self.x_c[-1])
        psi[1:] += f_base                      # base outflow (zero for a deep base)
        Q[1:] += np.cumsum(f_base * np.diff(th))
        psi_face = np.concatenate([[np.nan], f_surf])
        out = {"theta": th, "psi": psi, "Q": Q, "psi_face": psi_face, "prof": prof}
        if keep_profiles:
            out["profiles"] = P
        return out


def age_average(th, Q, s):
    """Eq. 7 integrated by parts onto the cumulative uptake ``Q``.

    ``s int_0^T e^{-s th} psi d th = s e^{-sT} Q(T) + s^2 int_0^T e^{-s th} Q d th``,
    with ``th = w^2`` and Simpson's rule on the uniform ``w`` grid. ``Q(0) = 0``
    and the integrand ``2 w e^{-s w^2} Q`` is smooth at the origin.

    Parameters
    ----------
    th : numpy.ndarray
        Ages, uniform in ``sqrt(th)``, an even number of intervals.
    Q : numpy.ndarray
        Cumulative uptake at those ages.
    s : float
        Renewal rate.

    Returns
    -------
    float
    """
    if (len(th) - 1) % 2:
        raise ValueError("Simpson's rule needs an even number of age intervals")
    w = np.sqrt(th)
    gw = 2 * w * np.exp(-s * th) * Q
    h = w[1] - w[0]
    simpson = h / 3 * (gw[0] + 4 * gw[1:-1:2].sum() + 2 * gw[2:-2:2].sum() + gw[-1])
    return s ** 2 * simpson + s * np.exp(-s * th[-1]) * Q[-1]


def rate(s=1.0, n_t=1600, *, D=1.0, r=0.0, depth=45.0, n_x=400, stretch=10.0,
         base="deep", **march_kw):
    """Age-averaged rate from :class:`Element` at ``n_t`` and ``2 n_t`` steps.

    Backward Euler is first order in the step, so ``2 R(2 n_t) - R(n_t)``
    removes the leading error.

    Parameters
    ----------
    s : float
        Renewal rate.
    n_t : int
        The coarser step count.
    D, r, depth, n_x, stretch, base
        :class:`Element` arguments.
    **march_kw
        :meth:`Element.march` arguments.

    Returns
    -------
    dict
        ``raw`` (``n_t``), ``fine`` (``2 n_t``) and ``extrapolated``.
    """
    el = Element(D=D, r=r, depth=depth, n_x=n_x, stretch=stretch, base=base)

    def one(nt):
        out = el.march(n_t=nt, **march_kw)
        return age_average(out["theta"], out["Q"], s)

    raw, fine = one(n_t), one(2 * n_t)
    return {"raw": raw, "fine": fine, "extrapolated": 2.0 * fine - raw}


class TwoSpeciesElement:
    """Danckwerts' case 5: dissolved A reacts with reagent B at rate ``r' c c'``.

    State ``(n_x, 2)``: space first, fields (A, B) last. The diffusion operator
    with a diffusivity per field is assembled once; each step solves

        (c - c_old)/dtheta = div(D grad c) + S(c),   S = (-q, -q),  q = r' c_A c_B

    with ``newton`` and the analytic Jacobian of ``S``, a 2 x 2 block per cell
    assembled with ``construct_coefficient_matrix``. With ``deplete=False`` the
    reagent is held at ``c_o'`` (``q = r' c_o' c_A``, B untouched), which is
    exactly the pseudo-first-order problem with ``r = r' c_o'``.

    Parameters
    ----------
    D, Dp : float, optional
        Diffusivities of A and B.
    rp : float, optional
        Second-order rate constant ``r'``.
    depth : float, optional
        Depth of the element.
    n_x : int, optional
        Number of cells.
    stretch : float, optional
        Grid stretch exponent.
    """

    def __init__(self, D=1.0, Dp=1.0, rp=1.0, depth=45.0, n_x=250, stretch=10.0):
        self.D, self.Dp, self.rp, self.n_x = float(D), float(Dp), float(rp), int(n_x)
        self.x_f = stretched_faces(depth, self.n_x, stretch)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.V = np.diff(self.x_f)
        self.shape = (self.n_x, 2)                       # (space, field): fields last

    def operators(self, cstar):
        """Diffusion matrix and boundary vector for a surface value ``cstar``."""
        # x = 0 (outward normal -x): A saturated, c_A = c*;
        #                            B non-volatile, -D' dc_B/dx = 0 (dc_B/dn = 0).
        # x = depth (outward normal +x): deep liquid, dc/dx = 0 for both.
        bc = ({"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[cstar, 0.0]]},
              {"a": [[1.0, 1.0]], "b": [[0.0, 0.0]], "d": [[0.0, 0.0]]})
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, bc)
        div = construct_div(self.shape, self.x_f, nu=0)    # nu = 0: Cartesian slab
        d_face = construct_coefficient_matrix(np.array([self.D, self.Dp]),
                                              shape=self.shape, axis=0)
        jac = (div @ d_face @ grad).tocsc()
        g = (div @ d_face @ grad_bc).toarray().ravel()
        return jac, g

    def source(self, c, c0p, deplete=True):
        """Pointwise source and its exact Jacobian.

        Parameters
        ----------
        c : numpy.ndarray
            State ``(n_x, 2)``.
        c0p : float
            Bulk reagent concentration (used when ``deplete`` is False).
        deplete : bool, optional
            Whether the reagent is consumed.

        Returns
        -------
        tuple
            ``S`` flattened ``(2 n_x,)`` and ``dS/dc`` as a sparse matrix.
        """
        cA = c[:, 0]
        cB = c[:, 1] if deplete else np.full(self.n_x, float(c0p))
        q = self.rp * cA * cB
        S = np.stack([-q, -q if deplete else np.zeros_like(q)], axis=1)
        dS = np.zeros((self.n_x, 2, 2))       # dS[i, a, b] = dS_a / dc_b in cell i
        dS[:, 0, 0] = -self.rp * cB
        if deplete:
            dS[:, 0, 1] = -self.rp * cA
            dS[:, 1, 0] = -self.rp * cB
            dS[:, 1, 1] = -self.rp * cA
        jac = construct_coefficient_matrix(dS, shape=((self.n_x, 2, 1), (self.n_x, 1, 2)))
        return S.ravel(), jac

    def march(self, cstar=1.0, c0p=50.0, th_max=30.0, n_t=400, deplete=True, tol=1e-10,
              frozen_jacobian=False):
        """March by backward Euler, one ``newton`` solve per step.

        Parameters
        ----------
        cstar, c0p : float, optional
            Surface value of A and bulk value of B.
        th_max : float, optional
            Last age.
        n_t : int, optional
            Steps, uniform in ``sqrt(theta)``.
        deplete : bool, optional
            Whether the reagent is consumed.
        tol : float, optional
            ``newton`` step tolerance on the unknowns scaled per field by
            ``(cstar, c0p)``.
        frozen_jacobian : bool, optional
            Give ``newton`` the source Jacobian of the frozen reagent while the
            reagent is consumed (a defect, for the break table).

        Returns
        -------
        dict
            ``theta``, ``Q`` (cumulative uptake of A from the discrete mass
            balance), ``c`` (last profile) and ``nit`` (Newton iterations per
            step). Raises ``RuntimeError`` if any step fails.
        """
        jac_diff, g_diff = self.operators(cstar)
        eye = eye_array(2 * self.n_x, format="csc")
        th = age_grid(th_max, n_t)
        scale = np.array([cstar, max(c0p, cstar)])       # per field, broadcast on (n_x, 2)
        c = np.zeros(self.shape)
        c[:, 1] = c0p
        Q, nits = np.zeros(n_t + 1), []
        for j in range(1, n_t + 1):
            dt = th[j] - th[j - 1]
            c_old = c.ravel().copy()

            def residual(cf):
                S, dS = self.source(cf.reshape(self.shape), c0p, deplete)
                if frozen_jacobian:
                    dS = self.source(cf.reshape(self.shape), c0p, False)[1]
                cf = cf.ravel()
                return ((cf - c_old) / dt - jac_diff @ cf - g_diff - S,
                        eye / dt - jac_diff - dS)

            sol = newton(residual, c, tol=tol, scale=scale, maxfev=60)
            if not sol.success:
                raise RuntimeError(f"newton failed at age {th[j]:.4g}: {sol.message}")
            new = sol.x.reshape(self.shape)
            nits.append(sol.nit)
            cB = new[:, 1] if deplete else c0p
            psi = (self.V @ (new[:, 0] - c[:, 0])) / dt + self.rp * (self.V @ (new[:, 0] * cB))
            Q[j] = Q[j - 1] + psi * dt
            c = new
        return {"theta": th, "Q": Q, "c": c, "nit": nits}


def rate_second_order(c0p, ratio, *, s=1.0, Dp=1.0, n_x=250, n_t=400, cstar=1.0,
                      deplete=True, th_max=30.0, frozen_jacobian=False):
    """Age-averaged rate of :class:`TwoSpeciesElement` at ``n_t`` and ``2 n_t``.

    Parameters
    ----------
    c0p : float
        Bulk reagent concentration over ``c*``.
    ratio : float
        ``r' c_o' / s``, the group the paper works in.
    s, Dp, n_x, n_t, cstar, deplete, th_max, frozen_jacobian
        As in :class:`TwoSpeciesElement` and its ``march``.

    Returns
    -------
    dict
        ``raw``, ``fine``, ``extrapolated`` and ``nit`` (all Newton counts).
    """
    el = TwoSpeciesElement(D=1.0, Dp=Dp, rp=ratio * s / c0p, depth=45.0, n_x=n_x,
                           stretch=10.0)
    outs = [el.march(cstar=cstar, c0p=c0p, th_max=th_max, n_t=nt, deplete=deplete,
                     frozen_jacobian=frozen_jacobian)
            for nt in (n_t, 2 * n_t)]
    raw, fine = (age_average(o["theta"], o["Q"], s) for o in outs)
    return {"raw": raw, "fine": fine, "extrapolated": 2.0 * fine - raw,
            "nit": outs[0]["nit"] + outs[1]["nit"]}


def lam_eq16(cstar, c0p, D, Dp, swap=False):
    """Solve eq. 16 for ``lambda = beta / sqrt(D)``.

    Written with ``erfcx(z) = e^{z^2} erfc(z)`` so that no exponential
    overflows.

    Parameters
    ----------
    cstar, c0p, D, Dp : float
    swap : bool, optional
        Read erf and erfc the wrong way round (a transcription defect).

    Returns
    -------
    float
    """
    def g(lam):
        mu = lam * np.sqrt(D / Dp)
        if swap:
            first = (cstar / np.sqrt(Dp)) * np.exp(mu * mu - lam * lam) * erf(mu)
            second = (c0p / np.sqrt(D)) * erfc(lam)
        else:
            first = (cstar / np.sqrt(Dp)) * np.exp(-lam * lam) * erfcx(mu)
            second = (c0p / np.sqrt(D)) * erf(lam)
        return first - second
    return brentq(g, 1e-12, 12.0, xtol=1e-15)


def eq15(cstar, c0p, D, Dp, s, swap=False):
    """Eq. 15, ``R = c* sqrt(Ds) / erf(beta/sqrt(D))``, with beta from eq. 16.

    Returns
    -------
    tuple
        ``(R, lambda)``.
    """
    lam = lam_eq16(cstar, c0p, D, Dp, swap)
    return cstar * np.sqrt(D * s) / erf(lam), lam
