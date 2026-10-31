"""Whitman (1923) two-film theory of gas absorption.

Model code of gallery page A3.1. Two stagnant films in series, a gas film
carrying the partial pressure ``p`` (mmHg) and a liquid film carrying the
concentration ``c`` (g/l), with local equilibrium ``p2 = f(c2)`` at the
interface between them. Whitman's coefficients ``k_p`` and ``k_c`` are
per-apparatus conductances; with films of unit thickness they are the
diffusivities ``D/delta`` directly.

Contents
--------
Equilibrium
    The straight line ``log10 p = A + B c`` that Whitman's Fig. 3 plots, fitted
    to equilibrium pairs.
two_film_algebraic
    Whitman's Eq. 4 closed with ``p2 = f(c2)``: a scalar root, no grid.
no_liquid_film_rate
    The ``k_c -> infinity`` limit of Eq. 4 in closed form.
reacting_film_rate
    The same interface with a first-order reaction in the liquid film, from
    the exact film profile: a scalar root, no grid.
TwoFilm
    The two films as one pymrm finite-volume system with the interface values
    as unknowns, solved by ``pymrm.newton``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import eye_array

from pymrm import BlockLayout, compute_boundary_values, construct_div, construct_grad, newton
from pymrm.checks import residual_check

__all__ = ["Equilibrium", "two_film_algebraic", "no_liquid_film_rate",
           "reacting_film_rate", "TwoFilm"]


@dataclass(frozen=True)
class Equilibrium:
    """Equilibrium partial pressure over the liquid, ``log10(p/mmHg) = A + B c``.

    Parameters
    ----------
    A : float
        Intercept, ``log10`` of the pressure at ``c = 0``.
    B : float
        Slope, per g/l.
    """

    A: float
    B: float

    @classmethod
    def fit(cls, c, p):
        """Least-squares straight line through ``(c, log10 p)``.

        Parameters
        ----------
        c, p : array_like
            Equilibrium pairs, g/l and mmHg.

        Returns
        -------
        Equilibrium
        """
        B, A = np.polyfit(np.asarray(c, float), np.log10(np.asarray(p, float)), 1)
        return cls(float(A), float(B))

    def p(self, c):
        """Equilibrium partial pressure (mmHg) at concentration ``c`` (g/l)."""
        return 10.0 ** (self.A + self.B * np.asarray(c, float))

    def c(self, p):
        """Concentration (g/l) in equilibrium with partial pressure ``p`` (mmHg)."""
        return (np.log10(np.asarray(p, float)) - self.A) / self.B

    def dpdc(self, c):
        """Slope ``dp/dc`` (mmHg per g/l) at ``c``."""
        return self.p(c) * np.log(10.0) * self.B

    def residuals(self, c, p):
        """Worst misfit of pairs ``(c, p)``: relative in ``p`` (%) and in ``c`` (g/l).

        Returns
        -------
        tuple of float
            ``(worst |p/f(c) - 1| in %, worst |c - f^-1(p)| in g/l)``.
        """
        r = np.log10(np.asarray(p, float)) - (self.A + self.B * np.asarray(c, float))
        return float(np.abs(10.0 ** r - 1.0).max() * 100.0), float(np.abs(r / self.B).max())


def two_film_algebraic(p1, c3, k_p, k_c, peq):
    """Whitman's Eq. 4, ``k_p (p1 - p2) = k_c (c2 - c3)``, closed with ``p2 = f(c2)``.

    A scalar root in ``c2``; no grid and no operators. The rate is returned
    from the gas side, ``k_p (p1 - p2)``. Eq. 4 makes the liquid side
    ``k_c (c2 - c3)`` equal to it, but at large ``k_c`` the difference
    ``c2 - c3`` is a difference of nearly equal floats, so that side is
    returned separately (``N_liquid_side``) and used for nothing else.

    Parameters
    ----------
    p1 : float
        Bulk-gas partial pressure (mmHg).
    c3 : float
        Bulk-liquid concentration (g/l).
    k_p, k_c : float
        Film conductances, g/(h mmHg) and g/(h (g/l)).
    peq : callable
        Equilibrium pressure ``f(c)``, increasing.

    Returns
    -------
    dict
        ``N`` (g/h), ``p2``, ``c2``, ``N_liquid_side``.
    """
    def f(x):
        return k_c * (x - c3) - k_p * (p1 - peq(x))

    lo, hi, step = c3, c3 + 1e-12, max(1.0, 0.05 * abs(c3) + 1.0)
    while f(hi) < 0.0 and hi < c3 + 1e7:
        hi += step
        step *= 1.6
    c2 = brentq(f, lo, hi, xtol=1e-13, rtol=4 * np.finfo(float).eps)
    p2 = float(peq(c2))
    return dict(N=k_p * (p1 - p2), p2=p2, c2=c2, N_liquid_side=k_c * (c2 - c3))


def no_liquid_film_rate(p1, c3, k_p, peq):
    """The ``k_c -> infinity`` limit of Eq. 4, in closed form.

    With no liquid-film resistance the interface concentration is the bulk
    one, so ``p2 = f(c3)`` and the rate is ``k_p (p1 - f(c3))``.

    Returns
    -------
    dict
        ``N`` (g/h), ``p2``, ``c2``.
    """
    p2 = float(peq(c3))
    return dict(N=k_p * (p1 - p2), p2=p2, c2=float(c3))


def reacting_film_rate(p1, c3, k_p, k_c, ka, peq, ceq):
    """Two films with a first-order reaction in the liquid film, from the exact profile.

    In a liquid film of unit thickness with diffusivity ``k_c``, rate constant
    ``ka``, ``c(0) = c2`` and ``c(1) = c3``, the profile is
    ``(c2 sinh(Ha (1 - x)) + c3 sinh(Ha x)) / sinh(Ha)`` with
    ``Ha = sqrt(ka / k_c)``, so the flux entering it is
    ``k_c Ha (c2 coth Ha - c3 / sinh Ha)``. Equated to ``k_p (p1 - f(c2))``
    this is a scalar root in ``c2`` on ``[0, f^-1(p1)]``.

    Parameters
    ----------
    p1, c3, k_p, k_c : float
        As in :func:`two_film_algebraic`.
    ka : float
        First-order rate constant in the liquid film (1/h), positive.
    peq, ceq : callable
        Equilibrium ``f(c)`` and its inverse.

    Returns
    -------
    dict
        ``N`` (g/h), ``p2``, ``c2``.
    """
    ha = np.sqrt(ka / k_c)
    coth = 1.0 / np.tanh(ha)
    csch = 2.0 * np.exp(-ha) / (1.0 - np.exp(-2.0 * ha))   # 1/sinh, no overflow

    def f(x):
        return k_c * ha * (x * coth - c3 * csch) - k_p * (p1 - peq(x))

    c2 = brentq(f, 0.0, float(ceq(p1)), xtol=1e-15, rtol=4 * np.finfo(float).eps)
    p2 = float(peq(c2))
    return dict(N=k_p * (p1 - p2), p2=p2, c2=c2)


class TwoFilm:
    """Steady diffusion through a gas film and a liquid film in series.

    The state is one column, ``(2 n + 2, 1)``, laid out with
    ``pymrm.BlockLayout`` along the axis that runs through the films: the gas
    film ``"p"`` ``(n, 1)``, the interface pressure ``"p2"`` ``(1, 1)``, the
    interface concentration ``"c2"`` ``(1, 1)`` and the liquid film ``"c"``
    ``(n, 1)``. The two interface values enter the films as external boundary
    values through ``shapes_d``, so every operator is constant and is built
    once. The row of ``"p2"`` is the equilibrium ``p2 - f(c2) = 0``, the row of
    ``"c2"`` the flux continuity across the interface face. Both rows are
    nonlinear only through ``f``, whose derivative is placed analytically.

    ``pymrm.construct_interface_matrices`` eliminates an interface whose
    conditions are linear (a partition coefficient); the equilibrium here is
    not, which is why the interface values are unknowns.

    Parameters
    ----------
    k_p, k_c : float
        Film conductances, g/(h mmHg) and g/(h (g/l)). With unit film
        thickness these are ``D/delta`` directly.
    peq, dpeq : callable
        Equilibrium pressure ``f(c)`` and its derivative.
    n : int, optional
        Cells per film.
    ka : float, optional
        First-order rate constant in the liquid film (1/h). ``0`` is Whitman.
    nu : int, optional
        Geometry index of ``construct_div``: ``0`` is a plane film. Exposed
        only so that a break row can set it wrongly.
    jump : bool, optional
        ``True`` applies ``p2 = f(c2)``. ``False`` replaces it by ``p2 = c2``,
        no equilibrium jump: a break row only.
    exact_jacobian : bool, optional
        ``False`` leaves ``-f'(c2)`` out of the Jacobian: a break row only.
    liquid_flux_factor : float, optional
        Multiplies the liquid flux in the flux-continuity row; ``1`` is the
        model, anything else a break row only.

    Attributes
    ----------
    layout : pymrm.BlockLayout
        Where each film and interface value sits in the state.
    jac_const : scipy.sparse.csc_array
        The constant part of the Jacobian.
    g_p1, g_c3 : numpy.ndarray
        The residual per unit bulk pressure and per unit bulk concentration.
    """

    def __init__(self, k_p, k_c, peq, dpeq, n=40, ka=0.0, nu=0, jump=True,
                 exact_jacobian=True, liquid_flux_factor=1.0):
        self.k_p, self.k_c, self.peq, self.dpeq = float(k_p), float(k_c), peq, dpeq
        self.n, self.ka, self.nu, self.jump = int(n), float(ka), nu, jump
        self.exact_jacobian = exact_jacobian
        self.liquid_flux_factor = float(liquid_flux_factor)
        self.x_f = np.linspace(0.0, 1.0, self.n + 1)   # unit thickness: D/delta = k
        self.layout = BlockLayout({"p": (self.n, 1), "p2": (1, 1), "c2": (1, 1),
                                   "c": (self.n, 1)}, axis=0)
        self.shape = self.layout.shape
        self._build_operators()

    def _build_operators(self):
        n, k_p, k_c, place = self.n, self.k_p, self.k_c, self.layout.place
        shape = self.layout.shapes["p"]
        # Both ends of both films are Dirichlet. Outward normal, a dphi/dn + b phi = d
        # with a = 0, b = 1, so phi = d at that end:
        #   gas film,    x = 0: p = p1 (bulk gas)      x = 1: p = p2 (interface)
        #   liquid film, x = 0: c = c2 (interface)     x = 1: c = c3 (bulk liquid)
        # With shapes_d, d = 1 is a coefficient on the external boundary values,
        # (p1, p2) and (c2, c3), so the operators do not change with them.
        self.bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
        grad, grad_lo, grad_hi = construct_grad(shape, self.x_f, None, self.bc, axis=0,
                                                shapes_d=((1, 1), (1, 1)))
        div = construct_div(shape, self.x_f, nu=self.nu, axis=0)  # nu = 0: plane film
        # face fluxes, n + 1 faces each: gas -k_p dp/dx, liquid -k_c dc/dx
        flux_p, flux_p_bulk, flux_p_int = -k_p * grad, -k_p * grad_lo, -k_p * grad_hi
        flux_c, flux_c_int, flux_c_bulk = -k_c * grad, -k_c * grad_lo, -k_c * grad_hi
        # the interface is the last face of the gas film and the first of the liquid film
        faces = BlockLayout({"faces": (n + 1, 1)}, axis=0)
        faces.subdomain("gas_interface", "faces", np.s_[n:])
        faces.subdomain("liquid_interface", "faces", np.s_[:1])
        at_gas, at_liq = faces.pick("gas_interface"), faces.pick("liquid_interface")
        one = eye_array(1, format="csc")
        lf = self.liquid_flux_factor
        jac = (
            # gas film: div(flux_p) = 0
            place(div @ flux_p, "p", "p") + place(div @ flux_p_int, "p", "p2")
            # liquid film: div(flux_c) + ka c = 0
            + place(div @ flux_c + self.ka * eye_array(n, format="csc"), "c", "c")
            + place(div @ flux_c_int, "c", "c2")
            # equilibrium row: p2 - f(c2) = 0 (f enters in residual())
            + place(one, "p2", "p2")
            # flux continuity at the interface face: gas flux - liquid flux = 0
            + place(at_gas @ flux_p, "c2", "p") + place(at_gas @ flux_p_int, "c2", "p2")
            - lf * place(at_liq @ flux_c, "c2", "c") - lf * place(at_liq @ flux_c_int, "c2", "c2"))
        if not self.jump:
            jac = jac - place(one, "p2", "c2")             # p2 - c2 = 0
        self.jac_const = jac.tocsc()
        self.g_p1, self.g_c3 = np.zeros(self.shape), np.zeros(self.shape)
        self.layout.put(self.g_p1, "p", (div @ flux_p_bulk).toarray())
        self.layout.put(self.g_p1, "c2", (at_gas @ flux_p_bulk).toarray())
        self.layout.put(self.g_c3, "c", (div @ flux_c_bulk).toarray())
        self.layout.put(self.g_c3, "c2", -lf * (at_liq @ flux_c_bulk).toarray())
        self.g_p1, self.g_c3 = self.g_p1.reshape(-1, 1), self.g_c3.reshape(-1, 1)

    def residual(self, u, p1, c3):
        """Residual of the films, the equilibrium and the flux continuity, and its Jacobian."""
        u = np.ascontiguousarray(u, dtype=float).reshape(self.shape)
        g = self.jac_const @ u.reshape(-1, 1) + self.g_p1 * p1 + self.g_c3 * c3
        if not self.jump:
            return g, self.jac_const
        c2 = float(self.layout.view(u, "c2")[0, 0])
        self.layout.add(g, "p2", -self.peq(c2))
        if not self.exact_jacobian:
            return g, self.jac_const
        return g, (self.jac_const + self.layout.place(
            np.array([[-float(self.dpeq(c2))]]), "p2", "c2")).tocsc()

    def initial_guess(self, p1, c3):
        """The deleted-liquid-film state: ``c2 = c3``, ``p2 = f(c3)``, linear profiles."""
        u = np.zeros(self.shape)
        c2 = float(c3)
        p2 = float(self.peq(c2)) if self.jump else c2
        x_c = 0.5 * (self.x_f[1:] + self.x_f[:-1])
        self.layout.put(u, "p", p1 + (p2 - p1) * x_c)
        self.layout.put(u, "p2", p2)
        self.layout.put(u, "c2", c2)
        self.layout.put(u, "c", c2 + (c3 - c2) * x_c)
        return u

    def scale(self, p1, c3):
        """Typical magnitude of each unknown: ``p1`` for pressures, ``max(c3, 1)`` g/l else."""
        s = np.empty(self.shape)
        for name in ("p", "p2"):
            self.layout.view(s, name)[...] = abs(p1)
        for name in ("c", "c2"):
            self.layout.view(s, name)[...] = max(abs(c3), 1.0)
        return s

    def solve(self, p1, c3, tol=1e-12, maxfev=50, require_success=True):
        """Solve with ``pymrm.newton`` from :meth:`initial_guess`.

        ``newton``'s step test is absolute, so ``scale`` gives every unknown its
        physical magnitude and ``tol`` acts on the step relative to it. The
        result is judged by its componentwise backward error
        (``pymrm.checks.residual_check``) at the returned state.

        Returns
        -------
        dict
            ``N`` (g/h, the gas-film flux at the interface face), ``p2``, ``c2``,
            the profiles ``p`` and ``c``, ``nit``, ``success``,
            ``backward_error``, and the two interface face fluxes read with
            ``compute_boundary_values`` (``N_gas_face``, ``N_liquid_face``).

        Raises
        ------
        RuntimeError
            When ``require_success`` and Newton does not converge.
        """
        def fun(u):
            return self.residual(u, p1, c3)

        res = newton(fun, self.initial_guess(p1, c3), tol=tol, maxfev=maxfev,
                     scale=self.scale(p1, c3))
        if require_success and not res.success:
            raise RuntimeError(f"TwoFilm at p1 = {p1:g}, c3 = {c3:g}: {res.message}")
        u = np.asarray(res.x, dtype=float).reshape(self.shape)
        lay = self.layout
        p, c = lay.view(u, "p").copy(), lay.view(u, "c").copy()
        p2, c2 = float(lay.view(u, "p2")[0, 0]), float(lay.view(u, "c2")[0, 0])
        # interface fluxes from the face gradients (gradients along +x)
        _, _, _, dp_hi = compute_boundary_values(
            p, self.x_f, None, ({"a": 0.0, "b": 1.0, "d": p1}, {"a": 0.0, "b": 1.0, "d": p2}))
        _, dc_lo, _, _ = compute_boundary_values(
            c, self.x_f, None, ({"a": 0.0, "b": 1.0, "d": c2}, {"a": 0.0, "b": 1.0, "d": c3}))
        n_gas, n_liq = -self.k_p * float(np.ravel(dp_hi)[0]), -self.k_c * float(np.ravel(dc_lo)[0])
        return dict(N=n_gas, N_gas_face=n_gas, N_liquid_face=n_liq, p2=p2, c2=c2,
                    p=p.ravel(), c=c.ravel(), nit=res.nit, success=bool(res.success),
                    backward_error=residual_check(fun, u)["backward_error"])

    def algebraic(self, p1, c3):
        """Whitman's Eq. 4 at this object's coefficients, with no grid."""
        return two_film_algebraic(p1, c3, self.k_p, self.k_c, self.peq)
