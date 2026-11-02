"""The shrinking core: Yagi and Kunii's conversion-time relation and its ash shell.

Model code of gallery page B3.1 (Yagi and Kunii 1955). The conversion-time
relation of the paper (eq. 6) holds three resistances in series; the page
checks it against the integral of the resistances (eq. 5) and against its three
classical limits, and solves the quasi-steady gas profile in the ash shell by
pymrm.

Contents
--------
eq6, groups, theta_from_resistances, flow_analytic
    The closed forms: eq. (6), the groups (omega, gamma) from physical
    coefficients, the same history by integrating the moving boundary, and the
    series-resistance flow to a core of radius ``r_c``.
LAWS, limit_gaps
    The three classical single-resistance laws and the distance of eq. (6) from
    each in its limit.
regime_map
    For a grid of (omega, gamma), the closest single-resistance law and its
    largest error in the conversion time.
single_law_error, crest, resistance_shares
    The same error at one point, its maximum over gamma at a given omega (the
    crest of the map, a kink found by a bounded 1D maximisation), and the
    shares of the three resistances.
Shell
    The gas profile in the shell ``r_c < r < R`` by pymrm, and the molar flows
    to the core and from the film, read on the faces.

Several functions take an argument whose default is the correct model and
whose alternative is a defect the page's break table injects.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import spsolve

from pymrm import compute_boundary_values, construct_div, construct_grad

__all__ = ["eq6", "groups", "theta_from_resistances", "flow_analytic", "LAWS", "limit_gaps",
           "regime_map", "single_law_error", "crest", "resistance_shares", "Shell"]


# ----------------------------------------------------------------- closed forms
def eq6(x, omega, gamma, cubic="printed"):
    """Yagi and Kunii eq. (6): ``theta / theta_B`` against ``x = r / R``.

    Parameters
    ----------
    x : array_like
        Core radius over particle radius, ``r / R``.
    omega, gamma : float
        ``k_c1 / k_d1`` and ``3 k_f1 / k_d1``.
    cubic : {"printed", "misread"}
        ``"printed"`` is ``omega (1 - 2 gamma)`` in the cubic term of the
        numerator; ``"misread"`` is ``omega (1 - gamma)``, a defect.

    Returns
    -------
    numpy.ndarray
        ``theta / theta_B``, zero at ``x = 1`` and one at ``x = 0``.
    """
    x = np.asarray(x, float)
    c3 = omega * (1.0 - 2.0 * gamma) if cubic == "printed" else omega * (1.0 - gamma)
    num = gamma * x + 3.0 * omega * gamma * x**2 + c3 * x**3
    return 1.0 - num / (omega + omega * gamma + gamma)


def groups(kc, kf, D, R, film_factor=3.0, kd_coef=12.0):
    """The groups ``(omega, gamma)`` from physical coefficients.

    Uses eq. (5) with ``k_d1 = kd_coef * D / D_p`` and ``D_p = 2 R``.

    Parameters
    ----------
    kc, kf, D, R : float
        Reaction and film coefficients, ash diffusivity, particle radius.
    film_factor : float
        The 3 in ``gamma = 3 k_f1 / k_d1``; another value is a defect.
    kd_coef : float
        The 12 in ``k_d1 = 12 D / D_p``; another value is a defect.

    Returns
    -------
    omega, gamma : float
    """
    kd = kd_coef * D / (2.0 * R)
    return kc / kd, film_factor * kf / kd


def theta_from_resistances(x, kc, kf, D, R):
    """``theta / theta_B`` by integrating the moving boundary over the resistances.

    Equates the consumption of the core to the series-resistance flow
    (:func:`flow_analytic`) and integrates in the core radius; it never uses
    eq. (6).

    Parameters
    ----------
    x : array_like
        ``r / R``.
    kc, kf, D, R : float
        Reaction and film coefficients, ash diffusivity, particle radius.

    Returns
    -------
    numpy.ndarray
        ``theta / theta_B``.
    """
    x = np.asarray(x, float)
    num = (1 - x) / kc + (R / D) * ((1 - x**2) / 2 - (1 - x**3) / 3) + (1 - x**3) / (3 * kf)
    den = 1 / kc + R / (6 * D) + 1 / (3 * kf)
    return num / den


def flow_analytic(rc, R, D, kc, kf, C1=1.0):
    """Molar flow to a core of radius ``rc``: three resistances in series.

    Parameters
    ----------
    rc, R, D, kc, kf, C1 : float
        Core and particle radius, ash diffusivity, reaction and film
        coefficients, bulk concentration.

    Returns
    -------
    float
        ``W = 4 pi C1 / (1/(kc rc^2) + (1/rc - 1/R)/D + 1/(kf R^2))``.
    """
    return 4.0 * np.pi * C1 / (1 / (kc * rc**2) + (1 / rc - 1 / R) / D + 1 / (kf * R**2))


LAWS = {
    "film": lambda X: X,
    "ash": lambda X: 1 - 3 * (1 - X) ** (2 / 3) + 2 * (1 - X),
    "reaction": lambda X: 1 - (1 - X) ** (1 / 3),
}
"""The classical conversion-time laws, ``theta / theta_B`` against ``X``."""


def limit_gaps(x, eps_film=1e-9, big_ash=1e7, eps_reaction=1e-9):
    """Largest distance of eq. (6) from each classical law in its limit.

    Parameters
    ----------
    x : array_like
        ``r / R`` values; the conversion is ``X = 1 - x**3``.
    eps_film : float
        ``gamma`` of the film limit (``omega = 1``).
    big_ash : float
        ``omega`` of the ash limit (``gamma = 3 omega``).
    eps_reaction : float
        ``omega`` of the reaction limit (``gamma = 1``).

    Returns
    -------
    dict
        ``{"film": ..., "ash": ..., "reaction": ...}``, each ``max |eq6 - law|``.
    """
    x = np.asarray(x, float)
    X = 1 - x**3
    return {
        "film": float(np.abs(eq6(x, 1.0, eps_film) - LAWS["film"](X)).max()),
        "ash": float(np.abs(eq6(x, big_ash, 3 * big_ash) - LAWS["ash"](X)).max()),
        "reaction": float(np.abs(eq6(x, eps_reaction, 1.0) - LAWS["reaction"](X)).max()),
    }


def regime_map(om_g, ga_g, n_x=300, laws=LAWS, cubic="printed"):
    """The closest single-resistance law and its error over a grid of (omega, gamma).

    Parameters
    ----------
    om_g, ga_g : array_like
        Values of omega (columns) and gamma (rows).
    n_x : int
        Points in ``x = r / R`` from 0.001 to 1 for the maximum over the history.
    laws : dict
        ``{name: f(X)}``; defaults to :data:`LAWS`, another dict is a defect.
    cubic : {"printed", "misread"}
        Passed to :func:`eq6`.

    Returns
    -------
    best : numpy.ndarray of int
        Index in ``laws`` of the closest law, shape ``(len(ga_g), len(om_g))``.
    err : numpy.ndarray
        Largest ``|eq6 - law|`` of that law, same shape.
    """
    xs = np.linspace(0.001, 1.0, n_x)
    Xs = 1 - xs**3
    ref = [f(Xs) for f in laws.values()]
    best = np.zeros((len(ga_g), len(om_g)), int)
    err = np.zeros(best.shape)
    for i, ga in enumerate(ga_g):
        for j, om in enumerate(om_g):
            full = eq6(xs, om, ga, cubic)
            e = [np.abs(full - r).max() for r in ref]
            best[i, j] = int(np.argmin(e))
            err[i, j] = min(e)
    return best, err


def single_law_error(omega, gamma, n_x=20001, laws=LAWS):
    """Error of the closest single-resistance law at one (omega, gamma).

    Parameters
    ----------
    omega, gamma : float
        The two groups of eq. (6).
    n_x : int
        Points in ``x = r / R`` from 0 to 1 for the maximum over the history.
    laws : dict
        ``{name: f(X)}``; defaults to :data:`LAWS`.

    Returns
    -------
    float
        ``min`` over the laws of ``max |eq6 - law|``.
    """
    xs = np.linspace(0.0, 1.0, n_x)
    Xs = 1 - xs**3
    full = eq6(xs, omega, gamma)
    return min(float(np.abs(full - f(Xs)).max()) for f in laws.values())


def crest(omega, gamma_range=(0.1, 10.0), laws=LAWS, method="bounded", n_scan=400, n_x=20001):
    """Largest error of the best single law over gamma at a given omega.

    The best-law error is a maximum of minima: where two laws tie it has a kink,
    which a grid in gamma undersamples. ``"bounded"`` scans ``n_scan`` values of
    gamma on a log scale, then maximises with
    ``scipy.optimize.minimize_scalar(method="bounded")`` between the neighbours
    of the best one; ``"grid"`` returns the scan's own maximum (a defect).

    Parameters
    ----------
    omega : float
        The group ``k_c1 / k_d1``.
    gamma_range : tuple of float
        Window in ``gamma`` to search.
    laws : dict
        ``{name: f(X)}``; defaults to :data:`LAWS`, another dict is a defect.
    method : {"bounded", "grid"}
        See above.
    n_scan : int
        Points of the scan.
    n_x : int
        Points in ``x`` of :func:`single_law_error`.

    Returns
    -------
    gamma, error : float
        Location and height of the maximum.
    """
    lg = np.linspace(np.log10(gamma_range[0]), np.log10(gamma_range[1]), n_scan)
    val = np.array([single_law_error(omega, 10.0**t, n_x, laws) for t in lg])
    k = int(np.argmax(val))
    if method == "grid":
        return float(10.0**lg[k]), float(val[k])
    res = minimize_scalar(lambda t: -single_law_error(omega, 10.0**t, n_x, laws),
                          bounds=(lg[max(k - 1, 0)], lg[min(k + 1, n_scan - 1)]),
                          method="bounded", options={"xatol": 1e-12})
    return float(10.0**res.x), float(-res.fun)


def resistance_shares(omega, gamma):
    """Shares of the three resistances of eq. (5) in their sum.

    Parameters
    ----------
    omega, gamma : float
        ``k_c1 / k_d1`` and ``3 k_f1 / k_d1``; the resistances are in the ratio
        ``1/omega : 1 : 1/gamma`` (reaction, ash layer, film).

    Returns
    -------
    reaction, ash, film : float
        Fractions summing to one.
    """
    r = np.array([1.0 / omega, 1.0, 1.0 / gamma])
    r = r / r.sum()
    return float(r[0]), float(r[1]), float(r[2])


# ------------------------------------------------------------------- ash shell
class Shell:
    """Quasi-steady gas in the ash shell ``r_c < r < R`` by pymrm finite volume.

    ``div(-D grad C) = 0`` on ``n`` equal cells in a spherical shell, with a
    Robin condition on each face, written with the outward normal of the shell.
    At the core face the outward normal points to the centre and the reaction
    removes gas: ``D dC/dn + k_c C = 0``. At the outer face the film delivers
    it: ``D dC/dn + k_f C = k_f C1``. The operators are assembled once.

    Parameters
    ----------
    rc, R : float
        Core and particle radius.
    D, kc, kf : float
        Ash diffusivity, reaction and film coefficients.
    n : int
        Cells.
    C1 : float
        Bulk concentration.
    nu : int
        Geometry of ``construct_div``: 2 is spherical; another value is a defect.
    inner_normal : {"outward", "inward"}
        ``"inward"`` writes the core condition with the inward normal, a defect.
    film : bool
        ``False`` imposes the bulk concentration at the outer face, dropping the
        film resistance (a defect).
    """

    def __init__(self, rc, R, D, kc, kf, n, C1=1.0, nu=2, inner_normal="outward", film=True):
        self.rc, self.R, self.D, self.kc, self.kf, self.n, self.C1 = rc, R, D, kc, kf, n, C1
        self.nu = nu
        self.x_f = np.linspace(rc, R, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        a_in = D if inner_normal == "outward" else -D
        self.bc = (
            # core face, outward normal toward the centre: D dC/dn + kc C = 0
            {"a": a_in, "b": kc, "d": 0.0},
            # outer face, outward normal: D dC/dn + kf C = kf C1 (film), or C = C1
            {"a": D, "b": kf, "d": kf * C1} if film else {"a": 0.0, "b": 1.0, "d": C1},
        )
        shape = (n, 1)
        grad, grad_bc = construct_grad(shape, self.x_f, self.x_c, self.bc)
        div = construct_div(shape, self.x_f, nu=nu)
        self._lhs = (div @ (-D * grad)).tocsc()
        self._rhs = -(div @ (-D * grad_bc)).toarray().ravel()

    def solve(self, permc_spec="COLAMD"):
        """Concentration in the cells.

        Parameters
        ----------
        permc_spec : str
            Column ordering of the sparse factorisation (``scipy.sparse.linalg.spsolve``).

        Returns
        -------
        numpy.ndarray
            ``C`` at the cell centres, shape ``(n,)``.
        """
        return spsolve(self._lhs, self._rhs, permc_spec=permc_spec)

    def faces(self, C):
        """Concentration and radial gradient on the two faces.

        Parameters
        ----------
        C : numpy.ndarray
            Cell concentrations.

        Returns
        -------
        C_core, dCdr_core, C_outer, dCdr_outer : float
            From ``compute_boundary_values`` with the shell's own boundary
            conditions; the gradients are along ``+r``.
        """
        v0, g0, v1, g1 = compute_boundary_values(C.reshape(-1, 1), self.x_f, self.x_c, self.bc)
        return float(v0.ravel()[0]), float(g0.ravel()[0]), float(v1.ravel()[0]), float(g1.ravel()[0])

    def flow_core(self, C, read="face"):
        """Molar flow reaching the core, ``4 pi r^2 D dC/dr`` at the core.

        Parameters
        ----------
        C : numpy.ndarray
            Cell concentrations.
        read : {"face", "cells", "cells_at_core"}
            ``"face"``: the discrete gradient on the core face. ``"cells"``:
            the difference of the first two cell centres, at their midpoint
            radius. ``"cells_at_core"``: the same difference at the core radius
            (a defect: the flow is read at the wrong radius).

        Returns
        -------
        float
        """
        area = 4.0 * np.pi
        if read == "face":
            return area * self.rc**2 * self.D * self.faces(C)[1]
        dCdr = (C[1] - C[0]) / (self.x_c[1] - self.x_c[0])
        r = 0.5 * (self.x_c[0] + self.x_c[1]) if read == "cells" else self.rc
        return area * r**2 * self.D * dCdr

    def flow_outer(self, C, read="face"):
        """Molar flow delivered by the film, ``4 pi R^2 k_f (C1 - C_surface)``.

        Parameters
        ----------
        C : numpy.ndarray
            Cell concentrations.
        read : {"face", "cell"}
            Surface concentration from the outer face value, or the last cell
            (a defect: the surface concentration read one half cell inside).

        Returns
        -------
        float
        """
        c_s = self.faces(C)[2] if read == "face" else float(C[-1])
        return 4.0 * np.pi * self.R**2 * self.kf * (self.C1 - c_s)
