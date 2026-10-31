"""The Robeson upper bound as reusable functions (page H1.8).

The bound is ``P_i = k * alpha_ij**n`` with ``P_i`` the fast-gas permeability in
barrers, ``alpha_ij = P_i / P_j`` the separation factor and ``n < 0``. Take
``(k, n)`` from the ``present`` block of ``robeson-2008-table12-upper-bounds``.

No pymrm operator acts on this relation: it is algebra between two scalars.
"""

from __future__ import annotations

import numpy as np

__all__ = ["P_bound", "alpha_bound", "excursion", "envelope_fit"]


def P_bound(alpha, k, n):
    """Upper-bound permeability at a separation factor.

    Parameters
    ----------
    alpha : float or array_like
        Separation factor ``P_i / P_j``, dimensionless.
    k : float
        Front factor, barrer.
    n : float
        Slope of the log-log bound, dimensionless, negative.

    Returns
    -------
    numpy.ndarray
        ``k * alpha**n``, barrer.
    """
    return k * np.asarray(alpha, float) ** n


def alpha_bound(P, k, n):
    """Upper-bound separation factor at a permeability; the inverse of :func:`P_bound`.

    Parameters
    ----------
    P : float or array_like
        Fast-gas permeability, barrer.
    k, n : float
        Front factor (barrer) and slope of the bound.

    Returns
    -------
    numpy.ndarray
        ``(P / k)**(1 / n)``, dimensionless.
    """
    return (np.asarray(P, float) / k) ** (1.0 / n)


def excursion(P, alpha, k, n):
    """Ratio of a point's permeability to the bound's at the point's selectivity.

    With ``n < 0`` a point lies above the bound (it beats it) exactly when
    ``r > 1``: ``log alpha > (log P - log k)/n`` multiplied by the negative
    ``n`` gives ``log k + n log alpha < log P``.

    Parameters
    ----------
    P, alpha : float or array_like
        Permeability (barrer) and separation factor of the point(s).
    k, n : float
        The bound.

    Returns
    -------
    numpy.ndarray
        ``r = P / (k * alpha**n)``.
    """
    return np.asarray(P, float) / P_bound(alpha, k, n)


def envelope_fit(x, y, tol=1e-10):
    """Tightest covering line ``y = c + s x`` with every point on or below it.

    Minimises ``sum (y_i - c - s x_i)**2`` subject to ``y_i <= c + s x_i``. This
    is a convex quadratic programme in two unknowns, so its optimum has 0, 1 or 2
    active constraints and is the best feasible member of a finite list: the
    unconstrained least-squares line, for each point the line through it with
    the least-squares slope about that point, and the line through each pair of
    points. All slopes are computed from centred differences, never from normal
    equations in raw ``x``.

    Parameters
    ----------
    x, y : array_like
        Abscissae and ordinates (here ``log10 alpha`` and ``log10 P``).
    tol : float
        Feasibility tolerance, relative to ``max(1, max|y|)``.

    Returns
    -------
    k : float
        ``10**c``, the front factor of the covering line.
    s : float
        Its slope.
    sse : float
        Its sum of squared residuals.

    Raises
    ------
    RuntimeError
        If no candidate is feasible, which cannot happen for a finite set of
        points with at least two distinct abscissae.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    dxm, dym = x - x.mean(), y - y.mean()
    if not np.any(dxm != 0):
        raise RuntimeError("envelope_fit needs at least two distinct abscissae")
    s_ols = float(np.sum(dxm * dym) / np.sum(dxm * dxm))
    cands = [(float(y.mean() - s_ols * x.mean()), s_ols)]
    for i in range(len(x)):
        dx, dy = x - x[i], y - y[i]
        ss = float(np.sum(dx * dx))
        s1 = float(np.sum(dx * dy) / ss) if ss > 0 else 0.0
        cands.append((y[i] - s1 * x[i], s1))                     # point i active
        for j in range(i + 1, len(x)):
            if x[j] != x[i]:
                s2 = (y[j] - y[i]) / (x[j] - x[i])
                cands.append((y[i] - s2 * x[i], s2))             # points i and j active
    scale = tol * max(1.0, float(np.max(np.abs(y))))
    best = None
    for c, s in cands:
        r = y - (c + s * x)
        if float(np.max(r)) <= scale:
            sse = float(np.sum(r * r))
            if best is None or sse < best[0]:
                best = (sse, c, s)
    if best is None:
        raise RuntimeError("no feasible covering line")
    return 10.0 ** best[1], best[2], best[0]
