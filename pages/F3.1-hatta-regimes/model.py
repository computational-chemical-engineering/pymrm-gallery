"""The stagnant-film model of gas absorption with reaction (page F3.1).

Gas A dissolves into a liquid film that holds a non-volatile reactant B,
``A + nu B -> products``, first order in each. In ``xi = x/delta``,
``a = A/A_i`` and ``b = B/B_b``::

    a'' = Ha^2 a b,        b'' = (Ha^2/q) a b,        E_i = 1 + q,
    a(0) = 1, b'(0) = 0,   a(1) = 0, b(1) = 1,

and the enhancement factor is ``E = -a'(0)``. With B in excess (``b = 1``)
the problem is linear and ``E = Ha/tanh(Ha)`` exactly.

Contents:

- :class:`Film`, the film equations by finite volumes (pymrm operators; one
  sparse solve for the linear pseudo-first-order case, ``newton`` with
  ``NumJac`` for the second-order case);
- :func:`vkh` and :func:`decoursey`, the two published approximations;
- :func:`unique_root` and :func:`argext_log`, the root and extremum helpers.
"""
import numpy as np
from scipy.sparse import eye_array
from scipy.sparse.linalg import spsolve
from pymrm import NumJac, construct_div, construct_grad, newton
from pymrm.checks import find_roots


# ------------------------------------------------------------- root helpers
def unique_root(f, lo, hi, n_scan=16, log=False, xtol=1e-12):
    """The single root of ``f`` on the physical range ``[lo, hi]``.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        The physical range, scanned in full by ``pymrm.checks.find_roots``.
    n_scan : int
        Scan points.
    log : bool
        Scan on a logarithmic grid.
    xtol : float
        Absolute tolerance of the refined root.

    Returns
    -------
    float

    Raises
    ------
    RuntimeError
        Unless the scan finds exactly one root, or if ``f`` cannot be
        evaluated at a scan point.
    """
    found = find_roots(f, lo, hi, n_scan=n_scan, log=log, xtol=xtol)
    if found["failed"]:
        raise RuntimeError(f"root scan could not evaluate f at {found['failed']}")
    if len(found["roots"]) != 1:
        raise RuntimeError(f"{len(found['roots'])} roots on [{lo}, {hi}], "
                           f"expected one: {found['message']}")
    return float(found["roots"][0])


def argext_log(f, lo, hi, n_scan=25, h=1e-3, refine=True):
    """The maximum of ``|f|`` on ``[lo, hi]``, refined from a scan of the window.

    A log-spaced scan of the whole window brackets the largest ``|f|``; the
    extremum is then the single root, between the scan neighbours, of the
    central difference ``|f(x e^h)| - |f(x e^-h)|`` in ``ln x``.

    Parameters
    ----------
    f : callable
        Scalar function.
    lo, hi : float
        The window, scanned on a logarithmic grid.
    n_scan : int
        Scan points.
    h : float
        Step of the central difference in ``ln x``.
    refine : bool
        False returns the largest sampled value (a break-table defect).

    Returns
    -------
    tuple
        ``(x, f(x), interior)``. ``interior`` is False when the largest
        sampled ``|f|`` sits on a window end (the maximum is at or beyond it,
        and that end is returned) or when ``refine`` is False.

    Raises
    ------
    RuntimeError
        If the bracket around an interior scan maximum holds other than one
        root of the central difference.
    """
    xs = np.geomspace(lo, hi, n_scan)
    fs = np.array([f(x) for x in xs])
    i = int(np.argmax(np.abs(fs)))
    if not refine or i in (0, n_scan - 1):
        return float(xs[i]), float(fs[i]), False

    def dfd(lx):
        x = np.exp(lx)
        return abs(f(x * np.exp(h))) - abs(f(x * np.exp(-h)))

    lx = unique_root(dfd, np.log(xs[i - 1]), np.log(xs[i + 1]), n_scan=8)
    return float(np.exp(lx)), float(f(np.exp(lx))), True


# ---------------------------------------------------- the two approximations
def vkh(ha, e_i, denom="q"):
    """Van Krevelen-Hoftijzer: ``E = Ha s/tanh(Ha s)``, ``s = sqrt((E_i-E)/(E_i-1))``.

    Solved for ``E`` as the single root on ``[1, E_i]``.

    Parameters
    ----------
    ha, e_i : float
        Hatta number and instantaneous enhancement factor.
    denom : {"q", "ei"}
        ``"q"`` divides by ``E_i - 1`` as published; ``"ei"`` by ``E_i`` (a
        break-table mistranscription).

    Returns
    -------
    float
    """
    d = e_i - 1.0 if denom == "q" else e_i

    def f(E):
        s = np.sqrt(max((e_i - E) / d, 0.0))
        return E - (ha * s / np.tanh(ha * s) if ha * s > 1e-12 else 1.0)

    return unique_root(f, 1.0, e_i)


def decoursey(ha, e_i, q_power=2):
    """DeCoursey's explicit approximation.

    ``E = -Ha^2/(2q) + sqrt(Ha^4/(4q^2) + E_i Ha^2/q + 1)`` with ``q = E_i - 1``.

    Parameters
    ----------
    ha, e_i : float
        Hatta number and instantaneous enhancement factor.
    q_power : int
        Power of ``q`` under ``Ha^4``: 2 as published; 1 is a break-table
        mistranscription.

    Returns
    -------
    float
    """
    q = e_i - 1.0
    return -ha**2 / (2 * q) + np.sqrt(ha**4 / (4 * q**q_power) + e_i * ha**2 / q + 1.0)


# ------------------------------------------------------- the film, in pymrm
class Film:
    """Stagnant liquid film on ``xi`` in ``[0, 1]``, uniform grid.

    Parameters
    ----------
    n_x : int
        Number of cells.
    n_c : {1, 2}
        1: ``a`` alone, with ``b = 1`` (pseudo-first-order). 2: ``(a, b)``.
    a_interface : float
        Interface value ``a(0)``; 1 is the model (a break-table parameter).
    flux_read : {"face", "cell"}
        ``"face"`` reads ``E`` from the face-flux operator at ``xi = 0``;
        ``"cell"`` from the first cell centre by a one-sided difference (a
        break-table defect).
    """

    def __init__(self, n_x=400, n_c=2, a_interface=1.0, flux_read="face"):
        self.n_x, self.n_c, self.flux_read = n_x, n_c, flux_read
        self.a_interface = a_interface
        self.shape = (n_x, n_c)                # spatial axis first, then fields
        self.x_f = np.linspace(0.0, 1.0, n_x + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # a dc/dn + b c = d, n the OUTWARD normal (n = -xi at xi = 0, +xi at 1)
        if n_c == 2:
            # xi = 0, interface: a = a_i (saturated) and -db/dxi = 0 (B non-volatile)
            left = {"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[a_interface, 0.0]]}
            # xi = 1, bulk edge: a = 0 (depleted bulk) and b = 1 (bulk value)
            right = {"a": [[0.0, 0.0]], "b": [[1.0, 1.0]], "d": [[0.0, 1.0]]}
        else:
            left = {"a": 0.0, "b": 1.0, "d": a_interface}   # xi = 0: a = a_i
            right = {"a": 0.0, "b": 1.0, "d": 0.0}          # xi = 1: a = 0
        grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c,
                                           (left, right), axis=0)
        div_mat = construct_div(self.shape, self.x_f, nu=0, axis=0)   # slab
        # the constant diffusion operator, assembled once
        self.flux_mat, self.flux_bc = -grad_mat, -grad_bc
        self.jac_diff = div_mat @ self.flux_mat
        self.g_diff = (div_mat @ self.flux_bc).toarray().ravel()
        # E = -da/dxi on the interface face: row 0 of the face-flux operator
        # (face 0, field a), boundary term included
        self.iface = np.asarray(self.flux_mat[[0], :].todense()).ravel()
        self.iface_bc = float(self.flux_bc.toarray().ravel()[0])
        # (n_x, n_c), never a bare (n_x,): the rate is pointwise in the fields
        self.numjac = NumJac((n_x, n_c))

    def enhancement(self, u):
        """``E``, the interfacial flux of ``a`` (the purely diffusive one is 1).

        Parameters
        ----------
        u : numpy.ndarray
            Solution, shape ``(n_x, n_c)``.

        Returns
        -------
        float
        """
        u = np.asarray(u).reshape(self.shape)
        if self.flux_read == "cell":
            return float((self.a_interface - u[0, 0]) / self.x_c[0])
        return float(self.iface @ u.ravel() + self.iface_bc)

    def solve_first_order(self, ha):
        """``b = 1``: the problem is linear, so one sparse solve, no Newton.

        Parameters
        ----------
        ha : float
            Hatta number.

        Returns
        -------
        tuple
            ``(a, E)``: ``a`` at the cell centres, shape ``(n_x,)``.
        """
        if self.n_c != 1:
            raise ValueError("solve_first_order needs a Film with n_c = 1")
        a = self.jac_diff + ha**2 * eye_array(self.n_x, format="csc")
        u = spsolve(a, -self.g_diff)
        return u, self.enhancement(u)

    def solve_second_order(self, ha, e_i, tol=1e-11):
        """Steady film by Newton from the pseudo-first-order profile.

        Parameters
        ----------
        ha, e_i : float
            Hatta number and instantaneous enhancement factor, ``q = E_i - 1``.
        tol : float
            Newton step tolerance on unknowns of order one.

        Returns
        -------
        tuple
            ``(u, E)``: ``u = (a, b)`` at the cell centres, shape ``(n_x, 2)``.

        Raises
        ------
        RuntimeError
            If Newton does not converge.
        """
        if self.n_c != 2:
            raise ValueError("solve_second_order needs a Film with n_c = 2")
        q = e_i - 1.0

        def rate(u):
            r = np.clip(u[:, 0], 0, None) * np.clip(u[:, 1], 0, None)
            return np.stack([ha**2 * r, (ha**2 / q) * r], axis=1)

        g = min(ha, 700.0)                         # sinh overflows past ~710
        u0 = np.stack([self.a_interface * np.sinh(g * (1.0 - self.x_c)) / np.sinh(g),
                       np.ones_like(self.x_c)], axis=1)

        def residual(u):
            g_r, jac_r = self.numjac(rate, u)
            return (self.jac_diff @ u.reshape((-1, 1)) + self.g_diff.reshape((-1, 1))
                    + g_r.reshape((-1, 1))), self.jac_diff + jac_r

        res = newton(residual, u0, maxfev=200, tol=tol, bounds=(0.0, np.inf),
                     stagnation=True)
        if not res.success:
            raise RuntimeError(f"film solve failed at Ha={ha}, E_i={e_i}: {res.message}")
        u = np.asarray(res.x).reshape(self.shape)
        return u, self.enhancement(u)


def err_pct(ha, e_i, film, approx):
    """Signed error of an approximation against the film solution, in percent.

    Parameters
    ----------
    ha, e_i : float
        Hatta number and instantaneous enhancement factor.
    film : Film
        A two-field film.
    approx : callable
        ``approx(ha, e_i) -> E``.

    Returns
    -------
    float
        ``100 (E_approx - E_film)/E_film``.
    """
    E = film.solve_second_order(ha, e_i)[1]
    return 100.0 * (approx(ha, e_i) - E) / E
