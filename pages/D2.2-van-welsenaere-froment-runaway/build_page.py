#!/usr/bin/env python3
"""Generate index.ipynb for page D2.2. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "Parametric sensitivity and runaway in a fixed bed"
description: "Two 1970 criteria for when a cooled tubular reactor stops being controllable — rebuilt, then swept over the whole operating plane the original could only sample."
categories: [sec:D, struct:S2, tier:T0, data:tier6, phase:gas-solid]
date: 2026-07-27
---

# Parametric sensitivity and runaway in a fixed bed

**Catalog ID:** `D2.2` · **Structures:** `S2` (plug flow with reaction) · **Tier:** T0

An exothermic reaction in a cooled tube always makes a hot spot. Raise the inlet
concentration a few percent and the hot spot may rise a few degrees — or it may
rise two hundred. There is a boundary, and Van Welsenaere and Froment found two
ways to locate it without integrating anything."""))

cells.append(md(r"""## Background

Put an exothermic reaction in a tube with a cooled wall and the temperature
profile has a maximum. Where that maximum sits is a design question with
consequences: too hot and the catalyst sinters, the selectivity collapses, or
the tube itself creeps. What makes it hard is that the hot spot is not merely
sensitive to the operating variables — over part of the range it is *violently*
sensitive, and a 1 % change in inlet partial pressure can move it by a hundred
degrees.

Before 1970 there were two ways to deal with this. Bilous and Amundson (1956)
perturbed a known steady state, which requires having integrated the reactor
first. Barkelew (1959) integrated a very large number of cases and condensed the
results into an empirical chart, which requires trial and error to use and is
tied to his modified rate expression.

Van Welsenaere and Froment did something different: they looked for **intrinsic**
features of the solution — characteristic points that exist in the temperature
profile whether or not anyone has computed it — and translated them into the
*p*–*T* phase plane, where they become geometry. The result is a set of closed-form
inequalities on the inlet conditions. No computer, which in 1970 was the point.

Note the terminology, which they are careful about: with gas and solid at the
same temperature, an ideal tubular reactor is always stable in the strict sense.
There is no bifurcation here. "Runaway" is a statement about *sensitivity*, not
about multiplicity — which is what distinguishes this page from
[`B1.1`](../B1.1-thiele-weisz-hicks/), where the multiplicity is real."""))

cells.append(md(r"""## The published model

**Reactor** (their Eqs. 3–4). One dimension, pseudo-homogeneous, constant wall
temperature, single irreversible pseudo-first-order reaction:

$$
\frac{\mathrm{d}p}{\mathrm{d}z} = -A\,p\,\mathrm{e}^{-a/T+b},
\qquad
\frac{\mathrm{d}T}{\mathrm{d}z} = B\,p\,\mathrm{e}^{-a/T+b} - C\,(T-T_w),
$$

$$
A = \frac{M P \rho_b}{\rho_g} p_B^0, \qquad
B = \frac{(-\Delta H)\rho_b}{c_p} p_B^0, \qquad
C = \frac{2U}{c_p R}, \qquad z = \frac{z'}{u},
$$

with $p = p^0$ and $T = T_0 = T_w$ at $z = 0$. The axial coordinate is a contact
time, not a length.

**The *p*–*T* plane.** Dividing one equation by the other removes $z$
altogether, so every reactor is a single trajectory in the $(T, p)$ plane. Two
loci live in that plane and neither depends on the inlet conditions:

- the **maxima curve** (Eq. 7), where $\mathrm{d}T/\mathrm{d}z = 0$:
  $p_m = (T_m - T_w) \big/ \left[\tfrac{B}{C}\mathrm{e}^{-a/T_m+b}\right]$;
- the **$p_s$ curve** (Eq. 18), the locus of the inflexion points that appear
  *before* the maximum: $p_s = T_i^2 \big/ \left[a\tfrac{B}{C}\mathrm{e}^{-a/T_i+b}\right]$.

**First criterion.** The maxima curve has a maximum of its own, at

$$
T_M = \tfrac{1}{2}\left[a - \sqrt{a\,(a - 4T_w)}\right] \qquad (8)
$$

A trajectory that crosses the maxima curve *beyond* its peak is one where a
further small increase in $p^0$ moves the hot spot sharply. So: **the trajectory
through the maximum of the maxima curve is critical.**

**Second criterion.** An inflexion point *before* the maximum is what makes the
hot spot grow. So: **the critical trajectory is the one tangent to the $p_s$
curve.** The tangency condition reduces to one implicit equation in
$t = a/T_i$ (Eq. 20),

$$
K = \frac{t-2}{t\,\mathrm{e}^{-t+20}}\left[1 - t\left(1 - \frac{t}{t_w}\right)\right],
\qquad K = \frac{A}{C}\mathrm{e}^{b-20}, \quad t_w = \frac{a}{T_w},
$$

which the paper solves *graphically*, through its Figs. 6 and 7, to avoid
needing a computer.

**From the critical point to the inlet.** Both criteria give a critical
$(T_{cr}, p_{cr})$ inside the reactor; getting back to the inlet needs an
extrapolation, and the paper gives a lower and an upper limit rather than one
value:

$$
(p^0)_l = \frac{T_{cr}^2}{a\frac{B}{C}\mathrm{e}^{-a/T_{cr}+b}} + \frac{A}{B}(T_{cr}-T_0)
\quad (27), \qquad
(p^0)_u = \frac{A}{B}(T_{cr}-T_w)\left[\frac{1}{\sqrt{X}} + 1\right]^2 \quad (28),
$$

with $X = \frac{A}{C}\mathrm{e}^{-a/T_{cr}+b}$. Their empirical observation — the
one that makes the method usable — is that the **mean of the two limits is an
excellent estimate of the true critical inlet pressure**, and for the first
criterion that mean has the closed form (30)

$$
(p^0)_m = \frac{A}{B}(T_M-T_w)\left[1 + \frac{1}{\sqrt{X}} + \frac{1}{X}\right].
$$"""))

cells.append(md(r"""## Parameters and assumptions

**Assumptions:** one dimension, no axial dispersion; pseudo-homogeneous, so gas
and solid share a temperature; constant wall temperature; constant fluid
density and velocity; a single irreversible reaction, pseudo-first-order in the
limiting reactant with the second reactant in excess at $p_B^0$; ideal gas.

**One trap in this paper, and it is a units trap.** $c_p$ is printed as
kcal/m³·°C — it is already a *volumetric* heat capacity. Multiplying by
$\rho_g$ "to fix the units" changes $C$ by a factor of 1.293 and moves the
runaway boundary with no other symptom. $\rho_g$ appears only in $A$.

**One trap in getting the numbers out.** This is a 1970 scan and the Elsevier
full-text API returns the publisher's OCR of it, which discards the mid-dot
decimal separator: `R = 00125 m` for 0.0125 m, `b = 19837` for 19.837,
`001353` for 0.01353 atm. The API is excellent for prose and useless for
numbers. Everything below was read from a 600 dpi render of the printed page,
and the validation section shows the check that the reading is right."""))

cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code('''import sys, urllib.request
from pathlib import Path

if not any("shared" in p for p in sys.path):
    local = Path.cwd()
    for _ in range(4):
        if (local / "shared" / "gallery_utils.py").is_file():
            sys.path.insert(0, str(local / "shared")); break
        local = local.parent
    else:
        url = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
               "pymrm-gallery/main/shared/gallery_utils.py")
        urllib.request.urlretrieve(url, "gallery_utils.py")
        sys.path.insert(0, ".")

import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse import SparseEfficiencyWarning
from scipy.sparse.linalg import MatrixRankWarning
from pymrm import construct_convflux_upwind, construct_div, NumJac, newton
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "D2.2-van-welsenaere-froment-runaway"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

cells.append(md(r"""## The data

Two datasets. The **parameters** are the base case of their Section 1. The
**examples** are every number printed in their Section 6 — four worked problems,
including the intermediate steps, so a reimplementation can be checked stage by
stage rather than only at the answer.

**These are not measurements.** The paper contains no experimental data at all;
its numbers are partly closed-form extrapolations and partly its own
fourth-order Runge–Kutta integrations. This page is therefore provenance tier 6,
validated against a published reference solution in the same sense as
[`B1.1`](../B1.1-thiele-weisz-hicks/) and [`F3.1`](../F3.1-hatta-regimes/). The
`kind` column separates the two, because a criterion value and an integration
value test different things and should not be pooled."""))

cells.append(code('''par = load_data("van-welsenaere-froment-1970-parameters.csv", page=PAGE)
ex = load_data("van-welsenaere-froment-1970-examples.csv", page=PAGE)
ex_meta = load_meta("van-welsenaere-froment-1970-examples.csv", page=PAGE)

P = dict(zip(par.symbol, par.value))
print(par.to_string(index=False))
print()
print(ex.groupby(["example", "kind"]).size().unstack(fill_value=0).to_string())
print(f"\\n{cite_data(ex_meta)}")
print(f"provenance tier {ex_meta['provenance_tier']['tier']} — "
      + " ".join(ex_meta["provenance_tier"]["note"].split())[:96] + " ...")'''))

cells.append(md(r"""## PyMRM implementation

Two pieces. The **reactor** is convection with a source in the contact time
$z$ — `construct_convflux_upwind` then `construct_div`, the same two calls as
[`C2.1`](../C2.1-xu-froment-smr/), with the state now $(p, T)$ instead of two
conversions. The **criteria** are algebra on top of it and touch no operator at
all.

Two implementation choices are worth stating.

*The exponential is guarded.* Beyond the critical inlet pressure the true
solution has no bounded hot spot, and a Newton step can walk $T$ somewhere that
overflows $\mathrm{e}^{-a/T+b}$. Clipping the exponent keeps the residual finite
so the solver reports failure instead of raising, and a failed solve is exactly
the signal that the case ran away.

*The critical inlet pressure is found by bisection on cold-started solves*, not
by walking a continuation curve. That is deliberate: a continuation sweep makes
the reported number depend on the path taken to reach it, which is the defect
CI caught on [`B1.1`](../B1.1-thiele-weisz-hicks/)."""))

cells.append(code('''A = P["M"] * P["P"] * P["rho_b"] / P["rho_g"] * P["p_B0"]
B = P["minus_dH"] * P["rho_b"] / P["c_p"] * P["p_B0"]
C_BASE = 2.0 * P["U"] / (P["c_p"] * P["R"])
A_EXP, B_EXP = P["a"], P["b"]


def C_of_R(R):
    """Eq. 4's cooling group. c_p is already volumetric -- see the units trap."""
    return 2.0 * P["U"] / (P["c_p"] * R)


def kexp(T):
    """exp(-a/T + b), guarded so a diverging Newton step cannot overflow."""
    return np.exp(np.clip(-A_EXP / np.maximum(T, 1.0) + B_EXP, -700.0, 50.0))


# ---- the two loci in the p-T plane, neither depending on the inlet ----------
def maxima_curve(T, Tw, C=C_BASE):
    """Eq. 7: where dT/dz = 0."""
    return (T - Tw) / (B / C * kexp(T))


def ps_curve(T, C=C_BASE):
    """Eq. 18: the locus of inflexion points before the maximum."""
    return T ** 2 / (A_EXP * B / C * kexp(T))


# ---- criterion 1 ------------------------------------------------------------
def T_M(Tw):
    """Eq. 8: the maximum of the maxima curve."""
    return 0.5 * (A_EXP - np.sqrt(A_EXP * (A_EXP - 4.0 * Tw)))


# ---- criterion 2 ------------------------------------------------------------
def T_cr_tangency(Tw, C=C_BASE):
    """Solve their Eq. 20 for t = a/(T_i)_t directly.

    The paper reaches this through Figs. 6 and 7 explicitly so that no computer
    is needed. With one, `brentq` on a scanned bracket is exact and their two
    figures become things to check rather than things to read.
    """
    K = A / C * np.exp(B_EXP - 20.0)
    tw = A_EXP / Tw

    def f(t):
        return (t - 2.0) / (t * np.exp(-t + 20.0)) * (1.0 - t * (1.0 - t / tw)) - K

    ts = np.linspace(A_EXP / T_M(Tw) * 0.80, tw - 1e-9, 4000)
    sign = np.nonzero(np.diff(np.sign(f(ts))) != 0)[0]
    if len(sign) == 0:
        return np.nan
    i = sign[-1]
    return A_EXP / brentq(f, ts[i], ts[i + 1])


# ---- from the critical point back to the inlet ------------------------------
def inlet_limits(T_cr, Tw, C=C_BASE, T0=None):
    """Eqs. 27, 28 and their mean.

    Eq. 27 is the general lower limit. At the first criterion's T_M it coincides
    with their Eq. 29, because Eq. 8 makes T_M^2 = a (T_M - T_w) an identity.
    """
    T0 = Tw if T0 is None else T0
    X = A / C * kexp(T_cr)
    lower = T_cr ** 2 / (A_EXP * B / C * kexp(T_cr)) + A / B * (T_cr - T0)
    upper = A / B * (T_cr - Tw) * (1.0 / np.sqrt(X) + 1.0) ** 2
    return lower, upper, 0.5 * (lower + upper)


def back_integrate(T_cr, p_cr, Tw, C=C_BASE, p_top=0.3):
    """The paper's own definition of the TRUE critical inlet pressure.

    Dividing Eq. 4 by Eq. 3 removes z and leaves dT/dp (their Eq. 5), so the
    trajectory through the critical point can be followed straight back to the
    inlet, which is where it reaches T = T_w. Note the direction: p RISES as we
    go back upstream, because p falls monotonically along the reactor.

    This is a phase-plane quadrature, not a reactor discretisation, so it is a
    genuinely independent check on the pymrm solve above rather than the same
    computation twice.
    """
    def dT_dp(p, T):
        return [-B / A + (C / A) * (T[0] - Tw) / (p * kexp(T[0]))]

    reached_inlet = lambda p, T: T[0] - Tw
    reached_inlet.terminal, reached_inlet.direction = True, -1
    sol = solve_ivp(dT_dp, (p_cr, p_top), [T_cr], events=reached_inlet,
                    rtol=1e-11, atol=1e-12, max_step=p_cr / 500)
    return float(sol.t_events[0][0]) if len(sol.t_events[0]) else np.nan


def critical_points(Tw, C=C_BASE):
    """Critical (T, p) inside the reactor for each criterion, and p0 for each."""
    TM = T_M(Tw)
    p1 = maxima_curve(TM, Tw, C)
    Tt = T_cr_tangency(Tw, C)
    p2 = ps_curve(Tt, C)
    return dict(T1=TM, p1=p1, p0_1=back_integrate(TM, p1, Tw, C),
                T2=Tt, p2=p2, p0_2=back_integrate(Tt, p2, Tw, C))'''))

cells.append(code('''class RunawayTube:
    """1-D pseudo-homogeneous tube with a constant wall temperature.

    State layout (n_z, 2): [..., 0] = p (atm), [..., 1] = T (K), against the
    contact time z = z'/u in hours, which is the paper's own coordinate.
    """

    def __init__(self, p0, Tw, T0=None, R=None, length=1.0, n_z=1500):
        self.p0, self.Tw = p0, Tw
        self.T0 = Tw if T0 is None else T0
        self.C = C_BASE if R is None else C_of_R(R)
        self.z_f = np.linspace(0.0, length / P["u"], n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.shape = (n_z, 2)
        # outward normal, so both dicts read a.dx/dn + b.x = d
        # inlet : p = p0, T = T0        -> a=0, b=1, d=(p0, T0)
        # outlet: dx/dn = 0, outflow    -> a=1, b=0, d=0
        self.bc = ({"a": 0.0, "b": 1.0, "d": np.array([p0, self.T0])},
                   {"a": 1.0, "b": 0.0, "d": 0.0})
        self.u = np.tile(np.array([p0, self.T0]), (n_z, 1))
        self.diverged = False
        self._build_operators()

    def _build_operators(self):
        # v = 1 by construction: z is the contact time, so d(v x)/dz = source.
        # nu = 0: the contact time is a Cartesian coordinate.
        conv_mat, conv_bc = construct_convflux_upwind(
            self.shape, self.z_f, self.z_c, self.bc, v=1.0, axis=0)
        div_mat = construct_div(self.shape, self.z_f, nu=0, axis=0)
        self.jac_const = div_mat @ conv_mat
        self.g_const = div_mat @ conv_bc
        self.numjac = NumJac(self.shape)   # pointwise source: last axis only

    def reaction(self, u):
        p, T = u[..., 0], u[..., 1]
        rate = p * kexp(T)
        return np.stack([-A * rate, B * rate - self.C * (T - self.Tw)], axis=-1)

    def residual(self, u):
        g_rxn, jac_rxn = self.numjac(self.reaction, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self, maxfev=200):
        try:
            result = newton(self.residual, self.u, maxfev=maxfev)
        except (ValueError, RuntimeError):
            self.diverged = True
            return None
        self.u = result.x.reshape(self.shape)
        self.diverged = not np.all(np.isfinite(self.u))
        return result

    @property
    def p(self):
        return self.u[:, 0]

    @property
    def T(self):
        return self.u[:, 1]


def hot_spot(p0, Tw, **kw):
    """Hot spot temperature, or +inf when the case runs away.

    A runaway case has no bounded solution, so the Newton solve is *expected*
    to fail there; the singular-Jacobian warning that goes with it is the
    signal, not a problem, and is silenced so it does not bury the output.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MatrixRankWarning)
        warnings.simplefilter("ignore", SparseEfficiencyWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        m = RunawayTube(p0, Tw, **kw)
        m.solve()
    return np.inf if m.diverged else float(m.T.max())


def critical_p0(Tw, target=None, lo=1e-4, hi=0.05, tol=1e-7, **kw):
    """Inlet pressure whose hot spot reaches `target` (default T_M).

    Plain bisection on cold-started solves. No continuation, so the answer does
    not depend on the path taken to reach it.
    """
    target = T_M(Tw) if target is None else target
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if hot_spot(mid, Tw, **kw) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)'''))

cells.append(md(r"""## Results

The paper's Figs. 1 and 2 are the sensitivity itself: a set of profiles at inlet
pressures that differ by a few percent. Below, the same sweep, with the critical
value marked."""))

cells.append(code('''Tw = 625.0
p_crit = critical_p0(Tw)
levels = [0.008, 0.012, 0.0145, 0.0158, 0.0164, p_crit, 0.0168, 0.0172]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
cmap = plt.cm.viridis(np.linspace(0.05, 0.92, len(levels)))
for p0, col in zip(levels, cmap):
    m = RunawayTube(p0, Tw)
    m.solve()
    if m.diverged:
        continue
    z_m = m.z_c * P["u"]                       # back to metres for the plot
    crit = abs(p0 - p_crit) < 1e-9
    kw = dict(color="tab:red" if crit else col, lw=2.2 if crit else 1.4,
              zorder=5 if crit else 2)
    axes[0].plot(z_m, m.p, **kw)
    axes[1].plot(z_m, m.T, label=f"{p0*1000:.2f}" + (" (critical)" if crit else ""),
                 **kw)
axes[0].set(xlabel="z' (m)", ylabel="p (atm)", xlim=(0, 1))
axes[1].set(xlabel="z' (m)", ylabel="T (K)", xlim=(0, 1))
axes[1].axhline(T_M(Tw), color="k", ls="--", lw=0.9)
axes[1].text(0.985, T_M(Tw) + 0.7, r"$T_M$ (Eq. 8)", fontsize=8, ha="right")
axes[1].legend(title=r"$p^0$ (10$^{-3}$ atm)", fontsize=8, loc="lower right",
               ncol=2)
fig.suptitle(f"Parametric sensitivity at $T_w$ = {Tw:.0f} K "
             f"(their Figs. 1 and 2, rebuilt)", fontsize=11)
fig.tight_layout()
plt.show()

print(f"critical inlet partial pressure at Tw = {Tw:.0f} K: {p_crit:.5f} atm")
print(f"paper, by back-integration from the critical point:  0.01651 atm")'''))

cells.append(md(r"""The whole argument is easier to see in the *p*–*T* plane,
where the trajectories are inlet-independent geometry against two fixed
curves."""))

cells.append(code('''fig, ax = plt.subplots(figsize=(6.6, 4.8))
Tgrid = np.linspace(Tw + 0.4, 700.0, 400)
ax.plot(Tgrid, maxima_curve(Tgrid, Tw), "k-", lw=1.6, label="maxima curve (Eq. 7)")
ax.plot(Tgrid, ps_curve(Tgrid), "k--", lw=1.4, label=r"$p_s$ curve (Eq. 18)")

TM = T_M(Tw)
ax.plot([TM], [maxima_curve(TM, Tw)], "ko", ms=7, zorder=6)
ax.annotate("criterion 1:\\nmaximum of the\\nmaxima curve",
            (TM, maxima_curve(TM, Tw)), textcoords="offset points",
            xytext=(26, -6), fontsize=8,
            arrowprops=dict(arrowstyle="-", lw=0.7, color="k"))
Tt = T_cr_tangency(Tw)
ax.plot([Tt], [ps_curve(Tt)], "s", color="tab:purple", ms=7, zorder=6)
ax.annotate("criterion 2:\\ntangency with $p_s$",
            (Tt, ps_curve(Tt)), textcoords="offset points",
            xytext=(-58, 42), fontsize=8, color="tab:purple",
            arrowprops=dict(arrowstyle="-", lw=0.7, color="tab:purple"))

for p0, col in zip(levels, cmap):
    m = RunawayTube(p0, Tw)
    m.solve()
    if m.diverged:
        continue
    crit = abs(p0 - p_crit) < 1e-9
    ax.plot(m.T, m.p, color="tab:red" if crit else col,
            lw=2.2 if crit else 1.1, zorder=5 if crit else 2)
ax.set(xlabel="T (K)", ylabel="p (atm)", xlim=(Tw, 690), ylim=(0, 0.022),
       title=r"Trajectories in the $p$–$T$ plane; red is critical")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
plt.show()'''))

cells.append(md("""## Validation

Four checks, in increasing order of how much of the page they exercise."""))

cells.append(code('''# 1. Does the parameter reading survive contact with the paper's own derived
#    quantities? Each of these depends on the whole table at once, and none was
#    used to obtain it, so a single mis-set decimal point breaks all of them.
print("1. Parameter reading, against derived values printed in the paper")
checks = [("A", A, None), ("B", B, None), ("C", C_BASE, None),
          ("ln K (Eq. 22)", np.log(A / C_BASE * np.exp(B_EXP - 20.0)), -2.055),
          ("t_w (Eq. 21b)", A_EXP / Tw, 21.818),
          ("T_M (Eq. 8)", T_M(Tw), 656.6)]
for name, val, paper in checks:
    tail = "" if paper is None else f"   paper {paper:>8}   " \\
                                    f"rel {abs(val - paper) / abs(paper) * 100:5.2f} %"
    print(f"   {name:<16} {val:>14.6g}{tail}")

# and the identity that makes Eqs. 27 and 29 agree at T_M
resid = T_M(Tw) ** 2 - A_EXP * (T_M(Tw) - Tw)
print(f"   identity T_M^2 - a(T_M - T_w) = {resid:.3e}  (exact, from Eq. 8)")'''))

cells.append(code('''# 2. Grid independence. First-order upwind, so successive differences halve.
print("2. Grid independence: hot spot at Tw = 625 K, p0 = 0.0164 atm")
prev, vals = None, []
for n in (200, 400, 800, 1600, 3200):
    v = hot_spot(0.0164, Tw, n_z=n)
    vals.append(v)
    print(f"   n_z = {n:5d}   T_max = {v:.5f} K"
          + ("" if prev is None else f"   change {v - prev:+.2e}"))
    prev = v
richardson = vals[-1] + (vals[-1] - vals[-2])
print(f"   Richardson estimate {richardson:.5f} K; "
      f"error at n_z = 1500 is about {abs(richardson - vals[-2]):.1e} K")'''))

cells.append(code('''# 3. The four worked examples, one by one. `kind` separates the values that
#    test the closed-form criteria from those that test the ODE solve.
rows = []


def rec(example, quantity, computed):
    hit = ex[(ex.example == example) & (ex.quantity == quantity)]
    if hit.empty:
        return
    paper = float(hit.value.iloc[0])
    rows.append({"example": example, "quantity": quantity, "kind": hit.kind.iloc[0],
                 "paper": paper, "computed": computed,
                 "rel_pct": abs(computed - paper) / abs(paper) * 100})


# Example 1(a): first criterion
TM = T_M(Tw)
lo1, up1, mean1 = inlet_limits(TM, Tw)
rec("1a", "T_cr", TM)
rec("1a", "p0_lower", lo1)
rec("1a", "p0_upper", up1)
rec("1a", "p0_mean", mean1)
cp = critical_points(Tw)
rec("1a", "p0_critical", cp["p0_1"])

# Example 1(b): second criterion, solving Eq. 20 instead of reading Figs. 6-7
K = A / C_BASE * np.exp(B_EXP - 20.0)
Tt = T_cr_tangency(Tw)
rec("1b", "t_w", A_EXP / Tw)
rec("1b", "ln_K", np.log(K))
rec("1b", "delta_t", A_EXP / Tw - 20.0)
rec("1b", "ln_K_r", np.log(K) - (A_EXP / Tw - 20.0))
rec("1b", "t", A_EXP / Tt)
rec("1b", "t_r", A_EXP / Tt - (A_EXP / Tw - 20.0))
rec("1b", "T_cr", Tt)
lo2, up2, mean2 = inlet_limits(Tt, Tw)
rec("1b", "p0_lower", lo2)
rec("1b", "p0_upper", up2)
rec("1b", "p0_mean", mean2)
rec("1b", "p0_critical", cp["p0_2"])

# Example 2: subcritical, an imposed hot spot of 640 K
Tm2 = 640.0
X2 = A / C_BASE * kexp(Tm2)
p_m2 = A / B * (Tm2 - Tw) * (1.0 + 1.0 / np.sqrt(X2) + 1.0 / X2)     # Eq. 30
rec("2", "p0_mean", p_m2)
rec("2", "p0_exact", back_integrate(Tm2, maxima_curve(Tm2, Tw), Tw))
rec("2", "T_max", hot_spot(p_m2, Tw))

# Example 3: given p0 and an imposed hot spot, find the wall temperature
p0_3, Tm3 = 0.0075, 675.0
dT_ad3 = B / A * p0_3                                                # Eq. 32
Q3 = 1.0 / np.sqrt(A / C_BASE * kexp(Tm3))                           # Eq. 34
dT_eff3 = dT_ad3 / (1.0 + Q3 + Q3 ** 2)                              # Eq. 31
Tw3 = Tm3 - dT_eff3                                                  # Eq. 33
rec("3", "dT_ad", dT_ad3)
rec("3", "Q", Q3)
rec("3", "dT_eff", dT_eff3)
rec("3", "T_w", Tw3)
rec("3", "T_max", hot_spot(p0_3, Tw3))

# Example 4: given p0 and Tw, find the radius that is critical
p0_4 = 0.0125
dT_ad4 = B / A * p0_4
dT_eff4 = T_M(Tw) - Tw
Q4 = brentq(lambda q: 1.0 + q + q * q - dT_ad4 / dT_eff4, 0.0, 1e3)
C4 = Q4 ** 2 * A * kexp(T_M(Tw))
rec("4", "T_cr", T_M(Tw))
rec("4", "dT_eff", dT_eff4)
rec("4", "dT_ad", dT_ad4)
rec("4", "Q", Q4)
rec("4", "C", C4)
rec("4", "R", 2.0 * P["U"] / (P["c_p"] * C4))

comp = pd.DataFrame(rows)
print("3. The four worked examples of Section 6")
print(comp.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
for kind, sub in comp.groupby("kind"):
    print(f"   {kind:<11}: n = {len(sub):2d}   mean |dev| = {sub.rel_pct.mean():.3f} %"
          f"   worst = {sub.rel_pct.max():.3f} % ({sub.loc[sub.rel_pct.idxmax(), 'quantity']})")
print(f"   overall    : n = {len(comp):2d}   mean |dev| = {comp.rel_pct.mean():.3f} %"
      f"   worst = {comp.rel_pct.max():.3f} %")'''))

cells.append(code('''# 4. Two independent routes to the same critical inlet pressure. The pymrm
#    reactor bisection discretises the reactor in z; the back-integration is a
#    quadrature in the p-T plane with an adaptive Runge-Kutta and never forms
#    the reactor grid at all. They should agree, and they are the only pair of
#    numbers on this page that does not involve the paper.
print("4. pymrm bisection vs phase-plane back-integration (criterion 1)")
Tws = np.arange(600.0, 701.0, 20.0)
worst_cross = 0.0
bracketed = True
for T in Tws:
    cpts = critical_points(T)
    fwd = critical_p0(T, hi=0.08)
    rel = abs(fwd - cpts["p0_1"]) / cpts["p0_1"] * 100
    worst_cross = max(worst_cross, rel)
    lo, up, _ = inlet_limits(cpts["T1"], T)
    inside = lo <= cpts["p0_1"] <= up
    bracketed &= bool(inside)
    print(f"   Tw = {T:.0f} K:  back-integration {cpts['p0_1']:.5f}   "
          f"pymrm {fwd:.5f}   {rel:5.2f} %   inside its bracket: {inside}")
print(f"   worst disagreement between the two methods: {worst_cross:.2f} %")
print(f"   criterion 1 brackets its own critical value everywhere: {bracketed}")

report_agreement("D2.2", {
    "examples_mean_dev_pct": comp.rel_pct.mean(),
    "examples_worst_dev_pct": comp.rel_pct.max(),
    "criterion_mean_dev_pct": comp[comp.kind == "criterion"].rel_pct.mean(),
    "integration_mean_dev_pct": comp[comp.kind == "integration"].rel_pct.mean(),
    "cross_method_worst_pct": worst_cross,
    "T_M_at_625K": T_M(625.0),
    "p0_critical_at_625K": critical_points(625.0)["p0_1"],
})'''))

cells.append(md(r"""## What pymrm adds

**Their Figs. 6 and 7 exist only because there was no computer.** Both are
graphical solutions of one implicit equation, Eq. 20, and the paper says as much
— the figures are there "to facilitate the solution of (20)". Solving it
directly turns them from *inputs* into *checks*: the curve below is their Fig. 7
recomputed rather than read, and their $t_r = 19.013$ is a value they took off
that graph by eye."""))

cells.append(code('''# Their Fig. 7: ln K_r against t_r, where K_r is K referred to t_w = 20.
tr = np.linspace(19.0, 20.0, 240)


def lnK_at(t, tw):
    return np.log((t - 2.0) / (t * np.exp(-t + 20.0))
                  * (1.0 - t * (1.0 - t / tw)))


fig, ax = plt.subplots(figsize=(6.0, 4.2))
ax.plot(tr, lnK_at(tr, 20.0), "k-", lw=1.8)
ax.plot([19.013], [np.log(K) - (A_EXP / Tw - 20.0)], "o", color="tab:red", ms=8,
        label=f"their Example 1(b):\\n$t_r$ = 19.013 read off the graph")
ax.set(xlabel=r"$t_r$", ylabel=r"$\\ln K_r$", xlim=(19.0, 20.0), ylim=(-5, 0),
       title="Their Fig. 7, recomputed from Eq. 20 rather than read")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
plt.show()

t_direct = A_EXP / T_cr_tangency(Tw)
print(f"t from solving Eq. 20 directly : {t_direct:.4f}   their graphical value 20.831")
print(f"t_r equivalent                 : {t_direct - (A_EXP / Tw - 20.0):.4f}"
      f"   their 19.013")
print(f"(T_i)_t                        : {T_cr_tangency(Tw):.2f} K   their 654.6 K")'''))

cells.append(md(r"""**And their Fig. 8 was four curves through a handful of
points.** Each exact value on it cost a numerical back-integration, so the paper
computed few. A bisection on cold-started solves costs nothing, so the boundary
can be drawn continuously — and the interesting quantity is not the boundary
itself but the *width of the bracket*, which is the price of not integrating."""))

cells.append(code('''Tw_sweep = np.arange(600.0, 701.0, 5.0)
lo1s, up1s, m1s, lo2s, up2s, m2s, ex1, ex2 = ([] for _ in range(8))
for T in Tw_sweep:
    cpts = critical_points(T)
    a1, b1, c1 = inlet_limits(cpts["T1"], T)
    lo1s.append(a1); up1s.append(b1); m1s.append(c1); ex1.append(cpts["p0_1"])
    a2, b2, c2 = inlet_limits(cpts["T2"], T)
    lo2s.append(a2); up2s.append(b2); m2s.append(c2); ex2.append(cpts["p0_2"])
lo1s, up1s, m1s, ex1 = map(np.array, (lo1s, up1s, m1s, ex1))
lo2s, up2s, m2s, ex2 = map(np.array, (lo2s, up2s, m2s, ex2))
# Each criterion has its OWN critical trajectory and therefore its own exact
# inlet pressure; comparing both against criterion 1's would flatter one of them.

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax = axes[0]
ax.fill_between(Tw_sweep, lo1s, up1s, color="tab:blue", alpha=0.16,
                label="criterion 1 bracket")
ax.fill_between(Tw_sweep, lo2s, up2s, color="tab:purple", alpha=0.16,
                label="criterion 2 bracket")
ax.plot(Tw_sweep, m1s, color="tab:blue", lw=1.3, ls="--", label="criterion 1 mean")
ax.plot(Tw_sweep, m2s, color="tab:purple", lw=1.3, ls="--", label="criterion 2 mean")
ax.plot(Tw_sweep, ex1, color="k", lw=2.0, label="criterion 1, exact")
ax.plot(Tw_sweep, ex2, color="k", lw=1.2, ls="-.", label="criterion 2, exact")
ax.set(xlabel=r"$T_w$ (K)", ylabel=r"critical $p^0$ (atm)", xlim=(600, 700),
       ylim=(0, 0.045), title="Their Fig. 8, drawn continuously")
ax.legend(fontsize=8)

ax = axes[1]
ax.plot(Tw_sweep, (up1s - lo1s) / ex1 * 100, color="tab:blue", lw=1.8,
        label="criterion 1, bracket width")
ax.plot(Tw_sweep, (up2s - lo2s) / ex2 * 100, color="tab:purple", lw=1.8,
        label="criterion 2, bracket width")
ax.plot(Tw_sweep, (m1s / ex1 - 1) * 100, color="tab:blue", lw=1.2, ls=":",
        label="criterion 1, error of the mean")
ax.plot(Tw_sweep, (m2s / ex2 - 1) * 100, color="tab:purple", lw=1.2, ls=":",
        label="criterion 2, error of the mean")
ax.axhline(0.0, color="k", lw=0.8)
ax.set(xlabel=r"$T_w$ (K)", ylabel="% of that criterion's exact $p^0$",
       xlim=(600, 700), title="Bracket width, and the error of taking the mean")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

print("bracket width, as a percentage of that criterion's own exact critical p0:")
for lbl, lo_, up_, ex_ in (("criterion 1", lo1s, up1s, ex1),
                           ("criterion 2", lo2s, up2s, ex2)):
    w = (up_ - lo_) / ex_ * 100
    print(f"   {lbl}: {w.min():.1f} % to {w.max():.1f} %")
print("error of the mean, as a percentage of that criterion's own exact value:")
for lbl, m_, ex_ in (("criterion 1", m1s, ex1), ("criterion 2", m2s, ex2)):
    e = (m_ / ex_ - 1) * 100
    print(f"   {lbl}: {e.min():+.2f} % to {e.max():+.2f} %")
print(f"criterion 2 is the more conservative: its exact critical p0 is lower "
      f"than criterion 1's at every wall temperature: {bool(np.all(ex2 <= ex1))}")
print(f"   by {(1 - ex2 / ex1).min() * 100:.1f} % to "
      f"{(1 - ex2 / ex1).max() * 100:.1f} %")'''))

cells.append(md(r"""Three things the original could only assert at a point are
now visible across the range.

**The bracket is wide** — 27 % to 66 % of the value it brackets for the first
criterion, 21 % to 37 % for the second — so the useful output of the method was
never the bracket itself but its midpoint. That is what makes the "rather
curious rule" the authors report worth reporting.

**The second criterion's bracket really is the narrower**, as they claim, and
its critical pressure really is the more conservative — lower than the first
criterion's at every wall temperature, by 0.3 % at 600 K widening to 15 % at
700 K.

**But the first criterion's midpoint is the better estimate**, not the second's:
+0.3 % to +1.2 % against its own exact value, where the second criterion's
midpoint runs +1.7 % to +4.7 % high. The paper compares the two only at
$T_w = 625$ K, where they are indistinguishable, and concludes they "entirely
confirm" each other. Over a range they do not quite: narrower bracket and better
midpoint are different virtues, and here they belong to different criteria.

That is the honest limit of this page. Nothing here contradicts the paper — the
criteria are reproduced, not improved — but it does quantify, across the
operating plane, claims the original supports at a single point."""))

cells.append(md(r"""## Reuse

**The criteria on their own.** `T_M`, `T_cr_tangency` and `inlet_limits` need
only the groups $A$, $B$, $C$ and the Arrhenius constants $a$, $b$ — no reactor
object, no pymrm. For a different chemistry, recompute $A$, $B$, $C$ from the
definitions in the model section and the criteria carry over unchanged, provided
the reaction is still pseudo-first-order and irreversible.

**The reactor.** `RunawayTube.residual` is the whole model. To add axial
dispersion, build a `construct_grad`/`construct_div` diffusion operator and add
it to `jac_const` — the source term does not change. To make the wall
temperature vary along the tube, pass an array where `Tw` is used in
`reaction`. To move to a two-phase (heterogeneous) description, the pellet
equation from [`B1.1`](../B1.1-thiele-weisz-hicks/) couples in as `D1.4` does.

**Watch the units.** $c_p$ is volumetric here. If you take $c_p$ from a data
book in kcal/(kg·°C), multiply by $\rho_g$ before using it in $B$ and $C$.

**Related pages.** [`C2.1`](../C2.1-xu-froment-smr/) (the same plug-flow
operators, real kinetics), [`B1.1`](../B1.1-thiele-weisz-hicks/) (multiplicity,
which this problem does *not* have), `C2.10` (o-xylene, the classic hot spot),
`D2.1`, `D2.3`.

**Cite the source, not this page:** Van Welsenaere, R. J. and Froment, G. F.,
*Parametric sensitivity and runaway in fixed bed catalytic reactors*, Chemical
Engineering Science **25**(10) 1503–1516 (1970),
[doi:10.1016/0009-2509(70)85073-4](https://doi.org/10.1016/0009-2509(70)85073-4)."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py", "mimetype": "text/x-python", "name": "python",
        "nbconvert_exporter": "python", "pygments_lexer": "ipython3",
        "version": "3.13.5"},
}
out = Path(__file__).with_name("index.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
