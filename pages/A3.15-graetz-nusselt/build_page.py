#!/usr/bin/env python3
"""Generate index.ipynb for page A3.15. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "The Graetz–Nusselt problem: the constants Graetz could not compute"
description: "The developing thermal boundary layer in a tube, solved as a 2-D PDE — and the two amplitudes Graetz declared out of reach and calibrated out of his own apparatus instead."
categories: [sec:A, struct:S6, tier:T0, data:tier1, phase:liquid]
date: 2026-08-05
---

# The Graetz–Nusselt problem: the constants Graetz could not compute

**Catalog ID:** `A3.15` · **Structure:** `S6` (2-D PDE) · **Tier:** T0

Leo Graetz did not set out to found a branch of convective heat transfer. He was
trying to measure the thermal conductivity of liquids, and the problem now named
after him is the calculation he needed in order to read his own instrument.

The calculation stopped one step short. He obtained the eigenvalues, and then
wrote that the *amplitudes* multiplying them would need "an exhaustive
investigation of the function $V(r,\beta,R)$, which is probably not simple" —
so for physical purposes he determined them **experimentally, as constants of the
apparatus**.

This page computes them, and then asks what that does to his measurements. The
answer is not the one the framing suggests."""))

cells.append(md(r"""## Background

*Everything quoted from the paper below is in the original German with a
translation. Nothing load-bearing rests on a translation alone.*

Graetz's apparatus (his §2) is a flow calorimeter. Liquid at temperature $T_1$
flows under a constant head through a narrow brass tube whose wall is held at
the temperature $T_0$ of running cooling water; the outflow temperature $U$ and
the weight $G$ collected in one minute are measured. From $T_0$, $T_1$, $U$, $G$
and the tube length $l$ the conductivity follows — *if* you can solve for the
temperature field inside the tube.

That is the Graetz problem. He solves it twice.

**First with a flat velocity profile** (§3), which gives Bessel functions and the
classical $J_0$ zeros. He then says plainly that this will not do:

> "Die Annahme jedoch, dass die Flüssigkeit in jedem Punkte der Röhre dieselbe
> mittlere Geschwindigkeit $\alpha$ habe, ist weder in der Natur erfüllt, noch
> wie die Berechnung zeigt, zur Vereinfachung erlaubt."

*(The assumption that the liquid has the same mean velocity $\alpha$ at every
point of the tube is neither fulfilled in nature nor, as the calculation shows,
permissible even as a simplification.)*

**Then with the Poiseuille profile**, which is the problem that carries his name.
The radial equation is no longer Bessel's, and the eigenfunctions are — as he
notes — an extension of the Bessel functions that he expects "may perhaps possess
independent mathematical interest".

Two things make this worth a page rather than a footnote:

1. It is the ancestor of every entrance-length correlation in heat and mass
   transfer, and the asymptotic Nusselt number of a tube at constant wall
   temperature is the first eigenvalue of this problem, halved.
2. **The paper contains measurements** — six liquids, three pressure heads — and
   it contains a calibration step whose weakness is invisible from inside the
   paper and obvious from outside it."""))

cells.append(md(r"""## The published model

**The equation** (journal page 87, read from the page image):

$$
2\alpha\left(1 - \frac{r^{2}}{R^{2}}\right)\frac{\partial u}{\partial z}
= a^{2}\left(\frac{\partial^{2}u}{\partial r^{2}} + \frac{1}{r}\frac{\partial u}{\partial r}\right),
\qquad a^{2} = \frac{k}{\varrho C},
$$

with $u$ measured from the wall temperature, $u = T_1$ at $z = 0$ and $u = 0$ at
$r = R$. $\alpha$ is the **mean** velocity, so $2\alpha(1-r^2/R^2)$ is Poiseuille
flow. Axial conduction is present in his plug-flow equation on page 83 and
dropped here; §7 of the validation checks that.

**The separation** $u = A\,e^{-(a^{2}/2\alpha)\beta^{2}z}V$ gives

$$
\frac{\mathrm{d}^{2}V}{\mathrm{d}r^{2}} + \frac{1}{r}\frac{\mathrm{d}V}{\mathrm{d}r}
+ \beta^{2}\left(1 - \frac{r^{2}}{R^{2}}\right)V = 0,
\qquad V(r,\beta,\infty) = J_{0}(\beta r),
$$

and the eigenvalues are the roots $\mu = \beta R$ of $V(R,\beta,R) = 0$. Graetz
expands $V$ in powers of $x = \beta r$ and prints the resulting polynomial in
$\mu$ and its first two roots.

**The measurable** (page 88) is the mean outlet temperature

$$
U = T_{1}\sum_{i} p_{i}\,
    \exp\!\left(-\frac{a^{2}\mu_{i}^{2}\pi l \varrho}{2G}\right),
$$

and since $a^{2} = k/(\varrho C)$ the density cancels, so **the observations give
$k/c$ directly** — which is what his tables report. The tube radius cancels too,
a point he makes explicitly on page 82.

**And then the step this page is about** (page 88):

> "Darin sind die $p_i$ Constanten, die ausser von den $\mu_i$ nur von $R$
> abhängen. Es wäre natürlich von wesentlichem Vortheil, wenn sich die $p_i$
> bestimmt durch $R$ und $\mu_i$ ausdrücken liessen. Indess gehört dazu eine
> eingehende Untersuchung der Function $V(r,\beta,R)$, die wahrscheinlich nicht
> einfach ist, die aber als Erweiterung der Bessel'schen Function vielleicht
> selbstständiges mathematisches Interesse besitzt. Für physikalische Zwecke
> lassen sich dagegen die $p_i$ sehr einfach als Constanten des Apparates
> experimentell bestimmen."

*(Therein the $p_i$ are constants which, besides on the $\mu_i$, depend only on
$R$. It would of course be of essential advantage if the $p_i$ could be expressed
definitely through $R$ and $\mu_i$. However, that requires an exhaustive
investigation of the function $V(r,\beta,R)$, which is probably not simple, but
which as an extension of the Bessel function may perhaps possess independent
mathematical interest. For physical purposes, on the other hand, the $p_i$ can
very simply be determined experimentally as constants of the apparatus.)*

He then keeps two terms, and fits $p_1$, $p_2$ **and** $k/c$ to three runs of one
liquid at three pressure heads. Three constants, three observations."""))

cells.append(md(r"""## Parameters and assumptions

**His:** steady state; fully developed Poiseuille flow from $z = 0$; constant
wall temperature; constant $k$, $\varrho$, $c$ over the run; no axial conduction
in the Poiseuille case; the measured outflow temperature equals the cross-section
mean of $u$ at $z = l$; and — for the working formula — two terms of the series
suffice, "Das zweite Glied macht noch etwa 3 Proc. in dem Werthe von $k$ aus"
(*the second term still accounts for about 3 per cent in the value of $k$*).

**Ours:** the same, written in dimensionless form. With $s = r/R$ and

$$
x \;=\; \frac{a^{2}z}{2\alpha R^{2}} \;=\; \frac{a^{2}\pi z}{2\dot V}
\;=\; \frac{(k/c)\,\pi z}{2G},
$$

the whole problem loses every parameter:

$$
(1-s^{2})\frac{\partial\theta}{\partial x}
= \frac{1}{s}\frac{\partial}{\partial s}\!\left(s\frac{\partial\theta}{\partial s}\right)
+ \frac{1}{\mathrm{Pe}^{2}}\frac{\partial^{2}\theta}{\partial x^{2}},
\qquad \theta(0,s) = 1,\quad \theta(x,1) = 0,
$$

with $\theta = (T-T_0)/(T_1-T_0)$ and $\mathrm{Pe} = 2\alpha R/a^{2}$ the axial
Péclet number. Setting $1/\mathrm{Pe}^{2} = 0$ is Graetz's equation exactly.
$x$ is twice the Graetz coordinate $x^{*} = z/(D_{h}\mathrm{Re}\,\mathrm{Pr})$
used in the modern literature.

The **cup-mean** (flow-weighted) outlet temperature is what a collected sample
measures:

$$
\theta_{\text{cup}}(x) = 4\!\int_{0}^{1}(1-s^{2})\,\theta\,s\,\mathrm{d}s
= \sum_{i} p_{i}\,e^{-\mu_{i}^{2}x}.
$$

Note that Graetz's own page-88 formula integrates $V$ against $r\,\mathrm{d}r$
and divides by $R^{2}\pi$ — that is the **area** mean, not the cup mean. §5
measures what the difference is worth."""))

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
from fractions import Fraction
from math import gamma

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, fsolve
from scipy.sparse.linalg import splu
from scipy.special import jn_zeros
from pymrm import (construct_grad, construct_div, construct_convflux_upwind,
                   interp_cntr_to_stagg_tvd, vanleer)
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A3.15-graetz-nusselt"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

cells.append(md(r"""## The data

Two files. `graetz-1882-observations.csv` is all 51 runs of his §4: six liquids,
the three temperatures, the collected weight, and the $k/c$ and $k$ he computes
from them. `graetz-1882-printed-constants.csv` is every constant printed in
§§2–4 with the journal page it is on, including the series coefficients and the
two fitted amplitudes.

**These are measurements**, which makes this one of the few tier-1 pages in the
gallery — but the $k/c$ column is *derived*, by his two-term formula with his
two fitted amplitudes. The raw measurement is the four temperature-and-mass
columns, and this page inverts those independently.

**Read the sidecars before using either file.** They record: three rows flagged
because the fit that uses them has zero degrees of freedom; the internal
consistency checks (`T = (T₁+U)/2`, `k = c·(k/c)`) and the one liquid that fails
the second; and three printed defects in the constants file."""))

cells.append(code('''obs = load_data("graetz-1882-observations.csv", page=PAGE)
con = load_data("graetz-1882-printed-constants.csv", page=PAGE)
obs_meta = load_meta("graetz-1882-observations.csv", page=PAGE)
con_meta = load_meta("graetz-1882-printed-constants.csv", page=PAGE)
C = dict(zip(con.quantity, con.value))
# journal page each constant is printed on, as an integer -- so that statements
# like "N pages earlier" are computed from the file rather than typed
PAGE_OF = {q: int(str(s).split()[-1]) for q, s in zip(con.quantity, con.source)}

obs["theta"] = (obs.T_out_C - obs.T_bath_C) / (obs.T_in_C - obs.T_bath_C)
CAL = obs[obs.calibration_row == 1].reset_index(drop=True)

print(f"{len(obs)} runs, {obs.liquid.nunique()} liquids, "
      f"{len(con)} printed constants")
print(obs.groupby("liquid", sort=False)
         .agg(runs=("theta", "size"), G_min=("mass_flow_g_min", "min"),
              G_max=("mass_flow_g_min", "max"), theta_min=("theta", "min"),
              theta_max=("theta", "max")).to_string())
print(f"\\n{cite_data(obs_meta)}")

# the two integrity checks the tables pay for -- both stated in the sidecar
print("\\nintegrity check 1:  T = (T_in + T_out)/2, which the paper states on p. 90")
dT = 0.5 * (obs.T_in_C + obs.T_out_C) - obs.T_mean_C_printed
print(f"   exact in {int((abs(dT) <= 0.02).sum())} of {len(obs)} rows; "
      f"worst residual {dT.abs().max():.2f} K "
      f"({int((obs.calibration_row[abs(dT) > 0.02] == 1).sum())} of the "
      "offenders are the calibration rows, which are means of five observations)")
print("\\nintegrity check 2:  k = c * (k/c)")
for liq, g in obs.groupby("liquid", sort=False):
    imp = g.k_printed / g.k_over_c_printed
    stated = g.specific_heat.iloc[0]
    flag = "" if (pd.isna(stated) or abs(imp.mean() / stated - 1) < 0.01) else "   <-- 4.1 % gap"
    print(f"   {liq:11s} stated c = {stated!s:>5}   implied {imp.mean():.4f} "
          f"+/- {imp.std():.4f}{flag}")
print("   turpentine is a printed slip in one of the two numbers, not scatter;")
print("   the page works in k/c throughout, so nothing here depends on it.")

# the three flagged rows, up front
print("\\nthe three rows Graetz fitted p1, p2 and k/c to (his stars, p. 91):")
print(CAL[["pressure_head", "T_in_C", "T_out_C", "T_bath_C", "mass_flow_g_min",
           "theta", "k_over_c_printed"]].to_string(index=False))
print(f"   flow range {CAL.mass_flow_g_min.max()/CAL.mass_flow_g_min.min():.2f}x, "
      "three constants, three observations -> zero degrees of freedom")'''))

cells.append(md(r"""## PyMRM implementation

The 2-D solve is the whole page: it produces the eigenvalue *and* the amplitude
without ever forming an eigenfunction, so it is an independent witness against
the series route in §3 of the validation.

Three things to note in the code.

- **`construct_div(..., nu=1)` on the radial axis** — cylindrical geometry, which
  is what puts the $\tfrac{1}{s}\partial_s(s\,\partial_s)$ in the equation. The
  axial axis is `nu=0`, Cartesian.
- **Every `bc` carries its physical equation in a comment**, because `a` is
  written on the *outward* normal and therefore means different things at the two
  ends of an axis.
- **First-order upwind is not good enough here.** The exact problem has no axial
  diffusion at all, so upwind numerical diffusion is pure error and it converges
  at first order — §3 and the break table measure it. A van Leer TVD deferred correction on the
  axial convection recovers second order at a few dozen extra back-substitutions
  through one stored factorisation. It needs under-relaxation, and the cell below
  measures that rather than asserting it: without it the limiter switches, the
  iteration limit-cycles, and every number on this page would depend on the
  iteration count."""))

cells.append(code('''def solve_graetz(X=1.0, n_z=600, n_r=128, nu_r=1, inv_pe2=0.0, plug=False,
                 z_grade=3.0, r_grade=2.0, tvd=True, tol=1e-11, maxit=200,
                 omega=0.6):
    """theta(x, s) for the Graetz problem on x in [0, X], s in [0, 1].

        (1-s^2) dtheta/dx = (1/s) d/ds(s dtheta/ds) + inv_pe2 * d2theta/dx2

    Graded faces: x = X t^z_grade puts cells where the boundary layer is thin,
    s = 1-(1-t)^r_grade puts them at the wall. `nu_r`, `plug` and `inv_pe2` exist
    so the validation section can break or switch them.
    """
    t = np.linspace(0.0, 1.0, n_z + 1)
    z_f = X * t ** z_grade
    z_c = 0.5 * (z_f[:-1] + z_f[1:])
    t = np.linspace(0.0, 1.0, n_r + 1)
    r_f = 1.0 - (1.0 - t) ** r_grade
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n_z, n_r)

    # u(r)/(2 * mean velocity); the 2 is absorbed into the definition of x
    w = np.ones_like(r_c) if plug else (1.0 - r_c ** 2)
    v = np.broadcast_to(w, (n_z + 1, n_r))

    # a dtheta/dn + b theta = d, n the OUTWARD normal.
    #   inlet  x=0 : theta = 1              -> a=0, b=1, d=1   (Dirichlet)
    #   outlet x=X : dtheta/dx = 0          -> a=1, b=0, d=0   (n = +x there)
    bc_z = ({"a": 0.0, "b": 1.0, "d": 1.0},
            {"a": 1.0, "b": 0.0, "d": 0.0})
    #   axis  s=0 : dtheta/ds = 0, symmetry -> a=1, b=0, d=0   (n = -s there,
    #                                          so a=1,d=0 is zero flux either way)
    #   wall  s=1 : theta = 0, wall at T_0  -> a=0, b=1, d=0   (Dirichlet)
    bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0},
            {"a": 0.0, "b": 1.0, "d": 0.0})

    conv, conv_bc = construct_convflux_upwind(shape, z_f, z_c, bc_z, v=v, axis=0)
    gz, gz_bc = construct_grad(shape, z_f, z_c, bc_z, axis=0)
    dz = construct_div(shape, z_f, nu=0, axis=0)        # nu=0: Cartesian, axial
    gr, gr_bc = construct_grad(shape, r_f, r_c, bc_r, axis=1)
    dr = construct_div(shape, r_f, nu=nu_r, axis=1)     # nu=1: cylindrical, radial

    # constant operators, assembled and factorised once
    A = (dz @ (conv - inv_pe2 * gz) - dr @ gr).tocsc()
    b = np.asarray((dz @ (conv_bc - inv_pe2 * gz_bc) - dr @ gr_bc).todense()).ravel()
    lu = splu(A)
    theta = lu.solve(-b).reshape(shape)

    nit, incr = 0, 0.0
    if tvd:                       # van Leer deferred correction, under-relaxed
        for nit in range(1, maxit + 1):
            _, dth = interp_cntr_to_stagg_tvd(theta, z_f, z_c, bc_z, v, vanleer,
                                              axis=0)
            corr = np.asarray(dz @ (v * dth).reshape(-1, 1)).ravel()
            new = theta + omega * (lu.solve(-b - corr).reshape(shape) - theta)
            incr = np.max(np.abs(new - theta)) / omega
            theta = new
            if incr < tol:
                break
        else:
            raise RuntimeError(f"deferred correction did not converge: {incr:.2e}")

    area_w = 2.0 * r_c * np.diff(r_f)                  # area weights, sum 1
    cup_w = area_w * w / np.sum(area_w * w)            # flow weights, sum 1
    return dict(x=z_c, x_f=z_f, s=r_c, s_f=r_f, theta=theta, nit=nit, incr=incr,
                cup=theta @ cup_w, area=theta @ area_w)


def decay(sol, xmin=0.35, key="cup"):
    """(mu_1, p_1) from the far-field slope and intercept of ln theta.

    `key` selects which cross-section mean is fitted: "cup" (flow-weighted, what
    a collected sample measures) or "area" (Graetz's own p. 88 formula).
    """
    x, cup = sol["x"], sol[key]
    m = x > xmin
    slope, icept = np.polyfit(x[m], np.log(cup[m]), 1)
    return np.sqrt(-slope), np.exp(icept)


def nusselt(sol):
    """Local Nu = -2 (dtheta/ds)|_{s=1} / theta_cup, from the wall gradient."""
    th, s_c, s_f = sol["theta"], sol["s"], sol["s_f"]
    dth_ds = (0.0 - th[:, -1]) / (s_f[-1] - s_c[-1])   # theta = 0 at the wall
    return sol["x"], -2.0 * dth_ds / sol["cup"]


# The under-relaxation is a claim about the solver, so measure it rather than
# assert it. omega = 1 is the textbook deferred correction.
for om in (1.0, 0.6):
    try:
        s_ = solve_graetz(n_z=300, n_r=64, omega=om, maxit=120)
        print(f"omega = {om:.1f}: converged in {s_['nit']:3d} iterations to "
              f"{s_['incr']:.1e}")
    except RuntimeError as e:
        print(f"omega = {om:.1f}: NOT converged in 120 iterations, stalled at "
              f"an increment of {str(e).split()[-1]}")'''))

cells.append(md(r"""And the series route, which exists only so that the 2-D solve
has something independent to be checked against. `eigen` shoots the radial ODE;
`amplitudes` projects the uniform inlet profile onto the eigenfunctions with the
weight $(1-s^{2})s$ that makes them orthogonal.

The cup and area amplitudes come out of the same integral with different
weights, which is exactly the difference between what Graetz wrote down and what
his thermometer measured."""))

cells.append(code('''def V_profile(mu, n=40001):
    """V(s) and V'(s) for V'' + V'/s + mu^2 (1-s^2) V = 0, V(0)=1, V'(0)=0."""
    def rhs(t, y):
        return [y[1], -(y[1] / t if t > 0 else 0.0) - mu ** 2 * (1 - t ** 2) * y[0]]
    s0 = 1e-10
    sol = solve_ivp(rhs, [s0, 1.0], [1 - mu ** 2 * s0 ** 2 / 4, -mu ** 2 * s0 / 2],
                    t_eval=np.linspace(s0, 1.0, n), rtol=1e-13, atol=1e-15)
    return sol.t, sol.y[0], sol.y[1]


def V_end(mu, rtol=1e-12):
    """V(1) alone: no t_eval, so the integrator is not asked to sample."""
    def rhs(t, y):
        return [y[1], -(y[1] / t if t > 0 else 0.0) - mu ** 2 * (1 - t ** 2) * y[0]]
    s0 = 1e-10
    return solve_ivp(rhs, [s0, 1.0],
                     [1 - mu ** 2 * s0 ** 2 / 4, -mu ** 2 * s0 / 2],
                     rtol=rtol, atol=1e-15).y[0, -1]


def eigen(n_modes=10, hi=40.0):
    """Roots of V(1; mu) = 0. Coarse cheap scan to bracket, accurate brentq."""
    grid = np.linspace(0.5, hi, 300)          # roots are ~4 apart; spacing <= 0.24
    val = np.array([V_end(m, rtol=1e-9) for m in grid])
    roots = [brentq(V_end, grid[i], grid[i + 1], xtol=1e-11)
             for i in np.nonzero(val[:-1] * val[1:] < 0)[0]]
    return np.array(roots[:n_modes])


def amplitudes(mus, n=40001):
    """Cup-mean and area-mean expansion coefficients for a uniform inlet."""
    p_cup, p_area, ident = [], [], []
    for mu in mus:
        s, V, dV = V_profile(mu, n)
        N = np.trapezoid((1 - s ** 2) * V * s, s)     # = -V'(1)/mu^2, by the ODE
        D = np.trapezoid((1 - s ** 2) * V ** 2 * s, s)
        M = np.trapezoid(V * s, s)
        p_cup.append(4 * N * N / D)
        p_area.append(2 * (N / D) * M)
        ident.append(abs(N / (-dV[-1] / mu ** 2) - 1))
    return np.array(p_cup), np.array(p_area), np.array(ident)


MU = eigen()
P_CUP, P_AREA, IDENT = amplitudes(MU)
SIGMA = jn_zeros(0, 10)          # the plug-flow case, Graetz's p. 83
P_PLUG = 4.0 / SIGMA ** 2        # and his own amplitudes for it, p. 84

print(f"{'i':>2}{'mu_i':>12}{'p_i (cup)':>12}{'p_i (area)':>12}"
      f"{'sigma_i (plug)':>16}{'4/sigma_i^2':>13}")
for i in range(6):
    print(f"{i+1:2d}{MU[i]:12.6f}{P_CUP[i]:12.6f}{P_AREA[i]:12.6f}"
          f"{SIGMA[i]:16.6f}{P_PLUG[i]:13.6f}")
_sl = np.polyfit(np.log(np.arange(1, 11)), np.log(P_CUP), 1)[0]
print(f"\\nsum of the first 10: cup {P_CUP.sum():.4f}, area {P_AREA.sum():.4f} "
      f"(both -> 1 as more modes")
print(f"are added, but only as p_i ~ i^{_sl:.2f}, so a truncated series is a real "
      "approximation)")
print(f"structural identity  N_i = -V'(1)/mu_i^2 holds to {IDENT.max():.1e} -- it "
      "follows from")
print("the ODE itself, so it checks the quadrature, not the physics.")'''))

cells.append(md(r"""## Results

Three views of the same field: the developing profile, the cup-mean decay that is
what Graetz measures, and the local Nusselt number that connects the two ends of
the tube to their classical limits."""))

cells.append(code('''sol = solve_graetz(X=1.0, n_z=600, n_r=128)
MU1_2D, P1_2D = decay(sol)
x_nu, Nu = nusselt(sol)
C_LEVEQUE = 2.0 / (gamma(4.0 / 3.0) * 9.0 ** (1.0 / 3.0))

fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))

xs_show = [1e-4, 1e-3, 1e-2, 0.05, 0.2]
for xt, col in zip(xs_show, plt.cm.viridis(np.linspace(0.1, 0.85, len(xs_show)))):
    j = np.argmin(abs(sol["x"] - xt))
    axes[0].plot(np.r_[sol["theta"][j], 0.0], np.r_[sol["s"], 1.0], color=col,
                 lw=1.7, label=f"$x$ = {sol['x'][j]:.0e}".replace("e-0", "e-"))
axes[0].set(xlabel=r"$\\theta$", ylabel=r"$s = r/R$",
            title="the boundary layer eating inwards")
axes[0].legend(fontsize=7.5, loc="lower left")

one = P_CUP[0] * np.exp(-MU[0] ** 2 * sol["x"])
ser = np.array([np.sum(P_CUP * np.exp(-MU ** 2 * xx)) for xx in sol["x"]])
ser2 = np.array([np.sum(P_CUP[:2] * np.exp(-MU[:2] ** 2 * xx)) for xx in sol["x"]])
axes[1].semilogx(sol["x"], sol["cup"] / one, color="tab:blue", lw=2.0,
                 label="2-D solve")
axes[1].semilogx(sol["x"], ser / one, "r:", lw=1.6, label="10-term series")
axes[1].semilogx(sol["x"], ser2 / one, "g-.", lw=1.2, label="2 terms")
axes[1].axhline(1.0, color="k", ls="--", lw=1.0)
gm = ((obs.k_over_c_printed * np.pi * obs.tube_length_cm)
      / (2 * obs.mass_flow_g_min)).values
axes[1].axvspan(gm.min(), gm.max(), color="tab:orange", alpha=0.16)
axes[1].text(gm.min() * 1.05, 1.45, "his runs", fontsize=8, color="tab:orange")
axes[1].set(xlabel=r"$x = (k/c)\\pi z/2G$", xlim=(3e-4, 1.0), ylim=(0.95, 1.7),
            ylabel=r"$\\theta_{\\rm cup}\\,/\\,p_1e^{-\\mu_1^2x}$",
            title="how many terms the thermometer needs")
axes[1].legend(fontsize=8)

m = x_nu > 3e-5
axes[2].loglog(x_nu[m], Nu[m], color="tab:blue", lw=2.0, label="2-D solve")
axes[2].loglog(x_nu[m], C_LEVEQUE * (x_nu[m] / 2) ** (-1 / 3), "k--", lw=1.1,
               label=r"Leveque $\\propto x^{-1/3}$")
axes[2].axhline(MU[0] ** 2 / 2, color="tab:red", lw=1.2, ls=":")
axes[2].text(2e-2, MU[0] ** 2 / 2 * 1.09, r"$\\mu_1^2/2$", color="tab:red", fontsize=9)
axes[2].set(xlabel=r"$x$", ylabel=r"local $\\mathrm{Nu}$", ylim=(2, 300),
            title="both ends have a closed form")
axes[2].legend(fontsize=8)
fig.tight_layout()
plt.show()

print(f"2-D solve, n_z = 600, n_r = 128, {sol['nit']} deferred-correction "
      f"iterations to {sol['incr']:.1e}")
print(f"  mu_1 from the far-field decay   {MU1_2D:.6f}")
print(f"  p_1  from the far-field decay   {P1_2D:.6f}")
print(f"  asymptotic Nu = mu_1^2/2        {MU1_2D**2/2:.4f}")
print(f"  Graetz's printed mu_1           {C['graetz_root_mu1']}")
print(f"  Graetz's fitted p_1             {C['p1_working_formula']}   <-- the "
      "constant he could not compute")'''))

cells.append(md("""## Validation

Seven checks and a defect-injection table. Two are against constants Graetz
prints, two against limits he never had, one is the arithmetic of his own
calibration re-run, one is the claim he makes about axial conduction, and one is
grid refinement on both axes and in the inlet region separately — followed by
§3b, which measures what the resulting agreement number is *not*."""))

cells.append(code('''# 1. His printed power series, against exact rational arithmetic.
#    V = sum a_k x^2k with a_k = sum_j c[k][j] mu^-2j; the recurrence follows from
#    the ODE. Setting x = mu collapses it to a polynomial in mu^2, which is what
#    he prints. Nothing here touches the pymrm solve.
def series_coeffs(N=60):
    c = [[Fraction(1)]]
    for m in range(N):
        row = [Fraction(0)] * (m + 2)
        for j in range(m + 2):
            t1 = c[m][j] if j < len(c[m]) else Fraction(0)
            t2 = c[m - 1][j - 1] if (m >= 1 and j >= 1 and j - 1 < len(c[m - 1])) \\
                else Fraction(0)
            row[j] = -(t1 - t2) / Fraction((2 * m + 2) ** 2)
        c.append(row)
    return [float(sum(c[k][k - p] for k in range(p, N + 1) if k - p < len(c[k])))
            for p in range(6)]


exact_poly = series_coeffs()
printed_poly = [1.0, -C["series_coeff_mu2"], C["series_coeff_mu4"],
                -C["series_coeff_mu6"], C["series_coeff_mu8"],
                -C["series_coeff_mu10"]]
print("1. The polynomial Graetz prints for V(mu) = 0 (journal p. 87)")
print(f"   {'power':>7}{'printed':>18}{'exact (rationals)':>22}{'rel. error':>13}")
poly_err = []
for k in range(6):
    e = abs(printed_poly[k] / exact_poly[k] - 1) if exact_poly[k] else 0.0
    poly_err.append(e)
    flag = "   <-- wrong" if e > 1e-4 else ""
    print(f"   mu^{2*k:<4d}{printed_poly[k]:18.10g}{exact_poly[k]:22.10g}"
          f"{e:13.2e}{flag}")
print("   Four coefficients are exact to every digit he prints, which is a strong")
print("   check on the page-image reading -- the text layer gives 0.1575 for the")
print("   0.1875 above. The last two are arithmetic slips of his: +0.58 % and a")
print(f"   factor of {1/(printed_poly[5]/exact_poly[5]):.2f}.")

# The one glyph on the page that had to be argued rather than simply read. The
# factor-of-2.65 claim turns on it, so state what the alternative would give.
MU10_ALT = -9.4938e-9            # the same digits with the broken group read as 9
MU10_ALT_ERR = abs(MU10_ALT / exact_poly[5] - 1)
print("\\n   ONE CAVEAT ON THE READING, because a headline number turns on it.")
print("   The third digit-group of the mu^10 coefficient is a BROKEN GLYPH on the")
print("   300 dpi bitmap. Magnified 16x it has an open left side, two left-hand")
print("   terminals and no closed upper bowl, so it is a 3 and the coefficient is")
print(f"   {printed_poly[5]:.4e}, which is the factor of "
      f"{1/(printed_poly[5]/exact_poly[5]):.2f} above. Had it been a 9:")
print(f"     read as 3 (taken here) {printed_poly[5]:12.4e}   "
      f"rel err {poly_err[5]:.3f}   <-- factor {1/(printed_poly[5]/exact_poly[5]):.2f}")
print(f"     read as 9 (rejected)   {MU10_ALT:12.4e}   "
      f"rel err {MU10_ALT_ERR:.3f}   <-- 2.4 %, an ordinary rounding slip")
print("   So the mu^8 slip stands whatever this glyph is, but the mu^10 FACTOR")
print("   rests on one damaged character read at 16x. The sidecar records it.")'''))

cells.append(code('''# 2. The eigenvalues. His mu_1 is excellent and his mu_2 is not, and the reason
#    is visible in his own polynomial.
print("2. Eigenvalues")
roots_printed = np.sort(np.roots(printed_poly[::-1]))
roots_printed = np.sqrt(roots_printed[(abs(roots_printed.imag) < 1e-9)
                                      & (roots_printed.real > 0)].real)
roots_exact_t = np.sort(np.roots(exact_poly[::-1]))
roots_exact_t = np.sqrt(roots_exact_t[(abs(roots_exact_t.imag) < 1e-9)
                                      & (roots_exact_t.real > 0)].real)
mu1_err = abs(C["graetz_root_mu1"] / MU[0] - 1)
mu2_err = abs(C["graetz_root_mu2"] / MU[1] - 1)
print(f"   mu_1  printed {C['graetz_root_mu1']:<9.4f} exact {MU[0]:.6f}"
      f"   rel {mu1_err:.2e}")
print(f"   mu_2  printed {C['graetz_root_mu2']:<9.2f} exact {MU[1]:.6f}"
      f"   rel {mu2_err:.2e}   <-- 2.7 % low")
print(f"   his own 10th-order polynomial has positive roots at "
      f"{np.array2string(roots_printed, precision=4)},")
print(f"   and with the coefficients CORRECTED it has "
      f"{np.array2string(roots_exact_t, precision=4)}. So the series can place")
print("   mu_1 to four figures and cannot place mu_2 at all: the 6.50 he prints")
print("   comes from neither polynomial, and truncation, not arithmetic, is why.")
print(f"\\n   plug-flow zeros he prints (p. 83), against J_0:")
for i, key in enumerate(["bessel_zero_1", "bessel_zero_2", "bessel_zero_3",
                         "bessel_zero_4"]):
    print(f"     sigma_{i+1}  printed {C[key]:<12.6f} exact {SIGMA[i]:.6f}"
          f"   rel {abs(C[key]/SIGMA[i]-1):.2e}")'''))

cells.append(code('''# 3. The 2-D pymrm solve against the series route, refined on each axis
#    separately. These share NO code: one is a finite-volume solve of the PDE,
#    the other shoots the radial ODE and projects. They share only the physics.
print("3. 2-D solve vs the eigenfunction route, refining one axis at a time")
LAD = {}
for tag, fixed, sweep in (("n_z", dict(n_r=256), (75, 150, 300, 600)),
                          ("n_r", dict(n_z=800), (16, 32, 64, 128))):
    print(f"   refine {tag} ({', '.join(f'{k} = {v}' for k, v in fixed.items())}):")
    print(f"     {tag:>6}{'mu_1':>12}{'err':>11}{'order':>8}"
          f"{'p_1':>12}{'err':>11}{'order':>8}")
    prev, errs = None, []
    for n in sweep:
        s_ = solve_graetz(**fixed, **{tag: n})
        mu_, p_ = decay(s_)
        e_mu, e_p = abs(mu_ / MU[0] - 1), abs(p_ / P_CUP[0] - 1)
        errs.append((e_mu, e_p))
        o = ("", "") if prev is None else (f"{np.log2(prev[0]/e_mu):8.2f}",
                                           f"{np.log2(prev[1]/e_p):8.2f}")
        print(f"     {n:6d}{mu_:12.6f}{e_mu:11.2e}{o[0]:>8}"
              f"{p_:12.6f}{e_p:11.2e}{o[1]:>8}")
        prev = (e_mu, e_p)
    # over the first 4x of the ladder, where the OTHER axis is not yet the floor
    LAD[tag] = 0.5 * np.log2(errs[0][0] / errs[2][0])
    print(f"     observed order in mu_1 over the first 4x: {LAD[tag]:.2f}")
print("   Both axes are second order, which is the deferred correction working:")
print("   plain first-order upwind gives order 1 and needs 16x the cells for the")
print("   same error (see the break table).")
print("   The order is read over the FIRST 4x of each ladder on purpose. Refining")
print("   one axis cannot take the error below what the other axis contributes, and")
print("   the last rung of each ladder is already at that floor -- which is why the")
print("   per-doubling column ends at 3.17 and why p_1 in the n_r ladder stops")
print("   moving at 1.5e-4 and then drifts. A ratio read there is not an order.")'''))

cells.append(md(r"""**3b. What the published `mu1_rel_err_2d` is, and what it is not.**

The orders above are the load-bearing result. The single agreement number at the
published grid is *not* an accuracy: the axial and radial errors carry opposite
signs and partially cancel there, so refining either axis **alone** makes the
reported number worse. The cell below measures that on three knobs the break
table has no row for — the grid, the fit window, and the domain length. A break
table cannot catch a baseline that is right by accident, so this has to be
measured directly and stated.

It also fits the far field of the solver's **area** mean, which is the amplitude
the 119 % row in §8f rests on and which otherwise has no independent witness."""))

cells.append(code('''print("3b. The published mu_1 agreement is a cancellation, not a converged error")
print(f"    {'n_z':>6}{'n_r':>6}{'mu_1':>12}{'rel err':>11}{'p_1 rel err':>13}")
GRID_SWEEP = {}
for nz, nr in ((600, 128), (600, 256), (1200, 128), (1200, 256)):
    s_ = solve_graetz(n_z=nz, n_r=nr)
    m_, p_ = decay(s_)
    GRID_SWEEP[(nz, nr)] = abs(m_ / MU[0] - 1)
    tag = "   <-- published" if (nz, nr) == (600, 128) else ""
    print(f"    {nz:6d}{nr:6d}{m_:12.6f}{abs(m_/MU[0]-1):11.2e}"
          f"{abs(p_/P_CUP[0]-1):13.2e}{tag}")
_pub = GRID_SWEEP[(600, 128)]
print(f"    Refining the RADIAL axis alone makes it "
      f"{GRID_SWEEP[(600, 256)]/_pub:.1f}x worse and the AXIAL axis")
print(f"    alone {GRID_SWEEP[(1200, 128)]/_pub:.1f}x worse; refining both together "
      "returns it to")
print(f"    {GRID_SWEEP[(1200, 256)]:.1e}. So the honest statement of the solver's "
      "accuracy at the")
print(f"    published grid is ~{max(GRID_SWEEP[(600, 256)], GRID_SWEEP[(1200, 128)]):.1e} "
      f"in mu_1; the {_pub:.1e} reported as `mu1_rel_err_2d`")
print("    is a partial cancellation, and reading it as an accuracy would")
print("    overstate the solve by an order of magnitude. What IS")
print("    load-bearing is the pair of observed orders above, which are read from")
print("    ratios and are insensitive to where the two error signs happen to meet.")

print(f"\\n    same message from two knobs with no break-table row:")
print(f"      {'fit window xmin':>22}{'mu_1 rel err':>14}")
for xm in (0.35, 0.5, 0.7):
    print(f"      {xm:22.2f}{abs(decay(sol, xmin=xm)[0]/MU[0]-1):14.2e}")
print(f"      {'domain length X':>22}{'mu_1 rel err':>14}")
for XX in (0.7, 1.0, 1.5, 2.0):
    print(f"      {XX:22.1f}"
          f"{abs(decay(solve_graetz(X=XX))[0]/MU[0]-1):14.2e}")

# The area amplitude, which nothing else on this page checks.
MU1_AREA_2D, P1_AREA_2D = decay(sol, key="area")
P1_AREA_ERR = abs(P1_AREA_2D / P_AREA[0] - 1)
print("\\n    and the AREA mean, the amplitude the 119 % row in 8f rests on:")
print(f"      2-D solve, far field of theta_area   p_1,area = {P1_AREA_2D:.6f}")
print(f"      series projection with weight s ds   p_1,area = {P_AREA[0]:.6f}")
print(f"      agree to {P1_AREA_ERR:.1e} -- the same two routes, no shared code, so")
print("      the counter-direction result in 8f has an independent witness too.")'''))

cells.append(code('''# 4. The inlet. The eigenfunction series is useless here -- it needs O(x^-1/2)
#    terms as x -> 0 -- but the boundary layer has its own similarity solution,
#    and its constant can be derived rather than quoted:
#       u ~ 4 u_bar y / R near the wall  =>  Nu -> 2/(Gamma(4/3) 9^(1/3)) x*^-1/3
print("4. The entrance region: local Nu against the Leveque limit")
print(f"   derived constant 2/(Gamma(4/3) 9^(1/3)) = {C_LEVEQUE:.6f}")
print(f"     {'n_z':>6}{'n_r':>6}{'extrapolated C0':>18}{'rel err':>11}{'order':>8}")
prev, LEV = None, None
for n_z, n_r in ((100, 40), (200, 80), (400, 160)):
    s_ = solve_graetz(X=1e-3, n_z=n_z, n_r=n_r, z_grade=6.0, r_grade=6.0)
    xx, nn = nusselt(s_)
    L = nn * (xx / 2) ** (1 / 3)
    i0, i1 = np.argmin(abs(xx - 1e-8)), np.argmin(abs(xx - 1e-7))
    c1 = (L[i1] - L[i0]) / (xx[i1] ** (1 / 3) - xx[i0] ** (1 / 3))
    c0 = L[i0] - c1 * xx[i0] ** (1 / 3)
    e = abs(c0 / C_LEVEQUE - 1)
    o = "" if prev is None else f"{np.log2(prev/e):8.2f}"
    print(f"     {n_z:6d}{n_r:6d}{c0:18.6f}{e:11.2e}{o:>8}")
    prev, LEV = e, e
LEVEQUE_ERR = LEV
print("   The extrapolation is in x^(1/3), the next term of the entrance")
print("   expansion; C0 is read from x = 1e-8 and 1e-7, deep inside the layer.")
print(f"\\n   and the other end: asymptotic Nu = mu_1^2/2 = {MU[0]**2/2:.5f}")
print("   -- the value every heat-transfer text prints for a tube at constant")
print("   wall temperature IS the first Graetz eigenvalue, halved.")'''))

cells.append(code('''# 5. Graetz's own calibration, re-run. Three equations (his I), three unknowns.
#    This tests the transcription and his arithmetic. It CANNOT test the model:
#    zero degrees of freedom means it fits whatever it is given.
def theta_series(kc, l, G, mus, ps, half=2.0):
    return float(np.sum(ps * np.exp(-kc * np.pi * l * mus ** 2 / (half * G))))


def invert(th, l, G, mus, ps, half=2.0):
    return brentq(lambda kc: theta_series(kc, l, G, mus, ps, half) - th,
                  1e-8, 20.0, xtol=1e-15, rtol=1e-15)


GMU = np.array([C["graetz_root_mu1"], C["graetz_root_mu2"]])
GP = np.array([C["p1_working_formula"], C["p2_working_formula"]])


def refit(mu2=None):
    m = np.array([GMU[0], GMU[1] if mu2 is None else mu2])
    def res(v):
        return [v[0] * np.exp(-v[2] * np.pi * l * m[0] ** 2 / (2 * G))
                + v[1] * np.exp(-v[2] * np.pi * l * m[1] ** 2 / (2 * G)) - t
                for t, l, G in zip(CAL.theta, CAL.tube_length_cm,
                                   CAL.mass_flow_g_min)]
    with warnings.catch_warnings():          # 8d feeds it a deliberate nonsense
        warnings.simplefilter("ignore")      # mu_2, which fsolve rightly dislikes
        return fsolve(res, [0.9, 0.012, 0.09])


p1f, p2f, kcf = refit()
LOGP1_ERR = abs(np.log10(p1f) + 10 - C["log_p1"])
print("5. His calibration re-solved from the three starred rows")
print(f"   {'':10}{'p_1':>12}{'p_2':>12}{'k/c':>10}{'log p_1':>12}{'log p_2':>10}")
print(f"   {'refitted':10}{p1f:12.5f}{p2f:12.5f}{kcf:10.5f}"
      f"{np.log10(p1f)+10:12.5f}{np.log10(p2f)+10:10.5f}")
print(f"   {'printed':10}{C['p1_working_formula']:12.5f}"
      f"{C['p2_working_formula']:12.5f}{0.0969:10.4f}"
      f"{C['log_p1']:12.5f}{C['log_p2']:10.5f}")
print(f"   log p_1 recovered to {LOGP1_ERR:.1e} in the logarithm; p_2 to "
      f"{abs(p2f/C['p2_working_formula']-1)*100:.1f} %, which is his own")
print("   iteration ('koennen nachher noch verbessert werden' -- can be improved")
print("   afterwards). The transcription and the reading of his formula are right.")
p1f2, p2f2, kcf2 = refit(mu2=MU[1])
print(f"   with the CORRECT mu_2 = {MU[1]:.4f} the same three rows give "
      f"p_1 = {p1f2:.5f}, p_2 = {p2f2:.5f}, k/c = {kcf2:.5f}")'''))

cells.append(code('''# 6. The misprint on page 90, settled against his own 51 rows.
print("6. The exponent printed in his working formula for k (p. 90)")
print(f"   {'mu_1 used':>28}{'mean |dev| from his k/c column':>34}")
MISPRINT = {}
for lab, m1 in ((f"{C['graetz_root_mu1']} (his p. 87 root)", C["graetz_root_mu1"]),
                (f"{C['mu1_working_formula']} (as printed on p. 90)",
                 C["mu1_working_formula"]),
                (f"{C['brass_corrected_root']} (his p. 86 gamma_1)",
                 C["brass_corrected_root"])):
    kc = np.array([invert(t, l, G, np.array([m1, GMU[1]]), GP)
                   for t, l, G in zip(obs.theta, obs.tube_length_cm,
                                      obs.mass_flow_g_min)])
    d_ = np.mean(np.abs(kc / obs.k_over_c_printed - 1))
    MISPRINT[round(m1, 4)] = d_
    print(f"   {lab:>28}{d_*100:32.2f} %")
FORMULA_DEV = MISPRINT[round(C["graetz_root_mu1"], 4)]
# how far apart the two printings are, taken from the CSV's own page column
PGAP = PAGE_OF["mu1_working_formula"] - PAGE_OF["brass_corrected_root"]
print(f"   Decisive: {C['mu1_working_formula']} is a misprint for "
      f"{C['graetz_root_mu1']}. {C['brass_corrected_root']} is the")
print(f"   brass-wall-corrected Bessel root {PGAP} pages earlier "
      f"(journal p. {PAGE_OF['brass_corrected_root']}")
print(f"   against p. {PAGE_OF['mu1_working_formula']}), a plausible "
      "compositor's substitution -- but the data,")
print("   not the guess, settle it.")'''))

cells.append(code('''# 7. Axial conduction. He keeps it in the plug-flow equation, then argues it is
#    negligible: alpha*beta/a^2 + beta^2 = sigma^2/R^2, and dropping beta^2
#    "beeinflusst den Werth von beta_1 erst in der sechsten Stelle".
#    That is a claim with a number in it, so check it at his own conditions.
R_TUBE = C["tube_diameter_cm"] / 2.0


def axial_effect(row):
    """|beta_1(approx)/beta_1(exact) - 1| from his own quadratic, one run."""
    V_dot = row.mass_flow_g_min / row.density_g_cm3       # cm^3/min
    a2 = row.k_printed / row.density_g_cm3                # c ~ 1 for water
    alpha = V_dot / (np.pi * R_TUBE ** 2)                 # cm/min
    q = alpha / a2
    b_ex = (-q + np.sqrt(q ** 2 + 4 * SIGMA[0] ** 2 / R_TUBE ** 2)) / 2
    b_ap = a2 * SIGMA[0] ** 2 / (alpha * R_TUBE ** 2)
    return a2, alpha, 2 * alpha * R_TUBE / a2, b_ex, b_ap, abs(b_ap / b_ex - 1)


# The effect FALLS as the flow rises, so his claim has to be tested at his
# SLOWEST run, not a convenient fast one. Both ends are printed.
wat = obs[obs.liquid == "water"]
row_slow = wat.loc[wat.mass_flow_g_min.idxmin()]
row_fast = wat.loc[wat.mass_flow_g_min.idxmax()]
print("7. Axial conduction, at the two extremes of his water runs")
print(f"   {'':22}{'G, g/min':>10}{'head':>6}{'alpha, cm/min':>15}{'Pe':>9}"
      f"{'rel effect':>13}")
for lab, r_ in (("slowest (worst case)", row_slow), ("fastest", row_fast)):
    a2, alpha, pe, b_ex, b_ap, rel = axial_effect(r_)
    print(f"   {lab:22}{r_.mass_flow_g_min:10.4f}{r_.pressure_head:>6}"
          f"{alpha:15.0f}{pe:9.0f}{rel:13.2e}")
a2, alpha, PE_AX, beta_exact, beta_approx, AXIAL_REL = axial_effect(row_slow)
print(f"   At the slowest run: a^2 = {a2:.4f} cm2/min, "
      f"beta_1 exact {beta_exact:.9f} /cm,")
print(f"   neglecting beta^2 {beta_approx:.9f} /cm, relative effect "
      f"{AXIAL_REL:.2e}.")
print(f"   He says 'erst in der sechsten Stelle' (only in the sixth place); at his")
print(f"   WORST run it reaches the {int(np.ceil(-np.log10(AXIAL_REL)))}th "
      "significant figure, which is exactly what he")
print("   claims, and at his fastest it is a further order down. The metric below")
print("   is the worst case, not a convenient one.")
print("\\n   And in the 2-D solve, where the term is 1/Pe^2 * d2theta/dx2. The")
print("   difference between two solves on the SAME grid isolates it: upwind")
print("   numerical diffusion is common to both and cancels.")
# first-order perturbation of the Poiseuille eigenvalue problem itself:
#   (1/s)(sV')' + lambda(1-s^2)V + (lambda^2/Pe^2)V = 0, lambda = mu^2
# gives  d(mu)/mu = -(mu^2/2Pe^2) * <V^2> / <(1-s^2)V^2>, both with weight s ds.
s1, V1, _ = V_profile(MU[0])
I_RATIO = (np.trapezoid(V1 ** 2 * s1, s1)
           / np.trapezoid((1 - s1 ** 2) * V1 ** 2 * s1, s1))
C_PRED = -MU[0] ** 2 / 2 * I_RATIO
print(f"     {'1/Pe^2':>12}{'mu_1':>12}{'shift':>12}{'shift x Pe^2':>14}")
base = decay(solve_graetz(n_z=300, n_r=64))[0]
C_OBS = None
for pe in (20.0, 40.0, 80.0):
    m_ = decay(solve_graetz(n_z=300, n_r=64, inv_pe2=1 / pe ** 2))[0]
    C_OBS = (m_ / base - 1) * pe ** 2
    print(f"     {1/pe**2:12.3e}{m_:12.6f}{m_/base-1:12.2e}{C_OBS:14.4f}")
AXIAL_PERT_ERR = abs(C_OBS / C_PRED - 1)
print(f"   The shift is C/Pe^2 with C converging to {C_OBS:.4f}. Perturbing the")
print(f"   eigenvalue problem analytically predicts C = -mu_1^2/2 * "
      f"<V^2>/<(1-s^2)V^2> = {C_PRED:.4f},")
print(f"   agreeing to {AXIAL_PERT_ERR*100:.1f} % -- a check that shares no code with the")
print("   finite-volume solve and that a wrong sign or a wrong Pe scaling breaks.")
print(f"   At his Pe = {PE_AX:.0f} that is {abs(C_OBS)/PE_AX**2:.1e} in mu_1, far below")
print("   the discretisation error of any affordable grid, and the same verdict as")
print("   his own quadratic reaches for the plug-flow case.")

# 7b. Is the flow laminar? He never says, and the model needs it. The three
#     printed pressure heights and the measured flows answer it with no external
#     constant at all: Poiseuille gives V proportional to head, turbulent pipe
#     flow roughly to head^0.5.
print("\\n7b. The flow regime, from his printed pressure heights alone")
heads = {"I": C["head_upper_cm"], "II": C["head_middle_cm"], "III": C["head_lower_cm"]}
print(f"   {'liquid':12s}{'heads used':>12}{'exponent n in V ~ head^n':>28}")
EXPO = []
for liq, g in obs.groupby("liquid", sort=False):
    m = g.groupby("pressure_head").apply(
        lambda d: (d.mass_flow_g_min / d.density_g_cm3).mean(), include_groups=False)
    if len(m) < 2:
        continue
    h = np.array([heads[k] for k in m.index])
    n_ = np.polyfit(np.log(h), np.log(m.values), 1)[0]
    EXPO.append(n_)
    print(f"   {liq:12s}{'/'.join(m.index):>12}{n_:28.2f}")
print(f"   mean exponent {np.mean(EXPO):.2f} (Poiseuille would give 1, fully")
print("   turbulent pipe flow about 0.5). The flow is laminar, so the Poiseuille")
print("   PROFILE is the right family -- but this says nothing about whether it is")
print("   fully DEVELOPED over the tube, which is a separate question the paper")
print("   gives no way to answer: it would need the viscosity, which he does not")
print("   print, and the entrance length is what a developing profile turns on.")'''))

cells.append(md(r"""### 8. Defect injection — what each published number moves for

One row per metric, plus the rows that show which checks are *structural* and
say so. The break table does not travel when a page directory is copied, and
this one was rebuilt from scratch for this physics — `A2.3`'s rows (mean
subtraction, the 192/48 factor, slug mass conservation) have no counterpart
here."""))

cells.append(code('''print("8a. The 2-D solve  (mu1_rel_err_2d, p1_rel_err_2d, p1_area_rel_err_2d,")
print("    nu_infinity, order_*)")
print(f"    {'injected defect':40s}{'mu_1':>10}{'rel err':>10}"
      f"{'p_1':>10}{'rel err':>10}{'p_1,area err':>13}")
BREAK = {}
for lab, kw in (("as published (nu_r=1, Poiseuille, TVD)", {}),
                ("construct_div nu_r = 1 -> 0 (slab)", dict(nu_r=0)),
                ("plug velocity instead of Poiseuille", dict(plug=True)),
                ("TVD off: first-order upwind", dict(tvd=False)),
                ("n_r = 128 -> 8 (wall layer unresolved)", dict(n_r=8))):
    s_ = solve_graetz(**{"n_z": 600, "n_r": 128, **kw})
    mu_, p_ = decay(s_)
    pa_ = decay(s_, key="area")[1]
    BREAK[lab] = (mu_, p_, pa_)
    print(f"    {lab:40s}{mu_:10.5f}{abs(mu_/MU[0]-1):10.2e}"
          f"{p_:10.5f}{abs(p_/P_CUP[0]-1):10.2e}"
          f"{abs(pa_/P_AREA[0]-1):13.2e}")
print("    Every row moves both numbers far outside the published error. The plug")
print("    row is not a defect so much as Graetz's own first model, and it lands")
print(f"    on his sigma_1 = {SIGMA[0]:.4f} and 4/sigma_1^2 = {P_PLUG[0]:.4f} instead -- the")
print("    sharpest confirmation that the velocity profile is what sets both.")
print("    And the order metrics (order_mu1_axial, order_mu1_radial) move too:")
e_ = [abs(decay(solve_graetz(n_z=n, n_r=256, tvd=False))[0] / MU[0] - 1)
      for n in (75, 300)]
print(f"      TVD off, n_z = 75 -> 300: err {e_[0]:.2e} -> {e_[1]:.2e}, "
      f"observed order {0.5*np.log2(e_[0]/e_[1]):.2f}")
print(f"      TVD on  (as published):                        "
      f"observed order {LAD['n_z']:.2f}")

print("\\n8b. The Leveque constant  (leveque_const_rel_err)")
print(f"    {'injected defect':46s}{'C0':>12}{'rel err':>11}")
for lab, kw in (("as published", {}),
                ("plug velocity: no wall shear at all", dict(plug=True)),
                ("uniform radial mesh (r_grade 6 -> 1)", dict(r_grade=1.0)),
                ("construct_div nu_r = 1 -> 0", dict(nu_r=0))):
    s_ = solve_graetz(**{"X": 1e-3, "n_z": 200, "n_r": 80, "z_grade": 6.0,
                        "r_grade": 6.0, **kw})
    xx, nn = nusselt(s_)
    L = nn * (xx / 2) ** (1 / 3)
    i0, i1 = np.argmin(abs(xx - 1e-8)), np.argmin(abs(xx - 1e-7))
    c1 = (L[i1] - L[i0]) / (xx[i1] ** (1 / 3) - xx[i0] ** (1 / 3))
    c0 = L[i0] - c1 * xx[i0] ** (1 / 3)
    print(f"    {lab:46s}{c0:12.6f}{abs(c0/C_LEVEQUE-1):11.2e}")
print("    Leveque exists because u ~ 4*u_bar*y/R at the wall. Plug flow has none,")
print("    so the x^-1/3 law is the wrong law there and the extrapolated constant")
print("    is meaningless -- 20x out. A uniform radial mesh cannot resolve a layer")
print("    that is 1e-3 of the radius thick, and misses by orders.")
print("    THE nu_r ROW IS POWERLESS, AND FOR A REASON: at x = 1e-8 the thermal")
print("    layer is ~1e-3 R thick, so the wall is locally flat and the leading")
print("    entrance constant cannot see the curvature. This metric tests the wall")
print("    shear and the near-wall resolution; it does NOT test the geometry, and")
print("    8a is where the geometry is tested.")'''))

cells.append(code('''print("8c. The series coefficients  (series_coeff_mu8_err, series_coeff_mu10_err)")
print("    These come from exact rational arithmetic, so a 'defect' has to be in")
print("    the recurrence. Dropping the (1 - x^2/mu^2) weight turns the equation")
print("    into Bessel's, and the polynomial must change:")
def bessel_poly():
    c, out = [Fraction(1)], []
    for m in range(30):
        c.append(-c[m] / Fraction((2 * m + 2) ** 2))
    return [float(c[k]) for k in range(6)]
bp = bessel_poly()
print(f"    {'power':>7}{'Graetz weight (1-s^2)':>24}{'weight dropped (Bessel)':>26}")
for k in range(4):
    print(f"    mu^{2*k:<4d}{exact_poly[k]:24.10g}{bp[k]:26.10g}")
print(f"    and the first root moves {MU[0]:.5f} -> {SIGMA[0]:.5f}. So the check")
print("    that his printed 0.1875 / 0.007921 / 0.00014404 are right can fail,")
print("    and does, for the one substitution that matters.")

print("\\n8d. His calibration refit  (calibration_logp1_abs_err)")
print(f"    {'injected defect':46s}{'p_1':>11}{'log p_1':>11}{'abs err':>11}")
for lab, m2 in (("as published (mu_2 = 6.50, his value)", GMU[1]),
                (f"mu_2 = {MU[1]:.4f} (the true second eigenvalue)", MU[1]),
                ("mu_2 = 3.00 (nonsense)", 3.0)):
    v = refit(mu2=m2)
    print(f"    {lab:46s}{v[0]:11.5f}{np.log10(v[0])+10:11.5f}"
          f"{abs(np.log10(v[0])+10-C['log_p1']):11.2e}")
print("    The refit recovers his log p_1 only when it is given HIS mu_2. That is")
print("    the check working: it is a check on the transcription and on his")
print("    arithmetic, and it is NOT a check on the model -- three constants")
print("    fitted to three rows fit anything. See 8f.")

print("\\n8e. The misprint test  (formula_reproduces_kc)")
for m1, d_ in sorted(MISPRINT.items()):
    print(f"    mu_1 = {m1:<8} mean |dev| from his k/c column {d_*100:8.2f} %")'''))

cells.append(code('''# 8f. THE STRUCTURAL ROW, and the point of the page.
print("8f. The spread of k/c across his three pressure heads")
print("    (the fitted row is STRUCTURAL -- three constants fitted to three rows --")
print("     so it is exactly 0 and is deliberately NOT reported as a metric)")
print(f"    {'route':44s}{'head III':>10}{'head II':>10}{'head I':>10}{'spread':>9}")
SPREAD = {}
routes = [("his fitted p_1, p_2 (nothing else known)", GMU, GP, 2.0),
          ("computed p_i, 1 term", MU[:1], P_CUP[:1], 2.0),
          ("computed p_i, 2 terms", MU[:2], P_CUP[:2], 2.0),
          ("computed p_i, 10 terms -- nothing fitted", MU, P_CUP, 2.0),
          ("computed p_i, 10 terms, AREA mean", MU, P_AREA, 2.0),
          ("plug flow, his own p. 84 result", SIGMA, P_PLUG, 1.0)]
order = CAL.mass_flow_g_min.values.argsort()
KC_ROUTE = {}
for lab, mus, ps, half in routes:
    v = [invert(CAL.theta[i], CAL.tube_length_cm[i], CAL.mass_flow_g_min[i],
                mus, ps, half) for i in order]
    sp = (max(v) - min(v)) / np.mean(v)
    SPREAD[lab] = (np.mean(v), sp)
    KC_ROUTE[lab] = v
    print(f"    {lab:44s}{v[0]:10.5f}{v[1]:10.5f}{v[2]:10.5f}{sp*100:8.2f}%")

EXACT_LAB = "computed p_i, 10 terms -- nothing fitted"
SPREAD_EXACT = SPREAD[EXACT_LAB][1]
KC10 = KC_ROUTE[EXACT_LAB]
X_CAL = [KC10[j] * np.pi * CAL.tube_length_cm[i] / (2 * CAL.mass_flow_g_min[i])
         for j, i in enumerate(order)]

# Two things the spread could be an artefact of, and is neither. Both are cheap
# and both were missing from the first version of this page.
print("\\n    Is it an artefact of truncating the series at 10 modes?")
MU_HI = eigen(n_modes=20, hi=70.0)
P_CUP_HI = amplitudes(MU_HI)[0]
v20 = [invert(CAL.theta[i], CAL.tube_length_cm[i], CAL.mass_flow_g_min[i],
              MU_HI, P_CUP_HI) for i in order]
SPREAD_20 = (max(v20) - min(v20)) / np.mean(v20)
print(f"      {'computed p_i, 20 terms':44s}{v20[0]:10.5f}{v20[1]:10.5f}"
      f"{v20[2]:10.5f}{SPREAD_20*100:8.2f}%")
print(f"      The three k/c agree with the 10-term values to "
      f"{max(abs(a / b - 1) for a, b in zip(v20, KC10)):.1e}. His three runs")
print(f"      sit at x = {', '.join(f'{xx:.3f}' for xx in X_CAL)}, where the 11th mode "
      f"is damped by")
print(f"      exp(-mu_11^2 x) <= {np.exp(-MU_HI[10]**2*min(X_CAL)):.0e} "
      f"(mu_11 = {MU_HI[10]:.2f}). Truncation is not the cause.")

print("\\n    Does it depend on the eigenfunction route at all? Invert the same")
print("    three rows through the 2-D SOLVE's own theta_cup(x) curve instead --")
print("    no eigenfunction anywhere, only the finite-volume field:")
_LOGCUP = CubicSpline(np.log(sol["x"]), np.log(sol["cup"]))


def invert_2d(th, l, G):
    """k/c from theta_cup(x) of the 2-D solve alone. x = (k/c) pi l / (2 G)."""
    lx = brentq(lambda t: float(_LOGCUP(t)) - np.log(th),
                np.log(sol["x"][0]), np.log(sol["x"][-1]),
                xtol=1e-15, rtol=1e-15)
    return np.exp(lx) * 2.0 * G / (np.pi * l)


v2d = [invert_2d(CAL.theta[i], CAL.tube_length_cm[i], CAL.mass_flow_g_min[i])
       for i in order]
SPREAD_2D = (max(v2d) - min(v2d)) / np.mean(v2d)
print(f"      {'2-D theta_cup(x), no series at all':44s}{v2d[0]:10.5f}"
      f"{v2d[1]:10.5f}{v2d[2]:10.5f}{SPREAD_2D*100:8.2f}%")
print(f"      against the series route's {SPREAD_EXACT*100:.2f} %. The central claim "
      "of this page does")
print("      not depend on the eigenfunction calculation.")

# how big is the thermometry error on that spread?
def amplify(th, l, G):
    k0 = invert(th, l, G, MU, P_CUP)
    k1 = invert(th * (1 + 1e-6), l, G, MU, P_CUP)
    return abs((k1 / k0 - 1) / 1e-6)


div = C["thermometer_division_inlet_outlet"]        # E1, E2: 1/10 degree, p. 82
div_bath = C["thermometer_division_bath"]          # E3:     1/5  degree, p. 82


def sigma_spread(coarsen=1.0, e3=False, mass_pct=0.0, average=True):
    """sigma on the head-to-head spread of k/c, for one error model.

    theta = (U - T0)/(T1 - T0), so a half-division on E1 (T1) and E2 (U) enters
    as 0.5*div*sqrt(1 + theta^2)/(T1 - T0); the bath thermometer E3 enters as
    0.5*div_bath*|theta - 1|/(T1 - T0). `mass_pct` is a relative error on the
    weighed-and-timed mass flow, which passes straight through because
    x = (k/c) pi l / (2 G) holds theta fixed: dln(k/c) = dln G exactly.
    `average` divides the READING errors by sqrt(5) -- each starred row is the
    mean of five observations (p. 91). A weighing or timing bias does not
    average, so `mass_pct` never gets the sqrt(5).
    """
    sig = []
    for i in range(3):
        r_ = CAL.iloc[i]
        A_ = amplify(r_.theta, r_.tube_length_cm, r_.mass_flow_g_min)
        d2 = (0.5 * div * coarsen) ** 2 * (1.0 + r_.theta ** 2)
        if e3:
            d2 += (0.5 * div_bath * coarsen) ** 2 * (r_.theta - 1.0) ** 2
        dth = np.sqrt(d2) / (r_.T_in_C - r_.T_bath_C)
        s_ = A_ * dth / r_.theta
        if average:
            s_ = s_ / np.sqrt(5)
        sig.append(np.hypot(s_, mass_pct / 100.0) if mass_pct else s_)
    return float(np.hypot(max(sig), min(sig)))


SIG_SPREAD = sigma_spread()
print(f"\\n    half a scale division ({div/2:.2f} K) on E1 and E2, five observations")
print(f"    per row, propagates to {SIG_SPREAD*100:.1f} % on that spread.")
print(f"    The 10-term spread is {SPREAD_EXACT*100:.1f} %, i.e. "
      f"{SPREAD_EXACT/SIG_SPREAD:.0f} sigma.")
print("    That is the NARROWEST defensible error budget -- it carries E1 and E2")
print("    and nothing else. 8g widens it, names what it leaves out and gives the")
print("    range, because the honest headline is a range and not the single")
print("    number 20.")
print("    The fitted route reads 0.00 % by construction and therefore proves")
print("    nothing; it is listed as the structural row it is.")

# 8g. The metrics that are not the 2-D solve.
print("8g. The remaining metrics")
print("    axial_conduction_rel_effect -- the cross-check in section 7 is the")
print("    perturbation constant C. It is sensitive to the weight it is projected")
print("    against, which is the thing that could be got wrong:")
_s1, _V1, _ = V_profile(MU[0])
_right = (np.trapezoid(_V1 ** 2 * _s1, _s1)
          / np.trapezoid((1 - _s1 ** 2) * _V1 ** 2 * _s1, _s1))
_wrong = 1.0            # projecting against the area weight, i.e. forgetting (1-s^2)
for lab, I_ in (("as published: weight (1-s^2) s ds", _right),
                ("weight s ds (the area weight, wrong)", _wrong)):
    print(f"      {lab:40s} C = {-MU[0]**2/2*I_:8.4f}   "
          f"vs the solve's {C_OBS:.4f}  ({abs((-MU[0]**2/2*I_)/C_OBS-1)*100:5.1f} %)")

print(f"\\n    kc_spread_thermometry_sigma -- the {SPREAD_EXACT*100:.1f} % spread in 8f "
      "against the error")
print("    budget. FIRST, what the published budget LEAVES OUT. It carries E1 and")
print("    E2 only. Two sources are missing, and adding either is legitimate:")
print(f"      {'error model':52s}{'sigma':>8}{'significance':>14}")
BUDGET = {}
for lab, kw in (
        ("published: E1,E2 half-division, 5 obs", {}),
        (f"+ E3, the {div_bath} K bath thermometer (p. 82)", dict(e3=True)),
        ("+ E3 + 1 % on the weighed-and-timed mass",
         dict(e3=True, mass_pct=1.0)),
        ("+ E3, and no sqrt(5) averaging at all",
         dict(e3=True, average=False))):
    sg = sigma_spread(**kw)
    BUDGET[lab] = sg
    print(f"      {lab:52s}{sg*100:7.2f} %{SPREAD_EXACT/sg:11.1f} sigma")
SIG_SPREAD_WIDE = max(BUDGET.values())
print(f"    E3's {div_bath} K division is in this page's own constants file and was")
print("    never in the budget; the mass flow -- a weighing over a timed minute,")
print(f"    p. 82 -- was not in it at all. So the honest headline is a RANGE: "
      f"{SPREAD_EXACT/SIG_SPREAD:.0f} sigma")
print(f"    on his stated thermometry down to {SPREAD_EXACT/SIG_SPREAD_WIDE:.0f} "
      "sigma on the widest budget that can")
print("    be constructed from what the paper says. The conclusion survives all of")
print("    them, and the single number 20 is the optimistic end.")
print("    NOTE WHAT THIS TABLE CANNOT DO. The rows below only make an ASSUMED")
print("    source coarser. No coarsening row can ADD a source, so a break table of")
print("    this shape is structurally incapable of detecting the understatement")
print("    above -- which is exactly why the omissions are named rather than")
print("    tested.")

print("\\n    SECOND, coarsening the assumed source, which is what a break row can do:")
for scale, avg, lab in (
        (1.0, True, "as published: 1/10 degree scale, 5 obs"),
        (10.0, True, "1 degree scale, still 5 obs (10x coarser)"),
        (1.0, False, "1/10 degree scale, a single observation"),
        (10.0, False, "1 degree scale AND a single observation")):
    sg = sigma_spread(coarsen=scale, average=avg)
    print(f"      {lab:52s}{sg*100:7.2f} %{SPREAD_EXACT/sg:11.1f} sigma")
print(f"    The 1-degree rows are a STRAWMAN and are printed only to bound the")
print(f"    table: p. 82 prints E1 and E2 as divided in 1/{int(round(1/div))} degree, and his own")
print("    tables quote temperatures to 0.01 K. They are not a scenario the paper")
print("    admits.")
print("    They are also where an earlier version of this page argued against")
print("    itself, so it is worth being plain: on the last row the spread would sit")
print(f"    at {SPREAD_EXACT/sigma_spread(coarsen=10.0, average=False):.1f} sigma, "
      "which is not survival -- it is the threshold of")
print("    insignificance. If that row were the truth the page would have no")
print("    result. It is not the truth, and the reason is p. 82, not the arithmetic.")
print("    What actually carries the conclusion is not sigma at all:")
print("      (i)  the spread is MONOTONE in flow rate (8f, low flow -> high k/c),")
print("           which reading noise has no reason to be; and")
print("      (ii) it reproduces OUT OF SAMPLE on liquids and runs that were never")
print("           fitted -- see 8h, which is the stronger of the two.")

print("\\n    head_I_over_III_amplification -- it is a property of the inversion,")
print("    so it changes with the constants inverted through:")
for lab, mus, ps in (("computed p_i, 10 terms (as published)", MU, P_CUP),
                     ("his fitted two-term formula", GMU, GP),
                     ("plug flow", SIGMA, P_PLUG)):
    a_ = []
    for h in ("I", "III"):
        g = obs[obs.pressure_head == h]
        half = 1.0 if lab == "plug flow" else 2.0
        vals = []
        for _, r in g.iterrows():
            k0 = invert(r.theta, r.tube_length_cm, r.mass_flow_g_min, mus, ps, half)
            k1 = invert(r.theta * (1 + 1e-6), r.tube_length_cm, r.mass_flow_g_min,
                        mus, ps, half)
            vals.append(abs((k1 / k0 - 1) / 1e-6))
        a_.append(np.mean(vals))
    print(f"      {lab:40s} head I / head III = {a_[0]/a_[1]:5.2f}")
print("    The conclusion -- head I is the worst-conditioned, not the best -- is")
print("    the same whichever constants are used, which is why it survives the")
print("    disagreement the rest of this page is about.")'''))

cells.append(md(r"""### 8h. The out-of-sample test — the one Graetz sets himself

Everything above is derived on the **three starred rows**, which are also the
three rows his $p_1$, $p_2$ and $k/c$ were fitted to. That is the weakest place
to make the argument, and it is not necessary: he states the test himself, twice.
On page 89 he asks whether it matters *"ob man nun die oberste, mittelste oder
unterste Druckhöhe benutzte"* (whether one uses the upper, middle or lower
pressure head), and on page 90 he asserts it does not.

So run it. For each liquid measured at more than one head, invert every run and
compare the head-to-head means — with his fitted constants, and with the computed
ones. **Forty-eight of the fifty-one rows were never fitted to anything.**"""))

cells.append(code('''print("8h. His own consistency test, run on every liquid he measured at more")
print("    than one pressure head. Spread = (max - min)/mean of the per-head mean")
print("    k/c. The three starred rows are 3 of the 51; everything else is out of")
print("    sample.")
print(f"    {'liquid':12s}{'runs':>6}{'heads':>12}"
      f"{'spread, his fitted p_i':>24}{'spread, computed p_i':>22}")
OOS = {}
for liq, g in obs.groupby("liquid", sort=False):
    heads = "/".join(sorted(g.pressure_head.unique()))
    if g.pressure_head.nunique() < 2:
        print(f"    {liq:12s}{len(g):6d}{heads:>12}"
              f"{'-- one head only, cannot test --':>46}")
        continue
    row_out = []
    for mus, ps in ((GMU, GP), (MU, P_CUP)):
        kc = np.array([invert(r.theta, r.tube_length_cm, r.mass_flow_g_min, mus, ps)
                       for _, r in g.iterrows()])
        m_ = pd.Series(kc, index=g.pressure_head.values).groupby(level=0).mean()
        row_out.append((m_.max() - m_.min()) / m_.mean())
    OOS[liq] = tuple(row_out)
    print(f"    {liq:12s}{len(g):6d}{heads:>12}"
          f"{row_out[0]*100:23.1f} %{row_out[1]*100:21.1f} %")
print(f"\\n    With his fitted constants the test passes, as he says it does: "
      f"{max(v[0] for v in OOS.values())*100:.1f} % at")
print(f"    worst. With the computed constants it fails on every liquid that has "
      "three")
print(f"    heads -- water {OOS['water'][1]*100:.1f} %, alcohol "
      f"{OOS['alcohol'][1]*100:.1f} %, and all {len(obs[obs.liquid=='cuso4'])} "
      f"copper-sulfate rows")
print(f"    {OOS['cuso4'][1]*100:.1f} % -- against {OOS['water'][0]*100:.1f} %, "
      f"{OOS['alcohol'][0]*100:.1f} % and {OOS['cuso4'][0]*100:.1f} % with his. "
      "The two-head")
print("    liquids span a much smaller flow range and separate the two routes much")
print("    less, which is itself the right direction.")
print(f"\\n    And the direction is the same every time: the lowest head gives the")
print("    highest k/c. Reading noise has no reason to do that on four liquids at")
print(f"    once, and this is not the fitted rows -- {len(obs)-len(CAL)} of the "
      f"{len(obs)} were never fitted.")
print("    This is the strongest evidence on the page that the 8f disagreement is")
print("    real, and it is stronger than any sigma in 8g, because it needs no error")
print("    model at all.")'''))

cells.append(md(r"""## What pymrm adds

Three things, and the third is not flattering to the exercise.

**1. The constants.** Graetz's $p_i$ are pure numbers — they do not depend on
$R$ at all, contrary to the sentence quoted above, because $R$ scales out of the
eigenvalue problem. They are computed here two independent ways.

**2. Two limits he did not have.** The entrance region and the far field both
have closed forms that the eigenfunction series cannot supply — the first
because it would need infinitely many terms, the second because it *is* the
first term, which he could see but not name. The asymptotic Nusselt number that
every heat-transfer text prints is $\mu_1^2/2$.

**3. What happens when his fitted constants are replaced by the true ones.** His
calibration fits three constants to three runs, so its perfect agreement is
guaranteed and is not evidence. Substituting the computed amplitudes removes the
fit — and the three runs then no longer agree on a single $k/c$. The
disagreement is monotone in flow rate, far outside his thermometry on any error
budget that can be built from the paper (§8g gives the range, and says what the
budget excludes), and — the part that needs no error budget at all — it
**reproduces on the forty-eight runs that were never fitted** (§8h). That last
test is Graetz's own, proposed on his page 89.

So the honest verdict runs the other way from the framing: **Graetz's empirical
calibration was doing real work, and it was not merely standing in for a
calculation he could not do.** It was absorbing a systematic error in the
apparatus that his own method could not see, and removing it makes his numbers
worse, not better.

### What this page cannot conclude

It does not identify the systematic effect. The spread is real — monotone in
flow rate and twenty times the propagated reading error — but *which* effect it
is stays open, and three candidates are worth writing down so the next reader
does not have to re-derive them:

- **A velocity profile that is not fully developed.** His printed heads and
  measured flows put $\dot V$ all but exactly proportional to head (§7b), so the
  flow is laminar and Poiseuille is the right family. Whether it is *developed* over
  the tube depends on the entrance length, which needs a viscosity the paper does
  not print, and settling it needs a momentum solve this page does not do. Note
  the direction: the other extreme, plug flow, moves the inferred $k/c$ further
  *down*, not up, so a flatter profile does not on its own close the gap.
- **Heat exchanged outside the cooled length** — in the corks, or in the glass
  tubes carrying the thermometers. Any flow-independent loss weighs more at low
  flow, which is the right direction, but the paper gives no dimensions for it.
- **The area-versus-cup mean.** Graetz's page-88 formula integrates over the
  cross-section; a collected sample is flow-weighted. §8f shows the two differ
  enormously — but the area mean makes the spread *worse*, so this is not the
  explanation either.

And the page does not claim Graetz's conductivities are wrong. It claims the
route by which he obtained them cannot be checked from inside the paper, that
checking it from outside makes it fail, and that his numbers survive because two
fitted constants absorbed the failure."""))

cells.append(code('''print("What replacing the two fitted constants does to his conductivities")
print("(k/c inverted run by run; k = c * k/c; his own reference values from p. 80")
print(" and p. 93 for scale, all in g-cal/(cm min K))")
obs["kc_exact"] = [invert(t, l, G, MU, P_CUP) for t, l, G in
                   zip(obs.theta, obs.tube_length_cm, obs.mass_flow_g_min)]
obs["kc_graetz"] = [invert(t, l, G, GMU, GP) for t, l, G in
                    zip(obs.theta, obs.tube_length_cm, obs.mass_flow_g_min)]
c_imp = obs.k_printed / obs.k_over_c_printed
summary = []
for liq, g in obs.groupby("liquid", sort=False):
    ci = (g.k_printed / g.k_over_c_printed).mean()
    summary.append((liq, len(g), g.k_printed.mean(), (g.kc_exact * ci).mean(),
                    (g.kc_exact * ci).mean() / g.k_printed.mean() - 1))
print(f"   {'liquid':12s}{'n':>3}{'his k':>10}{'k, unfitted':>13}{'change':>10}")
for liq, n, kh, ke, ch in summary:
    print(f"   {liq:12s}{n:3d}{kh:10.4f}{ke:13.4f}{ch*100:9.1f} %")
WATER_SHIFT = [s for s in summary if s[0] == "water"][0][4]
print(f"\\n   his own reference points for water, all printed in the paper:")
print(f"     Weber   (4 C)  {C['k_water_weber_4C']}      "
      f"Weber   (23 C) {C['k_water_weber_23C']}")
print(f"     Lorberg (23 C) {C['k_water_lorberg_23C']}     "
      f"Graetz  (30 C) {C['k_water_graetz_30C']}")
print(f"     Weber via Lorberg's temperature coefficient, carried to 30 C: "
      f"{C['k_water_lorberg_carried_30C']:.4f}")
print(f"   He notes that carried value, "
      "'was mit dem von mir bestimmten Werth gut uebereinstimmt'")
print(f"   -- {abs(C['k_water_lorberg_carried_30C']/C['k_water_graetz_30C']-1)*100:.1f} % "
      f"from his own {C['k_water_graetz_30C']}. The unfitted inversion gives "
      f"{[s for s in summary if s[0]=='water'][0][3]:.4f},")
print(f"   {abs(WATER_SHIFT)*100:.0f} % below both, and the agreement he had is gone.")'''))

cells.append(code('''# The other thing the 2-D solve settles: which of his three pressure heads is
# the good one. He states a preference, on p. 90, and it is the wrong way round.
print("Which pressure head gives the best-conditioned measurement?")
print("Graetz, p. 90: 'Die Beobachtungen bei der groessten Druckhoehe sind den")
print("anderen ueberlegen, weil die Ausflusstemperatur U sich bei diesen der")
print("Einflusstemperatur T_1 am meisten naehert.'  (The observations at the")
print("greatest pressure head are superior to the others, because the outflow")
print("temperature U there comes closest to the inflow temperature T_1.)")
print(f"\\n   {'head':>5}{'runs':>6}{'theta':>9}{'T_1 - U':>10}"
      f"{'|dln(k/c)/dln(theta)|':>23}{'sigma(k/c)':>12}")
AMP = {}
for h, g in obs.groupby("pressure_head"):
    a_ = np.mean([amplify(r.theta, r.tube_length_cm, r.mass_flow_g_min)
                  for _, r in g.iterrows()])
    dth = 0.5 * div * np.hypot(1.0, g.theta) / (g.T_in_C - g.T_bath_C)
    AMP[h] = a_
    print(f"   {h:>5}{len(g):6d}{g.theta.mean():9.4f}"
          f"{(g.T_in_C-g.T_out_C).mean():10.2f}{a_:23.2f}"
          f"{np.mean(a_*dth/g.theta)*100:11.1f} %")
print(f"\\n   His reason is sound as far as it goes: at head I the liquid changes")
print(f"   temperature by {(obs[obs.pressure_head=='I'].T_in_C-obs[obs.pressure_head=='I'].T_out_C).mean():.1f} K "
      f"rather than {(obs[obs.pressure_head=='III'].T_in_C-obs[obs.pressure_head=='III'].T_out_C).mean():.1f} K, so the T at which k is quoted is")
print("   better defined, and he is explicitly worried about the temperature")
print(f"   coefficient. But the inversion amplifies a temperature error by "
      f"{AMP['I']/AMP['III']:.1f}x")
print("   more at head I than at head III, and he does not weigh the two. The")
print("   trade-off is only visible once the inversion is done properly, which")
print("   needs the amplitudes he did not have.")

report_agreement("A3.15", {
    # the 2-D solve against the eigenfunction route -- no shared code
    # NOTE mu1_rel_err_2d is a partial cancellation of the axial and radial
    # errors at 600/128, not the accuracy of the solve -- see 3b. The orders are
    # the load-bearing pair.
    "mu1_rel_err_2d": float(abs(MU1_2D / MU[0] - 1)),
    "p1_rel_err_2d": float(abs(P1_2D / P_CUP[0] - 1)),
    "p1_area_rel_err_2d": float(P1_AREA_ERR),
    "nu_infinity": float(MU[0] ** 2 / 2),
    "order_mu1_axial": float(LAD["n_z"]),
    "order_mu1_radial": float(LAD["n_r"]),
    # the entrance region, its own limit and its own ladder
    "leveque_const_rel_err": float(LEVEQUE_ERR),
    # what the paper prints
    "series_coeff_mu8_printed_err": float(poly_err[4]),
    "series_coeff_mu10_printed_err": float(poly_err[5]),
    "mu1_printed_rel_err": float(mu1_err),
    "mu2_printed_rel_err": float(mu2_err),
    "calibration_logp1_abs_err": float(LOGP1_ERR),
    "formula_reproduces_kc_mean_dev": float(FORMULA_DEV),
    "axial_conduction_rel_effect": float(AXIAL_REL),
    # the result: the fitted rows, then the same test out of sample
    "kc_spread_calibration_unfitted": float(SPREAD_EXACT),
    "kc_spread_out_of_sample_water": float(OOS["water"][1]),
    "kc_spread_out_of_sample_cuso4": float(OOS["cuso4"][1]),
    "kc_spread_out_of_sample_alcohol": float(OOS["alcohol"][1]),
    # narrowest and widest defensible thermometry budgets -- the honest headline
    # is the range between them, not either endpoint
    "kc_spread_thermometry_sigma": float(SIG_SPREAD),
    "kc_spread_thermometry_sigma_wide": float(SIG_SPREAD_WIDE),
    "water_k_shift_unfitted": float(WATER_SHIFT),
    "head_I_over_III_amplification": float(AMP["I"] / AMP["III"]),
})'''))

cells.append(md(r"""## Reuse

**`solve_graetz` is the reusable object**, and almost nothing in it is about
heat. Replace $a^2$ by a diffusivity and it is the Graetz–Nusselt *mass* transfer
problem — a wall-coated monolith channel, a membrane tube, a
diffusion-limited wall reaction. The three switches are the ones you will want:
`nu_r` for the geometry, `plug` for the velocity profile, `inv_pe2` for axial
dispersion.

**Copy the deferred correction, and copy the under-relaxation with it.** A pure
convection–radial-diffusion problem has no physical axial diffusion, so
first-order upwind contributes error and nothing else, and it is first order. The
van Leer correction costs about 25 back-substitutions through one stored
factorisation and buys a factor of 16 in cells. Without `omega < 1` the limiter
switches and the iteration limit-cycles at the increment printed in the
implementation section — the solve still *returns*, and every number then depends
on the iteration count. The loop here raises rather than returning a
non-converged field.

**Read the amplitude table before lumping.** A "Nusselt number for a tube" of
$\mu_1^2/2$ is the $x\to\infty$ limit and it is wrong by an order of magnitude
in the entrance. The printed $p_i$ and $\mu_i$ above are all you need for the exact
answer at any $x$; the whole solution is five numbers.

**What to take from the calibration.** A fit with as many constants as
observations reproduces those observations exactly and tells you nothing, and
this is not a historical curiosity — it is `A2.3`'s lesson in another form and
the most common defect this repository finds. The test Graetz proposed for
himself was the right one, *"so mussten sich für die anderen Flüssigkeiten mit
diesen Constanten übereinstimmende Werthe von k ergeben"* (*then the other
liquids must yield agreeing values of k with these constants*) — the fit is
tested out of sample or not at all. What this page adds is a route that needs no
fit, and §8h actually runs his test with it, on every liquid he measured at more
than one pressure head. With his fitted constants it passes, as he says. With
the computed ones it fails on all three liquids that span three heads, in the
same direction every time, on runs that were never fitted to anything. That is
the strongest evidence on this page and it needs no error model at all — which
matters, because §8g shows the error model is the part that can be argued with.

**Related pages.** [`A2.3`](../A2.3-taylor-aris-dispersion/) is the gallery's
other `S6` page and the other half of this physics: Graetz resolves the radial
profile, Taylor–Aris asks when you may stop resolving it. `A2.1`, `A3.4`.

**Cite the source, not this page:** Graetz, L., *Ueber die
Wärmeleitungsfähigkeit von Flüssigkeiten. Erste Abhandlung*, Annalen der Physik
**254**(1) 79–94 (1882),
[doi:10.1002/andp.18822540106](https://doi.org/10.1002/andp.18822540106).
The catalogue records this as "Graetz (1883)"; the paper is signed *München,
15. Oct. 1882* and the publisher dates the issue 1882. The commonly cited 1883 is
the volume's other dating. A *zweite Abhandlung* followed in Ann. Phys. 261
(1885) and is **not** the paper read here."""))

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
