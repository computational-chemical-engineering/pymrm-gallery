#!/usr/bin/env python3
"""Generate index.ipynb for page A2.8. Run from the page directory.

Quoting convention, copied from A2.5/A2.6: markdown cells are raw triple-DOUBLE-
quoted strings and code cells are raw triple-SINGLE-quoted strings, so a code
cell may contain an ordinary Python docstring. Every one is RAW, so a single
backslash here is a single backslash in the notebook.

House rule this page follows strictly: no number that a cell computes is ever
retyped into a markdown cell. Anything with a computed number in it is emitted
by `display(Markdown(f"..."))` from the cell that computed it.
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# ---------------------------------------------------------------- title -----
cells.append(md(r"""---
title: "Zwietering's micromixing bounds: what a residence-time distribution cannot tell you"
description: "Complete segregation and maximum mixedness bracket the conversion for a given RTD. The bracket is reproduced against all 44 numbers Zwietering printed, and it is wide - a factor of 2.1 at kc0*tau = 50."
categories: [sec:A, struct:S1, struct:S3, struct:S5, tier:T0, data:tier6, phase:liquid]
date: 2026-08-05
---

# Zwietering's micromixing bounds: what a residence-time distribution cannot tell you

**Catalog ID:** `A2.8` · **Structures:** `S1` (batch element), `S3` (1-D steady
BVP), `S5` (convection-dominated, TVD deferred correction) · **Tier:** T0

Four gallery pages already treat the residence-time distribution as the answer.
[`A2.1`](../A2.1-danckwerts-boundary-conditions/) asks what boundary conditions
the axial-dispersion model needs, [`A2.3`](../A2.3-taylor-aris-dispersion/)
derives a dispersion coefficient from first principles,
[`A2.5`](../A2.5-edwards-richardson-dispersion/) measures one and
[`A2.6`](../A2.6-gunn-dispersion-correlations/) correlates it. All four end with
an RTD.

This page is about what you still do not know once you have it. Zwietering's
answer is a **pair of bounds**: for a given RTD and a reaction of order other
than one, the achievable exit concentration lies between the completely
segregated value and the maximum-mixedness value, and nothing about the RTD
narrows it further. The useful question is not whether the bounds are correct -
they are - but **how far apart they are**, because that gap is exactly the part
of a reactor's performance that an RTD measurement cannot reach."""))

# ----------------------------------------------------------- background -----
cells.append(md(r"""## Background

Danckwerts (1953, the paper behind [`A2.1`](../A2.1-danckwerts-boundary-conditions/))
introduced the *degree of segregation* $J$ and showed that two reactors with
**identical** residence-time distributions can give different conversions, because
the RTD says when molecules leave but not who they were mixed with on the way.
For a first-order reaction that does not matter - the rate is linear, so the
average of the rate is the rate of the average. For any other order it matters a
great deal.

Zwietering's contribution is to close the problem. Danckwerts had one extreme:
**complete segregation**, in which a fluid element entering the vessel keeps its
identity until the exit, so the reactor is a batch reactor whose contents are
mixed only in the outlet pipe. That is mixing *as late as possible*. Zwietering
constructs the opposite extreme - mixing *as early as possible*, compatible with
the same RTD - and calls it **maximum mixedness**. Between them, every physically
realisable reactor with that RTD must lie.

The two constructions are the same pipe run in opposite directions
(his Figs. 2 and 3):

- **Complete segregation** — a plug-flow tube with many *side exits*, tapped so
  that the total residence-time distribution is the given $f(t)$. Nothing mixes
  inside; all the mixing happens at the collecting header.
- **Maximum mixedness** — the same tube with the flow reversed, so the side taps
  become *side entrances*. Fresh feed is injected all along the tube, and every
  molecule is mixed on arrival with all the molecules that will leave at the same
  instant. Zwietering's own phrase: the second system "is simply constructed by
  reversal of the flow".

The natural coordinate for the second reactor is not the age of the fluid but its
**life expectation** $\lambda$ — the time still to be spent in the vessel. It runs
from $\infty$ at the entrance to $0$ at the exit, so the governing equation is
integrated *backwards*. That is the numerically interesting part of this page."""))

# ------------------------------------------------- the published model ------
cells.append(md(r"""## The published model

Write $F(t)$ for the residence-time distribution function, $f(t) = dF/dt$ for its
frequency function, $\tau = V/Q$ for the mean residence time, and
$w(t) = 1 - F(t)$ for the fraction of the feed still inside. Zwietering's eq. (8)
gives the internal age frequency function $\phi(\alpha) = w(\alpha)/\tau$, and his
eq. (12) the life-expectation frequency function $\psi(\lambda) = w(\lambda)/\tau$
— **the same function**, which is the formal statement that the two reactors are
one reactor read in two directions.

**Complete segregation.** Each element is a batch reactor. For $r = R(c)$ the
element leaving at time $t$ carries $c(t)$ from $dc/d\theta = -R(c)$, and the exit
is the flow-weighted mixture (his eq. 37)

$$ c_e = \int_0^\infty f(t)\, c(t)\, dt .$$

For $R = kc^2$ the element solution is his eq. (36), $c(t) = c_0/(1 + Kt/\tau)$
with $K = kc_0\tau$ (his eq. 35).

**Maximum mixedness.** A balance on the slice between $\lambda$ and
$\lambda + \Delta\lambda$ of the reversed tube gives his **eq. (31)**

$$ \frac{dc}{d\lambda} \;=\; R(c) \;+\; \frac{f(\lambda)}{1 - F(\lambda)}\,(c - c_0), $$

with $c_e = c(0)$. The three terms are, in order, the reaction, the dilution by
fresh feed entering through the side tubes, and the loss of the material that
leaves. The boundary condition at $\lambda \to \infty$ is $dc/d\lambda = 0$, which
turns the differential equation into an algebraic one there: since
$f/(1-F) \to n/\tau$ for $n$ equal tanks in series,

$$ K\gamma_\infty^{\,p} + n\,(\gamma_\infty - 1) = 0, \qquad \gamma = c/c_0 .$$

Zwietering prints exactly this for $n = 2$ ("$0 = K\gamma_\infty^2 + 2(\gamma_\infty - 1)$").

**The conservative form is the one to discretise.** Multiplying eq. (31) through
by $w$ and collecting,

$$ \frac{d}{d\lambda}\big[\,w(\lambda)\, c\,\big] \;=\; R(c)\,w(\lambda) \;-\; f(\lambda)\,c_0 ,$$

which is a plug-flow reactor with a **varying volumetric flow** $Qw(\lambda)$ and a
distributed side feed. In dimensionless form, with $x = \lambda/\tau$ and
$\gamma = c/c_0$,

$$ \operatorname{div}\big(v\,\gamma\big) \;=\; E(x) \;-\; K\,w(x)\,\gamma^{\,p},
\qquad v = -\,w(x), $$

where $E(x) = \tau f(\tau x)$. The negative velocity is the physics: the flow runs
toward **decreasing** $x$, entering where the life expectancy is large and leaving
at $x = 0$. Discretising the divergence of the flux rather than $v\,d\gamma/dx$
matters here for the same reason it mattered on
[`F2.3`](../F2.3-slurry-bubble-column-ft/): the "velocity" $w$ falls from 1 to 0
across the domain, and writing the convection non-conservatively loses that.

**The two RTDs Zwietering works out.** His eqs. (32) and (33),

$$ f_1(t) = \frac{4t}{\tau^2}e^{-2t/\tau}, \qquad
   f_2(t) = \frac{27}{2}\frac{t^2}{\tau^3}e^{-3t/\tau}, $$

are two and three equal well-mixed vessels in series. He evaluates **four**
states of mixing for each: complete segregation, $n$ real vessels each internally
segregated (mixing at every junction), $n$ real vessels each ideally mixed, and
maximum mixedness. The middle two are genuine intermediate reactors, and they are
what make the bounds testable rather than decorative."""))

# ------------------------------------------ parameters and assumptions ------
cells.append(md(r"""## Parameters and assumptions

There is **one** parameter, $K = kc_0\tau$, and it is dimensionless. Nothing on
this page is fitted; there is no constant to adjust and no property to look up.
That is worth stating plainly, because it fixes what the comparison against
Zwietering's tables is: every one of the 44 printed values is a **held-out test**
of a zero-parameter prediction, not a calibration.

Assumptions, all Zwietering's:

- Steady flow, one entrance and one exit, constant density, so $\tau = V/Q$
  (his eq. 7).
- The RTD is known exactly. Here it is the Erlang distribution of $n$ equal
  stirred tanks, which is analytic, so no measurement error enters.
- Isothermal, single reactant, rate a function of $c$ alone (his eq. 29).
- The reaction order is $p = 2$ throughout the tables; the page also sweeps
  $p$ to locate where the bounds close.

**Numerical parameters.** The maximum-mixedness domain is truncated at
$x = \lambda/\tau = X$ and the far-field root is imposed there; the grid is
graded as $x = X t^2$ so that cells pack toward the exit at $x = 0$, where the
solution turns fastest. Both $X$ and the cell count are refined below, and so is
the tolerance of the independent solver. The production values are printed by the
cell that sets them."""))

# --------------------------------------------------------------- data -------
cells.append(md(r"""## The data

There is **no experimental measurement anywhere in this paper.** It is a theory
paper: the three figures are schematics of pipe-and-tank arrangements with no data
on them, and nothing is digitised here. The data are Zwietering's own tabulated
numbers, read off 300 dpi page renders, and the page is **tier 6** by necessity.

Two files:

- `zwietering-1959-conversions.csv` — the 44 printed values of $c_e/c_0$ in
  Tables 1 and 2 (four mixing states x six values of $K$ x two RTDs, less the four
  dashes where Zwietering did not integrate eq. 31).
- `zwietering-1959-segregation-degree.csv` — the eight values of the degree of
  segregation $J$ in the same tables' last column.

**Resolution and the trap in this scan.** `pdfimages -list` reports CCITT-G4
bilevel images at 300 x 300 ppi native, so 300 dpi is the render resolution and
600 dpi would be interpolation. The text layer is worse than useless: Acrobat 3.0
Capture reads the Greek $\tau$ as the **digit 7** throughout, so eq. (II, 12)'s
$\xi(\lambda) = \tau$ comes back as "4th) = 7". Every equation and every number
here was read from an image, and each numeric cell was cropped and re-read at 2x,
with the four-decimal $J$ column re-read at 4x, because the 1959 typesetting sets
decimals as a raised mid dot and 3 and 8 are the confusable pair in this face.

**The page range is not what everyone cites.** The article is universally cited as
*Chem. Eng. Sci.* **11**, 1-11. The last page carries the printed folio **15**;
the PDF is fifteen pages, the folios run 1 to 15 without a gap, and page 15 holds
eqs. (II, 11) and (II, 12) and the reference list. Cite it as **11(1), 1-15**.

**No other page's dataset is loaded.** `A2.1` publishes no CSV at all
(`datasets: []` in its `meta.yaml`); `A2.3`, `A2.5` and `A2.6` publish dispersion
data that this page never uses, because the RTDs here are analytic Erlang
distributions rather than measured ones. Nothing in this page's numbers therefore
needs reconciling against a row in a sibling page's file. The one finding from a
sibling page that does bear on this one is `A2.6`'s: a quantity read at the last
**cell centre** instead of at the outlet **face** was 11.4 % low while looking
perfectly converged. That trap is live here too and is measured in the break
table below.

**And one thing about the tables themselves is wrong.** It is established from
Zwietering's own printed closed forms, in the first result cell."""))

# ------------------------------------------------------------ env cell ------
cells.append(code(r'''try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml'''))

cells.append(code(r'''import sys, urllib.request
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
from IPython.display import Markdown, display
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import gammaincc, gammaln, expi
from pymrm import (construct_convflux_upwind, construct_div, NumJac, newton,
                   interp_cntr_to_stagg_tvd, vanleer)
from gallery_utils import load_data, load_meta, cite_data, report_agreement

warnings.filterwarnings("ignore", category=UserWarning)
PAGE = "A2.8-zwietering-segregation"
rng = np.random.default_rng(20260805)      # nothing here is stochastic; seeded anyway
plt.rcParams.update({"figure.dpi": 110, "font.size": 9.5})'''))

# ------------------------------------------------ pymrm implementation ------
cells.append(md(r"""## PyMRM implementation

Three of the four mixing states need no PDE at all and are written out in closed
form or as a one-dimensional quadrature - which is the point of ranking the
validation before writing code. Only maximum mixedness needs a solver.

The RTD is the Erlang distribution of $n$ equal tanks. It is evaluated through
`scipy.special.gammaincc` rather than by summing $\sum (nx)^j/j!$, because that
sum overflows by $n \approx 170$ and the plug-flow limit below goes to
$n = 300$."""))

cells.append(code(r'''# ---- residence-time distribution of n equal stirred tanks, mean tau = 1 ----
def E_rtd(x, n):
    """Frequency function f(t) in x = t/tau. Erlang, shape n, scale 1/n."""
    x = np.asarray(x, float)
    if n == 1:
        return np.exp(-x)
    out = np.zeros_like(x)
    m = x > 0
    out[m] = np.exp(n * np.log(n) + (n - 1) * np.log(x[m]) - n * x[m] - gammaln(n))
    return out


def W_rtd(x, n):
    """w = 1 - F(x): fraction of the feed still inside. Regularised upper gamma."""
    return gammaincc(n, n * np.asarray(x, float))


def batch(x, K, p=2.0):
    """Batch element:  dc/dtheta = -k c^p,  gamma(0) = 1,  K = k c0^(p-1) tau."""
    x = np.asarray(x, float)
    if p == 1.0:
        return np.exp(-K * x)
    if p > 1.0:
        return (1.0 + (p - 1.0) * K * x) ** (-1.0 / (p - 1.0))
    # order below one: the element is exhausted in finite time and stays at zero
    return np.clip(1.0 - (1.0 - p) * K * x, 0.0, None) ** (1.0 / (1.0 - p))


# ---- state 1: complete segregation (Zwietering eqs. 36 and 37) -------------
def gamma_seg(K, n, p=2.0):
    f = lambda x: float(E_rtd(np.array([x]), n)[0] * batch(np.array([x]), K, p)[0])
    knots = [0.0, .5 / n, 1. / n, 2. / n, 4. / n, 10. / n, 30. / n]
    if p < 1.0:                      # split at the exhaustion time
        knots = sorted(set(knots + [1.0 / ((1.0 - p) * K)]))
    v = sum(quad(f, a, b, limit=400)[0] for a, b in zip(knots[:-1], knots[1:]))
    return v + quad(f, knots[-1], np.inf, limit=400)[0]


# ---- state 2: n real vessels, each internally segregated -------------------
def gamma_real_seg(K, n):
    """Batch elements inside each tank; the feed to tank i+1 is tank i's cup mix.
    This is exactly the chain Zwietering's eqs. (39) and (40) write out for n=2."""
    g = 1.0
    for _ in range(n):
        Kv = K * g
        g *= quad(lambda u: n * np.exp(-n * u) / (1.0 + Kv * u), 0, np.inf, limit=500)[0]
    return g


# ---- state 3: n real vessels, each ideally mixed (his eq. 41) --------------
def gamma_real_mixed(K, n):
    g = 1.0
    for _ in range(n):               # K g_i^2 = n (g_{i-1} - g_i)
        g = (-n + np.sqrt(n * n + 4 * K * n * g)) / (2 * K)
    return g


# ---- the far-field root, lambda -> infinity -------------------------------
def gamma_far(K, n, p=2.0):
    """K g^p + n (g - 1) = 0 : Zwietering's dc/dlambda = 0 condition."""
    return brentq(lambda g: K * max(g, 0.0) ** p + n * (g - 1.0), 1e-15, 1.0)'''))

cells.append(md(r"""The maximum-mixedness reactor is Zwietering's Fig. 3 written as a
pymrm convection-reaction problem. Everything constant is assembled once, in
`__init__`; the reaction is pointwise so `NumJac((ncell, 1), ...)` gives a
tridiagonal-blocked Jacobian rather than a dense one, and the shape is written
`(ncell, 1)` and never a bare `(ncell,)`.

Two details carry the physics and are worth reading before the code:

- **The face velocity is $v = -w(x)$**, not a constant. It is negative because the
  flow runs from large life expectancy to the exit, and it *decays to zero* at
  large $x$ because almost nothing is still in the vessel that far from leaving.
- **The exit value is read at the face $x = 0$**, through the same TVD
  reconstruction the deferred correction uses - not at the first cell centre. On a
  uniform grid the difference is first order in $h$ and is measured in the break
  table."""))

cells.append(code(r'''class MaxMixedness:
    """Zwietering eq. (31) in conservative form as a 1-D convection-reaction BVP.

        div(v gamma) = E(x) - K w(x) gamma^p ,   v = -w(x)  on the faces,

    x = lambda/tau, gamma = c/c0, w = 1-F. The flow runs toward DECREASING x and
    leaves the domain at the exit x = 0.

    The keyword arguments after `p` exist only to inject defects for the break
    table; the defaults are the model.
    """

    def __init__(self, K, n, X=6.0, ncell=800, grade=2.0, p=2.0, limiter=vanleer,
                 rxn_sign=+1.0, feed=True, v_const=False, far_field="root"):
        self.K, self.n, self.p = K, n, p
        self.limiter, self.rxn_sign, self.feed = limiter, rxn_sign, feed
        self.shape = (ncell, 1)                       # never a bare (ncell,)
        t = np.linspace(0.0, 1.0, ncell + 1)
        self.x_f = X * t ** grade                     # cells packed at the exit
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.g_inf = gamma_far(K, n, p) if far_field == "root" else 1.0
        # Outward normal, so both dicts read  a dgamma/dn + b gamma = d.
        #   x = 0  exit, pure outflow                 -> a=1, b=0, d=0
        #   x = X  entrance, lambda -> infinity       -> a=0, b=1, d=gamma_inf
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0},
                   {"a": 0.0, "b": 1.0, "d": self.g_inf})
        w_f = np.ones_like(self.x_f) if v_const else W_rtd(self.x_f, n)
        self.v = -w_f.reshape(-1, 1)                  # negative: flow toward x=0
        self.w_c = W_rtd(self.x_c, n).reshape(-1, 1)
        self.E_c = E_rtd(self.x_c, n).reshape(-1, 1)
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                  self.bc, v=self.v, axis=0)
        self.div = construct_div(self.shape, self.x_f, nu=0, axis=0)  # nu=0 Cartesian
        self.A = self.div @ conv
        self.b = self.div @ conv_bc
        self.numjac = NumJac(self.shape)              # pointwise source, last axis
        self.corr = np.zeros((ncell, 1))

    def source(self, g):
        """S = E(x) - K w gamma^p : side feed in, reaction out."""
        s = -self.rxn_sign * self.K * self.w_c * np.clip(g, 0.0, None) ** self.p
        return s + self.E_c if self.feed else s

    def residual(self, g):
        s, js = self.numjac(self.source, g)
        return (self.b + self.A @ g.reshape((-1, 1)) - s.reshape((-1, 1))
                + self.corr), self.A - js

    def solve(self, maxfev=80, max_it=120, tol=1e-13, strict=True):
        g = np.full(self.shape, self.g_inf)
        self.corr = np.zeros((self.shape[0], 1))
        res = newton(self.residual, g, maxfev=maxfev)
        g = res.x.reshape(self.shape)
        it, done = 0, self.limiter is None
        if self.limiter is not None:                  # van Leer deferred correction
            for it in range(1, max_it + 1):
                _, dg = interp_cntr_to_stagg_tvd(g, self.x_f, self.x_c, self.bc,
                                                 self.v, tvd_limiter=self.limiter,
                                                 axis=0)
                self.corr = np.asarray(self.div @ (self.v * dg.reshape(-1, 1)))
                res = newton(self.residual, g, maxfev=maxfev)
                gn = res.x.reshape(self.shape)
                done = np.max(np.abs(gn - g)) < tol
                g = gn
                if done:
                    break
        self.ok = bool(res.success and done)
        # A deferred correction that silently returns its cap is how an
        # unconverged number gets published. `strict=False` is used in exactly
        # one place - the break table, where failing to converge is itself a
        # second alarm.
        if strict:
            assert self.ok, "maximum-mixedness solve did not converge"
        self.g, self.n_it = g, it
        return self

    def exit_value(self):
        """gamma AT x = 0 - the reconstructed face value, not the cell centre."""
        up, dg = interp_cntr_to_stagg_tvd(self.g, self.x_f, self.x_c, self.bc,
                                          self.v, tvd_limiter=vanleer, axis=0)
        return float((up + dg).reshape(-1)[0])

    def exit_cell(self):
        return float(self.g[0, 0])


def gamma_mm(K, n, **kw):
    return MaxMixedness(K, n, **kw).solve().exit_value()


def mm_ode(K, n, X=6.0, p=2.0, rtol=1e-12, atol=1e-14, g_start=None,
           backward=True, ratio=None):
    """INDEPENDENT route: eq. (31) in its non-conservative form, integrated by an
    adaptive explicit Runge-Kutta (Dormand-Prince 8(5,3)). No pymrm operator, no
    grid, no flux, no Newton - and no implicit step either, which matters: the
    equation is contracting in the backward direction, so it is NOT stiff there
    and an explicit method is both legitimate and about ten times faster."""
    g0 = gamma_far(K, n, p) if g_start is None else g_start
    rr = ratio if ratio is not None else (
        lambda x: float(E_rtd(np.atleast_1d(x), n)[0] / W_rtd(x, n)))
    rhs = lambda x, y: K * max(y[0], 0.0) ** p + rr(x) * (y[0] - 1.0)
    s = solve_ivp(rhs, (X, 0.0) if backward else (0.0, X), [g0],
                  method="DOP853", rtol=rtol, atol=atol)
    return float(s.y[0, -1]) if s.success else np.nan


NCELL, XDOM = 800, 6.0
_m = MaxMixedness(10, 2, X=XDOM, ncell=NCELL).solve()
print(f"production settings: X = {XDOM}, {NCELL} cells graded as x = X t^2, "
      f"h(first cell) = {_m.x_f[1]:.2e}, h(last) = {np.diff(_m.x_f)[-1]:.3f}")
print(f"K = 10, n = 2 : gamma_inf = {_m.g_inf:.8f} -> gamma(0) = {_m.exit_value():.8f}"
      f"   ({_m.n_it} deferred-correction sweeps)")'''))

# ------------------------------------------------------------ results -------
cells.append(md(r"""## Results

### First, the tables are not what their captions say

Before any comparison can be scored, one thing has to be settled. Zwietering's
Table 1 is captioned *two* well-mixed vessels in series and Table 2 *three*. The
conversion bodies are the other way round, and this is established **from his own
printed closed forms**, with no model of ours entering:

- **Eq. (42)**, $c_e/c_0 = \left[-1 + \sqrt{-1 + 2\sqrt{1+2K}}\right]/K$, is derived
  three inches above the tables from eq. (41), whose $2/\tau$ makes it
  unambiguously the **two**-vessel ideally-mixed chain.
- **Eq. (38)** is derived by substituting eq. (32) - the **two**-vessel frequency
  function - into the segregation average.
- **Eqs. (39) and (40)** are the **two**-vessel segregated chain, with $2/K$ in the
  exponential integral.

All three are arithmetic on printed formulae. The next cell evaluates them and
asks which table's row they land on."""))

cells.append(code(r'''conv = load_data("zwietering-1959-conversions.csv", page=PAGE)
Jdat = load_data("zwietering-1959-segregation-degree.csv", page=PAGE)
meta = load_meta("zwietering-1959-conversions.csv", page=PAGE)
KS = sorted(conv["K"].unique())
ROWS = ["complete_segregation", "real_segregated", "real_ideally_mixed",
        "maximum_mixedness"]

def printed(table, row):
    d = conv[(conv.table == table) & (conv.mixing_state == row)]
    return dict(zip(d.K, d.ce_over_c0))

# Zwietering's own printed closed forms, all for TWO vessels
def eq38(K):  return 2 / K + (4 / K ** 2) * np.exp(2 / K) * expi(-2 / K)
def eq42(K):  return (-1 + np.sqrt(-1 + 2 * np.sqrt(1 + 2 * K))) / K
def eq3940(K):
    c1 = -np.exp(2 / K) / (K / 2) * expi(-2 / K)
    a = 2 / (K * c1)
    return c1 * (-a * np.exp(a) * expi(-a))

closed = {"complete_segregation": ("eq. (38)", eq38),
          "real_segregated":      ("eq. (39)+(40)", eq3940),
          "real_ideally_mixed":   ("eq. (42)", eq42)}

recs = []
for row, (name, fn) in closed.items():
    vals = {K: fn(K) for K in KS}
    d1 = max(abs(vals[K] - v) for K, v in printed(1, row).items())
    d2 = max(abs(vals[K] - v) for K, v in printed(2, row).items())
    recs.append({"printed formula (2 vessels)": name, "row": row,
                 "max |dev| vs Table 1 row": d1, "max |dev| vs Table 2 row": d2,
                 "lands on": "Table 1" if d1 < d2 else "Table 2"})
swap = pd.DataFrame(recs)
display(swap.style.format({"max |dev| vs Table 1 row": "{:.4f}",
                           "max |dev| vs Table 2 row": "{:.4f}"}).hide(axis="index"))

EQ38_T2 = max(abs(eq38(K) - v) for K, v in printed(2, "complete_segregation").items())
EQ38_T1 = max(abs(eq38(K) - v) for K, v in printed(1, "complete_segregation").items())
EQ42_T2 = max(abs(eq42(K) - v) for K, v in printed(2, "real_ideally_mixed").items())
EQ42_T1 = max(abs(eq42(K) - v) for K, v in printed(1, "real_ideally_mixed").items())
EQ3940_T2 = max(abs(eq3940(K) - v) for K, v in printed(2, "real_segregated").items())
EQ3940_T1 = max(abs(eq3940(K) - v) for K, v in printed(1, "real_segregated").items())

display(Markdown(
    f"All three two-vessel closed forms land on **Table 2**, which is captioned "
    f"three vessels. Eq. (42) reproduces Table 2's ideally-mixed row to "
    f"**{EQ42_T2:.4f}** and Table 1's only to {EQ42_T1:.4f}; eq. (38) gives "
    f"{EQ38_T2:.4f} against {EQ38_T1:.4f}; eqs. (39)+(40) give {EQ3940_T2:.4f} "
    f"against {EQ3940_T1:.4f}. The reading is therefore that **the conversion "
    f"bodies of the two tables are interchanged with respect to their captions**. "
    f"Everything below scores against the interchanged assignment and reports the "
    f"caption assignment beside it."))

BODY_N = {1: 3, 2: 2}      # the RTD the conversion body actually belongs to
CAP_N  = {1: 2, 2: 3}      # what the caption says'''))

cells.append(md(r"""The third piece of evidence needs no arithmetic at all. Three tanks
are a *narrower* residence-time distribution than two, and for a second-order
reaction a narrower distribution converts more, so the smaller $c_e/c_0$ belongs
to the three-tank table. Table 1, captioned two tanks, prints the smaller values
throughout.

The fourth is that the $J$ column does **not** move with the body - Table 1's $J$
values are the two-vessel ones and Table 2's the three-vessel ones, matching their
captions. That is why the finding is stated as *the conversion bodies are
interchanged* rather than *the two tables are swapped*: swapping two whole tables
is not an observable event, but a table whose numeric body disagrees with both its
own caption and its own last column is.

A fifth, found in verification, pins caption and $J$ column together from the
paper's own prose. Journal page 12 says: *"The results … are given in Table 1 …
The value $J = 0{\cdot}0275$ is the minimum value for this residence time
distribution"* — said of the **two-vessel** frequency function, and
$J = 0.0275$ *is* the two-vessel maximum-mixedness value
($(2e\,E_1(1)-1)/7 = 0.027528$). So the caption and the $J$ column belong
together, and only the conversion bodies can be the interchanged part.

### All 44 printed conversions, both readings"""))

cells.append(code(r'''MODEL = {"complete_segregation": gamma_seg,
         "real_segregated":      gamma_real_seg,
         "real_ideally_mixed":   gamma_real_mixed,
         "maximum_mixedness":    lambda K, n: gamma_mm(K, n, X=XDOM, ncell=NCELL)}

def score(assign):
    devs, per_row = [], {}
    for t in (1, 2):
        n = assign[t]
        for row in ROWS:
            d = [abs(MODEL[row](K, n) - val) for K, val in printed(t, row).items()]
            per_row[(t, row)] = max(d)
            devs += d
    return np.array(devs), per_row

DEV_FIX, PER_FIX = score(BODY_N)
DEV_CAP, PER_CAP = score(CAP_N)

tab = pd.DataFrame([
    {"reading": "conversion bodies interchanged", "N": DEV_FIX.size,
     "max |dev|": DEV_FIX.max(), "rms |dev|": np.sqrt(np.mean(DEV_FIX ** 2)),
     "within 0.001": int((DEV_FIX <= 0.001).sum())},
    {"reading": "as captioned", "N": DEV_CAP.size,
     "max |dev|": DEV_CAP.max(), "rms |dev|": np.sqrt(np.mean(DEV_CAP ** 2)),
     "within 0.001": int((DEV_CAP <= 0.001).sum())}])
display(tab.style.format({"max |dev|": "{:.4f}", "rms |dev|": "{:.5f}"}).hide(axis="index"))

per = pd.DataFrame([{"table": t, "caption": f"{CAP_N[t]} tanks",
                     "row": r, "max |dev| interchanged": PER_FIX[(t, r)],
                     "max |dev| as captioned": PER_CAP[(t, r)]}
                    for t in (1, 2) for r in ROWS])
display(per.style.format({"max |dev| interchanged": "{:.4f}",
                          "max |dev| as captioned": "{:.4f}"}).hide(axis="index"))

MM_ROW_MAX = max(PER_FIX[(1, "maximum_mixedness")], PER_FIX[(2, "maximum_mixedness")])
ALG_ROW_MAX = max(PER_FIX[(t, r)] for t in (1, 2) for r in ROWS[:3])

# Zwietering's own hand-integration bias, computed rather than asserted: the
# SIGNED relative deviation of his maximum-mixedness values from ours.
bias = {}
for t in (1, 2):
    n = BODY_N[t]
    bias[n] = np.array([(v - MODEL["maximum_mixedness"](K, n))
                        / MODEL["maximum_mixedness"](K, n)
                        for K, v in printed(t, "maximum_mixedness").items()])
MM_BIAS_2 = float(bias[2].mean())
MM_BIAS_3 = float(bias[3].mean())
MM_BIAS_2_ALLPOS = bool((bias[2] > 0).all())

# how many of the 44 round, at three decimals, to exactly what Zwietering printed
ROUND_EXACT = 0
for t in (1, 2):
    n = BODY_N[t]
    for row in ROWS:
        for K, val in printed(t, row).items():
            ROUND_EXACT += int(round(MODEL[row](K, n), 3) == round(val, 3))
display(Markdown(
    f"Under the interchanged reading the worst of the 44 is **{DEV_FIX.max():.4f}** "
    f"and {int((DEV_FIX <= 0.001).sum())} of {DEV_FIX.size} land inside 0.001. Under "
    f"the caption reading the worst is {DEV_CAP.max():.4f} and "
    f"{int((DEV_CAP <= 0.001).sum())} of {DEV_CAP.size} land inside 0.001 - the two "
    f"readings differ by a factor of {DEV_CAP.max()/DEV_FIX.max():.1f} in the worst "
    f"cell and {np.sqrt(np.mean(DEV_CAP**2))/np.sqrt(np.mean(DEV_FIX**2)):.0f} in "
    f"rms.\n\n"
    f"The residual structure is informative. The three rows that are algebra - the "
    f"two closed-form chains and the segregation quadrature - agree to "
    f"**{ALG_ROW_MAX:.4f}**, one unit in his last printed decimal. The "
    f"maximum-mixedness row, the only one he had to integrate by hand, is the "
    f"outlier at {MM_ROW_MAX:.4f}, and all of that sits in the two-tank case: his "
    f"three-tank maximum-mixedness values agree to "
    f"{PER_FIX[(1, 'maximum_mixedness')]:.4f}, which is rounding, while his two-tank "
    f"ones are off by up to {PER_FIX[(2, 'maximum_mixedness')]:.4f}. The two-tank "
    f"deviations are {'all positive' if MM_BIAS_2_ALLPOS else 'mixed in sign'} and "
    f"average **{MM_BIAS_2*100:+.2f} %** against {MM_BIAS_3*100:+.3f} % for the "
    f"three-tank row - a systematic positive bias of the kind a first-order hand "
    f"integration in finite steps produces, and the only place in either table "
    f"where anything sits outside his own rounding. Taken cell by cell, "
    f"**{ROUND_EXACT} of the 44** round at three decimals to exactly the digits he "
    f"printed."))'''))

cells.append(md(r"""### The bracket

The two bounds and the two intermediate reactors, against $K = kc_0\tau$. The
shaded band is the region Zwietering's theorem says every reactor with that RTD
must occupy; the markers are his printed numbers under the interchanged reading;
the dashed lines are the two ideal reactors, which are *not* bounds on this
problem - they belong to different RTDs entirely."""))

cells.append(code(r'''Kgrid = np.logspace(np.log10(0.3), np.log10(300), 26)
band = {}
for n in (2, 3):
    lo = np.array([gamma_seg(K, n) for K in Kgrid])
    hi = np.array([mm_ode(K, n, X=XDOM) for K in Kgrid])
    band[n] = (lo, hi,
               np.array([gamma_real_seg(K, n) for K in Kgrid]),
               np.array([gamma_real_mixed(K, n) for K in Kgrid]))

fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.3))
for ax, n in zip(axes, (2, 3)):
    lo, hi, rs, rm = band[n]
    ax.fill_between(Kgrid, lo, hi, color="tab:blue", alpha=0.18,
                    label="reachable with this RTD")
    ax.plot(Kgrid, lo, color="tab:blue", lw=2.0, label="complete segregation")
    ax.plot(Kgrid, hi, color="tab:red", lw=2.0, label="maximum mixedness")
    ax.plot(Kgrid, rs, color="tab:green", lw=1.3, ls="--",
            label=f"{n} real segregated vessels")
    ax.plot(Kgrid, rm, color="tab:purple", lw=1.3, ls="-.",
            label=f"{n} real ideally mixed vessels")
    ax.plot(Kgrid, 1 / (1 + Kgrid), color="k", lw=1.0, ls=":",
            label="plug flow (a different RTD)")
    ax.plot(Kgrid, (-1 + np.sqrt(1 + 4 * Kgrid)) / (2 * Kgrid), color="grey",
            lw=1.0, ls=":", label="one stirred tank (a different RTD)")
    t = 1 if BODY_N[1] == n else 2
    for row, mk in zip(ROWS, ("o", "s", "^", "D")):
        p = printed(t, row)
        ax.plot(list(p.keys()), list(p.values()), mk, ms=4.5, mfc="none",
                color="k", lw=0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$K = k c_0 \tau$"); ax.set_ylabel(r"$c_e/c_0$")
    ax.set_title(f"{n} equal tanks in series"
                 f"  (printed in Zwietering's Table {t})", fontsize=9.5)
    ax.grid(alpha=0.25, which="both")
axes[0].legend(fontsize=7.2, loc="lower left")
fig.suptitle("Zwietering's bounds. Open markers are the printed values; "
             "the band is what the RTD cannot resolve.", fontsize=10)
fig.tight_layout()
plt.show()
print(cite_data(meta))'''))

# ---------------------------------------------------------- validation ------
cells.append(md(r"""## Validation

The paper contains no measurement, so every check below is either against a
printed number, against an analytical limit, or against a second independent
computation. They are listed in descending order of what they can catch.

### 1. Two independent routes to the maximum-mixedness value

The break table further down tests *sensitivity* - whether the answer moves when
something is broken. It cannot test *correctness*: a defect present in both the
model and its own perturbation is invisible to it. The check that can catch that
is a second computation that shares no code.

`MaxMixedness` discretises the **conservative** form on a graded finite-volume
grid with an upwind flux and a van Leer deferred correction, and solves the whole
profile at once by Newton. `mm_ode` integrates the **non-conservative** form
(eq. 31 as Zwietering writes it) as an initial-value problem with an adaptive
explicit Runge-Kutta (Dormand-Prince 8(5,3)), forming no grid and no flux.
The equation is contracting in the backward direction, so it is not stiff there
and an explicit method is legitimate. They share the two RTD
functions and the far-field root, and nothing else."""))

cells.append(code(r'''recs = []
for n in (2, 3):
    for K in KS:
        a, b = gamma_mm(K, n, X=XDOM, ncell=NCELL), mm_ode(K, n, X=XDOM)
        recs.append({"n": n, "K": K, "pymrm (finite volume)": a,
                     "scipy DOP853 (IVP)": b, "difference": a - b})
two = pd.DataFrame(recs)
MM_TWO_ROUTES = two["difference"].abs().max()
display(two.style.format({"pymrm (finite volume)": "{:.9f}",
                          "scipy DOP853 (IVP)": "{:.9f}",
                          "difference": "{:+.2e}"}).hide(axis="index"))
display(Markdown(
    f"Worst disagreement over the twelve cases: **{MM_TWO_ROUTES:.1e}**, which is "
    f"the discretisation error of the {NCELL}-cell grid and not a modelling "
    f"difference - it falls at second order under refinement, below. That is "
    f"{MM_TWO_ROUTES/MM_ROW_MAX*100:.2f} % of the deviation being measured against "
    f"Zwietering's own maximum-mixedness numbers, so the comparison is limited by "
    f"his arithmetic and not by ours."))'''))

cells.append(md(r"""### 2. Both refinement axes, not one

The maximum-mixedness problem has **two** discretisation knobs, and refining only
one measures the wrong thing - the failure recorded on `H1.9` and `J3.1`, where in
both cases the unmeasured axis carried the larger error. Here they are the cell
count $h$ and the domain truncation $X$ at which the far-field root is imposed.
The independent solver adds a third, its tolerance."""))

cells.append(code(r'''# --- axis 1: grid ----------------------------------------------------------
grid_rows, ORDERS = [], {}
for K, n in ((10, 2), (30, 3)):
    ref = mm_ode(K, n, X=XDOM, rtol=1e-13, atol=1e-15)
    for lim, name in ((vanleer, "van Leer deferred correction"), (None, "bare upwind")):
        errs = []
        for nc in (100, 200, 400, 800):
            e = abs(gamma_mm(K, n, X=XDOM, ncell=nc, limiter=lim) - ref)
            errs.append(e)
            grid_rows.append({"K": K, "n": n, "scheme": name, "cells": nc,
                              "|error|": e})
        ORDERS[(K, n, name)] = float(np.log2(errs[-2] / errs[-1]))
g = pd.DataFrame(grid_rows)
piv = g.pivot_table(index=["K", "n", "scheme"], columns="cells", values="|error|")
piv["observed order"] = [ORDERS[(k, n, s)] for k, n, s in piv.index]
display(piv.style.format({c: "{:.2e}" for c in piv.columns[:-1]}
                         | {"observed order": "{:.2f}"}))

ORD_VL = min(ORDERS[(K, n, "van Leer deferred correction")] for K, n in ((10, 2), (30, 3)))
ORD_UW = min(ORDERS[(K, n, "bare upwind")] for K, n in ((10, 2), (30, 3)))

# --- axis 2: domain truncation --------------------------------------------
K0, N0 = 10, 2
ref_long = mm_ode(K0, N0, X=20.0, rtol=1e-13, atol=1e-15)
Xs = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
trunc = np.array([abs(mm_ode(K0, N0, X=X, rtol=1e-13, atol=1e-15) - ref_long)
                  for X in Xs])
sel = trunc > 1e-13
TRUNC_RATE = float(-np.polyfit(Xs[sel], np.log(trunc[sel]), 1)[0])
TRUNC_AT_PROD = float(trunc[np.argmin(np.abs(Xs - XDOM))]) if XDOM in Xs else \
    float(abs(mm_ode(K0, N0, X=XDOM, rtol=1e-13, atol=1e-15) - ref_long))
display(pd.DataFrame({"X": Xs, "|error| in gamma(0)": trunc})
        .style.format({"X": "{:.1f}", "|error| in gamma(0)": "{:.2e}"}).hide(axis="index"))

# The decay is not merely fitted, it is PREDICTED. A perturbation imposed at X is
# damped on the way to the exit by exp(-int_0^X mu dx) with
# mu = d(rhs)/d(gamma) = 2 K gamma + E/w, and the far-field error being corrected
# itself falls like 1/X, because E/w approaches its limit n algebraically
# (E/w = 4x/(1+2x) = 2 - 1/x + ... for two tanks). So
#     error(X1)/error(X2)  ~  (X2/X1) * exp( int_{X1}^{X2} mu dx ).
m = MaxMixedness(K0, N0, X=XDOM, ncell=NCELL).solve()
mu = 2 * K0 * m.g.ravel() + E_rtd(m.x_c, N0) / W_rtd(m.x_c, N0)
pred_rows = []
for a, b in zip(Xs[1:5], Xs[2:6]):
    s = (m.x_c >= a) & (m.x_c <= b)
    damp = np.exp(np.trapezoid(mu[s], m.x_c[s]))
    obs = trunc[list(Xs).index(a)] / trunc[list(Xs).index(b)]
    pred = damp * (b / a)
    pred_rows.append({"X window": f"{a:.1f} -> {b:.1f}", "observed ratio": obs,
                      "exp(int mu dx)": damp, "predicted ratio": pred,
                      "obs / predicted": obs / pred})
pt = pd.DataFrame(pred_rows)
TRUNC_MODEL_MAX = float(np.abs(pt["obs / predicted"] - 1).max())
display(pt.style.format({"observed ratio": "{:.2f}", "exp(int mu dx)": "{:.2f}",
                         "predicted ratio": "{:.2f}",
                         "obs / predicted": "{:.3f}"}).hide(axis="index"))

# --- axis 3: tolerance of the independent solver ---------------------------
ref_tol = mm_ode(K0, N0, X=XDOM, rtol=1e-13, atol=1e-15)
tolrow = {rt: abs(mm_ode(K0, N0, X=XDOM, rtol=rt, atol=rt * 1e-2) - ref_tol)
          for rt in (1e-6, 1e-8, 1e-10, 1e-12)}
TOL_AT_PROD = tolrow[1e-12]
display(Markdown(
    f"**Grid.** The van Leer deferred correction converges at observed order "
    f"**{ORD_VL:.2f}** and bare upwind at **{ORD_UW:.2f}**, so the correction is "
    f"buying a genuine order and not a constant.\n\n"
    f"**Truncation.** The error decays at **{TRUNC_RATE:.2f} per unit** of $X$ and "
    f"is **{TRUNC_AT_PROD:.1e}** at the production $X = {XDOM:g}$ - at the "
    f"double-precision floor, a factor {MM_TWO_ROUTES/max(TRUNC_AT_PROD,1e-300):.0e} "
    f"below the grid error. The reported numbers are therefore grid-limited and nothing "
    f"else. The decay is also *predicted* rather than merely fitted: a perturbation "
    f"imposed at $X$ is damped by $\\exp(-\\int_0^X \\mu\\,dx)$ with "
    f"$\\mu = 2K\\gamma + E/w$, and the far-field error itself falls like $1/X$ "
    f"because $E/w$ approaches $n$ algebraically. Observed over predicted lands "
    f"within **{TRUNC_MODEL_MAX*100:.0f} %** across the four windows.\n\n"
    f"**Tolerance.** The independent solver at rtol $10^{{-12}}$ is converged to "
    f"{TOL_AT_PROD:.1e}, three orders below the grid error, so it functions as an "
    f"exact reference at this resolution."))'''))

cells.append(md(r"""### 3. The analytical limits, each of which can fail

Four limits, in decreasing order of how much they test.

**(a) The exponential RTD.** For a single stirred tank ($n = 1$) Zwietering's
eq. (II, 12) states $\xi(\lambda) = \tau$: the function is constant, so the
maximum-mixedness solution must be constant and equal to the ordinary stirred-tank
root. This is a *degenerate* case for the solver - the domain over which the
solution is flat is the whole domain - and getting the coefficient $E/w$ wrong by
any factor breaks it immediately.

**(b) First order.** Zwietering's section 9b: for $R = kc$ every state of mixing
gives the same answer, $\int f(t)e^{-kt}dt$, which for $n$ tanks is
$(1 + K/n)^{-n}$. The segregation quadrature and the backward ODE integration are
different code paths and must land on the same number and on the analytical one.

**(c) Plug flow.** As $n \to \infty$ the RTD narrows to a delta function, which
leaves no freedom for a degree of segregation at all - Zwietering's own footnote:
"the plug flow reactor is at the same time completely segregated and in the state
of maximum mixedness". The bracket must close onto $1/(1+K)$.

**(d) Reaction order.** The bracket must *reverse* through order one - maximum
mixedness giving the least conversion above it and the most below it (his section
1, "the minimum conversion for a reaction of an order higher than one and equally
the maximum conversion for an order of reaction smaller than one")."""))

cells.append(code(r'''rows = []
for K in (1, 5, 20, 100):
    mm, root = mm_ode(K, 1, X=25.0), gamma_far(K, 1)
    rows.append({"K": K, "maximum mixedness, n=1": mm, "stirred-tank root": root,
                 "|difference|": abs(mm - root)})
lim_a = pd.DataFrame(rows)
LIM_CSTR = lim_a["|difference|"].max()
display(lim_a.style.format({"maximum mixedness, n=1": "{:.12f}",
                            "stirred-tank root": "{:.12f}",
                            "|difference|": "{:.1e}"}).hide(axis="index"))

rows = []
for n in (2, 3):
    for K in (5, 30):
        exact = (1 + K / n) ** (-n)
        s, mmv = gamma_seg(K, n, p=1.0), mm_ode(K, n, X=25.0, p=1.0)
        rows.append({"n": n, "K": K, "(1+K/n)^-n": exact,
                     "segregation - exact": s - exact, "max mixedness - exact": mmv - exact,
                     "bracket width": abs(mmv - s)})
lim_b = pd.DataFrame(rows)
FIRST_ORDER_WIDTH = lim_b["bracket width"].max()
FIRST_ORDER_VS_EXACT = max(lim_b["segregation - exact"].abs().max(),
                           lim_b["max mixedness - exact"].abs().max())
display(lim_b.style.format({"(1+K/n)^-n": "{:.10f}", "segregation - exact": "{:+.1e}",
                            "max mixedness - exact": "{:+.1e}",
                            "bracket width": "{:.1e}"}).hide(axis="index"))

rows = []
for n in (2, 3, 5, 10, 30, 100, 300):
    s, mmv = gamma_seg(30, n), mm_ode(30, n, X=max(3.0, 30.0 / n))
    rows.append({"tanks n": n, "segregation": s, "max mixedness": mmv,
                 "relative width": (mmv - s) / s,
                 "segregation - plug flow": s - 1 / 31})
lim_c = pd.DataFrame(rows)
PFR_WIDTH_300 = float(lim_c["relative width"].iloc[-1])
PFR_GAP_300 = float(lim_c["segregation - plug flow"].iloc[-1])
display(lim_c.style.format({"segregation": "{:.6f}", "max mixedness": "{:.6f}",
                            "relative width": "{:.4f}",
                            "segregation - plug flow": "{:.2e}"}).hide(axis="index"))

rows = []
for p in (0.5, 0.7, 0.9, 1.0, 1.1, 1.5, 2.0, 3.0):
    s, mmv = gamma_seg(10, 2, p=p), mm_ode(10, 2, X=XDOM, p=p)
    rows.append({"order p": p, "segregation": s, "max mixedness": mmv,
                 "mm - seg": mmv - s,
                 "who converts more": ("neither, to 1e-12" if abs(mmv - s) < 1e-12
                                      else ("maximum mixedness" if mmv < s
                                            else "segregation"))})
lim_d = pd.DataFrame(rows)
ORDER_MARGIN_BELOW = float(lim_d[lim_d["order p"] == 0.9]["mm - seg"].iloc[0])
ORDER_MARGIN_ABOVE = float(lim_d[lim_d["order p"] == 1.1]["mm - seg"].iloc[0])
display(lim_d.style.format({"order p": "{:.1f}", "segregation": "{:.7f}",
                            "max mixedness": "{:.7f}",
                            "mm - seg": "{:+.3e}"}).hide(axis="index"))

display(Markdown(
    f"**(a)** The exponential-RTD collapse holds to **{LIM_CSTR:.1e}** over two "
    f"decades of $K$.\n\n"
    f"**(b)** For a first-order reaction the two bounds coincide to "
    f"**{FIRST_ORDER_WIDTH:.1e}** and both sit on $(1+K/n)^{{-n}}$ to "
    f"{FIRST_ORDER_VS_EXACT:.1e}. This is a real check - the quadrature and the "
    f"backward integration share no code - but note it is **below "
    f"`check_agreement.py`'s `ABS_FLOOR` of 1e-12 and is therefore not compared by "
    f"CI at all**. It is reported here, and the break table shows what it catches.\n\n"
    f"**(c)** The bracket closes onto plug flow: at 300 tanks its relative width is "
    f"**{PFR_WIDTH_300*100:.1f} %**, down from {float(lim_c['relative width'].iloc[0])*100:.0f} % "
    f"at two, and the segregated value sits {PFR_GAP_300:.1e} above $1/(1+K)$.\n\n"
    f"**(d)** The bracket reverses through order one: at $p = 1.1$ maximum "
    f"mixedness lies {ORDER_MARGIN_ABOVE:+.1e} above segregation and at $p = 0.9$ it "
    f"lies {ORDER_MARGIN_BELOW:+.1e} below it, so which bound is the upper one is "
    f"decided by the order and not by the code. At $p = 0.5$ maximum mixedness "
    f"drives the exit concentration to "
    f"zero, to within {abs(float(lim_d[lim_d['order p']==0.5]['max mixedness'].iloc[0])):.0e} - complete "
    f"conversion, which for an order below one is reached in finite time - while "
    f"complete segregation leaves "
    f"{float(lim_d[lim_d['order p']==0.5]['segregation'].iloc[0]):.4f} unreacted."))'''))

cells.append(md(r"""### 4. The degree of segregation - a quantity the reactor solve never touches

$J = \operatorname{var}\alpha_P/\operatorname{var}\alpha$ depends on the RTD and the
state of mixing alone; no reaction enters it. Reproducing the eight printed values
therefore exercises the age and life-expectation machinery of sections 3-7, which
the conversion columns do not.

Four of the eight have exact closed forms. For $n$ equal tanks, the mean age in
tank $i$ is $i\tau/n$, so $n$ ideally mixed vessels give
$\operatorname{var}\alpha_P = \tau^2(n^2-1)/12n^2$ and, against
$\operatorname{var}\alpha = \tau^2(n+1)(n+5)/12n^2$ from his eq. (11),

$$ J_{\text{ideally mixed}} = \frac{n-1}{n+5} . $$

For $n$ internally segregated vessels a "point" carries a common sojourn in its
current tank but ages accumulated independently over the $i-1$ earlier ones, so
the variance *within* a point is $\tau^2(n-1)/2n^2$ and
$J = 1 - 6(n-1)/\big((n+1)(n+5)\big)$. The maximum-mixedness entries are
quadratures of his eq. (28)."""))

cells.append(code(r'''def var_alpha(n):     return (n + 1) * (n + 5) / (12 * n * n)
def J_real_mixed(n):  return (n - 1) / (n + 5)
def J_real_seg(n):    return 1 - 6 * (n - 1) / ((n + 1) * (n + 5))

def alpha_P(lam, n, U=60.0):
    """Mean age of the points with life expectation lam, Zwietering eq. (27)."""
    return quad(lambda s: W_rtd(s, n), lam, U, limit=500)[0] / W_rtd(lam, n)

def J_mm(n, abar=None, U=60.0):
    """Zwietering eq. (28)."""
    ab = (n + 1) / (2 * n) if abar is None else abar
    v = quad(lambda l: (alpha_P(l, n, U) - ab) ** 2 * W_rtd(l, n), 0, 30, limit=600)[0]
    return v / var_alpha(n)

def var_within_mm(n, U=60.0):
    """Variance of the ages WITHIN a point, from phi_P(a) = f(lam+a)/(1-F(lam))."""
    def inner(lam):
        m2 = quad(lambda a: a * a * E_rtd(np.array([lam + a]), n)[0] / W_rtd(lam, n),
                  0, U, limit=500)[0]
        return (m2 - alpha_P(lam, n, U) ** 2) * W_rtd(lam, n)
    return quad(inner, 0, 30, limit=400)[0]

rows = []
for t in (1, 2):
    n = CAP_N[t]                       # the J column DOES follow the caption
    got = {"complete_segregation": 1.0, "real_segregated": J_real_seg(n),
           "real_ideally_mixed": J_real_mixed(n), "maximum_mixedness": J_mm(n)}
    for r in ROWS:
        p = float(Jdat[(Jdat.table == t) & (Jdat.mixing_state == r)]["J"].iloc[0])
        rows.append({"table": t, "n (caption)": n, "row": r, "computed": got[r],
                     "printed": p, "deviation": got[r] - p})
Jtab = pd.DataFrame(rows)
J_MAX = Jtab["deviation"].abs().max()
J_MAX_EX = Jtab[~((Jtab.table == 2) & (Jtab.row == "maximum_mixedness"))]["deviation"].abs().max()
J_MM3 = float(Jtab[(Jtab.table == 2) & (Jtab.row == "maximum_mixedness")]["computed"].iloc[0])
J_MM3_P = float(Jtab[(Jtab.table == 2) & (Jtab.row == "maximum_mixedness")]["printed"].iloc[0])
display(Jtab.style.format({"computed": "{:.6f}", "printed": "{:.4f}",
                           "deviation": "{:+.5f}"}).hide(axis="index"))

VAR_RESID = max(abs(J_mm(n) * var_alpha(n) + var_within_mm(n) - var_alpha(n))
                for n in (2, 3))
display(Markdown(
    f"Seven of the eight reproduce: worst **{J_MAX_EX:.1e}**, and four of them as "
    f"exact rationals - 5/7, 1/7, 5/8 and 1/4 land on the printed 0.7143, 0.1429, "
    f"0.6250 and 0.2500. Each lands on the table whose **caption** matches, which "
    f"is what makes the interchange of the conversion bodies an interchange rather "
    f"than a relabelling.\n\n"
    f"The eighth does not. The three-tank maximum-mixedness entry computes to "
    f"**{J_MM3:.6f}** against the printed **{J_MM3_P:.4f}**, a "
    f"{abs(J_MM3/J_MM3_P-1)*100:.1f} % discrepancy, while the two-tank one "
    f"reproduces to four decimals. We cannot explain it and do not repair it; it is "
    f"reported as an unresolved deviation. The quadrature itself is not in doubt: "
    f"Zwietering's own Appendix I identity, var alpha = var between + var within, "
    f"closes to {VAR_RESID:.1e}. That identity is **structural** - his appendix "
    f"proves it must hold - and it has only partial power, since the two sides "
    f"share alpha_P and it therefore tests the *mean* consistency of alpha_P but "
    f"not its shape. The break table measures exactly how much it does catch."))'''))

cells.append(md(r"""### 5. Where the bounds bind, and how tight they are

A bound that is never approached is decoration. Two questions decide whether these
are worth anything: how wide is the bracket, and do real reactors actually spread
across it? The second is answerable because Zwietering supplies two *intermediate*
reactors with exactly the same RTD, so their position inside the bracket can be
measured rather than assumed. Define the micromixing position

$$ \mu \;=\; \frac{c_e - c_e^{\text{seg}}}{c_e^{\text{mm}} - c_e^{\text{seg}}}
\;\in\; [0,1] . $$"""))

cells.append(code(r'''rows = []
for n in (2, 3):
    for K in KS:
        s, hi = gamma_seg(K, n), mm_ode(K, n, X=XDOM)
        rs, rm = gamma_real_seg(K, n), gamma_real_mixed(K, n)
        rows.append({"n": n, "K": K, "segregation": s, "max mixedness": hi,
                     "relative width": (hi - s) / s,
                     "ratio hi/lo": hi / s,
                     "mu (real segregated)": (rs - s) / (hi - s),
                     "mu (real ideally mixed)": (rm - s) / (hi - s)})
bind = pd.DataFrame(rows)
display(bind.style.format({"segregation": "{:.5f}", "max mixedness": "{:.5f}",
                           "relative width": "{:.3f}", "ratio hi/lo": "{:.2f}",
                           "mu (real segregated)": "{:.3f}",
                           "mu (real ideally mixed)": "{:.3f}"}).hide(axis="index"))

mus = np.concatenate([bind["mu (real segregated)"].values,
                      bind["mu (real ideally mixed)"].values])
BRACKET_MIN_MARGIN = float(min(mus.min(), 1 - mus.max()))
REL_WIDTH_MAX = float(bind["relative width"].max())
RATIO_MAX = float(bind["ratio hi/lo"].max())
REL_WIDTH_K30_N2 = float(bind[(bind.n == 2) & (bind.K == 30)]["relative width"].iloc[0])
MU_SPAN = float(mus.max() - mus.min())

fig, axes = plt.subplots(1, 2, figsize=(11.4, 3.9))
ax = axes[0]
Kw = np.logspace(np.log10(0.3), np.log10(300), 16)
for n, c in ((2, "tab:blue"), (3, "tab:red"), (5, "tab:green"), (10, "tab:purple")):
    w = [(mm_ode(K, n, X=max(3.0, 30.0 / n)) - gamma_seg(K, n)) / gamma_seg(K, n)
         for K in Kw]
    ax.plot(Kw, w, color=c, lw=1.8, label=f"n = {n} tanks")
ax.set_xscale("log"); ax.set_xlabel(r"$K = k c_0 \tau$")
ax.set_ylabel("bracket width / segregated value")
ax.set_title("The bracket does not close at large K", fontsize=9.5)
ax.grid(alpha=0.25); ax.legend(fontsize=8)

ax = axes[1]
ps = np.array([0.6, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.3, 1.6, 2.0, 2.5, 3.0])
for n, c in ((2, "tab:blue"), (3, "tab:red")):
    d = [mm_ode(10, n, X=XDOM, p=pp) - gamma_seg(10, n, p=pp) for pp in ps]
    ax.plot(ps, d, "o-", color=c, ms=3.5, lw=1.6, label=f"n = {n} tanks")
ax.axhline(0.0, color="k", lw=0.8); ax.axvline(1.0, color="k", lw=0.8, ls=":")
ax.set_xlabel("reaction order p"); ax.set_ylabel(r"$c_e^{\rm mm} - c_e^{\rm seg}$")
ax.set_title("and it vanishes only at order one (K = 10)", fontsize=9.5)
ax.grid(alpha=0.25); ax.legend(fontsize=8)
fig.tight_layout(); plt.show()

display(Markdown(
    f"**The bounds bind hard, and they bind harder as the reaction gets faster.** "
    f"The widest of the twelve printed cases is the two-tank RTD at $K = 50$, where "
    f"the maximum-mixedness exit concentration is **{RATIO_MAX:.2f} times** the "
    f"segregated one - a relative bracket of {REL_WIDTH_MAX*100:.0f} %. At $K = 30$ "
    f"with two tanks it is {REL_WIDTH_K30_N2*100:.0f} %. In conversion terms that is "
    f"the difference between "
    f"{100*(1-float(bind[(bind.n==2)&(bind.K==50)]['segregation'].iloc[0])):.1f} % and "
    f"{100*(1-float(bind[(bind.n==2)&(bind.K==50)]['max mixedness'].iloc[0])):.1f} % "
    f"conversion for one measured RTD.\n\n"
    f"**And real reactors do occupy the bracket.** Zwietering's two intermediate "
    f"states run from mu = {mus.min():.2f} to mu = {mus.max():.2f} across the twelve "
    f"cases, spanning {MU_SPAN*100:.0f} % of the band, and the closest either comes "
    f"to a bound is **{BRACKET_MIN_MARGIN:.2f}** of the bracket width. So neither "
    f"bound is slack and neither is nearly attained: the band is genuinely the "
    f"reachable set, not a loose envelope. That is the number that makes the "
    f"bracketing check able to fail - if the maximum-mixedness solve were wrong by "
    f"more than {BRACKET_MIN_MARGIN*100:.0f} % of the bracket in the wrong "
    f"direction, an intermediate reactor would fall outside it."))'''))

cells.append(md(r"""### 6. Break table - what each reported number is sensitive to

Every metric on this page needs something that moves it. Defects are injected into
the maximum-mixedness solver at $K = 10$, $n = 2$ and the exit value is re-read;
the last column records whether the intermediate reactor with the same RTD still
lies inside the resulting bracket, which is the check the bounds themselves
provide."""))

cells.append(code(r'''K0, N0 = 10, 2
REF = mm_ode(K0, N0, X=XDOM, rtol=1e-13, atol=1e-15)
LO, RM = gamma_seg(K0, N0), gamma_real_mixed(K0, N0)
base = MaxMixedness(K0, N0, X=XDOM, ncell=NCELL).solve()

def brk(label, **kw):
    args = dict(X=XDOM, ncell=NCELL)
    args.update(kw)                      # so X= and ncell= can be overridden
    try:
        m = MaxMixedness(K0, N0, **args).solve(strict=False)
        v, ok = m.exit_value(), m.ok
    except Exception as exc:
        return {"injected defect": label, "gamma(0)": np.nan,
                "shift from reference": np.nan, "solver converged": False,
                "bracket still contains the real reactor": False,
                "note": type(exc).__name__}
    return {"injected defect": label, "gamma(0)": v, "shift from reference": v - REF,
            "solver converged": ok,
            "bracket still contains the real reactor": bool(LO <= RM <= v),
            "note": ""}

rows = [
    {"injected defect": "none (the model)", "gamma(0)": base.exit_value(),
     "shift from reference": base.exit_value() - REF, "solver converged": base.ok,
     "bracket still contains the real reactor": bool(LO <= RM <= base.exit_value()),
     "note": f"{base.n_it} sweeps"},
    brk("reaction sign flipped", rxn_sign=-1.0),
    brk("side feed E(x) deleted from the source", feed=False),
    brk("face velocity -1 instead of -(1-F)", v_const=True),
    brk("first-order reaction (p=1) against 2nd-order tables", p=1.0),
    brk("domain truncated at X = 0.5", X=0.5),
    brk("domain truncated at X = 1.0", X=1.0),
    brk("far-field value 1 instead of the root", far_field="one"),
    brk("far-field value 1 AND X = 0.5", far_field="one", X=0.5),
    brk("coarse grid, 25 cells", ncell=25),
    brk("uniform grid (grade = 1), 800 cells", grade=1.0),
    brk("bare upwind, no deferred correction", limiter=None),
]
uni = MaxMixedness(K0, N0, X=XDOM, ncell=NCELL, grade=1.0).solve()
rows.append({"injected defect": "exit read at the first cell centre (uniform grid)",
             "gamma(0)": uni.exit_cell(), "shift from reference": uni.exit_cell() - REF,
             "solver converged": uni.ok,
             "bracket still contains the real reactor": bool(LO <= RM <= uni.exit_cell()),
             "note": "A2.6's trap"})
rows.append({"injected defect": "exit read at the first cell centre (graded grid)",
             "gamma(0)": base.exit_cell(), "shift from reference": base.exit_cell() - REF,
             "solver converged": base.ok,
             "bracket still contains the real reactor": bool(LO <= RM <= base.exit_cell()),
             "note": "grading defuses it"})
fwd = mm_ode(K0, N0, X=XDOM, backward=False, g_start=REF)
rows.append({"injected defect": "eq. (31) integrated FORWARD from x = 0",
             "gamma(0)": fwd, "shift from reference": fwd - REF,
             "solver converged": bool(np.isfinite(fwd)),
             "bracket still contains the real reactor": False,
             "note": "the equation is unstable in this direction"})
cst = mm_ode(K0, N0, X=XDOM, ratio=lambda x: float(N0))
rows.append({"injected defect": "E/w replaced by its far-field constant n",
             "gamma(0)": cst, "shift from reference": cst - REF,
             "solver converged": True,
             "bracket still contains the real reactor": bool(LO <= RM <= cst),
             "note": "returns gamma_inf exactly"})

bt = pd.DataFrame(rows)
display(bt.style.format({"gamma(0)": "{:.6f}", "shift from reference": "{:+.2e}"})
        .hide(axis="index"))

PHYSICS_DEFECTS = ["reaction sign flipped",
                   "side feed E(x) deleted from the source",
                   "face velocity -1 instead of -(1-F)",
                   "first-order reaction (p=1) against 2nd-order tables",
                   "E/w replaced by its far-field constant n"]
DISCRETISATION_DEFECTS = ["domain truncated at X = 0.5", "domain truncated at X = 1.0",
                          "coarse grid, 25 cells", "bare upwind, no deferred correction",
                          "uniform grid (grade = 1), 800 cells"]
BREAK_MIN_REAL = float(bt[bt["injected defect"].isin(PHYSICS_DEFECTS)]
                       ["shift from reference"].abs().min())
BREAK_MIN_DISC = float(bt[bt["injected defect"].isin(DISCRETISATION_DEFECTS)]
                       ["shift from reference"].abs().min())
FARFIELD_SHIFT = float(bt[bt["injected defect"] ==
                          "far-field value 1 instead of the root"]["shift from reference"].iloc[0])
FARFIELD_SHIFT_SHORT = float(bt[bt["injected defect"] ==
                                "far-field value 1 AND X = 0.5"]["shift from reference"].iloc[0])
CELL_UNIFORM = float(bt[bt["injected defect"].str.contains("uniform grid\\)")]
                       ["shift from reference"].iloc[0])
CELL_GRADED = float(bt[bt["injected defect"].str.contains("graded grid\\)")]
                      ["shift from reference"].iloc[0])

# does the first-order collapse - the check CI cannot see - actually catch things?
fo_ref = (1 + K0 / N0) ** (-N0)
fo_rows = [
    ("none", abs(mm_ode(K0, N0, X=25.0, p=1.0) - gamma_seg(K0, N0, p=1.0))),
    ("E/w replaced by the constant n",
     abs(mm_ode(K0, N0, X=25.0, p=1.0, ratio=lambda x: float(N0)) - gamma_seg(K0, N0, p=1.0))),
    ("E/w for 3 tanks used with a 2-tank RTD",
     abs(mm_ode(K0, N0, X=25.0, p=1.0, ratio=lambda x: E_rtd(x, 3) / W_rtd(x, 3))
         - gamma_seg(K0, N0, p=1.0))),
]
fo = pd.DataFrame(fo_rows, columns=["injected defect", "|mm - segregation|"])
display(fo.style.format({"|mm - segregation|": "{:.2e}"}).hide(axis="index"))
FO_SENS = float(fo["|mm - segregation|"].iloc[1:].min())

display(Markdown(
    f"**What moves.** Every one of the five physics defects moves the answer by at "
    f"least **{BREAK_MIN_REAL:.1e}** and every one of the five discretisation "
    f"defects by at least {BREAK_MIN_DISC:.1e}, against a reported precision of "
    f"{MM_TWO_ROUTES:.1e} - two to eight orders of margin. Deleting the side feed "
    f"or dropping to first order pushes the answer *below* the segregation bound, "
    f"so the bracketing check catches both on its own. Integrating eq. (31) forward "
    f"instead of backward does not merely give a wrong answer, it fails to "
    f"integrate at all: the equation is unstable in that direction, which is the "
    f"reason Zwietering's construction starts at large lambda.\n\n"
    f"**The bracketing check has power in one direction only**, and the table says "
    f"so: the flipped reaction sign sends the answer to {bt[bt['injected defect']=='reaction sign flipped']['gamma(0)'].iloc[0]:.0f} "
    f"and the interior reactor is still 'inside' the resulting bracket, because an "
    f"upper bound that is too high cannot be caught by something it bounds. What "
    f"catches that defect is the solver refusing to converge, and the comparison "
    f"against the printed tables. The bracket catches errors in the other "
    f"direction, which is where it is used.\n\n"
    f"**What does not move, and must therefore not be claimed.** Replacing the "
    f"far-field root by gamma = 1 shifts the answer by **{abs(FARFIELD_SHIFT):.1e}** "
    f"- nothing. The backward integration is strongly contracting, so by $x = "
    f"{XDOM:g}$ every trace of the starting value is gone. This is Zwietering's own "
    f"remark ('the solution is very nearly independent of the starting point') "
    f"turned into a measurement, and it means **no check on this page has any power "
    f"over the far-field condition** at the production domain length. Shorten the "
    f"domain to X = 0.5 and the same defect moves the answer by "
    f"{abs(FARFIELD_SHIFT_SHORT):.1e}, which is what shows the insensitivity is a "
    f"property of the domain length and not of the check.\n\n"
    f"**The A2.6 trap is real here and was defused by the grid, not by luck.** On a "
    f"uniform 800-cell grid, reading the exit at the first cell centre instead of "
    f"the face costs {abs(CELL_UNIFORM):.1e} - {abs(CELL_UNIFORM)/MM_TWO_ROUTES:.0f} "
    f"times the reported precision, and it is first order in h so refining looks "
    f"like convergence. On the graded grid used here the first cell is "
    f"{base.x_f[1]:.1e} wide and the same mistake costs {abs(CELL_GRADED):.1e}. The "
    f"page reads the face anyway.\n\n"
    f"**The first-order collapse is below CI's floor but it is not powerless.** "
    f"Corrupting the coefficient E/w lifts it from {fo['|mm - segregation|'].iloc[0]:.1e} "
    f"to at least **{FO_SENS:.1e}** - eleven orders. So it is a real check on the "
    f"RTD machinery; it is simply one `check_agreement.py` will never compare."))'''))

# ------------------------------------------------------ what pymrm adds -----
cells.append(md(r"""## What pymrm adds

**To the physics, nothing.** Zwietering's bounds are right, his four mixing states
are correctly formulated, and his printed conversions are reproduced to the
precision measured in the Results section above. This page adds no term to
eq. (31) and no state of mixing he did not define.

What it adds is arithmetic he had to do by hand, and three things that follow from
being able to do it cheaply.

1. **The two dashes in each maximum-mixedness row are filled in.** Zwietering
   integrated eq. (31) only for $K = 5, 10, 20, 30$; at $K = 3$ and $K = 50$ the
   tables carry an em dash. Those four values are computed below. They are an
   extension, not agreement, and are labelled as such.
2. **The size of his own integration error is measurable.** Only one row in
   either table sits outside his own three-decimal rounding, and it is the one
   row he had to integrate by hand. The bias is computed above, not asserted.
3. **The interchange of the two tables' conversion bodies** is visible only once
   every value can be recomputed in a second.

**And the conservative discretisation is the pymrm-shaped part.** Eq. (31) is
usually presented, and usually solved, as an initial-value problem. Written as
$\operatorname{div}(v\gamma) = E - Kw\gamma^p$ with $v = -w$ it is an ordinary
convection-reaction boundary-value problem with a decaying velocity field and a
distributed feed - the same operator stack as a plug-flow reactor with side
injection, `construct_convflux_upwind` plus `construct_div` plus a pointwise
`NumJac` source. That reframing is what makes the life-expectation coordinate
behave like any other spatial axis, and it is why the van Leer deferred correction
applies unchanged and delivers second order."""))

cells.append(code(r'''fill = []
for t in (1, 2):
    n = BODY_N[t]
    for K in (3, 50):
        fill.append({"Zwietering's table": t, "caption says": f"{CAP_N[t]} tanks",
                     "RTD the body is": f"{n} tanks", "K": K,
                     "segregation (printed)": printed(t, "complete_segregation")[K],
                     "maximum mixedness (NOT printed)": gamma_mm(K, n, X=XDOM, ncell=NCELL),
                     "bracket width / segregated": None})
fdf = pd.DataFrame(fill)
fdf["bracket width / segregated"] = ((fdf["maximum mixedness (NOT printed)"]
                                      - fdf["segregation (printed)"])
                                     / fdf["segregation (printed)"])
display(fdf.style.format({"segregation (printed)": "{:.3f}",
                          "maximum mixedness (NOT printed)": "{:.4f}",
                          "bracket width / segregated": "{:.3f}"}).hide(axis="index"))
FILLED_K50_N2 = float(fdf[(fdf["RTD the body is"] == "2 tanks") & (fdf.K == 50)]
                      ["maximum mixedness (NOT printed)"].iloc[0])
display(Markdown(
    f"The em dash at $K = 50$ for the two-tank RTD is **{FILLED_K50_N2:.4f}**, "
    f"against the {printed(2, 'complete_segregation')[50]:.3f} printed for complete "
    f"segregation in the same column - the "
    f"widest bracket in either table, and the one Zwietering stopped short of."))'''))

# ----------------------------------------------------------------- reuse ----
cells.append(md(r"""## Reuse

**When this page's result matters.** Whenever a conversion is being predicted from
a measured RTD and the reaction is not first order. The bracket is the honest error
bar on that prediction, and the table in validation section 5 shows it is not
small - it exceeds the segregated value itself at the fast end of Zwietering's own
range, and it *widens* with $K$ rather than closing. A tracer experiment is not a
substitute for knowing the mixing state.

**How to reuse the solver.** `MaxMixedness` takes any $w(x) = 1 - F(x)$ and
$E(x)$; the Erlang pair here is only the example Zwietering worked. Substitute the
RTD you measured (or the one `A2.1`'s axial-dispersion model implies, or the one
`A2.3` predicts for a capillary) and the two bounds follow with no other change.
Three things to carry over:

- **Integrate backwards.** Forward integration of eq. (31) does not converge; the
  break table shows the adaptive solver failing outright. In the conservative
  formulation this is automatic - it is just the sign of the face velocity.
- **Read the exit at the face.** On a uniform grid the first cell centre is first
  order wrong and looks converged, which is `A2.6`'s finding in a different
  reactor.
- **Do not bother tuning the far-field boundary value, and do not claim it as a
  check.** Measured here, it makes no difference at all past $x \approx 2$; what
  *does* matter is that the domain is long enough, and that is the knob to refine.

**What this page cannot conclude.**

- **Nothing here is validated against experiment.** There is no measurement in
  Zwietering 1959. Every number is a comparison against his arithmetic or against
  an analytical limit; the page is tier 6 and calling it anything else would be
  wrong.
- **The bounds are not shown to be attainable by a real vessel.** Zwietering's two
  imaginary reactors are constructions. What the page shows is that two *ordinary*
  reactors with the same RTD sit well inside the band, not that anything reaches
  its edges. Whether a physical vessel can approach maximum mixedness is not
  addressed here or in the paper.
- **The three-tank degree of segregation under maximum mixedness is unresolved.**
  The computed and printed values are in validation section 4 and they disagree by
  more than three per cent. The two-tank companion reproduces exactly and
  Appendix I's variance identity closes, so the quadrature is probably right and
  the printed value probably a slip - but "probably" is as far as the evidence
  goes, and the page does not adjust either number to make them agree.
- **The interchange is an inference about typesetting, not about physics.** The
  evidence is that three printed two-vessel closed forms land on the table
  captioned three vessels, that the physical ordering agrees, and that the $J$
  column does not move. What actually happened in 1959 is not recoverable.
- **The far-field boundary condition is untested.** Not weakly tested - untested.
  See the break table.
- **`A2.4` (tanks in series) is not built**, so the RTDs used here have no page of
  their own to cross-check against. They are analytic and were verified in place:
  each integrates to one and has mean one."""))

cells.append(code(r'''metrics = {
    # --- the printed tables -------------------------------------------------
    "tables_max_dev_interchanged":     float(DEV_FIX.max()),
    "tables_rms_dev_interchanged":     float(np.sqrt(np.mean(DEV_FIX ** 2))),
    "tables_max_dev_as_captioned":     float(DEV_CAP.max()),
    "tables_algebraic_rows_max_dev":   float(ALG_ROW_MAX),
    "tables_max_mixedness_row_max_dev": float(MM_ROW_MAX),
    # --- the interchange, from Zwietering's own printed closed forms --------
    "eq38_vs_table2_max_dev":          float(EQ38_T2),
    "eq38_vs_table1_max_dev":          float(EQ38_T1),
    "eq42_vs_table2_max_dev":          float(EQ42_T2),
    "eq42_vs_table1_max_dev":          float(EQ42_T1),
    "eq3940_vs_table2_max_dev":        float(EQ3940_T2),
    "eq3940_vs_table1_max_dev":        float(EQ3940_T1),
    # --- two independent routes --------------------------------------------
    "mm_pymrm_vs_dop853_max":          float(MM_TWO_ROUTES),
    # --- refinement, both axes plus the reference solver's tolerance --------
    "mm_grid_order_vanleer":           float(ORD_VL),
    "mm_grid_order_upwind":            float(ORD_UW),
    "mm_truncation_decay_rate":        float(TRUNC_RATE),
    "mm_truncation_error_at_X6":       float(TRUNC_AT_PROD),
    "mm_truncation_model_max_rel_err": float(TRUNC_MODEL_MAX),
    "mm_reference_tolerance_error":    float(TOL_AT_PROD),
    # --- analytical limits --------------------------------------------------
    "limit_exponential_rtd_max_dev":   float(LIM_CSTR),
    "limit_first_order_bracket_width": float(FIRST_ORDER_WIDTH),
    "limit_first_order_vs_exact":      float(FIRST_ORDER_VS_EXACT),
    "limit_plugflow_rel_width_n300":   float(PFR_WIDTH_300),
    "limit_plugflow_seg_gap_n300":     float(PFR_GAP_300),
    "order_margin_p_1p1":              float(ORDER_MARGIN_ABOVE),
    "order_margin_p_0p9":              float(ORDER_MARGIN_BELOW),
    # --- degree of segregation ---------------------------------------------
    "J_max_dev_seven_of_eight":        float(J_MAX_EX),
    "J_max_dev_all_eight":             float(J_MAX),
    "J_mm_3tank_computed":             float(J_MM3),
    "J_variance_identity_residual":    float(VAR_RESID),
    # --- how tight the bounds are ------------------------------------------
    "bracket_max_ratio_hi_over_lo":    float(RATIO_MAX),
    "bracket_max_relative_width":      float(REL_WIDTH_MAX),
    "bracket_rel_width_K30_n2":        float(REL_WIDTH_K30_N2),
    "bracket_min_margin_of_interior":  float(BRACKET_MIN_MARGIN),
    # --- break-table sensitivities -----------------------------------------
    "break_min_physics_defect_shift":  float(BREAK_MIN_REAL),
    "break_farfield_root_shift":       float(abs(FARFIELD_SHIFT)),
    "break_farfield_root_shift_X0p5":  float(abs(FARFIELD_SHIFT_SHORT)),
    "break_exit_cell_uniform_grid":    float(abs(CELL_UNIFORM)),
    "break_first_order_collapse_sens": float(FO_SENS),
    # --- the extension ------------------------------------------------------
    "filled_dash_mm_K50_2tank":        float(FILLED_K50_N2),
}
report_agreement("A2.8", metrics)

FLOOR = 1e-12
below = {k: v for k, v in metrics.items() if abs(v) < FLOOR}
display(Markdown(
    "**Metrics below `check_agreement.py`'s `ABS_FLOOR` of 1e-12, which CI does "
    "not compare at all:** " +
    (", ".join(f"`{k}` ({v:.1e})" for k, v in below.items()) if below else "none") +
    ". These are reported for the reader, not protected by the regression suite; "
    "the break table above shows what each of them can and cannot catch."))'''))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                             "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
Path("index.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print(f"wrote index.ipynb with {len(cells)} cells")
