#!/usr/bin/env python3
"""Generate index.ipynb for page A2.4 (tanks in series). Run from the page directory.

Quoting convention, copied from A2.5/A2.6/A2.8: markdown cells are raw
triple-DOUBLE-quoted strings and code cells are raw triple-SINGLE-quoted strings,
so a code cell may contain an ordinary Python docstring. Every one is RAW, so a
single backslash here is a single backslash in the notebook.

House rule this page follows strictly: no number that a cell computes is ever
retyped into a markdown cell. Anything with a computed number in it is emitted by
`display(Markdown(f"..."))` from the cell that computed it.
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# ---------------------------------------------------------------- title -----
cells.append(md(r"""---
title: "Tanks in series: what two moments of an RTD do and do not buy you"
description: "Levenspiel's Ch. 14 fitted to its own worked examples, then pushed past them. A five-stage chain with exactly the same mean and variance as E14.4's four tanks leaves 23.2 % of the outlet concentration undetermined at the far end of the sweep, five printed defects are proved from the book's own arithmetic - including the inverted N-to-dispersion rule on Fig. 14.7 - and the model turns out to be the first-order upwind scheme with Pe = 2N."
categories: [sec:A, struct:S1, struct:S3, struct:S5, tier:T0, data:tier6, phase:liquid]
date: 2026-08-07
---

# Tanks in series: what two moments of an RTD do and do not buy you

**Catalog ID:** `A2.4` · **Structures:** `S1` (pointwise / staged algebra),
`S3` (1-D steady BVP), `S5` (convection-dominated) · **Tier:** T0

The tanks-in-series model has one adjustable number. You measure a
residence-time distribution, you take its mean and its variance, you divide, and
you get $N$. Levenspiel's Chapter 14 does exactly that three times, in worked
examples E14.1, E14.2 and E14.4, and it is the reason the model is in every
course: it is the cheapest possible description of a vessel that is neither plug
flow nor a stirred tank.

Fitting $N$ to an RTD and then reporting that the RTD is reproduced would be a
goodness of fit and nothing more. This page does not do that. It fits $N$ where
the book fits it — from **two moments** — and then asks what the fit is worth on
quantities that were held out of it: the *shape* of the exit curve, the
*conversion*, and the *additivity* property the model must satisfy exactly.

The answer is sharper than expected, and it is a limit rather than an agreement.
Levenspiel's E14.4 reads $\bar{t} = 60$ s and $\sigma^2 = 900$ s² off a pair of
tracer curves and concludes "4 tanks". A chain of **five unequal** tanks of
9, 9, 9, 9 and 24 s has exactly the same mean and exactly the same variance — the
two numbers the method uses — and leaves the **outlet concentration** open by
23.2 % at the far end of the range swept. That chain is only an example, not the
extreme: over *every* chain matching the same two moments the outlet
concentration is open by a factor that is computed in closed form below. Two
moments do not determine a conversion, and this page measures by how much.""" ))

# ----------------------------------------------------------- background -----
cells.append(md(r"""## Background

Four gallery pages already end at a residence-time distribution.
[`A2.1`](../A2.1-danckwerts-boundary-conditions/) settles what boundary
conditions the axial-dispersion model needs,
[`A2.3`](../A2.3-taylor-aris-dispersion/) derives a dispersion coefficient,
[`A2.5`](../A2.5-edwards-richardson-dispersion/) measures one and
[`A2.6`](../A2.6-gunn-dispersion-correlations/) correlates it.
[`A2.8`](../A2.8-zwietering-segregation/) asks what an RTD still leaves
undetermined once you have it, and answers with Zwietering's micromixing bounds.

This page is the *other* one-parameter flow model — the one Levenspiel opens
Chapter 14 by putting side by side with the dispersion model:

> "This model can be used whenever the dispersion model is used; and for not too
> large a deviation from plug flow both models give identical results, for all
> practical purposes. Which model you use depends on your mood and taste."

The tanks-in-series model replaces a partial differential equation by a chain of
$N$ equal, perfectly mixed tanks. It has three properties that make it worth a
page of its own:

1. **It is closed-form in $N$** — RTD, $F$ curve, first-order conversion,
   second-order conversion and the macrofluid conversion are all printed in
   Chapter 14 as explicit formulae.
2. **It is additive.** $M$ tanks followed by $N$ more tanks *is* $M+N$ tanks —
   Levenspiel's Eq. 4 — which is what lets E14.4 subtract a sloppy input signal
   from a sloppy output signal instead of deconvolving them.
3. **It is a discretisation.** $N$ equal stirred tanks in series are algebraically
   identical to the first-order upwind finite-volume scheme for plug flow on $N$
   cells. That identity is not in the chapter, it is derived below, and it is what
   makes the chapter's opening claim quantitative: the model's numerical diffusion
   is $u L/(2N)$, so the equivalent Péclet number is $\mathrm{Pe} = 2N$. The
   chapter *does* print a conversion rule of its own — inside the axes box of
   Fig. 14.7, book p. 327 — and it is inverted. The derivation below de-inverts
   it, and the book's own Chapter 13 settles which way round it goes.

**What was read, and from where.** Everything on this page comes from
Levenspiel, *Chemical Reaction Engineering*, 3rd edition (John Wiley & Sons,
1999), ISBN 0-471-25424-X, **Chapter 14 "The Tanks-In-Series Model", sections
14.1 (Pulse Response Experiments and the RTD) and 14.2 (Chemical Conversion)**,
book pp. 321–338, plus **book p. 303 in Chapter 13**, which Fig. 14.7's own
annotation sends the reader to by name and which supplies the two printed values
that settle the inversion. The book **attributes** the RTD, its means and variances, and
the $F$ curve of Fig. 14.7 to *MacMullin, R. B. and Weber, M., Jr., Trans.
AIChE **31**, 409 (1935)*. That paper is **not on disk and was not consulted**;
it is recorded as the origin and this monograph as the text actually read. No
equation on this page was written from memory.""" ))

# ------------------------------------------------- the published model ------
cells.append(md(r"""## The published model

Levenspiel defines two dimensionless times, on p. 321:

$$\theta_i = \frac{t}{\bar{t}_i}\ \text{(per tank)},\qquad
  \theta = \frac{t}{\bar{t}}\ \text{(all $N$ tanks)},\qquad
  \theta_i = N\theta .$$

### The RTD — his Eq. 3

Four equivalent forms, all boxed together as **Eq. 3** on book p. 323, with the
means and variances alongside:

$$\bar{t}\,\mathbf{E} = \left(\frac{t}{\bar t}\right)^{N-1}\frac{N^N}{(N-1)!}\,e^{-tN/\bar t}
  \quad\cdots\quad \bar t = N \bar t_i \quad\cdots\quad \sigma^2 = \frac{\bar t^{\,2}}{N}$$

$$\bar{t}_i\,\mathbf{E} = \left(\frac{t}{\bar t_i}\right)^{N-1}\frac{1}{(N-1)!}\,e^{-t/\bar t_i}
  \quad\cdots\quad \bar t_i = \frac{\bar t}{N} \quad\cdots\quad \sigma^2 = N \bar t_i^{\,2}$$

$$\mathbf{E}_{\theta_i} = \bar t_i \mathbf{E} = \frac{\theta_i^{\,N-1}}{(N-1)!}e^{-\theta_i}
  \quad\cdots\quad \sigma^2_{\theta_i} = N$$

$$\mathbf{E}_{\theta} = (N\bar t_i)\mathbf{E} = N\frac{(N\theta)^{N-1}}{(N-1)!}e^{-N\theta}
  \quad\cdots\quad \sigma^2_{\theta} = \frac{1}{N}$$

The chapter's own text says these "were first derived by MacMullin and Weber
(1935)".

### The properties of the curve — his Fig. 14.3

Figure 14.3 is a sketch, but it is annotated with five *algebraic* statements and
two numbers, and all seven are testable:

$$\mathbf{E}_{\theta,\max} = \frac{N(N-1)^{N-1}}{(N-1)!}e^{-(N-1)}
  \;\cong\; \frac{N}{\sqrt{2\pi(N-1)}},\ \text{"error} < 2\% \text{ for } N > 5\text{"}$$

$$\theta_{\max} = \frac{N-1}{N},\qquad
  \frac{\Delta\theta}{\theta_{\max}} = \frac{2}{\sqrt{N-1}},\qquad
  \text{Area} = 1,\qquad \sigma^2_\theta = \tfrac1N,$$

$$\mathbf{E}_{\theta,\mathrm{inf}} \cong 0.55\,\mathbf{E}_{\theta,\max}\ (N=4),
  \qquad \cong 0.61\,\mathbf{E}_{\theta,\max}\ (N\ge 10).$$

$\Delta\theta$ is the separation of the two inflection points.

### Independence — his Eq. 4

> "If $M$ tanks are connected to $N$ more tanks (all of the same size) then the
> individual means and variances (in ordinary time units) are additive"

$$\bar t_{M+N} = \bar t_M + \bar t_N \quad\text{and}\quad
  \sigma^2_{M+N} = \sigma^2_M + \sigma^2_N \tag{4}$$

with a footnote defining independence as "the fluid loses its memory as it passes
from vessel to vessel", and warning that laminar flow often does not satisfy it.

### The one-shot tracer relation — his Eq. 5

Under the heading *One-shot Tracer Input*, on book p. 324, Eq. 4 is turned into
the relation E14.4 actually uses:

$$\Delta\sigma^2 = \sigma^2_{\text{out}} - \sigma^2_{\text{in}}
  = \frac{(\Delta\bar t)^{\,2}}{N} \tag{5}$$

This is the most heavily exercised relation on the page: running it backwards is
the whole of E14.4's fit, and the headline in section 8 is a statement about it.
Eqs. 6a–c and 7a–c on p. 325, which superpose passes through a **closed
recirculation** loop, are printed but are not exercised anywhere here; nor is
Fig. 14.8's graphical construction for arbitrary kinetics (p. 328).

### The $F$ curve — his Eq. 8

$$\mathbf{F} = 1 - e^{-N\theta}\left[1 + N\theta + \frac{(N\theta)^2}{2!} + \cdots
  + \frac{(N\theta)^{N-1}}{(N-1)!}\right] \tag{8}$$

with a brace under the terms reading "For one tank use the first term / For
$N = 2$ ...". Figure 14.7, which plots it, is captioned "from MacMullin and Weber
(1935)" — and inside its axes box it carries the chapter's **only** rule for
turning $N$ into a dispersion number:

> "When $N > 50$ the curve becomes symmetrical in which case use fig. 13-11
> with $N = \tfrac12\left(\dfrac{\mathbf{D}}{uL}\right)$"

The ratio in that parenthesis is upside down. The book's own Chapter 13 proves
it, section 9 shows what the printed form costs, and it is reported as the most
consequential of the printed defects rather than repaired.

### Conversion — section 14.2

First order, one tank and then $N$ tanks (**Eq. 9**):

$$\frac{C_A}{C_{A0}} = \frac{1}{1+k\bar t_i} = \frac{1}{1+k\bar t}
  \qquad\longrightarrow\qquad
  \frac{C_A}{C_{A0}} = \frac{1}{(1+k\bar t_i)^N} = \frac{1}{\left(1+\frac{k\bar t}{N}\right)^N}$$

and, "for small deviations from plug flow (large $N$)",

$$\text{for same } C_{A,\text{final}}:\ \frac{V_{N\,\text{tanks}}}{V_p} = 1 + k\bar t_i
  = 1 + \frac{k\bar t}{2N},
\qquad
\text{for same volume } V:\ \frac{C_{A,N\text{tanks}}}{C_{Ap}} = 1 + \frac{(k\bar t)^2}{2N}.$$

**The first of those two is printed inconsistently** — the middle and right-hand
expressions differ by a factor of two — and the page proves which one is wrong
from the book's own Eq. 9 rather than from anywhere else.

Second order for a microfluid, $A\to R$ or $A+B\to R$ with $C_{A0}=C_{B0}$,
"Eq. 6.8 gives" (**Eq. 10**) a nest of $N$ radicals:

$$C_N = \frac{1}{4k\tau_i}\left(-2 + 2\sqrt{-1\cdots+2\sqrt{-1+2\sqrt{1+4C_0k\tau_i}}}\right)$$

and for a macrofluid, "combine Eq. 11.3 with Eq. 3" (**Eq. 11**):

$$\frac{C_A}{C_{A0}} = \frac{N^N}{(N-1)!\,\bar t^{\,N}}
  \int_0^{\infty}\left(\frac{C_A}{C_{A0}}\right)_{\text{batch}} t^{N-1}e^{-tN/\bar t}\,dt$$

The book sets that denominator as $(N-1)!\,\bar t_N$, with $N$ **lowered** as a
subscript. Only $\bar t^{\,N}$ makes the expression dimensionless and normalised,
and the book prints the raised form itself two pages later in E14.4, so the
subscript is a typesetting slip. It is reported below, not silently repaired.""" ))

# ------------------------------------------- parameters and assumptions -----
cells.append(md(r"""## Parameters and assumptions

There is **one parameter**, $N$, and it is an integer count of equal tanks. Two
groups set everything else: $k\bar t$ for a first-order reaction and
$kC_0\bar t$ for a second-order one. Nothing on this page is calibrated to
anything: every $N$ is either an integer chosen to sweep, or the value the book
itself computes from its own printed moments.

Assumptions the model makes, in the order they bite:

- **Each tank is perfectly mixed**, so its exit stream carries the tank's own
  concentration. This is not a boundary-condition choice; it is the model, and
  it is why the outflow below is *not* the second-order face reconstruction that
  pymrm's zero-gradient boundary condition would give. The difference is
  measured.
- **All tanks are equal in size.** Levenspiel says so explicitly in Eq. 4's
  sentence. The page's headline result is what goes wrong when this is assumed
  rather than checked.
- **Stages are independent** — no memory carried from tank to tank. Levenspiel's
  own footnote flags laminar flow as a case where this fails. Nothing on this
  page tests it; it is an assumption inherited whole.
- **Constant density and constant volumetric flow**, so $\bar t = V/v$ and
  $\bar t_i = \bar t/N$.

**Fit or test.** Every number this page compares against is a **held-out** one.
$N = 4$ in E14.4 is fitted — to two moments, by the book — and then the RTD
shape, the $F$ curve, the first- and second-order conversions and the macrofluid
conversion are all *consequences* of that fit that were not used to make it.
E14.1 and E14.2 fit nothing at all: they are pure arithmetic on printed inputs.
The Fig. 14.3 properties are fitted to nothing anywhere. There is **no experiment
in Chapter 14 at all**, so nothing here is experimental validation, and the page
says so wherever it matters.""" ))

# ----------------------------------------------------------- environment ----
cells.append(code(r'''# Colab: install pymrm if it is not already present.
try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm'''))

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

import textwrap
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from scipy.sparse import eye_array, identity
from scipy.sparse.linalg import splu
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln
from pymrm import (construct_convflux_upwind, construct_div, construct_grad,
                   construct_coefficient_matrix, compute_boundary_values,
                   interp_cntr_to_stagg_tvd, vanleer, NumJac, newton)
from gallery_utils import load_data, load_meta, cite_data, report_agreement

warnings.filterwarnings("ignore", category=UserWarning)
PAGE = "A2.4-tanks-in-series"
rng = np.random.default_rng(20260807)   # nothing here is stochastic; seeded anyway
plt.rcParams.update({"figure.dpi": 110, "font.size": 9.5})
np.set_printoptions(legacy="1.25")'''))

# -------------------------------------------------------------- the data ----
cells.append(md(r"""## The data

Chapter 14 contains **no measurement**. Its numbers are (a) the inputs and
intermediates of three worked examples, (b) the algebraic annotations of
Fig. 14.3, (c) two approximate ratios on that same figure, and (d) the threshold
and the coefficient of Fig. 14.7's dispersion-conversion annotation. Two further
rows come from **Chapter 13, book p. 303** — Fig. 13.12's printed
$\sigma^2_\theta = 2(\mathbf{D}/uL)$ and the printed validity limit
$\mathbf{D}/uL < 0.01$ of Fig. 13.11 — because Fig. 14.7's annotation names that
figure and those are the two numbers that settle whether its ratio is the right
way up. All of them are in the CSV below, transcribed from 600 dpi renders — the
scan's *native* resolution, checked with `pdfimages -list`, so anything higher
would be interpolation.

The text layer of this book is a specific kind of trap. The **prose** is
excellent and fully searchable; it is what located the chapter and its two
sections in the first place. The **equations and subscripts are destroyed**:
`pdftotext` returns Eq. 3 as `7. = (-&) NN e-tNli` and Eq. 4 as `%+, = tM+ t,`.
Not one numeric cell below came off the text layer.

Five of the printed items are **defective as printed**. They are carried in the
CSV exactly as they appear, with the defect named in the `printed_as` and `note`
columns, and each is proved from the book's own arithmetic in the *Printed
defects* section rather than quietly corrected.

**No other gallery page's dataset is loaded here**, so no cross-page
reconciliation applies. The one sibling page whose findings bear on this one is
[`A2.8`](../A2.8-zwietering-segregation/), which uses Erlang RTDs of $N$ equal
tanks and states in its own caveats that "`A2.4` is not built, so the RTDs used
here have no page of their own to cross-check against". It publishes no RTD
dataset, so there is nothing to borrow; what it establishes about these RTDs —
that each integrates to one and has mean one, and that for a reaction order other
than one the RTD alone does not fix the conversion — is re-established
independently below and is consistent.""" ))

cells.append(code(r'''PRINTED = load_data("levenspiel-1999-ch14-printed.csv", page=PAGE)
META = load_meta("levenspiel-1999-ch14-printed.csv", page=PAGE)


def P(item, quantity):
    """One printed value, by (worked example, quantity). Never retype a number."""
    row = PRINTED[(PRINTED["item"] == item) & (PRINTED["quantity"] == quantity)]
    if len(row) != 1:
        raise KeyError(f"{item}/{quantity}: {len(row)} rows")
    return float(row["value"].iloc[0])


display(Markdown(f"**Source:** {cite_data(META)}"))
display(PRINTED[["item", "quantity", "value", "units", "book_page", "printed_as"]])
display(Markdown(
    f"{len(PRINTED)} printed values, from book pp. "
    f"{PRINTED['book_page'].min()}-{PRINTED['book_page'].max()} "
    f"(PDF pp. {PRINTED['book_page'].min() + 16}-{PRINTED['book_page'].max() + 16}). "
    "Nothing is digitised: the page traces no curve anywhere."))'''))

# ------------------------------------------------ pymrm implementation ------
cells.append(md(r"""## PyMRM implementation

### The cascade is the first-order upwind scheme

A tracer balance on tank $i$ of volume $V_i = V/N$ with volumetric flow $v$ is

$$V_i\frac{dC_i}{dt} = v\,C_{i-1} - v\,C_i
  \qquad\Longleftrightarrow\qquad
  \frac{dC_i}{dt} = \frac{C_{i-1}-C_i}{\bar t_i},\quad \bar t_i = \frac{V_i}{v}.$$

Put the chain on a line $0\le z\le 1$ with $N$ uniform cells, superficial
velocity $u = 1/\bar t$ (so the transit time of the whole line is $\bar t$) and
$\Delta z = 1/N$. The first-order **upwind** convective flux on an interior face
is $u\,C_{i-1}$, and the divergence operator divides by the cell size, so

$$\left[\nabla\!\cdot(uC)\right]_i = \frac{u\,C_i - u\,C_{i-1}}{\Delta z}
 = \frac{C_i - C_{i-1}}{\bar t_i}.$$

The two are the same equation. `construct_convflux_upwind` plus `construct_div`
with `nu=0` **is** the tanks-in-series model — not an approximation to it — and
$N$ is a *physical parameter*, not a mesh resolution. That has a consequence the
rest of this page keeps running into: **there is no grid-refinement axis here.**
Refining the grid changes the model. The axes that do carry numerical error are
the time step, the quadrature panel count and (in the dispersion comparison at
the end) a genuine grid.

### The one place the standard boundary condition is wrong

A perfectly mixed tank discharges its own contents, so the flux leaving tank $N$
is exactly $v\,C_N$. pymrm's boundary machinery does not offer a "pure outflow"
condition: the nearest choice, a zero-gradient $\{a{=}1, b{=}0, d{=}0\}$, makes
the operator reconstruct the exit-face value to second order as
$\tfrac98 C_N - \tfrac18 C_{N-1}$. For a *discretised PDE* that is the right
thing to do and buys an order of accuracy. Here it is a **modelling error**,
because $N$ is not a mesh: it changes the last tank's balance and breaks the
identity with Levenspiel's Eq. 9. It is measured against $N$ in *Validation*, and
it decays like an ordinary first-order discretisation error — which is exactly
why it is easy to dismiss, and why it is in the break table.

So the exit face carries no operator flux (`b=1, d=0` suppresses it) and the
outflow enters as a cell-wise sink $v C_N / V_N$ on the last tank, built with
`construct_coefficient_matrix`. The resulting matrix, scaled by $\bar t_i$, is
printed below and is exactly the lower-bidiagonal cascade.""" ))

cells.append(code(r'''class Cascade:
    """N stirred tanks in series, assembled from the pymrm operators.

    The cell count IS Levenspiel's N, so nothing here is a discretisation of
    anything; `faces` exists only so that the unequal-tank chain in the Results
    section, and the break table, can use the same class.
    """

    def __init__(self, N, tbar=1.0, c_in=0.0, nu=0, outlet="stirred",
                 faces=None, v_sign=1.0):
        self.N, self.tbar, self.c_in = N, tbar, c_in
        self.shape = (N, 1)                      # (cells, fields) - never (N,)
        self.z_f = np.linspace(0.0, 1.0, N + 1) if faces is None else np.asarray(faces, float)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.dz = np.diff(self.z_f)
        self.u = v_sign / tbar                   # transit time of the whole line = tbar
        self.t_i = self.dz / abs(self.u)         # residence time of each tank

        # inlet: the feed concentration is imposed at the face -> a=0, b=1, d=c_in
        bc_in = {"a": 0.0, "b": 1.0, "d": c_in}
        if outlet == "stirred":
            # The exit stream of a stirred tank carries the tank's own value, so the
            # operator must put NO flux on the exit face: a=0,b=1,d=0 sets that face
            # value to zero and therefore its flux to zero.  The physical outflow
            # v*C_N is added below as a sink on the last tank alone.
            bc = (bc_in, {"a": 0.0, "b": 1.0, "d": 0.0})
            conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                      bc, v=self.u)
            div = construct_div(self.shape, self.z_f, nu=nu)   # nu=0: Cartesian
            coef = np.zeros((N, 1))
            coef[-1, 0] = abs(self.u) / self.dz[-1]
            self.A = (div @ conv
                      + construct_coefficient_matrix(coef, shape=self.shape)).tocsc()
        elif outlet == "zerograd":
            # The usual PDE outflow condition: dc/dn = 0 at z=1.  n = +z, so a=1,
            # b=0, d=0.  Correct for a discretised PDE, WRONG for this model.
            bc = (bc_in, {"a": 1.0, "b": 0.0, "d": 0.0})
            conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                      bc, v=self.u)
            div = construct_div(self.shape, self.z_f, nu=nu)
            self.A = (div @ conv).tocsc()
        else:
            raise ValueError(outlet)
        self.bc = bc
        self.b = np.asarray((div @ conv_bc).todense()).ravel()

    # ---- steady state ---------------------------------------------------
    def first_order(self, k):
        """Exit concentration for r = k C.  Linear: one factorisation, no Newton."""
        c = splu((self.A + k * eye_array(self.N, format="csc")).tocsc()).solve(-self.b)
        return c[-1], c

    def order_n(self, k, order=2.0, c_guess=None):
        """Exit concentration for r = k C^order, by Newton with a pointwise NumJac."""
        nj = NumJac(self.shape)          # (N,1): last axis is the field -> block-diagonal
        A, b = self.A, self.b

        def res_jac(c):
            c = np.asarray(c).reshape(self.shape)
            g, Jr = nj(lambda x: k * np.maximum(x, 0.0) ** order, c)
            return A @ c.ravel() + b + np.asarray(g).ravel(), (A + Jr).tocsc()

        guess = np.full(self.shape, self.c_in if c_guess is None else c_guess)
        sol = newton(res_jac, guess, tol=1e-14)
        c = np.asarray(sol.x).ravel()
        return c[-1], c, sol

    # ---- transient ------------------------------------------------------
    def pulse_response(self, tmax, nt, theta=0.5):
        """E(t) from a pulse into tank 1, by a fixed-step theta scheme.

        theta = 0.5 is Crank-Nicolson (2nd order), theta = 1 implicit Euler (1st).
        The injected amount is 1 mol into V_1, so E = u * C_N with no numerical
        normalisation anywhere - which is what lets 'Area = 1' be a real check.
        """
        c = np.zeros(self.N)
        c[0] = 1.0 / self.dz[0]
        dt = tmax / nt
        lu = splu((identity(self.N, format="csc") + theta * dt * self.A).tocsc())
        ts = np.empty(nt + 1)
        es = np.empty(nt + 1)
        ts[0], es[0] = 0.0, abs(self.u) * c[-1]
        for j in range(nt):
            c = lu.solve(c - (1 - theta) * dt * (self.A @ c) - dt * self.b)
            ts[j + 1], es[j + 1] = (j + 1) * dt, abs(self.u) * c[-1]
        return ts, es


_demo = Cascade(4, tbar=60.0)
print("A * t_i  for N = 4  (should be the lower-bidiagonal cascade):")
print(_demo.A.toarray() * _demo.t_i[0])
print(f"\ntank residence times t_i = {_demo.t_i} s, sum = {_demo.t_i.sum():g} s")'''))

cells.append(md(r"""### Levenspiel's own formulae, written once

Every closed form below is transcribed from a 600 dpi crop of the boxed equation
it is named after. They are kept as free functions so that the pymrm solve and
the algebra never share a line of code — that separation is what makes the
comparisons in *Validation* mean anything.

$\mathbf{E}$ is evaluated through `gammaln`, not by forming $N^N$ and $(N-1)!$,
because both overflow before $N = 200$ and the plug-flow sweep below goes to
$N = 400$. Eq. 8's $F$ is summed term by term exactly as the book writes it, so
that the truncation the book's brace prescribes is the one that is tested.""" ))

cells.append(code(r'''def E_time(t, N, tbar):
    """Levenspiel Eq. 3, second form: E(t) for N equal tanks of total mean tbar."""
    ti = tbar / N
    t = np.maximum(np.asarray(t, float), 1e-300)
    return np.exp((N - 1) * np.log(t / ti) - t / ti - gammaln(N)) / ti


def E_theta(th, N):
    """Levenspiel Eq. 3, fourth form: E_theta = N (N theta)^(N-1) e^(-N theta)/(N-1)!."""
    th = np.maximum(np.asarray(th, float), 1e-300)
    return np.exp(np.log(N) + (N - 1) * np.log(N * th) - N * th - gammaln(N))


def F_eq8(th, N):
    """Levenspiel Eq. 8, summed exactly as printed: 1 + Ntheta + ... + (Ntheta)^(N-1)/(N-1)!."""
    th = np.asarray(th, float)
    x = N * th
    term = np.ones_like(x)
    total = term.copy()
    for j in range(1, N):
        term = term * x / j
        total = total + term
    return 1.0 - np.exp(-x) * total


def eq9(N, ktbar):
    """Levenspiel Eq. 9: first order, N tanks."""
    return 1.0 / (1.0 + ktbar / N) ** N


def eq10(N, k, tau_i, C0):
    """Levenspiel Eq. 10: second-order microfluid, as the printed nest of N radicals."""
    x = np.sqrt(1.0 + 4.0 * C0 * k * tau_i)
    for _ in range(N - 1):
        x = np.sqrt(-1.0 + 2.0 * x)
    return (-2.0 + 2.0 * x) / (4.0 * k * tau_i)


def eq11(N, tbar, batch, npanel=40000, span=80.0):
    """Levenspiel Eq. 11: macrofluid conversion, by composite Simpson on [0, span*tbar].

    The printed denominator is (N-1)! * tbar_N with a LOWERED N.  Only tbar**N
    normalises, and E14.4 prints the raised form, so tbar**N is used and the
    subscript reading is reported as a defect, and exercised in the break table.
    """
    t = np.linspace(1e-12, span * tbar, npanel + 1)
    y = np.asarray(batch(t), float) * E_time(t, N, tbar)
    h = t[1] - t[0]
    return h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())


def max_brent(f, lo, hi, nscan=400):
    """max of a smooth scalar f on [lo, hi], BRACKETED by a coarse scan and then
    refined with Brent's method.  Every curve maximum on this page goes through
    here rather than through `.max()` of a sampled array: section 2's own
    sweep-versus-root-find margin is the reason."""
    xs = np.linspace(lo, hi, nscan)
    ys = np.array([f(x) for x in xs])
    j = int(np.argmax(ys))
    a, b = xs[max(j - 1, 0)], xs[min(j + 1, nscan - 1)]
    r = minimize_scalar(lambda x: -f(x), bounds=(a, b), method="bounded",
                        options={"xatol": 1e-10})
    return float(max(-r.fun, ys[j]))


def Emax_printed(N):
    """The closed form printed in Fig. 14.3 for E_theta,max."""
    return np.exp(np.log(N) + (N - 1) * np.log(N - 1) - (N - 1) - gammaln(N))


def Emax_stirling(N):
    """Fig. 14.3's approximation: N / sqrt(2 pi (N-1))."""
    return N / np.sqrt(2.0 * np.pi * (N - 1))


print("sanity: eq9(4, 2) =", eq9(4, 2.0),
      " eq10(1, 2, 1, 1) =", eq10(1, 2.0, 1.0, 1.0),
      " F_eq8(1e9, 3) =", F_eq8(1e9, 3))'''))

# ---------------------------------------------------------------- results ---
cells.append(md(r"""## Results

### 1. The three worked examples, reproduced from their printed inputs

E14.1 (*Modifications to a winery*) and E14.2 (*A fable on river pollution*) fit
nothing: they apply $\sigma^2 \propto N \propto L$ — Levenspiel's own sentence,
"for small deviations from plug flow, from Eq. 3 $\sigma^2 \propto N$ or
$\sigma^2 \propto L$" — to printed inputs and produce printed answers. They are
reproduced here because they are the only place in the chapter where the
model's *scaling with length* is exercised, and because getting the square wrong
is the single most likely way to misuse it.

E14.2's answer is obtained by **root-finding** the printed ratio
$14/10.5 = \sqrt{L/(L-119)}$ rather than by evaluating a rearrangement, and
checked against the closed form $L = 119/(1-(10.5/14)^2)$ — two routes that share
no algebra.""" ))

cells.append(code(r'''# ---- E14.1: the winery -------------------------------------------------
L1, S1, S1SQ, L2 = P("E14.1", "L_1"), P("E14.1", "sigma_1"), P("E14.1", "sigma_1_squared"), P("E14.1", "L_2")
E141_S2SQ_P, E141_S2_P = P("E14.1", "sigma_2_squared"), P("E14.1", "sigma_2")

assert S1SQ == S1 ** 2, "the book's own sigma_1^2 must be sigma_1 squared"
E141_S2SQ = L2 / L1 * S1SQ            # sigma^2 proportional to L
E141_S2 = np.sqrt(E141_S2SQ)
E141_S2SQ_DEV = abs(E141_S2SQ - E141_S2SQ_P)
E141_S2_DEV = abs(E141_S2 - E141_S2_P)
E141_WRONG_POWER = S1 * L2 / L1       # what scaling sigma (not sigma^2) would give

# ---- E14.2: the Ohio river --------------------------------------------
X, SP, SC = P("E14.2", "x_Portsmouth_to_Cincinnati"), P("E14.2", "spread_Portsmouth"), P("E14.2", "spread_Cincinnati")
E142_L_P = P("E14.2", "L")
E142_L = brentq(lambda L: SC / SP - np.sqrt(L / (L - X)), X * (1 + 1e-12), 1e6,
                xtol=1e-12, rtol=8.9e-16)
E142_L_CLOSED = X / (1.0 - (SP / SC) ** 2)          # second, independent route
E142_L_DEV = abs(E142_L - E142_L_P)
E142_TWO_ROUTES = abs(E142_L - E142_L_CLOSED)
E142_WRONG_POWER = X / (1.0 - SP / SC)              # spread ~ distance, not spread^2

display(Markdown(f"""
**E14.1** — pipeline {L1:g} m &rarr; {L2:g} m, spread {S1:g} bottles.
$\\sigma_2^2 = ({L2:g}/{L1:g})({S1SQ:g}) = {E141_S2SQ:.6g}$ against the printed
{E141_S2SQ_P:g} (dev {E141_S2SQ_DEV:.1e}), so
$\\sigma_2 = {E141_S2:.6g}$ against the printed {E141_S2_P:g}
(dev {E141_S2_DEV:.1e}) &mdash; **{E141_S2:.0f} bottles of vin rose**.
Scaling $\\sigma$ instead of $\\sigma^2$ would give {E141_WRONG_POWER:.4g} bottles,
{abs(E141_WRONG_POWER - E141_S2) / E141_S2 * 100:.0f} % high.

**E14.2** — root-finding the printed ratio {SC:g}/{SP:g} = $\\sqrt{{L/(L-{X:g})}}$
gives $L = {E142_L:.6f}$ miles against the printed {E142_L_P:g}
(dev {E142_L_DEV:.1e}); the closed form $L = {X:g}/(1-({SP:g}/{SC:g})^2)$ gives
{E142_L_CLOSED:.6f}, the two routes agreeing to {E142_TWO_ROUTES:.1e}. Using the
spread rather than its square would put the source at {E142_WRONG_POWER:.0f} miles
&mdash; past Marietta and Parkersburg, and a different answer to the question asked.
"""))'''))

cells.append(md(r"""### 2. E14.4: the fit, and everything it was not fitted to

E14.4 is the one place in the chapter where $N$ is obtained from data rather than
assumed, and it is the reason Eq. 4 exists. A sloppy tracer input with
$\bar t = 220$ s, $\sigma^2 = 100$ s² produces an output with $\bar t = 280$ s,
$\sigma^2 = 1000$ s². Rather than deconvolve, Levenspiel *subtracts moments*:

$$\Delta \bar t = 280-220 = 60\ \text{s},\qquad
  \Delta(\sigma^2) = 1000-100 = 900\ \text{s}^2,\qquad
  N = \frac{(\Delta\bar t)^2}{\Delta(\sigma^2)} = 4\ \text{tanks}.$$

That is **two numbers in, one number out.** Everything else in this section is a
consequence of those two numbers that was not used to produce them.

The prefactor of the resulting $E$ curve is computed two independent ways: from
$N^N/((N-1)!\,\bar t^{\,N})$, and by normalising $t^{3}e^{-4t/60}$ to unit area
with an independent quadrature. Neither route uses the other's algebra.""" ))

cells.append(code(r'''TB_IN, V_IN = P("E14.4", "tbar_in"), P("E14.4", "var_in")
TB_OUT, V_OUT = P("E14.4", "tbar_out"), P("E14.4", "var_out")
DTB_P, DV_P, N4_P = P("E14.4", "delta_tbar"), P("E14.4", "delta_var"), P("E14.4", "N")
PRE_P, EXP_P = P("E14.4", "E_prefactor"), P("E14.4", "E_exponent")
TICK_P = P("E14.4", "E_axis_top_tick")

# --- the fit: Eq. 4 run backwards -----------------------------------------
DTB, DV = TB_OUT - TB_IN, V_OUT - V_IN
E144_DTB_DEV, E144_DV_DEV = abs(DTB - DTB_P), abs(DV - DV_P)
E144_N = DTB ** 2 / DV
E144_N_DEV = abs(E144_N - N4_P)
N4, TB4 = int(round(E144_N)), DTB

# --- the E curve: two independent routes to the prefactor -----------------
E144_PRE = N4 ** N4 / (np.exp(gammaln(N4)) * TB4 ** N4)
_span, _npan = 200.0 * TB4, 400000
_t = np.linspace(1e-12, _span, _npan + 1)
_y = _t ** (N4 - 1) * np.exp(-N4 * _t / TB4)
_h = _t[1] - _t[0]
_area = _h / 3.0 * (_y[0] + _y[-1] + 4.0 * _y[1:-1:2].sum() + 2.0 * _y[2:-1:2].sum())
E144_PRE_QUAD = 1.0 / _area                     # second route: normalise numerically
E144_PRE_DEV = abs(E144_PRE - PRE_P) / PRE_P
E144_PRE_TWO_ROUTES = abs(E144_PRE - E144_PRE_QUAD) / E144_PRE
E144_EXP = N4 / TB4
E144_EXP_DEV = abs(E144_EXP - EXP_P)

# --- the peak: root-found, never sampled ----------------------------------
E144_TMAX = brentq(lambda t: (N4 - 1) / t - N4 / TB4, 1e-9, 10 * TB4, xtol=1e-13)
E144_EMAX = E144_PRE * E144_TMAX ** (N4 - 1) * np.exp(-N4 * E144_TMAX / TB4)
_grid = np.linspace(0, 4 * TB4, 1001)           # what a sampled max would have given
E144_EMAX_SAMPLED = float(E_time(_grid, N4, TB4).max())
# The SAME peak computed from the coefficients as the book ROUNDS AND PRINTS them
# (3.2922e-6 and 0.0667, both from the CSV), which is the curve the axis tick of
# Fig. E14.4b belongs to.  It is not the same number as E144_EMAX above, and the
# defect section must use this one.
E144_TMAX_PRINTED = brentq(lambda t: (N4 - 1) / t - EXP_P, 1e-9, 10 * TB4, xtol=1e-13)
E144_EMAX_PRINTED = PRE_P * E144_TMAX_PRINTED ** (N4 - 1) * np.exp(-EXP_P * E144_TMAX_PRINTED)
E144_TICK_ERR = abs(TICK_P - E144_EMAX_PRINTED) / E144_EMAX_PRINTED

display(Markdown(f"""
| step | this page | Levenspiel |
|---|---|---|
| $\\Delta\\bar t$ | {DTB:g} s | {DTB_P:g} s (dev {E144_DTB_DEV:.0e}) |
| $\\Delta(\\sigma^2)$ | {DV:g} s² | {DV_P:g} s² (dev {E144_DV_DEV:.0e}) |
| $N=(\\Delta\\bar t)^2/\\Delta(\\sigma^2)$ | {E144_N:.10g} | {N4_P:g} tanks (dev {E144_N_DEV:.0e}) |
| prefactor $N^N/((N-1)!\\,\\bar t^{{\\,N}})$ | {E144_PRE:.6e} s⁻⁴ | {PRE_P:.4e} (rel dev {E144_PRE_DEV:.1e}) |
| the same by normalising $t^3e^{{-4t/60}}$ | {E144_PRE_QUAD:.6e} s⁻⁴ | two routes agree to {E144_PRE_TWO_ROUTES:.1e} |
| exponent $N/\\bar t$ | {E144_EXP:.6g} s⁻¹ | {EXP_P:g} (dev {E144_EXP_DEV:.1e}, the book's rounding) |

The fitted $N$ is **exactly** {E144_N:.0f}: not rounded to it, equal to it, because
$60^2 = 3600$ and $4\\times 900 = 3600$.

The curve peaks at $t = {E144_TMAX:.6g}$ s with
$E = {E144_EMAX:.6g}$ s⁻¹, root-found from $dE/dt = 0$. Sampling the same curve on
a 1001-point grid over $[0, 4\\bar t]$ returns {E144_EMAX_SAMPLED:.6g} s⁻¹ instead
&mdash; close, and still a grid maximum. Evaluated instead with the coefficients
as the book *rounds and prints* them ({PRE_P:.4e} and {EXP_P:g}), the same curve
peaks at $t = {E144_TMAX_PRINTED:.6g}$ s with
$E = {E144_EMAX_PRINTED:.6g}$ s⁻¹ &mdash; a different number in the fifth digit, and
the one Fig. E14.4b's axis belongs to. Its top tick is printed
**{TICK_P:g}** s⁻¹, which is {E144_TICK_ERR * 100:.0f} % away from that maximum;
see *Printed defects*.
"""))'''))

cells.append(md(r"""### 3. Fig. 14.3's seven annotations

Figure 14.3 is drawn, not computed, but everything written on it is exact
algebra. Each statement below is obtained by **root-finding on Eq. 3** —
$dE_\theta/d\theta = 0$ for the maximum, $d^2E_\theta/d\theta^2 = 0$ for the two
inflection points — never by scanning a grid, and then compared with the printed
expression.

The Stirling annotation is the sharpest of them. "error $< 2\%$ for $N > 5$" is a
*threshold* claim, so it is tested by root-finding the $N$ at which the relative
error of $N/\sqrt{2\pi(N-1)}$ crosses 2 %, not by checking a few integers and
declaring victory.""" ))

cells.append(code(r'''NS = np.array([2, 3, 4, 5, 6, 8, 10, 16, 25, 50, 100, 200, 400])


def theta_max_root(N):
    return brentq(lambda th: (N - 1) / th - N, 1e-12, 50.0, xtol=1e-14, rtol=8.9e-16)


def inflections(N):
    """The two inflection points of Eq. 3, found by BRACKETING A SIGN CHANGE.

    d2E/dtheta2 = 0 is (N-1)(N-2)/u^2 - 2(N-1)/u + 1 = 0 with u = N theta.  The
    brackets come from a coarse scan, never from the analytic roots, so this is a
    root-find and not an evaluation.  At N = 2 the lower inflection sits exactly
    at theta = 0 (the density leaves the origin with zero curvature there), the
    scan finds only one interior sign change, and 0 is returned for the lower one.
    """
    g = lambda th: (N - 1) * (N - 2) / (N * th) ** 2 - 2 * (N - 1) / (N * th) + 1.0
    scan = np.geomspace(1e-9, 60.0, 20001)
    vals = np.array([g(s) for s in scan])
    idx = np.nonzero(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
    roots = [brentq(g, scan[i], scan[i + 1], xtol=1e-15, rtol=8.9e-16) for i in idx]
    if len(roots) == 1:
        return 0.0, roots[0]
    if len(roots) != 2:
        raise RuntimeError(f"N={N}: found {len(roots)} inflection points, expected 2")
    return roots[0], roots[1]


rows = []
for N in NS:
    tm = theta_max_root(N)
    lo, hi = inflections(N)
    em_root = E_theta(tm, N)
    rows.append(dict(N=N,
                     theta_max=tm, theta_max_printed=(N - 1) / N,
                     width_ratio=(hi - lo) / tm, width_printed=2 / np.sqrt(N - 1),
                     Emax_root=em_root, Emax_printed=Emax_printed(N),
                     stirling_rel_err=abs(Emax_stirling(N) - Emax_printed(N)) / Emax_printed(N),
                     Einf_lo=E_theta(lo, N) / em_root, Einf_hi=E_theta(hi, N) / em_root))
FIG143 = pd.DataFrame(rows)
FIG143["theta_max_dev"] = (FIG143.theta_max - FIG143.theta_max_printed).abs()
FIG143["width_dev"] = (FIG143.width_ratio - FIG143.width_printed).abs()
FIG143["Emax_dev"] = ((FIG143.Emax_root - FIG143.Emax_printed) / FIG143.Emax_printed).abs()

THETA_MAX_DEV = float(FIG143.theta_max_dev.max())
WIDTH_DEV = float(FIG143.width_dev.max())
EMAX_DEV = float(FIG143.Emax_dev.max())

# Area = 1 and sigma^2_theta = 1/N, by quadrature on Eq. 3 (independent of the algebra)
_th = np.linspace(1e-12, 60.0, 600001)
_hh = _th[1] - _th[0]
def _simp(y):
    return _hh / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())
AREA_DEV, MEAN_DEV, VAR_DEV = 0.0, 0.0, 0.0
for N in [2, 4, 10, 50]:
    e = E_theta(_th, N)
    a, m = _simp(e), _simp(_th * e)
    v = _simp((_th - 1.0) ** 2 * e)
    AREA_DEV = max(AREA_DEV, abs(a - 1.0))
    MEAN_DEV = max(MEAN_DEV, abs(m - 1.0))
    VAR_DEV = max(VAR_DEV, abs(v - 1.0 / N))

# the two printed sigma^2 forms of Eq. 3 are the same number: tbar^2/N == N tbar_i^2
_tb = 60.0
EQ3_VAR_IDENTITY = max(abs(_tb ** 2 / N - N * (_tb / N) ** 2) for N in [2, 4, 10, 50])

# the threshold, root-found
STIRLING_NCRIT = brentq(
    lambda N: abs(Emax_stirling(N) - Emax_printed(N)) / Emax_printed(N) - P("Fig14.3", "stirling_error_bound"),
    2.5, 40.0, xtol=1e-10)
STIRLING_ERR_N5 = float(FIG143.loc[FIG143.N == 5, "stirling_rel_err"].iloc[0])
STIRLING_ERR_N6 = float(FIG143.loc[FIG143.N == 6, "stirling_rel_err"].iloc[0])

EINF_LO_N4 = float(FIG143.loc[FIG143.N == 4, "Einf_lo"].iloc[0])
EINF_HI_N4 = float(FIG143.loc[FIG143.N == 4, "Einf_hi"].iloc[0])
EINF_LO_N10 = float(FIG143.loc[FIG143.N == 10, "Einf_lo"].iloc[0])
EINF_HI_N10 = float(FIG143.loc[FIG143.N == 10, "Einf_hi"].iloc[0])
EINF_ASYMPTOTE = np.exp(-0.5)
EINF_P4, EINF_P10 = P("Fig14.3", "Einf_over_Emax_N4"), P("Fig14.3", "Einf_over_Emax_Nge10")

display(FIG143[["N", "theta_max", "theta_max_printed", "width_ratio", "width_printed",
                "Emax_root", "Emax_printed", "stirling_rel_err"]].round(10))
display(Markdown(f"""
Over $N$ = {NS.min()} to {NS.max()}: $\\theta_{{\\max}}$ root-found against the printed
$(N-1)/N$ to **{THETA_MAX_DEV:.1e}**, $\\Delta\\theta/\\theta_{{\\max}}$ against the
printed $2/\\sqrt{{N-1}}$ to **{WIDTH_DEV:.1e}**, and $E_{{\\theta,\\max}}$ evaluated at
the root against the printed closed form to **{EMAX_DEV:.1e}** relative. By
quadrature, area $-1$ is at most **{AREA_DEV:.1e}**, mean $-1$ at most
**{MEAN_DEV:.1e}**, and $\\sigma^2_\\theta - 1/N$ at most **{VAR_DEV:.1e}**.

**The 2 % claim.** The relative error of $N/\\sqrt{{2\\pi(N-1)}}$ crosses 2 % at
$N = {STIRLING_NCRIT:.4f}$: it is {STIRLING_ERR_N5 * 100:.2f} % at $N=5$ and
{STIRLING_ERR_N6 * 100:.2f} % at $N=6$. The printed "error &lt; 2 % for $N > 5$" is
therefore correct, and correct at the tightest integer &mdash; one lower and it
would be false. (It is Stirling's series: the approximation is high by
$1/(12(N-1))$, which is 2 % at $N-1 = 4.17$.)

**The one annotation this page cannot settle.** Fig. 14.3 draws a single level
$E_{{\\theta,\\text{{inf}}}}$ and labels it $\\cong {EINF_P4:g}\\,E_{{\\theta,\\max}}$ for
$N=4$ and $\\cong {EINF_P10:g}\\,E_{{\\theta,\\max}}$ for $N\\ge 10$. The two inflection
points do **not** sit at the same height: at $N=4$ they are at
{EINF_LO_N4:.4f} and {EINF_HI_N4:.4f} of the maximum, and at $N=10$ at
{EINF_LO_N10:.4f} and {EINF_HI_N10:.4f}. The printed value lies between them in both
cases (their geometric mean is {np.sqrt(EINF_LO_N4 * EINF_HI_N4):.4f} at $N=4$,
{np.sqrt(EINF_LO_N10 * EINF_HI_N10):.4f} at $N=10$), and the $N\\ge 10$ figure is the
Gaussian limit $e^{{-1/2}} = {EINF_ASYMPTOTE:.4f}$. Which convention the figure
intends is not recoverable from the page, so **both heights are reported and
neither is scored as an agreement.**
"""))'''))

cells.append(md(r"""### 4. The RTD, three ways

The exit curve is computed by

1. **pymrm** — a pulse into tank 1 of the $N$-cell cascade, marched with
   Crank–Nicolson;
2. **Eq. 3** — the closed form, in the second of its four printed shapes;
3. **explicit convolution** — $E_M \star E_N$ evaluated by Simpson quadrature,
   which is the *definition* of two vessels in series and shares no assembly with
   either of the others.

Route 3 is what makes Eq. 4 testable rather than assumed, and it is the one used
in the additivity section below. Route 1 carries a genuine time-step error, and
it is refined.""" ))

cells.append(code(r'''NT_PROD, TMAX_PROD = 4000, 300.0


def rtd_error(N, tbar, nt, theta):
    ts, es = Cascade(N, tbar=tbar).pulse_response(TMAX_PROD, nt, theta=theta)
    return float(np.max(np.abs(es - E_time(ts, N, tbar)))), ts, es


NTS = np.array([500, 1000, 2000, 4000])
ERR_CN = np.array([rtd_error(N4, TB4, nt, 0.5)[0] for nt in NTS])
ERR_IE = np.array([rtd_error(N4, TB4, nt, 1.0)[0] for nt in NTS])
ORDER_CN = float(np.polyfit(np.log(1.0 / NTS), np.log(ERR_CN), 1)[0])
ORDER_IE = float(np.polyfit(np.log(1.0 / NTS), np.log(ERR_IE), 1)[0])
RTD_MAX_DEV, TS_P, ES_P = rtd_error(N4, TB4, NT_PROD, 0.5)

# area under the pymrm curve: 1 by physics, not by normalisation.  The window stops
# at TMAX_PROD, so the analytic tail of Eq. 3 beyond it is computed and reported
# separately - otherwise the truncation would be mistaken for solver error.
_h = TS_P[1] - TS_P[0]
RTD_AREA = float(_h / 3.0 * (ES_P[0] + ES_P[-1] + 4 * ES_P[1:-1:2].sum()
                             + 2 * ES_P[2:-1:2].sum()))
RTD_AREA_DEV = abs(RTD_AREA - 1.0)
_tail = np.linspace(TMAX_PROD, 40 * TB4, 200001)
_ht = _tail[1] - _tail[0]
_yt = E_time(_tail, N4, TB4)
RTD_TAIL = float(_ht / 3.0 * (_yt[0] + _yt[-1] + 4 * _yt[1:-1:2].sum() + 2 * _yt[2:-1:2].sum()))
RTD_AREA_RESIDUAL = abs(RTD_AREA_DEV - RTD_TAIL)


def convolve_E(t, M, N, tiM, tiN, npan=2000):
    """(E_M * E_N)(t) by composite Simpson on [0, t].  No cascade, no gamma identity."""
    out = np.empty_like(t)
    for j, tt in enumerate(t):
        if tt <= 0:
            out[j] = 0.0
            continue
        s = np.linspace(0.0, tt, npan + 1)
        y = E_time(s, M, M * tiM) * E_time(tt - s, N, N * tiN)
        h = s[1] - s[0]
        out[j] = h / 3.0 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())
    return out


_tc = np.linspace(1e-9, TMAX_PROD, 121)
CONV_EQ = convolve_E(_tc, 2, 2, TB4 / N4, TB4 / N4)
RTD_CONV_DEV = float(np.max(np.abs(CONV_EQ - E_time(_tc, N4, TB4))))
# ... and the same thing on 1/100 of the panels.  For EQUAL tank sizes the
# integrand E_M(s) E_N(t-s) collapses to s^(M-1)(t-s)^(N-1)e^(-t/t_i), a
# polynomial of degree M+N-2 = 2 times a constant-in-s exponential, and Simpson
# is EXACT on cubics.  If the two numbers below agree, the check is an identity
# and bounds nothing about the quadrature.
RTD_CONV_DEV_COARSE = float(np.max(np.abs(convolve_E(_tc, 2, 2, TB4 / N4, TB4 / N4, npan=20)
                                          - E_time(_tc, N4, TB4))))

fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.6))
for N in [1, 2, 4, 10, 30]:
    th = np.linspace(1e-6, 2.6, 800)
    ax[0].plot(th, E_theta(th, N), lw=1.4, label=f"N = {N}")
ax[0].set_xlabel(r"$\theta = t/\bar t$"); ax[0].set_ylabel(r"$E_\theta$")
ax[0].set_title("Eq. 3, dimensionless"); ax[0].legend(fontsize=7.5); ax[0].set_ylim(0, 2.4)

ax[1].plot(TS_P, E_time(TS_P, N4, TB4), "k-", lw=2.2, label="Eq. 3 (closed form)")
ax[1].plot(TS_P[::80], ES_P[::80], "o", ms=4.5, mfc="none", label="pymrm cascade (Crank-Nicolson)")
ax[1].plot(_tc[::4], CONV_EQ[::4], "x", ms=4.5, label=r"$E_2\star E_2$ (Simpson)")
ax[1].axvline(E144_TMAX, color="0.6", lw=0.8, ls=":")
ax[1].set_xlabel("t, s"); ax[1].set_ylabel(r"E, s$^{-1}$")
ax[1].set_title(f"E14.4: N = {N4}, $\\bar t$ = {TB4:g} s"); ax[1].legend(fontsize=7.5)
plt.tight_layout(); plt.show()

display(Markdown(f"""
At the production time step ($n_t = {NT_PROD}$ over {TMAX_PROD:g} s) the pymrm
cascade reproduces Eq. 3 to **{RTD_MAX_DEV:.2e}** s⁻¹, against a peak of
{E144_EMAX:.4g} s⁻¹. Refining the step: Crank-Nicolson at observed order
**{ORDER_CN:.3f}**, implicit Euler at **{ORDER_IE:.3f}** &mdash; both the expected
values, and the check is on the step because the *cell* count is the model.
The area under the pymrm curve falls **{RTD_AREA_DEV:.2e}** short of 1 with no
normalisation applied anywhere - the injected amount is one mole into $V_1$ and
$E = uC_N$ throughout - and {RTD_TAIL:.2e} of that is the analytic tail of Eq. 3
beyond $t = {TMAX_PROD:g}$ s, leaving **{RTD_AREA_RESIDUAL:.1e}** that the march is
responsible for.

The convolution route $E_2\\star E_2$ lands on the same curve to
**{RTD_CONV_DEV:.1e}** &mdash; and it does so on 2000 panels and on
{20} panels alike ({RTD_CONV_DEV_COARSE:.1e}), which is the tell. For *equal*
tank sizes the integrand is a cubic in $s$ times a constant, Simpson integrates
cubics exactly, and this number is therefore a **structural identity**, not an
independent measurement: it checks the transcription of Eq. 3 and the plumbing
of `convolve_E`, and it bounds nothing about the quadrature error on the
*unequal* chains that the next two sections actually integrate. It is labelled
structural in the coverage map for that reason, and the page's independent route
to the RTD is the Crank-Nicolson march against Eq. 3 above ({RTD_MAX_DEV:.2e}).
"""))'''))

cells.append(md(r"""### 5. Additivity, and the assumption it hides

Eq. 4 is a *structural* claim: for equal tanks it is the statement that a Gamma
distribution convolved with a Gamma distribution of the same scale is a Gamma
distribution, and the check below closes to round-off because it must. Reporting
that as evidence would be exactly the "check that cannot fail".

What the check *can* do is show what Eq. 4's parenthesis — "(all of the same
size)" — is worth. The means and variances of two **unequal** blocks still add:
that is true of any two independent stages and has nothing to do with tanks in
series. The *shape* does not follow. Below, two blocks with different tank sizes
are convolved, and the result is compared with the tanks-in-series curve carrying
the same mean and the same variance.""" ))

cells.append(code(r'''TI = TB4 / N4
# The window has to hold essentially the whole composite before the moment
# integrals below mean anything: the composite has tbar = 80 s, and a 340 s
# window truncates enough tail to dominate them (measured in the table below).
EQ4_WINDOW, EQ4_NPTS = 800.0, 331
_t = np.linspace(1e-9, EQ4_WINDOW, EQ4_NPTS)

# (a) equal tanks: structural, closes to round-off
EQ4_EQUAL = float(np.max(np.abs(convolve_E(_t, 2, 3, TI, TI) - E_time(_t, 5, 5 * TI))))

# (b) unequal blocks: 2 tanks of 10 s followed by 3 tanks of 20 s
tA, tB, MA, MB = 10.0, 20.0, 2, 3
CONV_UNEQ = convolve_E(_t, MA, MB, tA, tB)
MEAN_UNEQ = MA * tA + MB * tB
VAR_UNEQ = MA * tA ** 2 + MB * tB ** 2
N_MOMENT = MEAN_UNEQ ** 2 / VAR_UNEQ


def eq4_moment_devs(window, npts):
    """The two moment integrals on a given outer limit of the quadrature."""
    t = np.linspace(1e-9, window, npts)
    c = convolve_E(t, MA, MB, tA, tB)
    h = t[1] - t[0]
    sq = lambda y: h / 3.0 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())
    m = sq(t * c) / sq(c)
    v = sq((t - m) ** 2 * c) / sq(c)
    return abs(m - MEAN_UNEQ) / MEAN_UNEQ, abs(v - VAR_UNEQ) / VAR_UNEQ


EQ4_WIN = pd.DataFrame(
    [dict(window_s=w, npts=n, mean_reldev=eq4_moment_devs(w, n)[0],
          var_reldev=eq4_moment_devs(w, n)[1])
     for w, n in [(340.0, 141), (500.0, 207), (EQ4_WINDOW, EQ4_NPTS), (1200.0, 495)]])
EQ4_MEAN_DEV, EQ4_VAR_DEV = eq4_moment_devs(EQ4_WINDOW, EQ4_NPTS)
EQ4_MEAN_DEV_SHORT = float(EQ4_WIN.mean_reldev.iloc[0])
EQ4_VAR_DEV_SHORT = float(EQ4_WIN.var_reldev.iloc[0])

# The shape gap is an EXTREMUM, so it is bracketed and refined, not sampled.
_gap_fn = lambda t: abs(float(convolve_E(np.array([t]), MA, MB, tA, tB)[0])
                        - float(E_time(t, N_MOMENT, MEAN_UNEQ)))
CONV_UNEQ_PEAK = max_brent(
    lambda t: float(convolve_E(np.array([t]), MA, MB, tA, tB)[0]), 1e-9, EQ4_WINDOW)
EQ4_UNEQUAL_GAP = max_brent(_gap_fn, 1e-9, EQ4_WINDOW) / CONV_UNEQ_PEAK
EQ4_UNEQUAL_GAP_SAMPLED = float(np.max(np.abs(CONV_UNEQ - E_time(_t, N_MOMENT, MEAN_UNEQ)))
                                / CONV_UNEQ.max())

display(EQ4_WIN.round(10))
display(Markdown(f"""
**Equal tanks.** $E_2\\star E_3$ against $E_5$ with the same $\\bar t_i$:
**{EQ4_EQUAL:.2e}** s⁻¹. This is an identity and is labelled as one; it cannot
fail and it is not evidence for anything physical. It does not protect the
quadrature either: for equal tank sizes the integrand is a cubic and Simpson is
exact on cubics, which is why the same number comes back on 20 panels. The
quadrature is instead pinned by refining it directly, on the *unequal* integrand
where it is not exact.

**Unequal blocks.** {MA} tanks of {tA:g} s followed by {MB} tanks of {tB:g} s.
The first two moments of the composite are $\\sum t_i = {MEAN_UNEQ:g}$ s and
$\\sum t_i^2 = {VAR_UNEQ:g}$ s², and they add for *any* independent stages, exactly,
by a one-line argument &mdash; so the two deviations in the table above measure
**the outer limit of the quadrature and nothing else**. That is worth showing
rather than hiding: on the 340 s window a reader would see {EQ4_MEAN_DEV_SHORT:.1e}
and {EQ4_VAR_DEV_SHORT:.1e} and might take them for a test of additivity, whereas
they are the truncated tail of a distribution with $\\bar t = {MEAN_UNEQ:g}$ s. On
the {EQ4_WINDOW:g} s window used here they fall to {EQ4_MEAN_DEV:.1e} and
{EQ4_VAR_DEV:.1e}, and they stay there.

The *shape* is the part that can fail. The moment-matched tanks-in-series curve
($N = {N_MOMENT:.4f}$, $\\bar t = {MEAN_UNEQ:g}$ s) misses the true composite by
**{EQ4_UNEQUAL_GAP * 100:.4f} %** of its peak &mdash; bracketed and refined with
Brent's method rather than read off the plotting grid, which would have given
{EQ4_UNEQUAL_GAP_SAMPLED * 100:.4f} %. Eq. 4's parenthesis is doing work.
"""))'''))

cells.append(md(r"""### 6. Conversion: what the fitted $N$ predicts

Everything in this section is held out of the E14.4 fit, which used two moments
and nothing else. Three routes to the first-order exit concentration, sharing no
code:

- **pymrm** — one sparse factorisation of the cascade with a linear sink;
- **Eq. 9** — $\left(1+k\bar t/N\right)^{-N}$;
- **Eq. 11** — the macrofluid integral $\int_0^\infty e^{-kt}E(t)\,dt$ by Simpson
  quadrature, which is the Laplace transform of Eq. 3 computed numerically.

For a first-order reaction the last two *must* agree, because Levenspiel writes
"These equations apply to both micro- and macrofluids" — micromixing cannot
matter when the rate is linear — so route 3 is a real test of Eq. 11's
transcription, of the assumed $\bar t^{\,N}$ denominator, and of the quadrature,
against an algebraic answer none of them can see.

Second order splits the two: Eq. 10's nest of radicals is the microfluid,
Eq. 11 with the batch solution $C_0/(1+kC_0t)$ is the macrofluid, and they are
different numbers.""" ))

cells.append(code(r'''KTS = np.array([0.25, 0.5, 1.0, 2.0, 5.0, 10.0])
NN = [1, 2, 3, 4, 6, 10, 20, 40]

rows = []
for N in NN:
    for kt in KTS:
        pym, _ = Cascade(N, tbar=1.0, c_in=1.0).first_order(kt)
        alg = eq9(N, kt)
        mac = eq11(N, 1.0, lambda t, kk=kt: np.exp(-kk * t))
        rows.append(dict(N=N, ktbar=kt, pymrm=pym, eq9=alg, eq11=mac))
FO = pd.DataFrame(rows)
FO["pymrm_dev"] = (FO.pymrm - FO.eq9).abs()
FO["eq11_dev"] = (FO.eq11 - FO.eq9).abs()
PYMRM_EQ9_DEV = float(FO.pymrm_dev.max())
EQ11_EQ9_DEV = float(FO.eq11_dev.max())
EQ11_NORM_DEV = max(abs(eq11(N, 1.0, lambda t: np.ones_like(t)) - 1.0) for N in NN)

# Simpson refinement of Eq. 11 - the quadrature IS the refinable axis here
_pans = np.array([2500, 5000, 10000, 20000])
_qerr = np.array([abs(eq11(4, 1.0, lambda t: np.exp(-2.0 * t), npanel=p) - eq9(4, 2.0))
                  for p in _pans])
EQ11_QUAD_ORDER = float(np.polyfit(np.log(1.0 / _pans), np.log(_qerr), 1)[0])

# second order: pymrm Newton vs Eq. 10, and macro vs micro
rows = []
for N in NN:
    for kt in KTS:
        pym, _, sol = Cascade(N, tbar=1.0, c_in=1.0).order_n(kt, order=2.0)
        alg = eq10(N, kt, 1.0 / N, 1.0)
        mac = eq11(N, 1.0, lambda t, kk=kt: 1.0 / (1.0 + kk * t))
        rows.append(dict(N=N, kC0tbar=kt, pymrm=pym, eq10=alg, eq11_macro=mac,
                         nit=sol.nit))
SO = pd.DataFrame(rows)
SO["pymrm_dev"] = (SO.pymrm - SO.eq10).abs()
SO["macro_micro_gap"] = ((SO.eq11_macro - SO.eq10) / SO.eq10).abs()
PYMRM_EQ10_DEV = float(SO.pymrm_dev.max())
MACRO_MICRO_MAX = float(SO.macro_micro_gap.max())
MACRO_MICRO_ROW = SO.loc[SO.macro_micro_gap.idxmax()]
NEWTON_MAX_IT = int(SO.nit.max())

# the printed small-deviation asymptote, for equal volume
_ns = np.array([10, 20, 40, 80, 160, 320])
_kt = 1.0
_asym_err = np.array([abs((eq9(n, _kt) / np.exp(-_kt) - 1.0) / (_kt ** 2 / (2 * n)) - 1.0)
                      for n in _ns])
ASYM_ORDER = float(-np.polyfit(np.log(_ns), np.log(_asym_err), 1)[0])
ASYM_DEV_N320 = float(_asym_err[-1])

display(Markdown(f"""
**First order, three routes.** Over {len(NN)} values of $N$ and
{len(KTS)} values of $k\\bar t$ (to $k\\bar t = {KTS.max():g}$):
pymrm against Eq. 9 to **{PYMRM_EQ9_DEV:.1e}**, and the Eq. 11 macrofluid
quadrature against Eq. 9 to **{EQ11_EQ9_DEV:.1e}**. Eq. 11 with a batch solution
of 1 integrates to 1 to **{EQ11_NORM_DEV:.1e}**, which is the check that fixes the
$\\bar t^{{\\,N}}$ reading of its denominator. The Simpson error falls at observed
order **{EQ11_QUAD_ORDER:.2f}**.

**Second order.** pymrm's Newton solve (pointwise `NumJac`, at most
{NEWTON_MAX_IT} iterations) against Eq. 10's nest of $N$ radicals to
**{PYMRM_EQ10_DEV:.1e}**. Micro- and macrofluid are *not* the same here: the gap
reaches **{MACRO_MICRO_MAX * 100:.2f} %** at $N = {int(MACRO_MICRO_ROW.N)}$,
$kC_0\\bar t = {MACRO_MICRO_ROW.kC0tbar:g}$, and it is the macrofluid that converts
further. This is the same statement [`A2.8`](../A2.8-zwietering-segregation/)
makes with bounds; here it is the two specific states Levenspiel prints.

**The printed asymptote.** $C_{{A,N}}/C_{{Ap}} = 1 + (k\\bar t)^2/(2N)$ for the same
volume: the relative error of that expression falls at observed order
**{ASYM_ORDER:.2f}** in $N$ and is {ASYM_DEV_N320:.2e} at $N = {_ns[-1]}$, so the
printed formula is the correct first term. Its printed companion for the same
final conversion is *not*; see the next section.
"""))'''))

# -------------------------------------------------------- printed defects ---
cells.append(md(r"""### 7. Printed defects

Six items are wrong or ambiguous as printed — five defects and one typesetting
ambiguity. Each is established from the book's own numbers, and none is repaired
in the transcription: the CSV carries the glyphs as they appear and this section
says what they should be.

The first of them is the most consequential, because it is a **modelling rule**
and not a unit or a tick label: a reader who applies it gets a different reactor.
It is also the one this page nearly missed — an earlier draft asserted that the
chapter printed no such rule at all.""" ))

cells.append(code(r'''# (0) Fig. 14.7's dispersion-conversion rule, inverted.  Everything here is
# arithmetic on printed values; the COST of the inversion is measured in
# section 9, where a dispersion solve exists.
HALF = P("Fig14.7", "dispersion_rule_half")
N_SYM = P("Fig14.7", "N_symmetric_threshold")
SIG_COEF = P("Fig13.12", "sigma2_theta_over_DuL")     # sigma^2_theta = 2 (D/uL)
DUL_LIMIT = P("Fig13.11", "DuL_small_deviation_limit")
# Ch. 13's variance with Ch. 14's 1/N: 1/N = SIG_COEF (D/uL) -> N = (1/SIG_COEF)(uL/D)
DEFECT_FIG147_FROM_CH13 = abs(1.0 / SIG_COEF - HALF)
# the annotation's own threshold, pushed through the CORRECT rule
DEFECT_FIG147_THRESHOLD_DUL = HALF / N_SYM
DEFECT_FIG147_THRESHOLD_DEV = abs(DEFECT_FIG147_THRESHOLD_DUL - DUL_LIMIT)
# ... and through the rule AS PRINTED
DEFECT_FIG147_THRESHOLD_DUL_PRINTED = N_SYM / HALF

# (1) the factor of two in "for same C_A,final".
# The quantity is a LIMIT, x -> 0, of (V_N/V_p - 1)/x with x = k tbar_i, so it is
# extrapolated rather than fitted over an arbitrary window: a polyfit on a finite
# window returns 2.0003 or 2.03 depending only on how wide the window is, and the
# trailing digits would then be pinned in CI while meaning nothing.
def _vslope(x):
    """(V_N/V_p - 1)/x from the book's own Eq. 9 against plug flow."""
    return (x / np.log1p(x) - 1.0) / x


DEFECT_SLOPE_TABLE = pd.DataFrame(
    [dict(x=x, raw=_vslope(x), richardson=2 * _vslope(x / 2) - _vslope(x),
          polyfit_ratio=1.0 / float(np.polyfit(np.linspace(1e-6, x, 400),
                                               np.linspace(1e-6, x, 400)
                                               / np.log1p(np.linspace(1e-6, x, 400)) - 1.0,
                                               1)[0]))
     for x in [4e-1, 4e-2, 4e-3, 4e-4]])
DEFECT_EXACT_SLOPE = float(2 * _vslope(5e-5) - _vslope(1e-4))   # -> 1/2 - x^2/48
DEFECT_MID_SLOPE = 1.0                        # "1 + k tbar_i"       -> slope 1
DEFECT_RIGHT_SLOPE = 0.5                      # "1 + k tbar/(2N)"    -> slope 1/2
DEFECT_FACTOR = DEFECT_MID_SLOPE / DEFECT_EXACT_SLOPE
DEFECT_RIGHT_RATIO = DEFECT_RIGHT_SLOPE / DEFECT_EXACT_SLOPE
# direct numerical confirmation at a finite N: what volume ratio does Eq. 9 need?
_N, _ktb = 20, 1.0
_target = np.exp(-_ktb)                       # plug-flow conversion at k tbar = 1
_ratio = brentq(lambda r: eq9(_N, _ktb * r) - _target, 1.0, 5.0, xtol=1e-14)
DEFECT_RATIO_TRUE = _ratio
DEFECT_RATIO_MID = 1.0 + _ktb / _N
DEFECT_RATIO_RIGHT = 1.0 + _ktb / (2 * _N)

# (2) Eq. 11's denominator: tbar_N (subscript) against tbar**N (superscript)
_N11, _tb11 = 4, 60.0
NORM_SUPER = eq11(_N11, _tb11, lambda t: np.ones_like(t))
NORM_SUB = NORM_SUPER * _tb11 ** _N11 / _tb11          # what the subscript reading gives
DEFECT_EQ11_NORM_SUB = NORM_SUB

# (3) Fig. E14.4b's axis tick, against the peak of the curve AS PRINTED
DEFECT_TICK_RATIO = TICK_P / E144_EMAX_PRINTED

# (4) E14.2's "moles";  (5) E14.4's "900 s"
DEFECT_VAR_UNIT_VALUE_DEV = abs(DV - DV_P)     # the VALUE is right; the unit is not

display(DEFECT_SLOPE_TABLE.round(10))
display(Markdown(f"""
**(a) Fig. 14.7, book p. 327 — the dispersion-conversion rule is inverted.**
Inside the axes box: *"When $N > {N_SYM:g}$ the curve becomes symmetrical in which
case use fig. 13-11 with $N = {HALF:g}(\\mathbf{{D}}/uL)$"*. The ratio is upside
down, and the book supplies the proof twice over.

*From its own variance.* Chapter 13 prints
$\\sigma^2_\\theta = {SIG_COEF:g}(\\mathbf{{D}}/uL)$ (Fig. 13.12, book p. 303) and
Chapter 14 prints $\\sigma^2_\\theta = 1/N$ (Eq. 3). Together they give
$N = (1/{SIG_COEF:g})(uL/\\mathbf{{D}})$, whose coefficient matches the printed
{HALF:g} to **{DEFECT_FIG147_FROM_CH13:.0e}** &mdash; on the *reciprocal* ratio.

*From its own threshold.* Pushed through the corrected rule, the annotation's
$N > {N_SYM:g}$ is $\\mathbf{{D}}/uL < {DEFECT_FIG147_THRESHOLD_DUL:g}$, and
"Small Deviation from Plug Flow, $\\mathbf{{D}}/uL < {DUL_LIMIT:g}$" is the printed
heading of the section that owns Fig. 13.11: the two agree to
**{DEFECT_FIG147_THRESHOLD_DEV:.0e}**. As printed, the same threshold would send
the reader to $\\mathbf{{D}}/uL > {DEFECT_FIG147_THRESHOLD_DUL_PRINTED:g}$, off a
figure drawn between 0.00005 and 0.0128 and at the opposite end of the flow
regime. Section 9 measures what that costs in conversion.

**(b) Section 14.2, book p. 328 — a factor of two.** The line reads
$V_{{N\\,\\text{{tanks}}}}/V_p = 1 + k\\bar t_i = 1 + k\\bar t/(2N)$. With the book's own
$\\bar t = N\\bar t_i$ those two right-hand sides differ by exactly two. Expanding the
book's own Eq. 9 against plug flow gives $V_N/V_p = x/\\ln(1+x)$ with
$x = k\\bar t_i$, whose leading correction has slope
**{DEFECT_EXACT_SLOPE:.9f}** &mdash; the table above extrapolates it to $x = 0$
instead of fitting it, because the `raw` and `polyfit_ratio` columns show that a
finite window puts spurious digits on an exactly rational number. Against that,
the middle expression has slope {DEFECT_MID_SLOPE:g} (a factor
**{DEFECT_FACTOR:.9f}** too large, i.e. exactly two) and the right-hand one
{DEFECT_RIGHT_SLOPE:g} (factor {DEFECT_RIGHT_RATIO:.9f}, i.e. exactly one). At
$N = {_N}$, $k\\bar t = {_ktb:g}$ the volume ratio that Eq. 9 actually needs is
**{DEFECT_RATIO_TRUE:.6f}**, against {DEFECT_RATIO_MID:.6f} from the middle form
and {DEFECT_RATIO_RIGHT:.6f} from the right-hand one. **The middle expression is
missing a factor of $\\tfrac12$**; the right-hand form, and the companion
equal-volume formula $1+(k\\bar t)^2/(2N)$, are correct as printed. Reported, not
repaired.

**(c) Eq. 11, book p. 329 — a subscript where a superscript belongs.** The
prefactor denominator is set as $(N-1)!\\,\\bar t_N$ with $N$ lowered. With
$\\bar t^{{\\,N}}$ the kernel integrates to {NORM_SUPER:.10f}; with the printed
$\\bar t_N = \\bar t$ it integrates to {DEFECT_EQ11_NORM_SUB:.4g} at $N = {_N11}$,
$\\bar t = {_tb11:g}$ s &mdash; not dimensionless, let alone normalised. The book itself
prints the raised form two pages later in E14.4
($\\mathbf{{E}} = t^{{N-1}}/\\bar t^{{\\,N}}\\cdot N^N/(N-1)!\\cdot e^{{-tN/\\bar t}}$), which
settles the reading. A typesetting slip, not a modelling one.

**(d) Fig. E14.4b, book p. 334 — a transposed digit.** The top axis tick is
labelled {TICK_P:g} s⁻¹, on an axis whose other ticks are 0, 0.005 and 0.010 at
equal spacing. The book's own printed $E(t) = {PRE_P:.4e}\\,t^3e^{{-{EXP_P:g}t}}$
peaks at $t = {E144_TMAX_PRINTED:.6g}$ s with $E = {E144_EMAX_PRINTED:.7g}$ s⁻¹, so
the tick is {DEFECT_TICK_RATIO:.2f}$\\times$ the whole curve. It should read 0.015.
(Note which curve: with the *unrounded* coefficients the peak is
{E144_EMAX:.7g} s⁻¹ and the ratio {TICK_P / E144_EMAX:.2f}. The tick belongs to the
printed curve, so the printed coefficients are the ones used here.)

**(e) E14.2, book p. 331 — "$L = 272$ moles".** The quantity is a distance along
the Ohio River. The Comment two lines below reads "any location where
$L \\le 272$ **miles** is suspect", and every input to the calculation is in miles.
The *value* {E142_L_P:g} is exactly right, as the reproduction above shows; only
the unit is a slip.

**(f) E14.4, book p. 333 — "$\\Delta(\\sigma^2) = 1000 - 100 = 900$ s".** The same
class of slip as (e), in the same chapter: a variance carrying the unit of a
time. The inputs two lines above are read off Fig. E14.4a as {V_IN:g} s² and
{V_OUT:g} s², and the very next line divides $(\\Delta\\bar t)^2$ by this quantity to
get a dimensionless tank count, which only works in s². The *value* {DV_P:g} is
right &mdash; the reproduction above matches it to {DEFECT_VAR_UNIT_VALUE_DEV:.0e}
&mdash; and the CSV's `note` column carries the correction while `printed_as`
carries the glyphs.
"""))'''))

# ------------------------------------------------------- the headline -------
cells.append(md(r"""### 8. What two moments cannot buy: the headline

E14.4 uses exactly two numbers, $\Delta\bar t$ and $\Delta(\sigma^2)$, and
returns "4 tanks". The method assumes the vessel is $N$ **equal** tanks; with two
moments there is no way to check that assumption, because a chain of unequal
tanks has the same two moments as some equal chain by construction.

Here is a specific alternative. A five-stage chain with residence times
$9, 9, 9, 9, 24$ s has

$$\sum t_i = 60\ \text{s} = \Delta\bar t,\qquad
  \sum t_i^2 = 4(81) + 576 = 900\ \text{s}^2 = \Delta(\sigma^2),$$

**exactly** the two numbers E14.4 measured — the same mean, the same variance, to
the last digit. It is therefore indistinguishable from Levenspiel's four equal
tanks by his own method, and every gallery page that reuses a fitted $N$ inherits
the ambiguity. It is also *admissible*: E14.3's own Figs. E14.3b–d (book p. 333)
work with $V_1 = \tfrac13$, $V_2 = \tfrac23$ and with $N_1 = 12$, $N_2 = 24$, so
unequal blocks are not excluded by the chapter. And it is the natural choice —
9, 9, 9, 9, 24 is the second root of $a^2 - 24a + 135 = 0$, i.e. the *only*
non-equal member of the four-equal-plus-one family.

Three things have to be said carefully about the number that comes out, because
each of them is a way of overstating it:

- it is a gap in **outlet concentration**, which is not the same as a gap in
  conversion, and both are printed below;
- it is the value at the **right-hand edge** of the range swept, and it is
  monotone in $k\bar t$, so the sweep is extended past that edge rather than
  stopped at it;
- it is **one chain**, so it is a *lower* bound on what two moments leave open.
  The upper bound over the whole family is computed in closed form afterwards.

Both chains are built with the same `Cascade` class, differing only in `faces`.""" ))

cells.append(code(r'''T_EQ = np.full(4, 15.0)
T_UN = np.array([9.0, 9.0, 9.0, 9.0, 24.0])
assert np.isclose(T_UN.sum(), DTB) and np.isclose((T_UN ** 2).sum(), DV), "moments must match exactly"
AMBIG_MEAN_DEV = abs(T_UN.sum() - T_EQ.sum())
AMBIG_VAR_DEV = abs((T_UN ** 2).sum() - (T_EQ ** 2).sum())


def chain(ts, tbar_total):
    """A Cascade whose cells are sized in proportion to the given residence times."""
    faces = np.concatenate([[0.0], np.cumsum(ts) / ts.sum()])
    return Cascade(len(ts), tbar=tbar_total, c_in=1.0, faces=faces)


rows = []
for kt in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
    k = kt / DTB
    a1, _ = chain(T_EQ, DTB).first_order(k)
    b1, _ = chain(T_UN, DTB).first_order(k)
    a2, _, _ = chain(T_EQ, DTB).order_n(k, order=2.0)
    b2, _, _ = chain(T_UN, DTB).order_n(k, order=2.0)
    rows.append(dict(ktbar=kt, eq_first=a1, un_first=b1,
                     conc_gap=abs(b1 - a1) / a1,
                     eq_conv=1.0 - a1, un_conv=1.0 - b1,
                     conv_gap=abs((1.0 - b1) - (1.0 - a1)) / (1.0 - a1),
                     eq_second=a2, un_second=b2, second_gap=abs(b2 - a2) / a2))
AMB = pd.DataFrame(rows)
AMBIG_KT_HEADLINE = 10.0                      # the right-hand edge of the ORIGINAL sweep
_hl = AMB.loc[AMB.ktbar == AMBIG_KT_HEADLINE].iloc[0]
AMBIG_FIRST_MAX = float(_hl.conc_gap)         # the 5-stage chain AT k tbar = 10
AMBIG_FIRST_MAX_SWEPT = float(AMB.conc_gap.max())      # ... and at the new edge, 50
AMBIG_CONV_MAX = float(AMB.conv_gap.max())    # the same thing in CONVERSION
AMBIG_CONV_KT = float(AMB.loc[AMB.conv_gap.idxmax(), "ktbar"])
AMBIG_CONV_HL = float(_hl.conv_gap)
AMBIG_SECOND_MAX = float(AMB.loc[AMB.ktbar == AMBIG_KT_HEADLINE, "second_gap"].iloc[0])
AMBIG_FIRST_KT2 = float(AMB.loc[AMB.ktbar == 2.0, "conc_gap"].iloc[0])
# the pymrm unequal chain against the product formula it must satisfy
_k = 2.0 / DTB
_prod = 1.0 / np.prod(1.0 + _k * T_UN)
AMBIG_PYMRM_DEV = abs(chain(T_UN, DTB).first_order(_k)[0] - _prod)

_t = np.linspace(1e-9, 260.0, 400)
E_UN = convolve_E(_t, 4, 1, 9.0, 24.0)
# the RTD gap is an extremum too: bracketed and refined, and the denominator is
# the root-found peak of Eq. 3 at N = 4, tbar = 60 s computed in section 2.
AMBIG_RTD_GAP = max_brent(
    lambda t: abs(float(convolve_E(np.array([t]), 4, 1, 9.0, 24.0)[0])
                  - float(E_time(t, 4, DTB))), 1e-9, 260.0) / E144_EMAX

fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.6))
ax[0].plot(_t, E_time(_t, N4, DTB), lw=2.0, label=r"4 equal tanks of 15 s (Levenspiel)")
ax[0].plot(_t, E_UN, lw=2.0, ls="--", label=r"5 tanks of 9,9,9,9,24 s")
ax[0].set_xlabel("t, s"); ax[0].set_ylabel(r"E, s$^{-1}$")
ax[0].set_title(r"same $\bar t$ = 60 s, same $\sigma^2$ = 900 s$^2$")
ax[0].legend(fontsize=7.5)
ax[1].semilogy(AMB.ktbar, AMB.conc_gap * 100, "o-", label="first order, exit concentration")
ax[1].semilogy(AMB.ktbar, AMB.second_gap * 100, "s--", label="second order, exit concentration")
ax[1].semilogy(AMB.ktbar, AMB.conv_gap * 100, "^:", label="first order, CONVERSION")
ax[1].set_xlabel(r"$k\bar t$  (or $kC_0\bar t$)")
ax[1].set_ylabel("relative difference, %")
ax[1].set_title("what the two-moment fit does not fix"); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.show()
display(AMB.round(6))
display(Markdown(f"""
The two chains match E14.4's measured mean to **{AMBIG_MEAN_DEV:.0e}** s and its
measured variance to **{AMBIG_VAR_DEV:.0e}** s² &mdash; that is, exactly. Their exit
curves differ by **{AMBIG_RTD_GAP * 100:.1f} %** of the peak.

**In exit concentration** the two chains differ by
**{AMBIG_FIRST_MAX * 100:.1f} %** at $k\\bar t = {AMBIG_KT_HEADLINE:g}$, first order.
That is the number this page quotes, and the table shows what it is and is not:
it is monotone in $k\\bar t$, so it is not a maximum over the sweep but the value
at whichever $k\\bar t$ one stops at &mdash; {AMBIG_FIRST_KT2 * 100:.2f} % at
$k\\bar t = 2$, {AMB.conc_gap.iloc[-1] * 100:.1f} % at
$k\\bar t = {AMB.ktbar.iloc[-1]:g}$. Second order at
$k\\bar t = {AMBIG_KT_HEADLINE:g}$: {AMBIG_SECOND_MAX * 100:.1f} %.

**In conversion** it is a far smaller number, and saying "converts 23 %
differently" would be wrong. At $k\\bar t = {AMBIG_KT_HEADLINE:g}$ the two chains
convert {_hl.eq_conv * 100:.3f} % and {_hl.un_conv * 100:.3f} % &mdash;
{abs(_hl.un_conv - _hl.eq_conv) * 100:.2f} percentage points apart, or
{AMBIG_CONV_HL * 100:.2f} % relative &mdash; and over the whole sweep the relative gap in
conversion never exceeds **{AMBIG_CONV_MAX * 100:.2f} %**, at
$k\\bar t = {AMBIG_CONV_KT:g}$. The two statements are consistent: at high
$k\\bar t$ both chains convert almost everything, so a large relative difference in
what is *left* is a small difference in what is *converted*. Which of the two
matters is the reader's question, not the page's &mdash; for a purification duty,
or for anything sized on the outlet, it is the first.

The pymrm unequal chain reproduces the staged product $\\prod (1+kt_i)^{{-1}}$ to
{AMBIG_PYMRM_DEV:.1e}, so the difference is physics and not a solver artefact.

This is not a criticism of E14.4 &mdash; the method is exactly as good as its
inputs, and Levenspiel says plainly that the tanks are assumed equal. It is a
statement of what a two-moment fit is worth: **it fixes the RTD only within the
family it assumes, and the conversion it implies inherits that ambiguity.**
""" ))'''))

cells.append(md(r"""### 8b. The whole family, in closed form

One chain is an example. The honest question is the **range** of outlet
concentrations consistent with $\sum t_i = 60$ s and $\sum t_i^2 = 900$ s²,
because that is the set of vessels E14.4's two numbers cannot tell apart.

The answer is not a search. For a first-order reaction
$C/C_0 = \prod (1+kt_i)^{-1}$, so maximising the ambiguity means extremising
$\sum \ln(1+kt_i)$ under two linear-and-quadratic constraints. The Lagrange
condition is $k/(1+kt_i) = \lambda + 2\mu t_i$ for every $i$, and that equation
has **at most two roots** in $t$: the left side is increasing and concave, the
right side is a straight line. So *every* extremal chain consists of at most two
distinct tank sizes, and the extremum can be enumerated over
$(n, m)$ — $m$ tanks of one size, $n-m$ of another, with the two moments fixing
both sizes in closed form.

The supremum itself is closed form. Write $T = \sum t_i$, $S = \sum t_i^2$, and
note that $t_i \le \sqrt S$ for every $i$. With $h(t) = \ln(1+kt) - kt$, the
function $h(t)/t^2$ is increasing on $t>0$ (its numerator derivative reduces to
$u^2/(1+u)^2 \ge 0$ with $u = kt$), so $h(t_i) \le h(\sqrt S)\,t_i^2/S$ for every
$i$, and summing gives

$$\sum \ln(1+kt_i)\ \le\ \ln\!\left(1+k\sqrt S\right) + k\left(T - \sqrt S\right),$$

with equality only when every $t_i$ is $0$ or $\sqrt S$. The bound is therefore
attained only in the limit of **one tank of $\sqrt S = 30$ s followed by plug
flow of $T - \sqrt S = 30$ s** — a chain of infinitely many infinitesimal tanks.
That is the most-converting vessel E14.4's two numbers permit, and no chain
reaches it.""" ))

cells.append(code(r'''T_MOM, S_MOM = float(DTB), float(DV)          # the two numbers E14.4 measured
N_FAMILY_MAX = 400


def two_value_chains(n, m, T=T_MOM, S=S_MOM):
    """Residence times of a chain of m tanks of one size and n-m of another that
    matches BOTH printed moments.  Closed form: the constraints reduce to a
    quadratic in the first size."""
    p, q = m, n - m
    A = p * q + p * p
    B = -2.0 * T * p
    C = T * T - q * S
    disc = B * B - 4.0 * A * C
    if disc < 0.0:
        return []
    out = []
    for a in ((-B + np.sqrt(disc)) / (2.0 * A), (-B - np.sqrt(disc)) / (2.0 * A)):
        b = (T - p * a) / q
        if a > 0.0 and b > 0.0:
            out.append(np.concatenate([np.full(p, a), np.full(q, b)]))
    return out


def conc_of(ts, k):
    """prod (1 + k t_i)^-1, in logs so that 400-stage chains do not underflow."""
    return float(np.exp(-np.sum(np.log1p(k * np.asarray(ts)))))


def family_sup(k, T=T_MOM, S=S_MOM):
    """The closed-form supremum of the conversion: one sqrt(S) tank + plug flow."""
    return np.exp(-k * (T - np.sqrt(S))) / (1.0 + k * np.sqrt(S))


rows, FAM_CONV = [], {}
for kt in [2.0, AMBIG_KT_HEADLINE]:
    k = kt / T_MOM
    base = conc_of(T_EQ, k)
    lo, hi, arg = base, base, None
    conv = []
    for n in range(5, N_FAMILY_MAX + 1):
        best_n = 0.0
        for m in range(1, n):
            for ts in two_value_chains(n, m):
                c = conc_of(ts, k)
                if c < lo:
                    lo, arg = c, (n, m, ts[0], ts[-1])
                hi = max(hi, c)
                best_n = max(best_n, abs(c - base) / base)
        if n in (5, 6, 8, 12, 40, 100, N_FAMILY_MAX):
            conv.append(dict(stages=n, max_conc_gap=best_n))
    sup = family_sup(k)
    rows.append(dict(ktbar=kt, equal_chain=base,
                     most_converting=lo, least_converting=hi,
                     gap_below=abs(lo - base) / base, gap_above=abs(hi - base) / base,
                     closed_form_sup=sup, sup_gap=abs(sup - base) / base,
                     best_n=arg[0], best_a=arg[2], best_b=arg[3]))
    FAM_CONV[kt] = pd.DataFrame(conv)
FAM = pd.DataFrame(rows)

# the closed form must BOUND every enumerated chain - asserted, not assumed
for _r in rows:
    assert _r["most_converting"] >= family_sup(_r["ktbar"] / T_MOM) - 1e-14, _r

AMBIG_FAMILY_SUP_GAP = float(FAM.loc[FAM.ktbar == AMBIG_KT_HEADLINE, "sup_gap"].iloc[0])
AMBIG_FAMILY_SUP_GAP_KT2 = float(FAM.loc[FAM.ktbar == 2.0, "sup_gap"].iloc[0])
AMBIG_FAMILY_BEST_N = float(FAM.loc[FAM.ktbar == AMBIG_KT_HEADLINE,
                                    "gap_below"].iloc[0])
AMBIG_FAMILY_ABOVE = float(FAM.gap_above.max())

display(FAM.round(8))
display(FAM_CONV[AMBIG_KT_HEADLINE].round(6))
display(Markdown(f"""
At $k\\bar t = {AMBIG_KT_HEADLINE:g}$ the five-stage chain's
{AMBIG_FIRST_MAX * 100:.1f} % is a **lower bound**. Enumerating every two-value
chain up to {N_FAMILY_MAX} stages reaches {AMBIG_FAMILY_BEST_N * 100:.1f} %, and the
closed-form supremum &mdash; one {np.sqrt(S_MOM):g} s tank followed by
{T_MOM - np.sqrt(S_MOM):g} s of plug flow &mdash; is
**{AMBIG_FAMILY_SUP_GAP * 100:.1f} %**. At $k\\bar t = 2$ the same three numbers are
{AMBIG_FIRST_KT2 * 100:.2f} %, {float(FAM.loc[FAM.ktbar == 2.0, "gap_below"].iloc[0]) * 100:.2f} %
and {AMBIG_FAMILY_SUP_GAP_KT2 * 100:.2f} %. The ambiguity is **one-sided**: the
enumeration finds no chain at all that converts *less* than four equal tanks
(largest gap in that direction, {AMBIG_FAMILY_ABOVE * 100:.2f} %). Four equal
tanks already sit at the smallest $\\sum t_i^2$ their stage count allows, so every
other admissible chain must put more residence time into one stage and split the
rest more finely &mdash; and the enumeration says that always converts further.
Levenspiel's answer is therefore the *least*-converting vessel consistent with
his own two numbers, which is the conservative end and worth knowing.

So the correct statement of the headline is: E14.4's two numbers pin the outlet
concentration only to within a factor of
{1.0 / (1.0 - AMBIG_FAMILY_SUP_GAP):.2f} at $k\\bar t = {AMBIG_KT_HEADLINE:g}$; the
9, 9, 9, 9, 24 s chain realises {AMBIG_FIRST_MAX * 100:.1f} % of that ambiguity
with five tanks and is quoted because it is concrete, not because it is extreme.
"""))'''))

# ------------------------------------------------------ dispersion link -----
cells.append(md(r"""### 9. "Both models give identical results" — quantified

Chapter 14 opens by asserting that the tanks-in-series and dispersion models
agree "for not too large a deviation from plug flow ... for all practical
purposes". It also prints a rule for turning $N$ into a dispersion number — once,
inside the axes box of Fig. 14.7 on book p. 327:

> "When $N > 50$ the curve becomes symmetrical in which case use fig. 13-11
> with $N = \tfrac12\left(\dfrac{\mathbf{D}}{uL}\right)$"

**and it is inverted.** Three of the book's own printed values settle it, two of
them from the very figure the annotation points at:

1. Chapter 13 prints $\sigma^2_\theta = 2(\mathbf{D}/uL)$ for small deviations
   from plug flow (Fig. 13.12, book p. 303) and Chapter 14 prints
   $\sigma^2_\theta = 1/N$ (Eq. 3, book p. 323). Equating them gives
   $\mathbf{D}/uL = 1/(2N)$, i.e. $N = \tfrac12(uL/\mathbf{D})$ — the reciprocal
   of what is printed.
2. The annotation's own threshold. Under the correct form, $N > 50$ is exactly
   $\mathbf{D}/uL < 0.01$ — and "Small Deviation from Plug Flow,
   $\mathbf{D}/uL < 0.01$" is the printed heading of the section that owns
   Fig. 13.11. The threshold and the limit are the same statement, which they
   cannot be under the printed form.
3. Fig. 13.11 is drawn for $\mathbf{D}/uL$ between 0.00005 and 0.0128. As
   printed, a reader at $N > 50$ would enter it at $\mathbf{D}/uL > 100$ — four
   orders of magnitude off the axis, and at the mixed-flow end rather than the
   plug-flow end the annotation is about.

The derivation below reaches the same rule from a completely different direction,
the truncation error of the upwind scheme, and this section then measures what
the printed form costs. This is reported as printed defect (a), not repaired.

First-order upwind differencing of $u\,\partial C/\partial z$ on cells of width
$\Delta z$ has a leading truncation error $\tfrac12 u\,\Delta z\,\partial^2C/\partial z^2$
— it behaves like a diffusivity $D_{\text{num}} = u\Delta z/2$. With
$\Delta z = L/N$ that is $uL/(2N)$, so

$$\mathrm{Pe} = \frac{uL}{D_{\text{num}}} = 2N ,$$

which is $N = \tfrac12(uL/\mathbf{D})$ — Fig. 14.7's annotation the right way up.

The dispersion model at $\mathrm{Pe} = 2N$ is solved here with pymrm on a **real**
grid — Danckwerts inlet, zero-gradient outlet, van Leer deferred correction, the
construction published on [`A2.1`](../A2.1-danckwerts-boundary-conditions/) — so
this section, unlike everything above it, does have a grid-refinement axis, and
it is measured.

Then the question is turned around: rather than trusting $\mathrm{Pe} = 2N$, the
Péclet number that *reproduces* the $N$-tank conversion is root-found, and its
drift away from $2N$ is what "for all practical purposes" is worth.""" ))

cells.append(code(r'''def dispersion_exit(Pe, ktbar, n=1500, limiter=vanleer, max_it=200, tol=1e-13):
    """Closed vessel, Danckwerts inlet, first-order reaction.  u = 1, L = 1, tbar = 1."""
    shape = (n, 1)
    z_f = np.linspace(0.0, 1.0, n + 1)
    z_c = 0.5 * (z_f[:-1] + z_f[1:])
    D = 1.0 / Pe
    # inlet: u c* = u c - D dc/dz.  n = -z so dc/dn = -dc/dz  ->  a = D/u, b = 1, d = 1
    # outlet: dc/dz = 0.  n = +z  ->  a = 1, b = 0, d = 0
    bc = ({"a": D, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
    conv, conv_bc = construct_convflux_upwind(shape, z_f, z_c, bc, v=1.0)
    grad, grad_bc = construct_grad(shape, z_f, z_c, bc)
    div = construct_div(shape, z_f, nu=0)                     # nu=0: Cartesian
    A = (div @ (conv - D * grad) + ktbar * eye_array(n, format="csc")).tocsc()
    b0 = np.asarray((div @ (conv_bc - D * grad_bc)).todense()).ravel()
    lu = splu(A)
    c = lu.solve(-b0)
    if limiter is not None:
        for _ in range(max_it):
            _, dc = interp_cntr_to_stagg_tvd(c.reshape(shape), z_f, z_c, bc, 1.0,
                                             tvd_limiter=limiter, axis=0)
            cn = lu.solve(-b0 - np.asarray(div @ dc.reshape(-1, 1)).ravel())
            done = np.max(np.abs(cn - c)) < tol
            c = cn
            if done:
                break
        else:
            raise RuntimeError("deferred correction hit its iteration cap")
    return float(compute_boundary_values(c.reshape(shape), z_f, z_c, bc)[2].ravel()[0])


# grid refinement of the dispersion solve (this one IS a discretisation)
_ns = np.array([200, 400, 800, 1600, 3200])
_ex = np.array([dispersion_exit(8.0, 2.0, n=n) for n in _ns])
DISP_GRID_ORDER = float(-np.polyfit(np.log(_ns[:-1]), np.log(np.abs(_ex[:-1] - _ex[-1])), 1)[0])
_exu = np.array([dispersion_exit(8.0, 2.0, n=n, limiter=None) for n in _ns])
DISP_GRID_ORDER_UPWIND = float(-np.polyfit(np.log(_ns[:-1]),
                                           np.log(np.abs(_exu[:-1] - _exu[-1])), 1)[0])

rows = []
for N in [2, 4, 8, 16, 32, 64]:
    for kt in [0.5, 1.0, 2.0, 5.0]:
        tis = eq9(N, kt)
        dis = dispersion_exit(2.0 * N, kt)
        pe_eff = brentq(lambda Pe: dispersion_exit(Pe, kt, n=1200) - tis,
                        0.6 * N, 400.0 * N, xtol=1e-6, rtol=1e-10)
        rows.append(dict(N=N, ktbar=kt, tis=tis, disp_Pe2N=dis,
                         rel_diff=abs(dis - tis) / tis, Pe_eff=pe_eff,
                         Pe_ratio=pe_eff / (2.0 * N)))
DISP = pd.DataFrame(rows)
DISP_MAX_DIFF = float(DISP.rel_diff.max())
DISP_ROW_WORST = DISP.loc[DISP.rel_diff.idxmax()]
DISP_DIFF_N16_KT2 = float(DISP[(DISP.N == 16) & (DISP.ktbar == 2.0)].rel_diff.iloc[0])
DISP_DIFF_N64_KT05 = float(DISP[(DISP.N == 64) & (DISP.ktbar == 0.5)].rel_diff.iloc[0])
_sub = DISP[DISP.ktbar == 1.0]
DISP_DIFF_ORDER_IN_N = float(-np.polyfit(np.log(_sub.N), np.log(_sub.rel_diff), 1)[0])
PE_RATIO_N64 = float(DISP[(DISP.N == 64) & (DISP.ktbar == 0.5)].Pe_ratio.iloc[0])
PE_RATIO_N4_KT5 = float(DISP[(DISP.N == 4) & (DISP.ktbar == 5.0)].Pe_ratio.iloc[0])

# ---- what Fig. 14.7's rule costs AS PRINTED ------------------------------
# printed: N = HALF (D/uL)  ->  D/uL = N/HALF  ->  Pe = HALF/N
# derived: N = HALF (uL/D)  ->  Pe = uL/D = N/HALF = 2N
# (HALF and N_SYM come off the CSV in the printed-defects section above.)
_Nd, _ktd = 64, 2.0
_tis_d = eq9(_Nd, _ktd)
DEFECT_FIG147_PE_PRINTED = HALF / _Nd
DEFECT_FIG147_PE_DERIVED = _Nd / HALF
DEFECT_FIG147_PRINTED_RELDEV = abs(dispersion_exit(DEFECT_FIG147_PE_PRINTED, _ktd)
                                   - _tis_d) / _tis_d
DEFECT_FIG147_DERIVED_RELDEV = abs(dispersion_exit(DEFECT_FIG147_PE_DERIVED, _ktd)
                                   - _tis_d) / _tis_d

# ---- the third matching rule: variance, on the CLOSED vessel -------------
def sigma2_theta_closed(Pe):
    """Levenspiel's closed-vessel variance (Ch. 13), used here only to compare
    matching rules - the small-deviation form 2(D/uL) gives Pe = 2N exactly."""
    return 2.0 / Pe - 2.0 / Pe ** 2 * (1.0 - np.exp(-Pe))


PEVAR = pd.DataFrame([
    dict(N=N, Pe_var=brentq(lambda p: sigma2_theta_closed(p) - 1.0 / N, 1e-3, 1e5,
                            xtol=1e-12, rtol=8.9e-16))
    for N in [2, 4, 16, 64]])
PEVAR["two_N_minus_Pe"] = 2.0 * PEVAR.N - PEVAR.Pe_var
PE_VAR_OFFSET_N64 = float(PEVAR.two_N_minus_Pe.iloc[-1])
PE_EFF_OFFSET_N64 = float(2 * 64 - DISP[(DISP.N == 64) & (DISP.ktbar == 0.5)].Pe_eff.iloc[0])

display(DISP.round(6))
display(PEVAR.round(6))
display(Markdown(f"""
The dispersion solve converges at observed grid order
**{DISP_GRID_ORDER:.2f}** with the van Leer deferred correction and
**{DISP_GRID_ORDER_UPWIND:.2f}** with bare upwind &mdash; the usual pair, and the
only genuine grid axis on this page.

At $\\mathrm{{Pe}} = 2N$ the two models' first-order exit concentrations differ by
**{DISP_DIFF_N64_KT05 * 100:.3f} %** at $N = 64$, $k\\bar t = 0.5$, rising to
**{DISP_MAX_DIFF * 100:.1f} %** at $N = {int(DISP_ROW_WORST.N)}$,
$k\\bar t = {DISP_ROW_WORST.ktbar:g}$. At fixed $k\\bar t$ the gap falls as
$N^{{-{DISP_DIFF_ORDER_IN_N:.2f}}}$. Turned around: the Péclet number that exactly
reproduces $N$ tanks is {PE_RATIO_N64:.4f} of $2N$ at $N = 64$, $k\\bar t = 0.5$, but
only {PE_RATIO_N4_KT5:.4f} of it at $N = 4$, $k\\bar t = 5$.

So Levenspiel's "identical for all practical purposes" is a claim about the
*deviation from plug flow*, and it is a good one: at $N \\ge 16$ and
$k\\bar t \\le 2$ the two models are within {DISP_DIFF_N16_KT2 * 100:.2f} %, which no
tracer experiment resolves. It is **not** a claim you can carry to a strongly
converting reactor with a broad RTD, and the table above says by how much.

**What Fig. 14.7's rule costs as printed.** Read literally, $N = {HALF:g}
(\\mathbf{{D}}/uL)$ sends a reader at $N = {_Nd}$ to
$\\mathrm{{Pe}} = {DEFECT_FIG147_PE_PRINTED:.6g}$ instead of
{DEFECT_FIG147_PE_DERIVED:g} &mdash; essentially a single stirred tank. At
$k\\bar t = {_ktd:g}$ that is **{DEFECT_FIG147_PRINTED_RELDEV * 100:.0f} %** away from
the {_Nd}-tank conversion, against {DEFECT_FIG147_DERIVED_RELDEV * 100:.3f} % for
the de-inverted form. The de-inversion is not an opinion: the annotation's own
threshold $N > {N_SYM:g}$ maps through it to
$\\mathbf{{D}}/uL < {DEFECT_FIG147_THRESHOLD_DUL:g}$, which matches the printed
validity limit of Fig. 13.11 itself to {DEFECT_FIG147_THRESHOLD_DEV:.0e}; and
Ch. 13's own printed $\\sigma^2_\\theta = 2(\\mathbf{{D}}/uL)$ with Ch. 14's
$\\sigma^2_\\theta = 1/N$ gives the coefficient {HALF:g} to
{DEFECT_FIG147_FROM_CH13:.0e}, on the *other* ratio.

**Three matching rules, reconciled.** $\\mathrm{{Pe}} = 2N$ is what the upwind
truncation error gives and it is also what Ch. 13's small-deviation variance
gives, exactly. Matching the **closed-vessel** variance instead gives the
`Pe_var` column above, and matching the **conversion** gives the root-found
`Pe_eff` column. Both sit a little below $2N$ &mdash; at $N = 64$,
$2N - \\mathrm{{Pe}}$ is {PE_VAR_OFFSET_N64:.3f} on the variance rule and
{PE_EFF_OFFSET_N64:.3f} on the conversion rule at $k\\bar t = 0.5$ &mdash; that is, the
three rules differ by an $O(1)$ offset in $\\mathrm{{Pe}}$ and therefore by
$O(1/N)$ relatively. They are the same rule to leading order, which is the only
order the annotation claims.
"""))'''))

# ------------------------------------------------------------- validation ---
cells.append(md(r"""## Validation

### The $F$ curve, and the branch its brace prescribes

Eq. 8 is checked against a running integral of Eq. 3 — a different object
computed a different way — and at both ends of the brace the book draws under
it: the "for one tank use the first term" branch, and the general $N$ case.""" ))

cells.append(code(r'''# Eq. 8: the F curve, against a numerical integral of Eq. 3 - and both branches
_thF = np.linspace(1e-9, 6.0, 40001)
_hF = _thF[1] - _thF[0]
F_DEV = 0.0
for N in [1, 2, 3, 5, 10, 30]:
    e = E_theta(_thF, N)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * _hF)])
    F_DEV = max(F_DEV, float(np.max(np.abs(F_eq8(_thF, N) - cum))))
F_ENDPOINT_DEV = max(abs(F_eq8(1e6, N) - 1.0) for N in [1, 2, 5, 30])
F_ONE_TANK_DEV = float(np.max(np.abs(F_eq8(_thF, 1) - (1.0 - np.exp(-_thF)))))
EQ3_VAR_QUAD_DEV = VAR_DEV
FIG143_EMAX_N4 = float(FIG143.loc[FIG143.N == N4, "Emax_root"].iloc[0])

display(Markdown(f"""
Eq. 8 against the running integral of Eq. 3, for $N$ = 1, 2, 3, 5, 10 and 30:
**{F_DEV:.1e}** (trapezoid-limited). $F(\\infty) - 1$ is **{F_ENDPOINT_DEV:.1e}**.
The "for one tank use the first term" branch reduces to $1-e^{{-\\theta}}$ to
**{F_ONE_TANK_DEV:.1e}**, so the brace under Eq. 8 is exercised at its own first
term and not only in the general case.
"""))'''))

cells.append(md(r"""### The defect-injection table

Rebuilt for this physics; nothing in it travelled from another page. Each row
breaks one thing on purpose and reports what moves.

**What is measured and what is declared.** Every row carries a `metric` column
naming the entry of `agreement.json` whose *injected* value the `value` column
literally is; the row's `baseline` is then asserted, at the bottom of the
notebook, to equal that metric's clean value. A row cannot therefore claim to
cover something it did not recompute — the earlier draft of this page had one row
that reported a metric against zero, injected no defect at all, and was
nevertheless named as the mover of seven metrics. The wider `moves` column is the
row's *declared* scope: each name in it is checked to be a real key of
`agreement.json`, but the association itself is a judgement, and it is labelled as
one rather than as a measurement.

A row whose shift is smaller than the quantity it is supposed to protect is not
coverage, and is labelled.""" ))

cells.append(code(r'''BREAK = []


def brk(name, value, baseline, metric, moves=(), note=""):
    """One injected defect.

    `metric` is the agreement.json key whose injected value `value` IS, so that
    `baseline` can be asserted against the reported number; pass None only for a
    row whose value is a proxy, and say so in `note`.  `moves` is the declared
    scope, checked name-by-name against agreement.json.
    """
    names = tuple(dict.fromkeys(((metric,) if metric else ()) + tuple(moves)))
    BREAK.append(dict(defect=name, value=value, baseline=baseline,
                      shift_abs=abs(value - baseline), metric=metric,
                      moves=names, note=note))


# --- the printed-intermediate chains -----------------------------------
brk("E14.1: scale sigma, not sigma^2, with L",
    abs(E141_WRONG_POWER - E141_S2_P), E141_S2_DEV, "e14_1_sigma2_dev",
    ["e14_1_sigma2_sq_dev", "e14_1_wrong_power_gap"],
    "moves the answer by two and a half bottles")
brk("E14.2: spread ~ distance, not spread^2",
    abs(E142_WRONG_POWER - E142_L_P), E142_L_DEV, "e14_2_L_dev",
    ["e14_2_two_routes", "e14_2_wrong_power_gap"],
    "moves the source two hundred miles up the Ohio")
brk("E14.4: N from Delta(sigma), not Delta(sigma^2)",
    abs(DTB ** 2 / (np.sqrt(V_OUT) - np.sqrt(V_IN)) ** 2 - N4_P), E144_N_DEV,
    "e14_4_N_dev")
brk("E14.4: ignore the input signal (t_out^2/var_out)",
    abs(TB_OUT ** 2 / V_OUT - N4_P), E144_N_DEV, "e14_4_N_dev", (),
    "this is the point of Eq. 4")
brk("E14.4/Eq. 5: moments ADDED instead of subtracted",
    abs((TB_OUT + TB_IN) - DTB_P), E144_DTB_DEV, "e14_4_delta_tbar_dev",
    ["e14_4_delta_var_dev", "e14_4_N_dev", "e14_4_exponent_dev"])
brk("E14.4: prefactor with N! instead of (N-1)!",
    abs(N4 ** N4 / (np.exp(gammaln(N4 + 1)) * TB4 ** N4) - PRE_P) / PRE_P,
    E144_PRE_DEV, "e14_4_prefactor_reldev", ["e14_4_prefactor_two_routes"])
brk("E14.4: sampled maximum instead of root-found", E144_EMAX_SAMPLED, E144_EMAX,
    "e14_4_Emax", ["sweep_vs_rootfind_relerr"],
    "small, and still a grid max - that IS the finding")
brk("Fig E14.4b: tick judged against the UNROUNDED E(t)", TICK_P / E144_EMAX,
    DEFECT_TICK_RATIO, "defect_e144_tick_ratio", ["e14_4_Emax_printed"],
    "the tick belongs to the curve the book PRINTS, 3.2922e-6 and 0.0667")

# --- Eq. 3 itself, which everything downstream of it depends on ---------
def _E_theta_bad(th, N):
    """Eq. 3 with N! where the book writes (N-1)! - the classic transcription slip."""
    th = np.maximum(np.asarray(th, float), 1e-300)
    return np.exp(np.log(N) + (N - 1) * np.log(N * th) - N * th - gammaln(N + 1))


_ebad = _E_theta_bad(_th, 4)
brk("Eq. 3: (N-1)! written as N!", abs(float(_simp(_ebad)) - 1.0), AREA_DEV,
    "fig143_area_max_dev",
    ["fig143_mean_max_dev", "eq3_var_quadrature_dev", "fig143_Emax_max_reldev",
     "fig143_Emax_at_N4", "e14_4_Emax", "eq8_F_vs_integral_max_dev", "rtd_max_dev"])
brk("Eq. 3: theta_max read on the PER-TANK time theta_i = N theta",
    abs(float(N4 - 1) - (N4 - 1) / N4), THETA_MAX_DEV, "fig143_theta_max_max_dev",
    ["fig143_width_max_dev"])

# --- Fig. 14.3 ----------------------------------------------------------
brk("Fig 14.3: Stirling with sqrt(2 pi N), not sqrt(2 pi (N-1))",
    brentq(lambda N: abs(N / np.sqrt(2 * np.pi * N) - Emax_printed(N)) / Emax_printed(N)
           - P("Fig14.3", "stirling_error_bound"), 2.5, 400.0),
    STIRLING_NCRIT, "fig143_stirling_threshold_N",
    ["fig143_stirling_err_N5", "fig143_stirling_err_N6"])
_bad_infl = lambda th, N: N * (N - 2) / (N * th) ** 2 - 2 * (N - 1) / (N * th) + 1.0
_lo = brentq(lambda th: _bad_infl(th, N4), 1e-9, (N4 - 1) / N4, xtol=1e-15)
_hi = brentq(lambda th: _bad_infl(th, N4), (N4 - 1) / N4, 60.0, xtol=1e-14)
brk("Fig 14.3: inflection quadratic with N(N-2), not (N-1)(N-2)",
    abs((_hi - _lo) / ((N4 - 1) / N4) - 2 / np.sqrt(N4 - 1)), WIDTH_DEV,
    "fig143_width_max_dev", ["fig143_Einf_lo_N4", "fig143_Einf_hi_N4"])

# --- Eq. 8's truncation, which the book's own brace prescribes ----------
def _F_bad(th, N):
    """Eq. 8 summed to N terms instead of the printed N-1."""
    x = N * np.asarray(th, float)
    term = np.ones_like(x); total = term.copy()
    for j in range(1, N + 1):
        term = term * x / j
        total = total + term
    return 1.0 - np.exp(-x) * total


_thb = np.linspace(1e-9, 6.0, 40001)
_eb = E_theta(_thb, 5)
_cumb = np.concatenate([[0.0], np.cumsum(0.5 * (_eb[1:] + _eb[:-1]) * (_thb[1] - _thb[0]))])
brk("Eq. 8: series summed to N terms instead of N-1",
    float(np.max(np.abs(_F_bad(_thb, 5) - _cumb))), F_DEV,
    "eq8_F_vs_integral_max_dev", ["eq8_one_tank_branch_dev"])

# --- the printed small-deviation asymptote ------------------------------
brk("asymptote: (k tbar)^2/N instead of (k tbar)^2/(2N)",
    float(abs((eq9(320, 1.0) / np.exp(-1.0) - 1.0) / (1.0 / 320) - 1.0)),
    ASYM_DEV_N320, "asymptote_reldev_N320", ["asymptote_order_in_N"])

# --- the pymrm assembly -------------------------------------------------
_base_fo, _ = Cascade(N4, tbar=1.0, c_in=1.0).first_order(2.0)
_eq9_ref = eq9(N4, 2.0)
brk("pymrm: zero-gradient outflow bc (the PDE choice)",
    abs(Cascade(N4, tbar=1.0, c_in=1.0, outlet="zerograd").first_order(2.0)[0]
        - _eq9_ref), PYMRM_EQ9_DEV, "pymrm_eq9_max_dev",
    ["outflow_bc_relerr_N4", "outflow_bc_relerr_N2", "outflow_bc_decay_order_in_N"],
    "a MODELLING error, not an accuracy one - and the resolving row of the four")
brk("pymrm: geometrically graded cells (tanks unequal)",
    abs(Cascade(N4, tbar=1.0, c_in=1.0,
                faces=np.concatenate([[0.0], np.cumsum(1.6 ** np.arange(N4))
                                      / np.sum(1.6 ** np.arange(N4))])).first_order(2.0)[0]
        - _eq9_ref), PYMRM_EQ9_DEV, "pymrm_eq9_max_dev", (),
    "the second resolving row: the model still solves, it solves something else")
brk("pymrm: construct_div nu = 1 (cylindrical)",
    abs(Cascade(N4, tbar=1.0, c_in=1.0, nu=1).first_order(2.0)[0] - _eq9_ref),
    PYMRM_EQ9_DEV, "pymrm_eq9_max_dev", (),
    "COLLAPSE, not a shift: the exit goes to exactly 0, so this shows the metric "
    "reacting rather than resolving")
brk("pymrm: velocity reversed (downwind)",
    abs(Cascade(N4, tbar=1.0, c_in=1.0, v_sign=-1.0).first_order(2.0)[0] - _eq9_ref),
    PYMRM_EQ9_DEV, "pymrm_eq9_max_dev", (), "COLLAPSE, as above")
brk("pymrm: implicit Euler instead of Crank-Nicolson (time order)",
    ORDER_IE, ORDER_CN, "rtd_order_cn", ["rtd_order_ie", "rtd_max_dev"])
_zg_ts, _zg_es = Cascade(N4, tbar=TB4, outlet="zerograd").pulse_response(TMAX_PROD, NT_PROD)
_hz = _zg_ts[1] - _zg_ts[0]
brk("pymrm: RTD area with the zero-gradient outflow bc",
    abs(float(_hz / 3.0 * (_zg_es[0] + _zg_es[-1] + 4 * _zg_es[1:-1:2].sum()
                           + 2 * _zg_es[2:-1:2].sum())) - 1.0), RTD_AREA_DEV,
    "rtd_area_dev", ["rtd_area_residual"])
_nu_ts, _nu_es = Cascade(N4, tbar=TB4, nu=1).pulse_response(TMAX_PROD, NT_PROD)
_hn = _nu_ts[1] - _nu_ts[0]
brk("pymrm: RTD area with construct_div nu = 1 (cylindrical)",
    abs(float(_hn / 3.0 * (_nu_es[0] + _nu_es[-1] + 4 * _nu_es[1:-1:2].sum()
                           + 2 * _nu_es[2:-1:2].sum())) - 1.0), RTD_AREA_DEV,
    "rtd_area_dev", ["rtd_area_residual"])
_so = SO[(SO.N == N4) & (SO.kC0tbar == 2.0)].iloc[0]
brk("pymrm: Eq. 10 nested N-1 times instead of N",
    abs(_so.pymrm - eq10(N4 - 1, 2.0, 1.0 / N4, 1.0)), PYMRM_EQ10_DEV,
    "pymrm_eq10_max_dev",
    ["macro_micro_max_gap"])

# --- Eq. 11 and the quadrature -----------------------------------------
brk("Eq. 11: denominator tbar_N (as printed) instead of tbar^N",
    abs(DEFECT_EQ11_NORM_SUB - 1.0), EQ11_NORM_DEV, "eq11_normalisation_dev",
    ["eq11_vs_eq9_max_dev", "defect_eq11_subscript_norm"])
brk("Eq. 11: Simpson with 1/16 of the panels",
    abs(eq11(4, 1.0, lambda t: np.exp(-2.0 * t), npanel=2500) - eq9(4, 2.0)),
    EQ11_EQ9_DEV, "eq11_vs_eq9_max_dev", ["eq11_quad_order", "macro_micro_max_gap"])

# --- additivity: rows that INJECT a defect, not ones that report a metric ---
brk("Eq. 4: moment integrals truncated at 340 s instead of 800 s",
    EQ4_MEAN_DEV_SHORT, EQ4_MEAN_DEV, "eq4_mean_add_reldev",
    ["eq4_var_add_reldev"],
    "what those two metrics actually measure is this window")
brk("Eq. 4: moment-match the composite by TANK COUNT (N = 5), not by the moments",
    max_brent(lambda t: abs(float(convolve_E(np.array([t]), MA, MB, tA, tB)[0])
                            - float(E_time(t, MA + MB, MEAN_UNEQ))),
              1e-9, EQ4_WINDOW) / CONV_UNEQ_PEAK,
    EQ4_UNEQUAL_GAP, "eq4_unequal_shape_gap")

# --- the headline -------------------------------------------------------
_k10 = AMBIG_KT_HEADLINE / DTB
_c_eq10 = chain(T_EQ, DTB).first_order(_k10)[0]
brk("ambiguity: alternative chain sized by COUNT (5 equal tanks), so the "
    "variance no longer matches",
    abs(chain(np.full(5, DTB / 5), DTB).first_order(_k10)[0] - _c_eq10) / _c_eq10,
    AMBIG_FIRST_MAX, "ambiguity_first_order_max_gap",
    ["ambiguity_second_order_max_gap", "ambiguity_var_dev",
     "ambiguity_conc_gap_kt50", "ambiguity_conversion_max_gap"])
brk("ambiguity: RTD of the alternative chain taken as 5 EQUAL tanks",
    max_brent(lambda t: abs(float(E_time(t, 5, DTB)) - float(E_time(t, 4, DTB))),
              1e-9, 260.0) / E144_EMAX,
    AMBIG_RTD_GAP, "ambiguity_rtd_peak_gap")
brk("family bound: input variance not subtracted (S = var_out, not Delta var)",
    abs(np.exp(-_k10 * (DTB - np.sqrt(V_OUT))) / (1.0 + _k10 * np.sqrt(V_OUT))
        - _c_eq10) / _c_eq10,
    AMBIG_FAMILY_SUP_GAP, "ambiguity_family_sup_gap",
    ["ambiguity_family_sup_gap_kt2", "ambiguity_family_enum_gap"])

# --- the dispersion comparison and Fig. 14.7's printed rule -------------
brk("dispersion: no van Leer deferred correction (grid order)",
    DISP_GRID_ORDER_UPWIND, DISP_GRID_ORDER, "disp_grid_order",
    ["disp_grid_order_upwind"])
brk("dispersion: Pe = N instead of 2N",
    abs(dispersion_exit(16.0, 2.0) - eq9(16, 2.0)) / eq9(16, 2.0), DISP_DIFF_N16_KT2,
    "disp_reldiff_N16_kt2", ["disp_reldiff_max", "disp_reldiff_order_in_N",
                             "disp_reldiff_window_max", "disp_Pe_ratio_N64",
                             "disp_Pe_ratio_N4_kt5"])
brk("Fig 14.7: its rule AS PRINTED, N = 1/2 (D/uL), at N = 16",
    abs(dispersion_exit(HALF / 16.0, 2.0) - eq9(16, 2.0)) / eq9(16, 2.0),
    DISP_DIFF_N16_KT2, "disp_reldiff_N16_kt2",
    ["defect_fig147_printed_rule_reldev", "defect_fig147_derived_rule_reldev"],
    "printed defect (a), injected")
brk("Fig 14.7: threshold N > 50 mapped with the rule AS PRINTED",
    abs(DEFECT_FIG147_THRESHOLD_DUL_PRINTED - DUL_LIMIT),
    DEFECT_FIG147_THRESHOLD_DEV, "defect_fig147_threshold_dev")
brk("Fig 14.7: coefficient taken off the D/uL side of Ch. 13's variance",
    abs(SIG_COEF - HALF), DEFECT_FIG147_FROM_CH13, "defect_fig147_from_ch13")
brk("dispersion: variance matched on the open-vessel 2/Pe, not the closed vessel",
    abs(2.0 * 64 - brentq(lambda p: 2.0 / p - 1.0 / 64, 1e-3, 1e5, xtol=1e-12)),
    PE_VAR_OFFSET_N64, "disp_Pe_var_offset_N64", (),
    "the open-vessel rule gives Pe = 2N exactly, which is the point")

BT = pd.DataFrame(BREAK)
pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 58)
display(BT[["defect", "value", "baseline", "shift_abs", "metric", "note"]].round(8))
BREAK_MIN_SHIFT = float(BT[~BT.defect.str.contains("sampled maximum")].shift_abs.min())
display(Markdown(
    f"{len(BT)} rows, every one of which recomputes the metric named in its "
    f"`metric` column. Every row except the sampled-maximum one moves its target "
    f"by at least **{BREAK_MIN_SHIFT:.2e}**. The sampled-maximum row moves "
    f"{abs(E144_EMAX_SAMPLED - E144_EMAX):.1e} s⁻¹, which is *small* &mdash; and that is "
    f"the finding, not a pass: a 1001-point sweep of this curve gets within "
    f"{abs(E144_EMAX_SAMPLED - E144_EMAX) / E144_EMAX * 1e2:.4f} % of the true maximum, "
    f"so a swept extremum here would look converged while still being a grid value. "
    f"The page root-finds it."))'''))

cells.append(code(r'''# The outflow-boundary-condition penalty, across N - the one modelling error a
# reader of this page is most likely to reproduce.
rows = []
for N in [2, 4, 10, 40, 160]:
    good, _ = Cascade(N, tbar=1.0, c_in=1.0).first_order(2.0)
    bad, _ = Cascade(N, tbar=1.0, c_in=1.0, outlet="zerograd").first_order(2.0)
    rows.append(dict(N=N, model=good, zero_gradient=bad, rel_err=abs(bad - good) / good))
OUTFLOW = pd.DataFrame(rows)
OUTFLOW_N2 = float(OUTFLOW.loc[OUTFLOW.N == 2, "rel_err"].iloc[0])
OUTFLOW_N4 = float(OUTFLOW.loc[OUTFLOW.N == 4, "rel_err"].iloc[0])
OUTFLOW_ORDER = float(-np.polyfit(np.log(OUTFLOW.N), np.log(OUTFLOW.rel_err), 1)[0])

_asm = BT[BT.moves.apply(lambda ms: "pymrm_eq9_max_dev" in ms)]
_res = _asm[~_asm.note.str.startswith("COLLAPSE")]
ASSEMBLY_MIN_SHIFT = float(_res.shift_abs.min())
ASSEMBLY_MAX_SHIFT = float(_res.shift_abs.max())
ASSEMBLY_N_COLLAPSE = int(len(_asm) - len(_res))

_win = DISP[(DISP.N >= 16) & (DISP.ktbar <= 2.0)]
DISP_WINDOW_MAX = float(_win.rel_diff.max())
PE_RATIO_N64_WORST = float((1.0 - DISP[DISP.N == 64].Pe_ratio).max())
EMAX_SWEEP_RELERR = abs(E144_EMAX_SAMPLED - E144_EMAX) / E144_EMAX

display(OUTFLOW.round(8))
display(Markdown(f"""
**The outflow condition, quantified.** Using pymrm's zero-gradient outflow instead
of the stirred-tank one costs **{OUTFLOW_N2 * 100:.1f} %** at $N=2$ and
**{OUTFLOW_N4 * 100:.1f} %** at $N=4$, decaying as $N^{{-{OUTFLOW_ORDER:.2f}}}$ over
$N = 2$ to 160 &mdash; close enough to first order to pass for an ordinary
discretisation error, which is what makes it easy to dismiss. It is not one:
$N$ is the model, so this error does not go away by refining anything.

**What the assembly rows are worth.** Four rows break the pymrm assembly, and
they are not equal in value. {ASSEMBLY_N_COLLAPSE} of them &mdash; `nu=1` and a
reversed velocity &mdash; drive the exit concentration to exactly zero: the model
*collapses* rather than shifts, so they show the metric reacting, not resolving.
The two that resolve are the zero-gradient outflow and the graded cells, and they
carry the row on their own: they move `pymrm_eq9_max_dev` by
**{ASSEMBLY_MIN_SHIFT:.1e}** and **{ASSEMBLY_MAX_SHIFT:.1e}**, against a reported
{PYMRM_EQ9_DEV:.1e}. Those are the two numbers to quote.

**The sweep-versus-root-find margin.** A 1001-point sweep of E14.4's own curve
lands within **{EMAX_SWEEP_RELERR * 100:.4f} %** of the true maximum &mdash; close
enough to pass for converged, and still a grid value.
"""))'''))

cells.append(md(r"""### What the checks can and cannot do

**Structural, and labelled as such.**

- `pymrm_eq9_max_dev`. The $N$-cell upwind cascade is *algebraically identical*
  to Eq. 9; agreement at round-off is an identity, not evidence about the
  physics. What it can still catch is every way of assembling it wrongly, and
  four break rows do exactly that — the zero-gradient outflow, `nu=1`, a reversed
  velocity and graded cells — with the range of shifts printed just above.
- `eq4_equal_convolution_dev` **and `rtd_convolution_dev`**. These are the same
  identity at two stage counts: Gamma $\star$ Gamma at equal scale *is* Gamma.
  Neither can detect anything physical, and — this is the part an earlier draft
  got wrong — neither says anything about the quadrature either. For equal tank
  sizes the integrand $E_M(s)E_N(t-s)$ is $s^{M-1}(t-s)^{N-1}e^{-t/\bar t_i}$, a
  polynomial of degree $M+N-2 \le 3$ times a factor constant in $s$, and Simpson
  is exact on cubics: the same $10^{-18}$ comes back on 20 panels as on 2000, and
  it survives the $(N-1)!\to N!$ slip because the Beta integral absorbs it. Their
  above-floor companion is `eq4_unequal_shape_gap`, which is the same machinery
  applied where the identity does *not* hold and where the quadrature is genuinely
  refined.
- `eq3_variance_identity`. $\bar t^{\,2}/N$ and $N\bar t_i^{\,2}$ are the same
  number by the definition of $\bar t_i$; the check confirms the transcription of
  Eq. 3's two right-hand columns and nothing else. Companion:
  `eq3_var_quadrature_dev`, which integrates the density.
- `ambiguity_mean_dev` and `ambiguity_var_dev`. Zero by construction and asserted
  in code: the five-stage chain is *chosen* to match the two moments. They are
  reported so that a reader can see the match is exact rather than approximate,
  and the companion is `ambiguity_first_order_max_gap`.

**Second, independent computations.** Four headlines are computed twice with no
shared assembly. The pairings below are the ones the code actually computes, and
each names the metric that reports it:

| quantity | route 1 | route 2 | metric |
|---|---|---|---|
| E14.2's $L$ | `brentq` on the printed ratio | closed form $119/(1-(10.5/14)^2)$ | `e14_2_two_routes` |
| E14.4's prefactor | $N^N/((N-1)!\bar t^{\,N})$ | numerical normalisation of $t^3e^{-4t/60}$ | `e14_4_prefactor_two_routes` |
| first-order conversion | Eq. 11's quadrature $\int e^{-kt}E\,dt$ | Eq. 9's algebra | `eq11_vs_eq9_max_dev` |
| the $N$-tank RTD | pymrm Crank–Nicolson march | Eq. 3's closed form | `rtd_max_dev` |

The third is the strongest: a quadrature over Eq. 3 has no way of knowing that the
answer must be $(1+k\bar t/N)^{-N}$, and it would not agree if Eq. 3, Eq. 11's
denominator, or the quadrature were wrong. The Simpson convolution
$E_M\star E_N$ is **not** on this list: for equal tanks it is an identity that
Simpson integrates exactly (see above), so it is not a second route to anything.
`pymrm_eq9_max_dev` is not on it either — the cascade *is* Eq. 9, so that pair
shares its algebra rather than avoiding it.

**Refinement axes.** There is deliberately **no grid axis** for the cascade — the
cell count is Levenspiel's $N$, and refining it changes the model rather than the
error. The axes that do carry numerical error are all measured: the time step
(order 2 with Crank–Nicolson, 1 with implicit Euler), the Simpson panel count for
Eq. 11, and — in the dispersion comparison only, where there really is a
discretisation — the grid, at order 2 with the deferred correction and 1 without.

**Every printed branch is exercised somewhere it is live.** Eq. 9 and Eq. 10 are
evaluated at $N$ = 1, 2, 3, 4, 6, 10, 20 and 40, so Eq. 10's radical nest runs at
depths 0, 1, 2, 3, 5, 9, 19 and 39 — eight depths spanning that range, not every
integer in it. Eq. 8's series is summed at $N$ = 1 (its "for one tank use the
first term" branch), 2, 3, 5, 10 and 30. Eq. 11 is evaluated on both of its
branches: first order, where Levenspiel says it must equal Eq. 9, and second
order, where it must not. Eq. 4 is exercised both where it is exact (equal tanks)
and where its parenthesis matters (unequal blocks), and **Eq. 5** — the one-shot
tracer relation $\Delta\sigma^2 = (\Delta\bar t)^2/N$, book p. 324 — is the most
exercised relation on the page: run backwards it *is* E14.4's fit, and section 8
is a statement about what it leaves undetermined.

The printed relations **not** exercised anywhere are Fig. 14.8's graphical
construction for arbitrary kinetics (p. 328), which is a drawing rather than a
formula, and the closed-recirculation superpositions **Eqs. 6a–c and 7a–c**
(p. 325). Both are named in *The published model*, and neither is claimed as
tested. (Eqs. 6–7 are recirculation results; Eq. 5 is not, and an earlier draft of
this page said otherwise.)""" ))

# --------------------------------------------------------- what pymrm adds --
cells.append(md(r"""## What pymrm adds""" ))

cells.append(code(r'''display(Markdown(f"""
**To the physics, nothing.** Chapter 14 is right. Its three worked examples
reproduce with zero deviation, its Fig. 14.3 annotations reproduce to
{max(THETA_MAX_DEV, WIDTH_DEV):.0e}, its Eq. 8, Eq. 9 and Eq. 10 are correctly
transcribed and internally consistent, and its assertion that the model coincides
with the dispersion model for small deviations from plug flow is confirmed and
quantified. A 1999 textbook chapter is not the place to look for physics errors,
and none was found.

Four things are new here, and none of them is an improvement to the model:

**1. The identity with the upwind scheme, and the correction of the chapter's own
conversion rule.** $N$ equal tanks in series *is* `construct_convflux_upwind` +
`construct_div` on $N$ cells, exactly ({PYMRM_EQ9_DEV:.0e} against Eq. 9 over
{len(FO)} cases), and the equivalent Peclet number $2N$ follows from the scheme's
numerical diffusion $u\\Delta z/2$. Chapter 14 does print a rule of its own, on
Fig. 14.7 &mdash; $N = {HALF:g}(\\mathbf{{D}}/uL)$ &mdash; and it is **inverted**;
the derivation here is the same rule the right way up, and the book's own Ch. 13
settles it in two independent ways. Read as printed the rule is
{DEFECT_FIG147_PRINTED_RELDEV * 100:.0f} % wrong in conversion at $N = 64$,
$k\\bar t = 2$, against {DEFECT_FIG147_DERIVED_RELDEV * 100:.3f} % de-inverted. With
the rule fixed, the page then measures where it stops holding &mdash;
{DISP_WINDOW_MAX * 100:.2f} % at worst for $N\\ge 16,\\ k\\bar t\\le 2$, but
{DISP_MAX_DIFF * 100:.0f} % at $N = {int(DISP_ROW_WORST.N)}$,
$k\\bar t = {DISP_ROW_WORST.ktbar:g}$.

**2. The measurement of what a two-moment fit is worth.** The five-stage chain
9, 9, 9, 9, 24 s matches E14.4's $\\bar t$ and $\\sigma^2$ to
{max(AMBIG_MEAN_DEV, AMBIG_VAR_DEV):.0e} and differs in **outlet concentration**
by **{AMBIG_FIRST_MAX * 100:.1f} %** at $k\\bar t = {AMBIG_KT_HEADLINE:g}$ (in
*conversion*, {AMBIG_CONV_HL * 100:.2f} % &mdash; the two are not the same statement
and the page prints both). That chain is a lower bound: over every chain matching
the same two moments the closed-form supremum is
**{AMBIG_FAMILY_SUP_GAP * 100:.1f} %**, i.e. E14.4's two numbers pin the outlet
concentration only to within a factor of
{1.0 / (1.0 - AMBIG_FAMILY_SUP_GAP):.2f} there. Neither number is in the chapter,
and together they are the practical answer to "how much can I trust an $N$ fitted
this way".

**3. Five printed defects**, each proved from the book's own arithmetic: the
**inverted dispersion rule** on Fig. 14.7 (the consequential one), the factor of
{DEFECT_FACTOR:.0f} in $V_{{N\\text{{tanks}}}}/V_p$, the transposed digit on
Fig. E14.4b's axis ({TICK_P:g} against a printed-curve maximum of
{E144_EMAX_PRINTED:.5f}), "272 moles", and a variance printed in seconds. Plus one
typesetting ambiguity in Eq. 11's denominator, settled by the book's own E14.4.

**4. A caution about pymrm itself.** There is no pure-outflow boundary condition
in the library. The nearest, zero-gradient, reconstructs the exit face to second
order &mdash; right for a discretised PDE, wrong for a model whose cell count is a
physical parameter, and wrong by **{OUTFLOW_N4 * 100:.1f} %** at $N=4$ while decaying
like $N^{{-{OUTFLOW_ORDER:.2f}}}$, which is close enough to what ordinary
discretisation error looks like to be dismissed as one. The workaround used here (suppress the boundary flux, add the
outflow as a `construct_coefficient_matrix` sink on the last cell) is exact and
cheap, but it had to be worked out rather than looked up.
"""))'''))

# ------------------------------------------------------------------ reuse ---
cells.append(md(r"""## Reuse

**When to reach for this page.** Any staged process where a chain of well-mixed
compartments is the model rather than a discretisation: multistage extractors and
absorbers, a bank of crystallisers, an anaerobic digester train, compartment
models in pharmacokinetics, or the "$N$ tanks" that a measured tracer curve is
routinely collapsed to. `Cascade` takes any `faces`, so unequal stages cost
nothing extra, and `order_n` takes any pointwise rate.

**What to carry over.**

- **Do not use a zero-gradient outflow condition when the cells are tanks.** It
  reconstructs the exit face to second order and quietly changes the last stage's
  balance; the table in *Validation* gives the penalty against $N$. This is the
  mirror image of [`A2.6`](../A2.6-gunn-dispersion-correlations/)'s outlet trap —
  there the *cell-centre* read was wrong and the face read right; here the face
  reconstruction is wrong and the cell value right — and the rule that covers
  both is: ask what the model says the exit stream is, then read that.
- **`NumJac((N, 1))`, never `NumJac((N,))`.** With a bare 1-D shape the default
  stencil declares every tank coupled to every other and builds a dense
  Jacobian. It is in `AGENTS.md` and it is live in this page's `order_n`.
- **Root-find extrema and thresholds.** The margin printed in *Validation* is the
  reason: a 1001-point sweep of E14.4's own curve gets close enough to the true
  maximum to pass for converged while still being a grid value. Fig. 14.3's
  inflection points and the 2 % Stirling threshold are root-found for the same
  reason.
- **If you fit $N$ from two moments, say what it is worth.** Report the
  conversion the fit implies *and* the spread over chains that match the same two
  moments. Section 8 shows the cheap version — one alternative chain — and
  section 8b the complete one: the extremal chain has at most **two** distinct
  tank sizes (a two-line Lagrange argument), so the envelope is an enumeration
  rather than a search, and its supremum is the closed form
  $\ln(1+k\sqrt{S}) + k(T-\sqrt{S})$ with $T = \sum t_i$, $S = \sum t_i^2$. That
  bound transfers to any staged model fitted on a mean and a variance.
- **Say whether a gap is in concentration or in conversion.** They differ by more
  than an order of magnitude here, in opposite directions along the sweep, and
  quoting the larger one under the other one's name is the easiest overstatement
  on this page to make.

**What this page cannot conclude.**

- **Nothing here is validated against experiment.** Chapter 14 contains no
  measurement of any kind. Every comparison is against the book's own arithmetic
  or against an analytical limit; the page is tier 6, and calling it validation
  would be wrong.
- **Independence is assumed, not tested.** Levenspiel's own footnote says laminar
  flow often fails it. Nothing on this page can detect a vessel whose fluid
  carries memory from stage to stage, and the additivity check would pass anyway.
- **The equal-tank assumption is untestable from an RTD's first two moments** —
  that is the headline, and it cuts both ways: this page cannot say that
  Levenspiel's vessel *is* four equal tanks either, only that four equal tanks
  and the five-stage chain are indistinguishable by his method.
- **Fig. 14.3's $E_{\theta,\text{inf}}$ annotation is not settled.** The two
  inflection points have different heights and the figure draws one level; both
  are reported and neither is scored.
- **The dispersion comparison depends on a matching rule.** $\mathrm{Pe} = 2N$ is
  what the scheme's numerical diffusion gives and also what Levenspiel's own
  small-deviation $\sigma^2_\theta = 2(\mathbf{D}/uL)$ gives, exactly. Matching on
  the *closed-vessel* variance, or on the conversion, gives Péclet numbers a
  little below $2N$ — an $O(1)$ offset, printed as `disp_Pe_var_offset_N64` and
  as the root-found `Pe_eff` column, so the three rules agree to $O(1/N)$
  relatively and not exactly. Nothing here establishes that either model
  describes a real vessel.
- **Four of the five printed defects are typographic; the fifth is not.**
  Fig. 14.7's inverted rule is a modelling statement, and a reader who applies it
  gets a different reactor. All five are argued from the book's own numbers, and
  the page adjusts nothing: it re-derives every affected result from the correct
  form and reports the printed one beside it.""" ))

# ---------------------------------------------------- metrics and coverage --
cells.append(code(r'''metrics = {
    # --- printed intermediates, E14.1 and E14.2 -----------------------------
    "e14_1_sigma2_sq_dev":              float(E141_S2SQ_DEV),
    "e14_1_sigma2_dev":                 float(E141_S2_DEV),
    "e14_1_wrong_power_gap":            float(abs(E141_WRONG_POWER - E141_S2) / E141_S2),
    "e14_2_L_dev":                      float(E142_L_DEV),
    "e14_2_two_routes":                 float(E142_TWO_ROUTES),
    "e14_2_wrong_power_gap":            float(abs(E142_WRONG_POWER - E142_L) / E142_L),
    # --- E14.4: the fit and its consequences -------------------------------
    "e14_4_delta_tbar_dev":             float(E144_DTB_DEV),
    "e14_4_delta_var_dev":              float(E144_DV_DEV),
    "e14_4_N_dev":                      float(E144_N_DEV),
    "e14_4_prefactor_reldev":           float(E144_PRE_DEV),
    "e14_4_prefactor_two_routes":       float(E144_PRE_TWO_ROUTES),
    "e14_4_exponent_dev":               float(E144_EXP_DEV),
    "e14_4_Emax":                       float(E144_EMAX),
    "e14_4_Emax_printed":               float(E144_EMAX_PRINTED),
    # --- Fig. 14.3 ---------------------------------------------------------
    "fig143_theta_max_max_dev":         float(THETA_MAX_DEV),
    "fig143_width_max_dev":             float(WIDTH_DEV),
    "fig143_Emax_max_reldev":           float(EMAX_DEV),
    "fig143_area_max_dev":              float(AREA_DEV),
    "fig143_mean_max_dev":              float(MEAN_DEV),
    "fig143_stirling_threshold_N":      float(STIRLING_NCRIT),
    "fig143_stirling_err_N5":           float(STIRLING_ERR_N5),
    "fig143_stirling_err_N6":           float(STIRLING_ERR_N6),
    "fig143_Einf_lo_N4":                float(EINF_LO_N4),
    "fig143_Einf_hi_N4":                float(EINF_HI_N4),
    "fig143_Emax_at_N4":                float(FIG143_EMAX_N4),
    # --- Eq. 3 and Eq. 8 ---------------------------------------------------
    "eq3_variance_identity":            float(EQ3_VAR_IDENTITY),
    "eq3_var_quadrature_dev":           float(EQ3_VAR_QUAD_DEV),
    "eq8_F_vs_integral_max_dev":        float(F_DEV),
    "eq8_F_endpoint_dev":               float(F_ENDPOINT_DEV),
    "eq8_one_tank_branch_dev":          float(F_ONE_TANK_DEV),
    # --- Eq. 4, additivity -------------------------------------------------
    "eq4_equal_convolution_dev":        float(EQ4_EQUAL),
    "eq4_mean_add_reldev":              float(EQ4_MEAN_DEV),
    "eq4_var_add_reldev":               float(EQ4_VAR_DEV),
    "eq4_unequal_shape_gap":            float(EQ4_UNEQUAL_GAP),
    # --- the pymrm cascade -------------------------------------------------
    "pymrm_eq9_max_dev":                float(PYMRM_EQ9_DEV),
    "pymrm_eq10_max_dev":               float(PYMRM_EQ10_DEV),
    "rtd_max_dev":                      float(RTD_MAX_DEV),
    "rtd_area_dev":                     float(RTD_AREA_DEV),
    "rtd_area_residual":                float(RTD_AREA_RESIDUAL),
    "rtd_order_cn":                     float(ORDER_CN),
    "rtd_order_ie":                     float(ORDER_IE),
    "rtd_convolution_dev":              float(RTD_CONV_DEV),
    # --- Eq. 11, the macrofluid branch -------------------------------------
    "eq11_vs_eq9_max_dev":              float(EQ11_EQ9_DEV),
    "eq11_normalisation_dev":           float(EQ11_NORM_DEV),
    "eq11_quad_order":                  float(EQ11_QUAD_ORDER),
    "macro_micro_max_gap":              float(MACRO_MICRO_MAX),
    # --- the printed asymptote and the printed defect -----------------------
    "asymptote_order_in_N":             float(ASYM_ORDER),
    "asymptote_reldev_N320":            float(ASYM_DEV_N320),
    "defect_volume_ratio_factor":       float(DEFECT_FACTOR),
    "defect_volume_ratio_right_form":   float(DEFECT_RIGHT_RATIO),
    "defect_eq11_subscript_norm":       float(DEFECT_EQ11_NORM_SUB),
    "defect_e144_tick_ratio":           float(DEFECT_TICK_RATIO),
    "defect_fig147_printed_rule_reldev": float(DEFECT_FIG147_PRINTED_RELDEV),
    "defect_fig147_derived_rule_reldev": float(DEFECT_FIG147_DERIVED_RELDEV),
    "defect_fig147_threshold_dev":      float(DEFECT_FIG147_THRESHOLD_DEV),
    "defect_fig147_from_ch13":          float(DEFECT_FIG147_FROM_CH13),
    # --- the headline ------------------------------------------------------
    "ambiguity_mean_dev":               float(AMBIG_MEAN_DEV),
    "ambiguity_var_dev":                float(AMBIG_VAR_DEV),
    "ambiguity_first_order_max_gap":    float(AMBIG_FIRST_MAX),
    "ambiguity_conc_gap_kt50":          float(AMB.conc_gap.iloc[-1]),
    "ambiguity_conversion_max_gap":     float(AMBIG_CONV_MAX),
    "ambiguity_second_order_max_gap":   float(AMBIG_SECOND_MAX),
    "ambiguity_rtd_peak_gap":           float(AMBIG_RTD_GAP),
    "ambiguity_pymrm_vs_product_dev":   float(AMBIG_PYMRM_DEV),
    "ambiguity_family_sup_gap":         float(AMBIG_FAMILY_SUP_GAP),
    "ambiguity_family_sup_gap_kt2":     float(AMBIG_FAMILY_SUP_GAP_KT2),
    "ambiguity_family_enum_gap":        float(AMBIG_FAMILY_BEST_N),
    # --- the dispersion comparison -----------------------------------------
    "disp_grid_order":                  float(DISP_GRID_ORDER),
    "disp_grid_order_upwind":           float(DISP_GRID_ORDER_UPWIND),
    "disp_reldiff_N16_kt2":             float(DISP_DIFF_N16_KT2),
    "disp_reldiff_max":                 float(DISP_MAX_DIFF),
    "disp_reldiff_order_in_N":          float(DISP_DIFF_ORDER_IN_N),
    "disp_Pe_ratio_N64":                float(PE_RATIO_N64),
    "disp_Pe_ratio_N4_kt5":             float(PE_RATIO_N4_KT5),
    "disp_reldiff_window_max":          float(DISP_WINDOW_MAX),
    "disp_Pe_var_offset_N64":           float(PE_VAR_OFFSET_N64),
    "outflow_bc_relerr_N4":             float(OUTFLOW_N4),
    "outflow_bc_relerr_N2":             float(OUTFLOW_N2),
    "outflow_bc_decay_order_in_N":      float(OUTFLOW_ORDER),
    "sweep_vs_rootfind_relerr":         float(EMAX_SWEEP_RELERR),
    # --- break-table floor -------------------------------------------------
    "break_min_shift":                  float(BREAK_MIN_SHIFT),
}


def _rows(tag):
    """Break-table rows that declare this metric, with the shift each MEASURED.

    A row that recomputed the metric itself - so that its `baseline` is asserted
    below to be the reported value - is flagged [recomputed]; the rest declare it
    as scope.
    """
    hits = BT[BT.moves.apply(lambda ms: tag in ms)]
    if len(hits) == 0:
        raise KeyError(f"no break row names {tag!r}")
    return "break table: " + "; ".join(
        f"{r.defect} ({r.shift_abs:.2e})" + (" [recomputed]" if r.metric == tag else "")
        for r in hits.itertuples())


# Only the metrics that are NOT simply "moved by" something are written by hand.
# Everything else has its coverage entry BUILT from the break table's own
# `moves` lists, so a metric no row touches raises rather than getting a
# sentence, and a row that stops moving something shortens a list rather than
# leaving a stale claim behind.
LABELS = {
    "fig143_Einf_lo_N4":              "REPORTED, NOT SCORED - the figure draws one level for "
                                      "two inflection points of different height; both are "
                                      "printed and the ambiguity is stated rather than "
                                      "resolved.",
    "fig143_Einf_hi_N4":              "REPORTED, NOT SCORED - see fig143_Einf_lo_N4.",
    "eq3_variance_identity":          "STRUCTURAL - tbar^2/N and N tbar_i^2 are the same "
                                      "number by the definition of tbar_i; it checks the "
                                      "transcription of Eq. 3's two variance columns and "
                                      "nothing else. Above-floor companion: "
                                      "eq3_var_quadrature_dev",
    "eq8_F_endpoint_dev":             "STRUCTURAL - exp(-N theta) underflows to exactly 0 at "
                                      "theta = 1e6, so F(inf) = 1 in floating point whatever "
                                      "the series is. It checks that the series is finite and "
                                      "nothing else. Companion: eq8_F_vs_integral_max_dev",
    "eq4_equal_convolution_dev":      "STRUCTURAL - Gamma * Gamma at equal scale IS Gamma, and "
                                      "Simpson is EXACT on the resulting cubic integrand, so "
                                      "this bounds neither physics nor quadrature: the same "
                                      "number comes back on 20 panels. Above-floor companion: "
                                      "eq4_unequal_shape_gap",
    "rtd_convolution_dev":            "STRUCTURAL - the same identity as "
                                      "eq4_equal_convolution_dev at M = N = 2, and exact for "
                                      "the same reason. It is NOT an independent route to the "
                                      "RTD and is not cited as one. Above-floor companion: "
                                      "rtd_max_dev, which is the real second route",
    "pymrm_eq9_max_dev":              "STRUCTURAL - against the PHYSICS only: the N-cell upwind "
                                      "cascade is algebraically Eq. 9 - but NOT against the "
                                      "assembly. Two of the four rows below resolve it and two "
                                      "collapse the model instead; the note column says which.",
    "defect_volume_ratio_factor":     "DEFECT MEASUREMENT - it IS the ratio of the "
                                      "printed middle form's slope to the true one, and it is "
                                      "EXACTLY 2 (the true slope is the limit x/ln(1+x) - 1 "
                                      "over x -> 1/2, extrapolated, not fitted). Its companion "
                                      "defect_volume_ratio_right_form is exactly 1, which is "
                                      "what identifies WHICH of the two printed forms is wrong",
    "defect_volume_ratio_right_form": "DEFECT MEASUREMENT - companion of defect_volume_ratio_factor; if the "
                                      "expansion were wrong BOTH would move",
    "ambiguity_mean_dev":             "STRUCTURAL - by construction, and asserted in code: the "
                                      "alternative chain is chosen to match the moments "
                                      "exactly, and the assert fails if it does not. It is "
                                      "reported so a reader can see the match is exact, not "
                                      "approximate. Above-floor companion: "
                                      "ambiguity_first_order_max_gap",
    "ambiguity_var_dev":              "STRUCTURAL - see ambiguity_mean_dev. Companion: "
                                      "ambiguity_first_order_max_gap.",
    "ambiguity_pymrm_vs_product_dev": "STRUCTURAL - the unequal cascade is algebraically the "
                                      "staged product, so the two agree by construction. It "
                                      "guards the `faces` plumbing, and the graded-cells row "
                                      "shows that plumbing moving a real number. Companion: "
                                      "ambiguity_first_order_max_gap",
    "break_min_shift":                "it IS the break table's floor",
}

COVER = {}
for _k in metrics:
    _label = LABELS.get(_k)
    if _label is None:
        COVER[_k] = _rows(_k)                     # built, not written
    else:
        try:
            COVER[_k] = _label + " It is nevertheless moved: " + _rows(_k)
        except KeyError:
            COVER[_k] = _label

assert set(COVER) == set(metrics), set(metrics) ^ set(COVER)
assert all(v.strip() and not v.rstrip().endswith(":") for v in COVER.values())

# --- the break table's own bookkeeping, asserted rather than trusted --------
# 1. every name a row declares is a real key of agreement.json
for _r in BT.itertuples():
    for _m in _r.moves:
        assert _m in metrics, f"break row {_r.defect!r} names unknown metric {_m!r}"
# 2. every row recomputes the metric it is pinned to, and its baseline IS the
#    reported value of that metric - so no row can report a metric against an
#    unrelated number and still claim to cover it
for _r in BT.itertuples():
    assert _r.metric is not None, f"break row {_r.defect!r} pins no metric"
    _clean = metrics[_r.metric]
    assert abs(_r.baseline - _clean) <= 1e-12 * max(1.0, abs(_clean)), (
        f"break row {_r.defect!r}: baseline {_r.baseline!r} is not the reported "
        f"{_r.metric} = {_clean!r}")
# 3. and it actually moves it
    assert _r.shift_abs > 0.0, f"break row {_r.defect!r} moves nothing"

# Every metric below check_agreement.py's ABS_FLOOR is outside CI entirely, so each
# one is given an above-floor companion that IS inside CI and that the same defect
# would move.  Asserted, not claimed.
COMPANION = {
    "e14_1_sigma2_sq_dev": "e14_1_wrong_power_gap, the same chain's null baseline - what "
                           "scaling sigma instead of sigma^2 would give, 25 % away and "
                           "inside CI",
    "e14_1_sigma2_dev": "e14_1_wrong_power_gap",
    "e14_2_L_dev": "e14_2_wrong_power_gap, the same chain's null baseline - the source "
                   "location two hundred miles out, and inside CI",
    "e14_2_two_routes": "e14_2_wrong_power_gap",
    "e14_4_delta_tbar_dev": "e14_4_prefactor_reldev, which is downstream of it",
    "e14_4_delta_var_dev": "e14_4_prefactor_reldev",
    "e14_4_N_dev": "e14_4_prefactor_reldev",
    "e14_4_exponent_dev": "e14_4_prefactor_reldev",
    "e14_4_prefactor_two_routes": "e14_4_prefactor_reldev",
    "fig143_theta_max_max_dev": "fig143_Emax_at_N4",
    "fig143_width_max_dev": "fig143_Emax_at_N4",
    "fig143_Emax_max_reldev": "fig143_Emax_at_N4",
    "fig143_area_max_dev": "eq8_F_vs_integral_max_dev",
    "fig143_mean_max_dev": "eq8_F_vs_integral_max_dev",
    "eq3_variance_identity": "eq8_F_vs_integral_max_dev",
    "eq3_var_quadrature_dev": "eq8_F_vs_integral_max_dev",
    "eq8_F_endpoint_dev": "eq8_F_vs_integral_max_dev",
    "eq8_one_tank_branch_dev": "eq8_F_vs_integral_max_dev",
    "eq4_equal_convolution_dev": "eq4_unequal_shape_gap",
    "pymrm_eq9_max_dev": "outflow_bc_relerr_N4 and eq11_vs_eq9_max_dev",
    "pymrm_eq10_max_dev": "macro_micro_max_gap",
    "rtd_convolution_dev": "rtd_max_dev",
    "defect_fig147_threshold_dev": "defect_fig147_printed_rule_reldev",
    "defect_fig147_from_ch13": "defect_fig147_printed_rule_reldev",
    "ambiguity_mean_dev": "ambiguity_first_order_max_gap",
    "ambiguity_var_dev": "ambiguity_first_order_max_gap",
    "ambiguity_pymrm_vs_product_dev": "ambiguity_first_order_max_gap",
    "break_min_shift": "it is the floor of the table, not a physical quantity, and it is "
                       "itself above the floor: sweep_vs_rootfind_relerr",
}

ABS_FLOOR = 1e-12
below = [k for k, v in metrics.items() if abs(v) < ABS_FLOOR]
_missing = [k for k in below if k not in COMPANION]
assert not _missing, f"below-floor metrics with no named companion: {_missing}"
for _k, _c in COMPANION.items():
    _named = [w.strip(",.") for w in _c.split() if w.strip(",.") in metrics]
    # A companion string that names NO metric makes the loop below vacuous, which
    # is how four of these passed while asserting nothing.  Require a real name.
    assert _named, f"{_k}: companion string names no metric - {_c!r}"
    for _name in _named:
        assert abs(metrics[_name]) >= ABS_FLOOR, f"{_k}: companion {_name} is itself below floor"
print("COVERAGE MAP, asserted key-for-key against agreement.json.  The shifts are\n"
      "MEASURED and each row's baseline is asserted to be the reported metric; the\n"
      "wider association is the row's declared scope, checked name-by-name.\n")
for k in metrics:
    tag = next((t for t in ("STRUCTURAL", "REPORTED, NOT SCORED", "DEFECT MEASUREMENT")
                if COVER[k].startswith(t)), "moved by")
    body = COVER[k].removeprefix(tag + " - ")
    print(f"  {k}\n      {tag}: " + "\n      ".join(textwrap.wrap(body, 82)))
print(f"\n{len(below)} of {len(metrics)} metrics are below ABS_FLOOR = {ABS_FLOOR:g}, i.e. "
      f"OUTSIDE the CI\ncomparison entirely.  Each is an identity, an exact reproduction of "
      f"a printed integer,\nor an exactly-matched construction, and each is given an "
      f"above-floor companion that IS\ninside CI and that the same defect would move.  The "
      f"companion string must NAME a metric:\nan earlier draft had four that named none, so "
      f"their assert passed over an empty set.\n")
for k in below:
    print(f"  {k:<34} companion: {COMPANION[k]}")
print()
report_agreement("A2.4", metrics)'''))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                             "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
Path("index.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print(f"wrote index.ipynb with {len(cells)} cells")
