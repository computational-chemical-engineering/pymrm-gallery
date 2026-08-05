#!/usr/bin/env python3
"""Generate index.ipynb for page A3.3. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# ---------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Danckwerts surface renewal: the age distribution, and what the √D exponent can and cannot settle"
description: "Replace the stagnant film with a surface that is continually torn up and replaced, and k_L becomes √(Ds). The exponent on diffusivity is the falsifiable part — and this paper contains no measurement that could test it."
categories: [sec:A, struct:S4, tier:T0, data:tier6, phase:gas-liquid]
date: 2026-08-05
---

# Danckwerts surface renewal

**Catalog ID:** `A3.3` · **Structure:** `S4` (1-D transient PDE) · **Tier:** T0

The conventional picture of a gas dissolving into an agitated liquid puts a
stagnant film at the interface and lets the solute diffuse across it at steady
state. Danckwerts' objection in 1951 was that the conditions for such a film are
simply absent: turbulence reaches the surface, and the surface of liquid running
over a packing does not keep its identity past the next discontinuity.

His replacement is one sentence of physics. **Elements of liquid are brought to
the surface, absorb for a while by unsteady diffusion, and are swept back into
the bulk — and the chance of an element being replaced does not depend on how
long it has already been there.** That single assumption fixes the distribution
of surface ages, and averaging the penetration flux over it gives

$$k_L = \sqrt{Ds}.$$

Film theory gives $k_L \propto D$; this gives $k_L \propto \sqrt{D}$. That is a
discriminating prediction, and it is the reason the page exists.

**What this page does.** It reproduces the six numerical statements the paper
makes; it solves the transient element problem with pymrm and recovers the
paper's closed forms by the paper's own two routes; it maps where the $\sqrt{D}$
exponent survives and where it does not; it solves the one case Danckwerts says
he could not; and it says plainly that the paper contains no data capable of
testing any of it."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

**The film picture, as this paper states it.** Danckwerts opens by writing the
conventional result, attributing it to Sherwood's *Absorption and Extraction*
(1937):

$$R = \frac{D}{x_L}\,(c^* - c_o) \qquad\text{(eq. 1)},$$

with $x_L$ "the effective film thickness". Measured rates conform to

$$R = k_L (c^* - c_o) \qquad\text{(eq. 2)},$$

and his point is that eq. 2's proven usefulness "in no way provides support for
the stagnant-film hypothesis", because eq. 2 also follows from a picture he
regards as more realistic.

**The surface-renewal picture.** A turbulent liquid is a mass of eddies which
continually expose fresh surface and sweep away surface that has been in contact
with the gas for varying lengths of time. Take the exposed area as unity, let
$s$ be the mean rate of production of fresh surface, and — this is the whole
content — assume *the chance of an element of surface being replaced within a
given time is independent of its age*. Then the fractional replacement rate of
every age group is the same $s$, and the age distribution follows immediately.

**Where $s$ comes from.** Nowhere, in this paper. Danckwerts is explicit: "For
the present $s$ is regarded as a quantity which must be determined
experimentally for any given system." The expressions are therefore "likely to
be of use mainly in comparing rates of absorption in systems of equal $s$" — the
same packing, the same flow, a different solute or a different reaction. Keep
that in mind: it is what makes the $\sqrt{D}$ exponent the testable part and
$s$ itself untestable.

**Scope, and what this page is not.** Whitman's two-film theory is `A3.1` and
Higbie's fixed-contact-time penetration is `A3.2`; both are separate cases with
separate sources, and neither is built here. This page uses eq. 1 only as
Danckwerts prints it, and where his general age-distribution result reduces to a
fixed contact time (his "nonrandom packing" case) it says so *without* attaching
a name — Danckwerts cites Higbie in this paper only for evidence on interfacial
resistance, not for penetration theory. The gallery's other Danckwerts page,
[`A2.1`](../A2.1-danckwerts-boundary-conditions/), is the 1953 inlet/outlet
boundary conditions for a flow reactor: a different paper about a different
problem, and none of its provenance is borrowed here."""))

# ---------------------------------------------------------------- the model
cells.append(md(r"""## The published model

### The age distribution

An element in the age group $\theta \dots \theta + d\theta$ occupies area
$\phi(\theta)\,d\theta$. At steady state the area entering that group in a time
$d\theta$ equals the area in the group below it, less the fraction $s\,d\theta$
of it that is replaced:

$$\phi(\theta)\,d\theta = \phi(\theta - d\theta)\,d\theta\,(1 - s\,d\theta)
\;\Longrightarrow\; \frac{d\phi}{d\theta} = -s\phi,$$

and with $\int_0^\infty \phi\,d\theta = 1$,

$$\boxed{\;\phi(\theta) = s\,e^{-s\theta}\;}\qquad\text{(eq. 4)}.$$

### The average over ages

An element of age $\theta$ absorbs at the stagnant-liquid rate $\psi(\theta)$,
so the mean rate per unit area of surface is the Laplace transform of $\psi$:

$$R = s\int_0^\infty e^{-s\theta}\,\psi(\theta)\,d\theta \qquad\text{(eq. 7)},$$

which Danckwerts notes is much easier to evaluate than $\psi$ itself. For a
general age distribution the same average is

$$R = \int_0^\infty \phi(\theta)\,\psi(\theta)\,d\theta \qquad\text{(eq. 19)}.$$

### The stagnant element

For the simplest case — surface constantly saturated, no reaction — the element
obeys

$$\frac{\partial c}{\partial\theta} = D\frac{\partial^2 c}{\partial x^2},\qquad
c(x,0) = c_o,\quad c(0,\theta) = c^*,\quad c(\infty,\theta) = c_o,$$

with solution and surface flux

$$c = c_o + (c^*-c_o)\,\mathrm{erfc}\!\left[\frac{x}{2\sqrt{\theta D}}\right]
\quad\text{(eq. 27)},\qquad
\psi(\theta) = -D\left(\frac{\partial c}{\partial x}\right)_{x=0}
 = (c^*-c_o)\sqrt{\frac{D}{\pi\theta}} \quad\text{(eq. 3)}.$$

### The results this page reproduces

| # | System | Result |
|---|---|---|
| 8 | surface saturated, no reaction | $R = (c^*-c_o)\sqrt{Ds}$ |
| 6 | eddy-diffusion resistance beneath the surface | $1/k_L = 1/k_E + 1/\sqrt{Ds}$ |
| 9,10,11 | gas-film and/or surface resistance | $R = (c^*-c_o)\big/\left[\tfrac{1}{\sqrt{Ds}} + \tfrac{H}{k_G} + \tfrac{1}{k_S}\right]$ |
| 12,13,34 | first-order reaction, rate constant $r$ | $R = \left[c^* - c_o\frac{s}{r+s}\right]\sqrt{D(r+s)}$, and $= \left[\cdot\right]\sqrt{k_L^2 + Dr}$ |
| 15,16 | instantaneous reaction with a dissolved reagent | $R = \dfrac{c^*\sqrt{Ds}}{\mathrm{erf}[\beta/\sqrt{D}]}$, $\beta$ from eq. 16 |
| 17 | the same, when $D' = D$ | $R = (c^* + c_o')\sqrt{Ds}$ |
| 21 | any age distribution, no other resistance | $R = (c^*-c_o)\sqrt{D}\displaystyle\int_0^\infty\frac{\phi(\theta)}{\sqrt{\pi\theta}}d\theta$ |
| 22 | a patchwork of areas $a$ each with its own $s_a$ | $R = (c^*-c_o)\sqrt{D}\,\dfrac{\sum a\sqrt{s_a}}{\sum a}$ |
| 25 | undetermined $\phi$ **and** a gas film | $R = (c^*-c_o)\dfrac{k_G}{H}\displaystyle\int_0^\infty e^{k_G^2\theta/H^2D}\,\mathrm{erfc}\!\left[\tfrac{k_G}{H}\sqrt{\tfrac{\theta}{D}}\right]\phi(\theta)\,d\theta$ |

Eq. 29 gives the element flux under a surface resistance,
$\psi(\theta) = k(c^*-c_o)\,e^{k^2\theta/D}\,\mathrm{erfc}[k\sqrt{\theta/D}]$
with $k$ standing for $k_S$ or $k_G/H$; eq. 31 gives it under a first-order
reaction,
$\psi(\theta) = c^*\sqrt{Dr}\left[\mathrm{erf}\sqrt{r\theta} + e^{-r\theta}/\sqrt{\pi r\theta}\right]$.

**One printed defect, flagged rather than repaired.** Eq. 25 as printed carries
$e^{k_G^2\theta/HD}$ — $H$ **not** squared — while the $\mathrm{erfc}$ argument
in the same equation is $(k_G/H)\sqrt{\theta/D}$. Eq. 29, which eq. 25 explicitly
cites, fixes the exponent at $k^2\theta/D$ with $k = k_G/H$, so $H$ must be
squared. The table above writes $H^2D$, and *Validation* proves the correction
from the paper's own eq. 9 rather than from dimensional analysis alone."""))

# ---------------------------------------------------------------- parameters
cells.append(md(r"""## Parameters and assumptions

Everything on this page is done in the paper's own CGS units where a physical
value is involved, and in units of $D = s = c^* = 1$ where only a ratio matters.

| symbol | meaning | value used | source |
|---|---|---|---|
| $D$ | diffusivity of the absorbed gas | $10^{-5}$ cm²/s | printed, p. 1461 and p. 1462 |
| $s$ | fractional rate of renewal of surface | 5 s⁻¹ (illustration only) | printed, p. 1464 |
| $\theta$ | age of a surface element | — | — |
| $\mu, \rho$ | viscosity and density of water | 0.01 P, 1.0 g/cm³ | **not printed** — see below |
| $g$ | gravitational acceleration | 981 cm/s² | **not printed** |
| $r, r'$ | first- and second-order rate constants | swept | — |
| $k_S, k_G/H$ | surface and gas-film coefficients | swept | — |

**The two constants that are not in the paper.** The falling-film criterion
$L = lD\mu/(g\rho d^4)$ is stated with a worked value "about 0.05 (for water)",
and the phrase "for water" is where $\mu$ and $\rho$ enter. They are not printed
anywhere in the article, and neither is $g$. Using ordinary CGS water properties
is a reconstruction, and *Validation* reports it as one — but it is a
self-checking reconstruction, because the printed 0.05 is recovered to 2 % and
a wrong exponent on $d$ would move $L$ by four orders of magnitude.

**Assumptions the model rests on, as Danckwerts states them.**

1. The element absorbs exactly as a stagnant, infinitely deep liquid of the same
   age would. He bounds the error of this against a falling laminar film: less
   than 5 % provided $L < 0.1$.
2. The scale of turbulence is much larger than the depth of penetration, so
   velocity gradients within the penetrating layer may be ignored.
3. Replacement probability is independent of age (eq. 4). He devotes a section
   to the two ways this fails — nonrandom packing, and a nonuniform $s$.
4. Beneath every freshly formed surface is liquid at the bulk concentration
   $c_o$, and $c_o$ does not change with time.
5. Where a gas film is present, $k_G$ is constant and Henry's law holds. He
   calls this "almost certainly an oversimplification"."""))

# ---------------------------------------------------------------- the data
cells.append(md(r"""## The data

**There is none, and that is a finding rather than an inconvenience.** The
article has no figures and no tables. Danckwerts says so himself:

> No attempt has been made to compare the expressions derived here with
> published experimental measurements. In spite of the enormous number of these,
> very few are suitable for an analysis of this sort …

and closes the paper by listing two questions "which can be answered only by
experiment". The provenance tier is therefore **6**: everything here is checked
against the paper's own arithmetic and against exact solutions, and **nothing on
this page is validated against a measurement.**

The dataset shipped with the page is the set of sixteen numeric statements in
the running text — assumptions, inputs, definitions and computed results — each
with the sentence it comes from. It is a target list, not data.

**Provenance of the transcription.** `pdfimages -list` reports every page of this
scan as CCITT-G4 bilevel at **300 ppi native**, so the renders were made at that
resolution — rendering higher would be interpolation and would add nothing. The
two constants that the PDF text layer gets wrong were re-read on
nearest-neighbour zoomed crops at the same native resolution, which is where the
stacked "1/100" and the superscript $-5$ are unambiguous.

**No other page's dataset is loaded.** Nothing here is borrowed from `F3.1`,
`F3.5`, `A2.1` or anywhere else, so the cross-page reconciliation rule has
nothing to bite on. The only numbers on the page that are not computed by it are
the rows of this CSV.

**The one row that touches an experiment, and why it cannot help.** The renewal
rate $s \sim 5$ s⁻¹ is inferred from published $k_L$ for CO₂ and water in packed
towers, cited to Perry's *Chemical Engineers' Handbook* (1941). Those $k_L$
values are not printed here and the handbook is not on disk. **That $s$ comes
from inverting $k_L = \sqrt{Ds}$ — the very relation one would want to test — is
an inference, and is flagged as one here**, like $\mu$, $\rho$ and $g$ in the
falling-film check: the paper says only "tentative values of $s$, calculated from
published values of $k_L$", without writing the step down. Eq. 8 is the only
relation in the article connecting the two, and *What pymrm adds* runs the
inversion back to the $k_L$ it implies, which is the right order for a packed
tower — so the inference is safe. It is still an inference. *Validation* makes
the consequence quantitative."""))

# ---------------------------------------------------------------- env cell
cells.append(code(r"""# Colab: install pymrm if it is not already present.
try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm"""))

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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from scipy.linalg import solve_banded
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.sparse import eye_array
from scipy.special import erf, erfc, erfcx, erfcinv
from pymrm import construct_grad, construct_div, NumJac, newton
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A3.3-danckwerts-surface-renewal"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

M = {}                      # every reported metric lands here, once
printed = load_data("danckwerts-1951-printed-numbers.csv", page=PAGE).set_index("id")
meta = load_meta("danckwerts-1951-printed-numbers.csv", page=PAGE)
print(cite_data(meta))
display(printed[["journal_page", "quantity", "printed_value", "units", "kind"]])'''))

cells.append(code(r'''# Pull the paper's own values out of the CSV rather than retyping them.
P = {k: float(printed.loc[k, "printed_value"]) for k in printed.index}
D_w      = P["D_water"]            # cm2/s
th_1h    = 3600.0                  # 1 hour, the exposure the paper quotes
frac_pen = P["depth_fraction"]     # the 1/100 that defines "depth of penetration"
print(f"D  = {D_w:g} cm2/s      depth-of-penetration fraction = {frac_pen:g}")'''))

# ---------------------------------------------------------------- pymrm impl
cells.append(md(r"""## PyMRM implementation

The paper offers two routes to $R$ and uses both: solve the element PDE for
$\psi(\theta)$ and average it over eq. 4 (his cases 1, 2 and 4), or Laplace-
transform the PDE and read $R = -sD(\partial\bar c/\partial x)_{x=0}$ off the
transformed solution (his cases 3 and 6). The closed forms in the table above
are the second route. **This page takes the first**, numerically, so that the two
are independent: pymrm's operators and a time march on one side, `scipy.special`
on the other, sharing no code.

Three implementation choices carry the accuracy, and each is a place the check
could fail.

**The base of the element is a no-flux boundary, not a Dirichlet.** The element
is semi-infinite; truncating it at a finite depth with $c = c_o$ imposed would be
wrong in every case with a reaction, because the *unreacted* gas already
dissolved in a fresh element decays as $c_o e^{-r\theta}$, which is exactly what
Danckwerts' eq. 33 carries as its particular solution $c_o/(r+s)$. A closed base
reproduces that decay by itself; a Dirichlet base does not. (The same class is
used, with a Dirichlet base, to build the *film* limit under *What pymrm adds* —
there it is the physics, not a truncation.)

**The age grid is uniform in $w = \sqrt{\theta}$.** $\psi \sim \theta^{-1/2}$ at
short ages, so a uniform grid in $\theta$ cannot resolve the start and a
geometrically growing one cannot be refined honestly — on a geometric schedule
the step at a given age is set by the growth rate, not by the first step, so
sweeping the first step measures nothing. With $\theta_j = (j\,\Delta w)^2$,
halving $\Delta w$ halves the local step *at every age*, and the quadrature
integrand becomes smooth and finite at the origin.

**The surface flux comes from the discrete mass balance, not from a one-sided
difference.** Backward Euler satisfies accumulation + reaction = surface inflow
exactly, so $\psi$ read that way inherits the scheme's own conservation. The
boundary-gradient value is computed too, and the two are compared in
*Validation*."""))

cells.append(code(r'''def _banded(A):
    """Tridiagonal sparse matrix -> the (1,1) banded layout scipy wants."""
    A = A.tocsr(); n = A.shape[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = A.diagonal(1); ab[1] = A.diagonal(0); ab[2, :-1] = A.diagonal(-1)
    return ab


class Element:
    """One surface element of Danckwerts (1951): transient diffusion, with an
    optional first-order reaction, into liquid exposed to the gas at x = 0.

        dc/dtheta = D d2c/dx2 - r c,      c(x, 0) = c_o.

    Surface x = 0 : c = c*                              (eq. 27's condition), or
                    -D dc/dx = k (c* - c)               (eq. 29's condition).
    Base  x = depth : no flux  ('deep'  - truncated semi-infinite liquid), or
                      c = c_o  ('mixed' - a film over a well-mixed bulk).
    """

    def __init__(self, D=1.0, r=0.0, depth=45.0, n_x=400, stretch=10.0, base="deep"):
        self.D, self.r, self.depth, self.n_x, self.base = D, r, depth, n_x, base
        u = np.linspace(0.0, 1.0, n_x + 1)
        # Smooth exponential stretch: fine at the gas surface, coarse below.
        # Defined by a fixed map of u, so refining n_x refines the WHOLE grid
        # in proportion and second-order accuracy survives the grading.
        self.x_f = (depth * np.expm1(stretch * u) / np.expm1(stretch)
                    if stretch else depth * u)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.shape = (n_x, 1)          # (space, field) - never a bare (n,)
        self.V = np.diff(self.x_f)     # cell volumes per unit area
        self.div = construct_div(self.shape, self.x_f, nu=0)   # nu=0: Cartesian slab
        self.I = eye_array(n_x, format="csr")

    def _bc(self, cstar, c0, k=None, sign=1.0):
        # OUTWARD normal throughout: a dc/dn + b c = d.
        # Left face x = 0: n points in -x, so dc/dn = -dc/dx.
        #   Dirichlet  c = c*                                 -> a=0, b=1,   d=c*
        #   Robin  -D dc/dx = k(c*-c)  <=>  D dc/dn + k c = k c*
        #                                                     -> a=D, b=k,   d=k c*
        if k is None:
            left = {"a": 0.0, "b": 1.0, "d": cstar}
        else:
            left = {"a": self.D, "b": sign * k, "d": k * cstar}
        # Right face x = depth: n points in +x.
        #   no flux (deep liquid)     -> a=1, b=0, d=0
        #   well-mixed bulk (a film)  -> a=0, b=1, d=c_o
        right = ({"a": 1.0, "b": 0.0, "d": 0.0} if self.base == "deep"
                 else {"a": 0.0, "b": 1.0, "d": c0})
        return (left, right)

    def march(self, cstar=1.0, c0=0.0, k=None, th_max=30.0, n_t=800,
              bc_sign=1.0, drop_reaction_in_flux=False, profiles_at=()):
        """March in w = sqrt(theta) with uniform dw. Returns theta, psi, Q."""
        bc = self._bc(cstar, c0, k, bc_sign)
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, bc)
        A_diff = (self.div @ grad).tocsr() * self.D
        b_bc = np.asarray((self.div @ grad_bc).todense()).ravel() * self.D
        if self.base == "deep":
            base_out = lambda cc: 0.0
        else:
            hb = self.x_f[-1] - self.x_c[-1]
            base_out = lambda cc: self.D * (cc[-1] - c0) / hb

        w = np.linspace(0.0, np.sqrt(th_max), n_t + 1)
        th = w ** 2
        c = np.full(self.n_x, float(c0))
        psi = np.zeros(n_t + 1); psi_grad = np.zeros(n_t + 1); Q = np.zeros(n_t + 1)
        psi[0] = psi_grad[0] = np.nan if k is None else k * (cstar - c0)
        prof = {}
        for j in range(1, n_t + 1):
            dt = th[j] - th[j - 1]
            Mx = self.I / dt - A_diff + self.r * self.I
            rhs = c / dt + b_bc
            c_old, c = c, solve_banded((1, 1), _banded(Mx), rhs)
            # psi from the discrete mass balance: backward Euler conserves it
            # exactly, so this needs no boundary-gradient accuracy at all.
            rxn = 0.0 if drop_reaction_in_flux else self.r * (self.V * c).sum()
            psi[j] = (self.V * (c - c_old)).sum() / dt + rxn + base_out(c)
            # ...and, independently, from the boundary machinery of construct_grad.
            psi_grad[j] = -self.D * (grad @ c.reshape(-1, 1) + grad_bc)[0, 0]
            Q[j] = Q[j - 1] + psi[j] * dt
            for t in profiles_at:
                if th[j - 1] < t <= th[j]:
                    prof[t] = (th[j], c.copy())
        return th, psi, Q, psi_grad, prof


def simpson_w(w, g):
    """Simpson's rule on a uniform grid with an even number of intervals."""
    h = w[1] - w[0]
    return h / 3 * (g[0] + 4 * g[1:-1:2].sum() + 2 * g[2:-2:2].sum() + g[-1])


def age_average(th, Q, s):
    """Eq. 7, integrated by parts onto the cumulative uptake Q:

        s int_0^T e^{-s th} psi d th = s e^{-sT} Q(T) + s^2 int_0^T e^{-s th} Q d th,

    then substituting th = w^2. Q(0) = 0 exactly and the integrand 2w e^{-s w^2}Q
    is smooth at the origin, so no analytic small-age input is needed anywhere.
    The boundary term is 6e-13 of the total at the production T = 30/s and is kept
    only so that a DELIBERATELY truncated age integral still evaluates eq. 7
    exactly - which is what makes the corresponding break-table row meaningful.
    """
    w = np.sqrt(th)
    return (s ** 2 * simpson_w(w, 2 * w * np.exp(-s * th) * Q)
            + s * np.exp(-s * th[-1]) * Q[-1])


def R_pymrm(s=1.0, n_t=1600, richardson=True, D=1.0, r=0.0, depth=45.0,
            n_x=400, stretch=10.0, base="deep", **march_kw):
    """R from the pymrm element, with the first-order time error optionally
    Richardson-extrapolated out (the observed temporal order is 1.00)."""
    def one(nt):
        el_ = Element(D=D, r=r, depth=depth, n_x=n_x, stretch=stretch, base=base)
        th_, _, Q_, _, _ = el_.march(n_t=nt, **march_kw)
        return age_average(th_, Q_, s)
    R1 = one(n_t)
    return 2 * one(2 * n_t) - R1 if richardson else R1
'''))

# ---------------------------------------------------------------- results
cells.append(md(r"""## Results

### The element, and the age average

Below: the concentration profile inside one element at four ages (eq. 27), the
surface flux $\psi(\theta)$ it produces (eq. 3), and the age distribution eq. 4
that weights it. The third panel is the point of the whole construction — the
weight $s e^{-s\theta}$ falls away long before $\psi$ does, so the mean rate is
dominated by young surface."""))

cells.append(code(r'''D, s, cstar, c0 = 1.0, 1.0, 1.0, 0.0     # D = s = c* = 1: only ratios matter
el = Element(D=D, r=0.0, depth=45.0, n_x=400, stretch=10.0)
ages = (0.05, 0.25, 1.0, 4.0)
th, psi, Q, psi_grad, prof = el.march(cstar=cstar, c0=c0, th_max=30.0, n_t=1600,
                                      profiles_at=ages)

fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))
for a in ages:
    th_a, c_a = prof[a]
    ax[0].plot(el.x_c, c_a, lw=1.6, label=f"$s\\theta$ = {a:g}")
    ax[0].plot(el.x_c, c0 + (cstar - c0) * erfc(el.x_c / (2 * np.sqrt(D * th_a))),
               "k--", lw=0.8)
ax[0].set(xlim=(0, 6), xlabel="depth $x\\,\\sqrt{s/D}$", ylabel="$c/c^*$",
          title="element profiles (dashed: eq. 27)")
ax[0].legend(fontsize=8)

m = th > 0
ax[1].loglog(th[m], psi[m], lw=1.8, label="pymrm")
ax[1].loglog(th[m], (cstar - c0) * np.sqrt(D / (np.pi * th[m])), "k--", lw=1.0,
             label="eq. 3")
ax[1].set(xlabel="age $s\\theta$", ylabel="$\\psi$", title="element flux")
ax[1].legend(fontsize=8)

tt = np.linspace(0, 5, 400)
ax[2].plot(tt, s * np.exp(-s * tt), lw=1.8, label="$\\phi = s e^{-s\\theta}$ (eq. 4)")
ax[2].plot(tt[1:], s * np.exp(-s * tt[1:]) * np.sqrt(D / (np.pi * tt[1:])), lw=1.8,
           label="contribution to $R$")
ax[2].set(xlabel="age $s\\theta$", ylim=(0, 2.2), title="age distribution and its weight")
ax[2].legend(fontsize=8)
fig.tight_layout(); plt.show()'''))

cells.append(code(r'''# --- the three closed forms, by the pymrm route ---------------------------
rows = []

# eq. 8   surface saturated, no reaction
c0b = 0.3
R8_raw = R_pymrm(s=s, D=D, c0=c0b, cstar=cstar, th_max=30.0, richardson=False)
R8     = R_pymrm(s=s, D=D, c0=c0b, cstar=cstar, th_max=30.0)
e8     = (cstar - c0b) * np.sqrt(D * s)
rows.append(("8", "surface saturated, no reaction", e8, R8, R8 / e8 - 1))

# eq. 10 / 30   surface resistance k_S
kS = 0.7
R10 = R_pymrm(s=s, D=D, c0=c0b, cstar=cstar, k=kS, th_max=30.0)
e10 = (cstar - c0b) / (1 / np.sqrt(D * s) + 1 / kS)
rows.append(("10, 30", f"surface resistance, $k_S$ = {kS}", e10, R10, R10 / e10 - 1))

# eq. 12 / 34   first-order reaction
rr = 2.0
R34 = R_pymrm(s=s, D=D, r=rr, c0=c0b, cstar=cstar, th_max=30.0)
e34 = (cstar - c0b * s / (rr + s)) * np.sqrt(D * (rr + s))
rows.append(("12, 34", f"first-order reaction, $r/s$ = {rr:g}", e34, R34, R34 / e34 - 1))

# eq. 13 must be eq. 12 rewritten with k_L = sqrt(Ds)
kL = np.sqrt(D * s)
e13 = (cstar - c0b * s / (rr + s)) * np.sqrt(kL ** 2 + D * rr)

closed = pd.DataFrame(rows, columns=["eq.", "system", "closed form", "pymrm", "rel. dev."])
M["R_vs_eq8_rel"]  = float(R8 / e8 - 1)
M["R_vs_eq8_raw_rel"] = float(R8_raw / e8 - 1)
M["R_vs_eq30_rel"] = float(R10 / e10 - 1)
M["R_vs_eq34_rel"] = float(R34 / e34 - 1)
M["eq13_vs_eq12"]  = float(e13 / e34 - 1)
display(closed.style.format({"closed form": "{:.6f}", "pymrm": "{:.6f}",
                             "rel. dev.": "{:+.2e}"}))
display(Markdown(
    f"The pymrm element reproduces all three closed forms: "
    f"**{abs(M['R_vs_eq8_rel']):.1e}**, **{abs(M['R_vs_eq30_rel']):.1e}** and "
    f"**{abs(M['R_vs_eq34_rel']):.1e}** relative, with the first-order time error "
    f"Richardson-extrapolated out (raw, at the production grid, eq. 8 sits at "
    f"{abs(M['R_vs_eq8_raw_rel']):.1e}). Eq. 13 is eq. 12 rewritten with "
    f"$k_L = \\sqrt{{Ds}}$ and agrees with it to {abs(M['eq13_vs_eq12']):.1e} — "
    f"an algebraic identity, and it is listed as one below."))'''))

cells.append(md(r"""### Where the $\sqrt{D}$ exponent survives

This is the discriminating content of the paper, so it is worth being precise
about what the model actually predicts. Write

$$n \equiv \frac{d\ln k_L}{d\ln D},$$

so film theory (eq. 1) has $n = 1$ and eq. 8 has $n = \tfrac12$. The interesting
question is not "which is right" — the paper cannot answer that — but **how
fragile the $\tfrac12$ is**, because a prediction that dissolves under any
complication is not much of a prediction."""))

cells.append(code(r'''# n = d ln k_L / d ln D for each of the paper's system types.
def local_exponent(f, D0, h=1e-4):
    return (np.log(f(D0 * np.exp(h))) - np.log(f(D0 * np.exp(-h)))) / (2 * h)

s = 1.0
exp_rows = []
exp_rows.append(("1", "film theory as the paper prints it, $R = (D/x_L)\\Delta c$",
                 local_exponent(lambda D: D / 1.0, 1.0), "exactly 1"))
exp_rows.append(("8", "surface renewal, $\\phi = s e^{-s\\theta}$",
                 local_exponent(lambda D: np.sqrt(D * s), 1.0), "exactly 1/2"))

# eq. 21 with a UNIFORM age distribution (his 'nonrandom packing': vertically
# stacked rings, phi = 1/theta_c on 0 < theta < theta_c, 0 beyond).
thc = 1.0
R21_uniform = lambda D: 2 * np.sqrt(D / (np.pi * thc))
exp_rows.append(("21", "any age distribution (here $\\phi = 1/\\theta_c$, nonrandom packing)",
                 local_exponent(R21_uniform, 1.0), "exactly 1/2, whatever $\\phi$ is"))

# eq. 22, a patchwork of areas each with its own s_a
a_w, s_a = np.array([0.6, 0.3, 0.1]), np.array([0.5, 3.0, 20.0])
R22 = lambda D: np.sqrt(D) * (a_w * np.sqrt(s_a)).sum() / a_w.sum()
exp_rows.append(("22", "nonuniform $s$ over a patchwork of areas",
                 local_exponent(R22, 1.0), "exactly 1/2"))

# eq. 12: in D at fixed r
exp_rows.append(("12", "first-order reaction, $r/s$ = 2, at fixed $r$",
                 local_exponent(lambda D: np.sqrt(D * (2.0 + s)), 1.0), "exactly 1/2"))

# eq. 10/11: a series resistance at the surface
for kk in (0.1, 1.0, 10.0):
    exp_rows.append(("10, 11", f"series surface resistance, $k_S/\\sqrt{{Ds}}$ = {kk:g}",
                     local_exponent(lambda D, kk=kk: 1 / (1 / np.sqrt(D * s) + 1 / kk), 1.0),
                     "$\\tfrac12\\times$(liquid share of the resistance)"))

expo = pd.DataFrame(exp_rows, columns=["eq.", "system", "$n$", "what the algebra gives"])
display(expo.style.format({"$n$": "{:.5f}"}))
M["n_eq8"]  = float(exp_rows[1][2])
M["n_eq21_uniform_phi"] = float(exp_rows[2][2])
M["n_eq22_patchwork"]   = float(exp_rows[3][2])
M["n_eq10_kS1"] = float(exp_rows[6][2])'''))

cells.append(code(r'''# Eq. 25: an undetermined age distribution AND a gas film. Danckwerts states that
# such a system "will show no simple relationship between R and D". Quantify it.
def eq25(D, kG, H, phi, th_max=60.0, H_squared=True):
    kk = kG / H
    a2 = kG ** 2 / ((H ** 2 if H_squared else H) * D)
    f = lambda t: np.exp(min(a2 * t, 700.0)) * erfc(kk * np.sqrt(t / D)) * phi(t)
    return kk * quad(f, 0.0, th_max, limit=400)[0]

# H = 5, k = 0.7 is ONE point, chosen because the printed form diverges there.
# It is not representative: the printed integrand grows like exp{[k^2(H-1)/D-s]t},
# so it converges wherever k^2(H-1) < sD - at H = 5, k = 0.3 it converges to a
# perfectly finite (and wrong) number. The parameter-INDEPENDENT statement is the
# (H, k) table further down, and that is the one the argument rests on.
H = 5.0
kk = 0.7                                  # k = k_G/H, the velocity that matters
phi_exp = lambda t: s * np.exp(-s * t)
phi_uni = lambda t: (1.0 / thc) if t < thc else 0.0

# (i) With the exponential phi, eq. 25 must collapse onto eq. 9/10 exactly.
coll_sq = eq25(1.0, H * kk, H, phi_exp, th_max=200.0, H_squared=True)
eq9 = 1.0 / (1 / np.sqrt(1.0 * s) + H / (H * kk))
M["eq25_collapse_to_eq9"] = float(abs(coll_sq / eq9 - 1))

# ...and with the exponent as PRINTED it cannot collapse onto anything, because
# the integrand does not decay. erfc(k sqrt(th/D)) ~ e^{-k^2 th/D}, so the
# integrand of eq. 25 behaves as exp{[k^2(H-1)/D - s] th}/sqrt(th) when H is not
# squared, against exp(-s th)/sqrt(th) when it is.
growth = kk ** 2 * (H - 1) / 1.0 - s
M["eq25_printed_growth_rate"] = float(growth)
trunc = [(T, eq25(1.0, H * kk, H, phi_exp, th_max=T, H_squared=False)) for T in (10, 20, 40)]
M["eq25_printed_ratio_T40_over_T20"] = float(trunc[2][1] / trunc[1][1])
display(pd.DataFrame(
    [(T, v, v / eq9) for T, v in trunc],
    columns=["age cut-off $s\\theta_{max}$", "eq. 25 with the printed $HD$",
             "ratio to eq. 9"]).style.format({"eq. 25 with the printed $HD$": "{:.4g}",
                                              "ratio to eq. 9": "{:.4g}"}))

# (i-b) THE ARGUMENT THAT DOES NOT DEPEND ON THE POINT CHOSEN. Divergence is a
# property of (H, k); the collapse onto eq. 9 is not. Sweep both.
hk = []
for H_ in (1.5, 2.0, 3.0, 5.0, 10.0):
    for k_ in (0.3, 0.7, 1.5):
        e9 = 1.0 / (1 / np.sqrt(1.0 * s) + 1.0 / k_)
        sq = eq25(1.0, H_ * k_, H_, phi_exp, th_max=200.0, H_squared=True)
        p40 = eq25(1.0, H_ * k_, H_, phi_exp, th_max=40.0, H_squared=False)
        p80 = eq25(1.0, H_ * k_, H_, phi_exp, th_max=80.0, H_squared=False)
        hk.append((H_, k_, k_ ** 2 * (H_ - 1) - s, abs(sq / e9 - 1), p40 / e9, p80 / e9))
hkt = pd.DataFrame(hk, columns=["$H$", "$k_G/H$", "growth rate $k^2(H-1)/D - s$",
                                "$|H^2D$ form / eq. 9 $-1|$",
                                "printed $HD$ / eq. 9, cut at 40",
                                "same, cut at 80"])
display(hkt.style.format({"growth rate $k^2(H-1)/D - s$": "{:+.3f}",
                          "$|H^2D$ form / eq. 9 $-1|$": "{:.1e}",
                          "printed $HD$ / eq. 9, cut at 40": "{:.4g}",
                          "same, cut at 80": "{:.4g}"}).hide(axis="index"))
M["eq25_collapse_worst_over_Hk"] = float(hkt["$|H^2D$ form / eq. 9 $-1|$"].max())
M["eq25_printed_best_dev_over_Hk"] = float(
    (hkt["printed $HD$ / eq. 9, cut at 40"] - 1.0).abs().min())

# (ii) with a non-exponential phi, the exponent is no longer 1/2 and no longer
#      constant: it depends on D itself, which is what "no simple relationship"
#      means quantitatively.
Ds = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
rows25 = []
for kk2 in (0.3, 1.0, 3.0):
    ns = [local_exponent(lambda D, kk2=kk2: eq25(D, H * kk2, H, phi_uni, th_max=thc), Dv)
          for Dv in Ds]
    rows25.append([f"{kk2:g}"] + [f"{v:.3f}" for v in ns])
tab25 = pd.DataFrame(rows25, columns=["$k_G/H$"] + [f"$D$ = {d:g}" for d in Ds])
display(tab25)
M["n_eq25_uniform_phi_min"] = float(local_exponent(
    lambda D: eq25(D, H * 0.3, H, phi_uni, th_max=thc), 10.0))
M["n_eq25_uniform_phi_max"] = float(local_exponent(
    lambda D: eq25(D, H * 3.0, H, phi_uni, th_max=thc), 0.1))

display(Markdown(
    f"**The argument that settles the exponent is the second table, because it does "
    f"not depend on where it is evaluated.** Eq. 9 is eq. 25's own special case for "
    f"the exponential age distribution, so the corrected form must reproduce it at "
    f"*every* $(H, k_G/H)$ — and it does, worst case "
    f"**{M['eq25_collapse_worst_over_Hk']:.1e}** over the whole sweep, which unlike "
    f"the single-point {M['eq25_collapse_to_eq9']:.1e} is above "
    f"`check_agreement.py`'s floor and therefore inside CI. The printed form reproduces it "
    f"**nowhere**: the closest it ever comes is "
    f"{100*M['eq25_printed_best_dev_over_Hk']:.1f} % off, and that is at the mildest "
    f"corner of the sweep. That is a parameter-independent proof, and it is what the "
    f"reading of eq. 25 rests on.\n\n"
    f"The *divergence* is the dramatic symptom, and unlike the collapse it is "
    f"**conditional**. Since $\\mathrm{{erfc}}(k\\sqrt{{\\theta/D}}) \\sim "
    f"e^{{-k^2\\theta/D}}$, the printed integrand behaves as "
    f"$\\exp\\{{[k^2(H-1)/D - s]\\theta\\}}/\\sqrt\\theta$ where the corrected one "
    f"behaves as $e^{{-s\\theta}}/\\sqrt\\theta$, so the printed integral diverges "
    f"**only when $k^2(H-1) > sD$**. At $H = {H:g}$, $k = {kk:g}$ and $D = s = 1$ that "
    f"exponent is **{M['eq25_printed_growth_rate']:+.2f}** and it does diverge — the "
    f"first table shows it growing by a factor of "
    f"{M['eq25_printed_ratio_T40_over_T20']:.2g} when the cut-off is doubled from 20 to "
    f"40, where the corrected form is unchanged in its 13th digit. At $k_G/H$ = 0.3 the "
    f"same printed expression converges perfectly happily, to a finite and wrong "
    f"number; the last two columns of the second table separate the two behaviours. "
    f"So the slip is settled from the paper's own eq. 9 and its own eq. 29, "
    f"not from dimensional analysis and not from picking a divergent corner — and it "
    f"is reported, not silently repaired.\n\n"
    f"The last table is the quantitative version of Danckwerts' remark that a "
    f"system with an undetermined $\\phi$ **and** a gas film shows *no simple "
    f"relationship between $R$ and $D$*. The local exponent ranges from "
    f"{M['n_eq25_uniform_phi_min']:.3f} to {M['n_eq25_uniform_phi_max']:.3f} across the "
    f"sweep — never $\\tfrac12$, never constant, and never a power law at all."))'''))

cells.append(md(r"""**So the exponent is robust to the thing the paper worries about, and fragile
to the thing it does not dwell on.** Every complication *inside* the liquid — a
different age distribution, a patchwork of renewal rates, a first-order reaction
at fixed $r$ — leaves $n = \tfrac12$ untouched, because $\sqrt{D}$ factors
straight out of eqs. 21 and 22. What breaks it is a **resistance in series** with
the liquid: a gas film, an interfacial resistance, or (below) a liquid layer thin
enough to reach steady state. Those pull $n$ down towards 0 or up towards 1, and
in eq. 25 they destroy the power law entirely.

**None of these are identities in the interesting direction.** The rows marked
"exactly 1/2" in the table are algebra and are labelled as such — they cannot
fail, and they cannot be evidence for anything. The rows that move are the ones
that carry information, and what they say is that an exponent measurement is
only interpretable when the liquid-side resistance dominates."""))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

Five independent things are checked: the six numerical statements the paper
makes, the pymrm element against the paper's closed forms, the transcription of
the two equations that could not be checked any other way, the convergence
of the numerics in **both** axes, and — under *What pymrm adds* — a
**pure-discretisation control** on the nonlinear solver.

**What the defect-injection table below cannot do.** Every row of it perturbs an
input and checks that a number responds. That measures *sensitivity*, and a
metric can be sensitive to every defect anyone thinks to inject and still sit at
the wrong value, because the baseline it is measured from was never perturbed by
anything. Discretisation error is the standard way this happens: it is present in
the baseline and in every perturbation of it, so it cancels out of the *response*
and survives in the *value*. Two of this page's headline numbers were
contaminated exactly that way in an earlier draft. The remedy is not another
perturbation, it is a **control** — run the same solver on a case whose answer is
known exactly and read off what it invents. That is what the age-step and grid
refinements below do for the linear element, and what the $c_o'/c^* \to \infty$
control does for the second-order solver.

### 1. Every number the paper prints"""))

cells.append(code(r'''# --- (a) equivalent saturated layer absorbed by stagnant water in 1 hour ---
#     total uptake per area Q = 2 (c*-c0) sqrt(D theta / pi), i.e. a saturated
#     layer of that thickness.
layer_mm = 2 * np.sqrt(D_w * th_1h / np.pi) * 10.0
M["layer_1h_mm"] = float(layer_mm)

# --- (b) the depth-of-penetration coefficient, from eq. 27 -----------------
#     erfc(x / 2 sqrt(D theta)) = 1/100  ->  x = 2 erfcinv(1/100) sqrt(D theta)
coef = 2 * erfcinv(frac_pen)
depth_1h_mm = P["depth_coefficient"] * np.sqrt(D_w * th_1h) * 10.0
M["depth_coefficient"] = float(coef)
M["depth_coefficient_rel"] = float(coef / P["depth_coefficient"] - 1)
M["depth_1h_mm"] = float(depth_1h_mm)

# --- (c) the falling-film criterion L = l D mu / (g rho d^4) ---------------
#     mu, rho and g are NOT printed; ordinary CGS water properties are used and
#     the reconstruction is reported as one.
mu_w, rho_w, g_cgs = 0.01, 1.0, 981.0
L_val = P["L_example_l"] * P["L_example_D"] * mu_w / (g_cgs * rho_w * P["L_example_d"] ** 4)
M["L_falling_film"] = float(L_val)
M["L_falling_film_rel"] = float(L_val / P["L_example"] - 1)

# --- (d) truncating the age integral: sqrt(Ds) erf sqrt(s theta_c) ---------
defic = 1 - erf(np.sqrt(P["sthetac_10pct"]))
sthc_exact = brentq(lambda t: erf(np.sqrt(t)) - 0.90, 0.1, 5.0)
M["erf_deficit_at_sthetac_1p5"] = float(defic)
M["sthetac_for_exactly_10pct"] = float(sthc_exact)

# --- (e) the fraction of R contributed by surfaces older than 0.33 s -------
frac_old = 1 - erf(np.sqrt(P["s_packed_tower"] * P["age_cut_0p33"]))
M["frac_R_from_ages_over_0p33s"] = float(frac_old)

targets = pd.DataFrame([
    ("layer_1h",       "equivalent saturated layer at 1 h (mm)", "about 2", layer_mm),
    ("depth_coefficient", "coefficient in depth = coef sqrt(D theta)", "3.6", coef),
    ("depth_1h",       "depth of penetration at 1 h (mm)", "about 6", depth_1h_mm),
    ("L_example",      "L for the worked falling film", "about 0.05", L_val),
    ("sthetac_10pct",  "deficit of erf sqrt(s theta_c) at s theta_c = 1.5", "< 10 %", defic),
    ("age_cut_0p33",   "share of R from surfaces older than 0.33 s at s = 5", "< 1/10", frac_old),
], columns=["id", "quantity", "as printed", "recomputed"])
display(targets.style.format({"recomputed": "{:.4f}"}))'''))

cells.append(code(r'''display(Markdown(f"""
All six reproduce, and three of them are looser than they look — which is worth
stating rather than rounding away.

* The **saturated layer** comes out at **{layer_mm:.3f} mm** against "about 2 mm",
  {100*(layer_mm/2-1):+.0f} %.
* The **depth-of-penetration coefficient** is **{coef:.4f}** against the printed
  3.6, {100*M['depth_coefficient_rel']:+.2f} %. It is the number that confirms the
  fraction: the sentence defining it prints a stacked "1/100" that the PDF text
  layer drops entirely, and one-tenth would give {2*erfcinv(0.1):.3f} instead.
* The **depth at 1 hour** is **{depth_1h_mm:.2f} mm** using the paper's own 3.6
  ({coef*np.sqrt(D_w*th_1h)*10:.2f} mm using the exact coefficient) against "about
  6 mm" — {100*(depth_1h_mm/6-1):+.0f} %. Both of the paper's millimetre figures
  are **below** the value its own formula gives, so both are roundings *down*:
  {layer_mm:.3f} → "about 2" is ordinary rounding to the nearest millimetre,
  whereas {depth_1h_mm:.2f} → "about 6" is a truncation — rounding would have
  given {round(depth_1h_mm):.0f}. The whole of the {100*(depth_1h_mm/6-1):+.0f} %
  is that truncation; nothing about the calculation is in doubt, since
  {P['depth_coefficient']:g}·√(Dθ) with his own D and exposure leaves no freedom.
* The **falling-film criterion** gives **{L_val:.5f}** against "about 0.05",
  {100*M['L_falling_film_rel']:+.1f} %, with $\\mu$, $\\rho$ and $g$ supplied from
  outside the paper. That is a reconstruction, and it is self-checking twice over:
  the PDF text layer renders the printed $D$ as $10^{{-3}}$, which would give
  {P['L_example_l']*1e-3*mu_w/(g_cgs*rho_w*P['L_example_d']**4):.2f}, and reading the
  $d^{{-4}}$ as $d^{{-2}}$ would give
  {P['L_example_l']*P['L_example_D']*mu_w/(g_cgs*rho_w*P['L_example_d']**2):.2e}.
* The **truncated age average** is short of $\\sqrt{{Ds}}$ by
  **{100*defic:.2f} %** at $s\\theta_c = 1.5$, so "within 10 %" holds; the age at
  which the deficit is exactly 10 % is $s\\theta_c$ = **{sthc_exact:.4f}**, so the
  printed 1.5 is a conservative round number.
* At $s$ = 5 s⁻¹, surfaces older than 0.33 s carry **{100*frac_old:.2f} %** of the
  rate, so "less than one tenth" holds with room to spare.
"""))'''))

cells.append(md(r"""### 2. The one criterion that is a real calculation — eq. 36

Danckwerts' hardest case is the second-order reaction, for which he says "no
solution has so far been found". His way out is to substitute $r'c_o'$ for the
first-order $r$ and use eqs. 12–14, and eq. 36 is the condition under which he
claims that costs less than 10 %. He then illustrates it: with
$r'c_o'/s = 1$, "it is then found that $c_o'/c^* > 50$ if Equation 36 is to be
complied with".

That is a checkable arithmetic claim, and the paper prints the ratio in closed
form. Two independent routes are compared here: the printed closed form, and
direct numerical quadrature of eq. 31 — which shares nothing with it but the
paper."""))

cells.append(code(r'''def eq36_ratio_closed(a, s_, thc_):
    """The ratio printed on p. 1466: erf sqrt((a+s) thc) - sqrt(a/(a+s)) e^{-s thc}
    erf sqrt(a thc), with a = r' c_o' and thc = 0.05/(r' c*)."""
    return (erf(np.sqrt((a + s_) * thc_))
            - np.sqrt(a / (a + s_)) * np.exp(-s_ * thc_) * erf(np.sqrt(a * thc_)))


def eq36_ratio_quad(a, s_, thc_):
    """The same ratio, by integrating eq. 31 numerically. psi1 = c* sqrt(D a)
    [erf sqrt(a th) + e^{-a th}/sqrt(pi a th)], and R1 = c* sqrt(D(a+s))."""
    f = lambda t: s_ * np.exp(-s_ * t) * np.sqrt(a) * (
        erf(np.sqrt(a * t)) + np.exp(-a * t) / np.sqrt(np.pi * a * t))
    num = quad(f, 0.0, thc_, points=[0.0], limit=400)[0]
    return num / np.sqrt(a + s_)


a_, s_ = 1.0, 1.0                        # r' c_o' / s = 1, his illustration
thc_50 = 0.05 * P["so_ratio_threshold"] / a_     # theta_c = 0.05 c_o'/c* / (r'c_o')
r_closed = eq36_ratio_closed(a_, s_, thc_50)
r_quad   = eq36_ratio_quad(a_, s_, thc_50)
thc_star = brentq(lambda t: eq36_ratio_closed(a_, s_, t) - 0.95, 0.5, 8.0)
ratio_star = thc_star * a_ / 0.05

M["eq36_ratio_at_printed_50"] = float(r_closed)
M["eq36_closed_vs_quadrature"] = float(abs(r_closed / r_quad - 1))
M["eq36_threshold_c0p_over_cstar"] = float(ratio_star)
M["eq36_threshold_rel_to_printed"] = float(ratio_star / P["so_ratio_threshold"] - 1)

display(Markdown(f"""
The printed closed form and the direct quadrature of eq. 31 agree to
**{M['eq36_closed_vs_quadrature']:.1e}**, which is what certifies the transcription
of both equations — they were read from the same page but they are not the same
calculation.

At the printed $c_o'/c^* = 50$ the ratio is **{r_closed:.5f}**, and eq. 36 demands
more than 0.95. It is not met. The exact threshold is
$c_o'/c^*$ = **{ratio_star:.2f}** ({100*M['eq36_threshold_rel_to_printed']:+.0f} % on
the printed 50). So "$> 50$" is a rounded statement of "$> 53$", and 50 itself
sits {100*(0.95/r_closed-1):.2f} % on the wrong side of his own criterion. It changes
nothing downstream — the next section shows the approximation is far better than
eq. 36 makes it look — but it is his arithmetic, not ours, and it is recorded.
"""))'''))

cells.append(md(r"""### 3. The two equations that can only be checked by re-deriving them

Eq. 16 defines $\beta$ for the instantaneous-reaction case through a
transcendental equation with four error functions in it, and eq. 25's exponent is
printed inconsistently. Neither can be checked by substituting numbers into
itself. Both are settled here from something else the paper prints.

**Eq. 16 is the flux balance at the reaction plane.** Take the standard
moving-front construction — absorbed gas between the surface and a plane at
$x_f = 2\beta\sqrt{\theta}$, reagent beyond it, both concentrations zero at the
plane — and require that equivalents arrive there at equal rates. That balance,
written out, *is* eq. 16 character for character. This is a derivation, not a
test: it certifies that the equation on the page is the flux balance, and nothing
more. The check with power is in the next section, where a solver that never
forms $\beta$ at all is asked to reach the same answer.

**Eq. 17 follows from 15 and 16 when $D' = D$**, because eq. 16 then reduces to
$\mathrm{erf}(\beta/\sqrt{D}) = c^*/(c^*+c_o')$ and the erf cancels out of eq. 15.
That is an algebraic identity and is labelled one — but it is not powerless
against the *transcription*: swap the erf and erfc in eq. 16 and it stops
holding."""))

cells.append(code(r'''def lam_eq16(cstar_, c0p, D_, Dp_, swap=False):
    """Solve eq. 16 for lambda = beta/sqrt(D), scaled by exp(-lambda^2) so the
    two exponentials never overflow: erfcx(z) = e^{z^2} erfc(z)."""
    def g(lam):
        mu = lam * np.sqrt(D_ / Dp_)
        first  = (cstar_ / np.sqrt(Dp_)) * np.exp(-lam * lam) * erfcx(mu)
        second = (c0p / np.sqrt(D_)) * erf(lam)
        if swap:                       # the deliberate mis-reading: erf <-> erfc
            first  = (cstar_ / np.sqrt(Dp_)) * np.exp(-lam * lam) * np.exp(mu * mu) * erf(mu)
            second = (c0p / np.sqrt(D_)) * erfc(lam)
        return first - second
    return brentq(g, 1e-12, 12.0)


def eq15(cstar_, c0p, D_, Dp_, s_, swap=False):
    lam = lam_eq16(cstar_, c0p, D_, Dp_, swap)
    return cstar_ * np.sqrt(D_ * s_) / erf(lam), lam


ident = []
for c0p in (0.5, 1.0, 3.0, 10.0):
    v, lam = eq15(1.0, c0p, 1.0, 1.0, 1.0)
    ident.append((c0p, v, 1.0 + c0p, v / (1.0 + c0p) - 1, erf(lam), 1 / (1 + c0p)))
idf = pd.DataFrame(ident, columns=["$c_o'/c^*$", "eq. 15", "eq. 17", "rel. dev.",
                                   "erf$(\\beta/\\sqrt{D})$", "$c^*/(c^*+c_o')$"])
display(idf.style.format({"eq. 15": "{:.9f}", "eq. 17": "{:.9f}", "rel. dev.": "{:+.1e}",
                          "erf$(\\beta/\\sqrt{D})$": "{:.6f}",
                          "$c^*/(c^*+c_o')$": "{:.6f}"}))
M["eq17_from_eq15_16"] = float(max(abs(r[3]) for r in ident))
# The swap is INVISIBLE at c0' = c*, where the equation is symmetric in erf/erfc.
# Test it where it can be seen; the row below records the value that is used.
M["eq17_swapped_erf_dev"] = float(abs(
    eq15(1.0, 3.0, 1.0, 1.0, 1.0, swap=True)[0] / 4.0 - 1))
display(Markdown(
    f"Eq. 15 with $\\beta$ from eq. 16 returns eq. 17 to "
    f"**{M['eq17_from_eq15_16']:.1e}** — an identity, listed as one in the break "
    f"table, and blind to any error in the *physics* of eq. 16. It is not blind to "
    f"the transcription: swapping erf and erfc in eq. 16 moves eq. 15 by "
    f"**{100*M['eq17_swapped_erf_dev']:.0f} %** off eq. 17."))'''))

cells.append(md(r"""### 4. Convergence, in both axes

The brief for this repository is explicit that a transient page must refine the
grid *and* the time step and report observed orders, because refining one axis
alone measures the wrong error. Both are refined here on the eq. 8 case, where the
answer is known exactly.

The age grid is uniform in $w = \sqrt{\theta}$, so halving $\Delta w$ halves the
local step at every age — which is the property a geometrically growing schedule
does *not* have. Backward Euler should give first order in that step; the
cell-centred operators on a smoothly stretched grid should give second order in
space. Because the temporal error is the larger of the two at any affordable
setting, the spatial order is measured with the time error Richardson-removed."""))

cells.append(code(r'''ref = np.sqrt(1.0 * 1.0)          # exact sqrt(Ds) with D = s = 1

def R_raw(n_x, n_t):
    th_, _, Q_, _, _ = Element(D=1.0, r=0.0, depth=45.0, n_x=n_x,
                               stretch=10.0).march(th_max=30.0, n_t=n_t)
    return age_average(th_, Q_, 1.0)

rows_t, prev = [], None
for n_t in (200, 400, 800, 1600, 3200):
    e = abs(R_raw(1200, n_t) - ref)
    rows_t.append((n_t, e, np.nan if prev is None else np.log2(prev / e))); prev = e

rows_x, prev = [], None
for n_x in (50, 100, 200, 400, 800):
    e = abs((2 * R_raw(n_x, 1600) - R_raw(n_x, 800)) - ref)   # dt error extrapolated out
    rows_x.append((n_x, e, np.nan if prev is None else np.log2(prev / e))); prev = e

ct = pd.DataFrame(rows_t, columns=["$n_t$", "|err|", "order"])
cx = pd.DataFrame(rows_x, columns=["$n_x$", "|err|", "order"])
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
ax[0].loglog(ct["$n_t$"], ct["|err|"], "o-", label="pymrm")
ax[0].loglog(ct["$n_t$"], ct["|err|"].iloc[0] * (ct["$n_t$"].iloc[0] / ct["$n_t$"]),
             "k--", lw=0.8, label="first order")
ax[0].set(xlabel="$n_t$ (age steps)", ylabel="|R - $\\sqrt{Ds}$|",
          title="time refinement, $n_x$ = 1200"); ax[0].legend(fontsize=8)
ax[1].loglog(cx["$n_x$"], cx["|err|"], "o-", label="pymrm, $\\Delta t$ extrapolated")
ax[1].loglog(cx["$n_x$"], cx["|err|"].iloc[0] * (cx["$n_x$"].iloc[0] / cx["$n_x$"]) ** 2,
             "k--", lw=0.8, label="second order")
ax[1].set(xlabel="$n_x$ (cells)", title="grid refinement"); ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()

M["time_order"]  = float(np.nanmean(ct["order"].values[1:]))
M["space_order"] = float(np.nanmean(cx["order"].values[1:4]))
M["space_err_n400"] = float(cx["|err|"].iloc[3])
M["time_err_n3200"] = float(ct["|err|"].iloc[4])
display(pd.concat([ct, cx], axis=1).style.format({"|err|": "{:.3e}", "order": "{:.2f}"},
                                                  na_rep="—"))
display(Markdown(
    f"Observed orders: **{M['time_order']:.2f}** in the age step (backward Euler, "
    f"first order as expected) and **{M['space_order']:.2f}** in the grid "
    f"(second order, on a grid stretched by a factor $e^{{10}}$ from surface to base). "
    f"The last spatial row flattens because the extrapolated time error has stopped "
    f"being negligible, not because the discretisation has. At the production "
    f"settings the two errors are {M['space_err_n400']:.1e} and "
    f"{M['time_err_n3200']:.1e} — **the time error is the larger one**, which is why "
    f"the reported closed-form agreements have it extrapolated out."))'''))

cells.append(code(r'''# The surface flux, computed two ways, must agree wherever the age is resolved.
th_, psi_, Q_, psig_, _ = Element(D=1.0, r=0.0, depth=45.0, n_x=400,
                                  stretch=10.0).march(th_max=30.0, n_t=1600)
res = th_ > 1e-3                       # ages the grid can resolve at all
M["psi_massbalance_vs_gradient"] = float(np.max(np.abs(psi_[res] / psig_[res] - 1)))
exact_psi = np.sqrt(1.0 / (np.pi * th_[res]))
M["psi_vs_eq3_resolved"] = float(np.max(np.abs(psi_[res] / exact_psi - 1)))
res2 = th_ > 0.1
M["psi_vs_eq3_above_0p1"] = float(np.max(np.abs(psi_[res2] / np.sqrt(1.0 / (np.pi * th_[res2])) - 1)))
display(Markdown(
    f"The mass-balance and boundary-gradient readings of $\\psi$ agree to "
    f"**{M['psi_massbalance_vs_gradient']:.1e}** — an identity of backward Euler, "
    f"not evidence, and listed as one. Against eq. 3 the flux is good to "
    f"**{100*M['psi_vs_eq3_above_0p1']:.2f} %** for $s\\theta > 0.1$ and only "
    f"{100*M['psi_vs_eq3_resolved']:.0f} % once the first few unresolved steps at "
    f"$s\\theta \\sim 10^{{-3}}$ are included. **The pointwise flux at very short ages "
    f"is not converged and no claim rests on it**; what is converged is its age "
    f"average, because $s e^{{-s\\theta}}\\psi$ integrates to a finite value over that "
    f"region and the quadrature runs on the cumulative uptake, which is exactly zero "
    f"at $\\theta = 0$."))'''))

# ---------------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

Three things, in increasing order of how much they are worth.

### 1. The case Danckwerts says he could not solve

Of the second-order reaction — dissolved gas $A$ reacting with dissolved reagent
$B$ at rate $r'cc'$ — he writes: *"These equations are nonlinear and no solution
has so far been found."* He substitutes $r'c_o'$ for $r$, uses the first-order
result, and bounds the error with eq. 36.

pymrm solves the coupled nonlinear pair directly, so the bound can be measured
instead of assumed. `NumJac` handles the pointwise source term with the field
index as the last axis; the transport is the same operators as before."""))

cells.append(code(r'''class TwoSpeciesElement:
    """Danckwerts' case 5: A absorbed at the surface reacts with dissolved B,
    rate r' c c'. He states no solution has been found; this solves it."""

    def __init__(self, D=1.0, Dp=1.0, rp=1.0, depth=45.0, n_x=250, stretch=10.0):
        self.D, self.Dp, self.rp, self.n_x = D, Dp, rp, n_x
        u = np.linspace(0.0, 1.0, n_x + 1)
        self.x_f = depth * np.expm1(stretch * u) / np.expm1(stretch)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.shape = (n_x, 2)                    # (space, field): fields last
        self.V = np.diff(self.x_f)
        self.div = construct_div(self.shape, self.x_f, nu=0)     # nu=0: Cartesian
        self.numjac = NumJac(self.shape)         # pointwise source: last axis in full

    def march(self, cstar=1.0, c0p=50.0, th_max=30.0, n_t=400, tol=None,
              deplete=True):
        # deplete=False freezes the reagent at c_o', which turns the pair into
        # exactly the pseudo-first-order problem with r = r'c_o'. It is the
        # break-table row for the second-order error metric, and it doubles as
        # a pure-discretisation control: the answer is then eq. 12 exactly.
        # OUTWARD normal, per field: A saturated at the surface, B non-volatile;
        # both fields see a closed base, so the deep liquid keeps its reagent.
        bc = ({"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[cstar, 0.0]]},
              {"a": [[1.0, 1.0]], "b": [[0.0, 0.0]], "d": [[0.0, 0.0]]})
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, bc)
        Dv = np.tile([self.D, self.Dp], (self.n_x + 1, 1)).reshape(-1, 1)
        A_diff = (self.div @ grad.multiply(Dv)).tocsc()
        b_bc = np.asarray((self.div @ grad_bc.multiply(Dv)).todense())
        rp, I = self.rp, eye_array(self.n_x * 2, format="csc")
        # Newton converges on the size of its own update, so the tolerance
        # must be scaled by the concentration level or a 50:1 reagent
        # excess asks for more digits than the linear solve carries.
        tol = 1e-10 * max(1.0, cstar, c0p) if tol is None else tol

        def source(cf):
            cc = cf.reshape(self.shape)
            q = rp * cc[:, 0] * (cc[:, 1] if deplete else c0p)
            return np.stack([-q, -q if deplete else 0.0 * q],
                            axis=1).reshape(-1, 1)

        w = np.linspace(0.0, np.sqrt(th_max), n_t + 1)
        th = w ** 2
        c = np.zeros(self.shape); c[:, 1] = c0p
        Q, nits = np.zeros(n_t + 1), []
        for j in range(1, n_t + 1):
            dt = th[j] - th[j - 1]
            c_old = c.reshape(-1, 1).copy()

            def resid(cf):
                cf = cf.reshape(-1, 1)
                g, Jg = self.numjac(source, cf)
                return ((cf - c_old) / dt - (A_diff @ cf + b_bc) - g,
                        I / dt - A_diff - Jg)

            sol = newton(resid, c.reshape(-1, 1), tol=tol, maxfev=60)
            if not sol.success:            # never accept an unconverged step
                raise RuntimeError(f"Newton failed at age {th[j]:.4g}: {sol.message}")
            new = sol.x.reshape(self.shape)
            nits.append(sol.nit)
            # surface flux of A from the discrete balance (the base is closed)
            cB = new[:, 1] if deplete else c0p
            psi = ((self.V * (new[:, 0] - c[:, 0])).sum() / dt
                   + (self.V * rp * new[:, 0] * cB).sum())
            Q[j] = Q[j - 1] + psi * dt
            c = new
        return th, Q, c, nits


def R_second_order(c0p, ratio, s_=1.0, Dp=1.0, n_x=250, n_t=400, cstar=1.0,
                   deplete=True, richardson=True):
    """ratio = r' c_o' / s, the group Danckwerts works in.

    The march is backward Euler, so its error is FIRST ORDER in the age step and
    is Richardson-extrapolated out exactly as R_pymrm does it. Every number this
    function feeds to the page is extrapolated; the raw value at n_t is returned
    alongside so the size of the removed time error is visible. Skipping this
    step puts about a third of the reported second-order error into the age
    grid - see the discretisation control below.
    """
    def one(nt):
        el2 = TwoSpeciesElement(D=1.0, Dp=Dp, rp=ratio * s_ / c0p, depth=45.0,
                                n_x=n_x, stretch=10.0)
        th2, Q2, _, nits = el2.march(cstar=cstar, c0p=c0p, th_max=30.0, n_t=nt,
                                     deplete=deplete)
        return age_average(th2, Q2, s_), nits
    R1, nits1 = one(n_t)
    if not richardson:
        return R1, R1, nits1
    R2, nits2 = one(2 * n_t)
    return 2 * R2 - R1, R1, nits1 + nits2'''))

cells.append(code(r'''import time
t0 = time.time()
sweep, allnits = [], []
R1_pfo = np.sqrt(1.0 * (1.0 + 1.0))          # eq. 12 with r = r'c_o' = s, c_o = 0
for c0p in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0):
    R2v, R2raw, nits = R_second_order(c0p, 1.0)      # extrapolated in the age step
    allnits += list(nits)
    sweep.append((c0p, R2v, R1_pfo, 100 * (R1_pfo / R2v - 1),
                  100 * (R1_pfo / R2raw - 1)))
sw = pd.DataFrame(sweep, columns=["$c_o'/c^*$", "pymrm (exact 2nd order, extrapolated)",
                                  "eq. 12 with $r = r'c_o'$", "error of eq. 12 (%)",
                                  "same, raw at $n_t$ = 400 (%)"])
display(sw.style.format({"pymrm (exact 2nd order, extrapolated)": "{:.5f}",
                         "eq. 12 with $r = r'c_o'$": "{:.5f}",
                         "error of eq. 12 (%)": "{:+.3f}",
                         "same, raw at $n_t$ = 400 (%)": "{:+.3f}"}))

M["so_err_at_printed_50_pct"] = float(sw["error of eq. 12 (%)"].iloc[5])
M["so_err_at_ratio_1_pct"]    = float(sw["error of eq. 12 (%)"].iloc[0])
M["so_err_at_printed_50_raw_pct"] = float(sw["same, raw at $n_t$ = 400 (%)"].iloc[5])
M["so_newton_iters_max"]      = float(np.nanmax(allnits))
print(f"second-order sweep: {time.time()-t0:.1f} s, "
      f"Newton iterations {np.nanmin(allnits):.0f}-{np.nanmax(allnits):.0f} per step")'''))

cells.append(md(r"""**Before reading those numbers, the control that says how much of them is
physics.** Every row of the break table below perturbs an *input* and watches a
number move; none of them can catch a baseline that is wrong rather than
insensitive, because nothing was perturbed to produce it. The error of eq. 12 is
exactly such a baseline — a difference between two numbers, one of which comes
out of a time march.

So run the same solver where the answer is *known*: at $c_o'/c^* \to \infty$ the
reagent cannot be depleted, the problem is exactly pseudo-first-order, and the
answer is exactly $\sqrt{D(r+s)} = \sqrt2$. Whatever the solver reports there is
pure discretisation."""))

cells.append(code(r'''# Pure-discretisation control: c_o'/c* = 1e6, where the exact answer is sqrt(2).
t0 = time.time()
ctrl_rows = []
for n_t in (400, 800, 1600):
    _, Rraw, _ = R_second_order(1e6, 1.0, n_t=n_t, richardson=False)
    ctrl_rows.append((f"{n_t}", Rraw, 100 * (R1_pfo / Rraw - 1)))
Rc_ext, _, _ = R_second_order(1e6, 1.0, n_t=400)        # Richardson from 400/800
ctrl_rows.append(("400, 800 extrapolated", Rc_ext, 100 * (R1_pfo / Rc_ext - 1)))
ct2 = pd.DataFrame(ctrl_rows, columns=["$n_t$", "pymrm", "apparent error (%)"])
display(ct2.style.format({"pymrm": "{:.7f}", "apparent error (%)": "{:+.4f}"}).hide(axis="index"))
M["so_control_raw_n400_pct"] = float(ct2["apparent error (%)"].iloc[0])
M["so_control_extrap_pct"]   = float(ct2["apparent error (%)"].iloc[3])
print(f"discretisation control: {time.time()-t0:.1f} s")
display(Markdown(f"""
At the page's own production setting the solver invents
**{M['so_control_raw_n400_pct']:+.4f} %** of error out of the age step alone, and
extrapolation leaves **{M['so_control_extrap_pct']:+.4f} %**. Every second-order
number on this page is therefore quoted extrapolated: the raw column above is
kept only so the size of what was removed is visible. Read against this control,
{100 * M['so_control_raw_n400_pct'] / M['so_err_at_printed_50_raw_pct']:.0f} % of the
raw {M['so_err_at_printed_50_raw_pct']:.4f} % at $c_o'/c^* = 50$ was backward Euler,
not Danckwerts' substitution.
"""))'''))

cells.append(code(r'''display(Markdown(f"""
**Danckwerts' 10 % bound is correct and enormously conservative.** At his own
illustrative point — $r'c_o'/s = 1$ and $c_o'/c^* = 50$ — the pseudo-first-order
substitution is in error by **{M['so_err_at_printed_50_pct']:+.2f} %**, not 10 %
(raw at $n_t$ = 400 it reads {M['so_err_at_printed_50_raw_pct']:+.2f} %, of which
{M['so_control_raw_n400_pct']:.4f} points is the control's own discretisation).
The 10 % level is not reached until $c_o'/c^*$ falls all the way to about 1,
where the error is still only **{M['so_err_at_ratio_1_pct']:+.2f} %**.

That is not a defect in his reasoning: eq. 36 is a *sufficient* condition built by
bounding a truncated integral by a chain of inequalities, and a chain of bounds is
loose by construction. What the exact solve adds is the size of the looseness —
his criterion demands $c_o'/c^*$ above
**{M['eq36_threshold_c0p_over_cstar']:.0f}** where the exact 10 % point is near 1,
so it is conservative by that factor in the concentration ratio. Read
the other way: the substitution he offers as a last resort is accurate to better
than 1 % over essentially the whole range where "diffusion and reaction play
roughly equal parts", which is the regime he says interest centres on.
"""))'''))

cells.append(md(r"""### 2. An independent test of eqs. 15–17

The instantaneous-reaction result is the one place in the paper where the answer
depends on a transcendental equation with four error functions. The second-order
solver never forms $\beta$, never evaluates an error function, and does not know
eq. 15 exists — so driving $r'c_o'/s$ to infinity and landing on eq. 15 is a
check that **can** fail, unlike the eq. 15 → eq. 17 collapse above.

**Where the three rows are evaluated matters, and one of them had to be moved.**
At $c_o' = c^*$ *and* $D' = D$ eq. 16 is symmetric in erf and erfc: reading it
the wrong way round gives the identical target, so a row there cannot see the
one transcription error this test exists to catch — the same degeneracy the
break table records for its own erf/erfc row. The $D'/D = 1$ row is therefore
run at $c_o'/c^* = 3$ instead, and all three rows below now move under the
swap (the break table prints by how much)."""))

cells.append(code(r'''t0 = time.time()
inst, inst_R = [], {}
for Dp, c0p in ((0.5, 1.0), (1.0, 3.0), (2.0, 1.0)):
    tgt, lam = eq15(1.0, c0p, 1.0, Dp, 1.0)
    R2v, R2raw, _ = R_second_order(c0p, 3000.0, Dp=Dp, n_x=350, n_t=600)
    inst_R[(Dp, c0p)] = (R2v, tgt)
    inst.append((Dp, c0p, lam, tgt, R2v, 100 * (R2v / tgt - 1),
                 100 * (R2raw / tgt - 1)))
inf = pd.DataFrame(inst, columns=["$D'/D$", "$c_o'/c^*$", "$\\beta/\\sqrt{D}$ from eq. 16",
                                  "eq. 15", "pymrm at $r'c_o'/s = 3000$", "dev. (%)",
                                  "same, raw at $n_t$ = 600 (%)"])
display(inf.style.format({"$\\beta/\\sqrt{D}$ from eq. 16": "{:.5f}", "eq. 15": "{:.6f}",
                          "pymrm at $r'c_o'/s = 3000$": "{:.6f}", "dev. (%)": "{:+.3f}",
                          "same, raw at $n_t$ = 600 (%)": "{:+.3f}"}))
M["inst_vs_eq15_worst_pct"] = float(np.max(np.abs(inf["dev. (%)"])))
M["inst_vs_eq15_worst_raw_pct"] = float(np.max(np.abs(inf["same, raw at $n_t$ = 600 (%)"])))

# Is the residue the finite reaction rate, or a disagreement? Push the rate.
Dp_w, c0p_w = 1.0, 3.0                       # the worst row above
tgt_w = inst_R[(Dp_w, c0p_w)][1]
ratio_rows = []
for ratio in (300.0, 3000.0, 10000.0):
    Rv = (inst_R[(Dp_w, c0p_w)][0] if ratio == 3000.0 else
          R_second_order(c0p_w, ratio, Dp=Dp_w, n_x=350, n_t=600)[0])
    ratio_rows.append((ratio, Rv, 100 * (Rv / tgt_w - 1)))
rt = pd.DataFrame(ratio_rows, columns=["$r'c_o'/s$", "pymrm", "dev. from eq. 15 (%)"])
display(rt.style.format({"pymrm": "{:.6f}", "dev. from eq. 15 (%)": "{:+.3f}"}).hide(axis="index"))
M["inst_dev_at_ratio_300_pct"]   = float(rt["dev. from eq. 15 (%)"].iloc[0])
M["inst_dev_at_ratio_10000_pct"] = float(rt["dev. from eq. 15 (%)"].iloc[2])
print(f"{time.time()-t0:.1f} s")
display(Markdown(
    f"Worst deviation **{M['inst_vs_eq15_worst_pct']:.2f} %** across the three rows, "
    f"every value extrapolated in the age step (raw at $n_t$ = 600 the worst row reads "
    f"{M['inst_vs_eq15_worst_raw_pct']:.2f} %, so "
    f"{100*(1 - M['inst_vs_eq15_worst_pct']/M['inst_vs_eq15_worst_raw_pct']):.0f} % of the "
    f"unextrapolated figure was the time grid). The second table is the reason it may be called a "
    f"residue rather than a disagreement: at the worst row the deviation is one-signed "
    f"and falls from {abs(M['inst_dev_at_ratio_300_pct']):.2f} % at $r'c_o'/s$ = 300 to "
    f"{abs(M['inst_dev_at_ratio_10000_pct']):.2f} % at 10000 — it is the finite-rate "
    f"approach to an instantaneous reaction from below. "
    f"Since the two routes share no equation, this is what certifies the reading of "
    f"eq. 16; the eq. 17 identity above certifies only that eq. 16 was copied as a "
    f"flux balance."))'''))

cells.append(md(r"""### 3. Film theory and surface renewal are two ends of one calculation

Danckwerts states a condition for treating a liquid layer of restricted depth as
semi-infinite — "the time of exposure should be so short that the depth of
penetration is less than the depth of the liquid" — and leaves it there. The
element solver can just impose the other boundary: a layer of depth $d$ over a
well-mixed bulk held at $c_o$, with the surface still renewed at rate $s$.

Solving that in Laplace space gives a closed form the paper does not contain,

$$k_L = \sqrt{Ds}\,\coth\!\left(d\sqrt{s/D}\right)
 \;\Longrightarrow\;
 n = \frac{d\ln k_L}{d\ln D} = \frac12 + \frac{z}{\sinh 2z},\qquad z = d\sqrt{s/D},$$

whose two limits are **eq. 8** ($z\to\infty$, $k_L\to\sqrt{Ds}$, $n\to\tfrac12$)
and **eq. 1** ($z\to 0$, $k_L \to D/d$, $n\to 1$, with $x_L = d$). The film
picture is not a rival to surface renewal here; it is what surface renewal
becomes when the renewal is too slow to matter. The pymrm element reproduces the
closed form, and the exponent is the honest way to say where an experiment could
tell the two apart."""))

cells.append(code(r'''z_grid = np.logspace(-1.2, 1.0, 220)
n_closed = 0.5 + z_grid / np.sinh(2 * z_grid)
kL_closed = lambda z: 1.0 / np.tanh(z)          # in units of sqrt(Ds)

check = []
for dep in (0.3, 1.0, 2.0, 5.0):
    R = R_pymrm(s=1.0, D=1.0, depth=dep, base="mixed", stretch=3.0,
                n_x=400, n_t=800, cstar=1.0, c0=0.0, th_max=30.0)
    check.append((dep, R, np.sqrt(1.0) / np.tanh(dep), R / (1 / np.tanh(dep)) - 1))
ck = pd.DataFrame(check, columns=["$d\\sqrt{s/D}$", "pymrm $k_L/\\sqrt{Ds}$",
                                  "$\\coth z$", "rel. dev."])
display(ck.style.format({"pymrm $k_L/\\sqrt{Ds}$": "{:.6f}", "$\\coth z$": "{:.6f}",
                         "rel. dev.": "{:+.1e}"}))
M["kL_coth_worst_rel"] = float(np.max(np.abs(ck["rel. dev."])))

z_075 = brentq(lambda z: 0.5 + z / np.sinh(2 * z) - 0.75, 0.01, 10.0)
z_055 = brentq(lambda z: 0.5 + z / np.sinh(2 * z) - 0.55, 0.01, 10.0)
n_at_36 = 0.5 + P["depth_coefficient"] / np.sinh(2 * P["depth_coefficient"])
M["z_at_n_0p75"] = float(z_075)
M["z_at_n_0p55"] = float(z_055)
M["n_at_z_3p6"]  = float(n_at_36)

fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
ax[0].loglog(z_grid, kL_closed(z_grid), lw=1.8, label="$\\coth z$")
ax[0].loglog(z_grid, 1 / z_grid, "k--", lw=0.8, label="film, $k_L = D/d$")
ax[0].axhline(1.0, color="k", ls=":", lw=0.8)
ax[0].plot(ck["$d\\sqrt{s/D}$"], ck["pymrm $k_L/\\sqrt{Ds}$"], "o", ms=6, label="pymrm")
ax[0].set(xlabel="$z = d\\sqrt{s/D}$", ylabel="$k_L/\\sqrt{Ds}$",
          title="one calculation, two limits"); ax[0].legend(fontsize=8)
ax[1].semilogx(z_grid, n_closed, lw=1.8)
ax[1].axhline(0.5, color="k", ls=":", lw=0.8); ax[1].axhline(1.0, color="k", ls=":", lw=0.8)
ax[1].axvline(P["depth_coefficient"], color="C3", ls="--", lw=1.0)
ax[1].annotate("Danckwerts' own\ndepth-of-penetration rule", (P["depth_coefficient"], 0.78),
               fontsize=8, color="C3", ha="right")
ax[1].set(xlabel="$z = d\\sqrt{s/D}$", ylabel="$n = d\\ln k_L/d\\ln D$",
          ylim=(0.45, 1.05), title="the exponent between the two pictures")
fig.tight_layout(); plt.show()

display(Markdown(f"""
pymrm reproduces $\\coth z$ to **{M['kL_coth_worst_rel']:.1e}** over the whole
crossover. The exponent is half-way between the two theories at
$z$ = **{z_075:.3f}**, and within 10 % of $\\tfrac12$ only beyond
$z$ = **{z_055:.3f}**.

The vertical line is Danckwerts' own criterion, read as a number for the first
time: he defines the depth of penetration as {P['depth_coefficient']:g}$\\sqrt{{D\\theta}}$
and asks that it be smaller than the depth of the liquid, which at the mean age
$\\theta = 1/s$ is exactly $z > {P['depth_coefficient']:g}$. At that $z$ the exponent is
**{n_at_36:.4f}** — his rule of thumb puts $n$ within
{100*(n_at_36/0.5-1):.1f} % of $\\tfrac12$, so it is comfortably conservative for
the purpose he states it for.
"""))'''))

cells.append(md(r"""### And the thing pymrm cannot add: a test

The $\sqrt{D}$ exponent is falsifiable and the paper contains nothing that could
falsify it. It is worth being exact about why, because "no data" understates it.

Danckwerts obtains his one physical number, $s \sim 5$ s⁻¹, from published $k_L$
for CO₂ in water. He writes only that they are "calculated from published values
of $k_L$"; that the calculation is **the inversion of $k_L = \sqrt{Ds}$** is an
inference — eq. 8 is the only relation in the paper joining the two, and the
number below is what it implies — and it is labelled one rather than asserted.
Film theory fits the same
single measurement exactly as well, by inverting $k_L = D/x_L$. Both models have
one free parameter and one datum: **both residuals are zero, for any $k_L$
whatever.** No amount of care with a single system can separate them."""))

cells.append(code(r'''s_paper = P["s_packed_tower"]
kL_implied = np.sqrt(D_w * s_paper)                       # cm/s
xL_film_um = D_w / kL_implied * 1e4                       # micrometres
M["kL_implied_by_s5_cm_s"] = float(kL_implied)
M["xL_film_equivalent_um"] = float(xL_film_um)

# How far apart must two solutes' diffusivities be to separate n = 1 from n = 1/2?
# Two-point exponent estimate: sigma_n = sqrt(2) eps / ln(D2/D1) for a relative
# 1-sigma scatter eps on each k_L. Reject n = 1 in favour of n = 1/2 at 3 sigma.
def ratio_needed(eps, nsig=3.0):
    return float(np.exp(nsig * np.sqrt(2) * eps / 0.5))

need = pd.DataFrame([(e, ratio_needed(e)) for e in (0.02, 0.05, 0.10, 0.20)],
                    columns=["1σ scatter on $k_L$", "$D_2/D_1$ needed for a 3σ separation"])
display(need.style.format({"1σ scatter on $k_L$": "{:.0%}",
                           "$D_2/D_1$ needed for a 3σ separation": "{:.2f}"}))
M["D_ratio_needed_3sigma_5pct"] = ratio_needed(0.05)

display(Markdown(f"""
At $s$ = {s_paper:g} s⁻¹ and $D$ = {D_w:g} cm²/s the implied coefficient is
$k_L$ = **{kL_implied:.3e} cm/s**, and the film thickness that reproduces the very
same number is $x_L$ = **{xL_film_um:.1f} µm**. Two models, two fitted constants,
one measurement, zero residual each.

Separating them needs at least two solutes at the same hydrodynamics — which is
precisely the comparison Danckwerts recommends ("$k_L$ should be proportional to
$\\sqrt{{D}}$ for a series of solutes") and precisely the one he does not carry out.
The table says how demanding it is: at a 5 % scatter on each $k_L$ the two
diffusivities must differ by a factor of **{ratio_needed(0.05):.2f}** before the
exponents are three standard deviations apart. His eq. 6 gives the sharper version
of the same experiment — with an eddy-diffusion resistance beneath the surface,
$1/k_L$ is linear in $1/\\sqrt{{D}}$ under renewal and in $1/D$ under a film — and it
too is left for someone else.

**So this page reproduces a theory and tests its internal arithmetic. It does not
test the theory, and neither did the paper.** Danckwerts closes by saying so: the
two questions he raises "can be answered only by experiment".
"""))'''))

# ---------------------------------------------------------------- break table
cells.append(md(r"""### The break table

Every metric reported below is deliberately broken and the result recorded. Rows
that do not move are kept and labelled — an identity is worth keeping once it is
named as one, and the point of the table is to say what each number can and
cannot see."""))

cells.append(code(r'''breaks = []

# 1. flip the sign of b in the Robin surface condition
Rb = R_pymrm(s=1.0, D=1.0, c0=0.3, cstar=1.0, k=0.7, th_max=30.0, bc_sign=-1.0)
breaks.append(("surface Robin condition: sign of $b$ flipped", "R_vs_eq30_rel",
               f"{M['R_vs_eq30_rel']:+.1e}", f"{Rb/e10-1:+.3e}", "moves"))

# 2. give the element a Dirichlet base at c_o instead of a closed one, WITH reaction
Rd = R_pymrm(s=1.0, D=1.0, r=2.0, c0=0.3, cstar=1.0, base="mixed", depth=45.0,
             stretch=10.0, th_max=30.0)
breaks.append(("reacting element: Dirichlet base $c=c_o$ instead of no-flux",
               "R_vs_eq34_rel", f"{M['R_vs_eq34_rel']:+.1e}", f"{Rd/e34-1:+.3e}",
               "moves — the deep liquid must decay as $c_oe^{-r\\theta}$"))

# 3. forget the reaction term when reading the flux off the mass balance
th3, _, Q3, _, _ = Element(D=1.0, r=2.0, depth=45.0, n_x=400, stretch=10.0).march(
    cstar=1.0, c0=0.3, th_max=30.0, n_t=1600, drop_reaction_in_flux=True)
breaks.append(("surface flux read as accumulation only (reaction dropped)",
               "R_vs_eq34_rel", f"{M['R_vs_eq34_rel']:+.1e}",
               f"{age_average(th3,Q3,1.0)/e34-1:+.3e}", "moves"))

# 4. truncate the age integral at s theta = 1.5, the paper's own cut-off
th4, _, Q4, _, _ = Element(D=1.0, r=0.0, depth=45.0, n_x=400, stretch=10.0).march(
    cstar=1.0, c0=0.0, th_max=1.5, n_t=1600)
breaks.append(("age integral truncated at $s\\theta_c$ = 1.5", "R_vs_eq8_rel",
               f"{M['R_vs_eq8_rel']:+.1e}", f"{age_average(th4,Q4,1.0)-1:+.3e}",
               "moves, to the paper's own $-8.3$ %"))

# 5. shorten the element so the penetration front reaches the base
R5 = R_pymrm(s=1.0, D=1.0, c0=0.0, cstar=1.0, depth=1.5, stretch=3.0, th_max=30.0)
breaks.append(("element truncated at $x\\sqrt{s/D}$ = 1.5 instead of 45",
               "R_vs_eq8_rel", f"{M['R_vs_eq8_rel']:+.1e}", f"{R5-1:+.3e}", "moves"))

# 6. eq. 25 with the exponent as printed
breaks.append(("eq. 25 exponent read as printed, $k_G^2\\theta/HD$",
               "eq25_collapse_worst_over_Hk / eq25_printed_growth_rate",
               f"{M['eq25_collapse_worst_over_Hk']:.1e} / $-s$",
               f"{M['eq25_printed_best_dev_over_Hk']:.3f} at best / "
               f"{M['eq25_printed_growth_rate']:+.2f}",
               "moves at EVERY $(H, k)$ swept; the divergence on top of it is "
               "conditional on $k^2(H-1) > sD$"))

# 7b. eq. 31 with its second term dropped, against the printed eq. 36 ratio
def eq36_quad_broken(a, s_, thc_):
    f = lambda t: s_ * np.exp(-s_ * t) * np.sqrt(a) * erf(np.sqrt(a * t))
    return quad(f, 0.0, thc_, points=[0.0], limit=400)[0] / np.sqrt(a + s_)
breaks.append(("eq. 31 with its $e^{-r\\theta}/\\sqrt{\\pi r\\theta}$ term dropped",
               "eq36_closed_vs_quadrature", f"{M['eq36_closed_vs_quadrature']:.1e}",
               f"{abs(r_closed/eq36_quad_broken(a_, s_, thc_50)-1):.3e}", "moves"))

# 7. eq. 16 with erf and erfc swapped
breaks.append(("eq. 16 with erf and erfc swapped", "eq17_from_eq15_16",
               f"{M['eq17_from_eq15_16']:.1e}", f"{M['eq17_swapped_erf_dev']:.3e}", "moves"))

# 8. the depth-of-penetration fraction read as 1/10
breaks.append(("depth-of-penetration fraction read as 1/10 not 1/100",
               "depth_coefficient", f"{M['depth_coefficient']:.4f}",
               f"{2*erfcinv(0.1):.4f}", "moves"))

# 9. the falling-film D read as the text layer gives it
L_bad = P["L_example_l"] * 1e-3 * mu_w / (g_cgs * rho_w * P["L_example_d"] ** 4)
breaks.append(("$D$ in the $L$ example read as $10^{-3}$ (the PDF text layer)",
               "L_falling_film", f"{M['L_falling_film']:.5f}", f"{L_bad:.3f}",
               "moves by 100×"))

# 9b. the second-order error metric: drop the reagent depletion. Freezing c' at
#     c_o' turns the pair into exactly the pseudo-first-order problem the error
#     is measured against, so the metric must collapse to the numerical floor.
so_frozen, _, _ = R_second_order(50.0, 1.0, deplete=False)
breaks.append(("second-order pair with the reagent frozen at $c_o'$ (no depletion)",
               "so_err_at_printed_50_pct", f"{M['so_err_at_printed_50_pct']:+.4f} %",
               f"{100*(R1_pfo/so_frozen-1):+.4f} %",
               "moves — what is left is the discretisation control's own "
               f"{M['so_control_extrap_pct']:+.4f} %"))

# 9c. the instantaneous-limit test: read eq. 16 with erf and erfc swapped. This
#     is the metric's own target, so the deviation is what has to move.
sw_dev = []
for (Dp_b, c0p_b), (Rv_b, _) in inst_R.items():
    tgt_sw, _ = eq15(1.0, c0p_b, 1.0, Dp_b, 1.0, swap=True)
    sw_dev.append(abs(100 * (Rv_b / tgt_sw - 1)))
breaks.append(("eq. 16 with erf and erfc swapped, as the instantaneous-limit target",
               "inst_vs_eq15_worst_pct", f"{M['inst_vs_eq15_worst_pct']:.3f} %",
               f"{min(sw_dev):.1f}–{max(sw_dev):.1f} %",
               "moves on ALL THREE rows — the range is over the three, and the "
               "smallest of them is the point of moving the middle row off "
               "$c_o' = c^*$, where the swap moved nothing at all"))

# 9d. the coth crossover: give the layer a closed base instead of a well-mixed
#     bulk. Same depth, same solver, wrong physics for the film limit.
d_b = 1.0
R_closed_base = R_pymrm(s=1.0, D=1.0, depth=d_b, base="deep", stretch=3.0,
                        n_x=400, n_t=800, cstar=1.0, c0=0.0, th_max=30.0)
breaks.append((f"finite layer at $z$ = {d_b:g} given a closed base instead of a well-mixed bulk",
               "kL_coth_worst_rel", f"{M['kL_coth_worst_rel']:.1e}",
               f"{R_closed_base*np.tanh(d_b)-1:+.3f}", "moves"))

# 9e. eq. 25's exponent, against the exponents it produces
n25_bad = local_exponent(lambda D: eq25(D, H * 3.0, H, phi_uni, th_max=thc,
                                        H_squared=False), 0.1)
breaks.append(("eq. 25 exponent read as printed, in the $n$ sweep",
               "n_eq25_uniform_phi_max", f"{M['n_eq25_uniform_phi_max']:.4f}",
               f"{n25_bad:.4g}", "moves"))

# --- rows that do NOT move, kept and labelled -----------------------------
# 10. arithmetic instead of harmonic face mean for D
breaks.append(("arithmetic instead of harmonic face mean for $D$", "R_vs_eq8_rel",
               f"{M['R_vs_eq8_rel']:+.1e}", "identical",
               "STRUCTURAL — $D$ is uniform, so there is no jump for the rule to bite on"))
# 11. the age distribution normalisation
norm = quad(lambda t: 1.0 * np.exp(-1.0 * t), 0, np.inf)[0]
mean_age = quad(lambda t: t * 1.0 * np.exp(-1.0 * t), 0, np.inf)[0]
breaks.append(("$\\int\\phi\\,d\\theta$ = 1 and mean age = $1/s$",
               "(not reported)", f"{norm:.6f} / {mean_age:.6f}", "cannot be broken",
               "STRUCTURAL — eq. 4 is *constructed* to satisfy both; they test nothing"))
# 12. mass-balance vs boundary-gradient flux
breaks.append(("$\\psi$ from the mass balance vs from the gradient operator",
               "psi_massbalance_vs_gradient", f"{M['psi_massbalance_vs_gradient']:.1e}",
               "cannot be broken",
               "STRUCTURAL — backward Euler conserves mass exactly; a wrong BC moves both together"))
# 13. eq. 13 vs eq. 12
breaks.append(("eq. 13 against eq. 12", "eq13_vs_eq12", f"{M['eq13_vs_eq12']:.1e}",
               "cannot be broken",
               "STRUCTURAL — eq. 13 *is* eq. 12 with $k_L=\\sqrt{Ds}$ substituted"))
# 14. the four exponents that come out at exactly 1/2 (and the 1/4)
breaks.append(("the $n$ = $\\tfrac12$ exponents of eqs. 8, 21, 22 and the $\\tfrac14$ of eq. 10",
               "n_eq8, n_eq21_uniform_phi, n_eq22_patchwork, n_eq10_kS1",
               f"{M['n_eq8']:.4f} / {M['n_eq21_uniform_phi']:.4f} / "
               f"{M['n_eq22_patchwork']:.4f} / {M['n_eq10_kS1']:.4f}",
               "cannot be broken",
               "STRUCTURAL — these differentiate closed forms, not solutions: "
               "$\\sqrt D$ factors out of eqs. 21 and 22 identically. They record "
               "what the algebra says and are evidence for nothing"))

bt = pd.DataFrame(breaks, columns=["injected defect", "metric it should move",
                                   "as built", "when broken", "verdict"])
display(bt.style.hide(axis="index"))'''))

cells.append(md(r"""**Three warnings about this table, all of which apply to numbers on this page.**

**What a break table cannot do, and what caught it here.** Every row above
perturbs an input and checks that a number moves. That tests *sensitivity*, and
sensitivity is not correctness: a metric can respond correctly to every defect
you inject and still sit at the wrong value, because the baseline itself was
never perturbed by anything. That is exactly what happened to this page. The
second-order error and the instantaneous-limit deviation were first reported
without extrapolating the age step, and a large part of each of them — the two
sections above print how much — was backward Euler rather than physics. No row
in this table could have shown it, because every row moves the *same*
contaminated baseline.
What shows it is a **pure-discretisation control**: run the same solver where
the answer is known exactly and read off what it invents. That control is
printed above (`so_control_raw_n400_pct`, `so_control_extrap_pct`), it is a
reported metric, and the frozen-reagent row is its break-table twin — the two
agree to every digit printed, which is itself the check that the control is
measuring the discretisation of the *production* path and not of some other one.

`check_agreement.py` does not compare any metric whose magnitude is below
$10^{-12}$, so `eq25_collapse_to_eq9`, `eq36_closed_vs_quadrature`, `eq13_vs_eq12` and
`eq17_from_eq15_16` are **outside the regression suite entirely**. A break row is
not a metric and is not in `agreement.json`, so the rows above do not close that
gap either. What closes it is a companion metric well above the floor that comes
from the same calculation: `eq25_collapse_worst_over_Hk` (the same collapse swept
over $(H, k_G/H)$, and the stronger statement anyway) and
`eq25_printed_growth_rate` for eq. 25; `eq36_ratio_at_printed_50` and
`eq36_threshold_c0p_over_cstar`, both read off the same closed form, for eq. 36;
`eq17_swapped_erf_dev` for eq. 17. `eq13_vs_eq12` has no companion and needs
none — it is exactly zero by substitution and is labelled definitional.

And the five structural rows are structural in two different ways. The eq. 4
normalisation, the eq. 13/eq. 12 comparison and the exponents that come out at
exactly $\tfrac12$ are *definitional* — eq. 4 is constructed to normalise, eq. 13
is eq. 12 with $k_L$ substituted, and the exponents differentiate closed forms
that $\sqrt D$ factors out of, so no error of any kind can move any of them,
ever. What they cannot detect is anything at all about the numerics: they would
read the same on a page whose solver was broken. The harmonic face mean and the two readings of
$\psi$ are only *conditionally* structural: the face-mean rule has no jump to bite
on because $D$ is uniform here, and the two flux readings agree because backward
Euler conserves mass exactly. Both would come alive the moment this page acquired
a discontinuous diffusivity or a non-conservative scheme, and both are kept for
that reason."""))

cells.append(code(r'''report_agreement("A3.3", M)'''))

# ---------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**What to copy this page for.** Any problem of the form *"a transient
one-dimensional process runs for a random duration, and I want the average"*.
The machinery is eq. 7 — a Laplace transform of the transient response — and it
does not care what the transient is. `Element` is a general absorbing slab with
selectable surface and base conditions; `age_average` is the transform; the
$w=\sqrt{\theta}$ age grid is the part that makes a $\theta^{-1/2}$ response
integrable without any analytic input.

**Three things to carry across, and one not to.**

1. **Refine in $w$, not in $\theta$.** A geometrically growing time schedule
   cannot be refined honestly, because the step at a given age is fixed by the
   growth ratio rather than by the first step. $\theta_j = (j\Delta w)^2$ has one
   knob and it controls every step.
2. **Integrate the transform by parts onto the cumulative uptake.** $R = s^2\int
   e^{-s\theta}Q\,d\theta$ with $Q(0) = 0$ removes the singularity at the origin
   and needs no small-age asymptote — so the check against the closed form stays
   independent of it.
3. **Close the far boundary, do not pin it.** For a semi-infinite element with a
   volumetric sink, a no-flux base reproduces the correct decaying far field by
   itself; a Dirichlet base at the initial value does not, and the break table
   shows how far wrong it goes.
4. **Do not reuse $s$.** It is not a property of the fluid. Danckwerts says it
   must be measured per system, and everything on this page that touches a
   physical number goes through it.

**Related pages.** [`F3.1`](../F3.1-hatta-regimes/) solves the same
reaction–diffusion competition in a *steady film* and is the natural comparison
for the enhancement factors here; [`F3.5`](../F3.5-co2-amine-absorption/) does it
with a seven-species reversible network. [`A2.3`](../A2.3-taylor-aris-dispersion/)
and [`J1.5`](../J1.5-ldf-breakthrough/) are the other two pages built on a
transient-response-averaged-over-a-distribution structure.
[`A2.1`](../A2.1-danckwerts-boundary-conditions/) is the other Danckwerts page and
is unrelated to this one. `A3.1` (Whitman) and `A3.2` (Higbie) are the two
competitors named in the opening; neither is built, and this page deliberately
does not speak for either — the three-way comparison belongs on whichever of them
is built last, with all three sources in hand.

**Cite the source, not this page:** Danckwerts, P. V., *Significance of
liquid-film coefficients in gas absorption*, Industrial and Engineering Chemistry
**43**(6) 1460–1467 (1951),
[doi:10.1021/ie50498a055](https://doi.org/10.1021/ie50498a055).

Two earlier papers of his are cited by this one and are **not on disk and not
consulted**: *Research* **2**, 494 (1949), for the interfacial-resistance
hypothesis behind eq. 29, and *Trans. Faraday Soc.* **46**, 300 and 701 (1950),
for the solutions quoted as eqs. 31 and 35 and for the second-order error bound
behind eq. 36. Nothing on this page is attributed to them beyond what the 1951
paper itself prints, and where a result is quoted from them the 1951 text is the
transcription source."""))

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
