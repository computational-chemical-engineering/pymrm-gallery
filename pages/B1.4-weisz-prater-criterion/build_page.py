#!/usr/bin/env python3
"""Generate index.ipynb for page B1.4 (the Weisz-Prater criterion). Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "The Weisz-Prater criterion: what an observed rate can tell you"
description: "One inequality in measurable quantities decides whether a rate is kinetics or diffusion. This page computes what its threshold actually guarantees, tests the non-isothermal extension against a false-negative map, and finds the corner of the authors' own parameter range where their claim that the observable removes the multiplicity stops holding."
categories: [sec:B, struct:S3, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-01
---

# The Weisz-Prater criterion: what an observed rate can tell you

**Catalog ID:** `B1.4` · **Structures:** `S3` (1D steady BVP) · **Tier:** T0

Every other pellet page in this gallery runs the problem **forwards**: you know
the rate constant, you compute the effectiveness factor. An experimenter cannot
do that. They have a measured rate, a particle size, a concentration and a
diffusivity, and the question they need answered is whether the number they just
measured is chemistry or transport — *before* they know any kinetics at all.

The Weisz-Prater criterion answers it with one inequality in observables only:

$$\Phi \;=\; \frac{\mathrm{d}N/\mathrm{d}t}{c_0}\,\frac{R^2}{D} \;<\; 1 .$$

That is the inverse problem, and it is the whole content of this page. Three
things are done with it:

1. **What does the threshold actually guarantee?** "Order of magnitude" was all
   1954 could say. Here $\eta$ at $\Phi = 1$ is computed exactly, and it turns
   out to depend strongly on shape and on kinetics.
2. **Where does it fail?** A false-negative map over $\phi_0$, $\beta$ and
   $\gamma$, comparing the isothermal criterion against Weisz and Hicks' own
   non-isothermal extension, their eq. (10a).
3. **Does the observable really remove the multiplicity?** Weisz and Hicks state
   that it does, with one hedge. Over most of their range it does. Where it does
   not is located here, and it reaches **one of the curves their own Fig. 7
   plots**."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

A porous catalyst pellet consumes reactant faster than diffusion can supply it,
so the interior sees less than the surface concentration and the *observed* rate
falls below the intrinsic one. The ratio is the effectiveness factor $\eta$, and
the classical route to it — Thiele's — needs the intrinsic rate constant $k$,
through the Thiele modulus $\phi_0 = R\sqrt{k/D}$.

That is exactly backwards for an experimenter. As Weisz and Hicks put it, the
investigator "is in fact faced with the problem as to whether his catalyst
system reflects true kinetic behaviour": $k$ is what they are trying to measure,
so it cannot be an input to the test of whether the measurement is valid.

The way out is that the *product* $\phi_0^2\eta$ contains no rate constant. It is

$$\phi_0^2\,\eta \;=\; \frac{R^2k}{D}\cdot\frac{\text{observed rate}}{kc_0}
\;=\; \frac{\mathrm{d}N/\mathrm{d}t}{c_0}\,\frac{R^2}{D} \;\equiv\; \Phi ,$$

the observed rate per unit particle volume, the particle radius, the bulk
concentration and the effective diffusivity — four things a laboratory can
measure. Requiring $\Phi < 1$ is then a test that presupposes nothing about the
kinetics. It is the single most-used sanity check in heterogeneous catalysis.

**What it is not is exact.** It is a statement about a *mapping* — measured
$\Phi$ to unknown $\eta$ — asserted at order of magnitude. Two questions follow
immediately and neither is answered in the sources: how far from 1 can $\eta$ be
while $\Phi < 1$ still holds, and is $\eta$ even a function of $\Phi$? This page
answers both by sweeping the forward problem and asking what the criterion would
have concluded."""))

# ------------------------------------------------------------ the published model
cells.append(md(r"""## The published model

### Provenance: what was read, and what was not

**Weisz and Prater (1954), *Advances in Catalysis* 6, 143, is not on disk and
has no open-access route.** It is a book chapter that predates DOIs. It is cited
here as the origin of the result, exactly as Weisz and Hicks' own reference list
gives it on journal page 274 — with no chapter title, which is why none appears
anywhere on this page:

> [5] WEISZ P. B. and PRATER C. D., in *Advances in Catalysis*. Vol. 6, p. 143.
> Academic Press, New York 1954.

Everything transcribed below comes instead from **Weisz and Hicks (1962)**,
which is on disk and is already the source of the published
[`B1.1`/`B1.5`](../B1.1-thiele-weisz-hicks/) and
[`B1.6`](../B1.6-prater-relation/) pages. That paper does far more than mention
the criterion. It prints it as its **eq. (1)**; it states its purpose; it
devotes its **Section V** to the observable formulation; it extends it to the
exothermic non-isothermal case as its **eq. (10a)**; and it tests it against its
own Fig. 2. Every equation below was read off a 600 dpi render
(`pdftoppm -r 600 -f N -l N -png`) of the journal page named, because this scan's
text layer mangles equations badly — it renders eq. (1) as `$:$<I` and eq. (10a)
as `expCrBK1 + LOI`.

### An attribution the page has to be honest about

The inequality is universally called the *Weisz-Prater criterion*. Weisz and
Hicks do not attribute it to reference [5]. They write, on journal page 265:

> It was shown by WEISZ [9] that the conditions … offer a useful and general
> order of magnitude criterion

and their [9] is **WEISZ P. B., *Z. Phys. Chem.* 1957 **11** 1** — Weisz alone,
three years after the review. Reference [5], Weisz and Prater 1954, *is* named
and used in this paper, but for two other things: the insensitivity of $\eta(\phi)$
to particle shape, and the Arrhenius-slope relation reproduced below. So the
catalogue's citation is kept — the criterion carries both names in the
literature and Weisz and Prater is the origin the case names — while the
1957 attribution Weisz and Hicks themselves give is recorded alongside it.
Neither 1954 nor 1957 was consulted.

### The equations, as printed

Journal page 265, eq. (1), with $\mathrm{d}N/\mathrm{d}t$ "the reaction rate per
unit volume of a porous catalyst", $c_0$ the concentration "known externally of
the particle", $D$ the effective diffusivity and $R$ the radius:

$$\frac{\mathrm{d}N}{\mathrm{d}t}\,\frac{1}{c_0}\,\frac{R^2}{D} \;<\; 1 \tag{1}$$

Journal page 269 names the group, and the Fig. 7 caption on page 272 repeats it
verbatim:

$$\eta\phi_0^2 \;=\; \frac{\mathrm{d}N}{\mathrm{d}t}\,\frac{1}{c_0}\,\frac{R^2}{D}
\;=\; \Phi \tag{11}$$

The pellet itself is their eq. (9), journal page 267, written for the sphere:

$$\frac{\mathrm{d}^2y}{\mathrm{d}x^2} + \frac{2}{x}\frac{\mathrm{d}y}{\mathrm{d}x}
= \phi_0^2\,y\,\exp\!\left(\gamma\beta\,\frac{1-y}{1+\beta(1-y)}\right), \tag{9}$$

"where $\phi_0 = R\sqrt{(k_0/D)}$; $y = c/c_0$; $x = r/R$", subject to
"$y(1) = 1$, $(\mathrm{d}y/\mathrm{d}x)_{x=0} = 0$".

For an exothermic reaction the criterion is tightened. Journal page 269, in
terms of the Thiele modulus:

$$\phi_0^2\exp\!\left[\frac{\gamma\beta}{1+\beta}\right] \le 1,
\qquad\text{or}\qquad
\phi_0 \le \exp\!\left[-\tfrac{1}{2}\frac{\gamma\beta}{1+\beta}\right] \tag{10}$$

and then, "reformulated in terms of *observables* by eliminating $k_0$", the
same page gives the form this page is really about:

$$\frac{\mathrm{d}N}{\mathrm{d}t}\,\frac{1}{c_0}\,\frac{R^2}{D}\,
\exp\!\left[\frac{\gamma\beta}{1+\beta}\right] \;<\; 1. \tag{10a}$$

Two further printed statements are used as checks rather than as inputs.
Journal page 268, attributed to WEISZ [5] — that is, to Weisz and Prater 1954:

$$\frac{Q'}{Q} = 1 + \frac{1}{2}\frac{\mathrm{d}\ln\eta}{\mathrm{d}\ln\phi_0},$$

the apparent activation energy an Arrhenius plot returns. And journal page 269:
"For exothermic reactions, with severe thermal effects defined by
$\beta\gamma \gtrsim 5$, there exists a region of $\phi_0$-values … where
multiple solutions result for $\eta = \eta(\phi)$."

### The claim this page is built to test

Journal page 273, Section V, quoted in full so the hedge is not lost:

> The transformation of the abscissa to $\eta\phi_0^2$ has removed the
> multiplicity in the functions, in the sense that a given *observed rate* …
> defines the system uniquely. However, there will be a metastable region of
> rates, indicated by the dashed portions between arrows, which corresponds to
> the regions of multiple solutions in $\eta$ vs. $\phi_0$, which we believe
> cannot be realized by the steady-state catalytic system.

That is a strong, falsifiable statement about a mapping, and it is the reason
the criterion can be used at all: if two different pellet states could produce
the same $\Phi$, no threshold on $\Phi$ could decide between them. The hedge is
about *stability* — states they believe are not reachable — not about
uniqueness, and it is uniqueness that Section 5 below tests. The two are
different claims, and the numerical answer is different for each.

Fig. 7 has four panels, $\gamma = 10, 20, 30, 40$, and each carries its own set
of $\beta$ labels. The largest $\beta$ drawn on each is $0.8$, $0.8$, $0.4$ and
$0.3$ respectively — read off the printed curve labels on a 600 dpi render of
page 272, not digitised. Section 5 uses those four numbers to say exactly which
plotted curves its result reaches."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

The pellet model is Weisz and Hicks' eq. (9), transcribed above — the two
balances of their eqs. (4) and (5) reduced to one equation by the Prater
relation, which is proved and tested on
[`B1.6`](../B1.6-prater-relation/) and is not re-derived here. Written for a
general geometry index:

$$\nabla^2 y \;=\; \phi_0^2\,y\,
\exp\!\left[\frac{\gamma\beta(1-y)}{1+\beta(1-y)}\right],
\qquad \left.\frac{\partial y}{\partial n}\right|_{u=0}=0,\quad y(1)=1,$$

with $y = c/c_0$, $u$ the position scaled on $R$, $\gamma = Q/(R_gT_0)$ the
Arrhenius number, $\beta = (-\Delta H)Dc_0/(\lambda T_0)$ the Prater number, and
$\nabla^2$ carrying the geometry index $\nu = 0,1,2$ for slab, cylinder and
sphere. Weisz and Hicks solve the sphere; $\nu$ is swept here because the
threshold turns out to depend on it.

The effectiveness factor and the observable follow from the same surface flux:

$$\eta = \frac{\nu+1}{\phi_0^2}\left.\frac{\mathrm{d}y}{\mathrm{d}u}\right|_{u=1},
\qquad
\Phi = \phi_0^2\eta = (\nu+1)\left.\frac{\mathrm{d}y}{\mathrm{d}u}\right|_{u=1}.$$

Note what the second expression says: **$\Phi$ does not contain $\phi_0$ at
all.** It is the dimensionless surface flux, and a solver never needs the rate
constant to report it — which is the entire point of the criterion.

Parameter ranges follow the paper: $\gamma = 10$ to $40$, and $\beta$ from
$-0.8$ (endothermic) to $+0.8$ (exothermic); the sweeps here use
$-0.4 \le \beta \le 0.6$, which covers the whole qualitative range at a
solvable cost. Assumptions carried without test: steady state, one reaction,
constant $D$ and $\lambda$, no external film, and $\gamma$ held fixed while
$\phi_0$ varies (which is what the printed $Q'/Q$ relation assumes, since
strictly $\gamma = Q/R_gT_0$ also moves with temperature).

Four isothermal rate laws are used besides the paper's first-order one, all
normalised so that $\mathcal{R}(1) = 1$, which keeps $\phi_0^2$ equal to the
surface rate times $R^2/(Dc_0)$ and therefore keeps $\Phi$ exactly the printed
observable in every case:

$$\mathcal{R} = 1,\quad y^{1/2},\quad y,\quad y^2, \quad
\frac{36\,y}{(1+5y)^2}.$$

The last is a strongly inhibited Langmuir-Hinshelwood form, chosen because its
rate *decreases* with concentration above $y = 0.2$, which is the one isothermal
way to get $\eta > 1$."""))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

**There is none, and there does not need to be — provenance tier 6.**

Neither source reports a measurement that this page could be compared against.
Weisz and Hicks (1962) contains no experimental data at all (the same finding
recorded for [`B1.1`](../B1.1-thiele-weisz-hicks/) and
[`B1.6`](../B1.6-prater-relation/)), and Weisz and Prater (1954) is a review.
**Nothing on this page is validated against experiment**, and the page must not
be read as if it were: what follows is a set of internal consistency checks, a
comparison against exact closed forms, and a comparison between two independent
numerical methods.

No figure is digitised. The paper's Figs. 1-4 and 7 plot exactly the curves
computed here, but they are not needed: the criterion is an *inequality between
computable quantities*, so it can be tested by sweeping the forward problem, and
the two printed numerical statements used as checks — the worked
$\gamma = 20$, $\beta = 0.3$ case and the $\beta\gamma \gtrsim 5$ multiplicity
threshold — are in the text, not in a figure. No CSV is shipped.

Three things *are* read off the 600 dpi render of Fig. 7 on journal page 272,
and none of them is a data point. The $\beta$ labels printed beside the curves,
which bound which curves the paper drew; the pixel positions of panel (d)'s
five major gridlines; and the width in pixels of the printed curve stroke. The
last two exist solely to answer a resolving-power question in Section 5 — *could
the figure have shown the fold this page computes?* — and the answer they give
is negative, that it could not resolve it either way. No curve is traced and no
value is read off an axis, so the tier stays 6."""))

# ---------------------------------------------------------------- environment
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

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import i0, i1
from pymrm import construct_grad, construct_div, NumJac, newton, clip_approach
from gallery_utils import report_agreement

PAGE = "B1.4-weisz-prater-criterion"
np.seterr(all="ignore")
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

Two independent routes to the same pellet, because the page needs both and
because comparing them is the only check here with real resolving power.

**`Pellet`** is the pymrm finite-volume solver. It takes $\phi_0$ as an input
and returns $\eta$ and $\Phi$. Conventions, all of which matter:

- **Boundary conditions use the OUTWARD normal.** At the centre the outward
  normal points inward, so symmetry is $\partial y/\partial n = 0$
  (`{a:1, b:0, d:0}`); at the surface the Dirichlet condition $y = 1$ is
  `{a:0, b:1, d:1}`. Each carries its physical equation in a comment.
- **`nu` in `construct_div` is geometry**: `0` Cartesian slab, `1` cylindrical,
  `2` spherical. It is swept, and Section 3 shows the criterion's guarantee
  moves by a factor of nearly five across it - so getting it wrong is not cosmetic.
- **Constant operators are assembled once** in `__init__`. The Laplacian, its
  boundary contribution and the surface-gradient row never change; only the
  pointwise source block is rebuilt inside Newton.
- **Spatial axis first, field last**: the layout is `(n_u, 1)` even though there
  is only one field. This is not decoration. `NumJac(shape)` couples the **last**
  axis, so `NumJac((n_u,))` on a one-field problem declares every cell coupled to
  every other and builds a dense $n_u \times n_u$ numerical Jacobian — correct,
  but it costs $n_u$ function evaluations per Newton step. Writing the shape as
  `(n_u, 1)` makes the last axis the field axis and the stencil diagonal, for a
  bit-identical answer. The cell below measures the difference; following the
  house layout convention is the whole fix.
- The source clips $y$ into $[0,1]$, which the maximum principle says is where
  the solution lives, so that a Newton iterate overshooting the surface value
  cannot overflow the Arrhenius exponential.

**`shoot`** is the reference: integrate the ODE outward from the centre with
$y(0) = y_c$ prescribed, and stop where $y$ first reaches 1. Substituting
$s = \phi_0 u$ removes $\phi_0$ from the equation, so *the radius at which $y$
reaches 1 is $\phi_0$*. It discretises nothing and it forms no grid; the **only**
thing it shares with the pymrm route is the algebraic rate function, which is
deliberate — a difference there would be a transcription error, not a
discretisation error, and this comparison is meant to measure the latter.

The reason the page needs it is that parametrising by $y_c$ instead of by
$\phi_0$ traverses the **whole** solution branch in one pass, including the
unstable middle branch that no Newton iteration on $\phi_0$ can reach. Every
multiplicity result below rests on that, and it is deterministic: a fixed
parameter ladder, no warm-start continuation chain, so the numbers do not depend
on the path taken — the lesson `B1.1` learned when CI reproduced a different
$\eta$ on the ignited branch.

The branch parameter is $t$ with $y_c = 1 - e^{-t}$, one variable that resolves
$y_c \to 0$ and $y_c \to 1$ logarithmically at both ends."""))

cells.append(code('''NU = {"slab": 0, "cylinder": 1, "sphere": 2}


class Pellet:
    """Weisz-Hicks eq. (9) in pymrm.  Layout (n_u, 1): spatial axis first, field last."""

    def __init__(self, geom="sphere", n_u=400, nu_override=None):
        self.geom = geom
        self.nu = NU[geom] if nu_override is None else nu_override
        self.n_u, self.shape = n_u, (n_u, 1)
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)              # faces
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])        # centres
        # --- boundary conditions, OUTWARD normal ----------------------------
        #   u = 0 (centre) : symmetry,  dy/dn = 0
        #   u = 1 (surface): the external concentration is imposed,  y = 1
        bc = ({"a": 1.0, "b": 0.0, "d": 0.0},
              {"a": 0.0, "b": 1.0, "d": 1.0})
        g, gb = construct_grad(self.shape, self.u_f, self.u_c, bc)
        d = construct_div(self.shape, self.u_f, nu=self.nu)    # nu: 0 slab, 1 cyl, 2 sphere
        self.lap, self.lap_bc = d @ g, (d @ gb).toarray().reshape(-1, 1)
        self.grad, self.grad_bc = g, gb.toarray().reshape(-1, 1)
        # default stencil couples only the LAST axis, which here is the single
        # field, so the block is diagonal - which is right for a pointwise source
        self.numjac = NumJac(self.shape)
        vf = self.u_f ** (self.nu + 1) / (self.nu + 1.0)       # cell volumes for this nu
        self.dv, self.v_tot = np.diff(vf), vf[-1]

    def solve(self, phi, R, y_init=None, maxfev=80):
        p2 = phi * phi
        # clip to the interval the maximum principle guarantees, so that a Newton
        # iterate overshooting y = 1 cannot overflow the Arrhenius exponential
        src = lambda y: -p2 * R(np.clip(y, 0.0, 1.0))

        def residual(y):
            y = y.reshape(self.shape)
            g_s, j_s = self.numjac(src, y)
            return (self.lap @ y.reshape((-1, 1)) + self.lap_bc
                    + np.asarray(g_s).reshape((-1, 1))), self.lap + j_s

        y0 = np.ones(self.shape) if y_init is None else np.asarray(y_init, float)
        r = newton(residual, y0.reshape((-1, 1)), maxfev=maxfev, tol=1e-13,
                   callback=lambda x, g: clip_approach(x, g, 0.0, 1.0))
        y = r.x.ravel()
        res, _ = residual(y)
        return y, float(np.max(np.abs(res)))

    # ------------------------------------------------------------ observables
    def Phi(self, y):
        """The OBSERVABLE modulus, straight from the surface flux. No phi needed."""
        f = self.grad @ y.reshape((-1, 1)) + self.grad_bc
        return (self.nu + 1.0) * float(f[-1].item())

    def eta(self, y, phi):
        return self.Phi(y) / phi**2

    def Phi_from_volume(self, y, phi, R):
        """The same observable as a volume average of the intrinsic rate.

        NOTE this is NOT an independent route. Multiplying the converged residual
        by the cell volumes and summing telescopes the divergence, so
        Phi(flux) - Phi(volume) IS the volume-weighted mean Newton residual. It is
        therefore sensitive to a non-converged solve and to a mismatch between the
        rate in the model and the rate used here, and structurally blind to the
        geometry index and to the grid. Section 6 measures exactly that.
        """
        return phi * phi * float(np.sum(R(np.clip(y, 0.0, 1.0)) * self.dv) / self.v_tot)


# --------------------------------------------------------------- the rate laws
def wh_rate(beta, gamma, scale=1.0, sgn=1.0):
    """Weisz-Hicks eq. (9) kinetics.  `scale`/`sgn` exist only for Section 6."""
    return lambda y: scale * y * np.exp(sgn * gamma * beta * (1 - y) / (1 + beta * (1 - y)))


# All isothermal laws are normalised to R(1) = 1 so that phi_0^2 is the surface
# rate x R^2/(D c_0) in every case, and Phi is therefore exactly the printed
# observable of eq. (11) in every case.
KIN = {"zero order":         lambda y: np.ones_like(np.asarray(y, float)),
       "half order":         lambda y: np.sqrt(y),
       "first order":        lambda y: y,
       "second order":       lambda y: y * y,
       "inhibited LH (K=5)": lambda y: 36.0 * y / (1.0 + 5.0 * y) ** 2}


# ------------------------------------------------------------------- shooting
def shoot(nu, R, t, rtol=1e-10, dense=False):
    """y_c = 1 - exp(-t); integrate out to y = 1.  Returns (phi, eta, Phi[, sol]).

    Shares only the algebraic rate law with `Pellet` - no operator, no grid, no
    Newton iteration. Substituting s = phi_0 u removes phi_0 from the equation,
    so the radius at which y reaches 1 IS phi_0.
    """
    y_c = -np.expm1(-t)
    r0 = float(R(y_c))
    if not (0.0 < y_c < 1.0) or r0 <= 0.0:
        return None
    # start on the series solution y = y_c + R s^2 / (2(nu+1)), close enough that
    # the offset is 1e-10 of the distance to the nearer of the two bounds
    s0 = min(1e-3, np.sqrt(2 * (nu + 1) * 1e-10 * min(y_c, np.exp(-t)) / r0))
    rhs = lambda s, v: [v[1], float(R(np.clip(v[0], 0.0, 1.0))) - (nu / s) * v[1]]
    hit = lambda s, v: v[0] - 1.0
    hit.terminal, hit.direction = True, 1
    v0 = [y_c + r0 * s0**2 / (2 * (nu + 1)), r0 * s0 / (nu + 1)]
    sol = solve_ivp(rhs, (s0, 1e5), v0, events=hit, rtol=rtol, atol=1e-16,
                    method="DOP853", dense_output=dense)
    if not sol.t_events[0].size:
        return None
    s = float(sol.t_events[0][0])                 # s at y = 1  IS  phi_0
    yp = float(sol.y_events[0][0][1])
    out = (s, (nu + 1.0) * yp / s, (nu + 1.0) * s * yp)
    return out + (sol,) if dense else out


class Branch:
    """One solution branch, traced once and queried many times.

    Parametrised by the centre concentration, not by phi_0, so the trace covers
    the whole branch - upper, middle and lower - in one deterministic pass.
    """

    def __init__(self, nu, R, n=110, tmin=1e-9, tmax=24.0, rtol=1e-9):
        self.nu, self.R, self.rtol = nu, R, rtol
        ts = np.logspace(np.log10(tmin), np.log10(tmax), n)[::-1]   # mild -> severe
        rows = [(t, shoot(nu, R, t, rtol)) for t in ts]
        rows = [(t, o) for t, o in rows if o is not None and np.all(np.isfinite(o))]
        self.t = np.array([r[0] for r in rows])
        self.phi = np.array([r[1][0] for r in rows])
        self.eta = np.array([r[1][1] for r in rows])
        self.Phi = np.array([r[1][2] for r in rows])

    def all_where(self, f):
        """EVERY state along the branch where f(phi, eta, Phi) changes sign.

        Returning all of them, not the first, is deliberate: where the observable
        is multivalued there is more than one, and reporting only the first would
        hide exactly the failure this page is looking for.
        """
        v = np.array([f(p, e, P) for p, e, P in zip(self.phi, self.eta, self.Phi)])
        g = lambda t: (lambda o: np.nan if o is None else f(*o))(
            shoot(self.nu, self.R, t, self.rtol))
        out = []
        for i in range(len(v) - 1):
            if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i] * v[i + 1] < 0:
                ts = brentq(g, self.t[i], self.t[i + 1], xtol=1e-13, rtol=1e-14)
                out.append((ts,) + shoot(self.nu, self.R, ts, self.rtol))
        return out                       # list of (t, phi, eta, Phi)

    def where(self, f):
        """The first such state as (phi, eta, Phi), or None."""
        r = self.all_where(f)
        return r[0][1:] if r else None

    def seed(self, t, p):
        """A pymrm initial guess on `p`'s grid, taken from the branch state at t."""
        o = shoot(self.nu, self.R, t, self.rtol, dense=True)
        phi, sol = o[0], o[3]
        return phi, np.clip(sol.sol(np.clip(p.u_c * phi, sol.t[0], phi))[0], 0.0, 1.0)


def fold_depth(x):
    """Largest relative backward excursion of a sequence traced in increasing order.

    Zero means strictly monotone, i.e. the quantity determines the state uniquely.
    """
    run = np.maximum.accumulate(x)
    return float(((run - x) / run).max())


def fold_turning_points(br):
    """Every local extremum of Phi along the branch, LOCATED rather than described.

    A fold in the observable is a maximum of Phi followed by a minimum; between the
    two, three states share one measured Phi. Where those two values sit relative to
    Phi = 1 is the whole question of Section 5, so they are computed: each extremum
    bracketed on the traced branch is refined with a bounded Brent search on the
    shooting solver itself. The refined depth comes out slightly LARGER than the
    traced table's, which is the traced table's spacing, not a discrepancy.

    Returns a list of ("max"|"min", t, phi_0, eta, Phi) in branch order.
    """
    P, t = br.Phi, br.t
    out = []
    for i in range(1, len(P) - 1):
        if (P[i] - P[i - 1]) * (P[i + 1] - P[i]) >= 0:
            continue
        lo, hi = sorted((float(t[i - 1]), float(t[i + 1])))
        sgn = -1.0 if P[i] > P[i - 1] else 1.0          # a maximum: minimise -Phi
        r = minimize_scalar(lambda tt: sgn * shoot(br.nu, br.R, tt, br.rtol)[2],
                            bounds=(lo, hi), method="bounded",
                            options={"xatol": (hi - lo) * 1e-6})
        out.append(("max" if sgn < 0 else "min", float(r.x))
                   + shoot(br.nu, br.R, float(r.x), br.rtol))
    return out


def eta_exact(nu, phi):
    """Closed-form isothermal effectiveness factor."""
    if nu == 0:
        return np.tanh(phi) / phi
    if nu == 1:
        return 2 * i1(phi) / (phi * i0(phi))
    return (3.0 / phi) * (1.0 / np.tanh(phi) - 1.0 / phi)


print("operators are assembled once per Pellet instance; nothing is rebuilt inside Newton")'''))

cells.append(code('''# --- the layout convention, measured (timings are machine dependent) ----------
import time as _time


def _build_and_solve(shape):
    t0 = _time.time()
    u_f = np.linspace(0.0, 1.0, shape[0] + 1)
    u_c = 0.5 * (u_f[:-1] + u_f[1:])
    bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 0.0, "b": 1.0, "d": 1.0})
    g, gb = construct_grad(shape, u_f, u_c, bc)
    d = construct_div(shape, u_f, nu=2)
    lap, lap_bc = d @ g, (d @ gb).toarray().reshape(-1, 1)
    nj = NumJac(shape)
    t_build = _time.time() - t0
    src = lambda y: -np.clip(y, 0.0, 1.0)

    def res(y):
        y = y.reshape(shape)
        gs, js = nj(src, y)
        return (lap @ y.reshape((-1, 1)) + lap_bc + np.asarray(gs).reshape((-1, 1))), lap + js

    t0 = _time.time()
    r = newton(res, np.ones(shape).reshape((-1, 1)), maxfev=80, tol=1e-13,
               callback=lambda x, g: clip_approach(x, g, 0.0, 1.0))
    rr, _ = res(r.x.ravel())
    return t_build, _time.time() - t0, float(np.max(np.abs(rr))), r.x.ravel()


N_DEMO = 300
tb0, ts0, rr0, y0_ = _build_and_solve((N_DEMO,))
tb1, ts1, rr1, y1_ = _build_and_solve((N_DEMO, 1))
print(f"n_u = {N_DEMO}, phi_0 = 1, first order, sphere")
print(f"  shape (n_u,)  : build {tb0:7.3f} s  solve {ts0:7.3f} s  residual {rr0:.2e}")
print(f"  shape (n_u, 1): build {tb1:7.3f} s  solve {ts1:7.3f} s  residual {rr1:.2e}")
print(f"  max |y difference| between the two: {np.max(np.abs(y0_ - y1_)):.2e}")
print("Same answer; the last-axis stencil is what differs. Everything below uses (n_u, 1).")'''))

# --------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. The transcription, checked before anything is built on it

Three checks, in increasing strength. The first two are algebra and cannot fail
in any interesting way; they are here because a mis-read exponent would break
them, which is the only thing they are claimed to test."""))

cells.append(code('''# ---- check A: eq. (11) is a definition chase, and it must close --------------
# The observed rate per unit volume is the flux through the surface divided by the
# volume:  dN/dt = (S/V) D (dc/dr)|_R = (nu+1) D c_0 /R^2 * (dy/du)|_1.
# Multiply by R^2/(D c_0) and eq. (11) says that equals phi_0^2 eta.
nu_s, phi_s, dyd = sp.symbols("nu phi y1", positive=True)
Phi_from_rate = (nu_s + 1) * dyd                    # dN/dt * R^2/(D c_0)
eta_sym = (nu_s + 1) * dyd / phi_s**2               # the paper's eq. (2)/(3)
print("check A  Phi - phi^2 eta  =", sp.simplify(Phi_from_rate - phi_s**2 * eta_sym),
      "  (eq. 11 closes identically; this tests the bookkeeping, nothing else)")

# ---- check B: what the step from (10) to (10a) assumes -----------------------
# (10) constrains phi_0^2; (10a) constrains Phi = phi_0^2 eta. The two are the same
# inequality only where eta = 1, so (10a) is (10) tightened by a factor eta. The
# paper's own argument is a bootstrap: IF the criterion holds THEN eta ~ 1, so the
# substitution is self-consistent. The size of that assumption is what Section 4
# measures, and it is not zero.
g_s, b_s, eta_s = sp.symbols("gamma beta eta", positive=True)
lhs_10 = phi_s**2 * sp.exp(g_s * b_s / (1 + b_s))
lhs_10a = (phi_s**2 * eta_s) * sp.exp(g_s * b_s / (1 + b_s))
print("check B  eq.(10a) LHS / eq.(10) LHS =", sp.simplify(lhs_10a / lhs_10),
      " -> the two coincide exactly when eta = 1, and (10a) is stricter when eta > 1")

# ---- check C: the paper's own worked number ---------------------------------
# Journal page 269: 'a case like gamma = 20, beta = 0.3 becomes conservatively
# "safe" at phi_0 ~ 0.1 as predicted by (10) above.'
GAM_W, BET_W = 20.0, 0.3
PHI_W = float(np.exp(-0.5 * GAM_W * BET_W / (1 + BET_W)))
print(f"\\ncheck C  eq.(10) gives phi_0 <= exp[-0.5*gamma*beta/(1+beta)] = {PHI_W:.6f}")
print(f"         the paper states 'phi_0 ~ 0.1' for this case  ->  agrees to "
      f"{abs(PHI_W - 0.1) / 0.1 * 100:.1f} % of the quoted 0.1")
print("         A mis-read exponent breaks this: gamma*beta (no 1+beta) would give "
      f"{np.exp(-0.5 * GAM_W * BET_W):.4f}, and 0.5*gamma*beta/(1-beta) would give "
      f"{np.exp(-0.5 * GAM_W * BET_W / (1 - BET_W)):.4f}.")'''))

cells.append(md(r"""### 2. Two solvers, and what each is worth

The shooting reference is checked against the closed-form isothermal
effectiveness factor first, because everything downstream leans on it. Then
pymrm is checked against the same closed form, and against the shooting solution
in the non-isothermal case. **These are the only comparisons on this page with
resolving power**, and they are the ones that must converge."""))

cells.append(code('''# --- the shooting reference against the exact isothermal eta -----------------
shoot_ref = []
for nu in (0, 1, 2):
    for t in (0.01, 0.1, 0.5, 1.0, 3.0, 8.0, 15.0):
        o = shoot(nu, KIN["first order"], t, rtol=1e-11)
        if o:
            shoot_ref.append(abs(o[1] - eta_exact(nu, o[0])) / eta_exact(nu, o[0]))
SHOOT_EXACT = float(max(shoot_ref))
print(f"shooting vs the closed form, 3 geometries x 7 states: {SHOOT_EXACT:.2e} "
      f"worst relative")

# --- pymrm against the exact isothermal eta, under refinement ----------------
PHIS = (0.1, 0.5, 1.0, 3.0, 10.0)
pm_exact = {}
print(f"\\n{'n_u':>6} {'slab':>11} {'cylinder':>11} {'sphere':>11}")
for n_u in (100, 200, 400):
    row = []
    for geom, nu in (("slab", 0), ("cylinder", 1), ("sphere", 2)):
        p = Pellet(geom, n_u)
        e = [abs(p.eta(p.solve(f, KIN["first order"])[0], f) - eta_exact(nu, f))
             / eta_exact(nu, f) for f in PHIS]
        row.append(max(e))
    pm_exact[n_u] = max(row)
    print(f"{n_u:6d} " + " ".join(f"{v:11.3e}" for v in row))
PM_EXACT = float(pm_exact[400])
PM_ORDER = float(np.log2(pm_exact[200] / pm_exact[400]))
print(f"\\nworst at n_u = 400: {PM_EXACT:.3e}; observed order 200 -> 400: {PM_ORDER:.2f}")'''))

cells.append(code('''# --- pymrm against shooting, non-isothermal, over the whole branch -----------
# The shooting profile seeds Newton, so pymrm can be put on ANY branch state -
# including the unstable middle one. The comparison is then at equal phi_0, which
# is well conditioned; comparing at equal centre concentration is not, because on
# a steeply ignited branch eta is exponentially sensitive to y_c and the pymrm
# centre value sits at u = 1/(2 n_u), not at u = 0.
CASES = ((0.2, 20.0), (0.6, 20.0), (0.4, 30.0), (-0.4, 30.0))
TS = (0.02, 0.3, 2.0, 6.0, 12.0)
REF_STATES, SKIPPED = [], []
for beta, gamma in CASES:
    for t in TS:
        o = shoot(2, wh_rate(beta, gamma), t, rtol=1e-11, dense=True)
        if o is None:
            continue
        (REF_STATES if o[0] <= 30.0 else SKIPPED).append((beta, gamma, o))
n_states = len(REF_STATES)

pm_shoot, pm_shoot_rn = {}, {}
for n_u in (200, 400, 800):
    p = Pellet("sphere", n_u)
    worst, worst_rn = 0.0, 0.0
    for beta, gamma, (phi, eta_r, Phi_r, sol) in REF_STATES:
        y0 = np.clip(sol.sol(np.clip(p.u_c * phi, sol.t[0], phi))[0], 0.0, 1.0)
        y, rn = p.solve(phi, wh_rate(beta, gamma), y_init=y0)
        # an order-2 result is only meaningful if every state on it is CONVERGED;
        # a seeded solve that quietly failed to iterate would fake one
        worst_rn = max(worst_rn, rn)
        worst = max(worst, abs(p.eta(y, phi) - eta_r) / eta_r)
    pm_shoot[n_u], pm_shoot_rn[n_u] = worst, worst_rn
PM_SHOOT = float(pm_shoot[800])
PM_SHOOT_RN = float(max(pm_shoot_rn.values()))
PM_SHOOT_ORDER = float(np.log2(pm_shoot[400] / pm_shoot[800]))
assert PM_SHOOT_RN < 1e-8, PM_SHOOT_RN
for n_u, v in pm_shoot.items():
    print(f"  n_u = {n_u:4d}: worst relative deviation in eta over {n_states} branch states "
          f"{v:.3e}   worst Newton residual {pm_shoot_rn[n_u]:.1e}")
print(f"observed order 400 -> 800: {PM_SHOOT_ORDER:.2f}  (every state converged, so the")
print("order is a discretisation order and not a record of solves that stopped early)")
print(f"\\n{len(SKIPPED)} states were excluded for phi_0 > 30, and that exclusion is not")
print(f"cosmetic: they reach phi_0 = {max(s[2][0] for s in SKIPPED):.4g}, where the boundary layer is a")
print("fraction 1/phi_0 of the radius and no 800-cell grid resolves it. Saying so is")
print("part of the result: the criterion is never applied in that regime anyway,")
print("because Phi there is enormous and eq. (1) rejects it on sight.")'''))

cells.append(md(r"""### 3. What $\Phi < 1$ actually guarantees

The 1954 statement is an order of magnitude. With the forward problem in hand it
can be made exact: solve for the state at which $\Phi$ *equals* 1, and read off
$\eta$ there. That number is what an experimenter who just satisfies the
criterion is actually promised.

It is not one number. It depends on the shape and on the kinetics, and the
spread is larger than the criterion's own tolerance suggests."""))

cells.append(code('''def threshold_state(p, R, hi=25.0):
    """State where the observable Phi equals 1, by bisection on phi_0 in pymrm.

    Safe only where Phi(phi_0) is single valued, which is checked on the shooting
    branch before this is used. Every isothermal rate law on this page is: the
    kinetics cell below measures a fold depth of exactly zero in Phi for all five,
    including the inhibited Langmuir-Hinshelwood one. The non-isothermal cases of
    Section 5 are NOT, and no threshold there is located this way.

    Returns the Newton residual as its third value so the caller can assert that
    the state is a solution rather than infer it from the answer looking sensible.
    """
    f = lambda phi: p.Phi(p.solve(phi, R)[0]) - 1.0
    ph = brentq(f, 1e-3, hi, xtol=1e-12, rtol=1e-14)
    y, rn = p.solve(ph, R)
    return ph, p.eta(y, ph), rn


RN_TOL = 1e-8          # every Newton solve reported on this page must beat this


# --- geometry, first order: pymrm against the exact closed form ---------------
print("eta at Phi = 1, first-order isothermal kinetics")
print(f"{'geometry':>10} {'phi* pymrm':>11} {'eta* pymrm':>11} {'eta* exact':>11} {'rel':>9} "
      f"{'residual':>10}")
geo_rows, geo_dev, geo_rn = [], [], []
for geom, nu in (("slab", 0), ("cylinder", 1), ("sphere", 2)):
    ph, et, rn = threshold_state(Pellet(geom, 800), KIN["first order"])
    assert rn < RN_TOL, (geom, rn)
    phe = brentq(lambda f: f * f * eta_exact(nu, f) - 1.0, 1e-6, 25.0, xtol=1e-14, rtol=1e-15)
    ete = eta_exact(nu, phe)
    geo_rows.append((geom, ph, et, ete))
    geo_dev.append(abs(et - ete) / ete)
    geo_rn.append(rn)
    print(f"{geom:>10} {ph:11.6f} {et:11.6f} {ete:11.6f} {abs(et-ete)/ete:9.1e} {rn:10.1e}")
THRESH_EXACT_DEV = float(max(geo_dev))
THRESH_RN = float(max(geo_rn))
ETA_STAR_SPHERE = geo_rows[2][3]
ETA_STAR_SLAB = geo_rows[0][3]
print(f"\\npymrm reproduces the exact threshold to {THRESH_EXACT_DEV:.1e} relative, "
      f"worst Newton residual {THRESH_RN:.1e}.")
print(f"The guarantee spans {ETA_STAR_SLAB:.3f} (slab) to {ETA_STAR_SPHERE:.3f} (sphere):")
print(f"the SAME criterion admits a {100*(1-ETA_STAR_SLAB):.0f} % rate depression in one")
print(f"geometry and {100*(1-ETA_STAR_SPHERE):.1f} % in another, a factor "
      f"{(1-ETA_STAR_SLAB)/(1-ETA_STAR_SPHERE):.1f} in the error it lets through.")'''))

cells.append(code('''# --- the same, with the V/S length instead of R ------------------------------
# Phi is defined on R. On the volume-to-surface length L = R/(nu+1) it becomes
# Phi_VS = Phi/(nu+1)^2, which is the normalisation Aris' generalised modulus uses
# (see B1.2). Whether that makes the criterion's guarantee less shape-dependent is
# a question this page can answer rather than assert.
print("eta at Phi = 1 with each length scale")
print(f"{'geometry':>10} {'on R':>10} {'on V/S':>10}")
vs = []
for geom, nu in (("slab", 0), ("cylinder", 1), ("sphere", 2)):
    e_R = eta_exact(nu, brentq(lambda f: f * f * eta_exact(nu, f) - 1.0, 1e-6, 25.0,
                               xtol=1e-14, rtol=1e-15))
    e_V = eta_exact(nu, brentq(lambda f: f * f * eta_exact(nu, f) - (nu + 1.0) ** 2,
                               1e-6, 200.0, xtol=1e-14, rtol=1e-15))
    vs.append((e_R, e_V))
    print(f"{geom:>10} {e_R:10.4f} {e_V:10.4f}")
SPREAD_R = float(max(v[0] for v in vs) - min(v[0] for v in vs))
SPREAD_VS = float(max(v[1] for v in vs) - min(v[1] for v in vs))
print(f"\\nspread across geometries: {SPREAD_R:.3f} on R, {SPREAD_VS:.3f} on V/S "
      f"- a factor {SPREAD_R/SPREAD_VS:.1f} tighter")
print("So the shape sensitivity is mostly a choice of length scale, not physics.")
print("It is also why the criterion is quoted with different thresholds by different")
print("authors; the threshold and the length scale are not separable.")'''))

cells.append(code('''# --- kinetics, sphere ---------------------------------------------------------
# The branch is traced first, with the shooting reference, so that ALL states with
# Phi = 1 are found. A bisection on phi_0 would silently return one of them.
print("eta at Phi = 1, sphere, isothermal, all kinetics normalised to R(1) = 1")
print(f"{'kinetics':>20} {'roots':>6} {'phi*':>10} {'eta*':>10} {'pymrm eta*':>11} {'rel':>9} "
      f"{'residual':>10}")
KIN_ROOTS, kin_dev, kin_rn = {}, [], []
p800 = Pellet("sphere", 800)
for name, R in KIN.items():
    br = Branch(2, R, n=140)
    roots = br.all_where(lambda p, e, P: P - 1.0)
    KIN_ROOTS[name] = roots
    for j, (ts, ph, et, _) in enumerate(roots):
        # pymrm on the SAME branch state: seed Newton from the shooting profile
        _, y0 = br.seed(ts, p800)
        y, rn = p800.solve(ph, R, y_init=y0)
        assert rn < RN_TOL, (name, rn)
        dev = abs(p800.eta(y, ph) - et) / et
        kin_dev.append(dev)
        kin_rn.append(rn)
        print(f"{name if j == 0 else '':>20} {len(roots) if j == 0 else '':>6} "
              f"{ph:10.6f} {et:10.6f} {p800.eta(y, ph):11.6f} {dev:9.1e} {rn:10.1e}")
KIN_PYMRM_DEV = float(max(kin_dev))
KIN_RN = float(max(kin_rn))
first = {k: v[0][2] for k, v in KIN_ROOTS.items()}
KIN_MIN, KIN_MAX = float(min(first.values())), float(max(first.values()))
N_LH_ROOTS = len(KIN_ROOTS["inhibited LH (K=5)"])
print(f"\\npymrm agrees with the shooting branch on every one of these to "
      f"{KIN_PYMRM_DEV:.1e} relative, worst Newton residual {KIN_RN:.1e}.")
print(f"At the FIRST Phi = 1 crossing the guarantee runs {KIN_MIN:.4f} to {KIN_MAX:.4f}.")
print("Zero order gives exactly 1: eta = 1 identically until a dead core forms, which")
print("for a sphere is at Phi = 6, six times the threshold - the criterion is useless")
print("there, in the harmless direction.")
ETA_LH = float(KIN_ROOTS["inhibited LH (K=5)"][0][2])
print(f"\\nThe inhibited Langmuir-Hinshelwood rate gives eta = {ETA_LH:.4f} at Phi = 1, "
      f"ABOVE 1.")
print("Above y = 0.2 that rate falls with concentration, so depleting the pellet speeds")
print("it up. eta > 1 with no heat effect at all - and the criterion, which only sees")
print("Phi, cannot tell that apart from the exothermic enhancement of Section 4.")'''))

cells.append(code('''# --- why the thresholds are located on the branch and not by bisecting phi_0 ---
# This is not stylistic. A bisection on phi_0 that does not check the Newton
# residual returns a SECOND, spurious threshold for the inhibited rate law.
R_LH = KIN["inhibited LH (K=5)"]
br_lh = Branch(2, R_LH, n=200)
print(f"inhibited LH branch: fold depth in phi_0 {fold_depth(br_lh.phi):.2e}, "
      f"in Phi {fold_depth(br_lh.Phi):.2e}")
print("- single valued, so there is exactly one state with Phi = 1.\\n")
p = Pellet("sphere", 400)
print(f"{'phi_0':>8} {'Phi (pymrm)':>12} {'Newton residual':>16}")
for phi in (KIN_ROOTS["inhibited LH (K=5)"][0][1], 2.0, 5.5638):
    y, rn = p.solve(phi, R_LH)
    print(f"{phi:8.4f} {p.Phi(y):12.5f} {rn:16.1e}")
print("\\nThe third row also has Phi = 1 to four figures, and it is not a solution:")
print("its Newton residual is twelve orders of magnitude above the others. An earlier")
print("draft of this page bisected phi_0 for Phi = 1 and got that state, with a")
print("plausible-looking eta of 0.032. Locating the threshold on the shooting branch")
print("instead makes the failure impossible, and `threshold_state` returns the residual")
print("so it cannot be ignored. Assert that your solver converged; do not infer it")
print("from the answer looking reasonable.")'''))

cells.append(code('''# --- what the threshold costs in the quantity an experimenter reports ---------
# Weisz and Prater's printed relation, journal page 268:  Q'/Q = 1 + 0.5 dln eta/dln phi_0
def q_ratio(nu, phi, h=1e-5):
    return 1.0 + 0.5 * (np.log(eta_exact(nu, phi * (1 + h))) -
                        np.log(eta_exact(nu, phi * (1 - h)))) / np.log((1 + h) / (1 - h))


print("apparent activation energy at each geometry's own Phi = 1 threshold")
qr = {}
for geom, nu in (("slab", 0), ("cylinder", 1), ("sphere", 2)):
    phe = brentq(lambda f: f * f * eta_exact(nu, f) - 1.0, 1e-6, 25.0, xtol=1e-14, rtol=1e-15)
    qr[geom] = q_ratio(nu, phe)
    print(f"  {geom:>9}: Q'/Q = {qr[geom]:.4f}  ({100*(1-qr[geom]):.1f} % low)")
Q_SPHERE, Q_SLAB = float(qr["sphere"]), float(qr["slab"])

# the paper's own claim about the range of Q'/Q, checked
lim = np.array([q_ratio(2, f) for f in np.logspace(-3, 3, 400)])
Q_MIN, Q_MAX = float(lim.min()), float(lim.max())
Q_MONO = bool(np.all(np.diff(lim) <= 1e-4))     # tolerance = the truncation error below
print(f"\\nOver phi_0 = 1e-3 to 1e3, isothermal sphere: Q'/Q spans "
      f"[{Q_MIN:.4f}, {Q_MAX:.4f}], monotone decreasing: {Q_MONO}")
print(f"The paper states the measured value ranges 'between Q_0 and 1/2 Q_0'. It does.")
print(f"(The upper end overshoots 1 by {Q_MAX-1:.1e}; that is the central-difference")
print(" truncation error in the slope at phi_0 = 1e-3, not a property of eta.)")
print("This check can fail: a wrong factor in the printed relation, or a wrong")
print("eta(phi), would not land on 0.5 and 1.0 at the two ends. Substituting the")
print("factor 1 for 1/2 would give a lower limit of "
      f"{1 + 2*(q_ratio(2, 1e3) - 1):.4f} instead of 0.5.")'''))

cells.append(md(r"""### 4. The false-negative map: what the criterion would have concluded

Now the non-isothermal problem, which is where the isothermal criterion is
supposed to break and where Weisz and Hicks' eq. (10a) is supposed to rescue it.

The test sweeps the forward problem and asks what each criterion *would have
said*. For every $(\beta,\gamma)$, the branch is traced with the shooting
reference and three states are located exactly:

- the state where $\Phi = 1$, the isothermal criterion's own threshold — and
  $\eta$ there is what eq. (1) admits;
- the state where $\Phi = \exp[-\gamma\beta/(1+\beta)]$, eq. (10a)'s threshold —
  and $\eta$ there is what the extension admits;
- the state where $|\eta - 1|$ first reaches 5 %, which is the honest boundary
  of "no appreciable modification of chemical kinetics".

The ratio of the last to each threshold is the criterion's **safety factor**:
above 1 it is conservative, below 1 it is admitting a state it should not."""))

cells.append(code('''BETAS = (-0.4, -0.2, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6)
GAMMAS = (10.0, 20.0, 30.0, 40.0)
BR = {}                                   # cached branches, reused by every cell below
for gamma in GAMMAS:
    for beta in BETAS:
        BR[(beta, gamma)] = Branch(2, wh_rate(beta, gamma))
print(f"{len(BR)} branches traced, {sum(len(b.t) for b in BR.values())} shooting solves")'''))

cells.append(code('''MAP = {}
print("A '*' would mark a case where Phi = 1 has more than one solution, so that")
print("'eta there' is the first crossing only; Section 5 reports how many there are.")
print(f"\\n{'beta':>6} {'gamma':>6} {'b*g':>6} | {'eta @ Phi=1':>12} | {'Phi_c (10a)':>12} "
      f"{'eta there':>10} | {'Phi @ 5 %':>10} {'safety':>10}")
for gamma in GAMMAS:
    for beta in BETAS:
        br = BR[(beta, gamma)]
        r1 = br.all_where(lambda p, e, P: P - 1.0)
        s1 = r1[0][1:] if r1 else None
        Pc = float(np.exp(-gamma * beta / (1 + beta)))
        s2 = br.where(lambda p, e, P: P - Pc)
        s5 = br.where(lambda p, e, P: abs(e - 1.0) - 0.05)
        MAP[(beta, gamma)] = dict(eta_iso=None if s1 is None else s1[1], Phi_c=Pc,
                                  eta_nis=None if s2 is None else s2[1],
                                  Phi5=None if s5 is None else s5[2],
                                  n_roots=len(r1))
        f = None if (s5 is None) else s5[2] / Pc
        star = "*" if len(r1) > 1 else " "
        print(f"{beta:6.2f} {gamma:6.1f} {beta*gamma:6.1f} | "
              f"{s1[1]:11.5g}{star} | {Pc:12.4g} "
              f"{('%10.5g' % s2[1]) if s2 else '       n/a':>10} | "
              f"{('%10.4g' % s5[2]) if s5 else '       n/a':>10} "
              f"{('%10.3g' % f) if f else '       n/a':>10}")'''))

cells.append(code('''exo = [(k, v) for k, v in MAP.items() if k[0] > 0]
endo = [(k, v) for k, v in MAP.items() if k[0] < 0]
WORST_ISO_EXO = max(v["eta_iso"] for _, v in exo)
WORST_ISO_EXO_AT = max(exo, key=lambda kv: kv[1]["eta_iso"])[0]
WORST_NIS_EXO = max(abs(v["eta_nis"] - 1.0) for _, v in exo if v["eta_nis"])
WORST_ISO_ENDO = min(v["eta_iso"] for _, v in endo)
WORST_ISO_ENDO_AT = min(endo, key=lambda kv: kv[1]["eta_iso"])[0]
ETA_ISO_B0 = MAP[(0.0, 20.0)]["eta_iso"]
MOST_CONSERVATIVE = max((v["Phi5"] / v["Phi_c"], k) for k, v in exo if v["Phi5"])

print(f"eq. (1), Phi < 1, applied to an EXOTHERMIC pellet:")
print(f"  worst eta it calls safe : {WORST_ISO_EXO:.4g} at beta = {WORST_ISO_EXO_AT[0]}, "
      f"gamma = {WORST_ISO_EXO_AT[1]:.0f}")
print(f"  i.e. the observed rate is {WORST_ISO_EXO:.0f} times the intrinsic one and the")
print(f"  criterion reports no transport effect. This is the failure eq. (10a) exists for.")
print(f"\\neq. (10a), the same states:")
print(f"  worst |eta - 1| it calls safe : {WORST_NIS_EXO:.4f}")
print(f"  So over the whole exothermic sweep the extension holds eta to within "
      f"{100*WORST_NIS_EXO:.1f} % of 1.")
print(f"  It pays for that with conservatism: at beta = {MOST_CONSERVATIVE[1][0]}, "
      f"gamma = {MOST_CONSERVATIVE[1][1]:.0f} its threshold is a factor")
print(f"  {MOST_CONSERVATIVE[0]:.3g} below the true 5 % point, so it rejects conditions "
      f"that are in fact clean.")
print(f"\\neq. (1) at beta = 0, any gamma      : eta = {ETA_ISO_B0:.6f}  "
      f"({100*(1-ETA_ISO_B0):.1f} % depression admitted)")
print(f"eq. (1) applied to an ENDOTHERMIC pellet:")
print(f"  worst eta it calls safe : {WORST_ISO_ENDO:.4f} at beta = {WORST_ISO_ENDO_AT[0]}, "
      f"gamma = {WORST_ISO_ENDO_AT[1]:.0f}")
n_endo_na = sum(1 for _, v in endo if v["eta_nis"] is None)
print(f"eq. (10a) applied to an ENDOTHERMIC pellet: its threshold is "
      f"exp[-gamma*beta/(1+beta)] > 1,")
print(f"  reaching {max(v['Phi_c'] for _, v in endo):.3g}; the branch never gets that far in "
      f"{n_endo_na} of {len(endo)} cases,")
print("  so the criterion declares every attainable state safe. The paper says (10) is")
print("  'a sufficient criterion for exothermic reactions' and this is why the")
print("  restriction is load-bearing, not stylistic.")'''))

cells.append(code('''# --- a confusion count, over a stated grid of MEASUREMENTS --------------------
# The grid is in Phi, not in branch states, because Phi is what is measured. For
# each (beta, gamma) and each target Phi, every steady state consistent with that
# measurement is collected by linear interpolation along the traced branch - so
# where Phi is multivalued the experimenter's single number carries all of them,
# which is exactly the situation the criterion has to survive.
TOL, PHI_GRID = 0.05, np.logspace(-3, 2, 81)


def etas_at(br, target):
    """Every eta on the branch consistent with an observed Phi = target."""
    out, lp, lt = [], np.log(br.Phi), np.log(target)
    for i in range(len(lp) - 1):
        if (lp[i] - lt) * (lp[i + 1] - lt) < 0:
            w = (lt - lp[i]) / (lp[i + 1] - lp[i])
            out.append(br.eta[i] + w * (br.eta[i + 1] - br.eta[i]))
    return out


rows = []
for (beta, gamma), br in BR.items():
    for P in PHI_GRID:
        es = etas_at(br, P)
        if not es:
            continue
        # the measurement is "clean" only if EVERY state it is consistent with is
        truth = all(abs(e - 1.0) <= TOL for e in es)
        rows.append((beta, gamma, truth,
                     P < 1.0,
                     P * np.exp(gamma * beta / (1 + beta)) < 1.0,
                     len(es)))


def confuse(sel, idx):
    r = [v for v in rows if sel(v)]
    fn = sum(1 for v in r if v[idx] and not v[2])
    fp = sum(1 for v in r if (not v[idx]) and v[2])
    return fn, fp, len(r)


print(f"verdicts on {len(rows)} synthetic measurements ({len(PHI_GRID)} values of Phi "
      f"from 1e-3 to 1e2, x {len(BR)} (beta, gamma)),")
print(f"truth = every state consistent with that Phi has |eta - 1| <= {TOL}")
print(f"\\n{'subset':>14} {'criterion':>10} {'false neg':>10} {'false pos':>10} "
      f"{'measurements':>13}")
for label, sel in (("exothermic", lambda v: v[0] > 0), ("isothermal", lambda v: v[0] == 0),
                   ("endothermic", lambda v: v[0] < 0)):
    for name, idx in (("eq. (1)", 3), ("eq. (10a)", 4)):
        fn, fp, n = confuse(sel, idx)
        print(f"{label:>14} {name:>10} {fn:9d}  {fp:9d}  {n:13d}")
FN_ISO_EXO, FP_ISO_EXO, _ = confuse(lambda v: v[0] > 0, 3)
FN_NIS_EXO, FP_NIS_EXO, N_EXO = confuse(lambda v: v[0] > 0, 4)
FN_ISO_ENDO = confuse(lambda v: v[0] < 0, 3)[0]
FN_NIS_ENDO, _, N_ENDO = confuse(lambda v: v[0] < 0, 4)
FN_ISO_ISO, _, N_ISO = confuse(lambda v: v[0] == 0, 3)
N_AMBIG = sum(1 for v in rows if v[5] > 1)
print("\\nA false negative is the dangerous one: the criterion says the measurement is")
print("clean kinetics and it is not. A false positive only wastes an experiment.")
print(f"exothermic : eq. (1) gives {FN_ISO_EXO} false negatives out of {N_EXO}; "
      f"eq. (10a) gives {FN_NIS_EXO},")
print(f"             at the price of {FP_NIS_EXO} false positives against eq. (1)'s "
      f"{FP_ISO_EXO}.")
print(f"endothermic: eq. (1) gives {FN_ISO_ENDO} out of {N_ENDO}; eq. (10a) gives "
      f"{FN_NIS_ENDO} - it is WORSE than")
print("             doing nothing, which is why its stated restriction to exothermic")
print("             reactions matters.")
print(f"isothermal : eq. (1) gives {FN_ISO_ISO} out of {N_ISO} - the band "
      f"0.767 < Phi < 1 where")
print(f"             eta has already fallen below {1-TOL:.2f} but the criterion still passes.")
print(f"\\n{N_AMBIG} of the {len(rows)} measurements are consistent with more than one "
      f"steady state.")
print("Counts depend on the grid; the continuous quantities above do not, which is why")
print("they and not these are the headline.")'''))

cells.append(md(r"""### 5. Does the observable remove the multiplicity?

The paper's Section V says it does, and the criterion needs it to. The test is
direct: trace the branch, and ask whether $\Phi$ is a **monotone** function of
the state. If it is, a measured $\Phi$ picks out one $\eta$ and the criterion is
well posed. If it is not, the same observed rate belongs to more than one steady
state and no threshold on $\Phi$ can separate them.

`fold_depth` returns the largest relative backward excursion along the traced
branch: exactly zero means strictly monotone.

Whether a fold *matters* depends on where it sits relative to $\Phi = 1$, so the
turning points are located rather than described: a fold above the threshold
would only blur states eq. (1) has already rejected, while a fold below it is an
ambiguity inside the band eq. (1) certifies as safe. Every one found here is of
the second kind.

The same measurement applied to $\phi_0$ reproduces the paper's other printed
statement, $\beta\gamma \gtrsim 5$, and gives it a number."""))

cells.append(code('''print(f"{'beta':>6} {'gamma':>6} {'b*g':>6} {'fold in phi_0':>14} {'fold in Phi':>12} "
      f"{'eta max':>10}")
FOLD = {}
for gamma in GAMMAS:
    for beta in BETAS:
        br = BR[(beta, gamma)]
        fp, fP = fold_depth(br.phi), fold_depth(br.Phi)
        FOLD[(beta, gamma)] = (fp, fP)
        flag = "   <-- Phi is multivalued too" if fP > 1e-4 else ""
        print(f"{beta:6.2f} {gamma:6.1f} {beta*gamma:6.1f} {fp:14.3e} {fP:12.3e} "
              f"{br.eta.max():10.4g}{flag}")

bad = {k: v for k, v in FOLD.items() if v[1] > 1e-4}
WORST_PHI_FOLD = float(max(v[1] for v in FOLD.values()))
WORST_PHI_FOLD_AT = max(FOLD.items(), key=lambda kv: kv[1][1])[0]
N_PHI_FOLD = sum(1 for v in FOLD.values() if v[0] > 1e-4)
N_OBS_FOLD = len(bad)
print(f"\\nphi_0 is multivalued in {N_PHI_FOLD} of {len(FOLD)} cases; "
      f"Phi in {N_OBS_FOLD} of {len(FOLD)}.")
print(f"So the transformation removes the multiplicity in "
      f"{N_PHI_FOLD - N_OBS_FOLD} of the {N_PHI_FOLD} cases that have it - but not in all.")
print(f"Worst surviving fold: {WORST_PHI_FOLD:.3f} at beta = {WORST_PHI_FOLD_AT[0]}, "
      f"gamma = {WORST_PHI_FOLD_AT[1]:.0f}.")
print("Every one of those cases lies inside the range the paper states it computed,")
print("'gamma = 10, 20, 30 and 40, and beta in the range from 0 to +0.8'.")
FOLDED = [k for k, v in FOLD.items() if v[1] > 1e-4]
print("Specifically: " + ", ".join("(beta %g, gamma %g)" % k for k in sorted(FOLDED)) + ".")'''))

cells.append(code('''# --- WHERE the folds sit, on the very axis eq. (1) puts its threshold on -------
# This is the question that decides whether the fold matters. A fold ABOVE Phi = 1
# would only blur states the criterion has already rejected. A fold BELOW it is a
# failure inside the band the criterion certifies as safe. So the turning points are
# located, not described.
TURN = {k: fold_turning_points(BR[k]) for k in sorted(FOLDED)}
print(f"{'beta':>6} {'gamma':>6} | {'turning points of Phi':<44} | eta at each")
for k, tp in TURN.items():
    print(f"{k[0]:6.2f} {k[1]:6.1f} | "
          + f"{', '.join('%s %.4f' % (kind, P) for kind, _, _, _, P in tp):<44} | "
          + ", ".join(f"{e:.4g}" for _, _, _, e, _ in tp))
TURN_PHI = [P for tp in TURN.values() for *_, P in tp]
print(f"\\nAll {len(TURN_PHI)} turning points lie BELOW Phi = 1; the largest is "
      f"{max(TURN_PHI):.4f}.")

# every measurement on the Section 4 grid that more than one steady state fits
AMB = []
for (beta, gamma), br in BR.items():
    for P in PHI_GRID:
        es = etas_at(br, P)
        if len(es) > 1:
            AMB.append((beta, gamma, P, es))
assert len(AMB) == N_AMBIG, (len(AMB), N_AMBIG)
AMB_PHI_MAX = float(max(a[2] for a in AMB))
AMB_WORST = max(AMB, key=lambda a: max(a[3]) / min(a[3]))
AMB_DEEPEST = max(AMB, key=lambda a: a[2])
n_star = sum(1 for v in MAP.values() if v["n_roots"] > 1)

print(f"\\nThe equation Phi = 1 itself has a unique solution in {len(MAP) - n_star} of the "
      f"{len(MAP)} cases, so eq. (1)")
print("keeps a single, well defined verdict boundary everywhere here. That is not the")
print("reassurance it looks like, because of WHERE the folds are. They are not out")
print("beyond the threshold, among states the criterion already rejects. Every one of")
print("them is inside the band it certifies as safe - and so is every ambiguous")
print(f"measurement on the grid above: all {len(AMB)} of the {len(rows)} have Phi < 1, the "
      f"largest {AMB_PHI_MAX:.4f}.")
for tag, (b, g, P, es) in (("worst spread", AMB_WORST),
                           ("deepest into the band", AMB_DEEPEST)):
    print(f"\\n  {tag} - beta = {b}, gamma = {g:.0f}: an observed Phi = {P:.4f} passes eq. (1)")
    print(f"  by a factor of {1/P:.1f} and is consistent with eta = "
          + ", ".join(f"{e:.4g}" for e in es) + ".")
print("\\nSo what eq. (1) loses where Phi folds is not its verdict - it is the meaning of")
print("a PASS. That is a SECOND false-negative mechanism, independent of the loose")
print("threshold Section 4 measures: not a criterion that admits too much, but a")
print("measurement well inside the safe band that does not determine the state at all.")
print("An experimenter who clears the criterion threefold can still be sitting on an")
print("ignited pellet running two orders of magnitude above its intrinsic rate, and no")
print("amount of tightening the threshold would help, because the ambiguity is not at")
print("the boundary.")'''))

cells.append(code('''# --- where does each fold set in? --------------------------------------------
# A cheaper branch (n = 80) is used here because this is a bisection over branches;
# the values are cross-checked against the n = 110 table above at the end.
def depth(beta, gamma, which, n=80):
    br = Branch(2, wh_rate(beta, gamma), n=n)
    return fold_depth(br.phi if which == "phi" else br.Phi)


ONSET = {}
print(f"{'gamma':>6} {'phi_0 folds at':>16} {'b*g':>7} | {'Phi folds at':>14} {'b*g':>7}")
for gamma in (20.0, 30.0, 40.0):
    out = []
    for which in ("phi", "Phi"):
        f = lambda b: depth(b, gamma, which) - 1e-4
        out.append(brentq(f, 0.03, 1.3, xtol=4e-3) if f(1.3) > 0 else None)
    ONSET[gamma] = out
    a, b = out
    print(f"{gamma:6.0f} {a:16.3f} {a*gamma:7.2f} | "
          f"{('%14.3f' % b) if b else '           n/a'} "
          f"{('%7.2f' % (b*gamma)) if b else '    n/a'}")

PHI_ONSET = [ONSET[g][0] * g for g in (20.0, 30.0, 40.0)]
OBS_ONSET = [ONSET[g][1] * g for g in (20.0, 30.0, 40.0) if ONSET[g][1]]
print(f"\\nphi_0 multiplicity sets in at beta*gamma = {min(PHI_ONSET):.2f} to "
      f"{max(PHI_ONSET):.2f}.")
print("The paper states 'beta*gamma >~ 5'. That is reproduced, and the check can fail:")
print("a wrong geometry, a wrong exponent or a bad reference solver would not land")
print("there. The drift with gamma is real - the group beta*gamma is not by itself")
print("the controlling parameter, only nearly so - which is what '>~' is doing.")
print(f"\\nPhi multiplicity sets in at beta*gamma = {min(OBS_ONSET):.1f} to "
      f"{max(OBS_ONSET):.1f}, a factor "
      f"{min(OBS_ONSET)/max(PHI_ONSET):.1f} to {max(OBS_ONSET)/min(PHI_ONSET):.1f} higher.")
print("THAT is the quantitative content of Section V's claim. The observable does not")
print("remove the multiplicity; it postpones it by a factor of about two to four in")
print("beta*gamma. That covers every case the paper argues about in the text, and it")
print("does not cover the whole of the beta and gamma range it says it computed.")
# The largest beta labelled on each panel of Fig. 7, read off the printed curve
# labels on a 600 dpi render of journal page 272. Not digitised - these are the
# text labels beside the curves, and they bound which curves the paper drew.
FIG7_BETA_MAX = {10.0: 0.8, 20.0: 0.8, 30.0: 0.4, 40.0: 0.3}
print("\\nWhich curves the paper actually DREW are affected:")
n_drawn = 0
for g in GAMMAS:
    bmax = FIG7_BETA_MAX[g]
    if g not in ONSET:
        # gamma = 10 is outside the onset bisection above, which costs a branch per
        # iteration. The panel count needs a narrower question anyway - does the
        # LARGEST curve this panel draws fold? - so it is answered at that beta
        # directly, rather than left uncounted and then counted in the total.
        d = fold_depth(Branch(2, wh_rate(bmax, g)).Phi)
        hit = d > 1e-4
        n_drawn += int(hit)
        print(f"  Fig. 7 at gamma = {g:.0f}: fold in Phi at the largest beta drawn "
              f"({bmax}) is {d:.1e} -> {'AFFECTED' if hit else 'not reached'}")
        continue
    b = ONSET[g][1]
    if b is None:
        print(f"  Fig. 7 at gamma = {g:.0f}: no fold found for beta <= 1.3; "
              f"largest beta drawn is {bmax}")
        continue
    hit = b <= bmax
    n_drawn += int(hit)
    print(f"  Fig. 7 at gamma = {g:.0f}: Phi folds for beta > {b:.3f}; largest beta "
          f"drawn is {bmax} -> {'AFFECTED' if hit else 'not reached'}")
print(f"So {n_drawn} of the four panels contains a drawn curve in the folded regime.")
print(f"The large folds ({FOLD[(0.6, 40.0)][1]:.2f} at beta = 0.6, gamma = 40) are inside "
      f"the beta range the")
print("text says was computed but outside the beta range Fig. 7 draws. The claim being")
print("tested is uniqueness, not the stability hedge the paper attaches to it; a folded")
print("Phi means the observable does not determine the state, whether or not every")
print("state on the fold is stable.")

# --- could Fig. 7(d) have shown the one fold it draws? MEASURED, in pixels -----
# The resolving-power rule: compute the effect in pixels before quoting the figure
# as agreeing or disagreeing. Two quantities are measured off the same 600 dpi
# render of journal page 272 that the beta labels came from, and nothing else is
# taken from the image - no curve is digitised, no value is read off an axis:
#   * the Phi axis of panel (d). Its five major gridlines, Phi = 0.1, 1, 10, 100,
#     1000, have intensity-weighted centroids at x = 199.6, 560.1, 917.5, 1275.7
#     and 1631.3 px, so one decade is (1631.3 - 199.6)/4 px.
#   * the printed stroke. The horizontal width of the two near-vertical dashed
#     curve segments in that panel, over 163 scan rows clear of both gridlines,
#     has median 6 px at an ink threshold of 160/255 and 7 px at 200/255. (The
#     gridlines themselves are drawn much finer, 2 px, so this is the curve pen.)
FIG7D_DECADE_PX = (1631.3 - 199.6) / 4.0
FIG7D_STROKE_PX = (6.0, 7.0)                 # ink thresholds 160 and 200 of 255
tp_d = TURN[(0.3, 40.0)]
P_HI, P_LO = tp_d[0][4], tp_d[1][4]
FOLD_DECADES = float(np.log10(P_HI / P_LO))
FOLD_PX = float(FOLD_DECADES * FIG7D_DECADE_PX)
print(f"\\nCan Fig. 7(d) resolve the fold on the one drawn curve that has one?")
print(f"  beta = 0.3, gamma = 40: Phi turns back at {P_HI:.4f} and forward again at "
      f"{P_LO:.4f},")
print(f"  a width of log10({P_HI:.4f}/{P_LO:.4f}) = {FOLD_DECADES:.4f} decades.")
print(f"  Panel (d) prints {FIG7D_DECADE_PX:.1f} px per decade at 600 dpi, so the fold is")
print(f"  {FOLD_DECADES:.4f} x {FIG7D_DECADE_PX:.1f} = {FOLD_PX:.1f} px wide.")
print(f"  The curve is drawn {FIG7D_STROKE_PX[0]:.0f} to {FIG7D_STROKE_PX[1]:.0f} px wide, "
      f"so the fold spans "
      f"{FOLD_PX/FIG7D_STROKE_PX[1]:.1f} to {FOLD_PX/FIG7D_STROKE_PX[0]:.1f} line widths.")
print("  It is COMPARABLE to the stroke, not below it. The honest conclusion is that")
print("  Fig. 7 is silent: at that width the printed curve can neither show the fold")
print("  nor rule it out, so the figure is not evidence in either direction and the")
print("  claim here rests on the computation alone. (An earlier draft of this page")
print("  asserted the fold was 'below the line width' without measuring; it is not.)")
print(f"  Note the refined turning points give a slightly deeper fold, "
      f"{(P_HI - P_LO)/P_HI:.4f},")
print(f"  than the {FOLD[(0.3, 40.0)][1]:.4f} of the traced table above - that is the "
      f"branch spacing, and")
print("  the deeper value is the one used here, which is the conservative direction")
print("  for a claim that the figure cannot resolve it.")
print("\\nThe onset threshold (fold depth > 1e-4) is a convention, so its sensitivity")
print("is measured rather than assumed:")
for thr in (1e-5, 1e-4, 1e-3):
    b = brentq(lambda x: depth(x, 30.0, "Phi") - thr, 0.03, 1.3, xtol=4e-3)
    print(f"  gamma = 30, threshold {thr:.0e}: beta* = {b:.3f} (beta*gamma = {b*30:.2f})")
print("Branch resolution matters more than the threshold does: a fold shallower than")
print("the spacing between traced states cannot be seen at all, which biases every")
print("onset here slightly LATE. The reported onsets are therefore upper bounds.")'''))

# ------------------------------------------------------------------- figures
cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
cols = plt.cm.viridis(np.linspace(0.05, 0.85, 4))
for c, beta in zip(cols, (0.0, 0.2, 0.4, 0.6)):
    br = BR[(beta, 20.0)]
    axes[0].loglog(br.phi, br.eta, color=c, lw=1.6, label=rf"$\beta$ = {beta}")
    axes[1].loglog(br.Phi, br.eta, color=c, lw=1.6)
axes[0].set_xlabel(r"Thiele modulus $\phi_0$  (needs $k$)")
axes[1].set_xlabel(r"observable $\Phi = \phi_0^2\eta$  (needs no $k$)")
for a in axes:
    a.set_ylabel(r"$\eta$"); a.set_ylim(3e-2, 3e2); a.set_xlim(1e-2, 1e2)
    a.axhline(1.0, color="0.5", lw=0.8, ls=":")
axes[1].axvline(1.0, color="crimson", lw=1.0, ls="--")
axes[1].text(1.15, 4e-2, "eq. (1)", color="crimson", fontsize=9)
axes[0].legend(fontsize=8, loc="lower left")
axes[0].set_title(r"$\gamma = 20$: $\eta$ folds back in $\phi_0$")
axes[1].set_title(r"the same states, plotted against the observable")
fig.suptitle("The transformation Weisz and Hicks' Section V is about", y=1.01)
fig.tight_layout(); plt.show()'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
for c, beta in zip(plt.cm.plasma(np.linspace(0.05, 0.7, 3)), (0.2, 0.4, 0.6)):
    br = BR[(beta, 40.0)]
    axes[0].loglog(br.phi, br.eta, color=c, lw=1.6, label=rf"$\beta$ = {beta}")
    axes[1].loglog(br.Phi, br.eta, color=c, lw=1.6,
                   label=rf"$\beta$ = {beta}, fold {FOLD[(beta, 40.0)][1]:.2f}")
axes[0].set_xlabel(r"$\phi_0$"); axes[1].set_xlabel(r"$\Phi$")
for a in axes:
    a.set_ylabel(r"$\eta$"); a.axhline(1.0, color="0.5", lw=0.8, ls=":")
    a.legend(fontsize=8, loc="upper left"); a.set_ylim(0.5, 1e5)
axes[0].set_xlim(1e-3, 3.0)
axes[1].set_xlim(1e-3, 1e2)
axes[1].axvline(1.0, color="crimson", lw=1.0, ls="--")
axes[1].text(1.15, 0.7, "eq. (1)", color="crimson", fontsize=9)
axes[0].set_title(r"$\gamma = 40$ — inside the paper's Fig. 4 range")
axes[1].set_title(r"$\Phi$ folds back too: one $\Phi$, three $\eta$")
fig.tight_layout(); plt.show()'''))

cells.append(code(r'''fig, ax = plt.subplots(figsize=(7.2, 4.4))
for c, gamma in zip(plt.cm.viridis(np.linspace(0.05, 0.85, 4)), GAMMAS):
    b = [k[0] for k in MAP if k[1] == gamma]
    e1 = [MAP[(x, gamma)]["eta_iso"] for x in b]
    e2 = [MAP[(x, gamma)]["eta_nis"] for x in b]
    ax.semilogy(b, e1, "o-", color=c, lw=1.5, ms=4, label=rf"eq. (1), $\gamma$ = {gamma:.0f}")
    ax.semilogy([x for x, v in zip(b, e2) if v], [v for v in e2 if v], "s--",
                color=c, lw=1.0, ms=3, alpha=0.7)
ax.axhspan(0.95, 1.05, color="0.85", zorder=0)
ax.text(0.62, 1.0, r"$|\eta-1|\leq 5\,\%$", fontsize=8, va="center", ha="left")
ax.set_xlabel(r"Prater number $\beta$"); ax.set_ylabel(r"$\eta$ at the criterion's threshold")
ax.set_title("What each criterion admits (solid: eq. 1;  dashed: eq. 10a)")
ax.legend(fontsize=8, loc="upper left")
ax.set_xlim(-0.48, 0.85)
ax.annotate("eq. (10a) is stated for EXOTHERMIC reactions only;\n"
            "its endothermic branch is drawn to show why",
            xy=(-0.4, 4e-3), xytext=(-0.13, 3e-3), fontsize=7.5,
            arrowprops=dict(arrowstyle="->", lw=0.7))
fig.tight_layout(); plt.show()'''))

# --------------------------------------------------- results: defect injection
cells.append(md(r"""### 6. Breaking it on purpose

Every number above is a claim, and a claim is only worth what its check can
detect. This section injects defects and measures what moves.

It also settles the status of one check that looks like evidence and is not. The
observable can be computed two ways from the same converged solution: as the
**surface flux** $(\nu+1)\,\mathrm{d}y/\mathrm{d}u|_1$, or as the **volume
average of the intrinsic rate** $\phi_0^2\langle\mathcal{R}\rangle$. Eq. (11)
says they are the same thing, and they agree to roundoff.

That agreement is worth nothing on its own, and the reason is provable rather
than empirical. Multiply the discrete residual in cell $i$ by that cell's volume
and sum: `construct_div` is conservative, so the divergence term telescopes to
the surface flux and what is left is exactly

$$\Phi_{\text{flux}} - \Phi_{\text{volume}}
= \frac{\nu+1}{V}\sum_i V_i\,\mathrm{res}_i ,$$

**the volume-weighted mean Newton residual wearing a different hat.** So it can
detect a solve that has not converged and a mismatch between the rate in the
model and the rate used to interpret it, and it is structurally incapable of
detecting a wrong geometry index or a bad grid — a wrong $\nu$ changes both
sides identically. The table below measures each of those rather than asserting
them, and the conclusion is that this check adds nothing the Newton residual
does not already give."""))

cells.append(code('''R1 = KIN["first order"]
ref = Pellet("sphere", 400)
PH_REF, ET_REF, RN_REF = threshold_state(ref, R1)
y_ref, _ = ref.solve(PH_REF, R1)
FLUX_VOL = abs(ref.Phi(y_ref) - ref.Phi_from_volume(y_ref, PH_REF, R1))
print(f"reference: phi* = {PH_REF:.6f}, eta* = {ET_REF:.6f}, Newton residual {RN_REF:.1e}")
print(f"           Phi(flux) - Phi(volume) = {FLUX_VOL:.2e}\\n")

DEF = {}
# Each defect is a DIFFERENT MODEL solved correctly, not a failed solve, so the
# residual is asserted rather than discarded: without that, "the defect moved eta"
# and "the defect stopped Newton converging" are indistinguishable.
# (i) wrong geometry index
p = Pellet("sphere", 400, nu_override=0)
ph, et, rn = threshold_state(p, R1)
assert rn < RN_TOL, rn
y, _ = p.solve(ph, R1)
DEF["nu = 0 for a sphere"] = (et, abs(p.Phi(y) - p.Phi_from_volume(y, ph, R1)), rn)
# (ii) three-cell grid
p = Pellet("sphere", 3)
ph, et, rn = threshold_state(p, R1)
assert rn < RN_TOL, rn
y, _ = p.solve(ph, R1)
DEF["n_u = 3"] = (et, abs(p.Phi(y) - p.Phi_from_volume(y, ph, R1)), rn)
# (iii) rate 1 % too large in the model, interpreted with the true rate
p = Pellet("sphere", 400)
Rs = lambda y: 1.01 * y
ph, et, rn = threshold_state(p, Rs)
assert rn < RN_TOL, rn
y, _ = p.solve(ph, Rs)
DEF["rate scale 1 % off"] = (et, abs(p.Phi(y) - p.Phi_from_volume(y, ph, R1)), rn)
print(f"{'injected defect':>22} {'eta at Phi=1':>13} {'moved?':>8} {'flux-vol':>11} "
      f"{'residual':>10}")
print(f"{'(none)':>22} {ET_REF:13.6f} {'-':>8} {FLUX_VOL:11.2e} {RN_REF:10.1e}")
for k, (et, fv, rn) in DEF.items():
    mv = "yes" if abs(et - ET_REF) / ET_REF > 1e-3 else "NO"
    print(f"{k:>22} {et:13.6f} {mv:>8} {fv:11.2e} {rn:10.1e}")
print("(all three converged, so every move in eta below is the defect and not a solve")
print(" that gave up)")

# (iv) Newton stopped early on a genuinely hard, genuinely converged reference:
# an ignited state seeded from the shooting profile, so the "converged" comparison
# is known to be a real solution and not another failed solve.
br_h = BR[(0.6, 20.0)]
t_h = br_h.t[np.argmax(br_h.eta)]
hard = Pellet("sphere", 400)
PHI_H, y_seed = br_h.seed(t_h, hard)
R_H = wh_rate(0.6, 20.0)
y_ok, rn_ok = hard.solve(PHI_H, R_H, y_init=y_seed)
ETA_H_REF = shoot(2, R_H, t_h, rtol=1e-11)[1]
print(f"\\nignited reference state (beta = 0.6, gamma = 20, phi_0 = {PHI_H:.5f}):")
print(f"  converged  : eta = {hard.eta(y_ok, PHI_H):11.5f}  residual {rn_ok:.2e}  "
      f"vs shooting {ETA_H_REF:.5f}")
for mf in (1, 3):
    y_bad, rn_bad = hard.solve(PHI_H, R_H, y_init=np.ones(hard.shape), maxfev=mf)
    fv = abs(hard.Phi(y_bad) - hard.Phi_from_volume(y_bad, PHI_H, R_H))
    print(f"  maxfev = {mf}: eta = {hard.eta(y_bad, PHI_H):11.5f}  residual {rn_bad:.2e}  "
          f"flux-vol {fv:.2e}")
    if mf == 1:
        BLIND_NOCONV = float(fv)
        ETA_NOCONV = float(hard.eta(y_bad, PHI_H))
        RN_NOCONV = float(rn_bad)
    else:
        FV_BRANCH = float(fv)
        ETA_BRANCH = float(hard.eta(y_bad, PHI_H))
        RN_BRANCH = float(rn_bad)
ETA_H_OK = float(hard.eta(y_ok, PHI_H))
print("\\nRow 2 is what it looks like: Newton stopped after one step and neither")
print("diagnostic is fooled. The flux-volume identity DOES see a non-converged solve -")
print("it is the residual - but it tells you nothing the Newton residual beside it did")
print("not already say.")

# Row 3 is a different, and worse, failure - so it gets said rather than left in the
# table. Count the states the branch carries at this phi_0, interpolating on the
# traced branch (the exact branch point is hit, hence the a == 0 case).
n_states_here = 0
for i in range(len(br_h.phi) - 1):
    a_, b_ = br_h.phi[i] - PHI_H, br_h.phi[i + 1] - PHI_H
    n_states_here += int(a_ == 0.0 or a_ * b_ < 0.0)
print(f"\\nRow 3 is NOT a non-convergence, and it is the more instructive row. Its Newton")
print(f"residual is {RN_BRANCH:.2e} and its flux-volume gap {FV_BRANCH:.1e}: both diagnostics "
      f"are clean,")
print("because it is a genuine solution of the discrete system - just not the one that")
print(f"was asked for. The fold in phi_0 here is {FOLD[(0.6, 20.0)][0]:.2f} deep, so the "
      f"branch carries")
print(f"{n_states_here} steady states at phi_0 = {PHI_H:.5f}, and three Newton steps from the "
      f"y = 1 start")
print(f"landed on the lowest: eta = {ETA_BRANCH:.5f} against {ETA_H_OK:.5f} for the seeded "
      f"solve, a factor")
print(f"{ETA_H_OK/ETA_BRANCH:.0f} out. No residual can catch that, because nothing is wrong "
      f"with the residual.")
print("A converged solve tells you that you solved the equations; it does not tell you")
print("WHICH solution you found. Only the branch trace does, which is why every")
print("multiplicity result on this page is located on the shooting branch first - and")
print("it is the same failure, in miniature, as the one Section 5 reports for an")
print("experimenter reading a single number off a folded curve.")
DEF_SCALE = float(DEF["rate scale 1 % off"][1])
BLIND_GEOM = float(DEF["nu = 0 for a sphere"][1])
BLIND_COARSE = float(DEF["n_u = 3"][1])'''))

cells.append(code('''# --- do the checks that are supposed to catch these actually catch them? ------
print("the same defects, measured against the checks this page reports")
print(f"{'defect':>22} {'vs exact eta':>13} {'threshold eta*':>15}")
base = max(abs(Pellet("sphere", 400).eta(
    Pellet("sphere", 400).solve(f, R1)[0], f) - eta_exact(2, f)) / eta_exact(2, f)
    for f in PHIS)
print(f"{'(none)':>22} {base:13.2e} {ET_REF:15.6f}")
DEF_VS_EXACT = {}
for label, kw in (("nu = 0 for a sphere", dict(nu_override=0)), ("n_u = 3", dict(n_u=3))):
    q = Pellet("sphere", **({"n_u": 400} | kw))
    d = max(abs(q.eta(q.solve(f, R1)[0], f) - eta_exact(2, f)) / eta_exact(2, f) for f in PHIS)
    DEF_VS_EXACT[label] = float(d)
    print(f"{label:>22} {d:13.2e} {DEF[label][0]:15.6f}")
DEF_COARSE_VS_EXACT = DEF_VS_EXACT["n_u = 3"]

# --- the map's own sensitivity ------------------------------------------------
b_s = Branch(2, wh_rate(-0.3, 20.0)).where(lambda p, e, P: P - 1.0)[1]
b_r = MAP[(0.3, 20.0)]["eta_iso"]
Pc_ok = float(np.exp(-20.0 * 0.3 / 1.3))
Pc_bad = float(np.exp(-20.0 * 0.3))
s_bad = BR[(0.3, 20.0)].where(lambda p, e, P: P - Pc_bad)
print(f"\\nsign of beta flipped in the exponent (gamma = 20, beta = 0.3):")
print(f"  eta at Phi = 1 : {b_r:.4f} -> {b_s:.4f}")
print(f"eq. (10a) exponent mis-transcribed as gamma*beta instead of gamma*beta/(1+beta):")
print(f"  threshold Phi_c: {Pc_ok:.4e} -> {Pc_bad:.4e}, and eta there "
      f"{MAP[(0.3,20.0)]['eta_nis']:.5f} -> {s_bad[1]:.5f}")
print("reference solver degraded to rtol = 1e-4:")
for bb, gg in ((0.4, 40.0), (0.3, 40.0)):
    d_ok = FOLD[(bb, gg)][1]
    d_bad = fold_depth(Branch(2, wh_rate(bb, gg), n=110, rtol=1e-4).Phi)
    print(f"  beta = {bb}, gamma = {gg:.0f}: fold depth in Phi {d_ok:.4f} -> {d_bad:.4f}")
print("  (the fold survives a 1e-4 integrator, which is what it should do - a fold")
print("   that only exists at tight tolerance would be an integration artefact)")
DEF_BETA_SIGN = float(abs(b_s - b_r))
DEF_EXPONENT = float(abs(s_bad[1] - MAP[(0.3, 20.0)]["eta_nis"]))'''))

# ------------------------------------------------------------------ validation
cells.append(md(r"""## Validation

**Tier 6. Nothing on this page is compared with a measurement**, and none of
what follows should be called validation against experiment. Neither source
paper reports data. What is available is a closed form, a second numerical
method, and two printed numerical statements — and they are not of equal weight.

Ranked, strongest first:

1. **The shooting reference against the exact isothermal $\eta$** — three
   geometries, seven states each. This is what licenses every branch result,
   because the branch tracer is the only route to the unstable middle branch.
2. **pymrm against that same closed form**, refined, at the observed order.
3. **pymrm against the shooting solution in the non-isothermal case**, over the
   whole branch including its unstable part, refined. The two share only the
   algebraic rate law: one is a finite-volume BVP solve in $u$, the other an
   initial-value integration in $s = \phi_0 u$ that never forms a grid.
4. **The threshold itself, two ways**: bisection on the pymrm solver for
   $\Phi = 1$, against bisection on the closed form; and, for the four other
   kinetics, against every $\Phi = 1$ state the shooting branch contains.
5. **Three printed statements from the paper reproduced**: the worked
   $\gamma = 20$, $\beta = 0.3$ case giving $\phi_0 \sim 0.1$, the multiplicity
   onset at $\beta\gamma \gtrsim 5$, and $Q'/Q$ lying between 1 and $\tfrac12$.
   None was used as an input to anything.
6. **The measured sensitivity of all of the above** — Section 6 injects a
   series of defects and reports which numbers move and which do not.
7. **The flux/volume form of eq. (11)** — provably the volume-weighted Newton
   residual, so it detects non-convergence and a source mismatch and nothing
   else, and it is reported as such rather than as corroboration.

Not on that list, because it is not evidence for anything: the pixel arithmetic
in Section 5. It measures whether Fig. 7 *could* have shown the fold the page
computes, and concludes that it could not — a statement about the figure's
resolving power, not about the physics. It is printed because the alternative,
asserting that the fold is "below the line width", is what an earlier draft did
and it was wrong by a factor of about two.

Every Newton solve behind a reported number asserts its own residual against
`RN_TOL`; the worst anywhere on the page is $5.6\times10^{-10}$."""))

cells.append(code('''print(f"1. shooting vs exact isothermal eta       : {SHOOT_EXACT:.3e}  "
      f"(3 geometries x 7 states)")
print(f"2. pymrm vs exact isothermal eta          : {PM_EXACT:.3e} at n_u = 400, "
      f"order {PM_ORDER:.2f}")
print(f"3. pymrm vs shooting, non-isothermal      : {PM_SHOOT:.3e} at n_u = 800 "
      f"over {n_states} branch states,")
print(f"                                            order {PM_SHOOT_ORDER:.2f} "
      f"(the page's real discretisation test)")
print(f"   worst Newton residual over those states : {PM_SHOOT_RN:.1e} (all converged)")
print(f"4. threshold eta*, pymrm vs closed form   : {THRESH_EXACT_DEV:.3e} relative")
print(f"   threshold eta*, pymrm vs shooting, 5 kinetics : {KIN_PYMRM_DEV:.3e} relative")
print(f"   worst Newton residual over every threshold state : "
      f"{max(THRESH_RN, KIN_RN):.1e}")
print(f"5. eq.(10) worked case, phi_0             : {PHI_W:.6f} vs the paper's "
      f"'~ 0.1'")
print(f"   multiplicity onset in phi_0            : beta*gamma = {min(PHI_ONSET):.2f} - "
      f"{max(PHI_ONSET):.2f} vs the paper's '>~ 5'")
print(f"   Q'/Q over the full phi_0 range         : [{Q_MIN:.4f}, {Q_MAX:.4f}] vs the "
      f"paper's 'between Q and Q/2'")
print(f"6. defect sensitivity, eta at Phi = 1     : nu wrong -> "
      f"{DEF['nu = 0 for a sphere'][0]:.4f}, n_u = 3 -> {DEF['n_u = 3'][0]:.4f},")
print(f"                                            rate 1 % off -> "
      f"{DEF['rate scale 1 % off'][0]:.4f}  (reference {ET_REF:.4f})")
print(f"   beta sign flipped in the exponent      : eta at Phi = 1 moves by "
      f"{DEF_BETA_SIGN:.4f}")
print(f"   eq.(10a) exponent mis-transcribed      : eta at its threshold moves by "
      f"{DEF_EXPONENT:.5f}")
print(f"7. Phi(flux) - Phi(volume)                : {FLUX_VOL:.2e}  "
      f"[= the volume-weighted Newton residual]")
print(f"   it catches a 1 % rate-scale mismatch   : {DEF_SCALE:.2e}")
print(f"   it catches an unconverged Newton solve : {BLIND_NOCONV:.1e}  "
      f"(so does the residual, {RN_NOCONV:.1e})")
print(f"   it is BLIND to a wrong nu              : {BLIND_GEOM:.1e}  "
      f"(while eta* is "
      f"{100*abs(DEF['nu = 0 for a sphere'][0]-ET_REF)/ET_REF:.0f} % wrong)")
print(f"   it is BLIND to a 3-cell grid           : {BLIND_COARSE:.1e}  "
      f"(while eta vs exact is "
      f"{DEF_COARSE_VS_EXACT:.1e})")

print("\\n--- the results themselves, which are predictions rather than agreements ---")
print(f"eta at Phi = 1, isothermal : {ETA_STAR_SLAB:.4f} slab / "
      f"{geo_rows[1][3]:.4f} cylinder / {ETA_STAR_SPHERE:.4f} sphere")
print(f"                             {KIN_MIN:.4f} to {KIN_MAX:.4f} across five kinetics "
      f"on a sphere (first crossing)")
print(f"worst eta eq.(1) calls safe : {WORST_ISO_EXO:.4g} (exothermic), "
      f"{WORST_ISO_ENDO:.4f} (endothermic)")
print(f"worst |eta-1| eq.(10a) calls safe (exothermic) : {WORST_NIS_EXO:.4f}")
print(f"Phi multivalued in {N_OBS_FOLD} of {len(FOLD)} swept cases, worst fold "
      f"{WORST_PHI_FOLD:.3f}")
print(f"every fold turning point is BELOW Phi = 1 (largest {max(TURN_PHI):.4f}), and all "
      f"{len(AMB)} ambiguous")
print(f"measurements have Phi < 1 (largest {AMB_PHI_MAX:.4f}); worst case beta = "
      f"{AMB_WORST[0]}, gamma = {AMB_WORST[1]:.0f},")
print(f"Phi = {AMB_WORST[2]:.4f} fits eta = " + ", ".join(f"{e:.4g}" for e in AMB_WORST[3]))
print(f"the one drawn Fig. 7 curve that folds does so over {FOLD_PX:.1f} px at 600 dpi "
      f"against a")
print(f"{FIG7D_STROKE_PX[0]:.0f}-{FIG7D_STROKE_PX[1]:.0f} px printed stroke - the figure "
      f"cannot resolve it either way")

report_agreement("B1.4", {
    # checks that can fail
    "shooting_vs_exact_isothermal_eta": SHOOT_EXACT,
    "pymrm_vs_exact_isothermal_eta_n400": PM_EXACT,
    "pymrm_vs_exact_observed_order": PM_ORDER,
    "pymrm_vs_shooting_nonisothermal_n800": PM_SHOOT,
    "pymrm_vs_shooting_observed_order": PM_SHOOT_ORDER,
    "threshold_eta_pymrm_vs_closed_form": THRESH_EXACT_DEV,
    # printed statements reproduced
    "eq10_worked_case_phi0": PHI_W,
    "multiplicity_onset_betagamma_min": float(min(PHI_ONSET)),
    "multiplicity_onset_betagamma_max": float(max(PHI_ONSET)),
    "Q_ratio_min_over_phi_range": Q_MIN,
    "Q_ratio_max_over_phi_range": Q_MAX,
    # the page's own results
    "eta_at_Phi1_sphere_first_order": float(ETA_STAR_SPHERE),
    "eta_at_Phi1_slab_first_order": float(ETA_STAR_SLAB),
    "eta_at_Phi1_worst_exothermic": float(WORST_ISO_EXO),
    "eta_at_Phi1_worst_endothermic": float(WORST_ISO_ENDO),
    "eq10a_worst_eta_deviation_exothermic": float(WORST_NIS_EXO),
    "observable_fold_depth_worst": WORST_PHI_FOLD,
    "observable_fold_onset_betagamma_min": float(min(OBS_ONSET)),
    # where the folds sit, and whether the figure could have shown one
    "observable_fold_turning_point_Phi_max": float(max(TURN_PHI)),
    "ambiguous_measurement_Phi_max": AMB_PHI_MAX,
    "ambiguous_measurement_worst_eta_spread": float(max(AMB_WORST[3]) / min(AMB_WORST[3])),
    "fig7d_fold_width_decades_beta03_gamma40": FOLD_DECADES,
    "fig7d_fold_width_px_at_600dpi": FOLD_PX,
    "fig7d_printed_stroke_px_at_600dpi": float(FIG7D_STROKE_PX[0]),
    "fig7d_fold_in_line_widths": float(FOLD_PX / FIG7D_STROKE_PX[0]),
    # structural identity and its measured sensitivity
    "Phi_flux_minus_Phi_volume": FLUX_VOL,
    "Phi_identity_under_1pct_rate_error": DEF_SCALE,
    "Phi_identity_blind_wrong_nu": BLIND_GEOM,
    "Phi_identity_blind_three_cells": BLIND_COARSE,
    "Phi_identity_catches_unconverged": BLIND_NOCONV,
    "threshold_eta_pymrm_vs_shooting_kinetics": KIN_PYMRM_DEV,
})'''))

# --------------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Nothing to the criterion, which is one line of algebra**, and nothing to
Weisz and Hicks' extension of it. Both are reproduced here exactly as printed.
What is added is everything the 1954 and 1962 statements left as an order of
magnitude, because in 1962 answering it meant an IBM 704 and a Fortran
subroutine distributed on paper.

- **The threshold's guarantee, as a number.** $\Phi < 1$ admits a rate
  depression of a few per cent for a first-order sphere and around 30 % for a
  slab — a factor of nearly five in the error the same inequality lets through,
  and most of it is the choice of length scale rather than physics: measured on
  the volume-to-surface length instead of the radius, the spread across
  geometries roughly halves. Across kinetics on one shape the guarantee runs
  from *exactly* 1 (zero order, where $\eta = 1$ identically up to $\Phi = 6$)
  through 0.88 for second order, to **above** 1 for a self-inhibited
  Langmuir-Hinshelwood rate — which the criterion cannot tell apart from a
  thermal effect, because it only ever sees $\Phi$. Every one of those numbers
  is printed by the cells above; none is in either source.
- **Where the isothermal form fails, quantified.** For an exothermic pellet
  eq. (1) admits observed rates thousands of times the intrinsic one before it
  objects. Eq. (10a) removes every one of those false negatives in this sweep,
  holding $\eta$ within a couple of per cent of 1 everywhere — and it pays for
  that with a conservatism that reaches five orders of magnitude in $\Phi$ at
  the severe end, rejecting conditions that are in fact clean. Both halves of
  that trade are measured here; the paper could only state the safe direction.
- **The endothermic case, where the extension inverts.** Eq. (10a)'s threshold
  is $\exp[-\gamma\beta/(1+\beta)]$, which for $\beta < 0$ exceeds 1 and grows
  without bound. Over most of the endothermic sweep — 5 of the 8 cases, the
  count the cell prints — the branch never reaches it at all, so the criterion
  certifies every attainable state as safe while eq. (1) is itself already
  admitting a depression of more than 50 %. In the remaining three the
  threshold is reached, at an $\eta$ far below 1. Either way the paper's
  one-clause restriction of eq. (10) to exothermic reactions is load-bearing,
  and this page measures what it is worth.
- **A quantitative limit on the Section V claim.** The transformation to the
  observable does not remove the multiplicity; it postpones it, by a factor of
  two to four in $\beta\gamma$. That covers every case the paper argues about in
  its text, and it does not cover the whole of the $\gamma = 10$–$40$,
  $\beta = 0$–$0.8$ range it states it computed. Where $\Phi$ folds, one
  measured rate is consistent with three steady states differing by orders of
  magnitude in $\eta$, and no threshold on $\Phi$ can separate them.
- **And every fold sits *below* $\Phi = 1$ — inside the band the criterion
  certifies as safe.** $\Phi = 1$ still has a unique root in all 32 swept cases,
  so eq. (1) keeps a single verdict boundary; what a fold costs is the meaning
  of a *pass*. All 14 ambiguous measurements on this page's own grid have
  $\Phi < 1$, and the worst of them clears the threshold by more than a factor
  of three while admitting $\eta$ = 2.3, 5.1 **or 226**. That is a second
  false-negative mechanism, independent of the loose threshold above: not a
  criterion that admits too much, but a measurement well inside the safe band
  that does not determine the state at all, and one that tightening the
  threshold cannot fix. Section 5 prints the turning points that show it.
- **What Fig. 7 can and cannot settle, in pixels.** Section 5 locates the fold
  boundary against the $\beta$ labels Fig. 7 carries: one of its four panels
  contains a drawn curve in the folded regime ($\gamma = 40$, $\beta = 0.3$).
  Whether the printed figure contradicts the computed fold is a resolving-power
  question, so it is measured rather than asserted — the fold spans about 7 px
  on a 600 dpi render of panel (d), against a printed stroke of 6–7 px. It is
  *comparable to* the line width, not below it, so the figure is **silent**: it
  can neither show the fold nor rule it out, and its monotone appearance is not
  evidence against the computation. This is a limit on the claim, not a
  refutation of the paper: the criterion works everywhere the authors apply it.

**And the whole page is the inverse problem, which is what pymrm makes cheap.**
$\Phi$ is the dimensionless surface flux, so a finite-volume solver reports it
without ever being told the rate constant — `Pellet.Phi` takes no $\phi_0$
argument. Bisecting on $\phi_0$ until $\Phi$ hits a target is then a three-line
inversion of a measurement into an effectiveness factor, which is exactly what
the paper's Fig. 7 is for and what a laboratory actually needs."""))

# --------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

| Change | Where |
|---|---|
| Different geometry | `Pellet(geom=...)`, or a callable `nu` in `construct_div` |
| Different kinetics | pass any `R(y)` to `Pellet.solve`; `NumJac` picks up the derivatives |
| Your own measurement | `Branch(nu, R).all_where(lambda p, e, P: P - your_Phi)`; the inversion needs no $k$, and it returns *every* consistent state |
| Non-isothermal | `wh_rate(beta, gamma)`, which is the reduced eq. (9) |
| Two fields instead of one | [`B1.6`](../B1.6-prater-relation/)'s `TwoFieldPellet`, if you do not want the Prater relation substituted |
| The unstable branch | `Branch`, parametrised by $y_c$, not $\phi_0$ |

**How to actually use the criterion, given what is above.** Compute
$\Phi = r_{\mathrm{obs}}L^2/(D_ec_s)$ and **state the length scale** — the
threshold and the length are not separable, and quoting $\Phi < 1$ without
saying whether $L$ is $R$ or $V/S$ leaves a factor of nine on a sphere. If the
reaction is exothermic, estimate $\beta$ and $\gamma$ and use eq. (10a), not
eq. (1); eq. (1) is not conservative for an exothermic pellet in any useful
sense. If $\beta\gamma$ is above about 10, stop: the observable no longer
determines the state and the criterion has nothing to say. And note that a
passing $\Phi$ still admits several per cent in $\eta$ and a comparable error in
the apparent activation energy — Section 3 prints both, per geometry.

**Reading Weisz and Hicks.** The scan's text layer is unusable for equations —
eq. (1) OCRs as `$:$<I` and eq. (10a) as `expCrBK1 + LOI`. Render the journal
page at 600 dpi and read it. The prose OCRs acceptably; the numbers and symbols
do not. Same trap as `B3.1` and `B1.6`, and it is in `docs/pdf-findings.md`.

**Related pages.** [`B1.1`/`B1.5`](../B1.1-thiele-weisz-hicks/) — the forward
problem, $\eta(\phi_0)$ and the multiplicity this page views through the
observable; [`B1.6`](../B1.6-prater-relation/) — the Prater relation that
reduces the two balances to the one equation used here; `B1.2` Aris generalised
modulus and `B1.3` Bischoff generalised modulus — the length-scale question
Section 3 runs into; `B1.7` Mears criteria — the same observable-only logic
extended to external gradients.

## References

Weisz, P. B. and Prater, C. D. (1954). In *Advances in Catalysis*, Vol. 6,
p. 143. Academic Press, New York — **the origin the criterion is named for; not
consulted**, no copy on disk and no open-access route. Given exactly as Weisz
and Hicks' reference [5] prints it on journal page 274, which carries **no
chapter title**; none is supplied here, because the only document consulted for
this page does not contain one.

Weisz, P. B. (1957). *Z. Phys. Chem.* **11**, 1 — the reference Weisz and Hicks
themselves attach to eq. (1) ("It was shown by WEISZ [9] that the conditions
…"). **Not consulted.**

Weisz, P. B. and Hicks, J. S. (1962). *The behaviour of porous catalyst
particles in view of internal mass and heat diffusion effects.* Chemical
Engineering Science **17**(4), 265-275.
[doi:10.1016/0009-2509(62)85005-2](https://doi.org/10.1016/0009-2509(62)85005-2)
— **the source actually read.** Eqs. (1), (9), (10), (10a) and (11), the
$Q'/Q$ relation, the $\beta\gamma \gtrsim 5$ statement and the Section V claim
were all read off 600 dpi renders of journal pages 265, 267, 268, 269, 272 and
273; the two references above were read from the literature-cited block on
journal page 274, and the Fig. 7(d) axis calibration and stroke width of
Section 5 were measured on the render of page 272. Nothing on this page comes
from the scan's text layer, which is unusable for equations."""))

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
