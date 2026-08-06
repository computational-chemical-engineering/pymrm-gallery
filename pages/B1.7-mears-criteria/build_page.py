#!/usr/bin/env python3
"""Generate index.ipynb for page B1.7 (the Mears criteria). Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "The Mears criteria: what each transport-limit test actually guarantees"
description: "Mears (1971) collects every diagnostic inequality that decides whether a measured rate is chemistry or transport, each claiming a 5 % tolerance. This page solves the pellet, film and dispersion problems exactly and measures what each criterion's boundary really admits - the one that is exact, the ones the rounding broke, the window where one loses all power, and the limit where the combined criterion fails structurally."
categories: [sec:B, struct:S3, struct:S4, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-06
---

# The Mears criteria: what each transport-limit test actually guarantees

**Catalog ID:** `B1.7` · **Structures:** `S3` (1D steady BVP) + `S4` (axial dispersion) · **Tier:** T0

An experimenter measures a rate. Before that number can be called kinetics,
gradients must be ruled out in three domains — *intraparticle* (inside the
pellet), *interphase* (across the film around it), and *interparticle* (across
the reactor). Mears (1971) collects, and in three cases derives, one-line
inequalities in **observables** for each domain, almost all carrying the same
promise: satisfy this, and the measured rate deviates from the true one by less
than 5 %.

Those inequalities are still the standard screening tests fifty years later.
What none of them comes with is a measurement of what the promise is worth: the
derivations are first-order perturbation expansions, several right-hand sides
are rounded, and one is stated with an unquantified exclusion zone. The
published [`B1.4`](../B1.4-weisz-prater-criterion/) did this audit for the
Weisz-Prater criterion alone; this page does it for the rest of the family.
Every criterion the paper prints is solved against exactly, and the page
reports, per criterion, the deviation its boundary actually admits:

1. **One criterion is exact** — the interphase mass test (eq. 17) holds its 5 %
   to all digits for first-order kinetics, because the film balance makes the
   deviation *equal* to a third of the criterion's own left side.
2. **The intraparticle family is uniformly liberal by the same rounding.** The
   paper says "a numerical value of 0.75 … was rounded to 1". That rounding is
   exactly the 5 % claim: at the printed bounds the exact deviation is 6.3-6.8 %;
   at the unrounded 0.75 it is 4.8-5.1 %.
3. **The nonisothermal criterion (eq. 9) has a measurable window of no
   power.** Its bound diverges at the compensation point $n = \gamma\beta$;
   the exact 5 % threshold stays finite, so near compensation the criterion
   certifies states with arbitrarily large deviations — and this page computes
   the window of $\gamma\beta$ in which even restoring the unrounded 0.75
   cannot save it.
4. **The combined interphase-intraparticle criterion (eq. 21) fails
   structurally at low Biot number**, and even passing *every* pellet-scale
   criterion at once still admits ~11 % — the worst case puts two resistances
   at their boundaries at once, because 5 % tolerances stack.
5. **The axial-dispersion criterion (eq. 30) — Mears' own — is the sharpest in
   the paper**: it holds its 5 % at every conversion tested and approaches it
   asymptotically from below."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

The forward problem — given kinetics, compute the transport-corrupted rate — is
what [`B1.1`](../B1.1-thiele-weisz-hicks/) solves. The experimenter faces the
inverse problem: the intrinsic kinetics are what the experiment is *for*, so a
usable test must be stated in measured quantities only. The classic example is
Weisz and Prater's group $\Phi = \mathcal{R}\,r_p^2/(C_s D_e)$ — observed rate
per particle volume, particle radius, surface concentration, effective
diffusivity — with the criterion $\Phi < 1$.

Mears' 1971 review assembles the whole toolbox in that spirit and extends it.
Its route to each criterion is *perturbation analysis*: expand the rate in a
Taylor series in the concentration and temperature perturbations, and bound the
group of observables that keeps the first-order rate deviation below a stated
tolerance, usually 5 %. That construction has two systematic consequences this
page measures. First, a first-order expansion is only as good as the curvature
it drops, so each boundary sits slightly off the true 5 % contour — sometimes
on the safe side, sometimes not, and in the intraparticle family the deliberate
rounding of 0.75 to 1 moves it further. Second, wherever two effects cancel at
first order (the $n = \gamma\beta$ compensation case), the criterion's
denominator vanishes and the test silently loses all power — the paper says the
case "requires further consideration" and this page says how wide it is.

The paper contains no data, no figures and no tables: every claim in it is an
inequality plus prose. That makes it exactly the kind of source the gallery can
audit completely — the criteria are inequalities between computable quantities,
so the honest validation is to solve the underlying pellet, film and reactor
problems exactly and measure where each stated tolerance holds, where it fails,
and by how much."""))

# ------------------------------------------------------------ the published model
cells.append(md(r"""## The published model

### Provenance: what was read, and where each criterion actually comes from

**Mears, D. E. (1971), "Tests for Transport Limitations in Experimental
Catalytic Reactors", *Ind. Eng. Chem. Process Des. Develop.* 10(4), 541-547, is
on disk and is the only document consulted.** The scan is CCITT-G4 bilevel at
300 ppi native; every equation and every numeric threshold below was read off
cropped 300 dpi renders (`pdftoppm -r 300`) of the journal page named, because
the text layer eats exactly the load-bearing symbols: the effectiveness factor
extracts as a bare comma, eq. (25)'s right side `5.3` comes back as `°-`, and
eq. (3)'s $1/|n|$ loses its absolute-value bars. (The PDF also opens with the
tail of the *preceding* article — a catalytic-cracking study — so identity was
established from page 2 onward, the running head, and the ACS download stamp
`article-pdf/10/4/541/`.)

One attribution matter the page must be honest about, the same class of finding
as on [`B1.4`](../B1.4-weisz-prater-criterion/): **this paper is a review, and
the three criteria "derived by the author" are derived in two companion papers,
not here.** The interphase criteria (eqs. 14, 17) and the interparticle radial
criterion (eq. 23) cite "Mears, 1971b" — *J. Catal.* **20**, 127; the axial
criterion (eqs. 29-32) cites "Mears, 1971a" — *Chem. Eng. Sci.*, listed "in
press" in the reference list. **Neither companion was consulted**; what this
page audits are the criteria exactly as the review states them, which is also
how the literature cites them. The remaining criteria are attributed by the
review to Weisz & Prater (1954), Weisz (1957), Hudgins (1968), Anderson (1963),
Kubota & Yamanaka (1969), Weisz & Hicks (1962), Carberry (1961), Satterfield et
al. (1969), van den Bleek et al. (1969) and Petersen (1965b, 1968) — none
consulted either.

### The criteria, as printed

All symbols as in the paper's Nomenclature (journal page 547): $\mathcal{R}$ is
the **observed** rate per unit particle volume, $r_p$ the particle **radius**,
$d_p$ its diameter, $C_s$/$C_b$ the surface/bulk reactant concentration, $D_e$
the effective diffusivity, $\lambda$ the particle conductivity, $h$/$k_c$ the
film heat/mass transfer coefficients, $n$ the power-law order, $E$ the true
activation energy.

**Intraparticle, isothermal** (journal pages 542-543):

$$\frac{\mathcal{R}r_p^2}{C_sD_e} < 1 \tag{1}$$

(Weisz-Prater; first order, sphere, $\eta \ge 0.95$), with Weisz (1957)'s
linear-approximation values "conservatively 0.3 for second-order reactions and
6 for zero order", 0.6 for first order;

$$\frac{\mathcal{R}r_p^2}{C_sD_e} < \frac{\mathcal{R}_s}{C_s\mathcal{R}'_{cs}} \tag{2}$$

(Hudgins 1968, $\eta > 0.95$, "a numerical value of 0.75 … was rounded to 1",
extendable to negative orders "by taking the absolute value of the
derivative"), which for power-law kinetics becomes

$$\frac{\mathcal{R}r_p^2}{C_sD_e} < \frac{1}{|n|} \tag{3}$$

(eqs. 2 and 3 "fail for zero-order kinetics … the value of 6 on the right
obtained by Weisz is recommended").

**Intraparticle, thermal** (page 542): Anderson's quasi-isothermality test,
"the observed rate $\mathcal{R}$ must not differ from the rate that would
prevail at constant temperature by more than an acceptable amount, say 5 %",
with "again a numerical value of 0.75 on the right … rounded to 1", valid
"whether diffusional limitations exist in the particle or not":

$$\frac{|\Delta H|\mathcal{R}r_p^2}{\lambda T_s} < \frac{T_sR}{E} \tag{6}$$

**Intraparticle, combined concentration + temperature** (page 543): for
power-law kinetics (except zero order), $\eta = 1 \pm 0.05$ if

$$\frac{\mathcal{R}r_p^2}{C_sD_e} < \frac{1}{|n - \gamma\beta|},
\qquad \gamma = \frac{E}{RT_s},
\qquad \beta = \frac{(-\Delta H)D_eC_s}{\lambda T_s} \tag{9-11}$$

with the caveat that "the special case for which $|n-\gamma\beta|$ is close to
or equal to zero requires further consideration": there "the heat effect
compensates for the diffusion effect so that the rate stays nearly constant
until the concentration within the catalyst has dropped to less than 80 % of
the surface value", and Petersen's asymptotic method "gives a value of 13 for
the right side of Equation 9 when $n$ and $\gamma\beta$ both equal 1.0". The
isothermal-particle corollary is

$$|\gamma\beta| < 0.05\,n \tag{13}$$

and the comparison case is Weisz and Hicks (1962) for exothermic first-order
reactions, which "safely predicts positive deviations ($\eta > 1.05$) … but
fails to predict negative deviations for endothermic reactions":

$$\frac{\mathcal{R}r_p^2}{C_sD_e}\exp\!\left[\frac{\gamma\beta}{1+\beta}\right] < 1 \tag{12}$$

**Interphase** (pages 543-544): Mears' film-heat criterion, "if the observed
rate is to deviate by less than 5 %",

$$|\chi| = \left|\frac{-\Delta H\,\mathcal{R}\,r_p}{h\,T_b}\right| < 0.15\,\frac{RT_b}{E} \tag{14}$$

valid "whether transport limitations exist in the particle or not"; the film
heat resistance limits before the intraparticle one as long as

$$\mathrm{Bi}_p = h\,d_p/\lambda < 10 \tag{15}$$

Carberry (1961)'s first-order isothermal mass test $\eta k/k_ca < 0.1$
(eq. 16), "generalized to other reaction orders by the perturbation approach,
allowing a 5 % deviation as before":

$$\omega = \frac{\mathcal{R}\,r_p}{C_b\,k_c} < \frac{0.15}{n} \tag{17}$$

("the criterion also shows that zero-order reactions are not affected by
interphase transport"), and Satterfield et al.'s trickle-flow analogue
$\mathcal{R}r_p/C^*k_{LS} < 0.15$ (eq. 18).

**Combined interphase-intraparticle** (page 544): replacing the Dirichlet
surface conditions by film conditions (eqs. 19-20), "excluding $|n -
\gamma_b\beta_b|$ close to zero, the resulting criterion for $\eta = 1 \pm
0.05$ becomes"

$$\frac{\mathcal{R}r_p^2}{C_bD_e} <
\frac{1 + 0.33\,\gamma\chi}{|n-\gamma_b\beta_b|\,(1 + 0.33\,n\omega)} \tag{21}$$

with $\gamma_b, \beta_b$ at bulk conditions, and its isothermality corollary
$|\gamma_b\beta_b + 0.3\,n\gamma_b\chi| < 0.05\,n$ (eq. 22).

**Interparticle** (page 545): the radial-heat criterion at the hot spot,
$|\Delta H|\mathcal{R}_bR_o^2/(k_eT_w) < 0.4\,RT_w/E$ (eq. 23, wall-resistance
variant eq. 26, ordering rule eq. 25 with right side 5.3), the dilution
criterion $L/d_p > 250\,b/\delta$ (eq. 27), and the axial-dispersion family:
Petersen's asymptotic $\alpha = \sqrt{k_bD_a}/v < 1$ (eq. 28), against Mears'
perturbation criterion holding "the deviation in the required reactor length to
less than 5 %",

$$\alpha < 0.22 \tag{29}
\qquad\Longleftrightarrow\qquad
\frac{L}{d_p} > \frac{20\,n}{\mathrm{Pe}_a}\ln\frac{C_o}{C_f} \tag{30, 32}$$

with $\mathrm{Pe}_a = v d_p/D_a$, "about four and a half times more
conservative than the earlier one" because "asymptotic criteria actually
indicate the transition region … rather than the point at which transport
effects are just starting to become significant".

### What is audited, and what is stated only

The pellet-scale and axial criteria — eqs. 1-3, 6, 9, 12, 13, 14, 16, 17, 21
(isothermal part), 28-32 — are audited against exact solves below. Three groups
are **stated, not audited**, each for a reason given here rather than left to
be discovered: eqs. 4-5 (volume-change modulus) would need a pellet flux model
with net convection, and the paper's own position is that "in dilute cases the
volume-change effect becomes negligible" — with one printed wrinkle to carry:
the paper prints $\theta$ with opposite signs, $\theta = (\nu-1)Y_s$ in eq. 4
(page 542) against $\theta = (1-\nu)Y_s$ in the nomenclature (page 547), so a
builder of that criterion must pick a sign the source does not settle;
eqs. 23-27 (interparticle) test a
*reactor-scale* temperature field against $k_e$ correlations and a stochastic
bypassing model that live in other papers, so an exact audit has no
self-contained reference to solve; and eq. 22's thermal part couples film heat
to the pellet interior, whose honest audit belongs with eq. 21's thermal part —
both are beyond the isothermal reduction audited here, and the page says so
where it reports eq. 21. Eq. 18 is eq. 17 with $C^*k_{LS}$ in place of
$C_bk_c$ at $n = 1$, so its audit is the $n = 1$ column of eq. 17's."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

Everything is audited in the dimensionless variables the criteria themselves
use, so no physical parameter values enter any audit except the two
deliberately dimensional cross-checks in Section 3.

**Pellet (intraparticle).** Sphere, steady state, constant $D_e$ and $\lambda$,
single irreversible reaction, uniform activity — the paper's own assumptions.
With $y = c/C_s$, $u = r/r_p$, the two balances reduce by the Prater relation
(proved and stress-tested on [`B1.6`](../B1.6-prater-relation/), not re-derived
here) to

$$\nabla^2 y = \phi^2\,\mathcal{F}(y),\qquad
\mathcal{F}(y) = y^n\exp\!\left[\frac{\gamma\beta(1-y)}{1+\beta(1-y)}\right],
\qquad y'(0) = 0,\; y(1) = 1,$$

with $\nu = 2$ in the Laplacian. All rate laws are normalised to
$\mathcal{F}(1) = 1$, so the observable is exactly the printed group in every
case:

$$\Phi \equiv \frac{\mathcal{R}r_p^2}{C_sD_e} = \phi^2\eta = 3\,y'(1),
\qquad \eta = \frac{3\,y'(1)}{\phi^2}.$$

The **deviation** eqs. 1-3, 9, 12 bound is $|\eta - 1|$ — the observed rate
against the rate at surface conditions ("$\eta = 1 \pm 0.05$"). The deviation
eq. 6 bounds is different and is measured differently: the observed rate
against "the rate that would prevail at constant temperature", i.e.
$\eta/\eta_{\rm iso}(\phi) - 1$ at the same $\phi$.

**Film (interphase).** Resistances lumped at the surface (the paper's stated
assumption), sphere, so the external area per particle volume is $a = 3/r_p$
(nomenclature page: $a$ = "superficial (outside) surface area of catalyst
particle per unit particle volume"). The steady film balances then tie the
surface state to the criteria's own left sides **exactly**:

$$1 - \frac{C_s}{C_b} = \frac{\mathcal{R}r_p}{3\,k_cC_b} = \frac{\omega}{3},
\qquad
\frac{T_s - T_b}{T_b} = \frac{(-\Delta H)\mathcal{R}r_p}{3\,h\,T_b} = \frac{\chi}{3},$$

because $\mathcal{R}$ in $\omega$ and $\chi$ is the *observed* rate. The
deviation eqs. 14, 16, 17 bound is the observed rate against the rate at bulk
conditions, which is therefore a closed-form function of the observables — the
whole interphase audit is exact algebra, and the numerics only cross-check it.

**Combined (eq. 21).** The same pellet with the Dirichlet surface condition
replaced by the film condition $y'(1) = \mathrm{Bi}_m(1 - y(1))$,
$\mathrm{Bi}_m = k_cr_p/D_e$, everything now normalised on bulk:
$\Phi_b = \mathcal{R}r_p^2/(C_bD_e) = 3y'(1)$, $\eta_b$ against bulk
conditions, and $\omega = \Phi_b/\mathrm{Bi}_m$ — all observables again.

**Axial (eqs. 28-32).** Isothermal single reaction in a closed vessel with
axial dispersion, Danckwerts boundary conditions (the canonical treatment is
[`A2.1`](../A2.1-danckwerts-boundary-conditions/)), power-law rate. The
deviation eq. 30 bounds is in the **required reactor length**, and two
readings of "required length" exist — both are computed in Section 6, and the
one the headline numbers use is stated here so definition and implementation
match. At the criterion's own boundary
$L/d_p = (20n/\mathrm{Pe}_a)\ln(C_o/C_f)$ the vessel Peclet number is
$\mathrm{Pe}_L = vL/D_a = 20\,n\ln(C_o/C_f)$ — a function of conversion alone.
The headline numbers hold $\mathrm{Pe}_L$ at that boundary value and measure
the $k\tau$ deficit of that vessel against plug flow — the perturbation-theory
object, whose first-order term $\ln(C_o/C_f)/\mathrm{Pe}_L$ is exactly 5 % at
the boundary. The fully self-consistent reading — extra length at fixed $k$,
$v$, $D_a$, with $\mathrm{Pe}_L$ co-varying as $L$ grows — is printed beside
it in Section 6 and is slightly smaller everywhere.

**Parameter ranges.** Orders $n \in \{-1, \tfrac12, 1, 2, 3\}$ plus zero order
and one inhibited Langmuir-Hinshelwood rate; $\gamma = 10$-$40$ (the paper's
"typical values 10-40"); $\gamma\beta$ from $-2$ to $+2$, inside the paper's
$\beta$ range "$-0.1$-$0.1$" at these $\gamma$; $\mathrm{Bi}_m = 0.3$-$100$;
conversions 0.2-0.99. Carried without test: steady state, sphere geometry (the
criteria are all stated for $r_p$ of a sphere), constant properties, and the
Prater relation inside the pellet."""))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

**There is none, and there cannot be — provenance tier 6.**

The source paper prints no measurement, no figure and no table anywhere in its
seven pages (checked page by page on the renders; the only tabular matter is
the Nomenclature). Its claims are inequalities plus derivation prose, so
nothing on this page is validated against experiment and nothing here should be
read as if it were. What follows is: exact closed forms, two independent
numerical methods solving the same stated models, and the paper's own printed
numbers and internal identities used as transcription checks.

No dataset is shipped and no figure is digitised. No other page's CSV is
loaded either — the cross-page findings this page relies on
([`B1.4`](../B1.4-weisz-prater-criterion/)'s fold results, quoted in Section 4)
are cited from that page's published text, and this page's own branches
reproduce the one number the two pages share (the $\eta$ at $\Phi = 1$ for the
inhibited Langmuir-Hinshelwood rate), printed side by side where it occurs."""))

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

import time
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu
from pymrm import (construct_grad, construct_div, construct_convflux_upwind,
                   interp_cntr_to_stagg_tvd, vanleer, compute_boundary_values,
                   NumJac, newton, clip_approach)
from gallery_utils import report_agreement

T_START = time.time()
np.seterr(all="ignore")
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

Three solvers, and a deliberate division of labour between them.

**`shoot`** is the deterministic reference for every pellet threshold.
Substituting $s = \phi u$ removes $\phi$ from the pellet equation, so
integrating outward from a prescribed centre value $y_c$ until $y = 1$
(Dirichlet) — or until $s\,y'(s) = \mathrm{Bi}_m(1-y)$ (film surface) — yields
$\phi$ *as the stopping radius*. Parametrising by $y_c$ instead of $\phi$
traverses a solution branch in one deterministic pass with no continuation
chain, which is the lesson [`B1.1`](../B1.1-thiele-weisz-hicks/) paid for: CI
once reproduced a different ignited-branch $\eta$ because a feature had been
located on a warm-start path. Every threshold on this page is bracketed on a
fixed branch grid and refined with Brent's method on the shooting solver
itself.

**`Pellet`** is the pymrm finite-volume solver, the production tool the page is
really about, checked against `shoot` and against closed forms at measured
order. Its conventions all matter and all follow the house rules:

- **Boundary conditions use the OUTWARD normal.** Centre: symmetry
  $\partial y/\partial n = 0$ is `{a:1, b:0, d:0}`. Surface, Dirichlet:
  `{a:0, b:1, d:1}`. Surface, film: the flux condition
  $\mathrm{d}y/\mathrm{d}u = \mathrm{Bi}_m(1-y)$ has $n = +u$ there, so
  $\partial y/\partial n + \mathrm{Bi}_m\,y = \mathrm{Bi}_m$, i.e.
  `{a:1, b:Bim, d:Bim}` — the physical equation sits in a comment next to each.
- **Operators are assembled once** in `__init__`; only the pointwise source
  block is rebuilt inside Newton.
- **Layout `(n_u, 1)`, spatial axis first, field last** — `NumJac((n_u,))` on a
  one-field problem couples the last axis in full and builds a dense Jacobian
  (the trap measured on `B1.1`/`B1.6`/`F3.1`; bit-identical answers, 6× cost).
- The observable $\Phi = 3\,y'(1)$ is read from the boundary-respecting
  gradient row, never from the last cell centre.

**`axial_solve`** is the closed-vessel dispersion reactor on the pymrm
operators, following the published [`A2.1`](../A2.1-danckwerts-boundary-conditions/)
assembly: `construct_convflux_upwind` plus a van Leer TVD deferred correction
for second order, Danckwerts inlet `{a: D/v, b:1, d:1}` (outward normal points
against the flow, so $D/v\,\partial c/\partial n + c = c^\ast$), zero-gradient
outlet `{a:1, b:0, d:0}`. The outlet concentration is read through
`compute_boundary_values` — pymrm's zero-gradient outflow *extrapolates to the
face*, and reading the last cell centre instead disagrees at first order in the
grid (the `A3.7` lesson; the break table measures it here). The nonlinear
$n = 2$ case wraps the same assembly in a Picard loop with an asserted
convergence check."""))

cells.append(code('''# ----------------------------------------------------------- shooting reference
def shoot(R, t, nu=2, rtol=1e-10, bim=None, dense=False):
    """Integrate y'' + (nu/s) y' = R(y) outward from y(0) = y_c = 1 - exp(-t).

    Dirichlet surface (bim=None): stop at y = 1; the stopping radius IS phi.
    Film surface (bim=Bi_m)     : stop where s y'(s) = Bi_m (1 - y); y is c/C_b.
    Returns (phi, eta, Phi, omega, y_surf[, sol]) with eta and Phi on the
    surface (Dirichlet) or bulk (film) normalisation; omega = Phi/Bi_m or None.
    Shares NOTHING with the pymrm route but the algebraic rate law.
    """
    y_c = -np.expm1(-t)
    r0 = float(R(y_c))
    if not (0.0 < y_c < 1.0) or r0 <= 0.0:
        return None
    s0 = min(1e-3, np.sqrt(2 * (nu + 1) * 1e-10 * min(y_c, np.exp(-t)) / abs(r0)))
    rhs = lambda s, v: [v[1], float(R(np.clip(v[0], 0.0, 1.0))) - (nu / s) * v[1]]
    if bim is None:
        ev = lambda s, v: v[0] - 1.0
    else:
        ev = lambda s, v: s * v[1] - bim * (1.0 - v[0])
    ev.terminal, ev.direction = True, 1
    v0 = [y_c + r0 * s0 ** 2 / (2 * (nu + 1)), r0 * s0 / (nu + 1)]
    sol = solve_ivp(rhs, (s0, 1e5), v0, events=ev, rtol=rtol, atol=1e-16,
                    method="DOP853", dense_output=dense)
    if not sol.t_events[0].size:
        return None
    s = float(sol.t_events[0][0])
    y_s, yp = (float(x) for x in sol.y_events[0][0])
    out = (s, (nu + 1) * yp / s, (nu + 1) * s * yp,
           None if bim is None else (nu + 1) * s * yp / bim, y_s)
    return out + (sol,) if dense else out


def cross(R, val, lo=1e-6, hi=20.0, n_scan=90, nu=2, bim=None, rtol=1e-10,
          which="first"):
    """States on the branch where val(phi, eta, Phi, omega, y_s) crosses zero.

    Scans a FIXED log grid in the branch parameter t (deterministic - no
    continuation), brackets every sign change, refines each with brentq on the
    shooting solver itself. The scan runs from the DEEP end of the branch
    (small t = small y_c) toward the kinetic end, so which='first' returns the
    DEEPEST crossing - the largest-Phi boundary state, the deep-side entry
    into the condition's band. That is NOT the crossing a slowly worsening
    experiment meets first (the smallest-Phi one); wherever a condition has
    several crossings the page uses which='all' and prints the structure
    (Section 4 - it matters for eq. 9 above gamma*beta ~ 1.45).
    """
    ts = np.logspace(np.log10(lo), np.log10(hi), n_scan)
    f = lambda t: (lambda o: np.nan if o is None else val(*o))(
        shoot(R, t, nu=nu, bim=bim, rtol=rtol))
    v = np.array([f(t) for t in ts])
    out = []
    for i in range(len(v) - 1):
        if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i] * v[i + 1] < 0:
            tt = brentq(f, ts[i], ts[i + 1], xtol=1e-13, rtol=8.9e-16)
            out.append((tt,) + shoot(R, tt, nu=nu, bim=bim, rtol=rtol))
            if which == "first":
                return out[0]
    return (out if out else None) if which == "all" else None


# ------------------------------------------------------------------ rate laws
def power(n):
    return lambda y: np.clip(y, 1e-9 if n < 0 else 0.0, 1.0) ** n


def wh(n, gb, gamma):
    """Power-law Weisz-Hicks kinetics, parametrised by the PRODUCT gb = gamma*beta."""
    beta = gb / gamma
    return lambda y: (np.clip(y, 0.0, 1.0) ** n
                      * np.exp(gb * (1 - y) / (1 + beta * (1 - y))))


def eta_iso1(phi, nu=2):
    """Closed-form isothermal first-order effectiveness factor."""
    if nu == 0:
        return np.tanh(phi) / phi
    return (3.0 / phi) * (1.0 / np.tanh(phi) - 1.0 / phi)


# --------------------------------------------------------- pymrm pellet solver
class Pellet:
    """The pellet BVP on pymrm operators. Layout (n_u, 1): space first, field last."""

    def __init__(self, n_u=400, nu=2, bim=None):
        self.n_u, self.nu, self.bim, self.shape = n_u, nu, bim, (n_u, 1)
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])
        if bim is None:
            #   u = 1 (surface): Dirichlet  y = 1
            bc_s = {"a": 0.0, "b": 1.0, "d": 1.0}
        else:
            #   u = 1 (surface): film flux  dy/du = Bi_m (1 - y);  n = +u there,
            #   so  dy/dn + Bi_m y = Bi_m
            bc_s = {"a": 1.0, "b": bim, "d": bim}
        #   u = 0 (centre): symmetry  dy/dn = 0
        bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, bc_s)
        g, gb_ = construct_grad(self.shape, self.u_f, self.u_c, bc)
        d = construct_div(self.shape, self.u_f, nu=nu)   # nu = 2: sphere
        self.lap, self.lap_bc = d @ g, (d @ gb_).toarray().reshape(-1, 1)
        self.grad, self.grad_bc = g, gb_.toarray().reshape(-1, 1)
        self.numjac = NumJac(self.shape)                 # last axis = the field

    def solve(self, phi, R, y_init=None, maxfev=80):
        p2 = phi * phi
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

    def Phi(self, y):
        """The observable: (nu+1) x the surface flux, from the boundary-respecting
        gradient row - never the last cell centre."""
        f = self.grad @ y.reshape((-1, 1)) + self.grad_bc
        return (self.nu + 1.0) * float(f[-1].item())

    def eta(self, y, phi):
        return self.Phi(y) / phi ** 2


RN_TOL = 1e-8      # every pymrm solve behind a reported number must beat this
print("operators assembled once per Pellet; NumJac shape (n_u, 1); RN_TOL =", RN_TOL)'''))

cells.append(code(r'''# ------------------------------------------------ axial dispersion, pymrm route
def axial_solve(Pe_L, ktau, n_react=1, n_z=800, inlet="danckwerts",
                outlet_read="boundary", limiter=vanleer, max_pic=200):
    """Closed vessel, 0 <= z <= 1, v = 1, D = 1/Pe_L, source ktau * c^n_react.

    Follows the published A2.1 assembly: upwind convection + van Leer deferred
    correction, Danckwerts BCs on the outward normal. Returns the outlet
    concentration (c_out/c*), read at the boundary FACE via
    compute_boundary_values unless outlet_read='last-cell' (a deliberate defect
    for the break table). Nonlinear n_react solved by Picard with an asserted
    convergence check.
    """
    shape = (n_z, 1)
    z_f = np.linspace(0.0, 1.0, n_z + 1)
    z_c = 0.5 * (z_f[:-1] + z_f[1:])
    D = 1.0 / Pe_L
    if inlet == "danckwerts":
        # v c* = v c - D dc/dz at z=0; n = -z so dc/dn = -dc/dz:
        #   (D/v) dc/dn + c = c*  ->  a = D/v, b = 1, d = 1
        bc_in = {"a": D, "b": 1.0, "d": 1.0}
    else:                                   # the naive Dirichlet inlet (defect)
        bc_in = {"a": 0.0, "b": 1.0, "d": 1.0}
    # outlet dc/dz = 0; n = +z so dc/dn = dc/dz:  a = 1, b = 0, d = 0
    bc = (bc_in, {"a": 1.0, "b": 0.0, "d": 0.0})
    conv, conv_bc = construct_convflux_upwind(shape, z_f, z_c, bc, v=1.0)
    grad, grad_bc = construct_grad(shape, z_f, z_c, bc)
    div = construct_div(shape, z_f, nu=0)
    A0 = div @ (conv - D * grad)
    b0 = np.asarray((div @ (conv_bc - D * grad_bc)).todense()).ravel()
    c = np.full(n_z, 0.5)
    done = False
    for _ in range(max_pic):
        # linearise the source as ktau * c_old^(n-1) * c  (exact for n = 1)
        w = ktau * np.clip(c, 0.0, None) ** (n_react - 1)
        A = (A0 + eye_array(n_z, format="csc") * 0.0
             + __import__("scipy.sparse", fromlist=["diags"]).diags(w, format="csc"))
        lu = splu(A.tocsc())
        _, dc_f = interp_cntr_to_stagg_tvd(c.reshape(shape), z_f, z_c, bc, 1.0,
                                           tvd_limiter=limiter, axis=0)
        c_new = lu.solve(-b0 - np.asarray(div @ dc_f.reshape(-1, 1)).ravel())
        done = np.max(np.abs(c_new - c)) < 1e-13
        c = c_new
        if done:
            break
    assert done, "Picard/deferred-correction loop did not converge"
    if outlet_read == "boundary":
        _, _, c_out, _ = compute_boundary_values(c.reshape(shape), z_f, z_c, bc)
        return float(np.asarray(c_out).ravel()[0]), c, z_c
    return float(c[-1]), c, z_c            # the last-cell-centre defect


# ------------------------------------ closed form (Danckwerts / Wehner-Wilhelm)
def frac_unconverted(ktau, Pe_L):
    """First order, closed vessel: exit c/c*, the classical closed form."""
    q = np.sqrt(1.0 + 4.0 * ktau / Pe_L)
    den = (1 + q) ** 2 * np.exp(q * Pe_L / 2) - (1 - q) ** 2 * np.exp(-q * Pe_L / 2)
    return 4 * q * np.exp(Pe_L / 2) / den


def axial_dev(X, n_react=1, n_z=None):
    """Deviation in REQUIRED length at eq. 30/32's own boundary, conversion X.

    At the boundary Pe_L = 20 n ln(1/(1-X)). Find the ktau the dispersed vessel
    needs for conversion X, compare with plug flow. Closed form for n = 1;
    the pymrm reactor (n_z cells) for anything else or when n_z is given.
    """
    lnr = np.log(1.0 / (1.0 - X))
    Pe_L = 20.0 * n_react * lnr
    ktau_plug = lnr if n_react == 1 else X / (1.0 - X)   # plug flow, C0 = 1
    if n_z is None and n_react == 1:
        f = lambda k: frac_unconverted(k, Pe_L) - (1.0 - X)
    else:
        nz = n_z or 800
        f = lambda k: axial_solve(Pe_L, k, n_react=n_react, n_z=nz)[0] - (1.0 - X)
    ktau = brentq(f, ktau_plug, ktau_plug * 2.0, xtol=1e-13, rtol=8.9e-16)
    return ktau / ktau_plug - 1.0


print("axial machinery ready")'''))

# --------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. The transcription, checked against the paper's own identities

Every threshold was read off a native-resolution render, and a mis-read
threshold silently changes every verdict downstream — so before anything is
built on them, the criteria are checked against each other. The paper happens
to pay for four checks with real resolving power: printed numbers that are
*derivable from other printed numbers*. A mis-transcription of any of the six
constants involved (0.75, 0.15, 0.22, 20, 10, 0.1) breaks at least one of
them."""))

cells.append(code(r'''# ---- (a) eq. 29 and eq. 30 are the same criterion in two notations -----------
# alpha^2 = k_b D_a / v^2 = ktau / Pe_L.  Eq. 30 with tau = L/v reads
# 1 > 20 ktau/Pe_L = 20 alpha^2, i.e. alpha < 1/sqrt(20).
ALPHA_STAR = 1.0 / np.sqrt(20.0)
print(f"(a) eq. 30 in alpha form: alpha < 1/sqrt(20) = {ALPHA_STAR:.4f}; "
      f"the paper prints 0.22 (eq. 29) -> agrees to {abs(ALPHA_STAR-0.22)/0.22*100:.1f} %")
print(f"    and eq. 28 (alpha < 1) is 'four and a half times' looser: sqrt(20) = "
      f"{np.sqrt(20):.3f}")
print("    [STRUCTURAL: this is algebra relating two printed constants; it tests the")
print("     transcription of 0.22 and 20 against each other, and nothing else]")

# ---- (b) eq. 16 vs eq. 17 for a first-order sphere ---------------------------
# omega = R r_p/(C_b k_c) = 3 * (eta k/(k_c a)) for a = 3/r_p, first order,
# READING Carberry's eta as bulk-referenced. His 0.1 is then omega < 0.3 -
# exactly 2x Mears' 0.15. Under Mears' own surface-referenced eta (p. 542 and
# the nomenclature) the same 0.1 is omega < 3*0.1/1.1 = 0.273, factor 1.8;
# neither paper says which reading Carberry meant, so Section 3 prints the
# deviation both ways.
print(f"\n(b) Carberry eq. 16 (0.1) in eq. 17's variable: omega < 0.3 = "
      f"{0.3/0.15:.0f} x Mears' 0.15 (bulk-referenced")
print(f"    reading; surface-referenced: omega < {3*0.1/1.1:.3f}, factor "
      f"{3*0.1/1.1/0.15:.1f})")
print("    [STRUCTURAL algebra; Section 3 measures what each boundary admits, "
      "under both readings]")

# ---- (c) the 0.75 the paper says was 'rounded to 1' --------------------------
# Small-Phi expansion of the sphere: eta = 1 - n Phi/15 + O(Phi^2), so a 5 %
# deviation sits at Phi = 15 x 0.05 / n = 0.75/n. Measured from the exact
# branch at its KINETIC end (t large -> y_c -> 1 -> Phi -> 0; the small-t end
# is deep diffusion, the opposite limit), extrapolated linearly in Phi:
slope = {}
for n in (0.5, 1.0, 2.0):
    (e1, P1), (e2, P2) = [(o[1], o[2]) for o in (shoot(power(n), t)
                                                 for t in (6.2, 6.9))]
    d1, d2 = (1 - e1) / P1, (1 - e2) / P2
    slope[n] = d2 + (d1 - d2) * (0.0 - P2) / (P1 - P2)   # extrapolate to Phi = 0
    print(f"(c) n = {n}: measured d(1-eta)/dPhi at Phi->0 = {slope[n]:.5f}  "
          f"(n/15 = {n/15:.5f})")
SLOPE_ERR = max(abs(slope[n] * 15 / n - 1) for n in slope)
assert SLOPE_ERR < 5e-3, SLOPE_ERR
print(f"    perturbation slope confirmed to {SLOPE_ERR:.1e} relative: 15 x 0.05 = 0.75 is")
print("    the unrounded constant of eqs. 2/3, exactly as the paper says it rounded away.")

# ---- (d) eq. 15's threshold of 10 is derivable from 0.75 and 0.15 ------------
# Film heat binds before intraparticle heat iff gamma*chi/0.15 > gamma*beta*Phi/C6
# with C6 the eq. 6 bound. chi/(beta*Phi) = lambda/(h r_p) = 2/Bi_p, so the
# crossover is Bi_p = 2*C6/0.15: with the UNROUNDED C6 = 0.75 that is exactly 10;
# with the printed C6 = 1 it would be 13.3.
BI_STAR_UNROUNDED, BI_STAR_ROUNDED = 2 * 0.75 / 0.15, 2 * 1.0 / 0.15
print(f"\n(d) eq. 15 crossover derived from eqs. 6 + 14: Bi_p = 2 x 0.75/0.15 = "
      f"{BI_STAR_UNROUNDED:.0f} (unrounded)")
print(f"    vs 2 x 1/0.15 = {BI_STAR_ROUNDED:.1f} (rounded). The paper prints 10 - "
      f"consistent only with")
print("    the unrounded 0.75, which pins BOTH thresholds at once.")
print("    [This check can fail: mis-reading 0.15 as 0.5, or 0.75 as 0.5, breaks it.]")

# ---- (e) symbolic reductions the paper asserts -------------------------------
n_s, gb_s, om_s, chi_s, g_s = sp.symbols("n gammabeta omega gammachi gamma", positive=True)
Cs, k_s = sp.symbols("C_s k", positive=True)
R_pl = k_s * Cs ** n_s
hudgins = sp.simplify(R_pl / (Cs * sp.diff(R_pl, Cs)))     # eq. 2 -> eq. 3
eq21_iso = (1 + sp.Rational(33, 100) * chi_s * 0) / (n_s * (1 + sp.Rational(33, 100) * n_s * om_s))
print(f"\n(e) eq. 2 bound for R = k C^n: R_s/(C_s R'_cs) = {hudgins}  ->  eq. 3's 1/n")
print(f"    eq. 21 at omega = chi = 0, beta_b = 0: bound -> "
      f"{sp.simplify(eq21_iso.subs(om_s, 0))}  ->  eq. 3/9's isothermal 1/n")
print(f"    eq. 13's content: at |gamma beta| = 0.05 n, eq. 9's bound differs from")
print(f"    eq. 3's by |gb/(n-gb)| = {0.05/0.95:.4f} - a 5 % shift, self-consistent.")
print("    [STRUCTURAL: symbolic identities; they test the transcription only]")'''))

cells.append(md(r"""### 2. Intraparticle, isothermal: what eqs. 1-3 admit

The audit: trace the exact branch for each order, locate the state where the
observable $\Phi$ equals each printed bound, and read the deviation there; then
locate the exact 5 % state. The pymrm solver is run *at the same states* and
must agree at its measured order — and for $n = 1$ there is a third,
closed-form route, so the headline number is computed three independent
ways."""))

cells.append(code(r'''ISO = {}
print("isothermal sphere: what each printed bound actually admits")
print(f"{'n':>5} {'Phi at 5%':>10} {'bound 1/n':>10} {'dev at bound':>13} "
      f"{'dev at 0.75/n':>14} {'dev at Weisz':>13}")
WEISZ = {1.0: 0.6, 2.0: 0.3}
for n in (0.5, 1.0, 2.0, 3.0):
    R = power(n)
    s5 = cross(R, lambda p, e, P, o, ys: e - 0.95, lo=1e-3)
    dB = 1 - cross(R, lambda p, e, P, o, ys, tgt=1 / n: P - tgt, lo=1e-3)[2]
    d75 = 1 - cross(R, lambda p, e, P, o, ys, tgt=0.75 / n: P - tgt, lo=1e-3)[2]
    dW = (1 - cross(R, lambda p, e, P, o, ys, tgt=WEISZ[n]: P - tgt, lo=1e-3)[2]
          if n in WEISZ else None)
    ISO[n] = dict(Phi5=s5[3], phi5=s5[1], dev_bound=dB, dev_075=d75, dev_weisz=dW)
    print(f"{n:5.1f} {s5[3]:10.4f} {1/n:10.4f} {dB*100:12.3f}% {d75*100:13.3f}% "
          + (f"{dW*100:12.3f}%" if dW else "          n/a"))

DEV_AT_1_N1 = ISO[1.0]["dev_bound"]
DEV_AT_075_N1 = ISO[1.0]["dev_075"]
PHI5_N1 = ISO[1.0]["Phi5"]
print(f"\nThe pattern is uniform: at the printed bound 1/n the deviation is "
      f"{min(v['dev_bound'] for v in ISO.values())*100:.1f}-"
      f"{max(v['dev_bound'] for v in ISO.values())*100:.1f} %,")
print(f"at the unrounded 0.75/n it is "
      f"{min(v['dev_075'] for v in ISO.values())*100:.1f}-"
      f"{max(v['dev_075'] for v in ISO.values())*100:.1f} %. The 'rounded to 1' step IS "
      f"the 5 % claim.")
print(f"Weisz (1957)'s linear-approximation values 0.6 and 0.3 admit only "
      f"{ISO[1.0]['dev_weisz']*100:.1f} % and {ISO[2.0]['dev_weisz']*100:.1f} % -")
print("conservative, as the paper says ('conservatively 0.3 ... and 6 for zero order').")'''))

cells.append(code(r'''# --- the headline number three independent ways -------------------------------
# Route 1 (above): shooting.  Route 2: the closed form.  Route 3: pymrm FV.
phiB = brentq(lambda p: p * p * eta_iso1(p) - 1.0, 1e-3, 10, xtol=1e-14)
DEV_AT_1_CLOSED = 1 - eta_iso1(phiB)
p800 = Pellet(800)
f_pm = lambda p: p800.Phi(p800.solve(p, power(1.0))[0]) - 1.0
phiB_pm = brentq(f_pm, 0.5, 2.0, xtol=1e-12)
y_pm, rn_pm = p800.solve(phiB_pm, power(1.0))
assert rn_pm < RN_TOL, rn_pm
DEV_AT_1_PYMRM = 1 - p800.eta(y_pm, phiB_pm)
print("deviation admitted by eq. 1 (Phi < 1, first order, sphere) - three routes:")
print(f"  shooting      : {DEV_AT_1_N1*100:.4f} %")
print(f"  closed form   : {DEV_AT_1_CLOSED*100:.4f} %")
print(f"  pymrm (n=800) : {DEV_AT_1_PYMRM*100:.4f} %   (Newton residual {rn_pm:.1e})")
HEADLINE_SPREAD = max(abs(DEV_AT_1_N1 - DEV_AT_1_CLOSED),
                      abs(DEV_AT_1_PYMRM - DEV_AT_1_CLOSED))
print(f"  worst spread  : {HEADLINE_SPREAD:.2e}  "
      f"(pymrm and the closed form share nothing; shooting shares only the rate law)")

# --- pymrm at the 5 % states, refined -----------------------------------------
pm_err = {}
for n_u in (200, 400, 800):
    p = Pellet(n_u)
    worst = 0.0
    for n, v in ISO.items():
        y, rn = p.solve(v["phi5"], power(n))
        assert rn < RN_TOL, (n, rn)
        worst = max(worst, abs(p.eta(y, v["phi5"]) - 0.95) / 0.95)
    pm_err[n_u] = worst
PM_ISO = pm_err[800]
PM_ISO_ORDER = float(np.log2(pm_err[400] / pm_err[800]))
print(f"\npymrm vs shooting at the four 5 % states: "
      + "  ".join(f"n_u={k}: {v:.2e}" for k, v in pm_err.items()))
print(f"observed order 400 -> 800: {PM_ISO_ORDER:.2f}")

# --- zero order: the one exact threshold in the family ------------------------
# y = 1 - (phi^2/6)(1 - u^2) exactly, eta = 1 identically until the centre hits
# zero at phi^2 = 6. So Weisz's '6 for zero order' is not an approximation - it
# is the exact dead-core onset, and below it the deviation is exactly zero.
p400 = Pellet(400)
y0_, rn0 = p400.solve(np.sqrt(5.7), lambda y: np.ones_like(np.asarray(y, float)))
assert rn0 < RN_TOL
y_exact = 1 - (5.7 / 6) * (1 - p400.u_c ** 2)
ZERO_PROFILE_ERR = float(np.max(np.abs(y0_ - y_exact)))
ZERO_ETA_ERR = abs(p400.eta(y0_, np.sqrt(5.7)) - 1.0)
print(f"\nzero order at phi^2 = 5.7 (pymrm): max|y - exact parabola| = "
      f"{ZERO_PROFILE_ERR:.2e}, |eta - 1| = {ZERO_ETA_ERR:.2e}")
print("Weisz's 6 is exact: dev = 0 for Phi < 6, a dead core beyond. The one member")
print("of the family with nothing to audit - and eq. 17 will echo it at the film.")'''))

cells.append(code(r'''# --- the extensions the paper claims for eq. 2: negative order and inhibition --
# Negative order via |n|: n = -1, bound 1/|n| = 1. The deviation is UPWARD
# (depleting the pellet speeds the rate up), so the two-sided 5 % point is
# eta = 1.05.
R_neg = power(-1.0)
s_neg = cross(R_neg, lambda p, e, P, o, ys: e - 1.05, lo=1e-4, hi=6.0)
NEG_PHI5, NEG_ETA_AT_1 = s_neg[3], None
s_neg1 = cross(R_neg, lambda p, e, P, o, ys: P - 1.0, lo=1e-4, hi=6.0)
NEG_ETA_AT_1 = s_neg1[2]
print(f"n = -1: |eta-1| = 5 % at Phi = {NEG_PHI5:.4f}, but the bound 1/|n| = 1 is not")
print(f"reached until eta = {NEG_ETA_AT_1:.4f} - the absolute-value extension admits "
      f"{(NEG_ETA_AT_1-1)*100:.1f} %.")

# Hudgins eq. 2 with an inhibited Langmuir-Hinshelwood rate, R = y(1+K)^2/(1+Ky)^2,
# K = 5 - the same rate law B1.4 used. R'(1) < 0, so eq. 2's bound is
# |R_s/(C_s R'_cs)| = (1+K)/(K-1) = 1.5.
K = 5.0
R_lh = lambda y: np.clip(y, 0.0, 1.0) * (1 + K) ** 2 / (1 + K * np.clip(y, 0.0, 1.0)) ** 2
LH_BOUND = (1 + K) / (K - 1)
s_lh_b = cross(R_lh, lambda p, e, P, o, ys, tgt=LH_BOUND: P - tgt, lo=1e-4)
s_lh_5 = cross(R_lh, lambda p, e, P, o, ys: abs(e - 1) - 0.05, lo=1e-4)
s_lh_1 = cross(R_lh, lambda p, e, P, o, ys: P - 1.0, lo=1e-4)
LH_ETA_AT_BOUND, LH_PHI5, LH_ETA_AT_1 = s_lh_b[2], s_lh_5[3], s_lh_1[2]
print(f"\ninhibited LH (K = 5): eq. 2's bound = {LH_BOUND:.2f}; eta there = "
      f"{LH_ETA_AT_BOUND:.4f}.")
print(f"The criterion PASSES a rate that is {(LH_ETA_AT_BOUND-1)*100:.1f} % high - "
      f"because its promise is")
print("one-sided (eta > 0.95), and for rate-inhibited kinetics the deviation is upward.")
print(f"The two-sided 5 % point is Phi = {LH_PHI5:.4f}, {LH_BOUND/LH_PHI5:.2f}x below "
      f"the bound.")
print(f"\ncross-page reconciliation: eta at Phi = 1 on this branch = {LH_ETA_AT_1:.4f};")
print("the published B1.4 prints 1.0483 for the same rate law - same number, computed")
print("from an independently written branch tracer.")'''))

cells.append(md(r"""### 3. Interphase: the one criterion that is exact, and the one that is
liberal by construction

Because the film lumps each resistance at the surface, the steady balances make
the surface state an *exact* function of the observables:
$1 - C_s/C_b = \omega/3$ and $(T_s-T_b)/T_b = \chi/3$. The deviations the
criteria bound are therefore closed forms —

$$\mathrm{dev}_{\rm mass} = 1 - \left(1 - \frac{\omega}{3}\right)^{\!n},
\qquad
\mathrm{dev}_{\rm heat} = \exp\!\left[\frac{\gamma_b\,\chi/3}{1+\chi/3}\right] - 1,$$

and the entire audit reduces to evaluating them at each printed boundary. For
$n = 1$ the mass form is *linear*: $\mathrm{dev} = \omega/3$ identically, so
eq. 17's boundary $\omega = 0.15$ admits **exactly** 5.000 % — the only
criterion in the paper whose stated tolerance is exact rather than approximate.
Zero order gives $\mathrm{dev} \equiv 0$, which is the paper's own remark that
"zero-order reactions are not affected by interphase transport", here exact.

The closed forms are cross-checked by dimensional solves that never use them:
physical parameters, a Brent solve of the nonlinear film balance for the
surface state, observables computed from the solved state."""))

cells.append(code(r'''# --- eq. 17 (mass) at its boundary, per order ---------------------------------
FILM_MASS = {n: 1 - (1 - 0.05 / n) ** n for n in (0.5, 1.0, 2.0, 3.0)}
print("eq. 17 boundary omega = 0.15/n:")
for n, d in FILM_MASS.items():
    tag = "  <- EXACT" if n == 1.0 else ""
    print(f"  n = {n}: dev = {d*100:.3f} %{tag}")
print(f"  n = 0 : dev = 0 exactly ('zero-order reactions are not affected')")
print(f"  eq. 18 is eq. 17 at n = 1 with C* k_LS for C_b k_c: its 0.15 admits "
      f"exactly 5.000 %.")
# Carberry's eq. 16 group is eta*k/(k_c*a), and the deviation it admits
# depends on which rate his eta references. Under Mears' own eta definition
# (rate at SURFACE conditions - p. 542 and the nomenclature) the group is
# R/(k_c a C_s) = (C_b - C_s)/C_s, so 0.1 pins C_s/C_b = 1/1.1 and the
# deviation is exactly 1/11; read bulk-referenced (R/(k_c a C_b) = omega/3)
# it is omega = 0.3 and exactly 10 %. Neither this paper nor a consulted
# original settles the reading, so both are printed - either supports the
# ~2x comparison against eq. 17.
CARBERRY_DEV = 0.3 / 3                       # bulk-referenced reading
CARBERRY_DEV_SURF = 0.1 / 1.1                # surface-referenced (Mears' eta)
print(f"\nCarberry's eq. 16 boundary, first order: dev = "
      f"{CARBERRY_DEV_SURF*100:.2f} % (surface-referenced eta, Mears'")
print(f"own definition; omega = {3*0.1/1.1:.3f}) or {CARBERRY_DEV*100:.1f} % "
      f"(bulk-referenced; omega = 0.3) - the reading is")
print("Carberry's to settle, and either way 'negligible' in 1961 meant about twice")
print("what the 1971 generalisation allows.")

# --- eq. 14 (heat) at its boundary, per gamma ---------------------------------
FILM_HEAT = {}
for g in (10.0, 20.0, 30.0, 40.0):
    u = 0.05 / g                              # chi/3 at the boundary
    FILM_HEAT[g] = (np.exp(0.05 / (1 + u)) - 1, 1 - np.exp(-0.05 / (1 - u)))
print("\neq. 14 boundary gamma_b |chi| = 0.15:")
for g, (dx, dn) in FILM_HEAT.items():
    print(f"  gamma_b = {g:4.0f}: exothermic dev = {dx*100:.3f} %   "
          f"endothermic dev = {dn*100:.3f} %")
HEAT_EXO_MAX = max(v[0] for v in FILM_HEAT.values())
HEAT_EXO_MIN = min(v[0] for v in FILM_HEAT.values())
print(f"The exothermic side exceeds 5 % for EVERY gamma - "
      f"{HEAT_EXO_MIN*100:.2f}-{HEAT_EXO_MAX*100:.2f} % - because the")
print("perturbation linearises exp(): e^0.05 - 1 = 5.13 %. Liberal by construction,")
print("but only by ~0.1 percentage point; the endothermic side is conservative by the")
print("same curvature.")'''))

cells.append(code(r'''# --- dimensional cross-checks that never touch the closed forms ---------------
# Mass, n = 2: r_p = 2 mm, k_c = 0.01 m/s, C_b = 50 mol/m3. Brent-solve the
# physical balance k_c a (C_b - C_s) = k C_s^2 for C_s; tune k so the OBSERVED
# omega sits exactly on the boundary; read the deviation from the solved state.
rp, kc, Cb, n_ = 2e-3, 0.01, 50.0, 2.0
a_ = 3.0 / rp
def film_state(k):
    Cs = brentq(lambda C: kc * a_ * (Cb - C) - k * C ** n_, 0.0, Cb, xtol=1e-15)
    Robs = k * Cs ** n_
    return Cs, Robs, Robs * rp / (Cb * kc)          # omega from observables
k_star = brentq(lambda k: film_state(k)[2] - 0.15 / n_, 1e-9, 1e3, xtol=1e-18)
Cs_, Robs_, om_ = film_state(k_star)
DEV_MASS_NUM = 1 - Robs_ / (k_star * Cb ** n_)
print(f"dimensional mass check (n = 2): omega = {om_:.6f}, dev = "
      f"{DEV_MASS_NUM*100:.4f} % vs closed form {FILM_MASS[2.0]*100:.4f} %")
FILM_MASS_XCHK = abs(DEV_MASS_NUM - FILM_MASS[2.0])

# Heat, first order: T_b = 500 K, gamma_b = 30, r_p = 2 mm, h = 100 W/m2K,
# -dH = 1e5 J/mol. Solve the film energy balance h a (T_s - T_b) = dH R_obs on
# its LOW branch, parametrised by T_s (monotone there; tuning the prefactor
# instead runs into the ignition fold, where the inner solve loses its sign
# change); the prefactor follows from the balance. Tune T_s so the OBSERVED
# gamma_b*chi sits exactly on the boundary; mass film switched off (k_c -> inf)
# to isolate the heat criterion, as eq. 14's derivation assumes.
Tb, gam_b, h_, dH = 500.0, 30.0, 100.0, 1e5
E_R = gam_b * Tb
def heat_state(Ts):
    Robs = h_ * a_ * (Ts - Tb) / dH          # energy balance, lumped at surface
    kinf = Robs / (np.exp(-E_R / Ts) * Cb)   # the prefactor that puts T_s there
    chi = dH * Robs * rp / (h_ * Tb)         # chi from the OBSERVED rate
    return kinf, Robs, chi
Ts_ = brentq(lambda T: gam_b * heat_state(T)[2] - 0.15, Tb * (1 + 1e-9), 1.5 * Tb,
             xtol=1e-10)
k_star_h, Robs_h, chi_ = heat_state(Ts_)
# assert the state really solves the balance (the check on the check)
assert abs(h_ * a_ * (Ts_ - Tb) - dH * k_star_h * np.exp(-E_R / Ts_) * Cb) < 1e-8
DEV_HEAT_NUM = Robs_h / (k_star_h * np.exp(-E_R / Tb) * Cb) - 1
print(f"dimensional heat check (gamma_b = 30): gamma chi = {gam_b*chi_:.6f}, dev = "
      f"{DEV_HEAT_NUM*100:.4f} % vs closed form {FILM_HEAT[30.0][0]*100:.4f} %")
FILM_HEAT_XCHK = abs(DEV_HEAT_NUM - FILM_HEAT[30.0][0])
print("Both routes share the film model itself (they must), but not one line of")
print("algebra: a dropped 3, a wrong sign, or a slipped exponent moves them apart -")
print("the break table does exactly that and watches them separate.")'''))

cells.append(code(r'''# --- multiplicity cannot hide in the film observables -------------------------
# With BOTH films active and strong exothermicity the forward problem ignites:
# one prefactor, three steady states. The point measured here: unlike the
# distributed pellet (B1.4's published fold, where one observed Phi fits three
# eta), the film observables DETERMINE the deviation - dev is one closed-form,
# monotone function of (omega, gamma chi), whatever branch the state sits on.
kc2, h2, dH2, Cb2 = 0.004, 20.0, 3e5, 100.0
def both_films(kinf):
    """All steady states (C_s, T_s) for one prefactor, by dense scan + refine."""
    def F(R):
        Cs = Cb2 - R * rp / (3 * kc2)
        Ts = Tb + dH2 * R * rp / (3 * h2 * Tb) * Tb / 1.0
        if Cs <= 0:
            return np.nan
        return kinf * np.exp(-E_R / Ts) * Cs - R
    Rmax = 3 * kc2 * Cb2 / rp * (1 - 1e-12)
    Rs = np.linspace(1e-8, Rmax, 4001)
    vals = np.array([F(r) for r in Rs])
    roots = []
    for i in range(len(Rs) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            roots.append(brentq(F, Rs[i], Rs[i + 1], xtol=1e-12))
    return roots

KINF_MULTI = None
for kf in np.logspace(4, 9, 60):
    if len(both_films(kf)) == 3:
        KINF_MULTI = kf
        break
roots = both_films(KINF_MULTI)
print(f"one prefactor ({KINF_MULTI:.3e}), {len(roots)} steady states (S-curve):")
print(f"{'state':>8} {'omega':>9} {'gamma chi':>10} {'dev (model)':>12} "
      f"{'dev (closed form)':>18} {'eq17':>6} {'eq14':>6}")
COLLAPSE = 0.0
for lbl, R in zip(("low", "middle", "ignited"), roots):
    om = R * rp / (Cb2 * kc2)
    chi = dH2 * R * rp / (h2 * Tb)
    dev_model = R / (KINF_MULTI * np.exp(-E_R / Tb) * Cb2) - 1
    u = chi / 3
    dev_cf = (1 - om / 3) * np.exp(gam_b * u / (1 + u)) - 1
    COLLAPSE = max(COLLAPSE, abs(dev_model - dev_cf) / max(1.0, abs(dev_cf)))
    print(f"{lbl:>8} {om:9.4f} {gam_b*chi:10.3f} {dev_model:12.4g} {dev_cf:18.4g} "
          f"{'pass' if om < 0.15 else 'FAIL':>6} {'pass' if gam_b*chi < 0.15 else 'FAIL':>6}")
print(f"\nclosed-form collapse over all branches (relative): {COLLAPSE:.2e}")
print("[STRUCTURAL as a check - dev and the observables are computed from the same")
print(" solved state, so the collapse is an identity of the film model. What is NOT")
print(" structural is the finding it illustrates: every state beyond the low branch")
print(" carries observables that FAIL the criteria, so film multiplicity cannot")
print(" produce the false-pass that B1.4 measured for the pellet, where the fold sits")
print(" INSIDE the safe band. Lumped observables determine the state; distributed")
print(" ones need not.]")'''))

cells.append(md(r"""### 4. Intraparticle, nonisothermal: eq. 9's window of no power

Eq. 9 divides the isothermal bound by $|n - \gamma\beta|$. Near the
compensation point $n = \gamma\beta$ the bound diverges — but the exact 5 %
threshold stays finite, so there is a *window* of $\gamma\beta$ in which the
criterion certifies states it should reject. The window is measured against
the criterion's own honest constant: where even the **unrounded**
$0.75/|n-\gamma\beta|$ admits more than 5 %. Both bounds are measured here,
along with the paper's three printed statements about the compensation case:
Petersen's 13, the "80 %" concentration drop, and Weisz-Hicks eq. 12 as the
comparison.

One structural fact governs how every number here must be read: above
$\gamma\beta \approx 1.45$ the branch's $\eta$ overshoots past 1.05, so the
$\pm 5\,\%$ band has **three** crossings, and which one a number refers to
matters. The page's `cross` convention is the *deep-side* entry — the
largest-$\Phi$ boundary state — stated in its docstring; the cell after the
window prints the full crossing structure, and measures the second, mild
over-certification regime that opens past $\gamma\beta \approx 1.95$,
*outside* the window, where the unrounded bound overtakes the **kinetic**-side
crossing. So "outside the window the 0.75 constant restores the tolerance" is
true only up to that point, and the page says by how much beyond it.

Anderson's eq. 6 bounds a *different* deviation — observed rate vs the rate at
uniform temperature — so it is audited against exactly that:
$\eta/\eta_{\rm iso}(\phi) - 1$ at the state where its own observable
$\gamma|\beta|\Phi$ crosses 1."""))

cells.append(code(r'''GAMMA = 30.0
GBS = (-2.0, -1.0, -0.5, -0.2, 0.2, 0.5, 0.8, 1.0, 1.2, 2.0)
MAP9 = {}
print(f"eq. 9 audit, n = 1, gamma = {GAMMA:.0f} (beta = gb/gamma, inside the paper's "
      f"-0.1..0.1)")
print(f"{'g*b':>6} {'Phi at 5%':>10} {'bound eq9':>10} {'dev at eq9':>11} "
      f"{'bound eq12':>11} {'dev at eq12':>12}")
for gb in GBS:
    R = wh(1.0, gb, GAMMA)
    s5 = cross(R, lambda p, e, P, o, ys: abs(e - 1) - 0.05, lo=1e-4, hi=22.0,
               n_scan=110)
    b9 = 1 / abs(1 - gb) if abs(1 - gb) > 1e-12 else np.inf
    d9 = (abs(1 - cross(R, lambda p, e, P, o, ys, t=b9: P - t, lo=1e-4, hi=22.0,
                        n_scan=110)[2]) if np.isfinite(b9) else np.nan)
    beta = gb / GAMMA
    b12 = float(np.exp(-gb / (1 + beta)))
    d12 = abs(1 - cross(R, lambda p, e, P, o, ys, t=b12: P - t, lo=1e-4, hi=22.0,
                        n_scan=110)[2])
    MAP9[gb] = dict(Phi5=s5[3], bound9=b9, dev9=d9, bound12=b12, dev12=d12)
    print(f"{gb:6.1f} {s5[3]:10.4f} {b9:10.4f} "
          + (f"{d9*100:10.2f}%" if np.isfinite(d9) else "       inf")
          + f" {b12:11.4f} {d12*100:11.2f}%")

DEV9_GB05 = MAP9[0.5]["dev9"]
DEV9_GB08 = MAP9[0.8]["dev9"]
DEV12_ENDO = MAP9[-1.0]["dev12"]
PHI5_COMP = MAP9[1.0]["Phi5"]
endo_devs = [MAP9[g]["dev9"] for g in GBS if g < 0]
print(f"\nOn the endothermic side (gb < 0) the bound admits "
      f"{min(endo_devs)*100:.1f}-{max(endo_devs)*100:.1f} % - the familiar 0.75 "
      f"rounding.")
print(f"Approaching gb = n from below it climbs: {MAP9[0.2]['dev9']*100:.1f} % at "
      f"gb = 0.2, {DEV9_GB05*100:.1f} % at 0.5,")
print(f"{DEV9_GB08*100:.0f} % at 0.8, anything at 1. Past compensation it falls back "
      f"through {MAP9[1.2]['dev9']*100:.1f} % (gb = 1.2,")
print(f"whose single band crossing sits above even the printed bound) to "
      f"{MAP9[2.0]['dev9']*100:.1f} % (gb = 2).")
print("Convention: 'Phi at 5 %' is the DEEP-side band entry (the largest-Phi")
print("boundary state - see cross's docstring). For gb <= 1.2 in this table the")
print("crossing is unique, so the label is redundant there; for gb = 1.5 and 2.0 the")
print("band has three crossings and the next cells print all of them. 'dev at eq9'")
print("and 'dev at eq12' need no such label: Phi is monotone along these branches")
print("(no fold until gb ~ 11, B1.4's measured onset), so Phi = bound is one state.")
print(f"Weisz-Hicks eq. 12, exothermic: conservative everywhere here (its bound sits")
print(f"below the 5 % point), exactly the safe direction the paper claims. Endothermic,")
print(f"gb = -1: its bound is {MAP9[-1.0]['bound12']:.2f} and the deviation there is "
      f"{DEV12_ENDO*100:.0f} % - the failure the")
print("paper states ('fails to predict negative deviations'), now with a number.")'''))

cells.append(code(r'''# --- the window of no power, measured ------------------------------------------
# The PRINTED bound is liberal everywhere (that is Section 2's rounding), so
# "void" must be measured against the criterion's own honest constant: the
# window around gb = n where even Phi < 0.75/|n - gb| - the UNROUNDED bound -
# certifies deviations beyond 5 %. Inside it no constant fix works, because
# the bound's functional form diverges at compensation while the exact
# threshold stays finite. Outside it, restoring 0.75 restores the stated
# tolerance - but ONLY up to gb ~ 1.95: the next cell measures the second,
# mild regime beyond that, where the band has three crossings and the bound
# overtakes the kinetic-side one.
def window_edge(sign, n=1.0, c=0.75):
    def excess(gb):
        R = wh(n, gb, GAMMA)
        s5 = cross(R, lambda p, e, P, o, ys: abs(e - 1) - 0.05, lo=1e-4, hi=22.0,
                   n_scan=70)
        return c / abs(n - gb) - s5[3]
    gbs = n + sign * np.linspace(0.05, 0.95, 10)
    vals = [excess(g) for g in gbs]
    for i in range(len(gbs) - 1):
        if vals[i] * vals[i + 1] < 0:
            return brentq(excess, gbs[i], gbs[i + 1], xtol=2e-4)
    raise RuntimeError(f"no window edge inside the scanned range, side {sign:+d}")

W_LO, W_HI = window_edge(-1), window_edge(+1)
# what the unrounded bound admits at a point inside the window, for scale
R_in = wh(1.0, 1.0, GAMMA)
DEV_075_INSIDE = abs(1 - cross(R_in, lambda p, e, P, o, ys: P - 0.75 / 0.05,
                               lo=1e-4, hi=22.0, n_scan=110)[2])
print(f"even the UNROUNDED bound 0.75/|n - gb| over-certifies for gb in "
      f"({W_LO:.3f}, {W_HI:.3f})")
print(f"around n = 1 (gamma = {GAMMA:.0f}) - width {W_HI-W_LO:.3f} in gamma*beta, and "
      f"markedly ASYMMETRIC:")
print(f"{1-W_LO:.2f} below n against {W_HI-1:.2f} above. That matches the previous "
      f"cell's table: the")
print("liberality climbs steadily approaching n from below but drops back quickly")
print("past it, where the exact threshold stays high while the bound collapses like")
print("1/(gb - 1), turning even the printed bound conservative (gb = 1.2's single")
print("band crossing sits above it).")
print(f"At its edges the unrounded bound admits exactly 5 % - and at both edges the")
print(f"band crossing is unique (eta never reaches 1.05 below gb = 1.45), so that")
print(f"statement carries no crossing-selection ambiguity. At gb = 0.95 "
      f"(|n-gb| = 0.05,")
print(f"inside) it admits {DEV_075_INSIDE*100:.0f} % - the divergence is fast. That "
      f"is the measured")
print("content of the paper's unquantified 'close to or equal to zero requires")
print("further consideration'. The printed bound, 1/0.75 looser again, admits its")
print(f"{min(endo_devs)*100:.1f}-{max(endo_devs)*100:.1f} % on the endothermic side "
      f"(previous cell) and MORE than this inside")
print("the window.")

# --- the compensation case: Petersen's 13, and the 80 % statement -------------
print(f"\nAt gb = n = 1 exactly (gamma = {GAMMA:.0f}): the exact 5 % threshold is "
      f"Phi = {PHI5_COMP:.3f}.")
pet = {}
for g in (20.0, 30.0, 40.0):
    R = wh(1.0, 1.0, g)
    s5 = cross(R, lambda p, e, P, o, ys: abs(e - 1) - 0.05, lo=1e-3, hi=22.0)
    d13 = abs(1 - cross(R, lambda p, e, P, o, ys: P - 13.0, lo=1e-3, hi=22.0)[2])
    pet[g] = (s5[3], -np.expm1(-s5[0]), d13)
    print(f"  gamma = {g:.0f}: Phi at 5 % = {pet[g][0]:.3f}, centre y there = "
          f"{pet[g][1]:.3f}, dev at Petersen's 13 = {pet[g][2]*100:.1f} %")
PET13_DEV = pet[30.0][2]
YC_COMP = pet[30.0][1]
print(f"\nPetersen's 13 is not a 5 % criterion and never claimed to be: at Phi = 13 the")
print(f"deviation is {min(v[2] for v in pet.values())*100:.1f}-"
      f"{max(v[2] for v in pet.values())*100:.1f} % - it marks the transition region, "
      f"as the paper itself explains")
print("in the axial context. And the '80 %' statement is confirmed conservatively: the")
print(f"rate stays within 5 % until the centre concentration has fallen not to 0.8 but")
print(f"to {YC_COMP:.2f} of the surface value.")'''))

cells.append(code(r'''# --- the crossing structure above the window, and the second regime -----------
# Above the window eta OVERSHOOTS: once its peak exceeds 1.05 the +-5 % band
# has three crossings, and 'the 5 % threshold' stops being one number. Printed
# in full because the page's own convention (deep-side entry) must be
# distinguishable from the crossing a slowly worsening experiment meets first
# (the kinetic-side one) - the two differ by a factor ~15 at gb = 2.
CR = {}
for gb in (1.2, 1.5, 2.0):
    R = wh(1.0, gb, GAMMA)
    allc = cross(R, lambda p, e, P, o, ys: abs(e - 1) - 0.05, lo=1e-4, hi=22.0,
                 n_scan=200, which="all")
    pk = minimize_scalar(lambda t: -shoot(R, t)[1], bounds=(1e-4, 22.0),
                         method="bounded", options={"xatol": 1e-8})
    o_pk = shoot(R, pk.x)
    CR[gb] = (sorted(c[3] for c in allc), o_pk[1], o_pk[2])
    print(f"gb = {gb}: |eta-1| = 5 % crossings at Phi = "
          + ", ".join(f"{p:.4g}" for p in CR[gb][0])
          + f"; eta peaks at {CR[gb][1]:.4f} near Phi = {CR[gb][2]:.3g}")
MC_ONSET = brentq(lambda g: -minimize_scalar(
    lambda t: -shoot(wh(1.0, g, GAMMA), t)[1], bounds=(1e-4, 22.0),
    method="bounded", options={"xatol": 1e-6}).fun - 1.05, 1.2, 1.6, xtol=5e-4)
print(f"the overshoot reaches 1.05 - three crossings from there on - at gb = "
      f"{MC_ONSET:.3f}")

# The consequence for the window claim: past the point where the unrounded
# bound overtakes the KINETIC-side crossing, its certified band re-admits
# more than 5 % - mildly, from curvature, nothing like the window's
# divergence. Measured up to gb = 3, i.e. beta = 0.1 at gamma = 30, the edge
# of the paper's own beta range. (Phi is monotone along these branches, so
# 'dev at the bound' is one root-found state, not a grid max.)
def phi_kinetic(gb):
    R = wh(1.0, gb, GAMMA)
    allc = cross(R, lambda p, e, P, o, ys: abs(e - 1) - 0.05, lo=1e-4, hi=22.0,
                 n_scan=140, which="all")
    return min(c[3] for c in allc)

SECOND_ONSET = brentq(lambda g: 0.75 / abs(1 - g) - phi_kinetic(g), 1.6, 2.0,
                      xtol=5e-4)
OUT2 = {}
for gb in (2.0, 2.5, 3.0):
    R = wh(1.0, gb, GAMMA)
    b075 = 0.75 / abs(1 - gb)
    o = cross(R, lambda p, e, P, o_, ys, t=b075: P - t, lo=1e-4, hi=22.0,
              n_scan=140)
    OUT2[gb] = abs(o[2] - 1)
    print(f"gb = {gb}: unrounded bound {b075:.4f} vs kinetic-side crossing "
          f"{phi_kinetic(gb):.4f} -> dev at the bound {OUT2[gb]*100:.2f} %")

print(f"\nThe honest scope of the window claim, then: inside ({W_LO:.3f}, "
      f"{W_HI:.3f}) the criterion is")
print("beyond repair - the divergence certifies arbitrarily large deviations.")
print(f"From the upper edge to gb = {SECOND_ONSET:.2f} the unrounded 0.75 restores "
      f"the stated")
print(f"tolerance. Beyond that, up to the paper's own beta range edge (gb = 3 at")
print(f"gamma = 30), the bound re-admits {min(OUT2.values())*100:.2f}-"
      f"{max(OUT2.values())*100:.2f} % - the same mild curvature class")
print("as eq. 14's 5.10-5.12 %, reported the same way, and nothing like the window.")
print("On the endothermic side the crossing stays unique and the tolerance holds")
print("everywhere below the window's lower edge.")'''))

cells.append(code(r'''# --- Anderson's eq. 6, against the deviation it actually bounds ----------------
print("eq. 6 audit: dev_T = eta/eta_iso(phi) - 1 at the state where gamma|beta|Phi = 1")
print(f"{'gamma':>6} {'beta':>8} {'g*b':>6} | {'phi':>7} {'Phi':>8} {'dev_T':>9}")
AND6 = {}
for gamma, beta in ((20.0, 0.05), (30.0, 0.02), (40.0, 0.01), (30.0, 0.005),
                    (30.0, -0.02), (20.0, -0.05)):
    gb = gamma * beta
    R = wh(1.0, gb, gamma)
    sA = cross(R, lambda p, e, P, o, ys, gA=abs(gb): gA * P - 1.0,
               lo=1e-5, hi=22.0)
    devT = sA[2] / eta_iso1(sA[1]) - 1
    AND6[(gamma, beta)] = devT
    print(f"{gamma:6.0f} {beta:8.3f} {gb:6.1f} | {sA[1]:7.3f} {sA[3]:8.3f} "
          f"{devT*100:+8.3f}%")
AND6_MAX = max(abs(v) for v in AND6.values())
AND6_AT_GB1 = AND6[(20.0, 0.05)]
print(f"\nworst |dev_T| at the boundary: {AND6_MAX*100:.2f} % (strong heat, kinetic "
      f"regime); mild heat in")
print("deep diffusion is conservative. The claim 'valid whether diffusional limitations")
print("exist or not' holds to within ~1.1 percentage points - the same rounded-0.75")
print("looseness, softened because the criterion is usually met deep in the branch.")

# --- pymrm on the nonisothermal states, seeded from the branch ----------------
worst = {}
for n_u in (400, 800):
    p = Pellet(n_u)
    w = 0.0
    for gb in (-1.0, 0.5, 1.0):
        R = wh(1.0, gb, GAMMA)
        s5 = cross(R, lambda p_, e, P, o, ys: abs(e - 1) - 0.05, lo=1e-4, hi=22.0)
        o = shoot(R, s5[0], dense=True)
        phi, sol = o[0], o[5]
        y0 = np.clip(sol.sol(np.clip(p.u_c * phi, sol.t[0], phi))[0], 0.0, 1.0)
        y, rn = p.solve(phi, R, y_init=y0.reshape(-1, 1))
        assert rn < RN_TOL, (gb, rn)
        w = max(w, abs(p.eta(y, phi) - s5[2]) / abs(s5[2]))
    worst[n_u] = w
PM_WH = worst[800]
PM_WH_ORDER = float(np.log2(worst[400] / worst[800]))
print(f"\npymrm vs shooting at three nonisothermal 5 % states: n_u=400: "
      f"{worst[400]:.2e}, n_u=800: {worst[800]:.2e}, order {PM_WH_ORDER:.2f}")'''))

cells.append(code(r'''# --- where B1.4's fold finding lands on eq. 9 ---------------------------------
# B1.4 (published) measured that the observable Phi becomes MULTIVALUED for
# beta*gamma >~ 11-23, with every fold turning point BELOW Phi = 1 - so eq. 1
# certifies ignited states with enormous eta. Does the same fold defeat eq. 9?
# Trace one folded case (beta = 0.6, gamma = 40: B1.4's worst region) and score
# the SAME branch against both criteria's certified bands.
#
# The worst CERTIFIED deviation of a band Phi < b is a supremum whose sup sits
# AT the band edge, so a max over branch samples is biased low by exactly the
# grid spacing. The first staged version of this cell shipped that defect -
# 5.7 % from a 400-point grid max, ~20 % low; adversarial verification caught
# it, and it is the very grid-limited-extremum class this cell was recast to
# escape. The estimator now root-finds every band-edge state Phi = b on the
# branch and takes the max over edge states and the interior grid; the grid
# max at two densities is printed as the convergence evidence, and the break
# table re-injects the coarse-grid estimator as row (xi).
gb_f, gam_f = 24.0, 40.0
R_f = wh(1.0, gb_f, gam_f)

def trace_fold(n_pts):
    ts = np.logspace(-7, np.log10(24.0), n_pts)
    pairs = [(t, o) for t, o in ((t, shoot(R_f, t)) for t in ts) if o is not None]
    return (np.array([p[0] for p in pairs]),
            np.array([p[1][2] for p in pairs]),      # Phi
            np.array([p[1][1] for p in pairs]))      # eta

def band_worst(b, ts, Phi, eta):
    """sup |eta - 1| over the certified band Phi < b: every band-edge crossing
    Phi = b brentq-refined on the branch, plus the interior grid max.
    Returns (grid_max, worst, t_at_worst_edge_or_None)."""
    dev = np.abs(eta - 1)
    grid_max = float(dev[Phi < b].max())
    f = lambda t: shoot(R_f, t)[2] - b
    edge, d = [], Phi - b
    for i in range(len(ts) - 1):
        if d[i] * d[i + 1] < 0:
            tt = brentq(f, ts[i], ts[i + 1], xtol=1e-13, rtol=8.9e-16)
            edge.append((abs(shoot(R_f, tt)[1] - 1), tt))
    worst_edge = max(edge) if edge else (grid_max, None)
    return (grid_max, float(max(worst_edge[0], grid_max)),
            worst_edge[1] if worst_edge[0] >= grid_max else None)

ts_f, Phi_f, eta_f = trace_fold(400)
ts_d, Phi_d, eta_d = trace_fold(1500)
dev_f = np.abs(eta_f - 1)
run = np.maximum.accumulate(Phi_f)
FOLD_DEPTH = float(((run - Phi_f) / run).max())
b9_f = 1 / abs(1 - gb_f)
GM9_400, WORST_DEV_CERT_EQ9, T9_EDGE = band_worst(b9_f, ts_f, Phi_f, eta_f)
GM9_1500, w9_chk, _ = band_worst(b9_f, ts_d, Phi_d, eta_d)
GM1_400, WORST_DEV_CERT_EQ1, _ = band_worst(1.0, ts_f, Phi_f, eta_f)
GM1_1500, w1_chk, _ = band_worst(1.0, ts_d, Phi_d, eta_d)
YC9_EDGE = float(-np.expm1(-T9_EDGE))          # centre y at the eq. 9 edge state
# the edge estimate must be grid-independent, and the grid maxima must
# converge to it from below - a check that fails if the estimator is a mirage
assert abs(w9_chk - WORST_DEV_CERT_EQ9) < 1e-6, (w9_chk, WORST_DEV_CERT_EQ9)
assert GM9_400 <= GM9_1500 <= WORST_DEV_CERT_EQ9 + 1e-12
assert abs(w1_chk - WORST_DEV_CERT_EQ1) < 1e-3 * WORST_DEV_CERT_EQ1
assert GM1_400 <= GM1_1500 <= WORST_DEV_CERT_EQ1 * (1 + 1e-12)
# the ignited branch's lowest reach is a smooth turning point - refine it too
mask_i = dev_f > 1.0                                   # far-ignited: eta > 2
i0 = int(np.argmin(np.where(mask_i, Phi_f, np.inf)))
res_m = minimize_scalar(lambda t: shoot(R_f, t)[2],
                        bounds=(ts_f[i0 - 1], ts_f[i0 + 1]),
                        method="bounded", options={"xatol": 1e-10})
MIN_PHI_IGNITED = float(res_m.fun)

print(f"gb = {gb_f:.0f} (beta = {gb_f/gam_f}, gamma = {gam_f:.0f}): fold depth in Phi "
      f"= {FOLD_DEPTH:.3f} (multivalued, as B1.4 found)")
print(f"eq. 1's band (Phi < 1)          : worst certified deviation = "
      f"{WORST_DEV_CERT_EQ1:.0f} - i.e. eta ~ {WORST_DEV_CERT_EQ1+1:.0f},")
print("                                  the fold-class false pass B1.4 measured")
print(f"eq. 9's band (Phi < 1/{abs(1-gb_f):.0f} = {b9_f:.3f}): worst certified "
      f"deviation = {WORST_DEV_CERT_EQ9*100:.2f} %")
print(f"  estimator convergence: grid max {GM9_400*100:.2f} % (400 pts) -> "
      f"{GM9_1500*100:.2f} % (1500 pts) ->")
print(f"  {WORST_DEV_CERT_EQ9*100:.2f} % at the root-found band edge, identical on "
      f"both grids. The sup sits AT")
print(f"  Phi = 1/23 on the low branch (eta there = "
      f"{1+WORST_DEV_CERT_EQ9:.4f}, y_c = {YC9_EDGE:.3f} - near-kinetic, strongly")
print(f"  heated); eq. 1's edge behaves the same way "
      f"({GM1_400:.0f} -> {GM1_1500:.0f} -> {WORST_DEV_CERT_EQ1:.0f}).")
print(f"ignited branch reaches down to Phi = {MIN_PHI_IGNITED:.4f} = "
      f"{MIN_PHI_IGNITED/b9_f:.2f}x eq. 9's bound (refined turning point)")
print("The fold that defeats eq. 1 does not defeat eq. 9 - not because eq. 9 sees the")
print("fold, but because at fold severities |n - gamma beta| is large and its bound")
print("has already shrunk below the ignited branch's reach. What its band still")
print(f"admits is {WORST_DEV_CERT_EQ9*100:.1f} % at its own edge - the rounding-class "
      f"liberality again, a shade")
print("above eq. 3's 6.3-6.8 % because gb = 24 heats even a near-kinetic pellet -")
print("nothing ignited. Its weakness is the compensation window; its strength is")
print("exactly where eq. 1's weakness lives.")
print("(Convention note: this cell scores every traced state, folded ones included,")
print(" and its suprema are root-found band-edge states. Eq. 9 numbers elsewhere on")
print(" this page are deep-side band crossings, stated where they occur.)")'''))

cells.append(md(r"""### 5. The combined criterion (eq. 21), and what stacking tolerances costs

Eq. 21 repeats eq. 9 with bulk observables and two film-correction factors.
Audited here in its isothermal reduction ($\chi = 0$, $\beta_b = 0$):
$\Phi_b < 1/\bigl(n(1+0.33\,n\omega)\bigr)$, with $\omega = \Phi_b/\mathrm{Bi}_m$
an observable. The pellet is re-solved with the film boundary condition and the
5 % state located per Biot number; **the thermal part of eqs. 21-22 is not
audited** — it couples the film to the pellet interior and its honest audit
needs the two-field problem, which is out of scope here and said so.

The second question is the one an experimenter actually faces: after passing
*every* pellet-scale criterion in the paper at once, how wrong can the rate
still be? Each criterion allows its own 5 %, and deviations compound."""))

cells.append(code(r'''print("eq. 21 (isothermal part), n = 1 and 2: exact 5 % state vs the bound")
print(f"{'n':>3} {'Bi_m':>7} {'Phi_b at 5%':>12} {'omega there':>12} "
      f"{'bound(omega)':>13} {'ratio':>7}")
EQ21 = {}
for n in (1.0, 2.0):
    for Bim in (100.0, 10.0, 3.0, 1.0, 0.3):
        R = power(n)
        s5 = cross(R, lambda p, e, P, o, ys: e - 0.95, lo=1e-3, hi=18.0, bim=Bim)
        Phi_b, om = s5[3], s5[4]
        bound = 1 / (n * (1 + 0.33 * n * om))
        EQ21[(n, Bim)] = (Phi_b, om, bound, Phi_b / bound)
        print(f"{n:3.0f} {Bim:7.1f} {Phi_b:12.4f} {om:12.4f} {bound:13.4f} "
              f"{Phi_b/bound:7.3f}")
RATIO21_HI = EQ21[(1.0, 100.0)][3]
RATIO21_LO = EQ21[(1.0, 0.3)][3]
print(f"\nAt Bi_m = 100 the ratio is {RATIO21_HI:.2f} - Section 2's rounding again "
      f"(0.767 shifted by the")
print(f"bulk normalisation). As Bi_m falls the ratio collapses to {RATIO21_LO:.3f}: "
      f"the criterion's")
print("left side Phi_b vanishes with the pellet gradient while the deviation moves")
print("into the film, which Phi_b cannot see. Eq. 21 alone is structurally unable to")
print("detect film-dominated limitation - it must be applied WITH eq. 17, which the")
print("paper nowhere says.")

# --- and if you do apply everything at once? ----------------------------------
# The worst deviation among states passing eqs. 1, 17 and 21 at once is a
# supremum over a pass SET: its sup sits on the set's boundary (where the
# binding criterion crosses zero), at an INTERIOR optimum in Bi_m - so a
# coarse Bi_m x t grid is biased low by both spacings. The first staged
# version shipped 9.0 % from exactly that defect (a {0.3..100} Bi_m sweep,
# grid-only in t); adversarial verification caught it, and the break table
# keeps the defective estimator as row (xi). Here: per Bi_m, every crossing
# of the most-binding margin is brentq-refined and the sup taken over edge
# states and the passing grid; then the worst case is maximised over Bi_m,
# whose optimum is interior and smooth.
def joint_worst(Bim, n_scan=110, detail=False):
    state = lambda t: shoot(power(1.0), t, bim=Bim)
    def margin(o):                          # > 0 = some criterion fails
        phi, eta_b, Phi_b, om, ys = o
        return max(om - 0.15,               # eq. 17
                   Phi_b - 1 / (1 + 0.33 * om),   # eq. 21, isothermal part
                   Phi_b / ys - 1.0)        # eq. 1 on the surface observable
    best, arg, prev = 0.0, None, None
    for t in np.logspace(-4, 1.25, n_scan):
        o = state(t)
        if o is None:
            prev = None
            continue
        g = margin(o)
        if g < 0 and 1 - o[1] > best:
            best, arg = 1 - o[1], o
        if prev is not None and prev[1] * g < 0:   # pass-set edge: root-find it
            oo = state(brentq(lambda x: margin(state(x)), prev[0], t, xtol=1e-12))
            if 1 - oo[1] > best:
                best, arg = 1 - oo[1], oo
        prev = (t, g)
    return (best, arg) if detail else best

JW_SWEEP = {Bim: joint_worst(Bim) for Bim in (0.3, 1.0, 3.0, 5.0, 7.0, 10.0,
                                              30.0, 100.0)}
print("\nworst first-order deviation passing eqs. 1, 17 AND 21 at once, per Bi_m:")
print("  " + "  ".join(f"Bi_m={k:g}: {v*100:.2f}%" for k, v in JW_SWEEP.items()))
res_j = minimize_scalar(lambda lb: -joint_worst(10 ** lb),
                        bounds=(np.log10(3.0), np.log10(15.0)),
                        method="bounded", options={"xatol": 2e-3})
JOINT_BIM = float(10 ** res_j.x)
JOINT_WORST, o_j = joint_worst(JOINT_BIM, detail=True)
phi_j, eta_j, Phi_j, om_j, ys_j = o_j
print(f"\nthe maximiser is INTERIOR: Bi_m = {JOINT_BIM:.2f}, worst dev = "
      f"{JOINT_WORST*100:.2f} %")
print(f"  at the pass-set edge state Phi_b = {Phi_j:.4f}, omega = {om_j:.4f}: "
      f"eq. 17 binds exactly")
print(f"  (omega - 0.15 = {om_j-0.15:+.1e}) and eq. 1-on-the-surface is within "
      f"{abs(Phi_j/ys_j-1):.1e} of binding -")
print("  the worst case is two resistances at their boundaries at once, which is")
print("  exactly why the tolerances stack.")
JW_DENSE = joint_worst(JOINT_BIM, n_scan=300)
print(f"  scan-density check: n_scan 110 -> 300 moves the worst by "
      f"{abs(JW_DENSE-JOINT_WORST):.1e} (edge is root-found,")
print("  so the grid only brackets); the Bi_m optimum is a smooth interior peak")
print("  (sweep above).")
assert abs(JW_DENSE - JOINT_WORST) < 1e-6
# pymrm at the worst state - the new headline gets its own cross-route check
p_j = Pellet(800, bim=JOINT_BIM)
y_j, rn_j = p_j.solve(phi_j, power(1.0))
assert rn_j < RN_TOL, rn_j
JOINT_PM_GAP = abs((1 - p_j.eta(y_j, phi_j)) - JOINT_WORST)
assert JOINT_PM_GAP < 1e-5, JOINT_PM_GAP
print(f"  pymrm (film bc, n_u = 800) at this state: dev agrees with shooting to "
      f"{JOINT_PM_GAP:.1e}")
print(f"\nTolerances stack: two resistances, each individually inside its own 5 %")
print(f"criterion, compound to {JOINT_WORST*100:.1f} % - more than double the "
      f"promise. None of the")
print("criteria says so; the perturbation derivations treat each resistance alone.")

# --- pymrm with the film (Robin) boundary condition, refined ------------------
Bim_chk = 3.0
s5 = cross(power(1.0), lambda p, e, P, o, ys: e - 0.95, lo=1e-3, hi=18.0, bim=Bim_chk)
rob_err = {}
for n_u in (200, 400, 800):
    p = Pellet(n_u, bim=Bim_chk)
    y, rn = p.solve(s5[1], power(1.0))
    assert rn < RN_TOL, rn
    rob_err[n_u] = abs(p.eta(y, s5[1]) - s5[2]) / s5[2]
PM_ROBIN = rob_err[800]
PM_ROBIN_ORDER = float(np.log2(rob_err[400] / rob_err[800]))
print(f"\npymrm film-BC solve vs shooting at the Bi_m = 3 threshold state: "
      + "  ".join(f"n_u={k}: {v:.2e}" for k, v in rob_err.items())
      + f"; order {PM_ROBIN_ORDER:.2f}")
print("The film condition enters ONLY through the bc dict {a:1, b:Bim, d:Bim} on the")
print("outward normal - the operators are untouched.")'''))

cells.append(md(r"""### 6. Axial dispersion: the criterion that holds

Eq. 30's boundary fixes $\mathrm{Pe}_L = 20\,n\ln(C_o/C_f)$, so the audit is a
one-parameter map over conversion. First order has the classical closed form;
the pymrm reactor provides the independent route and the only route for
$n = 2$ (eq. 32)."""))

cells.append(code(r'''XS = (0.2, 0.3935, 0.5, 0.632, 0.8, 0.9, 0.95, 0.99)
AX1 = {X: axial_dev(X) for X in XS}
print("eq. 30 boundary, first order (closed form): deviation in required length")
for X, d in AX1.items():
    print(f"  X = {X:5.3f}: Pe_L = {20*np.log(1/(1-X)):6.2f}   dev = {d*100:.3f} %")
AX1_MIN, AX1_MAX = min(AX1.values()), max(AX1.values())
print(f"\nThe claim HOLDS at every conversion: {AX1_MIN*100:.2f}-{AX1_MAX*100:.2f} %, "
      f"approaching 5 % from below")
print("as X -> 1. Mears' own criterion is the sharpest in the paper - the first-order")
print("perturbation term ln(Co/Cf)/Pe_L is exactly 5 % at the boundary and the dropped")
print("curvature is conservative.")

# --- the second reading of 'required length', computed beside the first -------
# The numbers above hold Pe_L at its boundary value and measure the ktau
# deficit - the perturbation-theory object, as defined in Parameters and
# assumptions. The fully self-consistent reading fixes k, v and D_a, so the
# length ratio l scales ktau AND Pe_L together; solve for l directly.
def axial_dev_covary(X):
    lnr = np.log(1.0 / (1.0 - X))
    return brentq(lambda l: frac_unconverted(l * lnr, l * 20.0 * lnr) - (1.0 - X),
                  1.0, 2.0, xtol=1e-13) - 1.0

AXC = {X: axial_dev_covary(X) for X in (0.2, 0.632, 0.99)}
print("\nthe same boundary under the co-varying-length reading (fixed k, v, D_a):")
for X, d in AXC.items():
    print(f"  X = {X:5.3f}: dev = {d*100:.3f} %   (fixed-Pe_L reading above: "
          f"{AX1[X]*100:.3f} %)")
print("Both readings hold the claim everywhere; the co-varying one is smaller")
print("because the longer bed also disperses relatively less (Pe_L grows with L).")

# --- the pymrm route, refined, outlet read at the boundary face ---------------
X_CHK = 0.632
ax_err = {}
for n_z in (200, 400, 800, 1600):
    d = axial_dev(X_CHK, n_z=n_z)
    ax_err[n_z] = abs(d - AX1[X_CHK])
PM_AX = ax_err[1600]
PM_AX_ORDER = float(np.log2(ax_err[800] / ax_err[1600]))
print(f"\npymrm reactor vs closed form at X = {X_CHK} (dev = {AX1[X_CHK]*100:.3f} %):")
print("  " + "  ".join(f"n_z={k}: {v:.2e}" for k, v in ax_err.items())
      + f"   order {PM_AX_ORDER:.2f}")

# --- second order: eq. 32's n-scaling, pymrm only -----------------------------
AX2 = {X: axial_dev(X, n_react=2, n_z=800) for X in (0.5, 0.9)}
for X, d in AX2.items():
    print(f"eq. 32 boundary, n = 2, X = {X}: dev = {d*100:.3f} %  (pymrm, n_z = 800)")
AX2_MAX = max(AX2.values())
print(f"The factor n in eq. 32 does its job: the claim holds for second order too "
      f"({100*AX2_MAX:.2f} % worst).")

# --- Petersen's asymptotic alpha < 1, for contrast ----------------------------
PET_ALPHA = {}
for X in (0.393, 0.632, 0.9):
    lnr = np.log(1 / (1 - X))
    Pe = brentq(lambda p: frac_unconverted(p, p) - (1 - X), 1e-3, 60, xtol=1e-13)
    PET_ALPHA[X] = Pe / lnr - 1
    print(f"eq. 28 boundary (alpha = 1), X = {X}: dev = {PET_ALPHA[X]*100:.1f} %")
PET_ALPHA_MAX = max(PET_ALPHA.values())
print(f"Petersen's asymptotic boundary admits {min(PET_ALPHA.values())*100:.0f}-"
      f"{PET_ALPHA_MAX*100:.0f} % - it 'indicates the transition")
print("region', exactly as the paper says, and exactly like his 13 in Section 4.")
print("\n[Trickle-flow remark, structural: eq. 30's minimum length scales as 1/Pe_a,")
print(" so Pe_a = 0.1 vs 2 costs a factor 20 in L/d_p - the paper's 'order-of-")
print(" magnitude greater in trickle-flow' is the criterion's own algebra.]")'''))

# ------------------------------------------------------------------- figures
cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
# left: the intraparticle family - exact 5 % thresholds vs printed bounds
ns = np.array([0.5, 1.0, 2.0, 3.0])
axes[0].loglog(ns, [ISO[n]["Phi5"] for n in ns], "o-", color="#1f6feb",
               label=r"exact $\Phi$ at 5 %")
axes[0].loglog(ns, 1 / ns, "s--", color="crimson", label=r"printed bound $1/n$")
axes[0].loglog(ns, 0.75 / ns, "^:", color="#8250df", label=r"unrounded $0.75/n$")
axes[0].set_xlabel(r"reaction order $n$"); axes[0].set_ylabel(r"$\Phi$")
axes[0].set_title("eq. 3: the rounding is the 5 % claim")
axes[0].legend(fontsize=8)
# right: eq. 9 across gamma*beta
gbs = np.array(sorted(MAP9))
axes[1].semilogy(gbs, [MAP9[g]["Phi5"] for g in gbs], "o-", color="#1f6feb",
                 label=r"exact $\Phi$ at 5 %")
bs = np.array([MAP9[g]["bound9"] for g in gbs])
axes[1].semilogy(gbs, bs, "s--", color="crimson", label=r"eq. 9 bound $1/|n-\gamma\beta|$")
axes[1].semilogy(gbs, 0.75 * bs, "d--", color="#d29922", lw=1.0, ms=4,
                 label=r"unrounded $0.75/|n-\gamma\beta|$")
axes[1].semilogy(gbs, [MAP9[g]["bound12"] for g in gbs], "^:", color="#3fb950",
                 label="eq. 12 (Weisz-Hicks)")
axes[1].axvspan(W_LO, W_HI, color="crimson", alpha=0.12)
axes[1].text(1.0, 2e-2, "window of\nno power", ha="center", fontsize=8, color="crimson")
axes[1].set_xlabel(r"$\gamma\beta$  ($n = 1$, $\gamma = 30$)"); axes[1].set_ylabel(r"$\Phi$")
axes[1].set_title("eq. 9: liberal by 0.75, void near compensation")
axes[1].legend(fontsize=8, loc="upper left")
fig.tight_layout(); plt.show()'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
# left: film closed forms
om = np.linspace(0, 0.6, 200)
for n, c in zip((0.5, 1.0, 2.0), ("#8250df", "#1f6feb", "#3fb950")):
    axes[0].plot(om, (1 - (1 - om / 3) ** n) * 100, color=c, label=rf"$n$ = {n}")
    axes[0].axvline(0.15 / n, color=c, ls=":", lw=0.8)
axes[0].axhline(5, color="0.4", lw=0.8, ls="--")
axes[0].set_xlabel(r"$\omega = \mathcal{R}r_p/(C_bk_c)$")
axes[0].set_ylabel("rate deviation, %")
axes[0].set_title(r"eq. 17: dev = $1-(1-\omega/3)^n$; exact at $n=1$")
axes[0].legend(fontsize=8)
# right: axial dev map
Xs = np.linspace(0.05, 0.995, 60)
axes[1].plot(Xs, [axial_dev(X) * 100 for X in Xs], color="#1f6feb", lw=1.6,
             label="first order (closed form)")
axes[1].plot(list(AX2), [d * 100 for d in AX2.values()], "s", color="#3fb950",
             label="second order (pymrm)")
axes[1].axhline(5, color="0.4", lw=0.8, ls="--")
axes[1].set_xlabel("conversion $X$")
axes[1].set_ylabel("deviation in required length, %")
axes[1].set_ylim(0, 6)
axes[1].set_title("eq. 30/32 at their own boundary: the claim that holds")
axes[1].legend(fontsize=8)
fig.tight_layout(); plt.show()'''))

# --------------------------------------------------- results: defect injection
cells.append(md(r"""### 7. Breaking it on purpose

Every reported metric gets a defect that must move it; a check that nothing can
move is decoration, and this repository has caught that defect on its own pages
often enough to make the table mandatory. The rows below inject the defects;
the page's final cell then maps **each agreement metric** to the row that moves
it. Structural items are labelled as such above where they occur (Section 1's
identities, the film closed-form collapse); they are transcription checks and
illustrations, not evidence, and they are not re-broken here.

Two rows exist because adversarial verification caught **live instances of
defect classes no earlier row could move**: the first staged version's fold
bound (5.7 %) and joint-pass worst (9.0 %) were both grid-limited extrema —
a supremum read off a scan grid that stops short of the band edge — and no
row varied a scan grid or the *selection* among multiple crossings of a
condition. Rows (x) and (xi) now inject exactly those defects against the
recast, converged estimators, and both move their metrics by far more than
any tolerance.

What no row in this table can detect, stated plainly: a *convention* shared by
both routes. If "deviation" meant something else to Mears than to this page, or
$r_p$ meant a diameter, every route here would move together and agree. Those
readings are guarded only by the native-resolution transcription (the
nomenclature page prints "$r_p$ = particle radius, cm") and by Section 1's
paper-internal identities — which is precisely why that section exists."""))

cells.append(code(r'''print(f"{'defect':<46} {'metric':<28} {'reference':>10} {'broken':>10}")

# (i) wrong geometry index: sphere solved as slab - break BOTH the 5 % state
#     and the flagship headline (dev admitted at eq. 1's bound)
s5_slab = cross(power(1.0), lambda p, e, P, o, ys: e - 0.95, lo=1e-3, nu=0)
print(f"{'nu = 0 (slab) in the pellet':<46} {'Phi at 5 %, n = 1':<28} "
      f"{PHI5_N1:10.4f} {s5_slab[3]:10.4f}")
s1_slab = cross(power(1.0), lambda p, e, P, o, ys: P - 1.0, lo=1e-3, nu=0)
DEV_AT_1_SLAB = 1 - s1_slab[2]
print(f"{'  same defect, headline metric':<46} {'eq. 1 dev at bound':<28} "
      f"{DEV_AT_1_N1:10.4f} {DEV_AT_1_SLAB:10.4f}")
print(f"{'':<46} (eta {s1_slab[2]:.4f} - reconciles with B1.4's published slab "
      f"0.6948)")

# (ii) rate order mis-set: n = 2 solved where n = 1 is claimed
print(f"{'source exponent 2 where 1 is claimed':<46} {'Phi at 5 %, n = 1':<28} "
      f"{PHI5_N1:10.4f} {ISO[2.0]['Phi5']:10.4f}")

# (iii) film area factor: a = 2/r_p (cylinder) in the mass balance
DEV_WRONG_A = 1 - (1 - 0.15 / 2) ** 1
print(f"{'film a = 2/r_p instead of 3/r_p':<46} {'eq. 17 dev at boundary, n=1':<28} "
      f"{0.05:10.4f} {DEV_WRONG_A:10.4f}")

# (iv) sign of the heat of reaction flipped in the film balance
print(f"{'dH sign flipped (endo treated as exo)':<46} {'eq. 14 dev at boundary':<28} "
      f"{FILM_HEAT[30.0][0]:10.4f} {-FILM_HEAT[30.0][1]:10.4f}")

# (v) sign of beta flipped in the Weisz-Hicks source
s5_flip = cross(wh(1.0, -0.5, GAMMA), lambda p, e, P, o, ys: abs(e - 1) - 0.05,
                lo=1e-4, hi=22.0)
print(f"{'beta sign flipped in the source (gb = 0.5)':<46} {'Phi at 5 %':<28} "
      f"{MAP9[0.5]['Phi5']:10.4f} {s5_flip[3]:10.4f}")

# (vi) Danckwerts inlet replaced by the naive Dirichlet inlet
lnr = np.log(1 / (1 - X_CHK)); Pe_ = 20 * lnr
f_bad = lambda k: axial_solve(Pe_, k, n_z=800, inlet="dirichlet")[0] - (1 - X_CHK)
kt_bad = brentq(f_bad, lnr, 2.5 * lnr, xtol=1e-12)
DEV_DIRICHLET = kt_bad / lnr - 1
print(f"{'Dirichlet inlet instead of Danckwerts':<46} {'axial dev at X = 0.632':<28} "
      f"{AX1[X_CHK]:10.4f} {DEV_DIRICHLET:10.4f}")

# (vii) outlet read at the last cell centre instead of the boundary face
c_face, c_all, _ = axial_solve(Pe_, lnr * (1 + AX1[X_CHK]), n_z=100)
c_last = float(c_all[-1])
OUTLET_GAP = abs(c_face - c_last) / c_face
print(f"{'outlet read at last cell centre (n_z = 100)':<46} {'outlet c/c*':<28} "
      f"{c_face:10.6f} {c_last:10.6f}")
print(f"{'':<46} {'(relative gap':<28} {OUTLET_GAP:10.1e} {'- the A3.7 trap)':>10}")

# (viii) robustness: reference integrator degraded 1e-10 -> 1e-5
s5_loose = cross(power(1.0), lambda p, e, P, o, ys: e - 0.95, lo=1e-3, rtol=1e-5)
TOL_SHIFT = abs(s5_loose[3] - PHI5_N1) / PHI5_N1
print(f"{'shooting rtol 1e-10 -> 1e-5':<46} {'Phi at 5 %, n = 1':<28} "
      f"{PHI5_N1:10.4f} {s5_loose[3]:10.4f}")

# (ix) a defect in ONE route only - the pymrm pellet solved as a slab while the
#      closed form stays spherical. The cross-route agreement class (spread,
#      pymrm-vs-shooting, orders) must explode, or it is decoration.
p_slab400 = Pellet(400, nu=0)
y_sl, rn_sl = p_slab400.solve(phiB, power(1.0))
assert rn_sl < RN_TOL
SPREAD_BROKEN = abs((1 - p_slab400.eta(y_sl, phiB)) - DEV_AT_1_CLOSED)
print(f"{'nu = 0 in the pymrm route ONLY':<46} {'three-route spread':<28} "
      f"{HEADLINE_SPREAD:10.1e} {SPREAD_BROKEN:10.4f}")

# (x) crossing-selection convention: report the kinetic-side (smallest-Phi)
#     band crossing where the page's convention is the deep-side one. Below
#     gb ~ 1.45 the band crossing is unique and NO defect of this class can
#     exist - the convention only becomes load-bearing where the band has
#     three crossings, which is why the earlier table missed it.
CR2 = CR[2.0][0]                       # sorted crossings at gb = 2 (Section 4)
print(f"{'kinetic-side crossing selected (gb = 2)':<46} {'Phi at 5 %, gb = 2':<28} "
      f"{CR2[-1]:10.4f} {CR2[0]:10.4f}")

# (xi) scan-grid density: the two grid-limited extrema adversarial
#      verification caught on the first staged version, re-injected verbatim
#      against the recast, converged estimators.
print(f"{'fold-band worst from 400-pt grid max':<46} {'eq9 fold worst certified':<28} "
      f"{WORST_DEV_CERT_EQ9:10.4f} {GM9_400:10.4f}")
worst_coarse = 0.0                     # the first version's estimator, verbatim
for Bim_c in (0.3, 1.0, 3.0, 10.0, 30.0, 100.0):
    for t in np.logspace(-4, 1.25, 110):
        o = shoot(power(1.0), t, bim=Bim_c)
        if o is None:
            continue
        _, eta_b, Phi_b, om, ys = o
        if om < 0.15 and Phi_b < 1 / (1 + 0.33 * om) and Phi_b / ys < 1.0:
            worst_coarse = max(worst_coarse, 1 - eta_b)
JOINT_COARSE = worst_coarse
print(f"{'joint pass from the coarse Bi_m x t grid':<46} {'joint-pass worst dev':<28} "
      f"{JOINT_WORST:10.4f} {JOINT_COARSE:10.4f}")

print(f"\nEvery injected defect moves its metric except the deliberate robustness row")
print(f"(viii), which moves it by {TOL_SHIFT:.1e} - the thresholds are not artefacts "
      f"of solver tolerance.")
print("Row (vii) is the handoff's Neumann-outflow lesson measured in place: the")
print("zero-gradient outlet extrapolates to the face, and the last cell centre")
print(f"disagrees at first order ({OUTLET_GAP:.1e} at n_z = 100). It is small HERE")
print("because a Danckwerts outlet is flat - but 'small here' is a measured fact,")
print("not an assumption, and reading the face costs nothing.")'''))

# ------------------------------------------------------------------ validation
cells.append(md(r"""## Validation

**Tier 6. Nothing on this page is compared with a measurement** — the source
contains none. The evidence is ranked below, strongest first; every pymrm solve
behind a reported number asserts its Newton residual against `RN_TOL`.

1. **Three routes to the headline number** — the deviation admitted by
   eq. 1 — closed form, pymrm finite volume, and the shooting reference, which
   share respectively nothing and only the algebraic rate law.
2. **pymrm against the shooting reference at every audited threshold state**,
   Dirichlet and film boundary conditions, isothermal and not, at measured
   order ≈ 2; the axial reactor against the classical closed form at its
   measured order.
3. **The paper's own printed numbers reproduced from other printed numbers**
   (Section 1): 0.22 = $1/\sqrt{20}$, "four and a half" = $\sqrt{20}$,
   Carberry's 0.1 = 2 × Mears' 0.15, eq. 15's 10 = 2 × 0.75/0.15, and the
   measured perturbation slope $n/15$ that makes 0.75 the constant the paper
   says it rounded. None is used as an input to anything.
4. **The paper's prose statements measured**: Petersen's 13 admits ~38 %
   (a transition-region value, as the paper's own axial discussion explains);
   the "80 %" compensation statement holds conservatively (0.45); eq. 12 is
   safe exothermic and fails endothermic, now with numbers; zero order is
   exactly unaffected at the film and exactly bounded by 6 in the pellet.
5. **The defect-injection table**, Section 7 — including the deliberate
   null row (solver tolerance), the stated blind spot (shared conventions),
   and two rows injected after adversarial verification caught live
   grid-limited extrema on the first staged version (crossing-selection
   convention; scan-grid density — the fold bound and the joint-pass worst
   were both biased low by exactly those defects, and their cells were recast
   to root-found, grid-independent estimators). The final cell maps **every**
   reported agreement metric to the break row that moves it, and labels the
   structural ones.

The structural items — Section 1's algebraic identities and the film
closed-form collapse — are labelled where they occur and are counted as
transcription checks, not as evidence for any audit number. Three metrics sit
below `check_agreement`'s `ABS_FLOOR` of 1e-12 (the two dimensional film
cross-checks and the zero-order profile error): they are **outside the CI
regression suite** — protected only by this page's own asserts — and the final
cell names them."""))

cells.append(code(r'''RUNTIME_S = time.time() - T_START
print("=== the audit, in numbers ===")
print(f"eq. 1/3 (intraparticle): dev at printed bound {DEV_AT_1_N1*100:.2f} % (n=1; "
      f"{min(v['dev_bound'] for v in ISO.values())*100:.2f}-"
      f"{max(v['dev_bound'] for v in ISO.values())*100:.2f} % over n),")
print(f"  at unrounded 0.75/n {min(v['dev_075'] for v in ISO.values())*100:.2f}-"
      f"{max(v['dev_075'] for v in ISO.values())*100:.2f} %; Weisz 0.6/0.3 admit "
      f"{ISO[1.0]['dev_weisz']*100:.2f}/{ISO[2.0]['dev_weisz']*100:.2f} %; zero order "
      f"exact at 6")
print(f"eq. 2 extensions: n = -1 admits {(NEG_ETA_AT_1-1)*100:.1f} % at its bound; "
      f"inhibited LH passes eta = {LH_ETA_AT_BOUND:.3f}")
print(f"eq. 17 (interphase mass): EXACT 5.000 % at n = 1; "
      f"{FILM_MASS[0.5]*100:.2f} % (n=0.5) to {FILM_MASS[3.0]*100:.2f} % (n=3)")
print(f"eq. 14 (interphase heat): {HEAT_EXO_MIN*100:.2f}-{HEAT_EXO_MAX*100:.2f} % "
      f"exothermic (always > 5), {FILM_HEAT[30.][1]*100:.2f} % endothermic")
print(f"eq. 16 (Carberry): {CARBERRY_DEV_SURF*100:.1f} % (surface-referenced eta, "
      f"Mears' definition) to {CARBERRY_DEV*100:.0f} % (bulk)")
print(f"eq. 9: dev at bound {min(endo_devs)*100:.1f}-{max(endo_devs)*100:.1f} % on "
      f"the endothermic side; even the unrounded 0.75")
print(f"  form over-certifies for gb in ({W_LO:.2f}, {W_HI:.2f}) - the window where "
      f"no constant fix works -")
print(f"  and re-admits a mild {min(OUT2.values())*100:.2f}-"
      f"{max(OUT2.values())*100:.2f} % past gb = {SECOND_ONSET:.2f} (three-crossing "
      f"regime, to gb = 3)")
print(f"  fold case (gb = 24): eq. 1's band certifies {WORST_DEV_CERT_EQ1:.0f}; "
      f"eq. 9's band edge admits {WORST_DEV_CERT_EQ9*100:.1f} % (root-found)")
print(f"  Petersen's 13 admits {PET13_DEV*100:.0f} % (gamma = 30); compensation-case "
      f"5 % point Phi = {PHI5_COMP:.2f}, centre y = {YC_COMP:.2f}")
print(f"eq. 12: conservative exothermic; admits {DEV12_ENDO*100:.0f} % at gb = -1 "
      f"(the endothermic failure, quantified)")
print(f"eq. 6 (Anderson): worst {AND6_MAX*100:.2f} % at its boundary over the sweep")
print(f"eq. 21 (isothermal part): ratio to exact {RATIO21_HI:.2f} at Bi_m = 100 -> "
      f"{RATIO21_LO:.3f} at Bi_m = 0.3 (structural failure);")
print(f"  passing eqs. 1 + 17 + 21 at once still admits {JOINT_WORST*100:.1f} % "
      f"(interior optimum, Bi_m = {JOINT_BIM:.1f})")
print(f"eq. 30/32 (axial): {AX1_MIN*100:.2f}-{AX1_MAX*100:.2f} % (n=1), "
      f"{AX2_MAX*100:.2f} % worst (n=2) - the claim that HOLDS;")
print(f"  co-varying-length reading {AXC[0.2]*100:.2f}-{AXC[0.99]*100:.2f} % - both "
      f"readings hold it")
print(f"eq. 28 (Petersen axial): {PET_ALPHA_MAX*100:.0f} % worst at alpha = 1")
print(f"\nsolver cross-checks: headline three-route spread {HEADLINE_SPREAD:.1e}; "
      f"pymrm vs shooting {PM_ISO:.1e}")
print(f"  (iso, order {PM_ISO_ORDER:.2f}), {PM_WH:.1e} (nonisothermal, order "
      f"{PM_WH_ORDER:.2f}), {PM_ROBIN:.1e} (film bc, order {PM_ROBIN_ORDER:.2f});")
print(f"  axial vs closed form {PM_AX:.1e} at n_z = 1600 (order {PM_AX_ORDER:.2f}); "
      f"dimensional film checks {FILM_MASS_XCHK:.1e} / {FILM_HEAT_XCHK:.1e}")
print(f"\nruntime so far: {RUNTIME_S:.0f} s")

METRICS = {
    # what each criterion admits at its own printed boundary (the audit itself)
    "eq1_dev_at_bound_n1": float(DEV_AT_1_N1),
    "eq3_dev_at_bound_n05": float(ISO[0.5]["dev_bound"]),
    "eq3_dev_at_bound_n2": float(ISO[2.0]["dev_bound"]),
    "eq3_dev_at_bound_n3": float(ISO[3.0]["dev_bound"]),
    "eq3_dev_at_unrounded_075_n1": float(DEV_AT_075_N1),
    "weisz_06_dev_n1": float(ISO[1.0]["dev_weisz"]),
    "weisz_03_dev_n2": float(ISO[2.0]["dev_weisz"]),
    "neg_order_eta_at_bound": float(NEG_ETA_AT_1),
    "hudgins_LH_eta_at_bound": float(LH_ETA_AT_BOUND),
    "LH_eta_at_Phi1_vs_B14": float(LH_ETA_AT_1),
    "eq17_dev_at_bound_n1": float(0.05),
    "eq17_dev_at_bound_n05": float(FILM_MASS[0.5]),
    "eq17_dev_at_bound_n2": float(FILM_MASS[2.0]),
    "eq14_dev_exo_gamma30": float(FILM_HEAT[30.0][0]),
    "eq14_dev_endo_gamma30": float(FILM_HEAT[30.0][1]),
    "eq16_carberry_dev": float(CARBERRY_DEV),
    "eq16_carberry_dev_surface": float(CARBERRY_DEV_SURF),
    "eq9_dev_at_bound_gb05": float(DEV9_GB05),
    "eq9_dev_at_bound_gb08": float(DEV9_GB08),
    "eq9_nopower_window_lo": float(W_LO),
    "eq9_nopower_window_hi": float(W_HI),
    "eq9_Phi5_gb2_kinetic": float(CR[2.0][0][0]),
    "eq9_second_regime_onset": float(SECOND_ONSET),
    "eq9_outside_window_dev_gb2": float(OUT2[2.0]),
    "petersen13_dev_gamma30": float(PET13_DEV),
    "compensation_Phi5_gamma30": float(PHI5_COMP),
    "compensation_centre_y": float(YC_COMP),
    "eq12_dev_endo_gb_m1": float(DEV12_ENDO),
    "eq6_worst_devT": float(AND6_MAX),
    "eq21_ratio_Bim100": float(RATIO21_HI),
    "eq21_ratio_Bim03": float(RATIO21_LO),
    "joint_pass_worst_dev": float(JOINT_WORST),
    "joint_pass_worst_bim": float(JOINT_BIM),
    "eq30_dev_X0632": float(AX1[0.632]),
    "eq30_dev_X099": float(AX1[0.99]),
    "axial_dev_covarying_X0632": float(AXC[0.632]),
    "eq32_dev_n2_worst": float(AX2_MAX),
    "eq28_petersen_dev_X09": float(PET_ALPHA[0.9]),
    # the fold comparison (Section 4)
    "eq1_fold_worst_certified_dev": float(WORST_DEV_CERT_EQ1),
    "eq9_fold_worst_certified_dev": float(WORST_DEV_CERT_EQ9),
    "eq9_fold_ignited_margin": float(MIN_PHI_IGNITED / b9_f),
    # solver cross-checks (the checks that can fail)
    "headline_three_route_spread": float(HEADLINE_SPREAD),
    "pymrm_vs_shooting_iso": float(PM_ISO),
    "pymrm_vs_shooting_iso_order": float(PM_ISO_ORDER),
    "pymrm_vs_shooting_wh": float(PM_WH),
    "pymrm_vs_shooting_robin": float(PM_ROBIN),
    "pymrm_axial_vs_closed_form": float(PM_AX),
    "film_mass_dimensional_xcheck": float(FILM_MASS_XCHK),
    "film_heat_dimensional_xcheck": float(FILM_HEAT_XCHK),
    "zero_order_profile_err": float(ZERO_PROFILE_ERR),
    # structural / transcription identities, labelled as such
    "alpha_star_identity": float(ALPHA_STAR),
    "bi_star_identity": float(BI_STAR_UNROUNDED),
    "film_closed_form_collapse_rel": float(COLLAPSE),
}
report_agreement("B1.7", METRICS)

# ------- break-row coverage: every metric names the row that moves it ---------
# (i) slab nu, (ii) wrong exponent, (iii) film area factor, (iv) dH sign,
# (v) beta sign, (vi) Dirichlet inlet, (vii) last-cell outlet, (viii) null
# robustness row, (ix) single-route slab, (x) crossing-selection convention,
# (xi) scan-grid density. STRUCT = labelled structural where it occurs: an
# identity or transcription check, kept as such, not evidence.
PELLET, FILM_M, FILM_H, WH_ROWS, AXIAL, XROUTE = \
    "(i)+(ii)", "(iii)", "(iv)", "(v)", "(vi)+(vii)", "(ix)"
COVERAGE = {
    "eq1_dev_at_bound_n1": f"(i) directly ({DEV_AT_1_N1:.4f} -> {DEV_AT_1_SLAB:.4f})",
    "eq3_dev_at_bound_n05": PELLET, "eq3_dev_at_bound_n2": PELLET,
    "eq3_dev_at_bound_n3": PELLET, "eq3_dev_at_unrounded_075_n1": PELLET,
    "weisz_06_dev_n1": PELLET, "weisz_03_dev_n2": PELLET,
    "neg_order_eta_at_bound": PELLET,
    "hudgins_LH_eta_at_bound": PELLET,
    "LH_eta_at_Phi1_vs_B14": PELLET + " + reconciled with B1.4's printed 1.0483",
    "eq17_dev_at_bound_n1": FILM_M, "eq17_dev_at_bound_n05": FILM_M,
    "eq17_dev_at_bound_n2": FILM_M,
    "eq14_dev_exo_gamma30": FILM_H, "eq14_dev_endo_gamma30": FILM_H,
    "eq16_carberry_dev": FILM_M,
    "eq16_carberry_dev_surface": FILM_M,
    "eq9_dev_at_bound_gb05": WH_ROWS, "eq9_dev_at_bound_gb08": WH_ROWS,
    "eq9_nopower_window_lo": WH_ROWS, "eq9_nopower_window_hi": WH_ROWS,
    "eq9_Phi5_gb2_kinetic":
        f"(x) directly (the row swaps it with the deep-side {CR2[-1]:.4f})",
    "eq9_second_regime_onset": WH_ROWS + " + (x) class",
    "eq9_outside_window_dev_gb2": WH_ROWS + " + (x) class",
    "petersen13_dev_gamma30": WH_ROWS,
    "compensation_Phi5_gamma30": WH_ROWS, "compensation_centre_y": WH_ROWS,
    "eq12_dev_endo_gb_m1": WH_ROWS, "eq6_worst_devT": WH_ROWS,
    "eq21_ratio_Bim100": PELLET, "eq21_ratio_Bim03": PELLET,
    "joint_pass_worst_dev":
        f"(xi) directly ({JOINT_WORST:.4f} -> {JOINT_COARSE:.4f}); "
        + PELLET + "+" + FILM_M + " class too",
    "joint_pass_worst_bim": "(xi) (the coarse sweep has no interior optimum)",
    "eq30_dev_X0632": AXIAL, "eq30_dev_X099": AXIAL,
    "axial_dev_covarying_X0632": AXIAL,
    "eq32_dev_n2_worst": AXIAL, "eq28_petersen_dev_X09": AXIAL,
    "eq1_fold_worst_certified_dev": WH_ROWS + " + (xi) class",
    "eq9_fold_worst_certified_dev":
        f"(xi) directly ({WORST_DEV_CERT_EQ9:.4f} -> {GM9_400:.4f})",
    "eq9_fold_ignited_margin": WH_ROWS,
    "headline_three_route_spread":
        f"(ix) directly ({HEADLINE_SPREAD:.1e} -> {SPREAD_BROKEN:.2f})",
    "pymrm_vs_shooting_iso": XROUTE, "pymrm_vs_shooting_iso_order": XROUTE,
    "pymrm_vs_shooting_wh": XROUTE, "pymrm_vs_shooting_robin": XROUTE,
    "pymrm_axial_vs_closed_form": XROUTE + " class; (vi) moves it too",
    "film_mass_dimensional_xcheck": FILM_M,
    "film_heat_dimensional_xcheck": FILM_H,
    "zero_order_profile_err": PELLET,
    "alpha_star_identity": "STRUCT: transcription identity (0.22 vs 20)",
    "bi_star_identity": "STRUCT: transcription identity (10 vs 0.75/0.15)",
    "film_closed_form_collapse_rel": "STRUCT: identity of the film model",
}
missing, extra = set(METRICS) - set(COVERAGE), set(COVERAGE) - set(METRICS)
assert not missing and not extra, (missing, extra)
print("\nbreak-row coverage (which injected defect moves each reported metric):")
for k in METRICS:
    print(f"  {k:<30} {COVERAGE[k]}")
ABS_FLOOR = 1e-12
below = sorted(k for k, v in METRICS.items() if abs(v) < ABS_FLOOR)
print(f"\nmetrics below check_agreement's ABS_FLOOR = 1e-12, therefore OUTSIDE the")
print(f"CI regression suite and protected only by this page's own asserts: {below}")'''))

# --------------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Nothing to the criteria themselves** — each is one line of algebra,
reproduced exactly as printed. What is added is the part 1971 could not afford:
the exact solves behind every stated tolerance, which turn a family of rules of
thumb into a table of measured guarantees. Per criterion:

- **Eq. 17 is exact at first order** — the film balance makes the deviation
  equal to $\omega/3$ identically. The most-quoted "Mears criterion" is also
  the best-behaved one in the paper. Zero order is exactly unaffected, as the
  paper says.
- **The intraparticle family (eqs. 1-3, 9) is uniformly liberal by its own
  rounding.** The measured perturbation slope makes 0.75 the honest constant;
  at the printed bounds the exact deviation is 6.3-6.8 %. An experimenter who
  wants the stated 5 % should use $0.75/n$ — the number the paper had before
  rounding it away.
- **Eq. 9's exclusion zone has a width**: even with the unrounded 0.75
  restored, the bound certifies deviations beyond 5 % in a measured window of
  $\gamma\beta$ around $n$, and Petersen's 13 there is a transition-region
  marker admitting ~38 %, not a 5 % test. Past $\gamma\beta \approx 1.95$ —
  where the $\pm5\,\%$ band has three crossings and the bound overtakes the
  kinetic-side one — the unrounded bound re-admits a mild 5.0-5.2 %, the same
  curvature class as eq. 14. And eq. 9 is *immune* to the observable fold
  that `B1.4` showed defeats eq. 1: where $\Phi$ folds, $|n-\gamma\beta|$ is
  large and eq. 9's bound has already shrunk below the ignited branch's reach
  — its band edge still admits the rounding-class ≈7 % (root-found; a grid
  max under-reads it, which is row (xi) of the break table), nothing ignited.
- **Eq. 14 is liberal by construction but only just** (5.10-5.12 % for $\gamma$ = 10-40), and film-lumped observables *determine* the deviation — the
  multiplicity blind spot the pellet observable has does not exist at the film.
- **Eq. 21 fails structurally at low $\mathrm{Bi}_m$**: its left side vanishes
  with the pellet gradient while the deviation moves into the film. It must be
  applied together with eq. 17, which the paper nowhere states. And passing
  eqs. 1, 17 and 21 simultaneously still admits ~11 %, at an *interior*
  optimum $\mathrm{Bi}_m \approx 6$ where two criteria bind at once — 5 %
  tolerances stack, and no criterion in the paper says so.
- **Eq. 30/32 — Mears' axial criterion — holds its claim everywhere tested**,
  3.9-4.95 % across conversions and orders, asymptotically sharp; Petersen's
  $\alpha < 1$ admits 23-51 % at the same states.

The pymrm content: one `Pellet` class whose surface physics — Dirichlet or
film — is *entirely* a boundary-condition dict on the outward normal, the
operators untouched; and the `A2.1` reactor assembly reused with the outlet
read through `compute_boundary_values`. The deterministic shooting reference
and the closed forms give the finite-volume solver something to be checked
against at measured order, which is what makes the audit's numbers evidence
rather than output."""))

# --------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

| Change | Where |
|---|---|
| Audit your own criterion | write it as `val(phi, eta, Phi, omega, y_s)` and hand it to `cross` — the branch machinery finds the exact state where it binds; note the default returns the *deep-side* (largest-Φ) crossing, so use `which="all"` wherever the condition can cross more than once |
| Different kinetics | any `R(y)` into `shoot`/`Pellet.solve`; keep $\mathcal{F}(1)=1$ so $\Phi$ stays the printed observable |
| Film instead of Dirichlet surface | `Pellet(bim=...)` — one bc dict, `{a:1, b:Bim, d:Bim}` on the outward normal |
| Slab or cylinder | `nu=0/1` — but the printed thresholds are for the sphere, and the break table shows how much that matters |
| Different tolerance than 5 % | change the `0.05` in the `cross` calls; the bounds scale linearly only at first order, which is the whole point of the page |
| Axial criterion at your conversion and order | `axial_dev(X, n_react=...)` |

**How to actually screen a measurement, given what is above.** Use eq. 17 as
printed — it is exact at first order and near-exact otherwise. For the pellet,
use $0.75/n$, not $1/n$, if the 5 % is meant literally, and state $r_p$ is a
radius. If the reaction is exothermic, compute $\gamma\beta$ first: inside the
measured window around $n$ (Section 4) not even the unrounded constant makes
eq. 9 a 5 % test — measure at two particle sizes instead (the paper's own
empirical fallback). Apply eq. 17
*alongside* eq. 21, never eq. 21 alone, and remember that passing every test
still admits ~11 % when two resistances sit at their boundaries at once. The
axial criterion can be trusted as printed.

**Reading the paper.** The scan's text layer destroys the criteria (eq. 25's
5.3 comes back as `°-`; $\eta$ extracts as a comma); read every threshold off a
300 dpi render — the native resolution, so higher renders only interpolate.
The PDF opens with the tail of the preceding article; the paper starts partway
down page 1 of the file.

**Related pages.** [`B1.4`](../B1.4-weisz-prater-criterion/) — the same audit
for eq. 1 alone, including the fold analysis this page leans on;
[`B1.1`](../B1.1-thiele-weisz-hicks/) — the forward problem behind Section 4;
[`B1.6`](../B1.6-prater-relation/) — the reduction used in the pellet model;
[`A2.1`](../A2.1-danckwerts-boundary-conditions/) — the axial vessel and its
boundary conditions; `B1.2`/`B1.3` — the generalised moduli, i.e. the
length-scale question underneath every threshold here.

## References

Mears, D. E. (1971). *Tests for transport limitations in experimental
catalytic reactors.* Ind. Eng. Chem. Process Des. Develop. **10**(4), 541-547.
[doi:10.1021/i260040a020](https://doi.org/10.1021/i260040a020) — **the source
read**; every criterion and threshold transcribed from 300 dpi native renders
of journal pages 541-547.

Mears, D. E. (1971a). Chem. Eng. Sci. — cited by the source as "in press"; the
**origin of the axial criterion** (eqs. 29-32). **Not consulted.**

Mears, D. E. (1971b). J. Catal. **20**, 127 — the **origin of the interphase
and interparticle heat criteria** (eqs. 14, 23). **Not consulted.**

Cited by the source for the criteria it reviews, none consulted: Weisz &
Prater (1954); Weisz (1957); Weisz & Hicks (1962) — on disk as the source of
`B1.1`/`B1.4`/`B1.6`, not needed here; Hudgins (1968); Anderson (1963); Kubota
& Yamanaka (1969); Carberry (1961); Satterfield, Pelossof & Sherwood (1969);
Petersen (1965b, 1968); Stewart & Villadsen (1969); van den Bleek, van der
Wiele & van den Berg (1969)."""))

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
print(f"wrote {out}")
