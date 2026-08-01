#!/usr/bin/env python3
"""Generate index.ipynb for page B1.6 (the Prater relation). Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "The Prater relation: a pellet's temperature is slaved to its concentration"
description: "One linear identity fixes the temperature everywhere inside a catalyst particle from the concentration alone - whatever the kinetics, whatever the shape. Proved symbolically, confirmed numerically, measured against deliberately injected defects to establish what that confirmation is worth, and then broken on purpose."
categories: [sec:B, struct:S3, tier:T0, data:tier6, phase:gas-solid]
date: 2026-07-31
---

# The Prater relation: a pellet's temperature is slaved to its concentration

**Catalog ID:** `B1.6` · **Structures:** `S3` (1D steady BVP) · **Tier:** T0

Inside a working catalyst particle, reactant is consumed and heat is released at
the same points, in fixed proportion. Prater showed in 1958 that this makes the
temperature field a **linear function of the concentration field** — not
approximately, and not only for a sphere:

$$\lambda\,(T - T_s) \;=\; (-\Delta H)\,D_e\,(c_s - c).$$

No rate constant appears. No activation energy. No particle shape. Solve the
diffusion problem and the temperature comes free.

This page proves that, confirms it numerically in pymrm across three geometries
and three rate laws, establishes **what that numerical confirmation is and is
not sensitive to** by injecting defects on purpose, and then does the part the
1958 result could not: **finds where it stops being true, and by how much.**"""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

A porous catalyst pellet running an exothermic reaction is governed by two
balances that share one source term. Species is consumed at the local rate
$r$; heat is produced at $(-\Delta H)\,r$ at exactly the same place. The two
fields are therefore not independent, and Prater's observation is that they are
not even *approximately* independent — one determines the other exactly.

The practical consequence is the reason every non-isothermal pellet analysis
since is a **one-variable** problem. Weisz and Hicks (1962) say so in as many
words: the relation "enables us to write the reaction rate (6) in terms of one
variable instead of both, $c$ and $\Delta T$". Their whole paper, the
[`B1.1`/`B1.5`](../B1.1-thiele-weisz-hicks/) page's $\eta > 100$ and its three
coexisting steady states, all rest on it.

It also gives the single most useful number in non-isothermal catalysis for
free. Setting $c \to 0$ — a fully diffusion-starved centre — bounds the internal
temperature rise:

$$(\Delta T)_{\max} = \frac{(-\Delta H)\,D_e\,c_s}{\lambda},$$

which needs no kinetics at all, only transport properties and the surface
concentration. Divided by $T_s$ it is the **Prater number** $\beta$, the
parameter labelling the curves on the Weisz–Hicks diagram.

Why does an identity deserve a gallery page? Because it is an unusually clean
example of the second-highest kind of validation available to us — an
**internal identity the model must satisfy** — and because its assumptions are
usually left unstated. This page states them, and then breaks each one on
purpose."""))

# ------------------------------------------------------------ the published model
cells.append(md(r"""## The published model

### Provenance: what was read, and what was not

**Prater's own 1958 paper is not on disk and has no open-access route.** It is
cited here as the origin of the result, from Weisz and Hicks' own reference
list, which reads:

> [11] PRATER C. D., *Chem. Engng. Sci.* 1958 **8** 284.

Everything transcribed below comes instead from **Weisz and Hicks (1962)**,
which is on disk and is already the source for the published
[`B1.1`/`B1.5`](../B1.1-thiele-weisz-hicks/) page. They print the relation as
their **eq. (7)** and attribute it explicitly: Damköhler had it for the sphere,
and "this relationship was proven by PRATER [11] to apply to *any* particle
geometry."

So: **the equations on this page are Weisz and Hicks as printed; the result is
Prater's.** Nothing here is restated from memory or from a textbook. The scan's
text layer is badly mangled — it renders eq. (7) as `AT = T - To = - F` /
`[grad ~1,~~` — so every symbol below was read off a 600 dpi render of journal
page 266 (`pdftoppm -r 600 -f 2 -l 2 -png`).

### The equations, as printed

Weisz and Hicks define $K$ = thermal conductivity, $H$ = molar heat of reaction,
and $\mathrm{d}n/\mathrm{d}t$ = "the actual per unit volume reaction rate in a
volume element". Their transport equations are

$$D\nabla^2 c - \frac{\mathrm{d}n}{\mathrm{d}t} = 0 \tag{4}$$

$$K\nabla^2 T - H\frac{\mathrm{d}n}{\mathrm{d}t} = 0 \tag{5}$$

and the relation itself, printed exactly as

$$\Delta T = T - T_0 = -\frac{HD}{K}\,(c_0 - c), \tag{7}$$

"with $T_0$ and $c_0$ being the boundary values". The dimensionless groups
follow as

$$\gamma = \frac{Q}{RT_0}, \qquad
\beta = \frac{c_0 H D}{K T_0} = \left(\frac{\Delta T}{T_0}\right)_{\!\max},
\qquad y = \frac{c}{c_0}. \tag{8}$$

### One sign has to be resolved before anything is built

Eq. (7) is *exactly* what eqs. (4) and (5) give — the check below confirms it
symbolically. But eq. (7) at $c \to 0$ gives
$(\Delta T/T_0)_{\max} = -c_0HD/(KT_0)$, which is the **negative** of the middle
member of eq. (8). The two printed expressions for $\beta$ therefore disagree in
sign, and the disagreement is not cosmetic: it decides whether an exothermic
pellet is hotter or colder inside.

Three pieces of the paper's own content settle it, and all three agree:

1. Eq. (8)'s rate, $\mathrm{d}n/\mathrm{d}t = k_0c_0\,y\exp\{\gamma\beta(1-y)/[1+\beta(1-y)]\}$,
   is eq. (6) with $\Delta T/T_0$ replaced by $\beta(1-y)$. That substitution
   requires $\beta = +(\Delta T/T_0)_{\max}$, the right-hand member.
2. The paper sweeps "$\beta$ in the range from 0 to $+0.8$ (exothermic
   reaction) and 0 to $-0.8$ (endothermic)".
3. Its central result is $\eta \gg 1$, which needs the exponent to be
   *positive* where the pellet is depleted — i.e. hotter inside.

So $H$ in eqs. (5) and (7) is the molar **enthalpy** of reaction, negative when
heat is released, and $\beta$ is defined positive for an exothermic reaction. In
modern notation, with $(-\Delta H) > 0$ exothermic, $\lambda$ the effective
conductivity and $D_e$ the effective diffusivity, the two printed statements
become unambiguous and are what this page uses throughout:

$$\boxed{\;\lambda\,(T - T_s) = (-\Delta H)\,D_e\,(c_s - c)\;}
\qquad
\beta = \frac{(-\Delta H)\,D_e\,c_s}{\lambda\,T_s}
= \left(\frac{\Delta T}{T_s}\right)_{\!\max}.$$

This is the same $\beta$ that labels the curves on
[`B1.1`/`B1.5`](../B1.1-thiele-weisz-hicks/); the cross-check below confirms the
two pages use it identically."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

The relation holds under five assumptions, and this page tests each of them:

| # | assumption | what happens without it |
|---|---|---|
| 1 | steady state | fails; quantified below, exactly recovered when $Le = 1$ |
| 2 | one reaction, so both balances share one $r$ | fails (not tested here) |
| 3 | constant $D_e$ and $\lambda$ | the linear relation between $T$ and $c$ fails. If $D_e = D_e(c)$ and $\lambda = \lambda(T)$, the Kirchhoff potentials $\Phi_c = \int\! D_e\,\mathrm{d}c$ and $\Phi_T = \int\! \lambda\,\mathrm{d}T$ satisfy $\nabla^2[\Phi_T + (-\Delta H)\Phi_c] = 0$, so what survives is a linear relation between the *potentials*; between $T$ and $c$ it becomes nonlinear (not tested here) |
| 4 | the **same** boundary condition type for both fields | fails; quantified below |
| 5 | — nothing about the kinetics, and nothing about the geometry — | *these are the point: they are not assumptions* |

Everything on this page is dimensionless. With $u$ the position scaled on the
characteristic length $L$, $y = c/c_s$ and $\theta = T/T_s$, the steady pellet is

$$\nabla^2 y = \phi^2\,\mathcal{R}(y,\theta), \qquad
\nabla^2 \theta = -\beta\,\phi^2\,\mathcal{R}(y,\theta),$$

with $\phi = L\sqrt{k_s/D_e}$ the Thiele modulus, $\gamma = E/(R_gT_s)$ the
Arrhenius number, $\beta$ the Prater number, and $\nabla^2$ carrying the
geometry index $\nu = 0, 1, 2$ for slab, cylinder and sphere. **The two
equations are solved as two coupled fields and the Prater relation is never
substituted anywhere in the solver** — that is the whole design of the page.

Three rate laws are used, all pointwise in $(y,\theta)$:

$$\mathcal{R}_{1} = y\,e^{\gamma(1-1/\theta)}, \qquad
\mathcal{R}_{\mathrm{LH}} = \frac{y\,e^{\gamma(1-1/\theta)}}{\left[1 + K(\theta)\,y\right]^2},
\quad K(\theta) = K_s e^{\gamma_a(1/\theta-1)}, \qquad
\mathcal{R}_{0} = e^{\gamma(1-1/\theta)} .$$

$\mathcal{R}_1$ is Weisz and Hicks' own kinetics written in $\theta$; the other
two are chosen because the identity claims they cannot matter. The
Langmuir–Hinshelwood form carries a *second* temperature dependence, in the
adsorption constant, which is the sharpest available test of that claim."""))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

**There is none, and there does not need to be — provenance tier 6.**

This page compares a model with an identity, not with a measurement. Weisz and
Hicks (1962) report no experimental data at all (the same finding recorded for
`B1.1`), and Prater (1958) is a theoretical note. The validation route used here
is the second-ranked one in the gallery's brief — *an internal identity the
model must satisfy* — which needs no dataset and no digitised figure. Nothing on
this page is digitised and no CSV is shipped.

**This is a proof plus numerical confirmation. It is not experimental
validation, and the two must not be blurred.** A pellet whose measured centre
temperature matched $\beta T_s(1-y_c)$ would be validation; nothing below is.

The only inputs are the printed equations transcribed above, and the checks in
the next section are what stands in for a data comparison."""))

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
from scipy.integrate import solve_ivp, solve_bvp
from scipy.sparse import diags_array
from pymrm import (construct_grad, construct_div, construct_boundary_value_matrices,
                   NumJac, newton, clip_approach)
from gallery_utils import report_agreement

PAGE = "B1.6-prater-relation"
np.seterr(all="ignore")
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

One class, and the only thing that matters about it is what it does **not**
contain: nowhere does `TwoFieldPellet` use $\theta = 1 + \beta(1-y)$. It
discretises the two balances as an honest coupled system and lets the solver
find both fields.

Three conventions, all of which this page depends on:

- **Boundary conditions use the OUTWARD normal**, and are written per field. At
  the centre the outward normal points inward, so symmetry is
  $\partial y/\partial n = 0$; at the surface it points outward, so a film
  condition reads $\partial y/\partial n + Bi\,y = Bi$ — *the same sign at both
  ends of the film pair*, which is the trap. Each `bc` below carries its
  physical equation in a comment.
- **Spatial axis first, fields last**: the layout is `(n_u, 2)` with column 0
  the concentration and column 1 the temperature, everywhere.
- **`NumJac(shape)` with the default stencil**, which couples only the last
  axis. That is exactly right here: the rate is pointwise in $u$ and couples
  $y$ and $\theta$ *at the same point*, which is the last-axis block. No
  `axes_diagonals` is needed, because nothing in the source term reads a
  neighbouring cell.

Constant operators are assembled once in `__init__`; the only thing rebuilt per
Newton iteration is the pointwise source block."""))

cells.append(code('''NU = {"slab": 0, "cylinder": 1, "sphere": 2}


class TwoFieldPellet:
    """Species and energy as two coupled fields. The Prater relation is NOT used.

    Layout (n_u, 2): column 0 = y = c/c_s, column 1 = theta = T/T_s.
        grad^2 y     =  phi^2 * R(y, theta)
        grad^2 theta = -beta * phi^2 * R(y, theta)
    """

    def __init__(self, geom="sphere", n_u=200, bi_m=None, bi_h=None):
        self.geom, self.nu, self.n_u = geom, NU[geom], n_u
        self.shape = (n_u, 2)
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)          # faces
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])    # centres

        # --- boundary conditions, OUTWARD normal, one column per field -------
        # centre u=0: symmetry for both fields,  dy/dn = 0 and dtheta/dn = 0
        bc_l = {"a": np.array([[1.0, 1.0]]),
                "b": np.array([[0.0, 0.0]]),
                "d": np.array([[0.0, 0.0]])}
        # surface u=1, one column per field. Bi None -> Dirichlet, y = 1 / theta = 1.
        # With a film the flux INTO the pellet is supplied by the film, so with
        # the outward normal:
        #   D dc/dn = k_g (c_b - c_s)  ->  dy/dn     + Bi_m y     = Bi_m
        #   lam dT/dn = h  (T_b - T_s) ->  dtheta/dn + Bi_h theta = Bi_h
        # Both have the SAME sign pattern; writing the heat one with the inward
        # normal makes the pellet appear to absorb its own heat. The two fields
        # are set independently so that a MISMATCHED pair of boundary condition
        # types can be injected as a deliberate defect in section 4b.
        a_r = [0.0 if bi_m is None else 1.0, 0.0 if bi_h is None else 1.0]
        b_r = [1.0 if bi_m is None else bi_m, 1.0 if bi_h is None else bi_h]
        bc_r = {"a": np.array([a_r]), "b": np.array([b_r]), "d": np.array([b_r])}

        grad_mat, grad_bc = construct_grad(self.shape, self.u_f, self.u_c, (bc_l, bc_r))
        div_mat = construct_div(self.shape, self.u_f, nu=self.nu)  # nu: 0 slab, 1 cyl, 2 sphere
        self.lap = div_mat @ grad_mat
        self.lap_bc = (div_mat @ grad_bc).toarray().reshape(-1, 1)
        self.grad_mat = grad_mat
        self.grad_bc_v = grad_bc.toarray().reshape(-1, 1)
        s_mat, s_bc = construct_boundary_value_matrices(self.shape, self.u_f, self.u_c,
                                                        bc_r, axis=0, bound_id=1)
        self.surf_mat = s_mat
        self.surf_bc = np.asarray(s_bc.todense()).reshape(-1)
        self.numjac = NumJac(self.shape)          # default stencil: last-axis coupling
        # keep 1/theta finite during Newton; never active at a converged root
        self.lower = np.tile([0.0, 0.2], n_u).reshape(-1, 1)

    # ---------------------------------------------------------------- solving
    def solve(self, phi, beta, rate, v_init=None, tol=1e-12, maxfev=100, src_of=None):
        """`maxfev` and `src_of` exist only so that section 4b can stop Newton
        early and replace the source block by a deliberately defective one.
        `src_of(phi2, beta, rate)` returns the block, so it follows the ramp."""
        phi2 = phi * phi

        def src_correct(v):
            r = rate(v[..., 0], v[..., 1])
            return np.stack([-phi2 * r, beta * phi2 * r], axis=-1)

        src = src_correct if src_of is None else src_of(phi2, beta, rate)

        def residual(v):
            v = v.reshape(self.shape)
            g_s, jac_s = self.numjac(src, v)
            g = (self.lap @ v.reshape((-1, 1)) + self.lap_bc
                 + np.asarray(g_s).reshape((-1, 1)))
            return g, self.lap + jac_s

        v0 = np.ones(self.shape) if v_init is None else np.array(v_init, float).reshape(self.shape)
        res = newton(residual, v0.reshape((-1, 1)), maxfev=maxfev, tol=tol,
                     callback=lambda x, g: clip_approach(x, g, self.lower, None))
        v = res.x.reshape(self.shape)
        g, _ = residual(v)
        return v, float(np.max(np.abs(g)))

    # ------------------------------------------------------------ diagnostics
    def surface(self, v):
        """(y_s, theta_s) at u = 1 from pymrm's boundary-value operator."""
        s = (self.surf_mat @ v.reshape((-1, 1))).ravel() + self.surf_bc
        return float(s[-2]), float(s[-1])

    def eta(self, v, phi):
        """eta from the surface flux row: (nu+1)/phi^2 * dy/du at u = 1."""
        f = self.grad_mat @ v.reshape((-1, 1)) + self.grad_bc_v
        return (self.nu + 1.0) * float(f[-2].item()) / phi**2


def solve_ramp(pellet, phi, beta, rate, no_dead_core=True, **kw):
    """Deterministic continuation in beta from the isothermal solution.

    A fixed ladder of step counts, retried in order, so the path taken for a
    given case is identical on every machine - no warm-start chain along the
    swept parameter, whose convergence pattern would not be reproducible.

    `no_dead_core` additionally demands y > 0 everywhere. That matters only for
    zero-order kinetics, where the rate does not vanish with y and a dead core
    means there is no steady state at all; for the other rate laws y may
    underflow to 0 at the centre of a deeply ignited pellet, which is a
    perfectly good solution.

    Extra keywords are handed to `solve`, which is how section 4b ramps a
    defective source term.
    """
    for nstep in (8, 24, 72):
        v, rn = np.ones(pellet.shape), np.inf
        for f in np.linspace(0.0, 1.0, nstep + 1)[1:]:
            v, rn = pellet.solve(phi, f * beta, rate, v_init=v, **kw)
        ok = rn < 1e-8 and np.all(np.isfinite(v)) and v[:, 0].min() >= 0.0
        if ok and (not no_dead_core or v[:, 0].min() > 0.0):
            return v, rn, nstep
    return v, rn, -1


# --- the three rate laws, all pointwise in (y, theta) -----------------------
def r_first(gamma):
    """Weisz-Hicks kinetics written in theta: y exp[gamma(1 - 1/theta)]."""
    return lambda y, th: y * np.exp(gamma * (1.0 - 1.0 / th))


def r_lh(gamma, gamma_a=8.0, k_ads=3.0):
    """Langmuir-Hinshelwood, with a SECOND temperature dependence in K(theta)."""
    def rate(y, th):
        k = k_ads * np.exp(gamma_a * (1.0 / th - 1.0))
        return y * np.exp(gamma * (1.0 - 1.0 / th)) / (1.0 + k * y) ** 2
    return rate


def r_zero(gamma):
    """Zero order in the reactant; the rate does not vanish as y -> 0."""
    return lambda y, th: np.ones_like(y) * np.exp(gamma * (1.0 - 1.0 / th))


KIN = {"first order": r_first, "Langmuir-Hinshelwood": r_lh, "zero order": r_zero}
print("operators assembled per instance; nothing rebuilt inside Newton")'''))

# --------------------------------------------------------------- results: proof
cells.append(md(r"""## Results

### 1. The proof, in one line

Multiply the species balance by $\beta$ and add the energy balance. The rate
cancels *identically* — whatever function it is — leaving

$$\nabla^2\!\left(\theta + \beta y\right) = 0 .$$

So $w = \theta + \beta y$ is **harmonic**. It equals $1+\beta$ on the whole
boundary, because both fields are Dirichlet there with the same reference
values. A harmonic function equal to a constant on the boundary of a bounded
domain equals that constant everywhere (uniqueness of the Dirichlet problem, or
the maximum principle) — and that argument uses nothing about the shape of the
domain, which is exactly Prater's claim about "any particle geometry". Hence
$\theta + \beta y \equiv 1 + \beta$, which is the relation.

The cell below does this symbolically for a general geometry index $\nu$ and a
completely unspecified rate function."""))

cells.append(code('''u = sp.Symbol("u", positive=True)
nu_s, beta_s, phi_s, gamma_s = sp.symbols("nu beta phi gamma", real=True)
y_f, th_f, w_f = sp.Function("y"), sp.Function("theta"), sp.Function("w")
R_f = sp.Function("R")            # completely arbitrary kinetics


def laplacian(f):
    """(1/u^nu) d/du (u^nu df/du) - slab, cylinder, sphere and everything between."""
    return sp.diff(f, u, 2) + nu_s * sp.diff(f, u) / u


rate_sym = R_f(y_f(u), th_f(u))
res_species = laplacian(y_f(u)) - phi_s**2 * rate_sym            # = 0
res_energy = laplacian(th_f(u)) + beta_s * phi_s**2 * rate_sym   # = 0

combo = sp.simplify(res_energy + beta_s * res_species
                    - laplacian(th_f(u) + beta_s * y_f(u)))
print("energy + beta*species - laplacian(theta + beta*y)  =", combo)
assert combo == 0

# ...and a harmonic function with a symmetric centre is constant. dsolve's
# general form degenerates at nu = 1 (the cylinder), where the second solution
# is log u and not u^(1-nu), so that case is asked for separately.
gen = sp.dsolve(sp.Eq(laplacian(w_f(u)), 0), w_f(u)).rhs
print("general solution of laplacian(w) = 0, generic nu :", sp.simplify(gen))
for nu_val in (0, 1, 2):
    g_i = sp.dsolve(sp.Eq(laplacian(w_f(u)).subs(nu_s, nu_val), 0), w_f(u)).rhs
    print(f"   nu = {nu_val} ({['slab','cylinder','sphere'][nu_val]:8s})   "
          f"                    :", sp.simplify(g_i))
print("The second solution is singular at the centre for nu >= 1 (log u at")
print("nu = 1, u^(1-nu) for nu > 1) and has non-zero slope there for nu = 0,")
print("so regularity plus dw/du(0) = 0 leaves only the constant, and")
print("theta + beta*y = 1 + beta  everywhere, for ANY R and ANY nu.")'''))

cells.append(md(r"""### 2. Checking the transcription before building on it

A reading off a page image is a transcription and needs checking like any other.
Three checks, all using only the paper's own printed content."""))

cells.append(code('''# ---- check A: eq. (7) is exactly what eqs. (4) and (5) give ----------------
D, K, H, c0, T0, c, T, lap_c, lap_T, dndt = sp.symbols("D K H c_0 T_0 c T L_c L_T n", real=True)
eq4 = sp.Eq(D * lap_c - dndt, 0)                 # printed eq. (4)
eq5 = sp.Eq(K * lap_T - H * dndt, 0)             # printed eq. (5)
elim = sp.simplify(sp.solve([eq4, eq5], [dndt, lap_T], dict=True)[0][lap_T] - H * D * lap_c / K)
print("A. eliminating dn/dt from (4),(5):  K grad^2 T - H D grad^2 c = 0 ->", elim == 0)
# integrating twice with matched Dirichlet data: K(T-T0) = H D (c-c0)
dT_from_45 = H * D * (c - c0) / K
dT_printed = -H * D / K * (c0 - c)               # printed eq. (7)
print("   printed eq. (7) equals that integral      :",
      sp.simplify(dT_from_45 - dT_printed) == 0)

# ---- check B: the sign of beta, from eq. (7) at c -> 0 ---------------------
beta_from_eq7 = sp.simplify((dT_printed.subs(c, 0)) / T0)
beta_printed = c0 * H * D / (K * T0)             # middle member of eq. (8)
print(f"B. (dT/T0)_max from eq. (7)   : {beta_from_eq7}")
print(f"   beta as printed in eq. (8) : {beta_printed}")
print("   the two printed forms differ by a sign; the alternatives are")
print("   beta = +(dT/T0)_max  or  beta = -(dT/T0)_max.  Check C selects.")

# ---- check C: eq. (8)'s exponent is eq. (6) with dT/T0 = beta(1-y) ---------
y_s = sp.Symbol("y", positive=True)
dT_over_T0 = sp.Symbol("x", real=True)
exponent_eq6 = gamma_s * dT_over_T0 / (1 + dT_over_T0)          # from eq. (6)
exponent_eq8 = gamma_s * beta_s * (1 - y_s) / (1 + beta_s * (1 - y_s))  # printed eq. (8)
for name, sign in (("+", +1), ("-", -1)):
    got = sp.simplify(exponent_eq6.subs(dT_over_T0, sign * beta_s * (1 - y_s)) - exponent_eq8)
    print(f"C. beta = {name}(dT/T0)_max reproduces eq. (8)'s exponent : {got == 0}")
print()
print("Only beta = +(dT/T0)_max closes, which agrees with the paper sweeping")
print("beta from 0 to +0.8 for EXOTHERMIC reactions and with its central")
print("result eta >> 1. The page uses lam(T-T_s) = (-dH) D_e (c_s-c),")
print("beta = (-dH) D_e c_s / (lam T_s) > 0 exothermic.")'''))

cells.append(md(r"""### 3. What $\beta$ is worth in kelvin

$\beta$ is a transport group, so the maximum internal temperature rise can be
written down before any kinetics is known. The numbers below are **an
order-of-magnitude illustration with round values, taken from no source** — they
are here to show the size of the effect, not to characterise any real catalyst."""))

cells.append(code('''def prater_beta(dH_kJ_per_mol, D_e, c_s, lam, T_s):
    """beta = (-dH) D_e c_s / (lam T_s).  dH negative for an exothermic reaction."""
    return (-dH_kJ_per_mol * 1e3) * D_e * c_s / (lam * T_s)


print(f"{'(-dH) kJ/mol':>13} {'D_e m2/s':>10} {'c_s mol/m3':>11} "
      f"{'lam W/m/K':>10} {'T_s K':>7} {'beta':>9} {'dT_max K':>9}")
for dH, De, cs, lam, Ts in [(-100.0, 1e-6, 20.0, 0.3, 600.0),
                            (-100.0, 1e-6, 200.0, 0.3, 600.0),
                            (-800.0, 5e-6, 40.0, 0.5, 700.0),
                            (+180.0, 1e-6, 100.0, 0.3, 900.0)]:
    b = prater_beta(dH, De, cs, lam, Ts)
    print(f"{-dH:13.0f} {De:10.1e} {cs:11.0f} {lam:10.2f} {Ts:7.0f} "
          f"{b:9.4f} {b * Ts:9.1f}")
print()
print("A strongly exothermic oxidation can hold hundreds of kelvin inside a")
print("millimetre of catalyst; the last row is endothermic, so beta < 0.")'''))

# ------------------------------------------------- results: numerical sweep
cells.append(md(r"""### 4. Numerical confirmation, over everything the identity says is irrelevant

The sweep below solves the **two-field** system — never the reduced one — for
every combination of three geometries, three rate laws, three Arrhenius numbers
and a range of Prater numbers of both signs, and measures

$$\varepsilon = \frac{1}{|\beta|}\,\max_u\bigl|\,(\theta - 1) - \beta(1-y)\,\bigr| ,$$

the violation of the identity expressed as a fraction of $\beta$ itself. Every
case must converge; convergence is asserted separately, through the Newton
residual, because — as section 4b shows — $\varepsilon$ itself would not notice
if it had not.

**What $\varepsilon$ is.** The identity is *linear*, and one `construct_div`
serves both columns of the field array, so the discrete combination
$w = \theta + \beta y$ satisfies a closed linear subsystem: any root of the
discrete equations has $w \equiv 1+\beta$ whatever the grid, whatever the
geometry index, whatever the iterate it was reached from. $\varepsilon$ is
therefore a **floating-point measurement** — conditioning and roundoff — and not
a measure of the discretisation, the geometry, or the solve. Section 4b breaks
things on purpose to establish exactly which errors it does move, and by how
much.

Zero-order kinetics is restricted to small $\phi$: with no $y$ in the rate, a
dead core forms above a threshold and there is then *no* steady state to test —
which is physics, not a failure of the identity."""))

cells.append(code('''SWEEP = [("first order",           [0.3, 1.0, 2.0], (-0.4, -0.1, 0.1, 0.3, 0.6)),
         ("Langmuir-Hinshelwood",  [0.3, 1.0, 2.0], (-0.4, -0.1, 0.1, 0.3, 0.6)),
         ("zero order",            [0.2, 0.4],      (-0.4, -0.1, 0.1))]
GAMMAS = (10.0, 20.0, 30.0)

rows, failed = [], []
for kin, phis, betas in SWEEP:
    for geom in NU:
        for gamma in GAMMAS:
            for beta in betas:
                for phi in phis:
                    pel = TwoFieldPellet(geom, n_u=200)
                    v, rn, nstep = solve_ramp(pel, phi, beta, KIN[kin](gamma))
                    if nstep < 0:
                        failed.append((kin, geom, gamma, beta, phi, rn))
                        continue
                    y, th = v[:, 0], v[:, 1]
                    eps = np.max(np.abs((th - 1.0) - beta * (1.0 - y))) / abs(beta)
                    rows.append(dict(kin=kin, geom=geom, gamma=gamma, beta=beta, phi=phi,
                                     eps=eps, dT=np.abs(th - 1.0).max(), rn=rn,
                                     ymin=y.min(), eta=pel.eta(v, phi)))

assert not failed, f"cases that did not converge: {failed}"
eps_all = np.array([r["eps"] for r in rows])
worst = rows[int(np.argmax(eps_all))]
SWEEP_RESIDUAL = float(max(r["rn"] for r in rows))
print(f"{len(rows)} solves, all converged")
print(f"worst Newton residual     |g| = {SWEEP_RESIDUAL:.3e}   "
      f"(this, not eps, is what says the solves converged)")
print(f"worst identity violation  eps = {eps_all.max():.3e}   "
      f"(median {np.median(eps_all):.2e})")
print(f"  attained at: {worst['kin']}, {worst['geom']}, gamma={worst['gamma']:.0f}, "
      f"beta={worst['beta']:+.1f}, phi={worst['phi']}")
print(f"largest |theta-1| anywhere in the sweep : "
      f"{max(r['dT'] for r in rows):.4f}  (i.e. beta is reached)")
print(f"largest eta                             : {max(r['eta'] for r in rows):.1f}")
print(f"deepest depletion, min y                : {min(r['ymin'] for r in rows):.2e}")
SWEEP_WORST = float(eps_all.max())'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))

# left: profiles - theta plotted against the relation it is supposed to obey
demo = [("slab", 0), ("cylinder", 1), ("sphere", 2)]
BETA_D, GAMMA_D, PHI_D = 0.3, 20.0, 2.0
for geom, k in demo:
    pel = TwoFieldPellet(geom, n_u=400)
    v, rn, _ = solve_ramp(pel, PHI_D, BETA_D, KIN["first order"](GAMMA_D))
    c_ = plt.cm.viridis(0.1 + 0.35 * k)
    axes[0].plot(pel.u_c, v[:, 1], color=c_, lw=2, label=rf"{geom}: $\theta$ solved")
    axes[0].plot(pel.u_c[::20], 1 + BETA_D * (1 - v[::20, 0]), "o", color=c_,
                 ms=6, mfc="white", mew=1.5, label="_nolegend_")
axes[0].set_xlabel("$u$ [-]")
axes[0].set_ylabel(r"$\theta = T/T_s$ [-]")
axes[0].set_title(rf"Lines $\theta$ solved; markers $1+\beta(1-y)$"
                  rf" ($\beta$={BETA_D}, $\gamma$={GAMMA_D:.0f}, $\phi$={PHI_D})",
                  fontsize=10)
axes[0].legend(frameon=False, fontsize=9)

# right: the violation over the whole sweep
COLK = {"first order": "tab:blue", "Langmuir-Hinshelwood": "tab:orange",
        "zero order": "tab:green"}
MRK = {"slab": "o", "cylinder": "s", "sphere": "^"}
for kin in COLK:
    for geom in MRK:
        sel = [r for r in rows if r["kin"] == kin and r["geom"] == geom]
        if sel:
            axes[1].semilogy([abs(r["beta"]) for r in sel], [r["eps"] for r in sel],
                             MRK[geom], color=COLK[kin], ms=5, alpha=0.7,
                             mfc="none", label=f"{kin}, {geom}")
axes[1].axhline(SWEEP_WORST, color="k", ls=":", lw=0.9)
axes[1].set_xlabel(r"$|\beta|$ [-]")
axes[1].set_ylabel(r"$\varepsilon$, identity violation / $|\beta|$ [-]")
axes[1].set_ylim(1e-15, 3e-9)
axes[1].set_title(f"All {len(rows)} solves; nothing exceeds {SWEEP_WORST:.0e}",
                  fontsize=10)
axes[1].legend(frameon=False, fontsize=7, ncol=2, loc="upper center")
fig.tight_layout(); plt.show()'''))

# ------------------------------------------------- results: defect sensitivity
cells.append(md(r"""### 4b. Breaking it on purpose: what $\varepsilon$ can and cannot catch

An agreement number is worth exactly what it would cost to break. So one case
out of the 324 — the sphere at $\phi = 2$, $\beta = 0.3$, $\gamma = 20$ — is
solved again with **one implementation error injected at a time**, and
$\varepsilon$ is re-measured. Every number below is computed in the cells;
none is quoted.

The first table injects **inconsistencies between the two balances**: a sign, a
scale, a boundary condition type, a rate evaluated at the wrong state.
$\varepsilon$ moves by more than ten orders of magnitude, and it moves
*linearly* in the defect — a 1 % error in the heat source registers as
$\varepsilon = 10^{-2}$, which makes the number a measure of how big the mistake
is and not merely a flag.

The second table changes things that are wrong in the ordinary, expensive ways:
a Newton solve stopped before it converged, the wrong geometry index, a grid far
too coarse to resolve anything. **$\varepsilon$ does not move for any of them**,
because $w = \theta + \beta y$ is a closed linear subsystem of the discrete
equations: it is forced to $1+\beta$ at any root, on any grid, for any $\nu$.
One consequence is that the *severity* of the sweep — the extreme depletions and
efficiencies it reaches — adds no evidential weight to $\varepsilon$ at all. One
solve would have given the same number."""))

cells.append(code('''DPHI, DBETA, DGAMMA = 2.0, 0.3, 20.0          # one case out of the 324
d_rate = KIN["first order"](DGAMMA)


def eps_of(v, beta=DBETA):
    """Exactly the quantity the sweep measures."""
    y, th = v[:, 0], v[:, 1]
    return float(np.max(np.abs((th - 1.0) - beta * (1.0 - y))) / abs(beta))


def d_sign(phi2, beta, rate):
    """Heat released with the wrong sign: -beta where the energy balance wants +beta."""
    def src(v):
        r = rate(v[..., 0], v[..., 1])
        return np.stack([-phi2 * r, -beta * phi2 * r], axis=-1)
    return src


def d_scale(phi2, beta, rate):
    """A 1 % error in the heat of reaction - or in lambda, or in D_e."""
    def src(v):
        r = rate(v[..., 0], v[..., 1])
        return np.stack([-phi2 * r, 1.01 * beta * phi2 * r], axis=-1)
    return src


def d_state(phi2, beta, rate):
    """The energy balance evaluates the rate at theta = 1 instead of at theta."""
    def src(v):
        r = rate(v[..., 0], v[..., 1])
        r_cold = rate(v[..., 0], np.ones_like(v[..., 1]))
        return np.stack([-phi2 * r, beta * phi2 * r_cold], axis=-1)
    return src


DEFECTS = [("none - baseline", {}, {}),
           ("heat source uses -beta", {}, dict(src_of=d_sign)),
           ("heat source uses 1.01*beta", {}, dict(src_of=d_scale)),
           ("y Dirichlet but theta Robin, Bi=50", dict(bi_h=50.0), {}),
           ("rate at theta=1 in the energy eq", {}, dict(src_of=d_state))]

print(f"{'injected defect':>35} {'eps':>11} {'Newton |g|':>11} {'eta':>9}")
DEF = {}
for label, mkw, skw in DEFECTS:
    pel = TwoFieldPellet("sphere", n_u=200, **mkw)
    v, rn, ns = solve_ramp(pel, DPHI, DBETA, d_rate, no_dead_core=False, **skw)
    assert rn < 1e-8, (label, rn)     # every defective case still converges
    DEF[label] = eps_of(v)
    print(f"{label:>35} {DEF[label]:11.3e} {rn:11.2e} {pel.eta(v, DPHI):9.4f}")

DEF_BASE = DEF["none - baseline"]
DEF_SIGN = DEF["heat source uses -beta"]
DEF_SCALE = DEF["heat source uses 1.01*beta"]
DEF_BC = DEF["y Dirichlet but theta Robin, Bi=50"]
DEF_STATE = DEF["rate at theta=1 in the energy eq"]
print()
print(f"Every defect lifts eps from {DEF_BASE:.1e} to between {DEF_SCALE:.1e} and "
      f"{DEF_STATE:.1e},")
print("and each of these solves converged perfectly well - the Newton residual")
print("says nothing is wrong. A 1 % scale error registers as exactly "
      f"{DEF_SCALE:.1e}:")
print("the test is linear in the defect, so eps measures how big the mistake is.")'''))

cells.append(code('''print("blind (i): a Newton solve that never converged")
print(f"{'maxfev':>8} {'Newton |g|':>12} {'y at the centre':>16} {'eps':>11}")
blind_nc = []
for mf in (100, 3, 2, 1):
    pel = TwoFieldPellet("sphere", n_u=200)
    v, rn = pel.solve(DPHI, DBETA, d_rate, maxfev=mf, tol=1e-30)
    if mf < 100:
        blind_nc.append(eps_of(v))
    print(f"{mf:8d} {rn:12.3e} {v[0, 0]:16.6f} {eps_of(v):11.3e}")
BLIND_NOCONV = float(max(blind_nc))
print(f"The maxfev=1 iterate has y = {v[0, 0]:.2f} at the centre - more reactant")
print("inside the pellet than at its surface, which is impossible - on a residual")
print(f"of {rn:.1f}. Its eps is {eps_of(v):.1e}, inside the sweep's worst of "
      f"{SWEEP_WORST:.1e}.")
print()

print("blind (ii): the wrong geometry index. One construct_div serves both")
print("columns, so nu cancels out of the identity identically; a sphere solved")
print("with nu = 0 IS the slab operator, and returns a slab answer.")
print(f"{'nu':>4} {'geometry':>10} {'eta':>10} {'eps':>11}")
eta_g, blind_geom = {}, []
for g in ("sphere", "cylinder", "slab"):
    pel = TwoFieldPellet(g, n_u=200)
    v, rn, ns = solve_ramp(pel, DPHI, DBETA, d_rate)
    eta_g[g] = pel.eta(v, DPHI)
    blind_geom.append(eps_of(v))
    print(f"{pel.nu:4d} {g:>10} {eta_g[g]:10.4f} {eps_of(v):11.3e}")
BLIND_GEOM = float(max(blind_geom))
GEOM_ETA_ERR = abs(eta_g["slab"] - eta_g["sphere"]) / eta_g["sphere"]
print(f"Reporting the nu = 0 answer for a sphere is {GEOM_ETA_ERR:.0%} wrong in eta;")
print(f"eps stays below {BLIND_GEOM:.0e} for all three.")
print()

print("blind (iii): the grid. eps is compared with the error in eta, measured")
print("against the 200-cell solve of the same case.")
print(f"{'n_u':>6} {'eta':>10} {'error in eta':>14} {'eps':>11}")
ETA_REF, coarse = eta_g["sphere"], {}
for n in (200, 40, 12, 6, 3):
    pel = TwoFieldPellet("sphere", n_u=n)
    v, rn, ns = solve_ramp(pel, DPHI, DBETA, d_rate)
    e_n, eps_n = pel.eta(v, DPHI), eps_of(v)
    coarse[n] = (abs(e_n - ETA_REF) / ETA_REF, eps_n)
    print(f"{n:6d} {e_n:10.5f} {coarse[n][0]:13.1%} {eps_n:11.3e}")
BLIND_COARSE = float(coarse[3][1])
print(f"At n_u = 3 the efficiency is {coarse[3][0]:.0%} wrong while eps is "
      f"{BLIND_COARSE:.1e} -")
print(f"BETTER than the {coarse[200][1]:.1e} of the 200-cell solve. eps improves as")
print("the grid gets worse, because roundoff is all it ever measured.")
print()
print("So eps tests one thing: that the two balances see the same rate with the")
print("same sign and scale, and that the two fields carry the same kind of")
print("boundary condition. Convergence is established separately, by the Newton")
print("residual; the discretisation, by the shooting comparison in section 6.")'''))

# ------------------------------------------------- results: independent method
cells.append(md(r"""### 5. A second discretisation — what it settles and what it does not

The same problem is solved again with `scipy.solve_bvp`: collocation on an
adaptive mesh, a different method, a different solver, and no pymrm anywhere.
The continuation in $\beta$ is internal to this solver; the pymrm solution is
never used as a starting guess.

**This does not make the identity check independent.** `solve_bvp` puts the same
mesh and the same collocation nodes under every component of the first-order
system, and in that system $w = \theta + \beta y$ and $w' = \theta' + \beta y'$
again form a closed linear pair with the rate cancelled. The invariant is
inherited for the same structural reason as in Section 4, so $\varepsilon$
measured here carries no more weight than $\varepsilon$ measured there. The
table shows it plainly: several of these solves never reach their own tolerance
— `solve_bvp` exhausts `max_nodes` on the strongly depleted cases and returns
`status = 1` — and $\varepsilon$ does not distinguish them from the ones that
did.

**What the collocation run does contribute is the profiles.** Comparing $y(u)$
between the two solvers *is* a genuine test of the pymrm discretisation, and it
is one that can fail. The second cell makes that comparison under grid
refinement."""))

cells.append(code('''def bvp_pellet(nu, phi, beta, gamma, u0=1e-4, tol=1e-8, nstep=8, nodes=20000):
    """Independent collocation solve of the same two-field problem."""
    u_g = np.linspace(u0, 1.0, 300)
    w = np.vstack([np.ones_like(u_g), np.zeros_like(u_g),
                   np.ones_like(u_g), np.zeros_like(u_g)])
    sol = None
    for f in np.linspace(0.0, 1.0, nstep + 1)[1:]:
        b = f * beta

        def rhs(x, ww, b=b):
            yy, dy, tt, dt = ww
            r = yy * np.exp(gamma * (1.0 - 1.0 / tt))
            return np.vstack([dy, phi**2 * r - nu * dy / x,
                              dt, -b * phi**2 * r - nu * dt / x])

        def bcs(wa, wb):
            # dy/du = dtheta/du = 0 at the centre; y = theta = 1 at the surface
            return np.array([wa[1], wa[3], wb[0] - 1.0, wb[2] - 1.0])

        sol = solve_bvp(rhs, bcs, u_g, w, tol=tol, max_nodes=nodes)
        u_g, w = sol.x, sol.y
    uu = np.linspace(u0, 1.0, 4001)
    yy, tt = sol.sol(uu)[0], sol.sol(uu)[2]
    eps = np.max(np.abs((tt - 1.0) - beta * (1.0 - yy))) / abs(beta)
    return float(eps), float(yy.min()), float(tt.max()), sol


BVP_CASES = ((1.0, 0.3, 20.0), (2.0, -0.3, 20.0), (0.5, 0.6, 20.0), (3.0, 0.1, 30.0))
print(f"{'geometry':>10} {'phi':>5} {'beta':>6} {'gamma':>6} {'eps':>11} "
      f"{'y_min':>9} {'theta_max':>10} {'nodes':>7} {'rms res':>9} {'status':>7}")
bvp_eps, bvp_sols, bvp_bad = [], {}, []
for nu_i, gname in ((0, "slab"), (1, "cylinder"), (2, "sphere")):
    for phi, beta, gamma in BVP_CASES:
        e, ym, tm, sol = bvp_pellet(nu_i, phi, beta, gamma)
        bvp_eps.append(e)
        bvp_sols[(gname, phi, beta, gamma)] = sol
        rms = float(sol.rms_residuals.max())
        if sol.status != 0:
            bvp_bad.append(rms)
        print(f"{gname:>10} {phi:5.1f} {beta:+6.1f} {gamma:6.0f} {e:11.3e} "
              f"{ym:9.6f} {tm:10.6f} {sol.x.size:7d} {rms:9.1e} {sol.status:7d}")
BVP_WORST = float(max(bvp_eps))
BVP_NOT_CONVERGED = len(bvp_bad)
print(f"\\nworst eps over the collocation solves : {BVP_WORST:.3e}")
print(f"solves that did NOT reach their own tolerance (status != 0): "
      f"{BVP_NOT_CONVERGED} of {len(bvp_eps)},")
print(f"with rms residuals up to {max(bvp_bad):.1e}. Their eps is no larger than")
print("anyone else's: the collocation inherits the invariant just as pymrm does,")
print("so this table is not independent evidence for the identity.")'''))

cells.append(code('''# The profiles ARE an independent comparison of two discretisations.
print("pymrm vs collocation on y(u), under grid refinement")
print(f"{'case':>30} {'bvp status':>11} {'n_u=200':>10} {'n_u=400':>10} "
      f"{'n_u=800':>10} {'order':>7}")
prof_ok, prof_bad = [], []
for gname, phi, beta, gamma in (("slab", 1.0, 0.3, 20.0),
                                ("cylinder", 2.0, -0.3, 20.0),
                                ("sphere", 0.5, 0.6, 20.0),
                                ("sphere", 3.0, 0.1, 30.0)):
    sol = bvp_sols[(gname, phi, beta, gamma)]
    d = []
    for n_u in (200, 400, 800):
        pel = TwoFieldPellet(gname, n_u=n_u)
        v, rn, ns = solve_ramp(pel, phi, beta, KIN["first order"](gamma))
        d.append(float(np.max(np.abs(v[:, 0] - sol.sol(pel.u_c)[0]))))
    (prof_ok if sol.status == 0 else prof_bad).append(d[-1])
    lbl = "{} phi={:g} beta={:+.1f}".format(gname, phi, beta)
    print(f"{lbl:>30} {sol.status:11d} {d[0]:10.2e} {d[1]:10.2e} {d[2]:10.2e} "
          f"{np.log2(d[1] / d[2]):7.2f}")

PROFILE_MATCH = float(max(prof_ok))
PROFILE_STALLED = float(max(prof_bad))
print(f"\\nWhere the collocation solution is converged, pymrm approaches it at")
print(f"second order and agrees to {PROFILE_MATCH:.1e} at n_u = 800.")
print(f"Where it is not, the difference stalls at {PROFILE_STALLED:.1e} and does not")
print("fall under refinement - so it is not the pymrm discretisation error, and")
print("the status column says which of the two solutions to distrust. That is")
print("what a check with resolving power looks like.")'''))

# ------------------------------------------------- results: cross-page check
cells.append(md(r"""### 6. Cross-check against the published `B1.1`/`B1.5` page

The Weisz–Hicks page solves a **single** equation, obtained by substituting the
Prater relation into the species balance:

$$\nabla^2 y = \phi^2\,y\,\exp\!\left[\frac{\gamma\beta(1-y)}{1+\beta(1-y)}\right].$$

If the relation is exact, that reduced equation and this page's two-field system
must have the *same* solution, not merely a similar one. And its $\beta$ must be
this page's $\beta$: on `B1.1` it is called "the Prater number, the maximum
relative temperature rise", exactly the group defined here.

Three comparisons, at $\beta = 0.6$, $\gamma = 20$ — the parameters `B1.1` uses
for its multiplicity result — and they are **not of equal weight**, so each is
labelled for what it is:

1. $\beta$ **recovered** from the two-field solution as $(\theta-1)/(1-y)$,
   which should be the constant $0.6$ everywhere. This is the Section 4
   identity divided by $(1-y)$ and nothing more; it *cannot* fail, and its only
   content is that this page and `B1.1`/`B1.5` attach the same numerical value
   to the same symbol. A **convention check**, not corroboration.
2. The reduced equation solved on the same grid with the same operators. Given
   the exact discrete identity $\theta = 1 + \beta(1-y)$ and the algebraic
   identity $\gamma(1-1/\theta) \equiv \gamma\beta(1-y)/[1+\beta(1-y)]$, the
   two-field $y$ is an exact root of the reduced discrete system, so this
   agreement is guaranteed too. What it does test is that the two pages'
   $\gamma$ and $\beta$ enter the exponent the same way — an
   **implementation-consistency check** between pages.
3. The independent **shooting** reference from `B1.1`, which integrates the ODE
   and discretises nothing. This is the one comparison here that can fail, and
   the one that actually tests the pymrm discretisation: it must converge, and
   at second order."""))

cells.append(code('''class ReducedPellet:
    """The B1.1/B1.5 single equation: Prater already substituted."""

    def __init__(self, geom="sphere", n_u=400):
        self.nu, self.shape = NU[geom], (n_u,)
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])
        bc = ({"a": 1.0, "b": 0.0, "d": 0.0},    # centre: dy/dn = 0
              {"a": 0.0, "b": 1.0, "d": 1.0})    # surface: y = 1
        g, gb = construct_grad(self.shape, self.u_f, self.u_c, bc)
        d = construct_div(self.shape, self.u_f, nu=self.nu)   # 0/1/2 geometry
        self.lap, self.lap_bc = d @ g, (d @ gb).toarray().reshape(-1, 1)
        self.grad, self.grad_bc = g, gb.toarray().reshape(-1, 1)
        # (n_u, 1), NOT self.shape. NumJac's default stencil couples the LAST
        # axis in full; self.shape is (n_u,), whose last axis is SPACE, so
        # NumJac(self.shape) would declare every cell coupled to every other and
        # build a dense n_u x n_u Jacobian. `src` below is pointwise in y - the
        # tridiagonal coupling is in self.lap, which is added analytically - so
        # the block must be diagonal. Writing the shape with an explicit
        # length-1 field axis is the house layout and makes it so; the answer is
        # bit-identical, the Jacobian is n_u times cheaper to set up.
        self.numjac = NumJac((n_u, 1))

    def solve(self, phi, beta, gamma, y_init):
        p2 = phi * phi

        def src(y):
            yy = np.clip(y, 0.0, None)
            return -p2 * yy * np.exp(gamma * beta * (1 - yy) / (1 + beta * (1 - yy)))

        def residual(y):
            y = y.reshape(self.shape)
            g_s, j_s = self.numjac(src, y)
            return (self.lap @ y.reshape((-1, 1)) + self.lap_bc
                    + np.asarray(g_s).reshape((-1, 1))), self.lap + j_s

        r = newton(residual, np.asarray(y_init, float).reshape((-1, 1)), maxfev=60,
                   tol=1e-13, callback=lambda x, g: clip_approach(x, g, 0.0, None))
        return r.x.ravel()

    def eta(self, y, phi):
        f = self.grad @ y.reshape((-1, 1)) + self.grad_bc
        return (self.nu + 1.0) * float(f[-1].item()) / phi**2


def wh_shooting(beta, gamma, y_c):
    """B1.1's independent reference: shoot outward from the centre (sphere).

    Substituting s = phi*u removes phi, so the radius at which y reaches 1 IS
    phi. Reproduced here so the comparison uses no pymrm operator at all.
    """
    def rate(y):
        yy = np.clip(y, 0.0, None)
        return yy * np.exp(gamma * beta * (1 - yy) / (1 + beta * (1 - yy)))

    def rhs(s, v):
        return [v[1], rate(v[0]) - (2.0 / s) * v[1]]

    def hit(s, v):
        return v[0] - 1.0
    hit.terminal, hit.direction = True, 1
    s0 = 1e-6
    v0 = [y_c + rate(y_c) * s0**2 / 6.0, rate(y_c) * s0 / 3.0]
    sol = solve_ivp(rhs, (s0, 400.0), v0, events=hit, rtol=1e-10, atol=1e-12)
    if not sol.t_events[0].size:
        return None
    s_star = float(sol.t_events[0][0])
    return s_star, 3.0 * float(sol.y_events[0][0][1]) / s_star


BETA_X, GAMMA_X = 0.6, 20.0
print(f"{'geometry':>10} {'phi':>5} {'|beta-0.6|':>11} {'max|y2-yred|':>13} "
      f"{'eta 2-field':>12} {'eta reduced':>12}")
beta_err, red_err, eta_pairs = [], [], []
for geom in ("slab", "cylinder", "sphere"):
    for phi in (0.2, 0.4, 0.8, 2.0):
        pel = TwoFieldPellet(geom, n_u=400)
        v, rn, _ = solve_ramp(pel, phi, BETA_X, KIN["first order"](GAMMA_X))
        red = ReducedPellet(geom, n_u=400)
        y_red = red.solve(phi, BETA_X, GAMMA_X, v[:, 0])
        m = (1.0 - v[:, 0]) > 1e-10
        be = float(np.max(np.abs((v[m, 1] - 1.0) / (1.0 - v[m, 0]) - BETA_X)))
        dy = float(np.max(np.abs(v[:, 0] - y_red)))
        beta_err.append(be); red_err.append(dy)
        eta_pairs.append((pel.eta(v, phi), red.eta(y_red, phi)))
        print(f"{geom:>10} {phi:5.2f} {be:11.2e} {dy:13.3e} "
              f"{eta_pairs[-1][0]:12.5f} {eta_pairs[-1][1]:12.5f}")

BETA_RECOVERY = float(max(beta_err))
REDUCED_MATCH = float(max(red_err))
ETA_MATCH = float(max(abs(a - b) / b for a, b in eta_pairs))
print(f"\\nbeta recovered from the two-field solution : 0.6 to {BETA_RECOVERY:.2e}")
print("   (the identity restated - it cannot fail; what it checks is that both")
print("    pages mean the same number by the same symbol)")
print(f"two-field vs the B1.1/B1.5 reduced equation: {REDUCED_MATCH:.2e} in y, "
      f"{ETA_MATCH:.2e} relative in eta")
print("   (guaranteed by the exact discrete identity; what it checks is that")
print("    gamma and beta enter the two pages' exponents identically)")'''))

cells.append(code('''# The shooting reference discretises nothing; agreement with it must therefore
# improve under grid refinement, and at second order.
print(f"{'n_u':>6} {'phi':>5} {'eta 2-field':>12} {'eta shooting':>13} {'rel dev':>10}")
conv = {}
for n_u in (200, 400, 800):
    devs = []
    for phi in (0.2, 0.4, 0.8):
        pel = TwoFieldPellet("sphere", n_u=n_u)
        v, rn, _ = solve_ramp(pel, phi, BETA_X, KIN["first order"](GAMMA_X))
        ref = wh_shooting(BETA_X, GAMMA_X, float(v[0, 0]))
        e2, eref = pel.eta(v, phi), ref[1]
        devs.append(abs(e2 - eref) / eref)
        print(f"{n_u:6d} {phi:5.2f} {e2:12.5f} {eref:13.5f} {devs[-1]:10.2e}")
    conv[n_u] = max(devs)
SHOOT_DEV = float(conv[800])
order = np.log2(conv[400] / conv[800])
print(f"\\nworst deviation at n_u = 800: {SHOOT_DEV:.2e}; "
      f"observed order between 400 and 800: {order:.2f}")
print("The comparison is made at equal centre concentration, the parameter the")
print("two methods share, so no interpolation enters it.")
print("This is the page's real discretisation test: it can fail, and it falls at")
print(f"the expected rate. Its {SHOOT_DEV:.0e} carries more weight than any of the")
print("roundoff-level numbers above, all of which are inherited identities.")'''))

# ------------------------------------------------------- breakdown 1: film
cells.append(md(r"""### 7. Where it breaks (i): a film, with unequal Biot numbers

The proof needs $\theta + \beta y$ to be constant **on the boundary**. Two
Dirichlet conditions give that for free. A film does not.

With $Bi_m = k_gL/D_e$ and $Bi_h = hL/\lambda$, the steady external balances are
$Bi_m(1-y_s)$ moles in and $Bi_h(\theta_s-1)$ units of heat out, and these must
be in the same fixed proportion as inside, so

$$\theta_s - 1 = \beta\,\frac{Bi_m}{Bi_h}\,(1-y_s).$$

The **interior** relation still holds exactly, referenced to the *surface*
values — nothing inside the pellet changed. What fails is the relation people
actually write, referenced to the *bulk*. Subtracting the two gives a closed
form for the error, and it is independent of position:

$$\underbrace{(\theta - 1)}_{\text{true}} - \underbrace{\beta(1-y)}_{\text{bulk Prater}}
\;=\; \beta\,(1-y_s)\left(\frac{Bi_m}{Bi_h} - 1\right).$$

So the relation survives **if and only if $Bi_m = Bi_h$**, that is
$h/k_g = \lambda/D_e$. Nothing makes that happen:
$Bi_m/Bi_h = (k_g/h)\,(\lambda/D_e)$ multiplies a purely gas-side transfer ratio
by a purely pellet-side transport ratio, two numbers with no reason to be
reciprocal. The last column of the table below is the split it produces — how
much of the total temperature rise sits *outside* the pellet, in the film, where
the internal relation has no jurisdiction at all.

The cell below solves the film problem in pymrm and compares the measured
violation with the closed form above. Note the Robin conditions carry
$a = 1$ at the surface for both fields: same sign, outward normal.

**How to read the agreement with the closed form.** Both steps of the derivation
— that $\theta + \beta y$ is harmonic and that the surface balances fix its
boundary value — hold for the *discrete* system as well, so the numerical match
is forced and cannot fail. It confirms the algebra, not the solve; the second
cell shows the same match on an 8-cell grid that resolves nothing. Note also
that the closed form contains $y_s$, which is an output of the solve, so it is a
statement about the structure of the error rather than a prediction made in
advance. What is genuinely contributed here is the derivation itself, and the
magnitudes: how large the error is, and how much of the temperature rise moves
out of the pellet and into the film."""))

cells.append(code('''BETA_F, GAMMA_F, PHI_F = 0.05, 20.0, 3.0
FILM = [(100.0, 100.0), (100.0, 30.0), (100.0, 10.0), (100.0, 3.0), (10.0, 100.0)]

print(f"{'Bi_m':>7} {'Bi_h':>7} {'y_s':>9} {'theta_s':>9} {'theta_max':>10} "
      f"{'eps surface':>12} {'eps bulk':>11} {'closed form':>12} {'ext. dT':>9}")
film_rows, film_surf, film_pred_err = [], [], []
for bi_m, bi_h in FILM:
    pel = TwoFieldPellet("sphere", n_u=400, bi_m=bi_m, bi_h=bi_h)
    v, rn, ns = solve_ramp(pel, PHI_F, BETA_F, KIN["first order"](GAMMA_F),
                           no_dead_core=False)
    assert ns > 0, (bi_m, bi_h, rn)
    y, th = v[:, 0], v[:, 1]
    y_s, th_s = pel.surface(v)
    eps_surf = np.max(np.abs((th - th_s) - BETA_F * (y_s - y))) / BETA_F
    eps_bulk = np.max(np.abs((th - 1.0) - BETA_F * (1.0 - y))) / BETA_F
    pred = abs((1.0 - y_s) * (bi_m / bi_h - 1.0))
    film_rows.append((bi_m, bi_h, y_s, th_s, th.max(), eps_bulk, pred))
    film_surf.append(eps_surf)
    if pred > 1e-6:
        film_pred_err.append(abs(eps_bulk - pred) / pred)
    ext = (th_s - 1.0) / (th.max() - 1.0)      # fraction of the rise in the film
    print(f"{bi_m:7.0f} {bi_h:7.0f} {y_s:9.6f} {th_s:9.6f} {th.max():10.6f} "
          f"{eps_surf:12.2e} {eps_bulk:11.4e} {pred:12.4e} {ext:8.1%}")

FILM_SURF = float(max(film_surf))
FILM_PRED = float(max(film_pred_err))
print(f"\\nsurface-referenced Prater still exact to  : {FILM_SURF:.2e}")
print(f"bulk-referenced Prater vs the closed form : {FILM_PRED:.2e} relative")
print("The error is uniform in position, and its size is set by Bi_m/Bi_h alone.")'''))

cells.append(code('''# The match to the closed form is forced by the DISCRETE harmonicity, so it does
# not depend on resolution. Eight cells for phi = 3 resolve nothing at all, and
# the closed form is matched just as well - which is why this agreement confirms
# the algebra rather than the solve.
print(f"{'n_u':>6} {'y_s':>10} {'eps bulk':>12} {'closed form':>12} {'rel. diff':>11}")
grid_pred = []
for n in (8, 40, 400):
    pel = TwoFieldPellet("sphere", n_u=n, bi_m=100.0, bi_h=10.0)
    v, rn, ns = solve_ramp(pel, PHI_F, BETA_F, KIN["first order"](GAMMA_F),
                           no_dead_core=False)
    y_s, th_s = pel.surface(v)
    eb = np.max(np.abs((v[:, 1] - 1.0) - BETA_F * (1.0 - v[:, 0]))) / BETA_F
    pr = abs((1.0 - y_s) * (100.0 / 10.0 - 1.0))
    grid_pred.append(abs(eb - pr) / pr)
    print(f"{n:6d} {y_s:10.6f} {eb:12.4e} {pr:12.4e} {grid_pred[-1]:11.2e}")
FILM_PRED_COARSE = float(grid_pred[0])
print(f"\\nAt 8 cells the closed form is matched to {FILM_PRED_COARSE:.1e} - as well as")
print("at 400. The agreement is an algebraic identity of the discrete system, not")
print("a measure of how well the pellet was solved. The magnitudes in the table")
print("above (how much of the rise sits in the film) are the real content.")'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
show = [(100.0, 100.0), (100.0, 10.0), (100.0, 3.0)]
cols = plt.cm.plasma(np.linspace(0.1, 0.7, len(show)))
for (bi_m, bi_h), c_ in zip(show, cols):
    pel = TwoFieldPellet("sphere", n_u=400, bi_m=bi_m, bi_h=bi_h)
    v, rn, _ = solve_ramp(pel, PHI_F, BETA_F, KIN["first order"](GAMMA_F),
                          no_dead_core=False)
    y_s, th_s = pel.surface(v)
    lbl = rf"$Bi_m/Bi_h$ = {bi_m / bi_h:.0f}"
    axes[0].plot(pel.u_c, v[:, 1] - 1.0, color=c_, lw=2, label=lbl)
    axes[0].plot(pel.u_c, BETA_F * (1 - v[:, 0]), color=c_, lw=1.4, ls="--")
    axes[1].plot(pel.u_c, (v[:, 1] - th_s) - BETA_F * (y_s - v[:, 0]),
                 color=c_, lw=2, label=lbl)
axes[0].set_xlabel("$u$ [-]")
axes[0].set_ylabel(r"$\theta - 1$ [-]")
axes[0].set_yscale("log")
axes[0].set_title(r"Solid: true rise.  Dashed: what bulk Prater predicts",
                  fontsize=10)
axes[0].legend(frameon=False, fontsize=9)
axes[1].set_xlabel("$u$ [-]")
axes[1].set_ylabel(r"$(\theta-\theta_s) - \beta(y_s-y)$ [-]")
axes[1].set_ylim(-1e-11, 1e-11)
axes[1].set_title("Surface-referenced residual: zero to solver precision",
                  fontsize=10)
axes[1].legend(frameon=False, fontsize=9)
fig.tight_layout(); plt.show()'''))

# ------------------------------------------------------ breakdown 2: transient
cells.append(md(r"""### 8. Where it breaks (ii): a transient, unless $Le = 1$

The proof also needs steady state. With accumulation, scaling time on the
pellet's mass-diffusion time $\varepsilon L^2/D_e$,

$$\frac{\partial y}{\partial \tau} = \nabla^2 y - \phi^2\mathcal{R}, \qquad
\frac{\partial \theta}{\partial \tau} = Le\left[\nabla^2\theta + \beta\phi^2\mathcal{R}\right],
\qquad Le = \frac{\varepsilon\,\lambda}{\rho_p c_p D_e},$$

so that

$$\frac{\partial}{\partial\tau}(\theta + \beta y)
= \nabla^2(\theta + \beta y) + (Le - 1)\left[\nabla^2\theta + \beta\phi^2\mathcal{R}\right].$$

At $Le = 1$ the bracket disappears and $\theta + \beta y$ obeys a homogeneous
diffusion equation with constant data, so it stays at $1+\beta$ **at all times**:
the Prater relation is exact throughout the transient. That is a *linear
invariant* of the semi-discrete system, so any consistent integrator preserves
it exactly — which makes $Le = 1$ a control, and a stringent one.

Away from $Le = 1$ the relation is violated during the transient and recovers
only at steady state. The run below starts a pellet at bulk conditions
everywhere (so the identity holds at $\tau=0$) and lets it approach its steady
state."""))

cells.append(code('''BETA_T, GAMMA_T, PHI_T, N_T = 0.2, 20.0, 1.0, 100
pel_t = TwoFieldPellet("sphere", n_u=N_T)
rate_t = KIN["first order"](GAMMA_T)
PHI2_T = PHI_T**2


def transient(le, tau_end):
    """Startup from uniform bulk conditions. Operators assembled once, outside."""
    mult = np.tile([1.0, le], N_T).reshape(-1, 1)      # 1 on y, Le on theta
    M = diags_array(mult.ravel())

    def src(v):
        r = rate_t(v[..., 0], v[..., 1])
        return np.stack([-PHI2_T * r, le * BETA_T * PHI2_T * r], axis=-1)

    def rhs(t, w):
        s = src(w.reshape(pel_t.shape))
        return (mult * (pel_t.lap @ w.reshape((-1, 1)) + pel_t.lap_bc)
                + s.reshape((-1, 1))).ravel()

    def jac(t, w):
        _, j_s = pel_t.numjac(src, w.reshape(pel_t.shape))
        return M @ pel_t.lap + j_s

    t_eval = np.geomspace(1e-4, tau_end, 300)
    sol = solve_ivp(rhs, (0.0, tau_end), np.ones(2 * N_T), method="BDF", jac=jac,
                    rtol=1e-9, atol=1e-11, t_eval=t_eval)
    dev = np.array([np.max(np.abs((sol.y[:, k].reshape(pel_t.shape)[:, 1] - 1.0)
                                  - BETA_T * (1.0 - sol.y[:, k].reshape(pel_t.shape)[:, 0])))
                    / BETA_T for k in range(sol.y.shape[1])])
    return sol.t, dev


LEWIS = (0.01, 0.1, 1.0, 10.0)
fig, ax = plt.subplots(figsize=(7.8, 4.6))
cols = ["tab:blue", "tab:cyan", "k", "tab:red"]
tr = {}
print(f"{'Le':>7} {'max violation':>15} {'at tau':>9} {'at steady state':>17}")
for le, c_ in zip(LEWIS, cols):
    t, d = transient(le, 200.0 if le < 1.0 else 50.0)
    tr[le] = (t, d)
    ax.loglog(t, np.maximum(d, 1e-16), color=c_, lw=3 if le == 1.0 else 2,
              label=f"Le = {le:g}")
    print(f"{le:7.2f} {d.max():15.4e} {t[np.argmax(d)]:9.4f} {d[-1]:17.2e}")
ax.set_xlabel(r"$\\tau$, time on the mass-diffusion scale [-]")
ax.set_ylabel(r"$\\varepsilon$, identity violation / $\\beta$ [-]")
ax.set_title(r"Prater is exact at $Le=1$ at all times, and only there")
ax.legend(frameon=False)
fig.tight_layout(); plt.show()

TRANS_LE1 = float(tr[1.0][1].max())
TRANS_WORST = float(max(tr[le][1].max() for le in LEWIS if le != 1.0))
TRANS_PHYS = float(tr[0.1][1].max())
print(f"\\ncontrol, Le = 1 over the whole transient : {TRANS_LE1:.2e}")
print(f"worst violation away from Le = 1         : {TRANS_WORST:.4f} of beta")
print(f"at Le = 0.1, representative of a gas-filled porous pellet: "
      f"{TRANS_PHYS:.4f} of beta")
print("Every case returns to the identity as it approaches steady state;")
print("Le = 0.01 is the slowest and is still relaxing at the end of the window.")'''))

# ------------------------------------------------------------------ validation
cells.append(md(r"""## Validation

Seven items, and they are **not of equal weight**. Three of them are identities
of the discrete system that cannot fail; they are labelled as such where they
appear, and they are listed here the same way. None of them is a comparison with
a measurement — this page is **tier 6**, an identity proof with numerical
confirmation, and it must not be described as validated against experiment.

1. **The transcription.** Eq. (7) is exactly the integral of the printed eqs.
   (4) and (5); the printed $\beta$ in eq. (8) carries the opposite sign to the
   printed $(\Delta T/T_0)_{\max}$ next to it, and the paper's own rate
   expression, its stated sign convention for exothermic reactions and its
   $\eta > 1$ result all select the same reading.
2. **The identity, symbolically**, for arbitrary kinetics and arbitrary
   geometry index.
3. **The identity, numerically**, over the full sweep of geometries, rate laws
   and parameters — with convergence asserted separately, through the Newton
   residual, because $\varepsilon$ would not have noticed a solve that failed.
4. **The sensitivity of that test, measured** (Section 4b). Four injected
   defects lift $\varepsilon$ from roundoff to between $10^{-2}$ and $0.5$,
   linearly in the size of the defect; three ordinary errors — a non-converged
   solve, a wrong geometry index, a three-cell grid — leave it at roundoff. This
   is what fixes the scope of item 3.
5. **A second discretisation** (`scipy.solve_bvp` collocation) that shares no
   code with the first. It *inherits* the identity structurally, so it is not
   independent evidence for it; what it does give is an independent profile,
   and pymrm converges to that profile at second order.
6. **Cross-page.** $\beta$ recovered from the two-field solution and the match
   with `B1.1`/`B1.5`'s reduced equation are both guaranteed by the identity —
   they check conventions and implementation between pages, nothing more. The
   independent **shooting** reference, which discretises nothing, is the page's
   real discretisation test, and it converges at second order.
7. **The breakdowns are quantitative, not qualitative**: the film failure
   matches a closed form derived here (an algebraic identity of the discrete
   system, so it confirms the derivation rather than the solve), and the
   transient failure vanishes at $Le=1$ and is measured everywhere else."""))

cells.append(code('''print("1. transcription checks                    : symbolic, above (A, B, C)")
print("2. identity proved symbolically for any R, any nu")
print(f"3. worst violation over {len(rows)} pymrm solves      : {SWEEP_WORST:.3e}")
print(f"   worst Newton residual over those solves : {SWEEP_RESIDUAL:.3e}  "
      f"(this is what says they converged)")
print(f"4. eps under an injected defect            : sign {DEF_SIGN:.2e}, "
      f"1 % scale {DEF_SCALE:.2e},")
print(f"                                             BC-type mismatch {DEF_BC:.2e}, "
      f"wrong rate state {DEF_STATE:.2e}")
print(f"   eps where it is blind                   : unconverged solve "
      f"{BLIND_NOCONV:.1e}, wrong nu {BLIND_GEOM:.1e},")
print(f"                                             3-cell grid {BLIND_COARSE:.1e}")
print(f"5. worst violation, collocation (no pymrm) : {BVP_WORST:.3e}  "
      f"(inherited, and {BVP_NOT_CONVERGED} of those solves did not converge)")
print(f"   pymrm vs collocation profiles at n_u=800: {PROFILE_MATCH:.3e}  "
      f"(second order; this one can fail)")
print(f"6. beta recovered vs B1.1/B1.5 beta = 0.6  : {BETA_RECOVERY:.3e}  "
      f"(convention check; cannot fail)")
print(f"   two-field vs the reduced equation       : {REDUCED_MATCH:.3e} in y  "
      f"(guaranteed; cannot fail)")
print(f"   vs the independent shooting reference   : {SHOOT_DEV:.3e} at n_u = 800  "
      f"(the real discretisation test)")
print(f"7. film, surface-referenced still exact to : {FILM_SURF:.3e}")
print(f"   film, bulk-referenced vs closed form    : {FILM_PRED:.3e} relative  "
      f"({FILM_PRED_COARSE:.1e} on 8 cells - the match is algebraic)")
print(f"   transient control at Le = 1             : {TRANS_LE1:.3e}")
print(f"   transient worst away from Le = 1        : {TRANS_WORST:.4f} of beta")

report_agreement("B1.6", {
    # what the identity residual is, and what it is worth
    "sweep_worst_identity_violation": SWEEP_WORST,
    "sweep_worst_newton_residual": SWEEP_RESIDUAL,
    "defect_heat_source_sign_flip": DEF_SIGN,
    "defect_heat_source_1pct_scale": DEF_SCALE,
    "defect_bc_type_mismatch": DEF_BC,
    "defect_rate_at_inconsistent_state": DEF_STATE,
    "blind_eps_newton_not_converged": BLIND_NOCONV,
    "blind_eps_wrong_geometry_index": BLIND_GEOM,
    "blind_eps_three_cell_grid": BLIND_COARSE,
    # the checks that can fail
    "pymrm_vs_collocation_profile_n800": PROFILE_MATCH,
    "vs_shooting_reference_n800": SHOOT_DEV,
    # inherited identities, kept because they pin conventions between pages
    "collocation_worst_identity_violation": BVP_WORST,
    "two_field_vs_reduced_max_dy": REDUCED_MATCH,
    # the breakdowns
    "film_surface_referenced_residual": FILM_SURF,
    "film_bulk_error_vs_closed_form": FILM_PRED,
    "transient_control_Le1": TRANS_LE1,
    "transient_worst_violation": TRANS_WORST,
})'''))

# --------------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Nothing to the identity itself, and it would be dishonest to pretend
otherwise.** Prater's result is exact and its proof is two lines. Reproducing it
numerically confirms an implementation, not a theorem.

What pymrm adds is the part the 1958 result could not reach: **the size of the
error when the assumptions fail.** The original could only state the identity
under its hypotheses. Here each hypothesis is removed in turn and the
consequence is measured:

- **A film with $Bi_m \neq Bi_h$.** The interior relation survives, referenced
  to the surface. The bulk-referenced form — the one that gets written down —
  is wrong by $\beta(1-y_s)(Bi_m/Bi_h-1)$, uniformly in position. That closed
  form is derived here, not in either source. (The numerics agree with it to
  roundoff, but that agreement is forced: the derivation holds step for step in
  the discrete system, so it confirms the algebra and not the solve. What the
  computation supplies is the *magnitude* — how much of the temperature rise
  moves out of the pellet and into the film.) Since
  $Bi_m/Bi_h = (k_g/h)(\lambda/D_e)$ is not 1 for any real gas–solid catalyst,
  this is the case that matters in practice.
- **A transient with $Le \neq 1$.** The relation is a *linear invariant* of the
  time-dependent system exactly when $Le = 1$, and the computation confirms it
  to solver precision there while showing violations of tens of per cent of
  $\beta$ elsewhere, peaking early and decaying as the pellet reaches steady
  state. A quasi-steady pellet model that carries the Prater relation into a
  dynamic simulation is making an error of this size.
- **Geometry and kinetics as parameters, not rewrites.** `nu = 0, 1, 2` and a
  swapped rate function cover the whole sweep, so "any particle geometry, any
  kinetics" is exercised across nine combinations rather than argued.
  `construct_div` also accepts a callable `nu` for an arbitrary area profile;
  this page does not use one.

**And the identity is a cheap unit test for any non-isothermal pellet code —
for one specific class of bug.** Solve both balances, evaluate
$(\theta-1)-\beta(1-y)$, and a number that is not at roundoff means the two
balances disagree about the reaction: a wrong sign or a wrong scale on the heat
source, a rate evaluated at inconsistent states, or two fields carrying
different kinds of boundary condition. Section 4b measures how sharply it
responds — a 1 % scale error shows up as $\varepsilon \approx 10^{-2}$, and the
response is linear in the defect.

**It will not tell you anything else**, and a page recommending it has to say
so. $\varepsilon$ stays at roundoff for a Newton solve that never converged
(even one returning $y > 1$), for the wrong geometry index, and for a grid three
cells wide. Those need the Newton residual, a comparison against an independent
method, and grid refinement — all three of which are on this page precisely
because $\varepsilon$ cannot stand in for them."""))

# --------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

| Change | Where |
|---|---|
| Different geometry | `TwoFieldPellet(geom=...)`, or a callable `nu` in `construct_div` |
| Different kinetics | pass any `rate(y, theta)`; `NumJac` picks up the derivatives |
| External film | `TwoFieldPellet(..., bi_m=, bi_h=)`; both Robin rows have $a=1$ |
| Transient | `transient(le, tau_end)`; the operators are the steady ones |
| More species | widen the field axis to `(n_u, n_c+1)` and keep temperature last |
| Reversible or multiple reactions | the identity needs **one** shared rate; with several, one relation per independent reaction extent |

**The trap this page exists to flag.** The Prater relation is usually quoted
with the bulk concentration and bulk temperature, because that is what is
measurable. That version is only correct when the surface and bulk coincide, or
when $Bi_m = Bi_h$. Referenced to the surface it is exact; referenced to the
bulk it is not, and Section 7 gives the error in closed form.

**Reading Weisz and Hicks.** The scan's text layer is unusable for equations —
eq. (7) comes out as `AT = T - To = - F` and eq. (2) as `[grad ~1,~~`. Render
journal page 266 at 600 dpi and read it. The prose OCRs acceptably; the numbers
and symbols do not. This is the same trap recorded for `B3.1` and in
`docs/pdf-findings.md`.

**Related pages.** [`B1.1`/`B1.5`](../B1.1-thiele-weisz-hicks/) — the
non-isothermal pellet this relation makes solvable, where $\beta$ is a curve
label; `B1.2` Aris generalised modulus; `B1.4` Weisz–Prater criterion (the
observable-quantity counterpart of $\beta$); `D1.4` fixed bed with resolved
pellets, where the film conditions of Section 7 are supplied by the reactor;
`J4.7` immobilised enzyme particle, the same `S3` structure.

## References

Prater, C. D. (1958). *The temperature produced by heat of reaction in the
interior of porous particles.* Chemical Engineering Science **8**(3–4),
284–286. [doi:10.1016/0009-2509(58)85035-6](https://doi.org/10.1016/0009-2509(58)85035-6)
— **the origin of the result; not consulted for this page**, because no copy is
available. Cited as Weisz and Hicks cite it ("[11] PRATER C. D., *Chem. Engng.
Sci.* 1958 **8** 284"), with the DOI resolved and checked against that volume
and page.

Weisz, P. B. and Hicks, J. S. (1962). *The behaviour of porous catalyst
particles in view of internal mass and heat diffusion effects.* Chemical
Engineering Science **17**(4), 265–275.
[doi:10.1016/0009-2509(62)85005-2](https://doi.org/10.1016/0009-2509(62)85005-2)
— **the source actually read.** Eqs. (4), (5), (7) and (8) above are transcribed
from a 600 dpi render of journal page 266.

Damköhler, G. (1943). *Z. Phys. Chem.* **A193**, 16 — cited by Weisz and Hicks
as having obtained the relation for the spherical particle."""))

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
