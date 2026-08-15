#!/usr/bin/env python3
"""Generate index.ipynb for page I1.3. Run from the page directory."""
from __future__ import annotations

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# ----------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Heck, Wei & Katzer 1976: what the one-dimensional monolith model is adequate for"
description: "The paper's two claims, tested rather than illustrated - Table 1's asymptotic Nusselt numbers derived from scratch, the Nusselt singularity root-found and shown to sit 8.7 % downstream of the cause the text gives it, and the 1-D/2-D adequacy claim re-measured with the thresholds root-found instead of sampled."
categories: [sec:I, struct:S6, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-15
---

# Heck, Wei & Katzer 1976: what the one-dimensional monolith model is adequate for

**Catalog ID:** `I1.3` · **Structure:** `S6` (2-D PDE) · **Tier:** T0

The abstract makes two claims and the page exists to test both:

> "A two-dimensional model is shown to predict unusual behavior of the Nusselt
> number in the presence of rapid reaction. However, a simpler one-dimensional
> model is adequate for predicting monolith behavior."

*Unusual* compared with what — that is Table 1, the paper's only printed table,
and this page derives five of its eight cells from scratch (two of them in closed
form) rather than quoting them. *Adequate* for what, and where it stops — that is
a 1-D versus 2-D comparison the page can run itself, and the paper's own numbers
for it are **sampled inlet temperatures, not root-found thresholds**."""))

# ------------------------------------------------------------------- background
cells.append(md(r"""## Background

A monolith converter is a bundle of parallel ceramic channels with platinum on
the walls. CO oxidation over Pt is strongly exothermic *and* inhibited by CO, so
the rate is **negative order** in CO at high CO: the rate curve rises, peaks and
falls as the wall concentration increases. Put that curve against a straight
mass-transfer line and you get one or three intersections, and the jump from the
outer to the inner one is **light-off**.

The question the paper asks is a modelling question, not a chemistry one. A 1-D
model lumps the cross-section and closes the wall fluxes with a Nusselt and a
Sherwood number. A 2-D model resolves the radial profiles and *computes* the
Nusselt number from them. Young and Finlayson (1974) had reported that the
Nusselt number behaves strangely when the reaction is fast, and concluded the
expensive 2-D model was necessary. Heck, Wei and Katzer confirm the strange
behaviour and then argue the 1-D model is good enough anyway.

Both halves are checkable, and they are checkable **without touching a figure**,
because the paper states the numbers that matter in running text.

`I1.2` (Oh & Cavendish 1982) is the neighbouring page: the same reactor, the
same Voltz kinetics family, transient rather than steady, and reproduced against
a printed table. This page shares no data with it and loads none of its files."""))

# ------------------------------------------------------------- published model
cells.append(md(r"""## The published model

Read off the 300 ppi page bitmaps (`pdfimages -list` reports 300 ppi CCITT-G4
bilevel on all eight pages, so 300 is native and rendering higher only
interpolates). Nothing below came from the PDF text layer.

### The kinetics (p. 478), from Voltz et al. (1973)

$$r=\frac{k_r\,C\,C_{\mathrm{O_2}}}{[1+k_a C]^{2}},\qquad
k_r=k_r^{o}\exp\!\left[-\frac{E_r}{R\,(T_w+273)}\right],\qquad
k_a=k_a^{o}\exp\!\left[-\frac{E_a}{R\,(T_w+273)}\right]$$

with $k_r^{o}=4.14\times10^{8}$ kg mole/m²·s, $k_a^{o}=65.5$,
$E_r/R = 12{,}600$ °K and $E_a/R = -961$ °K. $C$ is a **mole fraction** (Notation,
p. 483). The paper writes no oxygen balance, so $C_{\mathrm{O_2}}$ is held at its
inlet value here and the alternative is a break-table row.

### The one-dimensional model (eqs. 8–14)

$$\frac{dC_G}{dX^{*}}=\frac{Sh}{Le}(C_w-C_G),\qquad
\frac{dT_G}{dX^{*}}=Nu\,(T_w-T_G),$$

closed by $-r\Delta H = h(T_w-T_G)$ and $r=\frac{\rho_G k_m}{M_G}(C_G-C_w)$, and
reduced with the adiabatic relation to a single algebraic wall condition. With
$Le\cong1$ and $Sh\cong Nu$ the paper's eq. (14) is

$$T_w \cong T_G^{o} + \Delta T_{\mathrm{AD}}\left(1-\frac{C_w}{C_G^{o}}\right),
\qquad \Delta T_{\mathrm{AD}}=\frac{\Delta H\,C_G^{o}}{M_G C_{PG}} .$$

**Equation (14) is not printed that way.** Verbatim, it reads

> $T_w \cong T_G{}^{o} + \Delta T_{AD} = \left(1 - \dfrac{C_w}{C_G{}^{o}}\right)$ 
> (14)

— an "=" where a multiplication is required. As printed it sets a temperature
equal to a dimensionless number, so it cannot be right, and the form above is what
eq. (13) gives when $Le\cong1$ and $Sh\cong Nu$ are imposed as the text instructs.
That is a derivation, not a guess, and it is labelled an inference. Reported, not
repaired.

### The two-dimensional model (eqs. 15–19)

$$\frac{(1-R^{2})}{2}\frac{\partial T_G}{\partial X^{*}}
=\frac{1}{R}\frac{\partial}{\partial R}\left(R\frac{\partial T_G}{\partial R}\right),
\qquad
\frac{(1-R^{2})}{2}\frac{\partial C_G}{\partial X^{*}}
=\frac{1}{Le\,R}\frac{\partial}{\partial R}\left(R\frac{\partial C_G}{\partial R}\right),$$

$$r\Delta H=\frac{2k_G}{D}\left.\frac{\partial T_G}{\partial R}\right|_w,
\qquad r=\frac{2\rho_G\mathcal{D}}{M_G D}\left.\frac{\partial C_G}{\partial R}\right|_w,
\qquad
Nu=\frac{2\,\partial T/\partial R|_w}{T_w-\bar T_G}. $$

Two things about this block are load-bearing for everything below.

**First, $X^{*}$ does not mean the same thing in eq. (9) and eq. (15).** The
Notation defines $X^{*}=4x/(D\,Re\,Pr)$. Equation (9) requires exactly that (it is
the 1-D energy balance with $Nu=hD/k_G$). Equation (15) requires
$X^{*}=x/(D\,Re\,Pr)$ — a factor of four smaller, and the figure abscissae are
labelled `X/(D RE PR)` without the four. The paper's own printed check value
settles which the 2-D results use, in §7a. A **third** scaling, $1000\,x/(D\,Re\,Pr)$,
is required by the Grigull–Tratz correlations the paper reprints on p. 482; §7d
shows that too, against a solve.

**Second, eqs. (17) and (18) cannot both have the sign they are printed with.**
The wall is a *sink* for CO and a *source* of heat, so $\partial C_G/\partial R|_w$
and $\partial T_G/\partial R|_w$ must have opposite signs; as printed both
equations set a positive $r$ proportional to a positive gradient. The same applies
to eq. (10) against eq. (12): $\Delta H$ cannot be positive in one and negative in
the other. These are sign-convention slips with an obvious intent, and they are
reported, not repaired: the physically consistent signs are used and stated."""))

# ------------------------------------------------------ parameters/assumptions
cells.append(md(r"""## Parameters and assumptions

Everything the paper prints is in `data/heck-1976-printed-values.csv`, with the
journal page and the location for each row: running text, a displayed equation, a
figure caption, or a character printed **inside** a figure's frame — an annotation
box, an abscissa tick label, a curve label or a case label. The data cell in *The
data* below enumerates every figure-sourced row and **counts them from the CSV
itself**, so the scope claim can be audited without diffing the file. **No curve was
digitised anywhere on this page.**

### What the paper does not print

The reacting model needs two dimensional groups the paper never states:

* $\Delta T_{\mathrm{AD}} = \Delta H\,C_G^{o}/(M_G C_{PG})$ — $\Delta H$ and
  $C_G^{o}$ are printed, $M_G C_{PG}$ is not;
* $\beta \equiv \rho_G k_m/M_G = (P/R_g T)\,Sh\,\mathcal{D}/D$ — $Sh$ and $D$
  are printed, the molar density and the diffusivity are not.

Neither can be recovered from the paper without circularity, so both are
**reconstructed from ordinary gas properties and stated as reconstructions**. The
important point for the page's conclusions is that the 1-D and the 2-D model use
*the same two numbers*, so the comparison between them — which is what the page is
about — is far less sensitive to the reconstruction than either model alone.
§7f measures exactly how much less.

### Scope: figures are out

Figures 1–13 carry the paper's profiles. Comparing against them would need
digitisation, which needs a maintainer review that is not available, so **it is
scoped out**: nothing is taken from a figure except characters *printed* in it, and
every such row is listed and counted below. What that costs is stated plainly in
*What this page cannot conclude*: the page cannot check any profile shape, any light-off *position*
against the paper's own, or the conversion curves of Figs. 4–6 and 12. What it
can do instead is check every number the paper states in words, and there are
more of those than the figure count suggests."""))

# ------------------------------------------------------------------- the data
cells.append(md(r"""## The data

Two files, both tier 6 (the paper's own printed values, no measurement).

* `heck-1976-table1-nusselt.csv` — Table 1, all eight cells. The paper's **only**
  printed table.
* `heck-1976-printed-values.csv` — every constant and every stated numerical
  result, with its journal page.

Table 1 is not the authors' computation: its footnote attributes it to Shah and
London (1971), which was **not consulted**. Two of its eight cells have exact
closed forms, derived here rather than quoted. Two cells cannot be checked from
this paper at all: the geometry column prints only a glyph, and for the
sinusoidal row the paper states no aspect ratio or wall profile, so that row is
under-determined by the source."""))

# --------------------------------------------------------------- cell 1: colab
cells.insert(1, code(r"""# Colab: install pymrm and fetch the shared helper if the repo is not checked out.
try:
    import pymrm
except ImportError:
    %pip install -q pymrm
import sys, pathlib
_here = pathlib.Path.cwd()
# published at pages/<id>-<slug>/ (repo root is parents[1]); staged one level
# deeper at queue_cases/<id>/page/ (parents[2]). Try both, then fall back to the
# raw URL for a Colab VM with no checkout.
_cands = [p / "shared" for p in list(_here.parents)[:3]] + [_here / "shared"]
for _p in _cands:
    if (_p / "gallery_utils.py").is_file():
        sys.path.insert(0, str(_p.parent))
        break
else:
    import urllib.request, os
    os.makedirs("shared", exist_ok=True)
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/computational-chemical-engineering/"
        "pymrm-gallery/main/shared/gallery_utils.py", "shared/gallery_utils.py")
    sys.path.insert(0, str(_here))
from shared.gallery_utils import load_data, load_meta, report_agreement
print("pymrm", getattr(pymrm, "__version__", "(dev)"))"""))

cells.append(code(r"""import numpy as np, pandas as pd, sympy as sp, time, warnings, re
np.set_printoptions(precision=6, suppress=False)
T1 = load_data("heck-1976-table1-nusselt.csv", page="I1.3-heck-wei-katzer-monolith")
PV = load_data("heck-1976-printed-values.csv", page="I1.3-heck-wei-katzer-monolith")
P = {r.quantity: r.value for r in PV.itertuples()}          # printed values, by name
def pf(name):                                                # printed value as float
    return float(P[name])
display(T1)
M = {}                                                        # every reported metric

from IPython.display import display_html
# Display a styled table as HTML ONLY. Two things about `display(df.style...)` are
# not reproducible and either would make two executions of this notebook differ byte
# for byte: the Styler's text/plain repr carries its MEMORY ADDRESS, and its HTML
# carries a RANDOM table uuid. This emits the HTML mime type alone, and every caller
# pins the uuid with `.set_uuid`.
def show(styler):
    display_html(styler.to_html(), raw=True)

# PROVENANCE AUDIT OF EVERY FIGURE-SOURCED ROW, generated from the CSV's own `where`
# column instead of being restated in prose, so that it cannot go stale when a row is
# added. Figures are scoped out on this page except for characters PRINTED inside
# their frames, and THIS is the list that promise has to be audited against - a
# maintainer checking the scope claim should not have to diff the CSV to find a row.
IN_FRAME = (("annotation", "annotation box"), ("tick label", "abscissa tick labels"),
            ("curve label", "curve label"), (r"case \d+ label", "case label"))
OUTSIDE  = (("caption", "figure caption"), ("running text", "running text"),
            ("displayed", "displayed equation"))
_kinds = lambda w, pats: ", ".join(lab for pat, lab in pats if re.search(pat, w)) or "-"
prov = pd.DataFrame([(r.quantity, r.where, _kinds(r.where, IN_FRAME),
                      _kinds(r.where, OUTSIDE))
                     for r in PV.itertuples() if re.search(r"Fig", r.where)],
                    columns=["row", "where (as recorded in the CSV)",
                             "read INSIDE the figure frame", "also printed outside it"])
show(prov.style.set_uuid("i13_prov").hide(axis="index"))
inside = prov["read INSIDE the figure frame"] != "-"
n_in = int(inside.sum())
n_only = int((inside & (prov["also printed outside it"] == "-")).sum())
print(f"{len(PV)} printed constants and stated results transcribed; "
      f"{len(T1)} rows in Table 1 (the paper's only printed table).")
print(f"{len(prov)} of those rows name a figure. {n_in} were read INSIDE a figure's "
      f"frame ({n_only} of them appear\nnowhere else in the paper, so the figure is "
      f"the only source); the remaining {len(prov)-n_in} are in\nrunning text or a "
      f"caption that merely names the figure.\nBy kind, counted from the CSV: "
      + ", ".join(f"{lab} {int(prov['read INSIDE the figure frame'].str.contains(lab, regex=False).sum())}"
                  for _, lab in IN_FRAME) + ".")
print("Every one is a printed CHARACTER. No curve is traced and no point is read off "
      "a plotted line.")"""))

# ------------------------------------------------------------ implementation
cells.append(md(r"""## PyMRM implementation

Three solvers, all built from `pymrm` operators and all reused by the validation
section.

1. **`duct(...)`** — fully developed flow and heat transfer in a duct
   cross-section, for Table 1. Cartesian 2-D for the square, cylindrical 1-D for
   the circle. `nu=1` in `construct_div` is the cylindrical geometry factor.
2. **`graetz(...)`** — the developing problem in a round tube, marched in the
   Graetz variable with `construct_convflux_upwind` + van Leer deferred
   correction, exactly as `A3.15` does. Used with an imposed wall flux (the
   paper's Figs. 8–9) and with a fixed wall temperature (its Fig. 10).
3. **`Radial` / `march_2d` / `march_1d`** — the reacting monolith channel, 2-D
   and 1-D, on a common axial coordinate.

Every switch a break-table row needs is a **function argument**, so a row cannot
silently perturb nothing: §8 asserts that each one changes the metric it claims
to move."""))

cells.append(code(r'''from scipy.sparse import diags, identity
from scipy.sparse.linalg import splu, eigsh
from scipy.optimize import brentq
from scipy.integrate import quad
from pymrm import (construct_grad, construct_div, construct_convflux_upwind,
                   compute_boundary_values, interp_cntr_to_stagg_tvd, vanleer)


def duct(nx, ny=None, geom="square", nu_r=1, plug=False, wall_neumann=False,
         weight_eigen=True):
    """Fully developed duct: Nu for constant wall flux (H1) and wall temperature (T).

    geom="square": quarter of a square of half-width 1, symmetry planes at x=y=0.
    geom="circle": radius 1, `nu_r=1` makes construct_div cylindrical.

    Solves lap(u) = -1 (u = 0 at the wall) for the velocity, then lap(th) = u for
    the H1 temperature, then the u-weighted eigenproblem for the T case.
    `plug`, `wall_neumann`, `nu_r` and `weight_eigen` exist for the break table.
    """
    if geom == "circle":
        f = np.linspace(0.0, 1.0, nx + 1); c = 0.5 * (f[:-1] + f[1:])
        shape, grids = (nx, 1), [(f, c, nu_r)]
        dv = (f[1:] ** 2 - f[:-1] ** 2) / 2.0                 # int r dr, 2*pi dropped
        A_full, P_full = np.pi, 2 * np.pi
    else:
        ny = ny or nx
        fx = np.linspace(0.0, 1.0, nx + 1); cx = 0.5 * (fx[:-1] + fx[1:])
        fy = np.linspace(0.0, 1.0, ny + 1); cy = 0.5 * (fy[:-1] + fy[1:])
        shape, grids = (nx, ny), [(fx, cx, 0), (fy, cy, 0)]
        dv = np.outer(np.diff(fx), np.diff(fy)).ravel()
        A_full, P_full = 4.0, 8.0                             # full square, half-width 1
    # a dphi/dn + b phi = d, n the OUTWARD normal.
    #   symmetry plane (lower) : dphi/dn = 0   -> a=1, b=0, d=0
    #   wall           (upper) : phi = 0       -> a=0, b=1, d=0   (Dirichlet)
    wall = {"a": 1.0, "b": 0.0, "d": 0.0} if wall_neumann else {"a": 0.0, "b": 1.0, "d": 0.0}
    bc = ({"a": 1.0, "b": 0.0, "d": 0.0}, wall)
    L, Lbc = 0, 0
    for ax, (f, c, nu) in enumerate(grids):
        g, gbc = construct_grad(shape, f, c, bc, axis=ax)
        d = construct_div(shape, f, nu=nu, axis=ax)            # nu=0 Cartesian, 1 cylindrical
        L = L + d @ g; Lbc = Lbc + d @ gbc
    Lbc = np.asarray(Lbc.todense()).ravel()
    lu = splu(L.tocsc())
    u = np.ones(int(np.prod(shape))) if plug else lu.solve(-1.0 - Lbc)
    th = lu.solve(u - Lbc)
    Aq = dv.sum()
    ub = (u * dv).sum() / Aq
    thb = (th * u * dv).sum() / (u * dv).sum()
    Nu_H1 = -4.0 * ub * A_full ** 2 / (P_full ** 2 * thb)
    w = dv * u if weight_eigen else dv
    # v0 fixed: ARPACK's default start vector is RANDOM, and two runs of this page
    # must produce byte-identical agreement.json.
    lam = eigsh((-diags(dv) @ L).tocsc(), k=1, M=diags(w).tocsc(),
                sigma=0.0, which="LM", v0=np.ones(int(np.prod(shape))))[0][0]
    Nu_T = 4.0 * lam * ub * A_full ** 2 / P_full ** 2
    fRe = (4.0 * A_full / P_full) ** 2 / (2.0 * ub)
    return dict(Nu_H1=Nu_H1, Nu_T=Nu_T, fRe=fRe, ub=ub)


def richardson(vals, ratio=2.0, order=2.0):
    """Richardson-extrapolate the last two rungs of a refinement ladder."""
    return vals[-1] + (vals[-1] - vals[-2]) / (ratio ** order - 1.0)


def observed_order(errs, ratio=2.0):
    return [float(np.log(errs[i] / errs[i + 1]) / np.log(ratio)) for i in range(len(errs) - 1)]


def fitted_order(ns, errs):
    """Least-squares slope of log(error) against log(n) over the whole ladder.

    Used where a rung-to-rung ratio is noisy - which it is on the axial ladder
    below, because the quantity being refined is LOCATED on that same axis.
    """
    return float(-np.polyfit(np.log(np.asarray(ns, float)),
                             np.log(np.asarray(errs, float)), 1)[0])'''))

cells.append(code(r'''class Graetz:
    """Developing temperature/concentration field in a round tube (Heck eqs. 15-16).

        half*(1 - R^2) dT/dX = (1/R) d/dR (R dT/dR)

    `half` is the coefficient AS PRINTED in eq. (15) (= 1/2). Passing half=2.0
    switches to the coefficient the Notation's X* = 4x/(D Re Pr) would require:
    that is the factor-of-four question, made a solver argument so section 7a can
    decide it by computation rather than by argument.
    """

    def __init__(self, X, n_x, n_r, half=0.5, x_grade=1.0, r_grade=2.0):
        t = np.linspace(0.0, 1.0, n_x + 1); self.x_f = X * t ** x_grade
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        t = np.linspace(0.0, 1.0, n_r + 1); self.r_f = 1.0 - (1.0 - t) ** r_grade
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        self.shape = (n_x, n_r)
        self.v = np.broadcast_to(half * (1.0 - self.r_c ** 2), (n_x + 1, n_r))
        aw = 2.0 * self.r_c * np.diff(self.r_f)
        self.area_w = aw
        self.cup_w = aw * (1.0 - self.r_c ** 2) / np.sum(aw * (1.0 - self.r_c ** 2))

    def solve(self, wall, inlet=0.0, tvd=True, tol=1e-10, maxit=600, omega=0.6,
              accept=1e-8):
        """`wall` is either ("flux", Qw(X)) or ("temp", value)."""
        n_x, n_r = self.shape
        # a dT/dn + b T = d, n the OUTWARD normal.
        #   inlet  X=0 : T = inlet            -> a=0, b=1, d=inlet   (Dirichlet)
        #   outlet X=X : dT/dX = 0            -> a=1, b=0, d=0       (n = +X there)
        bc_x = ({"a": 0.0, "b": 1.0, "d": inlet}, {"a": 1.0, "b": 0.0, "d": 0.0})
        #   axis  R=0 : dT/dR = 0, symmetry   -> a=1, b=0, d=0
        if wall[0] == "flux":
            q = np.asarray(wall[1](self.x_c), float).reshape(n_x, 1)
            #   wall R=1 : dT/dR = Qw(X)      -> a=1, b=0, d=Qw   (n = +R there)
            bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": q})
        else:
            q = None
            #   wall R=1 : T = value          -> a=0, b=1, d=value (Dirichlet)
            bc_r = ({"a": 1.0, "b": 0.0, "d": 0.0},
                    {"a": 0.0, "b": 1.0, "d": float(wall[1])})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                  bc_x, v=self.v, axis=0)
        dx = construct_div(self.shape, self.x_f, nu=0, axis=0)      # nu=0: Cartesian, axial
        gr, gr_bc = construct_grad(self.shape, self.r_f, self.r_c, bc_r, axis=1)
        dr = construct_div(self.shape, self.r_f, nu=1, axis=1)      # nu=1: cylindrical, radial
        A = (dx @ conv - dr @ gr).tocsc()
        b = np.asarray((dx @ conv_bc - dr @ gr_bc).todense()).ravel()
        lu = splu(A)
        T = lu.solve(-b).reshape(self.shape)
        nit, incr = 0, 0.0
        if tvd:                       # van Leer deferred correction, under-relaxed
            for nit in range(1, maxit + 1):
                _, dT = interp_cntr_to_stagg_tvd(T, self.x_f, self.x_c, bc_x,
                                                 self.v, vanleer, axis=0)
                corr = np.asarray(dx @ (self.v * dT).reshape(-1, 1)).ravel()
                new = T + omega * (lu.solve(-b - corr).reshape(self.shape) - T)
                incr = np.max(np.abs(new - T)) / omega
                T = new
                if incr < tol:
                    break
            else:
                # raise rather than return a half-converged field - but the bar is
                # `accept`, three orders below the coarsest discretisation error on
                # this page, not the tolerance the loop aims for.
                if incr > accept:
                    raise RuntimeError(f"deferred correction did not converge: {incr:.2e}")
        # the wall value is read with pymrm's own boundary reconstruction, not by
        # taking the last cell: the same call handles the Neumann and the
        # Dirichlet wall, and for the Dirichlet wall it must return the bc value.
        Tw, _ = compute_boundary_values(T, self.r_f, self.r_c, bc_r[1],
                                        axis=1, bound_id=1)
        return dict(x=self.x_c, r=self.r_c, T=T, Tw=np.asarray(Tw).ravel(),
                    cup=T @ self.cup_w, area=T @ self.area_w,
                    q=(None if q is None else q.ravel()), nit=nit, incr=incr)


def local_nusselt(sol):
    """Heck eq. (19): Nu = 2 (dT/dR)|_w / (T_w - Tbar), with dT/dR|_w the imposed flux."""
    return 2.0 * sol["q"] / (sol["Tw"] - sol["cup"])


def cross(x, y):
    """Root of y(x) at its LAST sign change, by interpolation - never a sampled index."""
    i = np.where(np.diff(np.sign(y)) != 0)[0]
    if len(i) == 0:
        return np.nan
    j = i[-1]
    return brentq(lambda t: np.interp(t, x, y), x[j], x[j + 1])'''))


cells.append(code(r'''# ---- the reacting monolith channel: kinetics, then the 1-D and 2-D models -----
KR0, KA0, ERR, EAR = pf("k_r0"), pf("k_a0"), pf("Er_over_R"), pf("Ea_over_R")
CG0, CO20 = pf("CG_CO_inlet"), pf("CG_O2_inlet")


def rate(Cw, Tw, o2="const", kr0=None, ka0=None, err=None, ear=None):
    """Voltz kinetics, eq. (5)-(7). Cw, Tw may be complex (complex-step derivative)."""
    kr0 = KR0 if kr0 is None else kr0
    ka0 = KA0 if ka0 is None else ka0
    err = ERR if err is None else err
    ear = EAR if ear is None else ear
    TK = Tw + 273.0
    kr = kr0 * np.exp(-err / TK)
    ka = ka0 * np.exp(-ear / TK)
    if o2 == "const":
        CO2 = CO20
    else:                                   # stoichiometric O2 depletion, 1/2 per CO
        CO2 = CO20 - 0.5 * (CG0 - Cw)
        CO2 = np.where(np.real(CO2) > 0.0, CO2, 0.0)
    return kr * Cw * CO2 / (1.0 + ka * Cw) ** 2


def r_of(Cw, TG0, dTad, **kw):
    """Rate on the adiabatic locus, eq. (14) as derived from eq. (13)."""
    return rate(Cw, TG0 + dTad * (1.0 - Cw / CG0), **kw)


def drdc(Cw, TG0, dTad, **kw):
    """dr/dCw by complex step - exact to machine precision, no step-size choice."""
    return np.imag(r_of(np.asarray(Cw, dtype=complex) + 1e-30j, TG0, dTad, **kw)) / 1e-30


def fold_cw(TG0, dTad, beta, **kw):
    """Wall concentration at the fold of the cold branch, ROOT-FOUND from r'(c) = -beta.

    F(c) = r(c) + beta c is the total sink; the cold branch is the largest root of
    F(c) = beta C_G, and it merges with the middle root at the local minimum of F.
    NaN when no fold exists - which is exactly the paper's own p. 480 condition,
    "when (-rho_G k_m / M_G) is greater than the maximum slope of (-dr/dC_w), there
    is no tangent point".
    """
    cs = np.geomspace(1e-8, CG0, 600)
    h = drdc(cs, TG0, dTad, **kw) + beta
    idx = np.where(np.diff(np.sign(h)) != 0)[0]
    if len(idx) == 0:
        return np.nan
    return brentq(lambda c: drdc(c, TG0, dTad, **kw) + beta,
                  cs[idx[-1]], cs[idx[-1] + 1], xtol=1e-18, rtol=8.9e-16)


def fold_CG(TG0, dTad, beta, **kw):
    """Bulk C_G at which the cold branch folds (light-off). NaN if no fold exists."""
    c2 = fold_cw(TG0, dTad, beta, **kw)
    if not np.isfinite(c2):
        return np.nan
    return (r_of(c2, TG0, dTad, **kw) + beta * c2) / beta


def Cw_cold(CG, TG0, dTad, beta, **kw):
    return brentq(lambda c: r_of(c, TG0, dTad, **kw) - beta * (CG - c),
                  1e-14, CG, xtol=1e-18, rtol=8.9e-16)


def lightoff_1d_quad(TG0, dTad, beta, Sh, npts=120, **kw):
    """Light-off position of the constant-Nu 1-D model, by quadrature in C_w.

    dC_G/dG = -4 Sh (C_G - C_w) with C_G = C_w + r(C_w)/beta on the cold branch, so

        G* = int_{c_fold}^{c_inlet}  [beta + r'(c)] / [4 Sh r(c)]  dc .

    Parametrising by C_w rather than by C_G matters: dC_G/dC_w vanishes at the
    fold, so the C_G form has a square-root endpoint and scipy's adaptive `quad`
    is 0.7 % out on it (measured in section 7e). This form has no singularity at
    all - the integrand VANISHES at the fold, because r'(c_fold) = -beta - and it
    is converged at 100 Gauss points. It shares no discretisation with march_1d.
    """
    cf = fold_cw(TG0, dTad, beta, **kw)
    if not np.isfinite(cf):
        return np.inf, np.nan
    CGs = (r_of(cf, TG0, dTad, **kw) + beta * cf) / beta
    if CGs >= CG0:
        return 0.0, CGs
    cin = Cw_cold(CG0, TG0, dTad, beta, **kw)
    x, w = np.polynomial.legendre.leggauss(npts)
    c = 0.5 * (cin - cf) * x + 0.5 * (cin + cf)
    f = (beta + drdc(c, TG0, dTad, **kw)) / (4.0 * Sh * r_of(c, TG0, dTad, **kw))
    return float(0.5 * (cin - cf) * np.sum(w * f)), CGs


def lightoff_1d_quad_CG(TG0, dTad, beta, Sh, limit=600, **kw):
    """The same integral parametrised by C_G - kept only to measure the trap."""
    CGs = fold_CG(TG0, dTad, beta, **kw)
    if not np.isfinite(CGs) or CGs >= CG0:
        return np.inf if not np.isfinite(CGs) else 0.0
    f = lambda CG: 1.0 / (4.0 * Sh * (CG - Cw_cold(CG, TG0, dTad, beta, **kw)))
    return quad(f, CGs * (1.0 + 1e-12), CG0, limit=limit)[0]'''))

cells.append(code(r'''# Grigull and Tratz (1965) as the paper reprints them on p. 482, verbatim.
NuQ_GT = lambda Xs: 4.364 + 8.68 * Xs ** -0.506 * np.exp(-0.041 * Xs)   # constant wall flux
NuT_GT = lambda Xs: 3.655 + 6.874 * Xs ** -0.488 * np.exp(-0.0572 * Xs)  # constant wall temp


def march_1d(TG0, dTad, b0, NuFun, dG, Gmax, o2="const", **kw):
    """1-D monolith on G = x/(D Re Pr).  dc/dG = -4 Nu(G) (c - c_w).

    b0 = beta/Sh is the property group; beta(G) = b0 Nu(G), so a variable Nusselt
    number moves the mass-transfer line as well as the heat balance, which is what
    the paper's Fig. 12 model does. Returns the light-off G (inf if none).
    """
    c, G, cw_prev = 1.0, 0.0, 1.0
    hist = []
    while G < Gmax - 1e-15:
        Nu = NuFun(G); beta = b0 * Nu
        Phi = lambda cw: c - r_of(cw * CG0, TG0, dTad, o2=o2, **kw) / (beta * CG0) - cw
        grid = np.linspace(1e-12, c, 500); v = Phi(grid)
        idx = np.where(np.diff(np.sign(v)) != 0)[0]
        if len(idx) == 0:
            return G, hist
        cw = brentq(Phi, grid[idx[-1]], grid[idx[-1] + 1], xtol=1e-16, rtol=8.9e-16)
        if cw < 0.5 * cw_prev and cw_prev - cw > 0.2 * cw_prev:
            return G, hist
        hist.append((G, c, cw))
        c = c - dG * 4.0 * Nu * (c - cw)
        cw_prev = cw; G += dG
    return np.inf, hist


class Radial:
    """Radial operators for the 2-D reacting channel, assembled once."""

    def __init__(self, n_r=48, grade=2.0, nu_r=1):
        t = np.linspace(0.0, 1.0, n_r + 1)
        self.r_f = 1.0 - (1.0 - t) ** grade
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        self.shape = (n_r, 1)
        # a dc/dn + b c = d, n the OUTWARD normal.
        #   R=0 axis : dc/dR = 0, symmetry  -> a=1, b=0, d=0
        #   R=1 wall : dc/dR = q, imposed   -> a=1, b=0, d=q  (n = +R there)
        self.bc0 = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        g, _ = construct_grad(self.shape, self.r_f, self.r_c, self.bc0, axis=0)
        d = construct_div(self.shape, self.r_f, nu=nu_r, axis=0)   # nu=1: cylindrical
        self.L = (d @ g).tocsc()
        bc1 = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 1.0})
        _, gbc1 = construct_grad(self.shape, self.r_f, self.r_c, bc1, axis=0)
        self.Lq = np.asarray((d @ gbc1).todense()).ravel()         # rhs per unit wall flux
        self.w = 0.5 * (1.0 - self.r_c ** 2)                       # eq. (15) coefficient
        aw = 2.0 * self.r_c * np.diff(self.r_f)
        self.cup_w = aw * (1.0 - self.r_c ** 2) / np.sum(aw * (1.0 - self.r_c ** 2))

    def wall_value(self, c, q):
        bc = {"a": 1.0, "b": 0.0, "d": q}
        return float(compute_boundary_values(c.reshape(self.shape), self.r_f,
                                             self.r_c, bc, axis=0,
                                             bound_id=1)[0].ravel()[0])


def march_2d(TG0, dTad, beta, Sh, dG, Gmax, n_r=48, o2="const", grade=2.0,
             nu_r=1, record=False, **kw):
    """2-D reacting channel, eqs. (15)-(18), reduced to one field by Le = 1.

    With Le = 1 the wall gradients of the scaled temperature and of (1 - c) are
    equal and both fields start uniform, so tau = 1 - c everywhere and eq. (14)
    holds pointwise at the wall with the 2-D wall concentration. The field is
    linear in the wall flux q, so each implicit step needs exactly two linear
    solves (q = 0 and q = 1) and then a scalar root-find for the wall value -
    which is the 2-D analogue of Fig. 1's mass-transfer line.
    """
    rad = Radial(n_r, grade=grade, nu_r=nu_r)
    gam = 2.0 * beta * CG0 / Sh                    # eq. (18) group; independent of Sh
    A = (identity(n_r, format="csc").multiply(rad.w[:, None] / dG) - rad.L).tocsc()
    lu = splu(A)
    c = np.ones(n_r); G = 0.0; cw_prev = 1.0; hist = []
    while G < Gmax - 1e-15:
        rhs0 = rad.w * c / dG
        c0 = lu.solve(rhs0); c1 = lu.solve(rhs0 + rad.Lq)
        w0 = rad.wall_value(c0, 0.0); S = rad.wall_value(c1, 1.0) - w0
        Phi = lambda cw: w0 - S * r_of(cw * CG0, TG0, dTad, o2=o2, **kw) / gam - cw
        grid = np.linspace(1e-12, max(w0, 1e-9), 500); v = Phi(grid)
        idx = np.where(np.diff(np.sign(v)) != 0)[0]
        if len(idx) == 0:
            return G, hist
        cw = brentq(Phi, grid[idx[-1]], grid[idx[-1] + 1], xtol=1e-16, rtol=8.9e-16)
        if cw < 0.5 * cw_prev and cw_prev - cw > 0.2 * cw_prev:
            return G, hist
        q = -r_of(cw * CG0, TG0, dTad, o2=o2, **kw) / gam
        c = c0 + (c1 - c0) * q
        G += dG; cw_prev = cw
        if record:
            cup = float(c @ rad.cup_w)
            hist.append((G, cw, cup, -2.0 * q / (cw - cup)))
    return np.inf, hist'''))

# --------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. Table 1, derived rather than quoted

The paper's only table is the "usual" behaviour against which the Nusselt number
is said to be unusual, so it is worth more than a citation. Two of its cells have
exact closed forms; the square pair does not, and is computed twice on
discretisations that share nothing."""))

cells.append(code(r'''# --- exact closed forms, by symbolic algebra (no grid anywhere) ---------------
r_, R0_ = sp.symbols("r R0", positive=True)
u_c = (R0_ ** 2 - r_ ** 2) / 4                                  # lap u = -1, u(R0) = 0
th_c = sp.integrate(sp.integrate(u_c * r_, r_) / r_, r_)
th_c = sp.simplify(th_c - th_c.subs(r_, R0_))                   # lap th = u, th(R0) = 0
I_c = lambda f: sp.integrate(f * 2 * sp.pi * r_, (r_, 0, R0_))
A_c, P_c = sp.pi * R0_ ** 2, 2 * sp.pi * R0_
Nu_circle_exact = sp.nsimplify(sp.simplify(
    -4 * (I_c(u_c) / A_c) * A_c ** 2 / (P_c ** 2 * (I_c(th_c * u_c) / I_c(u_c)))))

x_, y_ = sp.symbols("x y", real=True); s3 = sp.sqrt(3)
l1, l2, l3 = y_ + 1, (-y_ + s3 * x_) / 2 + 1, (-y_ - s3 * x_) / 2 + 1   # incircle radius 1
lap = lambda f: sp.diff(f, x_, 2) + sp.diff(f, y_, 2)
cc = sp.symbols("cc")
cc_val = sp.solve(sp.Poly(sp.expand(lap(cc * l1 * l2 * l3) + 1), x_, y_).coeffs(), cc)[cc]
u_t = sp.expand(cc_val * l1 * l2 * l3)
aa, bb = sp.symbols("aa bb")
th_t = sp.expand(l1 * l2 * l3 * (aa + bb * (x_ ** 2 + y_ ** 2)))
sol_t = sp.solve(sp.Poly(sp.expand(lap(th_t) - u_t), x_, y_).coeffs(), [aa, bb])
th_t = sp.expand(th_t.subs(sol_t))
assert sp.simplify(lap(th_t) - u_t) == 0 and sp.simplify(th_t.subs(y_, -1)) == 0
xl, xr = sp.solve(sp.Eq(l3, 0), x_)[0], sp.solve(sp.Eq(l2, 0), x_)[0]
I_t = lambda f: sp.integrate(sp.integrate(f, (x_, xr, xl)), (y_, -1, 2))
A_t, P_t = I_t(1), 6 * sp.sqrt(3)
Nu_tri_exact = sp.nsimplify(sp.simplify(
    -4 * (I_t(u_t) / A_t) * A_t ** 2 / (P_t ** 2 * (I_t(th_t * u_t) / I_t(u_t)))))

print(f"circular duct,   constant wall flux : Nu_H1 = {Nu_circle_exact} = "
      f"{float(Nu_circle_exact):.6f}   (Table 1 prints {T1.Nu_H[0]})")
print(f"equilateral tri, constant wall flux : Nu_H1 = {Nu_tri_exact} = "
      f"{float(Nu_tri_exact):.6f}   (Table 1 prints {T1.Nu_H[2]})")
print("\nBoth are derived here: the velocity is the exact Poisson solution on the")
print("cross-section (for the triangle, the product of the three wall distances),")
print("and the H1 temperature is an exact polynomial. Nothing is quoted.")
M["nu_h1_circle_exact"] = float(Nu_circle_exact)
M["nu_h1_triangle_exact"] = float(Nu_tri_exact)
M["nu_h1_circle_printed_relerr"] = abs(float(Nu_circle_exact) - float(T1.Nu_H[0])) / float(Nu_circle_exact)
M["nu_h1_triangle_printed_relerr"] = abs(float(Nu_tri_exact) - float(T1.Nu_H[2])) / float(Nu_tri_exact)'''))

cells.append(code(r'''# --- the same four numbers from pymrm finite volumes --------------------------
t0 = time.time()
lad_c = [duct(n, geom="circle") for n in (64, 128, 256, 512)]
lad_s = [duct(n, geom="square") for n in (32, 64, 128, 256)]
circ = lad_c[-1]; sq = lad_s[-1]
NuH1_c = richardson([d["Nu_H1"] for d in lad_c]); NuT_c = richardson([d["Nu_T"] for d in lad_c])
NuH1_s = richardson([d["Nu_H1"] for d in lad_s]); NuT_s = richardson([d["Nu_T"] for d in lad_s])
print(f"circle   fRe -> {richardson([d['fRe'] for d in lad_c]):.6f}   (exact 16)")
print(f"square   fRe -> {richardson([d['fRe'] for d in lad_s]):.6f}   (14.227 is the "
      f"textbook square-duct value; nothing on this page uses it, it checks the velocity field)")
rows = []
rows.append(("circular", "Nu_H", float(T1.Nu_H[0]), NuH1_c, float(Nu_circle_exact)))
rows.append(("circular", "Nu_T", float(T1.Nu_T[0]), NuT_c, np.nan))
rows.append(("square",   "Nu_H", float(T1.Nu_H[1]), NuH1_s, np.nan))
rows.append(("square",   "Nu_T", float(T1.Nu_T[1]), NuT_s, np.nan))
rows.append(("equilateral triangle", "Nu_H", float(T1.Nu_H[2]), np.nan, float(Nu_tri_exact)))
tab1 = pd.DataFrame(rows, columns=["geometry", "column", "printed", "pymrm (extrapolated)", "closed form"])
tab1["rel. dev."] = np.where(np.isfinite(tab1["pymrm (extrapolated)"]),
                             abs(tab1["printed"] - tab1["pymrm (extrapolated)"]) / tab1["pymrm (extrapolated)"],
                             abs(tab1["printed"] - tab1["closed form"]) / tab1["closed form"])
show(tab1.style.set_uuid("i13_table1").format({"printed": "{:.4f}", "pymrm (extrapolated)": "{:.6f}",
                           "closed form": "{:.6f}", "rel. dev.": "{:.2e}"}))
M["nu_h1_square_pymrm"] = NuH1_s
M["nu_t_square_pymrm"] = NuT_s
M["nu_t_circle_pymrm"] = NuT_c
M["nu_h1_circle_pymrm"] = NuH1_c
M["nu_h1_circle_solve_vs_exact"] = abs(NuH1_c - float(Nu_circle_exact)) / float(Nu_circle_exact)
M["nu_h1_square_printed_relerr"] = abs(float(T1.Nu_H[1]) - NuH1_s) / NuH1_s
M["nu_t_square_printed_relerr"] = abs(float(T1.Nu_T[1]) - NuT_s) / NuT_s
M["nu_t_circle_printed_relerr"] = abs(float(T1.Nu_T[0]) - NuT_c) / NuT_c
print(f"\n{time.time()-t0:.1f} s")'''))

cells.append(code(r'''# --- an independent route for the square: Fourier-Galerkin, no pymrm, no grid --
def square_spectral(N, nq=160):
    """Sine-series Galerkin on (0,2)^2. Shares no code with duct(): different basis,
    different quadrature, different eigen-formulation."""
    m = np.arange(1, N + 1); k = m * np.pi / 2.0
    b1 = 2.0 * (1.0 - np.cos(m * np.pi)) / (m * np.pi)
    K2 = k[:, None] ** 2 + k[None, :] ** 2
    U = np.outer(b1, b1) / K2                      # -lap u = 1
    ints = (1.0 - np.cos(m * np.pi)) / k
    IU = (U * np.outer(ints, ints)).sum(); ub = IU / 4.0
    TH = -U / K2                                   # lap th = u
    thb = (TH * U).sum() / IU
    A, Pm = 4.0, 8.0
    NuH1 = -4 * ub * A ** 2 / (Pm ** 2 * thb)
    xg, wg = np.polynomial.legendre.leggauss(nq); xg = xg + 1.0
    S = np.sin(np.outer(k, xg))
    ug = np.einsum("ab,ai,bj->ij", U, S, S, optimize=True)
    Wq = np.outer(wg, wg) * ug
    Amp = np.einsum("mi,pi->mpi", S, S).reshape(N * N, nq)
    # (Amp Wq Amp^T) is indexed ((m,p),(n,q)); the eigenproblem needs ((m,n),(p,q))
    Mm = np.transpose((Amp @ Wq @ Amp.T).reshape(N, N, N, N),
                      (0, 2, 1, 3)).reshape(N * N, N * N)
    import scipy.linalg as sl
    lam0 = sl.eigh(np.diag(K2.ravel()), Mm, eigvals_only=True,
                   subset_by_index=[0, 0])[0]
    return NuH1, 4 * lam0 * ub * A ** 2 / Pm ** 2

sp_lad = [square_spectral(N) for N in (12, 16, 20, 24)]
print("Fourier-Galerkin square duct (converges from below, slowly - the sine series")
print("of a constant is what limits it):")
for N, (a, b) in zip((12, 16, 20, 24), sp_lad):
    print(f"  N={N:3d}  Nu_H1={a:.6f}  Nu_T={b:.6f}")
print(f"\npymrm finite volume, extrapolated:  Nu_H1={NuH1_s:.6f}  Nu_T={NuT_s:.6f}")
print(f"""
NOT A BOUND, and an earlier draft of this page called it one. The sine-Galerkin
eigenproblem is a Rayleigh-Ritz restriction, so it bounds the eigenvalue from
ABOVE, while the truncated integral of u converges from BELOW; the product
4 mu ub A^2/P^2 carries no one-sided guarantee, and four monotone samples are not a
proof of monotonicity. What the sequence does show is that it is still increasing at
N = 24, where it has already passed {sp_lad[-1][1]:.6f}, above the interval [{float(T1.Nu_T[1])-5e-4:.4f}, {float(T1.Nu_T[1])+5e-4:.4f}]
that the printed {float(T1.Nu_T[1]):.4f} allows. The finite-volume extrapolation independently gives
{NuT_s:.7f}, and the two routes agree to {abs(sp_lad[-1][0]-NuH1_s)/NuH1_s:.1e} relative on Nu_H1.""")
M["nu_t_square_spectral_N24"] = float(sp_lad[-1][1])
M["nu_h1_square_two_routes_reldiff"] = abs(sp_lad[-1][0] - NuH1_s) / NuH1_s'''))

cells.append(md(r"""**Five of the eight cells are reproduced.** Four agree with the printed digits.
The fifth does not: the square duct at constant wall temperature is printed
**2.976** where both routes give **2.9775**. `2.976` implies the interval
[2.9755, 2.9765], which excludes 2.9775, so this is not a rounding difference. It
is also not a defect this page can pin on Heck et al.: Table 1 is attributed to
Shah and London (1971), which was **not consulted**, so the discrepancy could be
in the 1971 tabulation or in the transcription of it, and the page does not
choose. (The paper's own reference list prints that source as "**Shar**, R. K.,
and A. L. London" where Table 1's footnote says "Shah and London" — reported,
not repaired.)

**Three cells are not reproduced**, and the reasons differ. The equilateral
triangle at constant wall temperature needs a $u$-weighted eigenproblem on a
non-tensor domain, which this page does not build. Both sinusoidal cells are
**under-determined by the paper**: the geometry column prints a glyph and no
aspect ratio or wall profile anywhere, so there is no problem to solve. Saying
so is the honest outcome, not a gap to paper over."""))

# ---------------------------------------------------------- the Nusselt anomaly
cells.append(md(r"""### 2. The unusual Nusselt number — and what actually causes it

The paper demonstrates the anomaly on the *non-reacting* Graetz problem with an
imposed wall flux that rises through zero, and it prints both the flux functions
and the landmark positions. Everything needed is in running text, in Fig. 8's own
curve labels and in Fig. 9's caption — printed characters, no curve read:

* Fig. 8: $Q_w = -9.46 + (46.0X^{*})^{8}$
* Fig. 9: $Q_w = -3 + (40X^{*})^{8}$
* "the wall flux rises rapidly from a negative value and passes through zero near
  $X^{*}=0.029$. As the wall flux passes through zero, infinite negative and
  positive values for the Nusselt number occur."
* "When $X^{*}=0.032$ … the average gas temperature $T_G^{*}$ is higher than
  $T_w^{*}$; thus a large negative value results … At $X^{*}=0.033$,
  $(\partial T_G^{*}/\partial R)_w$ and $(T_w^{*}-\bar T_G^{*})$ are both
  positive."

So the paper states one location for the singularity in the first quotation and
brackets a different one in the second. That is decidable by computation."""))

cells.append(code(r'''t0 = time.time()
Q8 = lambda X: -9.46 + (46.0 * X) ** 8
Q9 = lambda X: -3.0 + (40.0 * X) ** 8
gz = Graetz(X=0.05, n_x=6000, n_r=200, r_grade=1.0)
res = {}
for nm, Q in (("Fig. 8", Q8), ("Fig. 9", Q9)):
    s = gz.solve(("flux", Q)); res[nm] = s
    zero_q = brentq(Q, 1e-4, 0.05)                              # exact, from the formula
    pole_cup = cross(s["x"], s["Tw"] - s["cup"])
    pole_area = cross(s["x"], s["Tw"] - s["area"])
    Nu = local_nusselt(s)
    pre = s["x"] < pole_cup                     # before the pole: Nu is continuous there
    zero_Nu = cross(s["x"][pre], Nu[pre])
    Nu_at_zero_q = float(np.interp(zero_q, s["x"], Nu))
    print(f"{nm}  Qw = 0 at X* = {zero_q:.6f}   (text says 'near 0.029')")
    print(f"          there Nu = {Nu_at_zero_q:+.3e}, i.e. it passes through ZERO, not infinity")
    print(f"          Nu = 0   at X* = {zero_Nu:.6f}   (the same point, root-found on Nu itself)")
    print(f"          Nu pole  at X* = {pole_cup:.6f} (mixing-cup mean), "
          f"{pole_area:.6f} (area mean)")
    print(f"          pole sits {100*(pole_cup/zero_q-1):.2f} % downstream of the flux zero\n")
    res[nm + "/marks"] = (zero_q, zero_Nu, pole_cup, pole_area)
zq, zN, pc, pa = res["Fig. 9/marks"]

# The paper's own TEN radial points, solved ONCE here and reused by section 7c, so
# that the two places which quote it quote the same computed number.
s10r = Graetz(X=0.05, n_x=3200, n_r=10, r_grade=1.0).solve(("flux", Q9))
pole_10r = cross(s10r["x"], s10r["Tw"] - s10r["cup"])

# NORMALISATION INVARIANCE, COMPUTED RATHER THAN ARGUED. The paper never states the
# magnitude of Q_w. The problem is linear in T, so halving the flux must halve the
# field and leave the zero of (Tw - Tbar) exactly where it is - which is why the
# landmark positions above do not depend on a convention the paper omits. Run twice
# on one grid rather than asserted from linearity:
gh = Graetz(X=0.05, n_x=3200, n_r=40, r_grade=1.0)
sh1 = gh.solve(("flux", Q9)); sh2 = gh.solve(("flux", lambda X: 0.5 * Q9(X)))
ph1 = cross(sh1["x"], sh1["Tw"] - sh1["cup"])
ph2 = cross(sh2["x"], sh2["Tw"] - sh2["cup"])
print(f"Q_w normalisation: pole at X* = {ph1:.12f} with Q_w, {ph2:.12f} with Q_w/2")
print(f"                   difference {abs(ph1-ph2):.2e} - the invariance is solved for, "
      f"not assumed\n")

M["fluxzero_X_fig9"] = zq
M["nu_zero_X_fig9"] = zN
M["nu_zero_minus_fluxzero_fig9"] = abs(zN - zq)
M["nu_pole_X_fig9"] = pc
M["nu_pole_over_fluxzero_fig9"] = pc / zq
M["nu_pole_X_fig8"] = res["Fig. 8/marks"][2]
print(f"{time.time()-t0:.1f} s")'''))

cells.append(code(r'''# The paper's own sign statements at X* = 0.025, 0.032, 0.033 (Fig. 9's flux).
s = res["Fig. 9"]
print("  X*      Qw        dT/dR|w   Tw - Tbar    Nu        paper says")
says = {0.025: "cooling, Nu positive and ordinary",
        0.032: "dT/dR|w > 0 but Tbar > Tw  ->  LARGE NEGATIVE Nu",
        0.033: "dT/dR|w > 0 and Tw > Tbar  ->  large positive Nu"}
ok = {}
for xq in (0.025, 0.032, 0.033):
    i = int(np.argmin(abs(s["x"] - xq)))
    d = s["Tw"][i] - s["cup"][i]; q = s["q"][i]
    print(f"{s['x'][i]:7.4f} {q:+9.4f} {q:+10.4f} {d:+11.5f} {2*q/d:+10.3f}   {says[xq]}")
    ok[xq] = (np.sign(q), np.sign(d))
M["fig9_sign_0p032_matches_paper"] = float(ok[0.032] == (1.0, -1.0))
M["fig9_sign_0p033_matches_paper"] = float(ok[0.033] == (1.0, 1.0))
print(f"""
The paper's X* = 0.033 statement is reproduced; its X* = 0.032 statement is not.
The pole is root-found at X* = {pc:.6f}, so by X* = 0.032 the solve has already
passed it and (Tw - Tbar) is positive. The paper's own two sentences bracket the
pole in (0.032, 0.033); this solve puts it {100*(1-pc/0.0325):.1f} % below the centre of that
bracket. The bracket and the pole disagree by more than any grid effect: on the
paper's own ten radial points the same solver gives X* = {pole_10r:.6f}, which is
{100*abs(1-pole_10r/pc):.2f} % from the converged value and still {100*(1-pole_10r/0.032):.1f} % short of 0.032. §7c refines the
radial axis from those 10 points to 320 and quotes this same number.""")'''))

cells.append(md(r"""**The cause the text gives is not the cause.** At the zero of the wall flux the
Nusselt number of eq. (19) has a vanishing *numerator*, so it passes through
**zero**, not infinity — which is exactly what the computation shows and what
Fig. 8's own case-2 curve does on its way down. The pole is the zero of the
*denominator* $T_w-\bar T_G$, and it lies **8.7 % further downstream** in both of
the paper's flux functions. The two are close but distinct, and they must be:
the bulk temperature integrates the flux history, so it keeps falling for a while
after the flux turns positive.

This is not a gap the authors missed. The very next paragraph — the Fig. 9
discussion quoted above — gives the correct mechanism explicitly, in terms of the
signs of the gradient and of $(T_w^{*}-\bar T_G^{*})$. What is wrong is the
one-sentence summary, and the numerical bracket that paragraph implies for the
pole, which this solve puts 4 % lower.

**The anomaly itself is entirely real and reproduces.** It is also
*normalisation-invariant*: the problem is linear in $T$, so scaling $Q_w$ scales
the field and leaves the pole where it is. Whatever convention the paper used for
the magnitude of $Q_w$ — which it never states — the landmark positions above are
unchanged."""))

# ------------------------------------------------- claim 2: 1-D versus 2-D
cells.append(md(r"""### 3. "A simpler one-dimensional model is adequate" — for what, and until where

The reacting comparison needs two groups the paper never prints, so they are
reconstructed from ordinary gas properties. **Nothing below is fitted to any
result of the paper.**

An earlier draft of this page said that "two printed facts then check the
reconstruction". They do not, and the cell below now measures how little they
check. One of the two is a **transcription check on the kinetics** — worth having,
but $k_a^{o}$ and $E_a/R$ were printed, not reconstructed — and it constrains
$\Delta T_{\mathrm{AD}}$ so weakly that the cell root-finds the whole range of
$\Delta T_{\mathrm{AD}}$ that passes it. The other has **no printed referent at
all**: the paper never states an adiabatic flame temperature in words anywhere in
pp. 477–483 (checked page by page on the 300 ppi bitmaps; the ~750 °C plateau
exists only as a plotted curve in Figs. 2, 3, 7 and 13, and figures are scoped
out), so it is a consequence of the reconstruction, not a test of it. And
$\beta$ is checked by neither. §7f is where the reconstruction is actually
held to account, by measuring what it moves."""))

cells.append(code(r'''# ---- reconstruction of the two groups the paper does not print ---------------
P_atm, Rg = 101325.0, 8314.0        # 1 atm; universal gas constant, J/kmol.K
M_G, C_PG = 29.0, 1093.0            # exhaust gas: kg/kmol and J/kg.K, ordinary values
D_ch = pf("channel_diameter") * 1e-3
Sh0 = pf("Nu_used_1D")              # 3.608, the paper's own square-channel value
T_props = pf("TG_inlet_fig1") + 273.0
D_AB = 0.208e-4 * (T_props / 300.0) ** 1.75      # CO in air, 0.208 cm2/s at 300 K
rho_over_M = P_atm / (Rg * T_props)
dTad = pf("delta_H") * 1e6 * CG0 / (M_G * C_PG)
b0 = rho_over_M * D_AB / D_ch                    # beta / Sh, a pure property group
beta = b0 * Sh0
print(f"reconstructed:  dT_AD = {dTad:.2f} K      (from the printed dH = "
      f"{pf('delta_H')} MJ/kmol, C_G0 = {CG0}, and M_G C_PG = {M_G*C_PG:.0f} J/kmol.K)")
print(f"                beta  = {beta:.5e} kmol/m2.s per unit mole fraction")
print(f"                       = (P/RgT) Sh D_AB / D  with D_AB = {D_AB:.4e} m2/s at "
      f"{T_props:.0f} K, D = {D_ch*1e3:.2f} mm\n")
from scipy.optimize import minimize_scalar
cw_at_max = lambda d: float(minimize_scalar(
    lambda c: -r_of(c, pf("TG_inlet_fig1"), d), bracket=(1e-5, 0.003, 0.04),
    method="brent", options=dict(xtol=1e-13)).x)
cmax = cw_at_max(dTad)
tgt = pf("Cw_at_rate_maximum")
print(f"THE ONE PRINTED FACT THAT TOUCHES THIS AT ALL: the paper states 'the rate")
print(f"        reaches a maximum at C_w = {tgt}' for T_G0 = {pf('TG_inlet_fig1'):.0f} C. Root-finding the")
print(f"        maximum of eq. (5) on the eq. (14) locus gives C_w = {cmax:.5f}, which")
print(f"        rounds to the one significant figure the paper prints.\n")
print(f"        BUT MEASURE ITS POWER BEFORE CALLING IT A CHECK. C_w at the maximum")
print(f"        barely depends on the reconstructed dT_AD:")
for f_ in (0.0, 0.5, 1.0, 2.0):
    print(f"          dT_AD = {f_*dTad:8.1f} K -> C_w,max = {cw_at_max(f_*dTad):.5f}")
half_ulp = 0.5 * 10.0 ** np.floor(np.log10(tgt))     # the paper prints ONE figure
d_hi = brentq(lambda d: cw_at_max(d) - (tgt + half_ulp), 0.5 * dTad, 8.0 * dTad, xtol=1e-6)
d_lo_ = minimize_scalar(cw_at_max, bracket=(0.1 * dTad, 0.5 * dTad, 1.0 * dTad),
                        method="brent", options=dict(xtol=1e-10))
ka_ = KA0 * np.exp(-EAR / (pf("TG_inlet_fig1") + 273.0))
print(f"""
        Root-found rather than sampled: C_w,max stays inside [{tgt-half_ulp}, {tgt+half_ulp}], i.e. still
        rounds to the printed {tgt}, for EVERY dT_AD from 0 up to {d_hi:.0f} K - its
        minimum over that whole range is {d_lo_.fun:.5f}, at dT_AD = {d_lo_.x:.0f} K, which
        clears the lower edge by {d_lo_.fun-(tgt-half_ulp):.1e}. The
        reconstructed {dTad:.0f} K sits in the middle of a range spanning a factor of {d_hi/dTad:.1f}.
        SO THIS DOES NOT CHECK dT_AD. What it does check is the transcription of the
        kinetics: at dT_AD = 0 the maximum is exactly 1/k_a(T_G0) = {1/ka_:.5f}, and
        every value above is a mild distortion of that by the adiabatic locus. k_a0
        and E_a/R were PRINTED, so this is a transcription check, not a check on a
        reconstructed group - and beta is checked by nothing here at all, since
        C_w,max does not contain it.

        The other number an earlier draft called a check: the same dT_AD puts a fully
        lit wall at T_G0 + dT_AD = {pf('TG_inlet_fig1')+dTad:.0f} C for a {pf('TG_inlet_fig1'):.0f} C inlet. The paper prints no
        adiabatic flame temperature in words anywhere in pp. 477-483, so there is
        nothing to compare it with; it is a CONSEQUENCE of the reconstruction and is
        reported as one.""")
M["dTad_reconstructed"] = dTad
M["beta_reconstructed"] = beta
M["rate_max_Cw"] = float(cmax)
M["rate_max_Cw_printed_absdiff"] = abs(float(cmax) - tgt)'''))

cells.append(md(r"""#### What "light-off threshold" has to mean

The paper quotes two windows, and each one describes a different figure:

* 1-D, p. 479, about **Fig. 2**: "For inlet gas temperatures below **304 °C**, the
  monolith temperature does not rise significantly … between 304° and **343 °C**,
  the reaction lights-off at some point in the channel."
* 2-D, p. 481, about **Fig. 11**: "no reaction light-off for low inlet gas
  temperature (case 1), $T^{o}_G <$ **293 °C**; light-off in the tube for
  intermediate inlet gas temperatures (case 2 and 3),
  293 $< T^{o}_G <$ **343 °C**."

**Neither is a root-found threshold.** 304, 316, 320 and 343 are the four inlet
temperatures of Fig. 2; 293, 316, 343 and 349 are the four cases of Fig. 11. Each
window is quoted from whichever sampled case happened to fall on the right side.

**And the two are quoted for channels of different length**, which is the whole
reason they look inconsistent with each other. "Lights off at some point in the
channel" is a statement about a channel of some stated length, and the two figures
do not share one: Fig. 2's abscissa is labelled 0.00 to **0.20** in
`X/(D RE PR)`, Fig. 11's 0 to **0.50** — tick labels, transcribed as printed
characters like every other number on this page, with no curve read. A longer
channel lights off at a lower inlet temperature, so a window read off a channel
2.5 times as long sits lower for that reason alone. The next cell runs each of the
paper's two statements at the length of the figure it is made about.

**A retraction, stated rather than quietly fixed.** An earlier draft of this page
read the two windows as thresholds on one channel, found the 2-D one 11 °C *below*
the 1-D one, and called that a contradiction with the mechanism the authors
themselves derive on p. 482 ("the one-dimensional model slightly underestimates
the Nusselt number; lower Nusselt number means earlier light-off"). **That claim
was wrong and it is withdrawn here and in every metadata file that carried it.**
Run at their own lengths, this page's own models reproduce *both* of the paper's
statements and reproduce the apparent 11 °C gap along with them — so the paper's
two numbers are consistent with each other and with its own mechanism, and the
gap chiefly measures the two figures' abscissae rather than the two models. *Chiefly*
and not *only*: the cell below decomposes it into a length term and a model term, and
the model term is several degrees in the direction p. 482 requires — it is what fixes
the *sign* that is length, not the whole of the gap. What survives
is the part that was always the contribution: the windows are **sampled**, and on
a **common** length the 2-D model sits *above* the 1-D, which is the direction
p. 482 requires.

Both thresholds are root-findable, and both are defined here without reference to
any figure:

* **upper** — light-off *at the inlet*: the inlet temperature at which the cold
  branch has already folded at $C_G = C_G^{o}$. Length-independent.
* **lower** — light-off leaves the channel: the inlet temperature at which the
  light-off position $G^{*}$ equals the channel length $L^{*}$. This one needs a
  length, which the paper never states; $L^{*}$ is swept."""))

cells.append(code(r'''Gmax_ = 0.40; dG_ = 1e-4
Gmax_long = 0.90          # long enough to contain Fig. 11's own channel and beyond
NuF_const = lambda Sh: (lambda G: Sh)
NuF_GT = lambda G: NuQ_GT(max(1000.0 * G, 1e-9))     # Grigull-Tratz, constant wall flux

def G_1d_quad(T, Sh=Sh0, **kw):
    g = lightoff_1d_quad(T, dTad, b0 * Sh, Sh, **kw)[0]
    return g if np.isfinite(g) else 10.0

def G_1d_march(T, NuFun=None, **kw):
    g = march_1d(T, dTad, b0, NuFun or NuF_const(Sh0), dG_, Gmax_, **kw)[0]
    return g if np.isfinite(g) else 10.0

def G_2d(T, Gmax=None, **kw):
    g = march_2d(T, dTad, beta, Sh0, dG_, Gmax or Gmax_, **kw)[0]
    return g if np.isfinite(g) else 10.0'''))

cells.append(code(r'''t0 = time.time()
# EACH OF THE PAPER'S TWO STATEMENTS, RUN AT THE LENGTH OF THE FIGURE IT IS ABOUT.
L_f2, L_f11 = pf("fig2_abscissa_end"), pf("fig11_abscissa_end")
T_f2, T_f11 = pf("lightoff_window_1D_low"), pf("lightoff_window_2D_low")
g1_low = G_1d_quad(T_f2)                       # the 1-D case the paper says does not light off
g2_low = G_2d(T_f11, Gmax=Gmax_long)           # the 2-D case it says does not light off
print(f"Fig. 2  (1-D) runs to L* = {L_f2:.2f}. The paper: below {T_f2:.0f} C 'the monolith")
print(f"        temperature does not rise significantly'. This model puts that case's")
print(f"        light-off at G* = {g1_low:.4f}, i.e. {'BEYOND' if g1_low > L_f2 else 'INSIDE'} the figure -> the statement holds.")
print(f"Fig. 11 (2-D) runs to L* = {L_f11:.2f}. The paper: 'no reaction light-off' below")
print(f"        {T_f11:.0f} C. This model puts that case's light-off at G* = {g2_low:.4f},")
print(f"        i.e. {'BEYOND' if g2_low > L_f11 else 'INSIDE'} the figure -> that statement holds too.\n")
t1_f2 = brentq(lambda T: G_1d_quad(T) - L_f2, 280.0, 380.0, xtol=1e-3)
t2_f11 = brentq(lambda T: G_2d(T, Gmax=Gmax_long) - L_f11, 280.0, 380.0, xtol=1e-3)
gap_model = t2_f11 - t1_f2
gap_paper = T_f11 - T_f2
# The gap is NOT length alone: it is length PLUS the difference between the two
# models, and the second term is this page's own result. Decomposed rather than
# asserted, by carrying the 1-D model out to Fig. 11's length as well.
t1_f11 = brentq(lambda T: G_1d_quad(T) - L_f11, 280.0, 380.0, xtol=1e-3)
d_len, d_mod = t1_f11 - t1_f2, t2_f11 - t1_f11
print(f"Root-found at each figure's own length:")
print(f"  1-D, Nu = {Sh0}, L* = {L_f2:.2f} (Fig. 2)   {t1_f2:7.2f} C")
print(f"  2-D, round tube, L* = {L_f11:.2f} (Fig. 11)  {t2_f11:7.2f} C")
M["G_lightoff_1d_at_window_low"] = g1_low
M["G_lightoff_2d_at_window_low"] = g2_low
M["T_low_1d_at_fig2_length"] = t1_f2
M["T_low_2d_at_fig11_length"] = t2_f11
M["T_windows_length_artefact_C"] = gap_model
print(f"""
so this page's own two models, on identical parameters, put the 2-D threshold {gap_model:+.2f} C
relative to the 1-D one WHEN EACH IS READ AT ITS OWN FIGURE'S LENGTH - against the
{gap_paper:+.0f} C between the paper's two quoted windows, which it reproduces to {abs(gap_model-gap_paper):.2f} C.

and that gap is ALMOST, BUT NOT ENTIRELY, the length. Taking the SAME 1-D model out to
Fig. 11's own length separates the two effects:

  length alone, 1-D from L* = {L_f2:.2f} to {L_f11:.2f} {d_len:+8.2f} C
  model difference at the common L* = {L_f11:.2f} {d_mod:+8.2f} C   <- this page's own result
  ------------------------------------------------------
  the gap between the two quoted windows   {d_len+d_mod:+8.2f} C

So the length term dominates and it is what FIXES THE SIGN; the two models still
differ by several degrees along the way, in the direction p. 482 requires, and that
term is not deleted by the retraction.

THE APPARENT INVERSION IS A CHANNEL-LENGTH ARTEFACT. Both of the paper's statements
are correct about their own figures, they are consistent with each other once the
lengths are accounted for, and they are consistent with the mechanism p. 482 derives.
An earlier draft of this page called them a contradiction; that is withdrawn. The
comparison that DOES test the two models is the next one, where the length is held
common - and there the 2-D sits above the 1-D, as p. 482 requires.
{time.time()-t0:.1f} s""")'''))

cells.append(code(r'''t0 = time.time()
# upper thresholds: length-independent, root-found on the fold condition
Tup = {}
Tup["1-D, Nu = 3.608 (square)"] = brentq(lambda T: fold_CG(T, dTad, b0 * Sh0) - CG0,
                                         280.0, 420.0, xtol=1e-8)
Tup["1-D, Nu = 4.364 (circle)"] = brentq(lambda T: fold_CG(T, dTad, b0 * 4.364) - CG0,
                                         280.0, 420.0, xtol=1e-8)
for k, v in Tup.items():
    print(f"upper threshold, light-off at the inlet:  {k:34s} {v:7.2f} C")
g349 = G_2d(pf("lightoff_window_2D_inlet_case"))
print(f"""
The paper's 1-D bracket for this, taken from its running text alone, is (304, 343] C:
the reaction "lights-off at some point in the channel" between 304 and 343, and above
343 "the monolith temperature is equal to the adiabatic flame temperature over its
entire length". Nothing here is read off a curve.

The 2-D model has NO such threshold, and this is the one place on the page where a
discretisation is load-bearing. A round tube has an infinite local Sherwood number at
X* = 0, so light-off always occurs at some X* > 0 and a STRICT inlet threshold does
not exist. What a march reports instead is the inlet temperature at which light-off
falls inside the FIRST AXIAL STEP - a number that depends on the step and on nothing
else. Refined rather than asserted:""")
lad_dG = []
for dGx in (4e-4, 2e-4, dG_, 5e-5):
    lad_dG.append(brentq(lambda T: (lambda g: g if np.isfinite(g) else 10.0)(
        march_2d(T, dTad, beta, Sh0, dGx, Gmax_)[0]) - 1e-9, 280.0, 470.0, xtol=1e-3))
    print(f"  dG = {dGx:.1e}   'light-off in the first step' at {lad_dG[-1]:7.2f} C"
          + ("   <- the dG used everywhere else on this page" if dGx == dG_ else ""))
div = lad_dG[-1] - lad_dG[-2]
print(f"""
  +{np.mean(np.diff(lad_dG)):.1f} C per halving and still rising, with no sign of a limit: the quantity
  DIVERGES as dG -> 0, exactly as the explanation above requires. It is reported as
  what it is - a resolution-dependent number at dG = {dG_:.0e} - and named for its step.
  Its break table carries a dG row for that reason; a row that moved it by a fraction
  of a percent would make a divergent number look converged.

What the paper means by "light-off at the tube inlet for high inlet gas temperatures
(case 4)" is therefore visual, not a threshold. Case 4 is {pf('lightoff_window_2D_inlet_case'):.0f} C, and this model puts
its light-off at G* = {g349:.5f} - nonzero, but a fiftieth of the abscissa of the figure it
is drawn on. The comparable, unambiguous quantity is the pair of 1-D values above,
which differ ONLY in the Nusselt number: {Tup["1-D, Nu = 4.364 (circle)"]-Tup["1-D, Nu = 3.608 (square)"]:.2f} C between the square-channel value
the paper's 1-D runs use and the round-tube value its 2-D model has. That is a SHAPE
effect and not dimensionality - and because it is measured on a different threshold
from the +4 C below, the like-for-like split is redone on the lower threshold two
cells down rather than carried across.""")
M["T_up_1d_square"] = Tup["1-D, Nu = 3.608 (square)"]
M["T_up_1d_circle"] = Tup["1-D, Nu = 4.364 (circle)"]
M["T_up_1d_shape_effect"] = Tup["1-D, Nu = 4.364 (circle)"] - Tup["1-D, Nu = 3.608 (square)"]
M["T_up_2d_first_step_dG1em4"] = lad_dG[2]
M["T_up_2d_first_step_divergence_C"] = div
M["G_lightoff_2d_at_case4"] = g349
print(f"{time.time()-t0:.1f} s")'''))

cells.append(code(r'''t0 = time.time()
rows = []
# L* is swept because the paper never prints a channel length. The two ends the paper
# DOES fix are read from the CSV rather than typed - Fig. 2's abscissa and Figs. 3-5's
# - with 0.30 as a margin past both, and every conclusion is checked over the sweep.
for L in (L_f2, pf("fig3_abscissa_end"), 0.30):
    a = brentq(lambda T: G_1d_quad(T) - L, 280.0, 380.0, xtol=1e-3)
    b = brentq(lambda T: G_1d_march(T, NuF_GT) - L, 280.0, 380.0, xtol=1e-3)
    c = brentq(lambda T: G_2d(T) - L, 280.0, 380.0, xtol=1e-3)
    rows.append((L, a, b, c, c - a, c - b))
low = pd.DataFrame(rows, columns=["L*", "1-D const Nu=3.608", "1-D Grigull-Tratz Nu(X*)",
                                  "2-D", "2-D minus 1-D const", "2-D minus 1-D G-T"])
show(low.style.set_uuid("i13_lowthr").format({"L*": "{:.2f}", "1-D const Nu=3.608": "{:.2f}",
                          "1-D Grigull-Tratz Nu(X*)": "{:.2f}", "2-D": "{:.2f}",
                          "2-D minus 1-D const": "{:+.2f}", "2-D minus 1-D G-T": "{:+.2f}"}))
# the row the *_L025 metrics are read from, located rather than indexed by hand
i25 = int(np.flatnonzero(np.isclose(low["L*"].to_numpy(), 0.25))[0])
M["T_low_1d_const_L025"] = low.iloc[i25, 1]
M["T_low_1d_GT_L025"] = low.iloc[i25, 2]
M["T_low_2d_L025"] = low.iloc[i25, 3]
M["T_low_2d_minus_1d_const_L025"] = low.iloc[i25, 4]
print(f"""
The 2-D model needs a HIGHER inlet temperature than either 1-D model, at every
channel length: +{low.iloc[i25,4]:.2f} C against the constant-Nu 1-D and +{low.iloc[i25,5]:.2f} C against the
paper's own Grigull-Tratz 1-D, at L* = {low.iloc[i25,0]:.2f}. That is the direction the paper
derives on p. 482. The sign does not depend on L*: it is the same at {low.iloc[0,0]:.2f} and at {low.iloc[-1,0]:.2f}.
The MAGNITUDE does depend on it - the threshold falls about {(low.iloc[0,1]-low.iloc[-1,1])/((low.iloc[-1,0]-low.iloc[0,0])/0.05):.1f} C per 0.05 of L* -
which is why the paper's two windows, read off figures whose abscissae differ by a
factor of {L_f11/L_f2:.1f}, cannot be compared with each other as thresholds.
The 2-D minus 1-D column shrinks with length too - +{low.iloc[0,4]:.2f} C here at L* = {low.iloc[0,0]:.2f}, and
+{d_mod:.2f} C at Fig. 11's L* = {L_f11:.2f} in the decomposition above - so the model difference the
retraction leaves standing is +{min(low.iloc[0,4], d_mod):.1f} to +{max(low.iloc[0,4], d_mod):.1f} C over the whole range of lengths at issue.

WHY THIS RANGE OF L*: the paper never prints a channel length, so the sweep is set by
the only lengths it does fix, its own one-dimensional figures' abscissae - {L_f2:.2f} for
Fig. 2 and {pf("fig3_abscissa_end"):.2f} for Figs. 3-5 - with {low.iloc[-1,0]:.2f} as a margin past both. Nothing above
depends on the choice.

NONE OF THESE ABSOLUTE TEMPERATURES IS THE PAPER'S NUMBER, and the {low.iloc[i25,1]:.1f} C in the
first column is not evidence for the {pf('lightoff_window_1D_low'):.0f} C it lands within a degree of. Two of the
model's groups are reconstructed (section 7f moves this column by tens of degrees),
and the paper's 304 is a sampled case rather than a threshold in the first place.
Only the DIFFERENCES between columns, which share those groups, are claimed.
{time.time()-t0:.1f} s""")'''))

cells.append(md(r"""#### Shape or dimensionality — split like-for-like on one threshold

The 1-D runs of the paper use $Nu = 3.608$, the *square*-channel value from its own
Table 1, while its 2-D model is a *round* tube. Two different things are therefore
bundled into "the 1-D/2-D difference", and they separate cleanly by running the 1-D
model at the round-tube Nusselt number: everything below the middle row is channel
shape, everything above it is the thermal entrance region plus dimensionality
proper. Done on **the same threshold and the same length** the +4 °C above was
measured on — an earlier draft attached a shape effect measured on the *upper*
threshold to a gap measured on the *lower* one, and the two are not interchangeable."""))

cells.append(code(r'''t0 = time.time()
L_ = float(low.iloc[i25, 0])
t_sq = float(low.iloc[i25, 1])                                   # 1-D, Nu = 3.608
t_ci = brentq(lambda T: G_1d_quad(T, Sh=4.364) - L_, 280.0, 380.0, xtol=1e-3)
t_2d = float(low.iloc[i25, 3])                                   # 2-D, round tube
print(f"at L* = {L_:.2f}, all three on identical parameters:")
print(f"  1-D, Nu = {Sh0} (square - what the paper's 1-D uses)  {t_sq:7.2f} C")
print(f"  1-D, Nu = 4.364 (round tube - shape matched to the 2-D) {t_ci:7.2f} C")
print(f"  2-D  (round tube, entrance region resolved)             {t_2d:7.2f} C")
M["T_low_1d_circle_L025"] = t_ci
M["T_low_shape_effect_L025"] = t_ci - t_sq
M["T_low_entrance_dim_L025"] = t_2d - t_ci
print(f"""
  channel SHAPE alone              {t_ci-t_sq:+.2f} C
  entrance region + dimensionality {t_2d-t_ci:+.2f} C
  ------------------------------------------
  the whole gap                    {t_2d-t_sq:+.2f} C

so the split is almost exactly half and half: {100*(t_ci-t_sq)/(t_2d-t_sq):.0f} % of what gets called the
"1-D/2-D difference" here is the Nusselt number of a SQUARE channel being compared
with a ROUND one, and has nothing to do with the number of dimensions. The paper's
own 1-D and 2-D runs carry that same mismatch, since its Figs. 2-6 use 3.608 and its
Fig. 11 is a cylinder.
{time.time()-t0:.1f} s""")'''))

cells.append(code(r'''t0 = time.time()
# the paper's own measure of adequacy: an equivalent shift in inlet temperature
rows = []
for T2 in (310.0, 316.0, 320.0, 330.0, 340.0):
    g2 = G_2d(T2)
    if g2 > 1.0:
        continue
    t1c = brentq(lambda T: G_1d_quad(T) - g2, 280.0, 420.0, xtol=1e-3)
    t1v = brentq(lambda T: G_1d_march(T, NuF_GT) - g2, 280.0, 420.0, xtol=1e-3)
    rows.append((T2, g2, T2 - t1c, T2 - t1v))
eq = pd.DataFrame(rows, columns=["2-D inlet T (C)", "2-D light-off G*",
                                 "equivalent shift, const-Nu 1-D (C)",
                                 "equivalent shift, Grigull-Tratz 1-D (C)"])
show(eq.style.set_uuid("i13_equiv").format({"2-D inlet T (C)": "{:.0f}", "2-D light-off G*": "{:.5f}",
                         "equivalent shift, const-Nu 1-D (C)": "{:+.2f}",
                         "equivalent shift, Grigull-Tratz 1-D (C)": "{:+.2f}"}))
j = int(np.argmin(abs(eq["2-D inlet T (C)"] - 316.0)))
M["equiv_shift_const_Nu_at316"] = float(eq.iloc[j, 2])
M["equiv_shift_GT_at316"] = float(eq.iloc[j, 3])
M["equiv_shift_GT_max"] = float(eq.iloc[:, 3].max())
M["equiv_shift_const_max"] = float(eq.iloc[:, 2].max())
print(f"""
The paper claims the 1-D/2-D difference is "equivalent to less than a 2 C
difference in inlet gas temperature", and it makes that claim for the
Grigull-Tratz 1-D model of its Fig. 12, not for the constant-Nu model of its
Figs. 2-5. Both are measured here:

  * Grigull-Tratz 1-D : {eq.iloc[j,3]:+.2f} C at a 316 C inlet, up to {eq.iloc[:,3].max():+.2f} C over the range.
  * constant-Nu 1-D   : {eq.iloc[j,2]:+.2f} C at a 316 C inlet, up to {eq.iloc[:,2].max():+.2f} C over the range.

So the claim is the right size for the model it was made about and roughly
{eq.iloc[j,2]/eq.iloc[j,3]:.1f}x too small for the simpler model that produced most of the paper's
figures. "Adequate" is a property of the Fig. 12 procedure, not of
one-dimensionality.
{time.time()-t0:.1f} s""")'''))

cells.append(md(r"""#### Where the difference actually comes from

The three sources are separable, and two of them are not about dimensionality at
all:

1. **Channel shape.** The paper's 1-D runs use $Nu = 3.608$, the *square*-channel
   value from its own Table 1; the 2-D model is a *round* tube, whose asymptote is
   4.364. The cell above measures it on the same threshold and the same length as
   the gap itself, and it comes out at **about half** of that gap.
2. **The thermal entrance region.** The 2-D model starts with an infinite local
   Nusselt number and relaxes towards the asymptote. A constant-$Nu$ 1-D model
   never has that, which is exactly why the paper's Fig. 12 procedure — the
   Grigull–Tratz correlation before light-off, restarted after it — closes most of
   the gap.
3. **Dimensionality proper**, i.e. what is left once 1 and 2 are matched.

The measurement above says the residue is small: with the entrance region put
back in, the two models differ by a couple of degrees of inlet temperature, which
is the paper's conclusion and it survives. What does not survive is the way the
conclusion is stated — as a property of *one-dimensional models*, when the model
that carries the paper's own parametric study (Figs. 2–6, constant $Nu$) is
several times further off than the claim allows."""))

# ------------------------------------------------------------------ validation
cells.append(md(r"""## Validation

### 7a. Which $X^{*}$ the paper's 2-D results use — decided by its own printed number

The Notation defines $X^{*} = 4x/(D\,Re\,Pr)$; eq. (15) requires
$X^{*} = x/(D\,Re\,Pr)$. The paper prints one number that separates them: for the
constant-wall-temperature Graetz problem with $\bar T_G^{*} = 1$ at the inlet and
$T_w^{*} = 0$, "at the point $X^{*} = 0.25$, $\bar T_G^{*}$ is 0.0224"."""))

cells.append(code(r'''t0 = time.time()
targ = pf("Tbar_at_X_0p25")
out = {}
for lab, half, X in (("eq. (15) as printed", 0.5, 0.30), ("Notation, 4x/(D Re Pr)", 2.0, 0.40)):
    g = Graetz(X=X, n_x=1600, n_r=200, half=half, x_grade=3.0)
    s = g.solve(("temp", 0.0), inlet=1.0)
    out[lab] = float(np.interp(0.25, s["x"], s["cup"]))
    print(f"{lab:26s}: Tbar_cup(X*=0.25) = {out[lab]:.6f}   "
          f"({out[lab]/targ:.2f} x the printed {targ})")
# the area mean, in case the paper's overbar meant that (its Notation offers both)
print(f"{'area mean, eq. (15)':26s}: {float(np.interp(0.25, s['x'], s['area'])):.6f}"
      if False else "", end="")
g = Graetz(X=0.30, n_x=1600, n_r=200, half=0.5, x_grade=3.0)
s15 = g.solve(("temp", 0.0), inlet=1.0)
area25 = float(np.interp(0.25, s15["x"], s15["area"]))
print(f"{'area mean, eq. (15)':26s}: {area25:.6f}   ({area25/targ:.2f} x printed)")
print(f"""
The Notation's scaling misses by a factor of {out['Notation, 4x/(D Re Pr)']/targ:.1f}; eq. (15)'s scaling misses by
{100*abs(out['eq. (15) as printed']/targ-1):.1f} %, and the AREA mean by {100*abs(area25/targ-1):.0f} %. So the 2-D results use eq. (15)'s
scaling and the overbar is the mixing-cup mean - both settled by the paper's own
printed number, neither assumed.""")
M["Tbar_X025_eq15"] = out["eq. (15) as printed"]
M["Tbar_X025_notation"] = out["Notation, 4x/(D Re Pr)"]
M["Tbar_X025_area"] = area25
M["Tbar_X025_notation_over_printed"] = out["Notation, 4x/(D Re Pr)"] / targ
print(f"{time.time()-t0:.1f} s")'''))

cells.append(code(r'''t0 = time.time()
# ... and the residual 6 % is the paper's own grid, not a transcription error.
print("radial grid  Tbar_cup(X* = 0.25)     (the paper states it used TEN radial points)")
lad = []
for n_r, grade in ((10, 1.0), (10, 2.0), (20, 1.0), (40, 1.0), (80, 2.0), (200, 2.0)):
    g = Graetz(X=0.30, n_x=1600, n_r=n_r, half=0.5, x_grade=3.0, r_grade=grade)
    sd_ = g.solve(("temp", 0.0), inlet=1.0)
    v = float(np.interp(0.25, sd_["x"], sd_["cup"]))
    lad.append((n_r, grade, v))
    print(f"  n_r = {n_r:3d}, grading {grade:.1f}   {v:.6f}")
conv = lad[-1][2]
n10 = [v for (n, gr, v) in lad if n == 10]
print(f"""
Converged: {conv:.6f}. At the paper's own ten radial points the answer lies between
{min(n10):.6f} (uniform) and {max(n10):.6f} (wall-clustered), and the printed 0.0224 sits
inside that interval - the wall-clustered ten-point solve is {100*(max(n10)/conv-1):.1f} % above converged
and the printed value is {100*(targ/conv-1):.1f} % above it. So this printed number is not a check on
the paper's 2-D model, it is a MEASUREMENT of that model's radial discretisation
error, and {100*(targ/conv-1):.1f} % is the size of it. That is the same
order as the 1-D/2-D difference the model is being used to adjudicate, and it is
why the adequacy question above is posed as an equivalent inlet-temperature shift
between two models run on the SAME grid rather than as an absolute accuracy.""")
M["Tbar_X025_converged"] = conv
M["Tbar_X025_paper_grid_error"] = targ / conv - 1.0
print(f"{time.time()-t0:.1f} s")'''))

cells.append(md(r"""### 7b–c. Refinement: every axis that carries error, and both axes of every 2-D solve"""))

cells.append(code(r'''t0 = time.time()
# --- square duct: the two Cartesian axes, refined separately ------------------
ref = duct(384, 384)
print("square duct, axes refined SEPARATELY against a (384,384) reference")
for axis in (0, 1):
    errs_H, errs_T = [], []
    for n in (24, 48, 96, 192):
        d = duct(n, 384) if axis == 0 else duct(384, n)
        errs_H.append(abs(d["Nu_H1"] - ref["Nu_H1"])); errs_T.append(abs(d["Nu_T"] - ref["Nu_T"]))
    oH, oT = observed_order(errs_H)[:2], observed_order(errs_T)[:2]
    print(f"  axis {axis}: Nu_H1 orders {oH[0]:.2f}, {oH[1]:.2f}   "
          f"Nu_T orders {oT[0]:.2f}, {oT[1]:.2f}")
    M[f"order_square_NuH1_axis{axis}"] = oH[1]
    M[f"order_square_NuT_axis{axis}"] = oT[1]
print("""
STRUCTURAL CAVEAT, stated rather than hidden: the square is symmetric under
swapping the axes, so these two columns are not independent measurements - they
must agree, and they do to three digits. The genuinely independent pair of axes
on this page is the Graetz solve below, where axis 0 carries convection and axis 1
carries diffusion.""")
print(f"{time.time()-t0:.1f} s")'''))

cells.append(code(r'''t0 = time.time()
# --- Graetz flux solve: axial (convection) and radial (diffusion) -------------
def pole_of(n_x, n_r, tvd=True, half=0.5, r_grade=1.0):
    g = Graetz(X=0.05, n_x=n_x, n_r=n_r, half=half, r_grade=r_grade)
    s = g.solve(("flux", Q9), tvd=tvd)
    return cross(s["x"], s["Tw"] - s["cup"])

ref_pole = pole_of(6400, 320)
print(f"reference pole (n_x=6400, n_r=320): X* = {ref_pole:.7f}")
for nm, ladder in (("radial", [(3200, n) for n in (10, 20, 40, 80)]),
                   ("axial ", [(n, 40) for n in (50, 100, 200, 400)])):
    base = ref_pole if nm == "radial" else pole_of(6400, 40)
    errs = [abs(pole_of(*g) - base) for g in ladder]
    ns = [g[1] if nm == "radial" else g[0] for g in ladder]
    o = observed_order(errs); fo = fitted_order(ns, errs)
    print(f"  {nm}: errors " + "  ".join(f"{e:.2e}" for e in errs) +
          "   rung orders " + ", ".join(f"{v:.2f}" for v in o) +
          f"   least-squares slope {fo:.2f}")
    M[f"order_pole_{nm.strip()}"] = fo
    if nm.strip() == "radial":
        radial_errs = errs                    # kept so the next rung's comparison is
    if nm.strip() == "axial":                 # computed, not retyped
        print(f"          The axial ladder is NOISY and its slope comes out above 2. "
              f"That is not\n          third-order accuracy: the pole is LOCATED on this "
              f"axis, so the locator's\n          own interpolation error is inside the "
              f"ladder, and the axial errors are\n          already {errs[0]/radial_errs[0]:.2f}x the "
              f"radial ones at the coarsest rung, so this is not\n          the limiting "
              f"axis. What the ladder DOES show cleanly is its break row:\n          with the "
              f"TVD deferred correction off the same fit gives first order,\n          which is "
              f"the upwind scheme it falls back to.")
M["nu_pole_X_fig9_refined"] = ref_pole
p10 = pole_10r                                # the same ten-point solve section 2 quotes
# THE PAPER'S OWN DISCRETISATION, not just its radial count: p. 480 states ten radial
# points AND dX* = 0.00025, and a 1955-vintage marching code is first order in X*, so
# the TVD correction goes off too. This is the strongest form of the grid question -
# it asks what the paper's OWN scheme would have produced.
s_pg = Graetz(X=0.05, n_x=int(round(0.05 / pf("axial_step_2D"))),
              n_r=int(pf("radial_points_2D")), r_grade=1.0).solve(("flux", Q9), tvd=False)
p_pg = cross(s_pg["x"], s_pg["Tw"] - s_pg["cup"])
M["nu_pole_X_fig9_paper_scheme"] = p_pg
print(f"""
Each axis is measured against a reference converged on the OTHER axis, so neither
order is read on the other axis's floor. The pole is remarkably insensitive to the
radial grid: at the paper's own TEN radial points it is X* = {p10:.6f} against
{ref_pole:.6f} converged, a difference of {100*abs(p10/ref_pole-1):.2f} %.

The stronger form of that question is what the paper's OWN scheme would have given,
not just its radial count. Its stated discretisation is ten radial points AND
dX* = {pf("axial_step_2D")}, and a 1955-vintage marching code is first order in X*, so the TVD
correction is switched off as well. That solve puts the pole at X* = {p_pg:.6f} -
{100*(p_pg/ref_pole-1):+.1f} % from converged, and {100*(1-p_pg/0.032):.1f} % SHORT of 0.032.

HOW SHORT, EXACTLY, IS SCHEME-DEPENDENT: {p_pg:.6f} is this page's own radial operator
run on the paper's stated counts, not the paper's 1955 code, which is not recoverable.
So the conclusion is stated as a MARGIN rather than as that number's last digits. The
two paper-count solves here span {min(p_pg,p10):.6f} to {max(p_pg,p10):.6f}, i.e. {100*(1-max(p_pg,p10)/0.032):.1f} to {100*(1-min(p_pg,p10)/0.032):.1f} % short,
and both sit within {max(abs(p_pg/ref_pole-1), abs(p10/ref_pole-1))*100:.1f} % of converged - while REACHING 0.032 takes {100*(0.032/ref_pole-1):+.1f} %, i.e.
{(0.032/ref_pole-1)/max(abs(p_pg/ref_pole-1), abs(p10/ref_pole-1)):.0f}x the discretisation error the paper's own stated grid produces here.
No discretisation tried here comes near its own (0.032, 0.033) bracket: the gap is
not a grid effect in either direction.
{time.time()-t0:.1f} s""")'''))

cells.append(md(r"""### 7d. The Grigull–Tratz correlations, against the solve — and a third $X^{*}$

The paper reprints two correlations from Grigull and Tratz (1965) and uses the
first of them as the Nusselt number of its Fig. 12 one-dimensional model. They
are printed without a definition of the $X^{*}$ they use, and it is neither of the
two already in play."""))

cells.append(code(r'''t0 = time.time()
g = Graetz(X=0.06, n_x=4000, n_r=200, r_grade=2.0, x_grade=3.0)
s = g.solve(("flux", lambda X: -np.ones_like(X)))          # constant wall flux
Nu_num = local_nusselt(s)   # q < 0 and (Tw - Tbar) < 0 here, so this is positive
print("  G = x/(D Re Pr)   solved Nu    G-T with X*=1000G   with X*=G    with X*=4G")
devs = []
for gq in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
    v = float(np.interp(gq, s["x"], Nu_num))
    devs.append(abs(NuQ_GT(1000 * gq) / v - 1.0))
    print(f"  {gq:14.4f} {v:12.4f} {NuQ_GT(1000*gq):16.4f} {NuQ_GT(gq):13.1f} {NuQ_GT(4*gq):12.1f}")
M["grigull_tratz_max_reldev"] = float(max(devs))
print(f"""
The correlation reproduces this page's own solve to {100*max(devs):.2f} % over that whole range
with X* = 1000 x/(D Re Pr) - which is how Grigull and Tratz wrote it - and is
wrong by one to two orders of magnitude under either of the paper's own two
meanings for X*. So the paper carries THREE incompatible definitions of one
symbol: 4x/(D Re Pr) in the Notation and eq. (9), x/(D Re Pr) in eqs. (15)-(16)
and on every figure abscissa, and 1000x/(D Re Pr) in the two reprinted
correlations. Reported, not repaired: the page uses each equation with the
scaling that equation requires, and says so at each use.
{time.time()-t0:.1f} s""")'''))

cells.append(md(r"""### 7e. Two independent routes to the light-off threshold

The constant-$Nu$ 1-D light-off position is computed two ways that share no
discretisation: a marched finite-volume integration of $dC_G/dG$, and a
quadrature of $dG/dC_G$ with the wall concentration root-found at each
$C_G$. The fold itself is root-found from $r'(C_w) = -\beta$ in both."""))

cells.append(code(r'''t0 = time.time()
print("  T_G0    march (dG=1e-4)   quadrature in C_w   diff/dG      quad in C_G")
ds, dq, n_giveup = [], [], 0
for T in (305.0, 310.0, 316.0, 322.0, 328.0):
    a = G_1d_march(T); b = G_1d_quad(T)
    with warnings.catch_warnings(record=True) as wlog:
        # scipy's subdivision-limit warning is EVIDENCE here, not noise - but its text
        # carries a kernel-local temp path, so it is captured and counted rather than
        # printed, which keeps two runs of this notebook byte-identical.
        warnings.simplefilter("always")
        c = lightoff_1d_quad_CG(T, dTad, beta, Sh0)
    n_giveup += any("subdivisions" in str(w.message) for w in wlog)
    ds.append(abs(a - b) / dG_); dq.append(abs(c - b) / b)
    print(f"  {T:6.1f} {a:16.6f} {b:18.6f} {ds[-1]:13.2f} {c:14.6f}")
print(f"  (scipy's adaptive quad exhausted its 600 subdivisions on {n_giveup} of the "
      f"{len(ds)} rows\n   of the last column, and still returned the wrong answer)")
M["lightoff_march_vs_quad_max_steps"] = float(max(ds))
M["quad_CG_endpoint_error"] = float(max(dq))
print(f"""
worst march-vs-quadrature difference {max(ds):.2f} axial steps - i.e. the march reports
the step at which the fold is crossed, and agrees with the quadrature to within one
step everywhere. Expressed as a relative error it looks worse at high inlet
temperature only because G* itself is small there; the dG row of the break table
confirms it scales with dG and not with anything else.

The last column is a TRAP worth recording. It is the same integral parametrised by
C_G instead of C_w, evaluated with scipy's adaptive `quad`, and it is up to
{100*max(dq):.2f} % out. dC_G/dC_w vanishes at the fold, so that form has a square-root
endpoint; the C_w form's integrand vanishes there instead and is converged at 100
Gauss points. Had the C_G form been the reference, a {100*max(dq):.2f} % error would have sat
under every threshold on this page and no break row would have shown it - which is
why the two routes here are two PARAMETRISATIONS as well as two discretisations.""")
print(f"{time.time()-t0:.1f} s")'''))

cells.append(md(r"""### 7f. How much of §3 survives the reconstruction

Two groups were reconstructed. The *thresholds themselves* depend on them
strongly; the *difference between the two models*, which is what the page
concludes about, does not. Measured rather than asserted."""))

cells.append(code(r'''t0 = time.time()
base = dict(dTad=dTad, beta=beta, b0=b0)
rows = []
for lab, f_dT, f_b in (("baseline", 1.00, 1.00), ("dT_AD +10 %", 1.10, 1.00),
                       ("dT_AD -10 %", 0.90, 1.00), ("beta x 1.5", 1.00, 1.50),
                       ("beta / 1.5", 1.00, 1 / 1.5)):
    dTad, b0 = base["dTad"] * f_dT, base["b0"] * f_b
    beta = b0 * Sh0
    tu1 = brentq(lambda T: fold_CG(T, dTad, beta) - CG0, 250.0, 500.0, xtol=1e-6)
    g2 = G_2d(316.0)
    t1c = brentq(lambda T: G_1d_quad(T) - g2, 250.0, 500.0, xtol=1e-3) if g2 < 1 else np.nan
    t1v = brentq(lambda T: G_1d_march(T, NuF_GT) - g2, 250.0, 500.0, xtol=1e-3) if g2 < 1 else np.nan
    rows.append((lab, tu1, g2 if g2 < 1 else np.nan, 316.0 - t1c, 316.0 - t1v))
dTad, b0, beta = base["dTad"], base["b0"], base["beta"]
sens = pd.DataFrame(rows, columns=["reconstruction", "1-D upper threshold (C)",
                                   "2-D light-off G* at 316 C",
                                   "shift, const-Nu 1-D (C)", "shift, G-T 1-D (C)"])
show(sens.style.set_uuid("i13_sens").format({"1-D upper threshold (C)": "{:.2f}",
                           "2-D light-off G* at 316 C": "{:.5f}",
                           "shift, const-Nu 1-D (C)": "{:+.2f}",
                           "shift, G-T 1-D (C)": "{:+.2f}"}))
sp_thr = float(sens["1-D upper threshold (C)"].max() - sens["1-D upper threshold (C)"].min())
sp_eq = float(sens["shift, const-Nu 1-D (C)"].max() - sens["shift, const-Nu 1-D (C)"].min())
lo_c, lo_v = float(sens["shift, const-Nu 1-D (C)"].min()), float(sens["shift, G-T 1-D (C)"].min())
M["recon_spread_upper_threshold_C"] = sp_thr
M["recon_spread_equiv_shift_C"] = sp_eq
M["recon_min_equiv_shift_const"] = lo_c
M["recon_min_equiv_shift_GT"] = lo_v
print(f"""
Read this honestly. A +/-10 % move in dT_AD and a factor 1.5 either way in beta
swing the ABSOLUTE upper threshold by {sp_thr:.1f} C, and the equivalent shift by {sp_eq:.2f} C -
so the shift is about {sp_thr/sp_eq:.0f}x less sensitive in degrees, but it is NOT insensitive,
and no number on this page should be quoted to better than that.

What survives every one of these variations is the two statements the page actually
makes: the constant-Nu 1-D model exceeds the paper's 2 C claim under EVERY variation
(smallest value {lo_c:+.2f} C), and the Grigull-Tratz 1-D model is at or above it under
every variation (smallest {lo_v:+.2f} C). Neither conclusion turns on the reconstruction;
the individual temperatures do, and they are labelled accordingly.

What no perturbation of this kind can test is whether the model STRUCTURE is the
paper's: constant oxygen, Le = 1 and Sh = Nu are the paper's own assumptions, and
the first of them is a break-table row rather than a certainty.
{time.time()-t0:.1f} s""")'''))

cells.append(md(r"""### 7g. The boundary read

Both marched models take their wall value from `compute_boundary_values`, not
from the last cell centre. For a Neumann wall the two differ at first order in the
near-wall spacing, and the wall value is precisely what the kinetics are evaluated
at, so the difference propagates straight into the rate."""))

cells.append(code(r'''t0 = time.time()
rad = Radial(48)
c0 = np.exp(-3.0 * (1.0 - rad.r_c))          # an arbitrary smooth profile
for q in (0.0, -0.7):
    bv = rad.wall_value(c0, q)
    naive = c0[-1]
    print(f"  q = {q:+.1f}:  compute_boundary_values -> {bv:.8f}   "
          f"last cell centre -> {naive:.8f}   difference {abs(bv-naive):.2e}")
# the Dirichlet case must return the boundary condition itself, exactly
gg = Graetz(X=0.05, n_x=200, n_r=40)
sd = gg.solve(("temp", 0.37), inlet=1.0)
print(f"  Dirichlet wall at 0.37: compute_boundary_values returns "
      f"{sd['Tw'][0]:.12f} (max deviation over the whole tube "
      f"{np.max(abs(sd['Tw']-0.37)):.2e})")
M["boundary_read_dirichlet_exact"] = float(np.max(abs(sd["Tw"] - 0.37)))
M["boundary_read_vs_lastcell"] = float(abs(rad.wall_value(c0, -0.7) - c0[-1]))
print(f"{time.time()-t0:.1f} s")'''))

cells.append(md(r"""### 8. Defect injection — and what it cannot reach

Every metric reported below gets a row that moves it. The rows are **built from
the solver arguments**, so a row cannot silently perturb nothing: the cell asserts
that each perturbation actually changes the number it names, and prints the
assertion result rather than trusting it."""))

cells.append(code(r'''t0 = time.time()
brk = []
COV = {}                       # metric -> list of (row, relative move)


def row(name, metrics, note):
    """Run a perturbation and record every metric it moves. `metrics` maps a
    metric name to its perturbed value."""
    for k, v in metrics.items():
        base = M[k]
        moved = abs(v - base) / max(abs(base), 1e-30) if np.isfinite(v) else np.inf
        brk.append((name, k, base, v, moved, note))
        COV.setdefault(k, []).append((name, moved))


# --- Table 1 solves ----------------------------------------------------------
d_ = duct(128, geom="square", plug=True)
row("plug velocity instead of Poiseuille",
    {"nu_h1_square_pymrm": d_["Nu_H1"], "nu_t_square_pymrm": d_["Nu_T"],
     "nu_h1_square_printed_relerr": abs(float(T1.Nu_H[1]) - d_["Nu_H1"]) / d_["Nu_H1"],
     "nu_t_square_printed_relerr": abs(float(T1.Nu_T[1]) - d_["Nu_T"]) / d_["Nu_T"],
     "nu_h1_square_two_routes_reldiff": abs(sp_lad[-1][0] - d_["Nu_H1"]) / d_["Nu_H1"]},
    "both duct problems are driven by u")
d_ = duct(64, geom="square", wall_neumann=True)
row("Neumann wall instead of Dirichlet",
    {"nu_h1_square_pymrm": d_["Nu_H1"]}, "removes the wall condition entirely")
d_ = duct(128, geom="square", weight_eigen=False)
row("unweighted eigenproblem (drop u)",
    {"nu_t_square_pymrm": d_["Nu_T"],
     "nu_t_square_spectral_N24": d_["Nu_T"]},
    "gives the plain Laplacian eigenvalue instead of the Graetz one")
d_ = duct(256, geom="circle", nu_r=0)
row("nu=0 (Cartesian) in the radial div",
    {"nu_t_circle_pymrm": d_["Nu_T"], "nu_h1_circle_pymrm": d_["Nu_H1"],
     "nu_t_circle_printed_relerr": abs(float(T1.Nu_T[0]) - d_["Nu_T"]) / d_["Nu_T"],
     "nu_h1_circle_solve_vs_exact": abs(d_["Nu_H1"] - float(Nu_circle_exact)) / float(Nu_circle_exact)},
    "the geometry factor of construct_div; it also breaks the agreement with 48/11")

# --- Graetz solves -----------------------------------------------------------
s10 = Graetz(X=0.05, n_x=6000, n_r=10, r_grade=1.0).solve(("flux", Q9))
p10_ = cross(s10["x"], s10["Tw"] - s10["cup"])
s10b = Graetz(X=0.05, n_x=6000, n_r=10, r_grade=1.0).solve(("flux", Q8))
row("the paper's own 10 radial points",
    {"nu_pole_X_fig9": p10_, "nu_pole_X_fig9_refined": p10_,
     "nu_pole_X_fig8": cross(s10b["x"], s10b["Tw"] - s10b["cup"]),
     "nu_pole_over_fluxzero_fig9": p10_ / M["fluxzero_X_fig9"]},
    "10 radial points is what p. 480 states")
s2_ = Graetz(X=0.20, n_x=2500, n_r=60, half=2.0, r_grade=1.0).solve(("flux", Q9))
row("half = 2 (the Notation's X*)",
    {"nu_pole_X_fig9": cross(s2_["x"], s2_["Tw"] - s2_["cup"])},
    "the factor-of-four question, made a solver argument")
def _pole_nt(nx):
    ss = Graetz(X=0.05, n_x=nx, n_r=40, r_grade=1.0).solve(("flux", Q9), tvd=False)
    return cross(ss["x"], ss["Tw"] - ss["cup"])
_b = pole_of(6400, 40)
_ns = (50, 100, 200, 400)
row("TVD deferred correction off",
    {"order_pole_axial": fitted_order(_ns, [abs(_pole_nt(n) - _b) for n in _ns])},
    "first-order upwind instead of van Leer, on the SAME ladder and the same fit")
sd10 = Graetz(X=0.30, n_x=1600, n_r=10, half=0.5, x_grade=3.0,
              r_grade=1.0).solve(("temp", 0.0), inlet=1.0)
v10 = float(np.interp(0.25, sd10["x"], sd10["cup"]))
row("the paper's own 10 radial points",
    {"Tbar_X025_eq15": v10, "Tbar_X025_converged": v10,
     "Tbar_X025_paper_grid_error": targ / v10 - 1.0,
     "Tbar_X025_area": float(np.interp(0.25, sd10["x"], sd10["area"]))},
    "this is the row that explains the printed 0.0224")
row("X* = G in Grigull-Tratz",
    {"grigull_tratz_max_reldev":
     float(max(abs(NuQ_GT(gq) / float(np.interp(gq, s["x"], Nu_num)) - 1.0)
               for gq in (0.001, 0.01, 0.05)))}, "the correlations' own X* scaling")

# --- the reacting model ------------------------------------------------------
def react_metrics(**kw):
    """Recompute every reacting metric under a perturbation of the kinetics."""
    tu_s = brentq(lambda T: fold_CG(T, dTad, b0 * Sh0, **kw) - CG0, 150.0, 900.0, xtol=1e-5)
    tu_c = brentq(lambda T: fold_CG(T, dTad, b0 * 4.364, **kw) - CG0, 150.0, 900.0, xtol=1e-5)
    g2 = (lambda g: g if g < 1 else 10.0)(march_2d(316.0, dTad, beta, Sh0, dG_, Gmax_, **kw)[0])
    cw = minimize_scalar(lambda c: -r_of(c, pf("TG_inlet_fig1"), dTad, **kw),
                         bracket=(1e-5, 0.003, 0.04), method="brent",
                         options=dict(xtol=1e-13)).x
    out = {"T_up_1d_square": tu_s, "T_up_1d_circle": tu_c,
           "T_up_1d_shape_effect": tu_c - tu_s, "rate_max_Cw": float(cw),
           "rate_max_Cw_printed_absdiff": abs(float(cw) - pf("Cw_at_rate_maximum"))}
    if 1e-4 < g2 < 1:      # skipped when the perturbed model lights off AT the inlet,
        t1c = brentq(lambda T: (lambda g: g if g < 1 else 10.0)(   # where no shift exists
            lightoff_1d_quad(T, dTad, b0 * Sh0, Sh0, **kw)[0]) - g2, 150.0, 900.0, xtol=1e-3)
        out["equiv_shift_const_Nu_at316"] = 316.0 - t1c
    return out

row("E_a sign flipped (+961)", react_metrics(ear=-EAR),
    "the negative activation energy IS the CO inhibition")
row("k_a0 = 0 (no CO inhibition term)", react_metrics(ka0=0.0),
    "removes the [1+k_a C]^2 denominator; the fold survives because T_w still "
    "rises as C_w falls, but it moves by over 100 C")
row("oxygen depleted stoichiometrically", react_metrics(o2="depleted"),
    "the paper writes no oxygen balance; this is the alternative reading")
_slab = lambda T, Gm: (lambda g: g if np.isfinite(g) else 10.0)(
    march_2d(T, dTad, beta, Sh0, dG_, Gm, nu_r=0)[0])
row("nu=0 (Cartesian) in the 2-D radial div",
    {"T_low_2d_L025": brentq(lambda T: _slab(T, Gmax_) - 0.25, 280.0, 380.0, xtol=1e-3),
     "T_up_2d_first_step_dG1em4": brentq(lambda T: _slab(T, Gmax_) - 1e-9,
                                         280.0, 520.0, xtol=1e-3),
     "G_lightoff_2d_at_case4": _slab(pf("lightoff_window_2D_inlet_case"), Gmax_),
     "T_low_2d_at_fig11_length": brentq(lambda T: _slab(T, Gmax_long) - L_f11,
                                        280.0, 380.0, xtol=1e-3),
     "G_lightoff_2d_at_window_low": _slab(T_f11, Gmax_long)},
    "cylindrical -> slab weighting of the radial diffusion")
row("dG x 2 (2e-4) in the 2-D march",
    {"T_up_2d_first_step_dG1em4": brentq(
        lambda T: (lambda g: g if np.isfinite(g) else 10.0)(
            march_2d(T, dTad, beta, Sh0, 2e-4, Gmax_)[0]) - 1e-9,
        280.0, 470.0, xtol=1e-3)},
    "THE POINT of this row: that metric has no continuum limit, so ONE step-size\n"
    "     change moves it by degrees. Every other 2-D number on the page is\n"
    "     step-independent to a hundredth of a degree.")
_g316 = G_2d(316.0)
row("Grigull-Tratz replaced by a constant Nu",
    {"T_low_1d_GT_L025": brentq(lambda T: G_1d_march(T, NuF_const(Sh0)) - 0.25,
                                280.0, 380.0, xtol=1e-3),
     "equiv_shift_GT_at316": 316.0 - brentq(
         lambda T: G_1d_march(T, NuF_const(4.364)) - _g316, 280.0, 400.0, xtol=1e-3)},
    "drops the thermal entrance region the correlation supplies: Nu -> 3.608 for\n     the threshold and -> 4.364 for the shift")
_m1 = lambda T, Nu: march_1d(T, dTad, b0, NuF_const(Nu), 1e-3, Gmax_)[0]
row("dG x 10 (1e-3), and the marched route in place of the quadrature",
    {"lightoff_march_vs_quad_max_steps":
     float(max(abs(_m1(T, Sh0) - G_1d_quad(T)) / 1e-3 for T in (310.0, 316.0, 322.0))),
     "T_low_1d_const_L025": brentq(lambda T: _m1(T, Sh0) - 0.25, 280.0, 380.0, xtol=1e-3),
     "T_low_1d_circle_L025": brentq(lambda T: _m1(T, 4.364) - 0.25, 280.0, 380.0, xtol=1e-3),
     "T_low_1d_at_fig2_length": brentq(lambda T: _m1(T, Sh0) - L_f2, 280.0, 380.0, xtol=1e-3),
     "G_lightoff_1d_at_window_low": _m1(T_f2, Sh0)},
    "the axial step of the marched route, and a route swap for the four thresholds\n"
    "     that the baseline gets from the C_w quadrature instead")
with warnings.catch_warnings():
    # scipy's "maximum number of subdivisions reached" is the POINT of this row, and
    # its message carries a kernel-local temp path, which would make two runs of this
    # notebook differ in their stream output. Silenced deliberately, not hidden.
    warnings.simplefilter("ignore")
    _q30 = float(max(abs(lightoff_1d_quad_CG(T, dTad, beta, Sh0, limit=30) / G_1d_quad(T) - 1.0)
                     for T in (310.0, 316.0, 322.0)))
row("quad subdivision limit 600 -> 30", {"quad_CG_endpoint_error": _q30},
    "the adaptive integrator's own effort on the square-root endpoint; it warns that\n"
    "     it gave up, which is the row working")
row("last cell instead of the boundary read",
    {"boundary_read_vs_lastcell": 0.0},
    "the metric IS that difference, so this row is the identity - it is a real row\n"
    "     (it drives the number to zero) and the metric is counted as moved, not\n"
    "     declared structural as well")

bt = pd.DataFrame(brk, columns=["injected defect", "metric", "baseline", "broken",
                                "rel. move", "what it perturbs"])
show(bt.style.set_uuid("i13_breaks").format({"baseline": "{:.6g}", "broken": "{:.6g}", "rel. move": "{:.2e}"}))
print(f"{len(brk)} (row, metric) pairs over {len(set(r[0] for r in brk))} injected defects; "
      f"{len(COV)} of {len(M)} metrics moved by at least one row.")
print(f"{time.time()-t0:.1f} s")'''))

cells.append(code(r'''# COVERAGE MAP, asserted key-for-key against what report_agreement will write.
STRUCTURAL = {
 "nu_h1_circle_exact": "48/11, a symbolic identity: no solver argument produces it, "
   "so none can move it. Live companion: nu_h1_circle_solve_vs_exact, which the "
   "nu=0 row moves by a factor of 2.",
 "nu_h1_triangle_exact": "28/9, the same. It has no finite-volume companion on this "
   "page - the triangle is solved only in closed form - so its live companion is "
   "nu_h1_triangle_printed_relerr, which is fixed for the same reason and is "
   "declared below.",
 "nu_h1_circle_printed_relerr": "exact value against a printed digit: both sides are "
   "fixed numbers, so nothing can move it. It is a transcription check, not a solve.",
 "nu_h1_triangle_printed_relerr": "the same.",
 "fluxzero_X_fig9": "the root of a polynomial the paper prints, so it is fixed. "
   "Live companion: nu_pole_over_fluxzero_fig9, which the 10-point row moves.",
 "nu_zero_X_fig9": "equals fluxzero_X_fig9 by construction - the Nusselt numerator IS "
   "the wall flux - and the next metric measures how well the solve reproduces that.",
 "nu_zero_minus_fluxzero_fig9": "BELOW ABS_FLOOR (1e-12) is not expected here, but it "
   "is an agreement between a root-found zero and a printed polynomial's root, so it "
   "is a solver accuracy check rather than a physical result.",
 "fig9_sign_0p032_matches_paper": "a boolean, and it is 0, so it sits AT the CI floor "
   "and is outside check_agreement's comparison. It is kept because it is the "
   "page's finding in its bluntest form; the above-floor carrier of the same "
   "information is nu_pole_X_fig9, which has two rows.",
 "fig9_sign_0p033_matches_paper": "a boolean, and it is 1. Same companion.",
 "dTad_reconstructed": "an INPUT, not a result. Section 7f is its perturbation study "
   "and moves it by +/-10 % explicitly.",
 "beta_reconstructed": "an INPUT. Section 7f moves it by a factor 1.5 either way.",
 "nu_pole_X_fig9_paper_scheme": "this metric IS a perturbation's output - the pole on "
   "the paper's own stated grid AND its own first-order scheme, i.e. a coarse-grid, "
   "TVD-off variant - so it is a break row rather than a baseline. Live companion: "
   "nu_pole_X_fig9, which two rows move.",
 "Tbar_X025_notation": "this metric IS the half=2 perturbation's own output - it is "
   "what the Notation's scaling gives - so it is a break row rather than a baseline.",
 "Tbar_X025_notation_over_printed": "the same ratio, restated.",
 "order_square_NuH1_axis0": "a measured convergence rate. `duct` has no first-order "
   "variant, so no injected defect on this page can change it; the quantities it "
   "qualifies (nu_h1_square_pymrm, nu_t_square_pymrm) each have two rows. The live "
   "order metric on this page is order_pole_axial, which the TVD row collapses.",
 "order_square_NuT_axis0": "the same.",
 "order_square_NuH1_axis1": "the same, AND structurally identical to axis 0: the "
   "square is symmetric under swapping the axes, so this is not an independent "
   "measurement and the page says so where it is reported.",
 "order_square_NuT_axis1": "the same.",
 "order_pole_radial": "a measured convergence rate with no first-order radial variant "
   "available. Companion: order_pole_axial, which the TVD row collapses from 2 to 1.",
 "recon_spread_upper_threshold_C": "this IS a perturbation result - the spread over "
   "section 7f's own variations - so it has no separate break row by construction.",
 "recon_spread_equiv_shift_C": "the same.",
 "recon_min_equiv_shift_const": "the same.",
 "recon_min_equiv_shift_GT": "the same.",
 "equiv_shift_GT_max": "the maximum over the inlet-temperature scan; its at-316 "
   "companion equiv_shift_GT_at316 carries a row.",
 "equiv_shift_const_max": "the same, with equiv_shift_const_Nu_at316 as the companion.",
 "boundary_read_dirichlet_exact": "EXACTLY ZERO, so it is below ABS_FLOOR = 1e-12 and "
   "outside the regression suite entirely. It is an identity - pymrm's boundary "
   "reconstruction must return a Dirichlet value unchanged - and it is kept as one. "
   "Above-floor companion: boundary_read_vs_lastcell.",
 "T_low_2d_minus_1d_const_L025": "a difference of two metrics that each carry a row.",
 "T_windows_length_artefact_C": "the same - the difference of the two figure-length "
   "thresholds, each of which carries a row.",
 "T_low_shape_effect_L025": "a difference of two metrics that each carry a row.",
 "T_low_entrance_dim_L025": "the same.",
 "T_up_2d_first_step_divergence_C": "this IS a refinement result - what one halving "
   "of dG does to a quantity that has no continuum limit - so it has no separate "
   "break row by construction, exactly like the section 7f spreads.",
}
missing = sorted(set(M) - set(COV) - set(STRUCTURAL))
# DEAD IS PER (ROW, METRIC) PAIR, not per metric. An earlier draft took the MAXIMUM
# over a metric's rows, which silenced a row that moved nothing whenever some other
# row moved the same metric - i.e. it was weaker than the sentence above it claims.
dead = sorted({(r, k) for k, v in COV.items() for r, m in v
               if m < 1e-6 and k not in STRUCTURAL})
overlap = sorted(set(COV) & set(STRUCTURAL))
extra = sorted(set(STRUCTURAL) - set(M))
# ... and the same hole on the OTHER side: a break row naming something that is not a
# metric at all would otherwise be counted as coverage and assert nothing.
stray = sorted(set(COV) - set(M))
pairs = [m for v in COV.values() for _, m in v]
print(f"metrics in agreement.json : {len(M)}")
print(f"  moved by >= 1 break row : {len(set(COV) & set(M))}")
print(f"  declared structural     : {len(set(STRUCTURAL) & set(M))}")
print(f"  the two must PARTITION the metric set: "
      f"{len(set(COV) & set(M))} + {len(set(STRUCTURAL) & set(M))} = "
      f"{len(set(COV) & set(M)) + len(set(STRUCTURAL) & set(M))} of {len(M)}")
print(f"  NEITHER (must be empty) : {missing}")
print(f"  BOTH   (must be empty)  : {overlap}")
print(f"  (row, metric) pairs that moved nothing : {dead}")
print(f"  structural keys that are not metrics : {extra}")
print(f"  perturbed keys that are not metrics  : {stray}")
print(f"  weakest single (row, metric) pair, over all {len(pairs)}: "
      f"{min(pairs):.2e} relative")
assert not missing, missing
assert not overlap, overlap
assert not dead, dead
assert not extra, extra
assert not stray, stray
assert len(set(COV) & set(M)) + len(set(STRUCTURAL) & set(M)) == len(M)
below_floor = sorted(k for k, v in M.items() if abs(v) < 1e-12)
print(f"\n  metrics below ABS_FLOOR = 1e-12, outside check_agreement entirely: {below_floor}")
print("  each is named in STRUCTURAL above with an above-floor companion.")'''))

cells.append(md(r"""**What the table cannot do, stated rather than left implicit.**

* `nu_h1_circle_exact` and `nu_h1_triangle_exact` are **symbolic identities**
  (48/11 and 28/9). No solver argument can move them, because no solver produces
  them. Their break rows perturb the *finite-volume solve that agrees with them*,
  which is the only thing a perturbation can reach; the companion metrics
  `nu_h1_circle_printed_relerr` and `nu_h1_triangle_printed_relerr` carry the
  comparison against the printed digits and are above the CI floor.
* `boundary_read_vs_lastcell` is a **difference between two reads of the same
  field**, so its break row is the identity — the row replaces the boundary
  reconstruction by the last cell centre and drives the metric to exactly zero.
  It is counted as *moved*, not as structural: an earlier draft listed it in both,
  which double-counted it in the coverage total. The same applied to the spectral
  square-duct value.
* No row perturbs the *model structure*. Constant oxygen, $Le = 1$, $Sh \cong Nu$
  and "the wall temperature follows the adiabatic locus" are the paper's
  assumptions, adopted here. The oxygen row is the one place where an assumption
  rather than a number is varied, and it moves the threshold by several degrees —
  which is the honest size of that particular unknown.
* **The reconstruction is not a defect and cannot be injected as one.** §7f
  measures its effect directly instead, and that measurement is why every
  conclusion in §3 is phrased as a difference between two models rather than as an
  absolute."""))

# ---------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Table 1 stops being a citation.** Five of its eight cells are derived here, two
of them in closed form — $48/11$ for the circular duct and $28/9$ for the
equilateral triangle, both obtained from exact polynomial solutions of the
cross-section problems rather than quoted from anywhere. The square pair is
computed on two discretisations that share no code. That matters because Table 1
is the baseline for the paper's first claim: you cannot call a Nusselt number
unusual without the usual one.

**The Nusselt singularity is root-found, and its cause is corrected.** The paper
attributes the infinite Nusselt numbers to the wall flux passing through zero. At
that point eq. (19) has a vanishing *numerator*, so the Nusselt number passes
through zero; the pole is the zero of $T_w - \bar T_G$ and sits **8.7 % further
downstream**, in both of the flux functions the paper prints. The paper's own next
paragraph gives the correct mechanism, so this is a wrong summary rather than a
missed effect — but the summary is the sentence a reader carries away.

**The 1-D/2-D thresholds are root-found instead of sampled, and on a common
length.** The paper's two light-off windows are quoted from the four inlet
temperatures each figure happened to run — 304, 316, 320, 343 in Fig. 2, and 293,
316, 343, 349 in Fig. 11 — so neither is a threshold. They are also quoted at
different channel lengths (Fig. 2's abscissa is labelled to 0.20, Fig. 11's to
0.50), and a longer channel lights off at a lower inlet temperature: run at each
figure's own length, this page's own models reproduce *both* windows and the
apparent 11 °C between them. **An earlier draft of this page called that apparent
gap a contradiction with the mechanism the authors derive on p. 482. It is not one,
and the claim is withdrawn** — it was a channel-length artefact, and the page now
measures it as such. Root-found on a *common* length and with identical parameters,
the 2-D model sits above the 1-D at every length tested, which is the direction
p. 482 requires.

**And about half of that "1-D/2-D difference" is not dimensionality at all.**
Split like-for-like on one threshold and one length: running the 1-D model at the
round-tube $Nu = 4.364$ instead of the square-channel 3.608 its own Table 1
supplies accounts for about half the gap, the entrance region and dimensionality
together for the other half. The paper's 1-D figures use the square value and its
2-D model is a cylinder, so its own comparison carries that mismatch.

**"Adequate" is given a number and a scope.** The paper's "less than a 2 °C
difference in inlet gas temperature" is measured here as an equivalent inlet-shift
between the two models. It holds, at about that size, for the variable-$Nu$
one-dimensional model of Fig. 12 — and it is several times larger for the
constant-$Nu$ model that produced Figs. 2–6, i.e. for most of the paper's
parametric conclusions. The adequacy belongs to the Fig. 12 *procedure*, not to
one-dimensionality.

**Three incompatible meanings of $X^{*}$ are separated**, each settled against a
computation rather than by inspection, and one of them is settled by the paper's
own printed check value.

**And that check value turns out to measure the paper's grid.** "At the point
$X^{*} = 0.25$, $\bar T_G^{*}$ is 0.0224" sits about 6 % above the converged
answer and is reproduced to a fraction of a percent by a solve on the ten radial
points p. 480 says were used. The paper's 2-D model therefore carries a radial
discretisation error of the same order as the 1-D/2-D difference it is being used
to settle — which is exactly why this page's adequacy statement is a difference
between two models on one grid, and not an absolute accuracy."""))

cells.append(md(r"""### What this page cannot conclude

* **Nothing about any figure.** Digitisation is out of scope here, so Figs. 1–13
  are used only for *characters printed inside their frames* — annotation boxes,
  abscissa tick labels, one curve label and one case label — each row labelled in
  the CSV with the figure it sits in, and **all of them enumerated and counted from
  the CSV** in *The data* above, which is what this promise is audited against.
  No profile shape, no light-off *position* against the paper's own, and none of
  the conversion curves of Figs. 4–6 and 12 is checked. The light-off positions computed here are
  therefore compared only against statements the paper makes in words.
* **The channel length is not printed anywhere.** The lower threshold needs one,
  so it is swept over 0.20–0.30 in Graetz units and every conclusion is checked to
  be independent of the choice. No single "the" lower threshold is claimed. The
  only lengths the paper does fix are its *figures'* abscissae — 0.20 for Fig. 2,
  0.25 for Figs. 3–5, 0.50 for Fig. 11 — and those are used **only** to set the
  swept range (the first two, which are the one-dimensional figures) and to read
  its two quoted windows at the lengths they were quoted for, never as *the*
  monolith's length.
* **There is no 2-D "light-off at the inlet" threshold to compute.** The local
  Sherwood number is infinite at $X^{*} = 0$, so light-off always happens at some
  $X^{*}>0$; what a march reports is the inlet temperature at which it falls inside
  the first step, and §3 shows that number diverging as the step shrinks. It is
  reported with its step size in its name and is not comparable with the paper's
  1-D bracket.
* **Two of the model's dimensional groups are reconstructions** from ordinary gas
  properties, not from the paper. §7f shows the absolute thresholds move by tens
  of degrees under plausible variations, so **no absolute threshold on this page
  should be read as the paper's number**; only differences between models run on
  the same parameters are.
* **The square-duct $Nu_T$ discrepancy is not attributed.** 2.976 is printed,
  2.9775 is computed, and Shah and London (1971) — the source Table 1 cites — was
  not consulted, so the page cannot say whether the last digit was lost in 1971 or
  in 1976.
* **Three of Table 1's eight cells are not reproduced**, and one of them cannot
  be: the sinusoidal geometry is printed as a glyph with no aspect ratio, so the
  paper does not determine the problem.
* **The transient model (eqs. 1–4, Fig. 7) is not built.** Its only quantitative
  statement — "the transience for this case is almost complete by 170 s" — needs
  the solid heat capacity and the monolith dimensions, none of which is printed.
  The neighbouring page `I1.2` (Oh & Cavendish 1982) does transient monolith
  light-off against a printed table."""))

cells.append(md(r"""## Reuse

`Graetz` is a general developing-flow solver for a round tube: pick the wall
condition with `("flux", f)` or `("temp", value)`, and the axial coefficient with
`half`. The wall flux may be **any function of the axial coordinate**, which is
what makes the paper's Figs. 8–9 reproducible; pymrm takes an array-valued `d` in
the boundary dictionary, so a varying Neumann value costs nothing and the operator
is still assembled once.

Copy the van Leer deferred correction **with its under-relaxation** (`omega=0.6`);
this is `A3.15`'s solver and the same caution applies — at `omega=1` the limiter
switches and the iteration limit-cycles.

`march_2d` is worth copying for any reacting-wall problem. The trick is that the
interior is *linear* in the wall flux, so each implicit step needs exactly two
linear solves (flux 0 and flux 1) and then a **scalar** root-find for the wall
value. That turns a 2-D nonlinear boundary problem into the same one-dimensional
algebra as Fig. 1's mass-transfer line, makes the fold structure explicit, and
lets the light-off point be root-found rather than stumbled over. It is also what
keeps the runtime in the tens of seconds.

Read wall values with `compute_boundary_values`, never from the last cell centre:
§7g measures the difference on this problem, and the wall value is the argument of
the kinetics.

## Cite the source, not this page

Heck, R. H., Wei, J. & Katzer, J. R., *Mathematical Modeling of Monolithic
Catalysts*, AIChE Journal **22**(3) 477–484 (1976), doi:10.1002/aic.690220310.
Received October 7, 1975; revision received January 16 and accepted January 17,
1976.

Table 1 is attributed by the paper to Shah and London (1971), Stanford University
Dept. of Mechanical Engineering Tech. Rept. No. 75 — **not consulted here**. The
kinetics are Voltz, Morgan, Liederman & Jacob, *Ind. Eng. Chem. Prod. Res.
Develop.* **12**, 295 (1973) — **not consulted here** either; the four constants
are used exactly as Heck et al. print them, which is the reprint route recorded in
`AGENTS.md`. The two Nusselt correlations are Grigull & Tratz, *Intern. J. Heat
Mass Transfer* **8**, 669 (1965) — not consulted; they are used as the paper
reprints them."""))

# ------------------------------------------------------------------- agreement
cells.append(code(r'''report_agreement("I1.3", M)'''))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.kernelspec = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata.language_info = {"name": "python", "version": "3"}
nbf.write(nb, "index.ipynb")
print(f"wrote index.ipynb with {len(cells)} cells")
