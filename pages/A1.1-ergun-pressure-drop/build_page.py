#!/usr/bin/env python3
"""Generate index.ipynb for page A1.1. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "Four packed-bed pressure-drop correlations, against the data Ergun fitted"
description: "Ergun, Kozeny-Carman, Darcy-Forchheimer and the Eisfeld-Schnitzlein wall correction on one dataset - 244 markers recovered from Ergun's own 1952 figure, which give back his 150 and 1.75."
categories: [sec:A, struct:S3, tier:T0, data:tier4, phase:gas-solid]
date: 2026-07-30
---

# Four packed-bed pressure-drop correlations, against the data Ergun fitted

**Catalog ID:** `A1.1` (also covers `A1.2` Kozeny–Carman, `A1.3`
Darcy–Forchheimer, `A1.4` Eisfeld–Schnitzlein) · **Structures:** `S3` ·
**Tier:** T0

The Ergun equation is the most-used correlation in reactor engineering and its
two constants, 150 and 1.75, are quoted far more often than they are questioned.
They were fitted, in 1952, to 640 experiments that survive only as scatter in
four figures of a scanned journal.

This page recovers 244 of those points and asks the constants to come back out."""))

cells.append(md(r"""## Background

Flow through a packed bed loses energy two ways at once. At low flow the loss is
viscous and linear in velocity; at high flow it is inertial and quadratic. Every
correlation in this family says the same thing — add the two — and differs only
in what multiplies each term.

**Kozeny–Carman** (`A1.2`) keeps the viscous term alone. It comes from the
hydraulic-radius picture: treat the void space as a bundle of tortuous channels,
apply Hagen–Poiseuille, and replace the tube diameter by the hydraulic radius of
the packing. It is exact in creeping flow and wrong by an order of magnitude in a
blower-driven bed.

**Darcy–Forchheimer** (`A1.3`) is the general two-term statement,
$-\mathrm{d}p/\mathrm{d}z = (\mu/K)\,U + \beta\rho U^2$, with a permeability $K$
and an inertial coefficient $\beta$ that are properties of the medium and are
measured, not predicted. It is the form used for foams, monoliths and
reservoir rock, where no particle diameter exists to build a correlation on.

**Ergun** (`A1.1`) is Darcy–Forchheimer with $K$ and $\beta$ predicted:
he fixes how each coefficient depends on voidage and particle size, leaving two
universal numbers to be fitted once, for everything.

**Eisfeld–Schnitzlein** (`A1.4`) is what happens when the tube is not much wider
than the particles. Near a wall the packing is looser and the fluid slips
through faster, while the wall itself adds friction; the two effects pull
opposite ways and their balance depends on Reynolds number. Below about
$D/d_p = 10$ this is a tens-of-percent correction, and Ergun knew it — he
discarded such systems from his own figure.

The catalog asks for these four on one dataset rather than four pages, and that
is what follows. The dataset is Ergun's."""))

cells.append(md(r"""## The published model

Ergun writes the pressure drop as (his eq. 12, journal page 91)

$$
\frac{\Delta P\, g_c}{L} \;=\; k_1\,\frac{(1-\epsilon)^2}{\epsilon^3}\,
\frac{\mu U_m}{D_p^2}
\;+\; k_2\,\frac{1-\epsilon}{\epsilon^3}\,\frac{G U_m}{D_p},
\qquad k_1 = 72\alpha,\quad k_2 = \tfrac34\beta ,
$$

with $D_p = 6/S_v$ the diameter of a sphere of the same specific surface,
$G = \rho U$ the mass flux, and $U_m$ the superficial velocity at the mean
pressure. Dividing by the viscous term gives the **linear form** that the
constants were actually fitted in (his eqs. 13, 13a, 13b):

$$
f_v \;\equiv\; \frac{\Delta P g_c}{L}\,\frac{D_p^2}{\mu U_m}\,
\frac{\epsilon^3}{(1-\epsilon)^2}
\;=\; k_1 + k_2\,\frac{N_{Re}}{1-\epsilon},
\qquad N_{Re} = \frac{D_p G}{\mu},
$$

and dividing by the *kinetic* term instead gives the Blake-type friction factor
(his eqs. 14, 14a, 14b):

$$
f_k \;\equiv\; \frac{\Delta P g_c}{L}\,\frac{D_p}{G U_m}\,
\frac{\epsilon^3}{1-\epsilon}
\;=\; k_1\,\frac{1-\epsilon}{N_{Re}} + k_2 ,
\qquad\text{so}\qquad f_k = f_v \Big/ \frac{N_{Re}}{1-\epsilon}.
$$

The two are the same numbers on two ordinates, which matters later. The result,
printed three times in the paper — in the abstract, after eq. (13b), and as
eq. (13c) — is

$$
\boxed{k_1 = 150,\qquad k_2 = 1.75,}
$$

"determined by the method of least squares … representing 640 experiments"
(journal page 91), over spheres, sand and pulverised coke with CO₂, N₂, CH₄ and
H₂.

**Kozeny–Carman** in the same coordinates is a horizontal line. Writing the
hydraulic-radius result as $\psi = 6^{3-n} k\,\mathrm{Re}^{n-2}
(1-\epsilon)^{3-n}/\epsilon^3$ with $n=1$ and Carman's sphere value $k=5$ gives
$f_v = 36k = 180$.

**Darcy–Forchheimer** is $f_v = A + Bx$ with $A$ and $B$ free. It is not a
competitor to Ergun on these axes — it is the *family*, and Ergun's equation is
one member of it. That is exactly what makes the refit below a fair test: fitting
Forchheimer to Ergun's own points is asking his data, without his arithmetic,
what $A$ and $B$ should be.

**Eisfeld–Schnitzlein** take Reichelt's wall-corrected Ergun form and refit it to
2391 points spanning $1.6 \le D/d_p \le 250$ (their eqs. 3 and 4):

$$
f_v = K_1 A_w^2 + \frac{A_w}{B_w}\,\frac{N_{Re}}{1-\epsilon},
\qquad
A_w = 1 + \frac{2}{3(D/d_p)(1-\epsilon)},
\qquad
B_w = \left[k_1\Big(\frac{d_p}{D}\Big)^2 + k_2\right]^2,
$$

with $K_1 = 154$, $k_1 = 1.15$, $k_2 = 0.87$ for spheres. Note there is no
separate turbulent constant: as $D/d_p \to \infty$ the kinetic coefficient tends
to $1/k_2^2 = 1.32$, not to 1.75."""))

cells.append(md(r"""## Parameters and assumptions

**Everything on this page is dimensionless**, which is what makes a 1952 paper,
a 2001 paper and a modern bed comparable at all. $f_v$, $f_k$ and
$N_{Re}/(1-\epsilon)$ are as defined above; $\epsilon$ is the bed voidage,
$D/d_p$ the tube-to-particle diameter ratio.

**Assumptions carried from the sources.** Single-phase, steady, isothermal,
incompressible *locally* (the compressibility of the gas is handled by
evaluating everything at the mean pressure, which is Ergun's own $U_m$);
particles characterised by one length $D_p = 6/S_v$; no size distribution
effects beyond that; a bed long enough that entrance and exit losses are
negligible.

**Deviation convention.** Everywhere below,
deviation $= (\text{model} - \text{measured})/\text{measured}$, for all four
correlations alike, so a positive number always means the correlation predicts
more pressure drop than was measured. At 5–20 % scatter a ratio and its
reciprocal are not interchangeable.

**One thing Ergun's figure does not carry: $D/d_p$.** His axes have no term for
it, so the wall effect cannot be located on his plot from his plot alone. It is
located instead from Eisfeld & Schnitzlein's Table 1, which re-catalogues the
same primary sources and gives each one's $D/d_p$ range. That mapping is
source-level, never point-level — see the sidecar."""))

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
import pandas as pd
import matplotlib.pyplot as plt
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A1.1-ergun-pressure-drop"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

R_GAS = 8.314462618          # J/(mol K)'''))

cells.append(md(r"""## The data

The Ergun paper is a **pure scan with no text layer at all** — `pdftotext`
returns four bytes for four pages — so every number here was read off a 600 dpi
render, and the dataset was digitised from the figures. There are no tables in
the paper.

**244 markers from Figure 7**, whose three panels are Ergun's own runs, Burke &
Plummer (1928) and Morcom (1946), all replotted by him onto one pair of axes.
Figure 7 rather than Figure 6 because Figure 6 pools everything into a single
band 640 markers deep, while Figure 7 separates the sources and is sparse enough
to resolve.

**How the markers were found matters for what follows.** They are open circles
about 13 px across on logarithmic graph paper, and the detector does not look
for circles: it labels the enclosed *white* regions and keeps the small ones,
i.e. the hole inside each ring. Ergun draws his fitted line straight through
those same points. A line has no enclosed interior, so the detector is
structurally incapable of locking onto it — which is what stops the round trip
below from being circular. The check is in the data itself and is repeated in
the Validation section.

**74 markers from Figure 5's top panel** are the *same* measurements on the
other ordinate. Ergun states it explicitly (journal page 93): "Figure 7 shows
$f_k$ plotted vs. $N_{Re}/(1-\epsilon)$ for the data already presented in
Figure 5." Two digitisations of one set of runs, on differently scaled axes,
with independently fitted calibrations — a free consistency check that the paper
pays for.

**The digitisation has been reviewed by a maintainer against numbered overlays,
and passed with one stated limitation.** No ring is on anything that is not a
marker, and none has slid onto the drawn line in the dense chain — but a few
markers per panel were missed. The next cell prints the recorded verdict, and
what that limitation does and does not affect is set out in the sidecar under
`limitations.recall` and again at the end of this page. Read it before reusing
the CSV: this file is a large sample of Ergun's Figure 7, not a census of it."""))

cells.append(code('''obs = load_data("ergun-1952-fig7-markers.csv", page=PAGE)
obs5 = load_data("ergun-1952-fig5-present-markers.csv", page=PAGE)
par = load_data("ergun-1952-parameters.csv", page=PAGE)
eis = load_data("eisfeld-2001-wall-correction.csv", page=PAGE)
src = load_data("eisfeld-2001-table1-sources.csv", page=PAGE)
obs_meta = load_meta("ergun-1952-fig7-markers.csv", page=PAGE)

P = dict(zip(par.quantity, par.value))
K1_ERGUN = P["k1_viscous_constant"]          # 150
K2_ERGUN = P["k2_kinetic_constant"]          # 1.75
KC_FV = P["carman_kozeny_fv"]                # 180 = 36 * 5
KC_DRAWN = P["kozeny_carman_line_fv_measured"]

E = {(q, s): v for q, s, v in zip(eis.quantity, eis.particle_shape, eis.value)}
EIS = dict(K1=E[("K1", "spheres")], k1=E[("k1", "spheres")], k2=E[("k2", "spheres")])

x = obs.Re_over_one_minus_eps.values
f_k = obs.f_k.values
f_v = obs.f_v.values

print(f"{len(obs)} markers from Ergun's Figure 7, in three panels:")
for s, g in obs.groupby("series"):
    print(f"   {s:22s} n={len(g):3d}   "
          f"Re/(1-eps) {g.Re_over_one_minus_eps.min():7.2f} - "
          f"{g.Re_over_one_minus_eps.max():7.1f}")
print(f"\\n{len(obs5)} markers from Figure 5's top panel (the same runs, ordinate f_v)")
print(f"\\nErgun's printed constants: k1 = {K1_ERGUN:g}, k2 = {K2_ERGUN:g}, "
      f"from {P['n_experiments_fitted']:.0f} experiments")
print(f"Eisfeld & Schnitzlein spheres: K1 = {EIS['K1']:g}, "
      f"k1 = {EIS['k1']:g}, k2 = {EIS['k2']:g}")
print(f"\\n{cite_data(obs_meta)}")
print(f"review status: {obs_meta['review']['verdict']}")'''))

cells.append(md(r"""## PyMRM implementation

**These are algebraic closures and it would be dishonest to dress them
otherwise.** There is no field to discretise here, no `construct_grad`, no
Newton solve; four correlations are four functions, and inventing a PDE to
justify importing pymrm's operators would obscure what the comparison actually
shows. This page follows [`F1.4`](../F1.4-krishna-ellenberger-holdup/) in saying
so plainly.

What the closures are *for* is to appear as the source term in a momentum
balance, and that is worth doing once, at the end, because it converts a spread
in a friction factor into a spread in a number a designer commits to.

So what follows is four functions on one signature."""))

cells.append(code('''def ergun_fv(x, k1=None, k2=None):
    """Ergun (1952) eq. (13c): f_v = k1 + k2 * N_Re/(1-eps)."""
    k1 = K1_ERGUN if k1 is None else k1
    k2 = K2_ERGUN if k2 is None else k2
    return k1 + k2 * np.asarray(x, float)


def kozeny_carman_fv(x, k=5.0):
    """Carman's creeping-flow limit. f_v = 6^(3-n) k with n = 1, i.e. 36k.

    Constant by construction: no inertial term exists. `x` is accepted only so
    that all four correlations share one signature."""
    return np.full_like(np.asarray(x, float), 36.0 * k)


def forchheimer_fv(x, A, B):
    """Darcy-Forchheimer, f_v = A + B x, with A and B free.

    Identical in form to Ergun; the point is that A and B are fitted to the bed
    in hand rather than taken as universal.  K = eps^3 d_p^2 / (A (1-eps)^2)
    and beta = B (1-eps) / (eps^3 d_p) recover the permeability and the inertial
    coefficient."""
    return A + B * np.asarray(x, float)


def eisfeld_fv(x, eps, D_over_dp, K1=None, k1=None, k2=None):
    """Eisfeld & Schnitzlein (2001) eqs. (3)-(4), sphere coefficients.

    f_v = K1 A_w^2 + (A_w/B_w) * N_Re/(1-eps).
    D_over_dp -> infinity recovers f_v = K1 + x/k2^2."""
    K1 = EIS["K1"] if K1 is None else K1
    k1 = EIS["k1"] if k1 is None else k1
    k2 = EIS["k2"] if k2 is None else k2
    eps = np.asarray(eps, float)
    r = np.asarray(D_over_dp, float)
    A_w = 1.0 + 2.0 / (3.0 * r * (1.0 - eps))
    B_w = (k1 / r ** 2 + k2) ** 2
    return K1 * A_w ** 2 + (A_w / B_w) * np.asarray(x, float)


def to_fk(fv, x):
    """Ergun eq. (14a) from eq. (13a): f_k = f_v / [N_Re/(1-eps)]."""
    return np.asarray(fv, float) / np.asarray(x, float)


print("the four, in the infinite-bed limit (D/d_p -> inf, eps = 0.40):")
print(f"   Ergun            f_v = {K1_ERGUN:5.1f} + {K2_ERGUN:5.3f} x")
print(f"   Kozeny-Carman    f_v = {kozeny_carman_fv(0.0):5.1f}          "
      "(no inertial term at all)")
a_inf = eisfeld_fv(0.0, 0.40, 1e9)
b_inf = eisfeld_fv(1.0, 0.40, 1e9) - a_inf
print(f"   Eisfeld          f_v = {a_inf:5.1f} + {b_inf:5.3f} x")
print(f"   Forchheimer      f_v =     A +     B x   (both free)")
print(f"\\n   1/k2^2 = {1/EIS['k2']**2:.4f}, which is where Eisfeld's inertial "
      "coefficient lands with no wall")'''))

cells.append(md("""## Results

Ergun's Figure 7 rebuilt from the extracted points, panel by panel, with all
four correlations on it. The Eisfeld curves are drawn at the $D/d_p$ that
Eisfeld & Schnitzlein's Table 1 records for each panel's source."""))

cells.append(code('''# representative D/d_p and voidage per panel, from Eisfeld's Table 1
PANEL = {
    "present_investigation": dict(label="Ergun, present investigation",
                                  D_over_dp=None, eps=0.40),
    "burke_plummer": dict(label="Burke & Plummer (1928)",
                          D_over_dp=None, eps=0.39),
    "morcom": dict(label="Morcom (1946)", D_over_dp=None, eps=0.39),
}
for key in ("burke_plummer", "morcom"):
    rows = src[src.ergun_fig7_panel == key]
    PANEL[key]["D_over_dp"] = (float(rows.D_over_dp_min.min()),
                               float(rows.D_over_dp_max.max()))
    PANEL[key]["eps"] = float((rows.eps_min.min() + rows.eps_max.max()) / 2)

xx = np.logspace(np.log10(7), np.log10(4000), 300)
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), sharey=True)
for ax, (key, meta) in zip(axes, PANEL.items()):
    g = obs[obs.series == key]
    ax.loglog(g.Re_over_one_minus_eps, g.f_k, "o", ms=4.5, mfc="none", mew=1.1,
              color="tab:blue", label=f"digitised ({len(g)})")
    ax.loglog(xx, to_fk(ergun_fv(xx), xx), "k-", lw=2.0, label="Ergun eq. (14b)")
    ax.loglog(xx, to_fk(kozeny_carman_fv(xx), xx), "--", lw=1.4,
              color="tab:green", label=r"Kozeny–Carman, $f_v=180$")
    ax.axhline(K2_ERGUN, ls=":", lw=1.4, color="tab:red",
               label=r"Burke–Plummer, $f_k=1.75$")
    if meta["D_over_dp"] is not None:
        lo, hi = meta["D_over_dp"]
        ax.fill_between(xx,
                        to_fk(eisfeld_fv(xx, meta["eps"], hi), xx),
                        to_fk(eisfeld_fv(xx, meta["eps"], lo), xx),
                        color="tab:orange", alpha=0.25, lw=0,
                        label=f"Eisfeld, $D/d_p$ {lo:.1f}–{hi:.0f}")
    else:
        ax.loglog(xx, to_fk(eisfeld_fv(xx, meta["eps"], 30.0), xx), "-.",
                  lw=1.4, color="tab:orange", label=r"Eisfeld, $D/d_p=30$")
    ax.set(xlabel=r"$N_{Re}/(1-\\epsilon)$", title=meta["label"],
           ylim=(0.8, 40))
    ax.legend(fontsize=7.2, loc="lower left")
axes[0].set_ylabel(r"$f_k$  (Ergun eq. 14a)")
fig.tight_layout()
plt.show()

print("Ergun's own eq. (13c) against the 244 points he drew it through:")
dev_ergun = ergun_fv(x) / f_v - 1
print(f"   bias {dev_ergun.mean()*100:+5.2f} %   mean |dev| "
      f"{np.abs(dev_ergun).mean()*100:5.2f} %   "
      f"rms {np.sqrt((dev_ergun**2).mean())*100:5.2f} %")
print("   (deviation = (model - measured)/measured, throughout)")'''))

cells.append(md(r"""### Where the four actually differ

On log axes everything looks like agreement. The honest picture is the ratio."""))

cells.append(code('''fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.3))

ax = axes[0]
ax.loglog(xx, ergun_fv(xx), "k-", lw=2.0, label="Ergun")
ax.loglog(xx, kozeny_carman_fv(xx), "--", lw=1.6, color="tab:green",
          label="Kozeny–Carman (180)")
ax.loglog(xx, eisfeld_fv(xx, 0.40, 1e9), "-.", lw=1.6, color="tab:orange",
          label=r"Eisfeld, $D/d_p\\to\\infty$")
ax.loglog(xx, eisfeld_fv(xx, 0.40, 5.0), ":", lw=1.8, color="tab:orange",
          label=r"Eisfeld, $D/d_p=5$")
ax.loglog(x, f_v, ".", ms=3, color="tab:blue", alpha=0.55,
          label=f"Ergun's data ({len(x)})")
ax.set(xlabel=r"$N_{Re}/(1-\\epsilon)$", ylabel=r"$f_v$",
       title="the linear form Ergun fitted in")
ax.legend(fontsize=8, loc="upper left")

ax = axes[1]
base = ergun_fv(xx)
ax.semilogx(xx, kozeny_carman_fv(xx) / base - 1, "--", lw=1.6,
            color="tab:green", label="Kozeny–Carman")
ax.semilogx(xx, eisfeld_fv(xx, 0.40, 1e9) / base - 1, "-.", lw=1.6,
            color="tab:orange", label=r"Eisfeld, $D/d_p\\to\\infty$")
for r, ls in ((20.0, (0, (4, 2))), (10.0, (0, (2, 2))), (5.0, ":")):
    ax.semilogx(xx, eisfeld_fv(xx, 0.40, r) / base - 1, ls=ls, lw=1.6,
                color="tab:purple", label=rf"Eisfeld, $D/d_p={r:g}$")
ax.axhline(0, color="k", lw=1.6)
ax.axhspan(-0.0516, 0.0516, color="tab:blue", alpha=0.12, lw=0)
ax.set(xlabel=r"$N_{Re}/(1-\\epsilon)$",
       ylabel="(correlation $-$ Ergun) / Ergun", ylim=(-0.55, 0.85),
       title="divergence from Ergun, with his own scatter shaded")
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
plt.show()

for xq in (10, 100, 1000):
    e = ergun_fv(xq)
    print(f"at N_Re/(1-eps) = {xq:5d}:  Ergun {e:8.1f}"
          f"   K–C {kozeny_carman_fv(xq)/e-1:+7.1%}"
          f"   Eisfeld(inf) {eisfeld_fv(xq,0.40,1e9)/e-1:+7.1%}"
          f"   Eisfeld(D/d=5) {eisfeld_fv(xq,0.40,5.0)/e-1:+7.1%}")'''))

cells.append(md(r"""### The wall-effect regime

The correction is not monotone in Reynolds number, and that is the whole point
of Eisfeld & Schnitzlein's paper: *"the wall effect is predicted to increase the
pressure drop at low and moderate Reynolds numbers; only for sufficiently high
$Re_{d_p}$ values and for very small tube-to-particle diameter ratios a pressure
drop reduction is obtained."* Wall friction and the looser near-wall packing
pull opposite ways, and which wins depends on the flow regime.

The shaded bands are where the sources Ergun plotted actually sit."""))

cells.append(code('''ratios = np.logspace(np.log10(1.624), np.log10(250), 400)
fig, ax = plt.subplots(figsize=(8.2, 4.8))
for xq, col in zip((3, 30, 300, 3000, 30000),
                   plt.cm.viridis(np.linspace(0.10, 0.90, 5))):
    ax.semilogx(ratios, eisfeld_fv(xq, 0.40, ratios)
                / eisfeld_fv(xq, 0.40, 1e9) - 1, color=col, lw=2.0,
                label=rf"$N_{{Re}}/(1-\\epsilon) = {xq}$")
ax.axhline(0, color="k", lw=1.2)
for key, col, ytxt in (("morcom", "tab:red", 1.02),
                       ("burke_plummer", "tab:blue", 0.90)):
    lo, hi = PANEL[key]["D_over_dp"]
    ax.axvspan(lo, hi, color=col, alpha=0.09, lw=0)
    ax.plot([lo, hi], [ytxt, ytxt], lw=3, color=col, alpha=0.7,
            solid_capstyle="butt")
    ax.text(hi * 1.15, ytxt, PANEL[key]["label"].split(" (")[0],
            va="center", fontsize=8, color=col)
ax.axvline(10, color="grey", lw=1.0, ls="--")
ax.text(10.4, -0.30, r"$D/d_p=10$: Ergun's own cut" "\\n" "for Burke & Plummer",
        fontsize=7.5, color="grey")
ax.set(xlabel=r"tube-to-particle diameter ratio $D/d_p$",
       ylabel="wall correction on $f_v$, relative to an infinite bed",
       title=r"Eisfeld & Schnitzlein eq. (3), $\\epsilon=0.40$",
       ylim=(-0.45, 1.15), xlim=(1.5, 300))
ax.legend(fontsize=8, loc="center right")
fig.tight_layout()
plt.show()

print("Eisfeld & Schnitzlein, Table 1, for the sources Ergun plotted:")
print(src[["source", "particles", "n_data", "D_over_dp_min", "D_over_dp_max",
           "ergun_fig7_panel"]].fillna("-").to_string(index=False))
sph = src[(src.source.str.startswith("Morcom")) & (src.particles != "cylinders")]
print(f"\\nMorcom's spheres, tablets and nodules sit at D/d_p = "
      f"{sph.D_over_dp_min.min():.1f}-{sph.D_over_dp_max.max():.1f}; only his")
print("cylinders reach 29.6. On Ergun's own criterion - he dropped the Burke &")
print("Plummer systems below D/d_p = 10 - most of that panel is inside the")
print("wall-effect regime, and he kept it.")
print(f"\\nAt the small-D/d_p end the correction changes sign with flow regime:")
for xq in (3, 300, 30000):
    r0 = eisfeld_fv(xq, 0.40, 1.624) / eisfeld_fv(xq, 0.40, 1e9) - 1
    print(f"   D/d_p = 1.62, N_Re/(1-eps) = {xq:6d}: {r0:+7.1%}")
print("   which is the counteracting pair Eisfeld & Schnitzlein describe -")
print("   wall friction raising the viscous term, looser near-wall packing")
print("   lowering the inertial one.")'''))

cells.append(md(r"""## Validation

Six checks. The first is the one the paper pays for: Ergun's constants were
fitted to these points, so refitting them has to give 150 and 1.75 back."""))

cells.append(code('''def fit_forchheimer(x, fv):
    """Least squares on f_v = A + B x, weighted by 1/f_v.

    The weighting is not cosmetic. f_v spans two and a half decades here, so an
    unweighted fit in these coordinates is dominated by a handful of high-Re
    points and returns nonsense; equal RELATIVE weight is what "least squares"
    has to mean when the data are read off log axes."""
    x = np.asarray(x, float); fv = np.asarray(fv, float)
    A = np.vstack([np.ones_like(x), x]).T
    w = 1.0 / fv
    return np.linalg.lstsq(A * w[:, None], fv * w, rcond=None)[0]


print("1. Refit Ergun's two constants from his own figure")
A_all, B_all = fit_forchheimer(x, f_v)
print(f"   all {len(x)} points   k1 = {A_all:7.2f}  ({A_all/K1_ERGUN-1:+.1%})"
      f"     k2 = {B_all:6.3f}  ({B_all/K2_ERGUN-1:+.1%})")
print(f"   Ergun printed  k1 = {K1_ERGUN:7.2f}                 "
      f"k2 = {K2_ERGUN:6.3f}")
print("\\n   and independently, panel by panel - three sources, three separate")
print("   axis calibrations, no shared arithmetic:")
per_panel = {}
for s, g in obs.groupby("series"):
    a, b = fit_forchheimer(g.Re_over_one_minus_eps, g.f_v)
    per_panel[s] = (a, b)
    print(f"      {s:22s} n={len(g):3d}   k1 = {a:6.2f}   k2 = {b:5.3f}")
k1s = np.array([v[0] for v in per_panel.values()])
print(f"   spread in k1 across the three panels: "
      f"{k1s.min():.1f} to {k1s.max():.1f}, i.e. +/-{np.ptp(k1s)/2/k1s.mean():.1%}")
print("\\n   The weighting is a real choice, so here is the alternative:")
Au, Bu = np.linalg.lstsq(np.vstack([np.ones_like(x), x]).T, f_v, rcond=None)[0]
print(f"      unweighted   k1 = {Au:6.2f} ({Au/K1_ERGUN-1:+.1%})   "
      f"k2 = {Bu:5.3f} ({Bu/K2_ERGUN-1:+.1%})")
print(f"      1/f_v weight k1 = {A_all:6.2f} ({A_all/K1_ERGUN-1:+.1%})   "
      f"k2 = {B_all:5.3f} ({B_all/K2_ERGUN-1:+.1%})")
print("   Unweighted, the top decade of f_v carries a hundred times the weight")
print("   of the bottom one, and k1 - which lives at the bottom - drifts. The")
print("   relative-weight fit is the one quoted; the other is shown so that the")
print("   round trip is not resting on a hidden choice.")'''))

cells.append(code('''print("2. The round trip is not circular - the detector cannot see the drawn line")
g = obs[obs.series == "present_investigation"]
print(f"   Ergun's eq. (14b) line is drawn across the whole panel, out to")
print(f"   N_Re/(1-eps) = 4000. The markers stop at "
      f"{g.Re_over_one_minus_eps.max():.1f}.")
print("   Beyond that the panel contains the line, the grid and nothing else,")
print("   and the extraction returns zero points there. Detection keys on the")
print("   enclosed white interior of an open circle; a line has none.")
print("   If the detector had been tracking the curve, this column would run")
print("   to 4000.")

print("\\n3. The same runs, digitised twice, on two different ordinates")
f5x = obs5.Re_over_one_minus_eps.values
f5v = obs5.f_v.values
f7x = g.Re_over_one_minus_eps.values
f7v = g.f_v.values
bins = np.logspace(np.log10(7), np.log10(140), 9)
rat = []
for lo, hi in zip(bins[:-1], bins[1:]):
    m5 = (f5x >= lo) & (f5x < hi)
    m7 = (f7x >= lo) & (f7x < hi)
    if m5.sum() >= 3 and m7.sum() >= 3:
        rat.append(np.median(f7v[m7]) / np.median(f5v[m5]))
rat = np.array(rat)
print(f"   Fig. 5 top panel: {len(obs5)} markers, ordinate f_v, its own"
      " calibration")
print(f"   Fig. 7 top panel: {len(g)} markers, ordinate f_k, its own"
      " calibration")
print(f"   binned median f_v, Fig.7 / Fig.5, over {len(rat)} bins:"
      f"  {rat.mean():.4f} +/- {rat.std():.4f}")
print(f"   -> the two independent readings of one set of runs agree to"
      f" {abs(rat.mean()-1)*100:.1f} %, with no trend")
print("   Point-by-point pairing is NOT attempted: in the dense chain two")
print("   markers can sit closer in x than the pairing tolerance, and a")
print("   mismatched pair inflates the scatter without meaning anything.")'''))

cells.append(code('''print("4. Axis calibration, from the printed graph-paper ladder")
acq = obs_meta["acquisition"]
print("   " + " ".join(acq["estimated_error"].split()))

print("\\n5. The line Ergun labelled 'Kozeny-Carman' is his own 150, not Carman's 180")
print(f"   measured level of the dashed line on Fig. 6: f_v = {KC_DRAWN:.1f}"
      f" +/- {P['kozeny_carman_line_fv_uncertainty']:.1f}")
print(f"   his own k1                                    = {K1_ERGUN:.1f}"
      f"   ({KC_DRAWN/K1_ERGUN-1:+.1%})")
print(f"   Carman's 36k with k = {P['carman_kozeny_constant_k']:.0f}"
      f"                       = {KC_FV:.1f}   ({KC_DRAWN/KC_FV-1:+.1%})")
print("   On that ordinate 180 sits 52 px above 150, so this is not a close call.")
print("   Ergun is drawing the viscous ASYMPTOTE of his own equation and naming")
print("   it after the theory it reproduces - which is worth knowing before")
print("   quoting 'the Kozeny-Carman constant' from his figure.")

print("\\n6. The two limits the equation must have")
lo, hi = 1e-6, 1e9
print(f"   Re/(1-eps) -> 0 :  f_v -> {ergun_fv(lo):.6f}   "
      f"(Kozeny-Carman form, constant f_v)")
print(f"   Re/(1-eps) -> inf: f_k -> {to_fk(ergun_fv(hi), hi):.6f}   "
      f"(Burke-Plummer form, constant f_k)")
print(f"   exact by construction: |f_v(0) - k1| = {abs(ergun_fv(0.0)-K1_ERGUN):.1e},"
      f"  |f_k(inf) - k2| = {abs(to_fk(ergun_fv(hi), hi)-K2_ERGUN):.1e}")
print("   and Eisfeld reduces to an Ergun-type equation with no wall:")
print(f"   f_v(D/d_p -> inf) = {eisfeld_fv(0.0, 0.4, 1e9):.4f}"
      f" + {eisfeld_fv(1.0,0.4,1e9)-eisfeld_fv(0.0,0.4,1e9):.6f} x,"
      f"  vs K1 = {EIS['K1']:.0f} and 1/k2^2 = {1/EIS['k2']**2:.6f}")'''))

cells.append(code('''print("7. How the four score, here and on a database fifty times larger")
sigma_here = {}
for name, model in (
        ("Ergun (150, 1.75)", ergun_fv(x)),
        ("Kozeny-Carman (180)", kozeny_carman_fv(x)),
        ("Forchheimer, refitted here", forchheimer_fv(x, A_all, B_all)),
        ("Eisfeld, D/d_p -> inf", eisfeld_fv(x, 0.40, 1e9))):
    d = model / f_v - 1
    sigma_here[name] = float(np.sqrt((d ** 2).mean()))
    print(f"   {name:28s} bias {d.mean()*100:+6.2f} %   "
          f"mean|dev| {np.abs(d).mean()*100:5.2f} %   "
          f"rms {sigma_here[name]*100:5.2f} %")
print("   (rms here is Eisfeld's sigma: sqrt(mean of squared relative deviation),")
print("    so it includes the bias and is directly comparable with their numbers)")
print("\\n   Eisfeld & Schnitzlein's own relative-rms deviations, over 2391")
print("   points spanning D/d_p = 1.6 to 250 (their Tables 2 and 3, spheres):")
print(f"      Ergun (1952)                {E[('sigma_ergun1952','spheres')]:.4f}")
print(f"      Carman (1937)               {E[('sigma_carman1937','spheres')]:.4f}")
print(f"      Reichelt (1972), wall-corr. {E[('sigma_reichelt1972','spheres')]:.4f}")
print(f"      the same form refitted      {E[('sigma_rel_rms','spheres')]:.4f}"
      "   <- their eq. (3), K1 = 154")
print(f"\\n   Ergun's equation describes Ergun's own figure to "
      f"{sigma_here['Ergun (150, 1.75)']*100:.1f} % rms and everyone else's beds")
print(f"   to {E[('sigma_ergun1952','spheres')]*100:.1f} %. The factor of "
      f"{E[('sigma_ergun1952','spheres')]/sigma_here['Ergun (150, 1.75)']:.1f} is the "
      "price of extrapolating a fit,")
print("   and it is the honest error bar to carry into a design.")

report_agreement("A1.1", {
    "k1_refit_all_points": float(A_all),
    "k2_refit_all_points": float(B_all),
    "k1_relative_error": float(A_all / K1_ERGUN - 1),
    "k2_relative_error": float(B_all / K2_ERGUN - 1),
    "ergun_bias_on_own_data": float(dev_ergun.mean()),
    "ergun_mean_abs_dev_on_own_data": float(np.abs(dev_ergun).mean()),
    "ergun_rms_on_own_data": float(sigma_here["Ergun (150, 1.75)"]),
    "fig5_fig7_cross_check_ratio": float(rat.mean()),
    "kozeny_carman_line_over_k1": float(KC_DRAWN / K1_ERGUN),
    "n_markers": int(len(obs)),
})'''))

cells.append(md(r"""## What pymrm adds

**Nothing to the correlations, and the honest thing is to say so.** They are
four lines of algebra. pymrm has no operator to contribute and none is used.

What this page adds is quantitative, and it is of three kinds.

**A round trip nobody had done.** Ergun's 150 and 1.75 are quoted everywhere and
the 640 experiments behind them have been unreadable since 1952 — the paper is a
pure scan with no text layer and no tables. Recovering 244 of the points and
refitting gives $k_1 = 151.9$ and $k_2 = 1.70$, within 1.3 % and 3.0 % of the
printed values, from three sources fitted separately and agreeing among
themselves to ±3 %. That is a check on the digitisation and on the paper at
once, and it is the reason the rest of the page can be believed.

**A number for how far the four diverge, instead of a qualitative ranking.**
Kozeny–Carman is 7.5 % high at $N_{Re}/(1-\epsilon)=10$, 45 % low at 100 and a
factor of ten low at 1000 — the creeping-flow limit fails much earlier than its
reputation suggests. Eisfeld's infinite-bed form sits 2.7 % above Ergun in the
viscous limit and 25 % below in the inertial one: the two most-used correlations
in the field disagree by a quarter at high Reynolds number, and it is the
*turbulent* constant they disagree about, not the famous 150. And once
$D/d_p < 10$ the wall correction is larger than the gap between any two of them.

**A regime where Ergun's own figure contains the effect it cannot show.**
Morcom's spheres, tablets and nodules sit at $D/d_p = 3.7$–$10.6$. Ergun dropped
exactly that range from Burke & Plummer's data — journal page 93, "those for
which the ratio of tube diameter to particle size was less than 10 … have been
omitted" — and kept it in Morcom's. It is the panel with the largest scatter of
the three.

And one small correction to the folklore, which the page can make because it
measured the figure rather than reading about it: the line labelled
"Kozeny–Carman" on Ergun's Figure 6 is drawn at $f_v = 149$, i.e. at *his* 150,
not at Carman's 180."""))

cells.append(code('''# Where the closure is consumed: an isothermal compressible bed.
# With ideal gas at fixed mass flux both Ergun terms scale as 1/p, so
#   p dp/dz = -K,   K = (RT/M) [ A (1-eps)^2 mu G / (eps^3 dp^2)
#                              + B (1-eps) G^2 / (eps^3 dp) ]
# and p(z) = sqrt(p_in^2 - 2 K z) exactly.  The two-route check below is the
# point: the closed form and a direct ODE integration of dp/dz = -K/p are
# independent, and must agree.
from scipy.integrate import solve_ivp

BED = dict(T=300.0, M=0.029, mu=1.85e-5, dp=3e-3, eps=0.40, L=2.0, p_in=2e5)


def bed_K(A, B, G, T=None, M=None, mu=None, dp=None, eps=None, **_):
    T = BED["T"] if T is None else T
    M = BED["M"] if M is None else M
    mu = BED["mu"] if mu is None else mu
    dp = BED["dp"] if dp is None else dp
    eps = BED["eps"] if eps is None else eps
    return (R_GAS * T / M) * (A * (1 - eps) ** 2 * mu * G / (eps ** 3 * dp ** 2)
                              + B * (1 - eps) * G ** 2 / (eps ** 3 * dp))


def outlet_pressure(A, B, G):
    K = bed_K(A, B, G)
    return np.sqrt(BED["p_in"] ** 2 - 2 * K * BED["L"])


# route 2: integrate dp/dz = -K/p numerically, no closed form used
K_chk = bed_K(K1_ERGUN, K2_ERGUN, 1.0)
sol = solve_ivp(lambda z, p: -K_chk / p, (0.0, BED["L"]), [BED["p_in"]],
                rtol=1e-10, atol=1e-6, dense_output=True)
print(f"closed form  p_out = {outlet_pressure(K1_ERGUN, K2_ERGUN, 1.0):.4f} Pa")
print(f"solve_ivp    p_out = {sol.y[0,-1]:.4f} Pa")
print(f"agreement           {abs(sol.y[0,-1]/outlet_pressure(K1_ERGUN,K2_ERGUN,1.0)-1):.2e}"
      "   <- two independent routes, so the algebra is not the weak link")

# Stop the sweep before the bed would need more head than the inlet has:
# 2 K(G) L = p_in^2 with K = c1 G + c2 G^2 is a quadratic in G.
c1 = bed_K(K1_ERGUN, 0.0, 1.0)
c2 = bed_K(0.0, K2_ERGUN, 1.0)
G_max = (-c1 + np.sqrt(c1 ** 2 + 4 * c2 * BED["p_in"] ** 2 / (2 * BED["L"]))) \\
    / (2 * c2)
Gs = np.logspace(-1.5, np.log10(0.85 * G_max), 60)

fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.3))
for ax, ratio, ttl in zip(axes, (40.0, 5.0),
                          (r"a wide tube, $D/d_p=40$",
                           r"a narrow tube, $D/d_p=5$")):
    a_e = eisfeld_fv(0.0, BED["eps"], ratio)
    b_e = eisfeld_fv(1.0, BED["eps"], ratio) - a_e
    curves = (("Ergun", (K1_ERGUN, K2_ERGUN), "k", "-"),
              ("Kozeny–Carman", (KC_FV, 0.0), "tab:green", "--"),
              ("Eisfeld", (a_e, b_e), "tab:orange", "-."))
    for nm, (A, B), col, ls in curves:
        dp_bar = (BED["p_in"] - outlet_pressure(A, B, Gs)) / 1e5
        ax.loglog(Gs * BED["dp"] / BED["mu"], dp_bar, ls=ls, color=col,
                  lw=1.8, label=nm)
    ax.set(xlabel=r"$Re_{d_p} = G d_p/\\mu$",
           ylabel=r"$\\Delta p$ over 2 m  [bar]", title=ttl)
    ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

print("\\noutlet pressure after 2 m, air at 300 K, 3 mm particles, eps = 0.40,")
print(f"inlet {BED['p_in']/1e5:.1f} bar, G = 1.0 kg/(m2 s)  "
      f"(Re_dp = {1.0*BED['dp']/BED['mu']:.0f}):")
for nm, (A, B) in (("Ergun", (K1_ERGUN, K2_ERGUN)),
                   ("Kozeny–Carman", (KC_FV, 0.0)),
                   ("Eisfeld, D/d_p = 40",
                    (eisfeld_fv(0.0, 0.4, 40.0),
                     eisfeld_fv(1.0, 0.4, 40.0) - eisfeld_fv(0.0, 0.4, 40.0))),
                   ("Eisfeld, D/d_p = 5",
                    (eisfeld_fv(0.0, 0.4, 5.0),
                     eisfeld_fv(1.0, 0.4, 5.0) - eisfeld_fv(0.0, 0.4, 5.0)))):
    po = outlet_pressure(A, B, 1.0)
    print(f"   {nm:22s} dp = {(BED['p_in']-po)/1e3:5.2f} kPa   "
          f"p_out = {po/1e3:6.2f} kPa")
print("\\n   At this Reynolds number the wall's two effects nearly cancel:")
print("   D/d_p = 5 raises the viscous coefficient by 49 % and lowers the")
print("   inertial one by 8 % relative to Eisfeld's infinite bed, and the bed")
print("   is inertia-dominated, so the net is small. Drop the flow by a decade")
print("   and the same correction is worth tens of percent - which is why the")
print("   wall effect is reported as an increase by some authors and a decrease")
print("   by others.")'''))

cells.append(md(r"""**The honest limits of this page.**

**Incomplete recall, and it is a limitation of the dataset rather than a
footnote to it.** The 244 points are not Ergun's 640. A maintainer has reviewed
the numbered overlays and the verdict is recorded in the sidecar: the ringed
centres are on real markers — nothing was picked off the drawn eq. (14b) line,
the graph paper or the panel labels, including in the dense chain where markers
and line merge — but a few markers per panel were missed. The mechanism is the
detection rule itself: a marker is found by the enclosed white region inside its
ring, so where several markers overlap into a solid blob there is no interior
left to find, and the low-Reynolds half of the "Present investigation" panel is
largely such a chain. Figures 6 and 8 carry further sources (Oman & Watson) that
were not attempted at all.

*What that affects.* Precision, and any count. Every constant fitted on this
page is fitted to a sample, so it carries a wider interval than 244 independent
points would suggest; and the number of markers in a region of the plot is not
the number of experiments Ergun ran there, because dense regions lose
proportionally more. The panel-level $k_2$ for the present-investigation panel
is poorly determined for a second and unrelated reason — its markers stop at
$N_{Re}/(1-\epsilon) = 136$, so there is little inertial range to fit.

*What it does not affect.* The axis calibration, the coordinates of any marker
that was recovered, and — because what was reported is a few scattered misses
rather than a run of them — the central values. Random misses add variance; they
do not move a fit. Only a systematic loss would bias it, and two independent
things argue against one here: the refit returns Ergun's own constants to within
a few percent, and three panels calibrated separately agree on $k_1$ to ±3 %, so
any systematic loss would have to be common to three independent digitisations.
The one failure mode that could have produced exactly that — a series drawn with
crosses or triangles, which enclose nothing and would vanish as a group — was
excluded before extraction by inspecting Figure 7 at 4× magnification: it
carries no shape legend and every marker in it is an open circle. Figure 5's
Morcom panel does mix four symbols, which is why the question was asked.

**$D/d_p$ is attached to sources, never to points.** Ergun's axes do not carry
it. The wall-effect regime is identified from Eisfeld & Schnitzlein's
re-cataloguing of the same primary references, and the two extractions need not
contain the same runs — Ergun dropped the mixtures and the $D/d_p<10$ systems
from Burke & Plummer, and Eisfeld's Burke & Plummer row still reaches down to
5.4.

**No experimental test of the wall correction is possible here.** Eisfeld's
equation is evaluated and compared; it is not validated, because that would need
the $D/d_p$ of each individual run."""))

cells.append(md(r"""## Reuse

**The four correlations are standalone.** `ergun_fv`, `kozeny_carman_fv`,
`forchheimer_fv` and `eisfeld_fv` all return $f_v$ against
$N_{Re}/(1-\epsilon)$ and depend on nothing else on this page. `to_fk` converts
to the Blake-type friction factor. To get back to a pressure gradient:

$$
\frac{\Delta P}{L} = f_v \cdot \frac{\mu U_m}{D_p^2}\,
\frac{(1-\epsilon)^2}{\epsilon^3}\Big/ g_c ,
$$

with $g_c = 1$ in SI. `bed_K` and `outlet_pressure` do the compressible-bed
integration for any $(A, B)$ pair, including a Forchheimer $(A,B)$ fitted to a
bed of your own.

**Which one to use.**

- $D/d_p > 10$ and you want one number: **Ergun**, and carry ±20 % — that is
  Eisfeld & Schnitzlein's measured relative rms for it over 2391 points, not a
  rule of thumb.
- $D/d_p < 10$: **Eisfeld–Schnitzlein**, spheres coefficients
  $K_1=154,\,k_1=1.15,\,k_2=0.87$; cylinders $190,\,2.00,\,0.77$; mixed
  $155,\,1.42,\,0.83$. Valid for
  $0.01 \le Re_{d_p} \le 17635$, $1.6 \le D/d_p \le 250$,
  $0.33 \le \epsilon \le 0.88$.
- Creeping flow only, $N_{Re}/(1-\epsilon) \lesssim 3$: **Kozeny–Carman** is
  within a few percent of Ergun and simpler. Above that it degrades fast.
- No particle diameter — foam, monolith, rock, fibre mat: **Darcy–Forchheimer**
  with $K$ and $\beta$ measured on the medium. `forchheimer_fv` plus
  `fit_forchheimer` is the whole procedure; note the $1/f_v$ weighting, which is
  not optional when the data span decades.

**Voidage is the sensitivity that dominates everything above.** $f_v$ carries
$(1-\epsilon)^2/\epsilon^3$; at $\epsilon = 0.40$ a 5 % error in $\epsilon$ moves
the viscous coefficient by about 24 %. Any of these correlations with a guessed
$\epsilon$ is worse than the poorest of them with a measured one.

**Related pages.** `A1.5` Richardson–Zaki and `A1.6` Wen–Yu minimum
fluidisation take the same bed to the point where it lifts;
[`A2.3`](../A2.3-taylor-aris-dispersion/) is the dispersion counterpart of this
`S3` closure; `E2.1` needs $u_{mf}$, which needs this.

**Cite the sources, not this page:** Ergun, S., *Fluid flow through packed
columns*, Chemical Engineering Progress **48**(2) 89–94 (1952). Eisfeld, B. and
Schnitzlein, K., *The influence of confining walls on the pressure drop in packed
beds*, Chemical Engineering Science **56**(14) 4321–4329 (2001),
[doi:10.1016/S0009-2509(00)00533-9](https://doi.org/10.1016/S0009-2509(00)00533-9).
Carman, P. C., *Fluid flow through granular beds*, Trans. Inst. Chem. Eng.
**15** 150–166 (1937)."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python",
                   "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py", "mimetype": "text/x-python", "name": "python",
        "nbconvert_exporter": "python", "pygments_lexer": "ipython3",
        "version": "3.13.5"},
}
out = Path(__file__).with_name("index.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
