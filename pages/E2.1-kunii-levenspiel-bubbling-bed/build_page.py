#!/usr/bin/env python3
"""Generate index.ipynb for page E2.1. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "The bubbling bed, and the conversion ceiling nobody can buy their way past"
description: "Kunii and Levenspiel's 1968 model of a fluidized bed reduces to one algebraic expression. Reproducing their three worked appendices verifies it, and then the same expression says a bed cannot exceed 99.1 % conversion however good the catalyst is."
categories: [sec:E, struct:S7, tier:T0, data:tier6, phase:gas-solid]
date: 2026-07-30
---

# The bubbling bed, and the conversion ceiling nobody can buy their way past

**Catalog ID:** `E2.1` · **Structures:** `S7` (multiphase coupling) · **Tier:** T0

A fluidized bed is the only common reactor that routinely converts **less than a
perfectly stirred tank**. Feed a catalytic bed hard enough and conversion stops
responding to the catalyst altogether: gas rides up inside bubbles, and bubbles
have almost no catalyst in them.

Kunii and Levenspiel's 1968 paper is where that behaviour got a quantitative
model with **one adjustable parameter** — the effective bubble diameter — instead
of the six that competing two-region models needed. This page reproduces it,
checks it against the three worked examples the authors printed, and then asks
the question the model answers but the paper never states: *how good can a
bubbling bed get?*"""))

cells.append(md(r"""## Background

Above minimum fluidization the excess gas travels as **bubbles**. A fast bubble
drags a thin shell of emulsion round with it — the **cloud** — plus a **wake** of
solids behind it, and the whole assembly rises faster than the gas percolating
through the emulsion. Fresh reactant therefore enters the bed almost entirely
inside bubbles, which contain hardly any catalyst, and it can only react after
crossing into the cloud and then into the emulsion.

That gives four resistances in series and parallel:

| step | rate coefficient |
|---|---|
| reaction on the few solids dispersed *in* the bubble | $\gamma_b K_r$ |
| interchange bubble $\to$ cloud+wake | $K_{bc}$ |
| reaction on the solids in the cloud+wake | $\gamma_c K_r$ |
| interchange cloud $\to$ emulsion | $K_{ce}$ |
| reaction on the solids in the emulsion | $\gamma_e K_r$ |

Every one of these follows from **one** measured quantity, the effective bubble
size $d_b$. That is the claim of the paper, and it is what makes the model
usable: the same $d_b$ has to explain gas–solid mass transfer, gas–solid heat
transfer *and* catalytic conversion, and the paper works all three.

Two of the older explanations are worth naming because the model has to beat
them. Tanks-in-series and axial-dispersion models can only produce conversions
**between** plug flow and backmix flow; measured beds routinely fall below
backmix. And a dozen two-region models with up to six parameters could each fit
one data set and no other."""))

cells.append(md(r"""## The published model

Everything below was read off 600 dpi renders of the printed pages. The scan has
the best text layer of the papers in this gallery — 12.4k characters per page —
and it is *still* useless for the equations: it renders equations 49, 50 and 56
as loose fragments and drops decimal points (appendix C's $13.2\,\mathcal{K}_m$
comes out as `1.32 K_m`). Equation numbers are the paper's.

**Flow.** With $u_0 > 2u_{mf}$ only large fast bubbles survive, and

$$
\delta = \frac{u_0-u_{mf}}{u_b},\qquad
u_b = u_0 - u_{mf} + u_{br},\qquad
u_{br} = 0.711\,(g d_b)^{1/2}
\tag{2, 3, 51}
$$

$$
1-\delta = \frac{L_{mf}}{L_f} = \frac{1-\varepsilon_f}{1-\varepsilon_{mf}},
\qquad
\frac{L_{mf}}{L_m} = \frac{1-\varepsilon_m}{1-\varepsilon_{mf}}
\tag{4}
$$

**Interchange**, per unit bubble volume per second:

$$
K_{bc} = 4.5\left(\frac{u_{mf}}{d_b}\right) + 5.85\,\frac{\mathfrak{D}^{1/2}g^{1/4}}{d_b^{5/4}},
\qquad
K_{ce} \cong 6.78\left(\frac{\varepsilon_{mf}\mathfrak{D}_e u_b}{d_b^{3}}\right)^{1/2}
\tag{10, 12}
$$

**Solids distribution**, as volume of solids per volume of bubble:

$$
\gamma_b \cong 0.001\text{–}0.01,
\qquad
\gamma_c = (1-\varepsilon_{mf})\left(\frac{3u_{mf}/\varepsilon_{mf}}{0.711(gd_b)^{1/2}-u_{mf}/\varepsilon_{mf}} + \alpha\right),
$$
$$
\delta(\gamma_b+\gamma_c+\gamma_e) = (1-\varepsilon_{mf})(1-\delta)
\tag{15, 16, 17}
$$

with $\alpha = V_w/V_b$ the wake-to-bubble volume ratio, 0.25–0.4 by equation 6.

**Reaction.** For a first-order reaction the three regions balance as

$$
-u_b\frac{dC_{Ab}}{dl} = (K_r)_b C_{Ab} = \gamma_b K_r C_{Ab} + K_{bc}(C_{Ab}-C_{Ac})
\tag{48a}
$$
$$
K_{bc}(C_{Ab}-C_{Ac}) \cong \gamma_c K_r C_{Ac} + K_{ce}(C_{Ac}-C_{Ae}),
\qquad
K_{ce}(C_{Ac}-C_{Ae}) \cong \gamma_e K_r C_{Ae}
\tag{48b, 48c}
$$

Equations 48b and 48c are **algebraic**: the cloud and the emulsion are assumed
to be in local steady state, with no axial transport of their own. Eliminating
$C_{Ac}$ and $C_{Ae}$ collapses the whole reactor to one number:

$$
\mathcal{K}_f = (K_r)_b\left(\frac{L_f}{u_b}\right)
= K_r\left(\frac{L_f}{u_b}\right)
\left[\gamma_b + \cfrac{1}{\cfrac{K_r}{K_{bc}} + \cfrac{1}{\gamma_c + \cfrac{1}{\cfrac{K_r}{K_{ce}}+\cfrac{1}{\gamma_e}}}}\right]
\tag{49}
$$

$$
K_r\left(\frac{L_f}{u_b}\right) = \frac{1}{1-\varepsilon_{mf}}\left(\frac{u_0}{u_{br}}\right)\mathcal{K}_m,
\qquad
\mathcal{K}_m = (1-\varepsilon_m)K_r L_m/u_0
\tag{50, 45}
$$

and the exit gas — taken to be bubble gas only — gives

$$
1 - X_A = e^{-\mathcal{K}_f}
\tag{54}
$$

against $e^{-\mathcal{K}_m}$ for plug flow (eq. 42) and $1/(1+\mathcal{K}_m)$ for
backmix flow (eq. 46).

**This is algebra, not a differential equation.** Say it plainly: the published
model is a continued fraction, and evaluating it needs numpy and nothing else.
The differential equation only comes back when the assumption behind equations
48b and 48c is relaxed, which is what the pymrm section below does."""))

cells.append(md(r"""### The same structure, three times over

The paper's own consistency argument is worth carrying: the mass-transfer result
(eq. 27), the heat-transfer result (eq. 38) and the fast-reaction limit (eq. 55)
are the *same* expression with $B_d$, $H_{bc}$ or $K_r$ in the reacting slot:

$$
(\mathrm{Sh})_{\text{over-all}} \cong \frac{\delta}{1-\varepsilon_f}
\left[\gamma_b(\mathrm{Sh}^*)_t + \frac{y\varphi_s d_p^2}{6\mathfrak{D}}K_{bc}\right],
\qquad
\mathrm{Nu}_{\text{apparent}} \cong \frac{1}{1-\varepsilon_f}
\left[\gamma_b(\mathrm{Nu}^*)_t + \frac{\varphi_s d_p^2}{6k_g}H_{bc}\right]
$$

The two prefactors differ deliberately — the text says equation 38 was obtained
after "replacing $u_b$ by $u_0/\delta$", which cancels the $\delta$. The local
coefficients come from the textbook single-sphere correlation, equation 1,
$\mathrm{Sh}^* = 2 + 0.6\,\mathrm{Sc}^{1/3}\mathrm{Re}^{1/2}$, evaluated at the
particle's terminal velocity."""))

cells.append(md(r"""## Parameters and assumptions

Three worked examples are printed as appendices, and between them they exercise
every equation above. Their inputs are in the sidecar of the dataset below; the
one used throughout the reaction sections is **appendix C**, Kobayashi and Arai's
ozone-decomposition bed:

| | |
|---|---|
| $\mathfrak{D} = \mathfrak{D}_e$ | 0.204 cm²/s |
| $d_t$, $L_m$ | 20 cm, 34 cm |
| $\varepsilon_m$, $\varepsilon_{mf}$ | 0.45, 0.50 |
| $u_{mf}$ | 2.1 cm/s |
| $u_0$ | 13.2 cm/s |
| $d_b$ | 3.7 cm (estimated — *the* fitted parameter) |
| $\gamma_b$ | 0 |
| $\alpha$ | 0.47, from Rowe and Partridge (1965) |

**Assumptions that matter.** The emulsion stays exactly at minimum fluidizing
conditions at every gas velocity; all bubbles have one size; bubbles are
solid-free apart from $\gamma_b$; the exit stream is bubble gas only. And
$\alpha = 0.47$ here sits *outside* the 0.25–0.40 range equation 6 quotes — the
paper takes it from a different source and says so.

**Two printed slips, recorded rather than repaired.** Appendix C prints
$u_0 = (6.6+9.9+13.2+20)/5 = 13.2$ — four terms over five. The value 13.2 is used
everywhere downstream and is confirmed by the paper's own
$u_b = 13.2-2.1+42.8 = 53.9$, so nothing is in doubt; a term was dropped in
typesetting. Appendix B prints $(1-\varepsilon_f)(u_0+5.8) = 8.70$ where its own
stated $\varepsilon_{mf} = 0.50$ gives 7.87. Both are carried into the dataset as
printed, and the validation section shows what each implies."""))

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
import matplotlib.pyplot as plt
from scipy.sparse import csc_array, identity, kron
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq
from pymrm import construct_div, construct_convflux_upwind
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "E2.1-kunii-levenspiel-bubbling-bed"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

G = 980.0                 # cm/s^2 - the paper works entirely in CGS
RHO_G, MU = 1.18e-3, 1.8e-4      # g/cc, g/(cm s): air, stated in appendices A and B

# The four numerical coefficients Kunii and Levenspiel print inside equations 10,
# 12 and 51. They are named rather than inlined so the validation section can
# perturb them one digit at a time and show which checks notice - see the
# defect-injection table. The defaults are exactly what the paper prints.
KL = dict(c_ubr=0.711, c_kbc_visc=4.5, c_kbc_diff=5.85, c_kce=6.78)

# Appendix C: Kobayashi and Arai's ozone bed
C = dict(u0=13.2, umf=2.1, emf=0.50, em=0.45, Lm=34.0, db=3.7, alpha=0.47,
         D=0.204, gamma_b=0.0)'''))

cells.append(md(r"""## The data

There is **no experimental data on this page**, and that is a deliberate choice
rather than an oversight.

The paper's experimental content is figures 3, 5, 7, 8 and 9 — Sherwood and
Nusselt numbers against Reynolds number, and conversion against
$\mathcal{K}_m$ — all of them scatter plots that would have to be digitised. On
figure 9, the sharpest of them, the markers at low $\mathcal{K}_m$ overlap into a
solid mass in which individual centres cannot be located, so the extraction is
staged for maintainer review rather than published here.

What *is* loaded is a tier-6 reference set: **the numbers Kunii and Levenspiel
printed while working their own three appendices**, read off 600 dpi renders.
Reproducing them tests this reimplementation and the reading of the equations. It
does not test whether the model describes a real bed, and the page does not claim
otherwise."""))

cells.append(code('''ref = load_data("kunii_levenspiel_1968_appendix_values.csv", page=PAGE)
meta = load_meta("kunii_levenspiel_1968_appendix_values.csv", page=PAGE)
print(cite_data(meta))
print()
print(ref.to_string(index=False))

def printed(appendix, symbol):
    """The value as printed in the paper."""
    row = ref[(ref.appendix == appendix) & (ref.symbol == symbol)]
    return float(row.value.iloc[0])'''))

cells.append(md(r"""## PyMRM implementation

Two layers, and they are honestly different in character.

**The published model is closed form.** `hydrodynamics` and `Kf_closed` below are
equations 2, 3, 10, 12, 15, 17, 45, 49, 50 and 51 transcribed. No solver is
involved and none is needed.

**The distributed model is where pymrm earns its place.** Equations 48b and 48c
freeze the cloud and the emulsion: they carry no gas anywhere, they only exchange
sideways. That is what makes the algebra close. But the paper itself flags the
limitation — for $u_0/u_{mf} < 6\sim11$ the emulsion gas flows *upward* and, in
its own words, "the conversion of emulsion gas will have to be found by an
extension of the present analysis". Kobayashi's runs are at
$u_0/u_{mf} = 3.1$–$9.6$, straddling exactly that boundary.

So the second layer solves all three regions as a **1-D convection–exchange
system on the fluidized bed height**, with

* bubble gas at superficial velocity $v_b = u_0 - u_{mf}$ (equation 2, exactly),
* cloud+wake gas carried up **with the bubble** at $u_b$, so
  $v_c = \delta(V_c/V_b+\alpha)\varepsilon_{mf}u_b$,
* emulsion gas at whatever is left, $v_e = u_{mf} - v_c$, which is **negative**
  for these conditions — the emulsion gas flows down.

The last line is a closure choice and is stated as one: the three superficial
velocities are required to sum to $u_0$. The paper's own equations 5 and 7 give
an emulsion velocity too, from the solids circulation, and it does *not* satisfy
that balance (it overshoots $u_0$ by 9.7 % here). The results section runs both.

Downflow makes this a genuinely two-point problem: the emulsion enters at the
**top**, carrying gas that has already been through the bed, and rejoins the feed
at the bottom. That recirculation loop is the flow pattern of the paper's own
figure 6, and it is precisely what the closed form cannot represent."""))

cells.append(code('''def hydrodynamics(u0, umf, emf, em, Lm, db, alpha, D, gamma_b=0.0, coef=None):
    """Equations 2, 3, 4, 10, 12, 15, 17 and 51, in CGS.

    ``coef`` overrides the four printed coefficients of equations 10, 12 and 51.
    The default, KL, is exactly what the paper prints; the validation section
    passes perturbed copies to measure which checks can see a mis-transcription.
    """
    k = KL if coef is None else coef
    ubr = k["c_ubr"] * np.sqrt(G * db)                               # eq. 51
    ub = u0 - umf + ubr                                              # eq. 3
    delta = (u0 - umf) / ub                                          # eq. 2
    Kbc = (k["c_kbc_visc"] * (umf / db)
           + k["c_kbc_diff"] * (D**0.5 * G**0.25 / db**1.25))        # eq. 10
    Kce = k["c_kce"] * np.sqrt(emf * D * ub / db**3)                 # eq. 12
    fc = 3.0 * (umf / emf) / (ubr - umf / emf) + alpha               # V_c/V_b + alpha, eq. 9 + 6
    gc = (1 - emf) * fc                                              # eq. 17
    ge = (1 - emf) * (1 - delta) / delta - gc - gamma_b              # eq. 15
    Lf = Lm * (1 - em) / ((1 - delta) * (1 - emf))                   # eq. 4, twice
    return dict(ubr=ubr, ub=ub, delta=delta, Kbc=Kbc, Kce=Kce, fc=fc,
                gb=gamma_b, gc=gc, ge=ge, Lf=Lf, emf=emf, u0=u0, umf=umf,
                em=em, Lm=Lm, db=db, alpha=alpha, D=D)


def Kr_of_Km(Km, h):
    """Equation 45 read backwards: the intrinsic rate constant behind K_m."""
    return Km * h["u0"] / ((1 - h["em"]) * h["Lm"])


def Kf_closed(Km, h):
    """Equations 49 + 50: the dimensionless rate group for the fluidized bed."""
    Kr = Kr_of_Km(Km, h)
    inner = Kr / h["Kce"] + 1.0 / h["ge"]
    mid = h["gc"] + 1.0 / inner
    brack = h["gb"] + 1.0 / (Kr / h["Kbc"] + 1.0 / mid)
    return (1.0 / (1 - h["emf"])) * (h["u0"] / h["ubr"]) * Km * brack


def conversion_closed(Km, h):
    """Equation 54."""
    return 1.0 - np.exp(-Kf_closed(Km, h))'''))

cells.append(code('''def solve_distributed(Km, h, n=800, convect=True, ve_mode="balance",
                     tol=1e-13, itmax=300):
    """Three-region axially distributed model on the fluidized bed height.

    Fields are ordered (bubble, cloud+wake, emulsion) on the LAST axis, spatial
    axis first: shape (n_z, 3), as the style guide requires.

    convect=False freezes the cloud and emulsion, which is exactly the paper's
    equations 48b and 48c, and the answer must then be equation 54.
    """
    Kr = Kr_of_Km(Km, h)
    d, Lf, emf = h["delta"], h["Lf"], h["emf"]
    vb = h["u0"] - h["umf"]                          # eq. 2, exactly
    vc = d * h["fc"] * emf * h["ub"] if convect else 0.0
    if not convect:
        ve = 0.0
    elif ve_mode == "balance":                       # v_b + v_c + v_e = u_0
        ve = h["umf"] - vc
    else:                                            # the paper's eqs. 5 and 7
        us = h["alpha"] * d * h["ub"] / (1 - d - h["alpha"] * d)
        ue = h["umf"] / emf - us
        ve = emf * (1 - d * (1 + h["fc"])) * ue
    v = np.array([vb, vc, ve])

    shape = (n, 3)
    x_f = np.linspace(0.0, Lf, n + 1)
    x_c = 0.5 * (x_f[:-1] + x_f[1:])
    div = construct_div(shape, x_f, nu=0, axis=0)    # nu=0: Cartesian, constant area

    # Outward normal at both ends. A phase that ENTERS at a face gets Dirichlet
    # there ({a:0, b:1, d:C_in}); a phase that LEAVES gets zero gradient
    # ({a:1, b:0, d:0}). Leaving an outlet as None makes the matrix singular.
    a_l, b_l = np.where(v > 0, 0.0, 1.0), np.where(v > 0, 1.0, 0.0)
    a_r, b_r = np.where(v < 0, 0.0, 1.0), np.where(v < 0, 1.0, 0.0)

    # Exchange + reaction, per unit BED volume: constant, so assemble once.
    S = np.zeros((3, 3))
    S[0, 0] = -d * (h["gb"] * Kr + h["Kbc"]);            S[0, 1] = d * h["Kbc"]
    S[1, 0] = d * h["Kbc"]
    S[1, 1] = -d * (h["Kbc"] + h["gc"] * Kr + h["Kce"]); S[1, 2] = d * h["Kce"]
    S[2, 1] = d * h["Kce"]
    S[2, 2] = -d * (h["Kce"] + h["ge"] * Kr)
    Sfull = kron(identity(n, format="csc"), csc_array(S), format="csc")

    Cin_bot, Ctop = 1.0, 1.0                         # feed concentration = 1
    for it in range(itmax):
        bc = ({"a": a_l, "b": b_l, "d": np.where(v > 0, Cin_bot, 0.0)},
              {"a": a_r, "b": b_r, "d": np.where(v < 0, Ctop, 0.0)})
        conv, conv_bc = construct_convflux_upwind(shape, x_f, x_c, bc=bc, v=v, axis=0)
        M = (div @ conv - Sfull).tocsc()
        rhs = -np.asarray((div @ conv_bc).todense()).ravel()
        Cf = spsolve(M, rhs).reshape(n, 3)
        if ve < 0:   # emulsion recirculates: top mixture in, bottom mixture back to feed
            Ct_new = (vb * Cf[-1, 0] + vc * Cf[-1, 1]) / (vb + vc)
            Cin_new = (h["u0"] * 1.0 + (-ve) * Cf[0, 2]) / (vb + vc)
        else:
            Ct_new, Cin_new = 1.0, 1.0
        done = abs(Ct_new - Ctop) < tol and abs(Cin_new - Cin_bot) < tol
        Ctop, Cin_bot = Ct_new, Cin_new
        if done:
            break

    up = np.maximum(v, 0.0)                          # mixing cup of the upward streams
    Cexit = float(up @ Cf[-1, :] / up.sum())
    return dict(Cexit=Cexit, C=Cf, x=x_c, v=v, iters=it, Cin=Cin_bot)'''))

cells.append(md(r"""## Results"""))

cells.append(code('''hC = hydrodynamics(**C)
print("Appendix C bed, recomputed from equations 2, 3, 10, 12, 15, 17, 51:")
for k in ("ubr", "ub", "delta", "Kbc", "Kce", "gc", "ge", "Lf"):
    print(f"  {k:6s} = {hC[k]:10.4f}")
print(f"  K_r    = {Kr_of_Km(1.0, hC):10.4f} x K_m")

Km = np.linspace(1e-3, 8.0, 400)
fig, ax = plt.subplots(figsize=(6.8, 4.8))
for db, style in ((3.7, "-"), (4.2, "--"), (5.0, "-.")):
    h = hydrodynamics(**{**C, "db": db})
    ax.semilogy(Km, np.exp(-Kf_closed(Km, h)), style, color="C0", lw=1.8,
                label=f"bubbling bed, $d_b$ = {db} cm")
ax.semilogy(Km, np.exp(-Km), ":", color="0.35", lw=1.6, label="plug flow (eq. 42)")
ax.semilogy(Km, 1.0 / (1.0 + Km), color="C3", lw=1.6, label="backmix flow (eq. 46)")
ceiling = hC["Kbc"] * hC["Lf"] / hC["ub"]
ax.axhline(np.exp(-ceiling), color="C1", lw=1.2, ls=(0, (1, 1)))
ax.text(3.4, np.exp(-ceiling) * 1.25, f"ceiling for $d_b$ = 3.7 cm: "
        f"{100 * (1 - np.exp(-ceiling)):.2f} % conversion", color="C1", fontsize=8.5)
ax.set_xlim(0, 8); ax.set_ylim(0.006, 1.05)
ax.set_xlabel(r"$\\mathcal{K}_m = (1-\\varepsilon_m)L_m K_r/u_0$")
ax.set_ylabel(r"$1-X_A$")
ax.set_title("Kobayashi's ozone bed: the three curves of the paper's figure 9")
ax.legend(fontsize=8.5, loc="lower left")
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""The three solid curves are the ones printed on the paper's figure 9, at the
bubble sizes the authors fitted. Two features are worth naming.

**The bed falls below backmix flow.** Not asymptotically — it crosses. Beyond a
finite $\mathcal{K}_m$ a bubbling bed converts less than a perfectly stirred
tank of the same catalyst, which no dispersion model can reproduce at all.

**And it flattens.** As $K_r \to \infty$ equation 49's bracket collapses to
$K_{bc}/K_r$, and $\mathcal{K}_f \to K_{bc}L_f/u_b$ — a **finite ceiling** that
does not contain the rate constant. Reaction stops being the limiting step and
bubble-to-cloud interchange takes over completely."""))

cells.append(code('''print("Where the bubbling bed becomes worse than a stirred tank, and where it stops:")
print(f"  {'d_b [cm]':>9}{'crosses backmix at K_m':>25}{'ceiling K_f':>14}{'max conversion':>17}")
rows = []
for db in (2.0, 3.0, 3.7, 4.2, 5.0, 7.0):
    h = hydrodynamics(**{**C, "db": db})
    ceil = h["Kbc"] * h["Lf"] / h["ub"]
    f = lambda K, h=h: np.exp(-Kf_closed(K, h)) - 1.0 / (1.0 + K)
    cross = brentq(f, 1e-8, 50.0) if f(1e-8) * f(50.0) < 0 else np.inf
    rows.append((db, cross, ceil, 1 - np.exp(-ceil)))
    txt = f"{cross:.3f}" if np.isfinite(cross) else "> 50"
    print(f"  {db:9.1f}{txt:>25}{ceil:14.3f}{100 * (1 - np.exp(-ceil)):16.4f} %")

Km_cross_37 = [r[1] for r in rows if r[0] == 3.7][0]
ceiling_37 = [r[3] for r in rows if r[0] == 3.7][0]
print(f"\\nAt the fitted d_b = 3.7 cm this bed cannot exceed {100 * ceiling_37:.2f} % conversion,")
print("however active the catalyst; the number is set by K_bc, L_f and u_b alone.")'''))

cells.append(md(r"""The ceiling has a closed form, and the interesting thing about it is what it
does *not* contain. Substituting equation 4 for $L_f$ and $(1-\delta)u_b = u_{br}$:

$$
\mathcal{K}_f^{\max} = \frac{K_{bc}L_f}{u_b}
= \frac{K_{bc}(d_b)\,(1-\varepsilon_m)L_m}{(1-\varepsilon_{mf})\,u_{br}(d_b)}
$$

Both $K_{bc}$ and $u_{br}$ depend on the bubble size and not on the gas
velocity, so **the ceiling is exactly independent of $u_0$**. Push more gas
through a bed and, at fixed bubble size, the best conversion it can ever reach
does not move at all."""))

cells.append(code('''db_grid = np.linspace(1.5, 8.0, 300)
fig, ax = plt.subplots(figsize=(6.6, 4.3))
for Lm, style in ((17.0, ":"), (34.0, "-"), (68.0, "--")):
    ceil = np.array([hydrodynamics(**{**C, "db": db, "Lm": Lm})["Kbc"]
                     * hydrodynamics(**{**C, "db": db, "Lm": Lm})["Lf"]
                     / hydrodynamics(**{**C, "db": db, "Lm": Lm})["ub"]
                     for db in db_grid])
    ax.plot(db_grid, 100 * (1 - np.exp(-ceil)), style, lw=1.9,
            label=f"$L_m$ = {Lm:.0f} cm" + ("  (appendix C)" if Lm == 34 else ""))
ax.plot(C["db"], 100 * ceiling_37, "o", ms=7, color="C3")
ax.annotate(f"{100 * ceiling_37:.2f} %", (C["db"], 100 * ceiling_37),
            textcoords="offset points", xytext=(8, -12), color="C3")
ax.set_xlabel(r"effective bubble diameter $d_b$ [cm]")
ax.set_ylabel("maximum attainable conversion [%]")
ax.set_title("The ceiling: what an infinitely fast catalyst would give")
ax.set_ylim(30, 101); ax.legend(fontsize=9)
fig.tight_layout(); plt.show()

# the u_0-independence is exact, not approximate
ceil_vs_u0 = [hydrodynamics(**{**C, "u0": u0})["Kbc"]
              * hydrodynamics(**{**C, "u0": u0})["Lf"]
              / hydrodynamics(**{**C, "u0": u0})["ub"] for u0 in (5.0, 13.2, 40.0, 200.0)]
print("ceiling K_f at u_0 = 5, 13.2, 40, 200 cm/s:",
      ", ".join(f"{c:.10f}" for c in ceil_vs_u0))
ceiling_u0_spread = float(np.ptp(ceil_vs_u0) / np.mean(ceil_vs_u0))
print(f"spread {ceiling_u0_spread:.2e} - independent of u_0 to rounding")
print("\\nThis is an ALGEBRAIC IDENTITY, not a check: K_bc(d_b) and u_br(d_b) contain")
print("no u_0 at all, so the spread is forced to rounding for any coefficient values.")
print("It confirms the claim in the paragraph above and nothing else. The Validation")
print("section injects a spurious u_0 dependence to show what it would take to move it.")'''))

cells.append(md(r"""This is a statement about equation 49 and about this bed's $L_m$,
$\varepsilon_m$, $\varepsilon_{mf}$ and $u_{mf}$ — not a universal chart. What it
shows is that the ceiling is governed by $d_b$ and by how deep the bed is, and by
nothing else. Doubling the bubble size from 3 to 6 cm costs more conversion than
any catalyst can give back, and the only lever left is bed depth, which enters
linearly in $\mathcal{K}_f^{\max}$ and therefore only logarithmically in the
conversion.

The caveat is the one the paper spends its conclusion on: real bubbles grow with
gas velocity, so $u_0$ does move the ceiling — through $d_b$, which the model
does not predict."""))

cells.append(code('''# The distributed model: first with the cloud and emulsion frozen (= the paper),
# then with the convection the paper leaves out.
Km_list = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
print("Effect of the convection equations 48b and 48c leave out")
print(f"  {'K_m':>6}{'1-X eq. 54':>13}{'1-X extended':>15}{'change':>10}{'  (v_e from overall balance)'}")
ext = []
for Km_ in Km_list:
    base = np.exp(-Kf_closed(Km_, hC))
    r = solve_distributed(Km_, hC, n=1600, convect=True, ve_mode="balance")
    ext.append(r["Cexit"])
    print(f"  {Km_:6.1f}{base:13.4f}{r['Cexit']:15.4f}{100 * (r['Cexit'] - base) / base:9.1f} %")

print(f"\\n  {'K_m':>6}{'1-X eq. 54':>13}{'1-X extended':>15}{'change':>10}{'  (v_e from eqs. 5 and 7)'}")
ext7 = []
for Km_ in Km_list:
    base = np.exp(-Kf_closed(Km_, hC))
    r = solve_distributed(Km_, hC, n=1600, convect=True, ve_mode="paper")
    ext7.append(r["Cexit"])
    print(f"  {Km_:6.1f}{base:13.4f}{r['Cexit']:15.4f}{100 * (r['Cexit'] - base) / base:9.1f} %")

r2 = solve_distributed(2.0, hC, n=1600, convect=True)
vb, vc, ve = r2["v"]
print(f"\\nSuperficial velocities [cm/s]: bubble {vb:.2f}, cloud+wake {vc:.2f}, "
      f"emulsion {ve:.2f}  (sum {vb + vc + ve:.2f} = u_0 {C['u0']})")
rp = solve_distributed(2.0, hC, n=400, convect=True, ve_mode="paper")
print(f"Equations 5 and 7 instead give emulsion {rp['v'][2]:.2f} cm/s, "
      f"and the three then sum to {rp['v'].sum():.2f}, "
      f"{100 * (rp['v'].sum() / C['u0'] - 1):.1f} % above u_0.")'''))

cells.append(code('''# What the profiles look like once the regions can carry gas.
r = solve_distributed(2.0, hC, n=800, convect=True)
r0 = solve_distributed(2.0, hC, n=800, convect=False)
fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.2), sharey=True)
labels = ("bubble", "cloud + wake", "emulsion")
for k in range(3):
    ax[0].plot(r0["C"][:, k], r0["x"], lw=1.8, label=labels[k])
    ax[1].plot(r["C"][:, k], r["x"], lw=1.8, label=labels[k])
ax[0].set_title("equations 48a-c as published\\n(cloud and emulsion frozen)")
ax[1].set_title("with cloud-wake carried up and\\nemulsion flowing down")
for a in ax:
    a.set_xlabel(r"$C_A/C_{Ai}$"); a.set_xlim(0, 1.02); a.legend(fontsize=9)
ax[0].set_ylabel("height above the distributor [cm]")
fig.suptitle(r"$\\mathcal{K}_m = 2$, appendix C bed", y=1.0)
fig.tight_layout(); plt.show()

print(f"Feed mixes with recirculated emulsion gas at the bottom: "
      f"C_in/C_Ai = {r['Cin']:.4f}")
print(f"Fixed point reached in {r['iters'] + 1} passes.")'''))

cells.append(md(r"""## Validation

Six checks, and they are **not of equal weight**. An earlier version of this page
ranked them the other way round, so the ranking is stated here explicitly and then
*measured* — every number below is followed by a defect-injection table showing
what moves it.

**The three that have power** are the appendix comparisons, checks 1–3. They are
the only things on this page that can see a wrong *number*: each one takes the
coefficients transcribed out of equations 10, 12 and 51 and puts them against
arithmetic the authors printed. A 22 % error in the `4.5` of equation 10 moves
check 1 by a factor of **23**.

**The three that do not** are checks 5 and 6 and the ceiling's $u_0$-independence.
Each is an **algebraic identity** of the model as implemented — guaranteed for any
parameter values — so each can catch a wrong *formula* and never a wrong *number*.
They are kept, because a wrong formula is a real error class and these are the only
checks that see it, but they are labelled below with what they cannot detect and
each carries its own injection row.

One coefficient has an external guard the appendices do not provide: `E1.2` rebuilds
the `4.5` that opens equation 10 from a Davidson-bubble flow field that never sees
this page, recovering it to 1.1e-6 as an *output*. The `5.85`, the `6.78` and the
`0.711` are pinned here by the appendix comparisons alone."""))

cells.append(code('''results = {}

print("1. Appendix C, every printed intermediate, recomputed from equations 2-51")
pairs = [("u_br", hC["ubr"], "u_br"), ("u_b", hC["ub"], "u_b"),
         ("delta", hC["delta"], "delta"), ("K_bc", hC["Kbc"], "K_bc"),
         ("K_ce", hC["Kce"], "K_ce"), ("gamma_c", hC["gc"], "gamma_c"),
         ("gamma_e", hC["ge"], "gamma_e"),
         ("K_r/K_m", Kr_of_Km(1.0, hC), "Kr_per_Km")]
print(f"   {'quantity':>10}{'printed':>12}{'recomputed':>14}{'deviation':>12}")
devC = []
for name, got, sym in pairs:
    want = printed("C", sym)
    dev = abs(got - want) / abs(want)
    devC.append(dev)
    print(f"   {name:>10}{want:12.4g}{got:14.6g}{100 * dev:11.2f} %")
results["appendixC_max_dev"] = max(devC)
print(f"   worst deviation {100 * max(devC):.2f} %, over eight numbers, nothing fitted")'''))

cells.append(code('''print("2. Appendix A: the same flow model, used for gas-solid MASS transfer")
A = dict(dp=0.028, phis=0.40, ut=41.0, D=0.065, y=0.9, Sc=2.39,
         umf=1.21, emf=0.50, gb=0.005, db=0.35)


def appendixA(coef=None):
    """Appendix A's nine printed quantities, from the same equations 1, 10, 27, 51."""
    k = KL if coef is None else coef
    Re_t = A["dp"] * RHO_G * A["ut"] / MU
    Sh_t = 2 + 0.6 * A["Sc"]**(1 / 3) * Re_t**0.5                    # eq. 1
    Kbc_A = (k["c_kbc_visc"] * A["umf"] / A["db"]
             + k["c_kbc_diff"] * (A["D"]**0.5 * G**0.25 / A["db"]**1.25))   # eq. 10
    ubr_A = k["c_ubr"] * np.sqrt(G * A["db"])                        # eq. 51
    # delta/(1-eps_f) = (u_0-u_mf) / [(1-eps_mf) u_br], since (1-delta) u_b = u_br
    brack_A = A["gb"] * Sh_t + A["y"] * A["phis"] * A["dp"]**2 / (6 * A["D"]) * Kbc_A
    slope_u0 = brack_A / ((1 - A["emf"]) * ubr_A)                    # eq. 27
    Re_per_u0 = A["dp"] * RHO_G / MU
    return {"Re_t": Re_t, "Sh_star_t": Sh_t, "K_bc": Kbc_A,
            "ub_minus_u0": ubr_A - A["umf"],
            "one_minus_ef_times_ub": (1 - A["emf"]) * ubr_A,
            "Re_per_u0": Re_per_u0, "Sh_slope_u0": slope_u0,
            "Sh_intercept": -slope_u0 * A["umf"], "Sh_slope_Re": slope_u0 / Re_per_u0}


def appendixA_dev(coef=None):
    got = appendixA(coef)
    return max(abs(v - printed("A", s)) / abs(printed("A", s)) for s, v in got.items())


got = appendixA()
print(f"   {'quantity':>22}{'printed':>12}{'recomputed':>14}{'deviation':>12}")
for sym, val in got.items():
    want = printed("A", sym)
    print(f"   {sym:>22}{want:12.4g}{val:14.6g}{100 * abs(val - want) / abs(want):11.2f} %")
results["appendixA_max_dev"] = appendixA_dev()
print(f"   (Sh)_over-all = {got['Sh_slope_Re']:.4f} Re - {-got['Sh_intercept']:.4f}"
      f"   vs the printed 0.045 Re - 0.0101")'''))

cells.append(code('''print("3. Appendix B: the same flow model again, for gas-solid HEAT transfer")
B = dict(dp=0.036, phis=0.806, ut=86.0, kg=6.25e-5, Cpg=0.24, Pr=0.69,
         umf=10.0, emf=0.50, gb=0.001, db=0.5)


def appendixB(coef=None):
    """Appendix B's chain, from equations 1, 11, 38, 51 - the same coefficients."""
    k = KL if coef is None else coef
    Re_tB = B["dp"] * RHO_G * B["ut"] / MU
    Nu_tB = 2 + 0.6 * B["Pr"]**(1 / 3) * Re_tB**0.5                  # eq. 1
    Hbc = (k["c_kbc_visc"] * B["umf"] * RHO_G * B["Cpg"] / B["db"]
           + k["c_kbc_diff"] * (B["kg"] * RHO_G * B["Cpg"])**0.5
           * G**0.25 / B["db"]**1.25)                                # eq. 11
    ubr_B = k["c_ubr"] * np.sqrt(G * B["db"])                        # eq. 51
    brack_B = B["gb"] * Nu_tB + B["phis"] * B["dp"]**2 / (6 * B["kg"]) * Hbc
    Re_per_u0B = B["dp"] * RHO_G / MU
    return dict(Re_t=Re_tB, Nu_star_t=Nu_tB, H_bc=Hbc, ubr=ubr_B,
                ub_minus_u0=ubr_B - B["umf"], Re_per_u0=Re_per_u0B,
                brack=brack_B, one_minus_ef_own=(1 - B["emf"]) * ubr_B)


def appendixB_dev(coef=None):
    b = appendixB(coef)
    pr = printed("B", "one_minus_ef_times_ub")      # 8.70, as printed
    devs = [abs(b[s] - printed("B", s)) / abs(printed("B", s))
            for s in ("Re_t", "Nu_star_t", "H_bc", "ub_minus_u0", "Re_per_u0")]
    slope = b["brack"] / pr / b["Re_per_u0"]
    inter = b["brack"] / pr * b["ub_minus_u0"]
    devs += [abs(slope - printed("B", "Nu_slope_Re")) / printed("B", "Nu_slope_Re"),
             abs(inter - printed("B", "Nu_intercept")) / printed("B", "Nu_intercept")]
    return max(devs)


b = appendixB()
one_minus_ef_own = b["one_minus_ef_own"]           # from its own stated eps_mf = 0.50
one_minus_ef_pr = printed("B", "one_minus_ef_times_ub")   # 8.70, as printed
print(f"   {'quantity':>22}{'printed':>12}{'recomputed':>14}{'deviation':>12}")
for sym in ("Re_t", "Nu_star_t", "H_bc", "ub_minus_u0", "Re_per_u0"):
    want = printed("B", sym)
    print(f"   {sym:>22}{want:12.4g}{b[sym]:14.6g}"
          f"{100 * abs(b[sym] - want) / abs(want):11.2f} %")

for tag, denom in (("as printed, 8.70", one_minus_ef_pr),
                   ("from its own eps_mf = 0.50", one_minus_ef_own)):
    slope = b["brack"] / denom / b["Re_per_u0"]
    inter = b["brack"] / denom * b["ub_minus_u0"]
    print(f"   (Nu)_apparent = {inter:.4f} + {slope:.4f} Re   [{tag}]")
print("   the printed answer is        0.0700 + 0.0510 Re")
results["appendixB_max_dev"] = appendixB_dev()
print(f"\\n   Only the printed 8.70 reproduces the printed answer, so 8.70 is what the")
print(f"   authors used. Their own stated eps_mf = 0.50 gives {one_minus_ef_own:.2f}, which")
print(f"   would move the line to {b['brack'] / one_minus_ef_own / b['Re_per_u0']:.4f} Re, "
      f"{100 * (one_minus_ef_pr / one_minus_ef_own - 1):.1f} % away.")
print(f"   8.70 corresponds to eps_mf = {1 - one_minus_ef_pr / b['ubr']:.3f}. Recorded, not repaired.")'''))

cells.append(code('''print("4. The distributed pymrm solve reproduces equation 54 when the cloud and")
print("   emulsion are frozen - two independent routes to the same number.")
print("   NOTE: convect=False is a DEGENERATE branch - v_c and v_e are identically")
print("   zero, so the recirculation loop is never entered and the fixed-point loop")
print("   never iterates. This ladder licenses the operator assembly, and it does NOT")
print("   license the convect=True numbers; those get their own ladder in check 7.")
print(f"   {'K_m':>6}{'eq. 54':>11}{'n=200':>12}{'n=400':>12}{'n=800':>12}{'n=1600':>12}{'order':>8}")
orders, finals = [], []
for Km_ in (0.5, 2.0, 8.0):
    base = np.exp(-Kf_closed(Km_, hC))
    errs = []
    for n in (200, 400, 800, 1600):
        e = abs(solve_distributed(Km_, hC, n=n, convect=False)["Cexit"] - base) / base
        errs.append(e)
    order = np.mean([np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)])
    orders.append(order); finals.append(errs[-1])
    print(f"   {Km_:6.1f}{base:11.5f}" + "".join(f"{e:12.2e}" for e in errs) + f"{order:8.2f}")
# Richardson-extrapolate the first-order upwind error away
rich = []
for Km_ in (0.5, 2.0, 8.0):
    base = np.exp(-Kf_closed(Km_, hC))
    c1 = solve_distributed(Km_, hC, n=800, convect=False)["Cexit"]
    c2 = solve_distributed(Km_, hC, n=1600, convect=False)["Cexit"]
    rich.append(abs((2 * c2 - c1) - base) / base)
results["distributed_vs_eq54_n1600"] = max(finals)
results["distributed_vs_eq54_richardson"] = max(rich)
results["upwind_order"] = float(np.mean(orders))
print(f"   first-order upwind, observed order {np.mean(orders):.3f};")
print(f"   Richardson-extrapolated agreement with equation 54: {max(rich):.1e}")'''))

cells.append(code('''print("5. The K_r -> 0 limit  [STRUCTURAL IDENTITY - see the injection table below]")
print("   Equation 57 asserts K_f = K_m; the algebra says the limit is")
print("   K_m u_0/(u_0-u_mf), because the exit stream is bubble gas only.")


def eq57_dev(h):
    """Deviation of lim_{K_r->0} K_f/K_m from u_0/(u_0-u_mf)."""
    return abs(Kf_closed(1e-8, h) / 1e-8 - C["u0"] / (C["u0"] - C["umf"])) \\
        / (C["u0"] / (C["u0"] - C["umf"]))


lim_num = Kf_closed(1e-8, hC) / 1e-8
lim_exact = C["u0"] / (C["u0"] - C["umf"])
results["eq57_limit_dev"] = eq57_dev(hC)
print(f"   numerical limit of K_f/K_m  = {lim_num:.6f}")
print(f"   u_0/(u_0-u_mf)              = {lim_exact:.6f}   deviation "
      f"{100 * results['eq57_limit_dev']:.2e} %")
print(f"   so equation 57's 'K_f = K_m' is {100 * (lim_exact - 1):.1f} % optimistic for this bed;")
print("   it is exact only as u_0/u_mf -> infinity.")
print("   WHAT THIS CANNOT DETECT: any value of the four printed coefficients. As")
print("   K_r -> 0 the bracket of eq. 49 collapses to gamma_b+gamma_c+gamma_e, which")
print("   eq. 15 DEFINES as (1-eps_mf)(1-delta)/delta, and with (1-delta)u_b = u_br")
print("   everything except u_0/(u_0-u_mf) cancels. The residual left is the O(K_r)")
print("   truncation from evaluating at K_m = 1e-8, not agreement.")

print("\\n6. Global mass balance on the extended solve")
print("   [STRUCTURAL IDENTITY - see the injection table below]")


def mass_balance_dev(h, n=2000, Km_=2.0, ve_mode="balance"):
    """Rate integrated over the bed against feed-minus-exit, on the extended solve."""
    r = solve_distributed(Km_, h, n=n, convect=True, ve_mode=ve_mode)
    Kr = Kr_of_Km(Km_, h)
    gam = np.array([h["gb"], h["gc"], h["ge"]])
    dx = h["Lf"] / r["C"].shape[0]
    consumed = float((h["delta"] * Kr * (r["C"] * gam).sum(axis=1)).sum() * dx)
    delivered = h["u0"] * (1.0 - r["Cexit"])
    return abs(consumed - delivered) / delivered, consumed, delivered


results["mass_balance_dev"], consumed, delivered = mass_balance_dev(hC)
print(f"   reactant consumed by integrating the rate over the bed : {consumed:.6f}")
print(f"   reactant lost between feed and exit  u_0 (1 - C_exit)  : {delivered:.6f}")
print(f"   closure {100 * results['mass_balance_dev']:.3f} %")
print("   WHAT THIS CANNOT DETECT: any value of the four printed coefficients. The")
print("   exchange block S has column sums -delta*gamma_k*K_r exactly - the K_bc and")
print("   K_ce terms cancel pairwise between rows - so summing the discrete equations")
print("   telescopes the convective flux, and v_b+v_c+v_e = u_0 by the closure choice.")
print("   The residual is first-order upwind truncation, and it is what it is for ANY")
print("   K_bc, K_ce, gamma or velocity that satisfies those two structural facts.")'''))

cells.append(md(r"""### Which of these six checks can see a mis-transcribed constant

The model has exactly four printed numerical coefficients — the `0.711` of
equation 51, the `4.5` and `5.85` of equation 10, and the `6.78` of equation 12.
Everything on this page is downstream of them, and the error that actually happens
when a number is read off a 1968 scan is a **transposed digit**, not a factor.

So each is perturbed by one digit and all six checks are re-evaluated. The table
below is the evidence for the ranking asserted at the top of this section."""))

cells.append(code('''INJECTIONS = [
    ("as published",                     dict(KL)),
    ("K_bc  4.5  -> 3.5   (22 % low)",   {**KL, "c_kbc_visc": 3.5}),
    ("K_bc  5.85 -> 5.58  (transposed)", {**KL, "c_kbc_diff": 5.58}),
    ("K_ce  6.78 -> 8.76  (transposed)", {**KL, "c_kce": 8.76}),
    ("u_br  0.711 -> 0.171 (transposed)", {**KL, "c_ubr": 0.171}),
]

print("Defect injection: one printed coefficient perturbed, all six checks re-run")
print(f"{'injection':<36}{'1 appC':>9}{'2 appA':>9}{'3 appB':>9}"
      f"{'5 eq57':>11}{'6 mass':>11}{'ceiling':>11}{'max X':>9}")
print(f"{'':36}{'%':>9}{'%':>9}{'%':>9}{'':>11}{'':>11}{'spread':>11}{'%':>9}")
inj_rows = []
for label, coef in INJECTIONS:
    h = hydrodynamics(**C, coef=coef)
    dC = max(abs(g - printed("C", s)) / abs(printed("C", s)) for _, g, s in
             [("u_br", h["ubr"], "u_br"), ("u_b", h["ub"], "u_b"),
              ("delta", h["delta"], "delta"), ("K_bc", h["Kbc"], "K_bc"),
              ("K_ce", h["Kce"], "K_ce"), ("gamma_c", h["gc"], "gamma_c"),
              ("gamma_e", h["ge"], "gamma_e"), ("K_r/K_m", Kr_of_Km(1.0, h), "Kr_per_Km")])
    dA, dB = appendixA_dev(coef), appendixB_dev(coef)
    mb = mass_balance_dev(h)[0]
    e57 = eq57_dev(h)
    cs = [hydrodynamics(**{**C, "u0": u}, coef=coef)["Kbc"]
          * hydrodynamics(**{**C, "u0": u}, coef=coef)["Lf"]
          / hydrodynamics(**{**C, "u0": u}, coef=coef)["ub"] for u in (5.0, 13.2, 40.0, 200.0)]
    spread = float(np.ptp(cs) / np.mean(cs))
    ceil37 = 1 - np.exp(-(h["Kbc"] * h["Lf"] / h["ub"]))
    inj_rows.append((label, dC, dA, dB, e57, mb, spread, ceil37))
    print(f"{label:<36}{100 * dC:9.2f}{100 * dA:9.2f}{100 * dB:9.2f}"
          f"{e57:11.2e}{mb:11.2e}{spread:11.2e}{100 * ceil37:9.2f}")

base = inj_rows[0]
# Rows 1-3 leave all three superficial velocities untouched, so they isolate the
# coefficient error. Row 4 (u_br) also changes u_b and therefore the whole flow
# pattern, which is why it is reported separately below.
same_flow = inj_rows[1:4]
results["injection_appendixC_worst_ratio"] = float(
    max(r[1] for r in inj_rows[1:]) / base[1])
results["injection_eq57_worst_ratio"] = float(
    max(abs(r[4]) for r in inj_rows[1:]) / abs(base[4]))
results["injection_mass_balance_worst_ratio_same_flow"] = float(
    max(abs(r[5]) for r in same_flow) / abs(base[5]))
print(f"\\nWorst-case movement, injected / as published:")
print(f"   check 1, appendix C   : x{results['injection_appendixC_worst_ratio']:8.1f}")
print(f"   check 5, eq. 57 limit : x{results['injection_eq57_worst_ratio']:8.1f}")
print(f"   check 6, mass balance : x"
      f"{results['injection_mass_balance_worst_ratio_same_flow']:8.1f}"
      "   (rows that leave the flow pattern alone)")
print("\\nThe appendix comparisons move by more than an order of magnitude on the very")
print("first row. Checks 5 and 6 stay inside their baseline band throughout, and the")
print("page's own printed verdict stays 'closure 0.001 %' on every one of them.")
print("\\nThe u_br row is the apparent exception: it moves the mass balance ~250x. It is")
print("not detection. Cutting u_br by a factor of four quintuples the emulsion downflow")
print("(v_e goes -2.3 to -12.0 cm/s), so the profile is far sharper and the FIRST-ORDER")
print("UPWIND TRUNCATION is larger. Refining n shows the residual is truncation in both:")
print(f"   {'n':>7}{'as published':>16}{'u_br 0.171':>16}")
mb_lad, mb_lad_u = [], []
h_u = hydrodynamics(**C, coef={**KL, "c_ubr": 0.171})
for n in (500, 1000, 2000, 4000):
    mb_lad.append(mass_balance_dev(hC, n=n)[0])
    mb_lad_u.append(mass_balance_dev(h_u, n=n)[0])
    print(f"   {n:7d}{mb_lad[-1]:16.3e}{mb_lad_u[-1]:16.3e}")
results["mass_balance_truncation_order"] = float(np.mean(
    [np.log2(mb_lad[i] / mb_lad[i + 1]) for i in range(len(mb_lad) - 1)]))
_ou = float(np.mean([np.log2(mb_lad_u[i] / mb_lad_u[i + 1]) for i in range(3)]))
print(f"   observed order {results['mass_balance_truncation_order']:.2f} as published and "
      f"{_ou:.2f} on the injected run - both")
print("   halve on every doubling and go to zero as the grid refines, so neither number")
print("   is an imbalance. Check 6 measures grid error, not physics.")'''))

cells.append(md(r"""### What the three identities *can* see

They are not deleted, because a wrong *formula* is a real error class and these are
the only checks on the page that catch it. Here is each one broken on purpose.

One more thing a reader should know about these three numbers: they are also the
ones the regression suite is weakest on. `scripts/check_agreement.py` ignores any
metric whose magnitude is below `ABS_FLOOR = 1e-12`, so **`ceiling_u0_spread` =
3.7e-16 is not compared by CI at all** — repo-wide, 74 of 920 published metrics sit
below that floor. A metric pinned at machine precision is not protected by the
regression suite, which is a second reason not to rank one first."""))

cells.append(code('''print("The identities, injected with the error class each actually guards against")
print(f"{'injected defect':<52}{'metric':>16}{'value':>12}")


def mb_asymmetric(scale, n=2000, Km_=2.0):
    """Mass balance with K_bc entering the two rows of S with different values."""
    h, Kr = hC, Kr_of_Km(Km_, hC)
    d, Lf, emf, n_ = h["delta"], h["Lf"], h["emf"], n
    vb = h["u0"] - h["umf"]; vc = d * h["fc"] * emf * h["ub"]; ve = h["umf"] - vc
    v = np.array([vb, vc, ve])
    x_f = np.linspace(0.0, Lf, n_ + 1); x_c = 0.5 * (x_f[:-1] + x_f[1:])
    div = construct_div((n_, 3), x_f, nu=0, axis=0)
    a_l, b_l = np.where(v > 0, 0.0, 1.0), np.where(v > 0, 1.0, 0.0)
    a_r, b_r = np.where(v < 0, 0.0, 1.0), np.where(v < 0, 1.0, 0.0)
    S = np.zeros((3, 3))
    S[0, 0] = -d * (h["gb"] * Kr + h["Kbc"]); S[0, 1] = d * h["Kbc"]
    S[1, 0] = d * h["Kbc"] * scale                    # <- the injected asymmetry
    S[1, 1] = -d * (h["Kbc"] + h["gc"] * Kr + h["Kce"]); S[1, 2] = d * h["Kce"]
    S[2, 1] = d * h["Kce"]; S[2, 2] = -d * (h["Kce"] + h["ge"] * Kr)
    Sfull = kron(identity(n_, format="csc"), csc_array(S), format="csc")
    Cin_bot, Ctop = 1.0, 1.0
    for _ in range(300):
        bc = ({"a": a_l, "b": b_l, "d": np.where(v > 0, Cin_bot, 0.0)},
              {"a": a_r, "b": b_r, "d": np.where(v < 0, Ctop, 0.0)})
        conv, conv_bc = construct_convflux_upwind((n_, 3), x_f, x_c, bc=bc, v=v, axis=0)
        Cf = spsolve((div @ conv - Sfull).tocsc(),
                     -np.asarray((div @ conv_bc).todense()).ravel()).reshape(n_, 3)
        Ct_new = (vb * Cf[-1, 0] + vc * Cf[-1, 1]) / (vb + vc)
        Cin_new = (h["u0"] + (-ve) * Cf[0, 2]) / (vb + vc)
        done = abs(Ct_new - Ctop) < 1e-13 and abs(Cin_new - Cin_bot) < 1e-13
        Ctop, Cin_bot = Ct_new, Cin_new
        if done:
            break
    up = np.maximum(v, 0.0)
    Cexit = float(up @ Cf[-1, :] / up.sum())
    gam = np.array([h["gb"], h["gc"], h["ge"]])
    consumed = float((d * Kr * (Cf * gam).sum(axis=1)).sum() * Lf / n_)
    delivered = h["u0"] * (1.0 - Cexit)
    return abs(consumed - delivered) / delivered


rows = []
rows.append(("(baseline)", "mass_balance_dev", results["mass_balance_dev"]))
rows.append(("K_bc asymmetric between the two rows of S (x1.1)",
             "mass_balance_dev", mb_asymmetric(1.1)))
rows.append(("emulsion velocity from eqs. 5 and 7 (sum != u_0)",
             "mass_balance_dev", mass_balance_dev(hC, ve_mode="paper")[0]))
rows.append(("(baseline)", "eq57_limit_dev", results["eq57_limit_dev"]))
h_ge = {**hC, "ge": hC["ge"] * 1.10}
rows.append(("gamma_e not from eq. 15's closure (x1.10)", "eq57_limit_dev", eq57_dev(h_ge)))
h_sw = {**hC, "gc": hC["ge"], "ge": hC["gc"]}
rows.append(("gamma_c and gamma_e swapped in eq. 49", "eq57_limit_dev", eq57_dev(h_sw)))
cs = [hydrodynamics(**{**C, "u0": u})["Kbc"] * (u / C["u0"])**0.1
      * hydrodynamics(**{**C, "u0": u})["Lf"] / hydrodynamics(**{**C, "u0": u})["ub"]
      for u in (5.0, 13.2, 40.0, 200.0)]
rows.append(("(baseline)", "ceiling_u0_spread", ceiling_u0_spread))
rows.append(("a spurious (u_0/13.2)^0.1 factor inside K_bc",
             "ceiling_u0_spread", float(np.ptp(cs) / np.mean(cs))))
for label, metric, val in rows:
    print(f"{label:<52}{metric:>16}{val:12.3e}")

results["identity_break_mass_balance_asymmetric_S"] = float(mb_asymmetric(1.1))
results["identity_break_eq57_gamma_e_off_closure"] = float(eq57_dev(h_ge))
results["identity_break_ceiling_u0_dependence"] = float(np.ptp(cs) / np.mean(cs))
print("\\nSo: check 6 tests that the exchange block is antisymmetric and that the three")
print("superficial velocities sum to u_0 - the bookkeeping of the assembly. Check 5")
print("tests that gamma_e obeys eq. 15 and that eq. 50's 1/(1-eps_mf) prefactor is")
print("present. The ceiling spread tests only the u_0-independence claim itself.")
print("Note the swap row: check 5 is blind even to gamma_c and gamma_e being")
print("exchanged, because only their SUM survives the K_r -> 0 limit.")'''))

cells.append(md(r"""### Refining the branch the results are actually reported on

Check 4 above refines the `convect=False` solve, which is the degenerate branch:
two of the three velocities are identically zero, the recirculation loop is never
entered, and the fixed-point iteration never runs. It establishes that the operator
assembly reproduces equation 54, and nothing more.

The number this page *reports* — the 7–15 % the extension adds, quoted in **What
pymrm adds** — comes from the `convect=True` branch, which the study above never
visits. It gets its own ladder here."""))

cells.append(code('''print("7. Grid refinement of the convect=True branch, at the conditions reported")
print(f"   {'n':>7}{'extension_max_change':>23}{'diff':>12}{'order':>8}")
ns = (200, 400, 800, 1600, 3200, 6400)
ladder = []
for n in ns:
    e = [solve_distributed(k, hC, n=n, convect=True)["Cexit"] for k in Km_list]
    ladder.append(float(max(abs(x - np.exp(-Kf_closed(k, hC))) / np.exp(-Kf_closed(k, hC))
                            for k, x in zip(Km_list, e))))
for i, (n, val) in enumerate(zip(ns, ladder)):
    if i == 0:
        print(f"   {n:7d}{val:23.6f}{'':12}{'':8}")
    else:
        d = ladder[i - 1] - val
        o = np.log2((ladder[i - 2] - ladder[i - 1]) / d) if i >= 2 else np.nan
        print(f"   {n:7d}{val:23.6f}{d:12.2e}" + (f"{o:8.2f}" if i >= 2 else f"{'':8}"))
ext_order = float(np.mean([np.log2((ladder[i - 1] - ladder[i]) / (ladder[i] - ladder[i + 1]))
                           for i in range(1, len(ladder) - 1)]))
ext_rich = float(2 * ladder[-1] - ladder[-2])          # first-order Richardson
results["extension_observed_order"] = ext_order
results["extension_max_change_n6400"] = ladder[-1]
results["extension_richardson"] = ext_rich
print(f"   observed order {ext_order:.2f} - first-order upwind, as expected")
print(f"   Richardson limit {ext_rich:.5f}; the reported n = 1600 value is "
      f"{100 * (ladder[3] / ext_rich - 1):.2f} % above it")
print("   The 7-15 % claim is a claim about the first digit and survives; the")
print("   evidence for it is this ladder, not check 4's frozen-branch ladder.")'''))

cells.append(md(r"""### Two claims the paper makes in passing, and what they are worth

Both are stated without support, and both are cheap to test once the model runs.
Neither survives the full range of the paper's own figure 9."""))

cells.append(code('''print("Claim 1 (appendix C): gamma_b anywhere in 0 to 0.01 changes the outlet")
print("concentration by less than 1 %")
print(f"   {'K_m':>6}{'gamma_b = 0':>14}{'gamma_b = 0.01':>16}{'change':>10}")
worst_gb = 0.0
for Km_ in (0.5, 1.0, 2.0, 4.0, 8.0):
    a = np.exp(-Kf_closed(Km_, hC))
    b = np.exp(-Kf_closed(Km_, hydrodynamics(**{**C, "gamma_b": 0.01})))
    worst_gb = max(worst_gb, abs(b - a) / a)
    print(f"   {Km_:6.1f}{a:14.5f}{b:16.5f}{100 * (b - a) / a:9.2f} %")
f = lambda K: (np.exp(-Kf_closed(K, hydrodynamics(**{**C, "gamma_b": 0.01})))
               - np.exp(-Kf_closed(K, hC))) / np.exp(-Kf_closed(K, hC)) + 0.01
Km_break = brentq(f, 0.05, 50.0)
results["gamma_b_claim_worst"] = worst_gb
print(f"   True up to K_m = {Km_break:.2f}. Figure 9 runs to K_m = 8, where the error is")
print(f"   {100 * worst_gb:.1f} % - five times the claim. It never changes a conclusion, but")
print("   it is not the 1 % stated.\\n")

print("Claim 2 (after eq. 53): adding the cloud-wake gas to the exit stream changes")
print("the fraction unconverted by 10 % at most")
vb = C["u0"] - C["umf"]
vc = hC["delta"] * hC["fc"] * hC["emf"] * hC["ub"]
print(f"   cloud+wake carries {100 * vc / (vb + vc):.1f} % of the upward gas flow")
print(f"   {'K_m':>6}{'bubble only':>14}{'+ cloud-wake':>15}{'change':>10}")
worst_cw = 0.0
for Km_ in (0.5, 1.0, 2.0, 4.0, 8.0):
    Kr = Kr_of_Km(Km_, hC)
    Cb = np.exp(-Kf_closed(Km_, hC))
    ge_eff = hC["Kce"] * hC["ge"] * Kr / (hC["Kce"] + hC["ge"] * Kr)   # eqs. 48b + 48c
    Cc = hC["Kbc"] * Cb / (hC["Kbc"] + hC["gc"] * Kr + ge_eff)
    mix = (vb * Cb + vc * Cc) / (vb + vc)
    worst_cw = max(worst_cw, abs(mix - Cb) / Cb)
    print(f"   {Km_:6.1f}{Cb:14.5f}{mix:15.5f}{100 * (mix - Cb) / Cb:9.1f} %")
results["cloud_wake_claim_worst"] = worst_cw
print(f"   True up to K_m near 4.9; at K_m = 8 it is {100 * worst_cw:.1f} %.")'''))

cells.append(code('''results["ceiling_u0_spread"] = ceiling_u0_spread
results["ceiling_conversion_db37"] = float(ceiling_37)
results["backmix_crossing_db37"] = float(Km_cross_37)
results["extension_max_change"] = float(max(
    abs(e - np.exp(-Kf_closed(k, hC))) / np.exp(-Kf_closed(k, hC))
    for k, e in zip(Km_list, ext)))
report_agreement("E2.1", results)'''))

cells.append(md(r"""## What pymrm adds

**Not a solver for the published model.** Equations 49, 50 and 54 are a continued
fraction and an exponential; they need numpy and nothing else, and this page says
so rather than dressing them up. What the reimplementation buys is that the three
appendices become an executable regression test — eight printed numbers in
appendix C reproduced to 0.46 %, nine in appendix A to 0.57 %, and appendix B
reproduced only when its printed 8.70 is used, which localises a slip that has
been in the literature since 1968.

**A ceiling the paper never states.** The same expression the authors fitted
contains $\mathcal{K}_f \to K_{bc}L_f/u_b$ as $K_r \to \infty$. For their own
appendix C bed that is 99.13 % conversion and no more, whatever catalyst is
loaded; it falls to 94.5 % at $d_b$ = 5 cm and 81.4 % at 7 cm; and it is
*exactly* independent of the gas velocity, because $K_{bc}$ and $u_{br}$ both
depend only on $d_b$. This follows from equation 49 alone; it is not a new model.

**And the extension the paper explicitly defers.** Equations 48b and 48c freeze
the cloud and the emulsion — they exchange gas sideways but carry none. That is
what closes the algebra, and the paper flags the cost: for $u_0/u_{mf} < 6\sim11$
the emulsion gas flows upward and "the conversion of emulsion gas will have to be
found by an extension of the present analysis". Kobayashi's runs sit at
$u_0/u_{mf}$ = 3.1 to 9.6.

The pymrm model puts the convection back: bubble gas up at $u_0-u_{mf}$,
cloud-and-wake gas carried up with the bubble at $u_b$, emulsion gas taking
whatever is left — which is *downward* here — and the two ends coupled by the
recirculation loop the paper's own figure 6 draws. With the cloud and emulsion
frozen it reproduces equation 54 to 1e-6 after Richardson extrapolation, which
licenses the operator assembly — and that ladder runs on a branch where two of the
three velocities are identically zero, so it is deliberately *not* offered as
evidence for the number below. The `convect=True` branch is refined separately, at
the conditions this number is reported at, in check 7: first order, and the
reported value sits 0.34 % above the Richardson limit. Switched on, the model
raises the unconverted fraction by 7 to 15 %, in the direction of *less*
conversion, because the cloud-wake stream is 28 % of the upward gas flow and it
bypasses the emulsion.

That number is not a correction to the paper. It is the size of the assumption,
and it is comparable to the effect of moving $d_b$ from 3.7 to 4.2 cm — which is
the one parameter the authors fitted. **The extension does not resolve anything;
it shows that the fitted bubble size is absorbing it.** That is the honest
reading, and it is the reason a page like this is worth more than a citation."""))

cells.append(md(r"""## Reuse

```python
from pymrm import construct_div, construct_convflux_upwind

h = hydrodynamics(u0=13.2, umf=2.1, emf=0.50, em=0.45, Lm=34.0,
                  db=3.7, alpha=0.47, D=0.204)
X = conversion_closed(Km=2.0, h=h)              # the published model
r = solve_distributed(Km=2.0, h=h, convect=True)  # with the convection put back
```

**What transfers.** `hydrodynamics` is the whole bubbling-bed flow model and is
reusable for any first-order process in a bubbling bed — the paper itself uses it
three times, for mass transfer, heat transfer and reaction. Change $K_r$ for
$B_d$ or $H_{bc}$ and equations 26, 38 and 49 are the same continued fraction.

**What does not.** $d_b$ is fitted, not predicted. The paper is explicit that
relating its "effective bubble size" to bubbles anyone has measured is unfinished
business, and it still largely is; do not read $d_b$ = 3.7 cm as a physical
diameter. First order only — equations 48a–c are linear in $C_A$ and the
continued fraction does not survive any other rate law, though the distributed
solve does, with `NumJac` on the last axis for a pointwise rate.

**Where else this structure appears.** A fast, catalyst-poor channel exchanging
with a slow, catalyst-rich one is the same `S7` problem as the large/small bubble
split in `F2.3`, and the same algebra as the film–ash–reaction series in `B3.1`.
If your bypassing phase carries most of the flow and almost none of the reaction,
this is your page."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.kernelspec = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata.language_info = {"name": "python"}
Path("index.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print(f"wrote index.ipynb ({len(cells)} cells)")
