#!/usr/bin/env python3
"""Generate index.ipynb for page A1.6. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Wen and Yu's minimum fluidisation velocity: where 33.7 and 0.0408 come from, and what they cost"
description: "The correlation everyone uses for u_mf, re-derived from the Ergun force balance, and tested against 21 measured minimum-fluidisation velocities it never saw."
categories: [sec:A, struct:S3, tier:T0, data:tier2, phase:gas-solid]
date: 2026-08-02
---

# Wen and Yu's minimum fluidisation velocity: where 33.7 and 0.0408 come from, and what they cost

**Catalog ID:** `A1.6` · **Structures:** `S3` · **Tier:** T0

$(N_{Re})_{mf} = \sqrt{(33.7)^2 + 0.0408\,N_{Ga}} - 33.7$ is the most-used
estimate of a minimum fluidisation velocity in the business. It is one line, it
needs no voidage and no shape factor, and those two properties are what it was
sold on.

They are also what it *costs*: the two constants exist only because two voidage
groups were replaced by numbers. This page recovers 33.7 and 0.0408 from the
Ergun force balance, works out what particle they secretly describe, and then
puts the correlation against 21 minimum fluidisation velocities measured seven
years later in another laboratory."""))

# ----------------------------------------------------------------- background
cells.append(md(r"""## Background

A packed bed lifts when the drag on it equals its buoyant weight. Set the
fixed-bed pressure gradient equal to the weight per unit volume,

$$
\frac{\Delta P}{L} \;=\; (1-\epsilon_{mf})\,(\rho_s - \rho_f)\,g ,
$$

put a two-term (Ergun-type) friction law on the left, and you have an equation
for the superficial velocity at incipient fluidisation. That is the whole
physics, and it has been known since the 1940s. The difficulty is entirely
practical: the friction law contains $\epsilon_{mf}$, the voidage *at the point
of lifting*, and $\phi_s$, the sphericity — and if you knew those you would
probably have measured $u_{mf}$ instead.

Every generalised $u_{mf}$ correlation is therefore an answer to one question:
**what do you do about $\epsilon_{mf}$ and $\phi_s$?**

- **Leva, Shirai and Wen (1956)**, and **Narsimhan (1965)** extending them, keep
  the two variables and prescribe them: Narsimhan takes $\epsilon_{mf} = 0.35$
  for spheres and makes it depend on particle diameter below 0.02 in. for
  non-spherical particles. His generalised expression is three equations.
- **Wen and Yu (1966)**, this page, do something different and much more useful.
  They observe that $\epsilon_{mf}$ and $\phi_s$ are correlated with each other,
  and that the friction law does not contain them separately — it contains two
  particular *groups* of them. Replace each group by a constant and both unknowns
  vanish together, leaving one equation with no adjustable input at all.

The communication reproduced here is three pages of *A.I.Ch.E. Journal* arguing
that the second answer beats the first. It is not the derivation — that is in a
companion paper, discussed under "The published model" below — but it prints the
result, the two approximations it rests on, the ranges they cover, and the
accuracy comparison. That is enough to take the correlation apart.

The interesting part is the *shape* of the approximation. It is not
"$\epsilon_{mf}$ is about 0.4"; it is a statement about combinations, and it
holds far better than either factor separately. What it costs is the subject of
this page, and the answer turns out to depend strongly on which end of the
Reynolds range you are at."""))

# ----------------------------------------------------------- published model
cells.append(md(r"""## The published model

### The correlation, exactly as printed

Journal page 610, equation (1):

$$
(N_{Re})_{mf} \;=\; \sqrt{(33.7)^2 + 0.0408\,N_{Ga}} \;-\; 33.7 \tag{1}
$$

with, from the notation list on journal page 612,

$$
N_{Ga} \;=\; \frac{d_p^3\,\rho_f(\rho_s-\rho_f)\,g}{\mu^2},
\qquad
N_{Re} \;=\; \frac{d_p\,\rho_f V}{\mu},
$$

$V$ the superficial fluid velocity and $(N_{Re})_{mf}$ "particle Reynolds number
at onset of fluidization". Wen and Yu call $N_{Ga}$ the **Galileo number**; it
is what most modern texts call the Archimedes number $Ar$, and this page keeps
the paper's symbol. For non-spherical particles $d_p$ "is defined as the
equivalent diameter of a spherical particle with the same volume", and the paper
adds that the geometric mean of two consecutive sieve openings is an acceptable
approximation to it.

### The two approximations the constants rest on

Journal page 611, equations (2) and (3), with the printed $\cong$:

$$
\frac{1-\epsilon_{mf}}{\phi_s^2\,\epsilon_{mf}^3} \;\cong\; 11 \tag{2}
$$

$$
\frac{1}{\phi_s\,\epsilon_{mf}^3} \;\cong\; 14 \tag{3}
$$

and the sentence that introduces them fixes their domain of validity: they
"cover the $d_p$ range from 0.002 to 1.97 in., $\epsilon_{mf}$ from 0.385 to
0.935, $\phi_s$ from 0.136 to 1.0, and with a particle diameter to column
diameter ratio ranging from 0.000807 to 0.25."

Those two groups are exactly the two that appear in the Ergun equation, which is
the point. The re-derivation is in [Results](#results).

### The accuracy claimed

Journal page 611: Equation (1) "gives an overall standard deviation of 34 % and
an average deviation of $\pm$ 25 % based on 284 points available in the
literature as shown in Figure 4", against 46 % and $\pm$ 34 % for Narsimhan's
correlations on 267 points. Table 1 breaks this down by particle class and is
transcribed below. The paper defines neither statistic.

### Two things this communication does *not* print

**It does not print the Ergun equation.** It says only, on journal page 610,
that the correlation was obtained "by employing the fixed-bed pressure drop
equation of Ergun", citing *Ergun, S., and A. A. Orning, Ind. Eng. Chem.* **41**,
1179 (1949) as reference 7. The constants of that friction law therefore have to
come from somewhere else, and on this page they come from
[`A1.1`](../A1.1-ergun-pressure-drop/), which reads them off Ergun's own 1952
paper — 150 and 1.75. **Which pair Wen and Yu actually used is not asserted here;
it is recovered, and the recovery is one of the page's checks.**

**It does not print the derivation.** Equation (1) is introduced as "the
correlation obtained by Wen and Yu (23)", and reference 23 is *Wen, C. Y., and
Y. H. Yu, Chem. Eng. Progr. Symposium Ser. No. 62, 62 (1966)* — a companion
paper, not on disk and not consulted. This communication states the result, its
two approximations, its ranges and its accuracy, and argues it against a rival.
The algebra that connects (2), (3) and the Ergun equation to (1) is *reconstructed*
here, not transcribed, and the page says so wherever it matters.

**A note on the catalogue citation.** The `A1.6` catalogue entry cites "Chem.
Eng. Prog. Symp. (1966)", which is reference 23 — the companion paper, not the
document read. The page is built from the *A.I.Ch.E. Journal* communication and
both are recorded separately in `meta.yaml`."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

Equation (1) is dimensionless and takes no parameters. Everything below it does.

| symbol | value | where it comes from |
|---|---|---|
| Ergun $k_1$ | 150 | Ergun (1952) eq. 13c, as printed — read on [`A1.1`](../A1.1-ergun-pressure-drop/) |
| Ergun $k_2$ | 1.75 | as above |
| Ergun $k_1$, refitted | 151.9 | [`A1.1`](../A1.1-ergun-pressure-drop/), from 244 markers digitised out of Ergun's own Figure 7 |
| Ergun $k_2$, refitted | 1.697 | as above |
| eq. (2) group | 11 | Wen and Yu, journal page 611 |
| eq. (3) group | 14 | Wen and Yu, journal page 611 |
| $g$ | 981 cm s$^{-2}$ | Geldart (1973), list of symbols — used only for the measured comparison |
| $\mu$ | $1.8\times10^{-4}$ poise | Geldart (1973), "for air", journal pages 289 and 290 |
| $\rho_f$ | $1.2\times10^{-3}$ g cm$^{-3}$ | **not printed by either paper.** Air at ambient; see the blind spots |

**The Ergun form used here**, in the shape that makes the two groups visible:

$$
\frac{\Delta P}{L} \;=\;
k_1\,\frac{(1-\epsilon)^2}{\phi_s^2\,\epsilon^3}\,\frac{\mu V}{d_p^2}
\;+\;
k_2\,\frac{1-\epsilon}{\phi_s\,\epsilon^3}\,\frac{\rho_f V^2}{d_p} .
$$

This is `A1.1`'s eq. (13c) with the sphericity written explicitly, i.e. with
Ergun's $D_p = 6/S_v$ replaced by $\phi_s d_p$. That substitution is a
*convention* and it is the one Wen and Yu's own notation forces, because their
$d_p$ is the volume-equivalent diameter and their $\phi_s$ is defined as the
surface-area ratio; it is not read from either paper. If it were wrong, the
recovery of 33.7 and 0.0408 below would fail, which is precisely what makes that
recovery worth doing.

**Deviation convention, used everywhere on this page without exception:**

$$
\text{deviation} \;=\; \frac{\text{model} - \text{measured}}{\text{measured}},
$$

so a negative number always means the correlation predicts a *lower* velocity
than was measured. At the 25–40 % scatter seen below, a ratio and its reciprocal
are not interchangeable.

**Assumptions carried from the sources.** Steady, isothermal, single-phase gas;
a bed deep enough that the distributor and the free surface do not matter; one
length scale per particle; incipient fluidisation treated as a single point
rather than a range, which is what a sieved narrow cut is for. Equation (1)
inherits everything the Ergun equation assumes, plus the two approximations."""))

# ------------------------------------------------------------ environment
cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code('''import sys, urllib.request
from pathlib import Path

if not any("shared" in p for p in sys.path):
    local = Path.cwd()
    for _ in range(5):
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
from scipy.optimize import brentq
from scipy.stats import spearmanr
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A1.6-wen-yu-minimum-fluidisation"
PAGE_A17 = "A1.7-geldart-classification"      # the measured u_mf come from here
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

# --- the printed constants of equation (1) ------------------------------
C1_PRINT, C2_PRINT = 33.7, 0.0408
# --- the printed voidage groups, equations (2) and (3) ------------------
G2_PRINT, G3_PRINT = 11.0, 14.0
# --- Ergun, from A1.1 ---------------------------------------------------
K1_ERGUN, K2_ERGUN = 150.0, 1.75            # as Ergun printed them
K1_REFIT, K2_REFIT = 151.93654, 1.6967056   # A1.1, refitted to 244 of his markers
# --- Geldart's stated gas properties (CGS throughout the measured test) -
G_CGS, MU_CGS, RHOF_CGS = 981.0, 1.8e-4, 1.2e-3'''))

# --------------------------------------------------------------------- data
cells.append(md(r"""## The data

Two datasets, and they do very different jobs.

### 1. Wen and Yu's Table 1 — the authors' own arithmetic, not measurement

The paper's only table. It is the standard deviation of each correlation against
literature data, by particle class. **Reproducing anything in it would be
reproduction, not validation**, and this page does not try; what it does with
Table 1 is check it against itself, which is the one thing the table can be asked
that it might fail.

**The 284 points behind it were not digitised.** They exist in the paper only as
markers in Figure 4, a log–log scatter three decades wide in
$(N_{Re})_{mf}$ and ten in $N_{Ga}$, and the routes in the Validation section are
all higher up the ladder in `docs/agent-brief.md` than a figure trace. Nothing on
this page depends on Figure 4.

### 2. Geldart's Table 1 — 21 measured minimum fluidisation velocities

The correlation needs to meet a measurement, and this paper contains none. The
gallery already holds one that fits exactly: **Geldart (1973), Table 1** —
22 narrow sieve cuts of three powders, fluidised by the author in a 5 cm column,
with $U_0$ (his symbol for $u_{mf}$) read from the pressure-drop/velocity curve.
21 of the 22 rows carry a $U_0$.

Three things make it a fair test rather than a convenient one.

- **It is genuinely held out.** Geldart's measurements are from 1973 and his own
  laboratory; Wen and Yu's 284 points are pre-1966 literature. No fitting of any
  kind connects them.
- **The independent variables are all printed**: surface/volume mean size from a
  microscope count of at least 650 particles per cut, particle density, and
  Geldart's own $\mu$ and $g$.
- **One of the three powders is spherical.** Geldart's section 4.2 describes
  Diakon as "a plastic moulding powder having spherical particles", so for its
  eight cuts $\phi_s = 1$ and Geldart's surface/volume diameter *is* Wen and Yu's
  volume-equivalent diameter. He prints no numerical sphericity — $\phi_s = 1$ is
  read off that sentence, and it is an assumption, not a measurement. Those eight
  rows are the cleanest subset *available* and are reported separately
  throughout — but "cleanest available" is not "exact", and Check 4 shows the
  headline number on them is a strong function of a $\phi_s$ nobody measured.

**The definitional mismatch on the other 13 rows is stated, not hidden.**
Geldart tabulates $d_{sv}$; Wen and Yu's $N_{Ga}$ and $N_{Re}$ are built on the
volume-equivalent $d_p$, and $d_{sv} = \phi_s d_p$. For the two cracking
catalysts $\phi_s$ is unknown and below 1, so substituting $d_{sv}$ for $d_p$
*understates* $d_p$ and therefore understates the predicted velocity. The
direction is known; the size is not, and no correction is applied.

Both catalyst densities are printed by Geldart as *approximate* ($\rho_s \simeq
1$ and $\simeq 1.5$ g cm$^{-3}$) against Diakon's exact 1.18 — another reason the
Diakon subset carries the headline."""))

cells.append(code('''wy = load_data("wen-yu-1966-table1.csv", page=PAGE)
wy_meta = load_meta("wen-yu-1966-table1.csv", page=PAGE)
gd = load_data("geldart_1973_table1.csv", page=PAGE_A17)
gd_meta = load_meta("geldart_1973_table1.csv", page=PAGE_A17)

print("Wen & Yu (1966) Table 1, journal page 611 - the authors' own statistics")
print(wy.to_string(index=False))
print("\\nprose, journal page 611: average deviation +/- "
      f"{wy_meta['notes']['prose_average_deviations']['equation_1_pct']} % for eq. (1), "
      f"+/- {wy_meta['notes']['prose_average_deviations']['narsimhan_pct']} % for Narsimhan")
print(f"\\n{cite_data(wy_meta)}")

meas = gd.dropna(subset=["U_0"]).copy()
print(f"\\n\\nGeldart (1973) Table 1: {len(gd)} rows, {len(meas)} with a measured U_0")
for p, g in meas.groupby("powder", sort=False):
    print(f"   {p:16s} n={len(g)}  d_sv {g.d_sv_um.min():3.0f}-{g.d_sv_um.max():3.0f} um   "
          f"U_0 {g.U_0.min():.2f}-{g.U_0.max():.2f} cm/s   rho_s={g.rho_s.iloc[0]:g}"
          f"{'' if g.rho_s_exact.iloc[0] else ' (approximate)'}")
print(f"\\n{cite_data(gd_meta)}")'''))

# ------------------------------------------------------ pymrm implementation
cells.append(md(r"""## PyMRM implementation

**There is no field to discretise on this page and it would be dishonest to
invent one.** Equation (1) is the positive root of a quadratic; the unapproximated
Ergun force balance is the positive root of a *different* quadratic; the
comparison between them is algebra from end to end. No `construct_grad`, no
`construct_div`, no Newton solve, and nothing here would run differently if pymrm
were absent. This page follows [`A1.1`](../A1.1-ergun-pressure-drop/) and
[`F1.4`](../F1.4-krishna-ellenberger-holdup/) in saying so in the section where a
reader would otherwise expect a solver.

What follows is therefore four functions. The only one worth reading twice is
`ergun_re_mf`, which is what equation (1) would be if nobody had approximated
anything."""))

cells.append(code('''def re_mf_wen_yu(N_Ga, C1=C1_PRINT, C2=C2_PRINT):
    """Wen & Yu (1966) eq. (1): (N_Re)_mf = sqrt(C1^2 + C2 N_Ga) - C1."""
    return np.sqrt(C1**2 + C2 * np.asarray(N_Ga, float)) - C1


def ergun_re_mf(N_Ga, phi_s, eps_mf, k1=K1_ERGUN, k2=K2_ERGUN):
    """The SAME force balance with no approximation: solve

        k2/(phi eps^3) Re^2  +  k1 (1-eps)/(phi^2 eps^3) Re  =  N_Ga

    for the positive root.  This is the reference equation (1) is measured
    against; it needs phi_s and eps_mf, which is the whole point.
    """
    a = k2 / (phi_s * eps_mf**3)
    b = k1 * (1.0 - eps_mf) / (phi_s**2 * eps_mf**3)
    return (np.sqrt(b**2 + 4.0 * a * np.asarray(N_Ga, float)) - b) / (2.0 * a)


def wen_yu_constants(k1=K1_ERGUN, k2=K2_ERGUN, g2=G2_PRINT, g3=G3_PRINT):
    """C1 and C2 of eq. (1) implied by an Ergun pair and the two voidage groups.

    Substituting eq. (2) and eq. (3) into the force balance gives
        (k2 g3) Re^2 + (k1 g2) Re = N_Ga,
    whose positive root is sqrt(C1^2 + C2 N_Ga) - C1 with
        C1 = k1 g2 / (2 k2 g3),   C2 = 1 / (k2 g3).
    """
    return k1 * g2 / (2.0 * k2 * g3), 1.0 / (k2 * g3)


def u_mf_from_re(Re, d, mu=MU_CGS, rho_f=RHOF_CGS):
    """Superficial velocity from a particle Reynolds number, CGS."""
    return Re * mu / (d * rho_f)'''))

# ------------------------------------------------------------------ results
cells.append(md(r"""## Results

### 1. Where 33.7 and 0.0408 come from

Set the Ergun pressure gradient equal to the buoyant bed weight,

$$
k_1\,\frac{(1-\epsilon_{mf})^2}{\phi_s^2\epsilon_{mf}^3}\,\frac{\mu V}{d_p^2}
+ k_2\,\frac{1-\epsilon_{mf}}{\phi_s\epsilon_{mf}^3}\,\frac{\rho_f V^2}{d_p}
= (1-\epsilon_{mf})(\rho_s-\rho_f)g ,
$$

divide through by $(1-\epsilon_{mf})$ and multiply by
$d_p^3\rho_f/\mu^2$. Every dimensional quantity collapses into $N_{Re}$ and
$N_{Ga}$, and what is left of the packing is exactly the two groups Wen and Yu
approximate:

$$
k_1\,\underbrace{\frac{1-\epsilon_{mf}}{\phi_s^2\epsilon_{mf}^3}}_{\text{eq. (2)}\;\cong\;11}
\,(N_{Re})_{mf}
\;+\;
k_2\,\underbrace{\frac{1}{\phi_s\epsilon_{mf}^3}}_{\text{eq. (3)}\;\cong\;14}
\,(N_{Re})_{mf}^2
\;=\; N_{Ga}.
$$

That is the substance of the paper in one line: **the friction law never needs
$\epsilon_{mf}$ and $\phi_s$ separately.** With the two numbers inserted it
becomes $24.5\,Re^2 + 1650\,Re = N_{Ga}$, and the positive root of that is
equation (1) with

$$
C_1 = \frac{k_1 \cdot 11}{2\,k_2\cdot 14},
\qquad
C_2 = \frac{1}{k_2 \cdot 14}.
$$

The next cell evaluates them."""))

cells.append(code('''C1_rec, C2_rec = wen_yu_constants()
d1 = 100 * (C1_rec - C1_PRINT) / C1_PRINT
d2 = 100 * (C2_rec - C2_PRINT) / C2_PRINT

print("Recovering equation (1) from the Ergun force balance")
print(f"  Ergun as printed on A1.1:  k1 = {K1_ERGUN:g}, k2 = {K2_ERGUN:g}")
print(f"  Wen & Yu eqs. (2), (3):    {G2_PRINT:g} and {G3_PRINT:g}")
print(f"  quadratic:                 {K2_ERGUN*G3_PRINT:g} Re^2 + {K1_ERGUN*G2_PRINT:g} Re = N_Ga")
print()
print(f"  C1 recovered = {C1_rec:.6f}   printed 33.7    deviation {d1:+.3f} %")
print(f"  C2 recovered = {C2_rec:.6f}   printed 0.0408  deviation {d2:+.3f} %")
print()
print("  Both printed values are the recovered ones rounded to three significant")
print(f"  figures: {C1_rec:.6f} -> {C1_rec:.3g}, and {C2_rec:.8f} -> {C2_rec:.3g}")

# what is actually identifiable from the two printed constants
print("\\nThe two identifiable groups (the constants are not independent):")
print(f"  viscous limit   Re/N_Ga     -> C2/(2 C1) = 1/(k1 g2) = {C2_PRINT/(2*C1_PRINT):.6e}"
      f"   vs 1/{K1_ERGUN*G2_PRINT:g} = {1/(K1_ERGUN*G2_PRINT):.6e}")
print(f"  turbulent limit Re^2/N_Ga   -> C2         = 1/(k2 g3) = {C2_PRINT:.6f}"
      f"   vs 1/{K2_ERGUN*G3_PRINT:g} = {1/(K2_ERGUN*G3_PRINT):.6f}")'''))

cells.append(md(r"""### 2. What particle do 11 and 14 describe?

Equations (2) and (3) are two equations in two unknowns, so *together* they pin
$\phi_s$ and $\epsilon_{mf}$ down completely. Dividing (2) by (3) gives
$\phi_s = 14(1-\epsilon_{mf})/11$, and putting that back into (3) leaves a single
quartic,

$$
(1-\epsilon_{mf})\,\epsilon_{mf}^3 \;=\; \frac{11}{14^2} .
$$

It has two roots in $(0,1)$, and one of them is outside the paper's own stated
$\phi_s$ range. The admissible one is the effective particle that equation (1)
describes — for every powder, whatever it is really made of."""))

cells.append(code('''rhs = G2_PRINT / G3_PRINT**2
f = lambda e: (1 - e) * e**3 - rhs
roots = []
for lo, hi in [(1e-9, 0.75), (0.75, 1 - 1e-12)]:
    if f(lo) * f(hi) < 0:
        e = brentq(f, lo, hi, xtol=1e-14)
        roots.append((e, G3_PRINT * (1 - e) / G2_PRINT))

PHI_RANGE = (0.136, 1.0)        # journal page 611
EPS_RANGE = (0.385, 0.935)      # journal page 611
print(f"(1 - eps) eps^3 = {G2_PRINT:g}/{G3_PRINT:g}^2 = {rhs:.8f}   ->  {len(roots)} roots in (0,1)")
for e, p in roots:
    inside = PHI_RANGE[0] <= p <= PHI_RANGE[1]
    print(f"   eps_mf = {e:.6f}   phi_s = {p:.6f}   "
          f"[{'inside' if inside else 'OUTSIDE'} the paper's own phi_s range "
          f"{PHI_RANGE[0]}-{PHI_RANGE[1]}]")
    print(f"      back-substitution: eq.(2) = {(1-e)/(p**2*e**3):.9f}, eq.(3) = {1/(p*e**3):.9f}")

EPS_EFF, PHI_EFF = [r for r in roots if PHI_RANGE[0] <= r[1] <= PHI_RANGE[1]][0]
print(f"\\nEffective particle behind 33.7 and 0.0408:  eps_mf = {EPS_EFF:.4f}, phi_s = {PHI_EFF:.4f}")
print("\\nFor comparison, the eps_mf values Wen & Yu themselves print for SPHERES")
print("(journal pages 610-611): a range 0.36 to 0.46, and reported averages")
print("0.386 (their ref. 21), 0.40 (ref. 4), 0.42 (ref. 23, their own data).")
print(f"A sphere has phi_s = 1, not {PHI_EFF:.3f}. The pair (2)+(3) is NOT a spherical packing.")

# --- and if you insist the particle IS a sphere, the two disagree ---------
EPS_FROM_2 = brentq(lambda e: (1 - e) / e**3 - G2_PRINT, 0.05, 0.95, xtol=1e-14)
EPS_FROM_3 = (1.0 / G3_PRINT) ** (1.0 / 3.0)
print(f"\\nSet phi_s = 1 and the two approximations stop agreeing with each other:")
print(f"   eq. (2) alone forces eps_mf = {EPS_FROM_2:.4f}   [(1-eps)/eps^3 = {G2_PRINT:g}]")
print(f"   eq. (3) alone forces eps_mf = {EPS_FROM_3:.4f}   [1/eps^3       = {G3_PRINT:g}]")
print(f"   they differ by {100*(EPS_FROM_3-EPS_FROM_2)/EPS_FROM_2:.1f} %, which is why (2)+(3) "
      f"together have to reach for phi_s = {PHI_EFF:.3f}.")
print(f"\\n{EPS_FROM_2:.4f} is the one that matters for small particles: eq. (2) sets the")
print(f"viscous term, so on SPHERES eq. (1) behaves like a packing at eps_mf = {EPS_FROM_2:.3f} -")
print(f"below the {EPS_RANGE[0]} lower bound of the paper's own stated eps_mf range.")'''))

cells.append(md(r"""That is the first substantive result and it is worth stating plainly.

The constants 33.7 and 0.0408 behave as though every powder were a packing of
sphericity $\phi_s \approx 0.67$ at a voidage $\epsilon_{mf} \approx 0.47$. Wen
and Yu spend the first column of their paper arguing that Narsimhan's
$\epsilon_{mf} = 0.35$ "seems to be too small" for spheres and citing 0.386, 0.40
and 0.42 instead — and then use a pair of approximations that, taken together,
describe something quite different from a sphere. There is no contradiction:
(2) and (3) are fitted across their whole population, in which most particles are
*not* spheres. But it means **the correlation should be expected to be least
accurate for spherical particles**, which is checkable, and is checked below.

**And forcing $\phi_s = 1$ shows the two approximations are not consistent with
each other for a sphere.** Equation (2) alone then demands
$\epsilon_{mf} = 0.383$ and equation (3) alone demands $\epsilon_{mf} = 0.415$ —
about 8 % apart, and only the first of them matters at small particle sizes,
because equation (2) is what multiplies the viscous term. That single number,
0.383, is the whole explanation of what happens to the fine powders in Check 3:
it sits *below* the 0.385 lower end of the range the paper itself states for its
own population, and well below the 0.40–0.42 that spherical packings actually
show.

The second root of the quartic, $\phi_s \approx 0.089$, is excluded by the
paper's own printed $\phi_s$ range, 0.136 to 1.0. It is recorded here because a
reader solving the same quartic will find it.

### 3. What the approximation costs, and where

Equation (1) against the same force balance solved without approximating
anything, for a sphere ($\phi_s=1$) at each of the voidages the paper itself
quotes."""))

cells.append(code('''N_Ga_grid = np.logspace(-2, 10, 600)
eps_quoted = [0.385, 0.386, 0.40, 0.42, 0.46]   # 0.385/0.46 the stated span, rest printed averages
re_wy_grid = re_mf_wen_yu(N_Ga_grid)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
rows = []
for eps in eps_quoted:
    ref = ergun_re_mf(N_Ga_grid, 1.0, eps)
    dev = 100 * (re_wy_grid - ref) / ref
    ax[1].semilogx(N_Ga_grid, dev, label=f"$\\\\epsilon_{{mf}}$ = {eps}")
    rows.append((eps, dev[0], dev[-1], dev.min(), dev.max()))
ax[0].loglog(N_Ga_grid, re_wy_grid, "k", lw=2, label="Wen & Yu eq. (1)")
for eps in (0.40, 0.42):
    ax[0].loglog(N_Ga_grid, ergun_re_mf(N_Ga_grid, 1.0, eps), "--",
                 label=f"Ergun, $\\\\phi_s$=1, $\\\\epsilon_{{mf}}$={eps}")
ax[0].set_xlabel("$N_{Ga}$"); ax[0].set_ylabel("$(N_{Re})_{mf}$")
ax[0].set_title("equation (1) and the unapproximated balance"); ax[0].legend(fontsize=8)
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set_xlabel("$N_{Ga}$"); ax[1].set_ylabel("deviation of eq. (1), %")
ax[1].set_title("cost of eqs. (2) and (3), spheres"); ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()

print("deviation of eq. (1) from the unapproximated Ergun balance, phi_s = 1")
print(f"{'eps_mf':>8s} {'viscous limit':>15s} {'turbulent limit':>17s}   {'over 1e-2 < N_Ga < 1e10':>26s}")
for eps, dv, dt, lo, hi in rows:
    print(f"{eps:8.3f} {dv:+14.1f} % {dt:+16.1f} %   {lo:+11.1f} % .. {hi:+.1f} %")

DEV_VISC_042 = [r[1] for r in rows if r[0] == 0.42][0]
DEV_TURB_042 = [r[3 if False else 2] for r in rows if r[0] == 0.42][0]
DEV_VISC_040 = [r[1] for r in rows if r[0] == 0.40][0]
print(f"\\nAt the eps_mf = 0.42 that Wen & Yu print for their own spherical data,")
print(f"eq. (1) is {abs(DEV_VISC_042):.1f} % LOW in the viscous limit and only "
      f"{abs(DEV_TURB_042):.1f} % {'low' if DEV_TURB_042 < 0 else 'high'} in the turbulent limit.")'''))

cells.append(md(r"""The asymmetry is the second substantive result, and it is large.

The two approximations do not degrade the correlation uniformly: they degrade
the **viscous** term far more than the inertial one. In the small-particle limit
equation (1) reduces to $(N_{Re})_{mf} = N_{Ga}/(k_1\!\times\!11) = N_{Ga}/1650$,
carrying the full error of equation (2). In the large-particle limit it reduces
to $(N_{Re})_{mf} = \sqrt{N_{Ga}/(k_2\!\times\!14)} = \sqrt{N_{Ga}/24.5}$, and
the square root halves whatever error equation (3) carries.

This matters for reading the rest of the page. Every measurement available here
is at the viscous end, which is where the correlation is at its worst — so the
deviations reported below are a lower bound on its quality, not a verdict on it.

### 4. Is any of this the Ergun constants' fault?

`A1.1` refitted Ergun's two constants to 244 markers digitised out of his own
Figure 7 and got 151.9 and 1.697 against his printed 150 and 1.75. Pushing that
uncertainty through the same derivation separates the two candidate explanations
for the deviations above."""))

cells.append(code('''C1_rf, C2_rf = wen_yu_constants(K1_REFIT, K2_REFIT)
re_rf = re_mf_wen_yu(N_Ga_grid, C1_rf, C2_rf)
re_rec = re_mf_wen_yu(N_Ga_grid, C1_rec, C2_rec)          # the UNROUNDED recovered pair
dev_ergun_only = 100 * (re_rf - re_rec) / re_rec           # refit vs recovered: no rounding folded in
dev_ergun_vs_printed = 100 * (re_rf - re_wy_grid) / re_wy_grid   # refit vs the PRINTED 33.7/0.0408

print("A1.1's refit of Ergun's own constants, pushed through the same derivation")
print(f"  k1: {K1_ERGUN:g} -> {K1_REFIT:.5g}   ({100*(K1_REFIT-K1_ERGUN)/K1_ERGUN:+.2f} %)")
print(f"  k2: {K2_ERGUN:g} -> {K2_REFIT:.5g}  ({100*(K2_REFIT-K2_ERGUN)/K2_ERGUN:+.2f} %)")
print(f"  C1: {C1_rec:.4f} -> {C1_rf:.4f}   ({100*(C1_rf-C1_rec)/C1_rec:+.2f} %)")
print(f"  C2: {C2_rec:.6f} -> {C2_rf:.6f} ({100*(C2_rf-C2_rec)/C2_rec:+.2f} %)")
print(f"\\n  effect on the PREDICTION over 1e-2 < N_Ga < 1e10: "
      f"{dev_ergun_only.min():+.2f} % .. {dev_ergun_only.max():+.2f} %")
print("  (measured against the recovered pair, not the printed 33.7/0.0408, so the")
print("   rounding of the printed constants is not folded into the comparison.)")

# --- LIKE FOR LIKE. Both costs vary across the Reynolds range, in OPPOSITE ---
# --- directions, so they must be compared at the same limit or not at all. --
ERG_VISC, ERG_TURB = float(dev_ergun_only[0]), float(dev_ergun_only[-1])
VOID_VISC, VOID_TURB = DEV_VISC_042, DEV_TURB_042
RATIO_VISC = abs(VOID_VISC) / abs(ERG_VISC)
RATIO_TURB = abs(VOID_TURB) / abs(ERG_TURB)
ERGUN_SHARE = max(abs(dev_ergun_only.min()), abs(dev_ergun_only.max()))
VOIDAGE_SHARE = abs(DEV_VISC_042)
print(f"\\n{'limit':>10s} {'eqs. (2),(3) at eps=0.42':>26s} {'Ergun constants (A1.1 refit)':>30s} "
      f"{'ratio':>8s}")
print(f"{'viscous':>10s} {VOID_VISC:+25.2f} % {ERG_VISC:+29.2f} % {RATIO_VISC:7.1f}x")
print(f"{'turbulent':>10s} {VOID_TURB:+25.2f} % {ERG_TURB:+29.2f} % {RATIO_TURB:7.2f}x")
print(f"\\nThe two costs are NOT one ratio. The voidage approximation outweighs the")
print(f"friction-law uncertainty by {RATIO_VISC:.1f}x in the VISCOUS limit and by only "
      f"{RATIO_TURB:.2f}x in the")
print(f"TURBULENT one, where they are the same size ({abs(VOID_TURB):.1f} % against "
      f"{abs(ERG_TURB):.1f} %). Quoting the worst")
print(f"voidage cost ({VOIDAGE_SHARE:.1f} %, viscous) against the worst Ergun cost "
      f"({ERGUN_SHARE:.2f} %, turbulent)")
print(f"would give {VOIDAGE_SHARE/ERGUN_SHARE:.0f}x, but that is two different limits and is not a")
print("like-for-like comparison. At the coarse end - exactly where this page tells a")
print("reader to expect a good answer - refitting Ergun's constants matters as much as")
print("the approximation does.")
print(f"\\n(The two limit-wise ratios are quoted against the RECOVERED pair. Against the")
print(f" printed 33.7/0.0408 instead they are "
      f"{abs(VOID_VISC)/abs(dev_ergun_vs_printed[0]):.0f}x viscous and "
      f"{abs(VOID_TURB)/abs(dev_ergun_vs_printed[-1]):.2f}x turbulent - the printed")
print(f" constants are rounded to three figures, and quoting against them folds Wen &")
print(f" Yu's rounding into a comparison that is supposed to be about Ergun's constants.")
print(f" Either way the answer is the same: ~20-25x in the viscous limit, ~1.2x in the")
print(f" turbulent one, and NOT one number.)")
print("\\nNote the constants move more than the prediction does. C1 and C2 are not")
print("independently identifiable from a Re: only 1/(k1 g2) and 1/(k2 g3) are, and")
print("a shift in k1 and k2 moves C1 and C2 partly in the same direction.")'''))

# --------------------------------------------------------------- validation
cells.append(md(r"""## Validation

Four checks, in the order of `docs/agent-brief.md`: an internal identity the
paper must satisfy, a second internal identity inside its own table, a
comparison against measurement, and a deliberate-break table for each. Nothing
here is a figure trace.

### Check 1 — the re-derivation, and whether it could have failed

The recovery of 33.7 and 0.0408 to a tenth of a per cent is only evidence if
plausible alternatives give something else. Each row below changes one input to
the derivation and re-reads the constants."""))

cells.append(code('''variants = [
    ("none - as built (Ergun 150/1.75, eqs 2 and 3)", K1_ERGUN, K2_ERGUN, G2_PRINT, G3_PRINT),
    ("eqs. (2) and (3) swapped",                      K1_ERGUN, K2_ERGUN, G3_PRINT, G2_PRINT),
    ("eq. (2) mis-read as 1.1",                       K1_ERGUN, K2_ERGUN, 1.1,      G3_PRINT),
    ("eq. (3) mis-read as 4",                         K1_ERGUN, K2_ERGUN, G2_PRINT, 4.0),
    ("Kozeny-Carman k1 = 180 instead of 150",         180.0,    K2_ERGUN, G2_PRINT, G3_PRINT),
    ("Eisfeld-Schnitzlein infinite bed, 154 / 1.32",  154.0,    1.32,     G2_PRINT, G3_PRINT),
    ("(not a defect) A1.1 refit, 151.94 / 1.6967",    K1_REFIT, K2_REFIT, G2_PRINT, G3_PRINT),
]
print(f"{'injected defect':47s} {'C1':>9s} {'dev %':>8s} {'C2':>9s} {'dev %':>8s}")
brk1 = []
for name, k1, k2, g2, g3 in variants:
    c1, c2 = wen_yu_constants(k1, k2, g2, g3)
    e1, e2 = 100*(c1-C1_PRINT)/C1_PRINT, 100*(c2-C2_PRINT)/C2_PRINT
    brk1.append((name, c1, e1, c2, e2))
    print(f"{name:47s} {c1:9.3f} {e1:+8.2f} {c2:9.6f} {e2:+8.2f}")
WORST_ALT = min(max(abs(e1), abs(e2)) for _, _, e1, _, e2 in brk1[1:-1])
BEST_READ = max(abs(d1), abs(d2))
REFIT_MOVE = max(abs(brk1[-1][2]), abs(brk1[-1][4]))
print(f"\\nSmallest disturbance any of the five wrong readings produces: {WORST_ALT:.2f} %,")
print(f"against {BEST_READ:.3f} % for the reading actually used - a factor "
      f"{WORST_ALT/BEST_READ:.0f}.")
print(f"\\nRESOLVING POWER. The last row is not a defect: it is A1.1's refit of Ergun's")
print(f"own constants to his own figure, and it moves the recovered constants by "
      f"{REFIT_MOVE:.1f} %.")
print(f"So this check resolves 150 from 180 or from 154; it does NOT resolve 150 from")
print(f"151.9. It identifies the FORM and the rounded pair, not the third digit.")
print("\\nWHAT THIS CHECK CANNOT DO:")
print("  * It cannot tell whether Ergun's 150 and 1.75 are RIGHT. It identifies which")
print("    friction-law constants Wen & Yu used, given that they used the Ergun form.")
print("  * It cannot separate a compensating pair of errors: k1 = 300 with eq.(2) = 5.5")
print(f"    gives C1 = {wen_yu_constants(300.0, K2_ERGUN, 5.5, G3_PRINT)[0]:.3f} and "
      f"C2 = {wen_yu_constants(300.0, K2_ERGUN, 5.5, G3_PRINT)[1]:.6f} - both exact.")
print("    Only the products k1*g2 and k2*g3 are identifiable.")'''))

cells.append(md(r"""### Check 2 — Wen and Yu's Table 1 against itself

A standard deviation over a pooled population is the root-mean-square of its
subgroup values weighted by their counts. Table 1 prints three subgroups, their
point counts, and an overall — so it can be asked whether the overall is the
pooling of the rows. Footnote * states the Narsimhan basis is 267 points, and the
three counts sum to exactly that, so for that column the question is
well posed. Footnote † says the Equation (1) overall is over all 284, so for
*that* column it is not, and the check reports the discrepancy rather than
excusing it."""))

cells.append(code('''rowsT = wy[wy.row < 4]
n = rowsT.n_points.to_numpy(float)
nar = rowsT.std_dev_narsimhan_pct.to_numpy(float)
eq1 = rowsT.std_dev_eq1_pct.to_numpy(float)
overall = wy[wy.row == 4].iloc[0]
N_NAR, N_EQ1 = 267, 284   # footnotes * and +

def pool(w, v):
    w, v = np.asarray(w, float), np.asarray(v, float)
    return float(np.sqrt((w * v**2).sum() / w.sum()))

pool_nar, pool_eq1 = pool(n, nar), pool(n, eq1)
print(f"row counts {n.astype(int)} sum to {int(n.sum())}; footnote * gives {N_NAR} for Narsimhan")
print(f"\\nNarsimhan  : rows pool to {pool_nar:.3f} %  vs printed {overall.std_dev_narsimhan_pct:g} %"
      f"   ({100*(pool_nar-overall.std_dev_narsimhan_pct)/overall.std_dev_narsimhan_pct:+.2f} %)")
print(f"Equation(1): rows pool to {pool_eq1:.3f} %  vs printed {overall.std_dev_eq1_pct:g} %"
      f"   ({100*(pool_eq1-overall.std_dev_eq1_pct)/overall.std_dev_eq1_pct:+.2f} %)")
print(f"\\nThe Narsimhan column closes to {abs(100*(pool_nar-46)/46):.1f} % on a printed integer,")
print("which tests all four of its numbers and the three point counts at once.")
print(f"\\nThe Equation (1) column does NOT close ON THE PRINTED COUNTS, and on that")
print(f"basis cannot be made to: adding the")
print(f"{N_EQ1-N_NAR} extra points would have to REDUCE a pooled rms from {pool_eq1:.2f} to "
      f"{overall.std_dev_eq1_pct:g},")
print("which is impossible for any deviation those points could carry -")
print(f"  required: ({N_EQ1} x {overall.std_dev_eq1_pct:g}^2 - {N_NAR} x {pool_eq1:.3f}^2) / "
      f"{N_EQ1-N_NAR} = "
      f"{(N_EQ1*overall.std_dev_eq1_pct**2 - N_NAR*pool_eq1**2)/(N_EQ1-N_NAR):.0f}, a NEGATIVE variance.")
print("So the Equation (1) row entries are computed on a different split of the 284")
print("than the printed counts, which the footnotes imply but do not spell out.")
print("This page reports that; it does not repair it.")
print("\\nSTATED PRECISELY, because the strong version is not proved: what cannot")
print("close is the printed row COUNTS as the basis of the Equation (1) column.")
print(f"With an unconstrained split of the {N_EQ1} among the three classes the printed "
      f"{overall.std_dev_eq1_pct:g}")
print(f"IS reachable - it needs sum(n v^2) = {N_EQ1*overall.std_dev_eq1_pct**2:.0f}, and with these three row")
print(f"values the achievable range is [{N_EQ1*min(eq1)**2:.0f}, {N_EQ1*max(eq1)**2:.0f}].")
print("The column is therefore not shown to be impossible; the stated basis is.")
print("The result is also robust to what 'standard deviation' means: read about the")
print("group means instead of about zero, pooling gives sum n (s^2 + (m-M)^2) / N,")
print("which is >= the value used here, so the contradiction only widens.")

print("\\n--- break table: does the pooling identity move? ---")
brk2 = [("none - as transcribed",            n,          nar),
        ("135.3 mis-read as 35.3",           n,          [43.4, 38.6, 35.3]),
        ("43.4 and 38.6 swapped",            n,          [38.6, 43.4, 135.3]),
        ("equal weights instead of counts",  np.ones(3), nar),
        ("the 9-point row dropped",          n[:2],      nar[:2]),
        ("count 203 mis-read as 20",         [55, 20, 9], nar)]
for name, w, v in brk2:
    print(f"  {name:35s} pooled = {pool(w, v):7.2f} %   vs printed 46")'''))

cells.append(md(r"""### Check 3 — against 21 measured minimum fluidisation velocities

This is the only check on the page that involves a measurement. Equation (1) is
evaluated on Geldart's printed sizes, densities, $\mu$ and $g$, and compared with
his measured $U_0$. Nothing is fitted, and Geldart's data postdates the
correlation by seven years."""))

cells.append(code('''d_cm = meas.d_sv_um.to_numpy() * 1e-4
drho = meas.rho_s.to_numpy() - RHOF_CGS
U_meas = meas.U_0.to_numpy()

N_Ga_meas = d_cm**3 * RHOF_CGS * drho * G_CGS / MU_CGS**2
Re_meas_pred = re_mf_wen_yu(N_Ga_meas)
u_wy = u_mf_from_re(Re_meas_pred, d_cm)
u_dr = 8e-4 * G_CGS * d_cm**2 * drho / MU_CGS       # Geldart eq. (3), Davies & Richardson

res = meas[["powder", "d_sv_um", "rho_s", "U_0", "U_MB"]].copy()
res["N_Ga"] = N_Ga_meas
res["Re_mf"] = Re_meas_pred
res["u_WenYu"] = u_wy
res["dev_WenYu_pct"] = 100 * (u_wy - U_meas) / U_meas
res["u_DavRich"] = u_dr
res["dev_DavRich_pct"] = 100 * (u_dr - U_meas) / U_meas
print(res.drop(columns=["U_MB"]).to_string(index=False, float_format=lambda x: f"{x:10.4g}"))

def summarise(sub, label):
    dw = sub.dev_WenYu_pct.to_numpy(); dd = sub.dev_DavRich_pct.to_numpy()
    return dict(label=label, n=len(sub), bias=dw.mean(), mad=np.abs(dw).mean(),
                rms=float(np.sqrt((dw**2).mean())), dr_bias=dd.mean(), dr_mad=np.abs(dd).mean())

groups = [summarise(res, "all rows"),
          summarise(res[res.powder == "Diakon"], "Diakon (spherical)"),
          summarise(res[res.powder == "Fresh catalyst"], "fresh catalyst"),
          summarise(res[res.powder == "Spent catalyst"], "spent catalyst")]
print(f"\\n{'subset':22s} {'n':>3s} {'WY bias':>9s} {'WY MAD':>8s} {'WY rms':>8s} |"
      f" {'D&R bias':>9s} {'D&R MAD':>8s}")
for g in groups:
    print(f"{g['label']:22s} {g['n']:3d} {g['bias']:+8.2f}% {g['mad']:7.2f}% {g['rms']:7.2f}% |"
          f" {g['dr_bias']:+8.2f}% {g['dr_mad']:7.2f}%")

ALL, DIA = groups[0], groups[1]
SD_ABOUT_MEAN = float(np.sqrt(ALL["rms"]**2 - ALL["bias"]**2))
print(f"\\nWen & Yu claim (journal page 611) a standard deviation of "
      f"{overall.std_dev_eq1_pct:g} % and an average")
print(f"deviation of +/- {wy_meta['notes']['prose_average_deviations']['equation_1_pct']} % "
      f"over their own 284 points. On these {ALL['n']} held-out rows the same")
print(f"statistics come out at {ALL['rms']:.1f} % and {ALL['mad']:.1f} %.")
print("\\nBUT THE PAPER DEFINES NEITHER STATISTIC, AND THE VERDICT TURNS ON THAT.")
print(f'  "standard deviation" read as rms about ZERO   : {ALL["rms"]:.1f} %  '
      f"vs their {overall.std_dev_eq1_pct:g} %  -> worse")
print(f'  "standard deviation" read as s.d. about the MEAN: {SD_ABOUT_MEAN:.1f} %  '
      f"vs their {overall.std_dev_eq1_pct:g} %  -> an exact match")
print(f"  [sqrt({ALL['rms']:.2f}^2 - {abs(ALL['bias']):.2f}^2) = {SD_ABOUT_MEAN:.2f}]")
print("Both readings are printed because nothing in the paper chooses between them.")
print("The rms-about-zero reading is the one used elsewhere on this page, and it is")
print("the less flattering of the two to this page's own conclusion; a reader who")
print("prefers the other reading should conclude that eq. (1) reproduces its authors'")
print("claimed scatter on data they never saw, with a systematic offset on top.")'''))

cells.append(code('''fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
mk = {"Diakon": ("o", "tab:blue"), "Fresh catalyst": ("s", "tab:orange"),
      "Spent catalyst": ("^", "tab:green")}
for p, g in res.groupby("powder", sort=False):
    m, c = mk[p]
    ax[0].loglog(g.U_0, g.u_WenYu, m, color=c, label=p)
    ax[1].semilogx(g.d_sv_um, g.dev_WenYu_pct, m, color=c, label=p)
lim = [0.05, 6]
ax[0].plot(lim, lim, "k-", lw=1, label="parity")
for f_ in (0.75, 1.25):
    ax[0].plot(lim, [f_*x for x in lim], "k:", lw=0.8)
ax[0].set_xlim(lim); ax[0].set_ylim(lim)
ax[0].set_xlabel("measured $U_0$, cm/s"); ax[0].set_ylabel("eq. (1) prediction, cm/s")
ax[0].set_title("Wen & Yu vs Geldart's measurements (dotted: $\\\\pm$25 %)")
ax[0].legend(fontsize=8)
ax[1].axhline(0, color="k", lw=0.8)
ax[1].axhline(DIA["bias"], color="tab:blue", ls="--", lw=1,
              label=f"Diakon mean {DIA['bias']:+.1f} %")
ax[1].axhline(DEV_VISC_042, color="k", ls="-.", lw=1,
              label=f"predicted from eqs (2),(3) at $\\\\epsilon_{{mf}}$=0.42: {DEV_VISC_042:+.1f} %")
ax[1].set_xlabel("$d_{sv}$, $\\\\mu$m"); ax[1].set_ylabel("deviation of eq. (1), %")
ax[1].set_title("deviation vs size"); ax[1].legend(fontsize=7)
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""The three powders behave differently and the split is worth reading before the
summary. The two catalysts are the rows carrying both the unknown sphericity and
the approximate density; the spent catalyst in particular has measured
velocities so low for its stated $\rho_s \simeq 1.5$ that *both* correlations
overpredict it badly, Davies and Richardson worse than Wen and Yu. That is a
property of those six rows, not a result about either correlation, and it is why
the headline below is the Diakon subset.

**The eight spherical rows are the ones worth reading in detail.** Diakon is the
only powder here that Geldart calls spherical, so it is the only one on which his
$d_{sv}$ can be identified with Wen and Yu's $d_p$ at all. Whether that
identification is *exact* depends on taking six words of prose at face value, and
Check 4 shows the headline number is a strong function of exactly that.

**Two of the five parts below are identities, not evidence, and both are labelled
as such.** This matters more than it sounds:

- Part **(b)** inverts each measured $U_0$ for the voidage the exact balance
  would need. That is a **change of variables on the same measurement**, not an
  independent quantity: the gap between $\epsilon_{mf} = 0.383$ and the inverted
  $\epsilon$ *is* the deviation re-expressed, by construction, for any data
  whatsoever. Both sides also share $\phi_s = 1$ and Ergun's $k_1 = 150$, which
  are needed to turn a constant into a voidage at all. The one thing in (b) that
  could genuinely have come out otherwise is whether the demanded voidage lands
  somewhere physically sensible — and part **(c)** shows that this is exactly the
  question the data do *not* settle.
- Part **(d)**, the viscous constant each row demands, is the same statement a
  third time: $u_{mf}$ is inversely proportional to that constant in this regime,
  so $\overline{K}/1650 - 1$ *is* the bias. Read (d) for its spread and its size
  trend, not for its centre.

**Part (c) is the contrary evidence, and it is in the page's own dataset.**
Geldart prints a bed voidage. Run the unapproximated balance at *that* voidage
and it is not a good reference at all — it is far high. The page reports this
rather than resolving it."""))

cells.append(code('''dia = res[res.powder == "Diakon"]
d_d = dia.d_sv_um.to_numpy() * 1e-4
drho_d = dia.rho_s.to_numpy() - RHOF_CGS
U_d = dia.U_0.to_numpy()
Ga_d = d_d**3 * RHOF_CGS * drho_d * G_CGS / MU_CGS**2
Re_d_meas = d_d * RHOF_CGS * U_d / MU_CGS

DEV_DIA_ROW = dia.dev_WenYu_pct.to_numpy()
print("(a0) eq. (1) row by row on the eight Diakon cuts - it is low on EVERY row,")
print("     but NOT by the same amount; the spread is a factor of ten.")
print("    " + "  ".join(f"{v:+.1f}" for v in DEV_DIA_ROW) + "  %")
print(f"    mean {DEV_DIA_ROW.mean():+.1f} %, range {DEV_DIA_ROW.min():+.1f} % to "
      f"{DEV_DIA_ROW.max():+.1f} %, rms {np.sqrt((DEV_DIA_ROW**2).mean()):.1f} %")
print("    No row is near the mean. 'low by 25 % on every row' would be wrong;")
print("    'low on every row, by 25 % on average' is what the data say. Part of")
print("    that spread is the systematic size trend quantified in (d) below.")

print("\\n(a) the exact Ergun balance, phi_s = 1, at voidages Wen & Yu themselves print")
ERGUN_AT_QUOTED = {}
for eps in (0.386, 0.40, 0.42):
    u = u_mf_from_re(ergun_re_mf(Ga_d, 1.0, eps), d_d)
    dv = 100 * (u - U_d) / U_d
    ERGUN_AT_QUOTED[eps] = float(dv.mean())
    print(f"    eps_mf = {eps:.3f}:  bias {dv.mean():+7.2f} %   MAD {np.abs(dv).mean():6.2f} %")
print(f"    eq. (1) itself      :  bias {DIA['bias']:+7.2f} %   MAD {DIA['mad']:6.2f} %")
print("  The unapproximated equation BRACKETS zero bias across the voidages Wen & Yu")
print("  print for spheres, and eq. (1) sits outside that bracket on the low side -")
print(f"  consistent with the eps_mf = {EPS_FROM_2:.3f} that eq. (2) forces at phi_s = 1.")
print(f"  But the bracket is WIDE and it is not centred: at 0.386 - Wen & Yu's own")
print(f"  ref. 21 average - the unapproximated balance is {ERGUN_AT_QUOTED[0.386]:+.2f} % biased, almost")
print(f"  as far off as eq. (1)'s {DIA['bias']:+.2f} %. Only the 0.42 end of the bracket makes")
print("  the Ergun equation look like a clean reference, and part (c) shows that the")
print("  one voidage Geldart himself reports does not support it at all.")

print("\\n(b) invert each row for the eps_mf the exact balance needs (phi_s = 1)"
      "\\n    [A REPARAMETERISATION OF (a0), NOT AN INDEPENDENT MEASUREMENT: the gap")
print("     between 0.383 and these numbers IS the deviation, by construction, and")
print("     both sides use the same phi_s = 1, the same d_sv/mu/g/drho and the same")
print("     Ergun k1 = 150. Its ONE falsifiable element is whether the demanded")
print("     voidage is physically sensible - which (c) puts in doubt.]")
eps_row = np.array([brentq(lambda e, A=A, R=R: ergun_re_mf(A, 1.0, e) - R, 0.20, 0.95, xtol=1e-12)
                    for A, R in zip(Ga_d, Re_d_meas)])
for dd, e in zip(dia.d_sv_um, eps_row):
    print(f"    d_sv = {dd:3.0f} um  ->  eps_mf = {e:.4f}")
print(f"    median {np.median(eps_row):.4f}, mean {eps_row.mean():.4f}, "
      f"range {eps_row.min():.4f} - {eps_row.max():.4f}")
print(f"    Wen & Yu print 0.36-0.46 for spheres, averages 0.386 / 0.40 / 0.42.")
print(f"    Eqs. (2)+(3) together imply {EPS_EFF:.4f} at phi_s = {PHI_EFF:.3f}.")
print(f"    Eq. (2) alone at phi_s = 1 forces {EPS_FROM_2:.4f}.")
print(f"    -> the demanded voidage and the paper's own literature values sit in the")
print(f"       same band ({np.median(eps_row):.3f} vs 0.386-0.42), BOTH above the "
      f"{EPS_FROM_2:.3f} the correlation uses.")
print(f"       That agreement is the one non-trivial content of (b) - and (c) is the")
print(f"       reason it cannot be called a closure between independent quantities.")

# --------------------------------------------------------------------------
print("\\n(c) CONTRARY EVIDENCE, from Geldart's own table: the voidage HE reports")
print("    Geldart prints eps_MB (bed voidage at minimum bubbling) and H_MB/H_0")
print("    (bed height at minimum bubbling over settled height) for every row, so")
print("    the settled voidage follows from a volume balance on the solids:")
print("        eps_0 = 1 - (1 - eps_MB) H_MB/H_0")
print("    and on the two coarsest Diakon cuts NO INFERENCE IS NEEDED AT ALL:")
print("    Geldart prints U_MB = U_0 and H_MB/H_0 = 1.000, so minimum bubbling IS")
print("    minimum fluidisation and the bed is still at its settled height. There")
print("    eps_mf = eps_MB = eps_0 exactly, as printed.")
eps_MB_d = meas.loc[dia.index, "eps_MB"].to_numpy()      # REPORTED by Geldart, method unstated
H_d = meas.loc[dia.index, "H_MB_over_H_0"].to_numpy()    # REPORTED by Geldart
eps_0_d = 1.0 - (1.0 - eps_MB_d) * H_d
u_at_MB = u_mf_from_re(ergun_re_mf(Ga_d, 1.0, eps_MB_d), d_d)
u_at_e0 = u_mf_from_re(ergun_re_mf(Ga_d, 1.0, eps_0_d), d_d)
dev_at_MB = 100 * (u_at_MB - U_d) / U_d
dev_at_e0 = 100 * (u_at_e0 - U_d) / U_d
NO_INFER = np.isclose(H_d, 1.000) & np.isclose(dia.U_MB.to_numpy(), U_d)
print(f"\\n{'d_sv':>6s} {'U_0':>6s} {'U_MB':>6s} {'H/H0':>6s} {'eps_MB':>7s} {'eps_0':>7s} "
      f"{'eps_inv':>8s} | {'eq.(1)':>8s} {'Ergun@eps_MB':>13s} {'Ergun@eps_0':>12s}")
for i in range(len(d_d)):
    tag = "  <- U_MB = U_0, H/H0 = 1.000: eps_mf = eps_MB, no inference" if NO_INFER[i] else ""
    print(f"{dia.d_sv_um.iloc[i]:6.0f} {U_d[i]:6.2f} {dia.U_MB.iloc[i]:6.2f} {H_d[i]:6.3f} "
          f"{eps_MB_d[i]:7.3f} {eps_0_d[i]:7.4f} {eps_row[i]:8.4f} | "
          f"{DEV_DIA_ROW[i]:+7.1f}% {dev_at_MB[i]:+12.1f}% {dev_at_e0[i]:+11.1f}%{tag}")
BIAS_AT_MB = float(dev_at_MB.mean())
BIAS_AT_E0 = float(dev_at_e0.mean())
print(f"\\n    mean bias, eight Diakon rows:")
print(f"      eq. (1)                                  {DIA['bias']:+8.2f} %")
print(f"      exact Ergun at Geldart's PRINTED eps_MB  {BIAS_AT_MB:+8.2f} %")
print(f"      exact Ergun at the settled eps_0         {BIAS_AT_E0:+8.2f} %")
print(f"      exact Ergun at the INVERTED eps of (b)     {0.0:+8.2f} %   (exact by construction)")
i_ni = np.flatnonzero(NO_INFER)
print(f"\\n    On the {len(i_ni)} rows where nothing is inferred, the 'unapproximated balance'")
print(f"    is {' and '.join(f'{dev_at_MB[i]:+.1f} %' for i in i_ni)} at eps_mf = "
      f"{eps_MB_d[i_ni][0]:.3f}, against eq. (1)'s "
      f"{' and '.join(f'{DEV_DIA_ROW[i]:+.1f} %' for i in i_ni)}.")
print("    THE EXACT ERGUN BALANCE IS WORSE THERE THAN THE CORRELATION IT IS BEING")
print("    USED TO JUDGE, in the opposite direction.")
print("\\n    Sphericity cannot rescue this. Expressed in d_sv, the exact Ergun u_mf is")
print("    EXACTLY phi_s-independent - the phi's cancel between d_p = d_sv/phi_s,")
print("    N_Ga, N_Re and both Ergun groups. Demonstrated, not asserted:")


def u_mf_exact_from_dsv(d_sv_cm, drho_, eps, phi):
    """Exact Ergun u_mf written in the SURFACE/VOLUME diameter Geldart prints."""
    dp = np.asarray(d_sv_cm, float) / phi
    Ga = dp**3 * RHOF_CGS * np.asarray(drho_, float) * G_CGS / MU_CGS**2
    return u_mf_from_re(ergun_re_mf(Ga, phi, eps), dp)


for phi_try in (1.0, 0.90, 0.60):
    uu = u_mf_exact_from_dsv(d_d, drho_d, eps_MB_d, phi_try)
    print(f"      phi_s = {phi_try:.2f}:  u_mf at Geldart's eps_MB = "
          + ", ".join(f"{v:.4f}" for v in uu[:3]) + ", ...  cm/s"
          + ("   <- identical" if phi_try != 1.0 else ""))
print(f"      max relative spread over those three phi_s: "
      f"{100*np.max(np.abs(u_mf_exact_from_dsv(d_d, drho_d, eps_MB_d, 0.60) / u_mf_exact_from_dsv(d_d, drho_d, eps_MB_d, 1.0) - 1)):.2e} %")
print("      The REFERENCE is phi_s-free; eq. (1) is not (Check 4). So sphericity")
print("      can move the correlation but cannot move this contradiction.")
print("\\n    SO THE PAGE DOES NOT CLAIM THAT THE BIAS BELONGS TO EQS. (2) AND (3)")
print("    RATHER THAN TO THE ERGUN EQUATION. What the data support is narrower:")
print("    the voidage eq. (2) forces at phi_s = 1 (0.383) is below the band that")
print("    reproduces Geldart's U_0 through the Ergun balance (median 0.409), and")
print("    that shortfall accounts for the bias arithmetically. But Geldart's own")
print(f"    reported voidage runs {eps_MB_d.min():.3f}-{eps_MB_d.max():.3f} (median "
      f"{np.median(eps_MB_d):.3f}) and is {eps_MB_d[NO_INFER][0]:.3f} on the two")
print(f"    rows that need no inference; at THOSE voidages the same balance")
print(f"    overpredicts his U_0 by {BIAS_AT_MB:.0f} % on average. Three readings are open, and")
print("    nothing on this page chooses between them:")
print("      (i)  eq. (2)'s voidage is too low and the Ergun balance is right at ~0.41;")
print("      (ii) Geldart's U_0 sit systematically BELOW the Ergun balance at his own")
print("           reported voidage, i.e. the reference itself is biased on this bed;")
print("      (iii) Diakon is not perfectly spherical, which moves eq. (1) but not the")
print("           reference at all (see Check 4).")
print("    CAVEAT ON eps_MB AND H_MB/H_0: Geldart's section 4.2 gives a measurement")
print("    method for U_0 and U_MB only. A1.7's data sidecar records eps_MB and")
print("    H_MB/H_0 as REPORTED, method unstated. They are not called measurements")
print("    here either - but they are the only bed voidages either paper prints.")
# --------------------------------------------------------------------------

print("\\n(d) the viscous constant each Diakon row demands"
      "   [ITS CENTRE IS THE BIAS RESTATED - read it for the spread and the trend]")
K_row = d_d**2 * drho_d * G_CGS / (MU_CGS * U_d)
K_WY = K1_ERGUN * G2_PRINT
K_DR = 1 / 8e-4
print(f"    per row: {np.round(K_row, 0)}")
print(f"    median {np.median(K_row):.0f}, spread {K_row.min():.0f}-{K_row.max():.0f} "
      f"(factor {K_row.max()/K_row.min():.2f})")
print(f"    Wen & Yu:               k1 x eq.(2) = {K_WY:g}")
print(f"    Davies & Richardson:    1 / 8e-4     = {K_DR:g}")

rho_s_, p_s = spearmanr(np.log(d_d), np.log(K_row))
print(f"\\n    residual independence: Spearman rank correlation of the demanded constant")
print(f"    with particle size is rho = {rho_s_:+.3f} (p = {p_s:.3f}) over {len(d_d)} rows.")
print("    The residuals are NOT independent of size, so no standard error is quoted")
print("    on the refitted constant below and none should be inferred from the spread.")

print("\\n(e) null baselines, on the same eight rows and the same deviation convention")
u_ref = d_d**2 * drho_d * G_CGS / (np.median(K_row) * MU_CGS)
nulls = [("constant u = mean of the measured column", np.full_like(U_d, U_d.mean())),
         ("constant u = median of the measured column", np.full_like(U_d, np.median(U_d))),
         ("one-parameter d^2 law, constant refit to these rows", u_ref),
         ("Davies & Richardson eq. (3), nothing fitted", 8e-4*G_CGS*d_d**2*drho_d/MU_CGS),
         ("Wen & Yu eq. (1), nothing fitted", u_mf_from_re(re_mf_wen_yu(Ga_d), d_d))]
for name, u in nulls:
    dv = 100 * (u - U_d) / U_d
    print(f"    {name:52s} bias {dv.mean():+7.2f} %  MAD {np.abs(dv).mean():6.2f} %")
MAD_REFIT = float(np.abs(100*(u_ref-U_d)/U_d).mean())
print(f"\\n    Reading: the two constant-velocity nulls are hopeless, so the d^2 scaling")
print(f"    is doing real work and a MAD in the teens is the floor here, not zero.")
print(f"    Eq. (1) carries {DIA['mad']:.1f} % MAD, all of it bias ({DIA['bias']:+.1f} %).")
print(f"    Refitting its single viscous constant to these very rows leaves {MAD_REFIT:.1f} % -")
print(f"    so {100-100*(DIA['mad']-MAD_REFIT)/DIA['mad']:.0f} % of eq. (1)'s error here survives ANY")
print("    one-parameter correction (part of it the size trend in (d)), and the rest")
print("    is what the approximation costs. Davies & Richardson, which is also a")
print("    one-parameter d^2 law but with a constant of 1250, is already at that floor.")
EPS_EQUIV_DR = brentq(lambda e: K1_ERGUN * (1 - e) / e**3 - K_DR, 0.10, 0.90, xtol=1e-12)
print(f"\\n    AND D&R IS NOT INDEPENDENT CORROBORATION OF THE VOIDAGE EITHER. Its hidden")
print(f"    constant {K_DR:g} is k1 (1-eps)/eps^3 at eps = {EPS_EQUIV_DR:.3f} - essentially the "
      f"{np.median(eps_row):.3f}")
print(f"    the inversion in (b) demands. It agrees because it encodes the same voidage,")
print(f"    not because two independent routes met. That is why it is presented here as a")
print(f"    null baseline already at the one-parameter floor, and not as a second witness.")'''))

cells.append(md(r"""### Check 4 — could the measured comparison have failed?

The comparison above is a chain of unit conversions on someone else's table, and
a chain of unit conversions is exactly the thing that produces a confident wrong
number. Each row below breaks one link and re-reads the Diakon bias.

**One of those links is not a unit conversion but an assumption, and it is the
one the headline rides on: $\phi_s = 1$.** Geldart prints no numerical sphericity
for Diakon; $\phi_s = 1$ is read off six words of his prose. Equation (1) is
built on the volume-equivalent $d_p$ and Geldart tabulates $d_{sv} = \phi_s d_p$,
so the prediction scales as $(d_{sv}/\phi_s)^2$ in this regime — while the
reference it is scored against is exactly $\phi_s$-free (part (c) above). All of
the sphericity risk therefore lands on the reported number, and none of it on the
reference. The break table below carries $\phi_s$ rows, and the cell solves for
the $\phi_s$ that would account for the whole bias."""))

cells.append(code('''def diakon_bias(d=None, drho_=None, mu=MU_CGS, g=G_CGS, rho_f=RHOF_CGS, U=None, fn=None):
    d = d_d if d is None else d
    dr = drho_d if drho_ is None else drho_
    U = U_d if U is None else U
    Ga = d**3 * rho_f * dr * g / mu**2
    Re = re_mf_wen_yu(Ga) if fn is None else fn(Ga)
    u = Re * mu / (d * rho_f)
    return float(np.mean(100 * (u - U) / U))

breaks = [
    ("none - as built",                              dict()),
    ("d_sv taken in mm instead of um",               dict(d=dia.d_sv_um.to_numpy()*0.1)),
    ("mu = 1.8e-5 Pa s instead of 1.8e-4 poise",     dict(mu=1.8e-5)),
    ("g = 9.81 (SI) instead of 981 cm/s2",           dict(g=9.81)),
    ("rho_f dropped from the density difference",    dict(drho_=dia.rho_s.to_numpy())),
    ("rho_f = 1.2e-4 (a decade wrong)",              dict(rho_f=1.2e-4)),
    ("measured U_MB used in place of U_0",           dict(U=dia.U_MB.to_numpy())),
    ("eq. (1) replaced by exact Ergun, eps=0.42",    dict(fn=lambda G: ergun_re_mf(G, 1.0, 0.42))),
    ("eq. (1) with C1 and C2 swapped in the root",   dict(fn=lambda G: re_mf_wen_yu(G, C2_PRINT, C1_PRINT))),
    # --- the assumption, not a unit: d_p = d_sv/phi_s -----------------------
    ("phi_s = 0.95 instead of 1 (d_p = d_sv/phi_s)",  dict(d=d_d/0.95)),
    ("phi_s = 0.90 instead of 1 (d_p = d_sv/phi_s)",  dict(d=d_d/0.90)),
]
base_bias = diakon_bias()
print(f"{'injected defect':46s} {'Diakon bias':>12s} {'moves by':>10s}")
for name, kw in breaks:
    b = diakon_bias(**kw)
    print(f"{name:46s} {b:+11.2f} % {b-base_bias:+9.2f}")

# --- how much sphericity would it take to erase the headline entirely? ------
PHI_ZERO_BIAS = brentq(lambda p: diakon_bias(d=d_d/p), 0.5, 1.0, xtol=1e-10)
print("\\nSPHERICITY SENSITIVITY - the assumption the headline rides on.")
print("Eq. (1) uses the volume-equivalent d_p; Geldart tabulates d_sv = phi_s d_p.")
print("In this regime u_mf ~ (d_sv/phi_s)^2, so the PREDICTION moves with phi_s while")
print("the reference of part (c) does not move at all.")
print(f"\\n{'assumed phi_s':>15s} {'Diakon bias':>13s}")
for p in (1.00, 0.95, 0.90, PHI_ZERO_BIAS):
    tag = "   <- the page's assumption, read off six words of prose" if p == 1.00 else (
        f"   <- the whole headline is gone" if p < 0.88 else "")
    print(f"{p:15.3f} {diakon_bias(d=d_d/p):+12.2f} %{tag}")
print(f"\\nphi_s = {PHI_ZERO_BIAS:.3f} would account for 100 % of the Diakon bias. A perfectly")
print(f"ordinary {0.95:.2f} for sieved suspension-polymerised beads removes "
      f"{100*(1-diakon_bias(d=d_d/0.95)/base_bias):.0f} % of it.")
print("So the measured verdict is a JOINT test of the voidage approximation AND of")
print("Diakon's sphericity, and the page cannot separate them: Geldart prints no")
print("numerical phi_s, only the phrase 'having spherical particles'.")
print("Under the plausible readings the headline is:")
for p, word in ((1.00, "spherical exactly"), (0.95, "very nearly spherical"),
                (0.90, "roughly spherical")):
    print(f"   phi_s = {p:.2f} ({word:22s}) -> eq. (1) low by {abs(diakon_bias(d=d_d/p)):.1f} % on average")

print("\\nWHAT THIS CHECK CANNOT DO, measured rather than asserted:")
print(f"  * rho_f. Dropping it from the density difference moves the bias by "
      f"{diakon_bias(drho_=dia.rho_s.to_numpy())-base_bias:+.2f} points,")
print("    and rho_f is not printed by either paper. The page makes no claim about it.")
print("    (It does NOT cancel from N_Ga and N_Re separately - only from the answer, in")
print(f"     the viscous limit, where u_mf = d^2 (rho_s-rho_f) g / (1650 mu).)")
print("  * The turbulent constant. See the resolving-power cell below.")
print(f"  * phi_s. It is not measured anywhere, and anything from {PHI_ZERO_BIAS:.3f} to 1 is")
print("    consistent with what Geldart writes. The headline is quoted at phi_s = 1")
print("    with that sensitivity stated, not with an error bar it does not have.")'''))

cells.append(md(r"""### What the measured test cannot resolve

Before quoting the agreement, the question `docs/agent-brief.md` insists on: is
the effect being tested larger than what the data can see? For half of equation
(1), the answer is no."""))

cells.append(code('''visc_only = N_Ga_meas / (2 * C1_PRINT / C2_PRINT)
inert_share = 100 * np.abs(Re_meas_pred - visc_only) / Re_meas_pred
print(f"(N_Re)_mf over Geldart's 21 rows: {Re_meas_pred.min():.2e} to {Re_meas_pred.max():.3f}")
print(f"Wen & Yu's Figure 4 spans (N_Re)_mf 0.001 to 4000 (journal page 611).")
print(f"\\nContribution of the inertial term of eq. (1) on these rows: "
      f"{inert_share.min():.4f} % to {inert_share.max():.2f} %")
print(f"Largest row: d_sv = {meas.d_sv_um.to_numpy()[np.argmax(inert_share)]:.0f} um.")
print("\\nSo this dataset tests 1/(k1 x eq.2) = 1/1650 and NOTHING ELSE. Setting the")
print("turbulent constant to any value at all changes the Diakon bias by:")
for c2 in (C2_PRINT/4, C2_PRINT*4):
    c1_ = C1_PRINT * (c2 / C2_PRINT)      # hold C2/(2C1) fixed, i.e. hold the viscous limit
    b = float(np.mean(100*(u_mf_from_re(re_mf_wen_yu(Ga_d, c1_, c2), d_d) - U_d)/U_d))
    print(f"   C2 x {c2/C2_PRINT:.2f} at fixed viscous limit -> bias {b:+.3f} % "
          f"({b - DIA['bias']:+.3f} points)")
print("\\nThat is the page's largest blind spot and it is structural: no group A powder")
print("can probe the inertial term. Section 3 says the approximation is at its BEST")
print("there, and nothing on this page tests that statement against a measurement.")

print("\\n--- and are these rows inside the ranges eqs. (2) and (3) are stated to cover? ---")
D_MIN_IN, D_MAX_IN = 0.002, 1.97          # journal page 611
DDC_MIN, DDC_MAX = 0.000807, 0.25         # journal page 611
COLUMN_CM = 5.0                           # Geldart section 4.2
d_in = meas.d_sv_um.to_numpy() / 25400.0
ddc = meas.d_sv_um.to_numpy() * 1e-4 / COLUMN_CM
print(f"stated d_p range      : {D_MIN_IN} to {D_MAX_IN} in. "
      f"= {D_MIN_IN*25400:.1f} um to {D_MAX_IN*25.4:.0f} mm")
print(f"Geldart's rows        : {meas.d_sv_um.min():.0f} to {meas.d_sv_um.max():.0f} um "
      f"= {d_in.min():.5f} to {d_in.max():.5f} in.")
print(f"stated d_p/D range    : {DDC_MIN} to {DDC_MAX}")
print(f"Geldart's rows (5 cm) : {ddc.min():.6f} to {ddc.max():.5f}")
out_d = meas[d_in < D_MIN_IN]
out_r = meas[ddc < DDC_MIN]
print(f"\\n{len(out_d)} of the {len(meas)} rows fall below the stated d_p bound, "
      f"{len(out_r)} below the stated d_p/D bound:")
for _, r in meas.assign(d_in=d_in, ddc=ddc).iterrows():
    flags = ("d_p" if r.d_in < D_MIN_IN else "   ") + " " + ("d_p/D" if r.ddc < DDC_MIN else "")
    if flags.strip():
        print(f"   {r.powder:16s} d_sv = {r.d_sv_um:3.0f} um   below: {flags}")
print(f"\\nAll of them are catalyst rows; every Diakon cut is inside both ranges")
print(f"({meas[meas.powder=='Diakon'].d_sv_um.min():.0f}-"
      f"{meas[meas.powder=='Diakon'].d_sv_um.max():.0f} um). Those bounds are the extremes of the")
print("population eqs. (2) and (3) were fitted over, not a proof of failure outside")
print("them - but the fresh-catalyst rows are the worst on the page, and they are")
print("also the ones at and below the fitted edge.")'''))

cells.append(md(r"""### Blind spots — claims this page does **not** make

Each of these is a thing a reader might reasonably assume the page has shown. It
has not.

1. **That equation (1) is 25 % low in general — or 25 % low on any particular
   row.** It is low on *every* one of eight spherical cuts of one plastic powder
   at $(N_{Re})_{mf} < 1$, by 25 % *on average*, and the per-row spread is a
   factor of ten (printed at the top of Check 3). No row is near the mean, and
   section 3 identifies that corner as the worst of the correlation's range. The
   paper's own claim over 284 points spanning $(N_{Re})_{mf}$ from 0.001 to 4000
   is $\pm$ 25 % *average* deviation, and this page neither confirms nor
   contradicts it there. Three of the 21 rows — all catalyst, none Diakon — are
   also finer than the smallest particle equations (2) and (3) were fitted over,
   which the last cell of Check 3 locates and which is worth knowing before
   reading the fresh-catalyst column.
2. **That Ergun's constants are right.** Check 1 identifies which friction-law
   constants Wen and Yu used; it says nothing about whether they are correct, and
   only the products $k_1\times 11$ and $k_2\times 14$ are identifiable at all.
   A compensating error in $k_1$ and in equation (2) is invisible to every check
   here.
3. **That the shape-factor claim holds.** Equation (1)'s selling point is that it
   needs no $\phi_s$. Testing that needs non-spherical powders with a *known*
   $\phi_s$, and neither paper on this page provides one — Geldart's two
   catalysts have no printed sphericity, which is exactly why substituting
   $d_{sv}$ for $d_p$ on those 13 rows is recorded as a definitional caveat and
   why the headline is the Diakon subset.
4. **That $\phi_s = 1$ for Diakon, or that the headline is insensitive to it.**
   This is the same gap as (3) pointed at the *reported number* rather than at
   the correlation's claim, and it is the strongest single dependence on the
   page. $\phi_s = 1$ rests on six words of Geldart's prose and no printed value.
   Because equation (1) is built on $d_p = d_{sv}/\phi_s$ and the reference of
   Check 3(c) is exactly $\phi_s$-free, the headline moves with $\phi_s^{-2}$ and
   the reference does not move at all: Check 4 prints the bias at
   $\phi_s = 0.95$ and $0.90$ and solves for the $\phi_s$ that would account for
   *all* of it. A perfectly ordinary departure from perfect sphericity removes a
   large fraction of the headline, and a 13–14 % departure removes the whole of
   it. The measured verdict is a joint test of the voidage approximation and of
   Diakon's sphericity, and nothing here separates the two.
5. **That the unapproximated Ergun balance is the right reference, or that the
   bias belongs to equations (2) and (3) rather than to it.** Check 3(c) is the
   page's own contrary evidence: at the bed voidage Geldart himself reports, the
   *unapproximated* balance overpredicts his $U_0$ badly — worse than equation (1)
   underpredicts it, and on two rows the voidage needs no inference at all
   because he prints $U_{MB} = U_0$ and $H_{MB}/H_0 = 1.000$. The inverted
   voidage of Check 3(b) explains the bias arithmetically; Geldart's reported
   voidage does not; and nothing on this page settles which is right.
6. **That $\rho_f$ matters.** It is not printed by either paper. The break table
   measures its influence and it is small; the page therefore asserts nothing
   about it rather than claiming insensitivity as a result.
7. **That the size trend in the Diakon residuals is understood.** The viscous
   constant each cut demands rises with size (Spearman $\rho = +0.833$) —
   *rises*, not monotonically; the printed sequence has two reversals. The page
   reports the rank correlation and its p-value and stops there. Sphericity,
   cohesion at the fine end, wall effects at
   $d_p/D = 318\,\mu\text{m}/5\,\text{cm}$ and the voidage of a sieved cut all
   vary along that axis and nothing here separates them.
8. **That anything in Wen and Yu's Figure 4 has been checked.** It was not
   digitised. The 284 points are quoted only through the statistics the authors
   computed from them.
9. **That the derivation is the authors'.** It is reconstructed. The companion
   paper that contains their own derivation was not consulted, and the evidence
   that the reconstruction is theirs is the recovery of 33.7 and 0.0408 — strong,
   but circumstantial."""))

cells.append(code('''metrics = dict(
    C1_recovered=C1_rec, C1_deviation_pct=d1,
    C2_recovered=C2_rec, C2_deviation_pct=d2,
    implied_eps_mf=EPS_EFF, implied_phi_s=PHI_EFF,
    dev_vs_exact_ergun_viscous_eps042_pct=DEV_VISC_042,
    dev_vs_exact_ergun_turbulent_eps042_pct=DEV_TURB_042,
    dev_vs_exact_ergun_viscous_eps040_pct=DEV_VISC_040,
    ergun_constant_share_pct=ERGUN_SHARE,
    voidage_approximation_share_pct=VOIDAGE_SHARE,
    # like-for-like, one limit at a time (the two above are DIFFERENT limits)
    ergun_constant_share_viscous_pct=ERG_VISC,
    ergun_constant_share_turbulent_pct=ERG_TURB,
    voidage_over_ergun_ratio_viscous=RATIO_VISC,
    voidage_over_ergun_ratio_turbulent=RATIO_TURB,
    table1_narsimhan_pooled=pool_nar,
    table1_narsimhan_printed=float(overall.std_dev_narsimhan_pct),
    table1_eq1_pooled=pool_eq1,
    table1_eq1_printed=float(overall.std_dev_eq1_pct),
    measured_rows=float(ALL["n"]),
    wenyu_bias_all_pct=ALL["bias"], wenyu_mad_all_pct=ALL["mad"], wenyu_rms_all_pct=ALL["rms"],
    wenyu_bias_diakon_pct=DIA["bias"], wenyu_mad_diakon_pct=DIA["mad"],
    davies_richardson_bias_diakon_pct=DIA["dr_bias"],
    davies_richardson_mad_diakon_pct=DIA["dr_mad"],
    refit_one_parameter_mad_diakon_pct=MAD_REFIT,
    inverted_eps_mf_median_diakon=float(np.median(eps_row)),
    inverted_viscous_constant_median_diakon=float(np.median(K_row)),
    residual_size_rank_correlation=float(rho_s_),
    residual_size_rank_pvalue=float(p_s),
    max_inertial_share_on_measured_rows_pct=float(inert_share.max()),
    eps_mf_forced_by_eq2_at_sphere=EPS_FROM_2,
    eps_mf_forced_by_eq3_at_sphere=EPS_FROM_3,
    exact_ergun_bias_diakon_eps042_pct=float(
        np.mean(100*(u_mf_from_re(ergun_re_mf(Ga_d, 1.0, 0.42), d_d) - U_d)/U_d)),
    rows_below_stated_dp_range=float(len(out_d)),
    # --- the contrary evidence in Geldart's own table (Check 3c) ---------
    exact_ergun_bias_diakon_at_reported_eps_MB_pct=BIAS_AT_MB,
    exact_ergun_bias_diakon_at_settled_eps_0_pct=BIAS_AT_E0,
    geldart_reported_eps_MB_no_inference_rows=float(int(NO_INFER.sum())),
    geldart_reported_eps_MB_no_inference_value=float(eps_MB_d[NO_INFER][0]),
    # --- per-row spread of the headline, and its sphericity dependence ----
    wenyu_dev_diakon_min_pct=float(DEV_DIA_ROW.min()),
    wenyu_dev_diakon_max_pct=float(DEV_DIA_ROW.max()),
    wenyu_rms_diakon_pct=float(np.sqrt((DEV_DIA_ROW**2).mean())),
    diakon_bias_at_phi_s_095_pct=diakon_bias(d=d_d/0.95),
    diakon_bias_at_phi_s_090_pct=diakon_bias(d=d_d/0.90),
    phi_s_that_erases_diakon_bias=float(PHI_ZERO_BIAS),
    # --- the other reading of "standard deviation" (Check 3) --------------
    wenyu_sd_about_mean_all_pct=SD_ABOUT_MEAN,
)
_ = report_agreement("A1.6", metrics)'''))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**To the correlation itself, nothing, and this page does not pretend otherwise.**
Equation (1) is a closed-form root and no operator, grid or solver appears
anywhere above. Sections 1 to 5 would run with pymrm uninstalled. What the
reimplementation adds is four things, none of which is a solver.

**The derivation, made explicit and made falsifiable.** The communication states
that equation (1) came from the Ergun equation and prints the two approximations,
but the algebra between them is in a companion paper nobody here has. Doing it
recovers 33.7 and 0.0408 to a tenth of a per cent, and the break table shows
every alternative reading of its inputs — a swapped pair, a Kozeny–Carman
constant, the Eisfeld–Schnitzlein infinite-bed constants — moving the answer by
tens of per cent. That fixes, from the outside, which friction-law constants the
correlation encodes, which the paper does not say. The same table states the
limit: it cannot separate 150 from `A1.1`'s refitted 151.9.

**The effective particle.** Equations (2) and (3) are usually quoted as two
independent conveniences. They are not independent: together they determine
$\phi_s$ and $\epsilon_{mf}$ uniquely, and the pair they determine is not a
sphere. Force $\phi_s = 1$ and they contradict each other by 8 % in
$\epsilon_{mf}$; the branch that survives into the viscous term is 0.383, which
is below the paper's own stated $\epsilon_{mf}$ floor. Anyone using equation (1)
on a fine spherical powder is using a correlation built around a packing looser
than any of the ones its authors cite, and this page is where that is computed
rather than asserted.

**The error budget, split — and split at one limit at a time.** The correlation's
deviation from the equation it came from is far larger than the uncertainty in
Ergun's own constants *in the viscous limit*, measured with `A1.1`'s refit of
those constants to Ergun's own figure — more than twenty times larger. In the
turbulent limit the two are the same size, within about 20 %. Quoting a single
ratio hides that, because the two costs move in opposite directions across the
range: the
voidage approximation is worst where the friction-law uncertainty is least and
best where it is largest. The asymmetry itself is not visible from the
correlation's form and is what a reader needs in order to set a design margin.
Small particles: expect a systematic underprediction, and the approximation is
the dominant term. Large particles: expect a better answer, but no better than
Ergun's own constants are known.

**A measurement the correlation never saw — and a contradiction inside it.** The
paper contains no data, only statistics computed from a figure. Geldart's Table 1
supplies 21 minimum fluidisation velocities from another laboratory seven years
later. Inverting the unapproximated balance row by row on the eight spherical
cuts asks what voidage would reproduce them; the answer sits in the same band as
the literature values Wen and Yu themselves quote for spheres and above the one
their correlation uses, and that shortfall accounts for the bias. **That
inversion is a reparameterisation of the bias, not an independent
corroboration** — it shares $\phi_s = 1$, the same sizes and gas properties and
Ergun's own $k_1$ with the thing it is compared to, and the gap *is* the
deviation by construction. Its one falsifiable element is whether the demanded
voidage is physically sensible, and there the same table supplies contrary
evidence the page reports rather than resolves: Geldart's own reported bed
voidage — needing no inference at all on the two coarsest cuts, where he prints
$U_{MB} = U_0$ and $H_{MB}/H_0 = 1.000$ — is *higher*, and at that voidage the
unapproximated balance overpredicts his velocities by more than equation (1)
underpredicts them. What this page adds is that both numbers are on the same
page, computed from the same table, with neither hidden."""))

# ------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**Use equation (1) when** you need $u_{mf}$ from size, density and gas
properties alone, and you do not have $\epsilon_{mf}$ or $\phi_s$. That is most
of the time, and it is why the correlation won. `re_mf_wen_yu` is two lines.

**Use `ergun_re_mf` instead when** you do have $\epsilon_{mf}$ *and* you know how
it was obtained — **and read Check 3(c) before you trust the result.** The
unapproximated balance is more faithful physics, but it is exquisitely sensitive
to the voidage you feed it, and on this page's own dataset the two available
voidages disagree by an amount that flips the verdict. At the voidage the eight
spherical cuts *demand*, it is unbiased by construction. At the voidage Geldart
*reports* for those same cuts — including two rows where no inference is needed —
it overpredicts his measured velocities by tens of per cent, worse in absolute
terms than equation (1) underpredicts them. A settled-bed height is exactly the
kind of number that produces the second case, not the first. So: `ergun_re_mf`
lets you use a measured voidage, but a measured voidage is not automatically the
right one, and $\epsilon_{mf}$ at the point of lifting is not the settled-bed
voidage unless the bed has not expanded. If you take one thing from this page
into a design, take the *sensitivity*: in the viscous-limit column of the
section-3 table equation (1) is fixed and only $\epsilon_{mf}$ moves, and between
0.385 and 0.46 the unapproximated balance's $u_{mf}$ changes by nearly a factor
of two.

**Expect the error to be size-dependent, not uniform, and do not expect a
constant offset.** Below $(N_{Re})_{mf}\approx 1$ the correlation carries the
full error of equation (2); above $(N_{Re})_{mf}\approx 100$ the square root has
halved the error of equation (3). Even within the eight-row spherical subset here
the per-row deviation spans a factor of ten and rises systematically with size.
Design margins should not be symmetric and should not be a single percentage.

**And if your particles are only nearly spherical, that matters more than it
looks.** Equation (1) wants the volume-equivalent $d_p$; most tabulated sizes are
surface/volume or sieve diameters, and $d_{sv} = \phi_s d_p$. In the viscous
regime the prediction goes as $(d_{sv}/\phi_s)^2$, so a 10 % shortfall in
sphericity moves the reported bias on these eight rows by 17 points — larger than
almost anything else on the page. Check 4 measures exactly that, and the same
table gives the $\phi_s$ at which the whole bias disappears.

**And note where $u_{mf}$ goes next.** It is an input, not an output, for most of
section E: the excess gas $u_0 - u_{mf}$ drives the bubble phase in the two-phase
theory, so a 25 % error in $u_{mf}$ at $u_0 = 3u_{mf}$ is a 12.5 % error in the
bubbling gas flow and a larger one close to incipient conditions. Nothing here
propagates that through a bed model; it is arithmetic, stated so a reader knows
which way it points.

**Related pages.** [`A1.1`](../A1.1-ergun-pressure-drop/) is the friction law
this correlation is a root of, and supplies both the printed and the refitted
constants used above. [`A1.7`](../A1.7-geldart-classification/) supplies the
measured dataset and uses a *different* $u_{mf}$ expression — Davies and
Richardson's, via Geldart's eq. (3) — inside its group A/B boundary; the two are
compared row by row in Check 3. `A1.5` (Richardson–Zaki) takes the same bed past
incipient fluidisation, and `E2.1` and `E1.2` consume $u_{mf}$ as a parameter.

**Cite the sources, not this page:** Wen, C. Y. and Yu, Y. H., *A generalized
method for predicting the minimum fluidization velocity*, A.I.Ch.E. Journal
**12**(3) 610–612 (1966),
[doi:10.1002/aic.690120343](https://doi.org/10.1002/aic.690120343); the
derivation is in Wen, C. Y. and Yu, Y. H., *Chem. Eng. Progr. Symposium Series*
No. 62, **62** (1966), which was not consulted. Geldart, D., *Types of gas
fluidization*, Powder Technology **7**(5) 285–292 (1973),
[doi:10.1016/0032-5910(73)80037-3](https://doi.org/10.1016/0032-5910(73)80037-3).
Ergun, S., *Fluid flow through packed columns*, Chemical Engineering Progress
**48**(2) 89–94 (1952)."""))

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
