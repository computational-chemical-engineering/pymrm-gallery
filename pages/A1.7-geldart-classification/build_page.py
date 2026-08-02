#!/usr/bin/env python3
"""Generate index.ipynb for page A1.7 (the Geldart classification).

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "Geldart's powder groups: the boundaries, and what his own table says about them"
description: "Geldart's 1973 classification is not a differential equation - it is two inequalities and a map. This page recomputes both boundary constants from the expressions he printed, tests the A/B boundary row by row against the 22 fractions he measured, measures what the minimum-fluidization correlation underneath it is actually worth against those same measurements, and shows that the B/D boundary is the Davidson cloud-existence condition written in other symbols."
categories: [sec:A, struct:S3, tier:T0, data:tier2, phase:gas-solid]
date: 2026-08-02
---

# Geldart's powder groups: the boundaries, and what his own table says about them

**Catalog ID:** `A1.7` · **Structures:** `S3` (the one pymrm solve on the page) · **Tier:** T0 · **Data tier:** 2 (a printed table of the author's own measurements)

Ask which fluidised-bed model applies to a powder and the first answer anyone
gives is a letter. **A** expands smoothly before it bubbles; **B** bubbles the
moment it fluidises; **C** will not fluidise at all; **D** spouts. The letters
come from one seven-page paper, and the thing the paper actually contributes is
not the letters but **two inequalities**:

$$(\rho_s-\rho_f)\,d' \le 225 \quad\text{(group A)},\qquad
  (\rho_s-\rho_f)\,(d')^2 \ge 10^6 \quad\text{(group D)},$$

with $\rho$ in g cm$^{-3}$ and $d'$ the mean particle size in µm.

Those two numbers are quantitative predictions, so they can be checked. This
page does four things with them and does not pretend to do a fifth:

1. **Recomputes 225 and $10^6$** from the expressions Geldart derives them from,
   using his own $g$, $\mu$, $K_{MB}$, $\varepsilon_0$ and $d_B$. Neither comes
   back exact.
2. **Tests the A/B boundary against the 22 size fractions Geldart measured**
   himself, 21 of which carry a measured $U_0$. The boundary uses only size and
   density; the verdict it is tested against uses only the two measured velocity
   columns. The two routes share no input column — one residual link, through
   the fit that supplies the 100 inside the 225, is named and measured in the
   validation section. Almost every fraction in the table is group A, so a
   whole-table score is nearly free: the page prints a **null baseline** beside
   it, names the handful of rows on which the criterion can possibly beat that
   baseline, and quotes its score on those rows as the result.
3. **Measures the minimum-fluidization correlation that the boundary is built
   on** — equation (3), which is Davies and Richardson's and was not fitted to
   these data — against the measured $U_0$ column. This is the weakest link in
   the boundary and the paper never quantifies it.
4. **Shows that the B/D boundary is a cloud-existence condition.** Geldart's
   equation (7) is, term for term, the statement that a 25 cm bubble rises more
   slowly than the interstitial gas — which is exactly the threshold the
   [`E1.2`](../E1.2-davidson-bubble/) page derives for the Davidson bubble to
   have no cloud at all. The two differ by the coefficient $1/\sqrt2$ against
   $0.711$, which is less than the rounding Geldart applied to his own answer.

**What this page does not do:** it does not reproduce the C/A boundary. Geldart
prints no expression for it — on his Figure 3 it is a shaded band drawn by hand
to separate the powders other authors called difficult from the ones they called
group A — and nothing here is digitised from any figure."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

Before 1973 the fluidisation literature was in the habit of generalising from
one powder. Geldart's opening complaint is exactly that: conclusions drawn from
cracking catalyst were being applied to powders of quite different size and
density, and this "is responsible for some of the apparent contradictions and
differences of opinion which appear in published papers."

His fix was to sort powders into four groups by two properties only — the
density difference $\rho_s-\rho_f$ and the mean size $d_{sv}$ — and then to be
strict about what may be carried across a group boundary. The summary is blunt
about the limit of its own usefulness:

> Generalizations concerning powders within a group can be made with reasonable
> confidence but conclusions drawn from observations made on a powder in one
> group should not in general be used to predict the behaviour of a powder in
> another group.

The four groups, in his words (section 3, journal pages 286–287):

| group | behaviour |
|---|---|
| **A** | beds "expand considerably before bubbling commences"; collapse slowly, typically 0.3–0.6 cm/s; gross powder circulation even with few bubbles; bubbles split and recoalesce; a maximum bubble size exists |
| **B** | "bubbles start to form ... at or only slightly above minimum fluidization velocity"; bed expansion small, collapse rapid; bubble size grows linearly with height and with $U-U_0$; no evidence of a maximum bubble size |
| **C** | "powders which are in any way cohesive"; lift as a plug or channel badly, because interparticle forces exceed what the fluid can exert; fixed by stirrers, vibration, fumed silica, or making the walls conducting |
| **D** | large and/or very dense; "all but the largest bubbles rise more slowly than the interstitial fluidizing gas, so that gas flows into the base of the bubble and out of the top"; can be made to spout |

Group D's definition is the one that matters later on this page: it is a
statement about a **bubble rising more slowly than the gas around it**, which is
the same object the [`E1.2`](../E1.2-davidson-bubble/) page solves for.

The classification is a `T0` case because it is upstream of everything else in
section E. Which of the gallery's fluidised-bed pages applies to a given powder
is, in practice, decided by its Geldart group: `E2.1`'s bubbling-bed model
presumes bubbles and a cloud, and `E1.2`'s cloud has a threshold below which it
does not exist. What this page adds to that chain is the arithmetic that says
which side of the threshold a powder is on."""))

# ------------------------------------------------------------ published model
cells.append(md(r"""## The published model

### Provenance

Everything below was read from **Geldart, D., "Types of gas fluidization",
*Powder Technology* **7**(5) 285–292 (1973)**, doi
[10.1016/0032-5910(73)80037-3](https://doi.org/10.1016/0032-5910(73)80037-3),
on 600 dpi renders of journal pages 285–292. The paper itself was read, not a
reprint of it.

The text layer of this Elsevier scan is not usable — it renders the paper's own
title as "Types of Gas Fhidization" — so **every equation, constant and table
value here was read visually off the rendered page image**, and nothing was
repaired by inference. Two places where the print is degraded are recorded as
judgements rather than silently fixed: the exponent in equation (7) is a
half-power, and the heading of Table 1's fifth column is $U_{MB}$ (the printed
subscript is damaged and scans as "MP"). Both are settled by the surrounding
text, and neither is a digit.

**One result on this page originates elsewhere.** Equation (3), the minimum
fluidization velocity, carries Geldart's reference 11: *L. Davies and J. F.
Richardson, Trans. Inst. Chem. Engrs., 44 (1966) T293*. That paper is **not on
disk and was not consulted**; equation (3) is used here exactly as Geldart
prints it, and every statement about "the correlation" is a statement about the
form printed on journal page 289. The same applies to equation (9)
(Verloop and Heertjes, ref. 10) and equation (10) (Oltrogge, ref. 36).

### The equations, as printed

The nomenclature is Geldart's own list of symbols, journal page 291: $g = 981$
cm s$^{-2}$; $d_{sv}$ the surface/volume diameter in **cm**; $d'$ the particle
size in **µm** ($=$ cm $\times 10^4$); $d_B$ the frontal diameter of a bubble in
cm; $U_0$ the superficial velocity at minimum fluidization, cm s$^{-1}$;
$U_{MB}$ the superficial velocity at minimum bubbling, cm s$^{-1}$;
$\varepsilon_0$ the bed voidage at minimum fluidization; $\mu$ the gas viscosity
in g cm$^{-1}$ s$^{-1}$; $\rho_f$, $\rho_s$ gas and particle density, g cm$^{-3}$.

**The bubble point** (journal page 288, from Figure 1):

$$U_{MB} = K_{MB}\,d_{sv} \tag{1}$$

$$U_{MB} = 100\,d_{sv} \tag{2}$$

$K_{MB}$ "has a value of 100 when $U_{MB}$ is in cm s$^{-1}$ and $d_{sv}$ in
cm". Geldart notes it is "curious that the bubble point can be correlated by
means of an equation which involves a term ($K_{MB}$) having units of
frequency", considers Hiby's 7–25 s$^{-1}$ bed oscillation as an explanation and
rejects it because 100 is much larger.

**Minimum fluidization** (journal page 289, his reference 11):

$$U_0 = \frac{8\times10^{-4}\,g\,d_{sv}^2\,(\rho_s-\rho_f)}{\mu} \tag{3}$$

**The A/B criterion** (journal page 289). A powder is group A when the bed
bubbles later than it fluidises:

$$\frac{U_{MB}}{U_0} \ge 1 \tag{4}$$

Substituting (1) and (3) into (4):

$$\frac{8\times10^{-4}\,g\,d_{sv}\,(\rho_s-\rho_f)}{K_{MB}\,\mu} \le 1 \tag{5}$$

and then, "for air at room temperature and pressure $K_{MB} = 100$ and
$\mu = 1.8\times10^{-4}$ poise", replacing $d_{sv}$ (cm) by $d'$ (µm):

$$(\rho_s-\rho_f)\,d' \le 225 \tag{6}$$

This is the line **XY** on his Figure 3.

**The B/D criterion** (journal pages 289–290). The density/size combinations for
which bubbles smaller than $d_B$ "would rise more slowly than the interstitial
gas velocity":

$$\left(g\,\frac{d_B}{2}\right)^{\tfrac12} \le
  \frac{8\times10^{-4}(\rho_s-\rho_f)\,g\,d_{sv}^2}{\mu\,\varepsilon_0} \tag{7}$$

with, for group D, "$d_B = 25$ cm. The choice is not critical since in eqn. (7)
we are considering $\sqrt{d_B}$. For large-particle systems $\varepsilon_0
\approx 0.4$ and for air $\mu = 1.8\times10^{-4}$ g cm$^{-1}$ s$^{-1}$":

$$(\rho_s-\rho_f)(d')^2 \ge 10^6 \tag{8}$$

Geldart calls this criterion "tentative" and says the use of (3) on the
right-hand side of (7) "is not strictly justified for these large particles
since the flow regime is transitional, not laminar."

**Two competing criteria**, for comparison (journal page 291). Verloop and
Heertjes (ref. 10), and Oltrogge (ref. 36), predict immediate bubbling when

$$\frac{(g\,d^3)^{\frac12}(\rho_s-\rho_f)}{\mu} > 5000 \tag{9}$$

$$\frac{(g\,d^3)^{\frac12}(\rho_s-\rho_f)}{\mu} > 400 \tag{10}$$

Geldart states that (9) "lies much too far to the right and is not shown on
Fig. 3", while (10) is drawn there as line O–O and "agrees with published
results about as well as eqn. (5)". Section [Validation](#validation) checks both
statements against the transcribed equations."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

Everything on this page runs in the paper's own CGS units, because every printed
constant — $8\times10^{-4}$, $K_{MB} = 100$, 225, $10^6$ — is dimensional and
carries them.

| symbol | value | where it comes from |
|---|---|---|
| $g$ | 981 cm s$^{-2}$ | list of symbols, journal page 291 |
| $\mu$ | $1.8\times10^{-4}$ g cm$^{-1}$ s$^{-1}$ | stated twice, journal pages 289 and 290, "for air" |
| $\rho_f$ | $1.2\times10^{-3}$ g cm$^{-3}$ | **not printed in this paper.** Air at ambient. See below. |
| $K_{MB}$ | 100 s$^{-1}$ | equation (2), journal page 288 |
| $\varepsilon_0$ | 0.4 | journal page 290, "for large-particle systems" |
| $d_B$ | 25 cm | journal page 290, "bubble sizes greater than 25 cm have rarely been reported" |

**The one number not printed in the paper is $\rho_f$**, and it does not matter.
Geldart writes every criterion in terms of $\rho_s-\rho_f$ and never says what he
takes $\rho_f$ to be; for air at ambient it is about $1.2\times10^{-3}$ g
cm$^{-3}$, which is 0.1 % of the smallest particle density in his own table.
The validation section measures how much of anything depends on it: the answer
is nothing, and that is recorded as a **blind spot** rather than as agreement.

**Assumptions carried by the criteria, all Geldart's own:**

- The A/B boundary inherits everything in equations (2) and (3). Equation (2) is
  a straight line fitted through the points of his Figure 1 with **no density
  term at all** — the bubble point is asserted to depend on size alone, for
  powders spanning 1 to 1.5 g cm$^{-3}$. Equation (3) is a laminar-regime
  expression whose single constant $8\times10^{-4}$ hides a voidage.
- The boundary is for **air at ambient**. Section 5.1 says so explicitly and
  predicts that raising the pressure moves line XY to the right, so that some
  group B powders become group A. That prediction is not tested here; no data
  for it exists in the paper.
- $\varepsilon_0 = 0.4$ and $d_B = 25$ cm enter only the B/D boundary, and the
  paper flags both as arbitrary-but-reasonable.

**Deviation convention, fixed once for the whole page.** Every percentage is a
signed relative deviation of the recomputed value from the printed one,

$$\delta = 100\,\frac{x_{\text{recomputed}} - x_{\text{printed}}}{x_{\text{printed}}}\ \%,$$

so a positive $\delta$ means Geldart printed a number smaller than his own
inputs give."""))

# ---------------------------------------------------------------- environment
cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

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
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
from pymrm import construct_grad, construct_div, construct_coefficient_matrix
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A1.7-geldart-classification"
pd.set_option("display.width", 130)
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

# ---- Geldart's own constants, CGS.  Every one is quoted in the table above.
P = dict(
    g       = 981.0,      # cm/s2      list of symbols, p. 291
    mu      = 1.8e-4,     # g/(cm s)   p. 289 and p. 290, "for air"
    K_MB    = 100.0,      # 1/s        eq. (2), p. 288
    C3      = 8.0e-4,     # -          the constant in eq. (3), p. 289
    eps_0   = 0.4,        # -          p. 290, "for large-particle systems"
    d_B     = 25.0,       # cm         p. 290
    rho_f   = 1.2e-3,     # g/cm3      NOT printed in the paper; air at ambient
)
PRINTED = dict(eq6=225.0, eq8=1.0e6, eq9=5000.0, eq10=400.0)
print({k: v for k, v in P.items()})'''))

# -------------------------------------------------------------------- data
cells.append(md(r"""## The data

**Table 1, journal page 288: "Experimental results on Group A powders".** 22
narrow size fractions of three powders — Diakon (a plastic moulding powder with
spherical particles), fresh cracking catalyst and spent cracking catalyst —
sieved from wide distributions, each sized by a microscope count of at least 650
particles, and fluidised as 200 g charges in a 5 cm glass column with a filter
paper distributor, giving beds about 20 cm deep.

The two columns this page's central test rests on are **measurements**, and that
is what lifts the page above the tier-6 pages elsewhere in the gallery. Section
4.2 states the method for both: $U_0$ from the pressure-drop / gas-velocity
curve, $U_{MB}$ by raising the air velocity until the first recognisable bubble
(about 0.5 cm) broke the surface, then lowering it until only one or two bubbles
remained, the two averaged over several repeats. Comparing a criterion with
those columns is validation against measurement, not reproduction of an author's
arithmetic.

**The other two columns are reported, not stated to be measured.** Section 4.2
says nothing about how $\varepsilon_{MB}$ and $H_{MB}/H_0$ were obtained.
$\varepsilon_{MB}$ in particular is not a directly readable quantity; the
natural route to it is the bed height, the 200 g charge and $\rho_s$, and
$\rho_s$ is printed as approximate for both catalysts. The page therefore calls
those two columns *reported* throughout, uses them only for the transcription
check in section 4, and claims nothing about their provenance.

Three limits on what they can test, and they are severe enough to state before
the numbers rather than after:

- **Every powder in the table is group A or borderline A/B.** The table's own
  title says so, the densities span 1 to 1.5 g cm$^{-3}$, and the sizes 25 to
  318 µm. Nothing here touches the C/A boundary or the B/D boundary. The A/B
  test below is the only boundary the data can reach.
- **Two of the three densities are printed as approximate** — "$\rho_s \simeq 1$"
  and "$\rho_s \simeq 1.5$" against Diakon's "$\rho_s = 1.18$". 14 of the 22 rows
  therefore carry a density uncertainty of unstated size, and any per-row
  quantity that divides by $\rho_s-\rho_f$ inherits it. Diakon is the clean
  series and is reported separately wherever that matters.
- **The $U_{MB}$ column is not independent of equation (2).** Figure 1 is these
  very points with the line $U_{MB} = 100\,d_{sv}$ drawn through them. Comparing
  equation (2) with the $U_{MB}$ column measures the quality of a fit, and this
  page labels it as such and never counts it as a test. The $U_0$ column carries
  no such circularity: equation (3) is Davies and Richardson's."""))

cells.append(code(r'''t1 = load_data("geldart_1973_table1.csv", page=PAGE)
meta = load_meta("geldart_1973_table1.csv", page=PAGE)
print(cite_data(meta))

t1["drho"] = t1["rho_s"] - P["rho_f"]                       # rho_s - rho_f
t1["prod6"] = t1["drho"] * t1["d_sv_um"]                    # LHS of eq. (6)
t1["ratio_meas"] = t1["U_MB"] / t1["U_0"]                   # LHS of eq. (4)
has_U0 = t1["U_0"].notna()

print(f"\n{len(t1)} rows, {has_U0.sum()} with a measured U_0 "
      f"({(~has_U0).sum()} printed as an em dash)")
print(f"d_sv  {t1.d_sv_um.min():.0f} to {t1.d_sv_um.max():.0f} um;  "
      f"rho_s {t1.rho_s.min()} to {t1.rho_s.max()} g/cm3;  "
      f"exact rho_s on {int(t1.rho_s_exact.sum())} of {len(t1)} rows")
display(t1[["powder", "nominal_range", "d_sv_um", "U_0", "U_MB",
            "eps_MB", "H_MB_over_H_0", "prod6", "ratio_meas"]].round(3))'''))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

**Be clear about what pymrm is doing here, because it is not much.** Geldart's
classification is two inequalities. It has no differential operator in it, it
does not need a solver, and section [Results](#results) evaluates it in closed
form. Building a PDE to produce a number that is three lines of arithmetic would
be padding, and this page does not do it.

The one solve on the page is there for a different reason. Equation (7) — the
B/D criterion — compares a **bubble rise velocity** with the **interstitial gas
velocity**. That comparison is the exact condition under which a rising bubble
has no cloud: gas that leaves the roof of the bubble does not come back, and the
Davidson recirculation cell does not close. The [`E1.2`](../E1.2-davidson-bubble/)
page solves that percolation problem in pymrm and derives

$$\left(\frac{R_c}{R}\right)^3 = \frac{u_{br}+2u_f}{u_{br}-u_f},$$

which has a **pole**, not a zero, at $u_{br} = u_f$. So Geldart's group D and
E1.2's no-cloud region are candidates for being the same set, and the way to
find out is to put a bubble in a Geldart-classified powder and look.

The solver below is E1.2's, condensed: the $\ell = 1$ radial mode of Laplace's
equation solved twice, once for the gas potential (Dirichlet on the bubble
surface: uniform pressure) and once for the solids potential (zero normal
gradient: no solids cross the surface), with the cloud read off as the radius
where the total radial velocity in the bubble frame changes sign. **E1.2 owns
that derivation and its validation; nothing about the operator is claimed as new
here.** It is used because it produces the cloud radius from a computed flow
field rather than from the formula being tested, and because the picture of a
25 cm bubble in a group B powder next to the same bubble in a group D powder is
the thing the classification implies and neither paper draws.

**The dependency that creates, stated plainly:** every cloud number on this page
— the computed $R_c/R$, the closed-form value it is checked against, and the
0.711 coefficient in the $d_b^{*}$ contour — belongs to `E1.2`. **If E1.2's
operator or its closed form changes, this page's cloud numbers move with it.**
Nothing in the Geldart reproduction proper (sections 1 to 4, the boundary
constants and the Table 1 classification) depends on E1.2 at all."""))

cells.append(code(r'''def _solve_mode(x_f, bc_bubble, A_far, nu=2):
    """Radial factor h(r) of the l = 1 harmonic potential  Phi = h(r) cos(theta).

    (1/r^nu) d/dr( r^nu dh/dr ) - nu h / r^2 = 0,  nu = 2 for a spherical bubble.
    This is E1.2's solve_mode with its sabotage hooks removed; see that page for
    the derivation, the grid convergence study and the operator break tests.
    """
    n = len(x_f) - 1
    x_c = 0.5 * (x_f[:-1] + x_f[1:])
    r_out = x_f[-1]

    # bc is (left, right) with a*dh/dn + b*h = d, n the OUTWARD normal:
    #   left  face r = R      : n = -r_hat, so dh/dn = -dh/dr
    #   right face r = r_out  : n = +r_hat, so dh/dn = +dh/dr
    # Right condition  dh/dr + (nu/r) h = (nu+1) A_far  is satisfied EXACTLY by
    # both radial modes h = A r and h = B r^-nu, so a finite r_out adds no
    # truncation error - only the grid does.
    bc = (bc_bubble, {"a": 1.0, "b": nu / r_out, "d": (nu + 1.0) * A_far})

    Grad, grad_bc = construct_grad((n, 1), x_f, x_c, bc)
    Div  = construct_div((n, 1), x_f, nu=nu)          # nu: 2 spherical
    Sink = construct_coefficient_matrix(-nu / x_c.reshape(-1, 1) ** 2)
    A = (Div @ Grad + Sink).tocsc()                   # constant operator
    h = spsolve(A, -(Div @ grad_bc)).reshape(-1, 1)
    return x_c, h.ravel(), (Grad @ h + grad_bc).ravel()


def davidson_cloud(d_b, u_f, u_br, n=800, r_out_mult=30.0):
    """Cloud radius around one spherical bubble, from the computed flow field.

    d_b  bubble diameter (cm);  u_f = U_0/eps_0 interstitial gas velocity (cm/s);
    u_br bubble rise velocity relative to the emulsion (cm/s).
    Returns R_c/R from the field, and nan when no dividing streamline exists.
    """
    R = 0.5 * d_b
    x_f = R * (r_out_mult ** np.linspace(0.0, 1.0, n + 1))     # geometric faces

    # GAS relative to the solids: u_rel = -grad(phi).  Uniform pressure over the
    # bubble surface -> phi = 0 there.   physical eq at r = R:  f = 0
    x_c, f, df = _solve_mode(x_f, {"a": 0.0, "b": 1.0, "d": 0.0}, -u_f)
    # SOLIDS: v = +grad(Phi_p).  No solids cross the surface.
    #   physical eq at r = R:  dg/dr = 0   (outward normal is -r_hat here)
    _,   g, dg = _solve_mode(x_f, {"a": 1.0, "b": 0.0, "d": 0.0}, -u_br)

    W = dg - df                       # total radial velocity amplitude at faces
    Rc = np.nan
    s = np.flatnonzero(np.diff(np.sign(W)))
    if s.size:
        Rc = brentq(CubicSpline(x_f, W), x_f[s[0]], x_f[s[0] + 1])
    return dict(R=R, Rc=Rc, RcR=Rc / R, x_f=x_f, x_c=x_c, f=f, g=g, W=W,
                T=(g - f) / x_c, u_f=u_f, u_br=u_br)


def cloud_closed_form(u_f, u_br):
    """E1.2's closed form.  Negative or nan means there is no cloud."""
    if u_br <= u_f:
        return np.nan
    return ((u_br + 2 * u_f) / (u_br - u_f)) ** (1 / 3)


# --- Geldart's own criteria, evaluated exactly as printed -------------------
def U_0_eq3(d_um, drho, p=P):
    """Equation (3), cm/s.  d_um in um, drho = rho_s - rho_f in g/cm3."""
    return p["C3"] * p["g"] * (d_um * 1e-4) ** 2 * drho / p["mu"]

def U_MB_eq2(d_um, p=P):
    """Equation (2), cm/s."""
    return p["K_MB"] * (d_um * 1e-4)

def u_br_davies_taylor(d_b, p=P):
    """0.711 sqrt(g d_b) - the coefficient E1.2 uses. NOT printed by Geldart."""
    return 0.711 * np.sqrt(p["g"] * d_b)

def u_br_geldart(d_b, p=P):
    """sqrt(g d_B / 2) - the left-hand side of Geldart's equation (7)."""
    return np.sqrt(p["g"] * d_b / 2.0)

print("solver and criteria defined")'''))

# ------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. The two boundary constants, recomputed

Equation (6) is equation (5) with $K_{MB} = 100$, $\mu = 1.8\times10^{-4}$ and
$d_{sv}$ in cm replaced by $d'$ in µm. Equation (8) is equation (7) with
$d_B = 25$ cm, $\varepsilon_0 = 0.4$, the same $\mu$, and the same substitution.
Both are pure arithmetic on printed inputs, so both can be redone."""))

cells.append(code(r'''def eq6_constant(p=P):
    """Constant on the RHS of eq. (6), from eq. (5).  Printed: 225.

    eq. (5):  C3 * g * d_sv * drho / (K_MB * mu) <= 1,  d_sv in cm
    with      d_sv = d' * 1e-4                          d'    in um
    =>        drho * d' <= K_MB * mu / (C3 * g * 1e-4)
    """
    return p["K_MB"] * p["mu"] / (p["C3"] * p["g"] * 1e-4)


def eq8_constant(p=P, u_br=None):
    """Constant on the RHS of eq. (8), from eq. (7).  Printed: 1e6.

    eq. (7):  sqrt(g d_B / 2) <= C3 * drho * g * d_sv^2 / (mu * eps_0)
    with      d_sv = d' * 1e-4
    =>        drho * d'^2 >= sqrt(g d_B/2) * mu * eps_0 / (C3 * g * 1e-8)
    `u_br` overrides the left-hand side, so the same routine gives the
    Davies-Taylor version used in section 5.
    """
    lhs = u_br_geldart(p["d_B"], p) if u_br is None else u_br
    return lhs * p["mu"] * p["eps_0"] / (p["C3"] * p["g"] * 1e-8)


C6, C8 = eq6_constant(), eq8_constant()
d6 = 100 * (C6 - PRINTED["eq6"]) / PRINTED["eq6"]
d8 = 100 * (C8 - PRINTED["eq8"]) / PRINTED["eq8"]

display(Markdown(
    f"Equation (6) recomputes to **{C6:.2f}** against the printed **225**, a "
    f"deviation of **{d6:+.2f} %**. Equation (8) recomputes to "
    f"**{C8:.4g}** against the printed **1e6**, **{d8:+.2f} %**. "
    f"Both are Geldart rounding his own arithmetic down to a memorable number, "
    f"and both round in the same direction. The question that matters is not "
    f"whether the roundings are defensible but whether they change any verdict; "
    f"section 2 answers that for the only boundary the data reach."))'''))

cells.append(md(r"""### 2. The A/B boundary against the 22 measured fractions

This is the test the page is for. Two verdicts per row, and they share no input
column:

- **the criterion**: group A when $(\rho_s-\rho_f)\,d' \le 225$. Uses the
  density and the size. Never sees a velocity.
- **the measurement**: group A when $U_{MB}/U_0 > 1$, equation (4) applied to
  the two measured velocity columns. Never sees a size or a density.

**The tie-break, stated up front, because it is load-bearing and the paper
prints the criterion both ways.** Two Diakon fractions print $U_0$ and $U_{MB}$
equal to the last digit — 3.11 against 3.11 at 263 µm, 4.11 against 4.11 at
318 µm — so their measured ratio is exactly 1.000 and the strict/non-strict
distinction alone decides them. *Displayed* as equation (4) on journal page 289
the criterion reads $U_{MB}/U_0 \ge 1$, which would put both in group A.
*Restated in words* one page earlier, inside the same section 4.3, it reads
"($U_{MB}/U_0 > 1$ for group A powders)" — and the sentence it sits inside says
those very fractions "should be classified as belonging to group B". **This page
takes the strict form**, `ratio > 1`, because it is the only reading that
reproduces the author's own verdict on his own rows. The choice is not free: the
non-strict reading changes every score below, and break table 2 reports by how
much.

**And a whole-table score is nearly worthless here, so it is not the headline.**
Almost all the comparable rows are measured group A, so a constant "everything
is group A" predictor already gets most of the table right without looking at
anything. The scoreboard below prints that null baseline beside the criterion's
score, names the rows on which the two can possibly differ, and scores the
criterion on those rows alone. That last number is the one to quote.

Nothing on the criterion side is fitted to the $U_0$ column, so a disagreement is
a real disagreement. There is one residual link — equation (2), which supplies
the 100 inside the 225, was fitted to the $U_{MB}$ column of this same table —
and the [Validation](#validation) section names it and reports equation (2)'s own
residual against that column."""))

cells.append(code(r'''def classify(df, threshold=None, size_col="d_sv_um", dens_col="drho", p=P):
    """Geldart group from eq. (6): 'A' when (rho_s - rho_f) d' <= threshold."""
    thr = PRINTED["eq6"] if threshold is None else threshold
    return np.where(df[dens_col] * df[size_col] <= thr, "A", "B")


# THE TIE-BREAK, fixed here and used everywhere on this page.
# Two Diakon rows print U_0 and U_MB equal to the last digit, so their measured
# ratio is exactly 1.000.  Equation (4) is DISPLAYED as U_MB/U_0 >= 1 (journal
# page 289) but RESTATED in words as "(U_MB/U_0 > 1 for group A powders)" in the
# same section 4.3 one page earlier - in the sentence that assigns those very
# fractions to group B.  We take the STRICT form, because it is the only reading
# that reproduces the author's verdict on his own rows.  Break table 2 reports
# what the non-strict form does to every score below.
t1["group_eq6"] = classify(t1)
t1["group_meas"]    = np.where(t1["ratio_meas"] >  1.0, "A", "B")   # strict, used
t1["group_meas_ge"] = np.where(t1["ratio_meas"] >= 1.0, "A", "B")   # eq. (4) as displayed
t1.loc[~has_U0, ["group_meas", "group_meas_ge"]] = "-"   # no U_0, no measured verdict

cmp_ = t1.loc[has_U0]
n_cmp = len(cmp_)
n_agree = int((cmp_.group_eq6 == cmp_.group_meas).sum())
disagree = cmp_.loc[cmp_.group_eq6 != cmp_.group_meas]
tie_rows = cmp_.loc[cmp_.ratio_meas == 1.0]

display(t1[["powder", "d_sv_um", "drho", "prod6", "group_eq6",
            "U_0", "U_MB", "ratio_meas", "group_meas"]].round(3))
print(f"\n{len(tie_rows)} rows print U_MB = U_0 exactly "
      f"({', '.join(f'{d:.0f} um' for d in tie_rows.d_sv_um)}), so the "
      f"tie-break alone decides them;\nthe strict reading of eq. (4) makes them "
      f"group B, which is what section 4.3 says they are")
print(f"\ncriterion agrees with the measurement on {n_agree} of {n_cmp} rows "
      f"with a measured U_0 - see the scoreboard below for what that is worth")
print("\nrows where they disagree:")
display(disagree[["powder", "nominal_range", "d_sv_um", "prod6",
                  "group_eq6", "ratio_meas", "group_meas"]].round(3))'''))

cells.append(code(r'''# --- What is that agreement worth?  Score it against a null model. ---------
# A predictor that says "group A" for everything sees no data at all.  Anything
# the criterion is worth is the margin over that, and the ONLY rows on which the
# two predictors can differ are the rows eq. (6) calls group B.
CMP    = has_U0.values
DIAKON = t1.rho_s_exact.values           # Diakon: the only series that CROSSES XY
ALLROWS = np.ones(len(t1), bool)

pred_eq6  = t1.group_eq6.values
pred_null = np.full(len(t1), "A")        # the null model
meas_gt   = t1.group_meas.values         # strict tie-break, this page's choice
meas_ge   = t1.group_meas_ge.values      # eq. (4) as displayed
DISC = CMP & (pred_eq6 != "A")           # the discriminating rows

def sc(pred, meas, mask):
    m = mask & CMP
    return int((pred[m] == meas[m]).sum()), int(m.sum())

def scs(pred, meas, mask):
    a, b = sc(pred, meas, mask); return f"{a}/{b}"

S = pd.DataFrame(
    [(lbl, scs(p, m, ALLROWS), scs(p, m, DIAKON), scs(p, m, DISC))
     for lbl, p, m in [
         ("criterion eq. (6),  tie-break >  (this page)", pred_eq6,  meas_gt),
         ('null model "all A", tie-break >',              pred_null, meas_gt),
         ("criterion eq. (6),  tie-break >= (as displayed)", pred_eq6,  meas_ge),
         ('null model "all A", tie-break >=',             pred_null, meas_ge)]],
    columns=["predictor", "all rows with a measured U_0", "Diakon only",
             "the rows eq. (6) calls B"])
display(S)

(a_all, _), (n_null, _)  = sc(pred_eq6, meas_gt, ALLROWS), sc(pred_null, meas_gt, ALLROWS)
(a_dk, n_dk), (d_null, _) = sc(pred_eq6, meas_gt, DIAKON), sc(pred_null, meas_gt, DIAKON)
(a_disc, n_disc)          = sc(pred_eq6, meas_gt, DISC)
(u_disc, _)               = sc(pred_null, meas_gt, DISC)
(ge_all, _), (ge_null, _) = sc(pred_eq6, meas_ge, ALLROWS), sc(pred_null, meas_ge, ALLROWS)
(ge_dk, _), (ge_disc, _)  = sc(pred_eq6, meas_ge, DIAKON), sc(pred_eq6, meas_ge, DISC)
n_measA = int((meas_gt[CMP] == "A").sum())
disc_sizes = ", ".join(f"{d:.0f}" for d in t1.d_sv_um[DISC])

display(Markdown(
    f"The criterion agrees with the measurement on **{a_all} of {n_cmp}** rows. "
    f"**That number should not be quoted on its own.** {n_measA} of the "
    f"{n_cmp} comparable rows are measured group A, so the null model — group A "
    f"for every powder, no arithmetic at all — already scores "
    f"**{n_null}/{n_cmp}**. The criterion is worth **{a_all - n_null} row** of "
    f"information over saying nothing.\n\n"
    f"The two predictors can differ only on the **{n_disc}** rows the criterion "
    f"calls group B: Diakon at {disc_sizes} µm. On those rows the criterion is "
    f"right **{a_disc} of {n_disc}**, against the null model's {u_disc} of "
    f"{n_disc}. ***{a_disc} of {n_disc} is the honest headline of this test.*** "
    f"The whole-table {a_all}/{n_cmp} is mostly a statement about how one-sided "
    f"the table is; Diakon alone, the only series that crosses the line, gives "
    f"{a_dk}/{n_dk} against the null model's {d_null}/{n_dk}.\n\n"
    f"All of that uses the strict tie-break. Read equation (4) as it is "
    f"displayed, $U_{{MB}}/U_0 \\ge 1$, and the two $U_{{MB}} = U_0$ rows become "
    f"group A: the criterion falls to {ge_all}/{n_cmp} (Diakon {ge_dk}/{n_dk}, "
    f"discriminating rows {ge_disc}/{n_disc}) while the null model rises to "
    f"**{ge_null}/{n_cmp}** — a perfect score for a predictor that has learned "
    f"nothing, which is by itself a reason to distrust that reading here. The "
    f"strict form is the one that reproduces Geldart's own sentence, and it is "
    f"what this page uses."))'''))

cells.append(md(r"""**The one disagreement is the row Geldart himself rules on — and he rules for
equation (6).** Section 4.3, journal page 288:

> It is interesting to note that the two largest sizes of Diakon bubbled at the
> incipient fluidization velocity and the 210–250-µm fraction at a velocity very
> close to $U_0$. According to the criterion ($U_{MB}/U_0 > 1$ for group A
> powders) these largest fractions should be classified as belonging to group B.

"These largest fractions" is all three of them: the two that bubbled at $U_0$
and the 210–250 µm fraction as well. So on that third row — actual size 220 µm —
there are three verdicts to compare, not two, and the cell below lays them out.
The short version is that **the strict ratio test is the odd one out**: equation
(6) says group B, Geldart's prose says group B, and only the mechanical
$U_{MB}/U_0 > 1$ says group A. He is reading the ratio loosely, and reasonably —
2.51 against 2.44 is "very close to $U_0$" for a bubble judged by eye — so the
single row on which the criterion and the measured ratio part company is a row
where the author of both had already said the ratio should not be read strictly,
and the boundary he kept is the one this page recomputes.

That is also, as section 5 of the paper says, *why* Diakon is in the table at
all — "it is particularly desirable to choose series of size fractions which
cross over the line representing eqn. (6). This was achieved with Diakon."
"""))

cells.append(code(r'''# The 220 um Diakon row, all three verdicts side by side.
row220 = t1.loc[t1.d_sv_um == 220].iloc[0]
V220 = pd.DataFrame([
    ("equation (6)", f"(rho_s - rho_f) d' = {row220.prod6:.1f} vs 225",
     row220.group_eq6),
    ("the strict ratio test, eq. (4)",
     f"U_MB/U_0 = {row220.U_MB}/{row220.U_0} = {row220.ratio_meas:.3f} vs 1",
     row220.group_meas),
    ("Geldart's prose, section 4.3",
     "the sentence quoted above, journal page 288", "B"),
], columns=["route", "what it evaluates", "verdict"])
display(V220)
display(Markdown(
    f"So on the one row where equation (6) and the measured ratio part company, "
    f"**Geldart's own prose sides with equation (6)**, not against it. The "
    f"disagreement this page reports is between the criterion and a strict "
    f"reading of the velocity columns — a reading its author explicitly declines "
    f"to make on this row. Two of the {n_disc} discriminating rows are the "
    f"$U_{{MB}} = U_0$ pair, where prose, criterion and strict ratio all agree on "
    f"group B; this is the third."))'''))

cells.append(code(r'''# How sharp is the boundary crossing?  Which rows sit closest to it, on
# either side, and does Geldart's 2 % rounding of 229 -> 225 move any of them?
C6_lo, C6_hi = min(PRINTED["eq6"], C6), max(PRINTED["eq6"], C6)
in_gap = t1.loc[(t1.prod6 > C6_lo) & (t1.prod6 <= C6_hi)]
near = t1.reindex((t1.prod6 - PRINTED["eq6"]).abs().sort_values().index).head(4)
near_A = t1.loc[t1.prod6 <= PRINTED["eq6"]].nlargest(1, "prod6").iloc[0]

# ... and the same question for eq. (8), which the page also claims no row
# falls into.  Claimed nowhere else on this page until now: compute it.
t1["prod8"] = t1.drho * t1.d_sv_um ** 2
C8_lo, C8_hi = min(PRINTED["eq8"], C8), max(PRINTED["eq8"], C8)
in_gap8 = t1.loc[(t1.prod8 > C8_lo) & (t1.prod8 <= C8_hi)]
d8_onset_1p5 = float(np.sqrt(PRINTED["eq8"] / (1.5 - P["rho_f"])))   # um, at rho_s = 1.5

display(near[["powder", "nominal_range", "d_sv_um", "prod6",
              "group_eq6", "ratio_meas", "group_meas"]].round(3))
display(Markdown(
    f"The last row on the group A side of the printed boundary sits at "
    f"**{near_A.prod6:.1f}**, which is "
    f"**{100*(PRINTED['eq6'] - near_A.prod6)/PRINTED['eq6']:.1f} %** "
    f"inside it, and it is measured group A "
    f"($U_{{MB}}/U_0 = {near_A.ratio_meas:.3f}$). "
    f"**{len(in_gap)}** rows of 22 fall in the gap between the printed 225 and "
    f"the recomputed {C6:.1f}, so Geldart's {d6:+.1f} % rounding changes no "
    f"verdict anywhere in his own table.\n\n"
    f"The same for equation (8), which is asserted elsewhere on this page and is "
    f"checked here rather than assumed: the largest "
    f"$(\\rho_s-\\rho_f)(d')^2$ in Table 1 is "
    f"**{t1.prod8.max():.3g}**, a factor {PRINTED['eq8']/t1.prod8.max():.1f} "
    f"below the printed $10^6$, so **{len(in_gap8)}** rows fall in the gap "
    f"between $10^6$ and the recomputed {C8:.4g} — and none is anywhere near it. "
    f"At the table's own densest powder, $\\rho_s = 1.5$ g cm$^{{-3}}$, equation "
    f"(8) does not begin until $d' = {d8_onset_1p5:.0f}$ µm, against a largest "
    f"measured fraction of {t1.d_sv_um.max():.0f} µm."))'''))

cells.append(md(r"""### 3. What equation (3) is worth against the measured $U_0$

The A/B boundary is the ratio of two correlations, and only one of them was
fitted to this table. Equation (3) — Davies and Richardson's, via Geldart's
reference 11 — is the load-bearing half, and the measured $U_0$ column is a
direct test of it that the paper never performs.

Rather than quote a scatter, invert it: solve equation (3) for the constant that
each measured row implies,

$$K_0^{\text{eff}} = \frac{U_0\,\mu}{g\,d_{sv}^2\,(\rho_s-\rho_f)},$$

and compare with the printed $8\times10^{-4}$. This is the same comparison as
$U_0^{\text{eq3}}$ against $U_0^{\text{meas}}$ but on a scale where the answer is
a single dimensionless number per row."""))

cells.append(code(r'''t1["U0_eq3"] = U_0_eq3(t1.d_sv_um, t1.drho)
t1["K0_eff"] = t1.U_0 * P["mu"] / (P["g"] * (t1.d_sv_um * 1e-4) ** 2 * t1.drho)
t1["K0_ratio"] = t1.K0_eff / P["C3"]

display(t1.loc[has_U0, ["powder", "d_sv_um", "U_0", "U0_eq3",
                        "K0_eff", "K0_ratio"]].round(4))

k = t1.loc[has_U0, "K0_ratio"]
kd = t1.loc[has_U0 & t1.rho_s_exact, "K0_ratio"]                 # Diakon only
kd_fine = t1.loc[has_U0 & t1.rho_s_exact & (t1.d_sv_um > 100), "K0_ratio"]
display(Markdown(
    f"Across all **{len(k)}** measured fractions the constant that equation (3) "
    f"would need runs from **{k.min():.2f}** to **{k.max():.2f}** times the "
    f"printed $8\\times10^{{-4}}$ — a spread of a factor **{k.max()/k.min():.1f}**, "
    f"with a median of {k.median():.2f}. Restricted to **Diakon**, the only "
    f"series whose density is printed exactly, it runs "
    f"**{kd.min():.2f}** to **{kd.max():.2f}**; dropping Diakon's finest "
    f"fraction (78 µm, the one outlier) leaves "
    f"**{kd_fine.min():.2f}** to **{kd_fine.max():.2f}**, i.e. within about "
    f"{100*max(1-kd_fine.min(), kd_fine.max()-1):.0f} % of the printed constant. "
    f"The one thing that comes out well is the *centre*: the median over all "
    f"{len(k)} fractions is {k.median():.3f} times the printed value, so "
    f"$8\\times10^{{-4}}$ is the right number for this collection of powders "
    f"to within {100*abs(k.median()-1):.0f} % — it is the row-to-row scatter, "
    f"not the constant, that is the problem. "
    f"So equation (3) is good to roughly ±{100*max(1-kd_fine.min(), kd_fine.max()-1):.0f} % "
    f"on a clean series of near-spherical particles above 100 µm, and much worse "
    f"below that and on the two powders whose density Geldart prints as "
    f"approximate."))'''))

cells.append(code(r'''fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
mk = {"Diakon": "s", "Fresh catalyst": "o", "Spent catalyst": "^"}
for name, grp in t1.loc[has_U0].groupby("powder"):
    ax[0].loglog(grp.U0_eq3, grp.U_0, mk[name], ms=6, label=name)
    ax[1].semilogx(grp.d_sv_um, grp.K0_ratio, mk[name], ms=6, label=name)
lim = [0.02, 6]
ax[0].plot(lim, lim, "k-", lw=1, label="parity")
ax[0].set(xlim=lim, ylim=lim, xlabel=r"$U_0$ from eq. (3)  (cm/s)",
          ylabel=r"$U_0$ measured, Table 1  (cm/s)",
          title="equation (3) against the measurement")
ax[1].axhline(1.0, color="k", lw=1)
ax[1].set(xlabel=r"$d_{sv}$  ($\mu$m)",
          ylabel=r"$K_0^{\rm eff} / (8\times10^{-4})$",
          title="the constant each row implies")
ax[1].set_xticks([25, 50, 100, 200, 300])
ax[1].set_xticklabels(["25", "50", "100", "200", "300"])
for a in ax:
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""### 4. The bed-expansion columns, and the voidage Geldart does not print

Table 1 carries two more columns, $\varepsilon_{MB}$ and $H_{MB}/H_0$. **The
paper does not say how either was obtained** — section 4.2 describes the
measurement of $U_0$ and $U_{MB}$ and stops there — so they are treated here as
*reported*, and what follows is a transcription check on them and nothing more.
They are not independent of each other. The solids in the bed
are conserved between the settled state and the minimum bubbling state, so

$$H_0\,(1-\varepsilon_0) = H_{MB}\,(1-\varepsilon_{MB})
\qquad\Longrightarrow\qquad
\varepsilon_0 = 1 - (1-\varepsilon_{MB})\,\frac{H_{MB}}{H_0}.$$

Geldart prints neither $\varepsilon_0$ nor $H_0$, so this recovers a number the
paper does not state — and it is simultaneously a check on the two columns,
because an $\varepsilon_0$ outside about 0.35–0.65 would mean one of them is
mis-transcribed. Two rows make it sharper still: the two coarsest Diakon
fractions print $H_{MB}/H_0 = 1.000$, which forces $\varepsilon_0 =
\varepsilon_{MB}$ exactly."""))

cells.append(code(r'''t1["eps_0_inferred"] = 1.0 - (1.0 - t1.eps_MB) * t1.H_MB_over_H_0
NOEXP = t1.index[t1.H_MB_over_H_0 == 1.000]      # fixed from the PRINTED table
resid_noexp = float(np.max(np.abs(t1.loc[NOEXP, "eps_0_inferred"]
                                  - t1.loc[NOEXP, "eps_MB"])))
dia = t1.loc[t1.rho_s_exact]

display(t1[["powder", "d_sv_um", "eps_MB", "H_MB_over_H_0",
            "eps_0_inferred"]].round(4))
display(Markdown(
    f"All 22 rows give $\\varepsilon_0$ between **{t1.eps_0_inferred.min():.3f}** "
    f"and **{t1.eps_0_inferred.max():.3f}**, every one of them a physically "
    f"admissible packed-bed voidage. Diakon, the near-spherical powder, is "
    f"tightest at **{dia.eps_0_inferred.min():.3f}–"
    f"{dia.eps_0_inferred.max():.3f}**: it falls steadily from "
    f"{dia.eps_0_inferred.iloc[0]:.3f} at {dia.d_sv_um.iloc[0]:.0f} µm to "
    f"{dia.eps_0_inferred.iloc[5]:.3f} at {dia.d_sv_um.iloc[5]:.0f} µm and is "
    f"flat at {dia.eps_0_inferred.iloc[-1]:.3f} for the two coarsest fractions, "
    f"which is what a powder that packs tighter as it coarsens should do. The "
    f"two catalysts, irregular and porous, sit higher at "
    f"**{t1.loc[~t1.rho_s_exact,'eps_0_inferred'].min():.3f}–"
    f"{t1.loc[~t1.rho_s_exact,'eps_0_inferred'].max():.3f}**. The "
    f"{len(NOEXP)} rows that print $H_{{MB}}/H_0 = 1.000$ return "
    f"$\\varepsilon_0 = \\varepsilon_{{MB}}$ to **{resid_noexp:.1e}** — but that "
    f"one is *structural*: with $H_{{MB}}/H_0$ exactly 1 the identity reduces to "
    f"$\\varepsilon_0 = \\varepsilon_{{MB}}$ algebraically, so it tests only that "
    f"those two cells were transcribed as 1.000 and nothing else.\n\n"
    f"Note what this does **not** say. Geldart uses $\\varepsilon_0 = 0.4$ in "
    f"equation (7); the powders in his own table sit at "
    f"{t1.eps_0_inferred.min():.2f}–{t1.eps_0_inferred.max():.2f}, well above "
    f"it. That is not a contradiction — 0.4 is stated for *large-particle* "
    f"systems and Table 1 contains none — but it does mean the table cannot be "
    f"used to check the $\\varepsilon_0$ that goes into the B/D boundary.\n\n"
    f"And note what it *rests* on. If $\\varepsilon_{{MB}}$ was itself obtained "
    f"from bed height, the 200 g charge and $\\rho_s$ — the natural route, though "
    f"the paper does not say — then for the "
    f"{int((~t1.rho_s_exact).sum())} catalyst rows it inherits a $\\rho_s$ "
    f"printed as approximate, and so does the $\\varepsilon_0$ recovered from it. "
    f"\"Physically admissible\" above therefore means admissible given the "
    f"printed numbers, not verified."))'''))

cells.append(md(r"""### 5. The B/D boundary is a cloud-existence condition

Read the two sides of equation (7) again:

$$\underbrace{\left(g\,\frac{d_B}{2}\right)^{1/2}}_{\text{bubble rise velocity}}
\;\le\;
\underbrace{\frac{8\times10^{-4}(\rho_s-\rho_f)\,g\,d_{sv}^2}{\mu\,\varepsilon_0}}_{
  \;=\;U_0/\varepsilon_0\;=\;\text{interstitial gas velocity}}$$

The right-hand side is equation (3) divided by $\varepsilon_0$ — the interstitial
velocity, which E1.2 calls $u_f$. The left-hand side is $\tfrac{1}{\sqrt2}
\sqrt{g\,d_B} = 0.7071\sqrt{g\,d_B}$, and the Davies–Taylor rise velocity that
E1.2 uses is $u_{br} = 0.711\sqrt{g\,d_b}$. **They are the same formula to half a
percent.** Geldart's group D condition is therefore

$$u_{br} \le u_f \qquad\text{for } d_B = 25\ \text{cm},$$

which is precisely E1.2's no-cloud condition — the pole in
$(R_c/R)^3 = (u_{br}+2u_f)/(u_{br}-u_f)$.

**What is new here, and what is not.** That a bubble rising more slowly than the
interstitial gas carries no cloud is standard — E1.2 says it in those words, and
Geldart's section 3.4 and journal page 289 describe group D in exactly the same
physical terms, gas entering the base of the bubble and leaving the top. Geldart
never mentions Davidson, and E1.2's source (Kunii and Levenspiel 1968) predates
the classification, but the *equivalence in words* is in both places. What is
absent from both is the arithmetic: identifying equation (7)'s left-hand side
with the Davies–Taylor rise velocity and its right-hand side with $u_f$, and so
reading the printed group D line as E1.2's pole. Read that way — which is a
re-reading, not a new result — **the whole Geldart diagram is also a map of the
smallest bubble that can carry a cloud**,

$$d_b^{*} = \frac{1}{g}\left(\frac{u_f}{0.711}\right)^2
          = \frac{1}{g}\left(\frac{U_0}{0.711\,\varepsilon_0}\right)^2,$$

whose 25 cm contour is Geldart's line for group D. Below is that identity
checked three ways: the two rise-velocity coefficients compared, the contour
compared with equation (8), and a bubble actually solved on each side of the
line."""))

cells.append(code(r'''# (a) the two coefficients, and what the difference does to the boundary
c_geldart = u_br_geldart(P["d_B"]) / np.sqrt(P["g"] * P["d_B"])   # = 1/sqrt(2)
c_davies  = 0.711
C8_dt = eq8_constant(u_br=u_br_davies_taylor(P["d_B"]))

display(Markdown(
    f"Geldart's coefficient is $1/\\sqrt2 = {c_geldart:.4f}$; Davies–Taylor is "
    f"{c_davies:.4f}, a difference of "
    f"**{100*(c_davies-c_geldart)/c_geldart:.2f} %** in velocity. Carried "
    f"through equation (7) that moves the boundary constant from "
    f"**{C8:.4g}** to **{C8_dt:.4g}**, "
    f"**{100*(C8_dt-C8)/C8:.2f} %** — against the **{d8:.2f} %** by which "
    f"Geldart rounded his own {C8:.4g} to $10^6$. *The difference between the "
    f"two theories is smaller than the difference between Geldart's arithmetic "
    f"and Geldart's printed answer.*"))'''))

cells.append(code(r'''# (b) two real powders, one on each side of the B/D line, same 25 cm bubble.
BEDS = {
    "sand, 100 um (B, far inside)":   dict(d_um=100.0, rho_s=2.6),
    "sand, 600 um (B, near the line)": dict(d_um=600.0, rho_s=2.6),
    "sand, 700 um (group D)":         dict(d_um=700.0, rho_s=2.6),
}
rows = []
for name, b in BEDS.items():
    drho = b["rho_s"] - P["rho_f"]
    U0   = U_0_eq3(b["d_um"], drho)
    u_f  = U0 / P["eps_0"]
    u_br = u_br_davies_taylor(P["d_B"])
    sol  = davidson_cloud(P["d_B"], u_f, u_br)
    rows.append(dict(bed=name,
                     prod_eq6=drho * b["d_um"], prod_eq8=drho * b["d_um"] ** 2,
                     group=("A" if drho * b["d_um"] <= PRINTED["eq6"] else
                            "D" if drho * b["d_um"] ** 2 >= PRINTED["eq8"] else "B"),
                     U0=U0, u_f=u_f, u_br=u_br, ratio=u_br / u_f,
                     RcR_pymrm=sol["RcR"], RcR_closed=cloud_closed_form(u_f, u_br)))
    BEDS[name]["sol"] = sol
bd = pd.DataFrame(rows).set_index("bed")
display(bd.round(4))

ok = bd.dropna(subset=["RcR_closed"])
dev_cloud = float(np.max(np.abs(ok.RcR_pymrm - ok.RcR_closed) / ok.RcR_closed))
display(Markdown(
    f"Both group B powders get a cloud and the group D powder does not, and the "
    f"pymrm solve reports the absence as a **missing dividing streamline**, not "
    f"as a negative number. The cloud is a thin shell far inside group B "
    f"($R_c/R = {bd.RcR_pymrm.iloc[0]:.3f}$ at 100 µm) and swells as the powder "
    f"approaches the line ($R_c/R = {bd.RcR_pymrm.iloc[1]:.2f}$ at 600 µm, where "
    f"$u_{{br}}/u_f$ is only {bd.ratio.iloc[1]:.2f}) — the pole, seen from the "
    f"inside. Where a cloud exists the computed $R_c/R$ matches E1.2's closed "
    f"form to **{dev_cloud:.1e}**; that is discretisation error against a "
    f"formula E1.2 already validated, not new evidence, and it is reported only "
    f"to show the port is faithful."))'''))

cells.append(code(r'''# (c) the classification map, drawn ONLY from the printed equations.
d_grid    = np.logspace(np.log10(10), np.log10(3000), 400)      # um
drho_grid = np.logspace(np.log10(0.1), np.log10(10), 400)       # g/cm3
DD, RR = np.meshgrid(d_grid, drho_grid)
db_star = (U_0_eq3(DD, RR) / (0.711 * P["eps_0"])) ** 2 / P["g"]   # cm

fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.0))

# --- left: the diagram itself
a = ax[0]
a.loglog(d_grid, PRINTED["eq6"] / d_grid, "k-", lw=2,
         label=r"eq. (6) as printed:  $\Delta\rho\,d' = 225$")
a.loglog(d_grid, C6 / d_grid, "k--", lw=1,
         label=rf"eq. (6) recomputed:  ${C6:.0f}$")
a.loglog(d_grid, PRINTED["eq8"] / d_grid ** 2, "-", color="tab:red", lw=2,
         label=r"eq. (8) as printed:  $\Delta\rho\,d'^2 = 10^6$")
a.loglog(d_grid, C8 / d_grid ** 2, "--", color="tab:red", lw=1,
         label=rf"eq. (8) recomputed:  ${C8:.3g}$")
C10 = PRINTED["eq10"] * P["mu"] / (np.sqrt(P["g"]) * 1e-6)
a.loglog(d_grid, C10 / d_grid ** 1.5, ":", color="tab:blue", lw=2,
         label=r"eq. (10), Oltrogge (line O-O)")
for name, grp in t1.groupby("powder"):
    for gletter, style in (("A", dict(mfc="w")), ("B", dict(mfc="tab:green")),
                           ("-", dict(mfc="0.8"))):
        s = grp[grp.group_meas == gletter]
        if len(s):
            a.loglog(s.d_sv_um, s.drho, mk[name], ms=7, mec="k", lw=0,
                     label=None, **style)
for name in mk:                                   # marker key only, no data
    a.loglog([], [], mk[name], ms=7, mec="k", mfc="w", lw=0,
             label=f"Table 1: {name.lower()}")
a.set(xlim=(10, 3000), ylim=(0.1, 10),
      xlabel=r"mean particle size  $d_{sv}$  ($\mu$m)",
      ylabel=r"density difference  $\rho_s-\rho_f$  (g/cm$^3$)",
      title="the boundaries Geldart prints, plus his own Table 1")
for lbl, xy in [("A", (30, 3.5)), ("B", (330, 2.2)), ("D", (1600, 3.5))]:
    a.text(*xy, lbl, fontsize=15, weight="bold", ha="center")
a.legend(fontsize=7.5, loc="lower left")
a.text(0.98, 0.97, "filled = measured group B\nopen = measured group A\n"
       "faint = no measured $U_0$", transform=a.transAxes, fontsize=7,
       va="top", ha="right")

# --- right: the same plane as a cloud map
a = ax[1]
lv = [0.1, 1.0, 5.0, 100.0]
cs = a.contour(DD, RR, db_star, levels=lv, colors="0.4", linewidths=1)
a.clabel(cs, fmt=lambda v: f"{v:g} cm", fontsize=7)
a.contour(DD, RR, db_star, levels=[P["d_B"]], colors="tab:purple", linewidths=2.5)
a.loglog(d_grid, PRINTED["eq8"] / d_grid ** 2, "-", color="tab:red", lw=1.5,
         label=r"eq. (8) as printed")
a.loglog(d_grid, C8_dt / d_grid ** 2, "--", color="tab:purple", lw=1.5,
         label=rf"$d_b^*={P['d_B']:.0f}$ cm contour  (${C8_dt:.3g}$)")
a.set(xscale="log", yscale="log", xlim=(10, 3000), ylim=(0.1, 10),
      xlabel=r"mean particle size  $d_{sv}$  ($\mu$m)",
      ylabel=r"density difference  $\rho_s-\rho_f$  (g/cm$^3$)",
      title=r"smallest bubble that carries a cloud, $d_b^*$")
a.legend(fontsize=8, loc="lower left")
fig.tight_layout()
plt.show()

print(f"eq. (10) in (drho, d') form:  drho * d'^1.5 = {C10:.0f}")'''))

cells.append(code(r'''# (d) the cloud around the same 25 cm bubble in each powder, from the field.
fig, ax = plt.subplots(1, 3, figsize=(13, 4.6))
th = np.linspace(0, 2 * np.pi, 400)
for a, (name, b) in zip(ax, BEDS.items()):
    sol = b["sol"]; R = sol["R"]
    rmax = 4.5 * R
    rr = np.linspace(R, rmax, 400)
    Wi = np.interp(rr, sol["x_f"], sol["W"])
    TH, RG = np.meshgrid(np.linspace(0, np.pi, 300), rr)
    psi = 0.5 * RG ** 2 * Wi[:, None] * np.sin(TH) ** 2
    # levels chosen from the equatorial value of psi, so the streamlines are
    # spread evenly around the bubble instead of bunching in the far field
    r_seed = np.linspace(1.02 * R, 0.95 * rmax, 16)
    lv = np.unique(np.sort(np.interp(r_seed, rr, 0.5 * rr ** 2 * Wi)))
    lv = np.unique(np.concatenate([lv, [0.0]]))
    for sgn in (1, -1):
        a.contour(sgn * RG * np.sin(TH), RG * np.cos(TH), psi,
                  levels=lv, colors="0.55", linewidths=0.7, linestyles="solid")
    a.fill(R * np.cos(th), R * np.sin(th), color="w", ec="k", lw=1.8, zorder=3)
    if np.isfinite(sol["Rc"]):
        a.plot(sol["Rc"] * np.cos(th), sol["Rc"] * np.sin(th), "-",
               color="tab:red", lw=2.2, zorder=4)
        sub = f"cloud at $R_c/R$ = {sol['RcR']:.3f}"
    else:
        sub = "no cloud: no dividing streamline"
    lim = 2.6 * R if not np.isfinite(sol["Rc"]) else max(2.6 * R, 1.2 * sol["Rc"])
    a.set(aspect="equal", xlim=(-lim, lim), ylim=(-lim, lim),
          title=f"{name}\n$u_{{br}}/u_f$ = {sol['u_br']/sol['u_f']:.2f}\n{sub}")
    a.grid(False)
fig.suptitle(f"the same {P['d_B']:.0f} cm bubble in three powders "
             "(streamlines in the bubble frame, pymrm; red = cloud boundary)",
             y=1.04)
fig.tight_layout()
plt.show()'''))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

Ranked as the gallery's brief ranks them, highest available first.

**Route 3 — stated numerical results in the text.** Geldart works one example on
journal page 289 with printed intermediates, and states a size at which a
1 g cm$^{-3}$ powder turns from A to B. Both are reproduced below.

**Route 3 — against measurement.** The A/B criterion against the 22 measured
fractions of Table 1, and equation (3) against the measured $U_0$ column. These
are the only comparisons on the page that involve a measurement.

**Route 2 — internal identities.** The two boundary constants recomputed from
the expressions they are derived from; the solids balance linking
$\varepsilon_{MB}$ and $H_{MB}/H_0$; the placement of equations (9) and (10)
relative to (6), against Geldart's prose about where they lie.

**Route 4 — digitised figures.** Not used. Nothing on this page is digitised."""))

cells.append(code(r'''# --- V1: the worked example, journal page 289 -------------------------------
d_ex, drho_ex = 100.0, 1.0
U0_ex, UMB_ex = U_0_eq3(d_ex, drho_ex), U_MB_eq2(d_ex)
U0_sand = U_0_eq3(d_ex, 2.7 - P["rho_f"])
d_cross = C6 / 1.0                       # drho = 1 g/cm3, from eq. (6)

V1 = pd.DataFrame([
    ("U_0 at d'=100 um, drho=1 g/cm3 (eq. 3)", U0_ex, 0.43,
     "read off Fig. 2 by Geldart, so printed to 2 s.f."),
    ("U_MB at d'=100 um (eq. 2)",              UMB_ex, 1.0, "exact by construction"),
    ("U_MB/U_0 at that point",                 UMB_ex / U0_ex, 2.33,
     "Geldart's 2.33 = 1/0.43, i.e. his graph-read U_0"),
    ("U_0 for 100 um sand, rho_s=2.7 (eq. 3)", U0_sand, 1.2, "printed to 2 s.f."),
    ("d' at which drho=1 turns A -> B",        d_cross, 250.0,
     'Geldart writes "about 250 um"'),
], columns=["quantity", "recomputed", "printed", "note"])
V1["dev_%"] = 100 * (V1.recomputed - V1.printed) / V1.printed
display(V1[["quantity", "recomputed", "printed", "dev_%", "note"]].round(4))

display(Markdown(
    f"Four of the five land within "
    f"{V1['dev_%'].abs().iloc[:4].max():.1f} % of the printed value, which is "
    f"all two significant figures can resolve: the worst of the four is "
    f"*{V1.quantity[V1['dev_%'].abs().iloc[:4].idxmax()]}* at "
    f"{V1['dev_%'].iloc[:4].loc[V1['dev_%'].abs().iloc[:4].idxmax()]:+.2f} %. "
    f"The ratio 2.33 is Geldart dividing by a number he read off a log-log "
    f"graph: $1/0.43 = {1/0.43:.3f}$ against $1/{U0_ex:.4f} = "
    f"{UMB_ex/U0_ex:.3f}$. "
    f"The fifth is different in kind. His \"about 250 µm\" is "
    f"{100*(d_cross-250)/250:+.1f} % from the {d_cross:.0f} µm his own equation "
    f"(6) gives, and $250 \\times 1 = 250$ is on the group B side of both 225 "
    f"and {C6:.0f}. It is a loose sentence, not a slip in an equation: nothing "
    f"downstream uses it."))'''))

cells.append(code(r'''# --- V2: where equations (9) and (10) sit, against Geldart's prose ---------
# (g d^3)^0.5 drho / mu > K  with d = d' * 1e-4 cm  =>  drho * d'^1.5 > K mu / (sqrt(g) 1e-6)
def const_9_10(K, p=P):
    return K * p["mu"] / (np.sqrt(p["g"]) * 1e-6)

C9, C10 = const_9_10(PRINTED["eq9"]), const_9_10(PRINTED["eq10"])
# where each crosses eq. (6):  drho d' = 225  and  drho d'^1.5 = C
cross9  = (C9 / PRINTED["eq6"]) ** 2
cross10 = (C10 / PRINTED["eq6"]) ** 2

# The LABELLED ticks of Fig. 3 are printed text, not digitised coordinates:
# 20, 50, 100, 200, 500, 1000 um on the abscissa; 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7
# g/cm3 on the ordinate.  Nothing else about the figure is used.
FIG3_D, FIG3_RHO = (20.0, 1000.0), (0.2, 7.0)
place = pd.DataFrame(
    [(r, (C9 / r) ** (2 / 3), (C10 / r) ** (2 / 3), PRINTED["eq6"] / r,
      ((C9 / r) ** (2 / 3)) / (PRINTED["eq6"] / r))
     for r in (FIG3_RHO[1], 1.0, FIG3_RHO[0])],
    columns=["drho (g/cm3)", "eq. (9)  d' (um)", "eq. (10)  d' (um)",
             "line XY  d' (um)", "eq.(9) / XY, in size"])
display(place.round(2))

d9_at_1 = (C9 / 1.0) ** (2 / 3)
d9_if_500 = (const_9_10(500.0) / 1.0) ** (2 / 3)      # eq. (9) with 5000 mis-read
display(Markdown(
    f"In $(\\Delta\\rho, d')$ coordinates the two competing criteria are "
    f"$\\Delta\\rho\\,d'^{{1.5}} = {C9:.0f}$ (Verloop and Heertjes, eq. 9) and "
    f"$\\Delta\\rho\\,d'^{{1.5}} = {C10:.0f}$ (Oltrogge, eq. 10). Both are "
    f"steeper in $d'$ than Geldart's line XY, so each meets XY exactly once:\n\n"
    f"* **eq. (10) crosses XY at $d' = {cross10:.0f}$ µm, $\\Delta\\rho = "
    f"{PRINTED['eq6']/cross10:.2f}$ g/cm³** — inside both of Figure 3's labelled "
    f"axis ranges ({FIG3_D[0]:.0f}–{FIG3_D[1]:.0f} µm, "
    f"{FIG3_RHO[0]}–{FIG3_RHO[1]:.0f} g/cm³), consistent with Geldart drawing it "
    f"there as line O–O and calling it as good as eqn. (5);\n"
    f"* **eq. (9) crosses XY at $d' = {cross9:.0f}$ µm, $\\Delta\\rho = "
    f"{PRINTED['eq6']/cross9:.4f}$ g/cm³** — far outside the diagram, and far "
    f"*below* its ordinate rather than beyond its abscissa. **That crossing "
    f"point is therefore not evidence about what is on the page**, and this "
    f"page previously mis-used it as such.\n\n"
    f"What the crossing does establish is that the two lines never meet inside "
    f"the diagram, so eq. (9) lies to the right of XY throughout it — by a "
    f"factor {place['eq.(9) / XY, in size'].iloc[0]:.1f} in size at the top of "
    f"the labelled ordinate and {place['eq.(9) / XY, in size'].iloc[-1]:.1f} at "
    f"the bottom. The line itself is **not off-scale**: at "
    f"$\\Delta\\rho = 1$ g/cm³ it sits at {d9_at_1:.0f} µm against XY's "
    f"{PRINTED['eq6']:.0f} µm, both inside the labelled abscissa. So Geldart's "
    f"\"lies much too far to the right\" is a statement about its distance from "
    f"XY and from his data points, not about it running off the axis; that it "
    f"\"is not shown on Fig. 3\" is his statement about his own drawing, and "
    f"nothing here checks it.\n\n"
    f"The placement has teeth: read the 5000 of eq. (9) as 500 and it lands at "
    f"{d9_if_500:.0f} µm at $\\Delta\\rho = 1$, i.e. to the *left* of XY, "
    f"inverting the conclusion entirely."))'''))

cells.append(md(r"""### Make sure the checks can fail

Every agreement number above is now attacked. Each row injects one defect that
the check is *supposed* to catch, and reports what the check does. A defect that
leaves a number unmoved is a blind spot and is named as one.

The three tables cover, in order: the boundary constants (V1 and section 1), the
Table 1 classification (section 2) and the two column-based checks (sections 3
and 4)."""))

cells.append(code(r'''# --- BREAK TABLE 1: the boundary constants ---------------------------------
def broken(**kw):
    p = dict(P); p.update(kw); return p

defects = [
    ("baseline",                                  dict()),
    ("mu read as 1.8e-3 (decimal slip)",          dict(mu=1.8e-3)),
    ("mu read as 1.8e-5",                         dict(mu=1.8e-5)),
    ("K_MB read as 10 instead of 100",            dict(K_MB=10.0)),
    ("the 8e-4 in eq. (3) read as 8e-3",          dict(C3=8.0e-3)),
    ("g in SI (9.81 m/s2) not converted",         dict(g=9.81)),
    ("eps_0 = 0.5 instead of 0.4",                dict(eps_0=0.5)),
    ("d_B = 50 cm instead of 25",                 dict(d_B=50.0)),
    ("rho_f dropped entirely (rho_f = 0)",        dict(rho_f=0.0)),
]
rows = []
for label, kw in defects:
    p = broken(**kw)
    c6, c8 = eq6_constant(p), eq8_constant(p)
    rows.append((label, c6, 100 * (c6 - PRINTED["eq6"]) / PRINTED["eq6"],
                 c8, 100 * (c8 - PRINTED["eq8"]) / PRINTED["eq8"]))
B1 = pd.DataFrame(rows, columns=["injected defect", "eq6 const", "dev_% vs 225",
                                 "eq8 const", "dev_% vs 1e6"])
display(B1.round(2))
print("\nBlind spot, and it is a real one: rho_f enters neither constant, "
      "because\nboth are written in the density DIFFERENCE and Geldart never "
      "substitutes a value\nfor rho_f. Setting rho_f = 0 moves nothing here. "
      "It moves the Table 1\nclassification by 0.1 % of one product, which is "
      "also nothing. This page therefore\nmakes NO claim that rho_f was read "
      "or used correctly - nothing on it depends on rho_f.")'''))

cells.append(code(r'''# --- BREAK TABLE 2: the Table 1 A/B classification -------------------------
def agreement_rows(threshold=None, dens_col="drho", size_col="d_sv_um",
                   swap_velocities=False, invert=False, df=None, subset=None,
                   strict=True, null_model=False):
    d = (t1 if df is None else df)
    if subset is not None:
        d = d.loc[subset]
    if null_model:
        pred = np.full(len(d), "A")
    else:
        pred = classify(d, threshold=threshold, dens_col=dens_col,
                        size_col=size_col)
    if invert:
        pred = np.where(pred == "A", "B", "A")
    ratio = (d["U_0"] / d["U_MB"]) if swap_velocities else (d["U_MB"] / d["U_0"])
    meas = np.where(ratio > 1.0 if strict else ratio >= 1.0, "A", "B")
    m = d["U_0"].notna().values
    return int((pred[m] == meas[m]).sum()), int(m.sum())

t1_cm = t1.copy()
t1_cm["d_sv_cm"] = t1_cm.d_sv_um * 1e-4          # the classic units slip

cases = [
    ("baseline: eq. (6) as printed, 225, tie-break >", dict()),
    ("TIE-BREAK: eq. (4) applied as displayed, >= not >", dict(strict=False)),
    ("REFERENCE: null model, 'group A' for everything", dict(null_model=True)),
    ("boundary recomputed, 229.4 instead of 225", dict(threshold=C6)),
    ("boundary mis-read as 150",                  dict(threshold=150.0)),
    ("boundary mis-read as 500",                  dict(threshold=500.0)),
    ("boundary mis-read as 2250 (decimal slip)",  dict(threshold=2250.0)),
    ("criterion inverted (A when product > 225)", dict(invert=True)),
    ("rho_s used in place of (rho_s - rho_f)",    dict(dens_col="rho_s")),
    ("d' left in cm, not converted to um",        dict(size_col="d_sv_cm", df=t1_cm)),
    ("U_0 and U_MB columns swapped",              dict(swap_velocities=True)),
]
rows = []
for label, kw in cases:                   # local names only: n_dk etc. are live
    b2_all = agreement_rows(**kw)
    kwd = dict(kw); kwd["subset"] = DIAKON
    if "df" in kwd:                       # keep the modified frame, subset it too
        kwd["df"] = kwd["df"]
    b2_dk = agreement_rows(**kwd)
    rows.append((label, f"{b2_all[0]}/{b2_all[1]}", f"{b2_dk[0]}/{b2_dk[1]}"))
B2 = pd.DataFrame(rows, columns=["injected defect", "all rows agreeing",
                                 "Diakon only"])
display(B2)
print("Two rows here are not defects and are the most important in the table.\n\n"
      "The TIE-BREAK row is the fork described in the results section: eq. (4)\n"
      "is DISPLAYED as U_MB/U_0 >= 1, and applying it that way costs the two\n"
      "rows that print U_MB = U_0 exactly. Section 4.3 says those fractions are\n"
      "group B, so the strict reading is the right one - but a reader who\n"
      "reimplements from the displayed inequality will get this row's numbers,\n"
      "not the baseline's, and would otherwise think the baseline invented.\n\n"
      "The REFERENCE row is a predictor that ignores the data entirely. Every\n"
      "score in this table has to be read against it, not against a notional\n"
      "perfect score - the classes here are badly unbalanced.\n\n"
      "The Diakon column is the one with power. 14 of the 22 fractions sit so\n"
      "far inside group A that no plausible mis-reading of the boundary moves\n"
      "them, so the whole-table score changes little; Diakon is the only series\n"
      "that CROSSES the line - which is why Geldart chose it, section 5:\n"
      '"It is particularly desirable to choose series of size fractions which\n'
      'cross over the line representing eqn. (6). This was achieved with Diakon."')'''))

cells.append(code(r'''# --- BREAK TABLE 3: the two column-based checks ----------------------------
def k0_span(exponent=2.0, dens_col="drho", df=None):
    d = (t1 if df is None else df)
    k = d.U_0 * P["mu"] / (P["g"] * (d.d_sv_um * 1e-4) ** exponent * d[dens_col])
    k = (k / P["C3"]).dropna()
    return k.median(), k.min(), k.max(), k.max() / k.min()

def eps0_span(df=None):
    """Range test + the H/H0 = 1 identity, on the rows the paper prints as 1.000."""
    d = (t1 if df is None else df)
    e = 1.0 - (1.0 - d.eps_MB) * d.H_MB_over_H_0
    n_out = int(((e < 0.35) | (e > 0.65)).sum())
    ident = float(np.max(np.abs(e[NOEXP] - d.eps_MB[NOEXP])))
    return e.min(), e.max(), n_out, ident

rows = []
for label, kw in [("baseline", dict()),
                  ("eq. (3) exponent read as d^1 not d^2", dict(exponent=1.0)),
                  ("eq. (3) exponent read as d^3",         dict(exponent=3.0)),
                  ("rho_s used in place of (rho_s-rho_f)", dict(dens_col="rho_s"))]:
    med, lo, hi, sp = k0_span(**kw)
    rows.append((label, med, lo, hi, sp))
B3a = pd.DataFrame(rows, columns=["injected defect", "K0_eff/8e-4 median",
                                  "min", "max", "spread factor"])
display(B3a.round(3))
print("Note which column has the power. The SPREAD barely moves under a wrong\n"
      "exponent (6.35 -> 5.69), so quoting a scatter alone would not catch it;\n"
      "the MEDIAN moves by orders of magnitude and does. The rho_s / drho row is\n"
      "the same rho_f blind spot as before: it moves nothing, and nothing on\n"
      "this page claims it would.")

t1_bad_eps  = t1.copy(); t1_bad_eps.loc[4, "eps_MB"] += 0.10          # 0.452 -> 0.552
t1_slip_eps = t1.copy(); t1_slip_eps.loc[4, "eps_MB"] = 0.0452        # decimal slip
t1_bad_H    = t1.copy(); t1_bad_H.loc[NOEXP, "H_MB_over_H_0"] = 1.05
t1_bad_H2   = t1.copy(); t1_bad_H2.loc[8, "H_MB_over_H_0"] = 4.3      # 1.43 -> 4.3
rows = []
for label, df in [("baseline", None),
                  ("one eps_MB out by +0.10 (0.452 -> 0.552)", t1_bad_eps),
                  ("one eps_MB decimal slip (0.452 -> 0.0452)", t1_slip_eps),
                  ("the two H/H0 = 1.000 rows read as 1.05",   t1_bad_H),
                  ("one H/H0 decimal slip (1.43 -> 4.3)",      t1_bad_H2)]:
    lo, hi, nout, ident = eps0_span(df)
    rows.append((label, lo, hi, nout, ident))
B3b = pd.DataFrame(rows, columns=["injected defect", "eps_0 min", "eps_0 max",
                                  "rows outside 0.35-0.65",
                                  "|eps_0 - eps_MB| where H/H0 printed 1.000"])
display(B3b.round(4))
print("The blind spot is row 2: a 0.10 error in one eps_MB lands INSIDE the\n"
      "range the other rows already span, and the identity column cannot see it\n"
      "either because that row does not print H/H0 = 1.000. So this check\n"
      "catches decimal slips and gross transcription errors, and nothing finer.\n"
      "It is a transcription check. It is not evidence about the physics.")

# --- and the cross-page identity
sol_ok = next(iter(BEDS.values()))["sol"]           # the 100 um group B sand
bad = davidson_cloud(P["d_B"], sol_ok["u_f"], u_br_davies_taylor(P["d_B"] / 2))
display(Markdown(
    f"**The cross-page identity, broken.** Building $u_{{br}}$ from the bubble "
    f"*radius* instead of its diameter — the commonest slip in that formula — "
    f"changes $u_{{br}}/u_f$ from {sol_ok['u_br']/sol_ok['u_f']:.2f} to "
    f"{bad['u_br']/bad['u_f']:.2f} and $R_c/R$ from {sol_ok['RcR']:.4f} to "
    f"{bad['RcR']:.4f}, and moves the eq. (8) constant by "
    f"{100*(eq8_constant(u_br=u_br_davies_taylor(P['d_B']/2))-C8)/C8:+.0f} %. "
    f"The 0.55 % agreement between $1/\\sqrt2$ and 0.711 is therefore a "
    f"statement about the coefficient and not an accident of a formula that "
    f"would agree anyway."))'''))

cells.append(md(r"""### What these checks cannot see

Stated as flatly as possible, because the failure this repository keeps finding
is a check presented as evidence for something it has no power over.

1. **Nothing here tests the C/A boundary.** Geldart prints no expression for it.
   His Figure 3 shows a shaded band PQ "drawn so as to separate empirically open
   and half-closed points", and every powder in Table 1 is far from it. The page
   plots no C/A line at all.
2. **Nothing here tests the B/D boundary against a measurement.** Table 1 stops
   at 318 µm and 1.5 g cm$^{-3}$; at that density equation (8) does not begin
   until 817 µm, computed in the boundary-gap cell above. The B/D work on this
   page is algebra — recomputing $10^6$, and
   identifying equation (7) with the cloud condition — and the cloud picture is
   a *consequence* of two models agreeing, not evidence that either is right.
   Geldart himself calls the criterion "tentative".
3. **Nothing here depends on $\rho_f$**, so nothing here checks it. Break table 1
   shows $\rho_f = 0$ moving no constant at all.
4. **The A/B agreement is weak in two separate ways, and both are quantified
   above rather than asserted here.** *First, the classes are lopsided*: almost
   every row of Table 1 is measured group A, so a null model that predicts A for
   everything scores nearly as well as the criterion, and the whole-table score
   is worth one row over it. The number to read is the score on the three rows
   the criterion calls group B. *Second, it is not fully independent of the data
   it is tested on*: the 100 inside the 225 is equation (2), fitted to this
   table's own $U_{MB}$ column. What is independent is the $U_0$ side — equation
   (3) was fitted elsewhere, and section 3 shows it is the sloppy half. Read the
   whole-table score as *"a boundary built from a fit to this table's $U_{MB}$
   and an outside correlation for $U_0$ sorts a table that is almost entirely
   group A"*, not as a blind prediction.
   **And the margin depends on a tie-break**: two rows print $U_{MB} = U_0$
   exactly, and reading equation (4) as it is displayed rather than as section
   4.3 applies it moves every score. Break table 2 reports both.
5. **The pymrm cloud number is not evidence about Geldart.** It is E1.2's
   operator reproducing E1.2's closed form; the deviation quoted is
   discretisation error. What carries content is the *qualitative* result — a
   cloud on one side of the line, no dividing streamline on the other — and the
   0.55 % coefficient comparison, which is arithmetic on printed constants.
6. **Equation (2)'s scatter is a fit residual, not a test.** It is reported below
   for completeness and labelled every time it appears."""))

cells.append(code(r'''# equation (2) against the U_MB column: a FIT RESIDUAL, reported as such.
grpA = t1.group_meas == "A"
res2 = 100 * (U_MB_eq2(t1.d_sv_um) - t1.U_MB) / t1.U_MB
display(Markdown(
    f"Equation (2), $U_{{MB}} = 100\\,d_{{sv}}$, against the $U_{{MB}}$ column "
    f"it was drawn through: mean absolute deviation "
    f"**{res2[grpA].abs().mean():.1f} %** over the {int(grpA.sum())} rows "
    f"measured group A, worst **{res2[grpA].abs().max():.1f} %** "
    f"(the {t1.d_sv_um[res2[grpA].abs().idxmax()]:.0f} µm "
    f"{t1.powder[res2[grpA].abs().idxmax()].lower()} fraction). Over all 22 rows "
    f"it is {res2.abs().mean():.1f} % mean, {res2.abs().max():.1f} % worst. "
    f"**This is the residual of a line to the points it was fitted through.** "
    f"It says the relation is simple and the scatter is real; it is not "
    f"evidence that equation (2) predicts anything, and it is not counted "
    f"anywhere on this page as agreement."))'''))

cells.append(code(r'''metrics = {
    "eq6_constant_recomputed":        C6,
    "eq6_constant_dev_pct":           d6,
    "eq8_constant_recomputed":        C8,
    "eq8_constant_dev_pct":           d8,
    "worked_U0_dev_pct":              float(V1["dev_%"].iloc[0]),
    "worked_ratio_dev_pct":           float(V1["dev_%"].iloc[2]),
    "worked_sand_U0_dev_pct":         float(V1["dev_%"].iloc[3]),
    # --- the A/B test.  The whole-table score is NOT the headline: it must be
    # read against the null baseline, and the informative score is the one on
    # the rows the criterion calls B.  The tie-break variants are exported too,
    # because they are what a reimplementation from the displayed eq. (4) gets.
    "AB_rows_agreeing":               n_agree,
    "AB_rows_compared":               n_cmp,
    "AB_rows_null_baseline":          n_null,
    "AB_rows_agreeing_diakon":        a_dk,
    "AB_rows_compared_diakon":        n_dk,
    "AB_rows_null_baseline_diakon":   d_null,
    "AB_discriminating_rows":         n_disc,
    "AB_discriminating_rows_agreeing": a_disc,
    "AB_tiebreak_is_strict":          1,   # 1 = ratio > 1 (section 4.3); 0 = >= 1
    "AB_rows_agreeing_nonstrict":     ge_all,
    "AB_rows_agreeing_diakon_nonstrict":       ge_dk,
    "AB_discriminating_rows_agreeing_nonstrict": ge_disc,
    "AB_rows_null_baseline_nonstrict": ge_null,
    "AB_rows_in_225_to_229_gap":      len(in_gap),
    "BD_rows_in_1e6_gap":             len(in_gap8),
    "BD_max_prod8_in_table1":         float(t1.prod8.max()),
    "BD_eq8_onset_um_at_rho_s_1p5":   d8_onset_1p5,
    "eq3_K0eff_ratio_min":            float(k.min()),
    "eq3_K0eff_ratio_max":            float(k.max()),
    "eq3_K0eff_ratio_min_diakon":     float(kd.min()),
    "eq3_K0eff_ratio_max_diakon":     float(kd.max()),
    "eps0_inferred_min":              float(t1.eps_0_inferred.min()),
    "eps0_inferred_max":              float(t1.eps_0_inferred.max()),
    "eps0_identity_residual":         resid_noexp,
    "eq2_fit_residual_mean_abs_pct":  float(res2[grpA].abs().mean()),
    "eq8_davies_taylor_vs_geldart_pct": float(100 * (C8_dt - C8) / C8),
    "cloud_pymrm_vs_closed_form":     dev_cloud,
}
report_agreement("A1.7", metrics)'''))

# ---------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**To the classification itself: nothing, and the page should say so plainly.**
Geldart's criteria are two inequalities in two variables. They need no solver,
no discretisation and no library; sections 1 to 4 above are arithmetic on
pandas columns and would run identically without pymrm installed. A page that
wrapped them in a PDE to justify the import would be dishonest about what the
1973 paper is.

What the reimplementation adds is four things, none of which is a solver.

**The boundary constants are Geldart's arithmetic rounded, and the rounding is
harmless.** 225 should be 229.4 and $10^6$ should be $1.016\times10^6$, low by
1.9 % and 1.6 %. No powder in his own Table 1 falls in the gap between 225 and
229.4, and none falls in the gap between $10^6$ and $1.016\times10^6$ either —
the largest $(\rho_s-\rho_f)(d')^2$ in the table is $1.19\times10^5$, a factor
8.4 below both. Both gaps are computed above rather than assumed, so the
roundings change no verdict there — but the recomputation is what licenses that
statement, and it also fixes the constants for anyone re-deriving the boundary
at a different $\mu$ or $K_{MB}$, which is exactly what section 5.1 of the paper
asks for when it discusses raised pressure.

**The A/B criterion sorted against a measurement, row by row — and scored
honestly.** Geldart's own verdict on this is one sentence of prose about Diakon.
Doing it arithmetically across all 21 fractions that carry a measured $U_0$ gives
20 agreements, but that figure is close to free: 19 of the 21 rows are measured
group A, so a null model that predicts group A for everything scores 19, and the
criterion is worth exactly one row over it. The three rows on which the two
predictors can differ are the Diakon fractions at 220, 263 and 318 µm, and the
criterion gets **two of those three** right. That is the number this page
reports. Two further things a bare score hides are made explicit: the margin
rests on a **tie-break** — two rows print $U_{MB} = U_0$ exactly, and equation
(4) is displayed as $\ge$ but applied in the text as $>$, so the strict reading
is used and the alternative is costed in break table 2 — and the boundary is
**sharp**, the last group A row sitting 0.5 % inside it and measured group A.

**And the one genuine conflict resolved the right way round.** On the 220 µm
fraction equation (6) says group B, the strict ratio test on the measured
velocities says group A, and section 4.3 says group B — so **Geldart's own prose
sides with equation (6)**, against a mechanical reading of his velocity columns.
The criterion does not disagree with its author; the strict ratio test does, on
the single row where the author had already said not to read it strictly.

**The weak half of the boundary identified and quantified.** The A/B line is
$U_{MB}/U_0$, and only $U_0$ is predicted by something not fitted to this table.
Inverting equation (3) row by row shows the constant it would need spans a factor
of 6.4 across Geldart's own powders — from 0.46 to 2.94 times the printed
$8\times10^{-4}$ — even though its median is 0.98 times that value. It holds to
about ±21 % on Diakon above 100 µm and is much worse on the fine catalyst
fractions. That is the number a
reader needs before carrying the boundary to a new powder, and it is not in the
paper. It is also the reason `A1.6` (Wen and Yu's minimum fluidization
correlation) is a separate case: replacing equation (3) is the obvious next move,
and this page does not make it.

**The B/D boundary re-read as a cloud condition.** Equation (7) is the statement
that a 25 cm bubble rises more slowly than the interstitial gas, and that is the
condition under which the Davidson bubble has no cloud —
[`E1.2`](../E1.2-davidson-bubble/)'s $(R_c/R)^3 = (u_{br}+2u_f)/(u_{br}-u_f)$ has
a pole exactly there. Geldart's $1/\sqrt2$ and Davies–Taylor's 0.711 differ by
0.55 %, less than the 1.6 % by which Geldart rounded his own answer to $10^6$,
so the two are the same boundary. The *physics* of that equivalence is standard
and is written in words in both documents; what is added here is the arithmetic
that closes it, and the re-reading it licenses: the classification diagram as
**a contour map of the smallest bubble a powder can wrap a cloud around**, with
group D as the region where 25 cm is not enough. This is the one place a solver
earns its keep, because the picture of the same bubble on either side of the
line comes from a computed flow field rather than from the formula being
examined — and it is therefore also the one part of the page that inherits
`E1.2`: if that page's operator or closed form changes, these cloud numbers
move.

**Honest accounting of the reverse direction.** Nothing on this page improves
Geldart's classification, and nothing on it validates the C/A or B/D boundary
against data, because the paper contains no data on either."""))

# ------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

The whole classification is four short functions, and they are worth lifting
directly:

```python
g, mu, K_MB, C3, rho_f = 981.0, 1.8e-4, 100.0, 8.0e-4, 1.2e-3   # CGS, air ambient

# Geldart (1973) group of a powder.  d_um in um, rho_s in g/cm3.
# Returns 'A', 'B' or 'D'.  Group C is NOT returned: Geldart prints no
# expression for the C/A boundary, only a hand-drawn band on his Fig. 3.
def geldart_group(d_um, rho_s, eps_0=0.4, d_B=25.0):
    drho = rho_s - rho_f
    if drho * d_um**2 >= np.sqrt(g * d_B / 2) * mu * eps_0 / (C3 * g * 1e-8):
        return "D"
    return "A" if drho * d_um <= K_MB * mu / (C3 * g * 1e-4) else "B"

def U_0(d_um, rho_s):                       # eq. (3), Davies & Richardson
    return C3 * g * (d_um * 1e-4)**2 * (rho_s - rho_f) / mu

def U_MB(d_um):                             # eq. (2)
    return K_MB * d_um * 1e-4

# Smallest bubble diameter (cm) that carries a Davidson cloud.
def min_cloud_bubble(d_um, rho_s, eps_0=0.4):
    return (U_0(d_um, rho_s) / (0.711 * eps_0))**2 / g
```

Note that `geldart_group` uses the constants **recomputed here**
($K_{MB}\mu/(C_3 g\,10^{-4})$ and the $\sqrt{gd_B/2}$ expression) rather than
the printed 225 and $10^6$, so it stays correct when $\mu$, $K_{MB}$,
$\varepsilon_0$ or $d_B$ are changed — which is the point of writing it out.
Substitute the printed constants if you want Geldart's diagram exactly; on his
own Table 1 the two give identical answers.

**Where this connects in the gallery.**

- [`E1.2`](../E1.2-davidson-bubble/) — the Davidson bubble. Owns the percolation
  solve reused here, the closed form for $R_c/R$, and the no-cloud threshold that
  section 5 identifies with Geldart's equation (7).
- [`E2.1`](../E2.1-kunii-levenspiel-bubbling-bed/) — the bubbling-bed model.
  Presumes a bubble phase with a cloud, i.e. presumes the powder is on the
  cloud-bearing side of the line this page draws.
- `A1.6` — Wen and Yu's minimum fluidization correlation. The obvious
  replacement for equation (3), and therefore for the weakest input to the A/B
  boundary. **Not built here.**
- [`A1.1`](../A1.1-ergun-pressure-drop/) — the Ergun equation, the pressure-drop
  relation whose laminar limit equation (3) is a special case of.
- `A1.8` — the fluidisation regime map, the next case downstream of this one.

**What to change first if you reuse this.** The $8\times10^{-4}$ in equation (3)
hides a voidage and is the single input most likely to be wrong for your powder;
section 3 shows what it costs. The $\varepsilon_0 = 0.4$ and $d_B = 25$ cm in the
B/D boundary are Geldart's own "arbitrary though reasonable" choices, and break
table 1 shows how far each moves the line."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.update({
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
})
nbf.write(nb, "index.ipynb")
print("wrote index.ipynb")
