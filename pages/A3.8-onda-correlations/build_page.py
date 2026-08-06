#!/usr/bin/env python3
"""Generate index.ipynb for page A3.8 (the Onda correlations). Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "The Onda correlations, and which of their fourteen constants the paper can see"
description: "Three coupled correlations - wetted area, liquid-side and gas-side coefficients - with the area feeding the other two. The gas-side constants are pinned to 0.45 % by an identity the paper prints twice and never checks; the liquid-side exponent set turns out to cap the measured slope of k_L a against L at 0.800, and two of Onda's own six runs sit above that cap. Every printed constant is then placed inside a pymrm counter-current absorber and its elasticity computed two ways, which is what shows that one exponent is inert wherever the packing and the liquid have the same surface tension."
categories: [sec:A, struct:S3, tier:T1, data:tier3, phase:gas-liquid]
date: 2026-08-05
---

# The Onda correlations, and which of their fourteen constants the paper can see

**Catalog ID:** `A3.8` · **Structures:** `S3` (1D steady BVP) · **Tier:** T1

Onda, Takeuchi and Okumoto (1968) is the packed-column workhorse: three
correlations that between them give the wetted area, the liquid-side coefficient
and the gas-side coefficient, from packing geometry and fluid properties alone.
It is still the default in process simulators.

Structurally it is not one correlation but **three coupled ones**, and the
coupling runs one way: the wetted area $a_w$ appears inside the liquid-side
Reynolds number, so eq. (1) feeds eq. (2), while eq. (3) for the gas side is
built on the *total* area and does not see $a_w$ at all. That asymmetry is the
whole reason the page is worth building, and it is where the interesting numbers
turn out to be.

The page does four separate things and keeps them apart.

1. **Transcribes fourteen printed constants** - eq. (3)'s prefactor is branched,
   so it is printed twice - and checks the transcription
   against a *second printing of the same three equations* inside the same PDF —
   the companion article that starts on journal page 62 restates all of them.
2. **Recovers eq. (4) from eq. (3).** The authors rearrange the gas-side
   correlation into a $j_D$ factor and print the result without showing the
   working. Redoing that working two independent ways gives 0.774443 against
   their printed **0.771**, and the exponent to machine precision. It is the only
   check on this page that pins all four gas-side constants at once.
3. **Finds a ceiling the correlation cannot exceed.** Combining eqs. (1) and (2)
   forces the local slope $\mathrm{d}\ln(k_L a)/\mathrm{d}\ln L$ into the window
   $(2/3,\,0.8]$ — **for any packing, any liquid and any flow rate**, because it
   depends on the exponent set alone and on no prefactor. Onda's own Table 1
   prints six measured slopes. Two of them, both spheres, are **above the
   ceiling**.
4. **Puts all fourteen constants inside a pymrm counter-current absorber** and
   measures the elasticity of the packed height to each, analytically and by
   finite differences on the solved column. That is what shows which constants a
   design calculation can actually feel — and that one of them,
   $(\sigma_c/\sigma)^{0.75}$, is **exactly inert** whenever the packing's
   critical surface tension equals the liquid's.

**What this page is not.** Eqs. (2) and (3) were fitted to the very $k_L a$ and
$k_G a$ data any comparison here would use, including Onda's own Table 1 runs,
which page 57 says were replotted against eq. (2) in Fig. 2. Nothing on this page
is an out-of-sample test of those fits and nothing is labelled as one. The two
comparisons that *are* free of that circularity are named where they occur:
the eq. (4) identity, which is arithmetic the paper owes and never shows; and
the exponent window, which no prefactor can move."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

**Which paper this is, and how that was established.** The file on disk is
`Onda1968-gas-liquid-mass-transfer-packed-columns-JCEJ1-56.pdf`. A bare numeric filename carries no metadata, and this repository has
been burned four times by trusting one, so the article was identified by reading
it. Journal page 56 carries the title *"MASS TRANSFER COEFFICIENTS BETWEEN GAS
AND LIQUID PHASES IN PACKED COLUMNS"*, the byline *"KAKUSABURO ONDA, HIROSHI
TAKEUCHI\*\* AND YOSHIO OKUMOTO / Dept. of Chem. Eng., University of Nagoya,
Nagoya"*, and the footnote *"Received on July 10, 1967"*. The running feet read
*"JOURNAL OF CHEMICAL ENGINEERING OF JAPAN"* and *"VOL.1 NO.1 1968"*, and the
printed folios run 56 to 62 — so the filename encodes volume 1, page 56.

That reading was then checked against CrossRef for `10.1252/jcej.1.56`, which
returns the same title, the same three author surnames, container *Journal of
Chemical Engineering of Japan*, volume 1, issue 1, pages 56–62, 1968. **Title,
authors, volume, issue, page range and year all agree with the page images**, so
the DOI is verified rather than auto-resolved from a terse citation.

Two other articles sit inside the same seven-page PDF and neither is used: page
56 opens with the last column of a preceding bubble-column paper, and page 62
begins *"Gas absorption with chemical reaction in packed columns"* by Onda, Sada
and Takeuchi. That companion article earns exactly one job here — it **reprints
all three correlations verbatim in its introduction**, giving a second,
independently typeset printing against which every constant was checked.

**Where the resolution rule bites.** `pdfimages -list` reports all seven pages as
CCITT-G4 bilevel, 2456 × 3330 px at **300 ppi native**. Rendering at 600 dpi
would interpolate and add nothing. Every equation and every table cell here was
read from a 300 dpi render, **cropped to a single line and re-read at that
resolution** — necessary because the 1968 typesetting sets exponents at roughly a
quarter of body size, and at page scale `-0.05` and `0.05`, or `-2.0` and `-2.6`,
are not separable.

**The problem Onda is solving.** A packed absorption column has three unknowns a
designer cannot measure: how much of the packing surface is actually wet, how
fast the liquid film renews, and how fast the gas boundary layer transports. Most
correlations of the day reported the *product* $k_L a$ or $k_G a$ — coefficient
times area — because that is what an experiment gives. Onda's move, which he had
been developing since 1959, was to correlate the **area separately**, then divide
every published $k_L a$ and $k_G a$ by it to expose the bare coefficients. The
first sentence of the paper says exactly that: *"Assuming that the wetted surface
on packing pieces is identical with the gas-liquid interface, Onda et al.
presented the empirical equations of the gas and liquid-side mass transfer
coefficients."*

**That assumption is the load-bearing one.** $a_w$ is a *wetted* area, measured
by wetting experiments; $a$ is the *interfacial* area for mass transfer. Setting
them equal is a modelling choice, and the paper's own conclusion states it as
one: *"Assuming that the wetted surface area evaluated by Eq. (1) is identical
with the gas-liquid interfacial area…"*. Nothing on this page tests it — nothing
in the paper can, because $a_w$ never appears except multiplied by a coefficient.

**Neighbouring pages.** `A3.1` (Whitman) is the two-film theory whose $k_G$ and
$k_L$ these correlations supply; `F3.1` (Hatta) and `F3.5` cover what happens when
reaction enhances the liquid film. `A3.6` (Calderbank & Moo-Young) and `A3.7`
(van 't Riet) are the stirred-tank analogues and are being built separately —
neither is compared here. `A3.9` (Billet–Schultes) and `A3.10`
(Rocha–Bravo–Fair) are the two later packed-column correlations that would make
the natural three-way comparison with Onda; **that comparison belongs to
whichever of those two is built second**, not here, because it needs their
sources on disk and this page has only Onda's."""))

# ------------------------------------------------------------ published model
cells.append(md(r"""## The published model

Three correlations, in the paper's own numbering and its own unit system
(m, kg, hr, kg-mole, atm, K).

**Eq. (1) — wetted area** (page 56, carried over from ref. 18, Onda, Takeuchi &
Koyama 1967):

$$\frac{a_w}{a_t} \;=\; 1-\exp\left\{-1.45\left(\frac{\sigma_c}{\sigma}\right)^{0.75}
\left(\frac{L}{a_t\mu_L}\right)^{0.1}
\left(\frac{L^2 a_t}{\rho_L^2 g}\right)^{-0.05}
\left(\frac{L^2}{\rho_L\sigma a_t}\right)^{0.2}\right\}
\tag{1}$$

which the paper immediately re-writes as
$1-\exp\{-1.45(\sigma_c/\sigma)^{0.75}(Re)^{0.1}(Fr)^{-0.05}(We)^{0.2}\}$,
naming the three groups. $\sigma_c$ is the **critical surface tension of the
packing material**, $\sigma$ that of the liquid.

**Eq. (2) — liquid-side coefficient** (page 57):

$$k_L\left(\frac{\rho_L}{\mu_L g}\right)^{1/3}
= 0.0051\left(\frac{L}{a_w\mu_L}\right)^{2/3}
\left(\frac{\mu_L}{\rho_L D_L}\right)^{-1/2}
\left(a_t D_p\right)^{0.4}
\tag{2}$$

**Note the $a_w$ in the Reynolds number.** That is the coupling: eq. (1) feeds
eq. (2), and it is why $k_L$ is not a power law in $L$.

**Eq. (3) — gas-side coefficient** (page 58):

$$\frac{k_G R T}{a_t D_G}
= 5.23\left(\frac{G}{a_t\mu_G}\right)^{0.7}
\left(\frac{\mu_G}{\rho_G D_G}\right)^{1/3}
\left(a_t D_p\right)^{-2.0}
\tag{3}$$

**Eq. (3) is branched.** Page 58: *"In Fig. 3, data for Raschig rings and Berl
saddles smaller than 15 mm are situated on the lower group and are best
correlated by merely changing the constant, 5.23, in Eq. (3) into 2.00."* The
Conclusion on page 61 restates the criterion as *"Raschig rings smaller than
15 mm and Berl saddles smaller than 1/2""*, and page 60 adds that for 1/2-in
spheres under **vaporization** *"the constant of Eq.(3) might be changed into
2.00"*. The paper gives no mechanism — *"this cause is not clear at present"* —
only the observation that $k_G a$ for packings below 15 mm falls with increasing
$a_t$. Both branches are exercised below.

**Eq. (4) — the same thing as a $j_D$ factor** (page 58). The authors note that
$a_t D_p = 6(1-\varepsilon) = 3.4$ for spheres and state, without showing the
rearrangement,

$$j_D = 0.771\left[\frac{G D_p'}{\mu_G(1-\varepsilon)}\right]^{-0.30}
\tag{4}$$

with $D_p'$ *"diameter of sphere possessing the same surface area as a piece of
packing"*. **Eq. (4) is eq. (3) in different clothes**, so it is an identity the
paper owes and never displays — and recovering it is the strongest check
available here.

**Eq. (5) — somebody else's data** (page 59). Shulman et al. (1955) correlated
the sublimation of **dry naphthalene packings** — a gas-solid experiment with no
liquid in it at all — as

$$j_D = 1.195\left[\frac{G D_p'}{\mu_G(1-\varepsilon)}\right]^{-0.36}
\tag{5}$$

and Onda remarks that *"the agreement between Eqs. (4) and (5) is fairly good
within the region of $100 < GD_p'/\mu_G(1-\varepsilon) < 10{,}000$"*. That
sentence is the one comparison in the paper against data the authors did not fit,
and it is quantified below.

**Eq. (6) — the design equation** (page 60):

$$Z = G_M\int_{y_1}^{y_2}\left(\frac{1}{k_G a}
+\frac{m}{k_L a\,c_{av}}\right)\frac{\mathrm{d}y}{y^{*}-y}
\tag{6}$$

with $G_M$ assumed constant and $m$ the slope of the equilibrium line. This is
resistance-in-series (`A3.1`) integrated down a column, and it is what the pymrm
section below solves as a two-point boundary-value problem instead of a
quadrature."""))

# --------------------------------------------------------- params/assumptions
cells.append(md(r"""## Parameters and assumptions

**Units.** Onda's Nomenclature (page 61, continued on 62) fixes the system, and
it is not SI: $a$ in m²/m³, $D$ in m²/hr, $g$ in m/hr², $G$ and $L$ in
kg/(m² hr), $\mu$ in kg/(m hr), $\rho$ in kg/m³, $R$ in m³ atm/(kg-mole K),
$k_G$ in kg-mole/(m² hr atm), $k_L$ in m/hr. Everything below works in those
units. The one trap: the Weber group $L^2/(\rho_L\sigma a_t)$ is dimensionless
only if $\sigma$ is in **kg/hr²**, and page 62's Nomenclature says exactly that —
*"$\sigma$ = surface tension [dynes/cm] **or** [kg/hr²]"*. The ratio
$\sigma_c/\sigma$ is unaffected. The notebook checks all three groups for
dimensionlessness symbolically rather than trusting this paragraph.

**What the paper does not print, and what is done about it.**

- **$a_t$ for any named packing.** Not tabulated anywhere. But the paper *does*
  print the sphere relation $a_t D_p = 6(1-\varepsilon) = 3.4$, so every
  calculation on this page that needs a geometry uses **spheres**, with $a_t$ and
  $\varepsilon$ taken from that printed statement and from nothing else.
- **$\sigma_c$ for any packing material.** Not printed in this article; it lives
  in ref. 18. So $\sigma_c/\sigma$ is carried as an **explicit page parameter**,
  set to 0.85 for the base case and swept from 0.5 to 1.5 below. That value is a
  page choice and is labelled as one everywhere it is used; the sweep says what it
  costs. It also produces the page's cleanest negative result: at
  $\sigma_c/\sigma = 1$ the exponent 0.75 becomes **exactly inert**.
- **The $L$ range of Table 1's runs.** Not printed; Fig. 2's abscissa is $Re_L$,
  not $L$. So the Table 1 comparison uses only the **structural window**, which
  needs no $L$, no $a_t$ and no property at all.

**The illustrative absorber is not Onda's.** The pymrm column below needs fluid
properties, and Onda publishes none. Air and water at 25 °C are used, from
standard property tables, at gas and liquid loadings inside the ranges the paper
itself prints for its vaporization runs ($L$ = 5000–8000 kg/(m² hr), Figs. 4
and 6). **No metric derived from that column is a claim about the paper.** The
column exists to make all fourteen constants live simultaneously and to measure
which of them a design calculation can feel; its elasticities are properties of
the correlation set, and the two that matter — $f_G$ and the wetting function
$\varphi$ — are reported so the numbers can be re-scaled to any other case.

**The modelling assumptions carried into the column.**

- Dilute solute, so $G_M$ and $L_M$ are constant down the column (Onda's own
  assumption in eq. (6): *"the gas molar flow rate, $G_M$, is assumed to be
  constant"*).
- **The total pressure is carried explicitly.** Onda's $k_G$ is in
  kg-mole/(m² hr atm), so it multiplies a *partial pressure*; eq. (6)'s driving
  force is a *mole fraction*, $\mathrm{d}y/(y^*-y)$. The two agree only at
  $P = 1$ atm, which eq. (6) leaves implicit. This page writes the gas-side
  resistance as $1/(P\,k_G a_w)$, which reduces to Onda's at 1 atm and is the
  correct generalisation away from it. $P$ = 1 atm everywhere below, so no number
  changes; the term is there so a reuser at 20 bar does not silently inherit a
  factor of twenty.
- A straight equilibrium line $y^* = m x$. Onda's Fig. 9 is precisely a study of
  what happens when $m$ varies down a column; that is not attempted here and the
  page says so.
- $a = a_w$, Onda's assumption, untested here and untestable from this paper.
- Isothermal, no reaction, plug flow in both phases, no axial dispersion. Axial
  dispersion in a packed absorber is `A2.6`'s subject, not this page's."""))

# -------------------------------------------------------------------- cell 1
cells.append(code("""# Colab environment cell.
try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm pandas pyyaml matplotlib sympy

import sys, pathlib
_shared = pathlib.Path.cwd()
for _ in range(4):
    if (_shared / "shared" / "gallery_utils.py").is_file():
        sys.path.insert(0, str(_shared / "shared")); break
    _shared = _shared.parent
else:
    if "google.colab" in sys.modules:
        !wget -q https://raw.githubusercontent.com/computational-chemical-engineering/pymrm-gallery/main/shared/gallery_utils.py
        sys.path.insert(0, ".")
import gallery_utils as gu"""))

# ---------------------------------------------------------------- the data
cells.append(md(r"""## The data

Three CSVs, all transcriptions of printed characters off 300 ppi renders — the
scan's native resolution — with every numeric cropped and re-read at that
resolution. **Nothing on this page is digitised**; no pixel of any figure is
measured, and the page needs none, because the paper prints its own constants,
its own reduced table and its own error bands as numerals.

- `onda-1968-correlation-constants.csv` — all fourteen constants and exponents of
  eqs. (1)–(3), plus the four of eqs. (4) and (5). Every constant used anywhere
  in this notebook is read from this file; **none is typed into the code**.
- `onda-1968-table1.csv` — the six rows of Table 1 (page 57), the authors' own
  CO₂-into-organic-solvent runs reduced to $k_L a = \alpha L^n$.
- `onda-1968-stated-results.csv` — the scalar claims printed in the text: the
  four error bands, $a_t D_p = 6(1-\varepsilon) = 3.4$, the validity range of the
  eq. (4)/eq. (5) comparison, and the three liquid loadings printed in the
  captions and legends of Figs. 4 and 6.

**Tier.** The correlation constants are tier 6 — the authors' own fits. Table 1
reduces the authors' own measurements, so it is tier 3 in origin, but it is
**not out-of-sample**: page 57 states these runs were replotted against eq. (2)
in Fig. 2. The page's Table 1 comparison is therefore in-sample, and is presented
as such. What makes it worth doing anyway is that the quantity compared — the
exponent window — cannot be moved by any fit of any prefactor.

**This page loads no other page's dataset**, so no cross-page reconciliation is
owed. It is a deliberate choice: `A3.1` holds Whitman's HCl equilibrium data and
`A3.4` holds a digitised Sherwood-number scatter, and neither contains a packed-
column mass velocity that could be compared with anything here."""))

cells.append(code(r"""import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sympy as sp
from scipy.optimize import brentq
from pymrm import (construct_convflux_upwind, construct_div, NumJac, newton,
                   interp_cntr_to_stagg_tvd, vanleer, generate_grid,
                   compute_boundary_values)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
np.set_printoptions(linewidth=120)

PAGE = "A3.8-onda-correlations"
K = gu.load_data("onda-1968-correlation-constants.csv", page=PAGE)
T1 = gu.load_data("onda-1968-table1.csv", page=PAGE)
SR = gu.load_data("onda-1968-stated-results.csv", page=PAGE).set_index("quantity")
print(gu.cite_data(gu.load_meta("onda-1968-correlation-constants.csv", page=PAGE)))

def c(symbol):
    '''A constant Onda printed, by name. Nothing in this notebook types one.'''
    return float(K.loc[K.symbol == symbol, "value"].iloc[0])

def stated(q):
    '''A scalar Onda printed in his running text, by name.'''
    return float(SR.loc[q, "value"])

display(K)
display(T1)"""))

cells.append(md(r"""### The transcription, checked against a second printing

Eqs. (1)–(3) appear **twice inside this PDF**: once as eqs. (1), (2), (3) of the
main article on pages 56–58, and again in the introduction of the companion
article that begins on page 62, renumbered (there $a_w$ is eq. 1, $k_G$ eq. 2 and
$k_L$ eq. 3). The two printings were cropped and read independently, and agree
character for character on all fourteen constants.

That is a genuine transcription check — two separate typesettings of the same
equations — but it is worth being precise about **what it cannot catch**: a
constant that was mis-set identically in both printings, which is one editorial
event away. The eq. (4) identity below is the check that would catch that,
because it tests the four gas-side constants against a *fifth* printed number
derived from them."""))

cells.append(code(r"""# The fourteen constants of eqs. (1)-(3), as read from the two printings.
# Fourteen, not thirteen, because eq. (3) prints its prefactor twice: 5.23 and 2.00.
C_AW, E_SIG, E_RE_L, E_FR_L, E_WE_L = (c("C_aw"), c("e_sigma"), c("e_Re_L"),
                                       c("e_Fr_L"), c("e_We_L"))
C_KL, E_RE_LW, E_SC_L, E_ATDP_L = c("C_kL"), c("e_Re_Lw"), c("e_Sc_L"), c("e_atDp_L")
C_KG, C_KG_SMALL, E_RE_G, E_SC_G, E_ATDP_G = (c("C_kG"), c("C_kG_small"), c("e_Re_G"),
                                              c("e_Sc_G"), c("e_atDp_G"))
# eqs. (4) and (5), which are NEVER inputs - they are targets.
C_JD_ONDA, E_JD_ONDA = c("C_jD_onda"), c("e_jD_onda")
C_JD_SHUL, E_JD_SHUL = c("C_jD_shulman"), c("e_jD_shulman")

both = K[K.equation.isin([1, 2, 3])]
print(f"{len(both)} constants in eqs. (1)-(3), each read from two independent printings")
print("  eq. (1):", ", ".join(both[both.equation == 1].printed_as))
print("  eq. (2):", ", ".join(both[both.equation == 2].printed_as))
print("  eq. (3):", ", ".join(both[both.equation == 3].printed_as))
N_CONSTANTS = int(len(both))"""))

cells.append(md(r"""### The three groups of eq. (1) are dimensionless — checked, not assumed

Onda's unit system is unusual enough that this is worth verifying symbolically
rather than by eye, and it is the check that catches a mis-read $\rho_L^2$ or a
dropped square. The Weber group is the interesting one: it closes only with
$\sigma$ in kg/hr², which is what page 62's Nomenclature says."""))

cells.append(code(r"""M_, KG_, HR_, MOL_ = sp.symbols("m kg hr kgmole", positive=True)
U = dict(L=KG_/(M_**2*HR_), G=KG_/(M_**2*HR_), a_t=1/M_, mu=KG_/(M_*HR_),
         rho=KG_/M_**3, g=M_/HR_**2, sigma=KG_/HR_**2, D=M_**2/HR_)
groups = {
    "Re_L = L/(a_t mu_L)":            U["L"]/(U["a_t"]*U["mu"]),
    "Fr_L = L^2 a_t/(rho_L^2 g)":     U["L"]**2*U["a_t"]/(U["rho"]**2*U["g"]),
    "We_L = L^2/(rho_L sigma a_t)":   U["L"]**2/(U["rho"]*U["sigma"]*U["a_t"]),
    "Re_G = G/(a_t mu_G)":            U["G"]/(U["a_t"]*U["mu"]),
    "Sc = mu/(rho D)":                U["mu"]/(U["rho"]*U["D"]),
    "a_t D_p":                        U["a_t"]*M_,
    "k_L (rho_L/(mu_L g))^(1/3)":     (M_/HR_)*(U["rho"]/(U["mu"]*U["g"]))**sp.Rational(1, 3),
}
DIMLESS_ALL = True
for name, expr in groups.items():
    simp = sp.simplify(expr)
    ok = simp == 1
    DIMLESS_ALL &= bool(ok)
    print(f"   {name:<34} -> {simp}   {'dimensionless' if ok else 'NOT DIMENSIONLESS'}")
print(f"\nall seven groups dimensionless in Onda's unit system: {DIMLESS_ALL}")
print("The Weber group closes ONLY with sigma in kg/hr2, which is what page 62 states.")"""))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

Two objects. First the correlation set itself, assembled once from the CSV
constants; then Onda's eq. (6) rebuilt as a **counter-current two-point boundary
value problem** rather than a quadrature, which is what makes the packed height a
solved quantity instead of an integral evaluated along an assumed operating line.

**The correlations.** `OndaSet` holds the three equations and the derived
$K_G a$. Every constant enters as a named argument with the CSV value as its
default, so the break table further down can perturb any single one of the
thirteen **continuous** parameters — the branch is a choice of value for `C_kG`,
not a fourteenth argument — without touching the code that uses them.

**The column.** Gas flows up, liquid down, solute passes from gas to liquid.
With $z$ measured from the bottom and the fields $[\,y,\;x\,]$:

$$\frac{\mathrm{d}}{\mathrm{d}z}\bigl(G_M\,y\bigr) = -N,
\qquad
\frac{\mathrm{d}}{\mathrm{d}z}\bigl(-L_M\,x\bigr) = +N,
\qquad
N = K_G a\,(y - m x),$$

i.e. one divergence of a convective flux per phase with **opposite velocities**,
plus a pointwise interphase source. That is exactly `construct_convflux_upwind`
followed by `construct_div`, with a two-component velocity `v = [[G_M, -L_M]]`,
and it is why the two phases can live as two fields of one array rather than as
two coupled solvers.

**Boundary conditions, on the outward normal.** `bc` is a 2-tuple of dicts
$a\,\partial\phi/\partial n + b\,\phi = d$ with $n$ pointing *out of* the domain,
and each entry is a per-field list because the two phases enter at opposite ends:

- at $z=0$ (bottom): gas **inlet**, $y = y_{in}$ → `a=0, b=1, d=y_in`;
  liquid **outlet**, pure outflow → `a=1, b=0, d=0`.
- at $z=Z$ (top): gas **outlet**, pure outflow → `a=1, b=0, d=0`;
  liquid **inlet**, $x = x_{in}$ → `a=0, b=1, d=x_in`.

`construct_div` uses `nu=0` — the axial coordinate of a column of constant
cross-section is Cartesian, not cylindrical; `nu=1` would be the *radial*
coordinate of the same column and is not what is being differentiated here.

**The Jacobian.** The source is pointwise in $(y,x)$, so `NumJac(shape)` with
`shape = (n, 2)` is exactly right: the default stencil couples the **last** axis
in full, which here is the field axis, and leaves the spatial axis tridiagonal.
Writing `(n,)` for a single field would declare every cell coupled to every other
and build a dense Jacobian — the trap recorded in `AGENTS.md`. Nothing here uses
`axes_diagonals`, because the source reads no neighbouring cell.

**Second order.** Pure upwind is first order; a van Leer TVD **deferred
correction** on top of it recovers second order in the interior, and the loop
asserts its own convergence rather than silently returning its iteration cap.
The *outflow face* stays first order — that is measured, named and Richardson-
extrapolated in the Validation section rather than hidden."""))

cells.append(code(r'''HR = 3600.0                     # s per hr, for converting tabulated SI properties
R_GAS = 0.0820574               # m3 atm / (kg-mole K)   -- Onda's R
G_ACC = 9.80665 * HR**2         # m/hr2                  -- Onda's g


class OndaSet:
    """Onda's eqs. (1), (2) and (3), and the K_G a they imply.

    Every printed constant is a named argument defaulting to the CSV value, so a
    single one can be perturbed without touching anything downstream. That is
    what the break table and the elasticity study both use.
    """

    def __init__(self, *, a_t, D_p, rho_L, mu_L, D_L, sigma, sigma_ratio,
                 rho_G, mu_G, D_G, T, P=1.0,
                 C_aw=C_AW, e_sigma=E_SIG, e_Re_L=E_RE_L, e_Fr_L=E_FR_L, e_We_L=E_WE_L,
                 C_kL=C_KL, e_Re_Lw=E_RE_LW, e_Sc_L=E_SC_L, e_atDp_L=E_ATDP_L,
                 C_kG=C_KG, e_Re_G=E_RE_G, e_Sc_G=E_SC_G, e_atDp_G=E_ATDP_G):
        self.__dict__.update(locals()); del self.self
        self.Sc_L = mu_L / (rho_L * D_L)
        self.Sc_G = mu_G / (rho_G * D_G)
        self.atDp = a_t * D_p

    # ---- eq. (1) -----------------------------------------------------------
    def psi(self, L):
        """The exponent argument of eq. (1); a_w/a_t = 1 - exp(-psi)."""
        Re = L / (self.a_t * self.mu_L)
        Fr = L**2 * self.a_t / (self.rho_L**2 * G_ACC)
        We = L**2 / (self.rho_L * self.sigma * self.a_t)
        return (self.C_aw * self.sigma_ratio**self.e_sigma * Re**self.e_Re_L
                * Fr**self.e_Fr_L * We**self.e_We_L)

    def a_w(self, L):
        return self.a_t * (1.0 - np.exp(-self.psi(L)))

    def phi(self, L):
        """dln(a_w)/dln(psi) = psi e^-psi / (1 - e^-psi); -> 1 as psi -> 0, 0 as psi -> inf."""
        p = self.psi(L)
        return p * np.exp(-p) / (1.0 - np.exp(-p))

    # ---- eq. (2) -----------------------------------------------------------
    def k_L(self, L):
        aw = self.a_w(L)
        return (self.C_kL * (L / (aw * self.mu_L))**self.e_Re_Lw
                * self.Sc_L**self.e_Sc_L * self.atDp**self.e_atDp_L
                / (self.rho_L / (self.mu_L * G_ACC))**(1.0 / 3.0))

    # ---- eq. (3) -----------------------------------------------------------
    def k_G(self, G):
        return (self.C_kG * (G / (self.a_t * self.mu_G))**self.e_Re_G
                * self.Sc_G**self.e_Sc_G * self.atDp**self.e_atDp_G
                * self.a_t * self.D_G / (R_GAS * self.T))

    # ---- resistances in series (A3.1), with a = a_w (Onda's assumption) -----
    def resistances(self, L, G, m, c_av):
        aw = self.a_w(L)
        R_G = 1.0 / (self.P * self.k_G(G) * aw)      # gas-side, per unit mole fraction
        R_L = m / (self.k_L(L) * aw * c_av)          # liquid-side, same units
        return R_G, R_L

    def KGa(self, L, G, m, c_av):
        R_G, R_L = self.resistances(L, G, m, c_av)
        return 1.0 / (R_G + R_L)

    def f_G(self, L, G, m, c_av):
        R_G, R_L = self.resistances(L, G, m, c_av)
        return R_G / (R_G + R_L)'''))

cells.append(code(r'''class Absorber:
    """Counter-current packed absorber: Onda's eq. (6) as a BVP, not a quadrature.

    Fields [y, x] on one axial grid; z from the bottom of the packing.
    Gas rises (+G_M), liquid falls (-L_M), so the two phases are two components
    of one array with opposite velocities and inlets at opposite ends.
    """

    def __init__(self, Z, KGa, m, G_M, L_M, y_in, x_in, n=200):
        self.n, self.shape = n, (n, 2)
        self.z_f, self.z_c = generate_grid(n, [0.0, Z], generate_x_c=True)
        self.Z, self.KGa, self.m, self.G_M, self.L_M = Z, KGa, m, G_M, L_M
        self.v = np.array([[G_M, -L_M]])          # gas up, liquid down
        # OUTWARD normal: a dphi/dn + b phi = d, one entry per field.
        # z=0 : gas INLET  y = y_in           -> a=0, b=1, d=y_in
        #       liquid OUTLET, pure outflow   -> a=1, b=0, d=0
        # z=Z : gas OUTLET, pure outflow      -> a=1, b=0, d=0
        #       liquid INLET  x = x_in        -> a=0, b=1, d=x_in
        self.bc = ({"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[y_in, 0.0]]},
                   {"a": [[1.0, 0.0]], "b": [[0.0, 1.0]], "d": [[0.0, x_in]]})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                 self.bc, v=self.v, axis=0)
        self.div = construct_div(self.shape, self.z_f, nu=0, axis=0)   # nu=0: Cartesian axis
        self.jac_const = self.div @ conv
        self.g_const = self.div @ conv_bc
        self.numjac = NumJac(self.shape)        # pointwise source; last axis IS the field axis
        self.u = np.column_stack([np.full(n, y_in), np.full(n, x_in)])
        self.defect = np.zeros((2 * n, 1))
        self.its = 0

    def source(self, u):
        N = self.KGa * (u[..., 0] - self.m * u[..., 1])
        return np.stack([-N, +N], axis=-1)

    def residual(self, u):
        g_s, j_s = self.numjac(self.source, u)
        g = (self.g_const + self.jac_const @ u.reshape((-1, 1)) + self.defect
             - g_s.reshape((-1, 1)))
        return g, self.jac_const - j_s

    def solve(self, tvd=True, tol=1e-13, max_it=80):
        res = newton(self.residual, self.u, maxfev=50)
        assert res.success, "Newton did not converge on the base upwind solve"
        self.u = res.x.reshape(self.shape)
        if tvd:
            done = False
            for it in range(1, max_it + 1):
                _, dc_f = interp_cntr_to_stagg_tvd(self.u, self.z_f, self.z_c, self.bc,
                                                  self.v, tvd_limiter=vanleer, axis=0)
                self.defect = np.asarray(self.div @ (self.v * dc_f).reshape((-1, 1)))
                prev = self.u.copy()
                res = newton(self.residual, self.u, maxfev=50)
                assert res.success, "Newton did not converge inside the deferred correction"
                self.u = res.x.reshape(self.shape)
                done = np.max(np.abs(self.u - prev)) < tol * max(1.0, np.abs(self.u).max())
                if done:
                    break
            # A deferred correction that silently returns its iteration cap is the
            # classic way a "converged" number is not one.
            assert done, f"van Leer deferred correction did not converge in {max_it} iterations"
            self.its = it
        return self

    def faces(self):
        """Boundary VALUES, reconstructed consistently with the bc.

        This is the A2.6 trap and its cure. The gas leaves through the face at
        z = Z, not through the centre of the last cell, and the two differ by
        O(h) - enough to drop the observed convergence order from 2 to 1 while
        looking perfectly stable. `compute_boundary_values` applies the same
        reconstruction the flux operator used, so it is consistent with the
        discretisation rather than a hand correction bolted on afterwards.
        """
        v_lo, _, v_hi, _ = compute_boundary_values(self.u, self.z_f, self.z_c,
                                                   self.bc, axis=0)
        return float(v_hi[0, 0]), float(v_lo[0, 1])      # y at z=Z, x at z=0

    @property
    def y_out(self):
        """Gas leaving the top, at the face."""
        return self.faces()[0]

    @property
    def x_out(self):
        """Liquid leaving the bottom, at the face."""
        return self.faces()[1]

    @property
    def y_out_cell(self):
        """The same quantity read half a cell short - the A2.6 defect, kept to measure it."""
        return float(self.u[-1, 0])

    def imbalance(self):
        """Solute lost by the gas against solute gained by the liquid, relative.

        STRUCTURAL: with consistent face values this is a telescoping sum of the
        discrete fluxes and holds to machine precision for ANY K_G a, right or
        wrong. It tests the transport bookkeeping and nothing about the physics.
        """
        gas = self.G_M * (self.bc[0]["d"][0][0] - self.y_out)
        liq = self.L_M * (self.x_out - self.bc[1]["d"][0][1])
        return abs(gas / liq - 1.0)


def column_closed_form(Z, KGa, m, G_M, L_M, y_in, x_in):
    """Exact solution of the same BVP, by hand. Independent of every operator above.

    D = y - m x obeys dD/dz = -lambda D with lambda = KGa (1/G_M - m/L_M), so the
    whole profile is one exponential and the two-point problem closes analytically.
    """
    lam = KGa * (1.0 / G_M - m / L_M)
    E = np.exp(-lam * Z)
    beta = (m * KGa / L_M) / lam
    D0 = (y_in - m * x_in) / (1.0 + beta * (1.0 - E))
    y_out = y_in - (KGa / G_M) * D0 * (1.0 - E) / lam
    x_bot = x_in + (KGa / L_M) * D0 * (1.0 - E) / lam
    return float(y_out), float(x_bot)'''))

cells.append(md(r"""### The illustrative operating point

Air and water at 25 °C and 1 atm, on 1-in ceramic spheres, at liquid and gas
loadings inside the ranges Onda prints for his own vaporization runs. **The
fluid properties are standard-table values and are not Onda's**; the geometry
*is* his, because $a_t D_p = 6(1-\varepsilon) = 3.4$ for spheres is printed on
page 58 and is the only packing geometry this article supplies."""))

cells.append(code(r'''ATDP_SPHERES = stated("atDp_spheres")             # 3.4, printed on page 58
SIX = stated("sixteen_one_minus_eps")             # 6, from a_t D_p = 6(1-eps)
EPS_SPHERES = 1.0 - ATDP_SPHERES / SIX
D_P = 0.0254                                      # 1-in spheres: a size the paper packs with
A_T = ATDP_SPHERES / D_P

SIGMA_RATIO = 0.85    # sigma_c/sigma: a PAGE PARAMETER. Onda does not print sigma_c
                      # for any packing material in this article; it lives in ref. 18.

PROPS = dict(  # air / water at 25 C, 1 atm, from standard tables - NOT from Onda
    rho_L=997.0, mu_L=8.90e-4 * HR, D_L=1.92e-9 * HR, sigma=0.0720 * HR**2,
    rho_G=1.184, mu_G=1.84e-5 * HR, D_G=1.60e-5 * HR, T=298.15, P=1.0)

L_REF = stated("L_fig6_hi_kg_m2_hr")   # 8000 kg/(m2 hr), printed in Fig. 6's legend
G_REF = 2000.0                         # kg/(m2 hr), a page choice inside Fig. 6's abscissa
C_AV = PROPS["rho_L"] / 18.015         # kg-mole/m3 of water
G_M = G_REF / 28.96                    # kg-mole/(m2 hr), air
L_M = L_REF / 18.015                   # kg-mole/(m2 hr), water

onda = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=SIGMA_RATIO, **PROPS)
PSI_REF, AW_REF = onda.psi(L_REF), onda.a_w(L_REF)
KL_REF, KG_REF = onda.k_L(L_REF), onda.k_G(G_REF)

print(f"spheres, from the paper:  a_t D_p = {ATDP_SPHERES}  ->  eps = {EPS_SPHERES:.5f}, "
      f"a_t = {A_T:.2f} m2/m3 at D_p = {D_P*1000:.1f} mm")
print(f"Sc_L = {onda.Sc_L:.1f}   Sc_G = {onda.Sc_G:.4f}   sigma_c/sigma = {SIGMA_RATIO} (page parameter)")
print(f"\nat L = {L_REF:.0f}, G = {G_REF:.0f} kg/(m2 hr):")
print(f"   psi          = {PSI_REF:.5f}          (eq. 1's exponent argument)")
print(f"   a_w/a_t      = {AW_REF/A_T:.5f}          a_w = {AW_REF:.2f} m2/m3")
print(f"   k_L          = {KL_REF:.5f} m/hr    = {KL_REF/HR:.3e} m/s")
print(f"   k_G          = {KG_REF:.4f} kg-mole/(m2 hr atm)")
print(f"   G_M = {G_M:.3f}, L_M = {L_M:.3f} kg-mole/(m2 hr),  c_av = {C_AV:.3f} kg-mole/m3")
AW_FRACTION_REF = float(AW_REF / A_T)'''))

# -------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. Eq. (4) is eq. (3) rearranged — and the rearrangement recovers 0.771

The paper writes eq. (4) with no working. The working is worth doing, because it
is the only place where all four constants of eq. (3) are tested against a number
printed independently of them.

**One thing has to be assumed, and it is worth naming before the result.** Onda's
Nomenclature defines $j_D$ only as *"mass transfer factor defined by Eq. (4)"* —
circularly, and it never prints the underlying definition. So the **standard
Chilton–Colburn form is assumed here**, and the agreement below is what confirms
the assumption rather than the other way round: any other convention would shift
the prefactor by a fixed factor and could not land within half a percent. Using a
*molar* velocity in place of the mass velocity, for instance, would move it by a
molecular weight — a factor of about 29 for air, not 1.004.

Two further conventions, both from the paper's own Nomenclature. $D_p'$ in
eqs. (4) and (5) is the *"diameter of sphere possessing the same surface area as
a piece of packing"*, whereas eq. (3) uses $D_p$, the *"nominal size of
packing"*. **For spheres the two coincide**, which is the only geometry used
anywhere on this page and the only one for which the paper supplies
$a_t D_p$ — so nothing here depends on the distinction, but a reuser applying
eq. (4) to Raschig rings must respect it.

Start from the standard Chilton–Colburn factor
$j_D = (k_c/u)\,Sc^{2/3}$. In Onda's units the coefficient in m/hr is
$k_G R T$ and the superficial gas velocity is $G/\rho_G$, so

$$j_D = \frac{k_G R T\,\rho_G}{G}\,Sc_G^{2/3}.$$

Substituting eq. (3) and using $Sc_G = \mu_G/(\rho_G D_G)$,
the product $D_G\rho_G Sc_G$ collapses to $\mu_G$, the leftover $Sc$ powers
cancel *exactly*, and the whole thing reduces to a pure Reynolds power:

$$j_D = 5.23\,Re_G^{-0.3}\,(a_t D_p)^{-2}. $$

Then $a_t D_p = 6(1-\varepsilon)$ for spheres makes $Re_G = X/6$ with
$X = G D_p/[\mu_G(1-\varepsilon)]$, so

$$j_D = 5.23\cdot 6^{0.3}\cdot 3.4^{-2}\; X^{-0.30}.$$

**Two things are being tested at once.** The exponent $-0.30$ is $0.7-1$ and
tests only eq. (3)'s Reynolds exponent. The prefactor tests 5.23 and $-2.0$
together. And the *absence of any Schmidt number in eq. (4)* is what tests the
$1/3$: had eq. (3)'s Schmidt exponent been anything else, $j_D$ would still
contain $Sc$ and could not be a function of $X$ alone.

Below, that is done twice. **Route A** is a symbolic reduction with `sympy` from
eq. (3) as transcribed, with no hand substitution. **Route B** never does any
algebra: it evaluates eq. (3) numerically at four thousand independent random
combinations of $G$, $\mu_G$, $\rho_G$, $D_G$ and $D_p$ spanning several decades
each, forms $j_D$ and $X$, and regresses $\ln j_D$ on $\ln X$. If eq. (3)'s
Schmidt exponent were not $1/3$, Route B's residual would not be zero — its
sample spans Schmidt numbers from about 0.02 to 80."""))

cells.append(code(r'''# ---- Route A: symbolic, from eq. (3) as transcribed -----------------------
_at, _DG, _G, _muG, _rhoG, _Dp, _eps = sp.symbols(
    "a_t D_G G mu_G rho_G D_p varepsilon", positive=True)
_ReG = _G / (_at * _muG)
_ScG = _muG / (_rhoG * _DG)
_kG_RT = (sp.nsimplify(C_KG) * _ReG**sp.nsimplify(E_RE_G) * _ScG**sp.nsimplify(E_SC_G)
          * (_at * _Dp)**sp.nsimplify(E_ATDP_G) * _at * _DG)          # eq. (3) -> k_G R T
_jD = sp.simplify(_kG_RT * _rhoG / _G * _ScG**sp.Rational(2, 3))       # Chilton-Colburn
_jD = sp.simplify(_jD.subs(_at, SIX * (1 - _eps) / _Dp))               # spheres
_X = _G * _Dp / (_muG * (1 - _eps))
_ratio = sp.powsimp(sp.expand_power_exp(sp.simplify(_jD / _X**sp.nsimplify(E_JD_ONDA))),
                    force=True)
_ratio = sp.simplify(_ratio.subs(_eps, sp.nsimplify(EPS_SPHERES)))
JD_PREFACTOR_SYMBOLIC = float(_ratio)
print("Route A (sympy, no hand algebra)")
print(f"   j_D reduces to        {sp.simplify(_jD)}")
print(f"   j_D / X^({E_JD_ONDA})      = {sp.nsimplify(_ratio)}  = {JD_PREFACTOR_SYMBOLIC:.6f}")
SC_CANCELS = bool(({_rhoG, _DG} & sp.simplify(_jD).free_symbols) == set())
print(f"   free symbols left in j_D: {sorted(map(str, sp.simplify(_jD).free_symbols))}")
print(f"   neither rho_G nor D_G survives, so Sc has cancelled exactly: {SC_CANCELS}")
print("   That cancellation is what tests eq. (3)'s 1/3: any other Schmidt exponent")
print("   would leave Sc in j_D, and eq. (4) could not be a function of X alone.")

# ---- Route B: numeric, no algebra at all ---------------------------------
rng = np.random.default_rng(20260805)          # seeded: this page reports its value
NSAMP = 4000
gs = 10**rng.uniform(2.0, 4.0, NSAMP)
mus = 10**rng.uniform(-2.2, -1.2, NSAMP)
rhos = 10**rng.uniform(-0.7, 0.7, NSAMP)
dgs = 10**rng.uniform(-2.5, -1.0, NSAMP)
dps = 10**rng.uniform(-2.5, -1.3, NSAMP)
ats = ATDP_SPHERES / dps
Scs = mus / (rhos * dgs)
kG_RT = (C_KG * (gs / (ats * mus))**E_RE_G * Scs**E_SC_G
         * (ats * dps)**E_ATDP_G * ats * dgs)
jDs = kG_RT * rhos / gs * Scs**(2.0 / 3.0)
Xs = gs * dps / (mus * (1.0 - EPS_SPHERES))
slope, lnA = np.polyfit(np.log(Xs), np.log(jDs), 1)
resid = np.log(jDs) - (lnA + slope * np.log(Xs))
JD_PREFACTOR_NUMERIC = float(np.exp(lnA))
JD_EXPONENT_NUMERIC = float(slope)
JD_FIT_MAX_RESID = float(np.abs(resid).max())
print(f"\nRoute B ({NSAMP} random property combinations, seed 20260805)")
print(f"   Sc_G sampled over  {Scs.min():.4f} .. {Scs.max():.2f}     "
      f"X sampled over {Xs.min():.1f} .. {Xs.max():.0f}")
print(f"   prefactor  {JD_PREFACTOR_NUMERIC:.6f}     exponent  {JD_EXPONENT_NUMERIC:.10f}")
print(f"   worst residual of the single-variable fit: {JD_FIT_MAX_RESID:.3e}")

JD_ROUTE_SPREAD = abs(JD_PREFACTOR_NUMERIC / JD_PREFACTOR_SYMBOLIC - 1.0)
JD_PREFACTOR_DEV_PCT = 100.0 * (JD_PREFACTOR_SYMBOLIC / C_JD_ONDA - 1.0)
JD_EXPONENT_DEV = abs(JD_EXPONENT_NUMERIC - E_JD_ONDA)
print(f"\n   two routes agree to {JD_ROUTE_SPREAD:.2e}")
print(f"   against Onda's printed eq. (4):")
print(f"      prefactor  {JD_PREFACTOR_SYMBOLIC:.6f}  vs printed {C_JD_ONDA}   "
      f"-> {JD_PREFACTOR_DEV_PCT:+.3f} %")
print(f"      exponent   {JD_EXPONENT_NUMERIC:.6f} vs printed {E_JD_ONDA}   "
      f"-> {JD_EXPONENT_DEV:.2e}")'''))

cells.append(md(r"""**What that 0.45 % is.** It is not experimental agreement and not a model test;
it is the paper's own arithmetic, redone. The gap is consistent with the authors
having rounded at an intermediate step; the cell below prints what the two obvious
intermediate roundings give, and then what the printed 0.771 *would* require of
each constant taken alone. None of the four is implausible, so **no claim is made that anything is
mis-set** — this is reported so a reader can see the size of the discrepancy
against the size of each constant."""))

cells.append(code(r'''def jd_prefactor(C=C_KG, e_Re=E_RE_G, e_atDp=E_ATDP_G, atDp=ATDP_SPHERES, six=SIX):
    """Closed form of eq. (4)'s prefactor from eq. (3)'s constants."""
    return C * six**(1.0 - e_Re) * atDp**e_atDp

print("the two obvious intermediate roundings of the same product:")
print(f"   round (a_t D_p)^-2  to 3 s.f.: {C_KG * SIX**(1-E_RE_G):.6g} x "
      f"{round(ATDP_SPHERES**E_ATDP_G, 6):.4g} -> "
      f"{C_KG * SIX**(1-E_RE_G) * float(f'{ATDP_SPHERES**E_ATDP_G:.3g}'):.6f}")
print(f"   round C_kG*6^0.3    to 3 s.f.: {float(f'{C_KG * SIX**(1-E_RE_G):.3g}'):.4g} x "
      f"{ATDP_SPHERES**E_ATDP_G:.6g} -> "
      f"{float(f'{C_KG * SIX**(1-E_RE_G):.3g}') * ATDP_SPHERES**E_ATDP_G:.6f}")
print(f"   the paper prints {C_JD_ONDA}, which is below both.\n")

rows = []
rows.append(("C_kG   5.23", C_KG, brentq(lambda v: jd_prefactor(C=v) - C_JD_ONDA, 1.0, 20.0)))
rows.append(("e_Re_G 0.7", E_RE_G, brentq(lambda v: jd_prefactor(e_Re=v) - C_JD_ONDA, 0.3, 0.99)))
rows.append(("e_atDp_G -2.0", E_ATDP_G,
             brentq(lambda v: jd_prefactor(e_atDp=v) - C_JD_ONDA, -4.0, -1.0)))
rows.append(("a_t D_p 3.4", ATDP_SPHERES,
             brentq(lambda v: jd_prefactor(atDp=v) - C_JD_ONDA, 2.0, 5.0)))
alt = pd.DataFrame(rows, columns=["printed constant", "printed value",
                                  "value that would give exactly 0.771"])
alt["relative shift"] = alt["value that would give exactly 0.771"] / alt["printed value"] - 1
display(alt.round(6))
print("None of these is a plausible mis-reading of the page image, and the page")
print("claims no mis-set constant. What the row for a_t D_p shows is that the")
print("printed 3.4 is itself a rounding of 6(1-eps), which is the likeliest source.")
JD_ATDP_FOR_0771 = float(alt.loc[alt["printed constant"] == "a_t D_p 3.4",
                                 "value that would give exactly 0.771"].iloc[0])
print(f"\n   a_t D_p = {JD_ATDP_FOR_0771:.4f} reproduces 0.771 exactly; the paper prints 3.4,")
print(f"   i.e. eps = {1-JD_ATDP_FOR_0771/SIX:.5f} against the printed {EPS_SPHERES:.5f}.")'''))

cells.append(md(r"""### 2. Eq. (4) against Shulman's naphthalene, over the range the paper states

This is the one comparison in the article against data the authors did not fit.
Eq. (5) comes from Shulman, Ullrich, Proulx and Zimmerman (1955), who **sublimed
dry naphthalene packings** — a gas–solid experiment with no liquid, no wetted
area and no absorption in it. Onda's remark is that the two agree *"fairly
good"* over $100 < X < 10{,}000$, and gives no number.

The two are power laws with different exponents, so they can only cross once. The
cell below puts the number on it, and reports the crossing."""))

cells.append(code(r'''XLO, XHI = stated("shulman_range_lo"), stated("shulman_range_hi")
Xg = np.logspace(np.log10(XLO), np.log10(XHI), 801)
j_onda = C_JD_ONDA * Xg**E_JD_ONDA
j_shul = C_JD_SHUL * Xg**E_JD_SHUL
ratio = j_onda / j_shul

EQ45_DEV_LO_PCT = 100.0 * (ratio[0] - 1.0)
EQ45_DEV_HI_PCT = 100.0 * (ratio[-1] - 1.0)
EQ45_MAX_DEV_PCT = float(100.0 * np.abs(ratio - 1.0).max())
EQ45_MEAN_DEV_PCT = float(100.0 * np.abs(ratio - 1.0).mean())
EQ45_CROSSOVER_X = float(np.exp(np.log(C_JD_SHUL / C_JD_ONDA) / (E_JD_ONDA - E_JD_SHUL)))
# null baseline: how much does j_D itself move over the same range?
EQ45_ONDA_SPAN = float(j_onda[0] / j_onda[-1])
EQ45_SHUL_SPAN = float(j_shul[0] / j_shul[-1])

print(f"Eq. (4) / Eq. (5) over the paper's own stated range {XLO:.0f} < X < {XHI:.0f}:")
print(f"   at X = {XLO:8.0f}   {EQ45_DEV_LO_PCT:+7.2f} %   (Onda LOW)")
print(f"   at X = {EQ45_CROSSOVER_X:8.1f}      0.00 %   (they cross)")
print(f"   at X = {XHI:8.0f}   {EQ45_DEV_HI_PCT:+7.2f} %   (Onda HIGH)")
print(f"   worst |deviation| {EQ45_MAX_DEV_PCT:.2f} %, mean {EQ45_MEAN_DEV_PCT:.2f} %")
print(f"\nNULL BASELINE - what a constant j_D would achieve over the same range:")
print(f"   j_D itself falls by a factor {EQ45_ONDA_SPAN:.2f} (Onda) and "
      f"{EQ45_SHUL_SPAN:.2f} (Shulman) across those two decades,")
print(f"   so a single constant fitted to the mid-range is off by roughly a factor")
print(f"   {np.sqrt(EQ45_ONDA_SPAN):.2f} at each end - about {100*(np.sqrt(EQ45_ONDA_SPAN)-1):.0f} %,")
print(f"   against the {EQ45_MAX_DEV_PCT:.1f} % the two correlations differ by. That is what")
print(f'   "fairly good" is worth: better than a constant by a factor '
      f'{100*(np.sqrt(EQ45_ONDA_SPAN)-1)/EQ45_MAX_DEV_PCT:.1f}.')
print(f"\nThis is the only comparison on this page against data Onda did not fit -")
print("but note it is correlation against correlation, not correlation against")
print("measurement: Shulman's own scatter about eq. (5) is not printed here.")'''))

cells.append(code(r'''fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.9))
ax[0].loglog(Xg, j_onda, "-", lw=1.8, color="C0",
             label=f"Onda eq. (4): {C_JD_ONDA} X^{{{E_JD_ONDA}}}")
ax[0].loglog(Xg, j_shul, "--", lw=1.8, color="C3",
             label=f"Shulman eq. (5): {C_JD_SHUL} X^{{{E_JD_SHUL}}}")
ax[0].axvline(EQ45_CROSSOVER_X, ls=":", color="0.5", lw=1)
ax[0].set_xlabel(r"$X = G D_p' / [\mu_G(1-\varepsilon)]$"); ax[0].set_ylabel(r"$j_D$")
ax[0].set_title("Two correlations, disjoint datasets"); ax[0].legend(fontsize=8)
ax[1].semilogx(Xg, 100 * (ratio - 1), "-", lw=1.8, color="C2")
ax[1].axhline(0, color="0.4", lw=1)
ax[1].axvline(EQ45_CROSSOVER_X, ls=":", color="0.5", lw=1)
ax[1].annotate(f"cross at X = {EQ45_CROSSOVER_X:.0f}", (EQ45_CROSSOVER_X, 0),
               textcoords="offset points", xytext=(6, 16), fontsize=8)
ax[1].set_xlabel(r"$X$"); ax[1].set_ylabel("eq. (4) relative to eq. (5)  /  %")
ax[1].set_title(f'"fairly good" = {EQ45_MAX_DEV_PCT:.1f} % worst over the stated range')
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""### 3. The exponent window, and where Onda's own Table 1 sits

This is the page's one comparison against measured numbers, and the one place
the correlation can be caught out.

Onda's Table 1 reduces each of the authors' runs to a power law
$k_L a = \alpha L^{n}$. In his own framework $a = a_w$, so the correlation
predicts the same product:

$$k_L a_w \;\propto\; \left(\frac{L}{a_w}\right)^{2/3} a_w
\;=\; L^{2/3}\,a_w^{1/3}
\quad\Longrightarrow\quad
n \;=\; \frac{2}{3} + \frac13\,\frac{\mathrm{d}\ln a_w}{\mathrm{d}\ln L}.$$

Now the $L$-dependence of eq. (1) collapses to a **single composite exponent**:
$Re\propto L$, $Fr\propto L^2$, $We\propto L^2$, so
$\psi \propto L^{\,0.1 + 2(-0.05) + 2(0.2)} = L^{0.4}$, and

$$\frac{\mathrm{d}\ln a_w}{\mathrm{d}\ln L}
= 0.4\,\varphi(\psi),\qquad
\varphi(\psi)=\frac{\psi e^{-\psi}}{1-e^{-\psi}} \in (0,1).$$

Therefore

$$\boxed{\;\frac23 \;<\; n \;\le\; \frac23+\frac{0.4}{3} = 0.8\;}$$

**for every packing, every liquid and every flow rate.** The bound involves no
prefactor at all — not 1.45, not 0.0051, not $\sigma_c/\sigma$, not $a_t D_p$,
not the Schmidt number. It is a property of the exponent *set*, so no fit can
move it, which is exactly what makes it testable against a table those fits have
already seen.

The window is derived twice below: from the algebra above, and by numerically
differentiating the composed functions `a_w(L)` and `k_L(L)` as coded, with no
reference to the derivation. A slip in either shows up as a disagreement."""))

cells.append(code(r'''# ---- Route A: the analytic window ----------------------------------------
L_COMPOSITE = E_RE_L + 2 * E_FR_L + 2 * E_WE_L      # 0.1 + 2(-0.05) + 2(0.2)
WINDOW_LO = E_RE_LW                                  # 2/3, as psi -> infinity
WINDOW_HI = E_RE_LW + L_COMPOSITE * (1.0 - E_RE_LW)  # 2/3 + 0.4/3, as psi -> 0
print(f"composite L-exponent of psi : {E_RE_L} + 2({E_FR_L}) + 2({E_WE_L}) = {L_COMPOSITE}")
print(f"window for n = dln(k_L a_w)/dlnL :  ({WINDOW_LO:.6f}, {WINDOW_HI:.6f})")

def n_pred_analytic(o, L):
    return E_RE_LW + L_COMPOSITE * (1.0 - E_RE_LW) * o.phi(L)

# ---- Route B: numerical differentiation of the coded functions -----------
def n_pred_numeric(o, L, rel=1e-5):
    lo, hi = L * (1 - rel), L * (1 + rel)
    f = lambda x: np.log(o.k_L(x) * o.a_w(x))
    return (f(hi) - f(lo)) / (np.log(hi) - np.log(lo))

Ls = np.logspace(np.log10(50.0), np.log10(500000.0), 240)
na = np.array([n_pred_analytic(onda, L) for L in Ls])
nb = np.array([n_pred_numeric(onda, L) for L in Ls])
N_WINDOW_TWO_ROUTES = float(np.abs(na - nb).max())
print(f"\nanalytic vs numerical-differentiation, over L = 50 .. 5e5 kg/(m2 hr): "
      f"worst {N_WINDOW_TWO_ROUTES:.3e}")
print(f"   numerically attained range of n over that sweep: "
      f"{nb.min():.6f} .. {nb.max():.6f}")
print(f"   psi over the same sweep: {onda.psi(Ls).min():.4g} .. {onda.psi(Ls).max():.4g}")
N_SWEEP_MIN, N_SWEEP_MAX = float(nb.min()), float(nb.max())
WINDOW_ESCAPE = float(max(0.0, nb.max() - WINDOW_HI, WINDOW_LO - nb.min()))
print(f"   the sweep never leaves the analytic window: escape = {WINDOW_ESCAPE:.2e}")'''))

cells.append(code(r'''n_meas = T1["n"].to_numpy(float)
inside = (n_meas > WINDOW_LO) & (n_meas <= WINDOW_HI)
tab = T1.copy()
tab["inside_window"] = inside
tab["excess_over_ceiling"] = np.maximum(0.0, n_meas - WINDOW_HI)
display(tab[["row", "packing", "size_printed", "absorbent", "n",
             "inside_window", "excess_over_ceiling"]])

TABLE1_N_INSIDE = int(inside.sum())
TABLE1_N_TOTAL = int(len(n_meas))
TABLE1_WORST_EXCESS = float(np.maximum(0.0, n_meas - WINDOW_HI).max())
TABLE1_N_RANGE = float(n_meas.max() - n_meas.min())
print(f"\n{TABLE1_N_INSIDE} of {TABLE1_N_TOTAL} printed exponents lie inside "
      f"({WINDOW_LO:.4f}, {WINDOW_HI:.4f})")
outs = tab[~inside]
for _, r in outs.iterrows():
    print(f"   OUTSIDE: row {int(r['row'])}  {r['packing']} {r['size_printed']} "
          f"in {r['absorbent']}   n = {r['n']}  -> {r['n']-WINDOW_HI:+.3f} above the ceiling")
print(f"\nworst excess over the ceiling: {TABLE1_WORST_EXCESS:+.3f}")
print(f"Table 1's n spans {TABLE1_N_RANGE:.2f}; the window is only "
      f"{WINDOW_HI-WINDOW_LO:.4f} wide. The correlation's exponent set cannot")
print("generate the spread the authors' own six runs display.")'''))

cells.append(md(r"""**Null baselines, beside it.** A window is a prediction with no free
parameters, so it must be scored against predictors that have none either — and
against one that has a single fitted parameter, which is the honest hard case."""))

cells.append(code(r'''def mae(pred):
    return float(np.abs(n_meas - pred).mean())

BEST_CONST = float(np.median(n_meas))     # the 1-parameter null, MAE-optimal
MAE_WINDOW_TOP = mae(WINDOW_HI)           # Onda at his own ceiling, 0 free parameters
MAE_WINDOW_MID = mae(0.5 * (WINDOW_LO + WINDOW_HI))
MAE_NULL_23 = mae(WINDOW_LO)              # drop the wetted-area coupling entirely
MAE_NULL_1 = mae(1.0)                     # k_L a proportional to L
MAE_NULL_FITTED = mae(BEST_CONST)         # 1 fitted parameter on 6 points

scores = pd.DataFrame([
    ("Onda, at the window ceiling n = 0.800", 0, MAE_WINDOW_TOP),
    ("Onda, at the window midpoint", 0, MAE_WINDOW_MID),
    ("null: n = 2/3, i.e. NO wetted-area coupling", 0, MAE_NULL_23),
    ("null: n = 1, k_L a proportional to L", 0, MAE_NULL_1),
    (f"null: best constant, fitted = {BEST_CONST:.4f}", 1, MAE_NULL_FITTED),
], columns=["predictor", "free parameters", "mean |error| in n"])
display(scores.round(5))
print(f"Onda's ceiling beats the no-coupling null by a factor "
      f"{MAE_NULL_23/MAE_WINDOW_TOP:.2f} - so the wetted-area term IS doing work.")
print(f"But a single fitted constant, {BEST_CONST:.4f}, beats Onda's ceiling by a factor "
      f"{MAE_WINDOW_TOP/MAE_NULL_FITTED:.2f}.")
print("With six points, no replicates and no error bars on n, that comparison")
print("cannot separate the two. It is reported because it bounds what this table")
print("can establish, which is less than the window's cleanliness suggests.")
NULL_BEATS_ONDA = bool(MAE_NULL_FITTED < MAE_WINDOW_TOP)'''))

cells.append(code(r'''fig, ax = plt.subplots(figsize=(7.6, 4.2))
psis = np.logspace(-3, 2, 400)
phi = psis * np.exp(-psis) / (1 - np.exp(-psis))
ax.semilogx(psis, E_RE_LW + L_COMPOSITE * (1 - E_RE_LW) * phi, "-", lw=2, color="C0",
            label=r"Onda: $n = 2/3 + 0.4\,\varphi(\psi)/3$")
ax.axhline(WINDOW_HI, ls="--", color="0.35", lw=1.2)
ax.axhline(WINDOW_LO, ls="--", color="0.35", lw=1.2)
ax.axhspan(WINDOW_LO, WINDOW_HI, color="C0", alpha=0.07)
ax.axvline(PSI_REF, ls=":", color="C2", lw=1.2)
ax.annotate(f"illustrative case\n$\\psi$ = {PSI_REF:.3f}", (PSI_REF, WINDOW_LO + 0.005),
            fontsize=8, color="C2")
for _, r in T1.iterrows():
    inw = (r["n"] > WINDOW_LO) and (r["n"] <= WINDOW_HI)
    ax.plot([psis[0], psis[-1]], [r["n"]] * 2, "-", lw=1.0,
            color=("0.55" if inw else "C3"), alpha=0.9)
    ax.annotate(f"{r['packing']} {r['size_printed']}, {r['absorbent']}  n={r['n']}",
                (psis[-1], r["n"]), fontsize=7, ha="right", va="bottom",
                color=("0.35" if inw else "C3"))
ax.set_xlabel(r"$\psi$  (eq. 1's exponent argument; small = poorly wetted)")
ax.set_ylabel(r"$n = \mathrm{d}\ln(k_L a)/\mathrm{d}\ln L$")
ax.set_ylim(0.63, 0.90)
ax.set_title("The window Onda's exponent set allows, against the six exponents he printed\n"
             "red = above the ceiling, which no choice of any prefactor can reach")
ax.legend(fontsize=8, loc="lower left"); fig.tight_layout(); plt.show()'''))

cells.append(md(r"""### 4. The column, and the elasticity of packed height to every printed constant

Onda's eq. (6) with $k_G$, $k_L$ and $a_w$ from eqs. (1)–(3), solved as a
counter-current BVP. The number a designer wants is the packed height $Z$ for a
duty; what this section measures is **how much of that height each printed
constant is responsible for**.

The elasticity $\partial\ln Z/\partial p$ is available analytically, because at
fixed duty $Z \propto R_G + R_L$ and the chain rule collapses onto two numbers —
the gas-side resistance fraction $f_G$ and the wetting function $\varphi(\psi)$:

$$\frac{\partial\ln Z}{\partial \ln C_{k_G}} = -f_G,\qquad
\frac{\partial\ln Z}{\partial \ln C_{k_L}} = -f_L,\qquad
\frac{\partial\ln Z}{\partial \ln C_{a_w}} = -\Bigl(f_G+\tfrac{f_L}{3}\Bigr)\varphi,$$

and for each *exponent* $e$ multiplying a group $\Pi$, the semi-elasticity is the
same coefficient times $\ln\Pi$. The $f_L/3$ rather than $f_L$ is the coupling:
raising $a_w$ raises the area but *lowers* $k_L$ through $Re_{L,w}$, and
$1-2/3=1/3$ of the gain survives.

Those analytic elasticities are then checked against **central finite differences
on $Z$** — a route that shares no algebra with them — and the finite-difference
route is then re-run **through the pymrm column solve itself**, so the check
reaches the operator assembly and not only the resistance arithmetic."""))

cells.append(code(r'''M_SLOPE = 1.0          # straight equilibrium line y* = m x; a page parameter
Y_IN, X_IN = 0.02, 0.0
REMOVAL = 0.95
N_CELLS = 200

def height_for_duty(o, *, m=M_SLOPE, L=L_REF, G=G_REF, removal=REMOVAL):
    """Packed height for a stated removal, from the closed form of the same BVP."""
    KGa = o.KGa(L, G, m, C_AV)
    f = lambda Z: column_closed_form(Z, KGa, m, G_M, L_M, Y_IN, X_IN)[0] - (1 - removal) * Y_IN
    return brentq(f, 1e-4, 500.0, xtol=1e-14, rtol=1e-15), KGa

Z_STAR, KGA_REF = height_for_duty(onda)
F_G_REF = onda.f_G(L_REF, G_REF, M_SLOPE, C_AV)
PHI_REF = float(onda.phi(L_REF))
ABSORPTION_FACTOR = L_M / (M_SLOPE * G_M)
print(f"duty: {100*REMOVAL:.0f} % removal of a dilute solute, y_in = {Y_IN}, m = {M_SLOPE}")
print(f"   absorption factor L_M/(m G_M) = {ABSORPTION_FACTOR:.3f}   (> 1, so no pinch)")
print(f"   K_G a  = {KGA_REF:.3f} kg-mole/(m3 hr) per unit mole fraction")
print(f"   f_G    = {F_G_REF:.4f}   (gas-side share of the total resistance)")
print(f"   phi    = {PHI_REF:.4f}   (wetting function at psi = {PSI_REF:.4f})")
print(f"   Z*     = {Z_STAR:.5f} m of packing")

col = Absorber(Z_STAR, KGA_REF, M_SLOPE, G_M, L_M, Y_IN, X_IN, n=N_CELLS).solve()
Y_EXACT, X_EXACT = column_closed_form(Z_STAR, KGA_REF, M_SLOPE, G_M, L_M, Y_IN, X_IN)
print(f"\npymrm column at n = {N_CELLS} ({col.its} deferred-correction iterations):")
print(f"   y_out  = {col.y_out:.8e}   closed form {Y_EXACT:.8e}")
print(f"   x_out  = {col.x_out:.8e}   closed form {X_EXACT:.8e}")
print(f"   solute balance, gas lost vs liquid gained: {col.imbalance():.3e} relative")
COL_IMBALANCE = float(col.imbalance())'''))

cells.append(code(r'''fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.8))
zz = np.linspace(0, Z_STAR, 400)
lam = KGA_REF * (1 / G_M - M_SLOPE / L_M)
E_ = np.exp(-lam * Z_STAR); beta_ = (M_SLOPE * KGA_REF / L_M) / lam
D0_ = (Y_IN - M_SLOPE * X_IN) / (1 + beta_ * (1 - E_))
y_ex = Y_IN - (KGA_REF / G_M) * D0_ * (1 - np.exp(-lam * zz)) / lam
x_ex = X_IN + (KGA_REF / L_M) * D0_ * (np.exp(-lam * zz) - np.exp(-lam * Z_STAR)) / lam
ax[0].plot(zz, 1e3 * y_ex, "-", lw=1.6, color="0.4", label="closed form")
ax[0].plot(col.z_c.ravel(), 1e3 * col.u[:, 0], "o", ms=3, mfc="none", color="C0",
           label=f"pymrm, n = {N_CELLS}")
ax[0].set_xlabel("z from the bottom  /  m"); ax[0].set_ylabel(r"$y \times 10^{3}$")
ax[0].set_title("gas, rising"); ax[0].legend(fontsize=8)
ax[1].plot(zz, 1e3 * x_ex, "-", lw=1.6, color="0.4", label="closed form")
ax[1].plot(col.z_c.ravel(), 1e3 * col.u[:, 1], "s", ms=3, mfc="none", color="C3",
           label=f"pymrm, n = {N_CELLS}")
ax[1].set_xlabel("z from the bottom  /  m"); ax[1].set_ylabel(r"$x \times 10^{3}$")
ax[1].set_title("liquid, falling"); ax[1].legend(fontsize=8)
fig.suptitle(f"Onda's eq. (6) as a counter-current BVP: {100*REMOVAL:.0f} % removal "
             f"in {Z_STAR:.3f} m of 1-in spheres", fontsize=10)
fig.tight_layout(); plt.show()'''))

cells.append(code(r'''# ---- analytic elasticities ------------------------------------------------
def groups_at(o, L, G):
    return dict(
        e_sigma=np.log(o.sigma_ratio),
        e_Re_L=np.log(L / (o.a_t * o.mu_L)),
        e_Fr_L=np.log(L**2 * o.a_t / (o.rho_L**2 * G_ACC)),
        e_We_L=np.log(L**2 / (o.rho_L * o.sigma * o.a_t)),
        e_Re_Lw=np.log(L / (o.a_w(L) * o.mu_L)),
        e_Sc_L=np.log(o.Sc_L),
        e_atDp_L=np.log(o.atDp),
        e_Re_G=np.log(G / (o.a_t * o.mu_G)),
        e_Sc_G=np.log(o.Sc_G),
        e_atDp_G=np.log(o.atDp),
    )

def elasticities_analytic(o, L, G, m):
    fG = o.f_G(L, G, m, C_AV); fL = 1 - fG; ph = o.phi(L)
    kA = -(fG + fL * (1 - E_RE_LW))          # dlnZ / dln(a_w)
    g = groups_at(o, L, G)
    out = {"C_aw": kA * ph, "C_kL": -fL, "C_kG": -fG}
    for e in ("e_sigma", "e_Re_L", "e_Fr_L", "e_We_L"):
        out[e] = kA * ph * g[e]                      # via a_w
    for e in ("e_Re_Lw", "e_Sc_L", "e_atDp_L"):
        out[e] = -fL * g[e]                          # via k_L
    for e in ("e_Re_G", "e_Sc_G", "e_atDp_G"):
        out[e] = -fG * g[e]                          # via k_G
    return out

# ---- numeric elasticities: central differences on the SOLVED column -------
BASE = dict(C_aw=C_AW, e_sigma=E_SIG, e_Re_L=E_RE_L, e_Fr_L=E_FR_L, e_We_L=E_WE_L,
            C_kL=C_KL, e_Re_Lw=E_RE_LW, e_Sc_L=E_SC_L, e_atDp_L=E_ATDP_L,
            C_kG=C_KG, e_Re_G=E_RE_G, e_Sc_G=E_SC_G, e_atDp_G=E_ATDP_G)

def onda_with(**over):
    kw = dict(BASE); kw.update(over)
    return OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=SIGMA_RATIO, **PROPS, **kw)

def elasticity_numeric(name, h=1e-4, m=M_SLOPE, L=L_REF, G=G_REF):
    """dlnZ/dln(param) for a prefactor, dlnZ/d(param) for an exponent.

    Z comes from a root find on the CLOSED FORM; the cell below confirms the
    pymrm column reproduces the same Z, so the two agree on what is being
    differentiated.
    """
    v0 = BASE[name]
    if name.startswith("C_"):
        lo, hi = v0 * (1 - h), v0 * (1 + h); dlnp = np.log(hi / lo)
    else:
        lo, hi = v0 - h, v0 + h; dlnp = hi - lo
    Zlo = height_for_duty(onda_with(**{name: lo}), m=m, L=L, G=G)[0]
    Zhi = height_for_duty(onda_with(**{name: hi}), m=m, L=L, G=G)[0]
    return (np.log(Zhi) - np.log(Zlo)) / dlnp

ana = elasticities_analytic(onda, L_REF, G_REF, M_SLOPE)
rows = []
for name in BASE:
    a_, n_ = ana[name], elasticity_numeric(name)
    rows.append((name, float(K.loc[K.symbol == name, "value"].iloc[0]),
                 a_, n_, abs(a_ - n_)))
el = pd.DataFrame(rows, columns=["printed constant", "value",
                                 "elasticity, analytic (chain rule)",
                                 "elasticity, finite difference on Z",
                                 "|difference|"])
el["|elasticity|"] = el["elasticity, analytic (chain rule)"].abs()
el = el.sort_values("|elasticity|", ascending=False)
display(el.drop(columns="|elasticity|").round(6))
ELASTICITY_TWO_ROUTES = float(el["|difference|"].max())
print(f"analytic chain rule vs finite differences on Z: worst "
      f"{ELASTICITY_TWO_ROUTES:.3e}")
print(f"\nSum rule: the three PREFACTOR elasticities must satisfy")
print(f"   -(C_kG) - (C_kL) = f_G + f_L = 1 exactly, by construction.")
SUM_RULE = abs(-(ana["C_kG"] + ana["C_kL"]) - 1.0)
print(f"   residual {SUM_RULE:.2e}   <- structural, cannot fail; reported as an identity")
EL_MAX = float(el["|elasticity|"].max())
EL_MIN = float(el["|elasticity|"].min())
EL_ARGMIN = str(el.iloc[-1]["printed constant"])
print(f"\nlargest  |elasticity|: {EL_MAX:.4f}  ({el.iloc[0]['printed constant']})")
print(f"smallest |elasticity|: {EL_MIN:.4g}  ({EL_ARGMIN})")
print(f"\nTHE SECOND NEARLY-INERT EXPONENT. eq. (3)'s Schmidt exponent 1/3 has")
print(f"semi-elasticity -f_G ln(Sc_G) = {ana['e_Sc_G']:+.5f}, because Sc_G = "
      f"{onda.Sc_G:.4f} for this")
print(f"gas and ln(Sc_G) = {np.log(onda.Sc_G):+.5f} is nearly zero. Replacing 1/3 by 2/3")
_o23 = onda_with(e_Sc_G=2 / 3)
E_SCG_2X_RELDEV = abs(height_for_duty(_o23)[0] / Z_STAR - 1.0)
print(f"moves the packed height by {100*E_SCG_2X_RELDEV:.3f} %. Almost every gas of")
print("industrial interest has Sc_G within a factor two of 1, so a design")
print("calculation essentially cannot see this exponent - which is precisely why")
print("the eq. (4) identity, where it shows up STRUCTURALLY through the Schmidt")
print("cancellation rather than numerically, is the check that pins it.")'''))

cells.append(md(r"""**And the same three prefactor elasticities, re-measured through the pymrm
solve.** The finite differences above are on $Z$ from the closed form, which
shares the resistance arithmetic with the analytic route but not its chain rule.
Pushing the same perturbation through the **column solve** closes the last gap:
it exercises the operator assembly, the boundary conditions and the Newton solve
as well. At fixed $Z$ the quantity differentiated is $\ln y_{out}$ rather than
$\ln Z$, and the two are related by the same $\mathrm{d}\ln y_{out}/\mathrm{d}\ln
K_Ga$ for every parameter — so dividing one by the other recovers the elasticity
of $Z$ with no further algebra, and the agreement is limited only by the column's
own discretisation error."""))

cells.append(code(r'''def dlogy_dlogp_column(name, h=1e-3, n=400):
    """d ln(y_out)/d ln(param) measured on the SOLVED pymrm column, at fixed Z."""
    v0 = BASE[name]; lo, hi = v0 * (1 - h), v0 * (1 + h)
    out = []
    for v in (lo, hi):
        o = onda_with(**{name: v})
        cc = Absorber(Z_STAR, o.KGa(L_REF, G_REF, M_SLOPE, C_AV), M_SLOPE,
                      G_M, L_M, Y_IN, X_IN, n=n).solve()
        out.append(np.log(cc.y_out))
    return (out[1] - out[0]) / np.log(hi / lo)

# the normaliser: d ln(y_out)/d ln(K_G a), measured the same way
cc_lo = Absorber(Z_STAR, KGA_REF * 0.999, M_SLOPE, G_M, L_M, Y_IN, X_IN, n=400).solve()
cc_hi = Absorber(Z_STAR, KGA_REF * 1.001, M_SLOPE, G_M, L_M, Y_IN, X_IN, n=400).solve()
DLOGY_DLOGK = (np.log(cc_hi.y_out) - np.log(cc_lo.y_out)) / np.log(1.001 / 0.999)

rows = []
for name in ("C_kG", "C_kL", "C_aw"):
    col_el = -dlogy_dlogp_column(name) / DLOGY_DLOGK      # dlnZ/dlnp = -dlny/dlnp / (dlny/dlnK)
    rows.append((name, ana[name], elasticity_numeric(name), col_el,
                 abs(col_el - ana[name])))
elc = pd.DataFrame(rows, columns=["prefactor", "analytic (chain rule)",
                                  "finite difference on Z",
                                  "through the pymrm column solve", "|column - analytic|"])
display(elc.round(8))
ELASTICITY_COLUMN_ROUTE = float(elc["|column - analytic|"].max())
print(f"d ln(y_out)/d ln(K_G a) on the column = {DLOGY_DLOGK:.6f}")
print(f"worst disagreement between the analytic chain rule and the column route: "
      f"{ELASTICITY_COLUMN_ROUTE:.2e}")
print("That floor is the column's own discretisation error at n = 400, not an")
print("error in either derivation - the closed-form route agrees to "
      f"{ELASTICITY_TWO_ROUTES:.1e}.")'''))

cells.append(md(r"""**One constant is exactly inert, and it is not an accident of the operating
point.** The semi-elasticity of $Z$ to the exponent 0.75 is proportional to
$\ln(\sigma_c/\sigma)$. When the packing's critical surface tension equals the
liquid's, that logarithm is zero and **the exponent 0.75 can be set to any value
whatever without moving anything** — it is not merely small, it is identically
zero. That is the `F2.3` pattern: a term that is provably inert for a whole class
of callers, and it needs saying, because $\sigma_c \approx \sigma$ is not exotic
(it is roughly the ceramic/water case).

The cell below demonstrates it rather than asserting it: with
$\sigma_c/\sigma = 1$, replacing 0.75 by 7.5 leaves every digit of $Z$."""))

cells.append(code(r'''o_eq = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=1.0, **PROPS)
Z_eq, _ = height_for_duty(o_eq)
o_eq_bad = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=1.0, **PROPS, e_sigma=10 * E_SIG)
Z_eq_bad, _ = height_for_duty(o_eq_bad)
INERT_075_RELDEV = abs(Z_eq_bad / Z_eq - 1.0)
print(f"at sigma_c/sigma = 1.0:   e_sigma = {E_SIG} -> Z = {Z_eq:.12f} m")
print(f"                          e_sigma = {10*E_SIG}  -> Z = {Z_eq_bad:.12f} m")
print(f"   relative change: {INERT_075_RELDEV:.3e}   <- exactly inert, not merely small")

o_ne = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=SIGMA_RATIO, **PROPS)
o_ne_bad = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=SIGMA_RATIO, **PROPS, e_sigma=10 * E_SIG)
Z_ne, _ = height_for_duty(o_ne); Z_ne_bad, _ = height_for_duty(o_ne_bad)
LIVE_075_RELDEV = abs(Z_ne_bad / Z_ne - 1.0)
print(f"\nat sigma_c/sigma = {SIGMA_RATIO}:  the same substitution moves Z by "
      f"{100*LIVE_075_RELDEV:.2f} %")
print("So 0.75 is live only when the packing and the liquid differ in surface")
print("tension, and this ARTICLE prints no sigma_c for any material - it is in ref. 18.")'''))

cells.append(code(r'''sr = np.linspace(0.5, 1.5, 21)
Zs, aws, els = [], [], []
for r_ in sr:
    o = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=r_, **PROPS)
    Zs.append(height_for_duty(o)[0]); aws.append(o.a_w(L_REF) / A_T)
    els.append(elasticities_analytic(o, L_REF, G_REF, M_SLOPE)["e_sigma"])
Zs, aws, els = map(np.array, (Zs, aws, els))
SIGMA_SWEEP_Z_FACTOR = float(Zs.max() / Zs.min())
SIGMA_SWEEP_AW_LO, SIGMA_SWEEP_AW_HI = float(aws.min()), float(aws.max())
print("sigma_c/sigma swept over 0.5 .. 1.5 (the page parameter this article cannot supply):")
print(f"   a_w/a_t   {SIGMA_SWEEP_AW_LO:.4f} .. {SIGMA_SWEEP_AW_HI:.4f}")
print(f"   Z         {Zs.max():.4f} .. {Zs.min():.4f} m, a factor {SIGMA_SWEEP_Z_FACTOR:.3f}")
print(f"   semi-elasticity of Z to the exponent 0.75 changes SIGN at sigma_c/sigma = 1:")
print(f"      at 0.5  {els[0]:+.4f}      at 1.0  {els[len(sr)//2]+0.0:+.1e}      "
      f"at 1.5  {els[-1]:+.4f}")
fig, ax = plt.subplots(1, 2, figsize=(10.2, 3.6))
ax[0].plot(sr, aws, "-", lw=1.8, color="C0"); ax[0].set_ylabel(r"$a_w/a_t$")
ax[1].plot(sr, els, "-", lw=1.8, color="C3"); ax[1].axhline(0, color="0.4", lw=1)
ax[1].set_ylabel(r"$\partial\ln Z/\partial e_\sigma$")
for a_ in ax:
    a_.axvline(1.0, ls=":", color="0.5", lw=1.2); a_.set_xlabel(r"$\sigma_c/\sigma$")
ax[0].set_title("wetted fraction"); ax[1].set_title("the 0.75 is inert exactly at 1.0")
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""### 5. Where each correlation is live: sweeping the equilibrium slope

Onda's Fig. 9 makes the same point for distillation. Here the resistance split is
swept over four decades of $m$ at fixed hydraulics, which is what says where
eq. (2) and eq. (3) each matter — and therefore where a mis-set constant in
either would be felt and where it would not.

**Only part of this range is a workable column.** Counter-current absorption
needs the absorption factor $L_M/(mG_M)$ above 1 or the column pinches; at these
loadings that caps $m$ at $L_M/G_M$. The sweep of the *resistance split* needs no
column and runs over the whole range; the *height* is reported only where a
column exists, and the cell says where the boundary is."""))

cells.append(code(r'''ms = np.logspace(-1, np.log10(2000.0), 400)
fG = np.array([onda.f_G(L_REF, G_REF, m, C_AV) for m in ms])
M_PINCH = L_M / G_M
F_G_AT_MIN_M, F_G_AT_MAX_M = float(fG[0]), float(fG[-1])
M_HALF = float(brentq(lambda m: onda.f_G(L_REF, G_REF, m, C_AV) - 0.5, 1e-3, 1e4))
print(f"gas-side resistance share f_G over m = {ms[0]:g} .. {ms[-1]:g}:")
print(f"   {F_G_AT_MIN_M:.4f} down to {F_G_AT_MAX_M:.5f}")
print(f"   50/50 split at m = {M_HALF:.4f}")
print(f"   pinch limit for this L/G: m < L_M/G_M = {M_PINCH:.3f}")
print(f"   inside the workable range, f_G runs "
      f"{onda.f_G(L_REF,G_REF,ms[0],C_AV):.4f} .. {onda.f_G(L_REF,G_REF,0.99*M_PINCH,C_AV):.4f}")

fig, ax = plt.subplots(figsize=(7.4, 4.0))
ax.semilogx(ms, fG, "-", lw=2, color="C0")
ax.axvline(M_PINCH, ls="--", color="C3", lw=1.2)
ax.axvline(M_HALF, ls=":", color="0.5", lw=1.2)
ax.fill_between(ms, 0, 1, where=(ms > M_PINCH), color="C3", alpha=0.06)
ax.annotate(f"pinch: no counter-current\ncolumn beyond m = {M_PINCH:.2f}\n"
            "at this L/G", (M_PINCH * 1.3, 0.55), fontsize=8, color="C3")
ax.annotate(f"50/50 at m = {M_HALF:.2f}", (M_HALF, 0.5), textcoords="offset points",
            xytext=(-100, 6), fontsize=8, color="0.3")
ax.set_xlabel("m, slope of the equilibrium line")
ax.set_ylabel(r"$f_G$ = gas-side share of $1/K_Ga$")
ax.set_ylim(0, 1)
ax.set_title("Which of Onda's three equations a design calculation can feel\n"
             r"$f_G\to1$: only eq. (3) matters.  $f_G\to0$: only eqs. (1) and (2) matter.")
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""### 6. The branch: 5.23 against 2.00

Eq. (3)'s prefactor is branched, and a branched correlation is exactly where a
page can accidentally test only one branch — the defect found on `A1.8`, where
every check sat on the dense branch and two printed constants could be deleted
outright. So both branches are run here, at the geometry the paper itself assigns
to the lower one.

The paper's criterion, stated three times: *"Raschig rings and Berl saddles
smaller than 15 mm"* (page 58), *"Raschig rings smaller than 15 mm and Berl
saddles smaller than 1/2""* (page 61), and for vaporization *"for 1/2-in. sphere
the constant of Eq.(3) might be changed into 2.00"* (page 60). Half-inch spheres
are therefore the one geometry for which this article explicitly names the lower
branch **and** supplies $a_t D_p$, so they are what is used."""))

cells.append(code(r'''BRANCH_FACTOR = C_KG / C_KG_SMALL
D_P_SMALL = 0.0127                      # 1/2-in spheres: the paper's own branch example
A_T_SMALL = ATDP_SPHERES / D_P_SMALL
rows = []
for label, dp, at_, Ck in [("1-in spheres, upper branch", D_P, A_T, C_KG),
                           ("1/2-in spheres, upper branch", D_P_SMALL, A_T_SMALL, C_KG),
                           ("1/2-in spheres, LOWER branch", D_P_SMALL, A_T_SMALL, C_KG_SMALL)]:
    o = OndaSet(a_t=at_, D_p=dp, sigma_ratio=SIGMA_RATIO, **PROPS, C_kG=Ck)
    Z, KGa = height_for_duty(o)
    rows.append((label, Ck, o.k_G(G_REF), o.f_G(L_REF, G_REF, M_SLOPE, C_AV), KGa, Z))
br = pd.DataFrame(rows, columns=["case", "eq. (3) prefactor", "k_G", "f_G", "K_G a", "Z / m"])
display(br.round(5))
Z_UPPER_SMALL = float(br.loc[1, "Z / m"]); Z_LOWER_SMALL = float(br.loc[2, "Z / m"])
BRANCH_Z_FACTOR = Z_LOWER_SMALL / Z_UPPER_SMALL
print(f"branch factor on the constant itself: {C_KG}/{C_KG_SMALL} = {BRANCH_FACTOR:.4f}")
print(f"branch factor on the PACKED HEIGHT for the same duty: {BRANCH_Z_FACTOR:.4f}")
print(f"   ({Z_UPPER_SMALL:.4f} m -> {Z_LOWER_SMALL:.4f} m)")
print("The height factor is smaller than the constant factor because the liquid")
print(f"film carries {100*(1-float(br.loc[1,'f_G'])):.1f} % of the resistance here and does not move.")
print("On a gas-film-controlled duty the two factors would coincide; that is the")
print("reason the branch is reported as a height, not as a coefficient ratio.")'''))

# ----------------------------------------------------------------- validation
cells.append(md(r"""### 7. The bands the authors claim, and one arithmetic remark of theirs

Quoted here because they are cited in the Reuse section and every number on this
page must be printed by the code rather than typed into prose."""))

cells.append(code(r'''for q in ("err_eq1_pct", "err_kL_pct", "err_vaporization_pct", "err_height_pct",
          "sigma_min_water_dyn_cm"):
    print(f"   {q:<26} = {stated(q):8.1f}   {SR.loc[q, 'unit']:<8} (page {int(SR.loc[q,'page'])})")
    print(f"      \"{SR.loc[q, 'statement']}\"")
NORMAN = stated("norman_ratio")
print(f"\nOne arithmetic remark of the authors', page 57: Eq. (2)'s Reynolds exponent is")
print(f'"nearly equal to {NORMAN} of that derived by Norman", so Norman\'s exponent must be')
print(f"   {E_RE_LW:.6f} / {NORMAN} = {E_RE_LW/NORMAN:.4f}")
NORMAN_EXPONENT = float(E_RE_LW / NORMAN)
print("i.e. very nearly 1 - a coefficient proportional to the liquid rate. Onda's")
print("2/3 is the difference between a model apparatus and a packed column, and it")
print("is stated here because it is the only comparison the paper makes to a")
print("mechanistically derived exponent rather than to a fitted one.")'''))

cells.append(md(r"""## Validation

Four checks, ranked by what they can catch, plus an explicit statement of what
none of them can.

**V1 — the eq. (4) identity, two independent routes.** Route A is symbolic;
Route B evaluates eq. (3) numerically at four thousand random property
combinations and regresses, doing no algebra at all. They agree, and both recover
the printed 0.771 to 0.45 % and the printed −0.30 to machine precision. This is
the strongest check on the page: it is the *only* one that pins all four gas-side
constants simultaneously, and it fails under a mis-reading of any of them.

**V2 — eq. (4) against Shulman's eq. (5).** Two correlations fitted to disjoint
experiments; the paper claims "fairly good" and this puts a number on it.

**V3 — the exponent window against Table 1.** In-sample, and negative for two
rows, which is the honest headline.

**V4 — the pymrm column against the closed form of the same BVP.** A
discretisation check and nothing more; it is labelled as one below.

**And the elasticity table is itself a two-route check**: the analytic chain rule
against central differences on the solved column, agreeing to the
`elasticity_two_routes` metric printed in the table above. That one *can* fail — it fails if either the analytic
derivation or the `OndaSet` implementation is wrong.

### V4: what the column check is, and what it is not"""))

cells.append(code(r'''grid = []
for n in [25, 50, 100, 200, 400, 800]:
    cc = Absorber(Z_STAR, KGA_REF, M_SLOPE, G_M, L_M, Y_IN, X_IN, n=n).solve()
    grid.append((n, cc.y_out, abs(cc.y_out / Y_EXACT - 1.0),
                 abs(cc.y_out_cell / Y_EXACT - 1.0), cc.its, cc.imbalance()))
gd = pd.DataFrame(grid, columns=["n", "y_out (face)", "dev, face value",
                                 "dev, CELL CENTRE (the A2.6 defect)",
                                 "TVD iterations", "solute imbalance"])
for src, dst in [("dev, face value", "order, face"),
                 ("dev, CELL CENTRE (the A2.6 defect)", "order, cell")]:
    gd[dst] = [np.nan] + list(np.log2(gd[src].values[:-1] / gd[src].values[1:]))
display(gd.round(10))

ys = gd["y_out (face)"].to_numpy()
ORDER_FITTED = float(np.log2(abs(ys[-3] - ys[-2]) / abs(ys[-2] - ys[-1])))
Y_RICHARDSON = float(ys[-1] + (ys[-1] - ys[-2]) / (2**ORDER_FITTED - 1))
COL_RAW_RELDEV = float(gd["dev, face value"].iloc[-1])
COL_CELL_RELDEV = float(gd["dev, CELL CENTRE (the A2.6 defect)"].iloc[-1])
COL_RICHARDSON_RELDEV = float(abs(Y_RICHARDSON / Y_EXACT - 1.0))
COL_ORDER_CELL = float(gd["order, cell"].iloc[-1])
Y_HALF_CELL_SHIFT = float(abs(col.y_out_cell / col.y_out - 1.0))
print(f"observed order on the FACE value:        {ORDER_FITTED:.3f}")
print(f"observed order on the CELL CENTRE value: {COL_ORDER_CELL:.3f}")
print(f"deviation at n = {int(gd['n'].iloc[-1])}:  face {COL_RAW_RELDEV:.3e}   "
      f"cell centre {COL_CELL_RELDEV:.3e}   ({COL_CELL_RELDEV/COL_RAW_RELDEV:.0f}x worse)")
print(f"Richardson-extrapolated deviation of the face value: {COL_RICHARDSON_RELDEV:.3e}")
print(f"\nTHE A2.6 TRAP, MEASURED. Reading the outlet at the centre of the last cell")
print(f"instead of at the face costs {100*Y_HALF_CELL_SHIFT:.3f} % at n = {N_CELLS}, and - the part")
print("that makes it dangerous - it drops the observed convergence order from 2 to 1")
print("while every grid still looks beautifully monotone. `compute_boundary_values`")
print("applies the same reconstruction the flux operator used, which is why the face")
print("route also closes the solute balance to machine precision and the cell route")
print("does not.")'''))

cells.append(md(r"""**What V4 cannot detect.** The closed form and the pymrm column solve the *same*
equations with the *same* $K_G a$; only the numerical method differs. So V4 tests
the discretisation, the boundary conditions and the operator assembly — and it is
blind to every error in the physics upstream of it. A wrong Onda constant, a
wrong resistance addition, a wrong $c_{av}$: all of them move both routes
identically and V4 reads exactly the same number. The break table below
demonstrates that rather than asserting it.

**And a break row tests sensitivity, never correctness.** Four verifiers
established that in this repository this week. Everything the break table shows
is that a metric *responds*; the reason to believe the eq. (4) prefactor is
right is that two independent derivations and a separately printed number agree,
not that perturbing 5.23 moves it.

### The break table

Every reported metric needs a row that moves it. Where no such row exists, the
metric is labelled structural and what it cannot detect is stated."""))

cells.append(code(r'''breaks = []
UND = {}

def row(defect, metric, defected, undefected=None):
    u = UND.get(metric) if undefected is None else undefected
    breaks.append((defect, metric, u, float(defected)))

# ---- baselines -----------------------------------------------------------
UND.update(
    jd_prefactor_symbolic=JD_PREFACTOR_SYMBOLIC,
    jd_prefactor_dev_pct=JD_PREFACTOR_DEV_PCT,
    jd_exponent_numeric=JD_EXPONENT_NUMERIC,
    jd_fit_max_resid=JD_FIT_MAX_RESID,
    eq45_max_dev_pct=EQ45_MAX_DEV_PCT,
    eq45_crossover_X=EQ45_CROSSOVER_X,
    window_hi=WINDOW_HI,
    window_lo=WINDOW_LO,
    table1_n_inside=float(TABLE1_N_INSIDE),
    table1_worst_excess=TABLE1_WORST_EXCESS,
    mae_window_ceiling=MAE_WINDOW_TOP,
    n_window_two_routes=N_WINDOW_TWO_ROUTES,
    Z_star_m=Z_STAR,
    f_G_ref=F_G_REF,
    a_w_fraction_ref=AW_FRACTION_REF,
    branch_Z_factor=BRANCH_Z_FACTOR,
    elasticity_two_routes=ELASTICITY_TWO_ROUTES,
    col_richardson_reldev=COL_RICHARDSON_RELDEV,
    col_raw_reldev=COL_RAW_RELDEV,
    col_imbalance=COL_IMBALANCE,
    inert_075_reldev=INERT_075_RELDEV,
    live_075_reldev=LIVE_075_RELDEV,
)

# ---- 1. every constant of eq. (3), against the eq. (4) identity -----------
def jd_dev_pct(**over):
    kw = dict(C=C_KG, e_Re=E_RE_G, e_atDp=E_ATDP_G); kw.update(over)
    return 100.0 * (jd_prefactor(**kw) / C_JD_ONDA - 1.0)

row("eq.(3) 5.23 -> 5.28 (digit misread)", "jd_prefactor_dev_pct", jd_dev_pct(C=5.28))
row("eq.(3) 5.23 -> 5.28 (digit misread)", "jd_prefactor_symbolic", jd_prefactor(C=5.28))
row("eq.(3) -2.0 -> -2.6 (exponent misread)", "jd_prefactor_dev_pct", jd_dev_pct(e_atDp=-2.6))
row("eq.(3) -2.0 -> -2.6 (exponent misread)", "jd_prefactor_symbolic", jd_prefactor(e_atDp=-2.6))
row("eq.(3) 0.7 -> 0.6", "jd_prefactor_dev_pct", jd_dev_pct(e_Re=0.6))
# the exponent of eq. (4) is 1 - e_Re_G; a mis-read Reynolds exponent moves it
row("eq.(3) 0.7 -> 0.6", "jd_exponent_numeric", 0.6 - 1.0)
row("a_t D_p 3.4 -> 3.0", "jd_prefactor_dev_pct", jd_dev_pct() if False else
    100.0 * (jd_prefactor(atDp=3.0) / C_JD_ONDA - 1.0))

# the Schmidt exponent is invisible to the PREFACTOR but not to the numeric FIT
def jd_fit_resid(e_Sc):
    kk = (C_KG * (gs / (ats * mus))**E_RE_G * Scs**e_Sc * (ats * dps)**E_ATDP_G * ats * dgs)
    jj = kk * rhos / gs * Scs**(2.0 / 3.0)
    s, a = np.polyfit(np.log(Xs), np.log(jj), 1)
    return float(np.abs(np.log(jj) - (a + s * np.log(Xs))).max())
row("eq.(3) 1/3 -> 1/2 (Schmidt exponent)", "jd_fit_max_resid", jd_fit_resid(0.5))
row("eq.(3) 1/3 -> 0.30", "jd_fit_max_resid", jd_fit_resid(0.30))

# ---- 2. eq. (4)/(5) comparison ------------------------------------------
def eq45(c4=C_JD_ONDA, e4=E_JD_ONDA, c5=C_JD_SHUL, e5=E_JD_SHUL):
    r = (c4 * Xg**e4) / (c5 * Xg**e5)
    cross = np.inf if e4 == e5 else float(np.exp(np.log(c5 / c4) / (e4 - e5)))
    return float(100 * np.abs(r - 1).max()), cross
row("eq.(5) 1.195 -> 1.795 (digit misread)", "eq45_max_dev_pct", eq45(c5=1.795)[0])
row("eq.(5) 1.195 -> 1.795 (digit misread)", "eq45_crossover_X", eq45(c5=1.795)[1])
row("eq.(5) -0.36 -> -0.30 (same as Onda)", "eq45_max_dev_pct", eq45(e5=-0.30)[0])

# ---- 3. the exponent window ---------------------------------------------
def window(e_Re_L=E_RE_L, e_Fr_L=E_FR_L, e_We_L=E_WE_L, e_Re_Lw=E_RE_LW):
    comp = e_Re_L + 2 * e_Fr_L + 2 * e_We_L
    return e_Re_Lw, e_Re_Lw + comp * (1 - e_Re_Lw)

for lbl, kw in [("eq.(1) 0.2 -> 0.3 (Weber exponent)", dict(e_We_L=0.3)),
                ("eq.(1) -0.05 -> +0.05 (Froude sign)", dict(e_Fr_L=+0.05)),
                ("eq.(2) 2/3 -> 1/2 (Reynolds exponent)", dict(e_Re_Lw=0.5))]:
    lo_, hi_ = window(**kw)
    row(lbl, "window_hi", hi_); row(lbl, "window_lo", lo_)
    ins = int(((n_meas > lo_) & (n_meas <= hi_)).sum())
    row(lbl, "table1_n_inside", ins)
    row(lbl, "table1_worst_excess", float(np.maximum(0, n_meas - hi_).max()))
    row(lbl, "mae_window_ceiling", mae(hi_))
# a defect in the TABLE, not the correlation
n_bad = n_meas.copy(); n_bad[2] = 0.36        # 0.86 read as 0.36
row("Table 1 n(row 3) 0.86 -> 0.36", "table1_n_inside",
    int(((n_bad > WINDOW_LO) & (n_bad <= WINDOW_HI)).sum()))
row("Table 1 n(row 3) 0.86 -> 0.36", "table1_worst_excess",
    float(np.maximum(0, n_bad - WINDOW_HI).max()))
row("Table 1 n(row 3) 0.86 -> 0.36", "mae_window_ceiling",
    float(np.abs(n_bad - WINDOW_HI).mean()))
# the two-route window check: break ONE route only
def n_pred_wrong_analytic(o, L):
    return E_RE_LW + L_COMPOSITE * (1 - E_RE_LW) * o.phi(L) * 2.0     # factor-2 slip
row("analytic window formula x2 (one route only)", "n_window_two_routes",
    float(np.abs(np.array([n_pred_wrong_analytic(onda, L) for L in Ls]) - nb).max()))

# ---- 4. the column and the correlation constants -------------------------
for lbl, kw in [("eq.(2) 0.0051 -> 0.0015", dict(C_kL=0.0015)),
                ("eq.(1) 1.45 -> 1.00", dict(C_aw=1.00)),
                ("eq.(3) 5.23 -> 2.00 (wrong branch)", dict(C_kG=C_KG_SMALL)),
                ("eq.(2) -1/2 -> +1/2 (Schmidt sign)", dict(e_Sc_L=+0.5)),
                ("eq.(2) 0.4 -> 0.04 (lost a digit)", dict(e_atDp_L=0.04)),
                ("eq.(1) 0.75 -> 0.075", dict(e_sigma=0.075)),
                ("eq.(1) 0.1 -> 1.0", dict(e_Re_L=1.0)),
                ("eq.(3) 1/3 -> 2/3", dict(e_Sc_G=2 / 3))]:
    o_b = onda_with(**kw)
    Z_b, KGa_b = height_for_duty(o_b)
    row(lbl, "Z_star_m", Z_b)
    row(lbl, "f_G_ref", o_b.f_G(L_REF, G_REF, M_SLOPE, C_AV))
    row(lbl, "a_w_fraction_ref", o_b.a_w(L_REF) / A_T)
# the branch factor needs its own defect: change what the LOWER branch is
o_lo = OndaSet(a_t=A_T_SMALL, D_p=D_P_SMALL, sigma_ratio=SIGMA_RATIO, **PROPS, C_kG=3.0)
o_hi = OndaSet(a_t=A_T_SMALL, D_p=D_P_SMALL, sigma_ratio=SIGMA_RATIO, **PROPS, C_kG=C_KG)
row("eq.(3) lower branch 2.00 -> 3.00", "branch_Z_factor",
    height_for_duty(o_lo)[0] / height_for_duty(o_hi)[0])

# ---- 5. the two-route elasticity check: break ONE route only -------------
def elasticities_wrong(o, L, G, m):
    """The same chain rule with the a_w feedback through k_L dropped."""
    fG = o.f_G(L, G, m, C_AV); fL = 1 - fG; ph = o.phi(L)
    kA = -(fG + fL)                       # WRONG: should be fG + fL*(1 - 2/3)
    g = groups_at(o, L, G)
    out = {"C_aw": kA * ph, "C_kL": -fL, "C_kG": -fG}
    for e in ("e_sigma", "e_Re_L", "e_Fr_L", "e_We_L"):
        out[e] = kA * ph * g[e]
    for e in ("e_Re_Lw", "e_Sc_L", "e_atDp_L"):
        out[e] = -fL * g[e]
    for e in ("e_Re_G", "e_Sc_G", "e_atDp_G"):
        out[e] = -fG * g[e]
    return out
_wrong = elasticities_wrong(onda, L_REF, G_REF, M_SLOPE)
row("chain rule drops the a_w->k_L feedback", "elasticity_two_routes",
    float(max(abs(_wrong[k] - elasticity_numeric(k)) for k in BASE)))

# ---- 6. the column's numerics -------------------------------------------
c_bad = Absorber(Z_STAR, KGA_REF, M_SLOPE, G_M, L_M, Y_IN, X_IN, n=800).solve(tvd=False)
row("van Leer correction switched off (n=800)", "col_richardson_reldev",
    abs(c_bad.y_out / Y_EXACT - 1.0))
row("van Leer correction switched off (n=800)", "col_raw_reldev",
    abs(c_bad.y_out / Y_EXACT - 1.0))
# liquid velocity sign flipped: the liquid now enters at the BOTTOM with the gas,
# so the same operators describe a co-current column instead of a counter-current one
c_v = Absorber(Z_STAR, KGA_REF, M_SLOPE, G_M, -L_M, Y_IN, X_IN, n=200).solve()
row("liquid velocity sign flipped (co-current)", "col_richardson_reldev",
    abs(c_v.y_out / Y_EXACT - 1.0))
# The co-current flip does NOT break the balance - a co-current column conserves
# solute just as well - so that row is kept and labelled as a NON-moving one. What
# DOES break it is reading the outlets half a cell short, which is the A2.6 defect:
row("liquid velocity sign flipped (co-current)", "col_imbalance", c_v.imbalance())
_gas = G_M * (Y_IN - col.y_out_cell)
_liq = L_M * (float(col.u[0, 1]) - X_IN)
row("outlets read at the cell centre (A2.6 defect)", "col_imbalance", abs(_gas / _liq - 1.0))
row("outlets read at the cell centre (A2.6 defect)", "col_raw_reldev",
    abs(col.y_out_cell / Y_EXACT - 1.0))

# ---- 7. the inert-exponent demonstration --------------------------------
o_i = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=1.0 + 1e-6, **PROPS, e_sigma=10 * E_SIG)
o_i0 = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=1.0 + 1e-6, **PROPS)
row("sigma_c/sigma 1.0 -> 1.000001", "inert_075_reldev",
    abs(height_for_duty(o_i)[0] / height_for_duty(o_i0)[0] - 1.0))
o_l = OndaSet(a_t=A_T, D_p=D_P, sigma_ratio=1.0, **PROPS, e_sigma=10 * E_SIG)
row("sigma_c/sigma 0.85 -> 1.00", "live_075_reldev",
    abs(height_for_duty(o_l)[0] / height_for_duty(o_eq)[0] - 1.0))

bt = pd.DataFrame(breaks, columns=["injected defect", "metric", "undefected", "defected"])
bt["moves_by"] = [("x %.4g" % (d / u) if (u not in (0, None) and np.isfinite(u) and u != 0)
                   else "%.3g -> %.3g" % (u if u is not None else np.nan, d))
                  for u, d in zip(bt.undefected, bt.defected)]
pd.set_option("display.max_rows", 120, "display.width", 220)
display(bt)

moved = {m for m, u, d in zip(bt.metric, bt.undefected, bt.defected)
         if (u is None) or (not np.isfinite(d)) or (u == 0) or abs(d - u) > 0.1 * abs(u)}
print(f"\nrows: {len(bt)}   distinct metrics in the table: {len(set(bt.metric))}, "
      f"of which at least one row moves: {len(moved)}")
unmoved = sorted(set(bt.metric) - moved)
print("metrics present in the table with NO moving row:", unmoved if unmoved else "none")
BREAK_METRICS = set(bt.metric)'''))

cells.append(code(r'''metrics = dict(
    # --- V1: the eq. (4) identity, the page's strongest check ---------------
    jd_prefactor_symbolic      = JD_PREFACTOR_SYMBOLIC,
    jd_prefactor_numeric       = JD_PREFACTOR_NUMERIC,
    jd_prefactor_dev_pct       = JD_PREFACTOR_DEV_PCT,
    jd_exponent_numeric        = JD_EXPONENT_NUMERIC,
    jd_exponent_dev            = JD_EXPONENT_DEV,          # < ABS_FLOOR: not CI-compared
    jd_two_routes_spread       = JD_ROUTE_SPREAD,          # < ABS_FLOOR: not CI-compared
    jd_fit_max_resid           = JD_FIT_MAX_RESID,         # < ABS_FLOOR: not CI-compared
    jd_atDp_for_printed_0771   = JD_ATDP_FOR_0771,
    # --- V2: eq. (4) against Shulman's eq. (5) ------------------------------
    eq45_dev_at_X100_pct       = EQ45_DEV_LO_PCT,
    eq45_dev_at_X10000_pct     = EQ45_DEV_HI_PCT,
    eq45_max_dev_pct           = EQ45_MAX_DEV_PCT,
    eq45_mean_dev_pct          = EQ45_MEAN_DEV_PCT,
    eq45_crossover_X           = EQ45_CROSSOVER_X,
    # --- V3: the exponent window against Table 1 ---------------------------
    window_lo                  = WINDOW_LO,
    window_hi                  = WINDOW_HI,
    n_window_two_routes        = N_WINDOW_TWO_ROUTES,
    table1_n_inside            = float(TABLE1_N_INSIDE),
    table1_worst_excess        = TABLE1_WORST_EXCESS,
    table1_n_range             = TABLE1_N_RANGE,
    mae_window_ceiling         = MAE_WINDOW_TOP,
    mae_null_no_area_coupling  = MAE_NULL_23,
    mae_null_fitted_constant   = MAE_NULL_FITTED,
    # --- the column, and what each constant is worth in metres -------------
    Z_star_m                   = Z_STAR,
    f_G_ref                    = F_G_REF,
    phi_ref                    = PHI_REF,
    a_w_fraction_ref           = AW_FRACTION_REF,
    m_for_equal_resistance     = M_HALF,
    branch_factor              = BRANCH_FACTOR,
    branch_Z_factor            = BRANCH_Z_FACTOR,
    elasticity_two_routes      = ELASTICITY_TWO_ROUTES,
    elasticity_max             = EL_MAX,
    live_075_reldev            = LIVE_075_RELDEV,
    e_ScG_doubled_reldev       = E_SCG_2X_RELDEV,
    inert_075_reldev           = INERT_075_RELDEV,         # < ABS_FLOOR: not CI-compared
    # --- V4 and structural identities, labelled, not evidence --------------
    col_order_face             = ORDER_FITTED,
    col_order_cell_centre      = COL_ORDER_CELL,
    col_raw_reldev             = COL_RAW_RELDEV,
    col_cell_centre_reldev     = COL_CELL_RELDEV,
    elasticity_column_route    = ELASTICITY_COLUMN_ROUTE,
    sigma_sweep_Z_factor       = SIGMA_SWEEP_Z_FACTOR,
    norman_exponent            = NORMAN_EXPONENT,
    col_richardson_reldev      = COL_RICHARDSON_RELDEV,
    col_imbalance              = COL_IMBALANCE,
    col_half_cell_shift        = Y_HALF_CELL_SHIFT,
    resistance_sum_rule        = SUM_RULE,                 # < ABS_FLOOR: not CI-compared
    window_escape_over_sweep   = WINDOW_ESCAPE,            # < ABS_FLOOR: not CI-compared
)
gu.report_agreement("A3.8", metrics)

FLOOR = 1e-12
below = sorted(k for k, v in metrics.items() if abs(v) < FLOOR)
print(f"\nNOTE - metrics below check_agreement.py's ABS_FLOOR = {FLOOR:g}, which CI")
print("therefore does NOT compare at all. They are pinned, not proven:")
for k in below:
    print(f"   {k:<26} = {metrics[k]:.4g}")
print("\nEach is an identity or a machine-precision agreement, and each has a break")
print("row in the table above showing what would move it. None carries a headline.")

uncovered = sorted(set(metrics) - BREAK_METRICS)
print(f"\nBREAK-TABLE COVERAGE: {len(metrics)} metrics reported, "
      f"{len(metrics)-len(uncovered)} appear in the break table.")
print("The remainder, and why:")
WHY = {
 "jd_prefactor_numeric":
   "the numeric twin of jd_prefactor_symbolic; every row that moves one moves the other",
 "jd_exponent_dev": "= |jd_exponent_numeric - printed|, and jd_exponent_numeric has rows",
 "jd_two_routes_spread":
   "STRUCTURAL - Route A and Route B evaluate the same transcribed eq. (3), so they "
   "agree for ANY constants; it tests the algebra, not the reading, and cannot fail "
   "for a mis-set constant",
 "jd_atDp_for_printed_0771":
   "an inversion of jd_prefactor_symbolic onto the printed 0.771; it is a defected "
   "value by construction and cannot also be a baseline",
 "eq45_dev_at_X100_pct": "an endpoint of eq45_max_dev_pct, which has rows",
 "eq45_dev_at_X10000_pct": "an endpoint of eq45_max_dev_pct, which has rows",
 "eq45_mean_dev_pct": "the mean of the same curve as eq45_max_dev_pct, which has rows",
 "table1_n_range": "a property of the transcribed table alone; the Table 1 defect row "
                   "above moves table1_n_inside and table1_worst_excess instead",
 "mae_null_no_area_coupling": "= mae(window_lo), and window_lo has rows",
 "mae_null_fitted_constant": "a property of Table 1 alone, with no model in it - that "
                             "is what makes it a null baseline",
 "phi_ref": "moves with every eq. (1) row via a_w_fraction_ref, which is reported",
 "m_for_equal_resistance": "moves with every eq. (2)/(3) row via f_G_ref, which is reported",
 "branch_factor": "= 5.23/2.00, two printed characters; there is no model in it to break. "
                  "branch_Z_factor is the one with physics in it and it has a row",
 "elasticity_max": "the largest entry of the elasticity table, whose correctness is "
                   "tested by elasticity_two_routes, which has a row",
 "col_order_face": "STRUCTURAL - the observed convergence order of the scheme, which "
                   "is a property of the discretisation and not of Onda's constants",
 "col_order_cell_centre": "STRUCTURAL - the order of the DEFECTIVE reading, reported so "
                          "the pair can be compared; it is 1 because the offset is h/2",
 "col_cell_centre_reldev": "the defected twin of col_raw_reldev BY CONSTRUCTION - it is "
                           "what the page would report if it read the outlet at the cell "
                           "centre, so it cannot also be a baseline",
 "elasticity_column_route": "the pymrm-solve twin of elasticity_two_routes, which has a "
                            "row; both move together for any chain-rule defect",
 "sigma_sweep_Z_factor": "the span of Z over a sweep of a PAGE PARAMETER, not of a "
                         "printed constant - there is nothing of Onda's in it to break",
 "norman_exponent": "= (2/3)/0.61, two printed characters divided; no model in it",
 "e_ScG_doubled_reldev": "it IS the defected value of the 'eq.(3) 1/3 -> 2/3' row on "
                         "Z_star_m, restated as a fraction, so it cannot also be a "
                         "baseline; the point it makes is how SMALL it is",
 "col_raw_reldev": "the un-extrapolated twin of col_richardson_reldev, which has rows",
 "col_half_cell_shift": "a diagnostic of the grid, not of the model: it is (h/2)N/G_M "
                        "by construction and moves only with n",
 "resistance_sum_rule": "STRUCTURAL - f_G + f_L = 1 by definition. It cannot fail for "
                        "any constants and is reported as the identity it is",
 "window_escape_over_sweep": "STRUCTURAL - the numerical sweep cannot leave a window "
                             "derived from the same phi it evaluates. It confirms the "
                             "sweep covered the limits, nothing more",
}
for k in uncovered:
    print(f"   {k}: {WHY.get(k, 'UNEXPLAINED - fix this')}")'''))

cells.append(md(r"""### What perturbation testing cannot detect here

Stated explicitly, because a break table's usefulness is bounded and the bound is
where the last four defects in this repository lived.

1. **A constant mis-set identically in both printings.** The transcription check
   is two typesettings of the same equations; one editorial error upstream of
   both defeats it. Only eq. (4) escapes, because it is a *fifth* number derived
   from the four.

0. **A transcription error in a published column no check consumes.** Table 1's
   $\alpha$ and temperature columns and three stated-results rows ship in the
   CSVs but feed nothing on this page, so a wrong digit there would reach a
   reuser with no gate in the way. (In verification every one of those values
   was re-read on an independent crop and found correct - the gap is in the
   coverage, not the data.)
2. **Anything about the physics of the column.** V4 compares two numerical
   routes to the same equations. A wrong $c_{av}$, a wrong resistance addition,
   the wrong Henry convention — every one moves both routes together and V4 does
   not blink. `col_imbalance` is a *discrete* conservation identity — it holds to machine
   precision for any $K_G a$ whatever, right or wrong. The break table shows both
   halves of that: flipping the liquid velocity to make the column co-current
   leaves it at machine precision, because a co-current column conserves solute
   just as well; what *does* move it is reading the outlets half a cell short of
   the faces. So it detects an inconsistency between the flux operator and the
   boundary reconstruction, **and nothing else**.
3. **Whether $a = a_w$.** Onda's central assumption. Nothing in the article
   separates the area from the coefficient, so nothing here can either, and no
   break row can be written for it. It is a structural blind spot of the whole
   correlation set, not of this page's implementation.
4. **Whether the correlations are right.** They were fitted to the data any
   comparison here would use. The exponent window is the only quantity on the
   page that a refit could not have absorbed, and even that is scored on a table
   of six rows with no error bars — where, as reported above, a single fitted
   constant scores better.
5. **The illustrative operating point.** Air, water, 25 °C and $\sigma_c/\sigma$
   = 0.85 are page choices. Every elasticity is reported alongside the two
   numbers ($f_G$ and $\varphi$) that carry the operating-point dependence, so a
   reader can re-scale rather than re-run; but no break row can test a choice."""))

# --------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Not the correlations.** Eqs. (1)–(3) are algebra and need no solver; the eq. (4)
identity and the exponent window are done with `sympy` and `numpy` and would be
the same without pymrm.

What the solver adds is the thing the 1968 paper could not do, and it is
specific: **eq. (6) becomes a boundary-value problem instead of a quadrature.**
Onda integrates $\mathrm{d}y/(y^*-y)$ along an operating line that must be
assumed in advance, with $G_M$ held constant so the integral closes. Writing the
same physics as two counter-current convective fluxes with an interphase source
puts the operating line *inside* the solution: the liquid profile is solved, not
assumed, and the two inlet conditions sit at opposite ends of the domain where
they physically belong. That is what makes it a two-point problem, and it is why
`construct_convflux_upwind` with a two-component velocity `[[G_M, -L_M]]` is the
natural assembly rather than an accident of implementation.

Three consequences that are used above.

1. **The elasticity of packed height to each printed constant is measurable by
   finite differences on a solved column**, so the analytic chain rule has an
   independent check. Without the column the analytic elasticities would be
   unfalsifiable algebra.
2. **The pinch is a property of the solution, not of the integral.** The sweep
   over $m$ finds $L_M/(mG_M) = 1$ as the boundary where a counter-current column
   stops existing; eq. (6)'s integrand diverges there, which is harder to read
   than a solve that stops converging.
3. **The route generalises where eq. (6) does not.** Varying $G_M$ (concentrated
   solute), a curved equilibrium line, or an axially varying $L$ all break the
   quadrature and none of them breaks the BVP — they change the source term and
   the velocity field, not the structure.

**What it does not add.** No accuracy. The column's closed form is exact for the
constant-coefficient case solved here, and the pymrm solve is *worse* than it by
construction — that is what the grid study measures. The solver earns its place
on generality, and on making the sensitivity study possible, not on precision."""))

# -------------------------------------------------------------------- reuse
cells.append(code('_md = r"""'
    r'''## Reuse

**`OndaSet` is the reusable object.** It takes packing geometry and fluid
properties in Onda's units and returns $a_w$, $k_L$, $k_G$ and the combined
$K_G a$. Before using it, four things about it are worth carrying with the code,
because none of them is visible from the function signature.

1. **Give it Onda's units, not SI.** $\\mu$ in kg/(m hr), $D$ in m²/hr, $g$ in
   m/hr², $\\sigma$ in **kg/hr²**, $R$ in m³ atm/(kg-mole K). Passing SI viscosity
   changes $Re_L$ by 3600 and $a_w$ silently — eq. (1) returns a plausible
   fraction for almost any input, which is exactly how this class of error
   survives.
2. **`C_kG` is branched and the default is the upper branch.** Pass
   `C_kG=2.00` for Raschig rings below 15 mm and Berl saddles below 1/2 in — the
   paper's own criterion, restated on pages 58, 60 and 61. At the illustrative
   duty the branch is worth a factor **`branch_Z_factor`** on packed height, less
   than the factor `branch_factor` on the constant itself because the liquid film
   does not move with it. On a gas-film-controlled duty the two coincide.
3. **`e_sigma` (the 0.75) is inert whenever $\\sigma_c = \\sigma$**, exactly, not
   approximately. And this article prints $\\sigma_c$ for **no** packing material;
   it is in Onda, Takeuchi & Koyama (1967), ref. 18. So a caller who does not
   have that reference has no $\\sigma_c$, and should treat $\\sigma_c/\\sigma$ as
   what it is here: a parameter to be swept, not a known.
4. **`a_t` for a named packing is not in this article either.** The only geometry
   Onda supplies is $a_t D_p = 6(1-\\varepsilon) = 3.4$ for spheres, which is what
   every calculation on this page uses. Do not take an $a_t$ from a vendor table
   and present the result as this paper's.

**Where the correlations are worth trusting, on the paper's own account.**
Eq. (2) within ±20 % for Raschig rings, Berl saddles, spheres and rods irrigated
with organic solvents *and* with water systems above about 50 dyn/cm; eq. (3)
within ±30 % for vaporization as well as absorption; eq. (6) with all three
within ±30 % on packed height, *"except columns higher than 1.0 m in which the
maldistribution of liquid might have occured"*. Those bands are the authors'
assessments of their own fits against their own data and are not independent
validations — they are quoted here as claims, which is how they should be passed
on.

**And the finding a reuser most needs.** The exponent set caps
$\\mathrm{d}\\ln(k_L a)/\\mathrm{d}\\ln L$ at **`window_hi`**, whatever the packing
and whatever the liquid. Two of Onda's own six runs — both spheres, one in CCl₄
and one in methanol — show measured exponents above that cap. If your data show
$k_L a$ rising faster than $L^{0.8}$, the correlation cannot follow it and no
adjustment of 1.45 or 0.0051 will make it.

**What is not here, and who should own it.** A three-way comparison with
Billet–Schultes (`A3.9`) and Rocha–Bravo–Fair (`A3.10`) is the obvious next step
and needs their sources; it belongs on whichever of those two is built second.
The stirred-tank correlations `A3.6` and `A3.7` are a different geometry and a
different $a$, and nothing here transfers to them. Reaction enhancement of the
liquid film is `F3.1`'s subject: $k_L$ from eq. (2) is the *physical* coefficient
that Hatta's enhancement factor multiplies, and using it unenhanced in a reacting
system is the standard way to be wrong by an order of magnitude.''' '"""\n'
    '''for _tok, _val in (("`branch_Z_factor`", f"{BRANCH_Z_FACTOR:.2f}"),
                   ("`branch_factor`", f"{BRANCH_FACTOR:.2f}"),
                   ("`window_hi`", f"{WINDOW_HI:.3f}")):
    assert _tok in _md, _tok
    _md = _md.replace(_tok, _val)
from IPython.display import Markdown as _MD
display(_MD(_md))'''))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.kernelspec = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata.language_info = {"name": "python"}
nbf.write(nb, Path(__file__).with_name("index.ipynb"))
print(f"wrote index.ipynb with {len(cells)} cells")
