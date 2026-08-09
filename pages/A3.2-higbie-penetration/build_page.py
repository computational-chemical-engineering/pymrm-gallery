#!/usr/bin/env python3
"""Generate index.ipynb for page A3.2 (Higbie penetration theory).

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells: list = []

# ---------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Higbie's penetration model: what a fixed contact time predicts, and the one measurement on disk that can test it"
description: "Replace the steady film by a surface element exposed for a fixed time and k_L becomes 2 sqrt(D / pi t_exp) - with the contact time supplied by the hydrodynamics, not fitted. Read from Bird, Stewart & Lightfoot section 18.5; cross-checked against Froment, De Wilde & Bischoff section 6.4; and finally set against Whitman's film and Danckwerts' surface renewal on the one printed table that carries all three."
categories: [sec:A, struct:S4, tier:T0, data:tier6, phase:gas-liquid]
date: 2026-08-09
---

# Higbie's penetration model

A liquid surface is exposed to a gas for a definite length of time, absorbs by
unsteady diffusion the whole while, and is then swept away and replaced. Nothing
reaches a steady state; there is no film thickness anywhere in the argument. The
mass-transfer coefficient that comes out is

$$k_L = 2\sqrt{\frac{\mathscr{D}_{AB}}{\pi\,t_{\rm exp}}},$$

so $k_L \propto \mathscr{D}^{1/2}$ where film theory gives $k_L \propto \mathscr{D}$.

**Where this page reads the model from, and where it does not.** R. Higbie,
*Trans. AIChE* **31**, 365-389 (1935) is the origin. It is pre-DOI, it is not on
disk, and **it was not consulted**. Everything below is read from
**Bird, Stewart & Lightfoot, *Transport Phenomena*, 2nd edn (Wiley, 2002),
section 18.5** (book pp. 558-561), which derives the result, names it the
"penetration model", attributes it to Higbie in the same paragraph, and prints
the full citation with a biographical footnote. A second book on disk,
**Froment, De Wilde & Bischoff, *Chemical Reactor Analysis and Design*, 3rd edn
(Wiley, 2011), section 6.4**, is used as a **cross-check only**: it treats
Higbie's picture as the uniform-surface-age special case of surface renewal and
prints the same constant and a closed-form enhancement factor, independently.
Where the two books disagree, that is reported as a finding, not adjudicated.

**What this page has that its two siblings did not.**
[`A3.1`](../A3.1-whitman-two-film/) (Whitman's two films) and
[`A3.3`](../A3.3-danckwerts-surface-renewal/) (Danckwerts' surface renewal) both
state that the three-way comparison between the pictures is open, and both
decline it - correctly, because neither source contains anything that could
separate them. This page is the last of the three, so it takes the comparison,
and it takes it with two things the others did not have: **one measured
mass-transfer coefficient** (Hammerton & Garner's 117 cm/hr for CO2 bubbles,
quoted in BSL's Problem 18A.7) and **one printed table carrying all three models
side by side** (Froment's Table 6.4.2.1). What those two can and cannot settle is
the subject of *Results*.
"""))

# ------------------------------------------------------------------ background
cells.append(md(r"""## Background

Whitman's film (1923) explains a liquid-side resistance by putting a stagnant
layer of thickness $\delta$ against the interface and letting it reach steady
state: $k_L = \mathscr{D}/\delta$. Section 18.5 sets up a different picture, and
BSL states it as a restriction rather than as an objection: the diffusion is to
take place "so slowly in the liquid film that $A$ will not 'penetrate' very far
into the film - that is, that the penetration distance will be small in
comparison with the film thickness." In the two systems the section and its
worked Example are about - liquid running down a wetted wall, liquid sliding
round a rising bubble - **the contact is too short for a steady state to be
established at all**: the solute penetrates a short distance into a surface that
is moving past, and then that piece of surface is gone.

**Whose argument that is, exactly.** All BSL prints about the 1935 paper is one
sentence - "This approach for studying gas absorption was apparently first
proposed by Higbie" - plus footnote 5, which is the citation and a biographical
note that also credits him with having "provided the basis for the 'penetration
model' of mass transfer". It does not say what Higbie's own argument was, which
geometries he treated, or what he compared his result against. The origin was not
consulted here, so **nothing on this page is attributed to the 1935 paper beyond
what BSL prints about it and, where Froment's cross-check is described, what
Froment prints about it** - his "Higbie's uniform age $\bar t$", printed
immediately above his Eq. 6.4.1-10, and his opening of section 6.4. Every
problem, example and contact time below is BSL's or Froment's, and is named as
such.

That changes what has to be solved. The film gives an ordinary differential
equation in one steady coordinate; penetration gives the transient diffusion
equation, with the time variable supplied by the flow. It also changes what is
*adjustable*. A film thickness is not observable and must be fitted. **A contact
time, in the two examples BSL applies the result to, is not fitted: it is
$L/v_{\max}$ for the falling film of section 18.5, and $D/v_t$ for the bubble of
diameter $D$ rising at $v_t$ in BSL's Example 18.5-1** - both set by the geometry
and the hydrodynamics. That is the difference this page is built around, because
it is the difference that makes one of the three pictures falsifiable on a single
measurement and the other two not.

**Where the sibling pages leave the question.** `A3.1` establishes that Whitman's
three runs cannot see the diffusivity exponent, and says the comparison is open.
`A3.3` goes further and proves *why* its own source cannot settle it: Danckwerts'
one physical number, $s \approx 5\ {\rm s^{-1}}$, is obtained by inverting
$k_L=\sqrt{\mathscr{D}s}$, and film theory fits the same single datum exactly as
well by inverting $k_L=\mathscr{D}/\delta$ - **one free constant each, one datum,
zero residual each**. That argument is `A3.3`'s and is not re-derived here; what
this page adds is a case in which the penetration model has *no* free constant to
absorb the datum with. Neither sibling's CSV is loaded here: nothing in them is an
input to anything below, and their numbers are referred to as theirs.

**A source that does not qualify.** Levenspiel's *Chemical Reaction Engineering*
(3rd edn, 1999) is on disk and mentions Higbie - "Alternatives to the film theory
are also in use. These models [Higbie (1935); Danckwerts (1950, 1955)] view that
the liquid at the interface is continually washed ..." - with the reference in the
chapter list. Named and attributed, but neither derived nor tested. It is not used
for this case.
"""))

# ------------------------------------------------------------ published model
cells.append(md(r"""## The published model

Everything in this section is transcribed from BSL section 18.5 on 300 ppi
renders of the scan at its native resolution. The equation numbers are BSL's.

**The problem.** Gas $A$ is absorbed by a laminar film of liquid $B$ running down
a wall, $A$ only slightly soluble so the velocity field is unaffected. The
velocity profile is Eq. 18.5-1,

$$v_z(x) = v_{\max}\left[1-\left(\frac{x}{\delta}\right)^{2}\right],$$

$x$ measured from the gas-liquid interface into the film of thickness $\delta$.
A shell balance and the two flux approximations - convection only in $z$,
diffusion only in $x$ - give Eq. 18.5-7,

$$v_{\max}\left[1-\left(\frac{x}{\delta}\right)^{2}\right]\frac{\partial c_A}{\partial z}
  = \mathscr{D}_{AB}\frac{\partial^{2}c_A}{\partial x^{2}},$$

with B.C. 1 $c_A=0$ at $z=0$ (Eq. 18.5-8), B.C. 2 $c_A=c_{A0}$ at $x=0$
(Eq. 18.5-9), and B.C. 3 $\partial c_A/\partial x=0$ at $x=\delta$
(Eq. 18.5-10) - $A$ cannot diffuse through the wall. BSL notes that this problem
"has been solved analytically in the form of an infinite series" (footnote 2:
R. L. Pigford, PhD thesis, University of Illinois, 1941) **but does not give that
solution**.

**The penetration approximation.** Instead BSL argues that if $A$ has penetrated
only a short distance, it "has the impression" that the film moves everywhere at
$v_{\max}$, and it does not "sense" the wall - so Eq. 18.5-7 and its boundary
conditions are replaced by Eq. 18.5-11,

$$v_{\max}\frac{\partial c_A}{\partial z} = \mathscr{D}_{AB}\frac{\partial^{2}c_A}{\partial x^{2}},$$

with $c_A=0$ at $z=0$, $c_A=c_{A0}$ at $x=0$, and $c_A=0$ at $x=\infty$
(Eqs. 18.5-12 to 14). **Two approximations are made at once here, and section
*What pymrm adds* separates them.** The solution is Eq. 18.5-16,

$$\frac{c_A}{c_{A0}} = \operatorname{erfc}\frac{x}{\sqrt{4\mathscr{D}_{AB}z/v_{\max}}},$$

the interfacial flux is Eq. 18.5-17,

$$N_{Ax}\big|_{x=0} = c_{A0}\sqrt{\frac{\mathscr{D}_{AB}v_{\max}}{\pi z}},$$

and the total absorbed by a film of length $L$ and width $W$ is Eq. 18.5-18,

$$W_A = WLc_{A0}\sqrt{\frac{4\mathscr{D}_{AB}v_{\max}}{\pi L}}.$$

BSL adds that "the same result is obtained by integrating the product
$v_{\max}c_A$ over the flow cross section at $z=L$ (see Problem 18C.3)" - a mass
balance rather than a surface-flux integral, and the route the numerical solve
below uses.

**The attribution, and the exposure time.** Immediately after Eq. 18.5-18:
"Equation 18.5-18 shows that the mass transfer rate is directly proportional to
the square root of the diffusivity and inversely proportional to the square root
of the 'exposure time,' $t_{\rm exp}=L/v_{\max}$. This approach for studying gas
absorption was apparently first proposed by Higbie.$^5$ The problem discussed in
this section illustrates the 'penetration model' of mass transfer. This model is
discussed further in Chapters 20 and 22." Footnote 5 is the full Higbie citation.
Writing Eq. 18.5-18 as $W_A = WL\,k_L c_{A0}$ gives the coefficient quoted at the
top of this page.

**The bubble, and the one place a measurement enters.** Example 18.5-1 applies the
same result to a gas bubble with internal (Rybczynski-Hadamard) circulation by
replacing $t_{\rm exp}=L/v_{\max}$ with $D/v_t$, $D$ the bubble diameter, giving
Eq. 18.5-19,

$$(N_A)_{\rm avg} = \sqrt{\frac{4\mathscr{D}_{AB}v_t}{\pi D}}\;c_{A0}.$$

BSL: this "turns out to be correct for potential flow of the liquid around the
bubble ... This equation has been approximately confirmed$^6$ for gas bubbles 0.3
to 0.5 cm in diameter rising through carefully purified water", footnote 6 being
D. Hammerton and F. H. Garner, *Trans. Inst. Chem. Engrs. (London)* **32**,
S18-S24 (1954). For **creeping** flow the corresponding result is Eq. 18.5-20,

$$(N_A)_{\rm avg} = \sqrt{\frac{4\mathscr{D}_{AB}v_t}{3\pi D}}\;c_{A0},$$

attributed to Levich, lower than Eq. 18.5-19 by exactly $\sqrt3$. BSL also records
that a surfactant skin, by killing the internal circulation, drops the exponent to
$\mathscr{D}^{1/3}$ "as for a solid sphere", and section 18.6 - a solid dissolving
into a falling film, where the velocity gradient at the interface is finite rather
than zero - gives $\mathscr{D}^{2/3}$. **Four exponents, then, in one chapter:**
$\mathscr{D}^{1}$ (film), $\mathscr{D}^{1/2}$ (penetration, mobile interface),
$\mathscr{D}^{2/3}$ (linear velocity gradient) and $\mathscr{D}^{1/3}$ (rigid
sphere).

**The cross-check.** Froment, De Wilde & Bischoff section 6.4 arrives at the same
constant from the other direction. Treating Higbie's uniform surface age $\bar t$
as a special case of the surface-age distribution, their Eq. 6.4.1-10 is followed
by: "since, for purely physical absorption, $N_A = 2C_{Ai}\sqrt{D_A/\pi\bar t}$,
but also $N_A = k_L C_{Ai}$." That is BSL's coefficient, from a different book,
a different derivation and a different starting point. They also print, as their
Eq. 6.4.2-8, the enhancement factor Higbie's uniform age gives for a
pseudo-first-order reaction,

$$F_A = \gamma\left[\left(1+\frac{\pi}{8\gamma^{2}}\right)
        \operatorname{erf}\!\left(\frac{2}{\sqrt\pi}\gamma\right)
        + \frac{1}{2\gamma}\exp\!\left(-\frac{4}{\pi}\gamma^{2}\right)\right],
  \qquad \gamma^{2}=\frac{D_Ak}{k_L^{2}},$$

against $\sqrt{1+\gamma^2}$ for Danckwerts' surface renewal (their Eq. 6.4.2-7)
and $\gamma/\tanh\gamma$ for the film. Their Table 6.4.2.1 tabulates all three.
"""))

# -------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

**Nothing on this page uses a transport property that is not printed in a source
on disk.** The two worked problems supply their own diffusivities, solubilities
and velocities; the one molecular weight needed comes from BSL's own Table E.1.

- **`18A.4`, chlorine into a falling water film.** $R=1.4$ cm and $L=13$ cm are
  text labels on Fig. 18A.4 (nothing is traced or digitised); $\langle v\rangle
  =17.7$ cm/s, $\mathscr{D}_{AB}=1.26\times10^{-5}$ cm$^2$/s and a saturation of
  0.823 g Cl$_2$ per 100 g water at 16 C are printed in the statement. The wetted
  width is $W=2\pi R$.
- **The factor $3/2$ is derived, not remembered.** Eq. 18.5-18 is written in
  $v_{\max}$ and the problem prints $\langle v\rangle$. Averaging the *printed*
  Eq. 18.5-1 across the film,
  $\langle v_z\rangle = \frac{1}{\delta}\int_0^{\delta}v_{\max}[1-(x/\delta)^2]\,dx
  = \tfrac23 v_{\max}$, so $v_{\max} = \tfrac32\langle v\rangle$.
- **One assumption, flagged.** Converting "0.823 g per 100 g water" to
  mol cm$^{-3}$ needs the density of water, which the problem does not print. It
  is taken as exactly 1.000 g/cm$^3$ and the solution treated as dilute. The
  break table contains a row using 0.999 instead. It **does** move the
  reproduction - by about a tenth of a percent, which is the whole size of the
  residual against BSL's printed three figures - and the two densities fall on
  **opposite sides** of a three-figure rounding step, so the choice is not
  harmless. What the printed answer still cannot do is **decide between them**:
  the rounding of the inputs 18A.4 *does* print - the diffusivity on its own -
  moves the answer by more than the gap between the two densities. *Results 1*
  shows that arithmetic rather than asserting it.
- **`18A.7`, CO$_2$ from a rising bubble.** $D=0.5$ cm, $v_t=22$ cm/s,
  $\mathscr{D}_{AB}=1.46\times10^{-5}$ cm$^2$/s and $c_{A0}=0.041$ g-mol/liter at
  18 C, 1 atm, all printed; and Hammerton & Garner's measured surface-averaged
  $k_c=117$ cm/hr, quoted by BSL from a paper that is **not on disk and was not
  consulted**.
- **$\delta$ is never invented.** Problem 18A.4 prints no film thickness, no
  viscosity and no flow rate, so the finite-thickness correction is reported as a
  threshold in the dimensionless group $\Lambda = \mathscr{D}_{AB}L/(v_{\max}
  \delta^{2})$ and as the $\delta$ at which that correction reaches 1 %, rather
  than as a number for that column.
- **Seeding.** Nothing on this page is stochastic. Two consecutive executions
  give identical content and an identical `agreement.json`.
"""))

# ------------------------------------------------------------------ colab cell
cells.append(code(r'''# Colab environment cell.
try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm pandas pyyaml matplotlib'''))

cells.append(code(r'''import sys, pathlib, urllib.request

_shared = pathlib.Path.cwd()
for _ in range(4):
    if (_shared / "shared" / "gallery_utils.py").is_file():
        sys.path.insert(0, str(_shared / "shared"))
        break
    _shared = _shared.parent
else:
    if "google.colab" in sys.modules:
        urllib.request.urlretrieve(
            "https://raw.githubusercontent.com/computational-chemical-engineering/"
            "pymrm-gallery/main/shared/gallery_utils.py", "gallery_utils.py")
        sys.path.insert(0, ".")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import eye_array
from scipy.linalg import solve_banded
from scipy.special import erf, erfc
from scipy.optimize import brentq
from scipy.integrate import quad
from IPython.display import display, Markdown

from pymrm import construct_grad, construct_div
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A3.2-higbie-penetration"
M: dict[str, float] = {}          # every reported metric lands here, once
np.set_printoptions(legacy="1.25")
plt.rcParams.update({"figure.dpi": 110, "font.size": 9})'''))

# ------------------------------------------------------------------- the data
cells.append(md(r"""## The data

Three CSVs, all of them printed numbers rather than measurements, which is why
the data tier is **6**. Two come from the pinned source (BSL), one from the
cross-check (Froment).

1. `bsl-2002-sec18p5-worked-problems.csv` - the inputs and printed answers of
   Problems 18A.4 and 18A.7, plus the molecular weight of Cl$_2$ from BSL's
   Table E.1. One row, `bub_kc_measured`, is a **measurement**: Hammerton &
   Garner's 117 cm/hr, quoted by BSL. That paper is not on disk and was not
   consulted; BSL prints no error bar for it and none is invented.
2. `bsl-2002-penetration-vs-film-tables.csv` - BSL's Table 20.1-1 (Arnold's
   $\varphi$ and $\psi$) and Table 22.8-2, headed *Comparison of Film and
   Penetration Models*. **That table is not about the exponent on diffusivity**;
   it compares the two models' corrections $\theta_x$ for a finite net
   mass-transfer rate. The page says so rather than letting it look like more
   than it is. The Arnold problem itself is case `A4.8`; $\varphi$ and $\psi$
   appear here only because Eq. 22.8-41 needs them.
3. `froment-2011-table6p4p2-1-three-models.csv` - Froment's Table 6.4.2.1,
   film / surface renewal / penetration enhancement factors at six $\gamma$.
   Froment prints it "After Beek [1968]"; **Beek 1968 is not on disk and was not
   consulted**, so every statement made about these numbers is a statement about
   Froment's printing.

**No other page's dataset is loaded.** `A3.1`'s Whitman table and `A3.3`'s
Danckwerts numbers are about different systems and are not inputs to anything
here; where those pages' *conclusions* are used, they are attributed in prose and
no number of theirs is retyped as if it were this page's own.
"""))

cells.append(code(r'''prob = load_data("bsl-2002-sec18p5-worked-problems.csv", page=PAGE).set_index("id")
tabs = load_data("bsl-2002-penetration-vs-film-tables.csv", page=PAGE)
fro  = load_data("froment-2011-table6p4p2-1-three-models.csv", page=PAGE)

print(cite_data(load_meta("bsl-2002-sec18p5-worked-problems.csv", page=PAGE)))
print(cite_data(load_meta("froment-2011-table6p4p2-1-three-models.csv", page=PAGE)))

# Pull BSL's own values out of the CSV rather than retyping them. Nothing in this
# notebook types a number that is a row of a file it loaded.
def P(key):
    """A value BSL prints, by name."""
    return float(prob.loc[key, "value"])

display(prob[["problem", "quantity", "value", "unit", "role"]])
display(tabs)
display(fro)'''))

cells.append(md(r"""### What the sources say about these rows, and what it costs

- BSL prints Problem 18A.7's answers to **three** significant figures and
  18A.4's to three as well. Every reproduction below is therefore at best a
  0.2 %-level statement - and 18A.4's residual is *at* that level, not under it.
  The break table carries the two unprinted inputs that live at the same scale,
  the density of water and the molecular weight of Cl$_2$; each moves the
  reproduction by about a tenth of a percent, which is the point rather than an
  embarrassment, because it means the printed answer cannot decide between them.
- Section 18.5 itself carries **no numerical result**. The only figures printed
  anywhere in it are the "0.3 to 0.5 cm" bubble-diameter range over which BSL
  says Eq. 18.5-19 was approximately confirmed, and the $\mathscr{D}^{1/3}$
  exponent of a surfactant-skinned bubble, which book p. 561 contrasts with the
  $\mathscr{D}^{1/2}$ of Eq. 18.5-19. (The $\mathscr{D}^{2/3}$ belongs to
  **section 18.6**, not to 18.5: book p. 563 is where BSL writes "in §18.5,
  $W_A\propto(\mathscr{D}_{AB}L)^{1/2}$, whereas in this section
  $W_A\propto(\mathscr{D}_{AB}L)^{2/3}$".) Everything else is derivation and
  attribution; the numbers live in the problems.
- The 117 cm/hr is a *surface-averaged* coefficient over a rising bubble, and
  BSL's own words for the agreement are "approximately confirmed". The page does
  not upgrade that to a validation.
- **Table 22.8-2's penetration column is not an independent computation.** It is
  Table 20.1-1's $\psi$ column multiplied by $(1-x_{A0})$, which is exactly what
  Eq. 22.8-41 says. The two reproductions below therefore carry the same
  information, to the last digit, and the notebook prints that coincidence
  instead of presenting them as two checks.
- Froment's Table 6.4.2.1 rows at $\gamma=0.01$ and $\gamma=10$ are examined
  below; they are inconsistent with the closed form printed at the head of their
  own column, and they are reported, not corrected.
"""))

# ----------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

The structure is `S4`: a one-dimensional transient (here, *marching*) diffusion
problem. The marching coordinate is the axial distance $z$, which plays the role
of time because Eq. 18.5-7 has no axial diffusion.

Non-dimensionalise with $\xi=x/\delta$, $\zeta=z/L$, $u=c_A/c_{A0}$ and

$$\Lambda \equiv \frac{\mathscr{D}_{AB}L}{v_{\max}\delta^{2}},
\qquad\text{so}\qquad
(1-\xi^{2})\frac{\partial u}{\partial\zeta} = \Lambda\frac{\partial^{2}u}{\partial\xi^{2}},$$

with $u(\xi,0)=0$, $u(0,\zeta)=1$ and $\partial u/\partial\xi|_{\xi=1}=0$. The
penetration depth relative to the film is $\sqrt{4\Lambda}$, so $\Lambda\to0$ is
BSL's short-contact-time limit. One class covers every case the page needs, by
switching two things independently:

- the **capacity** $(1-\xi^{2})$, the parabolic velocity profile, or $1$ for plug
  flow at $v_{\max}$ - this is BSL's first approximation;
- the **outer boundary**, a no-flux wall at $\xi=1$ or a domain long enough that
  the wall is never reached - this is BSL's second approximation.

Turning both off gives Eq. 18.5-11 exactly; turning both on gives Eq. 18.5-7
exactly. The dimensionless uptake reported at the outlet is the mixing cup

$$\Phi \equiv \int_{0}^{1}(1-\xi^{2})\,u(\xi,1)\,d\xi
  \qquad\Longrightarrow\qquad W_A = W\,\delta\,v_{\max}c_{A0}\,\Phi,$$

which is BSL's Problem 18C.3 route. In the penetration limit
$\Phi\to\sqrt{4\Lambda/\pi}$ and $\delta$ cancels out of $W_A$, as it must.

Two pymrm conventions matter here and both are commented in the code: the
boundary conditions use the **outward** normal, so the same dict means opposite
things at the two ends; and the shape passed around is `(n_x, 1)`, never a bare
`(n_x,)`.
"""))

cells.append(code(r'''def _banded(A):
    """Tridiagonal csr -> the (3, n) layout scipy.linalg.solve_banded wants."""
    A = A.tocsr()
    n = A.shape[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = A.diagonal(1)
    ab[1, :] = A.diagonal(0)
    ab[2, :-1] = A.diagonal(-1)
    return ab


class Film:
    """Absorption into a falling liquid film, marched in the axial coordinate.

        cap(xi) du/dzeta = Lam d2u/dxi2,   u(xi, 0) = 0

    xi = 0   gas-liquid interface, u = 1                    (BSL Eq. 18.5-9/13)
    xi = XI  no-flux wall (parabolic=True, XI=1)            (BSL Eq. 18.5-10)
             or a truncation far from the front (plug flow) (BSL Eq. 18.5-14)

    parabolic=True   cap = 1 - xi^2   -> BSL Eq. 18.5-7  (the FULL problem)
    parabolic=False  cap = 1          -> BSL Eq. 18.5-11 (the penetration form)
    """

    def __init__(self, Lam, XI=1.0, n_x=400, stretch=6.0, parabolic=True,
                 wall="noflux", cap_override=None):
        u = np.linspace(0.0, 1.0, n_x + 1)
        # Smooth exponential stretch, fine at the interface where the front is.
        # Fixed map of u, so refining n_x refines the WHOLE grid in proportion.
        self.x_f = (XI * np.expm1(stretch * u) / np.expm1(stretch)
                    if stretch else XI * u)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.shape = (n_x, 1)              # (space, field) - never a bare (n,)
        self.V = np.diff(self.x_f)
        self.n_x, self.Lam, self.XI = n_x, Lam, XI
        if cap_override is not None:
            self.cap = np.full(n_x, float(cap_override))
        else:
            self.cap = (1.0 - self.x_c ** 2) if parabolic else np.ones(n_x)

        # OUTWARD normal throughout: a dc/dn + b c = d.
        #   xi = 0   (n points in -xi): Dirichlet u = 1  -> a=0, b=1, d=1
        #   xi = XI  (n points in +xi): no flux          -> a=1, b=0, d=0
        #                               or u = 0 (sink)  -> a=0, b=1, d=0
        left = {"a": 0.0, "b": 1.0, "d": 1.0}
        right = ({"a": 1.0, "b": 0.0, "d": 0.0} if wall == "noflux"
                 else {"a": 0.0, "b": 1.0, "d": 0.0})
        self.div = construct_div(self.shape, self.x_f, nu=0)   # nu=0: Cartesian
        self.grad, self.grad_bc = construct_grad(self.shape, self.x_f, self.x_c,
                                                 (left, right))
        self.A = (self.div @ self.grad).tocsr() * Lam
        self.b = np.asarray((self.div @ self.grad_bc).todense()).ravel() * Lam
        self.C = eye_array(n_x, format="csr").multiply(self.cap[:, None]).tocsr()

    def march(self, n_t=800):
        """Backward Euler on a zeta grid uniform in sqrt(zeta).

        The front advances as sqrt(zeta), so uniform steps in w = sqrt(zeta)
        put the steps where the solution moves. Returns (u at the outlet,
        mixing-cup uptake, marching integral of the interfacial flux)."""
        w = np.linspace(0.0, 1.0, n_t + 1)
        z = w ** 2
        u = np.zeros(self.n_x)
        flux_int = 0.0
        for j in range(1, n_t + 1):
            dz = z[j] - z[j - 1]
            Mx = (self.C / dz - self.A).tocsr()
            rhs = (self.C @ u) / dz + self.b
            u = solve_banded((1, 1), _banded(Mx), rhs)
            # interfacial flux from the same boundary machinery that built A
            flux_int += dz * (-self.Lam
                              * (self.grad @ u.reshape(-1, 1) + self.grad_bc)[0, 0])
        return u, float((self.V * self.cap * u).sum()), float(flux_int)


def phi_num(Lam, n_t=800, richardson=True, **kw):
    """Mixing-cup uptake Phi(Lam). Backward Euler is first order in the marching
    step, so the production value is Richardson-extrapolated: 2*Phi(2n) - Phi(n)."""
    _, a, _ = Film(Lam, **kw).march(n_t)
    if not richardson:
        return a
    _, b, _ = Film(Lam, **kw).march(2 * n_t)
    return 2 * b - a


def phi_penetration(Lam):
    """BSL Eq. 18.5-16 integrated across the film: the exact plug-flow,
    semi-infinite answer, int_0^inf erfc(xi/sqrt(4 Lam)) dxi."""
    return np.sqrt(4.0 * Lam / np.pi)


def phi_wall_series(Lam, n_terms=20000):
    """EXACT series for plug flow with a no-flux wall at xi = 1: separation of
    variables on v = 1 - u with v(0)=0, dv/dxi(1)=0, v(xi,0)=1 gives
    v = sum_n (2/lam_n) sin(lam_n xi) e^{-Lam lam_n^2}, lam_n = (n+1/2) pi, so
    Phi = 1 - sum_n (2/lam_n^2) e^{-Lam lam_n^2}. Shares no code with the solver."""
    lam = (np.arange(n_terms) + 0.5) * np.pi
    return 1.0 - float((2.0 / lam ** 2 * np.exp(-Lam * lam ** 2)).sum())'''))

# ---------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. Problem 18A.4, two ways

The closed form is Eq. 18.5-18. The independent route never touches an error
function: it marches Eq. 18.5-11 on finite volumes and reads the outlet mixing
cup, which is BSL's Problem 18C.3 identity used as a *different* route to the
same number rather than as a check on itself.
"""))

cells.append(code(r'''R, L = P("film_R"), P("film_L")
v_avg, D_cl = P("film_v_avg"), P("film_D")
sol_g, M_Cl2 = P("film_solubility"), P("M_Cl2")
RHO_W = 1.000        # g/cm3, ASSUMED (not printed in the problem); see break table

W_wet = 2 * np.pi * R                       # wetted width of the column
v_max = 1.5 * v_avg                         # <v> = (2/3) v_max, from Eq. 18.5-1
cA0_cl = (sol_g / M_Cl2) / (100.0 / RHO_W)  # g-mol/cm3
t_exp_film = L / v_max
pen_depth_film = np.sqrt(4 * D_cl * t_exp_film)     # sqrt(4 D t_exp), cm
kL_film = 2 * np.sqrt(D_cl / (np.pi * t_exp_film))

# Eq. 18.5-18
WA_closed = W_wet * L * cA0_cl * np.sqrt(4 * D_cl * v_max / (np.pi * L))   # g-mol/s
WA_closed_hr = WA_closed * 3600.0

# Independent route: pymrm march of Eq. 18.5-11 (plug flow, deep domain).
# Lam is a free choice here because delta cancels out of W_A in this limit;
# XI = 10 sqrt(4 Lam) puts the truncation ten penetration depths away.
LAM_REF, NX_REF, NT_REF = 0.25, 600, 800
XI_REF = 10 * np.sqrt(4 * LAM_REF)
Phi_pymrm = phi_num(LAM_REF, n_t=NT_REF, XI=XI_REF, n_x=NX_REF, stretch=4.0,
                    parabolic=False)
delta_ref = np.sqrt(D_cl * L / (v_max * LAM_REF))
WA_pymrm_hr = W_wet * delta_ref * v_max * cA0_cl * Phi_pymrm * 3600.0

WA_printed = P("film_WA_printed")
M["eq18518_vs_printed_18A4_rel"] = WA_closed_hr / WA_printed - 1
M["pymrm_vs_printed_18A4_rel"] = WA_pymrm_hr / WA_printed - 1
M["pymrm_vs_eq18518_rel"] = WA_pymrm_hr / WA_closed_hr - 1

# Does the closed form actually round to the printed three figures? Three printed
# decimals allow +/- 0.0005, and the page must not claim more agreement than that.
# The two unprinted inputs are re-evaluated here, not asserted about.
WA_round3 = round(WA_closed_hr, 3)
round_half_rel = 0.0005 / WA_printed
_pref = W_wet * L * np.sqrt(4 * D_cl * v_max / (np.pi * L)) * 3600.0
WA_rho999 = _pref * (sol_g / M_Cl2) / (100.0 / 0.999)
WA_mw71 = _pref * (sol_g / 71.0) / (100.0 / RHO_W)

# The PRINTED inputs are three-figure too, so measure what one of them is worth:
# half a unit in the last printed place of D, through the same closed form.
def _WA_of_D(d):                       # Eq. 18.5-18, everything else held fixed
    return W_wet * L * cA0_cl * np.sqrt(4 * d * v_max / (np.pi * L)) * 3600.0

D_half = 0.5 * 10.0 ** (np.floor(np.log10(D_cl)) - 2)   # D is printed to 3 figures
D_window_rel = max(abs(_WA_of_D(D_cl + D_half) / WA_closed_hr - 1),
                   abs(_WA_of_D(D_cl - D_half) / WA_closed_hr - 1))
WA_excess = WA_closed_hr - (WA_printed + 0.0005)        # how far outside the window

display(Markdown(rf"""
| route | absorption rate (g-mol/hr) | vs BSL's printed answer |
|---|---|---|
| Eq. 18.5-18, closed form | {WA_closed_hr:.5f} | {M["eq18518_vs_printed_18A4_rel"]*100:+.3f} % |
| pymrm march of Eq. 18.5-11, outlet mixing cup | {WA_pymrm_hr:.5f} | {M["pymrm_vs_printed_18A4_rel"]*100:+.3f} % |
| **BSL, printed** | **{WA_printed:.3f}** | - |

The two routes share no assembly *where it matters* - one evaluates a square root,
the other solves {NX_REF} cells over {2*NT_REF} marching steps and integrates the
outlet profile - and they agree to
{abs(M["pymrm_vs_eq18518_rel"])*100:.2e} %. Both sit
{M["eq18518_vs_printed_18A4_rel"]*100:+.2f} % from the printed answer.

**That residual is a hair *outside* three-figure rounding, not inside it, and it
is worth saying so.** {WA_closed_hr:.5f} rounds to **{WA_round3:.3f}**, whereas
BSL prints **{WA_printed:.3f}**; three printed decimals allow
$\pm${round_half_rel*100:.3f} % and the deviation is
{M["eq18518_vs_printed_18A4_rel"]*100:+.3f} %. Two inputs Problem 18A.4 does
**not** print are each worth about that much. Taking the density of water as
0.999 g/cm$^3$ instead of the assumed {RHO_W:.3f} gives {WA_rho999:.5f}
({(WA_rho999/WA_printed-1)*100:+.3f} %); taking $M$(Cl$_2$) as 71.0 instead of
Table E.1's {M_Cl2:.3f} gives {WA_mw71:.5f}
({(WA_mw71/WA_printed-1)*100:+.3f} %); **each of those rounds to the printed
{WA_printed:.3f}**.

**It is not the size of those two alone, though, and the excess is thin.** The
amount to be explained is only {WA_excess:.2e} g-mol/hr,
{WA_excess/WA_closed_hr*100:.4f} % of the value - and the inputs 18A.4 *does*
print are three-figure as well. Half a unit in the last printed place of
$\mathscr{{D}}_{{AB}}$ moves $W_A$ by
$\pm${D_window_rel*100:.3f} % on its own, which is the whole
$\pm${round_half_rel*100:.3f} % window; and if BSL truncated to three decimals
rather than rounding, {WA_closed_hr:.5f} sits inside
[{WA_printed:.3f}, {WA_printed+0.001:.3f}) and there is nothing to explain at all.
So the page does not conclude that BSL used either unprinted value, and does not
claim the residual belongs to them specifically - it concludes that a three-figure
printed answer cannot tell, which is why both are break-table rows and why neither
is adopted.

**Be exact about which half of the calculation each number tests**, because the
two routes are *not* independent end to end. The wetted width $W=2\pi R$, the
factor $\tfrac32$, the solubility conversion and the 3600 are computed once and
enter both; and the reference thickness
$\delta_{{\rm ref}}=\sqrt{{\mathscr{{D}}L/(v_{{\max}}\Lambda)}}$ carries the same
dimensional group back out of the march, so $W$, $c_{{A0}}$, $v_{{\max}}$,
$\mathscr{{D}}$, $L$ and the 3600 all cancel identically out of the ratio.
`pymrm_vs_eq18518_rel` therefore **cannot see any of them**: the break table
injects $v_{{\max}}=\langle v\rangle$, $\rho=0.999$ and $M$(Cl$_2$)$=71.0$ into it
and it does not move. What the {abs(M["pymrm_vs_eq18518_rel"])*100:.2e} % *does*
test is the dimensionless solve - that the marched $\Phi(\Lambda)$ equals
$\sqrt{{4\Lambda/\pi}}$, an error function against {NX_REF} finite volumes, with
the 4 and the $\pi$ in it. Drop the 4 from Eq. 18.5-18 and this metric moves; that
row is in the break table too. **The dimensional prefactor is tested by the other
column of the table above**: the {M["eq18518_vs_printed_18A4_rel"]*100:+.3f} %
against BSL's printed {WA_printed:.3f} g-mol/hr is what says $W$, the $\tfrac32$,
the solubility conversion and the units are right, and it is the only thing on the
page that says so.

Exposure time $t_{{\rm exp}} = L/v_{{\max}}$ = **{t_exp_film:.4f} s**, penetration
depth $\sqrt{{4\mathscr{{D}}t_{{\rm exp}}}}$ = **{pen_depth_film*1e4:.1f} um**,
$k_L = 2\sqrt{{\mathscr{{D}}/\pi t_{{\rm exp}}}}$ = **{kL_film:.4e} cm/s**.
"""))'''))

cells.append(md(r"""### 2. Problem 18A.7(a): the instruction and the answer do not match

Part (a) says, in as many words, "Use Eq. 18.5-20" - the **creeping-flow** result
with $3\pi D$ in the denominator. Evaluating both printed equations on the
problem's own printed inputs settles which one produced the printed answer.
"""))

cells.append(code(r'''D_bub, DAB, cA0_bub, vt = (P("bub_D_bubble"), P("bub_DAB"),
                           P("bub_cA0") * 1e-3, P("bub_vt"))   # liter -> cm3
Na_a, Na_b = P("bub_Na_printed_a"), P("bub_Na_printed_b")
kc_meas = P("bub_kc_measured") / 3600.0                        # cm/hr -> cm/s

kL_19 = np.sqrt(4 * DAB * vt / (np.pi * D_bub))        # Eq. 18.5-19, potential flow
kL_20 = np.sqrt(4 * DAB * vt / (3 * np.pi * D_bub))    # Eq. 18.5-20, creeping flow
N19, N20, Nmeas = kL_19 * cA0_bub, kL_20 * cA0_bub, kc_meas * cA0_bub

M["eq18519_vs_printed_18A7a_rel"] = N19 / Na_a - 1
M["eq18520_vs_printed_18A7a_rel"] = N20 / Na_a - 1
M["eq18519_over_eq18520"] = kL_19 / kL_20
M["hg_measured_vs_printed_18A7b_rel"] = Nmeas / Na_b - 1

display(Markdown(rf"""
| equation | $(N_A)_{{\rm avg}}$ (g-mol cm$^{{-2}}$ s$^{{-1}}$) | vs the printed answer {Na_a:.2e} |
|---|---|---|
| Eq. 18.5-19 (potential flow, $\pi D$) | {N19:.4e} | {M["eq18519_vs_printed_18A7a_rel"]*100:+.3f} % |
| Eq. 18.5-20 (creeping flow, $3\pi D$) - **the one the problem names** | {N20:.4e} | {M["eq18520_vs_printed_18A7a_rel"]*100:+.2f} % |

**The printed answer to part (a) is Eq. 18.5-19's, not Eq. 18.5-20's.** The two
differ by exactly $\sqrt3$ = {M["eq18519_over_eq18520"]:.7f}, which is what
{N19:.4e}/{N20:.4e} comes to; nothing about this depends on a transcription, since
both equations and both numbers are printed on pages 561 and 570 of the same book.
This is reported, not repaired: the page does not decide whether the instruction or
the answer is the error, and every result below that uses a bubble flux states which
equation it came from.

Part (b) recomputes the rate from Hammerton & Garner's measured surface-averaged
$k_c$ = {P("bub_kc_measured"):.0f} cm/hr = {kc_meas:.6f} cm/s, giving
{Nmeas:.4e}, or {M["hg_measured_vs_printed_18A7b_rel"]*100:+.3f} % of the printed
{Na_b:.2e} - so part (b) reproduces cleanly and pins the units of the measurement.
"""))'''))

cells.append(md(r"""### 3. The only empirical test on disk, and how big it is

The penetration prediction for this bubble - BSL's Example 18.5-1 evaluated on
Problem 18A.7's printed inputs - carries **no adjustable constant**: the exposure
time is $D/v_t$, both printed. Put it against the measurement.
"""))

cells.append(code(r'''t_exp_bub = D_bub / vt
kL_higbie_hr, kL_creep_hr = kL_19 * 3600, kL_20 * 3600
kc_meas_hr = P("bub_kc_measured")

M["higbie_kL_vs_hg_measured_rel"] = kL_19 / kc_meas - 1
M["creeping_kL_vs_hg_measured_rel"] = kL_20 / kc_meas - 1

# What the other two pictures need in order to fit the same single datum.
delta_film_um = DAB / kc_meas * 1e4          # film theory: k_L = D/delta
s_renewal = kc_meas ** 2 / DAB               # surface renewal: k_L = sqrt(D s)
M["film_thickness_implied_um"] = delta_film_um
M["renewal_rate_implied_s_inv"] = s_renewal
M["higbie_texp_times_renewal_rate"] = s_renewal * t_exp_bub
# sqrt(D s)/(2 sqrt(D/(pi t))) - the ratio the two pictures would have to satisfy
M["renewal_over_higbie_kL_at_measured"] = np.sqrt(np.pi * s_renewal * t_exp_bub) / 2

display(Markdown(rf"""
Exposure time $t_{{\rm exp}} = D/v_t$ = **{t_exp_bub:.5f} s** - printed geometry
over printed velocity, nothing fitted.

| picture | free constant | $k_L$ (cm/hr) | vs measured {kc_meas_hr:.0f} cm/hr |
|---|---|---|---|
| penetration, Eq. 18.5-19 (potential flow) | **none** | {kL_higbie_hr:.2f} | {M["higbie_kL_vs_hg_measured_rel"]*100:+.1f} % |
| penetration, Eq. 18.5-20 (creeping flow) | **none** | {kL_creep_hr:.2f} | {M["creeping_kL_vs_hg_measured_rel"]*100:+.1f} % |
| Whitman film, $k_L=\mathscr{{D}}/\delta$ | $\delta$ | fits exactly | 0 by construction |
| Danckwerts renewal, $k_L=\sqrt{{\mathscr{{D}}s}}$ | $s$ | fits exactly | 0 by construction |

**This is the whole of the three-way comparison that a single measurement can
support, and it is worth being exact about what it says.** It does *not* say that
penetration theory is right and the other two wrong. It says that on this datum
only penetration makes a prediction at all: the film and renewal pictures each have
one free constant and there is one number to fit, so each reproduces it with zero
residual whatever it had been. That is `A3.3`'s argument, and this page inherits
it rather than re-deriving it. What is new is the other half - **the model with no
free constant is {abs(M["higbie_kL_vs_hg_measured_rel"])*100:.0f} % low**, and BSL's
own words for that agreement are "approximately confirmed".

Inverting the same measurement through the other two pictures gives the constants
they would need: a film **{delta_film_um:.2f} um** thick, or a renewal rate
**{s_renewal:.1f} s$^{{-1}}$** (mean surface age {1/s_renewal*1e3:.1f} ms). Those are
not independent facts about the liquid; they are the measurement, rewritten. The
renewal rate and the penetration exposure time come out at
$s\,t_{{\rm exp}}$ = {M["higbie_texp_times_renewal_rate"]:.4f}, and the two pictures
would predict the same $k_L$ when $\sqrt{{\pi s t_{{\rm exp}}}}/2 = 1$; the value
here is **{M["renewal_over_higbie_kL_at_measured"]:.4f}**, i.e. the surface-renewal
fit sits {(M["renewal_over_higbie_kL_at_measured"]-1)*100:.1f} % above the
parameter-free penetration prediction. That is not a second finding: it is the
{M["higbie_kL_vs_hg_measured_rel"]*100:+.1f} % of the table above, read backwards
($1/{M["renewal_over_higbie_kL_at_measured"]:.4f} - 1 =
{(1/M["renewal_over_higbie_kL_at_measured"]-1)*100:+.1f}$ %).
"""))'''))

cells.append(md(r"""### 4. Where the exponent question actually stands

The discriminating claim of penetration theory - as of surface renewal - is
$k_L\propto\mathscr{D}^{1/2}$ against the film's $\mathscr{D}^{1}$. **Nothing on
disk tests it for this case, and the reason is not "no data" but "one datum".**
Here is the search behind that statement, because a negative claim needs one:

- BSL section 18.5 (book pp. 558-561) carries **no numerical result**; the only
  figures in it are the "0.3 to 0.5 cm" bubble range over which Eq. 18.5-19 is
  said to have been approximately confirmed, and the competing exponents. Read
  in full at 300 ppi.
- **Its problems, all of them.** Every problem statement in Chapter 18 - 18A.1
  through 18D.2, book pp. 568-581 - was read on renders. Five touch
  Eqs. 18.5-16 to 20. **18A.3** propagates a $\pm5$ % solubility and $\pm15$ %
  diffusivity uncertainty through Eq. 18.5-18: a sensitivity question, not a
  test. **18A.4** prints one diffusivity and no measured rate at all - its
  "answer" is Eq. 18.5-18's own output. **18A.7** prints one diffusivity and the
  single measured $k_c$ used above. **18C.3** asks for Eq. 18.5-18 to be
  re-derived by integrating $v_{\max}c_A$ across the outlet - the identity this
  page's numerical route uses - and carries no numbers. **18D.2** asks for the
  analogue of Eq. 18.5-18 when $A$ reacts irreversibly in the film; its printed
  Answer is an algebraic expression, not a number. So exactly one measurement, at
  one diffusivity.
- BSL Table 22.8-2, headed *Comparison of Film and Penetration Models*, is
  reproduced below. It compares the models on the **high-net-flux correction**
  $\theta_x$, not on $\mathscr{D}^{n}$, and it contains no measurement.
- Froment section 6.4 is theory throughout; its Table 6.4.2.1 (reproduced below)
  is model-versus-model. Section 6.5 describes *how* $k_L$, $s$ and the
  interfacial area are measured - including Danckwerts et al.'s plot of
  $(N_AA_v/C_{Ai})^2$ against $k$, whose intercept gives $s$ - **but prints no
  measured values**.

One $(k_L, \mathscr{D})$ pair cannot resolve an exponent, so the honest reading is
that the exponent question is still exactly where `A3.3` left it, and this page
does not claim to have moved it. What has moved is the *magnitude* question, which
`A3.1` and `A3.3` could not touch at all.
"""))

cells.append(md(r"""### 5. The one printed table that carries all three models

Froment's Table 6.4.2.1 tabulates the film, surface-renewal and penetration
enhancement factors for a pseudo-first-order reaction at six values of
$\gamma=\sqrt{\mathscr{D}_Ak}/k_L$, each column headed by its own closed form. The
three closed forms are recomputed here and compared with the printed cells.
"""))

cells.append(code(r'''def FA_film(g):        return g / np.tanh(g)                    # Froment 6.3.2-11
def FA_renewal(g):     return np.sqrt(1.0 + g * g)                # Froment 6.4.2-7
def FA_penetration(g):                                            # Froment 6.4.2-8
    return g * ((1 + np.pi / (8 * g * g)) * erf(2 / np.sqrt(np.pi) * g)
                + 1 / (2 * g) * np.exp(-4 / np.pi * g * g))

# The film and renewal columns are printed to two decimals, so "outside the
# bracket" means outside by more than one unit in that last printed place.
BRACKET_TOL = 0.01

def bracket_violations(pen_printed, film_printed, renew_printed):
    lo = np.minimum(film_printed, renew_printed) - BRACKET_TOL
    hi = np.maximum(film_printed, renew_printed) + BRACKET_TOL
    return (pen_printed < lo) | (pen_printed > hi)

gam = fro["gamma"].to_numpy()
t3 = fro.copy()
t3["film_computed"] = FA_film(gam)
t3["renewal_computed"] = FA_renewal(gam)
t3["penetration_computed"] = FA_penetration(gam)
for col, pcol in (("film", "FA_film_printed"),
                  ("renewal", "FA_surface_renewal_printed"),
                  ("penetration", "FA_penetration_printed")):
    t3[f"{col}_rel"] = t3[f"{col}_computed"] / t3[pcol] - 1

t3["printed_pen_outside_printed_bracket"] = bracket_violations(
    t3["FA_penetration_printed"].to_numpy(), t3["FA_film_printed"].to_numpy(),
    t3["FA_surface_renewal_printed"].to_numpy())

M["froment_t6421_film_worst_rel"] = float(t3["film_rel"].abs().max())
M["froment_t6421_renewal_worst_rel"] = float(t3["renewal_rel"].abs().max())
M["froment_t6421_penetration_worst_rel"] = float(t3["penetration_rel"].abs().max())
M["froment_t6421_pen_at_gamma_10_rel"] = float(
    t3.loc[t3["gamma"] == 10, "penetration_rel"].iloc[0])
M["froment_t6421_pen_bracket_violations"] = float(
    t3["printed_pen_outside_printed_bracket"].sum())

display(t3[["gamma", "FA_film_printed", "film_computed", "film_rel",
            "FA_surface_renewal_printed", "renewal_computed", "renewal_rel",
            "FA_penetration_printed", "penetration_computed", "penetration_rel",
            "printed_pen_outside_printed_bracket"]].round(6))'''))

cells.append(code(r'''g = np.logspace(-2, 1.2, 400)
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
for a in ax:
    a.plot(g, FA_film(g), label=r"film, $\gamma/\tanh\gamma$")
    a.plot(g, FA_penetration(g), label="penetration (Higbie), Eq. 6.4.2-8")
    a.plot(g, FA_renewal(g), label=r"surface renewal, $\sqrt{1+\gamma^2}$")
    a.plot(t3["gamma"], t3["FA_penetration_printed"], "kx", ms=7,
           label="Table 6.4.2.1, penetration column as printed")
    a.set_xscale("log")
    a.set_xlabel(r"$\gamma$")
ax[0].set_yscale("log")
ax[0].set_ylabel(r"$F_A$")
ax[0].set_title("all six printed rows")
ax[1].set_xlim(5e-3, 0.6)
ax[1].set_ylim(0.9, 1.12)
ax[1].set_title(r"the small-$\gamma$ end, where $F_A \to 1$")
ax[0].legend(fontsize=7, loc="upper left")
fig.tight_layout()
plt.show()

display(Markdown(rf"""
**The surface-renewal column reproduces $\sqrt{{1+\gamma^2}}$ at all six
$\gamma$** (worst {M["froment_t6421_renewal_worst_rel"]*100:.2f} %, which is
rounding). **The film column is high by one unit in the last printed place at two
of six** ($\gamma=0.3$: 1.04 printed for {FA_film(0.3):.4f}; $\gamma=3$: 3.02
printed for {FA_film(3.0):.4f}), worst
{M["froment_t6421_film_worst_rel"]*100:.2f} %. **The penetration column disagrees
with the closed form printed at the head of its own column at every one of the
six**, worst {M["froment_t6421_penetration_worst_rel"]*100:.1f} %.

Two of those cells are wrong on the table's own evidence, with no computation
needed:

- at $\gamma=0.01$ the printed value is **0.94**, while the other two columns of
  the same row both read 1.00. An enhancement factor below 1 would mean a reaction
  that *suppresses* absorption;
- at $\gamma=10$ the printed value is **10.39**, above the surface-renewal
  column's 10.05 in the same row, while Froment's own text two paragraphs earlier
  says the three expressions "differ only by a few percent". The closed form gives
  {FA_penetration(10.0):.4f}, i.e. {M["froment_t6421_pen_at_gamma_10_rel"]*100:+.1f} %
  from what is printed, and it lies between the other two columns as it should.

`froment_t6421_pen_bracket_violations` =
{M["froment_t6421_pen_bracket_violations"]:.0f} counts exactly those two. The test
is: is the printed penetration cell outside the interval spanned by the *printed*
film and surface-renewal cells of the same row, by more than one unit
({BRACKET_TOL}) in the last decimal place those two columns are printed to? It uses
**the printed digits alone** and none of the closed forms.

Froment prints the table "After Beek [1968]"; Beek 1968 is not on disk and was not
consulted, so this is a statement about Froment's printing and nothing more.
**Nothing is repaired**: the CSV holds the cells as printed, and the curves above
are the closed forms, not the cells.
"""))'''))

cells.append(md(r"""### 6. BSL's own film-versus-penetration comparison, on a different axis

Chapter 22 compares the same two pictures on something else entirely: the
correction factor $\theta_x$ for a **finite net mass-transfer rate** across the
interface. Eq. 22.8-41 gives $\theta_x=(1-x_{A0})\psi(x_{A0})$ for the penetration
model, with $\psi$ from Arnold's problem (section 20.1, Table 20.1-1), and
Eq. 22.8-42 gives $\theta_x=\frac{1-x_{A0}}{x_{A0}}\ln\frac{1}{1-x_{A0}}$ for the
film. $\varphi$ is obtained here by root-finding BSL's Eq. 20.1-18, never by
reading the printed $\varphi$ back in, so both printed columns of Table 20.1-1 are
independent targets.
"""))

cells.append(code(r'''def xA0_of_phi(p):                                    # BSL Eq. 20.1-18
    return 1.0 / (1.0 + 1.0 / (np.sqrt(np.pi) * (1 + erf(p)) * p * np.exp(p * p)))

def phi_of_xA0(x, f=xA0_of_phi):
    return brentq(lambda p: f(p) - x, 1e-12, 10.0, xtol=1e-15, rtol=8.9e-16)

tt = tabs[(tabs["x_A0"] > 0) & (tabs["x_A0"] < 1)].copy()
xs = tt["x_A0"].to_numpy()
tt["phi_computed"] = [phi_of_xA0(x) for x in xs]
tt["psi_computed"] = tt["phi_computed"] * np.sqrt(np.pi) / xs        # table head
tt["theta_pen_computed"] = (1 - xs) * tt["psi_computed"]             # Eq. 22.8-41
tt["theta_film_computed"] = (1 - xs) / xs * np.log(1 / (1 - xs))     # Eq. 22.8-42
for a, b in (("phi", "phi_20p1_1"), ("psi", "psi_20p1_1"),
             ("theta_pen", "theta_penetration_22p8_2"),
             ("theta_film", "theta_film_22p8_2")):
    tt[f"{a}_rel"] = tt[f"{a}_computed"] / tt[b] - 1

M["bsl_t2011_phi_worst_rel"] = float(tt["phi_rel"].abs().max())
M["bsl_t2011_psi_worst_rel"] = float(tt["psi_rel"].abs().max())
M["bsl_t2282_theta_pen_worst_rel"] = float(tt["theta_pen_rel"].abs().max())
M["bsl_t2282_theta_film_worst_rel"] = float(tt["theta_film_rel"].abs().max())
M["bsl_t2282_theta_pen_minus_psi_dev"] = float(
    (tt["theta_pen_rel"] - tt["psi_rel"]).abs().max())
M["bsl_t2282_pen_film_gap_at_0p75"] = float(
    1 - tt.loc[tt["x_A0"] == 0.75, "theta_pen_computed"].iloc[0]
      / tt.loc[tt["x_A0"] == 0.75, "theta_film_computed"].iloc[0])

display(tt[["x_A0", "phi_20p1_1", "phi_computed", "phi_rel",
            "psi_20p1_1", "psi_computed", "psi_rel",
            "theta_penetration_22p8_2", "theta_pen_computed", "theta_pen_rel",
            "theta_film_22p8_2", "theta_film_computed", "theta_film_rel"]].round(6))

display(Markdown(rf"""
All four printed columns reproduce inside their own last printed digit: worst
{M["bsl_t2011_phi_worst_rel"]*100:.3f} % on $\varphi$,
{M["bsl_t2011_psi_worst_rel"]*100:.3f} % on $\psi$,
{M["bsl_t2282_theta_pen_worst_rel"]*100:.3f} % on the penetration $\theta_x$ and
{M["bsl_t2282_theta_film_worst_rel"]*100:.3f} % on the film $\theta_x$. Unlike
Froment's Table 6.4.2.1, **BSL's two tables are internally consistent throughout**.

**But the penetration $\theta_x$ is not an independent check.** Eq. 22.8-41 is
$(1-x_{{A0}})\psi$, so the row-by-row deviations of the $\theta_x$ column and of
the $\psi$ column are the same number: the worst difference between the two sets of
deviations is {M["bsl_t2282_theta_pen_minus_psi_dev"]:.2e}. One of the two metrics
is redundant, and the notebook says so rather than counting it twice.

The physics BSL draws from the table: "the penetration model predicts a stronger
correction $\theta_x$ for net mass transfer than does the film model. This is in
part because the net flow thickens the boundary layer, an effect that the film
model does not consider." The gap reaches
{M["bsl_t2282_pen_film_gap_at_0p75"]*100:.1f} % at $x_{{A0}}=0.75$. Note what this
does and does not say: the two rows at $x_{{A0}}=0$ and $x_{{A0}}=1$ are the two
analytic limits of both models and agree by construction, so only the four interior
rows carry information; and $\theta_x$ is a high-flux correction, **not** the
$\mathscr{{D}}^{{n}}$ question.
"""))'''))

cells.append(md(r"""### 7. When is the penetration approximation actually good?

BSL replaces Eq. 18.5-7 by Eq. 18.5-11 with two arguments at once: $A$ "has the
impression" the film moves at $v_{\max}$, and $A$ does not "sense" the wall. The
solver can switch those separately, which turns the sentence into two numbers.
Write

$$r_{\rm full}(\Lambda)=\frac{\Phi_{\rm full}}{\Phi_{\rm pen}},\qquad
r_{\rm wall}(\Lambda)=\frac{\Phi_{\rm plug,\,wall}}{\Phi_{\rm pen}},\qquad
r_{\rm profile}(\Lambda)=\frac{\Phi_{\rm full}}{\Phi_{\rm plug,\,wall}},$$

so that $r_{\rm full}=r_{\rm wall}\,r_{\rm profile}$ exactly. $r_{\rm wall}$ is
taken from the **exact eigenfunction series** of *Validation 2* rather than from
the solver, for a reason that is itself part of the answer: the wall effect turns
out to be smaller than the solver's own discretisation error, so the solver cannot
resolve it. $r_{\rm profile}$ is a ratio of two solves on the same grid with the
same step, so its discretisation error very largely cancels.
"""))

cells.append(code(r'''NX, NT, STR = 600, 800, 6.0

def r_profile(Lam):
    """Phi(full) / Phi(plug flow, same wall, same grid, same step)."""
    Ff = phi_num(Lam, n_t=NT, XI=1.0, n_x=NX, stretch=STR, parabolic=True)
    Fa = phi_num(Lam, n_t=NT, XI=1.0, n_x=NX, stretch=STR, parabolic=False)
    return Ff / Fa

def r_wall(Lam):
    """Exact: series answer with the wall, over the semi-infinite answer."""
    return phi_wall_series(Lam) / phi_penetration(Lam)

def r_full(Lam):
    return r_wall(Lam) * r_profile(Lam)

lams = np.logspace(-4, -0.7, 12)
rp = np.array([r_profile(l) for l in lams])
rw = np.array([r_wall(l) for l in lams])
rf = rw * rp

L_MARK = 0.1
M["profile_effect_at_lambda_0p1"] = 1 - r_profile(L_MARK)
M["wall_effect_series_at_lambda_0p1"] = 1 - r_wall(L_MARK)
M["profile_over_wall_effect_at_lambda_0p1"] = (
    M["profile_effect_at_lambda_0p1"] / M["wall_effect_series_at_lambda_0p1"])
# What the solver alone would say about the wall - kept to show it CANNOT say it.
Fa_solver = phi_num(L_MARK, n_t=NT, XI=1.0, n_x=NX, stretch=STR, parabolic=False)
M["wall_effect_solver_at_lambda_0p1"] = 1 - Fa_solver / phi_penetration(L_MARK)

# Root-find the 1 % threshold; never read it off the sweep above.
def lam_at_1pct(f=r_full, lo=1e-3, hi=0.5):
    return brentq(lambda l: f(l) - 0.99, lo, hi, xtol=1e-9)

lam_star = lam_at_1pct()
M["lambda_star_1pct"] = lam_star
M["pen_depth_over_delta_at_1pct"] = np.sqrt(4 * lam_star)
delta_star = np.sqrt(D_cl * L / (v_max * lam_star))
M["delta_star_18A4_um"] = delta_star * 1e4
# The same threshold read off a sweep by interpolation, for the break table:
# a sampled crossing is not a root, and how far off it is depends on the grid.
lam_star_swept = float(np.interp(-0.99, -rf, lams))
lams_coarse = np.logspace(-4, -0.7, 5)
rf_coarse = np.array([r_full(l) for l in lams_coarse])
lam_star_swept_coarse = float(np.interp(-0.99, -rf_coarse, lams_coarse))

fig, ax = plt.subplots(figsize=(5.8, 3.5))
ax.loglog(lams, 100 * (1 - rf), "o-", label="both (full Eq. 18.5-7)")
ax.loglog(lams, 100 * (1 - rp), "s--", label="parabolic velocity profile alone")
ax.loglog(lams, 100 * (1 - rw), "^:", label="finite film thickness alone (exact series)")
ax.axhline(1.0, color="0.6", lw=0.8)
ax.axvline(lam_star, color="0.6", lw=0.8)
ax.set_ylim(1e-12, 30)
ax.set_xlabel(r"$\Lambda = \mathscr{D}L/(v_{\max}\delta^2)$")
ax.set_ylabel("shortfall of Eq. 18.5-18, %")
ax.set_title("what the penetration approximation costs, and which half costs it")
ax.legend(fontsize=7)
fig.tight_layout()
plt.show()

display(Markdown(rf"""
**The two halves of BSL's argument are not equally load-bearing.** At
$\Lambda={L_MARK}$, where the penetration depth has reached
$\sqrt{{4\Lambda}}$ = {np.sqrt(4*L_MARK):.2f} of the film thickness, the finite
film costs {M["wall_effect_series_at_lambda_0p1"]*100:.5f} % while the parabolic
velocity profile costs {M["profile_effect_at_lambda_0p1"]*100:.2f} % - a factor
**{M["profile_over_wall_effect_at_lambda_0p1"]:.0f}** between them.

The reason is visible in the series: a reflecting wall only returns solute that a
semi-infinite film would have retained anyway, so what is lost is only the far tail
of the erfc profile - the part that would have reached beyond about *twice* the
film thickness - and the series puts that at
{M["wall_effect_series_at_lambda_0p1"]:.1e} here. What actually breaks
Eq. 18.5-18 is that the liquid *below* the surface moves more slowly than
$v_{{\max}}$ and therefore has longer to load up.

**The solver cannot see the wall effect at all**, and that is worth stating rather
than hiding: its own estimate,
{M["wall_effect_solver_at_lambda_0p1"]:.2e}, is a factor
{M["wall_effect_solver_at_lambda_0p1"]/M["wall_effect_series_at_lambda_0p1"]:.1f}
above the exact {M["wall_effect_series_at_lambda_0p1"]:.2e}, because the
discretisation error at these settings is of the same size. Being unable to resolve
an effect is the strongest available statement that it is negligible.

**The 1 % threshold, root-found rather than read off the sweep:**
$\Lambda^{{*}}$ = **{lam_star:.5f}**, i.e. Eq. 18.5-18 is within 1 % of the full
Eq. 18.5-7 for any film in which the penetration depth stays below
**{M["pen_depth_over_delta_at_1pct"]:.3f}$\,\delta$**. For Problem 18A.4's column
that means any water film thicker than **{delta_star*1e4:.0f} um**; the penetration
depth there is {pen_depth_film*1e4:.1f} um. Problem 18A.4 prints no film thickness,
no viscosity and no flow rate, so **this page does not state a correction for that
column** - it states the threshold, which is what the printed data support.

A sampled crossing is not a root, and how wrong it is depends entirely on the
sampling: reading the threshold off the twelve-point sweep plotted above gives
{lam_star_swept:.5f} ({abs(lam_star_swept/lam_star-1)*100:.1f} % away, which is
lucky), and off a five-point sweep over the same range
{lam_star_swept_coarse:.5f} ({abs(lam_star_swept_coarse/lam_star-1)*100:.0f} %
away). Both are in the break table; neither is what is reported.
"""))'''))

cells.append(code(r'''# Profiles: the full problem against Eq. 18.5-16, at a Lambda where they part.
LAM_P = 0.05
f_full = Film(LAM_P, XI=1.0, n_x=NX, stretch=STR, parabolic=True)
f_plug = Film(LAM_P, XI=1.0, n_x=NX, stretch=STR, parabolic=False)
u_full, cup_full, flux_full = f_full.march(1600)
u_plug, _, _ = f_plug.march(1600)
xi = f_full.x_c

fig, ax = plt.subplots(figsize=(5.4, 3.4))
ax.plot(xi, erfc(xi / np.sqrt(4 * LAM_P)), "k-", lw=2,
        label=r"Eq. 18.5-16, $\mathrm{erfc}$")
ax.plot(xi, u_plug, "--", label="pymrm, plug flow + no-flux wall")
ax.plot(xi, u_full, "-.", label="pymrm, full Eq. 18.5-7")
ax.set_xlabel(r"$\xi = x/\delta$")
ax.set_ylabel(r"$c_A/c_{A0}$ at $z = L$")
ax.set_title(rf"outlet profiles at $\Lambda$ = {LAM_P}")
ax.legend(fontsize=7)
fig.tight_layout()
plt.show()'''))

# ------------------------------------------------------------------ validation
cells.append(md(r"""## Validation

Five checks that can fail, one identity that cannot and is labelled as such, and a
break table.

### 1. The penetration branch against BSL's own closed forms

With the capacity set to 1 and the domain long, the solver is Eq. 18.5-11 and must
reproduce Eq. 18.5-16 pointwise and $\Phi=\sqrt{4\Lambda/\pi}$ - equivalently
Eq. 18.5-18 - in the integral. Neither shares any code with the solver.
"""))

cells.append(code(r'''LAM_V, NXV = 1e-2, 800
XIV = 12 * np.sqrt(4 * LAM_V)
fv = Film(LAM_V, XI=XIV, n_x=NXV, stretch=4.0, parabolic=False)
uv, Phi_raw, flux_raw = fv.march(1600)
_, Phi_half, _ = fv.march(800)
Phi_ext = 2 * Phi_raw - Phi_half

M["pen_profile_vs_eq18516_max_abs"] = float(
    np.abs(uv - erfc(fv.x_c / np.sqrt(4 * LAM_V))).max())
M["pen_uptake_vs_eq18518_rel"] = Phi_ext / phi_penetration(LAM_V) - 1
M["pen_uptake_vs_eq18518_raw_rel"] = Phi_raw / phi_penetration(LAM_V) - 1
M["mixing_cup_vs_flux_integral_rel"] = flux_raw / Phi_raw - 1

display(Markdown(rf"""
- outlet profile against Eq. 18.5-16, worst absolute deviation in $c_A/c_{{A0}}$:
  **{M["pen_profile_vs_eq18516_max_abs"]:.3e}**
- outlet uptake against $\sqrt{{4\Lambda/\pi}}$ (equivalently Eq. 18.5-18):
  **{M["pen_uptake_vs_eq18518_rel"]:+.3e}** extrapolated,
  {M["pen_uptake_vs_eq18518_raw_rel"]:+.3e} raw at $n_t$ = 1600
- outlet mixing cup against the marching integral of the surface flux
  (BSL's Problem 18C.3 identity): **{M["mixing_cup_vs_flux_integral_rel"]:+.3e}**

The third of those is **structural**: backward Euler on a conservative
finite-volume discretisation makes the accumulated surface flux equal the outlet
inventory to round-off, whatever the physics is, so it tests the assembly and the
boundary machinery and nothing else. It also sits **below** `ABS_FLOOR` = 1e-12, so
`check_agreement.py` will not compare it at all: it is not protected, it is merely
recorded. The two metrics above it are its above-floor companions and they are the
ones that can move.
"""))'''))

cells.append(md(r"""### 2. The wall branch against an exact series

Plug flow with a no-flux wall has a closed solution by separation of variables that
the books do not print, and which shares nothing with the solver: with $v=1-u$,
$v(0)=0$, $v_\xi(1)=0$, $v(\xi,0)=1$,

$$v = \sum_{n\ge0}\frac{2}{\lambda_n}\sin(\lambda_n\xi)\,e^{-\Lambda\lambda_n^{2}\zeta},
\qquad \lambda_n=(n+\tfrac12)\pi
\qquad\Longrightarrow\qquad
\Phi = 1-\sum_{n\ge0}\frac{2}{\lambda_n^{2}}e^{-\Lambda\lambda_n^{2}}.$$

Its $\Lambda\to0$ limit is $\sqrt{4\Lambda/\pi}$, which is a real check on both:
the series is built from sines and the closed form from an error function.
"""))

cells.append(code(r'''rows = []
for Lam in (1e-3, 1e-2, 1e-1):
    num = phi_num(Lam, n_t=NT, XI=1.0, n_x=NX, stretch=STR, parabolic=False)
    ser = phi_wall_series(Lam)
    rows.append((Lam, num, ser, num / ser - 1, 1 - ser / phi_penetration(Lam)))
wall = pd.DataFrame(rows, columns=["Lambda", "pymrm", "exact series",
                                   "rel dev", "wall effect (1 - series/pen)"])
M["wall_series_vs_pymrm_worst_rel"] = float(wall["rel dev"].abs().max())
display(wall)

display(Markdown(rf"""
Worst deviation between the finite-volume march and the exact series over three
decades of $\Lambda$: **{M["wall_series_vs_pymrm_worst_rel"]:.2e}** - the size of
the solver's marching error, not of the physics. The last column is why the wall
effect had to come from the series in *Results 7*: at $\Lambda=10^{{-3}}$ and
$10^{{-2}}$ it is at the double-precision floor, and at $\Lambda=0.1$ it is
{wall.iloc[-1, -1]:.2e}, still an order of magnitude below the solver's own error.
"""))'''))

cells.append(md(r"""### 3. Both axes refined, with observed orders

The marching step and the grid are refined separately, because the axis nobody
refines has repeatedly been the larger error. The two studies are deliberately
built so that each isolates its own axis: the marching study runs on a grid fine
enough that its space error is negligible and compares against the *exact* answer,
while the grid study holds the step fixed and compares against the finest grid at
that same step, so the marching error cancels out of the difference.
"""))

cells.append(code(r'''ref = phi_penetration(LAM_V)
NT_FIX, NX_FINE = 1600, 1600

tt_err = []
for n_t in (200, 400, 800, 1600):
    _, a, _ = Film(LAM_V, XI=XIV, n_x=3200, stretch=4.0, parabolic=False).march(n_t)
    tt_err.append((n_t, abs(a / ref - 1)))

_, a_fine, _ = Film(LAM_V, XI=XIV, n_x=NX_FINE, stretch=4.0,
                    parabolic=False).march(NT_FIX)
xx_err = []
for n_x in (100, 200, 400):
    _, a, _ = Film(LAM_V, XI=XIV, n_x=n_x, stretch=4.0,
                   parabolic=False).march(NT_FIX)
    xx_err.append((n_x, abs(a / a_fine - 1)))

def order(rows):
    (n1, e1), (n2, e2) = rows[-2], rows[-1]
    return float(np.log(e1 / e2) / np.log(n2 / n1))

M["time_order"] = order(tt_err)
M["space_order"] = order(xx_err)
M["time_err_n1600"] = float(tt_err[-1][1])
M["space_err_n400"] = float(xx_err[-1][1])
# Production runs n_x = 600, not the 400 the study's last row uses. Measure the
# space error there too, so the comparison below names settings that exist.
_, a_600, _ = Film(LAM_V, XI=XIV, n_x=NX, stretch=4.0,
                   parabolic=False).march(NT_FIX)
space_err_nx600 = abs(a_600 / a_fine - 1)

display(pd.DataFrame({"n_t": [r[0] for r in tt_err],
                      "rel err vs exact (n_x=3200)": [r[1] for r in tt_err]}))
display(pd.DataFrame({"n_x": [r[0] for r in xx_err],
                      f"rel diff vs n_x={NX_FINE} (n_t={NT_FIX})":
                          [r[1] for r in xx_err]}))
display(Markdown(rf"""
Observed orders: **{M["time_order"]:.3f}** in the marching step (backward Euler,
first order) and **{M["space_order"]:.3f}** in the grid (second order, on the
stretched mesh). **The marching step is the larger error**: at the finest settings
of these two studies, {M["time_err_n1600"]:.2e} at $n_t$ = {tt_err[-1][0]} (on
$n_x$ = 3200) against {M["space_err_n400"]:.2e} at $n_x$ = {xx_err[-1][0]} (at
$n_t$ = {NT_FIX}). Production is $n_x$ = {NX} with Richardson extrapolation over
$n_t$ = {NT}/{2*NT}; measured the same way, the space error at $n_x$ = {NX} is
smaller still, {space_err_nx600:.2e} - so the conclusion holds a fortiori. That is why every
$\Phi$ used for a physics statement above is Richardson-extrapolated in $n_t$ with
the raw value kept alongside - and why the wall effect, which is smaller than
either, is taken from the series instead.
"""))'''))

cells.append(md(r"""### 4. Froment's Eq. 6.4.2-8 re-derived from his Eq. 6.4.2-4

Froment prints both the instantaneous element flux for a pseudo-first-order
reaction at surface age $t$ (Eq. 6.4.2-4),

$$N_A(t)=\sqrt{kD_A}\,C_{Ai}\left[\operatorname{erf}\sqrt{kt}
      + \frac{e^{-kt}}{\sqrt{\pi kt}}\right],$$

and the enhancement factor Higbie's *uniform* age gives (Eq. 6.4.2-8). The second
must be the average of the first over $0<t<\bar t$, divided by
$k_LC_{Ai}=2C_{Ai}\sqrt{D_A/\pi\bar t}$ - and with $\gamma^2=D_Ak/k_L^2$ that means
$k\bar t=4\gamma^2/\pi$, which is exactly the argument of the exponential Froment
prints. Doing the quadrature checks both transcriptions at once, and it is what
makes the disagreement between Eq. 6.4.2-8 and the numbers printed underneath it a
defect in the table rather than in the reading.
"""))

cells.append(code(r'''def FA_penetration_quadrature(g):
    U = 4 * g * g / np.pi                       # = k tbar
    integrand = lambda v: erf(np.sqrt(v)) + np.exp(-v) / np.sqrt(np.pi * v)
    I, _ = quad(integrand, 0.0, U, limit=400)
    return (I / U) / (2.0 / np.sqrt(np.pi * U))

gg = np.array([0.05, 0.1, 0.3, 1.0, 2.0, 3.0, 5.0, 10.0])
qq = np.array([FA_penetration_quadrature(x) for x in gg])
cc = FA_penetration(gg)
M["froment_eq6428_closed_vs_quadrature_worst"] = float(np.abs(cc / qq - 1).max())
M["froment_penetration_vs_renewal_at_gamma_1"] = float(
    FA_penetration(1.0) / FA_renewal(1.0) - 1)
M["froment_penetration_vs_film_at_gamma_1"] = float(
    FA_penetration(1.0) / FA_film(1.0) - 1)

display(Markdown(rf"""
Worst relative deviation between Eq. 6.4.2-8 as printed and the quadrature of
Eq. 6.4.2-4 over eight $\gamma$ spanning {gg.min()} to {gg.max()}:
**{M["froment_eq6428_closed_vs_quadrature_worst"]:.2e}**. Both transcriptions are
therefore right.

**That number is pure round-off, and it is CI-brittle.** It sits at
{M["froment_eq6428_closed_vs_quadrature_worst"]/1e-12:.2f} x `ABS_FLOOR` = 1e-12,
so `check_agreement.py` *does* compare it - at a 5 % relative tolerance - but what
it is comparing is the accumulated float error of `scipy.integrate.quad` plus the
libm `erf`. Any scipy or libm change moves that by orders of magnitude and CI
reports a regression that is not one; and if a future build lands both sides under
1e-12, comparison stops silently instead. So it is carried by two companions well
above the floor, which are also the physics: at $\gamma=1$ the penetration
enhancement
sits **{M["froment_penetration_vs_renewal_at_gamma_1"]*100:+.2f} %** from surface
renewal's and **{M["froment_penetration_vs_film_at_gamma_1"]*100:+.2f} %** from the
film's. Those two numbers are the *whole* of the practical difference between the
three models on this axis, and they are the reason Froment concludes that "the
choice of the model matters little for design calculations: the predicted
differences are negligible with respect to the uncertainties of prediction of some
of the model or operating parameters".
"""))'''))

cells.append(md(r"""### 5. The break table

Every metric reported on this page needs something that moves it. Where nothing
does, the row is kept and the metric is labelled - an identity is worth having once
it is named as one. **This table was rebuilt for this page's physics; none of it
travelled from `A3.1` or `A3.3`.**
"""))

cells.append(code(r'''BR = []          # (defect, metric, baseline, broken, moved?)

def brk(label, **broken):
    for key, val in broken.items():
        base = M[key]
        denom = max(abs(base), 1e-30)
        moved = (abs(val - base) / denom > 0.10) or (np.sign(val) != np.sign(base))
        BR.append((label, key, base, val, "yes" if moved else "NO"))

def rel(a, b):
    return a / b - 1

def WA_hr(vmax=v_max, four=4.0, rho=RHO_W, mw=M_Cl2):
    c = (sol_g / mw) / (100.0 / rho)
    return W_wet * L * c * np.sqrt(four * D_cl * vmax / (np.pi * L)) * 3600

def WA_pymrm_hr_b(vmax=v_max, rho=RHO_W, mw=M_Cl2):
    """The SAME defects pushed through the marching route. delta_ref carries the
    dimensional group back out, so these cancel out of pymrm_vs_eq18518_rel - which
    is exactly what the paired rows below are here to show."""
    c = (sol_g / mw) / (100.0 / rho)
    d = np.sqrt(D_cl * L / (vmax * LAM_REF))
    return W_wet * d * vmax * c * Phi_pymrm * 3600.0

# --- the closed forms ---------------------------------------------------------
brk("v_max taken as the printed <v> (the 3/2 that follows from Eq. 18.5-1 forgotten)",
    eq18518_vs_printed_18A4_rel=rel(WA_hr(vmax=v_avg), WA_printed),
    pymrm_vs_eq18518_rel=rel(WA_pymrm_hr_b(vmax=v_avg), WA_hr(vmax=v_avg)))
brk("the 4 dropped from Eq. 18.5-18",
    eq18518_vs_printed_18A4_rel=rel(WA_hr(four=1.0), WA_printed),
    pymrm_vs_eq18518_rel=rel(WA_pymrm_hr_b(), WA_hr(four=1.0)))
brk("rho_water = 0.999 instead of the assumed 1.000 g/cm3",
    eq18518_vs_printed_18A4_rel=rel(WA_hr(rho=0.999), WA_printed),
    pymrm_vs_eq18518_rel=rel(WA_pymrm_hr_b(rho=0.999), WA_hr(rho=0.999)))
brk("M(Cl2) = 71.0 instead of BSL Table E.1's 70.905",
    eq18518_vs_printed_18A4_rel=rel(WA_hr(mw=71.0), WA_printed),
    pymrm_vs_eq18518_rel=rel(WA_pymrm_hr_b(mw=71.0), WA_hr(mw=71.0)))
brk("Eq. 18.5-20 used for 18A.7(a), as the problem statement instructs",
    eq18519_vs_printed_18A7a_rel=rel(N20, Na_a))
brk("the 3 dropped from Eq. 18.5-20, making it Eq. 18.5-19",
    eq18520_vs_printed_18A7a_rel=rel(N19, Na_a),
    eq18519_over_eq18520=1.0)
brk("k_c left in cm/hr (the /3600 forgotten)",
    hg_measured_vs_printed_18A7b_rel=rel(P("bub_kc_measured") * cA0_bub, Na_b),
    higbie_kL_vs_hg_measured_rel=rel(kL_19, P("bub_kc_measured")),
    creeping_kL_vs_hg_measured_rel=rel(kL_20, P("bub_kc_measured")),
    film_thickness_implied_um=DAB / P("bub_kc_measured") * 1e4,
    renewal_rate_implied_s_inv=P("bub_kc_measured") ** 2 / DAB)
brk("exposure time taken as 2D/v_t instead of D/v_t",
    higbie_kL_vs_hg_measured_rel=rel(2 * np.sqrt(DAB / (np.pi * 2 * t_exp_bub)),
                                     kc_meas),
    higbie_texp_times_renewal_rate=s_renewal * 2 * t_exp_bub,
    renewal_over_higbie_kL_at_measured=np.sqrt(np.pi * s_renewal
                                               * 2 * t_exp_bub) / 2)
brk("film thickness inverted from the PREDICTED k_L instead of the measured one",
    film_thickness_implied_um=DAB / kL_19 * 1e4)

# --- the reproduced printed tables --------------------------------------------
_phi_no_erf = lambda p: 1.0 / (1.0 + 1.0 / (np.sqrt(np.pi) * p * np.exp(p * p)))
_phi_b = np.array([phi_of_xA0(x, _phi_no_erf) for x in xs])
brk("the (1 + erf phi) factor dropped from Eq. 20.1-18",
    bsl_t2011_phi_worst_rel=float(np.abs(_phi_b / tt["phi_20p1_1"].to_numpy()
                                         - 1).max()),
    bsl_t2011_psi_worst_rel=float(np.abs(_phi_b * np.sqrt(np.pi) / xs
                                         / tt["psi_20p1_1"].to_numpy() - 1).max()))
brk("psi taken as phi/x_A0 (the sqrt(pi) of the Table 20.1-1 head dropped)",
    bsl_t2011_psi_worst_rel=float(np.abs(tt["phi_computed"].to_numpy() / xs
                                         / tt["psi_20p1_1"].to_numpy() - 1).max()),
    bsl_t2282_theta_pen_worst_rel=float(np.abs(
        (1 - xs) * tt["phi_computed"].to_numpy() / xs
        / tt["theta_penetration_22p8_2"].to_numpy() - 1).max()))
brk("Eq. 22.8-41 written as (1 - x) / psi instead of (1 - x) psi",
    bsl_t2282_theta_pen_worst_rel=float(np.abs(
        (1 - xs) / tt["psi_computed"].to_numpy()
        / tt["theta_penetration_22p8_2"].to_numpy() - 1).max()),
    bsl_t2282_theta_pen_minus_psi_dev=float(np.abs(
        ((1 - xs) / tt["psi_computed"].to_numpy()
         / tt["theta_penetration_22p8_2"].to_numpy() - 1)
        - tt["psi_rel"].to_numpy()).max()))
_theta_film_b = np.log(1 / (1 - xs))
brk("Eq. 22.8-42 written without the (1-x)/x prefactor",
    bsl_t2282_theta_film_worst_rel=float(np.abs(
        _theta_film_b / tt["theta_film_22p8_2"].to_numpy() - 1).max()),
    bsl_t2282_pen_film_gap_at_0p75=float(
        1 - tt.loc[tt["x_A0"] == 0.75, "theta_pen_computed"].iloc[0]
        / _theta_film_b[-1]))
brk("gamma/tanh(gamma) written upside down as tanh(gamma)/gamma",
    froment_t6421_film_worst_rel=float(np.abs(
        np.tanh(gam) / gam / t3["FA_film_printed"].to_numpy() - 1).max()),
    froment_penetration_vs_film_at_gamma_1=float(
        FA_penetration(1.0) / (np.tanh(1.0) / 1.0) - 1))
brk("sqrt(1 + gamma^2) written as 1 + gamma^2",
    froment_t6421_renewal_worst_rel=float(np.abs(
        (1 + gam ** 2) / t3["FA_surface_renewal_printed"].to_numpy() - 1).max()),
    froment_penetration_vs_renewal_at_gamma_1=float(
        FA_penetration(1.0) / 2.0 - 1))
_FA_pen_b = lambda x: x * ((1 + np.pi / (8 * x ** 2)) * erf(x)
                           + 1 / (2 * x) * np.exp(-4 / np.pi * x ** 2))
brk("erf(gamma) for erf(2 gamma/sqrt(pi)) in Froment Eq. 6.4.2-8",
    froment_t6421_penetration_worst_rel=float(np.abs(
        _FA_pen_b(gam) / t3["FA_penetration_printed"].to_numpy() - 1).max()),
    froment_t6421_pen_at_gamma_10_rel=float(
        _FA_pen_b(10.0) / t3["FA_penetration_printed"].to_numpy()[-1] - 1),
    froment_eq6428_closed_vs_quadrature_worst=float(np.abs(
        _FA_pen_b(gg) / qq - 1).max()))
brk("printed penetration column tested against a bracket ten times as wide",
    froment_t6421_pen_bracket_violations=float(bracket_violations(
        t3["FA_penetration_printed"].to_numpy(),
        t3["FA_film_printed"].to_numpy() - 0.09,
        t3["FA_surface_renewal_printed"].to_numpy() + 0.09).sum()))
# A TRANSCRIPTION defect, not a formula one: the two anomalous cells read as the
# values the closed form would give. This is the row the "worst over the table"
# metrics actually need, because one bad printed cell dominates them.
_pen_fixed = t3["FA_penetration_printed"].to_numpy().copy()
_pen_fixed[0], _pen_fixed[-1] = 0.99, 10.04
brk("the gamma = 0.01 and gamma = 10 penetration cells transcribed as 0.99 and 10.04",
    froment_t6421_penetration_worst_rel=float(np.abs(
        t3["penetration_computed"].to_numpy() / _pen_fixed - 1).max()),
    froment_t6421_pen_at_gamma_10_rel=float(
        t3["penetration_computed"].to_numpy()[-1] / _pen_fixed[-1] - 1),
    froment_t6421_pen_bracket_violations=float(bracket_violations(
        _pen_fixed, t3["FA_film_printed"].to_numpy(),
        t3["FA_surface_renewal_printed"].to_numpy()).sum()))

# --- the pymrm model ----------------------------------------------------------
_f = Film(LAM_V, XI=XIV, n_x=NXV, stretch=4.0, parabolic=False, wall="dirichlet")
_u, _P1, _ = _f.march(1600)
brk("outer boundary set to a sink (u = 0) instead of no flux",
    pen_uptake_vs_eq18518_raw_rel=rel(_P1, phi_penetration(LAM_V)),
    pen_profile_vs_eq18516_max_abs=float(
        np.abs(_u - erfc(_f.x_c / np.sqrt(4 * LAM_V))).max()),
    wall_series_vs_pymrm_worst_rel=abs(rel(
        Film(LAM_V, XI=1.0, n_x=NX, stretch=STR, parabolic=False,
             wall="dirichlet").march(2 * NT)[1], phi_wall_series(LAM_V))))
_f2 = Film(LAM_V, XI=np.sqrt(4 * LAM_V), n_x=NXV, stretch=4.0, parabolic=False)
_u2, _P2, _ = _f2.march(1600)
brk("domain truncated at one penetration depth instead of twelve",
    pen_uptake_vs_eq18518_raw_rel=rel(_P2, phi_penetration(LAM_V)),
    pen_profile_vs_eq18516_max_abs=float(
        np.abs(_u2 - erfc(_f2.x_c / np.sqrt(4 * LAM_V))).max()))
_f3 = Film(LAM_V, XI=XIV, n_x=NXV, stretch=4.0, parabolic=False)
_u3a, _P3a, _ = _f3.march(50)
_, _P3b, _ = Film(LAM_V, XI=XIV, n_x=NXV, stretch=4.0, parabolic=False).march(100)
brk("marching step coarsened to n_t = 50 (extrapolated over 50/100)",
    pen_uptake_vs_eq18518_raw_rel=rel(_P3a, phi_penetration(LAM_V)),
    pen_uptake_vs_eq18518_rel=rel(2 * _P3b - _P3a, phi_penetration(LAM_V)),
    pen_profile_vs_eq18516_max_abs=float(
        np.abs(_u3a - erfc(_f3.x_c / np.sqrt(4 * LAM_V))).max()),
    time_err_n1600=abs(rel(_P3a, phi_penetration(LAM_V))))
brk("convergence order read off the two COARSEST points instead of the two finest",
    time_order=order(tt_err[:2]), space_order=order(xx_err[:2]))
# The axis-isolation itself, broken: measure the grid error against the EXACT
# answer (so the marching error contaminates it) and the marching error on a
# grid too coarse to be innocent.
_xx_bad = [(n, abs(Film(LAM_V, XI=XIV, n_x=n, stretch=4.0,
                        parabolic=False).march(NT_FIX)[1] / ref - 1))
           for n in (100, 200, 400)]
_tt_bad = [(n, abs(Film(LAM_V, XI=XIV, n_x=25, stretch=4.0,
                        parabolic=False).march(n)[1] / ref - 1))
           for n in (800, 1600)]
brk("each axis measured without isolating it from the other",
    space_order=order(_xx_bad), time_order=order(_tt_bad))
_f25 = Film(LAM_V, XI=XIV, n_x=25, stretch=4.0, parabolic=False)
_u25, _P25, _ = _f25.march(NT_FIX)
brk("grid coarsened to n_x = 25",
    space_err_n400=abs(rel(_P25, a_fine)),
    pen_profile_vs_eq18516_max_abs=float(
        np.abs(_u25 - erfc(_f25.x_c / np.sqrt(4 * LAM_V))).max()))
_rp_flat = (phi_num(L_MARK, n_t=NT, XI=1.0, n_x=NX, stretch=STR, cap_override=1.0)
            / phi_num(L_MARK, n_t=NT, XI=1.0, n_x=NX, stretch=STR, parabolic=False))
_lam_star_flat = lam_at_1pct(r_wall, 0.05, 2.0)   # wall alone, wider bracket
brk("velocity profile flattened to plug flow inside the FULL model",
    profile_effect_at_lambda_0p1=1 - _rp_flat,
    profile_over_wall_effect_at_lambda_0p1=(
        (1 - _rp_flat) / M["wall_effect_series_at_lambda_0p1"]),
    lambda_star_1pct=_lam_star_flat,
    pen_depth_over_delta_at_1pct=np.sqrt(4 * _lam_star_flat),
    delta_star_18A4_um=np.sqrt(D_cl * L / (v_max * _lam_star_flat)) * 1e4)
brk("wall series truncated after 2 terms",
    wall_effect_series_at_lambda_0p1=1 - phi_wall_series(L_MARK, n_terms=2)
                                     / phi_penetration(L_MARK),
    wall_series_vs_pymrm_worst_rel=abs(rel(
        phi_num(1e-3, n_t=NT, XI=1.0, n_x=NX, stretch=STR, parabolic=False),
        phi_wall_series(1e-3, n_terms=2))))
brk("solver run at n_t = 100, so its own marching error swamps the wall effect",
    wall_effect_solver_at_lambda_0p1=1 - phi_num(L_MARK, n_t=100, XI=1.0, n_x=NX,
                                                 stretch=STR, parabolic=False)
                                     / phi_penetration(L_MARK))
brk("mixing cup taken without the (1 - xi^2) weight in the FULL model",
    mixing_cup_vs_flux_integral_rel=rel(float((f_full.V * u_full).sum()),
                                        flux_full))
brk("threshold read off the 12-point sweep by interpolation instead of root-found",
    lambda_star_1pct=lam_star_swept,
    pen_depth_over_delta_at_1pct=np.sqrt(4 * lam_star_swept),
    delta_star_18A4_um=np.sqrt(D_cl * L / (v_max * lam_star_swept)) * 1e4)
brk("threshold read off a 5-point sweep over the same range by interpolation",
    lambda_star_1pct=lam_star_swept_coarse,
    pen_depth_over_delta_at_1pct=np.sqrt(4 * lam_star_swept_coarse),
    delta_star_18A4_um=np.sqrt(D_cl * L / (v_max * lam_star_swept_coarse)) * 1e4)
brk("penetration solve for W_A truncated at one penetration depth",
    pymrm_vs_printed_18A4_rel=rel(
        W_wet * delta_ref * v_max * cA0_cl * 3600
        * phi_num(LAM_REF, n_t=NT_REF, XI=np.sqrt(4 * LAM_REF), n_x=NX_REF,
                  stretch=4.0, parabolic=False), WA_printed),
    pymrm_vs_eq18518_rel=rel(
        W_wet * delta_ref * v_max * cA0_cl * 3600
        * phi_num(LAM_REF, n_t=NT_REF, XI=np.sqrt(4 * LAM_REF), n_x=NX_REF,
                  stretch=4.0, parabolic=False), WA_closed_hr))
brk("arithmetic instead of harmonic face mean for the diffusivity",
    pen_uptake_vs_eq18518_raw_rel=M["pen_uptake_vs_eq18518_raw_rel"])

breaks = pd.DataFrame(BR, columns=["defect injected", "metric", "baseline",
                                   "with defect", "moved?"])
pd.set_option("display.max_colwidth", 74)
pd.set_option("display.max_rows", 100)
display(breaks)

# Coverage is GENERATED from the measured "moved?" column, not from the row
# labels: a metric whose only row leaves it bit-identical is NOT covered, and
# saying otherwise is how a headline number ends up with protection on paper only.
covered = set(breaks["metric"])
moving = set(breaks.loc[breaks["moved?"] == "yes", "metric"])
uncovered = [k for k in M if k not in covered]
no_mover = [k for k in M if k in covered and k not in moving]
print(f"\nbreak rows: {len(breaks)}   moved: {(breaks['moved?'] == 'yes').sum()}"
      f"   did not move: {(breaks['moved?'] == 'NO').sum()}")
print(f"metrics reported: {len(M)}   with at least one row: {len(covered & set(M))}"
      f"   with at least one MOVING row: {len(moving & set(M))}")
print("metrics with NO break row:", uncovered if uncovered else "none")
print("metrics with a row but NO MOVING row:", no_mover if no_mover else "none")
below_floor = [k for k, v in M.items() if abs(v) < 1e-12]
print("metrics below ABS_FLOOR = 1e-12 (outside CI):",
      below_floor if below_floor else "none")'''))

cells.append(md(r"""**What the rows that do not move are saying.** They are kept
deliberately, and they fall into five classes - each of which is a statement about
the page, not an oversight. Read every factor quoted below off the table itself;
nothing here is retyped.

Note which rows are **not** in this list. The water-density and $M$(Cl$_2$) rows
do move `eq18518_vs_printed_18A4_rel`, by about the size of the residual itself;
they are kept because that is what shows the printed three figures cannot decide
between the two unprinted inputs, not because they are inert. *Results 1* does
that arithmetic.

- **The metric is blind to the injected quantity by construction.** The
  $v_{\max}$, $\rho$ and $M$(Cl$_2$) rows each inject into *two* metrics:
  `eq18518_vs_printed_18A4_rel`, which they move, and `pymrm_vs_eq18518_rel`,
  which they cannot, because $\delta_{\rm ref}$ carries the whole dimensional
  prefactor back out of the marching route and it cancels out of the ratio. That
  pairing is the point of those rows: it is the demonstration behind the
  qualification in *Results 1*, and it is why the cross-route agreement is
  evidence about $\Phi(\Lambda)$ and not about $W$, $c_{A0}$ or the units. The row
  that drops the 4 from Eq. 18.5-18 moves the same metric, which is what fixes
  what it *does* see.
- **Structural.** The mixing cup and the accumulated surface flux are equal by
  construction under backward Euler on a conservative discretisation; and the
  diffusivity here is uniform, so the arithmetic and harmonic face means coincide -
  that row would matter only in a model with a jump in $\mathscr{D}$, and is kept
  as a marker.
- **The defect is genuinely harmless at the settings used, and other rows prove
  which setting is the assumption.** Changing the far boundary of the penetration
  solve from no-flux to a sink leaves the outlet profile **bit-identical**, because
  it sits twelve penetration depths out where $u$ is already zero. It would be
  wrong to read that as the profile check having nothing that can move it: three
  other rows inject into `pen_profile_vs_eq18516_max_abs` - *truncating the domain
  to one penetration depth*, *coarsening the marching step to $n_t=50$* and
  *coarsening the grid to $n_x=25$* - and all three move it, by the factors the
  table prints. The set says the domain **length**, the step and the mesh are the
  assumptions, and the far boundary-condition *type* is not. Similarly, reading the
  convergence orders off the two coarsest points instead of the two finest barely
  changes them - which is what "asymptotic" means - while measuring either axis
  without isolating it from the other destroys them.
- **A "worst over the table" metric can be pinned by one bad cell.** Breaking the
  penetration formula does not move
  `froment_t6421_penetration_worst_rel`, because that metric is dominated by the
  printed 0.94, and it does not move `froment_t6421_pen_at_gamma_10_rel` either,
  because $\operatorname{erf}$ is 1 to machine precision at both arguments when
  $\gamma=10$. The rows that do move them are the *transcription* row and, for the
  formula itself, `froment_eq6428_closed_vs_quadrature_worst`, which the same row
  lifts from round-off to a few percent. That is why the quadrature check exists.
- **A square root halves a defect.** `pen_depth_over_delta_at_1pct` and
  `delta_star_18A4_um` both scale as $\Lambda^{*\,\pm1/2}$, so the coarse-sweep row
  moves them by less than the 10 % the table calls a move even though it moves
  $\Lambda^{*}$ itself by more - compare the three rows in the table.
  $\Lambda^{*}$ is the metric that carries that row; the other two are reported
  because they are the readable form of it.

**A break table measures sensitivity, never correctness.** Nothing in it can catch
a wrong baseline. That is what the two independent routes are for: the
finite-volume march against Eq. 18.5-18 in *Results 1*, and the exact eigenfunction
series against the solver in *Validation 2*.

The cell above prints, rather than asserts, three things: which metrics have a
break row, which have a row that actually **moves** them - a row that leaves a
metric bit-identical is coverage on paper only, and the map is generated from the
measured `moved?` column rather than from the row labels - and which sit below
`ABS_FLOOR` = 1e-12 and are therefore outside the regression suite entirely.
"""))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**1. The equation BSL declines to solve.** Section 18.5 sets up Eq. 18.5-7 in full
and then says its analytic solution is "in the form of an infinite series" in
Pigford's 1941 thesis, which it does not print and which is not on disk. Solving it
directly turns BSL's two-sentence justification for Eq. 18.5-11 into a number: the
penetration formula is within 1 % of the full equation whenever the penetration
depth stays below $\sqrt{4\Lambda^{*}}$ of the film thickness, a threshold that is
root-found rather than read off a sweep.

**2. And it shows that the two halves of that justification are not equal.** BSL
offers two reasons for Eq. 18.5-11 - the solute does not sense the wall, and it
thinks the film moves at $v_{\max}$ - as though they were the same argument. They
are not: at the point where the front has crossed two thirds of the film, the
velocity profile is worth a couple of percent and the wall is so small that the
finite-volume solver cannot resolve it and the exact series has to be brought in to
measure it. The reason is visible once the wall case is solved exactly: a
reflecting wall returns solute that a semi-infinite film would have retained
anyway, whereas slower liquid below the surface has genuinely longer to load.
**That decomposition is not in either book.**

**3. The three-way comparison, taken rather than declined - with its limits
stated.** `A3.1` and `A3.3` both leave film-versus-penetration-versus-renewal open
because neither source can separate the pictures. Two things on disk let this page
close part of it: Froment's Table 6.4.2.1, which puts all three enhancement factors
side by side and whose penetration column turns out to disagree with its own
printed formula; and Hammerton & Garner's single measured $k_c$, on which the
penetration model is the only one of the three that makes a prediction without a
fitted constant. **The exponent question is not closed and this page says so** -
one $(k_L,\mathscr{D})$ pair cannot resolve $n$, exactly as `A3.3` established from
the other side.

**4. Two printed defects, proved from the sources' own numbers.** Problem 18A.7(a)
names Eq. 18.5-20 and prints Eq. 18.5-19's answer, a factor $\sqrt3$ apart; and two
cells of Froment's Table 6.4.2.1 are impossible against the printed digits of the
other two columns of their own rows. Both are reported and neither is repaired.

**What pymrm does not add.** The closed forms of section 18.5 are exact and need no
solver; the reproductions in *Results 1*, *5* and *6* are arithmetic, and the
notebook is honest that they are. The Arnold functions of Table 20.1-1 are `A4.8`'s
subject and are used here only as inputs to Eq. 22.8-41.
"""))

# ---------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**If you want a liquid-side coefficient for a short-contact device**, this is the
page: $k_L=2\sqrt{\mathscr{D}/\pi t_{\rm exp}}$ with $t_{\rm exp}$ from the
geometry - $L/v_{\max}$ for a wetted wall, $D/v_t$ for a circulating bubble. On the
one measurement BSL quotes, that runs about 12 % low with no adjustable constant.
Use Eq. 18.5-19 for a mobile interface and Eq. 18.5-20 for creeping flow, and note
that a surfactant skin moves you off this page entirely: BSL says the exponent then
becomes $\mathscr{D}^{1/3}$, as for a rigid sphere.

**If you are choosing between film, penetration and surface renewal**, read
*Results 3* and *4* before *Results 5*. On the practical axis - an enhancement
factor for a pseudo-first-order reaction - the three differ by a couple of percent
over the whole useful range, which is Froment's own conclusion and is reproduced
here; the choice is not worth agonising over for design. On the *scientific* axis -
the exponent on diffusivity - nothing on disk settles it, and the sibling pages
[`A3.1`](../A3.1-whitman-two-film/) and
[`A3.3`](../A3.3-danckwerts-surface-renewal/) explain from their own sources why.
What does distinguish the three is how many free constants they need: two of them
need one each, and penetration needs none.

**If you are building a falling-film or wetted-wall model**, the `Film` class here
is the skeleton: a capacity matrix times a marching derivative, `construct_grad`
and `construct_div` for the transverse operator, boundary conditions on the outward
normal. Substituting a reaction term or a second species is a change to the
residual, not to the assembly. `A3.3` marches the same structure in surface age
rather than axial distance, and `F3.1` (Hatta) is where reaction enhancement in
this geometry belongs.

**Do not take from this page**: a validated diffusivity exponent (none is
established); an error bar on the 12 % (BSL prints none for Hammerton & Garner's
number and none is invented); a film thickness for Problem 18A.4's column (not
printed, not inferred); or Froment's Table 6.4.2.1 penetration column as data (two
of its six cells are shown here to be inconsistent with the formula printed above
them).

**Cite the sources, not this page.** The origin is R. Higbie, *Trans. AIChE* **31**,
365-389 (1935), which **was not consulted**. The model as used here is read from
Bird, Stewart & Lightfoot, *Transport Phenomena*, 2nd edn (Wiley, 2002), section
18.5, with the worked problems 18A.4 and 18A.7 and Tables 20.1-1 and 22.8-2; the
cross-check is Froment, De Wilde & Bischoff, *Chemical Reactor Analysis and Design*,
3rd edn (Wiley, 2011), section 6.4 and Table 6.4.2.1, the last of which Froment
himself takes from Beek (1968), also not consulted. Hammerton & Garner,
*Trans. Inst. Chem. Engrs. (London)* **32**, S18-S24 (1954) supplies the one
measurement and is quoted only as BSL quotes it.
"""))

cells.append(code(r'''report_agreement("A3.2", M)'''))

# --------------------------------------------------------------------- write
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
