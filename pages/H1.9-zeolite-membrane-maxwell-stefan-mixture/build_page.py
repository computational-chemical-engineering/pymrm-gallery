#!/usr/bin/env python3
"""Generate index.ipynb for page H1.9. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# ---------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Maxwell–Stefan mixture diffusion in a zeolite membrane: a selectivity of 487 where the pure components predict 1.73"
description: "Two things break the pure-component prediction — the mixture isotherm is not multicomponent Langmuir, and the exchange coefficients are finite. Both are reproduced here, and against Bakker's and Funke's measurements the model is within 4 % on one case and 3.4x out on another the review passes over."
categories: [sec:H, struct:S9, struct:S4, tier:T1, data:tier2, phase:gas-solid]
date: 2026-08-03
---

# Maxwell–Stefan mixture diffusion in an MFI membrane

**Catalog ID:** `H1.9` · **Structures:** `S9` (multicomponent transport), `S4`
(1-D transient PDE) · **Tier:** T1

A methane–butane mixture is fed to a silicalite membrane at 95 and 5 kPa.
Methane diffuses a hundred times faster than butane and is nineteen times more
abundant, so the pure-component fluxes predict that butane will come through the
membrane 1.73 times more selectively than methane. It comes through **487**
times more selectively.

Krishna & Baur's explanation has two parts, and this page implements both: the
mixture isotherm is not the multicomponent Langmuir (size and configurational
entropy push the larger molecule out near saturation, so it needs IAST), and the
exchange coefficients $\eth_{ij}$ between the two sorbates are **finite**, so
the fast species is dragged along by the slow one. Take away the finite exchange
coefficients and the selectivity falls by a factor of forty; take away the
entropy effects IAST carries and it stops varying with feed composition at all —
where the model that keeps them falls steeply as the butane fraction rises, and
so do Bakker's two measurements.

The mixture side is where the review's experimental comparisons are, and they
are stated as numbers rather than drawn in a figure. Those are the only
measurements this page can be scored against, and they are used."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

Inside a zeolite pore every molecule is adsorbed all of the time, so "diffusion"
is a hopping problem on a lattice of sites and a "concentration" is an occupancy.
Page [`A4.7`](../A4.7-zeolite-micropore-maxwell-stefan/) works the
single-component case from the same review: the measured Fick coefficient is a
mobility $\eth$ times a thermodynamic factor $\Gamma$ read off the isotherm, and
the awkward loading dependence of $D$ lives in $\Gamma$.

**With two sorbates, three things appear that have no single-component
counterpart.**

1. $\Gamma$ becomes a **matrix**, and it is non-diagonal even for an ideal
   mixture: changing the loading of species 2 changes the fugacity of species 1,
   because they compete for the same sites.
2. There are **exchange coefficients** $\eth_{ij}$. A molecule that hops leaves a
   vacancy; whichever species fills it is correlated with the one that left. That
   friction between the sorbates slows a fast species and speeds up a slow one,
   and it is *not* present in the pure-component data.
3. The **mixture isotherm is not the pure isotherms side by side**. Near
   saturation the smaller molecule wins the remaining space on entropy grounds,
   whatever the adsorption energies say, and the multicomponent Langmuir cannot
   represent that at all.

This page takes Sections 3.1–3.4 and 4.2–4.3 of the review — the mixture theory
and the MFI membrane separations it was built for. The single-component flux
law, the isotherm machinery and the Kärger worked examples are `A4.7`'s and are
not repeated. The cut between the two cases is single-component / mixture rather
than micropore / membrane; the reason is under *Reuse* on both pages."""))

# ---------------------------------------------------------------- the model
cells.append(md(r"""## The published model

**The flux relations.** For $n$ sorbates in a stationary framework the review
writes the chemical potential gradients as linear functions of the fluxes
(eq. 42),

$$-\frac{\rho\,\theta_i}{RT}\nabla\mu_i
= \sum_{j\neq i}\frac{\Theta_jN_i-\Theta_iN_j}
{\Theta_{i,\mathrm{sat}}\Theta_{j,\mathrm{sat}}\eth_{ij}}
+ \frac{N_i}{\Theta_{i,\mathrm{sat}}\eth_i},
\qquad i=1\ldots n .$$

At $n=1$ the sum is empty and this is `A4.7`'s eq. (9). Defining
$u_i \equiv N_i/(\rho\,\Theta_{i,\mathrm{sat}})$ and using
$\nabla\mu_i/RT=\nabla\ln f_i$, eq. (42) becomes a small linear system

$$-\theta_i\,\nabla\ln f_i
= \sum_{j\neq i}\frac{\theta_ju_i-\theta_iu_j}{\eth_{ij}}
+ \frac{u_i}{\eth_i}
\;\equiv\;\sum_j B_{ij}u_j ,$$

whose matrix is the review's eq. (44),
$B_{ii}=1/\eth_i+\sum_{j\neq i}\theta_j/\eth_{ij}$,
$B_{ij}=-\theta_i/\eth_{ij}$. With the matrix thermodynamic factor of eq. (38),

$$\Gamma_{ij}\equiv\frac{\theta_i}{\theta_j}
\frac{\partial\ln f_i}{\partial\ln\theta_j},$$

the driving force $-\theta_i\nabla\ln f_i$ is $-\sum_j\Gamma_{ij}\nabla\theta_j$
and the flux relations close as eq. (45),

$$(N)=-\rho\,[\Theta_{\mathrm{sat}}]\,[B]^{-1}[\Gamma]\,(\nabla\theta),
\qquad [D]=[B]^{-1}[\Gamma] \quad\text{(eq. 46)} .$$

**The exchange coefficients** are not measurable from pure components. The
review estimates them with the Vignes-type interpolation of eq. (43),

$$\eth_{ij}=\eth_i^{\,\theta_i/(\theta_i+\theta_j)}\;
            \eth_j^{\,\theta_j/(\theta_i+\theta_j)} .$$

Letting $\eth_{ij}\to\infty$ removes the sorbate–sorbate friction entirely and
eq. (45) collapses to eq. (47), which is **Habgood's model** — the review is
explicit that this is a limiting case of the M–S approach, not a rival to it.

**The mixture isotherm** is ideal adsorbed solution theory (IAST), eqs. (62)–(65).
With $\varphi$ the surface potential and $P_i^0(\varphi)$ the pressure at which
pure $i$ has that potential,

$$Py_i = P_i^0(\pi)\,x_i,\qquad
\frac{\varphi}{k_BT}=\rho\!\int_0^{P_i^0}\!\Theta_i^0(P)\,\frac{\mathrm dP}{P},
\qquad
\Theta_{\mathrm{mix}}=\left[\sum_i\frac{x_i}{\Theta_i^0(P_i^0)}\right]^{-1}.$$

The alternative, eq. (39)'s multicomponent Langmuir
$\Gamma_{ij}=\delta_{ij}+\theta_i/(1-\sum_k\theta_k)$, is thermodynamically
consistent **only when all saturation capacities are equal**, which the review
states and which is exactly what fails here: C1 and nC4 hold 19 and 10 molecules
per unit cell.

**The membrane**, eqs. (57)–(59): loadings fixed at both faces by equilibrium
with the upstream and (swept) downstream partial pressures, and

$$\frac{\partial\theta_i}{\partial t}
=-\frac{1}{\rho\,\Theta_{i,\mathrm{sat}}}\frac{\partial N_i}{\partial z}.$$

The two selectivities, eqs. (60) and (61):

$$S_P=\frac{N_2/N_1}{p_{2,0}/p_{1,0}},
\qquad S=\frac{\Theta_2/\Theta_1}{p_2/p_1}.$$

*How those two were read, since it matters.* The only form of this paper on disk
is the Elsevier full text — there is no page image, so nothing here is read off a
render. That text flattens both equations in the same way and with the same
spacing, `(60) Sp = N2  N1   p20  p10` and `(61) S = Θ2  Θ1   p2  p1`: two spaces
inside each inner fraction, three between the two groups. That is the ordinary
flattening of a **compound fraction**, and both are read that way above —
numerator ratio *divided by* denominator ratio. The review's own worked value for
eq. (60), $3.1/34.0/(5/95)=1.73$, is the compound-fraction reading, so the
printed equation and the printed number agree and there is nothing to repair. The
rival reading — multiplying by the pressure ratio — is run in the break table
below and lands a factor of 361 from every printed selectivity, which is what
settles the reading. **That is a statement about which reading is meant, not a
claim that the paper prints eq. (60) wrongly**; the flattened text cannot support
the second claim, and an earlier draft of this page made it.

**Self-diffusivity**, eqs. (53)–(56): labelling species 1 and applying the
constraints $\nabla\theta_{2,3}=0$, $N_1=-N_{1^*}$ to eq. (45) gives

$$D_1^*=\left[\frac{1}{\eth_1}+\frac{\theta_1}{\eth_{11}}
        +\frac{\theta_2}{\eth_{12}}+\cdots\right]^{-1},$$

which contains no $[\Gamma]$ at all: sorption thermodynamics play no part in a
self-diffusivity.

**Two printed relations here are algebraic identities and are labelled as such
rather than presented as evidence.** Applying eq. (43) at $i=j$ returns
$\eth_{ii}=\eth_i$ — that is what "$x^a y^{1-a}$ with $y=x$" does, and it is not
a test of anything. And eq. (55), $D_1^*=\eth_1/(1+\theta_1)$, follows from
eq. (54) by substituting $\eth_{11}=\eth_1$ and dropping the other terms. Both
are stated because the review uses them; neither is scored."""))

# ---------------------------------------------------------------- parameters
cells.append(md(r"""## Parameters and assumptions

**Carried from the review.** Rigid framework; local equilibrium between the
sorbed phase and an ideal gas, so $f_i=p_i$; isothermal operation; no support
resistance and no membrane defects (Section 4.1 states both); downstream partial
pressures held at "vanishing values" by a sweep gas.

**Left as parameters, because the review does not fix them.**

- *The loading dependence of $\eth_i$ in the mixture.* Section 3.2 says the
  $\eth_i$ "depend on the loading following e.g. Eq. (17) or Eq. (19)" and never
  says which applies to the membrane examples. The production runs take
  eq. (17), $\eth_i=\eth_i(0)$, because that is what reproduces the
  single-component fluxes of Section 4.1 (`A4.7` establishes it, to 0.17 %).
  The sensitivity to eq. (19) instead is measured in *Validation*.
- *The downstream pressure.* "Vanishing" is implemented as a fixed small
  fraction of the upstream partial pressures. The sensitivity to that fraction
  is measured, and it is small.
- *The membrane thickness $\delta$.* Never printed for the C1/nC4 case. It
  cancels out of every selectivity, and — less obviously, so it is demonstrated
  rather than asserted below — out of the four *absolute* mixture fluxes as
  well, once it is reconstructed from a printed flux the way `A4.7` reconstructs
  it. **Exactly one quantity on this page carries it: the transient peak time,
  through $\delta^2$.** That is where `A4.7`'s $\approx 40\ \mu$m is tested, it
  is the only place it is tested, and it does not agree.

**What sets the scale, and what does not.** $S_P$ depends only on ratios: the
membrane thickness, the framework density and the overall magnitude of the
diffusivities all cancel. That is why the hexane results below can be computed
at all — the review prints $\eth_{n\mathrm{C6}}=5\,\eth_{3\mathrm{MP}}$ and
never prints either one."""))

# ---------------------------------------------------------------- env
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
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve, gmres, LinearOperator
from scipy.optimize import brentq
from scipy.integrate import quad
from pymrm import construct_grad, construct_div, NumJac
from gallery_utils import load_data, report_agreement

PAGE = "H1.9-zeolite-membrane-maxwell-stefan-mixture"
A47 = "A4.7-zeolite-micropore-maxwell-stefan"
np.random.seed(20260803)      # only the friction-system check draws numbers
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

# ---------------------------------------------------------------- the data
cells.append(md(r"""## The data

Three CSVs of this page's own, and two borrowed from `A4.7`.

**Own — `bakker-funke-mfi-membrane-selectivities.csv`, and it is why this page
is tier 2 rather than tier 6.** Four permeation selectivities *measured* on MFI
membranes: Bakker's 380 and 60 for methane–butane, and Funke et al.'s 24 for a
hexane-isomer mixture together with 1.3 for the same pair as pure components.
All four are stated in the review's **running text**, not drawn in a figure — so
on [`docs/data-strategy.md`](../../docs/data-strategy.md)'s ladder they are
tier 2, method `table`: numbers already in numeric form, with no digitisation
error at all. Neither Bakker's thesis nor Funke et al. is on disk; the review is
the document actually read, and both origins are recorded under
`origin_not_consulted` in the sidecar.

**Own — `krishna-baur-2003-mixture-printed-results.csv`**: every mixture result
the review prints as a number, which are the reproduction targets. **Own —
`krishna-baur-2003-mixture-parameters.csv`**: the parameters that appear only in
the mixture sections (the Sanborn & Snurr and Snurr & Kärger simulation
conditions, and Section 4.3's diffusivity ratios and time scale).

**Borrowed — `A4.7`'s two files**, for the DSL isotherm parameters of Tables 1
and 2 and for the Section 4.1 single-component results. Loading another page's
dataset means reading that page, so here is every finding `A4.7` states about
the rows used here and whether it bears on this one.

| `A4.7` finding | Rows | Bears on this page? |
|---|---|---|
| The two printed 3MP loadings are reproduced from Table 2's 3MP row to 0.19 % | Table 2, 3MP | **Yes** — that row drives every hexane result below, and it is certified |
| The printed single-component flux **ratio** is reproduced from the Table 1 C1 and nC4 rows to 0.17 % | Table 1, C1 and nC4 | **Yes, centrally.** Those two rows are the whole C1/nC4 mixture calculation, and the certification is used below as evidence in its own right |
| The review's printed $\eth_z(0)$ is exactly **half** what eq. (23) gives, confirmed out of sample on CH₄; and its two printed orientation-averaged diffusivities disagree by 4.1 % | Table 3, Section 2.2 | **No.** Kärger's relations produce zero-loading diffusivities for the *KMC* study of Section 3.3; no number on this page passes through eq. (23). `A4.7` declines to say whether the factor of two is Kärger's or the review's, and nothing here settles it either |
| $\Gamma$'s "two extrema" hold for three species the review does not name — but that survey has no printed target, and `A4.7`'s `GAMMA_IDENTITY` check is blind to a wrong Table 1 parameter by construction | Table 1, ten rows | **No.** This page quotes neither the sentence nor the survey, and the only Table 1 rows it uses are C1 and nC4 — the two `A4.7` certifies against a printed target |
| Three printed defects: the sign of eq. (25); $\nu_{zz}+\nu_{zz}$ in the denominator of eq. (23c); and $\nu_{\mathrm{str}}$ printed twice in the CH₄ sentence (a flagged repair by inference) | Sections 2.2, 2.4 | **No.** All three are single-component. Two further defects found on the mixture side are recorded below |
| The review's practice uses a halved $\eth_z$ in both worked applications of eq. (23) | Section 2.2, Table 3 | **No.** Same reason as above; stated so that it is not left unaddressed |

**And the rule that catches most of this class: no number that is a row in a
borrowed CSV is *stated as a result* without this page's own value being printed
beside it and the two reconciled.** The single-component fluxes 34 and 3.1
mmol m⁻² s⁻¹, the pure-component selectivity 1.73, the upstream pressures and
the two M–S diffusivities are all rows of `A4.7`'s results file. Every one of
them is loaded, printed beside this page's own value and reconciled — the
review's whole rhetorical setup ("we might expect 1.73; it is 487") runs
straight through them. Note the precise form of the claim: the *narrative*
markdown does name 1.73, 34, 3.1, 95 and 5 kPa, because a markdown cell cannot
interpolate a computed value. What it never does is assert one of them as an
outcome without the code printing this page's counterpart next to it.

**Deliberately scoped out.** Sections 4.4–4.7 — the N₂/CH₄ uptake in 4A
(Table 4), the O₂/N₂ chromatographic separation (Table 5), the C₅/C₆
breakthroughs (Table 6) and the pulsed hexane separation. Every one is a mixture
problem whose only validation is a figure; the review prints no numerical result
for any of them. Digitising a figure needs a maintainer review and none is
available, so none of them is transcribed and no claim here rests on one. Their
parameter tables are clean and would transcribe easily, so they are viable
future cases — a packed-bed chromatography page in section J is the natural
home. The same applies to Figs. 17–21, the review's MD/KMC verification of the
mixture theory: the theory is implemented and exercised below, but there is
nothing printed to score it against and the page says so where it happens."""))

cells.append(code('''par = load_data("krishna-baur-2003-sorption-parameters.csv", page=A47)
a47 = load_data("krishna-baur-2003-printed-results.csv", page=A47)
res = load_data("krishna-baur-2003-mixture-printed-results.csv", page=PAGE)
mpar = load_data("krishna-baur-2003-mixture-parameters.csv", page=PAGE)
expt = load_data("bakker-funke-mfi-membrane-selectivities.csv", page=PAGE)

PRINTED = dict(zip(res.result_id, res.value))       # this page's targets
A47_VAL = dict(zip(a47.result_id, a47.value))       # A4.7's rows - never retyped
MEAS = dict(zip(expt.measurement_id, expt.S_P_measured))

print(f"{len(res)} printed mixture targets, {len(mpar)} own parameter rows, "
      f"{len(expt)} MEASURED selectivities\\n")
print(res[["result_id", "where", "printed_as", "units"]].to_string(index=False))
print("\\nMeasured (tier 2, read through the review):")
print(expt[["measurement_id", "origin", "mixture", "composition",
            "upstream_pressure_Pa", "printed_as"]].to_string(index=False))


def params(component, source_table):
    s = par[(par.component == component) & (par.source_table == source_table)]
    if s.empty:
        raise KeyError(f"{component!r} not in {source_table!r}")
    return dict(zip(s.symbol, s.value))


def mprm(component, section, symbol):
    s = mpar[(mpar.component == component) & (mpar.source_section == section)
             & (mpar.symbol == symbol)]
    return float(s.value.iloc[0])'''))

# ---------------------------------------------------------------- implementation
cells.append(md(r"""## PyMRM implementation

Four objects. `DSL` is the dual-site Langmuir isotherm and its surface potential
— `A4.7`'s `Isotherm` with the potential promoted from a helper to the main
interface, because IAST is built on it. `IAST` solves eqs. (62)–(65) in **both**
directions: forward $(p)\to(\Theta)$ for the boundary conditions, and inverse
$(\Theta)\to(p)$ for the driving force, which is the direction the membrane
needs and the one the review does not write down. `B_matrix` is eq. (44), and
`Membrane` is the pymrm transient solve of eqs. (57)–(59).

Three implementation decisions worth stating, because each is a place a page can
quietly go wrong.

**The IAST inverse is a *scalar* root find, not an $n$-dimensional one.** Given
the loadings, the adsorbed-phase mole fractions $x_i=\Theta_i/\Theta_{\rm mix}$
and the total $\Theta_{\rm mix}$ are already known, so the only unknown left in
eqs. (62)–(65) is the surface potential — and $\Theta_{\rm mix}(\varphi)$ is
monotone. Log-bisection to bracket, then Newton to finish.

**$[\Gamma]$ comes from implicit differentiation of the IAST equations, not from
differencing the solve.** The chain rule through
$\ln f_i=\ln x_i+\ln P_i^0(\varphi)$ gives every element in closed form once
$\partial\varphi/\partial\Theta_j$ is had from the constraint. That is fast
enough to sit inside a Newton iteration — and it is a *different* route from
central differences, so the two can be compared, and they are.

**$[B]^{-1}$ is formed once from eq. (44) and once not at all.** Coding [B] twice
cannot detect an error inside [B] — that is `A4.2`'s finding, and its fix is
reused here: the friction system of eq. (42) is also solved matrix-free by GMRES,
applying the summation as written and never assembling anything. The two agree
to machine precision and disagree by hundreds of per cent when a term is dropped
from $B_{ii}$; both numbers are printed in *Validation*.

Boundary conditions use the outward normal, and both faces are Dirichlet
(eqs. 57–58), so the sign trap `AGENTS.md` warns about does not arise here — the
same `{a:0, b:1, d:Θ}` dict means the same thing at both ends. `nu=0` in
`construct_div`: the membrane is a slab."""))

cells.append(code('''class DSL:
    """Dual-site Langmuir, the review's eq. (32), plus the surface potential of
    eq. (64). Setting Theta_sat_B = 0 gives the single-site Langmuir."""

    def __init__(s, b_A, Th_sat_A, b_B=0.0, Th_sat_B=0.0):
        s.bA, s.tA, s.bB, s.tB = float(b_A), float(Th_sat_A), float(b_B), float(Th_sat_B)
        s.Th_sat = s.tA + s.tB
        s.H = s.tA * s.bA + s.tB * s.bB          # Henry coefficient, dTheta/dP at P -> 0
        s._lp = np.linspace(-45.0, 45.0, 9001)   # seed table for the inverse
        s._tab = s.psi(10.0 ** s._lp)

    def theta(s, P):
        r = s.tA * s.bA * P / (1 + s.bA * P)
        return r + (s.tB * s.bB * P / (1 + s.bB * P) if s.tB > 0 else 0.0)

    def dtheta(s, P):
        r = s.tA * s.bA / (1 + s.bA * P) ** 2
        return r + (s.tB * s.bB / (1 + s.bB * P) ** 2 if s.tB > 0 else 0.0)

    def psi(s, P):
        """int_0^P Theta^0(P') dP'/P' -- eq. (64) divided by rho k_B T."""
        r = s.tA * np.log1p(s.bA * P)
        return r + (s.tB * np.log1p(s.bB * P) if s.tB > 0 else 0.0)

    def P_of_psi(s, ps):
        """Invert the surface potential. Monotone, so a table seed + Newton."""
        ps = np.asarray(ps, float)
        P = 10.0 ** np.interp(ps, s._tab, s._lp)
        for _ in range(6):                       # dpsi/dlnP = Theta^0
            st = np.clip(-(s.psi(P) - ps) / np.maximum(s.theta(P), 1e-300), -50.0, 50.0)
            P = np.clip(P * np.exp(st), 1e-300, 1e300)
        return P


def _solve_psi(g, dg_dlnpsi, m, incr, seed=None):
    """Monotone scalar solve for the IAST surface potential, per grid point.

    Log-bisection to bracket (robust, cannot diverge) then Newton in ln(psi)
    (quadratic). Both stages are vectorised over the m grid points; a seed from
    the previous call narrows the bracket to six decades.
    """
    lo = np.full(m, 1e-14) if seed is None else np.maximum(seed * 1e-3, 1e-300)
    hi = np.full(m, 1e2) if seed is None else seed * 1e3
    for _ in range(40):
        bad = (g(lo) > 0) if incr else (g(lo) < 0)
        if not np.any(bad):
            break
        lo = np.where(bad, lo / 100, lo)
    for _ in range(40):
        bad = (g(hi) < 0) if incr else (g(hi) > 0)
        if not np.any(bad):
            break
        hi = np.where(bad, hi * 100, hi)
    for _ in range(12 if seed is not None else 30):
        mid = np.sqrt(lo * hi); gm = g(mid)
        up = (gm < 0) if incr else (gm > 0)
        lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    psi = np.sqrt(lo * hi)
    for _ in range(5):
        psi = np.clip(psi * np.exp(np.clip(-g(psi) / dg_dlnpsi(psi), -2.0, 2.0)), lo, hi)
    return psi


class IAST:
    """Ideal adsorbed solution theory, eqs. (62)-(65), both directions."""

    def __init__(s, isos):
        s.isos = isos; s.n = len(isos)

    def _P0(s, psi):
        return np.stack([s.isos[i].P_of_psi(psi) for i in range(s.n)], -1)

    def _T0(s, P0):
        return np.stack([s.isos[i].theta(P0[..., i]) for i in range(s.n)], -1)

    def forward(s, p, seed=None):
        """(p) -> (Theta). eq. (62) closes as sum_i p_i/P_i^0(psi) = 1."""
        p = np.atleast_2d(np.asarray(p, float))

        def g(psi):
            return np.sum(p / s._P0(psi), -1) - 1.0        # decreasing in psi

        def dg(psi):
            P0 = s._P0(psi)
            return -psi * np.sum(p / (P0 * s._T0(P0)), -1)

        psi = _solve_psi(g, dg, p.shape[0], incr=False, seed=seed)
        s.psi_last = psi
        P0 = s._P0(psi); x = p / P0
        return x * (1.0 / np.sum(x / s._T0(P0), -1))[..., None]     # eq. (65)

    def _psi_of(s, Th, seed=None):
        """(Theta) -> psi. x and Theta_mix are already known from (Theta), so
        eq. (65) is a SCALAR monotone equation for the surface potential."""
        Th = np.atleast_2d(np.asarray(Th, float))
        Tt = np.sum(Th, -1); x = Th / Tt[..., None]

        def g(psi):
            return 1.0 / np.sum(x / s._T0(s._P0(psi)), -1) - Tt    # increasing

        def dg(psi):
            P0 = s._P0(psi); T0 = s._T0(P0)
            dT0 = np.stack([s.isos[i].dtheta(P0[..., i]) * P0[..., i]
                            for i in range(s.n)], -1) / T0
            S = np.sum(x / T0, -1)
            return psi * np.sum(x * dT0 / T0 ** 2, -1) / S ** 2

        psi = _solve_psi(g, dg, Th.shape[0], incr=True, seed=seed)
        s.psi_last = psi
        return psi, x

    def inverse(s, Th, seed=None):
        psi, x = s._psi_of(Th, seed)
        return x * s._P0(psi)

    def inverse_and_gamma(s, Th, Th_sat, seed=None):
        """(Theta) -> (f) and [Gamma] of eq. (38), from ONE IAST solve.

        [Gamma] is obtained by implicit differentiation of eqs. (62)-(65):
        ln f_i = ln x_i + ln P_i^0(psi), so
            dln f_i/dTheta_j = delta_ij/Theta_i - 1/Theta_mix
                               + (1/Theta_i^0) dpsi/dTheta_j,
        and dpsi/dTheta_j comes from differentiating eq. (65). Compared against
        central differences in Validation - a slip in this algebra shows there.
        """
        Th = np.atleast_2d(np.asarray(Th, float)); n = s.n
        psi, x = s._psi_of(Th, seed)
        P0 = s._P0(psi); T0 = s._T0(P0)
        dT0 = np.stack([s.isos[i].dtheta(P0[..., i]) for i in range(n)], -1) * P0 / T0
        Tt = np.sum(Th, -1); S = np.sum(x / T0, -1)
        dF_dpsi = np.sum(x * dT0 / T0 ** 2, -1) / S ** 2
        inv = 1.0 / T0
        dF_dT = -(inv / Tt[:, None] - np.sum(x * inv, -1)[:, None] / Tt[:, None]) \\
            / S[:, None] ** 2 - 1.0
        dpsi_dT = -dF_dT / dF_dpsi[:, None]
        dlnf = (np.eye(n)[None] / Th[:, :, None] - 1.0 / Tt[:, None, None]
                + dpsi_dT[:, None, :] / T0[:, :, None])
        # Gamma_ij = (th_i/th_j) dln f_i/dln th_j = th_i (dln f_i/dTheta_j) Theta_j,sat
        return x * P0, (Th / Th_sat)[:, :, None] * dlnf * Th_sat[None, None, :]


def B_matrix(th, D_i, mode="vignes"):
    """eq. (44), with the exchange coefficients from eq. (43)."""
    m, n = th.shape
    B = np.zeros((m, n, n))
    for i in range(n):
        B[:, i, i] = 1.0 / D_i[:, i]
        if mode == "inf":                # eq. (47): the Habgood limit
            continue
        for j in range(n):
            if j == i:
                continue
            xi = th[:, i] / np.maximum(th[:, i] + th[:, j], 1e-300)
            Dij = D_i[:, i] ** xi * D_i[:, j] ** (1 - xi)          # eq. (43)
            B[:, i, i] += th[:, j] / Dij
            B[:, i, j] = -th[:, i] / Dij
    return B


def gamma_numerical(iast, Th, Th_sat, h=1e-5):
    """[Gamma] by central differences of the IAST inverse map, in ln(Theta).
    The independent route against which inverse_and_gamma is checked."""
    Th = np.atleast_2d(np.asarray(Th, float)); m, n = Th.shape
    th = Th / Th_sat; G = np.zeros((m, n, n))
    for j in range(n):
        Tp = Th.copy(); Tp[:, j] *= (1 + h)
        Tm = Th.copy(); Tm[:, j] *= (1 - h)
        d = (np.log(iast.inverse(Tp)) - np.log(iast.inverse(Tm))) \\
            / (np.log1p(h) - np.log1p(-h))
        for i in range(n):
            G[:, i, j] = th[:, i] / th[:, j] * d[:, i]
    return G


def friction_solve(th, D_i, d, mode="vignes"):
    """Solve eq. (42) for the velocities WITHOUT assembling [B] anywhere.

    The matrix-vector product IS the right-hand side of eq. (42) as printed, so
    this route cannot inherit an error in the B_ii / B_ij formulas of eq. (44).
    A4.2's finding, reused: coding [B] twice tests arithmetic, not physics.
    """
    n = len(th)

    def Dij(i, j):
        xi = th[i] / max(th[i] + th[j], 1e-300)
        return D_i[i] ** xi * D_i[j] ** (1 - xi)

    def mv(u):
        r = np.empty(n)
        for i in range(n):
            acc = u[i] / D_i[i]
            if mode != "inf":
                for j in range(n):
                    if j != i:
                        acc += (th[j] * u[i] - th[i] * u[j]) / Dij(i, j)
            r[i] = acc
        return r

    u, _ = gmres(LinearOperator((n, n), matvec=mv), d, rtol=1e-14, atol=0.0)
    return u'''))

cells.append(code('''class Membrane:
    """eqs. (57)-(59) with the flux relations eq. (45), on a pymrm slab grid.

    State is Theta (n_z, n_c) -- spatial axis first, field axis last. Fluxes are
    returned as N/rho in molecules m / (s . unit cell), so that a selectivity
    needs neither the framework density nor the membrane thickness.
    """

    def __init__(s, isos, Th_sat, D0, n_z=30, Dij="vignes", loading_dep=None,
                 iso_model="iast", delta=1.0, pd_frac=1e-6):
        s.iast = IAST(isos); s.isos = isos; s.n = len(isos)
        s.Th_sat = np.asarray(Th_sat, float); s.D0 = np.asarray(D0, float)
        s.n_z, s.Dij, s.iso_model, s.delta = n_z, Dij, iso_model, delta
        s.pd_frac = pd_frac; s.Th_floor = 1e-10
        s.loading_dep = loading_dep or (lambda th: np.ones_like(th))
        s.z_f = np.linspace(0.0, delta, n_z + 1)
        s.z_c = 0.5 * (s.z_f[:-1] + s.z_f[1:])
        s.shape = (n_z, s.n)                      # (space, field) - never (n,)
        s.div = construct_div(s.shape, s.z_f, nu=0)      # nu=0: Cartesian slab
        # The face operator is block diagonal with one n_c x n_c block per face,
        # and its sparsity never changes - only the numbers in it. Building the
        # CSR indices once and refilling `data` is ~10x faster than calling
        # scipy's block_diag inside every Newton iteration, and bit-identical.
        nf, nc = n_z + 1, s.n
        s._blk_cols = np.tile(
            (np.arange(nf)[:, None] * nc + np.arange(nc)[None, :]).repeat(nc, 0), 1).ravel()
        s._blk_ptr = np.arange(nf * nc + 1) * nc
        s._blk_shape = (nf * nc, nf * nc)

    # --- isotherm interface: IAST (eqs. 62-65) or multicomponent Langmuir (39)
    def load(s, p):
        p = np.atleast_2d(np.asarray(p, float))
        if s.iso_model == "iast":
            return s.iast.forward(p)
        den = 1 + np.sum([s.isos[i].bA * p[:, i] for i in range(s.n)], 0)
        return np.stack([s.Th_sat[i] * s.isos[i].bA * p[:, i] / den
                         for i in range(s.n)], -1)

    def gamma(s, Th):
        if s.iso_model == "iast":
            _, G = s.iast.inverse_and_gamma(Th, s.Th_sat, seed=getattr(s, "_seed", None))
            s._seed = s.iast.psi_last
            return G
        th = np.atleast_2d(Th) / s.Th_sat                              # eq. (39)
        return np.eye(s.n)[None] + th[:, :, None] / (1 - th.sum(-1))[:, None, None]

    def _M(s, Th, Th0, Thd):
        """Face-wise [Theta_sat][B]^-1[Gamma][Theta_sat]^-1, as a block matrix."""
        Thf = np.empty((s.n_z + 1, s.n))
        Thf[1:-1] = 0.5 * (Th[:-1] + Th[1:]); Thf[0] = Th0; Thf[-1] = Thd
        Thf = np.maximum(Thf, s.Th_floor)
        th = Thf / s.Th_sat
        G = s.gamma(Thf)
        B = B_matrix(th, s.D0 * s.loading_dep(th), s.Dij)
        M = np.linalg.solve(B, G) * s.Th_sat[None, :, None] / s.Th_sat[None, None, :]
        return sp.csr_array((M.ravel(), s._blk_cols, s._blk_ptr), shape=s._blk_shape)

    def _setup(s, p0, pd):
        p0 = np.asarray(p0, float)
        Th0 = s.load(p0[None, :])[0]
        # "the downstream partial pressures are maintained at vanishing values by
        # means of a sweep gas" (Section 4.1), as a fixed fraction of upstream.
        Thd = s.load((np.asarray(pd, float) if pd is not None else s.pd_frac * p0)[None, :])[0]
        # OUTWARD normal, a dTheta/dn + b Theta = d. Both faces are Dirichlet
        # (eqs. 57 and 58), so the same dict means the same thing at both ends.
        bc = ({"a": 0.0, "b": 1.0, "d": Th0.reshape(1, -1)},
              {"a": 0.0, "b": 1.0, "d": Thd.reshape(1, -1)})
        grad, gb = construct_grad(s.shape, s.z_f, s.z_c, bc)
        return Th0, Thd, grad, np.asarray(gb.todense()).reshape(-1)

    def flux(s, Th, Th0, Thd, grad, gbc):
        return -(s._M(Th, Th0, Thd) @ (grad @ Th.reshape(-1) + gbc)).reshape(-1, s.n)

    def residual(s, Th, Th_old, dt, grad, gbc, Th0, Thd):
        Th = np.maximum(Th.reshape(s.n_z, s.n), s.Th_floor)
        r = -(s.div @ (s._M(Th, Th0, Thd) @ (grad @ Th.reshape(-1) + gbc)))
        if dt is not None:
            r = (Th.reshape(-1) - Th_old.reshape(-1)) / dt - r
        return r.reshape(s.n_z, s.n)

    def _newton(s, Th, Th_old, dt, grad, gbc, Th0, Thd, jac, tol=1e-10, maxit=30):
        """Damped Newton on L = ln(Theta).

        Working in the logarithm does three things at once: the loadings stay
        positive with no clip; NumJac's perturbation becomes relative, which
        matters because the loadings span six decades across the membrane; and
        the update is multiplicative. The line search rejects any state with
        sum(theta) >= 1 -- IAST has no solution above saturation and its inverse
        map diverges there, which is how an undamped Newton produces NaN.
        """
        def flog(L):
            return s.residual(np.exp(L), Th_old, dt, grad, gbc, Th0, Thd) / s.Th_sat

        L = np.log(Th)
        for k in range(maxit):
            r, J = jac(flog, L)
            dL = np.clip(spsolve(J.tocsc(), -r.reshape(-1)).reshape(s.n_z, s.n), -2.0, 2.0)
            r0 = np.max(np.abs(r)); lam = 1.0; Ln = L
            for _ in range(14):
                Lt = L + lam * dL
                if np.max(np.sum(np.exp(Lt) / s.Th_sat, -1)) < 1 - 1e-9:
                    if np.max(np.abs(flog(Lt))) < r0 or lam < 1e-5:
                        Ln = Lt; break
                lam *= 0.5
            L = Ln
            s.newton_res = float(np.max(np.abs(flog(L))))
            if np.max(np.abs(lam * dL)) < tol:
                break
        s.newton_its = k + 1
        return np.exp(L)

    def _steady_once(s, p0, pd, tol, maxit, Th_start):
        Th0, Thd, grad, gbc = s._setup(p0, pd)
        w = (s.z_c / s.delta)[:, None]
        Th = Th_start if Th_start is not None else np.exp((1 - w) * np.log(Th0)
                                                          + w * np.log(Thd))
        # NumJac((n_z, n_c), axes_diagonals=[0]): the residual reads NEIGHBOURING
        # cells (through the face fluxes), so the space axis is tridiagonal; the
        # field axis is coupled in full by [B]^-1[Gamma]. ndims = 2, so
        # axes_diagonals is meaningful here (on a 1-D shape it would not be).
        jac = NumJac((s.n_z, s.n), axes_diagonals=[0])
        Th = s._newton(Th, None, None, grad, gbc, Th0, Thd, jac, tol=tol, maxit=maxit)
        J = s.flux(Th, Th0, Thd, grad, gbc)
        # at steady state div N = 0, so the flux must be z-independent. This is
        # THE convergence assertion; the Newton iteration count never is.
        s.flux_spread = float(np.max(np.abs(J / J.mean(0) - 1.0)))
        return Th, J

    def steady(s, p0, pd=None, tol=1e-9, maxit=25, Th_start=None,
               spread_tol=1e-7, n_cont=4):
        """Steady state, with div N = 0 asserted rather than assumed.

        A cold start is tried first. If the flux is not z-independent to
        `spread_tol` the solve is repeated along a FIXED pressure-continuation
        path (two decades below the target, logarithmically), which is
        deterministic and reproducible - no adaptive path, no warm start carried
        between reported cases. Validation checks on a case where both routes
        converge that they land on the same answer.
        """
        Th, J = s._steady_once(p0, pd, tol, maxit, Th_start)
        s.continued = False
        if s.flux_spread > spread_tol:
            s.continued = True
            Th = None
            for f_ in np.logspace(-2, 0, n_cont):
                Th, J = s._steady_once(np.asarray(p0, float) * f_, pd, tol, maxit, Th)
        return Th, J

    def transient(s, p0, t_out, pd=None, dt_frac=0.09, dt0=1e-6, tol=1e-11, maxit=40):
        """Backward Euler; inside each step the face coefficients are frozen and
        the linear system re-solved until the loadings stop moving. A step whose
        inner iteration does not converge is NEVER accepted - it is halved."""
        Th0, Thd, grad, gbc = s._setup(p0, pd)
        I = sp.eye_array(s.n_z * s.n, format="csr")
        Th = np.tile(Thd, (s.n_z, 1)); t = 0.0; rec = []
        s.nsteps = s.maxinner = s.cutbacks = 0
        for target in np.atleast_1d(np.asarray(t_out, float)):
            while t < target - 1e-15:
                step = min(max(dt0, dt_frac * max(t, dt0)), target - t)
                Th_old = Th.copy()
                while True:
                    Th_it = Th_old.copy(); ok = False
                    for k in range(1, maxit + 1):
                        Mf = s._M(Th_it, Th0, Thd)
                        A = (I / step - s.div @ Mf @ grad).tocsc()
                        b = Th_old.reshape(-1) / step + (s.div @ Mf) @ gbc
                        Tn = np.maximum(spsolve(A, b).reshape(s.n_z, s.n), s.Th_floor)
                        fs = np.sum(Tn / s.Th_sat, -1); bd = fs > 1 - 1e-9
                        if np.any(bd):
                            Tn[bd] *= (1 - 1e-9) / fs[bd, None]
                        d = np.max(np.abs(Tn - Th_it)) / max(1.0, np.max(Th_it))
                        Th_it = Tn
                        if d < tol:
                            ok = True; break
                    if ok:
                        break
                    step *= 0.25; s.cutbacks += 1
                    if step < 1e-14 * max(1.0, t):
                        raise RuntimeError("inner iteration failed")
                Th = Th_it; s.maxinner = max(s.maxinner, k); s.nsteps += 1; t += step
            rec.append((t, Th.copy(), s.flux(Th, Th0, Thd, grad, gbc)))
        return rec


def iso(component, table):
    g = params(component, table)
    return DSL(g["b_A"], g["Theta_sat_A"], g["b_B"], g["Theta_sat_B"])


C1, NC4 = iso("C1", "Table 1"), iso("nC4", "Table 1")
NC6, MP3 = iso("nC6", "Table 2"), iso("3MP", "Table 2")
DMB = iso("2,2DMB", "Table 2")
D_C1 = params("C1", "Section 4.1")["D_MS"]
D_NC4 = params("nC4", "Section 4.1")["D_MS"]
P_C1, P_NC4 = A47_VAL["p0_C1"], A47_VAL["p0_nC4"]
print(f"C1  : Theta_sat = {C1.Th_sat:4.1f}, Henry coefficient {C1.H:.5g} Pa^-1, "
      f"D_MS = {D_C1:.0e} m2/s, upstream {P_C1/1e3:.0f} kPa")
print(f"nC4 : Theta_sat = {NC4.Th_sat:4.1f}, Henry coefficient {NC4.H:.5g} Pa^-1, "
      f"D_MS = {D_NC4:.0e} m2/s, upstream {P_NC4/1e3:.0f} kPa")
print(f"nC6 / 3MP / 2,2DMB at 362 K: Theta_sat = {NC6.Th_sat:.1f} / "
      f"{MP3.Th_sat:.1f} / {DMB.Th_sat:.1f}")'''))

# ---------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. The mixture isotherm, and the review's Henry-limit selectivity

Section 4.2 states that for mixture loadings below 8 molecules per unit cell the
nC₄/C₁ sorption selectivity of the 95–5 mixture "is practically constant and
equals that calculated from the corresponding Henry coefficients, i.e. 2200".

That is two claims, and they can be checked separately: the **value**, which
follows from Table 1 alone, and the **plateau**, which needs IAST. The Henry
coefficient of a dual-site Langmuir is $H_i=\sum_s\Theta_{{\rm sat},s}b_s$, so
the value is arithmetic on four numbers."""))

cells.append(code('''S_HENRY = NC4.H / C1.H
S_PRINTED = PRINTED["sorp_S_henry"]
print("Henry coefficients from the Table 1 DSL rows, H_i = sum_s Theta_sat,s b_s:")
print(f"   H_C1  = {C1.H:.6g} Pa^-1      H_nC4 = {NC4.H:.6g} Pa^-1")
print(f"   S(Henry limit) = H_nC4/H_C1 = {S_HENRY:.1f}")
print(f"   the review prints                     {S_PRINTED:.0f}"
      f"    ({(S_HENRY/S_PRINTED-1)*100:+.1f} %)")

# The IAST solve must reproduce the Henry limit as P -> 0, from a completely
# different code path (eqs. 62-65 rather than two products and a quotient).
mix = IAST([C1, NC4])
y = np.array([0.95, 0.05])
Pt = np.logspace(-3, 8, 45)
Th = mix.forward(Pt[:, None] * y)
S_iast = (Th[:, 1] / Th[:, 0]) * (y[0] / y[1])
Thmix = Th.sum(1)
S_LIMIT_DEV = abs(S_iast[0] / S_HENRY - 1)
print(f"\\n   IAST at P_total = {Pt[0]:.0e} Pa gives S = {S_iast[0]:.1f}, "
      f"{S_LIMIT_DEV:.1e} from the Henry ratio")
print("   -- a genuinely different code path (eqs. 62-65) reaching the same limit,")
print("   so the 2200 discrepancy is not an arithmetic slip in the Henry formula.")

print("\\nAnd the plateau the review describes, computed:")
print(f"{'Theta_mix':>11}{'P_total/Pa':>13}{'S (IAST)':>11}{'/S_Henry':>10}")
for tgt in (0.01, 1.0, 4.0, 6.0, 8.0):
    k = np.searchsorted(Thmix, tgt)
    P_at = brentq(lambda lp: mix.forward((10 ** lp) * y)[0].sum() - tgt, -3, 8, xtol=1e-12)
    Th_at = mix.forward((10 ** P_at) * y)[0]
    S_at = (Th_at[1] / Th_at[0]) * (y[0] / y[1])
    print(f"{tgt:11.2f}{10**P_at:13.4g}{S_at:11.1f}{S_at/S_HENRY:10.3f}")
P_S8 = brentq(lambda lp: mix.forward((10 ** lp) * y)[0].sum()
              - PRINTED["sorp_Theta_plateau"], -3, 8, xtol=1e-12)
Th8 = mix.forward((10 ** P_S8) * y)[0]
S_AT_8 = float((Th8[1] / Th8[0]) * (y[0] / y[1]))
P_2200 = brentq(lambda lp: (lambda T: (T[1]/T[0])*(y[0]/y[1]))(
    mix.forward((10 ** lp) * y)[0]) - S_PRINTED, -3, 8, xtol=1e-12)
TH_AT_2200 = float(mix.forward((10 ** P_2200) * y)[0].sum())
print(f"\\nAt the review's own plateau edge, Theta_mix = "
      f"{PRINTED['sorp_Theta_plateau']:.0f}, IAST gives S = {S_AT_8:.0f} -")
print(f"{(1-S_AT_8/S_HENRY)*100:.0f} % below the Henry limit, though only "
      f"{np.log10(S_HENRY/S_AT_8):.2f} decades on an axis where S")
print(f"later falls by {np.log10(S_HENRY/S_iast[-1]):.1f} decades. So "
      f'"practically constant" is fair on a log scale;')
print(f"the printed 2200 is not the plateau VALUE - IAST passes through it at "
      f"Theta_mix = {TH_AT_2200:.2f},")
print("which is inside the plateau rather than at its Henry-limit end.")'''))

cells.append(md(r"""**So the printed 2200 does not follow from Table 1, and the next question is
whether Table 1 is the problem.** The obvious explanation would be a
transcription slip in one of the four C1/nC4 Table 1 parameters. That
explanation can be tested, because those two rows are the ones `A4.7` certifies
against a printed target: they reproduce the single-component flux ratio to
0.17 %. Solving for the parameter value that would make the Henry ratio 2200,
and then asking what that value does to the flux ratio, closes the question."""))

cells.append(code('''FLUX_C1, FLUX_NC4 = A47_VAL["flux_C1"], A47_VAL["flux_nC4"]
RATIO_PRINTED = FLUX_C1 / FLUX_NC4


def flux_ratio(c1=None, nc4=None):
    """A4.7's Section 4.1 check: the single-component flux ratio, in which the
    membrane thickness and the framework density both cancel exactly."""
    c1 = c1 or C1; nc4 = nc4 or NC4
    return (D_C1 * c1.psi(P_C1)) / (D_NC4 * nc4.psi(P_NC4))


BASE_RATIO = flux_ratio()
print(f"A4.7's certification of the two rows, recomputed here from the same CSV:")
print(f"   flux ratio from eqs. (12) and (32) : {BASE_RATIO:.4f}")
print(f"   from A4.7's stored 34 and 3.1      : {RATIO_PRINTED:.4f}   "
      f"({(BASE_RATIO/RATIO_PRINTED-1)*100:+.3f} %)")
print("   (the printed fluxes carry two significant figures, so the ratio itself")
print("    carries about 3 % of rounding - that bounds what this check can see)\\n")

FIELDS = [("b_A", 0), ("Theta_sat_A", 1), ("b_B", 2), ("Theta_sat_B", 3)]
BASE = {"C1": [C1.bA, C1.tA, C1.bB, C1.tB], "nC4": [NC4.bA, NC4.tA, NC4.bB, NC4.tB]}
print("Every single-parameter change to the C1 or nC4 Table 1 row that would make")
print("the Henry ratio equal the printed 2200, and what it does to that check:\\n")
print(f"{'row':>16}{'as printed':>13}{'needed for 2200':>17}{'factor':>9}"
      f"{'flux ratio dev':>16}")
ESCAPES = []
for sp_name in ("C1", "nC4"):
    for fname, k in FIELDS:
        v0 = BASE[sp_name][k]
        if v0 <= 0:
            continue

        def henry_with(val, sp_name=sp_name, k=k):
            w = list(BASE[sp_name]); w[k] = val
            a = DSL(*w) if sp_name == "C1" else C1
            b = NC4 if sp_name == "C1" else DSL(*w)
            return b.H / a.H

        try:
            val = 10 ** brentq(lambda lv: henry_with(10 ** lv) - S_PRINTED,
                               np.log10(v0) - 4, np.log10(v0) + 4, xtol=1e-14)
        except ValueError:
            print(f"{sp_name + ' ' + fname:>16}{v0:13.4g}"
                  f"{'no solution within 4 decades':>42}")
            continue
        w = list(BASE[sp_name]); w[k] = val
        fr = (flux_ratio(c1=DSL(*w)) if sp_name == "C1" else flux_ratio(nc4=DSL(*w)))
        dev = fr / RATIO_PRINTED - 1
        ESCAPES.append(abs(dev))
        print(f"{sp_name + ' ' + fname:>16}{v0:13.4g}{val:17.4g}{val/v0:9.3f}"
              f"{dev*100:+15.2f} %")
TIGHTEST = float(min(ESCAPES))
OBSERVED = abs(BASE_RATIO / RATIO_PRINTED - 1)
print(f"\\nThe least damaging of them still moves the flux ratio by "
      f"{TIGHTEST*100:.1f} %, which is")
print(f"{TIGHTEST/OBSERVED:.0f}x the {OBSERVED*100:.2f} % actually observed and outside "
      "the ~3 % the two-figure")
print("printed fluxes allow. So no single-character slip in the two Table 1 rows")
print("explains 2200 without breaking the one check that certifies them.")
print("\\nWHAT THIS DOES AND DOES NOT SETTLE. It settles that the value stored in")
print("this repository's copy of Table 1 is not the explanation. It does NOT say")
print("where 2200 comes from: the review's Fig. 24 is drawn against CBMC")
print("fugacities rather than against the DSL fit, so a Henry coefficient taken")
print("from the simulation data instead of from the Table 1 fit is a live")
print("possibility that nothing printed in the review can decide. Recorded as a")
print("printed-value discrepancy, not as an error.")'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2))
axes[0].loglog(Thmix, S_iast, lw=2.0, color="tab:blue", label="IAST (eqs. 62-65)")
axes[0].axhline(S_HENRY, color="tab:green", lw=1.3, ls="--",
                label=f"Henry ratio from Table 1 = {S_HENRY:.0f}")
axes[0].axhline(S_PRINTED, color="tab:red", lw=1.3, ls=":",
                label=f"the review prints {S_PRINTED:.0f}")
axes[0].axvline(PRINTED["sorp_Theta_plateau"], color="0.4", lw=1.0)
axes[0].annotate(r"$\Theta_{\rm mix}=8$", (8.4, 60), fontsize=8, color="0.4")
axes[0].set(xlabel=r"mixture loading $\Theta_{\rm mix}$ (molecules per unit cell)",
            ylabel=r"sorption selectivity $S$ (eq. 61)", xlim=(1e-3, 13),
            title="95-5 C1/nC4 in MFI at 300 K")
axes[0].legend(fontsize=8)

axes[1].semilogx(Pt, Th[:, 0], lw=2.0, color="tab:blue", label=r"C1, IAST")
axes[1].semilogx(Pt, Th[:, 1], lw=2.0, color="tab:red", label=r"nC4, IAST")
axes[1].semilogx(Pt, [C1.theta(0.95 * p) for p in Pt], lw=1.2, ls="--",
                 color="tab:blue", label="C1, pure at the same partial pressure")
axes[1].semilogx(Pt, [NC4.theta(0.05 * p) for p in Pt], lw=1.2, ls="--",
                 color="tab:red", label="nC4, pure at the same partial pressure")
axes[1].set(xlabel="total pressure (Pa)", ylabel=r"$\Theta$ (molecules per unit cell)",
            title="the size-entropy effect: nC4 is pushed out near saturation")
axes[1].legend(fontsize=8, loc="upper left")
fig.tight_layout(); plt.show()

k_peak = int(np.argmax(Th[:, 1]))
print(f"The nC4 loading peaks at {Th[k_peak,1]:.3f} molecules per unit cell at a total")
print(f"pressure of {Pt[k_peak]:.2g} Pa and then DECLINES as the pressure rises further,")
print(f"reaching {Th[-1,1]:.3f} at {Pt[-1]:.0g} Pa while C1 climbs from "
      f"{Th[k_peak,0]:.3f} to {Th[-1,0]:.3f}.")
print("The review describes exactly that ('increasing the total system pressure")
print("beyond 5 MPa leads to a decline in the loading of nC4!'), and it is the")
print("whole reason the multicomponent Langmuir cannot be used: no Langmuir")
print("isotherm has a non-monotonic component loading.")'''))

cells.append(md(r"""### 2. The hexane isomers sorb the other way round

The same machinery on Table 2's parameters at 362 K. Here the review states a
*threshold*: separation works when the upstream mixture loading exceeds 4
molecules per unit cell, "at 362 K this corresponds to a total upstream pressure
in excess of 1 kPa"."""))

cells.append(code('''hex_mix = IAST([NC6, MP3])
Ph = np.logspace(0, 6, 55)
Thh = hex_mix.forward(0.5 * Ph[:, None] * np.ones(2))
S_hex = Thh[:, 0] / Thh[:, 1]
P_TH4 = 10 ** brentq(lambda lp: hex_mix.forward(0.5 * (10 ** lp) * np.ones(2))[0].sum()
                     - 4.0, 0, 6, xtol=1e-12)
TH_AT_1KPA = float(hex_mix.forward(0.5 * PRINTED["hex_P_theta4"] * np.ones(2))[0].sum())
print(f"IAST puts Theta_mix = 4 at a total pressure of {P_TH4:.1f} Pa "
      f"({P_TH4/1e3:.3f} kPa).")
print(f"The review states 'in excess of {PRINTED['hex_P_theta4']/1e3:.0f} kPa'; at "
      f"exactly that pressure IAST gives Theta_mix = {TH_AT_1KPA:.2f}.")
print("The review's statement is therefore CORRECT but loose by a factor of "
      f"{PRINTED['hex_P_theta4']/P_TH4:.1f};")
print("it is a sufficient condition, not the threshold, and is scored as such.\\n")
print(f"{'P_total/kPa':>12}{'Theta_nC6':>11}{'Theta_3MP':>11}{'Theta_mix':>11}"
      f"{'S (nC6/3MP)':>13}")
for P in (0.1e3, 0.5e3, 1e3, 4e3, 15e3, 66e3):
    T = hex_mix.forward(0.5 * P * np.ones(2))[0]
    print(f"{P/1e3:12.2f}{T[0]:11.4f}{T[1]:11.4f}{T.sum():11.4f}{T[0]/T[1]:13.3f}")
print("\\nBelow Theta_mix = 4 the two isomers sorb almost equally; above it the")
print("branched one is squeezed out of the channel interiors. That is the")
print("configurational-entropy effect the separation is built on, and it is")
print("produced here by IAST from the pure-component fits alone.")'''))

# ------ section 3: the membrane
cells.append(md(r"""### 3. The membrane, and where the 1.73 goes

Now eqs. (57)–(59). Four runs: the 95–5 and 50–50 mixtures, each with the
exchange coefficients from eq. (43) and with $\eth_{ij}\to\infty$. Every
selectivity below is independent of the membrane thickness, the framework
density and the absolute size of the diffusivities."""))

cells.append(code('''def SP(J, p0):
    """eq. (60) as the review uses it: (N_2/N_1) divided by (p_2,0/p_1,0)."""
    Jm = J.mean(0)
    return float((Jm[1] / Jm[0]) / (p0[1] / p0[0]))


FEEDS = {"95-5": np.array([P_C1, P_NC4]), "50-50": np.array([50e3, 50e3])}
RUNS = {}
print(f"{'feed':>7}{'exchange':>26}{'Newton':>8}{'flux spread':>13}"
      f"{'S_P':>10}{'printed':>10}{'dev':>9}")
for feed, p0 in FEEDS.items():
    for mode, lab in (("vignes", "finite, eq. (43)"), ("inf", "infinite, eq. (47)")):
        m = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4],
                     n_z=40, Dij=mode)
        Th_m, J = m.steady(p0)
        s_p = SP(J, p0)
        key = f"mix_SP_{feed.replace('-', '')}_{'finite' if mode=='vignes' else 'inf'}"
        pr = PRINTED[key]
        RUNS[(feed, mode)] = (m, Th_m, J, s_p)
        print(f"{feed:>7}{lab:>26}{m.newton_its:8d}{m.flux_spread:13.1e}"
              f"{s_p:10.2f}{pr:10.1f}{(s_p/pr-1)*100:8.2f}%")

SP_955 = RUNS[("95-5", "vignes")][3]
SP_955_INF = RUNS[("95-5", "inf")][3]
SP_5050 = RUNS[("50-50", "vignes")][3]
SP_5050_INF = RUNS[("50-50", "inf")][3]
SP_PURE = A47_VAL["sp_pure"]
print(f"\\nA4.7's stored pure-component expectation is S_P = {SP_PURE}, and eq. (60)")
print(f"on its stored fluxes gives {(FLUX_NC4/FLUX_C1)/(P_NC4/P_C1):.3f} - arithmetic on two printed")
print("numbers, which tests the transcription and nothing else.")
print(f"The mixture value is {SP_955/SP_PURE:.0f}x larger, and removing the "
      "sorbate-sorbate")
print(f"friction alone takes it back down by a factor of "
      f"{SP_955/SP_955_INF:.0f}, to {SP_955_INF:.1f}.")'''))

cells.append(md(r"""**The absolute fluxes need a membrane thickness, which the review never prints.**
`A4.7` obtained $\delta\approx 40\ \mu$m by inverting each of the two
*single-component* steady fluxes with the only MFI density the review prints,
1800 kg m⁻³ from the caption of Fig. 40. It is reconstructed here from the same
stored rows and used to put the four mixture fluxes into mmol m⁻² s⁻¹.

**What that comparison tests, and what it cannot — stated before the numbers,
because an earlier version of this page got it wrong in exactly the way
[`docs/handoff.md`](../../docs/handoff.md) warns about.** That version claimed
the four mixture fluxes were an out-of-sample test of `A4.7`'s $\delta$ on the
$1/\delta$ power. **They are not. $\delta$ cancels identically.** Substituting
the reconstruction
$\delta=\rho_{\rm uc}\eth_{\rm C1}\psi_{\rm C1}(95\ {\rm kPa})/(N_A N_{\rm C1})$
into the conversion $N_i=\rho_{\rm uc}J_i/(N_A\delta)$ leaves

$$N_i \;=\; N_{\rm C1}^{\rm stored}\,\frac{J_i}{\eth_{\rm C1}\,\psi_{\rm C1}(95\ {\rm kPa})},$$

with no $\delta$, no $\rho_{\rm MFI}$, no unit-cell mass and no Avogadro
constant left in it. A four-fold error in the framework density would leave all
four "absolute" fluxes bit-identical, and the cell below demonstrates that rather
than asserting it. So this comparison **is** the $\delta$-free flux-ratio
comparison printed ten lines below it, rescaled by the stored 34 — the two
deviations differ only by the $n=1$ solver's discretisation error against the
analytic $\psi$, and both are printed so that the coincidence is visible.

What survives is still worth having, and it is about the *physics* rather than
the thickness: the pure C1 flux is **guaranteed** to come out at the stored 34
(that is what $\delta$ was fitted to) and the pure nC4 flux carries `A4.7`'s
0.17 % flux-ratio deviation and nothing more, whereas the four **mixture** fluxes
run through IAST, a matrix $[\Gamma]$ and finite $\eth_{ij}$, none of which
entered the reconstruction, and none of which is guaranteed to land anywhere.
The one quantity on this page that does carry $\delta$ is the transient peak
*time*, and it is taken up in the next section."""))

cells.append(code('''NA = 6.02214076e23
M_CELL = (96 * 28.0855 + 192 * 15.999) / NA * 1e-3      # kg per Si96O192 unit cell
RHO_UC = A47_VAL["rho_MFI"] / M_CELL                    # unit cells per m3
DELTA = float(RHO_UC * D_C1 * C1.psi(P_C1) / NA / (FLUX_C1 * 1e-3))
DELTA_NC4 = float(RHO_UC * D_NC4 * NC4.psi(P_NC4) / NA / (FLUX_NC4 * 1e-3))
print(f"A4.7's reconstruction, recomputed from the same stored rows:")
print(f"   from the stored C1 flux  : {DELTA*1e6:.3f} um")
print(f"   from the stored nC4 flux : {DELTA_NC4*1e6:.3f} um   "
      f"({abs(DELTA_NC4/DELTA-1)*100:.2f} % apart - this IS the flux-ratio check")
print("     inverted, not a second one)\\n")


def mmol(J, L=1.0):
    """Solver flux -> mmol m^-2 s^-1. `L` is the thickness the grid was built
    on: the delta-free runs use L = 1 and the dimensional ones use L = DELTA,
    and the physical flux scales as L/delta."""
    return RHO_UC * np.asarray(J) * (L / DELTA) / NA * 1e3


print(f"{'case':>34}{'computed':>11}{'printed':>10}{'dev':>9}")
FLUX_ROWS = []
for feed, mode, comp, key in (("95-5", "vignes", 0, "mix_N_C1_finite"),
                              ("95-5", "vignes", 1, "mix_N_nC4_finite"),
                              ("95-5", "inf", 0, "mix_N_C1_inf"),
                              ("95-5", "inf", 1, "mix_N_nC4_inf")):
    J = RUNS[(feed, mode)][2].mean(0)
    v = float(mmol(J)[comp]); pr = PRINTED[key]
    lab = f"{'C1' if comp == 0 else 'nC4'}, {feed}, " \\
          f"{'finite' if mode == 'vignes' else 'infinite'} Dij"
    FLUX_ROWS.append((lab, pr, v))
    print(f"{lab:>34}{v:11.4f}{pr:10.2f}{(v/pr-1)*100:8.2f}%")
FLUX_WORST = max(abs(v / p - 1) for _, p, v in FLUX_ROWS)
print(f"\\nWorst of the four mixture fluxes: {FLUX_WORST*100:.2f} %, against printed")
print("values carrying two or three significant figures. What that does NOT test")
print("is delta, and the demonstration is below.")

# ---- POWERLESS-CHECK DISCLOSURE: delta cancels out of all four fluxes --------
print("\\nDELTA, rho_MFI, THE UNIT-CELL MASS AND N_A ALL CANCEL. Swing the framework")
print("density over a factor of four - which moves the reconstructed thickness")
print("over the same factor - and every one of the four fluxes is BIT-IDENTICAL:")
print(f"\\n{'rho_MFI / kg m^-3':>19}{'delta / um':>13}{'N_C1 (finite)':>18}"
      f"{'N_nC4 (finite)':>18}{'worst of 4':>12}")
J_FIN = RUNS[("95-5", "vignes")][2].mean(0)
J_INF = RUNS[("95-5", "inf")][2].mean(0)
BASE4 = np.concatenate([mmol(J_FIN), mmol(J_INF)])
for rho in (A47_VAL["rho_MFI"], 1200.0, 3600.0, 900.0):
    ruc = rho / M_CELL
    dl = float(ruc * D_C1 * C1.psi(P_C1) / NA / (FLUX_C1 * 1e-3))
    four = np.concatenate([ruc * np.asarray(J) * (1.0 / dl) / NA * 1e3
                           for J in (J_FIN, J_INF)])
    print(f"{rho:19.0f}{dl*1e6:13.4f}{four[0]:18.12f}{four[1]:18.12f}"
          f"{np.max(np.abs(four - BASE4)):12.1e}")
CLOSED_FORM = FLUX_C1 * np.asarray(J_FIN) / (D_C1 * C1.psi(P_C1))
DELTA_CANCELS = float(np.max(np.abs(CLOSED_FORM / mmol(J_FIN) - 1)))
print(f"\\n   because the conversion collapses algebraically to "
      f"N_i = {FLUX_C1:g} J_i / (D_C1 psi_C1),")
print(f"   verified against the code path above to {DELTA_CANCELS:.1e}.")
print("   So the four-flux comparison is the delta-free flux-RATIO comparison")
print("   below, multiplied by the stored 34. It is a real test of the MIXTURE")
print("   physics - IAST, [Gamma] and finite Dij all enter J_i - and no test at")
print("   all of the membrane thickness. Kept and labelled rather than deleted,")
print("   because deleting it would hide that the two rows are one check.")

# the delta-free version of the same comparison, for readers who do not want to
# accept a reconstruction at all
pure_C1 = Membrane([C1], [C1.Th_sat], [D_C1], n_z=120)
pure_NC4 = Membrane([NC4], [NC4.Th_sat], [D_NC4], n_z=120)
JP_C1 = float(pure_C1.steady(np.array([P_C1]))[1].mean(0)[0])
JP_NC4 = float(pure_NC4.steady(np.array([P_NC4]))[1].mean(0)[0])
JM = RUNS[("95-5", "vignes")][2].mean(0)
R_C1 = float(JM[0]) / JP_C1
R_NC4 = float(JM[1]) / JP_NC4
print(f"\\nThe same statement with delta removed entirely, as flux ratios:")
print(f"   N_C1(mixture)/N_C1(pure)  : computed {R_C1:.5f}   printed "
      f"{PRINTED['mix_N_C1_finite']/FLUX_C1:.5f}   "
      f"({(R_C1/(PRINTED['mix_N_C1_finite']/FLUX_C1)-1)*100:+.2f} %)")
print(f"   N_nC4(mixture)/N_nC4(pure): computed {R_NC4:.5f}   printed "
      f"{PRINTED['mix_N_nC4_finite']/FLUX_NC4:.5f}   "
      f"({(R_NC4/(PRINTED['mix_N_nC4_finite']/FLUX_NC4)-1)*100:+.2f} %)")
print("   Methane's flux collapses by a factor of "
      f"{1/R_C1:.0f} in the mixture while butane's")
print(f"   is essentially unchanged ({(R_NC4-1)*100:+.1f} %). That asymmetry is the "
      "whole effect.")
RATIO_FREE_WORST = max(abs(R_C1 / (PRINTED["mix_N_C1_finite"] / FLUX_C1) - 1),
                       abs(R_NC4 / (PRINTED["mix_N_nC4_finite"] / FLUX_NC4) - 1))
print(f"\\n   And there is the coincidence, in the open: the four-flux worst")
print(f"   deviation is {FLUX_WORST:.7f} and this delta-free one is "
      f"{RATIO_FREE_WORST:.7f}. They")
print(f"   differ by {abs(FLUX_WORST-RATIO_FREE_WORST):.1e} in absolute terms, which is the "
      "n = 1 SOLVER")
print("   against the analytic psi and nothing else - i.e. they are the same")
print("   check reported twice, and both are kept.")
# the n = 1 collapse of eq. (45) to eq. (12): the solver run with one species
# must reproduce the closed-form integral A4.7 uses
ANA_C1 = D_C1 * C1.psi(P_C1)
ANA_NC4 = D_NC4 * NC4.psi(P_NC4)
COLLAPSE_1 = max(abs(JP_C1 / ANA_C1 - 1), abs(JP_NC4 / ANA_NC4 - 1))
print(f"\\n   (eq. 45 run with n = 1 against the closed-form eq. (12) integral: "
      f"{COLLAPSE_1:.1e})")'''))

cells.append(md(r"""#### The transient, and a peak that arrives before the butane does

The review reports that during the initial transience the methane flux
overshoots — 0.44 mmol m⁻² s⁻¹ at t = 0.73 s with finite exchange, 17.2 with
$\eth_{ij}\to\infty$ — because methane fills the membrane first and is then
displaced by the more strongly adsorbed butane.

**The peak time is the only quantity on this page that carries $\delta$ at all**,
and it carries it as $\delta^2$ where a steady flux carries it as $1/\delta$. The
four mixture fluxes above look like a second, independent handle on the
reconstruction and are not one — $\delta$ cancels out of them identically — so
this section is where `A4.7`'s $\approx 40\ \mu$m is tested, and it is the only
place it is tested.

**Which means the transient has to be resolved in space, and it is a different
resolution problem from the steady one.** The grid study in *Validation* refines
the steady selectivity, where the loading profile is smooth. The transient
carries a butane front, and its steepness is set by the site-A affinity at the
upstream face — printed below. Refining the steady solution says nothing about
it, so the grid is refined here, under the transient, before any number from it
is reported."""))

cells.append(code('''# Output times: a log sweep for the whole transient, plus a dense linear window
# through the peak so that its POSITION is not set by the output grid spacing.
t_out = np.unique(np.concatenate([np.logspace(-3, 3.3, 48),
                                  np.linspace(0.45, 1.15, 15)]))
t_pk = t_out[t_out <= 3.0]          # far enough past the peak to bracket it
NZ_PROD = 120                       # production grid, justified by the study below


def locate_peak(t, yv):
    """Peak of a curve, by a parabola in log t through the three points around
    the grid maximum, with the vertex clamped inside that bracket. Deterministic
    - no warm start, no continuation, no random seed."""
    i = int(np.argmax(yv))
    lt = np.log(t[i - 1:i + 2])
    c = np.polyfit(lt, yv[i - 1:i + 2], 2)
    tp = float(np.clip(np.exp(-c[1] / (2 * c[0])), t[i - 1], t[i + 1]))
    return tp, float(np.polyval(c, np.log(tp))), i


def run_transient(n_z, mode, dt_frac, tt):
    m = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4],
                 n_z=n_z, Dij=mode, delta=DELTA)
    rec = m.transient(np.array([P_C1, P_NC4]), tt, dt_frac=dt_frac)
    N = np.array([mmol(r[2][-1], DELTA) for r in rec])
    return m, rec, N


print("WHY THE STEADY GRID STUDY DOES NOT COVER THIS.\\n")
for nm, iso_, p0_ in (("C1 ", C1, P_C1), ("nC4", NC4, P_NC4)):
    bp = iso_.bA * p0_
    print(f"   {nm}: site-A affinity at the upstream face, b_A p = {bp:8.2f}"
          f"   -> that site is {100*bp/(1+bp):5.1f} % full")
print("   Both are essentially empty at the swept face. nC4 is the problem: a")
print("   site that is nearly saturated at one face and empty at the other has")
print("   an almost rectangular isotherm, and an almost rectangular isotherm")
print("   makes a SELF-SHARPENING front while it fills. The steady profile has")
print("   no such front; the transient does, and it has to be resolved.")

# --- production runs, full time range, at the grid the study below justifies
TRANS = {}
for mode in ("vignes", "inf"):
    m, rec, N = run_transient(NZ_PROD, mode, 0.05, t_out)
    TRANS[mode] = (m, N, rec)
    print(f"\\n   {mode:7s} at n_z = {NZ_PROD}: {m.nsteps} steps, at most "
          f"{m.maxinner} inner iterations, {m.cutbacks} cutbacks")

PEAKS = {}
for mode in ("vignes", "inf"):
    PEAKS[mode] = locate_peak(t_out, TRANS[mode][1][:, 0])

# how sharp that front actually is, on the production run. Measured as the
# largest cell-to-cell jump in the nC4 loading, as a fraction of its value at
# the upstream face - a grid-relative measure of how much structure one cell
# has to carry.
rec_v = TRANS["vignes"][2]
Th_face = float(TRANS["vignes"][0].load(np.array([[P_C1, P_NC4]]))[0][1])
JUMPS = [(t, float(np.max(np.abs(np.diff(Th[:, 1]))) / Th_face)) for t, Th, _ in rec_v]
T_WORST, FRONT_T = max(JUMPS, key=lambda r: r[1])
FRONT_PK = JUMPS[PEAKS["vignes"][2]][1]
FRONT_S = JUMPS[-1][1]
print(f"\\n   Measured on the n_z = {NZ_PROD} run, as the largest cell-to-cell jump in")
print(f"   the nC4 loading over the upstream face value:")
print(f"      worst over the whole transient (t = {T_WORST:.3g} s) : {FRONT_T*100:6.2f} %")
print(f"      at the peak time                            : {FRONT_PK*100:6.2f} %")
print(f"      at steady state                             : {FRONT_S*100:6.2f} %")
print(f"   {FRONT_T/FRONT_S:.0f}x at its worst. Early on one cell carries most of the whole")
print("   loading step, because the front starts AT the face; by the time of the")
print("   peak it has spread. The steady profile never has such a cell.")
print("   That is the mechanism. THE EVIDENCE is the grid study below, and it is")
print("   blunter: the steady selectivity is converged to 4e-05 at n_z = 40,")
print("   while the peak time is still 10 % out at n_z = 30.")'''))

cells.append(code('''print("GRID STUDY UNDER THE TRANSIENT. Only n_z changes; t_out, dt_frac = 0.05\\n"
      "and locate_peak are the ones used for the production run.\\n")
GRID_T = {}
for mode in ("vignes", "inf"):
    lvls = (30, 60, 240) if mode == "vignes" else (30, 60)
    rows = []
    for nz in lvls:
        _, _, N = run_transient(nz, mode, 0.05, t_pk)
        tp, yp, _ = locate_peak(t_pk, N[:, 0])
        rows.append((nz, tp, yp))
    rows.append((NZ_PROD, PEAKS[mode][0], PEAKS[mode][1]))   # the production run
    GRID_T[mode] = sorted(rows)

for mode, kf, kt in (("vignes", "mix_peak_C1_finite", "mix_peak_t_finite"),
                     ("inf", "mix_peak_C1_inf", None)):
    lab = "finite Dij (eq. 43)" if mode == "vignes" else "Dij -> infinity (eq. 47)"
    print(f"   {lab}")
    print(f"   {'n_z':>6}{'peak N_C1':>12}{'vs printed':>12}{'t_peak/s':>11}"
          f"{'vs printed':>12}")
    for nz, tp, yp in GRID_T[mode]:
        tcol = (f"{(tp/PRINTED[kt]-1)*100:+11.2f}%" if kt else f"{'-':>12}")
        star = "  <-- production" if nz == NZ_PROD else ""
        print(f"   {nz:6d}{yp:12.4f}{(yp/PRINTED[kf]-1)*100:+11.2f}%{tp:11.4f}"
              f"{tcol}{star}")
    print()

# Richardson on the finite-exchange peak time, which has four levels
gv = {nz: (tp, yp) for nz, tp, yp in GRID_T["vignes"]}
T_ORDER = float(np.log2(abs(gv[30][0] - gv[60][0]) / abs(gv[60][0] - gv[NZ_PROD][0])))
T_ORDER2 = float(np.log2(abs(gv[60][0] - gv[NZ_PROD][0]) / abs(gv[NZ_PROD][0] - gv[240][0])))
# f_inf = f_fine + (f_fine - f_coarse)/(r^p - 1), fine = 240, coarse = NZ_PROD
T_RICH = float(gv[240][0] + (gv[240][0] - gv[NZ_PROD][0]) / (2 ** T_ORDER2 - 1))
T_GRID_ERR = abs(gv[NZ_PROD][0] / T_RICH - 1)
print(f"   Observed order on t_peak: {T_ORDER:.2f} then {T_ORDER2:.2f}; Richardson from")
print(f"   the last pair extrapolates to {T_RICH:.4f} s, so the production n_z = "
      f"{NZ_PROD}")
print(f"   carries {T_GRID_ERR*100:.2f} % of spatial error on the peak time - "
      f"{abs(gv[NZ_PROD][0]/PRINTED['mix_peak_t_finite']-1)/T_GRID_ERR:.0f}x smaller than")
print(f"   its disagreement with the review, so it does not explain it.")
print(f"\\n   THE COARSE GRID FLATTERED BOTH NUMBERS, IN OPPOSITE DIRECTIONS. At")
print(f"   n_z = 30 the peak time is {gv[30][0]:.4f} s, "
      f"{(gv[30][0]/PRINTED['mix_peak_t_finite']-1)*100:+.1f} % from the printed "
      f"{PRINTED['mix_peak_t_finite']}, and the")
print(f"   peak height {gv[30][1]:.4f}, {(gv[30][1]/PRINTED['mix_peak_C1_finite']-1)*100:+.1f} %"
      f" from the printed {PRINTED['mix_peak_C1_finite']}. Refined, the time moves")
print(f"   AWAY to {(gv[NZ_PROD][0]/PRINTED['mix_peak_t_finite']-1)*100:+.1f} % and the height "
      f"CROSSES to {(gv[NZ_PROD][1]/PRINTED['mix_peak_C1_finite']-1)*100:+.1f} %. An earlier version of")
print("   this page reported the coarse pair and explained the height with")
print('   "backward Euler damps a maximum, so both are low". The dominant error')
print("   was spatial, both heights are now HIGH, and that explanation is")
print("   withdrawn rather than reworded.")'''))

cells.append(code('''print("TIME-STEP STUDY, AT THE PRODUCTION GRID. Three levels, not two - a")
print("two-level step study is what let the coarse grid pass unnoticed.\\n")
i_pk = PEAKS["vignes"][2]
GRID_SPACING = t_out[i_pk + 1] / t_out[i_pk] - 1
DT_ROWS = [(0.05, PEAKS["vignes"][0], PEAKS["vignes"][1])]
for dtf in (0.025, 0.0125):
    _, _, N_f = run_transient(NZ_PROD, "vignes", dtf, t_pk)
    DT_ROWS.append((dtf,) + locate_peak(t_pk, N_f[:, 0])[:2])
print(f"   output grid spacing at the peak: {GRID_SPACING*100:.1f} %\\n")
print(f"   {'dt_frac':>9}{'peak N_C1':>12}{'vs printed':>12}{'t_peak/s':>11}"
      f"{'vs printed':>12}")
for dtf, tp, yp in DT_ROWS:
    star = "  <-- production" if dtf == 0.05 else ""
    print(f"   {dtf:9.4f}{yp:12.4f}"
          f"{(yp/PRINTED['mix_peak_C1_finite']-1)*100:+11.2f}%{tp:11.4f}"
          f"{(tp/PRINTED['mix_peak_t_finite']-1)*100:+11.2f}%{star}")
TP_F, YP_F = DT_ROWS[1][1], DT_ROWS[1][2]
DT_PEAK_SHIFT = abs(YP_F / PEAKS["vignes"][1] - 1)
DT_TIME_SHIFT = abs(TP_F / PEAKS["vignes"][0] - 1)
RY = abs(DT_ROWS[0][2] - DT_ROWS[1][2]) / abs(DT_ROWS[1][2] - DT_ROWS[2][2])
RT = abs(DT_ROWS[0][1] - DT_ROWS[1][1]) / abs(DT_ROWS[1][1] - DT_ROWS[2][1])
print(f"\\n   successive-difference ratios: height {RY:.2f}, time {RT:.2f}. The height")
print(f"   is consistent with the first order backward Euler has; the TIME is not,")
print(f"   so no extrapolation in dt is offered for it - only the direction, which")
print("   is unambiguous. BOTH refinements push the same way: the height further")
print(f"   ABOVE the printed {PRINTED['mix_peak_C1_finite']} and the time further BELOW "
      f"the printed {PRINTED['mix_peak_t_finite']} s.")
print(f"   The reported values are the coarsest of the three, dt_frac = 0.05, so")
print("   the disagreement below is if anything understated.")

T_PEAK_DEV = PEAKS["vignes"][0] / PRINTED["mix_peak_t_finite"] - 1
DELTA_FOR_PEAK = float(np.sqrt(1 / (1 + T_PEAK_DEV)) - 1)
print(f"\\nWHAT THAT DOES TO A4.7'S RECONSTRUCTED THICKNESS. Converged, the peak")
print(f"time is {PEAKS['vignes'][0]:.4f} s against the printed "
      f"{PRINTED['mix_peak_t_finite']} s, {T_PEAK_DEV*100:+.1f} %. Since")
print(f"t_peak scales as delta^2, closing that gap needs delta "
      f"{DELTA_FOR_PEAK*100:+.1f} %, i.e.")
print(f"{DELTA*1e6:.1f} um -> {DELTA*(1+DELTA_FOR_PEAK)*1e6:.1f} um; at the finest step tried it "
      f"is {(np.sqrt(PRINTED['mix_peak_t_finite']/DT_ROWS[-1][1])-1)*100:+.1f} %. delta was")
print(f"inverted from a flux printed to two figures, and half a unit on "
      f"{FLUX_C1:g} is only")
print(f"{0.5/FLUX_C1*100:.1f} %, so rounding does not cover it either.")
print("\\nAND THE PAGE DOES NOT GET TO CALL THAT RECONCILED. The earlier version")
print("argued that moving delta would break the steady-flux agreement, so there")
print("was nothing to reconcile. That argument is void twice over: the required")
print(f"move is {DELTA_FOR_PEAK*100:+.1f} % and not a fraction of a per cent, and the four "
      "printed")
print("fluxes do not contain delta at all, so nothing about them constrains it.")
print("This is an OPEN TENSION and is reported as one. It sits in some")
print("combination of:")
print(f"  * rho_MFI = {A47_VAL['rho_MFI']:g} kg m^-3, read off a figure caption and")
print("    constrained by nothing else on this page or on A4.7 - delta scales")
print("    with it directly, so a 5 % error there is a 5 % error here;")
print("  * the Si96O192 unit-cell mass, hard-coded here from atomic weights and")
print("    traceable to nothing printed in the review;")
print(f"  * the printed {FLUX_C1:g} mmol m^-2 s^-1 that delta was inverted from, "
      "two figures;")
print("  * the review's transient itself, which may not be the calculation")
print("    reproduced here - it prints no delta, no grid and no time step.")
print("Nothing on this page decides between them, and nothing was adjusted.")'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2))
for ax, mode, ttl in ((axes[0], "vignes", r"finite $\eth_{ij}$ (eq. 43)"),
                      (axes[1], "inf", r"$\eth_{ij}\to\infty$ (eq. 47)")):
    N = TRANS[mode][1]
    ax.loglog(t_out, N[:, 0], lw=2.0, color="tab:blue", label="C1")
    ax.loglog(t_out, N[:, 1], lw=2.0, color="tab:red", label="nC4")
    tp, yp = PEAKS[mode][0], PEAKS[mode][1]
    ax.plot([tp], [yp], "o", ms=7, mfc="none", mew=1.6, color="tab:blue")
    key = "mix_peak_C1_finite" if mode == "vignes" else "mix_peak_C1_inf"
    ax.axhline(PRINTED[key], color="tab:blue", lw=1.0, ls=":")
    if mode == "vignes":
        ax.axvline(PRINTED["mix_peak_t_finite"], color="0.4", lw=1.0, ls="--")
        ax.annotate(f"printed peak\n{PRINTED[key]} at "
                    f"{PRINTED['mix_peak_t_finite']} s", (1.0, 0.6), fontsize=8,
                    color="0.35")
        ax.axhline(PRINTED["mix_N_C1_finite"], color="tab:blue", lw=0.8, ls="-.")
        ax.axhline(PRINTED["mix_N_nC4_finite"], color="tab:red", lw=0.8, ls="-.")
    else:
        ax.axhline(PRINTED["mix_N_C1_inf"], color="tab:blue", lw=0.8, ls="-.")
        ax.axhline(PRINTED["mix_N_nC4_inf"], color="tab:red", lw=0.8, ls="-.")
    ax.set(xlabel="time (s)", ylabel=r"$N$ (mmol m$^{-2}$ s$^{-1}$)",
           ylim=(1e-3, 1e2), title=ttl)
    ax.legend(fontsize=8, loc="lower left")
fig.tight_layout(); plt.show()
print("Dash-dot lines are the review's printed steady-state fluxes; the dotted line")
print("is its printed peak. Left and right differ ONLY in whether the two sorbates")
print("feel each other.")'''))

# ------ section 4: measurement
cells.append(md(r"""### 4. Against measurement

Bakker measured the nC₄/C₁ permeation selectivity of MFI at both feed
compositions. The review quotes both, calls the 95–5 comparison "quite close",
and says nothing about the 50–50 one. Sweeping the upstream composition at
constant total pressure puts all of it on one axis."""))

cells.append(code('''ys = np.array([0.99, 0.95, 0.80, 0.50, 0.20, 0.05])
SWEEP = {}
for mode in ("vignes", "inf"):
    Th_prev = None; row = []
    for yy in ys:
        m = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4],
                     n_z=30, Dij=mode)
        p0 = np.array([yy * 1e5, (1 - yy) * 1e5])
        Th_prev, J = m.steady(p0, Th_start=Th_prev)
        row.append(SP(J, p0))
    SWEEP[mode] = np.array(row)

print("Selectivity against MEASUREMENT. This is the only comparison on the page")
print("that is validation rather than reproduction.\\n")
print(f"{'feed':>7}{'measured':>10}{'finite Dij':>12}{'ratio':>8}"
      f"{'infinite Dij':>14}{'ratio':>8}")
MEAS_ROWS = [("95-5", MEAS["bakker_95_5"], SP_955, SP_955_INF),
             ("50-50", MEAS["bakker_50_50"], SP_5050, SP_5050_INF)]
for lab, meas, fin, inf in MEAS_ROWS:
    print(f"{lab:>7}{meas:10.0f}{fin:12.1f}{fin/meas:8.2f}{inf:14.2f}{inf/meas:8.3f}")
BAKKER_FIN = [f / m_ for _, m_, f, _ in MEAS_ROWS]
BAKKER_INF = [i / m_ for _, m_, _, i in MEAS_ROWS]
print(f"\\nThe model with finite exchange is {min(BAKKER_FIN):.1f}x to "
      f"{max(BAKKER_FIN):.1f}x the measurement;")
print(f"with the exchange removed it is {min(BAKKER_INF):.2f}x to "
      f"{max(BAKKER_INF):.2f}x, i.e. {1/max(BAKKER_INF):.0f}x to "
      f"{1/min(BAKKER_INF):.0f}x too LOW.")
print("So the measurement discriminates decisively between the two models - the")
print("gap between them is two orders of magnitude and the data sit on one side")
print(f"of it - while agreeing quantitatively with neither.")
print(f"\\nThe review's own words for the 95-5 case are 'quite close', and at "
      f"{BAKKER_FIN[0]:.2f}x that is")
print("defensible. It quotes the 50-50 measurement in the next sentence and does")
print(f"not comment on it; there the finite-exchange model is {BAKKER_FIN[1]:.1f}x high.")
print("Neither number is a fit: nothing in this calculation was tuned to either.")'''))

cells.append(code(r'''fig, ax = plt.subplots(figsize=(6.4, 4.4))
p_nC4 = (1 - ys) * 100.0
ax.semilogy(p_nC4, SWEEP["vignes"], "o-", lw=2.0, color="tab:blue",
            label=r"M-S, finite $\eth_{ij}$ (eq. 43)")
ax.semilogy(p_nC4, SWEEP["inf"], "s-", lw=2.0, color="tab:orange",
            label=r"M-S, $\eth_{ij}\to\infty$ (eq. 47)")
ax.semilogy([5.0, 50.0], [MEAS["bakker_95_5"], MEAS["bakker_50_50"]], "k*",
            ms=15, label="Bakker, measured")
ax.axhline(SP_PURE, color="tab:green", lw=1.3, ls="--",
           label=f"pure-component expectation, {SP_PURE}")
ax.set(xlabel="upstream nC4 partial pressure (kPa), total 100 kPa",
       ylabel=r"permeation selectivity $S_P$ (eq. 60)",
       title="nC4/C1 across MFI at 300 K")
ax.legend(fontsize=8); fig.tight_layout(); plt.show()
print("Every curve here is a prediction from pure-component data alone. The two")
print("stars are the only measured points the review states as numbers.")'''))

# ------ section 5: hexanes
cells.append(md(r"""### 5. The hexane isomers, and the review's best experimental case

Section 4.3 is the separation the entropy effect was found for. The review takes
$\eth_{n\mathrm{C6}}=5\,\eth_{3\mathrm{MP}}$ after Cavalcante & Ruthven, prints
$S_P = 43.2$ at 33 kPa upstream partial pressures, and reports Funke et al.'s
measured 24 at 15 kPa upstream — against 1.3 for the same pair measured as pure
components."""))

cells.append(code('''D_RATIO = mprm("nC6", "Section 4.3", "D_ratio_to_3MP")
P_I0 = mprm("mixture", "Section 4.3", "p_i0")
Ps = np.array([1e3, 2e3, 5e3, 15e3, 30e3, 2 * P_I0, 100e3])
Th_prev = None; SP_hex = []
for P in Ps:
    m = Membrane([NC6, MP3], [NC6.Th_sat, MP3.Th_sat], [D_RATIO, 1.0], n_z=30)
    Th_prev, J = m.steady(np.array([P / 2, P / 2]), Th_start=Th_prev)
    SP_hex.append(float(J.mean(0)[0] / J.mean(0)[1]))
SP_hex = np.array(SP_hex)

# the headline point, recomputed COLD so that no continuation chain can carry a
# warm-start artefact into a reported number, and grid-refined
SP_66, grid = None, []
Th_g = None
for nz in (15, 30, 60):
    m = Membrane([NC6, MP3], [NC6.Th_sat, MP3.Th_sat], [D_RATIO, 1.0], n_z=nz)
    Th_g, J = m.steady(np.array([P_I0, P_I0]), Th_start=None if nz == 15 else
                       np.exp(np.stack([np.interp(
                           np.linspace(0, 1, nz + 2)[1:-1],
                           np.linspace(0, 1, Th_g.shape[0] + 2)[1:-1],
                           np.log(Th_g[:, k])) for k in range(2)], -1)))
    grid.append((nz, float(J.mean(0)[0] / J.mean(0)[1]), m.flux_spread))
SP_66 = grid[-1][1]
print("nC6/3MP selectivity at the review's condition, grid-refined:")
for nz, v, spread in grid:
    print(f"   n_z = {nz:4d}   S_P = {v:.5f}   flux spread {spread:.1e}")
print(f"   the review prints {PRINTED['hex_SP_nC6_3MP']}   "
      f"({(SP_66/PRINTED['hex_SP_nC6_3MP']-1)*100:+.3f} %)")
SP_66_COLD = grid[0][1]
print(f"   from a cold start at n_z = {grid[0][0]}, with no continuation at all: "
      f"{SP_66_COLD:.5f}")

P_SP10 = 10 ** brentq(
    lambda lp: (lambda M: (lambda TJ: TJ[1].mean(0)[0] / TJ[1].mean(0)[1])(
        M.steady(np.array([10 ** lp / 2, 10 ** lp / 2]))))(
        Membrane([NC6, MP3], [NC6.Th_sat, MP3.Th_sat], [D_RATIO, 1.0], n_z=30)) - 10.0,
    np.log10(2e3), np.log10(5e3), xtol=4e-3)
print(f"\\nS_P crosses 10 at a total upstream pressure of {P_SP10/1e3:.2f} kPa;")
print(f"the review states 'in excess of {PRINTED['hex_P_SP10']/1e3:.0f} kPa' "
      f"({P_SP10/PRINTED['hex_P_SP10']:.2f}x), which it reads off its own")
print("Fig. 29(c). Scored as a figure-read statement, not as a printed number.")

k15 = int(np.where(Ps == 15e3)[0][0])
SP_15 = float(SP_hex[k15])
FUNKE = MEAS["funke_mixture"]
FUNKE_PURE = MEAS["funke_pure"]
print(f"\\nAGAINST MEASUREMENT. At the {Ps[k15]/1e3:.0f} kPa upstream pressure of "
      f"Funke et al.'s run:")
print(f"   computed S_P (M-S with finite exchange, IAST) : {SP_15:.2f}")
print(f"   measured by Funke et al.                     : {FUNKE:.1f}"
      f"   ({(SP_15/FUNKE-1)*100:+.1f} %)")
print(f"   measured by Funke et al. on PURE components  : {FUNKE_PURE}")
print(f"\\nThe mixture is {FUNKE/FUNKE_PURE:.0f}x more selective than the pure "
      f"components in the")
print(f"measurement and {SP_15/FUNKE_PURE:.0f}x in the model - and the model was "
      "given nothing but")
print("pure-component isotherm fits and one diffusivity ratio. This is the best")
print("experimental agreement in the review's mixture sections and, unlike")
print("Bakker's methane-butane pair, it is quantitative.")'''))

cells.append(md(r"""**Nothing was fitted — but three things were *chosen*, and the $+4.3\ \%$ means
nothing unless a reader can see them.** The review states the measurement in one
sentence, "for operation with an upstream pressure of 15 kPa, Funke et al. have
experimentally determined a value $S_P$ of 24", and states neither the feed
composition nor whether 15 kPa is the total or the per-component pressure. Each
of the three readings below is worth tens of per cent, which is an order of
magnitude more than the agreement being claimed. They are measured rather than
argued."""))

cells.append(code('''def sp_hex(P_total, y_nC6=0.5, ratio=None):
    """S_P of eq. (60) for the nC6/3MP pair: the flux ratio DIVIDED by the
    upstream partial-pressure ratio. At an equimolar feed the two coincide,
    which is why the cell above can use the bare flux ratio."""
    r = D_RATIO if ratio is None else ratio
    m = Membrane([NC6, MP3], [NC6.Th_sat, MP3.Th_sat], [r, 1.0], n_z=30)
    p0 = np.array([y_nC6 * P_total, (1 - y_nC6) * P_total])
    _, J = m.steady(p0)                      # cold, no continuation
    Jm = J.mean(0)
    return float(Jm[0] / Jm[1]) / (p0[0] / p0[1]), float(Jm[0] / Jm[1])


CHOICES = [
    ("AS BUILT: 15 kPa total, equimolar, x5", 15e3, 0.5, 5.0),
    ("  15 kPa EACH instead (30 kPa total)", 30e3, 0.5, 5.0),
    ("  feed y_nC6 = 0.4 instead of 0.5", 15e3, 0.4, 5.0),
    ("  feed y_nC6 = 0.6 instead of 0.5", 15e3, 0.6, 5.0),
    ("  D_nC6 = 4 D_3MP instead of 5", 15e3, 0.5, 4.0),
    ("  D_nC6 = 6 D_3MP instead of 5", 15e3, 0.5, 6.0),
]
print(f"THE THREE CHOICES BEHIND THAT {(SP_15/FUNKE-1)*100:+.1f} %, EACH MEASURED\\n")
print(f"{'reading':>40}{'S_P (eq. 60)':>14}{'vs 24':>9}{'flux ratio':>13}")
SENS = []
for lab, Pt, yy, rr in CHOICES:
    s60, sraw = sp_hex(Pt, y_nC6=yy, ratio=rr)
    SENS.append((lab, s60, sraw))
    print(f"{lab:>40}{s60:14.2f}{(s60/FUNKE-1)*100:8.0f}%{sraw:13.2f}")
FUNKE_CHOICE_SPREAD = max(s for _, s, _ in SENS) / min(s for _, s, _ in SENS)
print(f"\\n   The first row is a COLD solve and reproduces the {SP_15:.2f} reported above,")
print(f"   which came off a pressure sweep: {SENS[0][1]:.5f} against {SP_15:.5f}.")
print("   The last column is the bare flux ratio; it differs from S_P only for the")
print("   two non-equimolar rows, where eq. (60) divides out the feed ratio. The")
print(f"   measured {FUNKE:g} is an S_P, so the S_P column is the comparable one - but a")
print("   reader who thinks Funke reported a flux ratio should use the last.")
print(f"\\n   spread across the alternatives: {FUNKE_CHOICE_SPREAD:.1f}x, against the "
      f"{abs(SP_15/FUNKE-1)*100:.1f} % agreement claimed.")
print("\\nWHY THE READINGS USED ARE THE BETTER-SUPPORTED ONES, and they are still")
print("readings. The review gives the 43.2 point as 'keeping the upstream")
print("compartments at a TOTAL pressure of 66 kPa' while the Fig. 29(a) caption")
print("gives the same condition as p_i0 = 33 kPa, so in this paper an upstream")
print("pressure quoted in running text is a total; and the S_P > 10 threshold is")
print("stated as 'the pressure in the upstream compartment increases beyond 2 kPa'")
print(f"against the {P_SP10/1e3:.2f} kPa TOTAL computed above. Equimolar is the "
      "composition of every")
print("simulation in Section 4.3 and the curve the review places Funke's square")
print("on. Neither is stated for the MEASUREMENT, and the review does not say")
print(f"what Funke's feed was. The last two rows are not readings at all - the")
print(f"ratio of {D_RATIO:g} IS printed - but it is a single experimental number from a")
print("third paper (Cavalcante & Ruthven, recorded under origin_not_consulted and")
print("not opened here), and it drives the whole hexane result.")
print(f"\\nSo: nothing was tuned, and the {(SP_15/FUNKE-1)*100:+.1f} % should be read as one "
      "point inside")
print(f"a {FUNKE_CHOICE_SPREAD:.1f}x band of defensible readings, against a "
      "measurement printed to two")
print("significant figures with no stated uncertainty.")'''))

cells.append(code('''print("The other two hexane cases the review quotes, both read off ITS OWN")
print("figures rather than printed as computed results:\\n")
D32 = mprm("3MP", "Section 4.3", "D_ratio_to_22DMB")
D_TERN = mprm("ternary", "Section 4.3", "D_ratio_nC6_to_22DMB")
m = Membrane([MP3, DMB], [MP3.Th_sat, DMB.Th_sat], [D32, 1.0], n_z=30)
_, J = m.steady(np.array([30e3, 30e3]))
SP_BIN = float(J.mean(0)[0] / J.mean(0)[1])
Th_prev = None
for P in (2e3, 1e4, 9e4):
    m3 = Membrane([NC6, MP3, DMB], [NC6.Th_sat, MP3.Th_sat, DMB.Th_sat],
                  [D_TERN, D_TERN, 1.0], n_z=30)
    Th_prev, J3 = m3.steady(np.full(3, P / 3), Th_start=Th_prev)
SP_TERN = float(J3.mean(0)[1] / J3.mean(0)[2])
print(f"   3MP/2,2DMB, binary, 30 kPa each  : computed {SP_BIN:6.2f}   "
      f"the review says {PRINTED['hex_SP_3MP_DMB_binary']:.0f}"
      f"   ({(SP_BIN/PRINTED['hex_SP_3MP_DMB_binary']-1)*100:+.0f} %)")
print(f"   the same pair with nC6 present   : computed {SP_TERN:6.2f}   "
      f"the review says 'in excess of {PRINTED['hex_SP_3MP_DMB_ternary']:.0f}'"
      f"   (satisfied: {SP_TERN > PRINTED['hex_SP_3MP_DMB_ternary']})")
print(f"\\nAdding nC6 raises the 3MP/2,2DMB selectivity by "
      f"{SP_TERN/SP_BIN:.1f}x here, against the")
print(f"review's own factor of at least "
      f"{PRINTED['hex_SP_3MP_DMB_ternary']/PRINTED['hex_SP_3MP_DMB_binary']:.1f}x. "
      "The mechanism is reproduced; the two numbers")
print("are figure readings by the authors on a log axis and are not scored with")
print("the printed results.")'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2))
axes[0].loglog(Ph / 1e3, Thh[:, 0], lw=2.0, color="tab:blue", label=r"nC6")
axes[0].loglog(Ph / 1e3, Thh[:, 1], lw=2.0, color="tab:red", label=r"3MP")
axes[0].loglog(Ph / 1e3, Thh.sum(1), lw=1.3, ls="--", color="0.4",
               label=r"$\Theta_{\rm mix}$")
axes[0].axhline(4.0, color="0.4", lw=1.0, ls=":")
axes[0].axvline(P_TH4 / 1e3, color="tab:green", lw=1.0)
axes[0].annotate(f"$\\Theta_{{\\rm mix}}=4$ at {P_TH4/1e3:.2f} kPa",
                 (P_TH4 / 1e3 * 1.3, 0.02), fontsize=8, color="tab:green")
axes[0].set(xlabel="total upstream pressure (kPa)",
            ylabel=r"$\Theta$ (molecules per unit cell)", ylim=(1e-2, 20),
            title="50-50 nC6/3MP sorption at 362 K (IAST)")
axes[0].legend(fontsize=8, loc="upper left")

axes[1].loglog(Ps / 1e3, SP_hex, "o-", lw=2.0, color="tab:blue",
               label="M-S + IAST, this page")
axes[1].plot([15.0], [FUNKE], "k*", ms=15, label="Funke et al., mixture")
axes[1].plot([15.0], [FUNKE_PURE], "ko", ms=9, mfc="none", mew=1.8,
             label="Funke et al., pure components")
axes[1].plot([2 * P_I0 / 1e3], [PRINTED["hex_SP_nC6_3MP"]], "P", ms=11,
             color="tab:red", label=f"the review's printed "
             f"{PRINTED['hex_SP_nC6_3MP']}")
axes[1].axhline(10.0, color="0.5", lw=1.0, ls=":")
axes[1].set(xlabel="total upstream pressure (kPa)",
            ylabel=r"$S_P$ (nC6 over 3MP)",
            title="nC6/3MP permeation selectivity at 362 K")
axes[1].legend(fontsize=8, loc="upper left")
fig.tight_layout(); plt.show()'''))

# ------ section 6: the counterfactuals
cells.append(md(r"""### 6. Two counterfactuals the review states, and what they actually give

**The multicomponent Langmuir.** The review says that with eq. (39) and the
saturation capacities of C1 and nC4 set equal, $S_P$ "would be predicted to be
independent of the upstream composition and have a constant value of 800". That
is two claims again — a structural one and a numerical one — and they come apart.

**The zeolite 4A counter-example.** Section 4.4 is the case that goes the other
way, and a page showing only finite $\eth_{ij}$ winning would misreport the
review's own conclusion."""))

cells.append(code('''Th_eq = 10.0        # "the saturation capacities set equal to one another"
LANG = []
for yy in (0.95, 0.80, 0.50, 0.20, 0.05):
    m = Membrane([DSL(C1.bA, Th_eq), DSL(NC4.bA, Th_eq)], [Th_eq, Th_eq],
                 [D_C1, D_NC4], n_z=30, iso_model="langmuir")
    p0 = np.array([yy * 1e5, (1 - yy) * 1e5])
    _, J = m.steady(p0)
    LANG.append((yy, SP(J, p0)))
LANG_SPREAD = max(v for _, v in LANG) / min(v for _, v in LANG)
LANG_MEAN = float(np.mean([v for _, v in LANG]))
print("Multicomponent Langmuir, eq. (39), equal saturation capacities:")
print(f"{'upstream y_C1':>15}{'S_P':>12}")
for yy, v in LANG:
    print(f"{yy:15.2f}{v:12.2f}")
IAST_SPREAD = float(SWEEP["vignes"][1:].max() / SWEEP["vignes"][1:].min())
print(f"\\n   spread from 95 % down to 5 % methane in the feed: {LANG_SPREAD:.3f}x")
print(f"   IAST over exactly the same five feeds:            {IAST_SPREAD:.1f}x")
print("   -> THE STRUCTURAL CLAIM IS REPRODUCED. The multicomponent Langmuir does")
print("      make S_P composition-independent, to better than half a per cent,")
print("      and IAST does not. That is exactly the review's argument for IAST.")
print(f"\\n   the VALUE, however, is {LANG_MEAN:.0f} here against the printed "
      f"{PRINTED['mix_SP_langmuir']:.0f} "
      f"({LANG_MEAN/PRINTED['mix_SP_langmuir']:.2f}x).")
print("   The review does not say which saturation capacity it set the two equal")
print("   TO, nor which affinity it used once the dual-site fit is collapsed to a")
print("   single site. The first turns out not to matter - with equal capacities")
print("   the occupancies are independent of their common value and Theta_sat")
print("   cancels out of the flux ratio, checked below - but the second does. The")
print("   Table 1 site-A affinities are used here and the counterfactual is")
print("   reported as UNREPRODUCED rather than tuned until it lands on 800.")
LANG_CAP_CHECK = []
for cap in (4.0, 19.0):
    m = Membrane([DSL(C1.bA, cap), DSL(NC4.bA, cap)], [cap, cap],
                 [D_C1, D_NC4], n_z=30, iso_model="langmuir")
    p0 = np.array([50e3, 50e3]); _, J = m.steady(p0)
    LANG_CAP_CHECK.append(SP(J, p0))
LANG_CAP_DEV = abs(LANG_CAP_CHECK[0] / LANG_CAP_CHECK[1] - 1)
print(f"\\n   (Theta_sat = {4.0:.0f} vs {19.0:.0f} at 50-50: "
      f"{LANG_CAP_CHECK[0]:.3f} vs {LANG_CAP_CHECK[1]:.3f}, {LANG_CAP_DEV:.1e} apart)")'''))

cells.append(md(r"""**The review's conclusion 9 goes the other way for zeolite 4A**, and it is
quoted here rather than paraphrased into agreement. Section 4.4 compares both
implementations against Habgood's measured N₂/CH₄ uptake in 4A and finds that
"the M–S model, assuming $\eth_{ij}\to\infty$, does a very good job", concluding
that diffusion of N₂ and CH₄ in 4A "is essentially free from vacancy correlation
effects". The physical reason the review gives is topological: in 4A the
intracage hopping is fast and is not the limiting step, whereas in MFI the
exchange happens at the channel intersections and is not fast.

So the finding is not "finite $\eth_{ij}$ is right". It is that **$\eth_{ij}$ is
a property of the host topology, and MFI and 4A land at opposite limits of it**.
That case cannot be reproduced here: its only validation is Figs. 32–34, and the
review prints no numerical result for it. It is scoped out for exactly the reason
the rest of Sections 4.4–4.7 are.

### 7. Self-diffusivity, where the correlation effect is largest

Eqs. (53)–(56). The self-diffusivity contains no $[\Gamma]$ at all, so it
isolates the exchange coefficients from the sorption thermodynamics — which is
why the review uses it to argue that correlation effects are real. There is no
printed number to check against; the comparison is Fig. 21."""))

cells.append(code('''TS = np.array([mprm("CH4", "Section 3.4", "Theta_sat"),
               mprm("CF4", "Section 3.4", "Theta_sat")])
D0S = np.array([mprm("CH4", "Section 3.4", "D_MS_0"),
                mprm("CF4", "Section 3.4", "D_MS_0")])
TTOT = mprm("mixture", "Section 3.4", "Theta_total")


def D_self(Th, D0, Th_sat, facile=False):
    """eq. (56). D_11 = D_1 by eq. (43) applied at i = j - an algebraic
    identity, not a check, and it is used here only because the review does."""
    th = Th / Th_sat
    D_i = D0 * (1 - th.sum())                     # eq. (19) on the total occupancy
    out = []
    for i in (0, 1):
        j = 1 - i
        if facile:
            out.append(D_i[i])
            continue
        xi = th[i] / (th[i] + th[j])
        Dij = D_i[i] ** xi * D_i[j] ** (1 - xi)
        out.append(1.0 / (1 / D_i[i] + th[i] / D_i[i] + th[j] / Dij))
    return np.array(out)


print(f"Snurr & Karger's CH4/CF4 in MFI at 200 K, total loading {TTOT:.0f} "
      "molecules per unit cell")
print(f"(Theta_sat = {TS[0]:.0f} and {TS[1]:.0f}; "
      f"D(0) = {D0S[0]:.1e} and {D0S[1]:.1e} m2/s)\\n")
print(f"{'Theta_CH4':>10}{'D*_CH4':>12}{'facile':>12}{'ratio':>8}"
      f"{'D*_CF4':>12}{'facile':>12}{'ratio':>8}")
RATIOS = []
for T1 in (1.0, 3.0, 6.0, 9.0, 11.0):
    T = np.array([T1, TTOT - T1])
    a = D_self(T, D0S, TS); b = D_self(T, D0S, TS, facile=True)
    RATIOS.append(b[0] / a[0])
    print(f"{T1:10.1f}{a[0]:12.3e}{b[0]:12.3e}{b[0]/a[0]:8.2f}"
          f"{a[1]:12.3e}{b[1]:12.3e}{b[1]/a[1]:8.2f}")
SELF_RATIO = (float(min(RATIOS)), float(max(RATIOS)))
print(f"\\nRemoving the exchange raises the methane self-diffusivity by "
      f"{SELF_RATIO[0]:.1f}x to {SELF_RATIO[1]:.1f}x.")
print("The review states that the dotted facile-exchange lines of its Fig. 21(b)")
print("lie well above the MD points while the finite-exchange lines of Fig. 21(a)")
print("agree - but it prints no number for either, so THIS COMPARISON HAS NO")
print("TARGET. What is shown is the size and sign of the effect the review")
print("attributes to correlation, not evidence that the review is right about it.")
print("\\nONE CHOICE IS MADE SILENTLY BY EQ. (19) AND IS STATED HERE INSTEAD.")
print("Eq. (19) as printed is D_i = D_i(0)(1 - Theta_i/Theta_i,sat), the species'")
print("OWN occupancy. The cell above uses D0 * (1 - sum_k theta_k), the TOTAL")
print("occupancy, which is a generalisation to a mixture that the review states")
print("nowhere; Section 3.4 says only 'Eqs. (19), (43) and (56)'. The membrane")
print("results measure all three readings of the same ambiguity in Validation")
print("below; this section hard-codes one, and nothing on the page depends on it")
print("because there is no printed target here. What it costs, at the loadings")
print("above, is the difference between:")
for T1 in (1.0, 6.0, 11.0):
    T = np.array([T1, TTOT - T1]); th = T / TS
    print(f"      Theta_CH4 = {T1:4.1f}:  1 - theta_total = {1-th.sum():.4f}   "
          f"1 - theta_CH4 = {1-th[0]:.4f}   ({(1-th[0])/(1-th.sum()):.2f}x)")
print("\\nTwo relations used above are ALGEBRAIC IDENTITIES and are not checks:")
print("  eq. (43) at i = j returns D_ii = D_i, which is what x^a y^(1-a) does")
print("  when y = x; and eq. (55), D*_1 = D_1/(1+theta_1), is eq. (54) with")
print("  D_11 = D_1 substituted and the other terms dropped. Both are stated")
print("  because the review uses them; neither is scored.")'''))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

**This page has both kinds of comparison and keeps them apart in the words used.**
Four numbers are *measurements* — Bakker's two and Funke's two — and comparing
against them is **validation**. Everything else is a number Krishna & Baur
computed, and comparing against those is **reproduction**.

Six checks with something to fail against, then the break table."""))

cells.append(code('''print("1. The IAST solve, three ways it can be wrong\\n")

# (a) round trip: forward then inverse must return the pressures it started from
p_rt = np.array([[95e3, 5e3], [1e3, 1e3], [1e-2, 1e-4], [5e5, 5e5]])
Th_rt = mix.forward(p_rt)
IAST_ROUNDTRIP = float(np.max(np.abs(mix.inverse(Th_rt) / p_rt - 1)))
print(f"   (a) forward (p)->(Theta) then inverse (Theta)->(p), 4 states: "
      f"{IAST_ROUNDTRIP:.2e}")
print("       tests the inverse against the forward. It CANNOT see an error")
print("       common to both - the surface potential, which is (b).")

# (b) the surface potential, analytic against numerical quadrature
worst = 0.0
for it in (C1, NC4, NC6, MP3):
    for P in (1e-2, 1e2, 1e5, 1e7):
        num = quad(lambda x: it.theta(x) / x, 1e-14, P, limit=400)[0]
        worst = max(worst, abs(it.psi(P) / num - 1))
PSI_QUAD = float(worst)
print(f"\\n   (b) eq. (64) in closed form against adaptive quadrature of "
      f"Theta^0(P)/P,\\n       4 isotherms x 4 pressures: {PSI_QUAD:.2e}")
print("       This one CAN fail: the quadrature never sees the log form.")

# (c) IAST must reduce to the multicomponent Langmuir when the saturation
#     capacities are equal - the review's own explanation of why eq. (39)
#     predicts a constant selectivity
eq_a, eq_b = DSL(NC4.bA, 6.0), DSL(C1.bA, 6.0)
lang_mix = IAST([eq_a, eq_b])
worst_eq = worst_ne = 0.0
for Ptot in np.logspace(-2, 7, 19):
    for yy in (0.05, 0.5, 0.95):
        p = np.array([[yy * Ptot, (1 - yy) * Ptot]])
        T = lang_mix.forward(p)[0]
        den = 1 + eq_a.bA * p[0, 0] + eq_b.bA * p[0, 1]
        L = np.array([6.0 * eq_a.bA * p[0, 0] / den, 6.0 * eq_b.bA * p[0, 1] / den])
        worst_eq = max(worst_eq, float(np.max(np.abs(T / L - 1))))
        T2 = mix.forward(np.array([[p[0, 1], p[0, 0]]]))[0]   # the real, unequal pair
        den2 = 1 + C1.bA * p[0, 1] + NC4.bA * p[0, 0]
        L2 = np.array([C1.Th_sat * C1.bA * p[0, 1] / den2,
                       NC4.Th_sat * NC4.bA * p[0, 0] / den2])
        worst_ne = max(worst_ne, float(np.max(np.abs(T2 / L2 - 1))))
IAST_LANGMUIR = float(worst_eq)
IAST_UNEQUAL = float(worst_ne)
print(f"\\n   (c) IAST vs the multicomponent Langmuir with EQUAL saturation")
print(f"       capacities, 19 pressures x 3 compositions : {IAST_LANGMUIR:.2e}")
print(f"       the same comparison with the REAL Table 1 capacities (19 and 10):"
      f" {IAST_UNEQUAL*100:.0f} %")
print("       The first is the theorem the review relies on and it is a real test")
print("       of the IAST solve - a wrong surface potential breaks it. The second")
print("       says the theorem does not apply to this system, which is the whole")
print("       reason IAST is needed.")'''))

cells.append(code('''print("2. [Gamma] and [B], each by two routes that share no code\\n")
states = mix.forward(np.logspace(1, 6.5, 12)[:, None] * y)
G_an = mix.inverse_and_gamma(states, np.array([C1.Th_sat, NC4.Th_sat]))[1]
G_nu = gamma_numerical(mix, states, np.array([C1.Th_sat, NC4.Th_sat]))
GAMMA_TWO_ROUTES = float(np.max(np.abs(G_an / G_nu - 1)))
print(f"   [Gamma] by implicit differentiation of eqs. (62)-(65) against central")
print(f"   differences of the IAST inverse map, 12 states x 4 elements: "
      f"{GAMMA_TWO_ROUTES:.2e}")
print("   (the finite-difference floor; the break table moves it to O(1))")
print("   WHAT IT CANNOT SEE: both routes call the same IAST solve, so an error")
print("   inside eqs. (62)-(65) moves them together. That is what check 1 is for.")

rng = np.random.default_rng(20260803)
worst = worst_broken = 0.0
for _ in range(200):
    n = int(rng.integers(2, 4))
    th = rng.uniform(0.01, 0.30, n)
    D = 10.0 ** rng.uniform(-11, -9, n)
    d = rng.normal(size=n)
    u1 = friction_solve(th, D, d)
    B = B_matrix(th[None, :], D[None, :])[0]
    u2 = np.linalg.solve(B, d)
    worst = max(worst, float(np.max(np.abs(u1 / u2 - 1))))
    Bb = B.copy()
    for i in range(n):
        Bb[i, i] = 1.0 / D[i]                     # the exchange sum dropped
    worst_broken = max(worst_broken, float(np.max(np.abs(u1 / np.linalg.solve(Bb, d) - 1))))
FRICTION_TWO_ROUTES = worst
FRICTION_BROKEN = worst_broken
print(f"\\n   eq. (44)'s [B], assembled and inverted, against the raw friction")
print(f"   system of eq. (42) solved matrix-free by GMRES - 200 random 2- and")
print(f"   3-component states: {FRICTION_TWO_ROUTES:.2e}")
print(f"   with the exchange sum dropped from B_ii: {FRICTION_BROKEN*100:.0f} %")
print("   The GMRES route never forms [B], so it cannot inherit an error in the")
print("   B_ii / B_ij formulas. This is A4.2's replacement for a check that")
print("   compared [B] against itself, reused.")

# the n = 1 collapse: eq. (46) must give eq. (14), D = eth Gamma
th1 = np.array([[0.37]])
B1 = B_matrix(th1, np.array([[D_C1]]))[0]
G1 = np.array([[2.4]])
COLLAPSE_46 = float(abs((np.linalg.solve(B1, G1))[0, 0] / (D_C1 * G1[0, 0]) - 1))
print(f"\\n   eq. (46) at n = 1 against eq. (14), D = eth Gamma: {COLLAPSE_46:.1e}")
print("   ALGEBRAIC IDENTITY, labelled as one: at n = 1 the sum in eq. (44) is")
print("   empty so [B] = 1/eth by construction. It confirms the indexing of the")
print(f"   code and nothing else. The collapse that CAN fail is the one in")
print(f"   Results, where the n = 1 SOLVER is run against the closed-form")
print(f"   eq. (12) integral: {COLLAPSE_1:.1e}.")'''))

cells.append(code('''print("3. Discretisation, and the two modelling choices the review leaves open\\n")
GRID = []
Th_prev = None
for nz in (20, 40, 80):
    m = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=nz)
    _, J = m.steady(FEEDS["95-5"])
    GRID.append((nz, SP(J, FEEDS["95-5"]), m.flux_spread))
print(f"   {'n_z':>6}{'S_P':>12}{'rel. to n_z=80':>17}{'flux spread':>14}")
for nz, v, spr in GRID:
    print(f"   {nz:6d}{v:12.4f}{v/GRID[-1][1]-1:17.2e}{spr:14.1e}")
GRID_ERR = abs(GRID[1][1] / GRID[-1][1] - 1)
ORDER = float(np.log2(abs(GRID[0][1] - GRID[-1][1]) / abs(GRID[1][1] - GRID[-1][1])))
print(f"   observed order {ORDER:.2f}; the steady production runs use n_z = 40, so")
print(f"   they carry {GRID_ERR*100:.3f} % of spatial error - far below every steady")
print("   deviation reported above. The flux spread is div N at steady state and")
print("   is the convergence assertion; it is never inferred from the iteration")
print("   count.")
print("\\n   WHAT THIS STUDY DOES NOT COVER, stated because an earlier version of")
print("   this page let it imply otherwise: the TRANSIENT. The steady loading")
print("   profile is smooth, so it converges quickly and cheaply here; the")
print("   transient carries a self-sharpening n-butane front and needs a much")
print(f"   finer grid. That is refined separately, in Results, at n_z = {NZ_PROD}, and")
print("   the numbers move by an order of magnitude more than anything in this")
print("   table. A grid study is only evidence for the solution it was run on.")

print("\\n   Path independence. The 50-50 feeds need the continuation fallback; the")
print("   95-5 feed converges from a cold start, so both routes can be run on it")
print("   and compared. If they disagreed, every continued result would be suspect.")
m_a = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=40)
_, J_a = m_a.steady(FEEDS["95-5"])
m_b = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=40)
_, J_b = m_b.steady(FEEDS["95-5"], spread_tol=0.0)       # forces the continuation
PATH_DEV = abs(SP(J_b, FEEDS["95-5"]) / SP(J_a, FEEDS["95-5"]) - 1)
print(f"      cold start (continued={m_a.continued}) : S_P = {SP(J_a, FEEDS['95-5']):.6f}")
print(f"      forced continuation (continued={m_b.continued}) : "
      f"S_P = {SP(J_b, FEEDS['95-5']):.6f}   ({PATH_DEV:.1e})")
print("   Reported values are start-independent, which is the condition a")
print("   continuation chain has to meet before any number from it is quoted.")

print("\\n   The downstream pressure ('vanishing values by means of a sweep gas'):")
PD = []
for f_ in (1e-4, 1e-5, 1e-6, 1e-7):
    m = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=30,
                 pd_frac=f_)
    _, J = m.steady(FEEDS["95-5"])
    PD.append((f_, SP(J, FEEDS["95-5"])))
    print(f"      p_delta/p_0 = {f_:.0e}   S_P = {PD[-1][1]:.3f}   "
          f"(div N = {m.flux_spread:.0e})")
PD_SENS = abs(PD[0][1] / PD[-1][1] - 1)
print(f"   three decades of sweep pressure move S_P by {PD_SENS*100:.4f} %.")
print("   Well below that the Newton stops converging, which the div N column")
print("   would show; the production value of 1e-6 is clear of it.")

print("\\n   The loading dependence of eth, which Section 3.2 leaves as 'e.g.")
print("   eq. (17) or eq. (19)' and never fixes for the membrane:")
ALT = {}
for lab, dep in (("eq. (17), constant", None),
                 ("eq. (19), eth(0)(1-theta_i)", lambda th: 1.0 - th),
                 ("eq. (19) on the TOTAL occupancy",
                  lambda th: np.repeat((1.0 - th.sum(-1))[:, None], th.shape[-1], -1))):
    m = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=30,
                 loading_dep=dep)
    _, J = m.steady(FEEDS["95-5"])
    ALT[lab] = SP(J, FEEDS["95-5"])
    print(f"      {lab:35s} S_P = {ALT[lab]:8.2f}   "
          f"({ALT[lab]/PRINTED['mix_SP_955_finite']-1:+.1%} vs printed)")
DEP_SPREAD = max(ALT.values()) / min(ALT.values())
print(f"   {DEP_SPREAD:.2f}x between the readings. eq. (17) is the one that")
print("   reproduces the printed selectivity and the one A4.7 established for the")
print("   single components, so it is the production choice - but the review does")
print("   not say so, and the other readings are not far off either.")'''))

cells.append(code('''print("4. Summary of every printed target\\n")
SUMMARY = [
    ("S_P, 95-5, finite Dij (eq. 60)", PRINTED["mix_SP_955_finite"], SP_955),
    ("S_P, 95-5, infinite Dij (eq. 47)", PRINTED["mix_SP_955_inf"], SP_955_INF),
    ("S_P, 50-50, finite Dij", PRINTED["mix_SP_5050_finite"], SP_5050),
    ("S_P, 50-50, infinite Dij", PRINTED["mix_SP_5050_inf"], SP_5050_INF),
    ("S_P, nC6/3MP at 66 kPa, 362 K", PRINTED["hex_SP_nC6_3MP"], SP_66),
] + [(f"N ({lab})", pr, v) for lab, pr, v in FLUX_ROWS] + [
    ("peak N_C1, 95-5, finite Dij", PRINTED["mix_peak_C1_finite"], PEAKS["vignes"][1]),
    ("peak N_C1, 95-5, infinite Dij", PRINTED["mix_peak_C1_inf"], PEAKS["inf"][1]),
    ("t at that peak (finite Dij)", PRINTED["mix_peak_t_finite"], PEAKS["vignes"][0]),
    ("S (Henry limit), eq. (61)", PRINTED["sorp_S_henry"], S_HENRY),
    ("S_P, multicomponent Langmuir", PRINTED["mix_SP_langmuir"], LANG_MEAN),
]
print(f"{'target':>40}{'printed':>12}{'computed':>12}{'dev':>10}")
for name, pr, comp in SUMMARY:
    print(f"{name:>40}{pr:12.4g}{comp:12.4g}{(comp/pr-1)*100:9.2f}%")
CORE = SUMMARY[:9]                       # the five selectivities and four fluxes
WORST_CORE = max(abs(c / p - 1) for _, p, c in CORE)
PEAK_WORST = max(abs(c / p - 1) for _, p, c in SUMMARY[9:11])
print(f"\\nWorst deviation over the {len(CORE)} STEADY-STATE targets - the five printed")
print(f"selectivities and the four printed fluxes: {WORST_CORE*100:.2f} %, i.e. inside the")
print("printed precision throughout. Those are the results the model exists to")
print("produce and nothing in the calculation was tuned to them.")
print(f"\\nThe two TRANSIENT peak heights sit further out, at "
      f"{(SUMMARY[9][2]/SUMMARY[9][1]-1)*100:+.1f} % and "
      f"{(SUMMARY[10][2]/SUMMARY[10][1]-1)*100:+.1f} %,")
print(f"and both are HIGH. At n_z = 30, where an earlier version of this page")
print(f"stopped, the same two rows read "
      f"{(dict((n,y) for n,_,y in GRID_T['vignes'])[30]/PRINTED['mix_peak_C1_finite']-1)*100:+.1f} % and "
      f"{(dict((n,y) for n,_,y in GRID_T['inf'])[30]/PRINTED['mix_peak_C1_inf']-1)*100:+.1f} %, i.e. LOW, and")
print("that version attributed the sign to backward Euler damping a maximum.")
print(f"Refined in SPACE both cross over; refining the STEP as well moves them")
print(f"further up ({DT_PEAK_SHIFT*100:+.2f} % from one halving). The attribution was wrong "
      "and is")
print("withdrawn. They are reported as a separate row family rather than folded")
print("into the steady score.")
print("\\nThree further targets are reported separately, each for a stated reason:")
print(f"  * the peak TIME, {(PEAKS['vignes'][0]/PRINTED['mix_peak_t_finite']-1)*100:+.1f} %,"
      " AND IT IS THE ONE REAL DISAGREEMENT ON THIS")
print("    PAGE. It is the only quantity here that carries the reconstructed")
print(f"    membrane thickness at all, and it carries it as delta^2: closing the")
print(f"    gap needs delta {DELTA_FOR_PEAK*100:+.1f} %. Grid ({T_GRID_ERR*100:.2f} %) and step "
      f"({DT_TIME_SHIFT*100:.2f} %) refinement")
print("    both move it the wrong way, so it is not a resolution artefact. Left")
print("    open, and the candidates are listed where it is computed.")
print(f"  * the Henry-limit sorption selectivity, {(S_HENRY/PRINTED['sorp_S_henry']-1)*100:+.1f} %."
      " Established above to")
print("    follow from Table 1 by two independent routes and not to be")
print("    explicable by any single-parameter slip in the rows A4.7 certifies.")
print(f"  * the multicomponent-Langmuir counterfactual, {LANG_MEAN/PRINTED['mix_SP_langmuir']:.1f}x."
      " Its structural claim")
print("    is reproduced exactly; its value is underdetermined by what the review")
print("    prints, and nothing was tuned to reach it.")
print("\\nAnd two statements the review reads off its own figures, scored loosely:")
print(f"  * Theta_mix > 4 above {PRINTED['hex_P_theta4']/1e3:.0f} kPa: true, but the "
      f"actual threshold is {P_TH4/1e3:.2f} kPa")
print(f"  * S_P > 10 above {PRINTED['hex_P_SP10']/1e3:.0f} kPa: the crossing is at "
      f"{P_SP10/1e3:.2f} kPa")
print(f"  * 3MP/2,2DMB about {PRINTED['hex_SP_3MP_DMB_binary']:.0f} binary, above "
      f"{PRINTED['hex_SP_3MP_DMB_ternary']:.0f} ternary: computed {SP_BIN:.1f} and "
      f"{SP_TERN:.0f}")

print("\\n5. VALIDATION against measurement - the four tier-2 numbers\\n")
print(f"{'measurement':>28}{'measured':>10}{'this page':>12}{'model/measured':>17}")
VAL_ROWS = [("Bakker, 95-5 C1/nC4", MEAS["bakker_95_5"], SP_955),
            ("Bakker, 50-50 C1/nC4", MEAS["bakker_50_50"], SP_5050),
            ("Funke, nC6/3MP at 15 kPa", MEAS["funke_mixture"], SP_15)]
for lab, meas, comp in VAL_ROWS:
    print(f"{lab:>28}{meas:10.0f}{comp:12.2f}{comp/meas:17.2f}")
print(f"{'Funke, pure components':>28}{MEAS['funke_pure']:10.1f}"
      f"{SP_PURE:12.2f}{SP_PURE/MEAS['funke_pure']:17.2f}")
FUNKE_DEV = abs(SP_15 / MEAS["funke_mixture"] - 1)
BAKKER_WORST = max(abs(c / m_ - 1) for _, m_, c in VAL_ROWS[:2])
print(f"\\nThe last row is a DIFFERENT system's pure-component selectivity"
      f" ({SP_PURE}, C1/nC4)")
print(f"set beside Funke's ({MEAS['funke_pure']}, nC6/3MP): both are near unity and both are")
print("what the mixture calculation has to beat, but they are not the same")
print("quantity and are not scored against each other.")
print(f"\\nBest agreement with measurement: {FUNKE_DEV*100:.0f} % (Funke). Worst: "
      f"{BAKKER_WORST*100:.0f} % (Bakker, 50-50).")'''))

cells.append(md(r"""### Break table

Every number above is made to move by injecting a defect it should catch. Rows
that do not move are kept and labelled, because a check that cannot fail is not
evidence and deleting it hides that fact."""))

cells.append(code('''breaks = []


def sp_with(**kw):
    p0 = kw.pop("p0", FEEDS["95-5"])
    isos = kw.pop("isos", [C1, NC4])
    Ts = kw.pop("Th_sat", [C1.Th_sat, NC4.Th_sat])
    D = kw.pop("D0", [D_C1, D_NC4])
    m = Membrane(isos, Ts, D, n_z=30, **kw)
    _, J = m.steady(p0)
    return SP(J, p0)


BASE_SP = sp_with()

# 1. the mixture isotherm: IAST replaced by the multicomponent Langmuir, with
#    the REAL (unequal) saturation capacities. This is the review's central
#    claim and it must move the answer a long way.
lang_real = sp_with(iso_model="langmuir")
breaks.append(("IAST replaced by the multicomponent Langmuir, eq. (39)",
               "S_P, 95-5", f"{BASE_SP:.2f}", f"{lang_real:.2f}"))

# 2. the exchange coefficients removed
breaks.append(("exchange coefficients removed (eq. 47, Dij -> inf)", "S_P, 95-5",
               f"{BASE_SP:.2f}", f"{sp_with(Dij='inf'):.2f}"))

# 3. eq. (43) with the two exponents swapped - a plausible transcription slip
#    that leaves Dii = Di intact, so the identity in Results cannot see it
class BSwapped:
    pass


def B_swapped(th, D_i, mode="vignes"):
    m, n = th.shape
    B = np.zeros((m, n, n))
    for i in range(n):
        B[:, i, i] = 1.0 / D_i[:, i]
        for j in range(n):
            if j == i:
                continue
            xi = th[:, j] / np.maximum(th[:, i] + th[:, j], 1e-300)   # i and j swapped
            Dij = D_i[:, i] ** xi * D_i[:, j] ** (1 - xi)
            B[:, i, i] += th[:, j] / Dij
            B[:, i, j] = -th[:, i] / Dij
    return B


_B = B_matrix
try:
    globals()["B_matrix"] = B_swapped
    sp_swapped = sp_with()
finally:
    globals()["B_matrix"] = _B
breaks.append(("eq. (43): the two Vignes exponents swapped", "S_P, 95-5",
               f"{BASE_SP:.2f}", f"{sp_swapped:.2f}"))

# 4. [Gamma] replaced by the identity - i.e. the driving force taken as the
#    loading gradient instead of the chemical potential gradient
class NoGamma(Membrane):
    def gamma(s, Th):
        return np.repeat(np.eye(s.n)[None], np.atleast_2d(Th).shape[0], 0)


m_ng = NoGamma([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=30)
_, J_ng = m_ng.steady(FEEDS["95-5"])
breaks.append(("[Gamma] replaced by the identity matrix", "S_P, 95-5",
               f"{BASE_SP:.2f}", f"{SP(J_ng, FEEDS['95-5']):.2f}"))

# 5. [Gamma] kept but its off-diagonal elements zeroed: the coupling the review
#    says makes multicomponent diffusion in zeolites "a strongly coupled process"
class DiagGamma(Membrane):
    def gamma(s, Th):
        G = Membrane.gamma(s, Th)
        return G * np.eye(s.n)[None]


m_dg = DiagGamma([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=30)
_, J_dg = m_dg.steady(FEEDS["95-5"])
breaks.append(("[Gamma] off-diagonal elements zeroed", "S_P, 95-5",
               f"{BASE_SP:.2f}", f"{SP(J_dg, FEEDS['95-5']):.2f}"))

# 6. eq. (60) read as printed (multiplying by the pressure ratio) instead of as
#    the review uses it
breaks.append(("eq. (60) read as printed: x (p20/p10), not / ", "S_P, 95-5",
               f"{BASE_SP:.2f}",
               f"{BASE_SP*(P_NC4/P_C1)**2:.4f}"))

# 7. the surface potential: the second DSL site dropped from eq. (64)
class Psi1(DSL):
    def psi(s, P):
        return s.tA * np.log1p(s.bA * P)


bad_C1 = Psi1(C1.bA, C1.tA, C1.bB, C1.tB)
bad_NC4 = Psi1(NC4.bA, NC4.tA, NC4.bB, NC4.tB)
w = 0.0
for it in (bad_C1, bad_NC4):
    for P in (1e2, 1e5, 1e7):
        w = max(w, abs(it.psi(P) / quad(lambda x: it.theta(x) / x, 1e-14, P,
                                        limit=400)[0] - 1))
breaks.append(("eq. (64): the second DSL site dropped from the potential",
               "psi vs quadrature", f"{PSI_QUAD:.1e}", f"{w:.2f}"))

# 8. [Gamma]: implicit differentiation with the -1/Theta_mix term dropped
Ts_v = np.array([C1.Th_sat, NC4.Th_sat])


def gamma_missing(iastobj, Th, Th_sat):
    Th = np.atleast_2d(Th); n = iastobj.n
    psi, x = iastobj._psi_of(Th)
    P0 = iastobj._P0(psi); T0 = iastobj._T0(P0)
    dT0 = np.stack([iastobj.isos[i].dtheta(P0[..., i]) for i in range(n)], -1) * P0 / T0
    Tt = np.sum(Th, -1); S = np.sum(x / T0, -1)
    dF_dpsi = np.sum(x * dT0 / T0 ** 2, -1) / S ** 2
    inv = 1.0 / T0
    dF_dT = -(inv / Tt[:, None] - np.sum(x * inv, -1)[:, None] / Tt[:, None]) \\
        / S[:, None] ** 2 - 1.0
    dpsi_dT = -dF_dT / dF_dpsi[:, None]
    dlnf = np.eye(n)[None] / Th[:, :, None] + dpsi_dT[:, None, :] / T0[:, :, None]
    return (Th / Th_sat)[:, :, None] * dlnf * Th_sat[None, None, :]


G_bad = gamma_missing(mix, states, Ts_v)
breaks.append(("[Gamma]: the -1/Theta_mix term dropped from dln f/dTheta",
               "vs central differences", f"{GAMMA_TWO_ROUTES:.1e}",
               f"{float(np.max(np.abs(G_bad/G_nu-1))):.2f}"))

# 9. eq. (44): the exchange sum dropped from B_ii (measured above)
breaks.append(("eq. (44): the exchange sum dropped from B_ii",
               "vs the matrix-free eq. (42)", f"{FRICTION_TWO_ROUTES:.1e}",
               f"{FRICTION_BROKEN:.1f}"))

# 10. the slab solved as a cylinder
m_cyl = Membrane([C1, NC4], [C1.Th_sat, NC4.Th_sat], [D_C1, D_NC4], n_z=30)
m_cyl.div = construct_div(m_cyl.shape, m_cyl.z_f, nu=1)
_, J_cyl = m_cyl.steady(FEEDS["95-5"])
breaks.append(("membrane solved with nu = 1 (cylinder) instead of 0",
               "S_P, 95-5", f"{BASE_SP:.2f}", f"{SP(J_cyl, FEEDS['95-5']):.2f}"))

# 11. the two M-S diffusivities swapped
breaks.append(("the two M-S diffusivities swapped", "S_P, 95-5",
               f"{BASE_SP:.2f}", f"{sp_with(D0=[D_NC4, D_C1]):.2f}"))

# 12. a Table 1 slip: nC4's site-A affinity read one decade low
bad_row = DSL(NC4.bA / 10, NC4.tA, NC4.bB, NC4.tB)
breaks.append(("Table 1: nC4 b_A read a decade low", "S_P, 95-5",
               f"{BASE_SP:.2f}", f"{sp_with(isos=[C1, bad_row]):.2f}"))

print(f"{'injected defect':>54}{'quantity':>28}{'as built':>13}{'broken':>13}")
for a, b, c_, d_ in breaks:
    print(f"{a:>54}{b:>28}{c_:>13}{d_:>13}")

CYL_SHIFT = abs(SP(J_cyl, FEEDS["95-5"]) / BASE_SP - 1)
CYL_FLUX = float(J_cyl.mean(0)[1] / RUNS[("95-5", "vignes")][2].mean(0)[1])
CLAIMED = abs(SP_955 / PRINTED["mix_SP_955_finite"] - 1)
print(f"\\nThe geometry row is the weakest in the table and is reported as such.")
print(f"Solving the slab as a cylinder shifts S_P by {CYL_SHIFT*100:.1f} %, where every")
print("other row moves by a factor. The reason is structural: both species pass")
print("through the SAME construct_div, so an area profile multiplies the numerator")
print("and the denominator of a selectivity almost identically and largely cancels.")
print(f"It is not powerless - {CYL_SHIFT/CLAIMED:.0f}x the {CLAIMED*100:.2f} % deviation "
      "actually claimed - but a")
print("selectivity is a poor quantity to test a geometry with, and the ABSOLUTE")
print(f"fluxes are the right one: the same swap multiplies the nC4 flux by "
      f"{CYL_FLUX:.0f},")
print(f"against the {FLUX_WORST*100:.2f} % those fluxes are reproduced to.")
print(f"\\nThe eq. (60) row is of a different kind: it does not break the code, it")
print("breaks the READING of a printed equation. The Elsevier text flattens the")
print("compound fraction of eqs. (60) and (61) identically, so the reading has to")
print(f"be settled some other way. Multiplying by the pressure ratio instead of")
print(f"dividing moves S_P by a factor of {(P_C1/P_NC4)**2:.0f} and lands nowhere near any")
print("printed selectivity, and the review's own worked value for eq. (60),")
print(f"{FLUX_NC4}/{FLUX_C1} divided by {P_NC4/P_C1:g}, is the division too. Equation and worked")
print("value agree; this row settles WHICH READING IS MEANT and is not a claim")
print("that the paper prints eq. (60) wrongly.")'''))

cells.append(code('''print("PRINTED DEFECTS FOUND ON THE MIXTURE SIDE. Neither is repaired silently\\n")
print("and neither changes a number on this page.\\n")
print("  READ THIS FIRST. Both readings below come from the Elsevier full text,")
print("  which is the only form of this paper on disk. There is no page image and")
print("  nothing here was read off a 600 dpi render. A text dump is a poor witness")
print("  for anything involving a FRACTION - see the eq. (60) row of the break")
print("  table, where the flattening is genuinely ambiguous and no defect follows.")
print("  These two involve only subscripts, which the dump reproduces reliably.\\n")
print("  * eq. (37) is printed as 'L_ij = L_ij; i,j = 1,2,...n (i != j)' in the")
print("    sentence stating that the Onsager reciprocal relations demand [L] be")
print("    SYMMETRIC. The relation is L_ij = L_ji; as printed it is a tautology")
print("    and says nothing. Nothing on this page uses eq. (37): the Onsager")
print("    matrix is never formed, because eq. (45) goes straight from [B] and")
print("    [Gamma] to the fluxes.")
print("  * the caption of Fig. 21 prints 'Theta_1,sat = 22 and Theta_1,sat = 13';")
print("    the second subscript must be 2. This is NOT a repair by inference -")
print("    the running text of Section 3.4 gives the same pair correctly as")
print(f"    Theta_1,sat = {TS[0]:.0f} and Theta_2,sat = {TS[1]:.0f}, so the text")
print("    supplies the reading and the caption is the defective copy. Both are")
print("    recorded in the parameter CSV's row note.")
print("\\n  A4.7 records three further printed defects, all single-component: the")
print("  sign of eq. (25), the nu_zz + nu_zz denominator of eq. (23c), and")
print("  nu_str printed twice in the CH4 jump-frequency sentence. None of the")
print("  three touches anything computed here.")'''))

cells.append(code('''report_agreement("H1.9", {
    # --- reproduction: the review's own computed results
    "SP_955_finite_rel_dev": SP_955 / PRINTED["mix_SP_955_finite"] - 1,
    "SP_955_infinite_rel_dev": SP_955_INF / PRINTED["mix_SP_955_inf"] - 1,
    "SP_5050_finite_rel_dev": SP_5050 / PRINTED["mix_SP_5050_finite"] - 1,
    "SP_5050_infinite_rel_dev": SP_5050_INF / PRINTED["mix_SP_5050_inf"] - 1,
    "SP_hexane_66kPa_rel_dev": SP_66 / PRINTED["hex_SP_nC6_3MP"] - 1,
    "flux_worst_rel_dev": FLUX_WORST,
    "flux_ratio_delta_free_worst": RATIO_FREE_WORST,
    # the two above are ONE check: delta cancels out of the first. This is the
    # residual of that algebraic identity, and it is a disclosure, not evidence.
    "flux_delta_cancellation_residual": DELTA_CANCELS,
    "peak_flux_finite_rel_dev": PEAKS["vignes"][1] / PRINTED["mix_peak_C1_finite"] - 1,
    "peak_flux_infinite_rel_dev": PEAKS["inf"][1] / PRINTED["mix_peak_C1_inf"] - 1,
    "peak_time_rel_dev": T_PEAK_DEV,
    "worst_steady_target_rel_dev": WORST_CORE,
    "worst_peak_height_rel_dev": PEAK_WORST,
    "peak_height_dt_halved_shift": DT_PEAK_SHIFT,
    "peak_time_dt_halved_shift": DT_TIME_SHIFT,
    # the transient grid study: what the coarse grid used to report, what it
    # costs at the production grid, and what closing the peak time would need
    "peak_time_rel_dev_at_nz30": gv[30][0] / PRINTED["mix_peak_t_finite"] - 1,
    "peak_flux_finite_rel_dev_at_nz30": gv[30][1] / PRINTED["mix_peak_C1_finite"] - 1,
    "peak_time_grid_err_at_production": T_GRID_ERR,
    "peak_time_grid_order": T_ORDER2,
    "delta_shift_needed_for_peak_time": DELTA_FOR_PEAK,
    "transient_front_sharpening": FRONT_T / FRONT_S,
    # --- VALIDATION: against measurement, tier 2
    "funke_mixture_rel_dev": SP_15 / MEAS["funke_mixture"] - 1,
    # ... and the band of readings the review leaves open around it
    "funke_choice_spread": FUNKE_CHOICE_SPREAD,
    "bakker_955_model_over_measured": SP_955 / MEAS["bakker_95_5"],
    "bakker_5050_model_over_measured": SP_5050 / MEAS["bakker_50_50"],
    "bakker_955_infinite_over_measured": SP_955_INF / MEAS["bakker_95_5"],
    # --- the two printed numbers this page does NOT reproduce
    "henry_selectivity_rel_dev": S_HENRY / PRINTED["sorp_S_henry"] - 1,
    "henry_tightest_escape_flux_dev": TIGHTEST,
    "langmuir_SP_over_printed": LANG_MEAN / PRINTED["mix_SP_langmuir"],
    "langmuir_composition_spread": LANG_SPREAD,
    "iast_composition_spread": IAST_SPREAD,
    # --- internal checks
    "iast_roundtrip": IAST_ROUNDTRIP,
    "surface_potential_vs_quadrature": PSI_QUAD,
    "iast_vs_langmuir_equal_capacities": IAST_LANGMUIR,
    "gamma_two_routes": GAMMA_TWO_ROUTES,
    "friction_system_two_routes": FRICTION_TWO_ROUTES,
    "friction_system_broken": FRICTION_BROKEN,
    "n1_collapse_solver_vs_eq12": COLLAPSE_1,
    "grid_rel_err_at_nz40": GRID_ERR,
    "sweep_pressure_sensitivity": PD_SENS,
    "loading_dependence_spread": DEP_SPREAD,
    "continuation_path_dev": PATH_DEV,
    "delta_reconstructed_um": DELTA * 1e6,
})'''))

# ---------------------------------------------------------------- adds
cells.append(md(r"""## What pymrm adds

The review's mixture results are reproduced, not improved on. Six things here
are not in it.

**1. The review's mixture calculations are checked against measurement in the
words the distinction deserves.** Krishna & Baur state four measured
selectivities in running text — Bakker's 380 and 60, Funke's 24 and 1.3 — and
compare against them qualitatively. Put on one axis, the model with finite
exchange coefficients is 1.28× Bakker's 95–5 measurement and **3.4× the 50–50
one**, which the review quotes in the following sentence and does not comment on.
The infinite-exchange model is 31× and 11× too *low*. So the measurement
discriminates decisively between the two models — the gap between them is two
orders of magnitude and the data sit unambiguously on one side of it — while
agreeing quantitatively with neither. Funke's hexane-isomer case is different
and better: 25.0 computed against 24 measured, +4.3 %, with nothing fitted.
**That is the only quantitative validation against measurement in the review's
mixture sections**, and this page is where it is scored — with its error budget
printed beside it, because "nothing fitted" is not "nothing chosen". The review's
one sentence about Funke's run fixes neither the feed composition nor whether its
15 kPa is a total or a per-component pressure, and the diffusivity ratio of 5
comes from a third paper that was not consulted. All three readings are swept at
the point of comparison, and the alternatives span the 1.6× band printed there
against a +4.3 % agreement. That sweep is what makes the +4.3 % readable.

**2. The Henry-limit sorption selectivity of 2200 does not follow from Table 1,
and the obvious explanation is excluded.** The dual-site Henry coefficients give
2649.9 against the printed 2200, and the IAST solve reaches the same limit from
an entirely different code path — so it is not an arithmetic slip in the Henry formula. The
tempting explanation is a transcription error in one of the four C1/nC4 Table 1
parameters, and that can be tested rather than assumed: solving for the value
each parameter would need to make the ratio 2200 and feeding it back into the
one check that certifies those rows — `A4.7`'s single-component flux ratio,
observed at 0.17 % — moves that check by 4.2 % at best and 24 % at worst. No
single-character slip explains 2200 without breaking the certification. What
this does *not* settle is where 2200 comes from: the review's Fig. 24 is drawn
against CBMC fugacities, so a Henry coefficient taken from the simulation rather
than from the Table 1 fit remains possible and nothing printed decides it. The
subsidiary claim — that the selectivity is "practically constant" below
$\Theta_{\rm mix}=8$ — is fair on a log axis: IAST falls 32 % over that range
against 1.9 decades afterwards.

**3. `A4.7`'s reconstructed membrane thickness is tested on the one power of
$\delta$ that this page actually carries, and it fails by about 10 %.** The
reconstruction inverted the two *single-component* steady fluxes, where
$N\propto1/\delta$. An earlier version of this page claimed the four *mixture*
fluxes retest it out of sample on that same power. **They do not:** $\delta$,
$\rho_{\rm MFI}$, the unit-cell mass and Avogadro's number cancel identically
out of all four, the conversion collapsing to
$N_i=34\,J_i/(\eth_{\rm C1}\psi_{\rm C1})$, and swinging the framework density
over a factor of four leaves every one of the four bit-identical. That
comparison is the $\delta$-free flux-ratio check printed ten lines below it,
rescaled — the page now prints both deviations side by side so the coincidence
is visible, and keeps the check for what it does test, which is the mixture
physics. The transient peak *time* is the only quantity here that carries
$\delta$, as $\delta^2$. **Refined under the transient — which the steady grid
study never did, because the transient carries a self-sharpening butane front and
the steady profile does not — it lands about 10 % below the printed 0.73 s, and
both further grid and step refinement move it further away.** Closing that needs
$\delta$ about 5 % larger than 40 µm, which nothing else on this page forbids,
precisely because nothing else on this page contains $\delta$. It is reported as
an open tension with its four candidate homes named, not as a reconciliation.
The two peak **heights** move the same way: at $n_z=30$ they sat below the
printed values and the sign was attributed to backward Euler damping a maximum;
converged, both sit above, and that attribution is withdrawn.

**4. Two of the review's structural claims are separated from its numbers, and
they do not stand or fall together.** The multicomponent Langmuir with equal
saturation capacities really does make $S_P$ independent of upstream
composition — to 1.02× over feeds running from 95 % down to 5 % methane, where
IAST spans more than a factor of three over exactly the same five feeds —
which
is precisely the review's argument for IAST. Its printed *value* of 800 is
not reproduced: the computation gives about 1680, and the review does not say
which affinity it used once the dual-site fit is collapsed to a single site. The claim
is right and the number is underdetermined, and nothing was tuned to reach it.
The same separation applies to the hexane thresholds: "$\Theta_{\rm mix}>4$
above 1 kPa" is true but loose by 6.4×, and "$S_P>10$ above 2 kPa" crosses at
2.5 kPa.

**5. Every check on the coupled machinery has an independent second route, and
the two identities are labelled instead of scored.** $[\Gamma]$ is obtained by
implicit differentiation of eqs. (62)–(65) and checked against central
differences of the IAST inverse map; $[B]$ is assembled from eq. (44) and also
never assembled at all, the friction system of eq. (42) being solved matrix-free
by GMRES — `A4.2`'s fix for a check that compared $[B]$ against itself, reused,
and it moves by hundreds of per cent when a term is dropped from $B_{ii}$. The surface
potential is checked against adaptive quadrature, which never sees the closed
form. Against that, eq. (43) at $i=j$ giving $\eth_{ii}=\eth_i$ and eq. (55)
following from eq. (54) are **algebraic identities**, and so is eq. (46)
collapsing to eq. (14) at $n=1$; all three are stated as identities and none is
scored. The collapse that *can* fail — the $n=1$ solver against the closed-form
eq. (12) integral — is the one reported.

**6. The review's own counter-example is carried, not quietly dropped.** MFI
needs finite exchange coefficients; the review's conclusion for zeolite 4A is
the opposite, that Habgood's N₂/CH₄ uptake is better described by
$\eth_{ij}\to\infty$ because intracage hopping there is fast. The finding is not
that finite $\eth_{ij}$ wins — it is that $\eth_{ij}$ is a property of the host
topology and these two hosts sit at opposite limits.

**What the page does not add and cannot.** The 4A uptake, the O₂/N₂
chromatography and the C₅/C₆ breakthroughs of Sections 4.4–4.7 are all mixture
problems whose only validation is a figure, and the same is true of the MD/KMC
verification of Figs. 17–21: the mixture theory is implemented and exercised
here, but there is no printed number behind Sanborn & Snurr's Fick matrix or
Snurr & Kärger's self-diffusivities, so the self-diffusivity section shows the
size and sign of the correlation effect and claims nothing about whether the
review is right about it. Digitising any of those needs a maintainer review and
none is available. The loading dependence of $\eth_i$ inside the mixture is
never stated by the review; eq. (17) is used because it is what reproduces the
printed selectivities and what `A4.7` established for the single components, and
the spread across the three plausible readings — which turns out to be only a
few per cent — is measured and reported rather than assumed away."""))

# ---------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**`IAST` is the reusable object, and it is more general than this page.** It
takes any list of pure-component isotherms exposing $\Theta^0(P)$, its
derivative and the surface potential $\int\Theta^0\,{\rm d}P/P$, and solves
eqs. (62)–(65) in both directions for any number of components. Swap `DSL` for a
Toth or a Sips and nothing else changes. The inverse direction — loadings to
pressures — is the one a transport model needs and the one textbooks rarely
write down; here it is a *scalar* monotone root find, because the adsorbed-phase
mole fractions are already determined by the loadings.

**`inverse_and_gamma` is the piece worth lifting.** Almost every published IAST
transport model differentiates the isotherm numerically to get $[\Gamma]$, which
costs $2n$ nonlinear solves per evaluation point and puts a finite-difference
step size inside a Newton iteration. Implicit differentiation of the IAST
constraint gives every element in closed form from **one** solve, and the
numerical route is kept as the check rather than as the method.

**`Membrane` is a general $n$-component coupled-diffusion slab.** The only
zeolite-specific parts are `B_matrix` (the friction model) and the isotherm.
Replace $[B]^{-1}[\Gamma]$ with an ordinary Fick matrix and it is a
multicomponent membrane; replace the slab `nu=0` with `nu=2` and it is
[`A4.7`](../A4.7-zeolite-micropore-maxwell-stefan/)'s crystallite with mixtures
in it, which is what Sections 4.4–4.7 need.

**When the coupling matters, and when it does not.** If your saturation
capacities are equal, the multicomponent Langmuir is exact (checked here to
machine precision, at every loading) and IAST is wasted effort. If they differ —
19 against 10 here — the Langmuir is wrong by more than an order of magnitude at
high loading and the error is *qualitative*: no Langmuir isotherm can produce a component loading that
falls as the pressure rises, which is the effect the whole separation lives on.
If your sorbates have similar mobilities the exchange coefficients barely
matter; the further apart they are the more they do, and at a hundredfold ratio
they are worth a factor of 40 in selectivity.

**How to tell which limit you are in, without a mixture experiment.** The review
answers it topologically: exchange happens where the molecules meet, so a
structure with large cages and fast intracage hopping (LTA, FAU) tends to
$\eth_{ij}\to\infty$ while an intersecting-channel structure (MFI) does not.
Both limits are one argument apart in the code, and running both and reporting
the spread is cheaper than guessing.

**Where this page starts, and where `A4.7` ends.** Everything single-component —
the flux law eq. (12), the scalar $\Gamma$, the Kärger zero-loading relations,
the transient uptake and the LDF constant — is
[`A4.7`](../A4.7-zeolite-micropore-maxwell-stefan/)'s and is not repeated here.
The cut between the two cases is **single-component / mixture**, not micropore /
membrane: the mixture side needs a matrix $[\Gamma]$, exchange coefficients, an
$n\times n$ friction system and a non-linear IAST solve, none of which appears on
`A4.7`; it carries a different headline; and it is where the review's
experimental comparisons live, which should not be a subsection of somebody
else's page. `A4.7` keeps Section 4.1's two single-component fluxes because they
are the only place in the review where the single-component flux law produces a
printed number, and the flux-ratio check is the only one on that page exercising
$\eth$ and $\Gamma$ *inside a flux*.

**Related pages.** [`A4.2`](../A4.2-maxwell-stefan-vs-fick/) (the same $[B]$,
$[\Gamma]$, $[D]$ algebra in a bulk gas, and the source of the matrix-free
friction check), [`A4.3`](../A4.3-dusty-gas-model/) (the porous-medium form),
[`A4.9`](../A4.9-duncan-toor/) (what a non-diagonal $[D]$ does to a ternary
mixture), [`A4.7`](../A4.7-zeolite-micropore-maxwell-stefan/),
[`J1.5`](../J1.5-ldf-breakthrough/), and a packed-bed chromatography page in
section J once one exists — Sections 4.5–4.7 are waiting for it.

**Cite the source, not this page:** Krishna, R. & Baur, R., *Modelling issues in
zeolite based separation processes*, Separation and Purification Technology
**33**(3) 213–254 (2003),
[doi:10.1016/S1383-5866(03)00008-X](https://doi.org/10.1016/S1383-5866(03)00008-X).

The measured selectivities are **not** the review's own: they are Bakker's PhD
thesis and Funke et al., read here through the review, which is the document
actually consulted. Both origins are recorded under `origin_not_consulted` in
the dataset sidecar and neither was opened. So is **Cavalcante & Ruthven
(1995)**, which is not the origin of any measured row but *is* the origin of the
one diffusivity ratio the whole hexane result turns on; the review states the
ratio and the attribution, and nothing here checks either.

The catalogue's other named reference for this case, "Krishna & van den Broeke
(1995)", is not the source of this page: that paper is van den Broeke & Krishna,
*Experimental verification of the Maxwell–Stefan theory for micropore
diffusion*, Chem. Eng. Sci. **50**(16) 2507–2522, and its adsorbents are
activated carbon and a carbon molecular sieve — no membrane and no zeolite."""))

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
