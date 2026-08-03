#!/usr/bin/env python3
"""Generate index.ipynb for page A2.5. Run from the page directory.

Note on quoting: every code cell below is a RAW triple-quoted string, so a
single backslash in this file is a single backslash in the notebook. Do not
convert them to ordinary strings - `\\nu`, `\\times`, `\\bar`, `\\frac`,
`\\varepsilon` and `\\approx` all begin valid Python escapes and would be eaten.
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "Where the axial dispersion coefficient comes from: Edwards & Richardson's packed bed"
description: "A pulse of argon in air, two response curves, and the two-point moment method that turns them into D_L — reproduced in pymrm, against the points digitised from their Figure 9."
categories: [sec:A, struct:S4, tier:T1, data:tier2, phase:gas]
date: 2026-08-02
---

# Where the axial dispersion coefficient comes from

**Catalog ID:** `A2.5` · **Structures:** `S4` (1-D transient advection–dispersion) ·
**Tier:** T1

[`A2.1`](../A2.1-danckwerts-boundary-conditions/) takes the axial-dispersion
reactor as given and asks what its boundary conditions mean.
[`A2.3`](../A2.3-taylor-aris-dispersion/) derives a dispersion coefficient for a
tube from first principles. This page is about the third question, the
experimental one: **for a real packed bed, where does the number come from, and
how good is it?**"""))

cells.append(md(r"""## Background

Edwards and Richardson pushed air through beds of glass ballotini, sand and
powdered perspex, injected a pulse of argon, and watched it come out on a
$\beta$-ray ionisation detector. They did this over a Reynolds-number range of
**0.008 to 50** — wider than anyone before them — and produced the correlation
that is still quoted for gas-phase axial dispersion in packed beds.

Three things make this a page rather than a curve fit.

**The measurement is an inverse problem, and its forward model is the thing the
gallery already builds.** A dispersion coefficient is never observed. What is
observed is a concentration against time at two positions; $D_L$ is whatever
value of the parameter in

$$\frac{\partial c}{\partial t} + u\frac{\partial c}{\partial x}
  = D_L\frac{\partial^2 c}{\partial x^2}$$

makes the second curve follow from the first. So a page that can solve that
equation can *simulate the whole experiment*, invert it with the authors' own
formulae, and check that the number it gets back is the number it put in. That
check has real power, because the discretisation adds dispersion of its own and
the inversion cannot tell it apart from the physical kind.

**The correlation has structure, not just constants.** It is a molecular
diffusion floor plus an eddy-diffusion term, and each has a limit the model must
hit: at low Reynolds number $D_L$ must approach the molecular diffusivity times
a tortuosity factor, and at high Reynolds number the Péclet number must approach
2. In between, the data show a *maximum* in Péclet number that the sum of those
two mechanisms cannot produce at all — and that is what the paper's third term
is for.

**They also found where their own model fails.** With their finest particles the
moment analysis returned dispersion coefficients three times too large, and the
response curve computed from the fitted $D_L$ peaked *later* than the one they
measured. They attributed it to channelling. That diagnostic — the shape
disagreeing while the moments agree by construction — is reproduced in the last
section."""))

cells.append(md(r"""## The published model

### The measurement

For a pulse in an infinite uniform bed, Aris showed that two response curves
measured a distance $\Delta x$ apart are related by three exact statements,
which are Edwards & Richardson's Eqs. (5), (6) and (7):

$$A_1 = A_2, \qquad
\Delta\bar t = \bar t_2 - \bar t_1 = \frac{\Delta x}{u}, \qquad
\Delta\sigma^2 = \sigma_2^2 - \sigma_1^2 = \frac{2 D_L \Delta x}{u^3},$$

with area, mean and variance defined by their Eqs. (2)–(4),

$$A=\int_0^\infty c\,\mathrm dt,\qquad
\bar t=\frac1A\int_0^\infty c\,t\,\mathrm dt,\qquad
\sigma^2=\frac1A\int_0^\infty c\,t^2\,\mathrm dt-\bar t^{\,2}.$$

Everything that happens *before* the first detector — the shape of the injected
pulse, the distributor, the entrance section — cancels in the differences. That
is the whole point of measuring twice, and it is why the method needs no
knowledge of the input.

The three relations are exact for the equation above, not approximations. In
Laplace space the downstream transfer function is $\exp(\lambda\Delta x)$ with
$\lambda = \bigl(u-\sqrt{u^2+4D_Ls}\bigr)/2D_L = -s/u + D_Ls^2/u^3 + O(s^3)$,
and the first three coefficients of that expansion are exactly the three
statements above.

### The correlation

Two mechanisms, taken as additive (their Eqs. 13–15):

$$D_L = \gamma D_M + \tfrac12 u d_p .$$

The first term is molecular diffusion reduced by the tortuosity of the bed; the
second is the mixing-cell and turbulent-mixing limit $\mathrm{Pe} = u d_p/D_L =
2$ of Aris & Amundson and of Prausnitz. This form rises monotonically to
$\mathrm{Pe}=2$ and can never exceed it — but the measurements do. Radial mixing
inside the bed reduces axial spreading, and E&R put that in as an empirical
correction factor on the eddy term (their Eq. 16):

$$D_L = \gamma D_M + \frac{\tfrac12 u d_p}{1 + \beta D_M/(u d_p)} .$$

With $\gamma$ measured at 0.73 and $\beta$ fitted at 9.7 this is their Eq. (17),
quoted in the abstract as the result of the paper:

$$D_L = 0.73\,D_M + \frac{0.5\,u d_p}{1 + 9.7\,D_M/(u d_p)},
\qquad
\begin{aligned}&0.008 < \mathrm{Re} < 50\\ &0.0377 < d_p < 0.60\ \mathrm{cm}\end{aligned}$$

and, divided through by $u d_p$ with $\varepsilon = 0.37$ and $\mathrm{Sc}=0.72$,
their Eq. (18):

$$\mathrm{Pe}^{-1} = 0.38\,\mathrm{Re}^{-1} + \frac{0.5}{1+5.0\,\mathrm{Re}^{-1}} .$$

**Definitions matter here.** $u$ is the *interstitial* velocity throughout, so
$\mathrm{Re} = u d_p \varepsilon/\nu$ is built on the superficial velocity
$u\varepsilon$ while $\mathrm{Pe} = u d_p / D_L$ is not. Mixing the two is the
easiest way to be wrong by a factor $1/\varepsilon \approx 2.7$."""))

cells.append(md(r"""## Parameters and assumptions

**The system is one gas pair only.** Argon traced into air, $D_M = 0.205$ cm²/s.
The paper says plainly that "only one system (argon–air) has been studied and
thus the effect of variations in the kinematic viscosity term $\nu$ in the
Reynolds number has not been established". Nothing on this page tests it either.

**Assumptions carried into every calculation below.** One-dimensional flow with
no radial concentration gradient (the paper's porous distributor is what buys
this); constant $u$ and constant $D_L$; a non-porous, non-adsorbing solid, so no
intraparticle hold-up; and an infinite bed, in the sense that the pulse never
reaches either end while it is being observed.

**$\nu$ is not printed, but it is implied.** The paper gives $\mathrm{Sc} =
0.72$ and $D_M = 0.205$ cm²/s, so $\nu = \mathrm{Sc}\,D_M = 0.1476$ cm²/s. Every
velocity below comes from that. It is a reconstruction from two printed numbers,
not an outside value.

**Two runs are worked in detail**, both taken from figure legends: the
large-particle check of Figure 13 ($\mathrm{Re} = 24.6$, $d_p = 0.300$ cm,
$\Delta x = 100$ cm) and the fine-particle failure of Figure 15 ($d_p = 0.0097$
cm, $\mathrm{Re} = 0.0710$, apparent $D_L = 0.503$ cm²/s).

**Figure 15's legend prints no detector separation**, and the paper's test
section was one of two tubes, 21.3 cm and 100 cm. Where the Figure-15
reconstruction needs a $\Delta x$ it assumes the 100 cm one, which is what
Figure 13 prints; the section that uses it says so and reports the 21.3 cm
alternative beside it. Nothing else on the page depends on the choice."""))

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
from scipy.sparse import eye_array, diags_array, coo_array
from scipy.sparse.linalg import splu
from scipy.special import erfc
from scipy.optimize import least_squares
from pymrm import (construct_grad, construct_div, construct_convflux_upwind,
                   interp_cntr_to_stagg_tvd, vanleer)
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A2.5-edwards-richardson-dispersion"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
SEED = 20260802                # every resampling below seeds from this, explicitly


# Deviations are reported ONE way everywhere on this page:
#     dev = (model - measured) / measured
def dev(model, measured):
    return np.asarray(model, float) / np.asarray(measured, float) - 1.0'''))

cells.append(md(r"""## The data

Three files, of two very different kinds.

**Printed numbers.** The constants, ranges and worked-run conditions, and
Table 1's particle properties. These were read off 600 dpi page renders because
the text layer of this 1968 Pergamon scan drops the mid-dot decimal separator —
it renders $d_p = 0.607$ cm as `0607` and $d_p = 0.0097$ cm as `o%l97`. Nothing
was repaired by inference.

**Digitised measurements.** Edwards & Richardson print no table of dispersion
coefficients; the measurements exist only as Figures 7–9 and 14. Figure 9, the
Péclet–Reynolds plot for the five large bead sizes, was digitised — marker
positions only, and that is what makes this page tier 2 rather than tier 6. The
count and the Reynolds range come out of the file below, not out of this
sentence.

Two things made that digitisation cheap and safe. The dashed curve drawn through
the points *is* Eq. (18), whose constants are printed, so its pixel path could be
computed rather than traced and then erased before marker detection — and the
fact that the computed path lands on the printed dashes over 3.7 decades is
itself a check on both the axis calibration and the transcription. And Eq. (18)
contains no particle-size term, so the five marker shapes carry no information
this test needs: only positions were extracted and **no row carries a series
label**.

**Where the extraction needed help, and where it failed.** In four places the
glyphs overlap by more than half their width, and there the ink centroid of the
blob belongs to more than one marker. Those four clusters were redone as units:
every row inside them discarded, and one centre per glyph put back at the
best-fit position of that glyph's *own* template, seeded by a visual audit of
the 600 dpi render. Ten rows went in and eleven came out — six of them moved by
9 to 27 px, one came off a bare stretch of the dashed curve and onto the open
circle beside it, and one marker that had carried no row at all was added. The
sidecar carries the boxes, the old rows and the new centres, and the extraction
script reproduces the file from the page render. Three glyphs are still fused
into their neighbours badly enough that no centre can be fixed for them; they
are listed in the sidecar and printed below, and the sidecar records what
including the most nearly resolvable of them would do to the headline."""))

cells.append(code(r'''fig9 = load_data("edwards-richardson-1968-fig9.csv", page=PAGE)
stated = load_data("edwards-richardson-1968-stated.csv", page=PAGE)
tab1 = load_data("edwards-richardson-1968-table1.csv", page=PAGE)
P = dict(zip(stated.quantity, stated.value))          # printed values, by name

meta9 = load_meta("edwards-richardson-1968-fig9.csv", page=PAGE)
print(cite_data(meta9))
print(f"\nFigure 9: {len(fig9)} markers, "
      f"Re {fig9.reynolds.min():.5f} to {fig9.reynolds.max():.3g}, "
      f"Pe {fig9.peclet.min():.4f} to {fig9.peclet.max():.3f}")
print(f"the paper's own stated range: Re {P['reynolds_min']} to {P['reynolds_max']}"
      "  (never used in the axis calibration)")

# What the extraction could NOT do, read from the sidecar rather than retyped.
gone = meta9["missing"]["unrecovered_markers"]
print(f"\n{len(gone)} glyphs are fused into a neighbour and carry no row in the file,"
      " at approximately")
for m in gone:
    print(f"   Re {m['reynolds']:6.2f}  Pe {m['peclet']:5.2f}"
          f"   ({m['leftover_ink_px']} px of unexplained ink at"
          f" ({m['pixel_x']}, {m['pixel_y']}))")
print(f"   -> {meta9['missing']['effect_if_included']}")

print("\nTable 1, the five sizes Figures 7 and 9 cover:")
big = tab1[tab1.d_p_cm >= P["particle_diameter_min"]]
print(big[["material", "d_p_cm", "voidage", "material_density_g_cm3",
           "bulk_density_g_cm3"]].to_string(index=False))

# A check Table 1 pays for: 1 - bulk/material density must be the printed voidage.
# It tests three transcribed columns against each other on all eight rows.
void_from_rho = 1.0 - tab1.bulk_density_g_cm3 / tab1.material_density_g_cm3
resid_void = (void_from_rho - tab1.voidage).abs()
print(f"\nvoidage against 1 - bulk/material density, all {len(tab1)} rows:"
      f" worst absolute residual {resid_void.max():.4f}"
      f" (on d_p = {tab1.d_p_cm[resid_void.idxmax()]} cm),"
      f" mean {resid_void.mean():.4f}")

print("\nprinted constants used below:")
for k in ("tortuosity_factor_gamma", "correction_constant_beta",
          "molecular_diffusivity_argon_air", "voidage_used_in_eq18",
          "schmidt_number_used_in_eq18", "dispersion_coefficient_low_Re"):
    row = stated[stated.quantity == k].iloc[0]
    print(f"  {k:38s} {row.value:>9g} {row.unit:8s} (p. {row.page})")'''))

cells.append(md(r"""### Five checks the paper pays for, before any modelling

A page-image transcription is a transcription and needs checking like any other
(the `B3.1` lesson). This paper hands over five ways to do it, none of which
needs the figure or a solver, and each of which can fail on a single mis-read
digit. The validation section injects exactly such digits and shows the numbers
move. One of the five does *not* come out, and it is reported as it stands."""))

cells.append(code(r'''gamma = P["tortuosity_factor_gamma"]
beta = P["correction_constant_beta"]
D_M = P["molecular_diffusivity_argon_air"]
eps18 = P["voidage_used_in_eq18"]
Sc18 = P["schmidt_number_used_in_eq18"]
nu_air = Sc18 * D_M                       # reconstruction: nu = Sc * D_M
C_EDDY = P["eq18_eddy_coefficient"]       # 0.5, from Pe -> 2


def D_L_eq17(u, dp, D_M=D_M, gamma=gamma, beta=beta):
    """E&R Eq. (17): molecular floor plus corrected eddy term. u is INTERSTITIAL."""
    return gamma * D_M + 0.5 * u * dp / (1.0 + beta * D_M / (u * dp))


def Pe_eq18(Re, A=None, B=None, C=C_EDDY):
    """E&R Eq. (18): Pe^-1 = A/Re + C/(1 + B/Re)."""
    A = P["eq18_molecular_coefficient"] if A is None else A
    B = P["eq18_correction_coefficient"] if B is None else B
    return 1.0 / (A / Re + C / (1.0 + B / Re))


def Pe_max_eq18(A, B, C=C_EDDY):
    """Where d(Pe^-1)/dRe vanishes: C B/(Re+B)^2 = A/Re^2."""
    r = np.sqrt(A / (C * B))
    Re_s = r * B / (1.0 - r)
    return Re_s, Pe_eq18(Re_s, A=A, B=B, C=C)


A_printed = P["eq18_molecular_coefficient"]
B_printed = P["eq18_correction_coefficient"]

# --- 1. gamma is 0.150 / 0.205 -------------------------------------------
gamma_ratio = P["dispersion_coefficient_low_Re"] / D_M
gamma_ci = P["dispersion_coefficient_low_Re_ci95"] / D_M

# --- 2. Eq. (18)'s constants are Eq. (17)'s, scaled by eps/Sc ------------
#     D_M/(u d_p) = (nu/Sc)/(Re nu/eps) = eps/(Sc Re), so the two coefficients
#     of Re^-1 must be gamma*eps/Sc and beta*eps/Sc.
k18 = eps18 / Sc18
A_derived, B_derived = gamma * k18, beta * k18

# --- 3. where eps = 0.37 comes from --------------------------------------
eps_mean = float(big.voidage.mean())

# --- 4. the maximum of Eq. (18), which the paper only bounds -------------
Re_star, Pe_star = Pe_max_eq18(A_printed, B_printed)

print(f"1. gamma from the low-Re plateau: {P['dispersion_coefficient_low_Re']} / {D_M}"
      f" = {gamma_ratio:.4f} +/- {gamma_ci:.4f}, printed as {gamma}")
print(f"2. Eq.(18) coefficients from Eq.(17), scaled by eps/Sc = {k18:.5f}:")
print(f"     {gamma} x {k18:.5f} = {A_derived:.5f}   printed {A_printed}"
      f"   ({abs(A_derived - A_printed):.4f} absolute,"
      f" {abs(dev(A_derived, A_printed)):.2%})")
print(f"     {beta} x {k18:.5f} = {B_derived:.5f}   printed {B_printed}"
      f"   ({abs(B_derived - B_printed):.4f} absolute,"
      f" {abs(dev(B_derived, B_printed)):.2%})")
print(f"3. mean voidage of the five large sizes in Table 1 = {eps_mean:.4f},"
      f" printed as e = {eps18}")
print(f"4. maximum of Eq.(18): Pe = {Pe_star:.4f} at Re = {Re_star:.4f}")
print(f"     the conclusions say only 'a maximum, greater than"
      f" {P['peclet_maximum_claim']:.0f}'")
print(f"     Gunn (1993) Fig. 1 plots these data together with Gunn & Pryce's")
print(f"     frequency-response set and puts the maximum of the COMBINED gas-phase")
print(f"     results 'at a Reynolds number of ~4' - not of E&R's points alone")
print(f"   high-Re limit of Eq.(18): Pe -> 1/C = {1 / C_EDDY:.1f}, printed as"
      f" 'about {P['peclet_high_Re_limit']:.0f}'")
print(f"   low-Re limit of Eq.(17):  D_L -> gamma D_M = {gamma * D_M:.5f} cm2/s,"
      f" measured {P['dispersion_coefficient_low_Re']}")'''))

cells.append(code(r'''# --- 5. the crossover Reynolds number, and the one check that does NOT work.
#     Molecular diffusion "predominates" where gamma D_M exceeds the uncorrected
#     eddy term u d_p/2, i.e. at Re = 2 gamma eps D_M / nu.
def Re_crossover(D_mol, nu_fluid, gamma=gamma, eps=eps18):
    return 2.0 * gamma * eps * D_mol / nu_fluid


Re_gas = Re_crossover(P["typical_diffusivity_gas"], nu_air)
Re_gas_wrong_gamma = Re_crossover(P["typical_diffusivity_gas"], nu_air, gamma=0.63)
nu_implied = (2.0 * gamma * eps18 * P["typical_diffusivity_liquid"]
              / P["crossover_reynolds_liquid"])
NU_WATER_20C = 0.0100          # cm2/s, water at 20 C - an OUTSIDE number
Re_liq_water = Re_crossover(P["typical_diffusivity_liquid"], NU_WATER_20C)

display(Markdown(rf"""
**The gas crossover comes out to {abs(dev(Re_gas, P['crossover_reynolds_gas'])):.1%};
the liquid one does not come out at all.**
Using only numbers printed in the paper — $\gamma$ = {gamma}, $\varepsilon$ =
{eps18}, the "typical" gas diffusivity {P['typical_diffusivity_gas']} cm²/s and
$\nu = \mathrm{{Sc}}\,D_M$ = {nu_air:.4f} cm²/s — the Reynolds number at which
molecular diffusion stops dominating is **{Re_gas:.3f}**, against the
**{P['crossover_reynolds_gas']}** the paper prints: a deviation of
{dev(Re_gas, P['crossover_reynolds_gas']):+.1%}, which is what rounding "about
1.8" to two figures buys and is *not* an exact reproduction. Nothing was tuned,
though: a $\gamma$ of 0.63 instead of 0.73 would give {Re_gas_wrong_gamma:.2f},
{abs(dev(Re_gas_wrong_gamma, P['crossover_reynolds_gas'])):.0%} away.

The same expression with their liquid diffusivity $10^{{-5}}$ cm²/s reproduces
their printed $3\times10^{{-4}}$ only for a kinematic viscosity of
**{nu_implied:.4f} cm²/s**. Water at 20 °C is {NU_WATER_20C} cm²/s and would
give {Re_liq_water:.2e}, a factor
{Re_liq_water / P['crossover_reynolds_liquid']:.1f} larger. The liquid figure
therefore implies a liquid roughly {nu_implied / NU_WATER_20C:.1f} times as
viscous as water, or a different convention. **The page states that rather than
repairing it**; the number is printed, and it is used nowhere below.
"""))'''))

cells.append(md(r"""## PyMRM implementation

Two models, both `S4`, both assembled once and reused.

**`TracerBed`** is Eq. (1) on a bed long enough that the pulse never reaches
either end. `construct_convflux_upwind` for $u\,\partial c/\partial x$,
`construct_grad` and `construct_div` with `nu=0` for the Cartesian slab, a van
Leer TVD **deferred correction** on top of the first-order upwind flux, and a
$\theta$-method in time so that implicit Euler ($\theta=1$) and Crank–Nicolson
($\theta=\tfrac12$) can be compared. The operator is constant, so it is
factorised once and every time step is a back-substitution; the deferred
correction is warm-started from the previous step, so each step starts from the
previous step's correction rather than from scratch. How much that saves is not
measured here and no figure for it is quoted; what *is* printed, at the end of
the validation section, is the distribution of iteration counts actually needed
— median, the fraction that finish in five or fewer, and the worst case, which
is not small.

Every keyword after `limiter` exists only so that a defect can be injected in
the break table. The defaults are the model as published.

**`TwoZoneBed`**, defined later where it is first needed, is the same equation
on two mobile zones that exchange laterally. Layout `(n, 2)` — spatial axis
first, fields last — with the volume fractions in the mass matrix and in the
face velocities."""))

cells.append(code(r'''CORR_ITERS = []          # deferred-correction iteration counts, asserted at the end


class TracerBed:
    """dc/dt + u dc/dx = D d2c/dx2 on 0 <= x <= L, uniform grid.

    The pulse is released in the interior and must never reach either end, so
    the boundary conditions only have to be harmless: no tracer enters at x = 0,
    and the dispersive flux vanishes at x = L.
    """

    def __init__(self, L, n, u, D, limiter=vanleer, nu=0, D_sign=1.0,
                 outlet="zero-gradient", tol=1e-9, max_corr=40):
        self.n, self.u, self.D, self.L, self.h = n, u, D, L, L / n
        self.shape = (n, 1)                       # (cells, fields) - never (n,)
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # inlet, outward normal n = -x:   c = 0      ->  a = 0, b = 1, d = 0
        bc_in = {"a": 0.0, "b": 1.0, "d": 0.0}
        # outlet, outward normal n = +x:  dc/dx = 0  ->  a = 1, b = 0, d = 0
        bc_out = ({"a": 1.0, "b": 0.0, "d": 0.0} if outlet == "zero-gradient"
                  else {"a": 0.0, "b": 1.0, "d": 0.0})
        self.bc = (bc_in, bc_out)
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                  self.bc, v=u)
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc)
        div = construct_div(self.shape, self.x_f, nu=nu)      # nu=0: Cartesian slab
        self.div = div
        self.A = div @ (conv - D_sign * D * grad)
        self.b0 = np.asarray((div @ (conv_bc - D_sign * D * grad_bc)).todense()).ravel()
        self.limiter, self.tol, self.max_corr = limiter, tol, max_corr
        self._key, self._last_corr = None, np.zeros(n)

    def _factor(self, dt, theta):
        if self._key != (dt, theta):
            self._lu = splu((eye_array(self.n, format="csc") / dt
                             + theta * self.A).tocsc())
            self._key = (dt, theta)
        return self._lu

    def step(self, c, dt, theta=1.0, strict=True):
        lu = self._factor(dt, theta)
        rhs = c / dt - (1.0 - theta) * (self.A @ c) - self.b0
        if self.limiter is None:
            return lu.solve(rhs)
        c_new = lu.solve(rhs - self._last_corr)       # warm start
        done, it, corr = False, 0, self._last_corr
        for it in range(1, self.max_corr + 1):
            _, dc_f = interp_cntr_to_stagg_tvd(c_new.reshape(self.shape), self.x_f,
                                               self.x_c, self.bc, self.u,
                                               tvd_limiter=self.limiter, axis=0)
            corr = np.asarray(self.div @ (self.u * dc_f.reshape(-1, 1))).ravel()
            cc = lu.solve(rhs - corr)
            done = np.max(np.abs(cc - c_new)) <= self.tol * max(1.0, np.max(np.abs(cc)))
            c_new = cc
            if done:
                break
        # A deferred correction that silently returns its iteration cap is the
        # classic way to publish an unconverged number. Refuse to.
        if strict:
            assert done, f"deferred correction did not converge in {self.max_corr}"
        CORR_ITERS.append(it)
        self._last_corr = corr
        return c_new


def emg(x, x0, sigma, tail):
    """Gaussian with an exponential tail, which is what E&R's injected pulse is:
    their Fig. 5 shows 'a long tail' and their Eq. (8) fits it as an exponential
    decay. z > 0 lies behind the pulse, so it arrives later."""
    z = x0 - x
    return (0.5 / tail) * np.exp(0.5 * (sigma / tail) ** 2 - z / tail) * \
        erfc((sigma / tail - z / sigma) / np.sqrt(2))


def moments(t, c):
    """E&R Eqs. (2)-(4)."""
    A = np.trapezoid(c, t)
    m1 = np.trapezoid(c * t, t) / A
    s2 = np.trapezoid(c * (t - m1) ** 2, t) / A
    return A, m1, s2


def invert_two_point(t, c1, c2, dx):
    """E&R Eqs. (6) and (7): the whole measurement, in two lines."""
    A1, t1, s1 = moments(t, c1)
    A2, t2, s2 = moments(t, c2)
    u_hat = dx / (t2 - t1)
    D_hat = u_hat ** 3 * (s2 - s1) / (2.0 * dx)
    return u_hat, D_hat, A1, A2'''))

cells.append(code(r'''X0, X1, SIGMA, TAIL, L_TAIL = 30.0, 60.0, 4.0, 8.0, 80.0   # cm, the bed layout


def simulate_experiment(u, D, dx_probe=100.0, n=1200, cfl=0.5, theta=1.0,
                        limiter=vanleer, sigma=SIGMA, tail=TAIL, strict=True,
                        n_spread=7.0, l_tail=L_TAIL, t_end_mult=1.0, **bed_kw):
    """Run E&R's experiment: one pulse, two detectors a distance dx_probe apart.

    `t_end_mult` scales how long the record runs and `l_tail` is how far the bed
    continues past the second detector. Both are defects waiting to be injected.
    """
    x2 = X1 + dx_probe
    bed = TracerBed(x2 + l_tail, n, u, D, limiter=limiter, **bed_kw)
    c = emg(bed.x_c, X0, sigma, tail)
    i1 = int(np.argmin(np.abs(bed.x_c - X1)))
    i2 = int(np.argmin(np.abs(bed.x_c - x2)))
    dt = cfl * bed.h / u
    spread = np.sqrt(2.0 * D * (bed.x_c[i2] - X0) / u)
    t_end = t_end_mult * (bed.x_c[i2] - X0 + 6.0 * (sigma + tail)
                          + n_spread * spread) / u
    ns = int(np.ceil(t_end / dt))
    t = np.arange(ns + 1) * dt
    c1 = np.empty(ns + 1); c2 = np.empty(ns + 1)
    c1[0], c2[0] = c[i1], c[i2]
    for k in range(1, ns + 1):
        c = bed.step(c, dt, theta=theta, strict=strict)
        c1[k], c2[k] = c[i1], c[i2]
    return bed, t, c1, c2, bed.x_c[i2] - bed.x_c[i1]


# The run E&R print on Figure 13.
DP13 = P["fig13_particle_diameter"]
RE13 = P["fig13_reynolds"]
DX13 = P["fig13_separation"]
EPS13 = float(tab1.loc[np.isclose(tab1.d_p_cm, DP13), "voidage"].iloc[0])
U13 = RE13 * nu_air / (DP13 * EPS13)          # Re = u d_p eps / nu
D13 = D_L_eq17(U13, DP13)

print(f"Figure 13 conditions, as printed: Re = {RE13}, d_p = {DP13} cm, "
      f"Delta x = {DX13:.0f} cm")
print(f"  voidage from Table 1              {EPS13}")
print(f"  interstitial velocity  u        = {U13:8.4f} cm/s"
      f"   (superficial {U13 * EPS13:.4f} cm/s)")
print(f"  Eq. (17)               D_L      = {D13:8.4f} cm2/s")
print(f"  so Pe = u d_p / D_L             = {U13 * DP13 / D13:8.4f}")
print(f"  and Eq. (18) at Re = {RE13}      = {Pe_eq18(RE13):8.4f}")
print(f"  the two agree to {abs(dev(U13 * DP13 / D13, Pe_eq18(RE13))):.2%}, which is"
      " the rounding in Eq. (18)'s")
print("  printed constants and NOT an independent check - see the validation table.")'''))

cells.append(md(r"""## Results

### The correlation against the measurements

Left: the digitised markers of Figure 9 with Eq. (18) through them, and — the
point of the panel — the same model *without* the radial-mixing correction,
which is Eq. (15). Right: the dispersion coefficient itself against velocity for
the five bead sizes, showing the molecular floor they share and the linear
branches that separate them.

The right-hand panel is drawn against *velocity* rather than Reynolds number for
a reason. Eq. (17) depends on the particle size only through the product
$u d_p$, and $\mathrm{Re} = u d_p\varepsilon/\nu$ contains the same product, so
plotted against Re the five curves collapse onto one to within the spread of
their voidages. That is not a defect of the plot; it is why Figure 7 of the
paper carries a single dotted line through five series of markers."""))

cells.append(code(r'''Re_d, Pe_d = fig9.reynolds.to_numpy(), fig9.peclet.to_numpy()
Re_grid = np.logspace(np.log10(0.006), np.log10(80), 400)

# Eq. (15) is Eq. (16) with the correction factor removed, i.e. beta = 0.
Pe_eq15 = 1.0 / (A_printed / Re_grid + C_EDDY)

fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.4))
ax = axes[0]
ax.loglog(Re_d, Pe_d, "o", ms=4.5, mfc="none", color="tab:blue",
          label=f"Fig. 9, {len(Re_d)} digitised markers")
ax.loglog(Re_grid, Pe_eq18(Re_grid), "-", color="k", lw=1.8, label="Eq. (18)")
ax.loglog(Re_grid, Pe_eq15, "--", color="tab:red", lw=1.4,
          label="Eq. (15), no radial-mixing term")
ax.axhline(2.0, color="tab:green", lw=1.0, ls=":")
ax.text(0.012, 2.15, r"$\mathrm{Pe}\to 2$", color="tab:green", fontsize=8)
ax.plot([Re_star], [Pe_star], "*", ms=13, color="tab:orange", zorder=5,
        label=f"Eq. (18) maximum, Pe = {Pe_star:.2f}")
ax.set(xlabel=r"$\mathrm{Re}=u d_p \varepsilon/\nu$",
       ylabel=r"$\mathrm{Pe}=u d_p/D_L$", ylim=(0.015, 8),
       title="Figure 9, and what the correction factor is for")
ax.legend(fontsize=8, loc="upper left")

ax = axes[1]
u_grid = np.logspace(-2.2, 2.4, 300)          # interstitial velocity, cm/s
for _, row in big.iterrows():
    dp = row.d_p_cm
    ax.loglog(u_grid, D_L_eq17(u_grid, dp), lw=1.5, label=f"$d_p$ = {dp} cm")
ax.axhline(gamma * D_M, color="k", lw=1.2, ls=":")
ax.text(0.008, gamma * D_M * 1.2,
        r"$\gamma D_M$ = " + f"{gamma * D_M:.3f}" + r" cm$^2$/s", fontsize=8)
ax.set(xlabel=r"interstitial velocity $u$, cm/s", ylabel=r"$D_L$, cm$^2$/s",
       title="Eq. (17) against velocity: one floor, five branches")
ax.legend(fontsize=8, loc="upper left")
fig.tight_layout()
plt.show()

n_above2 = int((Pe_d > 2.0).sum())
display(Markdown(rf"""
**{n_above2} of the {len(Pe_d)} digitised points lie above $\mathrm{{Pe}} = 2$**,
the largest at {Pe_d.max():.2f}. Eq. (15) — molecular diffusion plus an
uncorrected eddy term — has a supremum of exactly {1 / C_EDDY:.0f} and can
therefore not produce any of them. That is not a matter of fitted constants: the
eddy coefficient is $\tfrac12$ *because* the mixing-cell limit is
$\mathrm{{Pe}} = 2$, so no choice of $\gamma$ moves the ceiling. The maximum is
the evidence for a third mechanism, and it is why the paper has an Eq. (16).
"""))'''))

cells.append(md(r"""### The experiment, simulated

The forward problem at the conditions of their Figure 13: two response curves
100 cm apart, the moments of each, and Eqs. (6) and (7) applied to the
differences. Note how much of the second curve's width is inherited from the
first — that is exactly what the subtraction removes."""))

cells.append(code(r'''bed, t, c1, c2, dx13 = simulate_experiment(U13, D13, dx_probe=DX13, n=1600,
                                           theta=0.5)
u_hat, D_hat, A1, A2 = invert_two_point(t, c1, c2, dx13)
A1m, t1m, s1m = moments(t, c1)
A2m, t2m, s2m = moments(t, c2)

fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.2))
ax = axes[0]
ax.plot(t, c1 / c1.max(), color="tab:blue", lw=1.8,
        label=f"first response curve, x = {X1:.0f} cm")
ax.plot(t, c2 / c2.max(), color="tab:red", lw=1.8,
        label=f"second response curve, x = {X1 + DX13:.0f} cm")
for tt, col in ((t1m, "tab:blue"), (t2m, "tab:red")):
    ax.axvline(tt, color=col, lw=0.9, ls=":")
ax.set(xlabel="time, s", ylabel="concentration (scaled)",
       title="one pulse, two detectors", xlim=(0, t2m + 6 * np.sqrt(s2m)))
ax.legend(fontsize=8)

ax = axes[1]
b2 = TracerBed(bed.L, 1600, U13, D13)
cs = emg(b2.x_c, X0, SIGMA, TAIL)
scale = cs.max()
dt2 = 0.5 * b2.h / U13
tt = 0.0
snaps = [0.0, 0.5, 1.5, 3.0, 4.5]
for target, col in zip(snaps, plt.cm.viridis(np.linspace(0.1, 0.85, len(snaps)))):
    while tt < target - 1e-12:
        cs = b2.step(cs, dt2, theta=0.5)
        tt += dt2
    ax.plot(b2.x_c, cs / scale, color=col, lw=1.5, label=f"t = {target:.1f} s")
for xx in (X1, X1 + DX13):
    ax.axvline(xx, color="k", lw=0.8, ls="--")
ax.set(xlabel="x, cm", ylabel="concentration (scaled)", xlim=(0, bed.L),
       title="the same run in space; dashed lines are the detectors")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

display(Markdown(rf"""
| quantity | first curve | second curve | difference | E&R |
|---|---|---|---|---|
| area $A$ | {A1m:.6g} | {A2m:.6g} | ratio {A2m / A1m:.8f} | Eq. (5), $A_1=A_2$ |
| mean, s | {t1m:.4f} | {t2m:.4f} | {t2m - t1m:.4f} | Eq. (6), $\Delta x/u$ = {dx13 / U13:.4f} |
| variance, s² | {s1m:.5f} | {s2m:.5f} | {s2m - s1m:.5f} | Eq. (7), $2D_L\Delta x/u^3$ = {2 * D13 * dx13 / U13 ** 3:.5f} |

Inverting the two differences with Eqs. (6) and (7) returns $u$ =
{u_hat:.5f} cm/s against the {U13:.5f} cm/s put in
({abs(dev(u_hat, U13)):.1e} relative) and $D_L$ = {D_hat:.5f} cm²/s against
{D13:.5f} ({abs(dev(D_hat, D13)):.1e}). The second curve's variance is
{s2m / s1m:.1f} times the first's, but only {(s2m - s1m) / s2m:.0%} of it is
what the bed between the detectors did; the rest is the injected pulse, and it
cancels.
"""))'''))

cells.append(md(r"""### What the grid measures instead of the bed

That recovery is not free. A first-order upwind convective flux adds a numerical
diffusivity $u h/2$, and a $\theta$-method in time adds $(\theta-\tfrac12)u^2
\Delta t$ — for implicit Euler, $u^2\Delta t/2$. Neither is distinguishable from
physical dispersion by *any* measurement, because they enter the equation in
exactly the same place. So a tracer experiment simulated on a coarse grid
**reports its own truncation error as a dispersion coefficient**.

`A2.1` measures the steady-state version of this ($u\,\mathrm dz/2$ masquerading
as a lower Péclet number). Here it is in the transient inverse problem, and it
is quantitative: the prediction is made before the runs, with no fitted
constant."""))

cells.append(code(r'''rows = []
for n in (400, 800, 1600, 3200):
    for tag, lim, th in (("upwind, implicit Euler", None, 1.0),
                         ("van Leer, implicit Euler", vanleer, 1.0),
                         ("van Leer, Crank-Nicolson", vanleer, 0.5)):
        bd, tt, a, b, dxx = simulate_experiment(U13, D13, dx_probe=DX13, n=n,
                                                cfl=0.5, theta=th, limiter=lim)
        uh, Dh, _, _ = invert_two_point(tt, a, b, dxx)
        dt_used = 0.5 * bd.h / U13
        pred = ((0.0 if lim is not None else U13 * bd.h / 2.0)
                + (th - 0.5) * U13 ** 2 * dt_used)
        rows.append(dict(scheme=tag, n=n, h=bd.h, D_hat=Dh, excess=Dh - D13,
                         predicted=pred))
tab = pd.DataFrame(rows)
tab["ratio"] = np.where(tab.predicted > 0, tab.excess / tab.predicted, np.nan)
tab["rel_err"] = np.abs(tab.D_hat / D13 - 1.0)

for tag, g in tab.groupby("scheme", sort=False):
    print(tag)
    for _, r in g.iterrows():
        extra = (f"predicted {r.predicted:8.5f}  ratio {r.ratio:6.4f}"
                 if r.predicted > 0 else "predicted 0 (both terms switched off)")
        print(f"   n = {int(r.n):5d}  h = {r.h:7.4f}  D_hat = {r.D_hat:9.5f}"
              f"  excess {r.excess:+9.5f}  {extra}")
    print()

cn = tab[tab.scheme == "van Leer, Crank-Nicolson"]
ord_cn = float(np.polyfit(np.log(cn.h.to_numpy()), np.log(cn.rel_err.to_numpy()), 1)[0])
worst_ratio = tab.dropna(subset=["ratio"]).ratio
coarse = tab.query("scheme == 'upwind, implicit Euler' and n == 400").D_hat.iloc[0]
display(Markdown(rf"""
The predicted excess is $u h/2$ for upwind plus $u^2\Delta t/2$ for implicit
Euler, with no fitted constant anywhere, and it accounts for the measured excess
to within **{abs(worst_ratio - 1).max():.1%}** over eight runs and a factor of
eight in grid spacing. Switch both off — van Leer plus Crank–Nicolson — and the
recovered $D_L$ converges on the input at observed order **{ord_cn:.2f}**,
reaching {cn[cn.n == 3200].rel_err.iloc[0]:.1e} relative at $n$ = 3200.

Read the first row the other way round and it is a warning about the
*experiment*, not about the code. A dispersion coefficient inferred from a pulse
is only as good as the resolution of whatever propagates it: at $n$ = 400 the
answer is {coarse / D13:.1f} times the truth, and the response curves it came
from look perfectly smooth.
"""))'''))

cells.append(md(r"""### Where the correlation says dispersion matters

Eq. (17) is usually wanted for one purpose: deciding whether a bed can be
treated as plug flow. The bed Péclet number is $\mathrm{Pe}_L = uL/D_L =
(L/d_p)\,\mathrm{Pe}$, so the correlation and the aspect ratio settle it between
them. What a given $\mathrm{Pe}_L$ then does to conversion is
[`A2.1`](../A2.1-danckwerts-boundary-conditions/)'s subject and is not rebuilt
here."""))

cells.append(code(r'''L_over_dp = np.logspace(0.7, 3.2, 220)
Re_ax = np.logspace(np.log10(P["reynolds_min"]), np.log10(P["reynolds_max"]), 220)
G_R, G_L = np.meshgrid(Re_ax, L_over_dp)
PeL = G_L * Pe_eq18(G_R)

fig, ax = plt.subplots(figsize=(7.0, 4.4))
cf = ax.contourf(G_R, G_L, np.log10(PeL), levels=np.linspace(-1, 4.5, 23),
                 cmap="viridis")
cl = ax.contour(G_R, G_L, PeL, levels=[1, 10, 100, 1000], colors="w", linewidths=1.1)
ax.clabel(cl, fmt="Pe_L = %g", fontsize=8)
ax.set(xscale="log", yscale="log", xlabel=r"$\mathrm{Re}$", ylabel=r"$L/d_p$",
       title="bed Peclet number from Eq. (18), argon in air")
fig.colorbar(cf, ax=ax, label=r"$\log_{10}\mathrm{Pe}_L$")
ax.grid(alpha=0.25, color="w")
fig.tight_layout()
plt.show()

need = 100.0
display(Markdown(rf"""
The map is read along a horizontal line. A bed 100 particles deep runs at
$\mathrm{{Pe}}_L$ = {100 * Pe_eq18(0.01):.1f} at $\mathrm{{Re}}$ = 0.01 but
{100 * Pe_eq18(10.0):.0f} at $\mathrm{{Re}}$ = 10 — the *same bed*, three orders
of magnitude apart in how nearly it is plug flow, because at low Reynolds number
$D_L$ stops falling with velocity and sits on the molecular floor
$\gamma D_M$ = {gamma * D_M:.3f} cm²/s. Reaching $\mathrm{{Pe}}_L$ =
{need:.0f} needs $L/d_p$ = {need / Pe_eq18(0.01):.0f} at $\mathrm{{Re}}$ = 0.01
and only {need / Pe_star:.0f} at the Péclet maximum, $\mathrm{{Re}}$ =
{Re_star:.1f}.

That is the practical content of the molecular-diffusion term, and it is why a
correlation carrying only the eddy branch would mislead badly in exactly the
regime — slow flow, small particles — where laboratory reactors live.
"""))'''))

cells.append(md(r"""## Validation

Five checks, in decreasing order of how much they can fail, and then a table
that breaks each of them on purpose.

**One deviation convention throughout the page**, including the sections above:
$\mathrm{dev} = (\text{model} - \text{measured})/\text{measured}$. At the
scatter of Figure 9 the reciprocal is not interchangeable with it, and the cell
below prints what taking it the other way round would do. The convention is
fixed once, in the helper `dev()`, and used everywhere."""))

cells.append(code(r'''# --- V1. Eq. (18) against the digitised measurements ---------------------
mod = Pe_eq18(Re_d)
d18 = dev(mod, Pe_d)
mad18 = float(np.abs(d18).mean())
bias18 = float(d18.mean())
logrms18 = float(np.std(np.log(mod / Pe_d)))

# nulls, so that "12 %" has something to be compared against
d15 = dev(1.0 / (A_printed / Re_d + C_EDDY), Pe_d)                  # Eq. (15)
Pe_const = float(np.exp(np.mean(np.log(Pe_d))))                     # best constant
d_const = dev(Pe_const, Pe_d)
d_mol = dev(Re_d / A_printed, Pe_d)                                 # molecular branch

print(f"V1  Eq. (18) vs {len(Pe_d)} markers.  deviation = (model - measured)/measured")
print(f"     mean |dev| {mad18:7.2%}   median |dev| {np.median(np.abs(d18)):7.2%}"
      f"   bias {bias18:+7.2%}   worst {np.abs(d18).max():7.2%}")
print(f"     log scatter (1 s.d.) {logrms18:.4f} in ln Pe ="
      f" {np.expm1(logrms18):.2%}")
print(f"\n     null baselines on the same {len(Pe_d)} points:")
print(f"       Eq. (15), no radial-mixing term        mean |dev|"
      f" {np.abs(d15).mean():8.2%}   bias {d15.mean():+7.2%}")
print(f"       molecular branch alone, Pe = Re/{A_printed}   mean |dev|"
      f" {np.abs(d_mol).mean():8.2%}")
print(f"       best single constant, Pe = {Pe_const:.3f}       mean |dev|"
      f" {np.abs(d_const).mean():8.2%}")
rev = np.asarray(Pe_d) / mod - 1.0        # the reciprocal convention
print(f"\n     taken the other way round, (measured - model)/model, the mean")
print(f"     absolute deviation is {np.abs(rev).mean():.2%} - barely different -")
print(f"     but the bias becomes {rev.mean():+.2%} instead of {bias18:+.2%}."
      "  Hence one convention.")

print("\n     Eq.(15) is the informative null: it has the right molecular branch")
print("     and the right high-Re limit, and still misses by a factor - because")
print(f"     its ceiling of {1 / C_EDDY:.0f} sits below {int((Pe_d > 2).sum())} of"
      f" the {len(Pe_d)} points.")'''))

cells.append(code(r'''# --- V2. refit the paper's own constants to the digitised points ---------
#     This checks the DIGITISATION, not the correlation: gamma and beta were
#     fitted by E&R to these very markers, so recovering them from an
#     independent reading of the same figure tests the reading.
def pe_model(Re, g, b):
    return 1.0 / (g * k18 / Re + C_EDDY / (1.0 + b * k18 / Re))


fit_both = least_squares(lambda p: np.log(pe_model(Re_d, *p) / Pe_d), [0.73, 9.7])
fit_beta = least_squares(lambda p: np.log(pe_model(Re_d, gamma, p[0]) / Pe_d), [9.7])

# The points are not independent: the series labels were not extracted, and
# neighbours in Re are often the same bead size. Ordinary standard errors would
# assume an independence the data do not have, so the interval is a moving-block
# bootstrap over contiguous Reynolds windows (the A3.4 lesson).
#
# AN INTERVAL BELONGS TO AN ESTIMATOR, NOT TO A QUANTITY. Two estimators are
# reported below - beta with gamma held at the printed 0.73, and both free - and
# each is bootstrapped THROUGH ITSELF. Quoting one estimator's interval beside
# the other's point estimate is how a page ends up claiming an exclusion it has
# not earned, and the two intervals here differ by about 0.5 in beta.
order = np.argsort(Re_d)
Re_s, Pe_s = Re_d[order], Pe_d[order]
BLOCK, NBOOT = 6, 1500


def block_bootstrap(block, seed, free_gamma, nboot=NBOOT):
    """Refit the SAME estimator on moving-block resamples of the digitised points."""
    r = np.random.default_rng(seed)
    n_blocks = int(np.ceil(len(Re_s) / block))
    out = []
    for _ in range(nboot):
        idx = np.concatenate([np.arange(s, min(s + block, len(Re_s)))
                              for s in r.integers(0, len(Re_s) - block + 1, n_blocks)])
        try:
            if free_gamma:
                f = least_squares(lambda p: np.log(pe_model(Re_s[idx], *p) / Pe_s[idx]),
                                  [0.73, 9.7])
                out.append(f.x)
            else:
                f = least_squares(lambda p: np.log(pe_model(Re_s[idx], gamma, p[0])
                                                   / Pe_s[idx]), [9.7])
                out.append([gamma, f.x[0]])
        except Exception:
            pass
    return np.array(out)


boot_free = block_bootstrap(BLOCK, SEED, True)
boot_fix = block_bootstrap(BLOCK, SEED, False)
gl, gh = np.percentile(boot_free[:, 0], [2.5, 97.5])
bl_free, bh_free = np.percentile(boot_free[:, 1], [2.5, 97.5])
bl_fix, bh_fix = np.percentile(boot_fix[:, 1], [2.5, 97.5])

print("V2  block bootstrap, each estimator resampled through itself")
print(f"    {'estimator':<34s} {'point':>8s} {'95 % interval':>20s}  contains"
      f" the printed {beta}?")
for lbl, est, lo, hi in (
        (f"beta, gamma held at {gamma}", fit_beta.x[0], bl_fix, bh_fix),
        ("beta, both free", fit_both.x[1], bl_free, bh_free)):
    print(f"    {lbl:<34s} {est:8.2f} {f'[{lo:.2f}, {hi:.2f}]':>20s}"
          f"   {'yes' if lo <= beta <= hi else 'NO'}")
print(f"    {'gamma, both free':<34s} {fit_both.x[0]:8.3f}"
      f" {f'[{gl:.3f}, {gh:.3f}]':>20s}"
      f"   printed {gamma}: {'yes' if gl <= gamma <= gh else 'NO'}")
print("\n    block-size dependence of the two beta intervals"
      " (the blocks are not a unit of independence,")
print("     only a dependence-robust device - see the caveat below):")
BLOCKS = (3, 6, 9, 12)
hold = {"gamma fixed": [], "both free": []}
for blk in BLOCKS:
    # BLOCK is already done above; reuse it so the row cannot disagree with the
    # interval quoted in the text.
    bfix = np.percentile((boot_fix if blk == BLOCK
                          else block_bootstrap(blk, SEED, False))[:, 1],
                         [2.5, 97.5])
    bfree = np.percentile((boot_free if blk == BLOCK
                           else block_bootstrap(blk, SEED, True))[:, 1],
                          [2.5, 97.5])
    hold["gamma fixed"].append((blk, bfix[0], bfix[1]))
    hold["both free"].append((blk, bfree[0], bfree[1]))
    print(f"      block {blk:2d}   gamma fixed [{bfix[0]:5.2f}, {bfix[1]:5.2f}]"
          f" {'contains' if bfix[0] <= beta <= bfix[1] else 'EXCLUDES':>8s} {beta}"
          f"    both free [{bfree[0]:5.2f}, {bfree[1]:5.2f}]"
          f" {'contains' if bfree[0] <= beta <= bfree[1] else 'EXCLUDES':>8s} {beta}")
n_hold_fix = sum(lo <= beta <= hi for _, lo, hi in hold["gamma fixed"])
n_hold_free = sum(lo <= beta <= hi for _, lo, hi in hold["both free"])
hi_fix = [hi for _, _, hi in hold["gamma fixed"]]
print(f"    -> the printed {beta} sits inside the gamma-fixed interval at"
      f" {n_hold_fix} of the {len(BLOCKS)} block sizes and inside the both-free"
      f" one at {n_hold_free};")
print(f"       the gamma-fixed upper limit runs from {min(hi_fix):.2f} to"
      f" {max(hi_fix):.2f}, i.e. it straddles {beta}. Whether {beta} is 'excluded'"
      " is decided")
print("       in the second decimal by a tuning constant, which is not a"
      " finding about the paper.")

display(Markdown(rf"""
**V2. The figure and the page image are independent witnesses, and they agree.**
$\gamma$ and $\beta$ were fitted by Edwards and Richardson to the markers of
Figure 9, so refitting them to a fresh digitisation of that figure does not test
the correlation — it tests the extraction and the transcription of Eq. (17)
against each other. Holding $\gamma$ at the printed {gamma}, least squares in
$\ln\mathrm{{Pe}}$ over the {len(Re_d)} points returns $\beta$ =
**{fit_beta.x[0]:.2f}** against the printed **{beta}**
({abs(dev(fit_beta.x[0], beta)):.1%}), with a 95 % block-bootstrap interval
**[{bl_fix:.2f}, {bh_fix:.2f}]**. Fitting both parameters is a *different*
estimator and returns $\gamma$ = {fit_both.x[0]:.3f}, $\beta$ =
{fit_both.x[1]:.2f}, with intervals [{gl:.3f}, {gh:.3f}] and
[{bl_free:.2f}, {bh_free:.2f}].

**Each interval is quoted beside the estimate it belongs to, and no exclusion is
claimed.** The printed {beta} sits inside the $\gamma$-fixed interval — the one
that matches the reported {fit_beta.x[0]:.2f} — at {n_hold_fix} of the
{len(BLOCKS)} block sizes tried, and inside the both-free interval at
{n_hold_free}. It is not comfortably inside either: the $\gamma$-fixed upper
limit runs from {min(hi_fix):.2f} to {max(hi_fix):.2f} across those block sizes,
so it straddles {beta}, and the both-free upper limit sits a little lower still.
**Whether {beta} falls inside or outside is therefore decided in the second
decimal by the block size**, which is a tuning constant of the interval and not a
fact about the paper. The printed {gamma} lies inside the both-free interval for
$\gamma$, [{gl:.3f}, {gh:.3f}]. Nothing here supports a claim that the digitisation
disagrees with the printed constants, and a uniform $\pm$5 px shift of the
extraction — one of the two systematics V3 bounds — moves $\beta$ by more than
the gap in either direction anyway.

**What the blocks are and are not.** Resampling contiguous blocks of {BLOCK}
points is a dependence-robust device, not a claim that a block is one
experiment: the five bead sizes interleave in Reynolds number above
$\mathrm{{Re}} \approx 1$, so a contiguous block there mixes series. The
interval is wider than an independence-assuming standard error, which is the
point; how much wider is not something this figure can settle.
"""))'''))

cells.append(code(r'''# --- V3. how much of V1 is the digitisation? ----------------------------
#     Two perturbations, both bounded by the sidecar: the worst ordinate tick
#     residual (3.8 px) and the ink-centroid-versus-bounding-box ambiguity of a
#     triangular glyph (about 5 px). The ordinate runs at 604.56 px per decade.
PX_PER_DECADE = 604.56
print("V3  sensitivity of the V1 headline to the two digitisation systematics")
print(f"     as extracted                          mean |dev| {mad18:7.2%}"
      f"   bias {bias18:+7.2%}")
v3 = []
for shift_px, why in ((+3.8, "worst ordinate tick residual"),
                      (-3.8, "worst ordinate tick residual"),
                      (+5.0, "triangle ink centroid vs box centre"),
                      (-5.0, "triangle ink centroid vs box centre")):
    Pe_shift = Pe_d * 10 ** (-shift_px / PX_PER_DECADE)
    ds = dev(Pe_eq18(Re_d), Pe_shift)
    v3.append((np.abs(ds).mean(), ds.mean()))
    print(f"     every point moved {shift_px:+5.1f} px in Pe      "
          f"mean |dev| {np.abs(ds).mean():7.2%}   bias {ds.mean():+7.2%}"
          f"   ({why})")
print("     -> the SCATTER is the measurement's, not the extraction's; the BIAS")
print("        is not resolved better than about +/- 2 %, so it is not claimed.")'''))

cells.append(code(r'''# --- V4. a second model, and an independent closed form for it -----------
class TwoZoneBed:
    """Two mobile zones exchanging laterally. Layout (n, 2).

        th_i dc_i/dt + d(th_i u_i c_i)/dx = d(th_i D0 dc_i/dx)/dx -/+ kappa (c_1 - c_2)
    """

    def __init__(self, L, n, th1, u1, u2, D0, kappa, limiter=vanleer,
                 tol=1e-9, max_corr=40, kappa_scale=1.0, mass="theta"):
        # kappa_scale and mass exist only so that V4b can break this class on
        # purpose. The defaults are the model.
        kappa = kappa * kappa_scale
        th = np.array([th1, 1.0 - th1]); uu = np.array([u1, u2])
        self.n, self.th, self.uu, self.h = n, th, uu, L / n
        self.shape = (n, 2)                        # spatial axis first, fields last
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.bc = ({"a": 0.0, "b": 1.0, "d": 0.0},      # c = 0 at the inlet
                   {"a": 1.0, "b": 0.0, "d": 0.0})      # dc/dx = 0 at the outlet
        self.v = np.broadcast_to(th * uu, (n + 1, 2))
        conv, conv_bc = construct_convflux_upwind(self.shape, self.x_f, self.x_c,
                                                  self.bc, v=self.v, axis=0)
        grad, grad_bc = construct_grad(self.shape, self.x_f, self.x_c, self.bc,
                                       axis=0)
        div = construct_div(self.shape, self.x_f, nu=0, axis=0)   # nu=0: Cartesian
        Df = diags_array(np.broadcast_to(th * D0, (n + 1, 2)).ravel())
        i = np.arange(n)
        rows_ = np.stack([2 * i, 2 * i, 2 * i + 1, 2 * i + 1]).T.ravel()
        cols_ = np.stack([2 * i, 2 * i + 1, 2 * i, 2 * i + 1]).T.ravel()
        vals_ = np.tile([kappa, -kappa, -kappa, kappa], n)
        E = coo_array((vals_, (rows_, cols_)), shape=(2 * n, 2 * n)).tocsc()
        self.div = div
        self.A = (div @ (conv - Df @ grad) + E).tocsc()
        self.b0 = np.asarray((div @ (conv_bc - Df @ grad_bc)).todense()).ravel()
        # accumulation carries th_i: the volume fraction multiplies dc_i/dt
        self.M = (diags_array(np.tile(th, n)) if mass == "theta"
                  else eye_array(2 * n, format="csc"))
        self.limiter, self.tol, self.max_corr = limiter, tol, max_corr
        self._key, self._last_corr = None, np.zeros(2 * n)

    def step(self, c, dt, theta=0.5):
        if self._key != (dt, theta):
            self._lu = splu((self.M / dt + theta * self.A).tocsc())
            self._key = (dt, theta)
        lu = self._lu
        rhs = self.M @ c / dt - (1.0 - theta) * (self.A @ c) - self.b0
        c_new = lu.solve(rhs - self._last_corr)
        done, it, corr = False, 0, self._last_corr
        for it in range(1, self.max_corr + 1):
            _, dc_f = interp_cntr_to_stagg_tvd(c_new.reshape(self.shape), self.x_f,
                                               self.x_c, self.bc, self.v,
                                               tvd_limiter=self.limiter, axis=0)
            corr = np.asarray(self.div @ (self.v * dc_f).reshape(-1, 1)).ravel()
            cc = lu.solve(rhs - corr)
            done = np.max(np.abs(cc - c_new)) <= self.tol * max(1.0, np.max(np.abs(cc)))
            c_new = cc
            if done:
                break
        assert done, f"two-zone deferred correction did not converge in {self.max_corr}"
        CORR_ITERS.append(it)
        self._last_corr = corr
        return c_new


def D_two_zone(th1, u1, u2, kappa):
    """Added dispersion of the two-zone model, from expanding c_i = c + a_i with
    th_1 a_1 + th_2 a_2 = 0:   D_add = th1^2 th2^2 (u1-u2)^2 / kappa.
    This shares no code with the solve - it comes from expanding the system, not
    from assembling it."""
    th2 = 1.0 - th1
    return (th1 * th2) ** 2 * (u1 - u2) ** 2 / kappa


def simulate_two_zone(th1, u1, u2, D0, kappa, dx_probe=100.0, n=1200, cfl=0.5,
                      sigma=SIGMA, tail=TAIL, **bed_kw):
    u = th1 * u1 + (1.0 - th1) * u2
    x2 = X1 + dx_probe
    bed = TwoZoneBed(x2 + L_TAIL, n, th1, u1, u2, D0, kappa, **bed_kw)
    c0 = emg(bed.x_c, X0, sigma, tail)
    c = np.column_stack([c0, c0]).ravel()
    i1 = int(np.argmin(np.abs(bed.x_c - X1)))
    i2 = int(np.argmin(np.abs(bed.x_c - x2)))
    dt = cfl * bed.h / max(u1, u2)
    Dg = D0 + D_two_zone(th1, u1, u2, kappa)
    spread = np.sqrt(2.0 * Dg * (bed.x_c[i2] - X0) / u)
    t_end = (bed.x_c[i2] - X0 + 6.0 * (sigma + tail) + 8.0 * spread) / u
    ns = int(np.ceil(t_end / dt))
    t = np.arange(ns + 1) * dt
    w = bed.th * bed.uu / u                     # the detector sees the flux average
    c1 = np.empty(ns + 1); c2 = np.empty(ns + 1)
    cr = c.reshape(n, 2); c1[0], c2[0] = cr[i1] @ w, cr[i2] @ w
    for k in range(1, ns + 1):
        c = bed.step(c, dt)
        cr = c.reshape(n, 2)
        c1[k], c2[k] = cr[i1] @ w, cr[i2] @ w
    return bed, t, c1, c2, bed.x_c[i2] - bed.x_c[i1], Dg


# fine-particle conditions, from the Figure 15 legend
DP15 = P["fig15_particle_diameter"]
RE15 = P["fig15_reynolds"]
DL15 = P["fig15_apparent_dispersion_coefficient"]
EPS15 = float(tab1.loc[np.isclose(tab1.d_p_cm, DP15), "voidage"].iloc[0])
UMF15 = float(tab1.loc[np.isclose(tab1.d_p_cm, DP15), "u_mf_cm_s"].iloc[0])
U15 = RE15 * nu_air / (DP15 * EPS15)
D15 = D_L_eq17(U15, DP15)

TH1, RATIO = 0.15, 2.0
u1_ref = RATIO * U15
u2_ref = (U15 - TH1 * u1_ref) / (1.0 - TH1)
kappa_ref = D_two_zone(TH1, u1_ref, u2_ref, 1.0) / (DL15 - D15)

print("V4  two-zone model: pymrm solve against an independent closed form")
prev, two_zone_rel = None, None
for n in (500, 1000, 2000):
    bd, tt, a, b, dxx, Dg = simulate_two_zone(TH1, u1_ref, u2_ref, D15,
                                              kappa_ref, n=n)
    uh, Dh, _, _ = invert_two_point(tt, a, b, dxx)
    two_zone_rel = abs(dev(Dh, Dg))
    rate = "" if prev is None else f"   ratio {prev / two_zone_rel:6.2f}"
    print(f"     n = {n:5d}  h = {bd.h:6.4f}  D_hat = {Dh:.6f}   "
          f"closed form {Dg:.6f}   rel {two_zone_rel:.2e}{rate}")
    prev = two_zone_rel

# V4b. A residual of 1e-04 is worth nothing until something is broken against
# it. Two injections into TwoZoneBed's assembly, both of which the closed form
# D_add = th1^2 th2^2 (u1-u2)^2 / kappa cannot know about.
print("\nV4b what that residual can see: the same comparison with the assembly"
      " broken on purpose")
print(f"     {'injected defect':<48s} {'D_hat':>10s} {'closed form':>12s} {'rel':>10s}")
v4_breaks = []
for lbl, kw in (("none - the model as published", {}),
                ("exchange coefficient doubled in the solve only",
                 dict(kappa_scale=2.0)),
                ("theta-weighted mass matrix replaced by the identity",
                 dict(mass="identity"))):
    bd, tt, a, b, dxx, Dg = simulate_two_zone(TH1, u1_ref, u2_ref, D15,
                                              kappa_ref, n=500, **kw)
    uh, Dh, _, _ = invert_two_point(tt, a, b, dxx)
    v4_breaks.append((lbl, Dh, Dg, abs(dev(Dh, Dg))))
    print(f"     {lbl:<48s} {Dh:10.4f} {Dg:12.6f} {abs(dev(Dh, Dg)):10.2e}")
v4_break_worst = max(r[3] for r in v4_breaks[1:])
v4_factors = [r[3] / v4_breaks[0][3] for r in v4_breaks[1:]]
print(f"     -> V4 is not decoration: the residual grows by factors of"
      f" {v4_factors[0]:.0f} and {v4_factors[1]:.0f}, and the")
print("        identity mass matrix drives the recovered coefficient"
      f" {'negative' if v4_breaks[2][1] < 0 else 'far from the closed form'}"
      f" ({v4_breaks[2][1]:.3f}).")'''))

cells.append(code(r'''# --- V5. deliberate defects. Every check above is broken on purpose and the
#         number reported. A check whose number does not move is decoration.
BASE = dict(u=U13, D=D13, dx_probe=DX13, n=800, cfl=0.5, theta=0.5)


def recovered(**kw):
    kw = {**BASE, **kw}
    try:
        bd, tt, a, b, dxx = simulate_experiment(**kw)
        uh, Dh, a1, a2 = invert_two_point(tt, a, b, dxx)
        if not np.isfinite(Dh):
            return np.nan, np.nan, np.nan
        return Dh, uh, a2 / a1
    except (AssertionError, ValueError, FloatingPointError):
        return np.nan, np.nan, np.nan


breaks = [
    ("none - the model as published", recovered()),
    ("sign of the dispersive flux flipped (D_sign = -1)",
     recovered(D_sign=-1.0)),
    ("nu = 1 (cylindrical) instead of nu = 0", recovered(nu=1)),
    ("Dirichlet outlet, bed 80 cm past the second detector",
     recovered(outlet="dirichlet")),
    ("Dirichlet outlet, bed 0.3 cm past the second detector",
     recovered(outlet="dirichlet", l_tail=0.3)),
    ("record cut off at 65 % of its length",
     recovered(t_end_mult=0.65)),
    ("bare upwind, n = 800 (numerical dispersion left in)",
     recovered(limiter=None)),
    ("bare upwind, n = 200", recovered(limiter=None, n=200)),
]
print("V5  defect sensitivity of the recovered dispersion coefficient")
print(f"    reference: D_L = {D13:.5f} cm2/s, u = {U13:.5f} cm/s")
print(f"    'nan' means the run refused to finish, which is also a signal\n")
print(f"    {'injected defect':<54s} {'D_hat':>9s} {'D_hat/D_L':>10s}"
      f" {'u_hat/u':>10s} {'A2/A1':>10s}")
for name, (Dh, uh, ar) in breaks:
    print(f"    {name:<54s} {Dh:9.4f} {Dh / D13:10.4f} {uh / U13:10.6f}"
          f" {ar:10.6f}")

# V5c. A defect NEITHER diagnostic catches, and the reason is structural: they
# test the solver, and this one is upstream of it. Table 1 has eight columns and
# the 0.203 cm row sits directly beside the 0.300 cm row used here, so reading
# the voidage off the wrong line is the likeliest transcription slip on the page.
print("\nV5c a defect the table above cannot see: the wrong row of Table 1")
print(f"    {'voidage used':<34s} {'u, cm/s':>9s} {'D_L Eq.(17)':>12s}"
      f" {'D_hat/D_L':>10s} {'u_hat/u':>10s} {'A2/A1':>10s}")
chain = []
for lbl, eps_used in ((f"{EPS13} - the 0.300 cm row (right)", EPS13),
                      (f"{float(tab1.loc[np.isclose(tab1.d_p_cm, 0.203), 'voidage'].iloc[0])}"
                       " - the 0.203 cm row beside it", float(
                           tab1.loc[np.isclose(tab1.d_p_cm, 0.203), "voidage"].iloc[0])),
                      (f"{eps18} - Eq. (18)'s rounded value", eps18)):
    u_e = RE13 * nu_air / (DP13 * eps_used)
    D_e = D_L_eq17(u_e, DP13)
    Dh, uh, arr = recovered(u=u_e, D=D_e)
    chain.append((eps_used, u_e, D_e, Dh / D_e, uh / u_e, arr))
    print(f"    {lbl:<34s} {u_e:9.4f} {D_e:12.4f} {Dh / D_e:10.5f}"
          f" {uh / u_e:10.6f} {arr:10.6f}")
chain_u_dev = abs(dev(chain[1][1], chain[0][1]))
chain_D_dev = abs(dev(chain[1][2], chain[0][2]))
print(f"    -> u moves {chain_u_dev:.1%} and D_L moves {chain_D_dev:.1%},"
      " and all three diagnostics are")
print("       identical to five decimals. They test the solver, not the"
      " parameter chain that feeds it.")
'''))

cells.append(code(r'''# and the transcription checks, broken by one digit each
print("V5b the transcription checks, one mis-read digit at a time")
print(f"    {'':<46s} {'derived':>9s} {'printed':>9s} {'|dev|':>8s}")
for name, g, b in (("as printed: gamma = 0.73, beta = 9.7", gamma, beta),
                   ("gamma mis-read as 0.63", 0.63, beta),
                   ("gamma mis-read as 0.78", 0.78, beta),
                   ("beta mis-read as 9.1", gamma, 9.1),
                   ("beta mis-read as 8.7", gamma, 8.7)):
    print(f"    Eq.(18) 1st coeff  {name:<28s} {g * k18:9.4f} {A_printed:9.4f}"
          f" {abs(dev(g * k18, A_printed)):8.2%}")
    print(f"    Eq.(18) 2nd coeff  {'':<28s} {b * k18:9.4f} {B_printed:9.4f}"
          f" {abs(dev(b * k18, B_printed)):8.2%}")
print("\n    the maximum of Eq. (18), which the paper bounds at 'greater than 3'")
for lbl, a_c, b_c in (("as printed", A_printed, B_printed),
                      ("1st coeff 0.30", 0.30, B_printed),
                      ("1st coeff 0.48", 0.48, B_printed),
                      ("2nd coeff 3.0", A_printed, 3.0),
                      ("2nd coeff 8.0", A_printed, 8.0)):
    rs, ps = Pe_max_eq18(a_c, b_c)
    print(f"    {lbl:<18s} Pe_max = {ps:6.3f} at Re = {rs:7.3f}"
          f"   consistent with '> 3': {'yes' if ps > 3 else 'NO'}")
print("\n    and the gas crossover Reynolds number, printed as 1.8")
for g in (0.63, 0.68, gamma, 0.78):
    rc = Re_crossover(P["typical_diffusivity_gas"], nu_air, gamma=g)
    print(f"    gamma = {g:4.2f}  ->  Re = {rc:5.3f}"
          f"   ({abs(dev(rc, P['crossover_reynolds_gas'])):6.2%} from printed)")'''))

cells.append(code(r'''ar = {n: v[2] for n, v in breaks}
dd = {n: v[0] / D13 for n, v in breaks}
uu_ = {n: v[1] / U13 for n, v in breaks}
NU1 = "nu = 1 (cylindrical) instead of nu = 0"
DIR_FAR = "Dirichlet outlet, bed 80 cm past the second detector"
DIR_NEAR = "Dirichlet outlet, bed 0.3 cm past the second detector"
CUT = "record cut off at 65 % of its length"
UPW = "bare upwind, n = 200"
display(Markdown(rf"""
### What each check can and cannot see, measured

| check | reference value | what moved it | what did not |
|---|---|---|---|
| $A_2/A_1$ | {ar["none - the model as published"]:.6f} | the wrong geometry ({ar[NU1]:.3f}) and a detector 0.3 cm from an absorbing outlet ({ar[DIR_NEAR]:.3f}) | the numerical dispersion that puts $\hat D_L$ out by a factor {dd[UPW]:.1f} — it stays at {ar[UPW]:.6f} |
| $\hat u/u$ | {uu_["none - the model as published"]:.6f} | almost nothing; the largest excursion in the whole table is {max(abs(v - 1) for v in uu_.values() if np.isfinite(v)):.1e} | every defect that matters. $\Delta\bar t = \Delta x/u$ holds in the discrete system too for a uniform velocity |
| $\hat D_L/D_L$ | {dd["none - the model as published"]:.5f} | numerical dispersion ({dd[UPW]:.3f} at $n$ = 200), a truncated record ({dd[CUT]:.3f}), a flipped dispersion sign (refuses to converge) | the wrong geometry ({dd[NU1]:.4f}) and the outlet ({dd[DIR_FAR]:.4f}) |
| Eq. (18) from Eq. (17) | {abs(dev(A_derived, A_printed)):.2%} / {abs(dev(B_derived, B_printed)):.2%} | one mis-read digit in $\gamma$ or $\beta$: {abs(dev(0.63 * k18, A_printed)):.1%} and {abs(dev(8.7 * k18, B_printed)):.1%} above | whether the printed constants are themselves right |
| Eq. (18) vs Figure 9 | {mad18:.1%} | removing the radial-mixing term ({np.abs(d15).mean():.0%}) | the constants, which were fitted to these very points |
| two-zone solve vs $D_\mathrm{{add}}$ closed form | {v4_breaks[0][3]:.1e} | doubling $\kappa$ in the solve ({v4_breaks[1][3]:.1e}) and replacing the $\theta$-weighted mass matrix with the identity ({v4_breaks[2][3]:.1e}, recovered $\hat D_L$ = {v4_breaks[2][1]:.3f}) | anything upstream of the assembly — it is handed $\theta_1$, the velocities and $\kappa$ |
| the parameter chain | — | **nothing in this table** | reading the voidage off the adjacent row of Table 1 moves $u$ by {chain_u_dev:.1%} and $D_L$ by {chain_D_dev:.1%} while $\hat D_L/D_L$, $\hat u/u$ and $A_2/A_1$ stay at {chain[1][3]:.5f}, {chain[1][4]:.6f}, {chain[1][5]:.6f} |

**Three things this table says that the numbers alone would not.**

$A_2/A_1$ and $\hat D_L/D_L$ are **complementary, not redundant**: each is blind
to exactly what the other catches. A wrong geometry leaves the recovered
dispersion coefficient right to {abs(dd[NU1] - 1):.2%} and destroys the area
ratio; the coarsest upwind grid does the reverse. Reporting either alone would
have hidden a defect the other catches.

**They are not jointly sufficient, and the last row is the proof.** Every
diagnostic on this page is computed *inside* the simulation, so all of them are
blind to what goes into it. Substituting the voidage of the 0.203 cm beads —
the row printed directly beside the 0.300 cm row this section uses, in a table
with eight columns — changes every physical number in the Figure-13 section by
{max(chain_u_dev, chain_D_dev):.1%} and leaves all three diagnostics identical
to five decimals. The pair tests the solver; it does not test the parameter
chain, and nothing here does. What guards that chain is elsewhere and is weaker:
V5b breaks the transcription of $\gamma$ and $\beta$ one digit at a time, the
Table 1 voidages are checked against the printed bulk and material densities in
*The data*, and the row itself is selected by matching $d_p$ rather than by
copying a number across. That last one is the only real defence against this
particular slip, and it is a coding convention, not a measurement.

$\hat u/u$ is **near-structural and is labelled as such**. It never leaves 1 by
more than {max(abs(v - 1) for v in uu_.values() if np.isfinite(v)):.1e} anywhere
in the table, including runs whose $D_L$ is out by a factor of five. It confirms
the bookkeeping and nothing else.

And the outlet row is a **blind spot stated as one**: with the bed continuing
80 cm past the second detector, changing the outlet condition from zero-gradient
to Dirichlet moves nothing at all — the two runs return the same $\hat D_L$ to
{max(1e-16, abs(dd[DIR_FAR] / dd["none - the model as published"] - 1)):.0e}
relative, which is to say they are the same numbers. Nothing on this page tests
the outlet condition, because in this geometry nothing can.
[`A2.1`](../A2.1-danckwerts-boundary-conditions/) is where that question lives.

Note finally what the {abs(dev(U13 * DP13 / D13, Pe_eq18(RE13))):.2%} agreement
between Eq. (17) and Eq. (18) at $\mathrm{{Re}}$ = {RE13} is *not*: those two
routes share $\gamma$, $\beta$, $\varepsilon$ and Sc, so all it measures is
the rounding of Eq. (18)'s printed constants to two significant figures. It is
on the page as arithmetic, not as evidence.

**Every deferred correction asserted its own convergence**: {len(CORR_ITERS)}
time steps across the whole page, median {int(np.median(CORR_ITERS))} iterations,
{np.mean(np.asarray(CORR_ITERS) <= 5):.0%} of them at 5 or fewer, and the worst
needing {max(CORR_ITERS)} of a cap of 40. The one defect that cannot converge —
the flipped dispersion sign, which makes the equation backwards-parabolic —
announces itself by failing that assertion rather than by returning a number.
"""))'''))

cells.append(code(r'''report_agreement("A2.5", {
    "eq18_vs_fig9_mean_abs_dev": mad18,
    "eq18_vs_fig9_bias": abs(bias18),
    "eq18_vs_fig9_log_scatter": logrms18,
    "eq15_vs_fig9_mean_abs_dev": float(np.abs(d15).mean()),
    "beta_refit_rel_dev": float(abs(dev(fit_beta.x[0], beta))),
    "eq18_coeff1_from_eq17_rel_dev": float(abs(dev(A_derived, A_printed))),
    "eq18_coeff2_from_eq17_rel_dev": float(abs(dev(B_derived, B_printed))),
    "gas_crossover_Re_rel_dev": float(abs(dev(Re_gas, P["crossover_reynolds_gas"]))),
    "gamma_from_low_Re_plateau_rel_dev": float(abs(dev(gamma_ratio, gamma))),
    "moment_recovery_rel_err_n3200": float(cn[cn.n == 3200].rel_err.iloc[0]),
    "numerical_dispersion_prediction_worst": float(abs(worst_ratio - 1).max()),
    "two_zone_vs_closed_form_rel_err": float(two_zone_rel),
    # the break tests, so that a check losing its power is itself a regression
    "two_zone_break_smallest_response": float(min(r[3] for r in v4_breaks[1:])),
    "parameter_chain_D_L_shift_unseen": float(chain_D_dev),
    "beta_ci_gamma_fixed_lo": float(bl_fix),
    "beta_ci_gamma_fixed_hi": float(bh_fix),
    "n_fig9_markers": float(len(Pe_d)),
})'''))

cells.append(md(r"""## What pymrm adds

### The fine-particle anomaly, and what their own diagnostic says about it

With the 0.0097 cm beads, the sand and the diakon, Edwards and Richardson's
moment analysis returned dispersion coefficients far above Eq. (17) — Figure 15
prints $D_L$ = 0.503 cm²/s at $\mathrm{Re}$ = 0.0710, where the correlation
gives about 0.15. They then applied their own consistency test, computing the
second response curve from the first with that $D_L$, and found the *measured*
curve peaked **earlier** than the computed one. They attributed it to
channelling: the materials "have a fair spread of particle sizes … it is thought
that channels are being formed, the fine particles being swept aside laterally".

The moment method cannot itself tell a genuinely larger $D_L$ from a velocity
distribution, because it only ever sees two moments. The shape can. So the
question a solver can answer is: **what velocity structure reproduces both their
$D_L$ and the sign of their shape discrepancy?**"""))

cells.append(code(r'''print(f"Figure 15 conditions: d_p = {DP15} cm, Re = {RE15}, "
      f"apparent D_L = {DL15} cm2/s")
print(f"  voidage (Table 1)                       {EPS15}")
print(f"  interstitial velocity                   {U15:.4f} cm/s")
print(f"  superficial velocity                    {U15 * EPS15:.4f} cm/s")
print(f"  minimum fluidisation velocity (Table 1) {UMF15} cm/s"
      f"   -> the run sits at {U15 * EPS15 / UMF15:.0%} of it")
print(f"  Eq. (17) prediction                     {D15:.5f} cm2/s")
print(f"  measured / predicted                    {DL15 / D15:.2f}"
      f"   (deviation {dev(D15, DL15):+.1%})")'''))

cells.append(code(r'''def inverse_gaussian_kernel(t, dx, u, D):
    """Impulse response of Eq. (1) over a distance dx in an infinite bed - E&R's
    Appendix B calculation in closed form. Its area is 1, its mean dx/u and its
    variance 2 D dx/u^3, i.e. exactly their Eqs. (5)-(7)."""
    g = np.zeros_like(t)
    m = t > 0
    g[m] = dx / np.sqrt(4 * np.pi * D * t[m] ** 3) * \
        np.exp(-(dx - u * t[m]) ** 2 / (4 * D * t[m]))
    return g


def peak_time(t, c):
    i = int(np.argmax(c))
    y0, y1, y2 = c[i - 1], c[i], c[i + 1]
    return t[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * (t[1] - t[0])


# DX_ASSUMED. Figure 15's legend prints the bead size, the Reynolds number and
# the apparent D_L, and NOT the detector separation. E&R p. 111 say the test
# section was "either of two tubes, one 21.3 cm and the other 100 cm long", and
# their conclusions quote 21.3 < Delta x < 121.3 cm. Everything below therefore
# ASSUMES the 100 cm section, which is the one Figure 13 prints; the alternative
# is run at the end of this section and it is not a small difference.
DX_ASSUMED = DX13


def channelling_case(flux_frac, th1=0.15, sigma=SIGMA, tail=TAIL, n=1000,
                     dx_probe=None):
    """A two-zone bed tuned so that the two-point moment method returns E&R's
    measured apparent D_L, with a given share of the total FLOW in the minority
    zone. Returns how far the computed second curve lags the simulated one, and
    the peak-to-peak of the shape residual between them."""
    dx_probe = DX_ASSUMED if dx_probe is None else dx_probe
    u1 = flux_frac * U15 / th1
    u2 = (1.0 - flux_frac) * U15 / (1.0 - th1)
    kappa = D_two_zone(th1, u1, u2, 1.0) / (DL15 - D15)
    bd, t, a, b, dxx, Dg = simulate_two_zone(th1, u1, u2, D15, kappa, n=n,
                                             sigma=sigma, tail=tail,
                                             dx_probe=dx_probe)
    uh, Dh, _, _ = invert_two_point(t, a, b, dxx)
    g = inverse_gaussian_kernel(t, dxx, uh, Dh)
    g = g / np.trapezoid(g, t)
    comp = np.convolve(a, g)[:len(t)] * (t[1] - t[0])
    p_sim, p_comp = peak_time(t, b), peak_time(t, comp)
    resid = float(np.ptp(comp / comp.max() - b / b.max()))
    return dict(flux_frac=flux_frac, th1=th1, u1=u1, u2=u2, kappa=kappa,
                D_hat=Dh, lag=(p_comp - p_sim) / p_sim, resid=resid,
                dx=dx_probe, t=t, sim=b, comp=comp, first=a)


scan = pd.DataFrame([{k: v for k, v in channelling_case(ff, n=800).items()
                      if k not in ("t", "sim", "comp", "first")}
                     for ff in (0.01, 0.03, 0.05, 0.08, 0.25, 0.40, 0.60, 0.80)])
scan["u1_over_u"] = scan.u1 / U15
scan["u2_over_u"] = scan.u2 / U15
print("The minority zone holds 15 % of the void. kappa is set in every row so")
print("that the two-point moment method returns E&R's measured 0.503 cm2/s.\n")
print(scan[["flux_frac", "u1_over_u", "u2_over_u", "kappa", "D_hat", "lag",
            "resid"]]
      .to_string(index=False,
                 formatters={"flux_frac": "{:.2f}".format,
                             "u1_over_u": "{:.3f}".format,
                             "u2_over_u": "{:.3f}".format,
                             "kappa": "{:.3f}".format,
                             "D_hat": "{:.4f}".format,
                             "lag": "{:+.2%}".format,
                             "resid": "{:.1%}".format}))'''))

cells.append(code(r'''# The scan above holds the minority VOLUME at 15 % and moves the flow split.
# That is one slice of the parameter space, and the peak displacement it bounds
# is a maximum over that slice only. Here the volume fraction moves too, with
# the velocity RATIO held fixed, so that every row is still tuned to the same
# measured D_L. Both directions are run: the minority zone slow, and fast.
def two_zone_split(th1, ratio, n=500, dx_probe=None):
    """Volume fraction th1 in zone 1, velocity ratio u1/u2 = ratio, kappa tuned
    so the moment method still returns E&R's 0.503 cm2/s."""
    u2 = U15 / (th1 * ratio + (1.0 - th1))
    u1 = ratio * u2
    ff = th1 * u1 / U15                       # share of the total flow in zone 1
    return channelling_case(ff, th1=th1, n=n, dx_probe=dx_probe)


SLOW, FAST = 0.2, 5.0                          # u1/u2 for the two directions
# th1 = 0.5 is deliberately not scanned: with equal volumes there is no minority
# zone, so the sign statement below has nothing to say about it.
th_scan = []
for th1 in (0.03, 0.05, 0.15, 0.35, 0.65, 0.85):
    for ratio in (SLOW, FAST):
        r = two_zone_split(th1, ratio)
        minority_slow = (th1 < 0.5) == (ratio < 1.0)
        th_scan.append(dict(th1=th1, ratio=ratio, minority_slow=minority_slow,
                            D_hat=r["D_hat"], lag=r["lag"], resid=r["resid"]))
th_scan = pd.DataFrame(th_scan)
th_scan["sign_ok"] = np.sign(th_scan.lag) == np.where(th_scan.minority_slow, 1, -1)
print("Volume split scanned at fixed velocity ratio; kappa retuned every row.")
print("'minority_slow' is whether the zone holding LESS of the void is the"
      " slower one;")
print("'sign_ok' is whether the computed peak then lags, as E&R report.\n")
print(th_scan.to_string(index=False,
                        formatters={"th1": "{:.2f}".format,
                                    "ratio": "{:.1f}".format,
                                    "D_hat": "{:.4f}".format,
                                    "lag": "{:+.3%}".format,
                                    "resid": "{:.1%}".format}))
sign_ok = int(th_scan.sign_ok.sum())
flux_ok = int((np.sign(scan.lag) == np.where(scan.u1 < scan.u2, 1, -1)).sum())
bad = th_scan[~th_scan.sign_ok]
LAG_FLOOR = 0.0005                    # 0.05 % of transit time
big = th_scan[th_scan.lag.abs() > LAG_FLOOR]
print(f"\nsign follows minority-by-volume-slow in {sign_ok} of {len(th_scan)}"
      f" rows here and {flux_ok} of {len(scan)} in the flux scan above.")
if len(bad):
    print(f"the {len(bad)} exceptions all sit at th1 = "
          + ", ".join(f"{v:.2f}" for v in bad.th1)
          + f", and every one has |lag| <= {bad.lag.abs().max():.3%} -")
    print("   near an even volume split the displacement is not just small but"
          " signless, so the")
    print("   rule holds where there is anything to have a sign:"
          f" {int(big.sign_ok.sum())} of {len(big)} rows with"
          f" |lag| > {LAG_FLOOR:.2%}.")
print(f"largest peak lag anywhere: {th_scan.lag.abs().max():.2%} of the transit"
      f" time, at th1 = {th_scan.loc[th_scan.lag.abs().idxmax(), 'th1']:.2f};"
      f"  largest shape residual: {th_scan.resid.max():.1%} of peak")
print(f"the th1 = 0.15 slice reported above tops out at"
      f" {scan.lag.abs().max():.2%} and {scan.resid.max():.1%}, so a maximum"
      " over one slice")
print("is not a maximum over the model class.")'''))

cells.append(code(r'''slow = channelling_case(0.03, n=800)
fast = channelling_case(0.40, n=800)

fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.4), sharey=True)
for ax, case, ttl in ((axes[0], slow, "minority zone SLOW (a stagnant pocket)"),
                      (axes[1], fast, "minority zone FAST (a channel)")):
    t_, sim, comp = case["t"], case["sim"], case["comp"]
    pk = sim.max()
    ax.plot(t_, sim / pk, color="k", lw=2.0,
            label="two-zone bed (the 'measurement')")
    ax.plot(t_, comp / comp.max(), "o", ms=3.2, mfc="none", color="tab:red",
            markevery=max(1, len(t_) // 60),
            label=r"computed from the first curve with the fitted $D_L$")
    ax.axvline(peak_time(t_, sim), color="k", lw=0.8, ls=":")
    ax.axvline(peak_time(t_, comp), color="tab:red", lw=0.8, ls=":")
    lo = peak_time(t_, sim)
    ax2 = ax.twinx()
    ax2.plot(t_, 100 * (comp / comp.max() - sim / pk), color="tab:blue", lw=1.3)
    ax2.axhline(0.0, color="tab:blue", lw=0.6, ls=":")
    ax2.set_ylim(-4.0, 4.0)
    ax2.set_ylabel("computed - simulated, % of peak", color="tab:blue",
                   fontsize=8)
    ax2.tick_params(axis="y", labelcolor="tab:blue", labelsize=8)
    ax2.grid(False)
    ax.set(xlabel="time, s", xlim=(lo - 8, lo + 16),
           title=f"{ttl}\ncomputed peak {case['lag']:+.2%} of transit time")
    ax.legend(fontsize=8, loc="upper left")
axes[0].set_ylabel("concentration (scaled)")
fig.tight_layout()
plt.show()

# is the sign an artefact of the injected pulse shape? repeat with a symmetric one
sym_slow = channelling_case(0.03, tail=0.4, n=800)
sym_fast = channelling_case(0.40, tail=0.4, n=800)


def sigma_over_tbar(dx):
    """sqrt(2 D dx/u^3) / (dx/u) at E&R's apparent D_L: falls as dx^-1/2."""
    return np.sqrt(2 * DL15 * dx / U15 ** 3) / (dx / U15)


sig_t = sigma_over_tbar(DX_ASSUMED)
DX_ALT = P["bed_length_min"]                  # the paper's other test section
sig_t_alt = sigma_over_tbar(DX_ALT)
alt = channelling_case(0.03, n=800, dx_probe=DX_ALT)
res_slow, res_fast = slow["resid"], fast["resid"]
display(Markdown(rf"""
**The sign of E&R's shape discrepancy picks out the slow zone, not the fast
one.** Constrained to reproduce their measured $D_L$ = {DL15} cm²/s, the
two-zone bed makes the computed curve peak *later* than the simulated one — the
sense they report — only when the zone holding the smaller share of the void is
the **slower** one. With a fast channel the sign is the wrong way round in every
case tried, between {scan.query("flux_frac >= 0.25").lag.min():+.2%} and
{scan.query("flux_frac >= 0.25").lag.max():+.2%}. Nor is it an artefact of the
injected pulse: repeating the two runs with a nearly symmetric injection instead
of the long-tailed one gives {sym_slow['lag']:+.2%} and {sym_fast['lag']:+.2%},
the same signs.

**How far that generalises, tested rather than assumed.** Scanning the volume
split as well — $\theta_1$ from {th_scan.th1.min():.2f} to
{th_scan.th1.max():.2f}, each direction, every row still tuned to the same
measured $D_L$ — the sign follows minority-by-volume-slow in
**{sign_ok} of {len(th_scan)}** rows, including the ones where the *majority*
zone is the slow one. The {len(bad)} exceptions are not a mixed verdict: they
sit at $\theta_1$ = {", ".join(f"{v:.2f}" for v in bad.th1)}, near an even
split, and none of them displaces the peak by as much as
{bad.lag.abs().max():.3%} of the transit time. Restricted to rows where the
displacement exceeds {LAG_FLOOR:.2%} — the ones where there is a sign to get
right — the rule holds in {int(big.sign_ok.sum())} of {len(big)}. So the claim
the page makes is the conditional one: *when* a two-zone bed at this $D_L$
displaces the peak measurably, the direction E&R report is the one a slow
minority zone gives.

**What is small is the peak displacement, not the asymmetry — and how small
depends on where you look.** At the published slice, $\theta_1$ = 0.15, the
largest lag is {scan.lag.abs().max():.2%} of the transit time, far below what a
chart recorder resolves. That is a maximum over one slice, not over the model
class: letting the volume fraction move as well, still tuned to the same
measured $D_L$, reaches **{th_scan.lag.abs().max():.2%}** at $\theta_1$ =
{th_scan.loc[th_scan.lag.abs().idxmax(), 'th1']:.2f}, and the shape residual
reaches {th_scan.resid.max():.1%} of peak. So the peak displacement bounds
nothing about the model class; what it says is that in the regime E&R's own
$\theta_1 \approx 0.15$ picture describes, the discrepancy is mostly *shape*.
The blue traces above are the computed curve minus the simulated one, spanning
{res_slow:.1%} of the peak in the slow case against {res_fast:.1%} in the fast
one.

**Is that visible to them?** Not answerable from the paper, and the page does
not claim it is. Appendix A gives percentage standard deviations of a five-pulse
average — {P['percent_sd_pulse_area']} % in area and
{P['percent_sd_pulse_variance']} % in variance — but those are the scatter of
*integrated moments*, and a residual of a few per cent of the peak is a
*pointwise* quantity. The two are not commensurable: a 1.4 % standard deviation
on an integral bounds pointwise noise neither above nor below. The only
pointwise statement in the paper is that below about $c_\mathrm{{max}}$/30 "the
noise level becomes too great for accurate measurements" (p. 114), which is a
statement about the tail, not the peak. So the residual is **comparable to, but
not demonstrably above, their reported reproducibility**, and what a two-zone
bed at their measured $D_L$ produces is a discrepancy of the size they would
have had to argue about — not one this page can show they must have seen.

**Δx is assumed, and it matters.** Figure 15's legend prints
{DP15} cm beads, $\mathrm{{Re}}$ = {RE15} and $D_L$ = {DL15} cm²/s, and **not**
the detector separation. Everything above takes $\Delta x$ =
{DX_ASSUMED:.0f} cm, which is what Figure 13 prints for the large-particle run;
the paper's test section was "either of two tubes, one {DX_ALT} cm and the other
{DX_ASSUMED:.0f} cm long". At {DX_ASSUMED:.0f} cm the response curve has
$\sigma/\bar t$ = {sig_t:.1%} — close to plug flow, and the argument made above
is that there is therefore not much shape for a velocity distribution to
distort. At {DX_ALT} cm it is {sig_t_alt:.1%}, {sig_t_alt / sig_t:.1f} times
larger, so that argument would be weaker at the shorter section — but running it
is not the same as arguing it. The same slow-zone case at {DX_ALT} cm gives a
lag of {alt['lag']:+.2%} against {slow['lag']:+.2%}, and a shape residual of
{alt['resid']:.1%} of peak against {res_slow:.1%}: the sign is unchanged and the
discrepancy gets **{"smaller" if alt['resid'] < res_slow else "larger"}**, not
larger, because the shorter section also gives the two zones less distance in
which to differentiate. The $\sigma/\bar t$ scaling alone predicts the wrong
direction here, which is the reason to run the alternative rather than reason
about it. What survives either way: the sign, and the fact that the shape
residual is much larger than the peak displacement.

A discrepancy large enough to see on a chart recorder therefore implies
**either** a slow region rather than a channel, **or** structure coarser than
any two-zone description tuned to {DL15} cm²/s can hold — for instance an
exchange so slow that the apparent $D_L$ is still growing with bed length, which
the two-point method would faithfully report as a length-dependent "dispersion
coefficient". The scan does not choose between those two branches, and the
metadata for this page must not either.

E&R checked for a length dependence and found none — but only for the *large*
particles, over {DX_ALT} to {P['bed_length_max']} cm. For the fine particles the
check is not reported, and it is exactly the measurement that would settle both
this and the $\Delta x$ assumption above.

**What this does not do.** It does not fit their response curves: those are not
tabulated and Figure 15's abscissa carries no scale. It reproduces the sign and
the order of magnitude of a diagnostic from their printed $D_L$ and their
printed conditions, and it says which of two channelling pictures is consistent
with it.
"""))'''))

cells.append(code(r'''display(Markdown(rf"""
### One number that bears on their explanation, from Table 1 alone — and the sentence that cuts against it

At $\mathrm{{Re}}$ = {RE15} the superficial velocity is {U15 * EPS15:.3f} cm/s,
and Table 1 gives the minimum fluidisation velocity of those same {DP15} cm
beads as {UMF15} cm/s. The anomalous run therefore sits at
**{U15 * EPS15 / UMF15:.0%} of incipient fluidisation** — a number the paper
never states, from two printed values on different pages, transcribed
independently.

It is offered as corroboration because of what E&R write on p. 120: "the
transition is sharper with the material with the narrowest size range, occurring
only at velocities near those giving rise to bed expansion". But **the sentence
immediately before that one, in the same paragraph, points the other way**:
"However, this phenomenon occurs at velocities considerably below those at which
either bed expansion or fluidisation takes place."

Both sentences are printed, three lines apart, and the paper is in tension with
itself here. {U15 * EPS15 / UMF15:.0%} of $u_\mathrm{{mf}}$ is what the *quoted*
sentence would lead you to expect and is hard to call "considerably below", so
the number corroborates one half of the paragraph and sits awkwardly against the
other. Quoting only the half it agrees with would be the easy version of this
paragraph and the wrong one.

Two further limits on what {U15 * EPS15 / UMF15:.0%} can carry. Table 1 prints
that $u_\mathrm{{mf}}$ with a tilde — it is the authors' own approximate value —
and the anomaly is reported for the sand and the diakon as well, whose
tabulated $u_\mathrm{{mf}}$ are
{float(tab1.loc[tab1.material == "sand", "u_mf_cm_s"].iloc[0])} and
{float(tab1.loc[tab1.material == "diakon", "u_mf_cm_s"].iloc[0])} cm/s against
Figure 15's single quoted run. So what the number establishes is narrower than
"the bed was near fluidisation": it is that this particular run was not in a
regime where incipient fluidisation is orders of magnitude away, so a mechanism
involving particle motion is not excluded on velocity grounds — while the
sentence above says the authors saw the effect where it was.
"""))'''))

cells.append(md(r"""### Blind spots — claims this page does not make

- **It does not test whether a packed bed obeys Eq. (1).** Every simulation here
  assumes it. The paper's evidence for the assumption is its Figure 13
  comparison, whose response curves are not tabulated and are not reproduced.
- **It does not test the Reynolds-number definition.** Only one gas pair was
  ever run, so the $\nu$ in Re is untested; the authors say so themselves and
  nothing here improves on it.
- **It does not test the particle-size independence.** Figure 9's five series
  were digitised without labels, on purpose, because Eq. (18) contains no
  particle-size term. The paper's claim that there is "no consistent trend with
  particle size" is therefore *not* checked here; checking it would need the
  labels and a maintainer-reviewed shape assignment.
- **The agreement against Figure 9 is a goodness of fit, not a prediction.**
  $\gamma$ and $\beta$ were fitted to these markers. The independent content of
  the comparison is the *shape* — that the maximum exists and sits where it does
  — and the recovery of $\beta$ from a fresh reading of the figure.
- **The moment-recovery check is about the solver, not the physics.** It shows
  that pymrm plus Eqs. (6) and (7) return what was put in, and it quantifies the
  numerical dispersion that would otherwise be mistaken for the real thing. It
  says nothing about whether $D_L$ = 4.07 cm²/s is what that bed actually did.
- **Neither that check nor the area ratio tests the parameter chain.** Both are
  computed inside the simulation, so both are blind to a wrong voidage, a wrong
  $\nu$ or the wrong row of Table 1 — the break table measures exactly how blind.
- **Three markers were not recovered**, one row was moved off a bare stretch of
  the dashed curve, seven merged-cluster centres were re-fitted and one missed
  marker was added; all of it is in the sidecar, with the effect on the headline
  of including the most nearly resolvable of the three. The review overlay is
  `queue_cases/A2.5/review/fig9-overlay.png`, and no conclusion turns on any
  single marker.
- **The liquid crossover Reynolds number is not reproduced**, and is left
  unexplained rather than repaired.
- **It does not test the outlet boundary condition**, and the break table proves
  it cannot: with 80 cm of bed past the second detector, a Dirichlet outlet and
  a zero-gradient outlet return the same numbers. That question belongs to
  [`A2.1`](../A2.1-danckwerts-boundary-conditions/).
- **The two-zone section is a demonstration of sufficiency, not an
  identification.** One number, an apparent $D_L$, cannot fix three parameters;
  the scan says what the *sign* of the shape discrepancy implies, and no more.
  In particular the peak displacement is **not** a bound on what a two-zone
  description can explain: it is a maximum over whichever slice of
  $(\theta_1, u_1/u_2)$ is scanned, and it more than quadruples when the volume
  fraction is allowed to move.
- **Figure 15's detector separation is assumed, not printed.** The 100 cm of
  Figure 13 is used; the paper's other test section is 21.3 cm, and the section
  above reports what that does."""))

cells.append(md(r"""## Reuse

**The measurement chain is the reusable part.** `simulate_experiment` plus
`invert_two_point` is a complete virtual tracer experiment: give it a velocity
and a dispersion coefficient and it returns what a two-detector rig would
report. Point it at any `S4` model — a reactor with reaction, an adsorbing bed,
a membrane module — and it answers "what apparent $D_L$ would a residence-time
experiment assign to this?" That is often the more useful number, because the
apparent one is what a fitted model inherits.

**`TwoZoneBed` generalises immediately.** Two mobile zones with lateral exchange
is also the structure of a mobile/stagnant adsorbent bed, of the bubble and
emulsion phases of a fluidised bed
([`E2.1`](../E2.1-kunii-levenspiel-bubbling-bed/)), and of a trickle bed's
wetted and dry regions. Only the exchange coefficient and the velocity split
change; the assembly does not.

**Take the numerical-dispersion budget with you.** $u h/2 + (\theta-\tfrac12)
u^2\Delta t$ is not specific to this problem. Any transient convective
calculation on a first-order upwind flux carries it, and in an *inverse* problem
it does not show up as a wiggle — it shows up as a plausible parameter value.

**Related pages.** [`A2.1`](../A2.1-danckwerts-boundary-conditions/) (what a bed
Péclet number does to conversion, and the steady-state form of the same
truncation error), [`A2.3`](../A2.3-taylor-aris-dispersion/) (where a dispersion
coefficient comes from when the velocity profile is known),
[`J1.5`](../J1.5-ldf-breakthrough/) (the same operator with adsorption), `A2.6`
(Gunn's stochastic model of the same quantity).

**On Gunn.** Gunn's 1993 note *On axial dispersion in fixed beds* replots these
very data — his Figure 1, "Edwards and Richardson [7]" — against his own
stochastic model. His Figure 1 carries **two** gas-phase data sets, Gunn &
Pryce's frequency-response measurements and E&R's pulse-response ones, and the
sentence usually quoted from it — "of particular interest is the maximum that is
found in the experimental results at a Reynolds number of $\approx$ 4" — is
about that combined set, not about E&R's points alone. Eq. (18)'s own maximum is
at Re = 3.20. Read as a statement about E&R's data specifically it would be an
over-reading of Gunn, so the comparison on this page is offered as "another
reader of the same figure put the maximum here", not as an independent
determination. Gunn also quotes the argon–air Schmidt group as 0.77 where
Edwards and Richardson use 0.72. That paper is
catalogued separately as `A2.6` and is deliberately *not* built into this page:
it prints its Eq. (42) but not the function $p(\mathrm{Re})$ that makes it
evaluable, which lives in Gunn's 1969 and 1987 papers.

**Cite the source, not this page:** Edwards, M. F. & Richardson, J. F., *Gas
dispersion in packed beds*, Chemical Engineering Science **23**(2) 109–123
(1968),
[doi:10.1016/0009-2509(68)87056-3](https://doi.org/10.1016/0009-2509(68)87056-3)."""))

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
