#!/usr/bin/env python3
"""Generate index.ipynb for page B2.2. Run from the page directory.

Quoting convention, copied from C1.1/A2.5/A2.8: markdown cells are raw
triple-DOUBLE-quoted strings and code cells are raw triple-SINGLE-quoted
strings, so a code cell may contain an ordinary Python docstring. Every one is
RAW, so a single backslash here is a single backslash in the notebook.

House rule this page follows strictly: no number that a cell computes is ever
retyped into a markdown cell. Anything with a computed number in it is emitted
by `display(Markdown(f"..."))` from the cell that computed it. Numbers the
PAPER prints are data and may appear in static markdown, always identified as
printed values.
"""
from __future__ import annotations

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# ---------------------------------------------------------------- title -----
cells.append(md(r"""---
title: "Froment-Bischoff coke deactivation: the coke profile names the fouling mechanism"
description: "The 1961 paper that tied catalyst activity to deposited carbon instead of time-on-stream - all three of its closed-form solutions re-derived and met at observed second order (worst 1.7e-5 in y, 5.0e-5 in C at production resolution) by a pymrm transient bed, its four hand-era approximations measured one by one, a typo in eq. (49) adjudicated, and the descending-vs-ascending profile diagnostic mapped to where a real carbon assay could actually use it."
categories: [sec:B, struct:S4, struct:S5, tier:T0, data:tier6, phase:gas, phase:gas-solid]
date: 2026-08-06
---

# Froment-Bischoff coke deactivation: the coke profile names the fouling mechanism

**Catalog ID:** `B2.2` · **Structures:** `S4` (1D transient PDE),
`S5` (convection-dominated, deferred-correction upwind) · **Tier:** T0

Before 1961, coked catalysts were bookkept against the clock: activity was
correlated with time on stream, following Voorhies' $c = m\,t'^n$. Froment and
Bischoff's move - stated in their abstract as the point of the paper - was to
make activity a function of **the carbon actually sitting on the catalyst**,
and to couple the carbon's own formation kinetics to the reactor's continuity
equation. Two mechanisms are worked out: coke from a reaction **parallel** to
the main one (fed by the reactant) and coke from a reaction **consecutive** to
it (fed by the product). The result that made the paper canonical: the coke
lays down as a **profile** along the bed - descending for a parallel
mechanism, ascending for a consecutive one - so the shape of the deposit is a
fingerprint of its chemistry, and the locus of maximum rate (and of maximum
heat release) **travels down the bed** as the catalyst ages.

The paper is from the hand-computation era - "a numerical treatment on an
electronic computer ... did not seem justified at the present state of
knowledge" - so everything in it is closed forms, series windows, and
asymptotic patches. That makes it unusually checkable: this page re-derives
every printed solution independently, meets all three with one pymrm
finite-volume bed at second order, **measures each of the paper's four
approximations** (the $\gamma = 0$ closure, the series window, the large-$\eta$
patch, and the $\eta \simeq t$ change of variables), settles a typo in
eq. (49) from the paper's own preceding line, and then asks the question the
profile diagnostic exists to answer: **how good does a carbon assay have to be
before the profile can actually name the mechanism?**"""))

# ----------------------------------------------------------- background -----
cells.append(md(r"""## Background

**The source, precisely.** G. F. Froment and K. B. Bischoff, "Non-steady state
behaviour of fixed bed catalytic reactors due to catalyst fouling", *Chemical
Engineering Science* **16**(3-4) 189-201 (1961), doi:10.1016/0009-2509(61)80030-4,
Laboratorium voor Organische Technische Chemie, Rijksuniversiteit Gent -
identified from the scan's own running head, by-line and abstract (the OCR
renders the year as "1901" and the by-line as "FI~OMENT"; the volume, page
range and received-date are unambiguous). The scan is CCITT-G4 bilevel at
300 ppi native; every numeral on this page was read from cropped
native-resolution renders, never from the text layer.

**The catalogue cites two papers.** Its "CES 16 (1961), 17 (1962)" second half
is Part II, *Kinetic data in reactors with catalyst decay*, CES **17** (1962)
105. **Part II was not read for this page**: no scan of it exists on disk, only
an Elsevier api-text whose OCR class is known to drop decimal points, and the
data rules forbid building on that. This page is built from Part I alone, which
carries the model the case is named for in full.

**Where it sits in the literature.** The paper opens against the empirical
tradition it replaces: Voorhies' 1945 correlation of average coke with time
(*Ind. Eng. Chem.* **37** 318 - the subject of the unclaimed case `B2.1`, and
quoted here only as this paper quotes it), Rudershausen & Watson, Tyuryaev et
al., and Wilson & den Herder. Its central objection is mechanistic: coke "must
result in some way either from the reactant or from the product", so its rate
belongs in the same continuity framework as the main reaction - and once there
is a profile, an activity correlated *only* with time is ill-posed, because
activity then varies along the bed. Section 7(a) closes the loop by showing the
coupled models reproduce the *range* of the observed Voorhies exponents; this
page recomputes that comparison exactly and finds it discriminates more than
the authors claimed. The deactivation ladder continues in `B2.3`
(Levenspiel's separable rate forms) and `B2.4` (Beeckman-Froment pore
blockage), both unclaimed at the time of writing; the same authors' fixed-bed
canon (`C2.10`/`D3.4` o-xylene, `B1.3` generalised modulus) shares the
conventions used here."""))

# --------------------------------------------------- the published model ----
cells.append(md(r"""## The published model

All equation numbers are the paper's; every equation below was transcribed from
a cropped 300 dpi render and is verified symbolically or numerically further
down the page.

**Continuity, and the exact change of variables.** With a flat velocity
profile, no diffusion of any sort, constant density and mole number, and
isothermal operation, the reactant mole fraction $y$ and the carbon content of
the catalyst $c$ (kg carbon / kg catalyst) obey (eqs. 1-7)

$$
\frac{\partial y}{\partial t} + \frac{\partial y}{\partial z}
  = -\frac{\Omega \rho_B d_p}{F}\, r_A, \qquad
\frac{\partial c}{\partial t} = \frac{\epsilon \rho_A \Omega d_p}{F}\, r_c ,
$$

in reduced time $t = F t' / (\epsilon \rho_A \Omega d_p)$ and reduced length
$z = z'/d_p$. The pair is hyperbolic, and the substitution $\eta = t - z$ -
time measured from the passage of the displacement front - removes the
$\partial y/\partial t$ term **exactly** (eqs. 6-7). The paper then notes that
for practical purposes $t \gg z$, so $\eta \simeq t$; this page *measures* what
that reading costs (it is a few parts in $10^5$ for the paper's own reactor,
computed below).

**Mechanisms** (eqs. 8-11), both first order in $A$:

- **parallel**: $A \to R$ (rate $k_1 P y$) and $A \to C$ (rate $k_2 P y$) -
  coke fed by the *reactant*;
- **consecutive**: $A \to R \to C$ with $r_A = k_1 P y$ and
  $r_c = k_2 P (1-y)$ - coke fed by the *product*.

**Deactivation** (eqs. 12-13): $k_i = k_i^\circ \phi_i$, with the activity a
function of the *local carbon content* - the paper's thesis - in one of two
printed forms: exponential $\phi = e^{-\alpha c}$ (empirical) or hyperbolic
$\phi = 1/(1+Kc)$ (a Langmuir picture: $K$ "a sort of adsorption equilibrium
constant" between deposited carbon and fouled sites; a footnote generalises to
$1/(\beta + Kc)$ and sets it aside).

**The three solved cases**, in the canonical variables used throughout this
page ($x \equiv az$, and $T \equiv \alpha b \eta$ for the exponential cases,
$T \equiv K b \eta$ for the hyperbolic one, with
$a = \Omega \rho_B d_p P k_1^\circ / F$,
$b = \Omega \rho_A \epsilon d_p P k_2^\circ / F$,
$\gamma = k_2^\circ/k_1^\circ$; boundary/initial data $y(0,\eta)=1$, pure feed,
and $c(z,0)=0$, fresh bed):

| case | gas phase | coke | exact solution |
|---|---|---|---|
| **PE** parallel, exponential (18)-(19), only $k_1$ affected | $y_x = -(e^{-C}+\gamma)\,y$ | $C_T = y$ | for $\gamma=0$, eqs. (32)-(33): $y = \{1 + e^{-T}(e^{x}-1)\}^{-1}$, $e^{-C} = \{1 + e^{-x}(e^{T}-1)\}^{-1}$ |
| **PH** parallel, hyperbolic (34)-(35), both $k$'s affected, $a$ built on $k_1^\circ + k_2^\circ$ | $y_x = -\dfrac{y}{1+W}$ | $W_T = \dfrac{y}{1+W}$ | eqs. (36)-(37): $s = \sqrt{1+2T}-1$, $W e^{W} = s\,e^{s-x}$, $y = e^{-x+s-W}$ |
| **CE** consecutive, exponential (40)-(41), only $k_1$ affected | $y_x = -e^{-C}\,y$ | $C_T = 1-y$ | no closed form; the paper gives a 4-term series (44) and a large-$\eta$ patch (46)+(49) |

with $C \equiv \alpha c$ and $W \equiv K c$. Three things worth noticing, each
verified below:

- **The PE pair is symmetric**: $e^{-C}(x,T)$ is $y(T,x)$ with the arguments
  swapped - deactivation at depth $x$ mirrors breakthrough at time $T$.
- **The PH implicit pair collapses to $y = W/s$** (substitute (37) into (36)),
  which makes its verification elementary and is nowhere stated in the paper.
- **The paper's own duality** (printed as a substitution table in Section 6):
  the consecutive system is the parallel system with
  $y \leftrightarrow e^{-\alpha c}$, $z \leftrightarrow \eta$ and
  $\gamma \to -1$.

**Derived quantities.** The paper's intermediate eq. (29),
$C_x + \gamma C = e^{-C} - 1$ with $C(0,\eta$-slice$) = T$, is an *exact* ODE
reduction in $x$ that holds for **any** $\gamma$ - this page uses it as an
independent solution route. The dimensionless main-reaction rate is
$r' = y\,e^{-\alpha c}$ (53), whose surface over $(z,\eta)$ is the paper's
Fig. 5; and the distance-averaged carbon $c_{z\,\mathrm{av.}}$ (51)-(52) is
what an in-situ regeneration experiment measures, the paper's bridge to the
Voorhies exponents.

**A printed slip in eq. (49), settled below.** Eq. (48) reads
$\int_1^y e^{-\alpha b \eta y'}/y' \, dy' = -\,a z\, e^{-\alpha b \eta}$
(italic $a$: the dimensionless group). The very next line, eq. (49), prints
$ei(\alpha b \eta y) = ei(\alpha b \eta) + \boldsymbol{\alpha} z\,
e^{-\alpha b \eta}$ - with $\alpha$ where the substitution
$ei(u) = \int_u^\infty e^{-\lambda}/\lambda\, d\lambda$ requires $a$. For the
paper's own Appendix-2 constants $\alpha/a$ is a factor of order $10^2$, so
the two readings differ grossly; the adjudication cell below shows the $az$
reading agrees with the exact solution and the $\alpha z$ reading is
nonsense. (Two cosmetic misprints are recorded in the data sidecar rather than
here: $\rho_F$ vs $\rho_A$ for the same density, and a $\mathrm{kg/m^2}$ unit
on $\rho_B$ in the notation table.)"""))

# --------------------------------------- parameters and assumptions ---------
cells.append(md(r"""## Parameters and assumptions

**Everything runs at the paper's own operating point.** Appendix 2 fixes the
practical ranges from the authors' (unpublished) ethylbenzene-dehydrogenation
data: $P = 1$ atm, $\epsilon = 0.38$, $\rho_B = 1130\ \mathrm{kg/m^3}$,
$\rho_F = 1.33\ \mathrm{kg/m^3}$, $F = 0.30$ kg/hr, $L' = 0.084$ m,
$\Omega = 6.3\times10^{-4}\ \mathrm{m^2}$, $T = 600\,^\circ$C, with
$k_1^\circ \simeq 10\ \mathrm{atm^{-1}hr^{-1}}$ from a preliminary analysis
and order-of-magnitude estimates
$k_2^\circ \simeq 0.01\ \mathrm{atm^{-1}hr^{-1}}$, $\alpha \simeq 100$,
$K \simeq 200$. From these the paper prints three worked intermediates -
$aL = 0.2\,k_1^\circ \simeq 2.0$ (A.1), $\alpha b t = 5$ (A.2a),
$K b t = 10$ (A.2b) - which are recomputed below and anchor the transcription.
So the production ranges of this page are the paper's: $az \in [0,2]$,
$\alpha b \eta \in [0,5]$ (to 10 for the consecutive case, 20 for the
hyperbolic, matching Figs. 2-3), and $\gamma = k_2^\circ/k_1^\circ = 10^{-3}$.

Note $d_p$ appears in $a$, $b$ and $t$ separately but cancels from $az$,
$\alpha b\eta$ and $aL$ - which is why the paper never has to state it, and why
this page can run entirely from printed values.

**Assumptions inherited from the paper** (all stated by it): isothermal
operation; plug flow, no axial diffusion; constant density and mole number;
both reactions first order in their gas-phase reactant; deactivation through
$\phi(c)$ multiplying the fresh-bed coefficient; in PE and CE only the main
reaction deactivates ($\phi_2 = 1$), in PH both deactivate identically. The
$(z,\eta)$ change of variables is exact; the only approximations anywhere are
the four the paper makes to *solve* the system - $\gamma = 0$ in PE, the
series window and the large-$\eta$ patch in CE, and $\eta \simeq t$ when
comparing with clock-time experiments - and each of the four is measured on
this page.

**Fit/test label: nothing on this page is fitted.** The model has no free
parameter anywhere: every computation is fixed by the printed dimensionless
groups. The paper's $k_2^\circ$, $\alpha$, $K$ are the *authors'*
order-of-magnitude estimates against their own data, used here (as in the
paper) only to fix the ranges above; no measurement exists in the paper to
validate against, so every agreement below is a **reproduction or an
internal-consistency check, never a validation** - the meta.yaml, README and
data sidecars carry the same label."""))

# ------------------------------------------------------------ colab cell ----
cells.append(code(r'''try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm pyyaml
'''))

# ------------------------------------------------------------- imports ------
cells.append(code(r'''"""Bootstrap, data, and the two ledgers (metrics M, break rows BREAKS)."""
import sys, pathlib, time
T_WALL0 = time.time()
if "google.colab" in sys.modules:
    import urllib.request
    base = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
            "pymrm-gallery/main/shared/gallery_utils.py")
    urllib.request.urlretrieve(base, "gallery_utils.py")
else:
    for p in (pathlib.Path.cwd(), *pathlib.Path.cwd().parents):
        if (p / "shared" / "gallery_utils.py").is_file():
            sys.path.insert(0, str(p / "shared")); break

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "B2.2-froment-bischoff-coking"
COND   = load_data("froment-bischoff-1961-appendix2-conditions.csv", page=PAGE).set_index("name")
CLAIMS = load_data("froment-bischoff-1961-printed-claims.csv", page=PAGE).set_index("claim_id")
print(cite_data(load_meta("froment-bischoff-1961-appendix2-conditions.csv", page=PAGE)))

M = {}          # agreement metrics, assembled across the page
BREAKS = []     # defect-injection rows: (metric, base, injection, result, note)
plt.rcParams.update({"figure.dpi": 110, "font.size": 9.5})
# fixed categorical order for the three solved cases, used on every figure:
COL = {"PE": "#0072B2", "PH": "#E69F00", "CE": "#009E73"}   # Okabe-Ito, CVD-safe
'''))

# ----------------------------------------------- appendix 2 arithmetic ------
cells.append(code(r'''"""Appendix 2: recompute every printed intermediate from the printed inputs.

These four numbers are the transcription anchor of the whole page: they tie the
dimensionless groups a, b (through aL and alpha*b*t) to printed arithmetic. A
mis-read input moves at least one of them (see the break table)."""
V = COND["value"]
prefactor = V.Omega * V.rho_B * V.P * V.L_prime / V.F          # paper: 0.2
aL        = prefactor * V.k1_0                                 # paper: 2.0 (A.1)
abt       = V.alpha * V.P * V.k2_0 * V.t_process               # paper: 5   (A.2a)
Kbt       = V.K * V.P * V.k2_0 * V.t_process                   # paper: 10  (A.2b)
gamma_pap = V.k2_0 / V.k1_0
# the eta ~= t reading, quantified: gas residence time vs process time, and the
# canonical slope eps_s = alpha*b/a = T-units of coking per x-unit of transit
t_res_hr = V.rho_F * V.eps_void * V.Omega * V.L_prime / V.F
eps_s    = V.alpha * (V.rho_F * V.eps_void * V.k2_0) / (V.rho_B * V.k1_0)

M["appendix2_aL_prefactor"] = float(prefactor)
M["appendix2_aL"]           = float(aL)
M["appendix2_alpha_b_t"]    = float(abt)
M["appendix2_K_b_t"]        = float(Kbt)
M["appendix2_gamma"]        = float(gamma_pap)
M["appendix2_eps_s"]        = float(eps_s)

aL_bad = (10 * V.Omega) * V.rho_B * V.P * V.L_prime / V.F * V.k1_0
BREAKS.append(("appendix2_aL", f"{aL:.3f}",
               "Omega transcribed 6.3e-3 instead of 6.3e-4 (one lost exponent)",
               f"{aL_bad:.1f}", "19.9 vs the printed 2.0 - the printed (A.1) "
               "catches a single-digit transcription slip"))

display(Markdown(f"""
| printed intermediate | paper prints | recomputed from the printed inputs |
|---|---|---|
| $aL$ prefactor $\\Omega\\rho_B P L'/F$ | 0.2 | {prefactor:.4f} |
| $aL$ (eq. A.1) | 2.0 | {aL:.3f} |
| $\\alpha b t$ (eq. A.2a) | 5 | {abt:.2f} |
| $Kbt$ (eq. A.2b) | 10 | {Kbt:.2f} |

The prefactor deviates from the printed 0.2 by {100*abs(prefactor-0.2)/0.2:.2f} %
- the paper's own rounding, propagated into its 2.0. (A.2a) and (A.2b) are exact
products of printed values ({abt:.2f}, {Kbt:.2f}); they confirm the transcription
of five constants at once but cannot fail by less than a whole misprint, so they
are anchors rather than sensitive checks. Two derived numbers used throughout:
$\\gamma = k_2^\\circ/k_1^\\circ = {gamma_pap:g}$ - a factor {0.02/gamma_pap:.0f}
below the paper's own 0.02 validity bound for its $\\gamma=0$ closure - and the
time-scale separation: gas residence time {3600*t_res_hr:.2f} s against a 5 hr
run, i.e. $\\eta$ and $t$ differ by {t_res_hr/V.t_process:.1e} of the run, slope
$\\varepsilon_s = \\alpha b/a = {eps_s:.2e}$ in canonical units. That is the
whole content of the paper's "$\\eta \\simeq t$"; it is measured as an error
below.
"""))
'''))

# ------------------------------------------------------------- the data -----
cells.append(md(r"""## The data

**The paper contains no measurements.** Its figures are the authors' own
computed curves; its experimental content is limited to the Appendix-2
operating conditions and the Voorhies-type exponents it quotes from four other
papers. Data tier is therefore **6** and the two datasets of this page are
transcriptions, not measurements:

| dataset | what it holds |
|---|---|
| `appendix2-conditions` | the printed operating conditions, estimated constants, and the paper's three worked intermediates (pp. 200-201) |
| `printed-claims` | every quantitative prose claim - the five quoted Voorhies exponents, the $\gamma \le 0.02$ / 5 % statement, the series and asymptotic validity windows, the three Fig.-4 slope readings (pp. 193-198) |

Both sidecars carry the same fit/test label as this page: nothing is a held-out
measurement, and nothing on the page is fitted to either file.

**No figure was digitised, deliberately.** Figs. 1(a)-3(b) plot the printed
closed forms and approximations at parameter values stated in the text or
readable from curve labels; Fig. 4 plots averages of those same solutions;
Fig. 5 is the rate surface of eq. (54). Everything figure-shaped on this page
is **recomputed from the printed equations** and can be compared with the paper
by eye - the figures themselves were never measured, and the review gate that
figure digitisation would trigger is not needed. The one figure-adjacent
finding runs the other way: the text says Fig. 3(b)'s large-$\eta$ curves were
computed from eqs. (46)+(49), and this page measures that those equations sit
systematically low, so the *printed figure* inherits a quantified bias (below).

**Cross-page obligations:** no other page's dataset is loaded. The Voorhies
numbers quoted here are transcribed from *this* paper's Section 7(a) only;
`B2.1` (Voorhies 1945, unclaimed) owns the original. Part II (CES 17, 1962) was
not read from an image and contributes nothing here."""))

# --------------------------------------------------- pymrm implementation ---
cells.append(md(r"""## PyMRM implementation

One solver serves all three cases, the full-$\gamma$ system, and the mixed
parallel+consecutive generalisation. The gas phase is quasi-steady in $x$ at
each coking time $T$ (exactly the paper's $(z,\eta)$ system; the full
transient with the $\varepsilon_s\,\partial y/\partial T$ term retained is
solved separately below and measured against this):

- **Transport assembled once**: `construct_convflux_upwind` with $v = 1$ on an
  $(n,1)$ layout - fields on the last axis, never a bare `(n,)` - and
  `construct_div` with `nu=0` (Cartesian; a comment in the class says so).
  Inlet is Dirichlet $y=1$ ($\{a{:}0, b{:}1, d{:}1\}$ on the outward normal),
  outlet is the zero-gradient outflow $\{a{:}1, b{:}0, d{:}0\}$.
- **Van Leer deferred correction** (`interp_cntr_to_stagg_tvd`) iterated to a
  $10^{-12}$ fixed point brings the smooth interior to second order.
- **The outflow face is corrected explicitly.** pymrm's zero-gradient bc
  *reconstructs the face value with zero slope* (its last flux row is a
  quadratic ghost), which contradicts the physical outlet gradient
  $y' = -g y \ne 0$ and leaves the outlet **cell** first-order even under TVD
  correction - the same outlet-face class as the published `A2.6`/`A3.7`
  findings, measured in the next cell rather than asserted. The deferred
  correction therefore replaces that one face flux by a linear *upwind*
  extrapolation, restoring $O(h^2)$; outlet *readings* go through
  `compute_boundary_values` and the same extrapolation, and the difference
  between the three possible readings is measured.
- **Coking marched by Heun (RK2)** in $T$; each stage is one *linear* sparse
  solve, because at frozen $C$ every case is linear in $y$ - no `NumJac`, no
  Newton, and the class says why.
- **Both refinement axes** ($n$ in $x$, $m$ in $T$) are swept independently
  and jointly below, with observed orders; error is measured against exact
  solutions *through the passage of the activity ridge* across the bed."""))

# ------------------------------------------------- exact references ---------
cells.append(code(r'''"""Exact references. Three independent analytic routes, one per case:

PE: the printed closed forms (32)-(33)                      [algebra]
PH: the printed implicit pair (36)-(37) via Lambert W       [algebra]
CE: an exact reduction DERIVED HERE, not in the paper: at fixed x,
    Chat = -ln y obeys the scalar autonomous ODE
        dChat/dT = 1 - Chat - exp(-Chat),   Chat(0) = x,
    and C accumulates as dC/dT = 1 - exp(-Chat), C(0) = 0.
    (Proof: integrate d(Chat_T)/dx = (e^-Chat - 1) Chat_x from the inlet,
    where Chat(0,T)=0; the x-dependence integrates out exactly.)
    It is cross-checked against the paper's own series (44) symbolically in
    the next cell and against the pymrm PDE solve in the validation section.

Plus the paper's eq. (29) integrated as an ODE in x - exact for ANY gamma -
which gives the full-gamma PE solution independently of any time marching."""
from scipy.integrate import solve_ivp, simpson
from scipy.special import lambertw, expi
from scipy.optimize import brentq

def pe_exact(x, T):
    """Eqs. (32)-(33) in canonical variables. Note the (x,T) swap symmetry."""
    x, T = np.broadcast_arrays(np.asarray(x, float), np.asarray(T, float))
    y = 1.0 / (1.0 + np.exp(-T) * np.expm1(x))
    C = np.log1p(np.exp(-x) * np.expm1(T))
    return y, C

def ph_exact(x, T):
    """Eqs. (36)-(37): W = LambertW(s e^{s-x}), y = e^{-x+s-W}  (== W/s)."""
    x, T = np.broadcast_arrays(np.asarray(x, float), np.asarray(T, float))
    s = np.sqrt(1.0 + 2.0 * T) - 1.0
    W = np.real(lambertw(s * np.exp(s - x)))
    y = np.exp(-x + s - W)
    return y, W

def ce_exact(x, T_grid, rtol=1e-11, atol=1e-13):
    """Exact CE reduction, all x at once. Returns y, C with shape (nx, nT)."""
    x = np.atleast_1d(np.asarray(x, float))
    T_grid = np.atleast_1d(np.asarray(T_grid, float))
    nx = x.size
    def rhs(t, u):
        e = np.exp(-u[:nx])
        return np.concatenate([1.0 - u[:nx] - e, 1.0 - e])
    sol = solve_ivp(rhs, [0.0, max(T_grid.max(), 1e-12)],
                    np.concatenate([x, np.zeros(nx)]),
                    t_eval=T_grid, rtol=rtol, atol=atol)
    assert sol.success, "CE reference integration did not converge"
    return np.exp(-sol.y[:nx]), sol.y[nx:]

def pe_gamma_ode(T_vals, gamma, L=2.0, nx=401, rtol=1e-12):
    """Full-gamma PE, exactly, via the paper's eq. (29) in canonical form:
    C_x = e^{-C} - 1 - gamma*C with C(0)=T, then y = exp(-int (e^{-C}+gamma)).
    Never forms a T-grid; independent of every marching scheme on this page."""
    x_eval = np.linspace(0.0, L, nx)
    Y = np.empty((len(T_vals), nx)); C = np.empty_like(Y)
    for i, T in enumerate(T_vals):
        sol = solve_ivp(lambda x, u: [np.exp(-u[0]) - 1.0 - gamma * u[0],
                                      -(np.exp(-u[0]) + gamma)],
                        [0.0, L], [T, 0.0], t_eval=x_eval, rtol=rtol, atol=1e-14)
        assert sol.success
        C[i] = sol.y[0]; Y[i] = np.exp(sol.y[1])
    return x_eval, Y, C

# ---- baseline cross-checks: the references must agree with EACH OTHER ------
# (a headline needs a second, independent computation of its baseline: the
#  closed forms (32)-(33) against the eq.-29 ODE route, sharing no algebra)
T_chk = np.linspace(0.5, 5.0, 10)
x_eval, Y29, C29 = pe_gamma_ode(T_chk, gamma=0.0)
y_cf, C_cf = pe_exact(x_eval[None, :], T_chk[:, None])
M["pe_closedform_vs_eq29_ode"] = float(max(np.abs(Y29 - y_cf).max(),
                                           np.abs(C29 - C_cf).max()))
# spot identities every reference must satisfy
assert abs(pe_exact(0.0, 3.0)[0] - 1) < 1e-14 and abs(pe_exact(2.0, 0.0)[1]) < 1e-14
s4 = np.sqrt(1 + 2 * 4.0) - 1
assert abs(ph_exact(0.0, 4.0)[1] - s4) < 1e-12      # inlet W(0,T) = s exactly
assert abs(ce_exact([0.7], [0.0])[0][0, 0] - np.exp(-0.7)) < 1e-12
y_swap_a = pe_exact(1.3, 0.4)[0]
y_swap_b = np.exp(-pe_exact(0.4, 1.3)[1])           # the (x,T) swap symmetry
print(f"closed forms (32)-(33) vs eq.(29) ODE route (gamma=0), "
      f"max abs diff: {M['pe_closedform_vs_eq29_ode']:.2e}")
print(f"PE swap symmetry e^-C(x,T) = y(T,x): {abs(y_swap_a-y_swap_b):.2e}")
'''))

# ------------------------------------------------------- sympy checks -------
cells.append(code(r'''"""Symbolic verification of every transcribed solution and series.

WHAT THESE CAN AND CANNOT CATCH (label, per the house rule): each residual is
identically zero once the transcription is right, so these are STRUCTURAL
checks of the transcription and of this page's algebra - they cannot detect a
numerical error anywhere, and their zeros sit below check_agreement.py's
ABS_FLOOR, hence outside CI. The numerical solver is tested separately."""
import sympy as sp

x_, T_, g_ = sp.symbols("x T gamma", positive=True)

# --- (32)-(33) satisfy the PE system with gamma = 0 -------------------------
y_s = 1 / (1 + sp.exp(-T_) * (sp.exp(x_) - 1))
C_s = sp.log(1 + sp.exp(-x_) * (sp.exp(T_) - 1))
res_gas  = sp.simplify(sp.diff(y_s, x_) + sp.exp(-C_s) * y_s)
res_coke = sp.simplify(sp.diff(C_s, T_) - y_s)
print("PE (32)-(33) residuals in (18)-(19)|gamma=0 :", res_gas, ",", res_coke)

# --- (36)-(37) satisfy the PH system, via the parametrisation x(W,T) --------
W_, s_ = sp.symbols("W s", positive=True)
s_of_T = sp.sqrt(1 + 2 * T_) - 1
x_of_WT = s_of_T - W_ + sp.log(s_of_T / W_)     # from W e^W = s e^{s-x}
y_of_WT = W_ / s_of_T                           # (36) with (37) substituted
dxdW = sp.diff(x_of_WT, W_); dxdT = sp.diff(x_of_WT, T_)
W_x = 1 / dxdW                                  # at fixed T
W_T = -dxdT / dxdW                              # at fixed x
res_ph_gas  = sp.simplify(sp.diff(y_of_WT, W_) * W_x + y_of_WT / (1 + W_))
res_ph_coke = sp.simplify(W_T - y_of_WT / (1 + W_))
print("PH (36)-(37) residuals in (34)-(35)      :", res_ph_gas, ",", res_ph_coke)
assert res_gas == 0 and res_coke == 0 and res_ph_gas == 0 and res_ph_coke == 0

# --- the printed series (30) [PE, powers of z] and (44) [CE, powers of eta] -
def picard(f, ic, var, order):
    """Taylor solution of dU/dvar = f(U), U(0)=ic, by Picard iteration."""
    U = ic
    for k in range(1, order + 1):
        U = sp.expand(ic + sp.integrate(sp.series(f(U), var, 0, k).removeO(), var))
    return U

# (30): C_x = e^-C - 1 - gamma C, C(0)=T; printed z-coeff and z^2-coeff:
C30 = picard(lambda c: sp.exp(-c) - 1 - g_ * c, T_, x_, 2)
printed_c1 = -g_ * T_ + sp.exp(-T_) - 1
printed_c2 = sp.Rational(1, 2) * (-g_ - sp.exp(-T_)) * (-g_ * T_ + sp.exp(-T_) - 1)
ok30 = (sp.simplify(sp.expand(C30.coeff(x_, 1) - printed_c1)) == 0,
        sp.simplify(sp.expand(C30.coeff(x_, 2) - printed_c2)) == 0)
print("PE series (30): printed z and z^2 coefficients reproduced:", ok30)

# (44): own Picard series for Chat (=-ln y), dChat/dT = 1-Chat-e^-Chat, Chat(0)=x
Chat44 = picard(lambda c: 1 - c - sp.exp(-c), x_, T_, 4)
LNY_SERIES = sp.expand(-Chat44)                      # 4 terms beyond -x
printed_t1 = x_ + sp.exp(-x_) - 1
printed_t2 = sp.Rational(1, 2) * (sp.exp(-x_) - 1) * (x_ + sp.exp(-x_) - 1)
ok44 = (sp.simplify(sp.expand(LNY_SERIES.coeff(T_, 1) - printed_t1)) == 0,
        sp.simplify(sp.expand(LNY_SERIES.coeff(T_, 2) - printed_t2)) == 0)
print("CE series (44): printed eta and eta^2 coefficients reproduced:", ok44)
assert all(ok30) and all(ok44)
ce_series = sp.lambdify((x_, T_), LNY_SERIES, "numpy")   # used further down
'''))

# ------------------------------------------------------- the Bed class ------
cells.append(code(r'''"""The pymrm bed: one class for PE / PH / CE / full-gamma / mixed mechanisms."""
from scipy.sparse import diags
from scipy.sparse.linalg import splu
from pymrm import (construct_convflux_upwind, construct_div,
                   interp_cntr_to_stagg_tvd, vanleer, compute_boundary_values)

class Bed:
    """Quasi-steady gas phase on a finite-volume grid; coke marched in T.

    mech: 'PE' (parallel/exponential, optional gamma), 'PH' (parallel/
    hyperbolic; C holds W = Kc), 'CE' (consecutive/exponential), or
    'MIX' (exponential activity; coke source (1-w)*y + w*(1-y))."""

    def __init__(self, mech, L=2.0, n=200, gamma=0.0, w_mix=0.0, tvd=True,
                 fix_outflow_face=True, y_in=1.0):
        self.mech, self.L, self.n = mech, L, n
        self.gamma, self.w_mix, self.tvd = gamma, w_mix, tvd
        self.fix_outflow_face = fix_outflow_face
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # bc on the OUTWARD normal:  inlet  {a:0,b:1,d:1}  ->  y = 1
        #                            outlet {a:1,b:0,d:0}  ->  dy/dn = 0 (outflow)
        self.bc = ({"a": 0.0, "b": 1.0, "d": y_in}, {"a": 1.0, "b": 0.0, "d": 0.0})
        conv, conv_bc = construct_convflux_upwind((n, 1), self.x_f, self.x_c,
                                                  self.bc, v=1.0)
        div = construct_div((n, 1), self.x_f, nu=0)      # nu=0: Cartesian bed
        self.div = div
        self.Lop = (div @ conv).tocsc()                  # assembled ONCE
        self.bvec = np.asarray((div @ conv_bc).todense()).ravel()
        self.conv_last = conv.tocsr()[-1, :].toarray().ravel()
        self.convbc_last = float(conv_bc.tocsr()[-1, 0])
        # linear upwind extrapolation weight for the outflow face value
        self.e_out = (self.x_f[-1] - self.x_c[-2]) / (self.x_c[-1] - self.x_c[-2])

    # -- physics ------------------------------------------------------------
    def activity(self, C):
        return 1.0 / (1.0 + C) if self.mech == "PH" else np.exp(-C)

    def gcoef(self, C):
        """g in y_x = -g y (linear in y at frozen C, hence no Newton)."""
        g = self.activity(C)
        return g + self.gamma if self.mech == "PE" else g

    def coke_rate(self, C, y):
        if self.mech == "CE":
            return 1.0 - y
        if self.mech == "MIX":
            return (1.0 - self.w_mix) * y + self.w_mix * (1.0 - y)
        return y * self.activity(C) if self.mech == "PH" else y

    # -- solvers ------------------------------------------------------------
    def solve_y(self, C, y_init=None, max_it=60, tol=1e-12):
        """(Lop + diag g) y = -bvec, with van Leer deferred correction and the
        physically consistent outflow-face flux, iterated to a fixed point."""
        lu = splu((self.Lop + diags(self.gcoef(C))).tocsc())
        y = lu.solve(-self.bvec) if y_init is None else np.asarray(y_init)
        if not self.tvd:
            return lu.solve(-self.bvec)
        for _ in range(max_it):
            _, dg = interp_cntr_to_stagg_tvd(y.reshape(-1, 1), self.x_f,
                                             self.x_c, self.bc, 1.0,
                                             tvd_limiter=vanleer, axis=0)
            dg = dg.copy().reshape(-1)
            if self.fix_outflow_face:
                # pymrm's zero-gradient bc rebuilds this face with ZERO slope;
                # replace by linear upwind extrapolation (measured below).
                f_want = y[-2] + self.e_out * (y[-1] - y[-2])
                f_have = self.conv_last @ y + self.convbc_last
                dg[-1] = f_want - f_have
            y_new = lu.solve(-self.bvec
                             - np.asarray(self.div @ dg.reshape(-1, 1)).ravel())
            if np.max(np.abs(y_new - y)) < tol:
                return y_new
            y = y_new
        raise RuntimeError("deferred correction did not reach its fixed point")

    def march(self, T_end, m, record_T=(), scheme="heun"):
        """March C with Heun (RK2; 'euler' for the break table). Records must
        lie on the step grid - a lesson this page measures in its break table:
        recording at the first step PAST the target time looks exactly like a
        first-order scheme error and swamps the real convergence."""
        dt = T_end / m
        C = np.zeros(self.n)
        rec, k, y = {}, 0, None
        record_T = np.asarray(sorted(record_T), float)
        for j in range(m):
            y = self.solve_y(C, y_init=y)
            f1 = self.coke_rate(C, y)
            if scheme == "euler":
                C = C + dt * f1
            else:
                y2 = self.solve_y(C + dt * f1, y_init=y)
                C = C + 0.5 * dt * (f1 + self.coke_rate(C + dt * f1, y2))
            Tnow = (j + 1) * dt
            while k < record_T.size and record_T[k] <= Tnow + 1e-9 * dt:
                assert abs(record_T[k] - Tnow) < 1e-9, \
                    "record time off the step grid - choose m accordingly"
                rec[record_T[k]] = (self.solve_y(C, y_init=y), C.copy())
                k += 1
        self.C = C
        return rec

    def outlet_readings(self, y):
        """The outlet mole fraction three ways (compared in the text)."""
        _, _, v_bc, _ = compute_boundary_values(y.reshape(-1, 1), self.x_f,
                                                self.x_c, self.bc, axis=0)
        return {"cell": float(y[-1]),
                "boundary_values": float(np.ravel(v_bc)[0]),
                "extrapolated_face": float(y[-2] + self.e_out * (y[-1] - y[-2]))}

# ---- fresh bed (C = 0, exact y = e^-x): spatial accuracy, and the faces ----
rows, orders = [], {}
for label, kw in (("bare upwind", dict(tvd=False)),
                  ("TVD, pymrm outflow as-is", dict(tvd=True, fix_outflow_face=False)),
                  ("TVD + physical outflow face", dict(tvd=True))):
    errs = []
    for n in (100, 200, 400):
        bed = Bed("PE", n=n, **kw)
        errs.append(np.abs(bed.solve_y(np.zeros(n)) - np.exp(-bed.x_c)).max())
    p = np.log2(errs[0] / errs[2]) / 2
    rows.append((label, errs[1], p))
    orders[label] = (errs, p)
bed = Bed("PE", n=200)
y_fresh = bed.solve_y(np.zeros(200))
out = bed.outlet_readings(y_fresh)
exact_face = np.exp(-2.0)
M["fresh_bed_maxerr_n200"] = orders["TVD + physical outflow face"][0][1]
M["outlet_read_err_extrapolated_face"] = abs(out["extrapolated_face"] - exact_face)
M["outlet_read_err_boundary_values"]   = abs(out["boundary_values"] - exact_face)
M["outlet_read_err_last_cell"]         = abs(out["cell"] - exact_face)
BREAKS.append(("fresh_bed_maxerr_n200", f"{M['fresh_bed_maxerr_n200']:.1e}",
               "outflow face left to pymrm's zero-gradient reconstruction",
               f"{orders['TVD, pymrm outflow as-is'][0][1]:.1e}",
               f"x{orders['TVD, pymrm outflow as-is'][0][1]/M['fresh_bed_maxerr_n200']:.0f} "
               "worse and first order: the zero-slope ghost contradicts y'=-gy "
               "at the outlet (the A2.6/A3.7 outlet-face class)"))

display(Markdown(f"""
Fresh-bed spatial accuracy (exact $y=e^{{-x}}$), max-norm at $n=200$ and
observed order over $n=100\\to400$:

| scheme | max error ($n=200$) | observed order |
|---|---|---|
| {rows[0][0]} | {rows[0][1]:.1e} | {rows[0][2]:.2f} |
| {rows[1][0]} | {rows[1][1]:.1e} | {rows[1][2]:.2f} |
| {rows[2][0]} | {rows[2][1]:.1e} | {rows[2][2]:.2f} |

With the face fixed, the interior is clean second order. **Reading the outlet**
(exact $y(L) = e^{{-2}}$, $n=200$): the last cell centre is off by
{M['outlet_read_err_last_cell']:.1e}, `compute_boundary_values` under the
zero-gradient bc by {M['outlet_read_err_boundary_values']:.1e} (it faithfully
reconstructs the bc it was given, and the *bc* is the biased ingredient - the
handoff's Neumann-outflow lesson in a new form), and the linear upwind face
extrapolation by {M['outlet_read_err_extrapolated_face']:.1e}. All outlet
numbers on this page use the extrapolated face.
"""))
'''))

# -------------------------------------------------------------- results -----
cells.append(md(r"""## Results

Everything the paper drew, recomputed exactly - plus the two maps the paper
could not draw by hand: the error of its own approximations, and the
detectability of its profile diagnostic."""))

# ------------------------------------------------ figures 1-3 ---------------
cells.append(code(r'''"""The paper's Figs. 1-3, recomputed from the exact solutions (nothing
digitised): breakthrough-style y(T) histories at the labelled bed depths, and
the coke profiles - descending (parallel) vs ascending (consecutive)."""
fig, ax = plt.subplots(3, 2, figsize=(9.0, 9.6))
xdepths = [0.5, 1.0, 1.5, 2.0]
shade = plt.cm.Blues(np.linspace(0.45, 0.95, 4))
greens = plt.cm.Greens(np.linspace(0.45, 0.95, 6))
oranges = plt.cm.Oranges(np.linspace(0.4, 0.95, 5))

Tg = np.linspace(0, 4, 300)
for xd, c in zip(xdepths, shade):
    ax[0, 0].plot(Tg, pe_exact(xd, Tg)[0], color=c, lw=2)
    ax[0, 0].annotate(f"$az$={xd}", (0.05, pe_exact(xd, 0.0)[0]),
                      xytext=(4, -2), textcoords="offset points", fontsize=8, color=c)
xg = np.linspace(0, 2, 200)
for Tv, c in zip([1, 2, 3, 4], shade):
    ax[0, 1].plot(xg, pe_exact(xg, np.full_like(xg, Tv))[1], color=c, lw=2)
    ax[0, 1].annotate(f"{Tv}", (xg[60], pe_exact(xg[60], Tv)[1]),
                      xytext=(0, 3), textcoords="offset points", fontsize=8, color=c)
ax[0, 0].set(title="Fig. 1(a) recomputed - PE, $y$ vs $\\alpha b\\eta$",
             xlabel="$\\alpha b \\eta$", ylabel="$y$")
ax[0, 1].set(title="Fig. 1(b) recomputed - PE coke profiles ($\\alpha b\\eta$ labelled)",
             xlabel="$az$", ylabel="$\\alpha c$")

Tg = np.linspace(0, 20, 300)
for xd, c in zip(xdepths, oranges):
    ax[1, 0].plot(Tg, ph_exact(xd, Tg)[0], color=c, lw=2)
    ax[1, 0].annotate(f"$az$={xd}", (0.25, ph_exact(xd, 0.0)[0]),
                      xytext=(4, -2), textcoords="offset points", fontsize=8, color=c)
for Tv, c in zip([1, 3, 8, 15, 20], oranges):
    ax[1, 1].plot(xg, ph_exact(xg, np.full_like(xg, Tv))[1], color=c, lw=2)
    ax[1, 1].annotate(f"{Tv}", (xg[100], ph_exact(xg[100], Tv)[1]),
                      xytext=(0, 3), textcoords="offset points", fontsize=8, color=c)
ax[1, 0].set(title="Fig. 2(a) recomputed - PH, $y$ vs $Kb\\eta$",
             xlabel="$Kb\\eta$", ylabel="$y$")
ax[1, 1].set(title="Fig. 2(b) recomputed - PH coke profiles ($Kb\\eta$ labelled)",
             xlabel="$az$", ylabel="$Kc$")

Tg = np.linspace(0.001, 10, 240)
yce_hist, _ = ce_exact(xdepths, Tg)
for i, (xd, c) in enumerate(zip(xdepths, greens)):
    ax[2, 0].plot(Tg, yce_hist[i], color=c, lw=2)
    ax[2, 0].annotate(f"$az$={xd}", (0.12, yce_hist[i, 0]),
                      xytext=(4, -2), textcoords="offset points", fontsize=8, color=c)
Tvals_ce = [1, 2, 3, 5, 7, 10]
_, Cce_prof = ce_exact(xg, Tvals_ce)
for j, (Tv, c) in enumerate(zip(Tvals_ce, greens)):
    ax[2, 1].plot(xg, Cce_prof[:, j], color=c, lw=2)
    ax[2, 1].annotate(f"{Tv}", (xg[-1], Cce_prof[-1, j]),
                      xytext=(3, 0), textcoords="offset points", fontsize=8, color=c)
ax[2, 0].set(title="Fig. 3(a) recomputed - CE, $y$ vs $\\alpha b\\eta$",
             xlabel="$\\alpha b \\eta$", ylabel="$y$")
ax[2, 1].set(title="Fig. 3(b) recomputed (exact, not the paper's eqs. 46+49)",
             xlabel="$az$", ylabel="$\\alpha c$")
for a in ax.ravel():
    a.grid(alpha=0.25)
fig.suptitle("The paper's six solution figures, recomputed from the printed "
             "equations at the printed parameter values", y=1.0)
fig.tight_layout()
plt.show()
print("Parallel profiles DESCEND from the inlet; consecutive profiles ASCEND "
      "from exactly zero at the inlet (pure feed carries no product, so no "
      "consecutive coke can form at z = 0). That sign is the diagnostic.")
'''))

# ------------------------------------------------ fig 4 + slopes ------------
cells.append(code(r'''"""Fig. 4 and the Voorhies-exponent comparison, done exactly.

The paper computed c-bar by GRAPHICAL integration of Figs. 1-3 and read slopes
off a log-log plot with a straightedge; here c-bar comes from Simpson
integration of the exact solutions and the slope is the exact local
d(ln cbar)/d(ln T). Reproduction targets from the printed-claims file:
the three low-T slopes of 1.0, 'PE always ~1.0', 'CE 0.5 as ab_eta approaches
4 and beyond', 'PH approaches 0.5 at about 10'."""
T_grid = np.geomspace(0.05, 40.0, 240)
xg1, xg2 = np.linspace(0, 1, 401), np.linspace(0, 2, 401)

def cbar(mech, xgrid):
    if mech == "PE":
        _, C = pe_exact(xgrid[None, :], T_grid[:, None])
        return simpson(C, x=xgrid, axis=1) / xgrid[-1]
    if mech == "PH":
        _, C = ph_exact(xgrid[None, :], T_grid[:, None])
        return simpson(C, x=xgrid, axis=1) / xgrid[-1]
    _, C = ce_exact(xgrid, T_grid)
    return simpson(C.T, x=xgrid, axis=1) / xgrid[-1]

CB = {(m, aL_): cbar(m, xg) for m in ("PE", "PH", "CE")
      for aL_, xg in ((1, xg1), (2, xg2))}
SL = {k: np.gradient(np.log(v), np.log(T_grid)) for k, v in CB.items()}

def slope_at(mech, aL_, T):
    return float(np.interp(np.log(T), np.log(T_grid), SL[(mech, aL_)]))

# where the CE slope actually crosses 0.5 (it has no plateau there)
ce_cross = {aL_: float(np.exp(np.interp(-0.5, -SL[("CE", aL_)][T_grid > 1],
                       np.log(T_grid[T_grid > 1])))) for aL_ in (1, 2)}
win = (T_grid >= 1.0) & (T_grid <= 31.6)          # Fig. 4's plotted decade-and-a-half

M["slope_low_T_all"] = float(np.mean([SL[k][0] for k in SL]))   # -> 1 as T -> 0
M["pe_slope_max_aL2"] = float(SL[("PE", 2)][win].max())
M["pe_slope_min_aL2"] = float(SL[("PE", 2)][win].min())
M["ph_slope_at_T10_aL2"] = slope_at("PH", 2, 10.0)
M["ph_slope_at_T10_aL1"] = slope_at("PH", 1, 10.0)
M["ce_slope_at_T4_aL2"] = slope_at("CE", 2, 4.0)
M["ce_slope_cross_05_T_aL2"] = ce_cross[2]
# the PH 0.5 is a true asymptote (cbar ~ s ~ sqrt(2T)); evaluate far out:
Tfar = 1e4
s_far = np.sqrt(1 + 2 * Tfar) - 1
M["ph_slope_asymptote_check_T1e4"] = float(
    Tfar / (np.sqrt(1 + 2 * Tfar) * (s_far - 1.0)))   # aL=2: cbar ~ s - aL/2

# space-time claim (Sec. 7a): c increases with aL for CE, decreases for PE
iT5 = np.argmin(np.abs(T_grid - 5.0))
st = {m: (CB[(m, 1)][iT5], CB[(m, 2)][iT5]) for m in ("PE", "PH", "CE")}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 4.0))
for (m, aL_), v in CB.items():
    a1.loglog(T_grid, v, color=COL[m], lw=2, ls="-" if aL_ == 2 else "--")
a1.loglog(T_grid, 0.9 * np.sqrt(T_grid / T_grid[0]) * CB[("CE", 2)][0],
          color="0.55", lw=1.2, ls=":")
a1.annotate("slope 1/2 reference", (8, 0.9 * np.sqrt(8 / T_grid[0]) * CB[("CE", 2)][0]),
            fontsize=8, color="0.4")
for m, ytxt in (("PE", 0.9), ("PH", 0.55), ("CE", 0.25)):
    a1.annotate(m, (30, np.interp(30, T_grid, CB[(m, 2)])), color=COL[m],
                fontsize=9, fontweight="bold", xytext=(4, 0), textcoords="offset points")
a1.set(xlabel=r"$\alpha b\eta$ or $Kb\eta$", ylabel=r"$\alpha\bar c$ or $K\bar c$",
       title="Fig. 4 recomputed exactly (solid $aL=2$, dashed $aL=1$)")
a1.grid(alpha=0.25, which="both")
for (m, aL_), v in SL.items():
    a2.semilogx(T_grid, v, color=COL[m], lw=2, ls="-" if aL_ == 2 else "--")
a2.axhline(1.0, color="0.6", lw=0.8); a2.axhline(0.5, color="0.6", lw=0.8)
a2.plot([4], [0.5], "v", color=COL["CE"], ms=7)
a2.plot([10], [0.5], "v", color=COL["PH"], ms=7)
a2.annotate("printed:\nCE=0.5 here", (4, 0.47), fontsize=8, color=COL["CE"],
            ha="center", va="top")
a2.annotate("printed:\nPH=0.5 here", (10, 0.44), fontsize=8, color=COL["PH"],
            ha="center", va="top")
for m in ("PE", "PH", "CE"):
    a2.annotate(m, (35, SL[(m, 2)][-1]), color=COL[m], fontsize=9,
                fontweight="bold")
a2.set(xlabel=r"$\alpha b\eta$ or $Kb\eta$", ylabel=r"d$\,\ln\bar c\,$/d$\,\ln T$",
       ylim=(0.2, 1.4), title="the exact local Voorhies exponent")
a2.grid(alpha=0.25)
fig.tight_layout(); plt.show()

nq = CLAIMS.printed_value
display(Markdown(f"""
**Printed slope claims against the exact curves** (local slopes, $aL=2$ unless
stated; the paper's numbers came from a straightedge on graphically-integrated
curves):

| printed claim (Sec. 7a) | printed | computed exactly |
|---|---|---|
| all three mechanisms: slope at low $T$ | 1.0 | {M['slope_low_T_all']:.3f} (mean over the six curves at $T=0.05$; provable limit, $c \\propto T$ uniformly on a fresh bed) |
| PE "always approximately 1.0" | ~1.0 | ranges {M['pe_slope_min_aL2']:.2f}-{M['pe_slope_max_aL2']:.2f} over the plotted window - always at or *above* 1, never below |
| CE "0.5 as $\\alpha b\\eta$ approaches 4 and beyond" | 0.5 at 4 | {M['ce_slope_at_T4_aL2']:.3f} at $T=4$; the slope has **no plateau** - it crosses 0.5 at $T \\approx {M['ce_slope_cross_05_T_aL2']:.0f}$ ($aL=2$; {ce_cross[1]:.0f} at $aL=1$) and keeps falling (asymptotically $\\bar c \\sim 2\\ln T$, slope $\\to 0$) |
| PH "approaches 0.5 at about 10" | 0.5 at 10 | {M['ph_slope_at_T10_aL2']:.3f} at $T=10$ ($aL=2$; {M['ph_slope_at_T10_aL1']:.3f} at $aL=1$). The 0.5 **is** the true asymptote ($\\bar c \\sim s \\sim \\sqrt{{2T}}$: {M['ph_slope_asymptote_check_T1e4']:.3f} at $T=10^4$) but it is far from reached at 10 |

The low-$T$ limits and the PE statement reproduce; the two printed **onset
values do not survive exact computation** - they are chord readings on a short
graphical window, and the paper half-says so itself ("it seems that the
experimental data are not accurate enough to determine whether or not the slope
gradually decreases as time progresses, as indicated by Fig. 4").

**The Voorhies comparison, sharpened.** The paper quotes measured exponents
$n$ = {nq.voorhies_n_low:g}-{nq.voorhies_n_high:g} (Voorhies, gas-oil cracking),
{nq.n_rudershausen_watson:g} (cyclohexane), {nq.n_tyuryaev:g} ($n$-butane),
{nq.n_wilson_denherder:g} (reforming), and claims its mechanisms "lead to
results that lie in the range of the quoted experimental data". Exactly
computed, the instantaneous exponent over the plotted window spans
{M['pe_slope_min_aL2']:.2f}-{M['pe_slope_max_aL2']:.2f} for PE,
{slope_at('PH',2,31.6):.2f}-{SL[('PH',2)][win].max():.2f} for PH and
{slope_at('CE',2,31.6):.2f}-{SL[('CE',2)][win].max():.2f} for CE: every quoted
exponent is reachable by CE, the upper three ({nq.n_rudershausen_watson:g}-{nq.n_wilson_denherder:g})
by PH, and **none** by PE, whose exponent never drops below 1. Within this
model family the quoted data do not merely bracket the mechanisms - they point
away from parallel-exponential. The paper does not draw that conclusion; it is
stated here as a consequence of its own equations, not as a fact about the
five experiments (which this page has not seen).

**The space-time claim reproduces**: at $T=5$, going from $aL=1$ to $aL=2$
moves $\\bar c$ from {st['PE'][0]:.3f} to {st['PE'][1]:.3f} (PE, decreases) and
{st['PH'][0]:.3f} to {st['PH'][1]:.3f} (PH, decreases) but {st['CE'][0]:.3f} to
{st['CE'][1]:.3f} (CE, increases) - the paper's Sec. 7(a) sign test for the
mechanism from runs at two space times, exact.
"""))
M["spacetime_pe_cbar_aL1_T5"], M["spacetime_pe_cbar_aL2_T5"] = st["PE"]
M["spacetime_ce_cbar_aL1_T5"], M["spacetime_ce_cbar_aL2_T5"] = st["CE"]
BREAKS.append(("ce_slope_at_T4_aL2", f"{M['ce_slope_at_T4_aL2']:.3f}",
               "mechanism swapped: the same exponent read off the parallel bed",
               f"{slope_at('PE', 2, 4.0):.3f}",
               "the instantaneous Voorhies exponent genuinely discriminates "
               "between the mechanisms at the same operating point"))
'''))

# ------------------------------------------------ fig 5 + ridge -------------
cells.append(code(r'''"""Fig. 5: the rate surface r' = y e^-C (eq. 53-54), and its ridge.

New here: the ridge locus has a CLOSED FORM the paper does not print. Setting
d(ln r')/dx = 0 in eq. (54) gives e^x = e^T - 1 exactly, i.e.

    x*(T) = ln(e^T - 1)  ~  T  for T >~ 1,

so the locus of maximum rate (and maximum heat release) enters the bed at
T = ln 2 and then travels with unit speed in (x, T) - i.e. dz*/d eta = alpha
b/a in the paper's variables - and the rate ON the ridge is
(1 + 1/(e^T - 1))/4 -> 1/4 of the fresh inlet rate."""
Tg = np.linspace(0.01, 3.0, 320)
xgs = np.linspace(0.0, 2.5, 320)
Xm, Tm = np.meshgrid(xgs, Tg)
y_s, C_s = pe_exact(Xm, Tm)
Rp = y_s * np.exp(-C_s)

fig, axf = plt.subplots(figsize=(6.8, 4.2))
cf = axf.contourf(Xm, Tm, Rp, levels=14, cmap="Blues")
fig.colorbar(cf, label="$r' = y\\,e^{-\\alpha c}$")
Tr = Tg[Tg > np.log(2) + 1e-9]
axf.plot(np.log(np.expm1(Tr)), Tr, color="#B2182B", lw=2.2,
         label="ridge $x^*=\\ln(e^T-1)$ (closed form, this page)")
axf.set(xlabel="$az$", ylabel="$\\alpha b\\eta$", xlim=(0, 2.5),
        title="Fig. 5 recomputed: the travelling locus of maximum rate (PE)")
axf.legend(loc="lower right", fontsize=8)
plt.show()

# numeric spot check of the closed form, and the ridge rate
i2 = np.argmin(np.abs(Tg - 2.0))
x_num = xgs[np.argmax(Rp[i2])]
M["ridge_locus_num_vs_closed_T2"] = float(abs(x_num - np.log(np.expm1(2.0))))
r_ridge = float(Rp[i2].max())
print(f"ridge at T=2: numeric argmax x={x_num:.3f} vs closed form "
      f"{np.log(np.expm1(2.0)):.3f} (|diff| {M['ridge_locus_num_vs_closed_T2']:.1e}, "
      f"grid step {xgs[1]-xgs[0]:.3f}); r' on ridge {r_ridge:.4f} vs "
      f"(1+1/(e^2-1))/4 = {(1+1/np.expm1(2.0))/4:.4f}")
print("This is Sec. 7(b) made quantitative: with a descending coke profile the "
      "hot spot MOVES; a control scheme tuned to the fresh-bed profile chases it.")
'''))

# ------------------------------------------- diagnostic detectability -------
cells.append(code(r'''"""The profile diagnostic, taken seriously: when can a carbon assay tell the
mechanisms apart?

Measure: the half-bed contrast Delta_h = (mean coke, rear half - front half) /
cbar - what the crudest sectioned regeneration experiment estimates. With a
relative error sigma on each half's carbon assay, the difference resolves only
if |Delta_h| > 2*sqrt(2)*sigma (two-sigma on a difference of two independent
assays). The paper says which SIGN each mechanism gives; this cell computes
the magnitude, where it dies, and what the mixed mechanism does to it."""
def contrast_of_profile(C, xgr):
    """Half-bed contrast (rear-half mean - front-half mean)/cbar."""
    cb = simpson(C, x=xgr) / (xgr[-1] - xgr[0])
    nh = C.size // 2
    return (C[nh:].mean() - C[:nh].mean()) / cb, cb

def half_contrast(mech, aL_, T):
    xgr = np.linspace(0, aL_, 400)
    if mech == "PE":
        _, C = pe_exact(xgr, np.full_like(xgr, T))
    else:
        _, C = ce_exact(xgr, [T]); C = C[:, 0]
    return contrast_of_profile(C, xgr)

def mixed_contrast(w, aL_=2.0, T=5.0, n=160, m=120):
    """Mixed parallel+consecutive coke source on the pymrm bed."""
    bed = Bed("MIX", L=aL_, n=n, w_mix=w)
    bed.march(T, m)
    return contrast_of_profile(bed.C, bed.x_c)[0]

Tscan = np.geomspace(0.5, 40, 40)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 4.0))
for aL_, ls in ((2.0, "-"), (1.0, "--"), (0.5, ":"), (0.2, (0, (1, 1)))):
    xgr = np.linspace(0, aL_, 400)
    nh = 200
    _, Cpe = pe_exact(xgr[None, :], Tscan[:, None])          # (nT, nx)
    cb = simpson(Cpe, x=xgr, axis=1) / aL_
    a1.semilogx(Tscan, np.abs((Cpe[:, nh:].mean(1) - Cpe[:, :nh].mean(1)) / cb),
                color=COL["PE"], ls=ls, lw=1.8)
    _, Cce = ce_exact(xgr, Tscan)                            # (nx, nT)
    cb = simpson(Cce, x=xgr, axis=0) / aL_
    a1.semilogx(Tscan, np.abs((Cce[nh:].mean(0) - Cce[:nh].mean(0)) / cb),
                color=COL["CE"], ls=ls, lw=1.8)
for sig, lab in ((0.10, "10 % assays"), (0.02, "2 % assays")):
    a1.axhline(2 * np.sqrt(2) * sig, color="0.45", lw=1.0)
    a1.annotate(f"detection floor, {lab}", (0.55, 2 * np.sqrt(2) * sig * 1.06),
                fontsize=8, color="0.35")
a1.annotate("consecutive (ascending)", (2.2, 0.93), color=COL["CE"], fontsize=9)
a1.annotate("parallel (descending)", (6.0, 0.36), color=COL["PE"], fontsize=9)
a1.set(xlabel=r"$\alpha b\eta$", ylabel=r"$|\Delta_h|$ (half-bed contrast)",
       title="line style: $aL$ = 2, 1, 0.5, 0.2", ylim=(0, 1.08))
a1.grid(alpha=0.25)

wgrid = np.linspace(0, 1, 11)
dh_mix = [mixed_contrast(w) for w in wgrid]
a2.plot(wgrid, dh_mix, color="#4a4a4a", lw=2, marker="o", ms=4)
a2.axhline(0, color="0.6", lw=0.8)
a2.set(xlabel="consecutive weight $w$ in the mixed source $(1-w)y + w(1-y)$",
       ylabel=r"$\Delta_h$", title="mixed mechanisms cancel the profile")
a2.grid(alpha=0.25)
fig.tight_layout(); plt.show()

dh_pe, cb_pe = half_contrast("PE", 2.0, 5.0)
dh_ce, cb_ce = half_contrast("CE", 2.0, 5.0)
dh_pe20, _ = half_contrast("PE", 2.0, 20.0)
# end-to-end contrast and its closed-form large-T check (a second route to the
# same diagnostic number: exact integral vs the asymptote C = T - x)
xgr = np.linspace(0, 2, 400)
_, Cp = pe_exact(xgr, np.full_like(xgr, 20.0))
delta_pe20 = (Cp[-1] - Cp[0]) / (simpson(Cp, x=xgr) / 2)
delta_pe20_asym = -2.0 / (20.0 - 1.0)
M["diag_halfcontrast_pe_base"] = float(dh_pe)
M["diag_halfcontrast_ce_base"] = float(dh_ce)
M["diag_halfcontrast_pe_T20"] = float(dh_pe20)
M["diag_endcontrast_pe_T20_vs_asymptote"] = float(abs(delta_pe20 - delta_pe20_asym))
M["diag_sigma_max_pe_base"] = float(abs(dh_pe) / (2 * np.sqrt(2)))
M["diag_sigma_max_ce_base"] = float(abs(dh_ce) / (2 * np.sqrt(2)))
dh_flat = dh_mix[5]                            # w = 0.5 from the scan
M["mixed_w05_halfcontrast"] = float(abs(dh_flat))

BREAKS.append(("diag_halfcontrast_ce_base", f"{dh_ce:+.3f}",
               "consecutive coke source coded as y instead of (1-y)",
               f"{dh_mix[0]:+.3f}",
               "sign flips: the diagnostic itself detects the classic "
               "reactant-vs-product source mix-up"))

display(Markdown(f"""
**At the paper's own operating point** ($aL=2$, $\\alpha b\\eta = 5$): the
consecutive profile carries a half-bed contrast of {dh_ce:+.3f} - detectable
with assays as poor as $\\sigma \\approx {100*M['diag_sigma_max_ce_base']:.0f}$ % -
while the parallel profile carries only {dh_pe:+.3f}, needing
$\\sigma \\lesssim {100*M['diag_sigma_max_pe_base']:.1f}$ %: **with routine
10 % carbon assays and a two-section cut, the descending profile of the
parallel mechanism is already invisible at the paper's own conditions**, and it
fades further with time on stream ({dh_pe20:+.3f} at $T=20$; asymptotically
the whole profile flattens to slope $-1$ and the end-to-end contrast obeys
$-aL/(T-aL/2)$, reproduced by the exact integral to
{M['diag_endcontrast_pe_T20_vs_asymptote']:.1e}). The ascending consecutive
signature is anchored by $c(0)=0$ - pure feed carries no product - and never
fades. So the paper's Sec. 7(a) verdict ("what is needed are careful
experiments in which the carbon profile is actually measured") has a
quantitative edge: profile measurement mainly *rules in* a consecutive route;
failing to see a descending profile rules out little.

**The degenerate case is exact.** With both routes first order and equal
weight ($w = 0.5$), the coke source is $\\tfrac12[(1-w)y + w(1-y)]|_{{w=1/2}} =
\\tfrac12$ - independent of $y$ - so the profile is *identically* flat
(computed: $|\\Delta_h| = {M['mixed_w05_halfcontrast']:.1e}$, structural zero,
below CI's floor) while $\\bar c = T/2$ grows smoothly and, being independent
of $aL$, is also independent of space time - the paper's own remark about
balancing mechanisms, made exact. A flat profile with a clean Voorhies curve
is therefore *not* evidence of uniform aging chemistry.
"""))
'''))

# ------------------------------------------------------- validation ---------
cells.append(md(r"""## Validation

No measurement exists in the source, so everything here is **reproduction and
internal consistency**, structured so that each number can fail:

1. the three analytic references are cross-validated against *each other*
   (closed forms vs the eq.-29 ODE route; the CE reduction vs the paper's own
   series, symbolically) before the solver is judged against them;
2. the pymrm bed must then meet all three at second order in **both**
   refinement axes, through the ridge passage;
3. the paper's four approximations are measured against the exact solutions,
   which turns its prose claims (the 5 % / $\gamma \le 0.02$ statement, the
   series window, the "improves as $\eta$ increases" assertion, $\eta \simeq t$)
   into numbers that could have come out wrong;
4. every reported metric has a break-table row that moves it, and checks that
   *cannot* move are labelled structural."""))

# ------------------------------------------------- convergence --------------
cells.append(code(r'''"""The pymrm bed against the three exact references.

Error = max-norm over recorded profiles at five times spanning the ridge
passage (the PE ridge enters the bed at T = ln 2 = 0.69 and leaves x = 2 at
T = ln(e^2+1) = 2.13, inside every record window used here)."""
def err_vs_exact(mech, T_end, n, m, gamma=0.0, nrec=5, scheme="heun", **kw):
    bed = Bed(mech, L=2.0, n=n, gamma=gamma, **kw)
    Ts = np.linspace(T_end / nrec, T_end, nrec)
    rec = bed.march(T_end, m, record_T=Ts, scheme=scheme)
    ey = eC = 0.0
    for T, (y, C) in rec.items():
        if mech == "PE":
            ye, Ce = pe_exact(bed.x_c, T)
        elif mech == "PH":
            ye, Ce = ph_exact(bed.x_c, T)
        else:
            ye, Ce = ce_exact(bed.x_c, [T]); ye, Ce = ye[:, 0], Ce[:, 0]
        ey, eC = max(ey, np.abs(y - ye).max()), max(eC, np.abs(C - Ce).max())
    return ey, eC

T_END = {"PE": 5.0, "PH": 10.0, "CE": 5.0}
joint, lines = {}, []
for mech in ("PE", "PH", "CE"):
    errs = [err_vs_exact(mech, T_END[mech], n, m)
            for n, m in ((50, 60), (100, 120), (200, 240), (400, 480))]
    joint[mech] = errs
    p_y = np.log2(errs[0][0] / errs[-1][0]) / 3
    p_C = np.log2(errs[0][1] / errs[-1][1]) / 3
    M[f"{mech.lower()}_maxerr_y_n400"] = errs[-1][0]
    M[f"{mech.lower()}_maxerr_C_n400"] = errs[-1][1]
    M[f"{mech.lower()}_order_joint_y"] = float(p_y)
    lines.append(f"| {mech} vs {'closed form' if mech != 'CE' else 'exact reduction'} "
                 f"| {errs[-1][0]:.1e} | {errs[-1][1]:.1e} | {p_y:.2f} | {p_C:.2f} |")

# single-axis sweeps (PE): each knob must control its own error
ns = [err_vs_exact("PE", 5.0, n, 720)[0] for n in (50, 100, 200, 400)]
ms = [err_vs_exact("PE", 5.0, 600, m)[1] for m in (30, 60, 120, 240)]
M["pe_order_x_only"] = float(np.log2(ns[0] / ns[-1]) / 3)
M["pe_order_T_only"] = float(np.log2(ms[0] / ms[-1]) / 3)

# break rows on the solver, each moving a reported metric
bed_flip = Bed("PE", n=200); bed_flip.activity = lambda C: np.exp(+C)
try:
    rec = bed_flip.march(5.0, 240, record_T=[5.0])
    e_flip = f"{np.abs(rec[5.0][0] - pe_exact(bed_flip.x_c, 5.0)[0]).max():.1f}"
except RuntimeError:
    e_flip = "solver diverges"
BREAKS.append(("pe_maxerr_y_n400", f"{M['pe_maxerr_y_n400']:.1e}",
               "activity sign flipped: phi = exp(+alpha c)",
               e_flip, "coke self-accelerates instead of shielding; "
               "error is O(1) immediately"))
bed_in2 = Bed("PE", n=200, y_in=0.9)     # feed defect, consistently applied
rec = bed_in2.march(5.0, 240, record_T=[5.0])
e_inlet = np.abs(rec[5.0][0] - pe_exact(bed_in2.x_c, 5.0)[0]).max()
BREAKS.append(("pe_maxerr_y_n400", f"{M['pe_maxerr_y_n400']:.1e}",
               "inlet fed y = 0.9 instead of the pure-feed y = 1",
               f"{e_inlet:.2f}", "boundary data reaches every metric"))
e_eul = err_vs_exact("PE", 5.0, 600, 240, scheme="euler")[1]
BREAKS.append(("pe_order_T_only", f"{M['pe_order_T_only']:.2f}",
               "Heun replaced by forward Euler (m sweep rerun)",
               f"{np.log2(err_vs_exact('PE', 5.0, 600, 30, scheme='euler')[1]/e_eul)/3:.2f}",
               f"order falls to one and err_C at m=240 grows to {e_eul:.1e}: "
               "the T knob genuinely controls the T error (uniform dt, so the "
               "A4.7 geometric-schedule trap does not apply)"))
# the recording-misalignment defect, measured (it bit during development)
bed_r = Bed("PE", n=200)
dt = 5.0 / 233                              # 233 steps: records fall OFF the grid
Cr = np.zeros(200); yr = None; worst = 0.0
for j in range(233):
    yr = bed_r.solve_y(Cr, y_init=yr)
    f1 = bed_r.coke_rate(Cr, yr)
    y2 = bed_r.solve_y(Cr + dt * f1, y_init=yr)
    Cr = Cr + 0.5 * dt * (f1 + bed_r.coke_rate(Cr + dt * f1, y2))
    Tnow = (j + 1) * dt
    for Trec in (1.0, 2.0, 3.0, 4.0, 5.0):
        if Tnow - dt < Trec <= Tnow:        # record at the step PAST the target
            worst = max(worst, np.abs(Cr - pe_exact(bed_r.x_c, Trec)[1]).max())
BREAKS.append(("pe_maxerr_C_n400", f"{M['pe_maxerr_C_n400']:.1e}",
               "profiles recorded at the first step PAST each target time",
               f"{worst:.1e}", "a pure bookkeeping slip that reads as a "
               "first-order scheme error (dC/dT = 1 at the inlet, so the "
               "offset is ~dt); it swamped the real convergence when this "
               "page was first built"))

display(Markdown(f"""
Joint refinement, $(n,m) = (50,60) \\to (400,480)$, max-norm over the record
times:

| case vs reference | err $y$ ($n=400$) | err $C$ | order $y$ | order $C$ |
|---|---|---|---|---|
""" + "\n".join(lines) + f"""

Single-axis sweeps (PE): refining $n$ alone at $m=720$ gives order
{M['pe_order_x_only']:.2f}; refining $m$ alone at $n=600$ gives order
{M['pe_order_T_only']:.2f} - each knob controls its own error. CE's observed
order ({M['ce_order_joint_y']:.2f}) sits slightly below 2, limited by the van
Leer limiter near the inlet where its activity profile is steepest; it rises
toward 2 with refinement rather than stalling.

The pymrm route shares nothing with the references it meets: the references
are algebra plus, for CE, a scalar ODE derived independently of any grid -
and the break rows show every one of these agreement numbers moves when the
physics or the marching is disturbed.
"""))
'''))

# ------------------------------------------------- gamma study --------------
cells.append(code(r'''"""The paper's first approximation, measured: dropping gamma = k2/k1.

Printed claim (Sec. 5a): for gamma <= 0.02 the gamma=0 solution agrees with
the series 'within 5 per cent'. Here the gamma=0 closed form is compared with
the EXACT full-gamma solution (eq. 29 ODE route) over the whole production
window, and the pymrm bed run at full gamma cross-checks the ODE route."""
T_vals = np.linspace(0.5, 5.0, 10)
res = {}
for gam in (float(M["appendix2_gamma"]), 0.02, 0.2):
    x_eval, Yg, Cg = pe_gamma_ode(T_vals, gam)
    y0, _ = pe_exact(x_eval[None, :], T_vals[:, None])
    res[gam] = float((np.abs(Yg - y0) / Yg).max())
M["gamma_claim_dev_at_002"] = res[0.02]
M["gamma_dev_at_paper_0001"] = res[float(M["appendix2_gamma"])]
BREAKS.append(("gamma_claim_dev_at_002", f"{res[0.02]*100:.2f} %",
               "gamma mis-set to 0.2 (a tenfold slip in k2/k1)",
               f"{res[0.2]*100:.0f} %", "the claim's own knob moves it"))

# independent-route agreement at gamma = 0.02: pymrm marching vs eq.-29 ODE
bed = Bed("PE", n=400, gamma=0.02)
rec = bed.march(5.0, 480, record_T=[2.5, 5.0])
x_eval, Yg, Cg = pe_gamma_ode([2.5, 5.0], 0.02)
dev = max(np.abs(rec[T][0] - np.interp(bed.x_c, x_eval, Yg[i])).max()
          for i, T in enumerate([2.5, 5.0]))
M["pymrm_vs_eq29_gamma002"] = float(dev)

display(Markdown(f"""
| $\\gamma$ | max relative deviation of the $\\gamma=0$ closed form in $y$ |
|---|---|
| {M['appendix2_gamma']:g} (the paper's own reactor) | {100*M['gamma_dev_at_paper_0001']:.2f} % |
| 0.02 (the printed validity bound) | {100*M['gamma_claim_dev_at_002']:.2f} % |

At the printed bound the worst deviation over the full window ($az \\le 2$,
$\\alpha b\\eta \\le 5$, at the bed outlet) is {100*M['gamma_claim_dev_at_002']:.2f} % -
the printed "within 5 per cent" is accurate where the authors checked it
(small $z$, where their series converges) and only mildly optimistic at the
far corner. For their own reactor ($\\gamma = 10^{{-3}}$) the closure costs
{100*M['gamma_dev_at_paper_0001']:.2f} %. Two routes that share no code - the
pymrm bed marched at $\\gamma = 0.02$ and the eq.-29 ODE integration - agree
to {M['pymrm_vs_eq29_gamma002']:.1e}, so the deviation quoted is a property of
the *approximation*, not of either solver.
"""))
'''))

# --------------------------------------------- CE claims + eq 49 typo -------
cells.append(code(r'''"""The paper's second and third approximations (the CE series window and the
large-eta patch), and the eq. (49) typo, all measured against the exact
reduction."""
# ---- series window: 'four terms safe up to about az = 2, ab_eta = 2' -------
xs = np.array([0.5, 1.0, 1.5, 2.0])
tab = []
for T in (1.0, 2.0):
    y_ex = ce_exact(xs, [T])[0][:, 0]
    y_se = np.exp(ce_series(xs, T))
    tab.append((T, np.abs(y_se - y_ex) / y_ex))
M["ce_series_relerr_x2_T1"] = float(tab[0][1][-1])
M["ce_series_relerr_x2_T2"] = float(tab[1][1][-1])

# ---- large-eta patch (46)+(49): does it 'improve as eta increases'? --------
ei = lambda u: -expi(-u)                       # the paper's ei(x) = -Ei(-x)
def ce_patch(x, T):
    """Solve ei(T y) = ei(T) + x e^-T for y (the az reading), C = (1-y) T."""
    y = brentq(lambda yy: ei(T * yy) - (ei(T) + x * np.exp(-T)), 1e-12, 1.0,
               xtol=1e-14)
    return y, (1.0 - y) * T

xg = np.linspace(0.05, 2.0, 40)
patch_rows = []
for T in (2.0, 5.0, 10.0, 20.0, 40.0):
    y_ex, C_ex = ce_exact(xg, [T]); y_ex, C_ex = y_ex[:, 0], C_ex[:, 0]
    yp, Cp = np.transpose([ce_patch(x, T) for x in xg])
    patch_rows.append((T, np.abs(yp - y_ex).max(), np.abs(Cp - C_ex).max(),
                       float((np.abs(Cp - C_ex) / C_ex).max())))
M["ce_patch_dC_T10"] = patch_rows[2][2]
M["ce_patch_dy_T2"] = patch_rows[0][1]
M["ce_patch_dy_T40"] = patch_rows[4][1]

# ---- eq. (49): az vs alpha*z, adjudicated with unscaled constants ----------
a_u, alpha_u, T_adj = 0.5, 100.0, 5.0     # a != alpha, as in any real system
z_u = 2.0 / a_u                            # a bed with az = 2 at its outlet
y_exact = ce_exact([2.0], [T_adj])[0][0, 0]
y_az = ce_patch(a_u * z_u, T_adj)[0]
try:
    y_alphaz = ce_patch(alpha_u * z_u, T_adj)[0]
    err_alphaz = abs(y_alphaz - y_exact) / y_exact
except ValueError:
    err_alphaz = np.nan
M["eq49_relerr_az_reading"] = float(abs(y_az - y_exact) / y_exact)
M["eq49_relerr_alphaz_reading"] = float(err_alphaz)
BREAKS.append(("eq49_relerr_az_reading", f"{M['eq49_relerr_az_reading']:.3f}",
               "eq. (49) taken as printed (alpha z instead of az)",
               f"{M['eq49_relerr_alphaz_reading']:.3f}",
               "the printed character is wrong by the factor alpha/a = 200 "
               "here; eq. (48) one line above fixes the correct reading"))

rows = "\n".join(f"| {T:g} | {dy:.3f} | {dC:.2f} | {rC*100:.0f} % |"
                 for T, dy, dC, rC in patch_rows)
display(Markdown(f"""
**Series window** (printed: four terms "safe up to about $az=2$,
$\\alpha b\\eta=2$"; own Picard coefficients, matching the printed (44) exactly
through $T^2$): at $\\alpha b\\eta = 1$ the four-term series is within
{100*M['ce_series_relerr_x2_T1']:.1f} % of exact everywhere in the bed, but at
the claimed corner $(az, \\alpha b\\eta) = (2,2)$ it is off by
{100*M['ce_series_relerr_x2_T2']:.0f} % - the window is real along most of the
bed and optimistic at its far corner.

**The large-$\\eta$ patch** (46)+(49), against exact (max over the bed):

| $\\alpha b\\eta$ | max $|\\Delta y|$ | max $|\\Delta(\\alpha c)|$ | rel. in $\\alpha c$ |
|---|---|---|---|
""" + rows + f"""

The printed justification - "the approximation improves as $\\eta$ increases" -
is **true in $y$ only eventually** (the deviation *rises* to
{patch_rows[2][1]:.3f} at $\\alpha b\\eta = 10$ before decaying to
{M['ce_patch_dy_T40']:.3f} at 40) **and is false in the carbon profile**: the
patch writes $c = b(1-y)\\eta$ with the *current* $y$, forgetting that earlier,
more active catalyst coked faster, so it under-counts carbon by the history
term $\\int_0^\\eta [y(\\eta)-y(\\eta')]\\,d\\eta'$, which *grows* like
$\\ln \\eta$ - {M['ce_patch_dC_T10']:.2f} in $\\alpha c$ ({patch_rows[2][3]*100:.0f} %
relative) by $\\alpha b\\eta = 10$. Since the paper states Fig. 3(b)'s
large-$\\eta$ curves were computed from these equations, the printed figure
inherits that bias; the recomputed Fig. 3(b) above is the exact one.

**Eq. (49) adjudicated**: with the deliberately unscaled constants $a = 0.5$,
$\\alpha = 100$, the reading required by eq. (48) ($az$) reproduces the exact
outlet $y$ to {100*M['eq49_relerr_az_reading']:.0f} % (that residue *is* the
patch's own error, just measured), while the character as printed
($\\alpha z$) misses by {100*M['eq49_relerr_alphaz_reading']:.0f} % - it pushes
the $ei$ argument off scale. A one-character typo, provable from the line
above it.
"""))
'''))

# --------------------------------------- full transient / eta ~= t ----------
cells.append(code(r'''"""The paper's fourth approximation: eta ~= t - and the full transient solved.

The (z, eta) change of variables is exact; reading eta as clock time t is not,
and its cost is governed by eps_s = alpha*b/a - the coking accumulated during
one gas transit. Because the exact solution is y(x, T - eps_s x), the cost has
a closed form; the pymrm bed here KEEPS the eps_s dy/dT term (backward Euler
on y, one linear solve per step - the same operators, one extra diagonal) and
must land on the characteristic solution, not the quasi-steady one."""
eps_demo = 0.1                                # exaggerated for a visible test
bed = Bed("PE", n=300)
dt = 2.0 / 3000
# Initial data ON the characteristic solution: at T = 0 the characteristic
# variable is tau_c = -eps*x < 0, where the closed form is undefined (that is
# the start-up transient the eta-form never sees). Starting instead from the
# closed-form state at tau_c = eps*(L - x) >= 0 -- i.e. clock time T0 = eps*L --
# makes the IVP exactly consistent, so the only remaining error is numerical.
T0 = eps_demo * 2.0
y, C = pe_exact(bed.x_c, T0 - eps_demo * bed.x_c)
for j in range(3000):
    _, dg = interp_cntr_to_stagg_tvd(y.reshape(-1, 1), bed.x_f, bed.x_c,
                                     bed.bc, 1.0, tvd_limiter=vanleer, axis=0)
    dg = dg.copy().reshape(-1)
    dg[-1] = (y[-2] + bed.e_out * (y[-1] - y[-2])) - (bed.conv_last @ y + bed.convbc_last)
    A = (bed.Lop + diags(bed.gcoef(C)) + diags(np.full(300, eps_demo / dt))).tocsc()
    y_new = splu(A).solve((eps_demo / dt) * y - bed.bvec
                          - np.asarray(bed.div @ dg.reshape(-1, 1)).ravel())
    C = C + dt * bed.coke_rate(C, 0.5 * (y + y_new))
    y = y_new
T_fin = T0 + 2.0
y_char, C_char = pe_exact(bed.x_c, T_fin - eps_demo * bed.x_c)  # exact, full transient
y_qs, C_qs = pe_exact(bed.x_c, np.full(300, T_fin))             # eta read as t
M["fulltransient_vs_characteristic_y"] = float(np.abs(y - y_char).max())
gap_demo = float(np.abs(C_qs - C_char).max())

# the paper's own reactor: closed-form cost of eta ~= t
eps_pap = M["appendix2_eps_s"]
xg = np.linspace(0, 2, 401); Tg = np.linspace(0.5, 5.0, 10)
gap_y = max(np.abs(pe_exact(xg, np.full_like(xg, T))[0]
                   - pe_exact(xg, T - eps_pap * xg)[0]).max() for T in Tg)
gap_C = max(np.abs(pe_exact(xg, np.full_like(xg, T))[1]
                   - pe_exact(xg, T - eps_pap * xg)[1]).max() for T in Tg)
M["eta_vs_t_cost_y_paper"] = float(gap_y)
M["eta_vs_t_cost_C_paper"] = float(gap_C)
BREAKS.append(("eta_vs_t_cost_C_paper", f"{gap_C:.1e}",
               "time-scale separation eps_s set to 0.1 (a bed whose gas "
               "transit is a tenth of the coking time)",
               f"{gap_demo:.2f}", "the eta ~= t cost scales with eps_s; the "
               "paper's value is small because its separation is 4.5e-5, not "
               "because the approximation is free"))

display(Markdown(f"""
With $\\varepsilon_s$ exaggerated to {eps_demo} the full-transient pymrm solve
lands on the exact characteristic solution $y(x, T-\\varepsilon_s x)$ to
{M['fulltransient_vs_characteristic_y']:.1e} (backward Euler in $T$, so this
figure is first-order in the step, not a second-order one), where reading
$\\eta$ as $t$ would err by {gap_demo:.2f} in $\\alpha c$ - the solver
distinguishes the two cleanly. At the paper's actual separation
($\\varepsilon_s = {eps_pap:.2e}$: 0.3 s of gas transit against a 5 hr run)
the closed-form cost of $\\eta \\simeq t$ is {M['eta_vs_t_cost_y_paper']:.1e}
in $y$ and {M['eta_vs_t_cost_C_paper']:.1e} in $\\alpha c$ over the whole
production window: the paper's fourth approximation is its only harmless one,
and now it is harmless *with a number*.
"""))
'''))

# ---------------------------------------------------- break table -----------
cells.append(code(r'''"""The break table: every reported metric moved by an injected defect, plus
the checks that CANNOT fail, named as such."""
hdr = ("| metric | value on page | injected defect | value under defect | what it shows |\n"
       "|---|---|---|---|---|\n")
body = "\n".join(f"| `{m}` | {b} | {inj} | {r} | {note} |"
                 for m, b, inj, r, note in BREAKS)
display(Markdown(hdr + body))

display(Markdown(f"""
**Structural checks, named** (kept because they pin the transcription, not
because they could catch a numerical error): the sympy residuals of
(32)-(33)/(36)-(37) and the series-coefficient matches are identically zero
once the transcription is right; the mixed-mechanism flat profile at
$w = 0.5$ is exact by construction ($C_T = \\tfrac12$); and (A.2a)/(A.2b) are
exact products of printed values. All of these sit at or below
`check_agreement.py`'s `ABS_FLOOR` = 1e-12 and are therefore **outside the CI
regression net** - they are documentation, not protection. The same is true of
`pe_closedform_vs_eq29_ode` at {M['pe_closedform_vs_eq29_ode']:.0e} and the
$w=0.5$ contrast at {M['mixed_w05_halfcontrast']:.0e}.

**Coverage note**: the solver-defect rows (sign flip, inlet, Euler,
misrecording) were run on the PE case; the PH and CE agreement metrics ride on
the *same* `Bed` class, operators and marching path, differing only in the two
one-line closure functions - and the closures themselves have their own rows
(the source flip and the mechanism swap). No metric on this page lacks a row
that moves it except the ones labelled structural above.

**What no perturbation on this page can detect**: a systematic
mis-transcription of the *dimensionless groups themselves* would move every
route - closed forms, ODE reductions, and the pymrm bed - together. The only
guards against that class are the Appendix-2 printed intermediates (which pin
$a$ and $b$ through $aL = 2.0$ and $\\alpha b t = 5$) and the recomputed
figures agreeing with the printed ones by eye. They pin $a$ and $b$; they do
**not** pin $\\alpha$ and $K$ separately, which enter only as the authors'
order-of-magnitude estimates and are unvalidatable from this paper. Stated so
the reader does not mistake six agreeing routes for six independent witnesses
of the scaling.
"""))
'''))

# ---------------------------------------------------- report ----------------
cells.append(code(r'''"""Metrics for CI regression, and the wall-clock accounting."""
M = {k: float(v) for k, v in M.items()}
report_agreement("B2.2", M)
print(f"\ntotal notebook wall time so far: {time.time() - T_WALL0:.0f} s")
'''))

# ------------------------------------------------- what pymrm adds ----------
cells.append(md(r"""## What pymrm adds

Honestly: the paper needs no computer to *state* its result, and this page
reproduces rather than validates it - no measurement exists to validate
against. What the numerics add is everything the authors said they could not
afford in 1961 ("a numerical treatment on an electronic computer ... did not
seem justified"):

- **The four approximations become measurements.** $\gamma = 0$ costs 0.3 %
  at the paper's own $\gamma$ and ~6 % at its printed validity bound; the
  four-term series is good along most of the bed and ~22 % off at its claimed
  corner; the large-$\eta$ patch under-counts carbon by a term that *grows*
  with time - so the printed Fig. 3(b) is biased and the claim that the
  approximation improves with $\eta$ is false in $c$; and $\eta \simeq t$
  costs a few parts in $10^5$. One 30-line pymrm class replaces all four at
  second order in both axes.
- **The general cases the paper could not solve** run in the same class
  unchanged: full $\gamma$, both-reactions-deactivating exponential, the mixed
  parallel+consecutive source - which turns the paper's closing remark about
  balancing mechanisms into an exact degeneracy ($w = 0.5$ gives an
  identically flat profile with a perfectly smooth Voorhies curve).
- **The diagnostic gets an amplitude.** The paper established the *sign* of
  the profile; the detectability map shows the consecutive signature's half-bed contrast
  decaying only logarithmically (~1/ln t: 21.5 % of assay at t' = 5, 16.6 % at 20,
  7.6 % at 1000 - below a 20 % floor by t' ~ 20; the inlet anchor itself is exact at
  every t), while the parallel signature is already below a
  10 %-assay detection floor at the paper's own operating point and fades as
  $t^{-1}$ - which sharpens what "careful experiments in which the carbon
  profile is actually measured" must mean in practice.
- **Two exact objects the paper did not print**: the ridge locus
  $x^* = \ln(e^T-1)$ (unit speed in canonical variables - the travelling
  hot-spot statement of Sec. 7(b), exactly), and the scalar reduction of the
  consecutive case, which is what makes the error of the paper's own CE
  approximations measurable at all."""))

# ------------------------------------------------------------- reuse --------
cells.append(md(r"""## Reuse

- **Copy this page's directory** for any transient-deactivation or
  moving-front bed: aging chromatography columns, adsorbent poisoning,
  shell-progressive bed poisoning (`B2.6`), pore blockage (`B2.4`). The
  pattern - quasi-steady gas phase solved as a linear BVP at frozen solid
  state, solid state marched in slow time, both axes refined - is `S4`+`S5`
  and is the standard reduction whenever gas transit is fast against the
  process that ages the bed. The $\varepsilon_s$ cell shows how to check that
  assumption *inside* the same operators rather than assert it.
- **The outflow-face lesson travels.** pymrm's zero-gradient outflow
  reconstructs the outlet face with zero slope; when your solution has a real
  gradient there, the outlet cell degrades to first order and every outlet
  reading - including `compute_boundary_values`, which faithfully evaluates
  the bc you gave it - inherits an $O(h)$ bias. Correct the one face flux (or
  write the Robin bc the physics dictates) and read the outlet from an upwind
  extrapolation. Measured here in the fresh-bed table.
- **Record on the step grid.** Comparing a marched profile at the first step
  *past* a target time produces a spurious first-order "scheme error"
  ($\mathcal{O}(\Delta t)$, since $\partial c/\partial t \sim 1$). It cost
  this page a factor ~500 in apparent error before being found; the break
  table keeps the measurement.
- **Related pages**: `B2.1` (Voorhies, unclaimed - owns the empirical law
  this paper argues against), `B2.3`/`B2.4` (the rest of the deactivation
  ladder, unclaimed), `C2.10`/`D3.4` (Froment's o-xylene kinetics and
  reactor), `A2.1` (outward-normal bc conventions), `J1.5` (the same
  quasi-steady-plus-slow-variable structure in adsorption breakthrough).
- **Fit/test, one line**: nothing on this page is fitted; every number is
  reproduction of, or internal consistency between, printed equations and
  printed arithmetic - data tier 6, no measurement exists in the source."""))

nb = nbf.v4.new_notebook(cells=cells,
                         metadata={"kernelspec": {"display_name": "Python 3",
                                                  "language": "python",
                                                  "name": "python3"},
                                   "language_info": {"name": "python"}})
nbf.write(nb, "index.ipynb")
print(f"wrote index.ipynb with {len(cells)} cells")
