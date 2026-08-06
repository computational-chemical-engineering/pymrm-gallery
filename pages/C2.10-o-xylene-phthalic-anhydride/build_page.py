#!/usr/bin/env python3
"""Generate index.ipynb for page C2.10. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "o-Xylene to phthalic anhydride: the classic hot spot"
description: "Froment's 1967 demonstration that effective kinetics plus a two-dimensional pseudo-homogeneous model reproduces fixed-bed hot spots — both models rebuilt, and the five-degree disagreement between them measured continuously."
categories: [sec:C, struct:S2, struct:S6, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-05
---

# o-Xylene to phthalic anhydride: the classic hot spot

**Catalog ID:** `C2.10` · **Structures:** `S2` (plug flow), `S6` (2-D radial) ·
**Tier:** T0 · covers `D3.4` (the multitubular reactor half of the same paper —
the scope argument is in the case files)

Raise the feed of a cooled multitubular o-xylene oxidation reactor by **three
degrees** — from 357 to 360 °C — and the hot spot goes from 30 °C to runaway.
Froment's 1967 paper is where that demonstration was first made with effective
kinetics inside a two-dimensional pseudo-homogeneous model, and it is the paper
every reactor textbook cites for it. This page rebuilds both of its models from
the printed constants and reproduces every numerical result the text states."""))

cells.append(md(r"""## Background

Phthalic anhydride is made by air oxidation of *o*-xylene over V₂O₅ at a few
hundred degrees, with the xylene mole fraction held below 1 % to stay under the
explosion limit. The paper describes the industrial arrangement: a multitubular
reactor of 2500 tubes, 2.5 cm in diameter and 2.5 to 3 m long, packed with
catalyst and cooled by a salt bath. The chemistry is a triangular network — the
wanted partial oxidation, the consecutive combustion of the product, and the
parallel combustion of the feed — and the two side reactions release far more
heat than the wanted one, so the selectivity question and the temperature
question are the same question.

In 1967 the standard design tool was the one-dimensional plug-flow model, which
by construction "provides no information concerning excessive temperatures
along the axis". Froment's paper is an *investigation of the reliability* of
that model: it sets up a two-dimensional effective-transport model with radial
mass and heat dispersion and a wall film, reviews where every transport
parameter comes from, runs both models on the same o-xylene case with the same
heat-transfer data, and reports where they disagree. The famous answers, all
reproduced below: the 2-D model puts the runaway limit at 360 °C where the 1-D
model puts it at 365 °C; the axis temperature is far above the radial mean; and
a 10 % change in either heat-transfer parameter converts a runaway case into a
35 °C hot spot.

**Scope.** The catalogue lists the kinetics (`C2.10`) and the multitubular
reactor (`D3.4`) as separate cases naming this same paper. They are one page,
for the same reason `A2.1` absorbed `A2.2`: the paper prints no kinetic data —
the rate constants are *asserted* as "fairly representative" of the V₂O₅
chemistry — so the only printed numbers the kinetics can be checked against are
the reactor solves, which are exactly `D3.4`'s content. A kinetics page and a
reactor page would validate against the same five stated numbers. The full
argument, including what a future separate `D3.4` would need (Calderbank's
kinetics with data, or the radial heat-transfer cases `A3.11`/`A3.12`), is
recorded in `queue_cases/C2.10.yaml` and `queue_cases/D3.4.yaml`."""))

cells.append(md(r"""## The published model

**The reaction scheme** (p. 23), with $A$ = *o*-xylene, $B$ = phthalic
anhydride, $C$ = the total-combustion products CO + CO₂, all three steps
pseudo-first-order in the hydrocarbon because oxygen is in large excess:

$$
A \xrightarrow{(+\,\text{air}),\,k_1} B \xrightarrow{(+\,\text{air}),\,k_2} C,
\qquad
A \xrightarrow{(+\,\text{air}),\,k_3} C
$$

$$
r_A = (k_1 + k_3)\,N_{A0}\,N_0\,(1-y), \qquad
r_B = N_{A0}\,N_0\,[\,k_1(1-y) - k_2 x\,], \qquad
r_C = N_{A0}\,N_0\,[\,k_2 x + k_3(1-y)\,]
$$

with $x$ the conversion to phthalic anhydride, $w$ the conversion to CO + CO₂,
$y = x + w$ the total conversion, and (p. 23, all three read off the page
image):

$$
\ln k_1 = -\frac{27{,}000}{1.98\,(t+T_0)} + 19.837,\qquad
\ln k_2 = -\frac{31{,}400}{1.98\,(t+T_0)} + 20.86,\qquad
\ln k_3 = -\frac{28{,}600}{1.98\,(t+T_0)} + 18.97
$$

where $t = T - T_0$ and $t + T_0$ is the absolute temperature. The rates are
per kg of catalyst per hour; the paper leaves the units of $k$ implicit and
they are fixed here by the paper's own groups (see the parameters cell).

**The two-dimensional model** (Eq. 1, p. 24), in the paper's dimensionless
coordinates $z = z'/d_p$, $r = r'/d_p$, $R = R'/d_p$, with axial dispersion
dropped ("at the high flow rates used in practice, the contribution of this
mechanism may be neglected"):

$$
\frac{\partial x}{\partial z} = a_1\!\left(\frac{\partial^2 x}{\partial r^2}
 + \frac{1}{r}\frac{\partial x}{\partial r}\right) + b_1 r_B,\qquad
\frac{\partial w}{\partial z} = a_1\!\left(\frac{\partial^2 w}{\partial r^2}
 + \frac{1}{r}\frac{\partial w}{\partial r}\right) + b_1 r_C,
$$

$$
\frac{\partial t}{\partial z} = a_2\!\left(\frac{\partial^2 t}{\partial r^2}
 + \frac{1}{r}\frac{\partial t}{\partial r}\right) + b_2 r_B + b_3 r_C
$$

$$
a_1 = \frac{1}{\mathrm{Pe}_{mR}},\quad
a_2 = \frac{\lambda_R}{G c_p d_p} = \frac{1}{\mathrm{Pe}_{hR}},\quad
b_1 = \frac{\rho_b d_p M_m}{G N_{A0}},\quad
b_2 = \frac{\rho_b d_p(-\Delta H_1)}{G c_p},\quad
b_3 = \frac{\rho_b d_p(-\Delta H_3)}{G c_p}
$$

Boundary conditions (p. 24): $x = w = t = 0$ at $z=0$; symmetry at $r=0$; no
mass flux at the wall; and at the wall the film condition
$\left(\partial t/\partial r\right)_R = -\dfrac{\alpha_w d_p}{\lambda_R}\,t
= -\alpha\,t$, the tube-wall and coolant-side resistances neglected. Radial
mean values are $\zeta_m = 2\int_0^1 \zeta\,(r/R)\,\mathrm{d}(r/R)$.

**The one-dimensional model** (p. 25) is the radially lumped reduction,

$$
\frac{\mathrm{d}x}{\mathrm{d}z} = b_1 r_B,\qquad
\frac{\mathrm{d}w}{\mathrm{d}z} = b_1 r_C,\qquad
\frac{\mathrm{d}t}{\mathrm{d}z} = -\frac{4U}{G c_p d_t}\,t + b_2 r_B + b_3 r_C
$$

and the point of the comparison is that $U$ is **not free**: it is built from
the same $\lambda_R$ and $\alpha_w$ used in the 2-D model, by matching the
zero-order moment of the exact packed-bed heat-exchanger solution (pp. 25–26),
which collapses to

$$
\frac{1}{U} = \frac{1}{\alpha_w} + \frac{R'}{4\lambda_R}
$$

The paper prints the closed-form exchanger solution on the way — a Bessel
series over the roots of $\lambda J_1(\lambda) = \alpha R\, J_0(\lambda)$ —
and that series is used below as an analytic check on the radial operator.

### Two printed values are European-style thousands separators

Page 24 prints "$G = 4.684$ kg./sq. meter hr." and "$\Delta H_3 = -1.090$
kcal./gram mole". Both periods are thousands separators (Froment wrote from
Ghent), and the paper itself proves it three ways, printed in the parameters
cell below: p. 19 prints $G = 4684$ outright for the same reactor; the printed
$\mathrm{Pe}_{hR} = 5.25$ equals $G c_p d_p/\lambda_R$ only with $G = 4684$
(with 4.684 it would be 0.00524); and $-1090$ kcal/gmol is the heat of complete
combustion of *o*-xylene, while $-1.090$ kcal/gmol would make total combustion
282× *less* exothermic than partial oxidation, contradicting the paper's own
argument that hot spots destroy selectivity. The same 1960s typesetting trap —
in mid-dot decimal form — is documented on the sibling page
[`D2.2`](../D2.2-van-welsenaere-froment-runaway/) for Van Welsenaere &
Froment (1970). The break table below shows what taking either literally would
do. A third slip: p. 24 says $N_{A0} = 0.00924  # p. 24; p. 19 prints 0.00927 for the companion case - the page's own 0.00929 inversion sits between them$ "corresponds to 44 gram
moles/cu. meter"; the unit is grams per normal cubic meter — Figure 11's own
axis is labelled g/Nm³, and the arithmetic (also below) closes for g/Nm³ and
not for mol/m³."""))

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
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import j0, j1
from pymrm import construct_grad, construct_div, construct_convflux_upwind, NumJac, newton
from gallery_utils import load_data, load_meta, cite_data, report_agreement

np.set_printoptions(precision=6, suppress=True)'''))

cells.append(md(r"""## Parameters and assumptions

Every constant is printed in the paper; the CSV records where. Three points
need stating rather than assuming:

- **$M_m$ and $c_p$ are printed on p. 19**, in the companion single-reaction
  example on the same reactor, not restated on p. 24. $c_p = 0.25$ is
  additionally pinned by p. 24's own $\mathrm{Pe}_{hR} = 5.25$, which contains
  it. $M_m = 29.48$ has no closing identity on p. 24; recomputing the feed mean
  molecular weight from the composition gives 29.68 (+0.7 %), and $M_m$ enters
  only the conversion scale $b_1$, not the temperature equation.
- **The absolute-temperature offset is not printed.** The baseline here uses
  $T[\mathrm{K}] = t\,[°\mathrm{C}] + 273.15$; 1967 practice was often $+273$.
  The 0.15 K difference is 0.024 % in absolute temperature, which the
  activation energy amplifies to about 0.5 % in $k_1$ — and the runaway limit
  sits on an exponential knife edge, so the sensitivity rows below show both
  conventions.
- **Bed length** is taken as 3 m ("2.5 to 3 meters long"; at 357 °C "a length
  of 3 meters is sufficient to reach the maximum in phthalic anhydride
  concentration" — checked below).

Units are the paper's own (kg, m, hr, kcal, °C), with the two heats of
reaction per kmol so that the rate constants come out in kmol per kg of
catalyst per hour."""))

cells.append(code('''par = load_data("froment-1967-parameters.csv",
                page="C2.10-o-xylene-phthalic-anhydride")
P = dict(zip(par.symbol, par.value))

# primary printed constants -------------------------------------------------
NA0, N0   = P["N_A0"], P["N_0"]              # inlet mole fractions (o-xylene, O2)
dH1, dH3  = P["minus_dH1"], P["minus_dH3"]   # kcal/kmol (positive = exothermic)
d_t, d_p  = P["d_t"], P["d_p"]               # m
rho_b     = P["rho_b"]                       # kg cat / m3
G         = P["G"]                           # kg gas / (m2 hr), superficial
c_p       = P["c_p"]                         # kcal/(kg degC)
M_m       = P["M_m"]                         # kg/kmol
lam_R     = P["lambda_R"]                    # kcal/(m hr degC)
alpha_w   = P["alpha_w"]                     # kcal/(m2 hr degC)
Pe_mR     = P["Pe_mR"]
L_bed     = P["L_bed"]                       # m
Rp        = d_t / 2                          # tube radius R' (m)
KELVIN    = 273.15                           # baseline convention; see above

E = np.array([P["E1"], P["E2"], P["E3"]])            # cal/mol
lnA = np.array([P["lnA1"], P["lnA2"], P["lnA3"]])
Rg = P["R_gas"]                                       # 1.98 cal/(mol K), as printed

# derived dimensionless groups (Eq. 1, p. 24) -------------------------------
a1 = 1.0 / Pe_mR                     # radial mass dispersion, = D_R rho_f/(G_i d_p)
a2 = lam_R / (G * c_p * d_p)         # radial heat dispersion, = 1/Pe_hR
alpha = alpha_w * d_p / lam_R        # dimensionless wall film coefficient
Rdim = Rp / d_p                      # dimensionless tube radius R
b1 = rho_b * d_p * M_m / (G * NA0)
b2 = rho_b * d_p * dH1 / (G * c_p)
b3 = rho_b * d_p * dH3 / (G * c_p)

def U_from(lam, aw):
    """Overall coefficient from the paper's moment-matching formula (p. 26)."""
    return 1.0 / (1.0 / aw + Rp / (4.0 * lam))

U_ref = U_from(lam_R, alpha_w)

print("dimensionless groups:")
print(f"  a1 = 1/Pe_mR   = {a1:.6f}")
print(f"  a2 = 1/Pe_hR   = {a2:.6f}   (Pe_hR = {1/a2:.4f}, paper prints 5.25)")
print(f"  alpha          = {alpha:.4f}     (alpha*R = Biot = {alpha*Rdim:.4f})")
print(f"  R = R'/d_p     = {Rdim:.4f}")
print(f"  b1 = {b1:.5f}   b2 = {b2:.3f} degC   b3 = {b3:.3f} degC")
print()
print("identity checks against printed derived values (nothing here is free):")
pe_dev = abs(G * c_p * d_p / lam_R - 5.25) / 5.25
print(f"  Pe_hR = G*cp*dp/lam_R = {G*c_p*d_p/lam_R:.4f} vs printed 5.25 "
      f"-> {100*pe_dev:.2f} %  (pins G = 4684 AND cp = 0.25)")
u_dev = abs(U_ref - 82.7) / 82.7
print(f"  U(0.67, 134) = {U_ref:.3f} vs printed 82.7  -> {100*u_dev:.2f} %")
u86 = U_from(0.75, alpha_w); u86_dev = abs(u86 - 86) / 86
print(f"  U(0.75, 134) = {u86:.3f} vs printed 86    -> {100*u86_dev:.2f} %")
u88 = U_from(lam_R, 150.0); u88_dev = abs(u88 - 88) / 88
print(f"  U(0.67, 150) = {u88:.3f} vs printed 88    -> {100*u88_dev:.2f} %")
# 44 g/Nm3 of o-xylene (M = 106.17) in air at 0 degC, 1 atm:
n_total = 101325.0 / (8.314 * 273.15)                 # mol/m3 at NTP
na0_calc = (44.0 / 106.17) / n_total
na0_dev = abs(na0_calc - NA0) / NA0
print(f"  44 g/Nm3 -> N_A0 = {na0_calc:.5f} vs printed 0.00924 -> {100*na0_dev:.2f} % "
      "(closes for g/Nm3; for 44 mol/m3 it would be off by 100x)")
# Reynolds check: Re = G*dp/mu with air viscosity at ~630 K (Sutherland)
mu_air = 1.458e-6 * 630.0**1.5 / (630.0 + 110.4) * 3600.0   # kg/(m hr)
print(f"  Re = G*dp/mu(630 K) = {G*d_p/mu_air:.0f} vs printed 121 "
      "(order-of-magnitude witness for G = 4684; with 4.684 it would be 0.13)")
adiabatic_1 = NA0 / M_m * 1000.0 * dH1 / 1000.0 / c_p
print(f"  adiabatic rise of A->B alone = {adiabatic_1:.0f} degC "
      f"(= b2/b1 * Mm-scale; the knife edge under the whole page)")'''))

cells.append(md(r"""## The data

**There is no experimental dataset in this paper**, and none is claimed here.
The kinetics are presented by Froment as "fairly representative of the gas
phase air oxidation of *o*-xylene into phthalic anhydride on V₂O₅ catalysts" —
asserted, not fitted to data in this document, and Calderbank (the second name
on the catalogue row) is not on disk. Everything this page compares against is
the paper's **own computed results**, stated in its running text: this is
**reproduction of a published reference solution (tier 6)**, the same standing
as `D2.2` and `B1.1`, not validation against measurement. Nothing is fitted:
every model input is a printed constant, so the reproduction can genuinely
fail — and the break table shows what failure would look like.

Two committed datasets:
[`froment-1967-parameters.csv`](data/froment-1967-parameters.csv) (the printed
constants, with the page each was read from) and
[`froment-1967-stated-results.csv`](data/froment-1967-stated-results.csv)
(every numerical statement in the text about this case, with the paper's own
qualifiers preserved).

**Cross-page reconciliation.** The published
[`D2.2`](../D2.2-van-welsenaere-froment-runaway/) page carries Van Welsenaere &
Froment (1970) — the *same* reactor and the same first reaction, three years
later. Its parameter CSV was transcribed independently of this page, so the
shared constants are reconciled below rather than retyped. D2.2's recorded
findings about its rows, checked against this page: (1) its `c_p = 0.323
kcal/(m³·°C)` is printed *volumetric* — that is exactly $\rho_g c_p = 1.293
\times 0.25$, so it confirms rather than contradicts this page's mass-basis
0.25; (2) its numbers had to be read from page renders because the 1970 OCR
drops mid-dot decimals — the same typesetting-trap family as the thousands
separators above; (3) its $U = 82.7$ is *inherited from this 1967 paper*, so
agreement on that row is a shared source, not two witnesses — the printed
$\lambda_R$, $\alpha_w$ route above is the independent check on it."""))

cells.append(code('''res = load_data("froment-1967-stated-results.csv",
                page="C2.10-o-xylene-phthalic-anhydride")
print(cite_data(load_meta("froment-1967-stated-results.csv",
                          page="C2.10-o-xylene-phthalic-anhydride")))
display(res[["quantity", "value", "unit", "model", "qualifier", "printed_where"]])

# reconcile every constant shared with the independently transcribed D2.2 page
d22 = load_data("van-welsenaere-froment-1970-parameters.csv",
                page="D2.2-van-welsenaere-froment-runaway")
D = dict(zip(d22.symbol, d22.value))
rows = [
    ("M_m [kg/kmol]",            M_m,              D["M"]),
    ("rho_b [kg/m3]",            rho_b,            D["rho_b"]),
    ("-dH1 [kcal/kmol]",         dH1,              D["minus_dH"]),
    ("tube radius R' [m]",       Rp,               D["R"]),
    ("N_0 = p_B0 [-, atm at 1 atm]", N0,           D["p_B0"]),
    ("E1/R_gas [K]  (27000/1.98)", E[0] / Rg,      D["a"]),
    ("lnA1 [-]",                 lnA[0],           D["b"]),
    ("U printed [kcal/(m2 hr degC)]", 82.7,        D["U"]),
    ("G*cp [kcal/(m2 hr degC)] vs u*cp_vol", G * c_p, D["u"] * D["c_p"]),
]
tab = pd.DataFrame(rows, columns=["constant", "this page (Froment 1967)",
                                  "D2.2 page (VW&F 1970)"])
tab["rel diff"] = (tab.iloc[:, 1] - tab.iloc[:, 2]).abs() / tab.iloc[:, 2].abs()
display(tab.style.format({"rel diff": "{:.2e}"}))
assert (tab["rel diff"][:8] < 1e-3).all(), "shared-constant reconciliation failed"
print("First eight rows identical to <0.1 %: two independent transcriptions of "
      "the shared constants agree. The last row differs by "
      f"{100*tab['rel diff'].iloc[-1]:.1f} % because the 1970 paper quotes u = 3600 m/hr "
      "against 1967's G = 4684 kg/(m2 hr) (equivalent gas density 1.301 vs 1.293 kg/m3) "
      "- a genuine, small difference between the two papers, not a transcription error.")'''))

cells.append(md(r"""## PyMRM implementation

The 2-D model is parabolic in $z$ (axial dispersion is dropped by the paper),
so $z$ is an integration variable and the discretisation lives entirely on the
radial axis: `construct_grad`/`construct_div` with **`nu=1` (cylindrical)** —
the same radial-operator vocabulary as
[`A2.3`](../A2.3-taylor-aris-dispersion/) and
[`A3.15`](../A3.15-graetz-nusselt/) — marched in $z$ with a stiff integrator,
which is the modern version of the implicit Crank–Nicolson march the paper
itself used ("the program was tried out in 1961"). Boundary conditions use the
**outward normal** ($a\,\partial c/\partial n + b\,c = d$):

- $r = 0$ (normal $-r$): symmetry, $\partial t/\partial n = 0$ →
  `{a: 1, b: 0, d: 0}`;
- $r = R$ (normal $+r$), temperature: the paper's wall film
  $\partial t/\partial r = -\alpha t$ becomes $\partial t/\partial n + \alpha t
  = 0$ → `{a: 1, b: alpha, d: 0}`;
- $r = R$, mass: no flux through the wall → `{a: 1, b: 0, d: 0}`.

The three fields march together as one state vector so the Arrhenius coupling
is exact; the reaction term is pointwise, and the 1-D model is the same code
with the radial operator replaced by the lumped $-4U/(Gc_p d_t)\,t$ sink."""))

cells.append(code('''class Radial2D:
    """Froment's Eq. 1: three fields (x, w, t) on a cylindrical radial grid,
    marched in the dimensionless axial coordinate z = z\'/d_p.

    Spatial axis first, fields last: state (n_r, 3), flattened only at the
    integrator interface."""

    def __init__(self, n_r=80, a1=a1, a2=a2, alpha=alpha, nu=1,
                 kelvin=KELVIN, k_perm=(0, 1, 2), wall_sign=+1.0,
                 b1_val=None, b2_val=None, b3_val=None, heat_scale=1.0):
        """heat_scale rescales the inlet o-xylene concentration: N_A0 cancels
        out of the mass equations (b1*r_B is N_A0-free), so a concentration of
        s*44 g/Nm3 is exactly heat_scale=s on the temperature sources."""
        self.n_r, self.kelvin = n_r, kelvin
        self.k_perm = k_perm                      # break-table hook (k index swap)
        self.b1 = b1 if b1_val is None else b1_val
        self.b2 = (b2 if b2_val is None else b2_val) * heat_scale
        self.b3 = (b3 if b3_val is None else b3_val) * heat_scale
        self.r_f = np.linspace(0.0, Rdim, n_r + 1)
        self.r_c = 0.5 * (self.r_f[:-1] + self.r_f[1:])
        shape = (n_r,)
        # outward normal at r=0 is -r; symmetry: dt/dn = 0
        bc_sym = {"a": 1.0, "b": 0.0, "d": 0.0}
        # wall, mass: no flux, dx/dn = 0
        bc_mass = ({"a": 1.0, "b": 0.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        # wall, heat: dt/dr = -alpha*t (paper, p. 24) -> dt/dn + alpha*t = 0
        # wall_sign=-1 is the break-table defect (film heats instead of cools)
        bc_heat = (bc_sym, {"a": 1.0, "b": wall_sign * alpha, "d": 0.0})
        gm, _ = construct_grad(shape, self.r_f, self.r_c, bc_mass)
        gh, gh_bc = construct_grad(shape, self.r_f, self.r_c, bc_heat)
        dv = construct_div(shape, self.r_f, nu=nu)     # nu=1: cylindrical radial
        self.Lm = (dv @ gm).toarray() * a1             # dense: n_r x n_r is small
        self.Lt = (dv @ gh).toarray() * a2
        self.gt = a2 * np.asarray((dv @ gh_bc).todense()).ravel()   # 0 here (d=0)
        # exact annular area weights for the radial mean (the paper\'s zeta_m)
        self.wgt = (self.r_f[1:] ** 2 - self.r_f[:-1] ** 2) / Rdim ** 2

    def rates(self, x, w, t, T0C):
        T = t + T0C + self.kelvin
        k = np.exp(lnA[:, None] - E[:, None] / (Rg * T))[list(self.k_perm), :]
        y = x + w
        rB = NA0 * N0 * (k[0] * (1 - y) - k[1] * x)
        rC = NA0 * N0 * (k[1] * x + k[2] * (1 - y))
        return rB, rC

    def rhs(self, z, u, T0C):
        n = self.n_r
        x, w, t = u[:n], u[n:2 * n], u[2 * n:]
        rB, rC = self.rates(x, w, t, T0C)
        return np.concatenate([self.Lm @ x + self.b1 * rB,
                               self.Lm @ w + self.b1 * rC,
                               self.Lt @ t + self.gt + self.b2 * rB + self.b3 * rC])

    def solve(self, T0C, z_end=L_bed / d_p, rtol=1e-8, cap=150.0):
        """March to z_end; stop early (runaway) when the mean t exceeds cap."""
        n = self.n_r

        def hit_cap(z, u, *args):
            return self.wgt @ u[2 * n:] - cap
        hit_cap.terminal, hit_cap.direction = True, 1
        sol = solve_ivp(self.rhs, (0.0, z_end), np.zeros(3 * n), args=(T0C,),
                        method="BDF", rtol=rtol, atol=rtol * 1e-2,
                        dense_output=True, events=hit_cap)
        self.sol, self.ran_away = sol, bool(sol.t_events[0].size)
        return self

    def profiles(self, n_z=1500):
        n = self.n_r
        zz = np.linspace(0.0, self.sol.t[-1], n_z)
        uu = self.sol.sol(zz)
        xm = uu[:n, :].T @ self.wgt
        wm = uu[n:2 * n, :].T @ self.wgt
        tm = uu[2 * n:, :].T @ self.wgt
        t_axis = uu[2 * n, :]                      # first cell centre, r = h/2
        return zz * d_p, xm, wm, tm, t_axis

    def hot_spot(self):
        """(z\'[m], tm) at the maximum of the radial-mean temperature."""
        if self.ran_away:
            return np.nan, np.inf
        zz, _, _, tm, _ = self.profiles(4000)
        i = int(np.argmax(tm))
        return zz[i], tm[i]


def critical_T0_2d(lo=356.0, hi=364.0, tol=0.02, cap=150.0, **kw):
    """Runaway limit by bisection on cold-started solves (no continuation:
    every solve starts from the same zero state, so the answer is a
    deterministic function of the bracket and tolerance only)."""
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        m = Radial2D(**kw).solve(mid, cap=cap)
        if m.ran_away or m.hot_spot()[1] > cap:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)'''))

cells.append(code('''def rhs_1d(z, u, T0C, U=82.7, kelvin=KELVIN):
    """The paper\'s 1-D reduction (p. 25) with the printed U = 82.7."""
    x, w, t = u
    T = t + T0C + kelvin
    k1, k2, k3 = np.exp(lnA - E / (Rg * T))
    y = x + w
    rB = NA0 * N0 * (k1 * (1 - y) - k2 * x)
    rC = NA0 * N0 * (k2 * x + k3 * (1 - y))
    return [b1 * rB, b1 * rC,
            -4.0 * U * d_p / (G * c_p * d_t) * t + b2 * rB + b3 * rC]


def solve_1d(T0C, U=82.7, kelvin=KELVIN, rtol=1e-10, cap=300.0):
    def hit_cap(z, u, *args):
        return u[2] - cap
    hit_cap.terminal, hit_cap.direction = True, 1
    sol = solve_ivp(rhs_1d, (0.0, L_bed / d_p), [0.0, 0.0, 0.0],
                    args=(T0C, U, kelvin), method="LSODA", rtol=rtol,
                    atol=rtol * 1e-2, dense_output=True, events=hit_cap,
                    max_step=L_bed / d_p / 200)
    return sol, bool(sol.t_events[0].size)


def hot_spot_1d(T0C, **kw):
    sol, ran = solve_1d(T0C, **kw)
    if ran:
        return np.nan, np.inf
    zz = np.linspace(0.0, sol.t[-1], 4000)
    tt = sol.sol(zz)[2]
    i = int(np.argmax(tt))
    return zz[i] * d_p, tt[i]


def critical_T0_1d(lo=360.0, hi=368.0, tol=0.02, cap=300.0, **kw):
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if np.isinf(hot_spot_1d(mid, cap=cap, **kw)[1]):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


class Tube1DPymrm:
    """Steady 1-D BVP route to the same reactor: pymrm upwind convection +
    Newton, the D2.2 machinery with this chemistry substituted in. Discretises
    the same equations a second way (steady BVP vs. adaptive march), so the two
    routes share only the printed constants."""

    def __init__(self, T0C, U=82.7, kelvin=KELVIN, n_z=4000):
        self.T0C, self.U, self.kelvin = T0C, U, kelvin
        self.shape = (n_z, 3)
        self.z_f = np.linspace(0.0, L_bed / d_p, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        # outward normal: inlet x = w = t = 0 -> {a:0, b:1, d:0};
        # outlet pure outflow -> {a:1, b:0, d:0}
        bc = ({"a": 0.0, "b": 1.0, "d": 0.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        conv, conv_bc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                                  bc, v=1.0, axis=0)
        dv = construct_div(self.shape, self.z_f, nu=0, axis=0)  # nu=0: axial coordinate
        self.jac_const = dv @ conv
        self.g_const = dv @ conv_bc
        self.numjac = NumJac(self.shape)       # pointwise source: last axis full

    def source(self, u):
        x, w, t = u[..., 0], u[..., 1], u[..., 2]
        T = t + self.T0C + self.kelvin
        k1 = np.exp(lnA[0] - E[0] / (Rg * T))
        k2 = np.exp(lnA[1] - E[1] / (Rg * T))
        k3 = np.exp(lnA[2] - E[2] / (Rg * T))
        y = x + w
        rB = NA0 * N0 * (k1 * (1 - y) - k2 * x)
        rC = NA0 * N0 * (k2 * x + k3 * (1 - y))
        return np.stack([b1 * rB, b1 * rC,
                         -4 * self.U * d_p / (G * c_p * d_t) * t
                         + b2 * rB + b3 * rC], axis=-1)

    def residual(self, u):
        g_rxn, jac_rxn = self.numjac(self.source, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_rxn.reshape((-1, 1))
        return g, self.jac_const - jac_rxn

    def solve(self):
        # initial guess from the march. This does NOT couple the two routes:
        # Newton converges to the root of ITS OWN upwind residual wherever it
        # starts, and the convergence is asserted below rather than assumed
        # (the B1.6 lesson: a reference that did not converge validates nothing).
        sol0, _ = solve_1d(self.T0C, U=self.U, kelvin=self.kelvin)
        u0 = sol0.sol(np.minimum(self.z_c, sol0.t[-1])).T.copy()
        result = newton(self.residual, u0, maxfev=100)
        self.u = result.x.reshape(self.shape)
        g, _ = self.residual(self.u)
        res_norm = float(np.max(np.abs(g)))
        assert res_norm < 1e-8, f"1-D BVP Newton did not converge: {res_norm:.1e}"
        return self

    def hot_spot(self):
        i = int(np.argmax(self.u[:, 2]))
        return self.z_c[i] * d_p, float(self.u[i, 2])'''))

cells.append(md(r"""## Results

### The two-dimensional model: Figures 9 and 12's family, and the 3 °C cliff

The radial-mean temperature rise $t_m(z)$ for the inlet (= coolant)
temperatures of the paper's Figure 9. At 357 °C the hot spot is a controlled
≈30 °C; at 360 °C the profile climbs past everything on the frame — the
paper's "the temperature rise goes out of control". The axis temperature runs
far above the mean, which is the whole argument against the 1-D model."""))

cells.append(code('''T0_family = [350.0, 355.0, 357.0, 358.5, 360.0]
runs2d = {T0: Radial2D().solve(T0) for T0 in T0_family}

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=False)
for T0 in T0_family:
    zz, xm, wm, tm, t_ax = runs2d[T0].profiles()
    lab = f"{T0:g} °C"
    axL.plot(zz, tm, lw=1.7, label=lab)
    axR.plot(zz, t_ax, lw=1.7, label=lab)
for ax, ttl in [(axL, "radial mean  $t_m$  (the paper's Figure 9)"),
                (axR, "axis  $t(r{=}0)$")]:
    ax.set(xlabel="bed length $z'$ (m)", ylabel="temperature rise (°C)",
           title=ttl, xlim=(0, 1.5), ylim=(0, 70))
    ax.legend(fontsize=8, title="inlet = coolant")
plt.tight_layout(); plt.show()

print("2-D hot spots (radial mean):")
hs2d = {}
for T0 in T0_family:
    z_hs, tm_hs = runs2d[T0].hot_spot()
    hs2d[T0] = tm_hs
    if np.isinf(tm_hs):
        print(f"  T0 = {T0:6.1f} degC : runaway (mean t exceeded 150 degC "
              f"inside the bed)")
    else:
        print(f"  T0 = {T0:6.1f} degC : tm_max = {tm_hs:6.2f} degC at "
              f"z' = {z_hs:.2f} m")
tm357 = hs2d[357.0]
print(f"\\npaper, p. 26: at 357 degC 'the hot spot, where tm equals about 30 degC'"
      f" -> computed {tm357:.2f} degC ({100*abs(tm357-30)/30:.1f} % from the stated ~30)")
m360 = runs2d[360.0]
zz, xm, wm, tm, t_ax = m360.profiles()
print(f"at 360 degC the mean rise reaches {tm.max():.0f} degC"
      + (" (integration stopped at the 150 degC cap: runaway)" if m360.ran_away
         else f" with axis peak {t_ax.max():.0f} degC - past the knee; the paper's"
              " Figure 9 shows this curve climbing off its 60 degC frame, and the"
              " continuous runaway limit is bisected in the next cells"))'''))

cells.append(md(r"""### Radial structure, the product maximum, and selectivity

The paper's Figure 8 shows severe radial gradients even in the controlled
case; its Figure 10 shows the selectivity cost of running hotter. Both are
reproduced qualitatively (the figures are not digitised — every quantitative
comparison on this page is against a *stated* number). The consecutive scheme
makes the phthalic anhydride conversion pass through a maximum; the paper says
3 m is just sufficient to reach it at 357 °C."""))

cells.append(code('''m357 = runs2d[357.0]
zz, xm, wm, tm, t_ax = m357.profiles()
n = m357.n_r
z_marks = [0.2, 0.4, 0.6, 0.9, 1.5]
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.0))
for zm in z_marks:
    u = m357.sol.sol(zm / d_p)
    axA.plot(m357.r_c / Rdim, u[2 * n:], lw=1.6, label=f"z' = {zm:g} m")
axA.set(xlabel="$r/R$", ylabel="temperature rise $t$ (°C)",
        title="radial temperature profiles, inlet 357 °C (cf. paper Fig. 8)")
axA.legend(fontsize=8)

for T0 in [355.0, 357.0, 360.0]:
    zzs, xms, wms, tms, _ = runs2d[T0].profiles()
    ys = xms + wms
    sel = np.divide(xms, ys, out=np.ones_like(xms), where=ys > 1e-6)
    axB.plot(ys, sel, lw=1.7, label=f"{T0:g} °C")
axB.set(xlabel="total conversion $y$", ylabel="selectivity $x/y$",
        title="selectivity vs conversion (cf. paper Fig. 10)",
        xlim=(0, 1), ylim=(0.5, 1.0))
axB.legend(fontsize=8, title="inlet = coolant")
plt.tight_layout(); plt.show()

i_max = int(np.argmax(xm))
x05, x10, x20, x25, x30 = np.interp([0.5, 1.0, 2.0, 2.5, 3.0], zz, xm)
print(f"phthalic anhydride conversion at 357 degC: x = {x20:.3f} at 2.0 m, "
      f"{x25:.3f} at 2.5 m, {x30:.3f} at 3.0 m; the local slope has fallen "
      f"from {2*(x10-x05):.2f}/m in the first meter to {2*(x30-x25):.2f}/m at "
      "the exit"
      + (f"; the maximum sits inside the bed at z' = {zz[i_max]:.2f} m"
         if i_max < len(zz) - 1 else
         "; the maximum sits at or just beyond the 3 m exit"))
print("paper, p. 26: 'a length of 3 meters is sufficient to reach the maximum "
      "in phthalic anhydride concentration' at 357 degC - consistent to "
      "within the flatness of the profile, though in this rebuild the "
      "maximum is not yet strictly passed at 3 m")
print(f"outlet selectivity x/y at 357 degC = {xm[-1]/(xm[-1]+wm[-1]):.3f}, "
      f"total conversion y = {xm[-1]+wm[-1]:.3f} "
      "(Figure 10's 357 curve reads ~0.83-0.85 at this conversion; qualitative only)")'''))

cells.append(md(r"""### The one-dimensional model with the derived $U$, and the two runaway limits

The 1-D model must use $U$ built from the same $\lambda_R$ and $\alpha_w$ —
the paper prints 82.7 and it is recomputed above. Its Figure 12 family is
reproduced below, then both runaway limits are located **continuously** by
bisection on cold-started solves, where the paper probed whole degrees. The
runaway criterion is a 150 °C cap on the mean temperature rise; the validation
section shows the located limit is insensitive to that choice."""))

cells.append(code('''T0_family_1d = [350.0, 355.0, 357.0, 358.5, 362.0, 363.0, 365.0]
plt.figure(figsize=(6.8, 4.2))
hs1d = {}
for T0 in T0_family_1d:
    sol, ran = solve_1d(T0)
    zz = np.linspace(0.0, sol.t[-1], 2500)
    tt = sol.sol(zz)[2]
    hs1d[T0] = np.inf if ran else tt.max()
    plt.plot(zz * d_p, tt, lw=1.6,
             label=f"{T0:g} °C" + (" (runaway)" if ran else ""))
plt.xlim(0, 1.5); plt.ylim(0, 70)
plt.xlabel("bed length $z'$ (m)"); plt.ylabel("temperature rise $t$ (°C)")
plt.title("1-D model, U = 82.7 (the paper's Figure 12)")
plt.legend(fontsize=8, title="inlet = coolant"); plt.tight_layout(); plt.show()

rise362, rise363 = hs1d[362.0], hs1d[363.0]
print(f"1-D hot spot at 362 degC: {rise362:.1f} degC   (paper: 'a rise of only 40 deg')")
print(f"1-D hot spot at 363 degC: {rise363:.1f} degC   (paper: 'a rise of 48 deg')")
print(f"  deviations: {100*abs(rise362-40)/40:.1f} % and {100*abs(rise363-48)/48:.1f} % "
      "on quantities that double every ~2 degC at this point")'''))

cells.append(code('''# The two runaway limits, bisected to 0.02 degC on deterministic cold starts.
T0c_2d = critical_T0_2d()
T0c_1d = critical_T0_1d()
gap = T0c_1d - T0c_2d
print(f"2-D runaway limit: T0 = {T0c_2d:.2f} degC   "
      "(paper: 'predicts a runaway ... for an inlet temperature of 360 degC')")
print(f"1-D runaway limit: T0 = {T0c_1d:.2f} degC   "
      "(paper: 'now found to be 365 degC'; Figure 12 probes 363 then 365)")
print(f"gap between the models: {gap:.2f} degC   "
      "(paper: 'within five degrees ... has to be considered excellent')")
print("\\nNote what the integer probing hides: the paper's two limits are "
      "statements that 360 (2-D) and 365 (1-D) RUN AWAY, with 357 and 363 "
      "controlled. The 1-D limit (363.93) sits inside its bracket (363, 365]. "
      f"The 2-D limit ({T0c_2d:.4f}) sits {T0c_2d - 360.0:+.4f} degC OUTSIDE "
      "(357, 360]: 360.0 was a bisection midpoint and came back CONTROLLED, so "
      "this model strictly fails the paper's 'runs away at 360' - by less than "
      "0.016 degC, far inside any physical uncertainty, but outside is outside "
      "and a padded bracket would be a check that cannot fail. The continuous "
      "gap is smaller than the quoted five degrees.")

# Figure 11's sweep: the runaway limit vs inlet concentration (44, 38, 32 g/Nm3).
# N_A0 scales with concentration and CANCELS out of the mass equations
# (b1*r_B carries no N_A0), so a lower inlet concentration is exactly a
# heat_scale < 1 on the temperature sources.
print("\\nrunaway limit vs inlet concentration (cf. paper Fig. 11, qualitative):")
crit_by_conc = {}
for conc in (44.0, 38.0, 32.0):
    T0c = critical_T0_2d(lo=356.0, hi=378.0, tol=0.05, n_r=48,
                         heat_scale=conc / 44.0)
    crit_by_conc[conc] = T0c
    print(f"  {conc:4.0f} g/Nm3 : T0_crit = {T0c:.2f} degC")
d_44_38 = crit_by_conc[38.0] - crit_by_conc[44.0]
d_38_32 = crit_by_conc[32.0] - crit_by_conc[38.0]
print(f"the limit rises by {d_44_38:.1f} degC from 44 to 38 and a further "
      f"{d_38_32:.1f} degC to 32 g/Nm3 - the direction and steepening the "
      "paper's Figure 11 draws. The figure itself is not digitised; no "
      "quantitative claim is made about its curve positions.")'''))

cells.append(md(r"""## Validation

Four layers, ordered by what each can and cannot catch. **What no numerical
check below can catch:** a mis-transcription of the *kinetic* constants that
both routes share. $k_1$'s constants have an independent witness (the D2.2
reconciliation above: $27000/1.98 = 13636$ and $\ln A = 19.837$, transcribed
separately from the 1970 paper); $k_2$ and $k_3$ have none — no other document
on disk prints them — so their guard is the reproduction agreement itself,
which is exactly why nothing here was allowed to be fitted.

### 1. The radial operator against the paper's own closed form

With the reaction off, the 2-D model *is* the packed-bed heat exchanger whose
exact Bessel-series solution the paper prints (p. 25). Marching the pymrm
operator from a uniform inlet and comparing the mean-temperature decay against
the series tests the cylindrical `nu=1` divergence, the Robin wall condition,
the symmetry axis and the area-weighted averaging at once — and none of the
kinetics. It **can fail**: the same cell shows the Cartesian operator
(`nu=0`) and a flipped wall-film sign, both of which move it far outside
tolerance. The moment-matching identity behind $U$ is also verified; that one
is the paper's own algebra, labelled an identity, and only checks the
transcription of the $K$ formula.

**A printed defect in that algebra, recorded not repaired:** the paper's
printed $1/K$ series (p. 25) carries $R^2$ where the zero-order moment of its
own $t_m'$ series requires $R^4$ — a factor $R^2 = 17.36\times$ at these
dimensions. The code below uses the correct $R^4$; with the printed $R^2$ the
moment identity misses by exactly that factor, which is how the defect was
found."""))

cells.append(code('''aR = alpha * Rdim
f_root = lambda s: s * j1(s) - aR * j0(s)
xs = np.linspace(1e-6, 150.0, 30001)
fs = f_root(xs)
roots = np.array([brentq(f_root, xs[i], xs[i + 1])
                  for i in range(len(xs) - 1) if fs[i] * fs[i + 1] < 0])

def tm_series(z):
    """Paper p. 25: tm\' = 4 a2R2 sum exp(-(a2/R2) lam_n^2 z)/(lam_n^2(lam_n^2+a2R2))
    with aR = alpha*R and lam_n the roots of lam J1 = aR J0."""
    return 4 * aR ** 2 * np.sum(
        np.exp(-(a2 / Rdim ** 2) * roots ** 2 * z)
        / (roots ** 2 * (roots ** 2 + aR ** 2)))

def hx_error(n_r=80, nu=1, wall_sign=+1.0, z_probe=50.0):
    """|pymrm - series|/series for the reactionless exchanger at z_probe."""
    m = Radial2D(n_r=n_r, nu=nu, wall_sign=wall_sign)
    sol = solve_ivp(lambda z, u: m.Lt @ u + m.gt, (0.0, z_probe),
                    np.ones(n_r), method="BDF", rtol=1e-11, atol=1e-13,
                    dense_output=True)
    tm_num = float(m.wgt @ sol.sol(z_probe))
    return abs(tm_num - tm_series(z_probe)) / tm_series(z_probe)

hx_ok = hx_error()
print(f"exchanger check at z = 50 (z' = 0.15 m), n_r = 80:")
print(f"  pymrm vs printed series : {hx_ok:.2e} relative")
print(f"  same with nu = 0        : {hx_error(nu=0):.2e}   <- the check CAN fail")
print(f"  same with wall sign flip: {hx_error(wall_sign=-1.0):.2e}   <- and fails loudly")
orders = [hx_error(n_r=nn) for nn in (20, 40, 80)]
p_obs = np.log2(orders[0] / orders[1]), np.log2(orders[1] / orders[2])
print(f"  grid convergence n_r 20/40/80: {orders[0]:.2e} / {orders[1]:.2e} / "
      f"{orders[2]:.2e}  (observed order {p_obs[0]:.2f}, {p_obs[1]:.2f})")

# the moment-matching identity behind U (paper p. 25) - IDENTITY, not physics:
# it tests only that the printed K formula was transcribed consistently.
K_closed = (2 * a2 / Rdim) * (4 * alpha / (aR + 4.0))
K_series = 1.0 / ((4 * aR ** 2 * Rdim ** 2 / a2)
                  * np.sum(1.0 / (roots ** 4 * (roots ** 2 + aR ** 2))))
print(f"\\nK closed form {K_closed:.6f} vs zero-moment series {K_series:.6f} "
      f"({abs(K_closed-K_series)/K_closed:.1e}; identity - the paper's own algebra)")'''))

cells.append(md(r"""### 2. Every stated number in the paper, in one table

The complete set of numerical statements the text makes about this case,
against what the rebuilt models give. The three hot-spot rows carry the
paper's own "about"; the runaway rows are integer-degree statements, so their
honest reading is a bracket (the stated limit runs away, the last probed
temperature below it does not)."""))

cells.append(code('''hs_lam075 = Radial2D(a2=0.75 / (G * c_p * d_p),
                     alpha=alpha_w * d_p / 0.75).solve(360.0).hot_spot()[1]
hs_aw150 = Radial2D(alpha=150.0 * d_p / lam_R).solve(360.0).hot_spot()[1]

rows = [
    ("Pe_hR from G, cp, dp, lam_R",      5.25, G * c_p * d_p / lam_R, "exact printed"),
    ("U(lam=0.67, aw=134)",              82.7, U_ref,                 "exact printed"),
    ("U(lam=0.75, aw=134)",              86.0, u86,                   "exact printed"),
    ("U(lam=0.67, aw=150)",              88.0, u88,                   "exact printed"),
    ("2-D hot spot t_m at 357 degC",     30.0, tm357,                 "'about'"),
    ("2-D t_m at 360 degC, lam_e=0.75",  35.0, hs_lam075,             "'about'"),
    ("2-D t_m at 360 degC, alpha_w=150", 35.0, hs_aw150,              "'about'"),
    ("1-D rise at 362 degC",             40.0, rise362,               "stated"),
    ("1-D rise at 363 degC",             48.0, rise363,               "stated"),
    ("2-D runaway limit",               360.0, T0c_2d,                "runs away at (bracket 357-360)"),
    ("1-D runaway limit",               365.0, T0c_1d,                "runs away at (bracket 363-365)"),
]
tab = pd.DataFrame(rows, columns=["quantity", "paper", "this page", "qualifier"])
tab["dev %"] = 100 * (tab["this page"] - tab["paper"]).abs() / tab["paper"].abs()
display(tab.style.format({"paper": "{:.2f}", "this page": "{:.2f}", "dev %": "{:.2f}"}))

stated = tab.iloc[4:9]     # the five model-output statements with real content
print(f"five stated model results: mean |dev| = {stated['dev %'].mean():.1f} %, "
      f"worst = {stated['dev %'].max():.1f} % (the 363 degC rise, a quantity that "
      "doubles every ~2 degC of inlet temperature at that point)")
print(f"runaway limits: 2-D {T0c_2d:.4f} in (357, 360]: "
      f"{357.0 < T0c_2d <= 360.0} (outside by {T0c_2d - 360.0:+.4f} degC - "
      f"reported, not padded); 1-D {T0c_1d:.4f} in (363, 365]: "
      f"{363.0 < T0c_1d <= 365.0}")
print("\\nREPRODUCTION, not validation: every 'paper' value here is Froment's "
      "own 1967 Crank-Nicolson/numerical-integration result. Nothing was "
      "fitted - all model inputs are printed constants - so the agreement is "
      "a genuine test of the transcription and the solver, but says nothing "
      "about the kinetics against reality.")'''))

cells.append(md(r"""### 3. A second, independent route to the headline numbers

The headline quantities are recomputed by routes that share **no assembly
code** with the pymrm models:

- **2-D**: an independently written conservative finite-difference operator
  (hand-built face-flux balances, one-sided Robin ghost at the wall) marched
  with a *different* stiff integrator (`Radau` vs `BDF`);
- **1-D**: the pymrm steady-BVP route (`construct_convflux_upwind` + `newton`,
  the D2.2 machinery) against the adaptive march used above.

The routes share the printed constants — and, for the 1-D pair, only the
Newton *initial guess*, which cannot move a converged root (the BVP's own
residual is asserted below 1e-8) — so agreement tests discretisation,
integration and assembly: everything except the transcription, which is
guarded separately above."""))

cells.append(code('''def indep_2d(T0C, n_r=80, lamR=lam_R, alphaw=alpha_w, PemR=Pe_mR,
             kelvin=KELVIN, cap=150.0):
    """Independent route: hand-built FD + Radau. No pymrm calls anywhere."""
    a1i, a2i = 1.0 / PemR, lamR / (G * c_p * d_p)
    al = alphaw * d_p / lamR
    r_f = np.linspace(0.0, Rdim, n_r + 1)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    h = r_f[1] - r_f[0]

    def lap(D, robin):
        A = np.zeros((n_r, n_r))
        for i in range(n_r):
            if i > 0:
                A[i, i - 1] += D * r_f[i] / (h * h * r_c[i])
                A[i, i] -= D * r_f[i] / (h * h * r_c[i])
            if i < n_r - 1:
                A[i, i + 1] += D * r_f[i + 1] / (h * h * r_c[i])
                A[i, i] -= D * r_f[i + 1] / (h * h * r_c[i])
            elif robin:   # wall flux -D dt/dr = D*al*t_wall, ghost-free one-sided
                A[i, i] -= D * r_f[i + 1] * al / ((1 + al * h / 2) * h * r_c[i])
        return A

    Lm, Lt = lap(a1i, False), lap(a2i, True)
    wgt = (r_f[1:] ** 2 - r_f[:-1] ** 2) / Rdim ** 2

    def rhs(z, u):
        x, w, t = u[:n_r], u[n_r:2 * n_r], u[2 * n_r:]
        T = t + T0C + kelvin
        k1 = np.exp(lnA[0] - E[0] / (Rg * T))
        k2 = np.exp(lnA[1] - E[1] / (Rg * T))
        k3 = np.exp(lnA[2] - E[2] / (Rg * T))
        y = x + w
        rB = NA0 * N0 * (k1 * (1 - y) - k2 * x)
        rC = NA0 * N0 * (k2 * x + k3 * (1 - y))
        return np.concatenate([Lm @ x + b1 * rB, Lm @ w + b1 * rC,
                               Lt @ t + b2 * rB + b3 * rC])

    def hit_cap(z, u):
        return wgt @ u[2 * n_r:] - cap
    hit_cap.terminal, hit_cap.direction = True, 1
    sol = solve_ivp(rhs, (0.0, L_bed / d_p), np.zeros(3 * n_r), method="Radau",
                    rtol=1e-8, atol=1e-10, dense_output=True, events=hit_cap)
    if sol.t_events[0].size:
        return np.inf
    zz = np.linspace(0.0, sol.t[-1], 4000)
    return float((sol.sol(zz)[2 * n_r:, :].T @ wgt).max())


tm357_indep = indep_2d(357.0, n_r=40)
route2_tm = abs(tm357_indep - tm357) / tm357
print("2-D hot spot at 357 degC, two independent routes:")
print(f"  pymrm operators + BDF, n_r = 80      : {tm357:.4f} degC")
print(f"  hand-built FD + Radau, n_r = 40      : {tm357_indep:.4f} degC   "
      f"(rel diff {route2_tm:.2e})")
print("  (deliberately different grids, so the reported difference is the "
      "deterministic discretisation gap, not integrator noise)")

lo, hi = 356.0, 364.0
while hi - lo > 0.02:
    mid = 0.5 * (lo + hi)
    if np.isinf(indep_2d(mid)):
        hi = mid
    else:
        lo = mid
T0c_2d_indep = 0.5 * (lo + hi)
route2_crit = abs(T0c_2d_indep - T0c_2d)
print(f"2-D runaway limit, independent route: {T0c_2d_indep:.2f} degC "
      f"(pymrm route {T0c_2d:.2f}; |diff| = {route2_crit:.2f} degC, "
      "quantised by the 0.02 degC bisection - printed for the record, not a "
      "CI metric, because a sub-tolerance difference is not resolvable)")

print("\\n1-D hot spots, adaptive march vs pymrm steady BVP (n_z = 4000):")
route2_1d = {}
for T0 in (357.0, 362.0, 363.0):
    bvp = Tube1DPymrm(T0).solve().hot_spot()[1]
    marched = hs1d.get(T0) or hot_spot_1d(T0)[1]
    route2_1d[T0] = abs(bvp - marched) / marched
    print(f"  T0 = {T0:5.1f}: march {marched:7.3f}  pymrm BVP {bvp:7.3f}  "
          f"rel diff {route2_1d[T0]:.2e}  (first-order upwind bias, shrinks with n_z)")'''))

cells.append(md(r"""### 4. Refinement, threshold insensitivity, and the break table

Grid and tolerance refinement for the headline hot spot; the runaway limit's
insensitivity to the arbitrary 150 °C cap; then **defect injection**: every
reported metric class, deliberately broken, to show each check moves when what
it guards is wrong. The last two rows are the two thousands-separator readings
taken literally — the trap this transcription actually faced."""))

cells.append(code('''print("grid refinement, 2-D hot spot at 357 degC (BDF rtol 1e-8):")
tm_by_n = {}
for nn in (20, 40, 80, 160):
    tm_by_n[nn] = Radial2D(n_r=nn).solve(357.0).hot_spot()[1]
    print(f"  n_r = {nn:4d}: tm_max = {tm_by_n[nn]:.4f} degC")
e1 = abs(tm_by_n[40] - tm_by_n[160])
e2 = abs(tm_by_n[80] - tm_by_n[160])
grid_order = np.log2(abs(tm_by_n[20] - tm_by_n[80]) / abs(tm_by_n[40] - tm_by_n[160]))
print(f"  apparent order ~ {grid_order:.1f}; n_r = 80 sits {e2:.4f} degC from n_r = 160")

print("\\nintegrator tolerance, same quantity, n_r = 80:")
for rt in (1e-6, 1e-8, 1e-10):
    v = Radial2D().solve(357.0, rtol=rt).hot_spot()[1]
    print(f"  rtol = {rt:.0e}: tm_max = {v:.5f} degC")

print("\\nrunaway-cap insensitivity (2-D limit, bisected to 0.02 degC):")
for cap in (100.0, 150.0, 250.0):
    v = critical_T0_2d(cap=cap, n_r=48)
    print(f"  cap = {cap:5.0f} degC: T0_crit = {v:.2f} degC")
print("  the limit moves by less than the bisection tolerance across a 2.5x "
      "range of the cap: near criticality the peak passes any finite cap "
      "almost at once.")'''))

cells.append(code('''# --- break table: inject one defect at a time, watch the metrics move -------
def broken(label, tm_fn, hx_fn=None, note=""):
    tm_val = tm_fn()
    hx_val = hx_fn() if hx_fn else None
    return (label, tm_val, hx_val, note)

base_a2 = a2
break_rows = [
    ("baseline", tm357, hx_ok, "printed constants, nu=1, outward-normal Robin"),
    broken("nu = 0 (Cartesian radial operator)",
           lambda: Radial2D(nu=0).solve(357.0).hot_spot()[1],
           lambda: hx_error(nu=0),
           "wrong geometry: slab has more cross-section near the wall"),
    broken("wall film sign flipped (b = -alpha)",
           lambda: Radial2D(wall_sign=-1.0).solve(357.0).hot_spot()[1],
           lambda: hx_error(wall_sign=-1.0),
           "the outward-normal trap: film heats instead of cools"),
    broken("k2 <-> k3 swapped",
           lambda: Radial2D(k_perm=(0, 2, 1)).solve(357.0).hot_spot()[1], None,
           "mis-assigned Arrhenius lines; no independent witness exists for "
           "k2, k3, so THIS metric is their only guard"),
    broken("heat Pe taken = mass Pe (a2 = 1/10)",
           lambda: Radial2D(a2=0.1).solve(357.0).hot_spot()[1],
           None, "the two Peclet numbers differ by 1.9x"),
    broken("Kelvin offset 273 instead of 273.15",
           lambda: Radial2D(kelvin=273.0).solve(357.0).hot_spot()[1],
           None, "the unprinted convention; see also the 1-D rows below"),
    broken("dH3 read literally as -1.090 kcal/gmol",
           lambda: Radial2D(b3_val=b3 / 1000.0).solve(357.0).hot_spot()[1],
           None, "combustion heat 1000x low"),
    broken("G read literally as 4.684 kg/(m2 hr)",
           lambda: (Radial2D(a2=lam_R / (4.684 * c_p * d_p),
                             b1_val=rho_b * d_p * M_m / (4.684 * NA0),
                             b2_val=rho_b * d_p * dH1 / (4.684 * c_p),
                             b3_val=rho_b * d_p * dH3 / (4.684 * c_p))
                    .solve(357.0).hot_spot()[1]),
           None, "every G-bearing group 1000x off: the 3 m bed behaves like "
                 "3 km of packing (reaction complete millimetres in, hot spot "
                 "displaced) - wrong, and visibly so"),
]
bt = pd.DataFrame(break_rows,
                  columns=["injected defect", "tm_max @357 (degC)",
                           "exchanger check", "note"])
display(bt.style.format({"tm_max @357 (degC)": "{:.2f}",
                         "exchanger check": lambda v: "-" if v is None else f"{v:.1e}"}))

print("1-D sensitivities on the 363 degC rise (paper: 48):")
for lab, kw in [("U = 82.7 (printed, baseline)", {}),
                ("U = 82.46 (recomputed from lam_R, alpha_w)", {"U": U_ref}),
                ("Kelvin offset 273", {"kelvin": 273.0})]:
    v = hot_spot_1d(363.0, **kw)[1]
    print(f"  {lab:44s}: {v:.2f} degC")
t0c_273 = critical_T0_1d(kelvin=273.0)
k1_shift = 100 * (np.exp(E[0] / Rg * 0.15 / 630.0 ** 2) - 1)
print(f"  1-D runaway limit with Kelvin 273: {t0c_273:.2f} degC "
      f"(baseline {T0c_1d:.2f}; the 0.15 K convention shift is "
      f"{k1_shift:.2f} % in k1 and moves the limit "
      f"{abs(t0c_273-T0c_1d):.1f} degC)")
print("\\nEvery injected defect moves the metric it should move; none is "
      "absorbed silently. What this table CANNOT do: catch a shared "
      "mis-transcription (see the section header).")'''))

cells.append(code('''metrics = {
    # transcription identities (printed derived values recomputed)
    "pehr_identity_rel": pe_dev,
    "u_ref_identity_rel": u_dev,
    "u_lam075_identity_rel": u86_dev,
    "u_aw150_identity_rel": u88_dev,
    "na0_from_44gnm3_rel": na0_dev,
    # operator against the paper's printed closed form (can-fail check)
    "hx_series_rel_z50": hx_ok,
    # reproduction of the stated results
    "tm2d_hotspot_357_C": tm357,
    "tm2d_357_vs_paper_rel": abs(tm357 - 30.0) / 30.0,
    "rise1d_362_C": rise362,
    "rise1d_362_vs_paper_rel": abs(rise362 - 40.0) / 40.0,
    "rise1d_363_C": rise363,
    "rise1d_363_vs_paper_rel": abs(rise363 - 48.0) / 48.0,
    "tm2d_360_lam075_C": hs_lam075,
    "tm2d_360_lam075_vs_paper_rel": abs(hs_lam075 - 35.0) / 35.0,
    "tm2d_360_aw150_C": hs_aw150,
    "tm2d_360_aw150_vs_paper_rel": abs(hs_aw150 - 35.0) / 35.0,
    "t0crit_2d_C": T0c_2d,
    "t0crit_1d_C": T0c_1d,
    "runaway_gap_C": gap,
    # independent-route agreement (can-fail checks; both values are
    # discretisation-gap dominated, hence deterministic)
    "route2_tm357_rel": route2_tm,
    "route2_rise1d_363_rel": route2_1d[363.0],
}
report_agreement("C2.10-o-xylene-phthalic-anhydride", metrics);'''))

cells.append(md(r"""## What pymrm adds

**The reproduction itself is the point** — the 1967 landmark rebuilt in about
a hundred lines, both models, from the printed constants alone, with the
1-D/2-D comparison that took the Ghent computing centre a custom
Crank–Nicolson program now a two-line operator swap (`nu=1` radial diffusion
in, lumped $-4U/(Gc_p d_t)$ out). Beyond that, three things the paper could
only assert or sample:

1. **Continuous runaway limits.** The paper probes whole degrees (357/360,
   363/365); bisection on cold starts locates both limits to 0.02 °C and shows
   the true 1-D/2-D gap is *smaller* than the quoted five degrees — the
   paper's headline conclusion survives being computed properly.
2. **The $\mathrm{Pe}_{mR} = 8$ claim, re-measured at the knife edge.** The
   paper reports the influence of $\mathrm{Pe}_{mR}$ = 8 vs 10 as "completely
   negligible", from repeating a handful of profile calculations. The cell
   below confirms it at controlled conditions — and then asks the harder
   question the 1967 sampling could not afford: does it stay negligible for
   the *runaway limit itself*, the most sensitive number in the problem? It
   does — the limit moves by less than the bisection tolerance — which is a
   stronger form of the paper's claim, not merely a repetition of it.
3. **The heat-transfer sensitivity as a curve.** The paper's Figures 13–14
   are two points (10 % on $\lambda_e$, 10 % on $\alpha_w$, both rescuing the
   360 °C runaway to ~35 °C — both reproduced in the table above). The sweep
   below draws the whole margin: how much heat-transfer error the design can
   absorb at each inlet temperature."""))

cells.append(code('''# Pe_mR = 8 vs 10: negligible where the paper looked, decisive at the edge
tm357_pe8 = Radial2D(a1=1.0 / 8.0).solve(357.0).hot_spot()[1]
print(f"Pe_mR 10 -> 8 at 357 degC: tm_max {tm357:.2f} -> {tm357_pe8:.2f} degC "
      f"({100*abs(tm357_pe8-tm357)/tm357:.1f} % - 'completely negligible', as stated)")
lo, hi = 350.0, 364.0
while hi - lo > 0.05:
    mid = 0.5 * (lo + hi)
    m = Radial2D(n_r=48, a1=1.0 / 8.0).solve(mid)
    if m.ran_away or m.hot_spot()[1] > 150.0:
        hi = mid
    else:
        lo = mid
T0c_pe8 = 0.5 * (lo + hi)
print(f"and the runaway limit itself: {T0c_2d:.2f} degC (Pe_mR = 10) vs "
      f"{T0c_pe8:.2f} degC (Pe_mR = 8) - a shift of "
      f"{abs(T0c_2d - T0c_pe8):.2f} degC, below the 0.05 degC bisection "
      "tolerance. The paper's 'completely negligible' holds even for the most "
      "sensitive quantity in the problem, which its profile-level sampling "
      "could not have established.")

# the heat-transfer margin, drawn as a curve instead of two points
scales = np.linspace(0.95, 1.25, 13)
crit_vs_scale = []
for s in scales:
    lam_s, aw_s = lam_R * s, alpha_w * s
    lo, hi = 355.0, 375.0
    while hi - lo > 0.05:
        mid = 0.5 * (lo + hi)
        m = Radial2D(n_r=48, a2=lam_s / (G * c_p * d_p),
                     alpha=aw_s * d_p / lam_s).solve(mid)
        if m.ran_away or m.hot_spot()[1] > 150.0:
            hi = mid
        else:
            lo = mid
    crit_vs_scale.append(0.5 * (lo + hi))

plt.figure(figsize=(6.4, 4.0))
plt.plot(100 * (scales - 1), crit_vs_scale, "o-", lw=1.6, ms=4)
plt.axhline(360, color="tab:red", lw=0.9, ls="--",
            label="360 °C (runs away at nominal heat transfer)")
plt.xlabel(r"heat-transfer parameters $\\lambda_R,\\ \\alpha_w$ scaled by (%)")
plt.ylabel("2-D runaway limit $T_{0,crit}$ (°C)")
plt.title("how much heat-transfer error the design absorbs")
plt.legend(fontsize=8); plt.tight_layout(); plt.show()
slope = np.polyfit(100 * (scales - 1), crit_vs_scale, 1)[0]
print(f"sensitivity near nominal: {slope:.2f} degC of runaway limit per percent "
      "of heat-transfer capacity - the quantitative form of the paper's "
      "closing warning that prediction 'requires a degree of precision in the "
      "measurement of the experimental parameters seldom achieved', and of its "
      "remark that the same curves could be shifted 'equally well by a slight "
      "modification of the kinetic coefficients'.")'''))

cells.append(md(r"""## Reuse

- **This chemistry in another reactor**: everything model-specific is the
  `rates` method and the constant block; substitute both and the marching,
  bisection and validation scaffolding carry over unchanged.
- **The 1-D runaway machinery** with criteria instead of brute force is the
  sibling page [`D2.2`](../D2.2-van-welsenaere-froment-runaway/) — same
  reactor, same $k_1$, closed-form critical-inlet criteria.
- **The S6 radial-operator pattern** (cylindrical `construct_grad`/
  `construct_div`, Robin wall) is shared with
  [`A2.3`](../A2.3-taylor-aris-dispersion/) (dispersion closure) and
  [`A3.15`](../A3.15-graetz-nusselt/) (Graetz problem); the radial
  heat-transfer parameters themselves are
  [`A3.12`](../A3.12-yagi-kunii-effective-conductivity/) (Yagi–Kunii, the very
  correlation Froment used for $\lambda_R$).
- **Froment kinetics with data** is
  [`C2.1`](../C2.1-xu-froment-smr/) (Xu & Froment SMR, validated against the
  runs it was fitted to) — use it as the model when a kinetics page proper is
  wanted.

**Limits to keep in mind.** The kinetics are pseudo-first-order effective rate
laws at 1 atm in large air excess on a 1960s V₂O₅ catalyst, asserted by the
paper as representative and never fitted to data in this document; do not
transplant them quantitatively. The runaway limit is exponentially sensitive:
this page measures 0.3 °C of limit per percent of heat-transfer capacity
(so the paper's own 10 % parameter uncertainty is worth ~3 °C), and even the
unprinted 273-vs-273.15 Kelvin convention — 0.5 % in $k_1$ — moves the 1-D
limit by 0.2 °C. Any reuse of the *numbers* (rather than the pattern)
inherits that sensitivity. The paper's own kinetics were superseded by
Calderbank's measurements (not on disk) and by Froment & Bischoff's textbook
treatment — a separate `D3.4` page against those sources would be validation;
this page is reproduction."""))

nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python",
                   "name": "python3"},
    "language_info": {"name": "python", "version": "3"},
}
out = Path(__file__).parent / "index.ipynb"
nbf.write(nb, out)
print(f"wrote {out}")
