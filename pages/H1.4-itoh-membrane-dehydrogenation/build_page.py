#!/usr/bin/env python3
"""Generate index.ipynb for page H1.4. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "Itoh Pd membrane reactor: shifting a dehydrogenation equilibrium"
description: "Hydrogen leaves through the wall, so a reaction capped at 18.7 % conversion runs to 99.7 %. The model puts that operating point on its thermodynamic ceiling and lands on the measured number — and the page quantifies how little that comparison can actually test."
categories: [sec:H, struct:S7, tier:T0, data:tier2, phase:gas-solid]
date: 2026-07-31
---

# Itoh Pd membrane reactor: shifting a dehydrogenation equilibrium

**Catalog ID:** `H1.4` · **Structures:** `S7` (two coupled 1-D domains) · **Tier:** T0

At 473 K and 1 atm, cyclohexane diluted in argon dehydrogenates to an
equilibrium conversion of **18.7 %** — and no catalyst, however good, can do
better in an ordinary tube. Itoh packed the catalyst inside a palladium tube
instead. Palladium passes hydrogen and *nothing else*, so the product hydrogen
leaves through the wall as fast as it forms, the reverse reaction is starved,
and the measured conversion reaches **99.7 %**.

This page builds the paper's model with nothing fitted and lands on that
number. It is equally careful about what the comparison can establish: at the
measured operating point the model sits **exactly on its own thermodynamic
ceiling**, so the agreement tests the equilibrium-plus-permeation description
and is completely insensitive to the kinetics and to the permeation constant.
The sensitivity sweep that shows this is printed in the validation section, and
the summary lines on this page are written to match it.

Two printed statements in the 3-page communication cannot be right as they
stand, and both are corrected below using only other printed statements of the
same paper. A third value that a reader might suspect — the headline purge rate
— turns out to be printed correctly, and is *confirmed* here rather than
corrected."""))

cells.append(md(r"""## Background

A reversible reaction stops at equilibrium; removing a product moves the
equilibrium. Raymont (1975) proposed exploiting this for H₂S decomposition, and
Kameyama et al. (1981) tried it with a porous Vycor glass membrane — but a
porous membrane passes *everything* to some degree, and Itoh's own simulations
(1984, 1985) showed that leakage of the feed caps what such a reactor can do.

A dense palladium membrane has no such cap. Hydrogen dissolves in the metal
(dissociatively — hence the square-root pressure law), diffuses through, and
recombines on the far side; cyclohexane and benzene cannot cross at all. The
1987 communication is the experimental demonstration: a double-tube reactor,
0.5 wt.% Pt/Al₂O₃ pellets packed inside a 200 µm palladium tube, cyclohexane in
argon fed through the bed, and an argon purge sweeping the shell side to carry
the permeated hydrogen away. Neither cyclohexane nor benzene was detected on
the shell side — the selectivity really is absolute.

The model reaction:

$$
\mathrm{C_6H_{12}} \;\rightleftharpoons\; \mathrm{C_6H_6} + 3\,\mathrm{H_2},
\qquad K_p = \frac{p_B\,p_H^3}{p_C} .
$$"""))

cells.append(md(r"""## The published model

Isothermal, isobaric plug flow on both sides, co-current, with molar flow rates
as the variables: $u_i$ on the reaction side, $v_i$ on the shell side, against
dimensionless length $L$. Hydrogen permeation obeys the half-power (Sieverts)
law, Eq. 1:

$$
Q_H = \alpha_H\left(\sqrt{p_H/P_0} - \sqrt{p'_H/P_0}\right),
\qquad
\alpha_H = \frac{2\pi l_0}{\ln(r_o/r_i)}\,D\,C_0 ,
$$

and the species balances (Eqs. 2–3) are

$$
\frac{du_C}{dL} = r_C V_r, \qquad
\frac{du_H}{dL} = -3\,r_C V_r
  - \alpha_H\!\left(\sqrt{\tfrac{u_H}{\sum u_i}} - \sqrt{\tfrac{v_H}{\sum v_i}}\right),
$$

(the square roots are of mole fractions because $P_{T_r} = P_{T_s} = P_0$ =
1 atm), with benzene and shell hydrogen closed algebraically (Eqs. 4–5), argon
constant on both sides (Eq. 6), and the Langmuir–Hinshelwood rate of Itoh et
al. (1985), Eq. 7:

$$
r_C = \frac{-k\,(K_p\,p_C/p_H^3 - p_B)}{1 + K_B K_p\,p_C/p_H^3},
\qquad
\begin{aligned}
k   &= 0.221\,e^{-4270/T} &&\mathrm{mol\,m^{-3}\,Pa^{-1}\,s^{-1}}\\
K_B &= 2.03\times10^{-10}\,e^{6270/T} &&\mathrm{Pa^{-1}}
\end{aligned}
$$

Inlet conditions (Eq. 8): everything enters at $L=0$; $u_B = u_H = v_H = 0$
there. The authors integrated Eqs. 2–7 by Runge–Kutta–Gill.

### Two misprints, and one printed value confirmed

**A note on how the paper was read.** The PDF's page images are 300 dpi bilevel
CCITT scans (`pdfimages -list`), so rendering at 600 dpi is a 2× upsample that
carries no extra information. Every marginal glyph quoted below was therefore
read on the *native* bitmap; the crops are kept with the case file.

Two printed statements cannot be right as they stand. Each is resolved **using
only other printed statements of the same paper** — nothing is fitted and
nothing is guessed.

**Misprint 1 — Eq. 5's sign.** The paper prints
$v_H = u_H - 3(u_C^0 - u_C)$. Integrating the paper's own Eqs. 2 and 3 together
with $dv_H/dL = Q_H$ gives the hydrogen balance
$u_H + v_H = 3(u_C^0 - u_C)$; substituting that into the printed Eq. 5 gives
$v_H = -v_H$, so **the printed equation forces $v_H \equiv 0$** — no hydrogen
in the permeate at all, which contradicts the permeation term of Eq. 3 outright.
The correct form is $v_H = 3(u_C^0 - u_C) - u_H$: a sign flip in print.

There is no sign-convention escape. Eq. 8 places *every* initial condition at
$L = 0$ ($u_H = 0$, $v_H = 0$ and $v_A = v_A^0$ there); Figure 3 draws $v_A^0$
entering the separation side at the left, alongside $u_C^0$ and $u_A^0$;
Figure 1 shows the purge gas entering at the feed end and leaving with the
product; and the sentence below Eq. 6 defines $u_i$ and $v_i$ as "the flow rates
of component $i$ in the reaction and separation sides", i.e. positive in $+L$.
The reactor is co-current in the apparatus, in the flow model and in the initial
conditions, so $v_H \geq 0$ throughout. Nothing computed here depends on the
call: the notebook integrates $dv_H/dL = +Q_H$ and uses the corrected identity
only as a check.

**Misprint 2 — the equilibrium constant.** The printed
$K_P = 4.89\times10^{35}\,e^{+3190/T}$ Pa³ gives $4.2\times10^{38}$ Pa³ at
473 K. Fed into Eq. 7's own equilibrium condition $K_p\,p_C/p_H^3 = p_B$, that
predicts essentially **complete** equilibrium conversion — contradicting the
paper's stated 18.7 %, the printed "Equilibrium conversion (0.187)" intercept
on Figure 4, and the premise of the entire study. Two independent checks show
this is not a reading error on our side:

* **A units error cannot explain it, in either direction.** Eq. 7 forces
  $[K_p] = \mathrm{Pa^3}$ ($K_p p_C/p_H^3$ is compared with $p_B$, and
  $K_B K_p p_C/p_H^3$ must be dimensionless given the printed
  $[K_B] = \mathrm{Pa^{-1}}$). Reading the printed number as atm³ instead makes
  it $1.04\times10^{15}$ times *larger* in Pa³ (table below); and expressed in
  atm³ the printed value is $4\times10^{23}$ atm³ against the
  $2.3\times10^{-4}$ atm³ the paper's own $x_{eq}$ requires — 27 orders out
  whichever way the conversion is applied.
* **The sign of the exponent is independently wrong.** Cyclohexane
  dehydrogenation is endothermic ($\Delta H^\circ \approx +206$ kJ/mol from the
  standard formation data used below), so $K_p$ must *rise* with temperature and
  the exponential argument must be negative. The printed $+3190/T$ describes an
  exothermic reaction.

The decision table below evaluates the printed expression and its candidate
repairs against the paper's own $x_{eq}$; none reproduces it. So
$K_p(473\,\mathrm{K})$ is **reconstructed from the stated equilibrium
conversion**: $x_{eq} = 0.187$ at $y_C^0 = 0.197$, 1 atm — all printed numbers.
The model runs at a single temperature, so the Arrhenius form is never needed.

**Confirmed, not corrected — the headline purge rate $v_A^0$.** The text prints
$v_A^0 = 11.8\times10^{-5}$ mol/s, and on the native bitmap that superscript
reads $-5$ cleanly (the same line carries an unambiguous $-7$ for comparison).
**This page therefore does not correct the paper here**; it confirms the printed
value, because everything else on the page would change if the exponent were
$-6$ and the confirmation happens to be available without any kinetics: with
$v_A^0 = 11.8\times10^{-6}$ mol/s, *no* kinetics and *no* permeance — both taken
infinitely fast — could push conversion past $X \approx 0.72$, so the measured
0.997 rules that reading out (numbers below). That ceiling is specific to the
**co-current** arrangement — a counter-current purge of the same size would do
much better — which is exactly why the three independent statements of
co-currency quoted above matter. What *is* illegible on the scan is
Figure 4's hand-lettered abscissa exponent ($v_A^0 \times 10^5$ or $10^6$);
consistency with the text makes it $\times 10^5$, and nothing on this page
depends on reading it.

**Deviation convention, used everywhere on this page:**
(model − measured)/measured."""))

cells.append(md(r"""## Parameters and assumptions

| Quantity | Value | Source |
|---|---|---|
| $T$ | 473 K | stated |
| $P_{T_r} = P_{T_s} = P_0$ | 1.013 × 10⁵ Pa | stated (isobaric, atmospheric) |
| Feed | 19.7 mol% cyclohexane in argon | stated |
| Tube | $r_i$ = 8.5 mm, $r_o$ = 8.7 mm, $l_0$ = 0.14 m | Figure 1 + α_H identity (below) |
| $D$ (H in Pd) | 9.23 × 10⁻¹⁰ m²/s | stated (Nagamoto & Inoue) |
| $C_0$ | 1280 mol/m³ | stated (Sieverts & Danz) |
| $\alpha_H$ | 4.47 × 10⁻⁵ mol/s | stated; recomputed to 0.03 % below |
| $k(473)$ | 2.653 × 10⁻⁵ mol m⁻³ Pa⁻¹ s⁻¹ | Eq. 7 |
| $K_B(473)$ | 1.160 × 10⁻⁴ Pa⁻¹ | Eq. 7 |
| $K_p(473)$ | 2.357 × 10¹¹ Pa³ | reconstructed from $x_{eq}=0.187$ (misprint 2) |
| $V_r$ | $\pi r_i^2 l_0$ = 3.18 × 10⁻⁵ m³ | "gross volume of the reaction section" |

**Note on the radii, and how little it matters.** The text says "17.0 mm OD",
but recomputing $\alpha_H$ from Eq. 1 reproduces the printed
4.47 × 10⁻⁵ mol/s only with $r_i$ = 8.5 mm, $r_o$ = 8.7 mm (0.03 % off; the OD
reading gives 2.4 %). Strictly, that identity shows what Itoh *substituted into
Eq. 1*, not what the tube physically was — "17.0 mm OD" may simply have been
mis-stated, or $\alpha_H$ computed with the wrong radii. Either way the model
uses the *printed* $\alpha_H$, and the only quantity that consumes $r_i$ is
$V_r$ (4.9 % between the two readings). $V_r$ enters the balances only as the
product $k\,V_r$, so the ×0.1–×10 variation of $k$ in the sensitivity sweep
below covers it many times over — and that sweep moves no number reported on
this page at all. The identity is recorded because it is a genuine internal
check on Eq. 1, not because anything hangs on it.

**Assumptions inherited from the paper:** isothermal, isobaric, plug flow on
both sides, negligible reaction on the palladium surface itself, Sieverts-law
permeation, and no mass-transfer resistance between bed and wall."""))

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
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from pymrm import construct_convflux_upwind, construct_div, NumJac
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "H1.4-itoh-membrane-dehydrogenation"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

# --- constants printed in the paper (all read off the native 300 dpi bitmap) ---
T = 473.0                        # K
P0 = 1.013e5                     # Pa, = P_Tr = P_Ts (isobaric, 1 atm)
Y_C0 = 0.197                     # feed mole fraction of cyclohexane in argon
L0 = 0.14                        # m, tube length
RI, RO = 8.5e-3, 8.7e-3          # m, tube radii (selected by the alpha_H identity)
D_PD = 9.23e-10                  # m2/s, H in Pd
C_0 = 1280.0                     # mol/m3, Sieverts reference concentration
ALPHA_H = 4.47e-5                # mol/s, printed permeation constant (Eq. 1)
K_RATE = 0.221 * np.exp(-4270.0 / T)      # mol m-3 Pa-1 s-1
K_B = 2.03e-10 * np.exp(6270.0 / T)       # Pa-1
V_R = np.pi * RI**2 * L0                  # m3, gross volume of the reaction section
N_AR = (1.0 - Y_C0) / Y_C0                # mol Ar fed per mol cyclohexane'''))

cells.append(md(r"""## The data

A 3-page communication with no tables: the validation targets live in the
running text. The one **measurement** is the headline conversion, 0.997, at
$u_C^0 = 2.90\times10^{-7}$ mol/s and $v_A^0 = 11.8\times10^{-5}$ mol/s; the
one **thermodynamic statement** is the equilibrium conversion, 0.187, for the
same feed without a membrane.

The measured conversion is quoted to three digits with no stated uncertainty,
so it means $X \in [0.9965, 0.9975]$ at best — which already bounds how finely
any deviation from it can be quoted. The validation section works that through.

Figure 4 also carries experimental markers at other purge rates. They are *not*
digitised here: the ranked-validation policy of this gallery takes a stated
numerical result over a digitised figure, and the stated result plus two
internal identities (below) close the model without any figure work. Digitising
Figure 4 would add breadth and could be done later without touching anything on
this page."""))

cells.append(code('''vals = load_data("itoh-1987-stated-values.csv", page=PAGE)
meta = load_meta("itoh-1987-stated-values.csv", page=PAGE)
print(vals.to_string(index=False, columns=["quantity", "value", "unit", "kind"]))
row = lambda q: float(vals.loc[vals.quantity == q, "value"].iloc[0])
X_MEASURED = row("measured_conversion")
X_EQ_STATED = row("equilibrium_conversion")
UC0_STATED = row("u_C0")
VA0_STATED = row("v_A0")
print(f"\\n{cite_data(meta)}")'''))

cells.append(md(r"""### Reconstructing $K_p(473\,\mathrm{K})$ from the stated equilibrium conversion

For the feed (1 mol C, $n_{Ar} = (1-y_C^0)/y_C^0$ mol Ar) at conversion $x$
and total pressure $P_0$, Eq. 7's equilibrium condition
$K_p p_C/p_H^3 = p_B$ rearranges to

$$
K_p(x) \;=\; \frac{27\,x^4\,P_0^3}{(1-x)\,\bigl(1 + n_{Ar} + 3x\bigr)^3},
$$

so the printed $x_{eq} = 0.187$ fixes $K_p$ at the one temperature the model
needs. The decision table then evaluates the printed expression and its
possible repairs against the paper's own statement, and a van 't Hoff estimate
from standard formation enthalpies and entropies (NIST WebBook values,
constant-$\Delta H$ approximation) checks the reconstruction's magnitude
independently.

The table also localises the corruption. The printed *prefactor* is of the
right thermodynamic family — $\exp(\Delta S^\circ/R) = 9.0\times10^{33}$ Pa³
against the printed $4.89\times10^{35}$, a factor 54 — whereas the printed
*argument* is wrong in sign and by 27 orders of magnitude. Repairing only the
argument, i.e. asking what slope makes the printed prefactor reproduce the
paper's own $x_{eq}$, gives $-26{,}480$ K, within 7 % of the van 't Hoff slope
$-\Delta H^\circ/R = -24{,}810$ K. That is the most economical diagnosis: one
number was mangled in typesetting, and repairing only that number reproduces
both the paper's own equilibrium conversion and standard thermodynamics. It is
a diagnosis, not a source: the page uses the reconstruction from $x_{eq}$,
which needs no such reasoning."""))

cells.append(code('''def kp_of_xeq(x):
    """Eq. 7's equilibrium condition solved for K_p at feed dilution Y_C0."""
    tot = 1.0 + N_AR + 3.0 * x
    return 27.0 * x**4 * P0**3 / ((1.0 - x) * tot**3)

def xeq_of_kp(Kp):
    """Inverse: equilibrium conversion for a given K_p (brentq, deterministic)."""
    f = lambda x: kp_of_xeq(x) - Kp
    if f(1.0 - 1e-9) < 0:          # even x -> 1 cannot reach this Kp
        return 1.0
    return brentq(f, 1e-9, 1.0 - 1e-9, xtol=1e-15)

K_P = kp_of_xeq(X_EQ_STATED)

# independent magnitude check: van 't Hoff from standard formation data
# (NIST WebBook: dHf298 benzene +82.9, cyclohexane -123.4 kJ/mol;
#  S298 benzene 269.2, cyclohexane 298.35, H2 130.68 J/mol/K)
DH298 = (82.9 - (-123.4)) * 1e3            # J/mol
DS298 = 269.2 + 3 * 130.68 - 298.35        # J/mol/K
lnK = -(DH298 - 298.15 * DS298) / (8.314 * 298.15) \
      + DH298 / 8.314 * (1 / 298.15 - 1 / T)
K_P_VH = np.exp(lnK) * (1e5) ** 3          # Pa^3 (1 bar reference)

SLOPE_REPAIR = T * np.log(K_P / 4.89e35)       # argument that fixes the printed prefactor
SLOPE_VH = -DH298 / 8.314                      # van 't Hoff slope, -dH/R
PREFAC_VH = np.exp(DS298 / 8.314) * (1e5) ** 3 # exp(dS/R) in Pa^3

print(f"K_p(473) reconstructed from x_eq = {X_EQ_STATED}: {K_P:.4e} Pa^3")
print(f"K_p(473) van 't Hoff estimate (independent):  {K_P_VH:.2e} Pa^3 "
      f"(ratio {K_P/K_P_VH:.2f} - consistent for a constant-dH estimate)")
print(f"standard reaction enthalpy from the same data: dH298 = {DH298/1e3:+.1f} kJ/mol")
print("  -> endothermic, so K_p must RISE with T: the exponential argument must")
print("     be negative, and the printed +3190/T cannot be right on its own.\\n")

print("decision table - each K_P candidate vs the paper's own x_eq = 0.187:")
print(f"  {'candidate':<40}{'K_p(473) [Pa^3]':>17}{'x_eq':>8}")
for label, kp in [
    ("printed: 4.89e35 exp(+3190/T)", 4.89e35 * np.exp(+3190 / T)),
    ("printed value read as atm^3, not Pa^3", 4.89e35 * np.exp(+3190 / T) * P0**3),
    ("sign repair: 4.89e35 exp(-3190/T)", 4.89e35 * np.exp(-3190 / T)),
    ("digit repair: 4.89e3 exp(+3190/T)", 4.89e3 * np.exp(+3190 / T)),
    (f"argument repair: 4.89e35 exp({SLOPE_REPAIR:.0f}/T)",
     4.89e35 * np.exp(SLOPE_REPAIR / T)),
    ("reconstructed from x_eq (this page)", K_P),
]:
    print(f"  {label:<40}{kp:>17.3e}{xeq_of_kp(kp):>8.3f}")
print("  -> the printed expression, a units re-reading, a sign flip and a digit")
print("     drop all fail against the paper's own 18.7 %. Only repairing the")
print("     exponential ARGUMENT works, and the slope it needs,")
print(f"     {SLOPE_REPAIR:.0f} K, is within "
      f"{abs(SLOPE_REPAIR - SLOPE_VH) / abs(SLOPE_VH) * 100:.0f} % of the van 't Hoff")
print(f"     slope -dH/R = {SLOPE_VH:.0f} K, while the printed prefactor is within a")
print(f"     factor {4.89e35 / PREFAC_VH:.0f} of exp(dS/R) = {PREFAC_VH:.2e} Pa^3.")
print("     Diagnosis only: the page uses the reconstruction from x_eq, which")
print("     assumes nothing about where the printed line came from.")'''))

cells.append(md(r"""## PyMRM implementation

The two coupled channels are one steady convection–reaction system on
$L \in [0,1]$ with state $(u_C, u_H, v_H)$ per cell — spatial axis first,
fields last, shape `(n_z, 3)`. Benzene and both argon flows are algebraic.
Convection is `construct_convflux_upwind` with $v = 1$ (the balances are
already in flow-rate form, so the "velocity" along dimensionless length is
unity) and `construct_div` with `nu=0` (Cartesian: $L$ is not a radial
coordinate). The constant operators are assembled **once**; only the pointwise
source is re-linearised, with `NumJac(shape)` coupling the last axis only.

Three numerical points, each of which cost real debugging time and is worth
stating:

- **Scale the unknowns to O(1).** `NumJac` floors its finite-difference
  perturbation at an *absolute* 10⁻⁶; molar flows of order 10⁻⁷ mol/s would be
  perturbed by more than their own value and the Jacobian becomes noise. The
  solver therefore works in $w = (u_C, u_H, v_H)/u_C^0$.
- **Regularise Eq. 7 algebraically, not with branches.** Multiplying through
  by $p_H^3$ gives $r_C = -k\,(K_p p_C - p_B\,p_H^3)/(p_H^3 + K_B K_p p_C)$,
  which is smooth at the $p_H \to 0$ inlet (limit $-k/K_B$: the reaction is
  zero-order when far from equilibrium) instead of the printed form's 0·∞.
- **Damp the Newton steps.** Near complete conversion the solution rides an
  equilibrium plateau ($u_C$ small, $p_H$ small) where the rate switches from
  zero-order to equilibrium-pinned over a very narrow range; undamped Newton
  overshoots and diverges. A backtracking line search on the residual norm
  (halve the step until the residual decreases) fixes it, starting from a
  deterministic explicit-march initial guess. Every Newton solve on this page
  records its iteration count and the observed range is printed in the
  validation section, so no iteration count in this prose can go stale.

Boundary conditions use the **outward normal** ($a\,\partial w/\partial n +
b\,w = d$): at the inlet face the physical condition is $w = (1, 0, 0)$
(Dirichlet, $a{=}0, b{=}1$); at the outlet, pure-convection outflow
$\partial w/\partial n = 0$ ($a{=}1, b{=}0$) — leaving the outlet `None` would
make the system singular."""))

cells.append(code('''NEWTON_ITERS = []          # every Newton solve on this page appends its count


def source(u, uC0, vA0, Kp=None, alphaH=ALPHA_H, kfac=1.0):
    """Pointwise source for state [..., (uC, uH, vH)] in mol/s.

    Clipping to >= 0 guards NumJac probe points; the converged solution is
    positive. Eq. 7 is used in the p_H^3-multiplied form (smooth at p_H = 0).
    `kfac` scales the rate constant for the sensitivity sweep; because the rate
    enters only as k*V_r, it equally represents an error in V_r.
    """
    Kp = K_P if Kp is None else Kp
    uC = np.maximum(u[..., 0], 0.0)
    uH = np.maximum(u[..., 1], 0.0)
    vH = np.maximum(u[..., 2], 0.0)
    uA = N_AR * uC0                      # argon, constant (Eq. 6)
    uB = uC0 - uC                        # benzene closure (Eq. 4)
    su = uC + uB + uH + uA               # total reaction-side flow
    sv = vA0 + vH                        # total shell-side flow
    pC, pB, pH = uC / su * P0, uB / su * P0, uH / su * P0
    # Eq. 7 times pH^3 / pH^3 -> regular everywhere:
    rC = -kfac * K_RATE * (Kp * pC - pB * pH**3) / (pH**3 + K_B * Kp * pC + 1e-30)
    # Eq. 1 with P_Tr = P_Ts = P0: sqrt of mole fractions
    perm = alphaH * (np.sqrt(uH / su) - np.sqrt(vH / sv))
    return np.stack([rC * V_R, -3.0 * rC * V_R - perm, perm], axis=-1)


class MembraneReactor:
    """Co-current Pd membrane reactor, Eqs. 1-8, as a pymrm FV/Newton solve.

    Unknowns are the dimensionless flows w = (uC, uH, vH)/uC0 so that NumJac's
    absolute perturbation floor (1e-6) stays small relative to the state.
    """

    def __init__(self, uC0, vA0, n_z=400, Kp=None, alphaH=ALPHA_H, kfac=1.0):
        self.uC0, self.vA0 = uC0, vA0
        self.Kp = K_P if Kp is None else Kp
        self.alphaH, self.kfac = alphaH, kfac
        self.z_f = np.linspace(0.0, 1.0, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.shape = (n_z, 3)
        # outward normal, a.dw/dn + b.w = d :
        #   inlet  (L=0): w = (1, 0, 0)   -> a=0, b=1  (everything fed at L=0, Eq. 8)
        #   outlet (L=1): dw/dn = 0       -> a=1, b=0  (convective outflow; None = singular)
        self.bc = ({"a": 0.0, "b": 1.0, "d": np.array([1.0, 0.0, 0.0])},
                   {"a": 1.0, "b": 0.0, "d": 0.0})
        self._build_operators()
        self.u = self._march()           # deterministic initial guess

    def _build_operators(self):
        # v = 1: the balances are in flow-rate form, d(w)/dL = source.
        # nu = 0: L is a Cartesian (axial) coordinate.
        conv, conv_bc = construct_convflux_upwind(
            self.shape, self.z_f, self.z_c, self.bc, v=1.0, axis=0)
        div = construct_div(self.shape, self.z_f, nu=0, axis=0)
        self.jac_const = div @ conv
        self.g_const = (div @ conv_bc).toarray().ravel()[:, None]
        self.numjac = NumJac(self.shape)     # pointwise source: last axis only

    def scaled_source(self, w):
        return source(w * self.uC0, self.uC0, self.vA0, self.Kp,
                      self.alphaH, self.kfac) / self.uC0

    def _march(self):
        """Explicit upwind march of the same grid: cheap, deterministic guess."""
        n_z = self.shape[0]
        u = np.zeros(self.shape)
        state = np.array([1.0, 0.0, 0.0])
        for i in range(n_z):
            state = np.maximum(state + self.scaled_source(state) / n_z, 0.0)
            state[0] = min(state[0], 1.0)
            u[i] = state
        return u

    def residual(self, w):
        g_s, jac_s = self.numjac(self.scaled_source, w)
        g = self.g_const + self.jac_const @ w.reshape((-1, 1)) - g_s.reshape((-1, 1))
        return g, self.jac_const - jac_s

    def solve(self, tol=1e-11, maxit=100):
        """Damped Newton: backtracking line search on the residual inf-norm."""
        w = self.u.copy()
        g, J = self.residual(w)
        ng = np.abs(g).max()
        for it in range(maxit):
            if ng < tol:
                self.u = w
                NEWTON_ITERS.append(it)
                return True
            dw = spsolve(J.tocsc(), -np.asarray(g).ravel()).reshape(self.shape)
            lam = 1.0
            while lam > 1e-6:
                w_try = w + lam * dw
                g2, J2 = self.residual(w_try)
                ng2 = np.abs(g2).max()
                if ng2 < (1.0 - 0.5 * lam) * ng + 1e-14:
                    break
                lam *= 0.5
            w, g, J, ng = w_try, g2, J2, ng2
        self.u = w
        NEWTON_ITERS.append(maxit)
        return bool(ng < tol)

    @property
    def conversion(self):
        return 1.0 - self.u[-1, 0]

    @property
    def h_balance_error(self):
        """Corrected Eq. 5 as an identity: 3(1 - uC) = uH + vH at the outlet."""
        uC, uH, vH = self.u[-1]
        return (uH + vH - 3.0 * (1.0 - uC)) / (3.0 * (1.0 - uC))


def conversion_fv(uC0, vA0, n_z=400, **kw):
    m = MembraneReactor(uC0, vA0, n_z=n_z, **kw)
    ok = m.solve()
    assert ok, f"Newton did not converge at uC0={uC0}, vA0={vA0}"
    return m


def conversion_ivp(uC0, vA0, Kp=None, alphaH=ALPHA_H, rtol=1e-10):
    """Independent reference: stiff IVP integration of Eqs. 2-8 (the paper
    used Runge-Kutta-Gill; Radau plays that role here)."""
    rhs = lambda L, y: source(np.asarray(y), uC0, vA0, Kp, alphaH)
    sol = solve_ivp(rhs, [0.0, 1.0], [uC0, 0.0, 0.0], method="Radau",
                    rtol=rtol, atol=1e-14)
    assert sol.success
    return 1.0 - sol.y[0, -1] / uC0'''))

cells.append(md(r"""### Two identities before any reactor is run

The permeation constant and the equilibrium limit can both be checked without
solving anything. The first shows which radii Itoh put into Eq. 1 (and nothing
more — see the note under Parameters); the second is the kinetics-free ceiling
that confirms the printed $v_A^0$."""))

cells.append(code('''print("identity 1: alpha_H from Eq. 1's own inputs")
for label, ri, ro in [("ID = 17.0 mm (r_i=8.5, r_o=8.7)", 8.5e-3, 8.7e-3),
                      ("OD = 17.0 mm (r_i=8.3, r_o=8.5)", 8.3e-3, 8.5e-3)]:
    a = 2 * np.pi * L0 * D_PD * C_0 / np.log(ro / ri)
    print(f"  {label}: {a:.4e} mol/s  vs printed {ALPHA_H:.2e}  "
          f"({(a - ALPHA_H) / ALPHA_H * 100:+.2f} %)")
ALPHA_H_RECOMPUTED = 2 * np.pi * L0 * D_PD * C_0 / np.log(RO / RI)
print("  -> the printed alpha_H was computed with r_i = 8.5 mm (0.03 %); this")
print("     says what Itoh substituted into Eq. 1, not what the tube physically")
print("     was, and it propagates to nothing: the model uses the printed")
print("     alpha_H, and r_i enters only through V_r, which appears only as the")
print("     product k*V_r (see the k sensitivity in check 1).\\n")

print("identity 2: the fast-permeation ceiling (Kp fixed, everything else infinite)")

def ceiling(uC0, vA0, Kp=None):
    """Max conversion with alpha_H -> inf and instantaneous kinetics: equal H2
    mole fraction y on both sides, reaction at equilibrium. Kinetics-free.
    Equilibrium Kp*pC = pB*pH^3 reduces to Kp(1-x) = x (y P0)^3."""
    Kp = K_P if Kp is None else Kp
    F0 = uC0 * (1.0 + N_AR)                     # non-H2 reaction-side flow
    def f(x):
        q = 3.0 * x * uC0 / (F0 + vA0)          # y/(1-y) for the common H2 fraction
        y = q / (1.0 + q)
        return Kp * (1.0 - x) - x * (y * P0) ** 3
    return brentq(f, 1e-12, 1.0 - 1e-12, xtol=1e-15)

for va, label in [(11.8e-6, "v_A0 = 11.8e-6 (the 10^-6 reading)"),
                  (11.8e-5, "v_A0 = 11.8e-5 (as printed)")]:
    print(f"  {label}: ceiling X = {ceiling(UC0_STATED, va):.4f}")
print(f"  measured X = {X_MEASURED}")
print(f"  same ceiling at 11.8e-6 with the van 't Hoff K_p instead: "
      f"{ceiling(UC0_STATED, 11.8e-6, K_P_VH):.4f}")
KP_NEEDED = np.exp(brentq(lambda lk: ceiling(UC0_STATED, 11.8e-6, np.exp(lk))
                          - X_MEASURED, np.log(1e10), np.log(1e18), xtol=1e-12))
print(f"  K_p that WOULD allow 0.997 at 11.8e-6: {KP_NEEDED:.2e} Pa^3 "
      f"({KP_NEEDED / K_P:.0f}x the reconstruction),")
print(f"  which implies a membrane-free equilibrium conversion of "
      f"{xeq_of_kp(KP_NEEDED):.3f} against the printed {X_EQ_STATED}.")
print("  -> at 11.8e-6 mol/s purge not even infinitely fast permeation and")
print("     kinetics reach the measured conversion: the co-current shell")
print("     accumulates enough H2 to stall the reaction (a pinch). This")
print("     CONFIRMS the printed 10^-5 exponent rather than correcting it, and")
print("     it needs only that standard thermodynamics is not wrong by two and")
print("     a half orders of magnitude. Figure 4's abscissa is then x 10^5.")'''))

cells.append(md(r"""## Results

The headline run first, then the operating map: conversion against purge rate
for the three feed rates whose curves Figure 4 draws (legend labels 0.29, 0.80,
1.64 × 10⁻⁶ mol/s — printed numbers). The sweep is this page's own computation
of the paper's calculated curves: a reproduction, not a validation, since the
figure's markers are not digitised."""))

cells.append(code('''m_stated = conversion_fv(UC0_STATED, VA0_STATED, n_z=800)
X_model = m_stated.conversion
DEV = (X_model - X_MEASURED) / X_MEASURED
DEV_1MX = ((1 - X_model) - (1 - X_MEASURED)) / (1 - X_MEASURED)
print(f"stated run: u_C0 = {UC0_STATED:.2e} mol/s, v_A0 = {VA0_STATED:.2e} mol/s")
print(f"  measured conversion  {X_MEASURED:.3f}")
print(f"  model conversion     {X_model:.4f}")
print(f"  deviation (model - measured)/measured = {DEV * 100:+.2f} %")
print(f"  ... but on the quantity that actually varies, the unconverted")
print(f"  fraction 1 - X: measured {1 - X_MEASURED:.4f}, model {1 - X_model:.5f},")
print(f"  deviation {DEV_1MX * 100:+.0f} %. Both numbers describe the same solve;")
print(f"  see validation check 1 for which one carries information.")
print(f"  equilibrium without the membrane: {X_EQ_STATED} -> the membrane")
print(f"  multiplies the attainable conversion by {X_model / X_EQ_STATED:.1f}")

L = m_stated.z_c
w = m_stated.u
su = (w[:, 0] + (1 - w[:, 0]) + w[:, 1] + N_AR)      # scaled total, tube side
sv = VA0_STATED / UC0_STATED + w[:, 2]
fig, ax = plt.subplots(1, 2, figsize=(11.8, 4.1))
ax[0].plot(L, w[:, 0], lw=2, label="$u_C/u_C^0$ cyclohexane")
ax[0].plot(L, 1 - w[:, 0], lw=2, label="$u_B/u_C^0$ benzene")
ax[0].plot(L, w[:, 1], lw=2, label="$u_H/u_C^0$ H$_2$, tube")
ax[0].plot(L, w[:, 2], lw=2, label="$v_H/u_C^0$ H$_2$, shell")
ax[0].set(xlabel="dimensionless length $L$", ylabel="flow / $u_C^0$",
          title="the stated run: H$_2$ is exported as fast as it forms")
ax[0].legend(fontsize=8)
ax[1].semilogy(L, w[:, 1] / su, lw=2, label="$y_{H}$ tube side")
ax[1].semilogy(L, w[:, 2] / sv, lw=2, label="$y'_{H}$ shell side")
ax[1].set(xlabel="dimensionless length $L$", ylabel="H$_2$ mole fraction",
          title="the membrane keeps $p_H$ low everywhere")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()'''))

cells.append(code('''UC0_CURVES = np.array([0.29e-6, 0.80e-6, 1.64e-6])   # printed legend labels
VA0_SWEEP = np.array([0.25, 0.5, 1, 2, 3, 4, 6, 8, 10, 11.8, 13.5, 15]) * 1e-5

fig, ax = plt.subplots(figsize=(7.2, 4.6))
for uC0, col in zip(UC0_CURVES, ("tab:blue", "tab:green", "tab:red")):
    X = [conversion_fv(uC0, va).conversion for va in VA0_SWEEP]
    ax.plot(VA0_SWEEP * 1e5, X, "-", lw=2, color=col,
            label=f"$u_C^0$ = {uC0 * 1e6:.2f} $\\\\times 10^{{-6}}$ mol/s")
ax.plot(VA0_STATED * 1e5, X_MEASURED, "k*", ms=16, zorder=5,
        label="measured, 99.7 % (stated in text)")
ax.axhline(X_EQ_STATED, color="0.4", ls="--", lw=1.4)
ax.text(9.0, X_EQ_STATED + 0.015, "equilibrium conversion 0.187", fontsize=8, color="0.3")
ax.set(xlabel="purge flow rate $v_A^0 \\\\times 10^5$ [mol/s]",
       ylabel="conversion $X$", ylim=(0, 1.02),
       title="the paper's Figure 4, recomputed (curves) — one stated measurement (star)")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); plt.show()
print("Every curve rises from the equilibrium conversion at zero purge toward")
print("complete conversion: the membrane turns a thermodynamic ceiling into a")
print("purge-rate knob. The curves reproduce the calculated curves of the")
print("paper's Figure 4 (a reproduction; the figure's markers are not digitised).")'''))

cells.append(md(r"""## Validation

Seven checks. The stated measured conversion is the one comparison with a
measurement, and check 1 spends most of its length establishing how much that
comparison can actually resolve — which is less than the headline deviation
suggests. The rest are identities the paper pays for, limits, and
method-independence checks: reproduction, not validation."""))

cells.append(code('''print("1. the measured conversion - and exactly how much it can resolve")
print(f"   measured {X_MEASURED:.3f}, model {X_model:.4f}, "
      f"deviation {DEV * 100:+.2f} %")
print("   Nothing was fitted: kinetics and alpha_H are printed constants, Kp")
print("   comes from the printed equilibrium conversion, geometry from Fig. 1.")
ceil_stated = ceiling(UC0_STATED, VA0_STATED)
print(f"   |model - fast-permeation ceiling| = {abs(X_model - ceil_stated):.1e}")
print("   At this operating point the model IS its own algebraic ceiling - a")
print("   closed-form function of (Kp, uC0, vA0, yC0, P0) alone. 'Nothing fitted'")
print("   is true, but almost nothing could have made it fail. Perturbing the")
print("   model shows what the comparison can and cannot see:\\n")
print(f"   {'variant':<32}{'X_model':>10}{'dev vs 0.997':>14}")
SENS = {}
for label, kw in [
    ("baseline", {}),
    ("alpha_H x 0.3", dict(alphaH=0.3 * ALPHA_H)),
    ("alpha_H x 10", dict(alphaH=10.0 * ALPHA_H)),
    ("k (== k*V_r) x 0.1", dict(kfac=0.1)),
    ("k (== k*V_r) x 10", dict(kfac=10.0)),
    ("Kp / 1.58 (the van 't Hoff value)", dict(Kp=K_P_VH)),
    ("Kp x 1.58", dict(Kp=1.58 * K_P)),
    ("Kp / 3", dict(Kp=K_P / 3)),
    ("Kp x 10", dict(Kp=10.0 * K_P)),
    ("Kp / 10", dict(Kp=K_P / 10)),
]:
    Xs = conversion_fv(UC0_STATED, VA0_STATED, n_z=800, **kw).conversion
    SENS[label] = Xs
    print(f"   {label:<32}{Xs:>10.5f}"
          f"{(Xs - X_MEASURED) / X_MEASURED * 100:>+13.2f} %")

SENS_ALPHAH = max(abs(SENS["alpha_H x 0.3"] - SENS["baseline"]),
                  abs(SENS["alpha_H x 10"] - SENS["baseline"]))
SENS_K = max(abs(SENS["k (== k*V_r) x 0.1"] - SENS["baseline"]),
             abs(SENS["k (== k*V_r) x 10"] - SENS["baseline"]))
print(f"\\n   alpha_H over a factor 33: X moves by {SENS_ALPHAH:.1e}")
print(f"   k (and hence V_r, and hence r_i) over a factor 100: "
      f"X moves by {SENS_K:.1e}")
print("   So this comparison tests NEITHER the Langmuir-Hinshelwood kinetics NOR")
print("   the Sieverts permeance. It bounds Kp to roughly an order of magnitude,")
print("   and it does not distinguish the reconstruction from the van 't Hoff")
print("   value - which in fact agrees slightly BETTER with the measurement.")
print("   What it does establish is that the reactor reached its co-current")
print("   fast-permeation asymptote and that the model puts that asymptote in")
print("   the right place.\\n")
print("   Two further honesty points on the same number:")
print(f"   (a) on the quantity that actually varies, the unconverted fraction:")
print(f"       1 - X measured {1 - X_MEASURED:.4f}, model {1 - X_model:.5f}, "
      f"deviation {DEV_1MX * 100:+.0f} %.")
print(f"       The model under-predicts unconverted cyclohexane by "
      f"{abs(DEV_1MX) * 100:.0f} %; that, not")
print(f"       {DEV * 100:+.2f} %, is the honest measure of the disagreement.")
DEV_LO = (X_model - 0.9975) / 0.9975
DEV_HI = (X_model - 0.9965) / 0.9965
print(f"   (b) '0.997' is three digits with no uncertainty, so X is only known to")
print(f"       lie in [0.9965, 0.9975] and the deviation only in "
      f"[{DEV_LO * 100:+.2f} %, {DEV_HI * 100:+.2f} %].\\n")

print("2. alpha_H identity: Eq. 1 recomputed from its own printed inputs")
print(f"   printed 4.47e-5, recomputed {ALPHA_H_RECOMPUTED:.4e}, "
      f"deviation {(ALPHA_H_RECOMPUTED - ALPHA_H) / ALPHA_H * 100:+.2f} %")
print("   (shows Itoh used r_i = 8.5 mm in Eq. 1; consequential for nothing on")
print("    this page - see the note under Parameters)\\n")

print("3. equilibrium limit: alpha_H = 0 recovers the stated 18.7 %")
m_eq = conversion_fv(2.9e-8, 1e-9, alphaH=0.0)   # slow feed, no permeation
print(f"   X(alpha_H=0, long residence) = {m_eq.conversion:.4f} vs stated {X_EQ_STATED}")
print(f"   (tests the reactor implementation against the closed-form Kp inversion;")
print(f"    Kp itself came from that statement, so this is a consistency check,")
print(f"    not independent evidence)\\n")

print("4. two independent methods and grid convergence")
X_ivp = conversion_ivp(UC0_STATED, VA0_STATED)
print(f"   stiff-IVP reference (Radau, the paper's own RKG role): X = {X_ivp:.6f}")
print(f"   FV at the stated run, every grid: X = "
      f"{conversion_fv(UC0_STATED, VA0_STATED, n_z=100).conversion:.6f}")
print(f"   |FV - IVP| = {abs(m_stated.conversion - X_ivp):.1e}, "
      f"|FV - ceiling| = {abs(m_stated.conversion - ceil_stated):.1e}")
print("   At the stated run there is no grid error to show: the outlet sits on")
print("   the fully-equilibrated fixed point (rate = 0, permeation = 0), which")
print("   IS the ceiling state - so FV at any n_z, the IVP, and the algebraic")
print("   ceiling all agree to ~1e-14. The ceiling claim made concrete.")
print("   Grid convergence is therefore shown mid-map, where the solution is")
print("   NOT equilibrated and the discretisation genuinely matters:")
X_ivp_mid = conversion_ivp(0.80e-6, 3e-5)
print(f"   (u_C0 = 0.8e-6, v_A0 = 3e-5 mol/s; IVP reference X = {X_ivp_mid:.6f})")
print(f"   {'n_z':>6}{'X_fv':>12}{'X_fv - X_ivp':>14}{'ratio':>8}")
errs = []
for n_z in (100, 200, 400, 800, 1600):
    Xn = conversion_fv(0.80e-6, 3e-5, n_z=n_z).conversion
    errs.append(abs(Xn - X_ivp_mid))
    ratio = f"{errs[-2] / errs[-1]:.2f}" if len(errs) > 1 else ""
    print(f"   {n_z:>6}{Xn:>12.6f}{Xn - X_ivp_mid:>+14.2e}{ratio:>8}")
X_FV_VS_IVP = errs[3]
print("   refinement ratio ~2: first-order upwind convergence toward the")
print("   independent reference, as it must be.\\n")

print("5. conservation (the corrected Eq. 5 as an identity) - and where it bites")
HBAL = abs(m_stated.h_balance_error)
m_mid_h = conversion_fv(0.80e-6, 3e-5)
HBAL_MID = abs(m_mid_h.h_balance_error)
print(f"   co-current, stated run: |uH + vH - 3(uC0-uC)|/3(uC0-uC) = {HBAL:.1e}")
print(f"   co-current, mid-map (X = {m_mid_h.conversion:.3f}, far from "
      f"equilibrium):  {HBAL_MID:.1e}")
print("   Both sit at rounding level - and that is a WEAK result, not a strong")
print("   one: co-current both streams are upwinded in the same direction, the")
print("   discrete sums telescope, and the identity therefore holds at any n_z.")
print("   It checks the implementation, not the discretisation. The identity is")
print("   NOT automatic counter-current, where the streams are upwinded in")
print("   opposite directions; that is the case where it proves something, and")
print("   it is measured (with its grid ladder) in 'What pymrm adds' below.")
print("   Carbon closure is algebraic (uB = uC0 - uC) and argon constant by")
print("   construction; the printed Eq. 5 (sign-flipped) would give vH < 0.\\n")

print("6. the model bounded by its own ceiling everywhere on the sweep")
worst = 0.0
for uC0 in UC0_CURVES:
    for va in VA0_SWEEP[::2]:
        worst = max(worst, conversion_fv(uC0, va).conversion - ceiling(uC0, va))
print(f"   max(X_model - X_ceiling) over the map = {worst:.2e}  (<= 0 up to")
print(f"   discretisation: finite kinetics and permeance never beat the bound)\\n")

print("7. Newton behaviour, instrumented rather than remembered")
print(f"   {len(NEWTON_ITERS)} damped-Newton solves up to this point on the page;")
print(f"   iteration counts {min(NEWTON_ITERS)} to {max(NEWTON_ITERS)} "
      f"(the whole-page range is in the CI metrics at the end)")
print("   (all from the same deterministic explicit-march guess; no warm start,")
print("    no continuation - every solve on this page is built from scratch)")'''))

cells.append(md(r"""## What pymrm adds

**The paper's own computation is an initial-value problem** — everything enters
at $L=0$, so Runge–Kutta marches it. Recasting it as a pymrm finite-volume
system with Newton solves reproduces it (check 4), which by itself adds only
robustness. The genuine addition is that the FV/Newton form does not care which
end a stream enters: **counter-current purge is the same assembly** with the
shell field given velocity $-1$ and its Dirichlet condition moved to $L=1$ —
per-field values in the same `bc` tuple — while a marching integrator would
need shooting. The comparison below is a *prediction*; Itoh's experiment was
co-current only, and no data exists to test it.

Unlike the co-current solve, the counter-current one is **not** grid-insensitive
— the co-current stated run sits on an algebraic fixed point and returns the
same number at every $n_z$, whereas the counter-current solution converges only
first order in $h$. It is therefore run on a converged grid, with its ladder and
its hydrogen closure printed below, and quoted to two decimals.

Two things the comparison shows, both mechanistic rather than obvious:

- **At moderate purge, counter-current wins, and the margin is large.** Fresh
  argon meets the reactor *outlet*, where the last of the cyclohexane needs the
  deepest hydrogen removal; the co-current reactor instead delivers its most
  hydrogen-loaded purge exactly there — the pinch that produced the 0.72
  ceiling above.
- **At very low purge the order reverses.** With little argon, the
  counter-current shell carries its full hydrogen load past the reactor
  *inlet*, back-permeating hydrogen into the fresh feed where the tube-side
  $p_H$ is still low. Co-current avoids that because both streams start clean
  together.

The archival contribution of this page is the misprint resolution: **two**
printed statements — Eq. 5's sign and the $K_P$ expression — are wrong as
printed, and both corrections are derived from the paper's own statements. The
third value a reader might suspect, $v_A^0 = 11.8\times10^{-5}$ mol/s, is
printed correctly and is confirmed here rather than corrected. All of it is
recorded so that nobody re-derives it from scratch."""))

cells.append(code('''class CounterCurrentReactor(MembraneReactor):
    """Same physics; purge enters at L = 1. Only the operators change:
    per-field velocity (+1, +1, -1) and per-field boundary values."""

    def _build_operators(self):
        # outward normal, per field (uC, uH, vH):
        #   L=0: uC, uH Dirichlet inlet (b=1); vH outflow d(vH)/dn = 0 (a=1)
        #   L=1: uC, uH outflow (a=1);        vH Dirichlet inlet vH = 0 (b=1)
        self.bc = ({"a": np.array([0.0, 0.0, 1.0]),
                    "b": np.array([1.0, 1.0, 0.0]),
                    "d": np.array([1.0, 0.0, 0.0])},
                   {"a": np.array([1.0, 1.0, 0.0]),
                    "b": np.array([0.0, 0.0, 1.0]),
                    "d": np.array([0.0, 0.0, 0.0])})
        conv, conv_bc = construct_convflux_upwind(
            self.shape, self.z_f, self.z_c, self.bc,
            v=np.array([1.0, 1.0, -1.0]), axis=0)
        div = construct_div(self.shape, self.z_f, nu=0, axis=0)
        self.jac_const = div @ conv
        self.g_const = (div @ conv_bc).toarray().ravel()[:, None]
        self.numjac = NumJac(self.shape)

    @property
    def h_balance_error(self):
        """Counter-current closure: uH leaves at L=1, vH leaves at L=0."""
        uC, uH = self.u[-1, 0], self.u[-1, 1]
        vH = self.u[0, 2]
        return (uH + vH - 3.0 * (1.0 - uC)) / (3.0 * (1.0 - uC))


# The counter-current solve is NOT grid-insensitive the way the co-current one
# is (which sits on an algebraic fixed point), so its grid has to be earned.
N_CC = 3200
CO_VS_CEIL = max(abs(conversion_fv(UC0_STATED, va, n_z=n).conversion
                     - ceiling(UC0_STATED, va))
                 for va in (0.5e-6, 8e-6) for n in (400, 3200))
print(f"co-current over this purge range sits on the algebraic ceiling to "
      f"{CO_VS_CEIL:.1e},")
print("so its numbers are the same at every n_z. The counter-current ones are")
print("not, and their grid has to be earned:")
print(f"   {'v_A0':>10}{'n_z':>7}{'X_cc':>10}{'H closure':>12}")
CC_LADDER = {}
for va in (0.5e-6, 8e-6):
    for n in (400, 800, 1600, 3200, 6400):
        m = CounterCurrentReactor(UC0_STATED, va, n_z=n)
        assert m.solve(), f"counter-current Newton failed at vA0={va}, n_z={n}"
        CC_LADDER[(va, n)] = (m.conversion, m.h_balance_error)
        print(f"   {va:>10.1e}{n:>7d}{m.conversion:>10.4f}"
              f"{m.h_balance_error:>+12.1e}")
print("   First order in h, in the solution and in the closure alike. The class")
print(f"   default n_z = 400 is NOT converged here: X = "
      f"{CC_LADDER[(0.5e-6, 400)][0]:.4f} against "
      f"{CC_LADDER[(0.5e-6, 6400)][0]:.4f} at n_z = 6400.")
print(f"   The comparison below is therefore run at n_z = {N_CC}, and quoted to")
print("   two decimals only.\\n")

VA0_CC = np.array([0.5, 1, 2, 3, 4, 5, 6, 7, 8]) * 1e-6
X_co, X_cc, HB_cc = [], [], []
for va in VA0_CC:
    X_co.append(conversion_fv(UC0_STATED, va, n_z=N_CC).conversion)
    mcc = CounterCurrentReactor(UC0_STATED, va, n_z=N_CC)
    ok = mcc.solve()
    assert ok, f"counter-current Newton failed at vA0={va}"
    X_cc.append(mcc.conversion)
    HB_cc.append(abs(mcc.h_balance_error))
HBAL_CC = max(HB_cc)

fig, ax = plt.subplots(figsize=(7.0, 4.4))
ax.plot(VA0_CC * 1e6, X_co, "o-", lw=2, color="tab:blue", label="co-current (the paper)")
ax.plot(VA0_CC * 1e6, X_cc, "s-", lw=2, color="tab:orange",
        label=f"counter-current (prediction, no data), $n_z$ = {N_CC}")
ax.axhline(X_EQ_STATED, color="0.4", ls="--", lw=1.2)
ax.text(6.2, 0.20, "equilibrium 0.187", fontsize=8, color="0.3")
ax.set(xlabel="purge flow rate $v_A^0 \\\\times 10^6$ [mol/s]",
       ylabel="conversion $X$",
       title=f"flow arrangement at $u_C^0$ = {UC0_STATED*1e7:.1f}e-7 mol/s")
ax.legend(fontsize=8, loc="center right")
fig.tight_layout(); plt.show()

print(f"at v_A0 = 8e-6 mol/s: co-current X = {X_co[-1]:.2f}, "
      f"counter-current X = {X_cc[-1]:.2f}")
print(f"at v_A0 = 0.5e-6 mol/s the order reverses: "
      f"co {X_co[0]:.2f} vs counter {X_cc[0]:.2f}")
print(f"worst counter-current hydrogen closure over the sweep at n_z = {N_CC}: "
      f"{HBAL_CC:.1e}")
print(" - this is the conservation check that is not automatic (check 5): the")
print("   two streams are upwinded in opposite directions, so the closure is a")
print("   genuine O(h) discretisation error and it is reported as such, not as")
print("   'machine precision'.")
print("(the counter-current Newton solve is carried over the purge range where")
print(" it converges cleanly, 0.5-8e-6 mol/s; at higher purge the damped solve")
print(" stalls, most likely because the counter-current shell inlet has v_H = 0")
print(" exactly, where the sqrt(v_H) term of Eq. 1 has an infinite derivative;")
print(" no result from that region is reported)")
X_CC_8 = X_cc[-1]'''))

cells.append(md(r"""**What this page does not establish.** The comparison with
measurement is a *single* stated conversion, at an operating point where the
model sits exactly on its thermodynamic ceiling. Check 1 measures what that
costs: varying $\alpha_H$ over a factor 33 and $k V_r$ over a factor 100 does
not move the modelled conversion at all, and $K_p$ can be divided by 3 or
multiplied by 10 with the deviation staying inside ±0.3 %. **So this page is not
experimental evidence for Itoh's Langmuir–Hinshelwood kinetics, nor for the
Sieverts permeance, and it does not select the reconstructed $K_p$ over the
van 't Hoff value** — the latter agrees marginally better. What the measurement
does establish is that the reactor reached its co-current fast-permeation
asymptote and that the model puts that asymptote in the right place.

On the quantity that actually varies, the unconverted fraction $1-X$, the model
is 45 % low. That is the honest measure of the residual disagreement, and it is
consistent with $K_p$ being somewhat under-reconstructed — but the measurement
cannot pin it down, because $1-X$ is a 0.3 % effect on a conversion quoted to
three digits.

The low-purge region, where the kinetics and the permeation resistance
genuinely shape the curves, is checked only as a reproduction of the paper's
own calculated curves, not against the experimental markers of Figure 4, which
were not digitised. The $K_p$ used is a reconstruction from the paper's stated
equilibrium conversion; it is consistent with independent thermodynamic data at
the factor-1.6 level, but the printed expression it replaces remains unexplained
— most plausibly a garbled transcription from the kinetics source (Itoh et al.
1985, in Japanese), which was not available to check. The counter-current
comparison is a pure prediction, on a converged grid but with no data behind
it."""))

cells.append(code('''# CI regression metrics. Deliberately includes the *resolving power* of the
# headline comparison, so that a future reader (or a future page citing this
# one) cannot take +0.13 % as evidence about the kinetics or the permeance.
report_agreement("H1.4", {
    "conversion_measured": X_MEASURED,
    "conversion_model": X_model,
    "conversion_dev_frac": DEV,
    "unconverted_measured": 1 - X_MEASURED,
    "unconverted_model": 1 - X_model,
    "unconverted_dev_frac": DEV_1MX,
    "sens_alphaH_x0.3_to_x10_maxabs_dX": SENS_ALPHAH,
    "sens_kVr_x0.1_to_x10_maxabs_dX": SENS_K,
    "sens_Kp_vanthoff_dev_frac": (SENS["Kp / 1.58 (the van 't Hoff value)"]
                                  - X_MEASURED) / X_MEASURED,
    "sens_Kp_div3_dev_frac": (SENS["Kp / 3"] - X_MEASURED) / X_MEASURED,
    "sens_Kp_x10_dev_frac": (SENS["Kp x 10"] - X_MEASURED) / X_MEASURED,
    "alphaH_recomputed_dev_frac": (ALPHA_H_RECOMPUTED - ALPHA_H) / ALPHA_H,
    "Kp_473_reconstructed_Pa3": K_P,
    "Kp_473_vanthoff_ratio": K_P / K_P_VH,
    "Kp_needed_for_0.997_at_11.8e-6_ratio": KP_NEEDED / K_P,
    "equilibrium_limit_X": m_eq.conversion,
    "ceiling_at_stated_run": ceil_stated,
    "ceiling_at_1e-6_reading": ceiling(UC0_STATED, 11.8e-6),
    "stated_run_fv_vs_ceiling_absdiff": abs(m_stated.conversion - ceil_stated),
    "fv_vs_ivp_absdiff_midmap_nz800": X_FV_VS_IVP,
    "h_balance_rel_error_cocurrent": HBAL,
    "h_balance_rel_error_countercurrent_nz3200": HBAL_CC,
    "cc_X_0p5e-6_nz3200": CC_LADDER[(0.5e-6, 3200)][0],
    "cc_X_8e-6_nz3200": CC_LADDER[(8e-6, 3200)][0],
    "cc_grid_shift_0p5e-6_nz400_to_nz6400": (CC_LADDER[(0.5e-6, 6400)][0]
                                             - CC_LADDER[(0.5e-6, 400)][0]),
    "newton_iters_min": min(NEWTON_ITERS),
    "newton_iters_max": max(NEWTON_ITERS),
})'''))

cells.append(md(r"""## Reuse

**To adapt to your membrane reactor:** replace the rate expression in
`source` (and its equilibrium closure `kp_of_xeq`), the permeation law (here
Sieverts, $\propto \sqrt{p}$ — a porous membrane would be linear in $p$ and
non-selective, adding terms to *every* balance), and the constants block.
The `(n_z, 3)` layout extends to more species by widening the last axis; the
operators never change.

**The traps, in the order they will bite:**

- `NumJac` floors its perturbation at an absolute 10⁻⁶ — **scale unknowns to
  O(1)** or the Jacobian of any small-magnitude problem is noise. This
  produced silent divergence, not an error message.
- Rate expressions with equilibrium groups like $K_p p_C/p_H^3$ must be
  **regularised algebraically** (multiply through) before a Newton solver sees
  them; branch-based guards leave kinks that stall damped Newton.
- A pure-convection outlet left as `None` makes the system **singular**; give
  the outlet $\partial w/\partial n = 0$ on the outward normal.
- Near-complete conversion puts the solution on an equilibrium plateau where
  undamped Newton diverges — **backtracking damping** (six lines) fixes it.
- Boundary conditions use the **outward normal**: the same `{a, b, d}` dict
  means different physics at the two ends, and per-field arrays in one dict
  are how two counter-flowing streams share one operator.

**Related pages.** [`H1.7`](../H1.7-solution-diffusion/) (membrane transport,
solution–diffusion), [`F2.3`](../F2.3-slurry-bubble-column-ft/) (plug-flow
reactor with printed-constant corrections established the same way),
[`D2.2`](../D2.2-van-welsenaere-froment-runaway/) (the same
convection–source FV/Newton assembly), [`C2.1`](../C2.1-xu-froment-smr/)
(equilibrium-limited kinetics against experiment).

**Cite the source, not this page:** Itoh, N., *A membrane reactor using
palladium*, AIChE Journal **33**(9) 1576–1578 (1987),
[doi:10.1002/aic.690330921](https://doi.org/10.1002/aic.690330921)."""))

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
