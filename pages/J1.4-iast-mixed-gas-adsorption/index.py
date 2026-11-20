# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ---
# title: "IAST 1965: mixture adsorption from pure-component isotherms alone, and a paper that prints no numbers"
# description: "Myers and Prausnitz derive the ideal adsorbed solution theory, the Raoult's-law analogue in which P y_i = P_i°(π) x_i at matched spreading pressure, and validate it on four gas pairs. Every one of those comparisons is a figure of external data, and even the pure-component inputs are never tabulated, so this page validates what the paper actually prints: the derivation chain (verified symbolically, which surfaces two printed equation defects), the graphical lever-rule construction, and the thermodynamic structure. The solver is proved against two closed forms it never touches: extended Langmuir, which the page derives symbolically as the exact IAST solution for equal capacities, and the Henry-law limit the paper itself proves rigorous. The inconsistency the paper reports for the Langmuir mixture rule is then made quantitative: integrating the Gibbs adsorption isotherm around a closed loop in (P, y) leaves a 12 % spreading-pressure mismatch for unequal-capacity extended Langmuir, computed twice, by line integral and by a symbolically derived curl, where IAST closes the same loop to round-off (below 1e-12)."
# categories: [sec:J, struct:S1, struct:S10, tier:T0, data:tier6, phase:gas-solid]
# date: 2026-08-08
# ---
#
# # IAST 1965: mixture adsorption from pure-component isotherms alone, and a paper that prints no numbers
#
# **Catalog ID:** `J1.4` · **Structures:** `S1` (pointwise algebra), `S10`
# (constrained algebraic system) · **Tier:** T0
#
# The ideal adsorbed solution theory is the `A2.1` of adsorption: the default
# mixture model everyone reaches for, sixty years on. Its content is one
# sentence, the paper's own abstract sentence: *the partial pressure of an
# adsorbed component is given by the product of its mole fraction in the adsorbed phase
# and the pressure which it would exert as a pure adsorbed component at the same
# temperature and **spreading pressure** as those of the mixture*,
#
# $$P y_i \;=\; P_i^{\circ}(\pi)\,x_i \qquad \text{(eq. 22)},$$
#
# Raoult's law with the pure-component reference taken at matched spreading
# pressure instead of matched temperature alone. Everything else follows from
# classical surface thermodynamics, and mixture equilibria follow from
# pure-component isotherms with **no mixture data and no mixture parameters**.
#
# **What this page can and cannot check, decided before any code.** The paper
# validates IAST against four systems, methane–ethane and ethylene–carbon
# dioxide on activated carbon, carbon monoxide–oxygen and propane–propylene on
# silica gel. **Every one of those comparisons is a figure** (Figs. 3–10), the
# experimental points in them belong to four *other* papers, and even the
# pure-component isotherms that feed the theory are never tabulated, Fig. 2's
# own legend box says its curves are *calculated from* the isotherms of Szepesy
# and Illés. There is not one table in the paper. So the experimental case is
# out of scope here (§4 gives the search and the reasoning), and the page
# validates what the paper actually prints:
#
# 1. **The derivation chain, symbolically.** Eqs. (9)–(35) and (45)–(46) are
#    re-derived with `sympy`; eleven identities close to exactly zero. Doing this
#    surfaces **two printed equation defects**: eq. (11) is missing the `+`
#    between its two $RT\ln$ terms, as printed it does not yield eq. (13), the
#    equilibrium relation the whole theory rests on, and the sentence
#    introducing eq. (21) cites eq. (21) itself where the derivation needs
#    eq. (20). Both are reported, neither repaired silently.
# 2. **The solver, against closed forms it never touches.** The page derives
#    symbolically that for two Langmuir isotherms of **equal capacity** the IAST
#    solution *is* the extended Langmuir rule, and the numerical solver
#    reproduces that closed form to its own tolerance (below 1e-12); in the
#    Henry-law limit, the one case the paper itself proves IAST rigorous
#    (eqs. 30–35), the solver converges to the closed-form selectivity at first
#    order in $P$, as it must.
# 3. **The thermodynamic consistency the paper's rivals lack, made
#    quantitative.** The paper reports (p. 125) that the Langmuir mixture model
#    is "not thermodynamically consistent". This page turns that sentence into a
#    number: integrating the Gibbs adsorption isotherm (eq. 15) around a closed
#    loop in $(P, y_1)$ leaves a spreading-pressure mismatch of **12 % of the
#    loop's own scale** for extended Langmuir with unequal capacities, computed
#    twice, by line integral and by a symbolically derived curl, two routes that
#    agree to round-off (below 1e-12), while IAST closes the same loop to
#    round-off (below 1e-12). With equal
#    capacities the curl is symbolically zero, which is exactly the case where
#    extended Langmuir *is* an IAST.
# 4. **One structural result the paper stops short of.** Under IAST, the
#    selectivity coefficient is *exactly* constant whenever the two pure
#    isotherms are the same curve shifted in $\log P$, the Henry limit of
#    eqs. (30)–(35) is the linear special case, so the composition dependence
#    the paper predicts in Fig. 10 is a measure of the two isotherms' difference
#    in **shape**, not in scale. Verified to the nested solver's tolerance
#    (relative spread below 1e-10) on a dual-site pair the closed forms never see.
#
# Nothing on this page is a fit (there is no data to fit), and nothing on it is
# experimental. **This page establishes that the theory is what the paper says
# it is, not that it is empirically right.** The empirical case lives in figures
# this page does not touch.

# %% [markdown]
# ## Background
#
# ### The problem in 1965
#
# Adsorption equilibria of gas *mixtures* are what separation design actually
# needs, and measuring them is far more work than measuring pure-component
# isotherms. The obvious programme, write the mixture isotherm's parameters as
# functions of the pure-component parameters, had, in the paper's words, "not
# been very successful; the predictions have not been in quantitative agreement
# with the experimental data ($2$) and often not even in qualitative agreement
# ($1$)". The one alternative on the table, Arnold's liquid entropy model, the
# paper dismisses on two grounds it states plainly: "it is thermodynamically
# inconsistent and … separate numerical integrations are required for each vapor
# composition."
#
# Myers and Prausnitz instead import the machinery that had organised
# vapour–liquid equilibria for decades: **treat the adsorbed phase as a
# solution, define activity coefficients for it, and let ideality be an
# empirical question.** The move that makes this possible is recognising the
# **spreading pressure** $\pi$, the two-dimensional analogue of pressure, the
# intensive variable conjugate to adsorbent area $A$, as a full thermodynamic
# variable. There is "no experimental technique for measuring the spreading
# pressure directly", the paper concedes, but it is computable from any
# measured isotherm through the Gibbs adsorption isotherm, and that is all the
# theory needs.
#
# Three assumptions carry the thermodynamics (p. 121): the adsorbent is inert;
# it "possesses a temperature-invariant area which is the same for all
# adsorbates", explicitly *not* valid for molecular sieves, where accessible
# area depends on the adsorbate; and adsorption is what the Gibbs definition
# says it is. On top of these the derivation assumes a perfect gas phase, which
# the paper defends for pressures "usually less than 1 atm."
#
# ### What the paper found, and why it was surprising
#
# The original plan, the paper says, was to *measure* activity coefficients
# from mixture data and then interpret them. "Surprisingly, the calculated
# activity coefficients were found to be equal to unity within the experimental
# error", so the ideal adsorbed solution, the $\gamma_i = 1$ case, predicts
# the mixture equilibria outright, from pure-component isotherms alone. The
# conclusions call it "very surprising" a second time and are careful about
# why: the adsorbents in question are highly heterogeneous, no statistical
# model of them existed to compare against, and the paper explicitly declines
# to claim the result must generalise, mixtures of "highly dissimilar
# components … may show appreciable nonideality".
#
# ### Where this page's siblings sit
#
# The pure-component isotherms IAST consumes are the subject of `J1.1`
# (Langmuir 1918) and `J1.3` (BET 1938); either drops into this page's
# $F_i(P_i^{\circ})$ slot. Adsorption *dynamics*, breakthrough, the
# linear-driving-force law, is `J1.5`, and an IAST equilibrium can serve as
# that page's mixture isotherm. This page is equilibrium only.

# %% [markdown]
# ## The published model
#
# All equation numbers are the paper's, read from a 300 ppi native render.
#
# ### Thermodynamics of the adsorbed phase (eqs. 1–18)
#
# Substituting $-\pi\,dA$ for $-P\,dV$ work makes the adsorbed phase a
# two-dimensional fluid:
#
# $$dU = T\,dS - \pi\,dA + \Sigma \mu_i\,dn_i \tag{1}$$
# $$dG = -S\,dT + A\,d\pi + \Sigma \mu_i\,dn_i \tag{2}$$
#
# so $G = G(T, \pi, n_i)$ and, by Euler's theorem, $G = \Sigma n_i \mu_i$ (3).
# Activity coefficients are then *defined* exactly as for a liquid mixture, at
# constant $T$ **and $\pi$** (5)–(8), giving the chemical potential
#
# $$\mu_i(T,\pi,x_1\ldots) = g_i^{\circ}(T,\pi) + RT \ln \gamma_i x_i, \tag{9}$$
#
# and, writing the pure-adsorbate reference through its equilibrium pressure
# $P_i^{\circ}(\pi)$ at the same spreading pressure (10), then equating with the
# perfect-gas chemical potential (12), the **equation of equilibrium for
# mixed-gas adsorption**:
#
# $$P y_i = P_i^{\circ}(\pi)\,\gamma_i x_i \quad (\text{constant } T). \tag{13}$$
#
# Eq. (14) restates (13) with fugacities for high pressure. The Gibbs adsorption
# isotherm $-A\,d\pi + \Sigma n_i\,d\mu_i = 0$ (15) yields the Gibbs–Duhem
# relation $\Sigma x_i\,d\ln\gamma_i = 0$ at constant $T$ and $\pi$ (17), so the
# activity coefficients are thermodynamically consistent by construction, the
# failure the paper charges Arnold's model with cannot occur here. The phase
# rule for adsorption (18) gives binary adsorption **3** degrees of freedom, one
# more than binary VLE, because area is an extra variable; the adsorbent itself
# is not counted ("The absorbent is not counted as a component in Equation (18)
# since it is assumed to be thermodynamically inert.", *sic*, "absorbent",
# p. 122; the paper spells "adsorbent" correctly everywhere else).
#
# ### Two printed equation defects, reported and not repaired
#
# **Eq. (11) is missing its `+`.** As printed (p. 122, verified on a crop of the
# 300 ppi native render):
#
# $$\mu_i(T,\pi,x_1\ldots) = g_i^{\circ}(T) + RT \ln P_i^{\circ}(\pi)\;RT\ln\gamma_i x_i
# \qquad [sic] \tag{11}$$
#
# with no operator between the two $RT\ln$ terms. The product form is dimensionally
# incoherent (energy²) and §2.1 shows symbolically that equating it with
# eq. (12) does **not** produce eq. (13), while the sum form produces it
# exactly. Since (11) is by construction "Substituting Equation (10) into
# Equation (9)", and both (10) and (9) are printed correctly, the defect is
# proved from the paper's own equations; the repair (a lost `+`) is an
# inference, labelled as one.
#
# **The sentence introducing eq. (21) cites eq. (21).** P. 123, verbatim: "With
# the pressure $P$ of the mixture held constant, Equation (21) becomes", and
# what follows *is* eq. (21). The derivation requires "Equation (20) becomes":
# at constant $P$, $d\ln P y_i = d\ln y_i$ turns (20) into (21), which §2.1
# also verifies. Again the repair is an inference; the defect is not.
#
# ### The ideal adsorbed solution (eqs. 19–29)
#
# For $\gamma_i = 1$, eq. (13) is the Raoult analogue (22) quoted in the title
# cell. The spreading pressure of a pure component comes from integrating (15):
#
# $$\pi(P_i^{\circ}) = \frac{RT}{A}\int_{t=0}^{P_i^{\circ}} n_i^{\circ}(t)\,d\ln t
# \quad (\text{constant } T), \tag{19}$$
#
# well defined at the lower limit because $n_i^{\circ} \propto P_i^{\circ}$ at
# low coverage. Mixing at constant $T,\pi$ has no enthalpy and no area change
# ((23), (24)), and the area balance gives the **total amount adsorbed**
#
# $$\frac{1}{n_t} = \frac{x_1}{n_1^{\circ}} + \frac{x_2}{n_2^{\circ}}
# \quad (\text{constant } T \text{ and } \pi) \tag{26}$$
#
# is a harmonic mix, not an arithmetic one; §7.5 shows the arithmetic version
# fails the Gibbs isotherm. Relative volatility and selectivity are defined by
# (27), (28) with $s_{1,2} = 1/\alpha_{1,2} = P_2^{\circ}/P_1^{\circ}$, and
# (22) + (41) give the Raoult total-pressure line
# $P = P_1^{\circ}x_1 + P_2^{\circ}x_2$ (29).
#
# ### The calculation scheme (eqs. 36–46) and the special case that anchors it
#
# The p. 124 scheme is a nine-unknown, seven-equation algebraic system:
# $\pi_1^{\circ} = \pi_2^{\circ}$ (40) is the matching condition, closed by
# specifying two quantities, e.g. $P$ and $y_1$; the paper's own count
# (9 − 7 = 2) agrees with its phase rule (3 degrees of freedom, one spent on
# $T$), which §4 checks from the transcribed counts. When the pure isotherms
# are analytic, (44) is one equation in one unknown; in general the paper gives
# the graphical construction of Fig. 1, whose lever rules (45), (46) §2.1
# proves and §6.2 draws. At very low coverage, Henry's law, $n_i^{\circ} =
# K_i P_i^{\circ}$ (30), the paper proves IAST **rigorous**: eqs. (32)–(35)
# reproduce the statistical-mechanical mixture result $n_i = K_i P y_i$ (31)
# exactly, with $\pi A/RT = K_i P_i^{\circ} = n_1^{\circ} = n_2^{\circ}$ (32)
# and constant selectivity $K_1/K_2$ (34). §7.3 uses this as a solver test.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.20.4"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "J1.4-iast-mixed-gas-adsorption"  # equals the page directory

import os, sys, urllib.request
from pathlib import Path
_shared = next((p / "shared" for p in (Path.cwd(), *Path.cwd().parents)
                if (p / "shared" / "gallery_bootstrap.py").is_file()), None)
if _shared is None:
    _shared = (Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
               / "pymrm-gallery" / GALLERY_REF)
    _shared.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/computational-chemical-engineering/"
        f"pymrm-gallery/{GALLERY_REF}/shared/gallery_bootstrap.py",
        _shared / "gallery_bootstrap.py")
sys.path.insert(0, str(_shared))
sys.modules.pop("gallery_bootstrap", None)
from gallery_bootstrap import setup
setup(PAGE, PYMRM, GALLERY_REF);  # the semicolon stops the cell echoing the path

# %%
# Imports. numpy, scipy, pandas, sympy and matplotlib are Colab defaults.
import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss
from pymrm.checks import observed_orders

from gallery_utils import Claims, cite_data, load_data, load_meta
from model import (TOL, circulation, dilute_selectivities, dualsite_pure, eqm_dev,
                   ext_langmuir, henry_pure, iast, langmuir_pure, lever_rule,
                   n1_dilute_dev, nt_dev_sup, path_rhs_gauss, path_rhs_trap, rel_spread,
                   z_quad)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
pd.set_option("display.width", 130)
claims = Claims(PAGE.split("-")[0])
EPS = np.finfo(float).eps

# %% [markdown]
# ### 2.1 The derivation chain, verified symbolically
#
# Eleven identities, every one an equation the paper prints (or, for the last
# three, a closed form this page derives and then leans on). A nonzero
# residual raises and fails the notebook. These checks verify the paper's
# algebra and this page's closed forms; they cannot detect a numerical defect,
# and the numerical machinery gets its own, breakable tests in §7. The count of
# identities has one break row of its own: eq. (11) read as printed.

# %%
# Symbols. All positive reals; gamma_i carried where the paper carries it.
RT, g0, P, y1, y2, x1, x2, P1o, P2o, gam, K1, K2, m, b, t, Pv = sp.symbols(
    "RT g0 P y1 y2 x1 x2 P1o P2o gamma K1 K2 m b t Pv", positive=True)
ZERO = []

# --- (11)+(12) -> (13): the equilibrium relation, with eq. (11) read as a SUM
mu_gas = g0 + RT*sp.log(P*y1)                        # eq. (12)


def eq11_to_13(form="sum"):
    """Residual of (13) from (11) and (12); eq. (11) read as a sum or as printed."""
    if form == "sum":
        mu_ads = g0 + RT*sp.log(P1o) + RT*sp.log(gam*x1)     # eq. (11), '+' restored
    else:
        mu_ads = g0 + RT*sp.log(P1o)*RT*sp.log(gam*x1)       # eq. (11) verbatim [sic]
    res13 = sp.simplify(sp.exp((mu_gas - g0)/RT) - sp.exp((mu_ads - g0)/RT))
    return sp.simplify(res13 - (P*y1 - P1o*gam*x1))


ZERO.append(("(11)+(12) -> (13):  P y1 - P1o gamma x1", eq11_to_13("sum")))

# --- and AS PRINTED (no operator: a product) it does NOT give (13)
mu_printed = g0 + RT*sp.log(P1o)*RT*sp.log(gam*x1)   # eq. (11) verbatim [sic]
Py_printed = sp.simplify(sp.exp((mu_printed - g0)/RT))
print("eq. (11) as printed (product) would give  P y1 =", Py_printed)
print("  -> not eq. (13), and dimensionally incoherent (the product term is energy^2).\n")

# --- (22)+(41)+(42) -> (29): the Raoult total-pressure line
sol29 = sp.solve([sp.Eq(P*y1, P1o*x1), sp.Eq(P*y2, P2o*x2),
                  sp.Eq(y1 + y2, 1), sp.Eq(x1 + x2, 1)], [P, y1, y2, x2], dict=True)
ZERO.append(("(22)+(41)+(42) -> (29):  P = P1o x1 + P2o x2",
             sp.simplify(sol29[0][P] - (P1o*x1 + P2o*(1 - x1)))))

# --- (20) -> (21) at constant P:  d ln(P y_i) = d ln y_i  (the misprinted cross-reference)
lnPy = sp.log(P*y1)
ZERO.append(("(20) -> (21):  d ln(P y1)/d y1 = d ln(y1)/d y1 at constant P",
             sp.simplify(sp.diff(lnPy, y1) - sp.diff(sp.log(y1), y1))))

# --- (19) on a Langmuir isotherm:  z = pi A/RT = m ln(1 + b P)
z_langmuir_sym = sp.integrate(m*b/(1 + b*t), (t, 0, Pv))   # integrand n(t)/t * t = n(t) d ln t
ZERO.append(("(19) on Langmuir:  z = m ln(1+bP)", sp.simplify(z_langmuir_sym - m*sp.log(1 + b*Pv))))

# --- (30) into (19) -> (32):  z = K_i P_i^o  in the Henry limit
z_henry_sym = sp.integrate(K1*t/t, (t, 0, Pv))
ZERO.append(("(30) into (19) -> (32):  z = K P", sp.simplify(z_henry_sym - K1*Pv)))

# --- (32)-(35) -> (31):  the Henry-limit IAST recovers n_i = K_i P y_i exactly
z_mix = P*(y1*K1 + (1 - y1)*K2)          # from x1+x2=1 with x_i = P y_i K_i / z
x1_h = P*y1*K1/z_mix
n1_h = x1_h*z_mix                        # n_t = z in this limit (eq. 32 + eq. 26)
ZERO.append(("(32)-(35) -> (31):  n1 = K1 P y1", sp.simplify(n1_h - K1*P*y1)))

for name, r in ZERO:
    print(f"  0 == {r}   [{name}]")
    if r != 0:
        raise RuntimeError(f"symbolic identity does not close: {name}")

# %%
# --- (45), (46): the lever rules of Fig. 1, from (22), (26), (29) alone.
# A = (0, z*), B = (P1o, z*), C = (P2o, z*) on the two pure curves at matched z*;
# D = (P, z*); F and E are the intersections of the vertical through D with the
# chords OB and OC.
zst, Ps = sp.symbols("zst Ps", positive=True)
x1_lr = (P2o - Ps)/(P2o - P1o)                 # eq. (29) solved for x1
y1_lr = P1o*x1_lr/Ps                           # eq. (22)
DE = zst - zst*Ps/P2o                          # D to E (E on chord OC)
FE = zst*Ps/P1o - zst*Ps/P2o                   # F to E (F on chord OB)
lever = []
lever.append(("(45):  y1 = DE/FE", sp.simplify(DE/FE - y1_lr)))
lever.append(("(46):  x1 = DC/BC", sp.simplify((P2o - Ps)/(P2o - P1o) - x1_lr)))

# --- equal-capacity Langmuir: the IAST solution IS extended Langmuir (this
#     page's closed form; not stated in the paper). Matching m ln(1+b1 P1o) =
#     m ln(1+b2 P2o) forces b1 P1o = b2 P2o = u, and x1+x2 = 1 gives
#     u = b1 P y1 + b2 P y2; then eq. (26) collapses:
b1s, b2s = sp.symbols("b1 b2", positive=True)
u = b1s*P*y1 + b2s*P*(1 - y1)
x1_eq = b1s*P*y1/u                              # x1 = P y1/P1o with P1o = u/b1
nt_eq = m*u/(1 + u)                             # n_i^o(P_i^o) = m u/(1+u), both components
n1_eq = sp.simplify(x1_eq*nt_eq)
extL = m*b1s*P*y1/(1 + b1s*P*y1 + b2s*P*(1 - y1))
lever.append(("equal-m IAST == extended Langmuir (n1)", sp.simplify(n1_eq - extL)))

# --- the curl of the extended-Langmuir Gibbs one-form, unequal capacities.
#     A d(pi)/RT = n_t d lnP + [n1/y1 - n2/y2] d y1 must be EXACT for pi to be
#     a state function. sympy derives the curl in closed form:
m1s, m2s, lnP = sp.symbols("m1 m2 lnP", real=True)
Pe = sp.exp(lnP)
u1, u2 = b1s*Pe*y1, b2s*Pe*(1 - y1)
D = 1 + u1 + u2
M = (m1s*u1 + m2s*u2)/D                         # coefficient of d lnP  (= n_t)
N = (m1s*b1s*Pe - m2s*b2s*Pe)/D                 # coefficient of d y1   (= n1/y1 - n2/y2)
curl = sp.simplify(sp.diff(N, lnP) - sp.diff(M, y1))
curl_closed = (m2s - m1s)*b1s*b2s*Pe**2/D**2
lever.append(("extended-Langmuir curl == (m2-m1) b1 b2 P^2 / D^2", sp.simplify(curl - curl_closed)))
lever.append(("... and == 0 at m1 = m2 (the consistent case)", sp.simplify(curl.subs(m2s, m1s))))

for name, r in lever:
    print(f"  0 == {r}   [{name}]")
    if r != 0:
        raise RuntimeError(f"symbolic identity does not close: {name}")
ZERO += lever
SYM_ZERO_COUNT = sum(r == 0 for _, r in ZERO)
claims.exact("sym_zero_identity_count", SYM_ZERO_COUNT,
             note="sympy residuals exactly zero; verifies algebra, not numerics")
print(f"\n{SYM_ZERO_COUNT} symbolic identities, all exactly zero. The curl closed form is the")
print("engine of §7.6: it vanishes iff m1 = m2, so the Langmuir mixture rule is")
print("thermodynamically consistent exactly when it happens to BE an ideal adsorbed solution.")
curl_fn = sp.lambdify((lnP, y1, m1s, b1s, m2s, b2s), curl, "numpy")

# %% [markdown]
# ## Parameters and assumptions
#
# **The paper's assumptions, inherited unchanged:** inert adsorbent;
# temperature-invariant area equal for all adsorbates; Gibbs definition of
# adsorption; perfect gas phase ($\gamma_i$ in eq. (13) set to 1, the ideal
# adsorbed solution; eq. (14)'s fugacity form is transcribed in §2 and **not
# exercised**, because the paper prints no real-gas model to put in it).
#
# **Illustrative parameters, and what they are not.** The paper tabulates no
# isotherm, not even the pure-component inputs of its own Fig. 2, so every
# number this page computes with is an **illustrative parameter set, chosen for
# no physical system, in arbitrary pressure and amount units**. No conclusion
# about any real gas pair is drawn from them anywhere on this page; they exist
# so that the theory's printed structure can be exercised and broken.
#
# | set | component 1 | component 2 | used for |
# |---|---|---|---|
# | **pair L** (equal capacities) | $m{=}1$, $b_1{=}10$ | $m{=}1$, $b_2{=}1$ | the closed-form solver test (§7.1); constant-selectivity demonstration |
# | **pair U** (unequal capacities) | $m_1{=}1$, $b_1{=}10$ | $m_2{=}0.5$, $b_2{=}2$ | everything else: Results, consistency loop, Henry limit |
# | **pair S** (dual-site, shape-shifted) | two-site Langmuir | the *same curve* at $\kappa P$, $\kappa = 1/4$ | the shape-translation result (§7.8) |
#
# with pure isotherms $n^{\circ}(P) = m\,bP/(1+bP)$ (single site), for which
# eq. (19) integrates to $z \equiv \pi A/RT = m\ln(1+bP)$, the closed form §2.1
# verified, and analytic inverses $P^{\circ}(z) = (e^{z/m}-1)/b$. All solver
# tolerances are $10^{-13}$ on the spreading-pressure match.

# %% [markdown]
# ## The data
#
# **There is none, and that is this page's central scope fact.** Stated
# precisely, with the search that backs each claim (all seven pages were read
# on 300 ppi native renders; the text layer was used only as a search index,
# never for a digit):
#
# - **The paper contains no tables.** The only tabular matter in its seven pages
#   is the NOTATION list and LITERATURE CITED (p. 126–127). No numerical result
#   of the theory, no predicted composition, amount, or selectivity, is
#   printed anywhere in the text.
# - **All four experimental comparisons are figures**: Figs. 3–4
#   (methane–ethane, activated carbon, 20 °C, 1 atm), Figs. 5–6 (CO–O₂, silica,
#   100 °C and 0 °C, 1 atm), Fig. 7 (propane–propylene, silica gel, 25 °C,
#   1 atm), Figs. 8–10 (ethylene–CO₂, activated carbon, 25.4 °C, 50–250 mm Hg).
#   The experimental points belong to Szepesy & Illés (1963), Markham & Benton
#   (1931), Lewis et al. (1950) and Bering & Serpenskii (1952) respectively:
#   none of which is on disk here.
# - **Even the pure-component inputs are not tabulated.** Fig. 2's legend box
#   reads "(CALCULATED FROM EXPERIMENTAL ADSORPTION ISOTHERMS OF SZEPESY &
#   ILLES (1963))": the paper's own inputs exist in it only as
#   already-integrated $\pi A/RT$ curves. Reproducing Fig. 3 would therefore
#   mean **differentiating a digitised integral** of data from a paper that is
#   not on disk, figure-derived inputs feeding figure-derived targets. That
#   double dependence, plus the repository rule that figure content needs a
#   human review gate, is why the experimental case is scoped out rather than
#   digitised. If the Szepesy & Illés tables (Acta Chim. Hung. **35**, 37, 53,
#   245) are ever acquired, the four comparisons become a natural companion
#   page.
#
# What the paper *does* print as numbers, conditions, the text's own counting,
# and the reference list, is transcribed in `myers-1965-printed-claims.csv`
# (the count is printed below, each scalar read off a digit-scale crop of the native render), and the
# page's computed IAST solutions are written to
# `iast-illustrative-reference.csv` so an independent implementation can
# regression-test against them. **No other page's dataset is loaded**, so none
# of the cross-page reconciliation obligations apply.

# %%
printed = load_data("myers-1965-printed-claims.csv", page=PAGE)
C = dict(zip(printed.key, pd.to_numeric(printed.printed)))
print(cite_data(load_meta("myers-1965-printed-claims.csv", page=PAGE)))
print(f"{len(printed)} printed scalars transcribed.\n")


# --- the paper's own counting, checked against itself (all numbers from the CSV)
def phase_rule_dof(components=2, phases=2):
    """Eq. (18): degrees of freedom of adsorption, components - phases + 3."""
    return components - phases + 3


DOF = phase_rule_dof()
if DOF != C["text_dof_binary"]:
    raise RuntimeError(f"eq. (18) gives {DOF}, the text prints {C['text_dof_binary']:.0f}")
free_after_T = C["text_dof_binary"] - 1      # T is fixed in the p.124 scheme
if not C["text_unknowns"] - C["text_equations"] == C["text_specified"] == free_after_T:
    raise RuntimeError("the p. 124 counts disagree with the phase rule")
claims.exact("eq18_dof_binary", DOF, note="eq. (18) for a binary on an inert adsorbent")
print(f"eq. (18) for a binary: 2 - 2 + 3 = {DOF} degrees of freedom (printed: {C['text_dof_binary']:.0f})")
print(f"p. 124 scheme: {C['text_unknowns']:.0f} unknowns - {C['text_equations']:.0f} equations = "
      f"{C['text_unknowns'] - C['text_equations']:.0f} to specify = printed {C['text_specified']:.0f}"
      f" = DOF {C['text_dof_binary']:.0f} with T already fixed.  Internally consistent.")

# --- two bibliographic defects, proved from the reference list's own numbers
print(f"\nFig. 6 legend attributes its points to 'MARKHAM & BENTON ({C['fig6_legend_year']:.0f})' [sic],")
print(f"but the paper's ONLY Markham & Benton reference (ref. 8) is J. Am. Chem. Soc. "
      f"{C['ref8_volume']:.0f}, 497 ({C['ref8_year']:.0f}),")
print(f"and the p. 125 text credits the 0 C data to '(8)'. {C['ref7_year']:.0f} is ref. 7, Lewis et al.")


def ref5_volume_rate(c):
    """Volumes per year J. Chem. Phys. would need between ref. 4 and ref. 5."""
    return (c["ref5_volume"] - c["ref4_volume"]) / (c["ref5_year"] - c["ref4_year"])


REF5_VOL_RATE = ref5_volume_rate(C)
claims.value("ref5_volumes_per_year", REF5_VOL_RATE, rtol=1e-9)
print(f"\nRef. 5 reads 'Ibid, {C['ref5_volume']:.0f}, 456 ({C['ref5_year']:.0f})' where Ibid = "
      f"J. Chem. Phys. (ref. 4: vol. {C['ref4_volume']:.0f} in {C['ref4_year']:.0f}).")
print(f"That requires {REF5_VOL_RATE:.1f} volumes/year from the same journal, and J. Chem. Phys. was at "
      f"vol. 37 in 1962;")
print("J. Phys. Chem.'s vol. 63 IS 1959. That ref. 5 points at the wrong journal is an inference;")
print("the arithmetic tension is not. Ref. 5 is also cited nowhere in the running text: every '(5)'")
print("in the text layer (searched as an index, then verified on crops) is equation (5).")


# %% [markdown]
# ## PyMRM implementation
#
# The p. 124 scheme, eqs. (36)–(42), reduced the way the paper's own eq. (44)
# reduces it: one unknown per state point, the reduced spreading pressure
# $z = \pi A/RT$, with the matching condition and mole-fraction closure
#
# $$R(z) \;=\; \frac{P y_1}{P_1^{\circ}(z)} + \frac{P y_2}{P_2^{\circ}(z)} - 1 \;=\; 0,$$
#
# solved simultaneously for a whole grid of $(P, y_1)$ states as a **block of
# independent scalar Newton problems**: `NumJac((K, 1))` (the shape is
# `(K, 1)`, never `(K,)`: a bare 1-D shape would declare all $K$ states coupled
# and build a dense Jacobian) and pymrm's `newton` on $w = \ln z$, which keeps
# the iterate positive, as $P^{\circ}(z)$ needs ($P^{\circ}(0) = 0$). Every
# solve raises unless `newton` reports success. Everything downstream is the
# paper's algebra: $x_1$ from (22), $n_t$ from (26), $s_{1,2}$ from (28).
#
# Two routes to $z(P^{\circ})$ exist and both are used: the **analytic** one for
# Langmuir forms ($z = m\ln(1+bP^{\circ})$, verified symbolically in §2.1), and
# the **general quadrature** of eq. (19) for isotherms known only pointwise:
# trapezoid on a logarithmic grid with the Henry-law segment $[0, t_0]$
# integrated analytically as $K t_0$, exactly the low-coverage care the paper's
# conclusions ask for ("the integration for spreading pressure is sensitive to
# this portion of the … isotherm"). §7.4 refines it and breaks it.

# %%
# The isotherm objects, the IAST solver and the rival rule live in model.py:
# langmuir_pure, henry_pure, dualsite_pure (its inverse a nested newton solve),
# iast (newton on NumJac((K, 1)) in w = ln z, raising unless it converges) and
# ext_langmuir (closed form, sharing nothing with the solver).

# --- the three illustrative pairs ---------------------------------------
mL, b1L, b2L = 1.0, 10.0, 1.0
m1U, b1U, m2U, b2U = 1.0, 10.0, 0.5, 2.0
LANG_L, LANG_U = (mL, b1L, mL, b2L), (m1U, b1U, m2U, b2U)   # extended-Langmuir constants
pL1, pL2 = langmuir_pure(mL, b1L), langmuir_pure(mL, b2L)
pU1, pU2 = langmuir_pure(m1U, b1U), langmuir_pure(m2U, b2U)
KAPPA = 0.25
pS1 = dualsite_pure(0.7, 20.0, 0.3, 0.5, 1.0)
pS2 = dualsite_pure(0.7, 20.0, 0.3, 0.5, KAPPA)   # the SAME curve, shifted in log P

demo = iast(1.0, np.array([0.5]), pU1, pU2)
print("pair U at P = 1, y1 = 0.5:")
print(f"  z = {demo['z'][0]:.6f}   x1 = {demo['x1'][0]:.6f}   n_t = {demo['nt'][0]:.6f}"
      f"   s12 = {demo['s12'][0]:.4f}")

# %%
# The construction the paper computes its Fig. 2 with: n(P) and z(P) for pair U.
# (Same CONSTRUCTION as the paper's figure; the curves are this page's
#  illustrative isotherms, not the paper's systems.)
Pplot = np.geomspace(1e-3, 20, 400)
fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
ax[0].semilogx(Pplot, pU1["n"](Pplot), "-",  color="C0", label="component 1  ($m_1{=}1,\\ b_1{=}10$)")
ax[0].semilogx(Pplot, pU2["n"](Pplot), "--", color="C1", label="component 2  ($m_2{=}0.5,\\ b_2{=}2$)")
ax[0].set(xlabel="$P$ (arbitrary units)", ylabel="$n^{\\circ}(P)$ (arbitrary units)",
          title="Pure isotherms, pair U (illustrative)")
ax[0].legend(fontsize=8)
ax[1].semilogx(Pplot, pU1["z"](Pplot), "-",  color="C0", label="component 1")
ax[1].semilogx(Pplot, pU2["z"](Pplot), "--", color="C1", label="component 2")
ax[1].set(xlabel="$P$ (arbitrary units)", ylabel="$z = \\pi A/RT$ (amount units)",
          title="Reduced spreading pressure, eq. (19)")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()

# %% [markdown]
# ## Results
#
# ### 6.1 The mixture equilibria (the construction of Figs. 3, 4, illustrative)
#
# The $(y{-}x)$ diagram and the total amount adsorbed for pair U, at three total
# pressures. These have the *structure* of the paper's Figs. 3 and 4, an
# azeotrope-free $y{-}x$ curve bowed toward the strongly adsorbed component,
# $n_t$ falling as the weakly adsorbed component takes over the gas, but they
# are computed from the illustrative parameters and **reproduce nothing in the
# paper**; the paper's own curves would need the Szepesy & Illés isotherms this
# page does not have.

# %%
yfine = np.linspace(1e-3, 1-1e-3, 399)
fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
for Pt, ls in [(0.2, ":"), (1.0, "-"), (5.0, "--")]:
    s = iast(Pt, yfine, pU1, pU2)
    ax[0].plot(s["x1"], yfine, ls, color="C0", label=f"$P = {Pt:g}$")
    ax[1].plot(yfine, s["nt"], ls, color="C0", label=f"$P = {Pt:g}$")
ax[0].plot([0, 1], [0, 1], lw=0.8, color="0.5")
ax[0].set(xlabel="$x_1$ (adsorbed phase)", ylabel="$y_1$ (gas phase)",
          title="$y{-}x$ diagram, pair U (illustrative)")
ax[0].legend(fontsize=8)
ax[1].set(xlabel="$y_1$ (gas phase)", ylabel="$n_t$ (amount units)",
          title="Total amount adsorbed, eq. (26)")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()
S_RESULTS = {Pt: iast(Pt, yfine, pU1, pU2) for Pt in (0.2, 1.0, 5.0)}

# %% [markdown]
# ### 6.2 The paper's graphical procedure, drawn and checked (Fig. 1's construction)
#
# The paper's general algorithm is geometric: at a chosen spreading pressure the
# horizontal line $ABC$ cuts the two pure $z(P)$ curves at $B = (P_1^{\circ}, z^{*})$
# and $C = (P_2^{\circ}, z^{*})$; a point $D$ on $BC$ picks the mixture pressure;
# the chords $OB$, $OC$ and the vertical through $D$ give $E$ and $F$; and
# eqs. (45), (46) read the equilibrium off the ruler:
# $y_1 = \mathrm{DE}/\mathrm{FE}$, $x_1 = \mathrm{DC}/\mathrm{BC}$. §2.1 proved
# both symbolically; here the construction is drawn for pair U and the lever
# rules evaluated numerically, the same numbers the Newton solver returns,
# which is eq. (44) and the ruler agreeing about eq. (40).

# %%
ZSTAR, X1STAR = 1.2, 0.85
LV = lever_rule(pU1, pU2, ZSTAR, X1STAR)          # eqs. (45), (46) against (22), (29)
P1s, P2s, Pd, Fh, Eh = LV["P1"], LV["P2"], LV["Pd"], LV["F"], LV["E"]
y1_lever, x1_lever = LV["y1"], LV["x1"]
LEVER_DEV = LV["dev"]
claims.bound("lever_rule_max_dev", LEVER_DEV, max=1e-12, basis="identity",
             note="lever readings against eqs. (22), (29); exact in exact arithmetic")
# and the defect that shows the check has power: swap the chords
LEVER_SWAPPED = lever_rule(pU1, pU2, ZSTAR, X1STAR, swap=True)["dev"]

Pc = np.geomspace(1e-3, P2s*1.08, 300)
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.plot(Pc, pU1["z"](Pc), color="C0", label="pure component 1")
ax.plot(Pc, pU2["z"](Pc), color="C1", label="pure component 2")
ax.plot([0, P2s], [ZSTAR, ZSTAR], "k--", lw=0.8)
ax.plot([0, Pd], [0, Fh], "k:", lw=0.8)                 # chord OB extended to F
ax.plot([0, P2s], [0, ZSTAR], "k:", lw=0.8)             # chord OC
ax.plot([Pd, Pd], [0, Fh], "k--", lw=0.8)
for x_, z_, nm, off in [(0, ZSTAR, "A", (-11, 4)), (P1s, ZSTAR, "B", (-3, 6)),
                        (P2s, ZSTAR, "C", (4, 4)), (Pd, ZSTAR, "D", (4, 5)),
                        (Pd, Eh, "E", (5, -11)), (Pd, Fh, "F", (5, 1))]:
    ax.plot(x_, z_, "ko", ms=4)
    ax.annotate(nm, (x_, z_), textcoords="offset points", xytext=off, fontsize=10)
ax.set(xlabel="$P$ (arbitrary units)", ylabel="$z = \\pi A/RT$",
       title=f"Fig. 1's construction, computed: $y_1 = $ DE/FE $ = {y1_lever:.6f}$,"
             f"  $x_1 = $ DC/BC $ = {x1_lever:.6f}$",
       xlim=(0, P2s*1.08), ylim=(0, Fh*1.12))
ax.legend(fontsize=8, loc="center right")
fig.tight_layout(); plt.show()
print(f"lever rules vs eqs. (22)/(29) directly: max deviation below 1e-12 (here {LEVER_DEV:.1e})")
print(f"with the chords deliberately swapped:   {LEVER_SWAPPED:.4f}  (the check can fail)")

# %% [markdown]
# ### 6.3 Selectivity: what varies, what cannot (the structure of Fig. 10's claim)
#
# The paper's Fig. 10 discussion is the one place it exhibits IAST predicting a
# *composition-dependent* selectivity, "a small composition dependence" the
# data of Bering and Serpenskii could not resolve. For pair U the dependence is
# anything but small, and both infinite-dilution ends have **closed forms** (the
# dilute-1 limit is $s \to P/P_1^{\circ}(z_2(P))$, the dilute-2 limit
# $s \to P_2^{\circ}(z_1(P))/P$; at $P = 1$ these evaluate to $5(1+\sqrt3) =
# 13.660254$ and exactly $60$), so the page reports limits, not sampled
# extrema. For pair L the selectivity is $b_1/b_2 = 10$ **identically**, under
# IAST, any two pure isotherms that are the same curve shifted in $\log P$ give
# exactly constant selectivity (§7.8 proves and breaks this), so composition
# dependence of $s$ under IAST measures the two isotherms' difference in
# *shape*, not in scale. The paper does not state this; its Henry-limit section
# (eqs. 30–35, constant $s = K_1/K_2$) is the linear special case of it.

# %%
# closed-form dilution limits at P = 1 (pair U)
S_DIL1, S_DIL2 = dilute_selectivities(1.0, pU1, pU2)   # y1 -> 0 and y1 -> 1
S_DIL_RATIO = S_DIL2/S_DIL1
claims.value("selectivity_dilute1_P1", S_DIL1, rtol=1e-9, note="closed form, 5(1+sqrt 3)")
claims.value("selectivity_dilute2_P1", S_DIL2, rtol=1e-9, note="closed form, exactly 60")
claims.value("selectivity_dilute_ratio_P1", S_DIL_RATIO, rtol=1e-9)
fig, ax = plt.subplots(figsize=(7.2, 4.0))
for Pt, ls in [(0.2, ":"), (1.0, "-"), (5.0, "--")]:
    s = S_RESULTS[Pt]
    ax.plot(yfine, s["s12"], ls, color="C0", label=f"pair U, $P = {Pt:g}$")
ax.plot(1e-3, S_DIL1, "o", color="C3", ms=6, zorder=5)
ax.plot(1-1e-3, S_DIL2, "o", color="C3", ms=6, zorder=5,
        label="closed-form dilution limits, $P{=}1$")
sL = iast(1.0, yfine, pL1, pL2)
ax.plot(yfine, sL["s12"], "-", color="C2", lw=1.2,
        label="pair L (equal shape): $s \\equiv b_1/b_2 = 10$")
ax.set(xlabel="$y_1$ (gas phase)", ylabel="selectivity $s_{1,2}$, eq. (28)", yscale="log",
       title="Selectivity under IAST: shape difference, not scale, drives the variation")
ax.legend(fontsize=8)
fig.tight_layout(); plt.show()
print(f"pair U at P = 1: dilute-in-1 limit s = {S_DIL1:.6f} = 5(1+sqrt(3)) = {5*(1+np.sqrt(3)):.6f}")
print(f"                 dilute-in-2 limit s = {S_DIL2:.6f} (exactly 60: e^(2 ln 11) = 121, (121-1)/2)")
print(f"                 composition-dependence factor {S_DIL_RATIO:.6f}")


def pair_L_selectivity(pure2=pL2):
    """Mean and relative spread of pair L's selectivity at P = 1 over the composition range."""
    s_ = iast(1.0, yfine, pL1, pure2)["s12"]
    return dict(mean=float(np.nanmean(s_)), spread=rel_spread(s_))


SEL_L = pair_L_selectivity()
claims.value("pair_L_selectivity_mean", SEL_L["mean"], rtol=1e-9, note="b1/b2 for equal capacities")
claims.bound("pair_L_selectivity_rel_spread", SEL_L["spread"], max=1e-10, basis="tolerance",
             note="constant by the shape-translation identity; newton tolerance 1e-13")
print(f"pair L: mean selectivity {SEL_L['mean']:.12f}, relative spread over the whole "
      f"composition range below 1e-10 (here {SEL_L['spread']:.1e})")

# %%
# The page's computed reference table, written for reuse (see the data sidecar:
# COMPUTED output, not data from any paper). Deterministic across runs.
rows = []
for pair, (q1, q2) in [("L", (pL1, pL2)), ("U", (pU1, pU2))]:
    for Pt in (0.2, 1.0, 5.0):
        yv = np.round(np.arange(0.1, 0.91, 0.1), 10)
        s = iast(Pt, yv, q1, q2)
        for i, yy in enumerate(yv):
            rows.append((pair, Pt, yy, s["z"][i], s["x1"][i], s["nt"][i],
                         s["n1"][i], s["n2"][i], s["s12"][i]))
ref = pd.DataFrame(rows, columns=["pair", "P", "y1", "z", "x1", "n_t", "n_1", "n_2", "s12"])
Path("data").mkdir(exist_ok=True)
hdr = ("# IAST reference solutions COMPUTED BY PAGE J1.4 on two illustrative Langmuir pairs.\n"
       "# NOT experimental data and NOT from Myers & Prausnitz (1965); see the .meta.yaml sidecar.\n"
       "# pair L: m=1, b1=10, b2=1.  pair U: m1=1, b1=10, m2=0.5, b2=2.  Arbitrary units.\n")
with open("data/iast-illustrative-reference.csv", "w", encoding="utf-8") as fh:
    fh.write(hdr)
    ref.to_csv(fh, index=False, float_format="%.12g")
print(f"wrote data/iast-illustrative-reference.csv  ({len(ref)} rows)")
ref.head(6)

# %% [markdown]
# ## Validation
#
# Eight numerical sections. Every number the page states is a registered
# claim, and §7.9's defect-injection table has a row that fails each one. Since
# nothing here is experimental, **every check below is a consistency test of
# the theory's printed structure and of this page's solver, none is evidence
# that IAST describes any real system.**
#
# ### 7.1 The solver against a closed form it never touches (equal capacities)
#
# §2.1 derived symbolically that for equal-capacity Langmuir pairs the IAST
# solution *is* the extended Langmuir rule. The Newton solver knows nothing of
# that closed form, it solves the spreading-pressure match by iteration, so
# agreement is a real test of the whole chain (matching, eq. (22), eq. (26)).

# %%
Pgrid, ygrid = np.meshgrid(np.geomspace(0.01, 10, 30), np.linspace(0.02, 0.98, 30))
GP, GY = Pgrid.ravel(), ygrid.ravel()

EQM_DEV = eqm_dev(pL1, pL2, LANG_L, GP, GY)
claims.bound("equal_m_closed_form_max_rel", EQM_DEV, max=1e-12, basis="tolerance",
             note="newton step tolerance 1e-13 on ln z")
print(f"IAST solver vs extended-Langmuir closed form, pair L, {GP.size} states: "
      f"max rel dev below 1e-12 (here {EQM_DEV:.1e})")
print("(the solver's tolerance; the break rows of 7.9 move it by orders of magnitude)")

# %% [markdown]
# ### 7.2 Where extended Langmuir stops being an IAST (unequal capacities)
#
# With $m_1 \ne m_2$ the two models genuinely differ, and the size of the
# difference is this page's first physical statement. Both extrema are
# computed honestly: the deviation grows monotonically with $P$ across the
# window (checked below: the study raises otherwise), so its supremum sits on
# the window's $P = 10$ edge and is a **window supremum by definition, not a
# sampled maximum**; the interior maximum over $y_1$ at that edge is a bounded
# maximisation refined from a scan of the whole $y_1$ range, raising when the
# scan maximum sits on a scan end (the 30-by-30 grid read of the same quantity
# is 0.7 % low, printed beside it); and the dilute-component
# deviation has a **closed form** at $y_1 \to 0$, where both models are linear
# in $y_1$ and the ratio of slopes is
# $\bigl[m_1 b_1 P/(1{+}b_2 P)\bigr] \big/ \bigl[(P/P_1^{\circ}(z_2(P)))\,n_2^{\circ}(P)\bigr]$.

# %%
Y_EDGE = 1e-3                                  # the y1 scan stays this far inside (0, 1)
N_SCAN = 41                                    # nodes of the y1 scan
P_WINDOW = np.array([0.5, 1.0, 2.0, 5.0, 10.0])   # pressures the growth with P is checked at


def window_sup(pure2=pU2, lang=LANG_U, refine=True):
    """Supremum over the (P, y1) window of the n_t deviation of extended Langmuir from IAST.

    The maximum over y1 at each pressure of P_WINDOW must grow with P, so the
    supremum is the maximum at the P = 10 edge; raises otherwise.
    """
    per_P = [nt_dev_sup(Pt, pU1, pure2, lang, N_SCAN, Y_EDGE, refine=refine) for Pt in P_WINDOW]
    sups = np.array([v for v, _ in per_P])
    if not np.all(np.diff(sups) > 0):
        raise RuntimeError("the maximum over y1 does not grow with P across the window")
    sup, ystar = max(per_P)                    # the P = 10 edge, by the check above
    grid = iast(GP, GY, pU1, pure2)
    n1e, n2e = ext_langmuir(GP, GY, *lang)
    grid_read = float(np.max(np.abs((n1e + n2e) - grid["nt"])/grid["nt"]))
    return dict(sup=sup, y=ystar, grid=grid_read, grid_low=1 - grid_read/sup)


WIN = window_sup()
NT_DEV_SUP, NT_DEV_Y, NT_DEV_GRID = WIN["sup"], WIN["y"], WIN["grid"]
claims.value("unequal_m_nt_dev_windowsup", NT_DEV_SUP, rtol=1e-6)
claims.value("unequal_m_nt_dev_y_at_sup", NT_DEV_Y, rtol=1e-5)
claims.value("unequal_m_nt_dev_grid_low_frac", WIN["grid_low"], rtol=1e-4,
             note="the 30 x 30 grid read against the refined supremum")

# dilute-component deviation, closed form at y1 -> 0, P = 10
N1_DILUTE_DEV = n1_dilute_dev(10.0, pU1, pU2, m1U, b1U, b2U)
claims.value("unequal_m_n1_dilute_limit_dev", N1_DILUTE_DEV, rtol=1e-9)

print(f"n_t deviation, window sup (P = 10 edge, y1 refined from a {N_SCAN}-node scan): "
      f"{NT_DEV_SUP:.6f} at y1* = {NT_DEV_Y:.6f}")
print(f"  the 30 x 30 grid read of the same quantity: {NT_DEV_GRID:.6f}  "
      f"({100*WIN['grid_low']:.2f} % low, why extrema are refined here)")
print(f"dilute-component-1 deviation at P = 10, closed form: {N1_DILUTE_DEV:.6f}")

# %% [markdown]
# ### 7.3 The Henry-law limit, the one case the paper proves rigorous
#
# Eqs. (30)–(35): at vanishing coverage IAST must reproduce $n_i = K_i P y_i$
# and the constant selectivity $K_1/K_2$. Two tests. **(a)** On *pure Henry
# isotherms* the solver must be exact at any pressure, and it is, to solver
# tolerance. **(b)** On pair U the selectivity must converge to
# $K_1/K_2 = m_1 b_1 / m_2 b_2 = 10$ as $P \to 0$, at first order in $P$ (the
# leading correction to both isotherms is $O(P)$). The order comes from
# successive differences of the selectivity on a ladder of pressures a decade
# apart, each order of the ladder printed, with the finest difference measured
# against the solver's round-off floor; the deviation at the smallest pressure
# is therefore a first-order approach to the limit, not a round-off figure.
# Pure isotherms without a Henry region (Sips exponent $t = 1/2$,
# $n = m(bP)^t/(1+(bP)^t)$) are the break row: the order drops to $t$.

# %%
# (a) pure Henry isotherms: closed case of eqs. (30)-(35)
KH1, KH2, P_H = 3.0, 0.7, 0.37
yv = np.linspace(0.05, 0.95, 11)


def henry_exact(match=True):
    """Largest relative deviation of the IAST n1 from n1 = K1 P y1 (eq. 31) on pure
    Henry isotherms; match=False skips the spreading-pressure match (x_i = y_i)."""
    if match:
        n1 = iast(P_H, yv, henry_pure(KH1), henry_pure(KH2))["n1"]
    else:
        n1 = yv/(yv/(KH1*P_H) + (1 - yv)/(KH2*P_H))
    return float(np.max(np.abs(n1 - KH1*P_H*yv)/(KH1*P_H*yv)))


HENRY_EXACT = henry_exact()
claims.bound("henry_model_solver_max_rel", HENRY_EXACT, max=1e-12, basis="tolerance",
             note="newton step tolerance 1e-13 on ln z")
print(f"(a) solver on pure Henry isotherms vs n1 = K1 P y1 (eq. 31): max rel below 1e-12 "
      f"(here {HENRY_EXACT:.1e})")

# (b) pair U selectivity -> K1/K2 at first order
S_INF = (m1U*b1U)/(m2U*b2U)
HENRY_LADDER = [10**k for k in range(3, 9)]       # 1/P: P = 1e-3 ... 1e-8, ratio 10


def sel_mid(Pt, p1=pU1, p2=pU2):
    """IAST selectivity of pair U at y1 = 0.5."""
    return float(iast(Pt, np.array([0.5]), p1, p2)["s12"][0])


def henry_dev(Pt=1e-8, s_inf=S_INF):
    """|s/s_inf - 1| at pressure Pt."""
    return abs(sel_mid(Pt)/s_inf - 1.0)


def henry_order(p1=pU1, p2=pU2, ladder=HENRY_LADDER):
    """Observed orders in P of s(P) from successive differences, and the round-off margin
    of the finest difference (against the solver tolerance times |s|)."""
    st = observed_orders(lambda n: sel_mid(1.0/n, p1, p2), ladder)
    q_prev, q_last = st["values"][-2:]
    orders = np.array(st["orders"])
    return dict(orders=orders, worst=float(max(orders, key=lambda v: abs(v - 1.0))),
                margin=abs(q_last - q_prev)/(max(TOL, EPS)*abs(q_last)))


HENRY_DEV = henry_dev()
HO = henry_order()
HENRY_ORDER = HO["worst"]
claims.value("henry_selectivity_dev_P1em8", HENRY_DEV, rtol=1e-3,
             note="first-order approach to K1/K2 at P = 1e-8")
claims.bound("henry_convergence_order", HENRY_ORDER, min=0.95, max=1.05, basis="theory",
             note="the order farthest from 1 on the ladder P = 1e-3 ... 1e-8")
claims.bound("henry_order_roundoff_margin", HO["margin"], min=1e3, basis="tolerance",
             note="finest successive difference over tol*|s|")
print(f"\n(b) |s/ (K1/K2) - 1| at P = 1e-8: {HENRY_DEV:.3e} (first order in P: "
      f"{HENRY_DEV/1e-8:.3f} P)")
print("    observed orders in P on the ladder P = 1e-3 ... 1e-8: "
      + "  ".join(f"{o_:.6f}" for o_ in HO["orders"]))
print(f"    min {HO['orders'].min():.6f}, max {HO['orders'].max():.6f}; round-off margin of "
      f"the finest difference {HO['margin']:.1e}")

# %% [markdown]
# ### 7.4 The general spreading-pressure quadrature, refined and broken
#
# Eq. (19) for isotherms known only pointwise: log-grid trapezoid with the
# Henry segment $[0, t_0]$ added analytically as $K t_0$, the low-coverage
# care the paper's conclusions single out. Tested against the analytic
# $z = m\ln(1+bP)$ on pair U's component 1. The order comes from successive
# differences on a ladder of 64 to 512 intervals and must be 2; the error
# against the closed form is reported at 257 nodes. The two design choices
# (log grid, Henry anchor) each get a break row. Dropping the anchor is the
# instructive one: successive refinements still converge at second order, to a
# value that is wrong by the missing Henry segment, so refinement looks
# converged while being wrong.

# %%
ZA = pU1["z"](1.0)                              # eq. (19) in closed form at P = 1
QUAD_LADDER = (64, 128, 256, 512)               # intervals; nodes = intervals + 1


def quad_study(ladder=QUAD_LADDER, **kw):
    """Error of z_quad at 257 nodes against the closed form, the observed orders from
    successive differences on the ladder, and the round-off margin of the finest
    difference (against n eps |z|, the summation round-off)."""
    st = observed_orders(lambda n: z_quad(1.0, pU1, n + 1, **kw), ladder)
    q_prev, q_last = st["values"][-2:]
    orders = np.array(st["orders"])
    return dict(err=abs(z_quad(1.0, pU1, 257, **kw)/ZA - 1.0), orders=orders,
                worst=float(max(orders, key=lambda v: abs(v - 2.0))),
                margin=abs(q_last - q_prev)/(max(ladder)*EPS*abs(q_last)))


QS = quad_study()
QUAD_ERR, QUAD_ORDER = QS["err"], QS["worst"]
claims.value("pi_quadrature_rel_err_257", QUAD_ERR, rtol=1e-3)
claims.bound("pi_quadrature_observed_order", QUAD_ORDER, min=1.9, max=2.1, basis="theory",
             note="trapezoid rule, ladder of 64 to 512 intervals")
claims.bound("pi_quadrature_roundoff_margin", QS["margin"], min=1e3, basis="tolerance",
             note="finest successive difference over n eps |z|")
print(f"rel err at 257 nodes {QUAD_ERR:.3e}; observed orders on {QUAD_LADDER} intervals: "
      + "  ".join(f"{o_:.4f}" for o_ in QS["orders"])
      + f"; round-off margin {QS['margin']:.1e}")

# the Henry anchor dropped, the grid starting at t0 = 1e-3
NA = quad_study(t0_frac=1e-3, anchor=False)
claims.value("pi_quadrature_no_anchor_rel", NA["err"], rtol=1e-3)
claims.bound("pi_quadrature_no_anchor_order", NA["worst"], min=1.9, max=2.1, basis="theory",
             note="refinement converges, at second order, to the wrong value")
print(f"Henry anchor dropped (t0 = 1e-3): rel err {NA['err']:.3e}, observed orders "
      + "  ".join(f"{o_:.4f}" for o_ in NA["orders"])
      + " (converged, and wrong by the missing segment)")
print(f"restoring the anchor at the same t0: rel err "
      f"{abs(z_quad(1.0, pU1, 257, t0_frac=1e-3)/ZA - 1.0):.3e}")

# %% [markdown]
# ### 7.5 The Gibbs adsorption isotherm along a path, eq. (21) as a solver test
#
# At constant $T$ and $P$, eq. (21) says
# $\;(A/RT)\,d\pi = n_1\,d\ln y_1 + n_2\,d\ln y_2$. The left side comes from
# the solver's $z(y_1)$, the *pure-component* route, eq. (19) at the matched
# $P_i^{\circ}$, while the right side integrates the *mixture* amounts along
# the composition path. The two routes share the solve but not the
# thermodynamics: if the solver returned self-consistent but wrong states
# (mismatched spreading pressures, or a wrong $n_t$ rule), the identity fails.
# Measured on 4001-node trapezoids, where the quadrature error of the healthy
# model is negligible: the **arithmetic-mean** $n_t$ (the plausible wrong
# reading of eq. (26)) misses by 1.5 %, and a 1 % mismatch between the two
# components' spreading pressures by 0.16 %.

# %%
YA, YB = 0.05, 0.95
PATH_LADDER = (64, 128, 256, 512)               # intervals; nodes = intervals + 1
N_FINE = 4001                                   # nodes of the defect measurements


def path_dz(**iast_kw):
    """Pure-route Delta z between the path ends, from the solver's z (defects passed on)."""
    s_ = iast(1.0, np.array([YA, YB]), pU1, pU2, **iast_kw)["z"]
    return float(s_[1] - s_[0])


DZ = path_dz()


def path_study(ladder=PATH_LADDER, rule="trapezoid", **iast_kw):
    """Residual of eq. (21) on 65 nodes and with Gauss-64, the observed orders of the
    trapezoid from successive differences, and the round-off margin of the finest
    difference (against max(tol, n eps) |q|)."""
    dz = path_dz(**iast_kw)
    st = observed_orders(lambda n: path_rhs_trap(n + 1, 1.0, YA, YB, pU1, pU2, rule=rule,
                                                 **iast_kw), ladder)
    q_prev, q_last = st["values"][-2:]
    orders = np.array(st["orders"])
    return dict(
        trap65=abs(path_rhs_trap(65, 1.0, YA, YB, pU1, pU2, rule=rule, **iast_kw) - dz)/abs(dz),
        gauss64=abs(path_rhs_gauss(64, 1.0, YA, YB, pU1, pU2, **iast_kw) - dz)/abs(dz),
        orders=orders, worst=float(max(orders, key=lambda v: abs(v - 2.0))),
        margin=abs(q_last - q_prev)/(max(TOL, max(ladder)*EPS)*abs(q_last)))


def path_fine(**iast_kw):
    """Residual of eq. (21) on N_FINE nodes: the effect of a defect, quadrature removed."""
    return abs(path_rhs_trap(N_FINE, 1.0, YA, YB, pU1, pU2, **iast_kw)
               - path_dz(**iast_kw))/abs(path_dz(**iast_kw))


PS = path_study()
PATH_TRAP65, PATH_GAUSS, PATH_TRAP_ORDER = PS["trap65"], PS["gauss64"], PS["worst"]
claims.value("gibbs_path_trap65_rel", PATH_TRAP65, rtol=1e-3)
claims.bound("gibbs_path_gauss64_rel", PATH_GAUSS, max=1e-12, basis="tolerance",
             note="Gauss-64 is converged; the gap is the newton tolerance")
claims.bound("gibbs_path_trap_order", PATH_TRAP_ORDER, min=1.9, max=2.1, basis="theory",
             note="trapezoid rule, ladder of 64 to 512 intervals")
claims.bound("gibbs_path_roundoff_margin", PS["margin"], min=1e3, basis="tolerance",
             note="finest successive difference over max(tol, n eps) |q|")
print(f"pure-route Delta z = {DZ:.10f}   (P = 1, y1: {YA} -> {YB})")
print(f"mixture-route integral: Gauss-64 residual below 1e-12 (here {PATH_GAUSS:.1e});")
print(f"  trapezoid-65 residual {PATH_TRAP65:.3e}; observed orders on {PATH_LADDER} intervals: "
      + "  ".join(f"{o_:.4f}" for o_ in PS["orders"]) + f"; round-off margin {PS['margin']:.1e}")

PATH_ARITH = path_fine(nt_mode="arith")
PATH_MISMATCH = path_fine(zfac2=1.01)
claims.value("gibbs_path_arith_mean_rel", PATH_ARITH, rtol=1e-3)
claims.value("gibbs_path_mismatch_rel", PATH_MISMATCH, rtol=1e-3)
print(f"\non {N_FINE}-node trapezoids (healthy model: {path_fine():.1e}):")
print(f"  n_t read as the ARITHMETIC mean of eq. (26)     -> residual {PATH_ARITH:.3e}")
print(f"  spreading pressures mismatched by 1 % (eq. 40)  -> residual {PATH_MISMATCH:.3e}")

# %% [markdown]
# ### 7.6 The closed loop: the paper's consistency charge, quantified
#
# The paper reports (p. 125) that Bering and Serpenskii "noted that the
# Langmuir equation for mixture adsorption was not thermodynamically
# consistent". Here is that sentence as a number. If a mixture model
# $n_i(P, y_1)$ respects the Gibbs adsorption isotherm, the one-form
# $n_t\,d\ln P + (n_1/y_1 - n_2/y_2)\,dy_1$ is exact and its integral around
# any closed loop vanishes: $\pi$ is a state function. Around the rectangle
# $P \in [0.2, 1]$, $y_1 \in [0.2, 0.8]$:
#
# - **IAST** closes the loop to round-off (below 1e-12), as its construction
#   from a single $\pi$ guarantees *in exact arithmetic*; the test can fail in
#   code: the same loop with the two components' spreading pressures
#   mismatched by 1 % (a solver bug this test would catch) fails at 2e-3.
# - **Extended Langmuir with unequal capacities does not close**: the deficit
#   is 12 % of the loop's own spreading-pressure scale. Computed **twice, by
#   routes sharing no assembly**, the boundary line integral of the model's
#   $n_i$, and the area integral of the curl that §2.1 derived symbolically
#   ($(m_2{-}m_1) b_1 b_2 P^2/D^2$), agreeing to round-off (below 1e-12). The curl's closed
#   form makes the deficit **exactly linear in $m_2 - m_1$** (its denominator
#   contains no $m$), and doubling $m_1 - m_2$ doubles the measured loop
#   deficit to round-off (below 1e-12), a third route to the same structure. At
#   $m_2 = m_1$ everything vanishes identically: the consistent case is
#   exactly the case where extended Langmuir *is* an IAST (§2.1, §7.1).

# %%
LOOP = dict(Pa=0.2, Pb=1.0, ya=0.2, yb=0.8)
Z_SCALE = float(pU1["z"](LOOP["Pb"]) - pU1["z"](LOOP["Pa"]))    # the loop's own z scale


def model_iast(**iast_kw):
    """IAST on pair U as a mixture model (P, y1) -> (n1, n2); defects passed on."""
    def f(Pt, yv_):
        s_ = iast(Pt, yv_, pU1, pU2, **iast_kw)
        return s_["n1"], s_["n2"]
    return f


def model_extl(m2=m2U):
    """Extended Langmuir on pair U's constants, capacity m2 as an argument."""
    return lambda Pt, yv_: ext_langmuir(Pt, yv_, m1U, b1U, m2, b2U)


def circulation_green(m2=m2U, ngauss=48):
    """Route 2: Green's theorem, the area integral of the symbolically derived curl."""
    xg, wg = leggauss(ngauss)
    hl = 0.5*(np.log(LOOP["Pb"]) - np.log(LOOP["Pa"]))
    lnv = hl*xg + 0.5*(np.log(LOOP["Pa"]) + np.log(LOOP["Pb"]))
    yv_ = 0.5*(LOOP["yb"] - LOOP["ya"])*xg + 0.5*(LOOP["ya"] + LOOP["yb"])
    LN, YV = np.meshgrid(lnv, yv_)
    W2 = np.outer(wg, wg)*hl*0.5*(LOOP["yb"] - LOOP["ya"])
    return float(np.sum(W2*curl_fn(LN, YV, m1U, b1U, m2, b2U)))


def loop_study(m2=m2U, ngauss=48, skip_edge4=False):
    """Extended-Langmuir loop deficit, its share of the z scale, the two-route gap and
    the deviation from exact linearity in m2 - m1 (m2 -> 0 doubles m1 - m2 here)."""
    ce = circulation(model_extl(m2), ngauss=ngauss, **LOOP)
    route = abs(ce/circulation_green(m2) - 1.0)
    lin = abs(circulation(model_extl(m1U - 2*(m1U - m2)), ngauss=ngauss, skip_edge4=skip_edge4,
                          **LOOP)
              / circulation(model_extl(m2), ngauss=ngauss, skip_edge4=skip_edge4, **LOOP) - 2.0)
    return dict(circ=ce, rel=ce/Z_SCALE, route=route, lin=lin)


CIRC_IAST = abs(circulation(model_iast(), **LOOP))
CIRC_IAST_MISMATCH = abs(circulation(model_iast(zfac2=1.01), **LOOP))
LS = loop_study()
CIRC_EXTL, CIRC_REL, ROUTE_AGREE, LINEARITY_DEV = LS["circ"], LS["rel"], LS["route"], LS["lin"]
claims.bound("loop_circulation_iast_abs", CIRC_IAST, max=1e-12, basis="identity",
             note="Gauss-48 on every edge is converged; zero for any model with a spreading pressure")
claims.value("loop_circulation_iast_mismatch_abs", CIRC_IAST_MISMATCH, rtol=1e-3,
             note="the same loop with a 1 % spreading-pressure mismatch injected")
claims.value("loop_circulation_ext_langmuir", CIRC_EXTL, rtol=1e-6)
claims.value("loop_circulation_over_scale", CIRC_REL, rtol=1e-6)
claims.value("loop_deficit_share_of_scale", abs(CIRC_REL), rtol=1e-6,
             note="the magnitude of the deficit as a fraction of the loop's z scale")
claims.bound("circulation_route_rel_agreement", ROUTE_AGREE, max=1e-12, basis="identity",
             note="Green's theorem; both routes on converged Gauss-48")
claims.bound("circulation_linearity_dev", LINEARITY_DEV, max=1e-12, basis="identity",
             note="the curl is linear in m2 - m1, its denominator holds no m")

print(f"loop deficit, IAST: below 1e-12 (here {CIRC_IAST:.1e}); with a 1 % spreading-pressure "
      f"mismatch injected: {CIRC_IAST_MISMATCH:.2e}")
print(f"loop deficit, extended Langmuir:  {CIRC_EXTL:.6f}  =  {100*CIRC_REL:.2f} % of the "
      f"loop's z-scale {Z_SCALE:.6f}")
print(f"same deficit via the symbolic curl (Green route): {circulation_green():.6f}   "
      f"routes agree to below 1e-12 (here {ROUTE_AGREE:.1e})")
print(f"doubling m1 - m2 multiplies the deficit by 2 to below 1e-12 (here {LINEARITY_DEV:.1e} off)")

# %% [markdown]
# ### 7.7 The lever rules, numerically (their symbolic proof is §2.1)
#
# An identity in exact arithmetic, the numeric check is kept because it is
# what §6.2's drawing rests on, and because its *break* (chords swapped, the
# plausible mis-drawing of Fig. 1) moves the answer by order one.

# %%
print(f"lever rules vs eqs. (22)/(29): max deviation below 1e-12 (here {LEVER_DEV:.1e})")
print(f"break, chords OB and OC swapped: {LEVER_SWAPPED:.4f}")

# %% [markdown]
# ### 7.8 The shape-translation result: when IAST selectivity cannot vary
#
# **Claim (this page's, not the paper's):** if the two pure isotherms are the
# same curve shifted in $\log P$, $n_2^{\circ}(P) = n_1^{\circ}(\kappa P)$:
# then eq. (19) gives $\pi_2(P^{\circ}) = \pi_1(\kappa P^{\circ})$, the match
# (40) forces $P_1^{\circ} = \kappa P_2^{\circ}$, and the selectivity (28) is
# **exactly** $1/\kappa$, independent of $P$ and composition. The Henry-limit
# constancy of eqs. (30)–(35) is the linear special case (any two lines are
# shape-identical, $\kappa = K_2/K_1$); §6.3's pair L is the equal-capacity
# Langmuir case ($\kappa = b_2/b_1$). Searched for in the paper and not found
# there: the Fig. 10 discussion (pp. 125–126) treats the composition dependence
# as a prediction to compare with data, and the conclusions do not connect
# selectivity constancy to isotherm shape; the closest printed statement is the
# Henry-limit special case itself. Verified here on pair S, dual-site
# isotherms the closed forms of §7.1 never see, inverted numerically, and
# broken by perturbing one site's affinity so the shapes genuinely differ.

# %%
Pg2, yg2 = np.meshgrid(np.geomspace(0.05, 5, 12), np.linspace(0.1, 0.9, 12))


def shape_spread(b_site2=0.5):
    """Relative spread of the IAST selectivity of pair S over the grid; b_site2 is the
    second-site affinity of component 2 (0.5 makes it the same curve as component 1)."""
    p2 = dualsite_pure(0.7, 20.0, 0.3, b_site2, KAPPA)
    return rel_spread(iast(Pg2.ravel(), yg2.ravel(), pS1, p2)["s12"])


sS = iast(Pg2.ravel(), yg2.ravel(), pS1, pS2)
SHAPE_SPREAD = shape_spread()
SHAPE_PERT_SPREAD = shape_spread(0.5*1.05)
claims.bound("shape_translation_selectivity_rel_spread", SHAPE_SPREAD, max=1e-10, basis="tolerance",
             note="exactly constant by the shape identity; nested newton tolerance 1e-13")
claims.value("shape_perturbed_selectivity_rel_spread", SHAPE_PERT_SPREAD, rtol=1e-4)
print(f"pair S (same shape, kappa = {KAPPA}): selectivity mean {np.mean(sS['s12']):.12f} "
      f"vs 1/kappa = {1/KAPPA}")
print(f"  relative spread over {Pg2.size} states spanning 100x in P: below 1e-10 "
      f"(here {SHAPE_SPREAD:.1e})")
print(f"one site's b perturbed 5 %: spread {SHAPE_PERT_SPREAD:.3e} "
      f"(the constancy is the shape identity, not the solver)")


# dilution-limit second route for 6.3's closed forms: the solver at small y1
def dilute1_route(y1=1e-8):
    """Relative gap of the solver's selectivity at y1 from the closed-form y1 -> 0 limit."""
    return abs(float(iast(1.0, np.array([y1]), pU1, pU2)["s12"][0])/S_DIL1 - 1.0)


DIL1_ROUTE_DEV, DIL1_ROUTE_DEV_7 = dilute1_route(1e-8), dilute1_route(1e-7)
claims.value("selectivity_dilute1_solver_vs_closed_rel", DIL1_ROUTE_DEV, rtol=1e-3,
             note="first-order approach in y1 at y1 = 1e-8")
claims.value("selectivity_dilute1_solver_vs_closed_rel_y1em7", DIL1_ROUTE_DEV_7, rtol=1e-3,
             note="the same at y1 = 1e-7: ten times larger, linear in y1")
print(f"\n6.3's dilute-in-1 closed form checked by the solver: rel dev {DIL1_ROUTE_DEV:.2e} at "
      f"y1 = 1e-8, {DIL1_ROUTE_DEV_7:.2e} at y1 = 1e-7 (linear in y1)")
S_DIL_L = dilute_selectivities(1.0, pL1, pL2)
print(f"composition-dependence factor, pair U: {S_DIL_RATIO:.6f};  pair L (equal shape): "
      f"{S_DIL_L[1]/S_DIL_L[0]:.6f} (collapses to 1)")

# %% [markdown]
# ### 7.9 The defect-injection table
#
# Every claim the page registers has a row here that fails its own check. Each
# row injects one defect through an argument of the study that computes the
# claim, everywhere the parameter enters, and recomputes the claim; a study
# that raises records nan for the claims it targets and stays in the table.
# The last row removes every injected defect from the four sensitivity studies
# (and the shape perturbation), so their values collapse: it is the row that
# shows those numbers are the defects' effects and nothing else.

# %%
def with_printed(key, v):
    """The transcribed scalars with one entry replaced."""
    c = dict(C)
    c[key] = v
    return c


pU2_same_m = langmuir_pure(m1U, b2U)              # pair U with m2 -> m1
ROWS = [
    ("eq. (11) read as printed (no operator: a product)", {
        "sym_zero_identity_count": lambda: sum(r == 0 for _, r in
                                               [("(11)", eq11_to_13("product"))] + ZERO[1:])}),
    ("adsorbent counted as a component in eq. (18)", {
        "eq18_dof_binary": lambda: phase_rule_dof(components=3)}),
    ("ref. 5 year 1959 read as 1969", {
        "ref5_volumes_per_year": lambda: ref5_volume_rate(with_printed("ref5_year", 1969.0))}),
    ("m2 = 0.6 in pair L (capacities made unequal)", {
        "equal_m_closed_form_max_rel": lambda: eqm_dev(pL1, langmuir_pure(0.6, b2L), LANG_L, GP, GY),
        "pair_L_selectivity_mean": lambda: pair_L_selectivity(langmuir_pure(0.6, b2L))["mean"],
        "pair_L_selectivity_rel_spread": lambda: pair_L_selectivity(langmuir_pure(0.6, b2L))["spread"]}),
    ("b1 wrong by 1 % in pair L", {
        "equal_m_closed_form_max_rel": lambda: eqm_dev(langmuir_pure(mL, 1.01*b1L), pL2, LANG_L,
                                                       GP, GY)}),
    ("newton stopped on a 1e-3 step alone (no backward-error test)", {
        "equal_m_closed_form_max_rel": lambda: eqm_dev(pL1, pL2, LANG_L, GP, GY, tol=1e-3,
                                                       backward_tol=None)}),
    ("m2 -> m1 in pair U (the models coincide again)", {
        # the scan maximum without the end guard: at round-off the argmax position is noise
        "unequal_m_nt_dev_windowsup": lambda: float(np.max(nt_dev(10.0, np.linspace(Y_EDGE, 1 - Y_EDGE, N_SCAN),
                                                                  pU1, pU2_same_m, (m1U, b1U, m1U, b2U)))),
        "unequal_m_n1_dilute_limit_dev": lambda: n1_dilute_dev(10.0, pU1, pU2_same_m, m1U, b1U, b2U),
        "loop_circulation_ext_langmuir": lambda: circulation(model_extl(m1U), **LOOP),
        "loop_circulation_over_scale": lambda: circulation(model_extl(m1U), **LOOP)/Z_SCALE,
        "loop_deficit_share_of_scale": lambda: abs(circulation(model_extl(m1U), **LOOP))/Z_SCALE}),
    ("window supremum read off the y1 scan, not refined", {
        "unequal_m_nt_dev_windowsup": lambda: window_sup(refine=False)["sup"],
        "unequal_m_nt_dev_y_at_sup": lambda: window_sup(refine=False)["y"],
        "unequal_m_nt_dev_grid_low_frac": lambda: window_sup(refine=False)["grid_low"]}),
    ("spreading-pressure matching skipped (x_i = y_i)", {
        "henry_model_solver_max_rel": lambda: henry_exact(match=False)}),
    ("Henry selectivity read at P = 0.1 instead of the limit", {
        "henry_selectivity_dev_P1em8": lambda: henry_dev(Pt=0.1)}),
    ("Henry constant taken as b (capacity forgotten)", {
        "henry_selectivity_dev_P1em8": lambda: henry_dev(s_inf=b1U/b2U)}),
    ("pure isotherms with Sips exponent t = 0.5 (no Henry region)", {
        "henry_convergence_order": lambda: henry_order(langmuir_pure(m1U, b1U, 0.5),
                                                       langmuir_pure(m2U, b2U, 0.5))["worst"]}),
    ("Henry ladder continued to P = 1e-13 (into the solver tolerance)", {
        "henry_order_roundoff_margin": lambda: henry_order(
            ladder=[10**k for k in range(8, 14)])["margin"]}),
    ("linear grid instead of logarithmic in eq. (19)", {
        "pi_quadrature_rel_err_257": lambda: quad_study(grid="linear")["err"]}),
    ("Henry anchor dropped (t0 = 1e-3)", {
        "pi_quadrature_rel_err_257": lambda: quad_study(t0_frac=1e-3, anchor=False)["err"]}),
    ("left rectangle rule instead of the trapezoid", {
        "pi_quadrature_observed_order": lambda: quad_study(rule="rectangle")["worst"],
        "pi_quadrature_no_anchor_order": lambda: quad_study(t0_frac=1e-3, anchor=False,
                                                            rule="rectangle")["worst"],
        "gibbs_path_trap_order": lambda: path_study(rule="rectangle")["worst"]}),
    ("quadrature ladder continued to 2^20 intervals (into round-off)", {
        "pi_quadrature_roundoff_margin": lambda: quad_study(
            ladder=(2**17, 2**18, 2**19, 2**20))["margin"]}),
    ("n_t read as the arithmetic mean of eq. (26)", {
        "gibbs_path_trap65_rel": lambda: path_study(ladder=(64, 128, 256), nt_mode="arith")["trap65"],
        "gibbs_path_gauss64_rel": lambda: path_study(ladder=(64, 128, 256), nt_mode="arith")["gauss64"]}),
    ("spreading pressures mismatched by 1 % (a broken eq. 40)", {
        "gibbs_path_trap65_rel": lambda: path_study(ladder=(64, 128, 256), zfac2=1.01)["trap65"],
        "loop_circulation_iast_abs": lambda: abs(circulation(model_iast(zfac2=1.01), **LOOP))}),
    ("path ladder continued to 2^16 intervals (into round-off)", {
        "gibbs_path_roundoff_margin": lambda: path_study(ladder=(2**14, 2**15, 2**16))["margin"]}),
    ("m2 -> 0 in extended Langmuir (m1 - m2 doubled)", {
        "loop_circulation_ext_langmuir": lambda: circulation(model_extl(0.0), **LOOP),
        "loop_circulation_over_scale": lambda: circulation(model_extl(0.0), **LOOP)/Z_SCALE,
        "loop_deficit_share_of_scale": lambda: abs(circulation(model_extl(0.0), **LOOP))/Z_SCALE}),
    ("line-integral route on 2-point Gauss", {
        "circulation_route_rel_agreement": lambda: loop_study(ngauss=2)["route"]}),
    ("edge 4 of the loop dropped", {
        "circulation_linearity_dev": lambda: loop_study(skip_edge4=True)["lin"]}),
    ("chords OB and OC swapped in Fig. 1's construction", {
        "lever_rule_max_dev": lambda: lever_rule(pU1, pU2, ZSTAR, X1STAR, swap=True)["dev"]}),
    ("one site's affinity of component 2 perturbed 5 %", {
        "shape_translation_selectivity_rel_spread": lambda: shape_spread(0.5*1.05)}),
    ("dilute selectivity read at y1 = 0.1 instead of the limit", {
        "selectivity_dilute1_solver_vs_closed_rel": lambda: dilute1_route(0.1),
        "selectivity_dilute1_solver_vs_closed_rel_y1em7": lambda: dilute1_route(0.1)}),
    ("pair L (equal shape) in place of pair U", {
        "selectivity_dilute1_P1": lambda: dilute_selectivities(1.0, pL1, pL2)[0],
        "selectivity_dilute2_P1": lambda: dilute_selectivities(1.0, pL1, pL2)[1],
        "selectivity_dilute_ratio_P1": lambda: (lambda d: d[1]/d[0])(dilute_selectivities(1.0, pL1, pL2))}),
    ("injected defects removed from the sensitivity studies", {
        "gibbs_path_arith_mean_rel": lambda: path_fine(nt_mode="harmonic"),
        "gibbs_path_mismatch_rel": lambda: path_fine(zfac2=1.0),
        "loop_circulation_iast_mismatch_abs": lambda: abs(circulation(model_iast(zfac2=1.0), **LOOP)),
        "pi_quadrature_no_anchor_rel": lambda: quad_study(t0_frac=1e-3, anchor=True)["err"],
        "shape_perturbed_selectivity_rel_spread": lambda: shape_spread(0.5)}),
]

claims.define_rows([r[0] for r in ROWS])
BREAK = []
for row, studies in ROWS:
    for name, fn in studies.items():
        try:
            out = fn()
        except Exception as exc:                 # noqa: BLE001 - a failed study is a result
            print(f"  {row}: {name} raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
            out = np.nan
        caught_ = claims.breaks(row, name, out)
        base_ = claims.metrics[name]["value"]
        if claims.metrics[name]["kind"] == "bound":   # round-off digits differ between machines
            base_, out = f"(here {base_:.1e})", f"(here {float(out):.1e})"
        BREAK.append((row, name, base_, out, caught_))

bt = pd.DataFrame(BREAK, columns=["injected defect", "metric", "baseline", "with defect", "caught"])
pd.set_option("display.max_colwidth", 60)
pd.set_option("display.max_rows", 300)
print(bt.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
idle = [r for r in bt["injected defect"].unique()
        if not bt[(bt["injected defect"] == r) & bt.caught].shape[0]]
print(f"\n{len(ROWS)} rows defined, {bt['injected defect'].nunique()} in the table; "
      f"rows that fail no claim: {idle or 'none'}")
missing = [k for k in claims.metrics if not bt[(bt.metric == k) & bt.caught].shape[0]]
print(f"claims no row fails: {missing or 'none'}")

# %% [markdown]
# ## What pymrm adds
#
# Nothing to the theory, and the page says so plainly: IAST is algebra plus one
# scalar root-find per state, and most of this notebook would run with pymrm
# uninstalled, exactly as on `J1.1`, `J1.3` and `A1.6`. What `newton` and
# `NumJac` do here is narrow and real:
#
# 1. **The whole $(P, y_1)$ grid as one Newton problem.** The spreading-pressure
#    match is solved for hundreds of states simultaneously with
#    `NumJac((K, 1))`, independent scalars, a diagonal Jacobian by
#    construction. The shape is `(K, 1)` and never `(K,)`: the bare 1-D shape
#    would couple every state to every other and build a dense $K \times K$
#    Jacobian for no change in the answer. The iteration runs in $\ln z$,
#    which keeps Newton on the feasible side of $P^{\circ}(0) = 0$, and every
#    solve raises unless `newton` reports success.
# 2. **The general-isotherm inverse.** For pair S the pure isotherms have no
#    analytic inverse; $P^{\circ}(z)$ is itself a vectorised `newton` solve
#    nested inside the IAST solve. That is the machinery a user with tabulated
#    isotherms (eq. 43) actually needs.
# 3. **Refined maxima, closed forms and measured orders where a sweep would
#    have been wrong.** The §7.2 supremum is a bounded maximisation refined
#    from a scan (the grid read is 0.7 % low, printed), the dilution limits are
#    closed forms rather than end-of-grid reads, and the quadrature and path
#    integrals carry observed orders from successive differences
#    (`pymrm.checks.observed_orders`), each with its round-off margin.
#
# The one thing on this page that is genuinely useful beyond reproduction is
# not pymrm at all: the **closed-loop consistency test of §7.6** (with its
# symbolically derived curl) as a portable diagnostic. It applies to *any*
# proposed mixture-adsorption model $n_i(P, y)$, needs no IAST solve, and
# returns a single number that is zero iff the model carries a consistent
# spreading pressure. The extended-Langmuir deficit, 12 % of scale on an
# ordinary parameter set, is what it reports on the most widely used shortcut
# model in the field, and the $m_2 = m_1$ collapse identifies *why* the
# shortcut is safe exactly when capacities match.

# %% [markdown]
# ## Reuse
#
# **Use this page for** a working, tested IAST solver (`iast` +
# `langmuir_pure`/`dualsite_pure`/`henry_pure`), the consistency-loop
# diagnostic (`circulation`), and the reference solutions in
# `data/iast-illustrative-reference.csv` for regression-testing your own
# implementation (parameters in the file header; equal-capacity pair L also has
# the closed form to check against).
#
# **Four things to carry away:**
#
# 1. **Extended Langmuir is safe exactly when capacities match.** Then it *is*
#    the ideal adsorbed solution (§2.1, symbolically). With unequal capacities
#    it is thermodynamically inconsistent: a 12 %-of-scale loop deficit on an
#    unremarkable parameter set (§7.6) and a 21.5 % total-amount error against
#    IAST at the illustrative window's edge (§7.2); the dilute component
#    is where it is worst (64 % there). If you must use it, check the capacity
#    ratio first.
# 2. **Eq. (26) is a harmonic mix.** The arithmetic reading is the natural slip
#    and it *fails the Gibbs isotherm* by 1.5 % on this page's mild example
#    (§7.5), a thermodynamic error, not a numerical one: refinement converges
#    to the wrong value, so no grid refinement will ever surface it.
# 3. **The low-pressure tail of the pure isotherms is load-bearing.** The
#    paper's own conclusions say so; §7.4 quantifies it: dropping the Henry
#    anchor at $t_0 = 10^{-3}$ costs 4e-3 in $z$ while successive refinements
#    still converge at second order, so refinement looks converged while being
#    wrong. Measure the Henry
#    constant, or extrapolate to it explicitly, before integrating eq. (19).
# 4. **Composition-dependent selectivity is a shape signal.** Under IAST,
#    isotherms differing only by a $\log P$ shift give exactly constant
#    $s_{1,2}$ (§7.8). If your IAST calculation shows strong composition
#    dependence, your two pure isotherms differ in shape; and if your *data*
#    show it where the fitted pure isotherms do not, IAST with those fits
#    cannot represent the system, whatever the parameters.
#
# **Do not use this page as evidence that IAST works.** It contains no
# experimental comparison (the paper's four validation systems are figures of
# external data, out of scope here), and every agreement on it is structural.
# The empirical literature on IAST is enormous; this page is about what the
# 1965 paper actually established on its own printed content. **Do not use it
# for dynamics**: breakthrough is `J1.5`, and this page's `iast` can serve as
# that page's mixture-equilibrium closure. **Do not read pairs L, U, S as any
# physical system**; they are illustrative parameter sets and nothing more.

# %% [markdown]
# ### Claims and the prose
#
# Every number the prose of this notebook and of the README states is checked
# below against the claim it quotes, so a number in the text cannot drift from
# the computation and a cut in its precision shows in CI (agreement.json
# "quoted"). Round-off-level results are quoted as bounds.

# %%
SAYS = {
    "sym_zero_identity_count": ["11"], "eq18_dof_binary": ["3"], "ref5_volumes_per_year": ["4.6"],
    "equal_m_closed_form_max_rel": ["below 1e-12"], "henry_model_solver_max_rel": ["below 1e-12"],
    "loop_circulation_iast_abs": ["below 1e-12"], "circulation_route_rel_agreement": ["below 1e-12"],
    "circulation_linearity_dev": ["below 1e-12"], "gibbs_path_gauss64_rel": ["below 1e-12"],
    "shape_translation_selectivity_rel_spread": ["below 1e-10"],
    "loop_deficit_share_of_scale": ["12 %", "11.9 %"], "loop_circulation_ext_langmuir": ["-0.1548"],
    "loop_circulation_iast_mismatch_abs": ["2e-3"],
    "unequal_m_nt_dev_windowsup": ["21.5 %", "0.2151"], "unequal_m_nt_dev_y_at_sup": ["0.0996"],
    "unequal_m_nt_dev_grid_low_frac": ["0.7 %", "0.67 %"],
    "unequal_m_n1_dilute_limit_dev": ["64 %", "64.2 %", "0.6417"],
    "selectivity_dilute1_P1": ["13.660254"], "selectivity_dilute2_P1": ["60"],
    "selectivity_dilute_ratio_P1": ["4.39", "4.3923"], "pair_L_selectivity_mean": ["10"],
    "gibbs_path_arith_mean_rel": ["1.5 %"], "gibbs_path_mismatch_rel": ["0.16 %"],
    "gibbs_path_trap65_rel": ["8.8e-5"], "pi_quadrature_no_anchor_rel": ["4e-3"],
    "pi_quadrature_rel_err_257": ["1.49e-5"],
    "henry_selectivity_dev_P1em8": ["2.75e-8"], "henry_convergence_order": ["between 0.95 and 1.05", "at least 0.95",
                                                                      "at most 1.05"],
    "pi_quadrature_observed_order": ["between 1.9 and 2.1"], "gibbs_path_trap_order": ["between 1.9 and 2.1"],
    "pi_quadrature_no_anchor_order": ["between 1.9 and 2.1", "at least 1.9", "at most 2.1"],
    "selectivity_dilute1_solver_vs_closed_rel": ["2.7e-8"],
    "selectivity_dilute1_solver_vs_closed_rel_y1em7": ["2.7e-7"],
    "shape_perturbed_selectivity_rel_spread": ["1.1e-2"],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)
print(f"prose literals: {sum(len(v) for v in SAYS.values())} checked against their claims with claims.says")
claims.report();   # the semicolon keeps the notebook from echoing the payload

# %% [markdown]
# ## References
#
# Myers, A. L. and Prausnitz, J. M. (1965). Thermodynamics of mixed-gas
# adsorption. *AIChE Journal* **11**(1), 121–127.
# [doi:10.1002/aic.690110125](https://doi.org/10.1002/aic.690110125), **the
# paper, and the only document read for content.** Identity confirmed from its
# own first page on a native-resolution render: the title, the by-line "A. L.
# MYERS and J. M. PRAUSNITZ / University of California, Berkeley, California",
# the footnote "A. L. Myers is at the University of Pennsylvania, Philadelphia,
# Pennsylvania", the running feet "Vol. 11, No. 1 / A.I.Ch.E. Journal /
# Page 121", and the Wiley PDF margin naming the same DOI; the Wiley Subject
# metadata independently gives "AIChE Journal 1965.11:121-127". Manuscript
# received April 24, 1964; revision received August 14, 1964; accepted
# August 19, 1964. `pdfimages -list` reports all seven pages as CCITT-G4
# bilevel at **300 ppi native**, so every numeric was read from a crop of a
# 300 ppi render at digit scale; the text layer was used only as a search
# index. (The seventh PDF page carries the start of the *next* article, Gabor,
# "Lateral Transport in a Fluidized-Packed Bed", below Myers & Prausnitz's
# reference list; nothing was read from it.)
#
# **Cited by the paper, not consulted, and no number here derives from them**:
# these four are the sources of every experimental point in Figs. 3–10, which
# is exactly why the experimental case is out of this page's scope:
#
# Szepesy, L. and Illés, V. (1963). *Acta Chim. Hung.* **35**, 37; *ibid.*
# p. 53; *ibid.* p. 245 (refs. 9–11), the methane, ethane and mixture
# isotherms behind Figs. 2–4. Markham, E. C. and Benton, A. F. (1931). *J. Am.
# Chem. Soc.* **53**, 497 (ref. 8), the CO–O₂/silica data of Figs. 5–6 (the
# Fig. 6 legend's "(1950)" contradicts this, the paper's only Markham & Benton
# reference; see §4). Lewis, W. K., Gilliland, E. R., Chertow, B. and
# Hoffman, W. H. (1950). *J. Am. Chem. Soc.* **72**, 1153 (ref. 7), the
# propane–propylene/silica-gel data of Fig. 7. Bering, B. P. and
# Serpenskii, V. V. (1952). *Zhur. Fiz. Khim.* **26**, 253 (ref. 2), the
# ethylene–CO₂ equilibrium surface behind Figs. 8–10 (the reference and text
# spell "Serpenskii"; the three figure legends spell "SERPINSKII", one
# surname, two transliterations, reported as an observation).
#
# Also cited, not consulted: Arnold (1949), *J. Am. Chem. Soc.* **71**, 104:
# the liquid entropy model the paper supersedes; de Boer (1953); Hill (1949),
# *J. Chem. Phys.* **17**, 520, and Hill (1960), *Introduction to Statistical
# Mechanics*, the adsorption thermodynamics and the low-coverage statistical
# mechanics behind eqs. (4), (30), (31); Young & Crowell (1962), the
# mixed-gas-adsorption review. **Ref. 5** ("Ibid, 63, 456 (1959)") is cited
# nowhere in the paper's running text and its "Ibid." (= *J. Chem. Phys.*)
# cannot carry volume 63 in 1959 (ref. 4's own pair puts that journal at
# vol. 17 in 1949, and it was at vol. 37 in 1962); *J. Phys. Chem.* **63**
# (1959) exists and carries Hill's adsorption-definitions paper at p. 456:
# that ref. 5 intends it is an inference, stated as one.
#
# `J1.1` (Langmuir 1918) and `J1.3` (BET 1938) are the pure-isotherm siblings;
# both were read for this page and neither is loaded, no dataset of either
# page is used and no number of theirs is retyped here.
