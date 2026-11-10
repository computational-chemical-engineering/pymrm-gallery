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
# title: "Fick's law with Wilke's mixture rule for the effective diffusivity"
# description: "Wilke's 1950 rule is exact when a species diffuses through gases that are genuinely stagnant, so reporting agreement there proves nothing. His own paper supplies the test that can fail: eight ternary film problems solved three ways, with the exact Stefan-Maxwell solution printed beside the rule. Recomputed here, the rule is wrong by tens of percent on the minority flux, and the failure is not in the harmonic average but in the composition it is frozen at, which his own second method largely repairs."
# categories: [sec:A, struct:S3, tier:T0, data:tier6, phase:gas]
# date: 2026-08-07
# ---
#
# # Fick's law with Wilke's mixture rule for the effective diffusivity
#
# **Catalog ID:** `A4.1` · **Structure:** `S3` (1-D steady BVP) · **Tier:** T0
#
# Every reactor model that carries one diffusion coefficient per species in a
# multicomponent gas is using some version of
#
# $$
# D'_A \;=\; \frac{1-y_A}{\dfrac{y_B}{D_{AB}}+\dfrac{y_C}{D_{AC}}+\dfrac{y_D}{D_{AD}}+\cdots}
# $$
#
# which is the equation printed in the abstract of Wilke (1950) and reprinted
# forty-three years later as eq. (6.1.14) of Taylor and Krishna, with
# `(Wilke, 1950)` inline. It is a weighted **harmonic** mean of the binary
# diffusivities, and Wilke derives it from the Stefan–Maxwell equations for
# "the effective diffusion coefficient of a gas with respect to a multicomponent
# mixture of **stagnant** gases".
#
# That last word is the whole problem with validating it. When species $i$ really
# does diffuse through $n-1$ gases that carry no flux, the rule is not an
# approximation at all: it is an algebraic identity, and any page that reports
# agreement there is reporting arithmetic. So this page does not stop there.
#
# Wilke's own paper contains the test that can fail. He sets up **eight ternary
# test examples**, solves each one three ways (his rule, Method 1; a refinement
# of it, Method 2; and the exact solution of the Stefan-Maxwell equations,
# Method 3) and prints all of it as Table 1. This page recomputes every cell of
# that table from the printed inputs, solves the exact problem three independent
# ways (a closed form derived here, a matrix exponential, and a pymrm
# finite-volume solve), and measures how the rule fails.
#
# **The result**, printed in the headline cell below: against the exact solution
# of his own eight problems the rule is wrong by tens of percent on the species
# that carries the minority of the net flux. But the failure is not in the
# harmonic average: replacing it with the weighted *arithmetic* mean that
# J. H. Arnold attributes to Hougen and Watson in the paper's own printed
# Discussion ("the latter is badly in error in some cases") leaves the median
# error almost where it was and makes the worst case far worse. What costs the
# accuracy is the *composition* the rule is frozen at: Wilke's Method 2 uses the
# identical harmonic rule at a flux-weighted film composition instead of the
# arithmetic mean, and its worst error is a fraction of Method 1's.
#
# Three printed defects are proved from the paper's own numbers along the way,
# two printed typos are quoted with `[sic]` rather than corrected, and one
# printed claim is checked and found to hold on the author's arithmetic and
# fail on the true solution.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "main"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A4.1-wilke-mixture-rule"  # equals the page directory

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
import pymrm
import re, textwrap, time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from pymrm.checks import observed_orders
from gallery_utils import load_data, Claims
from model import (unpack, Dmat, wilke_rule, arith_rule, exact_closed, exact_expm,
                   roots_both_signs, single, exact_solve, SMFilm, method1,
                   m2_pieces, method2, set_table,
                   P_MM, XLEN, T_K, R_GAS, P_ATM, C_TOT, SCALE)

T_START = time.time()
claims = Claims(PAGE.split("-")[0])

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
np.set_printoptions(precision=6, suppress=False)
print("numpy", np.__version__, "| pymrm", pymrm.__version__)

# %% [markdown]
# ## Background
#
# ### The paper, and which Wilke it is
#
# C. R. Wilke, *"Diffusional Properties of Multicomponent Gases"*, **Chemical
# Engineering Progress 46**(2) 95–104, February 1950, presented at the Ninth
# Regional Meeting, Los Angeles. The file used here is a microfilm scan whose text
# layer is one of the two worst in this library (182 characters extract from
# journal page 97; the by-line comes back as `C.R. WI E`), so **nothing here is
# read from the text layer**. `pdfimages -list` reports JPEG grayscale at
# 400 × 400 ppi native, and every number on this page was read off a 400 dpi
# render, cropped and re-read at digit scale. The identity was confirmed from the
# document itself, title page, by-line, affiliation, running feet
# (`Vol. 46, No. 2` / `CHEMICAL ENGINEERING PROGRESS` / `February, 1950`,
# pages 95–104), the 30-item Literature Cited, and the signed Discussion, and
# independently by Taylor and Krishna's reference list.
#
# The name attaches to at least four different results and they get confused. The
# one this page is about is the **gas-phase diffusion** mixture rule of
# *Chem. Eng. Prog.* **46** (1950). It is **not**:
#
# - the Buddenberg–Wilke **viscosity** mixture rule (*Ind. Eng. Chem.* **41**,
#   1345, 1949; *J. Chem. Phys.* **18**, 517, 1950), which is what Bird, Stewart
#   and Lightfoot and Chapman and Cowling carry when they cite "Wilke", and which
#   appears in *this* paper too, as eqs. (39)–(43) and Table 2, a separate topic
#   scoped out below;
# - Wilke and Lee's binary-diffusivity **estimation** correlation
#   (*Ind. Eng. Chem.* **47**, 1253, 1955), which is a different paper, journal,
#   year and equation;
# - the Wilke–Chang liquid-diffusivity correlation.
#
# ### The second witness
#
# Taylor and Krishna (1993), *Multicomponent Mass Transfer*, Wiley, is on disk and
# is used here for three things only: it prints the rule as eq. (6.1.14) with
# `(Wilke, 1950)` inline and states the condition, quoted as printed, stray
# article and all: *"When species $i$ diffuses through **a** [sic] $n-1$ stagnant
# gases and we have $N_j = 0$, $i \neq j$"*; its section 8.6
# records the pedigree ("the oldest, simplest, and still widely used methods,
# pioneered by Hougen and Watson (1947) and by Wilke (1950)"); and its
# Example 8.5.1, explicitly "based on experiments conducted by Fairbanks and
# Wilke (1950) with a view to assessing the validity of Wilke's effective
# diffusivity formula", supplies a fully worked problem *in the stagnant regime*
# with five printed intermediates. That example is reproduced here to show what
# "exact" buys and what it does not.
#
# ### Scope
#
# **In scope:** the diffusion half of the paper, eqs. (1)–(38), Table 1 in full,
# the Example-VI sample calculation, and the Example-IV diffusivity averages.
#
# **Out of scope, and why:**
#
# - **Figure 1** (toluene into hydrogen–argon; ethyl propionate into hydrogen–air)
#   is the paper's only experimental evidence, and it is a *figure*. Its data come
#   from Fairbanks's 1948 M.S. thesis (the paper's ref. 10) and the paper says the
#   experiments "are to be described in detail in a subsequent paper". Nothing on
#   this page rests on it and it is not digitised. The paper's own Comment section
#   is blunt about the consequence: validity of the Stefan–Maxwell equations "has
#   not been established experimentally in multicomponent systems for other than
#   the stagnant gas case cited here".
# - **Figure 3** is *recomputed*, not traced: both of its curves follow from the
#   paper's own eqs. (2) and (6) and the printed conditions of Example IV, and the
#   two scalars the text prints beside it (2.15 and 2.22 sq.cm./sec.) are
#   reproduced from that computation. Never trace a curve you can compute.
# - **Gas mixture viscosities** (eqs. 39–43, Table 2), the Buddenberg–Wilke
#   result, a different catalogue case.
# - **Equipment capacity coefficients** (eqs. 44–57), the transfer-unit algebra
#   that follows from the flux equations; it introduces no new closure.

# %% [markdown]
# ## The published model
#
# Notation is Wilke's own (journal p. 103): $p$ partial pressure in atm, $y$ mole
# fraction, $N$ rate of diffusion in g.mole/(sq.cm.)(sec.), $x$ the effective film
# thickness in cm, $R$ the gas constant in (cu.cm.)(atm.)/(g.mole)(°K).
#
# ### The rule: eq. (6)
#
# For diffusion of $A$ into a mixture of gases $B, C, D, \dots$, Wilke writes the
# Maxwell friction balance (his eq. 3), sets $U_B = U_C = \cdots = 0$ ("in the case
# of all gases being stagnant except $A$"), and reaches
#
# $$
# D'_A \;=\; \frac{1-y_A}{\dfrac{y_B}{D_{AB}}+\dfrac{y_C}{D_{AC}}+\dfrac{y_D}{D_{AD}}+\cdots}
# \tag{6}
# $$
#
# with the equivalent form (6a) in which the $y$'s are mole fractions "in the inert
# gas portion considered separately from diffusing gas $A$", so that $D'_A$ "is
# constant for a given inert gas composition and independent of the mole fraction
# of $A$". A footnote records that "this result was partially anticipated by
# Stefan".
#
# ### The flux equations: eqs. (9)–(17)
#
# The defining rate equation for component $A$ of a multicomponent system is
# Fick's law plus bulk flow,
#
# $$
# N_A \;=\; -\frac{D'_A P}{RT}\frac{\mathrm{d}y_A}{\mathrm{d}x} + N_t\,y_A ,
# \tag{9}
# $$
#
# with $N_t$ the total net rate of transfer. Wilke introduces
# $\phi_A = N_A/N_t$ (eq. 10) and integrates eq. (9) at constant $D''_A$ to
#
# $$
# N_A \;=\; \frac{\phi_A D''_A P}{RTx}\,\ln\frac{\phi_A - y_{A_2}}{\phi_A - y_{A_1}},
# \tag{14}
# $$
#
# which he rearranges into the film-pressure-factor form
#
# $$
# N_A \;=\; \frac{D''_A P}{RTx\,(p_f)_A}\,(p_{A_1}-p_{A_2}),
# \qquad
# (p_f)_A = \frac{P}{\phi_A}\,(\phi_A-y_A)_{lm}.
# \tag{15, 16}
# $$
#
# Eqs. (14) and (15) are the *same* equation; (14) is the form used below because
# it stays finite as $N_t \to 0$. The two familiar limits are printed: for one gas
# through a second stagnant gas $\phi_A = 1$ and $(p_f)_A = (P-p_A)_{lm}$, giving
# the classic Stefan equation (20); for equimolal counterdiffusion
# $\phi_A = \infty$ and $(p_f)_A = P$.
#
# ### Method 1: the rule, at the arithmetic-mean composition
#
# $$
# D''_A = \frac{1-y'_A}{\dfrac{y'_B}{D_{AB}}+\dfrac{y'_C}{D_{AC}}+\cdots},
# \qquad y'_A = \tfrac12\,(y_{A_1}+y_{A_2}),\ \ \text{etc.}
# \tag{25}
# $$
#
# ### Method 2: the same rule, at an *effective film* composition
#
# Wilke's own refinement. Partial film pressure factors are defined pairwise
# (eqs. 30–32), for instance
#
# $$
# (p_f)_{AB} = \frac{P}{\phi_{AB}}\,(y'_A+y'_B)\,(\phi_{AB}-y'_{AB}),
# \qquad
# y'_{AB} = \frac{y'_A}{y'_A+y'_B},\quad \phi_{AB} = \frac{N_A}{N_A+N_B},
# $$
#
# normalised into effective mole fractions (eqs. 27–29),
# $y''_B = (p_f)_{AB}/[(p_f)_{AB}+(p_f)_{AC}+\cdots]$, and fed back into the same
# harmonic average (eq. 26):
#
# $$
# D''_A = \frac{1}{\dfrac{y''_B}{D_{AB}}+\dfrac{y''_C}{D_{AC}}+\cdots}.
# \tag{26}
# $$
#
# One branch is stated in words and is live in the data: *"it is possible to
# obtain a negative partial film pressure factor in Equations (30)-(32). In such
# cases it is recommended that the presence of the component giving the negative
# value be neglected … and its partial film pressure factor assumed equal to zero
# in Equations (27)-(29)."* Example VIII takes that branch.
#
# ### Method 3: the exact solution
#
# For simultaneous diffusion of $A$ and $B$ through a third **stagnant** gas $C$
# Wilke quotes Sherwood's differential equations,
#
# $$
# -\frac{P}{RT}\frac{\mathrm{d}p_A}{\mathrm{d}x}
# =\left[\frac{1}{D_{AB}}-\frac{1}{D_{AC}}\right]N_A p_B
# -\left[\frac{N_B}{D_{AB}}+\frac{N_A}{D_{AC}}\right]p_A + \frac{N_A P}{D_{AC}},
# \tag{21}
# $$
#
# with eq. (22) the same expression under $A \leftrightarrow B$, and Gilliland's
# integral
#
# $$
# \frac{N_A}{D_{AC}}+\frac{N_B}{D_{BC}} = \frac{P}{RTx}\,\ln\frac{p_{C_2}}{p_{C_1}}
# \tag{23}
# $$
#
# together with eq. (24). Rates "calculated by Equations (21)–(24) are tabulated
# under method 3 in Table 1, and are considered to be correct solutions to these
# examples to be used for testing the approximate methods."
#
# ### These are the Stefan–Maxwell equations, and that is checked below
#
# Writing the standard Maxwell–Stefan balance for an ideal ternary at constant
# $T$, $P$, with $c = P/RT$,
#
# $$
# \frac{\mathrm{d}y_i}{\mathrm{d}x} \;=\; \sum_{j\neq i}\frac{y_i N_j - y_j N_i}{c\,D_{ij}},
# $$
#
# and imposing $N_C = 0$ reproduces eq. (21) term for term; the derivation is
# carried out in the notebook. **It reproduces eq. (22) only after one printed
# subscript is corrected**, and eq. (23) then follows in three lines as the exact
# integral of the $C$ equation. Both facts are used below.

# %% [markdown]
# ## Parameters and assumptions
#
# All eight examples share the conditions printed in Table 1's own footnote and
# repeated on journal p. 98: one film thickness, one total pressure, one
# temperature, the system H₂, H₂O and CO₂, and three distinct binary
# diffusivities. The cell below prints the values the page uses, from `model.py`.
#
# Assumptions, all Wilke's:
#
# 1. **Steady state, one dimension, constant $T$ and $P$.** $N_A$ and $N_B$ are
#    constants of the problem.
# 2. **Gas $C$ is stagnant in every example**: $N_C = 0$. This is what makes the
#    problem two-dimensional in the unknowns and what makes eq. (23) exact.
# 3. **Ideal gas.** $c = P/RT$ is uniform.
# 4. The diffusion coefficients are *not* real: the table's own dagger footnote
#    says they "are estimated and do not correspond to true values". Nothing on
#    this page is a statement about carbon dioxide, water or hydrogen; every
#    number is a statement about the *methods*, compared with each other on
#    identical inputs.
# 5. $R$ and $T$ are chosen: the paper prints neither. Every method-versus-method
#    error is a ratio of two fluxes computed at the same $c = P/RT$, so $c$
#    cancels from it; the break table measures what a misread $R$ does to the
#    reproduction of the printed table, which is not a ratio.
#
# **Fit or test?** Neither, and the distinction still matters. Nothing on this page
# is fitted to anything: there is no adjustable parameter anywhere in Methods 1, 2
# or 3. What the page reports is a **model-versus-model** comparison on the
# author's own eight problems, using the author's own exact solution as the
# reference: a *reproduction*, not a *validation*. The only experimental content
# in the paper is Figure 1, which is scoped out, so **no number on this page is
# tested against a measurement**, and the page never says otherwise.

# %%
print(f"film thickness x = {XLEN} cm,  P = {P_MM} mm Hg = {P_ATM:.6f} atm")
print(f"T = {T_K} K (the footnote's 40 degC),  R = {R_GAS} cc atm/(g.mol K), both chosen")
print(f"c = P/RT = {C_TOT:.6e} g.mol/cm3;  Table 1 rates are printed x {SCALE:.0e}")

# %% [markdown]
# ## The data
#
# Three tier-6 datasets, all of them the authors' own computations. No page
# image is reproduced and no figure is digitised.
#
# - `wilke-1950-table1.csv`: every cell of Table 1 (journal p. 99), which is
#   printed rotated and was rotated losslessly and read cell by cell.
# - `wilke-1950-printed-intermediates.csv`: the printed steps of the
#   Example-VI sample calculation (journal p. 101), the numbers of the
#   Example-VI film-pressure-factor paragraph (p. 99), the two Example-IV
#   diffusivity averages, and the paper's own accuracy claim for Method 2.
# - `taylor-krishna-1993-example-8-5-1.csv`: Taylor and Krishna's worked
#   evaporation example, printed pages 201 to 203.
#
# Every printed number is also loaded as the literal the paper prints
# (`dtype=str`), so that "within the printed precision" is a count in half-units
# of the last printed digit, never a guess from a float.
#
# **No other page's dataset is loaded**, so the cross-page reconciliation rule does
# not apply. The related page `A4.2` reaches a Wilke-versus-Maxwell-Stefan
# comparison from Krishna and Wesselingh (1997) on the Duncan-Toor two-bulb
# experiment; it shares no dataset, no solver and no number with this page.
#
# ### Three self-consistency checks the table imposes on itself
#
# Before anything is computed *from* the table, three constraints check what was
# read *off* it. All three are properties of the printed table alone.

# %%
T1 = load_data("wilke-1950-table1.csv", page=PAGE)
# the same files re-read as text, so that print-rounding bounds are taken from
# the LITERAL the paper prints rather than inferred from a float's repr (a
# printed trailing zero, "0.330", would otherwise be lost)
T1S = load_data("wilke-1950-table1.csv", page=PAGE, dtype=str)
PI = load_data("wilke-1950-printed-intermediates.csv", page=PAGE)
PIS = load_data("wilke-1950-printed-intermediates.csv", page=PAGE, dtype=str)
TK = load_data("taylor-krishna-1993-example-8-5-1.csv", page=PAGE)
TKS = load_data("taylor-krishna-1993-example-8-5-1.csv", page=PAGE, dtype=str)
tk = dict(zip(TK.quantity, TK.value))
tks = dict(zip(TKS.quantity, TKS.value))
pi = {(r.quantity, str(r.example), str(r.method)): r.value for r in PI.itertuples()}
pis = {(r.quantity, str(r.example), str(r.method)): r.value for r in PIS.itertuples()}

BASE = T1.drop_duplicates("example").set_index("example")
EXN = list(BASE.index)
set_table(BASE)


def _decimals(lit):
    """Digits after the decimal point in a printed literal."""
    lit = str(lit).strip()
    return len(lit.split(".")[1]) if "." in lit else 0


def half_unit(lit):
    """Half a unit in the last digit a printed literal shows."""
    return 0.5 * 10.0 ** -_decimals(lit)


def half_units(calc, printed, lit):
    """|calc - printed| in half-units of the literal's last printed digit."""
    return abs(calc - printed) / half_unit(lit)


s1 = BASE[["p1_A", "p1_B", "p1_C"]].sum(axis=1)
s2 = BASE[["p2_A", "p2_B", "p2_C"]].sum(axis=1)
print(f"check 1  {2 * len(BASE)} partial-pressure sums, max |sum - {P_MM:.0f} mm| = "
      f"{np.max(np.abs(np.r_[s1.values, s2.values] - P_MM)):.3e} mm Hg")
VALS = {}
for e, r in BASE.iterrows():
    for lab in ("AB", "AC", "BC"):
        pair = frozenset((r[f"gas_{lab[0]}"], r[f"gas_{lab[1]}"]))
        VALS.setdefault(pair, set()).add(r[f"D_{lab}"])
print(f"check 2  {3 * len(BASE)} diffusivity cells -> {len(VALS)} distinct gas pairs, "
      f"each with {[len(v) for v in VALS.values()]} distinct value(s)")
for pair, v in VALS.items():
    print(f"           {'-'.join(sorted(pair)):>10} : {sorted(v)[0]:.4f} cm2/s")
# check 3: the exponent glyph.  Both the Table 1 header multiplier and the
# exponent on journal p. 101 are illegible at native resolution.  Arithmetic,
# not pixel shape, fixes them: Method 1 for Example I and the p. 101 chain are
# recomputed below and land on the printed digits only for this multiplier.
print(f"check 3  multiplier set to {SCALE:.0e} by arithmetic (see Results); "
      "the printed glyph is not legible and is not used")

# %% [markdown]
# Check 1 is exact on every sum, which is what fixes the five-figure entries: a
# single mis-read digit in any of them would leave a residue. Check 2 collapses
# the diffusivity cells onto three numbers and passes on rows VI and VII, where
# $A$ and $C$ are permuted; a mis-assignment there is the most likely
# transcription error in the table, and the break table shows the check catching
# it. Check 3 is an ambiguous glyph settled by arithmetic rather than by pixel
# shape; the two recomputations that settle it are in Results.

# %% [markdown]
# ## PyMRM implementation
#
# ### The Stefan-Maxwell film, and Wilke's eqs. (21) to (23) derived from it
#
# With $c = P/RT$ constant, $N_C = 0$, and $y_C = 1-y_A-y_B$, the Maxwell-Stefan
# balances for $A$ and $B$ are **linear in $(y_A, y_B)$ with constant
# coefficients**:
#
# $$
# -c\,\frac{\mathrm{d}y_A}{\mathrm{d}x}
# = N_A\!\left(\frac{y_B}{D_{AB}}+\frac{y_C}{D_{AC}}\right) - N_B\frac{y_A}{D_{AB}},
# \qquad
# -c\,\frac{\mathrm{d}y_B}{\mathrm{d}x}
# = N_B\!\left(\frac{y_A}{D_{AB}}+\frac{y_C}{D_{BC}}\right) - N_A\frac{y_B}{D_{AB}} .
# $$
#
# Multiplying the first by $P$ and substituting $p = Py$ gives Wilke's eq. (21)
# **term for term**. The second gives his eq. (22) only if the bracket
# $[1/D_{AB} - 1/D_{AC}]$ is read as $[1/D_{AB} - 1/D_{BC}]$; that is proved to be
# a misprint below, from the paper's own eq. (23).
#
# Adding the two and using $y_A+y_B = 1-y_C$ collapses everything except a term in
# $y_C$:
#
# $$
# c\,\frac{\mathrm{d}y_C}{\mathrm{d}x} = \left(\frac{N_A}{D_{AC}}+\frac{N_B}{D_{BC}}\right) y_C
# \;\Longrightarrow\;
# \frac{N_A}{D_{AC}}+\frac{N_B}{D_{BC}} = \frac{c}{x}\ln\frac{y_{C_2}}{y_{C_1}},
# $$
#
# which is exactly eq. (23). It is used three ways: to reduce the exact problem to
# a **one-dimensional** root find, as a constraint on the *printed* Method-3
# numbers, and as the proof that eq. (22)'s subscript is wrong.
#
# ### Three routes to the exact answer
#
# | route | what it is | what it can catch |
# |---|---|---|
# | **closed form** | the linear system integrated analytically: $y_C$ is one eigenmode with rate $\lambda_C$, and $f = N_B y_A - N_A y_B$ obeys $f' = \lambda_2 f - (Q/c)\,y_C$ with $\lambda_2 = N_t/(cD_{AB})$ | everything except an error in this derivation |
# | **matrix exponential** | $\dot{y} = My + b$ propagated with a closed-form $2\times2$ $\exp(Mx)$, no eigen-decomposition of the problem structure | an error in the closed-form derivation |
# | **pymrm finite volume** | `construct_grad` + `construct_div`, face-wise $[B]^{-1}$ with the boundary faces from `compute_boundary_values`, `NumJac` + `pymrm.newton` | a discretisation or operator misuse; refines to the other two at second order |
#
# The first two share no algebra beyond the definition of $M$; the third shares
# only the matrix $[B]$ with either. The finite-volume residual is nonlinear in
# $y$ through $[B]^{-1}$, which is why it is a Newton solve; the two analytic
# routes need no solver beyond a root find on $N_A$.
#
# **Roots, and how many.** Every flux pair on this page is a root of one scalar
# equation, found by `pymrm.checks.find_roots` on a logarithmic scan of each sign:
# $N_A$ for the exact problem (eq. 23 eliminates $N_B$), the ratio $r = N_B/N_A$
# for Methods 1 and 2 (both film equations become functions of $r$ once
# $\phi_A + \phi_B = 1$). Each problem is expected to have exactly **one** root
# in the scanned range: one physical flux pair per example. The scan reports
# every sign change it finds, so a second root would be seen, and the count is a
# claim below.

# %%
IDENT = lambda t: t  # noqa: E731  (no defect: the binaries as printed)
_SOLVED = {}


def solve_all(rule=wilke_rule, comp="mean", clip=True, Dfix=IDENT, c=C_TOT, x=XLEN,
              lim=2e-2):
    """Root lists of the exact problem, Method 1, Method 2 and the arithmetic null
    for all eight examples, under one configuration (cached per configuration)."""
    key = (rule.__name__, comp, clip, id(Dfix), c, x, lim)
    if key not in _SOLVED:
        out = {}
        for e in EXN:
            b = Dfix(tuple(unpack(e)[3:]))
            out[e] = dict(ex=exact_solve(e, lim=lim, x=x, binaries=b, c=c),
                          m1=method1(e, rule=rule, comp=comp, x=x, binaries=b, c=c)[0],
                          m2=method2(e, clip=clip, x=x, binaries=b, c=c)[0],
                          nl=method1(e, rule=arith_rule, x=x, binaries=b, c=c)[0])
        _SOLVED[key] = out
    return _SOLVED[key]


def _one(roots):
    """The pair of a problem with exactly one root; nan when it has none or several."""
    return np.array(roots[0]) if len(roots) == 1 else np.full(2, np.nan)


def amax(v):
    """Largest |v|; nan if any entry is nan (a defect that left no solution) or
    the set is empty (a defect that emptied it)."""
    v = np.abs(np.asarray(v, float))
    return float(np.max(v)) if v.size else np.nan


def amin(v):
    """Smallest |v|; nan if any entry is nan or the set is empty."""
    v = np.abs(np.asarray(v, float))
    return float(np.min(v)) if v.size else np.nan


def amed(v):
    """Median |v|; nan if any entry is nan or the set is empty."""
    v = np.abs(np.asarray(v, float))
    return float(np.median(v)) if v.size else np.nan


def pf_A_of(NA, NB, y1, y2):
    """Film pressure factor (p_f)_A of eqs. (16) and (17), mm Hg."""
    phi = NA / (NA + NB)
    lmv = (y1[0] - y2[0]) / np.log((phi - y2[0]) / (phi - y1[0]))
    return P_MM / phi * lmv


# The two alternative Method-3 N_B readings are recorded in the CSV's own
# read_note column; they are parsed from there, not retyped here.
ALT = {}
for _, r in T1.iterrows():
    m = re.search(r"alternative reading (-?[\d.]+)", str(r.read_note))
    if m:
        ALT[r.example] = float(m.group(1))

GROSS = 0.01        # an eq.-(23) miss no hand-solution slop explains
PHI_NEAR = 0.05     # "phi close to 1": the dominant species carries the net flux


def study(T=T1, pi_=None, rule=wilke_rule, comp="mean", clip=True, Dfix=IDENT,
          c=C_TOT, x=XLEN, lim=2e-2):
    """Every Table-1 claim of the page, from the printed table T and printed
    intermediates pi_, under one injected defect (none by default).

    The solver inputs (terminal compositions) always come from the table as
    read; T and pi_ change only what the results are compared WITH, and the
    data checks.  Underscored keys are tables and margins for display."""
    pi_ = pi if pi_ is None else pi_
    A = {}
    base = T.drop_duplicates("example").set_index("example")
    # ---- checks 1 and 2: properties of the printed table alone
    sums = np.r_[base[["p1_A", "p1_B", "p1_C"]].sum(axis=1).values,
                 base[["p2_A", "p2_B", "p2_C"]].sum(axis=1).values]
    A["pp_sum_resid"] = float(np.max(np.abs(sums - P_MM)))
    vals = {}
    for e, r in base.iterrows():
        for lab in ("AB", "AC", "BC"):
            pair = frozenset((r[f"gas_{lab[0]}"], r[f"gas_{lab[1]}"]))
            vals.setdefault(pair, set()).add(r[f"D_{lab}"])
    A["pair_consistent"] = bool(all(len(v) == 1 for v in vals.values()))
    # ---- the four solutions of every example
    sol = solve_all(rule, comp, clip, Dfix, c, x, lim)
    A["examples_with_one_root_each_method"] = int(
        sum(len(v) == 1 for s in sol.values() for v in s.values()))
    EX, M1, M2, NL = ({e: _one(s[k]) for e, s in sol.items()} for k in ("ex", "m1", "m2", "nl"))
    A["_EX"], A["_M1"], A["_M2"], A["_NL"] = EX, M1, M2, NL
    # ---- Table 1, cell by cell: deviation and half-units of the printed digit
    rows = []
    for e in EXN:
        ib = T.index[T.example == e][0]
        ip = T.index[(T.example == e) & (T.m3_preferred == 1)][0]
        for i, spc in enumerate("AB"):
            row = dict(example=e, species=spc, gas=base.loc[e, f"gas_{spc}"])
            for m, idx, calc in (("m1", ib, M1[e][i]), ("m2", ib, M2[e][i]),
                                 ("m3", ip, EX[e][i])):
                prt, lit = T.loc[idx, f"{m}_N{spc}"], T1S.loc[idx, f"{m}_N{spc}"]
                row[f"{m}_print"], row[f"{m}_calc"] = prt, calc * SCALE
                row[f"{m}_dev"] = calc * SCALE / prt - 1.0
                row[f"{m}_hu"] = half_units(calc * SCALE, prt, lit)
            row["null_calc"] = NL[e][i] * SCALE
            rows.append(row)
    R = pd.DataFrame(rows)
    A["_R"] = R
    defect = ((R.example.isin(["I", "II"]) & (R.species == "B"))
              | ((R.example == "VI") & (R.species == "A")))
    A["_m1_defect"] = defect
    A["table1_method1_repro_max_reldev"] = amax(R.m1_dev[~defect])
    A["table1_exampleVI_method1_NA_dev"] = amax(R.m1_dev[(R.example == "VI") & (R.species == "A")])
    A["table1_method2_repro_max_reldev"] = amax(R.m2_dev)
    m3 = np.abs(R.m3_dev.values)
    order = np.argsort(-np.nan_to_num(m3, nan=np.inf))
    A["table1_method3_repro_max_reldev"] = float(m3[order[0]])
    A["table1_method3_repro_second_max"] = float(m3[order[1]])
    A["table1_method3_worst_cell"] = f"{R.example[order[0]]}-{R.species[order[0]]}"
    for m, mask in (("m1", ~defect), ("m2", R.m2_dev.notna() | True), ("m3", R.m3_dev.notna() | True)):
        hu = R[f"{m}_hu"][mask].values
        A[f"table1_{m}_cells_within_print"] = int(np.sum(hu <= 1.0))
        A[f"_{m}_margin"] = float(np.min(np.abs(hu - 1.0)))
    # ---- defect 1: the lost minus signs of Method 1 N_B
    sd = [e for e in ("I", "II", "III") if (base.loc[e, "m1_NB"] > 0) != (M1[e][1] > 0)]
    A["method1_sign_defect_rows"] = ",".join(sd)
    A["method1_sign_defect_magnitude_dev"] = (
        max(abs(abs(M1[e][1] * SCALE / base.loc[e, "m1_NB"]) - 1) for e in sd) if sd else np.nan)
    # ---- defect 2: Example VI, Method 1, against the paper's own p. 99 paragraph
    _, y1v, y2v, *_ = unpack("VI")
    pf_stated = pi_[("p_f_A", "VI", "1")]
    A["pf_A_exampleVI_method1"] = float(pf_A_of(*M1["VI"], y1v, y2v))
    A["pf_A_from_printed_cell"] = float(pf_A_of(base.loc["VI", "m1_NA"] / SCALE,
                                                base.loc["VI", "m1_NB"] / SCALE, y1v, y2v))
    A["pf_A_exampleVI_method1_reldev"] = abs(A["pf_A_exampleVI_method1"] / pf_stated - 1)
    A["pf_A_exampleVI_half_units"] = half_units(A["pf_A_exampleVI_method1"], pf_stated,
                                                pis[("p_f_A", "VI", "1")])
    ym = 0.5 * (y1v + y2v) * P_MM
    para = [half_units(v, pi_[(q, "VI", "1")], pis[(q, "VI", "1")])
            for v, q in ((ym[1], "p_B_mean"), (ym[2], "p_C_mean"), (ym[1] + ym[2], "p_BC_mean_sum"))]
    A["pf_paragraph_within_print"] = int(sum(h <= 1.0 for h in para))
    A["_para_margin"] = float(np.min(np.abs(np.array(para) - 1.0)))
    A["_para"] = (ym, para)
    # ---- the headline: each method against the exact solution
    H = []
    for e in EXN:
        ex = EX[e]
        for i, spc in enumerate("AB"):
            H.append(dict(example=e, species=spc, gas=base.loc[e, f"gas_{spc}"],
                          exact=ex[i] * SCALE, share=abs(ex[i]) / np.max(np.abs(ex)),
                          phi=ex[i] / ex.sum(), m1=M1[e][i] * SCALE, m2=M2[e][i] * SCALE,
                          null=NL[e][i] * SCALE))
    H = pd.DataFrame(H)
    for k in ("m1", "m2", "null"):
        H[f"{k}_err"] = H[k] / H["exact"] - 1.0
    H = H.sort_values("share").reset_index(drop=True)
    A["_H"] = H
    minor, major = H[H.share < 0.5], H[H.share >= 0.5]
    A["wilke_rule_max_err_vs_exact"] = amax(H.m1_err)
    A["wilke_rule_median_err_vs_exact"] = amed(H.m1_err)
    A["wilke_rule_max_err_majority_flux"] = amax(major.m1_err)
    A["wilke_rule_median_err_majority_flux"] = amed(major.m1_err)
    A["wilke_rule_max_err_minority_flux"] = amax(minor.m1_err)
    A["wilke_rule_median_err_minority_flux"] = amed(minor.m1_err)
    A["wilke_rule_min_err_minority_flux"] = amin(minor.m1_err.values)
    for e in ("IV", "VII"):
        A[f"wilke_rule_err_{e}_B"] = float(H[(H.example == e) & (H.species == "B")].m1_err.iloc[0])
    A["method2_max_err_vs_exact"] = amax(H.m2_err)
    A["method2_median_err_vs_exact"] = amed(H.m2_err)
    A["arithmetic_null_max_err_vs_exact"] = amax(H.null_err)
    A["arithmetic_null_median_err_vs_exact"] = amed(H.null_err)
    dom = H[H.share == 1.0]
    near = (dom.phi - 1.0).abs().values
    A["examples_dominant_phi_near_1"] = int(np.sum(near < PHI_NEAR))
    A["_near_margin"] = float(np.min(np.abs(near - PHI_NEAR)) / PHI_NEAR)
    A["_near"] = list(dom.example[near < PHI_NEAR])
    A["wilke_rule_max_err_dominant_phi_near_1"] = amax(dom.m1_err[near < PHI_NEAR])
    # ---- the paper's own claim for Method 2, both readings
    cl = pi_[("method2_max_deviation_claim", "-", "2")] / 100.0

    def own_max(alt=None):
        d, at = [], []
        for e in EXN:
            ib = T.index[T.example == e][0]
            ip = T.index[(T.example == e) & (T.m3_preferred == 1)][0]
            for spc in "AB":
                m3v = ALT[e] if (alt == e and spc == "B") else T.loc[ip, f"m3_N{spc}"]
                d.append(abs(T.loc[ib, f"m2_N{spc}"] / m3v - 1.0))
                at.append(f"{e}-{spc}")
        return float(np.max(d)), at[int(np.argmax(d))]

    own, own_at = own_max()
    alts = {a: own_max(a) for a in ALT}
    worst_read = max([own] + [v for v, _ in alts.values()])
    A["_own"] = (own, own_at, alts, cl)
    A["method2_max_dev_on_printed_columns"] = own
    A["method2_max_dev_on_printed_alt_reading"] = max(v for v, _ in alts.values())
    A["method2_claim_holds_on_printed"] = bool(own < cl)
    A["_m_holds_printed"] = abs(cl - own) / cl
    A["method2_claim_holds_every_reading"] = bool(worst_read < cl)
    A["_m_holds_every"] = abs(cl - worst_read) / cl
    A["method2_claim_holds_vs_exact"] = bool(A["method2_max_err_vs_exact"] < cl)
    A["_m_holds_exact"] = abs(cl - A["method2_max_err_vs_exact"]) / cl
    # ---- Example IV: eq. (25) and the path average along the exact profile
    _, y1_4, y2_4, *b4 = unpack("IV")
    b4 = Dfix(tuple(b4))
    yref4 = {"mean": 0.5 * (y1_4 + y2_4), "point1": y1_4, "point2": y2_4}[comp]
    d25 = float(rule(yref4, 1, Dmat(*b4)))
    NA4, NB4 = EX["IV"]

    def DpB(xx):
        yA, yB = exact_closed(NA4, NB4, y1_4, *b4, x=xx, c=c)
        yA, yB = float(np.ravel(yA)[0]), float(np.ravel(yB)[0])
        return (1 - yB) / (yA / b4[0] + (1 - yA - yB) / b4[2])     # eq. (6) for H2O

    dpath = (x / quad(lambda s: 1.0 / DpB(s), 0.0, x, limit=400, epsabs=1e-14,
                      epsrel=1e-12)[0]) if np.isfinite(NA4) else np.nan
    A["_DpB"] = DpB
    k25, kpath = ("D_pp_H2O_eq25", "IV", "1"), ("D_pp_H2O_graphical", "IV", "-")
    A["eq25_DppB_exampleIV"], A["path_average_DppB_exampleIV"] = d25, float(dpath)
    A["eq25_DppB_exampleIV_reldev"] = abs(d25 / pi_[k25] - 1)
    A["path_average_DppB_exampleIV_reldev"] = abs(dpath / pi_[kpath] - 1)
    A["eq25_DppB_exampleIV_half_units"] = half_units(d25, pi_[k25], pis[k25])
    A["path_average_DppB_exampleIV_half_units"] = half_units(dpath, pi_[kpath], pis[kpath])
    # ---- eq. (23) as a constraint on every printed Method-3 pair
    e23 = []
    for idx, r in T.iterrows():
        if pd.isna(r.m3_NA):
            continue
        _, y1, y2, *bb = unpack(r.example)
        DAB, DAC, DBC = Dfix(tuple(bb))
        rhs = c / x * np.log(y2[2] / y1[2])
        lhs = (r.m3_NA / DAC + r.m3_NB / DBC) / SCALE
        lit = T1S.loc[idx]
        # half a unit in the last printed digit of N_A and N_B, through eq. (23)
        bnd = (half_unit(lit.m3_NA) / DAC + half_unit(lit.m3_NB) / DBC) / SCALE / abs(rhs)
        # ...and of the two printed p_C, which enter through ln(p_C2/p_C1)
        h1, h2 = half_unit(lit.p1_C), half_unit(lit.p2_C)
        L0 = np.log(y2[2] / y1[2])
        lo = np.log((y2[2] * P_MM - h2) / (y1[2] * P_MM + h1))
        hi = np.log((y2[2] * P_MM + h2) / (y1[2] * P_MM - h1))
        pbnd = max(abs(lo / L0 - 1), abs(hi / L0 - 1))
        e23.append(dict(example=r.example, m3_set=int(r.m3_set), NA=r.m3_NA, NB=r.m3_NB,
                        resid=lhs / rhs - 1.0, round_bound=bnd, full_bound=bnd + pbnd,
                        rhs=rhs, DAC=DAC, DBC=DBC))
    E23 = pd.DataFrame(e23)
    A["_E23"] = E23
    ar = E23.resid.abs()
    A["eq23_pairs_inside_1pct"] = int(np.sum(ar <= GROSS))
    A["_e23_1pct_margin"] = float(np.min(np.abs(ar - GROSS)) / GROSS)
    A["eq23_pairs_inside_round_bound"] = int(np.sum(ar <= E23.round_bound))
    A["_e23_round_margin"] = float(np.min(np.abs(ar / E23.round_bound - 1)))
    A["eq23_pairs_inside_full_bound"] = int(np.sum(ar <= E23.full_bound))
    A["_e23_full_margin"] = float(np.min(np.abs(ar / E23.full_bound - 1)))
    A["eq23_max_resid_consistent_rows"] = amax(E23.resid[ar <= GROSS])
    r1 = E23[(E23.example == "I") & (E23.m3_set == 3)].iloc[0]
    r3 = E23[E23.example == "III"].iloc[0]
    A["eq23_exI_set3_resid"], A["eq23_exI_set3_round_bound"] = r1.resid, r1.round_bound
    A["eq23_exIII_resid"], A["eq23_exIII_round_bound"] = r3.resid, r3.round_bound
    A["eq23_exI_set3_NB_required"] = (r1.rhs - r1.NA / SCALE / r1.DAC) * r1.DBC * SCALE
    A["eq23_exI_set3_NA_required"] = (r1.rhs - r1.NB / SCALE / r1.DBC) * r1.DAC * SCALE
    return A


S0 = study()
print(f"{S0['examples_with_one_root_each_method']} of {4 * len(EXN)} problems (exact, "
      "Method 1, Method 2 and the null, for each example) have exactly one root in "
      "the scanned range")

# %% [markdown]
# ## Results
#
# ### The highest-ranked validation first: the sample calculation
#
# Journal p. 101 works Example VI by Method 2 and prints **every intermediate**.
# That is a rank-1 validation (a worked example with printed intermediates), so
# it comes before any comparison the page constructs for itself. The chain is
# reproduced with $\phi_A$ taken from the paper's own assumed pair ($N'_B$ and
# $N'_A$ of the first approximation, printed in the data file), so each printed
# line is checked against the same input the paper used, not against a
# converged answer.

# %%
_, y1_VI, y2_VI, DAB6, DAC6, DBC6 = unpack("VI")
yref6 = 0.5 * (y1_VI + y2_VI)
NpB, NpA = pi[("N_prime_B", "VI", "2")], pi[("N_prime_A_first", "VI", "2")]
phiA6 = NpA / (NpB + NpA)
SAMPLE_KEYS = ["phi_A", "p_f_A", "p_f_AB", "p_f_AC", "y_pp_B", "y_pp_C", "D_pp_A",
               "N_prime_A", "N_A"]


def sample_chain(Rg_=R_GAS, T_=T_K, x=XLEN, D_override=IDENT):
    """The printed steps of the p. 101 sample calculation, as a function of
    every input a misreading could damage."""
    DAB, DAC, DBC = D_override((DAB6, DAC6, DBC6))
    D = Dmat(DAB, DAC, DBC)
    lm_ = (y1_VI[0] - y2_VI[0]) / np.log((phiA6 - y2_VI[0]) / (phiA6 - y1_VI[0]))
    pf_A = P_MM / phiA6 * lm_                                  # eqs. (16) + (17)
    p = m2_pieces(yref6, phiA6, D)
    Npr = p[0] * (y1_VI[0] - y2_VI[0]) * P_MM / pf_A           # eq. (15), N' form
    Nabs = Npr * P_ATM / (Rg_ * T_ * x)
    return dict(zip(SAMPLE_KEYS, (phiA6, pf_A, p[("pf", 0)][1], p[("pf", 0)][2],
                                  p[("ypp", 0)][1], p[("ypp", 0)][2], p[0], Npr,
                                  Nabs * SCALE)))


def sample_study(**kw):
    """Largest relative deviation and the count within the printed digit."""
    ch = sample_chain(**kw)
    dev = [abs(ch[k] / pi[(k, "VI", "2")] - 1.0) for k in SAMPLE_KEYS]
    hu = np.array([half_units(ch[k], pi[(k, "VI", "2")], pis[(k, "VI", "2")])
                   for k in SAMPLE_KEYS])
    return {"sample_calc_max_reldev": float(np.max(dev)),
            "sample_calc_within_print": int(np.sum(hu <= 1.0)),
            "_margin": float(np.min(np.abs(hu - 1.0))), "_chain": ch, "_hu": hu, "_dev": dev}


SC = sample_study()
print(f"{'quantity':>12}{'printed':>12}{'recomputed':>14}{'rel dev':>11}{'half-units':>12}")
for k, d, h in zip(SAMPLE_KEYS, SC["_dev"], SC["_hu"]):
    print(f"{k:>12}{pis[(k, 'VI', '2')]:>12}{SC['_chain'][k]:14.6g}{d:11.2%}{h:12.2f}")
print(f"\n{SC['sample_calc_within_print']} of {len(SAMPLE_KEYS)} printed intermediates "
      "within half a unit of their last printed digit; max relative deviation "
      f"{SC['sample_calc_max_reldev']:.2%}")
print(f"The last line closes on {SC['_chain']['N_A'] / SCALE:.4g} g.mol/(s cm2), the printed "
      f"{pis[('N_A', 'VI', '2')]} only with the multiplier {SCALE:.0e}: that fixes the "
      "illegible exponent glyph.")

# %% [markdown]
# ### Table 1, recomputed cell by cell
#
# Deviations are relative; `hu` is the distance in half-units of the last digit
# the table prints, so `hu <= 1` means the recomputation rounds to the printed
# entry.

# %%
R = S0["_R"]
print("printed (as read) vs recomputed, all x 1e5 g.moles/(sec.)(sq.cm.)\n")
pd.set_option("display.width", 200, "display.max_columns", 40)
print(R[["example", "species", "gas", "m1_print", "m1_calc", "m1_dev", "m1_hu",
         "m2_print", "m2_calc", "m2_dev", "m2_hu", "m3_print", "m3_calc", "m3_dev",
         "m3_hu"]]
      .to_string(index=False,
                 formatters={**{c: (lambda v: f"{v:+.2%}") for c in ("m1_dev", "m2_dev", "m3_dev")},
                             **{c: (lambda v: f"{v:7.1f}") for c in ("m1_hu", "m2_hu", "m3_hu")}},
                 float_format=lambda v: f"{v:9.4f}"))
print(f"\nFew cells round to the printed entry: {S0['table1_m1_cells_within_print']} of "
      f"the {int((~S0['_m1_defect']).sum())} Method-1 cells (the three printed defects "
      "below left out), at the page's chosen R and T, which the paper does not print; "
      f"largest deviation {S0['table1_method1_repro_max_reldev']:.2%}")
for m, lab in (("m2", "Method 2"), ("m3", "Method 3")):
    print(f"{lab}: {S0[f'table1_{m}_cells_within_print']} of {len(R)} cells round to the "
          f"printed entry, at the page's chosen R and T; largest deviation "
          f"{S0[f'table1_method{m[1]}_repro_max_reldev']:.2%}, nothing left out")
print(f"Method 3's largest miss is cell {S0['table1_method3_worst_cell']}; the next "
      f"largest is {S0['table1_method3_repro_second_max']:.2%}")
print(f"\nMethod 1 for Example I lands on {S0['_M1']['I'][0]:.4g} g.mol/(s cm2) against a "
      f"printed {T1S.m1_NA.iloc[0]} x 1e5: the second recomputation fixing the multiplier.")

# %% [markdown]
# Three things in that table are not agreement, and each is established below
# on the paper's own numbers rather than asserted:
#
# 1. **Examples I and II, Method 1, $N_B$**: the printed magnitudes are right and
#    the **minus signs are missing**. Example III, which has the identical
#    geometry, prints its minus sign.
# 2. **Example VI, Method 1, $N_A$**: the printed entry is not what Method 1
#    gives, and the paper's own text says so.
# 3. **Example I, Method 3**: three solution sets are printed. Only one is a
#    solution of the exact problem.
#
# A fourth entry is *not* a defect and is not excluded from anything: Example VI's
# **Method-3** $N_A$, the largest Method-3 miss in the table. It is Wilke's own
# hand solution going soft on the **minority** flux, the same species the page's
# headline is about, and it is carried in the Method-3 reproduction at full
# weight.

# %%
# ---------------- defect 1: the two missing minus signs ---------------------
M1 = S0["_M1"]
print("defect 1: Method 1, N_B, Examples I and II\n")
for e in ("I", "II", "III"):
    _, y1, y2, *_ = unpack(e)
    print(f"  Example {e:>3}: p_B goes {y1[1]*P_MM:7.3f} -> {y2[1]*P_MM:7.3f} mm "
          f"(B moves in -x, so N_B < 0) | printed m1_NB = {BASE.loc[e, 'm1_NB']:+8.4f} | "
          f"recomputed {M1[e][1]*SCALE:+8.4f}")
print(f"\n  rows whose printed sign contradicts the recomputation: "
      f"{S0['method1_sign_defect_rows']}")
print(f"  magnitudes agree to {S0['method1_sign_defect_magnitude_dev']:.2%}, so this is a "
      "lost minus sign, not a wrong number.")

# ---------------- defect 2: Example VI, Method 1 ---------------------------
print("\n\ndefect 2: Method 1, N_A, Example VI, against the paper's own p. 99 text\n")
pf_stated = pi[("p_f_A", "VI", "1")]
print(f"  paper, journal p. 99: '(the) film pressure factor for the diffusion of A "
      f"as calculated\n  under Method 1 is {pis[('p_f_A', 'VI', '1')]} mm'\n")
print(f"  from the printed Table-1 pair ({BASE.loc['VI', 'm1_NA']:+.3f}, "
      f"{BASE.loc['VI', 'm1_NB']:+.2f}) : (p_f)_A = {S0['pf_A_from_printed_cell']:7.1f} mm")
print(f"  from Method 1 solved here  ({M1['VI'][0]*SCALE:+.3f}, {M1['VI'][1]*SCALE:+.2f}) : "
      f"(p_f)_A = {S0['pf_A_exampleVI_method1']:7.1f} mm, "
      f"{S0['pf_A_exampleVI_method1_reldev']:.2%} or {S0['pf_A_exampleVI_half_units']:.1f} "
      "half-units from the text")
ym, para = S0["_para"]
print(f"\n  the same paragraph's other three numbers, recomputed: "
      f"y'_B P = {ym[1]:.2f} mm (printed {pis[('p_B_mean', 'VI', '1')]}), "
      f"y'_C P = {ym[2]:.2f} mm (printed {pis[('p_C_mean', 'VI', '1')]}), "
      f"sum = {ym[1] + ym[2]:.2f} mm (printed {pis[('p_BC_mean_sum', 'VI', '1')]})")
print(f"  {S0['pf_paragraph_within_print']} of 3 round to the printed value: the paragraph "
      "is internally consistent and agrees with the recomputation; the TABLE CELL is the outlier.")
print("  Reported, not repaired: the CSV stores the cell as printed.")

# %%
# ---------------- defect 3: eq. (22)'s subscript ---------------------------
import sympy as smp
yA, yB, NAs, NBs, dAB, dAC, dBC = smp.symbols("y_A y_B N_A N_B D_AB D_AC D_BC")
yC = 1 - yA - yB
ms21 = yB*NAs/dAB - yA*NBs/dAB + yC*NAs/dAC           # Maxwell-Stefan, N_C = 0
ms22 = yA*NBs/dAB - yB*NAs/dAB + yC*NBs/dBC


def rhs21(d=dAC):
    """Wilke eq. (21) in mole fractions, its first bracket [1/D_AB - 1/d]."""
    return (1/dAB - 1/d)*NAs*yB - (NBs/dAB + NAs/dAC)*yA + NAs/dAC


def rhs22(d=dAC):
    """Wilke eq. (22) in mole fractions, its first bracket [1/D_AB - 1/d]."""
    return (1/dAB - 1/d)*NBs*yA - (NAs/dAB + NBs/dBC)*yB + NBs/dBC


def eq22_closes(bracket21=dAC):
    """Which reading of eq. (22)'s first bracket makes (21) + (22) integrate to eq. (23)."""
    ok = [tag for tag, d in (("D_AC", dAC), ("D_BC", dBC))
          if smp.simplify(smp.expand(rhs21(bracket21) + rhs22(d))
                          - (NAs/dAC + NBs/dBC)*yC) == 0]
    return ",".join(ok) or "none"


print("defect 3: eq. (22) as printed is inconsistent with eq. (23)\n")
print("  eq. (21) printed:  ... = [1/D_AB - 1/D_AC] N_A p_B - [N_B/D_AB + N_A/D_AC] p_A + N_A P/D_AC")
print("  eq. (22) printed:  ... = [1/D_AB - 1/D_AC] N_B p_A - [N_A/D_AB + N_B/D_BC] p_B + N_B P/D_BC")
print("  (21) minus Maxwell-Stefan :", smp.simplify(rhs21() - ms21))
print("  (22) as printed minus MS  :", smp.simplify(rhs22(dAC) - ms22))
print("  (22) with D_AC -> D_BC    :", smp.simplify(rhs22(dBC) - ms22))
EQ22 = eq22_closes()
print(f"  readings of (22)'s first bracket for which (21) + (22) integrate to eq. (23): {EQ22}")
print("\n  Eq. (23) is the equation Wilke says Method 3 was solved from, and most printed"
      " Method-3\n  pairs satisfy it (V2 below), so the free symbol is the subscript in (22)."
      "  Reported, not\n  repaired (the CSV holds no equations); every computation here uses"
      " the corrected form.")

# %% [markdown]
# ### The headline: where the rule is not exact, and how it fails
#
# The comparison is against the exact solution computed here, not against Wilke's
# Method-3 column; his column is used separately, as a transcription check.
# Species $C$ carries no flux in every example, so the departure from the rule's
# stated condition is measured by the **flux share** of each species,
# $|N_i|/\max(|N_A|,|N_B|)$, and by $\phi_i = N_i/N_t$.

# %%
H = S0["_H"]
print("Wilke's rule (Method 1), his refinement (Method 2) and the arithmetic-mean")
print("null, each against the exact solution of the SAME eight problems.\n")
print(H[["example", "species", "gas", "share", "phi", "exact", "m1", "m1_err",
         "m2", "m2_err", "null", "null_err"]].to_string(
    index=False,
    formatters={"m1_err": lambda v: f"{v:+.2%}", "m2_err": lambda v: f"{v:+.2%}",
                "null_err": lambda v: f"{v:+.2%}", "share": lambda v: f"{v:.4f}",
                "phi": lambda v: f"{v:8.4f}"},
    float_format=lambda v: f"{v:9.4f}"))
nmin, nmaj = int((H.share < 0.5).sum()), int((H.share >= 0.5).sum())
print(f"\n  Wilke rule (Method 1)   : max {S0['wilke_rule_max_err_vs_exact']:.1%}, "
      f"median {S0['wilke_rule_median_err_vs_exact']:.1%}")
print(f"  same rule, film comp.   : max {S0['method2_max_err_vs_exact']:.1%}, "
      f"median {S0['method2_median_err_vs_exact']:.1%}   (Method 2)")
print(f"  arithmetic-mean null    : max {S0['arithmetic_null_max_err_vs_exact']:.1%}, "
      f"median {S0['arithmetic_null_median_err_vs_exact']:.1%}")
print(f"\n  Method 1, minority-flux species (share < 0.5, n = {nmin}): median "
      f"{S0['wilke_rule_median_err_minority_flux']:.1%}, from "
      f"{S0['wilke_rule_min_err_minority_flux']:.1%} to {S0['wilke_rule_max_err_minority_flux']:.1%}")
print(f"  Method 1, majority-flux species (share >= 0.5, n = {nmaj}): median "
      f"{S0['wilke_rule_median_err_majority_flux']:.1%}, max {S0['wilke_rule_max_err_majority_flux']:.1%}")
print(f"  signed, minority water: Example IV {S0['wilke_rule_err_IV_B']:+.1%}, "
      f"Example VII {S0['wilke_rule_err_VII_B']:+.1%}")
print(f"\n  dominant species with |phi - 1| < {PHI_NEAR}: Examples {', '.join(S0['_near'])}; "
      f"Method 1 is within {S0['wilke_rule_max_err_dominant_phi_near_1']:.2%} on all of them")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
ax = axes[0]
for lab, col, c, mk in (("Wilke, eq. (25) at mean comp.", "m1_err", "C0", "o"),
                        ("Wilke, eq. (26) at film comp.", "m2_err", "C2", "s"),
                        ("arithmetic-mean null", "null_err", "C3", "^")):
    ax.loglog(np.maximum(H.share, 3e-3), 100 * H[col].abs().clip(lower=1e-3),
              mk, ms=6, color=c, label=lab, alpha=0.85)
ax.axhline(100 * S0["_own"][3], color="0.5", ls=":", lw=1)
ax.annotate("the paper's claimed ceiling for Method 2", xy=(3.2e-3, 12),
            fontsize=7.5, color="0.35")
ax.set_xlabel("flux share of the species,  $|N_i|\\,/\\,\\max(|N_A|,|N_B|)$")
ax.set_ylabel("| error against the exact solution |  [%]")
ax.set_ylim(3e-4, 1e3)
ax.set_title("A scalar effective diffusivity fails on the minority flux")
ax.legend(fontsize=8, loc="lower left")
ax.text(0.03, 0.965, "values below 1e-3 % are clipped to the axis floor",
        transform=ax.transAxes, fontsize=6.5, color="0.45", va="top")

ax = axes[1]
ax.axhline(0, color="0.6", lw=1)
ax.axvline(1.0, color="0.5", ls=":", lw=1)
for lab, col, c, mk in (("Method 1 (Wilke's rule)", "m1_err", "C0", "o"),
                        ("Method 2 (same rule, film comp.)", "m2_err", "C2", "s")):
    ax.plot(H.phi, 100 * H[col], mk, ms=6, color=c, label=lab, alpha=0.9)
for r in H.itertuples():                      # markers only: the points are not a curve
    if abs(r.m1_err) > 0.2:
        ax.annotate(f"{r.example}-{r.species}", (r.phi, 100 * r.m1_err),
                    textcoords="offset points", xytext=(6, -2), fontsize=7,
                    color="C0")
ax.set_xscale("symlog", linthresh=1.0)
ax.set_xlabel("$\\phi_i = N_i/N_t$   (symlog, linear within $|\\phi|<1$)")
ax.set_ylabel("signed error [%]")
ax.set_title("Signed error against the flux ratio")
ax.set_ylim(-55, 100)
ax.annotate("$\\phi_i = 1$: species $i$ carries\nthe whole net flux", xy=(1.3, 60),
            fontsize=7.5, color="0.35")
ax.legend(fontsize=8, loc="upper left")
plt.tight_layout(); plt.show()

# %% [markdown]
# ### The paper's own accuracy claim for Method 2
#
# Journal p. 100, quoted exactly as printed: *"Comparison with results of Method 3
# indicates satisfactory agreement in all examples with a maximum **derivation**
# [sic] of less than 10%."* The word is crisp at native resolution: it is the
# author's or the compositor's slip for *deviation*, and like every other printed
# defect on this page it is reported rather than repaired. That is a checkable
# statement, and it has two readings: against Wilke's own Method-3 column, or
# against the true exact solution.

# %%
own, own_at, alts, cl = S0["_own"]
ex_row = H.loc[H.m2_err.abs().idxmax()]
print(f"  claim, journal p. 100 : maximum deviation of Method 2 from Method 3 "
      f"< {pis[('method2_max_deviation_claim', '-', '2')]} %\n")
print(f"  on WILKE'S OWN two printed columns : max {own:.2%} (at {own_at})   -> claim "
      f"{'HOLDS' if S0['method2_claim_holds_on_printed'] else 'FAILS'}")
print(f"  against the EXACT solution         : max {S0['method2_max_err_vs_exact']:.2%} "
      f"(at {ex_row.example}-{ex_row.species})   -> claim "
      f"{'HOLDS' if S0['method2_claim_holds_vs_exact'] else 'FAILS'}")
print(f"  both entries at issue are minority fluxes (share {ex_row.share:.4f} for "
      f"{ex_row.example}-{ex_row.species})")
print(f"\n  The two Method-3 N_B cells with alternative readings in the CSV's read_note "
      f"({', '.join(f'{k} -> {v}' for k, v in ALT.items())}) move the first number:")
for e, (v, at) in alts.items():
    print(f"    Example {e} at its alternative reading : max {v:.2%} at {at}")
print(f"  The claim holds on every reading: "
      f"{'yes' if S0['method2_claim_holds_every_reading'] else 'no'}, with a smallest "
      f"margin of {100 * cl * S0['_m_holds_every']:.2f} points. The transcribed readings "
      "are what the page reports.")

# %% [markdown]
# ### Example IV and Figure 3, recomputed rather than traced
#
# For Example IV the paper obtains the CO₂ profile from eq. (2), plots $D'_{H_2O}$
# from eq. (6) along it, and evaluates
# $D''_{H_2O} = (x_2-x_1)\big/\!\int \mathrm{d}x/D'_{H_2O}$ **graphically**,
# "found to be 2.15 sq.cm./sec. compared with 2.22 sq.cm./sec. calculated by
# Equation (25)". Both numbers follow from the paper's own equations; no pixel is
# read from the figure.

# %%
_, y1_IV, y2_IV, DAB4, DAC4, DBC4 = unpack("IV")
DpB = S0["_DpB"]
EX = S0["_EX"]
H4 = float(H[(H.example == "IV") & (H.species == "B")].m1_err.iloc[0])
print(f"  D''_H2O from eq. (25), arithmetic-mean composition : "
      f"{S0['eq25_DppB_exampleIV']:.4f} (printed {pis[('D_pp_H2O_eq25', 'IV', '1')]}): "
      f"{S0['eq25_DppB_exampleIV_reldev']:.3%}, {S0['eq25_DppB_exampleIV_half_units']:.2f} half-units")
print(f"  D''_H2O from the path integral (harmonic average)  : "
      f"{S0['path_average_DppB_exampleIV']:.4f} (printed {pis[('D_pp_H2O_graphical', 'IV', '-')]}): "
      f"{S0['path_average_DppB_exampleIV_reldev']:.3%}, "
      f"{S0['path_average_DppB_exampleIV_half_units']:.2f} half-units (a graphical evaluation)")
print(f"  the paper calls the gap 'small': it is "
      f"{S0['eq25_DppB_exampleIV'] / S0['path_average_DppB_exampleIV'] - 1:+.2%}, while the "
      f"Method-1 error on Example IV's water flux is {H4:+.1%}, the paper's own conclusion "
      "that 'this small difference in D\" could not account for the large error'.")
print(f"\n  D'_H2O at the two film faces: {DpB(0.0):.4f} and {DpB(XLEN):.4f} cm2/s "
      f"(= D_BC = {DBC4:.4f} at the face where CO2 has vanished)")

xs = np.linspace(0, XLEN, 400)
fig, ax = plt.subplots(figsize=(7.0, 4.2))
ax.plot(xs, [float(np.ravel(exact_closed(*EX["IV"], y1_IV, DAB4, DAC4, DBC4, x=x)[0])[0]) * P_MM
             for x in xs], color="C0", label="$p_{CO_2}$  [mm Hg]  (left)")
ax.set_xlabel("$x$  [cm]"); ax.set_ylabel("$p_{CO_2}$  [mm Hg]")
ax2 = ax.twinx(); ax2.grid(False)
ax2.plot(xs, [DpB(x) for x in xs], color="C1", label="$D'_{H_2O}$  (right)")
ax2.axhline(S0["eq25_DppB_exampleIV"], color="C3", ls="--", lw=1.2,
            label=f"eq. (25): {S0['eq25_DppB_exampleIV']:.3f}")
ax2.axhline(S0["path_average_DppB_exampleIV"], color="C2", ls=":", lw=1.6,
            label=f"path average: {S0['path_average_DppB_exampleIV']:.3f}")
ax2.set_ylabel("$D'_{H_2O}$  [sq.cm./sec.]"); ax2.set_ylim(1.6, 4.0)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right")
ax.set_title("Example IV, recomputed from eqs. (2) and (6); Fig. 3 is not digitised")
plt.tight_layout(); plt.show()

# %% [markdown]
# ## Validation
#
# ### V1: the three exact routes, and the grid order
#
# The closed form and the matrix exponential are two integrations of the same
# $2\times2$ system, each root-found on its own; their gap is bounded by the root
# tolerance. The pymrm solve is a discretisation of it: its order is read from
# successive differences of a nested cell ladder (`pymrm.checks.observed_orders`),
# never against the exact value, with the round-off companion beside it, and the
# Richardson limit of the ladder is then compared with the closed form.
#
# The ladder's round-off floor is the flux non-uniformity $\max|N/\bar N - 1|$ at
# the finest grid: at steady state $\nabla\cdot N = 0$ makes $N$ the same on every
# face, so its spread is the solve's residual re-expressed. Examples II and III
# carry the smallest fluxes in the table, their exact profiles are close to
# linear, and their grid error sits at that floor; they are kept out of the order
# study and bounded on their own.

# %%
LADDER = (40, 80, 160, 320)
NEAR_LINEAR = ("II", "III")


def route_gap(Dfix_expm=IDENT):
    """Largest relative flux gap between the closed-form and matrix-exponential
    routes, each root-found on its own (Dfix_expm: a defect in the second only)."""
    gap = 0.0
    for e in EXN:
        a = np.array(single(exact_solve(e), f"closed form {e}"))
        b = np.array(single(exact_solve(e, route=exact_expm,
                                        binaries=Dfix_expm(tuple(unpack(e)[3:]))),
                            f"matrix exponential {e}"))
        gap = max(gap, float(np.max(np.abs(b / a - 1))))
    return gap


def film_study(near_linear=NEAR_LINEAR, face="mean", nu=0, ladder=LADDER):
    """The pymrm film on a nested ladder for all eight examples, against the
    closed form; orders from successive differences."""
    rows, orders, ext, over = [], [], [], []
    err_fin, near_err = {}, []
    for e in EXN:
        ex = np.array(EX[e])
        vals = {}
        for n in ladder:
            f = SMFilm(e, n_z=n, face=face, nu=nu)
            vals[n] = f.solve() / ex
            rows.append(dict(example=e, n_z=n, err=float(np.max(np.abs(vals[n] - 1))),
                             spread=f.flux_spread, newton_its=f.its))
        err_fin[e] = rows[-1]["err"]
        if e in near_linear:
            near_err.append(err_fin[e])
            continue
        st = observed_orders(lambda n: vals[n], ladder)
        orders += list(np.ravel(st["orders"]))
        ext.append(float(np.max(np.abs(st["extrapolated"] - 1))))
        over.append(float(np.min(np.abs(vals[ladder[-1]] - vals[ladder[-2]])))
                    / rows[-1]["spread"])
    G = pd.DataFrame(rows)
    return {"pymrm_vs_closed_form_reldiff_n320": max(err_fin.values()),
            "pymrm_grid_order_min": float(np.min(orders)),
            "pymrm_grid_order_max": float(np.max(orders)),
            "pymrm_extrapolated_vs_exact": max(ext),
            "pymrm_ladder_over_roundoff": min(over),
            "pymrm_flux_uniformity": float(G.spread.max()),
            "pymrm_err_near_linear_examples": max(near_err) if near_err else np.nan,
            "_G": G}


ROUTE = route_gap()
print(f"V1a  closed form vs matrix exponential, each root-found on its own: largest "
      f"relative flux gap below 1e-10 (here {ROUTE:.1e})")
FS = film_study()
print("\nV1b  pymrm finite volume against the closed form, nested ladder\n")
print(FS["_G"].to_string(index=False, formatters={
    "err": lambda v: f"{v:.3e}", "spread": lambda v: f"(here {v:.1e})"}))
print(f"\n     observed orders (successive differences, examples other than "
      f"{', '.join(NEAR_LINEAR)}): {FS['pymrm_grid_order_min']:.4f} to "
      f"{FS['pymrm_grid_order_max']:.4f}")
print(f"     finest successive difference over the flux non-uniformity, smallest: "
      f"above 1e3 (here {FS['pymrm_ladder_over_roundoff']:.1e})")
print(f"     Richardson limit of the ladder against the closed form: "
      f"{FS['pymrm_extrapolated_vs_exact']:.1e}")
print(f"     largest relative flux error at n_z = {LADDER[-1]}: "
      f"{FS['pymrm_vs_closed_form_reldiff_n320']:.3e}")
print(f"     Examples {', '.join(NEAR_LINEAR)}: error at n_z = {LADDER[-1]} "
      f"below 1e-9 (here {FS['pymrm_err_near_linear_examples']:.1e}), at the round-off floor")

# %% [markdown]
# ### V2: Wilke's own eq. (23) as a constraint on the printed table
#
# Eq. (23) is an exact integral of the equations Method 3 solves, so every printed
# Method-3 pair must satisfy it. This is a check on the **reading**, not on the
# page's computation, and it is the constraint that decides the ambiguous glyphs.
# `round_bound` propagates half a unit in the last printed digit of $N_A$ and
# $N_B$; `full_bound` adds the same for the two printed $p_C$, which enter through
# $\ln(p_{C_2}/p_{C_1})$. Check 1 is the argument that the printed pressures are
# not rounded (every sum is exact), so `round_bound` is the one the page uses and
# `full_bound` is carried beside it.

# %%
E23 = S0["_E23"]
print(E23[["example", "m3_set", "NA", "NB", "resid", "round_bound", "full_bound"]].to_string(
    index=False, formatters={"resid": lambda v: f"{v:+.4%}", "round_bound": lambda v: f"{v:8.4%}",
                             "full_bound": lambda v: f"{v:8.4%}"},
    float_format=lambda v: f"{v:9.4f}"))
print(f"\n  {S0['eq23_pairs_inside_round_bound']} of {len(E23)} pairs sit inside round_bound and "
      f"{S0['eq23_pairs_inside_full_bound']} of {len(E23)} inside full_bound.")
print(f"  {S0['eq23_pairs_inside_1pct']} of {len(E23)} sit inside {GROSS:.0%}, the worst of those "
      f"{S0['eq23_max_resid_consistent_rows']:.2%}: Wilke solved eqs. (23) and (24) by hand, so a\n"
      "  few tenths of a percent of solution error is expected and is not a reading error.")
for r in E23[E23.resid.abs() > GROSS].itertuples():
    print(f"  Example {r.example} set {r.m3_set} misses GROSSLY: {r.resid:+.1%} against a "
          f"{r.round_bound:.2%} rounding bound ({r.full_bound:.2%} with the pressures rounded too)")
r3 = E23[E23.example == "III"].iloc[0]
print(f"  Example III's p_C moves by {abs(BASE.loc['III', 'p2_C'] - BASE.loc['III', 'p1_C']):.3f} mm "
      f"out of {P_MM:.0f}, the most ill-conditioned row, but its miss is outside even full_bound by "
      f"a factor {abs(r3.resid) / r3.full_bound:.1f}. Holding N_A needs N_B = "
      f"{(r3.rhs - r3.NA / SCALE / r3.DAC) * r3.DBC * SCALE:+.4f}, holding N_B needs N_A = "
      f"{(r3.rhs - r3.NB / SCALE / r3.DBC) * r3.DAC * SCALE:+.4f}; the exact solution "
      f"({EX['III'][0]*SCALE:+.4f}, {EX['III'][1]*SCALE:+.4f}) sits between. Stored as printed.")

_, y1I, y2I, DABI, DACI, DBCI = unpack("I")
rI = E23[(E23.example == "I") & (E23.m3_set == 3)].iloc[0]
print(f"\n  Example I set 3, printed ({rI.NA}, {rI.NB}): holding N_A, eq. (23) requires "
      f"N_B = {S0['eq23_exI_set3_NB_required']:+.4f}; holding N_B it requires "
      f"N_A = {S0['eq23_exI_set3_NA_required']:+.4f}.")
print("    The constraint says one digit is wrong but not which; neither repair is adopted.")
print(f"\n  Do the printed sets solve the boundary-value problem?  y_A at x_2 should be "
      f"{y2I[0]:.6f}:")
for r in E23[E23.example == "I"].itertuples():
    yAx = float(np.ravel(exact_expm(r.NA / SCALE, r.NB / SCALE, y1I, DABI, DACI, DBCI)[0])[0])
    print(f"    set {r.m3_set} ({r.NA:+.3f}, {r.NB:+.3f}): y_A(x_2) = {yAx:+.6f}, "
          f"miss {yAx - y2I[0]:+.6f}")          # expm route: finite at N_t = 0
rhsI = C_TOT / XLEN * np.log(y2I[2] / y1I[2])
print(f"    exact solution found here : ({EX['I'][0]*SCALE:.4f}, {EX['I'][1]*SCALE:.4f})")
print(f"    the equimolal root of eq. (23) alone, N_A = -N_B = "
      f"{rhsI / (1/DACI - 1/DBCI) * SCALE:.4f} in closed form. The paper says so itself:\n"
      "    'Selection of N_A = -N_B always satisfies Equations (23) and (24) so that the "
      "equimolal\n    roots have no unique significance.'")

# %% [markdown]
# ### V3: Taylor and Krishna's worked example, in the regime where the rule is exact
#
# Example 8.5.1 of Taylor and Krishna (1993), printed pages 201 to 203: ethyl
# propionate (1) evaporating into air (2) and hydrogen (3), $N_2 = N_3 = 0$. This
# is the case Wilke's rule was derived for, and the case the paper's own Comment
# section says is the only one with experimental support. The exact flux solves
# $1 - y_{2,0}e^{N\ell/(cD_{12})} - y_{3,0}e^{N\ell/(cD_{13})} = y_{1,\delta}$,
# whose left side falls monotonically in $N$, so it has exactly one root on the
# whole line; the scan confirms it.

# %%
RG_SI, TK_K = 8.3143, 273.15 + tk["temperature"]
Pp = tk["pressure"] * 1e3
ell = tk["film_thickness"] * 1e-3
D12, D13 = tk["D_12"] * 1e-6, tk["D_13"] * 1e-6
y0 = np.array([tk["y_1_bulk"], tk["y_2_bulk"], tk["y_3_bulk"]])
yd = np.array([tk["y_1_interface"], tk["y_2_interface"], tk["y_3_interface"]])
TK_KEYS = ["c_t", "A_11", "Phi", "Xi", "N_1"]
TK_SWEEP = (yd[0], 0.2, 0.4, 0.6, 0.8)        # interface y_1: the printed one, then wider


def tk_study(Rg_=RG_SI, T_=TK_K, ell_=ell, D_override=IDENT):
    """T&K's five printed intermediates, and the frozen-composition error against
    the exact stagnant solution, under one misread input."""
    d12, d13 = D_override((D12, D13))
    ct_ = Pp / (Rg_ * T_)
    A11_ = 0.5 * (y0[1] + yd[1]) / d12 + 0.5 * (y0[2] + yd[2]) / d13
    Phi_ = np.log((1 - yd[0]) / (1 - y0[0]))
    Xi_ = 0.5 * Phi_ * (np.exp(Phi_) + 1) / (np.exp(Phi_) - 1)
    N1_ = ct_ * Xi_ * (y0[0] - yd[0]) / (A11_ * ell_)
    ch = dict(zip(TK_KEYS, (ct_, A11_, Phi_, Xi_, N1_)))
    dev = [abs(ch[k] / tk[k] - 1.0) for k in TK_KEYS]
    hu = np.array([half_units(ch[k], tk[k], tks[k]) for k in TK_KEYS])

    def exact_N(y1x):
        f = lambda N: (1 - y0[1] * np.exp(N * ell_ / (ct_ * d12))  # noqa: E731
                       - y0[2] * np.exp(N * ell_ / (ct_ * d13)) - y1x)
        return single(roots_both_signs(f, 1e-9, 50.0, n_scan=400), f"T&K y_1 = {y1x}")

    sweep = []
    for y1x in TK_SWEEP:
        Ne = exact_N(y1x)
        y2x, y3x = y0[1] * np.exp(Ne * ell_ / (ct_ * d12)), y0[2] * np.exp(Ne * ell_ / (ct_ * d13))
        A = 0.5 * (y0[1] + y2x) / d12 + 0.5 * (y0[2] + y3x) / d13
        Nw = ct_ * ((1 - 0.5 * (y0[0] + y1x)) / A) / ell_ * np.log((1 - y1x) / (1 - y0[0]))
        sweep.append(dict(y1_interface=y1x, exact=Ne, wilke=Nw, err=Nw / Ne - 1))
    S = pd.DataFrame(sweep)
    # Wilke's eq. (14) with phi_1 = 1 and D_1 = (1 - y_1)/A_11 at the mean composition
    N1w = ct_ * ((1 - 0.5 * (y0[0] + yd[0])) / A11_) / ell_ * Phi_
    return {"taylor_krishna_851_max_reldev": float(np.max(dev)),
            "tk_within_print": int(np.sum(hu <= 1.0)),
            "tk_wilke_identity_gap": abs(N1w / N1_ - 1),
            "tk_frozen_composition_err": float(S.err.iloc[0]),
            "stagnant_frozen_composition_worst_err": float(S.err.abs().max()),
            "_margin": float(np.min(np.abs(hu - 1.0))), "_ch": ch, "_dev": dev, "_hu": hu,
            "_S": S}


TKR = tk_study()
print(f"{'quantity':>8}{'printed':>12}{'recomputed':>14}{'rel dev':>11}{'half-units':>12}")
for k, d, h in zip(TK_KEYS, TKR["_dev"], TKR["_hu"]):
    print(f"{k:>8}{tks[k]:>12}{TKR['_ch'][k]:14.7g}{d:11.3%}{h:12.2f}")
print(f"\n  {TKR['tk_within_print']} of {len(TK_KEYS)} round to the printed value; max relative "
      f"deviation {TKR['taylor_krishna_851_max_reldev']:.1e}")
print(f"\n  IDENTITY, not agreement: Xi (y_10 - y_1d) = (Phi/2)(2 - y_10 - y_1d) = (1 - y_1,mean) Phi,"
      "\n  so T&K's flux formula and Wilke's eq. (14) with phi_1 = 1 are the SAME expression in this"
      f"\n  limit: they agree to round-off (here {TKR['tk_wilke_identity_gap']:.1e}), which proves nothing"
      " about the physics.")
print("\n  The live number: the frozen-composition flux against the exact stagnant solution,"
      "\n  at the printed interface and with the driving force widened (every N_j still zero):")
print(TKR["_S"].to_string(index=False, formatters={"err": lambda v: f"{v:+.3%}"},
                          float_format=lambda v: f"{v:10.5f}"))
print(f"\n  So 'exact for stagnant gases' is a statement about the POINTWISE coefficient,"
      f" eq. (6.1.14).\n  The flux still needs an integration, and freezing the rule at one"
      f" composition costs\n  {abs(TKR['tk_frozen_composition_err']):.3%} at the printed interface and "
      f"{TKR['stagnant_frozen_composition_worst_err']:.1%} at the widest driving force.")

# %% [markdown]
# ### Claims
#
# Every number the page states, with how it is checked. Root-found and closed-form
# results are held to `FIXED`; the finite-volume result to `SOLV`. Counts over a
# set are exact claims with their margin to the threshold that decided them.
# Residuals, two-route gaps and orders are bounds with their basis.

# %%
FIXED, SOLV = 1e-6, 1e-4
# A half-unit margin of 1e-4 is about 1e-7 of a three-figure entry, far above the
# root finder's 1e-14: a count closer to its threshold than that is decided by noise.
MARGINS = {   # count or verdict -> (margin to its threshold, smallest margin allowed)
    "table1_m1_cells_within_print": (S0["_m1_margin"], 1e-4),
    "table1_m2_cells_within_print": (S0["_m2_margin"], 1e-4),
    "table1_m3_cells_within_print": (S0["_m3_margin"], 1e-4),
    "pf_paragraph_within_print": (S0["_para_margin"], 1e-4),
    "examples_dominant_phi_near_1": (S0["_near_margin"], 0.05),
    "method2_claim_holds_on_printed": (S0["_m_holds_printed"], 1e-3),
    "method2_claim_holds_every_reading": (S0["_m_holds_every"], 1e-3),
    "method2_claim_holds_vs_exact": (S0["_m_holds_exact"], 1e-3),
    "eq23_pairs_inside_1pct": (S0["_e23_1pct_margin"], 0.01),
    "eq23_pairs_inside_round_bound": (S0["_e23_round_margin"], 0.01),
    "eq23_pairs_inside_full_bound": (S0["_e23_full_margin"], 0.01),
}
NOTES = {
    "examples_with_one_root_each_method": "find_roots scans of both signs; one root expected per problem",
    "pair_consistent": "check 2: every gas pair carries one diffusivity in all rows",
    "method1_sign_defect_rows": "rows whose printed Method-1 N_B sign contradicts Method 1",
    "table1_method3_worst_cell": "example-species of the largest Method-3 miss",
    "examples_dominant_phi_near_1": f"examples whose dominant species has |phi - 1| < {PHI_NEAR}",
}
for k, v in S0.items():
    if k.startswith("_"):
        continue
    if k == "pp_sum_resid":
        claims.bound(k, v, max=1e-9, basis="chosen",
                     note="check 1; a mis-read digit of a three-decimal pressure leaves at least 1e-3 mm")
    elif isinstance(v, (bool, str, int, np.integer)) and not isinstance(v, float):
        m = MARGINS.get(k)
        claims.exact(k, v, margin=None if m is None else m[0],
                     min_margin=None if m is None else m[1], note=NOTES.get(k))
    else:
        claims.value(k, v, rtol=FIXED)
claims.exact("eq22_bracket_that_gives_eq23", EQ22,
             note="symbolic: readings of eq. (22)'s first bracket for which (21) + (22) give eq. (23)")
claims.value("sample_calc_max_reldev", SC["sample_calc_max_reldev"], rtol=FIXED)
claims.exact("sample_calc_within_print", SC["sample_calc_within_print"],
             margin=SC["_margin"], min_margin=1e-4, note="printed steps within half a unit of their last digit")
claims.value("taylor_krishna_851_max_reldev", TKR["taylor_krishna_851_max_reldev"], rtol=FIXED)
claims.exact("tk_within_print", TKR["tk_within_print"], margin=TKR["_margin"], min_margin=1e-4,
             note="T&K intermediates within half a unit of their last digit")
claims.bound("tk_wilke_identity_gap", TKR["tk_wilke_identity_gap"], max=1e-12, basis="identity",
             structural=True, note="Xi (y_10 - y_1d) = (1 - y_1,mean) Phi: the same expression; no row can move it")
claims.value("tk_frozen_composition_err", TKR["tk_frozen_composition_err"], rtol=FIXED)
claims.value("stagnant_frozen_composition_worst_err", TKR["stagnant_frozen_composition_worst_err"],
             rtol=FIXED)
claims.bound("exact_two_routes_reldiff", ROUTE, max=1e-10, basis="tolerance",
             note="two integrations of one system, each root refined to 1e-14 in Table-1 units")
claims.value("pymrm_vs_closed_form_reldiff_n320", FS["pymrm_vs_closed_form_reldiff_n320"], rtol=SOLV)
claims.bound("pymrm_grid_order_min", FS["pymrm_grid_order_min"], min=1.99, basis="theory",
             note="second-order cells; successive differences of the nested ladder")
claims.bound("pymrm_grid_order_max", FS["pymrm_grid_order_max"], max=2.02, basis="theory",
             note="second-order cells; successive differences of the nested ladder")
claims.bound("pymrm_extrapolated_vs_exact", FS["pymrm_extrapolated_vs_exact"], max=1e-7, basis="theory",
             note="Richardson limit of a second-order ladder against the closed form; next-order terms")
claims.bound("pymrm_ladder_over_roundoff", FS["pymrm_ladder_over_roundoff"], min=1e3, basis="theory",
             note="finest successive difference over the flux non-uniformity of the finest solve")
claims.bound("pymrm_flux_uniformity", FS["pymrm_flux_uniformity"], max=1e-9, basis="tolerance",
             note="div N = 0 makes N uniform; the spread is the Newton solve's residual")
claims.bound("pymrm_err_near_linear_examples", FS["pymrm_err_near_linear_examples"], max=1e-9,
             basis="theory", note="Examples II and III: near-linear exact profiles, grid error at the round-off floor")
print(f"{len(claims.metrics)} claims registered")

# %% [markdown]
# ### The defect-injection table
#
# Every claim above needs a row that fails its own check. The Table-1 claims come
# from `study()`: each row injects one defect through a call argument (the mixture
# rule, the composition it is frozen at, the Method-2 clip, a binary diffusivity,
# the gas constant, the film thickness, the scan range) or hands it a table or
# printed intermediate with one cell misread, and reruns the whole study, so
# every claim the defect reaches is tested. A study that cannot be completed
# under a defect records `nan` for its claims, which fails them. The two printed
# chains, the symbolic test, the two exact routes and the pymrm film get rows of
# their own.

# %%
BREAK, DEFINED = [], []


def brk(row, values):
    """Record one row: every claim in ``values`` against its own check."""
    DEFINED.append(row)
    for name, v in values.items():
        if name in claims.metrics:
            BREAK.append((row, name, claims.breaks(row, name, v)))


def run_row(row, fn, keys, **kw):
    """Rerun one study under a defect; nan for all its claims if it cannot finish."""
    try:
        out = fn(**kw)
    except Exception as exc:                     # noqa: BLE001 - a failed study is a result
        print(f"  {row}: study raised {type(exc).__name__}, its claims record nan")
        out = {k: np.nan for k in keys}
    brk(row, {k: v for k, v in out.items() if not k.startswith("_")})


def misread(col, e, value, m3_set=None):
    """Table 1 with one cell (or one column of an example) read as ``value``."""
    T = T1.copy()
    mask = (T.example == e) if m3_set is None else (T.example == e) & (T.m3_set == m3_set)
    T.loc[mask, col] = value
    return T


def pi_misread(key, value):
    """The printed intermediates with one number read as ``value``."""
    return {**pi, key: value}


SK = [k for k in S0 if not k.startswith("_")]
SWAP = lambda t: (t[0], t[2], t[1])                     # noqa: E731
DAB15 = lambda t: (1.5 * t[0], t[1], t[2])              # noqa: E731
for row, kw in (("harmonic rule replaced by the arithmetic mean", dict(rule=arith_rule)),
                ("rule frozen at terminal 1 instead of the mean", dict(comp="point1")),
                ("D_AC and D_BC swapped", dict(Dfix=SWAP)),
                ("D_AB x 1.5", dict(Dfix=DAB15)),
                ("Method 2 negative-p_f clip disabled", dict(clip=False)),
                (f"R = {R_GAS} read as 83.0", dict(c=P_ATM / (83.0 * T_K))),
                (f"film thickness {XLEN} cm read as 0.11 cm", dict(x=0.11)),
                ("N_A scan cut to |N_A| < 1e-6 g.mol/(s cm2)", dict(lim=1e-6))):
    run_row(row, study, SK, **kw)
_v = T1.loc[T1.example == "V", "p1_C"].iloc[0]
run_row(f"Example V p1_C {_v} read as {_v + 0.01:.3f}", study, SK, T=misread("p1_C", "V", _v + 0.01))
_T = misread("D_AC", "VI", BASE.loc["VI", "D_BC"])
_T.loc[_T.example == "VI", "D_BC"] = BASE.loc["VI", "D_AC"]
run_row("Example VI D_AC and D_BC cells read in each other's place", study, SK, T=_T)
run_row("Example I Method-1 N_B read with a minus sign", study, SK,
        T=misread("m1_NB", "I", -abs(BASE.loc["I", "m1_NB"])))
for e, col, f in (("IV", "m2_NB", 1.1), ("VI", "m2_NA", 1.1), ("II", "m3_NA", 1.1),
                  ("III", "m3_NB", 1.1)):
    v0 = T1.loc[(T1.example == e) & (T1.m3_preferred == 1), col].iloc[0]
    run_row(f"Example {e} {col} {v0} read as {v0 * f:.4g}", study, SK,
            T=misread(col, e, v0 * f) if col[1] != "3" else
            misread(col, e, v0 * f, m3_set=int(T1.loc[(T1.example == e) & (T1.m3_preferred == 1), "m3_set"].iloc[0])))
run_row(f"Example VII Method-3 N_B at its alternative reading {ALT['VII']}", study, SK,
        T=misread("m3_NB", "VII", ALT["VII"]))
run_row(f"Example VI Method-1 N_A {BASE.loc['VI', 'm1_NA']} read as {BASE.loc['VI', 'm1_NA'] * 1.1:.4g}", study, SK,
        T=misread("m1_NA", "VI", BASE.loc["VI", "m1_NA"] * 1.1))
run_row(f"Example I set 3 N_B {rI.NB} read as {rI.NB * 0.9:.3f}", study, SK,
        T=misread("m3_NB", "I", rI.NB * 0.9, m3_set=3))
for key, f in ((("p_f_A", "VI", "1"), 0.98), (("p_C_mean", "VI", "1"), 1.04),
               (("D_pp_H2O_eq25", "IV", "1"), 1.02), (("D_pp_H2O_graphical", "IV", "-"), 1.02),
               (("method2_max_deviation_claim", "-", "2"), 1.5)):
    run_row(f"printed {key[0]} {pis[key]} read as {pi[key] * f:.4g}", study, SK,
            pi_=pi_misread(key, pi[key] * f))

# --- the printed chains --------------------------------------------------------------
for row, ks, kt in ((f"R + 1.1 % ({R_GAS} -> 83.0; {RG_SI} -> {RG_SI * 83.0 / R_GAS:.4f})",
                     dict(Rg_=83.0), dict(Rg_=RG_SI * 83.0 / R_GAS)),
                    ("T + 1 % on both chains", dict(T_=T_K * 1.01), dict(T_=TK_K * 1.01)),
                    ("film thickness x 1.1 on both chains", dict(x=XLEN * 1.1), dict(ell_=ell * 1.1)),
                    ("the two diffusivities each chain uses swapped",
                     dict(D_override=SWAP), dict(D_override=lambda t: (t[1], t[0]))),
                    ("D_AB / D_12 x 1.2", dict(D_override=lambda t: (1.2 * t[0], t[1], t[2])),
                     dict(D_override=lambda t: (1.2 * t[0], t[1])))):
    run_row(row + " (sample calculation)", sample_study, ["sample_calc_max_reldev",
                                                          "sample_calc_within_print"], **ks)
    run_row(row + " (Taylor and Krishna)", tk_study, [k for k in TKR if not k.startswith("_")], **kt)

# --- the symbolic test, the two exact routes and the pymrm film --------------------------
brk("eq. (21)'s first bracket read with D_BC", {"eq22_bracket_that_gives_eq23": eq22_closes(dBC)})
run_row("matrix-exponential route run with D_AC and D_BC swapped", lambda: {
    "exact_two_routes_reldiff": route_gap(SWAP)}, ["exact_two_routes_reldiff"])
FK = [k for k in FS if not k.startswith("_")]
run_row("film divergence with nu = 1 (cylinder) instead of 0", film_study, FK, nu=1)
run_row("interior face composition from the left cell (first order)", film_study, FK, face="left")
run_row("Examples II and III kept in the order study", film_study, FK, near_linear=())

bt = pd.DataFrame(BREAK, columns=["row", "claim", "caught"])
caught = bt[bt.caught].groupby("row", sort=False).claim.apply(list)
assert len(DEFINED) == len(set(DEFINED)) == bt.row.nunique(), "every defined row is printed once"
print(f"{len(DEFINED)} rows defined, {bt.row.nunique()} in the table; claims each row fails:")
for row, names in caught.items():
    print(f"  {row}:")
    print("      " + "\n      ".join(textwrap.wrap(", ".join(names), 96)))
idle = [r for r in bt.row.unique() if r not in caught.index]
print(f"\nrows that fail no claim: {idle if idle else 'none'}")
missing = [n for n in claims.metrics
           if not claims.metrics[n].get("structural") and not bt[(bt.claim == n) & bt.caught].shape[0]]
print(f"claims no row fails: {missing if missing else 'none'}")

# %% [markdown]
# ## What pymrm adds
#
# Honestly: **not the answer**. The exact problem here is a $2\times2$ linear
# system with constant coefficients, and it has a closed-form solution that this
# page derives and uses. A finite-volume solver is not needed to reproduce
# Wilke's Table 1 and the page does not pretend it is.
#
# What the pymrm route adds is three things the closed form cannot give:
#
# 1. **An independent check on the derivation.** The closed form rests on
#    spotting that $y_C$ is an eigenmode with rate $\lambda_C$ and that
#    $f = N_B y_A - N_A y_B$ obeys a forced first-order equation with rate
#    $\lambda_2 = N_t/(cD_{AB})$. That is exactly the kind of algebra that is
#    confidently wrong. `construct_grad`, `construct_div`, `NumJac` and
#    `pymrm.newton` know none of it (only $N = -c[B]^{-1}\nabla y$ and
#    $\nabla\!\cdot\! N = 0$), and the ladder in V1 converges onto the closed form
#    at second order, with a Richardson limit that meets it. That is the check
#    that could have failed, and it is why the closed form is trusted as the
#    reference everywhere else on the page.
# 2. **The route that survives when the linearity does not.** The closed form
#    exists only because $c$, $T$, $P$ and every $D_{ij}$ are constant and $N_C$
#    is exactly zero. Drop any one of those (a reacting film, a pressure
#    gradient, composition-dependent $D_{ij}$, a fourth species with its own flux)
#    and the analytic route disappears while the pymrm residual is unchanged in
#    form. The `Reuse` section below is written around that.
# 3. **The profiles, at no extra cost.** Figure 3's two curves, the CO₂ pressure
#    profile and $D'_{H_2O}(x)$ along it, are what Wilke had to construct
#    graphically and integrate by planimeter; the page recomputes both instead of
#    digitising them.
#
# `pymrm.checks.find_roots` does the root finding: every flux pair is a scanned
# root with its count reported, never the first root an iteration happens to
# reach.
#
# What the *page* adds beyond the paper is separable from the solver, and every
# number behind it is printed in the cells above:
#
# - **The rule is measured where it is not exact.** Wilke prints the comparison but
#   never reduces it to a number; his own summary of Method 1 is qualitative
#   ("good agreement when the component under consideration is diffusing at a
#   relatively high rate compared with that of the other component, and poor
#   agreement otherwise"). Quantified, the split is by **flux share**.
# - **The failure is attributed.** Swapping the harmonic average for the
#   arithmetic one barely moves the median and makes the worst case far worse;
#   keeping the harmonic average and moving the *composition* takes the worst
#   case down to the Method-2 column. The mixture rule's functional form is not
#   what limits it.
# - **A null baseline**, taken from the paper's own printed Discussion rather than
#   invented.
# - **The paper's own accuracy claim tested both ways**: it holds on Wilke's
#   arithmetic and fails on the true solution.
# - **Three printed defects proved from the paper's own numbers**, and one
#   ambiguity left unrepaired because the constraint does not decide it.
# - **A statement of what "exact for stagnant gases" is worth**: it is a statement
#   about the pointwise coefficient. Freezing that coefficient at one composition
#   costs little on Taylor and Krishna's example and a lot when the driving force
#   is widened, with every $N_j$ still exactly zero.

# %% [markdown]
# ## Reuse
#
# **When Wilke's rule is safe.** Species $i$ carries essentially all of the net
# flux ($\phi_i \to 1$, or every other species stagnant), the driving force is
# small, and you want the flux of *that* species. On this page's eight problems
# that regime is every example whose dominant species has $\phi$ close to one
# (listed in the headline cell), and the rule is within a few percent there.
#
# **When it is not.** Anything you intend to read off the *minority* flux. The
# minority-flux entries are wrong by tens of percent, and the sign of the error
# is not even consistent (compare Examples IV and VII in the headline cell). A
# scalar effective diffusivity applied to all species at once also violates the
# flux closure it is embedded in; page `A4.2` shows the same failure from the
# other direction, on measurements, and reaches it through Krishna and
# Wesselingh's review rather than through Wilke.
#
# **The cheap fix is not a better mixture rule.** It is a better composition. The
# harmonic-versus-arithmetic choice barely moves the median error; moving from
# the arithmetic-mean composition to Wilke's flux-weighted film composition cuts
# the worst case several-fold. If a scalar diffusivity has to be used, iterate it
# on the fluxes.
#
# **Reusing the code.**
#
# - `SMFilm` solves the exact ternary film with one stagnant species. To add a
#   fourth species: widen `shape` to `(n_z, 3)`, extend `_M` with the general
#   $B_{ii} = y_i/D_{in} + \sum_{k\neq i} y_k/D_{ik}$, $B_{ij} = -y_i(1/D_{ij} -
#   1/D_{in})$, and keep `NumJac(shape, axes_diagonals=[0])`: `axes_blocks`
#   defaults to the last axis, which is the species axis, and that is what you
#   want.
# - To make a species *react* at one face, change that face's `bc` dict rather
#   than the residual: the boundary condition is on the **outward** normal, so
#   $\{a: 1, b: 0, d: q\}$ means $\mathrm{d}y/\mathrm{d}x = -q$ at the left face
#   and $+q$ at the right. The face compositions come from
#   `compute_boundary_values` under the same dict, so the mobility follows.
# - `nu=0` in `construct_div` is the Cartesian film. A spherical pellet is
#   `nu=2` and needs nothing else changed in the residual.
# - `film_roots` is reusable for any two-species film-model closure: hand it a
#   callable returning $(K_A, K_B)$ as a function of $r = N_B/N_A$ and it returns
#   **every** root on the scanned range, not the first one a Newton iteration
#   happens to find. Wilke's own paper warns that "selection of the proper set of
#   solutions is not always obvious", and Example I is why.
#
# **What this page cannot conclude.**
#
# - Nothing about real carbon dioxide, water or hydrogen. The table's own footnote
#   says the diffusivities "are estimated and do not correspond to true values",
#   so every number here is a comparison *between methods* on identical inputs.
# - Nothing measured. The paper's only experimental content is Figure 1, which is
#   scoped out; the Comment section itself says the Stefan-Maxwell equations'
#   validity "has not been established experimentally in multicomponent systems
#   for other than the stagnant gas case". The exact solution is used here as the
#   *reference*, and if the Stefan-Maxwell equations are wrong for these mixtures
#   then so is every error bar on this page.
# - Nothing about mixtures of more than three species. Every example is ternary,
#   and Wilke says so himself of Method 2: "although the method has been tested
#   only with the present three-component examples, it seems reasonable to assume
#   that equally satisfactory results would be obtained in extending the method to
#   systems of any number of components." That extrapolation is not tested here.
# - Nothing about which digit is wrong in Example I's third Method-3 set. Eq. (23)
#   proves one of the two is wrong and does not say which.

# %% [markdown]
# ## References
#
# Wilke, C. R. (1950). *Diffusional Properties of Multicomponent Gases*.
# **Chemical Engineering Progress 46**(2) 95–104. The source. Read in full on
# renders at the scan's native 400 ppi; no DOI exists. Presented at the Ninth
# Regional Meeting, Los Angeles. Every equation, every cell of Table 1 and every
# printed intermediate used here was cropped and re-read at digit scale; the text
# layer was not used for any character.
#
# Arnold, J. H. (1950), and Sage, B. H. (1950). Discussion of the above,
# **Chemical Engineering Progress 46**(2) 104. Printed with the paper. Arnold's
# contribution is the source of this page's null baseline: Wilke's expression "is
# simply a weighted harmonic mean, as opposed to the weighted arithmetic mean
# suggested by Hougen and Watson; the latter is badly in error in some cases."
#
# Taylor, R. and Krishna, R. (1993). *Multicomponent Mass Transfer*. Wiley,
# ISBN 0-471-57417-1. The independent second witness, on disk and read at its
# native 600 ppi. Section 6.1.3 (printed p. 126) prints the rule as eq. (6.1.14)
# with "(Wilke, 1950)" inline and states its condition; section 8.6 (p. 204)
# records the pedigree; Example 8.5.1 (pp. 201–203) is reproduced here in full.
#
# Gilliland, E. R. and Sherwood, T. K. (1934). *Ind. Eng. Chem.* **26**, 516;
# Sherwood, T. K. (1937). *Absorption and Extraction*, McGraw-Hill;
# Hougen, O. A. and Watson, K. M. (1947). *Chemical Process Principles III*,
# Wiley; Fairbanks, D. F. (1948). M.S. thesis, University of California;
# Stefan, J. (1871, 1872, 1874). *Sitz. Akad. Wiss. Wien* **63**, **65**, **68**;
# Maxwell, J. C. (1866, 1890). **Origins not consulted.** These are Wilke's
# refs. 12, 24, 16, 10, 29–31, 20–21, read off the Literature Cited on journal
# p. 104. Note that Wilke's text cites "(25)" for both *"Sherwood (25) presents
# differential equations"* and *"Gilliland (25) has integrated (21) and (22)"*,
# but his ref. 25 is *Sherwood, T. K., Ibid., p. 11 (1937)*, the same book as
# ref. 24, while Gilliland & Sherwood (1934) is his ref. **12**. Both citations
# of "(25)" therefore point at Sherwood's book, not at the 1934 paper.
# Eqs. (21)–(24) reach this page **through**
# Wilke's transcription of them, and eqs. (21)–(23) are independently re-derived
# here from the Maxwell–Stefan balance rather than trusted. Hougen and Watson's
# arithmetic-mean rule reaches this page through **Arnold's one-sentence
# description of it** in the printed Discussion; their book was not opened and the
# page never cites a form of it beyond what that sentence names.
#
# Krishna, R. and Wesselingh, J. A. (1997). *Chemical Engineering Science* **52**,
# 861–911. Not used here. It is the source of sibling pages `A4.2`, `A4.3` and
# `A4.4`, which reach a related conclusion about scalar effective diffusivities by
# a wholly different route and on measured data.

# %%
print(f"Runtime (here {time.time() - T_START:.0f} s)")
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
