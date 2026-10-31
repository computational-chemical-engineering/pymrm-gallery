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
# title: "Five fixed-bed models on one reactor: what each rung of Froment's ladder buys"
# description: >
#   Chapter 11 of Froment, De Wilde & Bischoff classifies fixed-bed reactor models in
#   a 2x3 table and then works five of its six cells: basic one-dimensional plug flow
#   (11.5, S2), plus axial dispersion with Danckwerts inlet conditions (11.6, S4),
#   two-dimensional with radial dispersion and a wall coefficient (11.7, S6),
#   heterogeneous with interfacial gradients (11.8, S7) and heterogeneous with
#   intraparticle profiles as well (11.9, S8). This page builds all five on ONE
#   reactor - the hydrocarbon-oxidation design case of section 11.5.2 - and measures
#   what each rung buys, as the shift in the root-found inlet partial pressure at
#   which the hot spot reaches the top of the book's own stated operating range.
# categories: [sec:D, struct:S2, struct:S4, struct:S6, struct:S7, struct:S8, tier:T0, data:tier6, phase:gas-solid]
# date: 2026-08-17
# ---
#
# # Five fixed-bed models on one reactor: what each rung of Froment's ladder buys
#
# **Catalog ID:** `D1.1` (covering `D1.2`, `D1.3`, `D1.4`, `D1.5`) &middot;
# **Structures:** `S2` `S4` `S6` `S7` `S8` &middot; **Tier:** T0
#
# ## Background
#
# Table 11.4-1 of Froment, De Wilde & Bischoff (2011), on book page 505, is not a
# ladder. It is a **2 x 3 classification**: pseudohomogeneous models ($T = T_s$,
# $C = C_s$) against heterogeneous ones, crossed with one-dimensional basic,
# one-dimensional with mixing, and two-dimensional. Read from the book's own table:
#
# | | A - pseudohomogeneous, $T = T_s$; $C = C_s$ | B - heterogeneous, $T \ne T_s$; $C \ne C_s$ |
# |---|---|---|
# | one-dimensional | I sec. 11.5 &nbsp; basic, ideal | I sec. 11.8 &nbsp; + interfacial gradients |
# | | II sec. 11.6 &nbsp; + axial mixing | II sec. 11.9 &nbsp; + intraparticle gradients |
# | two-dimensional | III sec. 11.7 &nbsp; + radial mixing | III sec. 11.10 &nbsp; + radial mixing |
#
# **The ordering this page uses, and why.** The five cases the gallery catalogue
# calls `D1.1`-`D1.5` are five of those six cells. The book presents them in the
# order 11.5, 11.6, 11.7, 11.8, 11.9 - reading the table **down each column**, so
# that the two-dimensional pseudohomogeneous model (11.7) sits between the two
# one-dimensional heterogeneous ones. This page instead orders them by **how much
# each one adds to the transport description of the previous**:
#
# 1. **11.5** basic one-dimensional plug flow (`S2`) - convection only,
# 2. **11.6** + axial dispersion of heat and mass, Danckwerts inlet (`S4`),
# 3. **11.8** + a solid phase separated from the gas by a film (`S7`),
# 4. **11.9** + concentration profiles inside the particle (`S8`),
# 5. **11.7** + a radial coordinate in the tube (`S6`).
#
# That puts 11.7 last rather than third, for one reason that is arithmetic rather
# than taste: **rung 1's wall coefficient and rung 5's two radial parameters are
# the same measurement.** Section 11.7.4's equation (11.7.4-1) collapses the
# two-dimensional pair $(\alpha_w, \lambda_{er})$ into the single $\alpha_i$ that
# rung 1 uses, and the notebook reproduces section 11.5.2's printed
# $U = 0.096$ kJ/(m$^2$ s K) from section 11.7.3's printed $\alpha_w = 0.156$ and
# $\lambda_{er} = 0.78 \times 10^{-3}$ **to both of the digits it prints**.
#
# That is a consistency check, and it is stated as one. It is **not** evidence that
# the book derived the 0.096 that way: Van Welsenaere & Froment's own
# $U = 82.7$ kcal/(m$^2$ h K), converted to kJ/(m$^2$ s K), also rounds to 0.096
# at the two significant figures actually printed, and section 11.5.2 is otherwise
# their parameter set throughout. What the identity *does* decide, discriminatingly,
# is the tube radius: (11.7.4-1) returns 0.096 at $R_t = 0.0125$ m and a value that
# rounds to 0.095 at the 0.0127 m the printed 2.54 cm diameter implies (both computed
# in *Parameters and assumptions*).
# Rung 1 and rung 5 are therefore entangled whichever way the number was obtained,
# which is the reason for reading the table in this order rather than the book's.
#
# **The point of the page is the comparison.** Five restatements of five equation
# sets would be five-sixths of a chapter retyped. What is worth measuring is
# *where the simpler model already suffices and where it misleads*, and the
# comparison is only meaningful if all five rungs run on the **same reactor with the
# same chemistry**. Section 11.5.2's hydrocarbon-oxidation design case is the one
# the book itself carries furthest: 11.5.2 designs it, 11.5.3 works its runaway
# criteria, 11.7.3 puts a fuller kinetic scheme for it in two dimensions, and 11.8.2
# runs its startup with a solid phase. So that is the reactor here, with the
# pseudo-first-order rate of 11.5.2 held fixed across all five rungs.
#
# **What this page does not own.** The chapter overlaps six published pages, and
# this page loads their numbers rather than restating them:
#
# | the chapter's content | the page that owns it | what is done here |
# |---|---|---|
# | Van Welsenaere & Froment's runaway criteria, which 11.5.3 applies | `D2.2` | both its CSVs are loaded; Example 11.5.3.A is checked *against the original*, and the criteria are not rebuilt |
# | Froment (1967)'s o-xylene network and its 1-D/2-D runaway pair, which 11.7.3/11.7.4 restate | `C2.10` (which absorbed `D3.4`) | its parameter and stated-result CSVs are loaded and reconciled with the book's restatement; the 1967 case is **not** re-run |
# | Danckwerts' inlet and outlet conditions, which 11.6 needs | `A2.1` | its boundary-condition dicts are used as published, with its outlet-sign blind spot cited, not rediscovered |
# | the isothermal effectiveness factor, which 11.9.1 needs | `B1.1` | its exact-$\eta$ CSV is loaded as the reference for eq. (11.9.1-11) |
# | the Ergun friction factor of eq. (11.5.1-13) | `A1.1` | its printed-constant CSV is loaded; Example 11.5.1.A's *equivalent particle diameter* is the part `A1.1` does not carry |
# | Sh $= 2 + 1.1\,Sc^{1/3}Re^{0.6}$ for the film | `A3.4` | used as published for $k_g$, with the analogous Nu for $h_f$ |
# | Mears' interphase temperature criterion, printed on book p. 588 | `B1.7` | evaluated here on this reactor; the criterion family is `B1.7`'s |
# | Gunn's axial and radial dispersion correlations, cited by 11.6 and 11.7.1 | `A2.6` | the book's own "Pe$_a$ between 1 and 2" is used; Gunn's forms are `A2.6`'s |
#
# **Provenance.** Tier 6. Chapter 11 reports **no experiment**. Every comparison on
# this page is against a number the book computed itself or a sentence it asserts,
# so what follows is a reproduction of the book's arithmetic and of its claims about
# its own models - never a validation against data. Said again in *Validation*.


# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.10.31.5"  # replaced by the release tag at publication
PYMRM = ">=3.2.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "D1.1-froment-fixed-bed-ladder"  # equals the page directory

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
import inspect
from dataclasses import replace

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from gallery_utils import load_data, load_meta, load_claims, Claims
from model import (KCAL, Reactor, film, eta_sphere, phi_prime, Ladder, Pellet, hotspot,
                   solve_robust, T_max_safe, CTL_JUMP, CTL_STATS, NotACrossing,
                   crossing_control, p0_critical, cross_scan, ode_reference,
                   p0_critical_ode, richardson)
from model import _shift_or_none

# DETERMINISM: nothing on this page is stochastic - no sampling, no bootstrap, no
# random initial guess.  Every fixed-inlet solve starts from a guess that is a fixed
# function of the inputs, and every threshold is followed along its branch from one
# such start inside a single call, so no state crosses from one call to the next.
# Two consecutive executions give identical content and agreement.json.
# A pandas Styler's text/plain repr is a memory address and its HTML carries a
# random id, so Stylers go through to_html() with a PINNED uuid; and ipykernel
# flushes stdout on a timer, so the flush interval is disabled.
try:
    sys.stdout.flush_interval = 1e6
except Exception:
    pass
from IPython.display import HTML, display


def show(styler):
    display(HTML(styler.to_html()))


plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
C_BLUE, C_ORANGE, C_GREEN, C_PURPLE, C_RED = ("#0072B2", "#D55E00", "#009E73",
                                              "#CC79A7", "#8B0000")
np.set_printoptions(precision=8, suppress=False)

PAGE = "D1.1-froment-fixed-bed-ladder"
D22 = "D2.2-van-welsenaere-froment-runaway"
C210 = "C2.10-o-xylene-phthalic-anhydride"
B11 = "B1.1-thiele-weisz-hicks"
A11 = "A1.1-ergun-pressure-drop"

BK = load_data("froment-2011-ch11-printed.csv", page=PAGE).set_index("key")
BKMETA = load_meta("froment-2011-ch11-printed.csv", page=PAGE)


def bk(key, field="value"):
    """A number or sentence from Chapter 11, looked up - never retyped in a cell."""
    v = BK.loc[key, field]
    return float(v) if field == "value" else str(v)


print("this page's own dataset :", BKMETA["dataset_id"], "-",
      len(BK), "rows read from Chapter 11")
print("  ", BK.kind.value_counts().to_dict())
print("borrowed, and each source page read:")
VWF_PAR = load_data("van-welsenaere-froment-1970-parameters.csv", page=D22).set_index("symbol")
VWF_EX = load_data("van-welsenaere-froment-1970-examples.csv", page=D22)
F67_PAR = load_data("froment-1967-parameters.csv", page=C210).set_index("symbol")
F67_RES = load_data("froment-1967-stated-results.csv", page=C210)
ETA_REF = load_data("isothermal-exact.csv", page=B11)
ERG_PAR = load_data("ergun-1952-parameters.csv", page=A11)
for nm, pg, df in ((load_meta("van-welsenaere-froment-1970-parameters.csv", page=D22)["dataset_id"], D22, VWF_PAR),
                   (load_meta("van-welsenaere-froment-1970-examples.csv", page=D22)["dataset_id"], D22, VWF_EX),
                   (load_meta("froment-1967-parameters.csv", page=C210)["dataset_id"], C210, F67_PAR),
                   (load_meta("froment-1967-stated-results.csv", page=C210)["dataset_id"], C210, F67_RES),
                   (load_meta("isothermal-exact.csv", page=B11)["dataset_id"], B11, ETA_REF),
                   (load_meta("ergun-1952-parameters.csv", page=A11)["dataset_id"], A11, ERG_PAR)):
    print(f"   {nm:<34s} {len(df):5d} rows   from pages/{pg}/")

# what the four source pages establish about these rows, read from their own
# committed claims rather than retyped; The data, below, says what each affects
D22_C, C210_C, B11_C, A11_C = (load_claims(pg) for pg in (D22, C210, B11, A11))
C210_U_WORST = max(abs(C210_C[k]) for k in ("u_ref_identity_rel", "u_lam075_identity_rel",
                                             "u_aw150_identity_rel"))


# %% [markdown]
# ## The published model
#
# All five rungs are transcribed from **renders** of the book's own pages, not from
# its text layer. `docs/papers-inventory.yaml` records why: every Symbol-font
# operator in this file extracts as an unmappable Private-Use-Area glyph (U+F02D
# minus, U+F03D equals, U+F02B plus) in all three `pdftotext` modes, so an extracted
# equation *looks complete and has no operators at all*. Two things were caught on
# the render during this transcription that the text layer had silently dropped: the
# leading $\varepsilon$ of eq. (11.6-1), and the Haughey--Beveridge amplitude, which
# extracted as 0.0731 and is printed 0.073.
#
# **Rung 1 - section 11.5.1, basic one-dimensional pseudohomogeneous (`S2`).**
#
# $$-\frac{d(u_s C_A)}{dz} = r_A \rho_B \tag{11.5.1-1}$$
# $$u_s \rho_g c_p \frac{dT}{dz} = (-\Delta H) r_A \rho_B - 4\frac{U}{d_t}(T - T_r) \tag{11.5.1-2}$$
# $$-\frac{dp_t}{dz} = f\frac{\rho_g u_s^2}{d_p} \tag{11.5.1-3}$$
#
# with $C_A = C_{A0}$, $T = T_0$, $p_t = p_{t0}$ at $z = 0$. Note $r_A$ is the rate of
# *disappearance*, positive. Section 11.5.2 writes the same pair in partial
# pressures for the design case, which fixes the constant-density conversion the
# whole page uses, $C_A = \rho_g p_A / (M_m p_t)$:
#
# $$u_s\frac{dp_A}{dz} + \frac{M_m p_t}{\rho_g}\rho_B\,k\,p_B^0 p_A = 0 \tag{11.5.2-1}$$
# $$u_s\rho_g c_p\frac{dT}{dz} - (-\Delta H)\rho_B k p_B^0 p_A + \frac{4U}{d_t}(T - T_r) = 0 \tag{11.5.2-2}$$
#
# $$r_A = k\,p_B^0\,p_A, \qquad \ln k = 19.837 - \frac{13{,}636}{T}$$
#
# **Rung 2 - section 11.6, + axial mixing (`S4`).**
#
# $$\varepsilon D_{ea}\frac{d^2C_A}{dz^2} - u_s\frac{dC_A}{dz} - r_A\rho_B = 0 \tag{11.6-1}$$
# $$\lambda_{ea}\frac{d^2T}{dz^2} - \rho_g u_s c_P\frac{dT}{dz} + (-\Delta H)r_A\rho_B - \frac{4U}{d_t}(T-T_r) = 0 \tag{11.6-2}$$
#
# with the conditions the book calls "those generally used", verbatim:
#
# $$u_s(C_{A0} - C_A) = -\varepsilon D_{ea}\frac{dC_A}{dz}, \qquad
#   \rho_g u_s c_p (T_0 - T) = -\lambda_{ea}\frac{dT}{dz} \quad \text{for } z = 0$$
# $$\frac{dC_A}{dz} = \frac{dT}{dz} = 0 \quad \text{for } z = L$$
#
# These are the Danckwerts conditions of published `A2.1`, and its dict form is used
# unchanged. The notation list defines $\mathrm{Pe}_a = u_i d_p / D_{ea}$ on the
# *interstitial* velocity, so $\varepsilon D_{ea} = u_s d_p / \mathrm{Pe}_a$ and the
# dispersion coefficient of (11.6-1), divided by $u_s$, is simply the **length**
# $d_p/\mathrm{Pe}_{ma}$. Section 11.6 says $\mathrm{Pe}_a$ "may be considered to lie
# between 1 and 2", and of the thermal analogue that "little information is
# available".
#
# **Rung 3 - section 11.8.1, + interfacial gradients (`S7`).**
#
# $$-u_s\frac{dC}{dz} = k_g a_v (C - C_s^s) \tag{11.8.1-1}$$
# $$u_s\rho_g c_p\frac{dT}{dz} = h_f a_v (T_s^s - T) - 4\frac{U}{d_t}(T - T_r) \tag{11.8.1-2}$$
# $$\rho_B r_A = k_g a_v (C - C_s^s) \tag{11.8.1-3}$$
# $$(-\Delta H)\rho_B r_A = h_f a_v (T_s^s - T) \tag{11.8.1-4}$$
#
# with $r_A$ evaluated at $(C_s^s, T_s^s)$, and $C = C_0$, $T = T_0$ at $z = 0$.
# Section 11.8.1 states that it carries "the restrictions already mentioned in
# Section 11.5.1 for the basic case", which includes constant $u_s$, $\rho_g$, $c_p$.
#
# **Rung 4 - section 11.9.1, + intraparticle gradients (`S8`).** (11.9.1-1) and
# (11.9.1-2) repeat (11.8.1-1) and (11.8.1-2); the particle adds
#
# $$\frac{D_e}{\xi'^2}\frac{d}{d\xi'}\!\left(\xi'^2\frac{dC_s}{d\xi'}\right) - \rho_s r_A(C_s, T_s) = 0 \tag{11.9.1-3}$$
# $$\frac{\lambda_e}{\xi'^2}\frac{d}{d\xi'}\!\left(\xi'^2\frac{dT_s}{d\xi'}\right) + \rho_s(-\Delta H) r_A(C_s, T_s) = 0 \tag{11.9.1-4}$$
#
# $$\frac{dC_s}{d\xi'} = \frac{dT_s}{d\xi'} = 0 \ \text{ at } \xi' = 0; \qquad
#   k_g(C_s^s - C) = -D_e\frac{dC_s}{d\xi'}, \quad
#   h_f(T_s^s - T) = -\lambda_e\frac{dT_s}{d\xi'} \ \text{ at } \xi' = \tfrac{d_p}{2}$$
#
# The book's own footnote to (11.9.1-4) is worth quoting because it removes the
# commonest sign trap in this equation: *"The signs in (11.9.1-3) and (11.9.1-4) are
# those obtained when the rate is defined for first order as $r_A = kC_A$. Many books
# define $r_A = -kC_A$."* It then offers the reduction this page uses for the ladder,
# $\eta$ multiplying the surface rate,
#
# $$k_g a_v (C - C_s^s) = \eta\,\rho_B r_A(C_s^s, T_s^s) \tag{11.9.1-9}$$
# $$h_f a_v (T_s^s - T) = \eta\,\rho_B(-\Delta H) r_A(C_s^s, T_s^s) \tag{11.9.1-10}$$
# $$\eta = \frac{3}{\phi''^2}\left[\phi''\coth\phi'' - 1\right],\qquad
#   \phi'' = \frac{d_p}{2}\sqrt{\frac{k(T_s^s)\rho_s}{D_e}} \tag{11.9.1-11}$$
#
# **Rung 5 - section 11.7.2, two-dimensional pseudohomogeneous (`S6`).**
#
# $$(D_{er})_s\!\left(\frac{\partial^2 C}{\partial r^2} + \frac{1}{r}\frac{\partial C}{\partial r}\right) - u_s\frac{\partial C}{\partial z} - \rho_B r_A = 0 \tag{11.7.2-1}$$
# $$\lambda_{er}\!\left(\frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r}\right) - u_s\rho_g c_p\frac{\partial T}{\partial z} + \rho_B(-\Delta H) r_A = 0$$
#
# $$C = C_0,\ T = T_0 \ \text{ at } z=0; \qquad
#   \frac{\partial C}{\partial r} = 0 \ \text{ at } r = 0 \text{ and } r = R_t; \qquad
#   \frac{\partial T}{\partial r} = 0 \ \text{ at } r = 0$$
# $$\frac{\partial T}{\partial r} = -\frac{\alpha_w}{\lambda_{er}}(T_R - T_w) \ \text{ at } r = R_t$$
#
# There is **no** $4U/d_t$ term here: the wall coefficient $\alpha_w$ of eq.
# (11.7.1-1) replaces it, and eq. (11.7.4-1) is how the book converts between them,
#
# $$\frac{1}{\alpha_i} = \frac{1}{\alpha_w} + \frac{R_t}{4\lambda_{er}} \tag{11.7.4-1}$$
#
# **A printed inconsistency in rung 5's own groups, reported and not repaired.**
# Book page 574 prints
# $a_1 = (D_{er})_s/(u_i d_p) = 1/\mathrm{Pe}_{mr}$, and book page 572 prints, two
# pages earlier, that $(D_{er})_s$ "equals $\varepsilon D_{er}$, with $D_{er}$ based
# upon the interstitial flow velocity $u_i$". Those two statements are inconsistent
# by a factor $\varepsilon$: with $(D_{er})_s = \varepsilon D_{er}$ and
# $\mathrm{Pe}_{mr} = u_i d_p / D_{er}$,
# $(D_{er})_s/(u_i d_p) = \varepsilon/\mathrm{Pe}_{mr}$, while
# $(D_{er})_s/(u_s d_p) = 1/\mathrm{Pe}_{mr}$ exactly. The denominator should be
# $u_s$, matching $a_2 = \lambda_{er}/(G c_p d_p)$ which is built on $G = \rho_g u_s$
# in the very next line. Taken literally the printed group scales radial mass
# dispersion by $1/\varepsilon$ on this reactor (printed in the crossover cell). The
# page uses
# $(D_{er})_s/u_s = d_p/\mathrm{Pe}_{mr}$, states that as the repair-by-inference it
# is, and measures what the literal reading would cost - which is where the book's
# own sentence that "the computed results are not very sensitive with respect to
# $\mathrm{Pe}_{mr}$" gets tested rather than quoted.

# %% [markdown]
# ## Parameters and assumptions
#
# Section 11.5.2's parameters are **Van Welsenaere & Froment's (1970), restated in
# SI**, and the restatement is a unit conversion. Published `D2.2` carries their
# printed originals, so the conversion can be audited row by row rather than
# assumed - which is what the next cell does. Two of its findings decide numbers the
# rest of the page depends on:
#
# * **the tube radius.** Section 11.5.2 says the internal diameter is 2.54 cm, but
#   the printed $U = 0.096$ kJ/(m$^2$ s K) is eq. (11.7.4-1) evaluated at
#   $R_t = 0.0125$ m, and Example 11.5.3.A uses $R_t = 0.0125$ m throughout. At
#   0.0127 m the hot spot at the book's own $p_{A0} = 0.018$ is far outside the
#   operating range the same page states; at 0.0125 m it is just inside. Both
#   numbers are computed and printed in *Validation* below - none is typed here.
#   Arithmetic decides the radius, not the printed diameter.
# * **the specific heat.** 0.992 kJ/(kg K) is printed. Groups later in the same
#   chapter require more. The page uses the volumetric
#   $\rho_g c_p = 0.323$ kcal/(m$^3$ K) that they agree on, and reports the
#   printed value as a defect rather than repairing it. The deviation is computed
#   below.
#
# **Which pressure unit this page is in, and why it is not the one printed.**
# Section 11.5.2 prints *"the total pressure ... equal to 1 **bar**"* and
# *"$p_B^0 = 0.211$ **bar**"*, and its design point as *"0.018 **bar**"*. The model
# this page runs is in **atm**, and that is a deliberate reading with four reasons,
# each of which the notebook either prints or transcribes:
#
# 1. the rate constant. $\ln k = 19.837 - 13{,}636/T$ is Van Welsenaere & Froment's,
#    printed by them for $r_A = k\,p_B^0 p_A$ with the pressures in **atm**. Section
#    11.5.2 restates it *unchanged* while relabelling the pressures bar. Converting
#    the pressures without converting $k$ is a partial conversion, and the
#    self-consistent model is the unconverted one;
# 2. the book's own figure. Fig. 11.5.2-1's ordinate label reads **"$p_0$, atm"**
#    (read on a 300 ppi crop of book page 512; that unit word is the only thing this
#    page takes from any figure, and no curve, tick or point is read from any of
#    them). Page 513's own text ties the design pressures to that figure - it
#    discusses $p_{A0}$ = 0.0181, 0.0182 and 0.019 as "the upper part of the curves
#    ... (Figs. 11.5.2-1 and 11.5.2-2)";
# 3. Example 11.5.3.A's group $A$, printed in the same chapter, uses **0.208** -
#    the same oxygen pressure in atm, unconverted - where 11.5.2 prints 0.211 bar,
#    and the audit cell below shows 0.211 is 0.208 converted;
# 4. it is measurable, and the measurement is one-sided. Both readings are solved
#    in *Validation* below by the same grid-free Radau integration with one constant
#    changed: in atm the sensitivity sentence reproduces and its design point
#    sits just inside the ceiling the facing page states, while in the printed bar
#    the same sentence misses by tens of times more and its design point lands
#    hundreds of kelvin **above** that ceiling. Every one of those numbers is printed
#    by the code, none is typed here.
#
# So **every pressure this page reports is in atm**, including $p_{A0}^*$, the
# figure axes and section 11.7.3's 0.00924 - which book page 574 prints as a
# *mole fraction*, "corresponding to 44 g/m$^3$", and which is a partial pressure of
# 0.00924 atm only because $p_t = 1$ atm. It is Froment (1967)'s $N_{A0}$, and `C2.10`
# has found that the 44 is g per *normal* m$^3$; the cell below prints its value
# beside the printed one. The book's own strings are quoted verbatim with the unit it
# prints; nothing in the CSV is repaired.
#
# **What the chapter does not print, per rung.** This is the asymmetry the ladder
# turns out to have, and it is stated here before any of it is used:
#
# | rung | new coefficients it needs | printed for *this* reactor? |
# |---|---|---|
# | 1 `S2` | $U$ | **yes**, 0.096 kJ/(m$^2$ s K) - and reproduced by eq. (11.7.4-1) from rung 5's pair |
# | 2 `S4` | $D_{ea}$, $\lambda_{ea}$ | $\mathrm{Pe}_a$ "between 1 and 2"; of $\lambda_{ea}$, "little information is available", and $\mathrm{Pe}_{ha}$ is used in (11.6-3)/(11.6-4) without being defined anywhere in the book |
# | 3 `S7` | $k_g$, $h_f$ | **no** - the reader is sent to Chapter 3's correlation charts |
# | 4 `S8` | $D_e$, $\lambda_e$ | **no** - the reader is sent to Satterfield [1970] and Weisz & Hicks [1962], neither on disk |
# | 5 `S6` | $\lambda_{er}$, $\alpha_w$ | **yes**, both, to two figures |
#
# So the two ends of the ladder are pinned by the book and the middle is not.

# %%
# ---- Van Welsenaere & Froment's printed parameters, from D2.2's CSV -----------
V = dict(zip(VWF_PAR.index, VWF_PAR.value))
P = dict(M=V["M"], p_t=V["P"], rho_b=V["rho_b"], rho_g=V["rho_g"],
         minus_dH=V["minus_dH"], cp_vol=V["c_p"], u=V["u"], U=V["U"], R=V["R"],
         p_B0=V["p_B0"], a=V["a"], b=V["b"])
# c_p is ALREADY VOLUMETRIC, kcal/(m3 degC) - that trap is D2.2's finding, and
# D2.2's own unit column is what says so.  It is not rediscovered here.
assert VWF_PAR.loc["c_p", "unit"] == "kcal/(m3 degC)", VWF_PAR.loc["c_p", "unit"]

L_BED = bk("tube_length")
D_T_PRINTED = bk("tube_id")
D_P = bk("particle_diameter")
T_R0 = bk("coolant_temperature")                   # 352 C
T_PERM = bk("temp_range_high")                     # 415 C, the printed ceiling
R_T = P["R"]                                       # 0.0125 m; justified below
EPS = bk("voidage_a") + bk("voidage_b") * (
    1.0 + (D_T_PRINTED / D_P - 2.0) ** 2 / (D_T_PRINTED / D_P) ** 2)   # (11.5.1-16)
LAM_ER, ALPHA_W = bk("lam_er"), bk("alpha_w")
PE_MR, PE_HR = bk("pe_mr"), bk("pe_hr")

# the design case of model.py: every rung, pellet and reference solve reads it
RX = Reactor(**P, length=L_BED, d_p=D_P, T_r=T_R0, T_perm=T_PERM, eps=EPS,
             lam_er=LAM_ER, alpha_w=ALPHA_W, pe_mr=PE_MR, pe_hr=PE_HR,
             re_film=bk("Re_1173"))
A_G, B_G = RX.A, RX.B                              # the groups of (11.5.2-1/2)
THETA = RX.theta                                   # kmol/(m3 atm), from (11.5.2-1)
CP_MASS = RX.cp_mass                               # kJ/(kg K)
RHO_S = RX.rho_s
A_V = 6.0 * (1.0 - EPS) / D_P
C_G, kexp = RX.C, RX.kexp


# ================= the atm -> bar audit of section 11.5.2 =====================
ATM = 1.01325
rows = []


def add(what, book_key, orig, factor, note):
    b_ = bk(book_key)
    rows.append(dict(quantity=what, VWF_1970=orig, book_11_5_2=b_,
                     ratio=b_ / orig, vs_atm_to_bar=b_ / (orig * factor) - 1.0,
                     reading=note))


add("oxygen partial pressure p_B0", "oxygen_pressure", V["p_B0"], ATM, "converted")
add("heat of reaction (-dH)", "heat_of_reaction", V["minus_dH"] * KCAL, 1.0, "kcal -> kJ")
for k, q in (("vwf_p0_lower_625", "1a p0_lower"), ("vwf_p0_upper_625", "1a p0_upper"),
             ("vwf_p0_mean_625", "1a p0_mean"), ("vwf_p0_backint_625", "1a p0_critical")):
    ex = VWF_EX[(VWF_EX.example == "1a") & (VWF_EX.quantity == q.split()[1])]
    o = float(ex.value.iloc[0])
    add(f"Table 11.5.3.A-2 {q.split()[1]}", k, o, ATM,
        "converted" if "critical" not in q else "NOT converted")
AUDIT = pd.DataFrame(rows)
show(AUDIT.style.set_uuid("d11audit").hide(axis="index").format(
    {"VWF_1970": "{:.6g}", "book_11_5_2": "{:.6g}", "ratio": "{:.6f}",
     "vs_atm_to_bar": "{:+.2e}"}))

CONV_WORST = float(np.abs(AUDIT.vs_atm_to_bar.iloc[:-1]).max())
BACKINT_RESIDUAL = float(AUDIT.vs_atm_to_bar.iloc[-1])
print(f"\nThe first five rows are Van Welsenaere & Froment's atm values multiplied by"
      f"\n1 atm / 1 bar = {ATM}, reproduced to {CONV_WORST:.1e} relative - so section 11.5.2"
      f"\nIS their parameter set converted, and the conversion is not in doubt.")
print(f"The last row is not: the back-integrated critical inlet pressure is printed as"
      f"\n'{bk('vwf_p0_backint_625','as_printed')} bar' and is Van Welsenaere & Froment's"
      f" {VWF_EX[(VWF_EX.example=='1a')&(VWF_EX.quantity=='p0_critical')].value.iloc[0]} ATM"
      f"\nleft unconverted ({BACKINT_RESIDUAL:+.2e} against the converted value).")
_mean_bar = bk("vwf_p0_mean_625")
_bi_atm = float(VWF_EX[(VWF_EX.example == "1a") & (VWF_EX.quantity == "p0_critical")].value.iloc[0])
AGREE_AS_PRINTED = _mean_bar / _bi_atm - 1.0
AGREE_CONSISTENT = _mean_bar / (_bi_atm * ATM) - 1.0
print(f"\nThe book calls the two 'in excellent agreement'. As printed they differ by"
      f" {100*AGREE_AS_PRINTED:+.2f} %;"
      f"\nconverted consistently, by {100*AGREE_CONSISTENT:+.2f} % - the claim is BETTER than"
      f" the numbers beside it.")

# =================== what the chapter's own groups require ====================
route = []
route.append(("group B of Example 11.5.3.A = 257e6",
              bk("vwf_B"), P["minus_dH"] * P["rho_b"] * P["p_B0"] / P["cp_vol"],
              "p_B0 rho_b (-dH) / (rho_g c_p)"))
route.append(("group A of Example 11.5.3.A = 6150",
              bk("vwf_A"), A_G, "M p_t p_B0 rho_b / rho_g, the formula printed above it"))
_dTad = float(VWF_EX[(VWF_EX.example == "3") & (VWF_EX.quantity == "dT_ad")].value.iloc[0])
route.append(("A implied by dT_ad = 312.6 K at p0 = 0.0075 atm",
              bk("vwf_A"), bk("vwf_B") / (_dTad / 0.0075), "B / (dT_ad / p0)"))
GRP = pd.DataFrame(route, columns=["the book prints", "printed", "recomputed", "from"])
GRP["rel. dev."] = GRP.recomputed / GRP.printed - 1.0
show(GRP.style.set_uuid("d11groups").hide(axis="index").format(
    {"printed": "{:.6g}", "recomputed": "{:.7g}", "rel. dev.": "{:+.3e}"}))
A_PRINTED_DEV = float(GRP["rel. dev."].iloc[1])
B_PRINTED_DEV = float(GRP["rel. dev."].iloc[0])
print(f"\nB is reproduced to {abs(GRP['rel. dev.'].iloc[0]):.1e}, which pins (-dH), rho_b, p_B0"
      f"\nand the volumetric c_p together. A is not: the printed 6150 is"
      f" {100*A_PRINTED_DEV:+.3f} % from"
      f"\nits own printed formula, and the example's own dT_ad = 312.6 K at 0.0075 needs"
      f" {GRP.recomputed.iloc[2]:.1f}.")

# c_p: three independent routes against the printed 0.992 kJ/(kg K)
G_FLUX = P["rho_g"] * P["u"] / 3600.0                       # kg/(m2 s)
cp_routes = pd.DataFrame([
    ("volumetric c_p of D2.2's CSV / rho_g", P["cp_vol"] * KCAL / P["rho_g"]),
    ("B = 257e6 with (-dH), rho_b, p_B0 as printed",
     P["minus_dH"] * P["rho_b"] * P["p_B0"] / bk("vwf_B") * KCAL / P["rho_g"]),
    ("Pe_hr = 5.25 with lambda_er = 0.78e-3, G, d_p",
     PE_HR * LAM_ER / (G_FLUX * D_P)),
    ("Froment (1967) p. 19, loaded from C2.10", F67_PAR.loc["c_p", "value"] * KCAL),
], columns=["route", "c_p [kJ/(kg K)]"])
cp_routes["vs printed 0.992"] = cp_routes["c_p [kJ/(kg K)]"] / bk("specific_heat") - 1.0
show(cp_routes.style.set_uuid("d11cp").hide(axis="index").format(
    {"c_p [kJ/(kg K)]": "{:.6f}", "vs printed 0.992": "{:+.4f}"}))
CP_SPREAD = float(cp_routes["c_p [kJ/(kg K)]"].max() / cp_routes["c_p [kJ/(kg K)]"].min() - 1.0)
CP_PRINTED_DEV = float(cp_routes["vs printed 0.992"].mean())
print(f"\nFour routes, three of them inside Chapter 11 and the fourth Froment's own 1967"
      f"\npaper loaded from C2.10: they span {CP_SPREAD:.2e} and sit"
      f" {100*CP_PRINTED_DEV:+.1f} % from the 0.992 kJ/(kg K)"
      f"\nsection 11.5.2 prints. The printed value is recorded as a defect and is used"
      f" nowhere.")
print(f"\nBUT THE FOUR ARE NOT FOUR INDEPENDENT DETERMINATIONS, and the spread must not be"
      f"\nread as one. Route 2 reproduces route 1 to"
      f" {abs(cp_routes['c_p [kJ/(kg K)]'].iloc[1]/cp_routes['c_p [kJ/(kg K)]'].iloc[0]-1):.1e} BY CONSTRUCTION: B is printed as"
      f"\n(-dH) rho_b p_B0 / (rho_g c_p) and route 1 is the same c_p divided by the same"
      f" rho_g.\nRoute 4 is route 1 rounded -"
      f" {P['cp_vol']*KCAL/P['rho_g']/KCAL:.4f} kcal/(kg C) printed as"
      f" {F67_PAR.loc['c_p','value']} - and C2.10's own CSV note says its"
      f"\n0.25 is pinned on that page by the printed Pe_hR = 5.25, which is route 3. So the"
      f"\nfour trace to ONE printed determination plus route 3 evaluated exactly against the"
      f"\nsame number rounded. What survives is the CONCLUSION - the printed 0.992 is"
      f"\ninconsistent with the chapter's own groups - not an agreement among four witnesses."
      f"\nThat is why cp_route_spread carries a caveat in the break table.")

# the mean-density units error, provable from the same sentence
RHO_FROM_SENTENCE = bk("mass_flux") / bk("superficial_velocity")
print(f"\nSame sentence, same defect class: 'a mean fluid density of"
      f" {bk('mean_density','as_printed')} kg/m3'"
      f"\nwith G = {bk('mass_flux','as_printed')} kg/(m2 h) and u_s ="
      f" {bk('superficial_velocity','as_printed')} m/h, whose ratio is"
      f" {RHO_FROM_SENTENCE:.4f} kg/m3.\nD2.2's CSV gives Van Welsenaere & Froment's"
      f" {P['rho_g']} kg/m3; the printed value is out by 10^3.")

# ---- eq. (11.7.4-1): rung 1's U is rung 5's pair ----------------------------
ALPHA_I = 1.0 / (1.0 / ALPHA_W + R_T / (4.0 * LAM_ER))
U_VWF_SI = P["U"] * KCAL / 3600.0
ALPHA_I_PRINTED_R = 1.0 / (1.0 / ALPHA_W + (D_T_PRINTED / 2.0) / (4.0 * LAM_ER))
print(f"\neq. (11.7.4-1) at R_t = {R_T} m, alpha_w = {ALPHA_W}, lambda_er = {LAM_ER}:"
      f"\n  alpha_i = {ALPHA_I:.6f} kJ/(m2 s K)   against the printed"
      f" {bk('overall_U','as_printed')}, equal to round-off (within"
      f" 1e-12; here {ALPHA_I/bk('overall_U')-1:+.1e})"
      f"\n  Van Welsenaere & Froment's U = {P['U']} kcal/(m2 h K) ="
      f" {U_VWF_SI:.6f} kJ/(m2 s K)  -> {U_VWF_SI/ALPHA_I-1:+.3e}"
      f"\nAt R_t = {D_T_PRINTED/2} m (the printed 2.54 cm diameter) it would be"
      f" {ALPHA_I_PRINTED_R:.6f}, i.e. {ALPHA_I_PRINTED_R/bk('overall_U')-1:+.3e}"
      f"\nfrom the printed 0.096 - and it ROUNDS TO {ALPHA_I_PRINTED_R:.3f}, not to"
      f" {bk('overall_U'):.3f}.\nC2.10 recomputes U from this pair for the 1967 case, to"
      f" {100*C210_U_WORST:.2f} % at worst (its own claims);"
      f"\nthis is that check with the BOOK's SI values, and it is what fixes R_t = 0.0125 m."
      f"\nWHAT IT IS NOT: proof that the book DERIVED 0.096 from (11.7.4-1). VWF's own"
      f"\nU = {U_VWF_SI:.5f} rounds to the same 0.096 at the two significant figures printed,"
      f"\nso the identity discriminates the RADIUS, which is what it is used for here, and"
      f"\nnot the provenance of the number.")
U_FROM_2D_DEV = float(ALPHA_I / bk("overall_U") - 1.0)
U_PRINTED_R_DEV = float(ALPHA_I_PRINTED_R / bk("overall_U") - 1.0)

# ---- the unit the model is in, and the one the section prints ----------------
Y_A0 = bk("y_A0_1173")
M_IMPLIED = bk("y_A0_mass_1173") / 1000.0 / (Y_A0 * P["rho_g"] / P["M"])
print(f"\nUNITS. This page runs on Van Welsenaere & Froment's parameter set as D2.2 carries"
      f"\nit: p_t = {P['p_t']} ATM and p_B0 = {P['p_B0']} ATM, with ln k = {P['b']} - {P['a']:.0f}/T"
      f"\nas printed for those units. Section 11.5.2 prints"
      f" '{bk('total_pressure','as_printed')} bar' and"
      f"\n'{bk('oxygen_pressure','as_printed')} bar' for the same two quantities - and the audit"
      f" above shows the second\nIS the first converted - while restating k unchanged, which is"
      f" the partial\nconversion. Its own Fig. 11.5.2-1 labels the ordinate"
      f" '{bk('fig1152_1_ordinate_unit','as_printed')}' (book page"
      f" {int(BK.loc['fig1152_1_ordinate_unit','book_page'])};\nthat unit word is the only thing"
      f" this page takes from any figure in the chapter).\nEvery pressure reported below is"
      f" therefore in ATM. What the other reading costs is\nmeasured in Validation.")
print(f"\nAnd section 11.7.3's {Y_A0} is a MOLE FRACTION, not a pressure: page 574 prints it"
      f"\n'corresponding to {bk('y_A0_mass_1173','as_printed')} g/m3', which needs a molar mass of"
      f" {M_IMPLIED:.1f} kg/kmol - an aromatic C8,\nnot the {P['M']} kg/kmol mean. At p_t ="
      f" {P['p_t']} atm its partial pressure is numerically the\nsame number, which is why it can"
      f" be used as one.")
# the same 0.00924 is Froment (1967)'s N_A0, which C2.10 has already audited
NA0_F67 = float(F67_PAR.loc["N_A0", "value"])
assert NA0_F67 == Y_A0, "11.7.3's inlet mole fraction is not C2.10's N_A0"
NA0_C210 = NA0_F67 * (1.0 + C210_C["na0_from_44gnm3_rel"])     # C2.10's 44 g/Nm3 reading
C_NORMAL = 101325.0 / (8.314462618 * 273.15) / 1000.0         # kmol/m3 at 0 C and 1 atm
M_OXYLENE = 8 * 12.011 + 10 * 1.008                           # C8H10, kg/kmol
M_NORMAL = bk("y_A0_mass_1173") / 1000.0 / (Y_A0 * C_NORMAL)
print(f"\nC2.10 HAS READ THE SAME NUMBER, as Froment (1967)'s N_A0 = {NA0_F67}, and its CSV"
      f" records\nthat the 44 is g per NORMAL m3 (the paper's Figure 11 axis is g/Nm3):"
      f" 44 g/Nm3 of o-xylene\nis N_A0 = {NA0_C210:.5f}, {100*C210_C['na0_from_44gnm3_rel']:+.2f} %"
      f" from the printed {NA0_F67} (C2.10's claim).\nThe {M_IMPLIED:.1f} kg/kmol above divides by"
      f" rho_g/M = {P['rho_g']/P['M']:.5f} kmol/m3; the molar density at\n0 C and 1 atm is"
      f" {C_NORMAL:.5f}, and with it the molar mass is {M_NORMAL:.2f} kg/kmol, against"
      f" o-xylene's\n{M_OXYLENE:.2f}. The rest of the gap is the printed 0.00924 against"
      f" {NA0_C210:.5f}:\n{M_NORMAL:.2f} x {Y_A0}/{NA0_C210:.5f} = {M_NORMAL*Y_A0/NA0_C210:.2f}."
      f" The page uses the printed mole fraction, as the chapter does.")

print(f"\nderived: eps = {EPS:.6f} from (11.5.1-16) at d_t/d_p ="
      f" {D_T_PRINTED/D_P:.4f};  rho_s = {RHO_S:.1f} kg/m3;"
      f"  a_v = {A_V:.2f} m2/m3\n         theta = {THETA:.6f} kmol/(m3 atm);"
      f"  A = {A_G:.4f};  B = {B_G:.6e};  C = {C_G():.2f} 1/h")


# %% [markdown]
# ## The data
#
# There is no measurement anywhere on this page. Chapter 11 reports none, and the
# sources it draws on report none either, so the provenance tier is 6 throughout.
#
# **This page's own dataset.** `froment-2011-ch11-printed.csv` is every model
# constant, stated result and quantitative claim read out of Chapter 11 (the first
# code cell prints how many rows), each with the book page it came from and a `kind`
# of `parameter`, `stated_result` or `claim`. Its sidecar records that several rows are transcribed *because they are
# wrong* and lists them; none is repaired anywhere. `as_printed` keeps the character
# string, thousands separators and all, so that a defect can be quoted verbatim.
#
# **Six borrowed datasets from four pages, and what each source page already
# establishes about the rows used here.** The first code cell lists the six with
# their row counts; the cell below prints every number quoted here from the source
# page's own committed claims, so none is retyped.
#
# * `D2.2`, `van-welsenaere-froment-1970-parameters.csv` and `-examples.csv`.
#   `D2.2` establishes: **the printed `c_p` is already volumetric**, kcal/(m$^3$
#   $^\circ$C), so multiplying it by $\rho_g$ would shift the wall group by the gas
#   density - that finding is the reason the audit above can be run at all, and it
#   is asserted against `D2.2`'s own unit column rather than restated. `D2.2` also
#   establishes that the Elsevier full-text API drops the 1970 paper's mid-dot
#   decimal separators, how wide its criterion brackets are, and how much hot-spot
#   error its own first-order upwind leaves at $n_z = 1500$ - the last of which
#   **does** affect this page, and is why every number reported here is
#   Richardson-extrapolated. Its Section-6 values are used here only as the
#   *originals* against which the book's Example 11.5.3.A is checked; none of its
#   criteria is rebuilt.
# * `C2.10`, `froment-1967-parameters.csv` and `froment-1967-stated-results.csv`.
#   `C2.10` establishes: the 1967 paper's own gas constant is **1.98** cal/(mol K),
#   not 1.987, which is what makes 27,000/1.98 = 13,636 the number section 11.5.2
#   prints; that $G$ = 4684 is printed "4.684" on the paper's p. 24; that
#   $-\Delta H_3$ is printed "1.090"; where the 1-D and 2-D runaway limits of the
#   1967 case lie, located by bisection; and how closely $U$ recomputed from
#   $\lambda_R$ and $\alpha_w$ agrees with the printed one; and that the inlet mole
#   fraction $N_{A0}$ = 0.00924, the same number section 11.7.3 prints, is stated to
#   correspond to 44 g per *normal* m$^3$, which gives a different mole fraction
#   (affects this page only as a reading of the 44: the page uses the printed
#   0.00924, and *Parameters and assumptions* prints C2.10's value beside it). All
#   five matter here and all five are printed beside the book's restatement. `C2.10` also warns
#   that $U$ = 82.7 is a **shared source** with `D2.2`, not a second witness - so the
#   audit above never treats it as one.
# * `B1.1`, `isothermal-exact.csv` (`dataset_id: thiele-isothermal-exact`). `B1.1`
#   establishes the exact isothermal $\eta$ for slab, cylinder and sphere; how
#   closely its own pymrm pellet matches it; and - load-bearing here - that
#   **`NumJac` must be given `(n_r, 1)`, never a bare `(n_r,)`**. Its sphere rows are
#   the reference for eq. (11.9.1-11) below. Its fold at $\beta = 0.6$,
#   $\gamma = 20$ is *not* touched: the pellet here runs at a far smaller Prater
#   number, computed and printed.
# * `A1.1`, `ergun-1952-parameters.csv`. `A1.1` establishes that a refit of Ergun's
#   own markers gives constants different from 150 and 1.75, and what the line Ergun
#   labels Kozeny-Carman measures. Neither changes anything here, because Example
#   11.5.1.A is arithmetic *on the printed constants*; the refit is printed beside
#   them and is a break row. `A1.1` carries **no** McDonald and no Hicks constants and
#   does not discuss the equivalent-diameter definition, so those parts of Example
#   11.5.1.A are new here.
#
# No figure is digitised, no page image exists anywhere in this directory, and every
# figure-only content of Chapter 11 - Figs. 11.5.2-1/2, 11.6-1, 11.7.1-1/2/3, all six
# 11.7.3 figures, 11.8.1-1..4, 11.8.2-1 and the whole of Fig. 11.5.3-1's runaway band
# - is **out of scope** and never traced or read off. The one exception, and it is a
# unit rather than a datum: `fig1152_1_ordinate_unit` is the *label* on Fig. 11.5.2-1's
# ordinate, read on a 300 ppi crop as "$p_0$, atm". It decides which unit the whole
# page is in, so it is transcribed as a row like everything else and named here.

# %%
print("What the source pages establish about the rows used here, from their claims:")
print(f"   D2.2 : criterion-1 brackets {D22_C['bracket_width_c1_min_pct']:.0f} to"
      f" {D22_C['bracket_width_c1_max_pct']:.0f} % wide; first-order upwind leaves"
      f" {D22_C['hotspot_error_n1500_K']:.2f} K of hot-spot error at n_z = 1500")
print(f"   C2.10: runaway limits {C210_C['t0crit_1d_C']:.2f} C (1-D) and"
      f" {C210_C['t0crit_2d_C']:.2f} C (2-D); U from lambda_R and alpha_w within"
      f" {100*C210_U_WORST:.2f} %")
print(f"   B1.1 : its pymrm pellet against the exact table, worst"
      f" {B11_C['isothermal_max_rel_error_phi_le_30']:.1e} for phi <= 30 and"
      f" {B11_C['isothermal_max_rel_error_all_phi']:.1e} over all phi")
print(f"   A1.1 : refit k1 = {A11_C['k1_refit_all_points']:.1f},"
      f" k2 = {A11_C['k2_refit_all_points']:.3f}; Kozeny-Carman line"
      f" {bk('ergun_b')*A11_C['kozeny_carman_line_over_k1']:.1f}")

# %% [markdown]
# ## PyMRM implementation
#
# One class, five configurations: `Ladder` in this page's `model.py`, beside the
# pellet (`Pellet`), the threshold search (`p0_critical` with its
# `crossing_control`) and the grid-free rung-1 route (`ode_reference`). The state is
# `(n_z, n_f)` with the axial index first and fields last, except for rung 5 which
# is `(n_z, n_r, 2)` - spatial axes first, fields last, one layout throughout. Every
# constant operator is assembled once per solve in `_ops`, printed below, never
# inside a Newton iteration. Each solve is `pymrm.newton` with its default damping
# and direct sparse LU, on a Jacobian whose source part is `NumJac`'s finite
# differences.
#
# Three things are worth reading before the code:
#
# * **the divergence geometry differs between the two second axes.** Rung 5's radial
#   axis is a tube, so `construct_div(..., nu=1, axis=1)` is cylindrical; the pellet
#   of rung 4 is a sphere, so `nu=2`. Rung 1's and rung 2's axial axis is
#   `nu=0`, Cartesian. Getting these wrong is silent: the profile stays smooth.
# * **the two boundary conditions that carry physics** are written with the outward
#   normal spelt out in a comment, because the sign of `a` means opposite things at
#   the two ends. The Danckwerts inlet is `A2.1`'s published dict, scaled: the
#   dispersion coefficient divided by $u_s$ is a length, so `a` is
#   $d_p/\mathrm{Pe}$ and `b` is 1. The wall condition of rung 5 is
#   $\lambda_{er}\,\partial T/\partial n + \alpha_w T = \alpha_w T_w$ at $r = R_t$,
#   where the outward normal is $+r$.
# * **the film balances are scaled to the units of their own unknown.** Written as
#   the book writes them, (11.8.1-3) is in atm/h and (11.8.1-4) in
#   kcal/(m$^3$ h), rows of very different size, so the absolute residual test in
#   `Ladder.solve` would mean something different in each. Dividing each by its own
#   transfer coefficient puts the two algebraic rows in atm and K, the units of the
#   surface unknowns they determine.

# %%
# The five rungs are model.py's `Ladder`; its operator assembly, printed here, is
# where the rungs differ.  The bisection, the crossing control and the
# grid-free rung-1 route are `p0_critical`, `crossing_control` and
# `ode_reference` in the same file.
print(inspect.getsource(Ladder._ops))

# %%
P0_BOOK = bk("runaway_p0_safe")   # 0.018, section 11.5.2's own design point, in atm
P0_REF = bk("y_A0_1173")          # 0.00924, section 11.7.3's own inlet mole fraction
#                                 # for the SAME reactor, hence its partial pressure
#                                 # in atm at p_t = 1 atm.
NS_RUNG1 = (150, 300, 600, 1200, 2400, 4800)     # the ONE rung-1 refinement chain
_RUNG1_CACHE = {}
# The order fit is only meaningful if the finest grid's error is far above the
# thresholds' own noise floor (claims-design section 3).  That floor is NOT the 1e-12
# bisection bracket: a threshold moves by about 1e-8 relative when only its starting
# guess changes, so 1e-8 is the floor used here.
START_FLOOR = 1e-8


def rung1_chain(R=None, ceiling=None, ns=NS_RUNG1):
    """Rung 1's threshold and hot spot on a refinement chain, its observed orders,
    the Richardson limits, and the same two numbers from the grid-free ODE route.

    THIS IS THE ONLY PLACE rung 1's reported numbers are computed.  The Validation
    cell below DISPLAYS this function's return value, `h_rung1` RETURNS its metric
    subset, and the break rows call it with a perturbed R or ceiling - so the
    displayed number and the agreement.json number are the same object, not two
    computations that happen to agree.

    Memoised on its own arguments only, which is a pure function of them, so it
    changes nothing about determinism.
    """
    R = R_T if R is None else R
    ceiling = T_PERM if ceiling is None else ceiling
    key = (R, ceiling, tuple(ns))
    if key in _RUNG1_CACHE:
        return _RUNG1_CACHE[key]
    rows = []
    for n in ns:
        lad = Ladder(RX, rung=1, n_z=n, R=R)
        Tm, zm, _edge = hotspot(lad.z_c, lad.solve(P0_BOOK)[:, 1])
        rows.append(dict(n_z=n, T_hot=Tm, z_hot=zm,
                         p0_star=p0_critical(lad, ceiling=ceiling)))
    tab = pd.DataFrame(rows)
    ordT, limT = richardson(tab.T_hot, tab.n_z)
    ordP, limP = richardson(tab.p0_star, tab.n_z)
    ordZ, limZ = richardson(tab.z_hot, tab.n_z)
    T_ode, z_ode, _ = ode_reference(RX, P0_BOOK, R=R)
    p0_ode = p0_critical_ode(RX, R=R, ceiling=ceiling)
    out = dict(table=tab, order_z=ordZ, limit_z=limZ, limit_T=limT, limit_p0=limP,
               T_ode=T_ode, z_ode=z_ode, p0_ode=p0_ode,
               metrics=dict(rung1_p0_star_extrapolated=limP, rung1_order_p0=ordP,
                            rung1_order_T=ordT,
                            ode_vs_pymrm_p0=abs(limP / p0_ode - 1.0),
                            ode_vs_pymrm_T=abs(limT / T_ode - 1.0),
                            p0_star_vs_printed_0018=p0_ode / P0_BOOK - 1.0,
                            rung1_finest_err_over_start_floor=abs(
                                float(tab.p0_star.iloc[-1]) / limP - 1.0) / START_FLOOR))
    _RUNG1_CACHE[key] = out
    return out


_CROSS_R_CACHE = {}


def shift_rung5(R, n_z=600, n_r=10):
    """rung 5 against rung 1 AT THE SAME RADIUS - both thresholds move with R_t."""
    return 100.0 * (p0_critical(Ladder(RX, rung=5, n_z=n_z, n_r=n_r, R=R), hi=0.09)
                    / p0_critical(Ladder(RX, rung=1, n_z=n_z, R=R), hi=0.09) - 1.0)


R_SWEEP = (R_T, 0.0120, 0.0115, 0.0110, 0.0105, 0.0100, 0.0098, 0.0096, 0.0094,
           0.0090, 0.0085, 2.5 * D_P)
#          ^ from the chapter's own tube outward, down to d_t/d_p = 5, which is
#            where section 11.5.1 says the Ergun constants "drastically change".
#            The sweep runs until the target is straddled by two VALIDATED points
#            and stops there - each evaluation is two threshold bisections, one of
#            them two-dimensional, so the points past the bracket are not spent.
#            What happens on narrower tubes is measured separately, below the
#            crossover table, where it is the question being asked.


def cross_R(target, n_z=600, n_r=10, iters=4):
    """The tube radius at which rung 5's shift crosses `target`, from the same
    controlled scan as every other row.  Memoised on its own arguments so that the
    display cell and the break-table helper are ONE computation: each evaluation is
    two full threshold bisections, one of them two-dimensional, so this is the most
    expensive search on the page.

    ITERS IS A BRACKET TOLERANCE AND NOT AN ACCURACY: four refinement levels of the 0.5 mm straddling pair
    is about 2e-6 m, while the GRID moves the shift by about 0.1 percentage points,
    which on the local slope of the sweep is some 3e-5 m of R_t*.  The display cell
    computes that from the page's own two grids and prints R_t* with it."""
    key = (target, n_z, n_r, iters)
    if key not in _CROSS_R_CACHE:
        _CROSS_R_CACHE[key] = cross_scan(
            lambda R: shift_rung5(R, n_z=n_z, n_r=n_r), R_SWEEP, target,
            levels=iters)
    return _CROSS_R_CACHE[key]


def unit_reading(p_B0, p_t, R=R_T, ceiling=T_PERM):
    """Rung 1's threshold and its hot spot at the printed 0.018, on a stated
    reading of section 11.5.2's pressure unit: the design case with p_B0 and p_t
    replaced.  Grid-free (Radau) throughout, so the comparison between the two
    readings carries no discretisation at all."""
    rx = replace(RX, p_B0=p_B0, p_t=p_t)
    return dict(p0_star=p0_critical_ode(rx, R=R, ceiling=ceiling),
                T_at_book_p0=ode_reference(rx, P0_BOOK, R=R)[0])


G_FLUX = RX.mass_flux
FILM = film(RX)
print("film coefficients from Wakao & Funazkri (page A3.4), at the chapter's printed Re:")
for k in ("Re", "Sc", "Pr", "Sh", "Nu"):
    print(f"   {k:>3s} = {FILM[k]:10.4f}")
print(f"   k_g = {FILM['k_g']:.6f} m/s      h_f = {FILM['h_f']:.6f} kJ/(m2 s K)"
      f"      a_v = {FILM['a_v']:.2f} m2/m3")


# %% [markdown]
# ## Results
#
# The comparison metric is **the inlet partial pressure at which the hot spot
# reaches 415 $^\circ$C**, root-found on each rung. 415 $^\circ$C is not a choice:
# section 11.5.2 states the operating temperature range of this catalyst as
# "335$^\circ$ to 415$^\circ$C", so the ceiling is the book's own, and the quantity is
# a threshold rather than a sampled maximum. Call it $p_{A0}^*$.
#
# Reported alongside are the hot spot and its location at the book's own
# $p_{A0} = 0.018$ atm, so that both a fixed-input and a fixed-output comparison are
# on the page.

# %%
N_Z, N_R = 1200, 20
RUNGS = [(1, "11.5", "S2", "basic 1-D plug flow"),
         (2, "11.6", "S4", "+ axial dispersion"),
         (3, "11.8", "S7", "+ film resistance"),
         (4, "11.9", "S8", "+ intraparticle profiles"),
         (5, "11.7", "S6", "2-D radial dispersion")]
HI = {1: 0.0200, 2: 0.0200, 3: 0.0200, 4: 0.0600, 5: 0.0200}
# The profiles are compared at P0_REF rather than at 0.018 for a reason that is itself
# a result: 0.018 is above the threshold of some of the rungs (the count is COMPUTED and
# printed below, not asserted here), so their "hot spot" there is a runaway state rather
# than an operating point.  P0_REF is subcritical on all five, and the chapter prints it
# for this reactor.

SOL, LAD, ROWS = {}, {}, []
for r, sec, st, what in RUNGS:
    lad = Ladder(RX, rung=r, n_z=N_Z, n_r=N_R)
    LAD[r] = lad
    c = solve_robust(lad, P0_REF)
    SOL[r] = c
    Tm, zm, edge = hotspot(lad.z_c, lad.bulk_T(c))
    assert not edge, f"rung {r}: hot spot at a grid edge, not an interior root"
    pc = p0_critical(lad, hi=HI[r])
    ROWS.append(dict(rung=r, section=sec, structure=st, adds=what,
                     T_hot=Tm - 273.15, z_hot=zm,
                     p_exit=float(c[..., 0].mean(axis=-1)[-1] if r == 5 else c[-1, 0]),
                     p0_star=pc, newton_residual=lad.last_residual))
LADDER = pd.DataFrame(ROWS)
LADDER["vs_rung1_pct"] = 100.0 * (LADDER.p0_star / LADDER.p0_star.iloc[0] - 1.0)
# The hot spot at 0.018 atm is reported only for rungs still BELOW their threshold
# there.  Past it, the fixed-inlet solve starts cold from the rung-1 profile, which is
# then a runaway guess, and which state Newton reaches (or whether it converges) is
# decided by round-off (see the note on locating thresholds in Validation); branch
# following settles this for thresholds but not for a fixed-inlet solve above one.
LADDER["T_hot_at_0018"] = [(f"{T_max_safe(LAD[r], P0_BOOK) - 273.15:.1f}"
                            if LADDER.p0_star[i] > P0_BOOK else "past threshold")
                           for i, (r, *_) in enumerate(RUNGS)]
assert (LADDER.newton_residual < 1e-8).all(), "a ladder solve did not reach 1e-8"
show(LADDER.style.set_uuid("d11ladder").hide(axis="index").format(
    {"T_hot": "{:.2f}", "z_hot": "{:.4f}", "p_exit": "{:.3e}",
     "p0_star": "{:.7f}", "vs_rung1_pct": "{:+.2f}",
     "newton_residual": lambda v: "< 1e-8"}))

P0S = dict(zip(LADDER.rung, LADDER.p0_star))
SHIFT = dict(zip(LADDER.rung, LADDER.vs_rung1_pct))
R34_LIKE = 100.0 * (P0S[4] / P0S[3] - 1.0)      # rung 3 -> rung 4, against RUNG 3
R34_RATIO = abs(R34_LIKE) / abs(SHIFT[3])
N_SUPERCRIT = int((LADDER.p0_star < P0_BOOK).sum())
print(f"\nAll pressures are in ATM (see Parameters). T_hot and z_hot are at"
      f" p_A0 = {P0_REF} atm\n(section 11.7.3's own mole fraction for this reactor, at"
      f" p_t = {P['p_t']} atm); T_hot_at_0018 is\nat section 11.5.2's design point."
      f"\n\n'past threshold' MARKS A RUNG WHOSE OWN p0* IS BELOW 0.018 atm. Its hot spot"
      f" there is a\nrunaway state reached by a cold-started fixed-inlet solve, and which"
      f" state that\nsolve reaches is decided by round-off (see the note on locating"
      f" thresholds in Validation), so\nno temperature is printed for it. Newton residuals are all below"
      f" 1e-8 (asserted)."
      f"\n\nAT 0.018 ATM, {N_SUPERCRIT} OF THE FIVE RUNGS ARE ALREADY PAST THEIR OWN"
      f" THRESHOLD - rungs"
      f" {', '.join(str(r) for r in LADDER.rung[LADDER.p0_star < P0_BOOK])} - which is why"
      f" the\nprofile comparison is made at {P0_REF} atm and not there.")
print(f"\nEach rung, as a shift in the safe inlet partial pressure of the SAME reactor,"
      f"\nmeasured against RUNG 1:")
for r, sec, st, what in RUNGS[1:]:
    print(f"   {what:<26s} (sec. {sec}, {st})   {SHIFT[r]:+7.2f} %")
print(f"\nThe two heterogeneous rungs pull in OPPOSITE directions. Like for like - rung 4"
      f"\nagainst RUNG 3, which is the model it corrects, not against rung 1 -"
      f"\n   rung 1 -> rung 3 : {SHIFT[3]:+.2f} %   (of rung 1's threshold)"
      f"\n   rung 3 -> rung 4 : {R34_LIKE:+.2f} %   (of rung 3's threshold)"
      f"\nso the second correction is {R34_RATIO:.1f} times the size of the first. Against the"
      f" common\nrung-1 baseline the same pair reads {SHIFT[3]:+.2f} and {SHIFT[4]:+.2f} %, a"
      f" difference of {SHIFT[4]-SHIFT[3]:+.2f}\nPERCENTAGE POINTS - a different quantity, and"
      f" not a percentage of anything.\nEither way: stopping at rung 3 leaves you further from"
      f" rung 4"
      f" ({abs(P0S[4]-P0S[3]):.6f} atm)\nthan rung 1 was"
      f" ({abs(P0S[4]-P0S[1]):.6f} atm).")


# %%
fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.0))
COL = {1: C_BLUE, 2: C_GREEN, 3: C_ORANGE, 4: C_PURPLE, 5: C_RED}
for r, sec, st, what in RUNGS:
    lad, c = LAD[r], SOL[r]
    ax[0].plot(lad.z_c, lad.bulk_T(c) - 273.15, color=COL[r], lw=1.6,
               label=f"{r}. {sec} {st} {what}")
    ax[1].plot(lad.z_c, (c[..., 0].mean(axis=-1) if r == 5 else c[..., 0]),
               color=COL[r], lw=1.6)
ax[0].plot(LAD[5].z_c, LAD[5].axis_T(SOL[5]) - 273.15, color=C_RED, lw=1.0, ls=":",
           label="5. axis temperature")
ax[0].axhline(T_PERM - 273.15, color="0.35", lw=1.0, ls="--")
ax[0].annotate(f"{T_PERM-273.15:.0f} $^\\circ$C, top of the printed operating range",
               (1.05, T_PERM - 273.15 + 3), fontsize=8, color="0.25")
ax[0].set_xlabel("z [m]"); ax[0].set_ylabel("temperature [$^\\circ$C]")
ax[0].set_title(f"five models, one reactor, $p_{{A0}}$ = {P0_REF} atm")
ax[0].legend(fontsize=7.2, loc="upper right")
ax[1].set_xlabel("z [m]"); ax[1].set_ylabel("$p_A$ [atm]")
ax[1].set_title("hydrocarbon partial pressure")
fig.tight_layout(); plt.show()

fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.barh([f"{r}. sec. {s} ({t})" for r, s, t, _ in RUNGS], LADDER.vs_rung1_pct,
        color=[COL[r] for r, *_ in RUNGS], height=0.55)
ax.axvline(0.0, color="0.3", lw=1.0)
ax.set_xlabel("shift in the safe inlet partial pressure $p_{A0}^*$ [%] vs rung 1"
              "\n($p_{A0}^*$ in atm; see Parameters)")
ax.set_title("what each rung buys, on the same reactor")
ax.invert_yaxis(); fig.tight_layout(); plt.show()


# %% [markdown]
# ## Validation
#
# Ranked by strength, and the ranking is the honest one: the chapter reports no
# experiment, so nothing below is a validation against data.
#
# 1. **Example 11.5.1.A's pressure drop is fully printed arithmetic** - three
#    friction factors and three pressure drops from stated inputs. Reproducible to
#    the digit, and the strongest check available in the chapter. It is also where a
#    contradiction with eq. (11.5.1-9) shows up.
# 2. **Section 11.5.2's parametric-sensitivity sentence** is a root-findable
#    threshold with a printed ceiling, and it tests the whole of rung 1 - equations,
#    constants, wall group - in one number.
# 3. **Example 11.5.3.A** is Van Welsenaere & Froment's own worked example recopied,
#    so the book's transcription can be checked against `D2.2`'s CSV of the
#    originals. This checks the *transcription*, not the model - and, as the cell
#    below shows, a row-by-row transcription check is not enough on its own: it is
#    the **cross-column** test the table imposes on itself,
#    $\Delta T_{ad} = (B/A)\,p_0$, that finds the larger defect.
# 4. **Equation (11.9.1-11)** is an analytical solution of (11.9.1-3), so the pellet
#    grid solve can be checked against it and against `B1.1`'s exact table.
# 5. **Sections 11.8.2, 11.7.3 and 11.7.4 state results without arithmetic** - "of
#    the order of 1 $^\circ$C", "about 30", "less than 360", "365". Checkable to a
#    figure or two, no better.
# 6. **Grid convergence and a second solver family** are internal, not validation of
#    the book - but they are what makes the rung-to-rung differences meaningful,
#    because two of them are smaller than the discretisation error of a single
#    1200-cell solve.
#
# **Locating a threshold.** `p0_critical` makes one cold solve at the bottom of the
# bracket and starts every later solve from the converged state at the current lower
# end, so it follows the branch the reactor is on as the feed is raised. Cold-starting
# every solve from the rung-1 profile at the same $p_{A0}$ would not be neutral:
# between rung 1's threshold and a higher one that guess is already a runaway profile,
# and whether Newton then reaches the ignited state or fails (and the fallback ramp
# finds the lower state) is decided by round-off. Direct starts on rung 2 at
# $n_z = 900$ and $\mathrm{Pe}_a = 0.24$, at four inlet pressures $10^{-5}$ atm apart
# just below 0.0182 atm, reached the ignited state far above the ceiling, the lower
# state or no convergence, in no monotone order; nothing on this page recomputes
# them, because by this argument the outcome differs between machines. With the
# branch followed, the rung-2 crossover is 0.24883, inside the grid-set error bar the
# crossover cell prints, and every point of every scan passes the control, on three
# BLAS kernels.
# **What branch following does not cover:** a fixed-inlet solve *above* a threshold
# starts cold from the rung-1 profile. The page therefore prints no hot spot at
# 0.018 atm for a rung already past its threshold, states rung 5's coarse-grid state
# at 0.018 atm only as a bound, and does not say which rungs brentq survives.

# %%
# ============ 1. Example 11.5.1.A, and the two equivalent diameters ==========
E = dict(L=bk("ex1151a_length"), d=bk("ex1151a_cylinder"), eps=bk("ex1151a_voidage"),
         G=bk("ex1151a_mass_flux"), mu=bk("ex1151a_viscosity"),
         rho=bk("ex1151a_density"))
E["u"] = E["G"] / E["rho"]
ERG = dict(zip(ERG_PAR.quantity, ERG_PAR.value))
print("Ergun's constants, loaded from A1.1 rather than retyped:")
for q in ERG_PAR.quantity:
    print(f"   {q:<34s} {ERG[q]}")
print(f"   and A1.1's own refit of Ergun's markers gives"
      f" k1 = {A11_C['k1_refit_all_points']:.1f}, k2 = {A11_C['k2_refit_all_points']:.3f}"
      f" (its claims),\n   which is A1.1's finding and changes nothing here: Example 11.5.1.A is"
      f"\n   arithmetic ON the printed constants.")
assert ERG["k1_viscous_constant"] == bk("ergun_b") and ERG["k2_kinetic_constant"] == bk("ergun_a"), \
    "A1.1's Ergun constants and the book's printed a, b disagree"


def f_ergun(dp, a, b):
    Re = dp * E["G"] / E["mu"]
    return Re, (1.0 - E["eps"]) / E["eps"] ** 3 * (a + b * (1.0 - E["eps"]) / Re)


def f_hicks(dp):
    Re = dp * E["G"] / E["mu"]
    return Re, (bk("hicks_coefficient") * (1.0 - E["eps"]) ** bk("hicks_exponent_eps")
                / E["eps"] ** 3 * Re ** bk("hicks_exponent_Re"))


def dP(f, dp):
    return f * E["rho"] * E["u"] ** 2 / dp * E["L"]          # eq. (11.5.1-3)


DP_AREA = np.sqrt(1.5) * E["d"]      # sphere of equal surface AREA - what is printed
DP_SV = E["d"]                       # eq. (11.5.1-9), 6 V_p/S_p, for a cylinder d = h
DP_ROUNDED = bk("ex1151a_dp")        # the printed 0.0037, used by the book downstream
tab = []
for nm, dp in (("as printed, 0.0037 m", DP_ROUNDED),
               ("sqrt(1.5) x 0.003, unrounded", DP_AREA),
               ("eq. (11.5.1-9): 6 V/S = d", DP_SV)):
    Re, fe = f_ergun(dp, bk("ergun_a"), bk("ergun_b"))
    _, fm = f_ergun(dp, bk("mcdonald_a"), bk("mcdonald_b"))
    _, fh = f_hicks(dp)
    tab.append(dict(d_p=dp, basis=nm, Re=Re, f_Ergun=fe, dP_Ergun=dP(fe, dp),
                    f_McDonald=fm, dP_McDonald=dP(fm, dp),
                    f_Hicks=fh, dP_Hicks_bar=dP(fh, dp) / 1e5))
EX1 = pd.DataFrame(tab)
EX1.loc[len(EX1)] = dict(d_p=DP_ROUNDED, basis="THE BOOK PRINTS", Re=bk("ex1151a_Re"),
                         f_Ergun=bk("ex1151a_f_ergun"), dP_Ergun=bk("ex1151a_dp_ergun"),
                         f_McDonald=bk("ex1151a_f_mcdonald"),
                         dP_McDonald=bk("ex1151a_dp_mcdonald"),
                         f_Hicks=bk("ex1151a_f_hicks"), dP_Hicks_bar=bk("ex1151a_dp_hicks"))
show(EX1.style.set_uuid("d11ex1151a").hide(axis="index").format(
    {"d_p": "{:.6f}", "Re": "{:.2f}", "f_Ergun": "{:.4f}", "dP_Ergun": "{:.0f}",
     "f_McDonald": "{:.4f}", "dP_McDonald": "{:.0f}", "f_Hicks": "{:.4f}",
     "dP_Hicks_bar": "{:.4f}"}))

EX1151A_ERGUN_DEV = abs(EX1.dP_Ergun.iloc[0] / bk("ex1151a_dp_ergun") - 1.0)
EX1151A_MCD_DEV = abs(EX1.dP_McDonald.iloc[0] / bk("ex1151a_dp_mcdonald") - 1.0)
EX1151A_HICKS_DEV = abs(EX1.dP_Hicks_bar.iloc[0] / bk("ex1151a_dp_hicks") - 1.0)
EX1151A_WORST = max(EX1151A_ERGUN_DEV, EX1151A_MCD_DEV, EX1151A_HICKS_DEV)
DP_DEFINITION_COST = EX1.dP_Ergun.iloc[2] / EX1.dP_Ergun.iloc[0] - 1.0
print(f"\nOn the printed d_p = 0.0037 m the three pressure drops come back to"
      f" {EX1151A_WORST:.1e}\nrelative, i.e. to every digit the example prints."
      f"\n\nBut the printed d_p is NOT the chapter's own definition. Eq. (11.5.1-9) reads"
      f"\nd_p = 6(1-eps)/a_v, i.e. 6 V_p/S_p, and the text calls it 'the diameter of a"
      f"\nsphere with the same surface area PER UNIT VOLUME as the actual particle'. For a"
      f"\ncylinder with d = h = 3 mm, S_p/V_p = 6/d exactly, so eq. (11.5.1-9) gives"
      f"\nd_p = {DP_SV:.4f} m. The example instead uses sqrt(1.5) x 0.003 ="
      f" {DP_AREA:.6f} m, which is the\nsphere of equal surface AREA (pi d_s^2 = 1.5 pi d^2)"
      f" - a different definition.\nOn the chapter's own definition the Ergun pressure drop"
      f" is {EX1.dP_Ergun.iloc[2]/1e5:.4f} bar against\nthe printed"
      f" {bk('ex1151a_dp_ergun')/1e5:.3f} bar: {100*DP_DEFINITION_COST:+.1f} %. Reported, not"
      f" repaired.")

# and: no single Reynolds number produces all three printed friction factors
implied = {}
implied["Ergun"] = brentq(lambda R: (1 - E["eps"]) / E["eps"] ** 3 * (
    bk("ergun_a") + bk("ergun_b") * (1 - E["eps"]) / R) - bk("ex1151a_f_ergun"), 50, 500)
implied["McDonald"] = brentq(lambda R: (1 - E["eps"]) / E["eps"] ** 3 * (
    bk("mcdonald_a") + bk("mcdonald_b") * (1 - E["eps"]) / R) - bk("ex1151a_f_mcdonald"), 50, 500)
implied["Hicks"] = brentq(lambda R: bk("hicks_coefficient") * (1 - E["eps"]) ** bk(
    "hicks_exponent_eps") / E["eps"] ** 3 * R ** bk("hicks_exponent_Re")
    - bk("ex1151a_f_hicks"), 50, 500)
RE_IMPLIED_SPREAD = max(implied.values()) - min(implied.values())
print(f"\nA third, smaller printed inconsistency in the same example. Inverting each"
      f"\nprinted friction factor for the Reynolds number it needs, at the printed"
      f" eps = 0.38:")
for k, v in implied.items():
    print(f"   f_{k:<9s} = {bk('ex1151a_f_' + k.lower()):7.2f}  requires  Re = {v:.2f}")
print(f"   the example prints  Re = {bk('ex1151a_Re'):.0f}, and d_p = 0.0037 m gives"
      f" {DP_ROUNDED*E['G']/E['mu']:.2f}."
      f"\nThe Ergun and Hicks branches were evaluated at Re ~= {implied['Ergun']:.1f} and the"
      f"\nMcDonald branch at Re ~= {implied['McDonald']:.2f}; the spread is"
      f" {RE_IMPLIED_SPREAD:.2f}.\nA1.1 carries neither the McDonald nor the"
      f" Hicks constants, so this branch of the\nexample is not covered anywhere else.")


# %%
# ===== 2. the parametric-sensitivity sentence, and every axis that carries error
print("section 11.5.2, book page 513, verbatim:")
print("   " + bk("runaway_claim", "as_printed"))
print(f"   and the operating range on page 511 is {bk('temp_range_low','as_printed')}"
      f" to {bk('temp_range_high','as_printed')} C, which is what\n   'permissible limits'"
      f" means.\n")

R1 = rung1_chain()               # THE one rung-1 computation; agreement.json reads it
CONV1 = R1["table"]
ORD_T, LIM_T = R1["metrics"]["rung1_order_T"], R1["limit_T"]
ORD_P, LIM_P = R1["metrics"]["rung1_order_p0"], R1["limit_p0"]
ORD_Z, LIM_Z = R1["order_z"], R1["limit_z"]
show(CONV1.style.set_uuid("d11conv1").hide(axis="index").format(
    {"T_hot": "{:.6f}", "z_hot": "{:.6f}", "p0_star": "{:.8f}"}))
print(f"observed orders  T_hot {ORD_T:.3f}   z_hot {ORD_Z:.3f}   p0* {ORD_P:.3f}"
      f"   (first-order upwind)\nRichardson limits  T_hot {LIM_T:.5f} K ="
      f" {LIM_T-273.15:.3f} C   z_hot {LIM_Z:.5f} m   p0* {LIM_P:.8f} atm")
R1_FINEST_OVER_FLOOR = R1["metrics"]["rung1_finest_err_over_start_floor"]
print(f"The finest grid's threshold error against that limit is"
      f" {R1_FINEST_OVER_FLOOR:.0f} times the {START_FLOOR:.0e}\nby which a threshold moves"
      f" when only its starting guess changes, so the order fit is not\nfitted to that noise.")

T_ODE, Z_ODE = R1["T_ode"], R1["z_ode"]
P0_ODE = R1["p0_ode"]
ODE_T_DEV = R1["metrics"]["ode_vs_pymrm_T"]
ODE_P_DEV = R1["metrics"]["ode_vs_pymrm_p0"]
print(f"\nSECOND, INDEPENDENT ROUTE - adaptive Radau on the two ODEs, no grid, no pymrm"
      f"\noperator, no Newton, sharing only the constants:"
      f"\n   T_hot  {T_ODE:.6f} K = {T_ODE-273.15:.4f} C   vs the extrapolated pymrm"
      f" {LIM_T-273.15:.4f} C   ->  {ODE_T_DEV:.1e}"
      f"\n   z_hot  {Z_ODE:.6f} m                     vs {LIM_Z:.6f} m"
      f"\n   p0*    {P0_ODE:.8f} atm                 vs {LIM_P:.8f} atm     ->  {ODE_P_DEV:.1e}"
      f"\nA break row perturbs an input; this catches a wrong baseline, and it is the"
      f"\nreason the reported p0* is the extrapolated one and not the 1200-cell value"
      f" ({100*(P0S[1]/LIM_P-1):+.3f} %)."
      f"\nEVERY rung-1 number in agreement.json is an entry of the dict this cell is"
      f"\ndisplaying - one computation, not two that agree.")

T_AT_BOOK_P0 = T_ODE - 273.15
_p0_plus = P0_BOOK + bk("runaway_p0_increment")
print(f"\nSo, on the book's own constants with R_t = {R_T} m, IN ATM:"
      f"\n   p_A0 = {P0_BOOK:.4f} atm          ->  hot spot {T_AT_BOOK_P0:.2f} C, "
      f"{T_PERM-273.15 - T_AT_BOOK_P0:.2f} C BELOW the printed ceiling"
      f"\n   p_A0 = {_p0_plus:.4f} atm (+{bk('runaway_p0_increment'):.4f})  ->  hot spot"
      f" {ode_reference(RX, _p0_plus)[0]-273.15:.1f} C, far above it"
      f"\n   the ceiling is crossed at p_A0 = {P0_ODE:.6f} atm,"
      f" {100*(P0_ODE/P0_BOOK-1):+.3f} % from the printed 0.018."
      f"\nThe sentence is right, to one part in"
      f" {1.0/abs(P0_ODE/P0_BOOK-1):.0f}.")

# ---- WHAT THE SENTENCE IS AN AGREEMENT *WITH*: the two unit readings, measured ----
UNIT_ATM = unit_reading(P["p_B0"], P["p_t"])
UNIT_BAR = unit_reading(bk("oxygen_pressure"), bk("total_pressure"))
P0_STAR_BAR_DEV = UNIT_BAR["p0_star"] / P0_BOOK - 1.0
T_0018_BAR_C = UNIT_BAR["T_at_book_p0"] - 273.15
UNITS = pd.DataFrame([
    dict(reading=f"VWF's set, p_B0 = {P['p_B0']} atm, p_t = {P['p_t']} atm  (THIS PAGE)",
         p0_star=UNIT_ATM["p0_star"], vs_printed_0018=UNIT_ATM["p0_star"] / P0_BOOK - 1.0,
         T_hot_at_0018_C=UNIT_ATM["T_at_book_p0"] - 273.15),
    dict(reading=f"11.5.2 as printed, p_B0 = {bk('oxygen_pressure')} bar,"
                 f" p_t = {bk('total_pressure'):.0f} bar",
         p0_star=UNIT_BAR["p0_star"], vs_printed_0018=P0_STAR_BAR_DEV,
         T_hot_at_0018_C=T_0018_BAR_C),
])
show(UNITS.style.set_uuid("d11units").hide(axis="index").format(
    {"p0_star": "{:.6f}", "vs_printed_0018": "{:+.5f}", "T_hot_at_0018_C": "{:.1f}"}))
print(f"Both rows are the same grid-free Radau solve with one constant changed, so the"
      f"\ndifference is a PARAMETER-SET difference and carries no discretisation at all."
      f"\nThe +{100*UNIT_ATM['p0_star']/P0_BOOK-100:.3f} % reproduction of section 11.5.2's"
      f" sensitivity sentence is an agreement with"
      f"\nVan Welsenaere & Froment's 0.018 ATM, obtained with THEIR atm parameter set. Read"
      f"\nin the bar that section 11.5.2 itself prints, the same sentence is out by"
      f" {100*P0_STAR_BAR_DEV:.3f} %\nand its own design point sits at"
      f" {T_0018_BAR_C:.0f} C -"
      f" {T_0018_BAR_C-(T_PERM-273.15):.0f} K ABOVE the ceiling the facing page"
      f"\nstates, which no reading of the chapter can be right about. That is why the atm"
      f"\nreading is the one used, and it is stated rather than assumed:"
      f" Fig. 11.5.2-1's own\nordinate label reads"
      f" '{bk('fig1152_1_ordinate_unit','as_printed')}'.")

print(f"\nAt R_t = {D_T_PRINTED/2} m - the radius the printed 2.54 cm diameter implies -"
      f" the same\ninlet pressure gives a hot spot of"
      f" {ode_reference(RX, P0_BOOK, R=D_T_PRINTED/2)[0]-273.15:.0f} C, outside the operating range"
      f" the same page\nstates, and p0* falls to"
      f" {p0_critical_ode(RX, R=D_T_PRINTED/2):.6f} atm. The figures of section 11.5.2\ncannot have"
      f" been computed at that radius; R_t = {R_T} m is what reproduces them, and\nit is also"
      f" what eq. (11.7.4-1) needs to return the printed U. Two independent\nconstraints, same"
      f" answer.")

# ---- IS THE THRESHOLD A CEILING CROSSING OR A NEWTON-CONVERGENCE BOUNDARY? -------
print(f"\n\nCONTROL: p0_critical bisects and treats a non-converged solve as"
      f" supercritical, so\nthe number it returns COULD be a Newton-convergence boundary"
      f" rather than a crossing\nof the {T_PERM-273.15:.0f} C ceiling. It is checked, on"
      f" EVERY call this page makes - not\nonly on these five rungs:"
      f" `p0_critical` runs `crossing_control`"
      f" itself and RAISES if it fails, so a\nthreshold that is not a ceiling crossing"
      f" cannot reach a table, a sweep or a\nbisection. What follows is that same check,"
      f" displayed for the five reported rungs.")
ctl = []
for r, sec, st, what in RUNGS:
    lad, pc = LAD[r], P0S[r]
    lo_, hi_, ok_ = crossing_control(lad, pc)
    ctl.append(dict(rung=r, p0_star=pc, T_just_below_C=lo_ - 273.15,
                    T_just_above_C=hi_ - 273.15, jump_K=hi_ - lo_, crossing=ok_))
CTL = pd.DataFrame(ctl)
show(CTL.style.set_uuid("d11ctl").hide(axis="index").format(
    {"p0_star": "{:.7f}", "T_just_below_C": "{:.2f}", "T_just_above_C": "{:.2f}",
     "jump_K": "{:.3f}"}, na_rep="Newton did not converge"))
CROSSING_OK = bool(CTL.crossing.all())
assert CROSSING_OK, "a rung's threshold is NOT a ceiling crossing"
print(f"All five cross the ceiling, and cross it CONTINUOUSLY: the widest jump between"
      f" the\ntwo solves is {CTL.jump_K.abs().max():.3f} K, against the"
      f" {CTL_JUMP:.0f} K the control allows and the hundreds of\nkelvin an ignition"
      f" discontinuity gives. So no rung's p0* is set by a convergence\nfailure. Above"
      f" the threshold the Newton may give up - the situation bisection exists to\nsurvive -"
      f" and that is NOT an absence of a steady state: the grid-free Radau route\nreturns a"
      f" bounded one at the top of the bracket,"
      f" {ode_reference(RX, HI[1])[0]-273.15:.0f} C at p_A0 = {HI[1]}."
      f"\n\nBrentq on the identical bracket, on ALL FIVE rungs at the n_z = 600 the break"
      f" rows\nuse - this is what licenses bisection, so it is measured on every rung"
      f" rather than\non the two it is cheapest to show:")
_BRENT_OK, _BRENT_FAIL = [], []
for r in (1, 2, 3, 4, 5):
    _bs = p0_critical(Ladder(RX, rung=r, n_z=600, n_r=10), hi=HI[r])   # always, first
    try:
        _rt = brentq(lambda q: T_max_safe(Ladder(RX, rung=r, n_z=600, n_r=10), q) - T_PERM,
                     0.004, HI[r], xtol=1e-12)
        assert abs(_rt / _bs - 1.0) < 1e-10, f"rung {r}: brentq and bisection disagree"
        _BRENT_OK.append(r)
    except ValueError:
        _BRENT_FAIL.append(r)
    print(f"   rung {r}: bisection {_bs:.8f}")
assert _BRENT_FAIL, "brentq survived on every rung, so the text below is not what ran"
print(f"Brentq raises on a NaN on at least one rung, and on every rung where it survives it"
      f"\nreturns the bisection root to better than 1e-10 (asserted). WHICH rungs it survives"
      f"\nis decided by round-off, not by the rung: its points above the threshold are"
      f"\ncold-started fixed-inlet solves, whose convergence is the round-off-decided"
      f" event the\nnote on locating thresholds describes (here it fails on rungs"
      f" {', '.join(map(str, _BRENT_FAIL))}). So brentq is not wrong\nwhere it survives;"
      f" it simply cannot be relied on.")
_T5_COARSE = T_max_safe(Ladder(RX, rung=5, n_z=600, n_r=10), P0_BOOK)
assert np.isfinite(_T5_COARSE) and _T5_COARSE - 273.15 > 900.0, (
    "rung 5 at 0.018 atm on (600, 10) is not the bounded runaway state the text says")
print(f"\nAnd rung 5 at p_A0 = {P0_BOOK} atm, past its own threshold: on the coarser (600, 10)"
      f"\ngrid the cold-started solve reaches a bounded runaway state above 900 C (asserted),"
      f"\nmore than 480 K past the ceiling (here {_T5_COARSE-273.15:.0f} C), so the reactor HAS"
      f" a steady\nstate there. Its exact temperature, and whether the (1200, 20) grid"
      f" converges at all,\nis decided by round-off: this is a fixed-inlet solve above a"
      f" threshold, started\ncold, which branch following does not cover.")


# %%
# ---- rung 5 carries error on TWO axes: refine both, separately ---------------
r5 = []
for n_z in (300, 600, 1200):
    for n_r in (5, 10, 20, 40):
        lad = Ladder(RX, rung=5, n_z=n_z, n_r=n_r)
        c = solve_robust(lad, P0_REF)
        Tb, zb, _ = hotspot(lad.z_c, lad.bulk_T(c))
        Ta, za, _ = hotspot(lad.z_c, lad.axis_T(c))
        i = int(np.argmax(lad.axis_T(c)))
        r5.append(dict(n_z=n_z, n_r=n_r, T_bulk=Tb, T_axis=Ta,
                       dT_radial=float(c[i, 0, 1] - c[i, -1, 1])))
R5 = pd.DataFrame(r5)
show(R5.style.set_uuid("d11conv5").hide(axis="index").format(
    {"T_bulk": "{:.4f}", "T_axis": "{:.4f}", "dT_radial": "{:.4f}"}))
_z = R5[R5.n_r == 40]
_r = R5[R5.n_z == 1200]
ORD_Z5, LIM_Z5 = richardson(_z.T_axis, _z.n_z)
ORD_R5, LIM_R5 = richardson(_r.T_axis, _r.n_r)
T_AXIS_REPORTED = float(R5[(R5.n_z == N_Z) & (R5.n_r == N_R)].T_axis.iloc[0])
ERR_RADIAL = T_AXIS_REPORTED - LIM_R5
ERR_AXIAL = float(_z.T_axis.iloc[-1]) - LIM_Z5
AXIS_MINUS_BULK = T_AXIS_REPORTED - (LADDER.loc[4, "T_hot"] + 273.15)
print(f"axial axis (n_r = 40): observed order {ORD_Z5:.3f}, limit {LIM_Z5-273.15:.4f} C"
      f"\nradial axis (n_z = 1200): observed order {ORD_R5:.3f}, limit"
      f" {LIM_R5-273.15:.4f} C"
      f"\nThe radial axis converges at second order and the axial at first, so the"
      f"\ntwo-dimensional rung's error is set by the axis the model is NOT named after."
      f"\nBOTH ARE TINY. Comparing LIKE WITH LIKE - the AXIS temperature against an AXIS"
      f"\nlimit, on the same axis that was refined:"
      f"\n   radial discretisation error at ({N_Z}, {N_R}) : {ERR_RADIAL:+.4f} K"
      f"\n   axial  discretisation error at ({N_Z}, 40) : {ERR_AXIAL:+.4f} K"
      f"\ni.e. both below 0.03 K on a {LIM_R5-273.15:.0f} C hot spot (asserted), and the"
      f" axial one is\n{abs(ERR_AXIAL/ERR_RADIAL):.0f} times the radial one.")
assert max(abs(ERR_RADIAL), abs(ERR_AXIAL)) < 0.03, "a rung-5 axis error is not below 0.03 K"
print(f"\nWHAT IS *NOT* A GRID ERROR:"
      f"\nthe axis-minus-bulk difference at the reported grid,"
      f" {AXIS_MINUS_BULK:+.3f} K"
      f"\n({T_AXIS_REPORTED-273.15:.2f} C on the axis against"
      f" {LADDER.loc[4,'T_hot']:.2f} C for the radial mean). That is the radial"
      f"\ntemperature difference the two-dimensional model exists to resolve - the physics,"
      f"\nnot the discretisation - and it is {abs(AXIS_MINUS_BULK/ERR_RADIAL):.0f} times the"
      f" radial grid error beside it. The\nsame quantity measured wall-to-axis is the"
      f" dT_radial column above,"
      f" {float(R5[(R5.n_z==N_Z)&(R5.n_r==N_R)].dT_radial.iloc[0]):.2f} K.")
print(f"\nNOT REFINED ANYWHERE, and said so: rung 5's THRESHOLD, which is the quantity the"
      f"\nladder table reports. Its (600, 10) and ({N_Z}, {N_R}) values differ by about a"
      f"\ntenth of a percentage point - printed in the break-table cell below - and no"
      f"\nobserved order or Richardson limit is offered for it. The two refinements above"
      f"\nare of the hot spot at a fixed inlet, not of the threshold.")

# ---- rung 2: the difference it makes against the numerical diffusion of upwind
DISP_LEN = D_P / 1.5
print(f"\nRUNG 2 NEEDS THIS SAID OUT LOUD. First-order upwind adds a numerical"
      f"\ndiffusivity u dz/2, i.e. a dispersion LENGTH dz/2 - A2.1's finding, applied"
      f"\nhere. The physical dispersion length of (11.6-1) divided by u_s is"
      f"\nd_p/Pe_ma = {DISP_LEN*1e3:.3f} mm. At n_z = {N_Z} the numerical one is"
      f" {0.5*L_BED/N_Z*1e3:.3f} mm, i.e."
      f" {100*0.5*L_BED/N_Z/DISP_LEN:.0f} %\nof it, so a rung-1-versus-rung-2 comparison on a"
      f" single grid is mostly scheme.")
d2 = []
for n in (300, 600, 1200, 2400, 4800):
    l1, l2 = Ladder(RX, rung=1, n_z=n), Ladder(RX, rung=2, n_z=n)
    p1, p2 = p0_critical(l1), p0_critical(l2)
    d2.append(dict(n_z=n, numerical_mm=0.5 * L_BED / n * 1e3, p0_rung1=p1,
                   p0_rung2=p2, shift_pct=100.0 * (p2 / p1 - 1.0)))
D2 = pd.DataFrame(d2)
show(D2.style.set_uuid("d11r2conv").hide(axis="index").format(
    {"numerical_mm": "{:.4f}", "p0_rung1": "{:.8f}", "p0_rung2": "{:.8f}",
     "shift_pct": "{:+.4f}"}))
ORD_S2, LIM_S2 = richardson(D2.shift_pct, D2.n_z)
print(f"The shift itself converges (observed order {ORD_S2:.2f}) to {LIM_S2:+.4f} % as"
      f" h -> 0.\nIt is smaller than the {100*abs(P0S[1]/LIM_P-1):.3f} % discretisation"
      f" error of the rung-1 baseline at\nn_z = {N_Z}, so this rung is only resolvable if both"
      f" models are extrapolated -\nwhich is exactly why the ladder table above reports"
      f" root-found thresholds and\nthis cell reports their limit.")


# %%
# ===== 3. Example 11.5.3.A against Van Welsenaere & Froment's own numbers =====
# D2.2 owns the criteria; this checks the BOOK'S TRANSCRIPTION of their example.
tr = []
for bkey, ex_q in (("vwf_TM_625", "T_cr"), ("vwf_dT_625", "dT_eff"), ("vwf_Q_625", "Q"),
                   ("vwf_p0_lower_625", "p0_lower"), ("vwf_p0_upper_625", "p0_upper"),
                   ("vwf_p0_mean_625", "p0_mean"), ("vwf_p0_backint_625", "p0_critical")):
    src = VWF_EX[(VWF_EX.example.isin(["1a", "4"])) & (VWF_EX.quantity == ex_q)]
    o = float(src.value.iloc[0])
    # WHICH of Van Welsenaere & Froment's examples the row was found in, and at what
    # radius, is printed - because a row-by-row ratio of 1.000000 says nothing about
    # whether the two rows belong to the SAME case, and here two of them do not.
    tr.append(dict(quantity=ex_q, VWF_example=str(src.example.iloc[0]),
                   VWF_conditions=str(src.conditions.iloc[0])[:34],
                   VWF_1970=o, unit=str(src.unit.iloc[0]),
                   book=bk(bkey), ratio=bk(bkey) / o))
TRANS = pd.DataFrame(tr)
show(TRANS.style.set_uuid("d11trans").hide(axis="index").format(
    {"VWF_1970": "{:.6g}", "book": "{:.6g}", "ratio": "{:.6f}"}))
N_EXAMPLES = TRANS.VWF_example.nunique()
print(f"Every dimensionless and kelvin row transcribes unchanged and every ATM row is"
      f"\nmultiplied by {ATM} - except the last, which is not: the atm/bar audit again,"
      f"\nfrom the other side."
      f"\n\nBUT LOOK AT THE `VWF_example` COLUMN. These seven rows are looked up in"
      f" whichever of the\n1970 paper's worked examples happens to carry each quantity, and"
      f" they come from\n{N_EXAMPLES} of them, at two different tube radii. A ratio of"
      f" 1.000000 says the number was\ncopied correctly; it says nothing about whether the"
      f" rows belong to the SAME case.\nQ is the one where they do not, and the next block is"
      f" what this check misses.")

# Table 11.5.3.A-2's own arithmetic, on its own printed Q and dT
Q, DT = bk("vwf_Q_625"), bk("vwf_dT_625")
t2 = pd.DataFrame([
    ("Lower limit", "dT (1 + Q^2)", DT * (1 + Q ** 2), bk("vwf_dTad_lower_625")),
    ("Upper limit", "dT (1 + Q)^2", DT * (1 + Q) ** 2, bk("vwf_dTad_upper_625")),
    ("Mean", "dT (1 + Q + Q^2)", DT * (1 + Q + Q ** 2), bk("vwf_dTad_mean_625")),
], columns=["row", "formula printed on that row", "the formula gives", "the table prints"])
t2["ratio"] = t2["the table prints"] / t2["the formula gives"]
show(t2.style.set_uuid("d11t2").hide(axis="index").format(
    {"the formula gives": "{:.2f}", "the table prints": "{:.2f}", "ratio": "{:.4f}"}))
T2_UPPER_RATIO = float(t2.ratio.iloc[1])
_mean_of_printed = 0.5 * (bk("vwf_dTad_lower_625") + bk("vwf_dTad_upper_625"))
print(f"Row 1 checks. Row 2 does not: dT(1+Q)^2 = {DT*(1+Q)**2:.1f} C, and the"
      f" {bk('vwf_dTad_upper_625')} C\nprinted beside it is what row 3's formula gives"
      f" ({DT*(1+Q+Q**2):.1f} C). Row 3 then prints\n"
      f"{bk('vwf_dTad_mean_625')} C, which is the arithmetic mean of the two printed"
      f" numbers ({_mean_of_printed:.1f} C) and\nwhich no formula in the table produces."
      f" One slip, propagating once.\nNote that dT(1+Q+Q^2) IS identically the mean of"
      f" dT(1+Q^2) and dT(1+Q)^2 - the\nidentity holds for every Q - so the 'Mean' row's"
      f" formula is right and only its\nnumber is wrong. Table 11.5.3.A-1, at 635 K, is"
      f" internally consistent throughout:")
Q1, DT1 = bk("vwf_Q_635"), bk("vwf_TM_635") - 635.0
t1 = pd.DataFrame([
    ("Lower limit", DT1 * (1 + Q1 ** 2), bk("vwf_dTad_lower_635")),
    ("Upper limit", DT1 * (1 + Q1) ** 2, bk("vwf_dTad_upper_635")),
    ("Mean", DT1 * (1 + Q1 + Q1 ** 2), bk("vwf_dTad_mean_635")),
], columns=["row", "the formula gives", "the table prints"])
t1["rel. dev."] = t1["the formula gives"] / t1["the table prints"] - 1.0
show(t1.style.set_uuid("d11t1").hide(axis="index").format(
    {"the formula gives": "{:.2f}", "the table prints": "{:.2f}", "rel. dev.": "{:+.4f}"}))
T1_WORST = float(np.abs(t1["rel. dev."]).max())
print(f"worst {T1_WORST:.3f} on a table whose entries are printed to 3 figures.")

# ===== the LARGER defect in the same table, which no row-by-row check can see =====
# Table 11.5.3.A-2 has TWO columns, and the example prints the constant that ties
# them together: dT_ad = (B/A) p0.  ONE DIVISION on the example's own printed A and B.
BOA = bk("vwf_B") / bk("vwf_A")                 # per ATM, the units A and B are in
xc = []
for lab, pkey, dkey, form in (
        ("Lower limit", "vwf_p0_lower_625", "vwf_dTad_lower_625", lambda q: DT * (1 + q ** 2)),
        ("Upper limit", "vwf_p0_upper_625", "vwf_dTad_upper_625", lambda q: DT * (1 + q) ** 2),
        ("Mean", "vwf_p0_mean_625", "vwf_dTad_mean_625", lambda q: DT * (1 + q + q ** 2))):
    p0_atm = bk(pkey) / ATM                      # the column is printed in bar
    needed = p0_atm * BOA
    qi = brentq(lambda q: form(q) - needed, 0.1, 20.0, xtol=1e-13)
    xc.append(dict(row=lab, p0_printed_bar=bk(pkey), dTad_needed_K=needed,
                   dTad_printed_K=bk(dkey), dev=bk(dkey) / needed - 1.0, Q_implied=qi))
XC = pd.DataFrame(xc)
show(XC.style.set_uuid("d11t2cross").hide(axis="index").format(
    {"p0_printed_bar": "{:.5f}", "dTad_needed_K": "{:.1f}", "dTad_printed_K": "{:.1f}",
     "dev": "{:+.3f}", "Q_implied": "{:.4f}"}))
T2_DTAD_WORST_DEV = float(np.abs(XC.dev).max())
T2_Q_IMPLIED = float(XC.Q_implied.iloc[0])
T1_CROSS_WORST = max(abs(bk(dk) / (bk(pk) / ATM * BOA) - 1.0) for pk, dk in
                     (("vwf_p0_lower_635", "vwf_dTad_lower_635"),
                      ("vwf_p0_upper_635", "vwf_dTad_upper_635"),
                      ("vwf_p0_mean_635", "vwf_dTad_mean_635")))
print(f"THE TWO COLUMNS OF TABLE 11.5.3.A-2 DISAGREE. dT_ad = (B/A) p0 with the example's"
      f"\nown printed A = {bk('vwf_A'):.0f} and B = {bk('vwf_B'):.3g}, i.e. B/A ="
      f" {BOA:.0f} per atm, and the printed"
      f"\ndT_ad column is {100*T2_DTAD_WORST_DEV:.0f} % low at worst. Invert each row for the Q"
      f" it needs and all\nthree agree: Q = {XC.Q_implied.min():.4f} to"
      f" {XC.Q_implied.max():.4f}, against the printed {Q}."
      f"\nTABLE 11.5.3.A-1 PASSES THE SAME TEST - worst {100*T1_CROSS_WORST:.1f} % on the same"
      f" three rows,\nagainst {100*T2_DTAD_WORST_DEV:.0f} % - which rules out a wrong B/A, a"
      f" wrong conversion and a wrong\nreading of the formulae, since every one of those"
      f" would move BOTH tables.\nIt does NOT confine the defect to A-2:"
      f" the block below inverts each table's own printed Q"
      f" for the radius it\nimplies, and both of them point at the same 0.0175 m.")

# WHERE the printed Q comes from - the example's own second half, two paragraphs on
Q_SCALED = Q * np.sqrt(bk("vwf_R_critical_625") / R_T)
C_AT = lambda R: 2.0 * P["U"] / (P["cp_vol"] * R)
Q_FROM_C = lambda R, A, TM=None: np.sqrt(
    C_AT(R) / (A * np.exp(P["b"] - P["a"] / (bk("vwf_TM_625") if TM is None else TM))))
_r4 = VWF_EX[(VWF_EX.example == "4") & (VWF_EX.quantity == "R")]
print(f"\nWHY. The same Example 11.5.3.A asks, two paragraphs later, what radius would be"
      f"\ncritical, and answers - book page"
      f" {int(BK.loc['vwf_R_critical_625','book_page'])} -"
      f" '{bk('vwf_C_formula_claim','as_printed')}'"
      f"\nD2.2's CSV of the 1970 original carries the same answer as its Example 4:"
      f" R = {float(_r4.value.iloc[0])} m,\n\"{str(_r4.conditions.iloc[0])}\"."
      f" So Q = {Q} IS a printed value of this example - it is the Q of the"
      f"\nCRITICAL radius {bk('vwf_R_critical_625')} m, transplanted onto the"
      f" {R_T} m case that Table 11.5.3.A-2\ntabulates."
      f"\nThat printed C formula is itself the proof, because C = 2U/(c_p R_t) - with the"
      f"\nVOLUMETRIC c_p, D2.2's finding - is"
      f"\ninversely proportional to R_t at a fixed T_M, so Q scales as R_t^(-1/2):"
      f"\n   Q at R_t = {R_T} m, from Q({bk('vwf_R_critical_625')}) x sqrt(R_crit/R_t)   :"
      f" {Q_SCALED:.4f}"
      f"\n   Q at R_t = {R_T} m, from C = Q^2 A exp(b - E/RT_M) directly : "
      f"{Q_FROM_C(R_T, A_G):.4f}   (A = {A_G:.0f})"
      f"\n   Q at R_t = {bk('vwf_R_critical_625')} m by the same formula, against the printed"
      f" {Q}: {Q_FROM_C(bk('vwf_R_critical_625'), A_G):.4f}"
      f"\n   Q the p0 column of Table 11.5.3.A-2 needs                 : {T2_Q_IMPLIED:.4f}"
      f"\nFour routes, one answer: the case tabulated needs Q ~ {T2_Q_IMPLIED:.2f} and the"
      f" table was computed\nwith {Q}. The p0 column is Van Welsenaere & Froment's own,"
      f" converted;"
      f"\nthe dT_ad column beside it was recomputed with a Q from a different radius."
      f"\nThis is LARGER than the Upper/Mean slip above, and it is invisible to the"
      f"\ntranscription table: every one of those rows transcribes at ratio 1.000000,"
      f"\nbecause the number was copied correctly - from the wrong case. Reported, not"
      f"\nrepaired: nothing on this page uses Q.")

# ---- AND THE SAME TEST ON TABLE A-1, WHICH THE CROSS-COLUMN CHECK CANNOT SEE ----
# The C route is validated on Van Welsenaere & Froment's own printed numbers first:
# their Q at the critical radius, and their printed C, which pins R by itself.
_C4 = float(VWF_EX[(VWF_EX.example == "4") & (VWF_EX.quantity == "C")].value.iloc[0])
R_FROM_C_PRINTED = 2.0 * P["U"] / (P["cp_vol"] * _C4)
R_IMPLIED = lambda Qp, TM: R_T * (Q_FROM_C(R_T, A_G, TM) / Qp) ** 2
QINV = pd.DataFrame([
    dict(table="11.5.3.A-2", T_M_K=bk("vwf_TM_625"), Q_printed=Q,
         Q_at_0_0125=Q_FROM_C(R_T, A_G, bk("vwf_TM_625")),
         R_t_implied=R_IMPLIED(Q, bk("vwf_TM_625"))),
    dict(table="11.5.3.A-1", T_M_K=bk("vwf_TM_635"), Q_printed=Q1,
         Q_at_0_0125=Q_FROM_C(R_T, A_G, bk("vwf_TM_635")),
         R_t_implied=R_IMPLIED(Q1, bk("vwf_TM_635"))),
])
show(QINV.style.set_uuid("d11qinv").hide(axis="index").format(
    {"T_M_K": "{:.2f}", "Q_printed": "{:.4f}", "Q_at_0_0125": "{:.4f}",
     "R_t_implied": "{:.6f}"}))
R_A1_IMPLIED = float(QINV.R_t_implied.iloc[1])
print(f"The route is Van Welsenaere & Froment's own, and it is checked against their own"
      f"\nprinted numbers before it is inverted: their C = {_C4:.1f} 1/h with"
      f" C = 2U/(c_p R_t) alone\ngives R ="
      f" {R_FROM_C_PRINTED:.6f} m against the {bk('vwf_R_critical_625')} m the book prints"
      f" beside it, and\nQ({bk('vwf_R_critical_625')}) = {Q_FROM_C(bk('vwf_R_critical_625'), A_G):.4f}"
      f" against their printed {Q}."
      f"\n\nBOTH TABLES' Q BELONG TO R_t = {bk('vwf_R_critical_625')} m."
      f" Table A-1's printed Q = {Q1} implies\nR_t ="
      f" {R_A1_IMPLIED:.6f} m, the same radius as A-2's {Q}, and A-1 states the same"
      f" {R_T} m\nas A-2 does. So the transplant is in BOTH tables, and the reason only"
      f" A-2's columns\ndisagree is that A-1's two columns are not independent of each"
      f" other: they agree to\n{100*T1_CROSS_WORST:.1f} %, which is what a column"
      f" back-calculated from the other looks like, while\nA-2's p0 column is Van"
      f" Welsenaere & Froment's own - the transcription table above\nshows those three"
      f" entries at ratio {ATM} exactly - and its dT_ad column is not."
      f"\nA CROSS-COLUMN TEST CAN ONLY SEE A TRANSPLANT THAT MOVED ONE COLUMN. That is a"
      f"\nlimit of the test: A-1's pass does not confine the defect to A-2.")


# %%
# ======== 4. eq. (11.9.1-11) against a pellet grid solve and B1.1's table =====
# the pellet is model.py's `Pellet`: a spherical grid (`nu = 2`), surface values
# imposed, eta from the surface-flux row of the gradient operator.

# (a) the pellet grid solve against eq. (11.9.1-11), refined
# the pellet is examined at the rung-1 hot spot of section 11.5.2's OWN design point,
# which is the hottest state any rung reaches while still having a solution
_c1hs = Ladder(RX, rung=1, n_z=N_Z).solve(P0_BOOK)
_ihs = int(np.argmax(_c1hs[:, 1]))
T_HS = float(_c1hs[_ihs, 1])
P_HS = float(_c1hs[_ihs, 0])
PHI_HS = phi_prime(RX, T_HS, 1.0e-6)
pv = []
for n in (16, 32, 64, 128, 256):
    e = Pellet(RX, n_r=n).solve(THETA * P_HS, T_HS)["eta"]
    pv.append(dict(n_r=n, eta_numeric=e, eta_11_9_1_11=eta_sphere(PHI_HS),
                   rel_dev=e / eta_sphere(PHI_HS) - 1.0))
PV = pd.DataFrame(pv)
show(PV.style.set_uuid("d11pellet").hide(axis="index").format(
    {"eta_numeric": "{:.9f}", "eta_11_9_1_11": "{:.9f}", "rel_dev": "{:+.2e}"}))
ORD_ETA, LIM_ETA = richardson(PV.eta_numeric, PV.n_r)
ETA_DEV = float(abs(PV.rel_dev.iloc[-1]))
print(f"phi'' = {PHI_HS:.6f} at the rung-1 hot spot ({T_HS-273.15:.2f} C,"
      f" p_A = {P_HS:.6f} atm)."
      f"\nobserved order {ORD_ETA:.3f}, limit {LIM_ETA:.9f}, against eq. (11.9.1-11)'s"
      f" {eta_sphere(PHI_HS):.9f}\n  -> {abs(LIM_ETA/eta_sphere(PHI_HS)-1):.1e}. The analytical"
      f" solution and the grid solve agree.")

# (b) the same closed form against B1.1's exact table, which is an INDEPENDENT
#     transcription of it from Thiele (1939) rather than from Froment eq. (11.9.1-11)
SPH = ETA_REF[ETA_REF.geometry == "sphere"]
B11_DEV = float(np.abs(eta_sphere(SPH.phi.values) / SPH.eta.values - 1.0).max())
print(f"\nEq. (11.9.1-11) against B1.1's {len(SPH)} exact sphere rows, phi ="
      f" {SPH.phi.min():.2f} to {SPH.phi.max():.0f}:"
      f"\n   worst relative deviation {B11_DEV:.2e}. Same closed form, two independent"
      f"\n   transcriptions from two documents. B1.1's own pymrm pellet matches its table"
      f"\n   to {B11_C['isothermal_max_rel_error_phi_le_30']:.1e} for phi <= 30 and"
      f" {B11_C['isothermal_max_rel_error_all_phi']:.1e} over all phi - B1.1's claims, and"
      f"\n   this page's {ETA_DEV:.1e} at phi = {PHI_HS:.2f} is consistent with it.")

# (c) the book's claim that the particle is practically isothermal
print(f"\nsection 11.9.1, book page 598, verbatim:\n   "
      + bk("isothermal_pellet_claim", "as_printed"))
niso = Pellet(RX, n_r=128, isothermal=False).solve(THETA * P_HS, T_HS)
PELLET_DT = niso["T_centre"] - T_HS
BETA_PRATER = (P["minus_dH"] * KCAL) * 1.0e-6 * (THETA * P_HS) / (1.0e-3 * T_HS)
GAMMA_ARRH = P["a"] / T_HS
_c3hs = solve_robust(Ladder(RX, rung=3, n_z=N_Z), P0S[3])
FILM_DT_HS = float(np.max(_c3hs[:, 3] - _c3hs[:, 1]))
print(f"   intraparticle rise, centre minus surface : {PELLET_DT:+.4f} K"
      f"\n   interfacial rise across the film, at rung 3's own threshold: {FILM_DT_HS:+.4f} K"
      f"\n   ratio                                    : {FILM_DT_HS/PELLET_DT:.2f}"
      f"\n   Prater number beta = {BETA_PRATER:.3e}, Arrhenius gamma = {GAMMA_ARRH:.2f},"
      f" gamma*beta = {GAMMA_ARRH*BETA_PRATER:.4f}"
      f"\n   eta with the pellet energy balance on: {niso['eta']:.8f}, against"
      f" {eta_sphere(PHI_HS):.8f} isothermal"
      f"\nBoth halves of the claim hold on this reactor, and the second half is the"
      f"\nstronger: the film carries {FILM_DT_HS/PELLET_DT:.0f}x the temperature difference the"
      f" particle does. This\nsits far below B1.1's fold (beta = 0.6, gamma = 20), so none of"
      f" its multiplicity\nis in play here and none of it is borrowed.")


# %%
# ===== 5. the chapter's own criteria for deciding whether to climb a rung =====
print("Section 11.6 offers three yardsticks for whether rung 2 is needed. Two of them"
      "\ncannot be evaluated on the reactor section 11.5.2 designs.\n")
print("(i) the rule of thumb, book page 560:\n   " + bk("axial_50dp_claim", "as_printed"))
L_OVER_DP = L_BED / D_P
PE_L = L_OVER_DP * bk("pe_axial_low"), L_OVER_DP * bk("pe_axial_high")
print(f"   L/d_p = {L_OVER_DP:.0f}, which is {L_OVER_DP/bk('axial_50dp_claim'):.0f}x the"
      f" rule's 50, and the length Peclet number\n   Pe_a' = Pe_a L/d_p ="
      f" {PE_L[0]:.0f} to {PE_L[1]:.0f}. Page 563 says axial mixing would matter in"
      f"\n   methanol synthesis below Pe_a' = {bk('pe_ma_methanol'):.0f} and in ethylene"
      f" oxidation below {bk('pe_ma_ethylene'):.0f}."
      f"\n   MEASURED SHIFT, extrapolated: {LIM_S2:+.4f} % in p0*. The rule of thumb is"
      f" right.\n")

print("(ii) criterion (11.6-3), for a rate that decreases monotonically:")
r_A0 = kexp(T_R0) * P["p_B0"] * P0_BOOK                      # kmol/(kg cat h)
C0 = THETA * P0_BOOK
CRIT_MASS = r_A0 * P["rho_b"] * D_P / (P["u"] * C0)
print(f"      r_A0 rho_B d_p / (u_s C_0)  =  {CRIT_MASS:.3e}   <<  Pe_ma ="
      f" {bk('pe_axial_low'):.0f} to {bk('pe_axial_high'):.0f}   satisfied by"
      f" {bk('pe_axial_low')/CRIT_MASS:.1e}")
print(f"      (-dH) r_A0 rho_B d_p / [(T_0 - T_w) u_s rho_g c_p]  <<  Pe_ha")
print(f"      T_0 - T_w = {T_R0 - T_R0:.1f} K on this reactor, BY CONSTRUCTION: section"
      f" 11.5.2 sets\n      T = T_0 = T_r at z = 0, and section 11.5.3 builds its whole"
      f" runaway diagram\n      for 'the common situation where T_r = T_0'. The heat"
      f" criterion divides by\n      zero and cannot be evaluated for the chapter's own"
      f" design case.\n")

print("(iii) criterion (11.6-4), for a rate with an interior maximum, book page 561:")
lad1 = Ladder(RX, rung=1, n_z=4800)
c1 = lad1.solve(P0_BOOK)
dx_dzdp = float(np.max(np.abs(np.gradient(1.0 - c1[:, 0] / P0_BOOK, lad1.z_c))) * D_P)
dT_dzdp = float(np.max(np.abs(np.gradient(c1[:, 1], lad1.z_c))) * D_P)
print(f"      max |dx / d(z/d_p)|  =  {dx_dzdp:.4e}   <<  Pe_ma   satisfied by"
      f" {bk('pe_axial_low')/dx_dzdp:.0f}x")
print(f"      max |dT / d(z/d_p)|  =  {dT_dzdp:.4f} K   <<  Pe_ha")
print(f"      AS PRINTED THIS IS DIMENSIONALLY INHOMOGENEOUS: the left side is a"
      f" temperature\n      per particle diameter and the right side is a dimensionless"
      f" group. Its own\n      companion in (11.6-3) is normalised by (T_0 - T_w) and is"
      f" dimensionless;\n      (11.6-4)'s mass criterion is dimensionless too, because x"
      f" is. The heat one\n      is missing a temperature scale, and there is no candidate"
      f" in the chapter:\n      (T_0 - T_w) is zero here.")
print(f"\n      the book asserts, on the same page:\n      \"" +
      bk("axial_steep_claim", "as_printed") + "\"")
print(f"      Taking the natural repair - dividing by the adiabatic temperature rise,"
      f"\n      LABELLED AS AN INFERENCE and not as the book's - gives"
      f" {dT_dzdp/(B_G/A_G*P0_BOOK):.3e}, satisfied\n      by"
      f" {bk('pe_axial_low')/(dT_dzdp/(B_G/A_G*P0_BOOK)):.0f}x. So the assertion is"
      f" correct on any reading; it is the criterion,\n      not the conclusion, that cannot"
      f" be used.")
print(f"\n(iv) THE CHAPTER IS NOT THE ONLY PLACE THE BOOK DISCUSSES THESE. Section 12.7.2,"
      f"\n     book page {int(BK.loc['ch12_yf_claim','book_page'])},"
      f" {int(BK.loc['ch12_yf_claim','book_page']) - 561} pages later, says of the very same"
      f" criteria:\n     \"" + bk("ch12_yf_claim", "as_printed") + "\"")
print("     Read that carefully, because the page's gloss on it is an INFERENCE and is"
      "\n     labelled as one. What 12.7.2 STATES is that Mears found the criteria not"
      "\n     general and supplied alternates for equal feed and wall temperatures. That"
      "\n     those alternates exist BECAUSE of the two problems found above - a division"
      "\n     by (T_0 - T_w) that is zero when the feed and wall temperatures are equal,"
      "\n     and a dimensionally inhomogeneous heat criterion - is this page's reading of"
      "\n     the sentence, not something the sentence says. The one thing that is not an"
      "\n     inference is that the situation it singles out, equal feed and wall"
      "\n     temperatures, is exactly the situation section 11.5.2 builds. Section 11.6"
      "\n     carries no cross-reference to 12.7.2 and 12.7.2 none to 11.6; the nearest"
      "\n     thing is 12.7.2's own next paragraph, 'the latter has been thoroughly covered"
      "\n     in Chapter 11', which is about RADIAL dispersion, not about the criteria."
      "\n     (Searched: all 13 occurrences of 'Peclet number' in the 902-page text layer,"
      "\n     both 'Finlayson' hits in Ch. 11 and both in Ch. 12, and the notation list"
      "\n     entries between Pr and PN, which define only Pe_a = u_i d_p / D_ea and"
      "\n     Pe_a' = u_i L / D_ea - never Pe_ha or Pe_ma.)")

CRIT_HEAT_DENOM = T_R0 - T_R0


# %%
# ======== 6. rung 3's interfacial difference against the book's own claim =====
print("section 11.8.2, book page 590, on this same reactor and kinetic scheme:\n   \""
      + bk("film_dT_claim", "as_printed") + "\"")
lad3 = Ladder(RX, rung=3, n_z=N_Z, T_r=bk("startup_T0_1182"))
c3 = solve_robust(lad3, bk("y_A0_1173"))
DT_FILM_1182 = float(np.max(c3[:, 3] - c3[:, 1]))
DP_FILM_1182 = float(np.max(c3[:, 0] - c3[:, 2]))
print(f"\nAt the book's own startup conditions - T_0 = T_r ="
      f" {bk('startup_T0_1182','as_printed')} C, and the 11.7.3\ninlet mole fraction"
      f" y_A0 = {bk('y_A0_1173')} for the same reactor:"
      f"\n   max gas-to-solid temperature difference : {DT_FILM_1182:.4f} K   against the"
      f" printed 'of the order of 1 C'"
      f"\n   max partial-pressure drop across the film: {DP_FILM_1182:.3e} atm,"
      f" {100*DP_FILM_1182/bk('y_A0_1173'):.3f} % of the inlet"
      f"\nThe temperature claim holds; the chapter's other statement that 'the most likely"
      f"\ninterfacial gradient to occur is the temperature gradient' (page 585) is"
      f" confirmed\nby the same pair of numbers."
      f"\n\nTWO OF THE THREE INPUTS ARE NOT SECTION 11.8.2'S. It prints"
      f" T_0 = T_r = {bk('startup_T0_1182','as_printed')} C and\nnothing else quantitative for"
      f" this case: the inlet composition is borrowed from\nsection 11.7.3 and k_g, h_f come"
      f" from A3.4's correlation, which the chapter does not\nprint either. So this"
      f" reproduces 'of the order of 1 C' at the book's own startup\nTEMPERATURE with two"
      f" imported inputs, and it is worth exactly what an order-of-\nmagnitude claim with two"
      f" free inputs is worth. The crossover row below is what\nbounds the film's importance"
      f" without them.")

# Mears' criterion for the onset of interphase temperature gradients, page 588.
# The criterion family belongs to published B1.7; this is it evaluated here.
i_hs = int(np.argmax(c3[:, 1]))
r_obs = kexp(c3[i_hs, 3]) * P["p_B0"] * c3[i_hs, 2] / 3600.0          # kmol/(kg s)
CHI = ((P["minus_dH"] * KCAL) * r_obs * P["rho_b"] * D_P
       / (2.0 * FILM["h_f"] * c3[i_hs, 1]))
MEARS_RHS = bk("mears_constant") / (P["a"] / c3[i_hs, 1])
print(f"\nMears' interphase criterion, book page 588 (the criterion family is B1.7's;"
      f"\nthis is B1.7's eq. 14 evaluated on this reactor):"
      f"\n   |chi| = (-dH) r_A rho_B d_p / (2 h_f T) = {CHI:.5f}"
      f"\n   0.15 R T / E                            = {MEARS_RHS:.5f}"
      f"\n   -> {'SATISFIED' if CHI < MEARS_RHS else 'VIOLATED'}, by a factor"
      f" {MEARS_RHS/CHI:.2f}")
print(f"\nAND YET. At section 11.5.2's own design point the same {FILM_DT_HS:.2f} K film"
      f" difference moves\nthe safe inlet partial pressure by {SHIFT[3]:+.2f} %"
      f" - {abs(SHIFT[3])/(100*bk('runaway_p0_increment')/bk('runaway_p0_safe')):.1f} times"
      f" the {100*bk('runaway_p0_increment')/bk('runaway_p0_safe'):.2f} % margin\nthat section"
      f" 11.5.2's own sensitivity sentence hangs on. Mears' criterion asks whether\nthe"
      f" OBSERVED RATE deviates by 5 % from the intrinsic one, and it correctly says no."
      f"\nIt is not a statement about the runaway boundary, and near a runaway boundary a"
      f"\nrate that is right to 5 % is not enough. That is the ladder's lesson in one"
      f" line:\nthe criteria in the chapter test the RATE, and what a designer of this"
      f" reactor\nneeds tested is the THRESHOLD.")

print(f"\nAND THE CHAPTER SAYS SO ITSELF, in the sentence immediately after the one this"
      f"\npage uses as its anchor. Book page"
      f" {int(BK.loc['heat_mass_transfer_claim','book_page'])}, verbatim:\n   \""
      + bk("heat_mass_transfer_claim", "as_printed") + "\"")
print(f"That is section 11.5.2 stating, at exactly the operating points this page works"
      f"\nat - {P0_BOOK}, and the curves just above it - that the pseudohomogeneous model is"
      f"\nno longer adequate there and that interfacial transport has to be added. Which is"
      f"\nrung 3, qualitatively, PRINTED IN THE CHAPTER. So the contribution here is not"
      f"\nthe observation; it is the size and the sign - {SHIFT[3]:+.2f} % in the threshold,"
      f" against a\nMears criterion that passes by a factor {MEARS_RHS/CHI:.2f} - and the"
      f" point that the chapter's own\nyardsticks are rate criteria while its own warning is"
      f" about the threshold. The\nchapter noticed; it did not quantify, and it did not"
      f" reconcile the warning with\nthe criteria it prints three pages later.")


# %%
# ======== 7. the crossovers, root-found rather than sampled ==================
print("Where does each rung stop mattering? Every row of the table below is a CONTROLLED"
      "\nSCAN: a sweep outward from the chapter's own value in which each threshold must"
      "\npass the crossing control before its shift is used for anything, run until two"
      "\nVALIDATED points straddle the target, then refined between them with every"
      "\niterate controlled too. All four rows are ROOTS.\n"
      "\nThe thresholds follow the start-up branch (see the note on locating thresholds"
      "\nat the top of Validation), and the scan meets no rejected point at all.\n")
BASE600 = p0_critical(Ladder(RX, rung=1, n_z=600))
# ONE definition of "shift" in this table, and it is the same-grid one: rung 1 at the
# SAME n_z as the rung being moved.  Mixing two - some rows against the
# Richardson-extrapolated rung-1 threshold and others against a same-grid one -
# costs a large fraction of the target, because rung 1's own
# 600-cell discretisation error is a fifth of it.  Both readings of the film row are
# measured and printed at the end of this cell.  Same-grid is the right one: it is
# what cancels the discretisation between the two models being compared.


def shift_pct(lad, base=None):
    return 100.0 * (p0_critical(lad, hi=0.09) / (BASE600 if base is None else base)
                    - 1.0)


TARGET = 1.0                      # per cent shift in p0*, the crossover definition
cross, SCANS = [], {}


def _row(key, lbl, quantity, chapter, factor, res):
    """factor = how many times the chapter's own value the located one is."""
    SCANS[key] = res
    cross.append((lbl, quantity, res["kind"], res["value"], chapter, factor,
                  res["shift_pct"]))


_MARGIN_BEFORE_SCANS = CTL_STATS["margin_K"]     # the four scans' own margin, below
CTL_STATS["margin_K"] = np.inf
# (a) rung 2: how small must Pe_ma = Pe_ha be before axial dispersion moves p0* 1 %?
PE_SWEEP = (2.0, 1.5, 1.0, 0.6, 0.4, 0.32, 0.30, 0.25, 0.24, 0.23, 0.22, 0.21,
            0.20, 0.15, 0.10)
S2 = cross_scan(lambda pe: shift_pct(Ladder(RX, rung=2, n_z=600, pe_ma=pe, pe_ha=pe)),
                PE_SWEEP, TARGET)
_row("pe", "rung 2, sec. 11.6", "axial Peclet number Pe_a",
     f"{bk('pe_axial_low'):.0f} to {bk('pe_axial_high'):.0f} (page 560)",
     bk("pe_axial_low") / S2["value"], S2)

# (b) rung 3: how large must the film coefficients be before the film moves p0* 1 %?
FILM_SWEEP = (1.0, 2.0, 3.0, 4.0, 6.0, 10.0, 30.0, 100.0)
S3 = cross_scan(lambda fs: shift_pct(Ladder(RX, rung=3, n_z=600, film_scale=fs)),
                FILM_SWEEP, TARGET)
_row("film", "rung 3, sec. 11.8", "Wakao-Funazkri k_g and h_f, x",
     "1.0 (as correlated)", S3["value"] / 1.0, S3)

# (c) rung 4: how large must D_e be before intraparticle diffusion moves p0* 1 %?
#     measured against RUNG 3, not rung 1: as D_e -> infinity rung 4 collapses onto
#     rung 3 (eta -> 1), not onto rung 1, so rung 1 is the wrong baseline for THIS
#     rung's own contribution.  The other three rows each sit one rung above their
#     own baseline already.
BASE3 = p0_critical(Ladder(RX, rung=3, n_z=600))
DE_SWEEP = (1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4, 1.0e-3, 1.0e-2)
S4 = cross_scan(lambda de: shift_pct(Ladder(RX, rung=4, n_z=600, D_e=de), base=BASE3),
                DE_SWEEP, TARGET)
_row("De", "rung 4, sec. 11.9", "effective diffusivity D_e [m2/s], vs rung 3",
     "NOT PRINTED for this catalyst", S4["value"] / 1.0e-6, S4)

# (d) rung 5: how small must the tube radius be before 2-D moves p0* 1 %?
S5 = cross_R(TARGET)
SCAN_MARGIN_K = CTL_STATS["margin_K"]
CTL_STATS["margin_K"] = min(_MARGIN_BEFORE_SCANS, SCAN_MARGIN_K)
_row("R", "rung 5, sec. 11.7", "tube radius R_t [m]", f"{R_T} (this reactor)",
     R_T / S5["value"], S5)

SWEEP_LABEL = {"pe": "Pe_ma = Pe_ha", "film": "film coefficients, x",
               "De": "D_e [m2/s]", "R": "R_t [m]"}
for _k, _res in SCANS.items():
    _t = _res["sweep"].rename(columns={"x": SWEEP_LABEL[_k]})
    print(f"\nthe {SWEEP_LABEL[_k]} sweep: {_res['n_evaluated']} points evaluated,"
          f" {_res['n_controlled']} validated, the search uses only those,\nand the sweep"
          f" stops once two VALIDATED points straddle the target. Refinement:"
          f"\n{_res['levels']} levels took the {_res['bracket0']:.4g}-wide straddling pair"
          f" to {_res['bracket']:.3g}.")
    show(_t.style.set_uuid("d11sweep" + _k).hide(axis="index").format(
        {SWEEP_LABEL[_k]: "{:.6g}", "shift_pct": "{:+.4f}"},
        na_rep="the control REJECTS this point"))
PE_STAR, FILM_STAR, DE_STAR, R_STAR = (S2["value"], S3["value"], S4["value"],
                                       S5["value"])
CROSS = pd.DataFrame(cross, columns=["rung", "quantity", "kind", "crossover",
                                     "the chapter's value", "factor away",
                                     "shift there, %"])
# NOT asserted that list(CROSS.kind) == ["root"] * 4: cross_scan has a single
# `return`, always kind="root", so that assert could never fail - a toothless
# check on the one page that audits toothless checks.  The "kind" column stays,
# printed rather than guarded.
show(CROSS.style.set_uuid("d11cross").hide(axis="index").format(
    {"crossover": "{:.4g}", "factor away": "{:.3g}", "shift there, %": "{:+.4f}"}))

# the certificate the located rung-2 crossover itself carries, printed rather than
# asserted - it is the whole of the case that this row is a root and not an artefact
_L2S = Ladder(RX, rung=2, n_z=600, pe_ma=PE_STAR, pe_ha=PE_STAR)
_P2S = p0_critical(_L2S, hi=0.09)
_LO2, _HI2, _OK2 = crossing_control(_L2S, _P2S)
assert _OK2, "the located rung-2 crossover does not pass the crossing control"
SCAN_PROBED = sum(_r["probed"] for _r in SCANS.values())
SCAN_REFUSED = (sum(_r["refused"] for _r in SCANS.values())
                + sum(_r["n_evaluated"] - _r["n_controlled"] for _r in SCANS.values()))
print(f"\nALL FOUR ROWS ARE ROOTS. At Pe_a = {PE_STAR:.4f} the rung-2 shift is"
      f" {S2['shift_pct']:+.5f} %, and that\nthreshold carries the same certificate as the"
      f" five rungs of the ladder table -\np0* = {_P2S:.8f} atm with {_LO2-273.15:.2f} C just"
      f" below it and {_HI2-273.15:.2f} C just above, a\n{_HI2-_LO2:.3f} K jump against the"
      f" {CTL_JUMP:.0f} K the control allows. So axial dispersion has to be"
      f"\n{bk('pe_axial_low')/PE_STAR:.2f} times stronger than the chapter's own lower limit"
      f" before it moves the\nthreshold {TARGET:.0f} %. The four refinements probed"
      f" {SCAN_PROBED} interior points between them and,\nwith the sweep points, the control"
      f" rejected {SCAN_REFUSED} of all the thresholds the four scans used.")

RANK = CROSS.sort_values("factor away")
print(f"\nRead the 'factor away' column as a RANKING: it is how far each coefficient would"
      f"\nhave to move from the chapter's own value before its rung shifted the threshold"
      f"\n{TARGET:.0f} %, so SMALL means the rung is close to not mattering.")
for _, _row_ in RANK.iterrows():
    print(f"   {_row_['rung']:<20s} factor {_row_['factor away']:7.2f}"
          f"   {_row_['quantity']}")
print(f"\nThe ordering is DECIDABLE, all four of it: rung 5 is closest to irrelevance by a"
      f"\nwide margin, then rung 2 at {bk('pe_axial_low')/PE_STAR:.2f}, then rung 3 at"
      f" {FILM_STAR:.2f}, and rung 4 is furthest at\n{DE_STAR/1.0e-6:.2f}. Rungs 2 and 3 are"
      f" only {100*abs(FILM_STAR/(bk('pe_axial_low')/PE_STAR)-1):.0f} % apart, which is the"
      f" honest caveat on THAT pair.")

# THE CONTROL'S VERDICT ALONG Pe_a, on a regular 0.001 grid from 0.316 to 0.308.
_UNIFORM = []
for _p in np.round(np.arange(0.316, 0.3075, -0.001), 4):
    _UNIFORM.append((float(_p), *_shift_or_none(
        lambda p: shift_pct(Ladder(RX, rung=2, n_z=600, pe_ma=p, pe_ha=p)), float(_p))))
UNIFORM_OK = sum(1 for _, _, o in _UNIFORM if o)
assert UNIFORM_OK == len(_UNIFORM), "the control rejects a point of the fine Pe_a grid"
assert np.all(np.diff([v for _, v, _ in _UNIFORM]) > 0), (
    "the rung-2 shift is not monotone on the fine Pe_a grid")
print(f"\nTHE CONTROL'S VERDICT ALONG Pe_a, on a regular 0.001 grid:\n   "
      + "  ".join(f"{p:.3f} {v:+.4f} %" for p, v, _ in _UNIFORM)
      + f"\n{UNIFORM_OK} of {len(_UNIFORM)} validate, and the shift rises monotonically as"
      f" Pe_a falls.")

# THE REST OF THE Pe_a SWEEP, evaluated although the crossover did not need it. The
# scan above stops as soon as the target is bracketed, so "the largest shift the sweep
# shows" would otherwise be an artefact of where it stopped.  The tail is measured.
PE_TAIL = [(float(_p), *_shift_or_none(
    lambda x: shift_pct(Ladder(RX, rung=2, n_z=600, pe_ma=x, pe_ha=x)), float(_p)))
    for _p in PE_SWEEP[S2["n_evaluated"]:]]
_ALL2 = ([(float(r.x), float(r.shift_pct)) for r in S2["sweep"].itertuples() if r.controlled]
         + [(p, v) for p, v, ok in PE_TAIL if ok])
_MAX2 = max(_ALL2, key=lambda t: abs(t[1]))
PE_TAIL_OK = sum(1 for _, _, ok in PE_TAIL if ok)
print(f"\nWHAT THE WHOLE SWEEP REACHES, tail included:\n   "
      + "  ".join(f"{p:.2f} " + (f"{v:+.4f} %" if ok else "REFUSED") for p, v, ok in PE_TAIL)
      + f"\nAcross all {len(PE_SWEEP)} swept points the control validates {len(_ALL2)}, and the"
      f" largest |shift| any\nof them carries is {abs(_MAX2[1]):.4f} % at"
      f" Pe_a = {_MAX2[0]:.2f}. So a {2*TARGET:.0f} % crossover for rung 2 is bracketable on"
      f"\nthis sweep, and the break row that moves the crossover target moves it to"
      f" {2*TARGET:.0f} %.\nWhat lies below Pe_a = {min(p for p, _ in _ALL2):.2f} is not"
      f" reported here.")

# THE ACCURACY OF THE RUNG-2 ROOT, the same way rung 5's is done below: the same-grid
# shift at a validated sweep point on three grids, divided by the local slope.
PE_ANCHOR = S2["pair"][1]
PE_GRID = {}
for _n in (600, 900, N_Z):
    _b1 = BASE600 if _n == 600 else p0_critical(Ladder(RX, rung=1, n_z=_n))
    PE_GRID[_n] = 100.0 * (p0_critical(Ladder(RX, rung=2, n_z=_n, pe_ma=PE_ANCHOR,
                                              pe_ha=PE_ANCHOR), hi=0.09) / _b1 - 1.0)
PE_GAP = max(PE_GRID.values()) - min(PE_GRID.values())
PE_SLOPE = abs(S2["pair_shift"][1] - S2["pair_shift"][0]) / abs(S2["pair"][1]
                                                                - S2["pair"][0])
PE_STAR_UNC = PE_GAP / PE_SLOPE
print(f"\nAND THE DIGITS OF Pe_a* THAT MEAN ANYTHING. The refinement resolves the"
      f" straddling\npair to {S2['bracket']:.1e}, but that is a TOLERANCE. The accuracy is"
      f" set by the grid, and it is\nmeasured under the SAME-GRID definition this table uses"
      f" - which is the point, because\nagainst a fixed extrapolated baseline the same"
      f" comparison would carry rung 1's own\n600-to-{N_Z} movement instead of cancelling it."
      f" At Pe_a = {PE_ANCHOR:.4f}, a validated sweep point:"
      + "".join(f"\n   n_z = {n:5d}   shift = {v:+.5f} %" for n, v in PE_GRID.items())
      + f"\na spread of {PE_GAP:.4f} pp across a doubling of the grid. On the local slope of"
      f" the\nsweep, {PE_SLOPE*1e-3:.4f} pp per 0.001 in Pe_a, that is"
      f" Pe_a* = {PE_STAR:.4f} +/- {PE_STAR_UNC:.4f},\nfactor"
      f" {bk('pe_axial_low')/PE_STAR:.3f} +/-"
      f" {bk('pe_axial_low')*PE_STAR_UNC/PE_STAR**2:.3f}. Three figures on Pe_a*, not the"
      f" six its refinement\ntolerance would suggest - the same reading the tube radius gets"
      f" below.")
_sw5 = S5["sweep"]
_brk5 = abs(S5["pair"][1] - S5["pair"][0])     # the straddling pair the root sits in
_slope5 = abs(abs(S5["pair_shift"][1]) - abs(S5["pair_shift"][0])) / _brk5
_gap5 = abs(_sw5.shift_pct[0] - SHIFT[5])     # (600, 10) against (1200, 20), R_t as built
R_STAR_UNC = _gap5 / _slope5
assert np.all(np.diff(_sw5.shift_pct[_sw5.controlled].values) > 0), (
    "the rung-5 shift is not monotone over the validated sweep")
print(f"\nRUNG 5, THE ONE ROW WHOSE CROSSOVER IS INSIDE THIS REACTOR'S OWN RANGE."
      f" R_t* = {R_STAR:.6f} m,\nd_t/d_p = {2*R_STAR/D_P:.3f}, where |shift| ="
      f" {abs(S5['shift_pct']):.4f} %: a tube only {R_T/R_STAR:.2f}x narrower than this"
      f"\none - {R_STAR*1e3:.2f} mm rather than {R_T*1e3:.1f} mm - already brings the"
      f" two-dimensional rung\ninside {TARGET:.0f} %. THE DIGITS THAT MEAN ANYTHING:"
      f" the refinement resolves the bracket to\n{S5['bracket']*1e6:.1f}"
      f" micrometres, but that is a TOLERANCE and not an accuracy. The accuracy is set"
      f"\nby the grid: this page's own (600, 10) and ({N_Z}, {N_R}) values of the rung-5"
      f" shift differ by\n{_gap5:.3f} pp, and the local slope of the sweep is"
      f" {_slope5*1e-3:.1f} pp/mm, so R_t* = {R_STAR*1e3:.3f} +/- {R_STAR_UNC*1e3:.3f} mm,"
      f"\nd_t/d_p = {2*R_STAR/D_P:.2f} +/- {2*R_STAR_UNC/D_P:.2f}, factor"
      f" {R_T/R_STAR:.3f} +/- {R_T*R_STAR_UNC/R_STAR**2:.3f}. Three figures on R_t*, not six.")

# NARROWER TUBES THAN THE CROSSOVER, measured HERE and not inside the crossover scan:
# the scan stops as soon as the target is bracketed, and this is a different question.
R_NARROW = (0.0098, 0.00975, 0.0097, 0.0096, 0.0094, 0.0090)
_NARROW = [(R, *_shift_or_none(shift_rung5, R)) for R in R_NARROW]
assert all(ok for _, _, ok in _NARROW), "the control rejects a radius below the crossover"
_N_OK5 = [(R, v) for R, v, ok in _NARROW]
assert np.all(np.diff([v for _, v in _N_OK5]) > 0), (
    "the rung-5 shift is not monotone over the radii below the crossover")
_SIGN5 = [(R, v) for R, v in _N_OK5 if v > 0.0][0]
_NEG5 = [(R, v) for R, v in _N_OK5 if v < 0.0][-1]
assert _SIGN5[0] < _NEG5[0], "the rung-5 sign change is not where this text says"
_UNDER5 = [(R, v) for R, v in _N_OK5 if abs(v) < TARGET]
_OVER5 = [(R, v) for R, v in _N_OK5 if v > TARGET]
assert _OVER5 and min(R for R, _ in _UNDER5) > max(R for R, _ in _OVER5), (
    "the second 1 % crossing on narrow tubes is not where this text says")
print(f"\nAND ON TUBES NARROWER THAN THE CROSSOVER. On {len(R_NARROW)} radii below it the"
      f" control validates\nevery threshold:\n   "
      + "  ".join(f"{R*1e3:.2f} mm {v:+.4f} %" for R, v in _N_OK5)
      + f"\nThe signed shift rises monotonically as the tube narrows (checked, not assumed)"
      f" and\nCHANGES SIGN between {_NEG5[0]*1e3:.2f} mm ({_NEG5[1]:+.4f} %) and"
      f" {_SIGN5[0]*1e3:.2f} mm ({_SIGN5[1]:+.4f} %): the\ntwo-dimensional correction is"
      f" NEGATIVE on this tube and POSITIVE on a narrower one.\n|shift| stays under"
      f" {TARGET:.0f} % from R_t* = {R_STAR*1e3:.2f} mm down to"
      f" {min(R for R, _ in _UNDER5)*1e3:.2f} mm, and is past {TARGET:.0f} % again,"
      f"\non the positive side, by {max(R for R, _ in _OVER5)*1e3:.2f} mm"
      f" ({[v for R, v in _OVER5][0]:+.4f} %).")
print(f"\nNone of that contradicts section 11.7.3's own sentence, which reads, verbatim:"
      f"\n   \"" + bk("radial_gradient_claim", "as_printed") + "\"" + f"\nThat sentence is"
      f" about the radial temperature GRADIENT and this measure is about the"
      f"\nrunaway THRESHOLD: the gradient here is"
      f" {float(R5[(R5.n_z==N_Z)&(R5.n_r==N_R)].dT_radial.iloc[0]):.1f} K wall-to-axis while"
      f" the threshold moves {abs(SHIFT[5]):.1f} %. Severe\ngradients and a modest threshold"
      f" shift are the same reactor."
      f"\nThe film sits a factor {FILM_STAR:.1f} away, which is why that rung's size depends"
      f"\non which correlation you pick. The intraparticle rung is the one that is not"
      f"\ndecidable from the chapter at all: the crossover D_e is {DE_STAR:.2e} m2/s,"
      f"\nand the chapter prints no D_e for this catalyst - it sends the reader to two"
      f"\ndocuments ({bk('De_source_claim','as_printed')[:52]}...).")

# THE CONTROL ON A COARSE GRID: rung 4's anchor at n_z = 300, the grid of the
# coarse-grid break row.
_L300 = Ladder(RX, rung=4, n_z=300)
_P300 = p0_critical(_L300, hi=0.09)
_LO3, _HI3, _OK3 = crossing_control(_L300, _P300)
print(f"\nAND ON A COARSE GRID. At n_z = 300, rung 4's own anchor point - its default"
      f"\nD_e = {1.0e-6:.0e} m2/s, the first point of that sweep - returns p0* = {_P300:.10f}"
      f" with\n{_LO3-273.15:.3f} C just below it and {_HI3-273.15:.3f} C just above: a clean"
      f" crossing ({_OK3}), so the\ncoarse-grid break row runs at 300 cells.")

# WHAT THE OTHER SHIFT DEFINITION WOULD COST, measured on the row where it is largest
S3_EXTRAP = cross_scan(lambda fs: shift_pct(Ladder(RX, rung=3, n_z=600, film_scale=fs),
                                            base=LIM_P), FILM_SWEEP, TARGET)
print(f"\nTHE COST OF THE SHIFT DEFINITION, on the row where it is largest. Rung 1's"
      f" 600-cell\nthreshold is {BASE600:.8f} atm and its Richardson limit is {LIM_P:.8f},"
      f" a gap of\n{100*(BASE600/LIM_P-1):+.3f} percentage points - a fifth of the"
      f" {TARGET:.0f} % target. The film row against the"
      f"\nsame-grid baseline used above: {S3['value']:.4f}. Against the extrapolated"
      f" baseline: {S3_EXTRAP['value']:.4f}.\nOne definition, said once, in every row: the"
      f" same-grid one, because it is what\ncancels the discretisation between the two models"
      f" being compared.")

# what the printed a_1 group of page 574 would cost, given the chapter's own
# insensitivity claim for Pe_mr
print(f"\nAnd the printed a_1 = (D_er)_s/(u_i d_p) inconsistency of page 574, measured."
      f"\nThe chapter says: \"{bk('pe_mr_insensitive_claim','as_printed')}\"")
S_PEMR_AS_IS = shift_pct(Ladder(RX, rung=5, n_z=600, n_r=10))
S_PEMR_LITERAL = shift_pct(Ladder(RX, rung=5, n_z=600, n_r=10, pe_mr=PE_MR / EPS))
S_LAMER_5PCT = shift_pct(Ladder(RX, rung=5, n_z=600, n_r=10, lam_er=LAM_ER * 1.05,
                                pe_hr=PE_HR / 1.05))
print(f"   Pe_mr = {PE_MR:.0f} as this page reads the group : p0* shift"
      f" {S_PEMR_AS_IS:+.3f} %"
      f"\n   Pe_mr = {PE_MR/EPS:.1f}, the printed group read literally: p0* shift"
      f" {S_PEMR_LITERAL:+.3f} %  (delta"
      f" {S_PEMR_LITERAL-S_PEMR_AS_IS:+.3f} %)"
      f"\n   lambda_er raised 5 %                          : p0* shift"
      f" {S_LAMER_5PCT:+.3f} %  (delta {S_LAMER_5PCT-S_PEMR_AS_IS:+.3f} %)"
      f"\nThe chapter's sensitivity claim is confirmed on the threshold as well as on the"
      f"\nprofiles, and it is WHY the epsilon in a_1 is invisible: a factor 1/eps ="
      f" {1/EPS:.2f} in\nradial mass dispersion costs"
      f" {abs(S_PEMR_LITERAL-S_PEMR_AS_IS):.3f} % where 5 % on lambda_er costs"
      f" {abs(S_LAMER_5PCT-S_PEMR_AS_IS):.3f} %.\nA printed group can be wrong by a void"
      f" fraction and never be found, if the quantity\nit multiplies is the one the model is"
      f" insensitive to.")


# %%
# ===== 8. the book's restatement of Froment (1967), against C2.10's reading ====
# C2.10 owns the 1967 case. Nothing of it is re-run here; its numbers are printed
# beside the book's restatement of them, and the two are reconciled.
F67 = dict(zip(F67_PAR.index, F67_PAR.value))
rec = []
for lab, bkey, f67key, factor in (
        ("E_1/R", "E1_over_R_1173", None, None),
        ("E_2/R", "E2_over_R_1173", None, None),
        ("E_3/R", "E3_over_R_1173", None, None)):
    e = {"E_1/R": "E1", "E_2/R": "E2", "E_3/R": "E3"}[lab]
    rec.append(dict(quantity=lab, book_11_7_3=bk(bkey),
                    from_C2_10=F67[e] / F67["R_gas"],
                    ratio=bk(bkey) / (F67[e] / F67["R_gas"])))
rec.append(dict(quantity="E_1/R used in sec. 11.5.2", book_11_7_3=bk("rate_E_over_R"),
                from_C2_10=F67["E1"] / F67["R_gas"],
                ratio=bk("rate_E_over_R") / (F67["E1"] / F67["R_gas"])))
REC = pd.DataFrame(rec)
show(REC.style.set_uuid("d11rec").hide(axis="index").format(
    {"book_11_7_3": "{:.1f}", "from_C2_10": "{:.2f}", "ratio": "{:.6f}"}))
ER_INTERNAL = bk("rate_E_over_R") / bk("E1_over_R_1173") - 1.0
print(f"C2.10 establishes that the 1967 paper's own gas constant is"
      f" {F67['R_gas']} cal/(mol K), and\n{F67['E1']:.0f}/{F67['R_gas']} ="
      f" {F67['E1']/F67['R_gas']:.1f} is what section 11.5.2 prints. Section 11.7.3, fourteen"
      f"\npages earlier in the same chapter and for the same reaction on the same catalyst,"
      f"\nprints {bk('E1_over_R_1173','as_printed')}, which is {F67['E1']:.0f}/1.987."
      f" THE CHAPTER PRINTS THE SAME ACTIVATION\nENERGY TWO WAYS, {100*ER_INTERNAL:+.3f} %"
      f" apart, and neither section mentions the other's value.")

# what section 11.5.2's "first approximation" actually is
def k_1967(T, i):
    return np.exp(F67[f"lnA{i}"] - F67[f"E{i}"] / (F67["R_gas"] * T))


K1_HS, K3_HS = k_1967(T_HS, 1), k_1967(T_HS, 3)
assert P["minus_dH"] == F67["minus_dH1"], "11.5.2's heat of reaction is not C2.10's minus_dH1"
K_1152_HS = kexp(T_HS)
LUMP_DEFICIT = K_1152_HS / (K1_HS + K3_HS) - 1.0
print(f"\nAnd what the 'first approximation ... pseudo-first order' rate of section 11.5.2"
      f"\nIS, which the section does not say. At the rung-1 hot spot,"
      f" T = {T_HS-273.15:.1f} C:"
      f"\n   k of section 11.5.2      = {K_1152_HS:.6f}"
      f"\n   k_1 of section 11.7.3    = {K1_HS:.6f}     (o-xylene -> phthalic anhydride)"
      f"\n   k_3 of section 11.7.3    = {K3_HS:.6f}     (o-xylene -> CO, CO2)"
      f"\n   k_1 + k_3                = {K1_HS+K3_HS:.6f}     (total o-xylene consumption)"
      f"\nSection 11.5.2's single rate is k_1 alone to {abs(K_1152_HS/K1_HS-1):.1e}, not the"
      f" sum: it is the\nrate of the WANTED reaction used as the rate of total consumption,"
      f" {100*LUMP_DEFICIT:+.1f} % low.\nIts heat of reaction, {P['minus_dH']:,.0f} kcal/kmol, is"
      f" C2.10's minus_dH1 = {F67['minus_dH1']:.0f} exactly, so the\nheat release is"
      f" understated by the same factor and the {F67['minus_dH3']:.0f} kcal/kmol of the"
      f"\ncombustion route is absent. The ladder's base case is therefore a consistent"
      f"\nsingle-reaction reactor, and NOT a conservative approximation to the three-"
      f"\nreaction one - a point that matters for reading any of the shifts above as"
      f"\nadvice about a real o-xylene reactor.")

print(f"\nFinally, section 11.7.4's own one- versus two-dimensional comparison, book page"
      f" 579:\n   \"" + bk("equiv1d_claim", "as_printed") + "\"")
c210_2d = F67_RES[F67_RES.quantity.str.contains("runaway", case=False, na=False)]
show(c210_2d.style.set_uuid("d11c210").hide(axis="index"))
print(f"The book states {bk('runaway_1d_1174','as_printed')} C for the equivalent"
      f" one-dimensional model and"
      f" '{bk('runaway_2d_1174','as_printed')}' C for the\ntwo-dimensional one - a sharpened"
      f" pair. C2.10 located both\nby bisection on the 1967"
      f" paper's own kinetics and reports them above. Read the two\ntogether: the book turns a"
      f" range into a number. This page does not re-run that\ncase, and its own rung-1/rung-5"
      f" pair - {SHIFT[5]:+.2f} % in p0* on the pseudo-first-order\nchemistry - is a different"
      f" measurement of the same effect, in the same direction.")


# %%
# ============== defect injection: every metric needs a row that moves it ======
# The helpers below are the ONLY route to each baseline metric, so a row that
# perturbs a helper's argument provably reaches the code that produced the number.
# Each helper's default call is asserted against the baseline computed above.

def h_ex1151a(dp=None, a=None, b=None):
    """Calls the SAME f_ergun / f_hicks / dP the table above was built from, so a
    row that perturbs dp, a or b provably reaches the reported numbers."""
    dp = bk("ex1151a_dp") if dp is None else dp
    a = bk("ergun_a") if a is None else a
    b = bk("ergun_b") if b is None else b
    fe = f_ergun(dp, a, b)[1]
    fm = f_ergun(dp, bk("mcdonald_a"), bk("mcdonald_b"))[1]
    fh = f_hicks(dp)[1]
    d = dict(Ergun=dP(fe, dp) / bk("ex1151a_dp_ergun") - 1.0,
             McD=dP(fm, dp) / bk("ex1151a_dp_mcdonald") - 1.0,
             Hicks=dP(fh, dp) / 1e5 / bk("ex1151a_dp_hicks") - 1.0)
    ridx = dict(
        Ergun=brentq(lambda R: (1 - E["eps"]) / E["eps"] ** 3 * (
            a + b * (1 - E["eps"]) / R) - bk("ex1151a_f_ergun"), 20, 900),
        McD=brentq(lambda R: (1 - E["eps"]) / E["eps"] ** 3 * (
            bk("mcdonald_a") + bk("mcdonald_b") * (1 - E["eps"]) / R)
            - bk("ex1151a_f_mcdonald"), 20, 900),
        Hicks=brentq(lambda R: bk("hicks_coefficient") * (1 - E["eps"]) ** bk(
            "hicks_exponent_eps") / E["eps"] ** 3 * R ** bk("hicks_exponent_Re")
            - bk("ex1151a_f_hicks"), 20, 900))
    return dict(ex1151a_worst_dev=max(abs(v) for v in d.values()),
                dp_definition_cost=(dP(f_ergun(DP_SV, a, b)[1], DP_SV)
                                    / dP(f_ergun(DP_ROUNDED, a, b)[1], DP_ROUNDED)
                                    - 1.0),
                re_implied_spread=max(ridx.values()) - min(ridx.values()))


def h_units(p_B0_book=None, atm=None):
    p_B0_book = bk("oxygen_pressure") if p_B0_book is None else p_B0_book
    atm = ATM if atm is None else atm
    worst = max(abs(p_B0_book / (V["p_B0"] * atm) - 1.0),
                abs(bk("heat_of_reaction") / (P["minus_dH"] * KCAL) - 1.0))
    for q in ("p0_lower", "p0_upper", "p0_mean"):
        o = float(VWF_EX[(VWF_EX.example == "1a") & (VWF_EX.quantity == q)].value.iloc[0])
        worst = max(worst, abs(bk("vwf_p0_" + q.split("_")[1] + "_625") / (o * atm) - 1.0))
    o = float(VWF_EX[(VWF_EX.example == "1a") & (VWF_EX.quantity == "p0_critical")].value.iloc[0])
    return dict(conversion_worst_dev=worst,
                backint_residual=bk("vwf_p0_backint_625") / (o * atm) - 1.0)


def h_groups(cp_vol=None, minus_dH=None, p_B0=None):
    cp_vol = P["cp_vol"] if cp_vol is None else cp_vol
    minus_dH = P["minus_dH"] if minus_dH is None else minus_dH
    p_B0 = P["p_B0"] if p_B0 is None else p_B0
    Bc = minus_dH * P["rho_b"] * p_B0 / cp_vol
    Ac = P["M"] * P["p_t"] * P["rho_b"] / P["rho_g"] * p_B0
    cps = [cp_vol * KCAL / P["rho_g"],
           minus_dH * P["rho_b"] * p_B0 / bk("vwf_B") * KCAL / P["rho_g"],
           PE_HR * LAM_ER / (G_FLUX * D_P),
           F67_PAR.loc["c_p", "value"] * KCAL]
    return dict(group_B_dev=Bc / bk("vwf_B") - 1.0,
                group_A_printed_dev=Ac / bk("vwf_A") - 1.0,
                cp_printed_dev=float(np.mean(cps)) / bk("specific_heat") - 1.0,
                cp_route_spread=float(max(cps) / min(cps) - 1.0))


def h_equiv1d(R=None, alpha_w=None, lam_er=None):
    R = R_T if R is None else R
    alpha_w = ALPHA_W if alpha_w is None else alpha_w
    lam_er = LAM_ER if lam_er is None else lam_er
    ai = 1.0 / (1.0 / alpha_w + R / (4.0 * lam_er))
    ai_printed_R = 1.0 / (1.0 / alpha_w + (D_T_PRINTED / 2.0) / (4.0 * lam_er))
    # U_from_2d_dev lands at machine zero - the book's 0.096 IS eq. (11.7.4-1) at
    # R_t = 0.0125 m to every digit it prints - so it is registered as a bound.  The
    # companion is the same identity evaluated at the radius the printed tube
    # diameter implies, a value well away from zero.
    return dict(U_from_2d_dev=ai / bk("overall_U") - 1.0,
                U_from_2d_dev_at_printed_radius=ai_printed_R / bk("overall_U") - 1.0)


def h_tables(Q=None, dT=None, Q1=None, TM1=None, A=None, B=None, atm=None):
    """Table 11.5.3.A-1 and A-2, INCLUDING the cross-column test that a row-by-row
    transcription check cannot do: dT_ad = (B/A) p0 ties the two columns of A-2
    together, so A, B and the atm factor are arguments a break row can move."""
    Q = bk("vwf_Q_625") if Q is None else Q
    dT = bk("vwf_dT_625") if dT is None else dT
    Q1 = bk("vwf_Q_635") if Q1 is None else Q1
    DT1 = (bk("vwf_TM_635") if TM1 is None else TM1) - 635.0
    A = bk("vwf_A") if A is None else A
    B = bk("vwf_B") if B is None else B
    atm = ATM if atm is None else atm
    w = max(abs(DT1 * (1 + Q1 ** 2) / bk("vwf_dTad_lower_635") - 1.0),
            abs(DT1 * (1 + Q1) ** 2 / bk("vwf_dTad_upper_635") - 1.0),
            abs(DT1 * (1 + Q1 + Q1 ** 2) / bk("vwf_dTad_mean_635") - 1.0))
    boa, devs = B / A, []
    for pk, dk, form in (("vwf_p0_lower_625", "vwf_dTad_lower_625",
                          lambda q: dT * (1 + q ** 2)),
                         ("vwf_p0_upper_625", "vwf_dTad_upper_625",
                          lambda q: dT * (1 + q) ** 2),
                         ("vwf_p0_mean_625", "vwf_dTad_mean_625",
                          lambda q: dT * (1 + q + q ** 2))):
        needed = bk(pk) / atm * boa
        devs.append((bk(dk) / needed - 1.0, needed, form))
    q_imp = brentq(lambda q: devs[0][2](q) - devs[0][1], 0.1, 20.0, xtol=1e-13)
    return dict(table2_upper_ratio=bk("vwf_dTad_upper_625") / (dT * (1 + Q) ** 2),
                table1_worst_dev=w,
                table2_dTad_worst_dev=max(abs(d[0]) for d in devs),
                table2_Q_implied_by_p0=q_imp)


def h_eta(n_r=None, nu=None, phi=None, eta_col=None):
    n_r = 256 if n_r is None else n_r
    phi = PHI_HS if phi is None else phi
    nu = 2 if nu is None else nu                # another nu injects a wrong geometry
    e = Pellet(RX, n_r=n_r, nu=nu).solve(THETA * P_HS, T_HS)["eta"]
    ana = eta_sphere(phi)
    col = SPH.eta.values if eta_col is None else eta_col
    return dict(eta_grid_vs_11_9_1_11=abs(e / ana - 1.0),
                eta_vs_B11_worst=float(np.abs(eta_sphere(SPH.phi.values) / col - 1.0).max()))


def h_shifts(n_z=600, ceiling=None, **kw):
    """the four rung-to-rung shifts, on the coarse grid the break table uses"""
    ceiling = T_PERM if ceiling is None else ceiling
    base = p0_critical(Ladder(RX, rung=1, n_z=n_z, **{k: v for k, v in kw.items()
                                                  if k in ("R", "T_r")}),
                       ceiling=ceiling)
    out = {}
    for r, key, hi in ((2, "shift_rung2_pct", 0.02), (3, "shift_rung3_pct", 0.02),
                       (4, "shift_rung4_pct", 0.09), (5, "shift_rung5_pct", 0.02)):
        allowed = {k: v for k, v in kw.items() if k in Ladder.__init__.__code__.co_varnames}
        lad = Ladder(RX, rung=r, n_z=n_z, n_r=10, **allowed)
        out[key] = 100.0 * (p0_critical(lad, ceiling=ceiling, hi=hi) / base - 1.0)
    return out



def h_rung1(R=None, ceiling=None, ns=NS_RUNG1, extrapolate=True, ratio_read=None):
    """rung 1's seven metrics, from THE SAME `rung1_chain` call the Validation cell
    displayed - the default call is memoised, so these are the displayed numbers
    themselves rather than a second, coarser computation of them.  A row that
    perturbs R or the ceiling re-runs the full chain.  `extrapolate=False` compares
    the finest grid with the ODE route instead of the Richardson limit (so the
    finest grid is also the limit the order-fit guard measures against), and
    `ratio_read` fits the orders with the grid ratio misread; both reuse the chain."""
    out = rung1_chain(R=R, ceiling=ceiling, ns=ns)
    m = dict(out["metrics"])
    tab = out["table"]
    if not extrapolate:
        m["rung1_p0_star_extrapolated"] = float(tab.p0_star.iloc[-1])
        m["ode_vs_pymrm_p0"] = abs(float(tab.p0_star.iloc[-1]) / out["p0_ode"] - 1.0)
        m["ode_vs_pymrm_T"] = abs(float(tab.T_hot.iloc[-1]) / out["T_ode"] - 1.0)
        # the guard measured against the limit in use, which is now the finest grid
        m["rung1_finest_err_over_start_floor"] = abs(
            float(tab.p0_star.iloc[-1]) / m["rung1_p0_star_extrapolated"] - 1.0) / START_FLOOR
    if ratio_read is not None:
        ns_read = [ns[0] * ratio_read ** i for i in range(len(ns))]
        m["rung1_order_p0"] = richardson(tab.p0_star, ns_read)[0]
        m["rung1_order_T"] = richardson(tab.T_hot, ns_read)[0]
    return m


def h_units_model(p_B0=None, p_t=None, ceiling=None):
    """What section 11.5.2's sensitivity sentence gives on a stated reading of its
    own pressure unit.  Grid-free (Radau) on both readings, so a break row that
    flips the reading moves ONLY the parameter set."""
    p_B0 = bk("oxygen_pressure") if p_B0 is None else p_B0
    p_t = bk("total_pressure") if p_t is None else p_t
    ceiling = T_PERM if ceiling is None else ceiling
    u = unit_reading(p_B0, p_t, ceiling=ceiling)
    return dict(p0_star_bar_reading_dev=u["p0_star"] / P0_BOOK - 1.0,
                T_hot_0018_bar_reading_C=u["T_at_book_p0"] - 273.15)


_H_CROSS_DIAG = {}     # (target, n_z) -> the three cross_scan() dicts, kept rather
                       # than discarded, so a break-row scan's own refinement can be
                       # inspected the same way the reported rows' is - see the
                       # break-row-scan diagnostics printed after BREAKS below.


def h_cross(target=None, n_z=600):
    """The three one-dimensional crossover rows, from the same CONTROLLED SCAN as
    the display cell: same sweeps, same same-grid baseline, same bisection, so the
    two are identical by construction at the default arguments and are asserted so.
    A row that moves the target or the grid re-runs the whole scan, control and all;
    a scan that cannot straddle its target between two validated points RAISES
    rather than falling back on a bound, so no row can report an edge as an answer."""
    target = TARGET if target is None else target
    base = p0_critical(Ladder(RX, rung=1, n_z=n_z))
    def sh(lad, b=None):
        return 100.0 * (p0_critical(lad, hi=0.09) / (base if b is None else b) - 1.0)
    pe = cross_scan(lambda p: sh(Ladder(RX, rung=2, n_z=n_z, pe_ma=p, pe_ha=p)),
                    PE_SWEEP, target)
    fs = cross_scan(lambda f: sh(Ladder(RX, rung=3, n_z=n_z, film_scale=f)),
                    FILM_SWEEP, target)
    b3 = p0_critical(Ladder(RX, rung=3, n_z=n_z))       # rung 4's own baseline
    de = cross_scan(lambda d: sh(Ladder(RX, rung=4, n_z=n_z, D_e=d), b3), DE_SWEEP, target)
    _H_CROSS_DIAG[(target, n_z)] = dict(pe=pe, fs=fs, de=de)
    return dict(crossover_pe_axial=pe["value"], crossover_film_scale=fs["value"],
                crossover_D_e=de["value"])


def h_cross_R(target=None, n_z=600, n_r=10):
    """The tube-radius crossover, on the same sweep and the same grids as the
    display cell.  Memoised there, so a default call costs nothing twice; each
    evaluation is two full threshold bisections, one of them two-dimensional, which
    is what limits the refinement to four levels of the straddling pair - a
    TOLERANCE of a couple of micrometres against a grid-set ACCURACY of some
    3e-5 m, both printed above."""
    target = TARGET if target is None else target
    return dict(crossover_R_t=cross_R(target, n_z=n_z, n_r=n_r)["value"])


def h_film_dt(film_scale=None, n_r_pellet=None):
    """The three film/pellet temperature differences and Mears' margin, recomputed
    from the same solves the Validation cells displayed.

    They are covered by the film-scale and coarse-pellet rows."""
    fs = 1.0 if film_scale is None else film_scale
    n_p = 128 if n_r_pellet is None else n_r_pellet
    lad3 = Ladder(RX, rung=3, n_z=N_Z, n_r=N_R, film_scale=fs)
    # AT RUNG 3'S OWN THRESHOLD, recomputed for the film coefficients in force -
    # which is what the displayed number means, and what makes the helper safe to
    # perturb: at HALVED film coefficients the unperturbed threshold is past this
    # model's own, and the Newton diverges there rather than returning a number.
    hot = solve_robust(lad3, p0_critical(lad3, hi=HI[3]))
    dt_hot = float(np.max(hot[:, 3] - hot[:, 1]))
    start = solve_robust(Ladder(RX, rung=3, n_z=N_Z, T_r=bk("startup_T0_1182"),
                                film_scale=fs), P0_REF)
    dt_start = float(np.max(start[:, 3] - start[:, 1]))
    pdt = Pellet(RX, n_r=n_p, isothermal=False).solve(THETA * P_HS, T_HS)["T_centre"] - T_HS
    i = int(np.argmax(start[:, 1]))
    r_obs = kexp(start[i, 3]) * P["p_B0"] * start[i, 2] / 3600.0
    chi = ((P["minus_dH"] * KCAL) * r_obs * P["rho_b"] * D_P
           / (2.0 * film(RX, fs)["h_f"] * start[i, 1]))
    return dict(film_dT_at_hotspot_K=dt_hot, film_dT_at_startup_K=dt_start,
                film_dT_over_pellet_dT=dt_hot / pdt,
                mears_margin=(bk("mears_constant") / (P["a"] / start[i, 1])) / chi)


def h_control(p0s=None):
    """The crossing control on the five ladder rungs, as the Validation cell ran it:
    each rung's own threshold, solved either side from the branch state its own
    bisection stored.  `p0s` hands each rung a different threshold, which is how a
    break row tests that the control rejects one that is not the rung's own."""
    if p0s is None:
        for r, *_ in RUNGS:
            assert abs(LAD[r].branch[0] / P0S[r] - 1.0) < 1e-8, (
                f"rung {r}: the stored branch state is not at this rung's threshold")
    p0s = P0S if p0s is None else p0s
    out = [crossing_control(LAD[r], p0s[r], guess=LAD[r].branch[1]) for r, *_ in RUNGS]
    return dict(ladder_controls_rejected=sum(1 for *_, ok in out if not ok))


PE_FINE = tuple(float(x) for x in np.round(np.arange(0.316, 0.3075, -0.001), 4))


def h_rejections(jump=None):
    """The two rejection counts, recomputed over the ladder's five controls and the
    9-point fine Pe_a grid of the crossover cell, with the control's jump allowance
    as an argument.  The page-wide counts cover more thresholds than this; a row
    that injects a defect here recomputes what it can afford, and the counts it
    returns are rejections that defect produces.  CTL_STATS is restored, so a break
    row does not leak into the page-wide count."""
    saved = dict(CTL_STATS)
    try:
        lad_rej = sum(1 for r, *_ in RUNGS
                      if not crossing_control(LAD[r], P0S[r], guess=LAD[r].branch[1],
                                              jump=jump)[2])
        grid_rej = sum(1 for p in PE_FINE if not _shift_or_none(
            lambda x: 100.0 * (p0_critical(Ladder(RX, rung=2, n_z=600, pe_ma=x, pe_ha=x),
                                           hi=0.09, jump=jump) / BASE600 - 1.0), p)[1])
    finally:
        CTL_STATS.clear()
        CTL_STATS.update(saved)
    return dict(scan_points_rejected=grid_rej, thresholds_rejected=lad_rej + grid_rej)


def h_chapter_internal(E_over_R=None, E1_1173=None):
    """The two chapter-against-itself numbers.  Both hang on the same printed
    constant - section 11.5.2's own E/R - which is why one row moves both, and
    they were the other two agreement.json entries with no helper at all."""
    a_ = P["a"] if E_over_R is None else E_over_R
    e1 = bk("E1_over_R_1173") if E1_1173 is None else E1_1173
    k = float(np.exp(np.clip(P["b"] - a_ / T_HS, -700, 50)))
    return dict(E_over_R_internal_dev=a_ / e1 - 1.0,
                lump_deficit=k / (K1_HS + K3_HS) - 1.0)


BASE_M = {}
BASE_M.update(h_ex1151a()); BASE_M.update(h_units()); BASE_M.update(h_groups())
BASE_M.update(h_equiv1d()); BASE_M.update(h_tables()); BASE_M.update(h_eta())
BASE_M.update(h_rung1()); BASE_M.update(h_cross()); BASE_M.update(h_cross_R())
BASE_M.update(h_units_model())
BASE_M.update(h_film_dt()); BASE_M.update(h_chapter_internal()); BASE_M.update(h_control())
BASE_M.update(h_rejections())
BASE_SHIFTS = h_shifts()
BASE_M.update(BASE_SHIFTS)

# ---- the four shifts at BOTH resolutions, so the coarse ones are reported too ----
SH = pd.DataFrame([dict(rung=r,
                        displayed_n_z_1200_n_r_20=SHIFT[r],
                        break_rows_n_z_600_n_r_10=BASE_SHIFTS[f"shift_rung{r}_pct"],
                        gap_pp=BASE_SHIFTS[f"shift_rung{r}_pct"] - SHIFT[r])
                   for r in (2, 3, 4, 5)])
show(SH.style.set_uuid("d11shiftgrid").hide(axis="index").format(
    {"displayed_n_z_1200_n_r_20": "{:+.4f}", "break_rows_n_z_600_n_r_10": "{:+.4f}",
     "gap_pp": "{:+.4f}"}))
SHIFT_GRID_GAP = float(SH.gap_pp.abs().max())
print(f"agreement.json carries the RIGHT-HAND column, because every break row that"
      f" touches a\nshift has to re-run it. The gap to the displayed values is at most"
      f" {SHIFT_GRID_GAP:.3f} percentage\npoints, and it is printed here rather than left to"
      f" be discovered. Rung 5's row is\nthe reason rung 5's threshold is called unrefined"
      f" above.\n")

# the helpers ARE the route to the reported numbers, asserted rather than asserted-to
_R1M = R1["metrics"]
# ONE list, used BOTH to check each helper against its printed number and to check
# that no metric escapes the check.  Written once, so "checked as a set" cannot
# drift between two copies.
_ASSERTED = (("ex1151a_worst_dev", EX1151A_WORST),
                ("dp_definition_cost", DP_DEFINITION_COST),
                ("re_implied_spread", RE_IMPLIED_SPREAD),
                ("conversion_worst_dev", CONV_WORST),
                ("backint_residual", BACKINT_RESIDUAL),
                ("group_B_dev", B_PRINTED_DEV),
                ("group_A_printed_dev", A_PRINTED_DEV),
                ("cp_printed_dev", CP_PRINTED_DEV),
                ("cp_route_spread", CP_SPREAD),
                ("U_from_2d_dev", U_FROM_2D_DEV),
                ("U_from_2d_dev_at_printed_radius", U_PRINTED_R_DEV),
                ("table2_upper_ratio", T2_UPPER_RATIO),
                ("table1_worst_dev", T1_WORST),
                ("table2_dTad_worst_dev", T2_DTAD_WORST_DEV),
                ("table2_Q_implied_by_p0", T2_Q_IMPLIED),
                ("eta_grid_vs_11_9_1_11", ETA_DEV),
                ("eta_vs_B11_worst", B11_DEV),
                ("crossover_pe_axial", PE_STAR),
                ("crossover_film_scale", FILM_STAR),
                ("crossover_D_e", DE_STAR),
                ("crossover_R_t", R_STAR),
                ("p0_star_bar_reading_dev", P0_STAR_BAR_DEV),
                ("T_hot_0018_bar_reading_C", T_0018_BAR_C),
                ("rung1_p0_star_extrapolated", _R1M["rung1_p0_star_extrapolated"]),
                ("rung1_order_p0", _R1M["rung1_order_p0"]),
                ("rung1_order_T", _R1M["rung1_order_T"]),
                ("ode_vs_pymrm_p0", _R1M["ode_vs_pymrm_p0"]),
                ("ode_vs_pymrm_T", _R1M["ode_vs_pymrm_T"]),
                ("p0_star_vs_printed_0018", _R1M["p0_star_vs_printed_0018"]),
                ("rung1_finest_err_over_start_floor", R1_FINEST_OVER_FLOOR),
                ("shift_rung2_pct", SH.break_rows_n_z_600_n_r_10.iloc[0]),
                ("shift_rung3_pct", SH.break_rows_n_z_600_n_r_10.iloc[1]),
                ("shift_rung4_pct", SH.break_rows_n_z_600_n_r_10.iloc[2]),
                ("shift_rung5_pct", SH.break_rows_n_z_600_n_r_10.iloc[3]),
                ("film_dT_at_hotspot_K", FILM_DT_HS),
                ("film_dT_at_startup_K", DT_FILM_1182),
                ("film_dT_over_pellet_dT", FILM_DT_HS / PELLET_DT),
                ("mears_margin", MEARS_RHS / CHI),
                ("E_over_R_internal_dev", ER_INTERNAL),
                ("lump_deficit", LUMP_DEFICIT),
                ("ladder_controls_rejected", int((~CTL.crossing).sum())),
                ("scan_points_rejected", SCAN_REFUSED + (len(_UNIFORM) - UNIFORM_OK)),
                ("thresholds_rejected", CTL_STATS["rejected"]))
_MISMATCH = [(k, BASE_M[k], float(live)) for k, live in _ASSERTED
             if abs(BASE_M[k] - live) > 1e-12 * max(1.0, abs(live))]
assert not _MISMATCH, ("a helper does not reproduce the number reported above, so a "
                       "break row perturbing it would not reach that number: "
                       + repr(_MISMATCH))
_UNASSERTED = sorted(set(BASE_M) - {k for k, _ in _ASSERTED})
assert not _UNASSERTED, f"helper metrics not asserted against a printed number: {_UNASSERTED}"
print(f"ALL {len(BASE_M)} helper metrics are asserted equal, to 1e-12, to a number"
      f" printed earlier in\nthis notebook - none excepted. For the seven rung-1 metrics the"
      f" assertion is an\nidentity: the Validation cell displayed this helper's own return"
      f" value.")


# %%
# ---------------------------- the break rows ---------------------------------
def misread(key, swap=None, digit=None):
    """A printed number misread: the characters at `swap` and `swap + 1` of its
    printed string transposed, or one digit read as another."""
    t = bk(key, "as_printed")
    if swap is not None:
        t = t[:swap] + t[swap + 1] + t[swap] + t[swap + 2:]
    else:
        assert t.count(digit[0]) == 1, t
        t = t.replace(*digit)
    assert float(t) != bk(key), (key, t)
    return float(t)


def brk_dp_definition():
    """Example 11.5.1.A's d_p read as eq. (11.5.1-9) defines it, 6 V/S = d."""
    return h_ex1151a(dp=DP_SV)


def brk_ergun_refit():
    """A1.1's refit of Ergun's own 244 markers, k1 = 151.9, k2 = 1.697."""
    return h_ex1151a(a=A11_C["k2_refit_all_points"], b=A11_C["k1_refit_all_points"])


def brk_atm_factor():
    """the atm-to-bar factor set to 1, i.e. 11.5.2 read as a straight copy."""
    return dict(**h_units(atm=1.0), **h_tables(atm=1.0))


def brk_table_A():
    """Table 11.5.3.A-2's cross-column test on the A its own formula gives, 6165,
    instead of the printed 6150."""
    return h_tables(A=A_G)


def brk_unit_reading():
    """section 11.5.2's pressures read in VWF's atm rather than the printed bar."""
    return h_units_model(p_B0=P["p_B0"], p_t=P["p_t"])


def brk_target_R_2pct():
    """the tube-radius crossover defined as a 2 % shift instead of 1 %."""
    return h_cross_R(target=2.0)


def brk_pB0_bar():
    """p_B0 read as the 0.211 bar of 11.5.2 rather than the 0.208 atm of 11.5.3.A."""
    return h_units(p_B0_book=V["p_B0"])


def brk_cp_printed():
    """the volumetric c_p replaced by the printed 0.992 kJ/(kg K) x rho_g."""
    return dict(h_groups(cp_vol=bk("specific_heat") * P["rho_g"] / KCAL),
                **{})


def brk_dH_kcal_kJ():
    """(-dH) taken as the kJ/kmol of 11.5.2 without converting back to kcal."""
    return h_groups(minus_dH=bk("heat_of_reaction"))


def brk_radius_printed():
    """R_t = 0.0127 m, the radius the printed 2.54 cm diameter implies."""
    return dict(**h_equiv1d(R=D_T_PRINTED / 2.0), **h_shifts(R=D_T_PRINTED / 2.0),
                **h_rung1(R=D_T_PRINTED / 2.0))


def brk_lam_er_5pct():
    """lambda_er raised 5 %, the parameter 11.7.3 calls the sensitive one."""
    return dict(**h_equiv1d(lam_er=LAM_ER * 1.05),
                **h_shifts(lam_er=LAM_ER * 1.05, pe_hr=PE_HR / 1.05))


def brk_alpha_w_5pct():
    """alpha_w raised 5 %, the other sensitive one."""
    return h_equiv1d(alpha_w=ALPHA_W * 1.05)


def brk_table_Q():
    """Table 11.5.3.A-2's Q read as 3.4657 - two digits transposed."""
    return h_tables(Q=misread("vwf_Q_625", swap=4))


def brk_table_dT():
    """the permissible rise read as 31.9 K instead of 31.6 K."""
    return h_tables(dT=misread("vwf_dT_625", digit=("6", "9")))


def brk_pellet_nu():
    """construct_div nu = 1 (cylindrical) in the pellet, which is a sphere."""
    return h_eta(nu=1)


def brk_E_over_R_swap():
    """The chapter's two printed activation energies swapped: section 11.7.3's
    E_1/R = 13,588 K used where section 11.5.2's 13,636 K is, and the other way
    round.  Both chapter-against-itself metrics hang on that one constant."""
    return h_chapter_internal(E_over_R=bk("E1_over_R_1173"), E1_1173=P["a"])


def brk_pellet_coarse():
    """the pellet on 16 radial cells instead of 256 (and 16 instead of 128 for the
    non-isothermal solve the film/pellet ratio uses)."""
    return dict(**h_eta(n_r=16), **h_film_dt(n_r_pellet=16))


def brk_eta_csv():
    """B1.1's exact eta column perturbed in its last digit."""
    return h_eta(eta_col=SPH.eta.values * (1.0 + 1e-4))


def brk_pe_axial():
    """Pe_ma = Pe_ha at half the chapter's own lower limit, i.e. axial dispersion
    twice as strong; inside the validated region of the sweep the crossover row
    scans."""
    return h_shifts(pe_ma=bk("pe_axial_low") / 2.0, pe_ha=bk("pe_axial_low") / 2.0)


def brk_film_scale():
    """the Wakao-Funazkri coefficients halved."""
    return dict(**h_shifts(film_scale=0.5), **h_film_dt(film_scale=0.5))


def brk_De():
    """D_e raised a decade from the value rung 4 is evaluated at."""
    return h_shifts(D_e=10.0 * DE_SWEEP[0])


def brk_pe_mr_literal():
    """page 574's a_1 group read literally, i.e. Pe_mr divided by eps."""
    return h_shifts(pe_mr=PE_MR / EPS)


def brk_ceiling():
    """the permissible ceiling read as 410 C instead of the printed 415 C."""
    return dict(**h_shifts(ceiling=T_PERM - 5.0), **h_rung1(ceiling=T_PERM - 5.0))


def brk_pB0_group():
    """p_B0 read as the 0.211 bar section 11.5.2 prints, in the kcal groups."""
    return h_groups(p_B0=bk("oxygen_pressure"))


def brk_table1_Q():
    """Table 11.5.3.A-1's Q read as 2.9023 - two digits transposed."""
    return h_tables(Q1=misread("vwf_Q_635", swap=3))


def brk_eta_off():
    """rung 4 with eq. (11.9.1-11) switched off, i.e. eta forced to 1."""
    return h_shifts(eta_on=False)


def brk_target_up():
    """The crossover definition read as a 2 % shift instead of 1 %.  (The Pe_a
    sweep validates every point and reaches 2.5 %, printed in the crossover cell.)"""
    return h_cross(target=2.0)


def brk_target_half():
    """the crossover definition read as a 0.5 % shift instead of 1 %."""
    return h_cross(target=0.5)


def brk_pe_extrap_baseline():
    """The rung-2 shift measured against rung 1's Richardson-extrapolated
    threshold instead of the same-grid one: the mixed-baseline defect the
    crossover cell measures.  Rung 1's own 600-cell discretisation error is a
    fifth of the target, so the crossover moves."""
    sc = cross_scan(lambda p: 100.0 * (p0_critical(Ladder(RX, rung=2, n_z=600, pe_ma=p,
                                                          pe_ha=p), hi=0.09)
                                       / LIM_P - 1.0), PE_SWEEP, TARGET)
    return dict(crossover_pe_axial=sc["value"])


def brk_richardson_skipped():
    """rung 1's finest grid compared with the ODE route, no Richardson limit."""
    return h_rung1(extrapolate=False)


def brk_richardson_ratio():
    """the rung-1 orders fitted with the grid ratio read as 3 instead of 2."""
    return h_rung1(ratio_read=3.0)


def brk_wrong_rung_threshold():
    """each rung's crossing control run at the NEIGHBOURING rung's threshold (rung 1
    at rung 2's, every other rung at rung 1's): a threshold that is not the rung's
    own, which the control must reject."""
    return h_control(p0s={r: P0S[2 if r == 1 else 1] for r, *_ in RUNGS})


def brk_jump_allowance():
    """the crossing control's jump allowance read as 0.005 K instead of 1 K - below
    every jump the five rungs measure (the control table above), so a correct
    threshold is rejected and the rejection counts must leave zero."""
    return h_rejections(jump=0.005)


def brk_cross_coarse():
    """The crossovers located on a coarser grid, 300 cells instead of 600.  (The
    crossover cell prints rung 4's anchor point at 300 cells as a clean crossing.)"""
    return h_cross(n_z=300)


BREAK_FNS = [
    ("crossover definition read as 2 %, not 1 %", brk_target_up),
    ("crossover definition read as 0.5 %, not 1 %", brk_target_half),
    ("crossovers located on a 300-cell grid", brk_cross_coarse),
    ("rung-2 shift against rung 1's extrapolated threshold", brk_pe_extrap_baseline),
    ("rung 1: finest grid against the ODE, no Richardson limit", brk_richardson_skipped),
    ("rung 1: orders fitted with the grid ratio read as 3", brk_richardson_ratio),
    ("ladder controls run at the neighbouring rung's threshold", brk_wrong_rung_threshold),
    ("crossing control's jump allowance read as 0.005 K", brk_jump_allowance),
    ("tube-radius crossover read as 2 %, not 1 %", brk_target_R_2pct),
    (f"Table 11.5.3.A-2 cross-checked against A = {A_G:.0f}", brk_table_A),
    ("11.5.2's pressures read in atm, not the printed bar", brk_unit_reading),
    ("Ex 11.5.1.A: d_p from eq. (11.5.1-9), 6V/S = d", brk_dp_definition),
    (f"Ergun constants: A1.1's refit {A11_C['k1_refit_all_points']:.1f} /"
     f" {A11_C['k2_refit_all_points']:.3f}", brk_ergun_refit),
    ("units: atm-to-bar factor set to 1", brk_atm_factor),
    (f"units: p_B0 read as {V['p_B0']} bar, not {V['p_B0']} atm converted", brk_pB0_bar),
    (f"c_p: the printed {bk('specific_heat', 'as_printed')} kJ/(kg K) used instead", brk_cp_printed),
    ("(-dH): the kJ/kmol value used in the kcal groups", brk_dH_kcal_kJ),
    (f"R_t = {D_T_PRINTED / 2.0} m, from the printed {bk('tube_id', 'as_printed')} diameter",
     brk_radius_printed),
    ("lambda_er raised 5 %", brk_lam_er_5pct),
    ("alpha_w raised 5 %", brk_alpha_w_5pct),
    (f"Table 11.5.3.A-2: Q read {misread('vwf_Q_625', swap=4)}", brk_table_Q),
    (f"Table 11.5.3.A-2: dT read {misread('vwf_dT_625', digit=('6', '9'))} K", brk_table_dT),
    ("pellet: construct_div nu = 1, not 2", brk_pellet_nu),
    ("pellet: 16 radial cells, not 256", brk_pellet_coarse),
    ("B1.1's exact eta column shifted 1e-4", brk_eta_csv),
    (f"Pe_ma = Pe_ha = {bk('pe_axial_low') / 2.0} (2x the chapter's lower limit)", brk_pe_axial),
    ("film coefficients halved", brk_film_scale),
    ("D_e raised to 1e-5 m2/s", brk_De),
    ("Pe_mr from page 574's a_1 read literally", brk_pe_mr_literal),
    (f"permissible ceiling read as {T_PERM - 5.0 - 273.15:.0f} C", brk_ceiling),
    (f"p_B0 read as {bk('oxygen_pressure', 'as_printed')} bar in the kcal groups", brk_pB0_group),
    (f"Table 11.5.3.A-1: Q read {misread('vwf_Q_635', swap=3)}", brk_table1_Q),
    ("11.5.2's E/R and 11.7.3's E_1/R swapped", brk_E_over_R_swap),
    ("rung 4 with eq. (11.9.1-11) switched off", brk_eta_off),
]

METRICS = sorted(BASE_M)
MOVE_TOL = 1e-6
rows, COVERAGE, BREAK_GOT = [], {}, {}
for lbl, fn in BREAK_FNS:
    got = fn()
    BREAK_GOT[lbl] = got
    unknown = sorted(set(got) - set(METRICS))
    assert not unknown, f"break row '{lbl}' recomputes names that are not metrics: {unknown}"
    moved = []
    for k, v in got.items():
        base = BASE_M[k]
        rel = abs(v - base) / max(abs(base), 1e-12)
        if rel > MOVE_TOL:
            moved.append(k)
            COVERAGE.setdefault(k, []).append((lbl, rel))
    rows.append(dict(injected=lbl, metrics_touched=len(got), metrics_moved=len(moved),
                     worst_rel_move=max((abs(got[k] - BASE_M[k]) / max(abs(BASE_M[k]), 1e-12)
                                         for k in got), default=0.0)))
BREAKS = pd.DataFrame(rows)
show(BREAKS.style.set_uuid("d11breaks").hide(axis="index").format(
    {"worst_rel_move": "{:.3e}"}))

# ---- a break row's OWN cross_scan refinement is not announced above, unlike the ----
# ---- four reported rows' (printed earlier, in the crossover cell) - so print it ----
_DIAG = [(t, n, k, d) for (t, n), dd in sorted(_H_CROSS_DIAG.items()) for k, d in dd.items()]
print("\nTHE BREAK-ROW SCANS' OWN REFINEMENT, which the table above does not show:")
for _t, _n, _k, _d in _DIAG:
    print(f"   target {_t:4.1f} %  n_z = {_n:4d}  {_k:<3s}: {_d['levels']} levels,"
          f" {_d['probed']:3d} probed, {_d['refused'] + _d['n_evaluated'] - _d['n_controlled']}"
          f" rejected, bracket {_d['bracket']:.3g}"
          f" ({100.0 * _d['bracket'] / abs(_d['value']):.2g} % of the value)")


# %%
# ---- and the rows must COMPUTE what they return, parsed rather than trusted ----
# A row that returns a TYPED CONSTANT for a metric records a move of exactly 1
# whatever the reported value is: a hand-written coverage claim wearing the
# generator's clothes.  This guard is copied from published J4.2.  It parses each row's own source plus one level of the
# helpers that row names, and rejects any metric key bound to a numeric literal or
# to a local that is only ever assigned one.  What it cannot say is whether the
# expression is the RIGHT one; that is what the measured moves are for.
import ast, inspect, textwrap

_NUM_CTORS = {"float", "int", "complex", "float64", "float32"}


def _const(node, env, depth=0):
    """True if `node` is PROVABLY a numeric constant.  Constant-folds, so it also
    catches the laundering forms `0.1656 + 0.0001`, `float("15.688")`, `_T[0]` off
    a frozen tuple, and a helper whose every `return` is a literal; the negative
    controls below are exactly those four plus the bare literal."""
    if depth > 6:
        return False
    if isinstance(node, ast.Constant):
        return isinstance(node.value, (int, float)) and not isinstance(node.value, bool)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        return _const(node.operand, env, depth + 1)
    if isinstance(node, ast.BinOp):
        return _const(node.left, env, depth + 1) and _const(node.right, env, depth + 1)
    if isinstance(node, ast.Name):
        return node.id in env and _const(env[node.id], env, depth + 1)
    if isinstance(node, ast.Subscript):
        base = node.value
        if isinstance(base, ast.Name) and base.id in env:
            base = env[base.id]
        return (isinstance(base, (ast.Tuple, ast.List))
                and all(_const(e, env, depth + 1) for e in base.elts))
    if isinstance(node, ast.Call):
        fname = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
        if fname in _NUM_CTORS and len(node.args) == 1:
            a = node.args[0]
            if isinstance(a, ast.Constant) and isinstance(a.value, str):
                try:
                    float(a.value)
                    return True
                except ValueError:
                    return False
            return _const(a, env, depth + 1)
        obj = globals().get(getattr(node.func, "id", ""), None)
        if (inspect.isfunction(obj) and getattr(obj, "__module__", "") == "__main__"):
            try:
                sub = ast.parse(textwrap.dedent(inspect.getsource(obj)))
            except (OSError, TypeError):
                return False
            senv = _frozen(sub)
            rets = [n.value for n in ast.walk(sub)
                    if isinstance(n, ast.Return) and n.value is not None]
            return bool(rets) and all(_const(r, senv, depth + 1) for r in rets)
    return False


def _frozen(tree):
    """names assigned once, to something constant-foldable or to a frozen tuple.

    TUPLE UNPACKING IS INCLUDED: `a, b = 0.1657, 2.0` binds a local that only ever
    carries a literal, which is a form the guard's own positive claim names."""
    env = {}
    for n in ast.walk(tree):
        if isinstance(n, ast.Assign):
            for t in n.targets:
                if isinstance(t, ast.Name):
                    env[t.id] = n.value
                elif (isinstance(t, ast.Tuple) and isinstance(n.value, ast.Tuple)
                        and len(t.elts) == len(n.value.elts)):
                    for tt, vv in zip(t.elts, n.value.elts):
                        if isinstance(tt, ast.Name):
                            env[tt.id] = vv
    return env


def _key_bindings(src):
    tree, out = ast.parse(textwrap.dedent(src)), []
    for n in ast.walk(tree):
        if isinstance(n, ast.Dict):
            out += [(k.value, v) for k, v in zip(n.keys, n.values)
                    if isinstance(k, ast.Constant) and isinstance(k.value, str)]
        elif isinstance(n, ast.Call) and getattr(n.func, "id", "") == "dict":
            out += [(kw.arg, kw.value) for kw in n.keywords if kw.arg]
        elif isinstance(n, ast.Assign):
            for t in n.targets:
                if (isinstance(t, ast.Subscript) and isinstance(t.slice, ast.Constant)
                        and isinstance(t.slice.value, str)):
                    out.append((t.slice.value, n.value))
    return out


def _sources(fn):
    src = textwrap.dedent(inspect.getsource(fn))
    seen = [src]
    for n in ast.walk(ast.parse(src)):
        if isinstance(n, ast.Name):
            obj = globals().get(n.id)
            if (inspect.isfunction(obj) and getattr(obj, "__module__", "") == "__main__"
                    and obj is not fn):
                try:
                    seen.append(textwrap.dedent(inspect.getsource(obj)))
                except (OSError, TypeError):
                    seen.append(f"# UNREADABLE SOURCE: {n.id}")
    return seen


def literal_metrics(fn):
    bad = []
    for src in _sources(fn):
        tree = ast.parse(textwrap.dedent(src))
        env = _frozen(tree)
        bad += [k for k, v in _key_bindings(src) if k in METRICS and _const(v, env)]
    return sorted(set(bad))


def _nc_bare():
    """1. the bare literal."""
    return {"shift_rung2_pct": 0.0, "shift_rung3_pct": -4.7,
            "shift_rung4_pct": 15.7, "eta_grid_vs_11_9_1_11": 1e-5,
            "_diagnostic": np.log(2.0)}


def _nc_arith():
    """2. two literals added, and one multiplied by 1."""
    return {"shift_rung2_pct": 0.1656 + 0.0001, "shift_rung3_pct": -4.744 * 1.0}


def _nc_cast():
    """3. a literal laundered through float()."""
    return {"shift_rung4_pct": float("15.688")}


def _nc_tuple():
    """4. a frozen tuple, subscripted."""
    _T = (0.2488, 2.2665e-05)
    return {"crossover_pe_axial": _T[0], "crossover_D_e": _T[1]}


def _nc_helper_value():
    """the one-line laundering helper of control 5 - all returns are literals."""
    return 0.1657


def _nc_helper():
    """5. laundered through a __main__ helper - the case the guard ADVERTISES."""
    return {"shift_rung2_pct": _nc_helper_value()}


def _nc_unpack():
    """6. a local bound by TUPLE UNPACKING, which the guard's own positive claim -
    "a local that only ever carries one" - covers exactly.  `_frozen` reads tuple
    targets, so it is a negative control rather than an escape."""
    _u, _v = 0.2488, 2.0
    return {"crossover_pe_axial": _u}


_NC_GLOBAL = 0.1657          # a module-level constant, for the escape list below


def _nc_global():
    """a metric bound to a name assigned OUTSIDE the row's own source"""
    return {"shift_rung2_pct": _NC_GLOBAL}


def _nc_cond():
    return {"shift_rung2_pct": 0.1657 if True else 0.0}


def _nc_round():
    return {"shift_rung2_pct": round(0.1657, 4)}


def _nc_array():
    return {"shift_rung2_pct": np.array([0.1657])[0]}


def _nc_fstring():
    return {"shift_rung2_pct": float(f"{0.1657}")}


def _nc_sum():
    return {"shift_rung2_pct": sum([0.1657])}


def _nc_comprehension():
    return {k: 0.1657 for k in ("shift_rung2_pct",)}


def _nc_lambda():
    return {"shift_rung2_pct": (lambda: 0.1657)()}


def _nc_update():
    d = {}
    d.update(shift_rung2_pct=0.1657)
    return d


def _nc_zip():
    return dict(zip(["shift_rung2_pct"], [0.1657]))


def _nc_setdefault():
    d = {}
    d.setdefault("shift_rung2_pct", 0.1657)
    return d


def _nc_abs():
    """abs() of a literal: constant-FOLDABLE, but a call and not arithmetic, and
    `abs` is not one of the numeric constructors `_const` unwraps.  It is an
    escape, which is why the positive claim below does not say
    'constant-foldable'."""
    return {"shift_rung2_pct": abs(0.1657)}


def _nc_annassign():
    """13. an ANNOTATED assignment: `_a: float = 0.1657` is an ast.AnnAssign, which
    `_frozen` never reads at all (it only walks ast.Assign).  The guard's own
    positive claim about "a local name ... that only ever carries one" is a
    category, and this is one more form inside it that escapes."""
    _a: float = 0.1657
    return {"shift_rung2_pct": _a}


def _nc_list_target():
    """14. a LIST-target assignment: `[_a, _b] = 0.1657, 2.0` binds `_a` to a
    literal exactly the way tuple unpacking does, but `_frozen` only matches an
    ast.Tuple target, never an ast.List one, so it is missed."""
    [_a, _b] = 0.1657, 2.0
    return {"shift_rung2_pct": _a}


def _nc_nested_unpack():
    """15. NESTED tuple unpacking: `(_a, _b), _c = (0.1657, 2.0), 3.0` puts a
    Tuple where `_frozen` expects a bare Name inside the outer target's elements,
    so the recursion into the inner tuple never happens and `_a` is missed."""
    (_a, _b), _c = (0.1657, 2.0), 3.0
    return {"shift_rung2_pct": _a}


def _nc_walrus():
    """16. a WALRUS binding: `(_a := 0.1657)` is an ast.NamedExpr, a third kind of
    name-binding node `_frozen` never reads at all (it only walks ast.Assign)."""
    return {"shift_rung2_pct": (_a := 0.1657)}


def _fp_computed():
    """FALSE-POSITIVE control, and this one is not vacuous. It binds a METRIC key
    to a computed value; the guard must not flag it.  A NON-metric key would not
    do: `literal_metrics` filters it out before it looks at anything."""
    return {"shift_rung2_pct": float(np.log(2.0))}


ESCAPES = [("a module-level constant", _nc_global),
           ("0.1657 if True else ...", _nc_cond),
           ("round(0.1657, 4)", _nc_round),
           ("np.array([0.1657])[0]", _nc_array),
           ('float(f"{0.1657}")', _nc_fstring),
           ("sum([0.1657])", _nc_sum),
           ("a dict comprehension", _nc_comprehension),
           ("a lambda", _nc_lambda),
           ("d.update(key=0.1657)", _nc_update),
           ("dict(zip([key], [0.1657]))", _nc_zip),
           ("d.setdefault(key, 0.1657)", _nc_setdefault),
           ("abs(0.1657)", _nc_abs),
           ("an annotated assignment, `_a: float = 0.1657`", _nc_annassign),
           ("a list-target assignment, `[_a, _b] = 0.1657, 2.0`", _nc_list_target),
           ("nested tuple unpacking, `(_a, _b), _c = (0.1657, 2.0), 3.0`",
            _nc_nested_unpack),
           ("a walrus binding, `(_a := 0.1657)`", _nc_walrus)]

NEG_CONTROLS = [("bare literal", _nc_bare,
                 {"shift_rung2_pct", "shift_rung3_pct", "shift_rung4_pct",
                  "eta_grid_vs_11_9_1_11"}),
                ("literal arithmetic", _nc_arith,
                 {"shift_rung2_pct", "shift_rung3_pct"}),
                ("float(\"...\") cast", _nc_cast, {"shift_rung4_pct"}),
                ("frozen tuple, subscripted", _nc_tuple,
                 {"crossover_pe_axial", "crossover_D_e"}),
                ("laundered through a helper", _nc_helper, {"shift_rung2_pct"}),
                ("tuple unpacking", _nc_unpack, {"crossover_pe_axial"})]
LITERAL_ROWS = {lbl: bad for lbl, fn in BREAK_FNS for bad in [literal_metrics(fn)] if bad}
print(f"{len(BREAK_FNS)} rows parsed; {len(LITERAL_ROWS)} return a metric as a typed"
      f" constant.\nTHE GUARD'S TEETH, MEASURED on {len(NEG_CONTROLS)} laundering forms:")
_TEETH = []
for lbl, fn, want in NEG_CONTROLS:
    got = set(literal_metrics(fn))
    _TEETH.append(got == want)
    print(f"   {lbl:<28s} caught {sorted(got)}   -> {'PASS' if got == want else 'FAIL'}")
assert all(_TEETH), "the literal guard does not catch every negative control"
assert not literal_metrics(_fp_computed), "the guard flagged a COMPUTED metric value"
assert not LITERAL_ROWS, f"break rows returning typed constants: {LITERAL_ROWS}"
print(f"   false positives: a metric key bound to float(np.log(2.0)) is NOT flagged"
      f" ({literal_metrics(_fp_computed) == []}).")
print(f"\nAND WHAT IT DOES NOT CATCH: AT LEAST {len(ESCAPES)} more forms, measured rather"
      f" than left to be\ndiscovered - each of these binds a metric key to 0.1657 and walks"
      f" through. AT LEAST is\nthe honest quantifier: an independent sabotage set can"
      f" find a form this list\ndoes not name. The list:")
_ESCAPED = [(lbl, literal_metrics(fn)) for lbl, fn in ESCAPES]
print("   " + "; ".join(lbl for lbl, got in _ESCAPED if not got))
assert all(not got for _, got in _ESCAPED), (
    "an escape form is now caught - update the printed list: "
    + repr([lbl for lbl, got in _ESCAPED if got]))
print(f"So the accurate claim is a LIST OF FORMS and not a category: a metric key bound to"
      f"\na numeric literal, to +-*/** arithmetic or unary minus on literals, to a"
      f" float/int/\nnp.float64 cast of one, to a subscript of a frozen literal tuple, or to"
      f" a local name -\nin the row's own body or a one-level __main__ helper - bound by a"
      f" plain Name or\nflat-tuple-unpacking assignment (annotated, list-target, nested-tuple"
      f" and walrus\nbindings escape). NOT 'constant-foldable': abs(0.1657) is foldable,"
      f" but a call to\nsomething that is not a numeric constructor, and it escapes. And NOT 'any row"
      f" returning a metric as\na typed constant'. The realistic accidental case, a"
      f" module-level constant, is in\nthe escape list. What the guard still cannot say is"
      f" whether the expression is the\nRIGHT one; that is what the measured moves in the"
      f" table above are for.")

# ---------------- coverage map, GENERATED from the measured moves --------------
CAVEATS = {
    "cp_route_spread": ("A SPREAD among four routes to c_p, not an agreement among "
                        "four witnesses - the Parameters cell shows route 2 reproduces "
                        "route 1 by construction and route 4 is route 1 rounded, so "
                        "the four trace to one printed determination plus route 3. "
                        "Perturbing any one of them moves this number, and "
                        "brk_cp_printed and brk_dH_kcal_kJ both do; the caveat is "
                        "that it cannot detect a c_p that is wrong in the same "
                        "way in all four places, which is precisely what a shared "
                        "source looks like. C2.10 warns that U = 82.7 is such a "
                        "shared source between it and D2.2. The CONCLUSION it "
                        "supports - that the printed 0.992 is inconsistent with the "
                        "chapter's own groups - rests on cp_printed_dev, not on this."),
}
cov = []
for k in METRICS:
    movers = COVERAGE.get(k, [])
    cov.append(dict(metric=k, value=BASE_M[k], n_movers=len(movers),
                    best_mover=(max(movers, key=lambda t: t[1])[0] if movers else ""),
                    best_rel_move=(max(t[1] for t in movers) if movers else np.nan)))
COV = pd.DataFrame(cov)
show(COV.style.set_uuid("d11cov").hide(axis="index").format(
    {"value": "{:.6g}", "best_rel_move": "{:.2e}"}))
uncovered = sorted(set(METRICS) - set(COVERAGE))
assert not uncovered, f"metrics no row moves: {uncovered}"
print(f"{len(METRICS)} metrics, {len(COVERAGE)} moved by at least one row; none is"
      f" structural.")
for k, why in CAVEATS.items():
    print(f"\nCAVEAT - {k}:\n   {why}")
print(f"\nTHE CROSSING CONTROL, COUNTED OVER THE WHOLE PAGE. Every threshold this"
      f" notebook\nhas computed up to this point - the ladder table, the six-grid rung-1"
      f" chain, the\nrung-2 grid chain, the four crossover scans and all"
      f" {len(BREAK_FNS)} break rows -"
      f" has been\nchecked: {CTL_STATS['checked']} bisected thresholds, of which"
      f" {CTL_STATS['rejected']} were REJECTED as not being ceiling\ncrossings. The"
      f" smallest distance of any control solve from the ceiling, or of any\njump from"
      f" the {CTL_JUMP:.0f} K allowance, is {CTL_STATS['margin_K']:.2e} K, so no verdict"
      f" sits near round-off.")


# %%
# ============ claims: every published number says what it is ===================
# docs/claims-design.md.  One registration per metric.
# VALUES are physical, data and discretisation results; CI holds them to rtol of the
# committed baseline.  The default rtol is 0.05; it is tightened ONLY where a break
# row needs it, and only for numbers that are exact arithmetic on printed constants
# or bisected to 1e-12, so the tighter band is still far above any machine spread.
claims = Claims("D1.1")
RTOL = {
    "dp_definition_cost": 1e-3,          # printed-constant arithmetic; Ergun refit moves 0.8 %
    "table2_upper_ratio": 1e-3,          # printed-constant arithmetic; dT row moves 0.9 %
    "table2_dTad_worst_dev": 1e-3,       # printed-constant arithmetic; atm row moves 2.2 %
    "table2_Q_implied_by_p0": 1e-3,      # printed-constant arithmetic; atm row moves 0.7 %
    "rung1_p0_star_extrapolated": 1e-3,  # thresholds start-dependent at ~1e-8; R_t row moves 1.3 %
    "crossover_R_t": 1e-2,               # refinement bracket 2e-4 relative; 2 % row moves 2.8 %
}
# BOUNDS: identities, two-route gaps and convergence orders.  Each states its basis.
BOUNDS = {
    "U_from_2d_dev": dict(
        min=-1e-12, max=1e-12, basis="chosen",
        note="eq. (11.7.4-1) at R_t = 0.0125 m returns the printed U = 0.096 to "
             "round-off, so only a bound is stated"),
    "ode_vs_pymrm_p0": dict(
        max=1e-5, basis="chosen",
        note="two routes to rung 1's threshold: Richardson limit of the 150-4800 chain "
             "against grid-free Radau (rtol 1e-12); the limit's own error, measured "
             "4.0e-7; the bound is 25x that, and the finest grid alone is off by 2.6e-4"),
    "ode_vs_pymrm_T": dict(
        max=1e-4, basis="chosen",
        note="two routes to rung 1's hot spot at 0.018 atm, as for ode_vs_pymrm_p0; "
             "measured 3.3e-5, the bound 3x that; the 4800-cell value alone is off by "
             "about 17x the measured gap"),
    "rung1_order_p0": dict(
        min=0.9, max=1.1, basis="theory",
        note="first-order upwind on the 150-4800 chain"),
    "rung1_order_T": dict(
        min=0.85, max=1.15, basis="theory",
        note="first-order upwind; the hot spot is a root-found maximum, measured 0.94"),
    "eta_grid_vs_11_9_1_11": dict(
        max=1e-4, basis="theory",
        note="256-cell second-order pellet grid against the closed form; measured 4.6e-6"),
    "rung1_finest_err_over_start_floor": dict(
        min=1e3, basis="theory",
        note="|p0*(4800)/limit - 1| over the 1e-8 start-dependence floor of a "
             "threshold: the order fit is not fitted to that noise"),
    "eta_vs_B11_worst": dict(
        max=1e-5, basis="chosen",
        note="the same closed form against B1.1's tabulated exact eta, whose phi column "
             "is printed to 6 significant figures; measured 3.2e-6"),
}
EXACT = ("ladder_controls_rejected", "scan_points_rejected", "thresholds_rejected")
for k_ in METRICS:
    v_ = BASE_M[k_]
    if k_ in BOUNDS:
        claims.bound(k_, v_, **BOUNDS[k_])
    elif k_ in EXACT:
        continue
    else:
        claims.value(k_, v_, rtol=RTOL.get(k_, 0.05))

# THRESHOLDED RESULTS carry their margin.  The crossing control decides a verdict by
# comparing two temperatures with the ceiling and their difference with CTL_JUMP; the
# margin is the smallest of those distances, in K.  Round-off in T is about 1e-10 K.
_CTL_MARGIN = float(min(min(T_PERM - (lo + 273.15), (hi + 273.15) - T_PERM,
                            CTL_JUMP - abs(j))
                        for lo, hi, j in zip(CTL.T_just_below_C, CTL.T_just_above_C,
                                             CTL.jump_K)))
claims.exact("ladder_controls_rejected", BASE_M["ladder_controls_rejected"],
             margin=_CTL_MARGIN, min_margin=1e-6,
             note="margin in K: closest control solve to the 415 C ceiling, or jump to "
                  "the 1 K allowance, over the five ladder rungs")
claims.exact("thresholds_rejected", CTL_STATS["rejected"],
             margin=CTL_STATS["margin_K"], min_margin=1e-6,
             note="over every bisected threshold on the page, break rows included; "
                  "margin in K as for ladder_controls_rejected")
claims.exact("scan_points_rejected", SCAN_REFUSED + (len(_UNIFORM) - UNIFORM_OK),
             margin=SCAN_MARGIN_K, min_margin=1e-6,
             note="thresholds rejected in the four crossover scans (sweep points and "
                  "refinement probes) and the 9-point fine Pe_a grid; margin in K over "
                  "the four scans")

# every break row, measured above, is fed to the claims it touches
for lbl, fn in BREAK_FNS:
    for k_, v_ in BREAK_GOT[lbl].items():
        claims.breaks(lbl, k_, v_)

missing = sorted(set(METRICS) - set(claims.metrics))
assert not missing, f"metrics with a break row but no claim: {missing}"
_kinds = {kd: sorted(k_ for k_, c in claims.metrics.items() if c["kind"] == kd)
          for kd in ("value", "bound", "exact")}
print(f"{len(claims.metrics)} claims: {len(_kinds['value'])} values,"
      f" {len(_kinds['bound'])} bounds, {len(_kinds['exact'])} exact.")
for kd in ("bound", "exact"):
    for k_ in _kinds[kd]:
        c = claims.metrics[k_]
        print(f"   {kd:<5s} {k_:<34s} (here {c['value']:.2g})  {c.get('note', '')[:60]}")
print(f"Values tightened below rtol 0.05, each because a break row needs it: {RTOL}")

# PROSE: the numbers the Markdown cells quote, checked against their claims (one unit
# in the last written place).  Everything else this notebook prints is an f-string of
# the registered number itself.
SAYS = {
    "crossover_pe_axial": ["0.24883"],          # the note on locating thresholds
    "scan_points_rejected": ["0"],              # "every point of every scan passes"
}
for k_, lits in SAYS.items():
    for lit in lits:
        claims.says(lit, k_)


# %% [markdown]
# ## What pymrm adds
#
# **Honestly: not accuracy.** Every rung of this ladder is a system Froment, De Wilde
# & Bischoff solved in 1961-1967 with a Crank-Nicolson code, and section 11.7.4 says
# in as many words that "the possibilities of present-day computers are such that
# there is no longer any reason for not using the two-dimensional model". Nothing
# here is a numerical advance over the book.
#
# What pymrm adds is that **the five models become five arguments to one assembly
# rather than five programs**, and that is what makes the comparison possible at all.
# Concretely:
#
# * the ladder's five configurations differ by *which operators are built*, not by
#   which file is run: rung 2 adds a `construct_grad`/`construct_div` pair to the
#   same convection matrix, rung 5 adds a second `construct_grad`/`construct_div`
#   pair on a second axis with `nu=1`, and rungs 3 and 4 add two algebraic fields
#   masked out of the divergence by a `construct_coefficient_matrix`. Because the
#   five share the residual, the source term and the kinetics, a difference between
#   two rungs is a difference between two *models*, not between two codebases;
# * the effectiveness factor is read off **the surface-flux row of the gradient
#   operator**, which is the same quantity the film balance uses, so eq. (11.9.1-11)
#   can be checked against the equation it is a solution of rather than against a
#   volume-averaged reaction rate that would hide a boundary-condition error;
# * the geometry of a divergence is a parameter (`nu = 0, 1, 2`), so the pellet
#   sphere and the tube annulus are the same code, and the break table can inject a
#   wrong geometry into either;
# * `NumJac` makes the four-field algebraic-plus-differential system of rung 3 a
#   Newton solve of the same shape as rung 1's, and the film balances are two extra
#   rows rather than an inner iteration.
#
# Two things this page adds that the chapter does not have, both consequences of the
# above rather than of pymrm:
#
# * **the comparison is a threshold, not a profile.** The chapter compares its models
#   by plotting them; every figure in sections 11.5.2, 11.7.3 and 11.7.4 is a profile
#   or a sensitivity curve, and the one quantitative comparison it prints - "less
#   than 360" against "365" - is read off two such plots. Root-finding the inlet
#   partial pressure at which the hot spot reaches the chapter's own stated ceiling
#   turns each rung into a single number with a sign, which is what "what does this
#   rung buy" actually means to someone sizing a tube;
# * **the crossovers are located rather than asserted, and the solver's own limits
#   are located with them.** The chapter's yardsticks for climbing a rung are a rule
#   of thumb (50 particle diameters), two criteria that cannot be evaluated on its
#   own design case, and a sentence about sensitivity. Four root-found crossovers say
#   how far each coefficient would have to move before its rung mattered, which is the
#   same question with an answer - and they **rank** all four rungs, which is the part
#   a rule of thumb cannot do. Each root is located between two points the crossing
#   control validates, and each carries a grid-set error bar rather than its
#   refinement tolerance.

# %% [markdown]
# ## Reuse
#
# **What to lift.**
#
# * `Ladder` is the whole ladder. Substituting a different reaction means replacing
#   `source`; the five configurations, the boundary conditions and the diagnostics
#   carry over unchanged. `p0_critical` is the reusable part of the comparison.
# * **Locate a threshold by following the branch, not by cold starts.** Near a
#   runaway threshold there are several steady states, and Newton's damping does not
#   choose between them: which one a cold start reaches is decided by the guess, and
#   when the guess is itself a runaway profile, by round-off. The remedy is one cold solve at a
#   safe parameter value, then every scan point and every bisection midpoint started
#   from the converged state at the current lower end of the bracket. No state is
#   carried from one search to the next, so the result still cannot depend on the
#   direction a sweep was run in - the defect `B1.1` records for continuation chains.
# * **Put the control inside the root-find, not beside it.** A bisection that treats a
#   non-converged solve as "past the threshold" returns a *solver* boundary wherever
#   the solver gives out, and it does it silently. Run the $\pm 10^{-5}$ control on
#   **every** call and *raise* when it fails, so a bad threshold cannot become a
#   number anywhere, and write the searches so that a rejected point is recorded and
#   used for nothing and a scan that cannot straddle its target raises. On this page,
#   with the branch followed, the control rejects nothing; it is there for the case
#   where it would.
# * The **Danckwerts inlet as a length**: because $\varepsilon D_{ea}/u_s$ and
#   $\lambda_{ea}/(\rho_g c_p u_s)$ are both lengths, the inlet condition of rung 2
#   is `{"a": d_p/Pe, "b": 1.0, "d": c_in}` with no unit bookkeeping at all, and the
#   two fields differ only in their Péclet number. That is `A2.1`'s dict with the
#   scaling done once.
# * **Scale an algebraic residual to the units of its own unknown.** Written in the
#   units the book writes them in, (11.8.1-3) is in atm/h and (11.8.1-4) in
#   kcal/(m$^3$ h), so a residual test with one absolute threshold, like the one in
#   `Ladder.solve`, is strict for one row and loose for the other. Dividing each
#   balance by its own transfer coefficient puts both rows in the units of the
#   surface unknown they determine.
# * The **helper-plus-injection pattern** in the break table: each metric has exactly
#   one function that computes it, *every* metric's helper default is asserted equal
#   to a number printed earlier in the notebook to 1e-12, and every break row
#   perturbs one of its arguments. That is what makes "the row reaches the code it
#   perturbs" checkable rather than a claim: a row that reads a module constant
#   instead of the configuration does not reach the code it names. Two further rules:
#   * **an unasserted helper drifts.** A helper that runs its own coarser refinement
#     while the notebook displays a finer one carries a different number into
#     `agreement.json`, and nothing compares them. The remedy is structural rather
#     than a tighter tolerance - the display and the metric are the *same call*,
#     memoised - and the assertion list is checked against the metric set, so a
#     metric cannot be added without one.
#   * **the four shifts are the exception, and exceptions get printed.** The break
#     rows cannot afford the displayed grid once per row, so those four
#     metrics are computed on a coarser one. Both values and their gap are printed
#     in a table; an exception you can see is a caveat, an exception you cannot is a
#     discrepancy.
# * **The AST guard, but only if you measure its teeth on the forms you care about
#   AND publish the ones it misses.** This one parses each row plus one level of the
#   helpers it names and rejects a metric key bound to a numeric literal - and it
#   constant-folds, so `0.1656 + 0.0001`, `float("15.688")`, a subscript off a frozen
#   tuple, a one-line laundering helper and tuple unpacking are all caught. A
#   negative control that tests only the bare literal measures the guard's easiest
#   case rather than its coverage, so there is one control per laundering form, each
#   asserted to catch exactly its own keys - and a list of ESCAPE forms, printed and
#   asserted to escape (the guard cell prints both counts), so the notebook states
#   the guard's blind spots rather than implying it has none. A module-level
#   constant is one of them, and it is the realistic accidental case. The
#   false-positive control binds a *metric* key to a computed value; a non-metric
#   key is filtered out before the guard looks at anything, and would measure
#   nothing.
#
# **What not to lift.**
#
# * **The constants, without the audit.** Section 11.5.2's printed specific heat is
#   several per cent from what its own downstream groups require, its printed mean
#   density is out by $10^3$, group $A$ of Example 11.5.3.A misses its own printed
#   formula, its tube radius is not the one its own figures were computed at, and its
#   pressures are in a different unit from the one it prints. Every one of those was
#   found by *loading `D2.2`'s and `C2.10`'s CSVs and printing them beside the
#   book's*, which took an afternoon and is the single highest-yield thing on this
#   page. If you retype a constant from a textbook restatement of a paper, check the
#   restatement. (All the percentages are computed in *Parameters* and *Validation*;
#   none is typed here.)
# * **A row-by-row transcription check, on its own.** Every row of Table 11.5.3.A-2
#   transcribes from the 1970 paper at a ratio of exactly 1.000000, and the table is
#   still internally inconsistent by about a third, because one of those rows was
#   copied correctly *from a different worked example at a different tube radius*.
#   What caught it was a CROSS-COLUMN test the table imposes on itself -
#   $\Delta T_{ad} = (B/A)p_0$ with the example's own printed $A$ and $B$. If a table
#   has two columns computed from one quantity, check them against each other, not
#   only against the source. **And then check what the cross-column test cannot
#   see**: Table 11.5.3.A-1 *passes* it while carrying the same transplanted $Q$,
#   because both of its columns descend from that $Q$. A consistency test between two
#   numbers can only catch a defect that moved one of them.
# * **The shifts, as advice about a real o-xylene reactor.** The base case is section
#   11.5.2's deliberate simplification, and the notebook shows it is $k_1$ alone used
#   as the total consumption rate, with the combustion route's heat of reaction
#   (printed in the reconciliation cell) absent. The rung-to-rung *ordering* is what transfers; the magnitudes are for
#   this pseudo-first-order reactor.
# * **Any of it as validation.** Chapter 11 reports no experiment. Everything here is
#   a reproduction of the book's own arithmetic, of its transcription of a 1970 paper,
#   and of five claims it makes about its own models.
#
# **What the page cannot conclude.**
#
# * Whether rung 4 is needed on the real catalyst. The crossover $D_e$ is computed,
#   but the chapter prints no $D_e$ for V$_2$O$_5$ on silica and sends the reader to
#   two documents that are not on disk, so the sign of the rung-4 correction is
#   established and its size is not.
# * Whether $\lambda_{ea}$ makes rung 2 matter. The chapter says "little information
#   is available" and never defines $\mathrm{Pe}_{ha}$; the page sweeps the two Péclet
#   numbers together, so the crossover it root-finds is a crossover in *both* at once
#   and not in $\lambda_{ea}$ alone. Nor is the shift below the swept range reported:
#   the crossover cell prints the largest $|$shift$|$ on that sweep, and nothing is
#   claimed about smaller Péclet numbers.
# * Anything about section 11.10, the sixth cell of Table 11.4-1 - the
#   two-dimensional heterogeneous model. It is out of scope, and the Papageorgiou &
#   Froment [1995] kinetics it uses are a different reaction model.
# * Anything read off a figure. Figures 11.5.2-1/2, 11.6-1/2, 11.7.1-1/2/3, all six
#   11.7.3 figures, 11.7.4-1, 11.8.1-1..4, 11.8.2-1, 11.9.1.A-1..8 and the runaway
#   band of Fig. 11.5.3-1 are not digitised, not traced and not read. The **one**
#   thing taken from a figure anywhere on this page is the *unit word* on
#   Fig. 11.5.2-1's ordinate label, "$p_0$, atm", read on a 300 ppi crop; no curve,
#   tick or point is read from it or from any other.
# * **How rung 5's threshold converges.** Its hot spot is refined on both axes, but
#   the threshold - the quantity the ladder table actually reports for that rung - is
#   refined on neither. The (600, 10)-versus-(1200, 20) gap is printed; no observed
#   order and no Richardson limit is claimed for it, and the tube-radius crossover is
#   located at (600, 10) with that gap, divided by the local slope of the sweep,
#   as its error bar - three figures on $R_t^*$, not the six its bisection tolerance
#   would suggest.
# * **Whether a state exists that the Newton cannot reach.** `p0_critical` counts a
#   non-converged solve as supercritical. A control shows that every threshold this
#   page reports - not only the five default rungs - is a genuine 415 $^\circ$C crossing rather than
#   a convergence boundary, and thresholds that fail it are raised on rather than
#   returned. What the control cannot do is find a bounded state the Newton never
#   reaches: section 11.6's own three-steady-state discussion is a warning that a
#   computation can miss a state that exists. Each threshold here is, by
#   construction, the one on the branch followed up from a safe feed; other branches
#   are not searched for.

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
