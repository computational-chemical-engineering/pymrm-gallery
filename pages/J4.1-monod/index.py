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
# title: "Monod: what a growth curve can and cannot pick out of the family"
# description: "Froment, De Wilde & Bischoff section 1.5.2 prints the Monod equation as eq. (1.5.2-1), defines its two constants by exactly two properties - the asymptote r_m and the half-rate point K_S - and couples it to the biomass balance eq. (1.5.2-2); its Chapter 1 reference list carries the exact 1949 citation this case catalogues. Rawlings & Ekerdt p. 596 prints the same equation beside four alternatives - Blackman, Tessier, Moser, Contois - and calls it 'the simplest form of the Langmuir adsorption isotherm'. This page asks whether data can pick Monod out of that family, on the one batch run Levenspiel's Problem 29.18 reprints as Monod's own reported measurements, and answers with numbers. Two of the four alternatives are not distinguishable AT ALL by a constant-yield batch experiment: Moser contains Monod at n = 1, and CONTOIS IS EXACTLY MONOD REPARAMETERISED once the yield fixes X as an affine function of S - proved in sympy, and confirmed by two independent fits that land on the same curve, agreeing to better than 1e-3 in both constants (the exact residuals are optimiser exit tolerances, checked on every run as bounds rather than compared as values). Of the two that are, Tessier needs 10.3 replicates of this design and Blackman 1.30. And the ranking of the five on this run is decided by the DATA REDUCTION, not by the data: fit as the exercise instructs, and Monod is last of the four two-parameter laws and Moser's exponent n = 1.93 clears the 95 % F test - a verdict that is CO-CAUSED by the one printed cell the page flags, because at that cell's balance-implied value route A's F falls to 5.45, below the same 7.71. Do the interval average exactly instead - same seven growth rates, same models, no new parameters - and Monod's residual falls by 42 %, n drops to 1.51, and F falls from 12.3 to 1.44. A seeded parametric bootstrap prices that, with every draw's fit REFINED off the search grid rather than read from it: on data the Monod model itself generates, the reduction rejects Monod 36.5 % +- 1.1 % of the time against a nominal 5 %, and drives the median fitted exponent to 1.61 +- 0.0074 - while the same test on data drawn from route A's OWN model, the control, sits AT nominal size, 4.95 % +- 0.5 % against 5.00 %, and returns 1.01 +- 0.0049 (those six numbers are Monte-Carlo estimates at 2000 draws - FOUR FREQUENCIES AND TWO MEDIANS, not six frequencies - and every quotation of them carries its sampling error: at that draw count the four frequencies are seed-independent in their leading figure only, the two medians to at least two figures). So the small-sample F approximation at n = 7 is essentially exact and is not the culprit; the reduction is, and the observed F = 12.3 has p = 0.23 +- 0.0094 against a properly generated Monod null. One printed cell is reported and not repaired: the table's fourth interval-average substrate concentration is off the constant-yield line the other six rows fix to R^2 = 0.99996. Monod = Langmuir is proved as an identity; Monod and Blackman are shown to share BOTH of Froment's defining properties exactly while differing by exactly r_m/3 at the switch; and the experiment that does separate Monod from Contois is named and priced - a chemostat feed sweep, where the critical dilution rates cross exactly at the batch run's own substrate intercept and differ by a factor 1.83 at S_f = 20."
# categories: [sec:J, struct:S1, tier:T0, data:tier2, phase:liquid]
# date: 2026-08-13
# ---
#
# # Monod: what a growth curve can and cannot pick out of the family
#
# **Catalog ID:** `J4.1` · **Structures:** `S1` · **Tier:** T0
#
# ## Background
#
# **The origin of this case was not consulted.** Monod, J., *The growth of
# bacterial cultures*, **Annu. Rev. Microbiol. 3**, 371-394 (1949) is not on disk
# and could not be obtained. This page is built under the gallery's
# textbook-canonical-source rule from a monograph that states the model,
# attributes it, and carries it:
#
# > **Froment, G. F., De Wilde, J. & Bischoff, K. B.**, *Chemical Reactor Analysis
# > and Design*, 3rd edn, John Wiley & Sons (2011), ISBN 978-0-470-56541-4,
# > **section 1.5.2 "Microbial Kinetics", book pp. 26-29** (PDF 68-71).
#
# Section 1.5.2 prints the rate law as eq. (1.5.2-1), names $K_S$ *"the Monod- or
# saturation constant"*, records that *"in the original Monod equation, $r$ is
# represented by $\mu$, $r_m$ by $\mu_{max}$ and $C_A$ by $S$"*, and couples it to
# the biomass balance eq. (1.5.2-2) and the yield relation eq. (1.5.2-3). Its
# Chapter 1 reference list, on book p. 58, carries
#
# > Monod, J., *Annu. Rev. Microbiol.*, **3**, 371 (1949).
#
# - the exact citation this case catalogues, read at digit scale on a 300 ppi
# crop. So the origin is **identified** by the book rather than implied, which is
# what the gallery's E1.1 test asks for. **Everything this page says about Monod
# is bounded by what Froment prints about him.** There is no sentence of the form
# "Monod showed", "Monod assumed" or "Monod found" anywhere on this page, in
# `meta.yaml` or in `README.md`; that sweep was done deliberately, because
# attribution drift is the defect class that has cost this repository two pages.
#
# **The question this page is for** comes from a second book. **Rawlings, J. B. &
# Ekerdt, J. G.**, *Chemical Reactor Analysis and Design Fundamentals*, 2nd edn
# (2025 printing), **book p. 596**, prints the Monod equation as one of five:
#
# | | as printed on p. 596 |
# |---|---|
# | Monod | $\mu = \mu_m S/(K_s + S)$ |
# | Blackman | $\mu = \mu_m$ for $S \ge 2K_s$; $\mu = \mu_m S/(2K_s)$ for $S < 2K_s$ |
# | Tessier | $\mu = \mu_m(1 - e^{-K_s S})$ |
# | Moser | $\mu = \mu_m S^n/(K_s + S^n)$ |
# | Contois | $\mu = \mu_m S/(K_{sx} X + S)$ |
#
# and adds: *"We recognize the Monod equation [12] for cell growth rate as the
# simplest form of the Langmuir adsorption isotherm and the resulting
# Hougen-Watson kinetics for reaction rates on catalyst surfaces discussed in
# Chapter 5."* Its reference [12] is the same 1949 paper, with the page range.
# Rawlings & Ekerdt is used here as a **model catalogue and for its chemostat
# balances**, not as a second source for the rate law - the E1.1 target was pinned
# to Froment before this page was started, and nothing is adjudicated between the
# books.
#
# **So: can data pick Monod out of that family?** The answer below is in four
# parts - two of the four alternatives are not separable by a batch experiment at
# all and this is exact; the other two are separable and the page prices them; the
# apparent answer on the one real dataset available is set by the *data reduction*
# rather than by the data; and the experiment that would settle the remaining
# question is a different one.
#
# **The data.** Levenspiel, O., *Chemical Reaction Engineering*, 3rd edn (1999),
# **Problem 29.18, book p. 644**, prints seven time intervals of a batch run and
# introduces them as *"the reported data for one of his runs"*, attributing them
# to **"Monod, 1958, p. 74"** - the second edition of *Recherches sur la
# Croissance des Cultures Bacteriennes*, **not** the 1949 paper this case
# catalogues. That book is not on disk either and was not consulted. Levenspiel's
# own account, quoted as his, is that Monod proposed the equation *"in 1939, as
# part of his thesis"*, that *"the thesis was published as a book in 1948"*, and
# that it *"was later condensed and translated into English in 1949"*. **Which of
# those documents actually prints this table is not established here.**
#
# **Scope.** Substrate inhibition - Froment's eq. (1.5.2-4), attributed there to
# Graef and Andrews [1973], and Rawlings' $\mu = \mu_m S/(K_s + S + K_1S^2)$ - is
# **case J4.2** and gets no quantitative treatment on this page. Williams'
# structured model, eqs. (1.5.2-5) to (1.5.2-8) in the same section, is J4.9's.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.20"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "J4.1-monod"  # equals the page directory

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
# Imports. sympy is a Colab default; everything else is pymrm, numpy, scipy,
# pandas and matplotlib.
import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, minimize, minimize_scalar
from scipy.stats import f as f_dist, chi2

from pymrm.checks import find_roots, observed_orders
from gallery_utils import Claims, load_data, load_meta, cite_data
from model import (GROWTH, BatchMarch, PlugFlowFermenter, Chemostat, batch_reference,
                   integrate_interval, one_root)

claims = Claims(PAGE.split("-")[0])
plt.rcParams.update({"figure.dpi": 110, "font.size": 9.5})
# Okabe-Ito, assigned in fixed order and never cycled
C_BLUE, C_ORANGE, C_GREEN = "#0072B2", "#D55E00", "#009E73"
C_PURPLE, C_YELLOW, C_GREY = "#CC79A7", "#E69F00", "0.45"
COL = {"Monod": C_BLUE, "Blackman": C_ORANGE, "Tessier": C_GREEN,
       "Moser": C_PURPLE, "Contois": C_YELLOW}
np.set_printoptions(precision=6, suppress=False)

# DETERMINISM: two consecutive executions must give identical content, so
# (i) every stochastic step is seeded from SEED below and (ii) the whole of a
# long cell's output is batched into one message rather than flushed on a timer.
SEED = 20260813
try:
    sys.stdout.flush_interval = 1e6
except Exception:
    pass
from IPython.display import HTML, display


def show(styler):
    display(HTML(styler.to_html()))


M = {}          # -> agreement.json
BREAK = []      # defect-injection table
print("setup done")

# %% [markdown]
# ## The published model
#
# Froment eq. (1.5.2-1), verbatim, is
#
# $$r = \frac{r_m C_A}{K_S + C_A} \tag{1.5.2-1}$$
#
# with $r$ *"a fractional or specific growth rate of biomass"*, $r_m$ *"the
# maximum specific rate of biomass growth, i.e., the rate when the substrate
# concentration is not limiting ($C_A$ much larger than $K_S$)"*, and $K_S$
# *"the concentration of substrate $A$ at which the rate equals half of its
# maximum value, $r_m$"*. It is then coupled to the biomass balance
#
# $$r_B = \frac{{\rm d}C_B}{{\rm d}t} = r\,C_B = \frac{r_m C_B C_A}{K_S + C_A}
# \tag{1.5.2-2}$$
#
# and closed by a constant yield,
#
# $$-\frac{{\rm d}C_A}{{\rm d}t} = \frac{1}{y_B}\frac{{\rm d}C_B}{{\rm d}t}.
# \tag{1.5.2-3}$$
#
# **A note on the text layer, corrected.** This book is born-digital and its words
# survive extraction, but its equations do not: **every Symbol-font operator in
# them extracts as an unmappable Private-Use-Area glyph**, `U+F02D` for the minus,
# `U+F03D` for the equals, `U+F02B` for the plus. On eq. (1.5.2-3) all three modes
# of `pdftotext` emit the leading minus as `U+F02D`, which renders as nothing or
# as tofu - so the extracted line *looks* like `dC A dC B / dt y B dt`, an
# equation with no operators at all, and a reader who copies what is visible gets
# a sign error for free. The character is therefore **not dropped and the failure
# is not specific to the sign**: the `=` on the same line and the `=` and `+` of
# eq. (1.5.2-1) go the same way. The 300 ppi render carries all of them, and every
# transcription on this page was read there. This corrects the diagnosis the
# just-published `J4.6` page gives for eq. (1.5.1-17) on p. 26 of this same book -
# its *"the extractor dropped the leading term"* is the same phenomenon read the
# same wrong way, since `-layout` returns that line as
# `= C A <U+F02D> E <U+F02B> ...` with the leading term present. The **lesson** in
# both places is right and unchanged: read the render, never the text layer.
#
# **Five identities come first, before any data.** Three of them decide what an
# experiment can possibly do, and all five are proved symbolically below rather
# than asserted:
#
# 1. **Monod is the Langmuir isotherm**, which is what Rawlings & Ekerdt say on
#    p. 596. $\theta = KC/(1+KC) = C/(K^{-1}+C)$, so $K_S = 1/K$ and nothing else
#    changes. The gallery has that isotherm at `J1.1` and the enzyme form at
#    `J4.6`; they are one algebraic family.
# 2. **Monod and Blackman share *both* of the properties Froment uses to define
#    $r_m$ and $K_S$.** Blackman's law reaches $\mu_m$ exactly, and at $S = K_S$
#    it gives exactly $\mu_m/2$ - because $K_S < 2K_S$ puts that point on its
#    linear branch, where $\mu = \mu_m S/(2K_S)$. So the asymptote and the
#    half-rate point, the two properties by which the constants are defined, do
#    **not** distinguish the two laws at all. What does: at the branch switch
#    $S = 2K_S$ the two differ by **exactly $\mu_m/3$**, and on the lower branch
#    the largest gap is $\mu_m(3-2\sqrt2)/2 = 0.0858\,\mu_m$ at
#    $S = (\sqrt2-1)K_S$. Both are closed form and both are checked numerically
#    against a root-found maximum.
# 3. **Moser contains Monod** at $n = 1$. The two are nested, so no data can
#    *reject* Moser; the only meaningful question is whether $n$ differs from 1,
#    and that is an F test.
# 4. **Contois *is* Monod on any constant-yield batch run** - and this is the
#    result that decides most of the page. Constant yield makes the biomass an
#    affine function of the substrate, $X = \varphi(a - S)$, so
#    $$\mu = \frac{\mu_m S}{K_{sx}X + S}
#         = \frac{\mu_m S}{K_{sx}\varphi a + (1-K_{sx}\varphi)S}
#         = \frac{(\mu_m/\beta)\,S}{K_{sx}\varphi a/\beta + S},
#      \qquad \beta = 1 - K_{sx}\varphi,$$
#    which is Monod with $\mu_m' = \mu_m/\beta$ and $K_S' = K_{sx}\varphi a/\beta$.
#    The map is a bijection. **No batch experiment, of any length or precision,
#    separates Contois from Monod.** The symbolic residual below is exactly zero,
#    and two independent fits confirm it by landing on the same curve.
# 5. **Tessier is genuinely different** - it has no reparameterisation into
#    Monod - so it is the one alternative a batch run could in principle settle,
#    and the page prices what that would take.

# %%
# ---- five identities, proved not asserted
S_, X_, mum_, Ks_, Ksx_, phi_, a_, n_, K_ = sp.symbols(
    "S X mu_m K_s K_sx varphi a n K", positive=True)

monod    = mum_*S_/(Ks_ + S_)
blackman_lo = mum_*S_/(2*Ks_)                 # branch S < 2 K_s
tessier  = mum_*(1 - sp.exp(-Ks_*S_))
moser    = mum_*S_**n_/(Ks_ + S_**n_)
contois  = mum_*S_/(Ksx_*X_ + S_)
langmuir = K_*S_/(1 + K_*S_)                  # theta, J1.1's eq. (1)

I = {}
# 1. Monod == Langmuir with K_s = 1/K   (mu_m = 1 for the isotherm)
I["monod_is_langmuir"] = sp.simplify(langmuir - monod.subs({Ks_: 1/K_, mum_: 1}))
# 2a. Froment's two defining properties, for Monod
I["monod_asymptote"]   = sp.simplify(sp.limit(monod, S_, sp.oo) - mum_)
I["monod_half_rate"]   = sp.simplify(monod.subs(S_, Ks_) - mum_/2)
# 2b. the SAME two properties, for Blackman
I["blackman_asymptote"] = sp.simplify(mum_ - mum_)          # upper branch IS mu_m
I["blackman_half_rate"] = sp.simplify(blackman_lo.subs(S_, Ks_) - mum_/2)
# 2c. where they part: the gap at the switch, and the lower-branch extremum
gap_switch = sp.simplify((mum_ - monod).subs(S_, 2*Ks_))    # upper branch minus Monod
S_star     = sp.solve(sp.Eq(sp.diff(blackman_lo - monod, S_), 0), S_)[0]
gap_lower  = sp.simplify((blackman_lo - monod).subs(S_, S_star))
I["blackman_gap_at_switch"] = sp.simplify(gap_switch - mum_/3)
I["blackman_gap_lower"]     = sp.simplify(gap_lower - mum_*(2*sp.sqrt(2) - 3)/2)
# 3. Moser -> Monod at n = 1
I["moser_at_n1"] = sp.simplify(moser.subs(n_, 1) - monod)
# 4. Contois == Monod under a constant yield X = phi (a - S)
beta_    = 1 - Ksx_*phi_
contois_b = contois.subs(X_, phi_*(a_ - S_))
I["contois_is_monod"] = sp.simplify(
    contois_b - (mum_/beta_)*S_/(Ksx_*phi_*a_/beta_ + S_))
# 5. Tessier has no such reduction: its curvature ratio is not Monod's anywhere
#    (a POSITIVE statement, not an absence: the two second derivatives differ by
#     a factor that depends on S for every choice of the constants)
tess_ratio = sp.simplify(sp.diff(tessier, S_, 2)/sp.diff(tessier, S_))
mon_ratio  = sp.simplify(sp.diff(monod, S_, 2)/sp.diff(monod, S_))
I["tessier_not_monod"] = sp.simplify(sp.diff(tess_ratio - mon_ratio, S_))

for k, v in I.items():
    if k == "tessier_not_monod":
        continue
    if v != 0:
        raise RuntimeError(f"identity {k} does not close: residual {v}")
M["sym_monod_is_langmuir_residual"]    = float(abs(I["monod_is_langmuir"]))
M["sym_monod_half_rate_residual"]      = float(abs(I["monod_half_rate"]))
M["sym_blackman_half_rate_residual"]   = float(abs(I["blackman_half_rate"]))
M["sym_moser_at_n1_residual"]          = float(abs(I["moser_at_n1"]))
M["sym_contois_is_monod_residual"]     = float(abs(I["contois_is_monod"]))


def sym_gaps(switch=2):
    """Closed-form gaps over mu_m of a Blackman law switching at S = switch K_S:
    at the switch, and at the lower-branch extremum (sympy, exact)."""
    lo = mum_*S_/(switch*Ks_)
    g_sw = sp.simplify((mum_ - monod).subs(S_, switch*Ks_)/mum_)
    s_st = sp.solve(sp.Eq(sp.diff(lo - monod, S_), 0), S_)[0]
    return g_sw, sp.simplify(sp.Abs((lo - monod).subs(S_, s_st))/mum_)


GAPS2 = sym_gaps(2)
if sp.simplify(GAPS2[0] - sp.Rational(1, 3)) != 0 or sp.simplify(GAPS2[1] - (3 - 2*sp.sqrt(2))/2) != 0:
    raise RuntimeError(f"closed-form gaps {GAPS2} are not 1/3 and (3 - 2 sqrt 2)/2")
M["sym_blackman_gap_at_switch_over_mum"] = float(GAPS2[0])
M["sym_blackman_gap_lower_over_mum"]     = float(GAPS2[1])
S_STAR_OVER_KS = float(S_star/Ks_)
print("all five identities close exactly (sympy residual 0)")
print("  Monod == Langmuir with K_S = 1/K")
print("  Monod and Blackman: asymptote mu_m and half-rate at S = K_S, both exact")
print("  gap at the switch S = 2K_S            : mu_m *", sp.nsimplify(sp.Rational(1, 3)),
      "=", f"{M['sym_blackman_gap_at_switch_over_mum']:.6f} mu_m")
print("  largest gap on the lower branch       : mu_m *",
      sp.simplify((3 - 2*sp.sqrt(2))/2), "=",
      f"{M['sym_blackman_gap_lower_over_mum']:.6f} mu_m at S/K_S ="
      f" {S_STAR_OVER_KS:.6f}")
print("  Moser(n=1) == Monod;  Contois == Monod under a constant yield")
print("  Tessier: d/dS of the curvature ratio difference is not identically zero ->",
      sp.simplify(I["tessier_not_monod"]) != 0)

# %%
# The two Blackman gaps are CLOSED FORM above; here they are checked against a
# ROOT-FOUND extremum on the real line, because a closed form with an algebra
# slip in it looks exactly like a correct one.
mum_v, Ks_v = 1.0, 1.0


def blackman_numeric(switch=2.0, at=None, refine=True):
    """Blackman (switch at S = switch K_S) against Monod, numerically: the gap at
    S = at (the switch by default), and the lower-branch extremum, the one root of
    the analytic derivative on 0 < S < switch K_S (refine=False reads it off a
    64-point scan instead)."""
    f_lo = lambda S: mum_v*S/(switch*Ks_v) - mum_v*S/(Ks_v + S)
    d_lo = lambda S: mum_v/(switch*Ks_v) - mum_v*Ks_v/(Ks_v + S)**2
    if refine:
        s_hat = one_root(d_lo, 0.0, switch*Ks_v, xtol=1e-14, what="lower-branch extremum")
    else:
        grid = np.linspace(0.0, switch*Ks_v, 64)
        s_hat = float(grid[np.argmax(np.abs(f_lo(grid)))])   # the scan read-off (break row)
    S_at = switch*Ks_v if at is None else at
    gap_at = abs(mum_v*min(1.0, S_at/(switch*Ks_v)) - mum_v*S_at/(Ks_v + S_at))
    return s_hat, abs(f_lo(s_hat)), gap_at


S_hat, gap_lo_num, gap_sw_num = blackman_numeric()
M["num_blackman_gap_lower_relerr"] = abs(
    gap_lo_num - M["sym_blackman_gap_lower_over_mum"])/M["sym_blackman_gap_lower_over_mum"]
M["num_blackman_gap_switch_relerr"] = abs(
    gap_sw_num - M["sym_blackman_gap_at_switch_over_mum"])/M["sym_blackman_gap_at_switch_over_mum"]
M["num_blackman_Sstar_over_Ks_relerr"] = abs(S_hat - S_STAR_OVER_KS)/S_STAR_OVER_KS
# all three are round-off or the root's xtol, so each is printed against its bound
print(f"root-found S*/K_S = {S_hat:.12f} against the closed form {S_STAR_OVER_KS:.12f}"
      f"  (rel. below 1e-13, here {M['num_blackman_Sstar_over_Ks_relerr']:.1e})")
print(f"lower-branch gap  {gap_lo_num:.12f} mu_m   (rel. below 1e-14, here"
      f" {M['num_blackman_gap_lower_relerr']:.1e})")
print(f"gap at the switch {gap_sw_num:.12f} mu_m   (rel. below 1e-14, here"
      f" {M['num_blackman_gap_switch_relerr']:.1e})")
print("\nSo the two laws agree EXACTLY on both properties Froment defines the")
print("constants by, and differ by a third of the maximum rate at S = 2 K_S.")

# %% [markdown]
# ## Parameters and assumptions
#
# **Nothing on this page is a literature parameter value.** Section 1.5.2 of
# Froment prints **no value for $r_m$, $K_S$ or $K_i$ and no worked example for
# Monod kinetics** - no table, no fitted constant, nothing to reproduce. (It is
# not numberless: book p. 26 prints the selectivity identity
# $y_B + y_P + y_D = 1$, and the caption of Fig. 1.5.2-1 on book p. 29 prints
# $k_1 = 0.0125\,$s$^{-1}$, $k_2 = 0.025\,$s$^{-1}$, $M_0/D_0 = 1$ - but those
# belong to **Williams'** structured model, eqs. (1.5.2-5) to (1.5.2-8), which is
# `J4.9`'s, and none of them is a Monod constant. Both are transcribed into the
# printed-laws CSV so the claim is checkable.) Rawlings & Ekerdt's set
# ($\mu_m = 1$, $K_s = 1$, $S_f = 5$, $y = 1$) is printed on their book p. 597 as
# the values *"used to prepare"* Fig. 10.12 - body text, not a caption - and is
# used here only to reproduce their eq. (10.19) and nothing else. Every constant that describes the culture is
# fitted here, on this page, to the seven printed intervals, and is labelled a fit
# wherever it appears.
#
# **Units.** $\Delta t$ is printed in hours, so every growth rate here is in
# h$^{-1}$ - which is the dimension Froment prints for $\mu$ on p. 27, and the one
# independent unit check available. **The substrate and biomass columns carry no
# units at all** in the printed table. Every conclusion below is therefore stated
# either in h$^{-1}$ or as a ratio; the one quantity that inherits the unprinted
# scale is $K_S$ itself, and it is quoted as *"in the units of the printed
# $\bar C_A$ column"* every time.
#
# **The assumptions, and who prints each:**
#
# | assumption | printed by | where it bites |
# |---|---|---|
# | constant biomass yield | Froment eq. (1.5.2-3); Rawlings eq. (10.18); Levenspiel Ch. 29 | fixes $S$ as an affine function of $X$ - this is what makes Contois degenerate |
# | well mixed, isothermal batch | Froment "a closed volume with uniform composition and temperature" | the whole reduction |
# | exponential-phase growth only | Froment: *"The rate of biomass production, i.e., the increase in the number of cells during the exponential growth phase, is described by the empirical Monod-kinetics"* | the seven intervals are taken to be inside it |
# | no cell death | Froment prints $r_d = k_dC_B$ separately and it is not used here | would add a constant to every $\mu$ |
# | no product inhibition | Rawlings prints the form separately (p. 596); Levenspiel gives it a whole chapter (Ch. 30) | out of scope, and named |
#
# **One reduction, and it is exact.** For any growth law whatever,
# $\mathrm{d}(\ln C_C)/\mathrm{d}t = \mu$, so
#
# $$\mu_i \;=\; \frac{1}{\Delta t_i}\ln\frac{C_{C,i}^{\rm end}}{C_{C,i}^{\rm start}}$$
#
# is **exactly** the time-average of the specific growth rate over interval $i$.
# No approximation has been made at this point.
#
# **The approximation comes next, and it is the one this page is about.** The
# exercise invites you to compare that $\mu_i$ with $\mu(\bar C_{A,i})$ - the rate
# law evaluated at the interval's *average* substrate concentration. Since $\mu$
# is curved, $\langle\mu(C_A)\rangle \ne \mu(\langle C_A\rangle)$, and the two
# differ at second order in the substrate swing across the interval. That swing is
# not small here: the last two intervals consume more substrate than their own
# average. The page calls the first comparison **route A** and the exact one
# **route B**, and reports both.

# %% [markdown]
# ## The data
#
# Seven intervals, transcribed from a 600 ppi native render of Levenspiel p. 644
# and used exactly as printed. `data/levenspiel-p644-monod-batch.csv` carries them
# and its sidecar carries the provenance chain in full; the essentials are in
# *Background* above and are repeated in one line here: **these are attributed by
# Levenspiel to Monod (1958) p. 74, not to the 1949 paper, and neither document
# was consulted.**
#
# Two things the printed table says that the text layer does not:
#
# - the substrate column header carries an **overbar** - $\bar C_A$, an interval
#   average. `pdftotext` drops it in all three modes. It is load-bearing;
# - **no units** are printed for either concentration column.
#
# **A material balance the table must satisfy.** All three books assume a constant
# yield, which makes $C_A$ an affine function of $C_C$ *pointwise*, hence $\bar
# C_A$ an affine function of $\bar C_C$ *interval by interval*. Under exponential
# growth the exact time-average of $C_C$ is its **logarithmic** mean, so that line
# can be drawn with no fitted kinetics at all - and it is the only internal
# consistency check the table offers.

# %%
DATA = load_data("levenspiel-p644-monod-batch.csv", page=PAGE)
META = load_meta("levenspiel-p644-monod-batch.csv", page=PAGE)
PRINTED = load_data("printed-growth-laws.csv", page=PAGE)
P = {r.key: r.as_printed for r in PRINTED.itertuples()}

dt   = DATA["dt_hr"].to_numpy(float)
CAb  = DATA["CA_bar"].to_numpy(float)
Cc0  = DATA["CC_start"].to_numpy(float)
Cc1  = DATA["CC_end"].to_numpy(float)
NOBS = len(dt)
LN   = np.log(Cc1/Cc0)
mu   = LN/dt                                   # EXACT interval-average growth rate
Xlm  = (Cc1 - Cc0)/LN                          # exact time-average of C_C
Xam  = 0.5*(Cc0 + Cc1)                         # the other convention (break row)
RES_CA, RES_CC, RES_DT = 0.5, 0.05, 0.005      # printed half-ulps

print(cite_data(META))
print(DATA.to_string(index=False))
print("\n  mu_i = ln(Cc_end/Cc_start)/dt, 1/h :", np.round(mu, 5))
print("  log-mean biomass over the interval :", np.round(Xlm, 4))
print("\n  elapsed time of the whole run: %.2f h in %d intervals" % (dt.sum(), NOBS))


# %%
# ---- the constant-yield line, and the one cell that is not on it
def bal_line(drop=()):
    keep = np.array([i for i in range(NOBS) if i not in drop])
    A = np.vstack([np.ones(len(keep)), Xlm[keep]]).T
    c, *_ = np.linalg.lstsq(A, CAb[keep], rcond=None)
    r = CAb[keep] - A @ c
    R2 = 1 - np.sum(r**2)/np.sum((CAb[keep] - CAb[keep].mean())**2)
    return float(c[0]), float(c[1]), float(R2), r

a7, b7, R2_7, r7 = bal_line()
loo = {i: bal_line(drop=(i,)) for i in range(NOBS)}
worst = max(range(NOBS), key=lambda i: loo[i][2])          # the row whose removal helps most

A_BAL, B_BAL, R2_BAL, _ = loo[worst]
PHI_BAL = -1.0/B_BAL                     # yield: biomass formed per substrate consumed
implied = A_BAL + B_BAL*Xlm[worst]
end_of_interval = A_BAL - Cc1[worst]/PHI_BAL

lines = ["constant-yield line  CA_bar = a + b * Cc_logmean",
         f"  all {NOBS} rows      : a = {a7:10.4f}   b = {b7:8.5f}   R2 = {R2_7:.8f}",
         "  leave-one-out R2 (the check that says WHICH row, not merely that one is off):"]
for i in range(NOBS):
    lines.append(f"     without row {i+1}: R2 = {loo[i][2]:.8f}"
                 + ("   <-- by far the largest" if i == worst else ""))
lines += [f"\n  row {worst+1} is the one: dropping it takes R2 from {R2_7:.6f} to {R2_BAL:.8f};",
          f"  dropping any other row leaves R2 in "
          f"[{min(loo[i][2] for i in range(NOBS) if i != worst):.6f},"
          f" {max(loo[i][2] for i in range(NOBS) if i != worst):.6f}].",
          f"\n  the six-row line predicts CA_bar(row {worst+1}) = {implied:.4f};"
          f" the table prints {CAb[worst]:.0f}.",
          f"  residual {CAb[worst]-implied:+.4f}, against"
          f" {np.max(np.abs(loo[worst][3])):.4f} for the largest of the other six.",
          f"\n  yield from the six-row line : phi = {PHI_BAL:.6f} biomass per substrate",
          f"  substrate intercept         : a   = {A_BAL:.4f} (C_A when C_C -> 0)",
          f"  substrate at the run's start: {A_BAL - Cc0[0]/PHI_BAL:.4f}"
          f" ; at its end: {A_BAL - Cc1[-1]/PHI_BAL:.4f}"]
M["balance_R2_all_rows"]   = R2_7
M["balance_R2_six_rows"]   = R2_BAL
M["balance_row4_residual"] = float(CAb[worst] - implied)
M["balance_worst_other_residual"] = float(np.max(np.abs(loo[worst][3])))
M["balance_phi"] = PHI_BAL
M["balance_a"]   = A_BAL
BAD_ROW = worst

# A candidate explanation, LABELLED AN INFERENCE and then tested and refused
# as a systematic: 43 is close to the END-of-interval substrate concentration.
alt = {i: A_BAL - Cc1[i]/PHI_BAL for i in range(NOBS)}
lines += ["\n  CANDIDATE EXPLANATION (an inference, not a repair): the printed"
          f" {CAb[worst]:.0f} is within"
          f" {abs(CAb[worst]-end_of_interval)/CAb[worst]*100:.1f} % of that"
          f" interval's END substrate concentration, {end_of_interval:.2f}.",
          "  Tested as a systematic on the neighbouring rows and REFUSED:"]
for i in (worst-1, worst+1):
    lines.append(f"     row {i+1}: end-of-interval {alt[i]:.2f} against the printed"
                 f" {CAb[i]:.0f}  ({abs(alt[i]-CAb[i])/CAb[i]*100:.1f} % off)")
lines.append("  so the table is NOT systematically printing end values, and the"
             " cell is reported, not repaired.")
M["row4_end_of_interval_relgap"] = float(abs(CAb[worst]-end_of_interval)/CAb[worst])
print("\n".join(lines))

# %%
fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.2))
xx = np.linspace(0, Xlm.max()*1.05, 50)
ax[0].plot(xx, A_BAL + B_BAL*xx, color=C_GREY, lw=1.2,
           label=f"line from the other six rows\n$R^2$ = {R2_BAL:.6f}")
keep = [i for i in range(NOBS) if i != BAD_ROW]
ax[0].plot(Xlm[keep], CAb[keep], "o", color=C_BLUE, ms=6, label="printed rows")
ax[0].plot(Xlm[BAD_ROW], CAb[BAD_ROW], "s", color=C_ORANGE, ms=8,
           label=f"row {BAD_ROW+1}: printed {CAb[BAD_ROW]:.0f},\nline says {implied:.1f}")
ax[0].annotate("", xy=(Xlm[BAD_ROW], implied), xytext=(Xlm[BAD_ROW], CAb[BAD_ROW]),
               arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=1.2))
ax[0].set_xlabel(r"$\bar C_C$ over the interval (log mean)")
ax[0].set_ylabel(r"$\bar C_A$ as printed")
ax[0].set_title("the constant-yield balance the table\nmust satisfy, and one cell that does not",
                fontsize=9)
ax[0].legend(fontsize=7.5, loc="upper right")

ax[1].plot(CAb, mu, "o", color=C_BLUE, ms=6)
for i in range(NOBS):
    ax[1].annotate(str(i+1), (CAb[i], mu[i]), textcoords="offset points",
                   xytext=(5, -3), fontsize=7.5, color=C_GREY)
ax[1].plot(implied, mu[BAD_ROW], "s", color=C_ORANGE, ms=8, mfc="none",
           label=f"row {BAD_ROW+1} at the balance-implied $\\bar C_A$")
ax[1].set_xlabel(r"$\bar C_A$ as printed (units not printed)")
ax[1].set_ylabel(r"$\mu_i = \ln(C_C^{end}/C_C^{start})/\Delta t$,  h$^{-1}$")
ax[1].set_title("the seven intervals", fontsize=9)
ax[1].legend(fontsize=7.5, loc="lower right")
for a in ax:
    a.grid(alpha=.25)
fig.tight_layout()
plt.show()
print("Data: Levenspiel (1999) Problem 29.18 p. 644, attributed by him to")
print("Monod (1958) p. 74. Neither Monod publication was consulted.")

# %% [markdown]
# ## PyMRM implementation
#
# Four pieces, in increasing order of how much pymrm does:
#
# 1. **The five growth laws**, each written as $\mu = \mu_m\,g(\theta; S, X)$ so
#    that on route A $\mu_m$ can be profiled out in closed form. That is not
#    cosmetic: Blackman's objective has a kink wherever $2K_s$ crosses a data
#    point, and a single local start lands on the wrong minimum - measured, and
#    kept as a break-table row.
# 2. **Route B's forward prediction.** Constant yield collapses the batch to one
#    ODE, $\mathrm{d}C_C/\mathrm{d}t = \mu(a - C_C/\varphi,\,C_C)\,C_C$, whose
#    solution over each printed interval gives $\mu_i^{\rm pred}$ with no
#    averaging approximation anywhere. Rescaling $\tau = \mu_m t$ removes $\mu_m$
#    from that ODE entirely, so one dense-output solve per $(\theta,\,$interval$)$
#    prices every $\mu_m$ - which is what makes a **global** scan over $\theta$
#    affordable rather than a local start from a guess.
# 3. **A pymrm backward-Euler batch marcher**, `NumJac((1, 2))` + `newton` on the
#    two fields $(S, X)$, used as the **second, independent computation** of
#    route B's fit: different discretisation, different Jacobian, different linear
#    algebra, and it integrates $S$ as an unknown rather than substituting the
#    yield relation - so the yield line is *tested* rather than assumed, through
#    $S + X/\varphi = {\rm const}$.
# 4. **A plug-flow fermenter and a chemostat.** Levenspiel's Fig. 29.1 is captioned
#    *"Behavior of batch **or plug flow** reactors for Monod type microbial
#    fermentation"*, and that equivalence is what makes the plug-flow solve
#    testable: its outlet at residence time $\tau$ must equal the batch at
#    $t = \tau$, to the discretisation error and no further.
#    `construct_convflux_upwind` + `construct_div` (`nu=0`, Cartesian) assemble
#    $v\,\partial c/\partial z$ once; `NumJac((n, 2))` couples the two fields
#    pointwise, which is what a growth law is; `newton` solves. The chemostat is
#    Rawlings' eq. (10.18) with the same source term.
#
# **Boundary conditions use the outward normal.** At the inlet $z=0$ the outward
# normal points in $-z$, so the feed is `{"a": 0, "b": 1, "d": (S_f, X_f)}`; at
# the outlet the equation is pure convection, so the outlet is `{"outflow": True}`:
# the face carries the last cell's value, which is what first-order upwind
# transports. Outlet values come from `compute_boundary_values` under that same
# condition. Orders come from successive differences (`pymrm.checks.observed_orders`),
# each with its round-off companion: the finest difference over newton's floor.

# %%
# ---------------------------------------------------- the five printed laws
# mu = mu_m * g(theta; S, X).  theta is what is NOT linear.
G = {**GROWTH,      # the five laws and the zero-order null, from model.py
     # the first-order null, scaled by the largest printed CA_bar
     "first order": lambda th, S, X: S/np.max(CAb)}
NTH  = {"Monod": 1, "Blackman": 1, "Tessier": 1, "Moser": 2, "Contois": 1,
        "zero order": 0, "first order": 0}
LAWS  = ["Monod", "Blackman", "Tessier", "Moser", "Contois"]
NULLS = ["zero order", "first order"]
# GLOBAL scan grids.  A local start is not safe here - see the break table.
RANGE = {"Monod":    np.geomspace(1e-2, 1e5, 240),
         "Blackman": np.geomspace(1e-2, 1e5, 240),
         "Tessier":  np.geomspace(1e-5, 1e1, 240),
         "Contois":  np.geomspace(1e-4, 1e3, 240)}
MOSER_KS, MOSER_N = np.geomspace(1e-2, 1e12, 60), np.linspace(0.3, 6.0, 36)
MUM_LO, MUM_HI = 0.05, 5.0


def thetas(name, n_grid=None):
    """The global scan: RANGE[name], or the MOSER_KS x MOSER_N (or n_grid) grid."""
    if NTH[name] == 0:
        return [[]]
    if NTH[name] == 1:
        return [[t] for t in RANGE[name]]
    return [[k, n] for n in (MOSER_N if n_grid is None else n_grid) for k in MOSER_KS]


def check_interior(name, t, n_grid=None, what="scan"):
    """Raise when a scan minimum sits on an edge of its named grid: the optimum is
    then at or beyond the end of the grid, and a refinement from there is not a fit."""
    grids = ([RANGE[name]] if NTH[name] == 1 else
             [MOSER_KS, MOSER_N if n_grid is None else n_grid] if NTH[name] == 2 else [])
    for v, g in zip(t, grids):
        if len(g) > 1 and (v <= g[0] or v >= g[-1]):
            raise RuntimeError(f"{what} minimum of {name} on the grid edge"
                               f" ({v:g} in [{g[0]:g}, {g[-1]:g}])")


def Sof(Cc, a, phi):
    return a - Cc/phi


def profile_mum(gvec, y):
    """route A: mu_m enters linearly, so minimise ||mu_m g - y||^2 in closed form"""
    s = float(np.sum(y*gvec)/np.sum(gvec**2))
    return s, float(np.sum((s*gvec - y)**2))


def fit_A(name, y=mu, X=None, S=None, n_grid=None, n_fixed=None):
    """Route A: mu_m profiled out, a global scan on the named grid (raising at its
    edge), then a bounded refinement inside the winning cell. S is the substrate
    column (the printed CA_bar by default); n_grid replaces MOSER_N and n_fixed
    holds Moser's exponent at a value (a one-parameter fit in K)."""
    X = Xlm if X is None else X
    S = CAb if S is None else S
    if n_fixed is not None:                        # Moser with n held: scan K only
        gk = lambda v: G[name]([v, n_fixed], S, X)
        sse = np.array([profile_mum(gk(k), y)[1] for k in MOSER_KS])
        j = int(np.argmin(sse))
        if j in (0, len(MOSER_KS) - 1):
            raise RuntimeError(f"Moser (n = {n_fixed}) scan minimum on the MOSER_KS edge")
        r = minimize_scalar(lambda v: profile_mum(gk(np.exp(v)), y)[1],
                            bounds=(np.log(MOSER_KS[j - 1]), np.log(MOSER_KS[j + 1])),
                            method="bounded", options={"xatol": 1e-13})
        t = [float(np.exp(r.x)), float(n_fixed)]
    else:
        th = thetas(name, n_grid)
        sse = np.array([profile_mum(G[name](t, S, X), y)[1] for t in th])
        t = list(th[int(np.argmin(sse))])
        check_interior(name, t, n_grid)
    if n_fixed is None and len(t) == 1:            # refine inside the winning cell
        lo, hi = t[0]/1.3, t[0]*1.3
        r = minimize_scalar(lambda v: profile_mum(G[name]([v], S, X), y)[1],
                            bounds=(lo, hi), method="bounded",
                            options={"xatol": 1e-13})
        if min(r.x - lo, hi - r.x) < 1e-9*r.x:
            raise RuntimeError(f"{name} refinement ended on its bracket edge")
        t = [float(r.x)]
    elif n_fixed is None and len(t) == 2:
        r = minimize(lambda v: profile_mum(G[name]([np.exp(v[0]), v[1]], S, X), y)[1],
                     [np.log(t[0]), t[1]], method="Nelder-Mead",
                     options={"xatol": 1e-12, "fatol": 1e-16,
                              "maxiter": 20000, "maxfev": 20000})
        if not r.success:
            raise RuntimeError(f"{name} Nelder-Mead refinement: {r.message}")
        t = [float(np.exp(r.x[0])), float(r.x[1])]
    mum, sse = profile_mum(G[name](t, S, X), y)
    pred = mum*G[name](t, S, X)
    return dict(p=[mum] + t, sse=sse, rms=float(np.sqrt(sse/len(y))),
                npar=1 + NTH[name], pred=pred, resid=pred - y)


def psi(name, th, a, phi, tau_max=MUM_HI, rtol=1e-8):
    """mu_m-FREE trajectories: dCc/dtau = g(S(Cc)) Cc with tau = mu_m t.
    One dense-output solve per interval prices every mu_m at once."""
    out = []
    for i in range(NOBS):
        rhs = lambda t, y: [G[name](th, max(Sof(y[0], a, phi), 1e-14), y[0])*y[0]]
        out.append(integrate_interval(rhs, tau_max*dt[i], [Cc0[i]], rtol, rtol*1e-2,
                                      dense_output=True).sol)
    return out


def mu_pred_ode(name, mum, th, a, phi, rtol=1e-12, method="LSODA", dts=None):
    """route B's forward prediction: integrate to the printed dt (or dts), adaptively."""
    dts = dt if dts is None else dts
    out = np.empty(NOBS)
    for i in range(NOBS):
        rhs = lambda t, y: [mum*G[name](th, max(Sof(y[0], a, phi), 1e-14), y[0])*y[0]]
        s = integrate_interval(rhs, dts[i], [Cc0[i]], rtol, rtol*1e-2, method=method)
        out[i] = np.log(s.y[0, -1]/Cc0[i])/dts[i]
    return out


def fit_B(name, y=mu, a=None, phi=None, p0=None, tol=1e-13, refine=True):
    """global scan on route B's OWN objective, then an adaptive-ODE refinement."""
    a = A_BAL if a is None else a
    phi = PHI_BAL if phi is None else phi
    if p0 is None:
        best = (np.inf, None)
        for th in thetas(name):        # the same named grids as route A
            sols = psi(name, th, a, phi)
            r = minimize_scalar(
                lambda m: float(np.sum((np.array(
                    [np.log(float(sols[i](m*dt[i])[0])/Cc0[i])/dt[i]
                     for i in range(NOBS)]) - y)**2)),
                bounds=(MUM_LO, MUM_HI), method="bounded", options={"xatol": 1e-7})
            if r.fun < best[0]:
                best = (float(r.fun), [float(r.x)] + list(th))
        p0 = best[1]
        check_interior(name, p0[1:], what="route-B scan")
        if min(p0[0] - MUM_LO, MUM_HI - p0[0]) < 1e-6:
            raise RuntimeError(f"route-B scan of {name}: mu_m on [MUM_LO, MUM_HI] edge")
    p0 = np.asarray(p0, float)
    if not refine:                  # the scan optimum itself (a break row)
        res = mu_pred_ode(name, p0[0], p0[1:], a, phi) - y
        return dict(p=list(p0), sse=float(np.sum(res**2)), rms=float(np.sqrt(np.mean(res**2))),
                    npar=1 + NTH[name], pred=res + y, resid=res)
    lo = [1e-9]*len(p0)
    hi = [50.0] + ([1e12]*NTH[name] if name != "Moser" else [1e12, 10.0])
    # central differences with a relative step of 1e-5: the objective carries the
    # adaptive integrator's rtol 1e-12 as noise, which the default forward step
    # of 1.5e-8 turns into a Jacobian error near 1e-4 and an exit that moves with
    # the BLAS kernel; at this step the noise enters the Jacobian near 1e-7
    r = least_squares(lambda p: mu_pred_ode(name, p[0], p[1:], a, phi) - y, p0,
                      bounds=(lo, hi), xtol=tol, ftol=tol, gtol=tol,
                      jac="3-point", diff_step=1e-5)
    if r.status <= 0:
        raise RuntimeError(f"route-B least_squares for {name}: {r.message}")
    sse = float(np.sum(r.fun**2))
    return dict(p=list(r.x), sse=sse, rms=float(np.sqrt(sse/len(y))),
                npar=1 + NTH[name], pred=r.fun + y, resid=r.fun)


print("growth laws:", LAWS, "| nulls:", NULLS)
print("global scan:", {k: len(v) for k, v in RANGE.items()},
      "| Moser:", len(MOSER_KS), "x", len(MOSER_N))


# %%
# BatchMarch, PlugFlowFermenter, batch_reference and Chemostat are in model.py.
def mu_pred_pymrm(name, mum, th, a, phi, nstep=200, dts=None, newton_kw=None):
    """route B's prediction through the pymrm backward-Euler marcher."""
    dts = dt if dts is None else dts
    g = lambda S, X: mum*G[name](th, np.maximum(S, 1e-14), X)
    out = np.empty(NOBS)
    for i in range(NOBS):
        c = BatchMarch(g, phi, nstep, newton_kw).march([[Sof(Cc0[i], a, phi), Cc0[i]]], dts[i])
        out[i] = np.log(c[0, 1]/Cc0[i])/dts[i]
    return out


print("BatchMarch        : NumJac((1,2)) + newton, backward Euler on (S, X)")
print("PlugFlowFermenter : construct_convflux_upwind + construct_div(nu=0) + NumJac((n,2)) + newton")
print("Chemostat         : Rawlings eq. (10.18), same source term, marched to steady state")

# %% [markdown]
# ## Results
#
# ### Route A - the comparison the exercise invites
#
# Fit each law's $\mu(\bar C_A)$ to the seven $\mu_i$, with $\mu_m$ profiled out
# analytically and the remaining parameter found on a global scan and then
# refined. Contois needs a biomass too; the exact time average under exponential
# growth is the **logarithmic** mean of $C_C$, which is what is used (the
# arithmetic mean is a declared convention and a break row).

# %%
RA = {k: fit_A(k) for k in LAWS + NULLS}
rows = []
for k in LAWS + NULLS:
    r = RA[k]
    rows.append(dict(law=k, npar=r["npar"], mu_m=r["p"][0],
                     theta=("-" if NTH[k] == 0 else
                            ", ".join(f"{v:.5g}" for v in r["p"][1:])),
                     RMS=r["rms"], SSE=r["sse"]))
TA = pd.DataFrame(rows)
show(TA.style.format({"mu_m": "{:.6f}", "RMS": "{:.6f}", "SSE": "{:.6g}"})
     .hide(axis="index").set_uuid("routeA"))
for k in LAWS + NULLS:
    M[f"rmsA_{k.replace(' ', '_')}"] = RA[k]["rms"]
FA = ((RA["Monod"]["sse"] - RA["Moser"]["sse"])/(RA["Moser"]["sse"]/(NOBS - 3)))
FA_crit = float(f_dist.ppf(0.95, 1, NOBS - 3))
M["routeA_F_moser_vs_monod"] = float(FA)
M["routeA_n_moser"] = float(RA["Moser"]["p"][2])
M["F_crit_95"] = FA_crit
TWO_PAR = ["Monod", "Blackman", "Tessier", "Contois"]
print(f"\nMonod is LAST of the four two-parameter laws on this reduction:"
      f" RMS {RA['Monod']['rms']:.6f} 1/h against Contois {RA['Contois']['rms']:.6f},"
      f"\nTessier {RA['Tessier']['rms']:.6f} and Blackman {RA['Blackman']['rms']:.6f}"
      f" - a factor {RA['Monod']['rms']/RA['Blackman']['rms']:.2f} between the ends.")
print(f"Moser's exponent comes out n = {RA['Moser']['p'][2]:.4f}, and adding it to"
      f" Monod gives\nF = {FA:.4f} against F(1,{NOBS-3}) = {FA_crit:.4f} at 95 % -"
      f" SIGNIFICANT on this reduction.")

# ---- but that significance verdict is CO-CAUSED by the one flagged cell, and
#      this is where a reader meets the claim, so it is priced HERE and not only
#      in meta.yaml.  The imputed value is a labelled INFERENCE, not a repair.
worstA = max(TWO_PAR, key=lambda k: RA[k]["rms"])
CAB_IMPUTED = A_BAL + B_BAL*Xlm[BAD_ROW]
CAb_imp = CAb.copy()
CAb_imp[BAD_ROW] = CAB_IMPUTED                      # the balance-implied value
RA_IMP = {k: fit_A(k, S=CAb_imp) for k in LAWS}
FA_IMP = float((RA_IMP["Monod"]["sse"] - RA_IMP["Moser"]["sse"])
               / (RA_IMP["Moser"]["sse"]/(NOBS - 3)))
worstA_imp = max(TWO_PAR, key=lambda k: RA_IMP[k]["rms"])
print(f"\nIT IS NOT THE REDUCTION ALONE THAT CLEARS THAT TEST. Put row"
      f" {BAD_ROW+1}'s CA_bar at the\nvalue its own material balance implies,"
      f" {CAB_IMPUTED:.2f} instead of the printed {CAb[BAD_ROW]:.0f} - an"
      f"\nINFERENCE, not a repair, and nothing in the CSV is edited - and route A"
      f" gives\nF = {FA_IMP:.4f} against the same {FA_crit:.4f}: NOT significant."
      f" One printed cell,\nwhich the page flags as off the table's own balance"
      f" line, carries that verdict\nover the threshold on its own.")
print(f"The RANKING claim is unaffected: Blackman still beats Monod with the cell"
      f"\nimputed ({RA_IMP['Blackman']['rms']:.6f} against"
      f" {RA_IMP['Monod']['rms']:.6f}, a factor"
      f" {RA_IMP['Monod']['rms']/RA_IMP['Blackman']['rms']:.2f}\nwhere the printed"
      f" cell gives {RA['Monod']['rms']/RA['Blackman']['rms']:.2f}), and the"
      f" route-A-to-route-B reversal below is the\nreduction. What does move is"
      f" which law comes last of the four: {worstA} on the\nprinted cell,"
      f" {worstA_imp} on the imputed one, and they differ by"
      f" {abs(RA_IMP['Contois']['rms']-RA_IMP['Monod']['rms'])/RA_IMP['Monod']['rms']*100:.2f} %.")

print(f"\nEvery law beats both one-parameter nulls: zero order"
      f" {RA['zero order']['rms']:.6f}, first order\n{RA['first order']['rms']:.6f}"
      f" - Monod is {RA['first order']['rms']/RA['Monod']['rms']:.2f}x better than the"
      f" best of them, and Blackman {RA['first order']['rms']/RA['Blackman']['rms']:.2f}x.")

# %% [markdown]
# ### Route B - the same seven growth rates, averaged properly
#
# $\mu_i$ is *exactly* the interval average of $\mu$, so the honest comparison
# integrates the law across the interval instead of evaluating it at $\bar C_A$.
# Constant yield makes that a single ODE, and the two constants it needs -
# $\varphi$ and $a$ - come from the **material-balance line**, which is a property
# of the printed table and of no growth law. **No parameter is added and the models
# are the same five.**
#
# **What each route reads, exactly.** Both fit the same seven $\mu_i$, which are
# the whole of the kinetic information in the table. They differ in where the
# *substrate* comes from: route A reads all seven printed $\bar C_A$ directly,
# while route B never touches that column at all - it takes the substrate
# trajectory from $\varphi$ and $a$, and those come from the **six-row** balance
# line, the one that excludes the flagged row 4. So the datum route A's ranking
# leans on hardest is exactly the one route B's balance line leaves out. That is a
# real difference between the routes and not a nuisance: the break table prices it
# directly, refitting every route-B law - and route B's own $F$ - with the
# **all-seven** balance line instead.

# %%
# NO WALL-CLOCK TIMING IS PRINTED ANYWHERE ON THIS PAGE: elapsed seconds differ
# between two consecutive executions, and a page whose stream output is not
# reproducible cannot be diffed, which is how a real drift would be spotted.  The runtime lives in meta.yaml, where it belongs.
RB = {k: fit_B(k) for k in LAWS + NULLS}
rows = []
for k in LAWS + NULLS:
    r = RB[k]
    rows.append(dict(law=k, npar=r["npar"], mu_m=r["p"][0],
                     theta=("-" if NTH[k] == 0 else
                            ", ".join(f"{v:.5g}" for v in r["p"][1:])),
                     RMS=r["rms"], **{"RMS route A": RA[k]["rms"]},
                     **{"B / A": r["rms"]/RA[k]["rms"]}))
TB = pd.DataFrame(rows)
show(TB.style.format({"mu_m": "{:.6f}", "RMS": "{:.6f}", "RMS route A": "{:.6f}",
                      "B / A": "{:.4f}"}).hide(axis="index").set_uuid("routeB"))
for k in LAWS + NULLS:
    M[f"rmsB_{k.replace(' ', '_')}"] = RB[k]["rms"]
FB = ((RB["Monod"]["sse"] - RB["Moser"]["sse"])/(RB["Moser"]["sse"]/(NOBS - 3)))
M["routeB_F_moser_vs_monod"] = float(FB)
M["routeB_n_moser"] = float(RB["Moser"]["p"][2])
M["routeB_monod_mu_m"] = float(RB["Monod"]["p"][0])
M["routeB_monod_K_S"] = float(RB["Monod"]["p"][1])
M["routeA_over_routeB_monod_rms"] = RA["Monod"]["rms"]/RB["Monod"]["rms"]
worstA = max(TWO_PAR, key=lambda k: RA[k]["rms"])
worstB = max(TWO_PAR, key=lambda k: RB[k]["rms"])
print(f"\nThe ranking changes. Monod's residual falls"
      f" {100*(1 - RB['Monod']['rms']/RA['Monod']['rms']):.0f} %, from"
      f" {RA['Monod']['rms']:.6f} to {RB['Monod']['rms']:.6f} 1/h;"
      f"\nthe worst of the four two-parameter laws goes from {worstA} to {worstB}.")
print(f"Moser's exponent drops from {RA['Moser']['p'][2]:.4f} to"
      f" {RB['Moser']['p'][2]:.4f} and F from {FA:.4f} to {FB:.4f}, against the"
      f"\nsame F(1,{NOBS-3}) = {FA_crit:.4f}: NOT significant.")
print(f"\nMonod, route B, FITTED HERE: mu_m = {RB['Monod']['p'][0]:.6f} 1/h,"
      f" K_S = {RB['Monod']['p'][1]:.4f}\nin the units of the printed CA_bar column,"
      f" which the source does not print.")

# %%
# ---- Contois lands on the Monod curve, to the optimiser's own tolerance
mumC, KsxC = RB["Contois"]["p"]
beta = 1 - KsxC*PHI_BAL
map_mum, map_Ks = mumC/beta, KsxC*PHI_BAL*A_BAL/beta
# THE RULE, stated once here and applied over the WHOLE class in the claims
# cell: a quantity this page PROVES is exactly zero is an optimiser or
# root-finder exit tolerance, not a measurement, so its digits differ between
# machines.  It is registered as a BOUND (docs/claims-design.md), checked on
# every run against a limit that holds on any machine, and never compared as a
# value.  The printed digits below are this run's; the prose quotes the bound.
MAP_BOUND = 1e-3      # relative: the bound on the two map errors and on S_f/a - 1
RMS_GAP_BOUND = 1e-10  # 1/h: the bound on the Contois-minus-Monod RMS gap
MAP_MUM_RELERR = M["contois_map_mum_relerr"] = float(
    abs(map_mum - RB["Monod"]["p"][0])/RB["Monod"]["p"][0])
MAP_KS_RELERR = M["contois_map_KS_relerr"] = float(
    abs(map_Ks - RB["Monod"]["p"][1])/RB["Monod"]["p"][1])
CONTOIS_MINUS_MONOD_RMS = M["contois_minus_monod_rms"] = float(
    abs(RB["Contois"]["rms"] - RB["Monod"]["rms"]))
print("The Contois fit, pushed through the identity proved in *The published model*:")
print(f"   mu_m -> {map_mum:.5f}   against the Monod fit's"
      f" {RB['Monod']['p'][0]:.5f}   (rel. below {MAP_BOUND:g}) (here {MAP_MUM_RELERR:.1e})")
print(f"   K_S  -> {map_Ks:.3f}   against the Monod fit's"
      f" {RB['Monod']['p'][1]:.3f}   (rel. below {MAP_BOUND:g}) (here {MAP_KS_RELERR:.1e})")
print(f"   residual RMS differs by less than {RMS_GAP_BOUND:g} 1/h (here"
      f" {CONTOIS_MINUS_MONOD_RMS:.1e})")
print("   All three of those numbers are BOUNDS, not values: the page proves this")
print("   is one curve fitted twice, so they are two optimisers' exit tolerances,")
print("   and each is checked against a limit that holds on any machine.")
print("\nTwo optimisers, started from different places on two different")
print("parameterisations, landed on the same curve - because it IS the same curve.")

# %% [markdown]
# ### Why the two routes disagree: the reduction bias, measured
#
# $\mu$ is curved, so evaluating it at $\bar C_A$ is not the same as averaging it
# across the interval. The gap is second order in the substrate swing - and the
# swing is not small: **the last two intervals consume more substrate than their
# own printed average.**

# %%
pA = RA["Monod"]["p"]
point_est = pA[0]*G["Monod"](pA[1:], CAb, Xlm)
exact_est = mu_pred_ode("Monod", pA[0], pA[1:], A_BAL, PHI_BAL)
bias = point_est - exact_est
swing = (Cc1 - Cc0)/PHI_BAL                       # substrate consumed in the interval
TBias = pd.DataFrame({"interval": np.arange(1, NOBS+1), "CA_bar printed": CAb,
                      "substrate consumed": swing, "swing / CA_bar": swing/CAb,
                      "mu(CA_bar)": point_est, "<mu> exact": exact_est, "bias": bias})
show(TBias.style.format({"CA_bar printed": "{:.0f}", "substrate consumed": "{:.2f}",
                         "swing / CA_bar": "{:.3f}", "mu(CA_bar)": "{:.5f}",
                         "<mu> exact": "{:.5f}", "bias": "{:+.5f}"})
     .hide(axis="index").set_uuid("bias"))
i_worst = int(np.argmax(np.abs(np.where(np.arange(NOBS) != BAD_ROW, bias, 0.0))))
M["bias_max_abs"] = float(np.max(np.abs(bias)))
M["bias_max_abs_excl_bad_row"] = float(np.abs(bias[i_worst]))
M["bias_over_routeA_rms"] = M["bias_max_abs_excl_bad_row"]/RA["Monod"]["rms"]
M["max_swing_over_CAbar"] = float(np.max(swing/CAb))
print(f"\nLargest bias away from the flagged row: {bias[i_worst]:+.6f} 1/h at"
      f" CA_bar = {CAb[i_worst]:.0f} (interval {i_worst+1}),"
      f"\nwhich is {M['bias_over_routeA_rms']:.2f} times the residual scatter route A"
      f" ranks the five laws on.")
print(f"Row {BAD_ROW+1}'s entry, {bias[BAD_ROW]:+.6f}, is a different animal: it is"
      f" dominated by that row's\ndeparture from the material balance, not by curvature,"
      f" and is reported separately.")
print(f"Largest substrate swing relative to the printed average:"
      f" {M['max_swing_over_CAbar']:.3f} (interval {int(np.argmax(swing/CAb))+1}).")
print("\nThe bias has a SIGN: mu is concave, so mu(mean) exceeds the mean of mu, and")
print("the shortcut therefore makes every saturating law look too gradual near")
print("exhaustion - which is exactly the direction that inflates Moser's exponent.")


# %% [markdown]
# ### What that costs: a seeded bootstrap, with a control
#
# Generate data **from a Monod model**, add Gaussian noise at the observed residual
# scale, and run route A's F test on it. Do it twice, changing only what "a Monod
# model" means:
#
# - **from the forward model** - i.e. from data Monod would actually produce;
# - **from $\mu(\bar C_A)$** - route A's *own* model, for which route A is
#   correctly specified.
#
# The second is the control. If route A's size is inflated in both, the culprit is
# the small-sample F approximation; if only in the first, it is the reduction.
#
# **Every draw is fitted the way route A fits the real data - global scan, then
# refined.** This is not a detail. The observed $F$ the draws are compared against
# comes from a *refined* fit, and reading the draws' fits off the scan grid with no
# refinement is worth a factor two in the control's size: an unrefined Moser fit
# cannot reach its own optimum, so it under-states $F$ on every draw and makes the
# test look conservative when it is not. The refinement is a per-draw zoom on the winning grid cell, its convergence
# is printed below, and the grid-only numbers are kept as a **break row** because
# that is exactly the class of defect - a grid-limited extremum - that no
# perturbation test can see.

# %%
def gmat(name, X=None):
    X = Xlm if X is None else X
    th = thetas(name)
    return np.array([G[name](t, CAb, X) for t in th]), th


GA_M, _ = gmat("Monod")
GA_S, thA_S = gmat("Moser")
nlist = np.array([t[1] for t in thA_S])
# grid edges for the bootstrap: the upper end of n is the step limit
# n -> infinity and the Moser K's upper end the power-law limit K -> infinity, so a
# draw's optimum can lie beyond any finite grid; such draws are counted and printed,
# and each zoom continues from its grid winner as far as it reaches past the edge,
# so such a draw's F is a lower bound
EDGE_S = np.array([t[0] in (MOSER_KS[0], MOSER_KS[-1]) or t[1] in (MOSER_N[0], MOSER_N[-1])
                   for t in thA_S])
DLOG_M = np.log(RANGE["Monod"][1]/RANGE["Monod"][0])     # grid node spacings
DLOG_S, DN_S = np.log(MOSER_KS[1]/MOSER_KS[0]), MOSER_N[1] - MOSER_N[0]
ZOOM = 8                                # refinement rounds; convergence printed


def _sse_M(Y, K):
    """profiled SSE of Monod at a per-draw K (broadcasting over the leading axes)"""
    g = CAb/(K[..., None] + CAb)
    return (np.einsum("...i,...i->...", Y, Y)
            - np.einsum("...i,...i->...", Y, g)**2/np.einsum("...i,...i->...", g, g))


def _sse_S(Y, K, n):
    """profiled SSE of Moser at a per-draw (K, n)"""
    u = CAb**n[..., None]
    g = u/(K[..., None] + u)
    return (np.einsum("...i,...i->...", Y, Y)
            - np.einsum("...i,...i->...", Y, g)**2/np.einsum("...i,...i->...", g, g))


def _zoom1(Y, lk, half, rounds, m=9):
    """per-draw zoom on log K: m nodes across the bracket, halve, repeat"""
    for _ in range(rounds):
        cand = lk[:, None] + np.linspace(-1, 1, m)[None, :]*half[:, None]
        j = np.argmin(_sse_M(Y[:, None, :], np.exp(cand)), axis=1)
        lk, half = cand[np.arange(len(Y)), j], half*2.0/(m - 1)
    return np.exp(lk)


def _zoom2(Y, lk, n, hK, hN, rounds, m=7):
    """the same in two dimensions, on (log K_s, n)"""
    for _ in range(rounds):
        cK = lk[:, None] + np.linspace(-1, 1, m)[None, :]*hK[:, None]
        cN = np.maximum(n[:, None] + np.linspace(-1, 1, m)[None, :]*hN[:, None], 1e-6)
        s = _sse_S(Y[:, None, None, :], np.exp(cK)[:, :, None] + 0*cN[:, None, :],
                   cN[:, None, :] + 0*cK[:, :, None])
        jk, jn = np.unravel_index(np.argmin(s.reshape(len(Y), -1), axis=1), (m, m))
        lk, n = cK[np.arange(len(Y)), jk], cN[np.arange(len(Y)), jn]
        hK, hN = hK*2.0/(m - 1), hN*2.0/(m - 1)
    return np.exp(lk), n


def boot_stats(Y, rounds=None):
    """route-A F statistic and n_hat for a whole block of synthetic datasets.
    Global scan on the SAME grids route A uses, then the SAME refinement:
    rounds=0 reproduces a grid-only fit and is kept as a break row.  Also returns
    a diagnostic dict whose `nclip` is the refinement's convergence evidence,
    counted FROM THE GRID START ALONE - the note under the depth table says why
    it cannot be counted from the shipped pair of starts."""
    rounds = ZOOM if rounds is None else rounds
    nd = len(Y)
    y2 = (Y**2).sum(1)[:, None]
    q1 = y2 - (Y @ GA_M.T)**2/np.einsum("mi,mi->m", GA_M, GA_M)[None, :]
    q2 = y2 - (Y @ GA_S.T)**2/np.einsum("mi,mi->m", GA_S, GA_S)[None, :]
    i1, i2 = np.argmin(q1, axis=1), np.argmin(q2, axis=1)
    ntop = int(np.sum((i1 == 0) | (i1 == len(RANGE["Monod"]) - 1) | EDGE_S[i2]))
    if rounds == 0:
        s1, s2 = q1[np.arange(nd), i1], q2[np.arange(nd), i2]
        # nclip is measured from this row's own grid minima; nclip_pair is a
        # STRUCTURAL 0 - there is no pair of starts at depth 0 - and the sum that
        # reports it below excludes this row for exactly that reason.
        return ((s1 - s2)/(s2/(NOBS - 3)), nlist[i2],
                dict(nclip=int(np.sum(s2 > s1*(1 + 1e-12))), nclip_pair=0, ntop=ntop,
                     ident=0.0, nest_wins=0, best_gain=0.0, flips=0))
    K1 = _zoom1(Y, np.log(RANGE["Monod"][i1]), np.full(nd, DLOG_M), rounds)
    s1 = _sse_M(Y, K1)
    # TWO starts for the two-dimensional fit: the grid winner, and the NESTING
    # POINT n = 1 at the refined Monod K_S.  The second start is a MONOTONICITY
    # GUARANTEE, not a basin escape: at it the Moser objective IS the Monod
    # objective - the same quantity to the last bit, `ident` below - and _zoom2's
    # candidate set contains its own centre, so a fit started there can never end
    # above the Monod fit it is being tested against.  Measured below: it changes
    # no rejection decision on either arm and its largest SSE gain is a few parts
    # in 1e5.  What it buys is that the projection onto the nesting constraint
    # two blocks down is inert - and that is also why the count of draws where
    # that projection binds measures nothing at all.
    ident = float(np.max(np.abs(_sse_S(Y, K1, np.ones(nd)) - s1)))
    sols = []
    for lk0, n0 in ((np.log(np.array([thA_S[j][0] for j in i2])),
                     np.array([thA_S[j][1] for j in i2])),
                    (np.log(K1), np.ones(nd))):
        K2, n2 = _zoom2(Y, lk0, n0, np.full(nd, DLOG_S), np.full(nd, DN_S), rounds)
        sols.append((_sse_S(Y, K2, n2), n2))
    take = sols[0][0] <= sols[1][0]
    s2 = np.where(take, sols[0][0], sols[1][0])
    nhat = np.where(take, sols[0][1], sols[1][1])
    # Moser CONTAINS Monod, so its SSE cannot exceed Monod's, and a fit that ends
    # above it is not converged.  That count is informative for the GRID START
    # (nclip) and is identically zero for the shipped pair (nclip_pair), for the
    # reason given above; both are returned and both are printed.
    nclip = int(np.sum(sols[0][0] > s1*(1 + 1e-12)))
    nclip_pair = int(np.sum(s2 > s1*(1 + 1e-12)))
    gain = float(np.max((sols[0][0] - s2)/s2))
    sg = np.minimum(sols[0][0], s1)                 # the grid start on its own
    F_grid = (s1 - sg)/(sg/(NOBS - 3))
    s2 = np.minimum(s2, s1)
    F = (s1 - s2)/(s2/(NOBS - 3))
    return F, nhat, dict(nclip=nclip, nclip_pair=nclip_pair, ntop=ntop, ident=ident,
                         nest_wins=int(np.sum(sols[1][0] < sols[0][0])),
                         best_gain=gain,
                         flips=int(np.sum((F > FA_crit) != (F_grid > FA_crit))))


def prop_se(p, n=None):
    """binomial standard error of a bootstrap frequency, at n (default NBOOT) draws"""
    return float(np.sqrt(p*(1.0 - p)/(NBOOT if n is None else n)))


def median_se(v):
    """distribution-free standard error of a median, from the order statistics:
    half the width of the interval the binomial sign count puts one standard
    deviation either side of the middle rank.  No RNG, so it is deterministic."""
    w, N = np.sort(v), len(v)
    k = 0.5*np.sqrt(N)
    return float(0.5*(w[min(N - 1, int(np.ceil(N/2 + k)))]
                      - w[max(0, int(np.floor(N/2 - k)))]))


def boot_arm(yt, sd, seed=None, rounds=None, nboot=None):
    nb_ = NBOOT if nboot is None else nboot
    Y = yt[None, :] + np.random.default_rng(SEED if seed is None else seed
                                            ).normal(0.0, sd, (nb_, NOBS))
    Fb, nb, dg = boot_stats(Y, rounds)
    size = float((Fb > FA_crit).mean())
    pF, pn = float((Fb > FA).mean()), float((nb > RA["Moser"]["p"][2]).mean())
    # EVERY number below is a Monte-Carlo estimate at NBOOT draws - `size`, `pF`
    # and `pn` are frequencies, `med` is a median - so every one of them carries
    # a sampling error, printed beside it everywhere it is quoted: prop_se is the
    # binomial one and median_se the order-statistic one
    return dict(size=size, med=float(np.median(nb)), pF=pF, pn=pn,
                size_se=prop_se(size, nb_), pF_se=prop_se(pF, nb_), pn_se=prop_se(pn, nb_),
                med_se=median_se(nb), **dg)


sigma, NBOOT = RB["Monod"]["rms"], 2000
y_true = mu_pred_ode("Monod", RB["Monod"]["p"][0], RB["Monod"]["p"][1:],
                     A_BAL, PHI_BAL)
y_ctrl = point_est                                  # route A's own model
ARMS = (("forward", y_true), ("control", y_ctrl))

# ---- refinement evidence FIRST: the six numbers against the zoom depth
print(f"{NBOOT} seeded draws each, sigma = {sigma:.6f} 1/h, seed {SEED},"
      f" nominal size 0.0500")
print("\nrefinement of each draw's fit (rounds = 0 is the grid alone; the grid is"
      f"\n{len(RANGE['Monod'])} Monod nodes and {len(MOSER_KS)} x {len(MOSER_N)}"
      f" Moser nodes, an n spacing of {DN_S:.3f}):\n")
print(f"{'rounds':>7}{'fwd size':>10}{'fwd med n':>11}{'fwd p':>8}"
      f"{'ctrl size':>11}{'ctrl med n':>12}{'unconverged':>13}")
CONV = {}
for r in (0, 2, 4, 6, ZOOM):
    CONV[r] = {tag: boot_arm(yt, sigma, rounds=r) for tag, yt in ARMS}
    f_, c_ = CONV[r]["forward"], CONV[r]["control"]
    print(f"{r:>7}{f_['size']:>10.4f}{f_['med']:>11.4f}{f_['pF']:>8.4f}"
          f"{c_['size']:>11.4f}{c_['med']:>12.4f}"
          f"{f_['nclip'] + c_['nclip']:>13d}"
          + ("   <- the grid alone" if r == 0 else
             "   <- converged, used below" if r == ZOOM else ""))
# summed over the REFINED rows only: the rounds = 0 row is fitted from the grid
# winner alone, so it has no pair of starts and its entry is a structural 0 the
# code writes rather than a count it takes - folding that in would have made the
# claim below rest partly on a hard-coded zero.
N_REFINED = sum(1 for r in CONV if r > 0)
N_PAIR_CLIP = sum(CONV[r][t]["nclip_pair"]
                  for r in CONV if r > 0 for t in ("forward", "control"))
IDENT_MAX = max(CONV[r][t]["ident"] for r in CONV for t in ("forward", "control"))
print(f"\n  'unconverged' counts the draws, out of {2*NBOOT}, on which the"
      f" two-dimensional fit\n  STARTED FROM THE GRID WINNER ALONE still ends ABOVE"
      f" the one-dimensional one -\n  impossible for a converged fit, since Moser"
      f" contains Monod. Read it as the\n  convergence evidence for the zoom, and"
      f" read the six numbers' own trend beside\n  it: both say the same thing about"
      f" where this table settles.")
print(f"\n  WHAT IT CANNOT SHOW. The same count for the fit actually shipped - the"
      f" pair of\n  starts - is {N_PAIR_CLIP} summed over the {N_REFINED} REFINED rows"
      f" above (the grid-only row is\n  fitted from the grid winner alone, so it has no"
      f" pair of starts and no such\n  count), and it cannot be anything else at any"
      f" refined depth: the second start IS\n  the nesting point, where the Moser"
      f" objective and the Monod objective are the\n  same number to the last bit"
      f" (largest difference over all {2*NBOOT} draws and all\n  depths:"
      f" {IDENT_MAX:.1e}), and the zoom's candidate set contains its"
      f"\n  own centre, so that branch is monotone and can never end above it. A"
      f" column of\n  zeros there would be an identity of the code, not a"
      f" measurement, so it is not\n  what is printed. The shipped fit is at least"
      f" as converged as the column says.")
out = CONV[ZOOM]
GRID_ONLY = CONV[0]
M["boot_size_forward_truth"] = out["forward"]["size"]
M["boot_size_control_truth"] = out["control"]["size"]
M["boot_median_n_forward_truth"] = out["forward"]["med"]
M["boot_median_n_control_truth"] = out["control"]["med"]
M["boot_p_of_observed_F"] = out["forward"]["pF"]
M["boot_p_of_observed_n"] = out["forward"]["pn"]
print(f"draws whose grid winner is on a search-grid edge: {out['forward']['ntop']} forward,"
      f" {out['control']['ntop']} control, of {NBOOT} each.\nThe upper ends of the grids are"
      f" the limits n -> infinity (a step) and K -> infinity, which no finite grid"
      f"\ncontains. Each is refined from its grid winner as far as the zoom reaches past the"
      f" edge,\nso its F is a lower bound, and INCLUDED in every bootstrap statistic above;"
      f" the count is a claim.")

pm = lambda v, s: f"{v:.4f} +- {s:.4f}"
print(f"\n{'truth':<46}{'size':>18}{'median n':>18}{'p(F>obs)':>18}")
for tag, lab in (("forward", "the forward Monod model (what Monod produces)"),
                 ("control", "mu(CA_bar) - route A's own model [CONTROL]")):
    o = out[tag]
    print(f"{lab:<46}{pm(o['size'], o['size_se']):>18}"
          f"{pm(o['med'], o['med_se']):>18}{pm(o['pF'], o['pF_se']):>18}")
print(f"\nand on the forward arm, p for the OBSERVED exponent:"
      f" {pm(out['forward']['pn'], out['forward']['pn_se'])}")
print(f"\nEVERY NUMBER IN THAT TABLE AND IN THE LINE UNDER IT IS A MONTE-CARLO"
      f" ESTIMATE over\n{NBOOT} seeded draws, and they are NOT all of one kind: the"
      f" sizes and the p values\nare FREQUENCIES, the two exponents are MEDIANS. So"
      f" every one carries a sampling\nerror, printed beside it as +- one standard"
      f" error - binomial for the frequencies,\nand from the order statistics for the"
      f" medians.")


def sig_stable(v, s):
    """how many leading significant figures of v its own standard error s
    supports: the largest k with s below half a unit in the k-th figure."""
    k = 0
    while s < 0.5*10.0**(int(np.floor(np.log10(abs(v)))) - k) and k < 6:
        k += 1
    return k


SIX = [("forward size", out["forward"]["size"], out["forward"]["size_se"]),
       ("control size", out["control"]["size"], out["control"]["size_se"]),
       ("p(F>obs)", out["forward"]["pF"], out["forward"]["pF_se"]),
       ("p(n>obs)", out["forward"]["pn"], out["forward"]["pn_se"]),
       ("forward median n", out["forward"]["med"], out["forward"]["med_se"]),
       ("control median n", out["control"]["med"], out["control"]["med_se"])]
K_FIG = {lab: sig_stable(v, s) for lab, v, s in SIX}
FREQ6 = ("forward size", "control size", "p(F>obs)", "p(n>obs)")
MED6 = ("forward median n", "control median n")
print(f"\nAND HOW MANY FIGURES THAT BUYS, PER STATISTIC rather than in one sentence"
      f" - for\nthe six quoted off this page, the largest k whose k-th significant"
      f" figure one\nstandard error does not reach:\n")
for lab, v, s in SIX:
    half = 0.5*10.0**(int(np.floor(np.log10(abs(v)))) - K_FIG[lab])
    print(f"   {lab:<17}{v:>8.4f} +- {s:<7.4f} one SE is {s/half:>5.1f}x half a unit"
          f" in figure {K_FIG[lab] + 1}"
          f" ({half:.4g}) -> {K_FIG[lab]} figure" + ("s" if K_FIG[lab] != 1 else ""))
# what meta.yaml, README.md and models_entry.yaml say about these six - "the
# frequencies to their leading figure only, the medians to at least two" - is
# this table, and it is ASSERTED here rather than trusted to stay true
if max(K_FIG[lab] for lab in FREQ6) > 1 or min(K_FIG[lab] for lab in MED6) < 2:
    raise RuntimeError(f"supported figures changed: {K_FIG}")
print(f"\nSo the FREQUENCIES are seed-independent in their LEADING FIGURE ONLY"
      f" ({max(K_FIG[lab] for lab in FREQ6)} figure), and\nthe MEDIANS in at least"
      f" {min(K_FIG[lab] for lab in MED6)}. The four-figure values above are exact for"
      f" THIS seed and\nreproducible, not resolved - read the +-, not the digit count."
      f" Every quotation of\nthem on this page carries the same +-.")
print(f"\nThe control sits AT its nominal size ({out['control']['size']:.4f} +-"
      f" {out['control']['size_se']:.4f} against 0.0500) and\nrecovers the truth's"
      f" exponent ({out['control']['med']:.4f} +- {out['control']['med_se']:.4f}"
      f" against 1). So the small-sample F\napproximation at n = 7 is essentially EXACT"
      f" here - it is not conservative, and it\nis not the problem. Against data Monod"
      f" actually produces, the same test rejects\nit {out['forward']['size']:.1%} +-"
      f" {out['forward']['size_se']:.1%} of the time and puts the median exponent at"
      f" {out['forward']['med']:.4f} +- {out['forward']['med_se']:.4f}.\nThe whole of"
      f" the distortion is the reduction.")
print(f"\nWHAT THE SECOND START DOES, measured on the converged row rather than"
      f" asserted:\nit strictly wins on"
      f" {out['forward']['nest_wins'] + out['control']['nest_wins']} of the"
      f" {2*NBOOT} draws, but the largest SSE it ever gains is"
      f"\n{max(out['forward']['best_gain'], out['control']['best_gain']):.1e} relative,"
      f" and it changes"
      f" {out['forward']['flips'] + out['control']['flips']} rejection decisions on the"
      f" two arms combined.\nIt is a monotonicity guarantee - the fit can never end"
      f" above the Monod fit it is\nbeing tested against - and not an escape from a"
      f" wrong basin.")
print(f"\nRead the first row of the table above the same way: fitted on the grid"
      f" alone,\nthe control reads {GRID_ONLY['control']['size']:.4f} and its median"
      f" exponent {GRID_ONLY['control']['med']:.4f} - a coarse Moser grid\n(n spacing"
      f" {DN_S:.3f}) cannot reach its own optimum, so it deflates F on every draw and"
      f"\nmanufactures a conservative-looking control. That is a grid-limited"
      f" extremum,\nand no perturbation of the model would have revealed it.")
print(f"\nThe OBSERVED route-A verdict, F = {FA:.4f}, therefore has p ="
      f" {out['forward']['pF']:.4f} +- {out['forward']['pF_se']:.4f}\nagainst a"
      f" properly generated Monod null - not significant at all, and not close"
      f"\nenough to 0.05 for the sampling error to matter - and its exponent"
      f"\nn = {RA['Moser']['p'][2]:.4f} has p ="
      f" {out['forward']['pn']:.4f} +- {out['forward']['pn_se']:.4f}.")

# ---- the noise scale is a DECLARED CONVENTION for the forward arm, and is not
#      one for the control.  Both are printed rather than asserted.
SIG_ALT = {"route B RMS (used above)": sigma, "route A RMS": RA["Monod"]["rms"]}
SIG_TABLE = {lab: (boot_arm(y_true, sd)["size"], boot_arm(y_ctrl, sd)["size"])
             for lab, sd in SIG_ALT.items()}
print(f"\nWhich residual scatter supplies sigma is a declared convention, and it"
      f" matters\nfor one arm only:\n")
print(f"{'sigma':<26}{'value, 1/h':>12}{'forward size':>18}{'control size':>18}")
for lab, sd in SIG_ALT.items():
    print(f"{lab:<26}{sd:>12.6f}{pm(SIG_TABLE[lab][0], prop_se(SIG_TABLE[lab][0])):>18}"
          f"{pm(SIG_TABLE[lab][1], prop_se(SIG_TABLE[lab][1])):>18}")
print(f"\nThe FORWARD size is a bias-to-noise ratio, so it moves with sigma and the"
      f"\nquoted {out['forward']['size']:.1%} +- {out['forward']['size_se']:.1%} is"
      f" specific to route B's scatter. The"
      f" CONTROL's size is scale-free to\nMonte-Carlo error, which is exactly what a"
      f" correctly specified test should do,\nand is a second reason to read it as"
      f" the F approximation being exonerated.")

# %% [markdown]
# ### What this design *could* separate, if the data were clean
#
# Identifiability is a property of the experiment, not of the model. Take Monod as
# true, at the route-B fit; ask how close each alternative can get to it **at these
# seven design points**; and compare that distance with the scatter. $N^*$ is the
# number of independent repeats of the seven-interval run - Levenspiel's Problem
# 29.18, attributed by him to Monod (1958) - at which a $\chi^2$ separation at 95 %
# becomes possible.
#
# **Read $N^*$ as a threshold, not as a sample size.** It is the $N$ at which the
# *expected* excess sum of squares first **reaches** the 95 % critical value. An
# experiment run at exactly $N^*$ therefore lands above the critical value about
# half the time: $N^*$ is roughly a **50 %-power** point, not a 95 %-detection
# design. Three further choices in it - the confidence level, the degrees of
# freedom, and which route supplies the scatter - are declared conventions with
# break rows, and the separations $d$ that go into it are optimiser *upper* bounds,
# so $N^*$ errs low in the direction that matters.

# %%
s_obs = RB["Monod"]["rms"]
SEP_BOUND = 1e-9          # 1/h: the bound the two degenerate separations are held to
SEP_SCATTER_BOUND = 1e-7  # the same, over the route-B scatter
SEP = {}
sep_rows = []
for k in ["Blackman", "Tessier", "Moser", "Contois"]:
    r = fit_B(k, y=y_true, p0=RB[k]["p"])       # warm start: the manifold is the same
    d = r["rms"]
    SEP[k] = d
    Nstar = (chi2.ppf(0.95, NOBS - 2)*s_obs**2/(NOBS*d**2)) if d > 1e-9 else np.inf
    sep_rows.append(dict(alternative=k, d=d, **{"d / scatter": d/s_obs},
                         **{"N* repeats": Nstar},
                         verdict=("degenerate - see the identities" if d < 1e-9 else
                                  "within reach of one clean run, at ~50 % power"
                                  if Nstar <= 1.5 else "needs more data")))
    # every separation is a claim; the degenerate two are BOUNDS (the rule stated
    # where the Contois map is computed), registered as such in the claims cell
    M[f"sep_{k.lower()}"] = float(d)
    M[f"sep_over_scatter_{k.lower()}"] = float(d/s_obs)
    if np.isfinite(Nstar):     # only where a separation exists at all: for Moser
        M[f"Nstar_{k.lower()}"] = float(Nstar)   # and Contois it is not a number

TS = pd.DataFrame(sep_rows)
show(TS.style.format({"d": "{:.4e}", "d / scatter": "{:.5f}", "N* repeats": "{:.4g}"})
     .hide(axis="index").set_uuid("sep"))
print("\nMoser and Contois sit at the optimiser's own floor, and that is not an")
print("accident of this dataset - it is the algebra of *The published model*: Moser")
print("CONTAINS Monod at n = 1, and Contois IS Monod once a constant yield makes X an")
print("affine function of S. Their symbolic residuals are exactly zero; the numbers")
print(f"above (here {SEP['Moser']:.1e} and {SEP['Contois']:.1e}) measure"
      f" least_squares' xtol, not a\nseparation, so both are registered as BOUNDS,"
      f" below {SEP_BOUND:g} 1/h, and never compared as values.")
print(f"\nTessier and Blackman are genuinely separable. Blackman needs"
      f" {M['Nstar_blackman']:.2f} repeats of this")
print(f"design and Tessier {M['Nstar_tessier']:.1f}. So one clean run of it is"
      f" within {100*(M['Nstar_blackman']-1):.0f} % of enough for")
print(f"Blackman - not enough, and at N* the power is only about a half anyway,")
print(f"so 'one run would do' would be generous twice over.")
print("\nNote which way this cuts against route A: Blackman, the law route A liked")
print("best, is the ONLY one this design can actually resolve - so route A's ranking")
print("is not merely biased, it is confidently ordering laws it cannot tell apart.")

# %% [markdown]
# ### The experiment that *does* separate Contois from Monod
#
# The degeneracy above is a property of a **constant-yield batch run**, which ties
# $X$ to $S$. A chemostat does not: the dilution rate moves them independently.
# Rawlings' eq. (10.18) with a Monod law washes out at
# $D_c = \mu_m S_f/(K_s + S_f)$, which **depends on the feed concentration**. With
# a Contois law the same balance gives $\mu \to \mu_m$ as $X \to 0$, so
# $D_c = \mu_m$ **whatever the feed**. The two cross - and the crossing point is
# not arbitrary.

# %%
mumM, KsM = RB["Monod"]["p"]
Dc_monod = lambda Sf: mumM*Sf/(KsM + Sf)
Dc_contois = mumC                                     # independent of S_f, exactly
ratio = lambda Sf: Dc_contois/Dc_monod(Sf)
# ratio(S_f) falls monotonically in S_f; each level is the one root on a log scan
# of the feeds 1 to 1e5 (CA_bar units), raising unless exactly one is found
SF_LO, SF_HI, SF_XTOL = 1.0, 1e5, 1e-10
Sf_cross = one_root(lambda Sf: ratio(Sf) - 1.0, SF_LO, SF_HI, log=True, xtol=SF_XTOL,
                    what="D_c crossing")
Sf_10pct = one_root(lambda Sf: ratio(Sf) - 1.10, SF_LO, SF_HI, log=True, xtol=SF_XTOL,
                    what="10 % gap")
M["Sf_crossing"] = float(Sf_cross)
M["Sf_for_10pct_gap"] = float(Sf_10pct)
M["Dc_ratio_at_Sf20"] = float(ratio(20.0))
SF_CLOSED = mumC*KsM/(mumM - mumC)        # the same crossing, in closed form
SF_ROOTFIND_RELERR = abs(Sf_cross - SF_CLOSED)/SF_CLOSED
# the root finder's OWN guarantee (brentq inside find_roots), from the declared
# tolerances rather than
# from the residual that happened to come out: scipy returns an x0 with
# |x0 - x*| <= xtol + rtol|x0|, which is 5.1e-13 relative here.  The observed
# residual is three orders inside that, but it is a property of Brent's last
# interpolation, not of anything declared - it moves 40x under a change of
# BRACKET that changes nothing about the mathematics - so the claim
# Sf_rootfind_over_guarantee tests it against this bound, not against the value below.
SF_RTOL = 8.9e-16                         # brentq's default rtol
SF_ROOTFIND_BOUND = (SF_XTOL + SF_RTOL*abs(Sf_cross))/SF_CLOSED
CROSS_VS_A = M["Sf_crossing_vs_a_relerr"] = float(abs(Sf_cross - A_BAL)/A_BAL)
# The two exit tolerances of the Contois fit, SIGNED, and the sensitivity of the
# closed form to each.  S_f/a = (1+e_mu)/[(1+e_K)(1 - e_mu beta/(1-beta))], so to
# first order S_f/a - 1 = e_mu/(1-beta) - e_K: the mu_m error is amplified by
# 1/(1-beta) and the K_S error enters with coefficient EXACTLY 1.
E_MUM = (map_mum - mumM)/mumM
E_KS = (map_Ks - KsM)/KsM
AMP = 1/(1 - beta)
CROSS_DECOMP = E_MUM*AMP - E_KS
print("The crossing is EXACT, and it is the batch run's own substrate intercept:")
print("  putting mu_m^C = beta mu_m^M and K_S = K_sx phi a / beta into")
print("  mu_m^M S_f/(K_S + S_f) = mu_m^C gives S_f = a identically. Numerically:")
print(f"    root-found crossing        S_f = {Sf_cross:.2f}")
print(f"    material-balance intercept a   = {A_BAL:.2f}"
      f"   (rel. gap below {MAP_BOUND:g}) (here {CROSS_VS_A:.2e})")
print(f"\n  That gap is NOT the root-find's. Solving the same equation in closed form,"
      f"\n  S_f = mu_m^C K_S/(mu_m^M - mu_m^C) = {SF_CLOSED:.2f}, and the root find reproduces it"
      f" to\n  round-off (here {SF_ROOTFIND_RELERR:.1e}), inside its own"
      f" declared\n  guarantee (xtol + rtol|x| = {SF_ROOTFIND_BOUND:.1e} relative here),"
      f" which is what is asserted.")
print(f"\n  What the gap (here {CROSS_VS_A:.2e}) measures is the Contois fit's own two exit"
      f" tolerances, and\n  NOT in equal measure. Differentiating that closed form:"
      f" d ln S_f / d ln mu_m^C =\n  1/(1-beta) = {AMP:.3f}, while d ln S_f / d ln K_S"
      f" = 1 EXACTLY. So it is the mu_m error\n  that is amplified and the K_S error"
      f" that is not:")
print(f"      mu_m error (here {E_MUM:+.3e}), amplified x {AMP:.3f}, contributes"
      f" (here {E_MUM*AMP:+.3e})")
print(f"      K_S  error (here {E_KS:+.3e}), unamplified, other sign, contributes"
      f" (here {-E_KS:+.3e})")
print(f"      sum (here {CROSS_DECOMP:+.3e}), which is the gap above to within 0.5 %"
      f" (checked).")
# the decomposition is the justification for treating this quantity as a
# tolerance bound, so it is checked rather than left as prose: the first-order
# form is good to 1e-5
if abs(abs(CROSS_DECOMP) - CROSS_VS_A) > 5e-3*CROSS_VS_A:
    raise RuntimeError(f"first-order decomposition {CROSS_DECOMP} misses {CROSS_VS_A}")
print(f"\n  The two errors are exit tolerances, so their sizes and signs change from"
      f" one\n  machine to the next. The gap is the same exactly-zero\n  quantity as"
      f" the two map errors, and it is a BOUND for the same reason: below"
      f" {MAP_BOUND:g}\n  relative, checked on every run.")
print(f"\n  Away from it they separate: at S_f = 20 the two critical dilution rates"
      f"\n  differ by a factor {M['Dc_ratio_at_Sf20']:.4f}, and the gap exceeds 10 %"
      f" for every S_f\n  below {Sf_10pct:.2f} - in the units of the printed CA_bar"
      f" column, whose scale is set\n  by the run itself (its substrate ran from"
      f" {A_BAL - Cc0[0]/PHI_BAL:.1f} to {A_BAL - Cc1[-1]/PHI_BAL:.2f}).")

Sf_test, y_yield = 20.0, PHI_BAL
WASH, D_LO, D_HI = 1e-6, 0.05, 2.0    # washout level (of y S_f) and the D range scanned


def washout(g_, y_, Sf_, wash=None):
    """Critical dilution rate of the marched chemostat: the one D in [D_LO, D_HI]
    where the steady biomass falls to wash y S_f (a scan of 16 points, then brentq)."""
    w = WASH if wash is None else wash
    return one_root(lambda D: Chemostat(g_, y_, Sf_, D).steady()[1] - w*y_*Sf_,
                    D_LO, D_HI, n_scan=16, xtol=1e-9, what="washout")


gM_ = lambda S, X: mumM*S/(KsM + S)
gC_ = lambda S, X: mumC*S/(KsxC*X + S)
Dc_num = {}
for tag, g_ in (("Monod", gM_), ("Contois", gC_)):
    Dc_num[tag] = washout(g_, y_yield, Sf_test)
M["Dc_monod_closed_form"] = float(Dc_monod(Sf_test))
M["Dc_monod_transient_relerr"] = abs(Dc_num["Monod"] - Dc_monod(Sf_test))/Dc_monod(Sf_test)
M["Dc_contois_transient_relerr"] = abs(Dc_num["Contois"] - Dc_contois)/Dc_contois
print(f"\n  washout root-found on a MARCHED chemostat, S_f = {Sf_test:.0f}:")
print(f"    Monod   D_c = {Dc_num['Monod']:.5f} against the closed form"
      f" {Dc_monod(Sf_test):.5f}  (rel. {M['Dc_monod_transient_relerr']:.2e})")
print(f"    Contois D_c = {Dc_num['Contois']:.5f} against mu_m ="
      f" {Dc_contois:.5f}  (rel. {M['Dc_contois_transient_relerr']:.2e})")
print(f"\n  Two models that no batch run can tell apart differ by a factor"
      f" {M['Dc_ratio_at_Sf20']:.2f} in the\n  dilution rate at which the culture is"
      f" lost. That is the experiment to run.")

# %%
fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.3))
Sg = np.linspace(0.5, CAb.max()*1.08, 400)
Xg = PHI_BAL*(A_BAL - Sg)                       # the run's own yield line
for k in LAWS:
    p = RA[k]["p"]
    ax[0].plot(Sg, p[0]*G[k](p[1:], Sg, np.interp(Sg, CAb[::-1], Xlm[::-1])),
               color=COL[k], lw=1.3, label=f"{k} ({RA[k]['rms']:.4f})")
    p = RB[k]["p"]
    ax[1].plot(Sg, p[0]*G[k](p[1:], Sg, Xg), color=COL[k], lw=1.3,
               label=f"{k} ({RB[k]['rms']:.4f})")
for a in ax[:2]:
    a.plot(CAb, mu, "ko", ms=5, zorder=5)
    a.set_xlabel(r"substrate concentration (units not printed)")
    a.set_ylabel(r"$\mu$,  h$^{-1}$")
    a.legend(fontsize=7, loc="lower right")
    a.grid(alpha=.25)
    a.set_ylim(0, 1.02)
ax[0].set_title("route A: $\\mu$ evaluated at the interval\nAVERAGE (RMS, h$^{-1}$, in brackets)",
                fontsize=9)
ax[1].set_title("route B: the interval average done\nexactly - same data, same models",
                fontsize=9)
Sfg = np.geomspace(2, 1000, 300)
ax[2].plot(Sfg, Dc_monod(Sfg), color=C_BLUE, lw=1.4,
           label="Monod, $\\mu_m S_f/(K_S+S_f)$")
ax[2].axhline(Dc_contois, color=C_YELLOW, lw=1.4, label="Contois, $\\mu_m$ (no $S_f$)")
ax[2].plot([Sf_cross], [Dc_contois], "k*", ms=11, zorder=5,
           label=f"cross at $S_f = a$ = {Sf_cross:.1f}")
ax[2].axvline(Sf_10pct, color=C_GREY, ls=":", lw=1,
              label=f"> 10 % apart below $S_f$ = {Sf_10pct:.0f}")
ax[2].set_xscale("log")
ax[2].set_xlabel("$S_f$ (units of the printed $\\bar C_A$)")
ax[2].set_ylabel("critical dilution rate,  h$^{-1}$")
ax[2].set_title("the batch degeneracy broken:\nchemostat washout", fontsize=9)
ax[2].legend(fontsize=7, loc="lower right")
ax[2].grid(alpha=.25)
fig.tight_layout()
plt.show()
print("Left and centre: curves fitted on this page to Levenspiel p. 644's seven")
print("intervals. Right: those curves extrapolated to an experiment nobody ran.")

# %% [markdown]
# ## Validation
#
# Four kinds, kept apart, and none of them is the fit.
#
# 1. **A stated numerical result reproduced.** Rawlings & Ekerdt print a closed
#    form for the chemostat steady states, eq. (10.19), and the parameter set of
#    their Fig. 10.12 ($\mu_m = 1$, $K_s = 1$, $S_f = 5$, $y = 1$). Both are
#    reproduced here by marching eq. (10.18) to steady state - a different
#    computation entirely - and a break row changes one of their printed constants
#    and watches the reproduction fail.
# 2. **Two independent computations of route B's fit, and a statement of what that
#    does and does not establish.** The quoted constants come from an adaptive
#    LSODA solve; they are re-derived through a **pymrm backward-Euler marcher**
#    (`NumJac((1,2))` + `newton`) with Richardson extrapolation, which shares no
#    line of code, uses a different Jacobian and a different linear solve, and
#    integrates $S$ as an unknown rather than substituting the yield relation.
#    What that catches is a bug in *either integrator*. What it cannot catch is
#    anything the two share: both read $a$, $\varphi$, $\Delta t$ and $C_C^{\rm
#    start}$ from the same place and integrate the same ODE, so a wrong yield, a
#    wrong intercept or a mis-set interval moves **both** predictions together and
#    leaves their agreement where it was. That is measured below rather than
#    asserted.
#    $S + X/\varphi = {\rm const}$ belongs on the same list: `BatchMarch` builds
#    its source as $(-m/\varphi,\,m)$ from one $m$, so backward Euler conserves it
#    identically for *any* growth law. The measurement below shows it moving by
#    nothing under physics that is grossly wrong; only a different yield in the
#    source from the one in the invariant moves it, and that is its break row.
# 3. **An exact identity the discretisation must satisfy.** The plug-flow outlet at
#    residence time $\tau$ is the batch state at $t = \tau$ - Levenspiel's
#    Fig. 29.1 caption, and a theorem. Grid refined, order reported, Richardson
#    checked.
# 4. **A defect-injection table**, built for this page's physics.

# %%
# ---- 1. Rawlings eq. (10.19) and its printed parameter set, reproduced
mum_r, Ks_r, Sf_r, y_r = 1.0, 1.0, 5.0, 1.0
g_r = lambda S, X: mum_r*S/(Ks_r + S)
Dc_r_closed = mum_r*Sf_r/(Ks_r + Sf_r)
rows = []
for D in (0.25, 0.5, 0.75):
    Ss_c = D*Ks_r/(mum_r - D)                       # eq. (10.19)
    Xs_c = y_r*(Sf_r - Ss_c)
    Ss_n, Xs_n = Chemostat(g_r, y_r, Sf_r, D).steady()
    rows.append(dict(D=D, S_closed=Ss_c, S_marched=Ss_n, X_closed=Xs_c, X_marched=Xs_n,
                     **{"max rel. err": max(abs(Ss_n-Ss_c)/Ss_c, abs(Xs_n-Xs_c)/Xs_c)}))
TR = pd.DataFrame(rows)
show(TR.style.format({"D": "{:.2f}", "S_closed": "{:.9f}", "S_marched": "{:.9f}",
                      "X_closed": "{:.9f}", "X_marched": "{:.9f}",
                      "max rel. err": "{:.3e}"}).hide(axis="index").set_uuid("rawl"))
M["rawlings_10_19_max_relerr"] = float(TR["max rel. err"].max())
Dc_r_num = washout(g_r, y_r, Sf_r)
M["rawlings_Dc_relerr"] = abs(Dc_r_num - Dc_r_closed)/Dc_r_closed
print(f"\nRawlings' printed D_c = mu_m S_f/(K_s+S_f) = {Dc_r_closed:.10f}; root-found"
      f" on the marched\nmodel it is {Dc_r_num:.10f}  (rel. {M['rawlings_Dc_relerr']:.2e})."
      f"  Their eq. (10.19) is reproduced\nto below 1e-9 relative (here"
      f" {M['rawlings_10_19_max_relerr']:.1e}) at three dilution rates.")

# %%
# ---- 2. the SECOND, INDEPENDENT computation of route B's fit
pM = RB["Monod"]["p"]
ref_ode = mu_pred_ode("Monod", pM[0], pM[1:], A_BAL, PHI_BAL)
NS_BE = (50, 100, 200, 400)          # backward-Euler steps per interval, ratio 2
NEWTON_FLOOR = 1e-8                  # newton's backward_tol: the relative floor of a step


def be_study(ns=NS_BE, mum=None, floor=NEWTON_FLOOR, **kw):
    """Orders of the marcher from successive differences (observed_orders), per
    interval, with the round-off companion: the smallest finest successive
    difference over the newton floor NEWTON_FLOOR |mu|."""
    p_ = pM if mum is None else [mum, pM[1]]
    vals = {n: mu_pred_pymrm("Monod", p_[0], p_[1:], A_BAL, PHI_BAL, n, **kw) for n in ns}
    st = [observed_orders(lambda n, i=i: vals[n][i], ns) for i in range(NOBS)]
    fin = [o["orders"][-1] for o in st]
    diff = np.abs(vals[ns[-1]] - vals[ns[-2]])
    return dict(vals=vals, worst=float(max(fin, key=lambda q: abs(q - 1.0))),
                lo=float(min(min(o["orders"]) for o in st)),
                hi=float(max(max(o["orders"]) for o in st)),
                margin=float(np.min(diff/(floor*np.abs(vals[ns[-1]])))))


BES = be_study()
TM = pd.DataFrame([dict(nstep=n, **{"max rel. vs LSODA": float(np.max(np.abs(BES["vals"][n] - ref_ode)
                                                                   / np.abs(ref_ode)))})
                   for n in NS_BE])
show(TM.style.format({"max rel. vs LSODA": "{:.4e}"}).hide(axis="index").set_uuid("march"))
M["pymrm_batch_observed_order"] = BES["worst"]
M["pymrm_batch_order_roundoff_margin"] = BES["margin"]
print(f"orders from successive differences over {NS_BE}, every interval and triple:"
      f" {BES['lo']:.4f} to {BES['hi']:.4f};\nfinest triple, the one farthest from 1:"
      f" {BES['worst']:.4f}; smallest finest difference over the newton floor:"
      f" {BES['margin']:.2e}")

S0_b = Sof(Cc0[0], A_BAL, PHI_BAL)
c_end = BatchMarch(lambda S, X: pM[0]*S/(pM[1] + S), PHI_BAL, 400).march(
    [[S0_b, Cc0[0]]], dt[0])
inv0, inv1 = S0_b + Cc0[0]/PHI_BAL, c_end[0, 0] + c_end[0, 1]/PHI_BAL
M["pymrm_batch_conservation_relerr"] = float(abs(inv1 - inv0)/inv0)

pred_rich = lambda p: 2*mu_pred_pymrm("Monod", p[0], p[1:], A_BAL, PHI_BAL, 400) \
    - mu_pred_pymrm("Monod", p[0], p[1:], A_BAL, PHI_BAL, 200)
rr = least_squares(lambda p: pred_rich(p) - mu, pM, bounds=([1e-9, 1e-9], [50, 1e12]),
                   xtol=1e-12, ftol=1e-12, gtol=1e-12)
M["pymrm_refit_mum_relerr"] = float(abs(rr.x[0] - pM[0])/pM[0])
M["pymrm_refit_KS_relerr"] = float(abs(rr.x[1] - pM[1])/pM[1])
print(f"S + X/phi over interval 1: {inv0:.10g} -> {inv1:.10g}"
      f"   (rel. below 1e-13) (here {M['pymrm_batch_conservation_relerr']:.1e})")
print(f"\nroute B's Monod fit, re-derived through the pymrm marcher with Richardson")
print(f"extrapolation:  mu_m = {rr.x[0]:.5f}  K_S = {rr.x[1]:.3f}")
print(f"          LSODA: mu_m = {pM[0]:.5f}  K_S = {pM[1]:.3f}")
print(f"          relative: below 1e-4 for both (here"
      f" {M['pymrm_refit_mum_relerr']:.1e} and {M['pymrm_refit_KS_relerr']:.1e}),"
      f" the two fits' stopping rules")

# ---- WHAT THIS CROSS-CHECK CAN AND CANNOT SEE, measured.  Both routes read the
#      same a, phi, dt and Cc0 and integrate the same ODE, so every error in
#      THOSE moves the two together.  Three injections, each a failure mode
#      such a check might be thought to catch.
base_o = mu_pred_ode("Monod", pM[0], pM[1:], A_BAL, PHI_BAL)
blind = []
for lab, a_, ph_, dt_ in (("phi x 1.3   (a wrong yield)", A_BAL, PHI_BAL*1.3, dt),
                          ("a   x 1.2   (a wrong intercept)", A_BAL*1.2, PHI_BAL, dt),
                          ("dt  x 1.15  (a mis-set interval)", A_BAL, PHI_BAL, dt*1.15)):
    _sd, dt = dt, dt_                       # the marcher reads dt from the global
    try:
        o = mu_pred_ode("Monod", pM[0], pM[1:], a_, ph_)
        v = mu_pred_pymrm("Monod", pM[0], pM[1:], a_, ph_, 400)
    finally:
        dt = _sd
    blind.append(dict(injection=lab,
                      **{"prediction moves": float(np.max(np.abs(o - base_o)
                                                          / np.abs(base_o)))},
                      **{"LSODA vs pymrm": float(np.max(np.abs(v - o)/np.abs(o)))}))
TBL = pd.DataFrame([dict(injection="baseline (nothing injected)",
                         **{"prediction moves": 0.0},
                         **{"LSODA vs pymrm": float(np.max(np.abs(
                             mu_pred_pymrm("Monod", pM[0], pM[1:], A_BAL, PHI_BAL, 400)
                             - base_o)/np.abs(base_o)))})] + blind)
show(TBL.style.format({"prediction moves": "{:.3e}", "LSODA vs pymrm": "{:.3e}"})
     .hide(axis="index").set_uuid("blind"))
# the negative space of the cross-check: the claim blind_yield_gap_change_over_move
# checks the PROPERTY (the gap stays put while the prediction moves by orders)
# rather than pinning three tolerance-scale numbers.
BLIND_GAP = [float(v) for v in TBL["LSODA vs pymrm"]]
BLIND_MOVE = [float(v) for v in TBL["prediction moves"]]
print("\nThe two share no line of code, so a bug in EITHER INTEGRATOR moves them")
print("apart - which is what the order study and the 25-step break row test. They")
print("do NOT share their inputs' correctness: a yield wrong by 30 % moves the")
print(f"common prediction by {TBL['prediction moves'].iloc[1]:.1f} relative and"
      f" leaves their agreement at"
      f" {TBL['LSODA vs pymrm'].iloc[1]:.2e},\nagainst"
      f" {TBL['LSODA vs pymrm'].iloc[0]:.2e} with nothing injected. This check"
      f" measures the TIME\nDISCRETISATION and nothing else, and saying so is the"
      f" point of the table.")

# ---- the invariant holds for any growth law; a yield mismatch is what moves it
c_end_3x = BatchMarch(lambda S, X: 3*pM[0]*S/(pM[1] + S), PHI_BAL, 400).march(
    [[S0_b, Cc0[0]]], dt[0])
INV_3X_RELERR = float(abs((c_end_3x[0, 0] + c_end_3x[0, 1]/PHI_BAL) - inv0)/inv0)
print(f"\nS + X/phi is a TAUTOLOGY of the marcher, not a test of it: the source is")
print(f"(-m/phi, m) from one m, so backward Euler conserves it for ANY growth law.")
print(f"Tripling mu_m - grossly wrong physics, and it moves the biomass at the end"
      f"\nof interval 1 from {c_end[0, 1]:.2f} to {c_end_3x[0, 1]:.2f} - leaves the"
      f" invariant's relative change\nbelow 1e-13 (here {INV_3X_RELERR:.1e}), i.e. nowhere. What DOES break it is using a"
      f" different yield\nin the source from the one in the invariant, and that is"
      f" the break row.")

# %%
# ---- 3. plug flow == batch, and the grid order of the pymrm solve
S_START, X_START, TAU_VAL = Sof(Cc0[0], A_BAL, PHI_BAL), Cc0[0], 2.0
g_fit = lambda S, X: mumM*S/(KsM + S)
ref = batch_reference(g_fit, PHI_BAL, S_START, X_START, TAU_VAL)
NC_PFR = (200, 400, 800, 1600, 3200, 6400)     # cells, ratio 2


def pfr_study(ncs=NC_PFR, floor=NEWTON_FLOOR, **kw):
    """Outlet substrate (compute_boundary_values) on the ladder ncs, its orders from
    successive differences, the Richardson value against the batch reference, and
    the round-off companion: finest difference over the newton floor NEWTON_FLOOR S."""
    so = {n: PlugFlowFermenter(g_fit, PHI_BAL, S_START, X_START, TAU_VAL, ncell=n,
                               **kw).solve().outlet()[0] for n in ncs}
    st = observed_orders(lambda n: so[n], ncs)
    return dict(S=so, orders=st["orders"], order=float(st["orders"][-1]),
                relerr=float(abs(so[ncs[-1]] - ref[0])/ref[0]),
                rich=float(abs(st["extrapolated"] - ref[0])/ref[0]),
                margin=float(abs(so[ncs[-1]] - so[ncs[-2]])/(floor*abs(so[ncs[-1]]))))


PFS = pfr_study()
TP = pd.DataFrame([dict(ncell=n, S_outlet=PFS["S"][n], **{"rel. err": abs(PFS["S"][n] - ref[0])/ref[0]})
                   for n in NC_PFR])
show(TP.style.format({"S_outlet": "{:.9f}", "rel. err": "{:.4e}"})
     .hide(axis="index").set_uuid("pfr"))
M["pfr_batch_relerr_finest"] = PFS["relerr"]
M["pfr_observed_order"] = PFS["order"]
M["pfr_richardson_relerr"] = PFS["rich"]
M["pfr_order_roundoff_margin"] = PFS["margin"]
print(f"\nbatch reference at tau = {TAU_VAL}: S = {ref[0]:.3f}, X = {ref[1]:.3f}")
print(f"orders of the upwind plug-flow outlet from successive differences over {NC_PFR}:"
      f"\n  {', '.join(f'{q:.4f}' for q in PFS['orders'])}; finest {PFS['order']:.4f}")
print(f"Richardson value of the finest triple: {PFS['rich']:.3e} from the batch value;"
      f" finest difference\nover the newton floor: {PFS['margin']:.2e}.")
print("The outlet is a pure convective outflow ({'outflow': True}), so the face"
      " carries the last\ncell's value and compute_boundary_values reads exactly"
      " what the operator transported.")

# %%
# ---- the numbers README.md, meta.yaml and models.yaml quote that the cells
#      above print but do not register: each becomes a claim here, with a break
#      row in the table below, so every quoted number is checked by its kind.
def edge_fraction(sd, n_max=None):
    """fraction of the 2 NBOOT draws (both arms, on the seed boot_arm uses) whose
    scan optimum sits on an end of a search grid; n_max cuts the Moser n grid"""
    keep = np.array([n_max is None or t[1] <= n_max for t in thA_S])
    GS, TS = GA_S[keep], [t for t, k_ in zip(thA_S, keep) if k_]
    n_top = max(t[1] for t in TS)
    hits = 0
    for yt in (y_true, y_ctrl):
        Y = yt[None, :] + np.random.default_rng(SEED).normal(0.0, sd, (NBOOT, NOBS))
        y2 = (Y**2).sum(1)[:, None]
        i1 = np.argmin(y2 - (Y @ GA_M.T)**2/np.einsum("mi,mi->m", GA_M, GA_M)[None, :], axis=1)
        i2 = np.argmin(y2 - (Y @ GS.T)**2/np.einsum("mi,mi->m", GS, GS)[None, :], axis=1)
        on_s = np.array([TS[j][0] in (MOSER_KS[0], MOSER_KS[-1]) or TS[j][1] in (MOSER_N[0], n_top)
                         for j in i2])
        hits += int(np.sum((i1 == 0) | (i1 == len(RANGE["Monod"]) - 1) | on_s))
    return hits/(2*NBOOT)


def imputed(c4):
    """route A with row 4's CA_bar set to c4: the F test and the two RMS it quotes"""
    S_ = CAb.copy()
    S_[BAD_ROW] = c4
    R_ = {k: fit_A(k, S=S_) for k in ("Monod", "Blackman", "Moser")}
    F_ = float((R_["Monod"]["sse"] - R_["Moser"]["sse"])/(R_["Moser"]["sse"]/(NOBS - 3)))
    return {"balance_row4_imputed_CA": float(c4), "routeA_F_imputed": F_,
            "routeA_F_imputed_below_crit": FA_crit - F_,
            "rmsA_imputed_Monod": R_["Monod"]["rms"], "rmsA_imputed_Blackman": R_["Blackman"]["rms"],
            "routeA_imputed_blackman_margin": R_["Monod"]["rms"] - R_["Blackman"]["rms"]}


def end_of_interval(a_, b_):
    """the end-of-interval substrate a + b C_C,end of rows 3, 4 and 5, and the
    relative miss of rows 3 and 5 against their printed CA_bar"""
    e = {i: float(a_ + b_*Cc1[i]) for i in (2, 3, 4)}
    return {"eoi_row3": e[2], "eoi_row4": e[3], "eoi_row5": e[4],
            "eoi_row3_relgap": abs(e[2] - CAb[2])/CAb[2], "eoi_row5_relgap": abs(e[4] - CAb[4])/CAb[4]}


def boot_se(nboot=None):
    """the sampling errors the prose quotes beside the bootstrap estimates"""
    f_, c_ = boot_arm(y_true, sigma, nboot=nboot), boot_arm(y_ctrl, sigma, nboot=nboot)
    sA_ = RA["Monod"]["rms"]
    n_ = NBOOT if nboot is None else nboot
    return {"boot_size_forward_se": f_["size_se"], "boot_size_control_se": c_["size_se"],
            "boot_median_n_forward_se": f_["med_se"], "boot_median_n_control_se": c_["med_se"],
            "boot_p_of_observed_F_se": f_["pF_se"], "boot_p_of_observed_n_se": f_["pn_se"],
            "boot_size_forward_at_routeA_sigma_se": prop_se(boot_arm(y_true, sA_, nboot=nboot)["size"], n_),
            "boot_size_control_at_routeA_sigma_se": prop_se(boot_arm(y_ctrl, sA_, nboot=nboot)["size"], n_)}


NOM_SE = prop_se(0.05)                      # the binomial error of a nominal 5 % test
_sigB = "route B RMS (used above)"
if round(edge_fraction(sigma)*2*NBOOT) != out["forward"]["ntop"] + out["control"]["ntop"]:
    raise RuntimeError("edge_fraction does not reproduce the bootstrap's own edge count")
M.update({
    "routeB_monod_rms_drop": 1.0 - RB["Monod"]["rms"]/RA["Monod"]["rms"],
    "routeA_monod_last_margin": RA["Monod"]["rms"] - max(RA[k]["rms"] for k in TWO_PAR if k != "Monod"),
    "routeB_blackman_worst_margin": RB["Blackman"]["rms"] - max(RB[k]["rms"] for k in TWO_PAR if k != "Blackman"),
    "rmsA_null_over_monod": RA["first order"]["rms"]/RA["Monod"]["rms"],
    "rmsA_null_over_blackman": RA["first order"]["rms"]/RA["Blackman"]["rms"],
    "balance_R2_loo_min_other": min(loo[i][2] for i in range(NOBS) if i != BAD_ROW),
    "balance_R2_loo_max_other": max(loo[i][2] for i in range(NOBS) if i != BAD_ROW),
    "flagged_row": float(BAD_ROW + 1),
    "blackman_Nstar_shortfall": M["Nstar_blackman"] - 1.0,
    "rawlings_Dc_closed_form": float(Dc_r_closed),
    "boot_gridonly_size_forward": GRID_ONLY["forward"]["size"],
    "boot_gridonly_size_control": GRID_ONLY["control"]["size"],
    "boot_gridonly_median_n_forward": GRID_ONLY["forward"]["med"],
    "boot_gridonly_median_n_control": GRID_ONLY["control"]["med"],
    "boot_gridonly_p_of_observed_F": GRID_ONLY["forward"]["pF"],
    "boot_control_size_offset_se": abs(out["control"]["size"] - 0.05)/NOM_SE,
    "boot_size_forward_at_routeA_sigma": SIG_TABLE["route A RMS"][0],
    "boot_size_control_at_routeA_sigma": SIG_TABLE["route A RMS"][1],
    "boot_control_scale_free_se": abs(SIG_TABLE["route A RMS"][1] - SIG_TABLE[_sigB][1])
                                  / prop_se(SIG_TABLE[_sigB][1]),
    "boot_edge_fraction_main": edge_fraction(sigma),
    "boot_edge_fraction_routeA_sigma": edge_fraction(RA["Monod"]["rms"]),
    "blind_yield_gap_change_over_move": abs(BLIND_GAP[1] - BLIND_GAP[0])/BLIND_MOVE[1],
    "marcher_invariant_mu_tripled_relerr": INV_3X_RELERR,
    "Sf_rootfind_over_guarantee": SF_ROOTFIND_RELERR/SF_ROOTFIND_BOUND,
    **{k: v for k, v in imputed(CAB_IMPUTED).items()},
    **end_of_interval(A_BAL, B_BAL),
    **{k: out[a][q] for k, (a, q) in {
        "boot_size_forward_se": ("forward", "size_se"), "boot_size_control_se": ("control", "size_se"),
        "boot_median_n_forward_se": ("forward", "med_se"), "boot_median_n_control_se": ("control", "med_se"),
        "boot_p_of_observed_F_se": ("forward", "pF_se"), "boot_p_of_observed_n_se": ("forward", "pn_se")}.items()},
    "boot_size_forward_at_routeA_sigma_se": prop_se(SIG_TABLE["route A RMS"][0]),
    "boot_size_control_at_routeA_sigma_se": prop_se(SIG_TABLE["route A RMS"][1]),
})
if abs(M["routeA_F_imputed"] - FA_IMP) > 1e-12*FA_IMP:
    raise RuntimeError("imputed() does not reproduce route A's imputed F test")
print(f"route A F with row {BAD_ROW + 1} imputed at {CAB_IMPUTED:.2f}: {M['routeA_F_imputed']:.2f}"
      f" against {FA_crit:.2f}; RMS Monod {M['rmsA_imputed_Monod']:.5f},"
      f" Blackman {M['rmsA_imputed_Blackman']:.5f}")
print(f"route B: Monod's RMS {M['routeB_monod_rms_drop']:.0%} below route A's; the best"
      f" one-parameter null over Monod {M['rmsA_null_over_monod']:.2f}x, over Blackman"
      f" {M['rmsA_null_over_blackman']:.2f}x")
print(f"end-of-interval substrate, rows 3/4/5: {M['eoi_row3']:.2f}, {M['eoi_row4']:.2f},"
      f" {M['eoi_row5']:.2f}; rows 3 and 5 miss their printed CA_bar by"
      f" {M['eoi_row3_relgap']:.1%} and {M['eoi_row5_relgap']:.1%}")
print(f"bootstrap draws on a grid edge: {M['boot_edge_fraction_main']:.2e} of the draws at"
      f" route B's scatter, {M['boot_edge_fraction_routeA_sigma']:.2%} at route A's")

# %% [markdown]
# ### Defect injection
#
# Every metric that *can* be moved has a row that moves it. Metrics that cannot are
# named individually below and labelled structural. The table is built for this
# page's physics.

# %%
ROWS_L = []       # (row, note, study, targets): study() returns {metric: value}


def rec(name, metric, fn, note):
    """one row, one target: fn() recomputes the metric with the defect injected"""
    ROWS_L.append((name, note, lambda: {metric: fn()}, [metric]))


def rec_many(name, note, fn, targets):
    """one row, several targets: fn() returns {metric: value} from one study"""
    ROWS_L.append((name, note, fn, list(targets)))


def r2_line(X, keep):
    """coefficients, R^2 and residuals of the line CA_bar = a + b X on the rows keep"""
    A = np.vstack([np.ones(len(keep)), X[keep]]).T
    c = np.linalg.lstsq(A, CAb[keep], rcond=None)[0]
    r = CAb[keep] - A @ c
    return (float(c[0]), float(c[1]),
            float(1 - np.sum(r**2)/np.sum((CAb[keep] - CAb[keep].mean())**2)), r)


KEEP6 = np.array([i for i in range(NOBS) if i != BAD_ROW])


def f_moser(ra_monod, ra_moser):
    """F statistic of Moser against Monod, 1 and NOBS - 3 degrees of freedom"""
    return (ra_monod["sse"] - ra_moser["sse"])/(ra_moser["sse"]/(NOBS - 3))


# --- the identities, closed forms and the Blackman numerics
rec_many("Blackman switch placed at S = 4 K_S", "the closed forms are functions of the switch",
         lambda: dict(zip(("sym_blackman_gap_at_switch_over_mum", "sym_blackman_gap_lower_over_mum"),
                          map(float, sym_gaps(4)))),
         ["sym_blackman_gap_at_switch_over_mum", "sym_blackman_gap_lower_over_mum"])
rec("Blackman gap read at S = K_S instead of the switch S = 2 K_S",
    "num_blackman_gap_switch_relerr", lambda: abs(blackman_numeric(at=Ks_v)[2] - 1/3)/(1/3),
    "at S = K_S the two laws agree exactly; the gap lives at the switch")
rec_many("lower-branch extremum read off a 64-point scan, not root-found",
         "an extremum is the root of its derivative, not a grid read-off",
         lambda: (lambda o: {"num_blackman_Sstar_over_Ks_relerr": abs(o[0] - S_STAR_OVER_KS)/S_STAR_OVER_KS,
                             "num_blackman_gap_lower_relerr":
                                 abs(o[1] - M["sym_blackman_gap_lower_over_mum"])
                                 / M["sym_blackman_gap_lower_over_mum"]})(blackman_numeric(refine=False)),
         ["num_blackman_Sstar_over_Ks_relerr", "num_blackman_gap_lower_relerr"])
# --- the data reduction and the balance line
rec("balance line on the ARITHMETIC mean of C_C, not the log mean", "balance_R2_six_rows",
    lambda: r2_line(Xam, KEEP6)[2],
    "the exact time average under exponential growth is the LOG mean; declared convention")
rec("the flagged row kept in the balance line", "balance_R2_six_rows", lambda: bal_line()[2],
    "this IS the finding: one cell moves R^2")
rec_many("balance line drops row 1 instead of the flagged row",
         "the leave-one-out sweep says WHICH row, and it is not this one",
         lambda: (lambda b: {"balance_phi": -1.0/b[1], "balance_a": b[0]})(bal_line(drop=(0,))),
         ["balance_phi", "balance_a"])
rec("mu from a forward difference of C_C over the ARITHMETIC-mean biomass",
    "rmsA_Monod", lambda: fit_A("Monod", y=(Cc1 - Cc0)/dt/Xam)["rms"],
    "d(ln C)/dt = mu is exact; a difference quotient over the wrong mean is not")
rec("substrate swing taken as the biomass swing (yield omitted)", "max_swing_over_CAbar",
    lambda: float(np.max((Cc1 - Cc0)/CAb)),
    "the swing is substrate consumed, (C_C,end - C_C,start)/phi")
# --- route A / route B
rec("route B yield phi perturbed by 5 %", "routeB_n_moser",
    lambda: fit_B("Moser", phi=PHI_BAL*1.05, p0=RB["Moser"]["p"])["p"][2],
    "perturbs the computation, the row that could expose a wrong route-B baseline")
rec("Contois given the ARITHMETIC-mean biomass in route A", "rmsA_Contois",
    lambda: fit_A("Contois", X=Xam)["rms"], "declared convention")
rec("Blackman fitted from a single local start, no global scan", "rmsA_Blackman",
    lambda: (lambda r: float(np.sqrt(np.sum(r.fun**2)/NOBS)))(
        least_squares(lambda p: p[0]*np.minimum(1.0, CAb/(2*p[1])) - mu,
                      [0.8, 20.0], bounds=([1e-6, 1e-6], [10, 1e4]),
                      xtol=1e-15, ftol=1e-15, gtol=1e-15)),
    "the kinked objective has local minima; this is why the page scans globally")
rec("Moser's exponent frozen at n = 1", "routeA_F_moser_vs_monod",
    lambda: f_moser(RA["Monod"], fit_A("Moser", n_fixed=1.0)),
    "n = 1 IS Monod, so the F statistic must collapse to zero")
rec("route-A Moser scan restricted to n <= 1.5", "routeA_n_moser",
    lambda: fit_A("Moser", n_grid=MOSER_N[MOSER_N <= 1.5])["p"][2],
    "the optimum lies beyond the restricted grid, so the scan raises at its edge")
rec("F critical value with Moser's extra parameter not counted (1, 3 dof)", "F_crit_95",
    lambda: float(f_dist.ppf(0.95, 1, NOBS - 4)), "the residual degrees of freedom are NOBS - 3")
# --- the Contois map and the degenerate separations
rec_many("Contois map evaluated with the wrong yield (phi -> 2 phi)",
         "companion for the exactly-zero symbolic residual",
         lambda: {"contois_map_mum_relerr":
                      abs(mumC/(1 - KsxC*2*PHI_BAL) - RB["Monod"]["p"][0])/RB["Monod"]["p"][0],
                  "contois_map_KS_relerr":
                      abs(KsxC*2*PHI_BAL*A_BAL/(1 - KsxC*2*PHI_BAL) - RB["Monod"]["p"][1])
                      / RB["Monod"]["p"][1]},
         ["contois_map_mum_relerr", "contois_map_KS_relerr"])
rec_many("Contois separation measured on the ROUTE-A design",
         "the degeneracy is exact for the forward model only; the reduction manufactures a separation",
         lambda: (lambda d: {"sep_contois": d, "sep_over_scatter_contois": d/s_obs})(
             fit_A("Contois", y=RA["Monod"]["p"][0]*G["Monod"](RA["Monod"]["p"][1:], CAb, Xlm))["rms"]),
         ["sep_contois", "sep_over_scatter_contois"])
rec_many("Moser and Contois fits without their least_squares refinement",
         "the degenerate separations and the Contois-Monod gap are the fits' exit tolerances",
         lambda: (lambda dm, rc: {"sep_moser": dm, "sep_over_scatter_moser": dm/s_obs,
                                  "contois_minus_monod_rms": abs(rc - RB["Monod"]["rms"])})(
             fit_B("Moser", y=y_true, p0=RB["Moser"]["p"], refine=False)["rms"],
             fit_B("Contois", refine=False)["rms"]),
         ["sep_moser", "sep_over_scatter_moser", "contois_minus_monod_rms"])
# --- separations and the bootstrap
rec_many("N* computed at 99 % instead of 95 %", "declared convention",
         lambda: (lambda n_: {"Nstar_blackman": n_, "blackman_Nstar_shortfall": n_ - 1.0})(
             float(chi2.ppf(0.99, NOBS - 2)*s_obs**2/(NOBS*M["sep_blackman"]**2))),
         ["Nstar_blackman", "blackman_Nstar_shortfall"])
rec("scatter taken from route A instead of route B", "Nstar_tessier",
    lambda: float(chi2.ppf(0.95, NOBS-2)*RA["Monod"]["rms"]**2/(NOBS*M["sep_tessier"]**2)),
    "declared convention")
rec("bootstrap at half the noise level", "boot_size_forward_truth",
    lambda: boot_arm(y_true, sigma/2)["size"],
    "the size distortion is a bias-to-noise ratio, so it must grow as noise falls")
rec("bootstrap truth taken from route A's own model", "boot_median_n_forward_truth",
    lambda: boot_arm(y_ctrl, sigma)["med"], "the control, run through the forward arm's call")
# --- Rawlings, the chemostat and plug flow
rec("Rawlings' printed y = 1 replaced by 2", "rawlings_10_19_max_relerr",
    lambda: float(max(max(abs(v - c)/c for v, c in
                          zip(Chemostat(g_r, 2.0, Sf_r, D).steady(),
                              (D*Ks_r/(mum_r - D), y_r*(Sf_r - D*Ks_r/(mum_r - D)))))
                      for D in (0.25, 0.5, 0.75))),
    "reproducing a printed number must fail when the printed number is changed")
rec("Rawlings' printed S_f = 5 replaced by 8", "rawlings_10_19_max_relerr",
    lambda: float(max(max(abs(v - c)/c for v, c in
                          zip(Chemostat(g_r, y_r, 8.0, D).steady(),
                              (D*Ks_r/(mum_r - D), y_r*(5.0 - D*Ks_r/(mum_r - D)))))
                      for D in (0.25, 0.5, 0.75))),
    "the same, for the steady-state pair")
rec("Rawlings' printed K_s = 1 replaced by 2", "rawlings_Dc_relerr",
    lambda: abs(washout(lambda S, X: mum_r*S/(2.0 + S), y_r, Sf_r) - Dc_r_closed)/Dc_r_closed,
    "reproducing a printed number must fail when the printed number is changed")
rec_many("washout threshold loosened 1e-6 -> 1e-3 of y S_f",
         "declared convention, and the convergence evidence for the quoted D_c",
         lambda: {"Dc_monod_transient_relerr":
                      abs(washout(gM_, y_yield, Sf_test, wash=1e-3) - Dc_monod(Sf_test))/Dc_monod(Sf_test),
                  "Dc_contois_transient_relerr":
                      abs(washout(gC_, y_yield, Sf_test, wash=1e-3) - Dc_contois)/Dc_contois},
         ["Dc_monod_transient_relerr", "Dc_contois_transient_relerr"])
rec("plug flow with nu = 1 (cylindrical) instead of nu = 0", "pfr_batch_relerr_finest",
    lambda: pfr_study(ncs=(400, 800, 1600), nu=1)["relerr"],
    "nu is geometry, and a wrong one is invisible in a smooth profile")
rec("plug flow with the substrate sign flipped in the source", "pfr_batch_relerr_finest",
    lambda: pfr_study(ncs=(400, 800, 1600), substrate_sign=-1.0)["relerr"],
    "substrate must be consumed as biomass is made")
rec("plug-flow newton stopped on a 1e-4 step alone (floor 1e-4)", "pfr_order_roundoff_margin",
    lambda: pfr_study(newton_kw=dict(tol=1e-4, backward_tol=None), floor=1e-4)["margin"],
    "the order needs the finest difference far above the solver floor")
rec_many("every plug-flow number read off the three COARSEST grids",
         "an order quoted from an unconverged ladder is the classic way to publish a wrong one",
         lambda: (lambda o: {"pfr_observed_order": o["order"], "pfr_richardson_relerr": o["rich"]})(
             pfr_study(ncs=(25, 50, 100))),
         ["pfr_observed_order", "pfr_richardson_relerr"])
rec_many("the pymrm marcher's order read off 6/12/24 steps, and its refit at 25 steps",
         "the second route must be converged before it can confirm the first",
         lambda: {"pymrm_batch_observed_order": be_study(ns=(6, 12, 24))["worst"],
                  **(lambda r: {"pymrm_refit_mum_relerr": abs(r.x[0] - pM[0])/pM[0],
                                "pymrm_refit_KS_relerr": abs(r.x[1] - pM[1])/pM[1]})(
                      least_squares(lambda p: mu_pred_pymrm("Monod", p[0], p[1:], A_BAL, PHI_BAL, 25) - mu,
                                    pM, bounds=([1e-9, 1e-9], [50, 1e12]),
                                    xtol=1e-12, ftol=1e-12, gtol=1e-12))},
         ["pymrm_batch_observed_order", "pymrm_refit_mum_relerr", "pymrm_refit_KS_relerr"])
rec("marcher newton stopped on a 1e-3 step alone (floor 1e-3)", "pymrm_batch_order_roundoff_margin",
    lambda: be_study(newton_kw=dict(tol=1e-3, backward_tol=None), floor=1e-3)["margin"],
    "the order needs the finest difference far above the solver floor")
def _wrong_yield_invariant(mu_scale, phi_src):
    """Relative change of S + X/phi (phi = PHI_BAL) when the marcher's source uses phi_src."""
    c = BatchMarch(lambda S, X: mu_scale*pM[0]*S/(pM[1] + S), phi_src, 400).march(
        [[S0_b, Cc0[0]]], dt[0])
    return abs((c[0, 0] + c[0, 1]/PHI_BAL) - inv0)/inv0


rec_many("pymrm batch marcher with the yield used only in the source, not the state",
         "the conservation check holds only if BOTH sides use the same yield",
         lambda: {"pymrm_batch_conservation_relerr": _wrong_yield_invariant(1.0, PHI_BAL*1.01),
                  "marcher_invariant_mu_tripled_relerr": _wrong_yield_invariant(3.0, PHI_BAL*1.01)},
         ["pymrm_batch_conservation_relerr", "marcher_invariant_mu_tripled_relerr"])

# --- SWEEP ROWS: one injection, recorded against every metric it touches
RA_IMP_KEYS = ([f"rmsA_{k.replace(' ', '_')}" for k in LAWS + ["first order"]]
               + ["routeA_n_moser", "routeA_F_moser_vs_monod", "routeA_monod_last_margin",
                  "rmsA_null_over_monod", "rmsA_null_over_blackman"])
rec_many("every route-A fit with the flagged cell replaced by its balance-implied value",
         "the single flagged cell, propagated through all seven route-A fits and route A's F test;"
         " a labelled INFERENCE",
         lambda: (lambda R: {**{f"rmsA_{k.replace(' ', '_')}": R[k]["rms"] for k in LAWS + ["first order"]},
                             "routeA_n_moser": R["Moser"]["p"][2],
                             "routeA_F_moser_vs_monod": f_moser(R["Monod"], R["Moser"]),
                             "routeA_monod_last_margin": R["Monod"]["rms"]
                             - max(R[k]["rms"] for k in TWO_PAR if k != "Monod"),
                             "rmsA_null_over_monod": R["first order"]["rms"]/R["Monod"]["rms"],
                             "rmsA_null_over_blackman": R["first order"]["rms"]/R["Blackman"]["rms"]})(
             {k: fit_A(k, S=CAb_imp) for k in LAWS + ["first order"]}),
         RA_IMP_KEYS)
BOOT_KEYS = {"boot_size_forward_truth": ("forward", "size"), "boot_size_control_truth": ("control", "size"),
             "boot_median_n_forward_truth": ("forward", "med"),
             "boot_median_n_control_truth": ("control", "med"),
             "boot_p_of_observed_F": ("forward", "pF"), "boot_p_of_observed_n": ("forward", "pn")}
rec_many("the whole bootstrap fitted on the search grid alone, with no per-draw refinement",
         "an unrefined Moser fit cannot reach its own optimum, so it deflates F on every draw",
         lambda: {k: boot_arm(y_true if a == "forward" else y_ctrl, sigma, rounds=0)[q]
                  for k, (a, q) in BOOT_KEYS.items()},
         list(BOOT_KEYS))
_mu_fd = (Cc1 - Cc0)/dt/Xam            # a forward difference over the wrong mean
rec_many("mu itself from a forward difference over the ARITHMETIC-mean biomass",
         "d(ln C)/dt = mu is EXACT; the one-parameter nulls see only mu",
         lambda: {"rmsA_zero_order": fit_A("zero order", y=_mu_fd)["rms"],
                  "rmsB_zero_order": fit_B("zero order", y=_mu_fd)["rms"]},
         ["rmsA_zero_order", "rmsB_zero_order"])
RB7_KEYS = ([f"rmsB_{k.replace(' ', '_')}" for k in LAWS + ["first order"]]
            + ["routeB_monod_mu_m", "routeB_monod_K_S", "routeB_F_moser_vs_monod",
               "routeA_over_routeB_monod_rms", "routeB_monod_rms_drop", "routeB_blackman_worst_margin"])
rec_many("every route-B fit with the ALL-SEVEN balance line",
         "phi and a carry the flagged row into every route-B number at once",
         lambda: (lambda R: {**{f"rmsB_{k.replace(' ', '_')}": R[k]["rms"] for k in LAWS + ["first order"]},
                             "routeB_monod_mu_m": R["Monod"]["p"][0], "routeB_monod_K_S": R["Monod"]["p"][1],
                             "routeB_F_moser_vs_monod": f_moser(R["Monod"], R["Moser"]),
                             "routeA_over_routeB_monod_rms": RA["Monod"]["rms"]/R["Monod"]["rms"],
                             "routeB_monod_rms_drop": 1.0 - R["Monod"]["rms"]/RA["Monod"]["rms"],
                             "routeB_blackman_worst_margin": R["Blackman"]["rms"]
                             - max(R[k]["rms"] for k in TWO_PAR if k != "Blackman")})(
             {k: fit_B(k, a=a7, phi=-1.0/b7, p0=RB[k]["p"]) for k in LAWS + ["first order"]}),
         RB7_KEYS)
_idx3 = np.array([0, 3, 6])


def _sep3():
    """the two finite separations and their N*, on 3 of the 7 design points"""
    out_ = {}
    for k in ["Blackman", "Tessier"]:
        d2 = float(np.mean(fit_B(k, y=y_true, p0=RB[k]["p"])["resid"][_idx3]**2))
        out_[f"sep_{k.lower()}"] = np.sqrt(d2)
        out_[f"sep_over_scatter_{k.lower()}"] = np.sqrt(d2)/s_obs
        out_[f"Nstar_{k.lower()}"] = float(chi2.ppf(0.95, NOBS-2)*s_obs**2/(3*d2))
    return out_


rec_many("every separation measured at 3 of the 7 design points",
         "identifiability is a property of the DESIGN", _sep3,
         [f"{p_}_{k}" for k in ("blackman", "tessier") for p_ in ("sep", "sep_over_scatter", "Nstar")])
rec_many("the whole bootstrap re-drawn at half the noise level",
         "every bootstrap number is a property of the noise level and the seed",
         lambda: (lambda f_, c_: {"boot_p_of_observed_F": f_["pF"], "boot_p_of_observed_n": f_["pn"],
                                  "boot_size_control_truth": c_["size"],
                                  "boot_median_n_control_truth": c_["med"]})(
             boot_arm(y_true, sigma/2), boot_arm(y_ctrl, sigma/2, seed=SEED+1)),
         ["boot_p_of_observed_F", "boot_p_of_observed_n", "boot_size_control_truth",
          "boot_median_n_control_truth"])


def _bias_at(p_):
    """the reduction bias of a Monod curve with constants p_"""
    b_ = p_[0]*G["Monod"](p_[1:], CAb, Xlm) - mu_pred_ode("Monod", p_[0], p_[1:], A_BAL, PHI_BAL)
    ex = np.abs(np.where(np.arange(NOBS) != BAD_ROW, b_, 0.0))
    return {"bias_max_abs": float(np.max(np.abs(b_))), "bias_max_abs_excl_bad_row": float(np.max(ex)),
            "bias_over_routeA_rms": float(np.max(ex))/RA["Monod"]["rms"]}


rec_many("the reduction bias evaluated at the ROUTE-B constants",
         "the bias is a property of the curve it is evaluated on",
         lambda: _bias_at(RB["Monod"]["p"]),
         ["bias_max_abs", "bias_max_abs_excl_bad_row", "bias_over_routeA_rms"])
rec_many("the whole balance test on the ARITHMETIC-mean biomass",
         "declared convention: the exact time average under exponential growth is the LOG mean",
         lambda: (lambda b6, b7_: {"balance_R2_all_rows": b7_[2],
                                   "balance_row4_residual": float(CAb[BAD_ROW] - (b6[0] + b6[1]*Xam[BAD_ROW])),
                                   "balance_worst_other_residual": float(np.max(np.abs(b6[3]))),
                                   "row4_end_of_interval_relgap":
                                       float(abs(CAb[BAD_ROW] - (b6[0] + b6[1]*Cc1[BAD_ROW]))/CAb[BAD_ROW]),
                                   "balance_R2_loo_min_other": min(_loo_am),
                                   "balance_R2_loo_max_other": max(_loo_am)})(
             r2_line(Xam, KEEP6), r2_line(Xam, np.arange(NOBS))),
         ["balance_R2_all_rows", "balance_row4_residual", "balance_worst_other_residual",
          "row4_end_of_interval_relgap", "balance_R2_loo_min_other", "balance_R2_loo_max_other"])
_loo_am = [r2_line(Xam, np.array([j for j in range(NOBS) if j != i]))[2]
           for i in range(NOBS) if i != BAD_ROW]   # leave-one-out on the arithmetic mean


def _chemo_A():
    """the chemostat predictions built from the route-A constants"""
    DcA = lambda Sf: RA["Monod"]["p"][0]*Sf/(RA["Monod"]["p"][1] + Sf)
    ratA = lambda Sf: RA["Contois"]["p"][0]/DcA(Sf)
    return {"Dc_ratio_at_Sf20": float(ratA(20.0)),
            "Sf_crossing": one_root(lambda Sf: ratA(Sf) - 1.0, SF_LO, SF_HI, log=True, xtol=SF_XTOL),
            "Sf_for_10pct_gap": one_root(lambda Sf: ratA(Sf) - 1.10, SF_LO, SF_HI, log=True, xtol=SF_XTOL),
            "Dc_monod_closed_form": float(DcA(20.0))}


rec_many("the chemostat prediction built from the ROUTE-A constants",
         "which reduction you fit with changes what you predict for an experiment not run",
         lambda: (lambda d: {**d, "Sf_crossing_vs_a_relerr": abs(d["Sf_crossing"] - A_BAL)/A_BAL})(_chemo_A()),
         ["Dc_ratio_at_Sf20", "Sf_crossing", "Sf_for_10pct_gap", "Dc_monod_closed_form",
          "Sf_crossing_vs_a_relerr"])

# --- the quoted numbers registered in the cell above the table
rec("leave-one-out sweep picks the row whose removal helps LEAST", "flagged_row",
    lambda: float(min(range(NOBS), key=lambda i: loo[i][2]) + 1),
    "the flagged row is the one whose removal raises R^2 most")
IMP_KEYS = list(imputed(CAB_IMPUTED))
rec_many("the row-4 imputation left at the printed cell", "the inference is the balance-implied value",
         lambda: imputed(CAb[BAD_ROW]), IMP_KEYS)
EOI_KEYS = list(end_of_interval(A_BAL, B_BAL))
rec_many("end-of-interval values read off the ALL-SEVEN balance line",
         "the substitution uses the six-row line that leaves the flagged row out",
         lambda: end_of_interval(a7, b7), EOI_KEYS)
rec("Rawlings' printed K_s = 1 replaced by 2 in their D_c formula", "rawlings_Dc_closed_form",
    lambda: float(mum_r*Sf_r/(2.0 + Sf_r)), "the closed form is evaluated at the printed set")
rec("the control's distance from nominal read off the grid-only bootstrap",
    "boot_control_size_offset_se", lambda: abs(boot_arm(y_ctrl, sigma, rounds=0)["size"] - 0.05)/NOM_SE,
    "an unrefined fit deflates F, so the grid-only control looks conservative")
GO_KEYS = ["boot_gridonly_size_forward", "boot_gridonly_p_of_observed_F",
           "boot_gridonly_size_control", "boot_gridonly_median_n_control"]
rec_many("the grid-only bootstrap re-drawn at half the noise level",
         "the grid-only numbers are a property of the noise level and the seed",
         lambda: (lambda f_, c_: dict(zip(GO_KEYS, (f_["size"], f_["pF"], c_["size"], c_["med"]))))(
             boot_arm(y_true, sigma/2, rounds=0), boot_arm(y_ctrl, sigma/2, seed=SEED+1, rounds=0)),
         GO_KEYS)
rec("grid-only bootstrap truth taken from route A's own model", "boot_gridonly_median_n_forward",
    lambda: boot_arm(y_ctrl, sigma, rounds=0)["med"], "the control, run through the forward arm's call")
rec("the route-A-scatter forward arm at half its noise level", "boot_size_forward_at_routeA_sigma",
    lambda: boot_arm(y_true, RA["Monod"]["rms"]/2)["size"],
    "the forward size is a bias-to-noise ratio")
rec_many("the route-A-scatter control arm drawn from the forward truth",
         "only the correctly specified control is scale-free",
         lambda: (lambda fA_, fB_: {"boot_size_control_at_routeA_sigma": fA_,
                                    "boot_control_scale_free_se": abs(fA_ - fB_)/prop_se(fB_)})(
             boot_arm(y_true, RA["Monod"]["rms"])["size"], boot_arm(y_true, sigma)["size"]),
         ["boot_size_control_at_routeA_sigma", "boot_control_scale_free_se"])
SE_KEYS = ["boot_size_forward_se", "boot_size_control_se", "boot_median_n_forward_se",
           "boot_median_n_control_se", "boot_p_of_observed_F_se", "boot_p_of_observed_n_se",
           "boot_size_forward_at_routeA_sigma_se", "boot_size_control_at_routeA_sigma_se"]
rec_many("every bootstrap sampling error at 500 draws instead of 2000",
         "a sampling error scales as one over the square root of the draw count",
         lambda: boot_se(nboot=500), SE_KEYS)
rec_many("the Moser exponent grid cut at n <= 2",
         "a grid too narrow for the draws puts more of them on its edge",
         lambda: {"boot_edge_fraction_main": edge_fraction(sigma, n_max=2.0),
                  "boot_edge_fraction_routeA_sigma": edge_fraction(RA["Monod"]["rms"], n_max=2.0)},
         ["boot_edge_fraction_main", "boot_edge_fraction_routeA_sigma"])


def _blind_one_sided():
    """the cross-check with the wrong yield given to the LSODA route only"""
    o = mu_pred_ode("Monod", pM[0], pM[1:], A_BAL, PHI_BAL*1.3)
    v = mu_pred_pymrm("Monod", pM[0], pM[1:], A_BAL, PHI_BAL, 400)
    return abs(float(np.max(np.abs(v - o)/np.abs(o))) - BLIND_GAP[0])/BLIND_MOVE[1]


rec("the wrong yield given to the LSODA route only", "blind_yield_gap_change_over_move",
    _blind_one_sided, "the check is blind only to an error both routes share")

# Every row is declared before any runs; a study that raises records nan for
# every claim it targets and stays in the table.
claims.define_rows(list(dict.fromkeys(r[0] for r in ROWS_L)))
for name, note, study, targets in ROWS_L:
    try:
        vals = study()
        err = ""
    except Exception as exc:                 # noqa: BLE001 - a failed study is a result
        vals = {t: np.nan for t in targets}
        err = f"raised {type(exc).__name__}: {str(exc)[:60]}"
    for t in targets:
        base = M[t]
        v = float(vals[t])
        rel = abs(v - base)/max(abs(base), 1e-300) if np.isfinite(v) else np.inf
        BREAK.append(dict(row=name, metric=t, baseline=base, broken=v, **{"rel. move": rel},
                          note=(err or note)))

TBK = pd.DataFrame(BREAK)
show(TBK.style.format({"baseline": "{:.6g}", "broken": "{:.6g}", "rel. move": "{:.3e}"})
     .hide(axis="index").set_uuid("brk"))
print(f"\n{TBK['row'].nunique()} rows defined and run, {len(TBK)} (row, claim) entries;"
      f" rows whose study raised: {sorted(set(TBK.row[TBK.note.str.startswith('raised')])) or 'none'}")

# %%
# The quantities this page PROVES are exactly zero (identities, nesting, exact
# arithmetic, optimiser and root-finder exits).  Every one of them is a
# registered BOUND in the claims cell below, with a basis and a limit that holds
# on any machine, so none is ever compared as a value and the set of claims
# never depends on where a round-off number happens to land.
PROVED_ZERO = [
    "sym_monod_is_langmuir_residual", "sym_monod_half_rate_residual",
    "sym_blackman_half_rate_residual", "sym_moser_at_n1_residual",
    "sym_contois_is_monod_residual", "num_blackman_Sstar_over_Ks_relerr",
    "num_blackman_gap_lower_relerr", "num_blackman_gap_switch_relerr",
    "rawlings_10_19_max_relerr", "pymrm_batch_conservation_relerr",
    "sep_moser", "sep_over_scatter_moser", "sep_contois",
    "sep_over_scatter_contois", "contois_minus_monod_rms",
    "contois_map_mum_relerr", "contois_map_KS_relerr", "Sf_crossing_vs_a_relerr"]
if not set(PROVED_ZERO) <= set(M):
    raise RuntimeError(f"proved-zero quantities not computed: {set(PROVED_ZERO) - set(M)}")

# %%
# Every metric is registered ONCE with how it is checked (docs/claims-design.md):
# a value against the committed baseline within its rtol, a bound against its
# limit on every run and on any machine, an exact result by equality.  The set
# of names is fixed and never depends on a value.  Residuals, identities,
# two-route gaps and optimiser or root-finder exits are BOUNDS whatever their
# size, and the prose quotes the bound, never their round-off digits.
B = {   # name -> bound arguments; every name in PROVED_ZERO is here
    # -- symbolic identities: sympy simplifies each to exactly 0
    **{k: dict(max=0.0, basis="identity", structural=True,
               note="sympy residual, simplified to exactly 0")
       for k in ("sym_monod_is_langmuir_residual", "sym_monod_half_rate_residual",
                 "sym_blackman_half_rate_residual", "sym_moser_at_n1_residual",
                 "sym_contois_is_monod_residual")},
    # -- exact arithmetic and root-finds against closed forms
    "num_blackman_gap_switch_relerr": dict(
        max=1e-14, basis="identity",
        note="relative to 1/3; exact arithmetic, 1 - 2/3 against 1/3, measured 1.7e-16"),
    "num_blackman_gap_lower_relerr": dict(
        max=1e-14, basis="tolerance",
        note="relative to the closed-form gap; find_roots at xtol 1e-14 on a stationary"
             " point, so the error is second order in xtol; measured 1.6e-16"),
    "num_blackman_Sstar_over_Ks_relerr": dict(
        max=1e-13, basis="tolerance",
        note="relative to S*/K_S; the root's guarantee xtol + rtol|x| is 2.5e-14"
             " relative here; measured 0"),
    "rawlings_10_19_max_relerr": dict(
        max=1e-9, basis="tolerance",
        note="relative to the closed-form steady state; LSODA at rtol 1e-11 marched"
             " to t = 40000; measured exactly 0"),
    "pymrm_batch_conservation_relerr": dict(
        max=1e-13, basis="identity",
        note="relative to S + X/phi; backward Euler keeps a linear invariant of the"
             " source to round-off; measured 1.4e-16 to 2.9e-16 across BLAS kernels,"
             " the bound is about 350x the largest"),
    # -- optimiser exits on quantities the identities make exactly zero
    "sep_moser": dict(max=SEP_BOUND, basis="tolerance", unit="1/h",
                      note="Moser contains Monod at n = 1; least_squares at xtol 1e-13;"
                           " measured 2.8e-14 to 8.3e-14 across BLAS kernels and above 1e-12 on"
                           " other machines; 1e-9 is also the degeneracy threshold"),
    "sep_over_scatter_moser": dict(max=SEP_SCATTER_BOUND, basis="tolerance",
                                   note="sep_moser over the route-B scatter"),
    "sep_contois": dict(max=SEP_BOUND, basis="tolerance", unit="1/h",
                        note="Contois is Monod under a constant yield; least_squares"
                             " at xtol 1e-13; measured 1.2e-13 to 1.5e-13 across BLAS"
                             " kernels; 1e-9 is also the degeneracy threshold"),
    "sep_over_scatter_contois": dict(max=SEP_SCATTER_BOUND, basis="tolerance",
                                     note="sep_contois over the route-B scatter"),
    "contois_minus_monod_rms": dict(
        max=RMS_GAP_BOUND, basis="tolerance", unit="1/h",
        note="one curve fitted twice; two least_squares exits; measured 7.6e-15"
             " to 1.5e-13 across BLAS kernels, the bound is about 650x the largest"),
    "contois_map_mum_relerr": dict(
        max=MAP_BOUND, basis="tolerance",
        note="relative to the Monod mu_m; two least_squares exits; measured 2.7e-8"
             " to 3.7e-8 across BLAS kernels"),
    "contois_map_KS_relerr": dict(
        max=MAP_BOUND, basis="tolerance",
        note="relative to the Monod K_S; two least_squares exits; measured 1.1e-7"
             " to 2.9e-7 across BLAS kernels, the bound is about 3400x the largest"),
    "Sf_crossing_vs_a_relerr": dict(
        max=MAP_BOUND, basis="tolerance",
        note="relative to a; S_f = a identically, the gap is the Contois fit's two"
             " exit tolerances (decomposed on the page); measured 1.0e-7 to 1.7e-7"
             " across BLAS kernels, the bound is about 5800x the largest"),
    # -- convergence orders and the second computation of the fit
    "pymrm_batch_observed_order": dict(
        min=0.995, max=1.005, basis="theory",
        note="backward Euler is first order in dt; finest triple 100/200/400 of the"
             " successive differences, the interval farthest from 1"),
    "pfr_observed_order": dict(
        min=0.99, max=1.01, basis="theory",
        note="first-order upwind; finest triple 1600/3200/6400 of the successive"
             " differences of the outlet from compute_boundary_values"),
    "pfr_order_roundoff_margin": dict(
        min=1e2, basis="tolerance",
        note="finest successive difference of the outlet over newton's backward_tol"
             " 1e-8 times the outlet value: the order is read far above the solver floor"),
    "pymrm_batch_order_roundoff_margin": dict(
        min=1e2, basis="tolerance",
        note="smallest finest successive difference of the seven predictions over"
             " newton's backward_tol 1e-8 times the prediction"),
    "pfr_richardson_relerr": dict(
        max=1e-4, basis="theory",
        note="Richardson value of the finest triple (observed order) against the batch"
             " state: a discretisation residual"),
    "pymrm_refit_mum_relerr": dict(
        max=1e-4, basis="chosen",
        note="two-route gap, LSODA fit against Richardson-extrapolated backward Euler"
             " at 200/400 steps; set by the two fits' stopping rules, measured 7.0e-7 to"
             " 8.2e-7 across BLAS kernels; the bound is the prose claim and 1e2 below"
             " the 25-step row"),
    "pymrm_refit_KS_relerr": dict(
        max=1e-4, basis="chosen",
        note="as for mu_m; measured 2.5e-6 to 3.4e-6 across BLAS kernels; the bound"
             " is the prose claim and below the unconverged 25-step refit"),
    # -- verdicts and identities behind the quoted statements
    "routeA_monod_last_margin": dict(
        min=0.0, basis="chosen", unit="1/h",
        note="route-A Monod RMS minus the largest of the other three two-parameter laws:"
             " Monod is last of the four exactly when this is positive"),
    "routeA_F_imputed_below_crit": dict(
        min=0.0, basis="chosen",
        note="F(1,4) critical value minus route A's F with row 4 imputed: the imputation"
             " flips the verdict exactly when this is positive"),
    "boot_control_size_offset_se": dict(
        max=3.0, basis="chosen",
        note="distance of the control's size from 0.05 in binomial standard errors of a"
             " nominal test at NBOOT draws: 'at nominal' means within three"),
    "boot_control_scale_free_se": dict(
        max=3.0, basis="chosen",
        note="the control's size at route A's scatter minus at route B's, in binomial"
             " standard errors: 'scale-free' means within three"),
    "boot_edge_fraction_main": dict(
        max=2.5e-3, basis="chosen",
        note="fraction of the 2 NBOOT draws at route B's scatter whose scan optimum is on"
             " a grid edge; 0 to 3 draws over twelve seeds (SEED + 1000 to SEED + 1011),"
             " the bound is ten draws, about three times the largest"),
    "blind_yield_gap_change_over_move": dict(
        max=1e-3, basis="chosen",
        note="change of the LSODA-versus-pymrm gap under a 30 % yield error, over the"
             " prediction's own move: 'blind' means the gap moves by a small part of it"),
    "marcher_invariant_mu_tripled_relerr": dict(
        max=1e-13, basis="identity",
        note="S + X/phi after tripling mu_m, relative: backward Euler conserves the"
             " linear invariant of the source (-m/phi, m) for any growth law; a"
             " different yield in the source from the one in the invariant breaks it"),
    "Sf_rootfind_over_guarantee": dict(
        max=1.0, basis="tolerance", structural=True,
        note="a solver-contract bound: the S_f crossing's distance from S_f = a over"
             " find_roots' guaranteed accuracy xtol + rtol|x|, which no physics row can move"),
}
if not set(PROVED_ZERO) <= set(B):
    raise RuntimeError(f"proved-zero quantities without a bound: {set(PROVED_ZERO) - set(B)}")
# Value claims whose rtol is tightened so that their break row fails the claim's
# own check; each stays admissible (rtol*|x0| at least 100x its spread across
# BLAS kernels), and the row's relative move is in the comment.
RTOL = {"balance_R2_six_rows": 1e-3,          # flagged row kept: 1.5e-2
        "balance_R2_all_rows": 1e-4,          # arithmetic-mean biomass: 3.3e-4
        "balance_phi": 1e-3,                  # drop row 1 instead: 1.4e-2
        "balance_a": 1e-3,                    # drop row 1 instead: 2.5e-2
        "balance_row4_residual": 1e-3,        # arithmetic-mean biomass: 1.1e-2
        "rmsA_zero_order": 1e-3,              # forward-difference mu: 8.3e-3
        "rmsB_zero_order": 1e-3,              # forward-difference mu: 8.3e-3
        "rmsB_Monod": 1e-2,                   # all-seven balance line: 4.2e-2
        "rmsB_Contois": 1e-2,                 # all-seven balance line: 4.2e-2
        "rmsB_first_order": 1e-2,             # all-seven balance line: 3.9e-2
        "routeB_monod_mu_m": 1e-2,            # all-seven balance line: 2.9e-2
        "routeA_over_routeB_monod_rms": 1e-2,  # all-seven balance line: 4.3e-2
        "Sf_crossing": 5e-3,                  # route-A constants: 9.9e-3
        # bootstrap estimates: three times the relative spread over twelve
        # seeds (SEED + 1000 to SEED + 1011, 2000 draws each), measured
        # 0.033, 0.066, 0.0052, 0.0050, 0.045 and 0.056 in this order; the
        # control size is a frequency near 0.05 of 2000 draws, whose binomial
        # standard error is about 0.1 relative, so it takes three of those
        "boot_size_forward_truth": 0.10, "boot_size_control_truth": 0.30,
        "boot_median_n_forward_truth": 0.02, "boot_median_n_control_truth": 0.02,
        "boot_p_of_observed_F": 0.15, "boot_p_of_observed_n": 0.20,
        # grid-only, route-A-scatter and sampling-error numbers: three times the
        # relative spread over the same twelve seeds, measured 0.029, 0 (a grid
        # node), 0.059, 0.10, 0 (a grid node), 0.047, 0.070; sampling errors
        # 0.0077, 0.032, 0.10, 0.12, 0.017, 0.024, and 0.016 and 0.033 propagated
        # from the route-A sizes
        "boot_gridonly_size_forward": 0.10, "boot_gridonly_median_n_forward": 0.02,
        "boot_gridonly_p_of_observed_F": 0.18, "boot_gridonly_size_control": 0.30,
        "boot_gridonly_median_n_control": 0.02,
        "boot_size_forward_at_routeA_sigma": 0.15, "boot_size_control_at_routeA_sigma": 0.20,
        "boot_size_forward_se": 0.025, "boot_size_control_se": 0.10,
        "boot_median_n_forward_se": 0.30, "boot_median_n_control_se": 0.35,
        "boot_p_of_observed_F_se": 0.05, "boot_p_of_observed_n_se": 0.075,
        "boot_size_forward_at_routeA_sigma_se": 0.05, "boot_size_control_at_routeA_sigma_se": 0.10,
        # edge draws at route A's scatter: 28 to 54 of 4000 over the twelve
        # seeds, relative spread 0.19
        "boot_edge_fraction_routeA_sigma": 0.60,
        # deterministic quoted numbers; the break row's relative move in the comment
        "routeB_monod_rms_drop": 1e-2,           # all-seven balance line: 5.6e-2
        "routeB_blackman_worst_margin": 0.1,     # all-seven balance line: 0.45
        "rmsA_null_over_monod": 1e-2,            # flagged cell imputed: 2.1e-2
        "rmsA_null_over_blackman": 1e-2,         # flagged cell imputed: 0.42
        "balance_R2_loo_min_other": 1e-4,        # arithmetic-mean biomass: 5.5e-4
        "balance_R2_loo_max_other": 1e-4,        # arithmetic-mean biomass: 3.7e-4
        "balance_row4_imputed_CA": 1e-3, "routeA_F_imputed": 1e-2,
        "rmsA_imputed_Monod": 1e-3, "rmsA_imputed_Blackman": 1e-3,
        "eoi_row3": 1e-2, "eoi_row4": 1e-2, "eoi_row5": 1e-2,   # all-seven line: 3e-2 to 0.19
        "eoi_row3_relgap": 1e-2, "eoi_row5_relgap": 1e-2,
        "rawlings_Dc_closed_form": 1e-3}
EXACT = {"flagged_row"}
for k_, v in M.items():
    if k_ in EXACT:
        continue
    if k_ in B:
        claims.bound(k_, v, **B[k_])
    else:
        claims.value(k_, v, rtol=RTOL.get(k_, 0.05))
# The flagged row is decided by the leave-one-out sweep; its margin is the R^2
# gap from the row whose removal helps second most.
claims.exact("flagged_row", int(M["flagged_row"]),
             margin=R2_BAL - M["balance_R2_loo_max_other"], min_margin=1e-3,
             note="the row (1-based) whose removal raises the balance line's R^2 most;"
                  " margin: R^2 over the next best removal")
# The one count on the page: N* exists only for an alternative whose separation
# exceeds the 1e-9 1/h the degenerate two are bounded by, and the margin to that
# threshold is recorded so that the set of N* claims can never be decided by noise.
claims.exact("n_separable_alternatives", sum(1 for d in SEP.values() if d > SEP_BOUND),
             # the margin of the separable side, in decades above the bound; the
             # degenerate side is held below the same bound by the sep_moser and
             # sep_contois bound claims, whose round-off values differ by machine
             margin=min(np.log10(d/SEP_BOUND) for d in SEP.values() if d > SEP_BOUND),
             min_margin=1.0, structural=True,
             note="alternatives with a separation above the SEP_BOUND 1/h bound, so"
                  " that N* is defined: Blackman and Tessier; margin in decades to"
                  " that bound on the separable side, the degenerate side held by"
                  " its two bound claims; decided by the identities, which no row can move")
if claims["n_separable_alternatives"] != sum(1 for k in M if k.startswith("Nstar_")):
    raise RuntimeError("the N* claims do not match the separable alternatives")

# -- the break table, wired row by row to the claim it targets: a row breaks a
#    claim when the claim's own check (rtol or bound) fails on the injected value
for r in BREAK:
    claims.breaks(r["row"], r["metric"], r["broken"])
_kinds = {kd: [k_ for k_, c in claims.metrics.items() if c["kind"] == kd]
          for kd in ("value", "bound", "exact")}
print(f"{len(claims.metrics)} claims registered: {len(_kinds['value'])} values,"
      f" {len(_kinds['bound'])} bounds, {len(_kinds['exact'])} exact.  The bounds:")
for k_ in _kinds["bound"]:
    c = claims.metrics[k_]
    lo, hi = c.get("min"), c.get("max")
    rng = (f"[{lo:g}, {hi:g}]" if lo is not None and hi is not None
           else f">= {lo:g}" if lo is not None else f"<= {hi:g}")
    print(f"  {k_:36s} {rng:14s} [{c['basis']}]  (here {c['value']:.3g})")
# claims.report() runs at the end of the prose cell below, after every number
# the prose quotes has been checked against its claim with claims.says.

# %%
# Every number quoted in prose is checked here against its registered claim by
# claims.says; README.md, meta.yaml and models.yaml tag theirs with claim markers,
# which scripts/check_claims.py checks against agreement.json.
# BOUNDS, NOT PINS. Every quantity in the proved-exactly-zero class is an exit
# tolerance or a round-off residual, so it is a registered bound and is quoted
# as that bound in the prose, in meta.yaml and in models.yaml - never to two
# significant figures, whose last digit differs between machines.  The exact
# values are PRINTED in the cells above, beside their bound.  The limits join
# the sweep's pool below.  One bound is NOT a claim limit and is not in this
# list: the crossing's root-find is checked against the root finder's own xtol +
# rtol|x|, computed where the root is found, because a round number there
# would be a bound on an accident.
BOUNDS = sorted({c["max"] for c in claims.metrics.values()
                 if c["kind"] == "bound" and c.get("max")})

# -- numbers that ARE registered claims, as written in the prose: checked by
#    claims.says, within one unit of the last written place, or, for a bound
#    phrase such as "below 1e-9", at least as loose as the claim's own bound
SAYS = {
    "balance_R2_six_rows": ["0.99996"], "balance_R2_all_rows": ["0.98505"],
    "balance_row4_residual": ["-17.2"], "balance_phi": ["0.3174"],
    "balance_a": ["197.12"], "balance_worst_other_residual": ["0.65"],
    "row4_end_of_interval_relgap": ["3.1 %"],
    "rmsA_Monod": ["0.05469"], "rmsA_Blackman": ["0.02677"],
    "rmsA_Tessier": ["0.03871"], "rmsA_Contois": ["0.05153"],
    "rmsB_Monod": ["0.03147"], "rmsB_Blackman": ["0.03955"],
    "rmsB_Tessier": ["0.02631"],
    "routeA_n_moser": ["1.93"], "routeB_n_moser": ["1.51"],
    "routeA_F_moser_vs_monod": ["12.3"], "routeB_F_moser_vs_monod": ["1.44"],
    "F_crit_95": ["7.71"],
    "routeB_monod_mu_m": ["0.8809"], "routeB_monod_K_S": ["20.45"],
    "boot_size_forward_truth": ["0.37", "36.55 %"],
    "boot_size_control_truth": ["0.0495"],
    "boot_median_n_forward_truth": ["1.61"], "boot_median_n_control_truth": ["1.01"],
    "boot_p_of_observed_F": ["0.23"], "boot_p_of_observed_n": ["0.13"],
    "Nstar_tessier": ["10.3"], "Nstar_blackman": ["1.30"],
    "Dc_ratio_at_Sf20": ["1.83"], "Sf_for_10pct_gap": ["95.5"],
    "Sf_crossing": ["197.1"],
    "sym_blackman_gap_lower_over_mum": ["0.0858"],
    "bias_over_routeA_rms": ["1.21"], "bias_max_abs_excl_bad_row": ["0.0663"],
    "max_swing_over_CAbar": ["1.89"],
    "pfr_observed_order": ["1.00"],
    "pymrm_batch_observed_order": ["1.00"],
    "pfr_richardson_relerr": ["below 1e-4"],
    "rawlings_Dc_relerr": ["2.1e-4"],
    "rawlings_10_19_max_relerr": ["below 1e-9"],
    "Dc_monod_transient_relerr": ["below 1 %"],
    "Dc_contois_transient_relerr": ["below 1 %"],
    "pymrm_refit_mum_relerr": ["below 1e-4"],
    "pymrm_refit_KS_relerr": ["below 1e-4"],
    "pymrm_batch_conservation_relerr": ["below 1e-13"],
    # the proved-zero class, quoted only as the bounds they are checked against
    "contois_map_mum_relerr": [f"below {MAP_BOUND:g}"],
    "contois_map_KS_relerr": [f"below {MAP_BOUND:g}"],
    "Sf_crossing_vs_a_relerr": [f"below {MAP_BOUND:g}"],
    "contois_minus_monod_rms": [f"below {RMS_GAP_BOUND:g}"],
    "sep_moser": [f"below {SEP_BOUND:g}"], "sep_contois": [f"below {SEP_BOUND:g}"],
    # the quoted numbers registered next to the break table
    "routeB_monod_rms_drop": ["42 %"], "routeB_blackman_worst_margin": ["0.0081"],
    "routeA_monod_last_margin": ["above 0"], "routeA_F_imputed_below_crit": ["above 0"],
    "rmsA_null_over_monod": ["3.68"], "rmsA_null_over_blackman": ["7.53"],
    "balance_R2_loo_min_other": ["0.9757", "0.975738"],
    "balance_R2_loo_max_other": ["0.9848", "0.984806"],
    "blackman_Nstar_shortfall": ["30 %"], "rawlings_Dc_closed_form": ["0.8333"],
    "boot_gridonly_size_forward": ["0.308"], "boot_gridonly_size_control": ["2.1 %"],
    "boot_gridonly_median_n_forward": ["1.60"], "boot_gridonly_median_n_control": ["0.95"],
    "boot_gridonly_p_of_observed_F": ["0.18"],
    "boot_control_size_offset_se": ["below 3"], "boot_control_scale_free_se": ["below 3"],
    "boot_size_forward_at_routeA_sigma": ["0.2355"], "boot_size_control_at_routeA_sigma": ["0.0470"],
    "boot_size_forward_se": ["1.1 %", "0.0108"], "boot_size_control_se": ["0.5 %", "0.0049"],
    "boot_median_n_forward_se": ["0.0074"], "boot_median_n_control_se": ["0.0049"],
    "boot_p_of_observed_F_se": ["0.0094"], "boot_p_of_observed_n_se": ["0.0075"],
    "boot_size_forward_at_routeA_sigma_se": ["0.0095"],
    "boot_size_control_at_routeA_sigma_se": ["0.0047"],
    "boot_edge_fraction_main": ["below 2.5e-3"], "boot_edge_fraction_routeA_sigma": ["0.93 %"],
    "blind_yield_gap_change_over_move": ["below 1e-3"],
    "marcher_invariant_mu_tripled_relerr": ["below 1e-13"],
    "Sf_rootfind_over_guarantee": ["below 1"],
    "balance_row4_imputed_CA": ["60.18"], "routeA_F_imputed": ["5.45"],
    "rmsA_imputed_Monod": ["0.04687"], "rmsA_imputed_Blackman": ["0.04065"],
    "routeA_imputed_blackman_margin": ["0.0062"],
    "eoi_row3": ["74.89"], "eoi_row4": ["44.33"], "eoi_row5": ["13.46"],
    "eoi_row3_relgap": ["16.8 %"], "eoi_row5_relgap": ["53.6 %"],
    "flagged_row": ["4"],
}
for k_, lits in SAYS.items():
    for lit in lits:
        claims.says(lit, k_)
N_SAYS = sum(len(v) for v in SAYS.values())

print(f"{N_SAYS} prose literals checked against their registered claims by claims.says.")
claims.report();   # the semicolon keeps the notebook from echoing the payload
# %% [markdown]
# ## What pymrm adds
#
# **Honestly: the central argument is algebra and a seven-point fit, and scipy does
# most of it. Two things pymrm does here are load-bearing, and one is not.**
#
# Load-bearing:
#
# - **The second, independent computation of the fit.** `BatchMarch` is a
#   backward-Euler marcher built from `NumJac((1, 2))` and `newton`, and it
#   re-derives $\mu_m$ and $K_S$ through a different discretisation, a different
#   Jacobian and a different linear solve, agreeing with the LSODA fit to below
#   1e-4 in both constants. A break row confirms the check has teeth by running it unconverged. Be
#   precise about what it buys: it catches a bug in *either integrator*, and it is blind to everything the two
#   share - a wrong yield, a wrong intercept or a mis-set interval moves both
#   predictions together, which the Validation section measures. Carrying $S$
#   as an unknown instead of substituting the yield relation is worth having for
#   its own sake, but $S + X/\varphi$ is then a **tautology** of the scheme, not a
#   test of it: the source is $(-m/\varphi,\,m)$ built from one $m$, so backward
#   Euler conserves it for any growth law, and tripling $\mu_m$ moves it by
#   nothing at all.
# - **Turning a caption into a test.** `PlugFlowFermenter` makes Levenspiel's
#   *"batch **or plug flow**"* into a number: first-order upwind convergence at
#   observed order 1.00 from successive differences, and a Richardson value
#   within 1e-4 of the initial-value solution. The break table then prices the
#   usual mistakes on it, `nu=1` instead of `nu=0` and a flipped substrate sign,
#   both while the profile still looks smooth.
#
# Not load-bearing, and the page says so: a well-mixed batch and a chemostat are
# 0-D, and `solve_ivp` integrates them - as it does, throughout, as the reference.
#
# **What travels.** The growth law is written as a pointwise source over the last
# axis, which is the shape the rest of the gallery uses, so
# `PlugFlowFermenter.source` transplants unchanged into a pellet (`B1.1`), an
# immobilised-enzyme particle (`J4.7`) or a dispersed tubular reactor (`A2.1`,
# `A2.6`): the ingredients gain `construct_grad`, and the rate law is the same
# object.

# %% [markdown]
# ## Reuse
#
# **If you have batch growth data and want $\mu_m$ and $K_S$:** do not compare
# $\ln(C_C^{\rm end}/C_C^{\rm start})/\Delta t$ with $\mu$ evaluated at the
# interval-average substrate concentration unless the substrate barely moves across
# the interval. On these seven intervals that shortcut biases the rate by up to
# 0.0663 h$^{-1}$ - **1.21 times the residual scatter it is then judged against** -
# always in the direction that makes a saturating law look too gradual near
# exhaustion. Integrate the law across the interval instead; constant yield makes
# that a single ODE and it costs nothing. `fit_B` on this page is the whole of it.
#
# **If you are choosing between growth laws:** ask first whether your experiment
# *can* tell them apart. Two of the four alternatives Rawlings & Ekerdt print
# beside Monod cannot be separated from it by **any** constant-yield batch run.
# Moser contains Monod at $n = 1$. **Contois *is* Monod**, exactly, with
# $\mu_m' = \mu_m/(1-K_{sx}\varphi)$ and
# $K_S' = K_{sx}\varphi a/(1-K_{sx}\varphi)$ - so a batch fit that "prefers" one
# over the other is reporting its own optimiser and its own reduction. Of the two
# that are separable at this design, Blackman needs 1.30 repeats of the
# seven-interval run and Tessier 10.3 - and read those as ~50 %-power thresholds,
# so one clean run is *within 30 % of* enough for Blackman rather than enough.
#
# **If you need to separate Contois from Monod:** run a chemostat and sweep the
# *feed* concentration, not the dilution rate alone. Monod washes out at
# $D_c = \mu_m S_f/(K_S+S_f)$ and Contois at $D_c = \mu_m$ regardless of $S_f$, so
# the two predictions cross at exactly $S_f = a$ - the substrate intercept of your
# own batch run's material balance - and separate on either side. At the constants
# fitted here that is $S_f = 197.1$; the gap exceeds 10 % below $S_f = 95.5$ and
# reaches a factor 1.83 at $S_f = 20$.
#
# **If you are checking a printed table:** the constant-yield material balance is
# free, needs no kinetics and is sharp. Six of these seven rows sit on it at
# $R^2 = 0.99996$; the seventh is the whole reason the all-seven fit reads 0.98505.
# A leave-one-out sweep says *which* row, not merely that one is off - and that
# distinction matters, because "the table has scatter" and "one cell is wrong" lead
# to different analyses. Then price the cell before quoting anything that depends
# on it: on this table, putting that one cell at its balance-implied value takes
# route A's F from 12.27 to 5.45 and flips the only significance verdict the page
# could have reported.
#
# **If you want the neighbouring physics:** substrate inhibition - Froment's
# eq. (1.5.2-4) and Rawlings' $\mu_m S/(K_s+S+K_1S^2)$, the form whose rate has a
# *maximum* where Monod has an asymptote - is case **J4.2**, and gets no
# quantitative treatment here. Williams' structured model, eqs. (1.5.2-5) to
# (1.5.2-8) of the same section, is **J4.9**. The enzyme rate law Froment says
# Monod-kinetics are *"shaped after"* is **J4.6**, and the isotherm this page proves
# Monod identical to is **J1.1**.
#
# **What NOT to reuse.** $K_S$ from this page is in the units of a column whose
# units the source does not print, for an organism the source does not name, on a
# substrate identified only as lactose. It is a number about a table, not about a
# culture. And the seven rows are attributed by Levenspiel to Monod (1958) p. 74;
# neither that book nor the 1949 paper this case catalogues was consulted, so do
# not cite them as either.
