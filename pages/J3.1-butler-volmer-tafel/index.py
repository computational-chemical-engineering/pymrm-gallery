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
# title: "Butler–Volmer kinetics, and where its two limits stop being safe"
# description: "The interfacial charge-transfer law in the form a reactor model uses it, read off Doyle, Fuller and Newman (1993), the exchange current density derived from their own equations, and the overpotential at which the Tafel and linear approximations go wrong by 1, 5 and 10 %."
# categories: [sec:J, struct:S1, tier:T0, data:tier6, phase:liquid-solid]
# date: 2026-07-31
# ---
#
# # Butler–Volmer kinetics, and where its two limits stop being safe
#
# **Catalog ID:** `J3.1` · **Structures:** `S1` (algebraic constitutive law) ·
# **Tier:** T0
#
# Almost every electrochemical reactor model contains one nonlinear algebraic
# equation that turns a potential difference into a current density. Almost every
# one of them is then simplified, to a linear resistance, or to a single
# exponential, and the simplification is usually made without checking. This page
# does the checking, and puts numbers on it.

# %% [markdown]
# ## Background
#
# ### What this page is a page *of*, and what it is not
#
# The Butler–Volmer equation is conventionally attributed to Tafel (1905), to
# J. A. V. Butler (1924, 1932) and to Erdey-Grúz and Volmer (1930). **None of
# those papers was consulted for this page**, none is available here, and their
# bibliographic details are recorded from secondary attribution only, without a
# DOI, precisely because an unverified citation is worse than none. Nothing on
# this page is a restatement of a source that was not read.
#
# What *was* read is **Doyle, Fuller and Newman (1993)**, *J. Electrochem. Soc.*
# **140**(6) 1526–1533, which prints the law in exactly the form a reactor model
# uses it, twice: once as the general two-exponential kinetics of a flat electrode
# (their **Eq. 6**, with the exchange current density of **Eq. 8**), and once as a
# mass-action expression for an insertion reaction (their **Eq. 17**, with the
# exchange current density of **Eq. 30**). That paper is the source of the
# gallery's [`J3.4`](../J3.4-doyle-fuller-newman/) page, so its parameter set has
# already been extracted from 600 dpi page renders and reviewed, and its
# open-circuit expression has already been validated against the authors' own
# plotted curve to 3.3 mV. Everything algebraic here starts from those printed
# equations.
#
# A second paper on disk, **Marquis, Sulzer, Timms, Please and Chapman (2019)**
# (*J. Electrochem. Soc.* **166**(15) A3693–A3706, the source of
# [`J3.5`](../J3.5-single-particle-model/)), prints the same kinetics in its
# hyperbolic-sine form and prints the closed-form inversion. That gives a second,
# independent printed witness for two of the identities below.
#
# ### Why this is a separate page from `J3.4`
#
# `J3.4` *uses* Eq. 17 inside a full pseudo-two-dimensional cell; it never looks at
# the law itself. It does not derive the exchange current density, does not ask
# what happens when the two transfer coefficients do not sum to one, and does not
# touch either limit, its `arcsinh` is the exact inversion of the symmetric law,
# not an approximation to it. This page is about the constitutive law: what its
# parameters mean, which of its features are identities rather than choices, and,
# the part that is worth more than any restatement: **the overpotential at which
# each of its two standard approximations goes wrong by a stated amount**, checked
# against the overpotentials that the published `J3.4` cell actually operates at.
#
# Two honest qualifications, neither of which changes that.
#
# * **Part of the algebra is already in the repository, unremarked.** `J3.5`'s
#   `bv_invert` is a closed-form inversion of Eq. 17 whose structure silently
#   *contains* the reduction derived below: at zero current its root is Eq. 16
#   exactly, and its exchange prefactor is Eq. 30's. `J3.5` never says so, never
#   states the $\alpha_a+\alpha_c=1$ condition it depends on, and never uses it as
#   a check, but a reader who knows that page should know the algebra is not new
#   to the gallery, only unexamined in it.
# * **One number here already exists elsewhere on a different scale.** `J3.5`'s
#   `C_r_doyle` is the same interfacial-versus-superficial ratio this page uses to
#   classify the cathode, and reports it as $\sim10^{-4}$ against the $1.0\times
#   10^{-3}$ here. The two are not in conflict; they divide by different exchange
#   current densities, and the ratio between them is computed in the Reuse section
#   below.
#
# The answer for that cell is not the one most readers would guess, and it is
# different at its two electrodes.

# %% [markdown]
# ## The published model
#
# ### The law, as printed
#
# Doyle, Fuller and Newman write the kinetics of the lithium foil as (their
# **Eq. 6**, with the surface overpotential $\eta_{s1}=\Phi_1-\Phi_2-U_1$ of
# **Eq. 7**)
#
# $$
# I = i_{o1}\left[\exp\!\left(\frac{\alpha_{a1}F\eta_{s1}}{RT}\right)
# -\exp\!\left(-\frac{\alpha_{c1}F\eta_{s1}}{RT}\right)\right],
# \tag{6}
# $$
#
# $$
# i_{o1} = F(k_{a1})^{\alpha_{c1}}(k_{c1})^{\alpha_{a1}}
# (c_{\max}-c)^{\alpha_{c1}}(c)^{\alpha_{a1}}.
# \tag{8}
# $$
#
# That is the Butler–Volmer equation. The insertion reaction at the porous cathode
# is written in mass-action form instead (their **Eq. 17**, with $U'$ of Eq. 18 and
# $\eta=\Phi_1-\Phi_2$ of Eq. 19):
#
# $$
# i = Fk_2(c_{\max}-c)^{\alpha_c}c^{\alpha_A}
# \left[c_s\exp\!\left(\frac{\alpha_aF}{RT}(\eta-U')\right)
# -(c_T-c_s)\exp\!\left(-\frac{\alpha_cF}{RT}(\eta-U')\right)\right],
# \tag{17}
# $$
#
# $$
# U' = U_2^{\theta}-U_{\rm ref}^{\theta}+\frac{RT}{F}(\beta c_s+\zeta).
# \tag{18}
# $$
#
# Both were read off 600 dpi renders of journal pages 1527 and 1528; the scan's
# text layer mangles them beyond use (Eq. 17 comes out as
# `i = Fk2(emax - c)~~ ... ~R~01 - U'))`). The exponent printed on $c$ in Eq. 17 is
# typeset `α_A`; it is $\alpha_a$, and the check in the next section proves it,
# because no other reading reproduces the paper's own Eq. 30.
#
# Separately, the paper prints the **open-circuit potential** of the insertion
# material (their **Eq. 16**)
#
# $$
# U_2 = U_2^{\theta}-U_{\rm ref}^{\theta}
# +\frac{RT}{F}\left(\ln\frac{c_T-c_s}{c_s}+\beta c_s+\zeta\right)
# = U' + \frac{RT}{F}\ln\frac{c_T-c_s}{c_s},
# \tag{16}
# $$
#
# and the **exchange current density** of the insertion reaction (their **Eq. 30**)
#
# $$
# i_{o2} = F(k_2)(c_{\max}-c)^{\alpha_{c2}}(c)^{\alpha_{a2}}
# (c_T-c_s)^{\alpha_{a2}}(c_s)^{\alpha_{c2}}.
# \tag{30}
# $$
#
# Eqs. 16, 17 and 30 are printed in three different places, eleven journal pages
# apart in the case of the last two, and are **not** derived from one another in
# the paper. They therefore constitute an internal identity, and it is the
# strongest check available here.

# %% [markdown]
# ### The identity that ties them together
#
# Put $\eta_s = \eta - U_2$, the *surface* overpotential, measured from the local
# equilibrium potential rather than from the reference, and substitute Eq. 16 into
# Eq. 17. The bracket becomes
#
# $$
# c_s\left(\tfrac{c_T-c_s}{c_s}\right)^{\alpha_a}
# e^{\alpha_a F\eta_s/RT}
# -(c_T-c_s)\left(\tfrac{c_T-c_s}{c_s}\right)^{-\alpha_c}
# e^{-\alpha_c F\eta_s/RT},
# $$
#
# whose two prefactors are $c_s^{1-\alpha_a}(c_T-c_s)^{\alpha_a}$ and
# $c_s^{\alpha_c}(c_T-c_s)^{1-\alpha_c}$. **They are equal if and only if**
# $\alpha_a+\alpha_c=1$, and when they are, Eq. 17 collapses to
#
# $$
# i = i_{o2}\left[e^{\alpha_aF\eta_s/RT}-e^{-\alpha_cF\eta_s/RT}\right],
# \qquad
# i_{o2}=Fk_2(c_{\max}-c)^{\alpha_c}c^{\alpha_a}(c_T-c_s)^{\alpha_a}c_s^{\alpha_c},
# $$
#
# which is Doyle's Eq. 30 as printed, factor by factor and exponent by exponent,
# including which of $\alpha_a$ and $\alpha_c$ sits on which concentration, the
# detail a transcription is most likely to get wrong. It is also Eq. 6 with a
# different exchange current density, which is what makes Eq. 6 and Eq. 17 the same
# law.
#
# Read as a *check*, though, the agreement is narrower than the derivation. The
# $(c_{\max}-c)^{\alpha_c}c^{\alpha_a}$ prefactor is the same code expression on
# both sides and cancels in the ratio, so the numerical comparison tests only the
# $(c_T-c_s)^{\alpha_a}c_s^{\alpha_c}$ half against Eq. 17's bracket. And **at
# Doyle's own $\alpha_a=\alpha_c=\tfrac12$ every exponent swap is a no-op**, so
# none of it bites at the published parameters. What gives the check its power is
# the asymmetric rows of the sweep below. Both statements are measured, not
# asserted, in the injected-defect table that follows Check 2.
#
# Three consequences, all checked numerically below:
#
# 1. **The zero-current potential of Eq. 17 equals Eq. 16 only when
#    $\alpha_a+\alpha_c=1$.** Otherwise the kinetics and the thermodynamics of the
#    same paper disagree, by $\frac{RT}{F}\ln\frac{c_T-c_s}{c_s}
#    \left(\frac{1}{\alpha_a+\alpha_c}-1\right)$. This is a real constraint on the
#    model, not an accident of notation, and it is why $\alpha_a+\alpha_c$ appears
#    rather than $\alpha_a$ alone in the paper's own Eq. 29.
# 2. **An exchange current density exists only under the same condition.** For
#    $\alpha_a+\alpha_c\ne1$ the two branches have different prefactors and no
#    single $i_0$ can be factored out: Eq. 30 would then not be well defined.
# 3. Doyle's Table II has $\alpha_a=\alpha_c=0.5$, so the condition holds for the
#    published cell and Eq. 30 is what the paper says it is.
#
# ### The two limits
#
# With $\alpha_a+\alpha_c=1$ (write $n=\alpha_a+\alpha_c$ for the general case):
#
# **Linear.** Expanding about $\eta_s=0$,
# $i \to i_0\,\frac{nF}{RT}\,\eta_s$, i.e. a pure resistance
#
# $$
# R_{\rm ct} = \left(\frac{\partial i}{\partial\eta_s}\right)^{-1}_{\eta_s=0}
# = \frac{RT}{nFi_0}.
# $$
#
# The familiar textbook form $RT/(nFi_0)$ is therefore an identity only when the
# transfer coefficients sum to the electron number. That is an assumption, not a
# theorem, and the page states it rather than assuming it.
#
# **Tafel.** For $\eta_s$ large and positive the cathodic term dies and
# $i\to i_0e^{\alpha_aF\eta_s/RT}$, so a plot of $\eta_s$ against $\log_{10}i$ is a
# straight line of slope
#
# $$
# b = \frac{\ln 10\;RT}{\alpha_a F}\quad(\text{anodic}),\qquad
# b = -\frac{\ln 10\;RT}{\alpha_c F}\quad(\text{cathodic}).
# $$
#
# **Symmetric case.** When $\alpha_a=\alpha_c=\alpha$ the law is
# $i = 2i_0\sinh(\alpha F\eta_s/RT)$, which inverts in closed form,
#
# $$
# \eta_s = \frac{RT}{\alpha F}\,\operatorname{arcsinh}\!\frac{i}{2i_0}.
# $$
#
# Marquis et al. print both of these (their **Eqs. 1g–1h** for the hyperbolic-sine
# law and its exchange current density, and their **Eq. 18** for the inversion, in
# dimensionless form with $\alpha=\tfrac12$), and the `J3.4` page already uses the
# inversion for its lithium foil. It is checked below by feeding the closed form
# forward through the printed Eq. 17.

# %% [markdown]
# ## Parameters and assumptions
#
# **Assumptions built into the law itself**, all of them stated by Doyle et al.:
# a single elementary charge-transfer step; potential-independent transfer
# coefficients; concentrations at the reacting interface, not in the bulk (any
# mass-transfer resistance is outside the law); an ideal solution, "because
# activity coefficient data have not been reported"; isothermal operation; and no
# surface film. The insertion form additionally assumes site-limited kinetics on
# both sides of the interface, which is where the $(c_{\max}-c)$ and $(c_T-c_s)$
# factors come from.
#
# **Numerical inputs** are Doyle's Table II, Appendix A and Eq. 20, loaded
# **cross-page** from the reviewed `J3.4` dataset rather than re-extracted. The
# relevant ones are $k_2=10^{-10}$ m⁴/(mol s), $\alpha_a=\alpha_c=0.5$,
# $c_T=29\,000$ mol/m³, $c_{\max}=3920$ mol/m³, $c_0=1000$ mol/m³, $T=373.15$ K,
# and $i_{o1}=12.6$ A/m² for the lithium foil. Note the temperature: this cell runs
# at 100 °C, so $RT/F=32.15$ mV rather than the 25.69 mV of room temperature, and
# every threshold below scales with it. Both are reported.
#
# **One thing this page deliberately does not decide.** `J3.4` records an open
# question about the factor $\nu=2$ in the concentration term of Doyle's Eq. 5:
# thermodynamics requires it for a 1:1 salt, the printed equation omits it, and the
# paper's Figure 2 cannot separate the two. That question is about the electrolyte
# potential, not about the kinetics. Where the cell is solved below it is solved
# **both ways**, and the kinetic conclusions are reported for both, so nothing here
# depends on the answer.
#
# **Provenance tier 6.** Nothing on this page is compared with a measurement.
# Doyle's Figure 2 contains no experimental points at all, confirmed by the
# gallery maintainer on review of the original, and the checks here are algebraic
# identities plus the paper's own printed numbers. This page is *reproduced* and
# *internally verified*; it is **not** validated against experiment, and the two
# must not be blurred.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.20.4"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "J3.1-butler-volmer-tafel"  # equals the page directory

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
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Markdown
from scipy.sparse.linalg import MatrixRankWarning
from pymrm.checks import observed_orders
from gallery_utils import load_data, load_meta, cite_data, load_claims, Claims
from model import (R_GAS, F_CONST, Kinetics, Cell, kappa_from_delta, first_lobe_peak,
                   one_root, u_at_level)

warnings.filterwarnings("ignore", category=MatrixRankWarning)
PAGE_J34 = "J3.4-doyle-fuller-newman"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])
J34 = load_claims(PAGE_J34)          # J3.4's committed claims, quoted beside this page's

# %% [markdown]
# ## The data
#
# **This page ships no dataset of its own, and digitises no figure.** It has
# nothing to digitise: the result is an algebraic law and its limits, and the only
# numbers that enter come from Doyle's tables and running text.
#
# Two files are loaded cross-page from `J3.4`, whose provenance sidecars carry the
# full extraction record:
#
# * `doyle-fuller-newman-1993-parameters.csv`: Tables I and II, Appendix A and the
#   Eq. 20 constants, every value read off a 600 dpi render;
# * `doyle-fuller-newman-1993-stated-results.csv`: the scalar results Doyle prints
#   in running text, of which the printed $\delta$ and $\nu$ are used here as a
#   printed-number check on Eq. 30.
#
# What `J3.4` has established about these rows, and whether it matters here:
#
# * **The conductivity is not printed.** `J3.4` reconstructs it by inverting
#   Eqs. 28 and 29, and the two routes disagree (its claim
#   `kappa_route_spread_pct`, printed below). The kinetic law and its thresholds
#   do not contain it; the cell does, through the Eq. 28 route that `J3.4` uses,
#   and Check 3 below shows the disagreement is the same number as this page's
#   $\nu^2/\delta$ test, not a second one.
# * **The Bruggeman exponent of the salt flux** is read as $D\varepsilon^{0.5}$
#   on `J3.4`, which this page's cell inherits; it enters the cell's transport,
#   not the kinetics.
# * **The factor $\nu=2$ in Eq. 5** is left open on `J3.4`. The cell is run both
#   ways below and the kinetic conclusions do not depend on it.

# %%
par = load_data("doyle-fuller-newman-1993-parameters.csv", page=PAGE_J34)
stated = load_data("doyle-fuller-newman-1993-stated-results.csv", page=PAGE_J34)
par_meta = load_meta("doyle-fuller-newman-1993-parameters.csv", page=PAGE_J34)

P = {r.symbol: float(r.value) for r in par.itertuples()}
S = {r.quantity: float(r.value) for r in stated.itertuples()}
kin = Kinetics(P)                              # the law at Doyle's parameters
T = P["T"]
RTF = kin.rtf                                  # RT/F at Doyle's 100 C
f_RT = 1.0 / RTF
RTF_298 = R_GAS * 298.15 / F_CONST             # for the room-temperature column
C_T, C_MAX = P["c_T"], P["c_max"]
eps_c, dc, ds = P["epsilon"], P["delta_c"], P["delta_s"]
alpha, sigma = P["alpha_a"], P["sigma"]
a_spec = 3.0 * (1 - eps_c) / P["R_s"]          # interfacial area per unit volume
T_FULL = (1 - eps_c) * C_T * dc * F_CONST      # discharge time x current density

print(f"RT/F = {1e3*RTF:.3f} mV at Doyle's T = {T:.2f} K   "
      f"({1e3*RTF_298:.3f} mV at 298.15 K)")
print(f"alpha_a = {P['alpha_a']}, alpha_c = {P['alpha_c']}, "
      f"sum = {P['alpha_a'] + P['alpha_c']}")
print(f"k_2 = {P['k_2']:.3g} m4/(mol s), c_T = {C_T:.0f}, c_max = {C_MAX:.0f}, "
      f"c_0 = {P['c_0']:.0f} mol/m3, i_o1 = {P['i_0_1']} A/m2")
print(f"printed groups: delta = {S['delta']}, nu = {S['nu']:.0f}")
print(f"specific area a = {a_spec:.3g} 1/m")
print(f"J3.4's two routes to the unprinted conductivity differ by "
      f"{J34['kappa_route_spread_pct']:.1f} % (its claim kappa_route_spread_pct)")
print()
print(cite_data(par_meta))


# %% [markdown]
# ## PyMRM implementation
#
# ### The law is pointwise; the question is not
#
# The kinetic expression itself is algebra, it needs no grid, no operators and no
# solver, and manufacturing a spatial problem for it would be dishonest. So the
# first two thirds of this page are the law and its limits, evaluated directly.
#
# But *"is the approximation I am about to make safe?"* is not a pointwise
# question. The overpotential a reaction actually runs at is set by the whole
# reactor: by how the current distributes itself across a porous electrode, by how
# the exchange current density moves as the surface composition changes, and by
# what the transport does to it at the end of a discharge. That needs the reactor,
# and here the reactor is the published `J3.4` cell, a pymrm two-phase porous
# electrode. It is used as an **operating envelope**, not as a result to
# reproduce: the quantity taken from it is the distribution of $\eta_s$, which
# `J3.4` never reports.
#
# The `Cell` class in this page's `model.py` is `J3.4`'s cell with the kinetic
# law of the cathode selectable. Its time step carries a `dt_scale` that
# multiplies both the initial step and the ceiling, so that the temporal error
# can be refined at all; `dt_scale = 1` is the production schedule. The grid, the
# operators (`construct_grad`, `construct_div` with `nu=0` for the slab), the
# boundary conditions on the outward normal, the harmonic face means at the
# porosity jump and the salt balance in flux form are `J3.4`'s. Each implicit
# Euler step is one `newton` solve on a `NumJac((N, 3))` Jacobian whose stencil
# couples the three fields of a cell in full and neighbouring cells through the
# operators; a step whose `newton` does not report success is retried with a
# shorter step. The cell potential at a fixed utilisation is computed there, by
# one implicit-Euler step from the last state before the time at which Faraday's
# law puts that utilisation, not interpolated between steps. The four settings
# are
#
# | `kinetics` | what is solved |
# |---|---|
# | `eq17` | Doyle's Eq. 17 exactly as printed |
# | `bv` | the standard form $i_0[e^{\alpha_a f\eta_s}-e^{-\alpha_c f\eta_s}]$ with $i_0$ from Eq. 30 |
# | `linear` | $i_0\,nf\eta_s$ |
# | `tafel` | the cathodic branch alone, $-i_0e^{-\alpha_c f\eta_s}$ |
#
# `eq17` and `bv` must produce *identical* cells; that they do is one of the
# checks.

# %%
# The law at Doyle's parameters, from model.py: Eq. 16 (U_ocp), Eq. 18 (U_prime),
# Eq. 17 as printed (i_eq17, eta = Phi_1 - Phi_2 of Eq. 19), Eq. 30 (i0_eq30), the
# standard form (i_bv), its limits, R_ct and the validity thresholds.
U_ocp, U_prime, i_eq17, i0_eq30 = kin.U_ocp, kin.U_prime, kin.i_eq17, kin.i0_eq30
i_bv, i_linear, i_tafel, R_ct = kin.i_bv, kin.i_linear, kin.i_tafel, kin.R_ct
eta_equilibrium, i0_general = kin.eta_equilibrium, kin.i0_general
eta_tafel_threshold = kin.eta_tafel_threshold
eta_tafel_threshold_numeric = kin.eta_tafel_threshold_numeric
eta_linear_threshold = kin.eta_linear_threshold
print(f"i_o2 from Eq. 30 at c_0 and half utilisation: {float(i0_eq30(P['c_0'], 0.5 * C_T)):.3f} A/m2")


# %% [markdown]
# ### Check 1: Eq. 17, Eq. 16 and Eq. 30 are one equation
#
# The transcription is tested by feeding Eq. 17 as printed into the reduction
# derived above and comparing with Eq. 30 as printed. Nothing is fitted and the two
# readings are independent: Eq. 17 is on journal page 1528, Eq. 30 on page 1530.
#
# The sweep covers the whole physical range of both concentrations and four decades
# of current, and the transfer coefficients are swept **off** the symmetric value
# as well, so that the $\alpha_a+\alpha_c=1$ condition is tested rather than
# assumed.

# %%
rng = np.random.default_rng(20260731)
n_pts = 4000
c_r = rng.uniform(50.0, C_MAX - 50.0, n_pts)
cs_r = rng.uniform(0.005 * C_T, 0.995 * C_T, n_pts)
eta_r = rng.uniform(-0.30, 0.30, n_pts)          # surface overpotential, V

SUM1_PAIRS = [(0.5, 0.5), (0.3, 0.7), (0.7, 0.3), (0.15, 0.85)]
NOT1_PAIRS = [(0.5, 0.4), (0.5, 0.6), (0.25, 0.25), (1.0, 1.0)]
U2_r = U_ocp(cs_r / C_T)


def d_eq17_vs_eq30(aa, ac, e17pre=None, e17brk=None, e30pre=None, e30tail=None):
    """Check 1: Eq. 17 with U_2 from Eq. 16, against Eq. 30 in Eq. 6, max relative gap.
    The keyword arguments mis-set one exponent assignment each (defects)."""
    lhs = i_eq17(U2_r + eta_r, c_r, cs_r, aa, ac, pre=e17pre, brk=e17brk)
    rhs = i_bv(eta_r, i0_eq30(c_r, cs_r, aa, ac, pre=e30pre, tail=e30tail), aa, ac)
    sc = np.maximum(np.abs(lhs), np.abs(rhs))
    return float(np.max(np.abs(lhs - rhs) / np.where(sc > 0, sc, 1.0)))


_ROOTS = {}


def zero_current_roots(aa, ac, n=400):
    """eta_s at which Eq. 17 gives zero current, for the first n states: the single
    root on [-5, 5] V by pymrm.checks.find_roots (raises unless exactly one)."""
    if (aa, ac) not in _ROOTS:
        _ROOTS[(aa, ac)] = np.array([
            one_root(lambda e, c=c, s=s, u=u: float(i_eq17(u + e, c, s, aa, ac)), -5.0, 5.0,
                     "zero-current potential", n_scan=32, xtol=1e-15)
            for c, s, u in zip(c_r[:400], cs_r[:400], U2_r[:400])])
    return _ROOTS[(aa, ac)][:n]


def eq16_dev(pairs):
    """Largest |zero-current potential of Eq. 17 - Eq. 16| over 200 states and the pairs, V."""
    return float(max(np.max(np.abs(zero_current_roots(aa, ac, 200))) for aa, ac in pairs))


ident = pd.DataFrame([(aa, ac, aa + ac, d_eq17_vs_eq30(aa, ac), eq16_dev([(aa, ac)]))
                      for aa, ac in SUM1_PAIRS + NOT1_PAIRS],
                     columns=["alpha_a", "alpha_c", "sum", "max rel. dev. Eq.17 vs Eq.30",
                              "max |E_eq - Eq.16| / V"])
print(ident.to_string(index=False,
                      formatters={"max rel. dev. Eq.17 vs Eq.30": "{:.1e}".format,
                                  "max |E_eq - Eq.16| / V": "{:.1e}".format}))

IDENT_SYM = max(d_eq17_vs_eq30(aa, ac) for aa, ac in SUM1_PAIRS)
IDENT_EQ16 = eq16_dev(SUM1_PAIRS)
claims.bound("eq17_vs_eq30_max_rel", IDENT_SYM, max=1e-10, basis="identity",
             note="round-off of exponentials of order 10 over 4000 states")
claims.bound("eq17_zero_current_vs_eq16_V", IDENT_EQ16, max=1e-13, basis="tolerance",
             unit="V", note="find_roots xtol 1e-15 V")
print(f"\nwith alpha_a + alpha_c = 1 : Eq. 17 == Eq. 30 to below 1e-10 relative "
      f"(here {IDENT_SYM:.1e}),")
print(f"                             and Eq. 17's zero-current potential IS Eq. 16 "
      f"to below 1e-13 V (here {IDENT_EQ16:.1e})")
print(f"with alpha_a + alpha_c != 1 : the zero-current potential is off by up to "
      f"{1e3*eq16_dev(NOT1_PAIRS):.0f} mV.")

# %% [markdown]
# The failure for $\alpha_a+\alpha_c\ne1$ is not numerical noise, and it has a
# closed form: the zero-current potential of Eq. 17 sits at
# $U' + \frac{RT}{F(\alpha_a+\alpha_c)}\ln\frac{c_T-c_s}{c_s}$, so its offset from
# Eq. 16 is $\frac{RT}{F}\ln\frac{c_T-c_s}{c_s}\left(\frac{1}{\alpha_a+\alpha_c}-1\right)$.
# Checking that closed form against the numerically located root tests the
# statement rather than just observing it.

# %%
OFF_PAIRS = [(0.3, 0.7), (0.5, 0.4), (0.5, 0.6), (0.25, 0.25), (1.0, 1.0)]


def offset_check(aa, ac, drop_factor=False):
    """max |root - closed-form offset| over 400 states, V; drop_factor omits (1/n - 1)."""
    cs_s = cs_r[:400]
    closed = (np.zeros_like(cs_s) if drop_factor else
              RTF * np.log((C_T - cs_s) / cs_s) * (1.0 / (aa + ac) - 1.0))
    return float(np.max(np.abs(zero_current_roots(aa, ac) - closed)))


off = pd.DataFrame([(aa, ac, offset_check(aa, ac), offset_check(aa, ac, drop_factor=True))
                    for aa, ac in OFF_PAIRS],
                   columns=["alpha_a", "alpha_c", "max |numeric - closed form| / V",
                            "max |offset| / V"])
print(off.to_string(index=False, formatters={
    "max |numeric - closed form| / V": "{:.1e}".format,
    "max |offset| / V": "{:.4f}".format}))
OFFSET_ERR = max(offset_check(aa, ac) for aa, ac in OFF_PAIRS)
claims.bound("alpha_offset_closed_form_V", OFFSET_ERR, max=1e-13, basis="tolerance",
             unit="V", note="find_roots xtol 1e-15 V")
print(f"closed form vs located root: below 1e-13 V (here {OFFSET_ERR:.1e})")


# %% [markdown]
# ### Check 2, the things the law must do for *any* parameters
#
# Four statements the law must satisfy for every parameter set. They differ a great
# deal in how much they can detect, and the section after this one measures that
# rather than claiming it:
#
# * $i=0$ at $\eta_s=0$, exactly, for every $(\alpha_a,\alpha_c,i_0)$, this one is
#   sensitive, because it is evaluated at the closed-form equilibrium potential and
#   a wrong prefactor moves the root;
# * at equilibrium the two exponential terms are individually equal, the anodic
#   and cathodic partial currents balance, which is what makes $i_0$ an *exchange*
#   current density rather than a fitting constant. **This is an algebraic
#   consequence of $\alpha_a+\alpha_c=1$ and nothing else**: it is exact for any
#   concentration, it is the same statement as the first item at $n=1$, and the
#   only defect it can detect is $\alpha_a+\alpha_c\ne1$, which Check 1 already
#   tests. It is reported as the identity it is, not as evidence;
# * their common value is Eq. 30's $i_0$, evaluated **at asymmetric $\alpha$ as
#   well as at Doyle's symmetric pair**, because at $\alpha_a=\alpha_c$ the
#   comparison cannot see an exponent swap at all;
# * $R_{\rm ct}=RT/(nFi_0)$ with $n=\alpha_a+\alpha_c$, against a numerical
#   derivative of the printed Eq. 17;
# * the arcsinh inversion is exact for $\alpha_a=\alpha_c$, fed forward through
#   Eq. 17.

# %%
# eta_equilibrium (closed-form zero-current potential, U' + (RT/nF) ln((c_T-c_s)/c_s))
# and i0_general (the common value of the partial currents at equilibrium, any alpha)
# are Kinetics methods in model.py.
# (a) exactly zero current at the equilibrium potential, for every alpha pair
ZERO_PAIRS = [(0.5, 0.5), (0.3, 0.7), (0.15, 0.85), (0.9, 0.1), (0.25, 0.25)]


def zero_current(drop_n=False):
    """max |i(eta_eq)| / i0 over 5 alpha pairs and all states; drop_n omits 1/n in eta_eq."""
    return float(max(np.max(np.abs(i_eq17(eta_equilibrium(cs_r, aa, ac, drop_n), c_r, cs_r,
                                          aa, ac) / i0_general(c_r, cs_r, aa, ac)))
                     for aa, ac in ZERO_PAIRS))


zero_max = zero_current()
claims.bound("zero_current_at_zero_eta_rel", zero_max, max=1e-12, basis="identity")
print(f"(a) |i(eta_s=0)| / i0, worst over 5 alpha pairs x {n_pts} states : below 1e-12 "
      f"(here {zero_max:.1e})")

# (b) the two partial currents balance at equilibrium, and their common value is
#     Eq. 30's i_0. Evaluated at ASYMMETRIC alpha as well: at alpha_a = alpha_c the
#     comparison is blind to any swap of the two exponents (see the defect table).
ALPHA_PAIRS_N1 = SUM1_PAIRS


def d_partials(aa, ac, e17pre=None, e17brk=None, e30pre=None, e30tail=None):
    """Check 2(b): (|i_a - i_c| / i_a, |(i_a + i_c)/2 / Eq. 30 - 1|), max over states."""
    z = f_RT * (U2_r - U_prime(cs_r))
    pa, pc = (aa, ac) if e17pre is None else e17pre
    ba, bc = (aa, ac) if e17brk is None else e17brk
    p = F_CONST * P["k_2"] * (C_MAX - c_r) ** pc * c_r ** pa
    anod, cath = p * cs_r * np.exp(ba * z), p * (C_T - cs_r) * np.exp(-bc * z)
    ref = i0_eq30(c_r, cs_r, aa, ac, pre=e30pre, tail=e30tail)
    return (float(np.max(np.abs(anod - cath) / np.abs(anod))),
            float(np.max(np.abs(0.5 * (anod + cath) / ref - 1.0))))


def partials(**defect):
    """d_partials, worst over the four pairs with n = 1; a swap is given as True."""
    out = [d_partials(aa, ac, **{k: (ac, aa) for k in defect}) for aa, ac in ALPHA_PAIRS_N1]
    return max(o[0] for o in out), max(o[1] for o in out)


bal, i0_rel = partials()
claims.bound("partial_currents_balance_rel_structural", bal, max=1e-12, basis="identity",
             structural=True, note="exact for alpha_a + alpha_c = 1 at any state; blind to "
             "every exponent swap (row defect_worst_partial_currents_balance)")
claims.bound("i0_from_partials_vs_eq30_rel", i0_rel, max=1e-12, basis="identity")
print(f"(b) |i_a - i_c| / i_a at equilibrium, 4 alpha pairs with n = 1  : below 1e-12 "
      f"(here {bal:.1e})")
print("    (structural: exact whenever alpha_a + alpha_c = 1, for any concentration)")
print(f"    the common value IS Eq. 30's i_0, same 4 pairs, to below 1e-12 (here {i0_rel:.1e})")

# (c) R_ct against a central difference on Eq. 17, at the true equilibrium
h = 1e-7
RCT_PAIRS = [(0.5, 0.5), (0.3, 0.7), (0.15, 0.85), (0.5, 0.6), (1.0, 1.0)]


def rct_check(pairs=RCT_PAIRS, n_forced=None):
    """max |central dI/d(eta) RT / (n F i0) - 1|; n_forced replaces n = aa + ac (a defect)."""
    worst = 0.0
    for aa, ac in pairs:
        e_eq = eta_equilibrium(cs_r, aa, ac)
        di = (i_eq17(e_eq + h, c_r, cs_r, aa, ac) - i_eq17(e_eq - h, c_r, cs_r, aa, ac)) / (2 * h)
        n = aa + ac if n_forced is None else n_forced
        # d i / d eta at equilibrium = (alpha_a + alpha_c) F i_0 / RT, for any alpha pair
        worst = max(worst, float(np.max(np.abs(di * RTF / (n * i0_general(c_r, cs_r, aa, ac))
                                               - 1.0))))
    return worst


RCT_ID = rct_check()
claims.bound("Rct_identity_max_rel", RCT_ID, max=1e-7, basis="tolerance",
             note="central difference with h = 1e-7 V: round-off eps*|i|/(h di/deta)")
print(f"(c) numerical dI/d(eta) vs (n F i0)/RT, 5 alpha pairs: below 1e-7 (here {RCT_ID:.1e})")


# and for n = 1 the general i_0 is Eq. 30's
def i0_gen_check(aa=0.3, ac=0.7):
    """max |i0_general / Eq. 30 - 1| at one alpha pair."""
    return float(np.max(np.abs(i0_general(c_r, cs_r, aa, ac) / i0_eq30(c_r, cs_r, aa, ac) - 1.0)))


I0_GEN_VS_EQ30 = i0_gen_check()
claims.bound("i0_general_vs_eq30_rel", I0_GEN_VS_EQ30, max=1e-12, basis="identity")
print(f"    i_0 (general) vs Eq. 30 at alpha = (0.3, 0.7), sum = 1: below 1e-12 "
      f"(here {I0_GEN_VS_EQ30:.1e})")

# (d) the arcsinh inversion, checked forward through Eq. 17 (no root-finder
#     tolerance enters: the overpotential the closed form gives must reproduce
#     the current it was asked for)
i_targets = np.geomspace(1e-2, 1e3, 120)


def asinh_check(aa=None):
    """max |Eq. 17(arcsinh inversion) / i - 1|; aa != 1/2 (with ac = 1 - aa) is a defect:
    the inversion assumes alpha_a = alpha_c."""
    worst = 0.0
    for c_v, cs_v in [(1000.0, 290.0), (500.0, 14500.0), (3000.0, 28000.0)]:
        i0v = float(i0_eq30(c_v, cs_v, aa, None if aa is None else 1 - aa))
        U2v = float(U_ocp(cs_v / C_T))
        for sgn in (+1.0, -1.0):
            eta_closed = RTF / alpha * np.arcsinh(sgn * i_targets / (2 * i0v))
            back = i_eq17(U2v + eta_closed, c_v, cs_v, aa, None if aa is None else 1 - aa)
            worst = max(worst, float(np.max(np.abs(back / (sgn * i_targets) - 1.0))))
    return worst


ASINH_ID = asinh_check()
claims.bound("arcsinh_inversion_max_rel", ASINH_ID, max=1e-8, basis="identity",
             note="exact for alpha_a = alpha_c; the floor is the cancellation printed below")
print(f"(d) arcsinh inversion, fed back through Eq. 17, 5 decades of i: below 1e-8 "
      f"(here {ASINH_ID:.1e})")

# the floor on (d) is not the identity, it is cancellation: Eq. 17 as printed is a
# difference of two nearly equal exponentials, so it loses digits as eta_s -> 0.
canc = []
for scale in (1e0, 1e-2, 1e-4, 1e-6):
    e = scale * RTF
    lhs = np.exp(alpha * f_RT * e) - np.exp(-alpha * f_RT * e)
    canc.append((scale, abs(lhs / (2 * np.sinh(alpha * f_RT * e)) - 1.0)))
print("    difference-of-exponentials vs 2 sinh(.), as eta_s/(RT/F) shrinks:")
for s, v in canc:
    print(f"      eta_s = {s:8.0e} RT/F  ->  relative error {v:.2e}")
print("    -> evaluate the sinh form, not the printed difference, below ~1e-4 RT/F.")


# %% [markdown]
# ### What Checks 1 and 2 can and cannot catch, measured, by injecting defects
#
# An agreement number is worth exactly what it would cost to break, so the
# transcription is re-coded wrong, one mistake at a time, and the same two residuals
# are re-measured. Every number in the table below is computed in the cell; none is
# quoted. Each defect is a plausible mis-reading of a scanned equation, which
# $\alpha$ sits on which concentration, and each is applied to **one** of the two
# equations unless the row says otherwise.
#
# Two columns per check: at Doyle's own $\alpha_a=\alpha_c=\tfrac12$, and the worst
# over the sweep's asymmetric pairs with $\alpha_a+\alpha_c=1$. The result decides
# how the checks may be described:
#
# * **At $\alpha_a=\alpha_c$ every one of these defects is invisible.** Swapping two
#   exponents that are equal changes nothing, and both residuals stay at their
#   baseline for all six rows. A reader who assumes the check bites at the
#   published parameters would be wrong; all of its discriminating power comes from
#   the hypothetical asymmetric rows.
# * **Swapping the prefactor in *both* equations at once is invisible even there**,
#   because $(c_{\max}-c)^{\alpha_c}c^{\alpha_a}$ is written as the same code
#   expression in `i_eq17` and `i0_eq30` and cancels in the ratio. That half of
#   Eq. 30 is *not* tested by Check 1; only the $(c_T-c_s)^{\alpha_a}c_s^{\alpha_c}$
#   half is. Swapping it in either equation *alone* does break, which is what makes
#   the $\alpha_A$ reading of Eq. 17 decidable.
# * **The partial-current balance never moves at all**, under any of them. It is an
#   exact algebraic consequence of $\alpha_a+\alpha_c=1$; it is reported as an
#   identity and carries no evidential weight.

# %%
# d_eq17_vs_eq30 and d_partials above take each exponent assignment as an argument
DEFECTS = [
    ("none - as printed", {}),
    ("Eq. 30 tail exponents swapped", dict(e30tail="swap")),
    ("Eq. 30 prefactor exponents swapped", dict(e30pre="swap")),
    ("Eq. 17 prefactor read as c^alpha_c", dict(e17pre="swap")),
    ("the same swap in BOTH prefactors", dict(e17pre="swap", e30pre="swap")),
    ("Eq. 17 bracket exponents swapped", dict(e17brk="swap")),
]
ASYM_PAIRS = [(0.3, 0.7), (0.7, 0.3), (0.15, 0.85)]

def defect_row(spec):
    """The five columns of the defect table for one set of swapped exponent assignments."""
    def kw(aa, ac):
        return {k: ((ac, aa) if v == "swap" else v) for k, v in spec.items()}
    b_s, i_s = d_partials(0.5, 0.5, **kw(0.5, 0.5))
    ba_, ia_ = zip(*[d_partials(aa, ac, **kw(aa, ac)) for aa, ac in ASYM_PAIRS])
    return {"Eq.17 vs Eq.30, a=(.5,.5)": d_eq17_vs_eq30(0.5, 0.5, **kw(0.5, 0.5)),
            "Eq.17 vs Eq.30, asym": max(d_eq17_vs_eq30(aa, ac, **kw(aa, ac))
                                        for aa, ac in ASYM_PAIRS),
            "i0 vs Eq.30, a=(.5,.5)": i_s, "i0 vs Eq.30, asym": max(ia_),
            "partials balance, worst": max((b_s,) + ba_)}


def defect_metrics(defects=DEFECTS):
    """The defect-table summaries; defects=[] for every row removes the injected defects."""
    rows = {label: defect_row(spec) for label, spec in defects}
    sym = [r["Eq.17 vs Eq.30, a=(.5,.5)"] for r in rows.values()]
    get = lambda i, col: rows[DEFECTS[i][0]][col]          # noqa: E731
    return {"e30_tail": get(1, "Eq.17 vs Eq.30, asym"), "e17_pre": get(3, "Eq.17 vs Eq.30, asym"),
            "both_pre": get(4, "Eq.17 vs Eq.30, asym"), "sym_worst": max(sym) - min(sym),
            "bal_worst": max(r["partials balance, worst"] for r in rows.values()),
            "i0_asym": get(1, "i0 vs Eq.30, asym"), "rows": rows}


DEFECTS_REMOVED = [(label, {}) for label, _ in DEFECTS]
DM = defect_metrics()
dtab = pd.DataFrame([{"injected defect": k, **v} for k, v in DM["rows"].items()])
print(dtab.to_string(index=False, formatters={c: "{:.1e}".format for c in dtab.columns[1:]}))

DEF_BASE = DM["rows"][DEFECTS[0][0]]["Eq.17 vs Eq.30, asym"]
DEF_E30_TAIL, DEF_E17_PRE, DEF_BOTH_PRE = DM["e30_tail"], DM["e17_pre"], DM["both_pre"]
DEF_SYM_WORST, DEF_BAL_WORST, DEF_I0_ASYM = DM["sym_worst"], DM["bal_worst"], DM["i0_asym"]
claims.value("defect_eq30_tail_swap_eq17_vs_eq30", DEF_E30_TAIL, rtol=1e-6)
claims.value("defect_eq17_prefactor_swap_eq17_vs_eq30", DEF_E17_PRE, rtol=1e-6)
claims.value("defect_eq30_tail_swap_i0_from_partials", DEF_I0_ASYM, rtol=1e-6)
claims.bound("defect_both_prefactors_swap_eq17_vs_eq30", DEF_BOTH_PRE, max=1e-10,
             basis="identity", structural=True,
             note="a blind spot: the same swap on both sides cancels in the ratio")
claims.bound("defect_worst_at_symmetric_alpha_eq17_vs_eq30", DEF_SYM_WORST, max=1e-12,
             basis="identity", structural=True,
             note="a blind spot: swapping equal exponents changes nothing")
claims.bound("defect_worst_partial_currents_balance", DEF_BAL_WORST, max=1e-12,
             basis="identity", structural=True,
             note="a blind spot: no exponent swap moves the balance")
print()
print(f"at alpha = (0.5, 0.5): the six rows of Check 1 span below 1e-12 "
      f"(here {DEF_SYM_WORST:.1e})")
print("   -> every exponent swap is a no-op at Doyle's own parameters.")
print(f"at asymmetric alpha : baseline below 1e-10 (here {DEF_BASE:.1e}), and each "
      f"single-equation defect lifts it to O(1)")
print(f"   -> that is where the power is. Reading Eq. 17's prefactor as c^alpha_c "
      f"gives {DEF_E17_PRE:.2f},")
print(f"      which is why the printed alpha_A can be resolved as alpha_a at all.")
print(f"swapping BOTH prefactors together              : below 1e-10 (here {DEF_BOTH_PRE:.1e}; "
      f"cancels, undetectable)")
print(f"partial-current balance, worst over all rows   : below 1e-12 (here {DEF_BAL_WORST:.1e}; "
      f"never moves; structural)")

# %% [markdown]
# ### Check 3, a printed number that tests Eq. 30, with the conductivity divided out
#
# Doyle prints two dimensionless groups evaluated for his own cell:
#
# $$
# \delta = \frac{\alpha_aFI\delta_c}{RT}\Big(\tfrac1\kappa+\tfrac1\sigma\Big),
# \qquad
# \nu^2 = (\alpha_a+\alpha_c)\frac{Fai_o\delta_c^2}{RT}
# \Big(\tfrac1\kappa+\tfrac1\sigma\Big),
# $$
#
# with the printed values loaded from the stated-results file.
#
# Their ratio removes $\big(\tfrac1\kappa+\tfrac1\sigma\big)$ entirely:
#
# $$
# \frac{\nu^2}{\delta} = \frac{\alpha_a+\alpha_c}{\alpha_a}\,
# \frac{a\,i_{o2}\,\delta_c}{I}.
# $$
#
# This matters here because the conductivity is the one input the paper never
# prints: `J3.4` had to reconstruct it. The ratio needs none of it, so it is a
# **clean printed-number check on Eq. 30's exchange current density**, and it is
# independent of the reconstruction that `J3.4`'s own agreement rests on.
#
# **It is not, however, an independent second test.** `J3.4` inverts Eq. 28 for
# $\kappa$ and inverts Eq. 29 for $\kappa$ and reports the two answers differing
# (its `kappa_route_spread_pct`, printed below). That is the *same* disagreement,
# viewed from the other side: dividing the
# two equations is precisely what cancels $\kappa$, so
# $1/(1+\text{deviation})-1$ here and $\kappa_{28}/\kappa_{29}-1$ there must agree
# to whatever $\sigma$ contributes, which is nothing. The cell below computes both
# and shows they do. Read as $\kappa$, the disagreement is an ambiguity in a
# reconstructed conductivity; read as $\nu^2/\delta$, it is a statement about
# $i_{o2}$ and the specific area with no $\kappa$ in it. **One disagreement, two
# readings**, a reader should not count it twice.
#
# **The state the groups were evaluated at is not ambiguous.** The sentence
# immediately below Eq. 30 reads, verbatim off the 600 dpi render of page 1530:
#
# > For these calculations, the concentrations are taken to be at their initial
# > values. In our case, we find that $\delta$ = 1.95 and $\nu$ = 68 …
#
# and Figure 2's caption pins the initial solid state, which is the $u_0$ of
# Doyle's own Table II and of the `J3.4` dataset used here. So the
# comparison below is made at exactly the state the paper specifies, and "which
# state?" is *not* available as an explanation of the residual. Nor is $\kappa$: it
# cancels identically from $\nu^2/\delta$, as the algebra two paragraphs up shows.
# Both of the obvious candidates are ruled out by the paper's own text, and the
# residual is left standing and unexplained.

# %%
I_REF = 10.0                                  # A/m2, the current Doyle quotes these at


def nu2_study(u0=P["u_0"]):
    """Check 3 at initial solid state u0 (Doyle's u_0 as printed; another value is a
    defect): i_o2 there, nu^2/delta against Eq. 30, J3.4's two kappa routes, and the
    J3.5 scale ratio, all from the same i_o2."""
    i0_0 = float(i0_eq30(P["c_0"], u0 * C_T))
    lhs_ = S["nu"] ** 2 / S["delta"]
    rhs_ = (P["alpha_a"] + P["alpha_c"]) / P["alpha_a"] * a_spec * i0_0 * dc / I_REF
    dev = rhs_ / lhs_ - 1.0
    kd_ = 1.0 / (S["delta"] / (alpha * F_CONST * I_REF * dc / (R_GAS * T)) - 1.0 / sigma)
    kn_ = 1.0 / (S["nu"] ** 2 / (2 * alpha * F_CONST * a_spec * i0_0 * dc ** 2 / (R_GAS * T))
                 - 1.0 / sigma)
    pre0 = F_CONST * P["k_2"] * (C_MAX - P["c_0"]) ** alpha * P["c_0"] ** alpha
    return {"i0": i0_0, "lhs": lhs_, "rhs": rhs_, "dev": dev, "kd": kd_, "kn": kn_,
            "kappa_spread": kd_ / kn_ - 1.0, "as_kappa": 1.0 / (1.0 + dev) - 1.0,
            "scale_ratio": i0_0 / (pre0 * C_T)}


NU2 = nu2_study()
i0_init, lhs, rhs = NU2["i0"], NU2["lhs"], NU2["rhs"]
NU2_OVER_DELTA = NU2["dev"]
claims.value("nu2_over_delta_dev_frac", NU2_OVER_DELTA, rtol=1e-6)

print(f"evaluated at the state the paper specifies: c = c_0 = {P['c_0']:.0f} mol/m3, "
      f"u_0 = {P['u_0']:.2f}\n('the concentrations are taken to be at their initial "
      "values', p. 1530)")
print(f"i_o2 at the initial state (Eq. 30)      = {i0_init:.2f} A/m2")
print(f"paper:      nu^2 / delta = {S['nu']:.0f}^2 / {S['delta']:.2f} = {lhs:.1f}")
print(f"Eq. 30:  (n/alpha_a) a i_o2 delta_c / I = {rhs:.1f}")
print(f"deviation (model - paper)/paper         = {100*NU2_OVER_DELTA:+.1f} %")
print()

# could rounding of the printed groups account for it? nu^2/delta is quadratic in
# nu, so a half-unit either way is worth +-1.5 %, and delta to two decimals +-0.3 %.
lo = (S["nu"] - 0.5) ** 2 / (S["delta"] + 0.005)
hi = (S["nu"] + 0.5) ** 2 / (S["delta"] - 0.005)
nu_needed = np.sqrt(rhs * S["delta"])
print(f"rounding band on nu^2/delta (nu +-0.5, delta +-0.005) : "
      f"{lo:.0f} to {hi:.0f}  ({100*(lo/lhs-1):+.1f} % to {100*(hi/lhs-1):+.1f} %)")
print(f"nu that would close the gap                           : {nu_needed:.1f}, "
      f"{100*(nu_needed/S['nu']-1):+.1f} % from the printed {S['nu']:.0f}")
print("-> rounding does NOT account for it.")
print()

# what WOULD close it? Both candidates contradict a printed value, so neither is
# adopted; they are computed to show the size of the change that would be needed.
i0_needed = i0_init / (1.0 + NU2_OVER_DELTA)
U0_NEEDED = one_root(lambda u: float(i0_eq30(P["c_0"], u * C_T)) - i0_needed,
                     1e-6, 0.5, "u_0 that closes the gap", log=True, xtol=1e-14)
RS_NEEDED = P["R_s"] * (1.0 + NU2_OVER_DELTA)      # a = 3(1-eps)/R_s, so a i_o2 ~ 1/R_s
print(f"the ratio would close if a i_o2 rose by {100*(i0_needed/i0_init-1):+.1f} %, i.e. if")
print(f"  the initial solid state were u_0 = {U0_NEEDED:.4f}  (Table II and Fig. 2's "
      f"caption print {P['u_0']:.2f}), or")
print(f"  the particle radius were R_s = {1e6*RS_NEEDED:.3f} um       (Table II prints "
      f"{1e6*P['R_s']:.1f} um)")
print("Neither is available: both contradict a printed value, so neither is adopted.")
print("The residual is REPORTED, not explained.")
print()

# it is arithmetically J3.4's kappa disagreement, not a second independent test.
# Both inversions are J3.4's, reproduced here (Eq. 28 and Eq. 29 solved for kappa).
KAPPA_SPREAD, NU2_AS_KAPPA = NU2["kappa_spread"], NU2["as_kappa"]
claims.value("j34_kappa_route_spread", KAPPA_SPREAD, rtol=1e-6)
claims.value("nu2_over_delta_as_kappa_ratio", NU2_AS_KAPPA, rtol=1e-6)
claims.bound("j34_kappa_spread_vs_committed_rel",
             abs(100 * KAPPA_SPREAD / J34["kappa_route_spread_pct"] - 1.0), max=1e-6,
             basis="chosen", note="the same arithmetic as J3.4's committed claim, on another run")
print(f"J3.4 reads the same disagreement as two routes to the unprinted kappa:")
print(f"  kappa_eff from Eq. 28 / kappa_eff from Eq. 29 - 1 = {100*KAPPA_SPREAD:+.2f} %"
      f"  (J3.4's committed kappa_route_spread_pct: {J34['kappa_route_spread_pct']:.2f} %)")
print(f"  1/(1 + this page's deviation) - 1                 = {100*NU2_AS_KAPPA:+.2f} %")
print(f"  the two differ by {abs(KAPPA_SPREAD/NU2_AS_KAPPA - 1):.1e} relative, which is "
      f"the 1/sigma term.")
print("-> ONE disagreement seen from two sides, not two independent tests. This page's "
      "reading\n   is the more informative of the two only because kappa is absent "
      "from it.")
print()
print(f"R_ct on the interfacial area at that state = {R_ct(i0_init)*1e4:.4f} Ohm cm2")
print(f"and per unit of separator area             = "
      f"{R_ct(i0_init)/(a_spec*dc)*1e4:.3e} Ohm cm2")
print(f"1/kappa + 1/sigma is absent from the ratio; sigma alone would contribute "
      f"{1.0/sigma:.1e} against 1/kappa ~ 20, i.e. nothing.")


# %% [markdown]
# ## Results
#
# ### Where the Tafel approximation fails, in closed form
#
# Define the deviation the way a user meets it: the approximation is used to
# *predict a current* from a known overpotential, so the error is
#
# $$
# \varepsilon_T = \frac{i_{\rm Tafel}-i}{i}
# = \frac{e^{\alpha_a f\eta_s}}{e^{\alpha_a f\eta_s}-e^{-\alpha_c f\eta_s}}-1
# = \frac{1}{e^{(\alpha_a+\alpha_c)f\eta_s}-1},
# \qquad f=\frac{F}{RT},
# $$
#
# which inverts exactly:
#
# $$
# \boxed{\;\eta_s^{\rm Tafel}(\varepsilon)
# = \frac{RT}{(\alpha_a+\alpha_c)F}\,\ln\!\Big(1+\frac1\varepsilon\Big)\;}
# $$
#
# Two things follow that are worth stating, because they are not what is usually
# assumed.
#
# **The threshold does not depend on $\alpha_a$ and $\alpha_c$ separately, only on
# their sum.** A strongly asymmetric electrode with $\alpha_a=0.15$ reaches the
# Tafel regime at exactly the same overpotential as a symmetric one, provided
# $\alpha_a+\alpha_c$ is the same. The asymmetry changes the Tafel *slope*, not the
# Tafel *onset*. This is exact, not asymptotic.
#
# **The deviation convention matters.** If instead the approximation is used to
# *infer an overpotential* from a measured current, which is what reading a Tafel
# plot does, the error is in $\eta_s$, the $\alpha$-dependence returns, and the two
# conventions do not give the same number. Every threshold on this page uses the
# current convention.
#
# `Kinetics.eta_tafel_threshold` in `model.py` is the closed form;
# `eta_tafel_threshold_numeric` is the single root of $\varepsilon_T-\varepsilon$
# on the full law, located by `pymrm.checks.find_roots` on $[10^{-6}, 5]$ V and
# raising unless the scan finds exactly one.

# %%
levels = [0.10, 0.05, 0.01]
tab = []
for eps in levels:
    row = {"deviation": f"{100*eps:.0f} %"}
    row["closed form, 298.15 K / mV"] = 1e3 * eta_tafel_threshold(eps, 1.0, RTF_298)
    row["closed form, 373.15 K / mV"] = 1e3 * eta_tafel_threshold(eps, 1.0, RTF)
    row["x RT/F"] = eta_tafel_threshold(eps, 1.0, RTF) / RTF
    tab.append(row)
tafel_tab = pd.DataFrame(tab)
print("Tafel approximation, error in the predicted CURRENT at given eta_s "
      "(alpha_a + alpha_c = 1):")
print(tafel_tab.to_string(index=False, float_format=lambda v: f"{v:8.2f}"))

THR = {}                                     # thresholds by rtf, so every T row re-runs all


def thresholds(rtf=RTF):
    """The four headline thresholds in mV at one RT/F."""
    if rtf not in THR:
        THR[rtf] = {"tafel10": 1e3 * eta_tafel_threshold(0.10, 1.0, rtf),
                    "tafel1": 1e3 * eta_tafel_threshold(0.01, 1.0, rtf),
                    "lin5": 1e3 * eta_linear_threshold(0.05, 0.5, 0.5, +1, rtf),
                    "lin1": 1e3 * eta_linear_threshold(0.01, 0.5, 0.5, +1, rtf)}
    return THR[rtf]


claims.value("eta_tafel_10pct_mV", thresholds()["tafel10"], rtol=1e-6, unit="mV")
claims.value("eta_tafel_1pct_mV", thresholds()["tafel1"], rtol=1e-6, unit="mV")


# the closed form against root-finding on the full law, over a range of alpha
def tafel_check(form="log1p"):
    """max |numeric / closed form - 1|, 6 alpha_a x 3 levels; form="log" writes ln(1/eps)."""
    return float(max(abs(eta_tafel_threshold_numeric(e, aa, 1.0 - aa)
                         / eta_tafel_threshold(e, 1.0, form=form) - 1.0)
                     for aa in [0.15, 0.25, 0.35, 0.5, 0.65, 0.85] for e in levels))


TAFEL_CLOSED = tafel_check()
claims.bound("tafel_threshold_closed_vs_numeric", TAFEL_CLOSED, max=1e-12, basis="tolerance",
             note="find_roots xtol 1e-16 V on a threshold of order 0.1 V")
print(f"\nclosed form vs root-finding on the full law, 6 alpha_a x 3 levels: "
      f"below 1e-12 (here {TAFEL_CLOSED:.1e})")
print("-> the threshold is genuinely independent of how alpha is split.")


# %% [markdown]
# ### Where the linear approximation fails
#
# The linear limit has no closed-form threshold, because
# $\varepsilon_L = \dfrac{nf\eta_s}{e^{\alpha_af\eta_s}-e^{-\alpha_cf\eta_s}}-1$
# does not invert. For the symmetric case it reduces to
# $\varepsilon_L = u/\sinh u - 1$ with $u=\alpha f\eta_s$, which is negative, the
# linear form always **under**predicts the current, so it always **over**predicts
# the overpotential needed to drive a given current. That direction is the reason
# it is the conservative approximation of the two.
#
# Unlike the Tafel threshold, this one *does* depend on how $\alpha$ is split, and
# it is not symmetric in the sign of $\eta_s$ once $\alpha_a\ne\alpha_c$.
#
# One trap in computing it, which the code below is explicit about. For
# $\alpha_a\ne\alpha_c$ the deviation is **not monotonic** in $\eta_s$: the linear
# form first overshoots the exact law, then crosses back through it, then falls
# away. So $|\varepsilon_L|-\varepsilon$ has up to four roots, and a bracketing
# solver handed the whole interval returns one of the later ones. At
# $\alpha_a=0.3$ taking the last crossing reports the 5 % boundary more than an
# order of magnitude too far out (both are printed below), in the *unsafe*
# direction, which is exactly the sort of error that does not announce itself.
# The threshold is the **first** crossing: `find_roots` scans
# $[10^{-6}, 4]$ V on 1200 logarithmic points and the smallest root is kept.
#
# Expanding $\varepsilon_L$ gives two regimes, both of which are checked against
# the root-find below. With $b=\tfrac12(\alpha_a-\alpha_c)$ and
# $c=\tfrac16(\alpha_a^3+\alpha_c^3)$,
# $\varepsilon_L \simeq -b\,\hat\eta+(b^2-c)\hat\eta^2$ with
# $\hat\eta=F\eta_s/RT$, so
#
# * **symmetric** ($b=0$, $c=\tfrac1{24}$): $\hat\eta = \sqrt{24\varepsilon}$;
# * **strongly asymmetric** ($|b|\hat\eta$ dominant):
#   $\hat\eta = 2\varepsilon/|\alpha_a-\alpha_c|$.
#
# Between them the two error terms have opposite signs and partially cancel, and
# that cancellation produces a genuine, and genuinely fragile, anomaly, see
# below.

# %%
claims.value("eta_linear_5pct_mV", thresholds()["lin5"], rtol=1e-6, unit="mV")
claims.value("eta_linear_1pct_mV", thresholds()["lin1"], rtol=1e-6, unit="mV")
rows = []
for eps in levels:
    rows.append({
        "deviation": f"{100*eps:.0f} %",
        "298.15 K / mV": 1e3 * eta_linear_threshold(eps, 0.5, 0.5, +1, RTF_298),
        "373.15 K / mV": 1e3 * eta_linear_threshold(eps, 0.5, 0.5, +1, RTF),
        "x RT/F": eta_linear_threshold(eps, 0.5, 0.5, +1, RTF) / RTF,
        "sqrt(24 eps)": np.sqrt(24.0 * eps),
    })
lin_tab = pd.DataFrame(rows)
print("Linear approximation, symmetric alpha = 0.5, error in the predicted CURRENT:")
print(lin_tab.to_string(index=False, float_format=lambda v: f"{v:8.3f}"))
def sym_asymptote(k=24.0):
    """max |sqrt(k eps) / root-found threshold (in RT/F) - 1| over the three levels."""
    return float(max(abs(np.sqrt(k * e) / (eta_linear_threshold(e, 0.5, 0.5) / RTF) - 1.0)
                     for e in levels))


SYM_ASYMPT = sym_asymptote()
claims.value("linear_sym_asymptote_rel", SYM_ASYMPT, rtol=1e-6)
print(f"-> the small-overpotential expansion sqrt(24 eps) is good to "
      f"{100*SYM_ASYMPT:.1f} % at these levels.")

print("\nAsymmetry: the 5 % threshold, anodic and cathodic side, at 373.15 K")
asym = []
for aa in [0.5, 0.45, 0.42, 0.4, 0.3, 0.2, 0.1]:
    ac = 1.0 - aa
    asym.append({"alpha_a": aa, "alpha_c": ac,
                 "anodic +eta / mV": 1e3 * eta_linear_threshold(0.05, aa, ac, +1),
                 "cathodic -eta / mV": 1e3 * eta_linear_threshold(0.05, aa, ac, -1),
                 "2 eps/|da| x RT/F / mV": (np.nan if aa == ac else
                                           1e3 * RTF * 2 * 0.05 / abs(aa - ac)),
                 "Tafel 5 % / mV": 1e3 * eta_tafel_threshold(0.05, aa + ac)})
asym_tab = pd.DataFrame(asym)
print(asym_tab.to_string(index=False, float_format=lambda v: f"{v:10.2f}"))

LIN_5 = {aa: 1e3 * eta_linear_threshold(0.05, aa, 1 - aa, +1) for aa in (0.5, 0.42, 0.3)}
claims.value("linear_5pct_alpha0p42_mV", LIN_5[0.42], rtol=1e-6, unit="mV")
claims.value("linear_5pct_alpha0p3_mV", LIN_5[0.3], rtol=1e-6, unit="mV")
LIN_NAIVE_MV = 1e3 * eta_linear_threshold(0.05, 0.3, 0.7, crossing="last")
claims.value("break_linear_threshold_last_crossing_mV", LIN_NAIVE_MV, rtol=1e-6, unit="mV",
             note="the 5 % threshold at alpha_a = 0.3 taken as the last crossing: the trap")
print(f"\n5 % threshold at alpha_a = 0.3, first crossing: {LIN_5[0.3]:.2f} mV; "
      f"the last crossing: {LIN_NAIVE_MV:.2f} mV ({LIN_NAIVE_MV/LIN_5[0.3]:.0f}x)")


# free identity: swapping alpha_a and alpha_c must mirror the two sides exactly
def mirror(swap=True):
    """max |thr_anodic(aa) / thr_cathodic(1 - aa) - 1|; swap=False omits the alpha swap."""
    return float(max(abs(eta_linear_threshold(e, aa, 1 - aa, +1)
                         / eta_linear_threshold(e, *((1 - aa, aa) if swap else (aa, 1 - aa)), -1)
                         - 1.0)
                     for aa in (0.1, 0.2, 0.3, 0.42, 0.45) for e in levels))


MIRROR = mirror()
claims.bound("linear_sign_convention_rel", MIRROR, max=1e-12, basis="tolerance",
             note="the two sides are the same root-find with the sign handling reversed")
print(f"\nsign-convention check, thr_anodic(alpha_a) = thr_cathodic(1-alpha_a): below 1e-12 "
      f"(here {MIRROR:.1e})")
print("   (exact by construction IF the sign handling is right, and visibly "
      "non-zero if it is not:\n    a check on the implementation, not an "
      "independent identity)")


# the asymptote, away from the cancellation window
def asym_asymptote(aas=(0.1, 0.2, 0.3)):
    """max |2 eps RT/(|aa - ac| F) / threshold - 1| over aas and the three levels."""
    return float(max(abs(RTF * 2 * e / abs(aa - (1 - aa))
                         / eta_linear_threshold(e, aa, 1 - aa, +1) - 1.0)
                     for aa in aas for e in levels))


def narrowing(crossing="first"):
    """Anodic 5 % threshold at alpha_a = 0.5 over the one at alpha_a = 0.1."""
    return float(eta_linear_threshold(0.05, 0.5, 0.5, +1, crossing=crossing)
                 / eta_linear_threshold(0.05, 0.1, 0.9, +1, crossing=crossing))


ASYM_ASYMPT = asym_asymptote()
claims.value("linear_asym_asymptote_rel", ASYM_ASYMPT, rtol=1e-6)
print(f"strongly asymmetric asymptote  eta = 2 eps RT / (|alpha_a - alpha_c| F) : "
      f"good to {100*ASYM_ASYMPT:.1f} % for alpha_a <= 0.3")
LIN_ASYM_NARROW = narrowing()
claims.value("linear_window_narrowing_0p5_to_0p1", LIN_ASYM_NARROW, rtol=1e-6)
print(f"asymmetry narrows the linear window by {LIN_ASYM_NARROW:.1f}x between "
      f"alpha_a = 0.5 and 0.1, while the Tafel column does not move at all.")

# where the first-lobe overshoot reaches the tolerance: a root in alpha_a of the
# overshoot peak (model.first_lobe_peak, a bounded maximisation refined from a scan)
print("\nwhere the first-lobe overshoot reaches the tolerance - i.e. where the "
      "threshold jumps:")


def jump(eps, slope=None):
    """alpha_a in [0.30, 0.4999999] at which the first-lobe overshoot peak equals eps;
    slope replaces n = 1 in the linear form (a defect passes 2 alpha_a)."""
    return one_root(lambda a: first_lobe_peak(a, None if slope is None else slope(a)) - eps,
                    0.30, 0.4999999, "jump in alpha_a", xtol=1e-12)


JUMP = {eps: jump(eps) for eps in (0.01, 0.05)}
claims.value("linear_jump_alpha_a_1pct", JUMP[0.01], rtol=1e-6)
claims.value("linear_jump_alpha_a_5pct", JUMP[0.05], rtol=1e-6)
for eps in (0.01, 0.05):
    print(f"  {100*eps:4.0f} % tolerance : alpha_a = {JUMP[eps]:.4f}")
OVERSHOOT_045 = first_lobe_peak(0.45)
claims.value("linear_overshoot_peak_alpha0p45", OVERSHOOT_045, rtol=1e-6)
print(f"at alpha_a = 0.45 the overshoot peaks at only {100*OVERSHOOT_045:.2f} %, so a "
      f"5 % tolerance is never\nviolated on the first lobe there and the crossing "
      f"jumps to the far branch.")

# %% [markdown]
# ### The cancellation window, and why not to rely on it

# %%
GAP = {r: (1e3 * eta_linear_threshold(0.05, 0.5, 0.5, +1, r), 1e3 * eta_tafel_threshold(0.05, 1.0, r))
       for r in (RTF, RTF_298)}
display(Markdown(f"""
The table above is not monotonic in $\\alpha_a$, and that is not a numerical
artefact. On the anodic side with $\\alpha_a$ a little below $\\tfrac12$, the
linear form first runs *above* the exact law (the $-b\\hat\\eta$ term, $b<0$) and
then *below* it (the $-c\\hat\\eta^2$ term). The two cancel: at $\\alpha_a=0.45$ the
early overshoot peaks at only {100*OVERSHOOT_045:.2f} %, so a 5 % tolerance is never
violated there and the first crossing jumps out to the far branch. The 5 % boundary
is therefore **{LIN_5[0.42]:.0f} mV at $\\alpha_a=0.42$ against {LIN_5[0.5]:.0f} mV at
$\\alpha_a=0.50$**: the linear approximation is at its best not at the symmetric
point but slightly off it.

It is a true statement about the function and a bad thing to design on, because
the position of the jump is set entirely by the tolerance: the cell above locates
it as the root in $\\alpha_a$ of the overshoot peak, and it sits at
$\\alpha_a={JUMP[0.05]:.4f}$ for a 5 % tolerance but at $\\alpha_a={JUMP[0.01]:.4f}$
for a 1 % one. Either side of the window the threshold collapses to the
$2\\varepsilon/|\\alpha_a-\\alpha_c|$ asymptote. Quote the asymptote, not the anomaly.

One identity comes free and is checked: swapping $\\alpha_a$ and $\\alpha_c$ must
mirror the two branches exactly, $\\eta^{{+}}(\\alpha_a)=\\eta^{{-}}(1-\\alpha_a)$. It
holds to the root-finder's tolerance, which is a check on the sign handling as
much as on the algebra.

### The gap between the two limits

Between the linear threshold and the Tafel threshold lies a band where *neither*
approximation is within the stated tolerance and the full law is required. At the
5 % level and $\\alpha_a+\\alpha_c=1$ that band runs from
{GAP[RTF][0]/(1e3*RTF):.1f} $RT/F$ to {GAP[RTF][1]/(1e3*RTF):.1f} $RT/F$: from
{GAP[RTF][0]:.0f} to {GAP[RTF][1]:.0f} mV at 100 °C, {GAP[RTF_298][0]:.0f} to
{GAP[RTF_298][1]:.0f} mV at 25 °C. Nothing about that band depends on the electrode,
the reaction, or the exchange current density: it is a property of the functional
form and scales only with $RT/F$ and $1/(\\alpha_a+\\alpha_c)$.
"""))

# %%
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))

eh = np.linspace(-6.0, 6.0, 2400)                 # eta_s in units of RT/F
for aa, style in [(0.5, "-"), (0.25, "--"), (0.45, ":")]:
    ac = 1.0 - aa
    full = np.exp(aa * eh) - np.exp(-ac * eh)
    lin = (aa + ac) * eh
    ax[0].plot(eh, np.abs(lin / full - 1.0), "C0" + style, lw=1.4,
               label=rf"linear, $\alpha_a$={aa}")
# every alpha split with alpha_a + alpha_c = 1 gives the SAME Tafel deviation
taf_dev = np.abs(1.0 / np.expm1(np.abs(eh)))
ax[0].plot(eh, taf_dev, "C3-", lw=2.2, label="Tafel (all splits coincide)")
for eps in levels:
    ax[0].axhline(eps, color="0.6", lw=0.8)
    ax[0].text(-5.9, eps * 1.2, f"{100*eps:.0f} %", fontsize=8, color="0.4")
ax[0].set_yscale("log")
ax[0].set_ylim(1e-4, 3.0)
ax[0].set_xlabel(r"$\eta_s\,F/RT$")
ax[0].set_ylabel("relative deviation in current")
ax[0].set_title("Both limits, and the band where neither works")
ax[0].legend(fontsize=7, loc="lower left")

aa_grid = np.linspace(0.05, 0.95, 91)
for eps, col in zip(levels, ["C2", "C1", "C4"]):
    lin_hi = np.array([eta_linear_threshold(eps, a, 1 - a, +1) for a in aa_grid]) / RTF
    lin_lo = np.array([eta_linear_threshold(eps, a, 1 - a, -1) for a in aa_grid]) / RTF
    ax[1].plot(aa_grid, lin_hi, col + "-", lw=1.5, label=f"linear +, {100*eps:.0f} %")
    ax[1].plot(aa_grid, lin_lo, col + "--", lw=1.0, label=f"linear -, {100*eps:.0f} %")
    ax[1].axhline(eta_tafel_threshold(eps, 1.0) / RTF, color=col, ls=":", lw=1.8)
    ax[1].text(0.06, eta_tafel_threshold(eps, 1.0) / RTF * 1.06,
               f"Tafel {100*eps:.0f} %", fontsize=7, color=col)
ax[1].set_yscale("log")
ax[1].set_ylim(2e-2, 1e1)
ax[1].set_xlabel(r"$\alpha_a$   (with $\alpha_c = 1-\alpha_a$)")
ax[1].set_ylabel(r"$|\eta_s|\,F/RT$  at the threshold")
ax[1].set_title("Validity boundaries vs. the transfer coefficient")
ax[1].legend(fontsize=6.5, ncol=3, loc="lower center")
fig.tight_layout()
plt.show()

print("Left: the Tafel curve is drawn once because every split with "
      "alpha_a + alpha_c = 1 gives\nexactly the same deviation. The linear curves "
      "do depend on the split, and the dotted\none (alpha_a = 0.45) shows the "
      "cancellation dip that produces the anomaly.")
print("Right: Tafel boundaries (dotted) are horizontal lines - independent of the "
      "split.\nThe linear boundaries (solid = anodic, dashed = cathodic) fall away "
      "on both sides,\nwith the cancellation spike just off alpha_a = 0.5.")


# %% [markdown]
# ## The published cell as an operating envelope
#
# The thresholds above are properties of the equation. Whether they bite depends on
# where a real electrode sits, and that is a reactor question. Below is `J3.4`'s
# cell, ported with the kinetic law made selectable, run to extract the one thing
# that page never reports: the distribution of $\eta_s$.

# %%
BRUG = eps_c ** 1.5
KAPPA = kappa_from_delta(P, S["delta"]) / BRUG   # bulk value; kappa_eff = KAPPA * eps^1.5
print(f"kappa (bulk, reconstructed as on J3.4 from Eq. 28) = {KAPPA:.4f} S/m")
COLS = ["u", "V", "t", "salt", "max|eta_s|", "mean eta_s", "eta_s1", "mean|eta_s|"]
_RUNS = {}


def run(I=10.0, par=None, **kw):
    """Cell(...).march(I, ...) at Doyle's parameters (par overrides some), cached by its
    arguments: returns (out, dts, v_probe at u = 0.5). Constructor arguments of Cell
    (n_s, n_c, nu_factor, kinetics, face_mean, leak) and march arguments go in kw."""
    key = (I, tuple(sorted((par or {}).items())), tuple(sorted(kw.items())))
    if key not in _RUNS:
        k_ = kin if not par else Kinetics({**P, **par})
        ckw = {k: kw[k] for k in ("n_s", "n_c", "nu_factor", "kinetics", "face_mean", "leak",
                                  "newton_tol") if k in kw}
        mkw = {k: (np.asarray(v) if k == "dt_seq" else v) for k, v in kw.items() if k not in ckw}
        cell = Cell(k_, KAPPA, **ckw)
        out, dts = cell.march(I, u_probe=0.5, **mkw)
        _RUNS[key] = (out, dts, cell.v_probe)
    return _RUNS[key]


plt.rcParams["figure.max_open_warning"] = 0
print("Cell and run() ready.")

# %% [markdown]
# ### What the cell actually runs at
#
# The exact law (`eq17`) is solved at the three currents of Doyle's Figure 2 that
# reach the cutoff, and the surface overpotential is read out of the solution.
#
# **Where the readout stops is part of the result.** `Cell.march` is `J3.4`'s and
# halts at `v_stop = 1.55` V, below the cutoff $V_{\text{cutoff}}$ the paper
# states and the parameter file carries; that margin exists so the stepper has
# somewhere to go, not because the cell is meant to be operated there.
# Every envelope number below is therefore read at $V_{\text{cutoff}}$, which is
# what `J3.4` does with the same runs, and the march endpoint is printed beside it
# and labelled as the stepper artefact it is.

# %%
# Read at the paper's own cutoff, which is a loaded parameter, not a typed one.
V_CUT = P["V_cutoff"]


def envelope(I, par=None, level=V_CUT, **kw):
    """The run at current I read inside its envelope: up to whichever comes first in
    utilisation, the start of the transport-limited tail (u = 0.80) or the crossing of
    `level` (the root of the interpolant of V(u), unique or it raises)."""
    r = run(I, par, **kw)[0]
    u_cut = u_at_level(r[:, 0], r[:, 1], level)
    m = r[:, 0] <= min(0.80, u_cut)
    return {"u_cut": u_cut, "mean": float(r[m, 7].max()), "peak": float(r[m, 4].max()),
            "anode": float(r[m, 6].max()), "peak_full": float(r[:, 4].max()),
            "u_end": float(r[-1, 0]), "v_end": float(r[-1, 1]),
            "u19": u_at_level(r[:, 0], r[:, 1], 1.9) if r[:, 1].min() < 1.9 else np.nan}


runs = {I: run(I)[0] for I in (5.0, 10.0, 20.0)}
dtseq = {I: run(I)[1] for I in runs}
for I in runs:
    print(f"I = {I:4.0f} A/m2 : {len(runs[I]):4d} steps, "
          f"u_final = {runs[I][-1, 0]:.3f}, V_final = {runs[I][-1, 1]:.3f} V")
ENV = {I: envelope(I) for I in runs}
U_CUT = {I: e["u_cut"] for I, e in ENV.items()}
env_tab = pd.DataFrame([{
    "I / A m-2": I, "i_loc/i_0 if uniform": (I / (a_spec * dc)) / i0_init,
    f"u at {V_CUT} V cutoff": e["u_cut"], "cathode mean |eta_s| / mV": 1e3 * e["mean"],
    "cathode peak |eta_s| / mV": 1e3 * e["peak"], "anode eta_s1 / mV": 1e3 * e["anode"],
    "cathode peak PAST the cutoff / mV": 1e3 * e["peak_full"],
    "march endpoint u (v_stop=1.55)": e["u_end"]} for I, e in ENV.items()])
print()
print(env_tab.to_string(index=False, float_format=lambda v: f"{v:12.4g}"))

ETA_CATH_10, ETA_CATH_MEAN_10, ETA_AN_10 = ENV[10.0]["peak"], ENV[10.0]["mean"], ENV[10.0]["anode"]
ETA_AN_20, ETA_CATH_20 = ENV[20.0]["anode"], ENV[20.0]["peak"]
ETA_CATH_20_FULL, U_END_20, V_END_20 = ENV[20.0]["peak_full"], ENV[20.0]["u_end"], ENV[20.0]["v_end"]
U19_20 = ENV[20.0]["u19"]
claims.value("cell_cathode_mean_eta_s_mV", 1e3 * ETA_CATH_MEAN_10, rtol=1e-3, unit="mV")
claims.value("cell_cathode_peak_eta_s_mV", 1e3 * ETA_CATH_10, rtol=1e-3, unit="mV")
claims.value("cell_anode_eta_s_I10_mV", 1e3 * ETA_AN_10, rtol=1e-3, unit="mV")
claims.value("cell_anode_eta_s_I20_mV", 1e3 * ETA_AN_20, rtol=1e-3, unit="mV")
claims.value("cell_cathode_peak_I20_cutoff_mV", 1e3 * ETA_CATH_20, rtol=1e-3, unit="mV")
claims.value("u_cutoff_I20", U_CUT[20.0], rtol=1e-3)
claims.value("u_cutoff_I10", U_CUT[10.0], rtol=1e-3)
claims.value("u_1p9V_I20", U19_20, rtol=1e-3)
claims.bound("u_cutoff_I10_vs_J34_abs", abs(U_CUT[10.0] - J34["u_at_cutoff_I10"]), max=2e-3,
             basis="chosen", note="the same cell and cutoff on J3.4; its run and readout differ")
print()
print("The cathode PEAK sits in the first cell next to the separator, where the "
      "reaction crowds,\nand in the first TIME STEP. It is a boundary layer in both, "
      "it is converged in\nneither, and Validation quantifies each - together with "
      "the sustained value that\nthe cathode runs at once the initial transient has "
      "decayed, which is stable.\nThe cell-averaged value is not grid-sensitive, and "
      "all of them are reported.")

# Reconcile against the stated results loaded from J3.4, rather than restating
# the march endpoint: the CSV carries the paper's own number for this collapse.
print()
print(f"I = 20 A/m2, where the collapse is read:")
print(f"  V = 1.9 V crossing                 : u = {U19_20:.4f}")
print(f"  V = {V_CUT} V, the paper's cutoff       : u = {U_CUT[20.0]:.4f}")
print(f"  march endpoint, V = {V_END_20:.3f} V       : u = {U_END_20:.4f}   "
      f"<- a v_stop artefact, not a published condition")
print(f"  Doyle p. 1529, loaded from the CSV : u_at_sharp_drop = "
      f"{S['u_at_sharp_drop']:.2f}  (\"about 30%\")")
print(f"  I = 10 A/m2 at the same cutoff     : u = {U_CUT[10.0]:.4f}   "
      f"(J3.4's committed u_at_cutoff_I10: {J34['u_at_cutoff_I10']:.4f})")

# %% [markdown]
# ### The result: the two electrodes of the same cell are three decades apart
#
# This is the useful part, and it is not what a reader would guess from the fact
# that both electrodes obey the same equation.
#
# **State the tolerance first, because the verdict depends on it.** At $I=10$ both
# electrodes are inside the *5 %* linear window, the cathode by decades and the
# foil by a small factor. They are in genuinely *different* regimes only at the
# 1 % tolerance, where the foil is already out, or at $I=20$, where the foil is
# out at the 5 % level too. What does not depend on the tolerance is the
# **separation** in $|\eta_s|$ between the two electrodes of one cell, at one
# instant, under one equation. The cell below prints each electrode against both
# thresholds so the reader can apply their own.
#
# **The porous cathode is deep inside the 5 % linear window.** Its interfacial
# area $a\delta_c$ per unit separator area divides the local current density,
# $I/(a\delta_c)$, which is small against Eq. 30's exchange current density (both
# computed below under "Superficial or interfacial?"). The reaction is not
# uniform (Doyle's own $\nu$ says the ohmic drop dominates, so it crowds towards
# the separator), and the solve measures the real peak rather than assuming
# uniformity, but even the peak leaves $|\eta_s| \ll RT/F$. Doyle says as much in
# words, $k_2$ was "set equal to a value corresponding to a nearly reversible
# situation", and this puts a number on it.
#
# **The lithium foil is at the edge of it.** It is flat, so its current density is
# the cell's, $I$, against $i_{o1}$. At $I=20$ it is past the 5 % linear threshold
# and heading for the 10 % one.
#
# The cell after next substitutes each limit law in turn: linearising the cathode
# kinetics leaves the cell potential unchanged to within the solver tolerance,
# linearising the anode kinetics costs millivolts, and **Tafel is wrong at both**,
# catastrophically at the cathode, where $i\ll i_0$ is exactly the regime Tafel is
# not for.
#
# At $I=20$ the cell is transport-limited and collapses, and the cathode peak
# overpotential rises as the electrolyte depletes and $i_0$ falls with it. That
# collapse is read at Doyle's own cutoff; the cell below prints where it lands and
# what the peak is there.

# %%
print("Cathode, at the peak |eta_s| the solve finds (I = 10 A/m2):")
print(f"  |eta_s|                     = {1e3*ETA_CATH_10:.4f} mV "
      f"= {ETA_CATH_10/RTF:.3e} x RT/F")
print(f"  linear-limit error there    = "
      f"{abs((2*alpha*f_RT*ETA_CATH_10)/np.sinh(alpha*f_RT*ETA_CATH_10)/2 - 1):.3e}")
print(f"  Tafel-limit error there     = "
      f"{1.0/np.expm1(f_RT*ETA_CATH_10):.4g}  (i.e. {100/np.expm1(f_RT*ETA_CATH_10):.0f} %)")
print()
print(f"Anode (Doyle Eqs. 6-8, i_o1 = {P['i_0_1']} A/m2):")
an = []
for I in (5.0, 10.0, 20.0):
    i01 = P["i_0_1"]
    ex = 2 * RTF * np.arcsinh(I / (2 * i01))
    li = RTF * I / i01
    ta = 2 * RTF * np.log(I / i01)
    an.append({"I / A m-2": I, "I / i_o1": I / i01,
               "exact eta_s1 / mV": 1e3 * ex,
               "linear / mV": 1e3 * li, "linear error / mV": 1e3 * (li - ex),
               "Tafel / mV": 1e3 * ta, "Tafel error / mV": 1e3 * (ta - ex)})
an_tab = pd.DataFrame(an)
print(an_tab.to_string(index=False, float_format=lambda v: f"{v:10.3f}"))
LIN_ANODE_20 = float(an_tab.loc[an_tab["I / A m-2"] == 20.0,
                                "linear error / mV"].iloc[0])
TAF_ANODE_10 = float(an_tab.loc[an_tab["I / A m-2"] == 10.0,
                                "Tafel error / mV"].iloc[0])

# the verdict depends on the tolerance, so both are printed against both thresholds
print()
print("Where each electrode sits, against the linear thresholds (the verdict is "
      "tolerance-dependent):")
thr = {e: eta_linear_threshold(e, 0.5, 0.5) for e in (0.05, 0.01)}
FOIL_LIN_ERR = {}
for name, e_s in [("porous cathode, mean", ETA_CATH_MEAN_10),
                  ("porous cathode, peak", ETA_CATH_10),
                  ("lithium foil, I = 10", ETA_AN_10),
                  ("lithium foil, I = 20", ETA_AN_20)]:
    u = alpha * f_RT * e_s
    err = abs(u / np.sinh(u) - 1.0)
    FOIL_LIN_ERR[name] = err
    print(f"  {name:22s}: |eta_s| = {1e3*e_s:7.3f} mV = "
          f"{e_s/thr[0.05]:6.3f} x the 5 % threshold, {e_s/thr[0.01]:6.3f} x the 1 % one"
          f"   (linear error {100*err:.3g} %)")
def separation(an10=None, cath_mean10=None):
    """Decades between the foil and the cell-averaged cathode |eta_s| at I = 10."""
    return float(np.log10((ETA_AN_10 if an10 is None else an10)
                          / (ETA_CATH_MEAN_10 if cath_mean10 is None else cath_mean10)))


ETA_SEP_DECADES = separation()
ETA_SEP_DECADES_PEAK = float(np.log10(ETA_AN_10 / ETA_CATH_10))
claims.value("eta_separation_decades", ETA_SEP_DECADES, rtol=1e-3)
FOIL_RATIO = {(I, e): (ETA_AN_10 if I == 10 else ETA_AN_20) / thr[e]
              for I in (10, 20) for e in (0.05, 0.01)}
for (I, e), v in FOIL_RATIO.items():
    claims.value(f"foil_I{I}_over_linear_{int(100*e)}pct", v, rtol=1e-3)
print(f"  separation between the two electrodes at I = 10: "
      f"{ETA_SEP_DECADES:.2f} decades on the cell-averaged cathode value, "
      f"{ETA_SEP_DECADES_PEAK:.2f} on its peak")

display(Markdown(f"""
**The $I=20$ collapse, read at the cutoff.** The cell reaches
$V={V_CUT}$ V, Doyle's stated cutoff, loaded from the `J3.4` parameter file, at
$u={U_CUT[20.0]:.4f}$, and crosses 1.9 V at $u={U19_20:.4f}$. Doyle p. 1529 says
the potential "drops sharply when about 30 % of the cathode material is utilized",
and `u_at_sharp_drop = {S['u_at_sharp_drop']:.2f}` sits in the stated-results file
loaded at the top of this page. Inside that envelope the cathode peak
overpotential reaches **{1e3*ETA_CATH_20:.3f} mV**, still
{thr[0.01]/ETA_CATH_20:.1f} times inside the 1 % linear window and
{thr[0.05]/ETA_CATH_20:.1f} times inside the 5 % one. The march itself runs on to
$u={U_END_20:.4f}$ at $V={V_END_20:.3f}$ V, where `v_stop` stops the stepper,
{1e3*(V_CUT-V_END_20):.0f} mV past the published cutoff and outside the operating
envelope; the cathode peak over that whole run is {1e3*ETA_CATH_20_FULL:.1f} mV.
The largest kinetic overpotential in these runs is the lithium foil's at $I=20$,
{1e3*ETA_AN_20:.1f} mV.
"""))

# %% [markdown]
# ### Superficial or interfacial? The factor that decides the regime
#
# The classification above uses the **interfacial** current density,
# $I/(a\delta_c)$, and that is the physically correct one. It is worth putting a
# number on what the distinction is worth here, because the dramatic version of the
# warning, that ignoring $a\delta_c$ would put this cathode in the Tafel regime,
# is false, and a reader who believed it would mis-diagnose their own cell in
# exactly the direction the warning is meant to prevent.

# %%
print("Superficial vs interfacial current density, at Doyle's I = "
      f"{I_REF:.0f} A/m2 and i_0 = {i0_init:.2f} A/m2:")
i_sup, i_int = I_REF, I_REF / (a_spec * dc)
for name, i_loc in [("superficial I", i_sup), ("interfacial I/(a delta_c)", i_int)]:
    e_s = RTF / alpha * np.arcsinh(i_loc / (2 * i0_init))
    u = alpha * f_RT * e_s
    print(f"  {name:26s}: i = {i_loc:8.4f} A/m2, i/i_0 = {i_loc/i0_init:9.3e}, "
          f"eta_s = {1e3*e_s:8.4f} mV,\n"
          f"  {'':26s}  linear error {100*abs(u/np.sinh(u)-1):9.3g} %, "
          f"Tafel error {100/np.expm1(f_RT*e_s):9.3g} %")
print(f"  the factor between them is a delta_c = {a_spec*dc:.0f}")
print()
e_sup = RTF / alpha * np.arcsinh(i_sup / (2 * i0_init))
print(f"So even the SUPERFICIAL reading is inside the linear window: "
      f"{e_sup/thr[0.01]:.2f} x the 1 % threshold")
print(f"({1e3*thr[0.01]:.2f} mV) and {e_sup/thr[0.05]:.2f} x the 5 % one "
      f"({1e3*thr[0.05]:.2f} mV). It is not near Tafel:")
i_taf5 = 2 * i0_init * np.sinh(alpha * f_RT * eta_tafel_threshold(0.05))
print(f"the 5 % Tafel boundary sits at i/i_0 = {i_taf5/i0_init:.2f}, which is "
      f"{i_taf5/i_sup:.1f} x this cell's")
print("superficial current and " f"{i_taf5/i_int:.0f} x its interfacial one.")
print()
print("What a delta_c changes here is depth inside the linear window, not the "
      "regime:")
print(f"  linear error, superficial reading  = "
      f"{100*abs(alpha*f_RT*e_sup/np.sinh(alpha*f_RT*e_sup)-1):.3g} %")
e_int = RTF / alpha * np.arcsinh(i_int / (2 * i0_init))
print(f"  linear error, interfacial reading  = "
      f"{100*abs(alpha*f_RT*e_int/np.sinh(alpha*f_RT*e_int)-1):.3g} %")
print("Both are inside every tolerance on this page. The factor would decide the "
      "verdict only")
print(f"at about {i_taf5/i_sup:.0f} x this current, where the superficial reading "
      "reaches the Tafel boundary")
print("while the interface is still linear - which is exactly when getting it "
      "wrong costs something.")
A_DELTA_C = float(a_spec * dc)


def superficial(par=None):
    """Linear-limit errors of the superficial and interfacial readings at I_REF, and the
    Tafel headroom of the superficial one, with i_o2 of Eq. 30 at the initial state."""
    k_ = kin if not par else Kinetics({**P, **par})
    i0 = float(k_.i0_eq30(P["c_0"], P["u_0"] * C_T))
    err = {}
    for name, i_loc in (("sup", I_REF), ("int", I_REF / (a_spec * dc))):
        x = np.arcsinh(i_loc / (2 * i0))      # alpha f eta_s of the symmetric law
        err[name] = float(abs(x / np.sinh(x) - 1.0))
    i_t5 = 2 * i0 * np.sinh(alpha * f_RT * eta_tafel_threshold(0.05))
    return {"sup": err["sup"], "int": err["int"], "headroom": float(i_t5 / I_REF),
            "i0": i0, "e_sup": RTF / alpha * np.arcsinh(I_REF / (2 * i0))}


SUPER = superficial()
SUP_LIN_ERR, INT_LIN_ERR, SUP_TAFEL_HEADROOM = SUPER["sup"], SUPER["int"], SUPER["headroom"]
claims.value("superficial_reading_linear_err", SUP_LIN_ERR, rtol=1e-6)
claims.value("interfacial_reading_linear_err", INT_LIN_ERR, rtol=1e-6)
claims.value("superficial_tafel_headroom_factor", SUP_TAFEL_HEADROOM, rtol=1e-6)
claims.value("a_delta_c", A_DELTA_C, rtol=1e-9)
claims.value("superficial_eta_s_mV", 1e3 * SUPER["e_sup"], rtol=1e-6, unit="mV")
claims.value("superficial_i_over_i0", I_REF / i0_init, rtol=1e-6)
claims.value("i0_initial_A_m2", i0_init, rtol=1e-6, unit="A/m2")

# the same ratio appears on J3.5 as C_r, on a different i_0 scale. Reconcile it here
# so the two pages are not read as disagreeing.
PRE0 = F_CONST * P["k_2"] * (C_MAX - P["c_0"]) ** alpha * P["c_0"] ** alpha
C_R_J35 = I_REF / (a_spec * PRE0 * C_T * dc)      # J3.5's C_r_doyle, verbatim
C_R_HERE = i_int / i0_init
I0_SCALE_RATIO = NU2["scale_ratio"]
claims.value("j35_C_r_scale_ratio", I0_SCALE_RATIO, rtol=1e-6)
claims.value("interfacial_i_over_i0", C_R_HERE, rtol=1e-6)
print()
print("The same ratio on J3.5, and why the two pages print different numbers:")
print(f"  this page,  i/(a delta_c i_o2 from Eq. 30 at the initial state) = "
      f"{C_R_HERE:.3e}")
print(f"  J3.5 C_r,   i/(a delta_c * PRE0 * c_T)                          = "
      f"{C_R_J35:.3e}")
print(f"  ratio = {C_R_J35/C_R_HERE:.4f} = i_o2 / (PRE0 c_T) = "
      f"{I0_SCALE_RATIO:.4f} = sqrt(u_0 (1 - u_0)) = "
      f"{np.sqrt(P['u_0']*(1-P['u_0'])):.4f}")
print("J3.5 follows Marquis in scaling on a TYPICAL exchange current, with the solid")
print("concentration factors replaced by c_T; this page evaluates Eq. 30 at the state")
print(f"Doyle specifies, where the solid is only {100*P['u_0']:.0f} % lithiated. The two "
      "differ by exactly")
print("sqrt(u_0(1-u_0)) and are the same quantity. Neither is wrong; they are not "
      "comparable.")

# %% [markdown]
# ### Substituting the limit into the reactor, and paying for it
#
# The pointwise error is not the number an engineer needs; the number they need is
# what it does to the cell. So the same cell is re-solved with the limit laws in
# place of Eq. 17, and the cell potential is compared. The cathode substitution
# changes the field equations and needs a full re-solve; the anode substitution is
# explicit, so all three anode variants come from one solve.
#
# Two things have to be right before this comparison means anything.
#
# * **The two runs must step in time identically.** The reference run's step
#   sequence is recorded and replayed, so the comparison is step-for-step with no
#   interpolation, and a difference in the step schedule cannot enter it. (On
#   this cell the stepper takes the same steps for both laws over the compared window: a run without the replay
#   gives the same deviation, printed below.)
# * **Only one electrode may be substituted at a time.** The cell potential
#   contains the anode overpotential, so a `Cell` that switches both laws at once
#   reports the anode's error as the cathode's; the break table at the end of
#   Validation has that row. Here the anode always uses the
#   exact law and the switch acts on the porous cathode alone.
#
# The Tafel-substituted cathode is worth predicting before running it. Tafel with
# $i\ll i_0$ has no solution near $\eta_s=0$: the single exponential must be driven
# down to $i/i_0$ of the interfacial reading, which needs
# $\eta_s = \frac{RT}{\alpha_cF}\ln(i_0/|i|)$ of the *wrong sign* (printed below).
# The approximation does not merely lose accuracy: it reverses the sign of the
# polarisation.

# %%
pred_tafel = RTF / alpha * np.log(i0_init / (I_REF / (a_spec * dc)))
claims.value("tafel_cathode_predicted_eta_mV", 1e3 * pred_tafel, rtol=1e-6, unit="mV")
print(f"predicted Tafel cathode overpotential at I = 10 : {1e3*pred_tafel:+.1f} mV")

ref = runs[I_REF]
seq = tuple(dtseq[I_REF])    # replay the reference step sequence exactly


def sub_run(law, replay=True, anode_law="bv"):
    """The I_REF cell with `law` at the cathode: on the reference step sequence
    (replay=False steps adaptively, a defect), anode in `anode_law` (exact by default)."""
    kw = {"kinetics": law, "anode_law": anode_law}
    if replay:
        kw["dt_seq"] = seq
    return run(I_REF, **kw)[0]


sub = {law: sub_run(law) for law in ("bv", "linear", "tafel")}
for law in sub:
    print(f"  {law:7s}: {len(sub[law]):4d} steps (reference: {len(ref)}), "
          f"cathode mean eta_s at step 0 = {1e3*sub[law][0, 5]:+9.3f} mV, "
          f"V = {sub[law][0, 1]:.4f} V")


def tafel_pred(law="tafel"):
    """|first-step mean eta_s of the `law` run / the uniform-reaction Tafel prediction - 1|."""
    return float(abs(sub_run(law)[0, 5] / pred_tafel - 1.0))


TAFEL_PRED_DEV = tafel_pred()
claims.value("tafel_first_step_vs_prediction_frac", TAFEL_PRED_DEV, rtol=0.02)
print(f"\n  -> the Tafel run's first-step overpotential matches the "
      f"uniform-reaction prediction to {100*TAFEL_PRED_DEV:.2f} %")


def dv(law, lo=0.10, hi=0.70, **kw):
    """Cell-potential difference to the reference over lo < u < hi, step for step
    (the kw of sub_run inject the two comparison defects)."""
    s_ = sub_run(law, **kw)
    n = min(len(ref), len(s_))
    m = (ref[:n, 0] > lo) & (ref[:n, 0] < hi)
    return float(np.max(np.abs(s_[:n, 1][m] - ref[:n, 1][m])))


BV_VS_EQ17, LIN_CATH_DV, TAF_CATH_DV = dv("bv"), dv("linear"), dv("tafel")
claims.bound("bv_vs_eq17_cell_potential_mV", 1e3 * BV_VS_EQ17, max=1e-6, basis="tolerance",
             unit="mV", note="the same law written two ways; newton tol 1e-9 on scaled residuals")
claims.bound("linear_cathode_dV_mV", 1e3 * LIN_CATH_DV, max=1e-6, basis="tolerance", unit="mV",
             note="the linearisation error at |eta_s| below 1 mV is under the solver tolerance")
claims.value("tafel_cathode_dV_mV", 1e3 * TAF_CATH_DV, rtol=1e-3, unit="mV")
LIN_NOREPLAY = dv("linear", replay=False)
claims.bound("linear_cathode_dV_noreplay_mV", 1e3 * LIN_NOREPLAY, max=1e-6, basis="tolerance",
             unit="mV", note="the same comparison with both runs stepping adaptively")
print(f"  linear cathode, both runs stepping adaptively       : below 1e-6 mV "
      f"(here {1e3*LIN_NOREPLAY:.1e})")
print()
print(f"cell potential, max |deviation| over u = 0.10-0.70 at I = {I_REF:.0f} A/m2,")
print("compared step-for-step on the same time grid:")
print(f"  standard form (Eq. 30 + Eq. 6) vs Eq. 17 as printed : "
      f"below 1e-6 mV (here {1e3*BV_VS_EQ17:.1e})")
print(f"  linear cathode kinetics                             : "
      f"below 1e-6 mV (here {1e3*LIN_CATH_DV:.1e})")
print(f"  Tafel  cathode kinetics                             : "
      f"{1e3*TAF_CATH_DV:9.1f} mV")


# the anode substitution is explicit in I, so it needs no re-solve
def anode_shift(I, law, exact_law="bv", i01=None):
    """Cell-potential shift when the foil uses `law` instead of the exact law, V
    (exact_law=law makes the substitution exact by construction; i01 overrides i_o1)."""
    i01 = P["i_0_1"] if i01 is None else i01
    eta = {"bv": 2 * RTF * np.arcsinh(I / (2 * i01)), "linear": RTF * I / i01,
           "tafel": 2 * RTF * np.log(I / i01)}
    return float(-(eta[law] - eta[exact_law]))      # V = ... - eta_s1, so dV = -d(eta_s1)


print()
print("anode substitution, cell potential shift at I = 10 and 20 A/m2:")
AN_SHIFT = {(I, law): anode_shift(I, law) for I in (10.0, 20.0) for law in ("linear", "tafel")}
for (I, law), v in AN_SHIFT.items():
    print(f"  I = {I:4.0f}, {law:7s}: dV = {1e3*v:+9.2f} mV")
claims.value("linear_anode_dV_I20_mV", 1e3 * abs(AN_SHIFT[(20.0, "linear")]), rtol=1e-6, unit="mV")
claims.value("linear_anode_dV_I10_mV", 1e3 * abs(AN_SHIFT[(10.0, "linear")]), rtol=1e-6, unit="mV")

# %%
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))

for law, lab, st in [("eq17", "Eq. 17 as printed", "k-"),
                     ("linear", "linear kinetics", "C0--"),
                     ("tafel", "Tafel kinetics", "C3-.")]:
    r = ref if law == "eq17" else sub[law]
    ax[0].plot(r[:, 0], r[:, 1], st, lw=1.6, label=lab)
ax[0].set_xlabel("utilisation $u$")
ax[0].set_ylabel("cell potential / V")
ax[0].set_title(f"J3.4 cell at $I$ = {I_REF:.0f} A/m$^2$, kinetics substituted")
ax[0].legend(fontsize=8)

for I, col in zip((5.0, 10.0, 20.0), ["C0", "C1", "C3"]):
    r = runs[I]
    ax[1].plot(r[:, 0], 1e3 * r[:, 4], col + "-", lw=1.5,
               label=f"cathode peak, $I$ = {I:.0f}")
    ax[1].plot(r[:, 0], 1e3 * r[:, 6], col + "--", lw=1.2,
               label=f"anode, $I$ = {I:.0f}")
for eps, lab in [(0.05, "5 %"), (0.10, "10 %")]:
    thr = 1e3 * eta_linear_threshold(eps, 0.5, 0.5, +1)
    ax[1].axhline(thr, color="0.5", lw=0.9, ls=":")
    ax[1].text(0.02, thr * 1.06, f"linear {lab}", fontsize=7, color="0.35")
ax[1].set_yscale("log")
ax[1].set_ylim(1e-3, 3e2)
ax[1].set_xlabel("utilisation $u$")
ax[1].set_ylabel(r"$|\eta_s|$ / mV")
ax[1].set_title("Where each electrode operates")
ax[1].legend(fontsize=7, ncol=2, loc="lower right")
fig.tight_layout()
plt.show()


# %% [markdown]
# ## Validation
#
# Ranked before any code was written, and taken in order.
#
# **(2) Internal identities.** The strongest route available, because this is a
# constitutive law and the paper prints three equations that must be consistent.
# Eq. 17 (page 1528), Eq. 16 (page 1528) and Eq. 30 (page 1530) are printed
# separately and not derived from one another; reducing Eq. 17 with Eq. 16
# reproduces Eq. 30 exactly, and the reduction is only possible when
# $\alpha_a+\alpha_c=1$. Four further structural identities are checked above.
#
# **(3) A stated numerical result.** Doyle prints $\delta$ and $\nu$. Their ratio
# eliminates the conductivity, the one input the paper never gives, and becomes a
# direct test of Eq. 30's exchange current density.
#
# **(4) A digitised figure.** Not used, and not needed. Nothing was digitised for
# this page.
#
# Further checks below: the cell must be insensitive to the grid and to the time
# step, and its salt balance must close exactly, since the flux form makes
# conservation a property of the discretisation rather than of the solve.
#
# **Space and time are refined separately.** `Cell.march` sets `dt` and `dt_max`
# from the current alone: they are functions of $I$ and the reference discharge
# time and of *neither* $n_s$ nor $n_c$. So the grid ladder refines space at a
# **constant temporal error**, and the three `grid_*` numbers measure the spatial
# part of the discretisation error only. The `dt_scale` ladder, which multiplies
# `dt` *and* `dt_max`, measures the temporal part beside it, and on this cell it is
# the **larger** of the two.
#
# The salt balance is an identity of the flux form, exact by construction; the one
# thing that can move it is a salt-flux boundary condition, and the break table
# has that row.
#
# One number below is deliberately reported as **not converged**: the peak
# $|\eta_s|$. It sits in the first cathode cell, where the reaction crowds against
# the separator, and **in the very first time step**: the cell prints the step
# index and the time. It is therefore a boundary layer in $t$ as much as in $x$, it
# grows under *both* refinements, and a grid ladder at fixed `dt` cannot bound it.
# So the margin argument is made at the corner where **both** knobs are refined
# together. And the *sustained* peak, the largest $|\eta_s|$ once that initial
# transient has decayed, is printed beside it, because that one is stable across
# the whole refinement box and is the number a reader should carry away.

# %%
def v_at(r_):
    """Cell potential at u = 0.5 of a run() result, computed there (Cell.march u_probe)."""
    return float(r_[2])


GRIDS = [(12, 24), (24, 48), (48, 96), (96, 192)]


def salt_drift(r_):
    """Relative drift of the total salt over a run."""
    return float(abs(r_[:, 3].max() / r_[:, 3].min() - 1))


def ladder(runs_):
    """Spreads over a refinement ladder of run() results: V at u = 0.5 (V), the
    cell-averaged and the peak |eta_s| inside u < 0.8 (mV), salt drift."""
    vs_ = [v_at(x) for x in runs_]
    em_ = [1e3 * x[0][x[0][:, 0] < 0.8, 7].max() for x in runs_]
    ep_ = [1e3 * x[0][x[0][:, 0] < 0.8, 4].max() for x in runs_]
    return {"v": vs_, "em": em_, "ep": ep_, "V": float(max(vs_) - min(vs_)),
            "mean": float((max(em_) - min(em_)) / np.mean(em_)),
            "peak": float(max(ep_) / min(ep_)), "salt": max(salt_drift(x[0]) for x in runs_)}


def grid_study(par=None, **kw):
    """The grid ladder at the production time schedule (par: parameter overrides;
    kw: face_mean, leak)."""
    return ladder([run(I_REF, par, n_s=n_s, n_c=n_c, **kw) for n_s, n_c in GRIDS])


GS = grid_study()
gr = {g: run(I_REF, n_s=g[0], n_c=g[1])[0] for g in GRIDS}
for (n_s, n_c), v_, m_, p_ in zip(GRIDS, GS["v"], GS["em"], GS["ep"]):
    print(f"grid {n_s:3d}+{n_c:3d}: V at u=0.5 = {v_:.6f} V, mean|eta_s| = {m_:.5f} mV, "
          f"peak|eta_s| = {p_:.5f} mV, salt drift below 1e-12 "
          f"(here {salt_drift(gr[(n_s, n_c)]):.0e})")
vs, em, ep = GS["v"], GS["em"], GS["ep"]
GRID_V, GRID_ETA_MEAN, GRID_ETA_PEAK, SALT = GS["V"], GS["mean"], GS["peak"], GS["salt"]
claims.value("grid_V_spread_mV", 1e3 * GRID_V, rtol=0.02, unit="mV",
             note="spatial only: the time schedule is the same on every grid")
claims.value("grid_mean_eta_spread_frac", GRID_ETA_MEAN, rtol=0.02)
claims.value("grid_peak_eta_ratio", GRID_ETA_PEAK, rtol=1e-3)
claims.bound("salt_conservation", SALT, max=1e-12, basis="identity",
             note="telescoping divergence with zero salt flux on both boundary faces")
print(f"\ncell potential spread over an 8x refinement  : {1e3*GRID_V:.3f} mV")
print(f"cell-averaged |eta_s| spread, 8x refinement  : {100*GRID_ETA_MEAN:.2f} %")
print(f"peak |eta_s|, ratio finest/coarsest          : {GRID_ETA_PEAK:.2f} x")
print(f"total salt conserved to below 1e-12 (here {SALT:.1e})")
print()
print("The peak is a boundary-layer value in the first cathode cell and is NOT "
      "converged at\nthese grids. Its successive increments are "
      + ", ".join(f"{d:.3f}" for d in np.diff(ep)) + " mV, i.e. the sequence")
print("is converging but has not arrived; nothing is concluded from it alone.")
THR_1PCT = float(1e3 * eta_linear_threshold(0.01, 0.5, 0.5))
print(f"  finest-grid peak |eta_s|          = {max(ep):.4f} mV")
print(f"  1 % linear-limit threshold        = {THR_1PCT:.2f} mV")
print(f"  margin on THIS ladder alone       = {THR_1PCT/max(ep):.1f} x, against "
      f"{GRID_ETA_PEAK:.2f} x of grid movement")
print()
print("ALL THREE of these spreads are SPATIAL ONLY. march() sets dt and dt_max from")
print("the current alone, so every row above was computed at the same temporal")
print("error, and none of the three says anything about it. The next cell measures")
print("it, and the margin quoted here is re-made there under both knobs at once.")

# %%
# ---- the knob the ladder above holds fixed ----------------------------------
# march() sets dt = 2e-5 t_ref and dt_max = 6e-3 t_ref, with t_ref = t_full/|I|.
# Neither depends on n_s or n_c. dt_scale multiplies BOTH, which is the only way
# to refine a schedule that grows x1.3 per step: scaling the initial step alone
# leaves dt_max, and so the step actually taken over almost the whole discharge,
# exactly where it was.
DT_LEVELS = ((1.0, 2500), (0.5, 2500), (0.25, 3000), (0.125, 5000))
MAX_STEPS = dict(DT_LEVELS)
V_FLOOR = 1e-9 * RTF        # newton tol 1e-9 on the potential scaled by RT/F, V


def dt_run(dts, par=None, **kw):
    """The production-grid run at dt_scale dts (kw: scale_dt_max, newton_tol)."""
    return run(I_REF, par, dt_scale=dts, max_steps=MAX_STEPS[dts], **kw)


def dt_study(par=None, **kw):
    """The dt_scale ladder 1, 1/2, 1/4, 1/8: spreads, observed orders of V(u = 0.5) by
    pymrm.checks.observed_orders, the first-order Richardson temporal error at
    dt_scale = 1, and the finest difference over the solver floor."""
    lad = ladder([dt_run(d, par, **kw) for d, _ in DT_LEVELS])
    oo = observed_orders(lambda n: v_at(dt_run(1.0 / n, par, **kw)), (1, 2, 4, 8))
    v_rich = lad["v"][-1] - (lad["v"][-2] - lad["v"][-1])     # first-order Richardson
    return {**lad, "orders": oo["orders"], "err": float(abs(lad["v"][0] - v_rich)),
            "margin": float(abs(lad["v"][-1] - lad["v"][-2])
                            / (kw.get("newton_tol", 1e-9) * RTF))}


DS = dt_study()
dtr = {d: dt_run(d)[0] for d, _ in DT_LEVELS}
dv_, dm_, dp_ = DS["v"], DS["em"], DS["ep"]
ddif = [1e3 * abs(dv_[i + 1] - dv_[i]) for i in range(len(dv_) - 1)]

print("Time-step refinement at the production grid, I = 10 A/m2")
print(f"{'dt_scale':>9}{'steps':>7}{'V(u=0.5)/V':>13}{'change/mV':>11}{'order':>7}"
      f"{'mean|eta_s|/mV':>16}{'peak|eta_s|/mV':>16}")
for k, (dts, _) in enumerate(DT_LEVELS):
    chg = f"{ddif[k-1]:11.4f}" if k else " " * 11
    order = f"{np.log2(ddif[k-2]/ddif[k-1]):7.2f}" if k >= 2 else " " * 7
    print(f"{dts:9.4f}{len(dtr[dts]):7d}{dv_[k]:13.6f}{chg}{order}"
          f"{dm_[k]:16.5f}{dp_[k]:16.5f}")

DT_ORDER, DT_ORDER_MAX = min(DS["orders"]), max(DS["orders"])
DT_V, DT_ETA_MEAN, DT_ETA_PEAK, DT_V_ERR = DS["V"], DS["mean"], DS["peak"], DS["err"]
claims.bound("dt_observed_order", DT_ORDER, min=0.8, max=1.2, basis="theory",
             note="implicit Euler is first order; the smaller of the two ladder triples")
claims.bound("dt_observed_order_max", DT_ORDER_MAX, min=0.8, max=1.2, basis="theory")
claims.bound("dt_finest_over_roundoff", DS["margin"], min=1e3, basis="theory",
             note="finest V difference over the newton floor 1e-9 RT/F")
claims.value("dt_V_spread_mV", 1e3 * DT_V, rtol=0.02, unit="mV")
claims.value("dt_V_temporal_error_mV", 1e3 * DT_V_ERR, rtol=0.02, unit="mV")
claims.value("dt_mean_eta_spread_frac", DT_ETA_MEAN, rtol=0.02)
claims.value("dt_peak_eta_ratio", DT_ETA_PEAK, rtol=1e-3)
print(f"\nobserved orders {DT_ORDER:.2f} and {DT_ORDER_MAX:.2f} (pymrm.checks.observed_orders), "
      f"first order as implicit Euler")
print(f"gives, so the extrapolation below is the right one to use; the finest difference is "
      f"{DS['margin']:.0e} x the solver floor.")
print(f"Richardson-extrapolated temporal error at the production dt_scale = 1 : "
      f"{1e3*DT_V_ERR:.3f} mV")
print()
print("The same three quantities under each knob, like for like:")
print(f"{'':34}{'8x grid':>11}{'8x dt':>11}{'dt / grid':>12}")
print(f"{'cell potential spread / mV':34}{1e3*GRID_V:11.4f}{1e3*DT_V:11.4f}"
      f"{DT_V/GRID_V:12.2f}")
print(f"{'cell-averaged |eta_s| spread / %':34}{100*GRID_ETA_MEAN:11.3f}"
      f"{100*DT_ETA_MEAN:11.3f}{DT_ETA_MEAN/GRID_ETA_MEAN:12.2f}")
print(f"{'peak |eta_s|, coarsest -> finest':34}{GRID_ETA_PEAK:11.3f}{DT_ETA_PEAK:11.3f}"
      f"{'see below':>12}")
print()
print(f"So on the cell potential the temporal error is {DT_V_ERR/GRID_V:.1f}x the grid "
      f"spread once it is\nextrapolated, and {DT_V/GRID_V:.1f}x it measured the same "
      "way; on the cell-averaged")
print(f"overpotential it is {DT_ETA_MEAN/GRID_ETA_MEAN:.1f}x. On both quantities that "
      "converge, the larger of")
print("the two numerical errors is the temporal one. The peak is the exception and is")
print("not converged in either knob; it is taken apart below.")
print(f"J3.4, the same cell and schedule, commits a dt spread of {J34['dt_spread_mV']:.3f} mV "
      f"(its dt_spread_mV).")
print()
print("Every number on this page is taken at dt_scale = 1, and the argument they")
print(f"support (both electrodes inside the 5 % linear window, {ETA_SEP_DECADES:.1f} decades "
      f"apart) has room of")
print("orders of magnitude, not tenths of a millivolt. The three grid_* metrics are")
print("the spatial half of a two-part error, with the other half measured beside them.")

# %% [markdown]
# #### The peak overpotential is a *first-step* quantity
#
# The one number the two ladders disagree about is the peak. The cell below says
# where it actually sits, refines both knobs together, and separates the initial
# spike from the value the cathode runs at for the rest of the discharge.

# %%
prod = gr[(24, 48)]                  # the production run: default grid, dt_scale = 1
mprod = prod[:, 0] < 0.8


def peak_step(cut=0.0):
    """Step index of the largest |eta_s| of the production run inside u < 0.8 and
    t > cut: a named row of the run, not a crossing."""
    m = mprod & (prod[:, 2] > cut)
    return int(np.flatnonzero(m)[np.argmax(prod[m, 4])])


kpk = peak_step()
print("Where the peak |eta_s| sits in the I = 10 A/m2 discharge:")
print(f"  step {kpk} of {int(mprod.sum())} inside the envelope, "
      f"t = {prod[kpk, 2]:.4g} s of {prod[-1, 2]:.4g} s "
      f"({prod[kpk, 2]/prod[-1, 2]:.1e} of the run), u = {prod[kpk, 0]:.5f}")
claims.exact("peak_eta_at_first_step", kpk == 0)
claims.value("peak_eta_time_s", prod[kpk, 2], rtol=1e-6, unit="s")
claims.value("discharge_time_s", prod[-1, 2], rtol=1e-3, unit="s")
print("  -> it is the FIRST time step. The peak is the initial transient at the")
print("     separator face, and it is resolved by dt as much as by dx. That is why")
print("     it moves under both knobs and is converged under neither.")
print()


# the two knobs together, at the corner of the refinement box
def corner(n_s=96, n_c=192, dts=0.125):
    """Peak |eta_s| inside u < 0.8 with the grid and the time step refined together, mV."""
    r_ = run(I_REF, n_s=n_s, n_c=n_c, dt_scale=dts, max_steps=5000)[0]
    return float(1e3 * r_[r_[:, 0] < 0.8, 4].max())


def corner_margin(eps=0.01):
    """Margin of the corner peak to the linear threshold at tolerance eps."""
    return float(1e3 * eta_linear_threshold(eps, 0.5, 0.5) / corner())


CORNER = run(I_REF, n_s=96, n_c=192, dt_scale=0.125, max_steps=5000)[0]
CORNER_PEAK = corner()
PROD_PEAK = float(1e3 * prod[mprod, 4].max())
claims.value("corner_peak_eta_mV", CORNER_PEAK, rtol=1e-3, unit="mV")
claims.value("corner_peak_margin_factor", corner_margin(), rtol=1e-3)
claims.value("grid_peak_finest_mV", max(ep), rtol=1e-3, unit="mV")
claims.value("dt_peak_finest_mV", dp_[-1], rtol=1e-3, unit="mV")
print(f"{'peak |eta_s| / mV':38}{'value':>10}{'margin to the 1 % threshold':>30}")
for lab, val in [("production grid, production dt", PROD_PEAK),
                 ("8x grid, production dt", max(ep)),
                 ("production grid, 8x dt", dp_[-1]),
                 ("8x grid AND 8x dt together", CORNER_PEAK)]:
    print(f"{lab:38}{val:10.4f}{THR_1PCT/val:30.1f}")
CORNER_MARGIN = float(THR_1PCT / CORNER_PEAK)
print(f"  (1 % linear-limit threshold = {THR_1PCT:.3f} mV)")
print()
print(f"Refining both knobs together moves the peak by a factor "
      f"{CORNER_PEAK/PROD_PEAK:.2f} from the")
print(f"published value and takes the margin from {THR_1PCT/PROD_PEAK:.1f}x to "
      f"{CORNER_MARGIN:.1f}x. It is still an order of")
print("magnitude, so the classification of the cathode as deep inside the linear")
print("window survives, but the honest margin is the smaller one, and it is the")
print("number the caveats quote.")
print()

# the sustained value: what the cathode runs at once the initial spike has gone
SUS_CUT = 10.0                       # s


def sustained(r, cut=SUS_CUT):
    """Largest |eta_s| inside u < 0.8 after time `cut` (cut = 0 keeps the first step), mV."""
    m = (r[:, 0] < 0.8) & (r[:, 2] > cut)
    return float(1e3 * r[m, 4].max())


def sustained_study(cut=SUS_CUT):
    """The sustained peak at the production point and its spread over the refinement box."""
    box = list(gr.values()) + [dtr[d] for d, _ in DT_LEVELS] + [CORNER]
    s_ = [sustained(r, cut) for r in box]
    return {"peak": sustained(prod, cut), "spread": float((max(s_) - min(s_)) / np.mean(s_)),
            "all": s_}


SS = sustained_study()
sus, SUS_PEAK, SUS_SPREAD = SS["all"], SS["peak"], SS["spread"]
claims.value("sustained_peak_eta_mV", SUS_PEAK, rtol=1e-3, unit="mV")
claims.value("sustained_peak_box_spread_frac", SUS_SPREAD, rtol=0.02)
print(f"Sustained peak |eta_s| (t > {SUS_CUT:.0f} s, i.e. past the first-step spike):")
print(f"  at the production point                    : {SUS_PEAK:.4f} mV")
print(f"  over all {len(sus)} runs in the refinement box    : "
      f"{min(sus):.4f} to {max(sus):.4f} mV, spread {100*SUS_SPREAD:.1f} %")
print(f"  margin to the 1 % threshold                : {THR_1PCT/max(sus):.0f} x")
print(f"That is {100*SUS_SPREAD:.1f} % of movement over the whole refinement box, "
      f"against a factor {CORNER_PEAK/PROD_PEAK:.2f}")
print("for the spike, so it is the number to carry away. The first-step spike is")
print(f"real, but it is a transient of the initial condition: {PROD_PEAK/SUS_PEAK:.1f}x "
      f"the sustained\nvalue at the production point and {CORNER_PEAK/max(sus):.1f}x it "
      "at the corner, and resolved by neither knob.")


# %% [markdown]
# #### Salt conservation and the one thing that moves it
#
# `salt_conservation` is a telescoping identity: with zero salt flux on both
# boundary faces the discrete divergence sums to `n[N] - n[0] = 0` for **any**
# interior flux expression and any parameters, so the residual is exact by
# construction. The one error class it *can* see is a wrong salt-flux boundary
# condition; `Cell(leak=...)` puts a salt flux through the collector, as a fraction
# of $|I|/F$, and nothing else.

# %%
print("Salt conservation against a leaking collector boundary condition:")
print(f"{'n_salt[-1] as a fraction of I/F':34}{'steps':>7}{'salt drift':>13}"
      f"{'V(u=0.5)/V':>13}")
salt_break, leak_v = {}, {}
for leak in (0.0, 1e-3, 1e-2):
    rl = run(I_REF, leak=leak)
    salt_break[leak], leak_v[leak] = salt_drift(rl[0]), v_at(rl)
    drift_ = f"{salt_break[leak]:13.1e}" if salt_break[leak] > 1e-12 else f"{'below 1e-12':>13}"
    print(f"{leak:34.0e}{len(rl[0]):7d}{drift_}{leak_v[leak]:13.6f}")
SALT_BREAK = salt_break[1e-3]
claims.value("break_salt_collector_leak", SALT_BREAK, rtol=1e-3,
             note="salt drift with a collector leak of 1e-3 |I|/F")
claims.value("leak_cell_potential_shift_mV", 1e3 * abs(leak_v[1e-3] - leak_v[0.0]), rtol=1e-3,
             unit="mV")
print(f"\nA 0.1 % leak at the collector takes the drift from below 1e-12 to "
      f"{SALT_BREAK:.2e}, and moves the")
print(f"cell potential by {1e3*abs(leak_v[1e-3] - leak_v[0.0]):.2f} mV, more than the "
      "whole numerical error budget")
print("measured above, so this is an error class worth being able to see.")

# %% [markdown]
# #### What `Cell` does not carry
#
# Three differences from `J3.4`'s cell, all latent: nothing on this page depends
# on any of them.
#
# * **`march` has no `record` / `t_stop` snapshot machinery.** On `J3.4` that
#   interface is what makes its two checks with real discriminating power
#   runnable: the short-time semi-infinite comparison and the long-time pore-wall
#   flux. This page's cell checks are the grid and time-step ladders and the
#   salt identity.
# * **There is no charge branch in the stopping rule.** `J3.4` stops a charge at
#   a top voltage; this `march` does not. Every run on this page is a discharge,
#   so nothing here is affected, but a charge would run until `u > u_stop`.
# * **`eps` and `delta_c` come from the parameter dictionary, not constructor
#   arguments.** `J3.4` exposes them; this page passes a modified parameter
#   dictionary instead. *Reuse* below sends a reader to `J3.4` for the reactor.

# %% [markdown]
# ### The `J3.4` open question does not touch any of this
#
# `J3.4` records an unresolved question about the factor $\nu=2$ in the
# concentration term of Doyle's Eq. 5, and deliberately does not take a side. That
# term is in the electrolyte potential equation, not in the kinetics, so it should
# not move the overpotentials. "Should not" is not a check, so the cell is run both
# ways.

# %%
def nu_sens(nu_alt=1.0):
    """|peak cathode |eta_s| with nu_factor 2 / with nu_alt - 1| (nu_alt = 2 is no change)."""
    p2 = run(I_REF)[0]
    pa = run(I_REF, nu_factor=nu_alt)[0] if nu_alt != 2.0 else p2
    return float(abs(p2[p2[:, 0] < 0.8, 4].max() / pa[pa[:, 0] < 0.8, 4].max() - 1))


nu_tab = pd.DataFrame([{"nu factor": nu_f, "V at u = 0.5 / V": v_at(run(I_REF, **kw)),
                        "peak |eta_s| cathode / mV": 1e3 * (lambda r: r[r[:, 0] < 0.8, 4].max())(run(I_REF, **kw)[0]),
                        "eta_s anode / mV": 1e3 * (lambda r: r[r[:, 0] < 0.8, 6].max())(run(I_REF, **kw)[0])}
                       for nu_f, kw in ((2.0, {}), (1.0, {"nu_factor": 1.0}))])
print(nu_tab.to_string(index=False, float_format=lambda v: f"{v:12.5f}"))
NU_ETA = nu_sens()
NU_V = float(abs(nu_tab["V at u = 0.5 / V"].iloc[0] - nu_tab["V at u = 0.5 / V"].iloc[1]))
claims.value("nu_factor_eta_sensitivity_frac", NU_ETA, rtol=1e-3)
claims.value("nu_factor_V_shift_mV", 1e3 * NU_V, rtol=1e-3, unit="mV")
print(f"\nthe nu = 2 question moves the cell potential by {1e3*NU_V:.1f} mV "
      f"and the peak overpotential by {100*NU_ETA:.2f} %.")
print("The regime classification of both electrodes is unchanged either way, so "
      "nothing on this page depends on how that question is settled.")

# %% [markdown]
# ### The defect-injection table
#
# Every claim the page registers has a row here that fails its own check, apart
# from the structural identities, which are named as such where they are
# registered. Each row injects one defect through an argument of the study that
# computes the claim, or through Doyle's parameter dictionary for every run that
# uses it, and recomputes the claim; a study that raises records nan for the
# claims it targets and stays in the table. The row "injected defects removed"
# sets every studied defect or sensitivity to nothing, so those values collapse:
# it shows those numbers are the defects' effects and nothing else.

# %%
def safe_metrics(fn):
    """Evaluate a study once for several claims of one row; an exception gives nan."""
    try:
        return fn()
    except Exception as exc:                     # noqa: BLE001 - a failed study is a result
        print(f"  study raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
        return None


_ONCE = {}


def once(key, fn):
    """A study shared by several claims of one row, evaluated once."""
    if key not in _ONCE:
        _ONCE[key] = safe_metrics(fn)
    return _ONCE[key]


def pick(key, fn, item):
    """One item of a shared study, nan when the study raised."""
    out = once(key, fn)
    return np.nan if out is None else out[item]


K2_HALF = {"k_2": P["k_2"] / 2}
I01_HALF = {"i_0_1": P["i_0_1"] / 2}
D_HALF = {"D": P["D"] / 2}
i0_k2 = lambda: superficial(K2_HALF)                     # noqa: E731
env_k2 = lambda: envelope(10.0, K2_HALF)                 # noqa: E731
env_i01 = {I: (lambda I=I: envelope(I, I01_HALF)) for I in (10.0, 20.0)}
env_D = {I: (lambda I=I: envelope(I, D_HALF)) for I in (10.0, 20.0)}
thr_i01 = lambda e: 1e3 * eta_linear_threshold(e, 0.5, 0.5)    # noqa: E731
DEFECT_FREE = defect_metrics(DEFECTS_REMOVED)

ROWS = [
    ("Eq. 30 tail exponents swapped", {
        "eq17_vs_eq30_max_rel": lambda: max(d_eq17_vs_eq30(aa, ac, e30tail=(ac, aa))
                                            for aa, ac in SUM1_PAIRS),
        "i0_from_partials_vs_eq30_rel": lambda: partials(e30tail=True)[1]}),
    ("alpha_a + alpha_c != 1 in the zero-current check", {
        "eq17_zero_current_vs_eq16_V": lambda: eq16_dev(NOT1_PAIRS)}),
    ("(1/n - 1) factor dropped from the closed-form offset", {
        "alpha_offset_closed_form_V": lambda: max(offset_check(aa, ac, drop_factor=True)
                                                  for aa, ac in OFF_PAIRS)}),
    ("1/n dropped from the equilibrium potential", {
        "zero_current_at_zero_eta_rel": lambda: zero_current(drop_n=True)}),
    ("n forced to 1 in R_ct", {"Rct_identity_max_rel": lambda: rct_check(n_forced=1.0)}),
    ("arcsinh inversion used at alpha_a = 0.3", {
        "arcsinh_inversion_max_rel": lambda: asinh_check(0.3)}),
    ("ln(1/eps) written for ln(1 + 1/eps)", {
        "tafel_threshold_closed_vs_numeric": lambda: tafel_check("log")}),
    ("alpha swap omitted from the mirror", {"linear_sign_convention_rel": lambda: mirror(False)}),
    ("general i_0 compared with Eq. 30 at alpha = (0.5, 0.4)", {
        "i0_general_vs_eq30_rel": lambda: i0_gen_check(0.5, 0.4)}),
    ("T read as 298.15 K", {
        "eta_tafel_10pct_mV": lambda: thresholds(RTF_298)["tafel10"],
        "eta_tafel_1pct_mV": lambda: thresholds(RTF_298)["tafel1"],
        "eta_linear_5pct_mV": lambda: thresholds(RTF_298)["lin5"],
        "eta_linear_1pct_mV": lambda: thresholds(RTF_298)["lin1"],
        "linear_5pct_alpha0p42_mV": lambda: 1e3 * eta_linear_threshold(0.05, 0.42, 0.58, +1, RTF_298)}),
    ("sqrt(24 eps) written as sqrt(12 eps)", {"linear_sym_asymptote_rel": lambda: sym_asymptote(12.0)}),
    ("asymmetric asymptote applied at alpha_a = 0.42", {
        "linear_asym_asymptote_rel": lambda: asym_asymptote((0.1, 0.2, 0.3, 0.42))}),
    ("linear threshold taken at the last crossing", {
        "linear_window_narrowing_0p5_to_0p1": lambda: narrowing("last"),
        "linear_5pct_alpha0p3_mV": lambda: 1e3 * eta_linear_threshold(0.05, 0.3, 0.7,
                                                                       crossing="last")}),
    ("symmetric slope 2 alpha_a f eta at asymmetric alpha", {
        "linear_jump_alpha_a_1pct": lambda: jump(0.01, slope=lambda a: 2 * a),
        "linear_jump_alpha_a_5pct": lambda: jump(0.05, slope=lambda a: 2 * a),
        "linear_overshoot_peak_alpha0p45": lambda: first_lobe_peak(0.45, 0.9)}),
    ("initial solid state read as u_0 = 0.02", {
        "nu2_over_delta_dev_frac": lambda: nu2_study(0.02)["dev"],
        "nu2_over_delta_as_kappa_ratio": lambda: nu2_study(0.02)["as_kappa"],
        "j34_kappa_route_spread": lambda: nu2_study(0.02)["kappa_spread"],
        "j34_kappa_spread_vs_committed_rel": lambda: abs(100 * nu2_study(0.02)["kappa_spread"]
                                                         / J34["kappa_route_spread_pct"] - 1.0),
        "j35_C_r_scale_ratio": lambda: nu2_study(0.02)["scale_ratio"]}),
    ("k_2 halved", {
        "superficial_reading_linear_err": lambda: pick("k2s", i0_k2, "sup"),
        "interfacial_reading_linear_err": lambda: pick("k2s", i0_k2, "int"),
        "superficial_tafel_headroom_factor": lambda: pick("k2s", i0_k2, "headroom"),
        "superficial_eta_s_mV": lambda: 1e3 * pick("k2s", i0_k2, "e_sup"),
        "i0_initial_A_m2": lambda: pick("k2s", i0_k2, "i0"),
        "superficial_i_over_i0": lambda: I_REF / pick("k2s", i0_k2, "i0"),
        "interfacial_i_over_i0": lambda: I_REF / (a_spec * dc) / pick("k2s", i0_k2, "i0"),
        "tafel_cathode_predicted_eta_mV": lambda: 1e3 * RTF / alpha * np.log(
            pick("k2s", i0_k2, "i0") / (I_REF / (a_spec * dc))),
        "cell_cathode_mean_eta_s_mV": lambda: 1e3 * pick("k2e", env_k2, "mean"),
        "cell_cathode_peak_eta_s_mV": lambda: 1e3 * pick("k2e", env_k2, "peak"),
        "eta_separation_decades": lambda: separation(cath_mean10=pick("k2e", env_k2, "mean"))}),
    ("i_o1 of the foil halved", {
        "cell_anode_eta_s_I10_mV": lambda: 1e3 * pick("i10", env_i01[10.0], "anode"),
        "cell_anode_eta_s_I20_mV": lambda: 1e3 * pick("i20", env_i01[20.0], "anode"),
        "foil_I10_over_linear_5pct": lambda: 1e3 * pick("i10", env_i01[10.0], "anode") / thr_i01(0.05),
        "foil_I10_over_linear_1pct": lambda: 1e3 * pick("i10", env_i01[10.0], "anode") / thr_i01(0.01),
        "foil_I20_over_linear_5pct": lambda: 1e3 * pick("i20", env_i01[20.0], "anode") / thr_i01(0.05),
        "foil_I20_over_linear_1pct": lambda: 1e3 * pick("i20", env_i01[20.0], "anode") / thr_i01(0.01),
        "linear_anode_dV_I20_mV": lambda: 1e3 * abs(anode_shift(20.0, "linear", i01=P["i_0_1"] / 2)),
        "linear_anode_dV_I10_mV": lambda: 1e3 * abs(anode_shift(10.0, "linear", i01=P["i_0_1"] / 2)),
        "eta_separation_decades": lambda: separation(an10=pick("i10", env_i01[10.0], "anode"))}),
    ("electrolyte diffusivity halved", {
        "u_cutoff_I20": lambda: pick("D20", env_D[20.0], "u_cut"),
        "u_1p9V_I20": lambda: pick("D20", env_D[20.0], "u19"),
        "cell_cathode_peak_I20_cutoff_mV": lambda: 1e3 * pick("D20", env_D[20.0], "peak"),
        "u_cutoff_I10": lambda: pick("D10", env_D[10.0], "u_cut"),
        "u_cutoff_I10_vs_J34_abs": lambda: abs(pick("D10", env_D[10.0], "u_cut")
                                               - J34["u_at_cutoff_I10"]),
        "discharge_time_s": lambda: run(I_REF, D_HALF)[0][-1, 2]}),
    ("specific area from R_s doubled", {
        "a_delta_c": lambda: 3.0 * (1 - eps_c) / (2 * P["R_s"]) * dc}),
    ("Tafel cathode law in the standard-form comparison", {
        "bv_vs_eq17_cell_potential_mV": lambda: 1e3 * dv("tafel")}),
    ("linear law switched at both electrodes", {
        "linear_cathode_dV_mV": lambda: 1e3 * dv("linear", anode_law="linear"),
        "linear_cathode_dV_noreplay_mV": lambda: 1e3 * dv("linear", replay=False,
                                                          anode_law="linear")}),
    ("limit laws replaced by the exact law", {
        "tafel_cathode_dV_mV": lambda: 1e3 * dv("bv"),
        "linear_anode_dV_I20_mV": lambda: 1e3 * abs(anode_shift(20.0, "linear", "linear")),
        "linear_anode_dV_I10_mV": lambda: 1e3 * abs(anode_shift(10.0, "linear", "linear"))}),
    ("Tafel prediction tested on the exact-law run", {
        "tafel_first_step_vs_prediction_frac": lambda: tafel_pred("bv")}),
    ("harmonic face means replaced by arithmetic", {
        "grid_V_spread_mV": lambda: 1e3 * pick("ar", lambda: grid_study(face_mean="arithmetic"), "V"),
        "grid_mean_eta_spread_frac": lambda: pick("ar", lambda: grid_study(face_mean="arithmetic"), "mean"),
        "grid_peak_eta_ratio": lambda: pick("ar", lambda: grid_study(face_mean="arithmetic"), "peak"),
        "grid_peak_finest_mV": lambda: max(pick("ar", lambda: grid_study(face_mean="arithmetic"), "ep"))}),
    ("collector salt leak of 1e-3 |I|/F", {
        "salt_conservation": lambda: salt_drift(run(I_REF, leak=1e-3)[0])}),
    ("dt_max left unscaled in the dt ladder", {
        name: (lambda item=item, f=f: f(pick("dtm", lambda: dt_study(scale_dt_max=False), item)))
        for name, item, f in [
            ("dt_V_spread_mV", "V", lambda v: 1e3 * v), ("dt_V_temporal_error_mV", "err", lambda v: 1e3 * v),
            ("dt_observed_order", "orders", min), ("dt_observed_order_max", "orders", max),
            ("dt_finest_over_roundoff", "margin", float), ("dt_mean_eta_spread_frac", "mean", float),
            ("dt_peak_eta_ratio", "peak", float), ("dt_peak_finest_mV", "ep", lambda v: v[-1])]}),
    ("k_2 halved in the grid ladder", {
        "grid_mean_eta_spread_frac": lambda: pick("gk2", lambda: grid_study(K2_HALF), "mean"),
        "grid_peak_eta_ratio": lambda: pick("gk2", lambda: grid_study(K2_HALF), "peak"),
        "grid_peak_finest_mV": lambda: max(pick("gk2", lambda: grid_study(K2_HALF), "ep"))}),
    ("k_2 halved in the dt ladder", {
        "dt_mean_eta_spread_frac": lambda: pick("dk2", lambda: dt_study(K2_HALF), "mean"),
        "dt_peak_eta_ratio": lambda: pick("dk2", lambda: dt_study(K2_HALF), "peak"),
        "dt_peak_finest_mV": lambda: pick("dk2", lambda: dt_study(K2_HALF), "ep")[-1],
        "dt_V_spread_mV": lambda: 1e3 * pick("dk2", lambda: dt_study(K2_HALF), "V")}),
    ("newton tolerance 1e-4 in the dt ladder", {
        "dt_finest_over_roundoff": lambda: pick("dnt", lambda: dt_study(newton_tol=1e-4), "margin")}),
    ("corner taken at 4x instead of 8x", {"corner_peak_eta_mV": lambda: corner(48, 96, 0.25)}),
    ("margin taken to the 5 % threshold", {"corner_peak_margin_factor": lambda: corner_margin(0.05)}),
    ("first-step transient kept in the sustained peak", {
        "sustained_peak_eta_mV": lambda: sustained_study(0.0)["peak"],
        "sustained_peak_box_spread_frac": lambda: sustained_study(0.0)["spread"],
        "peak_eta_at_first_step": lambda: peak_step(SUS_CUT) == 0}),
    ("production time step halved", {"peak_eta_time_s": lambda: dt_run(0.5)[0][kpk, 2]}),
    ("injected defects removed", {
        "defect_eq30_tail_swap_eq17_vs_eq30": lambda: DEFECT_FREE["e30_tail"],
        "defect_eq17_prefactor_swap_eq17_vs_eq30": lambda: DEFECT_FREE["e17_pre"],
        "defect_eq30_tail_swap_i0_from_partials": lambda: DEFECT_FREE["i0_asym"],
        "break_linear_threshold_last_crossing_mV": lambda: 1e3 * eta_linear_threshold(0.05, 0.3, 0.7),
        "break_salt_collector_leak": lambda: salt_drift(run(I_REF, leak=0.0)[0]),
        "leak_cell_potential_shift_mV": lambda: 1e3 * abs(v_at(run(I_REF, leak=0.0)) - leak_v[0.0]),
        "nu_factor_eta_sensitivity_frac": lambda: nu_sens(2.0),
        "nu_factor_V_shift_mV": lambda: 1e3 * abs(v_at(run(I_REF)) - v_at(run(I_REF, nu_factor=2.0)))}),
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
        elif np.isfinite(float(out)) and abs(float(out)) < 1e-9:   # round-off: a bound
            out = "below 1e-9"
        BREAK.append((row, name, base_, out, caught_))

bt = pd.DataFrame(BREAK, columns=["injected defect", "metric", "baseline", "with defect", "caught"])
pd.set_option("display.max_colwidth", 60)
pd.set_option("display.max_rows", 300)
print(bt.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
idle = [r for r in bt["injected defect"].unique()
        if not bt[(bt["injected defect"] == r) & bt.caught].shape[0]]
print(f"\n{len(ROWS)} rows defined, {bt['injected defect'].nunique()} in the table; "
      f"rows that fail no claim: {idle or 'none'}")
missing = [k for k in claims.metrics if not claims.metrics[k].get("structural")
           and not bt[(bt.metric == k) & bt.caught].shape[0]]
print(f"claims no row fails: {missing or 'none'}")

# %% [markdown]
# ## What pymrm adds
#
# **Nothing to the law.** Butler–Volmer is one line of arithmetic. The
# identities, the closed-form Tafel threshold and the validity map need numpy and a
# root-finder; pymrm contributes only `pymrm.checks.find_roots`, which locates every
# threshold and every zero-current potential with a uniqueness check (or, for the
# linear threshold, keeps the first of several crossings it finds).
#
# What pymrm supplies is the **operating envelope**, the answer to *where does a
# real electrode actually sit*, which is the only reason the thresholds are
# actionable. That needs a coupled porous-electrode solve with electroneutrality:
# `J3.4`'s cell, built from `construct_grad` and `construct_div` with the
# boundary fluxes written on the outward-normal faces, harmonic face means at the
# porosity jump, and one `newton` solve per implicit-Euler step on a
# `NumJac((N, 3))` Jacobian, each step checked for success. The time-step ladder
# gets its orders from `pymrm.checks.observed_orders`. A cell borrowed from another
# page carries that page's evidence only for that page's physics, so the grid and
# time-step ladders and the salt identity are measured here again, each with a
# break row.
#
# That is the reusable pattern, and it is more interesting than the page: **a
# constitutive law and a reactor model in the same gallery let each check the
# other.** The law's limits say what the cell may safely assume; the cell says which
# part of the law's domain it lives in. Neither statement is available from one of
# them alone.
#
# Two things this page adds that are not in `J3.4` or in `J3.5`:
#
# * Eq. 30 is *derived* from Eq. 17 rather than used alongside it, which turns two
#   separately-transcribed equations into a check on both, and shows that Doyle's
#   Eqs. 16, 17 and 30 are mutually consistent only under $\alpha_a+\alpha_c=1$.
#   (`J3.5`'s `bv_invert` already *contains* this algebra without saying so; what is
#   new is stating it, testing it, and stating the condition it needs.)
# * the $\nu^2/\delta$ ratio, which tests the exchange current density with the
#   unprinted conductivity divided out. This is **not** a second, independent
#   disagreement: it is arithmetically identical to the spread between `J3.4`'s two
#   routes to $\kappa$, as the cell in Check 3 demonstrates. What the reading adds
#   is that the disagreement survives with $\kappa$ removed entirely, so it cannot
#   be attributed to the conductivity reconstruction.

# %% [markdown]
# ## Reuse
#
# **The thresholds transfer; the exchange current density does not.** Everything in
# the validity table is a property of the functional form: scale by $RT/F$ at your
# temperature, divide by $\alpha_a+\alpha_c$, and the Tafel numbers are exact. The
# linear ones need the root-find, which is three lines. No electrode property enters
# either.

# %%
display(Markdown(f"""
**What decides your regime is $i/i_0$, and $i$ means the *interfacial* current
density.** The single most common way to get this wrong is to use the superficial
current density of a porous electrode. This cell's factor between the two is
$a\\delta_c={A_DELTA_C:.0f}$, and it is worth being precise about what that factor
does and does not buy, because the obvious dramatic statement is false. **It does
not rescue this cathode from the Tafel regime, nothing does, because it was never
near it.** Even on the *superficial* current the cathode sits at
$i/i_0={I_REF/i0_init:.2f}$, which is $\\eta_s={1e3*SUPER['e_sup']:.2f}$ mV, a
{100*SUP_LIN_ERR:.2f} % linear error and already inside the 1 % linear window;
reaching the 5 % Tafel boundary would take {SUP_TAFEL_HEADROOM:.0f} times more
current. What $a\\delta_c$ buys is depth, not rescue: it moves the cathode from a
{100*SUP_LIN_ERR:.2f} % linear error to a {100*INT_LIN_ERR:.0e} % one. All of that is
computed in *Results*, under "Superficial or interfacial?".

The rule still matters, just not for the reason a dramatic version of it would
give. It matters at high rate: at {SUP_TAFEL_HEADROOM:.0f} times this current the
superficial reading would put the electrode at the Tafel boundary while its
interface was still deep in the linear window, and the two readings would then
disagree about which law to use. Compute $a$ before deciding anything; and if the
factor turns out not to change your verdict, say so rather than dressing it up.
"""))

# %% [markdown]
# **Take Eq. 30's structure, not its constants.** The composition dependence
# $(c_{\max}-c)^{\alpha_c}c^{\alpha_a}(c_T-c_s)^{\alpha_a}c_s^{\alpha_c}$ is what
# site-limited kinetics gives on both sides of the interface. Marquis et al. print
# the same structure for a liquid electrolyte with the $(c_{\max}-c)$ factor
# dropped, which is exactly the modification Doyle describes in words. The value of
# $k_2$ is not transferable at all: Doyle set it, in his own words, to "a value
# corresponding to a nearly reversible situation", because no exchange current
# density for TiS₂ had been reported. That choice is *why* this cathode is in the
# linear regime, and a measured $k_2$ would move it.
#
# **If you need the reactor, start from `J3.4`, not from here.** The `Cell` class on
# this page is a port with one switch added, and `J3.4` carries the provenance, the
# conductivity reconstruction and the digitised comparison. `J3.5` shows how the
# same cell reduces asymptotically, and its `bv_invert` is a ready-made closed-form
# inversion of Eq. 17, worth knowing about, but note that it is exact only under
# the $\alpha_a+\alpha_c=1$ condition derived here, which that page does not state.
# This page is the constitutive layer underneath both.
#
# **Two numbers in the gallery for one idea.** If you compare the interfacial and
# superficial current densities, `J3.5` prints the same ratio as `C_r`, smaller
# than this page's interfacial $i/i_0$ (both printed in *Results*). They are the
# same quantity divided by different exchange current densities: `J3.5` follows
# Marquis in scaling on a
# *typical* value with the solid factors replaced by $c_T$, this page evaluates
# Eq. 30 at Doyle's stated initial state, and the ratio between them is exactly
# $\sqrt{u_0(1-u_0)}$, computed in *Results*. Use whichever scale your reference
# uses, and do not compare the two numbers directly.
#
# **What not to reuse.** The absolute overpotentials. They belong to one cell, at
# 100 °C, with a rate constant its authors chose for convenience rather than
# measured. The *method*, solve the reactor, read out $\eta_s$, compare with the
# closed-form thresholds, is the part that transfers.

# %% [markdown]
# ### Claims and the prose
#
# Every number the README states is checked below against the claim it quotes, so
# a number in the text cannot drift from the computation and a cut in its
# precision shows in CI (agreement.json "quoted"). Round-off-level results are
# quoted as bounds.

# %%
SAYS = {
    "eq17_vs_eq30_max_rel": ["below 1e-10"], "eq17_zero_current_vs_eq16_V": ["below 1e-13"],
    "zero_current_at_zero_eta_rel": ["below 1e-12"],
    "partial_currents_balance_rel_structural": ["below 1e-12"],
    "i0_from_partials_vs_eq30_rel": ["below 1e-12"], "Rct_identity_max_rel": ["below 1e-7"],
    "arcsinh_inversion_max_rel": ["below 1e-8"], "tafel_threshold_closed_vs_numeric": ["below 1e-12"],
    "nu2_over_delta_dev_frac": ["-15.7 %"], "j34_kappa_route_spread": ["18.7 %"],
    "bv_vs_eq17_cell_potential_mV": ["below 1e-6"], "linear_cathode_dV_mV": ["below 1e-6"],
    "tafel_first_step_vs_prediction_frac": ["0.46 %"], "grid_V_spread_mV": ["0.149"],
    "dt_V_temporal_error_mV": ["0.273"], "dt_observed_order": ["at least 0.8"], "dt_observed_order_max": ["at most 1.2"],
    "salt_conservation": ["below 1e-12"], "break_salt_collector_leak": ["0.022"],
    "eta_tafel_10pct_mV": ["77.1"], "eta_linear_5pct_mV": ["35.9"], "eta_tafel_1pct_mV": ["148.4"],
    "eta_linear_1pct_mV": ["15.8"], "linear_5pct_alpha0p3_mV": ["8.29"],
    "break_linear_threshold_last_crossing_mV": ["203.4"],
    "cell_cathode_mean_eta_s_mV": ["0.032"], "cell_cathode_peak_eta_s_mV": ["0.74"],
    "cell_anode_eta_s_I10_mV": ["24.8"], "cell_anode_eta_s_I20_mV": ["46.6"],
    "foil_I10_over_linear_5pct": ["0.69"], "foil_I10_over_linear_1pct": ["1.57"],
    "foil_I20_over_linear_5pct": ["1.30"], "foil_I20_over_linear_1pct": ["2.95"],
    "corner_peak_eta_mV": ["1.519"], "corner_peak_margin_factor": ["10.4"],
    "sustained_peak_eta_mV": ["0.146"], "sustained_peak_box_spread_frac": ["10.6 %"],
    "interfacial_i_over_i0": ["0.0010"], "eta_separation_decades": ["2.9"],
    "tafel_cathode_dV_mV": ["519"], "u_cutoff_I20": ["0.3707"], "u_1p9V_I20": ["0.2638"],
    "cell_cathode_peak_I20_cutoff_mV": ["1.679"], "u_cutoff_I10": ["0.8312"],
    "defect_eq17_prefactor_swap_eq17_vs_eq30": ["0.95"],
    "defect_eq30_tail_swap_eq17_vs_eq30": ["0.97"], "peak_eta_time_s": ["0.39"],
    "discharge_time_s": ["16862"], "j35_C_r_scale_ratio": ["0.0995"],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)
print(f"prose literals: {sum(len(v) for v in SAYS.values())} checked against their claims with claims.says")
claims.report();   # the semicolon keeps the notebook from echoing the payload
