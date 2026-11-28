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
# title: "Michaelis-Menten and Briggs-Haldane: one rate law, two constants, and what data can decide between them"
# description: "Froment, De Wilde & Bischoff section 1.5.1 derives the same rate law twice - by rapid equilibrium, giving K_M = k_-1/k_1 (eq. 1.5.1-8), and by the pseudo steady state, giving K_M = (k_-1+k_2)/k_1 (eq. 1.5.1-13/14). The two differ by the factor 1 + k_2/k_-1, which can be anything. This page asks what an experiment could do about that, and the answer is bleak in a specific, quantifiable way: the steady-state rate law is EXACTLY independent of k_2/k_-1 - proved symbolically, and identical to double precision over eight decades of it - so no rate-versus-substrate dataset of any size or precision can separate the two readings of K_M; a batch transient does see the ratio, but only at order eps = C_E0/(c_A0+K_M), and most of that signal lives in the INDUCTION LAYER, whose peak is at t = 2.0e-7 min here, which is what Froment's own sentence about stopped-flow experiments is about and which no table sampled in minutes can reach; over the window the source table actually samples the spread is 0.242 eps, and with the other three constants free to compensate, telling k_2/k_-1 = 1e-4 from 1e4 at that table's own 5e-4 cmol/L print resolution needs eps = 5.80e-2 - a factor 18.2 more enzyme than the SAME statistic on the SAME grid with nothing free to compensate, and a loading at which the quasi-steady-state approximation that produces the Michaelis-Menten form is already wrong by 88.9 times that resolution. Ten symbolic identities close the printed chain; the seven-point batch table of Rawlings & Ekerdt Exercise 9.15 gives r_m = 0.162583 cmol/L/min and K_M = 0.503055 cmol/L at 4.02e-3 RMS against 3.90e-2 for a first-order null, and adding the two elementary-step parameters to that fit is not significant (F = 1.90 against F_crit,95 = 9.55). Two printed defects, reported not repaired: p. 24's high-concentration limit r_A = k_1 C_E^0, which the book's own eq. (1.5.1-15) and Problem 2.4(b) both contradict and which is dimensionally impossible; and p. 25's 'Lineweaver-Burke plot of 1/r versus C_A', where 1/r is strictly convex in C_A and exactly linear in 1/C_A. No straight line at all comes closer than 31.7 % of the curve's own range over 0.05-1.0 cmol/L - that is the Chebyshev floor, closed form, convention-free; the best LEAST-SQUARES line is further out still, at 50.3 % under the page's log-weighted convention and 66.9 % under a uniform one."
# categories: [sec:J, struct:S1, tier:T0, data:tier2, phase:liquid]
# date: 2026-08-13
# ---
#
# # Michaelis-Menten and Briggs-Haldane: one rate law, two constants, and what data can decide between them
#
# **Catalog ID:** `J4.6` · **Structures:** `S1` · **Tier:** T0
#
# Section 1.5.1 of Froment, De Wilde & Bischoff derives the enzyme rate law
# **twice**. The first route assumes the complex decomposition is rate
# determining, so the binding step reaches equilibrium; it ends at
#
# $$r_{\rm A} = r_{\rm P} = \frac{k_2 C_{\rm A} C_{\rm E}^0}{K_M + C_{\rm A}},
# \qquad K_M = \frac{k_{-1}}{k_1}, \tag{1.5.1-8}$$
#
# which the book's next line names *"the Michaelis-Menten equation for the rate of
# a simple enzymatic reaction"*. The second route - *"formulated by Briggs and
# Haldane"* - assumes only that the complex is at a pseudo steady state, and ends
# at the identical expression with
#
# $$K_M = \frac{k_{-1} + k_2}{k_1}. \tag{1.5.1-14}$$
#
# Same functional form. Different constant. The ratio between them is
# $1 + k_2/k_{-1}$, and nothing in the rate law bounds it.
#
# **So what can an experiment settle?** That is the question this page is for, and
# it has to be asked carefully, because the tempting answer - "fit a
# Michaelis-Menten curve, and if it fits, the mechanism is confirmed" - is wrong
# in both directions. The page answers it in three steps.
#
# **One: no amount of steady-state rate data can separate the two readings.**
# Write $\rho = k_2/k_{-1}$. **The result is algebraic, and it is proved that way
# here**: substituting the reparameterisation
# $k_2 = r_m/C_{\rm E}^0$, $k_{-1} = k_2/\rho$, $k_1 = (k_{-1}+k_2)/K_M$ into
# eq. (1.5.1-13) returns eq. (1.5.1-15) identically and
# $\partial r/\partial\rho \equiv 0$, in sympy, in the symbolic cell under *The published model*. The
# numerical sweep is a *check on the code*, not the finding: holding
# $r_m = k_2 C_{\rm E}^0$ and $K_M = (k_{-1}+k_2)/k_1$ fixed and sweeping $\rho$
# over **eight decades**, $r(C_{\rm A})$ evaluated from the elementary constants is
# identical to double precision over four decades of concentration - **maximum
# relative difference below 1e-14** on the nine swept $\rho$ at
# $C_{\rm E}^0 = 10^{-3}$ (where it comes out exactly zero, because the
# floating-point round trip through the elementary constants happens to be exact
# for those values), and **below 1e-14** over 2000 random $(\rho, C_{\rm E}^0)$
# draws, where that round trip is no longer exact and the difference is a few ulp.
# The exact zero is a property of those particular values; the *degeneracy* is a
# property of the algebra. The map from four elementary constants
# $(k_1, k_{-1}, k_2, C_{\rm E}^0)$ to the two observable ones $(r_m, K_M)$ has a
# **two-dimensional fibre**, and $\rho$ runs along it. Meanwhile
# $K_M^{\rm BH}/K_M^{\rm MM} = 1+\rho$ runs from 1.0001 to $10^4$ across that same
# sweep. The book says this itself, in one sentence on p. 25 that is easy to read
# past: *"If, for more insight into the process, the rate coefficients $k_1$,
# $k_{-1}$ and $k_2$ of the elementary steps themselves are required, only
# transient experimentation (stopped flow, or relaxation ...) can help."*
#
# **Two: a batch transient does see $\rho$ - at order $\varepsilon$, and mostly
# where nobody is looking.** With $(r_m, K_M, C_{\rm E}^0)$ fixed, the full
# mass-action solution $C_{\rm A}(t)$ moves across the same eight decades of
# $\rho$ by $0.977\,\varepsilon$, where
# $\varepsilon = C_{\rm E}^0/(c_{\rm A0}+K_M)$ - the group that organises the
# quasi-steady-state error, derived here rather than taken from anywhere. But
# **that maximum sits at $t = 2.0\times10^{-7}$ min, inside the induction layer**,
# whose duration is $1/[k_1(c_{\rm A0}+K_M)]$ and therefore differs by the same
# four decades: 2.1e-8 min at $\rho = 10^{-4}$ against 2.1e-4 min at
# $\rho = 10^{4}$. A uniform 181-point scan of the 18-minute run misses it by
# **66 %**, which is why the maximum here is root-found on a logarithmic scan and
# not sampled. **Froment's one sentence about stopped-flow and relaxation
# experiments is exactly right, and this is the number behind it.**
#
# Over the window an experiment sampled like the source table can actually see -
# $t \ge 3$ min, its first sample after zero - what is left is a persistent offset
# of $0.242\,\varepsilon$. Root-found against that table's 5e-4 cmol/L print
# resolution it needs $\varepsilon = 2.07\times10^{-3}$, and even that holds
# every other constant fixed, which no real fit does. **Let the other three float
# and the requirement moves to $\varepsilon = 5.80\times10^{-2}$** - at which
# loading the quasi-steady-state approximation, the thing that makes the rate law
# Michaelis-Menten in the first place, is already off by **88.9 times the
# resolution at which the discrimination becomes possible.** The regime where the
# two derivations differ observably in a slow measurement is the regime where
# neither of them is the right model.
#
# **How much of that is compensation, measured like for like.** The two
# $\varepsilon$'s above are *not* the same measurement, so the page decomposes the
# ratio between them instead of charging all of it to compensation. The first,
# $2.07\times10^{-3}$, is a **maximum over $t\in[3,18]$ min** of the $\rho$-spread;
# the second, $5.80\times10^{-2}$, is an **RMS over a 25-point grid starting at
# $t = 0$** of a *refitted* misfit. Repeat the first measurement with the second's
# statistic and grid, still with nothing free to compensate, and the threshold is
# $\varepsilon = 3.18\times10^{-3}$. So of the raw factor **28.0** between the two
# printed numbers, **1.54 is the change of statistic and window** and only
# **18.2 is the price of letting $r_m$, $K_M$ and $C_{\rm E}^0$ compensate.**
# 18.2 is the number this page means by the cost of compensation.
#
# **Three: on the actual data, the elementary constants are not merely poorly
# determined - they are free.** Rawlings & Ekerdt's Exercise 9.15 prints seven
# $c_{\rm S}$-versus-time points. They give $r_m = 0.162583$ cmol/L/min and
# $K_M = 0.503055$ cmol/L at 4.02e-3 RMS, against 3.90e-2 for the best first-order
# null - a real result, and the page's only empirical one. Pin $\rho$ anywhere from
# $10^{-4}$ to $10^{4}$, let $r_m$, $K_M$ and $C_{\rm E}^0$ re-optimise, and the
# achievable sum of squares changes by **8.6 %** across those eight decades,
# monotonically, with no interior optimum. Adding both elementary-step parameters
# to the two-parameter fit is not significant on any account: **F = 1.90 against
# F(2,3) = 9.55 at 95 %**.
#
# **What this page therefore does NOT claim.** It does not validate "the
# Michaelis-Menten mechanism". It establishes that a two-parameter saturating rate
# law describes one seven-point batch run far better than three simple-order
# alternatives, and that this is compatible with *every* mechanism in a
# two-parameter family that includes both of the book's derivations as end
# members. That is a much narrower claim than the name of the case suggests, and
# saying so is the point of the page.
#
# **Two printed defects, reported and not repaired**, both provable from the
# book's own equations. On p. 24: *"the rate levels off and becomes zero order
# with respect to the reactant, $r_A = k_1 C_E^0$"* [sic] - but eq. (1.5.1-8)'s
# own limit is $k_2 C_{\rm E}^0$, eq. (1.5.1-15) calls $k_2 C_{\rm E}^0$ "the
# maximum possible rate", the book's Problem 2.4(b) on p. 144 asks the student to
# *show* that limit, and $k_1 C_{\rm E}^0$ has units of s$^{-1}$, which a rate
# cannot have. On p. 25: *"the Lineweaver-Burke plot of $1/r$ versus $C_A$"*
# [sic] - $1/r$ is exactly linear in $1/C_{\rm A}$ (below 1e-13 off a straight line)
# and strictly convex in $C_{\rm A}$, $\mathrm{d}^2(1/r)/\mathrm{d}C_{\rm A}^2 =
# 2K_M/(r_m C_{\rm A}^3) > 0$, missing the best *least-squares* line by 50.3 % of
# its own range over 0.05-1.0 cmol/L. **That percentage is a convention and is
# reported as one**: it is the deviation from the line that minimises a
# *log-weighted* square error on the interval, which is how the page samples and
# plots it; weight the interval uniformly instead and the same curve misses by
# 66.9 %. Neither number answers *"how close can a straight line get"*, because
# the statistic is a **maximum** and the line that minimises a maximum is the
# **Chebyshev** line, not any least-squares one. That floor is **31.7 %** - a
# straight line does come that close - and it needs no weighting convention at
# all: for a convex $A/C + B$ it is closed form, with three-point equioscillation
# as its optimality certificate. What is convention-free is the sign of the second
# derivative, and that is what settles the defect.
#
# **Neither origin was consulted and neither is on disk.** Michaelis & Menten,
# *Biochem. Z.* **49**, 333-369 (1913) and Briggs & Haldane, *Biochem. J.* **19**,
# 338-339 (1925) are the origins recorded in the gallery catalogue; this page was
# built entirely from the 2011 Froment, De Wilde & Bischoff monograph, and says
# about the two 1913/1925 papers **only what that book prints about them** - which
# is two clauses and, as the page shows below, no citation at all.

# %% [markdown]
# ## Background
#
# ### The source, precisely
#
# **Gilbert F. Froment, Juray De Wilde and Kenneth B. Bischoff, *Chemical Reactor
# Analysis and Design*, 3rd edition, John Wiley & Sons (2011), ISBN
# 978-0-470-56541-4, Chapter 1 "Elements of Reaction Kinetics", section 1.5
# "Bio-kinetics", subsection 1.5.1 "Enzymatic Kinetics", book pp. 23-26
# (PDF pp. 65-68), equations (1.5.1-1) to (1.5.1-15).**
#
# Identified from the book's own title page ("Chemical Reactor Analysis and
# Design / 3rd Edition / Gilbert F. Froment, Texas A&M University / Kenneth B.
# Bischoff (deceased), University of Delaware / Juray De Wilde, Universite
# Catholique de Louvain, Belgium / John Wiley & Sons"), its imprint page and its
# Library of Congress CIP record. **The file's name says "Froment_Bischoff"; the
# document is the THIRD edition and De Wilde is a full author.** The file is
# born-digital with embedded TrueType subsets, so the text layer is clean -
# including equation numbers, Greek and sub/superscripts - and the only raster
# image in it is the 200 ppi cover JPEG.
#
# **The text layer was still not trusted for the two lines this page's defect
# findings rest on.** Both were rendered at 300 ppi and cropped to the line: p.
# 24's "$r_A = k_1C_E^0$" (the subscript is a **1**, upright, with no minus sign,
# and the line reads "...becomes zero order with respect to the reactant,
# $r_A = k_1C_E^0$. At low $C_A$ (1.5.1-8) degenerates into a first order rate
# equation") and p. 25's "Lineweaver-Burke plot of $1/r$ versus $C_A$" (the
# argument is $C_A$, not $1/C_A$). The extractor also silently mangled
# eq. (1.5.1-17) on p. 26 - it drops the leading "$= C_{A-E} +$" so that the
# equation reads as an identity for $C_{A-E}$ - which is a text-layer artefact,
# not a defect in the book; the render shows the equation is correct. That
# equation is outside this page's scope in any case (see *Scope*).
#
# ### This is a canonical-source page, and here is exactly what that means
#
# The catalogue reference for `J4.6` is "(1913); (1925)": Michaelis, L. & Menten,
# M. L., *Die Kinetik der Invertinwirkung*, **Biochemische Zeitschrift 49**,
# 333-369 (1913), and Briggs, G. E. & Haldane, J. B. S., **Biochemical Journal
# 19**, 338-339 (1925). **Neither is on disk and neither was consulted.** The
# page is built under the repository's textbook-canonical-source rule: the model
# is read from a monograph that states, attributes and carries the result, the
# monograph is named everywhere, and the origins go in `reference` while the text
# actually read goes in `reference_read_from`.
#
# The discipline that goes with that rule is *attribution*. **Nothing on this page
# of the form "Michaelis and Menten showed / assumed / argued", or "Briggs and
# Haldane objected", appears anywhere**, because no sentence of that form appears
# in the book. What the book prints about the origins is two clauses:
#
# > "This is the Michaelis-Menten equation for the rate of a simple enzymatic
# > reaction and $K_M = k_{-1}/k_1$ is known as the Michaelis-Menten constant."
# > (p. 24)
#
# > "In the second approach, formulated by Briggs and Haldane, the formation of
# > the complex A-E does not necessarily reach equilibrium, but its concentration
# > is eliminated by applying the pseudo steady state approximation" (p. 24)
#
# and, in Problem 2.4 on p. 144, the compound name *"The Michaelis-Menten
# (Briggs-Haldane) mechanism"* - under which the book prints
# $K_m = (k_2+k_3)/k_1$, i.e. the **Briggs-Haldane** constant. The book treats
# the two as one mechanism, which is exactly why this case's discriminating
# question is worth asking.
#
# **The book prints no citation for either origin.** That is a negative claim, so
# here is the search. (1) A case-insensitive full-text search of the book's own
# text layer, all 902 PDF pages, for `michaelis|menten|briggs|haldane` returns
# **fourteen lines**, the same fourteen in all three `pdftotext` modes (default,
# `-layout`, `-raw`): **five** in sections 1.5.1-1.5.2 (PDF 66 x3, 67, 69),
# **three** in Problem 2.4 (PDF 186), **three** in later running text and problem
# sets (PDF 436, 494, 503 - the last citing "Bischoff [1966]"), and **three** in
# the Subject Index (PDF 899: the Michaelis-Menten constant, equation and
# kinetics). None is a bibliographic entry. The index's fourth enzyme-related
# entry, "Lineweaver-Burke plot, 25", does not match that search pattern and is
# not counted here. (2) The
# **Chapter 1 REFERENCES list** (book pp. 58-59) was read in full at the
# alphabetical positions where the entries would fall - between "Boudart" and
# "Caddell" for Briggs, between "Mc Laughlin" and "McQuarrie"/"Monod" for
# Michaelis. Neither is there, **although Graef & Andrews [1973], Monod [1949] and
# Williams [1967] - the other named results of section 1.5 - all are.** (3) The
# Subject Index (book p. 857) carries page numbers only. So the two results this
# case is named for are the only named results in section 1.5 that the book
# attributes without citing. Everything this page says about 1913 and 1925 is
# bounded by the two clauses above.
#
# ### The second book, and what it is used for
#
# **Levenspiel, *Chemical Reaction Engineering*, 3rd edn (Wiley, 1999), Ch. 27, is
# NOT a second source for the derivation and is not adjudicated against Froment.**
# It is cited here for two specific things, both of which this page uses:
#
# 1. p. 615 prints the reciprocal plot's axes as *"$1/(-r_{\rm A})$ versus
#    $1/C_{\rm A}$ ... the Lineweaver plot"*, read on a 600 ppi native render.
#    That is independent corroboration of the second printed defect - though the
#    defect is settled by Froment's own eq. (1.5.1-15), not by Levenspiel.
# 2. the same page states, of the reciprocal plots, that fitting $C_{\rm A}$
#    versus $\tau$ directly *"is direct, is less prone to fiddling, and is more
#    reliable"*. That is a testable claim and the page measures it on the one
#    dataset it has.
#
# ### Where this sits in the gallery
#
# `J4.7` (immobilised enzyme particle) is this rate law inside a diffusion
# problem, which is `B1.1` with $r_m C/(K_M+C)$ substituted for the power law -
# the same `S3` machinery. `J4.8` (ASM1) and Monod growth are the microbial
# descendants; Froment's own section 1.5.2 says so in as many words, calling Monod
# kinetics *"shaped after the Michaelis-Menten kinetics for enzymatic
# reactions"*. `C1.1` (Langmuir-Hinshelwood-Hougen-Watson) is the heterogeneous
# catalysis twin, and Froment says that too: eq. (1.5.1-8), he writes, *"is
# entirely similar to the Hougen-Watson rate equations that will be derived in
# Chapter 2"*. The identifiability result on this page is the same shape as
# `B2.3`'s degeneracy families and `C1.2`'s discrimination floor: a rate form that
# several mechanisms produce, and a measurement that cannot tell them apart.
#
# ### Scope
#
# **In:** eqs. (1.5.1-1) to (1.5.1-15) - the scheme, both derivations, both
# constants, the saturating form and its half-max property. The batch
# integration of eq. (1.5.1-15), which the book does not print but which is what
# its own p. 25 sentence about "the rigorous solution, which requires numerical
# methods" is about. Rawlings & Ekerdt's Exercise 9.15 table.
#
# **Out:** competitive inhibition, eqs. (1.5.1-16) to (1.5.1-18) - a different
# case; and all of section 1.5.2 (Monod, Graef-Andrews, Williams' structured
# model), which belongs to `J4.9`/`J4.8`. No transport, no particle: this page is
# a well-mixed batch, structure `S1`.

# %% [markdown]
# ## The published model
#
# Every equation number below is the book's. Symbols are the book's:
# $C_{\rm A}$ is the reactant (substrate) concentration, $C_{\rm E}$ the free
# enzyme, $C_{A-E}$ the complex, $C_{\rm E}^0$ the total enzyme,
# $C_{\rm E}^0 = C_{\rm E} + C_{A-E}$. Rates are positive by the book's own
# convention ("a rate is positive; it is in the mass balance that the formation or
# disappearance leads to a quantity that is negative or positive", p. 23).
#
# ### The scheme and the three mass-action rates
#
# $$\mathrm{A} + \mathrm{E} \; \underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}} \;
# \mathrm{A}\text{-}\mathrm{E} \; \overset{k_2}{\longrightarrow} \; \mathrm{E} + \mathrm{P}
# \tag{1.5.1-1}$$
#
# $$r_{\rm A} = k_1 C_{\rm A} C_{\rm E} - k_{-1} C_{A-E} \tag{1.5.1-2}$$
# $$r_{\rm P} = k_2 C_{A-E} \tag{1.5.1-3}$$
# $$r_{A-E} = k_1 C_{\rm A} C_{\rm E} - (k_{-1} + k_2) C_{A-E} \tag{1.5.1-4}$$
#
# ### Route 1 - decomposition rate determining, so binding equilibrates
#
# $$C_{A-E} = K\, C_{\rm A} C_{\rm E}, \quad K = k_1/k_{-1} \tag{1.5.1-5}$$
# $$r_{\rm P} = k_2 C_{A-E} = k_2 K C_{\rm A} C_{\rm E} \tag{1.5.1-6}$$
# $$r_{\rm P} = \frac{k_2 K C_{\rm A} C_{\rm E}^0}{1 + K C_{\rm A}} \tag{1.5.1-7}$$
# $$r_{\rm A} = r_{\rm P} = \frac{k_2 C_{\rm A} C_{\rm E}^0}{K_M + C_{\rm A}},
# \qquad K_M = \frac{k_{-1}}{k_1} \tag{1.5.1-8}$$
#
# ### Route 2 - Briggs and Haldane: pseudo steady state on the complex
#
# $$\frac{\mathrm{d}C_{A-E}}{\mathrm{d}t} = 0 \tag{1.5.1-9}$$
# $$C_{A-E} = \frac{k_1 C_{\rm A} C_{\rm E}}{k_{-1}+k_2} \tag{1.5.1-10}$$
# $$r_{\rm A} = k_1 C_{\rm A} C_{\rm E}\left(1 - \frac{k_{-1}}{k_{-1}+k_2}\right) \tag{1.5.1-11}$$
# $$C_{\rm E} = \frac{C_{\rm E}^0}{1 + \dfrac{k_1 C_{\rm A}}{k_{-1}+k_2}} \tag{1.5.1-12}$$
# $$r_{\rm A} = r_{\rm P} = \frac{k_2 C_{\rm A} C_{\rm E}^0}{K_M + C_{\rm A}},
# \qquad K_M = \frac{k_{-1}+k_2}{k_1} \tag{1.5.1-13},(1.5.1-14)$$
#
# ### The common form
#
# $$r_{\rm A} = \frac{r_m C_{\rm A}}{K_M + C_{\rm A}}, \qquad r_m = k_2 C_{\rm E}^0
# \tag{1.5.1-15}$$
#
# *"where $r_m = k_2 C_E^0$ is the maximum possible rate. The rate reaches half
# this maximum value when $C_A$ equals $K_M$."*
#
# The first cell below closes this chain symbolically: ten identities, each one
# re-deriving a numbered equation from the numbered equations the book derives it
# from. These are **identities, not agreement** - what they protect is the
# transcription and the reading of the two routes, and they are registered as
# structural claims.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.28.3"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "J4.6-michaelis-menten-briggs-haldane"  # equals the page directory

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
from functools import lru_cache

import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.stats import f as f_dist

from pymrm.checks import find_roots, observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import (s_closed, ReducedCurve, elem, mass_action, max_gap,
                   BatchMassAction, BatchReduced)

plt.rcParams.update({"figure.dpi": 110, "font.size": 9.5})
# Okabe-Ito, assigned in fixed order and never cycled
C_BLUE, C_ORANGE, C_GREEN = "#0072B2", "#D55E00", "#009E73"
C_PURPLE, C_YELLOW, C_GREY = "#CC79A7", "#E69F00", "0.45"
np.set_printoptions(precision=6, suppress=False)

# DETERMINISM: two consecutive executions must give identical content.
# (i) ipykernel flushes stdout on a timer, so a cell that prints and then
#     computes for seconds gets split into a different number of stream outputs
#     each run. Batch the whole cell into one message instead.
# (ii) a pandas Styler's text/plain repr is a MEMORY ADDRESS and its HTML
#      carries a random CSS id prefix, so Stylers are rendered through
#      HTML(...to_html()) with a pinned uuid, never displayed directly.
try:
    sys.stdout.flush_interval = 1e6
except Exception:
    pass
from IPython.display import HTML, display


def show(styler):
    display(HTML(styler.to_html()))


claims = Claims(PAGE.split("-")[0])
M = {}          # every computed number, registered as a claim in the Claims cell
RES = 5e-4      # half-ulp of the source table's three printed decimals, cmol/L


def root1(f, lo, hi, n_scan, xtol, what):
    """The single root of f on [lo, hi] (pymrm.checks.find_roots); raises otherwise."""
    r = find_roots(f, lo, hi, n_scan=n_scan, xtol=xtol)
    if len(r["roots"]) != 1 or r["failed"]:
        raise RuntimeError(f"{what}: {len(r['roots'])} roots on [{lo:g}, {hi:g}], "
                           f"{len(r['failed'])} failed points: {r['message']}")
    return float(r["roots"][0])


print("setup done")

# %%
# ---- ten identities closing Froment's printed chain, eqs. (1.5.1-2)..(1.5.1-15)
k1, km1, k2 = sp.symbols("k1 k_m1 k2", positive=True)
CA, CE, CAE, CE0s = sp.symbols("C_A C_E C_AE C_E0", positive=True)
rm_s, KM_s, uu = sp.symbols("r_m K_M u", positive=True)

rA_2   = k1*CA*CE - km1*CAE                                  # (1.5.1-2)
rP_3   = k2*CAE                                              # (1.5.1-3)
rAE_4  = k1*CA*CE - (km1 + k2)*CAE                           # (1.5.1-4)
CAE_5  = (k1/km1)*CA*CE                                      # (1.5.1-5), K = k1/k-1
rP_7   = k2*(k1/km1)*CA*CE0s/(1 + (k1/km1)*CA)               # (1.5.1-7)
rA_8   = k2*CA*CE0s/(km1/k1 + CA)                            # (1.5.1-8), K_M = k-1/k1
CAE_10 = k1*CA*CE/(km1 + k2)                                 # (1.5.1-10)
rA_11  = k1*CA*CE*(1 - km1/(km1 + k2))                       # (1.5.1-11)
CE_12  = CE0s/(1 + k1*CA/(km1 + k2))                         # (1.5.1-12)
rA_13  = k2*CA*CE0s/((km1 + k2)/k1 + CA)                     # (1.5.1-13)+(1.5.1-14)
r_15   = rm_s*CA/(KM_s + CA)                                 # (1.5.1-15)

SYM = {
    "(1.5.1-4) = (1.5.1-2) - (1.5.1-3)":
        rAE_4 - (rA_2 - rP_3),
    "(1.5.1-7) from (5), (6) and C_E0 = C_E + C_A-E":
        k2*sp.solve(sp.Eq(CAE, CAE_5.subs(CE, CE0s - CAE)), CAE)[0] - rP_7,
    "(1.5.1-8) = (1.5.1-7) with K_M = k_-1/k_1":
        rP_7 - rA_8,
    "(1.5.1-10) from (1.5.1-4) = 0, i.e. (1.5.1-9)":
        sp.solve(sp.Eq(rAE_4, 0), CAE)[0] - CAE_10,
    "(1.5.1-11) from (1.5.1-2) and (1.5.1-10)":
        rA_2.subs(CAE, CAE_10) - rA_11,
    "(1.5.1-12) from C_E0 = C_E + C_A-E and (1.5.1-10)":
        sp.solve(sp.Eq(CE0s, CE + CAE_10), CE)[0] - CE_12,
    "(1.5.1-13)+(1.5.1-14) from (1.5.1-11) and (1.5.1-12)":
        rA_11.subs(CE, CE_12) - rA_13,
    "(1.5.1-8) is the k_2 -> 0 limit of (1.5.1-13)":
        sp.limit(rA_13, k2, 0) - sp.limit(rA_8, k2, 0),
    "(1.5.1-15) is exactly r_m/2 at C_A = K_M":
        r_15.subs(CA, KM_s) - rm_s/2,
    "1/r is exactly linear in u = 1/C_A":
        sp.diff((1/r_15).subs(CA, 1/uu), uu, 2),
}
SYM = {k: sp.simplify(v) for k, v in SYM.items()}
for k, v in SYM.items():
    print(f"  {v == 0}   {k}")
# the largest residual, evaluated at random positive values of every symbol: zero
# for an identity, nonzero for any residual sympy could not simplify away
_sym_rng = np.random.default_rng(20260813)
_vals = {s: float(_sym_rng.uniform(0.5, 2.0)) for s in (k1, km1, k2, CA, CE, CAE, CE0s,
                                                        rm_s, KM_s, uu)}
M["froment_chain_symbolic_max_residual"] = max(abs(float(v.subs(_vals))) for v in SYM.values())
M["froment_chain_identities_verified"] = sum(v == 0 for v in SYM.values())

CURV  = sp.simplify(sp.diff(1/r_15, CA, 2))
LIM8  = sp.limit(rA_8, CA, sp.oo)
LIM13 = sp.limit(rA_13, CA, sp.oo)
KMRAT = sp.simplify(((km1 + k2)/k1) / (km1/k1))
print()
print("  d2(1/r)/dC_A^2                      =", CURV, " (> 0 for K_M > 0)")
print("  lim_{C_A -> oo} of (1.5.1-8)        =", LIM8)
print("  lim_{C_A -> oo} of (1.5.1-13)       =", LIM13)
print("  K_M(Briggs-Haldane)/K_M(Michaelis-Menten) =", KMRAT, "= 1 + k2/k_-1")

# ---- THE DEGENERACY, PROVED rather than sampled. This is the page's central
# structural result; a grid can only *check* it.
rho_s = sp.symbols("rho", positive=True)
REPARAM = {k2: rm_s/CE0s, km1: (rm_s/CE0s)/rho_s,
           k1: ((rm_s/CE0s)/rho_s + rm_s/CE0s)/KM_s}
rA_13_rep = sp.simplify(rA_13.subs(REPARAM, simultaneous=True))
DEG_SAME  = sp.simplify(rA_13_rep - r_15)
DEG_DRHO  = sp.simplify(sp.diff(rA_13_rep, rho_s))
M["rho_degeneracy_symbolic_residual"] = max(
    abs(float(DEG_SAME.subs(_vals | {rho_s: 1.3}))), abs(float(DEG_DRHO.subs(_vals | {rho_s: 1.3}))))
print()
print("  SYMBOLIC PROOF of the rho-degeneracy (not a sweep):")
print("    (1.5.1-13) under k2 = r_m/C_E0, k_-1 = k2/rho, k1 = (k_-1+k2)/K_M")
print("      =", rA_13_rep, " = eq. (1.5.1-15), independent of rho")
print("    d/drho of that expression =", DEG_DRHO, " identically")

# %% [markdown]
# **The two printed defects, from the book's own equations.**
#
# *Defect 1, p. 24, quoted verbatim:* "At high reactant concentration ($C_A$ much
# larger than $K_M$), the rate levels off and becomes zero order with respect to
# the reactant, $r_A = k_1C_E^0$." [sic]
#
# The cell above prints the limit of the book's own eq. (1.5.1-8) as
# $C_{\rm A}\to\infty$: it is $k_2 C_{\rm E}^0$, and so is the limit of
# eq. (1.5.1-13). Two more of the book's own statements agree: eq. (1.5.1-15)
# introduces $r_m = k_2 C_{\rm E}^0$ as *"the maximum possible rate"*, and
# Problem 2.4(b) on p. 144 instructs the student to *"Show that the maximum initial
# rate is given by $-\mathrm{d}[S]/\mathrm{d}t|_{max} = k_3[E_0]$"*, where that
# problem's $k_3$ is section 1.5.1's $k_2$. And the printed expression is
# **dimensionally impossible**: $k_1$ is a second-order coefficient
# (m$^3$ mol$^{-1}$ s$^{-1}$), so $k_1 C_{\rm E}^0$ has units of s$^{-1}$, where a
# rate needs mol m$^{-3}$ s$^{-1}$. Reported, not repaired; the transcription CSV
# carries the printed form.
#
# *Defect 2, p. 25, quoted verbatim:* "One example is the Lineweaver-Burke plot of
# $1/r$ versus $C_A$." [sic]
#
# $1/r = (K_M/r_m)(1/C_{\rm A}) + 1/r_m$ is exactly affine in $1/C_{\rm A}$ - the
# cell above verifies $\mathrm{d}^2(1/r)/\mathrm{d}u^2 = 0$ identically in
# $u = 1/C_{\rm A}$ - while in $C_{\rm A}$ it is strictly convex,
# $\mathrm{d}^2(1/r)/\mathrm{d}C_{\rm A}^2 = 2K_M/(r_m C_{\rm A}^3) > 0$. The
# next-but-one cell measures how far from a straight line that is on the fitted
# constants. Levenspiel p. 615 prints the same plot with $1/C_{\rm A}$ on the
# abscissa. **He prints no surname at all**: his line reads
# *"$1/(-r_{\rm A})$ versus $1/C_{\rm A}$ ... the Lineweaver plot"*, read on a
# 600 ppi native render, and the string "Burk" does not occur anywhere in that
# book (case-insensitive full-text search of all 684 PDF pages, in all three
# `pdftotext` modes: two hits for "Lineweaver", p. 615 and the index; none for
# "Burk"). Froment prints "Lineweaver-Burke". Reported, not
# repaired - and nothing on this page turns on the spelling.
#
# **A defect this page does NOT report.** Extracting eq. (1.5.1-17) on p. 26 from
# the text layer produces "$C_{A-E} = (k_3 C_I/k_{-3} + 1)C_E$", which is wrong.
# The 300 ppi render shows the printed equation is
# "$= C_{A-E} + (k_3 C_I/k_{-3} + 1)C_E$" and is correct: the extractor dropped
# the leading term. Recorded here because it is exactly the failure mode the
# repository warns about - a clean, born-digital text layer is still not a page.

# %% [markdown]
# ## Parameters and assumptions
#
# There are no physical constants to transcribe: section 1.5.1 prints **no
# numbers at all**, no worked example and no table. Everything numeric on this page
# is one of three things, and each is labelled wherever it appears.
#
# 1. **Fitted** to the seven-point table of Rawlings & Ekerdt Exercise 9.15:
#    $r_m$ and $K_M$ (equivalently their $k$ and $K$). The agreement between the
#    Michaelis-Menten curve and those seven points is therefore a **goodness of
#    fit**, not a test, and a null baseline is computed beside it.
# 2. **Structural** - exact algebraic properties of the printed equations
#    (the ten identities, the degeneracy, the convexity of $1/r$ in $C_{\rm A}$).
#    These cannot fail for the right reason and are named as such.
# 3. **Computed properties** of the printed model: the cost of each reduction
#    against the full mass-action solution, and the enzyme loading an experiment
#    would need before $\rho = k_2/k_{-1}$ leaves a trace.
#
# ### The two dimensionless groups
#
# Both are defined here, not taken from anywhere:
#
# $$\rho \equiv \frac{k_2}{k_{-1}}, \qquad
# \varepsilon \equiv \frac{C_{\rm E}^0}{c_{\rm A0} + K_M}.$$
#
# $\rho$ is the whole difference between the two derivations:
# $K_M^{\rm BH} = (1+\rho)\,K_M^{\rm MM}$. $\varepsilon$ is the group that the
# quasi-steady-state error actually scales with - the page establishes that by
# sweeping $K_M$ over three decades at **fixed** $C_{\rm E}^0/c_{\rm A0}$ and
# watching which of the two candidate groups collapses the error. Froment's own
# printed criterion compares $C_{\rm E}^0$ with the substrate alone ("in
# particular if the enzyme concentration is relatively large compared to that of
# the substrate"), which is the $K_M \ll c_{\rm A0}$ corner of $\varepsilon$ and
# is right there; the $K_M$ dependence is the part this page adds.
#
# ### The reparameterisation used throughout
#
# Given $(r_m, K_M, C_{\rm E}^0, \rho)$ the elementary constants follow uniquely:
#
# $$k_2 = \frac{r_m}{C_{\rm E}^0}, \qquad k_{-1} = \frac{k_2}{\rho}, \qquad
# k_1 = \frac{k_{-1}+k_2}{K_M}.$$
#
# This is a bijection onto the four-parameter mass-action model, so sweeping
# $\rho$ at fixed $(r_m, K_M, C_{\rm E}^0)$ walks along the fibre of the
# observable map without changing anything a steady-state rate measurement can
# see. Route 1's constant is then $K_M^{\rm MM} = k_{-1}/k_1 = K_M/(1+\rho)$.
#
# ### Assumptions inherited from the book, stated
#
# - Isothermal, well mixed, constant volume, batch: structure `S1`.
# - Product formation is irreversible ($k_2$ has no reverse). The book's
#   eq. (1.5.1-1) prints a single arrow, and so does the page. Note that Problem
#   2.4 in Chapter 2 does *not*: it derives a reversible form with a
#   $[P]/K$ term. That is a different (larger) model and is out of scope.
# - Free enzyme and complex are the only enzyme species: $C_{\rm E}^0 =
#   C_{\rm E} + C_{A-E}$, the book's own balance.
# - No inhibitor ($C_I = 0$), which is what takes eqs. (1.5.1-16)-(1.5.1-18) out
#   of scope.

# %% [markdown]
# ## The data
#
# Section 1.5.1 publishes none, so the empirical half of this page comes from a
# **different book**: Rawlings & Ekerdt, *Chemical Reactor Analysis and Design
# Fundamentals*, 2nd edn (Nob Hill, 2025 printing), **Exercise 9.15 "Parameters of
# Michaelis-Menten kinetics"**, which prints seven $c_{\rm S}$-versus-time pairs.
# That book is a **data source here and nothing else** - the derivation, the
# attribution and both constants come from Froment, and no derivation is
# adjudicated between the two.
#
# **Provenance: unstated, and the page says so everywhere.** The exercise
# introduces the table with one sentence - *"The following measurements of cS
# versus time were taken in your laboratory."* - and gives no citation, no enzyme,
# no enzyme loading, no temperature, no pH and no solvent. The whole exercise was
# read and the book's bibliography searched; there is no reference attached to
# these numbers. "Taken in your laboratory" is a textbook framing addressed to a
# student, not a provenance claim, so **this page does not call the table a
# laboratory record.** Two things it *can* establish about it are computed below:
# the residual is far above the print resolution, so the table is not an
# exactly-rounded evaluation of the model it is set to fit; and the fitted
# constants land within a fraction of a percent of round values, which is
# consistent with - not proof of - a synthetic table.
#
# **$C_{\rm E}^0$ is not printed**, and for this case that is not a detail. Since
# $r_m = k_2 C_{\rm E}^0$, no fit can separate $k_2$ from the enzyme loading; and
# $\varepsilon$, the group that decides whether the quasi-steady-state
# approximation holds for *this* run, cannot be evaluated at all. Every
# $\varepsilon$ on this page is therefore a property of a *designed* experiment,
# never of Exercise 9.15's.
#
# The exercise's own hint - *"consider the transformation $1/r$ and formulate a
# linear least-squares problem that you can solve on a calculator"* - is the
# Lineweaver construction, and the page runs it as instructed alongside the direct
# fit, because that comparison is the one Levenspiel p. 615 makes a claim about.

# %%
DATA = load_data("rawlings-ekerdt-ex9.15-batch.csv", page=PAGE)
DMETA = load_meta("rawlings-ekerdt-ex9.15-batch.csv", page=PAGE)
PRINTED = load_data("froment-1.5.1-printed.csv", page=PAGE)
print(cite_data(DMETA))
print(cite_data(load_meta("froment-1.5.1-printed.csv", page=PAGE)))
t_d = DATA["t_min"].to_numpy(float)
c_d = DATA["c_S_cmol_per_L"].to_numpy(float)
cA0 = c_d[0]


def ma(p, CE0, te, a0=None, **kw):
    """model.mass_action started from the printed c_A0 unless a0 is given."""
    return mass_action(p, CE0, te, cA0 if a0 is None else a0, **kw)


display(DATA)
print(f"{len(t_d)} points, t = {t_d[0]:g}..{t_d[-1]:g} min, "
      f"c_S = {c_d.min():g}..{c_d.max():g} cmol/L")
print(f"conversion reached: {100*(1 - c_d[-1]/c_d[0]):.2f} %")
print()
print("provenance, from the dataset sidecar:")
_pv = " ".join(DMETA["notes"]["PROVENANCE_IS_NOT_STATED"].split())
print("  " + _pv[:_pv.index("So the table CANNOT")].strip())
print()
print(f"printed equations/quotes transcribed from Froment 1.5.1: {len(PRINTED)}")
print("  the two rows flagged as printed defects:")
for _, r in PRINTED[PRINTED["note"].str.startswith("PRINTED DEFECT")].iterrows():
    print(f"    p.{r['page']} [{r['key']}]  {r['as_printed']}")


# %% [markdown]
# ## PyMRM implementation
#
# Structure `S1`: an initial-value problem in time, in a single well-mixed cell.
# The pymrm way to march one is **backward Euler with `newton`, and `NumJac` for
# the pointwise reaction Jacobian** - which is what the two classes below do.
# Three models are needed, and keeping them separate is the point of the page:
#
# | model | state | what it is |
# |---|---|---|
# | `BatchMassAction` | $(C_{\rm A}, C_{\rm E}, C_{A-E}, C_{\rm P})$, shape `(1, 4)` | the scheme (1.5.1-1) with **no** approximation - the reference the two reductions are priced against |
# | `BatchReduced` | $C_{\rm A}$, shape `(1, 1)` | eq. (1.5.1-15) marched in time |
# | `s_closed` | - | the closed-form integral of eq. (1.5.1-15), $K_M\ln(c_{A0}/C_{\rm A}) + (c_{A0}-C_{\rm A}) = r_m t$, inverted with Lambert $W$ |
#
# Two `NumJac` details that matter and are easy to get wrong:
#
# - **`NumJac((1, 4))`, never `NumJac((4,))`.** The default stencil couples the
#   *last* axis in full, which is exactly right when the last axis is the field
#   index and the source term is pointwise. Passing a bare 1-D shape would make
#   space the last axis and build a dense Jacobian.
# - **No `axes_diagonals`.** The source term here reads only its own cell; there
#   is no spatial coupling at all in a batch.
#
# The reference solution for the mass-action system is `scipy`'s `Radau` with an
# analytic Jacobian; `LSODA` with the same Jacobian is used for the sweeps because
# it is ~20x faster on this stiff system, and the two are compared. The Lambert
# $W$ closed form and the pymrm marcher are two independent computations of the
# same reduced model, and the fit is done through **both**.

# %%
# The three models live in model.py: s_closed (Lambert W), ReducedCurve (the
# same curve with its time derivative), elem (rho, C_E0, r_m, K_M) -> (k_1, k_-1,
# k_2), mass_action (scipy reference, dense output with the analytic derivative),
# max_gap (largest difference in time, at a root of the analytic derivative), and
# the two pymrm marchers BatchMassAction and BatchReduced.
print("models:", BatchMassAction.__name__, BatchReduced.__name__, ReducedCurve.__name__,
      s_closed.__name__, mass_action.__name__, max_gap.__name__)


# %% [markdown]
# ## Results
#
# ### 1. The fit, and what it is worth against a null
#
# Fitting eq. (1.5.1-15) to the seven points means fitting its **integral**, not
# the rate law: the table is a batch decay curve, not a set of initial rates. That
# integral is closed-form, so the fit needs no ODE solver at all. Rawlings &
# Ekerdt's parameterisation is $r = k c_{\rm S}/(1+K c_{\rm S})$, which is
# eq. (1.5.1-15) with $k = r_m/K_M$ and $K = 1/K_M$; both are reported.
#
# **This is a fit, not a test.** The null baselines below are the point of the
# exercise: a saturating law has two parameters and a simple order has one, so
# "it fits" means nothing until the one-parameter alternatives have been tried.

# %%
def fit_mm(t, y, K_fixed=None, s0_free=False):
    """Multistart NLS of the integrated (1.5.1-15) to (t, y). Returns k, K, s0, resid."""
    def resid(q):
        i = 0
        kk = np.exp(q[i]); i += 1
        KK = K_fixed if K_fixed is not None else np.exp(q[i]); i += K_fixed is None
        ss = np.exp(q[i]) if s0_free else y[0]
        return s_closed(t, ss, kk/KK, 1/KK) - y

    best = None
    for k0 in (0.05, 0.2, 0.5, 1.0, 3.0):
        for K0 in (0.2, 1.0, 3.0, 10.0):
            g = [np.log(k0)] + ([] if K_fixed is not None else [np.log(K0)]) \
                + ([np.log(y[0])] if s0_free else [])
            r = least_squares(resid, g, xtol=1e-15, ftol=1e-15, gtol=1e-15)
            if best is None or r.cost < best.cost:
                best = r
    q = np.exp(best.x)
    kk = q[0]
    KK = K_fixed if K_fixed is not None else q[1]
    ss = q[-1] if s0_free else y[0]
    return kk, KK, ss, best.fun


k_fit, K_fit, _, res_fit = fit_mm(t_d, c_d)
rm_fit, KM_fit = k_fit/K_fit, 1.0/K_fit
rms_fit = float(np.sqrt(np.mean(res_fit**2)))
sse_fit = float(np.sum(res_fit**2))
M.update(mm_fit_k_per_min=k_fit, mm_fit_K_L_per_cmol=K_fit,
         mm_fit_rm_cmol_per_L_min=rm_fit, mm_fit_KM_cmol_per_L=KM_fit,
         mm_fit_rms_cmol_per_L=rms_fit,
         mm_fit_maxres_over_print_resolution=float(np.max(np.abs(res_fit))/RES))

print(f"Exercise 9.15 parameterisation : k   = {k_fit:.6f} 1/min, "
      f"K   = {K_fit:.6f} L/cmol")
print(f"Froment (1.5.1-15) equivalent  : r_m = {rm_fit:.6f} cmol/(L min), "
      f"K_M = {KM_fit:.6f} cmol/L")
print(f"RMS residual {rms_fit:.6e} cmol/L; max |residual| "
      f"{np.max(np.abs(res_fit)):.6e} = "
      f"{M['mm_fit_maxres_over_print_resolution']:.2f} x the 5e-4 print resolution")
print("residuals:", np.array2string(res_fit, precision=6))

# is the table an exactly-rounded model evaluation? and are the constants round?
_, _, _, res_K2 = fit_mm(t_d, c_d, K_fixed=2.0)
M["mm_fit_sse_penalty_K_pinned_at_2"] = float(np.sum(res_K2**2)/sse_fit - 1.0)
M["mm_fit_KM_rel_dev_from_half"] = abs(KM_fit/0.5 - 1)
print(f"\nK pinned at exactly 2 L/cmol (K_M = 0.5): SSE penalty "
      f"{100*M['mm_fit_sse_penalty_K_pinned_at_2']:.3f} % -- "
      f"K_M is within {100*abs(KM_fit/0.5 - 1):.2f} % of 0.5 cmol/L")

k3, K3, s03, res3 = fit_mm(t_d, c_d, s0_free=True)
M["mm_fit_KM_with_s0_free"] = 1.0/K3
print(f"initial concentration also free: k = {k3:.5f}, K = {K3:.5f} "
      f"(K_M = {1/K3:.6f}), c_A0 = {s03:.6f} vs printed {cA0:.3f}; "
      f"RMS {np.sqrt(np.mean(res3**2)):.6e}")

# --- null baselines: one-parameter simple orders, each best-fitted
def null_fit(kind, t=t_d, y=c_d):
    """Best one-parameter simple-order fit to (t, y): (RMS, rate constant)."""
    f = {"zero":   lambda q: np.maximum(y[0] - np.exp(q[0])*t, 0.0) - y,
         "first":  lambda q: y[0]*np.exp(-np.exp(q[0])*t) - y,
         "second": lambda q: y[0]/(1 + np.exp(q[0])*y[0]*t) - y}[kind]
    r = least_squares(f, [np.log(0.2)], xtol=1e-15, ftol=1e-15)
    return float(np.sqrt(np.mean(r.fun**2))), float(np.exp(r.x[0]))


NULLS = {k: null_fit(k) for k in ("zero", "first", "second")}
M["null_first_order_rms_ratio"] = NULLS["first"][0]/rms_fit
M["null_first_order_rms_cmol_per_L"] = NULLS["first"][0]
M["null_zero_order_rms_cmol_per_L"] = NULLS["zero"][0]
M["null_second_order_rms_cmol_per_L"] = NULLS["second"][0]
M["null_zero_order_rms_ratio"] = NULLS["zero"][0]/rms_fit
M["null_second_order_rms_ratio"] = NULLS["second"][0]/rms_fit
print("\nnull baselines (one parameter each, all best-fitted):")
for kind, (rr, kk) in NULLS.items():
    print(f"  {kind:6s} order  k = {kk:.5f}   RMS = {rr:.5e}  "
          f"({rr/rms_fit:6.2f} x the Michaelis-Menten fit)")

# %%
fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.1))
tt = np.linspace(0, 18, 400)
ax[0].plot(tt, s_closed(tt, cA0, rm_fit, KM_fit), color=C_BLUE, lw=1.6,
           label=f"eq. (1.5.1-15), fitted ($K_M$ = {KM_fit:.3f})")
ax[0].plot(tt, cA0*np.exp(-NULLS["first"][1]*tt), color=C_ORANGE, lw=1.4, ls="--",
           label="first-order null")
ax[0].plot(t_d, c_d, "o", ms=5, color="0.15", label="Exercise 9.15 (provenance unstated)")
ax[0].set_xlabel("t (min)"); ax[0].set_ylabel(r"$c_{\rm S}$ (cmol/L)")
ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=0.25, lw=0.5)

ax[1].axhline(0, color=C_GREY, lw=0.8)
ax[1].axhspan(-RES, RES, color=C_GREY, alpha=0.18,
              label="print resolution $\\pm 5\\times10^{-4}$")
ax[1].plot(t_d, res_fit, "o-", ms=5, lw=1.2, color=C_BLUE, label="Michaelis-Menten")
ax[1].plot(t_d, cA0*np.exp(-NULLS["first"][1]*t_d) - c_d, "s--", ms=4, lw=1.0,
           color=C_ORANGE, label="first-order null")
ax[1].set_xlabel("t (min)"); ax[1].set_ylabel("model - data (cmol/L)")
ax[1].legend(fontsize=7.5, frameon=False); ax[1].grid(alpha=0.25, lw=0.5)
fig.tight_layout(); plt.show()
print("Fit, not test: r_m and K_M were fitted to these same seven points. The "
      "residual band shows why the table cannot be an exactly-rounded model "
      "evaluation - it is 13.9x too wide for that.")

# %% [markdown]
# ### 2. The Lineweaver construction, run exactly as both books instruct
#
# Froment says to plot $1/r$ against $C_{\rm A}$; Levenspiel p. 615 says
# $1/(-r_{\rm A})$ against $1/C_{\rm A}$. The first cell settles which is right
# from Froment's own eq. (1.5.1-15) - by measuring, on the fitted constants, how
# far each is from a straight line. The second runs the exercise's hint on the
# actual table (rates by central difference, then unweighted linear least squares
# on the reciprocals) and compares what comes out with the direct fit, which is
# the comparison Levenspiel makes a claim about.

# %%
# "How far from a straight line" is a MAXIMUM OVER AN INTERVAL, and it is only
# well posed once WHICH straight line is named. Three answers are computed here,
# and the difference between them is larger than any grid effect:
#   (i) the CHEBYSHEV (minimax) line is the one that minimises this very
#       statistic, so it is the FLOOR - the honest answer to "how close can a
#       straight line get". It needs no weighting convention at all, and for a
#       convex A/C + B it is closed form with an equioscillation certificate.
#   (ii) a LEAST-SQUARES line minimises a different (square) norm, so it is
#       always further away in the max norm, and it additionally depends on how
#       the interval is WEIGHTED: a least-squares fit on a geomspace grid
#       weights it as dC/C, on a linspace grid as dC. Those two give 0.503 and
#       0.669; the page reports both and names the one it quotes.
#   (iii) in every case the maximum is ROOT-FOUND, not sampled: with the line
#       fixed, the deviation of A/C + B from q0 + q1 C has one interior
#       stationary point, C* = sqrt(-A/q1), in closed form.
# Moments and extrema are therefore all exact - there is no grid under any of
# these numbers. The grid sequence below is shown only as evidence that refining
# a 400-point grid converges to the closed form.
A_LW, B_LW = KM_fit/rm_fit, 1.0/rm_fit             # 1/r = A/C_A + B


def straightness_exact(a, b, weight="log", line=False):
    """max |1/r - best LEAST-SQUARES line in C_A| / range, on [a, b].

    Least squares minimises a square norm, NOT the maximum reported here, so
    this is an upper bound on the distance from the best straight line, never
    the distance itself - straightness_minimax gives that.
    `weight` = 'log' (dC/C, the geomspace convention) or 'uniform' (dC).
    `line=True` returns the line's (intercept, slope) instead, so that the
    figure can draw the same line the number is measured against.
    Closed-form weighted least squares, root-found extremum.
    """
    L = np.log(b/a)
    if weight == "log":
        m0, m1, m2 = L, b - a, (b*b - a*a)/2
        n0, n1 = A_LW*(1/a - 1/b) + B_LW*L, A_LW*L + B_LW*(b - a)
    else:
        m0, m1, m2 = b - a, (b*b - a*a)/2, (b**3 - a**3)/3
        n0, n1 = A_LW*L + B_LW*(b - a), A_LW*(b - a) + B_LW*(b*b - a*a)/2
    q0, q1 = np.linalg.solve(np.array([[m0, m1], [m1, m2]]), np.array([n0, n1]))
    if line:
        return float(q0), float(q1)
    dev = lambda x: A_LW/x + B_LW - q0 - q1*x
    xs = [a, b]
    xc = np.sqrt(-A_LW/q1) if q1 < 0 else None     # dev'(x) = -A/x^2 - q1 = 0
    if xc is not None and a < xc < b:
        xs.append(float(xc))
    return float(max(abs(dev(x)) for x in xs)/(A_LW*(1/a - 1/b)))


def dev_sup(q0, q1, a, b):
    """sup |1/r - (q0 + q1 C)| on [a, b], EXACTLY - no grid.

    1/r - line is convex in C, so its extrema on the interval are the two
    endpoints and the single interior stationary point sqrt(-A/q1).
    """
    d = lambda x: A_LW/x + B_LW - q0 - q1*x
    xs = [a, b]
    if q1 < 0 and a < np.sqrt(-A_LW/q1) < b:
        xs.append(float(np.sqrt(-A_LW/q1)))
    return max(abs(d(x)) for x in xs)


def straightness_minimax(a, b):
    """max |1/r - THE BEST STRAIGHT LINE in C_A| / range, on [a, b].

    Best in the CHEBYSHEV (minimax) sense: the line that minimises the maximum
    deviation, i.e. the statistic actually being reported. No straight line does
    better, so this is a convention-free FLOOR - it needs no weighting, and the
    normalisation by the range cancels r_m and K_M exactly, leaving a pure
    property of the interval.

    Closed form for a strictly convex f = A/C + B: the optimal slope is the
    chord slope m = (f(b)-f(a))/(b-a), the tangency is at sqrt(A/|m|) = sqrt(ab),
    and the deviation equioscillates + - + at a, sqrt(ab), b. Chebyshev's
    theorem makes that equioscillation an optimality CERTIFICATE - but by
    CONSTRUCTION here, not by an independent check: q0 is set to the midpoint
    of g(a) and g(xt), so d[0] = -d[1] identically, and m is the exact chord
    slope, so g(b) = g(a) and d[2] = d[0] identically. The equal-magnitude
    check below therefore compares each of those to itself up to round-off;
    what it actually earns is a < xt < b (tangency interior) and d[0] > 0
    (convexity, i.e. xt is a minimum of g, not a maximum). The 21x21
    perturbation scan two cells down is the part of this notebook that tests
    optimality from outside the construction.
    """
    f = lambda x: A_LW/x + B_LW
    m = (f(b) - f(a))/(b - a)                       # chord slope
    xt = float(np.sqrt(A_LW/abs(m)))                # = sqrt(ab); f'(xt) = m
    if not a < xt < b:
        raise RuntimeError(f"tangency {xt} outside [{a}, {b}]")
    g = lambda x: f(x) - m*x
    q0 = (g(a) + g(xt))/2                           # line midway between them
    d = [g(a) - q0, g(xt) - q0, g(b) - q0]
    if not (d[0] > 0 > d[1] and d[2] > 0):          # alternating signs, + - +
        raise RuntimeError(f"no equioscillation: {d}")
    if max(abs(abs(v) - abs(d[0])) for v in d) >= 1e-12*abs(d[0]):
        raise RuntimeError(f"deviations not equal in size: {d}")
    return float(abs(d[0])/(A_LW*(1/a - 1/b)))


def straightness(x, y):
    """The sampled version: max |y - lstsq line in x| / range(y) on a grid."""
    A = np.vstack([x, np.ones_like(x)]).T
    b = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(np.max(np.abs(y - A@b))/np.ptp(y))


LW_A, LW_B = 0.05, 1.0                              # spans the data's own range
c_ex = np.geomspace(LW_A, LW_B, 400)
r_ex = rm_fit*c_ex/(KM_fit + c_ex)
M["lineweaver_printed_axes_nonlinearity_minimax"] = straightness_minimax(LW_A, LW_B)
M["lineweaver_printed_axes_nonlinearity"] = straightness_exact(LW_A, LW_B, "log")
M["lineweaver_printed_axes_nonlinearity_uniform_measure"] = \
    straightness_exact(LW_A, LW_B, "uniform")
M["lineweaver_correct_axes_nonlinearity"] = straightness(1/c_ex, 1/r_ex)
print(f"exact eq. (1.5.1-15) on {LW_A}-{LW_B} cmol/L, distance from a straight "
      f"line,\nas a fraction of the range of 1/r over that interval. WHICH "
      f"straight line is\npart of the question, and it is the bigger lever of "
      f"the two conventions here:")
print(f"  1/r vs C_A, MINIMAX line -- THE best straight line   : "
      f"{M['lineweaver_printed_axes_nonlinearity_minimax']:.6f}"
      f"   <- floor, no convention")
print(f"  1/r vs C_A, best LEAST-SQUARES line, log-weighted    : "
      f"{M['lineweaver_printed_axes_nonlinearity']:.6f}   <- the page's quoted "
      f"convention")
print(f"  1/r vs C_A, best LEAST-SQUARES line, uniform         : "
      f"{M['lineweaver_printed_axes_nonlinearity_uniform_measure']:.6f}")
print(f"  1/r vs 1/C_A (Levenspiel p. 615), least squares      : "
      f"below 1e-13 (here {M['lineweaver_correct_axes_nonlinearity']:.1e})")
print(f"  and symbolically d2(1/r)/dC_A^2 = {CURV} > 0, so no straight line "
      f"exists\n  -- which is the convention-free part, and the part the defect "
      f"rests on.")
print(f"\n  THE STATISTIC IS A MAXIMUM, so the line that minimises it is the "
      f"CHEBYSHEV\n  line, not a least-squares one: a straight line does come "
      f"within "
      f"{100*M['lineweaver_printed_axes_nonlinearity_minimax']:.1f} % of the\n"
      f"  range, and the two least-squares numbers "
      f"({100*M['lineweaver_printed_axes_nonlinearity']:.1f} % and "
      f"{100*M['lineweaver_printed_axes_nonlinearity_uniform_measure']:.1f} %) "
      f"answer a\n  different question and do not even bracket it.")
_mmx = M["lineweaver_printed_axes_nonlinearity_minimax"]
_f = lambda x: A_LW/x + B_LW
_m = (_f(LW_B) - _f(LW_A))/(LW_B - LW_A)
_g = lambda x: _f(x) - _m*x
_q0 = (_g(LW_A) + _g(float(np.sqrt(A_LW/abs(_m)))))/2
_rng = A_LW*(1/LW_A - 1/LW_B)
_PN, _PD = 21, 0.05
_pert = [dev_sup(_q0*(1 + du), _m*(1 + dv), LW_A, LW_B)/_rng
         for du in np.linspace(-_PD, _PD, _PN)
         for dv in np.linspace(-_PD, _PD, _PN) if du or dv]
_worst = min(_pert)
if not _worst > _mmx:
    raise RuntimeError(f"a perturbed line beats the minimax floor: {_worst} <= {_mmx}")
print(f"  the three deviations at C = {LW_A:g}, {np.sqrt(LW_A*LW_B):.7f}, "
      f"{LW_B:g} are equal in size and\n  alternate in sign (checked), which is "
      f"Chebyshev's optimality certificate; and no\n  perturbed line in a "
      f"{_PN} x {_PN} grid of +-{100*_PD:g} % on its slope and intercept beats "
      f"it -\n  the best of those {len(_pert)} gets to {_worst:.6f} against the "
      f"floor's {_mmx:.6f}.")
print("\n  the two LEAST-SQUARES conventions are LIMITS of the two obvious "
      "grids, and a\n  400-point geomspace grid is "
      "not converged at four digits:")
for n in (400, 1600, 6400, 25600):
    cg, cl = np.geomspace(LW_A, LW_B, n), np.linspace(LW_A, LW_B, n)
    print(f"    n = {n:6d}   geomspace {straightness(cg, (KM_fit+cg)/(rm_fit*cg)):.6f}"
          f"   linspace {straightness(cl, (KM_fit+cl)/(rm_fit*cl)):.6f}")
print(f"    n -> inf     geomspace {M['lineweaver_printed_axes_nonlinearity']:.6f}"
      f"   linspace "
      f"{M['lineweaver_printed_axes_nonlinearity_uniform_measure']:.6f}"
      f"   (closed form)")


def lineweaver(tv, yv):
    """The exercise's hint: central-difference rates, then LS on 1/r vs 1/c."""
    cm = 0.5*(yv[:-1] + yv[1:])
    rr = -(yv[1:] - yv[:-1])/(tv[1:] - tv[:-1])
    A = np.vstack([1/cm, np.ones_like(cm)]).T
    sl, ic = np.linalg.lstsq(A, 1/rr, rcond=None)[0]
    return 1/sl, ic/sl, cm, rr                       # k = 1/slope, K = intercept*k


k_lb, K_lb, cm_d, r_d = lineweaver(t_d, c_d)
k_lbm, K_lbm, _, _ = lineweaver(t_d, s_closed(t_d, cA0, rm_fit, KM_fit))
M["lb_rm_bias_vs_nls"] = (k_lb/K_lb)/rm_fit - 1
M["lb_KM_bias_vs_nls"] = (1/K_lb)/KM_fit - 1
M["lb_KM_bias_finite_difference_only"] = (1/K_lbm)/KM_fit - 1
print(f"\nLineweaver on the seven printed points : k = {k_lb:.5f}, K = {K_lb:.5f}")
print(f"   -> r_m = {k_lb/K_lb:.5f} ({100*M['lb_rm_bias_vs_nls']:+.2f} %), "
      f"K_M = {1/K_lb:.5f} ({100*M['lb_KM_bias_vs_nls']:+.2f} %) vs the direct fit")
print(f"same construction on the FITTED CURVE sampled at the same seven times : "
      f"K_M = {1/K_lbm:.5f} ({100*M['lb_KM_bias_finite_difference_only']:+.2f} %)")
print("   so the differencing alone biases K_M UP by 14 %; the rest of the -47 %"
      " is the reciprocal transform amplifying the scatter of the last points.")
print("   Levenspiel p. 615 says the direct C_A-versus-tau fit 'is direct, is "
      "less prone to fiddling, and is more reliable'. On this dataset it is.")

# %%
fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.1))
# plotted over the data's own range; the METRICS above are read on 0.05-1.0,
# which is a stated choice and has its own break row
c_pl = np.geomspace(0.9*cm_d.min(), 1.0, 400)
r_pl = rm_fit*c_pl/(KM_fit + c_pl)
ax[0].plot(1/c_pl, 1/r_pl, color=C_BLUE, lw=1.6, label="eq. (1.5.1-15), exact")
A = np.vstack([1/c_pl, np.ones_like(c_pl)]).T
ax[0].plot(1/c_pl, A@np.linalg.lstsq(A, 1/r_pl, rcond=None)[0], color=C_GREY,
           lw=0.9, ls=":", label="straight line")
ax[0].plot(1/cm_d, 1/r_d, "o", ms=5, color="0.15", label="Exercise 9.15, differenced")
ax[0].set_xlabel(r"$1/C_{\rm A}$ (L/cmol)")
ax[0].set_ylabel(r"$1/r$ (L min/cmol)")
ax[0].set_title("Levenspiel p. 615 axes: exactly straight", fontsize=9)
ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=0.25, lw=0.5)

ax[1].plot(c_pl, 1/r_pl, color=C_ORANGE, lw=1.6, label="eq. (1.5.1-15), exact")
# BOTH lines are the ones fitted on 0.05-1.0 cmol/L - the interval the two
# reported percentages are measured on - drawn across the data's own range.
# The dotted one minimises a log-weighted SQUARE error (50.3 %); the dashed one
# minimises the MAXIMUM deviation, which is the statistic being reported, and is
# therefore the best straight line there is (31.7 %).
_qls = straightness_exact(LW_A, LW_B, "log", line=True)
ax[1].plot(c_pl, _qls[0] + _qls[1]*c_pl, color=C_GREY,
           lw=0.9, ls=":", label="best least-squares line")
ax[1].plot(c_pl, _q0 + _m*c_pl, color="0.35", lw=0.9, ls="--",
           label="best straight line (minimax)")
ax[1].plot(cm_d, 1/r_d, "o", ms=5, color="0.15", label="Exercise 9.15, differenced")
ax[1].set_xlabel(r"$C_{\rm A}$ (cmol/L)")
ax[1].set_ylabel(r"$1/r$ (L min/cmol)")
ax[1].set_title(f"Froment p. 25 axes, 0.05-1.0 cmol/L: no straight line "
                f"closer than "
                f"{100*M['lineweaver_printed_axes_nonlinearity_minimax']:.1f}"
                f" % of range", fontsize=9)
ax[1].legend(fontsize=7.5, frameon=False); ax[1].grid(alpha=0.25, lw=0.5)
fig.tight_layout(); plt.show()

# %% [markdown]
# ### 3. The degeneracy: what a steady-state rate measurement cannot see
#
# This is the case's own question, and it was **already settled symbolically** in
# the symbolic cell under *The published model*: under the reparameterisation, eq. (1.5.1-13) *is*
# eq. (1.5.1-15) and $\partial r/\partial \rho \equiv 0$. What follows is the
# numerical check on that algebra, not the evidence for it. Hold the two
# observable constants $(r_m, K_M)$ fixed, sweep $\rho = k_2/k_{-1}$ over eight
# decades, and evaluate the rate law **from the elementary constants** on a
# four-decade concentration grid. If the two derivations were separable by
# rate-versus-substrate data, the curves would move.
#
# The check returns exactly zero, and the cell then says what that zero is worth:
# it needs two floating-point round trips to be exact, and over 2000 random
# $(\rho, C_{\rm E}^0)$ draws they are not, so the same computation returns
# ulp-level noise instead. Both numbers are reported.

# %%
rho_dec = np.geomspace(1e-4, 1e4, 9)
cg = np.geomspace(1e-3, 10.0, 400)
base = rm_fit*cg/(KM_fit + cg)
CE0_probe = 1e-3



def degeneracy_rel(rho, CE0, k2_factor=1.0):
    """max |r(C_A) from the elementary constants / r from (1.5.1-15) - 1| on cg.

    k2_factor multiplies k_2 after the reparameterisation, with C_E0 left
    alone, so that r_m = k_2 C_E0 is no longer held fixed."""
    kf, kr, kc = elem(rho, CE0, rm_fit, KM_fit)
    kc = k2_factor*kc
    r_from_elem = kc*CE0*cg/((kr + kc)/kf + cg)              # (1.5.1-13)+(1.5.1-14)
    return float(np.max(np.abs(r_from_elem/base - 1)))


rows = []
for r_ in rho_dec:
    kf, kr, kc = elem(r_, CE0_probe, rm_fit, KM_fit)
    rows.append((r_, kf, kr, kc, kr/kf, degeneracy_rel(r_, CE0_probe)))
worst = max(r[-1] for r in rows)
M["rate_law_rho_degeneracy_max_rel"] = worst
M["KM_ratio_at_rho_1e4"] = 1.0 + 1e4

# The exact 0.0 above is a property of these nine rho at this C_E0: it needs the
# float round trips (k_-1+k_2)/k_1 == K_M and k_2 C_E0 == r_m to be exact, and
# they are here. Over random draws they are not, and the same computation returns
# ULP-level noise instead. THE ALGEBRA IS THE RESULT (proved in the first cell);
# this is the size of the arithmetic underneath it.
_rng = np.random.default_rng(20260813)             # seeded: determinism
_worst_rand, _inexact, _NDRAW = 0.0, 0, 2000
for _ in range(_NDRAW):
    r_, ce_ = 10**_rng.uniform(-4, 4), 10**_rng.uniform(-5, -1)
    kf, kr, kc = elem(r_, ce_, rm_fit, KM_fit)
    _worst_rand = max(_worst_rand, float(np.max(np.abs(
        kc*ce_*cg/((kr + kc)/kf + cg)/base - 1))))
    _inexact += not ((kr + kc)/kf == KM_fit and kc*ce_ == rm_fit)
M["rate_law_rho_degeneracy_random_draws_max_rel"] = _worst_rand

df = pd.DataFrame(rows, columns=["rho = k2/k-1", "k1 (L/cmol/min)", "k-1 (1/min)",
                                 "k2 (1/min)", "K_M(MM) = k-1/k1", "max |dr/r|"])
show(df.style.format({"rho = k2/k-1": "{:.0e}", "k1 (L/cmol/min)": "{:.4g}",
                      "k-1 (1/min)": "{:.4g}", "k2 (1/min)": "{:.4g}",
                      "K_M(MM) = k-1/k1": "{:.6f}", "max |dr/r|": "{:.2e}"})
     .set_caption("All nine mechanisms have r_m = %.6f and K_M(BH) = %.6f"
                  % (rm_fit, KM_fit)).set_uuid("j46deg"))
print(f"max relative difference in r(C_A) across eight decades of rho: below 1e-14 "
      f"(here {worst:.1e})")
print(f"K_M(Briggs-Haldane)/K_M(Michaelis-Menten) = 1 + rho over the same sweep: "
      f"{1+rho_dec[0]:.5f} .. {1+rho_dec[-1]:.4g}")
print()
print("WHAT THAT EXACT ZERO IS, AND IS NOT. The degeneracy is proved above in "
      "sympy;\nthe 0.0 here additionally requires the float round trips "
      "(k_-1+k_2)/k_1 == K_M and\nk_2 C_E0 == r_m to be exact, which they happen "
      f"to be for these nine rho at\nC_E0 = {CE0_probe:g}. Over {_NDRAW} random "
      f"(rho, C_E0) draws from the same ranges they are\ninexact in most "
      f"cases (here {_inexact}), and the same computation returns below 1e-14 (here up to "
      f"{M['rate_law_rho_degeneracy_random_draws_max_rel']:.1e})\n-- a few ulp of "
      "the rate, not a failure of the identity. Reported as such.")
print()
print("Companion that CAN move (the break table uses it): raise k2 by 1 % and "
      "leave C_E0 alone, so r_m is no longer held fixed --")
M["rate_law_break_k2_plus_1pct"] = degeneracy_rel(1.0, CE0_probe, k2_factor=1.01)
print(f"   max |dr/r| = {M['rate_law_break_k2_plus_1pct']:.4e}")

# %% [markdown]
# ### 4. What a batch transient sees, and at what order
#
# The rate law is blind to $\rho$; the full mass-action solution is not, because
# the quasi-steady-state approximation that produced the rate law is itself only
# approximate, and the size of its error depends on $\rho$. So the signal a batch
# run carries about $\rho$ is *the same size as the error in the model that made
# the question sensible*.
#
# Three measurements, and the order matters. **(a)** Over the whole run, where the
# maximum turns out to sit inside the induction layer - which is why it is
# root-found on a logarithmic scan and not read off a uniform one; the uniform
# scan is 66 % low. **(b)** Over the window an experiment sampled like the printed
# table can see, $t \ge 3$ min, where only about a quarter of that survives as a
# persistent offset. **(c)** The honest version of (b): the same question with
# $r_m$, $K_M$ and $C_{\rm E}^0$ all free to compensate, which is what a fit
# would do.
#
# **(c) is not measured the way (b) is, and the page does not pretend otherwise.**
# (b) is a *maximum over $t$* of the spread between two curves; (c) is an *RMS over
# a 25-point grid starting at $t = 0$* of a refitted misfit. Dividing (c) by (b)
# therefore charges to compensation a factor that is partly just the change of
# statistic and window. The second cell below measures the missing third quantity -
# the same RMS on the same grid with nothing free - so that the price of
# compensation can be quoted like for like.

# %%
te = np.linspace(0, 18, 181)
red = s_closed(te, cA0, rm_fit, KM_fit)


RHO_LO, RHO_HI = 1e-4, 1e4
eps_of = lambda CE0: CE0/(cA0 + KM_fit)


def ca_of_rho(rho, CE0):
    return ma(elem(rho, CE0, rm_fit, KM_fit), CE0, te, dense=True)


T_OBS = float(t_d[1])          # 3 min: the earliest time the printed table samples
WIN_ALL, WIN_OBS = (0.0, 18.0), (T_OBS, 18.0)


def rho_spread(CE0, lo=RHO_LO, hi=RHO_HI, window=WIN_OBS, check=False, with_t=False):
    """max |C_A(t; rho_lo) - C_A(t; rho_hi)| over `window`, ROOT-FOUND in t.

    check=True verifies on a 9-point rho grid that the ENDPOINTS really are the
    extremes, i.e. that the full-grid spread equals the endpoint spread. C_A(t)
    is NOT monotone in rho at every t - the curves cross late in the run, when
    they are all near zero - so this is checked rather than assumed.
    """
    if check:
        tw = np.linspace(*window, 181)
        Y = np.array([ma(elem(r_, CE0, rm_fit, KM_fit), CE0, tw)
                      for r_ in np.geomspace(lo, hi, 9)])
        full, ends = float(np.max(Y.max(0) - Y.min(0))), float(np.max(np.abs(Y[0] - Y[-1])))
        if abs(full - ends) > 1e-12 + 1e-9*full:
            raise RuntimeError(f"rho endpoints are not the extremes at C_E0 = {CE0}: "
                               f"{full} vs {ends}")
    return max_gap(ca_of_rho(lo, CE0), ca_of_rho(hi, CE0), window=window, with_t=with_t)


# (a) over the whole run, where the maximum sits INSIDE THE INDUCTION LAYER
v_all, t_all = rho_spread(1e-4, window=WIN_ALL, with_t=True)
M["rho_spread_induction_over_eps"] = v_all/eps_of(1e-4)
M["rho_spread_induction_peak_min"] = t_all
print(f"over the whole run at C_E0 = 1e-4: spread = {v_all:.6e} cmol/L "
      f"(= {v_all/eps_of(1e-4):.4f} eps), attained at t = {t_all:.3e} min")


def sampled_shortfall(CE0=1e-4, peak_by="root"):
    """1 - (max on the uniform 181-point grid te)/(the reported peak).

    peak_by="root" reports the root-found peak; "scan" reports the uniform-scan
    maximum itself, the defect the metric measures."""
    a, b = ca_of_rho(RHO_LO, CE0), ca_of_rho(RHO_HI, CE0)
    scan = float(np.max(np.abs(a(te) - b(te))))
    peak = rho_spread(CE0, window=WIN_ALL) if peak_by == "root" else scan
    return 1.0 - scan/peak


def induction_time(rho, CE0=1e-4, with_KM=True):
    """The induction time scale 1/(k_1 (c_A0 + K_M)), min; with_KM=False drops K_M."""
    return 1/(elem(rho, CE0, rm_fit, KM_fit)[0]*(cA0 + (KM_fit if with_KM else 0.0)))


sc_unif = float(np.max(np.abs(ca_of_rho(RHO_LO, 1e-4)(te) - ca_of_rho(RHO_HI, 1e-4)(te))))
M["sampled_max_shortfall_uniform_scan"] = sampled_shortfall()
M["induction_time_rho_lo_min"] = induction_time(RHO_LO)
M["induction_time_rho_hi_min"] = induction_time(RHO_HI)
print(f"  a UNIFORM 181-point scan of the same window reads {sc_unif:.6e}, "
      f"{100*M['sampled_max_shortfall_uniform_scan']:.1f} % low -- the whole "
      f"signal is in the first ~1e-3 min")
print(f"  induction timescale 1/(k1(c_A0+K_M)): "
      f"{M['induction_time_rho_lo_min']:.2e} min at rho = 1e-4, "
      f"{M['induction_time_rho_hi_min']:.2e} min at rho = 1e4")
print("  THIS IS INVISIBLE TO THE PRINTED TABLE, whose first sample after t = 0 "
      f"is at {T_OBS:g} min. It is what Froment p. 25 means by 'only transient\n"
      "  experimentation (stopped flow, or relaxation ...) can help'.")

# (b) over the window an experiment sampling like the printed table can see
SPREAD = {CE0: rho_spread(CE0, check=True) for CE0 in (1e-4, 1e-3, 1e-2, 1e-1)}
M["rho_spread_over_eps_trace_limit"] = SPREAD[1e-4]/eps_of(1e-4)
print(f"\nover the OBSERVABLE window t in [{T_OBS:g}, 18] min "
      f"(9-point rho grid confirms the endpoints are the extremes):")
for CE0, v in SPREAD.items():
    print(f"  C_E0 = {CE0:7.0e} cmol/L   eps = {eps_of(CE0):.4e}   "
          f"spread = {v:.6e} cmol/L   spread/eps = {v/eps_of(CE0):.4f}")



@lru_cache(maxsize=None)
def CE0_star_at(res=RES):
    """C_E0 at which the observable rho-spread reaches `res`, root-found in log10."""
    return 10**root1(lambda g: rho_spread(10**g) - res, -4, 0, n_scan=9, xtol=1e-8,
                     what="C_E0* (observable rho-spread at the resolution)")


CE0_star = CE0_star_at()
rho_spread(CE0_star, check=True)              # the endpoint check, at the root too
M["CE0_star_cmol_per_L"] = CE0_star
M["eps_star_rho_visible_at_print_resolution"] = eps_of(CE0_star)
print(f"\nROOT-FOUND in BOTH coordinates (the peak in t, then the threshold in "
      f"C_E0):\n  the observable spread reaches the {RES:g} cmol/L print "
      f"resolution at C_E0 = {CE0_star:.4e} cmol/L,\n  i.e. eps = "
      f"{eps_of(CE0_star):.4e}")
print("  -- and this holds r_m, K_M and C_E0 fixed, which no real fit does.")

# %%
# --- the honest version: can a fit at the WRONG rho reproduce the data?
tt_d = np.linspace(0, 18, 25)


def discriminate(eps, rho_true=1e4, rho_fit=1e-4, grid=None,
                 rtol=1e-11, atol=1e-14, starts=2):
    """Min RMS misfit of the rho_fit model to noise-free rho_true data,
    with r_m, K_M and C_E0 all free. An optimiser returns an UPPER BOUND on a
    minimum, so 'not separable' is the safe direction.

    rtol/atol are exposed so that a break row can loosen the ODE solve UNDERNEATH
    the root-find - a perturbation of the computation rather than of a declared
    convention, which is the only kind that can expose a wrong baseline."""
    grid = tt_d if grid is None else grid
    CE0 = eps*(cA0 + KM_fit)
    truth = ma(elem(rho_true, CE0, rm_fit, KM_fit), CE0, grid,
                        rtol=rtol, atol=atol)

    def rr(q):
        rmv, KMv, CE0v = np.exp(q)
        return ma(elem(rho_fit, CE0v, rmv, KMv), CE0v, grid,
                           rtol=rtol, atol=atol) - truth

    best = None
    for g in ([rm_fit, KM_fit, CE0], [1.2*rm_fit, 0.8*KM_fit, 0.5*CE0])[:starts]:
        r = least_squares(rr, np.log(g), xtol=1e-13, ftol=1e-13, gtol=1e-13)
        if best is None or r.cost < best.cost:
            best = r
    return float(np.sqrt(np.mean(best.fun**2))), np.exp(best.x)


print("misfit of a rho = 1e-4 model to rho = 1e4 data, everything else free:")
for eps in (1e-3, 1e-2, 3e-2, 1e-1):
    v, q = discriminate(eps)
    # the refitted constants are ill-determined (the misfit is flat in them), so
    # only the misfit is printed, at the digits that agree across BLAS kernels
    print(f"  eps = {eps:.0e}  RMS = {v:.2e} cmol/L")
rev = discriminate(1e-1, rho_true=1e-4, rho_fit=1e4)[0]
print(f"  reverse direction at eps = 1e-1: {rev:.2e} -- larger, so the direction "
      f"reported is the binding one")

# EVERY root-find on this threshold uses this one tolerance, in log10(eps):
# ln(10)*XTOL_LOG10 relative per root-find, twice that for a difference of two.
XTOL_LOG10 = 1e-6


@lru_cache(maxsize=None)
def eps_star_comp(res=RES, grid_n=25, rtol=1e-11, atol=1e-14, starts=2):
    """eps at which the compensated misfit reaches `res`, root-found in log10(eps)."""
    grid = np.linspace(0, 18, grid_n)
    return 10**root1(lambda g: discriminate(10**g, grid=grid, rtol=rtol, atol=atol,
                                            starts=starts)[0] - res,
                     -2, np.log10(0.2), n_scan=6, xtol=XTOL_LOG10,
                     what="eps* with compensation")


eps_star_c = eps_star_comp()
M["eps_star_discrimination_with_compensation"] = eps_star_c
print(f"\nROOT-FOUND: with the other three constants free, telling rho = 1e-4 from "
      f"rho = 1e4\n  at the {RES:g} cmol/L resolution needs eps = {eps_star_c:.4e}"
      f"  (C_E0 = {eps_star_c*(cA0+KM_fit):.4e} cmol/L,\n  root-found at xtol = "
      f"{XTOL_LOG10:g} in log10, i.e. {np.log(10)*XTOL_LOG10:.2e} relative).")

# ---- LIKE FOR LIKE. eps_star_rho_visible_at_print_resolution is a MAX over
# t in [3, 18] of the rho-spread; eps_star_discrimination_with_compensation is an
# RMS over tt_d (25 points, from t = 0) of a REFITTED misfit. Dividing one by the
# other charges the whole ratio to compensation, which is wrong: part of it is
# the change of statistic and window. So the SAME statistic on the SAME grid is
# measured here with rho pinned and NOTHING free, and that is the baseline the
# compensation penalty is quoted against.
def spread_rms_same_grid(eps, lo=RHO_LO, hi=RHO_HI):
    """RMS over tt_d of |C_A(rho=lo) - C_A(rho=hi)|; r_m, K_M, C_E0 all pinned."""
    CE0 = eps*(cA0 + KM_fit)
    a = ma(elem(lo, CE0, rm_fit, KM_fit), CE0, tt_d)
    b = ma(elem(hi, CE0, rm_fit, KM_fit), CE0, tt_d)
    return float(np.sqrt(np.mean((a - b)**2)))


@lru_cache(maxsize=None)
def eps_star_same_grid(res=RES, lo=RHO_LO, hi=RHO_HI):
    """eps at which the pinned RMS spread reaches `res`, root-found in log10(eps)."""
    return 10**root1(lambda g: spread_rms_same_grid(10**g, lo, hi) - res,
                     -5, -0.5, n_scan=12, xtol=1e-10,
                     what="eps* on the same grid, nothing free")


eps_star_ll = eps_star_same_grid()
M["eps_star_rms_25pt_grid_no_compensation"] = eps_star_ll
M["eps_star_compensation_penalty"] = eps_star_c/eps_star_ll
M["eps_star_statistic_and_window_factor"] = \
    eps_star_ll/M["eps_star_rho_visible_at_print_resolution"]
M["eps_star_raw_ratio_unlike_statistics"] = \
    eps_star_c/M["eps_star_rho_visible_at_print_resolution"]
print(f"\nTHE THREE EPSILONS, AND WHAT EACH ONE IS. All at the same {RES:g} cmol/L "
      f"resolution\nand the same rho pair [{RHO_LO:g}, {RHO_HI:g}]:")
print(f"  (a) max over t in [{T_OBS:g}, 18] of the spread, nothing free   "
      f"eps* = {M['eps_star_rho_visible_at_print_resolution']:.4e}")
print(f"  (b) RMS over the 25-point grid from t = 0, nothing free  "
      f"eps* = {eps_star_ll:.4e}")
print(f"  (c) the same RMS on the same grid, r_m/K_M/C_E0 free    "
      f"eps* = {eps_star_c:.4e}")
print(f"  (b)/(a) = {M['eps_star_statistic_and_window_factor']:.3f}  <- change of "
      f"STATISTIC and WINDOW only, no compensation")
print(f"  (c)/(b) = {M['eps_star_compensation_penalty']:.2f}  <- THE PRICE OF "
      f"COMPENSATION, like for like")
print(f"  (c)/(a) = {M['eps_star_raw_ratio_unlike_statistics']:.2f}  <- the raw "
      f"ratio of the two headline numbers; it is the product of the two above,\n"
      f"            and quoting it as the cost of compensation overstates that "
      f"cost by "
      f"{M['eps_star_statistic_and_window_factor']:.2f}x")

# %% [markdown]
# ### 5. The price of each reduction against the full mass-action solution
#
# Froment states when the pseudo steady state fails (*"in particular if the enzyme
# concentration is relatively large compared to that of the substrate"*, p. 25) and
# states the condition for route 1 only as words (*"the decomposition of the
# complex A-E into the product P and the enzyme E is the rate determining step"*).
# Both are turned into numbers here.

# %%
red_f = ReducedCurve(cA0, rm_fit, KM_fit)


def qssa_err(CE0, rho=1.0, route="BH", KM=None, total=False, rtol=1e-11, atol=1e-14):
    """max_t |C_A(full mass action) - C_A(reduced)|, cmol/L, at a root in t.

    route="BH" reduces with K_M = (k_-1 + k_2)/k_1 (eq. 1.5.1-14), route="MM" with
    K_M = k_-1/k_1 (eq. 1.5.1-8). KM sets the true K_M (default: the fitted one);
    total=True compares with C_A + C_A-E."""
    KMv = KM_fit if KM is None else KM
    p = elem(rho, CE0, rm_fit, KMv)
    KM_red = KMv if route == "BH" else p[1]/p[0]
    return max_gap(ma(p, CE0, te, dense=True, total=total, rtol=rtol, atol=atol),
                   ReducedCurve(cA0, rm_fit, KM_red), interp_tol=max(1e-8, 1e3*rtol*cA0))


QSSA_CE0 = (1e-5, 1e-4, 1e-3, 1e-2)          # the decades the order is read on


def qssa_order(route="BH", rtol=1e-11, atol=1e-14):
    """Order in eps of the quasi-steady-state error at rho = 1, decade to decade.

    The error is measured against the reduced model, which is the exact eps -> 0
    limit; the order is log10 of successive error ratios. margin is the smallest
    error over the integrator's floor rtol*c_A0 (the round-off companion)."""
    e = np.array([qssa_err(c, 1.0, route, rtol=rtol, atol=atol) for c in QSSA_CE0])
    o = np.log10(e[1:]/e[:-1])
    return dict(mean=float(np.mean(o)), min=float(o.min()), max=float(o.max()),
                margin=float(e.min()/(rtol*cA0)), err=e)


QSSA = {CE0: {r_: qssa_err(CE0, r_) for r_ in (1e-2, 1.0, 1e2)}
        for CE0 in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1)}
ks = sorted(QSSA)
QO = qssa_order()
M["qssa_error_order_in_eps"] = QO["mean"]
M["qssa_order_roundoff_margin"] = QO["margin"]
M["qssa_error_over_eps_rho1"] = QSSA[1e-5][1.0]/eps_of(1e-5)
print("max |C_A(full) - C_A(eq. 1.5.1-13)|, cmol/L, divided by eps:")
for CE0 in ks:
    print(f"  C_E0 = {CE0:7.0e}  eps = {eps_of(CE0):.4e}   " +
          "   ".join(f"rho={r_:<6g} {QSSA[CE0][r_]/eps_of(CE0):.4f}"
                     for r_ in (1e-2, 1.0, 1e2)))
print(f"observed order in eps at rho = 1, C_E0 = {QSSA_CE0[0]:g} to {QSSA_CE0[-1]:g}: "
      f"mean {QO['mean']:.4f}, decade by decade {QO['min']:.4f} to {QO['max']:.4f} "
      f"(first order, as the error being O(eps) requires); the smallest error is "
      f"{QO['margin']:.1e} times the integrator floor rtol*c_A0")

# the same error evaluated at the enzyme loading a DISCRIMINATING experiment needs
_ce = M["eps_star_discrimination_with_compensation"]*(cA0 + KM_fit)
qssa_at_eps_star = qssa_err(_ce)
M["qssa_error_at_discrimination_eps_over_resolution"] = qssa_at_eps_star/RES
M["qssa_error_at_discrimination_eps_cmol_per_L"] = qssa_at_eps_star
print(f"\nAT eps* = {M['eps_star_discrimination_with_compensation']:.4e}, where "
      f"rho first becomes visible with the other constants free, the pseudo "
      f"steady state\n  is itself off by {qssa_at_eps_star:.5f} cmol/L = "
      f"{M['qssa_error_at_discrimination_eps_over_resolution']:.1f} x the "
      f"{RES:g} resolution the discrimination is read at.")
# WHICH substrate the reduced model is compared against is a convention, and it
# moves this number, so both are printed. The reduced model has one substrate
# variable; the full model has free C_A and bound C_A-E.
qssa_at_eps_star_tot = qssa_err(_ce, total=True)
M["qssa_error_at_discrimination_eps_total_over_resolution"] = qssa_at_eps_star_tot/RES
print(f"  That compares eq. (1.5.1-15) with the FREE substrate C_A. Against the "
      f"TOTAL\n  substrate C_A + C_A-E - the other defensible convention - the "
      f"same error is\n  {qssa_at_eps_star_tot:.5f} cmol/L = "
      f"{qssa_at_eps_star_tot/RES:.1f} x the resolution. Both are far above 1, "
      f"which is the claim;\n  the factor quoted elsewhere on this page is the "
      f"free-C_A one.")

# which group? fix C_E0/c_A0 and sweep K_M over three decades
def qssa_group(KMs=(0.05, 0.5, 5.0, 50.0), ratio=1e-3, route="BH"):
    """At fixed C_E0/c_A0, the QSSA error over each candidate group, per K_M."""
    G = {}
    for KMv in KMs:
        CE0 = ratio*cA0
        e = qssa_err(CE0, KM=KMv, route=route)
        G[KMv] = (e, e/ratio, e/(CE0/(cA0 + KMv)))
    return dict(G=G, spread_cA0=max(g[1] for g in G.values())/min(g[1] for g in G.values()),
                spread_cA0_KM=max(g[2] for g in G.values())/min(g[2] for g in G.values()))


QG = qssa_group()
GROUP = QG["G"]
M["qssa_group_spread_CE0_over_cA0"] = QG["spread_cA0"]
M["qssa_group_spread_CE0_over_cA0_plus_KM"] = QG["spread_cA0_KM"]
M["qssa_error_over_eps_KM_0p05"] = GROUP[0.05][2]
print(f"\nat FIXED C_E0/c_A0 = 1e-3, sweeping K_M over three decades:")
for KMv, (e, a, b) in GROUP.items():
    print(f"  K_M = {KMv:6.2f}  error = {e:.4e}   error/(C_E0/c_A0) = {a:.4e}   "
          f"error/(C_E0/(c_A0+K_M)) = {b:.4f}")
print(f"  spread of the ratio: C_E0/c_A0 x{M['qssa_group_spread_CE0_over_cA0']:.1f}, "
      f"C_E0/(c_A0+K_M) x{M['qssa_group_spread_CE0_over_cA0_plus_KM']:.2f}")
print("  Froment's printed comparison is with the substrate alone, which is the "
      f"K_M << c_A0 corner of eps and is right there (coefficient "
      f"{M['qssa_error_over_eps_KM_0p05']:.4f} at K_M = 0.05); the K_M term is what "
      f"this page adds.")

# %%
# --- route 1 (rapid equilibrium, K_M = k_-1/k_1) priced against the full model
CE0t = 1e-5                       # trace enzyme, so the QSSA error is negligible


def eq_route_err(rho, **kw):
    """max_t |C_A(full) - C_A(route 1)|, the maximum at a root in t."""
    return qssa_err(CE0t, rho, route="MM", **kw)


@lru_cache(maxsize=None)
def rho_star_at(res=RES, rtol=1e-11, atol=1e-14, total=False):
    """rho at which route 1 mispredicts C_A(t) by `res`, root-found in log10(rho).

    total=True compares route 1 with C_A + C_A-E instead of the free substrate."""
    return 10**root1(lambda g: eq_route_err(10**g, rtol=rtol, atol=atol, total=total) - res,
                     -5, 0, n_scan=11, xtol=1e-10, what="rho* of route 1")


print("route 1 predicts C_A(t) from the elementary constants via K_M = k_-1/k_1:")
for r_ in (1e-3, 1e-2, 1e-1, 1.0, 10.0):
    p = elem(r_, CE0t, rm_fit, KM_fit)
    print(f"  rho = {r_:8.3g}   K_M(MM) = {p[1]/p[0]:.6f}  "
          f"max |dC_A| = {eq_route_err(r_):.4e}   (route 2: "
          f"{max_gap(ma(p, CE0t, te, dense=True), red_f):.4e})")

rho_star = rho_star_at()
M["rho_star_equilibrium_route"] = rho_star

# --- SECOND, INDEPENDENT ROUTE to the same number: analytic sensitivity of the
#     closed-form solution to K_M, with the extremum ROOT-FOUND, not sampled.
s_sym, s0_sym, KM_sym = sp.symbols("s s_0 K_M", positive=True)
g_expr = s_sym*sp.log(s0_sym/s_sym)/(s_sym + KM_sym)      # |ds/dK_M| along the curve
g_f = sp.lambdify((s_sym, s0_sym, KM_sym), g_expr, "numpy")
dg_f = sp.lambdify((s_sym, s0_sym, KM_sym), sp.diff(g_expr, s_sym), "numpy")


def dsdKM_max(sampled_on=None):
    """(s at the maximum, max) of |ds/dK_M| along the fitted curve.

    The maximum is the single root of the analytic d/ds on (0, c_A0); with
    sampled_on (a time grid) it is read off the closed-form curve there instead."""
    if sampled_on is not None:
        sv = np.clip(s_closed(sampled_on, cA0, rm_fit, KM_fit), 1e-12, None)
        gv = g_f(sv, cA0, KM_fit)
        return float(sv[np.argmax(gv)]), float(np.max(gv))       # the sampled defect
    s_ = root1(lambda s: dg_f(s, cA0, KM_fit), 1e-9, cA0*(1 - 1e-12), n_scan=64,
               xtol=1e-14, what="maximum of |ds/dK_M|")
    return s_, float(g_f(s_, cA0, KM_fit))


def rho_star_analytic(sampled_on=None, linear=False):
    """rho* from the analytic sensitivity: dK_M = K_M rho/(1+rho) meets RES.

    ``linear=True`` takes dK_M = K_M rho instead, dropping the 1/(1+rho).
    """
    g = dsdKM_max(sampled_on)[1]
    if linear:
        return RES/(g*KM_fit)
    return RES/(g*KM_fit - RES)


s_pk, gmax = dsdKM_max()
M["dsdKM_peak_s_cmol_per_L"] = s_pk
M["rho_star_equilibrium_analytic"] = rho_star_analytic()
M["rho_star_equilibrium_two_route_rel_diff"] = abs(
    M["rho_star_equilibrium_analytic"]/rho_star - 1)
s_samp = dsdKM_max(sampled_on=te)[1]
M["dsdKM_rootfound_over_sampled"] = gmax/s_samp

print(f"\nROOT-FOUND: route 1 reaches the {RES:g} cmol/L resolution at "
      f"rho = k_2/k_-1 = {rho_star:.6e}")
print(f"SECOND ROUTE (analytic sensitivity, no ODE solver, no mass-action model):")
print(f"  ds/dK_M = s ln(s_0/s)/(s + K_M); its maximum is root-found at "
      f"s = {s_pk:.6f}, value {gmax:.6f}")
print(f"  giving rho* = {M['rho_star_equilibrium_analytic']:.6e}, "
      f"{100*M['rho_star_equilibrium_two_route_rel_diff']:.3f} % from the "
      f"nonlinear root-find")
print(f"  (sampling that maximum on the 181-point t grid instead would give "
      f"{s_samp:.6f}, a factor {M['dsdKM_rootfound_over_sampled']:.6f} -- so on "
      f"THIS curve the sample happens to be right, which the root-find is what "
      f"establishes)")

# %%
fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.1))
epsv = np.array([eps_of(c) for c in ks])
for r_, col, mk in ((1e-2, C_BLUE, "o"), (1.0, C_ORANGE, "s"), (1e2, C_GREEN, "^")):
    ax[0].loglog(epsv, [QSSA[c][r_] for c in ks], mk + "-", ms=4, lw=1.4, color=col,
                 label=rf"$\rho$ = {r_:g}")
ax[0].loglog(epsv, epsv, ls=":", lw=0.9, color=C_GREY, label=r"slope 1")
ax[0].axhline(RES, color=C_PURPLE, ls="--", lw=1.0)
ax[0].text(epsv[0], RES*1.25, "print resolution", fontsize=7, color=C_PURPLE)
ax[0].set_xlabel(r"$\varepsilon = C_{\rm E}^0/(c_{\rm A0}+K_M)$")
ax[0].set_ylabel(r"max $|C_{\rm A}$(full) $- C_{\rm A}$(1.5.1-13)$|$")
ax[0].set_title("cost of the pseudo steady state", fontsize=9)
ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=0.25, lw=0.5, which="both")

rr_g = np.geomspace(1e-4, 1e1, 60)
ax[1].loglog(rr_g, [eq_route_err(r_) for r_ in rr_g], color=C_ORANGE, lw=1.6,
             label="route 1, $K_M = k_{-1}/k_1$")
ax[1].loglog(rr_g, [max_gap(ca_of_rho(r_, CE0t), red_f) for r_ in rr_g],
             color=C_BLUE, lw=1.6, label="route 2, $K_M = (k_{-1}+k_2)/k_1$")
ax[1].axhline(RES, color=C_PURPLE, ls="--", lw=1.0)
ax[1].axvline(rho_star, color=C_GREY, ls=":", lw=1.0)
ax[1].text(rho_star*1.15, 2e-5, rf"$\rho^* = {rho_star:.2e}$", fontsize=7.5)
ax[1].set_xlabel(r"$\rho = k_2/k_{-1}$")
ax[1].set_ylabel(r"max $|\Delta C_{\rm A}|$ (cmol/L)")
ax[1].set_title(r"at trace enzyme, $\varepsilon$ = %.1e" % eps_of(CE0t), fontsize=9)
ax[1].legend(fontsize=7.5, frameon=False); ax[1].grid(alpha=0.25, lw=0.5, which="both")
fig.tight_layout(); plt.show()


# %% [markdown]
# ### 6. What the actual seven points say about $\rho$
#
# Pin $\rho$, let $r_m$, $K_M$ and $C_{\rm E}^0$ re-optimise against the printed
# table, and read the best achievable sum of squares. If $\rho$ mattered to the
# data there would be an interior optimum and a sharp one.

# %%
def sse_pin_rho(rho, y=c_d):
    """Best SSE of the mass-action model to (t_d, y) with rho pinned, and its constants."""
    def rr(q):
        rmv, KMv, CE0v = np.exp(q)
        return ma(elem(rho, CE0v, rmv, KMv), CE0v, t_d, a0=y[0]) - y
    best = None
    for CE0g in (1e-4, 1e-2):
        r = least_squares(rr, np.log([rm_fit, KM_fit, CE0g]), jac="3-point",
                          diff_step=1e-5, xtol=1e-13, ftol=1e-13, gtol=1e-13)
        if best is None or r.cost < best.cost:
            best = r
    return 2*float(best.cost), np.exp(best.x)


RHO_PROFILE = (1e-4, 1e-2, 1.0, 1e2, 1e4)


def f_crit(level=0.95, n=len(t_d)):
    """The `level` point of F(2, n - 4)."""
    return float(f_dist.ppf(level, 2, n - 4))


def rho_profile(y=c_d, rhos=RHO_PROFILE):
    """SSE span over the pinned rho, and the F statistic of the four-parameter model."""
    sses = np.array([sse_pin_rho(r_, y)[0] for r_ in rhos])
    sse_mm = float(np.sum(fit_mm(t_d, y)[3]**2))
    best = float(sses.min())
    return dict(span=float(sses.max()/sses.min() - 1),
                F=float(((sse_mm - best)/2)/(best/(len(t_d) - 4))))


PROF = {r_: sse_pin_rho(r_) for r_ in RHO_PROFILE}
for r_, (s_, q_) in PROF.items():
    # the refitted constants are ill-determined; the SSE is printed at the
    # digits that agree across BLAS kernels
    print(f"  rho = {r_:8.0e}   SSE = {s_:.3e}")
sses = np.array([v[0] for v in PROF.values()])
M["rho_profile_sse_span_fraction"] = float(sses.max()/sses.min() - 1)
sse_best = float(sses.min())
M["full_model_F_statistic"] = float(((sse_fit - sse_best)/2)/(sse_best/(len(t_d) - 4)))
M["full_model_F_crit_95"] = f_crit()
print(f"\nSSE span across eight decades of rho: "
      f"{100*M['rho_profile_sse_span_fraction']:.2f} %, monotone, no interior optimum")
print(f"two-parameter Michaelis-Menten SSE {sse_fit:.6e}; best four-parameter "
      f"mass-action SSE {sse_best:.6e}")
print(f"F = {M['full_model_F_statistic']:.4f} against F(2, {len(t_d)-4}) = "
      f"{M['full_model_F_crit_95']:.4f} at 95 % -- the two extra parameters are "
      f"NOT justified by these data")
print("\nONE CAVEAT ON THE F TEST, STATED BECAUSE IT IS A REAL ONE. The models "
      "are nested\n(the reduced model is the C_E0 -> 0 limit at fixed r_m, "
      "K_M), but the nesting is\nNON-REGULAR: the null sits on the boundary "
      "C_E0 = 0 of the parameter space and rho\nis unidentified under it, so "
      "F(2, 3) is an approximation to the null distribution\nrather than the "
      "exact one. The direction is safe. A boundary null makes the true\nnull "
      "distribution stochastically no larger, so the real p-value is at least "
      "the\ntabulated one and 'not justified' holds a fortiori.")

# %% [markdown]
# ## Validation
#
# Five kinds, kept apart:
#
# 1. **Structural identities** (the ten symbolic ones, the symbolic proof of the
#    $\rho$-degeneracy, the exact linearity of $1/r$ in $1/C_{\rm A}$). These are
#    proofs; they cannot fail for the right reason, the numerical ones sit at
#    round-off and are registered as bounds, and each is named with a companion
#    that a break row can move. The *numerical* zero that goes with the degeneracy is reported
#    beside its own float spread over random draws, because an exact 0.0 in
#    floating point is a property of the values chosen, not of the algebra.
# 2. **Conservation** in the pymrm mass-action marcher: the two balances the
#    scheme implies, on a model that integrates all four fields independently, so
#    a wrong stoichiometric row would break them.
# 3. **Time-step refinement with observed orders** from successive differences
#    (`pymrm.checks.observed_orders`) and their round-off margin, plus Richardson
#    extrapolation checked against a different solver family.
# 4. **A second, independent computation of a headline** - two of them, in fact:
#    the fit through the pymrm marcher instead of the Lambert $W$ closed form,
#    and $\rho^*$ from an analytic sensitivity instead of a nonlinear root-find.
# 5. **A defect-injection table**: every non-structural claim has a row that
#    fails its own check, each defect injected through an argument of the study.
#    A single least-squares start instead of two is not a defect; its move of the
#    compensated $\varepsilon^*$ is a bound claim of its own.

# %%
# --- pymrm backward Euler: order in dt, conservation, and a second solver family
CE0h, n_out = 1e-2, 18
te_h = np.linspace(0, 18, n_out + 1)
p_h = elem(1.0, CE0h, rm_fit, KM_fit)
ref_lsoda = ma(p_h, CE0h, te_h)
ref_radau = ma(p_h, CE0h, te_h, method="Radau", rtol=1e-12, atol=1e-15)
M["lsoda_vs_radau_reference"] = float(np.max(np.abs(ref_lsoda - ref_radau)))

BE_LADDER = (180, 360, 720, 1440)
C_INIT = [cA0, CE0h, 0., 0.]


@lru_cache(maxsize=None)
def be_study(theta=1.0, tol=1e-13, ladder=BE_LADDER, extrapolate=True, **model_kw):
    """Time-step study of BatchMassAction: orders, Richardson, balances, round-off margin.

    The orders come from successive differences of C_A at the 19 output times
    (pymrm.checks.observed_orders); margin is the finest difference over the
    accumulated newton step tolerance n_steps*tol."""
    runs = {n: BatchMassAction(p_h, theta=theta, tol=tol, **model_kw).run(C_INIT, 18., n, n_out)
            for n in ladder}
    st = observed_orders(lambda n: runs[n][:, 0], ladder)
    o = np.array(st["orders"])
    f1, f2 = runs[ladder[-2]], runs[ladder[-1]]
    return dict(
        runs=runs, orders=o, mean=float(np.mean(o)),
        worst=float(max(o, key=lambda v: abs(v - 1.0))),
        margin=float(np.max(np.abs(f2[:, 0] - f1[:, 0]))/(ladder[-1]*tol)),
        err={n: float(np.max(np.abs(runs[n][:, 0] - ref_radau))) for n in ladder},
        richardson=float(np.max(np.abs((2*f2[:, 0] - f1[:, 0] if extrapolate
                                        else runs[ladder[0]][:, 0]) - ref_radau))),
        enzyme=float(np.max(np.abs(f2[:, 1] + f2[:, 2] - CE0h))),
        substrate=float(np.max(np.abs(f2[:, 0] + f2[:, 2] + f2[:, 3] - cA0))))


BE = be_study()
M["pymrm_be_time_order"] = BE["worst"]
M["pymrm_be_order_roundoff_margin"] = BE["margin"]
M["pymrm_richardson_vs_radau"] = BE["richardson"]
M["pymrm_enzyme_balance_max"] = BE["enzyme"]
M["pymrm_substrate_balance_max"] = BE["substrate"]
for i, n in enumerate(BE_LADDER):
    o = "" if i < 2 else f"   order from successive differences {BE['orders'][i-2]:.4f}"
    print(f"  n_steps = {n:5d}  dt = {18/n:.5f} min  max |dC_A - Radau| = {BE['err'][n]:.5e}{o}")
print(f"\n  observed orders {BE['orders'].min():.4f} to {BE['orders'].max():.4f} "
      f"(backward Euler: 1); finest difference {BE['margin']:.1e} times the "
      f"accumulated newton tolerance")
print(f"  Richardson(720, 1440) vs Radau      {M['pymrm_richardson_vs_radau']:.4e}")
print(f"  LSODA vs Radau on the same problem  below 1e-9 "
      f"(here {M['lsoda_vs_radau_reference']:.1e})")
print(f"  enzyme balance   max |C_E + C_A-E - C_E0|      below 1e-13 "
      f"(here {M['pymrm_enzyme_balance_max']:.0e})")
print(f"  substrate balance max |C_A + C_A-E + C_P - c_A0| below 1e-13 "
      f"(here {M['pymrm_substrate_balance_max']:.0e})")


# %%
# --- SECOND INDEPENDENT ROUTE TO THE FIT: pymrm backward Euler + Richardson,
#     sharing no algebra with the Lambert W closed form.
@lru_cache(maxsize=None)
def fit_via_pymrm(n_steps=1800, richardson=True):
    """(k, K) fitted through BatchReduced, Richardson-extrapolated from n and 2n steps."""
    def resid(q):
        kk, KK = np.exp(q)
        a = BatchReduced(kk/KK, 1/KK).run(c_d[0], 18., n_steps, t_d)
        if not richardson:
            return a - c_d
        b = BatchReduced(kk/KK, 1/KK).run(c_d[0], 18., 2*n_steps, t_d)
        return (2*b - a) - c_d
    r = least_squares(resid, np.log([0.3, 2.0]), xtol=1e-13, ftol=1e-13, gtol=1e-13)
    return np.exp(r.x)


k_alt, K_alt = fit_via_pymrm()
M["fit_second_route_k_rel_diff"] = abs(k_alt/k_fit - 1)
M["fit_second_route_K_rel_diff"] = abs(K_alt/K_fit - 1)
print(f"Lambert W closed form  : k = {k_fit:.6f}   K = {K_fit:.6f}")
print(f"pymrm BE + Richardson  : k = {k_alt:.6f}   K = {K_alt:.6f}")
print(f"relative difference    : below 1e-5 / below 1e-6 (here "
      f"{M['fit_second_route_k_rel_diff']:.1e} / "
      f"{M['fit_second_route_K_rel_diff']:.1e})")
print("\nWhat this catches that a break row cannot: an algebra error in the "
      "Lambert W inversion, or a wrong branch. Both routes integrate the same "
      "printed ODE; they share no line of code.")

# %% [markdown]
# ### The claims and the defect-injection table
#
# Every number the page states is registered once below with how it is checked
# (a value within its tolerance, a bound with its basis, an exact count). Then
# every non-structural claim gets a row that injects a defect through an
# argument of the study that computes it, and the row must fail that claim's own
# check. Structural claims (the symbolic identities, the float spread under an
# exact identity, a definitional ratio, a linear invariant of the source) are
# named as such.

# %%
# -- the fit quantities as one function of the data, so the perturbed-data row
#    recomputes every one of them through the same code
def lineweaver_printed_KM(y):
    """K_M read as -intercept/slope of 1/r against C_A (Froment's printed axes)."""
    _, _, cm, rr = lineweaver(t_d, y)
    sl, ic = np.linalg.lstsq(np.vstack([cm, np.ones_like(cm)]).T, 1/rr, rcond=None)[0]
    return float(-ic/sl)


def fit_metrics(y):
    """Every fit claim of Results 1 and 2 for the batch table y on t_d."""
    k, K, _, res = fit_mm(t_d, y)
    rms = float(np.sqrt(np.mean(res**2)))
    nulls = {kind: null_fit(kind, t_d, y)[0] for kind in ("zero", "first", "second")}
    k_l, K_l, _, _ = lineweaver(t_d, y)
    k_m, K_m, _, _ = lineweaver(t_d, s_closed(t_d, y[0], k/K, 1/K))
    return {
        "mm_fit_k_per_min": k, "mm_fit_K_L_per_cmol": K,
        "mm_fit_rm_cmol_per_L_min": k/K, "mm_fit_KM_cmol_per_L": 1/K,
        "mm_fit_rms_cmol_per_L": rms,
        "mm_fit_maxres_over_print_resolution": float(np.max(np.abs(res))/RES),
        "mm_fit_sse_penalty_K_pinned_at_2":
            float(np.sum(fit_mm(t_d, y, K_fixed=2.0)[3]**2)/np.sum(res**2) - 1.0),
        "mm_fit_KM_with_s0_free": 1.0/fit_mm(t_d, y, s0_free=True)[1],
        "mm_fit_KM_rel_dev_from_half": abs(1/K/0.5 - 1),
        "null_first_order_rms_ratio": nulls["first"]/rms,
        "null_zero_order_rms_ratio": nulls["zero"]/rms,
        "null_second_order_rms_ratio": nulls["second"]/rms,
        "null_first_order_rms_cmol_per_L": nulls["first"],
        "null_zero_order_rms_cmol_per_L": nulls["zero"],
        "null_second_order_rms_cmol_per_L": nulls["second"],
        "lb_rm_bias_vs_nls": (k_l/K_l)/(k/K) - 1,
        "lb_KM_bias_vs_nls": (1/K_l)/(1/K) - 1,
        "lb_KM_bias_finite_difference_only": (1/K_m)*K - 1,
    }


FIT_BASE = fit_metrics(c_d)
_drift = {k: (v, M[k]) for k, v in FIT_BASE.items() if not np.isclose(v, M[k], rtol=1e-12, atol=0)}
if _drift:
    raise RuntimeError(f"fit_metrics does not reproduce the cells above: {_drift}")

# -- the quantities only a sensitivity study computes, registered so that every
#    number quoted in the text is a claim
M["eps_star_comp_res_doubled"] = eps_star_comp(res=2*RES)
M["eps_star_comp_grid13"] = eps_star_comp(grid_n=13)
M["eps_star_rms_rho_narrow"] = eps_star_same_grid(lo=0.5, hi=2.0)
M["eps_star_single_start_rel_move"] = abs(eps_star_comp(starts=1)/eps_star_c - 1)
print(f"eps* with compensation: resolution doubled {M['eps_star_comp_res_doubled']:.3e}, "
      f"13-point grid {M['eps_star_comp_grid13']:.3e}; a single start moves it by below 1e-5 "
      f"(here {M['eps_star_single_start_rel_move']:.1e}), the order of the "
      f"{2*np.log(10)*XTOL_LOG10:.2e} two root-finds at xtol {XTOL_LOG10:g} in log10 resolve")
print(f"like-for-like eps* with the rho pair narrowed to [0.5, 2]: "
      f"{M['eps_star_rms_rho_narrow']:.4e}")

# %%
# Every number is registered ONCE with how it is checked (docs/claims-design.md).
# Round-off-level numbers are bounds; their digits differ between BLAS kernels.

# -- structure: proofs and their numerical companions
claims.bound("froment_chain_symbolic_max_residual", M["froment_chain_symbolic_max_residual"],
             max=0.0, basis="identity", structural=True,
             note="ten sympy identities simplify to 0; protects the transcription")
claims.exact("froment_chain_identities_verified", int(M["froment_chain_identities_verified"]),
             structural=True, note="the number of identities the code checks; it changes "
             "only if the code does")
claims.bound("rho_degeneracy_symbolic_residual", M["rho_degeneracy_symbolic_residual"],
             max=0.0, basis="identity", structural=True,
             note="sympy d r/d rho of eq. (1.5.1-13) under the reparameterisation")
claims.bound("rate_law_rho_degeneracy_max_rel", M["rate_law_rho_degeneracy_max_rel"],
             max=1e-14, basis="identity",
             note="relative to the rate; exactly 0 where the float round trips are exact")
claims.bound("rate_law_rho_degeneracy_random_draws_max_rel",
             M["rate_law_rho_degeneracy_random_draws_max_rel"], max=1e-14, basis="identity",
             structural=True, note="relative to the rate; a few ulp over 2000 draws, "
             "the arithmetic under an exact identity")
claims.exact("KM_ratio_at_rho_1e4", M["KM_ratio_at_rho_1e4"], structural=True,
             note="definitional: 1 + rho at rho = 1e4")
claims.value("rate_law_break_k2_plus_1pct", M["rate_law_break_k2_plus_1pct"],
             note="the above-floor companion of the degeneracy: r_m no longer held fixed")
claims.bound("lineweaver_correct_axes_nonlinearity", M["lineweaver_correct_axes_nonlinearity"],
             max=1e-13, basis="identity", note="fraction of the range of 1/r; 1/r is exactly "
             "affine in 1/C_A, measured 4.8e-16 to 1.2e-15 across BLAS kernels")

# -- the fit to the seven printed points. The tolerances are tight so that the
#    one-print-ulp row (a 1e-4 relative shift) fails them; the closed-form fit is
#    reproducible far below that.
claims.value("mm_fit_k_per_min", M["mm_fit_k_per_min"], rtol=1e-5, unit="1/min")
claims.value("mm_fit_K_L_per_cmol", M["mm_fit_K_L_per_cmol"], rtol=1e-5, unit="L/cmol")
claims.value("mm_fit_rm_cmol_per_L_min", M["mm_fit_rm_cmol_per_L_min"], rtol=1e-5,
             unit="cmol/L/min")
claims.value("mm_fit_KM_cmol_per_L", M["mm_fit_KM_cmol_per_L"], rtol=1e-5, unit="cmol/L")
claims.value("mm_fit_rms_cmol_per_L", M["mm_fit_rms_cmol_per_L"], rtol=1e-3, unit="cmol/L")
claims.value("mm_fit_maxres_over_print_resolution",
             M["mm_fit_maxres_over_print_resolution"], rtol=1e-3)
claims.value("mm_fit_sse_penalty_K_pinned_at_2", M["mm_fit_sse_penalty_K_pinned_at_2"],
             rtol=1e-3)
claims.value("mm_fit_KM_with_s0_free", M["mm_fit_KM_with_s0_free"], rtol=1e-5,
             unit="cmol/L")
claims.value("mm_fit_KM_rel_dev_from_half", M["mm_fit_KM_rel_dev_from_half"], rtol=1e-3,
             note="|K_M/0.5 - 1|")
claims.value("null_first_order_rms_ratio", M["null_first_order_rms_ratio"], rtol=1e-3)
for _k, _u in (("null_zero_order_rms_ratio", None), ("null_second_order_rms_ratio", None),
               ("null_first_order_rms_cmol_per_L", "cmol/L"),
               ("null_zero_order_rms_cmol_per_L", "cmol/L"),
               ("null_second_order_rms_cmol_per_L", "cmol/L")):
    claims.value(_k, M[_k], rtol=1e-5, unit=_u)
claims.value("lb_rm_bias_vs_nls", M["lb_rm_bias_vs_nls"], rtol=1e-5)
claims.value("lb_KM_bias_vs_nls", M["lb_KM_bias_vs_nls"], rtol=1e-5)
claims.value("lb_KM_bias_finite_difference_only", M["lb_KM_bias_finite_difference_only"],
             rtol=1e-4)

# -- straightness of the printed Lineweaver axes (closed form)
claims.value("lineweaver_printed_axes_nonlinearity_minimax",
             M["lineweaver_printed_axes_nonlinearity_minimax"])
claims.value("lineweaver_printed_axes_nonlinearity", M["lineweaver_printed_axes_nonlinearity"])
claims.value("lineweaver_printed_axes_nonlinearity_uniform_measure",
             M["lineweaver_printed_axes_nonlinearity_uniform_measure"])

# -- what a batch transient sees, and the epsilons
for _k in ("rho_spread_induction_over_eps", "rho_spread_induction_peak_min",
           "sampled_max_shortfall_uniform_scan", "rho_spread_over_eps_trace_limit",
           "induction_time_rho_lo_min", "induction_time_rho_hi_min",
           "CE0_star_cmol_per_L", "eps_star_rho_visible_at_print_resolution",
           "eps_star_discrimination_with_compensation",
           "eps_star_rms_25pt_grid_no_compensation", "eps_star_compensation_penalty",
           "eps_star_statistic_and_window_factor", "eps_star_raw_ratio_unlike_statistics",
           "eps_star_comp_res_doubled", "eps_star_comp_grid13", "eps_star_rms_rho_narrow",
           "qssa_error_over_eps_rho1", "qssa_error_at_discrimination_eps_over_resolution",
           "qssa_error_at_discrimination_eps_cmol_per_L",
           "qssa_error_at_discrimination_eps_total_over_resolution",
           "qssa_group_spread_CE0_over_cA0", "qssa_group_spread_CE0_over_cA0_plus_KM",
           "qssa_error_over_eps_KM_0p05"):
    claims.value(_k, M[_k])
claims.bound("rho_star_equilibrium_two_route_rel_diff",
             M["rho_star_equilibrium_two_route_rel_diff"], max=1e-3, basis="chosen",
             note="the analytic route linearises in dK_M = K_M rho/(1+rho); measured 8.0e-4, "
             "its fifth digit moves 2e-8 across BLAS kernels with the LSODA steps under "
             "the rho* root-find")
claims.bound("eps_star_single_start_rel_move", M["eps_star_single_start_rel_move"], max=1e-5,
             basis="tolerance", note="a single least-squares start against two, the order "
             "of the 4.6e-6 two root-finds at xtol 1e-6 in log10 resolve")
claims.bound("qssa_error_order_in_eps", M["qssa_error_order_in_eps"], min=0.9, max=1.1,
             basis="theory", note="the quasi-steady-state error is first order in eps "
             "(mean of the three decade-to-decade orders at rho = 1)")
claims.bound("qssa_order_roundoff_margin", M["qssa_order_roundoff_margin"], min=1e3,
             basis="tolerance", note="smallest error of the eps ladder over rtol*c_A0")
claims.value("rho_star_equilibrium_route", M["rho_star_equilibrium_route"], rtol=1e-4)
claims.value("rho_star_equilibrium_analytic", M["rho_star_equilibrium_analytic"], rtol=1e-5)
claims.value("dsdKM_rootfound_over_sampled", M["dsdKM_rootfound_over_sampled"], rtol=1e-5)
claims.value("dsdKM_peak_s_cmol_per_L", M["dsdKM_peak_s_cmol_per_L"], rtol=1e-6,
             unit="cmol/L")
claims.value("rho_profile_sse_span_fraction", M["rho_profile_sse_span_fraction"],
             rtol=2e-3, note="spreads 1e-8 relative across BLAS kernels (multi-start "
             "refits); the tolerance is 20 times that, below the one-print-ulp row's move")
claims.value("full_model_F_statistic", M["full_model_F_statistic"], rtol=1e-3)
claims.value("full_model_F_crit_95", M["full_model_F_crit_95"],
             note="definitional: the 95 % point of F(2, 3)")

# -- numerics of the pymrm marcher and the second routes
claims.bound("lsoda_vs_radau_reference", M["lsoda_vs_radau_reference"], max=1e-9,
             basis="tolerance", unit="cmol/L",
             note="LSODA at rtol 1e-11, atol 1e-14 against Radau at rtol 1e-12, atol 1e-15; "
             "measured 3.8e-11 to 5.8e-11 across BLAS kernels")
claims.bound("pymrm_be_time_order", M["pymrm_be_time_order"], min=0.9, max=1.1,
             basis="theory", note="backward Euler is first order in dt; the observed order "
             "farthest from 1 over successive differences on n = 180 to 1440")
claims.bound("pymrm_be_order_roundoff_margin", M["pymrm_be_order_roundoff_margin"], min=1e3,
             basis="tolerance", note="finest successive difference over n_steps*newton tol")
claims.value("pymrm_richardson_vs_radau", M["pymrm_richardson_vs_radau"], unit="cmol/L",
             note="the O(dt^2) remainder of Richardson(720, 1440), a discretisation error")
claims.bound("pymrm_enzyme_balance_max", M["pymrm_enzyme_balance_max"], max=1e-13,
             basis="identity", structural=True, unit="cmol/L",
             note="a linear invariant of the source, kept by every Newton step to round-off")
claims.bound("pymrm_substrate_balance_max", M["pymrm_substrate_balance_max"], max=1e-13,
             basis="identity", unit="cmol/L",
             note="a linear invariant of the source, kept by every Newton step to round-off")
claims.bound("fit_second_route_k_rel_diff", M["fit_second_route_k_rel_diff"], max=1e-5,
             basis="chosen", note="two-route gap, Lambert W against Richardson-extrapolated "
             "backward Euler at 1800/3600 steps; measured 5.9e-8 to 1.1e-7 across BLAS "
             "kernels, the bound is about 100x that and 1e2 below the no-Richardson row")
claims.bound("fit_second_route_K_rel_diff", M["fit_second_route_K_rel_diff"], max=1e-6,
             basis="chosen", note="as for k; measured 8e-10 to 7.4e-8 across BLAS kernels, "
             "the bound is about 13x the largest and 1e3 below the no-Richardson row")
_unreg = sorted(set(M) - set(claims.metrics))
if _unreg:
    raise RuntimeError(f"computed but not registered: {_unreg}")
print(f"{len(claims.metrics)} claims registered")

# %%
c_p = c_d.copy()
c_p[3] += RES                                # one print ulp on c_S(9 min)
_eps_obs = lambda: eps_of(CE0_star_at())
_t19 = np.linspace(0, 18, 19)
PERT = "c_S(9 min) moved by +5e-4 (one print ulp)"
MM8 = "eq. (1.5.1-14) -> (1.5.1-8): K_M = k_-1/k_1 at rho = 1"
RES2 = "resolution threshold doubled to 1e-3 cmol/L"
NARROW = "rho range narrowed from [1e-4, 1e4] to [0.5, 2]"
GRID13 = "discrimination misfit read on 13 grid points instead of 25"
SAMP19 = "analytic route's max |ds/dK_M| sampled on 19 points"
ROWS = [
    ("k2 raised 1 % with C_E0 left alone", {
        "rate_law_rho_degeneracy_max_rel":
            lambda: max(degeneracy_rel(r_, CE0_probe, k2_factor=1.01) for r_ in rho_dec)}),
    ("k2 raised 10 % with C_E0 left alone", {
        "rate_law_break_k2_plus_1pct": lambda: degeneracy_rel(1.0, CE0_probe, k2_factor=1.10)}),
    ("straightness of 1/r measured against C_A", {
        "lineweaver_correct_axes_nonlinearity": lambda: straightness(c_ex, 1/r_ex)}),
    ("Lineweaver run on Froment's printed axes (1/r vs C_A)", {
        "mm_fit_KM_cmol_per_L": lambda: lineweaver_printed_KM(c_d)}),
    ("c_A0 fitted instead of read off the table", {
        "mm_fit_KM_cmol_per_L": lambda: 1.0/fit_mm(t_d, c_d, s0_free=True)[1]}),
    (PERT, {**{k: (lambda k=k: fit_metrics_p()[k]) for k in FIT_BASE},
            "rho_profile_sse_span_fraction": lambda: rho_profile(c_p)["span"],
            "full_model_F_statistic": lambda: rho_profile(c_p)["F"]}),
    ("straightness measured on 0.4-0.6 cmol/L only (a tenth of the range)", {
        "lineweaver_printed_axes_nonlinearity": lambda: straightness_exact(0.4, 0.6, "log"),
        "lineweaver_printed_axes_nonlinearity_uniform_measure":
            lambda: straightness_exact(0.4, 0.6, "uniform"),
        "lineweaver_printed_axes_nonlinearity_minimax": lambda: straightness_minimax(0.4, 0.6)}),
    ("straightness weighted uniformly instead of logarithmically", {
        "lineweaver_printed_axes_nonlinearity": lambda: straightness_exact(LW_A, LW_B, "uniform")}),
    ("the floor measured against the LEAST-SQUARES line instead of the minimax one", {
        "lineweaver_printed_axes_nonlinearity_minimax":
            lambda: straightness_exact(LW_A, LW_B, "log")}),
    ("induction spread read over the observable window instead", {
        "rho_spread_induction_over_eps": lambda: rho_spread(1e-4, window=WIN_OBS)/eps_of(1e-4),
        "rho_spread_induction_peak_min": lambda: rho_spread(1e-4, window=WIN_OBS, with_t=True)[1]}),
    ("induction peak read on a uniform 181-point scan (sampled max)", {
        "sampled_max_shortfall_uniform_scan": lambda: sampled_shortfall(peak_by="scan")}),
    ("induction time scale taken as 1/(k_1 c_A0), K_M dropped", {
        "induction_time_rho_lo_min": lambda: induction_time(RHO_LO, with_KM=False),
        "induction_time_rho_hi_min": lambda: induction_time(RHO_HI, with_KM=False)}),
    (NARROW, {
        "rho_spread_over_eps_trace_limit": lambda: rho_spread(1e-4, 0.5, 2.0)/eps_of(1e-4),
        "eps_star_rms_25pt_grid_no_compensation": lambda: eps_star_same_grid(lo=0.5, hi=2.0)}),
    (NARROW + " (like-for-like eps* only)", {
        "eps_star_statistic_and_window_factor":
            lambda: eps_star_same_grid(lo=0.5, hi=2.0)/_eps_obs()}),
    (RES2, {
        "CE0_star_cmol_per_L": lambda: CE0_star_at(2*RES),
        "eps_star_rho_visible_at_print_resolution": lambda: eps_of(CE0_star_at(2*RES)),
        "eps_star_discrimination_with_compensation": lambda: eps_star_comp(res=2*RES),
        "eps_star_rms_25pt_grid_no_compensation": lambda: eps_star_same_grid(res=2*RES),
        "eps_star_compensation_penalty":
            lambda: eps_star_comp(res=2*RES)/eps_star_same_grid(res=2*RES)}),
    (GRID13, {
        "eps_star_discrimination_with_compensation": lambda: eps_star_comp(grid_n=13),
        "eps_star_raw_ratio_unlike_statistics": lambda: eps_star_comp(grid_n=13)/_eps_obs()}),
    ("single start compared with the two-start search on the 13-point grid", {
        "eps_star_single_start_rel_move":
            lambda: abs(eps_star_comp(starts=1)/eps_star_comp(grid_n=13) - 1)}),
    ("sensitivity studies run at the baseline settings", {
        "eps_star_comp_res_doubled": lambda: eps_star_comp(res=RES),
        "eps_star_comp_grid13": lambda: eps_star_comp(grid_n=25),
        "eps_star_rms_rho_narrow": lambda: eps_star_same_grid(lo=RHO_LO, hi=RHO_HI)}),
    (MM8, {
        "qssa_error_over_eps_rho1": lambda: eq_route_err(1.0)/eps_of(CE0t),
        "qssa_error_at_discrimination_eps_over_resolution": lambda: qssa_err(_ce, route="MM")/RES,
        "qssa_error_at_discrimination_eps_cmol_per_L": lambda: qssa_err(_ce, route="MM"),
        "qssa_error_at_discrimination_eps_total_over_resolution":
            lambda: qssa_err(_ce, route="MM", total=True)/RES,
        "qssa_error_order_in_eps": lambda: qssa_order(route="MM")["mean"],
        "qssa_error_over_eps_KM_0p05":
            lambda: qssa_group(KMs=(0.05,), route="MM")["G"][0.05][2]}),
    ("integrator at rtol 1e-6 under the eps ladder", {
        "qssa_order_roundoff_margin": lambda: qssa_order(rtol=1e-6, atol=1e-9)["margin"]}),
    ("K_M swept over 0.4-0.6 instead of 0.05-50", {
        "qssa_group_spread_CE0_over_cA0": lambda: qssa_group(KMs=(0.4, 0.6))["spread_cA0"],
        "qssa_group_spread_CE0_over_cA0_plus_KM":
            lambda: qssa_group(KMs=(0.4, 0.6))["spread_cA0_KM"]}),
    ("route 1 error measured on C_A + C_A-E instead of free C_A", {
        "rho_star_equilibrium_route": lambda: rho_star_at(total=True)}),
    ("max |ds/dK_M| sampled on a 19-point t grid", {
        "dsdKM_rootfound_over_sampled": lambda: gmax/dsdKM_max(sampled_on=_t19)[1]}),
    (SAMP19, {
        "rho_star_equilibrium_analytic": lambda: rho_star_analytic(sampled_on=_t19),
        "dsdKM_peak_s_cmol_per_L": lambda: dsdKM_max(sampled_on=_t19)[0]}),
    ("analytic route with dK_M = K_M rho, the 1/(1+rho) dropped", {
        "rho_star_equilibrium_analytic": lambda: rho_star_analytic(linear=True),
        "rho_star_equilibrium_two_route_rel_diff":
            lambda: abs(rho_star_analytic(linear=True)/rho_star - 1)}),
    ("F critical value read at 90 % instead of 95 %", {
        "full_model_F_crit_95": lambda: f_crit(0.90)}),
    ("LSODA at rtol 1e-6 instead of 1e-11", {
        "lsoda_vs_radau_reference":
            lambda: float(np.max(np.abs(ma(p_h, CE0h, te_h, rtol=1e-6, atol=1e-9) - ref_radau)))}),
    ("Crank-Nicolson (theta = 1/2) in place of backward Euler", {
        "pymrm_be_time_order": lambda: be_study(theta=0.5)["worst"]}),
    ("newton step tolerance loosened to 1e-6, backward-error test off", {
        "pymrm_be_order_roundoff_margin":
            lambda: be_study(tol=1e-6, backward_tol=None)["margin"]}),
    ("P row of the pymrm source doubled", {
        "pymrm_substrate_balance_max": lambda: be_study(p_row_factor=2.0)["substrate"]}),
    ("k_-1 term deleted from the pymrm source", {
        "pymrm_richardson_vs_radau": lambda: be_study(reverse=False)["richardson"]}),
    ("coarsest run (n = 180) reported instead of Richardson", {
        "pymrm_richardson_vs_radau": lambda: be_study(extrapolate=False)["richardson"]}),
    ("second fit route run without Richardson extrapolation", {
        "fit_second_route_k_rel_diff": lambda: abs(fit_via_pymrm(richardson=False)[0]/k_fit - 1),
        "fit_second_route_K_rel_diff": lambda: abs(fit_via_pymrm(richardson=False)[1]/K_fit - 1)}),
]


@lru_cache(maxsize=None)
def _fit_metrics_p():
    return fit_metrics(c_p)


fit_metrics_p = _fit_metrics_p
claims.define_rows([r[0] for r in ROWS])
BREAK = []
for row, studies in ROWS:
    for name, fn in studies.items():
        try:
            out = float(fn())
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
print(bt.to_string(index=False, float_format=lambda v: f"{v:.4g}"))   # kernel-stable digits
idle = [r for r in bt["injected defect"].unique()
        if not bt[(bt["injected defect"] == r) & bt.caught].shape[0]]
print(f"\n{len(ROWS)} rows defined, {bt['injected defect'].nunique()} in the table; "
      f"rows that fail no claim: {idle or 'none'}")
missing = [k for k in claims.metrics if not claims.metrics[k].get("structural")
           and not bt[(bt.metric == k) & bt.caught].shape[0]]
print(f"non-structural claims no row fails: {missing or 'none'}")

# %%
SAYS = {
    "mm_fit_KM_cmol_per_L": ["0.503055"],
    "mm_fit_rm_cmol_per_L_min": ["0.162583"],
    "mm_fit_k_per_min": ["0.323192"],
    "mm_fit_K_L_per_cmol": ["1.987853"],
    "mm_fit_rms_cmol_per_L": ["4.02e-3"],
    "mm_fit_maxres_over_print_resolution": ["13.9", "13.87"],
    "mm_fit_sse_penalty_K_pinned_at_2": ["0.108 %"],
    "mm_fit_KM_rel_dev_from_half": ["0.61 %"],
    "null_first_order_rms_ratio": ["9.71"],
    "null_zero_order_rms_ratio": ["13.29"],
    "null_second_order_rms_ratio": ["26.04"],
    "null_first_order_rms_cmol_per_L": ["3.90e-2"],
    "null_zero_order_rms_cmol_per_L": ["5.34e-2"],
    "null_second_order_rms_cmol_per_L": ["1.05e-1"],
    "rate_law_rho_degeneracy_max_rel": ["below 1e-14"],
    "rate_law_rho_degeneracy_random_draws_max_rel": ["below 1e-14"],
    "rate_law_break_k2_plus_1pct": ["9.8e-3"],
    "rho_spread_over_eps_trace_limit": ["0.242"],
    "rho_spread_induction_over_eps": ["0.977", "0.9773"],
    "rho_spread_induction_peak_min": ["2.0e-7", "2.011e-7"],
    "sampled_max_shortfall_uniform_scan": ["66 %", "66.3 %"],
    "induction_time_rho_lo_min": ["2.09e-8"],
    "induction_time_rho_hi_min": ["2.09e-4"],
    "eps_star_rho_visible_at_print_resolution": ["2.07e-3", "2.069e-3"],
    "eps_star_discrimination_with_compensation": ["5.80e-2", "5.798e-2"],
    "eps_star_rms_25pt_grid_no_compensation": ["3.18e-3", "3.185e-3"],
    "eps_star_compensation_penalty": ["18.2", "18.20"],
    "eps_star_statistic_and_window_factor": ["1.54", "1.539"],
    "eps_star_raw_ratio_unlike_statistics": ["28.0"],
    "eps_star_comp_res_doubled": ["8.757e-2"],
    "eps_star_comp_grid13": ["6.954e-2"],
    "eps_star_rms_rho_narrow": ["9.622e-3"],
    "eps_star_single_start_rel_move": ["below 1e-5"],
    "rho_profile_sse_span_fraction": ["8.6 %"],
    "full_model_F_statistic": ["1.90", "1.9009"],
    "full_model_F_crit_95": ["9.55", "9.5521"],
    "lineweaver_printed_axes_nonlinearity_minimax": ["31.7 %", "31.72 %", "31.7256 %"],
    "lineweaver_printed_axes_nonlinearity": ["50.3 %", "50.32 %"],
    "lineweaver_printed_axes_nonlinearity_uniform_measure": ["66.9 %", "66.88 %"],
    "lineweaver_correct_axes_nonlinearity": ["below 1e-13"],
    "lb_KM_bias_vs_nls": ["-0.4673", "-46.7 %"],
    "lb_rm_bias_vs_nls": ["-0.3406", "-34.1 %"],
    "lb_KM_bias_finite_difference_only": ["0.1402", "14.0 %"],
    "qssa_group_spread_CE0_over_cA0": ["94.9", "95.0"],
    "qssa_group_spread_CE0_over_cA0_plus_KM": ["1.91"],
    "qssa_error_over_eps_KM_0p05": ["0.9541"],
    "qssa_error_at_discrimination_eps_over_resolution": ["88.9"],
    "qssa_error_at_discrimination_eps_cmol_per_L": ["0.04443"],
    "qssa_error_at_discrimination_eps_total_over_resolution": ["38.8"],
    "rho_star_equilibrium_route": ["2.195e-3", "2.1951e-3"],
    "rho_star_equilibrium_analytic": ["2.1934e-3"],
    "rho_star_equilibrium_two_route_rel_diff": ["8.00e-4", "0.080 %"],
    "dsdKM_peak_s_cmol_per_L": ["0.228461"],
    "pymrm_be_time_order": ["between 0.9 and 1.1", "0.996"],
    "qssa_error_order_in_eps": ["between 0.9 and 1.1", "0.9988"],
    "pymrm_richardson_vs_radau": ["4.19e-7"],
    "pymrm_enzyme_balance_max": ["below 1e-13"],
    "pymrm_substrate_balance_max": ["below 1e-13"],
    "lsoda_vs_radau_reference": ["below 1e-9"],
    "fit_second_route_k_rel_diff": ["below 1e-5"],
    "fit_second_route_K_rel_diff": ["below 1e-6"],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)
print(f"prose literals: {sum(len(v) for v in SAYS.values())} checked against their claims "
      f"with claims.says")
claims.report();   # the semicolon keeps the notebook from echoing the payload

# %% [markdown]
# ## What pymrm adds
#
# **Honestly: for this page, not much, and the reason is worth stating.** A
# well-mixed batch of four species is 0-D. `newton` plus `NumJac((1, 4))` gives a
# compact, correct backward-Euler marcher and the conservation checks come out at
# round-off, but `solve_ivp` would have integrated the same system, and does, as
# the reference. What pymrm buys here is (i) a Jacobian assembled by the same
# mechanism the rest of the gallery uses, so the model drops into a spatial page
# unchanged, and (ii) the discipline of `NumJac((1, n_c))` rather than
# `NumJac((n_c,))`, which is the difference between a pointwise block and a dense
# matrix once the cell count stops being 1. Every Newton step checks `success`
# and raises; `pymrm.checks.find_roots` brackets every threshold on its whole
# range and refuses anything but one root; `pymrm.checks.observed_orders` gives
# the time-step orders from successive differences.
#
# **Point (i) is the real one.** `J4.7` is this rate law in a diffusing particle
# and `B1.1` is that same `S3` problem with a power law; the `source` method of
# `BatchReduced` is transplantable into either without change. The reason this
# page exists as `S1` and not as part of `J4.7` is that the question it settles -
# what data can decide between two derivations of the same constant - is a
# question about the *kinetics*, and putting it inside a pellet would only add a
# diffusional resistance that further hides the elementary constants.
#
# **What the page adds to the source**, as opposed to what pymrm adds:
#
# - the printed chain closed symbolically, ten identities, both routes;
# - the exact statement of what steady-state data cannot do - the rate law is
#   provably independent of $k_2/k_{-1}$, $\partial r/\partial\rho \equiv 0$ in
#   sympy, and identical to double precision over eight decades of it - which the
#   book gestures at in one sentence about stopped-flow experiments and does not
#   quantify;
# - a number for "the decomposition is the rate determining step":
#   $\rho \le 2.195\times10^{-3}$ for route 1 to predict the batch curve to the
#   source table's own resolution, arrived at two independent ways;
# - the group that organises the pseudo-steady-state error,
#   $\varepsilon = C_{\rm E}^0/(c_{\rm A0}+K_M)$, established by a $K_M$ sweep at
#   fixed $C_{\rm E}^0/c_{\rm A0}$ - Froment's printed comparison is with the
#   substrate alone and is the $K_M \ll c_{\rm A0}$ corner of it;
# - and the enzyme loading a discriminating experiment would need, with the
#   finding that at that loading the Michaelis-Menten form is already wrong by
#   88.9 times the resolution the discrimination is read at.
#
# **Two printed defects reported, neither repaired**, both settled from the book's
# own equations rather than from outside knowledge.

# %% [markdown]
# ## Reuse
#
# **If you have rate-versus-substrate data and want $K_M$:** fit
# eq. (1.5.1-15) - or, for a batch run, its integral - directly, and report
# $(r_m, K_M)$. Do **not** report $k_{-1}/k_1$ or $(k_{-1}+k_2)/k_1$ as though
# your data distinguished them. They do not, and this is exact, not approximate:
# the fibre of the observable map is two-dimensional and $k_2/k_{-1}$ runs along
# it. If you need the elementary constants, Froment p. 25 says what is required
# and this page prices it: transient experimentation, or a batch run at
# $\varepsilon = C_{\rm E}^0/(c_{\rm A0}+K_M)$ of order $10^{-1}$ - at which point
# you must fit the mass-action model, because the Michaelis-Menten form is no
# longer valid there.
#
# **If you are linearising:** the abscissa is $1/C_{\rm A}$, and even then
# consider not doing it. On the one dataset here, the reciprocal construction run
# exactly as both books instruct returns $K_M$ **47 % low** and $r_m$ **34 % low**
# against the direct fit of the same seven points; +14 points of that is the
# central-difference rate estimate and the rest is the reciprocal transform
# weighting the three smallest concentrations most heavily. Levenspiel p. 615 says
# the direct fit "is direct, is less prone to fiddling, and is more reliable"; on
# these data that is measurable and correct.
#
# **If you are checking whether the quasi-steady state is safe:** compute
# $\varepsilon = C_{\rm E}^0/(c_{\rm A0}+K_M)$, not $C_{\rm E}^0/c_{\rm A0}$. The
# maximum error in $C_{\rm A}(t)$ is $\approx 0.65$-$0.97\,\varepsilon$ (the
# coefficient falls with $\rho$), first order, over the four decades of
# $C_{\rm E}^0$ swept here (five values, $10^{-5}$ to $10^{-1}$). Comparing with
# the substrate alone mispredicts it by up to a factor 95 once $K_M$ is
# comparable to or larger than $c_{\rm A0}$.
#
# **If you want this rate law inside transport:** take `BatchReduced.source` to
# `J4.7` (immobilised enzyme particle, `S3`) or to `B1.1`'s pellet; the pymrm
# ingredients change from `newton` + `NumJac` alone to `construct_grad` +
# `construct_div` + `newton`, and the rate law is unchanged.
#
# **What NOT to reuse.** The seven-point table is an exercise dataset whose
# provenance the source book does not state. It is used here because the gallery
# is mostly tier-6 for want of any measurement at all, and because what the page
# concludes does not depend on it - the degeneracy, the $\varepsilon$ scaling and
# the two printed defects are all properties of the printed equations. Do not
# cite these numbers as a measured $K_M$ for any enzyme; no enzyme is named.
