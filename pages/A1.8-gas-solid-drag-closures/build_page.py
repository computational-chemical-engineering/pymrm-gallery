#!/usr/bin/env python3
"""Generate index.ipynb for page A1.8. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Two ways to close gas-solid drag: Syamlal-O'Brien against Ergun, and against the velocity-voidage law it stands in for"
description: "The MFIX Theory Guide names two kinds of data a drag law can be built from and picks one. This page builds both, measures where they agree and where they do not, and tests the closed form against the 1954 index it replaces and against 21 measured minimum-fluidisation velocities."
categories: [sec:A, struct:S3, tier:T1, data:tier2, phase:gas-solid]
date: 2026-08-02
---

# Two ways to close gas–solid drag

**Catalog ID:** `A1.8` · **Structures:** `S3` · **Tier:** T1

Every two-fluid simulation of a fluidised bed contains one scalar function that
decides almost everything the answer does: the drag coefficient $F_{gm}$ linking
the gas and the solids momentum equations. There is no first-principles form for
it, and the MFIX Theory Guide states the reason plainly — there are **two kinds
of experiment** you can build one from, and they do not overlap.

> Two types of experimental data can be used to develop fluid-solids drag
> formulas. One type, valid for high value of the solids volume fractions, is
> packed-bed pressure drop data expressed in the form of a correlation, such as
> the Ergun (1952) equation. Such a correlation must be supplemented with a drag
> correlation for low values of the solids volume fractions (Gidaspow 1986). The
> other type of data is available as correlations for the terminal velocity in
> fluidized or settling beds, expressed as a function of void fraction and
> Reynolds number (Richardson and Zaki 1954).
>
> — Syamlal, Rogers and O'Brien (1993), journal page 10

That paragraph is the page. It names three closures — an Ergun-based one, a
terminal-velocity-based one, and the blend rule that patches the first at low
solids fraction — and this page builds the two it prints in full, puts numbers on
where they agree and where they do not, and tests each against a measurement.

**Scope, stated before anything else.** The catalogue asks for
Gidaspow / Syamlal–O'Brien / Wen–Yu. The report on disk carries
**Syamlal–O'Brien complete** and **neither of the other two drag closures**: the
strings *"Wen and Yu"* and *"Wen & Yu"* do not occur anywhere in it, and it
prints **no Gidaspow drag law and no blend rule**. Gidaspow (1986) is cited
exactly once in §2.2.1, for the sentence quoted above.

*That last sentence is narrower than it looks, and an earlier draft of this page
got it wrong.* Gidaspow is **not** a marginal presence in this report — the name
runs through the introduction, the ill-posedness discussion, the solids-stress
section and the heat-transfer section, appears as first author on four separate
entries in the reference list, and the report **adopts a Ding & Gidaspow (1990)
expression as its own eq. (88)** for granular-energy transfer and a
Syamlal & Gidaspow (1985) model for solids conductivity. What is absent is
specifically a *drag* closure of that name, which is the only thing this page
needs to be true. Nothing about either missing closure is written here from
memory. What replaces them is the pair the report's own paragraph points at, both
on disk: **Ergun (1952)**, through the published page
[`A1.1`](../A1.1-ergun-pressure-drop/), and **Richardson & Zaki (1954)**, which
is the correlation Syamlal–O'Brien's closed form is a substitute for. The
[Reuse](#reuse) section says exactly what document would complete the
three-way comparison the catalogue names."""))

# ----------------------------------------------------------------- background
cells.append(md(r"""## Background

### Where the closure sits

The report's momentum equations, journal pages 7–9, are the standard two-fluid
pair. For the gas (eq. 7) and for the $m$-th solids phase (eq. 8),

$$
\frac{\partial}{\partial t}(\epsilon_g\rho_g\vec v_g) + \nabla\!\cdot(\epsilon_g\rho_g\vec v_g\vec v_g)
= \nabla\!\cdot\overline{S}_g + \epsilon_g\rho_g\vec g - \sum_m \vec I_{gm} + \vec f_g ,
$$

$$
\frac{\partial}{\partial t}(\epsilon_{sm}\rho_{sm}\vec v_{sm}) + \nabla\!\cdot(\epsilon_{sm}\rho_{sm}\vec v_{sm}\vec v_{sm})
= \nabla\!\cdot\overline{S}_{sm} + \epsilon_{sm}\rho_{sm}\vec g + \vec I_{gm} - \sum_{l\neq m} \vec I_{ml} ,
$$

coupled through the interaction force of eq. (9),

$$
\vec I_{gm} = -\epsilon_{sm}\nabla P_g \;-\; F_{gm}(\vec v_{sm}-\vec v_g)
\;+\; R_{0m}\left[\xi_{0m}\vec v_{sm} + \bar\xi_{0m}\vec v_g\right] ,
$$

whose three terms the report labels buoyancy, drag, and momentum transfer
accompanying mass transfer. **Everything on this page is about the single
scalar $F_{gm}$.**

Two consequences of that structure are used throughout and are worth writing
down now, because both are derived from the equations above rather than assumed.

**In a fixed bed** ($\vec v_{sm}=0$) the gas equation collapses to
$\;(-\mathrm{d}P/\mathrm{d}z)_{\text{friction}} = F_{gm}\,U/\epsilon_g^2$, with
$U$ the superficial velocity. That is what converts a pressure-drop correlation
into a drag coefficient, and it is how Ergun enters this page.

**In a uniformly fluidised or settling bed**, summing the two equations gives
$\mathrm{d}P/\mathrm{d}z = -(\epsilon_g\rho_g + \epsilon_s\rho_s)g$, and
substituting back into the solids equation gives

$$
F_{gm}\,|\vec v_g - \vec v_{sm}| \;=\; \epsilon_s\,\epsilon_g\,(\rho_s-\rho_g)\,g .
$$

Note the $\epsilon_g$ on the right: the buoyancy is on the *mixture* density, not
on the gas density. This relation holds for batch sedimentation and for
fluidisation alike, which is exactly the equivalence Richardson and Zaki
established experimentally, and it is what lets a terminal-velocity correlation
be converted into a drag law at all.

### Why there are two families and not one

A packed-bed pressure-drop correlation is measured where the particles touch. It
is at its best at $\epsilon_g \approx 0.4$ and it has nothing to say about an
isolated particle — as the report notes, it "must be supplemented" at low solids
fraction. A hindered-settling correlation is measured across the whole range from
a packed bed to infinite dilution, but it is fitted to a *velocity*, not to a
force, and the conversion needs a single-sphere drag curve to lean on.

The three closures the catalogue names are three answers:

- **Ergun-based**, extended into the dilute limit by a blend rule — the shape
  attributed to Gidaspow (1994), and the one this page cannot build.
- **Syamlal–O'Brien**, which converts a hindered-settling correlation, and so
  covers the whole voidage range with one continuous expression and no blend.
- **Wen–Yu**, a single-sphere drag curve with a voidage correction, which is the
  other half of the blend and is also not in the document on disk.

So this page can put the *first two families* against each other over the whole
$(\epsilon_g, Re)$ plane, and against measurement — which is the comparison the
report's own paragraph sets up."""))

# ----------------------------------------------------------- published model
cells.append(md(r"""## The published model

### Syamlal–O'Brien, exactly as printed

Section 2.2.1, journal pages 10–11, read on 400 dpi renders — which is the
embedded page images' own native resolution, since they are bilevel CCITT at
400 dpi. **The text layer of this report must not be used for any digit.** It
renders $0.06\,Re$ as `0.O6Re` with a capital letter O, $\epsilon_g^{2.65}$ as
`_g-2'6s`, "Ergun (1952)" as `Ergun (f952)`, and the equation numbers (12) and
(16) as `(121` and `(1''`.

$$
F_{gm} = \frac{3\,\epsilon_{sm}\,\epsilon_g\,\rho_g}{4\,V_{rm}^2\,d_{pm}}\;
C_{Ds}\!\left(\frac{Re_m}{V_{rm}}\right)\,\left|\vec v_{sm}-\vec v_g\right| \tag{11}
$$

$$
V_{rm} = 0.5\left(A - 0.06\,Re_m + \sqrt{(0.06\,Re_m)^2 + 0.12\,Re_m(2B-A) + A^2}\right) \tag{12}
$$

$$
A = \epsilon_g^{4.14} \tag{13}
\qquad\qquad
B = \begin{cases} 0.8\,\epsilon_g^{1.28} & \epsilon_g \le 0.85\\[2pt]
\epsilon_g^{2.65} & \epsilon_g > 0.85\end{cases} \tag{14}
$$

$$
Re_m = \frac{d_{pm}\,|\vec v_{sm}-\vec v_g|\,\rho_g}{\mu_g} \tag{15}
\qquad\qquad
C_{Ds}(Re) = \left(0.63 + \frac{4.8}{\sqrt{Re}}\right)^2 \tag{16}
$$

and the sentence under eq. (16), which is easy to miss and load-bearing: *"To use
this formula in equation (11), note that Re must be replaced with
$Re_m/V_{rm}$."*

**Attribution, as the report gives it**, and it is worth reading its own
reference list rather than the catalogue's. Eq. (12) with (13)–(14) is from
*Garside and Al-Dibouni (1977)*, `I&EC Proc. Des. Dev.` 16, 206–214; eq. (16) is
*Dalla Valle (1948)*, `Micromeritics`, Pitman. Neither is on disk and neither was
consulted; this page reads both from the 1993 report, which prints them in full
with attribution — the reprint route of `AGENTS.md`.

Eq. (11) is a different matter. The report attributes it to *Syamlal and O'Brien
(1987)*, and its reference list gives that as **"A Generalized Drag Correlation
for Multiparticle Systems," Unpublished report** — read on a 400 dpi render of
report page 44. So there is no origin paper to go and get: **this DOE report is
not a convenient source for the conversion, it is the citable published one.**
The catalogue's citation for the case, "Syamlal & O'Brien (1989)", is a
*different* item in the same list — "Computer Simulation of Bubbles in a
Fluidized Bed", AIChE Symposium Series No. 270, 85, 22–31 — which is a
bubble-simulation paper, not the drag derivation.

**There is no sphericity anywhere in this closure**, and the report says why. Its
paragraph "Four" on journal page 9 states that these formulations "deal with
uniform, smooth, spherical particles", that real systems do not, and that "there
are no well-accepted ways of treating such effects". That asymmetry against Ergun
matters later and is not a criticism: it is printed.

### Ergun, written as a drag coefficient

From [`A1.1`](../A1.1-ergun-pressure-drop/), Ergun's eq. (13c) with the two
constants that page recovers from his own figures. Substituting the fixed-bed
relation above, $(-\mathrm{d}P/\mathrm{d}z) = F_{gm}U/\epsilon_g^2$, and
$U = \epsilon_g v_{\text{slip}}$,

$$
F_{gm}^{\,\text{Ergun}}
= \frac{k_1\,\mu\,\epsilon_s^2}{\epsilon_g\,d_p^2}
\;+\; \frac{k_2\,\rho_g\,\epsilon_s\,|v_{\text{slip}}|}{d_p} .
$$

This is a *reading* of Ergun as a momentum-exchange coefficient, not something
Ergun wrote; the algebra is one line and is shown above. The constants are loaded
from `A1.1`'s dataset rather than typed.

### Richardson and Zaki, the correlation the closed form replaces

The report says $V_{rm}$ "can be calculated from the Richardson and Zaki (1954)
correlation only numerically; an explicit formula cannot be derived", and adopts
Garside and Al-Dibouni's closed form instead. That makes Richardson and Zaki the
natural reference for eqs. (12)–(14): **not an independent theory, but the thing
the substitution is a substitute for.**

Read from the verbatim Golden Jubilee reprint, `Trans IChemE` 75 (1997) S82–S100,
on 300 dpi renders (again the native resolution of the embedded bilevel images).
Their eq. (29), reprint page S93,

$$
\frac{V_c}{V_i} = \epsilon^{\,n}
$$

with $V_c$ "the falling velocity, relative to a fixed horizontal plane, of a
suspension of voidage $\epsilon$" and $V_i$ its value extrapolated to
$\epsilon = 1$; $V_i \to V_0$, the single-particle terminal velocity, as
$d/D \to 0$. Their Summary and Conclusions states the equivalence the whole thing
rests on:

> It has been shown that the falling velocity of a suspension relative to a fixed
> horizontal plane is equal to the upward velocity of liquid (based on the empty
> tube) required to maintain a suspension at the same concentration.

so $V_c$ is a **superficial** velocity in both processes. Two asymptotes are
printed on reprint page S94: $n = 4.65 + 19.5\,d/D$ for $Re < 0.2$ (their eq. 33)
and $n = 2.39$ for $Re$ above about 500 (their eq. 34); in between, their
Table VI tabulates $n_0$, the index extrapolated to $d/D = 0$, at eighteen
Reynolds numbers. That table is shipped with this page.

**The bridge between the two notations is one line, and it is not $V_{rm} =
\epsilon^n$.** Richardson and Zaki's $V_c$ is superficial, so
$V_c = \epsilon_g v_{\text{slip}}$; the force balance above shows (and the
Validation section proves numerically) that Syamlal–O'Brien's $V_{rm}$ is the
*slip* velocity ratio $v_{\text{slip}}/v_t$. Hence

$$
\epsilon_g\,V_{rm}(\epsilon_g, Re) \;\longleftrightarrow\; \epsilon_g^{\,n},
\qquad\text{i.e.}\qquad
n_{\text{implied}} = \frac{\ln\!\left(\epsilon_g V_{rm}\right)}{\ln \epsilon_g} .
$$

Getting that factor of $\epsilon_g$ wrong shifts every index on this page by
exactly 1 — an error the same size as the whole effect being measured — so it is
derived rather than asserted, and the derivation is checked in
[Validation](#validation)."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

The closures themselves take no parameters — they are dimensionless functions of
$\epsilon_g$ and $Re$ with printed constants. Everything numerical enters through
the two measured datasets.

| | |
|---|---|
| Ergun's constants $k_1$, $k_2$ | loaded from `A1.1`'s parameter file; the refit is loaded from its marker file and recomputed |
| gas properties for the $u_{mf}$ test | Geldart's own, as recorded in `A1.7`'s sidecar: $\mu = 1.8\times10^{-4}$ poise, $g = 981$ cm s$^{-2}$ |
| gas density for that test | $\rho_f = 1.2\times10^{-3}$ g cm$^{-3}$, air at room conditions. **Geldart prints no gas density**; this is an assumption, and the break table measures what it is worth |
| sphericity | $\phi_s = 1$ on the eight Diakon cuts, read off Geldart's phrase "a plastic moulding powder having spherical particles". See the warning below |
| voidage at minimum fluidisation | $\epsilon_0 = 1 - (1-\epsilon_{MB})H_{MB}/H_0$, the quantity `A1.7` derives from Geldart's two reported columns. Both that and $\epsilon_{MB}$ itself are reported on the page |

**Assumptions carried from the sources.** Steady, isothermal, single solids
phase, no mass transfer, no solids stress, uniform spheres of one size, and a bed
deep enough that the distributor and the free surface do not matter. The drag
comparison assumes the two closures are being asked the same question — a drag
force per unit volume at a given $(\epsilon_g, v_{\text{slip}})$ — which is
exactly what the two-fluid equations ask of them, and is *not* the same as
comparing a pressure drop with a settling velocity.

**One asymmetry that cannot be assumed away.** The Ergun form written in the
surface/volume diameter is exactly sphericity-free (`A1.6` demonstrates this),
while Syamlal–O'Brien is a sphere closure with no shape factor at all. So on any
non-spherical powder the two are not being given the same information, and this
page therefore takes its headline on Geldart's spherical cuts only."""))

# ------------------------------------------------------------ environment
cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code('''import sys, urllib.request
from pathlib import Path

if not any("shared" in p for p in sys.path):
    local = Path.cwd()
    for _ in range(6):
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
from matplotlib import colormaps
from IPython.display import display, Markdown
from scipy.optimize import brentq
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE     = "A1.8-gas-solid-drag-closures"
PAGE_A11 = "A1.1-ergun-pressure-drop"        # Ergun's constants and his 244 markers
PAGE_A17 = "A1.7-geldart-classification"     # the measured minimum-fluidisation velocities

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

# Okabe-Ito, assigned in a fixed order and never cycled:
#   0 Ergun (packed-bed route)  1 Syamlal-O'Brien  2 Richardson-Zaki  3 measured
COL = {"ergun": "#0072B2", "so": "#D55E00", "rz": "#009E73", "meas": "#E69F00"}'''))

# --------------------------------------------------------------------- data
cells.append(md(r"""## The data

Four datasets. Two are shipped with this page, two are borrowed — and borrowing
means reading the page that owns them, so what each of those pages has already
established about these exact rows is listed below and reconciled against this
page's own numbers where they overlap.

### 1. The eleven printed constants of eqs. (11)–(16) — shipped here

`mfix-1993-syamlal-obrien-constants.csv`. Not data in any evidential sense: they
are the coefficients of a correlation, and reproducing the correlation from them
would be reproduction, not validation. They are a dataset so that every constant
on this page has a single provenance record saying which equation it came from
and that it was read off an image rather than a text layer.

### 2. Richardson and Zaki's Table VI — shipped here

`richardson-zaki-1954-table6.csv`, 18 rows: the velocity–voidage index $n_0$ at
Reynolds numbers from 0.39 to 489, extrapolated to zero wall effect, plus the
paper's own two $\log_{10}$ columns.

**This is tier 2 but it is one step from a raw measurement, and the page says so
wherever it is used.** Each $n_0$ is a slope of $\log V_c$ against $\log\epsilon$
fitted by the authors to a run of measurements, and then extrapolated along a
*second* fitted line to $d/D = 0$. Two further limits: the runs are liquid–solid
while the report is a gas–solid code (nothing in eq. 29 is phase-specific, but
the data are not gas-fluidisation data); and the table states no voidage range,
so the window over which an index is compared is a choice made here, with its
sensitivity printed.

**And the two correlations are not independent witnesses.** Garside and
Al-Dibouni fitted their velocity–voidage relation to the same literature
Richardson and Zaki drew on. Comparing them measures **what the substitution the
report makes costs**, which is a real and falsifiable thing, and not "two
independent routes agreeing".

### 3. Ergun's constants and his 244 markers — borrowed from `A1.1`

What `A1.1` establishes about these exact rows, and whether it bears here:

- **The printed constants are $k_1$, $k_2$ in `ergun-1952-parameters.csv`.** They
  are printed below from that file and never typed into this page's prose.
  *Bears directly: they set the Ergun curve everywhere.*
- **Refitting the two constants to the 244 recovered markers gives a different
  pair**, with a stated $1/f_v$ weighting; unweighted the fit is dominated by a
  handful of high-Reynolds points and returns something else again. That refit is
  **recomputed here from the marker file** with the same weighting, and printed
  beside `A1.1`'s published value. *Bears: it bounds how much of any Ergun /
  Syamlal–O'Brien gap could be Ergun's own constants.*
- **The 244 markers are not Ergun's 640** — the dataset carries incomplete recall
  as a stated limitation. *Does not bear: nothing here counts markers, and
  scattered misses do not move a fit.*
- **`A1.1` measures Ergun's own equation against those markers at a small bias
  and a single-digit mean absolute deviation.** *Bears as context: that is the
  scatter of the correlation on its own data, and it is far smaller than the
  disagreements this page reports.*
- **Eisfeld and Schnitzlein's refit differs from Ergun mostly in the turbulent
  constant, not the famous 150.** *Bears as context and is quoted from `A1.1`'s
  own dataset in the closing discussion.*

### 4. Geldart's Table 1 — 21 measured minimum-fluidisation velocities, borrowed from `A1.7`

What `A1.7` and `A1.6` establish about these exact rows:

- **`U_0` and `U_MB` are measurements; $\epsilon_{MB}$ and $H_{MB}/H_0$ are
  *reported*, method unstated.** `A1.7`'s sidecar forbids calling the latter two
  measurements and this page does not. *Bears directly: the voidage this page
  needs comes from those two columns.*
- **`A1.7` derives $\epsilon_0 = 1-(1-\epsilon_{MB})H_{MB}/H_0$ from a solids
  balance and reports the range it takes over all 22 rows.** This page uses that
  derived quantity, recomputes it, and prints its own range beside `A1.7`'s
  reported one. *Bears directly.*
- **`A1.6` reports the unapproximated Ergun balance as strongly biased at the
  voidage Geldart himself reports on the eight spherical cuts, and calls it
  "worse than the correlation it is being used to judge".** That number is
  **recomputed independently here** and printed beside `A1.6`'s. *Bears directly
  — it is the reason this page does not present Ergun as a reference standard.*
- **`A1.6` reports Wen and Yu's eq. (1) on those same eight rows.** Also
  recomputed here, as a yardstick. *Bears as context.*
- **The one missing `U_0` is an em dash in the printed table**; 21 of 22 rows
  carry one. *Bears: the row counts on this page.*
- **All rows are group A or borderline A/B powders**, and $(N_{Re})_{mf} < 1$ on
  most of them. *Bears heavily: the inertial half of every closure here is
  essentially untested by this dataset, and the page says so rather than
  implying otherwise.*"""))

cells.append(code('''# ---- shipped with this page ------------------------------------------------
so_const = load_data("mfix-1993-syamlal-obrien-constants.csv", page=PAGE)
so_meta  = load_meta("mfix-1993-syamlal-obrien-constants.csv", page=PAGE)
rz       = load_data("richardson-zaki-1954-table6.csv", page=PAGE)
rz_meta  = load_meta("richardson-zaki-1954-table6.csv", page=PAGE)

# ---- borrowed --------------------------------------------------------------
erg_par  = load_data("ergun-1952-parameters.csv", page=PAGE_A11)
erg_mk   = load_data("ergun-1952-fig7-markers.csv", page=PAGE_A11)
erg_meta = load_meta("ergun-1952-parameters.csv", page=PAGE_A11)
gd       = load_data("geldart_1973_table1.csv", page=PAGE_A17)
gd_meta  = load_meta("geldart_1973_table1.csv", page=PAGE_A17)

print("Syamlal-O'Brien: the eleven printed constants of eqs. (11)-(16)")
print(so_const.to_string(index=False))
print(f"\\n{cite_data(so_meta)}")

print("\\n\\nRichardson & Zaki (1954) Table VI - the index the closed form stands in for")
print(rz.to_string(index=False))
print(f"\\n{cite_data(rz_meta)}")

# every constant is taken from the CSV, never retyped
C = {r.symbol: float(r.value) for r in so_const.itertuples()}
E = {r.quantity: float(r.value) for r in erg_par.itertuples()}
K1_PRINT, K2_PRINT = E["k1_viscous_constant"], E["k2_kinetic_constant"]
print(f"\\n\\nErgun's printed constants, read out of A1.1's parameter file: "
      f"k1 = {K1_PRINT:g}, k2 = {K2_PRINT:g}")'''))

cells.append(code('''# ---- A1.1's refit, recomputed here from A1.1's own marker file --------------
# A1.1 fits f_v = k1 + k2 x with EQUAL RELATIVE weight (1/f_v); it prints both
# that and the unweighted alternative, because the unweighted one is dominated
# by the top decade. The weighted one is what it quotes, so it is what is
# reproduced here.
_x, _fv = erg_mk.Re_over_one_minus_eps.to_numpy(), erg_mk.f_v.to_numpy()
_w = 1.0 / _fv
K1_REFIT, K2_REFIT = np.linalg.lstsq(
    np.vstack([np.ones_like(_x), _x]).T * _w[:, None], _fv * _w, rcond=None)[0]
display(Markdown(
    f"**Reconciliation with `A1.1`.** Refitting Ergun's two constants to the "
    f"**{len(erg_mk)}** markers in its dataset, with its stated $1/f_v$ weighting, gives "
    f"$k_1 = {K1_REFIT:.5f}$ and $k_2 = {K2_REFIT:.7f}$ against the printed "
    f"$k_1 = {K1_PRINT:g}$ and $k_2 = {K2_PRINT:g}$ read out of the same page's parameter "
    f"file, i.e. **{100*(K1_REFIT/K1_PRINT-1):+.2f} %** and "
    f"**{100*(K2_REFIT/K2_PRINT-1):+.2f} %**. `A1.1` publishes 151.9 and 1.697 for these; "
    f"the two agree to {100*abs(K1_REFIT/151.93654-1):.1e} % and "
    f"{100*abs(K2_REFIT/1.6967056-1):.1e} %, which confirms that this page is reading "
    f"that page's dataset the way that page did."))

# ---- Geldart, with A1.7's derived voidage ----------------------------------
gd["eps_0"] = 1.0 - (1.0 - gd.eps_MB) * gd.H_MB_over_H_0
meas = gd.dropna(subset=["U_0"]).copy()
dia  = meas[meas.powder == "Diakon"].copy()
display(Markdown(
    f"**Reconciliation with `A1.7`.** Geldart's table has **{len(gd)}** rows, "
    f"**{len(meas)}** with a measured $U_0$, of which **{len(dia)}** are the spherical "
    f"Diakon cuts. Recomputing `A1.7`'s solids-balance voidage "
    f"$\\\\epsilon_0 = 1-(1-\\\\epsilon_{{MB}})H_{{MB}}/H_0$ over all {len(gd)} rows gives "
    f"**{gd.eps_0.min():.3f} to {gd.eps_0.max():.3f}**; `A1.7` reports 0.441–0.570 for the "
    f"same quantity. On Diakon alone it is **{dia.eps_0.min():.3f} to {dia.eps_0.max():.3f}**, "
    f"and on the two coarsest cuts, where Geldart prints $H_{{MB}}/H_0 = 1.000$ and "
    f"$U_{{MB}} = U_0$, it equals the reported $\\\\epsilon_{{MB}}$ itself with no inference "
    f"at all: **{dia.eps_0.iloc[-1]:.3f}**."))
print(f"\\n{cite_data(gd_meta)}")'''))

# ------------------------------------------------------ pymrm implementation
cells.append(md(r"""## PyMRM implementation

**No pymrm operator appears on this page, and inventing one would obscure the
comparison rather than sharpen it.** Both closures are algebraic functions of
$(\epsilon_g, Re)$; the two comparisons against measurement are roots of scalar
equations; the limits are closed forms. There is no field to discretise. That is
the same call [`A1.1`](../A1.1-ergun-pressure-drop/) and
[`A1.6`](../A1.6-wen-yu-minimum-fluidisation/) make, and it is stated here where
a reader would look for a solver.

Where the closures do meet pymrm is one level up: $F_{gm}$ is the coefficient
that couples two momentum equations, so it is the term a two-fluid `S7`/`S10`
model carries. The [Reuse](#reuse) section says which pages that is, and the
functions below are written to be lifted straight into one — vectorised over
$(\epsilon_g, v_{\text{slip}})$ arrays, taking every constant from a dictionary
so that a defect can be injected without editing the model."""))

cells.append(code('''# ---------------------------------------------------------------------------
# The closures. Every constant comes from the CSV dictionary `C`, so the break
# tables below can perturb one and nothing else.
# ---------------------------------------------------------------------------
P0 = dict(pref=C["eq11_prefactor"], half=C["Vrm_half"], lin=C["Vrm_linear"],
          cross=C["Vrm_cross"], aexp=C["A_exponent"], bco=C["B_coefficient"],
          bde=C["B_exponent_dense"], bdi=C["B_exponent_dilute"],
          bsw=C["B_switch_voidage"], cd0=C["CD_offset"], cd1=C["CD_slope"])


def C_Ds(Re, p=P0):
    """Dalla Valle single-sphere drag function, eq. (16)."""
    return (p["cd0"] + p["cd1"] / np.sqrt(Re)) ** 2


def A_of(eps_g, p=P0):
    """Eq. (13)."""
    return eps_g ** p["aexp"]


def B_of(eps_g, p=P0):
    """Eq. (14) - two branches, switching at eps_g = 0.85."""
    return np.where(eps_g <= p["bsw"], p["bco"] * eps_g ** p["bde"], eps_g ** p["bdi"])


def V_rm(eps_g, Re_m, p=P0):
    """Garside & Al-Dibouni closed form, eq. (12). Re_m is built on the ACTUAL
    slip velocity, eq. (15) - not on the terminal velocity."""
    A, B = A_of(eps_g, p), B_of(eps_g, p)
    return p["half"] * (A - p["lin"] * Re_m
                        + np.sqrt((p["lin"] * Re_m) ** 2
                                  + p["cross"] * Re_m * (2 * B - A) + A * A))


def F_syamlal_obrien(eps_g, v_slip, d_p, rho_g, mu_g, p=P0):
    """Eq. (11), with Re replaced by Re_m/V_rm in C_Ds as the report instructs."""
    Re_m = d_p * v_slip * rho_g / mu_g                    # eq. (15)
    Vr = V_rm(eps_g, Re_m, p)
    return p["pref"] * (1 - eps_g) * eps_g * rho_g / (Vr ** 2 * d_p) \\
        * C_Ds(Re_m / Vr, p) * v_slip


def F_ergun(eps_g, v_slip, d_p, rho_g, mu_g, k1=None, k2=None):
    """Ergun's eq. (13c) read as a momentum-exchange coefficient. See the
    derivation in 'The published model': (-dP/dz) = F U / eps_g^2, U = eps_g v."""
    k1 = K1_PRINT if k1 is None else k1
    k2 = K2_PRINT if k2 is None else k2
    eps_s = 1 - eps_g
    return k1 * mu_g * eps_s ** 2 / (eps_g * d_p ** 2) + k2 * rho_g * eps_s * v_slip / d_p


def v_terminal(d_p, drho, rho_g, mu_g, g, p=P0):
    """Single-sphere terminal velocity on the SAME drag curve, eq. (16)."""
    f = lambda v: 0.75 * rho_g / d_p * C_Ds(d_p * v * rho_g / mu_g, p) * v * v - drho * g
    return brentq(f, 1e-14, 1e6, xtol=1e-18, rtol=1e-15)


def V_rm_at_balance(eps_g, Re_t, p=P0):
    """V_rm at the force balance. Implicit, because eq. (12) reads Re_m, which
    is built on the slip velocity, and the slip velocity is V_rm * v_t."""
    return brentq(lambda x: x - V_rm(eps_g, Re_t * x, p), 1e-16, 1 + 1e-9,
                  xtol=1e-16, rtol=1e-15)


print("closures defined; constants taken from the CSV:")
for k, v in P0.items():
    print(f"   {k:6s} = {v:g}")'''))

# ------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. What each closure reduces to

Both closures have closed asymptotes, and writing them out is what makes the
comparison legible. For Syamlal–O'Brien:

- **Viscous, $Re_m \to 0$.** $C_{Ds} \to 4.8^2/(Re_m/V_{rm})$ and
  $V_{rm} \to A = \epsilon_g^{4.14}$, so
  $F \to \tfrac{3}{4}\!\cdot\!4.8^2\,\epsilon_s\,\epsilon_g^{-3.14}\,\mu/d_p^2$.
- **Inertial, $Re_m \to \infty$.** $C_{Ds} \to 0.63^2$ and $V_{rm} \to B$, so
  $F \to \tfrac34 \cdot 0.63^2\,\epsilon_s\epsilon_g\rho_g v_{\text{slip}}/(B^2 d_p)$.
  The $V_{rm}\to B$ limit is not obvious from eq. (12) and is checked numerically
  below.
- **Dilute, $\epsilon_g \to 1$.** $A = B = 1$, so $V_{rm} = 1$ *exactly, at every
  Reynolds number* — the radicand becomes a perfect square — and eq. (11)
  collapses to the isolated-sphere drag times the number density. That is the
  correct dilute limit by construction, and it is the property the Ergun route
  lacks.

For Ergun the two limits are its own two terms. The ratio of the two closures in
each limit is therefore a closed function of voidage alone, printed below."""))

cells.append(code('''# --- the two limits, as closed forms, then checked against the full expression
def ratio_viscous(eps_g, k1=None, p=P0):
    k1 = K1_PRINT if k1 is None else k1
    return (p["pref"] * p["cd1"] ** 2) / k1 * eps_g ** (1 - p["aexp"]) * eps_g / (1 - eps_g)


def ratio_inertial(eps_g, k2=None, p=P0):
    k2 = K2_PRINT if k2 is None else k2
    return (p["pref"] * p["cd0"] ** 2) / k2 * eps_g / B_of(eps_g, p) ** 2


DP, RHOG, MUG = 300e-6, 1.2, 1.8e-5          # air, a 300 um particle: sets Re only
eps_probe = np.array([0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.99])
lim = pd.DataFrame({"eps_g": eps_probe})
lim["viscous, closed form"] = ratio_viscous(eps_probe)
lim["viscous, numeric"] = [F_syamlal_obrien(e, 1e-9, DP, RHOG, MUG)
                           / F_ergun(e, 1e-9, DP, RHOG, MUG) for e in eps_probe]
lim["inertial, closed form"] = ratio_inertial(eps_probe)
lim["inertial, numeric"] = [F_syamlal_obrien(e, 1e7, DP, RHOG, MUG)
                            / F_ergun(e, 1e7, DP, RHOG, MUG) for e in eps_probe]
print("F(Syamlal-O'Brien) / F(Ergun) in the two limits")
print(lim.round(4).to_string(index=False))
LIM_ERR = max(np.abs(lim["viscous, closed form"] / lim["viscous, numeric"] - 1).max(),
              np.abs(lim["inertial, closed form"] / lim["inertial, numeric"] - 1).max())

# single-sphere Stokes limit of the Dalla Valle curve
STOKES_RATIO = P0["cd1"] ** 2 / 24.0
display(Markdown(
    f"**The single-sphere curve is not Stokes, and by a known amount.** Eq. (16) gives "
    f"$C_{{Ds}}Re \\\\to {P0['cd1']**2:g}$ as $Re \\\\to 0$, against Stokes' 24 — "
    f"**{100*(STOKES_RATIO-1):+.2f} %**, so every terminal velocity on this page is "
    f"{100*(1/STOKES_RATIO-1):+.2f} % above the Stokes value. In the other limit it gives "
    f"$C_{{Ds}} \\\\to {P0['cd0']**2:.4f}$. Both are consequences of two printed constants "
    f"and neither is an error; they are stated because the {100*(1/STOKES_RATIO-1):.1f} % "
    f"propagates into every number in Check 4.\\n\\n"
    f"The closed forms above reproduce the full expression evaluated at "
    f"$Re = 10^{{-9}}$ and $10^{{7}}$ to **{100*LIM_ERR:.3f} %** at worst, which is the "
    f"residual of not having taken the limit exactly, not agreement between two routes."))

# where the limits cross
x_v = [brentq(lambda e: ratio_viscous(e) - 1, a, b) for a, b in [(0.25, 0.6), (0.6, 0.995)]]
x_i = brentq(lambda e: ratio_inertial(e) - 1, 0.25, P0["bsw"] - 1e-9)
from scipy.optimize import minimize_scalar
mv = minimize_scalar(ratio_viscous, bounds=(0.3, 0.99), method="bounded")
display(Markdown(
    f"In the **viscous** limit the two closures cross at $\\\\epsilon_g = {x_v[0]:.4f}$ and "
    f"again at ${x_v[1]:.4f}$, with Syamlal–O'Brien below Ergun in between (worst "
    f"{ratio_viscous(mv.x):.3f} at $\\\\epsilon_g = {mv.x:.3f}$) and above it outside. "
    f"In the **inertial** limit they cross once, at $\\\\epsilon_g = {x_i:.4f}$, and after "
    f"that Syamlal–O'Brien falls away monotonically — it is "
    f"{ratio_inertial(0.99):.3f} of Ergun at $\\\\epsilon_g = 0.99$.")) '''))

cells.append(md(r"""### 2. Where they agree, where they diverge, and by how much

The limits are the easy part. The question a simulation actually asks is what the
two closures do *in between*, and the answer is the least comfortable number on
this page."""))

cells.append(code('''# Both closures reduce to F d_p^2 / mu as a function of (eps_g, Re_m) alone.
def F_so_star(eps_g, Re_m, p=P0):
    Vr = V_rm(eps_g, Re_m, p)
    return p["pref"] * (1 - eps_g) * eps_g / Vr ** 2 * C_Ds(Re_m / Vr, p) * Re_m


def F_er_star(eps_g, Re_m, k1=None, k2=None):
    k1 = K1_PRINT if k1 is None else k1
    k2 = K2_PRINT if k2 is None else k2
    return k1 * (1 - eps_g) ** 2 / eps_g + k2 * (1 - eps_g) * Re_m


Re_grid = np.logspace(-4, 6, 400)
eps_fam = [0.40, 0.44, 0.50, 0.60, 0.70, 0.80]
shades = colormaps["Blues"](np.linspace(0.38, 0.95, len(eps_fam)))  # sequential: voidage IS a magnitude

fig, ax = plt.subplots(figsize=(7.2, 4.4))
for e, c in zip(eps_fam, shades):
    r = F_so_star(e, Re_grid) / F_er_star(e, Re_grid)
    ax.plot(Re_grid, r, color=c, lw=2, label=f"$\\\\epsilon_g$ = {e:.2f}")
    j = np.argmin(r)
    ax.plot(Re_grid[j], r[j], "o", ms=5, color=c, mec="white", mew=1.2)
ax.axhline(1.0, color="0.35", lw=1, ls="--")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("$Re_m$  (on the slip velocity)")
ax.set_ylabel("$F$(Syamlal–O'Brien) / $F$(Ergun)")
ax.set_title("Two closures for one coefficient; markers are the worst point on each curve")
ax.legend(frameon=False, ncol=2, fontsize=9)
plt.tight_layout(); plt.show()

rows = []
for e in eps_fam:
    r = F_so_star(e, Re_grid) / F_er_star(e, Re_grid)
    j = int(np.argmin(r))
    rows.append((e, ratio_viscous(e), ratio_inertial(e), r[j], Re_grid[j]))
worst = pd.DataFrame(rows, columns=["eps_g", "viscous limit", "inertial limit",
                                    "worst ratio", "at Re_m"])
worst["interior minimum"] = worst["at Re_m"] < 0.9 * Re_grid[-1]
print(worst.round(4).to_string(index=False))
W = worst[worst.eps_g == 0.44].iloc[0]
W50 = worst[worst.eps_g == 0.50].iloc[0]
V_AT_DIP = W["at Re_m"] * 1.8e-5 / (300e-6 * 1.2)      # slip velocity of a 300 um particle in air
# Ergun plots against Re/(1-eps_g), so put the dip on HIS abscissa and compare it
# with the span of the 244 markers this page loads, refits and quotes.
ERG_X = erg_mk.Re_over_one_minus_eps.to_numpy()
X_DIP44 = float(W["at Re_m"] / (1 - W["eps_g"]))
X_DIP50 = float(W50["at Re_m"] / (1 - W50["eps_g"]))
display(Markdown(
    f"**The headline.** At a packed-bed voidage the two closures agree in *both* asymptotic "
    f"limits to better than a factor {max(worst.loc[worst.eps_g<=0.44,'viscous limit'].max(), 1/worst.loc[worst.eps_g<=0.44,'inertial limit'].min()):.2f} — "
    f"and disagree by a factor **{1/W['worst ratio']:.2f}** in between. At "
    f"$\\\\epsilon_g = {W['eps_g']:.2f}$ the ratio bottoms out at **{W['worst ratio']:.3f}** at "
    f"$Re_m = {W['at Re_m']:.2f}$ — a 300 µm particle in air slipping at "
    f"{V_AT_DIP:.2f} m s$^{{-1}}$, which is not an exotic corner but the middle of a "
    f"bubbling bed. Agreement in both limits is therefore **not** evidence that two "
    f"closures agree, and a page that reported only the limits would have said the "
    f"opposite of the truth.\\n\\n"
    f"The dip is a genuine interior minimum only while the bed is dense: the last column "
    f"marks where the worst point is inside the swept range rather than at the inertial "
    f"asymptote, and it stops being interior above "
    f"$\\\\epsilon_g \\\\approx {worst.loc[worst['interior minimum'],'eps_g'].max():.2f}$. "
    f"Above that the closures simply diverge monotonically, to the inertial-limit ratios "
    f"in column three.\\n\\n"
    f"**And the qualification that belongs beside the headline, not in a footnote.** Ergun "
    f"plots against $Re/(1-\\\\epsilon_g)$, and on that abscissa the "
    f"$\\\\epsilon_g = {W['eps_g']:.2f}$ dip sits at **{X_DIP44:.2f}**. The "
    f"**{len(ERG_X)}** markers recovered from his Figure 7 — the ones this page loads, "
    f"refits and quotes — span **{ERG_X.min():.2f} to {ERG_X.max():.0f}**. So *every one of "
    f"Ergun's own recovered data points lies above the point where this page reports a "
    f"factor {1/W['worst ratio']:.2f} disagreement*, though only just: the gap is "
    f"{100*(ERG_X.min()/X_DIP44-1):.1f} %. At $\\\\epsilon_g = {W50['eps_g']:.2f}$ the dip "
    f"moves to {X_DIP50:.2f} and is comfortably inside the data. The obvious objection — "
    f"*of course they disagree, you have left the range Ergun measured* — is therefore "
    f"live at 0.44 and answered at 0.50, and the honest reading is that the disagreement "
    f"is real where Ergun has data and is an extrapolation of his correlation, by a "
    f"whisker, at the worst point. **No validity range for eq. (12) is available on "
    f"disk either**: the report gives none and Garside & Al-Dibouni (1977) is not here.")) '''))

cells.append(md(r"""### 3. The dilute limit — why the report says the Ergun route "must be supplemented"

That sentence is usually paraphrased as "Ergun is not valid in dilute flow",
which is true but not quantitative. It has a closed form. As $\epsilon_s \to 0$
the Ergun viscous term vanishes like $\epsilon_s^2$ while the drag on an isolated
sphere does not vanish at all, so the ratio of the Ergun drag *per particle* to
the isolated-sphere drag tends to $k_2/(\tfrac34 C_{Ds}(Re))$ — which goes to
zero in the Stokes regime and to a finite number above it."""))

cells.append(code('''Re_d = np.array([1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e6])
dil = pd.DataFrame({"Re": Re_d})
dil["Ergun / isolated sphere"] = K2_PRINT / (P0["pref"] * C_Ds(Re_d))
# and the same thing measured on the full expressions at a small but finite eps_s.
# NOT the same number: the closed form drops Ergun's viscous term, which survives
# at finite eps_s, so the last column also carries the share that term still has.
eps_dilute = 1 - 1e-4
dil["full expression, eps_s = 1e-4"] = [F_er_star(eps_dilute, r) / F_so_star(eps_dilute, r)
                                        for r in Re_d]
_visc = K1_PRINT * (1 - eps_dilute) ** 2 / eps_dilute
dil["Ergun viscous term, % of its total"] = 100 * _visc / (_visc + K2_PRINT * (1 - eps_dilute) * Re_d)
print("How badly the packed-bed route does where it was never measured")
print(dil.round(5).to_string(index=False))
DIL_GAP = float(dil["full expression, eps_s = 1e-4"].iloc[0] / dil["Ergun / isolated sphere"].iloc[0])
display(Markdown(
    f"**The two columns are not the same quantity at the dense-viscous corner, and the "
    f"third column says why.** The closed form is the $\\\\epsilon_s\\\\to0$ limit, in which "
    f"Ergun's viscous term has already vanished like $\\\\epsilon_s^2$; at a finite "
    f"$\\\\epsilon_s = 10^{{-4}}$ it has not, and at $Re = 10^{{-3}}$ it is still "
    f"**{dil['Ergun viscous term, % of its total'].iloc[0]:.0f} %** of Ergun's total, which is "
    f"why the second column sits **{DIL_GAP:.1f}×** above the first there. The two converge "
    f"as the viscous share dies: {dil['full expression, eps_s = 1e-4'].iloc[3]/dil['Ergun / isolated sphere'].iloc[3]:.2f}× "
    f"by $Re = 1$. The prose below quotes the **closed form**, i.e. the true dilute limit.\\n\\n"
    f"At $Re = 1$ the Ergun form gives **{100*K2_PRINT/(P0['pref']*C_Ds(1.0)):.1f} %** of the "
    f"drag on an isolated sphere; at $Re = 10^{{-2}}$ it gives "
    f"**{100*K2_PRINT/(P0['pref']*C_Ds(1e-2)):.2f} %**, and the shortfall is unbounded as "
    f"$Re \\\\to 0$. Above $Re \\\\approx {brentq(lambda r: K2_PRINT/(P0['pref']*C_Ds(r))-1, 1, 1e4):.0f}$ "
    f"it overshoots instead, reaching **{K2_PRINT/(P0['pref']*C_Ds(1e6)):.2f}×** at high "
    f"Reynolds number. **That two-sided failure is the whole reason a blend rule exists**, "
    f"and it is the piece of the catalogue's three-way comparison this page cannot supply: "
    f"the blend attributed to Gidaspow (1994) is in a monograph that is not on disk. The "
    f"Syamlal–O'Brien column needs no blend — by construction $V_{{rm}}(1) = 1$ and eq. (11) "
    f"becomes the isolated-sphere drag exactly, which is Check 1 below.")) '''))

cells.append(md(r"""### 4. The switch at $\epsilon_g = 0.85$

Eq. (14) is a two-branch definition, so the obvious question is how big the jump
is. It is worth asking precisely because the blend rule this page cannot build —
the one attributed to Gidaspow — is *known for* switching discontinuously, and it
would be easy to assume the same of any branched closure."""))

cells.append(code('''sw = P0["bsw"]
B_lo = P0["bco"] * sw ** P0["bde"]
B_hi = sw ** P0["bdi"]
jumps = [(r, 100 * (V_rm(sw + 1e-12, r) / V_rm(sw - 1e-12, r) - 1))
         for r in [1e-6, 1e-2, 1.0, 10.0, 100.0, 1e4, 1e8]]
print(f"B just below the switch : {B_lo:.8f}   (= {P0['bco']:g} * {sw:g}^{P0['bde']:g})")
print(f"B just above the switch : {B_hi:.8f}   (= {sw:g}^{P0['bdi']:g})")
print(f"relative jump in B      : {100*(B_hi/B_lo-1):+.5f} %\\n")
print("relative jump in V_rm across the switch, by Reynolds number:")
for r, j in jumps:
    print(f"   Re_m = {r:9.0e}   {j:+.6f} %")
JUMP_B = 100 * (B_hi / B_lo - 1)
JUMP_V = max(abs(j) for _, j in jumps)
display(Markdown(
    f"**It is a branch, not a discontinuity.** The two expressions meet at "
    f"$\\\\epsilon_g = {sw:g}$ to **{JUMP_B:+.3f} %** in $B$, and the largest jump this "
    f"produces in $V_{{rm}}$ anywhere in Reynolds number is **{JUMP_V:.3f} %** — smaller "
    f"than the third significant figure of either exponent, so eq. (14) is a change of "
    f"functional form and not a step. **This page therefore makes no claim about "
    f"any discontinuous drag closure**: the one the catalogue names switches at "
    f"$\\\\epsilon_g = 0.8$, is in a monograph that is not on disk, and the size of its jump "
    f"cannot be computed here.")) '''))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

Four checks and a transcription round trip, in the order of
`docs/agent-brief.md` — the two identities first because they are what licenses
the transcription, then the two that can actually fail.

**Read the first two as what they are.** They are algebraic identities of the
printed equations. They cannot fail for a correct transcription, and their value
is entirely in *which mis-readings they do and do not catch* — so each is
published with the break table that measures that, and each says explicitly what
it is blind to. Neither is evidence that the closure describes anything."""))

cells.append(md(r"""### Check 1 — $V_{rm}(\epsilon_g = 1) = 1$ at every Reynolds number

At $\epsilon_g = 1$ eqs. (13) and (14) both give 1, the radicand of eq. (12)
becomes $(0.06\,Re_m + 1)^2$, and $V_{rm}$ collapses to exactly 1. This is the
property that makes the closure need no blend rule, and it holds only for a
particular reading of the coefficients."""))

cells.append(code('''RE_ID = np.logspace(-8, 6, 80)   # capped at 1e6: above that the perfect square
                                 # (0.06 Re + 1)^2 loses the +1 to double-precision
                                 # cancellation, which is arithmetic, not physics


def check1(p=P0):
    return float(np.max(np.abs(V_rm(1.0, RE_ID, p) - 1.0)))


def V_rm_split_lin(eps_g, Re_m, lin_outer, lin_sq):
    """eq. (12) with the two appearances of 0.06 read differently."""
    A, B = A_of(eps_g), B_of(eps_g)
    return 0.5 * (A - lin_outer * Re_m
                  + np.sqrt((lin_sq * Re_m) ** 2 + 0.12 * Re_m * (2 * B - A) + A * A))


def V_rm_b_prefactor_both(eps_g, Re_m):
    """eq. (14) mis-read as carrying the 0.8 on BOTH branches."""
    A = A_of(eps_g)
    B = np.where(eps_g <= P0["bsw"], P0["bco"] * eps_g ** P0["bde"],
                 P0["bco"] * eps_g ** P0["bdi"])
    return 0.5 * (A - 0.06 * Re_m + np.sqrt((0.06 * Re_m) ** 2 + 0.12 * Re_m * (2 * B - A) + A * A))


B1 = [("as printed", check1()),
      ("outer 0.5 -> 0.55", check1(dict(P0, half=0.55))),
      ("0.12 -> 0.10 in the cross term", check1(dict(P0, cross=0.10))),
      ("0.06 read as 0.05 in the linear term only",
       float(np.max(np.abs(V_rm_split_lin(1.0, RE_ID, 0.05, 0.06) - 1)))),
      ("the 0.8 prefactor put on the dilute branch too",
       float(np.max(np.abs(V_rm_b_prefactor_both(1.0, RE_ID) - 1)))),
      ("A exponent 4.14 -> 4.65", check1(dict(P0, aexp=4.65))),
      ("dense B exponent 1.28 -> 2.28", check1(dict(P0, bde=2.28))),
      ("dilute B exponent 2.65 -> 2.39", check1(dict(P0, bdi=2.39))),
      ("B prefactor 0.8 -> 0.7", check1(dict(P0, bco=0.7))),
      ("switch voidage 0.85 -> 0.80", check1(dict(P0, bsw=0.80))),
      ("both 0.06 and 0.12 read a decimal place out",
       check1(dict(P0, lin=0.6, cross=1.2)))]
t1 = pd.DataFrame(B1, columns=["reading", "max |V_rm(1) - 1| over Re"])
print(t1.to_string(index=False, float_format=lambda v: f"{v:.3e}"))
CHK1 = B1[0][1]
TOL1 = max(CHK1, 1e-14) * 1.000001
blind1 = [n for n, v in B1[1:] if v <= TOL1]
display(Markdown(
    f"**Result: {CHK1:.1e}.** The radicand really is a perfect square, so the identity is "
    f"exact in exact arithmetic; what is left is the double-precision cancellation of "
    f"$(0.06Re+1)$ against $0.06Re$ at the top of the Reynolds range, and it shrinks with "
    f"the range. It is therefore read as a *threshold*, and anything at or below "
    f"{TOL1:.1e} counts as undetected.\\n\\n"
    f"**What it catches:** the outer 0.5, the pairing of the two 0.06's, the 0.12, and a "
    f"spurious prefactor on the dilute branch of eq. (14) — the last of which is precisely "
    f"the mis-reading the mangled text layer invites, since it renders that branch as "
    f"`_g-2'6s`.\\n\\n"
    f"**What it is blind to, measured:** {len(blind1)} of the {len(B1)-1} mis-readings above "
    f"leave it at or below that threshold — **every voidage exponent** (4.14, 1.28, 2.65), "
    f"the 0.8 itself, and the switch voidage. All are invisible for the same trivial reason: "
    f"$1^x = 1$, so eqs. (13) and (14) return 1 at $\\\\epsilon_g = 1$ whatever their "
    f"exponents are. A decimal-place slip that scales 0.06 and 0.12 *together* is very nearly "
    f"invisible too — {B1[-1][1]:.1e}, which is {B1[2][1]/B1[-1][1]:.0e} times below the "
    f"smallest genuine catch — because the perfect square survives any common scaling of the "
    f"pair. **The exponents are settled by Check 3 and by the page image, not by this.**")) '''))

cells.append(md(r"""### Check 2 — the force balance returns $v_{\text{slip}} = V_{rm}\,v_t$

This is the design property of eq. (11) and the reason $V_{rm}$ is called a
terminal-velocity ratio, but the report does not state it. Substituting eq. (11)
into the balance $F_{gm}v_{\text{slip}} = \epsilon_s\epsilon_g(\rho_s-\rho_g)g$
derived in the Background, the $\epsilon_s\epsilon_g$ cancels, $V_{rm}^2$ absorbs
$v_{\text{slip}}^2$, and what remains is the single-sphere terminal-velocity
equation in $v_{\text{slip}}/V_{rm}$.

**It is an algebraic identity of eq. (11), it does not test any constant inside
eq. (12) or eq. (16), and it is here because the $V_{rm}^2$ and the
$Re_m/V_{rm}$ argument are the two things in that equation most likely to be
mis-transcribed.**"""))

cells.append(code('''CASES = [(50e-6, 2500.0, 0.40), (300e-6, 2500.0, 0.60),
         (3e-3, 2500.0, 0.90), (1e-3, 1200.0, 0.50)]
RHO_A, MU_A, G_A = 1.2, 1.8e-5, 9.81


def check2(F=F_syamlal_obrien, p=P0):
    worst = 0.0
    for d_p, rho_s, e in CASES:
        drho = rho_s - RHO_A
        vt = v_terminal(d_p, drho, RHO_A, MU_A, G_A, p)
        v_bal = brentq(lambda v: F(e, v, d_p, RHO_A, MU_A, p) * v - (1 - e) * e * drho * G_A,
                       1e-16, 1e6, xtol=1e-18, rtol=1e-15)
        target = V_rm_at_balance(e, d_p * vt * RHO_A / MU_A, p) * vt
        worst = max(worst, abs(v_bal - target) / target)
    return worst


def _mk(power=2, cd_arg="ratio", pref=None, drop_eps_g=False):
    pf = P0["pref"] if pref is None else pref
    def F(eps_g, v, d_p, rho, mu, p=P0):
        Re = d_p * v * rho / mu
        Vr = V_rm(eps_g, Re, p)
        arg = Re / Vr if cd_arg == "ratio" else Re
        eg = 1.0 if drop_eps_g else eps_g
        return pf * (1 - eps_g) * eg * rho / (Vr ** power * d_p) * C_Ds(arg, p) * v
    return F


B2 = [("as printed", check2()),
      ("V_rm^2 -> V_rm in eq. (11)", check2(_mk(power=1))),
      ("V_rm^2 -> V_rm^3", check2(_mk(power=3))),
      ("C_Ds argument Re_m/V_rm -> Re_m", check2(_mk(cd_arg="plain"))),
      ("3/4 -> 1", check2(_mk(pref=1.0))),
      ("eps_g dropped from eq. (11)", check2(_mk(drop_eps_g=True))),
      ("A exponent 4.14 -> 4.65", check2(p=dict(P0, aexp=4.65))),
      ("C_D slope 4.8 -> 4.0", check2(p=dict(P0, cd1=4.0))),
      ("C_D offset 0.63 -> 0.44", check2(p=dict(P0, cd0=0.44))),
      ("outer 0.5 -> 0.55", check2(p=dict(P0, half=0.55)))]
t2 = pd.DataFrame(B2, columns=["reading", "worst relative residual"])
print(t2.to_string(index=False, float_format=lambda v: f"{v:.3e}"))
CHK2 = B2[0][1]
display(Markdown(
    f"**Result: {CHK2:.2e} over four decades of particle size and three voidages.**\\n\\n"
    f"**What it catches, measured:** the exponent on $V_{{rm}}$ (a factor of "
    f"{B2[1][1]:.1f} wrong if it is 1 instead of 2), the argument of $C_{{Ds}}$ "
    f"({100*B2[3][1]:.0f} %), the 3/4 ({100*B2[4][1]:.0f} %) and the $\\\\epsilon_g$ "
    f"({100*B2[5][1]:.0f} %).\\n\\n"
    f"**What it is blind to, measured:** every constant inside eq. (12) and eq. (16) — the "
    f"last four rows all sit at {max(v for _, v in B2[-4:]):.0e}, i.e. at the same machine "
    f"noise as the correct reading. They cancel because the identity uses the *same* "
    f"$C_{{Ds}}$ on both sides and the *same* $V_{{rm}}$ on both sides. **A page that "
    f"presented this residual as evidence that the drag law is right would be claiming "
    f"something it cannot see.**")) '''))

cells.append(md(r"""### Check 3 — against the velocity–voidage index the closed form replaces

This one can fail, and the report itself sets it up: eq. (12) is adopted
*because* Richardson and Zaki's correlation "can be calculated only numerically".
So the question with a number attached is **what the substitution costs**.

Using the bridge derived in *The published model*, the index implied by
Syamlal–O'Brien is $n_{\text{implied}} = \ln(\epsilon_g V_{rm})/\ln\epsilon_g$,
evaluated at the single-particle Reynolds number Table VI is indexed by. Two
honest complications are handled rather than hidden. First, eq. (12) is **not** a
power law, so its implied index depends on voidage — the table below reports it
at several, and the spread *is* part of the answer. Second, the two correlations
were fitted to overlapping literature, so this is not two independent witnesses;
it is the cost of one closed form standing in for another."""))

cells.append(code('''def implied_index(eps_g, Re_t, p=P0):
    return np.log(eps_g * V_rm_at_balance(eps_g, Re_t, p)) / np.log(eps_g)


EPS_PROBE = [0.45, 0.50, 0.60, 0.70, 0.85, 0.95]
tab = pd.DataFrame({"Re": rz.Re.to_numpy(), "n_0 (tabulated)": rz.n_0.to_numpy()})
for e in EPS_PROBE:
    tab[f"n at eps={e:.2f}"] = [implied_index(e, r) for r in rz.Re]
print("Richardson & Zaki's index, and the index eq. (12) implies")
print(tab.round(3).to_string(index=False))

dev = pd.DataFrame({"eps_g": EPS_PROBE})
dev["bias %"] = [100 * np.mean((tab[f"n at eps={e:.2f}"] - tab["n_0 (tabulated)"])
                               / tab["n_0 (tabulated)"]) for e in EPS_PROBE]
dev["MAD %"] = [100 * np.mean(np.abs((tab[f"n at eps={e:.2f}"] - tab["n_0 (tabulated)"])
                                     / tab["n_0 (tabulated)"])) for e in EPS_PROBE]
dev["rows above tabulated, of 18"] = [int((tab[f"n at eps={e:.2f}"] > tab["n_0 (tabulated)"]).sum())
                                     for e in EPS_PROBE]
print("\\n" + dev.round(2).to_string(index=False))
BIAS50 = float(dev.loc[dev.eps_g == 0.50, "bias %"].iloc[0])
BIAS95 = float(dev.loc[dev.eps_g == 0.95, "bias %"].iloc[0])
BIAS45 = float(dev.loc[dev.eps_g == 0.45, "bias %"].iloc[0])
d45 = 100 * (tab["n at eps=0.45"] - tab["n_0 (tabulated)"]) / tab["n_0 (tabulated)"]
d95 = 100 * (tab["n at eps=0.95"] - tab["n_0 (tabulated)"]) / tab["n_0 (tabulated)"]
N45_HI = int((d45 > 0).sum()); N95_HI = int((d95 > 0).sum())

fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.plot(rz.Re, rz.n_0, "o", ms=7, color=COL["meas"], mec="white", mew=1.2,
        label="Richardson & Zaki, Table VI (18 rows)", zorder=5)
for e, ls in zip([0.45, 0.60, 0.85], ["-", "--", ":"]):
    ax.plot(rz.Re, [implied_index(e, r) for r in rz.Re], ls, lw=2, color=COL["so"],
            label=f"eq. (12), $\\\\epsilon_g$ = {e:.2f}")
ax.axhline(4.65, color=COL["rz"], lw=1.2, ls="-.")
ax.text(0.45, 4.70, "their eq. (33), $Re<0.2$: $n=4.65$", fontsize=8, color=COL["rz"])
ax.axhline(2.39, color=COL["rz"], lw=1.2, ls="-.")
ax.text(0.45, 2.44, "their eq. (34), $Re>500$: $n=2.39$", fontsize=8, color=COL["rz"])
ax.set_xscale("log"); ax.set_xlabel("$Re$ on the single-particle terminal velocity")
ax.set_ylabel("velocity–voidage index $n$")
ax.set_title("The closed form against the correlation it stands in for")
ax.legend(frameon=False, fontsize=9)
plt.tight_layout(); plt.show()

display(Markdown(
    f"**Result, stated the way the numbers actually fall.** The mean bias is "
    f"**{BIAS45:+.1f} %** at $\\\\epsilon_g = 0.45$, **{BIAS50:+.1f} %** at 0.50 and "
    f"**{BIAS95:+.1f} %** at 0.95 — it **grows monotonically with voidage**, which it must, "
    f"because eq. (12) is not a power law and Richardson and Zaki's relation is. That "
    f"drift is the structural cost of the substitution and the report does not mention it.\\n\\n"
    f"But the sign is **not** uniform at the dense end, and saying otherwise would "
    f"misdescribe the table. At $\\\\epsilon_g = 0.45$ the implied index is above the "
    f"tabulated one on {N45_HI} of the 18 rows and below it on {18-N45_HI}, running from "
    f"{d45.min():+.1f} % to {d45.max():+.1f} % — a mean of {BIAS45:+.1f} % on a "
    f"{d45.max()-d45.min():.0f}-point spread, i.e. a **shape** mismatch as much as an "
    f"offset. The last column of the summary counts the rows on which it sits above; that "
    f"count first reaches all 18 at $\\\\epsilon_g = "
    f"{dev.loc[dev['rows above tabulated, of 18'] == 18, 'eps_g'].min():.2f}$ and stays "
    f"there ({N95_HI} of 18 at 0.95).\\n\\n"
    f"In the viscous limit the arithmetic is a one-liner: eq. (13) gives "
    f"$\\\\epsilon_gV_{{rm}} \\\\to \\\\epsilon_g^{{{P0['aexp']:g}+1}}$, so eq. (12) implies "
    f"$n = {P0['aexp']+1:g}$ where their eq. (33) gives 4.65 at $d/D = 0$ — "
    f"**{100*((P0['aexp']+1)/4.65-1):+.1f} %**, with no fitting and no window choice "
    f"anywhere in it. **Read the next cell before using that number**: the report's own "
    f"prose and its own eq. (11) put it one factor of $\\\\epsilon_g$ apart, and the other "
    f"reading reverses its sign.")) '''))

cells.append(md(r"""#### The one identification this rests on, and the reading that reverses its sign

Every number in Check 3 depends on a single question — **what is $V_{rm}$ a ratio
of?** — and the report answers it twice, inconsistently.

- **Eq. (11) answers it unambiguously.** Substituting it into the force balance
  gives $v_{\text{slip}} = V_{rm}v_t$ and nothing else; Check 2 measures that to
  machine precision. So *as implemented*, $V_{rm}$ is a **slip**-velocity ratio,
  Richardson and Zaki's superficial $V_c = \epsilon_g v_{\text{slip}}$, and the
  bridge carries the factor $\epsilon_g$: $n = \ln(\epsilon_g V_{rm})/\ln\epsilon_g$.
- **The report's prose answers it the other way.** It says $V_{rm}$ "can be
  calculated from the Richardson and Zaki (1954) correlation only numerically" —
  and *their* correlation is $V_c/V_i = \epsilon^n$ in the **superficial**
  velocity, which their own Summary states explicitly and this page's sidecar
  records in capitals. Read that way $V_{rm}$ *is* $\epsilon^n$, with no
  $\epsilon_g$ in the bridge and $n = \ln V_{rm}/\ln\epsilon_g$.

The two differ by exactly one unit of index, which is the error this page's
*Published model* section already names as "the same size as the whole effect
being measured". Nothing on disk settles it: Garside and Al-Dibouni (1977), whose
$V_r$ it actually is, is not here, and the report never writes down what its own
$V_{rm}$ is a ratio of. **So both are printed, and the page says which one it
adopts and why.**"""))

cells.append(code('''def alt_index(eps_g, Re_t, p=P0):
    """The index under the SUPERFICIAL reading: V_rm itself is eps^n, no eps_g
    in the bridge. Exactly one unit below implied_index, by construction."""
    return np.log(V_rm_at_balance(eps_g, Re_t, p)) / np.log(eps_g)


both = pd.DataFrame({"eps_g": EPS_PROBE})
both["page's reading: n = ln(eps V_rm)/ln eps"] = [
    100 * np.mean((np.array([implied_index(e, r) for r in rz.Re]) - rz.n_0) / rz.n_0)
    for e in EPS_PROBE]
both["alternative: n = ln(V_rm)/ln eps"] = [
    100 * np.mean((np.array([alt_index(e, r) for r in rz.Re]) - rz.n_0) / rz.n_0)
    for e in EPS_PROBE]
both.columns = ["eps_g", "bias %, slip reading", "bias %, superficial reading"]
print("Bias against Table VI's 18 indices, under both readings of V_rm")
print(both.round(2).to_string(index=False))

ALT45 = float(both.loc[both.eps_g == 0.45, "bias %, superficial reading"].iloc[0])
ALT50 = float(both.loc[both.eps_g == 0.50, "bias %, superficial reading"].iloc[0])
ALT95 = float(both.loc[both.eps_g == 0.95, "bias %, superficial reading"].iloc[0])
VISC_PAGE = 100 * ((P0["aexp"] + 1) / 4.65 - 1)
VISC_ALT = 100 * (P0["aexp"] / 4.65 - 1)
UNIT = max(abs(implied_index(e, r) - alt_index(e, r) - 1.0)
           for e in (0.45, 0.95) for r in rz.Re)
display(Markdown(
    f"""**The two readings differ by exactly one unit of index** — verified to {UNIT:.1e}
over the table — so they are a rigid offset of $\\\\overline{{1/n_0}}$ in per cent, and
every *sensitivity* measured in the break table below is identical under both. What is
not identical is **the sign of the answer**:

| | slip reading (adopted) | superficial reading |
|---|---|---|
| bias vs Table VI at $\\\\epsilon_g = 0.45$ | **{BIAS45:+.2f} %** | {ALT45:+.2f} % |
| at $\\\\epsilon_g = 0.50$ | **{BIAS50:+.2f} %** | {ALT50:+.2f} % |
| at $\\\\epsilon_g = 0.95$ | **{BIAS95:+.2f} %** | {ALT95:+.2f} % |
| viscous limit vs their tabulated 4.65 | **{P0['aexp']+1:g}, {VISC_PAGE:+.1f} %** | {P0['aexp']:g}, {VISC_ALT:+.1f} % |

**Which this page adopts, and why.** The slip reading, because it is the one
**eq. (11) implements** — it is forced by the equation MFIX actually evaluates, and
Check 2 proves it rather than assuming it. Everything Check 3 reports is therefore *the
index this closure imposes on a simulation*, which is the quantity a modeller needs.

**What the page does not claim.** That this is the index Garside and Al-Dibouni
*fitted*. The evidence there points the other way and is worth stating against the
page's own result: under the superficial reading the closed form lands within
**{abs(ALT95):.2f} %** of Table VI at $\\\\epsilon_g = 0.95$ and within
**{abs(float(both.loc[both.eps_g==0.85,'bias %, superficial reading'].iloc[0])):.2f} %** at
0.85 — near-exact agreement in exactly the dilute regime where a hindered-settling
correlation is best founded and where they would have fitted. That is suggestive, not
conclusive, and it cannot be settled without their 1977 paper.

**So the honest statement is conditional, and the conclusion depends on the reading.**
If eq. (11) means what it says, the substitution costs a systematically *larger* index,
growing with voidage. If instead Garside and Al-Dibouni's $V_r$ is the superficial ratio,
then their correlation reproduces Table VI well and **MFIX's eq. (11) is the thing that
carries the spare factor of $\\\\epsilon_g$** — a claim about the code, not about the
correlation, and one this page has no document to support. Both columns above are
CI-tracked so that neither can drift."""))'''))

cells.append(md(r"""#### The break table, and the defect that was inside it

**This table used to be run at $\epsilon_g = 0.50$ and nowhere else, and on the
strength of that this page claimed that "between the two checks every printed
constant of the closure is exercised by something". That claim was false**, and
it is the third time in this repository that a page's guard structure has itself
contained the class of defect it exists to catch. It is stated here rather than
quietly patched, because the failure is *structural* and will recur:

**eq. (14) is branched, so the voidage at which the check is evaluated decides
which of its constants are even live.** At $\epsilon_g = 0.50$ the dilute branch
is never reached, so the exponent 2.65 and the switch 0.85 could be set to
anything at all and this table would not twitch. Nor would any of the others:
Check 1 sits at $\epsilon_g = 1$, where $1^x = 1$; Check 2 is an algebraic
identity of eq. (11) that never reads eq. (14); and Check 4's Geldart voidages
are 0.441–0.498, every one of them on the dense branch. **Nothing on the page
evaluated the dilute branch against anything** — while the page's own headline
Check-3 number, $+29.2\,\%$ at $\epsilon_g = 0.95$, is computed *on* that branch
and is a CI-tracked metric. The check that would have caught the error was the
one the page publishes as its result.

The fix is to sweep the same table over voidage instead of fixing one. The
column headings below are the evaluation voidage; the entries are the *change*
in the bias, in points, against the as-printed row."""))

cells.append(code('''# --- break table for Check 3, swept over the evaluation voidage --------------
def check3_bias(p=P0, eps_g=0.50):
    d = np.array([implied_index(eps_g, r, p) for r in rz.Re])
    return 100 * np.mean((d - rz.n_0.to_numpy()) / rz.n_0.to_numpy())


B3 = [("A exponent 4.14 -> 4.65", dict(aexp=4.65)),
      ("A exponent 4.14 -> 4.00", dict(aexp=4.00)),
      ("dense B exponent 1.28 -> 1.39", dict(bde=1.39)),
      ("B prefactor 0.8 -> 1.0", dict(bco=1.0)),
      ("dilute B exponent 2.65 -> 2.39", dict(bdi=2.39)),
      ("dilute B exponent 2.65 -> 1.65", dict(bdi=1.65)),
      ("dilute branch deleted (B = 1 above the switch)", dict(bdi=0.0)),
      ("switch voidage 0.85 -> 0.80", dict(bsw=0.80)),
      ("switch voidage 0.85 -> 0.55", dict(bsw=0.55)),
      ("switch voidage 0.85 -> 0.95", dict(bsw=0.95)),
      ("BOTH former blind spots: 2.65 -> 1.0 and 0.85 -> 0.55", dict(bdi=1.0, bsw=0.55)),
      ("outer 0.5 -> 0.55", dict(half=0.55)),
      ("0.06 and 0.12 a decimal place out", dict(lin=0.6, cross=1.2)),
      ("C_D offset 0.63 -> 0.44", dict(cd0=0.44)),
      ("C_D slope 4.8 -> 4.0", dict(cd1=4.0))]

BASE3 = {e: check3_bias(P0, e) for e in EPS_PROBE}
t3 = pd.DataFrame([dict(reading=n,
                        **{f"{e:.2f}": check3_bias(dict(P0, **d), e) - BASE3[e]
                           for e in EPS_PROBE})
                   for n, d in B3])
print("as printed, bias vs Table VI in %:  "
      + "   ".join(f"eps={e:.2f}: {BASE3[e]:+.2f}" for e in EPS_PROBE))
print("\\nchange in that bias, in points, under each mis-reading:")
print(t3.to_string(index=False, float_format=lambda v: f"{v:+8.2f}"))

TOL3 = 0.01                       # points; below this the check saw nothing
_at50 = t3["0.50"].abs() > TOL3
_anywhere = t3[[f"{e:.2f}" for e in EPS_PROBE]].abs().max(axis=1) > TOL3
MISSED_BY_050 = int((_anywhere & ~_at50).sum())
NEVER_SEEN = list(t3.loc[~_anywhere, "reading"])
DIL_1_65 = float(t3.loc[t3.reading.str.contains("2.65 -> 1.65"), "0.95"].iloc[0])
DIL_DEL = float(t3.loc[t3.reading.str.contains("deleted"), "0.95"].iloc[0])
SW_055 = float(t3.loc[t3.reading.str.contains("0.85 -> 0.55"), "0.60"].iloc[0])
SW_095 = float(t3.loc[t3.reading.str.contains("0.85 -> 0.95"), "0.95"].iloc[0])
SW_080 = float(t3.loc[t3.reading.str.contains("0.85 -> 0.80")].iloc[0, 1:].abs().max())
JOINT = float(t3.loc[t3.reading.str.contains("BOTH")].iloc[0, 1:].abs().max())
A_UNIT = abs(check3_bias(dict(P0, aexp=5.14)) - check3_bias(dict(P0, aexp=4.14)))
CD_BLIND = float(t3.loc[t3.reading.str.startswith("C_D")].iloc[:, 1:].abs().to_numpy().max())
display(Markdown(
    f"""**What the sweep recovers.** **{MISSED_BY_050}** of the {len(B3)} mis-readings are
invisible in the $\\\\epsilon_g = 0.50$ column and visible somewhere else — that is the
measured size of the defect the single-voidage table had. The dilute exponent is now
exercised: 2.65 → 1.65 moves the 0.95 column by **{DIL_1_65:+.2f}** points and deleting
the dilute branch outright by **{DIL_DEL:+.2f}**. The switch voidage is exercised too,
though only when the mis-read switch **crosses an evaluation voidage**: 0.85 → 0.55 moves
the 0.60 column by **{SW_055:+.2f}** points and 0.85 → 0.95 moves the 0.95 column by
**{SW_095:+.2f}**. And the joint defect that defeated the old table — both former blind
spots injected at once — now moves it by up to **{JOINT:.1f}** points.

**What is still barely seen, and why that is a property of the closure rather than a gap
in the check.** A *small* displacement of the switch, 0.85 → 0.80, is worth at most
{SW_080:.3f} points anywhere in the table — detectable in principle, indistinguishable
from nothing in practice. It cannot be more: the two branches of eq. (14) meet at 0.85 to
**{JUMP_B:+.3f} %** (Results §4), so a switch moved a little way inside the region where
they nearly coincide is bounded by that continuity, whatever voidage you evaluate at. The
bound is measured, not assumed, and it is the honest reason this one constant is left to
the page image.

**And the check is exactly blind to eq. (16)**, at every voidage: the last two rows move
by **exactly zero** ({CD_BLIND:.1e} points), because the single-sphere drag curve cancels out of
an index comparing two hindered velocities to the same terminal velocity. Those two
constants are tested by Check 4 (the Dalla Valle slope is worth 63 points there) and by
nothing else on this page. So the corrected statement of joint sufficiency is:
**every printed constant of the closure is exercised by Check 3 or Check 4 except the
switch voidage, whose effect is bounded by the {JUMP_B:.3f} % continuity measured
directly.**

**Resolving power, unchanged and still modest.** One whole unit on the A exponent is
worth {A_UNIT:.1f} points at $\\\\epsilon_g = 0.50$; a 0.14 change in it is worth
{abs(float(t3.loc[t3.reading.str.contains("4.00"), "0.50"].iloc[0])):.2f} points, against a
row-to-row spread of {d45.max()-d45.min():.0f} points in the deviations themselves. This
check could not settle 4.14 against 4.00, and does not try to. The exponents are settled
by the page image; this measures the closure, not the transcription."""))'''))

cells.append(md(r"""### Check 4 — three closures against 21 measured minimum-fluidisation velocities

The strongest thing available: a measurement none of these closures ever saw.
Geldart fluidised 22 sieved cuts of three powders in a 5 cm column in 1973 and
read $U_0$ off the pressure-drop/velocity curve; 21 rows carry one.

At incipient fluidisation the force balance derived in the Background holds
exactly, so each closure gives a $u_{mf}$ with **no fitting of any kind**:

- **Ergun:** the positive root of its own two-term balance in $d_{sv}$, which is
  exactly sphericity-free (demonstrated on `A1.6`).
- **Syamlal–O'Brien:** $U_{mf} = \epsilon_{mf}\,V_{rm}\,v_t$, with $V_{rm}$ from
  the implicit balance and $v_t$ on the same eq. (16) curve.
- **Richardson–Zaki:** $U_{mf} = \epsilon_{mf}^{\,n_0(Re_t)}\,v_t$, with $n_0$
  interpolated in $\log Re$ from Table VI.

The headline is taken on the **eight Diakon cuts**: they are the only rows with a
sphericity (Geldart calls the powder spherical) and the only rows whose particle
density is printed exactly. That is `A1.6`'s choice and it is made for the same
reasons, one of which is sharper here — Syamlal–O'Brien has no shape factor at
all, so on a non-spherical powder it is not being asked a fair question."""))

cells.append(code('''G_CGS, MU_CGS, RHOF_CGS = 981.0, 1.8e-4, 1.2e-3        # Geldart's own g and mu; rho_f assumed


def u_mf_ergun(eps, d_sv, drho, k1=None, k2=None, rho=RHOF_CGS, mu=MU_CGS, g=G_CGS):
    k1 = K1_PRINT if k1 is None else k1
    k2 = K2_PRINT if k2 is None else k2
    a = k2 * rho * (1 - eps) / (eps ** 3 * d_sv)
    b = k1 * mu * (1 - eps) ** 2 / (eps ** 3 * d_sv ** 2)
    c = -(1 - eps) * drho * g
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)


def u_mf_so(eps, d_p, drho, p=P0, rho=RHOF_CGS, mu=MU_CGS, g=G_CGS):
    vt = v_terminal(d_p, drho, rho, mu, g, p)
    return eps * V_rm_at_balance(eps, d_p * vt * rho / mu, p) * vt


def n_table6(Re):
    """Table VI, interpolated in log Re, and CLIPPED to its own range [0.39, 489].
    Outside that range the index is held constant, which is an extrapolation the
    table does not license. The cell printing the 21 context rows below states how
    many rows it happens to, and on which."""
    return float(np.interp(np.log10(np.clip(Re, rz.Re.min(), rz.Re.max())),
                           np.log10(rz.Re.to_numpy()), rz.n_0.to_numpy()))


def u_mf_rz(eps, d_p, drho, p=P0, rho=RHOF_CGS, mu=MU_CGS, g=G_CGS):
    vt = v_terminal(d_p, drho, rho, mu, g, p)
    return eps ** n_table6(d_p * vt * rho / mu) * vt


def predict(sub, epscol="eps_0", p=P0, k1=None, k2=None,
            rho=RHOF_CGS, mu=MU_CGS, g=G_CGS, dscale=1e-4):
    d = sub.d_sv_um.to_numpy() * dscale
    drho = sub.rho_s.to_numpy() - rho
    U = sub.U_0.to_numpy()
    e = sub[epscol].to_numpy()
    out = pd.DataFrame({"d_sv um": sub.d_sv_um.to_numpy(), "eps_mf": e, "U_0 meas": U})
    out["v_t"] = [v_terminal(dd, dr, rho, mu, g, p) for dd, dr in zip(d, drho)]
    out["Re_t"] = d * out.v_t * rho / mu           # single particle, terminal
    out["Re_mf meas"] = d * (U / e) * rho / mu     # the bed itself, on the MEASURED slip
    out["Ergun"] = [u_mf_ergun(ee, dd, dr, k1, k2, rho, mu, g) for ee, dd, dr in zip(e, d, drho)]
    out["Syamlal-O'Brien"] = [u_mf_so(ee, dd, dr, p, rho, mu, g) for ee, dd, dr in zip(e, d, drho)]
    out["Richardson-Zaki"] = [u_mf_rz(ee, dd, dr, p, rho, mu, g) for ee, dd, dr in zip(e, d, drho)]
    for c in ["Ergun", "Syamlal-O'Brien", "Richardson-Zaki"]:
        out[c + " %"] = 100 * (out[c] - U) / U
    return out


CLOSURES = ["Ergun", "Syamlal-O'Brien", "Richardson-Zaki"]
res = predict(dia, "eps_0")
print("Eight spherical Diakon cuts, at the settled voidage A1.7 derives")
print(res.round(4).to_string(index=False))
summary = pd.DataFrame([(c, res[c + " %"].mean(), res[c + " %"].abs().mean(),
                         float(np.sqrt((res[c + " %"] ** 2).mean())))
                        for c in CLOSURES], columns=["closure", "bias %", "MAD %", "rms %"])
res_MB = predict(dia, "eps_MB")
summary["bias % at eps_MB"] = [res_MB[c + " %"].mean() for c in CLOSURES]
print("\\n" + summary.round(2).to_string(index=False))'''))

cells.append(code('''# --- reconciliation with A1.6, which owns this comparison for Ergun ---------
BIAS_ERG_MB = float(res_MB["Ergun %"].mean())


def re_mf_wen_yu(Ga, C1=33.7, C2=0.0408):
    return np.sqrt(C1 ** 2 + C2 * Ga) - C1


_d = dia.d_sv_um.to_numpy() * 1e-4
_drho = dia.rho_s.to_numpy() - RHOF_CGS
_Ga = _d ** 3 * RHOF_CGS * _drho * G_CGS / MU_CGS ** 2
_U_wy = re_mf_wen_yu(_Ga) * MU_CGS / (_d * RHOF_CGS)
DEV_WY = 100 * (_U_wy - dia.U_0.to_numpy()) / dia.U_0.to_numpy()

res_042 = predict(dia.assign(eps_0=0.42), "eps_0")
BIAS_ERG_042 = float(res_042["Ergun %"].mean())
BIAS_ERG_0 = float(res["Ergun %"].mean())
display(Markdown(
f"""**Reconciliation with `A1.6`, recomputed rather than quoted.** Three of the numbers on
this page also appear on that one, computed by different code. `A1.6`'s column is quoted
from its published metrics; this page's column is recomputed from the raw table:

| quantity | this page | `A1.6` publishes |
|---|---|---|
| exact Ergun balance vs the eight Diakon $U_0$, at Geldart's reported $\\\\epsilon_{{MB}}$ | **{BIAS_ERG_MB:+.2f} %** | 58.435 |
| the same, at the settled $\\\\epsilon_0$ derived from his two columns | **{BIAS_ERG_0:+.2f} %** | 42.060 |
| the same, at a textbook $\\\\epsilon_{{mf}} = 0.42$ | **{BIAS_ERG_042:+.2f} %** | 4.8185 |
| Wen and Yu's eq. (1) on the same eight rows | **{DEV_WY.mean():+.2f} %** | −25.214 |

All four reproduce. The third row is the one that matters, and an earlier version of this
page **left it out** while listing the other findings `A1.6` states about these exact
rows — which is precisely the omission `AGENTS.md` obligation 1 exists to prevent, since
it is the single strongest number arguing against this page's own conclusion.

**Read it as the argument it is.** At a commonly assumed $\\\\epsilon_{{mf}} = 0.42$ the
exact Ergun balance is only **{BIAS_ERG_042:+.1f} %** biased on these eight rows — for
practical purposes unbiased, and the reading that *would* crown a winner. It is put
beside this page's own numbers rather than buried.

**And the counter-argument, which is why the page still does not stop there.** That 0.42
is a textbook value, not this dataset's. `A1.7` derives
$\\\\epsilon_0 = 1-(1-\\\\epsilon_{{MB}})H_{{MB}}/H_0$ from Geldart's own two reported
columns and gets **{dia.eps_0.min():.3f}–{dia.eps_0.max():.3f}** on the Diakon cuts, with
{dia.eps_0.iloc[-1]:.3f} on the two rows that need no inference at all. So the +4.8 %
is bought by overriding the source's own voidage by
{100*(0.42/dia.eps_0.mean()-1):+.1f} % — and at *its* voidage the same balance is
{BIAS_ERG_0:+.1f} %. Every comparison on this page is at a stated voidage for exactly
this reason.

**Sphericity does not enter any of it.** Ergun written in $d_{{sv}}$ is exactly
$\\\\phi_s$-free (`A1.6` demonstrates it) and neither Syamlal–O'Brien nor
Richardson–Zaki has a shape factor to set. On these eight rows the *only* free
assumption separating the comparisons is the voidage.

The empirical $u_{{mf}}$ correlation, which is fitted to minimum-fluidisation data rather
than derived from a drag law, is **{DEV_WY.mean():+.1f} %** on these rows — and it is the
one comparison here with **no voidage in it at all**, Wen and Yu having absorbed the
voidage into their two fitted constants."""))'''))

cells.append(code('''fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.3))
for c, k in zip(CLOSURES, ["ergun", "so", "rz"]):
    ax[0].plot(res["U_0 meas"], res[c], "o-", ms=6, lw=1.6, color=COL[k], mec="white",
               mew=1.0, label=c)
lo, hi = 0.4, 13
ax[0].plot([lo, hi], [lo, hi], "--", color="0.35", lw=1)
ax[0].text(4.2, 3.4, "parity", fontsize=8, color="0.35")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("measured $U_0$  (cm s$^{-1}$)"); ax[0].set_ylabel("predicted $u_{mf}$  (cm s$^{-1}$)")
ax[0].set_title("Eight spherical Diakon cuts")
ax[0].legend(frameon=False, fontsize=9)

for c, k in zip(CLOSURES, ["ergun", "so", "rz"]):
    ax[1].plot(res["Re_t"], res[c + " %"], "o-", ms=6, lw=1.6, color=COL[k], mec="white",
               mew=1.0, label=c)
ax[1].axhline(0, color="0.35", lw=1, ls="--")
ax[1].set_xscale("log"); ax[1].set_xlabel("$Re_t$ of the particle")
ax[1].set_ylabel("(predicted $-$ measured) / measured, %")
ax[1].set_title("and how the deviation runs with particle Reynolds number")
ax[1].legend(frameon=False, fontsize=9)
plt.tight_layout(); plt.show()

rat_se = (res["Syamlal-O'Brien"] / res["Ergun"]).to_numpy()
rat_sr = (res["Syamlal-O'Brien"] / res["Richardson-Zaki"]).to_numpy()
# Is the S-O'/Ergun ratio actually monotone with size? Measure, do not assume.
D_SE = np.diff(rat_se)
SE_MONOTONE = bool(np.all(D_SE > 0))
SE_ARGMAX = int(np.argmax(rat_se))
SE_D = res["d_sv um"].to_numpy()

# The voidage each closure would need in order to be unbiased on these rows. This is
# an INVERSION of the measurement, not an agreement with it, and is labelled as one.
def _bias_at(eps, col):
    return float(predict(dia.assign(eps_0=eps), "eps_0")[col + " %"].mean())


EPS_IMPLIED = {c: brentq(lambda e: _bias_at(e, c), 0.20, 0.60, xtol=1e-10) for c in CLOSURES}
display(Markdown(
f"""**Result, and what it does and does not rank.** All three closures overpredict every
one of the eight measured velocities at the settled voidage: Ergun by
**{res['Ergun %'].mean():+.1f} %** on average, Syamlal–O'Brien by
**{res["Syamlal-O'Brien %"].mean():+.1f} %**, Richardson–Zaki by
**{res['Richardson-Zaki %'].mean():+.1f} %**.

**The blanket refusal to rank does not survive the reconciliation above, and is narrowed
rather than repeated.** What this dataset genuinely cannot settle is whether any closure
is *right*: the absolute bias is controlled almost entirely by a voidage Geldart never
measured directly, and moving it from the derived
{dia.eps_0.min():.3f}–{dia.eps_0.max():.3f} to a textbook 0.42 takes Ergun from
{BIAS_ERG_0:+.0f} % to {BIAS_ERG_042:+.1f} %. But the *ordering* is stable under every
voidage tried, and the cleanest way to see it is to invert the question — the voidage
each closure would need in order to be unbiased on these eight rows:

| closure | voidage that erases its bias |
|---|---|
| Ergun | **{EPS_IMPLIED['Ergun']:.4f}** |
| Richardson–Zaki | {EPS_IMPLIED['Richardson-Zaki']:.4f} |
| Syamlal–O'Brien | {EPS_IMPLIED["Syamlal-O'Brien"]:.4f} |

Ergun's is an ordinary minimum-fluidisation voidage; the other two are below the random
loose packing of spheres, which a bed at incipient fluidisation cannot be. **So on this
dataset Ergun is the least biased at every voidage tried, and that is a ranking.** What is
*not* supported is calling it validated: at the voidage this dataset's own columns imply
it is {BIAS_ERG_0:+.0f} % biased, `A1.6` finds it {BIAS_ERG_MB:+.0f} % at Geldart's
reported $\\\\epsilon_{{MB}}$, and an empirical $u_{{mf}}$ correlation with no drag law in
it beats all three. This is a numbered inversion of the measurement, not agreement with
it.

**And the gap between the two families — stated as the algebra it is, not as a
measurement.** The Syamlal–O'Brien / Ergun ratio runs from **{rat_se.min():.2f} to
{rat_se.max():.2f}**, rising over the first {SE_ARGMAX+1} cuts and **turning over on the
last**: {SE_D[SE_ARGMAX]:.0f} µm gives {rat_se.max():.4f} and {SE_D[-1]:.0f} µm gives
{rat_se[-1]:.4f}, a decrement of {D_SE[-1]:+.4f}. It is therefore **not monotone in
size or in $Re$** — an earlier version of this page said it was, and the eight numbers
printed two cells above contradicted it.

The Syamlal–O'Brien / Richardson–Zaki ratio stays inside
**{rat_sr.min():.2f}–{rat_sr.max():.2f}**, and **that one contains no measurement at
all.** $u_{{mf}}^{{SO}} = \\\\epsilon\\\\,V_{{rm}}v_t$ and
$u_{{mf}}^{{RZ}} = \\\\epsilon^{{n_0(Re_t)}}v_t$, so the terminal velocity **cancels
identically** out of their quotient, which is
$\\\\epsilon V_{{rm}}(\\\\epsilon,Re_t)/\\\\epsilon^{{n_0(Re_t)}}$ — a function of
$(\\\\epsilon_0, Re_t)$ alone, verified in the next cell. It is **Check 3's implied-index
comparison re-evaluated at Geldart's eight $(\\\\epsilon_0, Re_t)$ pairs**, in different
symbols, and it would read exactly the same if Geldart had never run the experiment. It
cannot detect an error in either terminal velocity, in the shared eq. (16), or in
Geldart's $U_0$ column. The S–O′/Ergun ratio is predicted/predicted in the same sense;
Geldart's $U_0$ enters only the per-cent columns. **Neither ratio is measurement-facing
evidence, and neither is offered as any.** What the measurement does carry is the three
bias columns above, and the ordering they imply."""))'''))

cells.append(code('''# --- the terminal velocity really does cancel out of S-O' / R-Z --------------
# Recompute the ratio from (eps_0, Re_t) alone, never forming v_t. If this agrees,
# the "0.74-1.12" above contains no velocity and cannot be a measured result.
rat_sr_no_vt = np.array([e * V_rm_at_balance(e, rt) / e ** n_table6(rt)
                         for e, rt in zip(res.eps_mf.to_numpy(), res.Re_t.to_numpy())])
VT_CANCELS = float(np.max(np.abs(rat_sr_no_vt / rat_sr - 1)))
print(f"S-O'/R-Z from the full expressions : {np.round(rat_sr, 6)}")
print(f"the same, with v_t never formed    : {np.round(rat_sr_no_vt, 6)}")
print(f"worst relative difference          : {VT_CANCELS:.3e}  on all {len(rat_sr)} rows")

# --- break table for Check 4 -------------------------------------------------
def bias_of(t): return {c: float(t[c + " %"].mean()) for c in CLOSURES}


B4 = [("as read", predict(dia, "eps_0")),
      ("d_sv taken as mm, not um -> cm  (x10)", predict(dia, "eps_0", dscale=1e-3)),
      ("mu a decimal place out (1.8e-5 poise)", predict(dia, "eps_0", mu=1.8e-5)),
      ("rho_f a decimal place out (1.2e-2)", predict(dia, "eps_0", rho=1.2e-2)),
      ("eps_mf = 0.383 for every row", predict(dia.assign(eps_0=0.383), "eps_0")),
      ("eps_mf = 0.42, the textbook value A1.6 uses", res_042),
      ("eps_mf = eps_MB, not the settled eps_0", predict(dia, "eps_MB")),
      ("eps_mf = 0.50 for every row", predict(dia.assign(eps_0=0.50), "eps_0")),
      ("A exponent 4.14 -> 4.65", predict(dia, "eps_0", p=dict(P0, aexp=4.65))),
      ("C_D slope 4.8 -> 4.0", predict(dia, "eps_0", p=dict(P0, cd1=4.0))),
      ("C_D offset 0.63 -> 0.44", predict(dia, "eps_0", p=dict(P0, cd0=0.44))),
      ("Ergun on A1.1's refitted constants", predict(dia, "eps_0", k1=K1_REFIT, k2=K2_REFIT))]
t4 = pd.DataFrame([dict(reading=n, **bias_of(t)) for n, t in B4])
print("bias in per cent, eight Diakon rows")
print(t4.to_string(index=False, float_format=lambda v: f"{v:+10.2f}"))

# index by the reading's name, never by row number - the table gains rows
R4 = t4.set_index("reading")


def d4(reading, closure):
    return abs(float(R4.loc[reading, closure]) - float(R4.loc["as read", closure]))


d_rho = d4("rho_f a decimal place out (1.2e-2)", "Syamlal-O'Brien")
d_cd = d4("C_D slope 4.8 -> 4.0", "Syamlal-O'Brien")
d_cd_rz = d4("C_D slope 4.8 -> 4.0", "Richardson-Zaki")
d_cd_er = d4("C_D slope 4.8 -> 4.0", "Ergun")
d_ref = d4("Ergun on A1.1's refitted constants", "Ergun")
d_e50 = d4("eps_mf = 0.50 for every row", "Syamlal-O'Brien")
d_e383 = d4("eps_mf = 0.383 for every row", "Syamlal-O'Brien")
d_route = abs(float(R4.loc["as read", "Syamlal-O'Brien"]) - float(R4.loc["as read", "Ergun"]))
display(Markdown(
f"""**The check has power, and the table says over what.** A unit slip moves it by
thousands of per cent. Replacing the derived per-row voidage
({dia.eps_0.min():.3f}–{dia.eps_0.max():.3f}) by a flat 0.50 moves it by
{d_e50:.0f} points, and by a flat 0.383 by {d_e383:.0f}
— **the single largest real uncertainty on this comparison**, which is why the
voidage is derived from Geldart's own two columns rather than assumed, and why every
alternative choice is reported rather than one. Crucially it moves under the constants
**Check 3 is blind to at every voidage**: the Dalla Valle slope is worth {d_cd:.0f}
points here.

**Joint sufficiency, restated correctly.** With the Check-3 table now swept over voidage,
every printed constant of the closure is exercised by Check 3 or by Check 4 **except the
switch voidage**, whose effect is bounded by the {JUMP_B:.3f} % continuity measured in
Results §4 and cannot exceed it for any switch displaced inside the region where the two
branches nearly coincide. An earlier version of this page asserted joint sufficiency
outright, on a Check-3 table run at one voidage; two constants were exercised by nothing,
and one of them carries this page's own headline.

Two rows are worth reading as controls rather than as defects, and the second needs a
correction. The Ergun column does not move at all when a Syamlal–O'Brien constant is
perturbed — {d_cd_er:.1e} points under the Dalla Valle slope — because the two closures
share no code. **The Richardson–Zaki column is a different matter: it moves
{d_cd_rz:.0f} points under the same perturbation**, because it takes its terminal velocity
from the same eq. (16) curve by construction. It is therefore *not* independent of
Syamlal–O'Brien, and the invariance control above should not be read as a general
independence claim. And the assumed gas density, the one input on this
page that Geldart does not print, is worth only {d_rho:.0f} points, so nothing here rests
on it. Ergun's own constants are worth less still: A1.1's refit moves the Ergun bias by
{d_ref:.2f} points, against a gap to Syamlal–O'Brien of {d_route:.0f} points — **the
disagreement between the two routes is {d_route/d_ref:.0f} times larger than what the
uncertainty in Ergun's own constants is worth.**"""))'''))

cells.append(code('''# --- Check 5: does Table VI survive its own round trip? ---------------------
r_Re = float(np.max(np.abs(np.log10(rz.Re) - rz.log_Re)))
r_n0 = float(np.max(np.abs(np.log10(rz.n_0) - rz.log_n0)))
bad = rz.copy()
bad.loc[4, "n_0"] = 41.7                       # a plausible slip: reading 4.17 as 41.7
r_bad = float(np.max(np.abs(np.log10(bad.n_0) - bad.log_n0)))
bad2 = rz.copy(); bad2.loc[10, "Re"] = 1.82    # 18.2 read as 1.82
r_bad2 = float(np.max(np.abs(np.log10(bad2.Re) - bad2.log_Re)))
display(Markdown(
f"""### Check 5 — Table VI against its own logarithm columns

Richardson and Zaki print $\\\\log_{{10}}Re$ and $\\\\log_{{10}}n_0$ beside $Re$ and $n_0$.
The columns are redundant, so recomputing them tests this page's transcription of all
{len(rz)} rows at once — and it is a check the *paper* pays for, in the sense of
`docs/handoff.md`.

| | |
|---|---|
| worst absolute residual on the $\\\\log_{{10}}Re$ column | **{r_Re:.1e}** |
| worst absolute residual on the $\\\\log_{{10}}n_0$ column | **{r_n0:.1e}** |
| the same, with 4.17 mis-read as 41.7 | {r_bad:.2f} |
| the same, with 18.2 mis-read as 1.82 | {r_bad2:.2f} |

Both residuals sit at the printed rounding, and a single-digit slip in either column moves
them by an order unity amount. **This is also how one cell was read at all**: the
$\\\\log n_0$ entry on the $Re = 2.02$ row is printed with a damaged 6, and
$\\\\log_{{10}}(4.17) = {np.log10(4.17):.4f}$ leaves 0.62 as the only two-decimal reading."""))'''))

cells.append(md(r"""### Blind spots — what this page does **not** claim

- **Nothing about Gidaspow's closure or about Wen–Yu drag.** Neither *drag law*
  is in the document on disk. The `ε = 0.8` discontinuity that a comparison of
  blended closures would be for is not computed here and is not asserted; the
  only switch measured on this page is Garside and Al-Dibouni's at 0.85, and it
  turns out to be continuous. (The report does cite Gidaspow extensively and
  adopts a Ding & Gidaspow expression elsewhere; what it does not print is a
  Gidaspow drag closure or a blend rule.)
- **The `B` switch voidage, 0.85, is the one printed constant no check on this
  page settles**, when it is displaced by only a little. The Check-3 sweep
  catches it when the mis-read switch crosses an evaluation voidage, but a small
  displacement is bounded by the +0.049 % continuity of eq. (14) at 0.85 and
  cannot be worth more than that. It is settled by the page image, and by the
  continuity measurement in Results §4 — not by any comparison against data.
- **Which of the two readings of `V_rm` the report meant.** The whole of Check 3
  is one unit of index away from its alternative, and that alternative reverses
  the sign of the result. Both are printed; the page adopts the one eq. (11)
  implements and says so, and Garside & Al-Dibouni (1977), which would settle it,
  is not on disk.
- **Neither of the two ratios in Check 4 is measurement-facing.** `v_t` cancels
  identically out of Syamlal–O'Brien / Richardson–Zaki, so that ratio is Check 3
  restated at Geldart's eight (ε₀, Re_t) pairs; Geldart's `U_0` enters only the
  per-cent columns. They are reported as algebra, not as evidence.
- **The Richardson–Zaki column is not independent of Syamlal–O'Brien.** Both take
  their terminal velocity from the same eq. (16), and the break table measures
  what that is worth: the Dalla Valle slope moves both, and Ergun neither.
- **`n_table6` holds Table VI's index constant outside Re 0.39–489.** Four of the
  21 context rows sit below the bottom; none of the eight Diakon rows does.
- **The inertial half of every closure is essentially untested against
  measurement.** Every Geldart row is a group A powder at a single-particle
  $Re_t$ and a bed $Re_{mf}$ that the Check 4 table prints, both small. Table VI reaches
  $Re = 489$, but it constrains eqs. (12)–(14) only — Check 3 is exactly blind to
  the single-sphere curve.
- **No claim that Ergun is the reference.** The `A1.6` reconciliation above is
  the reason. On this dataset every closure overpredicts, and which of the three
  readings `A1.6` leaves open is correct is not settled here either.
- **The 13 non-Diakon rows are reported but carry no sphericity**, and
  Syamlal–O'Brien has no shape factor to give them, so the whole-table numbers
  below are printed as context and are not the page's result.
- **Nothing on this page is digitised.** Both shipped datasets are transcribed
  tables; the one digitised dataset it loads is `A1.1`'s, and it is used only to
  bound Ergun's own constants.
- **`A1.1`'s 244 markers are not Ergun's 640**, and the refit inherits that
  page's stated recall limitation. It is used here only as a sensitivity, where a
  few scattered misses cannot matter."""))

cells.append(code('''res21 = predict(meas, "eps_0")
print("For completeness: all 21 rows with a measured U_0, 13 of which have no sphericity")
print(pd.DataFrame([(c, res21[c + " %"].mean(), res21[c + " %"].abs().mean())
                    for c in CLOSURES],
                   columns=["closure", "bias %", "MAD %"]).round(2).to_string(index=False))
print("\\nThe two cracking catalysts are porous, cohesive and printed with an approximate")
print("density, and their eps_0 runs far above Diakon's; the numbers above are context.")

# --- the clipping n_table6 does, declared rather than left in a docstring ----
_lo, _hi = float(rz.Re.min()), float(rz.Re.max())
CLIP21 = int(((res21.Re_t < _lo) | (res21.Re_t > _hi)).sum())
CLIP8 = int(((res.Re_t < _lo) | (res.Re_t > _hi)).sum())
print(f"\\nTable VI covers Re {_lo:g} to {_hi:g}, and n_table6 HOLDS THE INDEX CONSTANT "
      f"outside it.")
print(f"Of the 21 context rows, {CLIP21} fall below the bottom of the table "
      f"({', '.join(f'{v:.3f}' for v in sorted(res21.Re_t[res21.Re_t < _lo]))}), so their "
      f"Richardson-Zaki")
print(f"prediction uses an extrapolated index. Of the {len(res)} Diakon rows the page "
      f"takes its result on, {CLIP8} are clipped:")
print(f"their Re_t runs {res.Re_t.min():.3f} to {res.Re_t.max():.3f}, inside the table "
      f"throughout. The headline is unaffected;")
print("the context table is not, and that is why it is context.")

metrics = dict(
    check1_Vrm_unity_identity=CHK1,
    check2_force_balance_identity=CHK2,
    check3_index_bias_pct_eps050=BIAS50,
    check3_index_bias_pct_eps095=BIAS95,
    check3_viscous_index_vs_RZ_pct=100 * ((P0["aexp"] + 1) / 4.65 - 1),
    # the alternative (superficial) reading of V_rm - the sign of the result depends
    # on which one the report meant, so both are tracked
    check3_alt_reading_bias_pct_eps050=ALT50,
    check3_alt_reading_bias_pct_eps095=ALT95,
    check3_alt_reading_viscous_index_vs_RZ_pct=VISC_ALT,
    # the two constants the single-voidage break table could not see
    check3_break_dilute_exponent_2p65_to_1p65_pts_at_eps095=DIL_1_65,
    check3_break_switch_0p85_to_0p55_pts_at_eps060=SW_055,
    check4_ergun_bias_pct=float(res["Ergun %"].mean()),
    check4_syamlal_obrien_bias_pct=float(res["Syamlal-O'Brien %"].mean()),
    check4_richardson_zaki_bias_pct=float(res["Richardson-Zaki %"].mean()),
    check4_ergun_bias_pct_at_epsMB=BIAS_ERG_MB,
    check4_ergun_bias_pct_at_eps042=BIAS_ERG_042,
    check4_eps_erasing_ergun_bias=EPS_IMPLIED["Ergun"],
    check4_eps_erasing_syamlal_obrien_bias=EPS_IMPLIED["Syamlal-O'Brien"],
    check4_eps_erasing_richardson_zaki_bias=EPS_IMPLIED["Richardson-Zaki"],
    check4_vt_cancels_out_of_SO_over_RZ=VT_CANCELS,
    ratio_SO_over_Ergun_max_diakon=float(rat_se.max()),
    ratio_SO_over_Ergun_last_cut_diakon=float(rat_se[-1]),
    check5_tableVI_logRe_residual=r_Re,
    check5_tableVI_logn0_residual=r_n0,
    worst_ratio_SO_over_Ergun_at_eps044=float(W["worst ratio"]),
    B_branch_jump_pct_at_switch=JUMP_B,
    dalla_valle_vs_stokes_pct=100 * (STOKES_RATIO - 1),
)
_ = report_agreement("A1.8", metrics)'''))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**Nothing computational, and the page should say so plainly.** No operator is
assembled, no Newton solve is run, and every result above would be identical
without pymrm. That is the honest description of a closure-relation page, and it
is the same statement [`A1.1`](../A1.1-ergun-pressure-drop/) and
[`A1.6`](../A1.6-wen-yu-minimum-fluidisation/) make.

What the *page* adds, as against reading the report:

- **The two families are put on one axis.** The report names them in one
  paragraph and then implements one. Written as drag coefficients they can be
  divided, and the quotient is a function of $(\epsilon_g, Re_m)$ that behaves
  much worse in the middle than at either end.
- **$V_{rm}$ is identified as a slip-velocity ratio and proved to be one — and
  the report's own prose is shown to contradict that.** The report never says
  what $V_{rm}$ is a ratio of, eq. (11) forces one answer and the sentence
  introducing $V_{rm}$ implies the other, and the two are one factor of
  $\epsilon_g$ — one whole unit of index — apart. Both readings are printed with
  their numbers, because the sign of Check 3's result depends on which is right,
  and the document that would settle it is not on disk. Getting the factor from
  the momentum equations rather than by analogy is the reason that check means
  anything at all.
- **The substitution the report makes is costed, conditionally.** Eq. (12) exists
  because the Richardson–Zaki correlation is implicit, and nobody, including the
  report, says what that convenience costs. Under the reading eq. (11)
  implements it is a systematically larger index growing with voidage, whose
  viscous end is a one-line comparison of $4.14 + 1$ against 4.65.
- **Three closures meet a measurement none of them saw**, and the useful result
  is the shape of what it settles: not that any is right — the absolute bias is
  a statement about an unmeasured voidage — but that the *ordering* is stable,
  and that the voidage each would need in order to be unbiased separates them
  into one ordinary value and two below random loose packing.
- **Two published pages are audited for free.** Three of `A1.6`'s numbers on this
  dataset are recomputed here by independent code and all reproduce — including
  the one that argues against this page's own first conclusion — and `A1.1`'s
  refit is reproduced from its own marker file.
- **Two checks are shown to have been weaker than the page claimed.** The Check-3
  break table was run at one voidage and missed two of the eleven constants
  outright; it is now swept, the constants it recovers are measured, and the one
  it still cannot see is bounded rather than asserted."""))

# -------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

### Lifting the closure

`F_syamlal_obrien(eps_g, v_slip, d_p, rho_g, mu_g)` and
`F_ergun(...)` are vectorised over $(\epsilon_g, v_{\text{slip}})$ and take every
constant from a dictionary, so they drop into a two-fluid residual unchanged. In
pymrm terms $F_{gm}$ is a **pointwise coefficient on a phase-coupled layout** —
`(n_z, n_phase)` — and the momentum equations it couples are `S7`/`S10`, not
`S3`. Three things to carry across:

1. **$Re_m$ is built on the slip velocity, not the terminal velocity** (eq. 15).
   Eq. (12) is therefore implicit whenever the slip velocity is what you are
   solving for; `V_rm_at_balance` shows the one-dimensional root find, and in a
   full model it is simply part of the residual.
2. **The $\epsilon_g$ in the force balance is on the mixture buoyancy.**
   $F v_{\text{slip}} = \epsilon_s\epsilon_g(\rho_s-\rho_g)g$, not
   $\epsilon_s(\rho_s-\rho_g)g$. Dropping it scales the drag by
   $1/\epsilon_g$, and what that does to the balance is measured in Check 2's
   break table.
3. **Neither closure has a sphericity.** Ergun's can be given one by writing
   $d_p \to \phi_s d_p$, and then the balance in $d_{sv}$ is exactly
   $\phi_s$-free; Syamlal–O'Brien has nowhere to put one, and the report says
   there is no accepted way to add it.

### Choosing between them

The comparison on this page supports a narrow recommendation and no more:

- **In the two asymptotic limits at packed-bed voidage the two sit within tens of
  per cent of each other.** Results §1 prints the ratios; they are far larger
  than the uncertainty in Ergun's own constants, and small enough to live with.
- **In the transition regime at packed-bed voidage they differ by nearly a factor
  of two**, and that is where a fluidised powder sits. Which is right is not
  settled by anything on this page; that they differ, and by how much, is.
- **In dilute flow the packed-bed form is not usable at all** without a blend,
  and the numbers in Results §3 say how badly. Syamlal–O'Brien needs no blend.
- **Against a measured $u_{mf}$ all three overpredict at the voidage this
  dataset's own columns imply**, and an empirical $u_{mf}$ correlation with no
  drag law in it beats every one of them. If a minimum-fluidisation velocity is
  what you want, [`A1.6`](../A1.6-wen-yu-minimum-fluidisation/) is the page — and
  read its caveats, because they bound this page too.
- **On those rows Ergun is the least biased of the three at every voidage tried,
  and that is the one ordering this page will support.** Do not read it as
  validation: the absolute bias is controlled by a voidage nobody measured, and
  moving it from the derived value to a textbook 0.42 takes Ergun from +42 % to
  +4.8 % — `A1.6`'s number, recomputed here. The non-arbitrary form of the
  comparison is the inversion in Check 4: Ergun needs $\epsilon_{mf} = 0.415$ to
  be unbiased on these rows, Richardson–Zaki 0.371 and Syamlal–O'Brien 0.354, and
  the last two are below the random loose packing of spheres.
- **Do not lift either Check-4 ratio as measured corroboration.** The terminal
  velocity cancels identically out of Syamlal–O'Brien / Richardson–Zaki, so the
  0.74–1.12 band would read the same if the experiment had never been done. The
  measurement is in the bias columns, not in the ratios.

### What would complete the case

The catalogue asks for **Gidaspow / Syamlal–O'Brien / Wen–Yu** and this page
delivers one of the three plus two substitutes. Two documents would close it, and
neither can be replaced by memory:

- **Gidaspow, D. (1994), *Multiphase Flow and Fluidization*, Academic Press** —
  for the blend rule and the voidage at which it switches. Everything else needed
  to evaluate it is already on this page: the Ergun branch is built here, and the
  jump at the switch is a two-line computation once the rule is printed.
- **Wen, C. Y. and Yu, Y. H. (1966), *Mechanics of fluidization*, Chem. Eng.
  Progr. Symposium Series No. 62, 62, 100–111** — for the drag correlation. Note
  this is a **different paper** from the `A.I.Ch.E. Journal` communication on
  disk, which `A1.6` is built from and which contains no drag law; the Symposium
  Series paper is that communication's own reference 23.

### Related

[`A1.1`](../A1.1-ergun-pressure-drop/) (the packed-bed branch, and the source of
every Ergun constant here), [`A1.6`](../A1.6-wen-yu-minimum-fluidisation/) (the
$u_{mf}$ correlation, and the caveats that bound Check 4),
[`A1.7`](../A1.7-geldart-classification/) (the measured dataset and the derived
voidage), `A1.5` (Richardson–Zaki in its own right — the paper is on disk and
this page uses one of its tables), [`E1.2`](../E1.2-davidson-bubble/) and
[`E2.1`](../E2.1-kunii-levenspiel-bubbling-bed/) (fluidised-bed models that
consume $u_{mf}$ rather than a drag law)."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
Path("index.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print(f"wrote index.ipynb with {len(cells)} cells")
