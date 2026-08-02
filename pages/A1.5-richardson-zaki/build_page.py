#!/usr/bin/env python3
"""Generate index.ipynb for page A1.5 (Richardson and Zaki, 1954).

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "Richardson-Zaki: the exponent, the wall term, and what the paper's own tables say about both"
description: "u/u_t = eps^n is quoted everywhere; the 1954 paper it comes from is a data paper that prints 66 measured exponents, two wall laws, and a piecewise n(Re, d/D) fitted to them. This page transcribes eight tables off 600 dpi renders, quantifies the correlation as the fit residual it is, decides between the two wall laws on the authors' own numbers, and puts the flux function into a pymrm conservation law."
categories: [sec:A, struct:S1, tier:T0, data:tier2, phase:liquid-solid]
date: 2026-08-02
---

# Richardson-Zaki: the exponent, the wall term, and what the paper's own tables say about both

**Catalog ID:** `A1.5` · **Structures:** `S1` · **Tier:** T0 · **Data tier:** 2 (printed tables of the authors' own measurements)

Everybody quotes it as

$$\frac{u}{u_t} = \varepsilon^{\,n}.$$

The paper does not write that. It writes

$$\log V_c = n \log \varepsilon + \log V_i \tag{28}$$

with $V_i$ the intercept at $\varepsilon = 1$ — and then spends four pages
establishing that **$V_i$ is not $V_0$**. In sedimentation it is
($V_i = V_0$, eq. 40); in fluidisation it is not
($\log V_0 = \log V_i + d/D$, eq. 41). The famous form is the paper's
sedimentation special case with the wall term dropped.

Richardson and Zaki (1954) is a **data paper**. It prints:

- **Tables I and II** — 37 sedimentation runs, 15 particle/liquid pairs settled
  in two or three tubes each, with the measured slope $n$ and intercept.
- **Tables III and IV** — 29 fluidisation runs of spheres, in two columns.
- **Table V** — six non-spherical particles fluidised in water.
- **Tables VI and VII** — the extrapolation $n_0 = n|_{d/D \to 0}$ against $Re$,
  and the slope of $n$ against $d/D$.
- **Table VIII** — $n$, $K$ and $S$ for the non-spherical particles.

and only then the correlation everyone remembers, as **five** expressions
(eqs. 33, 34, 37, 38, 39) covering four Reynolds windows.

**What this page does.**

1. Transcribes all eight tables off 600 dpi renders and checks the transcription
   against four printed identities the tables themselves carry.
2. Measures the correlation against the **66 measured exponents** it was fitted
   to — and labels that a **fit residual, not a prediction**, with a null
   baseline beside it in every Reynolds window.
3. Tests the two wall laws, eq. (40) and eq. (41), against the intercept columns
   of the same tables — the one comparison on this page that is genuinely
   two-sided, because each law is falsifiable on the data set the other was
   fitted to.
4. Closes the internal identities the correlation must satisfy, including the
   one that could have failed and does not: eq. (39) at $Re = 500$ against the
   separately-stated eq. (34).
5. Puts $f(\phi) = V_0\,\phi\,(1-\phi)^n$ into a **pymrm** conservation law and
   shows what that solve does and does not test.

**What this page does not do.** It does not digitise a single figure. It makes
no claim about any curve in Figs. 5-24, and nothing here rests on a coordinate
read off a plot.""" ))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

Before 1954 the standard moves for hindered settling were to put the
*suspension's* density and viscosity into Stokes' law — Robinson (1926),
Steinour (1944), Hawksley (1951), all cited on the paper's first two pages.
Richardson and Zaki's opening argument is that this cannot be right:

> This cannot be true for a suspension of uniform particles all settling at the
> same rate, because each particle displaces its own volume of liquid as it
> settles. ... The effect of concentration on the resistance force encountered by
> a particle, for a given relative velocity, is attributable to the increase in
> the velocity gradient rather than to a change in viscosity.

So they drop the effective-property route entirely and do dimensional analysis
on the drag of a *constituent* particle, arriving at

$$\frac{V_c}{V_0} = f\!\left(\frac{V_0 d \rho}{\mu},\ \frac{d}{D},\ \varepsilon\right) \tag{22}$$

with two limits in which the Reynolds group drops out — streamline flow, where
$R$ is independent of $\rho$ (eq. 24), and fully turbulent flow, where $R$ is
independent of $\mu$ and $d$ (eq. 27). The correlation is then the experimental
determination of that function, and the tables are the experiment.

The other structural claim, which the sedimentation/fluidisation pairing is
designed to test, is that the two processes are the same problem:

$$V_s = \frac{V_c}{\varepsilon} \quad\text{and}\quad V_s = \frac{V_c'}{\varepsilon}
\ \Longrightarrow\ V_c = V_c'$$

— the settling velocity of a suspension relative to a fixed plane equals the
superficial liquid velocity needed to hold the same suspension at the same
voidage. Their Fig. 14 shows sedimentation and fluidisation of the same divinyl
benzene particles falling on one line, and Tables I-IV are the two halves of that
comparison.

`A1.5` is a `T0` case because $\varepsilon(u)$ is the closure every liquid-solid
fluidised-bed model needs, and because the flux function
$\phi\,V_0(1-\phi)^n$ is the standard constitutive input to batch-settling and
thickener models.""" ))

# ------------------------------------------------------------ published model
cells.append(md(r"""## The published model

### Provenance, and which pagination this page cites

Everything below was read from **Richardson, J. F. and Zaki, W. N.,
"Sedimentation and fluidisation: Part I", *Transactions of the Institution of
Chemical Engineers* **32**, 35-53 (1954)**.

The document on disk is **not the 1954 printing**. It is the verbatim reprint in
the ***Trans IChemE* Vol. 75, December 1997 Jubilee Supplement**, doi
[10.1016/S0263-8762(97)80006-8](https://doi.org/10.1016/S0263-8762(97)80006-8),
in which the paper occupies **reprint pages S82-S100** (19 pages). Every page
carries two page numbers: the reprint's own `S82`...`S100` in the running head,
and the original journal's footer `TRANS. INSTN CHEM. ENGRS, Vol. 32, 1954`.
**The original page numbers 35-53 are nowhere printed on the reprint.** They are
19 consecutive pages, which is consistent with the catalogue's "32, 35-53", but
this page has not verified that range from the document and does not assert it.

**Every "page" reference on this page is a reprint page, `S82`-`S100`**, and is
marked as such. The citation given is the 1954 original, because that is what the
reprint reproduces and what the running footer states; the reprint is recorded
separately as the text actually read.

### The scan, and how the exponents were settled

This is a pre-1980 scan of a 1954 letterpress printing. The mantissas render
cleanly at 600 dpi; **the superscript exponent glyphs do not** — `-2`, `-3` and
`-8` are frequently indistinguishable, which is the trap
[`docs/pdf-findings.md`](../../docs/pdf-findings.md) records for exactly this
class of document. The PDF text layer is worse still: it renders the paper's own
equation (1) as `F = 317 pVd` and Table I's title as
`Summal'lf of 1M RuulU Obtainellfrom SedifMl'ItalionEzperifMnts`. **It was not
used for any number.**

Rather than guess a damaged exponent, every one was settled from quantities
printed **in the same row**:

| identity | what it fixes |
|---|---|
| $d/D \times D = d$ | the exponent of $d$, and of $d/D$ |
| $\log_{10} V_0$, printed as its own column | the exponent of $V_0$ |
| $Re = V_0 d\rho/\mu$, from four printed columns | $V_0$, $d$, and the $Re$ exponent |
| Stokes' law $V_0 = d^2(\rho_s-\rho)g/18\mu$ — **which is what the Table I column is headed** | $V_0$ in Table I |

Section 1 runs all four as checks and reports where they fail. **No exponent was
repaired by inference from what looked plausible**, and every cell whose glyph
was not itself decisive is marked in the data sidecars.

### The correlation, as printed

**The form** (page S92, eq. 28), $n$ the slope and $\log V_i$ the intercept of
$\log V_c$ against $\log\varepsilon$:

$$\log V_c = n \log \varepsilon + \log V_i \tag{28}$$

**Viscous end** (page S94), from the runs with $Re < 0.2$:

$$n = 4{\cdot}65 + 19{\cdot}5\,\frac{d}{D} \tag{33}$$

**Turbulent end** (page S94), "at Reynolds numbers greater than about 500, $n$ is
independent both of $d/D$ and Reynolds number":

$$n = 2{\cdot}39 \tag{34}$$

**The extrapolated intercept** $n_0 = n|_{d/D\to 0}$, fitted to Table VI
(page S94):

$$n_0 = 4{\cdot}35\,Re^{-0{\cdot}03} \quad (0{\cdot}2 < Re < 1) \tag{35}$$
$$n_0 = 4{\cdot}45\,Re^{-0{\cdot}1} \quad (1 < Re < 500) \tag{36}$$

**The intermediate branches** (page S94-S95). Each is displayed twice: first as
$n_0/|(d/D)_{n=0}| \cdot Re^{-p}\,(d/D) + n_0$, then collected. The collected
coefficients are the authors' own division, and section 4 checks it:

$$n = \left(4{\cdot}35 + 17{\cdot}5\,\frac{d}{D}\right)Re^{-0{\cdot}03}
  \quad (0{\cdot}2 < Re < 1) \tag{37}$$
$$n = \left(4{\cdot}45 + 18\,\frac{d}{D}\right)Re^{-0{\cdot}1}
  \quad (1 < Re < 200) \tag{38}$$
$$n = 4{\cdot}45\,Re^{-0{\cdot}1} \quad (200 < Re < 500) \tag{39}$$

**Range of validity, stated by the authors on page S95**: "over a range of
Reynolds number from $2\times10^{-4}$ up to $7\times10^{3}$, and of $d/D$ from
zero up to $4\times10^{-2}$."

**The two wall laws.** For sedimentation (page S95):

$$V_i = V_0 \tag{40}$$

For fluidisation (page S95), "these lines cut the $\log V_i$ axis (at $d/D = 0$),
giving intercepts of approximately $\log V_0$":

$$\log V_0 = \log V_i + \frac{d}{D} \tag{41}$$

"This difference is attributable to the fact that in fluidisation a velocity
gradient is created in the liquid because of the drag exerted by the walls."

**Non-spherical particles** (page S96), $K$ Heywood's volumetric shape factor:

$$n = 2{\cdot}7\,K^{0{\cdot}16},\qquad
  K = \frac{\pi}{6}\cdot\frac{d_s^{\,3}}{d_p^{\,3}} \tag{42}$$

with $d_s$ the diameter of the sphere of equal volume and $d_p$ the diameter of
a circle of the same area as the particle's projected profile in its most stable
position.

### The recipe, in the authors' words (page S97)

> From the knowledge of the physical properties of the solid and liquid, the
> terminal falling velocity, $V_0$, of a single particle is calculated and,
> hence, the value of $V_i$ is obtained (equations (40), (41)). The Reynolds
> number $V_0 d\rho/\mu$ is then computed and the index $n$ is obtained from the
> appropriate equation (33, 34, 37, 38, 39). On substituting in the expression
> $V_c/V_i = \varepsilon^n$ the relation is obtained between $V_c$ ... and the
> porosity, $\varepsilon$.

Note that **$Re$ is built on $V_0$, not on the operating velocity**, so $n$ is a
property of the particle-liquid-tube combination and does not move as the bed
expands. Section 5 uses that.""" ))

# ----------------------------------------------------------- assumptions cell
cells.append(md(r"""## Parameters and assumptions

| quantity | value | where from |
|---|---|---|
| liquid temperature | 20 °C, "maintained" | equipment section, page S88 |
| fluidisation columns | 2.44 in. and 1.5 in. i.d. Perspex | page S88 |
| sedimentation tubes | 1.9, 2.8, 3.2, 4.8, 5 cm i.d. Pyrex, ~70 cm high | page S86 |
| particle sizes studied | $d > 100$ µm, uniformly sized spheres | summary, page S82 |
| $g$ | 981 cm s$^{-2}$ | CGS throughout; not printed |
| water at 20 °C | $\mu = 1$ cP, $\rho = 1$ g cm$^{-3}$ | **not printed**; used only to re-check Table III's $Re$ column, and flagged wherever used |

Four assumptions this page makes and states:

1. **The two fluidisation columns are 2.44 in. = 6.1976 cm and 1.5 in. =
   3.81 cm.** Not printed in cm anywhere. Fixed here by the fact that all 16
   printed $d/D$ values of Table III reproduce from the printed $d$ and those two
   diameters to better than 0.3 % — section 1 measures it.
2. **Table III's liquid is water at ambient**, $\mu = 1$ cP and $\rho = 1$. Used
   only for a transcription check, never for a result.
3. **$Re$ for branch selection is the printed $Re$ column**, not a recomputed
   one. The two agree to better than 0.5 % everywhere except two rows named in
   section 1.
4. **Deviations are $(\text{model} - \text{measured})/\text{measured}$**,
   as percentages, everywhere on this page. Where a quantity is a logarithm the
   residual is quoted in decades as well, because a 0.01 residual in
   $\log_{10} V$ is 2.3 % in $V$.""" ))

# ------------------------------------------------------------------ env cells
cells.append(code(r"""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code(r'''import sys, urllib.request
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
from IPython.display import Markdown, display
from pymrm import construct_div
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A1.5-richardson-zaki"
RNG = np.random.default_rng(20260802)

# cite_data() prints the sidecar's source.container, which carries the ORIGINAL 1954
# page range.  That range is inherited from the catalogue and was never read: the
# document on disk is the 1997 reprint, which numbers this paper S82-S100.  Every
# cite_data line on this page is printed with this flag under it, so the page never
# asserts a range it did not see.
CITE_NOTE = ("      [ '32, 35-53' above is the CATALOGUE's page range, inherited and NOT verified: "
             "the reprint\n        actually read numbers this paper S82-S100 and 35-53 appears "
             "nowhere on it. ]")
pd.set_option("display.width", 150)
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

G_CGS = 981.0          # cm/s2, CGS throughout the paper
D_244_CM = 2.44 * 2.54 # the 2.44 in. column, assumption 1
D_15_CM  = 1.50 * 2.54 # the 1.5 in. column,  assumption 1

# ---- the correlation, exactly as eqs. (33), (34), (37), (38), (39) print it.
#      Nothing here is fitted; every constant is transcribed.
def n_richardson_zaki(Re, dD):
    """Piecewise index n(Re, d/D).  Branch boundaries are the paper's own."""
    Re = np.asarray(Re, float); dD = np.asarray(dD, float)
    out = np.full(np.broadcast(Re, dD).shape, np.nan)
    out = np.where(Re <= 0.2,                    4.65 + 19.5 * dD,               out)  # (33)
    out = np.where((Re > 0.2) & (Re <= 1.0),    (4.35 + 17.5 * dD) * Re**-0.03,  out)  # (37)
    out = np.where((Re > 1.0) & (Re <= 200.0),  (4.45 + 18.0 * dD) * Re**-0.1,   out)  # (38)
    out = np.where((Re > 200.0) & (Re <= 500.0), 4.45 * Re**-0.1 + 0.0 * dD,     out)  # (39)
    out = np.where(Re > 500.0,                   2.39 + 0.0 * dD,                out)  # (34)
    return out

def rz_branch(Re):
    Re = np.asarray(Re, float)
    return np.select([Re <= 0.2, (Re > 0.2) & (Re <= 1), (Re > 1) & (Re <= 200),
                      (Re > 200) & (Re <= 500), Re > 500],
                     ["eq33", "eq37", "eq38", "eq39", "eq34"], default="?")

def dev_pct(model, measured):
    """The one deviation convention on this page: (model - measured)/measured."""
    return 100.0 * (np.asarray(model, float) - np.asarray(measured, float)) / np.asarray(measured, float)

metrics = {}
print("pymrm", pymrm.__version__ if hasattr(pymrm, "__version__") else "(version not exposed)")'''))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

Four datasets, all transcribed character by character from 600 dpi renders. The
sidecars carry the full transcription procedure, the list of cells whose exponent
came from a row identity, and every printed anomaly found.""" ))

cells.append(code(r'''sed = load_data("rz1954_sedimentation.csv",        page=PAGE)
flu = load_data("rz1954_fluidisation_spheres.csv", page=PAGE)
n0t = load_data("rz1954_n0_vs_re.csv",             page=PAGE)
nsp = load_data("rz1954_nonspherical.csv",         page=PAGE)
sed["flag"] = sed["flag"].fillna(""); flu["flag"] = flu["flag"].fillna("")
nsp["flag"] = nsp["flag"].fillna("")

print(cite_data(load_meta("rz1954_sedimentation.csv", page=PAGE)))
print(CITE_NOTE)
print()
print(f"Tables I + II  : {len(sed):3d} (run, tube) rows, "
      f"{sed.groupby(['particle','liquid','d_cm']).ngroups} particle/liquid groups, "
      f"Re {sed.Re_printed.min():.2e} to {sed.Re_printed.max():.2f}")
print(f"Tables III + IV: {len(flu):3d} (run, column) rows, "
      f"{flu.groupby(['particle','liquid']).ngroups} particle/liquid groups, "
      f"Re {flu.Re_printed.min():.2f} to {flu.Re_printed.max():.0f}")
print(f"Tables VI + VII: {n0t.n0_printed.notna().sum()} n0 values, "
      f"{n0t.slope_n_vs_dD_printed.notna().sum()} slope values")
print(f"Table V + VIII : {len(nsp)} particles, K {nsp.K_printed.min():.4f} to {nsp.K_printed.max():.3f}")
print(f"\nTOTAL measured exponents n: {len(sed) + len(flu)}")
print(f"d/D spans {min(sed.d_over_D_printed.min(), flu.d_over_D_printed.min()):.2e} to "
      f"{max(sed.d_over_D_printed.max(), flu.d_over_D_printed.max()):.3f}; "
      f"the paper claims validity to 4e-2.")'''))

# ------------------------------------------------- section 1 : transcription
cells.append(md(r"""### 1. Does the transcription hold up?

Four identities, all built from columns the paper prints side by side, so each is
a check on the *reading* rather than on the physics. They are **structural** in
the sense that they cannot tell us anything about sedimentation — but they are
not powerless, because a mis-read digit breaks them, and three of the four
locate printed errors, on named cells.

The Stokes check is the interesting one, because it is the only one that is
*supposed* to fail somewhere: Table I is headed "Reynolds Number $< 0.2$" and its
$V_0$ column is headed "Stokes' law", while Table II is headed "Reynolds Number
$> 0.2$" and its column is not. If Stokes reproduced both, the transcription
would be telling us the table titles are meaningless.""" ))

cells.append(code(r'''s = sed.copy()
s["dD_calc"]     = s.d_cm / s.D_cm
s["dD_dev"]      = dev_pct(s.d_over_D_printed, s.dD_calc)
s["Re_calc"]     = s.V0_cm_s * s.d_cm * s.rho_g_cm3 / (s.mu_cP * 1e-2)
s["Re_dev"]      = dev_pct(s.Re_printed, s.Re_calc)
s["logV0_res"]   = s.log10_V0_printed - np.log10(s.V0_cm_s)
s["V0_stokes"]   = s.d_cm**2 * (s.rho_s_g_cm3 - s.rho_g_cm3) * G_CGS / (18 * s.mu_cP * 1e-2)
s["stokes_dev"]  = dev_pct(s.V0_cm_s, s.V0_stokes)

f = flu.copy()
f["dD_calc"]   = f.d_cm / f.D_cm
f["dD_dev"]    = dev_pct(f.d_over_D_printed, f.dD_calc)
f["Re_calc"]   = f.V0_cm_s * f.d_cm * f.rho_g_cm3 / (f.mu_cP * 1e-2)
f["Re_dev"]    = dev_pct(f.Re_printed, f.Re_calc)
f["logV0_res"] = f.log10_V0_printed - np.log10(f.V0_cm_s)

bad_dD  = s.flag.str.contains("d_over_D_inconsistent")
clean   = s[~bad_dD]
runS    = clean[clean.dD_dev.abs() > 1]
rest    = clean[clean.dD_dev.abs() <= 1]
tI, tII = s[s.table == "I"], s[s.table == "II"]
rising  = tII.flag.str.contains("particles_rising")
tII_set, tII_ris = tII[~rising], tII[rising]

rows = [
    ["d/D x D = d  (sedimentation, the other 34 rows)", len(rest), f"{rest.dD_dev.abs().max():.2f} %"],
    ["d/D x D = d  (sedimentation, run S)",      len(runS), f"{runS.dD_dev.abs().max():.1f} %"],
    ["d/D x D = d  (sedimentation, runs H and N)", int(bad_dD.sum()), f"{s[bad_dD].dD_dev.abs().max():.1f} %"],
    ["d/D x D = d  (fluidisation, all rows)",    len(f), f"{f.dD_dev.abs().max():.2f} %"],
    ["Re = V0 d rho / mu  (sedimentation)",      len(s), f"{s.Re_dev.abs().max():.2f} %"],
    ["Re = V0 d rho / mu  (fluidisation)",       len(f), f"{f.Re_dev.abs().max():.2f} %"],
    ["log10 V0 column (sedimentation)",          len(s), f"{s.logV0_res.abs().max():.4f} decades"],
    ["log10 V0 column (fluidisation)",           len(f), f"{f.logV0_res.abs().max():.4f} decades"],
    ["Stokes' law vs Table I  V0 column",        len(tI), f"{tI.stokes_dev.abs().max():.2f} %"],
    ["Stokes' law vs Table II, settling rows",   len(tII_set),
     f"{tII_set.stokes_dev.abs().min():.0f} to {tII_set.stokes_dev.abs().max():.0f} %"],
    ["Stokes' law vs Table II, rising rows",     len(tII_ris),
     f"{tII_ris.stokes_dev.abs().max():.0f} % (WRONG SIGN)"],
]
print(pd.DataFrame(rows, columns=["identity", "rows",
                                  "|deviation| (worst, or range)"]).to_string(index=False))

metrics.update(
    transcription_dD_max_pct        = rest.dD_dev.abs().max(),
    transcription_Re_max_pct        = max(s.Re_dev.abs().max(), f.Re_dev.abs().max()),
    transcription_logV0_max_decades = max(s.logV0_res.abs().max(), f.logV0_res.abs().max()),
    stokes_tableI_max_pct           = tI.stokes_dev.abs().max(),
    stokes_tableII_settling_min_pct = tII_set.stokes_dev.abs().min(),
    stokes_tableII_settling_max_pct = tII_set.stokes_dev.abs().max(),
    column_244_in_cm                = D_244_CM,
    column_15_in_cm                 = D_15_CM,
)'''))

cells.append(code(r'''# The named failures, each printed with the number that exposes it.
print("PRINTED ANOMALIES FOUND BY THE FOUR IDENTITIES\n" + "-"*78)
for _, r in s[bad_dD].iterrows():
    print(f"  Table I run {r.runs:>2s}: d/D printed {r.d_over_D_printed:.3g} but d/D = "
          f"{r.d_cm:.3g}/{r.D_cm} = {r.dD_calc:.4g}  ({r.dD_dev:+.1f} %)")
odd = s[(s.dD_dev.abs() > 1) & ~bad_dD]
for _, r in odd.iterrows():
    print(f"  Table {r.table} run {r.runs:>2s}: d/D printed {r.d_over_D_printed:.3g} vs "
          f"{r.dD_calc:.4g}  ({r.dD_dev:+.1f} %)")
r12 = s[s.flag.str.contains("without_minus")].iloc[0]
print(f"\n  Table II run 12: intercept printed as {r12.intercept_log10_Vc_at_eps1:+.3f} while the other two "
      f"tubes of the same group and the log10 V0 column all read {r12.log10_V0_printed:+.3f}.")
print(f"    In the convention section 3 fits in, y = log10 V0 - intercept, the row reads "
      f"{r12.log10_V0_printed - r12.intercept_log10_Vc_at_eps1:+.3f} decades as printed and "
      f"{r12.log10_V0_printed + abs(r12.intercept_log10_Vc_at_eps1):+.3f} with the minus sign restored, "
      f"in line with every other Table II row. THE ROW IS KEPT AS PRINTED AND EXCLUDED FROM THE WALL FIT; "
      f"section 3 reports the fit with it excluded, as printed, and sign-restored.")
bad4 = f[f.flag.str.contains("V0_and_log10")].iloc[0]
print(f"\n  Table IV, ballotini {bad4.d_cm:g} cm in glycerol-water: V0 printed {bad4.V0_cm_s} with "
      f"log10 V0 printed {bad4.log10_V0_printed}, but log10({bad4.V0_cm_s}) = {np.log10(bad4.V0_cm_s):.4f}.")
print(f"    The printed Re column implies V0 = {bad4.Re_printed*bad4.mu_cP*1e-2/(bad4.d_cm*bad4.rho_g_cm3):.3f} "
      f"and the printed log implies {10**bad4.log10_V0_printed:.3f}; the V0 cell reads 7.35 at 600 dpi and is "
      f"recorded as printed. NOT REPAIRED - it enters no result on this page.")
g15 = f[f.flag.str.contains("d_printed_as_4.2x15", regex=False)]
grp15 = g15.groupby(["liquid", "Re_printed"])
print(f"\n  Table IV, glass spheres: d is printed '4.2 x 15^-1', a damaged '10'.  IT OCCURS "
      f"{grp15.ngroups} TIMES, once in each glass-sphere block of the table, affecting {len(g15)} "
      f"tabulated (row, column) entries:")
for (liq, Re_), gg in grp15:
    d_from_Re = Re_ * gg.mu_cP.iloc[0] * 1e-2 / (gg.V0_cm_s.iloc[0] * gg.rho_g_cm3.iloc[0])
    print(f"    {liq} at Re = {Re_:g} ({len(gg)} rows): the printed Re, V0, rho and mu return "
          f"d = {d_from_Re:.4f} cm")
print(f"    Table III prints 4.2e-1 for the same particle in water.  Recorded as 4.2e-1 on all "
      f"{len(g15)} rows and labelled a reconstruction on each.")
print(f"\n  Table I runs N and O: the printed Re is {s[s.runs=='N'].Re_printed.iloc[0]:.3g}, ABOVE the "
      f"table's own title of Re < 0.2.  V0 = {s[s.runs=='N'].V0_cm_s.iloc[0]:.3g} is confirmed both by the "
      f"printed log10 V0 and by Stokes' law to "
      f"{abs(s[s.runs=='N'].stokes_dev.iloc[0]):.2f} %, so Re really is ~0.3.  The row is assigned to the "
      f"eq. (37) branch by its printed Re, and section 2 reports what eq. (33) would have given.")'''))

cells.append(md(r"""**What the four identities establish, and what they cannot.**

They establish that the numbers on this page are the numbers on those pages, to
within the deviations printed above, and they locate five printed anomalies that
inspection alone would not have caught. They establish the two column diameters
in cm, which the paper never states.

They cannot establish anything about hindered settling. Stokes' law reproducing
Table I to a fraction of a per cent says only that the authors computed that
column the way its heading says. The same identity applied to Table II
over-predicts by 19 % to 66 % on the settling rows (the exact bounds are printed
in the table above), and on the three rows where
the particles *rise* ($\rho_s = 2.623$ in bromoform, $\rho = 2.89$) it returns
the wrong sign entirely — which is exactly what a table headed $Re > 0.2$ should
do. **Both outcomes are reported because the pair is the check**; either one
alone would be uninformative.""" ))

# --------------------------------------- section 2 : correlation vs the data
cells.append(md(r"""### 2. The correlation against the 66 measured exponents

Every row of Tables I-IV carries a measured slope $n$, a printed $Re$ and a
printed $d/D$. Eqs. (33), (34), (37), (38), (39) turn the last two into a
predicted $n$. That is 66 comparisons.

**This is a fit residual, not a prediction, and the page will not call it
anything else.** The five expressions were fitted to Fig. 18, and Fig. 18 was
built from these very slopes together with Steinour's and Lewis, Gilliland and
Bauer's — the paper says so on page S94. What the number measures is how tightly
the printed correlation reproduces the data it came from, which is worth knowing
(it bounds the transcription, it exposes the branch boundaries, and it can fail)
but is not evidence that the correlation predicts anything.

Two things on the same rows *are* informative:

- a **null baseline** in every Reynolds window — a constant $n$, fitted to the
  same rows, with no $Re$ and no $d/D$ in it;
- the rows the authors themselves **exclude**: "Some of the experimental results
  at higher values of $d/D$ do not follow the correlations given previously but
  they were obtained with oil, which caused the smaller glass spheres to
  agglomerate into flocs and to adhere to the tube walls" (page S95).""" ))

cells.append(code(r'''parts = []
for _, r in sed.iterrows():
    parts.append(dict(source=f"Table {r.table}", kind="sedimentation",
                      label=f"run {r.runs}: {r.particle} in {r.liquid}, D = {r.D_cm} cm",
                      group=f"{r.particle}|{r.liquid}|{r.d_cm}",
                      Re=r.Re_printed, dD=r.d_over_D_printed, n_meas=r.slope_n, is_oil=0,
                      rho_s=r.rho_s_g_cm3, rho=r.rho_g_cm3, mu=r.mu_cP, flag=r.flag))
for _, r in flu.iterrows():
    parts.append(dict(source=f"Table {r.table}", kind="fluidisation",
                      label=f"{r.particle} in {r.liquid}, D = {r.column_in} in",
                      group=f"{r.particle}|{r.liquid}",
                      Re=r.Re_printed, dD=r.d_over_D_printed, n_meas=r.slope_n, is_oil=int(r.is_oil),
                      rho_s=r.rho_s_g_cm3, rho=r.rho_g_cm3, mu=r.mu_cP, flag=r.flag))
A = pd.DataFrame(parts)
A["n_pred"] = n_richardson_zaki(A.Re.values, A.dD.values)
A["branch"] = rz_branch(A.Re.values)
A["dev"]    = dev_pct(A.n_pred, A.n_meas)

keep, oil = A[A.is_oil == 0], A[A.is_oil == 1]
print(f"{len(A)} comparisons: {len(keep)} in water or glycerol-water, {len(oil)} in the two oils.\n")
summary = []
for name, g in [("all 66", A), ("water / glycerol-water", keep), ("the two oils", oil)]:
    summary.append([name, len(g), f"{g.dev.abs().mean():.2f}", f"{np.sqrt((g.dev**2).mean()):.2f}",
                    f"{g.dev.abs().max():.2f}"])
print(pd.DataFrame(summary, columns=["subset", "rows", "mean |dev| %", "rms %", "max |dev| %"]).to_string(index=False))

per = []
for b, g in keep.groupby("branch"):
    const = g.n_meas.mean()
    dnull = dev_pct(np.full(len(g), const), g.n_meas.values)
    per.append([b, len(g), f"{g.Re.min():.3g}-{g.Re.max():.3g}", f"{g.dev.abs().mean():.2f}",
                f"{np.sqrt((g.dev**2).mean()):.2f}", f"{const:.3f}",
                f"{np.sqrt((dnull**2).mean()):.2f}"])
per = pd.DataFrame(per, columns=["branch", "rows", "Re range", "RZ mean|dev| %", "RZ rms %",
                                 "null constant n", "null rms %"])
print("\nPER-BRANCH, non-oil rows, with a constant-n null baseline fitted to the SAME rows:")
print(per.to_string(index=False))
c_all  = keep.n_meas.mean()
null_all = np.sqrt((dev_pct(np.full(len(keep), c_all), keep.n_meas.values)**2).mean())
print(f"\npooled over all {len(keep)} non-oil rows: one constant n = {c_all:.3f} gives "
      f"rms {null_all:.1f} % against the correlation's {np.sqrt((keep.dev**2).mean()):.2f} %.")

metrics.update(
    n_comparisons              = len(A),
    n_rows_nonoil              = len(keep),
    corr_nonoil_mean_abs_pct   = keep.dev.abs().mean(),
    corr_nonoil_rms_pct        = float(np.sqrt((keep.dev**2).mean())),
    corr_nonoil_max_abs_pct    = keep.dev.abs().max(),
    corr_oil_mean_abs_pct      = oil.dev.abs().mean(),
    corr_oil_max_abs_pct       = oil.dev.abs().max(),
    null_pooled_constant_n     = c_all,
    null_pooled_rms_pct        = null_all,
)'''))

cells.append(code(r'''# Where the null baseline is and is not beaten - the honest reading of the table above.
b38 = keep[keep.branch == "eq38"]; c38 = b38.n_meas.mean()
r38_rz   = np.sqrt((b38.dev**2).mean())
r38_null = np.sqrt((dev_pct(np.full(len(b38), c38), b38.n_meas.values)**2).mean())
b33 = keep[keep.branch == "eq33"]; c33 = b33.n_meas.mean()
r33_rz   = np.sqrt((b33.dev**2).mean())
r33_null = np.sqrt((dev_pct(np.full(len(b33), c33), b33.n_meas.values)**2).mean())
b34 = keep[keep.branch == "eq34"]; c34 = b34.n_meas.mean()
b37 = keep[keep.branch == "eq37"]; c37 = b37.n_meas.mean()
r37_rz   = np.sqrt((b37.dev**2).mean())
r37_null = np.sqrt((dev_pct(np.full(len(b37), c37), b37.n_meas.values)**2).mean())

display(Markdown(
    f"The correlation reproduces its own {len(keep)} non-oil rows to a mean "
    f"**{keep.dev.abs().mean():.2f} %** and a worst **{keep.dev.abs().max():.2f} %**. "
    f"Against a single constant $n = {c_all:.2f}$ that is a large margin — "
    f"{null_all:.0f} % rms — but the margin is **not evenly earned**:\n\n"
    f"- **Equation (38)**, the $1 < Re < 200$ branch and the largest at "
    f"{len(b38)} rows, is where the correlation does real work: rms "
    f"**{r38_rz:.2f} %** against the constant's **{r38_null:.1f} %**. Over that "
    f"window $n$ falls from about {b38.n_meas.max():.2f} to "
    f"{b38.n_meas.min():.2f}, and the $Re^{{-0.1}}$ tracks it.\n"
    f"- **Equation (33)**, the {len(b33)} viscous rows, is a real but much "
    f"smaller gain: rms **{r33_rz:.2f} %** against **{r33_null:.2f} %**. All the "
    f"information is in the $19.5\\,d/D$ term, and $d/D$ spans only "
    f"{b33.dD.min():.1e} to {b33.dD.max():.1e} here.\n"
    f"- **Equation (37)**, the {len(b37)} rows in $0.2 < Re < 1$, is **not "
    f"resolved either way**: rms {r37_rz:.2f} % against the constant's "
    f"{r37_null:.2f} %, a gap that reverses under all three of the sensitivities "
    f"below. Five rows over half a decade of $Re$ with an exponent of $-0.03$ "
    f"cannot be distinguished from a constant, and neither can be shown to beat "
    f"the other.\n"
    f"- **Equation (34) IS the null baseline.** It is the constant $n = 2.39$, "
    f"and the mean of the {len(b34)} measured slopes above $Re = 500$ is "
    f"{c34:.4f} — the same number to four decimals. There is nothing for it to "
    f"beat, and its {np.sqrt((b34.dev**2).mean()):.2f} % rms in that window is "
    f"the scatter of the six measurements, not a result.\n\n"
    f"So of the five printed expressions, **one (eq. 38) is clearly worth its "
    f"complexity on these data, one (eq. 33) is worth it modestly, one (eq. 34) "
    f"is a constant by construction, and one (eq. 37) is not resolved by the "
    f"data at all** — the next cells show exactly how unresolved."))'''))

# ------------------------------------- the eq. (37) window, three sensitivities
cells.append(md(r"""#### The eq. (37) window, and the three sensitivities that decide nothing

One branch, eq. (37) on $0.2 < Re < 1$, does not beat its own null. That gap is
0.02 percentage points on **five rows**, so before it is quoted as a finding it
gets three perturbations, each of which a referee would apply:

1. **Two of the five rows are Table I runs N and O**, whose printed $Re = 0.302$
   sits above the title of the table they are printed in ("Reynolds Number
   $< 0.2$"). §1 promised what eq. (33) would have given them; here it is, and
   the same rows re-placed with eq. (33) as the paper's own table placement
   implies.
2. **Leave one row out.**
3. **Penalise the null for the parameter it fits.** The constant is fitted
   in-sample on these five rows; eq. (37) has nothing fitted to them, so the
   comparison as it stands is generous to the null by one degree of freedom.""" ))

cells.append(code(r'''def rz_vs_null(g):
    """rms of the correlation and of a constant n fitted to the SAME rows."""
    c = g.n_meas.mean()
    return (float(np.sqrt((g.dev**2).mean())),
            float(np.sqrt((dev_pct(np.full(len(g), c), g.n_meas.values)**2).mean())), c)

b37 = keep[keep.branch == "eq37"].copy()
b37["n_eq33"] = 4.65 + 19.5 * b37.dD                       # eq. (33), for the two disputed rows
b37["dev33"]  = dev_pct(b37.n_eq33, b37.n_meas)
r37_rz, r37_null, c37 = rz_vs_null(b37)
print(f"THE {len(b37)} ROWS IN 0.2 < Re < 1, with what eq. (33) would give each")
print(b37[["source", "label", "Re", "dD", "n_meas", "n_pred", "dev", "n_eq33",
           "dev33"]].round(4).to_string(index=False))

NO   = b37[b37.flag.str.contains("Re_above_table_title")]
b37b = b37[~b37.flag.str.contains("Re_above_table_title")]
r3_rz, r3_null, c3 = rz_vs_null(b37b)

loo = []
for i in b37.index:
    g = b37.drop(i); rz_, nl_, _ = rz_vs_null(g)
    loo.append([b37.loc[i, "label"][:44], f"{rz_:.3f}", f"{nl_:.3f}",
                "RZ better" if rz_ < nl_ else "null better"])
n_flip = sum(1 for r in loo if r[3] == "RZ better")
pen = r37_null * np.sqrt(len(b37) / (len(b37) - 1))

print(f"\n(1) runs N and O are the {len(NO)} rows flagged Re_above_table_title; their printed Re is "
      f"{NO.Re.iloc[0]:.3g}.  Under eq. (37) they deviate by "
      f"{', '.join(f'{v:+.2f} %' for v in NO.dev)} - the two largest residuals in the window - and under "
      f"eq. (33) by {', '.join(f'{v:+.2f} %' for v in NO.dev33)}.")
print(f"\n(2) leave-one-out on the five rows:")
print(pd.DataFrame(loo, columns=["row dropped", "RZ rms %", "null rms %", "winner"]).to_string(index=False))

verdict = [["as placed by the printed Re", len(b37), f"{r37_rz:.3f}", f"{r37_null:.3f}",
            "null" if r37_null < r37_rz else "RZ"],
           ["N and O placed with eq. (33)", len(b37b), f"{r3_rz:.3f}", f"{r3_null:.3f}",
            "null" if r3_null < r3_rz else "RZ"],
           [f"null penalised x sqrt({len(b37)}/{len(b37)-1}) for its fitted constant", len(b37),
            f"{r37_rz:.3f}", f"{pen:.3f}", "null" if pen < r37_rz else "RZ"],
           [f"leave-one-out: RZ wins {n_flip} of {len(b37)} deletions", len(b37) - 1, "-", "-",
            "both" if 0 < n_flip < len(b37) else "consistent"]]
print()
print(pd.DataFrame(verdict, columns=["sensitivity", "rows", "RZ rms %", "null rms %",
                                     "which is lower"]).to_string(index=False))

metrics.update(
    eq37_rows                 = float(len(b37)),
    eq37_rz_rms_pct           = r37_rz,
    eq37_null_rms_pct         = r37_null,
    eq37_rz_rms_pct_no_NO     = r3_rz,
    eq37_null_rms_pct_no_NO   = r3_null,
    eq37_null_rms_pct_penalised = float(pen),
    eq37_loo_rz_wins          = float(n_flip),
)'''))

cells.append(code(r'''display(Markdown(
    f"**The eq. (37) window decides nothing, in either direction, and this page "
    f"does not let it.** As the rows are placed by their printed $Re$, the "
    f"correlation's rms is {r37_rz:.3f} % against the fitted constant's "
    f"{r37_null:.3f} % — a gap of {r37_rz - r37_null:+.3f} percentage points on "
    f"{len(b37)} rows. Every one of the three sensitivities above reverses it:\n\n"
    f"- placing runs N and O with **eq. (33)**, which is what the title of the "
    f"table they are printed in implies, leaves {len(b37b)} genuine rows and gives "
    f"RZ **{r3_rz:.3f} %** against the null's **{r3_null:.3f} %** — the "
    f"correlation now wins, and the two moved rows were the two largest residuals "
    f"in the window;\n"
    f"- **{n_flip} of {len(b37)}** leave-one-outs put the correlation ahead;\n"
    f"- charging the null the one degree of freedom it fits in-sample "
    f"({r37_null:.3f} $\\times\\sqrt{{{len(b37)}/{len(b37)-1}}}$ = {pen:.3f} %) "
    f"puts the correlation ahead.\n\n"
    f"So the honest statement is the weak one: **five rows spanning half a decade "
    f"of $Re$, two of them of disputed branch, cannot distinguish an exponent of "
    f"$-0.03$ from a constant** — not that eq. (37) is worse than a constant. An "
    f"earlier draft of this page said 'worse'; it survives none of the three "
    f"checks above and the claim is withdrawn."))'''))


cells.append(code(r'''# The rows the authors exclude - and the confound that makes their reason untestable.
tab = []
for nm, g in [("water / glycerol-water, d/D <= 0.04", keep[keep.dD <= 0.04]),
              ("water / glycerol-water, d/D >  0.04", keep[keep.dD >  0.04]),
              ("oil, d/D <= 0.04",                    oil[oil.dD <= 0.04]),
              ("oil, d/D >  0.04",                    oil[oil.dD >  0.04])]:
    tab.append([nm, len(g), f"{g.dev.abs().mean():.2f}", f"{g.dev.abs().max():.2f}"])
print(pd.DataFrame(tab, columns=["subset", "rows", "mean |dev| %", "max |dev| %"]).to_string(index=False))

win = A[(A.Re > 1) & (A.Re <= 200)]
hi  = win[win.dD > 0.04]
print(f"\nInside the eq. (38) window ({len(win)} rows), the rows above the stated d/D limit of 4e-2:")
print(hi[["label", "Re", "dD", "n_meas", "n_pred", "dev", "is_oil"]].to_string(index=False))
print(f"\nALL {len(hi)} of them are oil rows, and there are {len(win[(win.dD > 0.04) & (win.is_oil == 0)])} "
      f"non-oil rows above 4e-2 in that window.")

metrics.update(
    oil_highdD_max_pct   = oil[oil.dD > 0.04].dev.abs().max(),
    oil_lowdD_mean_pct   = oil[oil.dD <= 0.04].dev.abs().mean(),
    nonoil_highdD_in_eq38_window = float(len(win[(win.dD > 0.04) & (win.is_oil == 0)])),
)'''))

cells.append(code(r'''display(Markdown(
    f"The authors' own excluded rows behave exactly as they say: in the two oils "
    f"at $d/D > 0.04$ the correlation is out by a mean "
    f"**{oil[oil.dD > 0.04].dev.abs().mean():.1f} %** and a worst "
    f"**{oil[oil.dD > 0.04].dev.abs().max():.1f} %**, against "
    f"**{oil[oil.dD <= 0.04].dev.abs().mean():.2f} %** for the oil rows *below* "
    f"$d/D = 0.04$ and **{keep.dev.abs().mean():.2f} %** for everything not in "
    f"oil. Finding that is not a result — the paper predicts it.\n\n"
    f"**What is a result is that the paper's explanation cannot be tested on "
    f"this data set.** Two candidate causes are perfectly confounded. Inside the "
    f"eq. (38) window every one of the "
    f"{len(hi)} rows above the paper's own stated $d/D$ limit of "
    f"$4\\times10^{{-2}}$ is an oil row, and there are no non-oil rows there at "
    f"all. So 'the oil flocculated the spheres' and 'the correlation was never "
    f"claimed above $d/D = 0.04$' predict the same four failures, and nothing "
    f"in Tables I-VIII separates them. The paper asserts the first; this page "
    f"reports that the tables support neither over the other."))'''))

cells.append(code(r'''# Does the residual carry structure in a variable the correlation does not contain?
# The correlation is a function of Re and d/D ONLY - no fluid property, no density ratio.
w = keep.copy(); w["rho_ratio"] = w.rho_s / w.rho
gs = w.group.unique()
def cluster_boot(df, stat, nboot=2000):
    groups = df.group.unique()
    idx = {g: df.index[df.group == g] for g in groups}
    out = []
    for _ in range(nboot):
        pick = RNG.choice(groups, len(groups), replace=True)
        out.append(stat(df.loc[np.concatenate([idx[g] for g in pick])]))
    return np.percentile(out, [2.5, 97.5])

print(f"{len(w)} rows fall in {len(gs)} particle/liquid groups; rows within a group share a particle, a "
      f"liquid and a V0, so they are NOT independent.  Cluster bootstrap over groups, 2000 resamples.\n")
res = []
for v in ["rho_ratio", "mu", "rho", "rho_s"]:
    stat = lambda d, v=v: float(np.corrcoef(np.log10(d[v]), d.dev)[0, 1])
    lo, hi_ = cluster_boot(w, stat)
    res.append([f"log10 {v}", f"{stat(w):+.3f}", f"[{lo:+.3f}, {hi_:+.3f}]"])
for v, lab in [("Re", "log10 Re"), ("dD", "d/D")]:
    stat = lambda d, v=v: float(np.corrcoef(np.log10(d[v]) if v == "Re" else d[v], d.dev)[0, 1])
    lo, hi_ = cluster_boot(w, stat)
    res.append([lab, f"{stat(w):+.3f}", f"[{lo:+.3f}, {hi_:+.3f}]"])
print(pd.DataFrame(res, columns=["variable", "Pearson r with the residual", "cluster 95 % CI"]).to_string(index=False))

r_ratio = float(np.corrcoef(np.log10(w.rho_ratio), w.dev)[0, 1])
metrics.update(resid_corr_rho_ratio=r_ratio,
               resid_corr_mu=float(np.corrcoef(np.log10(w.mu), w.dev)[0, 1]))'''))

cells.append(code(r'''display(Markdown(
    f"Equations (33)-(39) contain **no fluid property and no density ratio** — "
    f"only $Re$ and $d/D$. The liquids behind these {len(w)} rows span "
    f"$\\mu$ from {w.mu.min():.2f} to {w.mu.max():.1f} cP and "
    f"$\\rho_s/\\rho$ from {w.rho_ratio.min():.2f} to {w.rho_ratio.max():.2f}, "
    f"so if a property were missing from the form the residual should know about "
    f"it. None of the six correlations above is distinguishable from zero once "
    f"the {len(gs)} particle/liquid clusters are respected.\n\n"
    f"**This is a null result and it is reported as one.** With {len(gs)} "
    f"clusters the confidence intervals are wide, so what the test excludes is a "
    f"*large* omitted-variable effect, not a small one. It does not show the "
    f"form is complete. Note also that the individual-row bootstrap would have "
    f"given intervals roughly half as wide and might have made one of these look "
    f"real; the clusters are the reason it does not."))'''))

# ------------------------------------------- section 3 : the two wall laws
cells.append(md(r"""### 3. The two wall laws, and the one test on this page that is genuinely two-sided

Equations (40) and (41) are the same statement with a different coefficient:

$$\log_{10} V_0 - \log_{10} V_i \;=\; \alpha \, \frac{d}{D},
\qquad \alpha = 0 \ \text{(eq. 40)},\qquad \alpha = 1 \ \text{(eq. 41)}.$$

Both sides of that are printed: the intercept column *is* $\log_{10} V_i$, and
$\log_{10} V_0$ and $d/D$ are their own columns. So $\alpha$ can be fitted, once
on the sedimentation tables and once on the fluidisation tables, and each fit is
a falsification test of the law the authors applied to the *other* data set.
**Neither law can be right on both sets** — that is the whole content of the two
equations — so this test has a real chance of embarrassing the paper.

The clusters matter here more than anywhere else on the page: a particle/liquid
group contributes two or three rows that share one $V_0$, so their residuals move
together. The interval below is a cluster bootstrap over groups, and the
individual-row interval is printed beside it.

**First, though: which rows can fail?** The quantity the fit sees is one number
per row, $y = \log_{10} V_0 - \log_{10} V_i$. If a table's intercept column
simply *is* its $\log_{10} V_0$ column, then $y \equiv 0$ on every row of it,
those rows vote for $\alpha = 0$ whatever the wall is doing, and they cannot
detect a wall effect of **any** size — while still counting as data and
narrowing the interval. That is exactly what Table II turns out to be, and it is
14 of the 36 sedimentation rows. The check that cannot fail, found inside the
page's strongest result. So the headline sedimentation fit is made on the rows
that can fail, the pooled fit is printed beside it, and the difference is stated
rather than averaged away.""" ))

cells.append(code(r'''def alpha_fit(df):
    x = df.d_over_D_printed.values
    y = (df.log10_V0_printed - df.intercept_log10_Vc_at_eps1).values
    return float((x * y).sum() / (x * x).sum())          # through the origin: alpha is the only parameter

def alpha_ci(df, gcol, nboot=2000):
    groups = df[gcol].unique(); idx = {g: df.index[df[gcol] == g] for g in groups}
    out = [alpha_fit(df.loc[np.concatenate([idx[g] for g in RNG.choice(groups, len(groups), replace=True)])])
           for _ in range(nboot)]
    return np.percentile(out, [2.5, 97.5])

sed_fit = sed[~sed.flag.str.contains("without_minus")].copy()
sed_fit["grp"] = sed_fit.particle + "|" + sed_fit.liquid + "|" + sed_fit.d_cm.astype(str)
sed_fit["y"] = sed_fit.log10_V0_printed - sed_fit.intercept_log10_Vc_at_eps1
flu_fit = flu.copy(); flu_fit["grp"] = flu_fit.particle + "|" + flu_fit.liquid
flu_fit["y"] = flu_fit.log10_V0_printed - flu_fit.intercept_log10_Vc_at_eps1
tI  = sed_fit[sed_fit.table == "I"]
tII = sed_fit[sed_fit.table == "II"]

# ---------- 0. CAN THE ROW FAIL?  y is the whole of what the fit sees.
print("CAN THE ROW FAIL?   y = log10 V0 - (intercept at eps = 1), one number per row")
inf = [[nm, len(d), f"{d.y.min():+.4f}", f"{d.y.max():+.4f}",
        f"{np.sqrt((d.y**2).mean()):.5f}", int((d.y.abs() <= 5e-4).sum())]
       for nm, d in [("sedimentation, Table I", tI), ("sedimentation, Table II", tII),
                     ("fluidisation, Tables III+IV", flu_fit)]]
print(pd.DataFrame(inf, columns=["subset", "rows", "min y", "max y", "rms y (decades)",
                                 "rows with |y| <= 0.0005"]).to_string(index=False))
print("\nTABLE II'S INTERCEPT COLUMN IS ITS log10 V0 COLUMN, group by group, as printed:")
tII_all = sed[sed.table == "II"]
for (p_, l_, d_), g in tII_all.groupby(["particle", "liquid", "d_cm"]):
    print(f"  {p_} {d_:g} cm in {l_[:34]:34s}  log10 V0 = {g.log10_V0_printed.iloc[0]:+.4f}"
          f"   intercepts {', '.join(f'{v:+.3f}' for v in g.intercept_log10_Vc_at_eps1)}"
          f"   at d/D {', '.join(f'{v:.2e}' for v in g.d_over_D_printed)}")
print(f"  -> {len(tII_all)} rows over {tII_all.groupby(['particle','liquid','d_cm']).ngroups} groups and "
      f"{tII_all.D_cm.nunique()} tube diameters; the intercept moves with the tube on NONE of them. The one "
      f"row that differs is run 12, whose 0.45 is the group's -0.457 with the minus sign lost - which is "
      f"itself evidence that the column is a copy, not a measurement.")

# ---------- 1. the fits.  The sedimentation headline is Table I: the rows that can fail.
a_tI  = alpha_fit(tI);  ci_tI  = alpha_ci(tI,  "grp")
a_tII = alpha_fit(tII); ci_tII = alpha_ci(tII, "grp")
a_sed = alpha_fit(sed_fit); ci_sed = alpha_ci(sed_fit, "grp")
a_flu = alpha_fit(flu_fit); ci_flu = alpha_ci(flu_fit, "grp")
flu_nooil = flu_fit[flu_fit.is_oil == 0]
a_fno = alpha_fit(flu_nooil); ci_fno = alpha_ci(flu_nooil, "grp")

out = [["sedimentation, TABLE I ONLY  <- headline", len(tI), tI.grp.nunique(), f"{a_tI:+.3f}",
        f"[{ci_tI[0]:+.3f}, {ci_tI[1]:+.3f}]"],
       ["sedimentation, Table II only (y == 0)", len(tII), tII.grp.nunique(), f"{a_tII:+.3f}",
        f"[{ci_tII[0]:+.3f}, {ci_tII[1]:+.3f}]"],
       ["sedimentation, Tables I+II pooled", len(sed_fit), sed_fit.grp.nunique(), f"{a_sed:+.3f}",
        f"[{ci_sed[0]:+.3f}, {ci_sed[1]:+.3f}]"],
       ["fluidisation, Tables III+IV", len(flu_fit), flu_fit.grp.nunique(), f"{a_flu:+.3f}",
        f"[{ci_flu[0]:+.3f}, {ci_flu[1]:+.3f}]"],
       ["fluidisation, non-oil only", len(flu_nooil), flu_nooil.grp.nunique(), f"{a_fno:+.3f}",
        f"[{ci_fno[0]:+.3f}, {ci_fno[1]:+.3f}]"]]
print("\n" + pd.DataFrame(out, columns=["data set", "rows", "groups", "fitted alpha",
                                        "cluster 95 % CI"]).to_string(index=False))
print(f"\nThe pooled interval is {np.ptp(ci_tI)/np.ptp(ci_sed):.1f}x NARROWER than the Table I interval, "
      f"on data half of which cannot move it.  The headline is the wide one.")

# individual-row bootstrap on the headline set, for the width comparison only.  SEEDED.
nb = [alpha_fit(tI.sample(len(tI), replace=True, random_state=np.random.default_rng(20260802 + k)))
      for k in range(2000)]
ci_naive = np.percentile(nb, [2.5, 97.5])
print(f"individual-row bootstrap on the Table I alpha: [{ci_naive[0]:+.3f}, {ci_naive[1]:+.3f}] - "
      f"the cluster interval is {np.ptp(ci_tI)/np.ptp(ci_naive):.2f}x wider.")

for nm, d in [("sedimentation Table I", tI), ("sedimentation Table II", tII),
              ("fluidisation", flu_fit)]:
    r40 = d.y
    r41 = r40 - d.d_over_D_printed
    print(f"\n{nm}: rms residual  with eq. (40) [alpha=0] {np.sqrt((r40**2).mean()):.5f} decades"
          f"   with eq. (41) [alpha=1] {np.sqrt((r41**2).mean()):.5f} decades")

metrics.update(wall_alpha_sed_tableI=a_tI, wall_alpha_sed_tableI_ci_lo=ci_tI[0],
               wall_alpha_sed_tableI_ci_hi=ci_tI[1],
               wall_alpha_sed_tableII=a_tII, wall_alpha_sed_tableII_ci_lo=ci_tII[0],
               wall_alpha_sed_tableII_ci_hi=ci_tII[1],
               wall_alpha_sed_pooled=a_sed, wall_alpha_sed_pooled_ci_lo=ci_sed[0],
               wall_alpha_sed_pooled_ci_hi=ci_sed[1],
               wall_alpha_sed_tableII_zero_rows=float((tII.y.abs() <= 5e-4).sum()),
               wall_alpha_fluidisation=a_flu, wall_alpha_flu_ci_lo=ci_flu[0], wall_alpha_flu_ci_hi=ci_flu[1],
               wall_alpha_flu_nonoil=a_fno)'''))

cells.append(code(r'''# ---------- 2. the two things that could have made the contrast an artefact.
# (a) run 12, the sign-corrupted cell.  It is a TABLE II row, so it does not enter
#     the headline fit at all; and no defensible handling of it moves anything.
sed_all = sed.copy(); sed_all["grp"] = sed_all.particle + "|" + sed_all.liquid + "|" + sed_all.d_cm.astype(str)
a_sed_all = alpha_fit(sed_all)
sed_sgn = sed_all.copy()
m12 = sed_sgn.flag.str.contains("without_minus")
sed_sgn.loc[m12, "intercept_log10_Vc_at_eps1"] = -sed_sgn.loc[m12, "intercept_log10_Vc_at_eps1"].abs()
a_sed_sgn = alpha_fit(sed_sgn)
print(f"RUN 12, the intercept printed without its minus sign, is a TABLE II row, so it is not in the "
      f"headline Table I fit ({a_tI:+.3f}) under any handling.  On the POOLED fit:")
print(f"  excluded (as above)                   alpha = {a_sed:+.3f}")
print(f"  kept exactly as printed, +0.45        alpha = {a_sed_all:+.3f}   (shift {a_sed_all - a_sed:+.3f})")
print(f"  minus sign restored, -0.45            alpha = {a_sed_sgn:+.3f}   (shift {a_sed_sgn - a_sed:+.3f})")
print(f"  -> the {a_sed_all - a_sed:+.2f} shift is what the cell AS PRINTED would do, and is a reason to "
      f"call it a printing error, not a fragility of the fit.")

# (b) the lever arm: the two data sets do not cover the same d/D.
print(f"\nd/D RANGE, because the fits are through the origin and the lever arm is not the same:")
print(f"  sedimentation d/D {sed_fit.d_over_D_printed.min():.2e} to {sed_fit.d_over_D_printed.max():.3f} "
      f"({(sed_fit.d_over_D_printed > 0.02).sum()} rows above 0.02)")
print(f"  fluidisation  d/D {flu_fit.d_over_D_printed.min():.2e} to {flu_fit.d_over_D_printed.max():.3f} "
      f"({(flu_fit.d_over_D_printed > 0.02).sum()} rows above 0.02)")
lo, hi_ = max(sed_fit.d_over_D_printed.min(), flu_fit.d_over_D_printed.min()), \
          min(sed_fit.d_over_D_printed.max(), flu_fit.d_over_D_printed.max())
s_ov  = sed_fit[(sed_fit.d_over_D_printed >= lo) & (sed_fit.d_over_D_printed <= hi_)]
s_ovI = s_ov[s_ov.table == "I"]
f_ov  = flu_fit[(flu_fit.d_over_D_printed >= lo) & (flu_fit.d_over_D_printed <= hi_)]
a_sov, ci_sov   = alpha_fit(s_ov),  alpha_ci(s_ov,  "grp")
a_sovI, ci_sovI = alpha_fit(s_ovI), alpha_ci(s_ovI, "grp")
a_fov, ci_fov   = alpha_fit(f_ov),  alpha_ci(f_ov,  "grp")
rows = [[f"sedimentation, Table I only", len(s_ovI), f"{a_sovI:+.3f}", f"[{ci_sovI[0]:+.3f}, {ci_sovI[1]:+.3f}]"],
        [f"sedimentation, pooled",       len(s_ov),  f"{a_sov:+.3f}",  f"[{ci_sov[0]:+.3f}, {ci_sov[1]:+.3f}]"],
        [f"fluidisation",                len(f_ov),  f"{a_fov:+.3f}",  f"[{ci_fov[0]:+.3f}, {ci_fov[1]:+.3f}]"]]
print(f"\nON THE COMMON d/D RANGE {lo:.2e} to {hi_:.3f}, WITH INTERVALS:")
print(pd.DataFrame(rows, columns=["data set", "rows", "fitted alpha", "cluster 95 % CI"]).to_string(index=False))

# (c) is alpha constant in d/D at all?  eq. (41) says it is.
SPLIT = 0.05
f_lo_, f_hi_ = flu_fit[flu_fit.d_over_D_printed < SPLIT], flu_fit[flu_fit.d_over_D_printed >= SPLIT]
a_flo, ci_flo = alpha_fit(f_lo_), alpha_ci(f_lo_, "grp")
a_fhi, ci_fhi = alpha_fit(f_hi_), alpha_ci(f_hi_, "grp")
print(f"\nIS THE FLUIDISATION alpha CONSTANT IN d/D, AS EQ. (41)'S LINEAR FORM REQUIRES?")
print(pd.DataFrame([[f"d/D < {SPLIT}", len(f_lo_), f"{a_flo:+.3f}", f"[{ci_flo[0]:+.3f}, {ci_flo[1]:+.3f}]"],
                    [f"d/D >= {SPLIT}", len(f_hi_), f"{a_fhi:+.3f}", f"[{ci_fhi[0]:+.3f}, {ci_fhi[1]:+.3f}]"]],
                   columns=["subset", "rows", "fitted alpha", "cluster 95 % CI"]).to_string(index=False))
print(f"  -> alpha falls by a factor {a_flo/a_fhi:.2f} across the range; a single linear term in d/D "
      f"cannot do that.")

metrics.update(wall_alpha_sed_with_run12=a_sed_all, wall_alpha_sed_run12_sign_restored=a_sed_sgn,
               wall_alpha_sed_common_range=a_sov, wall_alpha_sed_common_range_tableI=a_sovI,
               wall_alpha_sed_common_ci_lo=ci_sovI[0], wall_alpha_sed_common_ci_hi=ci_sovI[1],
               wall_alpha_flu_common_range=a_fov, wall_alpha_flu_common_ci_lo=ci_fov[0],
               wall_alpha_flu_common_ci_hi=ci_fov[1],
               wall_alpha_flu_dD_below_005=a_flo, wall_alpha_flu_dD_above_005=a_fhi)'''))

cells.append(code(r'''display(Markdown(
    f"**The tables pick the same two laws the paper does — on the rows that could "
    f"have said otherwise.**\n\n"
    f"On the {len(tI)} **Table I** rows $\\alpha = {a_tI:+.3f}$, 95 % CI "
    f"$[{ci_tI[0]:+.3f}, {ci_tI[1]:+.3f}]$ — consistent with eq. (40)'s zero and "
    f"**excluding eq. (41)'s $\\alpha = 1$**. On the fluidisation tables "
    f"$\\alpha = {a_flu:+.3f}$, CI $[{ci_flu[0]:+.3f}, {ci_flu[1]:+.3f}]$ — "
    f"consistent with eq. (41)'s one and **excluding zero**. Dropping the oil rows "
    f"moves it to {a_fno:+.3f}, CI $[{ci_fno[0]:+.3f}, {ci_fno[1]:+.3f}]$, so the "
    f"verdict is not an artefact of the rows the authors themselves distrust.\n\n"
    f"**Half the sedimentation data cannot fail, and that is a finding about the "
    f"paper, not an embarrassment for the page.** Table II's intercept column "
    f"reproduces its own $\\log_{{10}} V_0$ column on all "
    f"{len(tII_all)} of its rows: the {int((tII.y.abs() <= 5e-4).sum())} that "
    f"enter the fit agree to $|y| \\le 0.0005$ across two and three tube "
    f"diameters per group, and the one that does not is the one whose minus sign "
    f"is missing. The authors say as much on reprint page S95 (\"the "
    f"intercepts ... agree with those calculated for the terminal falling "
    f"velocity ... so that $V_i = V_0$\"), and whether they tabulated the "
    f"calculated value or measured that exactly cannot be told from the document. "
    f"Either way those rows carry no information about a wall effect: fitted "
    f"alone they give $\\alpha = {a_tII:+.3f}$, CI "
    f"$[{ci_tII[0]:+.3f}, {ci_tII[1]:+.3f}]$, an interval "
    f"{np.ptp(ci_tI)/np.ptp(ci_tII):.0f} times narrower than Table I's around a "
    f"number that is zero by construction. Pooling all 36 rows gives "
    f"{a_sed:+.3f} $[{ci_sed[0]:+.3f}, {ci_sed[1]:+.3f}]$ — the same verdict with "
    f"an interval {np.ptp(ci_tI)/np.ptp(ci_sed):.0f} times too narrow for the "
    f"evidence behind it. **The verdict survives the correction; the precision "
    f"does not, and the wide interval is the one quoted.**\n\n"
    f"**This is the strongest thing on the page, and it is still not independent "
    f"evidence.** Equation (41) was read off Fig. 21, and Fig. 21 was plotted "
    f"from these intercept columns; so the fluidisation $\\alpha$ is a refit, and "
    f"finding 1 is close to guaranteed. What is *not* guaranteed is the "
    f"sedimentation number: eq. (40) was justified on a different figure "
    f"(Fig. 20, $\\log V_0$ against $\\log V_i$, on which a wall term would show "
    f"as scatter rather than as a slope), and a wall effect of the fluidisation "
    f"size would have appeared here as $\\alpha \\approx 1$. It does not. **The "
    f"two data sets disagree about $\\alpha$ by more than either interval is "
    f"wide, and that disagreement is the paper's physical claim** — the wall "
    f"drags on the liquid in fluidisation and there is no through-flow to drag "
    f"on in sedimentation.\n\n"
    f"**The restricted range does more than qualify the contrast — it rejects "
    f"eq. (41)'s form.** The fits are through the origin, so the lever arm "
    f"matters, and the two data sets do not cover the same $d/D$: fluidisation "
    f"reaches {flu_fit.d_over_D_printed.max():.3f} where sedimentation stops at "
    f"{sed_fit.d_over_D_printed.max():.3f}. On the range they share the "
    f"fluidisation $\\alpha$ is {a_fov:+.3f}, CI "
    f"$[{ci_fov[0]:+.3f}, {ci_fov[1]:+.3f}]$ on {len(f_ov)} rows — an interval "
    f"that **excludes 1**, against sedimentation's {a_sovI:+.3f} "
    f"$[{ci_sovI[0]:+.3f}, {ci_sovI[1]:+.3f}]$ on the {len(s_ovI)} Table I rows "
    f"there. And $\\alpha$ is not constant in $d/D$ at all: "
    f"**{a_flo:+.3f}** below $d/D = 0.05$ ({len(f_lo_)} rows) against "
    f"**{a_fhi:+.3f}** above it ({len(f_hi_)} rows). A single linear term in "
    f"$d/D$, which is what eq. (41) is, cannot halve its own coefficient across "
    f"its own range. **So the full-range agreement with $\\alpha = 1$ is carried "
    f"by the large particles; on the small ones the wall term is more than twice "
    f"what eq. (41) writes.** The two-law verdict — a wall term in fluidisation, "
    f"none in sedimentation — is untouched by this. The *coefficient* 1 is not "
    f"supported outside the range that produced it.\n\n"
    f"And run 12, the intercept printed without its minus sign: it is a **Table "
    f"II** row, so it is not in the headline fit under any handling. On the "
    f"pooled fit, keeping it exactly as printed gives {a_sed_all:+.3f} and "
    f"restoring the minus sign gives {a_sed_sgn:+.3f}, against {a_sed:+.3f} with "
    f"it excluded. The {a_sed_all - a_sed:+.2f} is what the cell *as printed* "
    f"would do — which is a reason to call it a printing error, not a fragility "
    f"of the fit: no defensible handling of that cell moves $\\alpha$ by more "
    f"than {abs(a_sed_sgn - a_sed):.3f}."))'''))

cells.append(md(r"""#### The paired form, which removes $V_0$ entirely — and mostly cannot resolve anything

Tables III and IV run most particles in **both** columns. Subtracting the two
intercepts kills $\log V_0$:

$$\log_{10} V_i\big|_{2.44"} - \log_{10} V_i\big|_{1.5"}
 \;=\; \alpha\left(\frac{d}{D}\bigg|_{1.5"} - \frac{d}{D}\bigg|_{2.44"}\right).$$

That is a much cleaner statement — no terminal velocity, no Stokes law, nothing
but four printed numbers per particle. It is also, for most rows, **below the
resolution of the printed table**, and saying so is the point of including it.""" ))

cells.append(code(r'''pairs = []
for (t_, p_, l_), g in flu.groupby(["table", "particle", "liquid"]):
    if set(g.column_in) == {2.44, 1.5}:
        a = g[g.column_in == 2.44].iloc[0]; b = g[g.column_in == 1.5].iloc[0]
        pairs.append(dict(particle=p_, liquid=l_, is_oil=int(a.is_oil),
                          measured=a.intercept_log10_Vc_at_eps1 - b.intercept_log10_Vc_at_eps1,
                          predicted=b.d_over_D_printed - a.d_over_D_printed))
Pp = pd.DataFrame(pairs); Pp["ratio"] = Pp.measured / Pp.predicted

# The intercepts are printed to 2 or 3 decimals.  A difference of two such numbers
# carries a rounding uncertainty of order 0.01; below that the pair cannot vote.
QUANT = 0.01
resolvable = Pp[Pp.predicted > 2 * QUANT]
print(Pp.round(4).to_string(index=False))
print(f"\nthe intercept columns are printed to 2-3 decimals, so a DIFFERENCE of two of them is "
      f"uncertain by about +-{QUANT:.2f} decades.")
print(f"{len(resolvable)} of {len(Pp)} pairs have a predicted difference above twice that; "
      f"the other {len(Pp) - len(resolvable)} are predicted to differ by "
      f"{Pp[Pp.predicted <= 2*QUANT].predicted.min():.4f} to "
      f"{Pp[Pp.predicted <= 2*QUANT].predicted.max():.4f} decades and CANNOT TEST ANYTHING.")

slope_all = float((Pp.predicted * Pp.measured).sum() / (Pp.predicted**2).sum())
slope_res = float((resolvable.predicted * resolvable.measured).sum() / (resolvable.predicted**2).sum())
# SEEDED: pandas .sample() without random_state draws from numpy's GLOBAL RNG, which
# makes the interval move on every execution.  Each draw gets its own seeded Generator.
bs = [float((d.predicted * d.measured).sum() / (d.predicted**2).sum())
      for d in (resolvable.sample(len(resolvable), replace=True,
                                  random_state=np.random.default_rng(20260803 + k))
                for k in range(2000))]
ci = np.percentile(bs, [2.5, 97.5])
print(f"\nthrough-origin slope, all {len(Pp)} pairs           : {slope_all:.3f}")
print(f"through-origin slope, {len(resolvable)} resolvable pairs   : {slope_res:.3f}   95 % CI "
      f"[{ci[0]:.3f}, {ci[1]:.3f}]")
print(f"individual ratios on those pairs: {', '.join(f'{v:.2f}' for v in sorted(resolvable.ratio))}")

metrics.update(wall_paired_slope=slope_res, wall_paired_ci_lo=ci[0], wall_paired_ci_hi=ci[1],
               wall_paired_resolvable=float(len(resolvable)), wall_paired_total=float(len(Pp)))'''))

cells.append(code(r'''display(Markdown(
    f"The paired slope on the {len(resolvable)} pairs that can resolve it is "
    f"**{slope_res:.2f}**, 95 % CI $[{ci[0]:.2f}, {ci[1]:.2f}]$: consistent with "
    f"eq. (41)'s coefficient of 1, and excluding 0. Individual pairs run from "
    f"{resolvable.ratio.min():.2f} to {resolvable.ratio.max():.2f}, so what the "
    f"paired test supports is 'of order one', not 'one'.\n\n"
    f"**And {len(Pp) - len(resolvable)} of the {len(Pp)} pairs are pure noise for this "
    f"purpose.** Their two columns differ in $d/D$ by "
    f"{Pp[Pp.predicted <= 2*QUANT].predicted.max():.3f} decades or less, which is "
    f"at or below the rounding of the printed intercepts. If those pairs are "
    f"included the slope becomes {slope_all:.2f}, and it would be easy to quote "
    f"that as agreement — it is not, it is the large-$d/D$ pairs doing all the "
    f"work with four small-particle pairs contributing scatter. **A ratio "
    f"computed from a difference the table cannot represent is not evidence, "
    f"however close to 1 it lands.**"))'''))

# ------------------------------------- section 4 : internal identities
cells.append(md(r"""### 4. The identities the correlation has to satisfy

None of these needs a measurement. All of them can fail.""" ))

cells.append(code(r'''# (a) Tables VI and VII are related by an exact transformation.
#     The line of n against d/D at a given Re has slope s and reaches n = 0 at
#     d/D = -(d/D)_{n->0}, so its intercept is  n0 = s * |(d/D)_{n->0}|.
t = n0t.dropna(subset=["slope_n_vs_dD_printed", "minus_dD_at_n0_printed"]).copy()
t = t[t.slope_n_vs_dD_printed > 0]
t["n0_reconstructed"] = t.slope_n_vs_dD_printed * t.minus_dD_at_n0_printed
t["n0_reference"] = t.n0_printed
t.loc[t.Re == 0.2, "n0_reference"] = 4.65      # eq. (33)'s own intercept; Table VI has no Re = 0.2 row
t["dev"] = dev_pct(t.n0_reconstructed, t.n0_reference)
print("(a) TABLE VI x TABLE VII: n0 = slope x |(d/D)_{n->0}|")
print(t[["Re", "slope_n_vs_dD_printed", "minus_dD_at_n0_printed",
         "n0_reconstructed", "n0_reference", "dev"]].round(4).to_string(index=False))
print(f"    {len(t)} rows, worst {t.dev.abs().max():.3f} %, mean {t.dev.abs().mean():.3f} %.")
print(f"    The Re = 0.2 row is checked against eq. (33)'s intercept 4.65, which Table VI does not list, "
      f"and lands {t[t.Re == 0.2].dev.iloc[0]:+.3f} % away.")

mdd = n0t.minus_dD_at_n0_printed.dropna()
print(f"\n    the paper says |(d/D)_(n=0)| is 'approximately constant, with a mean value of -0.25'; "
      f"the {len(mdd)} tabulated values average {mdd.mean():.5f} "
      f"({dev_pct(mdd.mean(), 0.25):+.2f} % from the stated 0.25) and span "
      f"{mdd.min():.3f} to {mdd.max():.3f}.")

metrics.update(tableVI_VII_max_pct=t.dev.abs().max(), tableVI_VII_mean_pct=t.dev.abs().mean(),
               tableVII_mean_dD=mdd.mean())'''))

cells.append(code(r'''# (b) eqs (35) and (36) against the Table VI column they were fitted to.
v = n0t.dropna(subset=["n0_printed"]).copy()
v["n0_fit"] = np.where(v.Re < 1.0, 4.35 * v.Re**-0.03, 4.45 * v.Re**-0.1)
v["dev"] = dev_pct(v.n0_fit, v.n0_printed)
print("(b) eqs (35)/(36) vs the 18 tabulated n0 values (a FIT RESIDUAL)")
print(f"    mean |dev| {v.dev.abs().mean():.2f} %, worst {v.dev.abs().max():.2f} % at Re = "
      f"{v.loc[v.dev.abs().idxmax(), 'Re']:g}")
# the log columns of Table VI are a pure transcription identity - reported as such
lr = (v.log10_n0_printed - np.log10(v.n0_printed)).abs().max()
lre = (n0t.dropna(subset=['log10_Re_printed']).log10_Re_printed
       - np.log10(n0t.dropna(subset=['log10_Re_printed']).Re)).abs().max()
print(f"    Table VI's own log columns close to {lr:.2e} (log n0) and {lre:.2e} (log Re) on all 18 rows - "
      f"a TRANSCRIPTION check, structural, and it tests nothing physical.")

# (c) branch crossovers.
print("\n(c) BRANCH CROSSOVERS - the correlation is piecewise, so it need not be continuous")
def n_at(eq, Re, dD):
    return {"eq33": 4.65 + 19.5*dD, "eq37": (4.35 + 17.5*dD)*Re**-0.03,
            "eq38": (4.45 + 18.0*dD)*Re**-0.1, "eq39": 4.45*Re**-0.1, "eq34": 2.39}[eq]
cross = []
for Re, lo_eq, hi_eq in [(0.2, "eq33", "eq37"), (1.0, "eq37", "eq38"),
                         (200.0, "eq38", "eq39"), (500.0, "eq39", "eq34")]:
    row = [f"Re = {Re:g}", f"{lo_eq} -> {hi_eq}"]
    for dD in (0.0, 0.02, 0.04):
        row.append(f"{dev_pct(n_at(hi_eq, Re, dD), n_at(lo_eq, Re, dD)):+.3f} %")
    cross.append(row)
print(pd.DataFrame(cross, columns=["boundary", "branches", "d/D = 0", "d/D = 0.02",
                                   "d/D = 0.04"]).to_string(index=False))
jump500 = dev_pct(2.39, 4.45*500**-0.1)
jump200 = dev_pct(n_at("eq39", 200.0, 0.04), n_at("eq38", 200.0, 0.04))
metrics.update(crossover_Re500_pct=jump500, crossover_Re200_dD004_pct=jump200,
               crossover_Re1_pct=dev_pct(n_at("eq38", 1.0, 0.0), n_at("eq37", 1.0, 0.0)),
               crossover_Re02_pct=dev_pct(n_at("eq37", 0.2, 0.0), n_at("eq33", 0.2, 0.0)))'''))

cells.append(code(r'''# (d) the coefficients 17.5 and 18, which the paper derives on the line above.
print("(d) THE 17.5 AND THE 18")
for c, printed, eq in [(4.35, 17.5, "37"), (4.45, 18.0, "38")]:
    exact = c / 0.25
    at04  = dev_pct(c + printed*0.04, c + exact*0.04)
    print(f"    eq. ({eq}): the displayed line is {c}/(2.5e-1) = {exact:.3f}, collected as {printed} "
          f"({dev_pct(printed, exact):+.3f} %).  Effect on n at the stated d/D limit of 0.04: {at04:+.4f} %.")
print(f"    Using the ACTUAL mean of Table VII, {mdd.mean():.5f}, rather than the rounded 0.25:")
for c, printed in [(4.35, 17.5), (4.45, 18.0)]:
    print(f"      {c}/{mdd.mean():.5f} = {c/mdd.mean():.3f} against the printed {printed} "
          f"({dev_pct(printed, c/mdd.mean()):+.3f} %)")
metrics.update(coeff_175_vs_174_pct=dev_pct(17.5, 4.35/0.25),
               coeff_18_vs_178_pct=dev_pct(18.0, 4.45/0.25),
               coeff_rounding_effect_at_dD004_pct=dev_pct(4.45 + 18.0*0.04, 4.45 + 17.8*0.04))

# (e) HOW INDEPENDENT are eq. (39) and eq. (34)?  Ask what data each one sits on.
top = n0t[n0t.Re == n0t.Re.max()].iloc[0]
n0_top_fit = 4.45 * top.Re**-0.1
print("\n(e) WHERE eq. (39) AND eq. (34) MEET - AND WHAT THEY SHARE THERE")
print(f"    Table VI's LAST row is Re = {top.Re:g} with n0 printed {top.n0_printed:g} - "
      f"the same number eq. (34) states as the constant above Re = 500.")
print(f"    Table VII's row at the same Re prints the slope of n against d/D as "
      f"{top.slope_n_vs_dD_printed:g}, i.e. the d/D-independence eq. (39) asserts.")
print(f"    eq. (36)/(39) evaluated there: 4.45 x {top.Re:g}^-0.1 = {n0_top_fit:.4f}, "
      f"{dev_pct(n0_top_fit, top.n0_printed):+.2f} % from that same printed n0.")
print(f"    So the two expressions are not fitted to disjoint data at their boundary: they "
      f"share their top point, and meeting at Re = 500 to {jump500:+.3f} % is close to forced.")
metrics.update(eq39_at_tableVI_top=float(n0_top_fit),
               eq39_vs_tableVI_top_pct=float(dev_pct(n0_top_fit, top.n0_printed)),
               tableVI_top_Re=float(top.Re), tableVI_top_n0=float(top.n0_printed))'''))

cells.append(code(r'''display(Markdown(
    f"**A transcription check, not two independent fits meeting.** Equation (39) "
    f"is $4.45\\,Re^{{-0.1}}$, fitted on $1 < Re < 500$; equation (34) is the flat "
    f"$2.39$, stated separately from the runs above $Re = 500$. They meet at "
    f"$Re = 500$, where eq. (39) gives ${4.45*500**-0.1:.5f}$ against eq. (34)'s "
    f"$2.39$ — **{jump500:+.3f} %**, one part in seven thousand. An earlier draft "
    f"of this page called that an identity that could have failed. It could not, "
    f"or not by much: **Table VI's last row prints $n_0 = {top.n0_printed:g}$ at "
    f"$Re = {top.Re:g}$** — eq. (34)'s constant *is* that printed value — "
    f"eq. (36)/(39) already passes within {abs(dev_pct(n0_top_fit, top.n0_printed)):.2f} % "
    f"of it, and Table VII's row at the same $Re$ prints the $d/D$ slope as "
    f"{top.slope_n_vs_dD_printed:g}, which is eq. (39)'s dropped wall term. The "
    f"whole high-$Re$ end — those two table rows, eq. (36), eq. (39) and eq. (34) "
    f"— is one exercise on Figs. 18 and 19. What the agreement **does** test is "
    f"the reading: a mis-read $4{{\\cdot}}45$, $0{{\\cdot}}1$ or $2{{\\cdot}}39$ "
    f"would break it at once, and that is the whole of its value.\n\n"
    f"**The boundary that is genuinely discontinuous.** At $Re = 200$ the wall "
    f"term is simply dropped: eq. (39) has no $d/D$ in it. At $d/D = 0$ the two "
    f"branches match exactly; at the paper's own $d/D$ limit of 0.04 the index "
    f"jumps by **{jump200:+.1f} %** as $Re$ crosses 200. Anyone stepping $Re$ "
    f"through 200 in a solver will see it. The other two boundaries jump by "
    f"{metrics['crossover_Re02_pct']:+.2f} % (at $Re=0.2$) and "
    f"{metrics['crossover_Re1_pct']:+.2f} % (at $Re=1$) — small, but nonzero, and "
    f"the correlation is not continuous anywhere except at $Re = 500$.\n\n"
    f"**The 17.5 and the 18 are not what the line above them says.** The paper "
    f"displays eq. (37) as $\\frac{{4.35}}{{2.5\\times10^{{-1}}}}Re^{{-0.03}}"
    f"\\frac{{d}}{{D}} + 4.35\\,Re^{{-0.03}}$ and then collects it as "
    f"$(4.35 + 17.5\\,d/D)Re^{{-0.03}}$, but $4.35/0.25 = 17.400$. The same for "
    f"eq. (38): $4.45/0.25 = 17.800$, collected as 18. Both are rounded **up**, "
    f"by {dev_pct(17.5, 17.4):+.2f} % and {dev_pct(18.0, 17.8):+.2f} %. It costs "
    f"nothing — at the stated $d/D$ limit of 0.04 the two versions of eq. (38) "
    f"differ by {metrics['coeff_rounding_effect_at_dD004_pct']:+.3f} % in $n$ — "
    f"but it is worth knowing before anyone re-derives the coefficients and "
    f"finds 17.4."))'''))

# --------------------------------------- section 5 : non-spherical particles
cells.append(md(r"""### 5. Non-spherical particles: eq. (42), with the null baseline that nearly matches it

Table VIII gives $n$, $K$ and $S$ for seven particles, all in fully turbulent
flow ($Re$ from 2650 to 7150), where eq. (34) says $n$ should be a constant
2.39 for spheres. Equation (42), $n = 2.7\,K^{0.16}$, is the shape correction.

Before believing it: six of the seven particles have $K$ between 0.517 and 0.71
and $n$ between 2.38 and 2.55. **One point — the plates — carries the whole
correlation**, and the page says so with numbers.""" ))

cells.append(code(r'''x = nsp.copy()
x["n_eq42"]  = 2.7 * x.K_printed**0.16
x["dev"]     = dev_pct(x.n_eq42, x.n_table8)
x["K_recon"] = (np.pi / 6.0) * (x.d_s_cm / x.d_p_cm)**3          # the paper's own definition of K
x["K_dev"]   = dev_pct(x.K_recon, x.K_printed)
print(x[["particle", "Re_printed", "K_printed", "n_table8", "n_eq42", "dev"]].round(4).to_string(index=False))

c_ns  = x.n_table8.mean()
rms42 = float(np.sqrt((x.dev**2).mean()))
rmsc  = float(np.sqrt((dev_pct(np.full(len(x), c_ns), x.n_table8.values)**2).mean()))
y = x[~x.particle.str.contains("plates")]
c6    = y.n_table8.mean()
rms42_6 = float(np.sqrt((y.dev**2).mean()))
rmsc_6  = float(np.sqrt((dev_pct(np.full(len(y), c6), y.n_table8.values)**2).mean()))
print(f"\n            eq. (42) rms   constant-n null rms   null constant")
print(f"all 7        {rms42:8.2f} %        {rmsc:8.2f} %          {c_ns:.4f}")
print(f"drop plates  {rms42_6:8.2f} %        {rmsc_6:8.2f} %          {c6:.4f}")

n_sphere = 2.7 * (np.pi / 6.0)**0.16
print(f"\neq. (42) at the sphere's own K = pi/6 = {np.pi/6:.4f} gives n = {n_sphere:.4f}, against "
      f"eq. (34)'s 2.39: {dev_pct(n_sphere, 2.39):+.2f} %.")
print(f"Table VIII's K for the ball bearings is printed {x.K_printed.iloc[0]}, and pi/6 = {np.pi/6:.4f} "
      f"({dev_pct(np.pi/6, x.K_printed.iloc[0]):+.3f} %) - which is what fixes the definition of K.")
print(f"\nTable VIII's n column vs the 2.44 in. column slope of Tables V and III: "
      f"max |difference| = {(x.n_table8 - x.slope_244).abs().max():.3f}")
print(f"Table VIII's n column vs the 1.5 in. column slope:              "
      f"max |difference| = {(x.n_table8 - x.slope_15).abs().max():.3f}")

metrics.update(eq42_rms_pct=rms42, eq42_null_rms_pct=rmsc,
               eq42_rms_pct_no_plates=rms42_6, eq42_null_rms_pct_no_plates=rmsc_6,
               eq42_at_sphere_K_vs_eq34_pct=dev_pct(n_sphere, 2.39),
               tableVIII_n_equals_244_slope_maxdiff=(x.n_table8 - x.slope_244).abs().max())'''))

cells.append(code(r'''print("K RECOMPUTED FROM THE PAPER'S OWN PRINTED d_s, d_p AND DEFINITION")
print(x[["particle", "d_s_cm", "d_p_cm", "K_printed", "K_recon", "K_dev"]].round(4).to_string(index=False))
ok = x[x.K_dev.abs() < 0.1]; bad = x[x.K_dev.abs() >= 0.1]
print(f"\n{len(ok)} of the {len(x)} rows reproduce to better than 0.1 %: "
      f"{', '.join(ok.particle)}.")
print(f"{len(bad)} do not (three distinct shapes; the two 1/4 in. cylinder rows share d_s, d_p and K):")
for _, r in bad.iterrows():
    print(f"    {r.particle:24s} K printed {r.K_printed:.4f}, recomputed {r.K_recon:.4f}  ({r.K_dev:+.1f} %)")
metrics.update(K_identity_max_pct=x.K_dev.abs().max(), K_identity_rows_ok=float(len(ok)),
               K_identity_rows_bad=float(len(bad)))'''))

cells.append(code(r'''display(Markdown(
    f"**Equation (42) beats a constant, but only because of one particle.** Over "
    f"all seven the correlation's rms is **{rms42:.2f} %** against a constant "
    f"$n = {c_ns:.2f}$'s **{rmsc:.2f} %** — a factor "
    f"{rmsc/rms42:.1f}. Drop the steel plates, the only particle with $K$ far "
    f"from the rest ($K = {x.K_printed.min():.4f}$ against "
    f"{y.K_printed.min():.3f}-{y.K_printed.max():.3f} for the other six), and "
    f"the margin collapses to **{rms42_6:.2f} %** against **{rmsc_6:.2f} %**. "
    f"*The exponent 0.16 is fitted to seven points of which six are nearly "
    f"coincident in $K$.* That is the honest description of eq. (42) and it is "
    f"not what its ubiquity in textbooks suggests.\n\n"
    f"**A cross-check that closes.** Table VIII's $n$ column is, to the last "
    f"digit on all seven rows (max difference "
    f"{(x.n_table8 - x.slope_244).abs().max():.3f}), the **2.44 in. column** "
    f"slope from Tables V and III — and *not* the 1.5 in. slope, which differs "
    f"by up to {(x.n_table8 - x.slope_15).abs().max():.2f}. Equation (42) is "
    f"therefore a correlation of the large-column slopes only, at $d_p/D$ around "
    f"{nsp.dp_over_D_244.dropna().mean():.2f}, and the wall term of eq. (41) "
    f"has not been removed from it.\n\n"
    f"**And the abscissa it is fitted against is internally inconsistent on "
    f"{len(bad)} of the {len(x)} rows.** $K = (\\pi/6)(d_s/d_p)^3$ is the "
    f"paper's own printed definition and $d_s$, $d_p$ are its own printed "
    f"columns, so $K$ can be recomputed. Sphere, cubes and plates come back to "
    f"better than 0.1 %. The two 1/4 in. cylinder rows — which share $d_s$, "
    f"$d_p$ and $K$, so they are one shape counted twice — are out by "
    f"{x[x.particle.str.contains('cylinders')].K_dev.iloc[0]:+.1f} %, the "
    f"hexagonal prisms by {x[x.particle.str.contains('hexagonal')].K_dev.iloc[0]:+.1f} %, "
    f"and the 1/8 x 1/4 in. cylinder by "
    f"{x[x.particle == 'Steel cylinder'].K_dev.iloc[0]:+.1f} %. Note that Table "
    f"V's $d_s$ column is the one column of the whole paper that is *illegible* "
    f"at 600 dpi in this scan, so the page attributes the discrepancy to no "
    f"particular cell and draws no conclusion beyond this: **eq. (42)'s $x$-axis "
    f"cannot be reconstructed from the paper for {len(bad)} of its {len(x)} "
    f"points, three shapes out of six.**"))'''))

# ------------------------------------------------ pymrm implementation
cells.append(md(r"""## PyMRM implementation

The correlation itself needs no solver — it is five algebraic expressions, and
sections 1-5 above would run with pymrm uninstalled. What needs a solver is the
object the correlation becomes when it is used: **a flux function**.

A closed column of suspension, depth $x$ measured downwards from the free
surface, solids volume fraction $\phi = 1-\varepsilon$. A solids volume balance
over a slice is

$$\frac{\partial \phi}{\partial t} + \frac{\partial F}{\partial x} = 0,
\qquad F(\phi) \;=\; \phi\,V_0\,(1-\phi)^{\,n},$$

with $F = 0$ at both ends — no solids cross the free surface, none cross the
base. $F$ is eq. (28) written as a flux: the batch settling rate $V_c$ of a
suspension of voidage $\varepsilon$, times the solids fraction.

**This balance is not in Richardson and Zaki.** They measure the constant rate of
fall of the sludge line and say so; they never write a conservation law and never
model what happens when the settled bed grows up to meet the interface. The
balance is elementary and is derived here rather than cited. What the paper
contributes is $F$, and only $F$.

The discretisation: one constant `construct_div` with `nu = 0` — a column of
constant cross-section is Cartesian — built once, outside the time loop. $F$ is
nonlinear and non-monotone ($F$ rises to a maximum at $\phi^\star = 1/(n+1)$ and
falls back to zero), so a plain sign-of-velocity upwind is *wrong* above
$\phi^\star$; the face flux is the exact Godunov flux for a scalar law with one
interior maximum. Explicit stepping at CFL 0.4 — no Newton, no Jacobian, so none
of the `NumJac` shape traps arise.""" ))

cells.append(code(r'''def rz_flux(phi, V0, n):
    p = np.clip(phi, 0.0, 1.0)
    return V0 * p * (1.0 - p)**n                      # eq. (28) as a solids flux

def rz_flux_prime(phi, V0, n):
    p = np.clip(phi, 0.0, 1.0)
    return V0 * (1.0 - p)**(n - 1.0) * (1.0 - (n + 1.0) * p)

def settle(nx, t_end, phi0=0.10, V0=1.0, n=4.65, H=1.0, nu=0, sign=+1.0,
           flip_upwind=False, cfl=0.4):
    """Batch settling of a uniform suspension.  Returns cell centres and phi(t_end).

    nu = 0: Cartesian.  A settling column has constant cross-section, so the
    divergence operator carries no area profile.  Both boundaries are zero-flux,
    imposed directly on the face-flux vector rather than through a bc dict,
    because the flux is a nonlinear function of phi and not a linear form in it.
    """
    phi_star = 1.0 / (n + 1.0)
    x_f = np.linspace(0.0, H, nx + 1)
    x_c = 0.5 * (x_f[:-1] + x_f[1:])
    div = construct_div((nx, 1), x_f, nu=nu)          # CONSTANT: built once, outside the loop
    phi = np.full((nx, 1), float(phi0))
    t = 0.0; dx = H / nx
    while t < t_end:
        speed = np.nanmax(np.abs(rz_flux_prime(phi, V0, n)))
        if not np.isfinite(speed) or speed <= 0:
            return x_c, phi[:, 0], True
        dt = min(cfl * dx / speed, t_end - t)
        pl, pr = phi[:-1, 0], phi[1:, 0]
        if flip_upwind:
            F_int = rz_flux(pr, V0, n)                # deliberately wrong upwind direction
        else:
            lo, hi = np.minimum(pl, pr), np.maximum(pl, pr)
            f_lo, f_hi = rz_flux(lo, V0, n), rz_flux(hi, V0, n)
            F_min = np.minimum(f_lo, f_hi)
            F_max = np.maximum(f_lo, f_hi)
            F_max = np.where((lo <= phi_star) & (phi_star <= hi),
                             np.maximum(F_max, rz_flux(phi_star, V0, n)), F_max)
            F_int = np.where(pl <= pr, F_min, F_max)  # exact Godunov flux
        F = np.concatenate(([0.0], F_int, [0.0])).reshape(nx + 1, 1)
        phi = phi - dt * sign * (div @ F)
        t += dt
        if not np.all(np.isfinite(phi)):
            return x_c, phi[:, 0], True
    return x_c, phi[:, 0], False

PHI0, V0_REF, N_REF, H_REF, T_END = 0.10, 1.0, 4.65, 1.0, 0.35
u_rz = V0_REF * (1.0 - PHI0)**N_REF
print(f"phi0 = {PHI0}, n = {N_REF} (eq. 33 at d/D -> 0), V0 = {V0_REF}")
print(f"Richardson-Zaki settling velocity  V0 (1-phi0)^n            = {u_rz:.8f}")
print(f"Rankine-Hugoniot speed of the phi=0 / phi0 shock  F(phi0)/phi0 = "
      f"{rz_flux(PHI0, V0_REF, N_REF)/PHI0:.8f}")
print("\nThose two are the SAME NUMBER by construction: F(phi) = phi V0 (1-phi)^n, so")
print("F(phi0)/phi0 is V0(1-phi0)^n identically.  The identity is algebraic and is NOT a result.")'''))

# ------------------------------------------------------------------- results
cells.append(md(r"""## Results""" ))

cells.append(code(r'''x_c, phi, _ = settle(400, T_END)
xs = u_rz * T_END

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for nx, ls in [(50, ":"), (100, "--"), (400, "-")]:
    xc, ph, _ = settle(nx, T_END)
    ax[0].plot(1.0 - ph, xc, ls, lw=1.6, label=f"$n_x$ = {nx}")
ax[0].axhline(xs, color="k", lw=0.8)
ax[0].annotate(f"interface, $V_0\\varepsilon_0^{{\\,n}}t$ = {xs:.4f}", (0.35, xs),
               textcoords="offset points", xytext=(4, 6), fontsize=8)
ax[0].invert_yaxis(); ax[0].set_xlabel(r"voidage $\varepsilon$"); ax[0].set_ylabel("depth below the surface")
ax[0].set_title(f"batch settling at $t$ = {T_END}"); ax[0].legend(fontsize=8)

Re_grid = np.logspace(-4, 4, 800)
for dD, c in [(0.0, "k"), (0.02, "C0"), (0.04, "C1")]:
    ax[1].plot(Re_grid, n_richardson_zaki(Re_grid, np.full_like(Re_grid, dD)), c,
               lw=1.5, label=f"$d/D$ = {dD}")
ax[1].plot(A.Re, A.n_meas, "o", ms=4, mfc="none", color="C3", label="measured, Tables I-IV")
for r in (0.2, 1.0, 200.0, 500.0):
    ax[1].axvline(r, color="0.7", lw=0.7)
ax[1].set_xscale("log"); ax[1].set_xlabel("$Re = V_0 d\\rho/\\mu$"); ax[1].set_ylabel("$n$")
ax[1].set_title("eqs. (33), (37), (38), (39), (34)"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.show()
print(cite_data(load_meta("rz1954_fluidisation_spheres.csv", page=PAGE)))
print(CITE_NOTE)'''))

cells.append(code(r'''# What a bed of a real particle from Table III does, run through the paper's own recipe.
row = flu[(flu.particle == "Ballotini") & (flu.d_cm == 5.1e-2) & (flu.column_in == 2.44)].iloc[0]
n_row = n_richardson_zaki(row.Re_printed, row.d_over_D_printed)
V_i = 10**(np.log10(row.V0_cm_s) - row.d_over_D_printed)     # eq. (41)
eps = np.linspace(0.4, 1.0, 200)
print(f"Table III, {row.particle} d = {row.d_cm} cm in the {row.column_in} in. column")
print(f"  printed  Re = {row.Re_printed:g},  d/D = {row.d_over_D_printed:g},  V0 = {row.V0_cm_s} cm/s")
print(f"  eq. (38) n  = {n_row:.4f}   against the measured slope {row.slope_n} "
      f"({dev_pct(n_row, row.slope_n):+.2f} %)")
print(f"  eq. (41) V_i = {V_i:.4f} cm/s against the printed intercept "
      f"10^{row.intercept_log10_Vc_at_eps1} = {10**row.intercept_log10_Vc_at_eps1:.4f} cm/s "
      f"({dev_pct(V_i, 10**row.intercept_log10_Vc_at_eps1):+.2f} %)")
print(f"  the velocity to hold this bed at eps = 0.6:  V_c = V_i eps^n = "
      f"{V_i*0.6**n_row:.4f} cm/s   (measured n gives {10**row.intercept_log10_Vc_at_eps1*0.6**row.slope_n:.4f})")
metrics.update(worked_n_dev_pct=dev_pct(n_row, row.slope_n),
               worked_Vi_dev_pct=dev_pct(V_i, 10**row.intercept_log10_Vc_at_eps1))'''))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

### The pymrm solve: what it tests, and what it is blind to

The interface in a batch settler is a shock between $\phi = 0$ above and
$\phi_0$ below. Its speed follows from Rankine-Hugoniot, $F(\phi_0)/\phi_0$,
which for this flux is $V_0(1-\phi_0)^n$ **identically** — so comparing the two
is algebra, not evidence, and the page says that where it prints the number.

The check that can fail is grid convergence of the whole profile against the
exact shock solution, measured on a window well above the base so that the
growing sediment (which the correlation was never measured for) does not enter.""" ))

cells.append(code(r'''WINDOW = 0.6 * H_REF     # the L1 window: contains the interface, excludes the basal sediment
def L1_error(x_c, phi, u_ref, phi0=PHI0):
    m = x_c < WINDOW
    exact = np.where(x_c[m] < u_ref * T_END, 0.0, phi0)
    return float(np.abs(phi[m] - exact).sum() * (x_c[1] - x_c[0]))

print(f"L1 error against the exact shock solution on 0 < x < {WINDOW:g}")
prev, orders = None, []
for nx in [50, 100, 200, 400, 800, 1600]:
    xc, ph, _ = settle(nx, T_END)
    e = L1_error(xc, ph, u_rz)
    o = np.nan if prev is None else np.log2(prev / e)
    if prev is not None: orders.append(o)
    print(f"  nx = {nx:5d}   L1 = {e:.4e}" + ("" if prev is None else f"   order {o:.3f}"))
    prev = e
L1_400 = L1_error(*settle(400, T_END)[:2], u_rz)
print(f"\nmean observed order {np.mean(orders):.3f} - first-order upwind on a shock, which is what it is.")

xc, ph, _ = settle(400, T_END)
solids_res = abs(ph.sum() * (H_REF/400) - PHI0*H_REF) / (PHI0*H_REF)
print(f"solids balance residual at nx = 400: {solids_res:.2e}")
metrics.update(pymrm_L1_nx400=L1_400, pymrm_order_mean=float(np.mean(orders)),
               pymrm_solids_residual=solids_res)'''))

cells.append(md(r"""### Break table: inject a defect, confirm the number moves

Two checks are on trial here. `L1` is the profile error against the exact shock;
`solids residual` is the global mass balance. The point of the table is that the
second one **cannot fail** for almost anything, and would have looked like the
page's cleanest result if quoted alone.""" ))

cells.append(code(r'''def break_row(name, ref_u=None, nx=400, **kw):
    xc, ph, diverged = settle(nx, T_END, **kw)
    if diverged or not np.all(np.isfinite(ph)):
        return [name, "DIVERGED", "-", "-"]
    u = u_rz if ref_u is None else ref_u
    p0 = kw.get("phi0", PHI0)
    e = L1_error(xc, ph, u, p0)
    m = abs(np.nansum(ph) * (H_REF/nx) - p0*H_REF) / (p0*H_REF)
    return [name, f"{e:.3e}", f"x {e/L1_400:.1f}", f"{m:.1e}"]

brk = [break_row("none (reference)"),
       break_row("nu = 1: cylindrical instead of Cartesian", nu=1),
       break_row("flux sign flipped in the update", sign=-1.0),
       break_row("upwind direction reversed", flip_upwind=True),
       break_row("CFL = 1.6 (above the limit; phi is clipped, so this is a bad answer, not a blow-up)",
                 cfl=1.6),
       break_row("nx = 25 (far too coarse)", nx=25),
       break_row("n = 3.00 instead of 4.65", n=3.0),
       break_row("n = 3.00, judged against ITS OWN shock speed",
                 n=3.0, ref_u=V0_REF*(1.0-PHI0)**3.0),
       break_row("V0 1 % high, judged against V0 = 1", V0=1.01)]
print(pd.DataFrame(brk, columns=["injected defect", "L1 vs exact", "vs reference",
                                 "solids residual"]).to_string(index=False))
metrics.update(brk_nu1_ratio=L1_error(*settle(400, T_END, nu=1)[:2], u_rz)/L1_400,
               brk_upwind_ratio=L1_error(*settle(400, T_END, flip_upwind=True)[:2], u_rz)/L1_400,
               brk_n3_own_speed_ratio=L1_error(*settle(400, T_END, n=3.0)[:2],
                                              V0_REF*(1.0-PHI0)**3.0)/L1_400)'''))

cells.append(code(r'''display(Markdown(
    f"**The solids balance is a check that cannot fail.** It returns exactly "
    f"zero for a wrong exponent, a wrong $V_0$, an unstable time step and a "
    f"25-cell grid, and $10^{{-15}}$ for a flipped flux sign and a reversed "
    f"upwind direction. That is structural: `construct_div` telescopes, and both "
    f"boundary fluxes are set to zero, so the total is conserved whatever "
    f"nonsense is put into the interior faces. The **only** defect it catches is "
    f"`nu = 1`, and it catches that because a cylindrical divergence does not "
    f"telescope on a Cartesian grid. Quoting the balance residual as evidence "
    f"that the settling model is right would be exactly the defect the verifier "
    f"brief exists to catch.\n\n"
    f"**The L1 error is a real check on the discretisation and on nothing else.** "
    f"It moves by {metrics['brk_nu1_ratio']:.0f}x for `nu = 1`, "
    f"{metrics['brk_upwind_ratio']:.0f}x for a reversed upwind direction and "
    f"three orders of magnitude for a flipped sign. But set $n = 3.00$ and judge "
    f"the result against the shock speed *that value of $n$ implies* and the "
    f"error moves by only {metrics['brk_n3_own_speed_ratio']:.1f}x — **the solve "
    f"has no opinion about the value of $n$ at all.** A 1 % error in $V_0$ is "
    f"likewise almost invisible, because a 1 % shift of the interface is "
    f"comparable to the width the first-order scheme smears it over. Nothing in "
    f"this section is evidence for or against Richardson and Zaki's correlation; "
    f"it is evidence that the flux was integrated correctly."))'''))

cells.append(md(r"""### Every agreement number on this page, and whether it can fail

| # | what is compared | value | can it fail? |
|---|---|---|---|
| 1 | four printed identities against the transcription | see §1 | **yes** — three of them do, on five named cells |
| 2 | Stokes' law against Table I, and against Table II | §1 | **yes, and it must fail on Table II** — it does, by up to two orders |
| 3 | eqs. (33)-(39) against 58 non-oil measured slopes | §2 | yes, but it is a **fit residual** — these are the data the fit was made on |
| 4 | constant-$n$ null baseline per branch | §2 | it *is* the baseline; eq. (34) coincides with it |
| 4a | eq. (37) against its null on 5 rows | §2 | **no — it decides nothing.** The 0.02-point gap reverses under all three sensitivities printed there |
| 5 | residual against $\mu$, $\rho$, $\rho_s/\rho$ | §2 | **yes** — a null result, cluster-bootstrapped |
| 6 | wall coefficient $\alpha$, sedimentation vs fluidisation | §3 | **yes, on Table I and on Tables III+IV** — each law is falsifiable on the other's data |
| 6a | the 14 Table II rows inside that fit | §3 | **NO — their residual is identically zero**, because the intercept column reproduces $\log_{10}V_0$. Excluded from the headline; the pooled fit is printed beside it |
| 7 | paired wall test across the two columns | §3 | **yes, on 6 of 10 pairs**; the other 4 cannot resolve it |
| 8 | Table VI $\times$ Table VII transformation | §4 | **yes** — any transcription slip breaks it |
| 9 | eq. (39) at $Re = 500$ against eq. (34) | §4 | **as a transcription check, yes; as evidence of two independent fits, no** — Table VI's last row *is* eq. (34)'s 2.39 |
| 10 | eqs. (35)/(36) against Table VI | §4 | fit residual |
| 11 | Table VI's own log columns | §4 | **no** — a pure transcription identity, labelled as one |
| 12 | eq. (42) against Table VIII, with a null | §5 | **yes**, and the null nearly matches it |
| 13 | $K$ recomputed from printed $d_s$, $d_p$ | §5 | **yes** — it fails on 4 of the 7 rows (3 distinct shapes) |
| 14 | pymrm $L1$ against the exact shock | §Validation | **yes** — order 1, break table |
| 15 | pymrm solids balance | §Validation | **no** — structural; the break table proves it |
| 16 | RH speed against $V_0\varepsilon_0^n$ | §Validation | **no** — algebraically identical, and said so |""" ))

cells.append(code(r'''report_agreement("A1.5", metrics)'''))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**To the correlation: nothing.** Equations (33)-(42) are algebra and need no
solver, and this page does not pretend otherwise. What the reimplementation adds
is five things, none of which is a differential equation.

**The tables, machine-readable, with the transcription certified by identities
the paper itself prints.** Richardson and Zaki is cited thousands of times a year
and what is cited is one line of it. The 66 measured exponents, the 18 $n_0$
values and the two intercept columns have been unavailable in any usable form;
they are now four CSVs, and the four printed identities of §1 say how far they
can be trusted — including the five places the printing is wrong.

**The correlation costed against a null model, window by window.** The
whole-table agreement of 1.2 % is not the useful number. Of the five printed
expressions, eq. (38) beats a constant by a factor of seven and a half on 23
rows, eq. (33) beats it modestly, eq. (34) *is* a constant, and **eq. (37) is not
distinguishable from a constant in either direction on the five rows that fall in
its window** — the 0.02-percentage-point gap between them reverses if the two
rows of disputed branch are moved, under 2 of 5 leave-one-outs, and if the
in-sample null is charged for the parameter it fits. That last statement is not
in the paper and it is the one a reader needs before trusting the
$0.2 < Re < 1$ branch.

**The two wall laws decided on the authors' own numbers — on the rows that could
have decided otherwise.** Fitting the single coefficient $\alpha$ in
$\log V_0 - \log V_i = \alpha(d/D)$ separately to the sedimentation and
fluidisation tables gives intervals that exclude each other and that pick, in
each case, the equation the authors applied there. This is the page's strongest
result and its limits are stated: **Table II's intercept column reproduces its
own $\log_{10} V_0$ column on all 15 rows** — 14 of them exactly, the fifteenth
being the same number with its minus sign lost — so 14 of the 36
sedimentation rows carry an identically zero residual and cannot detect a wall
effect of any size — the sedimentation fit is therefore quoted on Table I alone,
at the several-times-wider interval §3 prints, and the verdict is unchanged.
The fluidisation number is a refit of the data Fig. 21 was drawn from; the
sedimentation number is not. And on the $d/D$ range the two data sets share, the
fluidisation interval **excludes** eq. (41)'s coefficient of 1: $\alpha$ is
roughly twice that coefficient below $d/D = 0.05$ and roughly equal to it above,
so the *linear* form of eq. (41) does not hold across its own range even though
the two-law verdict does. All four numbers are printed in §3.

**Three arithmetic findings, all small and all worth knowing.** The 17.5 and 18
of eqs. (37) and (38) are 17.4 and 17.8 rounded up, from the division displayed
on the line above them. The mean of Table VII is 0.2457, not the 0.25 the text
quotes. And the branch boundary at $Re = 200$ is discontinuous by up to 14 % in
$n$ because eq. (39) drops the wall term — which anyone sweeping $Re$ through 200
in a solver will meet, and which is invisible at $d/D = 0$ where it is usually
plotted.

**The flux function integrated, with the honest accounting of what that shows.**
$F(\phi) = \phi V_0(1-\phi)^n$ is what the correlation becomes downstream, in
every thickener and batch-settling model. The pymrm solve converges at order 1 on
the interface and its global solids balance is exactly zero — and the break table
shows that the balance is exactly zero for a wrong $n$, a wrong $V_0$, an
unstable step and a 16-fold-too-coarse grid as well. The solve tests the
integration, not the physics, and the page states that rather than banking the
$10^{-16}$.""" ))

# -------------------------------------------------------------------- caveats
cells.append(md(r"""## What this page does not claim

**It is not experimentally validated.** The data are the authors' own
measurements, so the *data tier* is 2 — but eqs. (33)-(39) were fitted to these
very slopes, and the 1.2 % is a fit residual. There is no hold-out set. The only
rows the authors excluded, the two oils at high $d/D$, are excluded for a reason
that is perfectly confounded with the correlation's own stated $d/D$ limit, so
they cannot serve as one either.

**It says nothing about $\varepsilon$ outside roughly 0.42 to 0.96.** Equation
(28) is a straight line on $\log V_c$ against $\log\varepsilon$ over whatever
range each run covered, and the tables record only its slope and intercept. The
most extreme printed abscissa tick this page read on any of those figures is
$\log\varepsilon = -0.38$ (Fig. 12, $\varepsilon = 0.42$), the least extreme
$-0.02$ ($\varepsilon = 0.96$) — those are **printed axis labels, not digitised
coordinates**, and no point was read off any figure. Everything the pymrm section
does below that voidage, which is the whole basal sediment, is outside the
correlation's measured range and is labelled so.

**It says nothing about non-uniform or non-spherical suspensions beyond
Table VIII.** The paper's own summary is explicit: "The present study has been
confined to relatively large particles ($d > 100$ microns)"; "The sedimentation
tests have all been carried out with uniform particles". Its mixed-size result is
a single figure (Fig. 15) with no tabulated numbers, and this page does not
reproduce it.

**The wall verdict is about $\alpha$, not about the mechanism.** The fits show a
$d/D$ term in the fluidisation intercepts and none in the sedimentation
intercepts. That the difference is "a velocity gradient created in the liquid
because of the drag exerted by the walls" is the authors' explanation, and
nothing in the tables tests it.

**And it is about the *presence* of the term, not the coefficient 1.** The
sedimentation half of that verdict rests on Table I alone, because Table II's
intercept column reproduces its own $\log_{10} V_0$ column and those 14 rows
cannot show a wall effect of any size; the interval on Table I is wide, and it is
the one quoted. The fluidisation half is a refit of the figure eq. (41) was read
off, and on the $d/D$ range the two data sets share it excludes 1 rather than
containing it. The page claims a wall term in fluidisation and none in
sedimentation. It does not claim that $\log V_0 - \log V_i$ equals $d/D$.

**The eq. (37) branch is not adjudicated.** Its five rows cannot separate the
correlation from a constant, in either direction, and §2 prints the three
sensitivities that show why. Nothing on this page says eq. (37) is worse than a
constant; an earlier draft did, and it was withdrawn.

**Table V's $d_s$ column is illegible in this scan** and no result on this page
depends on it. That is why §5's $K$ finding is stated as "eq. (42)'s abscissa
cannot be reconstructed for three of its seven points" rather than as "the
printed $K$ is wrong".

**Three quantities used here are not printed in the paper.** Water at 20 °C
($\mu = 1$ cP, $\rho = 1$), used only to re-check Table III's $Re$ column; the
two column diameters in cm, fixed to 0.3 % by the printed $d/D$ column; and
$g = 981$ cm s$^{-2}$, used only for Stokes' law in §1. None enters a result.

**The batch-settling balance is not Richardson and Zaki's.** They write no
conservation law. The balance in the pymrm section is derived on the page from a
solids volume balance; the kinematic-wave treatment of batch sedimentation is
standard and older than this paper, and no specific source is cited for it
because none was consulted. What is theirs is $F(\phi)$.

**No figure was digitised, and no claim here rests on one.** The only
figure-derived quantities anywhere on the page are printed axis tick labels,
quoted as printed text.""" ))

# --------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

```python
# The correlation, as printed.  Re is built on the TERMINAL velocity, so n is a
# property of the particle/liquid/tube combination and does not move with u.
def n_richardson_zaki(Re, dD):
    if Re <= 0.2:  return 4.65 + 19.5 * dD                    # eq. (33)
    if Re <= 1.0:  return (4.35 + 17.5 * dD) * Re**-0.03      # eq. (37)
    if Re <= 200:  return (4.45 + 18.0 * dD) * Re**-0.1       # eq. (38)
    if Re <= 500:  return 4.45 * Re**-0.1                     # eq. (39)
    return 2.39                                               # eq. (34)

# The authors' own recipe, page S97.  Do NOT drop the wall term in fluidisation.
V_i = V_0                              if sedimentation else V_0 * 10**(-d / D)   # eqs. (40), (41)
n   = n_richardson_zaki(V_0 * d * rho / mu, d / D)
V_c = V_i * eps**n                                                                # eq. (28)

# As a flux, which is what downstream models want:
F = lambda phi: phi * V_0 * (1.0 - phi)**n
```

**Three things to check before carrying this anywhere.**

1. **Which velocity you have.** The famous $u/u_t = \varepsilon^n$ is eq. (28)
   with $V_i$ replaced by $V_0$, i.e. eq. (40). In a column with $d/D = 0.04$ the
   wall term of eq. (41) is a factor $10^{-0.04} = 0.912$ — a 9 % error, larger
   than the correlation's own residual on this page.
2. **Where your $Re$ sits.** At $Re = 200$ the index jumps by up to 14 % because
   eq. (39) has no $d/D$ in it. At $Re = 0.2$ and $Re = 1$ there are smaller
   jumps. If you are sweeping $Re$, smooth the branches or expect kinks.
3. **How far outside $d/D \le 0.04$ you are.** That is the authors' stated limit
   and the only rows in the tables that exceed it inside the eq. (38) window are
   the ones the authors themselves discard.

**Where this connects in the gallery.**

- [`A1.1`](../A1.1-ergun-pressure-drop/) — Ergun. The fixed-bed end of the same
  $\varepsilon(u)$ curve; the bed stops obeying Ergun and starts obeying
  Richardson-Zaki at minimum fluidisation.
- `A1.6` — Wen and Yu's minimum fluidisation velocity. The other end of the
  matching condition, and **not built here**.
- [`A1.7`](../A1.7-geldart-classification/) — Geldart. Decides which fluidised-bed
  model applies at all; its group A is defined by *smooth expansion before
  bubbling*, which is the regime this correlation describes.
- `E1.1` — the two-phase theory. The gas-solid analogue of the expansion
  question; parked, see its case file.

**What to change first if you reuse this.** $V_0$. Every branch boundary and
every wall correction is computed from it, and the paper obtains it either by
measurement or from published drag data it does not reproduce. A 5 % error in
$V_0$ moves $Re$ by 5 %, which is enough to cross a branch boundary near
$Re = 1$, 200 or 500.""" ))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.update({
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
})
nbf.write(nb, "index.ipynb")
print("wrote index.ipynb")
