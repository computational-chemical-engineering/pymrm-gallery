# Findings from the priority-1 PDFs

Recorded 2026-07-26 after inspecting the three priority-1 papers. **No numbers
were transcribed into any dataset from these** — see the OCR warning below.

> **Location: `~/papers/pymrm-gallery/`** (moved there 2026-07-26; not in the
> repo, and must stay out — `.gitignore` blocks `*.pdf` and
> `scripts/check_metadata.py` errors on committed ones).

## Inventory and text-layer quality

Priorities 1 and 2 are all supplied. The characters-per-page figure is a good
predictor of whether tables will transcribe mechanically or need reading from a
page image: below roughly 5,000 expect trouble with sub/superscripts.

| chars/page | File | Paper | Page |
|---|---|---|---|
| 12,403 | `i260028a001.pdf` | Kunii & Levenspiel (1968) | `E2.1` |
| 8,987 | `i300005a006.pdf` | Oh & Cavendish (1982) | `I1.2` |
| 8,578 | `AIChE…1962…Duncan…pdf` | Duncan & Toor (1962) | `A4.9` ✔ done |
| 7,422 | `AIChE Journal - January 1989 - Xu.pdf` | Xu & Froment (1989) | `C2.1` |
| 6,773 | `AIChE…1987…Itoh…pdf` | Itoh (1987) | `H1.4` |
| 6,408 | `AIChE_Journal-1996-Krishna.pdf` | Krishna & Ellenberger (1996) | `F1.4` |
| 5,635 | `1-s2.0-0009250978851203-main.pdf` | Wakao & Funazkri (1978) | `A3.4` |
| 4,595 | `1-s2.0-0009250970850734-main.pdf` | Van Welsenaere & Froment (1970) | `D2.2` |
| 4,432 | `1-s2.0-0009250962850052-main.pdf` | Weisz & Hicks (1962) | `B1.1` ✔ published |

The two ACS scans have the best text layers by a wide margin, so **Kunii &
Levenspiel and Oh & Cavendish are the cheapest of the remaining pages to start**
— worth weighing against Xu & Froment's higher scientific value but harder
extraction. Note that a high character count does not guarantee correct
exponents; the Xu & Froment file scores 7,422 and still mangles them.

---

## 1. Xu & Froment (1989) — page `C2.1` ✔ built 2026-07-26

`AIChE Journal - January 1989 - Xu.pdf`

**Structure of the paper.** Table 1 experimental conditions · Table 2 the eleven
candidate reactions · Table 3 the retained mechanism · **Table 4 parameter
estimates per temperature** · Table 5 Arrhenius/van 't Hoff fits · **Table 6 the
final parameters for the fresh catalyst**.

**Table 6 is the one the page needs** — pre-exponential factors and activation
energies for the three reactions plus the adsorption constants.

**Blocker: the scan OCRs badly.** `pdftotext` renders exponents as literal
text — `8.664 lo-'` where the value is `8.664 × 10^-n` — so the parameter tables
cannot be transcribed mechanically without corrupting every number. Table 4 came
out partially legible; Table 6 did not resolve at all in the text layer.

**How to proceed:** render the relevant page at 600 dpi and read it as an image,
the same route used for Duncan & Toor Figure 2. Do **not** attempt to repair the
OCR by inference — a mis-read exponent is a silently wrong rate constant.

**Three subtleties already found in the text, worth carrying into the page:**

1. The steam-reforming runs used a **partially deactivated** catalyst. Table 6
   gives fresh-catalyst parameters obtained by extrapolating conversion to
   *t* = 0, which the authors themselves call "not too accurate".
2. The reverse water-gas-shift and methanation curves in Figures 4 and 5 were
   drawn after **multiplying the Table 6 rate parameters by 1.225**. Any
   reproduction of those figures must apply the same factor or it will not match.
3. Reference temperatures differ by parameter group: *T*ᵣ = 648 K for the rate
   coefficients and for *K*꜀ₒ and *K*ₕ₂, but 823 K for *K*_CH₄ and *K*_H₂O.
   Using one reference for all of them is a natural mistake and would bias the
   temperature dependence.

`K_CO2` was never statistically significant and correctly has **no** term in the
denominator of the rate equations.

### What the page actually needed, and how it went

All three subtleties above survived contact with the model and are on the page.
Two further things were found while building it:

4. **The paper never tabulates the equilibrium constants** *K*₁, *K*₂, *K*₃,
   although they appear in all three rate equations. They have to come from
   outside the paper. The check that makes this defensible is that
   d ln *K*/d(1/*T*) = −Δ*H*/*R*, so the correlation's slopes *are* reaction
   enthalpies and can be compared with the paper's own Table 7: 0.42 %, 1.93 %
   and 0.54 % for reactions I, II and III. This turned out to matter — see
   below.
5. **Table 5 ⇄ Table 6 is a free check on the page-image reading.** The two
   tables are related by *A* = value(*T*ᵣ)·exp(*E*/(*RT*ᵣ)). Recomputing all
   seven Table 6 pre-exponentials from Table 5 reproduces them to ≤1.44 %.
   Because *T*ᵣ is *in* the formula, the same check catches trap 3: forcing
   648 K on every parameter puts *K*_CH₄ out by a factor of 0.22 and *K*_H₂O by
   a factor of 33.

**Extraction, in the end.** Journal page 94 (PDF page 7) at 600 dpi was legible
without any image processing — both tables read directly, exponents and all.
Figures 2 and 3 (page 89) were much harder than Duncan & Toor: the markers are
~20 px glyphs on ~10 px curves, so morphological opening does not separate them.
What worked was a *local ink excess* — box-filtered ink density minus the
largest value explainable by a locally straight structure (the maximum over six
orientations of a grey-scale opening with a long line element) — followed by a
visual audit of every candidate at 600 dpi. The orientation maximum is the
essential part: a single horizontal element flags the whole steep near-origin
section as marker.

**Two independent checks made the digitisation trustworthy.** Figures 2 and 3
plot the same runs, so every Figure-2 point must have a Figure-3 partner at the
same space time — all 30 pair to within 0.0027 in *W*/*F*, from independently
fitted axis calibrations. And *x*_CO₂ ≤ *x*_CH₄ must hold: it does, everywhere.

**Result:** 0.0017 mean absolute deviation in conversion over 61 points, nothing
fitted. But the residual is not flat — it grows with space time and tracks the
approach to equilibrium, and a 5 % shift in the (external) equilibrium constants
moves the long-space-time predictions by the same amount. The bias there belongs
to finding 4, not to the kinetics.

## 2. Krishna & Ellenberger (1996) — page `F1.4` — reviewed, approach revised

`AIChE_Journal-1996-Krishna.pdf`

**Table 1** experimental setup, operating conditions and system properties —
this is the column-geometry table the page needs. **Table 2** physical
properties of the liquids used. **Table 3** average relative deviation of their
model against the Wilkinson correlation.

**Scale: 2,787 experiments, of which 1,735 in the churn-turbulent regime.**

**Table 1 has a free integrity check, and it passes.** The nineteen per-row
experiment counts — 346, 22, 60, 128, 460, 239, 209, 293, 37, 153, 63, 185, 99,
119, 70, 77, 64, 92, 71 — sum to exactly **2,787**, the total stated in the
text. Transcribe the column and add it up before doing anything else.

Watch three OCR slips in Table 1 that the page image corrects: `air-paraffin
oil (3)` is **(B)**, `SF, -tetradecane` is **SF₆**, and `0.00 1-0.249` is
**0.001–0.249**.

### The complete model, read from the page images

Wilkinson et al. (1992), their Eqs. 1–4:

- (1) ε_b = (*U* − *U*_trans)/*V*_b ; ε = ε_small + ε_b
- (2) *U*_trans = ε_trans *V*_small ; ε_trans = 0.5 exp(−193 ρ_G^−0.61 μ_L^0.5 σ^0.11)
- (3) *V*_small μ_L/σ = 2.25 (σ³ρ_L/(g μ_L⁴))^−0.273 (ρ_L/ρ_G)^0.03
- (4) *V*_b μ_L/σ = *V*_small μ_L/σ + 2.4 ((*U*−*U*_trans)μ_L/σ)^0.757 (σ³ρ_L/(g μ_L⁴))^−0.077 (ρ_L/ρ_G)^0.077

Reilly et al. (1994) for the transition, their Eq. 8, with **B = 3.85** (stated
in their text, not in the notation list):

- ε_trans = 0.59 *B*^1.5 √(ρ_G^0.96/ρ_L · σ^0.12)
- *V*_small = (1/2.84) ρ_G^−0.04 σ^0.12
- *U*_trans = *V*_small ε_trans (1 − ε_trans)

Their own large-bubble correlation, **Eq. 19** — the point of the paper:

> ε_b = 0.268 · *D*_T^−0.18 · (*U* − *U*_df)^−0.22 · (*U* − *U*_df)^4/5

i.e. ε_b = 0.268 *D*_T^−0.18 (*U* − *U*_df)^0.58, since 4/5 − 0.22 = 0.58. Valid
for *U* > 0.1 m/s and dispersion heights *H* > 1 m. Note it contains **no fluid
property at all** — that is the claim.

Derivation chain, if the page wants to show it: Eq. 12 (ε_b as an integral of
(*U*−*U*_df)/*V*_b over the dispersion height) → Eq. 13 (*V*_b = φ√(g d_b)) →
Eq. 14 (φ = φ₀ *D*_T^N) → Eq. 15 → Eq. 16 for *H* ≫ *h*\* → Eq. 17 with
*h*\* = b0(*U*−*U*_df)^b1 → Eq. 19 is Eq. 17 fitted.

### What is still missing, and it is the expensive part

**Table 3 cannot be reproduced without the raw data.** Its δ values (1735
points: Wilkinson 1.25 vs their model 0.16 for ε_b; 0.24 vs 0.23 for ε) are
averages over the 1,735 churn-turbulent measurements, and those measurements
exist only as scatter in Figures 7, 9 and 11. **Without digitising a figure this
page has no data and must not be published** — evaluating both correlations and
observing that one lies above the other demonstrates nothing.

Best digitisation targets, in order: **Figure 11** (large-bubble holdup vs
*U*−*U*_df at four gas densities, log–log, one column, the cleanest test of the
gas-density claim), then **Figure 7a/7b** (three column diameters, tests the
*D*_T^−0.18 exponent directly), then Figure 9. Expect the `C2.1` marker
extraction rather than the `A4.9` one, and expect worse: these are dense
overlapping scatter plots on log axes with several series per panel.

## 3. Weisz & Hicks (1962) — page `B1.1` (already published)

`1-s2.0-0009250962850052-main.pdf`

**No tables found in the text layer** — the η(φ) results are presented as
figures only, as expected for a 1962 computational paper.

This is **not blocking**: `B1.1` is already validated against exact isothermal
solutions and an independent shooting reference, which is a stronger test than
agreement with someone else's 1962 numerics. Digitising their curves would add a
third comparison and, more usefully, confirm that our β and γ conventions match
theirs — worth doing but low priority.

## 4. Van Welsenaere & Froment (1970) — page `D2.2` ✔ built 2026-07-27

`1-s2.0-0009250970850734-main.pdf`, and full text via the Elsevier PII endpoint.

**Structure.** Seven sections and an appendix, **no tables at all**, ten
figures. §1 model and base parameters · §2 the *p*–*T* phase plane · §3 critical
inlet conditions (the two criteria) · §4 subcritical conditions · §5 comparison
with Barkelew · §6 four worked numerical examples · §7 conclusion · Appendix I
on the locus of maxima.

**The API text is unusable for numbers.** It is the publisher's OCR of the same
scan and it discards the 1970 mid-dot decimal separator: `R = 00125 m` for
0.0125 m, `b = 19837` for 19.837, `(p°),,, = 001353 atm` for 0.01353 atm. It is
excellent for prose and section structure. Read every number off a 600 dpi page
render — the same discipline as `C2.1`.

**Base parameter set** (journal page 1504, read from the page image at 600 dpi,
not transcribed from the API text):

| | |
|---|---|
| *M* | 29.48 kg/kmol |
| *P* | 1 atm |
| ρ_b | 1300 kg/m³ |
| ρ_g | 1.293 kg/m³ |
| (−Δ*H*) | 307 000 kcal/kmol |
| *c*_p | 0.323 kcal/(m³·°C) — **volumetric**, as printed |
| *u* | 3600 m/hr |
| *U* | 82.7 kcal/(m²·hr·°C) |
| *R* | 0.0125 m |
| *p*_B⁰ | 0.208 atm |
| *a* | 13 636 K⁻¹ (= *E*/*R*_gas) |
| *b* | 19.837 |

Note *c*_p is printed as kcal/m³·°C, i.e. already ρ*c*_p, which is what makes
*C* = 2*U*/(*c*_p·*R*) come out in 1/hr against *z* = *z*′/*u* in hr. Do not
multiply by ρ_g again.

**Model** (their Eqs. 3–4), pseudo-first-order irreversible, one dimension,
constant wall temperature:

d*p*/d*z* = −*A p* e^(−a/T+b), d*T*/d*z* = *B p* e^(−a/T+b) − *C*(*T*−*T*_w),
with *A* = *M P* ρ_b/ρ_g · *p*_B⁰, *B* = (−Δ*H*)ρ_b/*c*_p · *p*_B⁰,
*C* = 2*U*/(*c*_p *R*), *z* = *z*′/*u*, and *p* = *p*⁰, *T* = *T*₀ = *T*_w at
*z* = 0.

**The first criterion, complete** (all read from the page images, journal pages
1506 and 1509–1511):

| Eq. | Formula |
|---|---|
| 6 | maxima curve: *B p*_m e^(−a/T_m+b) − *C*(*T*_m−*T*_w) = 0 |
| 7 | *p*_m = (*T*_m−*T*_w) / [(*B*/*C*) e^(−a/T_m+b)] |
| **8** | ***T*_M = ½[*a* − √(*a*(*a*−4*T*_w))]** — the critical hot spot |
| 28 | (*p*⁰)_u = (*A*/*B*)(*T*_M−*T*_w)·[1/√*X* + 1]², *X* = (*A*/*C*)e^(−a/T_M+b) |
| 29 | (*p*⁰)_l = (*T*_M−*T*_w)/[(*B*/*C*)e^(−a/T_M+b)] + (*A*/*B*)(*T*_M−*T*_w) |
| 30 | (*p*⁰)_m = (*A*/*B*)(*T*_M−*T*_w)·[1 + 1/√*X* + 1/*X*] — equals the mean of 28 and 29 |
| 31–34 | (Δ*T*)_ad/(Δ*T*)_eff = 1+*Q*+*Q*², (Δ*T*)_ad = (*B*/*A*)*p*⁰, (Δ*T*)_eff = *T*_M−*T*_w, *Q* = 1/√*X* |

**These were checked numerically and every published number comes back.** With
the parameter set above and *T*_w = 625 K:

| Quantity | Recomputed | Paper |
|---|---|---|
| ln *K* (Eq. 22) | −2.0568 | −2.055 |
| *t*_w (Eq. 21b) | 21.8176 | 21.818 |
| *T*_M (Eq. 8) | 656.62 K | 656.6 K |
| (*p*⁰)_l (Eq. 29) | 0.01353 atm | 0.01353 |
| (*p*⁰)_u (Eq. 28) | 0.01976 atm | 0.01976 |
| (*p*⁰)_m (Eq. 30) | 0.01664 atm | 0.01665 |
| Example 4: (Δ*T*)_ad, *Q*, *R* | 521.10 K, 3.4662, 0.01752 m | 521.09, 3.4675, 0.0175 m |
| Example 3: (Δ*T*)_ad | 312.66 K | 312.6 |

That is the parameter reading *and* all seven formulae confirmed at once, before
a line of the page is written. The script is in the session scratchpad but is
two dozen lines and trivial to recreate from the table above.

**The second criterion** needs their Eq. 20, *K* = (*t*−2)/(*t* e^(−t+20)) ·
[1 − *t*(1 − *t*/*t*_w)], with *t* = *a*/(*T*_i)_t and *K* = (*A*/*C*)e^(b−20).
The paper solves it graphically through Figs. 6 and 7 explicitly to avoid
needing a computer. Solve it directly with `brentq` instead and **Figs. 6 and 7
become validation targets rather than inputs** — that is the natural "what pymrm
adds" for this page, alongside a runaway-boundary sweep over (*T*_w, *p*⁰) that
the original could only sketch as Fig. 8.

**Validation, in order of value.** §6's four worked examples (above) test the
criteria and need no digitising at all. Fig. 8 — upper limit, lower limit, mean
and exact critical *p*⁰ against *T*_w over 600–700 K for both criteria, with
Barkelew's values marked ⊕ — is the richer target and would also let the page
reproduce the paper's comparison with Barkelew. **Read §6 from the page image**;
it is where the lost decimal points do the most damage.

**Provenance note.** These are the authors' own computed values, not
measurements, so `D2.2` is a tier-6 page like `B1.1` and `F3.1` unless Fig. 9
(Barkelew's diagram) is treated as an independent reference. Do not describe it
as experimentally validated.

---

## Recommended order for the next session

1. ~~**Xu & Froment Table 6** via page-image reading → build `C2.1`.~~ **Done
   2026-07-26.**
2. ~~**Van Welsenaere & Froment** `D2.2`.~~ **Done 2026-07-27** — and it cost
   the correction recorded above: the Elsevier API is *not* a way around OCR for
   numbers.
3. **Krishna & Ellenberger Tables 1–2** transcribe (clean text), then digitise
   the holdup figures → build `F1.4`. Expect the marker extraction to be closer
   to `C2.1` than to `A4.9`; `pages/C2.1-xu-froment-smr/extract_figures.py` is
   the starting point, not morphological opening. This is the next chance to
   add an experimentally validated page — 2,787 runs.
4. **Kunii & Levenspiel** `E2.1` — best text layer of the set.
5. Weisz & Hicks figures only if time allows.

---

## Batch triage, 2026-07-27 — ten papers arrived at once

Run `python scripts/probe_paper.py ~/papers/pymrm-gallery/*.pdf` to regenerate
this. The question it answers is *where do this paper's numbers live*, because
that, not page count, sets the cost of the page.

| Catalog | Paper | Text layer | Figures | Cost |
|---|---|---|---|---|
| `G1.8` | Herskowitz & Smith 1983 | 9 500 ch/pg, Tables 1–6 | 9 | tables transcribe directly |
| `J3.4` | Doyle–Fuller–Newman 1993 | 8 460 ch/pg | 8 | parameters in running text |
| `J4.8` | Henze et al. 1987 (ASM1) | 6 836 ch/pg, Tables 1–4 | 2 | the model *is* a table |
| `A4.2` | Krishna & Wesselingh 1997 | 6 883 ch/pg | 71 | a 51-page review, not one model |
| `B3.1` | Yagi & Kunii 1955 | 6 035 ch/pg, Tables 1–4 | 11 | tables transcribe directly |
| `H1.7` | Wijmans & Baker 1995 | 4 784 ch/pg | 12 | thin — read numbers off renders |
| `F2.3` | Maretto & Krishna 1999 | 4 420 ch/pg, Table 1 | 8 | thin — read numbers off renders |
| `J1.5` | Glueckauf 1955 | 2 469 ch/pg, Table 2 | 2 | thin — read numbers off renders |
| `A2.3` | Taylor 1953 | 2 469 ch/pg equivalent | few | thin, and hyphens for decimals |
| `A1.1` | Ergun 1952 | **none at all** | — | pure scan, everything off renders |

**The decimal-point trap is not confined to Elsevier.** Taylor's 1953 Royal
Society scan renders `48.0` as `48-0` and `59.8` as `59-8`. Same failure as Van
Welsenaere & Froment, different publisher, different decade. Assume it for any
pre-1980 scan and read the numbers off a render.

### Which of these can be built without a human in the loop

The gate is figure digitisation, so the split is:

- **`A2.3` Taylor–Aris** — validation is Taylor's own closed form
  *K* = *a*²*U*²/(48*D*) plus his measured table. Table, not figure. **No gate.**
- **`J1.5` Glueckauf** — the LDF coefficient 15*D*/*r*² against the exact
  spherical-diffusion series. Purely analytic. **No gate.**
- **`J4.8` ASM1** — the model is a stoichiometric matrix; validation is the
  continuity relations the matrix must satisfy. **No gate**, but a big model.
- **`A1.1` Ergun** — a one-line correlation, but its data is in a scanned figure
  and table. Needs the gate.
- Everything else — `F2.3`, `H1.7`, `G1.8`, `B3.1`, `J3.4` — validates against
  figures and needs the gate.

`A4.2` is a 51-page review with 71 figures and is not a single model; it is a
reference for the `S9` pages rather than a page of its own. Do not queue it as
one.

---

## `F1.4` Figure 11 — what the review changed

The digitisation was put to the maintainer as an overlay on the original figure
(private artifact, not committed). The verdict was **"positions fine, shapes
wrong — do not trust the gas labels"**, with these specifics:

- **SF₆ stops at about *U* − *U*_df = 0.05 m/s.** It is well separated and
  entirely on the left. Everything the detector called SF₆ to the right of that
  is squares and circles — 9 of its 21 SF₆ calls.
- **Most triangles (helium) lie below the eq. (19) line.** Above the line in the
  dense band the detector's triangles are wrong; those are circles and squares.
- Circles and squares are found where there is little overlap and missed in the
  dense region — hence only 6 argon and 4 air out of 63.
- One triangle is missed entirely: the second from the left, drawn inside a
  square.

**The fix is not a better classifier — it is to stop needing the labels.**
Eq. 19 contains no gas-density term at all:

> ε_b = 0.268 *D*_T^−0.18 (*U* − *U*_df)^0.58

so testing it against the extracted *positions*, ignoring which marker shape
each came from, is not a workaround — it is the correct test, and it uses all
63 points instead of a curated subset.

Doing that gives **13.8 % mean absolute deviation with a bias of +2.8 %**, which
sits alongside the δ = 0.16 the authors themselves report in their Table 3 for
this correlation over all 1,735 churn-turbulent runs.

**Two numbers in an earlier version of this section were wrong.** They are left
here corrected rather than deleted, because both mistakes generalise.

*Deviation direction.* The first pass computed measured/model for Eq. 19 and
model/measured for Wilkinson. At 14 % scatter a ratio and its reciprocal are not
interchangeable: the mean moved 13.3 → 13.8 % and the bias −0.2 → +2.8 %. Fix
one convention — the page uses (model − measured)/measured — and apply it to
every correlation being compared.

*The gas-density test was confounded.* The first pass split at 0.05 m/s and
reported the SF₆ group at −3.0 % against everything else at +0.5 %, a 3.5 %
difference across a 37-fold density range. That number measures **velocity, not
density**: on this figure the groups do not overlap at all in the abscissa —
every SF₆ point lies below 0.044 m/s and every other point above 0.051. If
Eq. 19's exponent of 0.58 is slightly off at the low end, the SF₆ group looks
displaced for a reason unrelated to gas.

The test that does work is extrapolation. Fit a free power law to the 51
helium/air/argon points only, then predict the 12 SF₆ points — 3.7 to 37 times
denser, and unseen by the fit. The result is a bias of **−7.7 % with 13.5 %
scatter**, the same scatter as the correlation itself, so SF₆ is not an outlier
group. Wilkinson would require 20–27 %. Part of the 7.7 % is extrapolation error
(the SF₆ window sits about half a decade below the fitted range), so it is an
upper bound on any real gas effect, and the page says so.

**Carry three things forward.**

1. Before building a shape classifier for a scatter figure, check whether the
   model being tested actually needs the series identity. Often the correlation
   under test has no term for whatever distinguishes the series, and then the
   labels are decoration.
2. Before reporting a difference between two groups, check the groups overlap in
   every *other* variable. A clean-looking group contrast can be measuring the
   confound.
3. Fix one deviation convention per page and state it. Reciprocals diverge as
   soon as the scatter is more than a few percent.
