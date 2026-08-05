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

### 600 dpi is not always 600 dpi — check the embedded image first

Itoh 1987 (`H1.4`) is a **300 dpi bilevel CCITT** scan. `pdftoppm -r 600` on it is a
2× upsample and carries **no information the 300 dpi bitmap does not already
have** — the extra pixels are interpolation. A verifier chasing a marginal
superscript had to read the native bitmap as an ASCII map instead, and that read
overturned a claimed misprint: the exponent the page said was mis-set is printed
correctly.

So "read it at 600 dpi" is shorthand for *read the image, not the text layer*, not
a guarantee of resolution. Before quoting a marginal glyph, check what is actually
in the file:

```bash
pdfimages -list f.pdf | head        # width, height, bpc, encoding per image
```

If `bpc` is 1 (bilevel) the anti-aliasing that makes a glyph legible at 600 dpi
does not exist, and a rendered upsample can invent apparent shapes. Render at the
image's own resolution and, for a genuinely ambiguous character, compare it with
an unambiguous instance of both candidates elsewhere on the same line — that is
what settled `-5` against `-7` here.

**Two more files in the same class, found 2026-08-02/03:** the MFIX Theory Guide
(`A1.8`) is bilevel CCITT at **400 dpi** on all 54 pages, and the Richardson &
Zaki Golden Jubilee reprint (`A1.5`, `A1.8`) carries a bilevel **300 dpi**
stencil over a 200 dpi grey background. Both pages read their constants at those
resolutions, correctly. Better than rendering at all: `pdfimages -png` extracts
the stored bitmap and bypasses the rasteriser.

### Look for the 1995 Golden Jubilee reprint before rendering pages

*Chemical Engineering Science* reprinted a set of its classics **verbatim** in the
1995 Golden Jubilee issue (vol. 50), re-typeset rather than re-scanned. Where the
original scan OCRs badly, the reprint OCRs cleanly, and the two printings can be
diffed against each other — which turns a guess into a decision. Aris 1957 renders
`0.698` as `898` in the 1957 scan and correctly in the 1995 reprint
(PII `0009250996818197`); an agent transcribed only numbers the two printings
agree on, and used a disagreement at *p* = 0.75 to work out which printing was
right by checking both against the paper's own stated limits.

Search the PII endpoint for the vol. 50 (1995) reprint of any pre-1970 CES paper
before falling back to 600 dpi renders. It is faster and it gives a second
independent witness for every constant.

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

---

## The Elsevier API text is not "safe after 1980" — measured on 15 articles, 2026-08-02

The correction recorded in [`handoff.md`](handoff.md#elsevier-full-text--use-this-instead-of-ocr)
says the API text drops decimal points "for pre-1980 scans". That framing is too
narrow and has already misled one dispatch. **The cut is not a date, it is
whether the article is a scan or born-digital**, and Elsevier was still scanning
into the 1990s.

Fifteen articles were fetched with the authorised key and inspected. The
diagnostic that separates them in one line is the count of well-formed decimals
(`\d+\.\d+`) against the character count:

| Article | chars | well-formed decimals | verdict |
|---|---|---|---|
| Krishna & Baur (2003), Sep. Purif. Technol. 33 | 218 k | many | **born-digital, usable for numbers** |
| Hulburt & Katz (1964), CES 19 | 84 k | 335 | theory paper, no tables; prose good |
| Uppal, Ray & Poore (1974), CES 29 | 67 k | 80 | body mangled, **figure labels clean** |
| Uppal, Ray & Poore (1976), CES 31 | 43 k | 79 | as above |
| Robeson (1991), J. Membr. Sci. 62 | 70 k | 26 | **decimals gone**, and this is 1991 |
| Szekely & Evans (1970), CES 25 | 52 k | 19 | prose only |
| Mars & van Krevelen (1954), CES 3 | 72 k | 17 | worst of the set, see below |
| Froment & Bischoff (1961), CES 16 | 46 k | 13 | prose only |
| van den Broeke & Krishna (1995), CES 50 | 73 k | — | scan; **Table 4 absent entirely** |
| van Deemter et al. (1956), CES 5 | 62 k | 7 | tables present, numerically destroyed |
| Whitman (1923 / IJHMT 1962) | 19 k | 7 | prose good, data are in figures |
| Nývlt (1968), J. Cryst. Growth 3 | 25 k | 7 | prose only |
| Calderbank & Moo-Young (1961), CES 16 | 81 k | **5** | 81 kB of text, five readable numbers |

**Robeson 1991 is the one to remember.** Sixty-nine kilobytes of clean-looking
English, and Table 2's Lennard-Jones kinetic diameters arrive as

> `Gas He H2 CO2 OZ NZ CH4 Kinetic diameter (A) 2 6 2 89 3 3 3 46 3 64 3 8`

for 2.6, 2.89, 3.3, 3.46, 3.64, 3.8 — the decimal point is rendered as a space,
so a six-value row becomes an eleven-value row and nothing looks wrong until you
count. Symbols degrade the same way: `P,` for *P*ₓ and `at,` for α.

**Mars & van Krevelen 1954 is the extreme.** The 1954 typesetting uses a mid-dot
decimal separator and the OCR picks a *different glyph for it within a single
table row*:

> `825 150 82·5 0·095 0·107 300 73'0 0-120 500 59·0 0'107`

`·`, `'` and `-` are all the same decimal point. Any rule of the form "replace ·
with ." silently corrupts two-thirds of the table.

**van Deemter 1956** loses them entirely: Table 1 reads `0 163 1100 0 2? 0 30`,
with `?` standing in for digits the OCR could not resolve. Both of its tables are
the experimental data a page would validate against, so the API route buys this
case nothing at all.

### There is no page-image fallback for these

The obvious remedy — render the page at 600 dpi, as for a Wiley scan — is not
available, because the API does not serve the page. Checked on every article
above:

- `GET /content/object/pii/<PII>` returns `{"choices":null}`. No page images, no
  figure objects.
- `GET /content/article/pii/<PII>` with `Accept: application/pdf` returns a
  **one-page entitlement preview**, and `content/article/entitlement/pii/<PII>`
  answers `AUTHENTICATION_ERROR — Requestor configuration settings insufficient
  for access to this resource`. This is the same behaviour already recorded for
  `B1.2` in [`papers-on-disk.yaml`](papers-on-disk.yaml).

So the key is a **prose-and-metadata** instrument, not an extraction instrument.

### What it is genuinely good for

Three things, and they are worth the request volume:

1. **Verifying identity.** The header carries title, journal, volume, issue and
   page range, and the body carries the author list. That is a title page in all
   but name, and it is how the `H1.9` mis-mapping and the `F1.6` / `F3.4`
   book-review DOIs were caught on 2026-08-02.
2. **Telling you which table a page needs**, so a PDF request can be ranked
   instead of guessed.
3. **Born-digital articles**, roughly 2000 onward, where it really is clean.
   Krishna & Baur (2003) returns Ð, Θ, θ, Γ, `Pa−1` and five tables with headers
   intact, and is being used as the source for `H1.9` and `A4.7`.

### The Elsevier text also drops whole tables

Not just characters. In van den Broeke & Krishna (1995) the body refers to
"Table 4" five times — it is the table of single-component Maxwell-Stefan
diffusivities, the one thing that case most needs — and **the table itself is not
in the returned text**. Tables 1, 2 and 3 are. So "the API returned 73 kB, the
paper is covered" is not a safe inference: check that each table you are counting
on has a *caption* in the text, not merely a cross-reference to it.

## The MFIX Theory Guide (`A1.8`) — good scan, bad text layer

`Syamlal-Rogers-OBrien-1993-MFIX-theory-guide-DOE-METC-94-1004.pdf`, 54 pp, from
OSTI. The document is legible on screen but its text layer is one of the worst
here: words are run together (`Fluid-SolidsMomentumTransfer`,
`ConservationofMass`), digits are substituted (`0.O6Re` with a capital O,
`Ergun (f952)` for 1952, `_g-2'6s` for ε_g^−2.65), and equation numbers come out
as `(121` and `(1''`. Read every constant off a page image — **at 400 dpi, not
600**, for the reason below.

**400 dpi is this file's ceiling, and 600 dpi would be interpolation.** Every one
of the 54 pages is a *single* CCITT-G4 image with `bpc = 1` (bilevel) at
**400 × 400 ppi** (a handful report 401 × 400 from rounding), 3520–3646 px wide by
4472–4572 tall:

```bash
pdfimages -list Syamlal-Rogers-OBrien-1993-MFIX-theory-guide-DOE-METC-94-1004.pdf
```

There is no higher native resolution to reach. `pdftoppm -r 600` upsamples a
1-bit image and invents grey that is not in the file. Better still, skip the
rasteriser entirely and pull the page image out as it is stored:

```bash
pdfimages -png -f 13 -l 14 <file.pdf> out    # eqs. (11)-(16), journal pages 10-11
```

**The same check on the Richardson & Zaki Golden Jubilee reprint**
(`1-s2.0-S0263876297800068-main.pdf`, `A1.5`'s source, whose Table VI `A1.8`
ships): each page carries a **`stencil` at 300 × 300 dpi, `bpc = 1`** over a
200 dpi grey background image. So 300 dpi is that file's native resolution too,
and `A1.8`'s "300 dpi renders" is right for the same reason.

*This is the `H1.4` lesson (below, "600 dpi is not always 600 dpi") turning up a
second and third time in one week. The `AGENTS.md` rule "read constants off
600 dpi renders" needs its precondition stated wherever it is written:* **or at
the embedded image's native resolution, whichever is lower — run `pdfimages
-list` first; for a bilevel scan, 600 dpi is interpolation.**

Its section 2.2.1 does carry the Syamlal–O'Brien drag law complete — eq. (11) for
the terminal-velocity-to-drag conversion, eq. (12) for Garside & Al-Dibouni's
closed-form *V*ᵣₘ with *A* and *B* in (13)–(14), and Dalla Valle's single-sphere
*C*_D. It does **not** contain Wen–Yu (the strings "Wen and Yu" and "Wen & Yu" do
not occur) and prints **no Gidaspow drag closure and no blend rule**, so it
sources one of the three drag laws the case names.

**Do not shorten that to "cites Gidaspow only in passing" — four files said so
and it is false.** Gidaspow (1986) is cited exactly once *in §2.2.1*, which is
the load-bearing part; but across the whole report the name occurs **28 times**
in the extracted text (18 body, 10 reference list) over ten reference entries,
four of them Gidaspow-first-author, and the report **adopts** a Ding & Gidaspow
(1990) expression as its own eq. (88) and a Syamlal & Gidaspow (1985)
conductivity model. The absence of a Gidaspow *drag* law is the claim; the
absence of Gidaspow is not. See `queue_cases/A1.8.yaml`.

---

## The 2026-08-05 second drop — sixteen files, per-file text-layer notes

Added by the source-mapping pass; the catalogue-ID map and the identification
evidence are in [`papers-on-disk.yaml`](papers-on-disk.yaml) and in each
`queue_cases/<ID>.yaml`. Full ranking in
[`source-map-2026-08-05.md`](source-map-2026-08-05.md).

### 600 dpi is wrong for every scanned file in this batch

`pdfimages -list` was run on all sixteen. **Every scanned file in this drop is
CCITT-G4 or JBIG2 bilevel at 300 ppi native**, except Fuller et al. (JPEG RGB,
also 300 ppi) and Li & Kwauk (JBIG2, ~285 ppi). So `pdftoppm -r 600` upsamples a
1-bit image on all of them and invents grey that is not in the file — the same
finding already recorded above for `A1.8`'s MFIX report (400 ppi) and for
Richardson & Zaki (300 ppi). **Render at 300 dpi, or pull the stored image out
directly with `pdfimages -png`.**

That now makes it three separate weeks running. The `AGENTS.md` rule "read
constants off 600 dpi renders" has never once been literally correct for a file
in this repository; the correct rule is *at the embedded image's native
resolution, or lower*.

### Per file

| File | Case | Pages | Text layer | Native | Verdict on numbers |
|---|---|---|---|---|---|
| Chapman & Cowling, *Non-Uniform Gases* 3rd edn | `A4.6` | 448 | ABBYY 8, ~2.0 k/pg | CCITT-G4 **bilevel 300 ppi** | prose good; **numbers no** |
| Li & Kwauk, *Particle-Fluid Two-Phase Flow* | `A1.9` | 214 | **NONE — 1 byte/page** | JBIG2 bilevel ~285 ppi | renders only, everything |
| Billet & Schultes 1999, Trans IChemE 77A | `A3.9` | 7 | **NONE — RC4, copy disabled** | CCITT-G4 bilevel 300 ppi | renders only; very legible |
| Dixon & Cresswell 1979, AIChE J 25 | `A3.11` | 14 | 4.8 k/pg | CCITT-G4 bilevel 300 ppi | prose ok; **digits destroyed** |
| Martin & Nilles 1993, CIT 65 (German) | `A3.14` | 10 | 4.5 k/pg | CCITT-G4 bilevel 300 ppi | see German section below |
| Westerterp et al. 1995, AIChE J 41 | `A2.7` | 16 | 5.4 k/pg | CCITT-G4 bilevel 300 ppi | prose excellent; **equations absent** |
| Zwietering 1959, CES 11 | `A2.8` | 15 | 3.2 k/pg | CCITT-G4 bilevel 300 ppi | **τ read as the digit 7 — see below** |
| Baldyga & Bourne 1989 I, Chem. Eng. J. 42 | `A2.9` | 10 | 3.4 k/pg | CCITT-G4 bilevel 300 ppi | prose good; **equation bodies dropped** |
| Baldyga & Bourne 1989 II, Chem. Eng. J. 42 | `A2.9` | 9 | 3.0 k/pg | CCITT-G4 bilevel 300 ppi | same |
| Rocha, Bravo & Fair 1993, IECR 32 | `A3.10` | 11 | 4.1 k/pg | CCITT-G4 bilevel 300 ppi | prose fair; **figure captions worst** |
| Fuller, Schettler & Giddings 1966, IEC 58 | `A4.5` | 10 | 6.0 k/pg | **JPEG RGB** 300 ppi | **best in batch — table survives** |
| Wilke & Lee 1955, IEC 47 | *none* | 5 | 5.1 k/pg | CCITT-G4 bilevel 300 ppi | formulae case-mangled |
| Kandula 2010, NASA KSC-2010-007 | `A3.13` 2nd | 15 | 1.7 k/pg | CCITT-G4 bilevel 300 ppi | **worst text layer in batch** |
| Mou et al. 2025, Powder Technol. 453 | `A3.13` | 17 | 3.2 k/pg | born-digital | clean, no render needed |
| Bezzo & Macchietto 2004, CACE 28 | `A2.10` | 13 | 3.8 k/pg | born-digital | clean, no render needed |
| Mills & Chang 2013, CES 90 | *(none — see `A4.8`)* | 7 | 4.2 k/pg | born-digital | clean, no render needed |

### Zwietering 1959 — the tau trap, and it is a new one

Worth its own entry because it fails *silently and numerically*, which is the
dangerous kind. Acrobat 3.0 Capture reads the Greek **τ as the digit 7**,
consistently, throughout the paper. On the page:

> ξ(λ) = τ  (II, 12)   …   τ/2 + (var t)/2τ   …   f(t) = (1/τ) exp(−t/τ)

comes back from the text layer as

> `4th) = 7`   …   `7/2 + (var t)/27`   …   `f(t) = (l/7) exp (- t/l)`

So `27` is *2τ*, not twenty-seven, and `(P/37)` is *(t̄²/3τ)*. Unlike the Robeson
missing-decimal-point case, nothing here looks malformed — every string is a
plausible number. Any parameter taken from this text layer is wrong. **Read the
whole paper on renders.**

Incidental correction found the same way: the article runs to a printed folio
**15**, not 11. It is `Chem. Eng. Sci. 11(1) (1959) 1–15`; the ubiquitous
"1–11" citation is short by four pages.

### Martin & Nilles 1993 — the second German paper in the repo

`A3.15`'s builder is working through Graetz 1882; the same care applies here, and
this file has four specific failure modes:

- **Umlauts and eszett are stripped.** `Wärmeleitung` → `Warmeleitung`,
  `durchströmten` → `durchstromten`, `Schüttungsrohren` → `Schuttungsrohren`,
  `Einfluß` → `EinfluB`. The document is correctly spelled; the OCR is not.
- **German decimal comma.** `L/D = 0,43 bis 2,5` is 0.43 to 2.5. A blind
  comma→point substitution is right here and wrong everywhere else in the repo;
  a blind point assumption reads `0,43` as two values.
- **Superscript powers of ten become the letters `lo`.** `10³ 10⁴ 10⁵` →
  `lo3 lo4 lo5`. Every axis in this paper is logarithmic.
- **Figures are `Abb.`, not `Fig.`** A grep for "Fig" returns nothing and will
  make you conclude the paper has no figures. It has eleven — and no tables at
  all, which is what makes `A3.14` a figure-only case.

### Kandula 2010 — a NASA report number is not a guarantee

`20100036467.pdf` is a genuine 15-page technical report (its own Title metadata
is `KSC-2010-007.pdf`), checked page by page — not an NTRS abstract record and
not a preview. But its Envision OCR is the worst in the batch at 1.7 k
chars/page: `eefcient` for *efficient*, `Zhner-Schlunder` for
*Zehner–Schlünder*, and the same φ = 0.6 appearing as `z, 0.6` and `^ 0.6` two
lines apart. Its Zehner–Schlünder equation comes back as
`ke _ k2/k1 / k1 (1— V) + V42/kt)`. Everything on renders.

### Fuller et al. 1966 — the one good scan, and still check it

The only non-bilevel scan in the drop (JPEG RGB, 8 bpc, 300 ppi), and the atomic
diffusion-volume table survives OCR **with its decimal points intact**:
`16.5 1.98 5.48 5.69 19.5 17.0 −20.2 7.07 6.70 2.88 17.9 16.6 20.1 5.59 16.1
22.8 37.9 18.9`. Equations still degrade (`MOW1'3 + CSX>1/3]2` for the
`(ΣV^⅓ + ΣV^⅓)²` denominator). **Verify the table against the image anyway** —
a table that *looks* right is exactly the Robeson 1991 failure mode, and this is
a 1966 paper with a published erratum that is not on disk.

### Chapman & Cowling — the mid-dot swallows the last row of Table 6

Worth one more line because it caught this pass in the act. Table 6 (book p. 185)
tabulates the collision integrals at sixteen values of *kT*/ε₁₂ ending at
**100·0**. Read off a full-page render at 200 dpi, the mid-dot vanishes and the
last row reads `1000`; the row count went wrong at the same time, and the first
draft of `queue_cases/A4.6.yaml` said "fifteen values … to 1000". Both errors
survived a full-page look and died on a **cropped** re-read at native 300 dpi.

The whole book uses the mid-dot decimal (`2·648`, `0·4291`), and the OCR also
reads leading zeros as the letter `o`. **Crop, go to 300 dpi, and read numerics
twice.** A full-page screenshot is for locating a table, not for transcribing it.
