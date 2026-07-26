# Findings from the priority-1 PDFs

Recorded 2026-07-26 after inspecting the three priority-1 papers. **No numbers
were transcribed into any dataset from these** — see the OCR warning below.

> **Move the PDFs somewhere durable.** They were delivered to a session
> scratchpad, which is not persistent. Suggested: `~/papers/pymrm-gallery/`.
> Still never inside this repository — `.gitignore` blocks `*.pdf` and
> `scripts/check_metadata.py` errors on committed ones.

---

## 1. Xu & Froment (1989) — page `C2.1`

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

## 2. Krishna & Ellenberger (1996) — page `F1.4`

`AIChE_Journal-1996-Krishna.pdf`

**Table 1** experimental setup, operating conditions and system properties —
this is the column-geometry table the page needs. **Table 2** physical
properties of the liquids used. **Table 3** average relative deviation of their
model against the Wilkinson correlation.

**Scale: 2,787 experiments, of which 1,735 in the churn-turbulent regime.** That
is a substantial dataset, and Tables 1–2 are text so the conditions transcribe
cleanly.

The holdup measurements themselves appear to be in figures (e.g. Figure 11,
large-bubble holdup vs *U* − *U*_trans for several gas densities), so expect
digitisation. Their large-bubble model is **Eq. 19**; the transition velocity
uses the **Reilly** correlation.

Good news for the page's framing: they report the Wilkinson correlation
consistently *over*predicting large-bubble holdup while *under*predicting the
small-bubble contribution — a concrete published claim the reproduction can test
rather than merely restate.

## 3. Weisz & Hicks (1962) — page `B1.1` (already published)

`1-s2.0-0009250962850052-main.pdf`

**No tables found in the text layer** — the η(φ) results are presented as
figures only, as expected for a 1962 computational paper.

This is **not blocking**: `B1.1` is already validated against exact isothermal
solutions and an independent shooting reference, which is a stronger test than
agreement with someone else's 1962 numerics. Digitising their curves would add a
third comparison and, more usefully, confirm that our β and γ conventions match
theirs — worth doing but low priority.

---

## Recommended order for the next session

1. **Xu & Froment Table 6** via page-image reading → build `C2.1`. Highest value:
   it is the most-used kinetics in the catalog, and it would be the gallery's
   second page validated against real measurements rather than analytics.
2. **Krishna & Ellenberger Tables 1–2** transcribe (clean text), then digitise
   the holdup figures → build `F1.4`.
3. Weisz & Hicks figures only if time allows.

Note the balance issue this addresses: of three published pages, only `A4.9` is
tested against experiment. `C2.1` and `F1.4` would make it three of five.
