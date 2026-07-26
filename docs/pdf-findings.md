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
