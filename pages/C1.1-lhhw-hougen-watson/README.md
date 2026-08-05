# C1.1 — Langmuir–Hinshelwood–Hougen–Watson kinetics

**Catalog ID:** `C1.1` · **Section:** C · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1`, `S2` · **Data tier:** 2 (printed experimental tables;
everything is *fit* data — see below)

Reproduces the founding worked example of the LHHW formalism: O. A. Hougen &
K. M. Watson, *Chemical Process Principles, Part Three: Kinetics and
Catalysis*, John Wiley & Sons (1947), Chapter XIX, Illustration 2 (book
pp. 943–958) — the vapor-phase hydrogenation of codimer over supported nickel,
with the complete 18-mechanism parameter estimation on the printed data of
Tschernitz, Bornstein, Beckmann & Hougen (*Trans. AIChE* **42**, 883, 1946,
not consulted; the book reprints the data in full).

The source is a book with no DOI, held as the combined volume (1157 PDF
pages, PDF page = book page + 16), identified from its own title and imprint
pages. Its scan is 150 ppi CCITT-G4 — the lowest native resolution in the
gallery — so every numeral was read from cropped native-resolution
enlargements and the text layer was never trusted for a digit.

## What the page shows

- **The transcription is proven, not asserted.** The book prints its own
  least-squares sums for the 200 °C block; the transcribed columns reproduce
  the three linear sums exactly (1.8e-15) and everything else to the book's
  own printed rounding.
- **The book's estimation reproduces.** Mechanism d's constants at 200 °C by
  two independent routes (transcribed data / the printed normal equations) to
  0.6 %; Σδ² = 3.002 and the famous ±8.44 % both recomputed.
- **A 1947 worksheet slip, diagnosed.** Eight of Table D's eighteen rows
  descend from one corrupted Σ(R·pS) entry — two independent printed rows
  imply the same wrong value to 1.1 %, and substituting that single number
  reproduces both rows to ~1 %. The acceptance verdicts survive it.
- **A misprint in the final equation, adjudicated.** Eq. (q) prints
  ΔS_S = −30.96 where Table G prints −30.46; the book's own K_S row decides
  for Table G at 10× separation. Typing the boxed equation costs ~22 % of K_S.
- **The discrimination question, answered on the founding dataset.** The
  joint sign-test decision {d, h} is threshold-robust and the *family* is
  genuinely selected, but molecular vs atomic hydrogen is not resolvable:
  0.4 percentage points apart in rate space, and even noiseless h-rates are
  fitted by mechanism d to 3.6 %. The book's fit-quality ground for choosing
  d does not survive; its chemical argument does.
- **The fitted law in a pymrm bed.** Isothermal plug flow at 200 °C, 3.5 atm
  (a demonstration — the book poses the problem but prints no answer):
  W/F = 106.6 for 99.8 % hydrogenation, two independent routes agreeing to
  1.6e-5 after Richardson extrapolation, and a measured +51 % bed-length
  penalty for ignoring the 42 % mole contraction.

## Fit vs test

Everything in the source is fitting data — every printed constant was fitted
by the book to the same 40 runs, and no held-out measurement exists. All
agreements on the page are labelled reproduction or internal-consistency
adjudication, never validation. The one synthetic computation (noiseless
h-rates refitted by d, to measure the design's discriminating power) is
labelled as a numerical experiment in the cell that runs it.

## Files

- `build_page.py` — generates `index.ipynb` (run `python build_page.py`, then
  execute the notebook).
- `data/hougen-watson-1947-*.csv` + `.meta.yaml` — the five transcribed
  tables with provenance, including both printings of every constant the book
  prints twice.
- `agreement.json` — 33 metrics for CI regression checking.

Runtime: about 5 seconds end to end.
