# A4.2 verification — adversarial review

Verifier pass on `queue_cases/A4.2/page/`, 2026-07-31.
Source re-read directly off 600 dpi renders (`pdftoppm -r 600`) of
`Krishna1997-maxwell-stefan-review-CES52-861.pdf`, journal pages 863, 870, 871, 872
(PDF pages 3, 9, 11, 12).

**Verdict: safe to publish after the fixes in F1–F3. Not a send-back.**

---

## Findings, by severity

### F1 — CONFIRMED. Stale number contradicted by the page's own output

`build_page.py:835` ("What pymrm adds") states the Wilke fluxes

> violate the closed-cell equimolar constraint by **~30%** of the largest flux

The Validation cell immediately above prints, from the code:

> Wilke closure defect: |sum_i J_i| reaches **104.8%** of the largest flux

`meta.yaml:35`, `README.md` and `queue_cases/A4.2.yaml:52` all say **105%**.
So the "~30%" in the notebook prose is the only stale copy and it is off by 3.5x.

*Failure scenario.* A reader scrolls from the Validation output (104.8%) to the
claim two cells later (~30%) and finds the page contradicting itself on a
headline number. This is the exact defect class the repo has seen before.

*Fix.* Interpolate the computed value, or set the prose to ~105%.

### F2 — CONFIRMED. The 8.9e-16 "MS = Fick" agreement is algebraically inevitable

`flux_ms` and `flux_fick` both call the **same** `build_b`. The only difference
is `np.linalg.solve` on a 2x2 versus the explicit adjugate/determinant of that
same 2x2. Agreement to 8.9e-16 is a fact about 2x2 linear algebra, not evidence
about eq. (27).

The page presents it otherwise, in three places:

- cell 9: "coded as an independent path so the equivalence is demonstrated, not assumed"
- cell 20 print: "the two formulations are the same physics to machine precision, as eq. (27) claims"
- cell 29: "both are *coded independently* ... and agree to ~10^-15"

*Failure scenario.* A transcription error inside `build_b` — a wrong `d_ms`
index, a dropped minus on the off-diagonal — propagates identically into both
closures, and `ms_fick_equivalence` still reads 8.9e-16. The check has **zero
power against exactly the error class it is presented as guarding**.

*What an independent check gives.* I assembled the full n=3 Maxwell–Stefan
friction system with the `sum_k J_k = 0` bootstrap, never forming `[B]`, and
recovered `[D]` column by column. Over 200 random ternary compositions it
reproduces the page's `fick_matrix` to **2.2e-13**. So the physics and the
closed-form elements are correct — but that is my check, not the page's.

*Fix.* Either relabel honestly ("a consistency check on the closed-form 2x2
inverse; the equivalence itself is algebraic, eqs. 25–27"), or replace
`flux_fick`'s matrix with a path that does not call `build_b` (~15 lines,
verified above to agree).

### F3 — CONFIRMED. The Wilke numbers depend on an unstated modelling choice

`n_i = 2` means the scalar law is imposed on H2 and N2, and **CO2 is the
dependent species that absorbs the closure**. Permuting which species is
eliminated — nothing else changed, same grid, same steps, same data — moves
every headline Wilke number:

| eliminated species | mean over all species | N2 bulb-difference | N2 difference crosses zero? |
|---|---|---|---|
| CO2 (the page's choice) | **3.14** mole % | **9.07** mole % | no |
| H2 | 4.31 mole % | 9.03 mole % | no |
| N2 | 1.44 mole % | 3.39 mole % | **yes**, peaks +0.1238 |
| Maxwell–Stefan (reference) | 0.59 mole % | 1.38 mole % | yes, peaks +0.1227 |

With N2 eliminated the scalar closure crosses zero and reaches +0.1238, which is
within 1% of the Maxwell–Stefan peak of +0.1227. That directly contradicts the
page's flat statements:

- cell 22 print: "no crossing, no barrier, **for ANY positive scalar diffusivities**"
- cell 29: "cannot cross zero for any positive diffusivities"
- title of the figure: "The scalar closure cannot cross zero"

The claim *is* defensible when the scalar law is imposed **on nitrogen** — both
such choices give ~9.05 mole % and no crossing — but the page never says which
species carry the scalar law, and the Reuse section actively invites the reader
to reorder species.

*Failure scenario.* A user adapts the page with a different species ordering
(explicitly offered in Reuse), gets a zero crossing from "the scalar closure",
and concludes the page's central structural claim is false.

*Fix.* (a) state that the scalar law is imposed on H2 and N2 with CO2 by
difference; (b) qualify the "cannot cross zero" claim as being about the species
the scalar law is imposed on; (c) either report the 1.44–4.31 mole % spread or
say the metric is conditioned on the elimination. Worth adding the flip side:
**Maxwell–Stefan gives 0.59 mole % under all three orderings** — it is
elimination-invariant, which is a genuine point in its favour the page misses.

### F4 — CONFIRMED, medium. "Identical to A4.9's published figure" is a tautology

`mean_abs_dev_ms` = `0.0058826829695680805` in A4.2's `agreement.json` is the
**bit-for-bit same double** as A4.9's `mean_abs_deviation`. So are
`mole_balance_error` (2.1632816939198067e-16) and `timestep_sensitivity`
(0.00014848833019553975). Same solver, same grid, same step count, same dataset,
same t=0 compositions — the number could not have come out otherwise.

The page does admit "the transport solver is the A4.9 page's, reused". But
`meta.yaml`'s `agreement:` ("identical to the published A4.9 page's own figure,
with nothing fitted") and the notebook print ("A4.9 reports 0.59 mole % for its
Maxwell-Stefan solve of the same data") both read as independent corroboration.
Fix: say "the same computation, hence necessarily the same number".

### F5 — CONFIRMED, medium. The scalar-shortcut bullet restates a published A4.9 result

A4.9's own `meta.yaml` `adds:` already says: *"wrong by 9.1 mole % in the N2
difference against 1.4 mole % for Maxwell-Stefan"*. A4.2 reports 9.07 vs 1.38 and
presents it under "What pymrm adds".

What A4.2 genuinely adds here is the **route** (a full local-composition Wilke
solve versus A4.9's frozen-D analytic exponential relaxation) plus the
equimolar-closure defect. That the two very different routes land on 9.07 and 9.1
is itself a good cross-check — claim that instead of the number.

### F6 — CONFIRMED, minor. The 105% closure defect is not CI-tracked

`wilke_closure` is printed and quoted in `README.md`, `meta.yaml` and the case
YAML, but is **not** in the `report_agreement` dict. The one number the "adds"
section leans on hardest is the one `check_agreement.py` cannot regress.

### F7 — CONFIRMED, minor, no numerical consequence. Latent bulb-label mismatch

Cell 6 comments `vol_1 = 77.99e-6  # bulb 1 (H2 + N2 charge, Duncan & Toor labelling)`,
but cell 26 passes `xr1` — the **review's** bulb 1, the N2+CO2 charge — as the
z=0 boundary. Volumes do not enter that cell (steady capillary solve only), and I
confirmed the reported J2/J1 is invariant under swapping the two boundary values
(−0.4675 either way). Harmless today; exactly the mismatch that bites on the next
edit.

### F8 — minor wording

"D1_eff ... within a few percent of its pair values (6.8-8.3e-5)": the minimum is
6.432e-5, which is **5.4% below** D_13. Defensible but at the edge of "a few".

---

## What survived attack (recorded so it is not re-litigated)

### The B12 sign typo — CONFIRMED, and provable more strongly than the page argues

I read journal p. 871 at 600 dpi myself. Printed:
`[B] = [[0.363, −0.036],[0.107, 0.495]] x 10^9`,
`[Gamma] = [[0.69, −0.13],[0.07, 1.05]]`,
`[D] = [[1.92, −0.58],[−0.28, 2.25]] x 10^-9`.
The minus glyph on −0.036 is unambiguous and identical in form to the minus on
−0.13 and −0.58 in the adjacent matrices. **Not an artefact, not a hyphen.**

Eq. (26), read at 600 dpi off p. 870, is exactly
`B_ij(i≠j) = −x_i (1/D_ij − 1/D_in)`. With D_12 = 3.4 > D_13 = 2.5 this forces
B_12 > 0.

*Does the "necessarily positive" argument assume the convention it establishes?*
**No.** The other worked example (p. 872, ideal gas) prints `B_12 = +0.007e5`
with the same ordering D_12 = 8.33 > D_13 = 6.8, and I reproduced all four of its
[B] elements to 0.09–3.5% and all four [D] elements to 0.02–1.03%. The two
printed examples are mutually inconsistent under any single sign convention, so
one is a typo — and the ideal one is the one that reproduces.

*Is the 1.65% closure circular — chosen after the correction that makes it work?*
**No.** The decisive route never touches the closure being reported. Back-solving
the unrounded diffusivities from the printed **B11, B21, B22 only** (B12 never
used) gives D = 3.375, 2.507, 1.661 e-9 — within ~1% of the printed 3.4/2.5/1.7 —
and eq. (26) then returns **B12 = +0.0359e9**: the printed magnitude 0.036
exactly, with the opposite sign. The page states this back-solve in prose
(3.38/2.50/1.66) without showing code; I verified it and it is exact.

Nor was the value tuned: +0.036 closes to 1.65%, +0.037 to 1.59%. The rival
hypothesis (B21 is the typo instead, +0.107 -> −0.107) gives **302%** and is
decisively rejected.

| hypothesis | max element deviation of [B]^-1[Gamma] vs printed [D] |
|---|---|
| as printed, B12 = −0.036e9 | 75.06% |
| B12 = +0.036e9 | **1.65%** |
| B12 = +0.037e9 (eq. 26) | 1.59% |
| B21 = −0.107e9 instead | 301.57% |
| both off-diagonals negated | 293.11% |

### Determinism — CONFIRMED clean

Fixed grid, fixed step count, fixed Newton tolerance, no continuation, no warm
start (each `simulate` call re-initialises from a linear profile). Re-executing
the notebook gives **bit-identical** stream outputs to the staged copy, and
`build_page.py` regenerates an identical notebook.

Feature times are robust to refinement in both directions:

| n_steps (n_z=40) | t_O | t_barrier | peak |
|---|---|---|---|
| 200 | 0.96 h | 7.60 h | +0.12240 |
| 400 (page) | 0.96 h | **7.52 h** | +0.12274 |
| 800 | 0.96 h | 7.52 h | +0.12291 |
| 1600 | 0.96 h | 7.51 h | +0.12300 |
| 3200 | 0.96 h | 7.51 h | +0.12304 |

n_z = 20 / 40 / 80 / 160 at n_steps = 400: t_O = 0.96 h, t_barrier = 7.52 h,
peak = +0.12274…5, mean deviation 0.588 mole % — identical to three figures.
Quoting 0.96 h and 7.52 h is defensible (barrier converges to 7.51 h, 0.1% away).

### Bulb labelling — CONFIRMED correct, no sign flip anywhere

Journal p. 863 at 600 dpi: eq. (3) reads
`Bulb 1: x1 = 0.00000, x2 = 0.50086, x3 = 0.49914` /
`Bulb 2: x1 = 0.50121, x2 = 0.49879, x3 = 0.00000`, transcribed exactly into the
CSV, and "The two bulbs were connected by means of a **86 mm** long capillary
tube". The next paragraph confirms the direction: "hydrogen diffuses from bulb 2
to bulb 1 ... Carbon dioxide diffuses from bulb 1 to bulb 2" — so the review's
bulb 1 **is** the N2+CO2 charge, opposite to Duncan & Toor. The page's
experimental figure labels bulb 1 "charged H2 + N2", matching the A4.9 dataset it
plots. Internally consistent; the uphill-diffusion signal is not flipped.

### The 6.7% frozen-matrix claim — survived a confounding attack

I suspected the 6.7% was a bookkeeping artefact: the review evaluates at
Δx = (−0.25, 0) (bulb-to-equilibrium, Δx2 exactly 0) while the full solve uses
bulb-to-bulb Δx = (0.50121, −0.00207). Evaluating the **frozen** printed [D] on
the full bulb-to-bulb driving force gives J2/J1 = −0.5001, against the review's
−0.4990 and the full solve's −0.4675. So the 6.7% really is the composition
dependence of [D] along the capillary, as claimed. The ratio is also invariant
under swapping the bulbs (−0.4675 both ways), so it is convention-free as stated.

### Traceability sweep

Every number I could find in the prose traces to a printed source, to code on the
page, or to a labelled reconstruction — **except** the "~30%" of F1. Spot-checked
and verified: `|D21/D22| ≈ 1.8` (1.774); the 3.38/2.50/1.66 back-solve (exact);
"experimental error 2.6 mole %" (traces to A4.9's `meta.yaml`); the ideal-example
matrices and eq. (3) compositions (read off the renders myself, CSV matches).

### Housekeeping

- Notebook executes clean in **15.4 s** against `runtime_seconds: 16`; all three
  closures converge; no errors.
- Nine required sections present in the required order; Colab cell is cell 1; no
  Quarto markdown; `load_data(..., page=...)` used throughout; `report_agreement`
  called.
- `models_entry.yaml` present, `slug`/`title` match `meta.yaml`, upgrade-in-place
  note included. Data tier 6 declared correctly; the tier-4 experimental set is
  loaded cross-page with provenance left in A4.9's sidecar, as claimed.
- `review/` holds no page images. Nothing was digitised for this page.

### Caveat honesty

Good overall. The page does not call reproduction "validation" — it says the
worked examples are the authors' own computations, tier 6, and separates them
from the tier-4 measurements. One thing it could add: the Maxwell–Stefan
0.59 mole % sits **at** the digitisation floor (A4.9's sidecar gives ±0.005 mole
fraction = 0.5 mole %), so that agreement is resolution-limited and should not be
read as a tighter result than it is. A4.9 says this; A4.2 does not repeat it.

---

## Fix list before publishing

1. **F1** — replace "~30%" with the computed value (105%).
2. **F2** — stop calling the 8.9e-16 check an independent implementation, or make
   it one.
3. **F3** — state the elimination choice; qualify "cannot cross zero for ANY
   positive scalar diffusivities"; note the 1.44–4.31 mole % spread (and that
   Maxwell–Stefan is elimination-invariant).
4. **F4/F5** — reword the A4.9 comparisons: same computation, not corroboration;
   the added value is the route, not the 9.07/1.38 numbers.
5. **F6** — add `wilke_closure` to `report_agreement`.
6. **F7** — fix the `vol_1` comment or the `xr1` usage so the labelling is
   consistent.
