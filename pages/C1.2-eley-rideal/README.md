# C1.2 — Eley-Rideal kinetics: the impact rate law, and the case that the name is wrong

**Catalog ID:** `C1.2` · **Section:** C · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1` · **Data tier:** 6 (analytical; the source prints no
data — see below)

Built from Prins, R., "Eley-Rideal, the Other Mechanism", *Topics in
Catalysis* **61** (2018) 714-721, doi:10.1007/s11244-018-0948-8 — a
historical and terminological review whose thesis is that the mechanism
universally called **Eley-Rideal** (a gas molecule colliding with a
chemisorbed one) was proposed by **Langmuir in 1922** and should be called
**Langmuir-Rideal**, while what Eley and Rideal actually studied (1940-41)
is a chemisorbed atom reacting with a *physisorbed* molecule. The page
reports that finding plainly (same class as `B1.4` and `H1.1`), builds the
impact rate law, and contrasts it with the Langmuir-Hinshelwood law
**imported from the published page `C1.1`** rather than re-derived.

## What the page shows

- **Every printed limit of Eqs. (3a)-(3c), derived.** All eleven
  limiting-order claims of p. 719 hold; derived symbolically and measured
  numerically on independently coded laws.
- **One printed defect, reported and not repaired.** "The rate of the L-H
  reaction should even reach a maximum for c_A = K_A." [sic] is
  dimensionally inconsistent with the paper's own Eq. (3a); the derived
  condition K_A·c_A = 1 + K_B·c_B is unit-invariant where the printed one
  moves by λ² under a unit rescale. The maximum is root-found by two routes
  sharing no algebra (4.8e-11 apart).
- **The p. 719 paragraph turned into a number.** The identifiability window:
  the chemisorbing species must span a factor w* = 11.8 (root-found, two
  optimizer routes to 1.8e-11) centred on the L-H maximum before the best
  impact-law fit misses noiseless L-H rates by more than the 1947 fit's own
  ±8.44 % — and centring the window where both laws are first order
  collapses the misfit by three orders: width buys nothing without
  curvature.
- **The impact law on the founding LHHW dataset (numerical experiment).**
  Hougen & Watson's own mechanism (o) *is* Prins's Eq. (3b) (with the
  product's adsorption term retained). Refit to noiseless mechanism-(d)
  rates over the 1947 design it misses by 27.4 % mean (48 % worst) — the
  mean-|Δ| at the least-squares-optimal fit, the declared protocol;
  directly minimising mean-|Δ| instead gives a smaller number, conclusion
  unchanged — against a 16.8 % noise floor — fit quality alone could have
  rejected it. But allow a hydrogen term in the first-power denominator
  (3.3 %, nonzero) or collapse the design's decade-wide hydrogen span
  (3.2 %) and it drops below the noise floor: the founding design
  discriminates only because it spans a decade in p_H (the design also
  straddles the L-H turnover, p_U* = 2.50 atm, inside the design range — an
  observation, not a second necessary condition: collapsing the hydrogen
  span alone already suffices). The denominator exponent (1 vs 2) is
  indistinguishable there at the design's noise floor even when the rates
  carry no noise.
- **Cross-page reconciliation that caught what it should.** The imported
  constants reproduce the book's printed Table C to 0.006 worst — except run
  25a at 0.503, independently reproducing `C1.1`'s documented digit slip
  (printed 5.40, row-consistent 5.90).

## No data, and what that means

The source prints **no measurement, no table and no figure** (checked by
full-text read, page-by-page render inspection, and the pdfimages
inventory). The page therefore claims **no experimental validation**:
everything is analytical or a labelled numerical experiment on `C1.1`'s
fitted constants — which are themselves fit data, as `C1.1`'s own metadata
records.

## Files

- `build_page.py` — generates `index.ipynb` (run `python build_page.py`,
  then execute the notebook).
- `data/prins-2018-printed-statements.csv` + `.meta.yaml` — every
  quantitative statement the paper prints, with page provenance (all
  numerics verified on 300 dpi render crops).
- Three `C1.1` datasets are loaded cross-page; that page's findings on the
  exact rows used are listed in the notebook's data section.
- `agreement.json` — 46 metrics for CI regression checking (two named
  structural identities sit below CI's ABS_FLOOR and carry above-floor
  companions).

Runtime: about 15 seconds end to end.
