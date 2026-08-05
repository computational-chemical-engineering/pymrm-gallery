# A3.3 — Danckwerts surface renewal

Replace the stagnant interfacial film with a surface that is continually torn up
and replaced, assume the chance of an element being replaced does not depend on
its age, and the age distribution is forced: `φ(θ) = s e^{−sθ}`. Averaging the
penetration flux over it gives `k_L = √(Ds)`, against film theory's `k_L ∝ D`.

**Source.** P. V. Danckwerts, *Significance of liquid-film coefficients in gas
absorption*, Ind. Eng. Chem. **43**(6) 1460–1467 (1951),
[doi:10.1021/ie50498a055](https://doi.org/10.1021/ie50498a055). Read from the ACS
PDF at its **native 300 ppi** (`pdfimages -list`: CCITT-G4 bilevel, 300 ppi —
rendering higher is interpolation), with the two constants the text layer
mangles re-read on nearest-neighbour crops at that resolution. The PDF is not in
this repository.

## What is on the page

- **Every number the paper prints** — all six reproduce, and three of them are
  loose enough to be worth stating rather than rounding away (the depth at 1 hour
  is 6.83 mm against his "about 6 mm").
- **A pymrm transient element** (`S4`) marched on an age grid uniform in
  `w = √θ`, age-averaged by his eq. 7, reproducing eqs. 8, 10/30 and 12/34 to
  ~4e-5. Both axes refined: observed orders 0.99 in the age step and 1.98 in the
  grid, with the time error the larger of the two — which is why **every**
  agreement on the page is Richardson-extrapolated in the age step, with the raw
  value printed beside it.
- **A printed defect in eq. 25**, proved from the paper's own eq. 9 and eq. 29:
  the exponent needs `H²`. The proof is that with `H²` the eq. 25 integral
  reproduces eq. 9 at every `(H, k_G/H)` swept (worst 4.1e-12) and as printed it
  reproduces it nowhere (closest 4.2 % off). The printed form also *diverges*,
  but only where `k²(H−1) > sD` — dramatic at `H = 5, k_G/H = 0.7`, and at
  `k_G/H = 0.3` the same expression converges to a finite wrong number.
- **The case he says he could not solve** — the second-order reaction, solved
  with `NumJac` + `newton`, which turns his 10 % error bound from an assumption
  into a measurement (+0.22 % at his own illustrative point, and conservative by
  a factor of 53 in the concentration ratio it demands). A **pure-discretisation
  control** on the same solver, at `c₀'/c* = 1e6` where the answer is exactly
  `√2`, says how much of that is the age grid: +0.1006 % raw at `n_t = 400`,
  +0.0104 % extrapolated.
- **Film theory and surface renewal as two ends of one calculation**:
  `k_L = √(Ds) coth(d√(s/D))`, whose exponent `n = ½ + z/sinh 2z` runs from 1 to
  ½ as the layer deepens.
- **Why the √D exponent cannot be tested from this paper**, stated as the page's
  main negative result.

## What it does not do

The paper has no data. `s` comes from inverting the very relation one would want
to test — an inference, flagged as one, since the paper says only "calculated
from published values of `k_L`" — so film theory and surface renewal fit its one
physical number equally well, with zero residual each. Nothing here is validated
against a measurement. Whitman (`A3.1`) and Higbie (`A3.2`) are separate cases
and are not built or spoken for here.

## Run it

```bash
python build_page.py          # regenerate index.ipynb from the builder script
```

Runtime about 60 s. Data: `data/danckwerts-1951-printed-numbers.csv` with its
provenance sidecar; no other page's dataset is loaded.
