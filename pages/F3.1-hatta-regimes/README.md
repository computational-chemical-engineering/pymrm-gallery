# F3.1 — Hatta number and the gas-liquid reaction regimes

Whether a dissolving gas reacts in the film or in the bulk decides which
absorber you should build. One number settles it.

- **Structure:** `S3` (1D steady BVP)
- **References:** Hatta (1932); Van Krevelen & Hoftijzer (1948)
  doi:10.1002/recl.19480670708; DeCoursey (1974) doi:10.1016/0009-2509(74)85003-7
- **Runtime:** ~30 s

## Results

Pseudo-first-order matches the exact `Ha/tanh(Ha)` to 6.3e-3 over
0.02 <= Ha <= 200, with the error concentrated at high Ha where the reaction
layer is thinner than a cell. The two-component solve reduces to the first-order
limit to 5.4e-6 when the co-reactant is in large excess, and reaches the
instantaneous ceiling E -> E_i at large Ha.

Accuracy of the standard approximations against the full film solution:

| approximation | max error | worst near |
|---|---|---|
| Van Krevelen-Hoftijzer | 2.1% | Ha ~ 9 |
| DeCoursey             | 8.7% | Ha ~ 1.5 |

DeCoursey is worst in the transition between the fast and instantaneous
regimes - the operating region of a real amine absorber.

## Data

**Provenance tier 6 - not experimental.** The sources are analytical: Hatta's
regimes and the first-order enhancement are exact, and VKH and DeCoursey are
published approximations. Validation is against exact solutions and limiting
cases.
