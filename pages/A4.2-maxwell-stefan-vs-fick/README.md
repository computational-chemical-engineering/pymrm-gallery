# A4.2 — Maxwell–Stefan vs Fick for multicomponent mixtures

Reproduces the bulk-fluid core of Krishna & Wesselingh (1997), *The
Maxwell–Stefan approach to mass transfer*, Chem. Eng. Sci. 52(6) 861–911
([doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7)):
the demonstration that the Maxwell–Stefan formulation and the generalized Fick
matrix `[D] = [B]^-1 [Gamma]` are two bookkeepings of the same physics, and
where the scalar effective-diffusivity simplification breaks (osmotic
diffusion, reverse diffusion, diffusion barrier).

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page: review's two worked matrix examples reproduced
  element by element (including a provable sign typo in the printed B12 of the
  non-ideal example), then the Duncan–Toor two-bulb cell solved under four
  flux closures in one pymrm skeleton — including the raw Maxwell–Stefan
  friction system, which shares no code with the `[B]` route and is what
  actually tests eqs. (26)–(27).
- `data/krishna-wesselingh-1997-worked-examples.csv` — the printed values of
  the two worked examples, transcribed from 600 dpi renders (tier 6, the
  authors' own computations; provenance in the sidecar).
- The experimental comparison loads the digitised Duncan & Toor dataset
  published with page `A4.9` cross-page; nothing was re-digitised.

Key numbers: printed matrices reproduced to ≤1.03% / ≤1.65% per element
(rounding-limited); the independent friction-system route agrees with the
`[B]` route to 7.8e-16 over the transient and 3.1e-12 on `[D](x)` across the
composition triangle; vs experiment 0.59 mole % mean (Maxwell–Stefan) against
3.14 mole % (Wilke with CO₂ eliminated), and 1.38 vs 9.07 mole % on the
nitrogen bulb-difference signal. Every Wilke number depends on which species
is eliminated (mean deviation 1.44–4.31 mole % across the three choices);
Maxwell–Stefan is invariant to 4e-16. Nothing fitted. Runtime ~33 s.

The mean deviation against the measurements is the number page `A4.9` already
publishes — same solver, grid and dataset — so it confirms a faithful port,
not an independent corroboration.
