# B1.1 + B1.5 — Thiele modulus and the Weisz-Hicks non-isothermal pellet

How much of a catalyst pellet actually works. Thiele (1939) answered it exactly
for the isothermal case; Weisz & Hicks (1962) showed that once the pellet can
heat up, the question has up to three answers at once and the pellet can run
more than a hundred times faster than its own surface conditions suggest.

- **Structures:** `S3` (1D steady BVP), `S10` (constrained continuation)
- **References:** Thiele (1939) I&EC 31(7) 916-920, doi:10.1021/ie50355a027;
  Weisz & Hicks (1962) Chem Eng Sci 17(4) 265-275, doi:10.1016/0009-2509(62)85005-2
- **Runtime:** ~5 min

## Results

Isothermal eta matches the exact closed form to 2.2e-4 for phi <= 30 across
slab, cylinder and sphere, with second-order convergence. At phi = 100 the error
grows to 6.5e-3 because the surface layer is thinner than a cell; refining to
n_r = 1600 brings it to 1.9e-5.

Non-isothermal: at beta = 0.6, gamma = 20 the fold spans 0.297 < phi < 0.319,
and at phi = 0.308 there are three steady states with eta = 1.08, 18.3 and 36.1
- a 33x spread at a single Thiele modulus. The pymrm continuation agrees with an
independent shooting solution to 1.0e-5.

## Data

**Provenance tier 6 - not experimental.** Thiele (1939) is analytical and
Weisz-Hicks (1962) is computational; neither reports measurements. Validation is
against exact solutions and an independent shooting reference. Experimental
effectiveness-factor data exists elsewhere and would make a good companion page.
