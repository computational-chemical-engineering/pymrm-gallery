# B1.2: Aris generalised shape modulus

From Aris, R., *On shape factors for irregular particles. I. The steady state
problem. Diffusion and reaction*, Chemical Engineering Science **6**(6) 262-268
(1957), doi:10.1016/0009-2509(57)85028-3.

- **Structures:** `S3` (1D steady BVP), `S13` (multi-dimensional geometry)
- **Source:** `index.py` (paired with `index.ipynb`); the solver in `model.py`,
  tested by `test_model.py`
- **Runtime:** about 55 s

One review blocker is open; see below and `meta.yaml`'s `blocker` block.

## Results

- Eqs. 15-17 reproduce 19 of the 21 Table 1 entries to
  0.00078 <!-- claim: table1_worst_19_entries -->, within one unit of the
  printed third decimal. The other two are the cylinder column at Lambda = 5 and
  10, printed 3.8 <!-- claim: cyl_printed_high_pct_Lambda5 --> % and
  2.6 <!-- claim: cyl_printed_high_pct_Lambda10 --> % above eq. 16, while the
  plate and sphere columns are right at those moduli to
  0.0003 <!-- claim: table1_plate_sphere_worst_Lambda_5_10 -->. No single
  rescaling of Lambda fits both cells (0.961 <!-- claim: cyl_printed_needs_Lambda_factor_5 -->x
  and 0.974 <!-- claim: cyl_printed_needs_Lambda_factor_10 -->x), so they are
  reported as slips, not corrected.
- The spherical-shell expression reproduces his five printed values to
  0.00041 <!-- claim: shell_formula_vs_printed -->, approaches eq. 17 as p^3 and
  eq. 15 linearly in 1 - p, and gives 0.7227 <!-- claim: shell_formula_p075 -->
  at p = 0.75, which decides between the two printings (1995: 0.723, 1957 OCR: 728).
- Aris's case (iv) series is right with p read as radius:length (his
  nomenclature): the 2-D solve matches it to
  4.2e-4 <!-- claim: fincyl_series_p_as_nomenclature --> on 160 radial cells,
  and that gap is grid error (the ladder's Richardson limit lies below
  1e-5 <!-- claim: fincyl_extrapolated_vs_series --> from the series). Read as
  in his running text, the series is off by
  47 % <!-- claim: fincyl_series_p_as_running_text --> of the solved value.
- Over fourteen shapes the band is
  0.1165 <!-- claim: band_all_shapes --> wide at its widest, against
  0.0918 <!-- claim: band_aris_plate_sphere_max --> for plate to sphere,
  1.27 <!-- claim: band_ratio_all_over_aris -->x. The plate is its upper edge and
  the cube its lower edge at all 28 <!-- claim: sweep_points_cube_lowest -->
  swept moduli.
- **The sphere is not the floor.** At the same v_p/s_x = a/3 the cube lies up to
  0.0261 <!-- claim: cube_below_sphere --> below it (3.6 <!-- claim: cube_below_sphere_pct_there --> %
  there), and up to 3.9 <!-- claim: cube_below_sphere_max_pct --> % in relative
  terms; the cylinder of length 2a lies between them at all
  28 <!-- claim: sweep_points_cyl2a_between_sphere_cube --> swept moduli. At equal volume the sphere
  is still the worst case (Aris's Leva argument); both hold.
- Eq. 18, eta Lambda -> 1, holds for five shapes including two with edges: the
  references extrapolate to 1 within
  1.2e-4 <!-- claim: eq18_reference_extrapolated_worst_dev --> (the line fit; the
  cylinder reference, a closed-form axial sum, is truncated by
  5e-8 <!-- claim: eq18_cyl_reference_truncation_Lambda50 --> at Lambda = 50),
  the pymrm solves within 0.0011 <!-- claim: eq18_extrapolated_worst_dev -->,
  mostly their grid error at Lambda = 50.
- Eqs. 19-20 reproduce his worked mixture example (mean modulus and eta at it);
  eq. 19 gives 0.4936 <!-- claim: mixture_eta_eq19 --> for the third number,
  printed 0.397 in the 1995 reprint and 0.897 in the 1957 scan.

pymrm against the references: 1-D to
2.4e-4 <!-- claim: pymrm_1d_worst_rel --> (n = 400, Lambda up to 20), square rod
8.2e-6 <!-- claim: square_rod_vs_series -->, cube on 32^3
7.5e-5 <!-- claim: cube_vs_series_n32 -->, the sweep's rod and finite cylinders
(graded grids) within 3.9e-3 <!-- claim: sweep_pymrm_vs_series_worst --> up to
Lambda = 20. The cube's curve in the sweep is its eigenfunction series at every
Lambda, checked against the 3-D solve at Lambda = 1 only: on 24^3 and 32^3
uniform grids the 3-D solve differs by 11 % <!-- claim: cube_n24_n32_gap_Lambda20 -->
at Lambda = 20.

## Data

Three small tables of the author's own computed values (Table 1, the section 3
shell table, the section 4 mixture example), each with a provenance sidecar.
**Provenance tier 6, not experimental**: Aris (1957) reports no measurements.

## Getting the paper: read the 1995 reprint, not the 1957 scan

Aris (1957) is closed access. It was retrieved through the Elsevier article
retrieval API (PII `0009250957850283`) under an institutional subscription. That
scan's text layer is the publisher's OCR of a mid-dot-decimal typesetting and is
useless for numbers: Table 1's 0.698 comes out as `898`, 0.432 as `482`.

**Chemical Engineering Science reprinted the paper verbatim in 1995** for its
Golden Jubilee, 50(24) 3899-3903, doi:`10.1016/0009-2509(96)81819-7`, PII
`0009250996818197`, and that re-typesetting OCRs cleanly. Both were retrieved
and every transcribed number is one on which the two agree, or the disagreement
is stated. **This is worth trying for any pre-1970 CES classic.**

## What could NOT be obtained

The entitlement returns a **one-page preview PDF only**, for both printings
(journal pages 262 and 3899). Table 1 is on page 265 / 3901, and the
article-object endpoint holds no page images for either PII. The 1957 PDF in
`~/papers/pymrm-gallery/` (`Aris1957-shape-generalised-modulus-CES6-262-PREVIEW1P.pdf`)
is that same one-page preview. The blocker needs a full PDF from the maintainer.

## The sign trap

The low end of every axis carries `{"a": 1, "b": 0, "d": 0}`: with the
**outward** normal that is dchi/dn = 0, symmetry for a solid shape and a
*sealed* face for the shell. Mark the shell's inner face exposed by mistake and
you get a perfectly convergent solution to a different problem, with no warning.

## Relationship to B1.1

`B1.1` plots the Aris collapse for Thiele's three geometries. This page is about
the shapes that have no closed form, and cross-links to `B1.1` rather than
restating it.

## Corrections

- 2026-11-01: the notebook said that reading p as in the running text puts the
  case (iv) series out by up to 40 %; it is
  47 % <!-- claim: fincyl_series_p_as_running_text --> of the solved value (as
  this README already said).
- 2026-11-01: said that every shape in the sweep other than the cube is solved
  with pymrm on a uniform grid, within 3e-3 of an exact reference at
  Lambda = 20. The plate, long cylinder, sphere and shells are closed forms, and
  on the uniform grids the finite cylinders were off by up to
  3.3 <!-- claim: sweep_disc_uniform_grid_worst_pct --> % at Lambda = 20 (the
  flat disc, l/a = 0.25). On grids graded toward
  the exposed faces they are within
  3.9e-3 <!-- claim: sweep_pymrm_vs_series_worst -->, and the ordering is read
  from the references.
- 2026-11-01: said eqs. 15-17 reproduce 19 of the 21 Table 1 entries to the
  printed rounding. Two of the nineteen (cylinder at Lambda = 0.5, sphere at
  Lambda = 2) are one unit off in the third decimal; 17 round to the printed
  digits.
- 2026-11-01: said the cube converges at second order, quoting an error ratio of
  2.29 for a 1.33x refinement. That ratio is order 2.9; from successive
  differences on 8^3, 16^3 and 32^3 the order is between 2 and 3
  (2.8 <!-- claim: cube_order_8_16_32 -->), the grids being pre-asymptotic.
- 2026-11-01: said that beyond Lambda of about 5 the cube's curve is its
  eigenfunction series rather than a pymrm solve. The sweep takes the cube from
  its series at every Lambda; the 3-D solve checks that series at Lambda = 1
  only.
- 2026-11-01: said the p -> 0 and p -> 1 limits are what catch a shell whose
  inner face is exposed by mistake. Those limits were checked on the formula,
  not on the solver, and p -> 0 cannot catch it: the exposed shell also tends to
  the sphere as p -> 0. The page's check is the solve against the shell
  expression at the printed p, which the exposed inner face moves by
  0.16 <!-- claim: shell_exposed_inner_face_dev --> (break table).

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        Pellet and the shape catalogue
test_model.py   fast checks of model.py and one stored claim
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
data/           Aris's printed tables, each with a provenance sidecar
```
