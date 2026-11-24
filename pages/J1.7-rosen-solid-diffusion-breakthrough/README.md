# J1.7: Rosen (1952), fixed-bed breakthrough with solid diffusion into spheres

From Rosen, J. B., *Kinetics of a fixed bed system for solid diffusion into spherical
particles*, *Journal of Chemical Physics* **20** 387-394 (1952),
[doi:10.1063/1.1700431](https://doi.org/10.1063/1.1700431).

- **Structures:** `S4` (1-D convection with exchange)
- **Runtime:** see `meta.yaml`.

## What the page does

Rosen solved the step response of a linear fixed bed with a liquid film and solid
(homogeneous particle) diffusion into spheres exactly, as a real integral (eq. 25), and gave
a long-bed error-function form (eq. 37) that he called "accurate to 1 percent" for
"γx≧50 and ν≦0.01" (p. 393). The paper does not say what the percent is of. The page
computes the exact curve by two inversions that share only the transform, reproduces
Tables II and IV, and tests eq. 37 over the whole stated range under two readings fixed
before any run: (a) the largest error in concentration, as a fraction of the feed; (b) the
largest relative error in the throughput at which a given fraction breaks through, over
a central band of breakthrough fractions (stated on the page). The catalogue lists the case as pore
diffusion; the paper's model is solid diffusion, and the page uses it.

## Findings

- **Reading (a), concentration: the statement does not hold at its lower end.** The
  largest error of eq. 37 is 0.0128 <!-- claim: r3_a_gx50_nu0 --> of the feed at the
  shortest stated bed without a film and 0.0140 <!-- claim: r3_a_T --> with the largest
  stated film; it falls to the stated "1 percent" at a bed length group of 81.3 <!-- claim: r3_gx_a_nu0 -->
  and 121.9 <!-- claim: r3_gx_a_nu0p01 --> respectively. It exceeds "1 percent" at
  12 <!-- claim: r3_a_nodes_fail --> of the 40 <!-- claim: r3_nodes_big_total --> ladder
  nodes inside the stated range. With the largest film its floor, however long
  the bed, is 0.0059 <!-- claim: r3_floor_nu0p01 -->.
- **Reading (b), throughput: the statement holds.** The largest error is
  0.59 % <!-- claim: r3_b_gx50_nu0 --> and 0.73 % <!-- claim: r3_b_T --> at the same two
  corners, within "1 percent" from a bed length group of 30.9 <!-- claim: r3_gx_b_nu0 --> and
  32.2 <!-- claim: r3_gx_b_nu0p01 -->, and above it at 0 <!-- claim: r3_b_nodes_fail -->
  nodes inside the stated range. On a wider band of breakthrough fractions (on the page) it
  is 1.46 % <!-- claim: r3_b_sens_T --> at the corner with the film, so this verdict
  depends on the band.
- **Cause.** Eq. 37 is eq. 26 with Table II's small-argument forms and no film, an
  identity to round-off. It carries the mean and variance of the curve and not its
  skewness. Adding the film's variance moves the largest error at the corner only to
  0.0127 <!-- claim: r4_filmvar_T -->; adding the third cumulant brings it to
  0.00060 <!-- claim: r4_edgeworth_T -->.
- **Table IV.** The exact column reproduces at all but
  1 <!-- claim: r2_u_cells_off_print --> of the 5 <!-- claim: r2_u_total --> finite entries; at the σθ row printed "1.0" the print "0.808" is
  12.8 <!-- claim: r2_u_st1p0_halfunits --> half-units from 0.814416 <!-- claim: r2_u_st1p0 -->.
  The saddle-point columns mostly do not reproduce (counts on the page), while Rosen's last
  column agrees with his own printed $u_0$ and $\delta_0$ in every row.

These are tests of the paper's formulas against its own exact solution; no measured
breakthrough curve is compared.

## Assumptions

- Linear isotherm, constant D and h, uniform spheres, a step into a clean bed, no axial
  dispersion (the paper's model). The verdicts cover the film groups of the stated range.
- The two readings of "1 percent" and the throughput band are the page's choice, fixed
  before computing.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/rosen-1952-*.csv` | Table II, Table IV and the statements tested, each with its provenance sidecar |
| `data/check_values.csv` | values of the page's functions at one test condition from a second route, asserted before any claim |

There is no `model.py`: the page is closed-form and its functions live in `index.py`.
