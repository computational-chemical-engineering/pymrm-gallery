# B1.3 - Bischoff's generalised modulus: measuring how narrow the narrow band is

Bischoff (1965) builds the effectiveness-factor modulus from the integral of
the rate, which makes the large-modulus ends of the curves for *every*
reaction rate form coincide exactly, and claims the rest of the curves stay
within "about 15%" (orders one-half to three) or "about 30%" (zero order
included). This page measures those claims, computes the LHHW and
volume-change collapses the paper only redrew or tabulated, and audits the
construction out of sample on a finite cylinder against a 1998 competitor.

- **Structure:** `S3` (1D steady BVP; one 2-D audit)
- **Reference:** Bischoff (1965) AIChE J 11(2) 351-355,
  doi:10.1002/aic.690110229; secondary comparison Pan & Zhu (1998)
  Chem Eng Sci 53(5) 933-946, doi:10.1016/S0009-2509(97)00385-0
- **Runtime:** ~1 min

## Results

Both printed claims hold, each under one reading, and the page pins which:
14.61% at m = 1.270 for orders one-half to three (his "about 15%", as
eta_max/eta_min - 1), and 40.55% max/min = 28.85% of the upper curve at the
zero-order kink m = 1 exactly (his "about 30%" is the fraction-of-upper
reading). Extrema are root-found, and both headlines are recomputed on an
independent route (closed forms vs quadrature, gap 5e-7). On the standard
modulus the same curves peak at 120.6% spread and never collapse (41.4%
forever). The LHHW family collapses to a 30.16% band ("again about the same"
is right against the ~30% figure, not the 15%); the volume-change family
stays within 3.18% of the first-order curve. A dead-zone onset formula the
paper does not print, m* = (n+1)/(1-n), is derived from its eq. (14).

Out of sample (finite cylinder Z = 1, second order, exact = 2-D pymrm):
the shape's own first-order curve read at Bischoff's generalised modulus
predicts eta to 0.030 absolute (4.4% relative - inside the advertised band);
reading Figure 1's slab curve costs 0.099 (the B1.2 shape band); Pan & Zhu's
1998 "arbitrary kinetics, < 1.5%" polynomial misses by 0.090, six times its
claimed bound.

## Printed defects found (reported, never repaired)

- **Bischoff Table 1, omega*Co = -0.25 / eq. (33):** prints 0.905; the
  closed-form eq. (33) gives 0.9107 - 0.0057 absolute, 11x the half-ULP
  rounding radius of the printed 3-decimal cell (and the same 11x with both
  sides taken relative). The other seven cells pin the transcription to
  within printed rounding.
- **Bischoff eqs. (27)-(28) (n = 1/2 elliptic reduction) as printed** are
  inconsistent with their own source eq. (20) by 29-41% in m: the where-block
  is reused from the n = 2 case, but the y = sqrt(t) substitution changes the
  prefactor exponent and the phi argument. Amended in those two places
  (labelled inference) the formula agrees to 2e-13 - and the paper's own
  printed 15% shows Figure 1 used the correct reduction, since the as-printed
  curve tops out at eta = 1/sqrt(2) and would give a ~41% band.
- **Bischoff eq. (14), denominator bracket exponent:** prints -1/2 where its
  own parents eqs. (8)+(13) force +1/2 (the numerator's inner bracket
  correctly prints -1/2, isolating the defect to one exponent). As printed
  the equation is dimensionally inconsistent and, in normalised variables,
  low by exactly the factor 1/(n+1) - half the correct m at first order,
  against the paper's own printed eq. (21). The page's Route B implements
  the forced +1/2 reading, labelled as such.
- **Pan & Zhu Table 5, mu2 at Z = 2.5:** prints 2.50; their eq. (30)
  gives 2.40 exactly (likely a copy of the Z = 2.0 cell).
- **Pan & Zhu eq. (36):** prints q3 = -54 - 405 mu1 Q1 + 9 mu2 Q2; the
  matching conditions they state force +405, and the printed sign makes
  their own Table 6 eta_a column negative.

## Data

**Provenance tier 6 - nothing here is a measurement.** All three transcribed
tables are the source authors' own computed values, read from digit-scale
crops of native-resolution renders (300 ppi Bischoff, 600 ppi Pan & Zhu).
The exact eta(m) curves computed by the page are exported as
`eta-generalised-modulus.csv`.
