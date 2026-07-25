# A4.9 — Duncan-Toor ternary gas diffusion

Two bulbs of gas joined by a capillary. Nitrogen starts 0.06 mole fraction apart,
closes that gap, then keeps going *against* its own gradient until it is twice as
far apart the other way, and finally stops dead at the moment its gradient is
largest. Fick's law can produce none of that.

- **Structures:** `S3` (1D steady BVP), `S9` (implicit multicomponent flux)
- **Reference:** Duncan & Toor (1962), AIChE J 8(1) 38-41, doi:10.1002/aic.690080112
- **Runtime:** ~11 s

## Agreement

0.59 mole % mean absolute deviation over 28 digitised points, nothing fitted.
The paper reports 0.45 mole % for its own Maxwell-Stefan predictions and quotes
2.6 mole % experimental error; about 0.5 mole % of our figure is digitisation
error rather than model error.

Point O (osmotic diffusion) at 0.96 h, diffusion barrier at 7.52 h.

## Data

Digitised from Figure 2 by programmatic marker extraction, validated by closure,
monotonicity and the closed-cell exchange balance. See
`data/duncan-toor-1962-run1.meta.yaml` for the method and its limits.

All model parameters were verified against the paper text. Note the paper uses an
*effective* (L/A) of 258.1 1/cm measured from binary runs, not the geometric
251.1 1/cm.
