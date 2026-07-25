# A4.9 — Duncan-Toor ternary gas diffusion

Two bulbs of gas joined by a capillary. Nitrogen has almost no gradient, so
Fick's law says it should not move. It moves a lot: uphill first, then it stops
against its own gradient, then it reverses.

- **Structures:** `S3` (1D steady BVP), `S9` (implicit multicomponent flux)
- **Reference:** Duncan & Toor (1962), AIChE J 8(1) 38-41, doi:10.1002/aic.690080112
- **Runtime:** ~12 s

## Status

Simulation complete and validated. **Experimental points not yet digitised** —
the CSV is a schema placeholder. See `data/duncan-toor-1962-run1.meta.yaml`.

Model parameters (geometry, initial compositions, MS diffusivities) are the
commonly tabulated values for this experiment and have not been verified against
the original article.
