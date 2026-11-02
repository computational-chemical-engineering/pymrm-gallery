# B1.4 - The Weisz-Prater criterion

**Catalog ID:** `B1.4` · **Section:** B · **Tier:** T0 · **Structures:** `S3`
· **Provenance tier:** 6 (no measurement anywhere in either source)

The inverse of every other pellet page in the gallery. `B1.1`/`B1.5` start from a
known rate constant and compute the effectiveness factor; the Weisz-Prater
criterion asks what an experimenter can conclude from **observables alone**
(measured rate, particle size, effective diffusivity, surface concentration)
without knowing the kinetics:

```
Phi = (dN/dt) * (1/c_0) * (R^2/D)  <  1
```

- **Source:** `index.py` (paired with `index.ipynb`); the rate laws, the pymrm
  pellet and the shooting branch tracer in `model.py`, tested by `test_model.py`
- **Runtime:** about 360 s

## Provenance: the reprint route

**Weisz & Prater (1954), *Adv. Catal.* 6, 143 is not on disk** and has no
open-access route; it is a pre-DOI book chapter. It is cited as the origin, from
Weisz & Hicks' own reference list, and **was not consulted**.

Everything transcribed is from **Weisz & Hicks (1962)**, *Chem. Eng. Sci.*
**17**(4), 265-275, `doi:10.1016/0009-2509(62)85005-2`, which is also the source
of `B1.1`/`B1.5` and `B1.6`. It prints the criterion as **eq. (1)**, states its
purpose, gives the observable group as **eq. (11)**, extends it to the
exothermic case as **eq. (10a)**, devotes its **Section V** to the observable
formulation, and tests it against its own Fig. 2.

Weisz & Hicks attach eq. (1) not to their reference [5] (Weisz & Prater 1954)
but to **[9] Weisz, *Z. Phys. Chem.* 1957 **11** 1**. Reference [5] is named and
used in the same paper, for the shape-insensitivity of eta(phi) and for the
printed Q'/Q relation. Both are recorded; neither was consulted.

Every equation was read off 600 dpi renders (journal pages 265, 267, 268, 269,
272, 273); both references come from the literature-cited block on page 274,
which prints **no chapter title** for [5]. The scan's text layer is unusable for
equations: eq. (1) OCRs as `$:$<I` and eq. (10a) as `expCrBK1 + LOI`.

## Results

1. **What Phi < 1 guarantees**, computed exactly: eta at Phi = 1 is
   0.9353 <!-- claim: eta_at_Phi1_sphere_first_order --> for a first-order
   sphere and 0.6948 <!-- claim: eta_at_Phi1_slab_first_order --> for a slab;
   across geometries the spread is 0.240 <!-- claim: eta_at_Phi1_spread_on_R -->
   on R and 0.132 <!-- claim: eta_at_Phi1_spread_on_VS --> on V/S. Across
   kinetics on a sphere it runs from exactly 1 (zero order) through
   0.8796 <!-- claim: eta_at_Phi1_second_order --> (second order) to
   1.0483 <!-- claim: eta_at_Phi1_inhibited_LH --> (self-inhibited
   Langmuir-Hinshelwood).
2. **A false-negative map** over beta and gamma: for an exothermic pellet
   eq. (1) admits eta up to 3136 <!-- claim: eta_at_Phi1_worst_exothermic -->;
   eq. (10a) holds eta within 0.0134 <!-- claim: eq10a_worst_eta_deviation_exothermic -->
   of 1 at its threshold, at a conservatism up to a factor
   1.01e5 <!-- claim: eq10a_max_conservatism --> in Phi. For an endothermic
   pellet eq. (10a) inverts: its threshold is reached only at eta of
   0.290 <!-- claim: eq10a_endothermic_eta_at_threshold_max --> or less.
3. **The paper's Section V claim** that the observable "has removed the
   multiplicity" is tested: it postpones it, from beta*gamma =
   4.91 <!-- claim: multiplicity_onset_betagamma_min --> to
   5.56 <!-- claim: multiplicity_onset_betagamma_max --> in phi_0 to
   10.84 <!-- claim: observable_fold_onset_betagamma_min --> to
   21.74 <!-- claim: observable_fold_onset_betagamma_max --> in Phi (each onset
   located as the cusp of the branch). Every fold sits below Phi = 1 (largest
   turning point 0.8055 <!-- claim: observable_fold_turning_point_Phi_max -->),
   inside the band the criterion certifies as safe. One of the four panels of
   Fig. 7 draws a curve in the folded regime; that fold spans
   6.9 <!-- claim: fig7d_fold_width_px_at_600dpi --> px on a 600 dpi render
   against a 6 to 7 px printed stroke, so the figure is silent: a limit on the
   claim, not a refutation of the paper.

## Validation

No experiment (tier 6). Ranked: the shooting reference against the exact
isothermal eta, below 1e-9 <!-- claim: shooting_vs_exact_isothermal_eta -->;
pymrm against the same closed form,
8.3e-6 <!-- claim: pymrm_vs_exact_isothermal_eta_n400 --> on 400 cells at the
observed order of every case; pymrm against the shooting solution over
non-isothermal branch states including unstable ones,
5.0e-5 <!-- claim: pymrm_vs_shooting_nonisothermal_n800 --> on 800 cells; the
threshold found by a root scan on pymrm and on the closed form; three printed
statements reproduced and used as inputs to nothing; a break table with one
injected defect per row. The flux/volume form of eq. (11) is reported as what it
is, the volume-weighted residual of the solve.

## Files

```
index.py         the page source (jupytext percent format)
index.ipynb      paired with index.py
model.py         rate laws, the pymrm pellet, the shooting branch tracer
test_model.py    import, fast solves, two stored claims
agreement.json   written by the Claims cell
meta.yaml        page metadata
```

No `data/`: nothing is digitised and no CSV is shipped.

## Corrections

- 2026-11-01: said the apparent activation energy Q'/Q overshoots 1 by 1.2e-5
  at small phi_0 because of the truncation error of the central difference,
  and called the curve monotone decreasing. The overshoot was round-off: the
  sphere's closed form cancels in coth(phi) - 1/phi at small phi, and the
  difference quotient magnifies it; with a series there, Q'/Q is
  1 - phi_0^2/15, below 1 and strictly decreasing, which the old values were
  not (they were monotone only within a tolerance of 1e-4).
- 2026-11-01: said eq. (10a) applied to an endothermic pellet certifies every
  attainable state, and that eq. (1) there already admits a depression of more
  than 50 %. Accepting every measurement passes all
  384 <!-- claim: not_clean_endothermic --> unclean endothermic measurements
  and eq. (10a) passes 369 <!-- claim: fn_eq10a_endothermic -->, so it is
  close to doing nothing; its threshold is always reached, at eta of
  0.290 <!-- claim: eq10a_endothermic_eta_at_threshold_max --> or less (the
  traced branch only stopped before it, in
  5 <!-- claim: endothermic_cases_beyond_ladder --> cases); and eq. (1) admits
  a depression above 50 % in only
  1 <!-- claim: endothermic_beyond_ladder_eta_below_half --> of those five
  cases (eta 0.4838 <!-- claim: eta_at_Phi1_worst_endothermic -->).
- 2026-11-01: said the fold onsets were upper bounds (beta*gamma 4.86 to 5.61
  in phi_0, 10.9 to 22.7 in Phi). They came from a fold-depth threshold on a
  coarse trace, bisected to 4e-3 in beta. Located as the cusp, they are
  4.91 <!-- claim: multiplicity_onset_betagamma_min --> to
  5.56 <!-- claim: multiplicity_onset_betagamma_max --> and
  10.84 <!-- claim: observable_fold_onset_betagamma_min --> to
  21.74 <!-- claim: observable_fold_onset_betagamma_max -->; the published
  phi_0 onsets at gamma = 30 and 40 were below the onset, not above it (the
  old Phi onsets were above the cusp).
- 2026-11-01: said the worst ambiguous measurement fits eta = 2.28, 5.09 or
  226. Those were interpolated along the traced branch; root-found on the
  shooting solution the three states spread over a factor
  100.7 <!-- claim: ambiguous_measurement_worst_eta_spread -->, and the
  smallest of them is 2.2 <!-- claim: ambiguous_measurement_worst_eta_min -->,
  not 2.3.
- 2026-11-01: the confusion counts (eq. (1) passes 261 unclean exothermic
  measurements, eq. (10a) rejects 397 clean ones) classified each state by eta interpolated along the traced
  branch, which put 14 <!-- claim: interpolated_class_flips --> measurements on
  the wrong side of the 5 % line. With
  each straddling state root-found the counts are
  252 <!-- claim: fn_eq1_exothermic --> and
  406 <!-- claim: fp_eq10a_exothermic -->.
