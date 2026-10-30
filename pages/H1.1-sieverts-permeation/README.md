# H1.1 — Sieverts-law hydrogen permeation through Pd membranes

The permeation closure of Itoh (1987), *A membrane reactor using palladium*,
AIChE Journal **33**(9) 1576–1578,
[doi:10.1002/aic.690330921](https://doi.org/10.1002/aic.690330921): his Eq. 1
and the constants beneath it.

This is the *closure*, not the reactor. The reactor built on it is
[`H1.4`](../H1.4-itoh-membrane-dehydrogenation/), from the same paper.

## What the source actually says

- Itoh names the law and uses it: hydrogen permeation "was assumed to obey the
  **half-power pressure law** (Bohmholdt and Wicke, 1967), i.e., it is
  proportional to the difference between **the roots** of the hydrogen partial
  pressures". Eq. 1 is
  `Q_H = alpha_H (sqrt(p_H/P_o) - sqrt(p'_H/P_o))`,
  `alpha_H = 2 pi l_o / ln(r_o/r_i) * D * C_o`.
- **The exponent is attributed to Bohmholdt & Wicke, not to Sieverts.**
  Sieverts & Danz (1936) is cited for exactly one number:
  `C_0 = 1280 <!-- claim: C_0 --> mol/m3`, the solubility of hydrogen in
  palladium at 473 K in equilibrium with 1.013e5 Pa, the amplitude of
  `c = C_0 sqrt(p/P_0)`.
- **No permeation measurement is printed.** No flux, no measured permeance, no
  pressure sweep. `alpha_H = 4.47e-5 <!-- claim: alpha_H_printed --> mol/s` is
  *computed* by Itoh from `D` and `C_0`, and it is the only external number this
  page can be checked against.
- **No temperature dependence is printed.** `D` and `C_0` are quoted at 473 K
  only, so the page says nothing about the activation energy of the permeance.
- **No misprint was found anywhere in the permeation law.** The two printing
  problems in this paper (Eq. 5's sign and the `K_P` expression) are elsewhere
  and are resolved on `H1.4`. The scan is a 300 dpi bilevel CCITT image, so a
  600 dpi render is a 2× upsample; every character was read on the native
  bitmap, where the three marginal superscripts (−10, −3, −5) all read cleanly.

## Validation route

Ranked as the builder brief requires; no figure was digitised and none was
needed.

1. **The one stated number**: `alpha_H`, *derived* by a pymrm radial
   finite-volume solve of the wall with Sieverts equilibrium at both faces (not
   by evaluating Eq. 1's formula) from Itoh's own printed `D`, `C_0` and `l_0`.
   Deviation -0.032 <!-- claim: alpha_H_dev_vs_printed_pct --> %, against a value
   printed to three digits. **Read its resolving power with it, and read what it
   is not.** It tests the transcription of `D`, `C_0`, `l_0` and the 200 µm
   wall, the structure of Eq. 1 and the radial assembly *given the geometry*. It
   does **not** test the geometry: the radii were selected by requiring this
   same identity to close, so the check chose them and cannot catch them being
   wrong. What it does do is discriminate the two candidate readings of the
   printed "17.0 mm" (-2.4 <!-- claim: alpha_H_od_reading_dev_pct --> % for the OD
   reading). It is also nearly blind to `nu`
   (1.16 <!-- claim: wall_nu0_break_thin_pct --> % against the baseline solve on
   a wall 2.4 <!-- claim: wall_over_radius_pct --> % of the tube radius) and to
   the grid (7.2e-7 <!-- claim: wall_grid_rel_err_n5 --> at `n_r = 5`). The same `nu`
   defect costs 45 <!-- claim: wall_nu0_break_thick_pct --> % on a wall with
   `r_o/r_i = 3`, which is how the page shows the assembly is right. And
   `alpha_H` constrains only `r_o/r_i`
   (= 1.023522 <!-- claim: ratio_locus_ro_over_ri -->), so the printed 200 µm
   wall is what turns a ratio into a radius.
2. **Internal identities.** The wall solve against its closed form, second order
   over seven grids; against a quadrature of the resistance integral on a wall
   with position-dependent `D`, where no closed form exists
   (2.9e-7 <!-- claim: graded_D_fv_vs_quadrature_rel -->). The channel solve
   against the exact inverted antiderivative `G` of the half-power law, first
   order, 2.4e-4 <!-- claim: channel_fv_vs_exact_maxabs_n3200 --> at
   `n_z = 3200`; the extinction-length invariant `L + G(w)/Pi = L_ext` to
   1.4e-4 <!-- claim: extinction_length_rel_err --> at every cell; the same
   ladder *above* `Pi_crit`, at Itoh's own `Pi = 51.4 <!-- claim: Pi_itoh -->`
   with the marched initial guess, first order with every value non-negative;
   and against an independent quadrature reference with a loaded permeate.
3. **Every one of them broken on purpose**, in a defect-injection table in which
   every claim fails at least one row. Of the five channel defects, four
   converge to a *wrong answer* and move the error by
   217 <!-- claim: channel_break_min_ratio -->× (`w` substituted for the mole
   fraction), 340 <!-- claim: channel_break_ratio_exponent -->× (exponent 1),
   1028 <!-- claim: channel_break_ratio_nu1 -->× (`nu = 1`) and
   2898 <!-- claim: channel_break_ratio_sink_sign -->× (sink sign); the
   outlet-Dirichlet variant stops Newton with a singular Jacobian.

**Nothing here is compared with a measurement**, because the source contains
none. The values are tier 2 (printed in a paper); the validation situation is
tier 6, and the page says so in its own Data section.

## The new content

Two numbers that are not in the source, and one that is:

- `Pi = alpha_H / u_H0`, the membrane's capacity divided by its duty, and
  `Pi_crit = G(1, beta) = 2.84 <!-- claim: Pi_crit -->`, below which permeation
  is the bottleneck;
- `y_H ~ Pi^(-1/n)`: under a half-power law spare capacity is worth its
  **square** (50× capacity gives a
  2500 <!-- claim: spare_capacity_50x_square -->× lower driving mole fraction);
- Itoh's own tube at `Pi = 51.4` against `Pi_crit = 2.84`: **permeation is
  18 <!-- claim: Pi_over_Pi_crit_itoh -->× from being the limiting resistance.**
  That is the condition (with equilibrium kinetics) that puts `H1.4`'s model on
  its co-current fast-permeation ceiling, and on that ceiling the conversion is
  a function of `K_p` and the purge split alone, so the permeance has no
  observable left to move. That *explains*, rather than restates, `H1.4`'s
  verified finding.
  **Two numbers that go with it, and they describe different problems.** The
  same membrane *run as a separator on a clean permeate* against the duty his
  reaction sets would strip it in
  5.5 % <!-- claim: L_ext_over_tube_itoh_clean_permeate --> of the tube at a
  driving mole fraction of 3.8e-4 <!-- claim: y_H_required_half_power_clean_permeate -->.
  Itoh's actual reactor does neither: his permeate is **loaded** with an argon
  purge, so at the ceiling both sides carry the same hydrogen mole fraction,
  `y = 7.2e-3 <!-- claim: y_H_both_sides_itoh_loaded_permeate -->`,
  19 <!-- claim: y_both_over_y_required -->× that driving fraction; the driving
  difference tends to zero, and with `y' > 0` there is no finite extinction
  length at all (check 5b).

The exponent is then given its consequences: a permeator with `n < 1` empties in
a **finite** length while `n = 1` never does (the logarithm diverges), and a
finite surface step in series with the wall produces an apparent exponent that
runs from ½ to 1 and is not even constant at fixed resistance ratio. That
section is a derivation from the mechanism with one dimensionless group,
compared with no data, and Itoh's membrane cannot be placed on its axis because
no surface kinetics are printed.

## Files

- `index.py`: the edited source of the page (jupytext percent format), paired
  with `index.ipynb`, the executed page
- `model.py`: the wall solve, the `Permeator` channel, the exact channel
  solutions and the surface-step series model; `test_model.py` checks it
- `data/itoh-1987-permeation.csv` + `.meta.yaml`: everything Itoh prints about
  permeation, with the attribution recorded per value
- `agreement.json`: the page's claims (format 2), written by its last cell

## Honesty notes

The reconstruction that the page depends on: `r_i = 8.5 mm`, `r_o = 8.7 mm`.
The running text says "200 µm thick, 17.0 mm OD"; only the *ID* reading
reproduces the printed `alpha_H`. **That comparison is circular as evidence
about the tube**: it picks between two readings of one sentence, and the radii
it selects are its own output, so the page says so and adds an independent
witness. **Figure 1 dimensions the same tube** and prints the values `0.2`,
`17ϕ`, `28ϕ` and `140` `[mm]`. Re-read on the native 300 dpi bitmap of journal
page 1577, the `17ϕ` arrowheads terminate on the **inner** faces of the two
membrane lines, while the wall is dimensioned separately as `0.2` (its tips
straddle the band from outside); the page notebook lists the pixel rows. Read
that way Figure 1 gives the **bore** as 17 mm independently of `alpha_H`, and in
contradiction to the text's "OD". The printed values are facts; the arrowhead
reading is an interpretation and is labelled as one, and a non-blocking
follow-up on `queue_cases/H1.1.yaml` asks the maintainer to confirm it. The
page's `alpha_H` argument stands either way, stated as a discrimination rather
than a validation.

Caravella et al. (2010), named in the catalog entry for the empirical pressure
exponent, is **not on disk and is not used**. Nothing on this page is a claim
about measured exponents in real palladium; a page that tests the exponent
against permeation data needs such a source and is recorded as a follow-up on
`queue_cases/H1.1.yaml`.

The channel model has a trap in the **choice of unknown**, and the page names it
as one rather than as a property of the law. Above `Pi_crit` the extinction
point falls inside the domain, and iterated in `w` from the flat state
`pymrm.newton` stops without converging. The discrete problem is fine: the cell
residual is continuous and strictly increasing, negative at `w = 0` and positive
at `w_{i-1}`, so a unique non-negative root exists at every cell and every `Pi`;
past the front it decays quadratically per cell, underflows to exactly 0 within
9 <!-- claim: tail_cells_to_zero --> cells of falling below 1e-3, and `w = 0`
satisfies the equation thereafter. What fails is the iteration in `w`: the
sink's slope grows like `w^(-1/2)`, every tail cell's Newton step overshoots past
the bound `w = 0`, and `newton`'s bound rule holds the step length at 0.5, so the
error only halves per iteration. The page solves for `s = sqrt(w)` instead, in
which the sink is linear near zero: `newton` then converges from the flat state
at every swept `Pi`, and an independent cell-by-cell march of the same equations
(`Permeator.march_guess`) agrees below 1e-9 <!-- claim: s_form_vs_march_maxabs -->.
The solve in `w` stays on the page, stalling at the
16 <!-- claim: sweep_points_w_form_stalls --> sweep points above `Pi_crit` and
only there, because the difference between "no solution" and "an unknown in
which Newton cannot reach the solution" is what a later reader needs.

## Corrections

- 2026-10-30: stated that the flat-start Newton stall above `Pi_crit` comes from
  `NumJac`'s absolute perturbation floor, with the line search parking on a
  slightly negative `w` where a `max(w, 0)` clip is non-smooth. The cause is the
  unknown: every tail cell's Newton step in `w` crosses `w = 0`, and the bound
  rule holds the step length at 0.5 (`w_form_step_at_half_fraction_min`). Solved
  for `s = sqrt(w)`, the same `newton` converges from the flat state
  (`s_form_flat_converges_everywhere`).
- 2026-10-30: the Reuse section said continuation in the capacity ratio works
  too. In `w`, `pymrm.continuation` stops just above `Pi_crit`
  (`continuation_last_pi_over_crit`) and names a fold the branch does not have;
  in `s` it reaches every value (`continuation_s_reaches_60`).
- 2026-10-30: labelled the channel breaks `w` substituted for the mole fraction
  (217×) and `nu = 1` (1028×) "not converged". Both converge, to wrong answers,
  and the numbers stand (`channel_break_ratio_mole_fraction`,
  `channel_break_ratio_nu1`).
- 2026-10-30: said that `H1.4` records the same flat-start trap; `H1.4` solves
  from its inlet state with `pymrm.newton` and no marched guess, so the
  cross-reference is removed.
