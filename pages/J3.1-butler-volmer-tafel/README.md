# J3.1: Butler–Volmer kinetics and its linear and Tafel limits

The interfacial charge-transfer law in the form a reactor model uses it, and the
overpotential at which each of its two standard approximations goes wrong by a
stated amount, checked against the overpotentials the published `J3.4` cell
actually runs at.

- **Structures:** `S1` (time-marched IVP, for the cell used as an operating
  envelope), `S10` (porous electrode with electroneutrality and a distributed
  interfacial reaction)
- **Reference:** Doyle, Fuller & Newman (1993), *J. Electrochem. Soc.* **140**(6)
  1526–1533, doi:10.1149/1.2221597, Eqs. 6, 8, 16, 17, 18, 28, 29, 30
- **Also:** Marquis, Sulzer, Timms, Please & Chapman (2019), *J. Electrochem.
  Soc.* **166**(15) A3693–A3706, doi:10.1149/2.0341915jes: the same kinetics in
  hyperbolic-sine form (Eqs. 1g–1h), and its printed closed-form inversion
  (Eq. 18)

## The source situation: read this first

The Butler–Volmer equation's historical sources, Tafel (1905), J. A. V. Butler
(1924, 1932), Erdey-Grúz & Volmer (1930), **were not consulted**. They are not
on disk, are not open access, and **no DOI is asserted for them**; they are
recorded as conventional attribution only.

The source that *was* read is Doyle, Fuller & Newman (1993), which prints the law
twice: as the two-exponential kinetics of a flat electrode (Eq. 6 with Eq. 8) and
as a mass-action expression for an insertion reaction (Eq. 17 with Eq. 30). Every
equation used here was read off 600 dpi renders of journal pages 1527, 1528 and
1530, because the scan's text layer mangles them.

## Provenance: tier 6, not experimental

Nothing here is compared with a measurement. Doyle's Figure 2 contains no
experimental points, confirmed by the maintainer on review of the original, and
every check is an algebraic identity or a number the paper prints. Do **not**
describe this page as validated against experiment.

## Why this is not covered by `J3.4`

`J3.4` uses Eq. 17 inside a full P2D cell and never examines the law: it does not
derive the exchange current density, does not ask what happens when the transfer
coefficients fail to sum to one, and does not touch either limit; its `arcsinh`
is the exact inversion of the symmetric law, not an approximation. This page is
about the constitutive law, and the two pages check each other: the law's limits
say what the cell may safely assume, the cell says which part of the law's domain
it lives in. `J3.5`'s `bv_invert` silently *contains* the reduction derived here
without stating the α_a + α_c = 1 condition it needs, and `J3.5`'s `C_r_doyle` is
the same interfacial/superficial ratio on a different *i*₀ scale; the factor
between the two is √(u₀(1−u₀)) = 0.0995 <!-- claim: j35_C_r_scale_ratio -->.

## What the page establishes

| what | result |
|---|---|
| Eq. 17 reduced with Eq. 16 vs the printed Eq. 30 | below 1e-10 <!-- claim: eq17_vs_eq30_max_rel --> relative, over random states and four α pairs with α_a + α_c = 1 |
| Eq. 17's zero-current potential vs Eq. 16 | below 1e-13 <!-- claim: eq17_zero_current_vs_eq16_V --> V, root-found; it fails when α_a + α_c ≠ 1 |
| *i* = 0 at zero surface overpotential | below 1e-12 <!-- claim: zero_current_at_zero_eta_rel --> |
| anodic and cathodic partial currents balance at equilibrium | below 1e-12 <!-- claim: partial_currents_balance_rel_structural -->, **structural**, and reported as such |
| that common value = Eq. 30's *i*₀ | below 1e-12 <!-- claim: i0_from_partials_vs_eq30_rel --> over four α pairs *including asymmetric ones* |
| *R*_ct = *RT*/(*nFi*₀), *n* = α_a+α_c | below 1e-7 <!-- claim: Rct_identity_max_rel --> against a central difference on Eq. 17 |
| arcsinh inversion fed forward through Eq. 17 | below 1e-8 <!-- claim: arcsinh_inversion_max_rel --> over five decades of current |
| closed-form Tafel threshold vs root-finding | below 1e-12 <!-- claim: tafel_threshold_closed_vs_numeric --> |
| Doyle's ν²/δ vs Eq. 30 (κ cancels) | -15.7 % <!-- claim: nu2_over_delta_dev_frac -->, unexplained, *and the same number as `J3.4`'s 18.7 % <!-- claim: j34_kappa_route_spread --> κ spread, not a second test* |
| standard form vs Eq. 17 as printed, in the cell | below 1e-6 <!-- claim: bv_vs_eq17_cell_potential_mV --> mV |
| Tafel-substituted first step vs prediction | 0.46 % <!-- claim: tafel_first_step_vs_prediction_frac --> |
| cell potential spread over an 8× grid refinement | 0.149 <!-- claim: grid_V_spread_mV --> mV, **spatial only** |
| temporal error at the production schedule | 0.273 <!-- claim: dt_V_temporal_error_mV --> mV, observed orders at least 0.8 <!-- claim: dt_observed_order --> and at most 1.2 <!-- claim: dt_observed_order_max --> |
| total salt conservation | below 1e-12 <!-- claim: salt_conservation -->, an identity; a collector leak of a thousandth of *I*/*F* takes it to 0.022 <!-- claim: break_salt_collector_leak --> |

### The validity thresholds at Doyle's 373.15 K (α_a + α_c = 1)

| deviation in current | Tafel | linear (α = 0.5) |
|---|---|---|
| 10 % | 77.1 <!-- claim: eta_tafel_10pct_mV --> mV | |
| 5 % | | 35.9 <!-- claim: eta_linear_5pct_mV --> mV |
| 1 % | 148.4 <!-- claim: eta_tafel_1pct_mV --> mV | 15.8 <!-- claim: eta_linear_1pct_mV --> mV |

The page prints the full table at 298.15 K and 373.15 K. The Tafel column has a
closed form, η = *RT* ln(1 + 1/ε)/((α_a+α_c)*F*), and depends on the transfer
coefficients **only through their sum**. The linear column does depend on the
split, and the linear deviation is not monotonic in η for α_a ≠ α_c: the
threshold is its **first** crossing. At α_a = 0.3 the first crossing is at
8.29 <!-- claim: linear_5pct_alpha0p3_mV --> mV and the last at
203.4 <!-- claim: break_linear_threshold_last_crossing_mV --> mV, so a threshold
taken at a later crossing is wrong in the unsafe direction.

### The published cell, at *I* = 10 A/m²

| | η_s | × the 5 % linear threshold | × the 1 % one |
|---|---|---|---|
| porous cathode | 0.032 <!-- claim: cell_cathode_mean_eta_s_mV --> mV mean, 0.74 <!-- claim: cell_cathode_peak_eta_s_mV --> mV peak¹ | far inside | far inside |
| lithium foil | 24.8 <!-- claim: cell_anode_eta_s_I10_mV --> mV | 0.69 <!-- claim: foil_I10_over_linear_5pct -->, inside | 1.57 <!-- claim: foil_I10_over_linear_1pct -->, outside |
| lithium foil at *I* = 20 | 46.6 <!-- claim: cell_anode_eta_s_I20_mV --> mV | 1.30 <!-- claim: foil_I20_over_linear_5pct -->, outside | 2.95 <!-- claim: foil_I20_over_linear_1pct -->, outside |

¹ The peak sits in the first cathode cell and in the first time step, and is
converged in neither the grid nor the time step: 1.519 <!-- claim: corner_peak_eta_mV -->
mV with both refined eightfold, a margin of
10.4 <!-- claim: corner_peak_margin_factor --> times to the 1 % threshold. The
sustained value past the initial transient is
0.146 <!-- claim: sustained_peak_eta_mV --> mV.

The cathode's interfacial current is a small fraction of its exchange current
density, *i*/*i*₀ = 0.0010 <!-- claim: interfacial_i_over_i0 -->. At *I* = 10
both electrodes are inside the **5 %** linear window; they are in genuinely
different regimes only at the 1 % tolerance, or at *I* = 20. What is
tolerance-independent is the 2.9 <!-- claim: eta_separation_decades --> decades
between them. Substituting the linear law in the cathode moves the cell potential
by below 1e-6 <!-- claim: linear_cathode_dV_mV --> mV. Substituting Tafel moves it
by 519 <!-- claim: tafel_cathode_dV_mV --> mV, of the wrong sign, because Tafel
with *i* ≪ *i*₀ has no solution near η_s = 0.

At *I* = 20 the cell reaches Doyle's cutoff at *u* = 0.3707 <!-- claim: u_cutoff_I20 -->
and crosses 1.9 V at *u* = 0.2638 <!-- claim: u_1p9V_I20 -->; the cathode peak
inside that envelope is 1.679 <!-- claim: cell_cathode_peak_I20_cutoff_mV --> mV.
At *I* = 10 the cutoff is at *u* = 0.8312 <!-- claim: u_cutoff_I10 -->, which
the page prints beside `J3.4`'s committed value for the same cell.

### What the identity checks can and cannot catch

Measured by injecting five mis-transcriptions one at a time, not asserted:

- at Doyle's own α_a = α_c = 0.5, **every** exponent swap is a no-op, so the
  whole discriminating power of the Eq. 17 vs Eq. 30 check comes from the
  sweep's hypothetical asymmetric rows, where a single-equation swap lifts the
  residual to 0.95 <!-- claim: defect_eq17_prefactor_swap_eq17_vs_eq30 --> (Eq. 17
  prefactor) or 0.97 <!-- claim: defect_eq30_tail_swap_eq17_vs_eq30 --> (Eq. 30
  tail);
- swapping the `(c_max−c)^α_c c^α_a` prefactor in **both** equations at once is
  undetectable even there, because it is the same code expression on both sides
  and cancels; only the `(c_T−c_s)^α_a c_s^α_c` half is actually tested;
- the partial-current balance never moves under any defect: it is an exact
  algebraic consequence of α_a + α_c = 1.

## Space and time, measured separately

`Cell.march` sets `dt` and `dt_max` from the current alone, so the grid ladder
refines space at a constant temporal error. The `dt_scale` ladder (it multiplies
`dt` and `dt_max` together) measures the temporal part beside it, and on this
cell the temporal part is the larger: 0.273 <!-- claim: dt_V_temporal_error_mV -->
mV extrapolated against a 0.149 <!-- claim: grid_V_spread_mV --> mV grid spread.
The peak |η_s| sits in the first time step, 0.39 <!-- claim: peak_eta_time_s -->
s into a discharge of 16862 <!-- claim: discharge_time_s --> s, which is why it
moves under both knobs and is converged under neither; the sustained peak moves
by 10.6 % <!-- claim: sustained_peak_box_spread_frac --> over the whole
refinement box.

## What pymrm adds

Nothing to the law: Butler–Volmer is one line of arithmetic, and the identities,
the closed-form Tafel threshold and the validity map need only numpy and a
root-finder (`pymrm.checks.find_roots`, which locates every threshold and every
zero-current potential with a uniqueness check, or keeps the first of several
crossings for the linear threshold). What pymrm supplies is the operating
envelope: `J3.4`'s porous-electrode cell, built from `construct_grad` and
`construct_div` with the boundary fluxes on the outward-normal faces, harmonic
face means at the porosity jump, and one `newton` solve per implicit-Euler step
on a `NumJac((N, 3))` Jacobian, each step checked for success. The time-step
orders come from `pymrm.checks.observed_orders`. The grid and time-step ladders
and the salt identity are measured on this page, each with a break row.

## Corrections

- 2026-08-02. The *I* = 20 A/m² envelope was read at the march endpoint, where
  `v_stop = 1.55` V stops the stepper, instead of at Doyle's cutoff. The page
  said "collapses at *u* = 0.40" and gave the cathode peak there as "the largest
  kinetic overpotential anywhere". Read at the cutoff the collapse is at
  *u* = 0.3707 <!-- claim: u_cutoff_I20 --> with a cathode peak of
  1.679 <!-- claim: cell_cathode_peak_I20_cutoff_mV --> mV, and the largest
  kinetic overpotential in these runs is the lithium foil's at *I* = 20,
  46.6 <!-- claim: cell_anode_eta_s_I20_mV --> mV.
- 2026-08-05. The margin of the cathode peak to the 1 % linear threshold was
  taken from the grid ladder alone, at a fixed time step; the page later said
  "it replaces the 15.7× taken from the grid ladder alone". With the grid and the
  time step refined together the margin is
  10.4 <!-- claim: corner_peak_margin_factor --> times.
- 2026-11-19. The page said "With independent adaptive stepping the linear and exact runs differ by 0.44 mV". Run without the replay, the linear-cathode cell differs from the exact one by below 1e-6 <!-- claim: linear_cathode_dV_noreplay_mV --> mV, the same as with it; the larger deviation comes from switching the law at both electrodes, a row of the break table.

## Data

No dataset of its own, and **no figure was digitised**. Two files are loaded
cross-page from `J3.4`: its reviewed parameter set and its stated-results table.
The page lists `J3.4`'s findings about those rows (the reconstructed
conductivity, the Bruggeman exponent of the salt flux, the open factor ν = 2 in
Eq. 5) and whether each matters here. `J3.4` must stay published for this page
to execute.
