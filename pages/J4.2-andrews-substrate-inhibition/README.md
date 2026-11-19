# J4.2: Andrews' inhibition function

**Two claims in one sentence of a 1968 Summary, tested instead of illustrated.**

Andrews puts Haldane's enzyme-inhibition function into a batch and a chemostat
balance and writes, on book p. 707:

> *"Simulation studies show that the primary result of inhibition by substrate
> in a batch culture is an increase in the lag time whereas in continuous
> culture inhibition by substrate may result in process instability."*

He also assumes, on p. 711, that there is **no lag phase**. So the first claim is
about a lag his kinetics manufacture, and it can be measured against a
no-inhibition control. The second is a statement about how many steady states a
chemostat has and which of them are stable, settled by root-finding and the
Jacobian of the two balances.

**The paper prints no table**: the extracted text layer has the word only inside
*"unstable"*. Every parameter set is typeset inside a figure legend, and reading
a typeset legend is transcription. **Nothing here is digitised**, so the page
does not establish empirical adequacy against Andrews' plotted curves. Two
things read *inside* figures are declared where they are used: Fig. 7's typeset
annotation, which is text and settles a sentence of his prose, and the fact that
both of Fig. 9's `S_1` curves start at the origin, which corroborates an initial
condition he never prints. No axis was calibrated, no point extracted, and no
number on the page comes from either curve.

## What it finds

**The batch claim survives, quantitatively.** Going from no inhibition to
Andrews' `K_i` multiplies the time to reach his target biomass by
5.2564 <!-- claim: batch_time_ratio_vs_no_inhibition -->, and of that increase
94.40 % <!-- claim: lag_share_of_delay_Ki2 --> is lag, against
5.65 % <!-- claim: lag_expo_share_Ki2 --> for a slower exponential phase. The
share grows with inhibition strength, from
87.78 % <!-- claim: lag_share_of_delay_Ki10 --> at `K_i = 10` g/l through
91.25 % <!-- claim: lag_share_of_delay_Ki5 --> at 5 g/l. Under Monod the lag is
identically zero by this definition, because the trajectory starts at its
steepest point: inhibition creates a lag where there was none.

**What that split does not say.** The peak specific growth rate changes by only
-19.43 % <!-- claim: peak_growth_rate_drop_rel -->, but the run reaches its peak
at 99.904 % <!-- claim: peak_time_fraction --> of its duration and is above half
the peak rate for 1.7490 % <!-- claim: above_half_peak_fraction --> of it; the
rate the culture actually realises, `ln(X_f/X_i)/t_f`, changes by
-80.98 % <!-- claim: batch_realised_mean_growth_drop_rel -->. During the lag the
biomass rises 114.67 <!-- claim: lag_biomass_fold_rise -->-fold. And the share
is a partition only near Andrews' own target: the lag does not depend on the
target and the delay does, so at a target of 1 g/l the same decomposition
returns 109.4003 % <!-- claim: lag_share_of_delay_at_Xf1 -->.

**Both of his stated numerical results come back out of a closed form.**
Constant yield collapses eqs. (4)-(5) to one autonomous ODE with a rational
integrand, so partial fractions give the batch time exactly (both sympy
residuals zero):

```
t = (1/mu_hat) [ (S_f - S_i)/K_i + (K_s/a) ln(S_i/S_f)
                 + (1 + a/K_i + K_s/a) ln(X_f/X_i) ],    a = S_i + X_i/Y
```

It gives 36.527476 <!-- claim: batch_time_direct_h --> h against Andrews' "36 hr"
(+1.4652 % <!-- claim: batch_time_direct_vs_printed36_rel -->) and
12.371597 <!-- claim: batch_time_seeded_h --> h against his "12.5 hr"
(-1.0272 % <!-- claim: batch_time_seeded_vs_printed12p5_rel -->), and agrees
with LSODA to below 1e-9 <!-- claim: batch_closedform_vs_lsoda_rel --> while
sharing no code with it. The second depends on reading one sentence of p. 715,
which supports three literal readings (reset, increment, or a 2/4/6 ladder) with
a raw spread of 11.5442 % <!-- claim: seed_readings_raw_spread_rel -->. The
typeset annotation inside Fig. 7, *"`S_i` ADDED IN INCREMENTS OF 2.0 GM/L"*, is
the increment reading, so that is the one reported; the reset reading is
-0.0419 % <!-- claim: seed_reset_vs_increment_rel --> from it, and the ladder
reading is the one the annotation excludes.

**The continuous claim is proved rather than illustrated.** At Andrews' `θ` and
`S_0` there are three steady states:

| | `S_1` (g/l) | `X_1` (g/l) | positive or slow eigenvalue (1/h) | |
|---|---|---|---|---|
| lower | 0.0150567 <!-- claim: steady_S_lower_g_per_l --> | 2.492472 <!-- claim: steady_X_lower_g_per_l --> | fast: -73.019142 <!-- claim: lower_node_fast_eigenvalue_per_h --> | stable node |
| upper | 3.984943 <!-- claim: steady_S_upper_g_per_l --> | 0.507528 <!-- claim: steady_X_upper_g_per_l --> | +0.0561790 <!-- claim: saddle_positive_eigenvalue_per_h --> | saddle, eq. (9)'s root |
| washout | `S_0` | 0 | -0.0481080 <!-- claim: washout_eigenvalue_haldane_per_h --> | stable node |

The saddle's sign is not numerical luck: with `μ(S*) = 1/θ`,
`det J = μ'(S*) X_1/(Yθ)`, so `det J` carries the sign of `μ'` exactly. Washout
is stable, so the reactor is bistable, which is what *"may result in process
instability"* means.

**Inhibition barely moves the steady state and completely changes its basin.**
At the same `mu_hat`, `K_s`, `θ` and `S_0`, dropping the inhibition term moves
the operating substrate by 0.3778 % <!-- claim: monod_vs_haldane_S_shift_rel -->
and the operating biomass by
-0.0011 % <!-- claim: monod_vs_haldane_X_shift_rel -->, and takes the washout
eigenvalue to +0.660702 <!-- claim: washout_eigenvalue_monod_per_h --> per hour,
from stable to unstable.

**There are two washout conditions and eq. (11) is only one of them.** Andrews'
`θ_w = 1/mu_hat_m` = 1.244949 <!-- claim: theta_washout_h --> h is the fold and
contains no `S_0`. The other, washout stable whenever `μ(S_0) < 1/θ`, i.e.
`θ < 1/μ(S_0)` = 3.506000 <!-- claim: theta_transcritical_h --> h, depends on
the feed and he does not print it. Andrews operates inside that bistable window.
Under Monod `μ` is monotone, the fold never happens, and only the transcritical
condition survives, as the `D_c = μ_m S_f/(K_S + S_f)` that J4.1 loads from
Rawlings & Ekerdt.

**All six of his stated continuous-culture outcomes reproduce**, and the page
root-finds the thresholds he does not print, on event-terminated integrations.
The critical inoculum is 0.1842466 <!-- claim: critical_inoculum_g_per_l --> g/l;
over the sampled `S_1(0)` Andrews never prints it varies
5.2120 % <!-- claim: critical_inoculum_S1_init_sensitivity_rel -->, and only
0.4071 % <!-- claim: critical_inoculum_S1_init_sensitivity_narrow_rel --> with
the largest sample left out, which is the flattering window and not the one
reported. The critical step is `S_0` =
17.28760 <!-- claim: critical_step_S0_g_per_l --> g/l, and the critical ramp
duration is 0.8720668 <!-- claim: critical_ramp_duration_h --> h, so his 1 hr
ramp clears it by only 14.67 % <!-- claim: andrews_ramp_margin_rel -->. The
critical inoculum and the critical step are each computed twice, by routes that
share only the right-hand side (a root of the outcome, and the saddle's stable
manifold or a separatrix side test), agreeing to below
1e-9 <!-- claim: critical_inoculum_two_routes_rel --> and below
1e-9 <!-- claim: critical_step_two_routes_rel -->.

## What the page adds that Andrews did not have

**Seed the feed and the instability disappears.** He defines `X_0` in eq. (7) and
never gives it a value. A feed biomass of
0.0213935 <!-- claim: feed_seed_removing_hysteresis_g_per_l --> g/l,
0.8583 % <!-- claim: feed_seed_share_of_operating_biomass --> of the operating
concentration, removes the fold, because washout stops being a steady state.
The reported value is the double root `G = G' = 0`, which has no grid in it;
counting steady states on a scan is grid-limited, and on its finest scan it is
still 1.2e-9 <!-- claim: feed_seed_fold_two_routes_rel --> away.

**Plug flow has none of this.** The same kinetics in a pymrm plug-flow
fermenter give an outlet biomass that increases with residence time for any
positive seed: no fold, no saddle, no washout, because plug flow is the batch
trajectory in space. That is the infinite-stage limit of Andrews' own
*"multistage operation with separate substrate supply to each reactor"*, and it
confirms his prediction in that limit only: a finite cascade is not solved
here, and his other two proposals are recycle problems and are untouched.

## What pymrm adds

Little to the 0-D physics: the batch and chemostat balances are two ODEs, and a
pymrm backward-Euler marcher (`NumJac((1, 2))` + `newton`) agrees with LSODA,
which shows the pymrm route is correct rather than improving on Andrews. The
batch closed form owes nothing to pymrm. What pymrm does add is the spatial
comparison: `construct_convflux_upwind` + `construct_div` + `NumJac` + `newton`,
with an `{"outflow": True}` outlet read through `compute_boundary_values`, solve
the plug-flow limit. Both refinement studies give observed orders from
successive differences of at least 0.95 <!-- claim: pymrm_march_time_step_order -->
and at most 1.05 <!-- claim: pymrm_march_time_step_order --> for the marcher, and
of at least 0.95 <!-- claim: pfr_outlet_grid_order --> and at most
1.05 <!-- claim: pfr_outlet_grid_order --> for the plug-flow grid, and the reported plug-flow outlet at `τ = 6` h is the
Richardson value 0.6076461 <!-- claim: pfr_outlet_X_at_tau6_richardson_g_per_l -->
g/l rather than a single grid's value.

## Seven printed defects, reported and none repaired

- **eq. (9) is printed with a bare `+`, not `±`**, read at digit scale. Its own
  paragraph says the quadratic has two roots. The root it returns is the
  **higher** one, which this same paper calls unstable. The page root-finds both
  and checks eq. (9) against the upper one; it does not correct it.
- **eq. (10) is printed `X_i = Y(S_0 − S_1)`** where `X_1` is meant: `X_i` is
  this paper's symbol for the *initial inoculum* (p. 713, Fig. 9's key). The text
  layer of this scan renders `X_1` as `XI`, so only a render could show it.
- **Figs. 4 and 7 print `Y = 0.5 GM/L`** for a yield the text defines as *"mass
  organisms produced/mass substrate utilized"*. What decides it is that printed
  definition plus the dimensions of eq. (5), not the count of legends. The
  numeral is not in question, only the unit.
- **`"Recherces"`** [sic] for *Recherches* in reference 1; **`"step chances"`**
  [sic] for *step changes* on p. 719; **`"innoculum"`** throughout; **`"a
  inhibition function"`** [sic] in the Summary.

## An approximation that is *not* filed as a defect

Andrews defines `K_s` and `K_i` as the low and high substrate concentrations at
which `μ = mu_hat/2`. Solving eq. (1) for that gives
`S = K_i(1 ∓ sqrt(1 − 4K_s/K_i))/2`, so both definitions are exact only as
`K_s/K_i → 0` and **neither has a solution unless `K_i ≥ 4K_s`**, which is
exactly the condition under which his own eq. (2) gives `mu_hat_m ≥ mu_hat/2`.
At the `K_i` of his chemostat figures the exact roots are
0.0304640 <!-- claim: half_rate_low_root_g_per_l --> and
1.969536 <!-- claim: half_rate_high_root_g_per_l --> g/l,
+1.5468 % <!-- claim: Ks_definition_rel_error --> and
-1.5232 % <!-- claim: Ki_definition_rel_error --> off the values he names them
for, and that is the mildest of the three `K_i` he prints: at his Fig. 4 it is
+3.1947 % <!-- claim: Ks_definition_rel_error_Ki1 --> and
-3.0958 % <!-- claim: Ki_definition_rel_error_Ki1 -->, at his Fig. 1
+6.8502 % <!-- claim: Ks_definition_rel_error_Ki0p5 --> and
-6.4110 % <!-- claim: Ki_definition_rel_error_Ki0p5 -->. The page prints all
three. The nearest thing to a qualification is on the same page, *"even when
`K_s` and `K_i` are well separated..."*, and the page **does not claim he was
unaware**.

## What J4.1 left here, and what goes back

`pages/J4.1-monod/` transcribed **Froment eq. (1.5.2-4)** and **Rawlings &
Ekerdt's `mu_m S/(K_s + S + K_1 S^2)`** into its CSV, flagged both *"OUT OF
SCOPE"* for J4.2, and never evaluated either. This page loads that CSV,
**retyping nothing**, and proves in sympy that both are Andrews' eq. (1), with
Froment's `r_m` = Andrews' `mu_hat` and Rawlings' `K_1` = `1/K_i`.

It also finds a trap in the crossing. J4.1's `froment_rm_definition` row records
Froment defining `r_m` as *"the maximum specific rate of biomass growth"*, true
for his Monod eq. (1.5.2-1) and **false** for his inhibition eq. (1.5.2-4), where
the maximum is `r_m/(1 + 2 sqrt(K_S/K_i))`. At Andrews' constants, carrying that
definition across overstates the maximum by
24.4949 % <!-- claim: rm_as_maximum_overstatement_rel -->.

Going the other way: Andrews' reference list carries the primary citations for
Moser (1958), Teissier (1936, spelled `Teissier` where Rawlings spells the law
`Tessier`; neither repaired) and Contois (1959), which J4.1 could not attribute.
**None is on disk, none was consulted, and nothing is attributed to any of
them.**

## Corrections

- 2026-11-19: the notebook said a non-zero-gradient outflow condition is "which is where the extra order would appear". No extra order appears on this first-order upwind scheme: under `{"outflow": True}` the outlet read converges with observed orders from successive differences of at least 0.95 <!-- claim: pfr_outlet_grid_order --> and at most 1.05 <!-- claim: pfr_outlet_grid_order -->.

## Sources

**Origin, on disk and read in full**: Andrews, J. F., *A Mathematical Model for
the Continuous Culture of Microorganisms Utilizing Inhibitory Substrates*,
**Biotechnol. Bioeng. 10**(6), 707-723 (1968), doi `10.1002/bit.260100602`.
Identity confirmed from its own first page. Every page was read on a render at
the file's native resolution and every numeral re-read on a digit-scale crop;
the landscape pages (Figs. 4, 7, 9) were rotated back before reading. **The text
layer was used only for searching**: it turns `X_1` into `XI`, which is the
failure that would have hidden eq. (10)'s defect.

**Cited by Andrews, not on disk, not consulted**: Haldane (1930), the origin of
the inhibition function itself; Dixon & Webb (1964), his source for fitting it;
Monod (1942); Moser (1958); Teissier (1936); Contois (1959); Koga & Humphrey
(1967); Brennan's PACTOLUS report (1964). Nothing is attributed to any of them
beyond what Andrews prints.

**Borrowed**: `pages/J4.1-monod/data/printed-growth-laws.csv`, six rows.
Froment, Rawlings & Ekerdt and Levenspiel were **not opened for this page**;
everything from them arrives through J4.1's transcription.

## Files

| file | what it is |
|---|---|
| `index.py` / `index.ipynb` | the page (jupytext pair) |
| `model.py` | growth laws, batch closed form, chemostat states and classifier, and the two pymrm models |
| `test_model.py` | fast checks of `model.py` |
| `data/andrews1968-printed-model.csv` | every equation, figure-legend parameter set, in-figure annotation and claim used, verbatim, with flags |
| `agreement.json` | the claims, written by the last code cell |

## Validation summary

- the paper's own identities, all exactly zero in sympy: eq. (1) ≡ eq. (6);
  eqs. (2) and (3) re-derived from eq. (1); the batch-time partial fraction and
  its antiderivative; Andrews ≡ Froment ≡ Rawlings;
- two stated numerical results reproduced from a closed form, labelled
  **reproduction, not validation**: they are Andrews' own PACTOLUS output, read
  by him off his own plot;
- six stated outcomes reproduced and the thresholds he omits root-found;
- four quantities computed twice by independent routes, and two transcription
  checks (eq. (9) against a root-find, the transcritical `θ` by two routes)
  labelled as such;
- both error-carrying axes refined, with orders from successive differences and
  a round-off companion;
- a break table in which every non-structural claim has a row that fails it.

## Reuse

Read `K_s` and `K_i` as the constants of eq. (1), not as the half-rate
concentrations they are named for. Do not carry a Monod textbook's definition of
`r_m` across to the inhibition form. Compute with eq. (6)'s form, not eq. (1)'s:
they are the same function but eq. (1) is `0/0`-shaped at `S = 0`, where every
batch and startup run begins. Check **both** washout conditions, not just the
printed one. And root-find your thresholds: a swept crossing would hide how
narrowly Andrews' ramp clears the critical one.
