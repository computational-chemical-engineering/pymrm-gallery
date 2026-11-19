# J4.4: Luedeking and Piret, lactic acid fermentation at controlled pH

**Six fitted pairs, one printed fermentation, and what that lets you check.**

Luedeking and Piret print their relation in the Summary itself:

> *"It was found that the instantaneous rate of acid formation dP/dt, could be
> related to the instantaneous rate of bacterial growth dN/dt, and to the
> bacterial density N, throughout a fermentation at a given pH, by the
> expression* dP/dt = α dN/dt + βN *where the constants α and β are determined
> by the pH of the fermentation."*

Luedeking, R. and Piret, E. L., *A Kinetic Study of the Lactic Acid
Fermentation. Batch Process at Controlled pH*, J. Biochem. Microbiol. Technol.
Eng. 1(4), 393-412 (1959), doi:10.1002/jbmte.390010406. The paper is on disk and
every page was read at the file's native resolution.

**The paper uses Greek α and β**; the catalogue writes *a* and *b*. The page keeps
the paper's symbols.

**Three printed tables and nothing digitised.** Table I (book p. 401, bacterial
density), Table II (p. 402: *N*, *P*, d*N*/d*t*, d*P*/d*t*, *k* and
(1/*N*)(d*P*/d*t*)) and Table III (p. 410, the six (pH, α, β) triples). Every
decimal point is a British mid-dot, and this scan's text layer relocates and drops
them (it renders book p. 406's `pH 5.4` as `pH 3.4`), so no number came off the
text layer: every cell was read on a crop enlarged to digit scale.

## The scope decision

Tables I and II tabulate one fermentation, the pH 6.0 run (book p. 401: *"As an
example, the rate information for the fermentation at pH 6·0 is tabulated in
Tables I and II"*). The other five appear only as plotted curves in Figs. 4-11 and
as fitted constants in Table III. A `pdftotext` search for "Table" finds the three
captions and two cross-references and nothing else, and the only ruled blocks are
the three tables. **So five of Table III's six rows are transcribed, used for
trends, and not refitted**, and no curve is digitised to manufacture them.

## What it finds

**The refit is a consistency check on the authors' own fit, not a validation.**
Book p. 410: *"The constants α and β are determined from the plots of Fig. 11"*,
and Fig. 11's axes are Table II's last two columns. Least squares on those
22 <!-- claim: n_rows_eq3 --> points gives α = 2.2666 <!-- claim: alpha_refit_eq3 -->
and β = 0.5287 <!-- claim: beta_refit_eq3 --> against the printed pair,
+3.03 % <!-- claim: alpha_refit_eq3_vs_printed --> in α and
-3.88 % <!-- claim: beta_refit_eq3_vs_printed --> in β, the printed values
0.51 <!-- claim: printed_alpha_in_se --> and
0.43 <!-- claim: printed_beta_in_se --> standard errors from the fit, so both
inside one (`printed_pair_inside_one_se`).

**The printed pair sits on the ridge.** The constants are correlated at
-0.89 <!-- claim: alpha_beta_corr_eq3 -->. Holding α at the printed value returns
β = 0.5513 <!-- claim: beta_given_printed_alpha -->, and the printed pair costs
0.65 % <!-- claim: eq3_printed_pair_rms_penalty --> in rms against the
least-squares optimum: their straight edge and least squares agree.

**Three routes.** Eq. (2), the same model weighted by *N*², and an integral form
that fits the *P* column with no differentiated column at all (α =
2.2727 <!-- claim: alpha_refit_integral -->, β = 0.5101 <!-- claim: beta_refit_integral -->).
The three routes spread by 1.03 % <!-- claim: alpha_route_spread --> in α and
3.64 % <!-- claim: beta_route_spread --> in β.

**Two nulls, both the paper's own rejected alternatives** (book p. 410,
*"Neither statement was found to hold true outside of the period of logarithmic
growth"*). Over all rows, acid production proportional to growth rate alone has
5.55 <!-- claim: eq2_null_growth_only_ratio --> times the two-term rms, and
proportional to density alone 4.92 <!-- claim: eq2_null_density_only_ratio --> times;
on the integrated acid curve 14.2 <!-- claim: integral_null_growth_only_ratio --> and
17.8 <!-- claim: integral_null_density_only_ratio --> times. The direction of those
ratios is algebra (each null is a nested submodel), so both halves of the paper's
sentence are measured: on the logarithmic plateau, the rows above the gap in the
printed *k* column, the one-term forms cost only
1.036 <!-- claim: eq2_logphase_null_growth_only_ratio --> and
1.033 <!-- claim: eq2_logphase_null_density_only_ratio --> times, and outside it
5.01 <!-- claim: eq2_outside_log_null_growth_only_ratio --> and
3.71 <!-- claim: eq2_outside_log_null_density_only_ratio --> times. A ratio near 1
on the plateau is degeneracy: there d*N*/d*t* = *k*<sub>c</sub>*N* and the two
regressors are nearly proportional.

**The integral test.** The printed constants in integrated eq. (2) reproduce the
measured acid curve to 1.70 % <!-- claim: P_curve_rms_rel --> rms and the final
printed acid to 3.12 % <!-- claim: P_final_vs_printed -->. Table II's d*P*/d*t* is
the *P* column differentiated graphically, so this is a coherence check between
two reductions of one run.

**The crossover the paper describes in words.** The two terms are equal at
*k*\* = β/α, crossed at *t* = 10.1 <!-- claim: crossover_time_h --> h, root-found
on the *k* column and again on the separate d*N*/d*t* and *N* columns (the gap
between the two is a bounded claim, `crossover_two_routes_rel`). It is quoted
to one decimal because a linear interpolant moves it by
0.16 % <!-- claim: crossover_interp_band_rel -->. Over the printed window the
growth-associated term makes 46.7 % <!-- claim: growth_assoc_share --> of the acid:
neither term is a correction to the other.

**Across pH, β collapses and α rises.** As the pH falls, β falls strictly over all
six levels by a factor of 5.0 <!-- claim: beta_range_ph6_over_ph45 -->, α rises by
1.61 <!-- claim: alpha_range_ph45_over_ph6 -->, so *k*\* falls
8.07 <!-- claim: kstar_range_ph6_over_ph45 -->-fold. Six fitted pairs are six
outputs of the Fig. 11 test, not the test, so the page claims neither that the
form holds at the other five levels nor a functional α(pH) or β(pH).

**Table II checks itself.** *k* and (1/*N*)(d*P*/d*t*) are quotients of its other
columns by the paper's definition (book p. 400). The *k* identity holds to
2.74 % <!-- claim: k_identity_max_rel --> at worst, and
87.4 % <!-- claim: digit_substitution_detection --> of the single-digit
substitutions available in the (*N*, d*N*/d*t*, *k*) triple break it by more.

## Two printed features, reported and not repaired

**The stray mark in Table III's pH 5·6 β cell.** At low resolution it renders as
`-0-49`. The value is 0·49, settled arithmetic first: read as a minus sign, the
strict fall of β with pH fails (a break row turns `beta_monotone_margin`
negative); the paper's Fig. 10 plots β on a zero-based axis with all six markers
above zero; and the mark is a connected component the size of the genuine
mid-dots in the same column, far smaller than the page's one hyphen (census in the
Table III sidecar). No value read off Fig. 10 enters a CSV, a fit or a claim. What
the mark is cannot be settled from one copy.

**Units.** Table III's footnotes are the only place the paper gives α and β units;
the Nomenclature (book p. 411) lists *"α, β Constants"* with none.

Table III's pH column is printed again in running text on book p. 400, and the
notebook checks the CSV against it (`ph_column_matches_text`).

## What pymrm adds

**Not much to the fit.** Refitting a two-parameter linear model is
`numpy.linalg.lstsq`. What pymrm adds is downstream:

- **The relation as a source term.** Eq. (2) drops unchanged into a batch marcher
  (backward Euler, `newton` + `NumJac((1, 2))`) and a steady plug-flow fermenter
  (`construct_convflux_upwind` + `construct_div`, `nu=0`, a `{"outflow": True}`
  outlet read with `compute_boundary_values`). Both residuals are linear in
  (*N*, *P*); every solve checks `success`.
- **A closed form to check against.** With *k*(*t*) prescribed from the printed
  *k* column (no growth law is imported), eq. (2) integrates exactly to
  *P* = *P*₀ + α(*N* − *N*₀) + β∫*N* d*t*. Both discretisations come back first
  order by `pymrm.checks.observed_orders`, observed orders
  between 0.9 and 1.1 <!-- claim: batch_marcher_order --> for backward Euler and
  between 0.9 and 1.1 <!-- claim: pfr_grid_order --> for donor-cell upwind, each with
  a round-off margin, and the first-order Richardson outlet is
  2.2e-5 <!-- claim: pfr_richardson_vs_closed_form --> from the closed form.
- **Crossings as checked roots** by `pymrm.checks.find_roots` over the whole
  printed window, raising unless there is exactly one.

The plug-flow model illustrates composability; nothing in this paper tests it.

## What the page cannot conclude

- That the relation is correct: one run, and the constants were fitted on it.
- That the form holds across pH: that needs the five other runs' raw data.
- Anything about α(pH) or β(pH) as functions.
- That the plug-flow extension describes a real fermenter.

## Cross-page

The notebook loads J4.1's `printed-growth-laws.csv` only to search it for
`luedeking`, `piret`, `lactic`, `product`, `alpha` and `beta`, and finds no hit
(`j41_search_hits`). That file holds growth laws μ(*S*); Luedeking and Piret
propose none, and eq. (2) contains no substrate concentration, so no row of it is
used and there is nothing to reconcile. J4.1 and J4.2 model growth, this page
models product formation given growth, and the two compose.

## Break table

Every claim fails at least one defect-injection row, each defect entering
through an argument of the study that computes the claim: single mis-transcribed
cells of Tables I, II and III (digit and mid-dot substitutions of the kind this
scan makes, and the pH 5·6 mark read as a sign), a dropped low-precision row, a
changed fit window, fit variable, interpolant or quadrature rule, a coarser grid
or step, orders taken against the finest level, loosened newton tolerances, and a
zero-gradient outlet in place of the outflow condition.

## Reuse

- **Fit α and β in the variables you will report them in**: eqs. (3) and (2) are
  one model and two least-squares problems.
- **Report *k*\* = β/α**: it has physical meaning, and the paper's qualitative
  account of the fermentation is a statement about it.
- **Check a printed rate column against the columns it was divided by.**
- **Look for arithmetic before pixels**: a column reprinted in running text, the
  same quantity on a zero-based axis, a monotone trend.
- **Test both halves of the sentence you are quoting.**
- **Fitted parameters at *n* conditions and raw data at one means one condition
  can be checked.** Check it properly and scope the rest out.
- To reuse the models, replace the prescribed *k*(*t*) by a growth law (J4.1,
  J4.2) coupled through *N*.

## Files

| file | what |
|---|---|
| `index.py`, `index.ipynb` | the page (paired by jupytext) |
| `model.py` | the fits, the closed form, the batch and plug-flow models |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | the claims, written by the last cell |
| `data/luedeking1959-table{1,2,3}-*.csv` | the three printed tables, with sidecars |
