# F1.4: Krishna–Ellenberger two-bubble-class gas holdup

Large-bubble gas holdup in a churn-turbulent bubble column, from Krishna &
Ellenberger, *AIChE Journal* **42**(9) 2627–2634 (1996),
[doi:10.1002/aic.690420923](https://doi.org/10.1002/aic.690420923).

The page tests Eq. 19, a correlation containing no fluid property at all,
against the markers of the paper's Figure 11, and against the Wilkinson
correlation it replaced.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/krishna-ellenberger-1996-fig11.csv` | the marker positions digitised from Figure 11 |
| `data/krishna-ellenberger-1996-parameters.csv` | Table 2 properties, gas densities, Reilly's *B*, the Table 3 deviations |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it. Read
`krishna-ellenberger-1996-fig11.meta.yaml` before using the `gas` column: it
records a maintainer review of the extraction and the reason most rows are
deliberately left `unassigned`.

## Two things to know before reusing this

**This is a correlation, not a discretised model.** Eqs. 8, 19 and 20 are
algebra; no pymrm operator appears on the page, and the page says so. The three
functions `eq19`, `wilkinson` and `reilly_transition` take SI arguments, the
dictionary of published constants and the parameters, and are standalone: lift
them into a reactor model. `F2.3` is the page that does.

**The gas labels are incomplete on purpose.** The four series on Figure 11
differ only by marker shape, drawn as overlapping open outlines. Automatic shape
classification failed in the dense band and a review confirmed it, so only the
SF₆ group is labelled. This costs nothing for the main test, since Eq. 19 has no
gas-density term and every marker can be used without labels, and the SF₆ group
alone carries the independence test.

## Headline result

Eq. 19 sits through the 63 <!-- claim: n_points --> markers with a mean absolute deviation of
13.8 % <!-- claim: eq19_mean_abs_dev --> and a bias of +2.8 % <!-- claim: eq19_bias -->. For the
33 <!-- claim: n_in_window --> markers that are certainly inside the window it was fitted over the
deviation is 11.1 % <!-- claim: eq19_mean_abs_dev_in_window --> with a bias of
-2.1 % <!-- claim: eq19_bias_in_window -->; the other 30 <!-- claim: n_below_window --> deviate by
16.9 % <!-- claim: eq19_mean_abs_dev_below_window --> with a bias of +8.3 % <!-- claim: eq19_bias_below_window -->,
and include every SF₆ marker. Wilkinson, evaluated at the plotted
abscissa and air density, is biased high by 64 % <!-- claim: wilkinson_bias -->, and deviates
4.6 <!-- claim: wilkinson_over_eq19_mean_abs_dev --> times as far as Eq. 19.

For gas independence, comparing group biases does not work on this figure: the SF₆
markers and the rest occupy disjoint velocity windows. The page extrapolates a free
power law fitted on the unlabelled markers onto the SF₆ markers, which sit
11.7 % <!-- claim: sf6_excess_over_extrapolation --> above it as a mean of point ratios, and
10.0 % <!-- claim: sf6_excess_geometric --> as a geometric mean; the second is the one the
1.2 <!-- claim: sf6_excess_z_vs_zero --> standard errors from zero refer to. Wilkinson's
fixed-velocity ratio is 5.9 % <!-- claim: sf6_wilkinson_excess_over_argon --> to
17.5 % <!-- claim: sf6_wilkinson_excess_over_helium --> there, depending on which light gas the
unlabelled markers are; a mix that changes with velocity, or a curve that bends below the
fitted range, moves that expectation by more than the standard error. The figure cannot
distinguish the two.

## Corrections

- 2026-11-03: said Wilkinson's gas-density dependence is "about six times larger than the
  measurements support" (`models.yaml`), and that it "is 2.6-3.5x larger than even that upper
  bound" (page). The comparison set Wilkinson's helium-to-SF₆ spread (old wording: "Wilkinson requires a
  helium-to-SF6 spread of 20 to 27 %") against the extrapolation as if the unlabelled markers were all helium. The SF₆ markers sit
  11.7 % <!-- claim: sf6_excess_over_extrapolation --> above the extrapolation, with a standard error of
  0.080 <!-- claim: sf6_log_excess_se --> in logs, and Wilkinson implies 5.9 % <!-- claim: sf6_wilkinson_excess_over_argon -->
  to 17.5 % <!-- claim: sf6_wilkinson_excess_over_helium --> if the unlabelled markers are one light
  gas throughout and the light-gas power law continues below the fitted range: the markers are inside that range, so the figure does not show Wilkinson's dependence to be exaggerated.
- 2026-11-03: said the SF₆ extrapolation bias of 7.7 % was "an upper bound on any real gas
  effect". It is a point estimate with no such property: stated as model over measured it is
  7.7 % <!-- claim: sf6_extrapolation_bias --> in absolute value, and as measured over extrapolated it is
  11.7 % <!-- claim: sf6_excess_over_extrapolation --> with a standard error of 0.080 <!-- claim: sf6_log_excess_se -->
  in logs.
- 2026-11-03: said the computed Reilly values "bracket the measured values" of Figure 6 (the
  dataset sidecar said it too). They do not at every gas: the argon dense-phase voidage read from
  Figure 6(a) is above the computed 0.1736 <!-- claim: reilly_eps_trans_argon -->. `F1.3` holds the reads.
- 2026-11-03: said "The 14 % scatter is a mixture of real experimental spread and digitisation error,
  and this page cannot separate the two." The sidecar's digitisation error moves the mean absolute
  deviation by only 0.11 % <!-- claim: digitisation_shift_mean_abs_dev --> (a shift in the percentage
  itself), so the scatter is the experiment's.
- 2026-11-03: gave as the reason the Wilkinson-to-Eq. 19 ratio is 4.6 <!-- claim: wilkinson_over_eq19_mean_abs_dev -->
  here and 7.8 <!-- claim: table3_eps_b_wilkinson_over_replacement --> in Table 3 that the paper's points
  "span six column diameters up to 0.63 m" (old wording: "its 1,735 points span six column
  diameters up to 0.63 m"). That is a conjecture, not a result of this page;
  `F1.3` works out the bias it would imply.
