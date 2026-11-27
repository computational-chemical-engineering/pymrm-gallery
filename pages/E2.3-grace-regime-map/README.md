# E2.3: Grace's regime map, Geldart's boundaries for any gas, and where his two B/D criteria part

The classification paper by Grace, *Can. J. Chem. Eng.* **64**(3) 353-363 (1986),
[doi:10.1002/cjce.5450640301](https://doi.org/10.1002/cjce.5450640301), writes Geldart's
four powder groups as boundaries in the dimensionless particle size $d_p^*$ as a function
of the density ratio $\Delta\rho/\rho$, for any gas. The page tests whether every
boundary number he prints follows, within its printed rounding, from the correlation and
constants he states for it, and measures where his two B/D criteria (the cloud criterion,
eqs 18-19, and the Ergun-term balance, eqs 20-21) part. A correlation page: nothing is
discretised, and pymrm contributes only `pymrm.checks.find_roots`.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/grace-1986-printed-numbers.csv` | every constant, table cell and boundary number the page uses, as printed (strings) |
| `data/check_values.csv` | the page's functions at one test condition, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it. The page borrows no dataset
from another page.

## Things to know before reusing this

- 36 <!-- claim: boundary_numbers_within --> of the 36
  <!-- claim: boundary_numbers_within_of --> comparisons with a number Grace prints are
  roundings of what his own equations and constants give, taking his printed constants
  as exact (eq 15's four with the page's air viscosity). The closest is Richardson's 53,
  computed 52.507 <!-- claim: t6_richardson_dBD -->, at
  0.985 <!-- claim: boundary_numbers_max_halfunits --> of a half-unit of its last digit;
  the rounding of his $C_1$, $C_2$ carries it across the rounding edge, so that 53 cannot be decided.
- Grace prints no gas properties. Eq 15 (C/A) reproduces only for an air viscosity
  between 18.19 <!-- claim: mu_window_lo_uPas_at_rhoT --> and
  18.41 <!-- claim: mu_window_hi_uPas_at_rhoT --> µPa s at the page's air density.
- Eq 19 (cloud criterion) and eq 21 (Ergun terms equal) are both built on eq 10 with the
  same $C_1$, $C_2$, and cross once at $\Delta\rho/\rho$ =
  437.1 <!-- claim: x_cross_19_21_exact -->, exactly twice eq 19's constant $b$, at the low
  end of what Grace calls "of order 500 to $10^4$". It reaches 500 for $\epsilon_{mf}$ =
  0.449 <!-- claim: eps_mf_for_cross_500 --> instead of the value he assumes, or $d_b$ =
  114 <!-- claim: db_over_dp_for_cross_500 --> $d_p$ instead of 100. Over 500 to $10^4$ they differ by a factor
  1.13 <!-- claim: ratio_21_over_19_500 --> to 4.44 <!-- claim: ratio_21_over_19_1e4 -->
  in $d_p^*$. For glass in air at the test condition, B/D is
  477 <!-- claim: eq19_dp_um_T --> µm by eq 19 and 1190 <!-- claim: dpBD21_um_T --> µm by
  eq 21. Grace recommends eq 21; the page measures the gap and does not choose. State
  which one you use.
- Eq 17 is eq 16 over $\Delta\rho/\rho$ = 1000 to 2000 only; use eq 16 elsewhere.
- Eq 19 has no group D below $\Delta\rho/\rho = b$: eq 18 has no root there.
- Fig 3 is not digitised, and Geldart's measured powders are not placed on the map.
