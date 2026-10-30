# H1.9 — Maxwell–Stefan mixture diffusion in a zeolite membrane

Mixture diffusion in an MFI (silicalite) membrane, from Krishna, R. & Baur, R.,
*Modelling issues in zeolite based separation processes*, **Separation and
Purification Technology 33**(3) 213–254 (2003),
[doi:10.1016/S1383-5866(03)00008-X](https://doi.org/10.1016/S1383-5866(03)00008-X).

A 95–5 methane/n-butane mixture is fed to the membrane. The pure-component
fluxes predict a butane-over-methane permeation selectivity of
1.73 <!-- claim: a47_sp_pure -->; the measured and computed values are two to
three orders of magnitude larger. This page implements the two reasons the
review gives, a mixture isotherm that needs IAST rather than the
multicomponent Langmuir and **finite** exchange coefficients between the two
sorbates, and scores both against the review's numbers and against the
measurements it quotes.

## What is here

| | |
|---|---|
| Sections of the review | 3.1–3.4 (mixture theory) and 4.2–4.3 (MFI membrane separations) |
| Equations implemented | (38)–(47) and (53)–(65) |
| Data tier | **2**: four *measured* selectivities, printed as numbers in the review's running text |
| Validation | Funke et al.'s measured selectivity reproduced to 4.3 % <!-- claim: funke_mixture_rel_dev -->, nothing fitted, inside a 1.57 <!-- claim: funke_choice_spread --> band of readings the review leaves open, all swept on the page |
| Reproduction | nine printed steady-state targets, worst 0.74 % <!-- claim: worst_steady_target_rel_dev -->; the transient peak **time** is 9.8 % <!-- claim: peak_time_shortfall --> early and is reported as an open tension |
| Files | `index.py` (edited source, paired with `index.ipynb`), `model.py` (IAST, the friction matrix, the membrane), `test_model.py` |

The single-component half of the same review is
[`A4.7`](../A4.7-zeolite-micropore-maxwell-stefan/) and is not repeated here.
The cut between the two cases is single-component / mixture, not micropore /
membrane; the reasoning is on both pages under *Reuse*.

## Data

| file | tier | what it is |
|---|---|---|
| `bakker-funke-mfi-membrane-selectivities.csv` | **2** | Four **measured** permeation selectivities, Bakker's two for CH₄/n-C₄H₁₀ and Funke et al.'s two for n-C₆/3MP, stated as numbers in the review's running text. Both origins are recorded under `origin_not_consulted`: neither is on disk and neither was opened. |
| `krishna-baur-2003-mixture-printed-results.csv` | 6 | Every mixture result the review prints as a number. The authors' own computed output; comparing against it is reproduction. |
| `krishna-baur-2003-mixture-parameters.csv` | 6 | The parameters that appear only in the mixture sections (3.3, 3.4, 4.3). |

Tables 1 and 2, the MFI density and the Section 4.1 single-component results are
loaded from `A4.7`. Every finding `A4.7` states about the rows used is listed in
*The data*, and every borrowed row that is compared with something is printed
beside this page's own value.

## Headline results

**Against measurement** (the only validation on the page):

| | model against measured |
|---|---|
| Funke, nC₆/3MP at 15 kPa | +4.3 % <!-- claim: funke_mixture_rel_dev --> |
| Bakker, 95–5 CH₄/nC₄ | 1.28 <!-- claim: bakker_955_model_over_measured --> times |
| Bakker, 50–50 CH₄/nC₄ | 3.4 <!-- claim: bakker_5050_model_over_measured --> times |

With the exchange coefficients removed the same model is
31 <!-- claim: bakker_955_infinite_too_low --> and
11 <!-- claim: bakker_5050_infinite_too_low --> times *too low* for Bakker's two
feeds. Both measurements lie between the two models, nearer the finite-exchange
one on a log axis, so the data favour finite exchange while agreeing
quantitatively with neither.

**The transient peak time and the membrane thickness.** The review never prints
the thickness; `A4.7` reconstructed it from a single-component flux. The four
mixture fluxes do not test it, because it cancels out of them identically. The
transient methane peak time carries it as δ², and refined in space and time it
lands 9.8 % <!-- claim: peak_time_shortfall --> before the printed value,
needing δ larger by 5.3 % <!-- claim: delta_shift_needed_for_peak_time -->.

**Two printed numbers that do not reproduce**, both investigated: the
Henry-limit sorption selectivity (2649.9 <!-- claim: S_henry --> from Table 1,
against the printed 2200 <!-- claim: printed_sorp_S_henry -->), and the multicomponent-Langmuir counterfactual
(1682 <!-- claim: langmuir_SP_mean --> against the printed
800 <!-- claim: printed_mix_SP_langmuir -->), whose
structural claim does hold.

## Numerics

The steady membrane is `pymrm.newton` with `damping="monotone"` on the
logarithm of the loadings, from a profile linear in the loadings; every
production case converges directly, and `pymrm.continuation` in the feed
pressure is the comparison route. The transient is backward Euler with a
frozen-coefficient inner iteration; `pymrm.pseudo_transient` does not apply,
because the page needs the physical transient. The loadings are one
`(n_z, n_c)` array with the species last, and the steady Jacobian is
`NumJac` with a tridiagonal stencil in space.

## Reproduce

```bash
cd pages/H1.9-zeolite-membrane-maxwell-stefan-mixture
jupyter nbconvert --execute --inplace index.ipynb      # a few minutes
GALLERY_OFFLINE=1 PYTHONPATH=../../shared python -m pytest -q test_model.py
```

Nothing on the page is stochastic; the one random draw (200 states for the
friction-system check) is seeded.

## Corrections

- 2026-10-30: said the finite- and infinite-exchange models are two orders of
  magnitude apart and that Bakker's measurements sit on one side of the gap.
  The models are a factor of about 40 <!-- claim: exchange_factor_955 --> apart,
  and both measurements lie between them
  (`bakker_955_log_position`, `bakker_5050_log_position`).
- 2026-10-30: said the multicomponent Langmuir makes S_P composition-independent
  to a fraction of a per cent; the spread over the five feeds is a factor of
  1.024 <!-- claim: langmuir_composition_spread -->.
- 2026-10-30: said the nine steady-state targets are reproduced "inside the
  printed precision throughout". Measured in half-units of the last printed
  digit, only some lie within the printed rounding
  (`steady_targets_within_printed_rounding`).
- 2026-10-30: the steady grid order was fitted to errors against the finest
  solution as if it were exact, and came out above the scheme's order. From
  successive differences on a further refined ladder it is second order
  (`steady_grid_order`), and the production grid error is measured against the
  Richardson limit (`grid_rel_err_at_nz40`).
