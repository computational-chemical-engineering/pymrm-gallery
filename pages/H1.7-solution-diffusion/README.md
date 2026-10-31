# H1.7 — Wijmans–Baker solution–diffusion model

From Wijmans, J. G. & Baker, R. W., *The solution-diffusion model: a review*,
*J. Membr. Sci.* **107**(1–2) 1–21 (1995),
[doi:10.1016/0376-7388(95)00102-I](https://doi.org/10.1016/0376-7388(95)00102-I).

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the edited source, paired with `index.ipynb` (jupytext percent format) |
| `index.ipynb` | the page; runs in about 37 s |
| `model.py` | the feed-side film: `film_wall_ratio` (pymrm solve), `film_analytic`, `polarised_fluxes` |
| `test_model.py` | import, one film solve, one stored claim via `load_claims` |
| `data/wijmans-baker-1995-fig5.csv` | 12 markers digitised from Figure 5, reviewed |
| `meta.yaml` | page metadata, validated against `models.yaml` |

## The argument

`A` is fitted to the water flux alone (Eq. 37) and `B` to the salt flux alone
(Eq. 40). Rejection is then not free: the permeate concentration is fixed by
the ratio of the two fluxes, so the third panel of Figure 5 is a **prediction**
from the first two, and that is what makes the figure worth reproducing.

Two results:

**The rejection panel cannot test the model it appears to test.** The prediction
is a rise of 0.31 <!-- claim: rejection_rise_pp --> percentage points across the
measured range: about 4 <!-- claim: rejection_rise_px --> px on a figure whose
printed line is 6 px thick. The markers sit
0.30 <!-- claim: marker_gap_min_pp --> to 0.54 <!-- claim: marker_gap_max_pp -->
percentage points above the prediction, about one line width.

**Concentration polarisation, assumed away in the derivation as transcribed
here, is quantified and excluded as the cause of the intercept gap.** The pymrm
film solve gives a wall enrichment of 29 <!-- claim: film_wall_50um_pct --> % for
a 50 µm boundary layer at the highest measured flux. On the water flux alone, a
98 <!-- claim: pol_delta_star_um --> µm film with a true permeability
1.8 <!-- claim: pol_A_true_ratio --> times the apparent one reproduces both the
measured slope and the fitted intercept (322 <!-- claim: delta_pi_psi --> psi
against the stated ~350). But the salt flux would follow the wall: B computed
from bulk would rise 29 % <!-- claim: pol_B_trend_pred --> across the range,
against a measured 2.4 % <!-- claim: B_trend_frac -->, 35
<!-- claim: pol_B_trend_gap_se --> standard errors. The measured trend allows a
film of at most 13 <!-- claim: salt_delta_max_um --> µm, which moves the
intercept by 0.3 <!-- claim: salt_dpi_shift_psi --> psi. The curvature of the
water flux points the same way only weakly (p = 0.09
<!-- claim: curvature_gap_p --> on one degree of freedom).

## The one trap

The film's wall boundary condition uses the **outward** normal, which at the
downstream face points in +x, so the rejected-salt term enters with a minus
sign. Get it wrong and the wall concentration comes out *below* bulk, salt
depleting against a membrane it cannot cross.

## Open

The fitted osmotic pressure sits 28 <!-- claim: intercept_gap_psi --> psi below
the caption's stated value, 3.8 <!-- claim: intercept_gap_se --> standard errors
of the fit with two degrees of freedom. Polarisation does not explain it unless
the intrinsic B falls by 25 % <!-- claim: pol_B_fall_needed --> across the range
in exact step with the wall; nothing is ruled in. The page says so rather
than explaining it away.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/H1.7-solution-diffusion/index.py   # after editing index.py
python scripts/run_kernels.py pages/H1.7-solution-diffusion
GALLERY_OFFLINE=1 PYTHONPATH=shared python -m pytest -q pages/H1.7-solution-diffusion/test_model.py
```

## Corrections

- 2026-10-30: stated that concentration polarisation is ruled out as the cause
  of the intercept gap, because a film thick enough would halve the apparent
  permeability "which the data does not show"; that test held the true
  permeability at the measured slope, which is not known, and the water flux
  fixes only the product of permeability and film thickness. With the
  permeability matched to the slope, a 98 <!-- claim: pol_delta_star_um --> µm
  film reproduces slope and intercept. What excludes it is the salt flux: B
  computed from bulk would rise 29 % <!-- claim: pol_B_trend_pred --> across the
  range, against a measured 2.4 % <!-- claim: B_trend_frac -->, and the measured
  trend allows at most 13 <!-- claim: salt_delta_max_um --> µm, worth
  0.3 <!-- claim: salt_dpi_shift_psi --> psi.
- 2026-10-30: stated the wall enrichment at the highest measured flux as 25 % and
  76 % (notebook) and 23 % and 70 % (README, `meta.yaml`) for 50 and 200 µm
  films; at the highest measured flux a 50 µm film gives
  29 <!-- claim: film_wall_50um_pct --> % and a 200 µm film
  175 <!-- claim: film_wall_200um_pct --> % (100 µm:
  66 <!-- claim: film_wall_100um_pct --> %). The profile figure labelled "at the
  highest measured flux" used a flux above the highest measured one.
- 2026-10-30: stated that B shows no trend with pressure; the trend is
  2.4 % <!-- claim: B_trend_frac --> across the range (standard error
  0.8 % <!-- claim: B_trend_se_frac -->), against the
  2.22 <!-- claim: B_convective_ratio -->-fold rise a salt leak carried by the
  water would give.
- 2026-10-30: stated that the rejection markers agree with the prediction and
  "would have agreed with any prediction above 99 %"; they sit
  0.30 <!-- claim: marker_gap_min_pp --> to 0.54 <!-- claim: marker_gap_max_pp -->
  percentage points above it, about one line width.
- 2026-10-30: read the film's wall value from the last cell centre and stated
  its error against the closed form as 4.3e-4 at n = 400; on the membrane face it
  is 8.1e-5 <!-- claim: film_err_n400 -->.
- 2026-10-30: stated that a systematic error in the pressure axis would be
  amplified in the intercept; an offset moves it one for one.
