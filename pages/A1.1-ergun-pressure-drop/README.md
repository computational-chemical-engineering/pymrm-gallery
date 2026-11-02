# A1.1: Ergun, Kozeny–Carman, Darcy–Forchheimer and the Eisfeld wall correction

One page, four correlations, one dataset: Ergun's own.

From Ergun, S., *Fluid flow through packed columns*, **Chem. Eng. Prog.**
48(2) 89–94 (1952), with the wall correction from Eisfeld, B. and Schnitzlein,
K., **Chem. Eng. Sci.** 56(14) 4321–4329 (2001),
[doi:10.1016/S0009-2509(00)00533-9](https://doi.org/10.1016/S0009-2509(00)00533-9).

This directory covers catalog entries **`A1.1`, `A1.2`, `A1.3` and `A1.4`**, which
is what `docs/catalog-A-foundations.md` asks for: *"One page overlaying Ergun,
Kozeny–Carman, Forchheimer, and the Eisfeld wall correction against a single Δp
dataset, with the wall-effect regime highlighted, is far more useful than four
pages."*

## Status: ready, the figure digitisation was reviewed and approved

The maintainer reviewed the numbered overlays on 2026-08-02: the ringed centres
are on real markers, including in the dense chain where markers merge into the
drawn eq. (14b) line, and nothing was picked up off the line, the graph paper or
the panel labels. A few markers per panel were missed. That is recorded as a
**stated limitation of the dataset**: `limitations.recall` in
`data/ergun-1952-fig7-markers.meta.yaml`, because it is a recall limitation:
it costs precision and distorts counts, and does not move the fitted constants
so long as the misses are scattered, which is what was reported.

## Contents

| file | what it is |
| --- | --- |
| `index.py`, `index.ipynb` | the page, as a jupytext pair; edit `index.py`. Runs in about 4 s |
| `data/ergun-1952-fig7-markers.csv` | **244** markers from Ergun's Figure 7, three panels |
| `data/ergun-1952-fig5-present-markers.csv` | 74 markers from Figure 5's top panel: the *same* runs, other ordinate, digitised independently as a cross-check |
| `data/ergun-1952-parameters.csv` | his constants, and the measured level of the line he labels "Kozeny–Carman" |
| `data/eisfeld-2001-wall-correction.csv` | Eisfeld & Schnitzlein Tables 2 and 3 |
| `data/eisfeld-2001-table1-sources.csv` | their Table 1 rows for the sources Ergun plotted; this is where D/dₚ comes from |
| `meta.yaml` | page metadata |

Every CSV has a `.meta.yaml` provenance sidecar. Read
`ergun-1952-fig7-markers.meta.yaml` before reusing the marker file.

## Three things to know before reusing this

**The paper has no text layer at all.** `pdftotext` returns four bytes for four
pages, and there are no tables anywhere in it. Every number on this page was read
off a 600 dpi render, and the dataset was digitised from figures. The catalog's
note that "Ergun's original data are tabulated in the 1952 paper" is wrong: it
is not.

**The round trip is the point, and it is not circular.** Ergun's constants were
least-squares fitted to these very points, so recovering them is the check the
paper pays for, but only if the extraction cannot be reading his drawn line
instead of his markers. It cannot: detection keys on the enclosed *white
interior* of an open circle, a line has none, and the half of the top panel that
contains the printed line and no markers returns zero points. That test is
repeated on the page from the CSV itself.

**"Least squares" needs a weight here.** f_v spans two and a half decades, so an
unweighted fit in Ergun's linear coordinates is dominated by a handful of
high-Reynolds points and returns k₁ = 168.5 <!-- claim: k1_unweighted --> (+12 % <!-- claim: k1_unweighted_relative_error -->). The page quotes the fit with
equal *relative* weight, k₁ = 151.9 <!-- claim: k1_refit_all_points -->, and prints both so the round trip is not
resting on a hidden choice.

## Headline results

Throughout, x = N_Re / (1 − ε), the Reynolds number on Ergun's axis, with
N_Re = Re_dp = G dₚ / μ.

| | |
| --- | --- |
| refit of k₁, k₂ from Ergun's own figure | **151.9 <!-- claim: k1_refit_all_points --> and 1.697 <!-- claim: k2_refit_all_points -->** vs his printed 150 and 1.75 |
| per-panel refit (3 sources, 3 calibrations) | k₁ = 146.6 <!-- claim: k1_panel_present_investigation --> / 154.6 <!-- claim: k1_panel_burke_plummer --> / 151.0 <!-- claim: k1_panel_morcom --> |
| Ergun's eq. vs the 244 <!-- claim: n_markers --> points | bias +0.35 % <!-- claim: ergun_bias_on_own_data -->, mean abs 5.16 % <!-- claim: ergun_mean_abs_dev_on_own_data -->, rms 7.13 % <!-- claim: ergun_rms_on_own_data --> |
| scatter about Ergun's line by panel | Morcom 10.2 % <!-- claim: ergun_rms_morcom -->, Burke and Plummer 5.1 % <!-- claim: ergun_rms_burke_plummer -->, his own runs 3.9 % <!-- claim: ergun_rms_present_investigation --> |
| Fig. 5 ⇄ Fig. 7 cross-check | 1.029 <!-- claim: fig5_fig7_cross_check_ratio --> on binned median f_v, no trend beyond its standard error |
| the line he calls "Kozeny–Carman" | f_v = 149.2 ± 0.5, i.e. 0.995 <!-- claim: kozeny_carman_line_over_k1 --> of his 150, not Carman's 180 (0.829 <!-- claim: kozeny_carman_line_over_carman -->) |
| Kozeny–Carman vs Ergun | +7.5 % <!-- claim: kc_vs_ergun_x10 --> at x = 10, -44.6 % <!-- claim: kc_vs_ergun_x100 --> at 100, -90.5 % <!-- claim: kc_vs_ergun_x1000 --> at 1000; +20 % <!-- claim: kc_vs_ergun_creeping --> as x → 0 |
| Eisfeld (infinite bed) vs Ergun | +2.7 % <!-- claim: eis_viscous_vs_ergun --> viscous, **-24.5 % <!-- claim: eis_inertial_vs_ergun --> inertial** |
| Eisfeld at D/dₚ = 5, vs Ergun | +46.0 % <!-- claim: eis_r5_vs_ergun_x10 --> on f_v at x = 10, -11.2 % <!-- claim: eis_r5_vs_ergun_x1000 --> at x = 1000 |
| the wall correction itself, D/dₚ = 1.624 | +178 % <!-- claim: wall_corr_smallest_x3 --> at x = 3, -24.5 % <!-- claim: wall_corr_smallest_x30000 --> at x = 30000, it changes sign |
| compressible bed, 2 m, G = 1 kg/(m² s) | Δp 6.29 kPa <!-- claim: bed_dp_kpa_ergun --> (Ergun); the wall at D/dₚ = 5 adds 22.4 % <!-- claim: wall_dp_effect_G1 --> to Eisfeld's infinite-bed value |

## Corrections

- 2026-11-02: said that at D/dₚ = 5 the wall "lowers the inertial coefficient by
  8 %" and that the wall's two effects nearly cancel in the compressible-bed
  example. The inertial coefficient rises by 10.3 % <!-- claim: wall_r5_inertial_coeff_change -->
  and the viscous one by 49.4 % <!-- claim: wall_r5_viscous_coeff_change -->, so the pressure
  drop at G = 1 kg/(m² s) is 22.4 % <!-- claim: wall_dp_effect_G1 --> above Eisfeld's
  infinite-bed value, not unchanged. It lands on Ergun's (6.30 <!-- claim: bed_dp_kpa_eis_r5 --> against 6.29 <!-- claim: bed_dp_kpa_ergun --> kPa)
  only because Eisfeld's infinite-bed form differs from Ergun's by -18.2 % <!-- claim: bed_dp_eis_inf_vs_ergun --> at this flow.
- 2026-11-02: said the wall correction "exceeds the gap between any two of them
  once D/dₚ < 10". At D/dₚ = 5 it exceeds the Ergun to Eisfeld gap only up to
  x of about 3.8e2 <!-- claim: wall_beats_eis_gap_x_max --> and the gap to Kozeny–Carman only up to
  x = 71 <!-- claim: wall_beats_kc_gap_x_max -->; at x = 1000 the wall correction is 14.3 % <!-- claim: wall_corr_r5_x1000 --> against an Ergun to Eisfeld gap of
  -22.4 % <!-- claim: eis_inf_vs_ergun_x1000 --> and an Ergun to Kozeny–Carman gap of -90.5 % <!-- claim: kc_vs_ergun_x1000 --> (both below Ergun).
- 2026-11-02: said Kozeny–Carman is "within a few percent of Ergun" for
  N_Re/(1−ε) up to about 3. Its constant is Carman's 180, not Ergun's 150: it is
  20 % <!-- claim: kc_vs_ergun_creeping --> above Ergun as x → 0 and 15.9 % <!-- claim: kc_vs_ergun_x3 --> at x = 3, and within 10 % only for
  x between 7.8 <!-- claim: kc_within_10pct_lo --> and 28.6 <!-- claim: kc_within_10pct_hi -->.
- 2026-11-02: said that on Fig. 6 the Carman line would sit 52 px above Ergun's
  150. At 643 px per decade the gap is 51 <!-- claim: carman_vs_k1_px --> px.
- 2026-11-02: called |f_k(∞) − k₂| "exact by construction". The
  limit is approached as k₁/x, so the gap at x = 10⁹ is not zero; the page now
  bounds it as a theoretical limit gap, below 1e-6 <!-- claim: fk_limit_gap_to_k2 -->,
  and checks separately that it equals k₁/x.
- 2026-11-02: said a 5 % error in ε moves the viscous coefficient by about
  24 %. That held only for ε 5 % low, +24.5 % <!-- claim: viscous_factor_eps_5pct_low -->; for ε 5 % high the
  change is -19.3 % <!-- claim: viscous_factor_eps_5pct_high -->.

## Regenerating

```bash
jupytext --sync index.py                   # after editing index.py
python - <<'EOF'
import nbformat; from nbclient import NotebookClient
nb=nbformat.read("index.ipynb",as_version=4)
NotebookClient(nb,timeout=1800,kernel_name="python3",resources={"metadata":{"path":"."}}).execute()
nbformat.write(nb,"index.ipynb"); print("OK")
EOF
```

The extraction itself is not in this directory: it needs the source scan, which
must never enter the repository. The method is described in full in the dataset
sidecars, in enough detail to reproduce from the PDF.
