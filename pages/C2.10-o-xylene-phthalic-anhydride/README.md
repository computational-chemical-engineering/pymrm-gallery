# C2.10: o-Xylene to phthalic anhydride, the classic hot spot

Froment's 1967 demonstration that effective kinetics inside a two-dimensional
pseudo-homogeneous model reproduces fixed-bed hot spots, the paper every
reactor textbook cites for it, rebuilt from the printed constants: both
models, every stated number, and the disagreement between their runaway limits
measured continuously.

- **Structures:** `S2` (1-D plug flow), `S6` (2-D radial dispersion, marched)
- **Reference:** Froment (1967), Ind. Eng. Chem. 59(2) 18–27,
  doi:10.1021/ie50686a006
- **Covers:** `D3.4` (the multitubular-reactor half of the same paper; the
  scope argument is recorded in `queue_cases/C2.10.yaml` and
  `queue_cases/D3.4.yaml`)
- **Runtime:** about 3 min
- **Source:** `index.py` (paired with `index.ipynb`); model code in `model.py`
  (`Tube2D`, `Tube1D`, `Tube1DSteady` take a pointwise source; they assume
  plug flow, no axial dispersion, constant dispersion and wall coefficients
  and a coolant constant along the tube; the o-xylene chemistry is
  `froment_source`), tested by `test_model.py`

## Agreement

Five stated model results reproduce to **4.0 <!-- claim: stated_mean_dev_pct --> % mean,
worst 8.6 <!-- claim: stated_worst_dev_pct --> %** (the 363 °C hot-spot rise, a
quantity that doubles every ~2 °C of inlet temperature there). The 2-D hot
spot at 357 °C comes out 29.53 <!-- claim: tm2d_hotspot_357_C --> °C against the
paper's "about 30". The 1-D runaway limit, 363.93 <!-- claim: t0crit_1d_C --> °C,
lands inside the paper's integer-probing bracket (363, 365]; the 2-D limit,
360.006 <!-- claim: t0crit_2d_C --> °C, is at the paper's 360 within 0.01 °C.
Which side of 360 it falls on is set by how the inputs are read, not by the
numerics: with the printed Pe_hR = 5.25 instead of the recomputed 5.2433 it is
359.964 <!-- claim: t0crit_2d_printed_pehr_C --> °C. The continuous gap
between the models is 3.93 <!-- claim: runaway_gap_C --> °C against the paper's
"within five degrees". Four printed transcription identities (Pe_hR, and U
for three λ_R/α_w combinations) close to 0.29 % or better. **Nothing is
fitted anywhere.**

The paper's "a length of 3 meters is sufficient to reach the maximum in
phthalic anhydride concentration" at 357 °C (p. 25, quoted in full in the
stated-results CSV) is not reproduced: the maximum is at
6.6 <!-- claim: z_pa_max_357_m --> m, and at 3 m the conversion is
10 % <!-- claim: pa_shortfall_at_3m_357 --> below it. At 360 °C it almost holds: the
maximum is at 3.2 <!-- claim: z_pa_max_360_m --> m, still past the exit.

## Refresh of 2026-10

- New layout: `index.py`, `model.py`, `test_model.py`, bootstrap cell, pymrm
  pin 3.1.0, every number registered with `Claims` (78 claims, one structural, the others each failed by
  at least one of 95 break rows that re-run its computation).
- `model.py` holds the tube classes, independent of the chemistry: state
  `(n_r, n_f)`, constant dispersion operator assembled once, `NumJac` on the
  pointwise source, and the exact sparse Jacobian handed to the BDF
  integrator (before, `solve_ivp` built it by finite differences; the two
  routes agree to below 1e-5 °C in the 357 °C hot spot, a comparison kept on the
  page).
- The runaway limits are bisected to 1e-4 °C instead of 0.02 °C: 2-D
  360.0078 to 360.0060 <!-- claim: t0crit_2d_C -->, 1-D 363.9297 to
  363.9333 <!-- claim: t0crit_1d_C -->, gap 3.9219 to 3.9274 <!-- claim: runaway_gap_C -->.
  Every other published metric moved by less than 1e-6 relative.
- The 1-D boundary-value route starts Newton from zero instead of from the
  march, and raises if it does not converge (it used to assert an absolute
  residual below 1e-8).
- pymrm's `pseudo_transient` and `continuation` do not apply: the page solves
  no steady problem with more than one solution, and the runaway limit is a
  threshold on a marched outcome, not a fold.
- Correction (2026-10-30), also on the page: this README and the page
  metadata said both limits land inside the paper's brackets, "360.01 in
  (357, 360]", while the notebook said the 2-D limit was outside. Neither is
  right as stated: the 2-D limit is at 360 within the input reading (printed
  or recomputed Pe_hR). The page quoted the 2-D margin as 0.0078 °C, the
  midpoint of a bisection bracket; it is
  0.006 <!-- claim: t0crit_2d_minus_360_C --> °C. It said the Pe_mR = 8 shift
  of the runaway limit, 0.04 °C, was below the 0.05 °C bisection tolerance;
  the two limits had been computed on different grids and tolerances, and on
  equal footing the shift is 0.043 <!-- claim: pe8_t0crit_shift_C --> °C,
  resolved and still small. It said the phthalic anhydride maximum at 357 °C
  sits "at or just beyond the 3 m exit"; it is at
  6.6 <!-- claim: z_pa_max_357_m --> m, with the exit
  10 % <!-- claim: pa_shortfall_at_3m_357 --> below it.

## Provenance: tier 6, reproduction not validation

The paper contains **no experimental data**. Its kinetics are asserted as
"fairly representative" of o-xylene oxidation on V₂O₅, never fitted to
measurements in this document, and every comparison target on this page is
Froment's own 1967 computed result. Calderbank (the catalogue row's second
name) is not on disk. Figures 7–14 are used qualitatively only; nothing was
digitised.

## Three typographic traps, all resolved from the paper's own numbers

1. **"G = 4.684 kg./sq. meter hr." is 4684**, a European thousands separator.
   Proof: p. 19 prints 4684 outright, and the printed Pe_hR = 5.25 equals
   G·c_p·d_p/λ_R only with 4684.
2. **"ΔH₃ = −1.090 kcal./gram mole" is −1090**, the same separator; −1.090
   would make total combustion 282× *less* exothermic than partial oxidation.
3. **"44 gram moles/cu. meter" is 44 g/Nm³**, Figure 11's own axis label; the
   arithmetic to N_A0 = 0.00924 closes for g/Nm³ and is off by 100× for mol/m³.

The break table runs all three taken literally.

## Cross-page

Shared constants are reconciled against the independently transcribed
parameter CSV of `D2.2` (Van Welsenaere & Froment 1970, the same reactor,
same k₁ three years later): eight of nine rows identical to <0.1 %. D2.2's
U = 82.7 is *inherited from this 1967 paper*, so that row is a shared source,
not a second witness; the independent check on U is the printed λ_R/α_w route.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/C2.10-o-xylene-phthalic-anhydride/index.py   # after editing index.py
python scripts/run_kernels.py pages/C2.10-o-xylene-phthalic-anhydride
GALLERY_OFFLINE=1 python -m pytest -q pages/C2.10-o-xylene-phthalic-anhydride/test_model.py
```

No PDF is needed to rebuild: both datasets are committed CSVs.
