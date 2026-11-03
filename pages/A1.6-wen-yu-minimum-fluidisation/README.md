# A1.6: Wen and Yu's minimum fluidisation velocity

One equation, two constants, and a promise: no voidage, no shape factor. This
page recovers the constants from the force balance they came from, computes the
particle they secretly describe, splits the resulting error between the
approximation and the friction law underneath it, and puts the correlation
against 21 <!-- claim: measured_rows --> minimum fluidisation velocities measured seven years later in another
laboratory.

- **Structures:** `S3`, algebraic closure. **The only pymrm call is `find_roots`**, and the
  page says so where a reader would look for a solver, as `A1.1` and `F1.4` do.
- **Runtime:** ~3 s
- **Data tier:** 2, a table printed in another paper, whose values are that
  author's own laboratory measurements.

## Source

**Wen, C. Y. and Yu, Y. H.**, *A Generalized Method for Predicting the Minimum
Fluidization Velocity*, A.I.Ch.E. Journal **12**(3) 610–612 (May 1966),
[doi:10.1002/aic.690120343](https://doi.org/10.1002/aic.690120343), read
directly, on 600 dpi renders of journal pages 610–612.

**The catalogue citation for `A1.6` was wrong, and usefully so.** It read
"Chem. Eng. Prog. Symp. (1966)". That is a real paper, *Wen, C. Y. and Yu, Y. H.,
Chem. Eng. Progr. Symposium Series No. 62, **62** (1966)*, which is this
communication's own reference 23 and contains the derivation, but it is **not on
disk and was not consulted**. Both are recorded in `meta.yaml`; they must not be
conflated.

**Two other origins cited but not consulted.** The communication attributes its
friction law to *Ergun, S. and Orning, A. A., Ind. Eng. Chem.* **41**, 1179
(1949) and **prints neither that equation nor the derivation**. The Ergun
constants used here are the 1952 pair, read on
[`A1.1`](../A1.1-ergun-pressure-drop/). Davies and Richardson (1966), used for
comparison, comes through Geldart, as on `A1.7`.

**A trap in this PDF.** Its first page of *extracted text* is a different
article, a communication on activity coefficients occupying the top of journal
page 610. Wen and Yu begin partway down, at "C. Y. WEN and Y. H. YU". Confirm the
page before transcribing anything.

## What the paper prints

```
(N_Re)_mf = sqrt( (33.7)^2 + 0.0408 N_Ga ) - 33.7                (1)

(1 - eps_mf) / (phi_s^2 eps_mf^3)  ~=  11                        (2)
1 / (phi_s eps_mf^3)               ~=  14                        (3)

N_Ga = d_p^3 rho_f (rho_s - rho_f) g / mu^2      ("Galileo number")
N_Re = d_p rho_f V / mu
```

covering, per journal page 611, `d_p` 0.002–1.97 in., `eps_mf` 0.385–0.935,
`phi_s` 0.136–1.0, `d_p/D` 0.000807–0.25 and `(N_Re)_mf` 0.001–4000. Table 1
compares equation (1) with Narsimhan's correlations by particle class; the 284
data points behind it appear only as scatter in Figure 4 and were **not
digitised**.

## What the page checks

| route | check | result |
|---|---|---|
| 2, identity | eqs. (2) and (3) put into the Ergun force balance | `24.5 Re² + 1650 Re = N_Ga`, root gives **C₁ = 33.6735 <!-- claim: C1_recovered --> (-0.079 <!-- claim: C1_deviation_pct --> %)** and **C₂ = 0.040816 <!-- claim: C2_recovered --> (+0.040 <!-- claim: C2_deviation_pct --> %)** |
| 2, identity | eqs. (2) and (3) solved *simultaneously* (`find_roots`) | one admissible root: **φ_s = 0.669 <!-- claim: implied_phi_s -->, ε_mf = 0.474 <!-- claim: implied_eps_mf -->**, not a sphere. Forced to φ_s = 1 they demand 0.383 <!-- claim: eps_mf_forced_by_eq2_at_sphere --> and 0.415 <!-- claim: eps_mf_forced_by_eq3_at_sphere --> |
| 2, identity | Table 1's rows pooled against its own overall | Narsimhan **46.24 <!-- claim: table1_narsimhan_pooled --> against a printed 46 <!-- claim: table1_narsimhan_printed -->**; Equation (1) 36.47 <!-- claim: table1_eq1_pooled --> against a printed 34 <!-- claim: table1_eq1_printed -->, which **cannot close on the printed row counts** (negative variance, -1572 <!-- claim: table1_eq1_extra_points_variance -->) and is reported, not repaired. The stronger claim, that the column is unreachable, is *not* made |
| 2, cross-page | A1.1's refit of Ergun's constants pushed through | the prediction moves **-1.27 <!-- claim: ergun_constant_share_viscous_pct --> % (viscous) to +1.56 <!-- claim: ergun_constant_share_turbulent_pct --> % (turbulent)**. Like for like at ε_mf = 0.42: **22.7 <!-- claim: voidage_over_ergun_ratio_viscous -->x at the viscous limit, 1.17 <!-- claim: voidage_over_ergun_ratio_turbulent -->x at the turbulent limit** |
| 2, cross-page | which end of the range is worse depends on the voidage | at 0.42 eq. (1) is -28.9 <!-- claim: dev_vs_exact_ergun_viscous_eps042_pct --> % (viscous) and -1.8 <!-- claim: dev_vs_exact_ergun_turbulent_eps042_pct --> % (turbulent) against the exact balance; the viscous end is the worse one at 3 <!-- claim: voidages_with_viscous_end_worse --> of the 5 voidages the paper quotes, the turbulent end at the lowest two |
| 3, measurement | eq. (1) vs Geldart's 21 <!-- claim: measured_rows --> measured `U_0` | bias **-17.4 <!-- claim: wenyu_bias_all_pct --> %**, MAD 33.3 <!-- claim: wenyu_mad_all_pct --> %, rms **38.1 <!-- claim: wenyu_rms_all_pct --> %** against their claimed 34 <!-- claim: table1_eq1_printed --> % (the paper defines neither statistic; as a standard deviation about the mean it is **33.9 <!-- claim: wenyu_sd_about_mean_all_pct --> %**). Both printed |
| 3, measurement | the spherical Diakon cuts, taking `d_sv = d_p` | **low on every row, by -25.2 <!-- claim: wenyu_bias_diakon_pct --> % on average (-54.5 <!-- claim: wenyu_dev_diakon_min_pct --> % to -5.5 <!-- claim: wenyu_dev_diakon_max_pct --> %)**; Davies and Richardson -0.75 <!-- claim: davies_richardson_bias_diakon_pct --> % with MAD 15.8 <!-- claim: davies_richardson_mad_diakon_pct --> % |
| 3, inversion | the voidage those Diakon cuts demand of the exact balance | median **0.409 <!-- claim: inverted_eps_mf_median_diakon -->**, a reparameterisation of the bias, not an independent closure |
| 3, contrary evidence | the voidage Geldart himself **reports** | ε_MB = **0.444 <!-- claim: geldart_reported_eps_MB_no_inference_value --> with no inference** on the two rows printing `U_MB = U_0`, `H_MB/H_0 = 1.000`; at the reported voidage the exact balance is **+58.4 <!-- claim: exact_ergun_bias_diakon_at_reported_eps_MB_pct --> %** over the Diakon rows, worse than eq. (1), in the other direction |
| 3, sensitivity | φ_s, assumed 1 from six words of prose | the headline is **-25.2 <!-- claim: wenyu_bias_diakon_pct --> % / -17.2 <!-- claim: diakon_bias_at_phi_s_095_pct --> % / -7.8 <!-- claim: diakon_bias_at_phi_s_090_pct --> %** at φ_s = 1 / 0.95 / 0.90, and **φ_s = 0.864 <!-- claim: phi_s_that_erases_diakon_bias --> erases it**. The reference is exactly φ_s-free |
| 4, digitised | none | **not used. Nothing on this page is digitised.** |

**What the page claims, and what it does not.** The algebra says equation (2)
forces a spherical packing to ε_mf = 0.383 <!-- claim: eps_mf_forced_by_eq2_at_sphere -->. Geldart's spherical cuts, inverted through the
unapproximated balance, demand 0.409 <!-- claim: inverted_eps_mf_median_diakon -->, and that shortfall accounts for the bias
arithmetically. **That inversion is not independent evidence**: it is a change of
variables on the same measured `U_0`. Its one falsifiable element is whether the
demanded voidage is physically sensible, and the page's own dataset says it may
not be: at the voidage Geldart reports, the unapproximated balance overpredicts
his velocities by more than eq. (1) underpredicts them. Sphericity cannot rescue
it, because the exact Ergun `u_mf` written in `d_sv` is φ_s-independent
(demonstrated in the notebook). **So the page does not claim that the bias belongs
to equations (2) and (3) rather than to the Ergun equation.** Three readings are
left open.

**Break tables.** Every claim has a break row built from a call that injects the
defect: wrong readings of the derivation, mis-transcribed Table 1 entries, unit errors, a wrong velocity column,
shifted rows, wrong constants and a wrong sphericity.

**Null baselines and residual independence.** A constant velocity scores 81 <!-- claim: constant_mean_velocity_mad_diakon_pct --> % MAD,
so a MAD in the teens is the floor; refitting equation (1)'s single viscous
constant to the Diakon rows leaves 14.7 <!-- claim: refit_one_parameter_mad_diakon_pct --> %, so 58 <!-- claim: eq1_error_surviving_refit_pct --> % of its error survives any
one-parameter correction. No standard error is quoted on that refit, because the
constant each row demands rises with size (Spearman ρ = +0.833 <!-- claim: residual_size_rank_correlation -->, p = 0.010 <!-- claim: residual_size_rank_pvalue -->)
though not monotonically.

## What is not here

- **The inertial half of equation (1) is untested against anything measured.**
  The largest (N_Re)_mf on any row is 0.82 <!-- claim: max_measured_re_mf -->, where the inertial term contributes 1.22 <!-- claim: max_inertial_share_on_measured_rows_pct --> %;
  a factor-four change in the turbulent constant moves the reported bias by
  0.89 <!-- claim: inertial_constant_scaled_max_shift_pts --> points at most.
- **No test of the shape-factor claim.** Geldart's two catalysts have no printed
  sphericity, so `d_sv` is substituted for `d_p` on those rows, which understates
  the predicted velocity by an unknown amount. The headline is taken on the
  Diakon rows.
- **No measurement of φ_s on the Diakon rows either; the headline rides on it.**
- **No resolution of the voidage contradiction.**
- **No digitisation of Figure 4.**
- **No pymrm operator.** The one pymrm call is `pymrm.checks.find_roots`.

## Corrections

- 2026-11-03: said that in the turbulent limit at ε_mf = 0.42 eq. (1) is "1.9 %
  low", that A1.1's refit moves the prediction "+1.55 %" there, and that the
  voidage approximation outweighs the friction-law uncertainty "1.21x" ("1.19x"
  against the printed constants). Those were values at the top of the grid "1e-2 < N_Ga < 1e10", not at the limit. The turbulent limit gives -1.8 <!-- claim: dev_vs_exact_ergun_turbulent_eps042_pct --> %,
  +1.56 <!-- claim: ergun_constant_share_turbulent_pct --> %, 1.17 <!-- claim: voidage_over_ergun_ratio_turbulent -->x and 1.16 <!-- claim: voidage_over_ergun_ratio_turbulent_vs_printed -->x.
- 2026-11-03: said the approximations "degrade the viscous term far more than
  the inertial one", that the viscous end "is where the correlation is at its worst", so that "the deviations reported below are a lower bound on its quality", and
  that the approximation is "at its best" at the turbulent end. That holds for
  ε_mf above the crossover at 0.392 <!-- claim: viscous_turbulent_flip_eps -->, not below it, where the turbulent end
  is worse: the viscous end is the worse one at 3 <!-- claim: voidages_with_viscous_end_worse --> of the 5 voidages quoted.
- 2026-11-03: said the correlation "should be expected to be least accurate for
  spherical particles, which is checkable, and is checked below". The per-powder comparison below cannot check it: the spherical Diakon rows have a MAD of 25.2 <!-- claim: wenyu_mad_diakon_pct --> % against 43.9 <!-- claim: wenyu_mad_fresh_pct --> % fresh and 31.6 <!-- claim: wenyu_mad_spent_pct --> % spent catalyst, but the catalysts carry an unknown
  sphericity and an approximate density, so the comparison neither supports nor refutes it.
- 2026-11-03: said ρ_f moves the reported bias "0.28 for a decade error".
  That run left ρ_f unchanged in the density difference. With ρ_f wrong
  everywhere it is used the shift is 0.35 <!-- claim: rho_f_decade_bias_shift_pts --> points; dropping it from the
  density difference shifts it 0.08 <!-- claim: rho_f_dropped_bias_shift_pts --> points.

## Related

[`A1.1`](../A1.1-ergun-pressure-drop/) (the friction law this is a root of, and
the source of both the printed and the refitted constants),
[`A1.7`](../A1.7-geldart-classification/) (the measured dataset, and the rival
`u_mf` expression inside Geldart's A/B boundary), `A1.5` (Richardson-Zaki, the
same bed past incipient fluidisation), `E1.2` and `E2.1` (which consume `u_mf`
as a parameter).
