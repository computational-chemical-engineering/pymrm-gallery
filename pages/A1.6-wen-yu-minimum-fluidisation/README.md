# A1.6 — Wen and Yu's minimum fluidisation velocity

One equation, two constants, and a promise: no voidage, no shape factor. This
page recovers the constants from the force balance they came from, computes the
particle they secretly describe, splits the resulting error between the
approximation and the friction law underneath it, and puts the correlation
against 21 minimum fluidisation velocities measured seven years later in another
laboratory.

- **Structures:** `S3` — algebraic closure. **No pymrm operator appears**, and the
  page says so where a reader would look for a solver, as `A1.1` and `F1.4` do.
  **Changed from the catalogue's `S1`**, deliberately, to match `A1.1` and `A1.7`;
  `docs/catalog-A-foundations.md` still says `S1` and needs the same edit.
  `check_metadata.py` will not catch this — it compares `meta.yaml` against
  `models.yaml` only.
- **Runtime:** ~4 s
- **Data tier:** 2 — a table printed in another paper, whose values are that
  author's own laboratory measurements.

## Source

**Wen, C. Y. and Yu, Y. H.**, *A Generalized Method for Predicting the Minimum
Fluidization Velocity*, A.I.Ch.E. Journal **12**(3) 610–612 (May 1966),
[doi:10.1002/aic.690120343](https://doi.org/10.1002/aic.690120343) — read
directly, on 600 dpi renders of journal pages 610–612.

**The catalogue citation for `A1.6` was wrong, and usefully so.** It read
"Chem. Eng. Prog. Symp. (1966)". That is a real paper — *Wen, C. Y. and Yu, Y. H.,
Chem. Eng. Progr. Symposium Series No. 62, **62** (1966)*, which is this
communication's own reference 23 and contains the derivation — but it is **not on
disk and was not consulted**. Both are recorded in `meta.yaml`; they must not be
conflated.

**Two other origins cited but not consulted.** The communication attributes its
friction law to *Ergun, S. and Orning, A. A., Ind. Eng. Chem.* **41**, 1179
(1949) and **prints neither that equation nor the derivation**. The Ergun
constants used here are the 1952 pair, read on
[`A1.1`](../A1.1-ergun-pressure-drop/). Davies and Richardson (1966), used for
comparison, comes through Geldart, as on `A1.7`.

**A trap in this PDF.** Its first page of *extracted text* is a different
article — a communication on activity coefficients occupying the top of journal
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

## What the page checks, in the brief's order

| route | check | result |
|---|---|---|
| 2, identity | eqs. (2) and (3) put into the Ergun force balance | `24.5 Re² + 1650 Re = N_Ga`, root gives **C₁ = 33.6735 (−0.079 %)** and **C₂ = 0.040816 (+0.040 %)** |
| 2, identity | eqs. (2) and (3) solved *simultaneously* | one admissible root: **φ_s = 0.669, ε_mf = 0.474** — not a sphere. Forced to φ_s = 1 they disagree: 0.383 vs 0.415 |
| 2, identity | Table 1's rows pooled against its own overall | Narsimhan **46.24 vs a printed 46**; Equation (1) 36.47 vs a printed 34, which **cannot close on the printed row counts** (negative variance, −1572) and is reported, not repaired. The stronger claim — that the column is unreachable — is *not* made: an unconstrained split of the 284 can reach 34 |
| 2, cross-page | A1.1's refit of Ergun's constants pushed through | prediction moves **−1.27 % (viscous) to +1.55 % (turbulent)**. Like for like against the voidage approximation: **22.7× at the viscous limit, 1.21× at the turbulent limit** — the two are the same size at the coarse end |
| 3, measurement | eq. (1) vs Geldart's 21 measured `U_0` | bias **−17.4 %**, MAD 33.3 %, rms **38.1 %** against their claimed 34 % — but as an s.d. *about the mean* it is **33.9 %, an exact match**, and the paper defines neither statistic. Both printed |
| 3, measurement | the 8 spherical Diakon cuts, taking `d_sv = d_p` | **low on every row, by 25.2 % on average (−5.5 % to −54.5 %)**; Davies & Richardson −0.75 % with MAD 15.8 % |
| 3, inversion | the voidage those 8 cuts demand of the exact balance | median **0.409** — in the same band as the 0.386/0.40/0.42 the paper quotes, above the 0.383 it uses. **A reparameterisation of the bias, not an independent closure** |
| 3, contrary evidence | the voidage Geldart himself **reports** for those cuts | ε_MB = **0.444 with no inference** on the two rows printing `U_MB = U_0`, `H_MB/H_0 = 1.000`. The exact balance there is **+47.5 % and +61.1 %** high (**+58.4 %** over all eight) — *worse* than eq. (1), in the other direction |
| 3, sensitivity | φ_s, assumed 1 from six words of prose | headline is **−25.2 % / −17.2 % / −7.8 %** at φ_s = 1 / 0.95 / 0.90, and **φ_s = 0.864 erases it entirely**. The reference is exactly φ_s-free |
| 4, digitised | — | **not used. Nothing on this page is digitised.** |

**What the page claims, and what it does not.** The algebra says equation (2)
forces a spherical packing to ε_mf = 0.383. Geldart's spherical cuts, inverted
through the unapproximated balance, demand 0.409, and the literature values Wen
and Yu themselves quote say 0.386–0.42; that shortfall accounts for the −25 %
bias arithmetically. **That inversion is not independent evidence** — it is a
change of variables on the same measured `U_0`, sharing φ_s = 1, the same sizes
and gas properties and Ergun's own k₁, so the gap *is* the deviation by
construction. Its one falsifiable element is whether the demanded voidage is
physically sensible, and **the page's own dataset says it may not be**: Geldart
reports a bed voidage of 0.444 on the two coarsest cuts with no inference
required, and at that voidage the unapproximated balance overpredicts his
velocities by +58 % — worse than the correlation it is being used to judge.
Sphericity cannot rescue it, because the exact Ergun `u_mf` written in `d_sv` is
exactly φ_s-independent (demonstrated in the notebook). **So the page does not
claim that the bias belongs to equations (2) and (3) rather than to the Ergun
equation.** Three readings are left open and none is resolved.

**Three break tables and one resolving-power cell.** Wrong readings of the
derivation move the constants 19.9–250 %, against 0.079 % for the reading used —
but A1.1's refit moves them 4.4 %, so the check resolves 150 from 180 or 154 and
**not** from 151.9. The pooling identity moves to 39.5/48.5/85.0/39.7/59.6 under
five mis-readings. Every link of the measured comparison's unit chain moves the
bias by hundreds of per cent — except ρ_f, which moves it by 0.08 points and
about which the page therefore claims nothing, and except **φ_s, which is not a
unit conversion but the page's load-bearing assumption and moves the headline by
8 points at 0.95 and 17 at 0.90**.

**Null baselines and residual independence are printed beside the deviations.** A
constant velocity scores 81 % MAD, so a MAD in the teens is the floor; refitting
equation (1)'s single viscous constant to the eight rows leaves 14.7 %, so 58 %
of its error survives any one-parameter correction. **No standard error is quoted
on that refit**, because the constant each row demands rises with size
(Spearman ρ = +0.833, p = 0.010, n = 8) — *rises*, not monotonically; the printed
sequence has two reversals.

## What is not here

- **The inertial half of equation (1) is untested against anything measured.**
  Every row is a group A powder at `(N_Re)_mf < 0.83`, where the inertial term
  contributes at most 1.22 %; a factor-four change in the turbulent constant
  moves the reported bias by under 0.9 points. Section 3's finding that the
  approximation is at its *best* there is algebra, not measurement.
- **No test of the shape-factor claim**, which is the correlation's main selling
  point. Geldart's two catalysts have no printed sphericity, so `d_sv` is
  substituted for `d_p` on those 13 rows; since `d_sv = φ_s d_p` this understates
  `d_p` and the predicted velocity, by an unknown amount. Direction stated, no
  correction applied, headline taken on the 8 Diakon rows.
- **And no measurement of φ_s on those 8 rows either — the headline rides on it.**
  φ_s = 1 for Diakon is read off Geldart's phrase "a plastic moulding powder
  having spherical particles"; he prints no number. The prediction scales as
  (d_sv/φ_s)², the reference does not move with φ_s at all, so the whole
  sphericity risk lands on the reported bias. The Check 4 break table prints it.
  The measured verdict is a joint test of the voidage approximation *and* of
  Diakon's sphericity, and nothing here separates them.
- **No resolution of the voidage contradiction.** The unapproximated Ergun
  balance is unbiased at the voidage the rows demand (by construction) and +58 %
  biased at the voidage Geldart reports. Either equation (2)'s voidage is too low,
  or Geldart's `U_0` sit below the Ergun balance at his own reported voidage, or
  Diakon is not perfectly spherical. The page prints all three and picks none.
- **No digitisation of Figure 4**, and so no independent contact with the 284
  points every statistic in the paper is computed from.
- **No pymrm operator.** Both routes are closed-form roots; inventing a PDE would
  obscure the comparison.

## Related

[`A1.1`](../A1.1-ergun-pressure-drop/) (the friction law this is a root of, and
the source of both the printed and the refitted constants),
[`A1.7`](../A1.7-geldart-classification/) (the measured dataset, and the rival
`u_mf` expression inside Geldart's A/B boundary — this page is the comparison
A1.7 named as its own obvious next move), `A1.5` (Richardson–Zaki, the same bed
past incipient fluidisation), `E1.2` and `E2.1` (which consume `u_mf` as a
parameter).
