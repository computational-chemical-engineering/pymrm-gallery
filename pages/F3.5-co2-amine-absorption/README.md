# F3.5 — CO2 absorption into amine-promoted carbonate solutions

Why adding 3 mol% of amine to hot potash multiplies the CO2 absorption flux —
answered by solving the film with every reversible reaction in it, the model
that ended the shuttle-versus-homogeneous-catalysis argument.

- **Structure:** `S3` (1D steady BVP) + `S10` (instantaneous-equilibrium
  constraints)
- **Reference:** Bosch, Versteeg & van Swaaij (1989), Chem. Eng. Sci. 44(11)
  2735-2743, doi:10.1016/0009-2509(89)85216-9. (The catalogue's original
  citation, "Versteeg & van Swaaij (1988)", is the group's kinetics paper —
  the source of this paper's correlation (11) — not the model paper.)
- **Runtime:** ~105 s

## Results

With Table 1 exactly as printed, all twelve enhancement factors come out
8.4–21.4% high and the twelve fluxes 10–92% high, systematically. The paper's
own identity J = E kL (ci − cb) implies a bulk dissolved CO2 1.42–1.49x above
what the printed constants give, and the model's structure localises the
discrepancy to the single constant K3 = Kw/Kc2 (whose rescaling
simultaneously fixes the bulk CO2 and, through the bulk OH−, the Hatta
number). One scalar, s = 1.4585, gives:

| comparison | deviation (model − paper)/paper |
|---|---|
| **unpromoted flux J_u (4 values) — headline** | **max 3.8%** |
| unpromoted E (4 values) | max 0.6% — see caveat |
| promoted E, DEA + HDA (8 values) | max 7.6% |
| promotion factors F (8 values) | max 5.6% |

**Read the E_u number with care.** The scalar is fitted to
cb = ci − J/(E kL), so it consumes the printed fluxes *and* the printed
enhancement factors. The identity does not admit one scalar exactly: the four
per-condition scalars span 1.419–1.490, and that residual lands in the flux,
because the fitted cb sits in E's denominator and makes E 2.6x less sensitive
to s than J is. The four E_u deviations (+0.2/+0.6/+0.6/+0.6%) are one offset
repeated, not four independent agreements; refitting to E_u alone gives
1.4735, 1.0% away. The flux is the honest metric.

**Evidence for, using no enhancement factor at all:** the unpromoted bulk CO2
read off Figures 2, 3 and 4 implies s = 1.419, 1.405 and 1.510 (mean 1.445),
and the α = 0.2 solve reproduces Figure 2's HDA profile (374 → 485 against
369 → 484 measured). Figure readings are a cross-check *pending maintainer
visual review*; no dataset is shipped. They were corrected on 2026-08-02
after a first review pass found two curve-identity errors (Figure 2's
promoted CO3²⁻ missed, Figure 3's two carbonate curves fused into one); the
corrected reading passes two independent charge balances — worst deviation
0.11% unpromoted, 0.20% promoted — and is itself awaiting confirmation. The three s
values above are unchanged by the correction: the curves that were wrong are
the promoted carbonates, which enter neither the loading identification nor
the bulk CO2 that gives s.

**Evidence against, stated plainly on the page:** the paper's only two
out-of-sample numbers both prefer the constants exactly as printed —
desorption F = 3.72 (−2.3%) printed against 3.35 (−12.0%) reconstructed, and
a driving-force reduction of 34.7% against 39.1% for the paper's "about one
third". The paper's own record is internally inconsistent here, since its
Figure 4 supports the reconstruction at the same condition.

Validation, ranked by what it can actually catch. Primary: the
pseudo-first-order closed form to 3.5e-5 over 0.1 ≤ Ha ≤ 300 (with grid
doubling 1.3e-5 and k_inst 2.2e-5), Van Krevelen–Hoftijzer within 0.1%, and a
Higbie penetration variant within **11.3 %** at `n_t = 640`, extrapolating to
11.6 %. Consistency checks, true by construction: physical-absorption limit
5e-11, a dimensionless cross-assembly sharing grid and operators 1e-8,
reversible→irreversible collapse 1e-5, carbon- and amine-flux closure 6e-10 and
1e-10. Electroneutrality (2.8e-11) and charge flux (1.1e-9) are algebraically
guaranteed by the equal ion diffusivities and are kept as regression guards only.

## The transient result was never time-refined, and the headline was 9 % low

Every steady result on this page is grid-refined. The Higbie penetration
variant — a transient solve, backward Euler on a geometric time grid — was run
once at `n_t = 160` and its result published. It was not converged:

| `n_t` | quadrature bias | film-vs-penetration max |
|---|---|---|
| **160** (published until 2026-08-05) | +1.42 % | **0.10488** |
| 320 | +0.72 % | 0.11032 |
| **640** (published now) | +0.37 % | **0.11301** |

Successive differences 5.44e-3 and 2.69e-3 give an **observed order of 1.02** —
clean first-order backward-Euler time error, not a quadrature property — and
Richardson extrapolates to **0.1156**. The withdrawn value was 9 % below that,
outside the 5 % tolerance `check_agreement.py` applies; the value now published
is 2.3 % below it. All three of the other knobs were already converged: `n_x`
400 → 800 and the domain factor 12 → 16 each move the number by 1–2 % of a
single `n_t` step, so `n_t` was both the knob that mattered and the only one
nobody had swept.

**The conclusion is unchanged and the number is not.** The penetration model
gives the *larger* enhancement, so the reproduction is still not an artefact of
the film idealisation and the film variant — which the paper's section 4.2
names — still matches Table 2 better. The size of the idealisation's effect was
understated by about one percentage point. The withdrawn value is published
beside the corrected one as `pen_vs_film_max_nt160_withdrawn`.

## A check that was labelled for the wrong reason

The carbon- and amine-flux closures were described as "genuine — they mix the
three different diffusivities". Section 9 of the notebook injects that and it is
not true: every reaction conserves carbon, so ∇·Σγⱼ**J**ⱼ = 0 identically for
*any* set of Dⱼ. Halving `D(HCO3⁻)` and doubling `D(CO2)` leave both closures at
1e-10 while E moves 11 % and 16 %. What they do catch is a wrong entry in the
**stoichiometry matrix**, which moves them by eight to ten orders. They are a
good check of a different thing than was claimed, and are relabelled.

## Data

**Provenance tier 6 — reproduction, not validation.** The paper tabulates no
measurements; its Table 2 is the authors' own numerics and its only
experimental numbers print as "±4" and "±6" (almost certainly "≈4" and "≈6").
Both tables are transcribed from 600 dpi page renders (the scan's OCR drops
exponents, and gives F_HDA = 3.74 where the print shows 3.76) with the
transcription checks described in the sidecars.
