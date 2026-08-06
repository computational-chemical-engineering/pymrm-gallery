# D2.1 — adversarial verification

Verifier pass, 2026-08-01. Staged page `queue_cases/D2.1/page/`.
Source re-read independently on 600 dpi renders of
`~/papers/pymrm-gallery/vanWelsenaere1970-parametric-sensitivity-runaway-CES25-1503.pdf`, PDF pages 2, 3, 9, 10,
11, 12 (journal pages 1504, 1505, 1511, 1512, 1513, 1514). Barkelew's own paper
was not available to me either.

**Verdict: safe to publish after the five fixes listed at the end (F1–F5).**
No fabricated number, no circular agreement, no powerless headline check. The
reconstruction the whole page rests on is *correct* — I proved it by a route the
page did not use. The defects are precision, two false statements about the
page's own code, and one break-table row that cannot fail.

---

## Verdicts on the four questions asked

**(1) Is the τ/N/S reconstruction uniquely determined by what is printed?
No — one step is a modelling choice, not a printed constraint. But it is right.**

What the printed text does pin, and I checked each independently:

- The *N/S* sentence (p. 1511, verbatim on the render) forces the cooling term to
  be `N·τ` and the generation to be `S·y·e^τ`: at τ=1 the transfer rate is N, at
  τ=0 with y=1 the generation rate is S, ratio N/S. A cooling term proportional
  to τ forces τ ∝ (T−T_w) and τ(T_w) = 0.
- The *S* sentence forces S to be the adiabatic Δτ from a feed at the wall
  temperature: with N=0, dτ/dy = −S, so τ(y=0) = S. Their Eq. 32 (printed on
  p. 1511) then gives S = a(ΔT)_ad/T_w².
- The rate must be **exactly** e^τ, and this is pinned by a *third* printed
  statement the page under-uses: their own τ_m = 1 for criterion 1. The maxima
  curve is y ∝ τ e^{−τ} — peak at τ = 1 — only if the exponent is τ itself. Any
  rescaling τ→λτ moves that peak to 1/λ.

What is **not** printed anywhere: the temperature at which the Arrhenius exponent
is linearised. `a(T−T_w)/T_w²` is the tangent at T_w, but `a(T−T_w)/T_M²` and
`a(T−T_w)/(T_w·T_M)` satisfy every printed sentence and every quotation above
equally well. That is a Frank-Kamenetskii convention, not a constraint from the
paper — and the page's own break table concedes this by injecting the T_M
linearisation as a *defect*. So the page's phrasing, that the boxed definition is
"pinned by the two definitions Section 5 *does* give in words", overstates by one
step.

The reconstruction is nevertheless **correct**, and decisively so — see finding 2.

**(2) Is "the collapse drifts" a sound refutation or an extrapolation artefact?
Sound. Not an artefact, and not a refutation of anything the paper asserts.**

- Reproduced by an independent formulation. I marched the reduced system in ζ
  with Radau (the page uses an LSODA quadrature in dτ/dy — a different
  independent variable, a different integrator, a separate transcription of the
  right-hand side) and located the envelope with a 5-point d/dlnS instead of a
  2-point difference. I get τ_m = 1.28791 / 1.27635 / 1.25468 at S = 4/8/16,
  against the page's 1.2868 / 1.2746 / 1.2521.
- The drift is physical, not a differencing artefact. Richardson-extrapolated to
  h→0 the tangency falls monotonically from 1.28791 (S=4) to 1.16217 (S=200),
  −9.76 %; the page's fixed h = 0.04 reports −10.3 %. The h-bias inflates the
  drift by about 0.6 of ~9.8 percentage points — 6 % of the effect, not the
  effect. (It is still a defect: finding 1.)
- The maximum near S ≈ 4 is confirmed at converged h: 1.26963 (S=2), 1.28158
  (2.5), 1.28625 (3), 1.28783 (3.5), **1.28790 (4)**, 1.28715 (4.5), 1.28296 (6).
- Not extrapolation. Van Welsenaere & Froment write "for all *S*-values", which
  is unbounded; and they write "very close to 1·275", which already declines to
  claim a constant. The page separates "the range Barkelew drew" from the wider
  range everywhere it quotes a number, never says the paper is wrong, and frames
  the result as quantifying "very close to". That is the right framing and it
  survives.
- Resolving power caveat, in the page's favour: the drift over S = 4–32 is 0.06
  in τ_m. Reading the crosses off the printed Fig. 9 by eye gives τ_m ≈ 1.27,
  1.27, 1.29, 1.22 with a scatter of about ±0.03 — so the figure *cannot* resolve
  the drift, which is exactly why the paper could only say "very close to". The
  page's claim is therefore new information, not a contradiction of a measurement.

**(3) Is the partial-cancellation result real? Yes — exactly reproduced.**

I recomputed the decomposition from the page's own three curves over
T_w = 600–700 K:

| T_w | rate law alone<br>(crit-1 exp vs crit-1 true) | criterion alone<br>(Barkelew vs crit-1, both exp) | total | product check |
|---|---|---|---|---|
| 600 | −5.07 % | +2.79 % | −2.43 % | −2.43 % |
| 625 | −5.57 % | +4.62 % | −1.25 % | −1.25 % |
| 650 | −6.23 % | +7.00 % | +0.32 % | +0.32 % |
| 675 | −7.12 % | +9.82 % | +2.11 % | +2.11 % |
| 700 | −7.55 % | +12.54 % | +4.04 % | +4.04 % |

The two have opposite sign at every wall temperature and the total is smaller in
magnitude than both components at every one of them. The criterion half is a
purely dimensionless quantity (the ratio of critical S at the tangency to
critical S at τ_m = 1, at the same N), so it does not inherit the reconstruction;
only the rate-law half does. The claim is real and it is the page's most
interesting result.

Two qualifications. §5 *does* say Barkelew's Fig. 8 values were obtained "using
the modified rate equation" (verbatim on the render of p. 1512), so the paper is
not hiding which rate law is on which side; the page's markdown says this
correctly ("mixes rate laws … the reader cannot tell how much of the difference
is which"), the case YAML's "conflates" is looser. And the numbers above appear
**nowhere on the page** — see finding 6.

**(4) Which checks survived my own break tests?**

Re-run by me, from the page's own functions: check 1 baseline 1.655 % ✓;
exp(τ)→1+τ, envelope gone ✓; cooling sign flipped, fails ✓ (mislabelled, see
finding 4); order 2 → 1.4759 ✓; check 3 baseline 1.7703 % ✓; T_M linearisation
→ 7.177 % ✓; square dropped → 99.845 % ✓ but guaranteed (finding 4); check 4
baseline 0.2593 % ✓; n_ζ = 25 → 13.359 % ✓; ν = 1 → degenerate ✓.

New break tests I ran that the page did not:

| injected defect | check 3 | comment |
|---|---|---|
| B × 1.1 | **1.7703 %, exactly unchanged** | check 3 is *exactly* invariant to B |
| U × 1.2 (C × 1.2) | 1.661 % | nearly blind |
| c_p × ρ_g (the page's own headline trap) | 2.008 % | nearly blind — but check 5 catches it, ln K → −1.800 vs printed −2.055 |
| a × 1.05 | 123.6 % | catches |

| injected defect (pymrm route only) | check 4 |
|---|---|
| cooling sign flipped in `BarkelewTube.reaction` | 0.259 % → **inf** (diverges) |
| S 10 % wrong in `BarkelewTube.reaction` | 0.259 % → **inf** |

The second table matters: the page says check 4 is blind to the model. It is not.

---

## Findings, by severity

### 1. Every tangency number carries an untested O(h²) truncation error of 0.07–0.66 %, biased in the direction of the headline drift — CONFIRMED

`tangency(S, tm, h=0.04)` locates the envelope from a two-point difference of
τ_m/S in ln S with `h` hard-coded and never varied. Richardson extrapolation from
h = 0.01 and 0.005 (LSODA rtol 1e-10; the answer is insensitive to rtol over
1e-7…1e-11) gives:

| S | page (h=0.04) | h=0.01 | h=0.005 | converged | page error |
|---|---|---|---|---|---|
| 3 | 1.28539 | 1.28621 | 1.28625 | 1.28627 | −0.068 % |
| 4 | 1.28680 | 1.28784 | 1.28790 | **1.28791** | −0.086 % |
| 8 | 1.27458 | 1.27624 | 1.27632 | **1.27635** | −0.138 % |
| 16 | 1.25210 | 1.25451 | 1.25463 | **1.25467** | −0.205 % |
| 32 | 1.22571 | 1.22906 | 1.22923 | **1.22929** | −0.291 % |
| 64 | 1.19828 | 1.20282 | 1.20305 | 1.20312 | −0.403 % |
| 200 | 1.15455 | 1.16167 | 1.16204 | **1.16217** | −0.655 % |

Consequences: the headline **1.655 % becomes 1.575 %**, the wide-range
**3.733 % becomes 3.512 %**, and the four values quoted on the README, in
`meta.yaml`, in `models_entry.yaml` and in `agreement.json` (1.2868, 1.2746,
1.2521, 1.2257, 1.1545) are wrong from the fourth digit.

Failure scenario: the tangency is an *exactly defined* mathematical object, and
the page presents it to five significant figures as the whole point of computing
rather than drawing it. A reader who re-derives it converged gets 1.2879, not
1.2868, and has no way to know why; and CI's `check_agreement.py` will fire the
day anybody touches `h`. The bias also grows monotonically with S, i.e. it points
the same way as the page's main new result, which is the worst possible direction
for an uncontrolled error to point.

### 2. The decisive test of the reconstruction is available, was not used, and passes — CONFIRMED

Their **Fig. 10** is the four cases p⁰ = 0.017 atm, T_w = 625/626/627/628 K
computed *with the simplified rate law*. It is the one printed object that can
discriminate the dimensional τ. The page already computes the blow-up abscissae
(0.7511, 0.5398, 0.4613, 0.4090 m) and then only says "all four reach the top of
the paper's own 800 K axis".

I measured them off the 600 dpi render (frame at x = 225.5 px and 1804.0 px for
z′ = 0 and 1.0 m; the four curves cross T = 800 K at x = 1417.5, 1084.5, 970.0,
868.5 px) and re-ran the ODE for each candidate reconstruction:

| T_w | figure | τ = a(T−T_w)/T_w² | dev | /T_M² | dev | /(T_w·T_M) | dev | /T_w |
|---|---|---|---|---|---|---|---|---|
| 628 | 0.4074 | **0.4090** | +0.40 % | 0.5608 | +37.7 % | 0.5285 | +29.7 % | no runaway |
| 627 | 0.4716 | **0.4613** | −2.17 % | no runaway | — | 0.7522 | +59.5 % | no runaway |
| 626 | 0.5442 | **0.5398** | −0.80 % | no runaway | — | no runaway | — | no runaway |
| 625 | 0.7551 | **0.7511** | −0.53 % | no runaway | — | no runaway | — | no runaway |

The page's reconstruction reproduces the printed figure to within the line width;
every alternative fails qualitatively. Two supporting checks in the same
direction: their **Fig. 3** (same four cases, true rate law) peaks at ≈778 K for
T_w = 628 K against the page's computed 777.57 K; and the tangency *abscissae*
the page computes land on the crosses printed on **Fig. 9** — read off the render
at N/S ≈ 0.93, 1.46, 1.86, 2.14 against the computed 0.9447, 1.4634, 1.8573,
2.1390, and the circles at τ_m = 1 at ≈1.17, 1.64, 2.00, 2.24 against 1.1738,
1.6458, 1.9899, 2.2309.

So the page is *more* right than it claims, and its self-assessment ("the 1.275
check is blind to the groups; check 3 is the only thing that sees them") is
pessimistic. The builder's reason for declining is defensible — measuring
positions on a curve is a figure digitisation and per `AGENTS.md` needs a
maintainer review — but the page should not leave a reader thinking the test was
impossible. Recommended, not blocking: add the Fig. 10 comparison with a numbered
overlay in `review/`, or say in one sentence that the comparison exists and was
declined so as not to digitise.

### 3. Two statements the page makes about its own code are false — CONFIRMED

(a) Validation table, row 4: check 4 cannot see "anything about the model, since
both routes call the same `reaction` function." **They do not share it.**
`tau_max_phase` transcribes dτ/dy in its own local `rhs`; `BarkelewTube.reaction`
transcribes the ζ-form separately. I flipped the cooling sign inside
`BarkelewTube.reaction` only, and separately made S 10 % wrong there only: check 4
went from 0.259 % to `inf` both times. The check has real power over the model in
the finite-volume route, and the page asserts a fact about its own implementation
that is not true.

(b) "What pymrm adds" says the method "turns the envelope into the **pointwise
minimum over S**". The code comment in the diagram cell says the opposite, and is
the one that matches the code: *"The envelope IS the locus of tangency points …
a pointwise minimum over a truncated S-range would be an artefact of where it was
truncated, not the envelope."* One of the two sentences has to go.

### 4. One break-table row cannot fail; another is mislabelled — CONFIRMED

`tau = a(T−T_w)/T_w, the square dropped` enters `check3` only through the factor
`T2(Tw)` in `p0 = Sc*A*T2(Tw)/(a*B)`. `kk`, `N` and the interpolated `Sc` are all
untouched, so p⁰ is multiplied by exactly 1/T_w ≈ 1/650 and the deviation is
1 − p⁰/(650·p_ref) ≈ 99.84 % **by construction**. It is the mirror image of the
defect `handoff.md` exists to catch: a break that cannot *not* break. Presented in
the same table as the informative rows, it reads as evidence that check 3
discriminates the reconstruction, which it is not.

The row that *is* informative is the T_M linearisation, 1.770 → 7.177 %. I
reproduced it. It is a genuine ~4× degradation from a genuinely plausible
alternative, and it is the only entry in the table that tests the one step of the
reconstruction the paper does not print.

Separately, `cooling term sign flipped` is reported as "no envelope exists". In my
run the failure is a `ValueError` raised by `u_at_tau`'s bracketing, which the
page's `except (ValueError, RuntimeError)` funnels into the same string. The
conclusion (the check catches it) stands; the label does not.

### 5. Check 3's blindness is stated only vaguely; it is *exactly* blind to B — CONFIRMED

The table says check 3 cannot see "an error common to both routes". Concretely:
p⁰(Barkelew) ∝ 1/B, and p⁰(VWF, back-integrated) ∝ 1/B exactly — substituting
p = q/B removes B from their Eq. 5 and from p_M. So check 3 is *identically*
invariant under any rescaling of B: I confirmed 1.7703 % unchanged at B×1.1. It
therefore cannot see (−ΔH), ρ_b or p_B⁰. It is also nearly blind to C (U×1.2 →
1.661 %) and to the page's own headline trap, c_p multiplied by ρ_g (→ 2.008 %,
inside the noise of "agreement").

Nothing is left unguarded — the c_p trap is caught by **check 5**: with c_p·ρ_g
the recomputed ln K is −1.800 against the printed −2.055, 12.4 % off, versus
0.09 % as built. But the README calls that error "silent", so the page should say
which check is not silent about it.

### 6. The page's most novel claim is never given numbers on the page — CONFIRMED

The rate-law/criterion decomposition (table under verdict 3) exists on the page
only as a dotted line in the right-hand panel. "5.1–7.6 %" appears in the case
YAML and nowhere on the page; the criterion half appears nowhere at all. A reader
cannot check the partial-cancellation claim from what is printed. One extra
`print` of the two columns fixes it.

### 7. "Pinned" overstates by one step, and one supporting sentence is wrong — CONFIRMED

See verdict 1 for the full argument. Two concrete corrections:

- *"the two definitions Section 5 does give in words, **both of which name the
  wall temperature as the reference**"* — only the S sentence names the wall
  temperature. The N/S sentence names τ = 1 and τ = 0.
- The one unprinted step is the linearisation temperature. Say so; the honest
  version is stronger than the current one, because the T_M break row is right
  there and shows the choice is testable.

Also: the case YAML says the reconstruction is "labelled as such on the page and
in the sidecar". It is prominent on the page, in `README.md` and in `meta.yaml` —
that satisfies `AGENTS.md`. The **sidecar does not mention it**; it says "Nothing
was reconstructed by inference" (of the data values, which is true and not a
contradiction, but the YAML's claim is inaccurate).

### 8. The scope decision overstates the D2.1/D2.2 independence — CONFIRMED, integrator-facing

The YAML: *"D2.1 re-derives ONE thing D2.2 also has … the two implementations
agree to 1.4e-8 relative at 625 K without either having seen the other."*
`p0_crit_vwf` and D2.2's `back_integrate` integrate the **same** ODE (Eq. 5) from
the **same** starting point (T_M from Eq. 8, p_M from Eq. 7) with the **same**
library, differing only in `method` and `rtol`. 1.4e-8 is scipy agreeing with
scipy to tolerance; it is not evidence of independence. Beyond that function,
D2.1's `Tube`, `T_M`, `C_of_R`, `kexp`, `A`, `B`, `C_BASE` and the whole
BC/operator block are copied from D2.2 down to the comments.

None of that is a defect — `AGENTS.md` explicitly asks builders to copy the
nearest page. But the integrator should judge scope on content, not on the
"one shared function" claim. **On content, BUILD is right**: the dimensionless
reduction, the computed envelope and tangency locus, the drift, the
criterion/rate-law decomposition and the reaction-order study are all genuinely
absent from D2.2, and D2.2 has no similarity claim to test.

### 9. Two prose numbers are not printed by any cell — CONFIRMED, minor

README's *"The two cross near 645 K"* — correct (I get 645.74 K by root-finding
on the deviation) but not printed on the page; the sweep only shows 600/625/650/
675/700. And "5.1–7.6 %" (finding 6). The notebook markdown is otherwise
commendably free of hard-coded numbers, and every number in it that I checked
(1.48 %, exactly 100 %, "a quarter of N/S at S = 4 to a percent at S = 200"
against 24.26 % and 1.22 %) matches the printed output.

---

## What I checked and could not break

- **Parameter reading, re-derived from the render.** A = 6165.0, C = 40965.3,
  ln K = −2.0570 against the printed −2.055; (B/A)·0.0125 = 521.09 K against the
  printed 521.09; T_M(625) = 656.62 against 656.6; a/625 = 21.818 against
  21.818. Confirmed on the 600 dpi render: `c_p = 0·323 kcal/m³ . °C` (volumetric
  — using c_p·ρ_g gives ln K = −1.800), `(−ΔH) = 307·000`, `u = 3·600`,
  `a = 13·636 °K⁻¹`. The mid-dot in this paper is *both* a decimal point and a
  thousands separator; `a = 13636` is forced by their own printed
  t_w = a/T_w = 21.818, and the printed unit `°K⁻¹` is the paper's own error.
- **The three §5 quotations** are verbatim on the renders of pp. 1511–1512,
  including *"very close to 1·275"* and *"No special meaning could be assigned to
  the corresponding trajectories, however."* The p. 1514 notation list gives only
  "*N* dimensionless group defined by Barkelew" — and the same one-liner for *S*
  and *τ* — so the reconstruction really is unavoidable. The page's README quotes
  only the *N* line; all three are missing definitions.
- **The 0.01651 reproduction is not circular.** Every input is a printed
  constant; T_M is their Eq. 8, p_M their Eq. 7; nothing anywhere is fitted to
  0.01651 or to anything downstream of it. It is, however, the same check D2.2
  already carries, so it is not new evidence.
- **τ_m = 1 for criterion 1** is a real derived identity (at a maximum
  y = (N/S)·τe^{−τ}, peak at τ = 1), matching their printed statement — a
  legitimate route-2 backing.
- **The ordering claim** holds for all 13 S, 34.3 % → 1.2 %, and is consistent
  with the circles sitting right of the crosses on the printed Fig. 9.
- **Determinism.** Re-executed end to end in 22.9 s. 156 output lines; exactly
  two differ from the staged notebook, both wall-clock timings. `agreement.json`
  rewrites byte-identical. No continuation chain, all root-finds bracketed with
  the sign change asserted.
- **Conventions.** `NumJac((n_z, 2))` — 2-D shape, last axis the field index,
  correct under the updated `AGENTS.md`; no bare 1-D shape and no
  `axes_diagonals` anywhere on the page. Boundary conditions are the
  outward-normal 2-tuple with the physical equation written in a comment beside
  each; `nu=0` in `construct_div` with the geometry stated. Deviation convention
  is (computed − reference)/reference everywhere I looked. Tier 6 is stated in
  `meta.yaml`, `README.md`, both sidecars, `models_entry.yaml` and on the page,
  and nothing is described as experimental or as validated against measurement.
- **`review/` empty by design** is correct: no curve was digitised by the
  builder. (My own Fig. 9/Fig. 10 measurements above are verification work, not
  page content.)

---

## Fixes required before publishing

- **F1.** Converge the tangency. Use `h = 0.005`, or Richardson-extrapolate two
  values of `h`, and re-report every τ_m, both mean deviations, `agreement.json`
  and the numbers echoed into `README.md`, `meta.yaml` and `models_entry.yaml`
  (1.575 % / 3.512 %; 1.2879, 1.2764, 1.2547, 1.2293; 1.1622 at S = 200). State
  the residual truncation. *(finding 1)*
- **F2.** Remove both false statements about the page's own code: the "pointwise
  minimum over S" sentence in "What pymrm adds", and check 4's "both routes call
  the same reaction function". Replace the latter with what is measured — a sign
  flip or a 10 % error in `BarkelewTube.reaction` alone sends check 4 to `inf`.
  *(finding 3)*
- **F3.** Relabel the `square dropped` break row as the guaranteed scale check it
  is (the deviation is 1 − 1/T_w by construction, so it tests nothing), and
  relabel the cooling-flip row as "bracketing fails", not "no envelope exists".
  *(finding 4)*
- **F4.** Replace "pinned by the two definitions" with an accurate account: the
  reduced system and the reference T_w are forced by print (and by their own
  τ_m = 1, which forces the exponent to be τ itself); the linearisation
  temperature is the one step that is a modelling choice. Fix "both of which name
  the wall temperature as the reference". *(finding 7)*
- **F5.** Print the rate-law and criterion deviations as numbers, not only as a
  plotted line — the partial-cancellation claim is the page's most novel result
  and currently cannot be checked from the page. *(finding 6)*

## Recommended, not blocking

- **R1.** Add the Fig. 10 comparison (finding 2) — it is the only decisive
  evidence for the reconstruction and it passes to 0.4–2.2 %. If that counts as
  digitisation, put a numbered overlay in `review/` for maintainer review; if
  not, at minimum say the comparison exists and was declined.
- **R2.** State check 3's exact invariance to B and its weak sensitivity to C,
  and note that the c_p × ρ_g trap is caught by check 5, not check 3.
  *(finding 5)*
- **R3.** Correct the case YAML for the integrator: the 1.4e-8 D2.2 agreement is
  scipy vs scipy, the code overlap is much larger than one function (and that is
  fine), and the reconstruction is not labelled in the sidecar. Scope call
  BUILD is nevertheless correct. *(findings 7, 8)*
- **R4.** Print the 645.7 K crossing rather than asserting it in the README.
  *(finding 9)*
