# B1.4 — adversarial verification

Verifier pass, 2026-08-01. Source re-read: `~/papers/pymrm-gallery/Weisz1962-nonisothermal-effectiveness-CES17-265.pdf`
at 600 dpi (journal pp. 265, 267, 268, 269, 272, 274). Weisz & Prater (1954) not
on disk, not consulted — same position as the builder.

Everything below was checked by re-deriving it, not by reading the notebook.
The page was also re-executed end to end (192.7 s): **every printed number
reproduces bit-for-bit**; the only diff across the whole stdout is the timing
line in cell 9 (2.938 s → 2.963 s). No prose/output drift found in the markdown.

---

## Findings

### F1 (HIGH, CONFIRMED) — "the folds all sit at Phi well above 1" is inverted

Cell 27 prints:

> The equation Phi = 1 itself has a unique solution in 32 of the 32 cases
> (the folds all sit at Phi well above 1) … What it loses where Phi folds is
> the ability to say what eta IS once the verdict is 'not safe'.

Both parentheses are wrong, in the dangerous direction. Turning points of Φ
along the branch, from the page's own `Branch` objects and reproduced with an
independent shooting formulation (`w = s·y`, LSODA, 400 states):

| case | Φ at the turning points | all < 1 ? |
|---|---|---|
| β=0.6, γ=30 | 0.5296, 0.4091 | yes |
| β=0.3, γ=40 | 0.8054, 0.7703 | yes |
| β=0.4, γ=40 | 0.5741, 0.3826 | yes |
| β=0.6, γ=40 | 0.3682, 0.2243, 0.2338, 0.1005 | yes |

Every fold lies **entirely inside the band eq. (1) certifies as safe**. So do all
14 of the page's own ambiguous synthetic measurements — 14 of 14 have Φ < 1:

```
beta=0.6 gamma=40 Phi=0.2054 etas = 1.545, 2.112, 71.95
beta=0.6 gamma=40 Phi=0.3652 etas = 3.352, 3.669, 323.5
beta=0.4 gamma=40 Phi=0.5623 etas = 3.142, 4.486, 34.18
```

A measured Φ = 0.37 passes the Weisz–Prater criterion by a factor of three and
is consistent with a steady state at η = 324. That is a *second, independent*
false-negative mechanism — non-uniqueness inside the safe band, not just a loose
threshold — and it is a stronger result than the page claims. The page has the
data and denies it in prose. The conclusion that Φ = 1 has a unique root in
32/32 is correct; only the stated reason and the stated consequence are wrong.

**Fix:** replace with something like "the folds all sit *below* Φ = 1, so eq. (1)
still has a single verdict boundary — but the ambiguity is entirely inside the
region it calls safe: at β = 0.6, γ = 40 a measured Φ = 0.37 is consistent with
η = 3.35, 3.67 and 323.5." One sentence, and it strengthens Section 5.

### F2 (MEDIUM, CONFIRMED) — "below the line width" is asserted, never measured

Cell 28 and the "What pymrm adds" section both say the γ=40, β=0.3 fold is
"below the figure's line width". `docs/agent-brief.md` §7 and the F1.4/H1.7
lesson require this in pixels. Measured on the 600 dpi render of Fig. 7(d):

- decade on the Φ axis = 357.5 px (gridlines at x = 250 / 609 / 967 / 1326 / 1680
  for Φ = 0.1 / 1 / 10 / 100 / 1000);
- fold = 0.8054 → 0.7703 = 0.0194 decades = **6.9 px** (6.6 px using a
  400-state branch, which gives the slightly deeper 0.0437);
- printed stroke width of the dashed β=0.3 curve, measured over 45 rows =
  **6.0 px median at threshold 160, 7.0 px at threshold 200**.

So the fold is **≈ 1.0 line widths, not below it**. The honest statement is
"comparable to the printed line width, so Fig. 7 cannot resolve it either way",
not "invisible". For completeness I traced the drawn curve: it is monotone in Φ
over η = 4.4 → 18.2 (Φ 0.803 → 0.984), i.e. the figure as drawn shows no fold,
and it disagrees with the computed branch by ~10 px at the fold minimum — which
is within what a 1962 hand-drawn curve through IBM 704 points can be trusted to.
The figure is **silent**, which is what the page should say.

### F3 (MEDIUM, CONFIRMED) — a reference title that is not in the source

Cell 38 (References):

> Weisz, P. B. and Prater, C. D. (1954). *Interpretation of measurements in
> experimental catalysis.* Advances in Catalysis **6**, 143. … Cited exactly as
> Weisz and Hicks' reference [5] gives it.

Reference [5] as printed on journal page 274 carries **no chapter title**:
"WEISZ P. B. and PRATER C. D., in *Advances in Catalysis*. Vol. 6, p. 143.
Academic Press, New York 1954." — which is exactly how cell 2 quotes it. The
italicised title is an untraceable input on a page whose whole provenance
argument is that nothing outside Weisz & Hicks was consulted, and the sentence
claiming it is verbatim is false. (The title happens to be the real one, so this
is a traceability defect, not a factual one.) Drop the title, or mark it as
supplied from outside the source.

### F4 (LOW, CONFIRMED) — a docstring asserting what the code contradicts

`threshold_state`: *"Safe only where Phi(phi_0) is monotone, which is checked
with the shooting branch before this is used - see the kinetics cell, where it
is NOT monotone."* Cell 19 prints `fold depth in phi_0 0.00e+00, in Phi
0.00e+00` for the inhibited LH branch — it **is** monotone; the spurious
φ₀ = 5.5638 root is a non-solution, not a second branch. Same class as the two
comments handoff.md records on `B1.6`/`E1.2`.

### F5 (LOW, CONFIRMED) — the `maxfev = 3` row is converged and on the wrong branch

```
converged  : eta = 44.55494  residual 6.13e-11
maxfev = 1 : eta =  1.13848  residual 9.71e-03  flux-vol 2.03e-03
maxfev = 3 : eta =  1.15529  residual 1.08e-10  flux-vol 6.40e-13
```

Row 3 is a **genuine** solution of the discrete system (φ₀ = 0.397 sits inside
the β=0.6, γ=20 fold in φ₀, depth 0.48, so three states exist) — Newton simply
landed on the low branch. Both diagnostics are clean and η is wrong by a factor
of 39. The surrounding prose only says "the flux-volume identity DOES see a
non-converged solve"; it never notes that the row below it shows a failure mode
*neither* diagnostic catches. Say so — it is the page's own subject.

### F6 (LOW, CONFIRMED) — "1 of the four panels" is computed over three

The Fig. 7 loop runs `for g in (20.0, 30.0, 40.0)`; γ = 10 is in
`FIG7_BETA_MAX` but never evaluated, yet the printed conclusion counts four
panels. I checked γ = 10 directly: no Φ fold at β = 0.8, 1.0 or 1.2 (900-state
branch, zero negative steps in Φ). **The conclusion is right**; the loop should
either include γ = 10 or the sentence should say "of the three tested".

I also re-confirmed the boundary cases at high resolution, since γ=30 is close:
β = 0.40, γ = 30 (the largest drawn on panel c) — no fold, min ΔΦ = +9.8e-05
over 900 states; β = 0.42 — no fold; β = 0.44 — fold 7.1e-03. So panel (c) is
genuinely not reached, by about 6 % in β.

### F7 (LOW, CONFIRMED) — endothermic overstatement in the markdown

Cell 37: eq. (10a) "certifies every attainable state as safe". True for 5 of the
8 endothermic cases. At β = −0.2, γ = 10 its threshold is Φ_c = 12.18, which the
branch does reach (η there 0.290). Cell 24's printed version is careful ("the
branch never gets that far in 5 of 8 cases"); the markdown generalises past it.

### F8 (LOW) — three thresholds discard the Newton residual

`ph, et, _ = threshold_state(...)` in cell 33 (ν=0, n_u=3, 1 % rate), and
cell 16 captures `rn` but never prints it. Given the page's own Section 6.6
lesson, print them. **No spurious threshold survives**: the geometry thresholds
are cross-checked against the exact closed form to 1.1e-07, which a
non-solution cannot satisfy, and every kinetics threshold is located on the
shooting branch. Priority-6 fix judged **complete**.

### F9 (LOW, CONFIRMED, touches a second published page) — the NumJac finding also applies to `B1.1`

See "NumJac" below. `pages/B1.1-thiele-weisz-hicks/index.ipynb` uses
`self.shape = (n_r,)` in `Pellet.__init__` and `numjac_m = NumJac((N_R,))` with
`N_R = 400` in its branch-continuation cell. The builder flagged only `B1.6`.

Also: the page's explanation ("it costs n_u function evaluations per Newton
step") names the smaller term. Measured at n = 400, the constructor is
4.378 s vs 0.0006 s and the whole Newton solve is 0.0103 s vs 0.0007 s — 99.9 %
of the cost is `NumJac(shape)` itself.

### F10 (trivial) — provenance page list omits p. 274

Reference [5]/[9] were read from the literature-cited block on journal page 274.
`meta.yaml`, README and cell 38 list 265/267/268/269/272/273.

---

## Verified correct

**Equations, all read off 600 dpi renders, all verbatim**

| eq. | page | as printed | page's transcription |
|---|---|---|---|
| (1) | 265 | `(dN/dt)(1/c₀)(R²/D) < 1` | correct |
| (9) | 267 | `d²y/dx² + (2/x)dy/dx = φ₀² y exp(γβ(1−y)/(1+β(1−y)))`, `φ₀=R√(k₀/D)`, `y=c/c₀`, `x=r/R` | correct |
| (10) | 269 | `φ₀² exp[γβ/(1+β)] ≤ 1`, or `φ₀ ≤ exp[−½γβ/(1+β)]` | correct |
| (10a) | 269 | `(dN/dt)(1/c₀)(R²/D) exp[γβ/(1+β)] < 1` | correct |
| (11) | 269 | `(dN/dt)(1/c₀)(R²/D) = φ₀²η = Φ` | correct |
| Q′/Q | 268 | `1 + ½ d ln η / d ln φ₀`, "WEISZ [5] has shown" | correct |

Also verified: "with severe thermal effects defined by βγ ≳ 5" (p. 269);
"a case like γ = 20, β = 0.3 becomes conservatively 'safe' at φ₀ ~ 0.1 as
predicted by (10)" (p. 269); the Section V quote (p. 273) including the
metastable hedge; the Fig. 7 caption (p. 272) repeating eq. (11) verbatim; the
computed range "γ = 10, 20, 30 and 40, and β from 0 to +0.8 … and 0 to −0.8"
(p. 267).

**Attribution** — p. 265 reads "It was shown by WEISZ **[9]** that the
conditions" immediately above eq. (1). [9] = *WEISZ P. B., Z. Phys. Chem. 1957
11 1*. [5] = *WEISZ P. B. and PRATER C. D., in Advances in Catalysis. Vol. 6,
p. 143. Academic Press, New York 1954*, cited on p. 265 for shape-insensitivity
("various geometric forms [5, 6]") and on p. 268 for Q′/Q. Builder's report
correct in every particular.

**Three printed statements reproduced, and the check discriminates**
- exp(−0.5·20·0.3/1.3) = 0.099491 vs "φ₀ ~ 0.1". The two mis-reads the page
  prints are the plausible ones and both are far away: γβ (no 1+β) → 0.0498,
  γβ/(1−β) → 0.0138.
- βγ onset 4.86–5.61 vs "βγ ≳ 5" — reproduced independently.
- Q′/Q ∈ [0.5005, 1.0000] vs "between Q and Q/2". Weak by construction (the
  limits of d ln η/d ln φ are 0 and −1 for tanh), but it does pin the factor ½,
  and the page labels it a reproduced statement rather than a validation.

**Headline results, re-derived independently** (own formulation `w = s·y`, which
removes the 2/s singularity, integrated with Radau/LSODA — shares no code with
the page):

| quantity | page | independent |
|---|---|---|
| η* sphere / cyl / slab, first order | 0.935278 / 0.880525 / 0.694816 | 0.9352778530 / 0.8805247011 / 0.6948165381 (mpmath, 30 digits) |
| η* on V/S length | 0.5633 / 0.5989 / 0.6948 | identical |
| η* inhibited LH | 1.048325 | 1.04832535 |
| η at Φ=1, β=0.6 γ=40 | 3136 | 3136 |
| … β=0.4 γ=40 / β=0.6 γ=30 / β=0.3 γ=40 | 125.87 / 134.98 / 18.594 | identical |
| fold in Φ, four folded cases | 0.7214 / 0.3325 / 0.2249 / 0.0423 | 0.7214 / 0.3325 / 0.2249 / 0.0423 |
| fold in φ₀, all cases | 0.963 … 0.176 | identical |

η > 1 isothermally is **real, not a solver artefact**: `36y/(1+5y)²` peaks at
y = 0.2 with R = 1.8 = 1.8·R(1), so depletion accelerates the pellet.

**Sweep ranges do not flatter the conclusion.** β is truncated at 0.6 against the
paper's 0.8. Extending to 0.8 would raise the worst η eq. (1) admits and would
add folded cases (β=0.8 folds at γ=30 and 40, not at 10 or 20 — so "4 of 32"
would become 6 of 36). Both truncations cut *against* the page's claims.

**"False negative" is defined symmetrically** in `confuse` for both criteria and
both signs of β; `truth` requires every state consistent with the measurement to
satisfy |η−1| ≤ 0.05, which is the right definition where Φ is multivalued.

**Endothermic branch is well posed**: 1 + β(1−y) ≥ 1 + β = 0.6 > 0 at β = −0.4,
no singularity; T_centre/T₀ = 0.6 is extreme but inside the paper's stated range.

**Conventions**: deviations are |computed − reference|/reference throughout;
`bc` is on the outward normal with the physical equation in a comment; `nu` is
commented at every `construct_div`; operators assembled once in `__init__`.

---

## Break tests

Every check the page presents as having resolving power moves hard when broken.
The third row is the one the page never breaks itself; I did.

| check | baseline | ν wrong | other defect |
|---|---|---|---|
| shooting vs exact isothermal η | 7.13e-11 | 8.37e-01 | — |
| pymrm vs exact isothermal η, n=400 | 8.31e-06 | 6.30e-01 | 3.47e-01 (n_u=3) |
| **pymrm vs shooting, non-isothermal, n=800** | **5.02e-05** | **4.70e+00** | **7.72e+02** (exponent sign flipped) |
| Φ=1 threshold η*, pymrm vs closed form | 1.11e-07 | 2.57e-01 | 8.18e-03 (n_u=3), 1.80e-04 (n_u=20) |
| Φ=1 threshold η*, pymrm vs shooting, 5 kinetics | 1.99e-07 | 4.37e-01 | 1.50e-02 (n_u=3) |

Newton residuals over the 18 non-isothermal branch states at n=800:
3.1e-10 to 5.6e-10, i.e. **all converged** — the order-2.00 result is genuine.
The seeded iterate moves 1e-13 to 2.6e-05 from the shooting profile, which is
just the exact/discrete gap at that resolution, not a failure to iterate.

**Flux/volume identity — the proof holds.** `construct_div` is conservative, so
Σᵢ Vᵢ resᵢ telescopes to A_N (dy/du)|₁ − A₀ (dy/du)|₀ − φ² Σ Vᵢ Rᵢ, and
A₀ (dy/du)|₀ = 0 for every ν (zero area for ν>0, zero gradient for ν=0). With
V_tot = 1/(ν+1) this gives exactly Φ_flux − Φ_vol = Σ Vᵢ resᵢ / V_tot — the
volume-weighted mean Newton residual, as the page states. A wrong ν changes
`construct_div` and `dv`/`v_tot` together, so the blindness is structural, not
empirical. Measurements reproduced exactly: 2.08e-12 (none), 4.23e-12 (ν=0),
2.11e-15 (n_u=3), 9.90e-03 (1 % rate mismatch), 2.03e-03 (maxfev=1).

**NumJac.** `NumJac((n,))` on a one-field problem builds a fully dense Jacobian:
nnz/n² = 1.000 at n = 300 and n = 400, against 300 and 400 non-zeros for
`(n,1)`. At n = 400: constructor 4.378 s → 0.0006 s, Newton solve 0.0103 s →
0.0007 s, and the answers are **byte-identical** (`np.array_equal` True,
`tobytes()` equal). `pages/B1.6-prater-relation/`'s `ReducedPellet`
(`self.shape = (n_u,)`, default n_u = 400) and `pages/B1.1-thiele-weisz-hicks/`
(`Pellet.shape = (n_r,)`, plus `NumJac((N_R,))` with N_R = 400) both carry it.
Performance only.

---

## Optional strengthening (not a defect)

The page's forward solution reproduces the *peaks* of Fig. 7 to about 1 %:
γ=30 gives η_max = 85.24 / 22.81 / 5.838 / 1.673 for β = 0.4 / 0.3 / 0.2 / 0.1
against roughly 85 / 22 / 5.8 / 1.65 drawn; γ=40 gives 128.8 (off the top of the
100 axis, and the drawn curve does exit the top) and 17.34 / 2.546 against
roughly 15 / 2.5. That would be a fourth reproduced-result check and it confirms
the fold analysis is on the right curves — but reading peak *values* off a
figure is a digitisation and would need the maintainer gate, unlike reading the
β labels. Mention it qualitatively or leave it out.
