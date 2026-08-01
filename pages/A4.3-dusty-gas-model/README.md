# A4.3 — The dusty gas model

Pore diffusion in its *n*-component form: the matrix that one effective
diffusivity per species is an approximation of, plus the viscous term that
page `A4.4` set to zero.

**Source.** Everything on the page is read from **Krishna, R. & Wesselingh,
J. A. (1997)**, "The Maxwell–Stefan approach to mass transfer", *Chemical
Engineering Science* **52**(6) 861–911,
[doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7),
which is on disk. Every equation and constant was read off a **600 dpi render**
of the printed page on 2026-08-01, never from the PDF text layer — that layer is
an Acrobat 3.0 Capture OCR and mangles eq. (86)'s prefactor and eq. (109)'s
stoichiometric subscript. The **origins** of the model — Mason & Malinauskas
(1983), Jackson (1977), Knudsen (1909), Bosanquet's 1944 report — are cited by
the review and were **not consulted**; none is on disk or reachable. So are
Kaza & Jackson (1980), Haynes (1978) and Remick & Geankoplis (1974), each known
only through a sentence in the review.

Files:

- `build_page.py` — writes `index.ipynb`; prose and code live here.
- `index.ipynb` — the page; executes end to end in about 25 s.
- `data/krishna-wesselingh-1997-dgm-constants.csv` — the structural constants of
  eqs. (83)–(85), (87), (91), (101), (103), (107), plus the four exponents the
  review states in words on its p. 887, with a provenance sidecar. Molar masses
  are IUPAC and flagged as not from the review.
- The pair diffusivities are **not** re-transcribed: they are loaded cross-page
  from the published `A4.2` dataset, which holds the review's own p. 872 values.

**The result.** Three things the matrix form says that no scalar effective
diffusivity can express, one that quantifies a sentence the review leaves
qualitative, and one theorem.

1. An **inert species** in a pellet running A → 2B carries exactly zero flux —
   conservation, not a prediction — with a concentration that varies by 22.5 % at
   Kn = 1 and 36.6 % under bulk control. A scalar closure **with a finite
   positive diffusivity** makes that variation exactly zero at every Knudsen
   number, which is run on the page (0.00e+00) rather than asserted. The premise
   matters: the review's own eq. (109) degenerates to D → 0 for an inert and
   predicts nothing there.
2. **Uphill diffusion has no critical pore size.** The pore diameter at which
   N₂'s flux reverses is a function of the near-end composition, and moves by a
   factor **4650** — from 4246 nm to 0.91 nm — across compositions that all lie
   inside the uphill window at 1 mm (x_N2 ∈ [0.4823, 0.5002]). At equal end
   compositions N₂ is uphill at **every** diameter from 1 nm to 1 mm. What pore
   size controls is the *width of the composition window*: its half-width falls
   from 1.77 × 10⁻² at 1 mm to 1.09 × 10⁻⁴ at 1 nm, a log–log slope of 0.98 in
   d₀, and it never closes. So uphill diffusion is not a macropore phenomenon
   that stops — it is one that needs an ever finer-tuned composition, and no
   scalar closure shows it at any pore size.
3. **Jackson's** *p*₀ = √ν_B *p* is the corner case of
   *p*₀/*p* = 1 + (ν_B √(M_B/M_A) − 1)·x_A,s, needing pure A at the surface,
   mass-conserving stoichiometry **and** Knudsen control. Under bulk control the
   rise nearly vanishes (1.0023 against 1.4141 at Kn = 0.01). The −5.0 × 10⁻⁵
   residual at Kn = 10⁴ is exactly 1/Kn (−5.0 × 10⁻⁶ at 10⁵, −5.0 × 10⁻⁷ at 10⁶);
   pushing φ does *not* shrink it (−5.09 × 10⁻⁵ at φ = 200).
4. **"Neglect of the viscous flow contribution is not very serious (Haynes,
   1978)"** becomes a screening group assembled from the review's own eqs. (85),
   (91) and (106): D_visc/Ðᵉ_iM = 3*p*d₀/(32 η v̄_i), linear in pore size — 2.64 ×
   10⁻⁴ for H₂ in a 1 nm pore and 0.26 in a 1 µm one. **That group cannot be
   swept alone**: eq. (85) makes Kn a function of the same d₀, running the other
   way (138.9 at 1 nm, 0.139 at 1 µm). Tying both to d₀, viscous flow changes the
   pellet-centre pressure by 0.009 % at 1 nm and 0.62 % at 1 µm, and **never by
   more than 0.63 %** anywhere from 1 nm to 100 µm — the maximum sits at 562 nm,
   Kn = 0.25. So the review's "not very serious" holds *everywhere*, not only for
   catalysts. Freezing Kn at 10⁴ instead — the natural but inconsistent sweep —
   gives 6.9 % at the 1 µm viscous group, an elevenfold overstatement of a
   pressure rise no 1 µm pore has; the page prices that mistake explicitly.
5. **[Bᵉ] is invertible for a reason, not by luck.** Σᵢ Bᵢⱼ = 0 identically for
   the wall-free matrix because Ð is symmetric, so (1,…,1) is an exact left null
   vector for any composition and any n ≥ 2. Adding the wall term makes every
   column sum 1/Ðᵉ_jM > 0 and the matrix strictly column-diagonally dominant,
   hence nonsingular by Levy–Desplanques. Checked on the assembled `b_matrix`
   for n = 2…6 (null vector to 1.3 × 10⁻¹⁶, dominance margin ≥ 0.92 of the
   diagonal), and V1's break table injects `knudsen=False` and gets an exception
   rather than a wrong answer.

**Validation.** **Tier 6** — nothing is compared with a measurement.

| check | result |
|---|---|
| independent collocation solution (shares no transport code) | 8.7 × 10⁻⁶ worst over six decades of pore size, converging at order 2.00 |
| eq. (103) recovered, viscous term included | 1.0 × 10⁻¹³ |
| collapse to `A4.4`'s published eq. (109) | 1.2 × 10⁻¹³ |
| Graham's law, predicted not imposed | 2.5 × 10⁻¹¹ (algebraic identity, labelled) |
| pressure uniformity across the capillary | 5.4 × 10⁻¹⁶ |
| free-molecule / Maxwell–Stefan asymptotes | 1.8 × 10⁻⁴ / 7.9 × 10⁻⁶ |
| the four p. 887 exponents | worst departure 0.021 |
| eqs. (105), (107), (103) derived from eq. (101) symbolically | identically 0 |
| grid order / 100-cell error | 2.03 / 3.9 × 10⁻⁷ |
| [Bᵉ] null vector / dominance margin | 1.3 × 10⁻¹⁶ / 0.92 of the diagonal |
| worst relative Newton residual over all 639 non-break solves | 1.8 × 10⁻⁸ (V4B at d₀ = 100 mm) |

The collocation route shares `knudsen_D` and `pair_matrix` with the solver — the
two helpers that turn eqs. (83) and (85) into numbers — so it cannot catch a
mistake inside those. Eq. (85)'s **mass exponent** is pinned separately by V3;
its **d₀/3 prefactor** and eq. (83)'s ε/τ are pinned by **nothing on the page**
and are transcription-only. That is measured, not asserted: breaking d₀/3 → d₀/2
in *every* route at once moves V1 only ×6.9 (which V1's own >100× rule prints as
BLIND), V2 from −1.0 × 10⁻¹³ to −2.8 × 10⁻¹³ and V3 from 5.7 × 10⁻¹⁴ to
8.2 × 10⁻¹⁴ — machine zero either side. The page carries that as a blind-spot
row, not a caveat.

Deliberate-break tables are printed for the first three. Flipping eq. (87)'s
off-diagonal sign moves the collocation deviation to 19.4, inverting eq. (85)'s
mass exponent to 1.0, and replacing the matrix by its diagonal to 2.9. Deleting
the Knudsen term does not move it — it makes [Bᵉ] **singular**, so there is no
answer to be wrong, which is the Background section's claim demonstrated rather
than asserted.

**Scope.** The scalar closure and the Bosanquet relation are `A4.4`, which
builds the binary 2×2 [Bᵉ] only as the yardstick for that closure; bulk
multicomponent diffusion with no wall is `A4.2`; uphill diffusion measured in
free space is `A4.9`. Not implemented here, though the review discusses all of
them: surface diffusion, pore-size distributions, viscous selectivity
(α′ᵢ ≠ 1), thermal diffusion, non-ideal mixtures.

**One open item, not blocking.** The review's Fig. 44 — the dusty gas model
against Remick & Geankoplis's (1974) He/Ne/Ar capillary measurements — has been
digitised, but the extraction has not been through the maintainer's visual
review, so it is parked in `queue_cases/A4.3/review/` (git-ignored) and no cell
on this page loads it. Adopting it would add one Results section and two
parameter-free experimental checks, and would take the data tier from 6 to 4. See
the `follow_up` block in `queue_cases/A4.3.yaml`.
