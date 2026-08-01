# A4.4 — Knudsen diffusion and the Bosanquet relation

What happens when the pore is small enough that molecule–wall collisions
compete with molecule–molecule collisions, and exactly where the Bosanquet
interpolation `1/D_eff = 1/D_K + 1/D_AB` is right.

**Source.** The origins are unreachable — Knudsen (1909) and Bosanquet's 1944
British wartime report BR-507 are not on disk and were **not consulted**.
Every equation was read from **Krishna & Wesselingh (1997)**, *The
Maxwell–Stefan approach to mass transfer*, Chem. Eng. Sci. 52(6) 861–911
([doi:10.1016/S0009-2509(96)00458-7](https://doi.org/10.1016/S0009-2509(96)00458-7)),
which is on disk, prints the dusty gas model in full (eqs. 82–87), derives the
Bosanquet formula as its eq. (110), *names* it, and states the three
conditions it needs. Equations 82, 83, 84, 85, 105, 107, 109 and 110 were read
off 600 dpi renders on 2026-08-01 — the PDF is an Acrobat Capture OCR whose
text layer renders eq. (86)'s prefactor as `c_t/RT` where the page prints
`1/RT`, and eq. (109)'s `nu_1` as `nu_{1i}`.

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page: eq. (85) and the transition curve; eq. (109)
  re-derived from eq. (82) with sympy; a binary dusty-gas solver in pymrm with
  three swappable flux closures; the isobaric counter-diffusion cell and the
  reacting pellet; the Bosanquet error mapped across Knudsen number.
- `data/krishna-wesselingh-1997-porous-media.csv` — the printed constants and
  stated results of the review's dusty-gas sections (tier 6; provenance in the
  sidecar). The review's pair diffusivities are loaded **cross-page** from the
  published `A4.2` dataset; nothing was re-transcribed.

**The result.** The error of the Bosanquet interpolation is
`D_Bos/D_1 - 1 = (rho - 1) x_1 / (1 + Kn)`, where `rho` is the flux ratio the
application imposes. It is therefore exact under Knudsen control and worst
under bulk control — the opposite of the instinct that an interpolation is
worst in the middle. Measured on a spherical pellet running `A -> 2B` **at a
Thiele modulus of 30**, Bosanquet overestimates the effectiveness factor by
27.0 % at Kn = 1e-3, 14.1 % at Kn = 1 and 3e-5 % at Kn = 1e6. The bulk-limit
figure is grid-converged: 26.97 % at the sweep's n = 600, 26.96 % at n = 4800.
It is also conditions-specific — the same error is 25.5 % at phi = 10 and
27.7 % at phi = 100, so phi is named wherever the number appears; what does not
depend on phi is the monotone decay with Kn, which is the result. In the
isobaric H₂/N₂ counter-diffusion cell Bosanquet underestimates the hydrogen
flux by 41.3 % at Kn = 1e-3 and 0.04 % at Kn = 1e3, and it gets the internal
pellet pressure wrong by a factor of two under bulk control.

**Validation** (tier 6 — no measurement anywhere on this page). Jackson's
pellet-centre pressure `p_0 = sqrt(nu_B) p`, which the review quotes on p. 892
as "a 40 % increase" for `nu_B = 2`, is reproduced by a solve that encodes it
nowhere: 41.42 % against `sqrt(2) - 1`, worst 1.1e-6 relative over four
stoichiometries and a spread of 3.6e-8 over 18 combinations of geometry,
Thiele modulus and grid. That leftover 4.2e-7 is a **finite-Kn model residue,
not solver accuracy**: it scales as 1/Kn (dev × Kn constant to 0.04 % from
Kn = 1e5 to 1e8 on a fixed grid). Graham's law (eq. 107) emerges from the
isobaric slab to 1.6e-12, and rises to 1.6e-1 when eq. (85)'s mass exponent is
inverted. The review's two scaling claims come out at slopes 0.993/0.996
against a stated 1 and 0.033/0.063 against a stated 0. Eq. (109) is re-derived
from eq. (82) exactly and reduces to eq. (110) exactly. The geometry index is
checked against the classical effectiveness factor of a slab, a cylinder and a
sphere (1.7e-6, against 4.9 % for the wrong index). The `sqrt(nu_B)` check is
broken on purpose eight ways: it catches a flipped mass exponent (100 %), a
dropped mass dependence (41 %) and a wrong stoichiometric coefficient (50 %),
and is **blind** to geometry, to a five-cell grid, to an unconverged solve and
to a sign flip in the `[B^e]` off-diagonals — all stated on the page, with the
sign flip instead caught by the closure comparison (1.0e0) and by Graham's law
(2.6e-2). Grid order 2.21; worst relative Newton residual 1.4e-10. Runtime
~9 s.

**Scope.** Binary and scalar (`S3`) by design. The n-component matrix form of
the dusty gas model, the viscous-flow term, and the review's Fig. 44
comparison against the Remick & Geankoplis capillary measurements are
catalogue entry `A4.3` and are not built here. Bulk multicomponent diffusion
with no pore is `A4.2`, built from the same review.
