# C2.10 — o-Xylene to phthalic anhydride: the classic hot spot

Froment's 1967 demonstration that effective kinetics inside a two-dimensional
pseudo-homogeneous model reproduces fixed-bed hot spots — the paper every
reactor textbook cites for it — rebuilt from the printed constants: both
models, every stated number, and the five-degree disagreement between them
measured continuously.

- **Structures:** `S2` (1-D plug flow), `S6` (2-D radial dispersion, marched)
- **Reference:** Froment (1967), Ind. Eng. Chem. 59(2) 18–27,
  doi:10.1021/ie50686a006
- **Covers:** `D3.4` (the multitubular-reactor half of the same paper; the
  scope argument is recorded in `queue_cases/C2.10.yaml` and
  `queue_cases/D3.4.yaml`)
- **Runtime:** ~2 min

## Agreement

Five stated model results reproduce to **4.0 % mean, worst 8.6 %** (the
363 °C hot-spot rise, a quantity that doubles every ~2 °C of inlet
temperature there). The 2-D hot spot at 357 °C comes out 29.53 °C against the
paper's "about 30". Both runaway limits land inside the paper's
integer-probing brackets — 360.01 °C in (357, 360] for the 2-D model, 363.93 °C
in (363, 365] for the 1-D — and the continuous gap between the models is
3.9 °C against the paper's "within five degrees". Four printed transcription
identities (Pe_hR, and U for three λ_R/α_w combinations) close to 0.29 % or
better. **Nothing is fitted anywhere.**

## Provenance: tier 6, reproduction not validation

The paper contains **no experimental data**. Its kinetics are asserted as
"fairly representative" of o-xylene oxidation on V₂O₅ — never fitted to
measurements in this document — and every comparison target on this page is
Froment's own 1967 computed result. Calderbank (the catalogue row's second
name) is not on disk. Figures 7–14 are used qualitatively only; nothing was
digitised.

## Three typographic traps, all resolved from the paper's own numbers

1. **"G = 4.684 kg./sq. meter hr." is 4684** — European thousands separator.
   Proof: p. 19 prints 4684 outright, and the printed Pe_hR = 5.25 equals
   G·c_p·d_p/λ_R only with 4684.
2. **"ΔH₃ = −1.090 kcal./gram mole" is −1090** — same separator; −1.090 would
   make total combustion 282× *less* exothermic than partial oxidation.
3. **"44 gram moles/cu. meter" is 44 g/Nm³** — Figure 11's own axis label;
   the arithmetic to N_A0 = 0.00924 closes for g/Nm³ and is off by 100× for
   mol/m³.

The break table runs the first two taken literally.

## Cross-page

Shared constants are reconciled against the independently transcribed
parameter CSV of `D2.2` (Van Welsenaere & Froment 1970 — the same reactor,
same k₁ three years later): eight of nine rows identical to <0.1 %. D2.2's
U = 82.7 is *inherited from this 1967 paper*, so that row is a shared source,
not a second witness; the independent check on U is the printed λ_R/α_w route.

## Rebuilding the page

```bash
python build_page.py                    # regenerate index.ipynb from the builder
python ../../../scripts/run_pages.py    # execute it (from a pages/ checkout: ../../scripts)
python ../../../scripts/check_agreement.py
```

No PDF is needed to rebuild: both datasets are committed CSVs.
