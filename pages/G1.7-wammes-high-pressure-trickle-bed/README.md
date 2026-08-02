# G1.7 — Wammes–Westerterp high-pressure trickle-bed hydrodynamics

Wammes, W. J. A. and Westerterp, K. R., *Hydrodynamics in a pressurized
cocurrent gas–liquid trickle-bed reactor*, Chemical Engineering & Technology
**14**(6) 406–413 (1991), doi:10.1002/ceat.270140608 (open access via the
University of Twente repository). Note: the gallery catalogue originally cited
the Chem. Eng. Sci. 46 (1991) 409–417 companion; the CET paper is the one that
derives the high-pressure correlations and is the one built here.

The page closes the paper's two mutually dependent correlations — pressure
gradient (Eq. 8, containing the hold-up) and dynamic hold-up (Eq. 9, containing
the pressure gradient) — as a coupled system, with nothing fitted, against 16
digitised hold-up points (Fig. 7) and 20 digitised pressure-gradient points
(Fig. 6). It also solves the paper's plug-flow CO₂-absorption balance (the
model behind every interfacial-area measurement) with pymrm operators and
verifies it against the paper's own closed form, Eq. (3).

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `data/` — digitised points and transcribed parameters, each with a
  provenance sidecar. **Digitisation reviewed and confirmed** by the maintainer
  on 2026-08-02 against the numbered overlays in `../review/` (git-ignored,
  they contain the page image): every cross on a real marker, one per marker,
  none missed, series correctly told apart by shape, and the four reconstructed
  centres in the touching pairs (fig7 at v_l = 2 and 4 mm/s; fig6 at
  ρ_g ≈ 0.45 and 1.6 kg/m³) checked individually. No CSV row changed.
- Headline numbers: hold-up mean |dev| 2.9 % (paper's own error: 8 %);
  pressure gradient 10.4 % (paper: 12 %). Deviation convention:
  (model − measured)/measured.
