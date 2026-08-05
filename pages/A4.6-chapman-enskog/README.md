# A4.6 — Chapman–Enskog binary diffusion

A binary diffusivity that is *derived* from the intermolecular potential rather
than fitted to diffusion data — and the honest question that makes it worth
building: how much of its accuracy is the kinetic theory, and how much is the
two-parameter fit underneath it.

**Source.** Chapman, S. and Cowling, T. G. (1970), *The Mathematical Theory of
Non-Uniform Gases: An Account of the Kinetic Theory of Viscosity, Thermal
Conduction and Diffusion in Gases*, 3rd edn, Cambridge University Press
(Cambridge Mathematical Library reissue 1990, reprinted 1993). This is the
**origin**, not a reprint — the primary account of the result the case is named
for, by one of the two people it is named after — so there is no
`reference_read_from`. The catalogue carried no citation at all for `A4.6` until
the monograph reached disk on 2026-08-05.

The PDF is a complete 448-page CCITT-G4 bilevel scan at **300 ppi native** with
an ABBYY FineReader text layer. The prose layer is good; the numbers are not —
the 1970 typesetting uses a raised mid-dot decimal separator that the OCR drops,
and leading zeros come back as the letter "o". Every numeral on the page was read
off a **cropped 300 dpi render**, twice. Rendering at 600 dpi would interpolate a
1-bit image.

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page.
- `data/chapman-cowling-1970-table22.csv` — **50 measured D₁₂ at S.T.P.** with
  the book's own literature references (**tier 2**).
- `data/chapman-cowling-1970-table6.csv` — the 16 rows of collision integrals.
- `data/chapman-cowling-1970-table11.csv` — molecular weights, viscosities and
  viscosity diameters.
- `data/chapman-cowling-1970-table17.csv` — the Lennard-Jones force constants,
  in both of the book's two published versions.

**No figure is digitised anywhere on this page, and no dataset is borrowed from
another page.**

## The result

Over the 39 pairs of Table 22 for which Table 17 supplies 12,6 constants for both
gases, the first approximation (eq. 14.2, 4) predicts the **measured** D₁₂ with a
**mean absolute deviation of 3.65 %**, a bias of **−2.74 %** and a worst case of
**−10.06 %** (N₂–CO₂); 29 of 39 inside 5 %, 38 of 39 inside 10 %. That is what
§14.4's "in general agree fairly well" is worth, and the book prints no such
number anywhere.

The comparison is **out of sample**: the force constants were fitted to pure-gas
*viscosity* (p. 235, explicitly), and no diffusion measurement enters the
prediction at any point. Table 23's constants, which *are* derived from D₁₂ data,
are deliberately not used.

**But most of that accuracy is the parameter fit, not the theory.** The book
prints a second set of force constants for the same molecules, obtained mainly
from virial coefficients, and warns that other authors' values "often differ
appreciably". Running the *identical* calculation on them gives **10.23 %** where
the viscosity constants give **4.13 %** on the same 27 pairs. Nothing about the
kinetic theory changed between those two numbers. So the page reports the
headline as what it is: *the 12,6 model transfers a two-parameter fit from
pure-gas viscosity to binary diffusion with a 3.6 % penalty, and from equilibrium
virial data with a 10 % penalty.*

The **null model** shows the collision integrals are nonetheless doing real work:
the same expression with 𝒲⁽¹⁾(1) deleted and Table 11's viscosity diameters gives
**17.24 %**, biased −17.24 % — every pair under-predicted, 4.7× worse. That bias
is predicted independently by the book's own statement that the viscosity
half-sum exceeds the diffusion σ₁₂ "by about 10 per cent", reproduced here at
10.46 % over all 50 rows.

## Validation

**The headline is computed a second, independent way.** Ω⁽¹'¹⁾\* and Ω⁽²'²⁾\* are
recomputed from the 12,6 potential by quadrature of the classical deflection
angle — no table consulted. It reproduces Table 6 to **0.18 %** (Ω⁽¹'¹⁾\*) and 0.30 % (Ω⁽²'²⁾\*) at worst over all
sixteen rows — both worst cases on the kT/ε = 100 row, where the quadrature is
itself converged to five digits under a 2.8× refinement, so that residual is a
difference between two numerical integrations rather than an error in this one.
Every other row agrees to better than 0.1 %. Substituting it for Table 6
entirely changes every predicted D₁₂ by at most **0.066 %** and moves the
headline from 3.647 % to 3.655 %. A break table
could not have established this: it shows the answer is *sensitive* to Table 6,
never that Table 6 was read *correctly*.

**One reading had to be reconstructed, and it corrects the source note.** Table
6's third column header carries a fraction whose glyph is destroyed in the scan.
It is **½, not ⅓** — established twice from printed material: eq. (9.8, 7) with
the rigid-sphere normalisation of eq. (10.2, 1) (2.9e-4 against ½, 0.50 against
⅓), and the independent quadrature (0.30 % against ½, 34 % against ⅓ — a factor
of 112). It changes
no diffusivity, but the prior source note for this book has it wrong.

**Four printed-table round trips.** Table 22's σ₁₂ from its D₁₂ through
eq. (14.2, 1) — worst 0.00509 on a column printed to 0.005, *and this pins the
state*: the mean signed residual is +0.00016 at 1 atm against a one-signed
+0.02257 at 1 bar (137× larger), so "S.T.P." is established from the book's own arithmetic
rather than assumed. Table 22's half-sum from **Table 11's** σ (worst 0.0075;
from Table 17's σ instead it is wrong by 0.740). Table 11's σ from its μ and M
through eq. (12.1, 6) (8.7e-4, rising to 8.8e-3 without the 1.016 correction).
Table 17's rₘ against 2^(1/6)σ (0.0052). Plus Table 6's A column against
eq. (9.8, 7) and its C column against eq. (14.4, 1).

**Four printed numerical results** of the second approximation, all recovered to
1e-12 or better: eq. (14.3, 2) at five mass ratios; Δ₁ → 1/13 in the Lorentz
limit; 1.083333 against the printed 1.083; 8.3333 % against "8⅓ per cent"; and
Kihara's B = ¾ giving 11.1111 % against "1/9, or 11·1 per cent". Six transcribed
equations have to be right at once for any of them to come out.

**Two "wrong" choices score better than the published ones**, and the page says
so: a geometric-mean σ₁₂ gives 3.26 % and a 1 bar reference 3.10 %. Neither is
adopted — eq. (14.2, 1) prints the arithmetic mean, and the σ₁₂ round trip pins
1 atm. The headline therefore sits on a shallow optimum and is quoted to one
decimal at most.

## What pymrm adds

§14.32 states that "the variation of D₁₂ with the proportions of a mixture is
ignored" when the data are discussed — but the apparatus behind many Table 22
entries is a closed tube whose composition runs from 0 to 1 along its own length.
Solving that cell in pymrm with D₁₂ = [D₁₂]₂(x₁) and inverting the constant-D
closed form for the apparent D, exactly as an experimenter would, gives
**+2.53 %** above [D₁₂]₁ for He–Xe, +1.80 % for H₂–N₂ and **+0.31 % for the
equal-mass control N₂–CO**. It is *not* the equimolar second approximation
(+2.68 / +1.89 / +0.31 %): the cell falls short by 0.15, 0.09 and 0.00 percentage
points, the shortfall growing with mass disparity and vanishing for the control.
A transient cell reports a weighted average over a composition history, not a
value at a stated composition. The page states plainly that this is a model of
the apparatus, not the apparatus, and is not a correction to be applied to
Table 22.

## Scope

Binary, dilute gases at S.T.P. The Fuller–Schettler–Giddings **empirical**
correlation for the same quantity is `A4.5`; the page argues that the joint
comparison — one measured-D₁₂ axis with a first-principles prediction and an
empirical correlation overlaid — should be built **there**, since `A4.5` must
transcribe its own source anyway and can load this page's Table 22 cross-page.
Multicomponent diffusion is `A4.2`/`A4.9`, pore transport `A4.3`/`A4.4`.

**This monograph does not source `A4.1`.** Its only Wilke is Buddenberg & Wilke
(1949) on the *viscosity* of a mixture, cited once in §12.43; the 1950 diffusion
mixture rule is not stated, named or attributed anywhere. Straight `E1.1`
failure.
