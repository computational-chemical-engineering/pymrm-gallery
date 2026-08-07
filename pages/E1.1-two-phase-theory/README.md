# E1.1 — Two-phase theory of fluidisation

From Toomey, R. D. & Johnstone, H. F., *Gaseous Fluidization of Solid Particles*,
Chemical Engineering Progress **48**(5) 220–226, May 1952 (presented at the
Forty-third Annual Meeting, AIChE, Columbus, Ohio).

**There is no DOI, and none was invented.** The article predates DOI
registration and AIChE has not retro-registered volume 48; CrossRef resolves
nothing for the title, the author pair or the page range. The citation was
verified from the scan's own title page and running footers.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 5 s (`runtime_seconds: 6` declared) |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata |
| `data/…_table1_beads.csv` | Table 1 — the five glass bead sizes |
| `data/…_table2_pressure_losses.csv` | Table 2 — 41 runs |
| `data/…_table3_wall_velocities.csv` | Table 3 — 13 runs of wall particle velocity from 16-mm Fastax film |
| `data/…_printed_results.csv` | the scalars printed in the running text |

Each CSV has a `.meta.yaml` provenance sidecar.

## Read this first: whose particles the headline is about

The four intercept ratios this page turns on — `u_mf/u_e` = 0.85, 0.69, 0.83,
0.87 — are printed by Toomey and Johnstone, but they are read off **Lewis,
Gilliland and Bauer's** four particle sizes, *not* off the five Scotchlite bead
sizes tabulated in this paper. The authors say so in the sentence before them:

> *"Only the data of Lewis, Gilliland and Bauer (12) extend to low enough
> velocities to obtain an accurate indication of the point of intersection."*

and their Figure 3 separates the two sets in its own legend — three filled
curves, *"AUTHORS' DATA"*, at `D_p` = 0.0148, 0.00861, 0.00418 in., and four
open curves, *"DATA OF LEWIS, GILLILAND AND BAUER"*, at 0.0224, 0.0178, 0.0112,
0.00881 in. Four open curves, four printed ratios.

So this page tests **the postulate**, on numbers this paper printed and read. It
is **not** a measurement of `u_e` for this apparatus, and the notebook says so in
its title cell, in Results §4 and in the Reuse section as well as here.

The page runs the refit on the authors' own beads anyway — the one it declines to
rely on. Fitting equation (11) on the legible cells, restricted to runs sharing a
bed height, gives `u_mf/u_e` = **2.05, 0.93, 0.47** for beads 2, 3 and 4: a
factor 4.4 spread **straddling 1**. Their own apparatus cannot fix even the sign
of the effect. That is what makes Lewis, Gilliland and Bauer's ratios
load-bearing rather than incidental.

## Why this page exists next to `E2.1`

`E2.1` reproduces Kunii and Levenspiel's bubbling-bed model, in which
`δ = (u₀ − u_mf)/u_b` is **assumed**. That page is tier 6 and says so: nothing on
it is compared with a measurement, and structurally nothing on it can be.

Toomey and Johnstone printed the postulate *and* the data. This page is the test.

## The reprint route was tried and failed

Kunii and Levenspiel (1968) state the two-phase relation as their equations 2
and 3 but never name it "two-phase theory", never attribute it to Toomey and
Johnstone (their whole Literature Cited list on journal page 492 was read), and
never test it. Froment, De Wilde and Bischoff's third edition does not contain
the string "Toomey" at all and routes the bubble–emulsion interaction through
Davidson and Harrison instead. `AGENTS.md` records this as the canonical example
of a reprint that carries the equation but not the case. Only the original would
do, and this page is built from it.

## The postulate has two halves and the printed numbers separate them

Equation (3), journal page 223:

```
(u_f − u_mf) A = V_g
```

Four columns later, on page 224, the authors report where the intercept of the
equation (11) lines actually lands — for Lewis, Gilliland and Bauer's particles,
as the section above explains:

> *"For the four particle diameters, the ratio u_mf/u_e is 0.85, 0.69, 0.83,
> 0.87."*

- **The half that survives**: the dense-phase flow really is constant,
  independent of `u_f`. That is what the straight lines demonstrate.
- **The half that fails**: that constant is not `u_mf`. It is about 1.23 `u_mf`,
  so the visible bubble flow is `Y = 0.765` of the postulate's at `u_f = 2 u_mf`,
  and `0.646` (200 °F) or `0.756` (70 °F) at the lowest of Toomey and Johnstone's
  own runs that robustly bubbles.

Everything after 1952 kept the second half.

## The three results, in order of how hard they are to argue with

**1. The bed-weight identity, and the claim Toomey makes without quantifying it.**
In the printed Discussion, Max Leva asks whether ΔP_ke is an increase over what
one would calculate from the weight of the bed, and Toomey replies that the
pressure drop at incipient fluidisation *"has been both calculated and measured
and the results checked very closely"*. He never says how closely.
`ΔP_mf/L_mf` from Table 2 must equal `(1−ε)(ρ_s−ρ_g)g` from Table 1 — a
manometer over a ruler against a pycnometer and a tapped-bed voidage. Worst
single run **0.49 %** over 25 runs and four independently printed voidages — and
that residual is not the voidages' fault, as the page's printed rounding budget
shows: ε buys 0.079–0.086 %, ρ_s about 0.28 %, and Table 2's own last digits
0.31–0.85 %. The worst run sits essentially on its own Table-2 rounding floor
(0.4884 % against 0.3782 %), so the identity holds to the resolution the printed
data allow.

The identity is also the tool that decided four otherwise ambiguous digits of
Table 2 and recovered bead 3's illegible void fraction (0.41929 → printed 0.420).
Bead 3 is therefore **excluded** from the test above; using it there would be
circular, and an earlier draft that did so is called out on the page.

**2. A printed constant that is wrong.** Table 3's derived column `u_p/D_p^0.5`
recomputes to within 0.79 % for bead sizes 3 and 4 across nine rows. Bead size
2's four rows are *all* off by the same **+3.96 %** — which is to say the
**printed** column is **low** by 3.81 %, the direction a too-large `D_p` in the
divisor requires — with a peak-to-peak spread about that offset of 0.089 %,
forty-five times smaller than the offset. Reading noise cannot do that; a single
wrong constant can. The column was evaluated with `D_p ≈ 0.0160 in.` where
**three independent printings** say 0.0148 in.: Table 1's inch column, Table 1's
micron column, and the label on the authors' own curve in Figure 3. So the error
is in Table 3's arithmetic, not in the particle size. **Reported, not repaired.**

Bead 4 carries the same signature six times weaker — offset +0.690 %, spread
0.180 %, implied `D_p` = 0.00424 in. against a printed 0.00418 in. — and the page
reports it as a second, weaker instance rather than leaving it silent. It is a
continuum, not a clean two-way split.

**3. An exponent that cannot be reconciled.** The paper says on p. 225 that *"the
velocity associated with the continuous phase is proportional to the square root
of the particle diameter"*; it says on p. 224 that `u_mf/u_e` shows no dependence
on particle diameter; and the five printed incipient-fluidisation Reynolds
numbers give `u_mf ~ D_p^1.51 ± 0.13` (and `D_p^1.47` from the two end points
with no fit at all). The gap is a full unit, 7.9 standard errors. All three
cannot be true. (Equation (13) itself is printed for `u_t`, the terminal
velocity, not `u_e` — the p. 226 Notation separates the subscripts — so the page
quotes the prose sentence, which is where the identification is actually made.)

**4. A cross-check on the column nothing else can see.** Table 3's `ΔP_ke/L` is
used in exactly one place, as Figure 5's ordinate, and `fig5_null_gain` is a poor
guard on it (a +30 % single-cell error drifts that metric 0.07–12 %, mostly
inside CI's 5 % tolerance). Dividing `ΔP_ke/L` by `ΔP_mf/L_mf` gives an
independent estimate of Table 2's `ΔP_ke/ΔP_mf` column: over the ten runs where
both are legible and the row match is unambiguous, median **5.25 %**, worst
24.5 %, RMS 9.41 %. That is `dPke_cross_check_rms`, it has its own break rows,
and it is the *partial recovery* the dataset sidecar refers to — computed, not
asserted.

## What pymrm does here, and what it does not

Results 1, 2 and 4 are arithmetic on printed tables. There is no PDE in a
bed-weight identity and none was manufactured.

The pymrm solve answers the one question the authors raise and drop: they say the
gas must expand as the pressure changes through the bed, and that this *"would
cause the data in Figure 3 to show a slight curvature near the abscissa"*. Their
deepest bed drops 44.3 in. H₂O against one atmosphere and the reference run used
here drops 40.7 in. over 24.6 in. of bed. The solve — a nonlinear, piecewise,
two-point problem with an unknown domain length, `construct_convflux_upwind`
+ `construct_div` (`nu=0`) + `newton`, with `L_f` root-found by `brentq` — gives:

* the excess velocity `u_f − u_mf` varies by **25.8 %** across the bed at the
  paper's stated 200 °F and by **20.5 %** at 70 °F; both exceed the shortfall the
  intercepts report on the additive reading, `1 − u_mf/u_e` = 19.0 % (on the
  multiplicative reading, `u_e/u_mf − 1` = 23.5 %, only the 200 °F figure does —
  the page states which it means);
* equation (1) applied with a single `u_f` over-predicts the expansion by 4.2 %;
* a band `1.235 ≤ u_f/u_mf ≤ 1.358`, exactly as wide as the bed's pressure ratio,
  in which the top bubbles and the bottom does not — straddling the threshold the
  Figure-3 intercepts are extrapolated to.

Reported metrics are the `n = 1600` values, not the Richardson limit; the
Richardson extrapolation is printed as evidence that the grid ladder is clean.

## No figure was digitised, because none needed to be

Figure 5 plots `ΔP_ke/L` against `u_p/D_p^0.5` — **both axes are columns of
Table 3**, so the thirteen points are printed. The page fits them, lets the
exponent on `D_p` go free, and reports the confidence interval: `[0.47, 0.77]`,
with the `n = 0` null baseline 55 % worse. The square root does real work; *the
square root specifically* is not identifiable from thirteen points over three
bead sizes.

## Reading a scan with no text layer at all

`pdftotext` returns **one byte per page**. Everything was read from 400 dpi
`pdftoppm` renders of a microfilm scan stored as heavily-compressed JPEG tiles.
The transferable method is the **arithmetic pin**: four digits of Table 2 were
decided not by looking harder at the pixels but by the constraint the table
imposes on itself. Where no such constraint exists — the whole `ΔP_ke/ΔP_mf`
column — cells are **empty**, not guessed. Seventeen of forty-one are empty, and
that is why this page takes the four intercept ratios the authors print instead
of refitting equation (11).
