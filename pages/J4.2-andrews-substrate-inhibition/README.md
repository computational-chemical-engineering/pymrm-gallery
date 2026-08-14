# J4.2 — Andrews' inhibition function

**Two claims in one sentence of a 1968 Summary, and what happens when you test
them instead of illustrating them.**

Andrews puts Haldane's enzyme-inhibition function into a batch and a chemostat
balance and writes, on book p. 707:

> *"Simulation studies show that the primary result of inhibition by substrate
> in a batch culture is an increase in the lag time whereas in continuous
> culture inhibition by substrate may result in process instability."*

He also assumes, on p. 711, that there is **no lag phase**. So the first claim is
about a lag his kinetics manufacture, and it can be measured against a
no-inhibition control. The second is a statement about how many steady states a
chemostat has and which of them are stable, and it is settled by root-finding and
a 2×2 Jacobian.

**The paper prints no table.** Every "table" hit in the extracted text layer is
inside the word *"unstable"* — two of them, and that is all. What it prints is
ten figures, and **every parameter set is typeset inside a legend**, which is why
this case is not figure-gated. Reading a typeset legend is transcription.
**Nothing here is digitised**, and the page says plainly that it therefore does
not establish empirical adequacy against Andrews' plotted curves. Two things
read *inside* figures are declared where they are used: Fig. 7's typeset
annotation, which is text and settles a sentence of his prose, and the fact that
both of Fig. 9's `S_1` curves start at the origin, which corroborates an initial
condition he never prints — no axis was calibrated, no point extracted, and no
number on the page comes from either curve.

## What it finds

**The batch claim survives, quantitatively, and it could have failed.** Going
from no inhibition to Andrews' `K_i = 2.0` g/l multiplies the time to reach
5 g/l of biomass by **5.256387**, from 6.949160 h to 36.527476 h. Of that
29.578315 h increase, **94.40 % is lag**, 5.65 % is a slower exponential phase
and −0.05 % is the tail. His word *"primary"* is right, and the share is monotone
in inhibition strength — 87.78 % at `K_i = 10`, 91.25 % at 5, 94.40 % at 2. Under
Monod the lag is **identically zero** by this definition, because the trajectory
starts at its steepest point: inhibition creates a lag where there was none.

**Two things about that split, and the page states both rather than leaving them
to be misread.** First, the peak specific growth rate falls only **19.4344 %**,
0.997009 → 0.803246 h⁻¹ — and that is not the effect of inhibition on growth.
The peak is reached at 36.492522 h of a 36.527476 h run, **99.90431 % of the way
through**; the trajectory is above half its peak rate for 0.638869 h,
**1.7490088 % of its duration**; and the rate the culture actually realises,
`ln(X_f/X_i)/t_f`, falls from 0.994042 to 0.189111 h⁻¹, **80.98 %**. During the
27.922318 h the split calls lag, biomass rises **114.67072**-fold, 0.005 →
0.573354 g/l. (Those three are written to five or more decimals deliberately:
unlike the 19.4344 % and the 80.98 %, they are **not** `agreement.json` metrics,
so CI never compares them, and at two decimals they fell outside the sweep
below as well. At these digits the sweep reads them.) Second, **the share is a partition only near Andrews' own target**: `λ` does
not depend on `X_f` and `Δt_f` does, so at `X_f = 1` g/l the same decomposition
returns **109.4003 %** and at 0.5 g/l 123.8027 %. The 94.40 % is quoted at his
5 g/l target and is defined there.

**Both of his stated numerical results come back out of a closed form.** Constant
yield collapses eqs. (4)–(5) to one autonomous ODE with a rational integrand, so
partial fractions give the batch time exactly:

```
t = (1/mu_hat) [ (S_f - S_i)/K_i + (K_s/a) ln(S_i/S_f)
                 + (1 + a/K_i + K_s/a) ln(X_f/X_i) ],    a = S_i + X_i/Y
```

Both sympy residuals — the decomposition and the antiderivative — are exactly
zero. It gives **36.527476 h against Andrews' 36** (+1.4652 %) and **12.371597 h
against his 12.5** (−1.0272 %), and it shares no code at all with the solver it
is checked against (6.7×10⁻¹² relative to LSODA).

**The second of those depends on reading one sentence, and Andrews' own figure
reads it for us.** *"The substrate concentration being increased in 2.0 g/l
steps after reduction in each case to 0.02 g/l"* supports three literal
readings — reset to 2.0, increment by 2.0, or a 2.0/4.0/6.0 ladder — giving
12.366419 h, 12.371597 h and 13.794022 h: a raw spread of **11.5442 %**. The
typeset annotation inside Fig. 7 itself, read at digit scale on a 300 ppi crop
of the rotated landscape page, is *"REPEATED ORGANISM SEPARATION & DILUTION.
`S_i` ADDED IN INCREMENTS OF 2.0 GM/L"*. That is the increment reading, so it is
the one reported; the reset reading is priced at 0.042 % away and the ladder
reading, at +10.3522 %, is the one the annotation excludes. **The page does not
call the sentence ambiguous, because his figure is not**, and the annotation is
transcribed into the CSV as `fig7_annotation`.

**The continuous claim is proved rather than illustrated.** At `θ = 3` h,
`S_0 = 5` g/l there are **three** steady states:

| | `S_1` (g/l) | `X_1` (g/l) | eigenvalues (h⁻¹) | |
|---|---|---|---|---|
| lower | 0.0150567 | 2.492472 | −73.019142, −1/3 | stable node |
| upper | 3.984943 | 0.507528 | **+0.0561790**, −1/3 | **saddle** — eq. (9)'s root |
| washout | 5 | 0 | **−0.0481080**, −1/3 | **stable node** |

The saddle's sign is not numerical luck. With `μ(S*) = 1/θ`,
`det J = μ'(S*) X_1/(Yθ)`, so `det J` carries the sign of `μ'` **exactly**: the
lower root is on the rising limb and stable, the upper is on the falling limb and
a saddle, whatever the constants. And washout is stable, so the reactor is
bistable — which is what *"may result in process instability"* means.

**Inhibition barely moves the steady state and completely changes its basin.**
Andrews says the first half himself (*"the inhibition function reduces to the
Monod function"*, p. 709) and understates it: at the same `mu_hat`, `K_s`, `θ`
and `S_0`, dropping the inhibition term moves the operating substrate
concentration by **0.3778 %** and the operating biomass by **0.0011 %** — and
takes the washout eigenvalue from **−0.0481080** to **+0.660702** h⁻¹, from
stable to unstable. Never conclude from *"the inhibition term is negligible at
the operating point"* that it is negligible.

**There are two washout conditions and eq. (11) is only one of them.** Andrews'
`θ_w = 1/mu_hat_m = 1.244949` h is the **fold**, where the two roots merge, and
it contains no `S_0`. The other — washout is stable whenever `μ(S_0) < 1/θ`, i.e.
`θ < 1/μ(S_0) = 3.506000` h — depends on the feed, and he does not print it. The
**exact bistable window** at `S_0 = 5` g/l is therefore

```
theta in (1.244949, 3.506000) h        and Andrews operates at 3 h, inside it,
```

which is why every one of his Fig. 8–10 experiments has something to show. The
two conditions also unify him with J4.1: under Monod `μ` is monotone, the fold
never happens, and only the transcritical condition survives — as
`D_c = μ_m S_f/(K_S + S_f)`, exactly the formula J4.1 loads from Rawlings &
Ekerdt. **One criterion, `D = μ(S)`, read at two different points.**

**All six of his stated continuous-culture outcomes reproduce, and the page
root-finds the thresholds he does not print.**

| Andrews, book p. 719 | outcome | this page |
|---|---|---|
| `S_0` 5 → 5.2 from the saddle | washout | washout ✔ |
| `S_0` 5 → 4.8 from the saddle | recovers | recovers ✔ |
| startup at `X_i = 0.10` | fails | fails ✔ |
| startup at `X_i = 0.50` | recovers | recovers ✔ |
| step `S_0` 5 → 20 | fails | fails ✔ |
| ramp `S_0` 5 → 20 over 1 hr | recovers | recovers ✔ |

The critical inoculum is **0.1842466 g/l** — his 0.10 is 45.7 % below it and his
0.50 is 171.4 % above; over the whole sampled range of the `S_1(0)` Andrews never
prints, `[0, 0.2]` g/l, it varies **5.2120 %** (0.4071 % over `[0, 0.05]` alone,
which is the flattering window and not the one reported), and `S_1(0) = 0` is
**corroborated by his own Fig. 9**, where both `S_1` curves emanate from the
origin. The critical step is **`S_0` = 17.28760 g/l**, and the
critical ramp duration is **0.8720668 h**, so **his 1 hr ramp clears it by only
14.67 %**. The Fig. 8 threshold, by contrast, is exactly `S_0 = 5.0` **by
construction** and is not a measurement; what those two excursions test is the
direction, and the page says so rather than reporting 5.0 as a result.

**Two of those thresholds are computed twice, by routes sharing nothing but the
right-hand side.** The critical inoculum by bisection on the outcome of the
initial-value problem *and* by backward-integrating the saddle's stable manifold
down to `S = 0`: agreement 2.4×10⁻¹². The critical step by bisection *and* by
asking which side of the new system's separatrix the old operating point sits
on: 8.1×10⁻¹⁴.

## What the page adds that Andrews did not have

**Seed the feed and the instability disappears.** He defines `X_0` in eq. (7) and
never gives it a value. A feed biomass of **0.0213935 g/l — 0.8583 % of the
operating concentration** — destroys the fold, because washout stops being a
steady state at all. **That number is where the page caught its own grid-limited
calculation**: bisecting on *how many* steady states there are counts sign
changes on a scan and can only see the two merging roots while they are more than
a grid cell apart, so it reports the fold too early. Solving for the **double
root** (`G = G' = 0`, which eliminates `X_0` and leaves
`μ(θμ − 1) = μ'(S_0 − S)`) has no grid in it. The double-root value is the one
reported, and refining the counting route drives it onto that value across six
orders of magnitude in three scans — 4.9×10⁻³, 5.3×10⁻⁶, 1.2×10⁻⁹.

**And one warning this page had to measure rather than inherit.** The standing
advice in this repository is that an outlet read off the last cell centre is
`O(h)` against a second-order `compute_boundary_values` read. Refining **both**
reads on the same four grids against the same batch reference says otherwise
here: orders **1.0011** and **1.0013** — both first order — and the last-cell
read is the *closer* of the two at every grid, by a ratio that hardly moves
under refinement (1.3215 at `n` = 100, **1.3226** at `n` = 800). That is
what a zero-gradient outflow condition does, and J4.1 measured the same thing on
its own problem with a different constant. The call is kept, and the reason
stated is the right one: `compute_boundary_values` returns the value the flux
operator transports, so a balance written on it closes — consistency, not
accuracy. For the same reason **no plug-flow outlet is reported off a single
grid**: the `τ = 6` h outlet is the Richardson value of `ncell` 1200 and 2400,
0.6076462 g/l, with the raw 1200-cell value printed beside it and its 0.1637 %
distance from the batch reference stated.

**Plug flow has none of this.** The same kinetics in a pymrm plug-flow fermenter
give an outlet biomass that increases strictly with residence time for any
positive seed — no fold, no saddle, no washout — because plug flow *is* the batch
trajectory in space and `dX/dt = μX > 0`. That is the infinite-stage limit of
Andrews' own *"multistage operation with separate substrate supply to each
reactor"*, and it confirms his prediction **in that limit only**: a finite cascade
is not solved here, and his other two proposals are recycle problems and are
untouched.

## Seven printed defects, reported and none repaired

- **eq. (9) is printed with a bare `+`, not `±`** — read at digit scale. Its own
  paragraph says the quadratic has two roots. The root it returns is the
  **higher** one, which this same paper calls unstable. The page root-finds both
  and checks eq. (9) against the upper one; it does not correct it.
- **eq. (10) is printed `X_i = Y(S_0 − S_1)`** where `X_1` is meant — `X_i` is
  this paper's symbol for the *initial inoculum* (p. 713, Fig. 9's key). The text
  layer of this scan renders `X_1` as `XI`, so only a render could show it.
- **Figs. 4 and 7 print `Y = 0.5 GM/L`** for a yield the text defines as *"mass
  organisms produced/mass substrate utilized"*. What decides it is that printed
  definition plus the dimensions of eq. (5), not the five-to-two count across the
  legends. The numeral is not in question, only the unit.
- **`"Recherces"`** [sic] for *Recherches* in reference 1; **`"step chances"`**
  [sic] for *step changes* on p. 719; **`"innoculum"`** throughout; **`"a
  inhibition function"`** [sic] in the Summary.

## An approximation that is *not* filed as a defect

Andrews defines `K_s` and `K_i` as the low and high substrate concentrations at
which `μ = mu_hat/2`. Solving eq. (1) for that gives
`S = K_i(1 ∓ sqrt(1 − 4K_s/K_i))/2`, so both definitions are exact only as
`K_s/K_i → 0` and **neither has a solution unless `K_i ≥ 4K_s`** — which is
exactly the condition under which his own eq. (2) gives `mu_hat_m ≥ mu_hat/2`. At
his own `K_i = 2.0` g/l the exact roots are **0.0304640** and **1.9695360** g/l,
**+1.5468 %** and **−1.5232 %** off the values he names them for — and that is
the mildest of the three `K_i` he prints: at the `K_i = 1.0` g/l of his Fig. 4
it is +3.1947 % and −3.0958 %, and at the `K_i = 0.50` g/l of his Fig. 1
**+6.8502 %** and **−6.4110 %**. The page prints all three. This is an
approximation, not a misprint, and the nearest thing to a qualification is on the
same page — *"even when `K_s` and `K_i` are well separated..."*. The page gives
the number and the exact condition and **does not claim he was unaware**.

## What J4.1 left here, and what goes back

`pages/J4.1-monod/` was published immediately before this page and left substrate
inhibition here on purpose: it transcribed **Froment eq. (1.5.2-4)** and
**Rawlings & Ekerdt's `mu_m S/(K_s + S + K_1 S^2)`** into its CSV, flagged both
*"OUT OF SCOPE — belongs to J4.2"*, and never evaluated either. This page loads
that CSV — **retyping nothing** — and proves in sympy that both are Andrews'
eq. (1), with Froment's `r_m` = Andrews' `mu_hat` and Rawlings' `K_1` = `1/K_i`.

And it finds a trap in the crossing. J4.1's `froment_rm_definition` row records
Froment defining `r_m` as *"the maximum specific rate of biomass growth"* — true
for his Monod eq. (1.5.2-1), **false** for his inhibition eq. (1.5.2-4), where
the maximum is `r_m/(1 + 2 sqrt(K_S/K_i))`. At Andrews' constants, carrying that
definition across **overstates the maximum by 24.4949 %**. J4.1's
`froment_asymptote` row already notes that the inhibition form *"exhibits a
maximum"* and that Froment prints no formula for it; Andrews' eq. (2) is that
formula.

Going the other way: J4.1 had to record that Rawlings & Ekerdt cite their five
growth laws to two textbooks, so *"nothing on this page is attributed to
Blackman, Tessier, Moser or Contois personally"*. **Andrews' reference list
carries the primary citations for three of those four** — Moser (1958), Teissier
(1936, spelled `Teissier` where Rawlings spells the law `Tessier`; neither
repaired) and Contois (1959). They are transcribed here for whoever builds those
cases. **None is on disk, none was consulted, and nothing is attributed to any of
them.**

## Sources

**Origin, on disk and read in full** — Andrews, J. F., *A Mathematical Model for
the Continuous Culture of Microorganisms Utilizing Inhibitory Substrates*,
**Biotechnol. Bioeng. 10**(6), 707-723 (1968), doi `10.1002/bit.260100602`.
Identity confirmed from its own first page. `pdfimages -list` reports every page
as CCITT-G4 bilevel at **300 ppi native**, so 300 ppi is native and higher would
be interpolation; all seventeen pages were read there and every numeral re-read
on a digit-scale crop. Three pages are landscape in the file (PDF 8, 12, 15 —
Figs. 4, 7, 9) and `pdftoppm` returns them rotated; those legends were rotated
back before reading. **The text layer was used only for searching**: it turns
`X_1` into `XI`, which is precisely the failure that would have hidden eq. (10)'s
defect.

**Cited by Andrews, not on disk, not consulted** — Haldane (1930), the origin of
the inhibition function itself; Dixon & Webb (1964), his source for fitting it,
which is one reason this page fits nothing; Monod (1942); Moser (1958); Teissier
(1936); Contois (1959); Koga & Humphrey (1967); Brennan's PACTOLUS report (1964).
Nothing is attributed to any of them beyond what Andrews prints.

**Borrowed** — `pages/J4.1-monod/data/printed-growth-laws.csv`, six rows.
Froment, Rawlings & Ekerdt and Levenspiel were **not opened for this page**;
everything from them arrives through J4.1's transcription.

## Files

| file | what it is |
|---|---|
| `index.ipynb` | the page; nine sections, executed clean, 56 s |
| `build_page.py` | regenerates `index.ipynb` |
| `data/andrews1968-printed-model.csv` | 51 rows — every equation, figure-legend parameter set, in-figure annotation and claim used, verbatim, 40 flagged |
| `agreement.json` | 63 metrics; 5 sit below CI's `ABS_FLOOR` and are named with above-floor companions |

## Validation summary

- **the paper's own identities, all exactly zero in sympy**: eq. (1) ≡ eq. (6);
  eq. (2) and eq. (3) re-derived from eq. (1); the batch-time partial fraction
  and its antiderivative; Andrews ≡ Froment ≡ Rawlings. Plus eq. (11) −
  1/eq. (2) = `0.0` in double precision. **All printed, none reported** — a
  symbolic zero is a proof, not a measurement;
- **two stated numerical results reproduced** to +1.4652 % and −1.0272 %, both
  from a closed form that shares no code with any solver, and **labelled
  reproduction, not validation**: they are Andrews' own PACTOLUS output, read by
  him off his own plot;
- **six stated outcomes reproduced** and the three thresholds he omits
  root-found on **event-terminated** integrations, never swept;
- **six quantities computed twice, of which four are genuinely independent** —
  the batch time (6.7×10⁻¹²), the critical inoculum (2.4×10⁻¹²), the critical
  step (8.1×10⁻¹⁴) and the `X_0` fold (1.2×10⁻⁹) — while the other two, eq. (9)
  against a Brent root (1.1×10⁻¹⁶) and the transcritical `θ` by its two routes
  (exactly `0.0`), are the same algebra written twice and are labelled
  **transcription checks**, not second computations; the break table
  mis-transcribes one term of each to show what they can catch;
- **both error-carrying axes refined**: marcher time step, observed order 1.0021,
  Richardson 338× better than the finest step; plug-flow grid, observed order
  1.0011, Richardson 637×;
- **both branches of the lag split and the switch between them** exercised, at
  `K_i = S_i²/K_s = 3333.3333` exactly;
- **twenty-four defect injections, with the coverage map generated from their
  measured moves** rather than written by hand: each row returns a dictionary
  keyed by metric name holding those metrics recomputed under the defect, a row
  covers a metric only if it moves it by more than 10⁻⁶ relative, and the
  notebook fails to execute if any metric is uncovered or any row moves nothing.
  All 63 metrics are covered, and covered *far enough to matter*: the weakest
  cover on the page is **10.00 %**, printed and asserted to clear the **5 %** at
  which `check_agreement` compares. That assertion is new, and it was written
  because one metric sat below it — the `S_1(0)` sensitivity, moved 0.23 % by its
  only mover, so no defect this page injects could have surfaced as a regression
  in it. It now has a row that truncates the sampled `S_1(0)` and re-root-finds,
  which moves it 92 %. The
  five below-floor metrics each carry a named above-floor companion, and the four
  genuinely structural quantities are named as identities — under the names the
  code prints them by, none of which is an `agreement.json` key;
- **and the rows are parsed, not trusted, to check they compute what they
  return**: a row that returned a *typed constant* for a metric would record a
  relative move of exactly 1.0 whatever the reported value was, so its coverage
  links could not fail — a hand-written coverage claim wearing the generator's
  clothes. This page shipped exactly that (five literal `0.0`s in the lag row,
  the only cover of four metrics), so a static guard now reads every row's own
  source plus one level of the helpers it names and rejects any `agreement.json`
  key bound to a numeric literal. Its teeth are **measured, not claimed**: the
  offending row is kept verbatim as a negative control and the guard is asserted
  to catch all five of its literal metrics — and to leave its `_`-prefixed
  diagnostic alone;
- **101 prose and metadata values** plus nine structural assertions checked
  against the live computation, **plus a mechanical sweep** of `meta.yaml`, this
  file, the data sidecar and `models_entry.yaml` that requires every number
  written to five or more decimal places to match a live value to **half an ulp
  of its own printed digits**. The token count depends on which shape the page is
  in and **both are pinned and asserted**: **93** across the four files in the
  case queue, **60** in the published page, where `models_entry.yaml` has been
  spliced into `models.yaml` and is counted as absent (on Colab all of them are).
  The notebook fails to execute if any token drifts, and the sweep **measures its
  own teeth** rather than claiming completeness: it rejects **100 %** of
  last-digit corruptions, and, hardened to every digit position by +1 and by +5,
  it rejects all but **one** substitution — the reported Richardson plug-flow
  outlet at `τ = 6` h, `0.6076462`, corrupted in its seventh decimal into the
  LSODA batch reference `0.60764725` that the page prints beside it, a millionth
  away. The code names that survivor by the live value it lands on, and a second
  distinct one would fail the build;
- **deterministic**: nothing stochastic, and the one sweep-shaped calculation
  uses a per-τ deterministic initial guess rather than a warm start, so no
  reported number can depend on sweep order. Two consecutive executions give a
  byte-identical `agreement.json` and an identical content hash over all cell
  sources and outputs, both figures included.

## Reuse

Read `K_s` and `K_i` as the constants of eq. (1), not as the half-rate
concentrations they are named for. Do not carry a Monod textbook's definition of
`r_m` across to the inhibition form. Compute with eq. (6)'s form, not eq. (1)'s —
they are the same function but eq. (1) is `0/0`-shaped at `S = 0`, where every
batch and startup run begins, and Andrews rewrote it for exactly that reason.
Check **both** washout conditions, not just the printed one. And root-find your
thresholds: a swept crossing would have agreed to three digits and hidden the
fact that Andrews' 1 hr ramp clears the critical ramp by only 14.67 %.
