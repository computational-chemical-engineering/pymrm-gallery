# D1.1 — Five fixed-bed models on one reactor

**Chapter 11 of Froment, De Wilde & Bischoff classifies fixed-bed reactor models in a
2×3 table and then works five of its six cells. This page builds all five on the same
reactor and measures what each one buys, as the shift in the root-found inlet partial
pressure at which the hot spot reaches the top of the book's own stated operating range.**

> *"If the partial pressure of the hydrocarbon were 0.018 bar, an increase of 0.0002 bar
> would raise the hot spot temperature beyond permissible limits."* — section 11.5.2,
> book page 513

That sentence is the page's anchor: the operating range is printed on the facing page
(335 °C to 415 °C), so "permissible limits" is a number, and the claim becomes a
root-findable threshold that tests the whole of rung 1 at once.

> **Units.** The sentence says *bar*. The model on this page is in **atm**, deliberately
> and throughout, and every pressure it reports is in atm. Section 11.5.2 restates Van
> Welsenaere & Froment's rate constant unchanged while relabelling their pressures bar;
> its own Fig. 11.5.2-1 labels the ordinate *"p₀, atm"*; and Example 11.5.3.A in the same
> chapter uses 0.208 where 11.5.2 prints 0.211 bar, which the notebook's audit shows is
> the same number converted. Measured both ways by one grid-free solve with one constant
> changed: in atm the sentence reproduces to **+0.028 %**; read in the printed bar it is
> out by **−1.105 %** and the hot spot at 0.018 is **627 °C**, 212 K *above* the ceiling
> the facing page states.

- **Structures:** `S2` `S4` `S6` `S7` `S8` — the widest of any page in the gallery
- **Covers:** `D1.1` (owning), `D1.2`, `D1.3`, `D1.4`, `D1.5`
- **Reference:** Froment, G. F., De Wilde, J. & Bischoff, K. B., *Chemical Reactor
  Analysis and Design*, 3rd ed., Wiley (2011), ISBN 978-0-470-56541-4, Chapter 11.
  Canonical-class source, read directly from 300–450 ppi renders of its own pages.
- **Provenance:** tier 6. **Chapter 11 reports no experiment.** Everything here
  reproduces the book's own arithmetic, its transcription of a 1970 paper, and five
  claims it makes about its own models.
- **Runtime:** about 37 minutes (measured: 2228 s and 2219 s for two concurrent executions)

## The five rungs, and why not in the book's order

Table 11.4-1 (book p. 505) is a 2×3 classification, not a ladder. The book presents the
five cases by reading it **down each column** — 11.5, 11.6, 11.7, 11.8, 11.9 — which puts
the two-dimensional pseudohomogeneous model between the two one-dimensional
heterogeneous ones. This page orders them by **how much each adds to the transport
description of the previous**, which puts 11.7 last:

| rung | section | structure | what it adds |
|---|---|---|---|
| 1 | 11.5 | `S2` | basic 1-D pseudohomogeneous plug flow |
| 2 | 11.6 | `S4` | + axial dispersion of heat and mass, Danckwerts inlet |
| 3 | 11.8 | `S7` | + a solid phase behind a film |
| 4 | 11.9 | `S8` | + concentration profiles inside the particle |
| 5 | 11.7 | `S6` | + a radial coordinate in the tube |

The reason is arithmetic, not taste. **Rung 1's wall coefficient and rung 5's two radial
parameters are the same measurement**: eq. (11.7.4-1) collapses (α_w, λ_er) into the α_i
that rung 1 uses, and the notebook reproduces section 11.5.2's printed
*U* = 0.096 kJ/(m² s K) from section 11.7.3's printed α_w = 0.156 and
λ_er = 0.78 × 10⁻³ **to both of the digits it prints**.

That is a consistency check and is stated as one — *not* evidence that the book derived
0.096 that way, since Van Welsenaere & Froment's own *U* = 82.7 kcal/(m² h K) = 0.09618
also rounds to 0.096 at two significant figures. What the identity does decide is the
tube **radius**: 0.096 at R_t = 0.0125 m against 0.0954 — which rounds to 0.095 — at the
0.0127 m the printed 2.54 cm diameter implies.

## What it finds

**1. Each rung, measured as a shift in the safe inlet partial pressure of the same
reactor.** The two heterogeneous rungs pull in *opposite* directions, and measured like
for like — rung 4 against **rung 3**, the model it corrects — the second correction is
**4.5 times** the first: adding a film makes the model 4.74 % more conservative, and then
resolving the particle it stands in front of moves +21.43 % the other way. (Against the
common rung-1 baseline the same pair reads −4.74 % and +15.67 %, a difference of
20.41 *percentage points* — a different quantity, and the page prints both.)
**Stopping at rung 3 leaves you further from rung 4 than rung 1 was.** Axial dispersion
moves the threshold by +0.1665 %, which is the same order as the 0.104 % discretisation
error of a single 1200-cell rung-1 solve — so it is only resolvable at all if both models
are Richardson-extrapolated, and the shift itself is refined over five grids.

**2. At section 11.5.2's own design point 0.018 atm, two of the five rungs are already
past their own runaway threshold** — rungs 3 and 5, whose "hot spot" there is 766 °C and
963 °C, not an operating point. The profile comparison is therefore made at section 11.7.3's own printed inlet
composition (0.00924, a *mole fraction*, hence 0.00924 atm at p_t = 1 atm), which is
subcritical on all five. The ladder table's one empty cell is a **solver** failure and is
labelled as one: at (n_z, n_r) = (1200, 20) rung 5's Newton does not converge at
0.018 atm, while on the coarser (600, 10) grid the break rows use it reaches a bounded —
and violently supercritical, 963 °C — state there. Nothing on the page rests on that
cell.

**3. Section 11.5.2's sensitivity sentence is right to 0.028 % in inlet partial pressure
— in atm, and as an agreement with Van Welsenaere & Froment's 0.018 atm** (see *Units*
above; read in the bar the section itself prints, the same sentence is out by −1.105 %).
Getting there settles two things the chapter leaves ambiguous. The tube radius is
0.0125 m, not the 0.0127 m its printed 2.54 cm diameter implies — at 0.0127 m the hot
spot at the book's own 0.018 is **662 °C**, far outside the range the same page states,
where at 0.0125 m it is **414.5 °C**, 0.5 K inside it — and 0.0125 m is *also* what
eq. (11.7.4-1) needs to return the printed *U*. Two independent constraints, one answer.

**4. The chapter's own criteria for deciding whether to climb from rung 1 to rung 2
cannot be evaluated on the reactor it designs.** (11.6-3)'s heat criterion divides by
(T₀ − T_w), which is zero here **by construction** — section 11.5.2 sets T = T₀ = T_r at
the inlet and section 11.5.3 builds its whole runaway diagram for "the common situation
where T_r = T₀". And (11.6-4)'s heat criterion is dimensionally inhomogeneous as printed:
a temperature per particle diameter compared with a dimensionless Péclet number, where
its own companion in (11.6-3) is normalised by (T₀ − T_w). **Section 12.7.2, 148 pages
later (book p. 709), says the criteria are not general** — *"Mears [1976] pointed out that
these criteria were not general, and he provided alternate ones for equal feed and wall
temperatures"* — with no cross-reference in either direction. That the alternates exist
*because of* the two problems above is this page's **inference** from that sentence, and
is labelled as one; what is not an inference is that the situation it singles out, equal
feed and wall temperatures, is exactly the one section 11.5.2 builds. The finding is the
missing cross-reference, not missing knowledge. The chapter's *conclusions* about axial
dispersion (the 50-particle-diameter rule, the Péclet thresholds on page 563) are all
confirmed by the measurement.

**5. Mears' interphase criterion passes by a factor 3.79 and the reactor is still 4.7 %
less safe than rung 1 says.** Mears' eq. 14 (published `B1.7`) asks whether the *observed
rate* deviates 5 % from the intrinsic one and correctly answers no. But the same film
difference — 0.85 K at section 11.8.2's printed startup temperature, and 3.31 K at rung
3's own threshold — moves the runaway threshold by 4.3 times the 1.11 % margin section
11.5.2's own sensitivity sentence hangs on. **The criteria in the chapter test the rate;
what a designer of this reactor needs tested is the threshold.**

**The chapter says the qualitative half of this itself, and the page now prints it
beside the finding.** Book p. 513, in the sentence immediately after the one quoted at
the top of this file: *"Note that for the upper part of the curves with
p_A0 = 0.0181, 0.0182, and 0.019 …, the model used here is no longer entirely adequate:
**heat and mass transfer effects would have to be taken into account**."* That is rung 3,
qualitatively, at exactly this operating point. The contribution here is the size and the
sign, and the observation that the chapter's own yardsticks three pages later are *rate*
criteria while its own warning is about the *threshold* — not the observation itself.

**6. As you climb the ladder the coefficients the book supplies for its own design case
run out — and then come back.** Rung 1's *U* is printed (and reproduced by
eq. (11.7.4-1) from rung 5's pair).
Rung 2's D_ea is "between 1 and 2" and its λ_ea is *"little information is available"*,
with Pe_ha used in two printed criteria and defined nowhere in the 902-page book. Rungs 3
and 4 are sent to Chapter 3's charts, to Satterfield [1970] and to Weisz & Hicks [1962].
Rung 5's λ_er and α_w are both printed. **The two ends of the ladder are pinned by the
book and the middle is not**, so the sign of the rung-4 correction is establishable and
its size is not — which is stated as a limit, not filled in.

So the page asks the question the other way round: *how far would each coefficient have
to move from the chapter's own value before its rung shifted the threshold by 1 %?* Each
answer is a controlled scan — swept outward from the chapter's value, every threshold
checked as a genuine 415 °C crossing, then refined between two points that passed that
check. **All four rows are roots and all four rungs order.** Rung 5 is closest to
irrelevance: a tube 1.25× narrower than this one (R_t* = 10.03 ± 0.03 mm, d_t/d_p = 6.69)
already brings the two-dimensional rung inside 1 %, which is the opposite of "never
avoidable". Then rung 2, at Pe_a* = 0.2488 ± 0.0002 — axial dispersion 4.02× stronger
than the chapter's own lower limit. Then rung 3, at 4.93× the correlated film
coefficients. Rung 4 is furthest, at 22.7× the D_e it was measured at. Rungs 2 and 3 are
23 % apart, which is the honest caveat on *that* pair.

**Two wrong versions of the rung-2 row are in the break table, and the second is the
interesting one.** The first located the crossover with the crossing control off and got
Pe_a = 0.097, where the shift is 2.8 % and both sides of the "threshold" are past the
ceiling. The repair put the control inside the scan — and then read its first refusal as
an edge, bisected it, and reported a *bound* of "at least 3.2×" with rungs 2 and 3
declared unorderable. That was the same class of error: **the control's verdict is
speckled in the swept parameter, not an interval.** On a regular 0.001 grid through the
region where the sweep first refuses, Pe_a = 0.316 and 0.314 validate while 0.315 and
0.313 do not, and the reported "edge" flipped its own verdict one ulp away. A refusal
means the bisection path to *that* threshold ran into a non-crossing; it says nothing
about the next point along. The crossover is on the same sweep, further out, on points
carrying the same ±1e-5 certificate as the five rungs of the ladder table — 414.98 °C
just below the threshold and 415.02 °C just above.

The same reading applies to the tube radius, with smaller stakes. Below the crossover the
*signed* shift passes through zero — **−0.169 % at 9.80 mm, +0.015 % at 9.75 mm** — so
the two-dimensional correction is negative on this tube and positive on a narrower one,
and |shift| stays under 1 % from 10.03 mm down to at least 9.675 mm. Rung 5 is speckled
too: 9.6875 mm is refused while 9.70 mm and 9.675 mm on either side of it validate, so
the previous version's "below 9.67 mm this solver says nothing" is withdrawn along with
the sign change "near 0.009 m" and the smallest |shift| of 0.093 % that an earlier
version read off points the control rejects.

**7. Fifteen printed defects, reported and none repaired.** Two are larger than the rest.

*The units.* **Section 11.5.2 is Van Welsenaere & Froment's parameter set converted from
atm to bar, and the conversion is done inconsistently.** Loading `D2.2`'s CSV of their
printed originals shows p_B0, the heat of reaction and all three of Table 11.5.3.A-2's
inlet limits multiplied by 1.01325 — the three table limits to five digits, the whole set
to 1.2e-3 relative — while the back-integrated critical value 0.01651 is left in atm and
labelled bar. The book calls the two *"in excellent agreement"*: as printed they differ
by +2.2 %, and **converting consistently improves the agreement to +0.84 %**. And the
*rate constant* was never converted at all, so the section's own design pressures are
atm wearing a bar label — which its own Fig. 11.5.2-1 confirms in its ordinate.

*Table 11.5.3.A-2's two columns disagree by up to 37 %, and this is the largest defect on
the page.* ΔT_ad = (B/A)p₀ with the example's **own** printed A = 6150 and B = 257 × 10⁶
requires 565.4 / 825.7 / 695.8 K where the table prints 411.5 / 521.1 / 466.3. Inverting
each row for the Q it needs gives 4.1102–4.1117 against the printed **Q = 3.4675** — which
is the Q of the *critical radius 0.0175 m* the same example computes two paragraphs
later, transplanted onto the 0.0125 m case tabulated. The book's own printed
C = Q²A·exp(−E/RT + b), with C = 2U/(c_p R_t), makes Q ∝ R_t^(−1/2) and closes the loop.
**A row-by-row transcription check cannot see this** — every row transcribes at ratio
1.000000, because the number was copied correctly from the wrong case — which is why the
page runs the cross-column test the table imposes on itself.

*And the cross-column test cannot see all of it either.* Table 11.5.3.A-1 passes that
test at 1.8 %, which the page used to read as evidence that the defect is *confined* to
A-2. It is not: inverting each table's own printed Q for the radius it implies — by the
same route, validated against Van Welsenaere & Froment's printed Q at two radii and
against their printed C, which pins R = 0.01753 m by itself — puts **both** tables' Q at
0.0175 m while both state 0.0125 m. A-1 passes only because its two columns are not
independent of each other, where A-2's p₀ column is Van Welsenaere & Froment's own. What
the A-1 pass does rule out is a wrong B/A, a wrong conversion or a misread formula, since
every one of those would move both tables.

The printed specific heat 0.992 kJ/(kg K) sits 5.7 % from what four routes require, three
of them inside the same chapter — though the notebook also shows those four are **not**
four independent determinations (two are the first rounded or rearranged), so the 0.94 %
"spread" is not a consensus and only the conclusion survives. Group *A* of Example
11.5.3.A is printed 6150 where its own formula gives 6165.

**8. Example 11.5.1.A contradicts the definition its own equation refers to.** It takes
d_p = √1.5 × 0.003 = 0.0037 m, the sphere of equal surface *area*. Eq. (11.5.1-9),
printed two pages earlier, is d_p = 6(1−ε)/a_v = 6 V_p/S_p, which for a cylinder with
d = h = 3 mm is exactly 0.003 m. On the chapter's own definition the Ergun pressure drop
rises by **+30.7 %**, from 0.319 bar to 0.417 bar. (`A1.1` owns Ergun's
constants, which are loaded and asserted equal to the printed ones; it carries neither
McDonald nor Hicks and does not discuss the equivalent-diameter ambiguity, so this branch
of the example is new.)

## What it does not own

The chapter overlaps six published pages, and this page loads their numbers rather than
restating them.

| the chapter's content | owner | what is done here |
|---|---|---|
| Van Welsenaere & Froment's runaway criteria (11.5.3) | `D2.2` | both CSVs loaded; Example 11.5.3.A checked *against the original*; criteria not rebuilt |
| Froment (1967)'s o-xylene network and its 1-D/2-D pair (11.7.3/11.7.4) | `C2.10` | both CSVs loaded and reconciled with the book's restatement; the 1967 case **not** re-run |
| Danckwerts' inlet/outlet conditions (11.6) | `A2.1` | its published dicts used unchanged; its upwind-dispersion finding applied, not rediscovered |
| the isothermal effectiveness factor (11.9.1) | `B1.1` | its exact-η CSV is the reference for eq. (11.9.1-11) |
| the Ergun friction factor (11.5.1-13) | `A1.1` | its printed-constant CSV loaded and asserted |
| Sh = 2 + 1.1 Sc^⅓ Re^0.6 | `A3.4` | used as published |
| Mears' criterion family | `B1.7` | its eq. 14 evaluated here and credited |

Also reconciled rather than annexed: the book states 365 °C for the equivalent
one-dimensional runaway limit and "less than 360" for the two-dimensional one, where
`C2.10` located both by bisection on the 1967 paper's own kinetics. **The book turns a
range into a number**; the two are printed side by side.

## Files

| file | what it is |
|---|---|
| `build_page.py` | generates `index.ipynb`; the physics lives here |
| `index.ipynb` | the executed page, nine sections in the required order |
| `data/froment-2011-ch11-printed.csv` | every constant, stated result and quantitative claim read out of Chapter 11, with the book page |
| `data/froment-2011-ch11-printed.meta.yaml` | provenance, the render-not-text-layer procedure, the defect list, and the search behind each negative claim |
| `meta.yaml` | page metadata |
| `agreement.json` | the metrics CI watches |

## Validation summary

- **Example 11.5.1.A** reproduced to 1.5e-3 relative on the printed d_p — three friction
  factors, three pressure drops, i.e. every digit it prints — and then shown to contradict eq. (11.5.1-9) by +30.7 %
  in Δp. A third, smaller inconsistency: no single Reynolds number produces all three
  printed friction factors (they need 155.5, 155.1 and 155.5 against a printed 155; the spread is 0.45).
- **Rung 1 computed a second, independent way** — adaptive Radau on the two ODEs, no
  grid, no pymrm operator, no Newton — agreeing with the Richardson-extrapolated pymrm
  threshold to **4.0e-07** and with the extrapolated hot spot to 3.3e-05, on the six-grid
  chain (150–4800) the notebook displays *and* `agreement.json` now carries: the display
  and the metric are the same memoised call. This is why the reported thresholds are the
  extrapolated ones: the 1200-cell value is 0.104 % off, comparable to the whole rung-2
  effect.
- **Both axes of the 2-D rung refined separately.** The radial axis converges at second
  order and the axial at first, so the two-dimensional model's error is set by the axis it
  is not named after — but both are tiny, and like for like (axis against an axis limit)
  the radial discretisation error at (1200, 20) is **−0.003 K** and the axial one
  **−0.029 K**. The 4.3 K axis-minus-bulk difference is the radial temperature profile the
  2-D model exists to resolve, not a grid error. Rung 5's *threshold*, by contrast, is
  refined on neither axis; its (600, 10)-vs-(1200, 20) gap of 0.107 percentage points is
  printed and no order is claimed for it.
- **Bisection, not brentq, and the reason is measured on all five rungs.** Above the
  threshold the upwind Newton stops converging and the hot-spot function returns NaN;
  brentq dies with *"The function value at x = … is NaN"* and, on the rungs where it
  survives, returns the bisection root to every digit. The notebook now runs that
  comparison on every rung rather than on the two it is cheapest to show.
- **The crossing control runs inside every threshold this page computes.** Because
  bisection counts a non-converged solve as supercritical, the number it returns could
  be a Newton-convergence boundary rather than a 415 °C crossing. `p0_critical` solves
  at p0*(1 ∓ 1e-5) and *raises* unless the hot spot is below the ceiling on one side,
  above it on the other, and continuous across it — so a threshold that is not a
  crossing cannot reach a table, a sweep or a bisection. **This is the repair the second
  review asked for**, and it was not cosmetic: the control used to run on the five
  default rungs only, and two published numbers came from thresholds it now rejects — a
  rung-2 "1 % crossover" located on an ignition jump, and two tube radii whose
  thresholds sat about 2 % below the real ones. The notebook counts them: **714 bisected
  thresholds checked, 93 rejected** as not being ceiling crossings, where the
  staged version checked five.
- **And a refusal is not an edge.** The third review found the repair above had created
  its own defect: the scans stopped at the first refusal and reported the stopping point.
  Validity is *speckled* in the swept parameter — the notebook prints a regular 0.001 grid
  in Pe_a on which accepted and refused points interleave, and a refused tube radius with
  validated radii on both sides of it — so a scan must sweep **through** refusals and
  refine by subdividing and keeping what validates, not by bisecting (roughly half the
  interior points of the rung-2 refinement are refused). A scan that cannot straddle its
  target between two validated points now *raises*; there is no bound fallback.
- **Eq. (11.9.1-11) checked twice** — against a spherical pymrm pellet solve (second-order
  convergent) and against `B1.1`'s 60 exact sphere rows, an independent transcription of
  the same closed form from a different document.
- **The book's own claim that the particle is practically isothermal** holds, and its
  second half holds more strongly: the film carries far more temperature difference than
  the particle does. The Prater number here is four orders below `B1.1`'s fold, so none of
  its multiplicity is in play and none is borrowed.
- **Group B = 257 × 10⁶ reproduced** from four other printed constants, which is what pins
  (−ΔH), ρ_B, p_B0 and the volumetric c_p together — and what convicts both the printed A
  and the ΔT_ad column of Table 11.5.3.A-2.
- **31 injected defects against 39 helper-computed metrics**, with the coverage map
  generated from the measured moves. Every metric has exactly one helper that computes it,
  and **every** helper default is asserted equal, to 1e-12, to a number printed earlier in
  the notebook — from one list, used both to check the numbers and to check that no metric
  escapes the check, so a metric cannot be added without an assertion. `agreement.json`
  and the break table now hold the **same** 39 metrics: six of them (three film/pellet
  temperature differences, Mears' margin and the two chapter-against-itself numbers) used
  to sit in `agreement.json` with no helper and no row, under a sentence asserting that
  rows moved them — which was false, and is replaced by two helpers and a new break row.
  The four rung-to-rung shifts are the single remaining exception, and it is printed
  rather than implicit: the break rows cannot afford the (1200, 20) grid 31 times, so
  `agreement.json` carries the (600, 10) values and a table gives both plus the gap (at
  most 0.107 pp). Two rows inject the two crossover-search defects this page was sent back
  for — the *uncontrolled* search, and the controlled one that *stopped at the first
  refusal* — verbatim, and measure how far each answer was. An AST guard (from published
  `J4.2`) rejects a metric key bound to a numeric literal, to arithmetic on literals, to a
  `float`/`int` cast of one, to a subscript of a frozen literal tuple, or to a local name —
  in the row's own body or a one-level helper — bound by a plain or tuple-unpacking
  assignment; its teeth are measured on **six** negative controls, one per laundering
  form, its false-positive control binds a *metric* key to a computed value, and
  **at least sixteen** forms it does *not* catch — a module-level constant among them —
  are printed and asserted to escape. The count and the wording have moved twice after
  verifiers wrote their own sabotage sets. The third pass found two escapes the page had
  not named: tuple unpacking, which the guard's own positive claim covered and which the
  guard now catches, and `abs(0.1657)`, which is constant-*foldable* but is a call rather
  than arithmetic — so the claim no longer says "constant-foldable" and the escape list is
  quantified as "at least". The fourth pass found that "a local name … that only ever
  carries one" was still a category rather than a list: annotated (`_a: float = 0.1657`),
  list-target (`[_a, _b] = 0.1657, 2.0`), nested-tuple (`(_a, _b), _c = (0.1657, 2.0),
  3.0`) and walrus (`(_a := 0.1657)`) bindings all escape the same way tuple unpacking used
  to, and are now named in the escape list instead of folded into that trailing clause.

## Honest limits

- Tier 6. There is no experimental comparison anywhere in Chapter 11, and none is invented.
- **Rung 4's magnitude is not decidable from the chapter.** D_e is not printed for this
  catalyst; the sign and the crossover are established, the size is not.
- Rung 2's magnitude (+0.1665 % extrapolated) depends on λ_ea, of which the chapter says
  "little information is available"; Pe_ma and Pe_ha are varied together, so the crossover
  at Pe_a* = 0.2488 ± 0.0002 is a crossover in *both* and not in λ_ea alone. Its accuracy
  is set by the grid, not by the refinement: the same-grid shift at a validated sweep
  point moves 0.0008 pp between 600 and 1200 cells, which on the local slope of the sweep
  is ±0.0002 in Pe_a — three figures, not the six the refinement tolerance would suggest.
- **What the rung-2 sweep does not reach.** The largest |shift| any control-validated
  point of that sweep carries is 1.19 %, so a 2 % rung-2 crossover is outside what this
  scan can locate and the break row that moves the target moves it to 1.1 %. That is a
  statement about the scan; nothing is claimed about what lies below it.
- The film coefficients come from `A3.4`'s correlation, not from the chapter, and the
  root-found crossover is at 4.93× the correlated coefficients — so that rung's size
  depends on which correlation you pick. (That number is 4.08 if the shift is measured
  against the Richardson-extrapolated rung-1 threshold instead of the same-grid one; the
  crossover table uses the same-grid definition in **every** row, which it did not
  before, and prints both readings of this row.)
- **The magnitudes of all four shifts belong to section 11.5.2's simplification.** The
  notebook shows its single "pseudo-first order" rate is k₁ alone — the wanted reaction —
  used as the rate of total o-xylene consumption, with the combustion route's much larger
  heat of reaction absent. Only the ordering and the signs transfer to a real reactor.
- Section 11.10 (the sixth cell of Table 11.4-1), sections 11.5.4–11.5.6, 11.7.5, 11.7.6
  and Examples 11.9.1.A/B/C are out of scope.
- **No figure is used**, with one stated exception that is a *unit* and not a datum: the
  ordinate label of Fig. 11.5.2-1, read on a 300 ppi crop as "p₀, atm", which is what
  decides the unit the whole page is in. No curve, tick or point is digitised, traced or
  read off any plot, and no page image exists in this directory.
- **Section 11.5.2 already says, qualitatively, what finding 5 measures** — that near the
  threshold "heat and mass transfer effects would have to be taken into account". The page
  prints that sentence beside the finding. What is new is the quantity and the sign.

## Reuse

- `Ladder` is the whole ladder; a different reaction means replacing `source`.
  `p0_critical` root-finds a threshold on cold-started solves whose initial guess is a
  fixed function of the inputs, so the number cannot depend on which direction a sweep was
  run in.
- **Put the control inside the root-find, not beside it.** A bisection that treats a
  non-converged solve as "past the threshold" returns a *solver* boundary wherever the
  solver gives out, and it does that silently. This page checked that boundary on the five
  configurations it displayed and not on the hundreds of thresholds its sweeps and
  bisections computed — so a crossover landed on an ignition jump and two swept rows were
  reported 2 % below their own thresholds. The repair is one line inside `p0_critical`: run
  the ±1e-5 control on **every** call and *raise* when it fails, so a bad threshold cannot
  become a number anywhere.
- **Then write the searches to expect refusals — and do not read the first one as an
  edge.** This is the part that is easy to get wrong, and this page got it wrong for a
  whole review cycle: the verdict is *speckled* in the swept parameter, so a scan that
  stops at the first refusal reports where it stopped as if it were a limit, and a
  bisection dies at its first refused midpoint. Sweep through refusals; refine by
  subdividing the bracket and keeping the points that validate; make the scan *raise* when
  it cannot straddle its target, because a bound is not an answer. Both wrong versions of
  the rung-2 search are break rows, and the table prints how far each was from the root.
- **The Danckwerts inlet as a length.** Because εD_ea/u_s and λ_ea/(ρ_g c_p u_s) are both
  lengths, rung 2's inlet condition is `{"a": d_p/Pe, "b": 1.0, "d": c_in}` with no unit
  bookkeeping, and the two fields differ only in their Péclet number.
- **Scale an algebraic residual to the units of its own unknown.** Written in the units the
  book writes them in, (11.8.1-3) is in bar/m and (11.8.1-4) in kcal/(m³ h) — 10⁷ apart —
  and the solve converges to a residual of 1e-10 *and a 123 K interfacial temperature
  difference*, with a smooth, plausible profile. Dividing each balance by its own transfer
  coefficient fixes it.
- **The helper-plus-injection pattern**: one function per metric, *every* reported number
  asserted against its default call as a checked set, every break row perturbing one of
  its arguments. The failure mode this page shipped once: an unasserted helper ran a
  coarser chain than the display, and `agreement.json` carried an observed order of 0.80
  where the page displayed 0.94.
- **Check a two-column table against itself.** Every row of Table 11.5.3.A-2 transcribes
  from its 1970 source at ratio 1.000000 and the table is still internally wrong by a
  third, because one row was copied correctly *from a different case*. The cross-column
  identity ΔT_ad = (B/A)p₀ is what found it.
- **Do not retype a constant from a textbook restatement of a paper.** Most of this page's
  design-case defects were found by loading `D2.2`'s and `C2.10`'s CSVs and printing them
  beside the book's — including the unit slip and the transplanted Q.
- **State which unit system a reproduction is in, and measure the alternative.** This page
  agrees with its source to 0.03 % in one unit and misses by 1.1 % in the other, and the
  source prints both.

## Regenerating

```bash
python build_page.py                    # regenerate index.ipynb from the builder
python ../../../scripts/run_pages.py    # execute it
python ../../../scripts/check_agreement.py
```
