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
- **Runtime:** about 19 minutes (1116 s per execution, single-threaded BLAS, measured 2026-10-31)

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
past their own runaway threshold**: rungs 3 and 5. The profile comparison is therefore
made at section 11.7.3's own printed inlet composition (0.00924, a *mole fraction*, hence
0.00924 atm at p_t = 1 atm), which is subcritical on all five. The ladder table prints no
hot spot at 0.018 atm for those two rungs, only "past threshold": that would be a
cold-started fixed-inlet solve above a threshold, whose answer (ignited state, or no
convergence) is decided by round-off and differs between BLAS kernels. On the coarser
(600, 10) grid rung 5 there reaches a bounded runaway state above 900 °C on every kernel
tried, which shows the reactor has a steady state there; nothing on the page rests on its
exact temperature.

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

**The chapter says the qualitative half of this itself, and the page prints it
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
irrelevance: a tube 1.25× narrower than this one (R_t* = 0.01003 m <!-- claim: crossover_R_t -->, ± 0.03 mm, d_t/d_p = 6.69)
already brings the two-dimensional rung inside 1 %, which is the opposite of "never
avoidable". Then rung 2, at Pe_a* = 0.2488 ± 0.0002 — axial dispersion 4.02× stronger
than the chapter's own lower limit. Then rung 3, at 4.93× the correlated film
coefficients. Rung 4 is furthest, at 22.7× the D_e it was measured at. Rungs 2 and 3 are
23 % apart, which is the honest caveat on *that* pair.

**How the thresholds are located.** `p0_critical` makes one cold solve at the bottom
of the bracket and starts every later scan point and bisection midpoint from the
converged state at the current lower end, so it follows the branch the reactor is on as
the feed is raised. A cold start from the rung-1 profile at the same p0 is not neutral:
between rung 1's threshold and a higher one that guess is already a runaway profile, and
Newton then lands on the ignited state (about 990 °C) or fails, by round-off. On one
machine (Haswell kernel), direct starts on rung 2 at n_z = 900, Pe_a = 0.24 and 0.01817,
0.01818, 0.01819 and 0.01820 atm gave 990.9 °C, 411.5 °C, no convergence and 992.1 °C;
the page does not recompute them, since by this argument they differ between machines.
With the branch followed, the rung-2 crossover is 0.24883 <!-- claim: crossover_pe_axial -->,
inside its ±0.0002; Pe_a = 0.30 is a clean crossing (414.98 / 415.02 °C, shift
+0.8288 %); and 0 <!-- claim: scan_points_rejected --> thresholds of the four scans
are rejected, on three BLAS kernels.

Below the crossover the *signed* shift passes through zero (**−0.169 % at 9.80 mm,
+0.015 % at 9.75 mm**), so the two-dimensional correction is negative on this tube and
positive on a narrower one. |shift| stays under 1 % from 10.03 mm down to 9.60 mm
(+0.572 %), and is past 1 % again, on the positive side, by 9.40 mm (+1.327 %); at
9.00 mm it is +2.879 %. Every one of those thresholds passes the control.

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
test at 1.8 %, which does not mean the defect is *confined* to A-2: inverting each table's own printed Q for the radius it implies — by the
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
| `index.py` | the edited source of the page (jupytext percent format), paired with `index.ipynb` |
| `index.ipynb` | the executed page, nine sections in the required order |
| `model.py` | the design case (`Reactor`), the five rungs (`Ladder`), the pellet, the threshold search with its crossing control, the crossover scan and the grid-free rung-1 route |
| `test_model.py` | import, fast solves of rung 1 and the pellet, and one stored claim, run offline in seconds |
| `data/froment-2011-ch11-printed.csv` | every constant, stated result and quantitative claim read out of Chapter 11, with the book page |
| `data/froment-2011-ch11-printed.meta.yaml` | provenance, the render-not-text-layer procedure, the defect list, and the search behind each negative claim |
| `meta.yaml` | page metadata |
| `agreement.json` | the page's claims in format 2 (values, bounds and exact counts), which CI checks by kind |

## Validation summary

- **Example 11.5.1.A** reproduced to 1.5e-3 relative on the printed d_p — three friction
  factors, three pressure drops, i.e. every digit it prints — and then shown to contradict eq. (11.5.1-9) by +30.7 %
  in Δp. A third, smaller inconsistency: no single Reynolds number produces all three
  printed friction factors (they need 155.5, 155.1 and 155.5 against a printed 155; the spread is 0.45).
- **Rung 1 computed a second, independent way** — adaptive Radau on the two ODEs, no
  grid, no pymrm operator, no Newton — agreeing with the Richardson-extrapolated pymrm
  threshold to **4.0e-07** and with the extrapolated hot spot to 3.3e-05, on the six-grid
  chain (150–4800) the notebook displays *and* `agreement.json` carries: the display
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
- **Bisection, not brentq.** Above the threshold the upwind Newton may stop converging
  and the hot-spot function returns NaN, on which brentq raises *"The function value at
  x = … is NaN"*. The notebook runs brentq on all five rungs: it fails on at least one,
  and where it survives it returns the bisection root to better than 1e-10. *Which* rungs
  it survives is decided by round-off (its points above the threshold are cold-started
  fixed-inlet solves) and differs between BLAS kernels, so the page does not say which.
- **Each threshold is followed along its branch, and the crossing control runs inside
  every one.** `p0_critical` makes one cold solve at the bottom of its bracket and starts
  every later solve from the converged state at the current lower end, so the threshold
  is the one on the branch the reactor follows as the feed is raised, and no state crosses
  from one call to the next. Because bisection counts a non-converged solve as
  supercritical, the number it returns could still be a Newton-convergence boundary, so
  `p0_critical` solves at p0*(1 ∓ 1e-5) from that branch state and *raises* unless the hot
  spot is below the ceiling on one side, above it on the other, and continuous across it.
  The notebook counts the checks: **0 rejected** over every bisected threshold on the
  page, the closest verdict about 5e-3 K from its threshold, registered as an exact claim
  with that margin.
- **The searches still tolerate a rejection, and meet none.** The crossover scans sweep
  through a rejected point, refine by subdividing and keeping what validates, and *raise*
  when they cannot straddle their target. On a regular 0.001 grid in Pe_a all 9 points
  validate and the shift rises monotonically.
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
- **34 injected defects against 43 helper-computed metrics**, with the coverage map
  generated from the measured moves, and **43 claims** registered with `Claims`: 32
  values, 8 bounds and 3 exact counts, none structural. Every metric has exactly one helper that computes
  it, and **every** helper default is asserted equal, to 1e-12, to a number printed
  earlier in the notebook, from one list used both to check the numbers and to check
  that no metric escapes the check. Every break row is fed to the claims it touches, and
  every claim has a row that fails its own check, including the guard on the rung-1
  order fit, which the row that drops the Richardson limit sends to zero. Residuals, identities, two-route gaps and
  convergence orders are bounds with a stated basis: the rung-1 ODE-versus-pymrm gaps
  below 1e-5 (threshold) and 1e-4 (hot spot), the orders in [0.9, 1.1] and [0.85, 1.15].
  The four rung-to-rung shifts are computed on (600, 10), because the break rows cannot
  afford the (1200, 20) grid once per row; a table gives both plus the gap (at most
  0.107 pp). Among the break rows: the rung-2 shift against rung 1's extrapolated
  threshold, the crossovers on a 300-cell grid, and two rung-1 rows (no Richardson
  limit; grid ratio misread) that give the bounds a row that fails them; a row running the ladder's controls at
  the neighbouring rung's threshold fails `ladder_controls_rejected`; and a row reading
  the control's jump allowance as 0.005 K (below every measured jump, 0.012 to 0.095 K)
  re-runs the five ladder controls and the 9-point fine Pe_a grid and fails the two
  page-wide rejection counts. An AST guard (from published
  `J4.2`) rejects a metric key bound to a numeric literal, to arithmetic on literals, to a
  `float`/`int` cast of one, to a subscript of a frozen literal tuple, or to a local name —
  in the row's own body or a one-level helper — bound by a plain or tuple-unpacking
  assignment; its teeth are measured on **six** negative controls, one per laundering
  form, its false-positive control binds a *metric* key to a computed value, and
  **at least sixteen** forms it does *not* catch — a module-level constant among them —
  are printed and asserted to escape. `abs(0.1657)` is constant-*foldable* but is a call
  rather than arithmetic and escapes, so the claim does not say "constant-foldable", and
  the escape list is quantified as "at least". Annotated (`_a: float = 0.1657`),
  list-target (`[_a, _b] = 0.1657, 2.0`), nested-tuple (`(_a, _b), _c = (0.1657, 2.0),
  3.0`) and walrus (`(_a := 0.1657)`) bindings all escape, and are named in the escape
  list.

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
- **What the rung-2 sweep does not reach.** The sweep ends at Pe_a = 0.10, where the
  shift is 2.52 %; nothing is claimed about smaller Péclet numbers.
- The film coefficients come from `A3.4`'s correlation, not from the chapter, and the
  root-found crossover is at 4.93× the correlated coefficients — so that rung's size
  depends on which correlation you pick. (That number is 4.08 if the shift is measured
  against the Richardson-extrapolated rung-1 threshold instead of the same-grid one; the
  crossover table uses the same-grid definition in **every** row, and prints both
  readings of this row.)
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
- **Locate a threshold by following the branch, not by cold starts.** Near a runaway
  threshold several steady states exist and Newton's damping does not choose between them;
  when the initial guess is itself a runaway profile, round-off does. One cold solve at a
  safe parameter value, then every scan point and midpoint started from the converged
  state at the current lower end of the bracket. No state crosses from one search to the
  next, so the number still cannot depend on the direction a sweep was run in.
- **Put the control inside the root-find, not beside it.** A bisection that treats a
  non-converged solve as "past the threshold" returns a *solver* boundary wherever the
  solver gives out, and it does that silently. Run the ±1e-5 control on **every** call and
  *raise* when it fails, record its margin, and write the searches so that a rejected
  point is used for nothing and a scan that cannot straddle its target raises.
- **The Danckwerts inlet as a length.** Because εD_ea/u_s and λ_ea/(ρ_g c_p u_s) are both
  lengths, rung 2's inlet condition is `{"a": d_p/Pe, "b": 1.0, "d": c_in}` with no unit
  bookkeeping, and the two fields differ only in their Péclet number.
- **Scale an algebraic residual to the units of its own unknown.** Written in the units the
  book writes them in, (11.8.1-3) is in atm/h and (11.8.1-4) in kcal/(m³ h), so one absolute
  residual threshold (the 1e-6 test in `Ladder.solve`) is strict for one row and loose for
  the other. Dividing each balance by its own transfer coefficient puts both in the units of
  the surface unknown they determine.
- **The helper-plus-injection pattern**: one function per metric, *every* reported number
  asserted against its default call as a checked set, every break row perturbing one of
  its arguments. The failure mode it prevents: an unasserted helper running a coarser
  chain than the display, so that `agreement.json` carries a number the page does not
  show.
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

From the repository root:

```bash
jupytext --sync pages/D1.1-froment-fixed-bed-ladder/index.py   # after editing index.py
python scripts/run_kernels.py pages/D1.1-froment-fixed-bed-ladder  # three BLAS kernels
python scripts/check_agreement.py
GALLERY_OFFLINE=1 PYTHONPATH=shared python -m pytest pages/D1.1-froment-fixed-bed-ladder/test_model.py
```

## Corrections

- 2026-10-29: the rung-2 crossover was stated as 0.24883362; it is 0.24883 (±0.0002).
- 2026-10-29: the crossing control's verdict was stated as "speckled" along the Pe_a
  sweep and the tube radius, with 13 of 24 refinement probes refused; every point of
  every scan is a clean crossing.
- 2026-10-29: Pe_a = 0.30 was stated as an "ignition jump" (414.2 against 999.3 °C); it
  is a clean crossing (414.98 / 415.02 °C).
- 2026-10-29: 9.6875 mm and 9.0 mm were stated as rejected thresholds; both are clean
  crossings.
- 2026-10-29: the largest validated rung-2 shift was stated as 1.19 %; the sweep
  reaches 2.52 %.
- 2026-10-31: the notebook said that without scaling the film balances the solve converges
  to a residual of 1e-10 and a 123 K interfacial temperature difference; with the rows
  unscaled `pymrm.newton` reaches the same state as with them scaled (0.5229 K at
  0.00924 atm on 1200 cells, gap 1e-13 K), as it must: multiplying a row by a constant
  does not change the Newton step. What the scaling decides is the absolute residual test
  after the solve.
- 2026-10-31: the notebook's Reuse said the AST guard has five negative controls and eleven
  escape forms; it runs six and sixteen.
- 2026-10-31: The data said five borrowed datasets; the page loads six, from four pages.
- 2026-10-31: `meta.yaml` gave R_t* = 10.028 mm; the located value is 10.029 mm
  (± 0.030 mm).
- 2026-10-31: The data said C2.10's `froment-1967-stated-results.csv` has 18 rows; it has 19.
- 2026-10-31: this README's Reuse advice said (11.8.1-3) is in bar/m; as the book writes it, it is
  in atm/h (rate of the surface reaction per unit bed volume, in partial-pressure units).
