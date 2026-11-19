# Refresh checklist

The rules every refreshed page meets: the refresh agent applies them, the reviewer
checks against them (citing rule numbers), the fixer must not break them. The contract
itself is `AGENTS.md` and `docs/claims-design.md`, cited by section. "(rows ...)" are
the evidence in `docs/refresh-learnings.md`. Generic pymrm rules point to the plugin
files under `~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/`: `pitfalls.md` (P1 to
P10), `numerics.md` and `style-guide.md` in `conventions/references/`, and
`verify-model/references/attack-list.md`.

## Claims and numbers

1. Every number in README, meta.yaml and models.yaml text carries
   `<!-- claim: name -->` right after it or sits inside an exact quote of old wording;
   only identifiers (catalog IDs, equation, table and figure numbers, dates,
   citations) are exempt. Notebook prose numbers are f-strings of computed values;
   refer to a printed count instead of restating it. Keep a tagged number clear of
   words ending in below, under, over or above. (AGENTS.md "Adding a page";
   claims-design 6; rows 23, 49, 55, 164)
2. Counts, labels and statements over a set ("inside every bracket", "all but two
   rows", "dominates only at the high end") are computed and registered as `count` or
   `exact` claims with their margin. A count typed into an f-string goes stale: compute
   it (`A4.5` said 2 of the free rows fail, 1 of 318 does; row 261). (claims-design 3;
   rows 23, 45, 91, 160, 162, 172, 261)
3. Classify by role, then by measured spread, never by size: residuals, identities,
   two-route gaps and tolerance effects are bounds with a basis. Print a
   round-off-level value as its bound with `(here ...)`; never divide by one. A
   "finite-difference floor" is truncation plus round-off `eps |f| / h`: estimate both
   before naming it (`A4.7`: gap 2.5e-6, truncation 2.2e-7, the rest round-off; row
   265). (claims-design 2, 3; AGENTS.md "Adding a page"; rows 1, 2, 6, 129, 265)
4. A seeded or sampled quantity gets enough draws that the seed moves it little, and a
   tolerance of about three times the spread measured over a dozen seeds; added draws
   get their own generator. A verdict ("excludes 1") is its own bound claim on the
   distance from the excluded value. (rows 157, 159)
5. Compare on the same basis: both sides of a ratio computed the same way (never one
   computed and one printed), on matched grids and tolerances; name the denominator;
   "off by a factor" divides by the bound. (rows 21, 52, 66, 73, 96, 106, 144)
6. A value labelled at a parameter value is computed there, not at the nearest plotting
   point. "Within the printed precision" of a source is a count in half-units of its
   last digit, read from the column as strings (`dtype=str`): count the cells that round
   to the printed entry, never quote a largest relative deviation (`A4.1` "to 5.80 %"
   was 1 of 16 cells; row 245). "Consistent with rounding" compares each row with the
   half units of every printed entry it combines (`A4.6` O2-Xe at 1.36 of its
   allowance; row 263); a mean of printed values also gets the movement rounding of its
   inputs can cause (`A4.5` N2 mean at 1.023 half-units against 0.1; row 261). (rows 44,
   46, 124, 143, 245, 261, 263)
7. Never retype a CSV value, another page's claim (use `load_claims`) or a hand read no
   CSV holds; a hand read belongs to the page that makes it. A threshold from another
   page is its claim, not its rounded prose (`C1.2` used 16.9 % for 16.856 %: window
   span 40.0009 against 39.768; row 250). (AGENTS.md "Data rules"; claims-design 5; rows
   174, 178, 250)

## Break rows

8. Every non-structural claim has a break row that fails that claim's own check, sized
   against its tolerance; a map of claims to defect classes is not a break table.
   Compute guards in the helper the rows call. `structural` is for a named identity no
   row can move. Measure which defect size crosses the bound before writing the row
   (`C1.2` step guard: 1e-2 moved the orders 2.6e-9 against a bound of 1e-8, 1.0 moves
   them 2.8e-5; row 249). A limiter of higher order (STOIC, orders 2.28 to 2.37 on
   `A3.8`) is the break row of an order's upper bound (row 240). (AGENTS.md "Adding a
   page"; claims-design 3, 5; rows 7, 54, 119, 240, 249)
9. Inject the defect through a call argument or physical parameter, everywhere the
   parameter enters, post-processing included; never edit a result after the solve or
   hand back the baseline value. Build labels from the call. A defect a mask can absorb
   (one bad cell) is too small. A looser `tol` is no defect for a `newton` cross-check
   on a well-conditioned root (`C1.2`: tol 1e-1 to 1e-8 all end at 4.8e-11); perturb
   what newton solves, or `maxfev` (row 248). (rows 28, 86, 154, 158, 191, 248)
10. A defect's effect is the difference of two solves that differ only by the defect.
    A comparison of an approximation with an exact solve has a row where the
    approximation is exact by construction. (rows 70, 78)
11. Convergence is a precondition, not a detector. A defect that leaves no solution
    records nan, with the non-existence shown by an independent route or a root scan
    over the whole physical range. (rows 61, 72, 156, 192)
    The row stays in the table: declare the rows with `claims.define_rows(labels)`
    before they run, and a study that raises returns nan for every claim the row
    targets, never an empty result. `claims.report()` fails for a declared row that
    recorded nothing, and accept_page.py prints "break rows: N defined, N recorded".
    A summary such as "rows that fail no claim: none" cannot see a row that recorded
    nothing. (A3.7 review, 9 Nov: 42 rows defined, 40 printed)

## Numerics and evidence

12. Crossings, thresholds and band edges are roots: `pymrm.checks.find_roots` on the
    physical range, raising unless it returns exactly one root; never a grid read-off,
    a sampled extreme, the nearest log-grid node or a capped scan. Raise when the
    domain is too short to hold the feature. Scan the fundamental domain of a symmetric
    model. An event of a time march (a flux that vanishes) is a root in time of an
    interpolant through the step values, then extrapolated in dt after its order is
    checked on three levels (`A4.2`: argmax 7.52 h, spline root 7.555 h, limit 7.5040
    h; row 247). (numerics.md "Bracket before you root-find"; rows 17, 33, 65, 81, 156,
    169, 182, 247)
13. An extremum is a root of an analytic derivative (pole multiplied out), or a bounded
    maximisation refined from a scan of the whole window; an argmax on a domain end
    means "at or beyond the end". Take limits in closed form, not at the grid end; at a
    finite point a limit is a rate, a theory bound. A maximum on the production grid is
    also extrapolated (rule 36; `A4.3`: sweep node 562 nm, maximum 744 nm; row 256).
    (rows 22, 96, 105, 131, 147, 173, 189, 256)
14. Orders come from successive differences of at least three levels
    (`pymrm.checks.observed_orders`), never against the finest grid or the exact value
    on a one-axis ladder. Report each case as a theory range with the ladder printed
    (minimum and maximum where they oscillate), plus the round-off companion bound.
    Check the companion case by case before pooling: a nearly linear exact profile
    measures round-off and can still give orders near 2 (`A4.1` Examples II and III,
    finest differences 1.5e-11 and 3.9e-12 against a non-uniformity of 1.9e-11 and
    5.3e-11; row 244). A refinement factor multiplies every term of a step schedule
    (`A4.7`: one unscaled term, order -13.8; row 266). (claims-design 3; rows 9, 19, 43,
    75, 95, 132, 134, 163, 244, 266)
15. Split an error before naming its cause: time against space (exact time by the
    eigendecomposition of the constant operator, or the other axis cancelled;
    observations on step times `(i + 1) * dt`), approximation against discretisation,
    each axis of a 2-D grid on its own ladder. Record every history of a march at the
    same times, initial state included, and pair by time, not index (`A4.2`: spread
    5.4 % paired by index, 4.5 % by time; row 246). (rows 48, 74, 82, 83, 89, 122, 126,
    133, 181, 246)
16. Name a cause only after removing it and watching the effect go, one ingredient at
    a time; a cause named for a gap is a claim, measured with the term added as a
    parameter. Never attribute an effect to one unmeasured cause. (rows 37, 79, 87,
    123, 166, 188)
17. Measure the reference first: scaled integrator `atol` and a second tolerance, a
    root-found reference checked at its root, series truncation at the largest
    parameter, `quad` tolerances set, cancelling closed forms against `mpmath`; evaluate
    it at cell centres, not through `np.interp`. Richardson removes discretisation, not
    an estimator's bias. A quadrature's error is the change when all its knobs refine
    together (`A4.6`: n_E alone 0.164 %, n_b alone 0.151 %, both 0.211 %; row 262).
    (rows 36, 63, 77, 86, 98, 117, 125, 135, 137, 262)
18. A statement over a range or a set of cases is evaluated at every case and point it
    covers, on the page's grids; root-find where a comparison flips. State conditions
    with the range ("for one light gas throughout"). (rows 94, 142, 143, 145, 170, 176,
    190)
19. A mean bias needs a standard error before it is called anything, and a reference
    group that is what it is said to be; split data where the abscissa certifies it.
    See rule 35 for the test a claim needs. (rows 176, 177, 179)
20. State what two routes share (end conditions, interpolant, factorisation, source
    function); routes sharing an interpolant agree to the coarser root finder. Measure
    a second `permc_spec` before claiming solver independence. (attack-list.md
    section 2; rows 59, 150, 153)
35. Test a statistical claim with the test that matches it. A claim about one
    predictor is a one-sample test of it, not a ratio of mean errors (`A3.8`: one-sample
    t 3.50, paired t 0.80; row 241); a fit gap between two models on the same runs needs
    the paired standard error (`C1.2`: 9.82 points at 1.57 standard errors; row 251); a
    count inside a band needs each estimate's standard error (`A4.5`: 20 of 26 inside,
    6 by two standard errors; row 260); "consistent with" a figure reading counts the
    points inside their sidecar reading error (`A3.12`: 4 of 14; row 243). A slope over
    ratios is printed with and without each named anomalous row, excluded by identity
    (`A4.5`: 1.74920 and 1.75001; row 259).
36. A quoted value at a condition, a difference between two parameter values on one
    grid, and an extremum of a sweep carry the grid error at that condition: extrapolate
    each on its own ladder before stating it (`A4.4`: 600-cell minus extrapolated 1.1e-5
    at phi = 10, 2.3e-3 at phi = 100; `A4.3`: phi = 200 against 60 on one grid looked
    larger, extrapolated both are -1/(2 Kn)). (rows 255, 256, 257)

## pymrm modelling (the plugin owns the rule)

21. A first-order convection-only outlet is `{"outflow": True}`; an outlet bc of `None`
    is a Dirichlet zero. Read outlet, wall and boundary values and fluxes with
    `compute_boundary_values` under the same bc, never from the last cell, which is
    half a cell short (`A4.3`: 24.43 % from faces on every grid, 22.5 to 24.2 % from
    cells; row 253). A van Leer profile with an outlet gradient, under outflow or zero
    gradient, has a first-order layer of 2 to 3 cells at the outlet while the interior
    stays second order: measure an order on the outlet face value or at a fixed x,
    never as a maximum over all cells. Under a per-field zero gradient the face value
    is second order and the last cell centre first order (`A3.8`: face orders 1.93 to
    1.99, last cell 1.14 to 1.02; row 239); this corrects the earlier wording that the
    last cell is a second-order outlet value. (P6, P8; rows 32, 35, 87, 121, 141, 148,
    149, 166, 187, 239, 253; pymrm #41)
22. Field shapes `(n, 1)`, `NumJac((n, 1))`, spatial axes first. A rate whose
    derivative is unbounded or varies fast gets its exact Jacobian through
    `construct_coefficient_matrix`; a state-dependent velocity is convected as a molar
    flux with a constant unit-velocity operator. (AGENTS.md "pymrm conventions you must
    follow"; P3, P9; rows 101, 110, 184)
23. Check every solve: `newton` `success` (`NewtonBreakdown` raises instead),
    `solve_ivp` `status`; `tol` above the round-off floor with `stagnation=True`;
    bounds rather than a custom step cap, tested apart from damping. An unchecked
    `success` can publish a number (`A4.4`: the 1.0e0 came from a solve that stopped on
    a bound; row 258). A fixed point stalling at an unchanged contraction rate has a
    step tolerance below round-off (`A3.12`: 1e-10 degC against a median step 7.0e-10;
    row 242). Measure the tolerance effect before quoting a small departure (`A4.3`:
    tenfold tighter moves 5.4e-5 by 2.2e-7; row 254). (P7; style-guide.md section 11;
    rows 13, 26, 47, 57, 104, 139, 242, 254, 258)
24. Newton converges to a nearby root, stable or not: locate a threshold by following
    the branch with warm starts, and say which branch a state is on. Choose
    `ConvFluxTVD` with `newton` or deferred correction by measuring both (`A3.8`: STOIC
    by deferred correction does not converge in 200 iterations, by newton it does; row
    240). (P10; rows 8, 18, 71, 76, 110, 120, 240)
25. Test a discrete identity on the scheme the page uses: exact cell volumes, a
    pure-Neumann source with zero mean on the grid, face gradients from
    `compute_boundary_values`. Its gap is the solve's residual, not a second route.
    (rows 31, 62, 69, 93, 183)
37. `construct_boundary_value_matrices` returns the inhomogeneous part with the bc's
    `d` in it; never multiply it by `d` again. A constant-coefficient validation cannot
    see the slip, so test a face value under a state-dependent coefficient (`A4.7`: a
    face at 0.9 read 0.81, Sh 35.47 against 35.83; row 264).

## Corrections

26. A Correction is for a published statement that was wrong. A change of definition,
    a fair rounding, a gain in precision or a true statement made more precise is not
    one; it goes in the commit message. A true old statement needs none. (AGENTS.md "A
    page shows current pymrm, not history"; rows 92, 97, 171)
27. Each Correction is dated, quotes the old wording character for character
    (emphasis markers may go; end a quote before an old em dash), and states what is
    true with claim markers. Old numbers appear only inside the quote. Every published
    number changed for any reason outside rule 26 has a Correction line; reproduce the
    old value first. (AGENTS.md
    "A page shows current pymrm, not history"; row 92)
28. A failure mode, trap or alternative result stated in prose needs the cell that
    produces it and a claim, or it goes (`C1.2` "a smaller number" from minimising
    mean-\|Δ\|, which no cell did: 22.76 %; row 252). (rows 53, 60, 252)

## Wording

29. "What pymrm adds" says what pymrm contributes on this page and no more; if it adds
    nothing to a correlation, say so. (AGENTS.md "Honesty requirements")
30. Say what a number is: what a sensitivity holds fixed (re-solve the coupled
    system), which statistic, which denominator, which uncertainties a margin covers. A
    statement about one grid ("cannot resolve") or one solver ("Newton cannot reach")
    is not one about the physics. (rows 18, 23, 29, 82, 171, 180)
31. No history narration and no pymrm version names in the notebook, comments,
    meta.yaml or models.yaml; a second route is described by what it computes. No em
    dashes anywhere. (AGENTS.md "A page shows current pymrm, not history")
32. A learnings row states only what was measured, with the numbers, and quotes the old
    page exactly. (row 92; ../SKILL.md "Learnings rows")
38. A learnings row has all 7 columns (# | Date | Page | Finding | Evidence | Target |
    Status), with literal pipes escaped as `\|`. A flush selects rows by Target and
    Status: on 2026-11-19, 15 rows (267 to 272, 275, 286 to 289, 298 to 301) lacked both
    and were skipped.

## Page structure and files

33. Files and layout as in AGENTS.md "Adding a page". Before an hour-long kernel run,
    check moved code against the old definitions on coarse grids. (row 56)
34. On a closed-form page every published constant sits in one dictionary that every
    function takes, so a break row per constant is one call; run each study separately
    inside the row loop and list the refusals. (row 175)

## Before you report

From the checkout, with `P=pages/<PAGE>`:

1. `python scripts/accept_page.py $P --only dashes,lint` runs three scans:
   Correction quotes against the old page, every quote found (rule 27); em dashes in
   the page, its models.yaml entry and its learnings rows, and pymrm version names
   in `$P/*.py` and `$P/meta.yaml` (rule 31); history words there, printed as LOOK
   lines: each is about the computation, never the page's history (rule 31); README numbers
   outside quotes without a claim marker, printed as LOOK lines: each hit is tagged,
   quoted or an identifier (rules 1, 27).
2. `grep -nE 'argmax|argmin|np\.interp|brentq|fsolve|minimize_scalar|\[-1\]' $P/*.py`:
   each hit is a checked `find_roots` result, a named data row, or an outlet read for
   `compute_boundary_values` (rules 12, 13, 21).
3. `grep -nE 'claims.breaks|newton\(|solve_ivp\(' $P/*.py`: each defect enters
   through an argument, each solve's `success` or `status` is checked, and each order
   claim has its round-off companion (rules 9, 14, 23).
4. Reread "What pymrm adds", the Corrections and the new learnings rows against the
   printed claims (rules 26 to 32, 38).
