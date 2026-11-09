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
   `exact` claims with their margin. (claims-design 3; rows 23, 45, 91, 160, 162, 172)
3. Classify by role, then by measured spread, never by size: residuals, identities,
   two-route gaps and tolerance effects are bounds with a basis. Print a
   round-off-level value as its bound with `(here ...)`; never divide by one.
   (claims-design 2, 3; AGENTS.md "Adding a page"; rows 1, 2, 6, 129)
4. A seeded or sampled quantity gets enough draws that the seed moves it little, and a
   tolerance of about three times the spread measured over a dozen seeds; added draws
   get their own generator. A verdict ("excludes 1") is its own bound claim on the
   distance from the excluded value. (rows 157, 159)
5. Compare on the same basis: both sides of a ratio computed the same way (never one
   computed and one printed), on matched grids and tolerances; name the denominator;
   "off by a factor" divides by the bound. (rows 21, 52, 66, 73, 96, 106, 144)
6. A value labelled at a parameter value is computed there, not at the nearest plotting
   point. "Within the printed precision" of a source is a count in half-units of its
   last digit, read from the column as strings (`dtype=str`). (rows 44, 46, 124, 143)
7. Never retype a CSV value, another page's claim (use `load_claims`) or a hand read no
   CSV holds; a hand read belongs to the page that makes it. (AGENTS.md "Data rules";
   claims-design 5; rows 174, 178)

## Break rows

8. Every non-structural claim has a break row that fails that claim's own check, sized
   against its tolerance; a map of claims to defect classes is not a break table.
   Compute guards in the helper the rows call. `structural` is for a named identity no
   row can move. (AGENTS.md "Adding a page"; claims-design 3, 5; rows 7, 54, 119)
9. Inject the defect through a call argument or physical parameter, everywhere the
   parameter enters, post-processing included; never edit a result after the solve or
   hand back the baseline value. Build labels from the call. A defect a mask can absorb
   (one bad cell) is too small. (rows 28, 86, 154, 158, 191)
10. A defect's effect is the difference of two solves that differ only by the defect.
    A comparison of an approximation with an exact solve has a row where the
    approximation is exact by construction. (rows 70, 78)
11. Convergence is a precondition, not a detector. A defect that leaves no solution
    records nan, with the non-existence shown by an independent route or a root scan
    over the whole physical range. (rows 61, 72, 156, 192)

## Numerics and evidence

12. Crossings, thresholds and band edges are roots: `pymrm.checks.find_roots` on the
    physical range, raising unless it returns exactly one root; never a grid read-off,
    a sampled extreme, the nearest log-grid node or a capped scan. Raise when the
    domain is too short to hold the feature. Scan the fundamental domain of a symmetric
    model. (numerics.md "Bracket before you root-find"; rows 17, 33, 65, 81, 156, 169,
    182)
13. An extremum is a root of an analytic derivative (pole multiplied out), or a bounded
    maximisation refined from a scan of the whole window; an argmax on a domain end
    means "at or beyond the end". Take limits in closed form, not at the grid end; at a
    finite point a limit is a rate, a theory bound. (rows 22, 96, 105, 131, 147, 173,
    189)
14. Orders come from successive differences of at least three levels
    (`pymrm.checks.observed_orders`), never against the finest grid or the exact value
    on a one-axis ladder. Report each case as a theory range with the ladder printed
    (minimum and maximum where they oscillate), plus the round-off companion bound.
    (claims-design 3; rows 9, 19, 43, 75, 95, 132, 134, 163)
15. Split an error before naming its cause: time against space (exact time by the
    eigendecomposition of the constant operator, or the other axis cancelled;
    observations on step times `(i + 1) * dt`), approximation against discretisation,
    each axis of a 2-D grid on its own ladder. (rows 48, 74, 82, 83, 89, 122, 126,
    133, 181)
16. Name a cause only after removing it and watching the effect go, one ingredient at
    a time; a cause named for a gap is a claim, measured with the term added as a
    parameter. Never attribute an effect to one unmeasured cause. (rows 37, 79, 87,
    123, 166, 188)
17. Measure the reference first: scaled integrator `atol` and a second tolerance, a
    root-found reference checked at its root, series truncation at the largest
    parameter, `quad` tolerances set, cancelling closed forms against `mpmath`; evaluate
    it at cell centres, not through `np.interp`. Richardson removes discretisation, not
    an estimator's bias. (rows 36, 63, 77, 86, 98, 117, 125, 135, 137)
18. A statement over a range or a set of cases is evaluated at every case and point it
    covers, on the page's grids; root-find where a comparison flips. State conditions
    with the range ("for one light gas throughout"). (rows 94, 142, 143, 145, 170, 176,
    190)
19. A mean bias needs a standard error before it is called anything, and a reference
    group that is what it is said to be; split data where the abscissa certifies it.
    (rows 176, 177, 179)
20. State what two routes share (end conditions, interpolant, factorisation, source
    function); routes sharing an interpolant agree to the coarser root finder. Measure
    a second `permc_spec` before claiming solver independence. (attack-list.md
    section 2; rows 59, 150, 153)

## pymrm modelling (the plugin owns the rule)

21. A first-order convection-only outlet is `{"outflow": True}`; an outlet bc of `None`
    is a Dirichlet zero. Read outlet, wall and boundary fluxes with
    `compute_boundary_values` under the same bc, never from the last cell. A van Leer
    profile with an outlet gradient, under outflow or zero gradient, has a first-order
    layer of 2 to 3 cells at the outlet while the last cell is a second-order outlet
    value and the interior stays second order: measure an order on the outlet value or
    at a fixed x, never as a maximum over all cells; an extrapolated outlet face
    removes the layer. (P6, P8; rows 32, 35, 87, 121, 141, 148, 149, 166, 187; pymrm
    #41)
22. Field shapes `(n, 1)`, `NumJac((n, 1))`, spatial axes first. A rate whose
    derivative is unbounded or varies fast gets its exact Jacobian through
    `construct_coefficient_matrix`; a state-dependent velocity is convected as a molar
    flux with a constant unit-velocity operator. (AGENTS.md "pymrm conventions you must
    follow"; P3, P9; rows 101, 110, 184)
23. Check every solve: `newton` `success` (`NewtonBreakdown` raises instead),
    `solve_ivp` `status`; `tol` above the round-off floor with `stagnation=True`;
    bounds rather than a custom step cap, tested apart from damping. (P7;
    style-guide.md section 11; rows 13, 26, 47, 57, 104, 139)
24. Newton converges to a nearby root, stable or not: locate a threshold by following
    the branch with warm starts, and say which branch a state is on. Choose
    `ConvFluxTVD` with `newton` or deferred correction by measuring both. (P10; rows 8,
    18, 71, 76, 110, 120)
25. Test a discrete identity on the scheme the page uses: exact cell volumes, a
    pure-Neumann source with zero mean on the grid, face gradients from
    `compute_boundary_values`. Its gap is the solve's residual, not a second route.
    (rows 31, 62, 69, 93, 183)

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
28. A failure mode or trap stated in prose needs the cell that produces it and a
    claim, or it goes. (rows 53, 60)

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
   printed claims (rules 26 to 32).
