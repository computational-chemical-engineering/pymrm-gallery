# Learnings from refreshing the published pages

One row per finding made while refreshing a page under
[`refresh-plan.md`](refresh-plan.md). A finding belongs here when it would change
what an agent writes or checks when building a pymrm model, and it carries evidence:
the page, and the numbers that showed it.

Targets in the pymrm repo (`plugins/pymrm/` unless stated):
`pitfall` (`skills/conventions/references/pitfalls.md`, with its test in
`pymrm/test/test_pitfalls.py`), `api-map`, `style-guide`, `numerics`, `exemplar`,
`eval` (a new case under `evals/`, usually built from the page), `verifier`
(`skills/verify-model/references/attack-list.md`), `library` (a pymrm bug or a missing
API; file it against `pymrm/src/pymrm/`).

Status: `open` (recorded), `filed` (in a pymrm PR, give the number), `merged`, or
`dropped` (with the reason).

Flush open rows to pymrm as one PR per batch of about ten refreshed pages.

| # | Date | Page | Finding | Evidence | Target | Status |
|---|---|---|---|---|---|---|
| 1 | 2026-10-29 | `J4.6`, `J5.3`, `F3.2` | Asserting a number at the round-off level in prose checks fails on another machine; CI was red from at least 2026-09-14 on this alone. State a bound or an order of magnitude, or compare at a tolerance. | CI runs 2026-09-14 to 2026-10-26: `9.7e-16` stated against `4.8e-16` computed (`J4.6`); runner-dependent prose drift on `F3.2` and `J5.3` | `style-guide`, `verifier` | open |
| 2 | 2026-10-29 | `J4.1`, `J5.3` | Deciding by a value's size whether it is a metric (reported only if below `ABS_FLOOR = 1e-12`) makes the metric set machine-dependent: an optimiser residual lands on different sides of the floor, and the page's own consistency checks fail. Classify by what the quantity is (zero by construction, tolerance-set, physical), never by its value. | `J4.1` `sep_moser` fails on every GitHub runner and under both alternative kernels; `J5.3` `mw_reconstruction_x_maxdiff` falls below the floor under `OPENBLAS_CORETYPE=Sandybridge` | `style-guide`, `verifier` | open |
| 3 | 2026-10-29 | all | Cross-machine reproducibility can be tested on one machine: `OPENBLAS_CORETYPE=Sandybridge` or `Nehalem` changes BLAS round-off the way a different runner does. It reproduced the CI failures of all four pages tried so far (`F3.2`, `J4.1`, `J4.6`, `J5.3`). A verifier can use it to tell a robust check from one that passes by luck. | `np.linalg.solve` on a 200x200 system differs in the 14th digit between kernels; all four failing pages fail locally under both alternatives | `verifier`, `numerics` | open |
| 4 | 2026-10-29 | `D1.1` | Near a threshold where several steady states exist, which branch Newton reaches from a deterministic cold start depends on round-off, so a page that locates the threshold on a different machine can land on the ignited branch below it. | `D1.1` rung 2 at `n_z = 1200` passes on some runners and raises `NotACrossing` on others (993.75 C below `p0*`, 415.02 C above) | `numerics`, `pitfall` | open |
| 5 | 2026-10-29 | CI | A regression check that runs only when every page passes protects nothing while any page fails: `check_agreement.py` was skipped on every CI run from at least 2026-09-14 to 2026-10-26. Run it with `if: always()`. | `execute.yml` before this branch | (gallery CI only) | fixed on `refresh/step0-ci-green` |
