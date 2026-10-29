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
