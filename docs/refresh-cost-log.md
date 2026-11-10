# Refresh cost log

Tokens and review findings per refreshed page, to judge the lean refresh cycle
(`.claude/skills/gallery-refresh/SKILL.md`, "Trial of the lean cycle"). One row per
page, added by the coordinator after the commit.

How to fill a row:

- tokens: `python scripts/agent_tokens.py <session dir> --page <ID>`, where the
  session directory is `~/.claude-tue/projects/-home-eajfpeters-Code-pymrm-suite-pymrm-gallery/<session-id>`
  of the coordinator that ran the page. Prompt tokens are uncached input plus cache
  reads plus cache writes, summed over the agent's calls; "last ctx" is the context
  of the agent's last call.
- must-fix: the number of MUST-FIX items in the review's numbered list (brief B).
- classes: the review's must-fix items by kind, in the words of SKILL.md "What the
  review keeps finding", plus any new kind spelled out.

## Lean cycle

| page | date | lane | port M | numerics or light M | refresh last ctx k | review M | fix round M | total M | must-fix | classes | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A4.2 | 2026-11-10 | model-solving | 0.9 | 15.8 | 279 | 1.5 | 4.6 | 22.8 | 5 | roots of step values published without the step error (extrapolate in dt); a moved published value without a Correction; a Correction quoting code instead of the printed line; a root helper that raises and crashes the study; a learnings row with the wrong lesson | lean cycle with call budget and break-row rule. Two Correction quotes are notebook outputs (checked by hand against origin/main); accept_page lint cannot see them |
| A4.1 | 2026-11-10 | model-solving | 2.4 | 20.9 | 337 | 1.9 | 0.9 | 26.1 | 4 | a break row that crashes the whole study instead of failing its claim; a moved published number not listed; an unmeasured cause; a learnings row saying more than measured | lean cycle with call budget and break-row rule: 89 calls |
| A3.12 | 2026-11-09 | model-solving | 1.7 | 14.8 | 282 | 2.0 | 1.7 | 20.2 | 4 | an effect asserted without measurement (an offset that would "move the zero just as much"; relaxation that "would not change that"); a stopping tolerance below round-off; a crossing claimed beyond the range of its input table | lean cycle with the call budget and the break-row rule (after the trial): 79 calls, 34 rows defined and printed |
| A3.8 | 2026-11-09 | model-solving | 1.9 | 17.8 | 289 | 1.7 | 2.5 | 24.0 | 4 | a Correction of a statement that was true; a boundary stated at the wrong parameter (a read-off where a root was needed); an unrefreshed models.yaml note with untagged numbers and a history reference; (plus the Correction rewording that followed) | lean cycle, trial page 3. Numerics brief with the call budget and the break-row count rule: 91 calls, 38 rows defined and printed |
| A3.7 | 2026-11-09 | model-solving | 1.8 | 17.2 | 324 | 1.5 | 0.5 | 20.9 | 4 | break rows that silently record nothing (an exception swallowed); a statement resting on such a row; a share stated without its measure; learnings rows with missing columns | lean cycle, trial page 2. The numerics brief added a call budget ("well under 100"): 78 calls against 127 on A3.5. Runtime 6 s -> 127 s (nine transient break rows) |
| A3.5 | 2026-11-09 | model-solving | 1.6 | 31.2 | 369 | 1.6 | 0.0 | 34.4 | 2 | an old number outside a Correction quote; overstated "What pymrm adds" (closed-form roots credited to the solve) | lean cycle, trial page 1. Fix round done by the coordinator (README and meta text). The numerics agent made 127 calls and ended at 369k context, about the old single-agent size; the port was cheap |
| A3.6 | 2026-11-09 | model-solving | - | 14.9 | 291 | 1.7 | 0.0 | 16.7 | 1 | a published wrong order statement not corrected | not the lean cycle: one Opus refresh agent on a hand-trimmed brief (checklist plus one worked page, the rest by grep, Edit and ranged reads); fix by the coordinator |
| A3.4 | 2026-11-09 | model-solving | - | 20.1 | 341 | 1.6 | 1.9 | 23.6 | 3 | a wrong published statement in a data sidecar not corrected; a verdict stated on a side its margin does not support; a computed number without a claim marker | not the lean cycle: one Opus refresh agent on a hand-trimmed brief |

Trial outcome (2026-11-10, the maintainer's decision relayed by pymrm-gallery): back to one Opus refresh agent per model-solving page, with the hand-trimmed brief of A3.6 and A3.4 plus the call budget; page_facts.py, accept_page.py and the break-row count rule stay. The split saved little: the port is cheap (1 to 2 M) but the numerics agent costs about what one trimmed agent does (15 to 31 M against 15 to 20 M), and must-fix counts stayed within the baseline.

## Baseline: before the lean cycle

Measured from the transcripts of the Team-account coordinator session `87e352fb`
(2026-10-31 to 2026-11-09) with the same script, one row per page whose agents named
it. Until 2026-11-03 the fix round went back to the refresh agent (its tokens are in
"refresh M" and the fix-round column is 0). "n/a" marks an agent whose transcript
ends without a usage record.

| page | refresh started | refresh calls | refresh M | refresh last ctx k | review M | fix round M | total M | must-fix |
|---|---|---|---|---|---|---|---|---|
| A3.15 | 2026-10-31 | 88 | 24.7 | 400 | 4.6 | 0.0 | 29.3 | 0 |
| A3.1 | 2026-10-31 | 104 | 32.3 | 443 | 4.0 | 0.0 | 36.3 | 5 |
| A3.2 | 2026-10-31 | 117 | 39.0 | 466 | 3.2 | 0.0 | 42.2 | 0 |
| A3.3 | 2026-10-31 | 95 | 27.0 | 402 | 2.6 | 0.0 | 29.6 | 2 |
| B1.2 | 2026-11-01 | 90 | 19.9 | 324 | 4.7 | 0.0 | 24.7 | 2 |
| B1.3 | 2026-11-01 | 119 | 38.1 | 484 | 6.4 | 0.0 | 44.5 | 4 |
| B1.4 | 2026-11-01 | 90 | 20.1 | 371 | 4.8 | 9.6 | 34.5 | 0 |
| B1.6 | 2026-11-01 | 111 | 32.8 | 443 | 5.6 | 2.6 | 41.0 | 2 |
| B1.7 | 2026-11-01 | 112 | 29.0 | 432 | 8.3 | 1.9 | 39.2 | 5 |
| B2.2 | 2026-11-02 | 92 | 25.3 | 396 | 1.6 | 0.4 | 27.4 | 1 |
| B2.3 | 2026-11-02 | 161 | 43.7 | 54 | 3.0 | 3.8 | 50.5 | 0 |
| B3.1 | 2026-11-02 | 38 | 4.9 | 193 | 0.9 | 1.4 | 7.2 | 2 |
| B3.2 | 2026-11-02 | 72 | 13.1 | 274 | 1.5 | 1.7 | 16.3 | 4 |
| B3.3 | 2026-11-02 | 76 | 21.5 | 412 | 2.0 | 0.8 | 24.3 | 2 |
| C2.1 | 2026-11-02 | 63 | 11.2 | 267 | 1.3 | 1.1 | 13.6 | 1 |
| A1.1 | 2026-11-02 | 48 | 7.3 | 232 | 1.1 | 1.6 | 10.0 | 2 |
| E1.1 | 2026-11-02 | 96 | 20.5 | 321 | 1.4 | 2.7 | 24.6 | 3 |
| E1.2 | 2026-11-02 | 67 | 12.1 | 264 | 1.0 | 0.7 | 13.8 | 1 |
| A1.5 | 2026-11-02 | 87 | 22.5 | 392 | 1.2 | 3.7 | 27.4 | 4 |
| E2.1 | 2026-11-03 | 68 | 13.9 | 291 | 1.2 | 2.2 | 17.3 | 2 |
| F1.3 | 2026-11-03 | 44 | 7.8 | 284 | 1.1 | 2.8 | 11.7 | 5 |
| F1.4 | 2026-11-03 | 36 | 4.1 | 170 | 0.9 | 1.4 | 6.4 | 6 |
| J1.5 | 2026-11-03 | 44 | 7.1 | 247 | 0.9 | 1.5 | 9.6 | 8 |
| F2.3 | 2026-11-03 | 94 | 23.7 | 391 | 1.4 | 3.1 | 28.2 | 6 |
| A1.6 | 2026-11-03 | 62 | 13.5 | 348 | 1.8 | 1.4 | 16.8 | 8 |
| C1.3 | 2026-11-03 | 90 | 23.9 | 382 | 1.7 | 1.0 | 26.6 | 3 |
| J3.4 | 2026-11-03 | 122 | 36.0 | 450 | 1.7 | 1.4 | 39.1 | 1 |
| B2.1 | 2026-11-04 | 97 | 24.8 | 379 | 1.4 | 1.7 | 27.9 | 2 |
| J3.5 | 2026-11-08 | 91 | 25.1 | 395 | 1.6 | 0.6 | 27.3 | 2 |
| A1.7 | 2026-11-09 | 83 | 23.2 | 439 | 1.3 | 2.6 | 27.0 | 4 |
| A1.8 | 2026-11-09 | 120 | 37.7 | 498 | 1.3 | 2.4 | 41.4 | 4 |
| A2.6 | 2026-11-09 | 113 | 37.1 | 500 | 1.7 | 2.2 | 41.0 | 2 |

32 pages: total per page median 27.3 M (range 6.4 to 50.5); must-fix median 2.0, range 0 to 8.
