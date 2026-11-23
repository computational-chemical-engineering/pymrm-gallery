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
| J1.1 | 2026-11-17 | model-solving | - | 19.5 | 309 | 2.7 | 3.3 | 25.5 | 1 | README counts (13 of 14, 17 of 20, ...) left as an audit list, two unchecked, instead of exact claims | one Opus refresh agent, F3.2 as worked page (98 calls); format-2 baseline with a quoted map (33 to 50 literals), no loosening entries; one Correction (newton stopped at maxfev on two Freundlich fits, unchecked; no published value moved) |
| J1.4 | 2026-11-17 | model-solving | - | 22.1 | 303 | 1.6 | 0 | 23.7 | 1 | Correction that was a change of definition (an order measured against the closed form, true at its own wording) | one Opus refresh agent, J1.3 as worked page (102 calls); format-1 baseline, 40 claims, 41 says literals, 28 break rows; orders now by successive differences; review fixes applied by the coordinator |
| J1.3 | 2026-11-17 | model-solving | - | 28.1 | 383 | 1.6 | 0 | 29.7 | 1 | Correction that was not one (an old sentence true at its own wording, counted as wrong) | one Opus refresh agent, J1.1 as worked page (115 calls); format-1 baseline, 87 claims, 70 says literals, 60 break rows; review fixes applied by the coordinator (Correction removed, Table I labels computed, a round-off print given its bound) |
| F3.5 | 2026-11-17 | model-solving | - | 19.5 | 316 | 2.2 | 0.9 | 22.6 | 2 | the claims.says list left empty (prose precision not recorded); a history phrase the dashes scan missed | one Opus refresh agent, F3.2 as worked page (90 calls); format-1 baseline, so no loosening entries; two Corrections (a break-row claim and the charge-flux reason) |
| F3.2 | 2026-11-16 | model-solving | - | 33.9 | 342 | 1.7 | 2.2 | 37.8 | 4 | a published asymptote that holds only for large E_i (Levenspiel's correction against the fixed-q limit, 0.9727 on the worst curve); a loosened.yaml reason calling a physical deficit round-off; a bracketed root left after rule 12; a validation line claiming a uniform-grid order no claim covers | one Opus refresh agent, A4.3 as worked page (149 calls, 2309-line old page with 19 break studies); the dropped claims.says list then failed CI's prose-literal loosening check (38 lines), parked for a Decision Desk card |
| C1.1 | 2026-11-16 | light (statistics, one small reactor) | - | 19.0 | 313 | 1.7 | 2.9 | 23.6 | 3 | a Correction that also retracted a true statement a nested F test would have confirmed (paired SE used where the comparison was nested); a duplicate-copy verdict decided below the page's own half-unit criterion | one Opus refresh agent (86 calls); review and fix in a fresh context after a /clear, at a fresh 5-hour window |
| A4.7 | 2026-11-16 | model-solving | - | 21.9 | 351 | 2.2 | 2.3 | 26.4 | 2 | production values compared with an exact series without their own discretisation error; a claim bound with no margin (0.99999 against 1) and a learnings row quoting a local value as the maximum | one Opus refresh agent, A4.6 as worked page (89 calls); no Corrections (the old values were true at their stated digits); check-in crons every 25 min kept the window watched |
| A4.6 | 2026-11-16 | model-solving | - | 20.8 | 332 | 1.6 | 1.8 | 24.2 | 5 | a Correction resting on an unconverged quadrature (single-knob doublings understated its error); a false interior-extremum finding; data sidecars left stating what the page revises; a learnings row without the measured lesson; a trend stated in the wrong direction | one Opus refresh agent, A4.4 as worked page (88 calls); started at 62 % of the window, ran past 90 % unwatched (hub-tue sees no readings while this session idles) |
| A4.5 | 2026-11-16 | model-solving | - | 15.3 | 284 | 1.6 | 1.7 | 18.6 | 4 | a false identity stated as the reason two routes agree; two Corrections for statements true by their own count or within input rounding (rule 26); a pymrm-versus-quadrature comparison not stated (rule 29) | one Opus refresh agent, A3.8 as worked page (76 calls); coordinator reworded the N2 rounding sentence |
| A4.4 | 2026-11-16 | model-solving | - | 18.7 | 321 | 1.7 | 0.4 | 20.8 | 1 | a Correction for a published figure that was true over what it covered (rule 26) | one Opus refresh agent, hand-trimmed brief with the sibling A4.3 as worked page (83 calls); review optional item 3 (a stronger round-off companion) not applied |
| A4.3 | 2026-11-16 | model-solving | - | 29.7 | 399 | 2.4 | 1.5 | 33.6 | 3 | a location called extrapolated but maximised on the 400-cell grid; a Correction for a check gap that was an upper bound, not wrong (rule 26); a learnings row stating a tighter constancy than measured | one Opus refresh agent, hand-trimmed brief plus call budget (104 calls); item 7 (value claim at rtol 1e-2) not applied: the kernel spread needs a tolerance above 1e-4, so the band claim stays |
| C1.2 | 2026-11-10 | light | - | 9.3 | 269 | 1.2 | 0.2 | 10.8 | 2 | a published verdict the page's own measurement contradicts, kept without a Correction; a noiseless misfit compared with a noisy scatter | one Sonnet agent (single-agent cycle after the trial), call budget 80: 52 calls |
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
| J4.4 | 2026-11-19 | 76 | 14.4 | 299 | 3.2 | 0.9 | 18.5 | 1 |
| J3.1 | 2026-11-19 | 131 | 37.0 | 419 | 2.3 | 0.9 | 40.2 | 2 |
| J4.2 | 2026-11-19 | 102 | 31.6 | 445 | 2.3 | 0.9 | 34.8 | 1 |
| J1.10 | 2026-11-19 | 137 | 35.4 | 370 | 2.4 | 2.4 | 40.2 | 4 |
| J4.8 | 2026-11-19 | 87 | 19.5 | 323 | 1.9 | 2.9 | 24.3 | 3 |
| F3.1 | 2026-11-19 | 61 | 9.2 | 210 | 1.4 | 3.6 | 14.2 | 5 |
| J4.1 | 2026-11-19 | 259 | 43.8 | 283 | 2.5 | 2.5 | 48.8 | 4 |
| J4.6 | 2026-11-19 | 144 | 37.4 | 361 | 2.4 | 2.7 | 42.5 | 4 |
| A4.9 | 2026-11-19 | 87 | 16.1 | 254 | 1.8 | 4.5 | 22.4 | 5 |

2026-11-19 sprint (personal account, pages in parallel, one worktree each; fix rounds on Opus where they changed code). J4.1's refresh took three agents on the call budget. Second, independent reviews of J1.10 and J4.8 found 2 and 1 further must-fix items (not in the must-fix column).

32 pages: total per page median 27.3 M (range 6.4 to 50.5); must-fix median 2.0, range 0 to 8.

## New pages (build cycle with a plan step, workflow C1-C4)

Tokens as prompt M (input plus cache) per agent from `scripts/agent_tokens.py`. Plan: `think` writes
PLAN.md, EQUATIONS.md and a check-value script; build: `implement` on Opus; review: `gallery-reviewer`
(Opus); fix: `gallery-fixer` (Sonnet). Reviewer test C5 (queue rule 8): the Fable review of the same
pre-fix commit runs on the personal account and is noted here, not added to the total.

| page | date | plan M | build M | build last ctx k | review M | fix round M | total M | must-fix | notes |
|---|---|---|---|---|---|---|---|---|---|
| F1.8 | 2026-11-23 | 4.3 | 9.8 | 232 | 1.6 | 6.3 | 22.0 | 5 distinct (Opus 4, Fable 3, all verified) | Team account. Fable review 0.662 M units against Opus about 0.42 M (C5 page 1: Fable -1 verified must-fix at about 1.6x; Fable's unique item, the Fig 6-a.1 difference unmatched in uG, moved the headline from 0.300 to 0.222). Coordinator: check values rewritten to 12 digits (6 digits failed rtol 1e-6), Styler display replaced (memory address in repr), an unverifiable handbook citation replaced by plain wording. Prose-diff after the fix: 0 new unsupported sentences. |
