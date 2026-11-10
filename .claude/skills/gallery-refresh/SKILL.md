---
name: gallery-refresh
description: Run the pymrm-gallery refresh cycle for one published page (refresh, independent review, fix round, spot-check, commit, PR, merge). Use when refreshing a page under docs/refresh-plan.md, briefing a refresh, review or fix agent for a gallery page, or committing and merging refreshed pages.
---

# Gallery refresh cycle

This skill is for the coordinating session. The page contract stays in `AGENTS.md`
and `docs/claims-design.md`; the evidence for every rule stays in
`docs/refresh-learnings.md`. What is here is the operational part: which agent runs
when, what each brief contains, the checks that decide a commit, and the
recurring defects condensed into one checklist.

## Files

| file | what it is |
| --- | --- |
| `references/refresh-checklist.md` | the rules every refreshed page must meet, distilled from the learnings rows; the refresh agent applies it up front, the reviewer checks against it |
| `references/briefs.md` | templates: A refresh (model-solving pages, Opus), A-light (light pages, Sonnet), B review, C fix round; the working-style block every brief carries |
| `references/commit-and-merge.md` | staging with shared files, spot-checks, commit and PR text, CI and merge rules, known checker behaviour |
| `.claude/agents/gallery-reviewer.md` | the review agent (Opus, xhigh, advice only, time-boxed) |
| `.claude/agents/gallery-fixer.md` | the fix-round agent (Sonnet, numbered list plus standard acceptance) |
| `scripts/page_facts.py` | the facts a brief needs about one page, from the old version (solver lines, metrics, sibling pages, learnings rows) |
| `scripts/accept_page.py` | every acceptance step in one command, about 30 lines of output, full log to a file |
| `scripts/agent_tokens.py` | tokens per agent and per page, from the session transcripts |
| `docs/refresh-cost-log.md` | tokens and must-fix count per page |

## Why the cycle is lean

Cost is roughly calls times context: whatever an agent reads early is paid again on
every later call. Measured on the Team account (31 Oct to 9 Nov), refresh agents made
about 90 % of the tokens: 80 to 120 calls ending at 400 to 500k context, 25 to 40 M
prompt tokens per page, while review and fix rounds took 1 to 3 M each. The context
grew from shell reads of whole files (42 %), script and notebook output (24 %), the
read-first list (16 %) and whole-file writes (10 %). Hence: page facts instead of
discovery, the checklist and the pitfalls in full and everything else by grep, one
worked page, Edit instead of Write, one quiet check command, a call budget, and few
coordinator check-ins. The review and the fix round are unchanged.

**Split refresh agents did not pay (trial, 9 Nov).** A Sonnet port followed by an
Opus numerics agent cost 34.4, 20.9 and 24.0 M on A3.5, A3.7 and A3.8 (must-fix 2, 4,
4), against 16.7 and 23.6 M for one Opus agent with the same reading discipline on
A3.6 and A3.4 (must-fix 1 and 3). The port was cheap (1.6 to 1.9 M) and exact, but
the numerics agent still read and rewrote the whole page (on A3.5, 110k characters of
page reads and 124k of cell-by-cell writes). So the cycle runs one refresh agent per
page. The trial also found one new must-fix class, break rows that silently recorded
nothing; `claims.define_rows` and the rows step of accept_page.py now catch it.

## The cycle, per page

1. **Page facts and lane.** `<VENV>/bin/python scripts/page_facts.py <ID> > <SCRATCH>/<id>_facts.md`.
   Read it (not build_page.py) and decide the lane: model-solving (a pymrm solve, a
   model.py) or light (closed forms, no pymrm operator). Pick ONE worked page of the
   same lane from its list (prefer a recent one with the same structure code) and the
   matching plugin exemplar, and write the lane B notes from the solver lines.
2. **Refresh**, one agent, brief written to the scratchpad: **A** (`implement`,
   `model: opus`) for a model-solving page, **A-light** (`implement`, Sonnet) for a
   light one. The agent writes a progress file after every checkpoint, so a
   usage-limit stop is recoverable.
3. **Review** with `gallery-reviewer` from brief B. Restate the implementer's
   findings as facts to verify and put the riskiest judgement call first (a
   reversed conclusion, a changed reading of an equation, a cited page).
4. **Fix round** with `gallery-fixer` from brief C: the reviewer's items, each with
   file:line and the exact replacement wording, plus anything the coordinator adds.
   Drop review items that are themselves wrong (a change of definition or a fair
   rounding proposed as a correction) and say why in the commit message.
5. **Spot-check and commit** (`references/commit-and-merge.md`): read the README
   Corrections yourself; fix small wording defects directly, then the acceptance
   subset in `references/commit-and-merge.md`.
6. **Log the cost**: one row in `docs/refresh-cost-log.md` from
   `scripts/agent_tokens.py <session dir> --page <ID>` and the review's must-fix count.
7. **Batch pages into one PR** so one CI run (about 2.5 h) covers several pages.
   Do not push to a branch whose PR CI is still running unless the new commit
   should be in that PR; merge when CI is green (standing authorisation for this
   sweep).

## Coordinator

- **Read little yourself.** The page facts file, the agents' reports and the
  acceptance summaries are what you need; open a page file only for the spot-check.
  Never paste a long output into the conversation when a file path will do.
- **Check-ins.** The harness notifies you when an agent finishes; do not poll. Ask
  `hub-tue` for the usage figure before starting each agent. Set a one-shot
  CronCreate check-in only for an agent expected to run longer than about 45
  minutes, and only one per agent: each check-in is a full-context call.

## Pacing (Team account)

- Work up to about 90 % of the 5-hour window; hub-tue warns at 80 % and sends the
  pause at 90 %. All development stops at 96 % weekly until the weekly reset
  (hub-tue sends the stop).
- One agent at a time.
- Measured costs before the reading discipline (points of the 5-hour window): Opus
  model-solving refresh 18 to 31 (long pages more), Opus review 5 to 10, Sonnet fix
  round 2 to 6, Sonnet light refresh 6 to 15. A full window is about 16 weekly
  points. The cost log has the tokens per page since 9 Nov.

## What the review keeps finding

Every refresh so far has needed a review: each found 0 to 8 must-fix items (34
reviews, 31 Oct to 9 Nov), almost all of the same classes. The checklist exists to
move them into the refresh itself: numbers without claim markers, corrections that
misquote the old wording or are a change of definition, grid read-offs where a root
is needed, orders without a round-off companion, effects attributed to one
unmeasured cause, overstated "What pymrm adds", learnings rows that say more than
was measured. Keep the review step anyway; it is the cheapest point of the cycle at
which a wrong published statement is caught.

## Learnings rows

A finding that would change what an agent writes or checks when building a pymrm
model goes to `docs/refresh-learnings.md` as a row with evidence. Gallery-only
rules belong in `references/refresh-checklist.md` as well; generic pymrm rules
(numerics, pitfalls, verifier attacks, library bugs) are flushed to the pymrm repo
(targets listed at the top of the learnings file). Mark a row `filed` or `merged`
with the pymrm PR number when that happens.
