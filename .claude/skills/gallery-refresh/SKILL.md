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
| `references/briefs.md` | templates: A refresh (model-solving lane), A-light (light lane), B review, C fix round |
| `references/commit-and-merge.md` | staging with shared files, spot-checks, commit and PR text, CI and merge rules, known checker behaviour |
| `.claude/agents/gallery-reviewer.md` | the review agent (Opus, xhigh, advice only, time-boxed) |
| `.claude/agents/gallery-fixer.md` | the fix-round agent (Sonnet, numbered list plus standard acceptance) |

## The cycle, per page

1. **Pick the page and the lane.** Model-solving pages (a pymrm solve, a model.py)
   get the refresh on Opus (`implement` with `model: opus`) early in a usage window;
   light pages (closed forms, no pymrm operator) get it on Sonnet (`implement`),
   which also fits near the top of a window. Read the old `build_page.py` briefly
   for the solver calls and root finds, and write the page-specific lane B notes.
2. **Refresh** from brief A or A-light, written to the scratchpad. The agent writes a
   progress file after every checkpoint, so a usage-limit stop is recoverable.
3. **Review** with `gallery-reviewer` from brief B. Restate the implementer's
   findings as facts to verify and put the riskiest judgement call first (a
   reversed conclusion, a changed reading of an equation, a cited page).
4. **Fix round** with `gallery-fixer` from brief C: the reviewer's items, each with
   file:line and the exact replacement wording, plus anything the coordinator adds.
   Drop review items that are themselves wrong (a change of definition or a fair
   rounding proposed as a correction) and say why in the commit message.
5. **Spot-check and commit** (`references/commit-and-merge.md`): read the README
   Corrections yourself; fix small wording defects directly, then re-run
   `jupytext --sync` and `scripts/check_claims.py`.
6. **Batch pages into one PR** so one CI run (about 2.5 h) covers several pages.
   Do not push to a branch whose PR CI is still running unless the new commit
   should be in that PR; merge when CI is green (standing authorisation for this
   sweep).

## Pacing (Team account)

- Ask `hub-tue` for the 5-hour figure before each agent. Work up to about 90 % of
  the window; hub-tue warns at 80 % and sends the pause at 90 %. All development
  stops at 96 % weekly until the weekly reset (hub-tue sends the stop).
- One agent at a time. For an agent expected to run longer than about 20 minutes,
  set a one-shot CronCreate check-in that asks hub-tue for the figure and sends
  "pause" to the agent near the limit.
- Measured costs (points of the 5-hour window): Opus model-solving refresh 18 to 31
  (long pages more), Opus review 5 to 10, Sonnet fix round 2 to 6, Sonnet light
  refresh 6 to 15. A full window is about 16 weekly points.

## What the review keeps finding

Every refresh so far has needed a review: each found 2 to 8 must-fix items, almost
all of the same classes. The checklist exists to move them into the refresh itself:
numbers without claim markers, corrections that misquote the old wording or are a
change of definition, grid read-offs where a root is needed, orders without a
round-off companion, effects attributed to one unmeasured cause, overstated "What
pymrm adds", learnings rows that say more than was measured. Keep the review step
anyway; it is the cheapest point of the cycle at which a wrong published
statement is caught.

## Learnings rows

A finding that would change what an agent writes or checks when building a pymrm
model goes to `docs/refresh-learnings.md` as a row with evidence. Gallery-only
rules belong in `references/refresh-checklist.md` as well; generic pymrm rules
(numerics, pitfalls, verifier attacks, library bugs) are flushed to the pymrm repo
(targets listed at the top of the learnings file). Mark a row `filed` or `merged`
with the pymrm PR number when that happens.
