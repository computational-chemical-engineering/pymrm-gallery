# Prompt for the next session

Copy everything in the block below as the first message of a fresh session started in
`~/Code/pymrm_suite/pymrm-gallery`.

---

```
Continue the refresh sweep of the pymrm gallery (repo ~/Code/pymrm_suite/pymrm-gallery).

Read, in this order: AGENTS.md; docs/refresh-plan.md, all of it, section 8 first (the
state on 2026-10-31, the next pages, the open items); docs/playbooks/refresh-page-brief.md
(the refresh brief, the review brief and the commit checklist); docs/claims-design.md
sections 6-7; the last 15 rows of docs/refresh-learnings.md; git log -10.

Then take the next page of section 8 through the cycle: implement agent (Opus) with
brief A; think review with brief B when a published statement or a metric beyond its
tolerance changes; the findings back to the same implementer; spot-checks; one commit
per page on a section branch; a PR every few pages; merge a case PR when CI is green.
Ask the maintainer before anything else outward-facing. Check first whether PR #11 has
been merged; if it has, confirm the site deployed.

Build a venv from requirements-lock.txt for the checks. Never git add -A; check
git diff --cached before every commit while an agent shares the checkout.
```

---

## Why each part is there

**Section 8 first** because it is the only place that records which pages are done and
what is waiting on the maintainer. The rest of `refresh-plan.md` holds the decisions
(lanes, order, review scope) that the brief assumes.

**The review is not optional where statements change.** On every page so far whose
published statements changed (`A2.1`, `A2.3`, `A2.5`, the reactor section), the
independent `think` review found must-fix defects: a correction resting on a
hand-picked subset of cases, a ratio of two values from different grids, the wrong
cause named for a broken identity.

**One page per commit, an index checked before committing,** because agents share the
checkout and its git index: on 2026-10-31 a commit swept in another agent's staged
deletion and needed a restore commit.

The build queue of 2026-08 (coordinator playbook, dashboards, one case at a time) is
paused while the sweep runs; `docs/handoff.md` and `docs/playbooks/coordinator.md`
describe it.
