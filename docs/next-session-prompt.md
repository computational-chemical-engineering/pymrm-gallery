# Prompt for the next session

Copy everything in the block below as the first message of a fresh session.

---

```
Continue the pymrm-gallery: work the 266-case catalogue with parallel agents.

Repo: ~/Code/pymrm_suite/pymrm-gallery   Python: ~/Code/pymrm_suite/.venv
Read first: docs/handoff.md (operating procedure), AGENTS.md,
docs/agent-brief.md, docs/pdf-findings.md.

State: 16 published, 3 need my judgement, 10 need papers, 232 unclaimed.
Every parked case carries a resume: block — use it, don't re-derive.

THE RULE: never block on me. A case needing input is parked with its
resume: block and you immediately start the next. Keep ~5 builders live
continuously. I read the dashboards when I choose — don't ping me per
case, don't idle waiting. I may be away from a screen for long stretches;
assume no answers are coming and keep working.

The loop, per docs/handoff.md: dispatch builder → on ready, dispatch an
adversarial verifier → integrate (splice models_entry.yaml, check_metadata,
run_pages --changed, check_agreement) → commit and push → refill to 5.
Blocked cases get parked and the dashboards regenerated and republished.

Dispatch policy: prefer cases whose paper is on disk AND whose validation
looks like a table, appendix or stated result rather than a figure — only
the figure route needs me, so only it stalls a case.

Keep the adversarial verifier on every ready page. Quality over throughput:
a missing page costs nothing, a page with a fabricated number costs the
credibility of all the others.

Publish everything that's finished.

Two decisions from me you're waiting on — carry them forward, ask once,
and don't block on either:
(1) purging three copyrighted overlays from git history (needs a
    force-push to a public repo);
(2) the Elsevier API key at ~/.config/elsevier/apikey is blocked by the
    Bash permission classifier — it's the highest-leverage unblock for
    section A.

Dashboards (republish these same URLs, don't mint new ones):
  papers    https://claude.ai/code/artifact/99b52225-ce54-4487-99ec-3a420f2ac4ad
  decisions https://claude.ai/code/artifact/237909a0-3432-4b82-8cfe-15b05019d036
```

---

## Why each part is there

**"Never block on me"** is the correction that prompted this file. Earlier
sessions stopped and reported after each blocked case. Throughput is set by how
rarely a case needs the maintainer, so a blocked case is parked — with a
`resume:` block recording what is already established, what each plausible
answer changes, and which files to touch — and the next case starts at once.

**The dispatch policy** is the main throughput lever, and it is about avoiding
the review gate rather than servicing it faster. Ranked: a worked example with
printed intermediates beats an internal identity beats a stated numerical result
beats a digitised figure. Only the last needs a human.

**The verifier** stays on the maintainer's explicit instruction. It roughly
doubles the cost per page and is aimed at the one failure this repository has
actually suffered: plausible, confident agreement that turned out to be circular.

**The two carried decisions** are genuinely blocked on the maintainer, but
neither blocks work — ask once, record the answer when it comes, continue.

## State at the time of writing (2026-07-31, commit `cdcfd4b`)

| | |
|---|---|
| Published pages | 15 (16 catalogue entries) |
| Awaiting judgement | `B1.2` (~2 min, no images), `A3.4` (3 overlays), `A1.1` (5 overlays) |
| Optional follow-ups on live pages | `E2.1`, `J3.4`, `J4.8` |
| Covered by another page | `A1.2`, `A1.3`, `A1.4` — all by `A1.1` |
| Needing a paper | 10, all section A, mostly pre-DOI classics |
| Not yet reached | 232 |

`B1.2` Aris is the cheapest win available: a two-minute yes publishes a finished
page unchanged. `A1.1` Ergun costs about fifteen minutes and unblocks four
catalogue entries at once.
